#include "backend.hpp"
#include "fftw_utils.hpp"
#include "host_stepper.hpp"
#include "parallel.hpp"
#include "spectral.hpp"

#include <fftw3-mpi.h>
#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
int rank = 0;
bool fftwThreadsInitialized = false;

struct LocalMode {
    std::size_t baseIndex, paddedIndex;
    Complex nonlinearFactor;
};

void mpiCheck(int status, const char *operation) {
    if (status == MPI_SUCCESS)
        return;
    char message[MPI_MAX_ERROR_STRING]{};
    int length = 0;
    MPI_Error_string(status, message, &length);
    throw std::runtime_error(std::string(operation) + ": " + std::string(message, length));
}

class MpiBackend final : public NonlinearBackend {
  public:
    explicit MpiBackend(const Parameters &parameters) : parameters_(parameters) {
        const std::size_t baseCount = parameters.nx * parameters.ny;
        if (baseCount > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("spectral field exceeds MPI count limit");
        allocLocal_ = fftw_mpi_local_size_2d(static_cast<ptrdiff_t>(parameters.my()),
                                             static_cast<ptrdiff_t>(parameters.mx()),
                                             MPI_COMM_WORLD, &localRows_, &firstRow_);
        if (allocLocal_ < localRows_ * static_cast<ptrdiff_t>(parameters.mx()))
            throw std::runtime_error("FFTW-MPI returned an invalid local size");
        if (parameters.nonlinearityCoefficient != 0.0) {
            for (FftwComplexField *buffer :
                 {&paddedSpectrum_, &psi_, &densitySpectrum_, &projectedDensity_})
                buffer->resize(static_cast<std::size_t>(allocLocal_));
            inversePsi_ = makePlan(paddedSpectrum_, psi_, FFTW_BACKWARD);
            forwardDensity_ = makePlan(projectedDensity_, densitySpectrum_, FFTW_FORWARD);
            inverseDensity_ = makePlan(densitySpectrum_, projectedDensity_, FFTW_BACKWARD);
            forwardNonlinear_ = makePlan(paddedSpectrum_, densitySpectrum_, FFTW_FORWARD);
            if (!inversePsi_ || !forwardDensity_ || !inverseDensity_ || !forwardNonlinear_)
                throw std::runtime_error("FFTW-MPI could not create dealiased plans");
        }
        const double scale = 1.0 / static_cast<double>(parameters.mx() * parameters.my());
        for (std::size_t i = 0; i < baseCount; ++i) {
            const std::size_t padded = paddedIndexForBaseMode(parameters, i);
            const ptrdiff_t paddedRow = static_cast<ptrdiff_t>(padded / parameters.mx());
            if (paddedRow < firstRow_ || paddedRow >= firstRow_ + localRows_)
                continue;
            const std::size_t localIndex =
                padded - static_cast<std::size_t>(firstRow_) * parameters.mx();
            const double gamma = ginzburgLandauFactor(
                parameters, waveNumberSquared(parameters, i % parameters.nx, i / parameters.nx));
            localModes_.push_back(
                {i, localIndex, parameters.nonlinearityCoefficient * scale / Complex(-gamma, 1.0)});
            if (i == 0)
                meanIndex_ = localIndex;
        }
        initializeOutputGather(baseCount);
    }

    void evaluate(const SpectralField &wavefunction, SpectralField &result) override {
        if (wavefunction.size() != parameters_.nx * parameters_.ny)
            throw std::runtime_error("invalid nonlinear input size");
        SpectralField localState(static_cast<std::size_t>(allocLocal_));
        for (const LocalMode &mode : localModes_)
            localState[mode.paddedIndex] = wavefunction[mode.baseIndex];
        SpectralField localResult;
        evaluateLocal(localState, localResult);
        SpectralField localBase(parameters_.nx * parameters_.ny);
        for (const LocalMode &mode : localModes_)
            localBase[mode.baseIndex] = localResult[mode.paddedIndex];
        result.resize(localBase.size());
        mpiCheck(MPI_Allreduce(localBase.data(), result.data(), static_cast<int>(localBase.size()),
                               MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD),
                 "MPI_Allreduce standalone nonlinear result");
    }

    NoiseLayout noiseLayout() const override { return NoiseLayout::forcedModes; }

    std::vector<std::size_t>
    ownedNoiseIndices(const std::vector<std::size_t> &indices) const override {
        std::vector<std::size_t> owned;
        for (const auto index : indices) {
            const auto row = static_cast<ptrdiff_t>(paddedIndexForBaseMode(parameters_, index) /
                                                    parameters_.mx());
            if (row >= firstRow_ && row < firstRow_ + localRows_)
                owned.push_back(index);
        }
        return owned;
    }

    void initializeTimeStepping(const IntegrationCoefficients &coefficients,
                                const SpectralField &deterministicForcing,
                                const SpectralField &state,
                                const std::vector<std::size_t> &noiseIndices) override {
        const std::size_t baseCount = parameters_.nx * parameters_.ny;
        if (state.size() != baseCount)
            throw std::runtime_error("invalid initial state size");
        state_.assign(static_cast<std::size_t>(allocLocal_), Complex{});
        if (!deterministicForcing.empty())
            deterministicForcing_.assign(state_.size(), Complex{});
        const auto globalFields = coefficients.fields();
        const auto localFields = coefficients_.fields();
        for (std::size_t field = 0; field < globalFields.size(); ++field)
            if (!globalFields[field]->empty())
                localFields[field]->assign(state_.size(), Complex{});
        for (const LocalMode &mode : localModes_) {
            state_[mode.paddedIndex] = state[mode.baseIndex];
            if (!deterministicForcing.empty())
                deterministicForcing_[mode.paddedIndex] = deterministicForcing[mode.baseIndex];
            for (std::size_t field = 0; field < globalFields.size(); ++field)
                if (!globalFields[field]->empty())
                    (*localFields[field])[mode.paddedIndex] =
                        (*globalFields[field])[mode.baseIndex];
        }
        noisePaddedIndices_.clear();
        localNoise_.clear();
        if (parameters_.forcingEnabled &&
            parameters_.forcingProfile != ForcingProfile::singleMode) {
            localNoise_.assign(state_.size(), Complex{});
            for (const auto index : noiseIndices) {
                const auto padded = paddedIndexForBaseMode(parameters_, index);
                const auto row = static_cast<ptrdiff_t>(padded / parameters_.mx());
                if (row < firstRow_ || row >= firstRow_ + localRows_)
                    throw std::runtime_error("MPI noise mode is not locally owned");
                noisePaddedIndices_.push_back(padded - static_cast<std::size_t>(firstRow_) *
                                                           parameters_.mx());
            }
        }
        workspace_.initialize(state_.size(), parameters_.nonlinearStageCount());
    }

    void advanceTimeStep(const SpectralField &noise) override {
        if (!noise.empty() &&
            (noise.size() != noisePaddedIndices_.size() || localNoise_.size() != state_.size()))
            throw std::runtime_error("invalid MPI stochastic-noise field");
        if (!noise.empty()) {
            forEachIndex(noise.size(),
                         [&](std::size_t i) { localNoise_[noisePaddedIndices_[i]] = noise[i]; });
        }
        const SpectralField &stepNoise = noise.empty() ? noise : localNoise_;
        const auto rhs = [&](const SpectralField &input, SpectralField &output) {
            evaluateLocal(input, output);
            if (!deterministicForcing_.empty())
                forEachIndex(output.size(),
                             [&](std::size_t i) { output[i] += deterministicForcing_[i]; });
        };
        advanceHostTimeStep(parameters_, coefficients_, state_, stepNoise, workspace_, rhs,
                            [&](SpectralField &localState) {
                                if (parameters_.hypoviscosity > 0.0 &&
                                    parameters_.hypoviscosityOrder < 0.0 &&
                                    meanIndex_ != std::numeric_limits<std::size_t>::max())
                                    localState[meanIndex_] = Complex{};
                            });
    }

    void downloadState(SpectralField &state) override { gatherLocalField(state_, state); }

    void evaluateCurrent(SpectralField &nonlinearTerm) override {
        SpectralField localResult;
        evaluateLocal(state_, localResult);
        gatherLocalField(localResult, nonlinearTerm);
    }

  private:
    void evaluateLocal(const SpectralField &localState, SpectralField &result) {
        if (localState.size() != static_cast<std::size_t>(allocLocal_))
            throw std::runtime_error("invalid distributed nonlinear input size");
        result.assign(static_cast<std::size_t>(allocLocal_), Complex{});
        if (parameters_.nonlinearityCoefficient == 0.0)
            return;
        std::copy(localState.begin(), localState.end(), paddedSpectrum_.begin());
        fftw_execute(inversePsi_);
        const std::size_t localCount = static_cast<std::size_t>(localRows_) * parameters_.mx();
        densityPointwise(psi_, projectedDensity_, localCount);
        fftw_execute(forwardDensity_);
        filterDensitySpectrum(densitySpectrum_, parameters_, static_cast<std::size_t>(firstRow_),
                              static_cast<std::size_t>(localRows_));
        fftw_execute(inverseDensity_);
        multiplyPointwise(projectedDensity_, psi_, paddedSpectrum_, localCount);
        fftw_execute(forwardNonlinear_);
        forEachIndex(localModes_.size(), [&](std::size_t i) {
            const LocalMode &mode = localModes_[i];
            result[mode.paddedIndex] = mode.nonlinearFactor * densitySpectrum_[mode.paddedIndex];
        });
    }

    void initializeOutputGather(std::size_t baseCount) {
        int ranks = 0;
        mpiCheck(MPI_Comm_size(MPI_COMM_WORLD, &ranks), "MPI_Comm_size");
        const int localCount = static_cast<int>(localModes_.size());
        if (rank == 0) {
            outputCounts_.resize(ranks);
            outputDisplacements_.resize(ranks);
        }
        mpiCheck(MPI_Gather(&localCount, 1, MPI_INT, outputCounts_.data(), 1, MPI_INT, 0,
                            MPI_COMM_WORLD),
                 "MPI_Gather output counts");
        if (rank == 0) {
            int total = 0;
            for (int i = 0; i < ranks; ++i) {
                outputDisplacements_[i] = total;
                total += outputCounts_[i];
            }
            if (static_cast<std::size_t>(total) != baseCount)
                throw std::runtime_error("distributed output modes do not cover the base grid");
            outputBaseIndices_.resize(baseCount);
            gatheredOutput_.resize(baseCount);
        }
        std::vector<std::uint64_t> indices(localModes_.size());
        for (std::size_t i = 0; i < indices.size(); ++i)
            indices[i] = localModes_[i].baseIndex;
        mpiCheck(MPI_Gatherv(indices.data(), localCount, MPI_UINT64_T, outputBaseIndices_.data(),
                             outputCounts_.data(), outputDisplacements_.data(), MPI_UINT64_T, 0,
                             MPI_COMM_WORLD),
                 "MPI_Gatherv output indices");
        packedOutput_.resize(localModes_.size());
    }

    void gatherLocalField(const SpectralField &localField, SpectralField &globalField) {
        forEachIndex(localModes_.size(), [&](std::size_t i) {
            packedOutput_[i] = localField[localModes_[i].paddedIndex];
        });
        mpiCheck(MPI_Gatherv(packedOutput_.data(), static_cast<int>(packedOutput_.size()),
                             MPI_CXX_DOUBLE_COMPLEX, gatheredOutput_.data(), outputCounts_.data(),
                             outputDisplacements_.data(), MPI_CXX_DOUBLE_COMPLEX, 0,
                             MPI_COMM_WORLD),
                 "MPI_Gatherv distributed spectral field");
        if (rank == 0)
            globalField.resize(gatheredOutput_.size());
        else
            globalField.clear();
        forEachIndex(gatheredOutput_.size(), [&](std::size_t i) {
            globalField[outputBaseIndices_[i]] = gatheredOutput_[i];
        });
    }

    fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output, int direction) {
        return fftw_mpi_plan_dft_2d(
            static_cast<ptrdiff_t>(parameters_.my()), static_cast<ptrdiff_t>(parameters_.mx()),
            fftwData(input), fftwData(output), MPI_COMM_WORLD, direction, fftwPlanningFlags());
    }

    Parameters parameters_;
    ptrdiff_t allocLocal_ = 0, localRows_ = 0, firstRow_ = 0;
    FftwComplexField paddedSpectrum_, psi_, densitySpectrum_, projectedDensity_;
    SpectralField state_, deterministicForcing_, localNoise_;
    std::vector<std::size_t> noisePaddedIndices_;
    IntegrationCoefficients coefficients_;
    HostIntegrationWorkspace workspace_;
    std::vector<LocalMode> localModes_;
    std::size_t meanIndex_ = std::numeric_limits<std::size_t>::max();
    std::vector<int> outputCounts_, outputDisplacements_;
    std::vector<std::uint64_t> outputBaseIndices_;
    SpectralField packedOutput_, gatheredOutput_;
    FftwPlan inversePsi_, forwardDensity_, inverseDensity_, forwardNonlinear_;
};
} // namespace

void backendInitialize(int &argc, char **&argv) {
    int provided = MPI_THREAD_SINGLE;
#ifdef _OPENMP
    constexpr int required = MPI_THREAD_FUNNELED;
#else
    constexpr int required = MPI_THREAD_SINGLE;
#endif
    mpiCheck(MPI_Init_thread(&argc, &argv, required, &provided), "MPI_Init_thread");
    mpiCheck(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "MPI_Comm_rank");
    if (provided < required)
        throw std::runtime_error("MPI runtime lacks required thread support");
#if defined(GP2D_HAVE_FFTW_THREADS) && defined(_OPENMP)
    if (!fftw_init_threads())
        throw std::runtime_error("FFTW thread initialization failed");
    fftwThreadsInitialized = true;
#endif
    fftw_mpi_init();
}

void backendFinalize(bool writeWisdom) {
    std::string wisdomError;
    int failed = 0;
    if (writeWisdom) {
        fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
        if (rank == 0) {
            try {
                saveFftwWisdom();
            } catch (const std::exception &error) {
                wisdomError = error.what();
                failed = 1;
            }
        }
        mpiCheck(MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD),
                 "MPI_Bcast wisdom export status");
    }
    fftw_mpi_cleanup();
#ifdef GP2D_HAVE_FFTW_THREADS
    if (fftwThreadsInitialized) {
        fftw_cleanup_threads();
        fftwThreadsInitialized = false;
    }
#endif
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
    if (failed)
        throw std::runtime_error(rank == 0 ? wisdomError
                                           : "rank zero could not export FFTW wisdom");
}
void backendAbort(int exitCode) { MPI_Abort(MPI_COMM_WORLD, exitCode); }
void backendBarrier() { mpiCheck(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier"); }
bool backendIsRoot() { return rank == 0; }
const char *backendName() {
#ifdef _OPENMP
    return "MPI/OpenMP";
#else
    return "MPI";
#endif
}
std::uint64_t backendSynchronizeSeed(std::uint64_t seed) {
    mpiCheck(MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD), "MPI_Bcast random seed");
    return seed;
}

std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &parameters) {
    int threads = 1;
#ifdef _OPENMP
    threads = parameters.threadCount > 0 ? parameters.threadCount : omp_get_max_threads();
    omp_set_num_threads(threads);
#endif
    if (threads > 1) {
#ifdef GP2D_HAVE_FFTW_THREADS
        fftw_plan_with_nthreads(threads);
#endif
    }
    configureFftw(parameters, rank == 0);
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);
    return std::make_unique<MpiBackend>(parameters);
}
