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
        for (FftwComplexField *buffer :
             {&paddedSpectrum_, &psi_, &squareSpectrum_, &projectedSquare_})
            buffer->resize(static_cast<std::size_t>(allocLocal_));
        inversePsi_ = makePlan(paddedSpectrum_, psi_, FFTW_BACKWARD);
        forwardSquare_ = makePlan(projectedSquare_, squareSpectrum_, FFTW_FORWARD);
        inverseSquare_ = makePlan(squareSpectrum_, projectedSquare_, FFTW_BACKWARD);
        forwardNonlinear_ = makePlan(paddedSpectrum_, squareSpectrum_, FFTW_FORWARD);
        if (!inversePsi_ || !forwardSquare_ || !inverseSquare_ || !forwardNonlinear_)
            throw std::runtime_error("FFTW-MPI could not create dealiased plans");
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
        }
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

    NoiseLayout noiseLayout() const override { return NoiseLayout::fullField; }

    void initializeTimeStepping(const IntegrationCoefficients &coefficients,
                                const SpectralField &deterministicForcing,
                                const SpectralField &state,
                                const std::vector<std::size_t> &) override {
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
        localNoise_.resize(state_.size());
        workspace_.initialize(state_.size(), parameters_.nonlinearStageCount());
    }

    void advanceTimeStep(const SpectralField &noise) override {
        if (!noise.empty() && noise.size() != parameters_.nx * parameters_.ny)
            throw std::runtime_error("invalid MPI stochastic-noise field");
        std::fill(localNoise_.begin(), localNoise_.end(), Complex{});
        if (!noise.empty())
            for (const LocalMode &mode : localModes_)
                localNoise_[mode.paddedIndex] = noise[mode.baseIndex];
        const auto rhs = [&](const SpectralField &input, SpectralField &output) {
            evaluateLocal(input, output);
            if (!deterministicForcing_.empty())
                forEachIndex(output.size(),
                             [&](std::size_t i) { output[i] += deterministicForcing_[i]; });
        };
        advanceHostTimeStep(
            parameters_, coefficients_, state_, noise.empty() ? SpectralField{} : localNoise_,
            workspace_, rhs, [&](SpectralField &localState) {
                if (parameters_.hypoviscosity > 0.0 && parameters_.hypoviscosityOrder < 0.0)
                    for (const LocalMode &mode : localModes_)
                        if (mode.baseIndex == 0)
                            localState[mode.paddedIndex] = Complex{};
            });
    }

    void downloadState(SpectralField &state) override { reduceLocalField(state_, state); }

    void evaluateCurrent(SpectralField &nonlinearTerm) override {
        SpectralField localResult;
        evaluateLocal(state_, localResult);
        reduceLocalField(localResult, nonlinearTerm);
    }

  private:
    void evaluateLocal(const SpectralField &localState, SpectralField &result) {
        if (localState.size() != static_cast<std::size_t>(allocLocal_))
            throw std::runtime_error("invalid distributed nonlinear input size");
        std::copy(localState.begin(), localState.end(), paddedSpectrum_.begin());
        fftw_execute(inversePsi_);
        const std::size_t localCount = static_cast<std::size_t>(localRows_) * parameters_.mx();
        squarePointwise(psi_, projectedSquare_, localCount);
        fftw_execute(forwardSquare_);
        filterPaddedSpectrum(squareSpectrum_, parameters_, static_cast<std::size_t>(firstRow_),
                             static_cast<std::size_t>(localRows_));
        fftw_execute(inverseSquare_);
        multiplyConjugatePointwise(projectedSquare_, psi_, paddedSpectrum_, localCount);
        fftw_execute(forwardNonlinear_);
        result.assign(static_cast<std::size_t>(allocLocal_), Complex{});
        forEachIndex(localModes_.size(), [&](std::size_t i) {
            const LocalMode &mode = localModes_[i];
            result[mode.paddedIndex] = mode.nonlinearFactor * squareSpectrum_[mode.paddedIndex];
        });
    }

    void reduceLocalField(const SpectralField &localField, SpectralField &globalField) const {
        SpectralField localBase(parameters_.nx * parameters_.ny);
        for (const LocalMode &mode : localModes_)
            localBase[mode.baseIndex] = localField[mode.paddedIndex];
        if (rank == 0)
            globalField.resize(localBase.size());
        else
            globalField.clear();
        mpiCheck(MPI_Reduce(localBase.data(), rank == 0 ? globalField.data() : nullptr,
                            static_cast<int>(localBase.size()), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0,
                            MPI_COMM_WORLD),
                 "MPI_Reduce distributed spectral field");
    }

    fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output, int direction) {
        return fftw_mpi_plan_dft_2d(
            static_cast<ptrdiff_t>(parameters_.my()), static_cast<ptrdiff_t>(parameters_.mx()),
            fftwData(input), fftwData(output), MPI_COMM_WORLD, direction, fftwPlanningFlags());
    }

    Parameters parameters_;
    ptrdiff_t allocLocal_ = 0, localRows_ = 0, firstRow_ = 0;
    FftwComplexField paddedSpectrum_, psi_, squareSpectrum_, projectedSquare_;
    SpectralField state_, deterministicForcing_, localNoise_;
    IntegrationCoefficients coefficients_;
    HostIntegrationWorkspace workspace_;
    std::vector<LocalMode> localModes_;
    FftwPlan inversePsi_, forwardSquare_, inverseSquare_, forwardNonlinear_;
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

void backendFinalize() {
    fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
    if (rank == 0)
        saveFftwWisdom();
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
