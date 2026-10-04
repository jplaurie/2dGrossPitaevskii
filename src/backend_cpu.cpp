#include "backend.hpp"
#include "fftw_utils.hpp"
#include "parallel.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
#ifdef GP2D_HAVE_FFTW_THREADS
bool fftwThreadsInitialized = false;
#endif

class CpuBackend final : public NonlinearBackend {
  public:
    explicit CpuBackend(const Parameters &parameters)
        : parameters_(parameters), paddedSpectrum_(parameters.mx() * parameters.my()),
          psi_(parameters.mx() * parameters.my()),
          squareSpectrum_(parameters.mx() * parameters.my()),
          projectedSquare_(parameters.mx() * parameters.my()),
          paddedIndexByBaseMode_(parameters.nx * parameters.ny),
          nonlinearMultiplier_(parameters.nx * parameters.ny) {
        inversePsi_ = makePlan(paddedSpectrum_, psi_, FFTW_BACKWARD);
        forwardSquare_ = makePlan(projectedSquare_, squareSpectrum_, FFTW_FORWARD);
        inverseSquare_ = makePlan(squareSpectrum_, projectedSquare_, FFTW_BACKWARD);
        forwardNonlinear_ = makePlan(paddedSpectrum_, squareSpectrum_, FFTW_FORWARD);
        if (!inversePsi_ || !forwardSquare_ || !inverseSquare_ || !forwardNonlinear_)
            throw std::runtime_error("FFTW could not create dealiased plans");
        const double scale = 1.0 / static_cast<double>(parameters.mx() * parameters.my());
        for (std::size_t i = 0; i < paddedIndexByBaseMode_.size(); ++i) {
            paddedIndexByBaseMode_[i] = paddedIndexForBaseMode(parameters, i);
            const double gamma = ginzburgLandauFactor(
                parameters, waveNumberSquared(parameters, i % parameters.nx, i / parameters.nx));
            nonlinearMultiplier_[i] =
                parameters.nonlinearityCoefficient * scale / Complex(-gamma, 1.0);
        }
    }

    void evaluate(const SpectralField &wavefunction, SpectralField &result) override {
        if (wavefunction.size() != parameters_.nx * parameters_.ny)
            throw std::runtime_error("invalid nonlinear input size");
        std::fill(paddedSpectrum_.begin(), paddedSpectrum_.end(), Complex{});
        forEachIndex(paddedIndexByBaseMode_.size(), [&](std::size_t i) {
            paddedSpectrum_[paddedIndexByBaseMode_[i]] = wavefunction[i];
        });
        fftw_execute(inversePsi_);
        const std::size_t paddedCount = parameters_.mx() * parameters_.my();
        squarePointwise(psi_, projectedSquare_, paddedCount);
        fftw_execute(forwardSquare_);
        filterPaddedSpectrum(squareSpectrum_, parameters_, 0, parameters_.my());
        fftw_execute(inverseSquare_);
        multiplyConjugatePointwise(projectedSquare_, psi_, paddedSpectrum_, paddedCount);
        fftw_execute(forwardNonlinear_);
        result.resize(parameters_.nx * parameters_.ny);
        forEachIndex(paddedIndexByBaseMode_.size(), [&](std::size_t i) {
            result[i] = nonlinearMultiplier_[i] * squareSpectrum_[paddedIndexByBaseMode_[i]];
        });
    }

  private:
    fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output, int direction) {
        return fftw_plan_dft_2d(static_cast<int>(parameters_.my()),
                                static_cast<int>(parameters_.mx()), fftwData(input),
                                fftwData(output), direction, FFTW_ESTIMATE);
    }

    Parameters parameters_;
    FftwComplexField paddedSpectrum_, psi_, squareSpectrum_, projectedSquare_;
    std::vector<std::size_t> paddedIndexByBaseMode_;
    SpectralField nonlinearMultiplier_;
    FftwPlan inversePsi_, forwardSquare_, inverseSquare_, forwardNonlinear_;
};
} // namespace

void backendInitialize(int &, char **&) {}
void backendFinalize() {
#ifdef GP2D_HAVE_FFTW_THREADS
    if (fftwThreadsInitialized) {
        fftw_cleanup_threads();
        fftwThreadsInitialized = false;
    }
#endif
}
void backendAbort(int) {}
void backendBarrier() {}
bool backendIsRoot() { return true; }
const char *backendName() {
#ifdef _OPENMP
    return "CPU/OpenMP";
#else
    return "CPU";
#endif
}
std::uint64_t backendSynchronizeSeed(std::uint64_t seed) { return seed; }

std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &parameters) {
    int threads = 1;
#ifdef _OPENMP
    threads = parameters.threadCount > 0 ? parameters.threadCount : omp_get_max_threads();
    omp_set_num_threads(threads);
#endif
    if (threads > 1) {
#ifdef GP2D_HAVE_FFTW_THREADS
        if (!fftw_init_threads())
            throw std::runtime_error("FFTW thread initialization failed");
        fftwThreadsInitialized = true;
        fftw_plan_with_nthreads(threads);
#endif
    }
    return std::make_unique<CpuBackend>(parameters);
}
