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
bool fftwThreadsInitialized = false;

class CpuBackend final : public NonlinearBackend {
public:
  explicit CpuBackend(const Parameters &p)
      : p_(p), paddedHat_(p.mx() * p.my()), psi_(p.mx() * p.my()),
        squareHat_(p.mx() * p.my()), square_(p.mx() * p.my()),
        baseToPadded_(p.nx * p.ny), nonlinearFactor_(p.nx * p.ny) {
    inversePsi_ = makePlan(paddedHat_, psi_, FFTW_BACKWARD);
    forwardSquare_ = makePlan(square_, squareHat_, FFTW_FORWARD);
    inverseSquare_ = makePlan(squareHat_, square_, FFTW_BACKWARD);
    forwardNonlinear_ = makePlan(paddedHat_, squareHat_, FFTW_FORWARD);
    if (!inversePsi_ || !forwardSquare_ || !inverseSquare_ ||
        !forwardNonlinear_)
      throw std::runtime_error("FFTW could not create dealiased plans");
    const double scale = 1.0 / static_cast<double>(p.mx() * p.my());
    for (std::size_t i = 0; i < baseToPadded_.size(); ++i) {
      baseToPadded_[i] = paddedIndexForBaseMode(p, i);
      const double gamma =
          ginzburgLandauFactor(p, waveNumberSquared(p, i % p.nx, i / p.nx));
      nonlinearFactor_[i] =
          p.nonlinearityCoefficient * scale / Complex(-gamma, 1.0);
    }
  }

  void evaluate(const SpectralField &wavefunction,
                SpectralField &result) override {
    if (wavefunction.size() != p_.nx * p_.ny)
      throw std::runtime_error("invalid nonlinear input size");
    std::fill(paddedHat_.begin(), paddedHat_.end(), Complex{});
    forEachIndex(baseToPadded_.size(), [&](std::size_t i) {
      paddedHat_[baseToPadded_[i]] = wavefunction[i];
    });
    fftw_execute(inversePsi_);
    const std::size_t paddedCount = p_.mx() * p_.my();
    squarePointwise(psi_, square_, paddedCount);
    fftw_execute(forwardSquare_);
    filterPaddedSpectrum(squareHat_, p_, 0, p_.my());
    fftw_execute(inverseSquare_);
    multiplyConjugatePointwise(square_, psi_, paddedHat_, paddedCount);
    fftw_execute(forwardNonlinear_);
    result.resize(p_.nx * p_.ny);
    forEachIndex(baseToPadded_.size(), [&](std::size_t i) {
      result[i] = nonlinearFactor_[i] * squareHat_[baseToPadded_[i]];
    });
  }

private:
  fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output,
                     int direction) {
    return fftw_plan_dft_2d(static_cast<int>(p_.my()),
                            static_cast<int>(p_.mx()), fftwData(input),
                            fftwData(output), direction, FFTW_ESTIMATE);
  }

  Parameters p_;
  FftwComplexField paddedHat_, psi_, squareHat_, square_;
  std::vector<std::size_t> baseToPadded_;
  SpectralField nonlinearFactor_;
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

std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &p) {
  int threads = 1;
#ifdef _OPENMP
  threads = p.threadCount > 0 ? p.threadCount : omp_get_max_threads();
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
  return std::make_unique<CpuBackend>(p);
}
