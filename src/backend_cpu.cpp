#include "backend.hpp"
#include "fftw_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
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
        squareHat_(p.mx() * p.my()), square_(p.mx() * p.my()) {
    inversePsi_ = makePlan(paddedHat_, psi_, FFTW_BACKWARD);
    forwardSquare_ = makePlan(square_, squareHat_, FFTW_FORWARD);
    inverseSquare_ = makePlan(squareHat_, square_, FFTW_BACKWARD);
    forwardNonlinear_ = makePlan(paddedHat_, squareHat_, FFTW_FORWARD);
    if (!inversePsi_ || !forwardSquare_ || !inverseSquare_ ||
        !forwardNonlinear_) {
      releasePlans();
      throw std::runtime_error("FFTW could not create dealiased plans");
    }
  }

  ~CpuBackend() override { releasePlans(); }

  void evaluate(const SpectralField &wavefunction,
                SpectralField &result) override {
    if (wavefunction.size() != p_.nx * p_.ny)
      throw std::runtime_error("invalid nonlinear input size");
    std::fill(paddedHat_.begin(), paddedHat_.end(), Complex{});
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (p_.nx * p_.ny >= 16384)
#endif
    for (std::ptrdiff_t rawY = 0; rawY < static_cast<std::ptrdiff_t>(p_.ny);
         ++rawY) {
      const auto y = static_cast<std::size_t>(rawY);
      const std::size_t py = paddedIndexForBase(y, p_.ny, p_.my());
      for (std::size_t x = 0; x < p_.nx; ++x) {
        const std::size_t px = paddedIndexForBase(x, p_.nx, p_.mx());
        paddedHat_[spectralIndex(px, py, p_.mx())] =
            wavefunction[spectralIndex(x, y, p_.nx)];
      }
    }
    fftw_execute(inversePsi_);
    const std::size_t paddedCount = p_.mx() * p_.my();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (paddedCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(paddedCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      square_[i] = psi_[i] * psi_[i];
    }
    fftw_execute(forwardSquare_);
    const double scale = 1.0 / static_cast<double>(paddedCount);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (paddedCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(paddedCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      const std::size_t px = i % p_.mx();
      const std::size_t py = i / p_.mx();
      if (retainedPaddedWave(signedWave(px, p_.mx()), p_.nx) &&
          retainedPaddedWave(signedWave(py, p_.my()), p_.ny))
        squareHat_[i] *= scale;
      else
        squareHat_[i] = Complex{};
    }
    fftw_execute(inverseSquare_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (paddedCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(paddedCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      paddedHat_[i] = square_[i] * std::conj(psi_[i]);
    }
    fftw_execute(forwardNonlinear_);
    result.resize(p_.nx * p_.ny);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (p_.nx * p_.ny >= 16384)
#endif
    for (std::ptrdiff_t rawY = 0; rawY < static_cast<std::ptrdiff_t>(p_.ny);
         ++rawY) {
      const auto y = static_cast<std::size_t>(rawY);
      const std::size_t py = paddedIndexForBase(y, p_.ny, p_.my());
      for (std::size_t x = 0; x < p_.nx; ++x) {
        const std::size_t px = paddedIndexForBase(x, p_.nx, p_.mx());
        const double gamma =
            ginzburgLandauFactor(p_, waveNumberSquared(p_, x, y));
        result[spectralIndex(x, y, p_.nx)] =
            p_.nonlinearityCoefficient * scale *
            squareHat_[spectralIndex(px, py, p_.mx())] / Complex(-gamma, 1.0);
      }
    }
  }

private:
  fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output,
                     int direction) {
    return fftw_plan_dft_2d(static_cast<int>(p_.my()),
                            static_cast<int>(p_.mx()), fftwData(input),
                            fftwData(output), direction, FFTW_ESTIMATE);
  }

  void releasePlans() noexcept {
    for (fftw_plan *plan :
         {&inversePsi_, &forwardSquare_, &inverseSquare_, &forwardNonlinear_}) {
      if (*plan)
        fftw_destroy_plan(*plan);
      *plan = nullptr;
    }
  }

  Parameters p_;
  FftwComplexField paddedHat_, psi_, squareHat_, square_;
  fftw_plan inversePsi_ = nullptr, forwardSquare_ = nullptr,
            inverseSquare_ = nullptr, forwardNonlinear_ = nullptr;
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
