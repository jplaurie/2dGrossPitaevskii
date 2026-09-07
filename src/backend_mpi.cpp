#include "backend.hpp"
#include "fftw_utils.hpp"
#include "spectral.hpp"

#include <fftw3-mpi.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
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

void mpiCheck(int status, const char *operation) {
  if (status == MPI_SUCCESS)
    return;
  char message[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  MPI_Error_string(status, message, &length);
  throw std::runtime_error(std::string(operation) + ": " +
                           std::string(message, length));
}

class MpiBackend final : public NonlinearBackend {
public:
  explicit MpiBackend(const Parameters &p) : p_(p) {
    const std::size_t baseCount = p.nx * p.ny;
    if (baseCount > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::runtime_error("spectral field exceeds MPI count limit");
    allocLocal_ = fftw_mpi_local_size_2d(
        static_cast<ptrdiff_t>(p.my()), static_cast<ptrdiff_t>(p.mx()),
        MPI_COMM_WORLD, &localRows_, &firstRow_);
    if (allocLocal_ < localRows_ * static_cast<ptrdiff_t>(p.mx()))
      throw std::runtime_error("FFTW-MPI returned an invalid local size");
    for (FftwComplexField *buffer : {&paddedHat_, &psi_, &squareHat_, &square_})
      buffer->resize(static_cast<std::size_t>(allocLocal_));
    local_.resize(baseCount);
    inversePsi_ = makePlan(paddedHat_, psi_, FFTW_BACKWARD);
    forwardSquare_ = makePlan(square_, squareHat_, FFTW_FORWARD);
    inverseSquare_ = makePlan(squareHat_, square_, FFTW_BACKWARD);
    forwardNonlinear_ = makePlan(paddedHat_, squareHat_, FFTW_FORWARD);
    if (!inversePsi_ || !forwardSquare_ || !inverseSquare_ ||
        !forwardNonlinear_) {
      releasePlans();
      throw std::runtime_error("FFTW-MPI could not create dealiased plans");
    }
  }

  ~MpiBackend() override { releasePlans(); }

  void evaluate(const SpectralField &wavefunction,
                SpectralField &result) override {
    if (wavefunction.size() != p_.nx * p_.ny)
      throw std::runtime_error("invalid nonlinear input size");
    std::fill(paddedHat_.begin(), paddedHat_.end(), Complex{});
    for (std::size_t y = 0; y < p_.ny; ++y) {
      const ptrdiff_t py =
          static_cast<ptrdiff_t>(paddedIndexForBase(y, p_.ny, p_.my()));
      if (py < firstRow_ || py >= firstRow_ + localRows_)
        continue;
      const std::size_t localY = static_cast<std::size_t>(py - firstRow_);
      for (std::size_t x = 0; x < p_.nx; ++x) {
        const std::size_t px = paddedIndexForBase(x, p_.nx, p_.mx());
        paddedHat_[spectralIndex(px, localY, p_.mx())] =
            wavefunction[spectralIndex(x, y, p_.nx)];
      }
    }
    fftw_execute(inversePsi_);
    const std::size_t localCount =
        static_cast<std::size_t>(localRows_) * p_.mx();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (localCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(localCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      square_[i] = psi_[i] * psi_[i];
    }
    fftw_execute(forwardSquare_);
    const double scale = 1.0 / static_cast<double>(p_.mx() * p_.my());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (localCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(localCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      const std::size_t px = i % p_.mx();
      const std::size_t localY = i / p_.mx();
      const std::size_t py = static_cast<std::size_t>(firstRow_) + localY;
      if (retainedPaddedWave(signedWave(px, p_.mx()), p_.nx) &&
          retainedPaddedWave(signedWave(py, p_.my()), p_.ny))
        squareHat_[i] *= scale;
      else
        squareHat_[i] = Complex{};
    }
    fftw_execute(inverseSquare_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (localCount >= 16384)
#endif
    for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(localCount);
         ++raw) {
      const auto i = static_cast<std::size_t>(raw);
      paddedHat_[i] = square_[i] * std::conj(psi_[i]);
    }
    fftw_execute(forwardNonlinear_);
    std::fill(local_.begin(), local_.end(), Complex{});
    for (std::size_t y = 0; y < p_.ny; ++y) {
      const ptrdiff_t py =
          static_cast<ptrdiff_t>(paddedIndexForBase(y, p_.ny, p_.my()));
      if (py < firstRow_ || py >= firstRow_ + localRows_)
        continue;
      const std::size_t localY = static_cast<std::size_t>(py - firstRow_);
      for (std::size_t x = 0; x < p_.nx; ++x) {
        const std::size_t px = paddedIndexForBase(x, p_.nx, p_.mx());
        const double gamma =
            ginzburgLandauFactor(p_, waveNumberSquared(p_, x, y));
        local_[spectralIndex(x, y, p_.nx)] =
            p_.nonlinearityCoefficient * scale *
            squareHat_[spectralIndex(px, localY, p_.mx())] /
            Complex(-gamma, 1.0);
      }
    }
    result.resize(local_.size());
    mpiCheck(MPI_Allreduce(local_.data(), result.data(),
                           static_cast<int>(local_.size()),
                           MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD),
             "MPI_Allreduce nonlinear result");
  }

private:
  fftw_plan makePlan(FftwComplexField &input, FftwComplexField &output,
                     int direction) {
    return fftw_mpi_plan_dft_2d(static_cast<ptrdiff_t>(p_.my()),
                                static_cast<ptrdiff_t>(p_.mx()),
                                fftwData(input), fftwData(output),
                                MPI_COMM_WORLD, direction, FFTW_ESTIMATE);
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
  ptrdiff_t allocLocal_ = 0, localRows_ = 0, firstRow_ = 0;
  FftwComplexField paddedHat_, psi_, squareHat_, square_;
  SpectralField local_;
  fftw_plan inversePsi_ = nullptr, forwardSquare_ = nullptr,
            inverseSquare_ = nullptr, forwardNonlinear_ = nullptr;
};
} // namespace

void backendInitialize(int &argc, char **&argv) {
  int provided = MPI_THREAD_SINGLE;
#ifdef _OPENMP
  constexpr int required = MPI_THREAD_FUNNELED;
#else
  constexpr int required = MPI_THREAD_SINGLE;
#endif
  mpiCheck(MPI_Init_thread(&argc, &argv, required, &provided),
           "MPI_Init_thread");
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
  mpiCheck(MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD),
           "MPI_Bcast random seed");
  return seed;
}

std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &p) {
  int threads = 1;
#ifdef _OPENMP
  threads = p.threadCount > 0 ? p.threadCount : omp_get_max_threads();
  omp_set_num_threads(threads);
#endif
  if (threads > 1) {
#ifdef GP2D_HAVE_FFTW_THREADS
    fftw_plan_with_nthreads(threads);
#endif
  }
  return std::make_unique<MpiBackend>(p);
}
