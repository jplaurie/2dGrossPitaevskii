#include "backend.hpp"
#include "fftw_utils.hpp"
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
        !forwardNonlinear_)
      throw std::runtime_error("FFTW-MPI could not create dealiased plans");
    const double scale = 1.0 / static_cast<double>(p.mx() * p.my());
    for (std::size_t i = 0; i < baseCount; ++i) {
      const std::size_t padded = paddedIndexForBaseMode(p, i);
      const ptrdiff_t paddedRow = static_cast<ptrdiff_t>(padded / p.mx());
      if (paddedRow < firstRow_ || paddedRow >= firstRow_ + localRows_)
        continue;
      const std::size_t localIndex =
          padded - static_cast<std::size_t>(firstRow_) * p.mx();
      const double gamma =
          ginzburgLandauFactor(p, waveNumberSquared(p, i % p.nx, i / p.nx));
      localModes_.push_back(
          {i, localIndex,
           p.nonlinearityCoefficient * scale / Complex(-gamma, 1.0)});
    }
  }

  void evaluate(const SpectralField &wavefunction,
                SpectralField &result) override {
    if (wavefunction.size() != p_.nx * p_.ny)
      throw std::runtime_error("invalid nonlinear input size");
    std::fill(paddedHat_.begin(), paddedHat_.end(), Complex{});
    forEachIndex(localModes_.size(), [&](std::size_t i) {
      const LocalMode &mode = localModes_[i];
      paddedHat_[mode.paddedIndex] = wavefunction[mode.baseIndex];
    });
    fftw_execute(inversePsi_);
    const std::size_t localCount =
        static_cast<std::size_t>(localRows_) * p_.mx();
    squarePointwise(psi_, square_, localCount);
    fftw_execute(forwardSquare_);
    filterPaddedSpectrum(squareHat_, p_, static_cast<std::size_t>(firstRow_),
                         static_cast<std::size_t>(localRows_));
    fftw_execute(inverseSquare_);
    multiplyConjugatePointwise(square_, psi_, paddedHat_, localCount);
    fftw_execute(forwardNonlinear_);
    std::fill(local_.begin(), local_.end(), Complex{});
    forEachIndex(localModes_.size(), [&](std::size_t i) {
      const LocalMode &mode = localModes_[i];
      local_[mode.baseIndex] =
          mode.nonlinearFactor * squareHat_[mode.paddedIndex];
    });
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

  Parameters p_;
  ptrdiff_t allocLocal_ = 0, localRows_ = 0, firstRow_ = 0;
  FftwComplexField paddedHat_, psi_, squareHat_, square_;
  SpectralField local_;
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
