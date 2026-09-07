#include <cuda_runtime.h>
#include <cufft.h>

__host__ __device__ inline cufftDoubleComplex operator+(cufftDoubleComplex a,
                                                        cufftDoubleComplex b) {
  return {a.x + b.x, a.y + b.y};
}
__host__ __device__ inline cufftDoubleComplex operator-(cufftDoubleComplex a,
                                                        cufftDoubleComplex b) {
  return {a.x - b.x, a.y - b.y};
}
__host__ __device__ inline cufftDoubleComplex operator*(cufftDoubleComplex a,
                                                        cufftDoubleComplex b) {
  return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
}
__host__ __device__ inline cufftDoubleComplex operator*(double a,
                                                        cufftDoubleComplex b) {
  return {a * b.x, a * b.y};
}

#include "backend.hpp"
#include "spectral.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>

static_assert(sizeof(Complex) == sizeof(cufftDoubleComplex));
static_assert(std::is_trivially_copyable_v<Complex>);

namespace {
void cudaCheck(cudaError_t status, const char *operation) {
  if (status != cudaSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             cudaGetErrorString(status));
}
void cufftCheck(cufftResult status, const char *operation) {
  if (status != CUFFT_SUCCESS)
    throw std::runtime_error(std::string(operation) + " failed (cuFFT status " +
                             std::to_string(static_cast<int>(status)) + ")");
}

template <class T> class DeviceBuffer {
public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t count) { allocate(count); }
  ~DeviceBuffer() { cudaFree(data_); }
  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;
  void allocate(std::size_t count) {
    cudaFree(data_);
    data_ = nullptr;
    count_ = count;
    if (count)
      cudaCheck(cudaMalloc(&data_, count * sizeof(T)),
                "allocate device buffer");
  }
  T *data() { return data_; }
  const T *data() const { return data_; }
  std::size_t size() const { return count_; }

private:
  T *data_ = nullptr;
  std::size_t count_ = 0;
};

__device__ long deviceWave(std::size_t index, std::size_t count) {
  return index <= (count - 1) / 2
             ? static_cast<long>(index)
             : static_cast<long>(index) - static_cast<long>(count);
}

__global__ void embedBase(const cufftDoubleComplex *input,
                          cufftDoubleComplex *padded, std::size_t nx,
                          std::size_t ny, std::size_t mx, std::size_t my) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= mx * my)
    return;
  const std::size_t px = i % mx, py = i / mx;
  const long kx = deviceWave(px, mx), ky = deviceWave(py, my);
  const bool retained =
      kx >= -static_cast<long>(nx / 2) && kx < static_cast<long>(nx / 2) &&
      ky >= -static_cast<long>(ny / 2) && ky < static_cast<long>(ny / 2);
  if (!retained) {
    padded[i] = {0.0, 0.0};
    return;
  }
  const std::size_t x =
      kx >= 0 ? static_cast<std::size_t>(kx)
              : static_cast<std::size_t>(static_cast<long>(nx) + kx);
  const std::size_t y =
      ky >= 0 ? static_cast<std::size_t>(ky)
              : static_cast<std::size_t>(static_cast<long>(ny) + ky);
  padded[i] = input[y * nx + x];
}

__global__ void squareField(const cufftDoubleComplex *input,
                            cufftDoubleComplex *output, std::size_t count) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < count)
    output[i] = input[i] * input[i];
}

__global__ void filterSquare(cufftDoubleComplex *field, std::size_t nx,
                             std::size_t ny, std::size_t mx, std::size_t my,
                             double scale) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= mx * my)
    return;
  const long kx = deviceWave(i % mx, mx), ky = deviceWave(i / mx, my);
  if (kx >= -static_cast<long>(nx / 2) && kx < static_cast<long>(nx / 2) &&
      ky >= -static_cast<long>(ny / 2) && ky < static_cast<long>(ny / 2)) {
    field[i].x *= scale;
    field[i].y *= scale;
  } else {
    field[i] = {0.0, 0.0};
  }
}

__global__ void cubicField(const cufftDoubleComplex *square,
                           const cufftDoubleComplex *psi,
                           cufftDoubleComplex *output, std::size_t count) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= count)
    return;
  const cufftDoubleComplex conjugate{psi[i].x, -psi[i].y};
  output[i] = square[i] * conjugate;
}

__global__ void extractBase(const cufftDoubleComplex *padded,
                            cufftDoubleComplex *output, std::size_t nx,
                            std::size_t ny, std::size_t mx, std::size_t my,
                            double lx, double ly, double coefficient,
                            double damping, double dampingCutoff,
                            double scale) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= nx * ny)
    return;
  const std::size_t x = i % nx, y = i / nx;
  const long kxIndex = deviceWave(x, nx), kyIndex = deviceWave(y, ny);
  const std::size_t px =
      kxIndex >= 0 ? static_cast<std::size_t>(kxIndex)
                   : static_cast<std::size_t>(static_cast<long>(mx) + kxIndex);
  const std::size_t py =
      kyIndex >= 0 ? static_cast<std::size_t>(kyIndex)
                   : static_cast<std::size_t>(static_cast<long>(my) + kyIndex);
  constexpr double twoPi = 6.283185307179586476925286766559;
  const double kx = twoPi * static_cast<double>(kxIndex) / lx;
  const double ky = twoPi * static_cast<double>(kyIndex) / ly;
  const double gamma =
      damping > 0.0 && sqrt(kx * kx + ky * ky) > dampingCutoff ? damping : 0.0;
  const cufftDoubleComplex value = padded[py * mx + px];
  const double denominator = gamma * gamma + 1.0;
  const double factor = coefficient * scale / denominator;
  // Divide by (-gamma + i).
  output[i] = {factor * (-gamma * value.x + value.y),
               factor * (-value.x - gamma * value.y)};
}

__global__ void addForcing(cufftDoubleComplex *rhs,
                           const cufftDoubleComplex *forcing,
                           std::size_t count) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < count)
    rhs[i] = rhs[i] + forcing[i];
}

struct DeviceStages {
  cufftDoubleComplex *n1{}, *n2{}, *n3{}, *n4{}, *a{}, *b{}, *c{};
};

__global__ void integrateStage(int stage, Integrator method, double h,
                               std::size_t count,
                               CoefficientPointers<cufftDoubleComplex> coeff,
                               cufftDoubleComplex *state, DeviceStages s,
                               const cufftDoubleComplex *noise) {
  const std::size_t i =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= count)
    return;
  if (stage == 0)
    s.a[i] = integrationStageA(method, h, i, coeff, state[i], s.n1[i]);
  else if (stage == 1)
    s.b[i] = integrationStageB(method, i, coeff, state[i], s.n1[i], s.n2[i]);
  else if (stage == 2)
    s.c[i] = integrationStageC(i, coeff, state[i], s.n1[i], s.n3[i]);
  else {
    const cufftDoubleComplex zero{0.0, 0.0};
    state[i] = integrationFinish(method, h, i, coeff, state[i], s.a[i], s.n1[i],
                                 s.n2[i], s.n3 ? s.n3[i] : zero,
                                 s.n4 ? s.n4[i] : zero);
    if (noise)
      state[i] = state[i] + noise[i];
  }
}

__global__ void suppressZeroMode(cufftDoubleComplex *state) {
  if (blockIdx.x == 0 && threadIdx.x == 0)
    state[0] = {0.0, 0.0};
}

class CudaBackend final : public NonlinearBackend {
public:
  explicit CudaBackend(const Parameters &p)
      : p_(p), baseCount_(p.nx * p.ny), paddedCount_(p.mx() * p.my()),
        input_(baseCount_), result_(baseCount_), paddedHat_(paddedCount_),
        psi_(paddedCount_), square_(paddedCount_), squareHat_(paddedCount_) {
    cufftCheck(cufftPlan2d(&plan_, static_cast<int>(p.my()),
                           static_cast<int>(p.mx()), CUFFT_Z2Z),
               "create dealiased transform plan");
  }

  ~CudaBackend() override {
    if (plan_)
      cufftDestroy(plan_);
  }

  bool deviceTimeStepping() const override { return true; }

  void initializeTimeStepping(const IntegrationCoefficients &coefficients,
                              const SpectralField &forcing,
                              const SpectralField &state) override {
    const auto hostFields = coefficients.fields();
    std::array<cufftDoubleComplex *, 10> pointers{};
    for (std::size_t i = 0; i < hostFields.size(); ++i) {
      coefficientBuffers_[i].allocate(hostFields[i]->size());
      if (!hostFields[i]->empty())
        cudaCheck(cudaMemcpy(coefficientBuffers_[i].data(),
                             hostFields[i]->data(),
                             hostFields[i]->size() * sizeof(Complex),
                             cudaMemcpyHostToDevice),
                  "upload integration coefficients");
      pointers[i] = coefficientBuffers_[i].data();
    }
    coefficients_ = {pointers[0], pointers[1], pointers[2], pointers[3],
                     pointers[4], pointers[5], pointers[6], pointers[7],
                     pointers[8], pointers[9]};
    for (const std::size_t index : {0UL, 1UL, 4UL})
      stageBuffers_[index].allocate(baseCount_);
    if (p_.integrator == Integrator::etd3 || p_.integrator == Integrator::etd4)
      for (const std::size_t index : {2UL, 5UL})
        stageBuffers_[index].allocate(baseCount_);
    if (p_.integrator == Integrator::etd4)
      for (const std::size_t index : {3UL, 6UL})
        stageBuffers_[index].allocate(baseCount_);
    stages_ = {stageBuffers_[0].data(), stageBuffers_[1].data(),
               stageBuffers_[2].data(), stageBuffers_[3].data(),
               stageBuffers_[4].data(), stageBuffers_[5].data(),
               stageBuffers_[6].data()};
    if (!forcing.empty()) {
      forcing_.allocate(forcing.size());
      cudaCheck(cudaMemcpy(forcing_.data(), forcing.data(),
                           forcing.size() * sizeof(Complex),
                           cudaMemcpyHostToDevice),
                "upload deterministic forcing");
    }
    if (p_.forcingEnabled && p_.forcingProfile != ForcingProfile::singleMode)
      noise_.allocate(baseCount_);
    uploadState(state);
  }

  void advance(const SpectralField &noise) override {
    if (!noise.empty()) {
      if (noise.size() != baseCount_ || !noise_.data())
        throw std::runtime_error("invalid device noise field");
      cudaCheck(cudaMemcpy(noise_.data(), noise.data(),
                           baseCount_ * sizeof(Complex),
                           cudaMemcpyHostToDevice),
                "upload stochastic increment");
    }
    rightHandSide(input_.data(), stages_.n1);
    launchStage(0);
    rightHandSide(stages_.a, stages_.n2);
    if (p_.integrator == Integrator::etd3 ||
        p_.integrator == Integrator::etd4) {
      launchStage(1);
      rightHandSide(stages_.b, stages_.n3);
    }
    if (p_.integrator == Integrator::etd4) {
      launchStage(2);
      rightHandSide(stages_.c, stages_.n4);
    }
    launchStage(3);
    if (p_.hypoviscosity > 0.0 && p_.hypoviscosityOrder < 0.0) {
      suppressZeroMode<<<1, 1>>>(input_.data());
      cudaCheck(cudaGetLastError(), "suppress device zero mode");
    }
  }

  void downloadState(SpectralField &state) override {
    state.resize(baseCount_);
    cudaCheck(cudaMemcpy(state.data(), input_.data(),
                         baseCount_ * sizeof(Complex), cudaMemcpyDeviceToHost),
              "download integrated wavefunction");
  }

  void evaluate(const SpectralField &state, SpectralField &output) override {
    uploadState(state);
    evaluateDevice(input_.data(), result_.data());
    output.resize(baseCount_);
    cudaCheck(cudaMemcpy(output.data(), result_.data(),
                         baseCount_ * sizeof(Complex), cudaMemcpyDeviceToHost),
              "download nonlinear term");
  }

private:
  static constexpr int threads = 256;
  int blocks(std::size_t count) const {
    return static_cast<int>((count + threads - 1) / threads);
  }

  void uploadState(const SpectralField &state) {
    if (state.size() != baseCount_)
      throw std::runtime_error("invalid nonlinear input size");
    cudaCheck(cudaMemcpy(input_.data(), state.data(),
                         baseCount_ * sizeof(Complex), cudaMemcpyHostToDevice),
              "upload wavefunction");
  }

  void launchStage(int stage) {
    integrateStage<<<blocks(baseCount_), threads>>>(
        stage, p_.integrator, p_.timeStep, baseCount_, coefficients_,
        input_.data(), stages_, noise_.data());
    cudaCheck(cudaGetLastError(), "integrate device Runge-Kutta stage");
  }

  void rightHandSide(const cufftDoubleComplex *state,
                     cufftDoubleComplex *output) {
    evaluateDevice(state, output);
    if (forcing_.data()) {
      addForcing<<<blocks(baseCount_), threads>>>(output, forcing_.data(),
                                                  baseCount_);
      cudaCheck(cudaGetLastError(), "add deterministic device forcing");
    }
  }

  void evaluateDevice(const cufftDoubleComplex *state,
                      cufftDoubleComplex *output) {
    embedBase<<<blocks(paddedCount_), threads>>>(
        state, paddedHat_.data(), p_.nx, p_.ny, p_.mx(), p_.my());
    cudaCheck(cudaGetLastError(), "embed device spectral field");
    cufftCheck(
        cufftExecZ2Z(plan_, paddedHat_.data(), psi_.data(), CUFFT_INVERSE),
        "execute wavefunction inverse transform");
    squareField<<<blocks(paddedCount_), threads>>>(psi_.data(), square_.data(),
                                                   paddedCount_);
    cudaCheck(cudaGetLastError(), "square device wavefunction");
    cufftCheck(
        cufftExecZ2Z(plan_, square_.data(), squareHat_.data(), CUFFT_FORWARD),
        "execute square forward transform");
    const double scale = 1.0 / static_cast<double>(paddedCount_);
    filterSquare<<<blocks(paddedCount_), threads>>>(
        squareHat_.data(), p_.nx, p_.ny, p_.mx(), p_.my(), scale);
    cudaCheck(cudaGetLastError(), "filter square convolution");
    cufftCheck(
        cufftExecZ2Z(plan_, squareHat_.data(), square_.data(), CUFFT_INVERSE),
        "execute filtered-square inverse transform");
    cubicField<<<blocks(paddedCount_), threads>>>(
        square_.data(), psi_.data(), paddedHat_.data(), paddedCount_);
    cudaCheck(cudaGetLastError(), "form cubic device nonlinearity");
    cufftCheck(cufftExecZ2Z(plan_, paddedHat_.data(), squareHat_.data(),
                            CUFFT_FORWARD),
               "execute nonlinear forward transform");
    extractBase<<<blocks(baseCount_), threads>>>(
        squareHat_.data(), output, p_.nx, p_.ny, p_.mx(), p_.my(), p_.lx(),
        p_.ly(), p_.nonlinearityCoefficient, p_.ginzburgLandauDamping,
        p_.ginzburgLandauCutoff, scale);
    cudaCheck(cudaGetLastError(), "extract device nonlinear modes");
  }

  Parameters p_;
  std::size_t baseCount_, paddedCount_;
  DeviceBuffer<cufftDoubleComplex> input_, result_, paddedHat_, psi_, square_,
      squareHat_;
  std::array<DeviceBuffer<cufftDoubleComplex>, 10> coefficientBuffers_;
  std::array<DeviceBuffer<cufftDoubleComplex>, 7> stageBuffers_;
  DeviceBuffer<cufftDoubleComplex> forcing_, noise_;
  CoefficientPointers<cufftDoubleComplex> coefficients_;
  DeviceStages stages_;
  cufftHandle plan_ = 0;
};
} // namespace

void backendInitialize(int &, char **&) {
  cudaCheck(cudaFree(nullptr), "initialize CUDA");
}
void backendFinalize() {}
void backendAbort(int) {}
void backendBarrier() {}
bool backendIsRoot() { return true; }
const char *backendName() { return "CUDA"; }
std::uint64_t backendSynchronizeSeed(std::uint64_t seed) { return seed; }
std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &p) {
  return std::make_unique<CudaBackend>(p);
}
