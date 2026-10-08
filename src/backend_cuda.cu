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
__host__ __device__ inline cufftDoubleComplex operator*(double a, cufftDoubleComplex b) {
    return {a * b.x, a * b.y};
}
#ifdef GP2D_CUDA_MIXED
__host__ __device__ inline cufftComplex operator*(cufftComplex a, cufftComplex b) {
    return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
}
#endif

#include "backend.hpp"
#include "fftw_utils.hpp"
#include "spectral.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

static_assert(sizeof(Complex) == sizeof(cufftDoubleComplex));
static_assert(std::is_trivially_copyable_v<Complex>);

namespace {
#ifdef GP2D_CUDA_MIXED
using TransformComplex = cufftComplex;
using TransformReal = float;
constexpr cufftType transformType = CUFFT_C2C;
#else
using TransformComplex = cufftDoubleComplex;
using TransformReal = double;
constexpr cufftType transformType = CUFFT_Z2Z;
#endif

void cudaCheck(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
void cufftCheck(cufftResult status, const char *operation) {
    if (status != CUFFT_SUCCESS)
        throw std::runtime_error(std::string(operation) + " failed (cuFFT status " +
                                 std::to_string(static_cast<int>(status)) + ")");
}

bool environmentEnabled(const char *name) {
    const char *value = std::getenv(name);
    return value && value[0] != '\0' && std::string_view(value) != "0";
}

template <class T> class DeviceBuffer {
  public:
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t count) { allocate(count); }
    ~DeviceBuffer() { cudaFree(data_); }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    void allocate(std::size_t count) {
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::runtime_error("device buffer size overflow");
        cudaFree(data_);
        data_ = nullptr;
        count_ = count;
        if (count)
            cudaCheck(cudaMalloc(&data_, count * sizeof(T)), "allocate device buffer");
    }
    T *data() { return data_; }
    const T *data() const { return data_; }
    std::size_t size() const { return count_; }

  private:
    T *data_ = nullptr;
    std::size_t count_ = 0;
};

__device__ long deviceWave(std::size_t index, std::size_t count) {
    return index <= (count - 1) / 2 ? static_cast<long>(index)
                                    : static_cast<long>(index) - static_cast<long>(count);
}

__global__ void embedBaseSpectrum(const cufftDoubleComplex *input, TransformComplex *padded,
                                  std::size_t nx, std::size_t ny, std::size_t mx, std::size_t my) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= mx * my)
        return;
    const std::size_t px = i % mx, py = i / mx;
    const long kx = deviceWave(px, mx), ky = deviceWave(py, my);
    const bool retained = kx >= -static_cast<long>(nx / 2) && kx < static_cast<long>(nx / 2) &&
                          ky >= -static_cast<long>(ny / 2) && ky < static_cast<long>(ny / 2);
    if (!retained) {
        padded[i] = {0, 0};
        return;
    }
    const std::size_t x = kx >= 0 ? static_cast<std::size_t>(kx)
                                  : static_cast<std::size_t>(static_cast<long>(nx) + kx);
    const std::size_t y = ky >= 0 ? static_cast<std::size_t>(ky)
                                  : static_cast<std::size_t>(static_cast<long>(ny) + ky);
    const cufftDoubleComplex value = input[y * nx + x];
    padded[i] = {static_cast<TransformReal>(value.x), static_cast<TransformReal>(value.y)};
}

__global__ void squareField(const TransformComplex *input, TransformComplex *output,
                            std::size_t count) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        output[i] = input[i] * input[i];
}

__global__ void filterProjectedSquare(TransformComplex *field, std::size_t nx, std::size_t ny,
                                      std::size_t mx, std::size_t my, TransformReal scale) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= mx * my)
        return;
    const long kx = deviceWave(i % mx, mx), ky = deviceWave(i / mx, my);
    if (kx >= -static_cast<long>(nx / 2) && kx < static_cast<long>(nx / 2) &&
        ky >= -static_cast<long>(ny / 2) && ky < static_cast<long>(ny / 2)) {
        field[i].x *= scale;
        field[i].y *= scale;
    } else {
        field[i] = {0, 0};
    }
}

__global__ void formCubicProduct(const TransformComplex *square, const TransformComplex *psi,
                                 TransformComplex *output, std::size_t count) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    const TransformComplex conjugate{psi[i].x, -psi[i].y};
    output[i] = square[i] * conjugate;
}

__global__ void extractNonlinearTerm(const TransformComplex *padded, cufftDoubleComplex *output,
                                     std::size_t nx, std::size_t ny, std::size_t mx, std::size_t my,
                                     double lx, double ly, double coefficient, double damping,
                                     double dampingCutoff, double scale) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= nx * ny)
        return;
    const std::size_t x = i % nx, y = i / nx;
    const long kxIndex = deviceWave(x, nx), kyIndex = deviceWave(y, ny);
    const std::size_t px = kxIndex >= 0 ? static_cast<std::size_t>(kxIndex)
                                        : static_cast<std::size_t>(static_cast<long>(mx) + kxIndex);
    const std::size_t py = kyIndex >= 0 ? static_cast<std::size_t>(kyIndex)
                                        : static_cast<std::size_t>(static_cast<long>(my) + kyIndex);
    constexpr double twoPi = 6.283185307179586476925286766559;
    const double kx = twoPi * static_cast<double>(kxIndex) / lx;
    const double ky = twoPi * static_cast<double>(kyIndex) / ly;
    const double gamma = damping > 0.0 && sqrt(kx * kx + ky * ky) > dampingCutoff ? damping : 0.0;
    const TransformComplex value = padded[py * mx + px];
    const double denominator = gamma * gamma + 1.0;
    const double factor = coefficient * scale / denominator;
    // Divide by (-gamma + i).
    const double real = static_cast<double>(value.x);
    const double imaginary = static_cast<double>(value.y);
    output[i] = {factor * (-gamma * real + imaginary), factor * (-real - gamma * imaginary)};
}

__global__ void addDeterministicForcing(cufftDoubleComplex *rhs, const cufftDoubleComplex *forcing,
                                        std::size_t count) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        rhs[i] = rhs[i] + forcing[i];
}

__global__ void scatterStochasticNoise(cufftDoubleComplex *noise,
                                       const cufftDoubleComplex *compactNoise,
                                       const std::size_t *indices, std::size_t count) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        noise[indices[i]] = compactNoise[i];
}

struct DeviceStages {
    cufftDoubleComplex *nonlinearAtStart{}, *nonlinearAtStageA{}, *nonlinearAtStageB{},
        *nonlinearAtStageC{}, *stageA{}, *stageB{}, *stageC{};
};

__global__ void advanceIntegrationStage(int stage, Integrator method, double h, std::size_t count,
                                        CoefficientPointers<cufftDoubleComplex> coeff,
                                        cufftDoubleComplex *state, DeviceStages s,
                                        const cufftDoubleComplex *noise) {
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    if (stage == 0)
        s.stageA[i] = integrationStageA(method, h, i, coeff, state[i], s.nonlinearAtStart[i]);
    else if (stage == 1)
        s.stageB[i] = integrationStageB(method, i, coeff, state[i], s.nonlinearAtStart[i],
                                        s.nonlinearAtStageA[i]);
    else if (stage == 2)
        s.stageC[i] =
            integrationStageC(i, coeff, state[i], s.nonlinearAtStart[i], s.nonlinearAtStageB[i]);
    else {
        const cufftDoubleComplex zero{0.0, 0.0};
        state[i] = integrationFinish(method, h, i, coeff, state[i], s.stageA[i],
                                     s.nonlinearAtStart[i], s.nonlinearAtStageA[i],
                                     s.nonlinearAtStageB ? s.nonlinearAtStageB[i] : zero,
                                     s.nonlinearAtStageC ? s.nonlinearAtStageC[i] : zero);
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
    explicit CudaBackend(const Parameters &parameters)
        : parameters_(parameters), baseCount_(parameters.nx * parameters.ny),
          paddedCount_(parameters.mx() * parameters.my()), wavefunction_(baseCount_),
          nonlinearTerm_(baseCount_), paddedSpectrum_(paddedCount_), psi_(paddedCount_),
          projectedSquare_(paddedCount_), squareSpectrum_(paddedCount_),
          compactNoiseEnabled_(!environmentEnabled("GP2D_CUDA_FULL_NOISE")),
          profileEnabled_(environmentEnabled("GP2D_CUDA_PROFILE")),
          graphEnabled_(parameters.cudaGraphEnabled) {
        cudaCheck(cudaStreamCreate(&stream_), "create CUDA integration stream");
        cufftCheck(cufftPlan2d(&plan_, static_cast<int>(parameters.my()),
                               static_cast<int>(parameters.mx()), transformType),
                   "create dealiased transform plan");
        cufftCheck(cufftSetStream(plan_, stream_), "attach cuFFT plan to integration stream");
        if (profileEnabled_) {
            cudaCheck(cudaEventCreate(&profileStart_), "create CUDA profile event");
            cudaCheck(cudaEventCreate(&profileNoiseReady_), "create CUDA noise profile event");
            cudaCheck(cudaEventCreate(&profileEnd_), "create CUDA profile event");
        }
    }

    ~CudaBackend() override {
        if (profileEnabled_ && profileSteps_ > 0) {
            const double averageStep = profileStepMilliseconds_ / profileSteps_;
            const double averageNoise = profileNoiseMilliseconds_ / profileSteps_;
            const double percentage =
                profileStepMilliseconds_ > 0.0
                    ? 100.0 * profileNoiseMilliseconds_ / profileStepMilliseconds_
                    : 0.0;
            std::cout << "CUDA profile: noise_mode=" << (compactNoiseEnabled_ ? "compact" : "full")
                      << " steps=" << profileSteps_ << " average_step_ms=" << averageStep
                      << " average_noise_prepare_ms=" << averageNoise
                      << " noise_percentage=" << percentage << '\n';
        }
        if (profileStart_)
            cudaEventDestroy(profileStart_);
        if (profileNoiseReady_)
            cudaEventDestroy(profileNoiseReady_);
        if (profileEnd_)
            cudaEventDestroy(profileEnd_);
        if (graphExec_)
            cudaGraphExecDestroy(graphExec_);
        if (graph_)
            cudaGraphDestroy(graph_);
        if (plan_)
            cufftDestroy(plan_);
        if (stream_)
            cudaStreamDestroy(stream_);
    }

    NoiseLayout noiseLayout() const override {
        return compactNoiseEnabled_ ? NoiseLayout::forcedModes : NoiseLayout::fullField;
    }

    void initializeTimeStepping(const IntegrationCoefficients &coefficients,
                                const SpectralField &forcing, const SpectralField &state,
                                const std::vector<std::size_t> &noiseIndices) override {
        if (graphExec_) {
            cudaGraphExecDestroy(graphExec_);
            graphExec_ = nullptr;
        }
        if (graph_) {
            cudaGraphDestroy(graph_);
            graph_ = nullptr;
        }
        const auto hostFields = coefficients.fields();
        std::array<cufftDoubleComplex *, 10> pointers{};
        for (std::size_t i = 0; i < hostFields.size(); ++i) {
            coefficientBuffers_[i].allocate(hostFields[i]->size());
            if (!hostFields[i]->empty())
                cudaCheck(cudaMemcpy(coefficientBuffers_[i].data(), hostFields[i]->data(),
                                     hostFields[i]->size() * sizeof(Complex),
                                     cudaMemcpyHostToDevice),
                          "upload integration coefficients");
            pointers[i] = coefficientBuffers_[i].data();
        }
        coefficients_ = {pointers[0], pointers[1], pointers[2], pointers[3], pointers[4],
                         pointers[5], pointers[6], pointers[7], pointers[8], pointers[9]};
        const std::size_t stageCount = parameters_.nonlinearStageCount();
        for (std::size_t i = 0; i < stageCount; ++i)
            stageBuffers_[i].allocate(baseCount_);
        for (std::size_t i = 1; i < stageCount; ++i)
            stageBuffers_[i + 3].allocate(baseCount_);
        stages_ = {stageBuffers_[0].data(), stageBuffers_[1].data(), stageBuffers_[2].data(),
                   stageBuffers_[3].data(), stageBuffers_[4].data(), stageBuffers_[5].data(),
                   stageBuffers_[6].data()};
        if (!forcing.empty()) {
            forcing_.allocate(forcing.size());
            cudaCheck(cudaMemcpy(forcing_.data(), forcing.data(), forcing.size() * sizeof(Complex),
                                 cudaMemcpyHostToDevice),
                      "upload deterministic forcing");
        }
        if (parameters_.forcingEnabled &&
            parameters_.forcingProfile != ForcingProfile::singleMode) {
            noise_.allocate(baseCount_);
            if (compactNoiseEnabled_) {
                if (noiseIndices.empty())
                    throw std::runtime_error("compact stochastic noise has no modes");
                noiseIndices_.allocate(noiseIndices.size());
                compactNoise_.allocate(noiseIndices.size());
                cudaCheck(cudaMemcpy(noiseIndices_.data(), noiseIndices.data(),
                                     noiseIndices.size() * sizeof(std::size_t),
                                     cudaMemcpyHostToDevice),
                          "upload stochastic-noise indices");
                cudaCheck(cudaMemset(noise_.data(), 0, baseCount_ * sizeof(cufftDoubleComplex)),
                          "initialize device stochastic increment");
            }
        }
        uploadState(state);
        if (graphEnabled_) {
            // cuFFT may perform one-time setup on its first execution, which is
            // not legal during stream capture. Warm the plan without changing state.
            evaluateDevice(wavefunction_.data(), nonlinearTerm_.data());
            cudaCheck(cudaStreamSynchronize(stream_), "warm CUDA graph transform plan");
        }
    }

    void advanceTimeStep(const SpectralField &noise) override {
        if (profileEnabled_)
            cudaCheck(cudaEventRecord(profileStart_, stream_), "record CUDA profile start");
        if (!noise.empty()) {
            const std::size_t expectedSize =
                compactNoiseEnabled_ ? noiseIndices_.size() : baseCount_;
            if (noise.size() != expectedSize || !noise_.data())
                throw std::runtime_error("invalid device noise field");
            if (compactNoiseEnabled_) {
                cudaCheck(cudaMemcpyAsync(compactNoise_.data(), noise.data(),
                                          noise.size() * sizeof(Complex), cudaMemcpyHostToDevice,
                                          stream_),
                          "upload compact stochastic increment");
                scatterStochasticNoise<<<blocks(noise.size()), threads, 0, stream_>>>(
                    noise_.data(), compactNoise_.data(), noiseIndices_.data(), noise.size());
                cudaCheck(cudaGetLastError(), "scatter device stochastic increment");
            } else {
                cudaCheck(cudaMemcpyAsync(noise_.data(), noise.data(), baseCount_ * sizeof(Complex),
                                          cudaMemcpyHostToDevice, stream_),
                          "upload stochastic increment");
            }
        }
        if (profileEnabled_)
            cudaCheck(cudaEventRecord(profileNoiseReady_, stream_),
                      "record CUDA noise profile event");
        if (!noise.empty())
            // The host noise buffer is reused immediately by the solver. Ensure
            // the asynchronous copy/scatter has consumed it before returning.
            cudaCheck(cudaStreamSynchronize(stream_), "prepare CUDA stochastic increment");
        if (graphEnabled_) {
            if (!graphExec_)
                captureStepGraph();
            cudaCheck(cudaGraphLaunch(graphExec_, stream_), "launch CUDA timestep graph");
        } else {
            executeStep();
        }
        if (profileEnabled_) {
            cudaCheck(cudaEventRecord(profileEnd_, stream_), "record CUDA profile end");
            cudaCheck(cudaEventSynchronize(profileEnd_), "synchronize CUDA profile event");
            float stepMilliseconds = 0.0F, noiseMilliseconds = 0.0F;
            cudaCheck(cudaEventElapsedTime(&stepMilliseconds, profileStart_, profileEnd_),
                      "measure CUDA step time");
            cudaCheck(cudaEventElapsedTime(&noiseMilliseconds, profileStart_, profileNoiseReady_),
                      "measure CUDA noise time");
            profileStepMilliseconds_ += stepMilliseconds;
            profileNoiseMilliseconds_ += noiseMilliseconds;
            ++profileSteps_;
        }
    }

    void executeStep() {
        rightHandSide(wavefunction_.data(), stages_.nonlinearAtStart);
        launchStage(0);
        rightHandSide(stages_.stageA, stages_.nonlinearAtStageA);
        if (parameters_.nonlinearStageCount() >= 3) {
            launchStage(1);
            rightHandSide(stages_.stageB, stages_.nonlinearAtStageB);
        }
        if (parameters_.nonlinearStageCount() == 4) {
            launchStage(2);
            rightHandSide(stages_.stageC, stages_.nonlinearAtStageC);
        }
        launchStage(3);
        if (parameters_.hypoviscosity > 0.0 && parameters_.hypoviscosityOrder < 0.0) {
            suppressZeroMode<<<1, 1, 0, stream_>>>(wavefunction_.data());
            cudaCheck(cudaGetLastError(), "suppress device zero mode");
        }
    }

    void downloadState(SpectralField &state) override {
        cudaCheck(cudaStreamSynchronize(stream_), "synchronize integrated wavefunction");
        state.resize(baseCount_);
        cudaCheck(cudaMemcpy(state.data(), wavefunction_.data(), baseCount_ * sizeof(Complex),
                             cudaMemcpyDeviceToHost),
                  "download integrated wavefunction");
    }

    void evaluate(const SpectralField &state, SpectralField &output) override {
        uploadState(state);
        evaluateDevice(wavefunction_.data(), nonlinearTerm_.data());
        cudaCheck(cudaStreamSynchronize(stream_), "synchronize nonlinear evaluation");
        output.resize(baseCount_);
        cudaCheck(cudaMemcpy(output.data(), nonlinearTerm_.data(), baseCount_ * sizeof(Complex),
                             cudaMemcpyDeviceToHost),
                  "download nonlinear term");
    }

    void evaluateCurrent(SpectralField &output) override {
        evaluateDevice(wavefunction_.data(), nonlinearTerm_.data());
        cudaCheck(cudaStreamSynchronize(stream_), "synchronize current nonlinear evaluation");
        output.resize(baseCount_);
        cudaCheck(cudaMemcpy(output.data(), nonlinearTerm_.data(), baseCount_ * sizeof(Complex),
                             cudaMemcpyDeviceToHost),
                  "download current nonlinear term");
    }

  private:
    static constexpr int threads = 256;
    int blocks(std::size_t count) const {
        return static_cast<int>((count + threads - 1) / threads);
    }

    void uploadState(const SpectralField &state) {
        if (state.size() != baseCount_)
            throw std::runtime_error("invalid nonlinear input size");
        cudaCheck(cudaMemcpy(wavefunction_.data(), state.data(), baseCount_ * sizeof(Complex),
                             cudaMemcpyHostToDevice),
                  "upload wavefunction");
    }

    void executeTransform(TransformComplex *input, TransformComplex *output, int direction,
                          const char *operation) {
#ifdef GP2D_CUDA_MIXED
        cufftCheck(cufftExecC2C(plan_, input, output, direction), operation);
#else
        cufftCheck(cufftExecZ2Z(plan_, input, output, direction), operation);
#endif
    }

    void launchStage(int stage) {
        advanceIntegrationStage<<<blocks(baseCount_), threads, 0, stream_>>>(
            stage, parameters_.integrator, parameters_.timeStep, baseCount_, coefficients_,
            wavefunction_.data(), stages_, noise_.data());
        cudaCheck(cudaGetLastError(), "integrate device Runge-Kutta stage");
    }

    void rightHandSide(const cufftDoubleComplex *state, cufftDoubleComplex *output) {
        evaluateDevice(state, output);
        if (forcing_.data()) {
            addDeterministicForcing<<<blocks(baseCount_), threads, 0, stream_>>>(
                output, forcing_.data(), baseCount_);
            cudaCheck(cudaGetLastError(), "add deterministic device forcing");
        }
    }

    void evaluateDevice(const cufftDoubleComplex *state, cufftDoubleComplex *output) {
        embedBaseSpectrum<<<blocks(paddedCount_), threads, 0, stream_>>>(
            state, paddedSpectrum_.data(), parameters_.nx, parameters_.ny, parameters_.mx(),
            parameters_.my());
        cudaCheck(cudaGetLastError(), "embed device spectral field");
        executeTransform(paddedSpectrum_.data(), psi_.data(), CUFFT_INVERSE,
                         "execute wavefunction inverse transform");
        squareField<<<blocks(paddedCount_), threads, 0, stream_>>>(
            psi_.data(), projectedSquare_.data(), paddedCount_);
        cudaCheck(cudaGetLastError(), "square device wavefunction");
        executeTransform(projectedSquare_.data(), squareSpectrum_.data(), CUFFT_FORWARD,
                         "execute square forward transform");
        const double scale = 1.0 / static_cast<double>(paddedCount_);
        filterProjectedSquare<<<blocks(paddedCount_), threads, 0, stream_>>>(
            squareSpectrum_.data(), parameters_.nx, parameters_.ny, parameters_.mx(),
            parameters_.my(), static_cast<TransformReal>(scale));
        cudaCheck(cudaGetLastError(), "filter square convolution");
        executeTransform(squareSpectrum_.data(), projectedSquare_.data(), CUFFT_INVERSE,
                         "execute filtered-square inverse transform");
        formCubicProduct<<<blocks(paddedCount_), threads, 0, stream_>>>(
            projectedSquare_.data(), psi_.data(), paddedSpectrum_.data(), paddedCount_);
        cudaCheck(cudaGetLastError(), "form cubic device nonlinearity");
        executeTransform(paddedSpectrum_.data(), squareSpectrum_.data(), CUFFT_FORWARD,
                         "execute nonlinear forward transform");
        extractNonlinearTerm<<<blocks(baseCount_), threads, 0, stream_>>>(
            squareSpectrum_.data(), output, parameters_.nx, parameters_.ny, parameters_.mx(),
            parameters_.my(), parameters_.lx(), parameters_.ly(),
            parameters_.nonlinearityCoefficient, parameters_.ginzburgLandauDamping,
            parameters_.ginzburgLandauCutoff, scale);
        cudaCheck(cudaGetLastError(), "extract device nonlinear modes");
    }

    void captureStepGraph() {
        cudaCheck(cudaStreamBeginCapture(stream_, cudaStreamCaptureModeThreadLocal),
                  "begin CUDA timestep graph capture");
        executeStep();
        cudaCheck(cudaStreamEndCapture(stream_, &graph_), "end CUDA timestep graph capture");
        cudaCheck(cudaGraphInstantiate(&graphExec_, graph_, nullptr, nullptr, 0),
                  "instantiate CUDA timestep graph");
    }

    Parameters parameters_;
    std::size_t baseCount_, paddedCount_;
    DeviceBuffer<cufftDoubleComplex> wavefunction_, nonlinearTerm_;
    DeviceBuffer<TransformComplex> paddedSpectrum_, psi_, projectedSquare_, squareSpectrum_;
    std::array<DeviceBuffer<cufftDoubleComplex>, 10> coefficientBuffers_;
    std::array<DeviceBuffer<cufftDoubleComplex>, 7> stageBuffers_;
    DeviceBuffer<cufftDoubleComplex> forcing_, noise_, compactNoise_;
    DeviceBuffer<std::size_t> noiseIndices_;
    CoefficientPointers<cufftDoubleComplex> coefficients_;
    DeviceStages stages_;
    bool compactNoiseEnabled_ = true, profileEnabled_ = false, graphEnabled_ = false;
    cudaEvent_t profileStart_ = nullptr, profileNoiseReady_ = nullptr, profileEnd_ = nullptr;
    cudaStream_t stream_ = nullptr;
    cudaGraph_t graph_ = nullptr;
    cudaGraphExec_t graphExec_ = nullptr;
    std::uint64_t profileSteps_ = 0;
    double profileStepMilliseconds_ = 0.0, profileNoiseMilliseconds_ = 0.0;
    cufftHandle plan_ = 0;
};
} // namespace

void backendInitialize(int &, char **&) { cudaCheck(cudaFree(nullptr), "initialize CUDA"); }
void backendFinalize() { saveFftwWisdom(); }
void backendAbort(int) {}
void backendBarrier() { cudaCheck(cudaDeviceSynchronize(), "synchronize CUDA device"); }
bool backendIsRoot() { return true; }
const char *backendName() {
#ifdef GP2D_CUDA_MIXED
    return "CUDA mixed (FP64 state / FP32 FFT)";
#else
    return "CUDA";
#endif
}
std::uint64_t backendSynchronizeSeed(std::uint64_t seed) { return seed; }
std::unique_ptr<NonlinearBackend> makeBackend(const Parameters &parameters) {
    configureFftw(parameters);
    return std::make_unique<CudaBackend>(parameters);
}
