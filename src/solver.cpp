#include "solver.hpp"
#include "parallel.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
SpectralField makeSpectralField(const Parameters &parameters, bool needed = true) {
    return needed ? SpectralField(parameters.nx * parameters.ny) : SpectralField{};
}

void requireFinite(const SpectralField &values, const char *description) {
    if (!std::all_of(values.begin(), values.end(), [](Complex value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag());
        }))
        throw std::runtime_error(std::string(description) + " contains a non-finite coefficient");
}

std::uint64_t resolveRandomSeed(std::uint64_t configuredSeed) {
    if (configuredSeed == 0)
        configuredSeed = static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
    return backendSynchronizeSeed(configuredSeed);
}
} // namespace

Solver::Solver(Parameters parameters, std::unique_ptr<NonlinearBackend> backend)
    : parameters_(std::move(parameters)), backend_(std::move(backend)), baseTransform_(parameters_),
      linearOperator_(makeSpectralField(parameters_)),
      noise_(makeSpectralField(parameters_,
                               parameters_.forcingEnabled &&
                                   parameters_.forcingProfile != ForcingProfile::singleMode)),
      diagnosticNonlinearTerm_(makeSpectralField(parameters_)),
      deterministicForcing_(makeSpectralField(parameters_, parameters_.forcingEnabled &&
                                                               parameters_.forcingProfile ==
                                                                   ForcingProfile::singleMode)),
      forcingAmplitude_(parameters_.nx * parameters_.ny), stochasticNoiseScale_(noise_.size()),
      random_(parameters_.randomSeed = resolveRandomSeed(parameters_.randomSeed)) {
#ifdef _OPENMP
    if (parameters_.threadCount > 0)
        omp_set_num_threads(parameters_.threadCount);
#endif
    std::filesystem::create_directories(parameters_.dataDirectory);
    std::filesystem::create_directories(parameters_.outputDirectory);
    if (!backend_->supportsDeviceTimeStepping()) {
        const std::size_t count = parameters_.nonlinearStageCount();
        for (std::size_t i = 0; i < count; ++i)
            nonlinearStages_[i] = makeSpectralField(parameters_);
        for (std::size_t i = 1; i < count; ++i)
            stageStates_[i - 1] = makeSpectralField(parameters_);
    }
    buildLinearOperator();
    buildIntegrationCoefficients();
    SpectralField{}.swap(linearOperator_);
    if (parameters_.forcingEnabled) {
        buildForcing();
        if (backend_->compactStochasticNoise() &&
            parameters_.forcingProfile != ForcingProfile::singleMode)
            noise_.resize(forcedIndices_.size());
    }
}

void Solver::buildLinearOperator() {
    forEachIndex(linearOperator_.size(), [&](std::size_t index) {
        const std::size_t x = index % parameters_.nx;
        const std::size_t y = index / parameters_.nx;
        const double k2 = waveNumberSquared(parameters_, x, y);
        const double k = std::sqrt(k2);
        const double gamma = ginzburgLandauFactor(parameters_, k2);
        const double hamiltonian =
            -parameters_.dispersionCoefficient * k2 + parameters_.chemicalPotential;
        Complex value = hamiltonian / Complex(-gamma, 1.0);
        if (parameters_.hyperviscosity > 0.0 &&
            (!parameters_.hyperviscosityCutoffEnabled || k > parameters_.hyperviscosityCutoff))
            value -= parameters_.hyperviscosity * std::pow(k2, parameters_.hyperviscosityOrder);
        if ((k2 > 0.0 || parameters_.hypoviscosityOrder >= 0.0) &&
            parameters_.hypoviscosity > 0.0 &&
            (!parameters_.hypoviscosityCutoffEnabled || k < parameters_.hypoviscosityCutoff))
            value -= parameters_.hypoviscosity * std::pow(k2, parameters_.hypoviscosityOrder);
        linearOperator_[index] = value;
    });
    if (parameters_.hypoviscosity > 0.0 && parameters_.hypoviscosityOrder < 0.0)
        linearOperator_[0] = Complex{};
    requireFinite(linearOperator_, "linear operator; reduce dissipation coefficients or orders");
}

void Solver::buildIntegrationCoefficients() {
    const auto phi = [](Complex z, int order) {
        double factorial = 1.0;
        for (int k = 2; k <= order; ++k)
            factorial *= k;
        if (std::abs(z) < 2.0) {
            Complex term = 1.0 / factorial;
            Complex sum = term;
            for (int k = 1; k < 64; ++k) {
                term *= z / static_cast<double>(order + k);
                sum += term;
                if (std::abs(term) < 1.e-17 * std::abs(sum))
                    break;
            }
            return sum;
        }
        Complex value = (std::exp(z) - 1.0) / z;
        double previousFactorial = 1.0;
        for (int k = 2; k <= order; ++k) {
            value = (value - 1.0 / previousFactorial) / z;
            previousFactorial *= k;
        }
        return value;
    };
    IntegrationCoefficients &coefficients = coefficients_;
    coefficients.e1 = makeSpectralField(parameters_);
    if (parameters_.usesEtd()) {
        coefficients.q1 = makeSpectralField(parameters_);
        coefficients.f1 = makeSpectralField(parameters_);
        if (parameters_.integrator != Integrator::etd2) {
            coefficients.e2 = makeSpectralField(parameters_);
            coefficients.q2 = makeSpectralField(parameters_);
            coefficients.f2 = makeSpectralField(parameters_);
            coefficients.f3 = makeSpectralField(parameters_);
        }
        if (parameters_.integrator == Integrator::etd4) {
            coefficients.q3 = makeSpectralField(parameters_);
            coefficients.q4 = makeSpectralField(parameters_);
            coefficients.q5 = makeSpectralField(parameters_);
        }
    }
    const double h = parameters_.timeStep;
    forEachIndex(linearOperator_.size(), [&](std::size_t i) {
        const Complex z = h * linearOperator_[i];
        coefficients.e1[i] = std::exp(z);
        if (parameters_.usesEtd()) {
            const Complex p1 = phi(z, 1), p2 = phi(z, 2);
            if (parameters_.integrator == Integrator::etd2) {
                coefficients.q1[i] = h * p1;
                coefficients.f1[i] = h * p2;
            } else {
                const Complex p3 = phi(z, 3), halfP1 = phi(0.5 * z, 1);
                coefficients.e2[i] = std::exp(0.5 * z);
                coefficients.q1[i] = (0.5 * h) * halfP1;
                coefficients.q2[i] = h * p1;
                coefficients.f1[i] = h * (p1 - 3.0 * p2 + 4.0 * p3);
                coefficients.f2[i] = h * (p2 - 2.0 * p3);
                coefficients.f3[i] = h * (-p2 + 4.0 * p3);
                if (parameters_.integrator == Integrator::etd4) {
                    const Complex halfP2 = phi(0.5 * z, 2);
                    coefficients.q2[i] = h * (0.5 * halfP1 - halfP2);
                    coefficients.q3[i] = h * halfP2;
                    coefficients.q4[i] = h * (p1 - 2.0 * p2);
                    coefficients.q5[i] = (2.0 * h) * p2;
                }
            }
        }
        if (!stochasticNoiseScale_.empty()) {
            const double a = linearOperator_[i].real();
            const double ah = a * h;
            const double variance =
                std::abs(ah) < 1.e-8 ? h * (1.0 + ah) : std::expm1(2.0 * ah) / (2.0 * a);
            stochasticNoiseScale_[i] = std::sqrt(variance);
        }
    });
    for (const auto *values : coefficients.fields())
        requireFinite(*values, "time integration coefficients");
    if (!stochasticNoiseScale_.empty() &&
        !std::all_of(stochasticNoiseScale_.begin(), stochasticNoiseScale_.end(),
                     [](double x) { return std::isfinite(x) && x >= 0.0; }))
        throw std::runtime_error("stochastic integration coefficients are non-finite");
}

void Solver::buildForcing() {
    const long mode = static_cast<long>(parameters_.forcingWavenumber);
    const double maximumLog = std::log(std::numeric_limits<double>::max());
    for (std::size_t y = 0; y < parameters_.ny; ++y) {
        const long kyIndex = signedWave(y, parameters_.ny);
        for (std::size_t x = 0; x < parameters_.nx; ++x) {
            const long kxIndex = signedWave(x, parameters_.nx);
            const double k2 = waveNumberSquared(parameters_, x, y);
            const double k = std::sqrt(k2);
            double amplitude = 0.0;
            if (k > 0.0 && parameters_.forcingProfile == ForcingProfile::annulus &&
                std::abs(k - parameters_.forcingWavenumber) < parameters_.forcingWidth)
                amplitude = parameters_.forcingAmplitude;
            else if (k > 0.0 && parameters_.forcingProfile == ForcingProfile::gaussian)
                amplitude = parameters_.forcingAmplitude *
                            std::exp(-0.5 * std::pow((k - parameters_.forcingWavenumber) /
                                                         parameters_.forcingWidth,
                                                     2.0));
            else if (k > 0.0 && parameters_.forcingProfile == ForcingProfile::exponential) {
                const double logRatio =
                    parameters_.forcingShapeOrder * std::log(k / parameters_.forcingWavenumber);
                if (std::isfinite(logRatio) && logRatio <= maximumLog) {
                    const double ratio = std::exp(logRatio);
                    amplitude = parameters_.forcingAmplitude * std::exp(logRatio - ratio);
                }
            } else if (k > 0.0 && parameters_.forcingProfile == ForcingProfile::logNormal)
                amplitude = parameters_.forcingAmplitude *
                            std::exp(-0.5 * std::pow(std::log(k / parameters_.forcingWavenumber) /
                                                         parameters_.forcingLogWidth,
                                                     2.0));
            else if (parameters_.forcingProfile == ForcingProfile::singleMode &&
                     std::abs(kxIndex) == mode && std::abs(kyIndex) == mode)
                amplitude = parameters_.forcingAmplitude;
            const std::size_t index = spectralIndex(x, y, parameters_.nx);
            forcingAmplitude_[index] = amplitude;
            if (amplitude != 0.0) {
                forcedIndices_.push_back(index);
                ++forcedModeCount_;
                waveActionInjectionCoefficient_ += amplitude * amplitude;
                quadraticEnergyInjectionCoefficient_ +=
                    (-parameters_.dispersionCoefficient * k2 + parameters_.chemicalPotential) *
                    amplitude * amplitude;
            }
        }
    }
    if (forcedModeCount_ == 0)
        throw std::runtime_error(
            "forcingEnabled is true, but the selected profile contains no modes");
    if (parameters_.targetWaveActionInjectionRate > 0.0) {
        const double scale =
            std::sqrt(parameters_.targetWaveActionInjectionRate / waveActionInjectionCoefficient_);
        for (double &value : forcingAmplitude_)
            value *= scale;
        waveActionInjectionCoefficient_ *= scale * scale;
        quadraticEnergyInjectionCoefficient_ *= scale * scale;
    }
    if (parameters_.forcingProfile != ForcingProfile::singleMode &&
        parameters_.nonlinearityCoefficient != 0.0) {
        const long minimumX = -static_cast<long>(parameters_.nx / 2);
        const long maximumX = static_cast<long>((parameters_.nx - 1) / 2);
        const long minimumY = -static_cast<long>(parameters_.ny / 2);
        const long maximumY = static_cast<long>((parameters_.ny - 1) / 2);
        const std::size_t stride = parameters_.nx + 1;
        // A summed-area table evaluates the accepted forcing power for every
        // shifted spectral overlap without a nested convolution over all modes.
        std::vector<double> prefix((parameters_.ny + 1) * stride, 0.0);
        for (std::size_t y = 0; y < parameters_.ny; ++y)
            for (std::size_t x = 0; x < parameters_.nx; ++x) {
                const std::size_t centeredX =
                    static_cast<std::size_t>(signedWave(x, parameters_.nx) - minimumX);
                const std::size_t centeredY =
                    static_cast<std::size_t>(signedWave(y, parameters_.ny) - minimumY);
                const double amplitude = forcingAmplitude_[spectralIndex(x, y, parameters_.nx)];
                prefix[(centeredY + 1) * stride + centeredX + 1] = amplitude * amplitude;
            }
        for (std::size_t y = 0; y < parameters_.ny; ++y)
            for (std::size_t x = 0; x < parameters_.nx; ++x) {
                const std::size_t i = (y + 1) * stride + x + 1;
                prefix[i] += prefix[i - 1] + prefix[i - stride] - prefix[i - stride - 1];
            }
        stochasticQuarticInjectionWeight_.resize(parameters_.nx * parameters_.ny);
        for (std::size_t y = 0; y < parameters_.ny; ++y)
            for (std::size_t x = 0; x < parameters_.nx; ++x) {
                const long waveX = signedWave(x, parameters_.nx);
                const long waveY = signedWave(y, parameters_.ny);
                const long firstX = std::max(minimumX, minimumX - waveX);
                const long lastX = std::min(maximumX, maximumX - waveX);
                const long firstY = std::max(minimumY, minimumY - waveY);
                const long lastY = std::min(maximumY, maximumY - waveY);
                const std::size_t x0 = static_cast<std::size_t>(firstX - minimumX);
                const std::size_t x1 = static_cast<std::size_t>(lastX - minimumX + 1);
                const std::size_t y0 = static_cast<std::size_t>(firstY - minimumY);
                const std::size_t y1 = static_cast<std::size_t>(lastY - minimumY + 1);
                const double acceptedForcingPower =
                    prefix[y1 * stride + x1] - prefix[y0 * stride + x1] - prefix[y1 * stride + x0] +
                    prefix[y0 * stride + x0];
                stochasticQuarticInjectionWeight_[spectralIndex(x, y, parameters_.nx)] =
                    2.0 * parameters_.nonlinearityCoefficient * acceptedForcingPower;
            }
    }
    if (!deterministicForcing_.empty())
        forEachIndex(deterministicForcing_.size(),
                     [&](std::size_t i) { deterministicForcing_[i] = forcingAmplitude_[i]; });
    if (backendIsRoot())
        std::cout << "number of forced modes = " << forcedModeCount_
                  << "\nwave-action injection coefficient = " << waveActionInjectionCoefficient_
                  << "\nquadratic-energy injection coefficient = "
                  << quadraticEnergyInjectionCoefficient_ << '\n';
}

void Solver::generateNoise(SpectralField &noise) {
    const bool compact = backend_->compactStochasticNoise();
    if (compact) {
        if (noise.size() != forcedIndices_.size())
            throw std::logic_error("invalid compact stochastic-noise buffer");
    } else {
        std::fill(noise.begin(), noise.end(), Complex{});
    }
    constexpr double circularScale = 0.7071067811865475244;
    for (std::size_t position = 0; position < forcedIndices_.size(); ++position) {
        const std::size_t i = forcedIndices_[position];
        const Complex value = forcingAmplitude_[i] * circularScale * stochasticNoiseScale_[i] *
                              Complex(normal_(random_), normal_(random_));
        noise[compact ? position : i] = value;
    }
}

void Solver::rightHandSide(const SpectralField &input, SpectralField &output) {
    backend_->evaluate(input, output);
    if (parameters_.forcingEnabled && parameters_.forcingProfile == ForcingProfile::singleMode)
        forEachIndex(output.size(), [&](std::size_t i) { output[i] += deterministicForcing_[i]; });
}

void Solver::step(SpectralField &wavefunction) {
    if (!noise_.empty())
        generateNoise(noise_);
    if (backend_->supportsDeviceTimeStepping()) {
        backend_->advanceDeviceState(noise_);
        return;
    }
    const auto coefficientPointers = coefficients_.pointers();
    SpectralField &nonlinearAtStart = nonlinearStages_[0];
    SpectralField &nonlinearAtStageA = nonlinearStages_[1];
    SpectralField &nonlinearAtStageB = nonlinearStages_[2];
    SpectralField &nonlinearAtStageC = nonlinearStages_[3];
    SpectralField &stageA = stageStates_[0];
    SpectralField &stageB = stageStates_[1];
    SpectralField &stageC = stageStates_[2];

    rightHandSide(wavefunction, nonlinearAtStart);
    forEachIndex(wavefunction.size(), [&](std::size_t i) {
        stageA[i] = integrationStageA(parameters_.integrator, parameters_.timeStep, i,
                                      coefficientPointers, wavefunction[i], nonlinearAtStart[i]);
    });
    rightHandSide(stageA, nonlinearAtStageA);
    if (!nonlinearAtStageB.empty()) {
        forEachIndex(wavefunction.size(), [&](std::size_t i) {
            stageB[i] =
                integrationStageB(parameters_.integrator, i, coefficientPointers, wavefunction[i],
                                  nonlinearAtStart[i], nonlinearAtStageA[i]);
        });
        rightHandSide(stageB, nonlinearAtStageB);
    }
    if (!nonlinearAtStageC.empty()) {
        forEachIndex(wavefunction.size(), [&](std::size_t i) {
            stageC[i] = integrationStageC(i, coefficientPointers, wavefunction[i],
                                          nonlinearAtStart[i], nonlinearAtStageB[i]);
        });
        rightHandSide(stageC, nonlinearAtStageC);
    }
    forEachIndex(wavefunction.size(), [&](std::size_t i) {
        wavefunction[i] =
            integrationFinish(parameters_.integrator, parameters_.timeStep, i, coefficientPointers,
                              wavefunction[i], stageA[i], nonlinearAtStart[i], nonlinearAtStageA[i],
                              nonlinearAtStageB.empty() ? Complex{} : nonlinearAtStageB[i],
                              nonlinearAtStageC.empty() ? Complex{} : nonlinearAtStageC[i]);
        if (!noise_.empty())
            wavefunction[i] += noise_[i];
    });
    enforceStateConstraints(wavefunction, parameters_);
}

void Solver::writeState(const RestartState &state) {
    writeWavefunction(parameters_, baseTransform_, state.wavefunction, state.frame);
    std::ostringstream randomState, distributionState;
    randomState << random_;
    distributionState << normal_;
    writeRestart(parameters_, state.time, state.frame, state.wavefunction, randomState.str(),
                 distributionState.str());
}

void Solver::validateRunBounds(const RestartState &state) const {
    const long double finalTime =
        static_cast<long double>(state.time) +
        static_cast<long double>(parameters_.timeStep) * parameters_.numberOfSteps;
    if (finalTime > static_cast<long double>(std::numeric_limits<double>::max()))
        throw std::runtime_error("requested run would overflow simulation time");

    const std::uint64_t scheduledOutputs =
        parameters_.numberOfSteps / parameters_.outputIntervalSteps +
        (parameters_.numberOfSteps % parameters_.outputIntervalSteps != 0 ? 1 : 0);
    if (state.frame > std::numeric_limits<std::uint64_t>::max() - scheduledOutputs)
        throw std::runtime_error("requested run would overflow output frame count");
}

void Solver::restoreRandomState(const RestartState &state) {
    if (!state.randomEngineState.empty()) {
        std::istringstream savedEngine(state.randomEngineState);
        if (!(savedEngine >> random_))
            throw std::runtime_error("cannot restore random-generator state");
        parameters_.randomSeed = state.randomSeed;
    }
    if (!state.randomDistributionState.empty()) {
        std::istringstream savedDistribution(state.randomDistributionState);
        if (!(savedDistribution >> normal_))
            throw std::runtime_error("cannot restore normal-distribution state");
    } else if (state.restarting) {
        normal_.reset();
    }
}

void Solver::initializeDeviceTimeStepping(const SpectralField &wavefunction) {
    if (!backend_->supportsDeviceTimeStepping())
        return;
    backend_->initializeDeviceState(coefficients_, deterministicForcing_, wavefunction,
                                    forcedIndices_);
    coefficients_ = {};
}

RestartState Solver::prepareRun() {
    const bool recoveredFresh = backendIsRoot() && recoverOutputTransaction(parameters_);
    backendBarrier();
    RestartState state = readRestart(parameters_, baseTransform_, backendIsRoot());
    validateRunBounds(state);
    restoreRandomState(state);

    if (backendIsRoot()) {
        prepareOutputFiles(parameters_, state.restarting || recoveredFresh, state.frame);
        writeRunRecords(parameters_, backendName(), state.time, state.frame, forcingAmplitude_,
                        forcedModeCount_, waveActionInjectionCoefficient_,
                        quadraticEnergyInjectionCoefficient_);
        if (!state.restarting) {
            beginOutputTransaction(parameters_, state.frame);
            writeState(state);
            finishOutputTransaction(parameters_);
        }
        std::cout << "backend = " << backendName() << "\nnx = " << parameters_.nx
                  << " ny = " << parameters_.ny << " timeStep = " << parameters_.timeStep
                  << " numberOfSteps = " << parameters_.numberOfSteps
                  << " outputIntervalSteps = " << parameters_.outputIntervalSteps
                  << "\nrandom seed = " << parameters_.randomSeed << '\n';
    }
    initializeDeviceTimeStepping(state.wavefunction);
    return state;
}

double Solver::writeOutputFrame(const RestartState &state, DiagnosticsAverages &averages) {
    beginOutputTransaction(parameters_, state.frame);
    const double energy = writeDiagnostics(
        parameters_, baseTransform_, state.time, state.frame, state.wavefunction,
        diagnosticNonlinearTerm_, forcingAmplitude_, stochasticQuarticInjectionWeight_, averages);
    writeState(state);
    finishOutputTransaction(parameters_);
    return energy;
}

void Solver::run() {
    RestartState state = prepareRun();
    DiagnosticsAverages averages;
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t stepIndex = 0; stepIndex < parameters_.numberOfSteps; ++stepIndex) {
        const std::uint64_t stepNumber = stepIndex + 1;
        step(state.wavefunction);
        state.time += parameters_.timeStep;
        if (stepNumber % parameters_.outputIntervalSteps == 0 ||
            stepNumber == parameters_.numberOfSteps) {
            ++state.frame;
            if (backend_->supportsDeviceTimeStepping())
                backend_->downloadDeviceState(state.wavefunction);
            backend_->evaluate(state.wavefunction, diagnosticNonlinearTerm_);
            if (backendIsRoot()) {
                const double energy = writeOutputFrame(state, averages);
                std::cout << "time = " << state.time << " file = " << state.frame
                          << " Energy = " << energy << '\n';
            }
        }
    }
    if (backendIsRoot()) {
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        std::cout << "time taken for code is = " << elapsed.count() << '\n';
    }
}
