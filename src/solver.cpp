#include "solver.hpp"
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
[[maybe_unused]] constexpr std::size_t parallelThreshold = 16384;

SpectralField field(const Parameters &p, bool needed = true) {
  return needed ? SpectralField(p.nx * p.ny) : SpectralField{};
}

template <class Operation>
void forEachIndex(std::size_t count, Operation operation) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count >= parallelThreshold)
#endif
  for (std::ptrdiff_t raw = 0; raw < static_cast<std::ptrdiff_t>(count); ++raw)
    operation(static_cast<std::size_t>(raw));
}

void requireFinite(const SpectralField &values, const char *description) {
  if (!std::all_of(values.begin(), values.end(), [](Complex value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
      }))
    throw std::runtime_error(std::string(description) +
                             " contains a non-finite coefficient");
}

std::uint64_t resolveRandomSeed(std::uint64_t seed) {
  if (seed == 0)
    seed = static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
  return backendSynchronizeSeed(seed);
}
} // namespace

Solver::Solver(Parameters p, std::unique_ptr<NonlinearBackend> backend)
    : p_(std::move(p)), backend_(std::move(backend)), baseTransform_(p_),
      linear_(field(p_)),
      noise_(field(p_, p_.forcingEnabled &&
                           p_.forcingProfile != ForcingProfile::singleMode)),
      n1_(field(p_, !backend_->deviceTimeStepping())),
      n2_(field(p_, !backend_->deviceTimeStepping())),
      n3_(field(p_, !backend_->deviceTimeStepping() &&
                        (p_.integrator == Integrator::etd3 ||
                         p_.integrator == Integrator::etd4))),
      n4_(field(p_, !backend_->deviceTimeStepping() &&
                        p_.integrator == Integrator::etd4)),
      stageA_(field(p_, !backend_->deviceTimeStepping())),
      stageB_(field(p_, !backend_->deviceTimeStepping() &&
                            (p_.integrator == Integrator::etd3 ||
                             p_.integrator == Integrator::etd4))),
      stageC_(field(p_, !backend_->deviceTimeStepping() &&
                            p_.integrator == Integrator::etd4)),
      diagnosticNonlinear_(field(p_)),
      deterministicForcing_(
          field(p_, p_.forcingEnabled &&
                        p_.forcingProfile == ForcingProfile::singleMode)),
      forcingAmplitude_(p_.nx * p_.ny), noiseScale_(noise_.size()),
      random_(p_.randomSeed = resolveRandomSeed(p_.randomSeed)) {
#ifdef _OPENMP
  if (p_.threadCount > 0)
    omp_set_num_threads(p_.threadCount);
#endif
  std::filesystem::create_directories(p_.dataDirectory);
  std::filesystem::create_directories(p_.outputDirectory);
  buildLinearOperator();
  buildIntegrationCoefficients();
  SpectralField{}.swap(linear_);
  if (p_.forcingEnabled)
    buildForcing();
}

void Solver::buildLinearOperator() {
  forEachIndex(linear_.size(), [&](std::size_t index) {
    const std::size_t x = index % p_.nx;
    const std::size_t y = index / p_.nx;
    const double k2 = waveNumberSquared(p_, x, y);
    const double k = std::sqrt(k2);
    const double gamma = ginzburgLandauFactor(p_, k2);
    const double hamiltonian =
        -p_.dispersionCoefficient * k2 + p_.chemicalPotential;
    Complex value = hamiltonian / Complex(-gamma, 1.0);
    if (p_.hyperviscosity > 0.0 &&
        (!p_.hyperviscosityCutoffEnabled || k > p_.hyperviscosityCutoff))
      value -= p_.hyperviscosity * std::pow(k2, p_.hyperviscosityOrder);
    if ((k2 > 0.0 || p_.hypoviscosityOrder >= 0.0) && p_.hypoviscosity > 0.0 &&
        (!p_.hypoviscosityCutoffEnabled || k < p_.hypoviscosityCutoff))
      value -= p_.hypoviscosity * std::pow(k2, p_.hypoviscosityOrder);
    linear_[index] = value;
  });
  if (p_.hypoviscosity > 0.0 && p_.hypoviscosityOrder < 0.0)
    linear_[0] = Complex{};
  requireFinite(linear_,
                "linear operator; reduce dissipation coefficients or orders");
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
  auto &c = coefficients_;
  c.e1 = field(p_);
  if (p_.usesEtd()) {
    c.q1 = field(p_);
    c.f1 = field(p_);
    if (p_.integrator != Integrator::etd2) {
      c.e2 = field(p_);
      c.q2 = field(p_);
      c.f2 = field(p_);
      c.f3 = field(p_);
    }
    if (p_.integrator == Integrator::etd4) {
      c.q3 = field(p_);
      c.q4 = field(p_);
      c.q5 = field(p_);
    }
  }
  const double h = p_.timeStep;
  forEachIndex(linear_.size(), [&](std::size_t i) {
    const Complex z = h * linear_[i];
    c.e1[i] = std::exp(z);
    if (p_.usesEtd()) {
      const Complex p1 = phi(z, 1), p2 = phi(z, 2);
      if (p_.integrator == Integrator::etd2) {
        c.q1[i] = h * p1;
        c.f1[i] = h * p2;
      } else {
        const Complex p3 = phi(z, 3), halfP1 = phi(0.5 * z, 1);
        c.e2[i] = std::exp(0.5 * z);
        c.q1[i] = (0.5 * h) * halfP1;
        c.q2[i] = h * p1;
        c.f1[i] = h * (p1 - 3.0 * p2 + 4.0 * p3);
        c.f2[i] = h * (p2 - 2.0 * p3);
        c.f3[i] = h * (-p2 + 4.0 * p3);
        if (p_.integrator == Integrator::etd4) {
          const Complex halfP2 = phi(0.5 * z, 2);
          c.q2[i] = h * (0.5 * halfP1 - halfP2);
          c.q3[i] = h * halfP2;
          c.q4[i] = h * (p1 - 2.0 * p2);
          c.q5[i] = (2.0 * h) * p2;
        }
      }
    }
    if (!noiseScale_.empty()) {
      const double a = linear_[i].real();
      const double ah = a * h;
      const double variance = std::abs(ah) < 1.e-8
                                  ? h * (1.0 + ah)
                                  : std::expm1(2.0 * ah) / (2.0 * a);
      noiseScale_[i] = std::sqrt(variance);
    }
  });
  for (const auto *values : c.fields())
    requireFinite(*values, "time integration coefficients");
  if (!noiseScale_.empty() &&
      !std::all_of(noiseScale_.begin(), noiseScale_.end(),
                   [](double x) { return std::isfinite(x) && x >= 0.0; }))
    throw std::runtime_error(
        "stochastic integration coefficients are non-finite");
}

void Solver::buildForcing() {
  const long mode = static_cast<long>(p_.forcingWavenumber);
  const double maximumLog = std::log(std::numeric_limits<double>::max());
  for (std::size_t y = 0; y < p_.ny; ++y) {
    const long kyIndex = signedWave(y, p_.ny);
    for (std::size_t x = 0; x < p_.nx; ++x) {
      const long kxIndex = signedWave(x, p_.nx);
      const double k2 = waveNumberSquared(p_, x, y);
      const double k = std::sqrt(k2);
      double amplitude = 0.0;
      if (k > 0.0 && p_.forcingProfile == ForcingProfile::annulus &&
          std::abs(k - p_.forcingWavenumber) < p_.forcingWidth)
        amplitude = p_.forcingAmplitude;
      else if (k > 0.0 && p_.forcingProfile == ForcingProfile::gaussian)
        amplitude = p_.forcingAmplitude *
                    std::exp(-0.5 * std::pow((k - p_.forcingWavenumber) /
                                                 p_.forcingWidth,
                                             2.0));
      else if (k > 0.0 && p_.forcingProfile == ForcingProfile::exponential) {
        const double logRatio =
            p_.forcingShapeOrder * std::log(k / p_.forcingWavenumber);
        if (std::isfinite(logRatio) && logRatio <= maximumLog) {
          const double ratio = std::exp(logRatio);
          amplitude = p_.forcingAmplitude * std::exp(logRatio - ratio);
        }
      } else if (k > 0.0 && p_.forcingProfile == ForcingProfile::logNormal)
        amplitude =
            p_.forcingAmplitude *
            std::exp(-0.5 * std::pow(std::log(k / p_.forcingWavenumber) /
                                         p_.forcingLogWidth,
                                     2.0));
      else if (p_.forcingProfile == ForcingProfile::singleMode &&
               std::abs(kxIndex) == mode && std::abs(kyIndex) == mode)
        amplitude = p_.forcingAmplitude;
      const std::size_t index = spectralIndex(x, y, p_.nx);
      forcingAmplitude_[index] = amplitude;
      if (amplitude != 0.0) {
        forcedIndices_.push_back(index);
        ++forcedModeCount_;
        waveActionInjectionCoefficient_ += amplitude * amplitude;
        quadraticEnergyInjectionCoefficient_ +=
            (-p_.dispersionCoefficient * k2 + p_.chemicalPotential) *
            amplitude * amplitude;
      }
    }
  }
  if (forcedModeCount_ == 0)
    throw std::runtime_error(
        "forcingEnabled is true, but the selected profile contains no modes");
  if (p_.targetWaveActionInjectionRate > 0.0) {
    const double scale = std::sqrt(p_.targetWaveActionInjectionRate /
                                   waveActionInjectionCoefficient_);
    for (double &value : forcingAmplitude_)
      value *= scale;
    waveActionInjectionCoefficient_ *= scale * scale;
    quadraticEnergyInjectionCoefficient_ *= scale * scale;
  }
  if (!deterministicForcing_.empty())
    forEachIndex(deterministicForcing_.size(), [&](std::size_t i) {
      deterministicForcing_[i] = forcingAmplitude_[i];
    });
  if (backendIsRoot())
    std::cout << "number of forced modes = " << forcedModeCount_
              << "\nwave-action injection coefficient = "
              << waveActionInjectionCoefficient_
              << "\nquadratic-energy injection coefficient = "
              << quadraticEnergyInjectionCoefficient_ << '\n';
}

void Solver::generateNoise(SpectralField &noise) {
  std::fill(noise.begin(), noise.end(), Complex{});
  constexpr double circularScale = 0.7071067811865475244;
  for (const std::size_t i : forcedIndices_)
    noise[i] = forcingAmplitude_[i] * circularScale * noiseScale_[i] *
               Complex(normal_(random_), normal_(random_));
}

void Solver::rightHandSide(const SpectralField &input, SpectralField &output) {
  backend_->evaluate(input, output);
  if (p_.forcingEnabled && p_.forcingProfile == ForcingProfile::singleMode)
    forEachIndex(output.size(),
                 [&](std::size_t i) { output[i] += deterministicForcing_[i]; });
}

void Solver::step(SpectralField &w) {
  if (!noise_.empty())
    generateNoise(noise_);
  if (backend_->deviceTimeStepping()) {
    backend_->advance(noise_);
    return;
  }
  const auto c = coefficients_.pointers();
  rightHandSide(w, n1_);
  forEachIndex(w.size(), [&](std::size_t i) {
    stageA_[i] =
        integrationStageA(p_.integrator, p_.timeStep, i, c, w[i], n1_[i]);
  });
  rightHandSide(stageA_, n2_);
  if (!n3_.empty()) {
    forEachIndex(w.size(), [&](std::size_t i) {
      stageB_[i] = integrationStageB(p_.integrator, i, c, w[i], n1_[i], n2_[i]);
    });
    rightHandSide(stageB_, n3_);
  }
  if (!n4_.empty()) {
    forEachIndex(w.size(), [&](std::size_t i) {
      stageC_[i] = integrationStageC(i, c, w[i], n1_[i], n3_[i]);
    });
    rightHandSide(stageC_, n4_);
  }
  forEachIndex(w.size(), [&](std::size_t i) {
    w[i] = integrationFinish(p_.integrator, p_.timeStep, i, c, w[i], stageA_[i],
                             n1_[i], n2_[i], n3_.empty() ? Complex{} : n3_[i],
                             n4_.empty() ? Complex{} : n4_[i]);
    if (!noise_.empty())
      w[i] += noise_[i];
  });
  enforceStateConstraints(w, p_);
}

void Solver::run() {
  const bool recoveredFresh = backendIsRoot() && recoverOutputTransaction(p_);
  backendBarrier();
  RestartState state = readRestart(p_, baseTransform_, backendIsRoot());
  const long double finalTime =
      static_cast<long double>(state.time) +
      static_cast<long double>(p_.timeStep) * p_.numberOfSteps;
  if (finalTime > static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::runtime_error("requested run would overflow simulation time");
  const std::uint64_t scheduledOutputs =
      p_.numberOfSteps / p_.outputIntervalSteps +
      (p_.numberOfSteps % p_.outputIntervalSteps != 0 ? 1 : 0);
  if (state.frame >
      std::numeric_limits<std::uint64_t>::max() - scheduledOutputs)
    throw std::runtime_error("requested run would overflow output frame count");
  if (!state.randomEngineState.empty()) {
    std::istringstream saved(state.randomEngineState);
    if (!(saved >> random_))
      throw std::runtime_error("cannot restore random-generator state");
    p_.randomSeed = state.randomSeed;
  }
  if (!state.randomDistributionState.empty()) {
    std::istringstream saved(state.randomDistributionState);
    if (!(saved >> normal_))
      throw std::runtime_error("cannot restore normal-distribution state");
  } else if (state.restarting) {
    normal_.reset();
  }
  if (backendIsRoot()) {
    prepareOutputFiles(p_, state.restarting || recoveredFresh, state.frame);
    writeRunRecords(p_, backendName(), state.time, state.frame,
                    forcingAmplitude_, forcedModeCount_,
                    waveActionInjectionCoefficient_,
                    quadraticEnergyInjectionCoefficient_);
    if (!state.restarting) {
      beginOutputTransaction(p_, state.frame);
      writeWavefunction(p_, baseTransform_, state.wavefunction, state.frame);
      std::ostringstream randomState, distributionState;
      randomState << random_;
      distributionState << normal_;
      writeRestart(p_, state.time, state.frame, state.wavefunction,
                   randomState.str(), distributionState.str());
      finishOutputTransaction(p_);
    }
    std::cout << "backend = " << backendName() << "\nnx = " << p_.nx
              << " ny = " << p_.ny << " timeStep = " << p_.timeStep
              << " numberOfSteps = " << p_.numberOfSteps
              << " outputIntervalSteps = " << p_.outputIntervalSteps
              << "\nrandom seed = " << p_.randomSeed << '\n';
  }
  if (backend_->deviceTimeStepping()) {
    backend_->initializeTimeStepping(coefficients_, deterministicForcing_,
                                     state.wavefunction);
    coefficients_ = {};
  }
  DiagnosticsAverages averages;
  const auto start = std::chrono::steady_clock::now();
  for (std::uint64_t stepIndex = 0; stepIndex < p_.numberOfSteps; ++stepIndex) {
    const std::uint64_t stepNumber = stepIndex + 1;
    step(state.wavefunction);
    state.time += p_.timeStep;
    if (stepNumber % p_.outputIntervalSteps == 0 ||
        stepNumber == p_.numberOfSteps) {
      ++state.frame;
      if (backend_->deviceTimeStepping())
        backend_->downloadState(state.wavefunction);
      backend_->evaluate(state.wavefunction, diagnosticNonlinear_);
      if (backendIsRoot()) {
        beginOutputTransaction(p_, state.frame);
        const double energy = writeDiagnostics(p_, baseTransform_, state.time,
                                               state.frame, state.wavefunction,
                                               diagnosticNonlinear_, averages);
        writeWavefunction(p_, baseTransform_, state.wavefunction, state.frame);
        std::ostringstream randomState, distributionState;
        randomState << random_;
        distributionState << normal_;
        writeRestart(p_, state.time, state.frame, state.wavefunction,
                     randomState.str(), distributionState.str());
        finishOutputTransaction(p_);
        std::cout << "time = " << state.time << " file = " << state.frame
                  << " Energy = " << energy << '\n';
      }
    }
  }
  if (backendIsRoot()) {
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - start;
    std::cout << "time taken for code is = " << elapsed.count() << '\n';
  }
}
