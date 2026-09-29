#pragma once

#include "backend.hpp"
#include "fftw_utils.hpp"
#include "output.hpp"

#include <array>
#include <memory>
#include <random>
#include <vector>

class Solver {
  public:
    Solver(Parameters parameters, std::unique_ptr<NonlinearBackend> backend);
    void run();

  private:
    void buildLinearOperator();
    void buildIntegrationCoefficients();
    void buildForcing();
    void generateNoise(SpectralField &noise);
    void rightHandSide(const SpectralField &input, SpectralField &output);
    void step(SpectralField &wavefunction);
    RestartState prepareRun();
    void validateRunBounds(const RestartState &state) const;
    void restoreRandomState(const RestartState &state);
    void initializeDeviceTimeStepping(const SpectralField &wavefunction);
    void writeState(const RestartState &state);
    double writeOutputFrame(const RestartState &state, DiagnosticsAverages &averages);

    Parameters parameters_;
    std::unique_ptr<NonlinearBackend> backend_;
    BaseTransform baseTransform_;
    SpectralField linearOperator_;
    IntegrationCoefficients coefficients_;
    SpectralField noise_;
    std::array<SpectralField, 4> nonlinearStages_;
    std::array<SpectralField, 3> stageStates_;
    SpectralField diagnosticNonlinearTerm_, deterministicForcing_;
    std::vector<double> forcingAmplitude_, stochasticNoiseScale_;
    std::vector<double> stochasticQuarticInjectionWeight_;
    std::vector<std::size_t> forcedIndices_;
    std::size_t forcedModeCount_ = 0;
    double waveActionInjectionCoefficient_ = 0.0;
    double quadraticEnergyInjectionCoefficient_ = 0.0;
    std::mt19937_64 random_;
    std::normal_distribution<double> normal_{0.0, 1.0};
};
