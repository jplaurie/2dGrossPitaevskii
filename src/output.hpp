#pragma once

#include "backend.hpp"
#include "fftw_utils.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct RestartState {
  double time = 0.0;
  std::uint64_t frame = 0;
  std::uint64_t randomSeed = 0;
  SpectralField wavefunction;
  std::string randomEngineState;
  std::string randomDistributionState;
  bool restarting = false;
};

struct DiagnosticsAverages {
  std::vector<double> waveActionSpectrum;
  std::vector<double> quadraticEnergySpectrum;
  std::vector<double> waveActionFlux;
  std::vector<double> fullEnergyFlux;
  std::size_t count = 0;
};

bool recoverOutputTransaction(const Parameters &parameters);
void beginOutputTransaction(const Parameters &parameters, std::uint64_t frame);
void finishOutputTransaction(const Parameters &parameters);
void writeRunRecords(const Parameters &parameters, const std::string &backend,
                     double time, std::uint64_t frame,
                     const std::vector<double> &amplitude,
                     std::size_t forcedModes,
                     double waveActionInjectionCoefficient,
                     double quadraticEnergyInjectionCoefficient);

RestartState readRestart(const Parameters &parameters, BaseTransform &transform,
                         bool isRoot);
void prepareOutputFiles(const Parameters &parameters, bool restarting,
                        std::uint64_t restartFrame);
void writeWavefunction(const Parameters &parameters, BaseTransform &transform,
                       const SpectralField &wavefunction, std::uint64_t frame);
double writeDiagnostics(const Parameters &parameters, BaseTransform &transform,
                        double time, std::uint64_t frame,
                        const SpectralField &wavefunction,
                        const SpectralField &nonlinear,
                        const std::vector<double> &forcingAmplitude,
                        const std::vector<double>
                            &stochasticQuarticInjectionWeight,
                        DiagnosticsAverages &averages);
void writeRestart(const Parameters &parameters, double time,
                  std::uint64_t frame, const SpectralField &wavefunction,
                  const std::string &randomEngineState,
                  const std::string &randomDistributionState);
void writeForcingFiles(const Parameters &parameters,
                       const std::vector<double> &amplitude,
                       std::size_t forcedModes,
                       double waveActionInjectionCoefficient,
                       double quadraticEnergyInjectionCoefficient);
