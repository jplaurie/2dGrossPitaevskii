#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

enum class Integrator { etd2, etd3, etd4, integratingFactorRk2 };
enum class ForcingProfile {
  annulus,
  gaussian,
  exponential,
  logNormal,
  singleMode
};

struct Parameters {
  std::size_t nx = 512;
  std::size_t ny = 512;
  double aspectRatio = 1.0;
  double timeStep = 1.0e-5;
  std::uint64_t numberOfSteps = 9'999'999'999ULL;
  std::uint64_t outputIntervalSteps = 100;

  Integrator integrator = Integrator::etd4;
  double dispersionCoefficient = -1.0;
  double nonlinearityCoefficient = 1000.0;
  double chemicalPotential = 0.0;

  double hyperviscosity = 1.0e-36;
  double hyperviscosityOrder = 8.0;
  bool hyperviscosityCutoffEnabled = false;
  double hyperviscosityCutoff = 1024.0;
  double hypoviscosity = 5.0e3;
  double hypoviscosityOrder = -2.0;
  bool hypoviscosityCutoffEnabled = false;
  double hypoviscosityCutoff = 4.0;
  double ginzburgLandauDamping = 0.0;
  double ginzburgLandauCutoff = 0.5;

  bool forcingEnabled = true;
  ForcingProfile forcingProfile = ForcingProfile::logNormal;
  double forcingWavenumber = 32.0;
  double forcingWidth = 2.0;
  double forcingAmplitude = 0.1;
  double forcingShapeOrder = 4.0;
  double forcingLogWidth = 0.05;
  double targetWaveActionInjectionRate = 0.0;
  std::uint64_t randomSeed = 0;

  bool writeModeDiagnostics = false;
  int threadCount = 0;
  bool overwriteOutput = false;
  std::filesystem::path initialConditionFile;
  std::filesystem::path dataDirectory = "data";
  std::filesystem::path outputDirectory = "output";

  [[nodiscard]] std::size_t mx() const { return 3 * nx / 2; }
  [[nodiscard]] std::size_t my() const { return 3 * ny / 2; }
  [[nodiscard]] double lx() const;
  [[nodiscard]] double ly() const;
  [[nodiscard]] std::size_t spectrumBins() const;
  [[nodiscard]] bool usesEtd() const {
    return integrator != Integrator::integratingFactorRk2;
  }
  [[nodiscard]] std::size_t nonlinearStageCount() const {
    return integrator == Integrator::etd4   ? 4
           : integrator == Integrator::etd3 ? 3
                                            : 2;
  }
};

[[nodiscard]] const char *integratorName(Integrator integrator);
[[nodiscard]] const char *forcingProfileName(ForcingProfile profile);
Parameters readParameters(const std::filesystem::path &path);
void validateParameters(const Parameters &parameters);
void writeParameterRecord(const Parameters &parameters,
                          const std::string &backend,
                          const std::filesystem::path &recordDirectory = {});
