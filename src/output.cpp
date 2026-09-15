#include "output.hpp"
#include "io_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
constexpr std::size_t checkpointChunkComplexValues = 4096;

bool finiteField(const SpectralField &field) {
  return std::all_of(field.begin(), field.end(), [](Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  });
}

std::filesystem::path wavefunctionPath(const Parameters &p,
                                       std::uint64_t frame) {
  std::ostringstream name;
  name << "wavefunction_" << std::setw(8) << std::setfill('0') << frame
       << ".dat";
  return p.dataDirectory / name.str();
}

std::filesystem::path checkpointPath(const Parameters &p, std::uint64_t frame) {
  std::ostringstream name;
  name << "checkpoint_" << std::setw(8) << std::setfill('0') << frame << ".bin";
  return p.dataDirectory / name.str();
}

std::ofstream numericOutput(const std::filesystem::path &path,
                            std::ios::openmode mode = std::ios::out) {
  std::ofstream out(path, mode);
  if (!out)
    throw std::runtime_error("cannot write output file: " + path.string());
  out << std::scientific << std::setprecision(12);
  return out;
}

SpectralField readCheckpoint(const Parameters &p, std::uint64_t frame) {
  const auto path = checkpointPath(p, frame);
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("cannot open spectral checkpoint: " +
                             path.string());
  char magic[8]{};
  std::uint64_t nx = 0, ny = 0, count = 0;
  input.read(magic, sizeof(magic));
  input.read(reinterpret_cast<char *>(&nx), sizeof(nx));
  input.read(reinterpret_cast<char *>(&ny), sizeof(ny));
  input.read(reinterpret_cast<char *>(&count), sizeof(count));
  if (!input || std::string(magic, 7) != "GP2DCP1" || nx != p.nx ||
      ny != p.ny || count != p.nx * p.ny)
    throw std::runtime_error("invalid spectral checkpoint header: " +
                             path.string());
  SpectralField field(static_cast<std::size_t>(count));
  std::vector<double> buffer(2 * checkpointChunkComplexValues);
  for (std::size_t offset = 0; offset < field.size();) {
    const std::size_t chunk =
        std::min(checkpointChunkComplexValues, field.size() - offset);
    input.read(reinterpret_cast<char *>(buffer.data()),
               static_cast<std::streamsize>(2 * chunk * sizeof(double)));
    if (!input)
      break;
    for (std::size_t i = 0; i < chunk; ++i)
      field[offset + i] = Complex(buffer[2 * i], buffer[2 * i + 1]);
    offset += chunk;
  }
  if (!input || input.peek() != std::ifstream::traits_type::eof())
    throw std::runtime_error("invalid spectral checkpoint payload: " +
                             path.string());
  if (!finiteField(field))
    throw std::runtime_error("spectral checkpoint contains non-finite values");
  return field;
}

std::vector<Complex> readPhysicalField(const std::filesystem::path &path,
                                       const Parameters &p) {
  const std::size_t count = p.nx * p.ny;
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open wavefunction field: " +
                             path.string());
  std::vector<double> values;
  double value = 0.0;
  while (input >> value) {
    if (!std::isfinite(value))
      throw std::runtime_error("wavefunction contains a non-finite value: " +
                               path.string());
    values.push_back(value);
  }
  if (!input.eof())
    throw std::runtime_error("invalid value in wavefunction field: " +
                             path.string());
  if (values.size() != 2 * count)
    throw std::runtime_error(
        "wavefunction field must contain two values per grid point: " +
        path.string());
  std::vector<Complex> field(count);
  for (std::size_t i = 0; i < count; ++i)
    field[i] = Complex(values[2 * i], values[2 * i + 1]);
  return field;
}

void initializeCsv(const std::filesystem::path &path, const std::string &header,
                   bool append, bool overwrite) {
  if (append && std::filesystem::exists(path)) {
    std::ifstream input(path);
    std::string existingHeader;
    std::getline(input, existingHeader);
    if (existingHeader != header)
      throw std::runtime_error(
          "CSV header does not match this solver version: " + path.string());
    return;
  }
  if (!append && std::filesystem::exists(path) && !overwrite)
    throw std::runtime_error("refusing to overwrite existing output: " +
                             path.string());
  std::ofstream output(path, std::ios::trunc);
  if (!output)
    throw std::runtime_error("cannot initialize CSV file: " + path.string());
  output << header << '\n';
}

std::uint64_t lastCsvFrame(const std::filesystem::path &path) {
  std::ifstream input(path);
  std::string line, last;
  while (std::getline(input, line))
    if (!line.empty())
      last = line;
  if (last.empty() || last.starts_with("time,frame,"))
    return 0;
  const auto first = last.find(','), second = last.find(',', first + 1);
  if (first == std::string::npos || second == std::string::npos)
    throw std::runtime_error("malformed CSV: " + path.string());
  return std::stoull(last.substr(first + 1, second - first - 1));
}

bool containsRunData(const std::filesystem::path &directory) {
  if (!std::filesystem::exists(directory))
    return false;
  for (const auto &entry : std::filesystem::directory_iterator(directory)) {
    const std::string name = entry.path().filename().string();
    if (name == "restart_state.txt" || name.starts_with("wavefunction_") ||
        name.starts_with("checkpoint_"))
      return true;
  }
  return false;
}

bool containsSolverOutput(const std::filesystem::path &directory) {
  constexpr std::array names{"diagnostics.csv",
                             "spectra.csv",
                             "fluxes.csv",
                             "modes.csv",
                             "forcing_summary.csv",
                             "forcing_spectrum.csv",
                             "resolved_parameters.txt"};
  return std::any_of(names.begin(), names.end(), [&](const char *name) {
    return std::filesystem::exists(directory / name);
  });
}

} // namespace

RestartState readRestart(const Parameters &p, BaseTransform &transform,
                         bool isRoot) {
  std::filesystem::create_directories(p.dataDirectory);
  std::filesystem::create_directories(p.outputDirectory);
  RestartState state;
  std::vector<Complex> physical;
  const auto restartPath = p.dataDirectory / "restart_state.txt";
  if (std::filesystem::exists(restartPath)) {
    std::ifstream input(restartPath);
    std::string format, key;
    std::getline(input, format);
    if (format != "gp2d_restart_v1")
      throw std::runtime_error("unsupported restart format: " +
                               restartPath.string());
    std::size_t savedNx = 0, savedNy = 0;
    double savedAspect = 0.0;
    if (!(input >> key >> state.time) || key != "time" ||
        !(input >> key >> state.frame) || key != "frame" ||
        !(input >> key >> savedNx) || key != "nx" ||
        !(input >> key >> savedNy) || key != "ny" ||
        !(input >> key >> savedAspect) || key != "aspectRatio" ||
        !(input >> key >> state.randomSeed) || key != "randomSeed")
      throw std::runtime_error("malformed restart metadata: " +
                               restartPath.string());
    input >> std::ws;
    std::getline(input, key, ' ');
    if (key != "randomEngine")
      throw std::runtime_error("restart is missing randomEngine state");
    std::getline(input, state.randomEngineState);
    std::getline(input, key, ' ');
    if (key != "randomDistribution")
      throw std::runtime_error("restart is missing randomDistribution state");
    std::getline(input, state.randomDistributionState);
    if (!std::isfinite(state.time) || state.time < 0.0 || savedNx != p.nx ||
        savedNy != p.ny || savedAspect != p.aspectRatio)
      throw std::runtime_error(
          "restart grid, aspect ratio, or time is invalid");
    state.wavefunction = readCheckpoint(p, state.frame);
    state.restarting = true;
    if (isRoot)
      std::cout << "restarting at time " << state.time << " from frame "
                << state.frame << '\n';
  } else {
    if (p.initialConditionFile.empty()) {
      physical.assign(p.nx * p.ny, Complex{});
      if (isRoot)
        std::cout << "starting from zero wavefunction\n";
    } else {
      physical = readPhysicalField(p.initialConditionFile, p);
      if (isRoot)
        std::cout << "starting from " << p.initialConditionFile << '\n';
    }
  }
  if (state.wavefunction.empty()) {
    transform.forward(physical, state.wavefunction);
    enforceStateConstraints(state.wavefunction, p);
  }
  if (!finiteField(state.wavefunction))
    throw std::runtime_error(
        "initial wavefunction contains non-finite coefficients");
  return state;
}

void prepareOutputFiles(const Parameters &p, bool restarting,
                        std::uint64_t restartFrame) {
  constexpr std::array csvNames{"diagnostics.csv", "spectra.csv", "fluxes.csv",
                                "modes.csv"};
  if (!restarting && !p.overwriteOutput &&
      (containsRunData(p.dataDirectory) ||
       containsSolverOutput(p.outputDirectory)))
    throw std::runtime_error("solver output already exists; use new "
                             "directories or set overwriteOutput true");
  if (restarting) {
    for (const char *name : csvNames) {
      const auto path = p.outputDirectory / name;
      if (std::filesystem::exists(path) && lastCsvFrame(path) > restartFrame)
        throw std::runtime_error(path.string() +
                                 " contains frames newer than the restart");
    }
  }
  initializeCsv(
      p.outputDirectory / "diagnostics.csv",
      "time,frame,total_energy,kinetic_energy,potential_energy,nonlinear_"
      "energy,"
      "wave_action,wave_action_dissipation_hypo,wave_action_dissipation_hyper,"
      "quadratic_energy_dissipation_hypo,quadratic_energy_dissipation_hyper,"
      "total_energy_dissipation_hypo,total_energy_dissipation_hyper,"
      "expected_full_energy_injection",
      restarting, p.overwriteOutput);
  initializeCsv(
      p.outputDirectory / "spectra.csv",
      "time,frame,wavenumber,wave_action_spectrum,quadratic_energy_spectrum,"
      "segment_mean_wave_action_spectrum,segment_mean_quadratic_energy_"
      "spectrum",
      restarting, p.overwriteOutput);
  initializeCsv(p.outputDirectory / "fluxes.csv",
                "time,frame,wavenumber,wave_action_flux,full_energy_flux,"
                "segment_mean_wave_action_flux,segment_mean_full_energy_flux",
                restarting, p.overwriteOutput);
  if (p.writeModeDiagnostics)
    initializeCsv(p.outputDirectory / "modes.csv",
                  "time,frame,psi_1_0_real,psi_1_0_imag,psi_0_1_real,"
                  "psi_0_1_imag,psi_1_1_real,psi_1_1_imag,psi_2_1_real,"
                  "psi_2_1_imag,psi_0_3_real,psi_0_3_imag",
                  restarting, p.overwriteOutput);
}

void writeWavefunction(const Parameters &p, BaseTransform &transform,
                       const SpectralField &wavefunction, std::uint64_t frame) {
  const auto path = wavefunctionPath(p, frame);
  if (std::filesystem::exists(path) && !p.overwriteOutput)
    throw std::runtime_error("refusing to overwrite wavefunction snapshot: " +
                             path.string());
  std::vector<Complex> physical;
  transform.inverse(wavefunction, physical);
  const auto temporary = std::filesystem::path(path.string() + ".tmp");
  auto out = numericOutput(temporary);
  for (std::size_t y = 0; y < p.ny; ++y) {
    for (std::size_t x = 0; x < p.nx; ++x) {
      const Complex value = physical[spectralIndex(x, y, p.nx)];
      out << value.real() << ' ' << value.imag();
      out << (x + 1 == p.nx ? '\n' : ' ');
    }
  }
  closeChecked(out, "failed while writing wavefunction snapshot");
  std::filesystem::rename(temporary, path);
}

double writeDiagnostics(const Parameters &p, BaseTransform &transform,
                        double time, std::uint64_t frame,
                        const SpectralField &w, const SpectralField &nonlinear,
                        const std::vector<double> &forcingAmplitude,
                        const std::vector<double>
                            &stochasticQuarticInjectionWeight,
                        DiagnosticsAverages &avg) {
  if (forcingAmplitude.size() != w.size() ||
      (!stochasticQuarticInjectionWeight.empty() &&
       stochasticQuarticInjectionWeight.size() != w.size()))
    throw std::runtime_error("invalid forcing diagnostics size");
  const std::size_t bins = p.spectrumBins();
  if (avg.waveActionSpectrum.empty()) {
    avg.waveActionSpectrum.assign(bins, 0.0);
    avg.quadraticEnergySpectrum.assign(bins, 0.0);
    avg.waveActionFlux.assign(bins, 0.0);
    avg.fullEnergyFlux.assign(bins, 0.0);
  }
  std::vector<double> waveSpectrum(bins), quadraticSpectrum(bins),
      waveShell(bins), quadraticShell(bins);
  SpectralField hamiltonianRate(w.size());
  const double area = p.lx() * p.ly();
  const double binWidth = std::min(2.0 * gpPi / p.lx(), 2.0 * gpPi / p.ly());
  double waveAction = 0.0, kineticEnergy = 0.0, potentialEnergy = 0.0;
  double waveHypo = 0.0, waveHyper = 0.0;
  double quadraticHypo = 0.0, quadraticHyper = 0.0;
  double totalEnergyHypo = 0.0, totalEnergyHyper = 0.0;
  double expectedFullEnergyInjection = 0.0;
  const bool stochasticForcing =
      p.forcingEnabled && p.forcingProfile != ForcingProfile::singleMode;
  const bool deterministicForcing =
      p.forcingEnabled && p.forcingProfile == ForcingProfile::singleMode;
  for (std::size_t y = 0; y < p.ny; ++y) {
    for (std::size_t x = 0; x < p.nx; ++x) {
      const std::size_t index = spectralIndex(x, y, p.nx);
      const double k2 = waveNumberSquared(p, x, y), k = std::sqrt(k2);
      const double n = std::norm(w[index]);
      const double kineticWeight = -p.dispersionCoefficient * k2;
      const double weight = kineticWeight + p.chemicalPotential;
      waveAction += area * n;
      kineticEnergy += area * kineticWeight * n;
      potentialEnergy += area * p.chemicalPotential * n;
      double hyper = 0.0, hypo = 0.0;
      if (p.hyperviscosity > 0.0 &&
          (!p.hyperviscosityCutoffEnabled || k > p.hyperviscosityCutoff))
        hyper = p.hyperviscosity * std::pow(k2, p.hyperviscosityOrder);
      if ((k2 > 0.0 || p.hypoviscosityOrder >= 0.0) && p.hypoviscosity > 0.0 &&
          (!p.hypoviscosityCutoffEnabled || k < p.hypoviscosityCutoff))
        hypo = p.hypoviscosity * std::pow(k2, p.hypoviscosityOrder);
      waveHyper += 2.0 * area * hyper * n;
      waveHypo += 2.0 * area * hypo * n;
      quadraticHyper += 2.0 * area * weight * hyper * n;
      quadraticHypo += 2.0 * area * weight * hypo * n;
      const std::size_t bin = static_cast<std::size_t>(k / binWidth);
      if (bin >= bins)
        throw std::logic_error("spectrum bin count is too small");
      waveSpectrum[bin] += area * n;
      quadraticSpectrum[bin] += area * weight * n;
      const double gamma = ginzburgLandauFactor(p, k2);
      const Complex conservativeNonlinear =
          Complex(1.0, gamma) * nonlinear[index];
      hamiltonianRate[index] =
          Complex(0.0, -weight) * w[index] + conservativeNonlinear;
      const Complex hamiltonianGradient =
          weight * w[index] + Complex(0.0, 1.0) * conservativeNonlinear;
      if (stochasticForcing) {
        const double amplitude = forcingAmplitude[index];
        expectedFullEnergyInjection += area * weight * amplitude * amplitude;
        if (!stochasticQuarticInjectionWeight.empty())
          expectedFullEnergyInjection +=
              area * stochasticQuarticInjectionWeight[index] * n;
      } else if (deterministicForcing) {
        expectedFullEnergyInjection +=
            2.0 * area *
            std::real(std::conj(hamiltonianGradient) * forcingAmplitude[index]);
      }
      totalEnergyHyper += 2.0 * area * hyper *
                          std::real(std::conj(hamiltonianGradient) * w[index]);
      totalEnergyHypo += 2.0 * area * hypo *
                         std::real(std::conj(hamiltonianGradient) * w[index]);
      const double transfer =
          -2.0 * area * std::real(std::conj(w[index]) * conservativeNonlinear);
      waveShell[bin] += transfer;
      quadraticShell[bin] += weight * transfer;
    }
  }
  SpectralField square, conservativeSquareRate;
  transform.projectedSquareSpectra(w, hamiltonianRate, square,
                                   conservativeSquareRate);
  // The two-pass cubic backend is the gradient of
  // (g/2) * ||P(psi^2)||^2, where P retains the base spectral band.
  // Using base-grid |psi|^4 here would describe a different Hamiltonian.
  double nonlinearEnergy = 0.0;
  for (const Complex value : square)
    nonlinearEnergy += std::norm(value);
  nonlinearEnergy *= 0.5 * area * p.nonlinearityCoefficient;
  const double totalEnergy = kineticEnergy + potentialEnergy + nonlinearEnergy;
  std::vector<double> fullEnergyShell = quadraticShell;
  for (std::size_t y = 0; y < p.ny; ++y)
    for (std::size_t x = 0; x < p.nx; ++x) {
      const std::size_t index = spectralIndex(x, y, p.nx);
      const std::size_t bin = static_cast<std::size_t>(
          std::sqrt(waveNumberSquared(p, x, y)) / binWidth);
      fullEnergyShell[bin] -=
          area * p.nonlinearityCoefficient *
          std::real(std::conj(square[index]) * conservativeSquareRate[index]);
    }
  std::vector<double> waveFlux(bins), fullEnergyFlux(bins);
  double cumulativeWave = 0.0, cumulativeEnergy = 0.0;
  for (std::size_t i = 0; i < bins; ++i) {
    cumulativeWave += waveShell[i];
    cumulativeEnergy += fullEnergyShell[i];
    waveFlux[i] = cumulativeWave;
    fullEnergyFlux[i] = cumulativeEnergy;
  }
  const auto finiteVector = [](const std::vector<double> &values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
  };
  if (!std::isfinite(totalEnergy) || !std::isfinite(waveAction) ||
      !std::isfinite(totalEnergyHypo) || !std::isfinite(totalEnergyHyper) ||
      !std::isfinite(expectedFullEnergyInjection) ||
      !finiteVector(waveSpectrum) || !finiteVector(quadraticSpectrum) ||
      !finiteVector(waveFlux) || !finiteVector(fullEnergyFlux))
    throw std::runtime_error(
        "diagnostics became non-finite; reduce the time step or coefficients");
  ++avg.count;
  for (std::size_t i = 0; i < bins; ++i) {
    avg.waveActionSpectrum[i] += waveSpectrum[i];
    avg.quadraticEnergySpectrum[i] += quadraticSpectrum[i];
    avg.waveActionFlux[i] += waveFlux[i];
    avg.fullEnergyFlux[i] += fullEnergyFlux[i];
  }
  auto diagnostics = numericOutput(p.outputDirectory / "diagnostics.csv",
                                   std::ios::out | std::ios::app);
  diagnostics << time << ',' << frame << ',' << totalEnergy << ','
              << kineticEnergy << ',' << potentialEnergy << ','
              << nonlinearEnergy << ',' << waveAction << ',' << waveHypo << ','
              << waveHyper << ',' << quadraticHypo << ',' << quadraticHyper
              << ',' << totalEnergyHypo << ',' << totalEnergyHyper << ','
              << expectedFullEnergyInjection << '\n';
  closeChecked(diagnostics, "failed while writing diagnostics.csv");
  auto spectra = numericOutput(p.outputDirectory / "spectra.csv",
                               std::ios::out | std::ios::app);
  auto fluxes = numericOutput(p.outputDirectory / "fluxes.csv",
                              std::ios::out | std::ios::app);
  for (std::size_t i = 0; i < bins; ++i) {
    const double k = static_cast<double>(i) * binWidth;
    spectra << time << ',' << frame << ',' << k << ',' << waveSpectrum[i] << ','
            << quadraticSpectrum[i] << ','
            << avg.waveActionSpectrum[i] / static_cast<double>(avg.count) << ','
            << avg.quadraticEnergySpectrum[i] / static_cast<double>(avg.count)
            << '\n';
    fluxes << time << ',' << frame << ',' << k << ',' << waveFlux[i] << ','
           << fullEnergyFlux[i] << ','
           << avg.waveActionFlux[i] / static_cast<double>(avg.count) << ','
           << avg.fullEnergyFlux[i] / static_cast<double>(avg.count) << '\n';
  }
  closeChecked(spectra, "failed while writing spectra.csv");
  closeChecked(fluxes, "failed while writing fluxes.csv");
  if (p.writeModeDiagnostics) {
    auto modes = numericOutput(p.outputDirectory / "modes.csv",
                               std::ios::out | std::ios::app);
    modes << time << ',' << frame;
    for (const auto [x, y] :
         {std::pair{1UL, 0UL}, std::pair{0UL, 1UL}, std::pair{1UL, 1UL},
          std::pair{2UL, 1UL}, std::pair{0UL, 3UL}}) {
      if (x >= p.nx || y >= p.ny)
        modes << ",,";
      else {
        const Complex value = w[spectralIndex(x, y, p.nx)];
        modes << ',' << value.real() << ',' << value.imag();
      }
    }
    modes << '\n';
    closeChecked(modes, "failed while writing modes.csv");
  }
  return totalEnergy;
}

void writeRestart(const Parameters &p, double time, std::uint64_t frame,
                  const SpectralField &wavefunction,
                  const std::string &randomEngineState,
                  const std::string &randomDistributionState) {
  if (!finiteField(wavefunction))
    throw std::runtime_error("refusing to checkpoint non-finite coefficients");
  const auto binaryPath = checkpointPath(p, frame);
  if (std::filesystem::exists(binaryPath) && !p.overwriteOutput)
    throw std::runtime_error("refusing to overwrite spectral checkpoint: " +
                             binaryPath.string());
  const auto temporary = std::filesystem::path(binaryPath.string() + ".tmp");
  std::ofstream checkpoint(temporary, std::ios::binary | std::ios::trunc);
  if (!checkpoint)
    throw std::runtime_error("cannot write spectral checkpoint: " +
                             binaryPath.string());
  const char magic[8] = {'G', 'P', '2', 'D', 'C', 'P', '1', '\0'};
  const std::uint64_t nx = p.nx, ny = p.ny, count = wavefunction.size();
  checkpoint.write(magic, sizeof(magic));
  checkpoint.write(reinterpret_cast<const char *>(&nx), sizeof(nx));
  checkpoint.write(reinterpret_cast<const char *>(&ny), sizeof(ny));
  checkpoint.write(reinterpret_cast<const char *>(&count), sizeof(count));
  std::vector<double> buffer(2 * checkpointChunkComplexValues);
  for (std::size_t offset = 0; offset < wavefunction.size();) {
    const std::size_t chunk =
        std::min(checkpointChunkComplexValues, wavefunction.size() - offset);
    for (std::size_t i = 0; i < chunk; ++i) {
      buffer[2 * i] = wavefunction[offset + i].real();
      buffer[2 * i + 1] = wavefunction[offset + i].imag();
    }
    checkpoint.write(reinterpret_cast<const char *>(buffer.data()),
                     static_cast<std::streamsize>(2 * chunk * sizeof(double)));
    offset += chunk;
  }
  closeChecked(checkpoint, "failed while writing spectral checkpoint");
  std::filesystem::rename(temporary, binaryPath);
  const auto metadata = p.dataDirectory / "restart_state.txt";
  const auto metadataTemporary = p.dataDirectory / "restart_state.tmp";
  std::ofstream out(metadataTemporary, std::ios::trunc);
  if (!out)
    throw std::runtime_error("cannot write restart metadata");
  out << std::setprecision(17) << "gp2d_restart_v1\ntime " << time << "\nframe "
      << frame << "\nnx " << p.nx << "\nny " << p.ny << "\naspectRatio "
      << p.aspectRatio << "\nrandomSeed " << p.randomSeed << "\nrandomEngine "
      << randomEngineState << "\nrandomDistribution " << randomDistributionState
      << '\n';
  closeChecked(out, "cannot write restart metadata");
  std::filesystem::rename(metadataTemporary, metadata);
}

void writeForcingFiles(const Parameters &p,
                       const std::vector<double> &amplitude,
                       std::size_t forcedModes, double waveCoefficient,
                       double quadraticEnergyCoefficient) {
  std::ofstream summary(p.outputDirectory / "forcing_summary.csv",
                        std::ios::trunc);
  summary << "enabled,profile,temporal_type,forced_modes,wave_action_injection_"
             "coefficient,quadratic_energy_injection_coefficient\n"
          << std::boolalpha << p.forcingEnabled << ','
          << forcingProfileName(p.forcingProfile) << ',';
  if (!p.forcingEnabled)
    summary << "disabled";
  else if (p.forcingProfile == ForcingProfile::singleMode)
    summary << "deterministic";
  else
    summary << "stochastic";
  summary << ',' << forcedModes << ',';
  if (p.forcingEnabled && p.forcingProfile != ForcingProfile::singleMode)
    summary << std::setprecision(17) << waveCoefficient << ','
            << quadraticEnergyCoefficient;
  else
    summary << ',';
  summary << '\n';
  closeChecked(summary, "failed while writing forcing_summary.csv");
  auto spectrum = numericOutput(p.outputDirectory / "forcing_spectrum.csv");
  spectrum << "kx,ky,amplitude\n";
  if (p.forcingEnabled)
    for (std::size_t y = 0; y < p.ny; ++y)
      for (std::size_t x = 0; x < p.nx; ++x) {
        const double kx =
            2.0 * gpPi * static_cast<double>(signedWave(x, p.nx)) / p.lx();
        const double ky =
            2.0 * gpPi * static_cast<double>(signedWave(y, p.ny)) / p.ly();
        spectrum << kx << ',' << ky << ','
                 << amplitude[spectralIndex(x, y, p.nx)] << '\n';
      }
  closeChecked(spectrum, "failed while writing forcing_spectrum.csv");
}
