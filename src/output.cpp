#include "output.hpp"
#include "hdf5_io.hpp"
#include "io_utils.hpp"
#include "spectral.hpp"
#include "vortex_diagnostics.hpp"

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

std::filesystem::path wavefunctionPath(const Parameters &parameters, std::uint64_t frame) {
    std::ostringstream name;
    name << "wavefunction_" << std::setw(8) << std::setfill('0') << frame << ".dat";
    return parameters.dataDirectory / name.str();
}

std::filesystem::path hdf5WavefunctionPath(const Parameters &parameters, std::uint64_t frame) {
    std::ostringstream name;
    name << "wavefunction_" << std::setw(8) << std::setfill('0') << frame << ".h5";
    return parameters.dataDirectory / name.str();
}

std::filesystem::path checkpointPath(const Parameters &parameters, std::uint64_t frame) {
    std::ostringstream name;
    name << "checkpoint_" << std::setw(8) << std::setfill('0') << frame << ".bin";
    return parameters.dataDirectory / name.str();
}

std::ofstream numericOutput(const std::filesystem::path &path,
                            std::ios::openmode mode = std::ios::out) {
    std::ofstream out(path, mode);
    if (!out)
        throw std::runtime_error("cannot write output file: " + path.string());
    out << std::scientific << std::setprecision(12);
    return out;
}

SpectralField readCheckpoint(const Parameters &parameters, std::uint64_t frame) {
    const auto path = checkpointPath(parameters, frame);
    std::ifstream input(path, std::ios::binary);
    if (!input)
        throw std::runtime_error("cannot open spectral checkpoint: " + path.string());
    char magic[8]{};
    std::uint64_t nx = 0, ny = 0, count = 0;
    input.read(magic, sizeof(magic));
    input.read(reinterpret_cast<char *>(&nx), sizeof(nx));
    input.read(reinterpret_cast<char *>(&ny), sizeof(ny));
    input.read(reinterpret_cast<char *>(&count), sizeof(count));
    if (!input || std::string(magic, 7) != "GP2DCP1" || nx != parameters.nx ||
        ny != parameters.ny || count != parameters.nx * parameters.ny)
        throw std::runtime_error("invalid spectral checkpoint header: " + path.string());
    SpectralField field(static_cast<std::size_t>(count));
    std::vector<double> buffer(2 * checkpointChunkComplexValues);
    for (std::size_t offset = 0; offset < field.size();) {
        const std::size_t chunk = std::min(checkpointChunkComplexValues, field.size() - offset);
        input.read(reinterpret_cast<char *>(buffer.data()),
                   static_cast<std::streamsize>(2 * chunk * sizeof(double)));
        if (!input)
            break;
        for (std::size_t i = 0; i < chunk; ++i)
            field[offset + i] = Complex(buffer[2 * i], buffer[2 * i + 1]);
        offset += chunk;
    }
    if (!input || input.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("invalid spectral checkpoint payload: " + path.string());
    if (!finiteField(field))
        throw std::runtime_error("spectral checkpoint contains non-finite values");
    return field;
}

std::vector<Complex> readPhysicalField(const std::filesystem::path &path,
                                       const Parameters &parameters) {
    const std::size_t count = parameters.nx * parameters.ny;
    if (path.extension() == ".h5" || path.extension() == ".hdf5") {
        Hdf5Field field = readHdf5Field(path);
        if (field.nx != parameters.nx || field.ny != parameters.ny)
            throw std::runtime_error("HDF5 initial-condition dimensions do not match nx and ny");
        const auto differs = [](double left, double right) {
            return std::abs(left - right) >
                   1.e-12 * std::max({1.0, std::abs(left), std::abs(right)});
        };
        if (differs(field.lengthX, parameters.lx()) ||
            differs(field.lengthY, parameters.ly()))
            throw std::runtime_error(
                "HDF5 initial-condition domain lengths do not match the configured domain");
        return std::move(field.wavefunction);
    }
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open wavefunction field: " + path.string());
    std::vector<double> values;
    double value = 0.0;
    while (input >> value) {
        if (!std::isfinite(value))
            throw std::runtime_error("wavefunction contains a non-finite value: " + path.string());
        values.push_back(value);
    }
    if (!input.eof())
        throw std::runtime_error("invalid value in wavefunction field: " + path.string());
    if (values.size() != 2 * count)
        throw std::runtime_error("wavefunction field must contain two values per grid point: " +
                                 path.string());
    std::vector<Complex> field(count);
    for (std::size_t i = 0; i < count; ++i)
        field[i] = Complex(values[2 * i], values[2 * i + 1]);
    return field;
}

void initializeCsv(const std::filesystem::path &path, const std::string &header, bool append,
                   bool overwrite) {
    if (append && std::filesystem::exists(path)) {
        std::ifstream input(path);
        std::string existingHeader;
        std::getline(input, existingHeader);
        if (existingHeader != header)
            throw std::runtime_error("CSV header does not match this solver version: " +
                                     path.string());
        return;
    }
    if (!append && std::filesystem::exists(path) && !overwrite)
        throw std::runtime_error("refusing to overwrite existing output: " + path.string());
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
    constexpr std::array names{
        "diagnostics.csv",     "spectra.csv",          "fluxes.csv",
        "modes.csv",           "vortices.csv",         "vortex_positions.csv",
        "forcing_summary.csv", "forcing_spectrum.csv", "resolved_parameters.txt"};
    return std::any_of(names.begin(), names.end(),
                       [&](const char *name) { return std::filesystem::exists(directory / name); });
}

} // namespace

RestartState readRestart(const Parameters &parameters, BaseTransform &transform, bool isRoot) {
    std::filesystem::create_directories(parameters.dataDirectory);
    std::filesystem::create_directories(parameters.outputDirectory);
    RestartState state;
    std::vector<Complex> physical;
    const auto restartPath = parameters.dataDirectory / "restart_state.txt";
    if (std::filesystem::exists(restartPath)) {
        std::ifstream input(restartPath);
        std::string format, key;
        std::getline(input, format);
        if (format != "gp2d_restart_v1")
            throw std::runtime_error("unsupported restart format: " + restartPath.string());
        std::size_t savedNx = 0, savedNy = 0;
        double savedAspect = 0.0;
        if (!(input >> key >> state.time) || key != "time" || !(input >> key >> state.frame) ||
            key != "frame" || !(input >> key >> savedNx) || key != "nx" ||
            !(input >> key >> savedNy) || key != "ny" || !(input >> key >> savedAspect) ||
            key != "aspectRatio" || !(input >> key >> state.randomSeed) || key != "randomSeed")
            throw std::runtime_error("malformed restart metadata: " + restartPath.string());
        input >> std::ws;
        std::getline(input, key, ' ');
        if (key != "randomEngine")
            throw std::runtime_error("restart is missing randomEngine state");
        std::getline(input, state.randomEngineState);
        std::getline(input, key, ' ');
        if (key != "randomDistribution")
            throw std::runtime_error("restart is missing randomDistribution state");
        std::getline(input, state.randomDistributionState);
        if (!std::isfinite(state.time) || state.time < 0.0 || savedNx != parameters.nx ||
            savedNy != parameters.ny || savedAspect != parameters.aspectRatio)
            throw std::runtime_error("restart grid, aspect ratio, or time is invalid");
        state.wavefunction = readCheckpoint(parameters, state.frame);
        state.restarting = true;
        if (isRoot)
            std::cout << "restarting at time " << state.time << " from frame " << state.frame
                      << '\n';
    } else {
        if (parameters.initialConditionFile.empty()) {
            physical.assign(parameters.nx * parameters.ny, Complex{});
            if (isRoot)
                std::cout << "starting from zero wavefunction\n";
        } else {
            physical = readPhysicalField(parameters.initialConditionFile, parameters);
            if (isRoot)
                std::cout << "starting from " << parameters.initialConditionFile << '\n';
        }
    }
    if (state.wavefunction.empty()) {
        transform.forward(physical, state.wavefunction);
        enforceStateConstraints(state.wavefunction, parameters);
    }
    if (!finiteField(state.wavefunction))
        throw std::runtime_error("initial wavefunction contains non-finite coefficients");
    return state;
}

void prepareOutputFiles(const Parameters &parameters, bool restarting, std::uint64_t restartFrame) {
    constexpr std::array csvNames{"diagnostics.csv", "spectra.csv",  "fluxes.csv",
                                  "modes.csv",       "vortices.csv", "vortex_positions.csv"};
    if (!restarting && !parameters.overwriteOutput &&
        (containsRunData(parameters.dataDirectory) ||
         containsSolverOutput(parameters.outputDirectory)))
        throw std::runtime_error("solver output already exists; use new "
                                 "directories or set overwriteOutput true");
    if (restarting) {
        for (const char *name : csvNames) {
            const auto path = parameters.outputDirectory / name;
            if (std::filesystem::exists(path) && lastCsvFrame(path) > restartFrame)
                throw std::runtime_error(path.string() + " contains frames newer than the restart");
        }
    }
    initializeCsv(parameters.outputDirectory / "diagnostics.csv",
                  "time,frame,total_energy,kinetic_energy,potential_energy,nonlinear_"
                  "energy,"
                  "wave_action,wave_action_dissipation_hypo,wave_action_dissipation_hyper,"
                  "quadratic_energy_dissipation_hypo,quadratic_energy_dissipation_hyper,"
                  "total_energy_dissipation_hypo,total_energy_dissipation_hyper,"
                  "expected_full_energy_injection",
                  restarting, parameters.overwriteOutput);
    initializeCsv(parameters.outputDirectory / "spectra.csv",
                  "time,frame,wavenumber,wave_action_spectrum,quadratic_energy_spectrum,"
                  "segment_mean_wave_action_spectrum,segment_mean_quadratic_energy_"
                  "spectrum",
                  restarting, parameters.overwriteOutput);
    initializeCsv(parameters.outputDirectory / "fluxes.csv",
                  "time,frame,wavenumber,wave_action_flux,full_energy_flux,"
                  "segment_mean_wave_action_flux,segment_mean_full_energy_flux",
                  restarting, parameters.overwriteOutput);
    if (parameters.writeModeDiagnostics)
        initializeCsv(parameters.outputDirectory / "modes.csv",
                      "time,frame,psi_1_0_real,psi_1_0_imag,psi_0_1_real,"
                      "psi_0_1_imag,psi_1_1_real,psi_1_1_imag,psi_2_1_real,"
                      "psi_2_1_imag,psi_0_3_real,psi_0_3_imag",
                      restarting, parameters.overwriteOutput);
    prepareVortexDiagnostics(parameters, restarting);
}

void writeWavefunction(const Parameters &parameters, BaseTransform &transform,
                       const SpectralField &wavefunction, double time, std::uint64_t frame) {
    std::vector<Complex> physical;
    transform.inverse(wavefunction, physical);
    if (parameters.fieldOutputFormat != FieldOutputFormat::hdf5) {
        const auto path = wavefunctionPath(parameters, frame);
        if (std::filesystem::exists(path) && !parameters.overwriteOutput)
            throw std::runtime_error("refusing to overwrite wavefunction snapshot: " +
                                     path.string());
        const auto temporary = std::filesystem::path(path.string() + ".tmp");
        auto out = numericOutput(temporary);
        for (std::size_t y = 0; y < parameters.ny; ++y) {
            for (std::size_t x = 0; x < parameters.nx; ++x) {
                const Complex value = physical[spectralIndex(x, y, parameters.nx)];
                out << value.real() << ' ' << value.imag();
                out << (x + 1 == parameters.nx ? '\n' : ' ');
            }
        }
        closeChecked(out, "failed while writing wavefunction snapshot");
        std::filesystem::rename(temporary, path);
    }
    if (parameters.fieldOutputFormat != FieldOutputFormat::text) {
        const auto path = hdf5WavefunctionPath(parameters, frame);
        if (std::filesystem::exists(path) && !parameters.overwriteOutput)
            throw std::runtime_error("refusing to overwrite HDF5 wavefunction snapshot: " +
                                     path.string());
        const auto temporary = std::filesystem::path(path.string() + ".tmp");
        writeHdf5Field(temporary, parameters, time, frame, physical);
        std::filesystem::rename(temporary, path);
    }
}

void writeRestart(const Parameters &parameters, double time, std::uint64_t frame,
                  const SpectralField &wavefunction, const std::string &randomEngineState,
                  const std::string &randomDistributionState) {
    if (!finiteField(wavefunction))
        throw std::runtime_error("refusing to checkpoint non-finite coefficients");
    const auto binaryPath = checkpointPath(parameters, frame);
    if (std::filesystem::exists(binaryPath) && !parameters.overwriteOutput)
        throw std::runtime_error("refusing to overwrite spectral checkpoint: " +
                                 binaryPath.string());
    const auto temporary = std::filesystem::path(binaryPath.string() + ".tmp");
    std::ofstream checkpoint(temporary, std::ios::binary | std::ios::trunc);
    if (!checkpoint)
        throw std::runtime_error("cannot write spectral checkpoint: " + binaryPath.string());
    const char magic[8] = {'G', 'P', '2', 'D', 'C', 'P', '1', '\0'};
    const std::uint64_t nx = parameters.nx, ny = parameters.ny, count = wavefunction.size();
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
    const auto metadata = parameters.dataDirectory / "restart_state.txt";
    const auto metadataTemporary = parameters.dataDirectory / "restart_state.tmp";
    std::ofstream out(metadataTemporary, std::ios::trunc);
    if (!out)
        throw std::runtime_error("cannot write restart metadata");
    out << std::setprecision(17) << "gp2d_restart_v1\ntime " << time << "\nframe " << frame
        << "\nnx " << parameters.nx << "\nny " << parameters.ny << "\naspectRatio "
        << parameters.aspectRatio << "\nrandomSeed " << parameters.randomSeed << "\nrandomEngine "
        << randomEngineState << "\nrandomDistribution " << randomDistributionState << '\n';
    closeChecked(out, "cannot write restart metadata");
    std::filesystem::rename(metadataTemporary, metadata);
}
