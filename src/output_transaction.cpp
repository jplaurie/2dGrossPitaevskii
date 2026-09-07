#include "output.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {
constexpr std::array csvNames{"diagnostics.csv", "spectra.csv", "fluxes.csv",
                              "modes.csv"};

std::array<std::filesystem::path, 2> frameFiles(const Parameters &p,
                                                std::uint64_t frame) {
  std::ostringstream suffix;
  suffix << std::setw(8) << std::setfill('0') << frame;
  return {p.dataDirectory / ("wavefunction_" + suffix.str() + ".dat"),
          p.dataDirectory / ("checkpoint_" + suffix.str() + ".bin")};
}

std::filesystem::path appended(const std::filesystem::path &path,
                               const char *suffix) {
  return path.string() + suffix;
}

std::optional<std::uint64_t> committedFrame(const Parameters &p) {
  const auto path = p.dataDirectory / "restart_state.txt";
  if (!std::filesystem::exists(path))
    return std::nullopt;
  std::ifstream in(path);
  std::string format, key;
  double time = 0.0;
  std::uint64_t frame = 0;
  if (!(in >> format) || format != "gp2d_restart_v1" || !(in >> key >> time) ||
      key != "time" || !(in >> key >> frame) || key != "frame")
    throw std::runtime_error(
        "cannot recover output with malformed restart metadata");
  return frame;
}

struct Journal {
  std::uint64_t frame{};
  bool previousMetadata{};
  std::uint64_t previousFrame{};
  std::array<bool, 4> csvExisted{};
  std::array<std::uintmax_t, 4> csvSizes{};
  std::array<bool, 2> frameExisted{};
};

Journal readJournal(const Parameters &p) {
  std::ifstream in(p.dataDirectory / "output_transaction.txt");
  std::string format, directory;
  Journal journal;
  if (!(in >> format >> std::quoted(directory) >> journal.frame >>
        journal.previousMetadata >> journal.previousFrame) ||
      format != "gp2d_output_transaction_v1")
    throw std::runtime_error("malformed output transaction journal");
  if (std::filesystem::canonical(p.outputDirectory) !=
      std::filesystem::path(directory))
    throw std::runtime_error(
        "recover the interrupted run using its original outputDirectory: " +
        directory);
  for (std::size_t i = 0; i < csvNames.size(); ++i)
    if (!(in >> journal.csvExisted[i] >> journal.csvSizes[i]))
      throw std::runtime_error("malformed CSV offsets in output journal");
  for (auto &existed : journal.frameExisted)
    if (!(in >> existed))
      throw std::runtime_error("malformed frame records in output journal");
  if (!(in >> std::ws).eof())
    throw std::runtime_error("unexpected data in output transaction journal");
  return journal;
}

void removeJournal(const Parameters &p) {
  std::filesystem::remove(p.dataDirectory / "output_transaction.txt");
  std::filesystem::remove(p.dataDirectory / "output_transaction.tmp");
}
} // namespace

bool recoverOutputTransaction(const Parameters &p) {
  if (!std::filesystem::exists(p.dataDirectory / "output_transaction.txt"))
    return false;
  const Journal journal = readJournal(p);
  const auto committed = committedFrame(p);
  if (committed && *committed == journal.frame) {
    finishOutputTransaction(p);
    return false;
  }
  if (committed.has_value() != journal.previousMetadata ||
      (committed && *committed != journal.previousFrame))
    throw std::runtime_error(
        "restart metadata does not match the interrupted output transaction");
  for (std::size_t i = 0; i < csvNames.size(); ++i) {
    const auto path = p.outputDirectory / csvNames[i];
    if (journal.csvExisted[i] &&
        (!std::filesystem::exists(path) ||
         std::filesystem::file_size(path) < journal.csvSizes[i]))
      throw std::runtime_error("committed CSV data is missing: " +
                               path.string());
  }
  for (std::size_t i = 0; i < csvNames.size(); ++i) {
    const auto path = p.outputDirectory / csvNames[i];
    if (journal.csvExisted[i])
      std::filesystem::resize_file(path, journal.csvSizes[i]);
    else
      std::filesystem::remove(path);
  }
  const auto files = frameFiles(p, journal.frame);
  for (std::size_t i = 0; i < files.size(); ++i) {
    const auto backup = appended(files[i], ".previous");
    if (std::filesystem::exists(backup))
      std::filesystem::rename(backup, files[i]);
    else if (!journal.frameExisted[i])
      std::filesystem::remove(files[i]);
    std::filesystem::remove(appended(files[i], ".tmp"));
  }
  std::filesystem::remove(p.dataDirectory / "restart_state.tmp");
  removeJournal(p);
  std::cout << "recovered interrupted output frame " << journal.frame << '\n';
  return journal.frame == 0 && !committed;
}

void beginOutputTransaction(const Parameters &p, std::uint64_t frame) {
  if (std::filesystem::exists(p.dataDirectory / "output_transaction.txt"))
    throw std::runtime_error("an output transaction is already active");
  Journal journal;
  journal.frame = frame;
  const auto previous = committedFrame(p);
  journal.previousMetadata = previous.has_value();
  journal.previousFrame = previous.value_or(0);
  const auto files = frameFiles(p, frame);
  for (std::size_t i = 0; i < files.size(); ++i) {
    journal.frameExisted[i] = std::filesystem::exists(files[i]);
    if (journal.frameExisted[i] &&
        (!p.overwriteOutput || !std::filesystem::is_regular_file(files[i])))
      throw std::runtime_error("refusing to overwrite output frame file: " +
                               files[i].string());
    for (const char *suffix : {".tmp", ".previous"})
      if (std::filesystem::exists(appended(files[i], suffix)))
        throw std::runtime_error("untracked temporary output file: " +
                                 appended(files[i], suffix).string());
  }
  for (std::size_t i = 0; i < csvNames.size(); ++i) {
    const auto path = p.outputDirectory / csvNames[i];
    journal.csvExisted[i] = std::filesystem::exists(path);
    if (journal.csvExisted[i])
      journal.csvSizes[i] = std::filesystem::file_size(path);
  }
  const auto temporary = p.dataDirectory / "output_transaction.tmp";
  std::ofstream out(temporary);
  out << "gp2d_output_transaction_v1\n"
      << std::quoted(std::filesystem::canonical(p.outputDirectory).string())
      << '\n'
      << journal.frame << ' ' << journal.previousMetadata << ' '
      << journal.previousFrame << '\n';
  for (std::size_t i = 0; i < csvNames.size(); ++i)
    out << journal.csvExisted[i] << ' ' << journal.csvSizes[i] << '\n';
  for (const bool existed : journal.frameExisted)
    out << existed << '\n';
  out.close();
  if (!out)
    throw std::runtime_error("cannot write output transaction journal");
  std::filesystem::rename(temporary,
                          p.dataDirectory / "output_transaction.txt");
  for (std::size_t i = 0; i < files.size(); ++i)
    if (journal.frameExisted[i])
      std::filesystem::rename(files[i], appended(files[i], ".previous"));
}

void finishOutputTransaction(const Parameters &p) {
  const Journal journal = readJournal(p);
  if (committedFrame(p) != std::optional{journal.frame})
    throw std::runtime_error("cannot finish an uncommitted output transaction");
  for (const auto &path : frameFiles(p, journal.frame)) {
    std::filesystem::remove(appended(path, ".previous"));
    std::filesystem::remove(appended(path, ".tmp"));
  }
  removeJournal(p);
}

void writeRunRecords(const Parameters &p, const std::string &backend,
                     double time, std::uint64_t frame,
                     const std::vector<double> &amplitude,
                     std::size_t forcedModes, double waveCoefficient,
                     double quadraticEnergyCoefficient) {
  const auto history = p.outputDirectory / "segments";
  std::filesystem::create_directories(history);
  std::filesystem::path segment;
  for (std::uint64_t index = 1;; ++index) {
    std::ostringstream name;
    name << "segment_" << std::setw(8) << std::setfill('0') << index;
    segment = history / name.str();
    if (std::filesystem::create_directory(segment))
      break;
  }
  writeParameterRecord(p, backend, segment);
  Parameters segmentParameters = p;
  segmentParameters.outputDirectory = segment;
  writeForcingFiles(segmentParameters, amplitude, forcedModes, waveCoefficient,
                    quadraticEnergyCoefficient);
  std::ofstream manifest(segment / "segment.txt");
  manifest << std::setprecision(17) << "startTime " << time << "\nstartFrame "
           << frame << "\nstochasticUpdate exact_linear_covariance_v1\n";
  manifest.close();
  if (!manifest)
    throw std::runtime_error("cannot write run segment record");
  for (const char *name : {"resolved_parameters.txt", "forcing_summary.csv",
                           "forcing_spectrum.csv"}) {
    const auto temporary = p.outputDirectory / (std::string(name) + ".tmp");
    std::filesystem::copy_file(
        segment / name, temporary,
        std::filesystem::copy_options::overwrite_existing);
    std::filesystem::rename(temporary, p.outputDirectory / name);
  }
}
