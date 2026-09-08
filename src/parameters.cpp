#include "parameters.hpp"
#include "io_utils.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

bool parseBool(const std::string &text, const std::string &key) {
  if (text == "true" || text == "1")
    return true;
  if (text == "false" || text == "0")
    return false;
  throw std::runtime_error(key + " must be true or false, got: " + text);
}

template <class T>
T parseNumber(const std::string &text, const std::string &key) {
  if constexpr (std::is_unsigned_v<T>)
    if (!text.empty() && text.front() == '-')
      throw std::runtime_error(key + " cannot be negative: " + text);
  std::istringstream input(text);
  T value{};
  input >> value;
  if (!input || !(input >> std::ws).eof())
    throw std::runtime_error("invalid value for " + key + ": " + text);
  return value;
}

std::string trim(const std::string &value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};
  return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

template <class T, std::size_t N>
bool parseMember(
    std::string_view key, const std::string &value, Parameters &parameters,
    const std::pair<std::string_view, T Parameters::*> (&members)[N]) {
  const auto entry = std::find_if(
      std::begin(members), std::end(members),
      [key](const auto &candidate) { return candidate.first == key; });
  if (entry == std::end(members))
    return false;
  if constexpr (std::is_same_v<T, bool>)
    parameters.*entry->second = parseBool(value, std::string(key));
  else if constexpr (std::is_same_v<T, std::filesystem::path>)
    parameters.*entry->second = value;
  else
    parameters.*entry->second = parseNumber<T>(value, std::string(key));
  return true;
}

Integrator parseIntegrator(const std::string &text) {
  if (text == "etd2")
    return Integrator::etd2;
  if (text == "etd3")
    return Integrator::etd3;
  if (text == "etd4")
    return Integrator::etd4;
  if (text == "rk2")
    return Integrator::integratingFactorRk2;
  throw std::runtime_error("integrator must be etd2, etd3, etd4, or rk2");
}

ForcingProfile parseForcingProfile(const std::string &text) {
  if (text == "annulus")
    return ForcingProfile::annulus;
  if (text == "gaussian")
    return ForcingProfile::gaussian;
  if (text == "exponential")
    return ForcingProfile::exponential;
  if (text == "logNormal")
    return ForcingProfile::logNormal;
  if (text == "singleMode")
    return ForcingProfile::singleMode;
  throw std::runtime_error(
      "forcingProfile must be annulus, gaussian, exponential, logNormal, or "
      "singleMode");
}

bool isFinite(double value) { return std::isfinite(value); }
} // namespace

double Parameters::lx() const { return 2.0 * pi * aspectRatio; }
double Parameters::ly() const { return 2.0 * pi; }

std::size_t Parameters::spectrumBins() const {
  const double binWidth = std::min(2.0 * pi / lx(), 2.0 * pi / ly());
  const double maximumKx = pi * static_cast<double>(nx) / lx();
  const double maximumKy = pi * static_cast<double>(ny) / ly();
  const double count =
      std::floor(std::hypot(maximumKx, maximumKy) / binWidth) + 1.0;
  if (!isFinite(count) || count < 1.0 ||
      count >= static_cast<double>(std::numeric_limits<std::size_t>::max()))
    throw std::runtime_error("aspectRatio produces an invalid spectrum grid");
  return static_cast<std::size_t>(count);
}

const char *integratorName(Integrator integrator) {
  switch (integrator) {
  case Integrator::etd2:
    return "etd2";
  case Integrator::etd3:
    return "etd3";
  case Integrator::etd4:
    return "etd4";
  case Integrator::integratingFactorRk2:
    return "rk2";
  }
  throw std::logic_error("unknown integrator");
}

const char *forcingProfileName(ForcingProfile profile) {
  switch (profile) {
  case ForcingProfile::annulus:
    return "annulus";
  case ForcingProfile::gaussian:
    return "gaussian";
  case ForcingProfile::exponential:
    return "exponential";
  case ForcingProfile::logNormal:
    return "logNormal";
  case ForcingProfile::singleMode:
    return "singleMode";
  }
  throw std::logic_error("unknown forcing profile");
}

Parameters readParameters(const std::filesystem::path &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open parameter file: " + path.string());
  Parameters p;
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (const auto comment = line.find('#'); comment != std::string::npos)
      line.erase(comment);
    line = trim(line);
    if (line.empty())
      continue;
    std::replace(line.begin(), line.end(), '=', ' ');
    std::istringstream fields(line);
    std::string key, value, extra;
    fields >> key >> value;
    if (key.empty() || value.empty() || (fields >> extra))
      throw std::runtime_error("invalid parameter line " +
                               std::to_string(lineNumber));

    static constexpr std::pair<std::string_view, std::size_t Parameters::*>
        sizes[]{{"nx", &Parameters::nx}, {"ny", &Parameters::ny}};
    static constexpr std::pair<std::string_view, std::uint64_t Parameters::*>
        counts[]{{"numberOfSteps", &Parameters::numberOfSteps},
                 {"outputIntervalSteps", &Parameters::outputIntervalSteps},
                 {"randomSeed", &Parameters::randomSeed}};
    static constexpr std::pair<std::string_view, double Parameters::*> reals[]{
        {"aspectRatio", &Parameters::aspectRatio},
        {"timeStep", &Parameters::timeStep},
        {"dispersionCoefficient", &Parameters::dispersionCoefficient},
        {"nonlinearityCoefficient", &Parameters::nonlinearityCoefficient},
        {"chemicalPotential", &Parameters::chemicalPotential},
        {"hyperviscosity", &Parameters::hyperviscosity},
        {"hyperviscosityOrder", &Parameters::hyperviscosityOrder},
        {"hyperviscosityCutoff", &Parameters::hyperviscosityCutoff},
        {"hypoviscosity", &Parameters::hypoviscosity},
        {"hypoviscosityOrder", &Parameters::hypoviscosityOrder},
        {"hypoviscosityCutoff", &Parameters::hypoviscosityCutoff},
        {"ginzburgLandauDamping", &Parameters::ginzburgLandauDamping},
        {"ginzburgLandauCutoff", &Parameters::ginzburgLandauCutoff},
        {"forcingWavenumber", &Parameters::forcingWavenumber},
        {"forcingWidth", &Parameters::forcingWidth},
        {"forcingAmplitude", &Parameters::forcingAmplitude},
        {"forcingShapeOrder", &Parameters::forcingShapeOrder},
        {"forcingLogWidth", &Parameters::forcingLogWidth},
        {"targetWaveActionInjectionRate",
         &Parameters::targetWaveActionInjectionRate}};
    static constexpr std::pair<std::string_view, bool Parameters::*> booleans[]{
        {"hyperviscosityCutoffEnabled",
         &Parameters::hyperviscosityCutoffEnabled},
        {"hypoviscosityCutoffEnabled", &Parameters::hypoviscosityCutoffEnabled},
        {"forcingEnabled", &Parameters::forcingEnabled},
        {"writeModeDiagnostics", &Parameters::writeModeDiagnostics},
        {"overwriteOutput", &Parameters::overwriteOutput}};
    static constexpr std::pair<std::string_view,
                               std::filesystem::path Parameters::*>
        paths[]{{"initialConditionFile", &Parameters::initialConditionFile},
                {"dataDirectory", &Parameters::dataDirectory},
                {"outputDirectory", &Parameters::outputDirectory}};

    const bool recognized = parseMember(key, value, p, sizes) ||
                            parseMember(key, value, p, counts) ||
                            parseMember(key, value, p, reals) ||
                            parseMember(key, value, p, booleans) ||
                            parseMember(key, value, p, paths);
    if (key == "integrator")
      p.integrator = parseIntegrator(value);
    else if (key == "forcingProfile")
      p.forcingProfile = parseForcingProfile(value);
    else if (key == "threadCount")
      p.threadCount = parseNumber<int>(value, key);
    else if (!recognized)
      throw std::runtime_error("unknown parameter key on line " +
                               std::to_string(lineNumber) + ": " + key);
  }
  validateParameters(p);
  return p;
}

void validateParameters(const Parameters &p) {
  if (p.nx < 4 || p.ny < 4 || p.nx % 2 != 0 || p.ny % 2 != 0)
    throw std::runtime_error("nx and ny must be even and at least four");
  const auto fftLimit =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (p.nx > 2 * (fftLimit / 3) || p.ny > 2 * (fftLimit / 3))
    throw std::runtime_error(
        "3/2-rule grid dimensions exceed FFT library limits");
  const auto allocationLimit = std::min(
      static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()),
      std::numeric_limits<std::size_t>::max() / sizeof(std::complex<double>));
  const auto validateProduct = [&](std::size_t a, std::size_t b,
                                   const char *description) {
    if (a && b > allocationLimit / a)
      throw std::runtime_error(std::string(description) +
                               " exceeds addressable array limits");
  };
  validateProduct(p.nx, p.ny, "base grid");
  validateProduct(p.mx(), p.my(), "dealiased grid");
  if (!(p.aspectRatio > 0.0) || !isFinite(p.aspectRatio))
    throw std::runtime_error("aspectRatio must be finite and positive");
  (void)p.spectrumBins();
  if (!(p.timeStep > 0.0) || !isFinite(p.timeStep))
    throw std::runtime_error("timeStep must be finite and positive");
  if (p.numberOfSteps == 0 || p.outputIntervalSteps == 0)
    throw std::runtime_error(
        "numberOfSteps and outputIntervalSteps must be positive");
  if (p.threadCount < 0)
    throw std::runtime_error("threadCount cannot be negative");
  for (const auto [value, name] :
       {std::pair{p.dispersionCoefficient, "dispersionCoefficient"},
        std::pair{p.nonlinearityCoefficient, "nonlinearityCoefficient"},
        std::pair{p.chemicalPotential, "chemicalPotential"},
        std::pair{p.hyperviscosity, "hyperviscosity"},
        std::pair{p.hyperviscosityOrder, "hyperviscosityOrder"},
        std::pair{p.hypoviscosity, "hypoviscosity"},
        std::pair{p.hypoviscosityOrder, "hypoviscosityOrder"},
        std::pair{p.ginzburgLandauDamping, "ginzburgLandauDamping"}})
    if (!isFinite(value))
      throw std::runtime_error(std::string(name) + " must be finite");
  if (p.hyperviscosity < 0.0 || p.hypoviscosity < 0.0 ||
      p.ginzburgLandauDamping < 0.0)
    throw std::runtime_error("dissipation coefficients cannot be negative");
  if (p.hyperviscosityOrder < 0.0)
    throw std::runtime_error("hyperviscosityOrder cannot be negative");
  if (p.hyperviscosityCutoff < 0.0 || p.hypoviscosityCutoff < 0.0 ||
      p.ginzburgLandauCutoff < 0.0 || !isFinite(p.hyperviscosityCutoff) ||
      !isFinite(p.hypoviscosityCutoff) || !isFinite(p.ginzburgLandauCutoff))
    throw std::runtime_error(
        "dissipation cutoffs must be finite and nonnegative");
  if (p.forcingEnabled) {
    if (!(p.forcingWavenumber > 0.0) || p.forcingWidth < 0.0 ||
        p.forcingAmplitude < 0.0 || p.forcingShapeOrder <= 0.0 ||
        p.forcingLogWidth <= 0.0 || p.targetWaveActionInjectionRate < 0.0 ||
        !isFinite(p.forcingWavenumber) || !isFinite(p.forcingWidth) ||
        !isFinite(p.forcingAmplitude) || !isFinite(p.forcingShapeOrder) ||
        !isFinite(p.forcingLogWidth) ||
        !isFinite(p.targetWaveActionInjectionRate))
      throw std::runtime_error("forcing parameters are invalid or non-finite");
    if (p.forcingProfile == ForcingProfile::singleMode &&
        (p.forcingWavenumber != std::floor(p.forcingWavenumber) ||
         p.forcingWavenumber >=
             static_cast<double>(std::min(p.nx, p.ny)) / 2.0))
      throw std::runtime_error(
          "singleMode forcingWavenumber must be an integer "
          "below both Nyquist modes");
    if (p.forcingProfile == ForcingProfile::gaussian && !(p.forcingWidth > 0.0))
      throw std::runtime_error(
          "forcingWidth must be positive for gaussian forcing");
    if (p.forcingProfile == ForcingProfile::singleMode &&
        p.targetWaveActionInjectionRate > 0.0)
      throw std::runtime_error(
          "targetWaveActionInjectionRate applies only to stochastic forcing");
  }
  if (!p.initialConditionFile.empty() &&
      !std::filesystem::exists(p.initialConditionFile))
    throw std::runtime_error("initialConditionFile does not exist: " +
                             p.initialConditionFile.string());
  if (p.dataDirectory.empty() || p.outputDirectory.empty())
    throw std::runtime_error("output directories cannot be empty");
}

void writeParameterRecord(const Parameters &p, const std::string &backend,
                          const std::filesystem::path &recordDirectory) {
  const auto path =
      (recordDirectory.empty() ? p.outputDirectory : recordDirectory) /
      "resolved_parameters.txt";
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("cannot write parameter record: " + path.string());
  out << std::boolalpha << std::setprecision(17) << "backend " << backend
      << '\n'
      << "nx " << p.nx << '\n'
      << "ny " << p.ny << '\n'
      << "aspectRatio " << p.aspectRatio << '\n'
      << "domainLengthX " << p.lx() << '\n'
      << "domainLengthY " << p.ly() << '\n'
      << "timeStep " << p.timeStep << '\n'
      << "numberOfSteps " << p.numberOfSteps << '\n'
      << "outputIntervalSteps " << p.outputIntervalSteps << '\n'
      << "integrator " << integratorName(p.integrator) << '\n'
      << "dispersionCoefficient " << p.dispersionCoefficient << '\n'
      << "nonlinearityCoefficient " << p.nonlinearityCoefficient << '\n'
      << "chemicalPotential " << p.chemicalPotential << '\n'
      << "hyperviscosity " << p.hyperviscosity << '\n'
      << "hyperviscosityOrder " << p.hyperviscosityOrder << '\n'
      << "hyperviscosityCutoffEnabled " << p.hyperviscosityCutoffEnabled << '\n'
      << "hyperviscosityCutoff " << p.hyperviscosityCutoff << '\n'
      << "hypoviscosity " << p.hypoviscosity << '\n'
      << "hypoviscosityOrder " << p.hypoviscosityOrder << '\n'
      << "hypoviscosityCutoffEnabled " << p.hypoviscosityCutoffEnabled << '\n'
      << "hypoviscosityCutoff " << p.hypoviscosityCutoff << '\n'
      << "ginzburgLandauDamping " << p.ginzburgLandauDamping << '\n'
      << "ginzburgLandauCutoff " << p.ginzburgLandauCutoff << '\n'
      << "forcingEnabled " << p.forcingEnabled << '\n'
      << "forcingProfile " << forcingProfileName(p.forcingProfile) << '\n'
      << "forcingWavenumber " << p.forcingWavenumber << '\n'
      << "forcingWidth " << p.forcingWidth << '\n'
      << "forcingAmplitude " << p.forcingAmplitude << '\n'
      << "forcingShapeOrder " << p.forcingShapeOrder << '\n'
      << "forcingLogWidth " << p.forcingLogWidth << '\n'
      << "targetWaveActionInjectionRate " << p.targetWaveActionInjectionRate
      << '\n'
      << "randomSeed " << p.randomSeed << '\n'
      << "writeModeDiagnostics " << p.writeModeDiagnostics << '\n'
      << "threadCount " << p.threadCount << '\n'
      << "overwriteOutput " << p.overwriteOutput << '\n'
      << "initialConditionFile " << p.initialConditionFile.string() << '\n'
      << "dataDirectory " << p.dataDirectory.string() << '\n'
      << "outputDirectory " << p.outputDirectory.string() << '\n';
  closeChecked(out, "failed while writing parameter record: " + path.string());
}
