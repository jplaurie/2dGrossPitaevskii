#include "parameters.hpp"
#include "io_utils.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
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

template <class T> T parseNumber(const std::string &text, const std::string &key) {
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
    throw std::runtime_error("forcingProfile must be annulus, gaussian, exponential, logNormal, or "
                             "singleMode");
}

struct ParameterSetting {
    std::string key;
    std::string value;
    std::size_t lineNumber;
};

std::optional<ParameterSetting> parseParameterLine(std::string line, std::size_t lineNumber) {
    if (const auto comment = line.find('#'); comment != std::string::npos)
        line.erase(comment);
    line = trim(line);
    if (line.empty())
        return std::nullopt;

    std::replace(line.begin(), line.end(), '=', ' ');
    std::istringstream fields(line);
    ParameterSetting setting{{}, {}, lineNumber};
    std::string extra;
    fields >> setting.key >> setting.value;
    if (setting.key.empty() || setting.value.empty() || (fields >> extra))
        throw std::runtime_error("invalid parameter line " + std::to_string(lineNumber));
    return setting;
}

void applyParameter(const ParameterSetting &setting, Parameters &parameters) {
    const std::string &key = setting.key;
    const std::string &value = setting.value;
    const auto readNumber = [&]<class T>(T &destination) {
        destination = parseNumber<T>(value, key);
    };

    if (key == "nx")
        readNumber(parameters.nx);
    else if (key == "ny")
        readNumber(parameters.ny);
    else if (key == "aspectRatio")
        readNumber(parameters.aspectRatio);
    else if (key == "timeStep")
        readNumber(parameters.timeStep);
    else if (key == "numberOfSteps")
        readNumber(parameters.numberOfSteps);
    else if (key == "outputIntervalSteps")
        readNumber(parameters.outputIntervalSteps);
    else if (key == "integrator")
        parameters.integrator = parseIntegrator(value);
    else if (key == "dispersionCoefficient")
        readNumber(parameters.dispersionCoefficient);
    else if (key == "nonlinearityCoefficient")
        readNumber(parameters.nonlinearityCoefficient);
    else if (key == "chemicalPotential")
        readNumber(parameters.chemicalPotential);
    else if (key == "hyperviscosity")
        readNumber(parameters.hyperviscosity);
    else if (key == "hyperviscosityOrder")
        readNumber(parameters.hyperviscosityOrder);
    else if (key == "hyperviscosityCutoffEnabled")
        parameters.hyperviscosityCutoffEnabled = parseBool(value, key);
    else if (key == "hyperviscosityCutoff")
        readNumber(parameters.hyperviscosityCutoff);
    else if (key == "hypoviscosity")
        readNumber(parameters.hypoviscosity);
    else if (key == "hypoviscosityOrder")
        readNumber(parameters.hypoviscosityOrder);
    else if (key == "hypoviscosityCutoffEnabled")
        parameters.hypoviscosityCutoffEnabled = parseBool(value, key);
    else if (key == "hypoviscosityCutoff")
        readNumber(parameters.hypoviscosityCutoff);
    else if (key == "ginzburgLandauDamping")
        readNumber(parameters.ginzburgLandauDamping);
    else if (key == "ginzburgLandauCutoff")
        readNumber(parameters.ginzburgLandauCutoff);
    else if (key == "forcingEnabled")
        parameters.forcingEnabled = parseBool(value, key);
    else if (key == "forcingProfile")
        parameters.forcingProfile = parseForcingProfile(value);
    else if (key == "forcingWavenumber")
        readNumber(parameters.forcingWavenumber);
    else if (key == "forcingWidth")
        readNumber(parameters.forcingWidth);
    else if (key == "forcingAmplitude")
        readNumber(parameters.forcingAmplitude);
    else if (key == "forcingShapeOrder")
        readNumber(parameters.forcingShapeOrder);
    else if (key == "forcingLogWidth")
        readNumber(parameters.forcingLogWidth);
    else if (key == "targetWaveActionInjectionRate")
        readNumber(parameters.targetWaveActionInjectionRate);
    else if (key == "randomSeed")
        readNumber(parameters.randomSeed);
    else if (key == "writeModeDiagnostics")
        parameters.writeModeDiagnostics = parseBool(value, key);
    else if (key == "threadCount")
        readNumber(parameters.threadCount);
    else if (key == "overwriteOutput")
        parameters.overwriteOutput = parseBool(value, key);
    else if (key == "initialConditionFile")
        parameters.initialConditionFile = value;
    else if (key == "dataDirectory")
        parameters.dataDirectory = value;
    else if (key == "outputDirectory")
        parameters.outputDirectory = value;
    else
        throw std::runtime_error("unknown parameter key on line " +
                                 std::to_string(setting.lineNumber) + ": " + key);
}

bool isFinite(double value) { return std::isfinite(value); }
} // namespace

double Parameters::lx() const { return 2.0 * pi * aspectRatio; }
double Parameters::ly() const { return 2.0 * pi; }

std::size_t Parameters::spectrumBins() const {
    const double binWidth = std::min(2.0 * pi / lx(), 2.0 * pi / ly());
    const double maximumKx = pi * static_cast<double>(nx) / lx();
    const double maximumKy = pi * static_cast<double>(ny) / ly();
    const double count = std::floor(std::hypot(maximumKx, maximumKy) / binWidth) + 1.0;
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

Parameters readParameters(const std::filesystem::path &path, bool requireExistingInitialCondition) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open parameter file: " + path.string());
    Parameters parameters;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (const auto setting = parseParameterLine(line, lineNumber))
            applyParameter(*setting, parameters);
    }
    validateParameters(parameters, requireExistingInitialCondition);
    return parameters;
}

void validateParameters(const Parameters &parameters, bool requireExistingInitialCondition) {
    if (parameters.nx < 4 || parameters.ny < 4 || parameters.nx % 2 != 0 || parameters.ny % 2 != 0)
        throw std::runtime_error("nx and ny must be even and at least four");
    const auto fftLimit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (parameters.nx > 2 * (fftLimit / 3) || parameters.ny > 2 * (fftLimit / 3))
        throw std::runtime_error("3/2-rule grid dimensions exceed FFT library limits");
    const auto allocationLimit =
        std::min(static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()),
                 std::numeric_limits<std::size_t>::max() / sizeof(std::complex<double>));
    const auto validateProduct = [&](std::size_t a, std::size_t b, const char *description) {
        if (a && b > allocationLimit / a)
            throw std::runtime_error(std::string(description) +
                                     " exceeds addressable array limits");
    };
    validateProduct(parameters.nx, parameters.ny, "base grid");
    validateProduct(parameters.mx(), parameters.my(), "dealiased grid");
    if (!(parameters.aspectRatio > 0.0) || !isFinite(parameters.aspectRatio))
        throw std::runtime_error("aspectRatio must be finite and positive");
    (void)parameters.spectrumBins();
    if (!(parameters.timeStep > 0.0) || !isFinite(parameters.timeStep))
        throw std::runtime_error("timeStep must be finite and positive");
    if (parameters.numberOfSteps == 0 || parameters.outputIntervalSteps == 0)
        throw std::runtime_error("numberOfSteps and outputIntervalSteps must be positive");
    if (parameters.threadCount < 0)
        throw std::runtime_error("threadCount cannot be negative");
    for (const auto [value, name] :
         {std::pair{parameters.dispersionCoefficient, "dispersionCoefficient"},
          std::pair{parameters.nonlinearityCoefficient, "nonlinearityCoefficient"},
          std::pair{parameters.chemicalPotential, "chemicalPotential"},
          std::pair{parameters.hyperviscosity, "hyperviscosity"},
          std::pair{parameters.hyperviscosityOrder, "hyperviscosityOrder"},
          std::pair{parameters.hypoviscosity, "hypoviscosity"},
          std::pair{parameters.hypoviscosityOrder, "hypoviscosityOrder"},
          std::pair{parameters.ginzburgLandauDamping, "ginzburgLandauDamping"}})
        if (!isFinite(value))
            throw std::runtime_error(std::string(name) + " must be finite");
    if (parameters.hyperviscosity < 0.0 || parameters.hypoviscosity < 0.0 ||
        parameters.ginzburgLandauDamping < 0.0)
        throw std::runtime_error("dissipation coefficients cannot be negative");
    if (parameters.hyperviscosityOrder < 0.0)
        throw std::runtime_error("hyperviscosityOrder cannot be negative");
    if (parameters.hyperviscosityCutoff < 0.0 || parameters.hypoviscosityCutoff < 0.0 ||
        parameters.ginzburgLandauCutoff < 0.0 || !isFinite(parameters.hyperviscosityCutoff) ||
        !isFinite(parameters.hypoviscosityCutoff) || !isFinite(parameters.ginzburgLandauCutoff))
        throw std::runtime_error("dissipation cutoffs must be finite and nonnegative");
    if (parameters.forcingEnabled) {
        if (!(parameters.forcingWavenumber > 0.0) || parameters.forcingWidth < 0.0 ||
            parameters.forcingAmplitude < 0.0 || parameters.forcingShapeOrder <= 0.0 ||
            parameters.forcingLogWidth <= 0.0 || parameters.targetWaveActionInjectionRate < 0.0 ||
            !isFinite(parameters.forcingWavenumber) || !isFinite(parameters.forcingWidth) ||
            !isFinite(parameters.forcingAmplitude) || !isFinite(parameters.forcingShapeOrder) ||
            !isFinite(parameters.forcingLogWidth) ||
            !isFinite(parameters.targetWaveActionInjectionRate))
            throw std::runtime_error("forcing parameters are invalid or non-finite");
        if (parameters.forcingProfile == ForcingProfile::singleMode &&
            (parameters.forcingWavenumber != std::floor(parameters.forcingWavenumber) ||
             parameters.forcingWavenumber >=
                 static_cast<double>(std::min(parameters.nx, parameters.ny)) / 2.0))
            throw std::runtime_error("singleMode forcingWavenumber must be an integer "
                                     "below both Nyquist modes");
        if (parameters.forcingProfile == ForcingProfile::gaussian &&
            !(parameters.forcingWidth > 0.0))
            throw std::runtime_error("forcingWidth must be positive for gaussian forcing");
        if (parameters.forcingProfile == ForcingProfile::singleMode &&
            parameters.targetWaveActionInjectionRate > 0.0)
            throw std::runtime_error(
                "targetWaveActionInjectionRate applies only to stochastic forcing");
    }
    if (requireExistingInitialCondition && !parameters.initialConditionFile.empty() &&
        !std::filesystem::exists(parameters.initialConditionFile))
        throw std::runtime_error("initialConditionFile does not exist: " +
                                 parameters.initialConditionFile.string());
    if (parameters.dataDirectory.empty() || parameters.outputDirectory.empty())
        throw std::runtime_error("output directories cannot be empty");
}

void writeParameterRecord(const Parameters &parameters, const std::string &backend,
                          const std::filesystem::path &recordDirectory) {
    const auto path = (recordDirectory.empty() ? parameters.outputDirectory : recordDirectory) /
                      "resolved_parameters.txt";
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("cannot write parameter record: " + path.string());
    output << std::boolalpha << std::setprecision(17);
    const auto write = [&](const char *name, const auto &value) {
        output << name << ' ' << value << '\n';
    };
    write("backend", backend);
    write("nx", parameters.nx);
    write("ny", parameters.ny);
    write("aspectRatio", parameters.aspectRatio);
    write("domainLengthX", parameters.lx());
    write("domainLengthY", parameters.ly());
    write("timeStep", parameters.timeStep);
    write("numberOfSteps", parameters.numberOfSteps);
    write("outputIntervalSteps", parameters.outputIntervalSteps);
    write("integrator", integratorName(parameters.integrator));
    write("dispersionCoefficient", parameters.dispersionCoefficient);
    write("nonlinearityCoefficient", parameters.nonlinearityCoefficient);
    write("chemicalPotential", parameters.chemicalPotential);
    write("hyperviscosity", parameters.hyperviscosity);
    write("hyperviscosityOrder", parameters.hyperviscosityOrder);
    write("hyperviscosityCutoffEnabled", parameters.hyperviscosityCutoffEnabled);
    write("hyperviscosityCutoff", parameters.hyperviscosityCutoff);
    write("hypoviscosity", parameters.hypoviscosity);
    write("hypoviscosityOrder", parameters.hypoviscosityOrder);
    write("hypoviscosityCutoffEnabled", parameters.hypoviscosityCutoffEnabled);
    write("hypoviscosityCutoff", parameters.hypoviscosityCutoff);
    write("ginzburgLandauDamping", parameters.ginzburgLandauDamping);
    write("ginzburgLandauCutoff", parameters.ginzburgLandauCutoff);
    write("forcingEnabled", parameters.forcingEnabled);
    write("forcingProfile", forcingProfileName(parameters.forcingProfile));
    write("forcingWavenumber", parameters.forcingWavenumber);
    write("forcingWidth", parameters.forcingWidth);
    write("forcingAmplitude", parameters.forcingAmplitude);
    write("forcingShapeOrder", parameters.forcingShapeOrder);
    write("forcingLogWidth", parameters.forcingLogWidth);
    write("targetWaveActionInjectionRate", parameters.targetWaveActionInjectionRate);
    write("randomSeed", parameters.randomSeed);
    write("writeModeDiagnostics", parameters.writeModeDiagnostics);
    write("threadCount", parameters.threadCount);
    write("overwriteOutput", parameters.overwriteOutput);
    write("initialConditionFile", parameters.initialConditionFile.string());
    write("dataDirectory", parameters.dataDirectory.string());
    write("outputDirectory", parameters.outputDirectory.string());
    closeChecked(output, "failed while writing parameter record: " + path.string());
}
