#include "vortex_field.hpp"

#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
using Complex = std::complex<double>;

double minimumImage(double difference, double length) { return std::remainder(difference, length); }

std::string uncomment(std::string line) {
    if (const auto comment = line.find('#'); comment != std::string::npos)
        line.erase(comment);
    return line;
}

std::vector<std::string> splitCsv(const std::string &line) {
    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ','))
        fields.push_back(field);
    return fields;
}

void readPointVortexRunRecord(const std::filesystem::path &path,
                              PointVortexMetadata &metadata) {
    std::ifstream input(path);
    if (!input)
        return;
    std::string key;
    while (input >> key) {
        if (key == "boundaryCondition") {
            std::string value;
            if (!(input >> value))
                throw std::runtime_error("invalid boundaryCondition in " + path.string());
            metadata.geometry = value;
        } else if (key == "boxLengthX" || key == "boxLengthY") {
            double value = 0.0;
            if (!(input >> value) || !std::isfinite(value) || value <= 0.0)
                throw std::runtime_error("invalid " + key + " in " + path.string());
            if (key == "boxLengthX")
                metadata.lengthX = value;
            else
                metadata.lengthY = value;
        } else {
            std::string remainder;
            std::getline(input, remainder);
        }
    }
    if (input.bad())
        throw std::runtime_error("failed while reading PointVortex run record: " +
                                 path.string());
}

template <class T> T parseField(const std::string &text, const std::string &description) {
    std::istringstream input(text);
    T value{};
    if (!(input >> value) || !(input >> std::ws).eof())
        throw std::runtime_error("invalid " + description + ": " + text);
    return value;
}

int windingFromCirculation(double circulation, std::size_t lineNumber) {
    if (!std::isfinite(circulation) || circulation == 0.0)
        throw std::runtime_error("invalid zero or non-finite circulation on line " +
                                 std::to_string(lineNumber));
    return std::signbit(circulation) ? -1 : 1;
}

Complex scaledThetaOne(Complex z, double aspect) {
    // theta_1(z,q) with its positive real factor 2*q^(1/4) removed.  Only its
    // phase is needed. q = exp(-pi*aspect), where aspect = Ly/Lx.
    const double logQ = -gpPi * aspect;
    Complex sum{};
    constexpr std::size_t maximumTerms = 100000;
    constexpr double tolerance = 2.0e-15;
    const double imaginaryMagnitude = std::abs(z.imag());
    for (std::size_t n = 0; n < maximumTerms; ++n) {
        const double logCoefficient = logQ * static_cast<double>(n) * static_cast<double>(n + 1);
        const double coefficient = std::exp(logCoefficient);
        const double harmonic = static_cast<double>(2 * n + 1);
        Complex term = coefficient * std::sin(harmonic * z);
        if (n % 2 != 0)
            term = -term;
        sum += term;
        // Do not use the current term alone as the stopping test: an individual
        // sine can vanish while later harmonics do not. This envelope bounds
        // |sin(harmonic*z)| and is already decreasing once the second condition
        // holds.
        const double logEnvelope = logCoefficient + harmonic * imaginaryMagnitude;
        const bool envelopeDecreasing =
            logQ * static_cast<double>(2 * n + 2) + 2.0 * imaginaryMagnitude < 0.0;
        const double logThreshold = std::log(tolerance * std::max(1.0, std::abs(sum)));
        if (n >= 2 && envelopeDecreasing && logEnvelope <= logThreshold)
            return sum;
    }
    throw std::runtime_error("periodic phase series did not converge; check the domain aspect "
                             "ratio");
}

void validateFiniteField(const std::vector<Complex> &field, const std::string &description) {
    if (!std::all_of(field.begin(), field.end(), [](Complex value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag());
        }))
        throw std::runtime_error(description + " contains a non-finite value");
}
} // namespace

PointVortexMetadata readPointVortexMetadata(const std::filesystem::path &path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open point-vortex file: " + path.string());

    PointVortexMetadata metadata;
    std::string line;
    while (std::getline(input, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            continue;
        if (line[first] != '#') {
            metadata.trajectory =
                line.find("time,frame,index,x,y,circulation", first) == first;
            break;
        }

        std::istringstream fields(line.substr(first + 1));
        std::string field;
        while (fields >> field) {
            const auto equals = field.find('=');
            if (equals == std::string::npos)
                continue;
            const std::string key = field.substr(0, equals);
            const std::string value = field.substr(equals + 1);
            if (key == "geometry")
                metadata.geometry = value;
            else if (key == "box_length") {
                const double length = parseField<double>(value, "PointVortex box_length");
                if (!std::isfinite(length) || length <= 0.0)
                    throw std::runtime_error("PointVortex box_length must be finite and positive");
                metadata.lengthX = length;
                metadata.lengthY = length;
            }
        }
    }
    if (input.bad())
        throw std::runtime_error("failed while reading point-vortex file: " + path.string());

    if (metadata.trajectory)
        readPointVortexRunRecord(path.parent_path() / "resolved_parameters.txt", metadata);
    return metadata;
}

void validatePointVortexDomain(const PointVortexMetadata &metadata,
                               const Parameters &parameters, double coordinateScale) {
    if (!(coordinateScale > 0.0) || !std::isfinite(coordinateScale))
        throw std::invalid_argument("coordinate scale must be finite and positive");
    if (metadata.geometry && *metadata.geometry != "periodic")
        throw std::invalid_argument("PointVortex input uses geometry=" + *metadata.geometry +
                                    "; the GP initial-condition tools require a doubly "
                                    "periodic PointVortex state");

    const auto significantlyDifferent = [](double left, double right) {
        return std::abs(left - right) >
               1.0e-10 * std::max({1.0, std::abs(left), std::abs(right)});
    };
    const auto checkLength = [&](const std::optional<double> &source, double target,
                                 const char *axis) {
        if (source && significantlyDifferent(coordinateScale * *source, target)) {
            std::ostringstream message;
            message << "scaled PointVortex " << axis << " length "
                    << coordinateScale * *source << " does not match GP " << axis << " length "
                    << target << "; use matching domains or set --coordinate-scale";
            throw std::invalid_argument(message.str());
        }
    };
    checkLength(metadata.lengthX, parameters.lx(), "x");
    checkLength(metadata.lengthY, parameters.ly(), "y");
}

std::vector<PointVortex> readPointVortices(const std::filesystem::path &path,
                                           std::optional<std::uint64_t> trajectoryFrame) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open point-vortex file: " + path.string());

    std::vector<std::pair<std::uint64_t, PointVortex>> trajectory;
    std::vector<PointVortex> plain;
    bool trajectoryFormat = false;
    bool formatKnown = false;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        line = uncomment(std::move(line));
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
            continue;
        if (!formatKnown) {
            trajectoryFormat = line.find("time,frame,index,x,y,circulation") != std::string::npos;
            formatKnown = true;
            if (trajectoryFormat)
                continue;
        }

        if (trajectoryFormat) {
            const auto fields = splitCsv(line);
            if (fields.size() < 6)
                throw std::runtime_error("invalid trajectory row on line " +
                                         std::to_string(lineNumber));
            const double time = parseField<double>(fields[0], "trajectory time");
            const auto frame = parseField<std::uint64_t>(fields[1], "trajectory frame");
            (void)parseField<std::uint64_t>(fields[2], "trajectory vortex index");
            const double x = parseField<double>(fields[3], "vortex x coordinate");
            const double y = parseField<double>(fields[4], "vortex y coordinate");
            const double circulation = parseField<double>(fields[5], "vortex circulation");
            if (!std::isfinite(time) || !std::isfinite(x) || !std::isfinite(y))
                throw std::runtime_error("non-finite vortex coordinate on line " +
                                         std::to_string(lineNumber));
            trajectory.push_back({frame, {x, y, windingFromCirculation(circulation, lineNumber)}});
        } else {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream fields(line);
            double x = 0.0, y = 0.0, circulation = 0.0;
            std::string trailing;
            if (!(fields >> x >> y >> circulation) || (fields >> trailing) || !std::isfinite(x) ||
                !std::isfinite(y))
                throw std::runtime_error("invalid point-vortex row on line " +
                                         std::to_string(lineNumber));
            plain.push_back({x, y, windingFromCirculation(circulation, lineNumber)});
        }
    }
    if (input.bad())
        throw std::runtime_error("failed while reading point-vortex file: " + path.string());

    if (!trajectoryFormat) {
        if (trajectoryFrame)
            throw std::runtime_error("--frame applies only to trajectory.csv input");
        if (plain.empty())
            throw std::runtime_error("point-vortex file is empty: " + path.string());
        return plain;
    }
    if (trajectory.empty())
        throw std::runtime_error("point-vortex trajectory is empty: " + path.string());
    const std::uint64_t selectedFrame =
        trajectoryFrame.value_or(std::max_element(trajectory.begin(), trajectory.end(),
                                                  [](const auto &left, const auto &right) {
                                                      return left.first < right.first;
                                                  })
                                     ->first);
    std::vector<PointVortex> selected;
    for (const auto &[frame, vortex] : trajectory)
        if (frame == selectedFrame)
            selected.push_back(vortex);
    if (selected.empty())
        throw std::runtime_error("trajectory contains no rows for frame " +
                                 std::to_string(selectedFrame));
    return selected;
}

double padeVortexDensity(double radius) {
    if (!std::isfinite(radius) || radius < 0.0)
        throw std::invalid_argument("vortex-profile radius must be finite and nonnegative");
    const double a1 = 0.34010790700196714760;
    const double b1 = (2304.0 * std::pow(a1, 3) + 656.0 * std::pow(a1, 2) - 421.0 * a1 - 28.0) /
                      (7680.0 * std::pow(a1, 2) - 1680.0 * a1 - 330.0);
    const double a2 = a1 * (b1 - 0.25);
    const double b2 =
        ((737280.0 * std::pow(a1, 3) + 209920.0 * std::pow(a1, 2) - 134720.0 * a1 - 8960.0) * b1 -
         364544.0 * std::pow(a1, 3) + 70144.0 * std::pow(a1, 2) + 18256.0 * a1 + 393.0) /
        (2457600.0 * std::pow(a1, 2) - 537600.0 * a1 - 105600.0);
    const double a3 = a1 * (192.0 * b2 - 48.0 * b1 + 16.0 * a1 + 5.0) / 192.0;
    const double b3 =
        ((61440.0 * a1 - 9600.0) * b2 + (-30720.0 * std::pow(a1, 2) + 640.0 * a1 + 560.0) * b1 +
         8448.0 * std::pow(a1, 2) - 1056.0 * a1 - 21.0) /
        (368640.0 * a1 - 92160.0);
    const double a4 =
        (4608.0 * a1 * b3 - 1152.0 * a1 * b2 + b1 * (384.0 * std::pow(a1, 2) + 120.0 * a1) -
         128.0 * std::pow(a1, 2) - 7.0 * a1) /
        4608.0;
    const double r2 = radius * radius;
    const double numerator = r2 * (a1 + r2 * (a2 + r2 * (a3 + r2 * a4)));
    const double denominator = 1.0 + r2 * (b1 + r2 * (b2 + r2 * (b3 + r2 * a4)));
    const double density = numerator / denominator;
    if (!std::isfinite(density) || density < -1.0e-14)
        throw std::runtime_error("Padé vortex profile produced an invalid density");
    return std::clamp(density, 0.0, 1.0);
}

Complex periodicVortexPhaseFactor(double x, double y, const Parameters &parameters,
                                  const std::vector<PointVortex> &vortices, int phaseWindingX,
                                  int phaseWindingY) {
    const double lx = parameters.lx(), ly = parameters.ly();
    const bool useXAsThetaPeriod = lx >= ly;
    Complex factor{1.0, 0.0};
    double windingMoment = 0.0;
    int totalWinding = 0;
    for (const PointVortex &vortex : vortices) {
        if (vortex.winding != -1 && vortex.winding != 1)
            throw std::invalid_argument("only singly quantized vortices are supported");
        // Keep a centered representative. Moving one member of a neutral set
        // by a whole period leaves its position unchanged but changes the
        // uniform phase-winding sector of the theta-function construction.
        // Current PointVortex output is centered, so remainder() preserves its
        // intended zero-background-flow representative.
        const double vx = std::remainder(vortex.x, lx);
        const double vy = std::remainder(vortex.y, ly);
        // Use the longer side as the theta function's real period. This keeps
        // the imaginary argument bounded and avoids overflow for elongated
        // domains. y-i*x is a rotation of x+i*y and preserves vortex sign.
        const Complex argument = useXAsThetaPeriod ? gpPi * Complex((x - vx) / lx, (y - vy) / lx)
                                                   : gpPi * Complex((y - vy) / ly, -(x - vx) / ly);
        const double thetaAspect = useXAsThetaPeriod ? ly / lx : lx / ly;
        const Complex theta = scaledThetaOne(argument, thetaAspect);
        if (std::abs(theta) > 32.0 * std::numeric_limits<double>::min()) {
            const Complex unit = theta / std::abs(theta);
            factor *= vortex.winding > 0 ? unit : std::conj(unit);
        }
        windingMoment += static_cast<double>(vortex.winding) * (useXAsThetaPeriod ? vx : vy);
        totalWinding += vortex.winding;
    }
    if (totalWinding != 0)
        throw std::invalid_argument("a doubly periodic GP field requires zero net vortex winding");

    const double periodicityCorrection = useXAsThetaPeriod
                                             ? -2.0 * gpPi * y * windingMoment / (lx * ly)
                                             : 2.0 * gpPi * x * windingMoment / (lx * ly);
    const double correction =
        periodicityCorrection + 2.0 * gpPi *
                                    (static_cast<double>(phaseWindingX) * x / lx +
                                     static_cast<double>(phaseWindingY) * y / ly);
    factor *= std::exp(Complex(0.0, correction));
    const double magnitude = std::abs(factor);
    return magnitude == 0.0 ? Complex{1.0, 0.0} : factor / magnitude;
}

std::vector<Complex> imprintPointVortices(const Parameters &parameters,
                                          std::vector<PointVortex> vortices,
                                          const VortexImprintOptions &options) {
    if (!(options.backgroundDensity > 0.0) || !std::isfinite(options.backgroundDensity))
        throw std::invalid_argument("background density must be finite and positive");
    if (!(options.healingLength > 0.0) || !std::isfinite(options.healingLength))
        throw std::invalid_argument("healing length must be finite and positive");
    if (!(options.coordinateScale > 0.0) || !std::isfinite(options.coordinateScale))
        throw std::invalid_argument("coordinate scale must be finite and positive");
    if (vortices.empty())
        throw std::invalid_argument("at least one point vortex is required");

    const double lx = parameters.lx(), ly = parameters.ly();
    for (PointVortex &vortex : vortices) {
        vortex.x = std::remainder(options.coordinateScale * vortex.x, lx);
        vortex.y = std::remainder(options.coordinateScale * vortex.y, ly);
    }

    std::vector<Complex> field(parameters.nx * parameters.ny);
    const double backgroundAmplitude = std::sqrt(options.backgroundDensity);
    for (std::size_t yIndex = 0; yIndex < parameters.ny; ++yIndex) {
        const double y = ly * static_cast<double>(yIndex) / static_cast<double>(parameters.ny);
        for (std::size_t xIndex = 0; xIndex < parameters.nx; ++xIndex) {
            const double x = lx * static_cast<double>(xIndex) / static_cast<double>(parameters.nx);
            double amplitude = backgroundAmplitude;
            for (const PointVortex &vortex : vortices) {
                const double dx = minimumImage(x - vortex.x, lx);
                const double dy = minimumImage(y - vortex.y, ly);
                const double scaledRadius = std::hypot(dx, dy) / options.healingLength;
                amplitude *= std::sqrt(padeVortexDensity(scaledRadius));
            }
            field[spectralIndex(xIndex, yIndex, parameters.nx)] =
                amplitude * periodicVortexPhaseFactor(x, y, parameters, vortices,
                                                      options.phaseWindingX, options.phaseWindingY);
        }
    }
    validateFiniteField(field, "imprinted wavefunction");
    return field;
}

std::vector<Complex> readWavefunction(const std::filesystem::path &path,
                                      const Parameters &parameters) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open wavefunction: " + path.string());
    std::vector<Complex> field;
    field.reserve(parameters.nx * parameters.ny);
    double real = 0.0, imaginary = 0.0;
    while (input >> real) {
        if (!(input >> imaginary) || !std::isfinite(real) || !std::isfinite(imaginary))
            throw std::runtime_error("invalid real/imaginary pair in wavefunction: " +
                                     path.string());
        field.emplace_back(real, imaginary);
    }
    if (!input.eof())
        throw std::runtime_error("invalid value in wavefunction: " + path.string());
    if (field.size() != parameters.nx * parameters.ny)
        throw std::runtime_error("wavefunction must contain exactly 2*nx*ny numbers: " +
                                 path.string());
    return field;
}

void writeWavefunctionFile(const std::filesystem::path &path,
                           const std::vector<Complex> &wavefunction, const Parameters &parameters,
                           bool overwrite) {
    if (wavefunction.size() != parameters.nx * parameters.ny)
        throw std::invalid_argument("wavefunction has the wrong grid size");
    validateFiniteField(wavefunction, "wavefunction");
    if (std::filesystem::exists(path) && !overwrite)
        throw std::runtime_error("refusing to overwrite wavefunction: " + path.string());
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    std::ofstream output(temporary, std::ios::trunc);
    if (!output)
        throw std::runtime_error("cannot write wavefunction: " + path.string());
    output << std::scientific << std::setprecision(12);
    for (std::size_t y = 0; y < parameters.ny; ++y) {
        for (std::size_t x = 0; x < parameters.nx; ++x) {
            const Complex value = wavefunction[spectralIndex(x, y, parameters.nx)];
            output << value.real() << ' ' << value.imag();
            output << (x + 1 == parameters.nx ? '\n' : ' ');
        }
    }
    output.close();
    if (!output)
        throw std::runtime_error("failed while writing wavefunction: " + path.string());
    if (overwrite && std::filesystem::exists(path))
        std::filesystem::remove(path);
    std::filesystem::rename(temporary, path);
}
