#include "vortex_diagnostics.hpp"
#include "io_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr const char *header = "time,frame,positive_vortices,negative_vortices,total_vortices";
constexpr const char *positionHeader =
    "time,frame,index,x,y,circulation,cell_x,cell_y,core_residual";

double phaseDifference(Complex left, Complex right) {
    return std::remainder(std::arg(right) - std::arg(left), 2.0 * gpPi);
}

void prepareCsv(const std::filesystem::path &path, const char *expectedHeader, bool restarting,
                bool overwrite) {
    if (restarting && std::filesystem::exists(path)) {
        std::ifstream input(path);
        std::string existingHeader;
        std::getline(input, existingHeader);
        if (existingHeader != expectedHeader)
            throw std::runtime_error("CSV header does not match this solver version: " +
                                     path.string());
        return;
    }
    if (!restarting && std::filesystem::exists(path) && !overwrite)
        throw std::runtime_error("refusing to overwrite existing output: " + path.string());
    std::ofstream output(path, std::ios::trunc);
    output << expectedHeader << '\n';
    closeChecked(output, "cannot initialize vortex diagnostics: " + path.string());
}

std::pair<Complex, std::pair<double, double>> locateBilinearZero(Complex a, Complex b, Complex c,
                                                                 Complex d) {
    const Complex mixed = a - b - d + c;
    double u = 0.5, v = 0.5;
    for (int iteration = 0; iteration < 12; ++iteration) {
        const Complex value = a + u * (b - a) + v * (d - a) + u * v * mixed;
        const Complex derivativeU = b - a + v * mixed;
        const Complex derivativeV = d - a + u * mixed;
        const double determinant =
            derivativeU.real() * derivativeV.imag() - derivativeV.real() * derivativeU.imag();
        if (std::abs(determinant) < 1.e-15)
            break;
        const double deltaU =
            (-value.real() * derivativeV.imag() + derivativeV.real() * value.imag()) / determinant;
        const double deltaV =
            (-derivativeU.real() * value.imag() + value.real() * derivativeU.imag()) / determinant;
        u += deltaU;
        v += deltaV;
        if (!std::isfinite(u) || !std::isfinite(v) || u < -0.25 || u > 1.25 || v < -0.25 ||
            v > 1.25) {
            u = v = 0.5;
            break;
        }
        if (std::abs(deltaU) + std::abs(deltaV) < 1.e-12)
            break;
    }
    u = std::clamp(u, 0.0, 1.0);
    v = std::clamp(v, 0.0, 1.0);
    return {a + u * (b - a) + v * (d - a) + u * v * mixed, {u, v}};
}
} // namespace

std::vector<DetectedVortex> detectVortices(const Parameters &parameters,
                                           const std::vector<Complex> &physical) {
    if (physical.size() != parameters.nx * parameters.ny)
        throw std::runtime_error("invalid physical field for vortex detection");
    std::vector<DetectedVortex> vortices;
    for (std::size_t y = 0; y < parameters.ny; ++y) {
        const std::size_t nextY = (y + 1) % parameters.ny;
        for (std::size_t x = 0; x < parameters.nx; ++x) {
            const std::size_t nextX = (x + 1) % parameters.nx;
            const Complex a = physical[spectralIndex(x, y, parameters.nx)];
            const Complex b = physical[spectralIndex(nextX, y, parameters.nx)];
            const Complex c = physical[spectralIndex(nextX, nextY, parameters.nx)];
            const Complex d = physical[spectralIndex(x, nextY, parameters.nx)];
            const double winding = phaseDifference(a, b) + phaseDifference(b, c) +
                                   phaseDifference(c, d) + phaseDifference(d, a);
            const int charge = winding > gpPi ? 1 : winding < -gpPi ? -1 : 0;
            if (charge == 0)
                continue;
            const auto [coreValue, location] = locateBilinearZero(a, b, c, d);
            const auto [u, v] = location;
            double physicalX = parameters.lx() * (static_cast<double>(x) + u) / parameters.nx;
            double physicalY = parameters.ly() * (static_cast<double>(y) + v) / parameters.ny;
            physicalX = std::fmod(physicalX, parameters.lx());
            physicalY = std::fmod(physicalY, parameters.ly());
            vortices.push_back({charge, physicalX, physicalY, std::abs(coreValue), x, y});
        }
    }
    return vortices;
}

void prepareVortexDiagnostics(const Parameters &parameters, bool restarting) {
    if (!parameters.writeVortexDiagnostics)
        return;
    prepareCsv(parameters.outputDirectory / "vortices.csv", header, restarting,
               parameters.overwriteOutput);
    prepareCsv(parameters.outputDirectory / "vortex_positions.csv", positionHeader, restarting,
               parameters.overwriteOutput);
}

void writeVortexDiagnostics(const Parameters &parameters, BaseTransform &transform, double time,
                            std::uint64_t frame, const SpectralField &wavefunction) {
    if (!parameters.writeVortexDiagnostics)
        return;
    std::vector<Complex> physical;
    transform.inverse(wavefunction, physical);
    const std::vector<DetectedVortex> vortices = detectVortices(parameters, physical);
    std::uint64_t positive = 0, negative = 0;
    for (const DetectedVortex &vortex : vortices) {
        if (vortex.charge > 0)
            ++positive;
        else
            ++negative;
    }
    const auto path = parameters.outputDirectory / "vortices.csv";
    std::ofstream output(path, std::ios::app);
    output << std::scientific << std::setprecision(12) << time << ',' << frame << ',' << positive
           << ',' << negative << ',' << positive + negative << '\n';
    closeChecked(output, "failed while writing vortices.csv");
    auto positions =
        std::ofstream(parameters.outputDirectory / "vortex_positions.csv", std::ios::app);
    positions << std::scientific << std::setprecision(12);
    for (std::size_t index = 0; index < vortices.size(); ++index) {
        const DetectedVortex &vortex = vortices[index];
        positions << time << ',' << frame << ',' << index << ',' << vortex.x << ',' << vortex.y
                  << ',' << vortex.charge << ',' << vortex.cellX << ',' << vortex.cellY << ','
                  << vortex.coreResidual << '\n';
    }
    closeChecked(positions, "failed while writing vortex_positions.csv");
}
