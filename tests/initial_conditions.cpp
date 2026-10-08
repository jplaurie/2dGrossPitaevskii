#include "vortex_diagnostics.hpp"
#include "vortex_field.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        require(padeVortexDensity(0.0) == 0.0, "Padé density is nonzero at the vortex core");
        require(padeVortexDensity(0.1) > 0.0 && padeVortexDensity(0.1) < 1.0,
                "Padé density has an invalid core value");
        require(std::abs(padeVortexDensity(100.0) - 1.0) < 1.0e-4,
                "Padé density does not approach the background");

        Parameters parameters;
        parameters.nx = 16;
        parameters.ny = 12;
        parameters.aspectRatio = 1.4;
        const std::vector<PointVortex> vortices{{1.1, 2.2, 1}, {4.4, 3.3, -1}};
        const auto value = periodicVortexPhaseFactor(0.73, 1.19, parameters, vortices);
        const auto periodicX =
            periodicVortexPhaseFactor(0.73 + parameters.lx(), 1.19, parameters, vortices);
        const auto periodicY =
            periodicVortexPhaseFactor(0.73, 1.19 + parameters.ly(), parameters, vortices);
        require(std::abs(value - periodicX) < 2.0e-12, "vortex phase is not periodic in x");
        require(std::abs(value - periodicY) < 2.0e-12, "vortex phase is not periodic in y");
        require(std::abs(std::abs(value) - 1.0) < 1.0e-14,
                "vortex phase factor is not unit magnitude");

        // A centered horizontal (+,-) dipole moves toward +y in the standard
        // c=-1/2 convention. Remove the singular + vortex contribution from a
        // numerical phase gradient and check the sign of the regular flow.
        Parameters dipoleParameters = parameters;
        dipoleParameters.aspectRatio = 1.0;
        const std::vector<PointVortex> horizontalDipole{{-0.7, 0.0, 1}, {0.7, 0.0, -1}};
        const double offset = 1.0e-2, difference = 1.0e-6;
        const auto phaseBelow = periodicVortexPhaseFactor(
            horizontalDipole[0].x + offset, -difference, dipoleParameters, horizontalDipole);
        const auto phaseAbove = periodicVortexPhaseFactor(
            horizontalDipole[0].x + offset, difference, dipoleParameters, horizontalDipole);
        const double totalGradient =
            std::arg(std::conj(phaseBelow) * phaseAbove) / (2.0 * difference);
        const double singularGradient = std::atan(difference / offset) / difference;
        require(totalGradient - singularGradient > 0.0,
                "centered dipole phase has the wrong drift direction");

        Parameters tallParameters = parameters;
        tallParameters.aspectRatio = 0.6;
        const auto tallValue = periodicVortexPhaseFactor(0.31, 1.27, tallParameters, vortices);
        const auto tallPeriodicX =
            periodicVortexPhaseFactor(0.31 + tallParameters.lx(), 1.27, tallParameters, vortices);
        const auto tallPeriodicY =
            periodicVortexPhaseFactor(0.31, 1.27 + tallParameters.ly(), tallParameters, vortices);
        require(std::abs(tallValue - tallPeriodicX) < 2.0e-12,
                "tall-domain vortex phase is not periodic in x");
        require(std::abs(tallValue - tallPeriodicY) < 2.0e-12,
                "tall-domain vortex phase is not periodic in y");

        VortexImprintOptions options;
        options.backgroundDensity = 1.0;
        options.healingLength = 0.2;
        const auto field = imprintPointVortices(parameters, vortices, options);
        require(field.size() == parameters.nx * parameters.ny,
                "imprinted field has the wrong size");
        require(std::all_of(field.begin(), field.end(),
                            [](std::complex<double> z) {
                                return std::isfinite(z.real()) && std::isfinite(z.imag());
                            }),
                "imprinted field contains a non-finite value");
        const auto detected = detectVortices(parameters, field);
        require(detected.size() == vortices.size(),
                "phase-winding detector found the wrong number of imprinted vortices");
        for (const PointVortex &expected : vortices) {
            const auto match =
                std::find_if(detected.begin(), detected.end(), [&](const auto &found) {
                    const double dx = std::remainder(found.x - expected.x, parameters.lx());
                    const double dy = std::remainder(found.y - expected.y, parameters.ly());
                    return found.charge == expected.winding &&
                           std::hypot(dx, dy) < 1.5 * std::max(parameters.lx() / parameters.nx,
                                                               parameters.ly() / parameters.ny);
                });
            require(match != detected.end(), "detected vortex position or charge is incorrect");
        }

        bool rejectedNonNeutral = false;
        try {
            (void)periodicVortexPhaseFactor(0.2, 0.3, parameters, {{1.0, 1.0, 1}});
        } catch (const std::invalid_argument &) {
            rejectedNonNeutral = true;
        }
        require(rejectedNonNeutral, "non-neutral periodic vortex state was accepted");

        PointVortexMetadata metadata;
        metadata.geometry = "periodic";
        metadata.lengthX = parameters.lx();
        metadata.lengthY = parameters.ly();
        validatePointVortexDomain(metadata, parameters, 1.0);
        bool rejectedWrongGeometry = false;
        try {
            metadata.geometry = "disk";
            validatePointVortexDomain(metadata, parameters, 1.0);
        } catch (const std::invalid_argument &) {
            rejectedWrongGeometry = true;
        }
        require(rejectedWrongGeometry, "non-periodic PointVortex metadata was accepted");
        metadata.geometry = "periodic";
        metadata.lengthX = 2.0;
        bool rejectedWrongLength = false;
        try {
            validatePointVortexDomain(metadata, parameters, 1.0);
        } catch (const std::invalid_argument &) {
            rejectedWrongLength = true;
        }
        require(rejectedWrongLength, "mismatched PointVortex domain length was accepted");

        const std::filesystem::path trajectoryDirectory = "gp2d_test_point_vortex_run";
        std::filesystem::remove_all(trajectoryDirectory);
        std::filesystem::create_directory(trajectoryDirectory);
        const std::filesystem::path trajectoryPath = trajectoryDirectory / "trajectory.csv";
        {
            std::ofstream trajectory(trajectoryPath);
            trajectory << "time,frame,index,x,y,circulation,u,v\n"
                       << "0,0,0,1,2,1,0,0\n"
                       << "0,0,1,3,4,-1,0,0\n"
                       << "1,2,0,5,6,1,0,0\n"
                       << "1,2,1,7,8,-1,0,0\n";
        }
        {
            std::ofstream record(trajectoryDirectory / "resolved_parameters.txt");
            record << std::setprecision(17) << "POINT_VORTEX_RUN_RECORD 1\n"
                   << "boundaryCondition periodic\n"
                   << "domainLengthX " << parameters.lx() << '\n'
                   << "domainLengthY " << parameters.ly() << '\n';
        }
        const auto latest = readPointVortices(trajectoryPath);
        const auto first = readPointVortices(trajectoryPath, 0);
        const auto trajectoryMetadata = readPointVortexMetadata(trajectoryPath);
        std::filesystem::remove_all(trajectoryDirectory);
        require(latest.size() == 2 && latest[0].x == 5.0 && latest[1].winding == -1,
                "latest trajectory frame was not selected correctly");
        require(first.size() == 2 && first[0].x == 1.0 && first[1].y == 4.0,
                "explicit trajectory frame was not selected correctly");
        require(trajectoryMetadata.trajectory &&
                    trajectoryMetadata.geometry == std::optional<std::string>("periodic") &&
                    trajectoryMetadata.lengthX && trajectoryMetadata.lengthY,
                "trajectory run metadata was not identified");
        validatePointVortexDomain(trajectoryMetadata, parameters, 1.0);

        std::cout << "periodic vortex phase and Padé profile agree with their invariants\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
