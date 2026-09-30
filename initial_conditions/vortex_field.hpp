#pragma once

#include "parameters.hpp"

#include <complex>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

struct PointVortex {
    double x = 0.0;
    double y = 0.0;
    int winding = 0;
};

struct VortexImprintOptions {
    double backgroundDensity = 1.0;
    double healingLength = 1.0;
    double coordinateScale = 1.0;
    int phaseWindingX = 0;
    int phaseWindingY = 0;
};

struct PointVortexMetadata {
    std::optional<std::string> geometry;
    std::optional<double> lengthX;
    std::optional<double> lengthY;
    bool trajectory = false;
};

std::vector<PointVortex>
readPointVortices(const std::filesystem::path &path,
                  std::optional<std::uint64_t> trajectoryFrame = std::nullopt);
PointVortexMetadata readPointVortexMetadata(const std::filesystem::path &path);
void validatePointVortexDomain(const PointVortexMetadata &metadata,
                               const Parameters &parameters, double coordinateScale);

double padeVortexDensity(double radiusInHealingLengths);

std::complex<double> periodicVortexPhaseFactor(double x, double y, const Parameters &parameters,
                                               const std::vector<PointVortex> &vortices,
                                               int phaseWindingX = 0, int phaseWindingY = 0);

std::vector<std::complex<double>> imprintPointVortices(const Parameters &parameters,
                                                       std::vector<PointVortex> vortices,
                                                       const VortexImprintOptions &options);

std::vector<std::complex<double>> readWavefunction(const std::filesystem::path &path,
                                                   const Parameters &parameters);
void writeWavefunctionFile(const std::filesystem::path &path,
                           const std::vector<std::complex<double>> &wavefunction,
                           const Parameters &parameters, bool overwrite);
