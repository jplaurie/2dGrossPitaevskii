#pragma once

#include "backend.hpp"

#include <cstdint>
#include <vector>

struct DetectedVortex {
    int charge = 0;
    double x = 0.0, y = 0.0;
    double coreResidual = 0.0;
    std::size_t cellX = 0, cellY = 0;
};

std::vector<DetectedVortex> detectVortices(const Parameters &parameters,
                                           const std::vector<Complex> &physical);

void prepareVortexDiagnostics(const Parameters &parameters, bool restarting);
void writeVortexDiagnostics(const Parameters &parameters, double time, std::uint64_t frame,
                            const std::vector<Complex> &physical);
