#include "spectral.hpp"

#include <cmath>
#include <stdexcept>

double waveNumberSquared(const Parameters &parameters, std::size_t x, std::size_t y) {
    const double kx =
        2.0 * gpPi * static_cast<double>(signedWave(x, parameters.nx)) / parameters.lx();
    const double ky =
        2.0 * gpPi * static_cast<double>(signedWave(y, parameters.ny)) / parameters.ly();
    return kx * kx + ky * ky;
}

double ginzburgLandauFactor(const Parameters &parameters, double k2) {
    return parameters.ginzburgLandauDamping > 0.0 && std::sqrt(k2) > parameters.ginzburgLandauCutoff
               ? parameters.ginzburgLandauDamping
               : 0.0;
}

void enforceStateConstraints(SpectralField &field, const Parameters &parameters) {
    if (field.size() != parameters.nx * parameters.ny)
        throw std::runtime_error("invalid spectral field size");
    if (parameters.hypoviscosity > 0.0 && parameters.hypoviscosityOrder < 0.0)
        field[0] = Complex{};
}
