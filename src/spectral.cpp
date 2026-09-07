#include "spectral.hpp"

#include <cmath>
#include <stdexcept>

double waveNumberSquared(const Parameters &p, std::size_t x, std::size_t y) {
  const double kx =
      2.0 * gpPi * static_cast<double>(signedWave(x, p.nx)) / p.lx();
  const double ky =
      2.0 * gpPi * static_cast<double>(signedWave(y, p.ny)) / p.ly();
  return kx * kx + ky * ky;
}

double ginzburgLandauFactor(const Parameters &p, double k2) {
  return p.ginzburgLandauDamping > 0.0 && std::sqrt(k2) > p.ginzburgLandauCutoff
             ? p.ginzburgLandauDamping
             : 0.0;
}

void enforceStateConstraints(SpectralField &field, const Parameters &p) {
  if (field.size() != p.nx * p.ny)
    throw std::runtime_error("invalid spectral field size");
  if (p.hypoviscosity > 0.0 && p.hypoviscosityOrder < 0.0)
    field[0] = Complex{};
}
