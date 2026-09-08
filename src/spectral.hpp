#pragma once

#include "backend.hpp"

#include <cstddef>

constexpr double gpPi = 3.141592653589793238462643383279502884;

inline std::size_t spectralIndex(std::size_t x, std::size_t y,
                                 std::size_t width) {
  return y * width + x;
}

inline long signedWave(std::size_t index, std::size_t count) {
  return index <= (count - 1) / 2
             ? static_cast<long>(index)
             : static_cast<long>(index) - static_cast<long>(count);
}

inline bool retainedPaddedWave(long wave, std::size_t baseCount) {
  return wave >= -static_cast<long>(baseCount / 2) &&
         wave < static_cast<long>(baseCount / 2);
}

inline std::size_t paddedIndexForBase(std::size_t index, std::size_t baseCount,
                                      std::size_t paddedCount) {
  const long wave = signedWave(index, baseCount);
  return wave >= 0
             ? static_cast<std::size_t>(wave)
             : static_cast<std::size_t>(static_cast<long>(paddedCount) + wave);
}

inline std::size_t paddedIndexForBaseMode(const Parameters &p,
                                          std::size_t index) {
  const std::size_t x = index % p.nx, y = index / p.nx;
  return spectralIndex(paddedIndexForBase(x, p.nx, p.mx()),
                       paddedIndexForBase(y, p.ny, p.my()), p.mx());
}

double waveNumberSquared(const Parameters &parameters, std::size_t x,
                         std::size_t y);
double ginzburgLandauFactor(const Parameters &parameters, double k2);
void enforceStateConstraints(SpectralField &field,
                             const Parameters &parameters);
