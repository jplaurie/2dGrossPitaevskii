#pragma once

#include <cmath>
#include <complex>
#include <cstdint>

inline std::uint64_t mixNoiseBits(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// A Box--Muller pair keyed by step and global Fourier-mode index. One shared
// mt19937_64 draw supplies each step key; ranks only evaluate their own modes.
// Midpoint uniforms exclude both endpoints, keeping log(u) finite.
inline std::complex<double> indexedComplexNormal(std::uint64_t stepKey, std::uint64_t mode) {
    const auto key = stepKey ^ mixNoiseBits(mode);
    const double u = std::ldexp(static_cast<double>(mixNoiseBits(key) >> 12) + 0.5, -52);
    const double v =
        std::ldexp(static_cast<double>(mixNoiseBits(key + 0x9e3779b97f4a7c15ULL) >> 12) + 0.5, -52);
    const double radius = std::sqrt(-std::log(u)); // Circular complex variance one.
    const double phase = 6.283185307179586476925286766559 * v;
    return radius * std::complex<double>(std::cos(phase), std::sin(phase));
}
