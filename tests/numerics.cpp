#include "backend.hpp"
#include "fftw_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
Complex phase(double angle) { return std::exp(Complex(0.0, angle)); }

SpectralField directNonlinear(const Parameters &p, const SpectralField &w) {
  const std::size_t paddedCount = p.mx() * p.my();
  std::vector<Complex> psi(paddedCount), squareSpectrum(paddedCount),
      cubic(paddedCount);
  for (std::size_t py = 0; py < p.my(); ++py)
    for (std::size_t px = 0; px < p.mx(); ++px) {
      Complex sum{};
      for (std::size_t y = 0; y < p.ny; ++y)
        for (std::size_t x = 0; x < p.nx; ++x) {
          const double angle =
              2.0 * gpPi *
              (static_cast<double>(signedWave(x, p.nx)) * px / p.mx() +
               static_cast<double>(signedWave(y, p.ny)) * py / p.my());
          sum += w[spectralIndex(x, y, p.nx)] * phase(angle);
        }
      psi[spectralIndex(px, py, p.mx())] = sum;
    }
  for (std::size_t ky = 0; ky < p.my(); ++ky)
    for (std::size_t kx = 0; kx < p.mx(); ++kx) {
      Complex sum{};
      for (std::size_t py = 0; py < p.my(); ++py)
        for (std::size_t px = 0; px < p.mx(); ++px) {
          const double angle =
              -2.0 * gpPi *
              (static_cast<double>(signedWave(kx, p.mx())) * px / p.mx() +
               static_cast<double>(signedWave(ky, p.my())) * py / p.my());
          const Complex value = psi[spectralIndex(px, py, p.mx())];
          sum += value * value * phase(angle);
        }
      if (retainedPaddedWave(signedWave(kx, p.mx()), p.nx) &&
          retainedPaddedWave(signedWave(ky, p.my()), p.ny))
        squareSpectrum[spectralIndex(kx, ky, p.mx())] =
            sum / static_cast<double>(paddedCount);
    }
  for (std::size_t py = 0; py < p.my(); ++py)
    for (std::size_t px = 0; px < p.mx(); ++px) {
      Complex sum{};
      for (std::size_t ky = 0; ky < p.my(); ++ky)
        for (std::size_t kx = 0; kx < p.mx(); ++kx) {
          const double angle =
              2.0 * gpPi *
              (static_cast<double>(signedWave(kx, p.mx())) * px / p.mx() +
               static_cast<double>(signedWave(ky, p.my())) * py / p.my());
          sum += squareSpectrum[spectralIndex(kx, ky, p.mx())] * phase(angle);
        }
      const std::size_t i = spectralIndex(px, py, p.mx());
      cubic[i] = sum * std::conj(psi[i]);
    }
  SpectralField result(p.nx * p.ny);
  for (std::size_t y = 0; y < p.ny; ++y)
    for (std::size_t x = 0; x < p.nx; ++x) {
      Complex sum{};
      for (std::size_t py = 0; py < p.my(); ++py)
        for (std::size_t px = 0; px < p.mx(); ++px) {
          const double angle =
              -2.0 * gpPi *
              (static_cast<double>(signedWave(x, p.nx)) * px / p.mx() +
               static_cast<double>(signedWave(y, p.ny)) * py / p.my());
          sum += cubic[spectralIndex(px, py, p.mx())] * phase(angle);
        }
      result[spectralIndex(x, y, p.nx)] = p.nonlinearityCoefficient * sum /
                                          static_cast<double>(paddedCount) /
                                          Complex(0.0, 1.0);
    }
  return result;
}
} // namespace

int main(int argc, char **argv) {
  try {
    backendInitialize(argc, argv);
    Parameters p;
    p.nx = 4;
    p.ny = 6; // Its 3/2-padded dimension is odd, exercising FFT indexing.
    p.nonlinearityCoefficient = 1.7;
    p.hyperviscosity = 0.0;
    p.hypoviscosity = 0.0;
    p.ginzburgLandauDamping = 0.0;
    p.threadCount = 1;
    if (signedWave(4, 9) != 4 || signedWave(5, 9) != -4)
      throw std::runtime_error("odd-length Fourier indexing is incorrect");
    SpectralField w(p.nx * p.ny);
    w[spectralIndex(0, 0, p.nx)] = {0.7, -0.2};
    w[spectralIndex(1, 0, p.nx)] = {0.12, 0.05};
    w[spectralIndex(0, 1, p.nx)] = {-0.08, 0.04};
    w[spectralIndex(1, 1, p.nx)] = {0.03, -0.07};
    w[spectralIndex(3, 2, p.nx)] = {0.02, 0.01};
    auto backend = makeBackend(p);
    SpectralField actual;
    backend->evaluate(w, actual);
    const SpectralField expected = directNonlinear(p, w);
    double error = 0.0, magnitude = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
      error = std::max(error, std::abs(actual[i] - expected[i]));
      magnitude = std::max(magnitude, std::abs(expected[i]));
    }
    if (error > 2.e-11 * std::max(1.0, magnitude))
      throw std::runtime_error(
          "dealiased cubic term disagrees with direct DFT: " +
          std::to_string(error));

    BaseTransform transform(p);
    SpectralField conservativeRate(w.size()), dampingRate(w.size()),
        zeroRate(w.size());
    for (std::size_t y = 0; y < p.ny; ++y)
      for (std::size_t x = 0; x < p.nx; ++x) {
        const std::size_t i = spectralIndex(x, y, p.nx);
        const double k2 = waveNumberSquared(p, x, y);
        const double weight =
            -p.dispersionCoefficient * k2 + p.chemicalPotential;
        conservativeRate[i] = Complex(0.0, -weight) * w[i] + actual[i];
        dampingRate[i] = -(0.2 + 0.01 * k2) * w[i];
      }
    SpectralField square, conservativeSquareRate;
    transform.projectedSquareSpectra(w, conservativeRate, square,
                                     conservativeSquareRate);
    SpectralField dampingSquareRate;
    transform.projectedSquareSpectra(w, dampingRate, square, dampingSquareRate);
    const double area = p.lx() * p.ly();
    const auto quadraticRate = [&](const SpectralField &rate) {
      double value = 0.0;
      for (std::size_t y = 0; y < p.ny; ++y)
        for (std::size_t x = 0; x < p.nx; ++x) {
          const std::size_t i = spectralIndex(x, y, p.nx);
          const double weight =
              -p.dispersionCoefficient * waveNumberSquared(p, x, y) +
              p.chemicalPotential;
          value += 2.0 * area * weight * std::real(std::conj(w[i]) * rate[i]);
        }
      return value;
    };
    const auto quarticRate = [&](const SpectralField &rate) {
      double value = 0.0;
      for (std::size_t i = 0; i < square.size(); ++i)
        value += area * p.nonlinearityCoefficient *
                 std::real(std::conj(square[i]) * rate[i]);
      return value;
    };
    const double conservativeEnergyRate =
        quadraticRate(conservativeRate) + quarticRate(conservativeSquareRate);
    if (std::abs(conservativeEnergyRate) > 2.e-11)
      throw std::runtime_error(
          "projected Hamiltonian transfer does not close: " +
          std::to_string(conservativeEnergyRate));

    const auto projectedEnergy = [&](const SpectralField &state) {
      SpectralField stateSquare, unused;
      transform.projectedSquareSpectra(state, zeroRate, stateSquare, unused);
      double value = 0.0;
      for (std::size_t y = 0; y < p.ny; ++y)
        for (std::size_t x = 0; x < p.nx; ++x) {
          const std::size_t i = spectralIndex(x, y, p.nx);
          const double weight =
              -p.dispersionCoefficient * waveNumberSquared(p, x, y) +
              p.chemicalPotential;
          value += area * weight * std::norm(state[i]);
        }
      for (const Complex squareMode : stateSquare)
        value += 0.5 * area * p.nonlinearityCoefficient * std::norm(squareMode);
      return value;
    };
    constexpr double epsilon = 1.e-7;
    SpectralField plus = w, minus = w;
    for (std::size_t i = 0; i < w.size(); ++i) {
      plus[i] += epsilon * dampingRate[i];
      minus[i] -= epsilon * dampingRate[i];
    }
    const double finiteDifferenceLoss =
        -(projectedEnergy(plus) - projectedEnergy(minus)) / (2.0 * epsilon);
    const double analyticLoss =
        -quadraticRate(dampingRate) - quarticRate(dampingSquareRate);
    if (std::abs(finiteDifferenceLoss - analyticLoss) >
        2.e-8 * std::max(1.0, std::abs(analyticLoss)))
      throw std::runtime_error(
          "projected-Hamiltonian damping rate is inconsistent: " +
          std::to_string(finiteDifferenceLoss - analyticLoss));
    backend.reset();
    backendFinalize();
    std::cout << "dealiased cubic term and projected Hamiltonian agree\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
