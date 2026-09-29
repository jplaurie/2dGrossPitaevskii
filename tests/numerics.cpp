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

SpectralField directNonlinear(const Parameters &parameters, const SpectralField &w) {
    const std::size_t paddedCount = parameters.mx() * parameters.my();
    std::vector<Complex> psi(paddedCount), squareSpectrum(paddedCount), cubic(paddedCount);
    for (std::size_t py = 0; py < parameters.my(); ++py)
        for (std::size_t px = 0; px < parameters.mx(); ++px) {
            Complex sum{};
            for (std::size_t y = 0; y < parameters.ny; ++y)
                for (std::size_t x = 0; x < parameters.nx; ++x) {
                    const double angle =
                        2.0 * gpPi *
                        (static_cast<double>(signedWave(x, parameters.nx)) * px / parameters.mx() +
                         static_cast<double>(signedWave(y, parameters.ny)) * py / parameters.my());
                    sum += w[spectralIndex(x, y, parameters.nx)] * phase(angle);
                }
            psi[spectralIndex(px, py, parameters.mx())] = sum;
        }
    for (std::size_t ky = 0; ky < parameters.my(); ++ky)
        for (std::size_t kx = 0; kx < parameters.mx(); ++kx) {
            Complex sum{};
            for (std::size_t py = 0; py < parameters.my(); ++py)
                for (std::size_t px = 0; px < parameters.mx(); ++px) {
                    const double angle = -2.0 * gpPi *
                                         (static_cast<double>(signedWave(kx, parameters.mx())) *
                                              px / parameters.mx() +
                                          static_cast<double>(signedWave(ky, parameters.my())) *
                                              py / parameters.my());
                    const Complex value = psi[spectralIndex(px, py, parameters.mx())];
                    sum += value * value * phase(angle);
                }
            if (retainedPaddedWave(signedWave(kx, parameters.mx()), parameters.nx) &&
                retainedPaddedWave(signedWave(ky, parameters.my()), parameters.ny))
                squareSpectrum[spectralIndex(kx, ky, parameters.mx())] =
                    sum / static_cast<double>(paddedCount);
        }
    for (std::size_t py = 0; py < parameters.my(); ++py)
        for (std::size_t px = 0; px < parameters.mx(); ++px) {
            Complex sum{};
            for (std::size_t ky = 0; ky < parameters.my(); ++ky)
                for (std::size_t kx = 0; kx < parameters.mx(); ++kx) {
                    const double angle = 2.0 * gpPi *
                                         (static_cast<double>(signedWave(kx, parameters.mx())) *
                                              px / parameters.mx() +
                                          static_cast<double>(signedWave(ky, parameters.my())) *
                                              py / parameters.my());
                    sum += squareSpectrum[spectralIndex(kx, ky, parameters.mx())] * phase(angle);
                }
            const std::size_t i = spectralIndex(px, py, parameters.mx());
            cubic[i] = sum * std::conj(psi[i]);
        }
    SpectralField result(parameters.nx * parameters.ny);
    for (std::size_t y = 0; y < parameters.ny; ++y)
        for (std::size_t x = 0; x < parameters.nx; ++x) {
            Complex sum{};
            for (std::size_t py = 0; py < parameters.my(); ++py)
                for (std::size_t px = 0; px < parameters.mx(); ++px) {
                    const double angle =
                        -2.0 * gpPi *
                        (static_cast<double>(signedWave(x, parameters.nx)) * px / parameters.mx() +
                         static_cast<double>(signedWave(y, parameters.ny)) * py / parameters.my());
                    sum += cubic[spectralIndex(px, py, parameters.mx())] * phase(angle);
                }
            result[spectralIndex(x, y, parameters.nx)] = parameters.nonlinearityCoefficient * sum /
                                                         static_cast<double>(paddedCount) /
                                                         Complex(0.0, 1.0);
        }
    return result;
}
} // namespace

int main(int argc, char **argv) {
    try {
        backendInitialize(argc, argv);
        Parameters parameters;
        parameters.nx = 4;
        parameters.ny = 6; // Its 3/2-padded dimension is odd, exercising FFT indexing.
        parameters.nonlinearityCoefficient = 1.7;
        parameters.hyperviscosity = 0.0;
        parameters.hypoviscosity = 0.0;
        parameters.ginzburgLandauDamping = 0.0;
        parameters.threadCount = 1;
        if (signedWave(4, 9) != 4 || signedWave(5, 9) != -4)
            throw std::runtime_error("odd-length Fourier indexing is incorrect");
        SpectralField w(parameters.nx * parameters.ny);
        w[spectralIndex(0, 0, parameters.nx)] = {0.7, -0.2};
        w[spectralIndex(1, 0, parameters.nx)] = {0.12, 0.05};
        w[spectralIndex(0, 1, parameters.nx)] = {-0.08, 0.04};
        w[spectralIndex(1, 1, parameters.nx)] = {0.03, -0.07};
        w[spectralIndex(3, 2, parameters.nx)] = {0.02, 0.01};
        auto backend = makeBackend(parameters);
        SpectralField actual;
        backend->evaluate(w, actual);
        const SpectralField expected = directNonlinear(parameters, w);
        double error = 0.0, magnitude = 0.0;
        for (std::size_t i = 0; i < actual.size(); ++i) {
            error = std::max(error, std::abs(actual[i] - expected[i]));
            magnitude = std::max(magnitude, std::abs(expected[i]));
        }
        if (error > 2.e-11 * std::max(1.0, magnitude))
            throw std::runtime_error("dealiased cubic term disagrees with direct DFT: " +
                                     std::to_string(error));

        BaseTransform transform(parameters);
        SpectralField conservativeRate(w.size()), dampingRate(w.size()), zeroRate(w.size());
        for (std::size_t y = 0; y < parameters.ny; ++y)
            for (std::size_t x = 0; x < parameters.nx; ++x) {
                const std::size_t i = spectralIndex(x, y, parameters.nx);
                const double k2 = waveNumberSquared(parameters, x, y);
                const double weight =
                    -parameters.dispersionCoefficient * k2 + parameters.chemicalPotential;
                conservativeRate[i] = Complex(0.0, -weight) * w[i] + actual[i];
                dampingRate[i] = -(0.2 + 0.01 * k2) * w[i];
            }
        SpectralField square, conservativeSquareRate;
        transform.projectedSquareSpectra(w, conservativeRate, square, conservativeSquareRate);
        SpectralField dampingSquareRate;
        transform.projectedSquareSpectra(w, dampingRate, square, dampingSquareRate);
        const double area = parameters.lx() * parameters.ly();
        const auto quadraticRate = [&](const SpectralField &rate) {
            double value = 0.0;
            for (std::size_t y = 0; y < parameters.ny; ++y)
                for (std::size_t x = 0; x < parameters.nx; ++x) {
                    const std::size_t i = spectralIndex(x, y, parameters.nx);
                    const double weight =
                        -parameters.dispersionCoefficient * waveNumberSquared(parameters, x, y) +
                        parameters.chemicalPotential;
                    value += 2.0 * area * weight * std::real(std::conj(w[i]) * rate[i]);
                }
            return value;
        };
        const auto quarticRate = [&](const SpectralField &rate) {
            double value = 0.0;
            for (std::size_t i = 0; i < square.size(); ++i)
                value += area * parameters.nonlinearityCoefficient *
                         std::real(std::conj(square[i]) * rate[i]);
            return value;
        };
        const double conservativeEnergyRate =
            quadraticRate(conservativeRate) + quarticRate(conservativeSquareRate);
        if (std::abs(conservativeEnergyRate) > 2.e-11)
            throw std::runtime_error("projected Hamiltonian transfer does not close: " +
                                     std::to_string(conservativeEnergyRate));

        const auto projectedEnergy = [&](const SpectralField &state) {
            SpectralField stateSquare, unused;
            transform.projectedSquareSpectra(state, zeroRate, stateSquare, unused);
            double value = 0.0;
            for (std::size_t y = 0; y < parameters.ny; ++y)
                for (std::size_t x = 0; x < parameters.nx; ++x) {
                    const std::size_t i = spectralIndex(x, y, parameters.nx);
                    const double weight =
                        -parameters.dispersionCoefficient * waveNumberSquared(parameters, x, y) +
                        parameters.chemicalPotential;
                    value += area * weight * std::norm(state[i]);
                }
            for (const Complex squareMode : stateSquare)
                value += 0.5 * area * parameters.nonlinearityCoefficient * std::norm(squareMode);
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
        const double analyticLoss = -quadraticRate(dampingRate) - quarticRate(dampingSquareRate);
        if (std::abs(finiteDifferenceLoss - analyticLoss) >
            2.e-8 * std::max(1.0, std::abs(analyticLoss)))
            throw std::runtime_error("projected-Hamiltonian damping rate is inconsistent: " +
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
