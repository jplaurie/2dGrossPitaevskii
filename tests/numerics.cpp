#include "backend.hpp"
#include "fftw_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace {
// Exact, non-wrapping Fourier convolutions: independent of FFT padding.
std::optional<std::size_t> modeIndex(const Parameters &p, long x, long y) {
    if (!retainedPaddedWave(x, p.nx) || !retainedPaddedWave(y, p.ny))
        return std::nullopt;
    return spectralIndex(x < 0 ? p.nx + x : x, y < 0 ? p.ny + y : y, p.nx);
}

SpectralField directDensity(const Parameters &p, const SpectralField &w) {
    SpectralField density(w.size());
    for (std::size_t d = 0; d < w.size(); ++d) {
        const long dx = signedWave(d % p.nx, p.nx), dy = signedWave(d / p.nx, p.ny);
        if (!retainedDensityWave(dx, p.nx) || !retainedDensityWave(dy, p.ny))
            continue;
        for (std::size_t j = 0; j < w.size(); ++j) {
            const auto other =
                modeIndex(p, signedWave(j % p.nx, p.nx) + dx, signedWave(j / p.nx, p.ny) + dy);
            if (other)
                density[d] += w[*other] * std::conj(w[j]);
        }
    }
    return density;
}

SpectralField directNonlinear(const Parameters &p, const SpectralField &w) {
    const auto density = directDensity(p, w);
    SpectralField result(w.size());
    for (std::size_t k = 0; k < w.size(); ++k)
        for (std::size_t d = 0; d < w.size(); ++d) {
            const auto other = modeIndex(p, signedWave(k % p.nx, p.nx) - signedWave(d % p.nx, p.nx),
                                         signedWave(k / p.nx, p.ny) - signedWave(d / p.nx, p.ny));
            if (other)
                result[k] += p.nonlinearityCoefficient * density[d] * w[*other] / Complex(0., 1.);
        }
    return result;
}

void checkMomentumAndPlaneWaves(const Parameters &p, NonlinearBackend &backend) {
    SpectralField w(p.nx * p.ny), rate;
    for (std::size_t i = 0; i < w.size(); ++i)
        w[i] = {0.12 * std::cos(0.71 * i), 0.09 * std::sin(0.39 * i)};
    backend.evaluate(w, rate);
    const auto exact = directNonlinear(p, w);
    double massRate = 0., momentumXRate = 0., momentumYRate = 0.;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (std::abs(rate[i] - exact[i]) > 1.e-12)
            throw std::runtime_error("density convolution disagrees at a spectral boundary");
        const double transfer = 2. * std::real(std::conj(w[i]) * rate[i]);
        massRate += transfer;
        momentumXRate += signedWave(i % p.nx, p.nx) * transfer;
        momentumYRate += signedWave(i / p.nx, p.ny) * transfer;
    }
    if (std::max({std::abs(massRate), std::abs(momentumXRate), std::abs(momentumYRate)}) > 1.e-12)
        throw std::runtime_error("density-projected nonlinearity does not conserve mass/momentum");
    for (std::size_t mode = 0; mode < w.size(); ++mode) {
        std::fill(w.begin(), w.end(), Complex{});
        w[mode] = {0.7, 0.2};
        backend.evaluate(w, rate);
        for (std::size_t i = 0; i < w.size(); ++i) {
            const Complex expected =
                p.nonlinearityCoefficient * std::norm(w[mode]) * w[i] / Complex(0., 1.);
            if (std::abs(rate[i] - expected) > 1.e-12)
                throw std::runtime_error("plane-wave nonlinear phase is incorrect");
        }
    }
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
        checkMomentumAndPlaneWaves(parameters, *backend);
        {
            Parameters edgeParameters = parameters;
            edgeParameters.nx = edgeParameters.ny = 8;
            auto edgeBackend = makeBackend(edgeParameters);
            checkMomentumAndPlaneWaves(edgeParameters, *edgeBackend);
        }
        SpectralField actual;
        backend->evaluate(w, actual);
        const SpectralField expected = directNonlinear(parameters, w);
        double error = 0.0, magnitude = 0.0;
        for (std::size_t i = 0; i < actual.size(); ++i) {
            error = std::max(error, std::abs(actual[i] - expected[i]));
            magnitude = std::max(magnitude, std::abs(expected[i]));
        }
        if (error > 2.e-11 * std::max(1.0, magnitude))
            throw std::runtime_error(
                "dealiased cubic term disagrees with exact Fourier convolution: " +
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
        SpectralField density, conservativeDensityRate;
        transform.projectedDensitySpectra(w, conservativeRate, density, conservativeDensityRate);
        SpectralField dampingDensityRate;
        transform.projectedDensitySpectra(w, dampingRate, density, dampingDensityRate);
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
            for (std::size_t i = 0; i < density.size(); ++i)
                value += area * parameters.nonlinearityCoefficient *
                         std::real(std::conj(density[i]) * rate[i]);
            return value;
        };
        const double conservativeEnergyRate =
            quadraticRate(conservativeRate) + quarticRate(conservativeDensityRate);
        if (std::abs(conservativeEnergyRate) > 2.e-11)
            throw std::runtime_error("projected Hamiltonian transfer does not close: " +
                                     std::to_string(conservativeEnergyRate));

        const auto projectedEnergy = [&](const SpectralField &state) {
            SpectralField stateDensity, unused;
            transform.projectedDensitySpectra(state, zeroRate, stateDensity, unused);
            double value = 0.0;
            for (std::size_t y = 0; y < parameters.ny; ++y)
                for (std::size_t x = 0; x < parameters.nx; ++x) {
                    const std::size_t i = spectralIndex(x, y, parameters.nx);
                    const double weight =
                        -parameters.dispersionCoefficient * waveNumberSquared(parameters, x, y) +
                        parameters.chemicalPotential;
                    value += area * weight * std::norm(state[i]);
                }
            for (const Complex densityMode : stateDensity)
                value += 0.5 * area * parameters.nonlinearityCoefficient * std::norm(densityMode);
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
        const double analyticLoss = -quadraticRate(dampingRate) - quarticRate(dampingDensityRate);
        if (std::abs(finiteDifferenceLoss - analyticLoss) >
            2.e-8 * std::max(1.0, std::abs(analyticLoss)))
            throw std::runtime_error("projected-Hamiltonian damping rate is inconsistent: " +
                                     std::to_string(finiteDifferenceLoss - analyticLoss));

        // Independently differentiate the Hamiltonian along both real noise
        // directions. Circular complex noise has variance f_j^2/2 per direction.
        constexpr double noiseEpsilon = 1.e-3;
        const double energy = projectedEnergy(w);
        double trace = 0., predictedTrace = 0.;
        for (std::size_t j = 0; j < w.size(); ++j) {
            const double power = .01 * (1. + .1 * j);
            double laplacian = 0.;
            for (const Complex direction : {Complex(1., 0.), Complex(0., 1.)}) {
                auto upper = w, lower = w;
                upper[j] += noiseEpsilon * direction;
                lower[j] -= noiseEpsilon * direction;
                laplacian += (projectedEnergy(upper) - 2. * energy + projectedEnergy(lower)) /
                             (noiseEpsilon * noiseEpsilon);
            }
            trace += .25 * power * laplacian;
            const double weight =
                -parameters.dispersionCoefficient *
                    waveNumberSquared(parameters, j % parameters.nx, j / parameters.nx) +
                parameters.chemicalPotential;
            predictedTrace += area * power * weight;
            for (std::size_t i = 0; i < w.size(); ++i) {
                const long dx = signedWave(j % parameters.nx, parameters.nx) -
                                signedWave(i % parameters.nx, parameters.nx);
                const long dy = signedWave(j / parameters.nx, parameters.ny) -
                                signedWave(i / parameters.nx, parameters.ny);
                const double accepted =
                    retainedDensityWave(dx, parameters.nx) && retainedDensityWave(dy, parameters.ny)
                        ? 1.
                        : 0.;
                predictedTrace += area * parameters.nonlinearityCoefficient * power *
                                  (1. + accepted) * std::norm(w[i]);
            }
        }
        if (std::abs(trace - predictedTrace) > 2.e-6 * std::max(1., std::abs(predictedTrace)))
            throw std::runtime_error("density Hamiltonian Ito trace is inconsistent");
        backend.reset();
        backendFinalize();
        std::cout << "density convolution, plane waves, momentum, and Hamiltonian agree\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
