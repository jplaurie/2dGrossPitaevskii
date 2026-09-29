#include "io_utils.hpp"
#include "output.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
std::ofstream numericOutput(const std::filesystem::path &path,
                            std::ios::openmode mode = std::ios::out) {
    std::ofstream output(path, mode);
    if (!output)
        throw std::runtime_error("cannot write output file: " + path.string());
    output << std::scientific << std::setprecision(12);
    return output;
}
} // namespace

double writeDiagnostics(const Parameters &parameters, BaseTransform &transform, double time,
                        std::uint64_t frame, const SpectralField &wavefunction,
                        const SpectralField &nonlinear, const std::vector<double> &forcingAmplitude,
                        const std::vector<double> &stochasticQuarticInjectionWeight,
                        DiagnosticsAverages &averages) {
    if (forcingAmplitude.size() != wavefunction.size() ||
        (!stochasticQuarticInjectionWeight.empty() &&
         stochasticQuarticInjectionWeight.size() != wavefunction.size()))
        throw std::runtime_error("invalid forcing diagnostics size");
    const std::size_t bins = parameters.spectrumBins();
    if (averages.waveActionSpectrum.empty()) {
        averages.waveActionSpectrum.assign(bins, 0.0);
        averages.quadraticEnergySpectrum.assign(bins, 0.0);
        averages.waveActionFlux.assign(bins, 0.0);
        averages.fullEnergyFlux.assign(bins, 0.0);
    }
    std::vector<double> waveSpectrum(bins), quadraticSpectrum(bins), waveShell(bins),
        quadraticShell(bins);
    SpectralField hamiltonianRate(wavefunction.size());
    const double area = parameters.lx() * parameters.ly();
    const double binWidth = std::min(2.0 * gpPi / parameters.lx(), 2.0 * gpPi / parameters.ly());
    double waveAction = 0.0, kineticEnergy = 0.0, potentialEnergy = 0.0;
    double waveHypo = 0.0, waveHyper = 0.0;
    double quadraticHypo = 0.0, quadraticHyper = 0.0;
    double totalEnergyHypo = 0.0, totalEnergyHyper = 0.0;
    double expectedFullEnergyInjection = 0.0;
    const bool stochasticForcing =
        parameters.forcingEnabled && parameters.forcingProfile != ForcingProfile::singleMode;
    const bool deterministicForcing =
        parameters.forcingEnabled && parameters.forcingProfile == ForcingProfile::singleMode;
    for (std::size_t y = 0; y < parameters.ny; ++y) {
        for (std::size_t x = 0; x < parameters.nx; ++x) {
            const std::size_t index = spectralIndex(x, y, parameters.nx);
            const double k2 = waveNumberSquared(parameters, x, y), k = std::sqrt(k2);
            const double modalWaveAction = std::norm(wavefunction[index]);
            const double kineticWeight = -parameters.dispersionCoefficient * k2;
            const double weight = kineticWeight + parameters.chemicalPotential;
            waveAction += area * modalWaveAction;
            kineticEnergy += area * kineticWeight * modalWaveAction;
            potentialEnergy += area * parameters.chemicalPotential * modalWaveAction;
            double hyperviscousRate = 0.0, hypoviscousRate = 0.0;
            if (parameters.hyperviscosity > 0.0 &&
                (!parameters.hyperviscosityCutoffEnabled || k > parameters.hyperviscosityCutoff))
                hyperviscousRate =
                    parameters.hyperviscosity * std::pow(k2, parameters.hyperviscosityOrder);
            if ((k2 > 0.0 || parameters.hypoviscosityOrder >= 0.0) &&
                parameters.hypoviscosity > 0.0 &&
                (!parameters.hypoviscosityCutoffEnabled || k < parameters.hypoviscosityCutoff))
                hypoviscousRate =
                    parameters.hypoviscosity * std::pow(k2, parameters.hypoviscosityOrder);
            waveHyper += 2.0 * area * hyperviscousRate * modalWaveAction;
            waveHypo += 2.0 * area * hypoviscousRate * modalWaveAction;
            quadraticHyper += 2.0 * area * weight * hyperviscousRate * modalWaveAction;
            quadraticHypo += 2.0 * area * weight * hypoviscousRate * modalWaveAction;
            const std::size_t bin = static_cast<std::size_t>(k / binWidth);
            if (bin >= bins)
                throw std::logic_error("spectrum bin count is too small");
            waveSpectrum[bin] += area * modalWaveAction;
            quadraticSpectrum[bin] += area * weight * modalWaveAction;
            const double gamma = ginzburgLandauFactor(parameters, k2);
            const Complex conservativeNonlinear = Complex(1.0, gamma) * nonlinear[index];
            hamiltonianRate[index] =
                Complex(0.0, -weight) * wavefunction[index] + conservativeNonlinear;
            const Complex hamiltonianGradient =
                weight * wavefunction[index] + Complex(0.0, 1.0) * conservativeNonlinear;
            if (stochasticForcing) {
                const double amplitude = forcingAmplitude[index];
                expectedFullEnergyInjection += area * weight * amplitude * amplitude;
                if (!stochasticQuarticInjectionWeight.empty())
                    expectedFullEnergyInjection +=
                        area * stochasticQuarticInjectionWeight[index] * modalWaveAction;
            } else if (deterministicForcing) {
                expectedFullEnergyInjection +=
                    2.0 * area *
                    std::real(std::conj(hamiltonianGradient) * forcingAmplitude[index]);
            }
            totalEnergyHyper += 2.0 * area * hyperviscousRate *
                                std::real(std::conj(hamiltonianGradient) * wavefunction[index]);
            totalEnergyHypo += 2.0 * area * hypoviscousRate *
                               std::real(std::conj(hamiltonianGradient) * wavefunction[index]);
            const double transfer =
                -2.0 * area * std::real(std::conj(wavefunction[index]) * conservativeNonlinear);
            waveShell[bin] += transfer;
            quadraticShell[bin] += weight * transfer;
        }
    }
    SpectralField square, conservativeSquareRate;
    transform.projectedSquareSpectra(wavefunction, hamiltonianRate, square, conservativeSquareRate);
    // The two-pass cubic backend is the gradient of
    // (g/2) * ||P(psi^2)||^2, where P retains the base spectral band.
    // Using base-grid |psi|^4 here would describe a different Hamiltonian.
    double nonlinearEnergy = 0.0;
    for (const Complex value : square)
        nonlinearEnergy += std::norm(value);
    nonlinearEnergy *= 0.5 * area * parameters.nonlinearityCoefficient;
    const double totalEnergy = kineticEnergy + potentialEnergy + nonlinearEnergy;
    std::vector<double> fullEnergyShell = quadraticShell;
    for (std::size_t y = 0; y < parameters.ny; ++y)
        for (std::size_t x = 0; x < parameters.nx; ++x) {
            const std::size_t index = spectralIndex(x, y, parameters.nx);
            const std::size_t bin =
                static_cast<std::size_t>(std::sqrt(waveNumberSquared(parameters, x, y)) / binWidth);
            fullEnergyShell[bin] -=
                area * parameters.nonlinearityCoefficient *
                std::real(std::conj(square[index]) * conservativeSquareRate[index]);
        }
    std::vector<double> waveFlux(bins), fullEnergyFlux(bins);
    double cumulativeWave = 0.0, cumulativeEnergy = 0.0;
    for (std::size_t i = 0; i < bins; ++i) {
        cumulativeWave += waveShell[i];
        cumulativeEnergy += fullEnergyShell[i];
        waveFlux[i] = cumulativeWave;
        fullEnergyFlux[i] = cumulativeEnergy;
    }
    const auto finiteVector = [](const std::vector<double> &values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
    };
    if (!std::isfinite(totalEnergy) || !std::isfinite(waveAction) ||
        !std::isfinite(totalEnergyHypo) || !std::isfinite(totalEnergyHyper) ||
        !std::isfinite(expectedFullEnergyInjection) || !finiteVector(waveSpectrum) ||
        !finiteVector(quadraticSpectrum) || !finiteVector(waveFlux) ||
        !finiteVector(fullEnergyFlux))
        throw std::runtime_error(
            "diagnostics became non-finite; reduce the time step or coefficients");
    ++averages.count;
    for (std::size_t i = 0; i < bins; ++i) {
        averages.waveActionSpectrum[i] += waveSpectrum[i];
        averages.quadraticEnergySpectrum[i] += quadraticSpectrum[i];
        averages.waveActionFlux[i] += waveFlux[i];
        averages.fullEnergyFlux[i] += fullEnergyFlux[i];
    }
    auto diagnostics = numericOutput(parameters.outputDirectory / "diagnostics.csv",
                                     std::ios::out | std::ios::app);
    diagnostics << time << ',' << frame << ',' << totalEnergy << ',' << kineticEnergy << ','
                << potentialEnergy << ',' << nonlinearEnergy << ',' << waveAction << ',' << waveHypo
                << ',' << waveHyper << ',' << quadraticHypo << ',' << quadraticHyper << ','
                << totalEnergyHypo << ',' << totalEnergyHyper << ',' << expectedFullEnergyInjection
                << '\n';
    closeChecked(diagnostics, "failed while writing diagnostics.csv");
    auto spectra =
        numericOutput(parameters.outputDirectory / "spectra.csv", std::ios::out | std::ios::app);
    auto fluxes =
        numericOutput(parameters.outputDirectory / "fluxes.csv", std::ios::out | std::ios::app);
    for (std::size_t i = 0; i < bins; ++i) {
        const double k = static_cast<double>(i) * binWidth;
        spectra << time << ',' << frame << ',' << k << ',' << waveSpectrum[i] << ','
                << quadraticSpectrum[i] << ','
                << averages.waveActionSpectrum[i] / static_cast<double>(averages.count) << ','
                << averages.quadraticEnergySpectrum[i] / static_cast<double>(averages.count)
                << '\n';
        fluxes << time << ',' << frame << ',' << k << ',' << waveFlux[i] << ',' << fullEnergyFlux[i]
               << ',' << averages.waveActionFlux[i] / static_cast<double>(averages.count) << ','
               << averages.fullEnergyFlux[i] / static_cast<double>(averages.count) << '\n';
    }
    closeChecked(spectra, "failed while writing spectra.csv");
    closeChecked(fluxes, "failed while writing fluxes.csv");
    if (parameters.writeModeDiagnostics) {
        auto modes =
            numericOutput(parameters.outputDirectory / "modes.csv", std::ios::out | std::ios::app);
        modes << time << ',' << frame;
        for (const auto [x, y] : {std::pair{1UL, 0UL}, std::pair{0UL, 1UL}, std::pair{1UL, 1UL},
                                  std::pair{2UL, 1UL}, std::pair{0UL, 3UL}}) {
            if (x >= parameters.nx || y >= parameters.ny)
                modes << ",,";
            else {
                const Complex value = wavefunction[spectralIndex(x, y, parameters.nx)];
                modes << ',' << value.real() << ',' << value.imag();
            }
        }
        modes << '\n';
        closeChecked(modes, "failed while writing modes.csv");
    }
    return totalEnergy;
}

void writeForcingFiles(const Parameters &parameters, const std::vector<double> &amplitude,
                       std::size_t forcedModes, double waveCoefficient,
                       double quadraticEnergyCoefficient) {
    std::ofstream summary(parameters.outputDirectory / "forcing_summary.csv", std::ios::trunc);
    summary << "enabled,profile,temporal_type,forced_modes,wave_action_injection_"
               "coefficient,quadratic_energy_injection_coefficient\n"
            << std::boolalpha << parameters.forcingEnabled << ','
            << forcingProfileName(parameters.forcingProfile) << ',';
    if (!parameters.forcingEnabled)
        summary << "disabled";
    else if (parameters.forcingProfile == ForcingProfile::singleMode)
        summary << "deterministic";
    else
        summary << "stochastic";
    summary << ',' << forcedModes << ',';
    if (parameters.forcingEnabled && parameters.forcingProfile != ForcingProfile::singleMode)
        summary << std::setprecision(17) << waveCoefficient << ',' << quadraticEnergyCoefficient;
    else
        summary << ',';
    summary << '\n';
    closeChecked(summary, "failed while writing forcing_summary.csv");
    auto spectrum = numericOutput(parameters.outputDirectory / "forcing_spectrum.csv");
    spectrum << "kx,ky,amplitude\n";
    if (parameters.forcingEnabled)
        for (std::size_t y = 0; y < parameters.ny; ++y)
            for (std::size_t x = 0; x < parameters.nx; ++x) {
                const double kx = 2.0 * gpPi * static_cast<double>(signedWave(x, parameters.nx)) /
                                  parameters.lx();
                const double ky = 2.0 * gpPi * static_cast<double>(signedWave(y, parameters.ny)) /
                                  parameters.ly();
                spectrum << kx << ',' << ky << ',' << amplitude[spectralIndex(x, y, parameters.nx)]
                         << '\n';
            }
    closeChecked(spectrum, "failed while writing forcing_spectrum.csv");
}
