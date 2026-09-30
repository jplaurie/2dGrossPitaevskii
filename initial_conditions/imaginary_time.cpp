#include "backend.hpp"
#include "fftw_utils.hpp"
#include "integrator.hpp"
#include "parallel.hpp"
#include "spectral.hpp"
#include "vortex_field.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {
using Complex = std::complex<double>;

struct CommandLine {
    std::filesystem::path parameterFile;
    std::filesystem::path inputFile;
    std::filesystem::path outputFile;
    std::filesystem::path diagnosticsFile;
    double driftX = 0.0;
    double driftY = 0.0;
    double tolerance = 1.0e-10;
    std::uint64_t minimumSteps = 0;
    bool overwrite = false;
};

void usage(std::ostream &out) {
    out << "Usage: gp2d_relax --parameters FILE [--input FILE] --output FILE [options]\n"
           "\n"
           "Relax a GP wavefunction by CPU imaginary-time integration. Grid, physical,\n"
           "time-step, integrator, step-count, output interval, and thread settings come\n"
           "from the normal solver parameter file. --input defaults to its\n"
           "initialConditionFile. Forcing and real-time dissipation settings are ignored.\n"
           "\n"
           "Options:\n"
           "  --drift-x U          comoving-frame x velocity (default 0)\n"
           "  --drift-y U          comoving-frame y velocity (default 0)\n"
           "  --tolerance EPS      relative stationary residual; 0 disables early stop\n"
           "                       (default 1e-10)\n"
           "  --minimum-steps N    do not test convergence before N steps (default 0)\n"
           "  --diagnostics FILE   convergence CSV (default OUTPUT.relaxation.csv)\n"
           "  --overwrite          replace existing output and diagnostics files\n"
           "  -h, --help           show this help\n";
}

const char *nextArgument(int &index, int argc, char **argv) {
    if (++index >= argc)
        throw std::runtime_error(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

template <class T> T number(const char *text, const std::string &option) {
    const std::string value(text);
    std::size_t used = 0;
    try {
        if constexpr (std::is_same_v<T, double>) {
            const double parsed = std::stod(value, &used);
            if (used != value.size() || !std::isfinite(parsed))
                throw std::invalid_argument("invalid");
            return parsed;
        } else {
            if (!value.empty() && value.front() == '-')
                throw std::invalid_argument("negative");
            const unsigned long long parsed = std::stoull(value, &used);
            if (used != value.size())
                throw std::invalid_argument("invalid");
            return static_cast<T>(parsed);
        }
    } catch (const std::exception &) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
}

CommandLine parseCommandLine(int argc, char **argv) {
    CommandLine command;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "-h" || argument == "--help") {
            usage(std::cout);
            std::exit(0);
        } else if (argument == "--parameters")
            command.parameterFile = nextArgument(i, argc, argv);
        else if (argument == "--input")
            command.inputFile = nextArgument(i, argc, argv);
        else if (argument == "--output")
            command.outputFile = nextArgument(i, argc, argv);
        else if (argument == "--diagnostics")
            command.diagnosticsFile = nextArgument(i, argc, argv);
        else if (argument == "--drift-x")
            command.driftX = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--drift-y")
            command.driftY = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--tolerance")
            command.tolerance = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--minimum-steps")
            command.minimumSteps = number<std::uint64_t>(nextArgument(i, argc, argv), argument);
        else if (argument == "--overwrite")
            command.overwrite = true;
        else
            throw std::runtime_error("unknown option: " + argument);
    }
    if (command.parameterFile.empty() || command.outputFile.empty())
        throw std::runtime_error("--parameters and --output are required");
    if (command.tolerance < 0.0)
        throw std::runtime_error("--tolerance cannot be negative");
    if (command.diagnosticsFile.empty())
        command.diagnosticsFile = command.outputFile.string() + ".relaxation.csv";
    return command;
}

Complex phi(Complex z, int order) {
    double factorial = 1.0;
    for (int k = 2; k <= order; ++k)
        factorial *= k;
    if (std::abs(z) < 2.0) {
        Complex term = 1.0 / factorial;
        Complex sum = term;
        for (int k = 1; k < 64; ++k) {
            term *= z / static_cast<double>(order + k);
            sum += term;
            if (std::abs(term) < 1.0e-17 * std::abs(sum))
                break;
        }
        return sum;
    }
    Complex value = (std::exp(z) - 1.0) / z;
    double previousFactorial = 1.0;
    for (int k = 2; k <= order; ++k) {
        value = (value - 1.0 / previousFactorial) / z;
        previousFactorial *= k;
    }
    return value;
}

IntegrationCoefficients buildCoefficients(const Parameters &parameters,
                                          const std::vector<Complex> &linearOperator) {
    IntegrationCoefficients result;
    const std::size_t count = linearOperator.size();
    result.e1.resize(count);
    if (parameters.usesEtd()) {
        result.q1.resize(count);
        result.f1.resize(count);
        if (parameters.integrator != Integrator::etd2) {
            result.e2.resize(count);
            result.q2.resize(count);
            result.f2.resize(count);
            result.f3.resize(count);
        }
        if (parameters.integrator == Integrator::etd4) {
            result.q3.resize(count);
            result.q4.resize(count);
            result.q5.resize(count);
        }
    }
    forEachIndex(count, [&](std::size_t i) {
        const Complex z = parameters.timeStep * linearOperator[i];
        result.e1[i] = std::exp(z);
        if (!parameters.usesEtd())
            return;
        const Complex p1 = phi(z, 1), p2 = phi(z, 2);
        if (parameters.integrator == Integrator::etd2) {
            result.q1[i] = parameters.timeStep * p1;
            result.f1[i] = parameters.timeStep * p2;
            return;
        }
        const Complex p3 = phi(z, 3), halfP1 = phi(0.5 * z, 1);
        result.e2[i] = std::exp(0.5 * z);
        result.q1[i] = 0.5 * parameters.timeStep * halfP1;
        result.q2[i] = parameters.timeStep * p1;
        result.f1[i] = parameters.timeStep * (p1 - 3.0 * p2 + 4.0 * p3);
        result.f2[i] = parameters.timeStep * (p2 - 2.0 * p3);
        result.f3[i] = parameters.timeStep * (-p2 + 4.0 * p3);
        if (parameters.integrator == Integrator::etd4) {
            const Complex halfP2 = phi(0.5 * z, 2);
            result.q2[i] = parameters.timeStep * (0.5 * halfP1 - halfP2);
            result.q3[i] = parameters.timeStep * halfP2;
            result.q4[i] = parameters.timeStep * (p1 - 2.0 * p2);
            result.q5[i] = 2.0 * parameters.timeStep * p2;
        }
    });
    for (const auto *field : result.fields())
        if (!std::all_of(field->begin(), field->end(), [](Complex value) {
                return std::isfinite(value.real()) && std::isfinite(value.imag());
            }))
            throw std::runtime_error("non-finite imaginary-time integration coefficient");
    return result;
}

class Relaxer {
  public:
    Relaxer(Parameters parameters, double driftX, double driftY,
            std::unique_ptr<NonlinearBackend> backend)
        : parameters_(std::move(parameters)), driftX_(driftX), driftY_(driftY),
          backend_(std::move(backend)), transform_(parameters_),
          linearOperator_(parameters_.nx * parameters_.ny),
          stages_{SpectralField(linearOperator_.size()), SpectralField(linearOperator_.size()),
                  SpectralField(linearOperator_.size()), SpectralField(linearOperator_.size())},
          states_{SpectralField(linearOperator_.size()), SpectralField(linearOperator_.size()),
                  SpectralField(linearOperator_.size())},
          work_(linearOperator_.size()), zeroRate_(linearOperator_.size()) {
        for (std::size_t y = 0; y < parameters_.ny; ++y)
            for (std::size_t x = 0; x < parameters_.nx; ++x) {
                const std::size_t i = spectralIndex(x, y, parameters_.nx);
                const double kx = 2.0 * gpPi * static_cast<double>(signedWave(x, parameters_.nx)) /
                                  parameters_.lx();
                const double ky = 2.0 * gpPi * static_cast<double>(signedWave(y, parameters_.ny)) /
                                  parameters_.ly();
                const double k2 = kx * kx + ky * ky;
                linearOperator_[i] = parameters_.dispersionCoefficient * k2 -
                                     parameters_.chemicalPotential + driftX_ * kx + driftY_ * ky;
            }
        coefficients_ = buildCoefficients(parameters_, linearOperator_);
    }

    BaseTransform &transform() { return transform_; }

    void stationaryResidual(const SpectralField &wavefunction, SpectralField &residual) {
        nonlinear(wavefunction, residual);
        forEachIndex(residual.size(),
                     [&](std::size_t i) { residual[i] += linearOperator_[i] * wavefunction[i]; });
    }

    void step(SpectralField &wavefunction) {
        const auto c = coefficients_.pointers();
        nonlinear(wavefunction, stages_[0]);
        forEachIndex(wavefunction.size(), [&](std::size_t i) {
            states_[0][i] = integrationStageA(parameters_.integrator, parameters_.timeStep, i, c,
                                              wavefunction[i], stages_[0][i]);
        });
        nonlinear(states_[0], stages_[1]);
        if (parameters_.integrator == Integrator::etd3 ||
            parameters_.integrator == Integrator::etd4) {
            forEachIndex(wavefunction.size(), [&](std::size_t i) {
                states_[1][i] = integrationStageB(parameters_.integrator, i, c, wavefunction[i],
                                                  stages_[0][i], stages_[1][i]);
            });
            nonlinear(states_[1], stages_[2]);
        }
        if (parameters_.integrator == Integrator::etd4) {
            forEachIndex(wavefunction.size(), [&](std::size_t i) {
                states_[2][i] =
                    integrationStageC(i, c, wavefunction[i], stages_[0][i], stages_[2][i]);
            });
            nonlinear(states_[2], stages_[3]);
        }
        forEachIndex(wavefunction.size(), [&](std::size_t i) {
            wavefunction[i] = integrationFinish(parameters_.integrator, parameters_.timeStep, i, c,
                                                wavefunction[i], states_[0][i], stages_[0][i],
                                                stages_[1][i], stages_[2][i], stages_[3][i]);
        });
    }

    struct Diagnostic {
        double comovingEnergy = 0.0;
        double waveAction = 0.0;
        double relativeResidual = 0.0;
    };

    Diagnostic diagnostic(const SpectralField &wavefunction) {
        stationaryResidual(wavefunction, work_);
        double residual2 = 0.0, state2 = 0.0, quadratic = 0.0;
        for (std::size_t y = 0; y < parameters_.ny; ++y)
            for (std::size_t x = 0; x < parameters_.nx; ++x) {
                const std::size_t i = spectralIndex(x, y, parameters_.nx);
                const double kx = 2.0 * gpPi * static_cast<double>(signedWave(x, parameters_.nx)) /
                                  parameters_.lx();
                const double ky = 2.0 * gpPi * static_cast<double>(signedWave(y, parameters_.ny)) /
                                  parameters_.ly();
                const double k2 = kx * kx + ky * ky;
                residual2 += std::norm(work_[i]);
                state2 += std::norm(wavefunction[i]);
                quadratic += (-parameters_.dispersionCoefficient * k2 +
                              parameters_.chemicalPotential - driftX_ * kx - driftY_ * ky) *
                             std::norm(wavefunction[i]);
            }
        SpectralField square, unused;
        transform_.projectedSquareSpectra(wavefunction, zeroRate_, square, unused);
        double quartic = 0.0;
        for (const Complex value : square)
            quartic += std::norm(value);
        const double area = parameters_.lx() * parameters_.ly();
        return {area * (quadratic + 0.5 * parameters_.nonlinearityCoefficient * quartic),
                area * state2,
                std::sqrt(residual2 / std::max(state2, std::numeric_limits<double>::min()))};
    }

  private:
    void nonlinear(const SpectralField &wavefunction, SpectralField &result) {
        backend_->evaluate(wavefunction, result);
        // The real-time CPU backend returns g|psi|^2 psi / i = -i g|psi|^2 psi.
        // Multiplication by -i converts it to the imaginary-time term -g|psi|^2 psi.
        forEachIndex(result.size(), [&](std::size_t i) { result[i] *= Complex(0.0, -1.0); });
    }

    Parameters parameters_;
    double driftX_, driftY_;
    std::unique_ptr<NonlinearBackend> backend_;
    BaseTransform transform_;
    SpectralField linearOperator_;
    IntegrationCoefficients coefficients_;
    std::array<SpectralField, 4> stages_;
    std::array<SpectralField, 3> states_;
    SpectralField work_, zeroRate_;
};

void ensureWritable(const std::filesystem::path &path, bool overwrite) {
    if (std::filesystem::exists(path) && !overwrite)
        throw std::runtime_error("refusing to overwrite file: " + path.string());
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
}
} // namespace

int main(int argc, char **argv) {
    bool backendReady = false;
    try {
        const CommandLine command = parseCommandLine(argc, argv);
        Parameters parameters = readParameters(command.parameterFile, false);
        const std::filesystem::path input =
            command.inputFile.empty() ? parameters.initialConditionFile : command.inputFile;
        if (input.empty())
            throw std::runtime_error(
                "provide --input or initialConditionFile in the parameter file");
        if (command.outputFile.lexically_normal() == command.diagnosticsFile.lexically_normal())
            throw std::runtime_error(
                "wavefunction output and diagnostics must use different files");
        ensureWritable(command.outputFile, command.overwrite);
        ensureWritable(command.diagnosticsFile, command.overwrite);

        backendInitialize(argc, argv);
        backendReady = true;
        {
            parameters.ginzburgLandauDamping = 0.0;
            auto backend = makeBackend(parameters);
            Relaxer relaxer(parameters, command.driftX, command.driftY, std::move(backend));
            auto physical = readWavefunction(input, parameters);
            SpectralField wavefunction;
            relaxer.transform().forward(physical, wavefunction);

            std::ofstream diagnostics(command.diagnosticsFile, std::ios::trunc);
            if (!diagnostics)
                throw std::runtime_error("cannot write diagnostics: " +
                                         command.diagnosticsFile.string());
            diagnostics << std::scientific << std::setprecision(12)
                        << "imaginary_time,step,comoving_energy,wave_action,relative_residual\n";

            const auto report = [&](std::uint64_t step) {
                const Relaxer::Diagnostic values = relaxer.diagnostic(wavefunction);
                const double time = parameters.timeStep * static_cast<double>(step);
                diagnostics << time << ',' << step << ',' << values.comovingEnergy << ','
                            << values.waveAction << ',' << values.relativeResidual << '\n';
                diagnostics.flush();
                std::cout << "imaginary time = " << time << " step = " << step
                          << " comoving energy = " << values.comovingEnergy
                          << " residual = " << values.relativeResidual << '\n';
                return values.relativeResidual;
            };

            double residual = report(0);
            std::uint64_t completedSteps = 0;
            for (std::uint64_t iteration = 0; iteration < parameters.numberOfSteps; ++iteration) {
                const std::uint64_t step = iteration + 1;
                relaxer.step(wavefunction);
                completedSteps = step;
                const bool scheduled =
                    step % parameters.outputIntervalSteps == 0 || step == parameters.numberOfSteps;
                if (scheduled) {
                    residual = report(step);
                    if (command.tolerance > 0.0 && step >= command.minimumSteps &&
                        residual <= command.tolerance)
                        break;
                }
            }
            diagnostics.close();
            if (!diagnostics)
                throw std::runtime_error("failed while writing diagnostics: " +
                                         command.diagnosticsFile.string());

            relaxer.transform().inverse(wavefunction, physical);
            writeWavefunctionFile(command.outputFile, physical, parameters, command.overwrite);
            std::cout << "wrote relaxed wavefunction to " << command.outputFile << " after "
                      << completedSteps << " steps; final residual = " << residual << '\n';
        }
        backendFinalize();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        if (backendReady)
            backendFinalize();
        return 1;
    }
}
