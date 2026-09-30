#include "vortex_field.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {
struct CommandLine {
    std::filesystem::path parameterFile;
    std::filesystem::path vortexFile;
    std::filesystem::path outputFile;
    std::optional<std::uint64_t> frame;
    std::optional<double> backgroundDensity;
    std::optional<double> healingLength;
    double coordinateScale = 1.0;
    int phaseWindingX = 0;
    int phaseWindingY = 0;
    bool overwrite = false;
};

void usage(std::ostream &out) {
    out << "Usage: gp2d_vortex_imprint --parameters FILE --vortices FILE --output FILE [options]\n"
           "\n"
           "Generate a periodic Gross--Pitaevskii wavefunction from point-vortex positions.\n"
           "Current 2dPointVortex initial-condition and trajectory.csv files are accepted.\n"
           "Circulation magnitudes are mapped to singly quantized windings by sign.\n"
           "\n"
           "Options:\n"
           "  --frame N                 trajectory.csv frame (default: last frame)\n"
           "  --background-density RHO override the automatic value -mu/g\n"
           "  --healing-length XI       override the automatic value sqrt(c/mu)\n"
           "  --coordinate-scale S      multiply input positions before wrapping (default 1)\n"
           "  --winding-x N             add N whole-domain phase windings in x\n"
           "  --winding-y N             add N whole-domain phase windings in y\n"
           "  --overwrite               replace an existing output file\n"
           "  -h, --help                show this help\n";
}

const char *nextArgument(int &index, int argc, char **argv) {
    if (++index >= argc)
        throw std::runtime_error(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

template <class T> T number(const char *text, const std::string &option) {
    std::string value(text);
    std::size_t used = 0;
    try {
        if constexpr (std::is_same_v<T, double>) {
            const double parsed = std::stod(value, &used);
            if (used != value.size() || !std::isfinite(parsed))
                throw std::invalid_argument("trailing");
            return parsed;
        } else if constexpr (std::is_same_v<T, int>) {
            const long parsed = std::stol(value, &used);
            if (used != value.size() || parsed < std::numeric_limits<int>::min() ||
                parsed > std::numeric_limits<int>::max())
                throw std::invalid_argument("range");
            return static_cast<int>(parsed);
        } else {
            if (!value.empty() && value.front() == '-')
                throw std::invalid_argument("negative");
            const unsigned long long parsed = std::stoull(value, &used);
            if (used != value.size())
                throw std::invalid_argument("trailing");
            return static_cast<T>(parsed);
        }
    } catch (const std::exception &) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
}

CommandLine parseCommandLine(int argc, char **argv) {
    CommandLine options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "-h" || argument == "--help") {
            usage(std::cout);
            std::exit(0);
        } else if (argument == "--parameters")
            options.parameterFile = nextArgument(i, argc, argv);
        else if (argument == "--vortices")
            options.vortexFile = nextArgument(i, argc, argv);
        else if (argument == "--output")
            options.outputFile = nextArgument(i, argc, argv);
        else if (argument == "--frame")
            options.frame = number<std::uint64_t>(nextArgument(i, argc, argv), argument);
        else if (argument == "--background-density")
            options.backgroundDensity = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--healing-length")
            options.healingLength = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--coordinate-scale")
            options.coordinateScale = number<double>(nextArgument(i, argc, argv), argument);
        else if (argument == "--winding-x")
            options.phaseWindingX = number<int>(nextArgument(i, argc, argv), argument);
        else if (argument == "--winding-y")
            options.phaseWindingY = number<int>(nextArgument(i, argc, argv), argument);
        else if (argument == "--overwrite")
            options.overwrite = true;
        else
            throw std::runtime_error("unknown option: " + argument);
    }
    if (options.parameterFile.empty() || options.vortexFile.empty() || options.outputFile.empty())
        throw std::runtime_error("--parameters, --vortices, and --output are required");
    return options;
}
} // namespace

int main(int argc, char **argv) {
    try {
        const CommandLine command = parseCommandLine(argc, argv);
        const Parameters parameters = readParameters(command.parameterFile, false);
        VortexImprintOptions options;
        options.coordinateScale = command.coordinateScale;
        options.phaseWindingX = command.phaseWindingX;
        options.phaseWindingY = command.phaseWindingY;

        if (command.backgroundDensity)
            options.backgroundDensity = *command.backgroundDensity;
        else {
            if (!(parameters.nonlinearityCoefficient > 0.0) ||
                !(parameters.chemicalPotential < 0.0))
                throw std::runtime_error("automatic background density requires "
                                         "nonlinearityCoefficient > 0 and chemicalPotential < 0; "
                                         "supply --background-density otherwise");
            options.backgroundDensity =
                -parameters.chemicalPotential / parameters.nonlinearityCoefficient;
        }
        if (command.healingLength)
            options.healingLength = *command.healingLength;
        else {
            if (!(parameters.dispersionCoefficient < 0.0) || !(parameters.chemicalPotential < 0.0))
                throw std::runtime_error("automatic healing length requires "
                                         "dispersionCoefficient < 0 and chemicalPotential < 0; "
                                         "supply --healing-length otherwise");
            options.healingLength =
                std::sqrt(parameters.dispersionCoefficient / parameters.chemicalPotential);
        }

        const PointVortexMetadata metadata = readPointVortexMetadata(command.vortexFile);
        validatePointVortexDomain(metadata, parameters, options.coordinateScale);
        auto vortices = readPointVortices(command.vortexFile, command.frame);
        const std::size_t count = vortices.size();
        auto field = imprintPointVortices(parameters, std::move(vortices), options);
        writeWavefunctionFile(command.outputFile, field, parameters, command.overwrite);
        std::cout << "wrote " << command.outputFile << " on a " << parameters.nx << 'x'
                  << parameters.ny << " grid from " << count << " vortices\n"
                  << "background density = " << options.backgroundDensity
                  << "\nhealing length = " << options.healingLength << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
