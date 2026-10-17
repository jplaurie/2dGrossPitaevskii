#include "fftw_utils.hpp"
#include "parallel.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <filesystem>
#include <utility>

namespace {
unsigned planningFlags = FFTW_ESTIMATE;
std::filesystem::path wisdomPath;
} // namespace

void configureFftw(const Parameters &parameters, bool importWisdom) {
    switch (parameters.fftwPlanning) {
    case FftwPlanning::estimate:
        planningFlags = FFTW_ESTIMATE;
        break;
    case FftwPlanning::measure:
        planningFlags = FFTW_MEASURE;
        break;
    case FftwPlanning::patient:
        planningFlags = FFTW_PATIENT;
        break;
    }
    wisdomPath = parameters.fftwWisdomFile;
    if (importWisdom && !wisdomPath.empty() && std::filesystem::exists(wisdomPath) &&
        !fftw_import_wisdom_from_filename(wisdomPath.c_str()))
        throw std::runtime_error("cannot import FFTW wisdom: " + wisdomPath.string());
}

void saveFftwWisdom() {
    // Consume the request even if export fails: error cleanup must not retry it.
    const auto path = std::exchange(wisdomPath, {});
    if (!path.empty() && !fftw_export_wisdom_to_filename(path.c_str()))
        throw std::runtime_error("cannot export FFTW wisdom: " + path.string());
}

unsigned fftwPlanningFlags() { return planningFlags; }

class BaseTransform::DensityTransform {
  public:
    explicit DensityTransform(const Parameters &parameters)
        : parameters_(parameters), field_(parameters.mx() * parameters.my()),
          rate_(parameters.mx() * parameters.my()), work_(parameters.mx() * parameters.my()) {
        inverse_ = makePlan(field_, FFTW_BACKWARD);
        forward_ = makePlan(work_, FFTW_FORWARD);
        if (!inverse_ || !forward_)
            throw std::runtime_error("FFTW could not create density-transform plans");
    }

    void evaluate(const SpectralField &spectral, const SpectralField &spectralRate,
                  SpectralField &density, SpectralField &densityRate) {
        load(spectral, field_);
        transformProduct(field_, field_, 1.0, density);
        transformRate(spectralRate, densityRate);
    }

  private:
    fftw_plan makePlan(FftwComplexField &field, int direction) {
        return fftw_plan_dft_2d(static_cast<int>(parameters_.my()),
                                static_cast<int>(parameters_.mx()), fftwData(field),
                                fftwData(field), direction, fftwPlanningFlags());
    }

    void load(const SpectralField &input, FftwComplexField &padded) {
        if (input.size() != parameters_.nx * parameters_.ny)
            throw std::runtime_error("invalid density-transform input size");
        std::fill(padded.begin(), padded.end(), Complex{});
        forEachIndex(input.size(), [&](std::size_t i) {
            padded[paddedIndexForBaseMode(parameters_, i)] = input[i];
        });
        fftw_execute_dft(inverse_, fftwData(padded), fftwData(padded));
    }

    void transformRate(const SpectralField &input, SpectralField &output) {
        load(input, rate_);
        transformProduct(field_, rate_, 2.0, output);
    }

    void transformProduct(const FftwComplexField &left, const FftwComplexField &right,
                          double factor, SpectralField &output) {
        forEachIndex(work_.size(), [&](std::size_t i) {
            work_[i] = factor * std::real(std::conj(left[i]) * right[i]);
        });
        fftw_execute(forward_);
        extract(work_, output);
    }

    void extract(const FftwComplexField &padded, SpectralField &output) {
        output.resize(parameters_.nx * parameters_.ny);
        const double scale = 1.0 / static_cast<double>(parameters_.mx() * parameters_.my());
        forEachIndex(output.size(), [&](std::size_t i) {
            const auto x = i % parameters_.nx, y = i / parameters_.nx;
            output[i] = retainedDensityWave(signedWave(x, parameters_.nx), parameters_.nx) &&
                                retainedDensityWave(signedWave(y, parameters_.ny), parameters_.ny)
                            ? scale * padded[paddedIndexForBaseMode(parameters_, i)]
                            : Complex{};
        });
    }

    Parameters parameters_;
    FftwComplexField field_, rate_, work_;
    FftwPlan inverse_, forward_;
};

void densityPointwise(const FftwComplexField &input, FftwComplexField &output, std::size_t count) {
    forEachIndex(count, [&](std::size_t i) { output[i] = std::norm(input[i]); });
}

void filterDensitySpectrum(FftwComplexField &field, const Parameters &parameters,
                           std::size_t firstRow, std::size_t rowCount) {
    const std::size_t count = rowCount * parameters.mx();
    const double scale = 1.0 / static_cast<double>(parameters.mx() * parameters.my());
    forEachIndex(count, [&](std::size_t i) {
        const std::size_t x = i % parameters.mx(), y = firstRow + i / parameters.mx();
        if (retainedDensityWave(signedWave(x, parameters.mx()), parameters.nx) &&
            retainedDensityWave(signedWave(y, parameters.my()), parameters.ny))
            field[i] *= scale;
        else
            field[i] = Complex{};
    });
}

void multiplyPointwise(const FftwComplexField &left, const FftwComplexField &right,
                       FftwComplexField &output, std::size_t count) {
    forEachIndex(count, [&](std::size_t i) { output[i] = left[i] * right[i]; });
}

BaseTransform::BaseTransform(const Parameters &parameters)
    : parameters_(parameters), planningInput_(parameters.nx * parameters.ny),
      planningOutput_(parameters.nx * parameters.ny) {
    forward_ = fftw_plan_dft_2d(static_cast<int>(parameters.ny), static_cast<int>(parameters.nx),
                                fftwData(planningInput_), fftwData(planningOutput_), FFTW_FORWARD,
                                fftwPlanningFlags());
    inverse_ = fftw_plan_dft_2d(static_cast<int>(parameters.ny), static_cast<int>(parameters.nx),
                                fftwData(planningInput_), fftwData(planningOutput_), FFTW_BACKWARD,
                                fftwPlanningFlags());
    if (!forward_ || !inverse_)
        throw std::runtime_error("FFTW could not create base-grid plans");
}

BaseTransform::~BaseTransform() = default;

void BaseTransform::forward(const std::vector<Complex> &physical, SpectralField &spectral) {
    if (physical.size() != parameters_.nx * parameters_.ny)
        throw std::runtime_error("invalid physical field size");
    planningInput_.assign(physical.begin(), physical.end());
    fftw_execute(forward_);
    spectral.assign(planningOutput_.begin(), planningOutput_.end());
    const double scale = 1.0 / static_cast<double>(parameters_.nx * parameters_.ny);
    for (Complex &value : spectral)
        value *= scale;
}

void BaseTransform::inverse(const SpectralField &spectral, std::vector<Complex> &physical) {
    if (spectral.size() != parameters_.nx * parameters_.ny)
        throw std::runtime_error("invalid spectral field size");
    planningInput_.assign(spectral.begin(), spectral.end());
    fftw_execute(inverse_);
    physical.assign(planningOutput_.begin(), planningOutput_.end());
}

void BaseTransform::projectedDensitySpectra(const SpectralField &spectral,
                                            const SpectralField &spectralRate,
                                            SpectralField &density, SpectralField &densityRate) {
    if (!densityTransform_)
        densityTransform_ = std::make_unique<DensityTransform>(parameters_);
    densityTransform_->evaluate(spectral, spectralRate, density, densityRate);
}
