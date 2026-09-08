#include "fftw_utils.hpp"
#include "parallel.hpp"
#include "spectral.hpp"

#include <algorithm>

class BaseTransform::SquareTransform {
public:
  explicit SquareTransform(const Parameters &p)
      : p_(p), field_(p.mx() * p.my()), rate_(p.mx() * p.my()),
        work_(p.mx() * p.my()) {
    inverse_ = makePlan(field_, FFTW_BACKWARD);
    forward_ = makePlan(work_, FFTW_FORWARD);
    if (!inverse_ || !forward_)
      throw std::runtime_error("FFTW could not create square-transform plans");
  }

  void evaluate(const SpectralField &spectral,
                const SpectralField &spectralRate, SpectralField &square,
                SpectralField &squareRate) {
    load(spectral, field_);
    transformProduct(field_, field_, 1.0, square);
    transformRate(spectralRate, squareRate);
  }

private:
  fftw_plan makePlan(FftwComplexField &field, int direction) {
    return fftw_plan_dft_2d(static_cast<int>(p_.my()),
                            static_cast<int>(p_.mx()), fftwData(field),
                            fftwData(field), direction, FFTW_ESTIMATE);
  }

  void load(const SpectralField &input, FftwComplexField &padded) {
    if (input.size() != p_.nx * p_.ny)
      throw std::runtime_error("invalid square-transform input size");
    std::fill(padded.begin(), padded.end(), Complex{});
    forEachIndex(input.size(), [&](std::size_t i) {
      padded[paddedIndexForBaseMode(p_, i)] = input[i];
    });
    fftw_execute_dft(inverse_, fftwData(padded), fftwData(padded));
  }

  void transformRate(const SpectralField &input, SpectralField &output) {
    load(input, rate_);
    transformProduct(field_, rate_, 2.0, output);
  }

  void transformProduct(const FftwComplexField &left,
                        const FftwComplexField &right, double factor,
                        SpectralField &output) {
    forEachIndex(work_.size(), [&](std::size_t i) {
      work_[i] = factor * left[i] * right[i];
    });
    fftw_execute(forward_);
    extract(work_, output);
  }

  void extract(const FftwComplexField &padded, SpectralField &output) {
    output.resize(p_.nx * p_.ny);
    const double scale = 1.0 / static_cast<double>(p_.mx() * p_.my());
    forEachIndex(output.size(), [&](std::size_t i) {
      output[i] = scale * padded[paddedIndexForBaseMode(p_, i)];
    });
  }

  Parameters p_;
  FftwComplexField field_, rate_, work_;
  FftwPlan inverse_, forward_;
};

void squarePointwise(const FftwComplexField &input, FftwComplexField &output,
                     std::size_t count) {
  forEachIndex(count, [&](std::size_t i) { output[i] = input[i] * input[i]; });
}

void filterPaddedSpectrum(FftwComplexField &field, const Parameters &p,
                          std::size_t firstRow, std::size_t rowCount) {
  const std::size_t count = rowCount * p.mx();
  const double scale = 1.0 / static_cast<double>(p.mx() * p.my());
  forEachIndex(count, [&](std::size_t i) {
    const std::size_t x = i % p.mx(), y = firstRow + i / p.mx();
    if (retainedPaddedWave(signedWave(x, p.mx()), p.nx) &&
        retainedPaddedWave(signedWave(y, p.my()), p.ny))
      field[i] *= scale;
    else
      field[i] = Complex{};
  });
}

void multiplyConjugatePointwise(const FftwComplexField &left,
                                const FftwComplexField &right,
                                FftwComplexField &output, std::size_t count) {
  forEachIndex(
      count, [&](std::size_t i) { output[i] = left[i] * std::conj(right[i]); });
}

BaseTransform::BaseTransform(const Parameters &p)
    : p_(p), planningInput_(p.nx * p.ny), planningOutput_(p.nx * p.ny) {
  forward_ = fftw_plan_dft_2d(
      static_cast<int>(p.ny), static_cast<int>(p.nx), fftwData(planningInput_),
      fftwData(planningOutput_), FFTW_FORWARD, FFTW_ESTIMATE);
  inverse_ = fftw_plan_dft_2d(
      static_cast<int>(p.ny), static_cast<int>(p.nx), fftwData(planningInput_),
      fftwData(planningOutput_), FFTW_BACKWARD, FFTW_ESTIMATE);
  if (!forward_ || !inverse_)
    throw std::runtime_error("FFTW could not create base-grid plans");
}

BaseTransform::~BaseTransform() = default;

void BaseTransform::forward(const std::vector<Complex> &physical,
                            SpectralField &spectral) {
  if (physical.size() != p_.nx * p_.ny)
    throw std::runtime_error("invalid physical field size");
  planningInput_.assign(physical.begin(), physical.end());
  fftw_execute(forward_);
  spectral.assign(planningOutput_.begin(), planningOutput_.end());
  const double scale = 1.0 / static_cast<double>(p_.nx * p_.ny);
  for (Complex &value : spectral)
    value *= scale;
}

void BaseTransform::inverse(const SpectralField &spectral,
                            std::vector<Complex> &physical) {
  if (spectral.size() != p_.nx * p_.ny)
    throw std::runtime_error("invalid spectral field size");
  planningInput_.assign(spectral.begin(), spectral.end());
  fftw_execute(inverse_);
  physical.assign(planningOutput_.begin(), planningOutput_.end());
}

void BaseTransform::projectedSquareSpectra(const SpectralField &spectral,
                                           const SpectralField &spectralRate,
                                           SpectralField &square,
                                           SpectralField &squareRate) {
  if (!squareTransform_)
    squareTransform_ = std::make_unique<SquareTransform>(p_);
  squareTransform_->evaluate(spectral, spectralRate, square, squareRate);
}
