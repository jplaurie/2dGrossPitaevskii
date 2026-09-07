#include "fftw_utils.hpp"

BaseTransform::BaseTransform(const Parameters &p)
    : p_(p), planningInput_(p.nx * p.ny), planningOutput_(p.nx * p.ny) {
  forward_ = fftw_plan_dft_2d(
      static_cast<int>(p.ny), static_cast<int>(p.nx), fftwData(planningInput_),
      fftwData(planningOutput_), FFTW_FORWARD, FFTW_ESTIMATE);
  inverse_ = fftw_plan_dft_2d(
      static_cast<int>(p.ny), static_cast<int>(p.nx), fftwData(planningInput_),
      fftwData(planningOutput_), FFTW_BACKWARD, FFTW_ESTIMATE);
  if (!forward_ || !inverse_) {
    if (forward_)
      fftw_destroy_plan(forward_);
    if (inverse_)
      fftw_destroy_plan(inverse_);
    throw std::runtime_error("FFTW could not create base-grid plans");
  }
}

BaseTransform::~BaseTransform() {
  if (forward_)
    fftw_destroy_plan(forward_);
  if (inverse_)
    fftw_destroy_plan(inverse_);
}

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
