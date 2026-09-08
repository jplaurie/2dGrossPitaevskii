#pragma once

#include "backend.hpp"

#include <fftw3.h>

#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

template <class T> struct FftwAllocator {
  using value_type = T;
  FftwAllocator() = default;
  template <class U> FftwAllocator(const FftwAllocator<U> &) {}
  T *allocate(std::size_t count) {
    if (count > std::numeric_limits<std::size_t>::max() / sizeof(T))
      throw std::bad_array_new_length();
    auto *memory = static_cast<T *>(fftw_malloc(count * sizeof(T)));
    if (!memory)
      throw std::bad_alloc();
    return memory;
  }
  void deallocate(T *memory, std::size_t) { fftw_free(memory); }
  template <class U> bool operator==(const FftwAllocator<U> &) const {
    return true;
  }
};

using FftwComplexField = std::vector<Complex, FftwAllocator<Complex>>;

class FftwPlan {
public:
  FftwPlan() = default;
  ~FftwPlan() { reset(); }
  FftwPlan(const FftwPlan &) = delete;
  FftwPlan &operator=(const FftwPlan &) = delete;
  FftwPlan &operator=(fftw_plan plan) noexcept {
    if (plan != plan_)
      reset(plan);
    return *this;
  }
  explicit operator bool() const { return plan_ != nullptr; }
  operator fftw_plan() const { return plan_; }

private:
  void reset(fftw_plan plan = nullptr) noexcept {
    if (plan_)
      fftw_destroy_plan(plan_);
    plan_ = plan;
  }
  fftw_plan plan_ = nullptr;
};

void squarePointwise(const FftwComplexField &input, FftwComplexField &output,
                     std::size_t count);
void filterPaddedSpectrum(FftwComplexField &field, const Parameters &parameters,
                          std::size_t firstRow, std::size_t rowCount);
void multiplyConjugatePointwise(const FftwComplexField &left,
                                const FftwComplexField &right,
                                FftwComplexField &output, std::size_t count);

template <class Allocator>
inline fftw_complex *fftwData(std::vector<Complex, Allocator> &field) {
  static_assert(sizeof(Complex) == sizeof(fftw_complex));
  return reinterpret_cast<fftw_complex *>(field.data());
}

class BaseTransform {
public:
  explicit BaseTransform(const Parameters &parameters);
  ~BaseTransform();
  BaseTransform(const BaseTransform &) = delete;
  BaseTransform &operator=(const BaseTransform &) = delete;

  void forward(const std::vector<Complex> &physical, SpectralField &spectral);
  void inverse(const SpectralField &spectral, std::vector<Complex> &physical);
  void projectedSquareSpectra(const SpectralField &spectral,
                              const SpectralField &spectralRate,
                              SpectralField &square, SpectralField &squareRate);

private:
  class SquareTransform;
  Parameters p_;
  FftwPlan forward_, inverse_;
  FftwComplexField planningInput_, planningOutput_;
  std::unique_ptr<SquareTransform> squareTransform_;
};
