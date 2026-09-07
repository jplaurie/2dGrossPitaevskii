#pragma once

#include "backend.hpp"

#include <fftw3.h>

#include <limits>
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

private:
  Parameters p_;
  fftw_plan forward_ = nullptr;
  fftw_plan inverse_ = nullptr;
  FftwComplexField planningInput_, planningOutput_;
};
