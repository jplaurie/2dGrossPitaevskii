cmake_minimum_required(VERSION 3.20)
project(GrossPitaevskii2D VERSION 0.5.0 LANGUAGES C CXX)

find_package(Git QUIET)
set(GP2D_GIT_COMMIT "unknown")
set(GP2D_GIT_DIRTY "unknown")
if(Git_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GP2D_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GP2D_GIT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(GP2D_GIT_STATUS STREQUAL "")
    set(GP2D_GIT_DIRTY "false")
  else()
    set(GP2D_GIT_DIRTY "true")
  endif()
endif()

option(GP2D_OPENMP "Enable OpenMP and threaded FFTW" ON)
option(GP2D_MPI "Build the hybrid MPI/OpenMP backend" ON)
option(GP2D_CUDA "Build the CUDA/cuFFT backend" ON)
option(GP2D_HDF5 "Enable optional HDF5 field snapshots and exporter" ON)
option(GP2D_BACKEND_TESTS "Run MPI/CUDA integration tests (requires runtime access)" OFF)
set(GP2D_CUDA_ARCHITECTURES "" CACHE STRING "Optional CUDA architecture override")

find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3 REQUIRED)
find_library(FFTW_THREADS_LIBRARY NAMES fftw3_threads)
if(GP2D_HDF5)
  find_package(HDF5 COMPONENTS C QUIET)
endif()

set(GP2D_COMMON_SOURCES
    src/parameters.cpp
    src/fftw_utils.cpp
    src/spectral.cpp
    src/diagnostics.cpp
    src/vortex_diagnostics.cpp
    src/hdf5_io.cpp
    src/output.cpp
    src/output_transaction.cpp
    src/solver.cpp)
add_library(gp2d_common ${GP2D_COMMON_SOURCES})
target_include_directories(gp2d_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
target_link_libraries(gp2d_common PUBLIC ${FFTW_LIBRARY})
target_compile_features(gp2d_common PUBLIC cxx_std_20)
target_compile_definitions(gp2d_common PUBLIC
  GP2D_VERSION="${PROJECT_VERSION}" GP2D_GIT_COMMIT="${GP2D_GIT_COMMIT}"
  GP2D_GIT_DIRTY="${GP2D_GIT_DIRTY}")
if(GP2D_HDF5 AND HDF5_FOUND)
  target_include_directories(gp2d_common PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(gp2d_common PUBLIC ${HDF5_LIBRARIES})
  target_compile_definitions(gp2d_common PUBLIC GP2D_HAVE_HDF5=1)
elseif(GP2D_HDF5)
  message(WARNING "HDF5 not found; HDF5 field output and exporter will not be available")
endif()

if(GP2D_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(gp2d_common PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

set(GP2D_SERIAL_COMMON gp2d_common)
if(GP2D_OPENMP AND OpenMP_CXX_FOUND)
  add_library(gp2d_common_serial ${GP2D_COMMON_SOURCES})
  target_include_directories(gp2d_common_serial
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
  target_link_libraries(gp2d_common_serial PUBLIC ${FFTW_LIBRARY})
  target_compile_features(gp2d_common_serial PUBLIC cxx_std_20)
  target_compile_definitions(gp2d_common_serial PUBLIC
    GP2D_VERSION="${PROJECT_VERSION}" GP2D_GIT_COMMIT="${GP2D_GIT_COMMIT}"
    GP2D_GIT_DIRTY="${GP2D_GIT_DIRTY}")
  if(GP2D_HDF5 AND HDF5_FOUND)
    target_include_directories(gp2d_common_serial PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(gp2d_common_serial PUBLIC ${HDF5_LIBRARIES})
    target_compile_definitions(gp2d_common_serial PUBLIC GP2D_HAVE_HDF5=1)
  endif()
  set(GP2D_SERIAL_COMMON gp2d_common_serial)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(gp2d_common PRIVATE -Wall -Wextra -Wpedantic)
  if(TARGET gp2d_common_serial)
    target_compile_options(gp2d_common_serial PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

add_executable(gross_pitaevskii_cpu src/main.cpp src/backend_cpu.cpp)
target_link_libraries(gross_pitaevskii_cpu PRIVATE gp2d_common)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(gross_pitaevskii_cpu PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(gross_pitaevskii_cpu_serial src/main.cpp src/backend_cpu.cpp)
target_link_libraries(gross_pitaevskii_cpu_serial PRIVATE ${GP2D_SERIAL_COMMON})
add_executable(gp2d_benchmark_cpu_serial src/benchmark.cpp src/backend_cpu.cpp)
target_link_libraries(gp2d_benchmark_cpu_serial PRIVATE ${GP2D_SERIAL_COMMON})
add_executable(gp2d_benchmark_cpu src/benchmark.cpp src/backend_cpu.cpp)
target_link_libraries(gp2d_benchmark_cpu PRIVATE gp2d_common)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(gross_pitaevskii_cpu_serial PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(gp2d_benchmark_cpu_serial PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(gp2d_benchmark_cpu PRIVATE -Wall -Wextra -Wpedantic)
endif()
if(FFTW_THREADS_LIBRARY)
  target_link_libraries(gp2d_benchmark_cpu PRIVATE ${FFTW_THREADS_LIBRARY})
  target_compile_definitions(gp2d_benchmark_cpu PRIVATE GP2D_HAVE_FFTW_THREADS=1)
endif()
if(FFTW_THREADS_LIBRARY)
  target_link_libraries(gross_pitaevskii_cpu PRIVATE ${FFTW_THREADS_LIBRARY})
  target_compile_definitions(gross_pitaevskii_cpu PRIVATE GP2D_HAVE_FFTW_THREADS=1)
else()
  message(WARNING "FFTW threads library not found; CPU FFT plans remain single-threaded")
endif()

add_library(gp2d_initial_conditions initial_conditions/vortex_field.cpp)
target_include_directories(gp2d_initial_conditions PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/initial_conditions)
target_link_libraries(gp2d_initial_conditions PUBLIC gp2d_common)
target_compile_features(gp2d_initial_conditions PUBLIC cxx_std_20)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(gp2d_initial_conditions PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(gp2d_vortex_imprint initial_conditions/vortex_imprint.cpp)
target_link_libraries(gp2d_vortex_imprint PRIVATE gp2d_initial_conditions)

if(GP2D_HDF5 AND HDF5_FOUND)
  add_executable(gp2d_hdf5_export src/hdf5_export.cpp)
  target_link_libraries(gp2d_hdf5_export PRIVATE gp2d_common)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gp2d_hdf5_export PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

# Imaginary-time relaxation deliberately uses the production CPU nonlinear
# backend, so its dealiasing and normalization stay identical to real-time runs.
add_executable(gp2d_relax initial_conditions/imaginary_time.cpp src/backend_cpu.cpp)
target_link_libraries(gp2d_relax PRIVATE gp2d_initial_conditions)
if(FFTW_THREADS_LIBRARY)
  target_link_libraries(gp2d_relax PRIVATE ${FFTW_THREADS_LIBRARY})
  target_compile_definitions(gp2d_relax PRIVATE GP2D_HAVE_FFTW_THREADS=1)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(gp2d_vortex_imprint PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(gp2d_relax PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(GP2D_MPI)
  find_package(MPI COMPONENTS CXX)
  find_path(FFTW_MPI_INCLUDE_DIR fftw3-mpi.h)
  find_library(FFTW_MPI_LIBRARY NAMES fftw3_mpi)
  if(MPI_CXX_FOUND AND FFTW_MPI_INCLUDE_DIR AND FFTW_MPI_LIBRARY)
    add_executable(gross_pitaevskii_mpi src/main.cpp src/backend_mpi.cpp)
    target_include_directories(gross_pitaevskii_mpi PRIVATE ${FFTW_MPI_INCLUDE_DIR})
    target_link_libraries(gross_pitaevskii_mpi PRIVATE gp2d_common MPI::MPI_CXX
                                                        ${FFTW_MPI_LIBRARY})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(gross_pitaevskii_mpi PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    if(FFTW_THREADS_LIBRARY)
      target_link_libraries(gross_pitaevskii_mpi PRIVATE ${FFTW_THREADS_LIBRARY})
      target_compile_definitions(gross_pitaevskii_mpi PRIVATE GP2D_HAVE_FFTW_THREADS=1)
    endif()
    add_executable(gp2d_benchmark_mpi src/benchmark.cpp src/backend_mpi.cpp)
    target_include_directories(gp2d_benchmark_mpi PRIVATE ${FFTW_MPI_INCLUDE_DIR})
    target_link_libraries(gp2d_benchmark_mpi PRIVATE gp2d_common MPI::MPI_CXX
                                                      ${FFTW_MPI_LIBRARY})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(gp2d_benchmark_mpi PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    if(FFTW_THREADS_LIBRARY)
      target_link_libraries(gp2d_benchmark_mpi PRIVATE ${FFTW_THREADS_LIBRARY})
      target_compile_definitions(gp2d_benchmark_mpi PRIVATE GP2D_HAVE_FFTW_THREADS=1)
    endif()
  else()
    message(WARNING "MPI and/or FFTW-MPI not found; gross_pitaevskii_mpi will not be built")
  endif()
endif()

if(GP2D_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_executable(gross_pitaevskii_cuda src/main.cpp src/backend_cuda.cu)
    add_executable(gross_pitaevskii_cuda_mixed src/main.cpp src/backend_cuda.cu)
    add_executable(gp2d_benchmark_cuda src/benchmark.cpp src/backend_cuda.cu)
    add_executable(gp2d_benchmark_cuda_mixed src/benchmark.cpp src/backend_cuda.cu)
    target_link_libraries(gross_pitaevskii_cuda PRIVATE gp2d_common CUDA::cudart CUDA::cufft)
    target_link_libraries(gross_pitaevskii_cuda_mixed PRIVATE gp2d_common CUDA::cudart CUDA::cufft)
    target_link_libraries(gp2d_benchmark_cuda PRIVATE gp2d_common CUDA::cudart CUDA::cufft)
    target_link_libraries(gp2d_benchmark_cuda_mixed PRIVATE gp2d_common CUDA::cudart CUDA::cufft)
    target_compile_definitions(gross_pitaevskii_cuda_mixed PRIVATE GP2D_CUDA_MIXED=1)
    target_compile_definitions(gp2d_benchmark_cuda_mixed PRIVATE GP2D_CUDA_MIXED=1)
    set(GP2D_CUDA_TARGETS gross_pitaevskii_cuda gross_pitaevskii_cuda_mixed
      gp2d_benchmark_cuda gp2d_benchmark_cuda_mixed)
    set_target_properties(${GP2D_CUDA_TARGETS}
      PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
    if(NOT GP2D_CUDA_ARCHITECTURES STREQUAL "")
      set_target_properties(${GP2D_CUDA_TARGETS}
        PROPERTIES CUDA_ARCHITECTURES "${GP2D_CUDA_ARCHITECTURES}")
    endif()
  else()
    message(WARNING "CUDA toolkit not found; gross_pitaevskii_cuda will not be built")
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(gp2d_parameters_tests tests/parameters.cpp)
  target_link_libraries(gp2d_parameters_tests PRIVATE gp2d_common)
  add_test(NAME parameter_parsing COMMAND gp2d_parameters_tests)

  add_executable(gp2d_numerics_tests tests/numerics.cpp src/backend_cpu.cpp)
  target_link_libraries(gp2d_numerics_tests PRIVATE gp2d_common)
  if(FFTW_THREADS_LIBRARY)
    target_link_libraries(gp2d_numerics_tests PRIVATE ${FFTW_THREADS_LIBRARY})
    target_compile_definitions(gp2d_numerics_tests PRIVATE GP2D_HAVE_FFTW_THREADS=1)
  endif()
  add_test(NAME nonlinear_direct_dft COMMAND gp2d_numerics_tests)
  add_executable(gp2d_initial_condition_tests tests/initial_conditions.cpp)
  target_link_libraries(gp2d_initial_condition_tests PRIVATE gp2d_initial_conditions)
  add_test(NAME initial_condition_math COMMAND gp2d_initial_condition_tests)
  set(GP2D_SMOKE_DIR "${CMAKE_CURRENT_BINARY_DIR}/smoke")
  file(MAKE_DIRECTORY "${GP2D_SMOKE_DIR}")
  configure_file(tests/smoke.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/smoke.cmake" @ONLY)
  add_test(NAME cpu_output_smoke COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/smoke.cmake")
  set(GP2D_INITIAL_SMOKE_DIR "${CMAKE_CURRENT_BINARY_DIR}/initial_condition_smoke")
  configure_file(tests/initial_conditions_smoke.cmake.in
                 "${CMAKE_CURRENT_BINARY_DIR}/initial_conditions_smoke.cmake" @ONLY)
  add_test(NAME initial_condition_workflow
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/initial_conditions_smoke.cmake")
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME regression_cpu COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.py cpu
      --cpu $<TARGET_FILE:gross_pitaevskii_cpu>)
    add_test(NAME temporal_convergence COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/convergence.py
      --cpu $<TARGET_FILE:gross_pitaevskii_cpu>)
    if(TARGET gp2d_hdf5_export)
      add_test(NAME hdf5_output COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/hdf5_output.py
        --cpu $<TARGET_FILE:gross_pitaevskii_cpu>
        --exporter $<TARGET_FILE:gp2d_hdf5_export>)
    endif()
    if(GP2D_BACKEND_TESTS)
      foreach(backend IN ITEMS mpi cuda cuda_mixed)
        if(TARGET gross_pitaevskii_${backend})
          set(launcher)
          if(backend STREQUAL "mpi")
            list(APPEND launcher --mpiexec ${MPIEXEC_EXECUTABLE}
              --mpi-numproc-flag=${MPIEXEC_NUMPROC_FLAG})
          endif()
          if(backend STREQUAL "cuda_mixed")
            list(APPEND launcher --tolerance 5e-8)
          endif()
          add_test(NAME regression_${backend} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.py backend
            --cpu $<TARGET_FILE:gross_pitaevskii_cpu>
            --candidate $<TARGET_FILE:gross_pitaevskii_${backend}> ${launcher})
          set_tests_properties(regression_${backend} PROPERTIES LABELS "backend" TIMEOUT 300)
        endif()
      endforeach()
    endif()
  endif()
endif()
