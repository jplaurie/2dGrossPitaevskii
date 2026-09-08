# 2D Gross–Pitaevskii pseudo-spectral solver

A C++20 solver for a complex Gross–Pitaevskii field on a doubly periodic
domain. A shared numerical model is available through OpenMP, hybrid
MPI/OpenMP, and CUDA backends.

Three executables share one model, parameter format, set of time integrators,
and output format:

| Executable | Backend | Use case |
| --- | --- | --- |
| `gross_pitaevskii_cpu` | FFTW, with optional OpenMP | Shared-memory runs |
| `gross_pitaevskii_mpi` | FFTW-MPI, with optional OpenMP | Distributed-memory runs |
| `gross_pitaevskii_cuda` | CUDA and cuFFT | NVIDIA GPU runs |

The solver provides dealiased pseudo-spectral nonlinear evaluation, four
fixed-step exponential or integrating-factor schemes, reproducible stochastic
forcing, atomic checkpoints, automatic restart, and cross-backend regression
tests.

## Equation and discretization

The solver evolves normalized Fourier coefficients of the complex field
`psi` on `Lx = 2*pi*aspectRatio`, `Ly = 2*pi`:

```math
\partial_t\psi_k =
\frac{(-c|k|^2+\mu)\psi_k+g\,\widehat{|\psi|^2\psi}_k}
     {i-\Gamma_k}
-\left[\nu |k|^{2p}+\alpha |k|^{2q}\right]\psi_k+F_k.
```

Here `c`, `g`, and `mu` are `dispersionCoefficient`,
`nonlinearityCoefficient`, and `chemicalPotential`. The optional
Ginzburg–Landau coefficient `Gamma_k` is applied above its configured cutoff.
Hyper- and hypoviscosity can each optionally be limited to one side of a
wavenumber cutoff.

The cubic term uses a two-pass 3/2-rule treatment:
`psi` is embedded on the padded grid, `psi^2` is transformed and truncated to
the retained band, and that filtered product is multiplied by `conj(psi)`
before the final transform. This is important for a cubic nonlinearity and
provides a momentum-conserving discretization.

Fixed-step integrators are ETDRK2 (`etd2`), ETDRK3 (`etd3`), ETDRK4-B
(`etd4`), and second-order integrating-factor Runge–Kutta (`rk2`). The full
linear operator is integrated analytically. Stochastic forcing uses the exact
linear covariance over one step, and its generator state is checkpointed.

## Requirements

The CPU build requires CMake 3.20+, a C++20 compiler, and FFTW3 development
headers and libraries. OpenMP and FFTW's threads library are optional. The
hybrid executable also needs MPI and FFTW-MPI. The CUDA executable needs the
NVIDIA CUDA Toolkit and cuFFT, plus an NVIDIA GPU at run time. Python 3 enables
the end-to-end regression tests.

For example, the required packages can be installed on Arch Linux with:

```bash
sudo pacman -S cmake gcc fftw openmpi fftw-openmpi cuda python
```

## Build and test

Start with a portable CPU build:

```bash
cmake -S . -B build/cpu -DCMAKE_BUILD_TYPE=Release \
  -DGP2D_MPI=OFF -DGP2D_CUDA=OFF
cmake --build build/cpu -j
ctest --test-dir build/cpu --output-on-failure
```

To build every backend supported by the local toolchain and enable the MPI and
CUDA regression tests:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
  -DGP2D_BACKEND_TESTS=ON
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
```

CMake omits MPI if MPI or FFTW-MPI is unavailable and omits CUDA if no CUDA
compiler is found. Useful options are:

```text
-DGP2D_OPENMP=OFF
-DGP2D_MPI=OFF
-DGP2D_CUDA=OFF
-DGP2D_CUDA_ARCHITECTURES=<CUDA architecture>
-DGP2D_BACKEND_TESTS=ON
```

Backend tests are opt-in because they require a working MPI launcher and, for
CUDA, a visible GPU. Convenience targets are `make cpu`, `make mpi`,
`make cuda`, and `make test`; set `BUILD_DIR` if desired.

## Quick start

[`examples/quickstart.params`](examples/quickstart.params) is a small,
repeatable forced run:

```bash
./build/cpu/gross_pitaevskii_cpu examples/quickstart.params
```

It writes snapshots and checkpoints under `data/quickstart/` and CSV
diagnostics under `output/quickstart/`. Remove those directories or change the
paths before starting a fresh run.

All executables accept an optional parameter-file path and otherwise read
`params.txt`:

```bash
./build/release/gross_pitaevskii_cpu run.params

mpirun -n 2 ./build/release/gross_pitaevskii_mpi run.params

./build/release/gross_pitaevskii_cuda run.params
```

`threadCount` selects OpenMP and threaded-FFTW threads per process; zero uses
the OpenMP runtime default, including `OMP_NUM_THREADS` when it is set. For
MPI, plan for `ranks * threadCount` CPU cores. Only rank zero writes files.
CUDA keeps time-integration stages and nonlinear FFTs on the device; host
transfers occur for stochastic increments and output.

## Parameter files

Each nonempty line is `key value` or `key = value`; `#` begins a comment. Keys
are case-sensitive, and invalid keys, values, or extra fields stop the run.

| Key | Purpose |
| --- | --- |
| `nx`, `ny` | Even physical-grid dimensions, each at least four |
| `aspectRatio` | Positive `Lx/(2*pi)` |
| `timeStep`, `numberOfSteps` | Fixed step and additional steps this invocation |
| `outputIntervalSteps` | Save cadence; the final step is always saved |
| `integrator` | `etd2`, `etd3`, `etd4`, or `rk2` |
| `dispersionCoefficient` | `c` in the equation |
| `nonlinearityCoefficient` | Cubic coefficient `g` |
| `chemicalPotential` | Linear coefficient `mu` |
| `hyperviscosity`, `hyperviscosityOrder` | `nu` and positive-scale spectral power `p` |
| `hyperviscosityCutoffEnabled`, `hyperviscosityCutoff` | Apply hyperviscosity only above the cutoff |
| `hypoviscosity`, `hypoviscosityOrder` | `alpha` and spectral power `q` (negative is allowed) |
| `hypoviscosityCutoffEnabled`, `hypoviscosityCutoff` | Apply hypoviscosity only below the cutoff |
| `ginzburgLandauDamping`, `ginzburgLandauCutoff` | `Gamma_k` strength and lower wavenumber threshold |
| `forcingEnabled` | Enable or disable forcing |
| `forcingProfile` | `annulus`, `gaussian`, `exponential`, `logNormal`, or `singleMode` |
| `forcingWavenumber` | Profile center; an integer Fourier-mode index for `singleMode` |
| `forcingWidth` | Annulus half-width or Gaussian standard deviation |
| `forcingAmplitude` | Spectral forcing amplitude |
| `forcingShapeOrder` | Exponent for `exponential` forcing |
| `forcingLogWidth` | Log-space sigma for `logNormal` forcing |
| `targetWaveActionInjectionRate` | Normalize stochastic forcing when positive |
| `randomSeed` | Reproducible 64-bit seed; zero chooses and records a time-based seed |
| `writeModeDiagnostics` | Write selected complex Fourier modes |
| `threadCount` | Host threads per process; zero uses the runtime default |
| `overwriteOutput` | Permit frame replacement; it does not disable restart detection |
| `initialConditionFile` | Optional physical complex field |
| `dataDirectory`, `outputDirectory` | State and diagnostic locations |

Booleans accept `true`/`false` or `1`/`0`. `singleMode` is deterministic at
the four modes with `|kx index| = |ky index| = forcingWavenumber`; the other
profiles use circular complex Gaussian, white-in-time forcing.
Forcing-specific numeric constraints are checked only when `forcingEnabled` is
true; disabled forcing parameters are parsed but otherwise ignored.
For negative `hypoviscosityOrder`, the hypoviscous multiplier is singular at
`k=0`. The mean mode is therefore explicitly set to zero on initialization and
after every time step.

An initial-condition file contains `2*nx*ny` whitespace-delimited numbers,
ordered as row-major `real imag` pairs. A `wavefunction_NNNNNNNN.dat` snapshot
has this format and can be used directly as an initial condition. Relative
paths are resolved from the directory in which the executable is launched.

## Output and restart

Fresh runs save frame zero, then the requested cadence and final step.

| Location | Contents |
| --- | --- |
| `dataDirectory/wavefunction_NNNNNNNN.dat` | Physical `real imag` pairs (`ny` by `2*nx`) |
| `dataDirectory/checkpoint_NNNNNNNN.bin` | Normalized complex spectral state |
| `dataDirectory/restart_state.txt` | Latest time, frame, grid identity, and RNG state |
| `outputDirectory/diagnostics.csv` | Hamiltonian components, wave action, and damping rates |
| `outputDirectory/spectra.csv` | Wave-action and quadratic-energy shell spectra |
| `outputDirectory/fluxes.csv` | Wave-action and full Hamiltonian-energy fluxes |
| `outputDirectory/modes.csv` | Optional selected complex Fourier modes |
| `outputDirectory/forcing_*.csv` | Forcing summary and spectrum |
| `outputDirectory/segments/` | Per-invocation resolved parameters and forcing records |

Here `quadratic_energy` means the kinetic-plus-chemical-potential part,
`area * sum_k (-c*|k|^2 + mu)*|psi_k|^2`. The diagnostics keep its kinetic
and potential components in separate columns; total energy additionally
includes the quartic nonlinear contribution. For consistency with the
two-pass de-aliased cubic term, quartic energy is evaluated from the retained
spectrum of `psi^2`. The energy flux includes both the quadratic and quartic
transfers and uses the legacy low-wavenumber cumulative convention. The
`total_energy_dissipation_hypo` and
`total_energy_dissipation_hyper` columns include each operator's effect on
both the quadratic and quartic Hamiltonian energy; the similarly named
`quadratic_energy_*` columns retain the quadratic-only values.

When `restart_state.txt` exists, the solver resumes automatically.
`numberOfSteps` means additional steps, CSV files append, and frame numbering
continues. Matching `nx`, `ny`, and `aspectRatio` are required; physical
parameters may change between invocations. Output frames are journaled and
committed atomically, so a partially written frame is rolled back on restart.
Use new data and output directories for an independent run. Existing restart
metadata always resumes that run; `overwriteOutput` only permits replacement
of colliding output files.

## Code structure

```text
src/
  main.cpp                    shared executable entry point
  parameters.cpp/.hpp         parse, validate, and record settings
  spectral.cpp/.hpp           Fourier indexing and retained-band rules
  fftw_utils.cpp/.hpp         base-grid complex FFTW transforms
  solver.cpp/.hpp             linear operator, forcing, and time stepping
  integrator.hpp              shared CPU/CUDA stage formulas
  output.cpp/.hpp             diagnostics, snapshots, and checkpoints
  output_transaction.cpp      atomic output recovery and run history
  backend.hpp                 common nonlinear-backend interface
  backend_cpu.cpp             FFTW/OpenMP cubic backend
  backend_mpi.cpp             FFTW-MPI/OpenMP cubic backend
  backend_cuda.cu             CUDA/cuFFT backend and device time stepping
tests/
  numerics.cpp                direct-DFT nonlinear verification
  regression.py               restart and cross-backend comparisons
```
