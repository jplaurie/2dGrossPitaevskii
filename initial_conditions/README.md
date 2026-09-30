# Point-vortex initial conditions and imaginary-time relaxation

This directory contains the CPU workflow for turning a state produced by the
current [`2dPointVortex`](https://github.com/jplaurie/2dPointVortex) code into
an initial condition for this Gross--Pitaevskii (GP) solver:

1. `gp2d_vortex_imprint` reads a periodic PointVortex state, constructs a
   periodic GP phase, adds finite vortex cores, and writes a GP wavefunction.
2. `gp2d_relax` evolves that field in imaginary time until it is closer to a
   stationary GP solution. A uniform drift velocity can be included for a
   translating state such as a vortex dipole.
3. The relaxed file can be supplied to any real-time GP backend through
   `initialConditionFile`.

The relaxation executable is CPU-only by design. It uses the production
FFTW/OpenMP nonlinear backend and the same twice-dealiased cubic evaluation as
the real-time CPU solver. MPI or CUDA would add complexity without much benefit
to this one-off preparation stage; the resulting file is still accepted by the
CPU, MPI, and CUDA real-time executables.

## Requirements and build

The minimum requirements are CMake 3.20, a C++20 compiler, and FFTW3
development files. OpenMP and FFTW's threads library are optional. From the
`2dGrossPitaevskii` repository root, a CPU-only build is:

```bash
cmake -S . -B build/cpu -DCMAKE_BUILD_TYPE=Release \
  -DGP2D_MPI=OFF -DGP2D_CUDA=OFF
cmake --build build/cpu --target \
  gp2d_vortex_imprint gp2d_relax gross_pitaevskii_cpu --parallel
```

`make initial` is a convenience alternative, using `build/release` by default.
To generate a new PointVortex input file, also build the adjacent repository:

```bash
cmake -S ../2dPointVortex -B ../2dPointVortex/build/release \
  -DCMAKE_BUILD_TYPE=Release
cmake --build ../2dPointVortex/build/release --target \
  point_vortex_initial point_vortex_cpu --parallel
```

Run `ctest --test-dir build/cpu --output-on-failure` to test the complete
imprint, relaxation, and real-time hand-off.

## Quick start: a translating dipole

Run these commands from the `2dGrossPitaevskii` repository root. They generate
the same horizontal dipole stored in
[`examples/dipole.dat`](examples/dipole.dat): a positive vortex at the left, a
negative vortex at the right, and separation `d=pi/8` in a `2*pi` square.

```bash
mkdir -p runs/dipole

../2dPointVortex/build/release/point_vortex_initial \
  --geometry periodic \
  --case dipole \
  --box-length 6.283185307179586 \
  --ring-radius 0.392699081698724 \
  --output runs/dipole/vortices.dat
```

Imprint the two vortices on the GP grid described by
[`examples/relax.params`](examples/relax.params):

```bash
./build/cpu/gp2d_vortex_imprint \
  --parameters initial_conditions/examples/relax.params \
  --vortices runs/dipole/vortices.dat \
  --output runs/dipole/imprinted.dat
```

For the example coefficients, the point-vortex speed estimate is
`1/d=8/pi=2.546479089470326`. This dipole moves in the positive `y` direction,
so relax it in that moving frame:

```bash
./build/cpu/gp2d_relax \
  --parameters initial_conditions/examples/relax.params \
  --input runs/dipole/imprinted.dat \
  --output runs/dipole/relaxed.dat \
  --drift-y 2.546479089470326
```

Finally, start a conservative real-time run with the supplied example:

```bash
./build/cpu/gross_pitaevskii_cpu \
  initial_conditions/examples/realtime.params
```

All paths in these commands and parameter files are resolved relative to the
directory in which the executable is launched. Existing outputs are protected;
add `--overwrite` to the preparation command only when replacement is intended.

## Accepted `2dPointVortex` inputs

### Generated initial-condition files

The preferred input is the file written by `point_vortex_initial`. It contains
metadata followed by one `x y circulation` row per vortex:

```text
# PointVortex initial condition
# geometry=periodic pattern=dipole seed=1234567
# box_length=6.2831853071795862 disk_radius=1 infinite_half_width=1
# minimum_separation=0.39269908169872397
# x y circulation
-0.19634954084936199 0 1
 0.19634954084936199 0 -1
```

Only `geometry=periodic` is compatible with the GP solver's doubly periodic
domain. The importer checks the geometry and box length metadata and reports a
clear error if they disagree with the GP parameter file.

### PointVortex trajectories

A `trajectory.csv` written by the current PointVortex solver is also accepted:

```bash
./build/cpu/gp2d_vortex_imprint \
  --parameters my_relax.params \
  --vortices ../2dPointVortex/runs/my_run/trajectory.csv \
  --frame 25 \
  --output runs/my_state/imprinted.dat
```

The CSV header must be the current
`time,frame,index,x,y,circulation,u,v` format. `--frame N` selects frame `N`;
without it, the numerically largest frame is used. When
`resolved_parameters.txt` is beside `trajectory.csv`, its boundary condition
and box lengths are checked automatically. If the CSV was copied without that
run record, the user must verify the domain manually.

PointVortex coordinates are centered around zero, while the GP grid is stored
on `[0,Lx) x [0,Ly)`. The importer maps these periodically onto the GP grid
while retaining the centered representative when it chooses the background
phase-winding sector. It maps the sign of PointVortex circulation to GP winding
`+1` or `-1`; its magnitude is not used because this tool creates singly
quantized GP vortices.

## Domain, neutrality, and units

The GP domain is

```math
L_x=2\pi\,\texttt{aspectRatio},\qquad L_y=2\pi.
```

The current PointVortex periodic generator uses a square. Direct transfer
therefore normally uses `aspectRatio 1` and `--box-length 2*pi`, as in the quick
start. If the two codes use different but proportional length units, apply one
uniform conversion factor:

```bash
--coordinate-scale GP_LENGTH/POINT_VORTEX_LENGTH
```

For example, a PointVortex square of side `2` maps to the GP square of side
`2*pi` with `--coordinate-scale 3.141592653589793`. Metadata validation is
performed after this scale is applied. A single scale cannot map a square
PointVortex state to a rectangular GP domain.

A doubly periodic phase requires zero total winding, so the numbers of positive
and negative singly quantized vortices must balance. The tool rejects a
non-neutral state. It also accepts `--winding-x N` and `--winding-y N` to add an
integer whole-domain background phase winding; both default to zero.

The generated phase uses a Jacobi-theta representation plus the uniform phase
gradient required to make it exactly periodic in both directions. Each core
uses the order-eight
[Caliari--Zuccher Padé density profile](https://doi.org/10.1016/j.cpc.2017.09.013),
and the wavefunction modulus is the square root of that density.

## GP parameters used while imprinting

Both tools use the normal GP parameter-file syntax. For the defocusing,
finite-background convention

```math
i\partial_t\psi=c\Delta\psi+g|\psi|^2\psi+\mu\psi,
\qquad c<0,\quad g>0,\quad \mu<0,
```

the importer infers

```math
\rho_\infty=-\mu/g,
\qquad
\xi=\sqrt{c/\mu},
```

where `rho_infinity` is background density and `xi` is the core/healing-length
scale used by the Padé profile. Override these only when using another
nondimensionalization:

```text
--background-density RHO
--healing-length XI
```

The grid must resolve `xi` and the distance between nearby vortices. Here
`dx=Lx/nx` and `dy=Ly/ny`; repeat a calculation at higher `nx,ny` and verify
that the relaxed residual and observables no longer change materially. The
imprinted file contains exactly `nx*ny` real/imaginary pairs, so relaxation and
the subsequent real-time run must use the same grid dimensions and aspect
ratio.

The full imprint command-line interface is available with:

```bash
./build/cpu/gp2d_vortex_imprint --help
```

## Imaginary-time relaxation

The relaxer integrates

```math
\partial_\tau\psi
=-\left(c\Delta\psi+g|\psi|^2\psi+\mu\psi
        +i\boldsymbol{U}\cdot\nabla\psi\right).
```

A converged field therefore solves the GP equation in a frame translating at
`U=(drift-x,drift-y)`. Use zero drift for a stationary state. For a translating
dipole, start with its point-vortex velocity and adjust the relevant component
to minimize the final residual. Finite cores and periodic images mean that the
best GP drift can differ from the elementary `1/d` estimate.

The orientation and signs matter. In the quick-start state, `+` is left of `-`,
so the motion and `--drift-y` are positive. Reversing the circulation signs
reverses the drift; rotating that pair by 90 degrees moves the drift into the
`x` component.

The following normal parameter keys control relaxation:

| Parameter key | Meaning during relaxation |
| --- | --- |
| `nx`, `ny`, `aspectRatio` | Grid and periodic domain; must match the imprinted field |
| `dispersionCoefficient`, `nonlinearityCoefficient`, `chemicalPotential` | GP equation being solved |
| `timeStep` | Imaginary-time step |
| `numberOfSteps` | Maximum number of imaginary-time steps |
| `outputIntervalSteps` | Diagnostic and convergence-check interval |
| `integrator` | `etd2`, `etd3`, `etd4`, or `rk2` |
| `threadCount` | OpenMP/FFTW thread count; `0` uses the runtime default |

Forcing, hyperviscosity, hypoviscosity, and Ginzburg--Landau damping do not
participate in relaxation. Set them to zero/false anyway so the same file
clearly describes the intended subsequent conservative run. The algorithm is
grand-canonical: the chemical potential is fixed and wave action is not
renormalized after each step.

Useful relaxer options are:

| Option | Meaning |
| --- | --- |
| `--input FILE` | Imprinted field; defaults to `initialConditionFile` |
| `--output FILE` | Required relaxed wavefunction |
| `--drift-x U`, `--drift-y U` | Translating-frame velocity components |
| `--tolerance EPS` | Stop when the relative stationary residual is at most `EPS`; default `1e-10` |
| `--minimum-steps N` | Do not stop before step `N` |
| `--tolerance 0` | Disable early stopping and run exactly `numberOfSteps` |
| `--diagnostics FILE` | Override the default `OUTPUT.relaxation.csv` path |
| `--overwrite` | Replace both output and diagnostic files if present |

At step zero and every `outputIntervalSteps`, the console and diagnostics CSV
report:

- `comoving_energy`: the Hamiltonian including the translating-frame term;
- `wave_action`: the integral of `|psi|^2`;
- `relative_residual`: the norm of the stationary GP equation divided by the
  field norm.

Successful relaxation normally lowers the residual before it levels off at a
value set by timestep, resolution, drift accuracy, or the fact that the chosen
vortex arrangement is not a relative equilibrium. Energy is useful as a trend,
but the residual is the direct convergence criterion. Imaginary time does not
pin vortex coordinates: a nonstationary configuration may move, merge, or
annihilate.

## Choosing numerical settings

Start with `etd4` and the example timestep. If the field or diagnostics become
non-finite, or if the residual oscillates or grows, reduce `timeStep`. If the
residual decreases smoothly but has not plateaued, increase `numberOfSteps`.
Use a shorter `outputIntervalSteps` while tuning, then increase it for long
runs. Check the result at a finer spatial grid and a smaller imaginary-time
step before treating it as converged.

For a moving dipole, repeat relaxation over a small range of drift values and
compare the final residuals. A poor drift often causes the cores to translate
during relaxation and produces a residual plateau. More complicated
PointVortex snapshots generally are not exact uniformly translating GP states;
relaxation is then a way to remove core-profile radiation, not a guarantee that
all prescribed positions remain fixed.

## Using the result in a real-time run

Set

```text
initialConditionFile runs/dipole/relaxed.dat
```

in the real-time parameter file. Keep `nx`, `ny`, and `aspectRatio` identical.
Keep the GP coefficients identical if the relaxed field is meant to remain a
stationary or steadily translating solution. Choose new `dataDirectory` and
`outputDirectory` paths so the real-time run does not collide with an existing
run. The example [`examples/realtime.params`](examples/realtime.params) shows a
complete conservative configuration.

## Troubleshooting

- **“requires a doubly periodic PointVortex state”**: regenerate with
  `point_vortex_initial --geometry periodic`, or run the PointVortex simulation
  with `boundaryCondition periodic`.
- **“scaled PointVortex ... length does not match”**: use the same box length in
  both programs, or supply the uniform `--coordinate-scale` conversion.
- **“zero net vortex winding”**: use equal total counts of positive and negative
  vortices. A single vortex cannot exist by itself on this periodic torus.
- **“trajectory contains no rows for frame”**: inspect the second CSV column and
  choose an existing `--frame`, or omit the option to use the last frame.
- **“wavefunction must contain exactly 2*nx*ny numbers”**: the imprint and
  relaxation parameter files use different grids, or the input file is
  incomplete.
- **An output already exists**: select a new path, or explicitly pass
  `--overwrite`. This protection also applies to the relaxation CSV.
- **Residual stalls well above tolerance**: refine the grid, reduce the
  timestep, tune the drift, or accept that the configuration is not stationary
  in a uniformly moving frame.
- **Vortices disappear during relaxation**: check that the initial separation
  is resolved and physically viable; annihilation can also be the correct
  unconstrained imaginary-time evolution.

Use `gp2d_vortex_imprint --help` and `gp2d_relax --help` for the authoritative
option lists. The main [solver README](../README.md) documents every real-time
parameter and output file.
