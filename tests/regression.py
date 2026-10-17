#!/usr/bin/env python3
import argparse
import cmath
import csv
import math
import struct
import subprocess
import tempfile
from pathlib import Path


def write_params(folder, **overrides):
    folder.mkdir(parents=True, exist_ok=True)
    values = dict(nx=12, ny=14, aspectRatio=1.5, timeStep=2e-4,
                  numberOfSteps=4, outputIntervalSteps=2, integrator="etd4",
                  dispersionCoefficient=-1, nonlinearityCoefficient=2,
                  chemicalPotential=.1, hyperviscosity=1e-6,
                  hyperviscosityOrder=2, hypoviscosity=0,
                  forcingEnabled="true", forcingProfile="annulus",
                  forcingWavenumber=2, forcingWidth=.7,
                  forcingAmplitude=.01, randomSeed=7654, threadCount=1,
                  dataDirectory=folder / "data", outputDirectory=folder / "output")
    values.update(overrides)
    path = folder / "run.params"
    path.write_text("".join(f"{key} {value}\n" for key, value in values.items()))
    return path


def run(command, params):
    result = subprocess.run([*command, str(params)], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(f"solver failed:\n{result.stdout}\n{result.stderr}")
    return result


def run_error(command, params, expected):
    result = subprocess.run([*command, str(params)], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode != 0 and expected in result.stderr, (
        f"expected failure containing {expected!r}:\n{result.stdout}\n{result.stderr}")


def checkpoint(folder, frame=2):
    payload = (folder / "data" / f"checkpoint_{frame:08d}.bin").read_bytes()
    magic, nx, ny, count = struct.unpack_from("8sQQQ", payload)
    assert magic.startswith(b"GP2DCP1") and count == nx * ny
    return struct.unpack_from(f"{2*count}d", payload, 32)


def relative_difference(a, b):
    return max(abs(x-y) for x, y in zip(a, b)) / max(1.0, max(map(abs, a)))


def mode(values, nx, x, y):
    offset = 2 * (y * nx + x)
    return complex(values[offset], values[offset + 1])


def signed_wave(index, count):
    return index if index <= (count - 1) // 2 else index - count


def initial_field(path, nx, ny):
    with path.open("w") as stream:
        for y in range(ny):
            row = []
            for x in range(nx):
                angle_x = 2 * math.pi * x / nx
                angle_y = 2 * math.pi * y / ny
                value = (0.7 + 0.15 * complex(math.cos(angle_x), math.sin(angle_x))
                         + 0.1 * complex(math.cos(2*angle_y), -math.sin(2*angle_y))
                         + 0.08 * complex(math.cos(2*angle_x-angle_y),
                                          math.sin(2*angle_x-angle_y)))
                row.extend((f"{value.real:.17g}", f"{value.imag:.17g}"))
            stream.write(" ".join(row) + "\n")


def restart_constraints_checks(root, command, tolerance=5e-11):
    initial = root / "restart_initial.dat"
    initial_field(initial, 12, 14)
    changed = root / "restart_changed_hypo"
    common = dict(numberOfSteps=1, outputIntervalSteps=1, forcingEnabled="false",
                  initialConditionFile=initial, threadCount=2)
    run(command, write_params(changed, **common))
    saved = changed / "data/wavefunction_00000001.dat"
    constrained = common | dict(hypoviscosity=.03, hypoviscosityOrder=-1)
    run(command, write_params(changed, **constrained))
    fresh = root / "fresh_changed_hypo"
    run(command, write_params(fresh, **(constrained | {"initialConditionFile": saved})))
    error = relative_difference(checkpoint(changed, 2), checkpoint(fresh, 1))
    assert error < tolerance, f"restart initialization constraints differ: {error}"

    full = root / "missing_initial_full"
    split = root / "missing_initial_split"
    run(command, write_params(full, **(common | {"numberOfSteps": 2})))
    params = write_params(split, **common)
    run(command, params)
    initial.unlink()
    run(command, params)
    assert checkpoint(full, 2) == checkpoint(split, 2), "restart requires its original IC"


def normalized_forcing_checks(root, cpu):
    for profile in ("annulus", "gaussian", "exponential", "logNormal"):
        reference = None
        for name, amplitude in (("unit", 1), ("large", "1e200"), ("small", "1e-200")):
            case = root / f"normalized_{profile}_{name}"
            run(cpu, write_params(case, forcingProfile=profile, forcingAmplitude=amplitude,
                                  targetWaveActionInjectionRate=1))
            values = checkpoint(case)
            if reference is None:
                reference = values
            assert values == reference, "normalized forcing depends on the common amplitude"
            assert max(map(abs, values)) > 0, "normalized forcing silently produced zero noise"
            with (case / "output/forcing_summary.csv").open() as stream:
                summary = next(csv.DictReader(stream))
            assert math.isclose(float(summary["wave_action_injection_coefficient"]), 1, rel_tol=1e-14)
            assert math.isfinite(float(summary["quadratic_energy_injection_coefficient"]))
    overflowing = write_params(root / "forcing_overflow", forcingAmplitude="1e200")
    run_error(cpu, overflowing, "forcing injection coefficients overflow or underflow")


def linear_checks(root, command):
    initial = root / "linear_initial.dat"
    initial_field(initial, 12, 14)
    for method in ("etd2", "etd3", "etd4", "rk2"):
        case = root / f"linear_{method}"
        run(command, write_params(case, integrator=method, initialConditionFile=initial,
                                  forcingEnabled="false", nonlinearityCoefficient=0,
                                  numberOfSteps=5, outputIntervalSteps=5))
        before, after = checkpoint(case, 0), checkpoint(case, 1)
        for y in range(14):
            for x in range(12):
                k2 = (signed_wave(x, 12) / 1.5) ** 2 + signed_wave(y, 14) ** 2
                linear = complex(-1e-6 * k2 ** 2, -(k2 + .1))
                expected = mode(before, 12, x, y) * cmath.exp(linear * 5 * 2e-4)
                assert abs(mode(after, 12, x, y) - expected) < 5e-13


def cpu_checks(root, cpu):
    full = root / "full"
    run(cpu, write_params(full, writeVortexDiagnostics="true"))
    with (full / "output/diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [row["frame"] for row in rows] == ["1", "2"]
    assert all(math.isfinite(float(row["total_energy"])) for row in rows)
    with (full / "output/vortices.csv").open() as stream:
        vortex_rows = list(csv.DictReader(stream))
    assert [row["frame"] for row in vortex_rows] == ["1", "2"]
    assert all(int(row["positive_vortices"]) >= 0 and
               int(row["negative_vortices"]) >= 0 for row in vortex_rows)
    with (full / "output/vortex_positions.csv").open() as stream:
        position_rows = list(csv.DictReader(stream))
    for vortex_row in vortex_rows:
        frame_positions = [row for row in position_rows
                           if row["frame"] == vortex_row["frame"]]
        assert len(frame_positions) == int(vortex_row["total_vortices"])
    assert all(int(row["circulation"]) in (-1, 1) and
               0 <= float(row["x"]) < 1.5 * 2 * math.pi and
               0 <= float(row["y"]) < 2 * math.pi and
               math.isfinite(float(row["core_residual"]))
               for row in position_rows)
    assert "kinetic_energy" in rows[0]
    assert all(math.isfinite(float(row["expected_full_energy_injection"]))
               for row in rows)
    values = checkpoint(full, 1)
    nx, ny, aspect = 12, 14, 1.5
    minimum_x, maximum_x = -(nx // 2), (nx - 1) // 2
    minimum_y, maximum_y = -(ny // 2), (ny - 1) // 2
    forcing = {}
    quadratic_injection = 0.0
    for y in range(ny):
        wave_y = signed_wave(y, ny)
        for x in range(nx):
            wave_x = signed_wave(x, nx)
            k2 = (wave_x / aspect) ** 2 + wave_y ** 2
            amplitude = (.01 if math.sqrt(k2) > 0 and
                         abs(math.sqrt(k2) - 2) < .7 else 0)
            forcing[wave_x, wave_y] = amplitude
            quadratic_injection += (k2 + .1) * amplitude ** 2
    quartic_injection = 0.0
    total_power = sum(amplitude ** 2 for amplitude in forcing.values())
    for iy in range(minimum_y, maximum_y + 1):
        for ix in range(minimum_x, maximum_x + 1):
            accepted_power = 0.0
            for (jx, jy), amplitude in forcing.items():
                if abs(jx - ix) < nx / 2 and abs(jy - iy) < ny / 2:
                    accepted_power += amplitude ** 2
            quartic_injection += (total_power + accepted_power) * abs(
                mode(values, nx, ix % nx, iy % ny)) ** 2
    area = 1.5 * (2 * math.pi) ** 2
    expected_injection = area * (quadratic_injection + 2 * quartic_injection)
    assert math.isclose(float(rows[0]["expected_full_energy_injection"]),
                        expected_injection, rel_tol=2e-11)
    for quantity in ("total_energy_dissipation_hypo",
                     "total_energy_dissipation_hyper"):
        assert all(math.isfinite(float(row[quantity])) for row in rows)
    for name in ("spectra.csv", "forcing_summary.csv"):
        header = (full / "output" / name).read_text().splitlines()[0]
        assert "quadratic_energy" in header
    flux_path = full / "output/fluxes.csv"
    with flux_path.open() as stream:
        fluxes = csv.DictReader(stream)
        assert fluxes.fieldnames == [
            "time", "frame", "wavenumber", "wave_action_flux",
            "full_energy_flux", "segment_mean_wave_action_flux",
            "segment_mean_full_energy_flux"]
        flux_rows = list(fluxes)
    checked_flux_columns = ("wave_action_flux", "full_energy_flux")
    assert all(math.isfinite(float(row[quantity]))
               for row in flux_rows for quantity in checked_flux_columns)
    for diagnostic in rows:
        frame_rows = [row for row in flux_rows
                      if row["frame"] == diagnostic["frame"]]
        scale = max(1.0, *(abs(float(row["full_energy_flux"]))
                           for row in frame_rows))
        assert abs(float(frame_rows[-1]["full_energy_flux"])) < 2e-10 * scale

    # For a spatially constant field and order-zero damping, the full rate has
    # the closed form 2*A*nu*(mu*|psi|^2 + g*|psi|^4).
    damping = root / "constant_damping"
    damping.mkdir()
    damping_initial = damping / "initial.dat"
    damping_initial.write_text((" ".join(["0.4 0.2"] * 12) + "\n") * 14)
    damping_params = write_params(
        damping, numberOfSteps=1, outputIntervalSteps=1,
        initialConditionFile=damping_initial, forcingEnabled="false",
        hyperviscosity=.2, hyperviscosityOrder=0, hypoviscosity=0)
    run(cpu, damping_params)
    with (damping / "output/diagnostics.csv").open() as stream:
        damping_row = next(csv.DictReader(stream))
    assert float(damping_row["expected_full_energy_injection"]) == 0
    amplitude = mode(checkpoint(damping, 1), 12, 0, 0)
    density = abs(amplitude) ** 2
    area = 1.5 * (2 * math.pi) ** 2
    expected = 2 * area * .2 * (.1 * density + 2 * density ** 2)
    assert math.isclose(float(damping_row["total_energy_dissipation_hyper"]),
                        expected, rel_tol=2e-11)
    assert (float(damping_row["total_energy_dissipation_hyper"])
            > float(damping_row["quadratic_energy_dissipation_hyper"]))
    split = root / "split"
    params = write_params(split, numberOfSteps=2)
    run(cpu, params)
    run(cpu, params)
    assert checkpoint(full) == checkpoint(split)
    assert len(list((split / "output/segments").glob("segment_*"))) == 2

    # An uncommitted frame is rolled back to its exact CSV offsets before the
    # step is repeated from the committed checkpoint.
    recovery = root / "recovery"
    recovery_params = write_params(recovery, numberOfSteps=2)
    run(cpu, recovery_params)
    csv_names = ("diagnostics.csv", "spectra.csv", "fluxes.csv", "modes.csv",
                 "vortices.csv", "vortex_positions.csv")
    records = []
    for name in csv_names:
        path = recovery / "output" / name
        exists = path.exists()
        size = path.stat().st_size if exists else 0
        records.append((exists, size))
        if exists:
            with path.open("ab") as stream:
                stream.write(b"interrupted row\n")
    (recovery / "data/wavefunction_00000002.dat").write_text("partial\n")
    (recovery / "data/checkpoint_00000002.bin").write_bytes(b"partial")
    journal = ["gp2d_output_transaction_v3",
               f'"{(recovery / "output").resolve()}"', "2 1 1"]
    journal.extend(f"{int(exists)} {size}" for exists, size in records)
    journal.extend(("0", "0", "0"))
    (recovery / "data/output_transaction.txt").write_text("\n".join(journal) + "\n")
    recovered = run(cpu, recovery_params)
    assert "recovered interrupted output frame 2" in recovered.stdout
    assert checkpoint(full) == checkpoint(recovery)

    # Expensive FFTW planning can be amortized across runs through wisdom.
    wisdom = root / "wisdom"
    wisdom_file = wisdom / "fftw.wisdom"
    wisdom_params = write_params(
        wisdom, numberOfSteps=1, outputIntervalSteps=1, forcingEnabled="false",
        fftwPlanning="measure", fftwWisdomFile=wisdom_file)
    run(cpu, wisdom_params)
    assert wisdom_file.exists() and wisdom_file.stat().st_size > 0
    run(cpu, wisdom_params)

    failed_wisdom = write_params(
        root / "failed_wisdom", numberOfSteps=1, outputIntervalSteps=1,
        forcingEnabled="false", fftwWisdomFile=root / "absent/wisdom")
    result = subprocess.run([*cpu, str(failed_wisdom)], capture_output=True, text=True)
    assert result.returncode == 1 and "cannot export FFTW wisdom" in result.stderr
    assert "terminate called" not in result.stderr

    # Reject a restarted timestep that cannot change the saved double timestamp.
    stagnant = root / "stagnant_time"
    stagnant_params = write_params(stagnant, numberOfSteps=1, outputIntervalSteps=1,
                                   forcingEnabled="false")
    run(cpu, stagnant_params)
    metadata = stagnant / "data/restart_state.txt"
    lines = metadata.read_text().splitlines()
    lines[1] = "time 10000000000000000"
    metadata.write_text("\n".join(lines) + "\n")
    before = metadata.read_bytes()
    run_error(cpu, stagnant_params, "timeStep is too small")
    assert metadata.read_bytes() == before
    assert not (stagnant / "data/checkpoint_00000002.bin").exists()

    # Frame collisions are detected before any diagnostic append is committed.
    collision = root / "collision"
    collision_params = write_params(collision, numberOfSteps=2)
    run(cpu, collision_params)
    (collision / "data/wavefunction_00000002.dat").write_text("collision\n")
    run_error(cpu, collision_params, "refusing to overwrite output frame file")

    # The inverse-power operator is singular at k=0, so suppress the condensate
    # explicitly while damping nonzero low modes.
    mean = root / "inverse_hypo_mean"
    mean_initial = mean / "initial.dat"
    mean.mkdir()
    initial_field(mean_initial, 12, 14)
    mean_params = write_params(
        mean, numberOfSteps=1, outputIntervalSteps=1, forcingEnabled="false",
        dispersionCoefficient=0, nonlinearityCoefficient=0,
        chemicalPotential=0, hyperviscosity=0, hypoviscosity=.5,
        hypoviscosityOrder=-1, initialConditionFile=mean_initial)
    run(cpu, mean_params)
    initial_checkpoint = checkpoint(mean, 0)
    final_checkpoint = checkpoint(mean, 1)
    assert mode(initial_checkpoint, 12, 0, 0) == 0
    assert mode(final_checkpoint, 12, 0, 0) == 0
    assert (abs(mode(final_checkpoint, 12, 1, 0))
            < abs(mode(initial_checkpoint, 12, 1, 0)))
    restart_constraints_checks(root, cpu)
    normalized_forcing_checks(root, cpu)
    linear_checks(root, cpu)
    print("CPU outputs and exact stochastic restart passed")


def backend_checks(root, cpu, candidate, tolerance):
    initial = root / "initial.dat"
    initial_field(initial, 12, 14)
    boundary_initial = root / "boundary_initial.dat"
    with boundary_initial.open("w") as stream:
        for y in range(14):
            row = []
            for x in range(12):
                value = sum((complex(.02 * math.cos(.71 * i), .015 * math.sin(.39 * i))
                             * cmath.exp(2j * math.pi * (signed_wave(i % 12, 12) * x / 12
                                                        + signed_wave(i // 12, 14) * y / 14)))
                            for i in range(12 * 14))
                row.extend((f"{value.real:.17g}", f"{value.imag:.17g}"))
            stream.write(" ".join(row) + "\n")
    cases = [
        ("spectral_boundaries", dict(forcingEnabled="false", initialConditionFile=boundary_initial)),
        ("etd2_gaussian", dict(integrator="etd2", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("etd3_gaussian", dict(integrator="etd3", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("etd4_gaussian", dict(integrator="etd4", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("etd4_cuda_graph", dict(integrator="etd4", forcingProfile="gaussian",
                                 forcingWidth=.8, cudaGraphEnabled="true")),
        ("rk2_gaussian", dict(integrator="rk2", forcingProfile="gaussian",
                              forcingWidth=.8)),
        ("annulus", dict(forcingProfile="annulus")),
        ("exponential", dict(forcingProfile="exponential",
                             forcingShapeOrder=3)),
        ("log_normal", dict(forcingProfile="logNormal",
                            forcingLogWidth=.3)),
        ("single_mode", dict(forcingProfile="singleMode")),
        ("unforced", dict(forcingEnabled="false")),
        ("linear_cuda_graph", dict(forcingEnabled="false", nonlinearityCoefficient=0,
                                   cudaGraphEnabled="true")),
        ("linear_stochastic_cuda_graph", dict(nonlinearityCoefficient=0,
                                              cudaGraphEnabled="true")),
        ("linear_single_mode", dict(nonlinearityCoefficient=0, forcingProfile="singleMode")),
        ("ginzburg_hypo", dict(forcingEnabled="false",
                                ginzburgLandauDamping=.2,
                                ginzburgLandauCutoff=1,
                                hypoviscosity=.03,
                                hypoviscosityOrder=-1,
                                hypoviscosityCutoffEnabled="true",
                                hypoviscosityCutoff=3)),
    ]
    largest_difference = 0.0
    for name, settings in cases:
        reference = root / f"cpu_{name}"
        comparison = root / f"candidate_{name}"
        common = dict(numberOfSteps=3, outputIntervalSteps=3, threadCount=2,
                      initialConditionFile=initial)
        common.update(settings)
        run(cpu, write_params(reference, **common))
        run(candidate, write_params(comparison, **common))
        error = relative_difference(checkpoint(reference, 1),
                                    checkpoint(comparison, 1))
        largest_difference = max(largest_difference, error)
        assert error < tolerance, f"{name} backend mismatch: {error}"

    sustained_initial = root / "sustained_initial.dat"
    initial_field(sustained_initial, 96, 128)
    sustained = dict(nx=96, ny=128, aspectRatio=1.5, integrator="etd4",
                     timeStep=1e-4, numberOfSteps=64, outputIntervalSteps=64,
                     initialConditionFile=sustained_initial, forcingEnabled="false",
                     threadCount=2)
    reference = root / "cpu_sustained"
    comparison = root / "candidate_sustained"
    run(cpu, write_params(reference, **sustained))
    run(candidate, write_params(comparison, **sustained))
    error = relative_difference(checkpoint(reference, 1), checkpoint(comparison, 1))
    largest_difference = max(largest_difference, error)
    assert error < tolerance, f"sustained rectangular backend mismatch: {error}"

    # A backend must also restore its own spectral and random-generator state
    # exactly when a stochastic run is split across invocations.
    full = root / "candidate_restart_full"
    split = root / "candidate_restart_split"
    restart_settings = dict(numberOfSteps=4, outputIntervalSteps=2,
                            initialConditionFile=initial,
                            forcingProfile="annulus", threadCount=2)
    run(candidate, write_params(full, **restart_settings))
    split_params = write_params(split, **(restart_settings | {"numberOfSteps": 2}))
    run(candidate, split_params)
    run(candidate, split_params)
    assert checkpoint(full) == checkpoint(split), "candidate restart is not exact"
    restart_constraints_checks(root, candidate, tolerance)
    linear_checks(root, candidate)
    print("Integrators, forcing profiles, damping, and restart agree across backends "
          f"(maximum relative difference {largest_difference:.3g})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["cpu", "backend"])
    parser.add_argument("--cpu", required=True)
    parser.add_argument("--candidate")
    parser.add_argument("--mpiexec")
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpi-numproc-flag", default="-n")
    parser.add_argument("--tolerance", type=float, default=5e-11)
    args = parser.parse_args()
    if args.mode == "backend" and not args.candidate:
        parser.error("backend mode requires --candidate")
    if args.mpi_ranks < 1:
        parser.error("--mpi-ranks must be positive")
    cpu = [args.cpu]
    candidate = [args.candidate] if args.candidate else None
    if args.mpiexec:
        candidate = [args.mpiexec, args.mpi_numproc_flag, str(args.mpi_ranks), *candidate]
    with tempfile.TemporaryDirectory(prefix="gp2d-regression-") as temporary:
        if args.mode == "cpu":
            cpu_checks(Path(temporary), cpu)
        else:
            backend_checks(Path(temporary), cpu, candidate, args.tolerance)


if __name__ == "__main__":
    main()
