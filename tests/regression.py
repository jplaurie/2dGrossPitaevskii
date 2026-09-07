#!/usr/bin/env python3
import argparse
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


def cpu_checks(root, cpu):
    full = root / "full"
    run(cpu, write_params(full))
    with (full / "output/diagnostics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [row["frame"] for row in rows] == ["1", "2"]
    assert all(math.isfinite(float(row["total_energy"])) for row in rows)
    assert "kinetic_energy" in rows[0]
    for name in ("spectra.csv", "fluxes.csv", "forcing_summary.csv"):
        header = (full / "output" / name).read_text().splitlines()[0]
        assert "quadratic_energy" in header
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
    csv_names = ("diagnostics.csv", "spectra.csv", "fluxes.csv", "modes.csv")
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
    journal = ["gp2d_output_transaction_v1",
               f'"{(recovery / "output").resolve()}"', "2 1 1"]
    journal.extend(f"{int(exists)} {size}" for exists, size in records)
    journal.extend(("0", "0"))
    (recovery / "data/output_transaction.txt").write_text("\n".join(journal) + "\n")
    recovered = run(cpu, recovery_params)
    assert "recovered interrupted output frame 2" in recovered.stdout
    assert checkpoint(full) == checkpoint(recovery)

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
    print("CPU outputs and exact stochastic restart passed")


def backend_checks(root, cpu, candidate):
    initial = root / "initial.dat"
    initial_field(initial, 12, 14)
    cases = [
        ("etd2_gaussian", dict(integrator="etd2", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("etd3_gaussian", dict(integrator="etd3", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("etd4_gaussian", dict(integrator="etd4", forcingProfile="gaussian",
                               forcingWidth=.8)),
        ("rk2_gaussian", dict(integrator="rk2", forcingProfile="gaussian",
                              forcingWidth=.8)),
        ("annulus", dict(forcingProfile="annulus")),
        ("exponential", dict(forcingProfile="exponential",
                             forcingShapeOrder=3)),
        ("log_normal", dict(forcingProfile="logNormal",
                            forcingLogWidth=.3)),
        ("single_mode", dict(forcingProfile="singleMode")),
        ("unforced", dict(forcingEnabled="false")),
        ("ginzburg_hypo", dict(forcingEnabled="false",
                                ginzburgLandauDamping=.2,
                                ginzburgLandauCutoff=1,
                                hypoviscosity=.03,
                                hypoviscosityOrder=-1,
                                hypoviscosityCutoffEnabled="true",
                                hypoviscosityCutoff=3)),
    ]
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
        assert error < 5e-11, f"{name} backend mismatch: {error}"

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
    print("Integrators, forcing profiles, damping, and restart agree across backends")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["cpu", "backend"])
    parser.add_argument("--cpu", required=True)
    parser.add_argument("--candidate")
    parser.add_argument("--mpiexec")
    parser.add_argument("--mpi-numproc-flag", default="-n")
    args = parser.parse_args()
    if args.mode == "backend" and not args.candidate:
        parser.error("backend mode requires --candidate")
    cpu = [args.cpu]
    candidate = [args.candidate] if args.candidate else None
    if args.mpiexec:
        candidate = [args.mpiexec, args.mpi_numproc_flag, "2", *candidate]
    with tempfile.TemporaryDirectory(prefix="gp2d-regression-") as temporary:
        if args.mode == "cpu":
            cpu_checks(Path(temporary), cpu)
        else:
            backend_checks(Path(temporary), cpu, candidate)


if __name__ == "__main__":
    main()
