#!/usr/bin/env python3
"""Check temporal order against the exact spatially uniform GP solution."""

import argparse
import cmath
import math
import struct
import subprocess
import tempfile
from pathlib import Path


def checkpoint_mode(path):
    payload = path.read_bytes()
    magic, nx, ny, count = struct.unpack_from("8sQQQ", payload)
    assert magic.startswith(b"GP2DCP1") and count == nx * ny
    real, imaginary = struct.unpack_from("2d", payload, 32)
    return complex(real, imaginary)


def run_case(executable, root, method, steps, final_time, amplitude):
    case = root / f"{method}-{steps}"
    case.mkdir()
    initial = case / "initial.dat"
    row = " ".join([f"{amplitude.real:.17g} {amplitude.imag:.17g}"] * 8)
    initial.write_text((row + "\n") * 8)
    parameters = {
        "nx": 8,
        "ny": 8,
        "aspectRatio": 1,
        "timeStep": f"{final_time / steps:.17g}",
        "numberOfSteps": steps,
        "outputIntervalSteps": steps,
        "integrator": method,
        "dispersionCoefficient": -1,
        "nonlinearityCoefficient": 1.3,
        "chemicalPotential": 0.2,
        "hyperviscosity": 0,
        "hypoviscosity": 0,
        "ginzburgLandauDamping": 0,
        "forcingEnabled": "false",
        "threadCount": 1,
        "initialConditionFile": initial,
        "dataDirectory": case / "data",
        "outputDirectory": case / "output",
    }
    parameter_file = case / "run.params"
    parameter_file.write_text("".join(f"{key} {value}\n" for key, value in parameters.items()))
    result = subprocess.run([executable, parameter_file], text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"{method}/{steps} failed:\n{result.stdout}\n{result.stderr}")
    actual = checkpoint_mode(case / "data/checkpoint_00000001.bin")
    frequency = 0.2 + 1.3 * abs(amplitude) ** 2
    exact = amplitude * cmath.exp(-1j * frequency * final_time)
    return abs(actual - exact)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", required=True)
    args = parser.parse_args()
    expected_orders = {"rk2": 2, "etd2": 2, "etd3": 3, "etd4": 4}
    amplitude = complex(0.7, 0.2)
    steps = (4, 8, 16, 32)
    with tempfile.TemporaryDirectory(prefix="gp2d-convergence-") as temporary:
        root = Path(temporary)
        for method, expected_order in expected_orders.items():
            errors = [run_case(args.cpu, root, method, n, 0.8, amplitude) for n in steps]
            observed = math.log(errors[-2] / errors[-1], 2)
            minimum = expected_order - 0.35
            if not math.isfinite(observed) or observed < minimum:
                raise AssertionError(
                    f"{method} observed order {observed:.3f}, expected at least {minimum:.2f}; "
                    f"errors={errors}"
                )
            print(f"{method}: order={observed:.3f}, errors={errors}")


if __name__ == "__main__":
    main()
