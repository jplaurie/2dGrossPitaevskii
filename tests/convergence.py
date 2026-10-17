#!/usr/bin/env python3
"""Check temporal order for uniform and constrained multimode GP solutions."""

import argparse
import cmath
import math
import struct
import subprocess
import tempfile
from pathlib import Path


def checkpoint_field(path):
    payload = path.read_bytes()
    magic, nx, ny, count = struct.unpack_from("8sQQQ", payload)
    assert magic.startswith(b"GP2DCP1") and count == nx * ny
    values = struct.unpack_from(f"{2 * count}d", payload, 32)
    return [complex(values[2 * i], values[2 * i + 1]) for i in range(count)]


def run_case(executable, root, method, steps, final_time, amplitude, constrained=False):
    case = root / f"{'constrained' if constrained else 'uniform'}-{method}-{steps}"
    case.mkdir()
    initial = case / "initial.dat"
    with initial.open("w") as stream:
        for y in range(8):
            row = []
            for x in range(8):
                value = amplitude
                if constrained:
                    value = (0.4 * cmath.exp(2j * math.pi * x / 8)
                             + 0.3 * cmath.exp(2j * math.pi * y / 8)
                             + 0.2 * cmath.exp(2j * math.pi * (x + y) / 8))
                row.extend((f"{value.real:.17g}", f"{value.imag:.17g}"))
            stream.write(" ".join(row) + "\n")
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
        "hypoviscosity": 0.1 if constrained else 0,
        "hypoviscosityOrder": -1,
        # Suppress the singular mean without damping any nonzero mode.
        "hypoviscosityCutoffEnabled": "true",
        "hypoviscosityCutoff": 0.1,
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
    actual = checkpoint_field(case / "data/checkpoint_00000001.bin")
    if constrained:
        if actual[0] != 0:
            raise AssertionError("constrained mean mode is nonzero")
        return actual
    frequency = 0.2 + 1.3 * abs(amplitude) ** 2
    exact = amplitude * cmath.exp(-1j * frequency * final_time)
    return abs(actual[0] - exact)


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

            fields = [run_case(args.cpu, root, method, n, 0.5, amplitude, constrained=True)
                      for n in (10, 20, 40, 80)]
            differences = [max(abs(a - b) for a, b in zip(coarse, fine))
                           for coarse, fine in zip(fields, fields[1:])]
            observed = math.log(differences[-2] / differences[-1], 2)
            if not math.isfinite(observed) or observed < minimum:
                raise AssertionError(
                    f"{method} constrained order {observed:.3f}, expected at least {minimum:.2f}; "
                    f"differences={differences}"
                )
            print(f"{method} constrained: order={observed:.3f}, differences={differences}")


if __name__ == "__main__":
    main()
