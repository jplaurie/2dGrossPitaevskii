#!/usr/bin/env python3
"""Exercise HDF5 snapshots, export formats, and HDF5 initial conditions."""

import argparse
import math
import subprocess
import tempfile
from pathlib import Path


def run(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"command failed: {' '.join(map(str, command))}\n"
                           f"{result.stdout}\n{result.stderr}")


def numbers(path):
    return [float(value) for value in path.read_text().split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", required=True)
    parser.add_argument("--exporter", required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="gp2d-hdf5-") as temporary:
        root = Path(temporary)
        initial = root / "initial.dat"
        with initial.open("w") as stream:
            for y in range(6):
                row = []
                for x in range(8):
                    value = complex(0.4 + 0.03 * math.cos(2 * math.pi * x / 8),
                                    0.2 * math.sin(2 * math.pi * y / 6))
                    row.extend((f"{value.real:.17g}", f"{value.imag:.17g}"))
                stream.write(" ".join(row) + "\n")
        params = root / "run.params"
        params.write_text(
            "nx 8\nny 6\naspectRatio 1.3\ntimeStep 0.001\nnumberOfSteps 1\n"
            "outputIntervalSteps 1\nintegrator etd4\ndispersionCoefficient -1\n"
            "nonlinearityCoefficient 1\nchemicalPotential 0\nhyperviscosity 0\n"
            "hypoviscosity 0\nforcingEnabled false\nthreadCount 1\n"
            "fieldOutputFormat both\nhdf5CompressionLevel 1\n"
            f"initialConditionFile {initial}\ndataDirectory {root / 'data'}\n"
            f"outputDirectory {root / 'output'}\n"
        )
        run([args.cpu, params])
        for frame in (0, 1):
            stem = f"wavefunction_{frame:08d}"
            text_field = root / "data" / f"{stem}.dat"
            hdf5_field = root / "data" / f"{stem}.h5"
            exported = root / f"{stem}_exported.dat"
            assert text_field.exists() and hdf5_field.exists()
            run([args.exporter, hdf5_field, exported])
            expected, actual = numbers(text_field), numbers(exported)
            assert len(expected) == len(actual) == 2 * 8 * 6
            assert max(abs(a - b) for a, b in zip(expected, actual)) < 1e-12

        table = root / "gnuplot.dat"
        run([args.exporter, root / "data/wavefunction_00000001.h5", table,
             "--format", "gnuplot"])
        rows = [line.split() for line in table.read_text().splitlines()
                if line and not line.startswith("#")]
        assert len(rows) == 8 * 6 and all(len(row) == 6 for row in rows)
        run([args.cpu, params])
        assert (root / "data/wavefunction_00000002.h5").exists()

        restart_root = root / "from_hdf5"
        restart_params = root / "from_hdf5.params"
        restart_params.write_text(
            "nx 8\nny 6\naspectRatio 1.3\ntimeStep 0.001\nnumberOfSteps 1\n"
            "outputIntervalSteps 1\nintegrator etd2\ndispersionCoefficient -1\n"
            "nonlinearityCoefficient 0\nchemicalPotential 0\nhyperviscosity 0\n"
            "hypoviscosity 0\nforcingEnabled false\nthreadCount 1\n"
            "fieldOutputFormat hdf5\n"
            f"initialConditionFile {root / 'data/wavefunction_00000001.h5'}\n"
            f"dataDirectory {restart_root / 'data'}\n"
            f"outputDirectory {restart_root / 'output'}\n"
        )
        run([args.cpu, restart_params])
        assert (restart_root / "data/wavefunction_00000000.h5").exists()
        assert not (restart_root / "data/wavefunction_00000000.dat").exists()
        print("HDF5 snapshots, exports, compression, and initial-condition input passed")


if __name__ == "__main__":
    main()
