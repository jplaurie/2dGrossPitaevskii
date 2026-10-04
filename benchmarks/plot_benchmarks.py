#!/usr/bin/env python3
"""Plot complete-timestep timing and single-core-normalized speedup."""

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


ORDER = ("CPU serial", "CPU/OpenMP", "MPI", "CUDA FP64", "CUDA mixed")
COLORS = {"CPU serial": "#577590", "CPU/OpenMP": "#277da1", "MPI": "#f8961e",
          "CUDA FP64": "#43aa8b", "CUDA mixed": "#d1495b"}
MARKERS = {"CPU serial": "D", "CPU/OpenMP": "o", "MPI": "s",
           "CUDA FP64": "^", "CUDA mixed": "v"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, nargs="+",
                        default=[Path("benchmarks/results.csv")],
                        help="one or more benchmark CSV files")
    parser.add_argument("--metadata", type=Path, default=Path("benchmarks/system.json"))
    parser.add_argument("--output-prefix", type=Path, default=Path("benchmarks/backend_scaling"))
    args = parser.parse_args()
    rows = []
    for input_path in args.input:
        with input_path.open(newline="") as stream:
            rows.extend(csv.DictReader(stream))
    metadata = json.loads(args.metadata.read_text())
    samples = {}
    for row in rows:
        samples.setdefault(row["backend"], {}).setdefault(int(row["resolution"]), []).append(
            float(row["seconds_per_step"]))
    labels = {"CPU serial": "CPU serial (1 core)",
              "CPU/OpenMP": f"CPU/OpenMP ({metadata['cpu_openmp_threads']} threads)",
              "MPI": f"MPI ({metadata['mpi_ranks']} ranks)", "CUDA FP64": "CUDA FP64",
              "CUDA mixed": "CUDA mixed (FP64 state, FP32 FFT)"}
    figure, (timing, speedup) = plt.subplots(1, 2, figsize=(10.8, 4.3))
    medians = {}
    for backend in ORDER:
        if backend not in samples:
            continue
        resolutions = sorted(samples[backend])
        values = [median(samples[backend][n]) for n in resolutions]
        medians[backend] = dict(zip(resolutions, values))
        lower = [value - min(samples[backend][n]) for n, value in zip(resolutions, values)]
        upper = [max(samples[backend][n]) - value for n, value in zip(resolutions, values)]
        timing.errorbar(resolutions, values, yerr=[lower, upper], label=labels[backend],
                        color=COLORS[backend], marker=MARKERS[backend], linewidth=2, capsize=3)
    timing.set(xscale="log", yscale="log", xlabel="Grid resolution, N (for N x N)",
               ylabel="Time per ETD4 step (s)", title="Complete timestep time")
    timing.grid(True, which="both", alpha=.25)
    timing.legend(frameon=False, fontsize=8.5)
    baseline = medians.get("CPU serial", {})
    for backend in ORDER:
        if backend not in medians:
            continue
        resolutions = sorted(set(baseline) & set(medians[backend]))
        speedup.plot(resolutions, [baseline[n] / medians[backend][n] for n in resolutions],
                     label=labels[backend], color=COLORS[backend], marker=MARKERS[backend],
                     linewidth=2)
    speedup.axhline(1, color=".45", linewidth=1, linestyle="--")
    speedup.set(xscale="log", yscale="log", xlabel="Grid resolution, N (for N x N)",
                ylabel="Speedup over one CPU core", title="Backend speedup")
    speedup.grid(True, which="both", alpha=.25)
    resolutions = sorted({n for backend in samples.values() for n in backend})
    for axis in (timing, speedup):
        axis.set_xscale("log", base=2)
        axis.set_xticks(resolutions, labels=[f"{n:,}" for n in resolutions])
        axis.xaxis.set_minor_locator(NullLocator())
    figure.suptitle(f"2D Gross--Pitaevskii ETD4 benchmark\n{metadata['cpu_model']} · "
                    f"{metadata['gpu']}", fontsize=11)
    figure.tight_layout()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "png"):
        path = args.output_prefix.with_suffix(f".{extension}")
        figure.savefig(path, dpi=180, bbox_inches="tight")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
