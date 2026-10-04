#!/usr/bin/env python3
"""Run calibrated CPU, MPI, and CUDA 2D Gross--Pitaevskii benchmarks."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys


BACKENDS = {
    "serial": ("CPU serial", "gp2d_benchmark_cpu_serial"),
    "cpu": ("CPU/OpenMP", "gp2d_benchmark_cpu"),
    "mpi": ("MPI", "gp2d_benchmark_mpi"),
    "cuda": ("CUDA FP64", "gp2d_benchmark_cuda"),
    "cuda_mixed": ("CUDA mixed", "gp2d_benchmark_cuda_mixed"),
}


def output(command):
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def physical_cores():
    cores = set()
    for cpu in Path("/sys/devices/system/cpu").glob("cpu[0-9]*"):
        try:
            cores.add(((cpu / "topology/physical_package_id").read_text().strip(),
                       (cpu / "topology/core_id").read_text().strip()))
        except OSError:
            pass
    return len(cores) or (os.cpu_count() or 1)


def parse_result(stdout):
    line = next((item for item in reversed(stdout.splitlines())
                 if item.startswith("BENCHMARK,")), "")
    fields = next(csv.reader([line]), [])
    if len(fields) != 7:
        raise RuntimeError(f"could not parse benchmark output:\n{stdout}")
    return dict(zip(("record", "reported_backend", "nx", "ny", "steps", "seconds",
                     "seconds_per_step"), fields))


def run_once(backend, executable, resolution, steps, warmup, mpi_ranks, cpu_threads):
    threads = cpu_threads if backend == "cpu" else 1
    command = [str(executable), str(resolution), str(steps), str(warmup), str(threads)]
    if backend == "mpi":
        command = ["mpirun", "--bind-to", "core", "--map-by", "core", "-n",
                   str(mpi_ranks), *command]
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": str(threads), "OMP_PROC_BIND": "close",
                        "OMP_PLACES": "cores"})
    completed = subprocess.run(command, text=True, capture_output=True, env=environment)
    if completed.returncode:
        raise RuntimeError(f"benchmark failed ({' '.join(command)}):\n"
                           f"{completed.stdout}{completed.stderr}")
    return parse_result(completed.stdout)


def calibrate(backend, executable, resolution, warmup, mpi_ranks, cpu_threads,
              target_seconds, minimum_steps):
    steps = 1
    while True:
        result = run_once(backend, executable, resolution, steps, warmup,
                          mpi_ranks, cpu_threads)
        seconds = float(result["seconds"])
        if seconds >= .1 or steps >= 1_000_000:
            break
        steps *= 4
    return max(minimum_steps, min(1_000_000, math.ceil(target_seconds * steps / seconds)))


def system_metadata(repo, args):
    cpu_model = "unknown"
    try:
        match = re.search(r"^model name\s*:\s*(.+)$", Path("/proc/cpuinfo").read_text(),
                          re.MULTILINE)
        if match:
            cpu_model = match.group(1).strip()
    except OSError:
        pass
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(), "platform": platform.platform(),
        "cpu_model": cpu_model, "logical_cpus": os.cpu_count(),
        "physical_cores_detected": physical_cores(),
        "cpu_openmp_threads": args.cpu_threads, "mpi_ranks": args.mpi_ranks,
        "gpu": output(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]),
        "gpu_driver": output(["nvidia-smi", "--query-gpu=driver_version",
                              "--format=csv,noheader"]),
        "compiler": output(["c++", "--version"]).splitlines()[0],
        "mpi": output(["mpirun", "--version"]).splitlines()[0],
        "cuda_compiler": output(["nvcc", "--version"]).splitlines()[-1],
        "git_commit": output(["git", "-C", str(repo), "rev-parse", "HEAD"]),
        "git_dirty": output(["git", "-C", str(repo), "status", "--porcelain"])
                     not in ("", "unknown"),
        "resolutions": args.resolutions, "trials": args.trials,
        "target_timed_seconds_per_trial": args.target_seconds,
        "minimum_timed_steps": args.minimum_timed_steps,
        "warmup_steps": args.warmup_steps, "integrator": "ETD4-B",
        "dealiasing": "two-pass 3/2 padding",
        "mixed_precision": "FP64 state/integration; FP32 padded FFTs and cubic products",
        "timing_excludes": ["process and MPI startup", "CUDA context creation", "FFT planning",
                            "allocation", "coefficient construction and upload",
                            "state initialization and upload", "warm-up steps", "state download",
                            "diagnostics", "file output"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build/bench"))
    parser.add_argument("--output", type=Path, default=Path("benchmarks/results.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("benchmarks/system.json"))
    parser.add_argument("--backends", nargs="+", choices=BACKENDS, default=list(BACKENDS))
    parser.add_argument("--resolutions", nargs="+", type=int,
                        default=[64, 128, 256, 512, 1024])
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--target-seconds", type=float, default=.75)
    parser.add_argument("--warmup-steps", type=int, default=2)
    parser.add_argument("--minimum-timed-steps", type=int, default=3)
    parser.add_argument("--cpu-threads", type=int, default=physical_cores())
    parser.add_argument("--mpi-ranks", type=int, default=physical_cores())
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if (any(value <= 0 for value in [*args.resolutions, args.trials, args.warmup_steps,
                                     args.minimum_timed_steps, args.cpu_threads, args.mpi_ranks])
            or args.target_seconds <= 0):
        parser.error("resolutions, counts, and durations must be positive")

    repo = Path(__file__).resolve().parents[1]
    build = (repo / args.build_dir).resolve() if not args.build_dir.is_absolute() else args.build_dir
    executables = {}
    for backend in args.backends:
        executable = build / BACKENDS[backend][1]
        if not executable.is_file():
            parser.error(f"missing {executable}; build the benchmark targets first")
        if backend == "mpi" and not shutil.which("mpirun"):
            parser.error("mpirun is not available")
        executables[backend] = executable
    result_path = (repo / args.output).resolve() if not args.output.is_absolute() else args.output
    metadata_path = ((repo / args.metadata).resolve()
                     if not args.metadata.is_absolute() else args.metadata)
    if result_path.exists() and not args.overwrite:
        parser.error(f"{result_path} exists; pass --overwrite to replace it")
    result_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(system_metadata(repo, args), indent=2) + "\n")

    columns = ["backend", "resolution", "grid_points", "trial", "steps", "warmup_steps",
               "seconds", "seconds_per_step", "grid_points_per_second",
               "cpu_threads", "mpi_ranks"]
    with result_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for resolution in args.resolutions:
            for backend in args.backends:
                steps = calibrate(backend, executables[backend], resolution, args.warmup_steps,
                                  args.mpi_ranks, args.cpu_threads, args.target_seconds,
                                  args.minimum_timed_steps)
                print(f"{BACKENDS[backend][0]:>12} {resolution}x{resolution:<5} "
                      f"timed_steps={steps} warmup_steps={args.warmup_steps}", flush=True)
                for trial in range(1, args.trials + 1):
                    result = run_once(backend, executables[backend], resolution, steps,
                                      args.warmup_steps, args.mpi_ranks, args.cpu_threads)
                    step_time = float(result["seconds_per_step"])
                    writer.writerow({"backend": BACKENDS[backend][0], "resolution": resolution,
                                     "grid_points": resolution * resolution, "trial": trial,
                                     "steps": result["steps"], "warmup_steps": args.warmup_steps,
                                     "seconds": result["seconds"], "seconds_per_step": step_time,
                                     "grid_points_per_second": resolution * resolution / step_time,
                                     "cpu_threads": args.cpu_threads, "mpi_ranks": args.mpi_ranks})
                    stream.flush()
                    print(f"  trial {trial}: {step_time:.6g} s/step", flush=True)
    print(f"wrote {result_path}\nwrote {metadata_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
