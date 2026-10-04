# Backend benchmarks

This benchmark times complete ETD4-B timesteps on square grids. Each step performs four cubic
nonlinear evaluations, and each evaluation performs four complex FFTs on a 3/2-padded grid.

Construction, FFT planning, allocation, coefficient/state upload, and warm-up happen before the
timed interval. Synchronization is included at both timing boundaries. Diagnostics, output, and
the final GPU download are excluded.

```bash
cmake -S . -B build/bench -DCMAKE_BUILD_TYPE=Release
cmake --build build/bench --parallel
python3 benchmarks/run_benchmarks.py --overwrite
python3 benchmarks/plot_benchmarks.py
```

The checked-in 2,048 x 2,048 results were collected separately because a single serial timestep
already takes several seconds:

```bash
python3 benchmarks/run_benchmarks.py --resolutions 2048 --trials 3 \
  --minimum-timed-steps 1 --warmup-steps 1 \
  --output benchmarks/results_large.csv --metadata benchmarks/system_large.json --overwrite
python3 benchmarks/plot_benchmarks.py \
  --input benchmarks/results.csv benchmarks/results_large.csv
```

CPU/OpenMP and MPI default to one worker per detected physical core. MPI uses single-threaded
ranks; CPU serial is compiled without OpenMP or threaded FFTW. The mixed backend keeps spectral
state, integration, coefficients, forcing, and noise in FP64 while using FP32 for padded FFT fields
and pointwise cubic products. The plot reports absolute step time and speedup normalized by the
single-core CPU serial backend.

Raw trials are in `results.csv` and `results_large.csv`, with matching system/build metadata and
timing exclusions in `system.json` and `system_large.json`. The runner uses the Python standard
library; plot generation also requires Matplotlib.
