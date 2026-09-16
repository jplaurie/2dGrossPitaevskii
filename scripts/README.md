# Plotting and movies

These notebooks and scripts read the current solver output format directly.
They use Matplotlib with serif fonts and external LaTeX text rendering, in the
same visual spirit as the legacy `2DGrossPitaevski` plotting scripts.

## Notebooks

- `wavefunction.ipynb` writes separate intensity, real-part, imaginary-part,
  and phase PDFs.
- `spectra.ipynb` writes separate wave-action and quadratic-energy spectrum
  PDFs against either wavenumber or angular frequency.
- `fluxes.ipynb` plots wave-action and full-energy fluxes.
- `diagnostics.ipynb` writes one PDF for conserved-quantity contributions and a
  second PDF for damping/injection rates.

Edit the clearly marked **Configuration** cell near the top of each notebook.
`MODE` can be `"single"`, `"multiple"`, or `"average"`; `FRAMES=None` uses all
available frames, while negative entries such as `FRAMES=[-1]` select from the
end.  Every notebook writes PDF files below `figures/` by default.

The frequency axis is derived from the solver's linear dispersion relation,
`omega(k) = abs(-c*k**2 + mu)`.  In `spectra.ipynb`, `DENSITY_IN_FREQUENCY=True`
also divides the spectral density by `abs(d omega/d k)`; the zero mode is
omitted because this Jacobian vanishes there.

Launch Jupyter from either the repository root or this directory:

```bash
jupyter notebook scripts/wavefunction.ipynb
```

The dependencies are Python 3, NumPy, Matplotlib, Jupyter, and a working LaTeX
installation.  Set `USE_TEX=False` in a notebook (or pass `--no-tex` to a movie
script) when LaTeX is unavailable.

## Movies

The movie scripts use the same readers and styling.  MP4 output requires
`ffmpeg`; GIF output uses Matplotlib's Pillow writer.

```bash
python scripts/movie_wavefunction.py --quantity all --frames 0 1 2 3 \
    --output figures/wavefunction.mp4

python scripts/movie_spectra.py --start 1 --stop 100 --stride 2 \
    --output figures/spectra.mp4

python scripts/movie_fluxes.py --start 1 --stop 100 --stride 2 \
    --codec h265 --output figures/fluxes.mp4
```

Omit `--frames`, `--start`, and `--stop` to animate every available frame.  Run
any script with `--help` for axis, scale, frame-rate, codec, and input-path
options.  MP4 defaults to H.264 (`--codec h264`); choose `--codec h265` for
HEVC/H.265 output.  The codec option is ignored when the output is a GIF.
Physical-state movies use one global color scale by default (requiring one
pre-scan of the selected snapshots); `--color-scale first` avoids that scan,
and `--color-scale dynamic` rescales every frame.
