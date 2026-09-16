#!/usr/bin/env python3
"""Create an MP4 or GIF of wave-action and full-energy fluxes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from gp2d_plotting import (
    available_frames,
    parameter_float,
    read_csv,
    read_parameters,
    repository_root,
    rows_for_frame,
    select_frames,
    spectral_abscissa,
    use_plot_style,
    write_animation,
)


def parse_arguments() -> argparse.Namespace:
    root = repository_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "output/fluxes.csv")
    parser.add_argument(
        "--parameters", type=Path, default=root / "output/resolved_parameters.txt"
    )
    parser.add_argument("--output", type=Path, default=root / "figures/fluxes.mp4")
    parser.add_argument("--frames", nargs="*", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument(
        "--x-axis", choices=["wavenumber", "frequency"], default="wavenumber"
    )
    parser.add_argument("--y-scale", choices=["linear", "symlog"], default="symlog")
    parser.add_argument("--fps", type=float, default=12.0)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument(
        "--codec",
        choices=["h264", "h265"],
        default="h264",
        help="ffmpeg video codec for MP4 output (ignored for GIF)",
    )
    parser.add_argument("--no-tex", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    use_plot_style(not args.no_tex, font_size=14)
    table = read_csv(args.input)
    frames = select_frames(
        available_frames(table),
        args.frames if args.frames else None,
        start=args.start,
        stop=args.stop,
        stride=args.stride,
    )
    if not frames:
        raise ValueError("the frame selection is empty")
    parameters = read_parameters(args.parameters) if args.parameters.exists() else {}
    coefficient = parameter_float(parameters, "dispersionCoefficient", -1.0)
    chemical_potential = parameter_float(parameters, "chemicalPotential", 0.0)

    def values(frame: int):
        rows = rows_for_frame(table, frame)
        x, x_label = spectral_abscissa(
            np.asarray(rows["wavenumber"]), args.x_axis, coefficient, chemical_potential
        )
        wave = np.asarray(rows["wave_action_flux"])
        energy = np.asarray(rows["full_energy_flux"])
        mask = np.isfinite(x) & np.isfinite(wave) & np.isfinite(energy) & (x > 0.0)
        order = np.argsort(x[mask])
        return (
            x[mask][order],
            wave[mask][order],
            energy[mask][order],
            x_label,
            float(rows["time"][0]),
        )

    datasets = [values(frame) for frame in frames]
    x_all = np.concatenate([entry[0] for entry in datasets])
    wave_all = np.concatenate([entry[1] for entry in datasets])
    energy_all = np.concatenate([entry[2] for entry in datasets])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    wave_line, = axes[0].plot([], [], color="C0")
    energy_line, = axes[1].plot([], [], color="C1")
    coordinate = r"\omega" if args.x_axis == "frequency" else "k"
    for axis, ylabel, data in zip(
        axes,
        [rf"$\Pi_N({coordinate})$", rf"$\Pi_H({coordinate})$"],
        [wave_all, energy_all],
    ):
        axis.axhline(0.0, color="0.25", linewidth=0.8)
        axis.set_xscale("log")
        axis.set_yscale(args.y_scale)
        axis.set_xlim(x_all.min(), x_all.max())
        limit = max(1.05 * np.max(np.abs(data)), np.finfo(float).eps)
        axis.set_ylim(-limit, limit)
        axis.set_xlabel(datasets[0][3])
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.2)
    title = fig.suptitle("")

    def update(index: int):
        x, wave, energy, _, time = datasets[index]
        wave_line.set_data(x, wave)
        energy_line.set_data(x, energy)
        title.set_text(rf"frame {frames[index]}, $t={time:.6g}$")
        return wave_line, energy_line, title

    movie = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=1000.0 / args.fps, blit=False
    )
    destination = write_animation(movie, args.output, args.fps, args.dpi, args.codec)
    plt.close(fig)
    print(f"wrote {destination}")


if __name__ == "__main__":
    main()
