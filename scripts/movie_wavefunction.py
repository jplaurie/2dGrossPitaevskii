#!/usr/bin/env python3
"""Create an MP4 or GIF from physical-space wavefunction snapshots."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from gp2d_plotting import (
    discover_wavefunctions,
    parameter_float,
    read_csv,
    read_parameters,
    read_wavefunction,
    repository_root,
    select_frames,
    use_plot_style,
    write_animation,
)


QUANTITIES = {
    "intensity": (lambda psi: np.abs(psi) ** 2, r"$|\psi|^2$", "magma"),
    "real": (np.real, r"$\mathrm{Re}(\psi)$", "RdBu_r"),
    "imaginary": (np.imag, r"$\mathrm{Im}(\psi)$", "RdBu_r"),
    "phase": (np.angle, r"$\arg(\psi)$", "twilight"),
}


def parse_arguments() -> argparse.Namespace:
    root = repository_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "data")
    parser.add_argument(
        "--parameters", type=Path, default=root / "output/resolved_parameters.txt"
    )
    parser.add_argument(
        "--diagnostics", type=Path, default=root / "output/diagnostics.csv"
    )
    parser.add_argument("--output", type=Path, default=root / "figures/wavefunction.mp4")
    parser.add_argument("--quantity", choices=[*QUANTITIES, "all"], default="intensity")
    parser.add_argument("--frames", nargs="*", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--fps", type=float, default=12.0)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument(
        "--codec",
        choices=["h264", "h265"],
        default="h264",
        help="ffmpeg video codec for MP4 output (ignored for GIF)",
    )
    parser.add_argument("--vmin", type=float, help="override scale for a single quantity")
    parser.add_argument("--vmax", type=float, help="override scale for a single quantity")
    parser.add_argument(
        "--color-scale",
        choices=["global", "first", "dynamic"],
        default="global",
        help="derive fixed limits from all frames (default), the first frame, or each frame",
    )
    parser.add_argument("--no-tex", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    use_plot_style(not args.no_tex, font_size=14)
    files = discover_wavefunctions(args.data_dir)
    requested = args.frames if args.frames else None
    frames = select_frames(
        files, requested, start=args.start, stop=args.stop, stride=args.stride
    )
    if not frames:
        raise ValueError("the frame selection is empty")

    parameters = read_parameters(args.parameters) if args.parameters.exists() else {}
    aspect_ratio = parameter_float(parameters, "aspectRatio", 1.0)
    lx = parameter_float(parameters, "domainLengthX", 2.0 * np.pi * aspect_ratio)
    ly = parameter_float(parameters, "domainLengthY", 2.0 * np.pi)
    extent = (0.0, lx, 0.0, ly)
    times: dict[int, float] = {}
    if args.diagnostics.exists():
        diagnostics = read_csv(args.diagnostics)
        times = {
            int(row["frame"]): float(row["time"])
            for row in diagnostics
        }

    first = read_wavefunction(files[frames[0]])
    names = list(QUANTITIES) if args.quantity == "all" else [args.quantity]
    if args.quantity == "all":
        fig, axes_grid = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
        axes = list(axes_grid.flat)
    else:
        fig, axis = plt.subplots(figsize=(7, 6))
        axes = [axis]

    real_limit = max(float(np.max(np.abs(first.real))), float(np.max(np.abs(first.imag))))
    intensity_limit = float(np.max(np.abs(first) ** 2))
    if args.color_scale == "global":
        for frame in frames[1:]:
            psi = read_wavefunction(files[frame])
            real_limit = max(
                real_limit,
                float(np.max(np.abs(psi.real))),
                float(np.max(np.abs(psi.imag))),
            )
            intensity_limit = max(intensity_limit, float(np.max(np.abs(psi) ** 2)))
    real_limit = max(real_limit, np.finfo(float).eps)
    intensity_limit = max(intensity_limit, np.finfo(float).eps)
    images = []
    for axis, name in zip(axes, names):
        transform, label, cmap = QUANTITIES[name]
        if name == "phase":
            limits = (-np.pi, np.pi)
        elif name in {"real", "imaginary"}:
            limits = (-real_limit, real_limit)
        else:
            limits = (0.0, intensity_limit)
        if len(names) == 1:
            limits = (
                limits[0] if args.vmin is None else args.vmin,
                limits[1] if args.vmax is None else args.vmax,
            )
        image = axis.imshow(
            transform(first),
            origin="lower",
            extent=extent,
            interpolation="bilinear",
            cmap=cmap,
            vmin=limits[0],
            vmax=limits[1],
        )
        axis.set_title(label)
        axis.set_xlabel(r"$x$")
        axis.set_ylabel(r"$y$")
        axis.set_aspect("equal")
        fig.colorbar(image, ax=axis)
        images.append(image)

    title = fig.suptitle("")

    def update(index: int):
        frame = frames[index]
        psi = read_wavefunction(files[frame])
        for image, name in zip(images, names):
            image.set_data(QUANTITIES[name][0](psi))
            if args.color_scale == "dynamic" and name != "phase":
                if name == "intensity":
                    limit = max(float(np.max(np.abs(psi) ** 2)), np.finfo(float).eps)
                    limits = (0.0, limit)
                else:
                    component = psi.real if name == "real" else psi.imag
                    limit = max(float(np.max(np.abs(component))), np.finfo(float).eps)
                    limits = (-limit, limit)
                if len(names) == 1:
                    limits = (
                        limits[0] if args.vmin is None else args.vmin,
                        limits[1] if args.vmax is None else args.vmax,
                    )
                image.set_clim(*limits)
        time_text = rf", $t={times[frame]:.6g}$" if frame in times else ""
        title.set_text(rf"frame {frame}{time_text}")
        return [*images, title]

    movie = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=1000.0 / args.fps, blit=False
    )
    destination = write_animation(movie, args.output, args.fps, args.dpi, args.codec)
    plt.close(fig)
    print(f"wrote {destination}")


if __name__ == "__main__":
    main()
