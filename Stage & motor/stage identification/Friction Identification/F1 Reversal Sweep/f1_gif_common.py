#!/usr/bin/env python3
"""Shared parts of the three plain-language F1 GIFs (P, M and Q).

The briefs (F1_P_loops_gif_simple.md, F1_M_Q_gifs_simple.md) ask for one
visual grammar across all three: the same colours, the same font, the same
canvas, the same pace and the same footer. That only holds if the three
renderers share the pieces rather than each carrying its own copy, so the
palette, the canvas, the frame capture and the GIF writer live here.

Nothing in this module knows anything about friction; the physics comes from
render_f1_companions, which all three renderers import.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "companion"

# Colours, pace and canvas, from the briefs.
TEAL, ORANGE = "#1a6478", "#c48a33"
GHOST, INK, MUTED = "#d0d4d8", "#2b3640", "#59636b"
DOTTED = "#9aa0a6"           # the grey "what if" alternatives
FPS = 20
SIZE_PX = (1000, 560)
MAX_BYTES = 5_000_000
FOOTER = "illustration, not measured data"

RC = {
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
}


def new_figure():
    plt.rcParams.update(RC)
    return plt.figure(figsize=(SIZE_PX[0] / 100, SIZE_PX[1] / 100), dpi=100,
                      facecolor="white")


def footer(fig, extra: str = "") -> None:
    fig.text(0.075, 0.022, FOOTER + (f"  ·  {extra}" if extra else ""),
             fontsize=9, color=MUTED, va="center")


def lighten(colour: str, amount: float) -> tuple[float, float, float]:
    """Mix a colour towards white. 0 keeps it, 1 is white."""
    rgb = matplotlib.colors.to_rgb(colour)
    return tuple(c + (1.0 - c) * amount for c in rgb)


def capture(fig) -> Image.Image:
    fig.canvas.draw()
    return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())


def save_gif(frames: list[Image.Image], path: Path, fps: int = FPS) -> None:
    """One shared palette, so colours stay put and frame differences stay
    small. Falls back to narrower frames only if the file would be large;
    the briefs ask for about 1000 px wide."""
    path.parent.mkdir(parents=True, exist_ok=True)
    shades = [(v, v, v) for v in np.linspace(0, 255, 40).astype(int)]
    for base in ((26, 100, 120), (196, 138, 51), (208, 212, 216), (43, 54, 64)):
        shades.extend(tuple(round(c + (255 - c) * t) for c in base)
                      for t in np.linspace(0, 1, 24))
    flat = [int(c) for rgb in shades for c in rgb]
    flat += flat[-3:] * (256 - len(shades))      # a palette is exactly 256
    palette = Image.new("P", (1, 1))
    palette.putpalette(flat)
    for width in (SIZE_PX[0], 900, 800, 700, 600):
        size = (width, round(frames[0].height * width / frames[0].width))
        resized = (frames if width == frames[0].width else
                   [f.resize(size, Image.Resampling.LANCZOS) for f in frames])
        indexed = [f.convert("RGB").quantize(palette=palette,
                                             dither=Image.Dither.NONE)
                   for f in resized]
        indexed[0].save(path, save_all=True, append_images=indexed[1:],
                        duration=round(1000 / fps), loop=0, optimize=True)
        if path.stat().st_size <= MAX_BYTES:
            return
    raise RuntimeError(f"GIF stays above {MAX_BYTES / 1e6:.0f} MB")


def report(path: Path, frames: list, fps: int = FPS) -> None:
    print(f"{path.relative_to(HERE)}")
    print(f"  {len(frames)} frames, {len(frames) / fps:.1f} s at {fps} fps, "
          f"{path.stat().st_size / 1e6:.2f} MB")
