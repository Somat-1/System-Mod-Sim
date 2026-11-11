#!/usr/bin/env python3
"""Shared figure style for the friction chapter.

Section 1 of the F0b results-plot instructions, in one place so every figure
in the chapter looks the same: serif type, one meaning per colour, spines
left and bottom, findings in the caption rather than the plot.

    from thesis_style import apply, panel_letter, C, save

Colours carry meaning per block, and every distinction a reader needs is
also carried by line style or marker fill, so the figures survive being
printed in greyscale.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

FULL_WIDTH = 7.2        # in, full text width
HALF_WIDTH = 3.5

# One meaning per colour, across the whole chapter.
C = {
    "N": "#2b3640",          # noise floor, near-black
    "L": "#1a6478",          # ladder, teal
    "S": "#c48a33",          # slow ramps, ochre
    "D": "#4f8a6f",          # dwell, green
    "V": "#b5473a",          # speed, red
    "R": "#6c5a8e",          # reference pairs, purple
    "delta": "#2b3640",      # pair difference, dark and thick
    "ref": "#8a9096",        # reference lines, grey
    "grid": "#e5e7ea",
    "ink": "#2b3640",
    "muted": "#59636b",
}

RC = {
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.4,
    "lines.linewidth": 1.4, "lines.markersize": 4,
    "figure.facecolor": "white",
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.facecolor": "white",
}

# The exact axis labels, so the chapter uses one wording and one unit style.
LABEL = {
    "travel": r"commanded travel $s$ (µm)",
    "lost": r"lost motion $\ell = s - \Delta x$ (µm)",
    "lag": r"lag $e = x - u$ (µm)",
    "stroke": r"stroke $s$ (µm)",
    "fraction": r"fraction arriving $\Delta x / s$",
    "tau": r"averaging time $\tau$ (s)",
    "precision": r"precision of a window mean (nm)",
}


def apply() -> None:
    mpl.rcParams.update(RC)


def panel_letter(ax, s: str, dx: float = -0.10, dy: float = 1.03) -> None:
    """Bold (a), (b) … outside the axes at the top left. No panel titles."""
    ax.text(dx, dy, f"({s})", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def finish(ax, grid: bool = True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(grid, color=C["grid"], lw=0.4)
    ax.tick_params(length=3)


def ref_line(ax, value: float, text: str, axis: str = "y", side: str = "right",
             style=(0, (1, 2)), pad: float = 0.012, at: float | None = None) -> None:
    """A reference line with its label offset from it, inside the axes.

    Labels never sit on their line and never cross the data: horizontal
    lines are labelled at one end, vertical lines at the top.
    """
    if axis == "y":
        ax.axhline(value, color=C["ref"], lw=0.7, ls=style, zorder=1)
        x = 1 - pad if side == "right" else pad
        ax.text(x, value, text, transform=ax.get_yaxis_transform(),
                ha="right" if side == "right" else "left", va="bottom",
                fontsize=8, color=C["ref"])
    else:
        ax.axvline(value, color=C["ref"], lw=0.7, ls=style, zorder=1)
        ax.text(value, pad if at is None else at, text,
                transform=ax.get_xaxis_transform(), ha="left", va="bottom",
                fontsize=8, color=C["ref"], rotation=90)


def save(fig, out_dir: Path, name: str) -> list[Path]:
    """Vector PDF for the thesis, 300 dpi PNG for slides."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in (".pdf", ".png"):
        path = out_dir / f"{name}{suffix}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
        written.append(path)
    plt.close(fig)
    return written
