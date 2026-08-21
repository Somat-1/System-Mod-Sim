#!/usr/bin/env python3
"""Figure for the Co-MAC Stage 3 section of comac_mode_extraction.md: per-DOF
COMAC bars (rotational vs. translational sub-systems), computed over the
filtered/sign-aligned paired modes from plot_comac_mode_extraction.py. Reads
the already-saved comac_mode_extraction_data.npz -- no re-solve.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

LUGRE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Split Method"

PALETTE = ["#454040", "#605B51", "#D8D365", "#E6F082"]
THRESHOLD_COLOR = "#605B51"


def main() -> None:
    d = np.load(DATA_DIR / "npz" / "comac_mode_extraction_data.npz")
    comac_rot, comac_trans = d["comac_rot"], d["comac_trans"]
    rot_labels = [str(x) for x in d["rot_labels"]]
    trans_labels = [str(x) for x in d["trans_labels"]]

    vals_rot = comac_rot[:, 3]
    vals_trans = comac_trans[:, 3]

    cmap = LinearSegmentedColormap.from_list("comac_palette", PALETTE)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 5.3),
                              gridspec_kw={"width_ratios": [4, 2]})

    panels = [(rot_labels, vals_rot, axes[0], "Rotational COMAC"),
              (trans_labels, vals_trans, axes[1], "Translational COMAC")]

    for labels, vals, ax, title in panels:
        x = np.arange(len(labels))
        colors = [cmap(v) for v in vals]
        bars = ax.bar(x, vals, color=colors, width=0.6, zorder=3)

        ax.axhline(0.9, color=THRESHOLD_COLOR, linewidth=1.3, linestyle="--", zorder=2)

        for xi, v in zip(x, vals):
            text_color = PALETTE[0] if v > 0.55 else PALETTE[2]
            va = "bottom" if v < 0.9 else "top"
            offset = 0.02 if v < 0.9 else -0.02
            ax.text(xi, v + offset, f"{v:.2f}", ha="center", va=va, fontsize=9,
                     color=PALETTE[0])

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylim(0.0, 1.05)
        ax.set_yticks(np.arange(0.0, 1.01, 0.2))
        ax.grid(axis="y", color="#cccccc", linewidth=0.8, zorder=0)
        ax.text(0.5, -0.14, title, transform=ax.transAxes, ha="center", va="top", fontsize=10)

    axes[0].set_ylabel("COMAC", fontsize=9.5)
    axes[1].set_yticklabels([])

    fig.suptitle("Coordinate Modal Assurance Criterion (COMAC) per Degree of Freedom",
                 fontsize=12)
    fig.tight_layout(rect=[0.0, 0.06, 1.0, 0.93])

    out_path = DATA_DIR / "comac_stage3.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
