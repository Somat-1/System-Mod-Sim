#!/usr/bin/env python3
"""Figure for the Kinematically Scaled Co-MAC report: unified 6-DOF COMAC bar
chart, reading from the already-saved comac_kinematic_scaled_data.npz -- no
re-solve. Same visual conventions as Split Method/comac_stage3.png (palette,
0.9 reference line, title below), but one unified 6-bar chart instead of two
separate rotational/translational panels.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

LUGRE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Kinematically Scaled"

PALETTE = ["#454040", "#605B51", "#D8D365", "#E6F082"]
THRESHOLD_COLOR = "#605B51"


def main() -> None:
    d = np.load(DATA_DIR / "npz" / "comac_kinematic_scaled_data.npz")
    labels = [str(x) for x in d["state_labels"]]
    vals = d["comac"][:, 3]

    cmap = LinearSegmentedColormap.from_list("comac_palette", PALETTE)

    fig, ax = plt.subplots(figsize=(7.5, 5.3))
    x = np.arange(len(labels))
    colors = [cmap(v) for v in vals]
    ax.bar(x, vals, color=colors, width=0.6, zorder=3)

    ax.axhline(0.9, color=THRESHOLD_COLOR, linewidth=1.3, linestyle="--", zorder=2)

    for xi, v in zip(x, vals):
        va = "bottom" if v < 0.9 else "top"
        offset = 0.02 if v < 0.9 else -0.02
        ax.text(xi, v + offset, f"{v:.2f}", ha="center", va=va, fontsize=9, color=PALETTE[0])

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0.0, 1.05)
    ax.set_yticks(np.arange(0.0, 1.01, 0.2))
    ax.set_ylabel("COMAC", fontsize=9.5)
    ax.grid(axis="y", color="#cccccc", linewidth=0.8, zorder=0)
    ax.text(0.5, -0.12, "Unified 6-DOF COMAC (Kinematically Scaled)",
             transform=ax.transAxes, ha="center", va="top", fontsize=10)

    fig.suptitle("Coordinate Modal Assurance Criterion (COMAC) -- Kinematically Scaled", fontsize=12)
    fig.tight_layout(rect=[0.0, 0.06, 1.0, 0.93])

    out_path = DATA_DIR / "comac_kinematic_comac.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
