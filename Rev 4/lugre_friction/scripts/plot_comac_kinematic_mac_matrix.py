#!/usr/bin/env python3
"""Figure for the Kinematically Scaled Co-MAC report: single-panel heatmap of
MAC_scaled (6-DOF, Set A x Set B), reading from the already-saved
comac_kinematic_scaled_data.npz -- no re-solve. Same visual conventions as
Split Method/comac_mac_matrices.png (palette, white gridlines, green
row-argmax border, subplot title below), but one unified 6x6 matrix instead
of two 4x4/2x2 sub-matrices.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

LUGRE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Kinematically Scaled"

PALETTE = ["#454040", "#605B51", "#D8D365", "#E6F082"]
BEST_MATCH_COLOR = "#48A111"


def main() -> None:
    d = np.load(DATA_DIR / "npz" / "comac_kinematic_scaled_data.npz")
    mac = d["mac_scaled"]
    n = mac.shape[0]

    a_labels = [f"A-mode{i+1}" for i in range(n)]
    b_labels = [f"B-mode{j+1}" for j in range(n)]

    cmap = LinearSegmentedColormap.from_list("comac_palette", PALETTE)

    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(mac, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")

    ax.set_xticks(np.arange(-0.5, n, 1.0), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1.0), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.6)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)

    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(b_labels, rotation=45, ha="left", fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(a_labels, fontsize=8)

    for i in range(n):
        for j in range(n):
            val = mac[i, j]
            text_color = PALETTE[0] if val > 0.6 else PALETTE[3]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5, color=text_color)

    for i in range(n):
        j = int(np.argmax(mac[i, :]))
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, fill=False,
                                edgecolor=BEST_MATCH_COLOR, linewidth=2.5,
                                zorder=5, clip_on=False))

    ax.text(0.5, -0.08, "Kinematically Scaled Global MAC (6-DOF)",
             transform=ax.transAxes, ha="center", va="top", fontsize=10)

    fig.suptitle("Global MAC Matrix -- Kinematically Scaled (6-DOF, unified)", fontsize=12)
    fig.tight_layout(rect=[0.0, 0.06, 0.88, 0.93])
    cbar_ax = fig.add_axes([0.90, 0.15, 0.03, 0.68])
    fig.colorbar(im, cax=cbar_ax)

    out_path = DATA_DIR / "comac_kinematic_mac_matrix.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
