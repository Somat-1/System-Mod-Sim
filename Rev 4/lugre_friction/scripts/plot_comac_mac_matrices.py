#!/usr/bin/env python3
"""Figure for the two MAC tables in comac_mode_extraction.md (Co-MAC Stage 2):
heatmaps of MAC_rot and MAC_trans (Set A x Set B), with the per-row best-matching
cell outlined rather than shown as a separate panel. Reads the already-saved
comac_mode_extraction_data.npz -- no re-solve of the eigenproblem.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

LUGRE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = LUGRE_DIR / "rendered_assets" / "temp" / "Co-MAC" / "Split Method"

PALETTE = ["#454040", "#605B51", "#D8D365", "#E6F082"]
BEST_MATCH_COLOR = "#48A111"

# Flagged slots (0-indexed [A-mode row, B-mode col]) -- degenerate/aliased matches called
# out by row: low-frequency subspace leakage (A-mode2 vs B-mode1, both matrices) and, in
# the translational matrix only, the high-frequency block where x_s/x_n carry too little
# modal energy to distinguish shapes (A-mode3..6 vs B-mode5..6).
FLAGGED_ROT = [(1, 0)]
FLAGGED_TRANS = [(1, 0), (2, 4), (2, 5), (3, 4), (3, 5), (4, 4), (4, 5), (5, 4), (5, 5)]


def draw_flagged_cell(ax, i: int, j: int, val: float) -> None:
    """Blank the cell to black, print the value flanked by a big X on each side
    (same light text color used elsewhere for dark/grey cells), and box the
    value/X group -- in place of the heatmap color/plain text used for the
    other (unflagged) cells."""
    text_color = PALETTE[3]
    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor="black",
                            edgecolor="none", zorder=1))
    ax.text(j - 0.30, i, "X", ha="center", va="center", fontsize=11,
             fontweight="bold", color=text_color, zorder=10)
    ax.text(j + 0.30, i, "X", ha="center", va="center", fontsize=11,
             fontweight="bold", color=text_color, zorder=10)
    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5,
             color=text_color, zorder=10)
    ax.add_patch(Rectangle((j - 0.44, i - 0.19), 0.88, 0.38, fill=False,
                            edgecolor=text_color, linewidth=1.0, zorder=10))


def main() -> None:
    d = np.load(DATA_DIR / "npz" / "comac_mode_extraction_data.npz")
    mac_rot, mac_trans = d["mac_rot"], d["mac_trans"]
    n = mac_rot.shape[0]

    a_labels = [f"A-mode{i+1}" for i in range(n)]
    b_labels = [f"B-mode{j+1}" for j in range(n)]

    cmap = LinearSegmentedColormap.from_list("comac_palette", PALETTE)

    panels = [
        (mac_rot, "Torsional Global MAC"),
        (mac_trans, "Linear/Translational Global MAC"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.3))
    im = None
    for k, (ax, (mat, title)) in enumerate(zip(axes.flat, panels)):
        im = ax.imshow(mat, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")
        ax.text(0.5, -0.06, title, transform=ax.transAxes, ha="center", va="top", fontsize=10)

        ax.set_xticks(np.arange(-0.5, n, 1.0), minor=True)
        ax.set_yticks(np.arange(-0.5, n, 1.0), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.6)
        ax.tick_params(which="minor", length=0)
        ax.tick_params(which="major", length=0)

        ax.set_xticks(np.arange(n))
        ax.set_xticklabels(b_labels, rotation=45, ha="left", fontsize=8)
        ax.xaxis.tick_top()
        if k == 0:
            ax.set_yticks(np.arange(n))
            ax.set_yticklabels(a_labels, fontsize=8)
        else:
            ax.set_yticks([])

        flagged = set(FLAGGED_ROT if k == 0 else FLAGGED_TRANS)

        for i in range(n):
            for j in range(n):
                if (i, j) in flagged:
                    continue
                val = mat[i, j]
                text_color = PALETTE[0] if val > 0.6 else PALETTE[3]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                         fontsize=7.5, color=text_color)

        for i in range(n):
            j = int(np.argmax(mat[i, :]))
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, fill=False,
                                    edgecolor=BEST_MATCH_COLOR, linewidth=2.5,
                                    zorder=5, clip_on=False))

        for i, j in flagged:
            draw_flagged_cell(ax, i, j, mat[i, j])

    fig.suptitle("Global MAC Matrices (Rotational vs. Translational Sub-Vectors)",
                 fontsize=12)
    fig.tight_layout(rect=[0.0, 0.06, 0.9, 0.93])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.68])
    fig.colorbar(im, cax=cbar_ax)

    out_path = DATA_DIR / "comac_mac_matrices.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
