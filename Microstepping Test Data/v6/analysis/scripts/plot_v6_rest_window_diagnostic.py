"""Diagnose the resting-level bug that produced the spurious ~0.07 mm travel
scatter (and the apparent ~240 um MRES-1-only leg-to-leg swing) in an earlier
version of `travel_and_drift.png`.

That version computed each leg's resting level as the median of every
"still" (low-slope) sample in the WHOLE gap between the previous leg's end
and this leg's start. That gap is not a plain dwell: it also contains the
marker that precedes every leg (a fast move to a leg-specific offset, held
for its own ~1 s dwell, then a fast return) and, at the start of each MRES
group, the CONFIG/OSCILLATION markers plus the full 30 s oscillation block.
Depending on how many still samples fell in the marker's dwell versus the
genuine pre-ramp rest, the median landed on whichever cluster happened to be
larger -- essentially at random from one leg to the next.

The fix (now in `plot_v6_exp_run.py`) uses each leg's own `ramp_extent`
window instead: bounded to where motion has genuinely stopped on either side
of the ramp, so it cannot reach back into an unrelated marker dwell.

This script renders both windows on the same axes for the six MRES 1 legs
(where the bug was most visible), so the difference is visible in the raw
signal rather than asserted. Legs 3, 4 and 6 happen to show ~0 um
old-vs-new difference (their gaps are short enough that the marker dwell
doesn't dominate); legs 1, 2 and 5 show the effect clearly, leg 5 by a full
239 um.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ids_common import (ANALYSIS_DIR, apply_style, load_ids, still_mask,
                        ramp_extent, C_MEASURED, C_INK, C_MUTED, C_GRID)
from plot_v6_exp_run import CSV, find_trajectories

OUT_DIR = ANALYSIS_DIR / "exp_run"
N_LEGS_SHOWN = 6           # the MRES 1 group: where the bug was most visible.

C_OLD = "#d94801"          # old (whole-gap) window and its resulting median.
C_NEW = "#178a5f"          # new (ramp-local) window and its resulting median.


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    pos_mm = pos / 1000.0
    spans = find_trajectories(t, pos, dt_s)
    still = still_mask(pos, dt_s)

    fig, axes = plt.subplots(2, 3, figsize=(16, 7.5))
    labels = ["D-slow", "D-mod", "D-fast", "I-slow", "I-mod", "I-fast"]

    for i, ax in enumerate(axes.ravel()):
        s0, _apex, _e0 = spans[i]
        e0_prev = spans[i - 1][2] if i else 0
        a, _b = ramp_extent(still, s0, spans[i][2], dt_s, pad_s=1.0)

        old_seg, old_still = pos[e0_prev:s0], still[e0_prev:s0]
        old_rest = (float(np.median(old_seg[old_still]))
                   if old_still.sum() >= 50 else float(np.median(old_seg))) / 1000.0
        new_seg, new_still = pos[a:s0], still[a:s0]
        new_rest = (float(np.median(new_seg[new_still]))
                   if new_still.sum() > 20 else float(np.median(new_seg))) / 1000.0

        lo, hi = e0_prev, s0
        tt, yy = t[lo:hi] - t[lo], pos_mm[lo:hi]
        ax.plot(tt, yy, color=C_MEASURED, lw=0.8, zorder=3)

        ax.axvspan(0, t[s0] - t[lo], color=C_OLD, alpha=0.08, lw=0,
                  label="Old window (whole inter-leg gap)" if i == 0 else None)
        ax.axvspan(t[a] - t[lo], t[s0] - t[lo], color=C_NEW, alpha=0.16, lw=0,
                  label="New window (ramp-local dwell only)" if i == 0 else None)

        ax.axhline(old_rest, color=C_OLD, lw=1.3, ls="--", zorder=4,
                  label="Old rest (median of whole gap)" if i == 0 else None)
        ax.axhline(new_rest, color=C_NEW, lw=1.3, ls="-", zorder=4,
                  label="New rest (median of ramp-local dwell)" if i == 0 else None)

        diff_um = (old_rest - new_rest) * 1000.0
        ax.set_title(f"#{i+1}  MRES 1  {labels[i]}"
                     f"\nold-new = {diff_um:+.1f} µm", loc="left", fontsize=9.5)
        ax.set_xlabel("Time since previous leg ended [s]")
        if i % 3 == 0:
            ax.set_ylabel("Position [mm]")
        ax.tick_params(labelsize=8)

    axes.ravel()[0].legend(loc="upper right", fontsize=7, frameon=True,
                           framealpha=0.9).get_frame().set_edgecolor(C_GRID)

    fig.suptitle("EXPrun — resting-level window diagnostic (MRES 1 group)"
                "\nshaded = sample window used; dashed/solid = resulting median; "
                "large old-new gaps show the old method landing on the marker's "
                "own dwell instead of the true rest",
                x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUT_DIR / "rest_window_diagnostic.png", dpi=150,
               bbox_inches="tight")
    plt.close(fig)
    print("wrote rest_window_diagnostic.png")


if __name__ == "__main__":
    main()
