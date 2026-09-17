"""Marker return-to-rest residual: does the axis land back where it started
after each leg's preceding marker, or does the reversal leave it somewhere
else?

Each leg is preceded by a marker: a fast move to a leg-specific negative
offset, held ~1 s, then a fast return, held ~0.5 s (see the firmware's
MARKER_REVERSE_DWELL_MS / MARKER_SETTLE_MS). This is a small, self-contained
out-and-back, structurally identical to the 10 mm trajectory legs themselves
but far smaller and with a direction reversal in the middle -- and direction
reversal is exactly where backlash / detent-torque effects should show up.

For each pair of consecutive legs, this compares:
  - rest_after[i]   -- the genuine dwell right after leg i's own ramp ends
                        (measured in leg i's own ramp_extent window), to
  - rest_before[i+1] -- the genuine dwell right before leg i+1's ramp starts
                        (measured in leg i+1's own ramp_extent window).
Both are local, ramp-bounded measurements (see leg_rest_levels in
plot_v6_exp_run.py) -- neither reaches into the marker itself, so the gap
between them is entirely what happened during that leg's own marker.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ids_common import (ANALYSIS_DIR, apply_style, load_ids, still_mask,
                        C_MEASURED, C_COMMANDED, C_INK, C_MUTED, C_GRID)
from plot_v6_exp_run import CSV, find_trajectories, leg_rest_levels, MRES_VALUES

OUT_DIR = ANALYSIS_DIR / "exp_run"


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t, pos, meta = load_ids(CSV)
    dt_s = meta["dt_s"]
    spans = find_trajectories(t, pos, dt_s)
    still = still_mask(pos, dt_s)

    _idx_b, rest_b, _idx_a, rest_a = leg_rest_levels(t, pos, spans, still, dt_s)
    residual = rest_b[1:] - rest_a[:-1]        # um, one per leg-to-leg gap
    transitions = np.arange(1, len(spans))     # gap i is "leg i -> leg i+1"
    mres_of_gap = [MRES_VALUES[i // 6] for i in range(len(spans) - 1)]

    fig, ax = plt.subplots(figsize=(13, 5))
    colors = [C_COMMANDED if m == 1 else C_MEASURED for m in mres_of_gap]
    ax.bar(transitions, residual, color=colors, width=0.7, zorder=3)
    ax.axhline(0.0, color=C_INK, lw=0.8, zorder=2)
    for g in range(6, len(spans), 6):
        ax.axvline(g + 0.5, color=C_MUTED, lw=0.8, ls=":", zorder=1)

    for x, m in zip([3.5, 9.5, 15.5, 21.5], MRES_VALUES):
        ax.text(x, ax.get_ylim()[1] * 0.92, f"MRES {m}", ha="center",
               fontsize=9, color=C_INK)

    ax.set_xlabel("Transition (leg i -> leg i+1)")
    ax.set_ylabel("Marker return residual [µm]\n(rest before next leg) - (rest after this leg)")
    ax.set_title("EXPrun — does each leg's marker return to where it started?",
                loc="left", pad=28)
    ax.text(0.0, 1.03,
            "orange = MRES 1 transitions, blue = MRES 4/16/32; "
            "large bars mean the reversal inside that marker left the axis "
            "somewhere else",
            transform=ax.transAxes, fontsize=8.5, color=C_MUTED, va="bottom")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "marker_return_residual.png", dpi=150,
               bbox_inches="tight")
    plt.close(fig)
    print("wrote marker_return_residual.png")

    print("\ngap   MRES  residual[um]")
    for i, (r, m) in enumerate(zip(residual, mres_of_gap)):
        print(f"{i+1:2d}->{i+2:<2d}  {m:4d}  {r:+8.2f}")
    mres1 = np.array([r for r, m in zip(residual, mres_of_gap) if m == 1])
    rest = np.array([r for r, m in zip(residual, mres_of_gap) if m != 1])
    print(f"\nMRES 1 transitions:      mean {mres1.mean():+.1f} um, "
          f"sd {mres1.std(ddof=1):.1f} um, range [{mres1.min():+.1f}, {mres1.max():+.1f}]")
    print(f"MRES 4/16/32 transitions: mean {rest.mean():+.2f} um, "
          f"sd {rest.std(ddof=1):.2f} um, range [{rest.min():+.2f}, {rest.max():+.2f}]")


if __name__ == "__main__":
    main()
