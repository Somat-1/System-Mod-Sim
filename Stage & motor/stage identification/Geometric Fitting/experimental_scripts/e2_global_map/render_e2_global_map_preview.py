#!/usr/bin/env python3
"""A first look at what E2 was run for: the error map over the whole travel.

Not the finished analysis - no drift correction, no weighting, no uncertainty
budget. Just the three views that say whether the 400-target seeded grid did
its job: the direction-averaged deviation along the travel, what is left after
a single linear gain, and that remainder folded onto one turn of the screw.

The fold is the reason the grid is irregular. 400 targets spread over 50 screw
phase bins give every part of one revolution five or more visits at different
places along the travel, so a periodic screw term separates from whatever else
varies with position.

Writes rendered/e2_global_map_preview.png.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e2_run

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parents[3] / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, style   # noqa: E402

OUT = Path(os.environ.get("E2_OUT", _HERE / "rendered"))
AVG, FIT, FOLD = "#161616", "#99483e", "#4c78a8"
LEAD_MM = 2.0
MINUS = "−"


def num(value, places=2):
    return f"{value:.{places}f}".replace("-", MINUS)


def thesis_style():
    style()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "figure.facecolor": "#ffffff",
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",
    })


def main():
    thesis_style()
    run = e2_run.load()
    plan = run.plan
    hold = run.hold()
    resid = (hold - plan.target_mm) * 1000.0
    ok = run.matched
    up = ok & (plan.state == "UP")
    dn = ok & (plan.state == "DOWN")

    common, iu, idn = np.intersect1d(plan.index[up], plan.index[dn],
                                     return_indices=True)
    x = plan.target_mm[up][iu]
    avg = 0.5 * (resid[up][iu] + resid[dn][idn])
    order = np.argsort(x)
    x, avg = x[order], avg[order]

    slope, intercept = np.polyfit(x - x.mean(), avg, 1)
    line = np.polyval((slope, intercept), x - x.mean())
    rest = avg - line

    fig = plt.figure(figsize=(11.0, 8.0))
    gs = fig.add_gridspec(3, 1, hspace=0.42)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_a)
    ax_c = fig.add_subplot(gs[2])

    ax_a.plot(x, avg, ".", ms=3.4, color=AVG)
    ax_a.plot(x, line, "-", lw=1.6, color=FIT,
              label=f"linear gain {num(slope * 1000, 1)} ppm")
    ax_a.axhline(0.0, color=INK2, lw=0.8)
    ax_a.set_ylabel("deviation [µm]")
    ax_a.legend(loc="lower left", fontsize=9, framealpha=1.0)
    ax_a.tick_params(labelbottom=False)
    ax_a.set_title(f"(a) direction-averaged deviation, {x.size} targets "
                   f"— {np.ptp(avg):.2f} µm peak to peak",
                   fontsize=11, loc="left")

    ax_b.plot(x, rest, ".", ms=3.4, color=AVG)
    ax_b.axhline(0.0, color=INK2, lw=0.8)
    ax_b.set_xlabel("commanded position [mm]")
    ax_b.set_ylabel("after the gain [µm]")
    ax_b.set_title(f"(b) what the gain does not explain — RMS "
                   f"{rest.std():.3f} µm, {np.ptp(rest):.2f} µm "
                   f"peak to peak", fontsize=11, loc="left")

    phase = np.mod(x, LEAD_MM) / LEAD_MM
    ax_c.plot(phase, rest, ".", ms=3.4, color=FOLD, alpha=0.75)
    bins = np.linspace(0, 1, 26)
    which = np.clip(np.digitize(phase, bins) - 1, 0, bins.size - 2)
    centre = 0.5 * (bins[:-1] + bins[1:])
    mean = np.array([rest[which == k].mean() if (which == k).any() else np.nan
                     for k in range(bins.size - 1)])
    sem = np.array([rest[which == k].std(ddof=1) / np.sqrt((which == k).sum())
                    if (which == k).sum() > 1 else np.nan
                    for k in range(bins.size - 1)])
    ax_c.errorbar(centre, mean, yerr=2 * sem, fmt="o-", ms=4.5, lw=1.4,
                  color=AVG, mfc="#ffffff", mec=AVG, mew=1.2, ecolor=AVG,
                  elinewidth=0.9, capsize=2,
                  label="binned mean ±2σ of the mean")
    ax_c.axhline(0.0, color=INK2, lw=0.8)
    ax_c.set_xlabel("position within one screw turn [revolutions]")
    ax_c.set_ylabel("after the gain [µm]")
    ax_c.legend(loc="upper right", fontsize=9, framealpha=1.0)
    amp = np.nanmax(mean) - np.nanmin(mean)
    ax_c.set_title(f"(c) the same points folded on the {LEAD_MM:.0f} mm lead "
                   f"— binned swing {amp:.3f} µm", fontsize=11,
                   loc="left")

    for ax in (ax_a, ax_b, ax_c):
        ax.set_yticklabels([num(v, 1) for v in ax.get_yticks()])

    fig.text(0.012, 0.012,
             "A first look only: no drift correction, no weighting, no "
             "uncertainty budget. Reversal is removed by the direction "
             "average, which is why it does not appear here.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "e2_global_map_preview.png", dpi=200)
    print("wrote", OUT / "e2_global_map_preview.png")
    print(f"gain {slope * 1000:+.2f} ppm; pk-pk {np.ptp(avg):.3f} um; "
          f"after gain RMS {rest.std():.3f} um; fold swing {amp:.3f} um")


if __name__ == "__main__":
    main()
