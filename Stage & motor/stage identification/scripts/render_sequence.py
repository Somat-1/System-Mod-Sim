#!/usr/bin/env python3
"""Draw the stage identification sequence: the whole run, then each phase close up.

Reads the plan from `sequence_design.py` (not from the CSV, so the staircase can
be drawn as real microsteps rather than inferred from a label) and writes
`../sequence/stage_identification_sequence.{png,svg}`.

The top panel is the 42-minute run end to end — what the stage does, and when.
Below it each phase is drawn on its own time axis, because the phases differ in
scale by three orders of magnitude: a staircase microstep is 0.6 µm against a
20 mm gain-calibration move, and on one shared axis four of the five phases are
a flat line.

Palette is Okabe-Ito, checked rather than assumed: every pair clears ΔE 8 in
OKLab under simulated protanopia, deuteranopia and tritanopia (worst pair 8.2),
and every hue clears 3:1 against white. Phases are also separated by position and
by a direct label on each panel, so identity never rests on colour alone.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sequence_design as sd  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

INK = {
    "gain": "#0072B2",
    "staircase": "#D55E00",
    "settling": "#009E73",
    "repeatability": "#CC79A7",
    "drift": "#666666",
}
GRID = "#D7D7D7"
TEXT = "#222222"
MUTED = "#6B6B6B"

PHASE_ORDER = ["gain", "staircase", "settling", "repeatability"]
PHASE_SHORT = {
    "gain": "1 · Gain calibration",
    "staircase": "2 · Microstep staircase",
    "settling": "3 · Settling matrix",
    "repeatability": "4 · Bidirectional repeatability",
    "drift": "Drift check",
}


def trace(rows):
    """Expand the plan into a position-vs-time curve, one segment per row.

    A move is a ramp at its commanded speed then a hold; a staircase leg is the
    individual microsteps with their dwells, which is the whole shape of that
    phase and is lost if it is drawn as a ramp.
    """
    t_pts, p_pts, seg = [], [], []
    for r in rows:
        t0, start, target = r["t_start"], r["start"], r["target"]
        n_micro = int(r.get("micro_steps", 0) or 0)
        if n_micro > 0:
            dwell = float(r["micro_dwell_s"])
            delta = (target - start) / n_micro
            t, p = t0, start
            for k in range(n_micro):
                p_next = start + delta * (k + 1)
                t_pts += [t, t]              # the step itself: vertical
                p_pts += [p, p_next]
                t_pts += [t, t + dwell]      # then the dwell: horizontal
                p_pts += [p_next, p_next]
                t += dwell
                p = p_next
        else:
            move_t = abs(target - start) / r["speed"]
            t_pts += [t0, t0 + move_t, r["t_end"]]
            p_pts += [start, target, target]
        seg.append((r["phase"], t0, r["t_end"]))
    return np.array(t_pts), np.array(p_pts), seg


def mm(steps):
    return np.asarray(steps) / sd.FULL_STEPS_PER_MM


def shade_phases(ax, rows, alpha=0.08, label=False):
    """Tint the background by phase so the structure reads without the legend."""
    runs = []
    for r in rows:
        if runs and runs[-1][0] == r["phase"]:
            runs[-1][2] = r["t_end"]
        else:
            runs.append([r["phase"], r["t_start"], r["t_end"]])
    seen = set()
    for phase, t0, t1 in runs:
        ax.axvspan(t0 / 60, t1 / 60, color=INK[phase],
                   alpha=alpha if phase != "drift" else alpha * 2.2, lw=0)
        # Drift checks chop a phase into several runs; label only the first,
        # or the overview reads as though the phase restarts each time.
        if label and (t1 - t0) > 120 and phase not in seen:
            seen.add(phase)
            ax.text((t0 + t1) / 120, 0.965, PHASE_SHORT[phase].split(" · ")[0],
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=8.5, color=INK[phase], weight="bold")
    return runs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--mres", type=int, default=16)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    out = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    phases, rows, total = sd.build_sequence(args.seed, args.mres)
    summary = sd.summarise(phases, rows, total)
    t, p, _ = trace(rows)

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#999999",
                         "axes.labelcolor": TEXT, "text.color": TEXT,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig = plt.figure(figsize=(15.5, 10.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.45, 1, 1],
                          hspace=0.72, wspace=0.17,
                          left=0.065, right=0.985, top=0.9, bottom=0.105)

    # ---------------------------------------------------------- whole run
    ax = fig.add_subplot(gs[0, :])
    shade_phases(ax, rows, label=True)
    ax.plot(t / 60, mm(p), lw=0.9, color="#1A1A1A")
    for r in rows:
        if r["phase"] == "drift":
            ax.plot(r["t_end"] / 60, mm(r["target"]), "v", ms=6,
                    color=INK["drift"], zorder=5)
    ax.set_xlim(0, total / 60)
    ax.set_ylabel("position  [mm from centre]")
    ax.set_xlabel("time  [min]")
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_title(f"Stage identification sequence — {total/60:.0f} min, "
                 f"{len(rows)} commanded moves",
                 fontsize=13, pad=26, weight="bold", loc="left")
    ax.axhline(0, color=INK["drift"], lw=0.8, ls="--", alpha=0.7)

    handles = [Line2D([], [], color=INK[k], lw=6, alpha=0.45,
                      label=PHASE_SHORT[k]) for k in PHASE_ORDER]
    handles.append(Line2D([], [], color=INK["drift"], marker="v", ls="--",
                          lw=0.8, ms=6,
                          label=f"Drift check — reference every "
                                f"{sd.DRIFT_INTERVAL_S/60:.0f} min"))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.015),
              ncol=5, frameon=False, fontsize=8.6, handlelength=1.6,
              columnspacing=1.5)

    # ------------------------------------------------------ phase close-ups
    windows = {}
    for r in rows:
        if r["phase"] in PHASE_ORDER:
            lo, hi = windows.get(r["phase"], (r["t_start"], r["t_end"]))
            windows[r["phase"]] = (min(lo, r["t_start"]), max(hi, r["t_end"]))

    detail = {
        "gain": dict(zoom=None,
                     caption="whole electrical periods AND whole screw "
                             "revolutions,\nso phase error and screw "
                             "eccentricity both cancel"),
        "staircase": dict(zoom=None,
                          caption="one microstep at a time, forward then back;\n"
                                  "the gap between passes is the hysteresis"),
        "settling": dict(zoom=260.0,
                         caption="distance × speed × direction, randomised;\n"
                                 "dwell after each stop is where it is measured"),
        "repeatability": dict(zoom=200.0,
                              caption="each target approached from both sides,\n"
                                      "with a run-up so backlash resolves the same way"),
    }

    for i, key in enumerate(PHASE_ORDER):
        axd = fig.add_subplot(gs[1 + i // 2, i % 2])
        t0, t1 = windows[key]
        zoom = detail[key]["zoom"]
        if zoom:
            t1 = min(t1, t0 + zoom)
        if key == "staircase":
            # The staircase spans 8 full steps - 80 um - against a +/-5 mm
            # transit. On a shared mm axis it is a flat line, so this panel is
            # drawn relative to its own start, in micrometres.
            legs = [r for r in rows if r["phase"] == "staircase" and r["micro_steps"]][:2]
            base = legs[0]["start"]
            tt0 = legs[0]["t_start"]
            tt1 = legs[-1]["t_end"]
            sel = (t >= tt0) & (t <= tt1)
            axd.plot(t[sel] - tt0, (p[sel] - base) * 10.0, lw=1.0, color=INK[key])
            axd.set_xlim(0, tt1 - tt0)
            axd.set_ylabel("position  [µm from leg start]")
            for k in range(1, 3):
                axd.axhline(k * sd.ELECTRICAL_PERIOD_STEPS * 10.0, color=GRID,
                            lw=0.8, ls=":")
                axd.text(2, k * sd.ELECTRICAL_PERIOD_STEPS * 10.0,
                         f" electrical period {k}", fontsize=7, color=MUTED,
                         va="bottom")
            # Each microstep is 0.625 um on an 80 um axis - true, but sub-pixel,
            # so the staircase reads as a smooth ramp. An inset over the first
            # few seconds shows the discreteness the phase is named for.
            ins = axd.inset_axes([0.52, 0.12, 0.45, 0.4])
            w = 12.0
            m2 = (t >= tt0) & (t <= tt0 + w)
            ins.plot(t[m2] - tt0, (p[m2] - base) * 10.0, lw=1.1, color=INK[key])
            ins.set_title(f"first {w:.0f} s", fontsize=7, color=MUTED, pad=2)
            ins.tick_params(labelsize=6.5, colors=MUTED)
            ins.grid(True, color=GRID, lw=0.5, alpha=0.7)
            ins.set_xlabel("s", fontsize=6.5, labelpad=1)
            ins.set_ylabel("µm", fontsize=6.5, labelpad=1)
            t1 = tt1
            t0 = tt0
        else:
            sel = (t >= t0) & (t <= t1)
            axd.plot(t[sel] - t0, mm(p[sel]), lw=1.1, color=INK[key])
            for r in [r for r in rows if t0 <= r["t_start"] <= t1]:
                if r["micro_steps"]:
                    continue
                axd.plot((r["t_end"] - t0), mm(r["target"]), ".", ms=3.5,
                         color=INK[key], alpha=0.65)
            axd.set_xlim(0, t1 - t0)
        axd.grid(True, color=GRID, lw=0.6, alpha=0.8)
        axd.set_xlabel("time within phase  [s]")
        if key != "staircase":
            axd.set_ylabel("position  [mm]")
        info = summary["phases"][key]
        axd.set_title(f"{PHASE_SHORT[key]}   ·   {info['moves']} moves, "
                      f"{info['minutes']:.1f} min"
                      + (f"   (first {zoom:.0f} s)" if zoom else
                         ("   (one forward + reverse leg)" if key == "staircase" else "")),
                      fontsize=10, color=INK[key], weight="bold", loc="left")
        axd.text(0.0, -0.30, detail[key]["caption"], transform=axd.transAxes,
                 ha="left", va="top", fontsize=7.8, color=MUTED,
                 linespacing=1.45)

    fig.text(0.065, 0.004,
             f"{sd.FULL_STEPS_PER_MM} full steps/mm · 1 full step = 10 µm · "
             f"electrical period = {sd.ELECTRICAL_PERIOD_STEPS} full steps = 40 µm · "
             f"screw revolution = {sd.FULL_STEPS_PER_REV} full steps = 2 mm · "
             f"working range ±{sd.WORKING_HALF_SPAN_MM:.0f} mm of {sd.TRAVEL_MM:.0f} mm travel",
             fontsize=7.8, color=MUTED)

    png = out / "stage_identification_sequence.png"
    fig.savefig(png, dpi=170)
    fig.savefig(png.with_suffix(".svg"))
    plt.close(fig)
    print(f"wrote {png}")
    print(f"wrote {png.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
