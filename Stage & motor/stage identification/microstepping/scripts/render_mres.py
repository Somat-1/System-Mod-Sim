#!/usr/bin/env python3
"""Draw the MRES test: the session, one block end to end, and the staircase design.

Reads the plan from `mres_design.py` (not from the CSVs, so staircase tiers are
drawn as their real increments) and writes `../sequence/mres_test_sequence.png`
at thesis print size - 160 mm across - in the campaign's thesis style.

  (a) the session: part 1 (warm-up, not yet defined) and part 2 (stepping)
  (b) one 1/256 block, commanded position against time
  (c) the B3 two-tier staircase at every MRES, first forward + reverse pair
  (d) microstep size against the interferometer floor, with fine-tier counts
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mres_design as md  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

STYLE = HERE.parents[1] / "analysis" / "thesis.mplstyle"

# Phase hues follow the campaign's categorical order (analysis/vizstyle.py);
# the reference visits share one neutral so the measured phases stand out.
NEUTRAL = "#b4b3ad"
PHASE_COLOR = {
    "REF": NEUTRAL, "E": NEUTRAL,
    "B1": "#86857f",
    "B2": "#eda100",
    "B3": "#2a78d6",
}
PHASE_NAME = {
    "REF": "x_ref baseline / end reference",
    "B1": "B1 configure + marker",
    "B2": "B2 condition",
    "B3": "B3 two-tier staircase",
}
# MRES is a magnitude, so it gets one hue light -> dark (vizstyle SEQUENTIAL).
MRES_COLOR = dict(zip(md.MRES_TESTED,
                      ["#9ec5f4", "#6da7ec", "#3987e5", "#184f95", "#0d366b"][-len(md.MRES_TESTED):]))
MUTED = "#6b6a66"
PART1_W = 2.6                     # drawn width of the undefined part 1 [min]


def runs(rows, key=lambda r: r.phase):
    """Collapse consecutive rows with the same key into (key, t0, t1)."""
    out = []
    for r in rows:
        k = key(r)
        if out and out[-1][0] == k:
            out[-1][2] = r.t_end
        else:
            out.append([k, r.t_start, r.t_end])
    return out


def trace(rows):
    """Commanded position against time; staircase tiers as real increments."""
    t, p = [], []
    for r in rows:
        if r.kind == "stair":
            step = (r.target - r.start) / r.n_inc
            dt = r.inc_usteps * md.B3_PULSE_S + r.inc_dwell_s
            for k in range(r.n_inc):
                tk = r.t_start + k * dt
                t += [tk, tk, tk + dt]
                p += [r.start + k * step, r.start + (k + 1) * step,
                      r.start + (k + 1) * step]
        else:
            move = md.motion_duration_s(r.target - r.start, r.speed_fs) if r.speed_fs else 0
            t += [r.t_start, r.t_start + move, r.t_end]
            p += [r.start, r.target, r.target]
    return np.array(t), np.array(p, dtype=float)


def um(u256):
    return np.asarray(u256) / md.U * md.FULL_STEP_UM


def panel(ax, letter):
    ax.annotate(f"({letter})", (-0.01, 1.02), xycoords="axes fraction",
                ha="right", va="bottom", fontweight="bold")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = Path(args.out) if args.out else HERE.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    plt.style.use(str(STYLE))
    order, rows = md.build(args.seed)

    fig = plt.figure(figsize=(160 / 25.4, 158 / 25.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.5, 1.0, 1.05],
                          hspace=0.62, wspace=0.3,
                          left=0.09, right=0.975, top=0.91, bottom=0.075)

    # ------------------------------------------------------- (a) session
    ax = fig.add_subplot(gs[0, :])
    t_max = rows[-1].t_end / 60
    for y in (0,):
        # Part 1 has no rows yet: a hatched placeholder of nominal width.
        ax.add_patch(Rectangle((-PART1_W - 0.25, y - 0.32), PART1_W, 0.64,
                               fc="white", ec=NEUTRAL, lw=0.6, hatch="////"))
        ax.text(-PART1_W / 2 - 0.25, y + 0.38, "part 1", ha="center",
                va="bottom", fontsize=7, color=MUTED)
        for ph, t0, t1 in runs(rows, key=lambda r: "REF" if r.phase == "E" else r.phase):
            ax.broken_barh([(t0 / 60, (t1 - t0) / 60)], (y - 0.32, 0.64),
                           color=PHASE_COLOR[ph], lw=0)
        for seg, t0, t1 in runs(rows, key=lambda r: r.segment):
            if seg.startswith("block"):
                mres = next(r.mres for r in rows if r.segment == seg)
                ax.text((t0 + t1) / 120, y + 0.38, f"1/{mres}",
                        ha="center", va="bottom", fontsize=7.5,
                        fontweight="bold" if mres == md.ANCHOR else "normal")
                ax.plot([t0 / 60] * 2, [y - 0.4, y + 0.4], color="black", lw=0.6)
        ax.text(rows[-1].t_end / 60 + 0.3, y, f"{rows[-1].t_end / 60:.1f} min",
                va="center", fontsize=7.5, color=MUTED)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.45, 0.85)
    ax.set_xlim(-PART1_W - 0.5, t_max + 2.2)
    ax.set_xticks(np.arange(0, t_max + 1, 5))
    ax.set_xlabel("time in part 2, stepping sequence  [min]")
    handles = [Patch(fc="white", ec=NEUTRAL, hatch="////", lw=0.6,
                     label="Part 1: homing + warm-up (to be defined)")]
    handles += [Patch(color=PHASE_COLOR[k], label=PHASE_NAME[k])
                for k in ("REF", "B1", "B2", "B3")]
    ax.annotate(f"anchor 1/{md.ANCHOR} (bold) first and last", (1.0, 1.02),
                xycoords="axes fraction", ha="right", va="bottom",
                fontsize=6.5, color=MUTED)
    ax.legend(handles=handles, ncol=3, loc="lower left",
              bbox_to_anchor=(0.0, 1.16), fontsize=7, handlelength=1.2,
              columnspacing=1.2)
    panel(ax, "a")

    # ------------------------------------------------- (b) one 1/256 block
    seg = f"block{order.index(md.ANCHOR) + 1}"
    blk = [r for r in rows if r.segment == seg]
    tb0 = blk[0].t_start
    t, p = trace(blk)
    ax = fig.add_subplot(gs[1, :])
    for ph, t0, t1 in runs(blk):
        ax.axvspan(t0 - tb0, t1 - tb0, color=PHASE_COLOR[ph], alpha=0.16, lw=0)
        ax.text((t0 + t1) / 2 - tb0, 1.02, ph,
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=7.5, color=PHASE_COLOR[ph], fontweight="bold")
    ax.plot(t - tb0, um(p), color="black", lw=0.6)
    # +/-2 mm conditioning and a 40 um staircase on one axis: symlog keeps the
    # staircase linear and compresses the strokes and the marker.
    ax.set_yscale("symlog", linthresh=50, linscale=1.6)
    ax.set_yticks([-1000, -100, 0, 20, 40, 100, 1000],
                  ["$-10^3$", "$-10^2$", "0", "20", "40", "$10^2$", "$10^3$"])
    ax.minorticks_off()
    ax.set_xscale("function", functions=(lambda x: np.sqrt(np.maximum(x, 0)),
                                         lambda x: x**2))
    ax.set_xticks([0, 5, 20, 40, 100, 200, 300, 400])
    ax.set_xlim(0, blk[-1].t_end - tb0)
    ax.set_xlabel(r"time within block  [s, $\sqrt{\ }$ scale]")
    ax.set_ylabel(r"commanded position  [$\mu$m]")
    ax.annotate(f"1/{md.ANCHOR}, block 1 - {(blk[-1].t_end - tb0) / 60:.1f} min",
                (0.995, 0.04), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=7.5, color=MUTED)
    panel(ax, "b")

    # -------------------------------------------- (c) B3 two-tier staircase
    ax = fig.add_subplot(gs[2, 0])
    by_mres = {}
    for r in rows:
        if r.phase == "B3" and r.kind == "stair" and r.repeat == 0:
            by_mres.setdefault(r.mres, []).append(r)
    fine_um = um(md.fine_span_u256())
    span_um = um(md.B3_SPAN)
    tmax = 0
    for mres in md.MRES_TESTED:
        legs = by_mres[mres][:4]
        t, p = trace(legs)
        t -= legs[0].t_start
        tmax = max(tmax, t[-1])
        ax.plot(t, um(p - legs[0].start), color=MRES_COLOR[mres], lw=0.7,
                label=f"1/{mres}")
    ax.axhspan(0, fine_um, color=PHASE_COLOR["B3"], alpha=0.10, lw=0)
    ax.axhspan(span_um - fine_um, span_um, color=PHASE_COLOR["B3"], alpha=0.10, lw=0)
    ax.text(0.7 * tmax, fine_um / 2, "shaded: B3a fine tiers", ha="center",
            va="center", fontsize=6.5, color=MUTED)
    ax.set_xlim(0, tmax)
    ax.set_ylim(-1.5, span_um + 1.5)
    ax.set_xlabel("time  [s]")
    ax.set_ylabel(r"commanded position  [$\mu$m]")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=6.5,
              ncol=5, handlelength=0.9, columnspacing=0.7, handletextpad=0.4,
              borderaxespad=0.2)
    panel(ax, "c")

    # ------------------------------------ (d) resolvability and fine tier
    ax = fig.add_subplot(gs[2, 1])
    xs = np.arange(len(md.MRES_TESTED))
    deltas = [md.delta_nm(m) for m in md.MRES_TESTED]
    fine = md.fine_span_u256()
    ax.plot(xs, deltas, color="black", lw=0.6, zorder=1)
    for x, m, d in zip(xs, md.MRES_TESTED, deltas):
        ax.plot(x, d, "o", ms=4, color=MRES_COLOR[m], mec="black", mew=0.4,
                zorder=2)
        ax.annotate(f"{fine * m // md.U}", (x, d), xytext=(0, 5),
                    textcoords="offset points", ha="center", fontsize=6.5)
    floor = 3 * md.SIGMA_FLOOR_NM
    ax.axhline(floor, color="black", lw=0.6, ls=(0, (3, 2)))
    ax.text(xs[0] - 0.2, floor * 1.25,
            rf"$3\sigma_\mathrm{{floor}}$ = {floor:g} nm (placeholder)",
            fontsize=6.5, color=MUTED, va="bottom")
    ax.set_yscale("log")
    ax.set_xticks(xs, [f"1/{m}" for m in md.MRES_TESTED])
    ax.set_xlim(-0.4, len(xs) - 0.6)
    ax.set_ylim(floor / 2, max(deltas) * 3)
    ax.set_xlabel("microstep resolution")
    ax.set_ylabel(r"microstep size $\delta$  [nm]")
    ax.annotate(f"labels: B3a microsteps over $S_\\mathrm{{fine}}$ = "
                f"{fine_um:.3f} $\\mu$m",
                (0.98, 0.96), xycoords="axes fraction", ha="right", va="top",
                fontsize=6.5, color=MUTED)
    panel(ax, "d")

    png = out / "mres_test_sequence.png"
    fig.savefig(png, dpi=600)
    fig.savefig(out / "mres_test_sequence.pdf")
    plt.close(fig)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
