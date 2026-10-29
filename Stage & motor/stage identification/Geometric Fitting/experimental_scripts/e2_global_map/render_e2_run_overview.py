#!/usr/bin/env python3
"""The E2 run as it was actually measured, end to end.

The capture opens part-way through the homing sweep that preceded the
sequence, so the first six minutes run the full travel rather than the plan.
That is left in: it is what the raw file contains, and it is the only part of
the run that approaches the switches.

State bands come from the plan, placed at the times of the measured holds they
were matched to. This scope-derived alignment is intentionally independent of
the available controller log.

Writes rendered/e2_run_overview.png.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e2_run
from e2_run import X_REF_MM, TURN_LOW_MM, TURN_HIGH_MM

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parents[3] / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, style   # noqa: E402

OUT = Path(os.environ.get("E2_OUT", _HERE / "rendered"))
MINUS = "−"
COLOURS = {"HOMING": "#68717a", "DATUM": "#7a6f91", "UP": "#4c78a8",
           "TURN_REF": "#817765", "DOWN": "#598a76", "END_REF": "#817765"}
V_CMD = 2.0


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


def blocks(run):
    """(start, end, state) in minutes, from the matched holds."""
    ok = np.flatnonzero(run.matched)
    state = run.plan.state[ok]
    start = run.t[run.hold_start[ok]] / 60.0
    end = run.t[np.minimum(run.hold_end[ok], run.t.size - 1)] / 60.0
    out = [(0.0, start[0], "HOMING")]
    a, cur = start[0], state[0]
    for i in range(1, state.size):
        if state[i] != cur:
            out.append((a, start[i], cur))
            a, cur = start[i], state[i]
    out.append((a, end[-1], cur))
    return out


def main():
    thesis_style()
    run = e2_run.load()
    tm = run.t / 60.0
    spans = blocks(run)
    ok = run.matched

    fig = plt.figure(figsize=(11.0, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=(0.42, 2.3, 1.0), hspace=0.52)
    ax_s = fig.add_subplot(gs[0])
    ax_p = fig.add_subplot(gs[1], sharex=ax_s)
    ax_v = fig.add_subplot(gs[2], sharex=ax_s)

    ax_s.grid(False)
    for a, b, st in spans:
        ax_s.axvspan(a, b, color=COLOURS.get(st, "#999999"), lw=0)
        if b - a > 2.0:
            ax_s.text(0.5 * (a + b), 0.5, st, ha="center", va="center",
                      color="#ffffff", fontsize=9)
    for n, (a, b, st) in enumerate(s for s in spans if s[1] - s[0] <= 2.0):
        ax_s.annotate(st, (0.5 * (a + b), 0.0),
                      xycoords=("data", "axes fraction"),
                      textcoords="offset points",
                      xytext=(0, -8 if n % 2 == 0 else -19), fontsize=7.5,
                      ha="right" if a > 0.93 * spans[-1][1] else "center",
                      va="top", color=INK2,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=INK2,
                                      shrinkA=0, shrinkB=1))
    ax_s.set_ylim(0, 1)
    ax_s.set_yticks([])
    ax_s.tick_params(labelbottom=False, bottom=False)
    ax_s.set_title("(a) phase of the run, from the plan matched to the holds",
                   fontsize=11, loc="left", pad=8)

    for a, b, st in spans:
        k = slice(*np.searchsorted(tm, (a, b)))
        ax_p.plot(tm[k], run.pos[k], lw=1.0, color=COLOURS.get(st, "#999999"),
                  rasterized=True)
    for y, label in ((TURN_LOW_MM, f"lower turn {TURN_LOW_MM:.0f} mm"),
                     (TURN_HIGH_MM, f"upper turn {TURN_HIGH_MM:.3f} mm")):
        ax_p.axhline(y, color="#99483e", lw=0.8, ls="--", zorder=2)
    ax_p.axhline(X_REF_MM, color=INK2, lw=0.8, ls=":", zorder=2)
    ax_p.text(tm[-1], TURN_HIGH_MM + 2.5, f"upper turn {TURN_HIGH_MM:.3f} mm",
              ha="right", fontsize=8.5, color="#99483e")
    ax_p.text(tm[-1], TURN_LOW_MM - 5.5, f"lower turn {TURN_LOW_MM:.1f} mm",
              ha="right", fontsize=8.5, color="#99483e")
    ax_p.text(tm[0] + 0.3, X_REF_MM + 2.0,
              r"$x_{\mathrm{ref}}$ = 61.35 mm", fontsize=8.5, color=INK2)
    lo = run.pos[run.hold_start[ok].min():run.hold_end[ok].max()]
    ax_p.set_ylabel("stage position [mm]")
    ax_p.set_ylim(-6, 130)
    ax_p.tick_params(labelbottom=False)
    ax_p.set_title(f"(b) measured position, 1 kHz "
                   f"({run.t.size:,} samples) — measurement stays within "
                   f"{lo.min():.3f}–{lo.max():.3f} mm", fontsize=11,
                   loc="left")

    v = np.abs(np.diff(run.pos)) * 1000.0
    v10 = np.convolve(v, np.ones(10) / 10, mode="valid")
    idx = np.flatnonzero(ok)
    pk, when = [], []
    for a, b in zip(run.hold_end[idx[:-1]], run.hold_start[idx[1:]]):
        if b > a + 20:
            pk.append(v10[a:b].max())
            when.append(run.t[b] / 60.0)
    pk, when = np.array(pk), np.array(when)
    for st, colour in (("UP", COLOURS["UP"]), ("DOWN", COLOURS["DOWN"])):
        lo_t, hi_t = [(a, b) for a, b, s in spans if s == st][0]
        m = (when >= lo_t) & (when <= hi_t)
        ax_v.plot(when[m], pk[m], ".", ms=3.0, color=colour, rasterized=True)
    other = np.ones(pk.size, bool)
    for st in ("UP", "DOWN"):
        lo_t, hi_t = [(a, b) for a, b, s in spans if s == st][0]
        other &= ~((when >= lo_t) & (when <= hi_t))
    ax_v.plot(when[other], pk[other], ".", ms=3.0, color=COLOURS["DATUM"],
              rasterized=True)
    ax_v.axhline(V_CMD, color="#99483e", lw=0.9, ls="--", zorder=3)
    ax_v.text(tm[0] + 0.3, V_CMD - 0.04, "commanded 2.0 mm/s", fontsize=8.5,
              color="#99483e", va="top")
    ax_v.set_ylim(1.5, 2.35)
    ax_v.set_xlabel("time in the capture [min]")
    ax_v.set_ylabel("peak speed [mm/s]")
    ax_v.set_title(f"(c) peak speed of each move between dwells — median "
                   f"{np.median(pk):.2f} mm/s; AquadB ceiling 8.333 mm/s",
                   fontsize=11, loc="left")
    ax_v.set_xlim(0, tm[-1])

    fig.text(0.012, 0.012,
             "The datum "
             f"({run.datum_mm:.4f} mm) and every window edge come from the "
             "trace, with the plan supplying stop order; the serial log is "
             "cross-checked separately in the findings analysis.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "e2_run_overview.png", dpi=200)
    print("wrote", OUT / "e2_run_overview.png")


if __name__ == "__main__":
    main()
