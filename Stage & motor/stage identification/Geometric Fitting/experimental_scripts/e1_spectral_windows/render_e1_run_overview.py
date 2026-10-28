#!/usr/bin/env python3
"""The E1 run as it was actually measured, end to end.

One figure over the whole 1 kHz capture: which state the firmware was in, where
the stage was, and how fast it moved. The capture starts inside the homing move
that preceded the measurement, so the first minutes are the approach to x_ref
rather than part of the sequence - that is left in deliberately, since it is
what the raw file contains.

Writes rendered/e1_run_overview.{png,pdf}.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e1_run
from e1_run import WINDOWS, OVERRUN_MM, X_REF_MM

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parents[3] / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, style   # noqa: E402

OUT = Path(os.environ.get("E1_OUT", _HERE / "rendered"))
MINUS = "−"
COLOURS = {
    "HOME": "#68717a", "INIT": "#68717a", "DATUM": "#7a6f91",
    "PILOT": "#a9783f", "W1": "#4c78a8", "W2": "#598a76", "W3": "#9d6770",
    "REF": "#817765", "END_REF": "#817765", "PARK": "#68717a",
    "DONE": "#68717a",
}
V_CMD = 2.0


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


def state_spans(run):
    """Contiguous (start, end, state) in minutes, clipped to the capture."""
    m = run.events.kind == "state_enter"
    ms, st = run.events.ms[m] / 1000.0, run.events.state[m]
    keep = np.concatenate(([True], st[1:] != st[:-1]))
    ms, st = ms[keep], st[keep]
    lo, hi = run.t[0], run.t[-1]
    spans = []
    for i, s in enumerate(st):
        a = ms[i]
        b = ms[i + 1] if i + 1 < ms.size else hi
        a, b = max(a, lo), min(b, hi)
        if b > a:
            spans.append((a / 60.0, b / 60.0, s))
    return spans


def main():
    thesis_style()
    run = e1_run.load()
    tm = run.t / 60.0
    spans = state_spans(run)

    fig = plt.figure(figsize=(11.0, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=(0.42, 2.3, 1.0), hspace=0.52)
    ax_s = fig.add_subplot(gs[0])
    ax_p = fig.add_subplot(gs[1], sharex=ax_s)
    ax_v = fig.add_subplot(gs[2], sharex=ax_s)

    # (a) what the firmware was doing
    ax_s.grid(False)
    for a, b, st in spans:
        ax_s.axvspan(a, b, color=COLOURS.get(st, "#999999"), lw=0)
        if b - a > 1.2:
            ax_s.text(0.5 * (a + b), 0.5, st, ha="center", va="center",
                      color="#ffffff", fontsize=9)
    short = [(a, b, st) for a, b, st in spans if b - a <= 1.2]
    for n, (a, b, st) in enumerate(short):
        if st in ("INIT", "HOME", "DONE"):
            continue
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
    ax_s.set_title("(a) firmware state, from the ESP log", fontsize=11,
                   loc="left", pad=8)

    # (b) the measured trajectory
    for name, a, b in WINDOWS:
        ax_p.axhspan(a, b, color=COLOURS[name], alpha=0.10, lw=0)
        ax_p.text(tm[-1], 0.5 * (a + b), f" {name}", va="center", fontsize=8.5,
                  color=COLOURS[name])
    ax_p.axhline(X_REF_MM, color=INK2, lw=0.8, ls="--", zorder=2)
    ax_p.text(tm[0] + 0.3, X_REF_MM + 1.2, r"$x_{\mathrm{ref}}$ = 61.35 mm",
              fontsize=8.5, color=INK2)
    for a, b, s in spans:
        k = slice(*np.searchsorted(tm, (a, b)))
        ax_p.plot(tm[k], run.pos[k], lw=1.0, color=COLOURS.get(s, "#999999"),
                  rasterized=True, solid_joinstyle="round")

    for name, a, b in WINDOWS:
        t_pre = run.events.when(name, "arrive", "lower_pre")[0] / 60.0
        t_ovr = run.events.when(name, "arrive", "upper_overrun")[0] / 60.0
        for tt, y in ((t_pre, a - OVERRUN_MM), (t_ovr, b + OVERRUN_MM)):
            ax_p.plot([tt], [y], marker="v" if y < a else "^", ms=6,
                      mfc="#ffffff", mec=COLOURS[name], mew=1.2, zorder=5)
    ax_p.plot([], [], "v", ms=6, mfc="#ffffff", mec=INK2, mew=1.2,
              label=f"{OVERRUN_MM:.0f} mm pre-approach below the window")
    ax_p.plot([], [], "^", ms=6, mfc="#ffffff", mec=INK2, mew=1.2,
              label=f"{OVERRUN_MM:.0f} mm overrun above the window")
    ax_p.legend(loc="center left", fontsize=8.5, framealpha=1.0,
                bbox_to_anchor=(0.012, 0.33))
    ax_p.set_ylabel("stage position [mm]")
    ax_p.set_ylim(-4, 128)
    ax_p.tick_params(labelbottom=False)
    ax_p.set_title("(b) measured position, 1 kHz interferometer "
                   f"({run.t.size:,} samples, no gaps)",
                   fontsize=11, loc="left")

    # (c) speed, one point per commanded move
    v = np.diff(run.pos) * 1000.0
    k = 10
    v10 = np.convolve(np.abs(v), np.ones(k) / k, mode="valid")
    tv = run.t[k // 2: k // 2 + v10.size]

    def peak(a, b):
        i, j = np.searchsorted(tv, (a, b))
        return v10[i:j].max() if j > i else np.nan

    st = run.stops
    grid = np.array([peak(a, b) for a, b in zip(st.move_start_s, st.arrive_s)])
    for name, _, _ in WINDOWS:
        m = st.window == name
        ax_v.plot(st.arrive_s[m] / 60.0, grid[m], ".", ms=3.0,
                  color=COLOURS[name], rasterized=True)
    m = np.isin(st.window, ("PILOT", "REF"))
    ax_v.plot(st.arrive_s[m] / 60.0, grid[m], ".", ms=3.0, color=COLOURS["PILOT"],
              rasterized=True)
    long_moves = [(name, run.events.when(name, "move_start", lab)[0],
                   run.events.when(name, "arrive", lab)[0])
                  for name, _, _ in WINDOWS
                  for lab in ("lower_pre", "upper_overrun")]
    ax_v.plot([0.5 * (a + b) / 60.0 for _, a, b in long_moves],
              [peak(a, b) for _, a, b in long_moves], "s", ms=5,
              mfc="#ffffff", mec=INK, mew=1.1, zorder=5,
              label="pre-approach and overrun moves")
    ax_v.axhline(V_CMD, color="#99483e", lw=0.9, ls="--", zorder=3)
    ax_v.text(tm[0] + 0.3, V_CMD - 0.09, "commanded 2.0 mm/s", fontsize=8.5,
              color="#99483e", va="top")
    ax_v.plot([], [], ".", ms=8, color=INK2, label="grid and reference moves")
    ax_v.legend(loc="lower right", fontsize=8.5, framealpha=1.0, ncol=2)
    ax_v.set_ylim(1.5, 2.45)
    ax_v.set_xlabel("time on the ESP clock [min]")
    ax_v.set_ylabel("peak speed [mm/s]")
    ax_v.set_title("(c) peak speed of every commanded move — all "
                   f"{np.isfinite(grid).sum()} reach the commanded 2.0 mm/s "
                   f"(median {np.nanmedian(grid):.2f}); ceiling 8.333 mm/s",
                   fontsize=11, loc="left")
    ax_v.set_xlim(tm[0], tm[-1])

    fig.text(0.012, 0.012,
             "Capture put on the ESP clock by fitting its rate "
             f"({run.rate_ppm:+.1f} ppm) and offset to all 676 logged moves; "
             .replace("-", MINUS) +
             f"position offset from the 654 holds, fit MAD "
             f"{run.fit_mad_um:.2f} µm.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    OUT.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"e1_run_overview.{ext}", dpi=200)
    print("wrote", OUT / "e1_run_overview.png")


if __name__ == "__main__":
    main()
