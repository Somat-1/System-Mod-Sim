#!/usr/bin/env python3
"""Where the v3 record windows actually landed, and what it cost.

The ISO evaluation took each 0.5 s averaging window from the board log and
mapped it onto the interferometer trace with a fitted clock. The clock was
fitted against a stop detector that smoothed velocity with a centred
51-sample boxcar, which reports rest about 34 ms late, so the whole board
timeline came out ~38 ms late. The firmware starts the next move at the
instant the window closes, so there was no margin to absorb that: every
window ran past rest and into the following move.

Four panels, following the chain:

    a   one hold, zoomed: the window closes 38 ms after the stage has left
    b   all 240 windows, aligned on their closing edge - it is systematic,
        and it is signed with the pass direction
    c   what that does to each window mean: up biased one way, down the
        other, so B = xbar_up - xbar_dn adds the two
    d   mean B_i against the assumed clock offset: a flat plateau at the
        true -1.60 um, and the cliff the published -0.365 um sits on

Writes plots/temp/v3_window_overrun.{png,pdf}.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
for p in (str(_MOTOR), str(_HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)
from staircase_figures import INK, INK2, INK3, style, titleblock    # noqa: E402
from iso230_2_evaluate import (load_stops, load_scope, fit_clock,   # noqa: E402
                               evaluate, summarise)

OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))
CACHE = Path(os.environ.get("ISO_CACHE", OUT))

UP, DOWN = "#5f8bc5", "#c88e5c"      # categorical: pass direction
BAD, GOOD = "#99483e", "#356449"     # the defective window, the clear one

MOVE_GATE = 0.10        # um/ms; the traverse runs at ~2 um/ms, rest noise <0.03
MOVE_RUN = 20           # consecutive samples required, so noise cannot trip it
CLEAR_BACK = 50         # a clear window ends this far before the stage leaves

MU = "µm"


def num(v, dp=3, sign=True):
    """A value label with a typographic minus, as the other figures use."""
    fmt = "%+." + str(dp) + "f" if sign else "%." + str(dp) + "f"
    return (fmt % v).replace("-", "−")


def load(log, csv):
    """Stops from the log, trace from the CSV, with the trace cached."""
    stops, _ = load_stops(log)
    cache = CACHE / "v3_pos_mm.npy"
    if cache.exists():
        pos = np.load(cache)
    else:
        pos = load_scope(csv)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache, pos)
    return stops, pos


def onsets(pos, a_, b_):
    """For each window, the sample at which the next move really begins."""
    sp = np.convolve(np.abs(np.diff(pos)), np.ones(5) / 5, mode="same") * 1000
    out = []
    for a, b in zip(a_, b_):
        seg = (sp[a:b + 4000] > MOVE_GATE).astype(int)
        c = np.convolve(seg, np.ones(MOVE_RUN, int), mode="valid")
        out.append(a + int(np.argmax(c == MOVE_RUN)) if (c == MOVE_RUN).any()
                   else -1)
    return np.array(out)


def bracket(ax, x0, x1, y, colour, label, below=False):
    """A span marker with end caps, for drawing a window on the time axis."""
    h = 0.75
    ax.plot([x0, x1], [y, y], color=colour, lw=1.5, zorder=7)
    for x in (x0, x1):
        ax.plot([x, x], [y - h, y + h], color=colour, lw=1.5, zorder=7)
    ax.text(0.5 * (x0 + x1), y - 1.4 if below else y + 1.4, label,
            color=colour, size=8.6, ha="center",
            va="top" if below else "bottom", zorder=7)


def panel_a(ax, pos, a, b, onset, direction):
    """One hold: the record window, and the move that starts inside it."""
    t0, t1 = a - 110, b + 12
    t = np.arange(t0, t1) - b
    rest = pos[onset - 550:onset - CLEAR_BACK].mean()
    y = (pos[t0:t1] - rest) * 1000.0

    ax.axvspan(onset - b, 0, color=BAD, alpha=0.20, lw=0, zorder=1)
    ax.axvline(onset - b, color=BAD, lw=1.2, zorder=4)
    ax.axhline(0.0, color=INK3, lw=0.9, zorder=2)
    ax.plot(t, y, color=INK, lw=1.5, zorder=5, solid_capstyle="round")

    bracket(ax, a - b, 0, -3.4, BAD,
            "the 0.5 s window the evaluation averaged")
    bracket(ax, onset - 550 - b, onset - CLEAR_BACK - b, -7.2, GOOD,
            "the same 0.5 s, clear of the move", below=True)

    m_bad = (pos[a:b].mean() - rest) * 1000.0
    ax.annotate("the next move starts here,\n%d ms before the window closes"
                % (b - onset),
                xy=(onset - b, 2.2), xytext=(onset - b - 55, 11.0),
                color=BAD, size=8.8, ha="right", va="center",
                arrowprops=dict(arrowstyle="-", color=BAD, lw=0.9,
                                shrinkA=4, shrinkB=4))
    ax.text(-58, 24.6, "%s %s of travel\nfalls inside the window,\nand its "
            "mean comes out\n%s %s off the rest level"
            % (num((pos[b] - pos[onset]) * 1000.0, 1), MU, num(m_bad), MU),
            color=BAD, size=8.8, ha="right", va="top", zorder=7)

    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-9.5, 26.0)
    ax.set_xlabel("time relative to the window closing (ms)")
    ax.set_ylabel("position, from the settled level (%s)" % MU)
    ax.set_title("a   one hold, as the published evaluation saw it (%s pass)"
                 % ("up" if direction > 0 else "down"),
                 loc="left", color=INK, size=11, pad=10)


def panel_b(ax, pos, b_, onset, dir_):
    """Every window, aligned on its closing edge, levelled on its own rest."""
    span, lead = 560, 70
    t = np.arange(-span, lead)
    med = {}
    for d, colour in ((+1, UP), (-1, DOWN)):
        sel = np.where(dir_ == d)[0]
        stack = np.empty((sel.size, t.size))
        for r, k in enumerate(sel):
            b, o = b_[k], onset[k]
            rest = pos[o - 550:o - CLEAR_BACK].mean()
            stack[r] = (pos[b - span:b + lead] - rest) * 1000.0
        ax.plot(t, stack.T, color=colour, lw=0.5, alpha=0.18, zorder=2)
        med[d] = np.median(stack, axis=0)

    over = b_ - onset
    ax.axvspan(-span, 0, color=BAD, alpha=0.055, lw=0, zorder=0)
    ax.axvspan(-np.median(over), 0, color=BAD, alpha=0.16, lw=0, zorder=1)
    ax.axvline(0, color=INK3, lw=1.0, zorder=3)

    ax.plot(t, med[+1], color=UP, lw=2.0, zorder=5,
            label="up pass, median of %d" % int((dir_ == 1).sum()))
    ax.plot(t, med[-1], color=DOWN, lw=2.0, zorder=5,
            label="down pass, median of %d" % int((dir_ == -1).sum()))

    ax.text(-95, 25.5,
            "every one of the %d windows closes\n%.0f ms into the next move "
            "(%d–%d ms)"
            % (over.size, np.median(over), over.min(), over.max()),
            color=BAD, size=8.8, ha="right", va="top")
    ax.text(-span + 16, 4.0, "flat while the stage is at rest\n"
            "— the first 460 ms of every window are sound",
            color=INK2, size=8.8, ha="left", va="bottom")

    ax.set_xlim(-span, lead)
    ax.set_ylim(-30, 30)
    ax.set_xlabel("time relative to the window closing (ms)")
    ax.set_ylabel("position, from the settled level (%s)" % MU)
    ax.set_title("b   all 240 windows, levelled on their own rest",
                 loc="left", color=INK, size=11, pad=10)
    ax.legend(loc="lower left", ncol=1, handlelength=2.0)


def panel_c(ax, bias, dir_):
    """The bias each window mean carries, by direction, and what B does."""
    rng = np.random.default_rng(7)
    stats = {}
    for row, (d, colour, label) in enumerate(((+1, UP, "up pass"),
                                              (-1, DOWN, "down pass"))):
        v = bias[dir_ == d]
        y = 1 - row + rng.uniform(-0.17, 0.17, v.size)
        ax.scatter(v, y, s=13, color=colour, alpha=0.55, lw=0, zorder=3)
        m = float(np.median(v))
        stats[d] = m
        ax.plot([m, m], [1 - row - 0.30, 1 - row + 0.30], color=colour,
                lw=2.4, zorder=5)
        ax.text(m, 1 - row + 0.40, "%s   median %s %s" % (label, num(m), MU),
                color=colour, size=9.0, va="bottom",
                ha="right" if m > 0 else "left")

    ax.axvline(0.0, color=INK3, lw=1.0, zorder=1)
    shift = stats[+1] - stats[-1]
    ax.annotate("", xy=(stats[+1], -0.50), xytext=(stats[-1], -0.50),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2,
                                shrinkA=0, shrinkB=0), zorder=6)
    ax.text(0.0, -0.66,
            "B = x̄↑ − x̄↓ subtracts them, so the two "
            "biases add:   %s %s" % (num(shift), MU),
            color=INK, size=9.4, ha="center", va="top")
    ax.text(0.0, -1.06,
            "(x̄↑ + x̄↓)/2 sums them, so they cancel "
            "— which is why E, M\nand the gain all came out looking "
            "normal", color=INK2, size=8.8, ha="center", va="top")

    ax.set_yticks([])
    ax.set_xlim(-0.95, 0.95)
    ax.set_ylim(-1.95, 1.75)
    ax.set_xlabel("window mean minus the settled level (%s)" % MU)
    ax.set_title("c   the bias on each of the 240 window means",
                 loc="left", color=INK, size=11, pad=10)
    ax.grid(axis="y", visible=False)


def panel_d(ax, stops, pos, rate, off, grid):
    """Mean B_i as a function of the assumed clock offset."""
    vals = []
    for g in grid:
        per, *_ = evaluate(stops, pos, rate, off + g)
        vals.append(summarise(per)["B_mean"])
    vals = np.array(vals)

    ax.plot(grid, vals, color=INK, lw=1.6, zorder=4)
    ax.scatter(grid, vals, s=11, color=INK, lw=0, zorder=5)

    plateau = float(np.median(vals[grid <= -60]))
    ax.axhline(plateau, color=GOOD, lw=1.4, ls=(0, (4, 2)), zorder=2,
               label="true reversal, %s %s (plateau, and the blind result)"
                     % (num(plateau), MU))
    pub = float(vals[np.argmin(np.abs(grid))])
    ax.scatter([0], [pub], s=78, facecolor="none", edgecolor=BAD, lw=1.8,
               zorder=6, label="as published, %s %s" % (num(pub), MU))

    ax.axvspan(grid[0], -45, color=GOOD, alpha=0.09, lw=0, zorder=0)
    ax.text(grid[0] + 5, plateau + 0.55,
            "window clear of the move:\nflat, insensitive to the offset",
            color=GOOD, size=8.8, ha="left", va="bottom")
    ax.annotate("+25 ms of offset is +1.2 %s of\nreversal — and "
                "fit_clock searched\non a 25 ms grid" % MU,
                xy=(0, pub), xytext=(-30, 4.3), color=BAD, size=8.8,
                ha="right", va="top",
                arrowprops=dict(arrowstyle="-", color=BAD, lw=0.9,
                                shrinkA=2, shrinkB=5))

    ax.set_xlim(grid[0], grid[-1])
    ax.set_ylim(-2.6, 8.6)
    ax.set_xlabel("shift applied to the fitted clock offset (ms)")
    ax.set_ylabel("mean Bᵢ over the 24 targets (%s)" % MU)
    ax.set_title("d   the cliff the published value sits on",
                 loc="left", color=INK, size=11, pad=10)
    ax.legend(loc="upper left", ncol=1, handlelength=2.0)


def main() -> int:
    log = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        _HERE.parent / "data" / "terminal_20261021_183115.log")
    csv = Path(sys.argv[2]) if len(sys.argv) > 2 else (
        _HERE.parents[1] / "ISO230-2_v3.csv")

    stops, pos = load(log, csv)
    rate, off, resid = fit_clock(stops, pos)
    iso = [s for s in stops if s["state"] == "ISO"]
    a_ = np.array([int(round(s["t0"] * rate + off)) for s in iso])
    b_ = np.array([int(round(s["t1"] * rate + off)) for s in iso])
    dir_ = np.array([s["dir"] for s in iso])
    onset = onsets(pos, a_, b_)
    over = b_ - onset

    clear = np.array([pos[o - 550:o - CLEAR_BACK].mean() for o in onset])
    bias = (np.array([pos[a:b].mean() for a, b in zip(a_, b_)]) - clear) * 1000
    travel = np.abs(pos[b_] - pos[onset]) * 1000.0
    shift = np.median(bias[dir_ == 1]) - np.median(bias[dir_ == -1])

    print("clock: rate %.6f, offset %.0f ms, resid %.1f ms" % (rate, off, resid))
    print("windows: %d, all cutting into the next move; overrun median %.0f ms "
          "(%d-%d ms)" % (over.size, np.median(over), over.min(), over.max()))
    print("bias: up %+.3f, down %+.3f um; B shift %+.3f um"
          % (np.median(bias[dir_ == 1]), np.median(bias[dir_ == -1]), shift))

    cand = np.where((dir_ == +1) & (over == np.median(over)))[0]
    k = int(cand[cand.size // 2])

    style()
    fig = plt.figure(figsize=(13.2, 10.4))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.22,
                  top=0.805, bottom=0.058, left=0.058, right=0.985)
    panel_a(fig.add_subplot(gs[0, 0]), pos, a_[k], b_[k], onset[k], dir_[k])
    panel_b(fig.add_subplot(gs[0, 1]), pos, b_, onset, dir_)
    panel_c(fig.add_subplot(gs[1, 0]), bias, dir_)
    panel_d(fig.add_subplot(gs[1, 1]), stops, pos, rate, off,
            np.arange(-120, 41, 5.0))

    titleblock(
        fig,
        "ISO 230-2 v3: the averaging windows ran into the next move",
        "Record windows taken from the board log through the fitted clock, "
        "drawn on the interferometer trace. The capture is sound; the window "
        "placement is not.\n"
        "%s  ·  %s  ·  240 ISO stops, 24 targets × 2 directions "
        "× 5 runs  ·  motion gate 0.10 %s/ms held 20 ms"
        % (Path(csv).name, Path(log).name, MU),
        "Every window closes %.0f ms after the stage has left the target "
        "— %.0f %s of travel, signed with the pass — which shifts B "
        "by %+.2f %s and turns −1.60 %s into −0.37 %s."
        % (np.median(over), np.median(travel), MU, shift, MU, MU, MU))

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        p = OUT / ("v3_window_overrun.%s" % ext)
        fig.savefig(p)
        print("wrote %s" % p)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
