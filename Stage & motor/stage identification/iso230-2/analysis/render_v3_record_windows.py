#!/usr/bin/env python3
"""The v3 record windows drawn on the trace, over three steps of the sequence.

Each shaded band is the half-second the published evaluation averaged to get
one measured position - the window edges come from the board log, mapped onto
the interferometer trace by fit_clock() exactly as iso230_2_evaluate.py does
it. The point of the figure is to show which samples that actually selected.

The inset is the closing edge of the middle window at full sample resolution:
the window does not end at rest, it ends inside the following move.

Writes plots/temp/v3_record_windows.{png,pdf}.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
for p in (str(_MOTOR), str(_HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)
from staircase_figures import INK, INK2, INK3, style, titleblock   # noqa: E402
from iso230_2_evaluate import (load_stops, load_scope,             # noqa: E402
                               fit_clock)

OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))
CACHE = Path(os.environ.get("ISO_CACHE", OUT))

WINDOW = "#99483e"      # the half-second that was averaged
PLATEAU = "#356449"     # the stage actually at rest
MU = "µm"

MOVE_GATE = 0.10        # um/ms; the traverse runs at ~2 um/ms
MOVE_RUN = 20           # consecutive samples, so noise cannot trip the gate
NSTEP = 3               # how many steps to show


def load(log, csv):
    stops, _ = load_stops(log)
    cache = CACHE / "v3_pos_mm.npy"
    if cache.exists():
        pos = np.load(cache)
    else:
        pos = load_scope(csv)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache, pos)
    return stops, pos


def moving_mask(pos):
    """True where the stage is really moving.

    The gate has to be held for MOVE_RUN samples before it counts, so a stray
    noise sample at rest cannot open it - but the mask must still start at the
    first sample of that run, not MOVE_RUN/2 before it. So: mark the starts of
    every full run, then spread each start forward over its own run. Smoothing
    the gate instead would drag the edge earlier and overstate the overrun.
    """
    sp = np.convolve(np.abs(np.diff(pos)), np.ones(5) / 5, mode="same") * 1000
    raw = (sp > MOVE_GATE).astype(int)
    full = (np.convolve(raw, np.ones(MOVE_RUN, int), "valid") == MOVE_RUN)
    held = np.convolve(full.astype(int), np.ones(MOVE_RUN, int),
                       "full")[:raw.size] > 0
    return np.append(held, held[-1])


def rest_spans(moving, lo, hi, min_ms=800):
    """Contiguous at-rest intervals inside [lo, hi)."""
    m = moving[lo:hi]
    edges = np.flatnonzero(np.diff(np.concatenate(([1], m.astype(int), [1]))))
    out = []
    for s, e in zip(edges[0::2], edges[1::2]):
        if e - s >= min_ms:
            out.append((lo + s, lo + e))
    return out


def main() -> int:
    log = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        _HERE.parent / "data" / "terminal_20261021_183115.log")
    csv = Path(sys.argv[2]) if len(sys.argv) > 2 else (
        _HERE.parents[1] / "ISO230-2_v3.csv")

    stops, pos = load(log, csv)
    rate, off, _ = fit_clock(stops, pos)
    iso = [s for s in stops if s["state"] == "ISO"]
    a_ = np.array([int(round(s["t0"] * rate + off)) for s in iso])
    b_ = np.array([int(round(s["t1"] * rate + off)) for s in iso])
    dir_ = np.array([s["dir"] for s in iso])
    moving = moving_mask(pos)

    # the one constant that ties the counter to machine coordinates, exactly
    # as evaluate() derives it: no gain, no per-target offset
    meas = np.array([pos[a:b].mean() for a, b in zip(a_, b_)])
    cmd = np.array([s["cmd_mm"] for s in iso])
    pos = pos - float(np.mean(meas - cmd))

    onset_all = np.array([a + int(np.argmax(moving[a:b + 3000]))
                          for a, b in zip(a_, b_)])
    over_med = float(np.median(b_ - onset_all))

    # three consecutive stops of one ascending pass, from the middle of the run
    k0 = next(k for k in range(len(iso) // 2, len(iso) - NSTEP)
              if set(dir_[k:k + NSTEP]) == {+1})
    ks = range(k0, k0 + NSTEP)
    lo, hi = a_[k0] - 1800, b_[k0 + NSTEP - 1] + 2600
    t = (np.arange(lo, hi) - lo) / 1000.0
    y = pos[lo:hi]

    style()
    fig = plt.figure(figsize=(12.6, 6.4))
    ax = fig.add_axes([0.072, 0.108, 0.910, 0.675])
    ax.set_axisbelow(True)

    for n, (s, e) in enumerate(rest_spans(moving, lo, hi)):
        ax.axvspan((s - lo) / 1000.0, (e - lo) / 1000.0, color=PLATEAU,
                   alpha=0.10, lw=0, zorder=0,
                   label="the stage at rest on a target" if n == 0 else None)
    for n, k in enumerate(ks):
        ax.axvspan((a_[k] - lo) / 1000.0, (b_[k] - lo) / 1000.0,
                   color=WINDOW, alpha=0.28, lw=0, zorder=1,
                   label="the 0.5 s averaged for this stop" if n == 0 else None)

    ax.plot(t, y, color=INK, lw=1.5, zorder=4, solid_capstyle="round")

    for k in ks:
        ax.annotate("target %.3f mm" % iso[k]["cmd_mm"],
                    xy=(0.5 * (a_[k] + b_[k] - 2 * lo) / 1000.0, pos[a_[k]]),
                    xytext=(0, 9), textcoords="offset points",
                    color=INK2, size=8.8, ha="center", va="bottom")

    ax.set_xlim(t[0], t[-1])
    ax.set_xlabel("time (s)")
    ax.set_ylabel("measured position (mm)")
    ax.legend(loc="upper left", handlelength=1.6, labelspacing=0.5)

    # ---- inset: the closing edge of the middle window, sample by sample ----
    k = k0 + NSTEP // 2
    seg = moving[a_[k]:b_[k] + 3000]
    onset = a_[k] + int(np.argmax(seg)) if seg.any() else b_[k]
    i0, i1 = b_[k] - 90, b_[k] + 22
    rest = pos[onset - 550:onset - 50].mean()

    axi = ax.inset_axes([0.660, 0.118, 0.318, 0.390])
    ti = np.arange(i0, i1) - b_[k]
    yi = (pos[i0:i1] - rest) * 1000.0
    still = ti <= (onset - b_[k])
    axi.axvspan(onset - b_[k], 0, color=WINDOW, alpha=0.28, lw=0, zorder=1)
    axi.plot(ti, yi, color=INK3, lw=0.9, zorder=3)
    axi.scatter(ti[still], yi[still], s=11, color=INK, lw=0, zorder=4,
                label="samples at rest")
    bad = ~still & (ti <= 0)
    axi.scatter(ti[bad], yi[bad], s=15, color=WINDOW, lw=0, zorder=5,
                label="moving samples, inside the window")
    after = ti > 0
    axi.scatter(ti[after], yi[after], s=11, facecolor="none",
                edgecolor=INK3, lw=0.7, zorder=4, label="after the window")
    axi.axvline(0, color=INK3, lw=0.9, zorder=2)

    axi.set_xlim(ti[0], ti[-1])
    axi.set_xlabel("ms, relative to the window closing", size=8.5)
    axi.set_ylabel("%s from rest" % MU, size=8.5)
    axi.tick_params(labelsize=8)
    axi.set_title("the middle window's closing edge, every 1 ms sample",
                  loc="left", size=9.0, color=INK, pad=5)
    axi.legend(loc="upper left", fontsize=7.8, handlelength=1.1,
               labelspacing=0.35, borderpad=0.3)
    ax.indicate_inset([(i0 - lo) / 1000.0, pos[i0] - 0.012,
                       (i1 - i0) / 1000.0, 0.060],
                      inset_ax=axi, edgecolor=INK3, lw=0.9, alpha=0.9)

    titleblock(
        fig,
        "ISO 230-2 v3: which samples each measured position was averaged from",
        "Three consecutive stops of one pass. The shaded half-second is the "
        "window the published evaluation averaged, placed on the trace by the "
        "same fitted clock it used.\n"
        "%s  ·  %s  ·  stops %d–%d of 240, %s pass"
        % (Path(csv).name, Path(log).name, k0 + 1, k0 + NSTEP,
           "up" if dir_[k0] > 0 else "down"),
        "This window ends %d ms after the stage has already left the target, "
        "so its last samples are mid-move — across all 240 stops the "
        "median is %.0f ms." % (b_[k] - onset, over_med))

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        p = OUT / ("v3_record_windows.%s" % ext)
        fig.savefig(p)
        print("wrote %s" % p)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
