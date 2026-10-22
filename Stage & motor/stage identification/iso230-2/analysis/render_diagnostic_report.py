#!/usr/bin/env python3
"""Three figures on the open-loop diagnostic (IsoTest.csv, 2026-10-21).

    1  the experiment   what was run, and what question it was built to answer
    2  the results      the recorded data for all four profiles, side by side
    3  the limit        2.500 mm/s, where it comes from, what it does

Written to plots/temp/. No terminal log was saved for this run, so the profile
boundaries are recovered from the trace rather than assumed: the transfers are
found by their length, grouped four to a profile in run order, and the three
plain transfers in each group kept (the fourth has the home approach appended
and runs much longer).
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.ticker import MaxNLocator

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa

SRC = Path(os.environ.get("ISO_TEST_CSV", _HERE.parent.parent / "IsoTest.csv"))
OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))

HEADER_ROWS = 23
WRAP = 2 ** 32
NM_PER_COUNT = 1.0            # Resolution AquadB = 1000 pm
CLOCK_NS = 400.0              # AquadB clock
CEILING = 1e9 / CLOCK_NS * NM_PER_COUNT / 1e6      # 2.500 mm/s
CEILING_NM = CEILING * 1e3                         # nm per 1 ms sample
L_MM = 59.338                 # the transfer distance
HOLD_S = 10.0

PROFILES = [("A", 5.0, 40.0, "the setting v2 ran"),
            ("B", 5.0, 10.0, "same speed, 4x less acceleration"),
            ("C", 2.0, 40.0, "same acceleration, 2.5x less speed"),
            ("D", 1.0, 5.0, "both well derated")]
COLOUR = {"A": "#99483e", "B": "#c88e5c", "C": "#5f8bc5", "D": "#356449"}
RED, GREEN, BLUE = "#b3261e", "#1f7a4d", "#2a78d6"


# ------------------------------------------------------------------ loading

def load(path=SRC):
    idx, val = [], []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i < HEADER_ROWS:
                continue
            line = line.strip()
            if not line or line.startswith("EOF"):
                continue
            p = line.split("\t")
            if len(p) == 2:
                idx.append(int(p[0]))
                val.append(int(p[1]))
    val = np.array(val, np.int64)
    dv = np.diff(val)
    dv[dv > WRAP // 2] -= WRAP          # UINT32 rollover, not motion
    dv[dv < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(dv)
    return np.array(idx, np.int64), pos


def find_transfers(pos, min_mm=50.0):
    d = np.diff(pos)
    sp = np.convolve(np.abs(d).astype(float), np.ones(101) / 101,
                     mode="same") / 1e3
    moving = sp > 0.03
    out, i = [], 0
    while i < len(moving):
        if moving[i]:
            j = i
            while j < len(moving) and moving[j]:
                j += 1
            if abs(pos[j] - pos[i]) / 1e6 > min_mm:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def measure(pos, runs):
    d = np.diff(pos).astype(float)
    per = len(runs) // len(PROFILES)
    rows = []
    for k, (name, v, acc, aim) in enumerate(PROFILES):
        block = runs[k * per:(k + 1) * per]
        clean = sorted(block, key=lambda ab: ab[1] - ab[0])[:3]
        a, b = clean[0]
        seg = np.abs(d[a + 2000:b - 2000])
        # deceleration present at the end of the recorded move
        sm = np.convolve(np.abs(d[b - 1500:b + 100]), np.ones(15) / 15,
                         mode="valid")
        pk = float(np.median(sm[:400]))
        hi = np.where(sm >= 0.9 * pk)[0]
        lo = np.where(sm <= 0.10 * pk)[0]
        lo = lo[lo > hi[-1]] if len(hi) else lo
        # stationary time recorded inside the following hold
        tail = np.convolve(np.abs(d[b:b + 12000]), np.ones(51) / 51,
                           mode="same")
        n = 0
        for s in tail < 20.0:
            if s:
                n += 1
            elif n > 200:
                break
            else:
                n = 0
        rows.append(dict(
            name=name, v=v, acc=acc, aim=aim, slice=(a, b),
            block=(block[0][0], block[-1][1]),
            dist=float(np.mean([abs(pos[y] - pos[x]) / 1e6 for x, y in clean])),
            recorded=float(np.mean([(y - x) / 1000.0 for x, y in clean])),
            physics=L_MM / v + v / acc,
            drain=L_MM / CEILING,
            clip=100.0 * float(np.mean(seg == CEILING_NM)),
            decel=int(lo[0] - hi[-1]) if len(lo) and len(hi) else -1,
            decel_needed=1000.0 * min(v, CEILING) / acc,
            still=n / 1000.0))
    return rows


def head(fig, n, title, sub):
    fig.text(0.045, 0.962, f"{n}.  {title}", size=20, weight="bold",
             color=INK, va="top")
    fig.text(0.045, 0.908, sub, size=12, color=INK2, va="top")


# ================================================================= figure 1
def fig_experiment(pos, rows):
    """Full recorded trace, then each phase at a common scale.

    No annotation beyond panel labels, settings and axes: this one is built
    to sit in a paper, where the caption carries the explanation.
    """
    fig = plt.figure(figsize=(7.6, 5.8))
    t = np.arange(len(pos)) / 1000.0
    mm = (pos - pos.min()) / 1e6

    ax = fig.add_axes([0.093, 0.590, 0.890, 0.360])
    ax.plot(t, mm, color=INK, lw=0.7)
    for r in rows:
        a, b = r["block"]
        ax.axvspan(a / 1000.0, b / 1000.0, color=COLOUR[r["name"]],
                   alpha=0.17, lw=0)
        ax.text((a + b) / 2000.0, 134, r["name"], ha="center", va="center",
                size=9.5, weight="bold", color=COLOUR[r["name"]])
    ax.set_xlim(0, t[-1])
    ax.set_ylim(-6, 146)
    ax.set_xlabel("$t$ (s)", size=9)
    ax.set_ylabel("$x$ (mm)", size=9)
    ax.tick_params(labelsize=8)
    ax.grid(False)

    for i, r in enumerate(rows):
        bx = fig.add_axes([0.093 + 0.2255 * i, 0.098, 0.196, 0.330])
        a, b = r["block"]
        seg = mm[a:b]
        # Each phase fills its own axes. They differ in duration, so one
        # shared time axis would strand the fast phases in white space;
        # the durations are read off the tick labels instead.
        bx.plot(np.arange(len(seg)) / 1000.0, seg, color=COLOUR[r["name"]],
                lw=0.7)
        bx.set_xlim(0, (b - a) / 1000.0)
        bx.set_ylim(-4, 70)
        bx.xaxis.set_major_locator(
            MaxNLocator(nbins=4, steps=[1, 2, 5, 10], integer=True))
        bx.tick_params(labelsize=8)
        bx.set_xlabel("$t$ (s)", size=9)
        if i == 0:
            bx.set_ylabel("$x$ (mm)", size=9)
        else:
            bx.set_yticklabels([])
        bx.grid(False)
        bx.text(0.0, 1.20, r["name"], transform=bx.transAxes, size=9.5,
                weight="bold", color=COLOUR[r["name"]], va="top")
        bx.text(0.0, 1.075,
                f"$v$ = {r['v']:g} mm s$^{{-1}}$,  "
                f"$a$ = {r['acc']:g} mm s$^{{-2}}$",
                transform=bx.transAxes, size=7.4, color=INK2, va="top")
    return fig


# ================================================================= figure 2
def fig_results(pos, rows):
    """The four recorded transfers as a shared-axes matrix.

    One time axis and one speed axis serve all four panels, so the traces
    compare directly and only the outer edges carry labels. Everything the
    panels used to spell out in text belongs in the caption.
    """
    d = np.diff(pos).astype(float)
    fig = plt.figure(figsize=(7.6, 5.0))

    left, right, bottom, top = 0.105, 0.985, 0.115, 0.965
    gx, gy = 0.018, 0.098
    w = (right - left - gx) / 2
    h = (top - bottom - gy) / 2
    axes = []
    for i in range(4):
        col, row = i % 2, i // 2
        x0 = left + col * (w + gx)
        y0 = bottom + (1 - row) * (h + gy)
        axes.append(fig.add_axes([x0, y0, w, h]))

    for i, (r, ax) in enumerate(zip(rows, axes)):
        a, b = r["slice"]
        sm = np.convolve(np.abs(d[a:b + 13000]), np.ones(101) / 101,
                         mode="same") / 1e3
        t = np.arange(len(sm)) / 1000.0
        ax.axhline(CEILING, color=INK3, lw=0.8, ls=(0, (4, 3)))
        # The dwell is the stretch the trace spends at zero after the
        # move. Marking its length makes the cost of the lag visible:
        # the same 10 s hold was commanded in every profile.
        x0 = r["recorded"]
        x1 = x0 + r["still"]
        ax.axvspan(x0, x1, ymin=0, ymax=0.40, color=INK3, alpha=0.18,
                   lw=0)
        ax.plot(t, sm, color=COLOUR[r["name"]], lw=0.8)
        ax.annotate("", xy=(x0, 0.70), xytext=(x1, 0.70),
                    arrowprops=dict(arrowstyle="<->", color=INK,
                                    lw=0.8))
        ax.text(0.5 * (x0 + x1), 0.84, f"{r['still']:.1f} s",
                ha="center", va="bottom", size=7.6, color=INK)
        # x fits each record, y stays common: the traces then fill their
        # panels instead of trailing off into white space, and the
        # durations are read from the tick labels.
        ax.set_xlim(0, len(sm) / 1000.0)
        ax.xaxis.set_major_locator(
            MaxNLocator(nbins=5, steps=[1, 2, 5, 10], integer=True))
        ax.set_ylim(0, 3.1)
        ax.tick_params(labelsize=8)
        ax.grid(False)
        if i % 2 == 1:
            ax.set_yticklabels([])
        # One line along the top: the band above 2.6 mm/s is the only strip
        # free of data in all four panels once each fills its own width.
        ax.text(0.018, 0.965, r["name"], transform=ax.transAxes,
                va="top", size=9.5, weight="bold",
                color=COLOUR[r["name"]])
        ax.text(0.080, 0.955,
                f"$v$ = {r['v']:g} mm s$^{{-1}}$,  "
                f"$a$ = {r['acc']:g} mm s$^{{-2}}$",
                transform=ax.transAxes, va="top", size=7.4,
                color=INK2)
        if i == 0:
            ax.text(0.26, CEILING + 0.07, "2.500 mm s$^{-1}$",
                    transform=ax.get_yaxis_transform(), size=7.2,
                    color=INK3)

    fig.text((left + right) / 2, 0.028, "$t$ (s)", size=9, ha="center")
    fig.text(0.022, (bottom + top) / 2, "$|\\dot{x}|$ (mm s$^{-1}$)", size=9,
             va="center", rotation="vertical")
    return fig


# ================================================================= figure 3
def fig_limit(pos, rows):
    """Saturation, closure, and the lag the saturation causes.

    (a) frequency distribution of the per-sample displacement, (b) the
    closure of every transfer against the band they all fall inside, and
    (c) the recorded speed of the fastest profile, with the intervals marked.
    No prose: the caption carries it.
    """
    d = np.diff(pos).astype(float)
    fig = plt.figure(figsize=(7.6, 6.3))
    ax = fig.add_axes([0.093, 0.610, 0.375, 0.330])
    bx = fig.add_axes([0.610, 0.610, 0.375, 0.330])
    cx = fig.add_axes([0.093, 0.078, 0.892, 0.425])

    # --- (a) per-sample displacement, one consistent binning --------------
    steady = np.concatenate([np.abs(d[r["slice"][0] + 2000:
                                     r["slice"][1] - 2000]) for r in rows])
    bins = np.arange(0, CEILING_NM + 26, 25)
    ax.hist(steady, bins=bins, color=INK3, lw=0)
    ax.axvline(CEILING_NM, color=INK, lw=0.9, ls=(0, (4, 3)))
    ax.set_yscale("log")
    ax.set_xlim(0, CEILING_NM + 120)
    ax.set_ylim(3e2, 3e5)
    ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
    ax.set_xlabel("displacement per 1 ms sample (nm)", size=9)
    ax.set_ylabel("samples", size=9)
    ax.tick_params(labelsize=8)
    ax.grid(False)
    ax.text(0.015, 0.965, "(a)", transform=ax.transAxes, size=9.5,
            weight="bold", color=INK, va="top")

    # --- (b) closure of every transfer ------------------------------------
    dists = np.array([r["dist"] for r in rows])
    lo, hi = dists.min(), dists.max()
    band = 1000.0 * (hi - lo)
    xs = np.arange(len(rows))
    bx.axhspan(lo, hi, color=INK3, alpha=0.22, lw=0)
    bx.plot(xs, dists, "o", ms=7, mfc=[COLOUR[r["name"]] for r in rows][0],
            mec="none", alpha=0)        # placeholder to fix autoscale order
    for i, r in enumerate(rows):
        bx.plot([i], [dists[i]], "o", ms=7, mfc=COLOUR[r["name"]],
                mec=SURFACE, mew=0.8)
    mid = 0.5 * (lo + hi)
    bx.annotate("", xy=(3.78, lo), xytext=(3.78, hi),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
    bx.text(3.90, mid, f"{band:.1f} µm", ha="left", va="center", size=8,
            color=INK)
    bx.set_xlim(-0.5, 4.75)
    bx.set_ylim(mid - 0.0035, mid + 0.0035)
    bx.set_xticks(xs)
    bx.set_xticklabels([r["name"] for r in rows], size=9)
    bx.set_ylabel("transfer distance (mm)", size=9)
    bx.tick_params(labelsize=8)
    bx.grid(False)
    bx.text(0.015, 0.965, "(b)", transform=bx.transAxes, size=9.5,
            weight="bold", color=INK, va="top")

    # --- (c) the lag the saturation causes --------------------------------
    r = rows[0]
    a, b = r["slice"]
    sm = np.convolve(np.abs(d[a:a + 34000]), np.ones(101) / 101,
                     mode="same") / 1e3
    t = np.arange(len(sm)) / 1000.0
    arrive, done = r["physics"], r["recorded"]
    still_end = done + r["still"]

    cx.axvspan(arrive, arrive + HOLD_S, color=INK3, alpha=0.14, lw=0)
    cx.plot(t, sm, color=COLOUR["A"], lw=1.0)
    cx.axhline(CEILING, color=INK, lw=0.9, ls=(0, (4, 3)))
    for xv in (arrive, done, still_end):
        cx.axvline(xv, color=INK3, lw=0.7, ls=(0, (2, 3)))

    def interval(x0, x1, y, label, colour=INK):
        cx.annotate("", xy=(x0, y), xytext=(x1, y),
                    arrowprops=dict(arrowstyle="<->", color=colour, lw=0.9))
        cx.text(0.5 * (x0 + x1), y + 0.075, label, ha="center", va="bottom",
                size=8.5, color=colour)

    interval(arrive, done, 2.86, f"{done - arrive:.1f} s  recording lag")
    interval(arrive, arrive + HOLD_S, 0.33,
             f"{HOLD_S:.0f} s  commanded hold")
    cx.annotate(f"{r['still']:.1f} s  stationary",
                xy=(0.5 * (done + still_end), 0.10),
                xytext=(still_end + 2.2, 0.72), size=8.5, color=INK,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=INK3, lw=0.7))
    cx.set_xlabel("$t$ (s)", size=9)
    cx.set_ylabel("recorded speed (mm s$^{-1}$)", size=9)
    cx.set_xlim(0, 33)
    cx.set_ylim(0, 3.25)
    cx.tick_params(labelsize=8)
    cx.grid(False)
    cx.text(arrive - 0.35, 2.17, "stage stops", ha="right", va="center",
            size=8.5, color=INK2)
    cx.text(done + 0.35, 2.17, "record stops", ha="left", va="center",
            size=8.5, color=INK2)
    cx.text(0.008, 0.965, "(c)", transform=cx.transAxes, size=9.5,
            weight="bold", color=INK, va="top")
    cx.text(0.047, 0.965, r["name"], transform=cx.transAxes, size=9.5,
            weight="bold", color=COLOUR[r["name"]], va="top")
    cx.text(0.076, 0.955,
            f"$v$ = {r['v']:g} mm s$^{{-1}}$,  "
            f"$a$ = {r['acc']:g} mm s$^{{-2}}$",
            transform=cx.transAxes, size=7.6, color=INK2, va="top")
    return fig


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    style()
    idx, pos = load()
    rows = measure(pos, find_transfers(pos))
    print(f"{'':2} {'cmd':>11} {'dist mm':>8} {'physics':>8} {'recorded':>9} "
          f"{'clip%':>6} {'decel':>7} {'needed':>7} {'still':>7}")
    for r in rows:
        print(f"{r['name']:2} {r['v']:5.1f}/{r['acc']:<5.0f} {r['dist']:8.3f} "
              f"{r['physics']:8.2f} {r['recorded']:9.2f} {r['clip']:6.1f} "
              f"{r['decel']:7d} {r['decel_needed']:7.0f} {r['still']:7.1f}")
    OUT.mkdir(parents=True, exist_ok=True)
    for stem, fn in [("01_experiment", fig_experiment),
                     ("02_results", fig_results),
                     ("03_velocity_limit", fig_limit)]:
        fig = fn(pos, rows)
        fig.savefig(OUT / f"{stem}.png")
        fig.savefig(OUT / f"{stem}.pdf")
        plt.close(fig)
        print(f"  wrote {stem}.png / .pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
