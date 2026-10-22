#!/usr/bin/env python3
"""Why v3 still saturates while moving, and why its results survive it.

v3 commands 2 mm/s, under the 2.500 mm/s AquadB ceiling, yet 46 % of its
moving samples still sit exactly on the ceiling. The reason is that the
commanded speed is a mean, not an instantaneous value: the firmware stalls
its step-pulse train every 20 ms to write a telemetry line, the rotor rings
around the commanded position during the stall, and the velocity swings the
full +-2.500 mm/s range. The link clips those peaks and queues the remainder.

The queue is tiny and transient, and the 2 s settle window drains it long
before the ISO record window opens, which is why none of the 221 record
windows contains a saturated sample.

Writes plots/temp/04_v3_saturation_impact.{png,pdf}.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa

SRC = Path(os.environ.get("ISO_V3_CSV", _HERE.parent.parent / "ISO230-2_v3.csv"))
OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))
HEADER_ROWS = 23
WRAP = 2 ** 32
CEIL = 2500.0            # nm per 1 ms sample = 2.500 mm/s
REC_MS = 500             # t_record_s = 0.5
TELEM_MS = 20            # SAMPLE_INTERVAL_MS in the firmware
BLUE, RED, GREEN = "#5f8bc5", "#99483e", "#356449"


def load(path=SRC):
    val = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i < HEADER_ROWS:
                continue
            line = line.strip()
            if not line or line.startswith("EOF"):
                continue
            p = line.split("\t")
            if len(p) == 2:
                val.append(int(p[1]))
    val = np.array(val, np.int64)
    d = np.diff(val)
    d[d > WRAP // 2] -= WRAP
    d[d < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(d)
    return pos


def stops(pos):
    d = np.diff(pos)
    sp = np.convolve(np.abs(d).astype(float), np.ones(51) / 51,
                     mode="same") / 1e3
    still = sp < 0.02
    out, i = [], 0
    while i < len(still):
        if still[i]:
            j = i
            while j < len(still) and still[j]:
                j += 1
            if 1500 < j - i < 5000:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    style()
    pos = load()
    d = np.diff(pos).astype(float)
    sat = np.abs(d) == CEIL
    runs = stops(pos)

    # margin from the last saturated sample to the window opening
    margins = []
    for a, b in runs:
        w0 = b - REC_MS
        prev = np.where(sat[:w0])[0]
        margins.append((w0 - prev[-1]) / 1000.0 if len(prev) else np.nan)
    margins = np.array(margins)
    in_window = sum(1 for a, b in runs if sat[b - REC_MS:b].any())

    # saturation rate against time before the window opens
    span = 5000
    acc = np.zeros(span + REC_MS)
    for a, b in runs:
        w0 = b - REC_MS
        s, e = w0 - span, w0 + REC_MS
        if s < 0:
            continue
        acc += sat[s:e]
    rate = 100.0 * acc / len(runs)
    tax = (np.arange(-span, REC_MS)) / 1000.0

    sp = np.convolve(np.abs(d), np.ones(51) / 51, mode="same") / 1e3
    mov = sp > 0.5
    print(f"v3: {len(runs)} record stops")
    print(f"  saturated samples: {100*sat.mean():.1f}% overall, "
          f"{100*sat[mov].mean():.1f}% while moving, "
          f"{100*sat[~mov].mean():.3f}% while stopped")
    print(f"  windows containing saturation: {in_window} of {len(runs)}")
    print(f"  margin to the window: min {np.nanmin(margins):.2f} s, "
          f"median {np.nanmedian(margins):.2f} s")

    # One representative ISO target move and the dwell that follows it.
    MOVE_T0, DWELL_T0, DWELL_T1 = 410_110, 413_440, 415_870
    ZOOM_T0, ZOOM_W = 411_500, 140

    fig = plt.figure(figsize=(7.6, 7.7))
    ax = fig.add_axes([0.093, 0.722, 0.890, 0.218])    # overview
    bx = fig.add_axes([0.093, 0.428, 0.890, 0.208])    # the zoom
    cx = fig.add_axes([0.093, 0.072, 0.400, 0.250])    # saturation vs window
    dx = fig.add_axes([0.610, 0.072, 0.375, 0.250])    # margin

    # --- (a) one move and its dwell, with the zoom located ----------------
    s0, s1 = MOVE_T0 - 200, DWELL_T1 + 250
    seg = d[s0:s1]
    tt = (np.arange(len(seg)) + s0 - MOVE_T0) / 1000.0
    dw0 = (DWELL_T0 - MOVE_T0) / 1000.0
    dw1 = (DWELL_T1 - MOVE_T0) / 1000.0
    ax.axvspan(dw0, dw1, color=INK3, alpha=0.13, lw=0)
    ax.axvspan(dw1 - REC_MS / 1000.0, dw1, color=GREEN, alpha=0.22, lw=0)
    for s in (+1, -1):
        ax.axhline(s * CEIL, color=INK, lw=0.9, ls=(0, (4, 3)))
    ax.axhline(0, color=INK3, lw=0.7)
    ax.plot(tt, seg, color=BLUE, lw=0.5)
    z0 = (ZOOM_T0 - MOVE_T0) / 1000.0
    ax.axvspan(z0, z0 + ZOOM_W / 1000.0, color=RED, alpha=0.30, lw=0)
    ax.annotate("(b)", xy=(z0 + 0.07, 2750), xytext=(z0 + 0.52, 3250),
                size=9, weight="bold", color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    ax.text(0.5 * (0 + dw0), -3800, "move", ha="center", size=8.6,
            color=INK2)
    ax.text(0.5 * (dw0 + dw1), -3800, f"dwell {(dw1-dw0):.2f} s",
            ha="center", size=8.6, color=INK2)
    ax.annotate(f"record window {REC_MS/1000:.1f} s",
                xy=(dw1 - 0.25, 1500), xytext=(dw1 - 1.55, 3100),
                size=8.6, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.9))
    ax.set_xlim(tt[0], tt[-1])
    ax.set_ylim(-CEIL * 1.78, CEIL * 1.78)
    ax.set_xlabel("time from the start of the move (s)", size=9)
    ax.set_ylabel("displacement\nper 1 ms (nm)", size=9)
    ax.tick_params(labelsize=8)
    ax.grid(False)
    ax.text(0.008, 0.965, "(a)", transform=ax.transAxes, size=9.5,
            weight="bold", color=INK, va="top")
    ax.text(0.052, 0.965, "one 5.1 mm ISO target move, commanded "
                          "2 mm s$^{-1}$",
            transform=ax.transAxes, size=8.6, color=INK2, va="top")

    # --- (b) 140 ms from inside that move ---------------------------------
    seg = d[ZOOM_T0:ZOOM_T0 + ZOOM_W]
    t = np.arange(len(seg))
    bx.axhline(0, color=INK3, lw=0.7)
    for s in (+1, -1):
        bx.axhline(s * CEIL, color=INK, lw=0.9, ls=(0, (4, 3)))
    first = int(np.argmin(seg[:TELEM_MS]))
    for k in range(0, 8):
        xv = first + TELEM_MS * k
        if xv <= t[-1]:
            bx.axvline(xv, color=INK3, lw=0.7, ls=(0, (2, 3)))
    bx.plot(t, seg, color=BLUE, lw=1.1, marker="o", ms=2.2)
    bx.set_xlim(0, ZOOM_W - 1)
    bx.set_ylim(-CEIL * 1.34, CEIL * 1.34)
    bx.set_xlabel("time within the move (ms)", size=9)
    bx.set_ylabel("displacement\nper 1 ms (nm)", size=9)
    bx.tick_params(labelsize=8)
    bx.grid(False)
    bx.text(0.008, 0.965, "(b)", transform=bx.transAxes, size=9.5,
            weight="bold", color=INK, va="top")
    bx.text(0.052, 0.965, "every 1 ms sample of the shaded strip above",
            transform=bx.transAxes, size=8.6, color=INK2, va="top")
    bx.text(0.008, 0.038, f"dotted lines {TELEM_MS} ms apart",
            transform=bx.transAxes, size=8.2, color=INK2)
    bx.text(0.998, CEIL + 210, "+2.500 mm s$^{-1}$",
            transform=bx.get_yaxis_transform(), ha="right", size=7.6,
            color=INK)
    bx.text(0.998, -CEIL - 760, "−2.500 mm s$^{-1}$",
            transform=bx.get_yaxis_transform(), ha="right", size=7.6,
            color=INK)

    # --- (c) where saturation sits relative to the record window ----------
    cx.axvspan(-REC_MS / 1000.0, 0, color=GREEN, alpha=0.22, lw=0)
    cx.plot(tax, np.convolve(rate, np.ones(15) / 15, mode="same"),
            color=BLUE, lw=1.2)
    cx.set_xlim(-5, 0.5)
    cx.set_ylim(-2, 62)
    cx.set_xlabel("time to the end of the record window (s)", size=9)
    cx.set_ylabel("stops with a saturated sample (%)", size=9)
    cx.tick_params(labelsize=8)
    cx.grid(False)
    cx.text(0.015, 0.965, "(c)", transform=cx.transAxes, size=9.5,
            weight="bold", color=INK, va="top")
    cx.annotate("record window", xy=(-0.25, 4), xytext=(-1.95, 17),
                size=8.4, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.9))
    cx.text(0.97, 0.90, f"{in_window} of {len(runs)} windows\ncontain one",
            transform=cx.transAxes, ha="right", va="top", size=8.4,
            color=GREEN)

    # --- (d) the margin that protects the measurement ---------------------
    good = margins[~np.isnan(margins)]
    dx.hist(good, bins=np.arange(1.0, 2.35, 0.05), color=INK3, lw=0)
    dx.axvline(good.min(), color=RED, lw=1.2)
    dx.set_xlabel("margin before the window (s)", size=9)
    dx.set_ylabel("stops", size=9)
    dx.tick_params(labelsize=8)
    dx.grid(False)
    dx.text(0.015, 0.965, "(d)", transform=dx.transAxes, size=9.5,
            weight="bold", color=INK, va="top")
    dx.annotate(f"min {good.min():.2f} s",
                xy=(good.min(), dx.get_ylim()[1] * 0.45),
                xytext=(good.min() + 0.12, dx.get_ylim()[1] * 0.62),
                size=8.4, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "04_v3_saturation_impact.png"
    fig.savefig(png)
    fig.savefig(OUT / "04_v3_saturation_impact.pdf")
    plt.close(fig)
    print(f"\nwrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
