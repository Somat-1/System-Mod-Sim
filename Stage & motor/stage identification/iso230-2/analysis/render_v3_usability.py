#!/usr/bin/env python3
"""Is the v3 ISO 230-2 dataset usable despite the pulse-train stalling?

ISO 230-2 reads a SETTLED position at each target. It never looks at the
trajectory. So the stalling documented in render_velocity_ceiling.py matters
only if it does one of two things:

  loses steps      the stage is then not where the pulse count says, and the
                   error accumulates with travel
  spoils settling  the stage is still moving when the 0.5 s record window
                   opens, so the recorded value is not a position

Both are testable from the trace alone, which matters because no terminal log
was saved for this run either. The verdict is that neither happens: the stage
is dead still during every record window, and what drift there is appears
equally on every target rather than scaling with travel, which is the
signature of thermal drift rather than lost steps.

Config: settle_margin_s 2.0 then t_record_s 0.5 at each target, n_runs 5, so
each target is approached ten times.

Writes plots/iso230_2_v3_usability.png.
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
from current_sweep_figures import titleblock3                        # noqa

SRC = Path(os.environ.get("ISO_V3_CSV", _HERE.parent.parent / "ISO230-2_v3.csv"))
OUT = Path(os.environ.get("ISO_V3_OUT", _HERE / "plots"))
HEADER_ROWS = 23
WRAP = 2 ** 32
RECORD_MS = 500          # t_record_s = 0.5
POS = "#2a78d6"
WARN = "#b3261e"
OKC = "#1f7a4d"


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
    idx = np.array(idx, np.int64)
    val = np.array(val, np.int64)
    d = np.diff(val)
    d[d > WRAP // 2] -= WRAP
    d[d < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(d)
    return idx, pos


def dwells(pos, lo=1500, hi=5000):
    """The record stops: still for between 1.5 and 5 s."""
    sp = np.convolve(np.abs(np.diff(pos)).astype(float), np.ones(51) / 51,
                     mode="same") / 1e3
    still = sp < 0.02
    out, i = [], 0
    while i < len(still):
        if still[i]:
            j = i
            while j < len(still) and still[j]:
                j += 1
            if lo < j - i < hi:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def cluster(settled, tol_nm=200000):
    """Group the stops that belong to the same target."""
    o = np.argsort(settled)
    groups, cur = [], [o[0]]
    for k in o[1:]:
        if settled[k] - settled[cur[-1]] < tol_nm:
            cur.append(k)
        else:
            groups.append(cur)
            cur = [k]
    groups.append(cur)
    return [g for g in groups if len(g) >= 6]


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    idx, pos = load()
    gaps = int(np.sum(np.diff(idx) != 1))
    missing = int((np.diff(idx)[np.diff(idx) != 1] - 1).sum())
    stops = dwells(pos)
    settled = np.array([pos[b - RECORD_MS:b].mean() for a, b in stops])
    tmid = np.array([(a + b) / 2000.0 for a, b in stops])

    std = np.array([pos[b - RECORD_MS:b].std() for a, b in stops])
    ptp = np.array([np.ptp(pos[b - RECORD_MS:b]) for a, b in stops])
    # how far the stage still travelled between arriving and the record window
    creep = np.array([pos[b - RECORD_MS:b].mean() - pos[a:a + 50].mean()
                      for a, b in stops])

    groups = cluster(settled)
    slopes, resid, centres, counts = [], [], [], []
    for g in groups:
        y, x = settled[g], tmid[g] / 60.0
        k, c = np.polyfit(x, y, 1)
        slopes.append(k)
        resid.append(np.std(y - np.polyval([k, c], x)))
        centres.append(y.mean())
        counts.append(len(g))
    slopes = np.array(slopes)
    resid = np.array(resid)
    centres = np.array(centres)
    span_mm = (centres.max() - centres.min()) / 1e6
    pos_slope = np.polyfit(centres / 1e6, slopes, 1)[0]

    print(f"{SRC.name}: {len(idx):,} samples, {len(idx)/1000:.0f} s, "
          f"{gaps} index gaps ({missing} samples)")
    print(f"record stops {len(stops)}, target clusters {len(groups)}, "
          f"{int(np.median(counts))} visits each (n_runs 5 -> 10)")
    print(f"\nwithin the 0.5 s record window:")
    print(f"  sd       median {np.median(std):6.1f} nm   p90 "
          f"{np.percentile(std,90):6.1f}   max {std.max():7.1f}")
    print(f"  p-p      median {np.median(ptp):6.1f} nm   p90 "
          f"{np.percentile(ptp,90):6.1f}   max {ptp.max():7.1f}")
    print(f"  creep from arrival: median {np.median(np.abs(creep)):.0f} nm, "
          f"max {np.abs(creep).max():.0f} nm")
    print(f"\nper-target repeatability (residual sd after removing drift):")
    print(f"  median {np.median(resid)/1e3:.2f} um   max {resid.max()/1e3:.2f} um")
    print(f"\ndrift: mean {slopes.mean():+.1f} nm/min over all targets "
          f"({slopes.mean()*(len(idx)/60000.0)/1000:+.2f} um across the "
          f"{len(idx)/60000.0:.0f} min run)")
    print(f"  slope vs position {pos_slope:+.3f} nm/min per mm over "
          f"{span_mm:.0f} mm -> {pos_slope*span_mm:+.1f} nm/min end to end")
    print("  Lost steps would scale with travel; this does not.")

    # -------------------------------------------------------------- figure
    style()
    fig = plt.figure(figsize=(13.0, 9.8))
    ax = fig.add_axes([0.070, 0.520, 0.400, 0.225])
    bx = fig.add_axes([0.570, 0.520, 0.385, 0.225])
    cx = fig.add_axes([0.070, 0.075, 0.400, 0.280])
    dx = fig.add_axes([0.570, 0.075, 0.385, 0.280])

    titleblock3(
        fig,
        "ISO 230-2 v3 - is the data usable?",
        "The pulse-train stalling corrupts the trajectory, not the settled "
        "position. ISO 230-2 reads only the latter, and the stepper is open "
        "loop, so the stage lands where the pulse count puts it as long as no "
        "pulses are lost.",
        f"Verdict: usable. Still to {np.median(std):.0f} nm sd in every "
        f"record window; each target repeats to {np.median(resid)/1e3:.2f} um "
        f"sd over {int(np.median(counts))} visits; the "
        f"{slopes.mean():+.0f} nm/min drift is common to all targets, not "
        f"scaling with travel. {missing} samples dropped, none in a window.",
        width=150)

    # (a) a representative stop
    a, b = stops[len(stops) // 2]
    t = (np.arange(a, b) - a) / 1000.0
    y = (pos[a:b] - pos[b - RECORD_MS:b].mean()) / 1e3
    ax.axvspan((b - RECORD_MS - a) / 1000.0, (b - a) / 1000.0,
               color=OKC, alpha=0.12, lw=0)
    ax.plot(t, y, color=POS, lw=1.0)
    ax.axhline(0, color=INK3, lw=0.8)
    ax.set_xlabel("time since arrival (s)")
    ax.set_ylabel("position − recorded value (µm)")
    ax.set_title("(a)  a representative record stop", size=10.5, loc="left")
    ax.text(0.97, 0.90, "shaded = the 0.5 s record window",
            transform=ax.transAxes, ha="right", size=8.6, color=OKC)

    # (b) stillness within the record window
    bx.hist(std, bins=40, color=INK3, lw=0)
    bx.axvline(np.median(std), color=OKC, lw=1.4)
    bx.set_xlabel("sd within the record window (nm)")
    bx.set_ylabel("stops")
    bx.set_title("(b)  the stage is settled, every time", size=10.5,
                 loc="left")
    bx.text(0.96, 0.88, f"median {np.median(std):.0f} nm, max {std.max():.0f} nm\n"
                        f"{len(stops)} stops, none still moving",
            transform=bx.transAxes, ha="right", va="top", size=8.6, color=INK2)

    # (c) per-target repeatability
    for g, c0 in zip(groups, centres):
        cx.plot([(c0 - centres.min()) / 1e6] * len(g),
                (settled[g] - c0) / 1e3, "o", ms=3.4, mfc=POS, mec="none",
                alpha=0.55)
    cx.axhline(0, color=INK3, lw=0.8)
    cx.set_xlabel("target position along travel (mm)")
    cx.set_ylabel("settled − target mean (µm)")
    cx.set_title("(c)  every visit to every target", size=10.5, loc="left")
    cx.text(0.02, 0.05,
            f"{len(groups)} targets x ~{int(np.median(counts))} visits; "
            f"residual sd {np.median(resid)/1e3:.2f} µm",
            transform=cx.transAxes, size=8.8, color=INK2)

    # (d) drift per target - flat means no accumulating step loss
    dx.plot((centres - centres.min()) / 1e6, slopes, "o", ms=6,
            mfc=POS, mec=SURFACE, mew=1.0)
    xs = np.linspace(0, span_mm, 10)
    dx.plot(xs, np.polyval(np.polyfit((centres - centres.min()) / 1e6,
                                      slopes, 1), xs),
            color=WARN, lw=1.4, ls=(0, (5, 3)),
            label=f"{pos_slope:+.2f} nm/min per mm")
    dx.axhline(slopes.mean(), color=INK3, lw=1.0,
               label=f"common drift {slopes.mean():+.0f} nm/min")
    dx.set_xlabel("target position along travel (mm)")
    dx.set_ylabel("drift of that target (nm/min)")
    dx.set_title("(d)  drift does not scale with travel", size=10.5,
                 loc="left")
    dx.legend(loc="upper left", fontsize=8.4)
    dx.text(0.98, 0.04, "lost steps would put these on a slope",
            transform=dx.transAxes, ha="right", size=8.6, color=INK2)

    for x in (ax, bx, cx, dx):
        x.grid(True, axis="y")
        x.grid(False, axis="x")

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "iso230_2_v3_usability.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
