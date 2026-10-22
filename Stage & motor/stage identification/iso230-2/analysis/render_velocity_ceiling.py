#!/usr/bin/env python3
"""How the 2.500 mm/s ceiling was observed, and what it is not.

The open-loop diagnostic (IsoTest.csv, 2026-10-21) repeats the same 59.338 mm
transfer under four commanded motion profiles. The four were chosen so that
the two candidate causes move independently:

    A  5.0 mm/s  40 mm/s^2     the profile v2 ran
    B  5.0 mm/s  10 mm/s^2     same speed, acceleration cut 4x
    C  2.0 mm/s  40 mm/s^2     same acceleration, speed cut 2.5x
    D  1.0 mm/s   5 mm/s^2     both derated

That is the whole experiment: if the limit were torque or resonance it must
respond to A -> B; if it were a speed or data-rate limit it must respond to
A -> C. Only A -> C moves it, and even C does not clear it.

No terminal log with DIAG markers was saved, so the earlier figure had to
assume which block was which. This script does not: the transfers are found
from the trace, and the four groups are then identified by their own duration
and ramp time, which reproduce the commanded order independently.

Three things are measured, each answering one question:

  distance per transfer   is the COUNTER losing counts? (no - all 12
                          transfers close to within 1 um, while the clipped
                          fraction ranges from 2% to 99.6%)
  ramp time A vs B        did the acceleration setting take effect at all?
                          (yes - B's ramp is ~2x A's, and the ceiling is
                          unchanged, so acceleration is not the cause)
  step quantisation       what IS the ceiling? (exactly 4 microsteps per 1 ms
                          sample: 4 x 625 nm = 2500 nm = 2.500 mm/s)

Writes plots/iso_velocity_ceiling_explained.png.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa
from current_sweep_figures import titleblock3                        # noqa

SRC = Path(os.environ.get("ISO_TEST_CSV", _HERE.parent.parent / "IsoTest.csv"))
OUT = Path(os.environ.get("ISO_CEILING_OUT", _HERE / "plots"))

HEADER_ROWS = 23
WRAP = 2 ** 32
DT = 1e-3                       # s, the export's SampleTime

# Geometry, from config/iso_config.json: lead 2 mm, 200 steps/rev, 1/16.
FS_PER_MM = 100.0
MRES = 16
NM_PER_USTEP = 1e6 / (FS_PER_MM * MRES)     # 625 nm
CEILING_NM = 4 * NM_PER_USTEP               # 2500 nm per 1 ms sample

# commanded (speed mm/s, acceleration mm/s^2), in the order the config runs them
PROFILES = (("A", 5.0, 40.0), ("B", 5.0, 10.0),
            ("C", 2.0, 40.0), ("D", 1.0, 5.0))
COLOUR = {"A": "#b3261e", "B": "#e08214", "C": "#2a78d6", "D": "#1f7a4d"}


# ------------------------------------------------------------------- loading

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
            if len(p) != 2:
                continue
            idx.append(int(p[0]))
            val.append(int(p[1]))
    return np.array(idx, np.int64), np.array(val, np.int64)


def unwrap_u32(v):
    d = np.diff(v)
    d[d > WRAP // 2] -= WRAP
    d[d < -WRAP // 2] += WRAP
    out = np.empty_like(v)
    out[0] = v[0]
    out[1:] = v[0] + np.cumsum(d)
    return out


# ------------------------------------------------------- finding the transfers

def find_transfers(pos, min_mm=50.0):
    """The long constant-velocity transfers, found from the trace itself."""
    d = np.diff(pos)
    speed = np.convolve(np.abs(d).astype(float), np.ones(101) / 101,
                        mode="same") / 1e3                       # mm/s
    moving = speed > 0.03
    runs, i = [], 0
    while i < len(moving):
        if moving[i]:
            j = i
            while j < len(moving) and moving[j]:
                j += 1
            if abs(pos[j] - pos[i]) / 1e6 > min_mm:
                runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def measure(pos, a, b):
    """Per-transfer numbers. The steady part excludes 2 s of ramp each end."""
    d = np.diff(pos)
    seg = np.abs(d[a + 2000:b - 2000])
    m = {}
    m["dist_mm"] = abs(pos[b] - pos[a]) / 1e6
    m["dur_s"] = (b - a) * DT
    m["mean_mm_s"] = m["dist_mm"] / m["dur_s"]
    # Signed, along the direction of travel. |step| hides retrograde motion
    # by counting it as progress, and that is the whole story for C and D.
    sd = d[a + 2000:b - 2000].astype(float)
    # The sign comes from the net travel of the whole run, never from one
    # sample: with the rotor ringing, the first non-zero sample is often a
    # reversal, which would invert the run and report 82% "backward".
    sgn = np.sign(float(pos[b] - pos[a])) or 1.0
    fwd = sd * sgn
    m["back_pct"] = 100.0 * float(np.mean(fwd < 0))
    m["back_mm"] = float(fwd[fwd < 0].sum()) / 1e6
    m["net_mm_s"] = float(fwd.mean()) / 1e3
    m["fwd"] = fwd
    moving = seg[seg > 0.5 * seg.max()]
    m["burst_mm_s"] = float(np.median(moving)) / 1e3
    m["clip_pct"] = 100.0 * float(np.mean(seg == CEILING_NM))
    # ramp: time from motion start to first reaching 90% of the burst rate
    head = np.abs(d[a:a + 4000])
    thr = 0.9 * np.median(moving)
    m["ramp_ms"] = int(np.argmax(head >= thr)) if (head >= thr).any() else -1
    m["slice"] = (a, b)
    return m


def group(runs, pos):
    """Four groups, in run order, three clean transfers kept from each.\n
    Each profile contributes four long runs: three plain transfers and a
    fourth that the home approach is appended to, which runs much longer and
    would drag the plateau statistics down. Keeping the three shortest in
    each group drops it without having to know where homing began.
    """
    rows = [measure(pos, a, b) for a, b in runs]
    rows = [r for r in rows if r["dur_s"] > 5.0]
    per = len(rows) // len(PROFILES)
    out = {}
    for k, (name, v, acc) in enumerate(PROFILES):
        block = rows[k * per:(k + 1) * per]
        sel = sorted(block, key=lambda r: r["dur_s"])[:3]
        out[name] = {"commanded_mm_s": v, "commanded_mm_s2": acc,
                     "runs": sel,
                     "dur_s": float(np.mean([r["dur_s"] for r in sel])),
                     "mean_mm_s": float(np.mean([r["mean_mm_s"] for r in sel])),
                     "burst_mm_s": float(np.mean([r["burst_mm_s"] for r in sel])),
                     "clip_pct": float(np.mean([r["clip_pct"] for r in sel])),
                     "ramp_ms": float(np.mean([r["ramp_ms"] for r in sel])),
                     "dist_mm": float(np.mean([r["dist_mm"] for r in sel])),
                     "back_pct": float(np.mean([r["back_pct"] for r in sel])),
                     "back_mm": float(np.mean([r["back_mm"] for r in sel])),
                     "net_mm_s": float(np.mean([r["net_mm_s"] for r in sel]))}
    return out


# ------------------------------------------------------------------ rendering

def render(pos, g, rows_all):
    style()
    d = np.diff(pos)
    fig = plt.figure(figsize=(13.4, 12.6))
    ax = fig.add_axes([0.075, 0.565, 0.880, 0.155])   # recorded velocity
    bx = fig.add_axes([0.075, 0.340, 0.400, 0.170])   # distance vs clipping
    cx = fig.add_axes([0.555, 0.340, 0.400, 0.170])   # commanded vs achieved
    dx = fig.add_axes([0.075, 0.055, 0.400, 0.185])   # step histogram
    ex = fig.add_axes([0.555, 0.055, 0.400, 0.185])   # 120 ms detail

    titleblock3(
        fig,
        "Where the 2.500 mm/s ceiling comes from",
        "One 59.338 mm transfer repeated under four commanded profiles "
        "(IsoTest.csv, open-loop diagnostic). A->B changes acceleration "
        "alone, A->C changes speed alone.",
        f"2.500 mm/s caps the FIRMWARE's pulse rate, not the stage or the "
        f"interferometer. A and B run at exactly 4 microsteps per 1 ms, 100% "
        f"forward, no jitter. Commanding below it does not avoid it: C and D "
        f"burst at the same rate but stall every 20 ms for the telemetry "
        f"write, reversing on {g['C']['back_pct']:.0f}-"
        f"{g['D']['back_pct']:.0f}% of samples.",
        width=152)

    # --- (a) recorded velocity, one representative transfer per profile ---
    for name, v, acc in PROFILES:
        r = g[name]["runs"][0]
        a, b = r["slice"]
        seg = d[a:b].astype(float)
        sm = np.abs(np.convolve(seg, np.ones(101) / 101, mode="same")) / 1e3
        t = np.arange(len(sm)) * DT
        ax.plot(t, sm, color=COLOUR[name], lw=1.6,
                label=f"{name}  commanded {v:g} mm/s, {acc:g} mm/s² "
                      f"→ recorded {g[name]['mean_mm_s']:.3f} mm/s")
    ax.axhline(CEILING_NM / 1e3, color=INK, lw=1.1, ls=(0, (5, 3)))
    ax.text(0.998, CEILING_NM / 1e3 - 0.19,
            f"ceiling {CEILING_NM / 1e3:.3f} mm/s = 4 microsteps per 1 ms",
            transform=ax.get_yaxis_transform(), ha="right", size=9, color=INK)
    ax.set_xlabel("time from the start of the transfer (s)")
    ax.set_ylabel("recorded speed (mm/s)\n100 ms mean")
    ax.set_xlim(0, 80)
    ax.set_ylim(0, 3.1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.10), ncol=2,
              fontsize=8.8, handlelength=1.8)
    ax.text(0.004, 0.97, "(a)  A and B sit on the ceiling together although "
            "their accelerations differ 4x", size=9.6, color=INK2,
            transform=ax.transAxes, va="top")

    # --- (b) distance closes regardless of how hard it clips --------------
    for name, _, _ in PROFILES:
        xs = [r["clip_pct"] for r in g[name]["runs"]]
        ys = [r["dist_mm"] for r in g[name]["runs"]]
        bx.plot(xs, ys, "o", ms=7, mfc=COLOUR[name], mec=SURFACE, mew=1.0,
                label=name)
    dists = [r["dist_mm"] for name, _, _ in PROFILES for r in g[name]["runs"]]
    bx.axhline(float(np.mean(dists)), color=INK3, lw=0.9, ls=(0, (4, 3)))
    bx.set_xlabel("samples at exactly 2500 nm (%)")
    bx.set_ylabel("transfer distance (mm)")
    bx.set_ylim(np.mean(dists) - 0.004, np.mean(dists) + 0.004)
    bx.set_title("(b)  not the interferometer", size=10.5, loc="left")
    bx.text(0.03, 0.10,
            f"all 12 transfers close within "
            f"{1000 * (max(dists) - min(dists)):.1f} µm while clipping runs "
            f"from {min(dists and [r['clip_pct'] for n, _, _ in PROFILES for r in g[n]['runs']]):.0f}% "
            f"to {max([r['clip_pct'] for n, _, _ in PROFILES for r in g[n]['runs']]):.0f}%.\n"
            f"A counter dropping counts would record short.",
            transform=bx.transAxes, size=8.6, color=INK2)
    bx.legend(loc="lower right", fontsize=8.4, ncol=4, handletextpad=0.2,
              columnspacing=0.8)

    # --- (c) where the retrograde motion is -------------------------------
    names = [p[0] for p in PROFILES]
    xs = np.arange(len(names))
    back = [g[n]["back_pct"] for n in names]
    cx.bar(xs, back, width=0.55, color=[COLOUR[n] for n in names], lw=0)
    for i, n in enumerate(names):
        cx.text(i, back[i] + 0.6,
                f"{back[i]:.1f}%" + ("" if back[i] > 0.05 else "  none"),
                ha="center", size=8.8, color=INK)
        if g[n]["back_mm"] < -0.05:
            cx.text(i, back[i] + 2.6, f"{g[n]['back_mm']:+.1f} mm",
                    ha="center", size=8.0, color=INK2)
    cx.set_xticks(xs)
    cx.set_xticklabels([f"{n}\n{g[n]['commanded_mm_s']:g} mm/s" for n in names],
                       size=8.8)
    cx.set_ylabel("samples moving backward (%)")
    cx.set_ylim(0, 24)
    cx.set_title("(c)  commanding slower makes it worse", size=10.5, loc="left")
    cx.text(0.03, 0.80,
            "A and B never reverse. C and D give back\n"
            "7 mm of a 53 mm transfer and take it again.",
            transform=cx.transAxes, size=8.6, color=INK2)

    # --- (d) the step distribution, in microstep units --------------------
    steady = []
    for n in names:
        for r in g[n]["runs"]:
            a, b = r["slice"]
            steady.append(np.abs(d[a + 2000:b - 2000]))
    steady = np.concatenate(steady)
    dx.hist(steady, bins=np.arange(0, CEILING_NM + 25, 25), color=INK3, lw=0)
    for k in range(1, 5):
        dx.axvline(k * NM_PER_USTEP, color=INK, lw=0.8, ls=(0, (2, 2)),
                   alpha=0.55)
        dx.text(k * NM_PER_USTEP, 2.6e5, f"{k}", ha="center", size=8.2,
                color=INK)
    dx.set_yscale("log")
    dx.set_xlabel("|step| per 1 ms sample (nm)")
    dx.set_ylabel("samples")
    dx.set_xlim(0, CEILING_NM + 260)
    dx.set_ylim(top=6e5)
    dx.set_title("(d)  microsteps per sample", size=10.5, loc="left")
    dx.annotate(f"{int(np.sum(steady == CEILING_NM)):,} samples at exactly "
                f"{CEILING_NM:.0f} nm,\n0 beyond it",
                xy=(CEILING_NM, np.sum(steady == CEILING_NM)),
                xytext=(-18, -64), textcoords="offset points", ha="right",
                size=8.6, color="#b3261e",
                arrowprops=dict(arrowstyle="-", color="#b3261e", lw=0.9))

    # --- (e) signed detail of profile C, against the telemetry period -----
    a, b = g["C"]["runs"][0]["slice"]
    fwd = g["C"]["runs"][0]["fwd"]
    z = slice(10000, 10120)
    tz = np.arange(120)
    ex.axhline(0, color=INK3, lw=0.9)
    ex.axhline(CEILING_NM, color=INK, lw=1.1, ls=(0, (5, 3)))
    ex.plot(tz, fwd[z], color=COLOUR["C"], lw=1.2, marker="o", ms=2.4)
    ex.fill_between(tz, 0, np.minimum(fwd[z], 0), color="#b3261e", alpha=0.30,
                    lw=0)
    for k in range(0, 7):
        t0 = 20 * k + float(np.argmin(fwd[z][:20]))
        if t0 <= 118:
            ex.axvline(t0, color=INK3, lw=0.8, ls=(0, (2, 3)))
    ex.set_xlabel("time (ms)")
    ex.set_ylabel("signed step per 1 ms (nm)")
    ex.set_ylim(-CEILING_NM * 1.15, CEILING_NM * 1.25)
    ex.set_title("(e)  profile C, signed - reversals every 20 ms",
                 size=10.5, loc="left")
    ex.text(0.03, 0.05,
            "dotted lines are 20 ms apart: SAMPLE_INTERVAL_MS.\n"
            "The pulse train stalls for the telemetry write and the rotor rings.",
            transform=ex.transAxes, size=8.4, color=INK2)

    for x in (ax, bx, cx, dx, ex):
        x.grid(True, axis="y")
        x.grid(False, axis="x")

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "iso_velocity_ceiling_explained.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    idx, val = load()
    pos = unwrap_u32(val)
    runs = find_transfers(pos)
    g = group(runs, pos)

    print(f"{SRC.name}: {len(idx):,} samples, {len(idx) * DT:.1f} s, "
          f"{len(runs)} transfers over 50 mm\n")
    print(f"{'':2} {'cmd v':>6} {'cmd a':>6} {'dist':>8} {'dur':>7} "
          f"{'mean v':>7} {'burst':>7} {'clip%':>6} {'ramp ms':>8}")
    for name, v, acc in PROFILES:
        p = g[name]
        print(f"{name:2} {v:6.1f} {acc:6.1f} {p['dist_mm']:8.3f} "
              f"{p['dur_s']:7.2f} {p['mean_mm_s']:7.3f} "
              f"{p['burst_mm_s']:7.3f} {p['clip_pct']:6.1f} "
              f"{p['ramp_ms']:8.0f}")
    print(f"\n1 microstep = {NM_PER_USTEP:.0f} nm; the ceiling is 4 of them "
          f"per 1 ms = {CEILING_NM:.0f} nm = 2.500 mm/s")
    print(f"= {4 / DT / 1000:.0f} kpulses/s at 1/{MRES}, i.e. "
          f"{4 / DT / MRES:.0f} full steps/s")
    print(f"\nwrote {render(pos, g, runs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
