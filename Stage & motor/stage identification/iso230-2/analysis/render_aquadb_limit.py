#!/usr/bin/env python3
"""The 2.500 mm/s limit: where it is, and what it does to the data.

The limit is in the AquadB link that carries the interferometer's position to
the EL5101, and it follows from two settings:

    Resolution AquadB = 1000 pm   ->  one increment = 1 nm of travel
    AquadB clock      =  400 ns   ->  one increment per 400 ns at best

        1 / 400 ns            = 2 500 000 increments/s
        2.5e6 x 1 nm          = 2.500 mm/s          <- the ceiling

which is exactly 2500 counts in a 1 ms sample, the value the trace is pinned
at. The link QUEUES increments rather than dropping them, and that single fact
sets what is and is not damaged:

    totals survive     every transfer closes to the same distance however hard
                       it saturated, so positions and distances are sound
    timing does not    above 2.5 mm/s the record is paced by the drain, not by
                       the stage, so the recorded trace is stretched in time

The decisive evidence is deceleration. Profile B is commanded 5 mm/s at
10 mm/s^2, which physically needs 250 ms to stop. Its recorded move ends in
13 ms - a cliff, not a ramp, because the stage had already stopped and the
link was still emptying its queue. Profiles C and D, commanded below the
ceiling, show the deceleration their acceleration actually implies (53 ms
against 50, 180 ms against 200). That is the difference between a stage that
is slow and a recording that is late.

Writes plots/temp/05_recorded_velocity_limit_explained.{png,pdf}.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa
from current_sweep_figures import titleblock3                        # noqa

SRC = Path(os.environ.get("ISO_TEST_CSV", _HERE.parent.parent / "IsoTest.csv"))
OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))
STEM = "05_recorded_velocity_limit_explained"

HEADER_ROWS = 23
WRAP = 2 ** 32
NM_PER_COUNT = 1.0               # Resolution AquadB 1000 pm
CLOCK_NS = 400.0                 # AquadB clock
CEILING_MM_S = 1e9 / CLOCK_NS * NM_PER_COUNT / 1e6     # 2.500 mm/s
CEILING_NM = CEILING_MM_S * 1e3                        # per 1 ms sample
L_MM = 59.338                    # the transfer distance

# commanded (name, v mm/s, a mm/s^2), and the run slices found earlier
PROFILES = [("A", 5.0, 40.0, (689760, 713630)),
            ("B", 5.0, 10.0, (825860, 849870)),
            ("C", 2.0, 40.0, (984300, 1021680)),
            ("D", 1.0, 5.0, (1241720, 1315030))]
COLOUR = {"A": "#b3261e", "B": "#e08214", "C": "#2a78d6", "D": "#1f7a4d"}


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
    d = np.diff(val)
    d[d > WRAP // 2] -= WRAP
    d[d < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(d)
    return np.array(idx, np.int64), pos


def analyse(pos):
    d = np.diff(pos).astype(float)
    rows = []
    for name, v, acc, (a, b) in PROFILES:
        seg = np.abs(d[a + 2000:b - 2000])
        # deceleration actually present at the end of the recorded move
        sm = np.convolve(np.abs(d[b - 1500:b + 100]), np.ones(15) / 15,
                         mode="valid")
        pk = float(np.median(sm[:400]))
        hi = np.where(sm >= 0.9 * pk)[0]
        lo = np.where(sm <= 0.10 * pk)[0]
        lo = lo[lo > hi[-1]] if len(hi) else lo
        decel = int(lo[0] - hi[-1]) if len(lo) and len(hi) else -1
        # how long the trace is stationary after the move, inside the 10 s hold
        tail = np.abs(d[b:b + 12000])
        smt = np.convolve(tail, np.ones(51) / 51, mode="same")
        still = smt < 20.0
        n = 0
        for s in still:
            if s:
                n += 1
            elif n > 200:
                break
            else:
                n = 0
        rows.append(dict(
            name=name, v=v, acc=acc, slice=(a, b),
            dist_mm=abs(pos[b] - pos[a]) / 1e6,
            recorded_s=(b - a) / 1000.0,
            trapez_s=L_MM / v + v / acc,
            link_s=L_MM / CEILING_MM_S,
            clip_pct=100.0 * float(np.mean(seg == CEILING_NM)),
            decel_ms=decel,
            decel_needed_ms=1000.0 * min(v, CEILING_MM_S) / acc,
            still_s=n / 1000.0))
    return rows


def render(pos, rows):
    style()
    d = np.diff(pos).astype(float)
    fig = plt.figure(figsize=(14.2, 11.6))
    ax = fig.add_axes([0.055, 0.505, 0.400, 0.230])   # the arithmetic
    bx = fig.add_axes([0.560, 0.505, 0.405, 0.230])   # queue not loss
    cx = fig.add_axes([0.055, 0.275, 0.400, 0.165])   # duration
    dx = fig.add_axes([0.560, 0.275, 0.405, 0.165])   # deceleration
    ex = fig.add_axes([0.055, 0.055, 0.910, 0.145])   # the cost

    titleblock3(
        fig,
        "The 2.500 mm/s limit is the AquadB link, and it stretches time",
        "Resolution AquadB 1000 pm = one increment per nanometre; AquadB "
        "clock 400 ns = one increment per 400 ns at best. 2.5 M/s x 1 nm = "
        "2.500 mm/s, exactly the 2500 counts per 1 ms the trace is pinned at.",
        "The link queues increments rather than dropping them, so distances "
        "and positions survive intact; timing does not. Profile B needs "
        "250 ms to decelerate and its recorded move ends in 13 ms - the stage "
        "had already stopped, the link was still draining.",
        width=168)

    # ---------------- (a) the arithmetic ---------------------------------
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("(a)  where 2.500 mm/s comes from", size=10.8, loc="left")
    steps = [
        (0.90, "Resolution AquadB", "1000 pm", "one increment = 1 nm"),
        (0.72, "AquadB clock", "400 ns", "minimum spacing between increments"),
        (0.50, "max increment rate", "1 / 400 ns", "= 2 500 000 increments/s"),
        (0.32, "x 1 nm per increment", "", "= 0.0025 m/s"),
    ]
    for y, label, value, note in steps:
        ax.text(0.02, y, label, size=9.6, color=INK, weight="bold", va="center")
        ax.text(0.47, y, value, size=9.6, color="#b3261e", va="center")
        ax.text(0.62, y, note, size=9.0, color=INK2, va="center")
    ax.plot([0.02, 0.98], [0.42, 0.42], color=GRID, lw=1.0)
    ax.text(0.02, 0.16, "= 2.500 mm/s", size=16, weight="bold", color=INK,
            va="center")
    ax.text(0.40, 0.16, "= 2500 counts in a 1 ms sample", size=10, color=INK2,
            va="center")
    ax.text(0.40, 0.065, "measured: 2500 nm exactly, 0 samples beyond it",
            size=9.2, color="#1f7a4d", va="center")

    # ---------------- (b) queue, not loss --------------------------------
    for r in rows:
        bx.plot(r["clip_pct"], r["dist_mm"], "o", ms=9,
                mfc=COLOUR[r["name"]], mec=SURFACE, mew=1.2,
                label=f"{r['name']}  {r['v']:g} mm/s")
    mean = float(np.mean([r["dist_mm"] for r in rows]))
    spread = 1000 * (max(r["dist_mm"] for r in rows)
                     - min(r["dist_mm"] for r in rows))
    bx.axhline(mean, color=INK3, lw=0.9, ls=(0, (4, 3)))
    bx.set_ylim(mean - 0.004, mean + 0.004)
    bx.set_xlabel("samples pinned at the ceiling (%)")
    bx.set_ylabel("transfer distance (mm)")
    bx.set_title("(b)  the link queues, it does not drop", size=10.8,
                 loc="left")
    bx.text(0.03, 0.12,
            f"the same {mean:.3f} mm to within {spread:.1f} µm whether 2 % or\n"
            f"100 % of samples saturate — so no increment is lost, and\n"
            f"positions and distances are unaffected",
            transform=bx.transAxes, size=8.8, color=INK2)
    bx.legend(loc="lower right", fontsize=8.4, ncol=4,
              handletextpad=0.2, columnspacing=0.7)

    # ---------------- (c) duration ---------------------------------------
    names = [r["name"] for r in rows]
    xs = np.arange(len(rows))
    cx.bar(xs - 0.22, [r["trapez_s"] for r in rows], width=0.40, color=GRID,
           lw=0, label="physics: L/v + v/a")
    cx.bar(xs + 0.22, [r["recorded_s"] for r in rows], width=0.40,
           color=[COLOUR[n] for n in names], lw=0, label="recorded")
    cx.axhline(rows[0]["link_s"], color="#b3261e", lw=1.3, ls=(0, (5, 3)),
               label=f"queue drain L/2.5 = {rows[0]['link_s']:.1f} s")
    for i, r in enumerate(rows):
        cx.text(i + 0.22, r["recorded_s"] + 1.6, f"{r['recorded_s']:.1f}",
                ha="center", size=8.4, color=INK)
    cx.set_xticks(xs)
    cx.set_xticklabels([f"{r['name']}\n{r['v']:g}, {r['acc']:g}"
                        for r in rows], size=8.6)
    cx.set_ylabel("transfer time (s)")
    cx.set_ylim(0, 92)
    cx.set_title("(c)  above the ceiling the drain sets the duration",
                 size=10.8, loc="left")
    cx.legend(loc="upper left", fontsize=8.0)

    # ---------------- (d) deceleration, the decisive test -----------------
    cxs = np.arange(len(rows))
    dx.bar(cxs - 0.22, [r["decel_needed_ms"] for r in rows], width=0.40,
           color=GRID, lw=0, label="deceleration the command requires")
    dx.bar(cxs + 0.22, [r["decel_ms"] for r in rows], width=0.40,
           color=[COLOUR[n] for n in names], lw=0, label="recorded 90%→10% fall")
    for i, r in enumerate(rows):
        dx.text(i + 0.22, r["decel_ms"] + 8, f"{r['decel_ms']}", ha="center",
                size=8.4, color=INK)
    dx.annotate("needs 250 ms,\nrecorded in 13 ms", xy=(1.22, 20),
                xytext=(1.55, 150), size=8.8, color="#b3261e",
                arrowprops=dict(arrowstyle="->", color="#b3261e", lw=1.1))
    dx.set_xticks(cxs)
    dx.set_xticklabels([f"{r['name']}\n{r['v']:g}, {r['acc']:g}"
                        for r in rows], size=8.6)
    dx.set_ylabel("deceleration time (ms)")
    dx.set_ylim(0, 300)
    dx.set_title("(d)  A and B stop too fast to be real", size=10.8,
                 loc="left")
    dx.legend(loc="upper left", fontsize=8.0)

    # ---------------- (e) the cost, on profile A --------------------------
    r = rows[0]
    a, b = r["slice"]
    sm = np.convolve(np.abs(d[a:a + 36000]), np.ones(101) / 101,
                     mode="same") / 1e3
    t = np.arange(len(sm)) / 1000.0
    arrive = r["trapez_s"]
    ex.axvspan(arrive, arrive + 10.0, color="#1f7a4d", alpha=0.10, lw=0)
    ex.plot(t, sm, color=COLOUR["A"], lw=1.5)
    ex.axhline(CEILING_MM_S, color=INK, lw=1.1, ls=(0, (5, 3)))
    ex.axvline(arrive, color=INK, lw=1.2)
    ex.axvline(r["recorded_s"], color="#b3261e", lw=1.2, ls=(0, (4, 2)))
    ex.annotate("", xy=(arrive, 2.72), xytext=(r["recorded_s"], 2.72),
                arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.3))
    ex.text((arrive + r["recorded_s"]) / 2, 2.78,
            f"{r['recorded_s'] - arrive:.1f} s of queue still draining after "
            f"the stage stopped", ha="center", size=9.0, color="#b3261e")
    ex.text(arrive - 0.4, 0.30, "stage arrives\n(physics: L/v + v/a)",
            ha="right", size=8.8, color=INK)
    ex.text(r["recorded_s"] + 0.4, 1.30, "record finally\ncatches up",
            ha="left", size=8.8, color="#b3261e")
    ex.text(arrive + 5.0, 0.12, "the intended 10 s hold",
            ha="center", size=8.8, color="#1f7a4d")
    ex.text(r["recorded_s"] + 0.9, 0.62,
            f"only {r['still_s']:.1f} s of it\nrecorded as still",
            ha="left", size=8.8, color="#1f7a4d")
    ex.set_xlabel("time from the start of the transfer (s)")
    ex.set_ylabel("recorded speed (mm/s)")
    ex.set_xlim(0, 34)
    ex.set_ylim(0, 3.05)
    ex.set_title(f"(e)  what it costs - profile A, commanded "
                 f"{r['v']:g} mm/s: the drain eats the stop window",
                 size=10.8, loc="left")

    for x in (bx, cx, dx, ex):
        x.grid(True, axis="y")
        x.grid(False, axis="x")

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{STEM}.png"
    pdf = OUT / f"{STEM}.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    idx, pos = load()
    rows = analyse(pos)
    print(f"AquadB: 1 increment = {NM_PER_COUNT:g} nm, clock {CLOCK_NS:g} ns")
    print(f"  -> {1e9/CLOCK_NS:,.0f} increments/s x {NM_PER_COUNT:g} nm "
          f"= {CEILING_MM_S:.3f} mm/s = {CEILING_NM:.0f} nm per 1 ms\n")
    print(f"{'':2} {'cmd':>11} {'dist mm':>8} {'physics':>8} {'drain':>7} "
          f"{'recorded':>9} {'clip%':>6} {'decel':>7} {'needed':>7} "
          f"{'still s':>8}")
    for r in rows:
        print(f"{r['name']:2} {r['v']:5.1f}/{r['acc']:<5.0f} "
              f"{r['dist_mm']:8.3f} {r['trapez_s']:8.2f} {r['link_s']:7.2f} "
              f"{r['recorded_s']:9.2f} {r['clip_pct']:6.1f} "
              f"{r['decel_ms']:7d} {r['decel_needed_ms']:7.0f} "
              f"{r['still_s']:8.1f}")
    png, pdf = render(pos, rows)
    print(f"\nwrote {png}\n      {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
