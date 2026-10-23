#!/usr/bin/env python3
"""Why the stage runs slow: the step-timing budget, measured.

The board logs, for every move, how many pulses it requested and how long it
took. Dividing gives the pulse period it actually achieved, against the one
the velocity profile asked for. Across 2, 5 and 8 mm/s those points fall on a
straight line,

    T_actual = a * T_nominal + b

and the two coefficients are the two defects, separated:

    a - 1   a fraction of every interval is lost to the telemetry write that
            blocks inside the pulse loop; it scales with the period, so it
            costs the same percentage at any speed
    b       a fixed cost added to every interval, because the scheduler takes
            its deadline from micros() AFTER the pulse rather than advancing
            an absolute one; it only bites once the period approaches it

Writes plots/temp/05_step_timing_shortfall.{png,pdf}.
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
sys.path.insert(0, str(_HERE))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa
from aquadb_bandwidth_check import load_log                           # noqa

OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))
FS_PER_MM = 100.0
MRES = 16
MIN_MM = 5.0
COLOUR = {2.0: "#5f8bc5", 5.0: "#c88e5c", 8.0: "#99483e"}
RED, GREEN = "#99483e", "#356449"


def measure(log):
    """One row per long move: nominal and achieved pulse period, in us."""
    _, moves, _ = load_log(log)
    rows = []
    for m in moves:
        if abs(m["x_to"] - m["x_from"]) < MIN_MM or m["requested"] == 0:
            continue
        rows.append(dict(v=m["v"],
                         t_nom=1e6 / (m["v"] * FS_PER_MM * MRES),
                         t_act=m["elapsed_ms"] * 1000.0 / m["requested"]))
    return rows


def main() -> int:
    log = sys.argv[1] if len(sys.argv) > 1 else str(
        _HERE.parent / "data" / "terminal_20261022_162627.log")
    rows = measure(log)
    if not rows:
        print(f"no long moves found in {log}")
        return 2
    nom = np.array([r["t_nom"] for r in rows])
    act = np.array([r["t_act"] for r in rows])
    vel = np.array([r["v"] for r in rows])
    a, b = np.polyfit(nom, act, 1)
    f = 1.0 - 1.0 / a

    speeds = sorted(set(vel))
    print(f"{Path(log).name}: {len(rows)} long moves")
    print(f"  fit  T_act = {a:.4f} * T_nom + {b:.2f} us")
    print(f"  blocked fraction {100*f:.1f}%   per-pulse overhead {b:.2f} us\n")
    print(f"{'v cmd':>6} {'T_nom':>8} {'T_act':>8} {'v achieved':>11} "
          f"{'shortfall':>10}")
    ach = {}
    for v in speeds:
        m = vel == v
        ta = act[m].mean()
        ach[v] = v * nom[m][0] / ta
        print(f"{v:6.1f} {nom[m][0]:8.1f} {ta:8.1f} {ach[v]:11.3f} "
              f"{100*(1-ach[v]/v):9.1f}%")

    style()
    fig = plt.figure(figsize=(12.4, 4.3))
    ax = fig.add_axes([0.055, 0.165, 0.255, 0.70])
    bx = fig.add_axes([0.395, 0.165, 0.255, 0.70])
    cx = fig.add_axes([0.735, 0.165, 0.245, 0.70])

    # --- (a) what the stage actually does ---------------------------------
    lim = [0, 9]
    ax.plot(lim, lim, color=GREEN, lw=1.3, ls=(0, (5, 3)),
            label="after the fix")
    for v in speeds:
        ax.plot([v], [ach[v]], "o", ms=9, mfc=COLOUR[v], mec=SURFACE, mew=1.2,
                zorder=4)
        ax.annotate(f"−{100*(1-ach[v]/v):.0f}%", (v, ach[v]), xytext=(6, -11),
                    textcoords="offset points", size=8.6, color=COLOUR[v])
    ax.plot(speeds, [ach[v] for v in speeds], color=INK3, lw=1.0, zorder=3,
            label="measured")
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 9)
    ax.set_xlabel("commanded speed (mm/s)")
    ax.set_ylabel("achieved speed (mm/s)")
    ax.set_title("(a)  the stage never reaches the command", size=10,
                 loc="left", color=INK)
    ax.legend(frameon=False, fontsize=8.4, loc="upper left")
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")

    # --- (b) the straight line that separates the two causes --------------
    xs = np.linspace(0, 340, 50)
    bx.plot(xs, a * xs + b, color=INK, lw=1.2, zorder=3)
    bx.plot(xs, xs, color=GREEN, lw=1.2, ls=(0, (5, 3)), zorder=2,
            label="ideal  $T_{act}=T_{nom}$")
    for v in speeds:
        m = vel == v
        bx.plot(nom[m], act[m], "o", ms=8, mfc=COLOUR[v], mec=SURFACE,
                mew=1.1, zorder=4, label=f"{v:g} mm/s")
    bx.annotate(f"slope {a:.3f}\n→ {100*f:.0f}% of every\n     interval blocked",
                (250, a * 250 + b), xytext=(-6, -58),
                textcoords="offset points", ha="right", size=8.2, color=INK2)
    bx.annotate(f"intercept {b:.1f} µs\n→ fixed per-pulse cost",
                (0, b), xytext=(16, 26), textcoords="offset points",
                size=8.2, color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK3, lw=0.8))
    bx.set_xlim(0, 340)
    bx.set_ylim(0, 420)
    bx.set_xlabel("pulse period asked for (µs)")
    bx.set_ylabel("pulse period achieved (µs)")
    bx.set_title("(b)  two defects, separated", size=10, loc="left", color=INK)
    bx.legend(frameon=False, fontsize=8.0, loc="upper left")
    bx.grid(True, axis="y")
    bx.grid(False, axis="x")

    # --- (c) where each microsecond goes ----------------------------------
    idx = np.arange(len(speeds))
    base = np.array([1e6 / (v * FS_PER_MM * MRES) for v in speeds])
    prop = base * (a - 1.0)
    fixed = np.full(len(speeds), b)
    cx.bar(idx, base, width=0.6, color=INK3, lw=0, label="useful (commanded)")
    cx.bar(idx, prop, width=0.6, bottom=base, color="#c88e5c", lw=0,
           label=f"blocked by telemetry ({100*f:.0f}%)")
    cx.bar(idx, fixed, width=0.6, bottom=base + prop, color=RED, lw=0,
           label=f"scheduler overhead ({b:.0f} µs)")
    for i, v in enumerate(speeds):
        tot = base[i] + prop[i] + fixed[i]
        cx.text(i, tot + 9, f"{100*(tot/base[i]-1):.0f}% over",
                ha="center", size=8.4, color=INK)
    cx.set_xticks(idx)
    cx.set_xticklabels([f"{v:g} mm/s" for v in speeds], size=9)
    cx.set_ylabel("pulse period (µs)")
    cx.set_ylim(0, 470)
    cx.set_title("(c)  the period budget", size=10, loc="left", color=INK)
    cx.legend(frameon=False, fontsize=7.8, loc="upper right")
    cx.grid(True, axis="y")
    cx.grid(False, axis="x")

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "05_step_timing_shortfall.png"
    fig.savefig(png)
    fig.savefig(OUT / "05_step_timing_shortfall.pdf")
    plt.close(fig)
    print(f"\nwrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
