#!/usr/bin/env python3
"""Why E1's 0.12 mm grid cancels the motor, and what it does instead.

The microstep error is periodic in the motor's electrical cycle: 4 full steps,
64 microsteps at MRES 16, 0.04 mm of stage travel. E1 steps 192 microsteps,
which is exactly three of those cycles, so every target in a window sits at the
same electrical phase and the motor contributes the *same* error to all of
them. It stops being a position-dependent error and becomes one constant offset
per window, which the per-window offset subtraction removes.

The curve here is measured, not assumed: the direction-averaged static error
over one electrical cycle from the motor-identification staircase, one sample
per microstep.

The same spacing deliberately does not divide the 2 mm screw lead. 0.12 mm is
0.06 of a revolution, so the screw phase advances every target and closes only
after 50 of them - the screw is swept while the motor is frozen. That is the
whole design in one sentence, and panel (d) is the picture of it.

Writes rendered/e1_motor_cycle_rejection.png.
"""
from pathlib import Path
import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent.parent / "motor identification"
if str(_MOTOR / "scripts") not in sys.path:
    sys.path.insert(0, str(_MOTOR / "scripts"))
from staircase_figures import INK, INK2, INK3, style   # noqa

CYCLE_CSV = _MOTOR / "analysis" / "staircase" / "fig2_electrical_cycle.csv"
CONFIG = _HERE / "e1_config.json"
OUT = Path(os.environ.get("E1_OUT", _HERE / "rendered"))

GOOD, BAD = "#356449", "#99483e"       # the chosen grid, and the trap
MOTOR, SCREW = "#99483e", "#5f8bc5"
NAIVE_U = 200                          # 0.125 mm: the round number to avoid
UM_PER_FULL_STEP = 10.0                # 2 mm lead / 200 steps
DEG_PER_FULL_STEP = 1.8


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


def load():
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    with CYCLE_CSV.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    # Measured shaft error over one electrical cycle, one sample per microstep,
    # converted to the stage travel it causes: 1.8 deg of shaft = 10 um.
    err = np.array([float(r["mean_deg"]) for r in rows])
    err = (err - err.mean()) / DEG_PER_FULL_STEP * UM_PER_FULL_STEP
    return cfg, err


def sample(err, index):
    """Error at an absolute microstep index, the cycle repeating forever."""
    return err[np.mod(index, err.size).astype(int)]


def main() -> int:
    cfg, err = load()
    up = float(cfg["ustep_per_mm"])          # 1600 microsteps / mm
    cyc = int(cfg["motor_cycle_usteps"])     # 64
    grid = int(cfg["grid_spacing_usteps"])   # 192
    lead = float(cfg["lead_mm"])             # 2.0
    n_win = int(cfg["targets_per_window"]) if "targets_per_window" in cfg else 101
    pp = float(np.ptp(err))

    print(f"electrical cycle {cyc} usteps = {cyc/up*1000:.0f} um; "
          f"grid {grid} usteps = {grid/up*1000:.0f} um = {grid//cyc} cycles")
    print(f"measured motor error {pp:.3f} um peak-to-peak")

    span = 4 * grid                                   # 768 usteps = 0.48 mm
    u = np.arange(span + 1)
    curve = sample(err, u)
    good = np.arange(0, span + 1, grid)
    naive = np.arange(0, span + 1, NAIVE_U)

    thesis_style()
    fig = plt.figure(figsize=(9.4, 8.4))
    ax = fig.add_axes([0.088, 0.740, 0.885, 0.185])
    bx = fig.add_axes([0.088, 0.512, 0.885, 0.185])
    cx = fig.add_axes([0.088, 0.248, 0.885, 0.175])
    dx = fig.add_axes([0.088, 0.060, 0.885, 0.098])

    # --- (a) and (b): the measured cycle, sampled two ways -----------------
    for axis, pts, colour, lab, title in (
            (ax, good, GOOD, f"E1 target, every {grid} µsteps",
             f"(a)  E1 grid: {grid} µsteps = {grid//cyc} × {cyc} = "
             f"{grid//cyc} electrical cycles → every target at one phase"),
            (bx, naive, BAD, f"target, every {NAIVE_U} µsteps",
             f"(b)  A round {NAIVE_U/up*1000:.0f} µm instead: "
             f"{NAIVE_U} µsteps is not a multiple of {cyc} → the phase walks")):
        for k in range(0, span + 1, cyc):
            axis.axvline(k, color=INK3, lw=0.5, alpha=0.55, zorder=0)
        axis.plot(u, curve, color=INK2, lw=1.2, zorder=2,
                  label="measured motor error")
        axis.plot(pts, sample(err, pts), "o", color=colour, ms=8,
                  mec="white", mew=1.3, zorder=3, label=lab)
        axis.axhline(0, color=INK3, lw=0.6, zorder=0)
        axis.set_xlim(0, span)
        axis.set_ylim(-0.72 * pp, 1.18 * pp)
        axis.set_ylabel("motor error (µm)")
        axis.set_title(title, size=9.6, loc="left", pad=6)
        axis.legend(frameon=False, fontsize=8.0, ncol=2, loc="upper right")
        axis.grid(False)
    # One stage-travel axis for the pair, on top of (a); (b) shares the scale.
    top = ax.secondary_xaxis(
        "top", functions=(lambda v: v / up, lambda v: v * up))
    top.set_xlabel("stage travel (mm)", size=9, labelpad=3)
    ax.tick_params(labelbottom=False)
    ax.annotate("", (0, 0.70 * pp), xytext=(cyc, 0.70 * pp),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.1))
    ax.text(cyc + 14, 0.70 * pp, f"one electrical cycle = {cyc} µsteps "
            f"= {cyc/up*1000:.0f} µm", ha="left", va="center", size=8.0,
            color=INK2)
    ax.annotate("", (0, 0.42 * pp), xytext=(grid, 0.42 * pp),
                arrowprops=dict(arrowstyle="<->", color=GOOD, lw=1.3))
    ax.text(grid / 2, 0.46 * pp, f"E1 step = {grid} µsteps = "
            f"{grid/up*1000:.0f} µm", ha="center", va="bottom", size=8.0,
            color=GOOD)
    bx.set_xlabel("microsteps from the start of the window")

    # --- (c) what each target in a real window actually sees ---------------
    idx = np.arange(n_win)
    e_good = sample(err, idx * grid)
    e_naive = sample(err, idx * NAIVE_U)
    cx.axhline(0, color=INK3, lw=0.6, zorder=0)
    cx.plot(idx, e_naive, "o", color=BAD, ms=3.4, alpha=0.85,
            label=f"{NAIVE_U} µsteps: spans {np.ptp(e_naive):.3f} µm, "
                  f"{len(set(np.round(e_naive, 6)))} distinct values")
    cx.plot(idx, e_good, "o", color=GOOD, ms=4.0,
            label=f"{grid} µsteps: constant {e_good[0]:+.3f} µm — "
                  f"one offset, removed per window")
    cx.set_xlim(-1, n_win)
    cx.set_xlabel(f"target index within a window (all {n_win})")
    cx.set_ylabel("motor error (µm)")
    cx.set_title("(c)  The error each target actually samples", size=9.6,
                 loc="left")
    cx.legend(frameon=False, fontsize=8.0, loc="upper right", ncol=1)
    cx.grid(True, axis="y")
    cx.grid(False, axis="x")

    # --- (d) the other half: the screw is swept, not frozen ----------------
    motor_phase = np.mod(idx * grid, cyc) / cyc
    screw_phase = np.mod(idx * grid / up, lead) / lead
    dx.plot(screw_phase, np.full(n_win, 1.0), "|", color=SCREW, ms=13, mew=1.1)
    dx.plot(motor_phase, np.full(n_win, 0.0), "|", color=MOTOR, ms=13, mew=1.6)
    dx.set_yticks([0, 1])
    dx.set_yticklabels([f"motor\n{len(set(motor_phase))} phase",
                       f"screw\n{len(set(np.round(screw_phase, 9)))} phases"],
                       size=8.2)
    dx.set_ylim(-0.6, 1.6)
    dx.set_xlim(-0.02, 1.02)
    dx.set_xlabel("phase within the cycle (motor: 64 µsteps · screw: 2 mm lead)")
    dx.set_title("(d)  The same spacing freezes the motor and sweeps the screw",
                 size=9.6, loc="left")
    dx.grid(True, axis="x")
    dx.grid(False, axis="y")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "e1_motor_cycle_rejection.png")
    plt.close(fig)
    print(f"E1 samples {len(set(motor_phase))} motor phase, "
          f"{len(set(np.round(screw_phase, 9)))} screw phases")
    print(f"wrote {OUT / 'e1_motor_cycle_rejection.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
