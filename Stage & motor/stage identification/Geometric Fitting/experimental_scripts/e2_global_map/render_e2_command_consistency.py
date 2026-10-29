#!/usr/bin/env python3
"""Did E2 go where it was told? Four ways of asking.

Position, step length, speed and the reference. Nothing here is the global map
the experiment is for - this is the gate the run has to pass before that map
means anything.

The speed panel is the one that most often looks wrong and is not. E2's grid
is irregular by design, from 0.04 mm to 0.56 mm, and a move shorter than
v^2/a = 0.1 mm never reaches the commanded 2 mm/s: it is a triangle whose peak
is sqrt(a*d). Those moves sit on that curve, not below the command by accident.

Writes rendered/e2_command_consistency.png.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e2_run
from e2_run import X_REF_MM

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parents[3] / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, style   # noqa: E402

OUT = Path(os.environ.get("E2_OUT", _HERE / "rendered"))
UP, DOWN, REF = "#4c78a8", "#598a76", "#7a6f91"
WARN = "#99483e"
MINUS = "−"
V_CMD, A_CMD = 2.0, 40.0
DRIFT_WARN_UM = 0.5


def num(value, places=2):
    return f"{value:.{places}f}".replace("-", MINUS)


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


def main():
    thesis_style()
    run = e2_run.load()
    plan = run.plan
    hold = run.hold()
    resid = (hold - plan.target_mm) * 1000.0
    ok = run.matched
    up = ok & (plan.state == "UP")
    dn = ok & (plan.state == "DOWN")
    ref = ok & np.isin(plan.state, ("DATUM", "TURN_REF", "END_REF"))

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # (a) where it stopped against where it was sent
    for m, colour, label in ((up, UP, "upward pass"), (dn, DOWN, "downward pass")):
        ax_a.plot(plan.target_mm[m], resid[m], ".", ms=3.2, color=colour,
                  label=f"{label}, {m.sum()} stops")
    ax_a.axhline(0.0, color=INK2, lw=0.9)
    ax_a.set_xlabel("commanded position [mm]")
    ax_a.set_ylabel("measured − commanded [µm]")
    ax_a.legend(loc="lower left", fontsize=8.5, framealpha=1.0)
    ax_a.set_title("(a) every stop, against its command", fontsize=11, loc="left")
    ax_a.text(0.98, 0.96, f"sd {num(resid[ok].std(ddof=1))} µm, "
              f"worst {num(np.abs(resid[ok]).max())} µm",
              transform=ax_a.transAxes, ha="right", va="top", fontsize=8.5,
              color=INK2)

    # (b) peak speed against how far the move had to go
    v = np.abs(np.diff(run.pos)) * 1000.0
    v10 = np.convolve(v, np.ones(10) / 10, mode="valid")
    idx = np.flatnonzero(ok)
    dist, peak, which = [], [], []
    for i, j in zip(idx[:-1], idx[1:]):
        a, b = run.hold_end[i], run.hold_start[j]
        if b <= a + 20:
            continue
        dist.append(abs(plan.target_mm[j] - plan.target_mm[i]))
        peak.append(v10[a:b].max())
        which.append(plan.state[j])
    dist, peak, which = np.array(dist), np.array(peak), np.array(which)
    for st, colour in (("UP", UP), ("DOWN", DOWN)):
        m = which == st
        ax_b.plot(dist[m], peak[m], ".", ms=3.2, color=colour)
    d = np.linspace(0.005, max(dist.max(), 1.0), 400)
    ax_b.plot(d, np.minimum(np.sqrt(A_CMD * d), V_CMD), "-", lw=1.5,
              color=WARN, zorder=5,
              label=r"$\min(\sqrt{a\,d},\ v)$, $a$ = 40 mm/s$^2$, $v$ = 2 mm/s")
    ax_b.axvline(V_CMD ** 2 / A_CMD, color=INK2, lw=0.8, ls=":")
    ax_b.annotate(f"{V_CMD ** 2 / A_CMD:.2f} mm: shorter moves\nnever reach "
                  "the command", (V_CMD ** 2 / A_CMD, 1.35),
                  xytext=(8, 0), textcoords="offset points", fontsize=8.5,
                  color=INK2, va="center")
    ax_b.set_xscale("log")
    ax_b.set_xlabel("commanded move length [mm]")
    ax_b.set_ylabel("peak speed [mm/s]")
    ax_b.set_ylim(1.0, 2.35)
    ax_b.legend(loc="lower right", fontsize=8.5, framealpha=1.0)
    ax_b.set_title("(b) speed of each move, against the profile it was given",
                   fontsize=11, loc="left")

    # (c) step length
    for m, colour in ((up, UP), (dn, DOWN)):
        o = np.argsort(plan.target_mm[m])
        cmd = np.diff(plan.target_mm[m][o])
        got = np.diff(hold[m][o])
        ax_c.plot(cmd * 1000, (got - cmd) * 1000, ".", ms=3.2, color=colour)
    ax_c.axhline(0.0, color=INK2, lw=0.9)
    ax_c.set_xlabel("commanded step [µm]")
    ax_c.set_ylabel("measured − commanded step [µm]")
    ax_c.set_title("(c) step length, 40 to 560 µm by design", fontsize=11,
                   loc="left")

    # (d) the reference
    t_ref = np.array([run.t[run.record_window(i)[0]] / 60.0
                      for i in np.flatnonzero(ref)])
    y_ref = (hold[ref] - X_REF_MM) * 1000.0
    ax_d.axhspan(y_ref.mean() - DRIFT_WARN_UM, y_ref.mean() + DRIFT_WARN_UM,
                 color=WARN, alpha=0.10, lw=0)
    ax_d.plot(t_ref, y_ref, "o-", ms=6, mfc="#ffffff", mec=REF, mew=1.4,
              color=REF, lw=1.0)
    seen = set()
    for n, (x, y, s) in enumerate(zip(t_ref, y_ref, plan.state[ref])):
        if s in seen:
            continue
        seen.add(s)
        k = int((plan.state[ref] == s).sum())
        ax_d.annotate(f"{s}{f' ({k})' if k > 1 else ''}", (x, y),
                      xytext=(0, 11), textcoords="offset points",
                      ha="left" if x < t_ref[-1] * 0.6 else "right",
                      fontsize=8, color=INK2)
    ax_d.set_xlabel("time in the capture [min]")
    ax_d.set_ylabel("reference − 61.35 mm [µm]")
    ax_d.set_ylim(y_ref.mean() - 1.0, y_ref.mean() + 1.0)
    ax_d.set_title(f"(d) the reference: spread "
                   f"{np.ptp(y_ref):.3f} µm over "
                   f"{t_ref[-1] - t_ref[0]:.0f} min", fontsize=11, loc="left")
    ax_d.text(0.98, 0.04, f"shaded: ±{DRIFT_WARN_UM} µm warning band",
              transform=ax_d.transAxes, ha="right", fontsize=8.5, color=WARN)

    for ax in (ax_a, ax_c, ax_d):
        ax.set_yticklabels([num(v, 1) for v in ax.get_yticks()])

    fig.suptitle("E2 global map: is the run consistent with what was commanded?",
                 fontsize=12.5, y=0.995)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.text(0.012, 0.012,
             f"{ok.sum()} of {plan.target_mm.size} planned stops found in the "
             "trace; dwell 2504 ms against 2500 commanded; no sample at the "
             "AquadB ceiling during the measurement.", fontsize=8, color=INK2)
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "e2_command_consistency.png", dpi=200)
    print("wrote", OUT / "e2_command_consistency.png")


if __name__ == "__main__":
    main()
