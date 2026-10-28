#!/usr/bin/env python3
"""Is the reversal a property of the stage or of the processing?

The three runs are put through one code path that never touches a board clock:
dwells are found in the interferometer trace, matched to the planned stop order
by position, and averaged inside each measured hold. The log contributes only
the ordered list of (target, direction) - plan information - and for E2 there
is no log at all.

Panel (a) varies the averaging window over six recipes. All three runs sit
between -1.60 and -1.72 um and move by less than 0.05 um across the recipes.
The published ISO value of -0.365 um is not reachable by any of them; it needs
the board-timestamped window, whose tail falls on the start of the next move.

Panel (b) slides a 500 ms window through the dwell. The settling tail is worth
about 0.07 um over the whole 2.5 s - more than an order below the 1.2 um the
ISO windowing defect contributed - so where the window sits does not decide B.

Needs the raw captures; results are cached in data/_window_independence.npz.

Writes rendered/reversal_window_independence.{png,pdf}.
"""
from pathlib import Path
import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_STAGE = _HERE.parents[2]
for extra in (_HERE.parent / "e2_global_map",
              _HERE.parents[3] / "motor identification" / "scripts",
              _STAGE / "iso230-2" / "analysis"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))
import e2_run                                              # noqa: E402
from staircase_figures import INK, INK2, style             # noqa: E402
from render_raw_measurement import read_scope              # noqa: E402
from iso230_2_evaluate import load_scope, load_stops       # noqa: E402

OUT = Path(os.environ.get("E1_OUT", _HERE / "rendered"))
CACHE = _HERE / "data" / "_window_independence.npz"
ISO_RESULTS = _STAGE / "iso230-2" / "analysis" / "iso230_2_v3_results.json"
EDGE_MARGIN_MS = 50
MINUS = "−"

RECIPES = ("firmware window,\nmeasured edges", "last 500 ms\nof the hold",
           "middle 50 %", "middle 20 %", "second half",
           "whole hold,\n50 ms trimmed")
SLIDE = (100, 400, 700, 1000, 1300, 1600, 1900)
RUNS = ("E1", "E2", "ISO v3")
COLOURS = {"E1": "#4c78a8", "E2": "#598a76", "ISO v3": "#9d6770"}
PUBLISHED = "#99483e"


def num(value, places=3):
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


def recipe_edges(hs, he):
    """The six averaging windows, every edge taken from the measured hold."""
    cap = he - EDGE_MARGIN_MS
    span = he - hs
    return ((hs + 2000, np.minimum(hs + 2500, cap)),
            (cap - 500, cap),
            (hs + span // 4, he - span // 4),
            (hs + 2 * span // 5, he - 2 * span // 5),
            (hs + span // 2, cap),
            (hs + 50, cap))


def bmean(pos, ok, target, direction, lo, hi):
    """B_i over the targets measured in both directions, for one window."""
    h = np.full(target.size, np.nan)
    for i in np.flatnonzero(ok):
        if hi[i] - lo[i] >= 100:
            h[i] = pos[lo[i]:hi[i]].mean()
    bi = []
    for t in np.unique(target[ok]):
        at = ok & (target == t) & np.isfinite(h)
        u, d = h[at & (direction == 1)], h[at & (direction == -1)]
        if u.size and d.size:
            bi.append(u.mean() - d.mean())
    return np.array(bi) * 1000.0


def blind(pos_raw, target, direction, tol, eligible=None):
    """Holds from the trace, matched to the plan order. No board clock.

    `eligible` marks the plan rows that belong in the reversal: the grid stops
    only. Reference and pilot visits sit on top of grid targets in E1, and
    averaging them into one direction of such a target shifts its B_i.
    """
    e2_run.MATCH_TOL_MM = tol
    a, b = e2_run.find_holds(pos_raw)
    hold_mm = np.array([pos_raw[x + (y - x) // 4: y - (y - x) // 4].mean()
                        for x, y in zip(a, b)])
    assign, datum = e2_run.match(hold_mm, target)
    ok = assign >= 0
    if eligible is not None:
        ok &= eligible
    hs = np.where(ok, a[np.clip(assign, 0, None)], -1)
    he = np.where(ok, b[np.clip(assign, 0, None)], -1)
    pos = pos_raw - datum

    recipes = [bmean(pos, ok, target, direction, lo, hi)
               for lo, hi in recipe_edges(hs, he)]
    slide = []
    for off in SLIDE:
        lo = hs + off
        hi = np.where(lo + 500 <= he - EDGE_MARGIN_MS, lo + 500, lo)
        slide.append(bmean(pos, ok, target, direction, lo, hi))
    return dict(matched=int(ok.sum()), total=int(target.size),
                holds=int(a.size),
                recipe_mean=np.array([x.mean() for x in recipes]),
                recipe_sd=np.array([x.std(ddof=1) for x in recipes]),
                recipe_n=np.array([x.size for x in recipes]),
                slide_mean=np.array([x.mean() for x in slide]))


def e1_inputs():
    rows = [r for r in csv.DictReader(
        (_HERE / "sequence" / "e1_measurement_plan.csv").open())
        if r["action"] == "stop"]
    target = np.array([float(r["target_mm"]) for r in rows])
    direction = np.array([1 if r["direction"] == "up" else -1
                          if r["direction"] == "down" else 0 for r in rows])
    window = np.array([r["window"] for r in rows])
    pos = read_scope(_HERE / "data" / "E1_Experimental_Run.csv").displacement_mm
    return pos, target, direction, 0.06, np.isin(window, ("W1", "W2", "W3"))


def e2_inputs():
    plan = e2_run.load_plan()
    pos = read_scope(_HERE.parent / "e2_global_map" / "data" /
                     "E2_Experimental_Run.csv").displacement_mm
    return (pos, plan.target_mm, plan.direction, 0.06,
            np.isin(plan.state, ("UP", "DOWN")))


def iso_inputs():
    stops, _ = load_stops(_STAGE / "iso230-2" / "data" /
                          "terminal_20261021_183115.log")
    iso = [s for s in stops if s["state"] == "ISO"]
    target = np.array([s["cmd_mm"] for s in iso])
    direction = np.array([1 if s["dir"] > 0 else -1 for s in iso])
    return (load_scope(_STAGE / "ISO230-2_v3.csv"), target, direction, 1.0,
            np.ones(target.size, bool))


def gather(use_cache=True):
    if use_cache and CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        return {k: z[k].item() for k in RUNS}
    out = {}
    for name, fn in zip(RUNS, (e1_inputs, e2_inputs, iso_inputs)):
        pos, target, direction, tol, eligible = fn()
        r = out[name] = blind(pos, target, direction, tol, eligible)
        print(f"{name}: {r['holds']} holds, matched {r['matched']}/"
              f"{r['total']}, B {r['recipe_mean'][1]:+.3f} um over "
              f"{r['recipe_n'][1]} targets, recipe spread "
              f"{np.ptp(r['recipe_mean']):.3f} um")
    np.savez_compressed(
        CACHE, **{k: np.array(v, dtype=object) for k, v in out.items()})
    return out


def published_b():
    s = json.loads(ISO_RESULTS.read_text())
    return float(s.get("summary", s)["B_mean"])


def render(res):
    thesis_style()
    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(11.8, 4.6),
        gridspec_kw=dict(width_ratios=(1.55, 1.0), wspace=0.28))

    y = np.arange(len(RECIPES))[::-1]
    pub = published_b()
    lo_band = min(res[n]["recipe_mean"].min() for n in RUNS) - 0.03
    hi_band = max(res[n]["recipe_mean"].max() for n in RUNS) + 0.03
    ax.axvspan(lo_band, hi_band, color="#eceff3", zorder=0)
    ax.axvline(pub, color=PUBLISHED, lw=1.6, ls=(0, (5, 2)), zorder=2)
    ax.annotate(f"published ISO 230-2 v3,  {num(pub)} µm\n"
                f"(board-timestamped window)",
                xy=(pub, y[-1] - 0.1), xytext=(pub + 0.06, y[-1] - 0.40),
                color=PUBLISHED, fontsize=8.5, ha="left", va="center",
                clip_on=False)
    for k, name in enumerate(RUNS):
        r = res[name]
        ax.errorbar(r["recipe_mean"], y + (k - 1) * 0.22,
                    xerr=r["recipe_sd"] / np.sqrt(r["recipe_n"]),
                    fmt="o", ms=6.5, lw=1.4, capsize=2.5,
                    color=COLOURS[name], mec="#ffffff", mew=0.9,
                    zorder=3,
                    label=f"{name}   {num(r['recipe_mean'][1])} µm"
                          f"   ({r['recipe_n'][1]} targets)")
    ax.set_yticks(y)
    ax.set_yticklabels(RECIPES, fontsize=8.5)
    ax.set_xlim(-2.00, 0.06)
    ax.set_ylim(y[-1] - 0.72, y[0] + 0.52)
    ax.set_xlabel("mean reversal  $\\bar{B}$   (µm)")
    ax.xaxis.set_major_formatter(lambda v, _: num(v, 1))
    ax.annotate(f"every run, every window:\n{num(hi_band - 0.03)} to "
                f"{num(lo_band + 0.03)} µm",
                xy=(hi_band + 0.05, y[3] + 0.5), ha="left", va="center",
                fontsize=8.5, color=INK2)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5,
              handletextpad=0.4, borderaxespad=0.4, labelspacing=0.5)
    ax.set_title("(a)  six averaging windows, every edge from the measured hold",
                 loc="left", fontsize=10, color=INK)

    for name in RUNS:
        bx.plot(SLIDE, res[name]["slide_mean"], "o-", ms=5, lw=1.6,
                color=COLOURS[name], mec="#ffffff", mew=0.9, label=name)
    bx.axhline(pub, color=PUBLISHED, lw=1.6, ls=(0, (5, 2)))
    bx.annotate(f"published,  {num(pub)} µm", xy=(SLIDE[0], pub),
                xytext=(0, 7), textcoords="offset points", ha="left",
                fontsize=8.5, color=PUBLISHED)
    swing = max(np.ptp(res[n]["slide_mean"]) for n in RUNS)
    bx.annotate("settling tail across the whole dwell:\n"
                f"at most {swing:.3f} µm", xy=(0.04, 0.40),
                xycoords="axes fraction", ha="left", va="bottom",
                fontsize=8.5, color=INK2)
    bx.set_xlabel("start of a 500 ms window, from the hold start   (ms)")
    bx.set_ylabel("mean reversal  $\\bar{B}$   (µm)")
    bx.yaxis.set_major_formatter(lambda v, _: num(v, 1))
    bx.set_ylim(-2.00, 0.06)
    bx.legend(loc="upper right", frameon=False, fontsize=8.5,
              borderaxespad=1.6)
    bx.set_title("(b)  where in the dwell the window sits",
                 loc="left", fontsize=10, color=INK)

    fig.suptitle("The reversal does not depend on how it is windowed, and no "
                 "window reproduces the published value",
                 x=0.010, ha="left", fontsize=11.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"reversal_window_independence.{ext}", dpi=220)
    print(f"wrote {OUT / 'reversal_window_independence.png'}")


if __name__ == "__main__":
    render(gather(use_cache="--fresh" not in sys.argv))
