#!/usr/bin/env python3
"""Re-derive E1's reversal with the board log removed from the chain.

e1_run.py puts the capture into stage coordinates by fitting the ESP's clock
rate and offset to its own logged move intervals, then averages each stop over
the window the firmware says it recorded. That is a reasonable thing to do and
it is also the one step a sceptic should want removed, because it is the step
that failed in the ISO 230-2 v3 evaluation.

So this script never opens the log. It finds the dwells in the trace, matches
them to the planned stop order by position alone, and averages inside each
measured hold. Several window choices are reported, including ones that use no
nominal settle time at all, because a number that moves with the window is a
windowing artefact and a number that does not is the stage.
"""
from __future__ import annotations

from pathlib import Path
import csv
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_E2 = _HERE.parent / "e2_global_map"
_ISO = _HERE.parents[2] / "iso230-2" / "analysis"
for p in (str(_E2), str(_ISO)):
    if p not in sys.path:
        sys.path.insert(0, p)
from render_raw_measurement import read_scope            # noqa: E402
from e2_run import find_holds, match                     # noqa: E402

SCOPE = _HERE / "data" / "E1_Experimental_Run.csv"
PLAN = _HERE / "sequence" / "e1_measurement_plan.csv"
EDGE_MARGIN_MS = 50


def load_plan():
    rows = [r for r in csv.DictReader(PLAN.open()) if r["action"] == "stop"]
    return (np.array([r["state"] for r in rows]),
            np.array([r["window"] for r in rows]),
            np.array([r["direction"] for r in rows]),
            np.array([int(r["index"] or 0) for r in rows]),
            np.array([float(r["target_mm"]) for r in rows]))


def windows(a, b):
    """Window edge pairs, per named recipe, from the measured hold alone."""
    hi_cap = b - EDGE_MARGIN_MS
    return {
        "firmware window, measured edges": (a + 2000, np.minimum(a + 2500, hi_cap)),
        "last 500 ms of the hold":         (hi_cap - 500, hi_cap),
        "middle 50 %":                     (a + (b - a) // 4, b - (b - a) // 4),
        "middle 20 %":                     (a + 2 * (b - a) // 5, b - 2 * (b - a) // 5),
        "second half":                     (a + (b - a) // 2, hi_cap),
        "whole hold, 50 ms trimmed":       (a + 50, hi_cap),
    }


def reversal(hold_mm, target_mm, direction, keep):
    """Mean B_i over the targets measured in both directions."""
    up = direction == "up"
    down = direction == "down"
    bi, per_window = [], {}
    for t in np.unique(target_mm[keep]):
        at = keep & (target_mm == t)
        u, d = hold_mm[at & up], hold_mm[at & down]
        if u.size and d.size and np.isfinite(u).all() and np.isfinite(d).all():
            bi.append(u.mean() - d.mean())
    return np.array(bi) * 1000.0, per_window


def main():
    scope = read_scope(SCOPE)
    d = scope.displacement_mm
    state, window, direction, index, target = load_plan()

    a, b = find_holds(d)
    hold_mm = np.array([d[x + (y - x) // 4: y - (y - x) // 4].mean()
                        for x, y in zip(a, b)])
    assign, datum = match(hold_mm, target)
    pos = d - datum
    ok = assign >= 0
    print(f"{d.size} samples, {scope.time_s[-1] / 60:.2f} min")
    print(f"holds found {a.size}, matched {ok.sum()}/{target.size}, "
          f"datum {datum:.6f} mm")

    hs = np.where(ok, a[np.clip(assign, 0, None)], -1)
    he = np.where(ok, b[np.clip(assign, 0, None)], -1)
    dwell = (he - hs)[ok]
    print(f"measured dwell ms: median {np.median(dwell):.0f}  "
          f"5-95% {np.percentile(dwell, 5):.0f}-{np.percentile(dwell, 95):.0f}")

    recipes = windows(hs, he)
    print()
    hdr = f"{'window recipe':<32}" + "".join(f"{w:>12}" for w in ("W1", "W2", "W3", "all"))
    print(hdr)
    print("-" * len(hdr))
    for name, (lo, hi) in recipes.items():
        h = np.full(target.size, np.nan)
        for i in np.flatnonzero(ok):
            if hi[i] - lo[i] >= 100:
                h[i] = pos[lo[i]:hi[i]].mean()
        row = f"{name:<32}"
        for w in ("W1", "W2", "W3"):
            keep = ok & (window == w) & np.isfinite(h)
            bi, _ = reversal(h, target, direction, keep)
            row += f"{bi.mean():>+12.3f}"
        keep = ok & np.isin(window, ("W1", "W2", "W3")) & np.isfinite(h)
        bi, _ = reversal(h, target, direction, keep)
        row += f"{bi.mean():>+12.3f}"
        print(row)

    # Spread and trend for the headline recipe, so the mean is not quoted bare.
    lo, hi = recipes["last 500 ms of the hold"]
    h = np.full(target.size, np.nan)
    for i in np.flatnonzero(ok):
        h[i] = pos[lo[i]:hi[i]].mean()
    keep = ok & np.isin(window, ("W1", "W2", "W3")) & np.isfinite(h)
    tg = np.unique(target[keep])
    bi, _ = reversal(h, target, direction, keep)
    slope = np.polyfit(tg, bi, 1)[0] * 1000.0
    print()
    print(f"last 500 ms: n={bi.size} targets, mean {bi.mean():+.3f} um, "
          f"sd {bi.std(ddof=1):.3f}, max|Bi| {np.abs(bi).max():.3f}, "
          f"slope {slope:+.2f} nm/mm")

    # Does the window edge matter at all? Sweep where the 500 ms sits.
    print()
    print("500 ms window slid through the hold (offset from hold start):")
    for off in (100, 400, 700, 1000, 1300, 1600, 1900):
        h = np.full(target.size, np.nan)
        for i in np.flatnonzero(ok):
            x, y = hs[i] + off, hs[i] + off + 500
            if y <= he[i] - EDGE_MARGIN_MS:
                h[i] = pos[x:y].mean()
        keep = ok & np.isin(window, ("W1", "W2", "W3")) & np.isfinite(h)
        bi, _ = reversal(h, target, direction, keep)
        if bi.size:
            print(f"  +{off:>4d} ms  n={bi.size:>3d} targets  "
                  f"B {bi.mean():+.3f} um")


if __name__ == "__main__":
    main()
