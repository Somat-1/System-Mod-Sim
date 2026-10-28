#!/usr/bin/env python3
"""The control for check_e1_no_board_clock.py: the same code on the ISO v3 run.

If a reversal near -1.6 um were an artefact of how E1 is processed here, then
feeding the ISO 230-2 v3 capture through the identical path - dwells found in
the trace, windows anchored to the measured hold - would reproduce the
published -0.365 um. It does not.

The log is read for the ordered list of (target, direction) only. That is plan
information, the same thing e1_measurement_plan.csv is. No board timestamp
reaches any window edge.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_E2 = _HERE.parent / "e2_global_map"
_ISO = _HERE.parents[2] / "iso230-2" / "analysis"
for p in (str(_E2), str(_ISO)):
    if p not in sys.path:
        sys.path.insert(0, p)
import e2_run                                               # noqa: E402
from iso230_2_evaluate import load_stops, load_scope        # noqa: E402

LOG = _HERE.parents[2] / "iso230-2" / "data" / "terminal_20261021_183115.log"
SCOPE = _HERE.parents[2] / "ISO230-2_v3.csv"
EDGE_MARGIN_MS = 50

# The ISO grid is 4.9 mm apart, so the matcher can be given a wide tolerance;
# the dwell is the same 2.0 s settle + 0.5 s record as E1.
e2_run.MATCH_TOL_MM = 1.0


def main():
    stops, _ = load_stops(LOG)
    iso = [s for s in stops if s["state"] == "ISO"]
    target = np.array([s["cmd_mm"] for s in iso])
    direction = np.array(["up" if s["dir"] > 0 else "down" for s in iso])
    pos_raw = load_scope(SCOPE)

    a, b = e2_run.find_holds(pos_raw)
    hold_mm = np.array([pos_raw[x + (y - x) // 4: y - (y - x) // 4].mean()
                        for x, y in zip(a, b)])
    assign, datum = e2_run.match(hold_mm, target)
    pos = pos_raw - datum
    ok = assign >= 0
    print(f"{pos_raw.size} samples, {pos_raw.size / 60000:.2f} min")
    print(f"holds found {a.size}, matched {ok.sum()}/{target.size}, "
          f"datum {datum:.6f} mm")

    hs = np.where(ok, a[np.clip(assign, 0, None)], -1)
    he = np.where(ok, b[np.clip(assign, 0, None)], -1)
    dwell = (he - hs)[ok]
    print(f"measured dwell ms: median {np.median(dwell):.0f}  "
          f"5-95% {np.percentile(dwell, 5):.0f}-{np.percentile(dwell, 95):.0f}")

    hi_cap = he - EDGE_MARGIN_MS
    recipes = {
        "firmware window, measured edges": (hs + 2000, np.minimum(hs + 2500, hi_cap)),
        "last 500 ms of the hold":         (hi_cap - 500, hi_cap),
        "middle 50 %":                     (hs + (he - hs) // 4, he - (he - hs) // 4),
        "middle 20 %":                     (hs + 2 * (he - hs) // 5, he - 2 * (he - hs) // 5),
        "second half":                     (hs + (he - hs) // 2, hi_cap),
        "whole hold, 50 ms trimmed":       (hs + 50, hi_cap),
    }

    print()
    print(f"{'window recipe':<32}{'B mean':>10}{'max|Bi|':>10}{'n':>6}")
    print("-" * 58)
    for name, (lo, hi) in recipes.items():
        h = np.full(target.size, np.nan)
        for i in np.flatnonzero(ok):
            if hi[i] - lo[i] >= 100:
                h[i] = pos[lo[i]:hi[i]].mean()
        bi = []
        for t in np.unique(target[ok]):
            at = ok & (target == t) & np.isfinite(h)
            u, d = h[at & (direction == "up")], h[at & (direction == "down")]
            if u.size and d.size:
                bi.append(u.mean() - d.mean())
        bi = np.array(bi) * 1000.0
        print(f"{name:<32}{bi.mean():>+10.3f}{np.abs(bi).max():>10.3f}"
              f"{bi.size:>6d}")

    print()
    print("published iso230_2_v3_results.json, for comparison:")
    import json
    r = json.loads((_ISO / "iso230_2_v3_results.json").read_text())
    s = r.get("summary", r)
    print(f"{'board-timestamped window':<32}{s['B_mean']:>+10.3f}"
          f"{s['B']:>10.3f}{'24':>6}")


if __name__ == "__main__":
    main()
