#!/usr/bin/env python3
"""Load the 30 Oct repeat run and tie it to the plan, using the capture alone.

Everything substantive is `e2_run`: the hold detector, the dynamic-program
alignment of holds to planned stops, the mode-of-differences datum, and the
record window anchored to each hold's own measured edges. Only the capture, the
cache and the serial log differ. Processing both sessions through one code path
is the point - a day-to-day comparison is worthless if the two days were
reduced by two implementations.

The plan comes from e2_global_map/sequence. The firmware was re-run unmodified,
so the accepted seed and all 400 targets are the same; see sequence/PROVENANCE.md.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_E2 = _HERE.parent / "e2_global_map"
if str(_E2) not in sys.path:
    sys.path.insert(0, str(_E2))
import e2_run                                      # noqa: E402
from e2_run import (Plan, Run, find_holds, match,   # noqa: E402,F401
                    load_plan, SETTLE_MS, RECORD_MS, EDGE_MARGIN_MS,
                    MATCH_TOL_MM, X_REF_MM)

_ISO = _HERE.parents[2] / "iso230-2" / "analysis"
if str(_ISO) not in sys.path:
    sys.path.insert(0, str(_ISO))
from render_raw_measurement import read_scope       # noqa: E402

DATA = _HERE / "data"
SCOPE = DATA / "E3_Experimental_Run.csv"
LOG = DATA / "e3_terminal_20261030_135726.log"
CACHE = DATA / "e3_run_aligned.npz"
PLAN = _E2 / "sequence" / "e2_measurement_plan.csv"

SESSION = "E3"
SESSION_LABEL = "E3 (repeat of E2 targets, 30 Oct)"


def load(use_cache: bool = True) -> Run:
    """The same Run object e2_run.load() returns, for the 30 Oct capture."""
    plan = load_plan(PLAN)
    if use_cache and CACHE.exists():
        z = np.load(CACHE)
        return Run(t=z["t"], pos=z["pos"], plan=plan,
                   hold_start=z["hold_start"], hold_end=z["hold_end"],
                   datum_mm=float(z["datum_mm"]))
    scope = read_scope(SCOPE)
    d = scope.displacement_mm
    a, b = find_holds(d)
    hold_mm = np.array([d[x + (y - x) // 4: y - (y - x) // 4].mean()
                        for x, y in zip(a, b)])
    assign, datum = match(hold_mm, plan.target_mm)
    pos = d - datum
    hs = np.where(assign >= 0, a[np.clip(assign, 0, None)], -1)
    he = np.where(assign >= 0, b[np.clip(assign, 0, None)], -1)
    np.savez_compressed(CACHE, t=scope.time_s, pos=pos, hold_start=hs,
                        hold_end=he, datum_mm=datum)
    return Run(t=scope.time_s, pos=pos, plan=plan, hold_start=hs,
               hold_end=he, datum_mm=datum)


if __name__ == "__main__":
    run = load(use_cache=False)
    n = run.matched.sum()
    h = run.hold()
    resid = (h - run.plan.target_mm) * 1000.0
    print(f"{run.t.size} samples, {run.t[-1] / 60:.2f} min, datum "
          f"{run.datum_mm:.6f} mm")
    print(f"matched {n}/{run.plan.target_mm.size} planned stops")
    print(f"residual um: median {np.nanmedian(resid):+.3f}  "
          f"sd {np.nanstd(resid):.3f}  |max| {np.nanmax(np.abs(resid)):.3f}")
