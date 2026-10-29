#!/usr/bin/env python3
"""Load the E2 run and tie it to the plan, using the capture alone.

The alignment is deliberately independent of the ESP log: holds are found in
the trace, matched to the planned stop sequence by position and order, and the
record window inside each hold is anchored to that hold's own measured start
and end—the firmware's 2.0 s settle and 0.5 s record, but with edges taken
from the data. The serial log is available and is checked separately by the
findings analysis for completion, step phase and reference MSCNT consistency.

That last point is deliberate. The ISO 230-2 v3 evaluation took its window
edges from the board timestamps, which put the last 38 ms of every window on
the start of the next move and biased the reversal by 1.2 um. Anchoring to the
measured hold cannot make that mistake.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_ISO = _HERE.parents[2] / "iso230-2" / "analysis"
if str(_ISO) not in sys.path:
    sys.path.insert(0, str(_ISO))
from render_raw_measurement import read_scope   # noqa: E402

DATA = _HERE / "data"
SCOPE = DATA / "E2_Experimental_Run.csv"
PLAN = _HERE / "sequence" / "e2_measurement_plan.csv"
CACHE = DATA / "e2_run_aligned.npz"

MOVING_NM = 60.0          # smoothed increment that counts as motion
DISTURBANCE_UM = 5.0      # an excursion smaller than this is not a move
MIN_HOLD_MS = 1500        # the commanded dwell is 2500 ms
SETTLE_MS, RECORD_MS = 2000, 500
EDGE_MARGIN_MS = 50
MATCH_TOL_MM = 0.06       # well under the 0.04 mm minimum target spacing... of
                          # the grid, but over the stage's own ~8 um of error
X_REF_MM = 61.35
TURN_LOW_MM, TURN_HIGH_MM = 2.0, 120.695


@dataclass
class Plan:
    order: np.ndarray
    state: np.ndarray
    direction: np.ndarray
    index: np.ndarray
    target_mm: np.ndarray
    phase_bin: np.ndarray


@dataclass
class Run:
    t: np.ndarray             # seconds from the first sample
    pos: np.ndarray           # stage position, mm
    plan: Plan
    hold_start: np.ndarray    # index into t, per matched plan stop (-1 if none)
    hold_end: np.ndarray
    datum_mm: float

    @property
    def matched(self):
        return self.hold_start >= 0

    def record_window(self, i):
        """The firmware's window, placed inside the measured hold."""
        a, b = self.hold_start[i], self.hold_end[i]
        lo = a + SETTLE_MS
        hi = min(a + SETTLE_MS + RECORD_MS, b - EDGE_MARGIN_MS)
        if hi - lo < 100:                       # a short hold: use what is left
            hi = b - EDGE_MARGIN_MS
            lo = max(a + 200, hi - RECORD_MS)
        return lo, hi

    def hold(self):
        out = np.full(self.plan.target_mm.size, np.nan)
        for i in np.flatnonzero(self.matched):
            lo, hi = self.record_window(i)
            out[i] = self.pos[lo:hi].mean()
        return out


def load_plan(path: Path = PLAN) -> Plan:
    rows = [r for r in csv.DictReader(path.open()) if r["action"] == "stop"]
    return Plan(order=np.array([int(r["order"]) for r in rows]),
                state=np.array([r["state"] for r in rows]),
                direction=np.array([1 if r["direction"] == "up" else -1
                                    if r["direction"] == "down" else 0
                                    for r in rows]),
                index=np.array([int(r["index"]) if r["index"] else 0
                                for r in rows]),
                target_mm=np.array([float(r["target_mm"]) for r in rows]),
                phase_bin=np.array([int(r["screw_phase_bin"]) for r in rows]))


def _runs(flag):
    edge = np.diff(flag.astype(np.int8))
    a = np.flatnonzero(edge == 1) + 1
    b = np.flatnonzero(edge == -1) + 1
    if flag[0]:
        a = np.r_[0, a]
    if flag[-1]:
        b = np.r_[b, flag.size]
    return a, b


def find_holds(pos):
    """Stationary stretches, with sub-micron disturbances not counted as moves."""
    step = np.abs(np.diff(pos)) * 1e6
    still = np.convolve(step, np.ones(5) / 5, mode="same") < MOVING_NM
    a, b = _runs(~still)
    # Judge a moving stretch by how far it strays, not by where it ends: the
    # reference visits leave 61.35 for the 58.35 pre-approach and come back,
    # so their net displacement is nil while the excursion is 3 mm.
    span = np.array([np.ptp(pos[x:min(y + 1, pos.size)]) for x, y in zip(a, b)])
    small = span * 1000.0 < DISTURBANCE_UM
    for x, y in zip(a[small], b[small]):
        still[x:y] = True
    a, b = _runs(still)
    keep = (b - a) > MIN_HOLD_MS
    return a[keep], b[keep]


def _align(hold_mm, target_mm, datum, tol=MATCH_TOL_MM, gap=1.0):
    """Sequence alignment of holds to planned stops, allowing gaps on both sides.

    A greedy walk is no good here: one unmatched stop early on and the pointer
    runs away down the list. This is the usual dynamic program - match, skip a
    hold, or skip a stop - so an extra hold from homing or a stop the detector
    merged costs one gap and nothing after it.
    """
    n, m = hold_mm.size, target_mm.size
    cost = np.full((n + 1, m + 1), np.inf)
    cost[0, 0] = 0.0
    cost[1:, 0] = np.arange(1, n + 1) * gap
    cost[0, 1:] = np.arange(1, m + 1) * gap
    back = np.zeros((n + 1, m + 1), np.int8)
    for i in range(1, n + 1):
        d = np.abs(hold_mm[i - 1] - datum - target_mm)
        sub = np.where(d < tol, d, np.inf)
        row, prev = cost[i], cost[i - 1]
        for j in range(1, m + 1):
            c0 = prev[j - 1] + sub[j - 1]
            c1 = prev[j] + gap
            c2 = row[j - 1] + gap
            if c0 <= c1 and c0 <= c2:
                row[j], back[i, j] = c0, 0
            elif c1 <= c2:
                row[j], back[i, j] = c1, 1
            else:
                row[j], back[i, j] = c2, 2
    out = np.full(m, -1)
    i, j = n, m
    while i > 0 and j > 0:
        k = back[i, j]
        if k == 0:
            out[j - 1] = i - 1
            i -= 1
            j -= 1
        elif k == 1:
            i -= 1
        else:
            j -= 1
    return out


def match(hold_mm, target_mm):
    """Datum and hold-to-stop assignment, from positions and order alone.

    The datum is the mode of every hold-minus-target difference: most holds do
    sit on a planned target, so the true offset is a spike a micron wide in a
    sea of unrelated differences.
    """
    diff = (hold_mm[:, None] - target_mm[None, :]).ravel()
    lo, hi = np.percentile(diff, [1, 99])
    counts, edges = np.histogram(diff, bins=int((hi - lo) / 0.001) + 1,
                                 range=(lo, hi))
    datum = 0.5 * (edges[counts.argmax()] + edges[counts.argmax() + 1])
    for _ in range(3):
        a = _align(hold_mm, target_mm, datum)
        ok = a >= 0
        if not ok.any():
            break
        datum = float(np.median(hold_mm[a[ok]] - target_mm[ok]))
    return _align(hold_mm, target_mm, datum), float(datum)


def load(use_cache=True) -> Run:
    plan = load_plan()
    if use_cache and CACHE.exists():
        z = np.load(CACHE)
        return Run(t=z["t"], pos=z["pos"], plan=plan, hold_start=z["hold_start"],
                   hold_end=z["hold_end"], datum_mm=float(z["datum_mm"]))
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
    return Run(t=scope.time_s, pos=pos, plan=plan, hold_start=hs, hold_end=he,
               datum_mm=datum)


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
