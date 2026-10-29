#!/usr/bin/env python3
"""ISO 230-2:2014 evaluation from the capture alone, with no board log.

iso230_2_evaluate.py needs the ESP log: it takes the commanded target, the
direction and - crucially - the record-window edges from it. The v4 run has no
log, and in any case taking window edges from board timestamps is what put the
last tens of milliseconds of every v3 window on the start of the next move.

So this reads only the interferometer trace. Dwells are found in the data, each
one's direction comes from the sign of the move that led into it, and the
target it belongs to is the nearest one on the planned grid. The averaging
window sits inside the measured hold. The ISO quantities are then computed by
the same summarise() as the logged path, so the two are directly comparable.

    python3 iso230_2_evaluate_blind.py --csv ../../ISO230-2_v4.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_E2 = (_HERE.parents[1] / "Geometric Fitting" / "experimental_scripts" /
       "e2_global_map")
if str(_E2) not in sys.path:
    sys.path.insert(0, str(_E2))
import e2_run                                            # noqa: E402
from iso230_2_evaluate import load_scope, summarise      # noqa: E402

SETTLE_MS, RECORD_MS = 2000, 500
EDGE_MARGIN_MS = 50
MATCH_TOL_MM = 1.0        # the ISO grid is 4.891 mm apart


def record_window(a, b):
    """The firmware's 2.0 s settle + 0.5 s record, placed in the measured hold."""
    cap = b - EDGE_MARGIN_MS
    lo = a + SETTLE_MS
    hi = min(a + SETTLE_MS + RECORD_MS, cap)
    if hi - lo < 100:
        hi = cap
        lo = max(a + 200, hi - RECORD_MS)
    return lo, hi


APPROACH_MS = 150


def classify(pos, a, b, target_mm, datum):
    """Per hold: mean position, the target it sits on, and how it was entered.

    Direction is the sign of the last APPROACH_MS of travel before the hold
    begins - the approach direction ISO 230-2 means. Net travel from the
    previous hold will not do: at the two turning points the pass ends and the
    next one starts on the same target, so consecutive holds differ only by the
    stage's own error and the sign is noise.
    """
    hold = np.array([pos[x + (y - x) // 4: y - (y - x) // 4].mean()
                     for x, y in zip(a, b)])
    d = np.abs(hold[:, None] - datum - target_mm[None, :])
    idx = d.argmin(axis=1)
    near = d[np.arange(d.shape[0]), idx] < MATCH_TOL_MM
    approach = pos[a] - pos[np.maximum(a - APPROACH_MS, 0)]
    direction = np.where(np.abs(approach) * 1000.0 < 1.0, 0,
                         np.sign(approach)).astype(int)
    return hold, idx, near, direction


def evaluate(meas, cmd, idx, dirs):
    """ISO 230-2 per-target table, from measured means and commanded targets."""
    datum = float(np.mean(meas - cmd))
    dev_um = (meas - datum - cmd) * 1000.0
    out = {}
    for i in sorted(set(int(v) for v in idx)):
        row = {}
        for dkey, dval in (("up", +1), ("dn", -1)):
            v = dev_um[(idx == i) & (dirs == dval)]
            row[f"n_{dkey}"] = int(v.size)
            row[f"xbar_{dkey}"] = float(v.mean())
            row[f"s_{dkey}"] = float(v.std(ddof=1))
        row["B"] = row["xbar_up"] - row["xbar_dn"]
        row["R_up"] = 4 * row["s_up"]
        row["R_dn"] = 4 * row["s_dn"]
        row["R"] = max(2 * row["s_up"] + 2 * row["s_dn"] + abs(row["B"]),
                       row["R_up"], row["R_dn"])
        row["xbar"] = 0.5 * (row["xbar_up"] + row["xbar_dn"])
        row["cmd_mm"] = float(cmd[idx == i][0])
        out[i] = row
    return out, dev_um, datum


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=_HERE.parents[1] / "ISO230-2_v4.csv")
    ap.add_argument("--targets", default=_HERE.parent / "sequence" /
                    "iso_targets_seed2302_rom122.695.csv")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    target_mm = np.array([float(line.split(",")[1])
                          for line in Path(a.targets).read_text()
                          .splitlines()[1:] if line.strip()])
    pos_raw = load_scope(Path(a.csv))
    hs, he = e2_run.find_holds(pos_raw)

    # Datum from the mode of every hold-minus-target difference: most holds do
    # sit on a target, so the true offset is a narrow spike among the rest.
    hold0 = np.array([pos_raw[x + (y - x) // 4: y - (y - x) // 4].mean()
                      for x, y in zip(hs, he)])
    diff = (hold0[:, None] - target_mm[None, :]).ravel()
    lo, hi = np.percentile(diff, [1, 99])
    counts, edges = np.histogram(diff, bins=int((hi - lo) / 0.002) + 1,
                                 range=(lo, hi))
    datum0 = 0.5 * (edges[counts.argmax()] + edges[counts.argmax() + 1])
    for _ in range(3):
        _, idx, near, _ = classify(pos_raw, hs, he, target_mm, datum0)
        datum0 = float(np.median(hold0[near] - target_mm[idx[near]]))

    hold0, idx, near, dirs = classify(pos_raw, hs, he, target_mm, datum0)
    keep = near & (dirs != 0)
    print(f"csv {Path(a.csv).name}   targets {Path(a.targets).name}")
    print(f"{pos_raw.size} samples, {pos_raw.size / 60000:.2f} min, "
          f"{hs.size} holds found, {keep.sum()} on the grid")
    dwell = (he - hs)[keep]
    print(f"dwell ms: median {np.median(dwell):.0f}  "
          f"5-95% {np.percentile(dwell, 5):.0f}-{np.percentile(dwell, 95):.0f}")

    meas = np.array([pos_raw[slice(*record_window(x, y))].mean()
                     for x, y in zip(hs[keep], he[keep])])
    cmd = target_mm[idx[keep]]
    per, dev_um, datum = evaluate(meas, cmd, idx[keep] + 1, dirs[keep])

    n_up = np.array([p["n_up"] for p in per.values()])
    n_dn = np.array([p["n_dn"] for p in per.values()])
    print(f"datum {datum:.6f} mm   {len(per)} targets x 2 directions x "
          f"{n_up.min()}-{n_up.max()} up / {n_dn.min()}-{n_dn.max()} down\n")

    print(f"{'i':>3} {'target mm':>10} {'xbar_up':>9} {'xbar_dn':>9} "
          f"{'s_up':>7} {'s_dn':>7} {'B_i':>8} {'R_i':>7}")
    for i, p in per.items():
        print(f"{i:3d} {p['cmd_mm']:10.4f} {p['xbar_up']:9.3f} "
              f"{p['xbar_dn']:9.3f} {p['s_up']:7.3f} {p['s_dn']:7.3f} "
              f"{p['B']:+8.3f} {p['R']:7.3f}")

    S = summarise(per)
    print(f"\n{'quantity':34} {'value (um)':>11}")
    for k, label in (("B", "B   reversal value (max |B_i|)"),
                     ("R_up", "R^  unidirectional repeatability, up"),
                     ("R_dn", "Rv  unidirectional repeatability, down"),
                     ("R", "R   bidirectional repeatability"),
                     ("E_up", "E^  systematic deviation, up"),
                     ("E_dn", "Ev  systematic deviation, down"),
                     ("E", "E   systematic deviation, bidirectional"),
                     ("M", "M   mean bidirectional deviation range"),
                     ("A", "A   accuracy")):
        print(f"{label:34} {S[k]:11.3f}")
    print(f"\n  mean B_i (signed) {S['B_mean']:+.3f} um   "
          f"median s_i {S['s_med']:.3f} um   max s_i {S['s_max']:.3f} um")

    if a.json:
        Path(a.json).write_text(json.dumps(
            {"summary": S,
             "per_target": {str(k): v for k, v in per.items()},
             "datum_mm": datum, "csv": Path(a.csv).name,
             "targets": Path(a.targets).name,
             "method": "measured holds, no board log"},
            indent=2))
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
