#!/usr/bin/env python3
"""ISO 230-2:2014 evaluation of an ISO run, from the log and the scope CSV.

The board log supplies what was commanded and when each 0.5 s record window
opened and closed; the scope CSV supplies what the interferometer measured.
The two clocks are tied together by fitting scope = rate*board + offset on the
arrival instants, which lands every stop to within a few tens of ms - far
inside the 500 ms window.

One constant datum ties the counter's absolute value to machine coordinates.
Nothing is gain-fitted: a scale error is a real systematic and belongs in E,
not in the fit. Every quantity the standard defines is a difference or a
spread, so the datum choice cannot change any of them.

    python analysis/iso230_2_evaluate.py \\
        --log data/terminal_20261021_183115.log \\
        --csv ../ISO230-2_v3.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

U256 = 25600.0
HEADER_ROWS = 23
WRAP = 2 ** 32


def load_stops(path):
    stops, targets = [], {}
    for line in Path(path).open(encoding="utf-8", errors="replace"):
        if "# STOP," in line:
            f = line.split("# STOP,")[1].strip().split(",")
            stops.append(dict(state=f[0], run=int(f[1]), dir=int(f[2]),
                              i=int(f[3]), cmd_mm=float(f[5]),
                              t_arr=int(f[6]), t0=int(f[7]), t1=int(f[8])))
        elif "# TARGET," in line:
            f = line.split("# TARGET,")[1].strip().split(",")
            targets[int(f[0])] = float(f[2])
    return stops, targets


def load_scope(path):
    val = []
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i < HEADER_ROWS:
                continue
            line = line.strip()
            if not line or line.startswith("EOF"):
                continue
            q = line.split("\t")
            if len(q) == 2:
                val.append(int(q[1]))
    val = np.array(val, np.int64)
    d = np.diff(val)
    d[d > WRAP // 2] -= WRAP
    d[d < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(d)
    return pos / 1e6                      # mm


def fit_clock(stops, mx):
    """scope_ms = rate*board_ms + offset, fitted on arrival instants."""
    t_arr = np.array([s["t_arr"] for s in stops], float)
    d = np.diff(mx)
    sp = np.convolve(np.abs(d), np.ones(51) / 51, mode="same") * 1000
    moving = sp > 0.02
    ends = np.where(moving[:-1] & ~moving[1:])[0] + 1

    def cost(rate, off):
        pred = t_arr * rate + off
        k = np.clip(np.searchsorted(ends, pred), 1, len(ends) - 1)
        return np.median(np.minimum(np.abs(ends[k] - pred),
                                    np.abs(ends[k - 1] - pred)))

    best = (1e18, 1.0, 0.0)
    for off in np.arange(-5000, 120000, 250.0):
        c = cost(1.0, off)
        if c < best[0]:
            best = (c, 1.0, off)
    for rate in np.arange(0.99990, 1.00011, 0.00001):
        for off in np.arange(best[2] - 2000, best[2] + 2000, 25.0):
            c = cost(rate, off)
            if c < best[0]:
                best = (c, rate, off)
    return best[1], best[2], best[0]


def evaluate(stops, mx, rate, off):
    iso = [s for s in stops if s["state"] == "ISO"]
    meas = []
    for s in iso:
        a = int(round(s["t0"] * rate + off))
        b = int(round(s["t1"] * rate + off))
        a, b = max(0, a), min(len(mx), b)
        meas.append(mx[a:b].mean())
    meas = np.array(meas)
    cmd = np.array([s["cmd_mm"] for s in iso])
    datum = float(np.mean(meas - cmd))        # one constant, no gain fit
    dev_um = (meas - datum - cmd) * 1000.0

    idx = np.array([s["i"] for s in iso])
    dirs = np.array([s["dir"] for s in iso])
    out = {}
    for i in sorted(set(int(v) for v in idx)):
        row = {}
        for dkey, dval in (("up", +1), ("dn", -1)):
            m = (idx == i) & (dirs == dval)
            v = dev_um[m]
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
    return out, dev_um, idx, dirs, datum


def summarise(per):
    up = np.array([p["xbar_up"] for p in per.values()])
    dn = np.array([p["xbar_dn"] for p in per.values()])
    s_up = np.array([p["s_up"] for p in per.values()])
    s_dn = np.array([p["s_dn"] for p in per.values()])
    mid = np.array([p["xbar"] for p in per.values()])
    B = np.array([p["B"] for p in per.values()])
    hi = np.maximum(up + 2 * s_up, dn + 2 * s_dn)
    lo = np.minimum(up - 2 * s_up, dn - 2 * s_dn)
    return {
        "B":    float(np.max(np.abs(B))),
        "B_mean": float(np.mean(B)),
        "R_up": float(np.max([p["R_up"] for p in per.values()])),
        "R_dn": float(np.max([p["R_dn"] for p in per.values()])),
        "R":    float(np.max([p["R"] for p in per.values()])),
        "E_up": float(up.max() - up.min()),
        "E_dn": float(dn.max() - dn.min()),
        "E":    float(max(up.max(), dn.max()) - min(up.min(), dn.min())),
        "M":    float(mid.max() - mid.min()),
        "A":    float(hi.max() - lo.min()),
        "s_med": float(np.median(np.concatenate([s_up, s_dn]))),
        "s_max": float(np.max(np.concatenate([s_up, s_dn]))),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--log", default=here.parent /
                    "data" / "terminal_20261021_183115.log")
    ap.add_argument("--csv", default=here.parent.parent / "ISO230-2_v3.csv")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    stops, targets = load_stops(a.log)
    mx = load_scope(a.csv)
    rate, off, resid = fit_clock(stops, mx)
    per, dev, idx, dirs, datum = evaluate(stops, mx, rate, off)
    S = summarise(per)

    n_iso = sum(1 for s in stops if s["state"] == "ISO")
    print(f"log {Path(a.log).name}   csv {Path(a.csv).name}")
    print(f"clock: scope = {rate:.6f}*board + {off:.0f} ms "
          f"(median arrival residual {resid:.1f} ms)")
    print(f"datum: {datum:.6f} mm   stops: {n_iso} = "
          f"{len(per)} targets x 2 directions x "
          f"{per[min(per)]['n_up']} runs\n")

    print(f"{'i':>3} {'target mm':>10} {'xbar_up':>9} {'xbar_dn':>9} "
          f"{'s_up':>7} {'s_dn':>7} {'B_i':>8} {'R_i':>7}")
    for i, p in per.items():
        print(f"{i:3d} {p['cmd_mm']:10.4f} {p['xbar_up']:9.3f} "
              f"{p['xbar_dn']:9.3f} {p['s_up']:7.3f} {p['s_dn']:7.3f} "
              f"{p['B']:+8.3f} {p['R']:7.3f}")

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
            {"clock": {"rate": rate, "offset_ms": off, "resid_ms": resid},
             "datum_mm": datum, "per_target": per, "summary": S}, indent=2))
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
