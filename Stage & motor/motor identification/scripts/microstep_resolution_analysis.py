#!/usr/bin/env python3
"""Analyse a microstep-resolution recording: how finely does the motor step?

Reads the CSV written by run_campaign.py for sequence/microstep_resolution_plan.csv
and, for every staircase and +-1 action, finds where the rotor settled after
each commanded microstep (median of the last SETTLE_MS of each hold). Reports
per MRES:

  step       mean measured microstep vs the ideal, in encoder counts and
             arc-seconds, and its spread (std as % of ideal)
  dead       steps that moved the rotor by less than half the ideal step
  err p-p    peak-to-peak position error over the staircase (settled angle
             minus commanded angle, start offset removed)
  hyst       mean forward-minus-back difference at the same commanded angles
  pm1        mean settled response to one microstep, as % of the ideal,
             and the largest transient overshoot

Writes <recording>_steps.csv (every settled step) and
<recording>_summary.csv (one row per MRES) beside the recording.

    python scripts/microstep_resolution_analysis.py recordings/<file>.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics as st
from collections import defaultdict
from pathlib import Path

FULL_STEPS_PER_REV = 200
SETTLE_MS = 50.0


def load(path: Path, cpr: int):
    """Samples per action index, and each action's tag/kind/mres."""
    samples = defaultdict(list)          # index -> [(t_ms, event, cmd_counts, enc)]
    meta = {}
    scale = cpr / (256 * FULL_STEPS_PER_REV)   # u256 -> encoder counts
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("D,"):
                p = line.strip().split(",")
                if len(p) != 10:
                    continue
                try:
                    samples[int(p[3])].append(
                        (int(p[2]) / 1000.0, int(p[4]), int(p[5]) * scale, int(p[6])))
                except ValueError:
                    continue
            elif line.startswith("# ACTION_BEGIN"):
                m = re.search(r"index=(\d+),tag=([^,]+),kind=(\d+),current_ma=(\d+),"
                              r"mres=(\d+)", line)
                if m:
                    meta[int(m.group(1))] = dict(tag=m.group(2), kind=int(m.group(3)),
                                                 mres=int(m.group(5)))
    return samples, meta


def settled(seg):
    """Median encoder reading over the last SETTLE_MS of a hold segment."""
    if not seg:
        return None
    t_end = seg[-1][0]
    # Last SETTLE_MS, or the last half of a shorter hold (2 ms at 1/256).
    window = min(SETTLE_MS, 0.5 * (t_end - seg[0][0]))
    tail = [s[3] for s in seg if s[0] >= t_end - window] or [seg[-1][3]]
    return st.median(tail)


def segments(rows):
    """Split an action's samples at each marked step edge (event != 0).

    Returns [(command_after_edge, samples_until_next_edge)]."""
    out, cur, cmd = [], None, None
    for r in rows:
        if r[1] != 0:
            if cur is not None:
                out.append((cmd, cur))
            cur, cmd = [r], r[2]
        elif cur is not None:
            cur.append(r)
    if cur is not None:
        out.append((cmd, cur))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("recording", type=Path)
    ap.add_argument("--cpr", type=int, default=20000)
    args = ap.parse_args()
    samples, meta = load(args.recording, args.cpr)
    counts_per_fs = args.cpr / FULL_STEPS_PER_REV
    arcsec_per_count = 360 * 3600 / args.cpr

    step_rows, by_mres = [], defaultdict(dict)
    for idx in sorted(meta):
        info, rows = meta[idx], samples.get(idx, [])
        m, kind = info["mres"], info["kind"]
        ideal = counts_per_fs / m
        if kind == 2 and rows:                                  # staircase
            prev_rows = samples.get(idx - 1, [])
            base = settled(prev_rows)
            base_cmd = prev_rows[-1][2] if prev_rows else 0.0
            if base is None:
                continue
            pts = [(base_cmd, base)] + [(c, settled(seg)) for c, seg in segments(rows)]
            direction = "fwd" if pts[-1][0] > pts[0][0] else "back"
            for k in range(1, len(pts)):
                d = pts[k][1] - pts[k - 1][1]
                step_rows.append(dict(mres=m, direction=direction, step=k,
                                      command_counts=round(pts[k][0], 3),
                                      settled_counts=pts[k][1],
                                      delta_counts=d, ideal_counts=ideal))
            by_mres[m][direction] = pts
        elif kind == 3 and rows:                                # +-1 response
            amps, overs = [], []
            segs = segments(rows)
            for k in range(1, len(segs)):
                before, after = settled(segs[k - 1][1]), settled(segs[k][1])
                sign = 1 if segs[k][0] > segs[k - 1][0] else -1
                amps.append(sign * (after - before))
                peak = max(sign * (s[3] - before) for s in segs[k][1])
                overs.append(peak - sign * (after - before))
            if amps:
                by_mres[m]["pm1"] = (st.mean(amps), max(overs))

    summary = []
    print(f"\n{'MRES':>6} {'ideal':>8} {'mean step':>10} {'spread':>7} {'dead':>11} "
          f"{'err p-p':>9} {'hyst':>8} {'+-1 resp':>9} {'overshoot':>9}")
    print(f"{'':>6} {'[cnt]':>8} {'[cnt]':>10} {'[%]':>7} {'steps':>11} "
          f"{'[arcsec]':>9} {'[arcsec]':>8} {'[% ideal]':>9} {'[cnt]':>9}")
    for m in sorted(by_mres):
        d = by_mres[m]
        ideal = counts_per_fs / m
        deltas = [r["delta_counts"] for r in step_rows if r["mres"] == m]
        signed = [r["delta_counts"] * (1 if r["direction"] == "fwd" else -1)
                  for r in step_rows if r["mres"] == m]
        mean_step = st.mean(signed) if signed else float("nan")
        spread = (st.pstdev(signed) / ideal * 100) if len(signed) > 1 else float("nan")
        # Only meaningful where one microstep is at least ~2 encoder counts.
        dead = sum(1 for x in signed if x < 0.5 * ideal) if ideal >= 2 else None
        err_pp = hyst = float("nan")
        if "fwd" in d:
            e = [enc - cmd for cmd, enc in d["fwd"]]
            e = [x - e[0] for x in e]
            err_pp = (max(e) - min(e)) * arcsec_per_count
        if "fwd" in d and "back" in d:
            back = {round(c, 3): enc for c, enc in d["back"]}
            diffs = [enc - back[round(c, 3)] for c, enc in d["fwd"] if round(c, 3) in back]
            if diffs:
                hyst = st.mean(diffs) * arcsec_per_count
        resp, over = d.get("pm1", (float("nan"), float("nan")))
        row = dict(mres=m, ideal_counts=ideal, mean_step_counts=mean_step,
                   step_spread_pct=spread, dead_steps=dead, steps=len(deltas),
                   error_pp_arcsec=err_pp, hysteresis_arcsec=hyst,
                   pm1_response_pct=resp / ideal * 100 if ideal else float("nan"),
                   pm1_overshoot_counts=over)
        summary.append(row)
        dead_txt = f"{dead}/{len(deltas)}" if dead is not None else "n/a"
        print(f"{'1/' + str(m):>6} {ideal:8.2f} {mean_step:10.3f} {spread:7.0f} "
              f"{dead_txt:>11} {err_pp:9.0f} {hyst:8.0f} "
              f"{row['pm1_response_pct']:9.0f} {over:9.1f}")
    print(f"\n1 count = {arcsec_per_count:.1f} arcsec; steps under ~1 count "
          f"(1/128, 1/256) are below the encoder's resolution individually.")

    stem = args.recording.with_suffix("")
    for name, rows in ((f"{stem}_steps.csv", step_rows), (f"{stem}_summary.csv", summary)):
        if rows:
            with open(name, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
            print(f"wrote {name}")


if __name__ == "__main__":
    main()
