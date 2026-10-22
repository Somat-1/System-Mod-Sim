#!/usr/bin/env python3
"""Does the EL5101 keep up? Closure test for the AquadB bandwidth diagnostic.

Three things can go wrong on a move, and the log plus the trace separate all
three:

    firmware did not emit    emitted_pulses < requested_pulses, from the board
    stage did not follow     cannot be seen directly - the axis is open loop
    metrology did not record recorded distance < commanded distance

Only the third is what this test is about. A transmitter that clips and queues
loses nothing, so the recorded distance still matches; a receiver that cannot
decode the edge rate drops counts, so it comes up short and the error grows
with travel. Closure is therefore the criterion, not saturation.

Profiles and the per-move boundaries come from the board log's DIAG records,
so nothing is inferred from the shape of the trace. The board clock is tied to
the scope by fitting scope = rate*board + offset on the move-end instants.

    python analysis/aquadb_bandwidth_check.py \\
        --log data/terminal_YYYYMMDD_HHMMSS.log \\
        --csv ../AquadbTest.csv --ceiling 8.333
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

HEADER_ROWS = 23
WRAP = 2 ** 32
MIN_MM = 5.0            # ignore short setup hops when judging closure


def load_log(path):
    """DIAG_PROFILE declarations and every DIAG_MOVE_START/END pair."""
    profiles, starts, moves, holds = {}, {}, [], []
    kv = lambda s: dict(p.split("=", 1) for p in s.strip().split(",") if "=" in p)
    for line in Path(path).open(encoding="utf-8", errors="replace"):
        if "# DIAG_PROFILE," in line:
            f = kv(line.split("# DIAG_PROFILE,")[1])
            profiles[f["name"]] = (float(f["v_mm_s"]), float(f["a_mm_s2"]))
        elif "# DIAG_MOVE_START," in line:
            f = kv(line.split("# DIAG_MOVE_START,")[1])
            starts[(f["profile"], f["cycle"], f["leg"])] = f
        elif "# DIAG_MOVE_END," in line:
            f = kv(line.split("# DIAG_MOVE_END,")[1])
            s = starts.get((f["profile"], f["cycle"], f["leg"]))
            if s is None:
                continue
            moves.append(dict(
                profile=f["profile"], cycle=int(f["cycle"]), leg=f["leg"],
                t0=int(s["t_ms"]), t1=int(f["t_ms"]),
                elapsed_ms=int(f["elapsed_ms"]),
                x_from=float(s["x_from_mm"]), x_to=float(s["x_to_mm"]),
                v=float(s["v_mm_s"]), a=float(s["a_mm_s2"]),
                requested=int(f["requested_pulses"]),
                emitted=int(f["emitted_pulses"]),
                complete=int(f["complete"])))
        elif "# DIAG_HOLD_END," in line:
            f = kv(line.split("# DIAG_HOLD_END,")[1])
            holds.append(dict(profile=f["profile"], cycle=int(f["cycle"]),
                              point=f["point"], t=int(f["t_ms"]),
                              x=float(f["x_mm"])))
    return profiles, moves, holds


def scope_start_unix(path):
    """The export header's FILETIME, as a unix epoch (100 ns ticks since 1601)."""
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i > 8:
                break
            if "Starttime of export" in line:
                ft = int(line.split("\t")[1])
                return ft / 1e7 - 11644473600.0
    return None


def log_wall_anchor(path):
    """Median (wall_unix, board_ms) pair from records that carry both."""
    pairs = []
    for line in Path(path).open(encoding="utf-8", errors="replace"):
        m = re.match(r"^(\d{10}\.\d+)\s+# DIAG_\w+,.*?t_ms=(\d+)", line)
        if m:
            pairs.append((float(m.group(1)), int(m.group(2))))
    return pairs


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
    return pos / 1e6


def fit_clock(seed_off, rate=1.0):
    """scope = rate*board + offset, from the two files' own wall clocks.

    Motion-based anchoring cannot be used here: a move END is late by exactly
    the backlog under test, and a move START is clean only where the queue
    drained during the preceding hold, which is false on the fast profiles.
    The export's FILETIME header against the terminal's line timestamps gives
    the offset outright. A rate error of 1e-4 is 0.1 s over a 20 min run,
    well inside the 0.5 s windows, so the rate is left at 1 unless given.
    """
    return rate, seed_off


def main() -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--rate", type=float, default=1.0,
                    help="board-to-scope clock rate; 1.0 is right to 1e-4")
    ap.add_argument("--ceiling", type=float, default=8.333,
                    help="metrology ceiling in mm/s (400 ns -> 2.5, 120 ns -> 8.333)")
    a = ap.parse_args()

    profiles, moves, holds = load_log(a.log)
    t_csv = scope_start_unix(a.csv)
    pairs = log_wall_anchor(a.log)
    seed = (float(np.median([(w - t_csv) * 1000.0 - b for w, b in pairs]))
            if t_csv and pairs else None)
    mx = load_scope(a.csv)
    rate, off = fit_clock(seed, a.rate)
    d = np.diff(mx)
    ceil_nm = round(a.ceiling * 1000.0)

    print(f"log {Path(a.log).name}   csv {Path(a.csv).name}")
    print(f"ceiling {a.ceiling:.3f} mm/s = {ceil_nm} nm per 1 ms sample")
    print(f"clock: scope = {rate:.6f}*board + {off:.0f} ms "
          f"(from the export FILETIME and the log wall clock)")
    print(f"profiles declared: " +
          ", ".join(f"{k} ({v[0]:g} mm/s, {v[1]:g} mm/s^2)"
                    for k, v in profiles.items()) + "\n")

    # Closure is measured settled-to-settled: the position at the end of
    # one hold against the end of the next. Both are drain-free whenever
    # the hold outlasts the backlog, and the commanded displacement
    # between them is exact.
    holds.sort(key=lambda h: h["t"])
    REC = 500
    rows = {}
    for h0, h1 in zip(holds, holds[1:]):
        if h0["profile"] != h1["profile"]:
            continue
        cmd = h1["x"] - h0["x"]
        if abs(cmd) < MIN_MM:
            continue
        w = []
        for h in (h0, h1):
            k = int(round(h["t"] * rate + off))
            seg = mx[max(0, k - REC):min(len(mx), k)]
            if seg.size == 0:
                w = None
                break
            w.append((seg.mean(), seg.std() * 1000.0))
        if w is None:
            continue
        i0 = int(round(h0["t"] * rate + off))
        i1 = int(round(h1["t"] * rate + off))
        seg = np.abs(d[max(0, i0):min(len(d), i1)])
        rows.setdefault(h0["profile"], []).append(dict(
            cmd=abs(cmd), rec=abs(w[1][0] - w[0][0]),
            clip=100.0 * float(np.mean(seg == ceil_nm)) if seg.size else 0.0,
            settled_nm=max(w[0][1], w[1][1])))

    hdr = (f"{'profile':13} {'n':>2} {'cmd mm':>8} {'rec mm':>8} "
           f"{'closure um':>11} {'pulses':>7} {'v_cmd':>6} {'v_real':>7} "
           f"{'clip%':>6} {'settled nm':>11}")
    print(hdr)
    print("-" * len(hdr))
    ok_all = True
    for name, (v, acc) in profiles.items():
        seg = rows.get(name, [])
        mv = [m for m in moves if m["profile"] == name
              and abs(m["x_to"] - m["x_from"]) > MIN_MM]
        if not seg or not mv:
            continue
        pulses_ok = all(m["emitted"] == m["requested"] and m["complete"] == 1
                        for m in mv)
        vreal = np.mean([abs(m["x_to"] - m["x_from"]) /
                         (m["elapsed_ms"] / 1000.0) for m in mv])
        cmd = np.array([r["cmd"] for r in seg])
        rec = np.array([r["rec"] for r in seg])
        closure = 1000.0 * float(np.mean(rec - cmd))
        settled = float(np.max([r["settled_nm"] for r in seg]))
        bad = abs(closure) > 10.0 or not pulses_ok or settled > 200.0
        ok_all &= not bad
        print(f"{name:13} {len(cmd):2d} {cmd.mean():8.3f} {rec.mean():8.3f} "
              f"{closure:+11.1f} {'ok' if pulses_ok else 'SHORT':>7} "
              f"{v:6.1f} {vreal:7.3f} "
              f"{np.mean([r['clip'] for r in seg]):6.1f} {settled:11.0f}"
              f"{'   <-- CHECK' if bad else ''}")

    print("")
    print("  cmd / rec mm   commanded vs recorded displacement, hold to hold")
    print("  closure        recorded - commanded; lost counts show here")
    print("  pulses         the board's own emitted == requested check")
    print("  v_real         commanded distance / firmware move time")
    print("  settled nm     worst sd in the 0.5 s before a hold ended; large")
    print("                 means the queue had not drained, so that row's")
    print("                 closure is not yet trustworthy")
    print("")
    print("verdict: " + ("closure holds on every profile" if ok_all
                        else "FLAGGED - see the marked rows"))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
