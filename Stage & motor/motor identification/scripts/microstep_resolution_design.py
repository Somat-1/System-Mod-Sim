#!/usr/bin/env python3
"""Generate the microstep-resolution test: how finely can the motor step?

A dedicated companion to the full identification campaign, under 5 minutes.
It runs on the same firmware and runner (scripts/run_campaign.py --plan ...)
and turns one full shaft revolution, one microstep at a time, at each of the
fine resolutions 1/32, 1/64 and 1/128 (change with --mres):

  align     1 s at MRES 1, so the driver's microstep counter starts on a
            full step (it survives ESP32 resets; earlier moves can leave it
            between full steps, and an MRES change then snaps the rotor)
  settle    1 s hold after each MRES change
  fwd       one revolution forward in single microsteps (200 x MRES of them)
  back      the same revolution in reverse (hysteresis = fwd - back)
  pm1       +1 / -1 microstep, repeated, 200 ms each side: does one microstep
            move the rotor, by how much, and how does it settle (2 kHz)

Each revolution takes about 42 s each way, so the hold after each microstep is
7 ms at 1/32, 3 ms at 1/64 and 2 ms at 1/128: a slow quasi-static rotation
(about 1.4 rpm) sampled at 2 kHz through every microstep. Every block returns
to the start angle, a full step, which is on every MRES grid.

Encoder: 20000 counts/rev = 100 counts per full step. One microstep is
100/MRES counts: 6.25 at 1/16, 1.56 at 1/64, 0.39 at 1/256. Individual steps
below about 1/64 are therefore under the encoder's resolution; their average
size over the staircase still is not.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from sequence_design import (DEFAULT_CURRENT_MA, DRIVER_MODE, ELECTRICAL_PERIOD_FS,
                             FULL_STEPS_PER_REV, Action, timeline)

MRES_LIST = (32, 64, 128)         # the fine end, where the question is
MAX_SECONDS = 300.0                # the whole test fits in 5 minutes
REV_SECONDS = 42.0                 # time per revolution, each direction
COUNTS_PER_REV = 20000             # measured at bring-up, 2026-10-13


def stair_dwell_ms(m: int, rev_seconds: float = REV_SECONDS) -> int:
    return max(1, round(rev_seconds * 1000 / (FULL_STEPS_PER_REV * m)))


def stair_sample_hz(dwell: int) -> int:
    # Several readings inside even the shortest hold.
    return 2000 if dwell < 10 else (1000 if dwell < 50 else 500)


def build_sequence(current_ma: int = DEFAULT_CURRENT_MA, dwell_ms: int = 200,
                   pm1_repeats: int = 10, mres_list=MRES_LIST,
                   rev_seconds: float = REV_SECONDS) -> list[Action]:
    a: list[Action] = []
    # The driver keeps its own microstep counter (MSCNT) while powered, and
    # earlier moves can leave it between full steps; switching MRES then snaps
    # the rotor to the new grid. One second at MRES 1 puts it on a full step,
    # so every block below starts from the same, aligned angle.
    a.append(Action("align", "hold", "full-step alignment", current_ma, 1,
                    dwell_ms=1000, sample_hz=500,
                    note="snap MSCNT to a full step; not analysed"))
    for m in mres_list:
        n = FULL_STEPS_PER_REV * m
        sd = stair_dwell_ms(m, rev_seconds)
        hz = stair_sample_hz(sd)
        tag = f"m{m:03d}"
        a.append(Action(f"{tag}_settle", "hold", f"1/{m} settle", current_ma, m,
                        dwell_ms=1000, sample_hz=500,
                        note="hold after the MRES change"))
        a.append(Action(f"{tag}_fwd", "stair", f"1/{m} rev +{n}", current_ma, m,
                        microsteps=n, dwell_ms=sd, sample_hz=hz,
                        note=f"one revolution forward, {sd} ms per microstep"))
        a.append(Action(f"{tag}_back", "stair", f"1/{m} rev -{n}", current_ma, m,
                        microsteps=-n, dwell_ms=sd, sample_hz=hz,
                        note="same revolution in reverse"))
        a.append(Action(f"{tag}_pm1", "step_response", f"1/{m} +-1 x{pm1_repeats}",
                        current_ma, m, microsteps=1, dwell_ms=dwell_ms,
                        repeats=pm1_repeats, sample_hz=2000,
                        note="single-microstep response, both directions"))
    return a


def validate(actions: list[Action]) -> None:
    position = 0
    for i, row in enumerate(actions):
        if not 1 <= row.mres <= 256 or row.mres & (row.mres - 1):
            raise ValueError(f"row {i}: invalid MRES")
        if not 10 <= row.sample_hz <= 5000:
            raise ValueError(f"row {i}: invalid sample rate")
        if position % (256 // row.mres):
            raise ValueError(f"row {i}: start is off the 1/{row.mres} grid")
        position += row.net_u256()
    if position:
        raise ValueError(f"sequence ends at {position} u256, not zero")
    total = sum(x.duration_s() for x in actions)
    if total > MAX_SECONDS:
        raise ValueError(f"sequence takes {total:.0f} s, over {MAX_SECONDS:.0f} s")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--current-ma", type=int, default=DEFAULT_CURRENT_MA)
    ap.add_argument("--mres", type=lambda t: tuple(int(x) for x in t.split(",")),
                    default=MRES_LIST, help="comma-separated, default 32,64,128")
    ap.add_argument("--dwell-ms", type=int, default=200,
                    help="hold each side of the +-1 test (default 200)")
    ap.add_argument("--rev-seconds", type=float, default=REV_SECONDS,
                    help="time per revolution per direction (default 120)")
    ap.add_argument("--pm1-repeats", type=int, default=10)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent.parent / "sequence")
    args = ap.parse_args()
    actions = build_sequence(args.current_ma, args.dwell_ms, args.pm1_repeats,
                             mres_list=args.mres, rev_seconds=args.rev_seconds)
    validate(actions)
    rows = timeline(actions)
    args.out.mkdir(parents=True, exist_ok=True)
    plan = args.out / "microstep_resolution_plan.csv"
    with plan.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    total = rows[-1]["t_end_s"]
    meta = {
        "total_seconds": round(total, 2), "total_minutes": round(total / 60, 2),
        "actions": len(rows), "mres": list(args.mres), "pm1_dwell_ms": args.dwell_ms,
        "rev_seconds": args.rev_seconds,
        "stair_dwell_ms": {str(m): stair_dwell_ms(m, args.rev_seconds) for m in args.mres},
        "pm1_repeats": args.pm1_repeats, "current_ma": args.current_ma,
        "counts_per_rev": COUNTS_PER_REV,
        "counts_per_microstep": {str(m): COUNTS_PER_REV / FULL_STEPS_PER_REV / m
                                 for m in args.mres},
        "driver_mode_invariants": DRIVER_MODE,
    }
    (args.out / "microstep_resolution_summary.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"{len(rows)} actions, {total / 60:.2f} min")
    for m in args.mres:
        sd = stair_dwell_ms(m, args.rev_seconds)
        print(f"  1/{m:<3d}  {FULL_STEPS_PER_REV * m:6d} microsteps/rev, hold {sd:4d} ms, "
              f"{stair_sample_hz(sd)} Hz, 1 microstep = "
              f"{COUNTS_PER_REV / FULL_STEPS_PER_REV / m:6.2f} counts")
    print(f"wrote {plan}")


if __name__ == "__main__":
    main()
