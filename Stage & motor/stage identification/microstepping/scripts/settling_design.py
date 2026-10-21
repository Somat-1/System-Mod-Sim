#!/usr/bin/env python3
"""Settling matrix across microstep sizes - the campaign, defined as data.

The question: does microstep size change how far the stage actually goes when
it is told to go a given distance?

So the measured quantity is a DISPLACEMENT, not a position. For each commanded
distance D the stage settles, moves by D, settles again, and the achieved
displacement is the difference between the two settled values:

    delta  = x_settled_after - x_settled_before
    offset = delta - D                      <- what this campaign is for

Nothing here needs an absolute datum, a drift correction or a position
reference, because every number is a difference taken over a few seconds. That
is the main reason this is a separate sequence from the ISO 230-2 baseline,
which measures absolute position against a datum and answers a different
question.

Structure, per microstep size (one block):

    M1  set MRES at a full step, read back, settle
    M2  condition: 3 strokes of +-2 mm, then back to x_start
    M3  the matrix: for each distance D, walk out n steps and back n steps,
        settling at every stop

Walking out and back rather than returning to x_start between moves gives a
run of same-direction displacements, so the reversal at the turn is isolated
to one move instead of contaminating all of them. n is capped per distance so
the walk stays inside the soft limits.

    python scripts/settling_design.py
    python scripts/settling_design.py --mres 16 32 --repeats 6
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, field
from pathlib import Path
import argparse
import csv
import json
import random

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "sequence"

U256 = 256                     # u256 per full step
FULL_STEPS_PER_MM = 100        # 2 mm lead, 200 full steps per revolution
U256_PER_MM = U256 * FULL_STEPS_PER_MM
FULL_STEP_UM = 10.0

# The microstep sizes under test. 1/256 is carried as a deliberate "too fine"
# control: the staircase data already put the stick-slip lump near 700 nm, so
# 1/256 should degrade rather than improve, and the matrix should show it.
MRES_SET = (16, 32, 64, 256)

# PROCEDURE.md B4 distances, in u256. The smallest is one microstep at the
# coarsest size under test, which is what makes the set comparable: every D
# is a whole number of microsteps at every MRES.
DISTANCES = (
    ("1/16 FS", 16),
    ("1/4 FS", 64),
    ("1 FS", 256),
    ("4 FS", 1024),
    ("16 FS", 4096),
    ("64 FS", 16384),
    ("5 mm", 5 * U256_PER_MM),
    ("20 mm", 20 * U256_PER_MM),
)

REPEATS = 4                    # steps per direction, capped by travel
WALK_LIMIT_U256 = 55 * U256_PER_MM     # keep the walk inside +-55 mm

# Step rate ceiling the firmware will generate reliably. It binds hard here:
# the finest size in the set fixes the speed for ALL of them, because settling
# depends on approach velocity, and letting each MRES run at its own top speed
# would compare speeds rather than microstep sizes.
MAX_PULSE_HZ = 20000
V_TARGET_FS_S = 500            # 5 mm/s, the operational speed, if allowed


def speed_for(mres_set) -> int:
    return min(V_TARGET_FS_S, MAX_PULSE_HZ // max(mres_set))


V_FS_S = speed_for(MRES_SET)
ACCEL_MM_S2 = 40.0
DWELL_S = 3.0                  # settle at every stop
SETTLE_READ_S = 1.0            # settled value = mean of the final second
COND_STROKE_U256 = 2 * U256_PER_MM
MARKER_FS = {16: 32, 32: 36, 64: 40, 128: 44, 256: 48}
SEED = 4201

DRIVER = dict(chopper="spreadcycle", en_spreadCycle=1, toff=5, tbl=24,
              irun_ma=400, ihold_multiplier=1.0, iholddelay=0, tpowerdown=0,
              intpol=0, r_sense_ohm=0.10)


@dataclass
class Row:
    block: int
    mres: int
    phase: str
    kind: str
    label: str
    d_label: str = ""
    d_u256: int = 0            # the commanded displacement this row measures
    direction: int = 0
    index: int = 0             # which step of the walk
    measured: int = 0          # 1 if the analysis reads a displacement here
    start_u256: int = 0
    target_u256: int = 0
    speed_fs: int = V_FS_S
    dwell_s: float = DWELL_S
    t_start: float = 0.0
    t_end: float = 0.0
    note: str = ""


def move_seconds(delta_u256: int, speed_fs: int) -> float:
    """Trapezoidal time for a move, with the ramp included."""
    d_mm = abs(delta_u256) / U256_PER_MM
    if d_mm == 0:
        return 0.0
    v = speed_fs / FULL_STEPS_PER_MM          # mm/s
    t_ramp = v / ACCEL_MM_S2
    d_ramp = v * t_ramp                       # accelerate + decelerate
    if d_mm <= d_ramp:
        return 2.0 * (d_mm / ACCEL_MM_S2) ** 0.5
    return 2.0 * t_ramp + (d_mm - d_ramp) / v


class Builder:
    def __init__(self, repeats: int, seed: int, speed: int = V_FS_S):
        self.rows: list[Row] = []
        self.pos = 0
        self.t = 0.0
        self.repeats = repeats
        self.speed = speed
        self.rng = random.Random(seed)

    def add(self, **kw) -> None:
        kw.setdefault("speed_fs", self.speed)
        r = Row(**kw)
        r.start_u256 = self.pos
        r.t_start = round(self.t, 4)
        dur = (move_seconds(r.target_u256 - self.pos, r.speed_fs)
               if r.kind == "move" else 0.0)
        self.t += dur + r.dwell_s
        r.t_end = round(self.t, 4)
        self.pos = r.target_u256
        self.rows.append(r)

    # ---------------------------------------------------------------- blocks
    def configure(self, block: int, mres: int) -> None:
        self.add(block=block, mres=mres, phase="M1", kind="set_mres",
                 label=f"MRES 1/{mres}", target_u256=self.pos, dwell_s=2.0,
                 note="write CHOPCONF.MRES at a full step, read back, log MSCNT")
        amp = MARKER_FS[mres] * U256
        self.add(block=block, mres=mres, phase="M1", kind="move",
                 label=f"marker {MARKER_FS[mres]} FS", target_u256=-amp,
                 dwell_s=1.0, note="separator signature for the recording")
        self.add(block=block, mres=mres, phase="M1", kind="move",
                 label="marker return", target_u256=0, dwell_s=0.5)

    def condition(self, block: int, mres: int) -> None:
        for k in range(3):
            self.add(block=block, mres=mres, phase="M2", kind="move",
                     label=f"+2 mm {k + 1}/3", target_u256=COND_STROKE_U256,
                     dwell_s=0.2)
            self.add(block=block, mres=mres, phase="M2", kind="move",
                     label=f"-2 mm {k + 1}/3", target_u256=-COND_STROKE_U256,
                     dwell_s=0.2)
        self.add(block=block, mres=mres, phase="M2", kind="move",
                 label="x_start from below", target_u256=0, dwell_s=2.0,
                 note="defined friction state before the matrix")

    def matrix(self, block: int, mres: int) -> None:
        order = list(DISTANCES)
        self.rng.shuffle(order)
        for d_label, d in order:
            n = max(1, min(self.repeats, int(WALK_LIMIT_U256 // d)))
            # Straddle x_start rather than walking out one way from it. The
            # measurement is unaffected - every number is a difference - but
            # the stage then needs half the travel. floor(n/2) keeps the
            # offset a whole number of microsteps at every MRES.
            offset = -(n // 2) * d
            if offset:
                self.add(block=block, mres=mres, phase="M3", kind="move",
                         label=f"{d_label} walk start", target_u256=offset,
                         dwell_s=1.0, note="straddle x_start")
            # Walk out, then walk back. The first move of each run carries the
            # reversal; the rest are same-direction displacements.
            for sign in (+1, -1):
                for k in range(n):
                    self.add(block=block, mres=mres, phase="M3", kind="move",
                             label=f"{d_label} {'+' if sign > 0 else '-'}"
                                   f"{k + 1}/{n}",
                             d_label=d_label, d_u256=d, direction=sign,
                             index=k + 1, measured=1,
                             target_u256=self.pos + sign * d,
                             note="reversal" if k == 0 else "")
            if self.pos != offset:
                raise AssertionError(f"{d_label} walk did not close")
            if self.pos != 0:
                self.add(block=block, mres=mres, phase="M3", kind="move",
                         label=f"{d_label} close", target_u256=0, dwell_s=1.0,
                         note="walk did not close; returning to x_start")

    def build(self, mres_set) -> list[Row]:
        for i, mres in enumerate(mres_set, start=1):
            self.configure(i, mres)
            self.condition(i, mres)
            self.matrix(i, mres)
        return self.rows


def validate(rows: list[Row]) -> list[str]:
    problems = []
    pos = 0
    for i, r in enumerate(rows):
        step = U256 // r.mres
        if r.kind == "move" and (r.target_u256 - pos) % step:
            problems.append(f"row {i}: {r.label} is not a whole microstep "
                            f"at 1/{r.mres}")
        if r.kind == "set_mres" and pos % U256:
            problems.append(f"row {i}: MRES change at {pos} u256, "
                            "not a full step")
        if abs(r.target_u256) > WALK_LIMIT_U256 + COND_STROKE_U256:
            problems.append(f"row {i}: {r.label} reaches "
                            f"{r.target_u256 / U256_PER_MM:.2f} mm")
        if r.speed_fs * r.mres > MAX_PULSE_HZ:
            problems.append(f"row {i}: {r.speed_fs} FS/s at 1/{r.mres} is "
                            f"{r.speed_fs * r.mres} Hz, over {MAX_PULSE_HZ}")
        pos = r.target_u256
    if pos:
        problems.append(f"sequence ends at {pos} u256, not x_start")
    for d_label, d in DISTANCES:
        for m in MRES_SET:
            if d % (U256 // m):
                problems.append(f"{d_label} is not a whole microstep at 1/{m}")
    return problems


def summarise(rows: list[Row], mres_set, repeats: int, seed: int) -> dict:
    measured = [r for r in rows if r.measured]
    per_block: dict[str, dict] = {}
    for r in rows:
        b = per_block.setdefault(f"block{r.block}", dict(
            mres=r.mres, rows=0, measured=0, minutes=0.0))
        b["rows"] += 1
        b["measured"] += r.measured
    for r in rows:
        per_block[f"block{r.block}"]["minutes"] = round(
            max(x.t_end for x in rows if x.block == r.block) / 60
            - min(x.t_start for x in rows if x.block == r.block) / 60, 2)
    return {
        "mres_set": list(mres_set),
        "distances": [{"label": lab, "u256": d,
                       "um": d / U256 * FULL_STEP_UM,
                       "microsteps": {f"1/{m}": d // (U256 // m)
                                      for m in mres_set}}
                      for lab, d in DISTANCES],
        "repeats_per_direction": repeats,
        "seed": seed,
        "rows": len(rows),
        "measured_displacements": len(measured),
        "total_minutes": round(rows[-1].t_end / 60, 2),
        "dwell_s": DWELL_S,
        "settle_read_s": SETTLE_READ_S,
        "speed_fs_s": V_FS_S,
        "pulse_rate_hz": {f"1/{m}": V_FS_S * m for m in mres_set},
        "driver": DRIVER,
        "blocks": per_block,
        "measures": "displacement = settled_after - settled_before; "
                    "offset = displacement - commanded. Relative throughout, "
                    "so no datum or drift correction is required.",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mres", type=int, nargs="+", default=list(MRES_SET))
    ap.add_argument("--repeats", type=int, default=REPEATS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--speed", type=int,
                    help="full-step/s; default is the fastest the finest "
                         "size in the set allows")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    speed = args.speed or speed_for(args.mres)
    if speed * max(args.mres) > MAX_PULSE_HZ:
        raise SystemExit(f"{speed} FS/s at 1/{max(args.mres)} needs "
                         f"{speed * max(args.mres)} Hz, over {MAX_PULSE_HZ}")
    rows = Builder(args.repeats, args.seed, speed).build(args.mres)
    problems = validate(rows)
    info = summarise(rows, args.mres, args.repeats, args.seed)
    info["speed_fs_s"] = speed
    info["speed_mm_s"] = speed / FULL_STEPS_PER_MM
    info["pulse_rate_hz"] = {f"1/{m}": speed * m for m in args.mres}

    args.out.mkdir(parents=True, exist_ok=True)
    plan = args.out / "settling_plan.csv"
    with plan.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(rows[0])))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    (args.out / "settling_summary.json").write_text(
        json.dumps(info | {"checks": problems}, indent=2) + "\n",
        encoding="utf-8")

    print(f"{info['rows']} rows, {info['measured_displacements']} measured "
          f"displacements, {info['total_minutes']:.1f} min")
    print(f"microstep sizes: {', '.join(f'1/{m}' for m in args.mres)}")
    print(f"{'distance':>9} {'um':>10} " +
          " ".join(f"{'1/'+str(m):>8}" for m in args.mres) + "   microsteps")
    for d in info["distances"]:
        print(f"{d['label']:>9} {d['um']:10.3f} " +
              " ".join(f"{d['microsteps']['1/'+str(m)]:8d}" for m in args.mres))
    print()
    for name, b in info["blocks"].items():
        print(f"  {name}  1/{b['mres']:<4d} {b['rows']:4d} rows  "
              f"{b['measured']:4d} measured  {b['minutes']:6.2f} min")

    print(f"\nspeed {speed} FS/s = {speed / FULL_STEPS_PER_MM:.2f} mm/s, the "
          f"same for every size - set by 1/{max(args.mres)} against the "
          f"{MAX_PULSE_HZ // 1000} kHz ceiling")
    print("  pulse rate: " + ", ".join(f"1/{m} = {speed * m / 1000:.1f} kHz"
                                       for m in args.mres))
    if speed < V_TARGET_FS_S:
        coarser = [m for m in args.mres if m != max(args.mres)]
        print(f"  NOTE: below the {V_TARGET_FS_S / FULL_STEPS_PER_MM:.1f} mm/s "
              f"operational speed. Without 1/{max(args.mres)} the set would "
              f"allow {speed_for(coarser)} FS/s.")
    if problems:
        print("\nCHECKS FAILED:")
        for p in problems:
            print(f"  {p}")
    else:
        print("\nchecks passed")
    print(f"wrote {plan}")


if __name__ == "__main__":
    main()
