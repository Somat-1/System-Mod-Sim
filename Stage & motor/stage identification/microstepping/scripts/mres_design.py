#!/usr/bin/env python3
"""The microstep-resolution (MRES) test, defined as data.

No plotting and no hardware here: this module builds the ordered list of moves,
staircase legs and dwells for the test in `../PROCEDURE.md`, and everything
else reads it. `render_mres.py` draws it; `run_mres.py` executes it.

One session, in two parts:

  Part 1  warm-up            homing to the two end switches, then warm-up.
                             NOT DEFINED YET - see `warmup()`.
  Part 2  stepping sequence  x_ref baseline, six MRES blocks, end reference

and each block is

  B1  configure and align     MRES write at a full step, sync marker   ~0.1 min
  B2  condition               +/-2 mm strokes, fixed-side approach     ~0.6 min
  B3  two-tier staircase      fine tier then coarse, fwd + rev, x2     ~2-7 min

1/256 is the anchor: it runs first and last, bracketing the other four in
shuffled order, so the two anchor runs show whatever changed over the session.
The procedure split the blocks over two sessions only because B4 and B5 made
each block ~16 min; without them all six fit in one ~30 min session.

The procedure's B4 (settling matrix) and B5 (creep / hold) are not part of
this test. The final positions sit on a common grid for every MRES, so B4
could only detect an MRES effect in the approach dynamics, and the ISO
baseline at the chosen MRES covers the distance range later anyway. Creep
depends on friction and hold current rather than MRES, and is measured once
at the chosen MRES during the friction identification.

Conventions worth knowing before reading the numbers:

* Positions are integers in **u256 = 1/256 full step** (39.0625 nm), measured
  from x_start. Every MRES tested divides 256, so every target in every block
  lies on every block's grid, and spans are identical across MRES by
  construction - the design principle of the procedure.
* 100 full steps = 1 mm, so one full step (FS) is 10 um, and one electrical
  period is 4 FS.
* Speeds are full steps/s and are the **same for every MRES**.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass, asdict

# --- rig constants (shared with ../../scripts/sequence_design.py) ----------
FULL_STEPS_PER_MM = 100
SCREW_LEAD_MM = 2.0               # 200 full steps per revolution
FULL_STEP_UM = 1000.0 / FULL_STEPS_PER_MM
ELECTRICAL_PERIOD_FS = 4
TRAVEL_MM = 44.0
U = 256                           # u256 per full step
ACCEL_FULL_STEPS_S2 = 4000        # firmware trapezoid, 40 mm/s^2
MAX_PULSE_RATE_HZ = 20_000        # firmware step-generator ceiling

# --- section 1 parameters --------------------------------------------------
# 1/16 is left out of this session (2026-10-20): its coarse microsteps were
# not a candidate, and dropping it keeps the session to five blocks.
MRES_TESTED = (32, 64, 128, 256)
ANCHOR = 256

# Reversal lost motion R. The ISO 230-2 reversal value from the stage
# identification repeatability phase (analysis/stage_results.json:
# repeatability.reversal_mean = -1.71 um, sd 0.18 um over 25 targets).
# The gain-pair estimate is 2.25 +/- 0.49 um; the repeatability value is used
# because it is the direct, better-determined measurement of the same thing.
REVERSAL_R_UM = 1.71

# Interferometer floor at a 300 ms average. Not yet measured on this rig - the
# procedure expects low single-digit nm. Replace with the Allan minimum.
SIGMA_FLOOR_NM = 3.0              # PLACEHOLDER

SPEED_FS = 70                     # v = 0.7 mm/s, the rig's MODERATE
X_START = 0                       # centre of travel (ZERO)
X_REF = -100 * U                  # -1 mm: drift reference
APPROACH = 10 * U                 # fixed-side run-up, 0.1 mm >> R
SIDE = -1                         # every fixed-side approach comes from below
RUNUP_DWELL_S = 0.3

# --- block structure ---------------------------------------------------------
B1_SETTLE_S = 2.0                 # wait after the MRES write
MARKER_DWELL_S, MARKER_SETTLE_S = 1.0, 0.5

B2_STROKE = 200 * U               # +/-2 mm
B2_STROKES = 3
B2_TURN_S = 0.2
B2_FINAL_DWELL_S = 2.0

B3_COARSE_INC = U // 16           # 1/16 FS, the same spacing for every MRES
B3_SPAN = ELECTRICAL_PERIOD_FS * U   # 4 FS = one electrical cycle
B3_DWELL_S = 0.5
B3_REPEATS = 2
B3_PULSE_S = 0.001                # one pulse per ms within an increment

REF_DWELL_S = 3.0
REF_VISITS = 3

SESSION_BUDGET_S = 60 * 60

# --- fixed driver settings (section 2), written and read back ---------------
# Hold current equals run current and TPOWERDOWN is 0, as in the stage
# identification campaign: the dwells are the measurements, and a current
# reduction part-way through one changes what is measured. If the operating
# configuration reduces hold current, set it here; `check_parameters()` then
# confirms every B3 averaging window starts after the reduction has completed.
DRIVER = dict(
    chopper="spreadcycle", en_spreadCycle=1, toff=5, tbl=24,
    hstrt="firmware default - read back CHOPCONF", hend="firmware default",
    irun_ma=360, ihold_multiplier=1.0, iholddelay=0, tpowerdown=0,
    intpol=0, r_sense_ohm=0.10,
)
TPOWERDOWN_UNIT_S = 2**18 / 12e6  # ~21.8 ms per count at the 12 MHz clock
B3_WINDOW_S = 0.3                 # averaging window at the end of each dwell

# One session: the anchor first and last, the rest shuffled between them with
# a logged seed.
SESSION_ORDER_SEED_OFFSET = 1

# IDS separator signatures: a rapid negative leap, 1 s dwell, return, settle.
# The amplitude identifies what follows, so the interferometer record can be
# segmented without a shared clock. Every one is larger than anything the
# stepping sequence commands as a single excursion below 2 mm.
MARKER_SESSION_FS = 60
MARKER_END_FS = 56


def marker_fs(mres: int) -> int:
    """32, 36, 40, 44, 48 FS for 1/16 ... 1/256."""
    return 16 + 4 * int(math.log2(mres))


# --- derived quantities ------------------------------------------------------

def delta_nm(mres: int) -> float:
    return FULL_STEP_UM * 1000.0 / mres


def fine_span_u256() -> int:
    """S_fine = max(1/4 FS, 3R), rounded UP to the 1/16 FS grid.

    Rounding to 1/16 FS keeps the fine tier an integer number of microsteps at
    every MRES (including 1/16 itself) and lets the coarse tier continue on the
    same 1/16 FS spacing without a remainder.
    """
    s_um = max(0.25 * FULL_STEP_UM, 3 * REVERSAL_R_UM)
    return math.ceil(s_um / FULL_STEP_UM * 16) * B3_COARSE_INC


def resolvable(mres: int) -> bool:
    return delta_nm(mres) >= 3 * SIGMA_FLOOR_NM


def motion_duration_s(distance_u256: int, speed_fs: int) -> float:
    """Duration of the firmware's symmetric accelerate-cruise-decelerate profile."""
    d = abs(distance_u256) / U
    if d == 0:
        return 0.0
    ramp = speed_fs**2 / ACCEL_FULL_STEPS_S2
    if d >= ramp:
        return 2 * speed_fs / ACCEL_FULL_STEPS_S2 + (d - ramp) / speed_fs
    return 2 * math.sqrt(d / ACCEL_FULL_STEPS_S2)


# --- the plan ----------------------------------------------------------------

@dataclass
class Row:
    """One firmware action. `kind` is move, stair or set_mres.

    A stair row is `n_inc` increments of `inc_usteps` microsteps, each followed
    by `inc_dwell_s`: the whole staircase tier in one row, with its count
    recorded so the renderer and the firmware both see the real steps.
    """
    part: str                     # warmup / stepping
    segment: str                  # upload unit: reference, block1..3, end
    block: int                    # 0 outside blocks
    mres: int
    phase: str                    # REF, B1..B3, E
    tier: str                     # B3a / B3b / marker / runup / ...
    kind: str
    label: str
    t_start: float
    t_end: float
    start: int                    # u256
    target: int                   # u256
    speed_fs: int
    dwell_s: float
    n_inc: int = 0
    inc_usteps: int = 0
    inc_dwell_s: float = 0.0
    direction: int = 0
    repeat: int = 0
    note: str = ""


class Builder:
    def __init__(self):
        self.rows: list[Row] = []
        self.t = 0.0
        self.pos = X_START
        self.part = "warmup"
        self.segment = ""
        self.block = 0
        self.mres = 0

    def _add(self, **kw) -> Row:
        kw.setdefault("start", self.pos)
        r = Row(part=self.part, segment=self.segment,
                block=self.block, mres=self.mres, **kw)
        self.rows.append(r)
        self.t = r.t_end
        self.pos = r.target
        return r

    def set_mres(self, mres: int, phase: str):
        if self.pos % U:
            raise ValueError("MRES may only change at a full-step position")
        self.mres = mres
        self._add(phase=phase, tier="config", kind="set_mres",
                  label=f"MRES 1/{mres}", t_start=self.t,
                  t_end=self.t + B1_SETTLE_S, target=self.pos, speed_fs=0,
                  dwell_s=B1_SETTLE_S,
                  note="write CHOPCONF.MRES, read back, log MSCNT; "
                       "interferometer zero = mean of this dwell")

    def move(self, target: int, dwell_s: float, phase: str, tier: str = "",
             label: str = "", **kw) -> Row:
        dur = motion_duration_s(target - self.pos, SPEED_FS) + dwell_s
        return self._add(phase=phase, tier=tier, kind="move", label=label,
                         t_start=self.t, t_end=self.t + dur, target=target,
                         speed_fs=SPEED_FS, dwell_s=dwell_s, **kw)

    def stair(self, n_inc: int, inc_u256: int, sign: int, phase: str,
              tier: str, label: str, **kw) -> Row:
        usteps = inc_u256 * self.mres // U
        dur = n_inc * (usteps * B3_PULSE_S + B3_DWELL_S)
        return self._add(phase=phase, tier=tier, kind="stair", label=label,
                         t_start=self.t, t_end=self.t + dur,
                         target=self.pos + sign * n_inc * inc_u256,
                         speed_fs=0, dwell_s=0.0, n_inc=n_inc,
                         inc_usteps=usteps, inc_dwell_s=B3_DWELL_S, **kw)

    def approach(self, target: int, dwell_s: float, phase: str, label: str,
                 **kw):
        """Arrive at `target` from the fixed side, with a run-up if needed."""
        if (target - self.pos) * SIDE >= 0 and self.pos != target:
            self.move(target + SIDE * APPROACH, RUNUP_DWELL_S, phase,
                      "runup", f"run-up {label}", **kw)
        return self.move(target, dwell_s, phase, "", label, **kw)

    def marker(self, amp_fs: int, phase: str, label: str):
        home = self.pos
        self.move(home - amp_fs * U, MARKER_DWELL_S, phase, "marker",
                  f"marker {amp_fs} FS: {label}")
        self.move(home, MARKER_SETTLE_S, phase, "marker", "marker return")

    def ref_visits(self, phase: str):
        for k in range(REF_VISITS):
            if k:   # leave x_ref below so the next arrival is a fresh approach
                self.move(X_REF + SIDE * APPROACH, RUNUP_DWELL_S, phase,
                          "runup", "leave x_ref")
            self.approach(X_REF, REF_DWELL_S, phase, f"x_ref {k + 1}/{REF_VISITS}")


# --- part 1: warm-up -----------------------------------------------------------

def warmup(b: Builder):
    """Part 1 - homing and warm-up. Intentionally empty for now.

    To be defined once the two end switches are fitted: home to both ends to
    find the travel, set x_start at its centre, then the warm-up of
    PROCEDURE.md (~4 min of +/-5 mm shuttling at v, until x_ref is stable).
    Whatever goes here must leave the stage at x_start, on a full step, so
    part 2 can start from there.
    """
    b.part = "warmup"


# --- part 2: stepping sequence ----------------------------------------------------

def reference(b: Builder, first_mres: int):
    """The session's drift baseline, at the first block's MRES."""
    b.part, b.segment, b.block = "stepping", "reference", 0
    b.set_mres(first_mres, "REF")
    b.marker(MARKER_SESSION_FS, "REF", "session start")
    b.ref_visits("REF")
    b.approach(X_START, 1.0, "REF", "x_start")


def block(b: Builder, index: int, mres: int):
    b.part, b.segment, b.block = "stepping", f"block{index}", index

    # B1 - configure at a full step, then the IDS separator signature.
    b.set_mres(mres, "B1")
    b.marker(marker_fs(mres), "B1", f"block MRES 1/{mres}")

    # B2 - erase the friction history; finish from below.
    for _ in range(B2_STROKES):
        b.move(X_START + B2_STROKE, B2_TURN_S, "B2", "stroke", "+2 mm")
        b.move(X_START - B2_STROKE, B2_TURN_S, "B2", "stroke", "-2 mm")
    b.move(X_START, B2_FINAL_DWELL_S, "B2", "", "x_start from below")

    # B3 - fine tier at single microsteps, then 1/16 FS increments to 4 FS;
    # the reverse sweep starts with its own fine tier, which is where the
    # reversal dead zone is.
    fine = fine_span_u256()
    n_fine = fine * mres // U
    n_coarse = (B3_SPAN - fine) // B3_COARSE_INC
    for rep in range(B3_REPEATS):
        if rep:
            b.approach(X_START, B2_FINAL_DWELL_S, "B3", "x_start, fixed side",
                       repeat=rep)
        for sign, name in ((+1, "fwd"), (-1, "rev")):
            b.stair(n_fine, U // mres, sign, "B3", "B3a",
                    f"{name} fine {n_fine} x 1 ustep", direction=sign, repeat=rep)
            b.stair(n_coarse, B3_COARSE_INC, sign, "B3", "B3b",
                    f"{name} coarse {n_coarse} x 1/16 FS", direction=sign,
                    repeat=rep)


def end_reference(b: Builder):
    b.part, b.segment, b.block = "stepping", "end", 0
    b.marker(MARKER_END_FS, "E", "end of session")
    b.ref_visits("E")
    b.approach(X_START, 1.0, "E", "x_start")


def session_order(seed: int) -> list[int]:
    middle = [m for m in MRES_TESTED if m != ANCHOR]
    random.Random(seed * 10 + SESSION_ORDER_SEED_OFFSET).shuffle(middle)
    return [ANCHOR, *middle, ANCHOR]


def build_warmup() -> list[Row]:
    """Part 1 alone. Empty until homing is defined."""
    b = Builder()
    warmup(b)
    return b.rows


def build_stepping(seed: int = 7):
    """Part 2: (order, rows). Assumes the stage is at x_start on a full step."""
    order = session_order(seed)
    b = Builder()
    reference(b, order[0])
    for i, mres in enumerate(order, 1):
        block(b, i, mres)
    end_reference(b)
    return order, b.rows


def build_single_block(mres: int) -> list[Row]:
    """B1-B3 at one MRES, for bring-up or a repeat. Starts and ends at x_start."""
    b = Builder()
    block(b, 1, mres)
    return b.rows


def build(seed: int = 7):
    return build_stepping(seed)


def check_parameters():
    """Section 1/2 checks that can be made before anything moves."""
    issues = []
    f_step = SPEED_FS * max(MRES_TESTED)
    if f_step > MAX_PULSE_RATE_HZ:
        issues.append(f"step rate {f_step} Hz at 1/256 exceeds {MAX_PULSE_RATE_HZ} Hz")
    reach_mm = (B2_STROKE + APPROACH) / U / FULL_STEPS_PER_MM
    if reach_mm > TRAVEL_MM / 2 - 1.0:
        issues.append(f"B2 strokes reach {reach_mm:.1f} mm from x_start")
    hold_s = DRIVER["tpowerdown"] * TPOWERDOWN_UNIT_S
    if hold_s > B3_DWELL_S - B3_WINDOW_S:
        issues.append("TPOWERDOWN ends inside the B3 averaging window")
    return issues


def summarise(seed: int, order, rows) -> dict:
    fine = fine_span_u256()
    out = {
        "seed": seed,
        "units": "positions in u256 = 1/256 full step = 39.0625 nm; "
                 "100 full steps/mm",
        "parts": {"warmup": "not defined yet (homing to end switches + warm-up)",
                  "stepping": "x_ref baseline, B1-B3 per MRES block, end reference"},
        "sessions": 1,
        "dropped_phases": ["B4 settling matrix", "B5 creep / hold"],
        "parameters": {
            "screw_lead_mm": SCREW_LEAD_MM,
            "s_fs_um": FULL_STEP_UM,
            "delta_nm": {f"1/{m}": delta_nm(m) for m in MRES_TESTED},
            "reversal_R_um": REVERSAL_R_UM,
            "S_fine_um": fine / U * FULL_STEP_UM,
            "S_fine_fs": fine / U,
            "sigma_floor_nm": SIGMA_FLOOR_NM,
            "sigma_floor_is_placeholder": True,
            "v_mm_s": SPEED_FS / FULL_STEPS_PER_MM,
            "f_step_hz_at_256": SPEED_FS * 256,
            "x_start_mm": X_START / U / FULL_STEPS_PER_MM,
            "x_ref_mm": X_REF / U / FULL_STEPS_PER_MM,
            "approach_mm": APPROACH / U / FULL_STEPS_PER_MM,
            "resolvable": {f"1/{m}": resolvable(m) for m in MRES_TESTED},
            "fine_tier_counts": {f"1/{m}": fine * m // U for m in MRES_TESTED},
            "coarse_tier_increments": (B3_SPAN - fine) // B3_COARSE_INC,
        },
        "driver": DRIVER,
        "markers_fs": {"session": MARKER_SESSION_FS, "end": MARKER_END_FS,
                       **{f"1/{m}": marker_fs(m) for m in MRES_TESTED}},
        "checks": check_parameters(),
    }
    segments = {}
    for r in rows:
        d = segments.setdefault(r.segment, {"mres": r.mres, "rows": 0,
                                            "seconds": 0.0, "phases": {}})
        d["rows"] += 1
        d["seconds"] += r.t_end - r.t_start
        d["phases"][r.phase] = d["phases"].get(r.phase, 0.0) + r.t_end - r.t_start
    for d in segments.values():
        d["minutes"] = round(d.pop("seconds") / 60, 2)
        d["phases"] = {k: round(v / 60, 2) for k, v in d["phases"].items()}
    out["order"] = [f"1/{m}" for m in order]
    out["stepping_minutes"] = round(rows[-1].t_end / 60, 2)
    out["segments"] = segments
    return out


def main():
    import argparse
    import csv
    from pathlib import Path

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None, help="output folder (default ../sequence)")
    args = ap.parse_args()

    out = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    order, rows = build(args.seed)
    with (out / "mres_plan.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in asdict(r).items()})
    summary = summarise(args.seed, order, rows)
    (out / "mres_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                          encoding="utf-8")

    p = summary["parameters"]
    print(f"S_fine = {p['S_fine_um']:.3f} um ({p['S_fine_fs']:g} FS); "
          f"fine counts {p['fine_tier_counts']}; "
          f"coarse increments {p['coarse_tier_increments']}")
    print("part 1 (warm-up): not defined yet")
    print(f"session: {' -> '.join(summary['order'])}  "
          f"stepping {summary['stepping_minutes']:.1f} min")
    for seg, d in summary["segments"].items():
        print(f"  {seg:<9} {'1/' + str(d['mres']):>6} {d['rows']:4d} rows "
              f"{d['minutes']:6.1f} min  {d['phases']}")
    for issue in summary["checks"]:
        print(f"CHECK: {issue}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
