#!/usr/bin/env python3
"""The stage identification sequence, defined as data.

No plotting and no hardware here: this module builds the ordered list of moves
and dwells that make up the campaign, and everything else reads it.
`render_sequence.py` draws it; the firmware executes the same structure.

Five phases, about 45 minutes plus interleaved drift checks:

  1  gain calibration        the true step -> distance gain            ~3 min
  2  microstep staircase     static error vs phase, and hysteresis    ~10 min
  3  settling matrix         the main metric                          ~20 min
  4  bidirectional repeat.   ISO 230-2 style accuracy numbers         ~10 min
  D  drift checks            separate drift from the effects          interleaved

Two conventions worth knowing before reading the numbers:

* Positions are **full steps** from the centre of travel. 100 full steps = 1 mm,
  so one full step is 10 um and the working range below is +/-10 mm.
* An **electrical period is 4 full steps** (40 um) and a **screw revolution is
  200 full steps** (2 mm). 200 is a whole number of both, which is why every
  gain-calibration move is a multiple of 200 — see `gain_calibration()`.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass, asdict, field

# --- rig constants ---------------------------------------------------------
FULL_STEPS_PER_MM = 100
FULL_STEPS_PER_REV = 200          # one screw revolution = 2 mm
ELECTRICAL_PERIOD_STEPS = 4       # one electrical period of the 2-phase motor
TRAVEL_MM = 44.0                  # measured total travel
WORKING_HALF_SPAN_MM = 10.0       # +/- about centre; well inside the travel

# --- speeds, as full steps per second --------------------------------------
SPEED_SLOW = 25
SPEED_MODERATE = 70               # the rig's "MODERATE", 0.7 mm/s
SPEED_FAST = 200
SPEED_QUASISTATIC = 10
ACCEL_FULL_STEPS_S2 = 4000       # must match stage_identification.ino

DRIFT_INTERVAL_S = 300.0          # return to the reference every ~5 min
DRIFT_DWELL_S = 3.0
REFERENCE_POSITION = 0            # centre of travel


@dataclass
class Move:
    """One commanded motion, or one dwell when `target` repeats the last one."""
    phase: str
    label: str
    target: int                   # full steps from centre
    speed: int                    # full steps/s
    dwell_s: float                # hold after arriving
    mres: int = 16                # microstep resolution for this move
    note: str = ""
    # A staircase leg is one row here but many discrete microsteps on the stage.
    # Recording the count makes the row self-describing, so the renderer can
    # draw the real steps instead of parsing them back out of `note`.
    micro_steps: int = 0
    micro_dwell_s: float = 0.0
    # Moves sharing a group must not be separated by a drift check. A staircase
    # is forward and reverse over the same ground, and a repeatability approach
    # is a run-up followed by its target: in both cases a detour to the
    # reference in the middle inserts a backlash reversal and measures
    # something else. Empty means the move stands alone.
    group: str = ""

    def distance_from(self, previous: int) -> int:
        return abs(self.target - previous)

    def duration_s(self, previous: int) -> float:
        if self.micro_steps:
            # Firmware emits each staircase pulse over 1 ms, then holds it.
            return self.micro_steps * (0.001 + self.micro_dwell_s)
        distance = self.distance_from(previous)
        if distance == 0:
            return self.dwell_s
        return motion_duration_s(distance, self.speed) + self.dwell_s


@dataclass
class Phase:
    key: str
    title: str
    purpose: str
    method: str
    moves: list = field(default_factory=list)


def mm(steps: float) -> float:
    return steps / FULL_STEPS_PER_MM


def steps(millimetres: float) -> int:
    return int(round(millimetres * FULL_STEPS_PER_MM))


def motion_duration_s(distance_steps: int, speed_steps_s: int) -> float:
    """Duration of the symmetric triangular/trapezoidal firmware profile."""
    if distance_steps <= 0:
        return 0.0
    accel_and_decel_distance = speed_steps_s**2 / ACCEL_FULL_STEPS_S2
    if distance_steps >= accel_and_decel_distance:
        return (2 * speed_steps_s / ACCEL_FULL_STEPS_S2
                + (distance_steps - accel_and_decel_distance) / speed_steps_s)
    return 2 * math.sqrt(distance_steps / ACCEL_FULL_STEPS_S2)


# ---------------------------------------------------------------------------
# 1. gain calibration
# ---------------------------------------------------------------------------

def gain_calibration() -> Phase:
    """Long moves of whole electrical periods *and* whole screw revolutions.

    Both conditions matter and they are different errors:

    * A move that is not a whole number of **electrical periods** starts and
      ends at different points in the motor's electrical cycle, so the static
      phase error does not cancel and leaks straight into the fitted slope.
    * A move that is not a whole number of **screw revolutions** starts and ends
      at different screw angles, so once-per-revolution eccentricity does not
      cancel either. On this rig that term was measured at 159.5 um amplitude -
      far larger than the gain error being looked for, so it has to cancel by
      construction rather than by averaging.

    200 full steps is 50 electrical periods and exactly one screw revolution, so
    every multiple of 200 satisfies both at once.
    """
    ph = Phase(
        key="gain",
        title="1 · Gain calibration",
        purpose="the true step → distance gain",
        method="long moves of whole electrical periods and whole screw "
               "revolutions, both directions; fit the slope",
    )
    # Several lengths, each done out and back, so direction is balanced.
    lengths_rev = [1, 2, 3, 4, 5]        # screw revolutions = 2,4,6,8,10 mm
    position = 0
    for _pass in range(2):           # two passes: replicates for the slope fit
      for rev in lengths_rev:
        span = rev * FULL_STEPS_PER_REV
        for sign in (+1, -1):
            target = sign * span // 2
            ph.moves.append(Move(
                "gain", f"{rev} rev {'+' if sign > 0 else '-'}", target,
                SPEED_MODERATE, 2.0,
                note=f"{span} full steps = {span // ELECTRICAL_PERIOD_STEPS} "
                     f"electrical periods = {rev} screw rev"))
            position = target
    ph.moves.append(Move("gain", "return", 0, SPEED_MODERATE, 2.0))
    return ph


# ---------------------------------------------------------------------------
# 2. microstep staircase
# ---------------------------------------------------------------------------

def microstep_staircase(mres: int = 16, periods: int = 2,
                        dwell_s: float = 0.75,
                        locations_mm=(-5.0, 0.0, +5.0)) -> Phase:
    """One microstep at a time through >=2 electrical periods, then back.

    Forward and back over the same ground is the whole point: the difference
    between the two passes at the same commanded position *is* the hysteresis,
    and it cannot be recovered from a one-way sweep.

    Repeating at several locations separates error that follows the motor's
    electrical phase (which repeats every 4 full steps wherever you are) from
    error that follows the screw (which does not).
    """
    ph = Phase(
        key="staircase",
        title="2 · Microstep staircase",
        purpose="static error map vs electrical phase, and hysteresis",
        method=f"quasi-static staircase through {periods} electrical periods "
               f"forward then back, {dwell_s:g} s per microstep, at "
               f"{len(locations_mm)} locations",
    )
    span_full = periods * ELECTRICAL_PERIOD_STEPS      # full steps
    n_micro = span_full * mres                         # microsteps each way
    for loc in locations_mm:
        base = steps(loc)
        grp = f"staircase@{loc:+.0f}"
        ph.moves.append(Move("staircase", f"to {loc:+.0f} mm", base,
                             SPEED_MODERATE, 2.0, mres=mres, group=grp))
        # The staircase is recorded as its two endpoints plus the microstep
        # count; expanding 128 individual dwells into this list would bury the
        # structure without adding information.
        ph.moves.append(Move(
            "staircase", f"up {n_micro} usteps", base + span_full,
            SPEED_QUASISTATIC, dwell_s * n_micro, mres=mres,
            micro_steps=n_micro, micro_dwell_s=dwell_s, group=grp,
            note=f"{n_micro} microsteps, {dwell_s:g} s dwell each, forward"))
        ph.moves.append(Move(
            "staircase", f"down {n_micro} usteps", base,
            SPEED_QUASISTATIC, dwell_s * n_micro, mres=mres,
            micro_steps=n_micro, micro_dwell_s=dwell_s, group=grp,
            note=f"{n_micro} microsteps, {dwell_s:g} s dwell each, reverse"))
    ph.moves.append(Move("staircase", "return", 0, SPEED_MODERATE, 2.0))
    return ph


# ---------------------------------------------------------------------------
# 3. point-to-point settling matrix
# ---------------------------------------------------------------------------

def settling_matrix(seed: int = 7, repeats: int = 5,
                    dwell_s: float = 5.0) -> Phase:
    """Randomised full-factorial over distance x speed x direction.

    Randomised because the confounders here are all monotonic in time - the
    motor warms, the lubricant redistributes, the room drifts. Run in a tidy
    nested order, any of those would masquerade as an effect of whichever
    factor happened to vary slowest.
    """
    ph = Phase(
        key="settling",
        title="3 · Point-to-point settling matrix",
        purpose="the main metric: where it actually stops",
        method="acceleration-limited moves in randomised order over distance x speed x "
               "direction, 2-5 s dwell after each stop",
    )
    distances_mm = [0.1, 1.0, 10.0]
    speeds = [("slow", SPEED_SLOW), ("moderate", SPEED_MODERATE),
              ("fast", SPEED_FAST)]
    directions = [+1, -1]

    cells = [(d, sn, sv, sg)
             for d in distances_mm
             for sn, sv in speeds
             for sg in directions
             for _ in range(repeats)]
    random.Random(seed).shuffle(cells)

    position = 0
    for dist_mm, sname, sval, sign in cells:
        delta = sign * steps(dist_mm)
        target = position + delta
        # Keep inside the working range by reflecting rather than clipping, so
        # the commanded distance stays exactly as the factor level says.
        limit = steps(WORKING_HALF_SPAN_MM)
        if abs(target) > limit:
            target = position - delta
        ph.moves.append(Move(
            "settling", f"{dist_mm:g} mm {sname} {'+' if delta > 0 else '−'}",
            target, sval, dwell_s,
            note=f"distance={dist_mm:g}mm speed={sname} dir={'+' if delta > 0 else '-'}"))
        position = target
    ph.moves.append(Move("settling", "return", 0, SPEED_MODERATE, 2.0))
    return ph


# ---------------------------------------------------------------------------
# 4. bidirectional repeatability
# ---------------------------------------------------------------------------

def bidirectional_repeatability(n_targets: int = 5, repeats: int = 5,
                                approach_mm: float = 0.5) -> Phase:
    """ISO 230-2 style: each target approached from both sides, several times.

    The approach run-up is what makes it bidirectional rather than merely
    repeated: the stage is taken `approach_mm` past the target in the direction
    of travel and brought back, so the backlash is always taken up the same way
    within one direction and the two directions can be compared.
    """
    ph = Phase(
        key="repeatability",
        title="4 · Bidirectional repeatability",
        purpose="standardised accuracy and repeatability numbers",
        method=f"ISO 230-2 style: {n_targets} targets along travel, approached "
               f"from both sides, {repeats} repeats",
    )
    span = WORKING_HALF_SPAN_MM
    targets_mm = [-span + 2 * span * i / (n_targets - 1) for i in range(n_targets)]
    run_up = steps(approach_mm)
    for rep in range(repeats):
        for t_mm in targets_mm:
            t = steps(t_mm)
            for sign in (+1, -1):
                grp = f"approach@{t_mm:+.1f}/{rep}/{sign}"
                ph.moves.append(Move(
                    "repeatability", f"run-up {t_mm:+.1f}", t - sign * run_up,
                    SPEED_MODERATE, 0.5, group=grp,
                    note=f"approach from {'below' if sign > 0 else 'above'}"))
                ph.moves.append(Move(
                    "repeatability", f"target {t_mm:+.1f}", t,
                    SPEED_MODERATE, 2.5, group=grp,
                    note=f"rep {rep+1}/{repeats}, "
                         f"{'forward' if sign > 0 else 'reverse'} approach"))
    ph.moves.append(Move("repeatability", "return", 0, SPEED_MODERATE, 2.0))
    return ph


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

def build_sequence(seed: int = 7, mres: int = 16):
    """All phases in order, with drift checks interleaved by elapsed time."""
    phases = [
        gain_calibration(),
        microstep_staircase(mres=mres),
        settling_matrix(seed=seed),
        bidirectional_repeatability(),
    ]

    timeline = []
    t = 0.0
    position = 0
    next_drift = DRIFT_INTERVAL_S
    drift_index = 0

    open_group = ""
    for ph in phases:
        for mv in ph.moves:
            # A drift check is due: go to the reference, hold, come back. It is
            # inserted between moves rather than replacing one, so the phase
            # under way is not altered by the measurement of its own drift -
            # but never inside a group, where the detour would be measured as
            # part of the effect.
            # An ungrouped move is always a safe boundary; a grouped one only
            # when it starts a different group from the one just finished.
            at_boundary = (mv.group == "") or (mv.group != open_group)
            if t >= next_drift and at_boundary:
                drift_index += 1
                move_t = motion_duration_s(
                    abs(REFERENCE_POSITION - position), SPEED_MODERATE)
                timeline.append(dict(
                    phase="drift", label=f"drift check {drift_index}",
                    t_start=t, t_end=t + move_t + DRIFT_DWELL_S,
                    start=position, target=REFERENCE_POSITION,
                    speed=SPEED_MODERATE, dwell_s=DRIFT_DWELL_S, mres=mres,
                    micro_steps=0, micro_dwell_s=0.0, group="",
                    note="return to the fixed reference"))
                t += move_t + DRIFT_DWELL_S
                position = REFERENCE_POSITION
                # Measure the interval from this check, not from the grid. A
                # check deferred past a long atomic block would otherwise leave
                # the schedule behind real time and fire several in a burst at
                # the next boundary.
                next_drift = t + DRIFT_INTERVAL_S

            dur = mv.duration_s(position)
            timeline.append(dict(
                phase=mv.phase, label=mv.label, t_start=t, t_end=t + dur,
                start=position, target=mv.target, speed=mv.speed,
                dwell_s=mv.dwell_s, mres=mv.mres,
                micro_steps=mv.micro_steps, micro_dwell_s=mv.micro_dwell_s,
                group=mv.group, note=mv.note))
            t += dur
            position = mv.target
            open_group = mv.group

    return phases, timeline, t


def summarise(phases, timeline, total_s):
    by_phase = {}
    for row in timeline:
        d = by_phase.setdefault(row["phase"], {"moves": 0, "seconds": 0.0})
        d["moves"] += 1
        d["seconds"] += row["t_end"] - row["t_start"]
    return {
        "total_seconds": round(total_s, 1),
        "total_minutes": round(total_s / 60.0, 2),
        "moves": len(timeline),
        "working_range_mm": [-WORKING_HALF_SPAN_MM, WORKING_HALF_SPAN_MM],
        "travel_mm": TRAVEL_MM,
        "full_steps_per_mm": FULL_STEPS_PER_MM,
        "electrical_period_full_steps": ELECTRICAL_PERIOD_STEPS,
        "full_steps_per_rev": FULL_STEPS_PER_REV,
        "acceleration_full_steps_s2": ACCEL_FULL_STEPS_S2,
        "phases": {
            p.key: {
                "title": p.title,
                "purpose": p.purpose,
                "method": p.method,
                "moves": by_phase.get(p.key, {}).get("moves", 0),
                "minutes": round(by_phase.get(p.key, {}).get("seconds", 0) / 60, 2),
            } for p in phases
        },
        "drift": {
            "interval_s": DRIFT_INTERVAL_S,
            "checks": by_phase.get("drift", {}).get("moves", 0),
            "minutes": round(by_phase.get("drift", {}).get("seconds", 0) / 60, 2),
        },
    }


def main():
    import argparse
    import csv
    from pathlib import Path

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--mres", type=int, default=16)
    ap.add_argument("--out", default=None, help="output folder (default ../sequence)")
    args = ap.parse_args()

    out = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    phases, timeline, total = build_sequence(args.seed, args.mres)
    summary = summarise(phases, timeline, total)

    with open(out / "sequence_plan.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(timeline[0].keys()))
        w.writeheader()
        for row in timeline:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in row.items()})
    (out / "sequence_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"{len(timeline)} moves, {total/60:.1f} min total")
    for key, info in summary["phases"].items():
        print(f"  {info['title'].replace(chr(183), '-'):<34} {info['moves']:4d} moves  "
              f"{info['minutes']:6.1f} min")
    print(f"  {'D - Drift checks':<34} {summary['drift']['checks']:4d} checks "
          f"{summary['drift']['minutes']:6.1f} min")
    print(f"wrote {out/'sequence_plan.csv'} and sequence_summary.json")


if __name__ == "__main__":
    main()
