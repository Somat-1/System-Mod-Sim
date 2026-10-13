#!/usr/bin/env python3
"""Generate the native-microstep rotor identification campaign."""

from __future__ import annotations
import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

FULL_STEPS_PER_REV = 200
ELECTRICAL_PERIOD_FS = 4
DEFAULT_MRES = 16
DEFAULT_CURRENT_MA = 360
SAMPLE_HZ = 2000
DRIVER_MODE = {
    "spreadcycle": True,
    "microplyer_intpol": False,
    "stealthchop": False,
    "coolstep": False,
    "standstill_current_reduction": False,
}


@dataclass
class Action:
    block: str
    action: str
    label: str
    current_ma: int = DEFAULT_CURRENT_MA
    mres: int = DEFAULT_MRES
    microsteps: int = 0
    speed_fs_s: float = 0.0
    dwell_ms: int = 0
    repeats: int = 1
    sample_hz: int = SAMPLE_HZ
    note: str = ""

    def duration_s(self) -> float:
        n, dwell = abs(self.microsteps), self.dwell_ms / 1000
        if self.action == "hold":
            return dwell
        if self.action in {"move", "velocity"}:
            return n / (self.mres * self.speed_fs_s) + dwell
        if self.action == "stair":
            return n * dwell
        if self.action == "step_response":
            return 2 * n * self.repeats * dwell
        if self.action == "reversal":
            return 4 * n * self.repeats / (self.mres * self.speed_fs_s) + dwell
        raise ValueError(f"unknown action {self.action}")

    def net_u256(self) -> int:
        if self.action in {"move", "velocity", "stair"}:
            return self.microsteps * (256 // self.mres)
        return 0


def hold(block: str, label: str, seconds: float, current: int,
         note: str = "", sample_hz: int = 500) -> Action:
    return Action(block, "hold", label, current_ma=current,
                  dwell_ms=round(1000 * seconds), sample_hz=sample_hz, note=note)


def build_sequence(nominal_current_ma: int = DEFAULT_CURRENT_MA,
                   stair_dwell_ms: int = 150,
                   include_reversal: bool = True) -> list[Action]:
    """Build a <=1 h plan that starts and ends at command zero."""
    if nominal_current_ma <= 0:
        raise ValueError("nominal current must be positive")
    if stair_dwell_ms < 20:
        raise ValueError("stair dwell must be >=20 ms")
    a: list[Action] = []
    mres = DEFAULT_MRES

    # 0/1: opening reference, run-in, and electrical-phase reference.
    a.append(hold("warmup", "energized hold", 60, nominal_current_ma,
                  "opening standstill; log encoder mean and temperature"))
    warm_n = 6 * FULL_STEPS_PER_REV * mres
    for rep in range(1, 7):
        a.append(Action("warmup", "velocity", f"warmup + rep {rep}",
                        nominal_current_ma, mres, +warm_n, 80, 500,
                        sample_hz=500,
                        note="six shaft revolutions"))
        a.append(Action("warmup", "velocity", f"warmup - rep {rep}",
                        nominal_current_ma, mres, -warm_n, 80, 500,
                        sample_hz=500,
                        note="same six revolutions in reverse"))
    align_n = ELECTRICAL_PERIOD_FS * mres
    a.extend([
        Action("alignment", "stair", "electrical cycle +", nominal_current_ma,
               mres, +align_n, dwell_ms=120, sample_hz=500),
        Action("alignment", "stair", "electrical cycle -", nominal_current_ma,
               mres, -align_n, dwell_ms=120, sample_hz=500),
    ])

    # 2: every MRES=16 table entry over a complete shaft revolution.
    rev_n = FULL_STEPS_PER_REV * mres
    a.extend([
        Action("static_map", "stair", "full revolution +", nominal_current_ma,
               mres, +rev_n, dwell_ms=stair_dwell_ms,
               sample_hz=500,
               note="average settled tail of each dwell"),
        Action("static_map", "stair", "full revolution -", nominal_current_ma,
               mres, -rev_n, dwell_ms=stair_dwell_ms,
               sample_hz=500,
               note="identical command positions in reverse"),
    ])

    # 3: current scaling separates detent torque from waveform angle error.
    currents = sorted({round(nominal_current_ma * x) for x in (0.50, 0.75, 1.0)})
    current_n = 4 * ELECTRICAL_PERIOD_FS * mres
    for current in currents:
        a.append(hold("current_sweep", f"settle at {current} mA", 4, current))
        a.extend([
            Action("current_sweep", "stair", f"{current} mA +", current,
                   mres, +current_n, dwell_ms=120,
                   sample_hz=500,
                   note="four electrical cycles"),
            Action("current_sweep", "stair", f"{current} mA -", current,
                   mres, -current_n, dwell_ms=120,
                   sample_hz=500,
                   note="four electrical cycles in reverse"),
        ])

    # 4: twelve out-and-back responses at four mechanical angles.
    position_fs = 0
    for rotor_fs in (0, 50, 100, 150):
        delta = rotor_fs - position_fs
        if delta:
            a.append(Action("dynamics", "move", f"move to {rotor_fs} fs",
                            nominal_current_ma, mres, delta * mres, 40, 1000,
                            sample_hz=1000))
            position_fs = rotor_fs
        for test_mres, name in ((1, "full step"), (16, "one microstep")):
            for direction in (+1, -1):
                sign = "+" if direction > 0 else "-"
                a.append(Action("dynamics", "step_response",
                                f"{name} {sign} at {rotor_fs} fs",
                                nominal_current_ma, test_mres, direction,
                                dwell_ms=350, repeats=12,
                                sample_hz=2000,
                                note="retain each excitation and reset edge"))
    a.append(Action("dynamics", "move", "return from dynamics",
                    nominal_current_ma, mres, -position_fs * mres, 40, 1000,
                    sample_hz=1000))

    # 5: paired log-spaced constant-velocity plateaus.
    for speed in (0.25, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 400):
        span_fs = max(2 * ELECTRICAL_PERIOD_FS,
                      min(3 * FULL_STEPS_PER_REV, round(speed * 12)))
        for direction in (+1, -1):
            sign = "+" if direction > 0 else "-"
            a.append(Action("velocity", "velocity", f"{speed:g} fs/s {sign}",
                            nominal_current_ma, mres,
                            direction * span_fs * mres, speed, 1500,
                            sample_hz=2000,
                            note=f"constant-rate span {span_fs} full steps"))

    # 6: optional native-MRES triangular reversals for pre-sliding friction.
    if include_reversal:
        for frequency in (0.2, 1.0):
            for amplitude in (1, 2, 4, 8, 16, 32):
                speed = 4 * amplitude / mres * frequency
                a.append(Action("reversal", "reversal",
                                f"A={amplitude} us f={frequency:g} Hz",
                                nominal_current_ma, mres, amplitude, speed,
                                1000, 8,
                                sample_hz=1000,
                                note="eight symmetric native-microstep cycles"))

    # 7: closing repeat and drift reference.
    close_n = 2 * ELECTRICAL_PERIOD_FS * mres
    a.extend([
        Action("close", "stair", "repeat two cycles +", nominal_current_ma,
               mres, +close_n, dwell_ms=stair_dwell_ms, sample_hz=500),
        Action("close", "stair", "repeat two cycles -", nominal_current_ma,
               mres, -close_n, dwell_ms=stair_dwell_ms, sample_hz=500),
        hold("close", "final reference standstill", 30, nominal_current_ma,
             "compare with opening reference for drift"),
    ])
    validate_sequence(a)
    return a


def validate_sequence(actions: list[Action]) -> None:
    allowed = {"hold", "move", "stair", "step_response", "velocity", "reversal"}
    position = 0
    for i, row in enumerate(actions):
        if row.action not in allowed:
            raise ValueError(f"row {i}: invalid action")
        if not 1 <= row.mres <= 256 or row.mres & (row.mres - 1):
            raise ValueError(f"row {i}: invalid MRES")
        if not 0 < row.current_ma <= 2000:
            raise ValueError(f"row {i}: invalid current")
        if not 10 <= row.sample_hz <= 5000:
            raise ValueError(f"row {i}: invalid sample rate")
        if row.action in {"move", "velocity", "reversal"}:
            if row.speed_fs_s <= 0 or row.speed_fs_s * row.mres > 20000:
                raise ValueError(f"row {i}: invalid pulse rate")
        if row.duration_s() <= 0:
            raise ValueError(f"row {i}: zero duration")
        position += row.net_u256()
    if position:
        raise ValueError(f"campaign ends at {position} u256, not zero")
    if sum(x.duration_s() for x in actions) > 3600:
        raise ValueError("campaign exceeds one hour")
    static = [x.microsteps for x in actions if x.block == "static_map"]
    expected = FULL_STEPS_PER_REV * DEFAULT_MRES
    if static != [expected, -expected]:
        raise ValueError("static map is not one full revolution out and back")


def timeline(actions: list[Action]) -> list[dict]:
    rows, elapsed, position = [], 0.0, 0
    for index, action in enumerate(actions):
        duration, end = action.duration_s(), position + action.net_u256()
        row = asdict(action)
        row.update(index=index, t_start_s=round(elapsed, 6),
                   t_end_s=round(elapsed + duration, 6),
                   start_u256=position, end_u256=end)
        rows.append(row)
        elapsed, position = elapsed + duration, end
    return rows


def summary(actions: list[Action]) -> dict:
    blocks: dict[str, dict] = {}
    for action in actions:
        item = blocks.setdefault(action.block, {"actions": 0, "seconds": 0.0})
        item["actions"] += 1
        item["seconds"] += action.duration_s()
    total = sum(x.duration_s() for x in actions)
    return {
        "total_seconds": round(total, 3), "total_minutes": round(total / 60, 2),
        "actions": len(actions), "full_steps_per_revolution": FULL_STEPS_PER_REV,
        "electrical_period_full_steps": ELECTRICAL_PERIOD_FS,
        "native_mres": DEFAULT_MRES, "sample_hz": SAMPLE_HZ,
        "driver_mode_invariants": DRIVER_MODE,
        "blocks": {k: {"actions": v["actions"],
                        "minutes": round(v["seconds"] / 60, 2)}
                   for k, v in blocks.items()},
    }


def write_plan(actions: list[Action], out: Path) -> tuple[Path, Path]:
    out.mkdir(parents=True, exist_ok=True)
    rows = timeline(actions)
    plan = out / "motor_identification_plan.csv"
    meta = out / "motor_identification_summary.json"
    with plan.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    meta.write_text(json.dumps(summary(actions), indent=2) + "\n", encoding="utf-8")
    return plan, meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nominal-current-ma", type=int, default=DEFAULT_CURRENT_MA)
    parser.add_argument("--stair-dwell-ms", type=int, default=150)
    parser.add_argument("--no-reversal", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    out = args.out or Path(__file__).resolve().parent.parent / "sequence"
    actions = build_sequence(args.nominal_current_ma, args.stair_dwell_ms,
                             not args.no_reversal)
    plan, meta = write_plan(actions, out)
    info = summary(actions)
    print(f"{info['actions']} actions, {info['total_minutes']:.2f} min")
    for name, item in info["blocks"].items():
        print(f"  {name:<14} {item['actions']:>3} actions  {item['minutes']:>6.2f} min")
    print(f"wrote {plan}")
    print(f"wrote {meta}")


if __name__ == "__main__":
    main()
