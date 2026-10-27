#!/usr/bin/env python3
"""Generate, validate, and export the E1 measurement-only sequence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "e1_config.json"
SEQUENCE_DIR = ROOT / "sequence"
FIRMWARE_HEADER = ROOT / "firmware" / "e1_spectral_windows" / "e1_config.h"


@dataclass(frozen=True)
class Action:
    order: int
    state: str
    action: str
    pass_number: int | None = None
    window: str = ""
    direction: str = ""
    index: int | None = None
    target_usteps: int | None = None
    target_mm: float | None = None
    settle_s: float = 0.0
    record_s: float = 0.0
    record: bool = False
    note: str = ""


def load_config(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def nominal_stop_s(c: dict[str, Any]) -> float:
    distance = c["nominal_spacing_mm"]
    velocity = c["velocity_mm_s"]
    acceleration = c["acceleration_mm_s2"]
    ramp_distance = velocity * velocity / acceleration
    move = (2.0 * math.sqrt(distance / acceleration)
            if distance <= ramp_distance else
            2.0 * velocity / acceleration + (distance - ramp_distance) / velocity)
    return move + c["settle_s"] + c["record_s"]


def validate_config(c: dict[str, Any]) -> None:
    required = {
        "lead_mm", "steps_per_rev", "mres", "ustep_per_mm", "velocity_mm_s",
        "acceleration_mm_s2", "settle_s", "record_s", "overrun_mm",
        "x_ref_mm", "ref_preapproach_mm", "soft_min_mm", "soft_max_mm",
        "windows_mm", "nominal_spacing_mm", "grid_spacing_usteps",
        "motor_cycle_usteps", "motor_cycles_per_spacing",
        "pilot_targets", "datum_visits", "end_ref_visits", "time_limit_min",
        "pilot_pause_budget_s", "measurement_passes", "aquadb_clock_ns",
        "rmt_resolution_hz", "rmt_chunk_pulses", "rmt_queue_depth",
        "rmt_ring_buffers", "sync_gpio",
    }
    missing = sorted(required - c.keys())
    if missing:
        raise ValueError(f"missing config keys: {', '.join(missing)}")
    expected_upmm = c["steps_per_rev"] * c["mres"] / c["lead_mm"]
    if not math.isclose(c["ustep_per_mm"], expected_upmm):
        raise ValueError(
            f"ustep_per_mm={c['ustep_per_mm']} but mechanics imply {expected_upmm:g}"
        )
    if c["mres"] <= 0 or c["mres"] & (c["mres"] - 1):
        raise ValueError("mres must be a positive power of two")
    if c["settle_s"] < c.get("t_hold_s", 0.0) + 0.5:
        raise ValueError("settle_s must be at least t_hold_s + 0.5 s")
    if c["hold_current_ma"] != c["run_current_ma"] and c.get("t_hold_s", 0.0) <= 0:
        raise ValueError("t_hold_s must be derived when hold current differs from run current")
    if len(c["windows_mm"]) != 3:
        raise ValueError("E1 requires exactly three windows")
    spacing_u = round(c["nominal_spacing_mm"] * c["ustep_per_mm"])
    if spacing_u != c["grid_spacing_usteps"]:
        raise ValueError("nominal_spacing_mm and grid_spacing_usteps disagree")
    if spacing_u != c["motor_cycle_usteps"] * c["motor_cycles_per_spacing"]:
        raise ValueError("grid spacing is not the configured whole number of motor cycles")
    if c["measurement_passes"] not in (1, 2):
        raise ValueError("measurement_passes must be 1 or 2")
    if c["measurement_passes"] == 2 and nominal_stop_s(c) > 1.9:
        raise ValueError("two passes are permitted only when nominal t_stop <= 1.9 s")
    if c["rmt_ring_buffers"] <= c["rmt_queue_depth"]:
        raise ValueError("rmt_ring_buffers must exceed rmt_queue_depth")


def generate_targets(c: dict[str, Any]) -> dict[str, list[int]]:
    """Generate a plain grid whose spacing is an integer motor-cycle count."""
    upmm = int(c["ustep_per_mm"])
    step_u = int(c["grid_spacing_usteps"])
    result: dict[str, list[int]] = {}
    for number, (a_mm, b_mm) in enumerate(c["windows_mm"], 1):
        a_u, b_u = round(a_mm * upmm), round(b_mm * upmm)
        values = list(range(a_u, b_u + 1, step_u))
        if len(values) != 101:
            raise ValueError(f"W{number} has {len(values)} targets, expected 101")
        if any(b - a != step_u for a, b in zip(values, values[1:])):
            raise ValueError(f"W{number} is not a regular {step_u}-microstep grid")
        if values[0] != a_u or values[-1] != b_u:
            raise ValueError(f"W{number} edges are not pinned")
        result[f"W{number}"] = values
    return result


class PlanBuilder:
    def __init__(self, c: dict[str, Any]):
        self.c = c
        self.rows: list[Action] = []

    def add(self, state: str, action: str, *, pass_number: int | None = None,
            window: str = "", direction: str = "", index: int | None = None,
            target_u: int | None = None, settle_s: float = 0.0,
            record_s: float = 0.0, record: bool = False, note: str = "") -> None:
        target_mm = None if target_u is None else target_u / self.c["ustep_per_mm"]
        self.rows.append(Action(len(self.rows) + 1, state, action, pass_number,
                                window, direction, index, target_u, target_mm,
                                settle_s, record_s, record, note))

    def transition(self, state: str, note: str = "", pass_number: int | None = None) -> None:
        self.add(state, "state_enter", pass_number=pass_number, note=note)

    def move(self, state: str, target_mm: float, *, pass_number: int | None = None,
             window: str = "", note: str = "") -> None:
        self.add(state, "move", pass_number=pass_number, window=window,
                 target_u=round(target_mm * self.c["ustep_per_mm"]), note=note)

    def stop(self, state: str, window: str, direction: str, index: int,
             target_u: int, *, pass_number: int | None = None, note: str = "") -> None:
        self.add(state, "stop", pass_number=pass_number, window=window,
                 direction=direction, index=index, target_u=int(target_u),
                 settle_s=self.c["settle_s"], record_s=self.c["record_s"],
                 record=True, note=note)

    def ref_visit(self, state: str, index: int, *, pass_number: int | None = None) -> None:
        self.move(state, self.c["ref_preapproach_mm"], pass_number=pass_number,
                  window="REF", note="fixed-side pre-approach")
        self.stop(state, "REF", "up", index,
                  round(self.c["x_ref_mm"] * self.c["ustep_per_mm"]),
                  pass_number=pass_number, note="reference visit")

    def window(self, name: str, targets: list[int], pass_number: int) -> None:
        a_mm, b_mm = self.c["windows_mm"][int(name[1:]) - 1]
        self.transition(name, f"measurement pass {pass_number}", pass_number)
        self.move(name, a_mm - self.c["overrun_mm"], pass_number=pass_number,
                  window=name, note="lower pre-approach")
        for i, target in enumerate(targets, 1):
            self.stop(name, name, "up", i, int(target), pass_number=pass_number)
        self.move(name, b_mm + self.c["overrun_mm"], pass_number=pass_number,
                  window=name, note="upper overrun")
        for i in range(len(targets), 0, -1):
            self.stop(name, name, "down", i, int(targets[i - 1]),
                      pass_number=pass_number)


def build_measurement_plan(c: dict[str, Any], targets: dict[str, list[int]]) -> list[Action]:
    b = PlanBuilder(c)
    b.transition("INIT", "precondition: homed, warm, at x_ref, interferometer recording")
    b.add("INIT", "hardware_checks",
          note="verify driver, position precondition, signal, AquadB, limits, and logs")
    b.transition("DATUM")
    for visit in range(1, int(c["datum_visits"]) + 1):
        b.ref_visit("DATUM", visit)

    b.transition("PILOT", "first 20 W1 targets in both directions")
    pilot = targets["W1"][: int(c["pilot_targets"])]
    first_mm = int(pilot[0]) / c["ustep_per_mm"]
    last_mm = int(pilot[-1]) / c["ustep_per_mm"]
    b.move("PILOT", first_mm - c["overrun_mm"], window="PILOT",
           note="lower pre-approach")
    for i, target in enumerate(pilot, 1):
        b.stop("PILOT", "PILOT", "up", i, int(target))
    b.move("PILOT", last_mm + c["overrun_mm"], window="PILOT",
           note="upper overrun before downward pilot")
    for i in range(len(pilot), 0, -1):
        b.stop("PILOT", "PILOT", "down", i, int(pilot[i - 1]))
    b.move("PILOT", c["ref_preapproach_mm"], window="REF",
           note="return to reference from below")
    b.move("PILOT", c["x_ref_mm"], window="REF", note="park for pilot review")
    b.add("PILOT", "operator_gate",
          note="continue only after all Section 7 checks pass; otherwise ABORT")

    ref_index = 0
    for pass_number in range(1, int(c["measurement_passes"]) + 1):
        for name in ("W1", "W2", "W3"):
            b.window(name, targets[name], pass_number)
            if name in ("W1", "W2"):
                ref_index += 1
                b.transition("REF", f"post-{name}", pass_number)
                b.ref_visit("REF", ref_index, pass_number=pass_number)
            if name == "W2":
                b.add("REF", "time_watchdog", pass_number=pass_number,
                      note="skip remaining W3 work if projected end exceeds 60 min")

    b.transition("END_REF", "after final W3: three visits and environment snapshot")
    for visit in range(1, int(c["end_ref_visits"]) + 1):
        b.ref_visit("END_REF", visit)
    b.add("END_REF", "environment_snapshot")
    b.transition("PARK", "remain at x_ref; stop recording; write summary")
    b.add("PARK", "write_summary")
    b.transition("DONE")
    return b.rows


def move_time_s(distance_mm: float, c: dict[str, Any]) -> float:
    distance_mm = abs(distance_mm)
    velocity = c["velocity_mm_s"]
    acceleration = c["acceleration_mm_s2"]
    ramp_distance = velocity * velocity / acceleration
    if distance_mm <= ramp_distance:
        return 2.0 * math.sqrt(distance_mm / acceleration)
    return 2.0 * velocity / acceleration + (distance_mm - ramp_distance) / velocity


def plan_time_s(rows: list[Action], c: dict[str, Any]) -> float:
    elapsed = 0.0
    position = c["x_ref_mm"]
    for row in rows:
        if row.target_mm is not None:
            elapsed += move_time_s(row.target_mm - position, c)
            position = row.target_mm
        elapsed += row.settle_s + row.record_s
    return elapsed


def validate_plan(c: dict[str, Any], targets: dict[str, list[int]],
                  measurement: list[Action]) -> dict[str, Any]:
    upmm = c["ustep_per_mm"]
    lo_u, hi_u = round(c["soft_min_mm"] * upmm), round(c["soft_max_mm"] * upmm)
    outside = [r for r in measurement if r.target_usteps is not None
               and not lo_u <= r.target_usteps <= hi_u]
    if outside:
        r = outside[0]
        raise ValueError(f"action {r.order} target {r.target_mm:.6f} mm is outside soft limits")
    window_stops = [r for r in measurement if r.state in ("W1", "W2", "W3")
                    and r.action == "stop"]
    pilot_stops = [r for r in measurement if r.state == "PILOT" and r.action == "stop"]
    expected_window_stops = 2 * sum(len(values) for values in targets.values()) * c["measurement_passes"]
    if len(window_stops) != expected_window_stops:
        raise ValueError(f"measurement has {len(window_stops)} window stops, expected {expected_window_stops}")
    if len(pilot_stops) != 2 * c["pilot_targets"]:
        raise ValueError("pilot stop count is wrong")
    refs = [r for r in measurement if r.action == "stop" and r.window == "REF"]
    estimate_s = plan_time_s(measurement, c) + c["pilot_pause_budget_s"]
    aquadb_limit = 1000.0 / c["aquadb_clock_ns"]
    summary = {
        "targets_per_window": {name: len(values) for name, values in targets.items()},
        "total_unique_targets": sum(len(v) for v in targets.values()),
        "measurement_passes": c["measurement_passes"],
        "window_measurement_stops": len(window_stops),
        "pilot_stops": len(pilot_stops),
        "datum_stops": c["datum_visits"],
        "post_window_ref_stops": 2 * c["measurement_passes"],
        "end_ref_stops": c["end_ref_visits"],
        "measurement_reference_stops_total": len(refs),
        "grid_move_peak_velocity_mm_s": min(
            c["velocity_mm_s"],
            math.sqrt(c["acceleration_mm_s2"] * c["nominal_spacing_mm"]),
        ),
        "aquadb_velocity_limit_mm_s": aquadb_limit,
        "estimated_nominal_stop_s": nominal_stop_s(c),
        "estimated_sequence_without_operator_pause_min": plan_time_s(measurement, c) / 60.0,
        "estimated_recording_with_pilot_budget_min": estimate_s / 60.0,
        "checks": {
            "measurement_only_no_home_or_warmup": not any(
                r.state in ("HOME", "WARMUP") for r in measurement
            ),
            "all_targets_integer_microsteps": True,
            "all_planned_positions_inside_soft_limits": True,
            "each_window_has_101_targets": all(len(values) == 101 for values in targets.values()),
            "all_window_targets_strictly_increasing": True,
            "window_edges_pinned": True,
            "grid_spacing_is_three_motor_cycles": (
                c["grid_spacing_usteps"] == 3 * c["motor_cycle_usteps"]
            ),
            "motor_cycle_phase_is_constant_within_each_window": all(
                len({value % c["motor_cycle_usteps"] for value in values}) == 1
                for values in targets.values()
            ),
            "all_targets_share_one_motor_cycle_phase": len({
                value % c["motor_cycle_usteps"]
                for values in targets.values() for value in values
            }) == 1,
            "velocity_not_above_aquadb_limit": c["velocity_mm_s"] <= aquadb_limit,
            "nominal_stop_not_above_3p5_s": nominal_stop_s(c) <= 3.5,
            "estimated_recording_not_above_60_min": estimate_s <= c["time_limit_min"] * 60.0,
        },
    }
    if not all(summary["checks"].values()):
        failed = [name for name, passed in summary["checks"].items() if not passed]
        raise ValueError("failed plan checks: " + ", ".join(failed))
    return summary


def write_actions(path: Path, rows: list[Action]) -> None:
    fields = list(Action.__dataclass_fields__)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            values = asdict(row)
            if values["target_mm"] is not None:
                values["target_mm"] = f"{values['target_mm']:.6f}"
            writer.writerow(values)


def write_targets(path: Path, c: dict[str, Any], targets: dict[str, list[int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["window", "index", "target_usteps", "target_mm"])
        for name, values in targets.items():
            for i, target in enumerate(values, 1):
                writer.writerow([name, i, int(target), f"{target / c['ustep_per_mm']:.6f}"])


def write_firmware_header(c: dict[str, Any], digest: str,
                          targets: dict[str, list[int]]) -> None:
    def array(values: list[int]) -> str:
        ints = [int(v) for v in values]
        return ",\n".join(
            "  " + ", ".join(str(v) for v in ints[start:start + 12])
            for start in range(0, len(ints), 12)
        )

    lines = [
        "/* GENERATED by e1_sequence.py from e1_config.json. Do not edit. */",
        "#pragma once", "#include <stdint.h>", "",
        f'constexpr char CFG_HASH[] = "{digest[:12]}";',
        f"constexpr double CFG_LEAD_MM = {c['lead_mm']:.6f};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {c['steps_per_rev']};",
        f"constexpr uint16_t CFG_MRES = {c['mres']};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {c['run_current_ma']};",
        f"constexpr uint16_t CFG_HOLD_CURRENT_MA = {c['hold_current_ma']};",
        f"constexpr uint8_t CFG_TPOWERDOWN = {c['tpowerdown']};",
        f"constexpr uint8_t CFG_IHOLDDELAY = {c['iholddelay']};",
        f"constexpr double CFG_V_MM_S = {c['velocity_mm_s']:.6f};",
        f"constexpr double CFG_A_MM_S2 = {c['acceleration_mm_s2']:.6f};",
        f"constexpr uint32_t CFG_SETTLE_MS = {round(c['settle_s'] * 1000)};",
        f"constexpr uint32_t CFG_RECORD_MS = {round(c['record_s'] * 1000)};",
        f"constexpr uint32_t CFG_T_HOLD_MS = {round(c['t_hold_s'] * 1000)};",
        f"constexpr uint32_t CFG_NOMINAL_STOP_MS = {round(nominal_stop_s(c) * 1000)};",
        f"constexpr int32_t CFG_OVERRUN_U = {round(c['overrun_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_XREF_U = {round(c['x_ref_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_REF_PRE_U = {round(c['ref_preapproach_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {round(c['soft_min_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {round(c['soft_max_mm'] * c['ustep_per_mm'])};",
        f"constexpr uint16_t CFG_GRID_SPACING_U = {c['grid_spacing_usteps']};",
        f"constexpr uint8_t CFG_MOTOR_CYCLES_PER_SPACING = {c['motor_cycles_per_spacing']};",
        f"constexpr uint8_t CFG_PILOT_TARGETS = {c['pilot_targets']};",
        f"constexpr uint8_t CFG_DATUM_VISITS = {c['datum_visits']};",
        f"constexpr uint8_t CFG_END_REF_VISITS = {c['end_ref_visits']};",
        f"constexpr uint8_t CFG_MEASUREMENT_PASSES = {c['measurement_passes']};",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(c['time_limit_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {c['aquadb_clock_ns']};",
        f"constexpr float CFG_EXPECTED_NOISE_UM = {c['expected_noise_um']:.6f}F;",
        f"constexpr float CFG_REF_DRIFT_WARNING_UM = {c['reference_drift_warning_um']:.6f}F;",
        f"constexpr uint32_t CFG_RMT_HZ = {c['rmt_resolution_hz']}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {c['rmt_chunk_pulses']};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {c['rmt_queue_depth']};",
        f"constexpr uint8_t CFG_RMT_RING = {c['rmt_ring_buffers']};",
        f"constexpr int8_t CFG_SYNC_GPIO = {c['sync_gpio']};",
        f"constexpr uint16_t CFG_TARGETS_PER_WINDOW = {len(targets['W1'])};",
    ]
    for name in ("W1", "W2", "W3"):
        lines.extend(["", f"constexpr int32_t CFG_{name}_TARGETS[CFG_TARGETS_PER_WINDOW] = {{",
                      array(targets[name]), "};"])
    FIRMWARE_HEADER.parent.mkdir(parents=True, exist_ok=True)
    FIRMWARE_HEADER.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    c, digest = load_config(args.config.resolve())
    validate_config(c)
    targets = generate_targets(c)
    measurement = build_measurement_plan(c, targets)
    summary = validate_plan(c, targets, measurement)
    summary.update({
        "experiment": c["experiment"], "config_sha256": digest, "config": c,
        "target_grid": "regular_192_usteps_three_motor_cycles",
        "ustep_per_mm": c["ustep_per_mm"],
        "microstep_um": 1000.0 / c["ustep_per_mm"],
        "electrical_period_usteps": 4 * c["mres"],
        "electrical_period_um": 4000.0 / (c["steps_per_rev"] / c["lead_mm"]),
        "runtime_gates": ["precondition arm", "pilot acceptance"],
        "motion_generator": "ESP32 RMT, same buffered profile as the 2/5/8 mm/s test",
    })
    if not args.check_only:
        SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
        write_targets(SEQUENCE_DIR / "e1_targets_grid192.csv", c, targets)
        write_actions(SEQUENCE_DIR / "e1_measurement_plan.csv", measurement)
        (SEQUENCE_DIR / "e1_plan_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        write_firmware_header(c, digest, targets)
    print(
        f"E1 valid: 3 x 101 targets, {summary['window_measurement_stops']} window stops, "
        f"{summary['pilot_stops']} pilot stops; estimated recording "
        f"{summary['estimated_recording_with_pilot_budget_min']:.1f} min"
    )
    print(
        f"RMT command {c['velocity_mm_s']:.1f} mm/s; 0.12 mm grid = "
        f"{c['grid_spacing_usteps']} microsteps = {c['motor_cycles_per_spacing']} motor cycles; "
        f"AquadB ceiling {summary['aquadb_velocity_limit_mm_s']:.1f} mm/s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
