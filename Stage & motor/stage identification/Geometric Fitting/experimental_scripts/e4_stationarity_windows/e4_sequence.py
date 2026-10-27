#!/usr/bin/env python3
"""Generate, validate, and export the E4 stationarity-window sequence."""
from __future__ import annotations

import argparse, csv, hashlib, json, math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "e4_config.json"
SEQUENCE_DIR = ROOT / "sequence"
FIRMWARE_HEADER = ROOT / "firmware" / "e4_stationarity_windows" / "e4_config.h"

@dataclass(frozen=True)
class Action:
    order: int
    state: str
    action: str
    window: str = ""
    direction: str = ""
    index: int | None = None
    target_usteps: int | None = None
    target_mm: float | None = None
    screw_phase_bin: int | None = None
    settle_s: float = 0.0
    record_s: float = 0.0
    record: bool = False
    note: str = ""

def load_config(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()

def move_time_s(distance_mm: float, c: dict[str, Any]) -> float:
    d, v, a = abs(distance_mm), c["velocity_mm_s"], c["acceleration_mm_s2"]
    return 2 * math.sqrt(d / a) if d <= v * v / a else 2 * v / a + (d - v * v / a) / v

def validate_config(c: dict[str, Any]) -> None:
    expected = c["steps_per_rev"] * c["mres"] / c["lead_mm"]
    if not math.isclose(c["ustep_per_mm"], expected):
        raise ValueError("ustep_per_mm disagrees with mechanics")
    if len(c["window_centres_mm"]) != 7:
        raise ValueError("E4 requires exactly seven windows")
    if c["screw_rev_usteps"] != c["steps_per_rev"] * c["mres"]:
        raise ValueError("screw revolution disagrees with mechanics")
    if c["motor_cycle_usteps"] != 4 * c["mres"]:
        raise ValueError("motor cycle must equal four full steps")
    if round(c["grid_spacing_mm"] * c["ustep_per_mm"]) != c["grid_spacing_usteps"]:
        raise ValueError("grid spacing in mm and microsteps disagree")
    if c["grid_spacing_usteps"] != 3 * c["motor_cycle_usteps"]:
        raise ValueError("0.12 mm spacing must equal three motor cycles")
    intervals = round(c["window_span_mm"] / c["grid_spacing_mm"])
    if intervals + 1 != c["targets_per_window"]:
        raise ValueError("window span/grid does not give 51 targets")
    if c["settle_s"] < c["t_hold_s"] + 0.5:
        raise ValueError("settle time is shorter than hold transition + 0.5 s")
    if c["rmt_ring_buffers"] <= c["rmt_queue_depth"]:
        raise ValueError("RMT ring must exceed queue depth")

def generate_targets(c: dict[str, Any]) -> dict[str, list[int]]:
    upmm, step = c["ustep_per_mm"], c["grid_spacing_usteps"]
    half = c["window_span_mm"] / 2
    result = {}
    for number, centre in enumerate(c["window_centres_mm"], 1):
        start = round((centre - half) * upmm)
        values = [start + i * step for i in range(c["targets_per_window"])]
        result[f"W{number}"] = values
    return result

class PlanBuilder:
    def __init__(self, c: dict[str, Any]):
        self.c, self.rows = c, []

    def add(self, state: str, action: str, *, window: str = "",
            direction: str = "", index: int | None = None,
            target_u: int | None = None, settle_s: float = 0.0,
            record_s: float = 0.0, record: bool = False, note: str = "") -> None:
        mm = None if target_u is None else target_u / self.c["ustep_per_mm"]
        phase = None if target_u is None else (
            target_u % self.c["screw_rev_usteps"]) // self.c["motor_cycle_usteps"]
        self.rows.append(Action(len(self.rows) + 1, state, action, window,
                                direction, index, target_u, mm, phase,
                                settle_s, record_s, record, note))

    def transition(self, state: str, note: str = "") -> None:
        self.add(state, "state_enter", note=note)

    def move(self, state: str, target_u: int, window: str, note: str) -> None:
        self.add(state, "move", window=window, target_u=int(target_u), note=note)

    def stop(self, state: str, window: str, direction: str,
             index: int, target_u: int, note: str = "") -> None:
        self.add(state, "stop", window=window, direction=direction,
                 index=index, target_u=int(target_u), settle_s=self.c["settle_s"],
                 record_s=self.c["record_s"], record=True, note=note)

    def reference(self, state: str, index: int) -> None:
        pre = round(self.c["ref_preapproach_mm"] * self.c["ustep_per_mm"])
        xref = round(self.c["x_ref_mm"] * self.c["ustep_per_mm"])
        self.move(state, pre, "REF", "fixed-side pre-approach")
        self.stop(state, "REF", "up", index, xref, "reference visit")

    def window(self, name: str, targets: list[int]) -> None:
        overrun = round(self.c["overrun_mm"] * self.c["ustep_per_mm"])
        self.transition(name, f"stationarity window {name[1:]}")
        self.move(name, targets[0] - overrun, name, "lower pre-approach")
        for i, target in enumerate(targets, 1):
            self.stop(name, name, "up", i, target)
        self.move(name, targets[-1] + overrun, name, "upper overrun")
        for i in range(len(targets), 0, -1):
            self.stop(name, name, "down", i, targets[i - 1])

def build_plan(c: dict[str, Any], targets: dict[str, list[int]]) -> list[Action]:
    b = PlanBuilder(c)
    b.transition("INIT", "homed, warm, at x_ref, interferometer recording")
    b.add("INIT", "hardware_checks", note="driver, switches, coordinate arm, header")
    b.transition("DATUM")
    for visit in range(1, c["datum_visits"] + 1):
        b.reference("DATUM", visit)
    ref_index = 0
    for number in range(1, 8):
        name = f"W{number}"
        b.window(name, targets[name])
        if number < 7:
            ref_index += 1
            b.transition("REF", f"post-{name}")
            b.reference("REF", ref_index)
        if number == c["watchdog_after_window"]:
            b.add("REF", "time_watchdog",
                  note="abort remaining windows if projected finish exceeds 58 min")
    b.transition("END_REF")
    for visit in range(1, c["end_ref_visits"] + 1):
        b.reference("END_REF", visit)
    b.add("END_REF", "environment_snapshot")
    b.transition("PARK", "remain at x_ref; write summary")
    b.add("PARK", "write_summary")
    b.transition("DONE")
    return b.rows

def plan_time_s(rows: list[Action], c: dict[str, Any]) -> float:
    elapsed, position = 0.0, c["x_ref_mm"]
    for row in rows:
        if row.target_mm is not None:
            elapsed += move_time_s(row.target_mm - position, c)
            position = row.target_mm
        elapsed += row.settle_s + row.record_s
    return elapsed

def validate_plan(c: dict[str, Any], targets: dict[str, list[int]],
                  rows: list[Action]) -> dict[str, Any]:
    grid, rev = c["motor_cycle_usteps"], c["screw_rev_usteps"]
    overrun = round(c["overrun_mm"] * c["ustep_per_mm"])
    lo = round(c["soft_min_mm"] * c["ustep_per_mm"])
    hi = round(c["soft_max_mm"] * c["ustep_per_mm"])
    phases, bounds = {}, {}
    for name, values in targets.items():
        phase = [(value % rev) // grid for value in values]
        phases[name] = phase
        bounds[name] = [values[0] / c["ustep_per_mm"], values[-1] / c["ustep_per_mm"]]
    stops = [r for r in rows if r.action == "stop" and r.state.startswith("W")]
    refs = [r for r in rows if r.action == "stop" and r.window == "REF"]
    planned = [r.target_usteps for r in rows if r.target_usteps is not None]
    estimate = plan_time_s(rows, c)
    checks = {
        "seven_windows": len(targets) == 7,
        "fifty_one_targets_per_window": all(len(v) == 51 for v in targets.values()),
        "all_targets_on_64_ustep_grid": all(v % grid == 0 for values in targets.values() for v in values),
        "constant_192_ustep_spacing": all(all(b - a == c["grid_spacing_usteps"]
            for a, b in zip(values, values[1:])) for values in targets.values()),
        "all_50_screw_phases_once_per_50_intervals": all(len(set(p[:-1])) == 50 for p in phases.values()),
        "endpoint_phase_repeats": all(p[0] == p[-1] for p in phases.values()),
        "window_starts_on_64_ustep_grid": all(v[0] % grid == 0 for v in targets.values()),
        "approaches_and_overruns_inside_soft_limits": all(
            lo <= values[0] - overrun and values[-1] + overrun <= hi
            for values in targets.values()),
        "all_planned_motion_inside_soft_limits": min(planned) >= lo and max(planned) <= hi,
        "seven_hundred_fourteen_analysis_stops": len(stops) == 714,
        "twelve_reference_stops": len(refs) == 12,
        "estimated_recording_below_60_min": estimate <= c["time_limit_min"] * 60,
        "aquadb_velocity_limit_respected": c["velocity_mm_s"] <= 1000 / c["aquadb_clock_ns"],
    }
    if not all(checks.values()):
        raise ValueError("failed checks: " + ", ".join(k for k, v in checks.items() if not v))
    return {
        "experiment": c["experiment"], "windows": 7,
        "window_centres_mm": c["window_centres_mm"], "window_bounds_mm": bounds,
        "targets_per_window": 51, "unique_phases_per_window": 50,
        "analysis_stops": len(stops), "reference_stops": len(refs),
        "total_record_stops": len(stops) + len(refs),
        "estimated_recording_min": estimate / 60,
        "phase_bins_by_window": phases, "checks": checks,
    }

def write_targets(path: Path, c: dict[str, Any], targets: dict[str, list[int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["window", "index", "target_usteps", "target_mm", "screw_phase_bin"])
        for name, values in targets.items():
            for index, target in enumerate(values, 1):
                phase = (target % c["screw_rev_usteps"]) // c["motor_cycle_usteps"]
                writer.writerow([name, index, target,
                                 f"{target / c['ustep_per_mm']:.6f}", phase])

def write_actions(path: Path, rows: list[Action]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Action.__dataclass_fields__))
        writer.writeheader()
        for row in rows:
            values = asdict(row)
            if values["target_mm"] is not None:
                values["target_mm"] = f"{values['target_mm']:.6f}"
            writer.writerow(values)

def write_header(c: dict[str, Any], digest: str,
                 targets: dict[str, list[int]]) -> None:
    def array(values: list[int]) -> str:
        return ",\n".join("  " + ", ".join(map(str, values[i:i + 12]))
                           for i in range(0, len(values), 12))
    stop_ms = round((move_time_s(c["grid_spacing_mm"], c) +
                     c["settle_s"] + c["record_s"]) * 1000)
    lines = [
        "/* GENERATED by e4_sequence.py from e4_config.json. Do not edit. */",
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
        f"constexpr uint32_t CFG_NOMINAL_STOP_MS = {stop_ms};",
        f"constexpr int32_t CFG_OVERRUN_U = {round(c['overrun_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_XREF_U = {round(c['x_ref_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_REF_PRE_U = {round(c['ref_preapproach_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {round(c['soft_min_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {round(c['soft_max_mm'] * c['ustep_per_mm'])};",
        f"constexpr uint16_t CFG_GRID_SPACING_U = {c['grid_spacing_usteps']};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {c['motor_cycle_usteps']};",
        f"constexpr uint16_t CFG_SCREW_REV_U = {c['screw_rev_usteps']};",
        f"constexpr uint8_t CFG_DATUM_VISITS = {c['datum_visits']};",
        f"constexpr uint8_t CFG_END_REF_VISITS = {c['end_ref_visits']};",
        f"constexpr uint8_t CFG_WATCHDOG_AFTER_WINDOW = {c['watchdog_after_window']};",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['watchdog_projection_min'] * 60000)}UL;",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(c['time_limit_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {c['aquadb_clock_ns']};",
        f"constexpr float CFG_REF_DRIFT_WARNING_UM = {c['reference_drift_warning_um']:.6f}F;",
        f"constexpr uint32_t CFG_RMT_HZ = {c['rmt_resolution_hz']}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {c['rmt_chunk_pulses']};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {c['rmt_queue_depth']};",
        f"constexpr uint8_t CFG_RMT_RING = {c['rmt_ring_buffers']};",
        f"constexpr int8_t CFG_SYNC_GPIO = {c['sync_gpio']};",
        "constexpr uint8_t CFG_WINDOW_COUNT = 7;",
        f"constexpr uint8_t CFG_TARGETS_PER_WINDOW = {c['targets_per_window']};",
        f"constexpr uint16_t CFG_ANALYSIS_STOP_COUNT = {7 * 2 * c['targets_per_window']};",
        f"constexpr uint16_t CFG_TOTAL_RECORD_STOP_COUNT = "
        f"{7 * 2 * c['targets_per_window'] + c['datum_visits'] + 6 * c['inter_window_ref_visits'] + c['end_ref_visits']};",
    ]
    for name in targets:
        lines += ["", f"constexpr int32_t CFG_{name}_TARGETS[CFG_TARGETS_PER_WINDOW] = {{",
                  array(targets[name]), "};"]
    lines += ["", "constexpr const int32_t *CFG_WINDOWS[CFG_WINDOW_COUNT] = {",
              "  CFG_W1_TARGETS, CFG_W2_TARGETS, CFG_W3_TARGETS, CFG_W4_TARGETS,",
              "  CFG_W5_TARGETS, CFG_W6_TARGETS, CFG_W7_TARGETS", "};", ""]
    FIRMWARE_HEADER.parent.mkdir(parents=True, exist_ok=True)
    FIRMWARE_HEADER.write_text("\n".join(lines), encoding="utf-8")

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    c, digest = load_config(args.config.resolve())
    validate_config(c)
    targets = generate_targets(c)
    rows = build_plan(c, targets)
    summary = validate_plan(c, targets, rows)
    summary.update({
        "config_sha256": digest, "config": c,
        "target_design": "seven 6 mm windows; 0.12 mm grid; 50 intervals/51 targets",
        "phase_design": "spacing advances three of 50 phase bins; all 50 visited before endpoint repeat",
        "motion_generator": "ESP32 RMT buffered 1 MHz trapezoidal/triangular profiles",
        "runtime_gates": ["precondition arm"],
    })
    if not args.check_only:
        SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
        write_targets(SEQUENCE_DIR / "e4_targets_grid192.csv", c, targets)
        write_actions(SEQUENCE_DIR / "e4_measurement_plan.csv", rows)
        (SEQUENCE_DIR / "e4_plan_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        write_header(c, digest, targets)
    print(f"E4 valid: 7 x 51 targets x 2 directions = {summary['analysis_stops']} stops; "
          f"50 screw phases/window; {summary['reference_stops']} references; "
          f"estimated {summary['estimated_recording_min']:.1f} min")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
