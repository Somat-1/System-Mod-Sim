#!/usr/bin/env python3
"""Generate, validate, and export the RMT ISO 230-2 sequence.

The targets come from place_targets() in ../scripts/iso_design.py - the same
seed-2302 placement the v3 baseline used - and are checked against the TARGET
lines of the v3 controller log, so this run measures v3's positions exactly.
The move pattern is v3's: DATUM x3, then per run an x_ref visit, an ascending
pass from first-target - overrun, a descending pass from last-target + overrun,
and the return to first-target - overrun; then END_REF x3 and PARK.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, re, sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "iso_rmt_config.json"
SEQUENCE_DIR = ROOT / "sequence"
FIRMWARE_HEADER = ROOT / "firmware" / "iso230_2_rmt" / "iso_rmt_config.h"

sys.path.insert(0, str(ROOT.parent / "scripts"))
import iso_design  # noqa: E402

@dataclass(frozen=True)
class Action:
    order: int
    state: str
    action: str
    run: int | None = None
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

def move_time_s(distance_mm: float, c: dict[str, Any]) -> float:
    d, v, a = abs(distance_mm), c["velocity_mm_s"], c["acceleration_mm_s2"]
    return 2 * math.sqrt(d / a) if d <= v * v / a else 2 * v / a + (d - v * v / a) / v

def validate_config(c: dict[str, Any]) -> None:
    expected = c["steps_per_rev"] * c["mres"] / c["lead_mm"]
    if not math.isclose(c["ustep_per_mm"], expected):
        raise ValueError("ustep_per_mm disagrees with mechanics")
    if c["settle_s"] < c["t_hold_s"] + 0.5:
        raise ValueError("settle time is shorter than hold transition + 0.5 s")
    if c["rmt_ring_buffers"] <= c["rmt_queue_depth"]:
        raise ValueError("RMT ring must exceed queue depth")
    if c["watchdog_projection_min"] >= c["time_limit_min"]:
        raise ValueError("watchdog threshold must precede hard limit")
    for key in ("ref_preapproach_mm", "x_ref_mm"):
        if not c["soft_min_mm"] <= c[key] <= c["soft_max_mm"]:
            raise ValueError(f"{key} outside proven soft limits")

def generate_targets(c: dict[str, Any]) -> list[int]:
    """iso_design's placement, in microsteps of the configured MRES."""
    design = {key: c[key] for key in ("lead_mm", "steps_per_rev", "mres", "end_margin_mm",
                                      "overrun_mm", "p_nom_mm", "target_end_quantum_mm",
                                      "seed", "x_ref_mm")}
    plan = iso_design.place_targets(c["target_rom_mm"], design)
    per_u = iso_design.U256 // c["mres"]
    if any(u % per_u for u in plan["targets_u256"]):
        raise ValueError("iso_design target is off the microstep grid")
    return [u // per_u for u in plan["targets_u256"]]

def v3_targets(c: dict[str, Any]) -> list[int] | None:
    """The TARGET lines of the v3 log, u256 -> microsteps; None if absent."""
    path = (ROOT / c["v3_target_log"]).resolve()
    if not path.exists():
        return None
    found = re.findall(r"# TARGET,\d+,(\d+),", path.read_text(encoding="utf-8"))
    per_u = iso_design.U256 // c["mres"]
    n = len(found) // 2 or len(found)   # the log lists them once per header
    return [int(u) // per_u for u in found[:n]]

class PlanBuilder:
    def __init__(self, c: dict[str, Any]):
        self.c, self.rows = c, []

    def add(self, state: str, action: str, *, run: int | None = None,
            direction: str = "", index: int | None = None,
            target_u: int | None = None, settle_s: float = 0.0,
            record_s: float = 0.0, record: bool = False, note: str = "") -> None:
        mm = None if target_u is None else target_u / self.c["ustep_per_mm"]
        self.rows.append(Action(len(self.rows) + 1, state, action, run, direction,
                                index, target_u, mm, settle_s, record_s, record, note))

    def transition(self, state: str, run: int | None = None, note: str = "") -> None:
        self.add(state, "state_enter", run=run, note=note)

    def move(self, state: str, target_u: int, note: str, run: int | None = None) -> None:
        self.add(state, "move", run=run, target_u=int(target_u), note=note)

    def stop(self, state: str, run: int, direction: str, index: int,
             target_u: int, note: str = "") -> None:
        self.add(state, "stop", run=run, direction=direction, index=index,
                 target_u=int(target_u), settle_s=self.c["settle_s"],
                 record_s=self.c["record_s"], record=True, note=note)

    def reference(self, state: str, run: int, index: int) -> None:
        upm = self.c["ustep_per_mm"]
        self.move(state, round(self.c["ref_preapproach_mm"] * upm),
                  "fixed-side pre-approach", run)
        self.stop(state, run, "up", index, round(self.c["x_ref_mm"] * upm),
                  "reference visit")

def build_plan(c: dict[str, Any], targets: list[int]) -> list[Action]:
    b = PlanBuilder(c)
    ov = round(c["overrun_mm"] * c["ustep_per_mm"])
    b.transition("HOME", note="'g'/'h': home both ends, check travel, centre at x_ref")
    b.transition("INIT", note="driver, switches, header")
    b.transition("DATUM")
    for visit in range(1, int(c["datum_visits"]) + 1):
        b.reference("DATUM", 0, visit)
    for run in range(1, int(c["n_runs"]) + 1):
        b.transition("ISO_REF", run)
        b.reference("ISO_REF", run, run)
        b.transition("ISO", run, "ascending then descending pass")
        b.move("ISO", targets[0] - ov, "lower overrun", run)
        for index, target in enumerate(targets, 1):
            b.stop("ISO", run, "up", index, target)
        b.move("ISO", targets[-1] + ov, "upper overrun", run)
        for index in range(len(targets), 0, -1):
            b.stop("ISO", run, "down", index, targets[index - 1])
        b.move("ISO", targets[0] - ov, "lower overrun", run)
        if run < int(c["n_runs"]):
            b.add("ISO", "time_watchdog", run=run,
                  note="skip the remaining runs above the watchdog projection")
    b.transition("END_REF")
    for visit in range(1, int(c["end_ref_visits"]) + 1):
        b.reference("END_REF", int(c["n_runs"]), visit)
    b.add("END_REF", "environment_snapshot")
    b.transition("PARK", note="remain at x_ref; write summary")
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

def validate_plan(c: dict[str, Any], targets: list[int], v3: list[int] | None,
                  rows: list[Action]) -> dict[str, Any]:
    upm = c["ustep_per_mm"]
    planned = [r.target_usteps for r in rows if r.target_usteps is not None]
    lo_u, hi_u = round(c["soft_min_mm"] * upm), round(c["soft_max_mm"] * upm)
    stops = [r for r in rows if r.action == "stop" and r.state == "ISO"]
    refs = [r for r in rows if r.action == "stop" and r.state != "ISO"]
    spacing = [b - a for a, b in zip(targets, targets[1:])]
    phase = sorted((t % c["screw_rev_usteps"]) / c["screw_rev_usteps"] for t in targets)
    gaps = [b - a for a, b in zip(phase, phase[1:] + [phase[0] + 1])]
    estimate = plan_time_s(rows, c)
    n_ref = int(c["datum_visits"]) + int(c["n_runs"]) + int(c["end_ref_visits"])
    checks = {
        "unique_and_increasing": all(s > 0 for s in spacing),
        "matches_v3_targets": v3 is not None and v3 == targets,
        "all_motion_inside_proven_soft_limits": min(planned) >= lo_u and max(planned) <= hi_u,
        "two_directional_stops_per_target_per_run":
            len(stops) == 2 * len(targets) * int(c["n_runs"]),
        "reference_stops": len(refs) == n_ref,
        "estimated_recording_below_watchdog": estimate <= c["watchdog_projection_min"] * 60,
        "aquadb_velocity_limit_respected": c["velocity_mm_s"] <= 1000 / c["aquadb_clock_ns"],
    }
    if not all(checks.values()):
        raise ValueError("failed checks: " + ", ".join(k for k, v in checks.items() if not v))
    return {
        "experiment": c["experiment"], "seed": c["seed"], "targets": len(targets),
        "runs": c["n_runs"], "analysis_stops": len(stops), "reference_stops": len(refs),
        "target_min_mm": targets[0] / upm, "target_max_mm": targets[-1] / upm,
        "spacing_min_mm": min(spacing) / upm, "spacing_max_mm": max(spacing) / upm,
        "spacing_mean_mm": sum(spacing) / len(spacing) / upm,
        "lead_phase_max_gap": max(gaps),
        "estimated_run_min": (plan_time_s(rows, c) - plan_time_s(
            [r for r in rows if r.run in (None, 0) or r.state == "END_REF"], c))
            / int(c["n_runs"]) / 60,
        "estimated_recording_min": estimate / 60,
        "checks": checks,
    }

def write_targets(path: Path, c: dict[str, Any], targets: list[int]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "target_usteps", "target_mm", "lead_phase"])
        for index, target in enumerate(targets, 1):
            writer.writerow([index, target, f"{target / c['ustep_per_mm']:.6f}",
                             f"{(target % c['screw_rev_usteps']) / c['screw_rev_usteps']:.4f}"])

def write_actions(path: Path, rows: list[Action]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Action.__dataclass_fields__))
        writer.writeheader()
        for row in rows:
            values = asdict(row)
            if values["target_mm"] is not None:
                values["target_mm"] = f"{values['target_mm']:.6f}"
            writer.writerow(values)

def write_header(c: dict[str, Any], digest: str, targets: list[int]) -> None:
    upm = c["ustep_per_mm"]
    array = ",\n".join("  " + ", ".join(map(str, targets[i:i + 8]))
                       for i in range(0, len(targets), 8))
    spacing = (targets[-1] - targets[0]) / (len(targets) - 1) / upm
    stop_ms = round((move_time_s(spacing, c) + c["settle_s"] + c["record_s"]) * 1000)
    lines = [
        "/* GENERATED by iso_rmt_sequence.py from iso_rmt_config.json. Do not edit. */",
        "#pragma once", "#include <stdint.h>", "",
        f'constexpr char CFG_HASH[] = "{digest[:12]}";',
        f"constexpr uint32_t CFG_SEED = {c['seed']}UL;",
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
        f"constexpr int32_t CFG_XREF_U = {round(c['x_ref_mm'] * upm)};",
        f"constexpr int32_t CFG_REF_PRE_U = {round(c['ref_preapproach_mm'] * upm)};",
        f"constexpr int32_t CFG_OVERRUN_U = {round(c['overrun_mm'] * upm)};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {round(c['soft_min_mm'] * upm)};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {round(c['soft_max_mm'] * upm)};",
        f"constexpr uint8_t CFG_N_RUNS = {c['n_runs']};",
        f"constexpr uint8_t CFG_DATUM_VISITS = {c['datum_visits']};",
        f"constexpr uint8_t CFG_END_REF_VISITS = {c['end_ref_visits']};",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(c['time_limit_min'] * 60000)}UL;",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['watchdog_projection_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {c['aquadb_clock_ns']};",
        f"constexpr float CFG_REF_DRIFT_WARNING_UM = {c['reference_drift_warning_um']:.6f}F;",
        f"constexpr uint32_t CFG_RMT_HZ = {c['rmt_resolution_hz']}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {c['rmt_chunk_pulses']};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {c['rmt_queue_depth']};",
        f"constexpr uint8_t CFG_RMT_RING = {c['rmt_ring_buffers']};",
        f"constexpr int8_t CFG_SYNC_GPIO = {c['sync_gpio']};",
        f"constexpr uint8_t CFG_TARGET_COUNT = {len(targets)};", "",
        "constexpr int32_t CFG_TARGETS[CFG_TARGET_COUNT] = {", array, "};", "",
    ]
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
    summary = validate_plan(c, targets, v3_targets(c), rows)
    summary.update({
        "config_sha256": digest, "config": c,
        "generator": "iso_design.place_targets (seed 2302), checked against the v3 log",
        "motion_generator": "ESP32 RMT buffered 1 MHz trapezoidal/triangular profiles",
        "telemetry": "no serial output during motion; MSCNT read after each record window",
    })
    if not args.check_only:
        SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
        write_targets(SEQUENCE_DIR / "iso_rmt_targets.csv", c, targets)
        write_actions(SEQUENCE_DIR / "iso_rmt_measurement_plan.csv", rows)
        (SEQUENCE_DIR / "iso_rmt_plan_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        write_header(c, digest, targets)
    print(f"ISO RMT valid: {len(targets)} targets x 2 directions x {c['n_runs']} runs = "
          f"{summary['analysis_stops']} stops + {summary['reference_stops']} references; "
          f"targets match v3: {summary['checks']['matches_v3_targets']}; "
          f"estimated {summary['estimated_recording_min']:.1f} min")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
