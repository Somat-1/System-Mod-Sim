#!/usr/bin/env python3
"""Generate, validate, and export the E3 global-map sequence.

E3 repeats the E2 design on a different day at *different* positions, which is
what makes it a test of the map rather than of its repeatability. Two things
therefore differ from e2_sequence.py and nothing else does:

  * a different seed, so the stratified draw lands elsewhere; and
  * an explicit exclusion of the E2 targets.

The second is not redundant. Each stratum holds only about seven grid points,
so a fresh seed on its own still reproduces roughly one E2 target in seven -
some fifty of them - and those contribute nothing to a test at new positions.

Targets stay on the 64-microstep grid, as in E1/E2/E4. That is one complete
motor electrical cycle - four full steps at mres = 16 - so every target is
commanded at the same point in the microstep table. The motor's own microstep
error is then a constant offset common to all targets instead of a varying
term aliased into the screw harmonics, which is what keeps the geometry
identifiable. It also keeps this session on the same electrical phase as E2,
so the two are directly comparable.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
import numpy as np

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "e3_config.json"
SEQUENCE_DIR = ROOT / "sequence"
FIRMWARE_HEADER = ROOT / "firmware" / "e3_global_map" / "e3_config.h"

@dataclass(frozen=True)
class Action:
    order: int
    state: str
    action: str
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
    if c["motor_cycle_usteps"] != 4 * c["mres"]:
        raise ValueError("motor cycle must equal four full steps")
    if c["screw_rev_usteps"] != c["steps_per_rev"] * c["mres"]:
        raise ValueError("screw revolution disagrees with mechanics")
    if c["settle_s"] < c["t_hold_s"] + 0.5:
        raise ValueError("settle time is shorter than hold transition + 0.5 s")
    if c["rmt_ring_buffers"] <= c["rmt_queue_depth"]:
        raise ValueError("RMT ring must exceed queue depth")
    if c["watchdog_projection_min"] >= c["time_limit_min"]:
        raise ValueError("watchdog threshold must precede hard limit")
    for key in ("lower_approach_mm", "upper_turn_mm", "ref_preapproach_mm", "x_ref_mm"):
        if not c["soft_min_mm"] <= c[key] <= c["soft_max_mm"]:
            raise ValueError(f"{key} outside proven soft limits")

def load_iso_targets(c: dict[str, Any]) -> tuple[np.ndarray, Path]:
    source = (ROOT / c["iso_target_source"]).resolve()
    if not source.exists():
        raise FileNotFoundError(f"E0 exclusion source not found: {source}")
    with source.open(newline="", encoding="utf-8") as handle:
        values = [round(float(row["cmd_mm"]) * c["ustep_per_mm"])
                  for row in csv.DictReader(handle)]
    if not values:
        raise ValueError("E0 exclusion source contains no targets")
    return np.asarray(values, dtype=np.int64), source

def load_prior_targets(c: dict[str, Any]) -> tuple[np.ndarray, Path]:
    """The E2 target list, to be kept clear of."""
    source = (ROOT / c["e2_target_source"]).resolve()
    if not source.exists():
        raise FileNotFoundError(f"E2 exclusion source not found: {source}")
    with source.open(newline="", encoding="utf-8") as handle:
        values = [int(row["target_usteps"]) for row in csv.DictReader(handle)]
    if len(values) != c["target_count"]:
        raise ValueError(f"E2 source has {len(values)} targets, expected "
                         f"{c['target_count']}")
    return np.asarray(values, dtype=np.int64), source


def candidate_pools(c: dict[str, Any], iso_u: np.ndarray,
                    prior_u: np.ndarray) -> list[np.ndarray]:
    grid, n = int(c["motor_cycle_usteps"]), int(c["target_count"])
    lo = round(c["target_min_mm"] * c["ustep_per_mm"])
    hi = round(c["target_max_mm"] * c["ustep_per_mm"])
    exclusion = int(c["iso_exclusion_grid_points"]) * grid
    prior_exclusion = int(c["e2_exclusion_grid_points"]) * grid
    edges = np.linspace(lo, hi, n + 1)
    pools = []
    for k in range(n):
        g0 = int(np.ceil(edges[k] / grid))
        g1 = int(np.ceil(edges[k + 1] / grid)) - (0 if k == n - 1 else 1)
        candidates = np.arange(g0, g1 + 1, dtype=np.int64) * grid
        candidates = candidates[
            np.min(np.abs(candidates[:, None] - iso_u[None, :]), axis=1) > exclusion]
        if not len(candidates):
            raise ValueError(f"stratum {k + 1} empty after E0 exclusion")
        kept = candidates[
            np.min(np.abs(candidates[:, None] - prior_u[None, :]), axis=1)
            > prior_exclusion]
        # A stratum with nothing left after both exclusions would otherwise
        # silently reopen to E2's own position. Fail instead: an empty stratum
        # means the two exclusions are too wide for the stratum width, which
        # is a design question, not something to paper over here.
        if not len(kept):
            raise ValueError(
                f"stratum {k + 1} empty after the E2 exclusion; "
                f"{len(candidates)} candidate(s) survived E0 alone")
        pools.append(kept)
    return pools

def generate_targets(c: dict[str, Any], iso_u: np.ndarray,
                     prior_u: np.ndarray) -> tuple[np.ndarray, int]:
    """Restricted stratified design; increment seed until phase-balanced."""
    pools = candidate_pools(c, iso_u, prior_u)
    rev, grid = int(c["screw_rev_usteps"]), int(c["motor_cycle_usteps"])
    seed = int(c["seed_start"])
    while seed < 2**32 - 1:
        rng = np.random.default_rng(seed)
        targets = np.asarray([rng.choice(pool) for pool in pools], dtype=np.int64)
        phase = np.bincount((targets % rev) // grid, minlength=rev // grid)
        if phase.min() >= int(c["phase_min_count"]):
            return np.sort(targets), seed
        seed += 1
    raise RuntimeError("no phase-balanced seed found")

class PlanBuilder:
    def __init__(self, c: dict[str, Any]):
        self.c, self.rows = c, []

    def add(self, state: str, action: str, *, direction: str = "",
            index: int | None = None, target_u: int | None = None,
            settle_s: float = 0.0, record_s: float = 0.0,
            record: bool = False, note: str = "") -> None:
        mm = None if target_u is None else target_u / self.c["ustep_per_mm"]
        phase = None if target_u is None else (
            target_u % self.c["screw_rev_usteps"]) // self.c["motor_cycle_usteps"]
        self.rows.append(Action(len(self.rows) + 1, state, action, direction,
                                index, target_u, mm, phase, settle_s,
                                record_s, record, note))

    def transition(self, state: str, note: str = "") -> None:
        self.add(state, "state_enter", note=note)

    def move_mm(self, state: str, target_mm: float, note: str) -> None:
        self.add(state, "move", target_u=round(target_mm * self.c["ustep_per_mm"]),
                 note=note)

    def stop(self, state: str, direction: str, index: int,
             target_u: int, note: str = "") -> None:
        self.add(state, "stop", direction=direction, index=index,
                 target_u=int(target_u), settle_s=self.c["settle_s"],
                 record_s=self.c["record_s"], record=True, note=note)

    def reference(self, state: str, index: int) -> None:
        self.move_mm(state, self.c["ref_preapproach_mm"], "fixed-side pre-approach")
        self.stop(state, "up", index,
                  round(self.c["x_ref_mm"] * self.c["ustep_per_mm"]),
                  "reference visit")

def build_plan(c: dict[str, Any], targets: np.ndarray) -> list[Action]:
    b = PlanBuilder(c)
    b.transition("INIT", "homed, warm, at x_ref, interferometer recording")
    b.add("INIT", "hardware_checks", note="driver, switches, coordinate arm, header")
    b.transition("DATUM")
    for visit in range(1, int(c["datum_visits"]) + 1):
        b.reference("DATUM", visit)
    b.transition("UP", "ascending global-map pass")
    b.move_mm("UP", c["lower_approach_mm"], "safe turn; requested -1 mm is outside limits")
    for index, target in enumerate(targets, 1):
        b.stop("UP", "up", index, int(target))
    b.move_mm("UP", c["upper_turn_mm"], "safe turn; requested 123 mm is outside limits")
    b.transition("TURN_REF")
    b.reference("TURN_REF", 1)
    b.move_mm("TURN_REF", c["upper_turn_mm"], "return to upper turn point")
    b.add("TURN_REF", "time_watchdog", note="abort down pass above 58 min projection")
    b.transition("DOWN", "descending mirror pass")
    for index in range(len(targets), 0, -1):
        b.stop("DOWN", "down", index, int(targets[index - 1]))
    b.transition("END_REF")
    for visit in range(1, int(c["end_ref_visits"]) + 1):
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

def validate_plan(c: dict[str, Any], targets: np.ndarray, iso_u: np.ndarray,
                  prior_u: np.ndarray, seed: int,
                  rows: list[Action]) -> dict[str, Any]:
    grid, rev = int(c["motor_cycle_usteps"]), int(c["screw_rev_usteps"])
    phase = np.bincount((targets % rev) // grid, minlength=rev // grid)
    spacing = np.diff(targets)
    planned = [r.target_usteps for r in rows if r.target_usteps is not None]
    lo_u = round(c["soft_min_mm"] * c["ustep_per_mm"])
    hi_u = round(c["soft_max_mm"] * c["ustep_per_mm"])
    exclusion = int(c["iso_exclusion_grid_points"]) * grid
    prior_exclusion = int(c["e2_exclusion_grid_points"]) * grid
    to_prior = np.min(np.abs(targets[:, None] - prior_u[None, :]), axis=1)
    stops = [r for r in rows if r.action == "stop" and r.state in ("UP", "DOWN")]
    refs = [r for r in rows if r.action == "stop" and
            r.state in ("DATUM", "TURN_REF", "END_REF")]
    estimate = plan_time_s(rows, c)
    checks = {
        "target_count": len(targets) == c["target_count"],
        "unique_and_increasing": len(np.unique(targets)) == len(targets) and bool(np.all(spacing > 0)),
        "targets_on_64_ustep_grid": bool(np.all(targets % grid == 0)),
        "phase_bins_at_least_five": int(phase.min()) >= c["phase_min_count"],
        "no_e0_target_within_one_grid_point": bool(
            np.min(np.abs(targets[:, None] - iso_u[None, :])) > exclusion),
        "no_e2_target_within_one_grid_point": bool(to_prior.min() > prior_exclusion),
        "no_target_shared_with_e2": bool(not np.intersect1d(targets, prior_u).size),
        "seed_differs_from_e2": seed != 5315,
        "all_motion_inside_proven_soft_limits": min(planned) >= lo_u and max(planned) <= hi_u,
        "two_directional_stops_per_target": len(stops) == 2 * len(targets),
        "seven_reference_stops": len(refs) == 7,
        "estimated_recording_below_60_min": estimate <= c["time_limit_min"] * 60,
        "aquadb_velocity_limit_respected": c["velocity_mm_s"] <= 1000 / c["aquadb_clock_ns"],
    }
    if not all(checks.values()):
        raise ValueError("failed checks: " + ", ".join(k for k, v in checks.items() if not v))
    return {
        "experiment": c["experiment"], "seed_start": c["seed_start"],
        "seed_used": seed, "targets": len(targets), "analysis_stops": len(stops),
        "reference_stops": len(refs), "phase_counts": phase.tolist(),
        "phase_count_min": int(phase.min()), "phase_count_max": int(phase.max()),
        "target_min_mm_actual": float(targets.min() / c["ustep_per_mm"]),
        "target_max_mm_actual": float(targets.max() / c["ustep_per_mm"]),
        "spacing_min_mm": float(spacing.min() / c["ustep_per_mm"]),
        "spacing_max_mm": float(spacing.max() / c["ustep_per_mm"]),
        "spacing_mean_mm": float(spacing.mean() / c["ustep_per_mm"]),
        "distance_to_nearest_e2_target_mm": {
            "min": float(to_prior.min() / c["ustep_per_mm"]),
            "median": float(np.median(to_prior) / c["ustep_per_mm"]),
            "max": float(to_prior.max() / c["ustep_per_mm"]),
        },
        "e2_targets_reproduced": int(np.intersect1d(targets, prior_u).size),
        "estimated_recording_min": estimate / 60,
        "requested_turn_points_mm": [-1.0, 123.0],
        "implemented_safe_turn_points_mm": [c["lower_approach_mm"], c["upper_turn_mm"]],
        "endpoint_clearance_note": "Nominal -1/123 mm exceed proven E1 limits.",
        "checks": checks,
    }

def write_targets(path: Path, c: dict[str, Any], targets: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "target_usteps", "target_mm", "screw_phase_bin"])
        for index, target in enumerate(targets, 1):
            phase = int((target % c["screw_rev_usteps"]) // c["motor_cycle_usteps"])
            writer.writerow([index, int(target), f"{target / c['ustep_per_mm']:.6f}", phase])

def write_actions(path: Path, rows: list[Action]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Action.__dataclass_fields__))
        writer.writeheader()
        for row in rows:
            values = asdict(row)
            if values["target_mm"] is not None:
                values["target_mm"] = f"{values['target_mm']:.6f}"
            writer.writerow(values)

def write_header(c: dict[str, Any], digest: str, targets: np.ndarray, seed: int) -> None:
    values = [int(value) for value in targets]
    array = ",\n".join("  " + ", ".join(map(str, values[i:i + 12]))
                       for i in range(0, len(values), 12))
    nominal_spacing = (c["target_max_mm"] - c["target_min_mm"]) / c["target_count"]
    stop_ms = round((move_time_s(nominal_spacing, c) + c["settle_s"] + c["record_s"]) * 1000)
    lines = [
        "/* GENERATED by e3_sequence.py from e3_config.json. Do not edit. */",
        "#pragma once", "#include <stdint.h>", "",
        f'constexpr char CFG_HASH[] = "{digest[:12]}";',
        f"constexpr uint32_t CFG_SEED = {seed}UL;",
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
        f"constexpr int32_t CFG_XREF_U = {round(c['x_ref_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_REF_PRE_U = {round(c['ref_preapproach_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_LOWER_APPROACH_U = {round(c['lower_approach_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_UPPER_TURN_U = {round(c['upper_turn_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {round(c['soft_min_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {round(c['soft_max_mm'] * c['ustep_per_mm'])};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {c['motor_cycle_usteps']};",
        f"constexpr uint16_t CFG_SCREW_REV_U = {c['screw_rev_usteps']};",
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
        f"constexpr uint16_t CFG_TARGET_COUNT = {len(values)};", "",
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
    iso_u, iso_source = load_iso_targets(c)
    prior_u, prior_source = load_prior_targets(c)
    targets, seed = generate_targets(c, iso_u, prior_u)
    rows = build_plan(c, targets)
    summary = validate_plan(c, targets, iso_u, prior_u, seed, rows)
    summary.update({
        "config_sha256": digest, "config": c,
        "e0_exclusion_source": str(iso_source),
        "e2_exclusion_source": str(prior_source),
        "generator": "restricted stratified random on 64-ustep grid, excluding "
                     "the E0 and E2 target neighbourhoods; fixed seed search",
        "motion_generator": "ESP32 RMT buffered 1 MHz trapezoidal/triangular profiles",
    })
    if not args.check_only:
        SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
        write_targets(SEQUENCE_DIR / "e3_targets_seeded_grid64.csv", c, targets)
        write_actions(SEQUENCE_DIR / "e3_measurement_plan.csv", rows)
        (SEQUENCE_DIR / "e3_plan_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        write_header(c, digest, targets, seed)
    nearest = summary["distance_to_nearest_e2_target_mm"]
    print(f"E3 valid: {len(targets)} targets x 2 directions; seed {seed}; "
          f"phase counts {summary['phase_count_min']}..{summary['phase_count_max']}; "
          f"estimated {summary['estimated_recording_min']:.1f} min")
    print(f"New positions: 0 of {len(targets)} shared with E2; nearest E2 target "
          f"{nearest['min'] * 1000:.0f} um away, median {nearest['median'] * 1000:.0f} um")
    print(f"All targets on the {c['motor_cycle_usteps']}-ustep grid: one motor "
          f"electrical cycle, so the microstep error stays a constant")
    print("Safe endpoint adaptation: requested -1/123 mm -> implemented 2.000/120.695 mm")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
