#!/usr/bin/env python3
"""Generate E5 procedure-v2 commands and the ESP32 flash-resident action table.

The host evaluates the frozen E2-only H=2 model.  Firmware only executes the
resulting integer microstep commands; it contains no model code.

    python e5_make_table.py --check-only
    python e5_make_table.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "e5_config.json"
SEQ = HERE / "sequence"
FIRMWARE_HEADER = HERE / "firmware" / "e5_compensated_validation" / "e5_config.h"
sys.path.insert(0, str(HERE))
from e5_geometry import FrozenModel, theta_of  # noqa: E402


@dataclass
class Row:
    order: int
    state: str
    action: str
    block: str = ""
    run: int | None = None
    mode: str = ""
    dir: int | None = None
    target_id: int | None = None
    row: int | None = None
    u_nom: int | None = None
    u_cmd: int | None = None
    x_nom_mm: float | None = None
    G_rel_um: float | None = None
    corr_um: float | None = None
    q_um: float | None = None
    theta_deg: float | None = None
    el_off: int | None = None
    edge: int | None = None
    note: str = ""


TABLE_COLUMNS = ["row", "block", "run", "mode", "dir", "target_id", "u_nom",
                 "u_cmd", "corr_um", "q_um", "theta_deg", "el_off", "edge",
                 "x_nom_mm", "G_rel_um"]
STATE_CODES = {"HOME": 1, "INIT": 2, "R1": 3, "B2A": 4, "R2": 5,
               "B2B": 6, "R3": 7, "B3": 8, "R4": 9, "PARK": 10}
BLOCK_CODES = {"": 0, "R1": 1, "R2": 2, "R3": 3, "R4": 4,
               "B2": 5, "B3": 6}
MODE_CODES = {"": 0, "U": 1, "C": 2}
ACTION_CODES = {"state_enter": 0, "move": 1, "stop": 2}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_config(path: Path) -> tuple[dict, str]:
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def read_usteps(path: Path) -> list[int]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [int(r["target_usteps"]) for r in csv.DictReader(fh)]


def move_time_s(distance_mm: float, c: dict) -> float:
    d, v, a = abs(distance_mm), c["velocity_mm_s"], c["acceleration_mm_s2"]
    return 2 * math.sqrt(d / a) if d <= v * v / a else 2 * v / a + (d - v * v / a) / v


def signed_el_off(u: int) -> int:
    return int((u + 32) % 64) - 32


def b3_targets(c: dict, exclude: np.ndarray) -> tuple[list[int], list[int], int]:
    """One seeded target per stratum; nominal offsets follow -32,-30,...,+30."""
    upm, n = c["ustep_per_mm"], c["b3_target_count"]
    lo, hi = (round(v * upm) for v in c["b3_span_mm"])
    edges = np.linspace(lo, hi, n + 1)
    offsets = [-32 + c["b3_offset_step_usteps"] * k for k in range(n)]
    guard = c["b3_exclusion_grid_points"] * c["motor_cycle_usteps"]
    for seed in range(c["b3_seed_start"], c["b3_seed_start"] + 1000):
        rng, chosen = random.Random(seed), []
        for k, offset in enumerate(offsets):
            a, b = edges[k], edges[k + 1]
            grid = np.arange(math.ceil((a - offset) / 64) * 64,
                             math.floor((b - 1 - offset) / 64) * 64 + 1, 64)
            candidates = [int(g + offset) for g in grid]
            candidates = [u for u in candidates if a <= u < b and
                          np.min(np.abs(exclude - u)) > guard]
            if not candidates:
                break
            chosen.append(rng.choice(candidates))
        if len(chosen) == n:
            return chosen, offsets, seed
    raise RuntimeError("no valid B3 seed found within 1000 attempts")


class Plan:
    def __init__(self, c: dict, model: FrozenModel):
        self.c, self.model, self.rows, self.stops = c, model, [], 0

    def add(self, state: str, action: str, **kw) -> Row:
        row = Row(len(self.rows) + 1, state, action, **kw)
        self.rows.append(row)
        return row

    def move(self, state: str, u_cmd: int, note: str, block: str, run=None) -> None:
        self.add(state, "move", block=block, run=run, u_nom=u_cmd, u_cmd=u_cmd,
                 x_nom_mm=u_cmd / self.c["ustep_per_mm"], note=note)

    def stop(self, state: str, block: str, run: int, mode: str, direction: int,
             target_id: int, u_nom: int, note: str = "") -> None:
        upm = self.c["ustep_per_mm"]
        x = u_nom / upm
        if mode == "C":
            theta_nom = theta_of(np.array([u_nom]))
            g = float(self.model.G_rel(np.array([x]), theta_nom)[0])
            x_c = x - g / 1000.0 - direction * self.model.beta0_um / 2000.0
            u_cmd = int(round(x_c * upm))
            q = (u_cmd / upm - x_c) * 1000.0
        else:
            g, u_cmd, q = None, u_nom, 0.0
        self.stops += 1
        lo, hi = self.c["edge_span_mm"]
        self.add(state, "stop", block=block, run=run, mode=mode,
                 dir=direction, target_id=target_id, row=self.stops,
                 u_nom=u_nom, u_cmd=u_cmd, x_nom_mm=x, G_rel_um=g,
                 corr_um=(u_cmd - u_nom) / upm * 1000.0, q_um=q,
                 theta_deg=float(np.degrees(theta_of(np.array([u_cmd]))[0])),
                 el_off=signed_el_off(u_cmd), edge=int(x < lo or x > hi), note=note)

    def reference_block(self, number: int) -> None:
        state, block = f"R{number}", f"R{number}"
        upm = self.c["ustep_per_mm"]
        approach = round(self.c["reference_preapproach_mm"] * upm)
        self.add(state, "state_enter", block=block, run=number, mode="U")
        for target_id, x in enumerate(self.c["reference_positions_mm"], 1):
            u = round(x * upm)
            self.move(state, u - approach, "reference pre-approach below", block, number)
            self.stop(state, block, number, "U", +1, target_id, u, "reference up")
            self.move(state, u + approach, "reference overrun above", block, number)
            self.stop(state, block, number, "U", -1, target_id, u, "reference down")

    def bidirectional(self, state: str, block: str, run: int, mode: str,
                      targets: list[int]) -> None:
        overrun = round(self.c["overrun_mm"] * self.c["ustep_per_mm"])
        self.move(state, targets[0] - overrun, "lower overrun", block, run)
        for target_id, u in enumerate(targets, 1):
            self.stop(state, block, run, mode, +1, target_id, u)
        self.move(state, targets[-1] + overrun, "upper overrun", block, run)
        for target_id in range(len(targets), 0, -1):
            self.stop(state, block, run, mode, -1, target_id, targets[target_id - 1])
        self.move(state, targets[0] - overrun, "lower overrun", block, run)


def build(c: dict, model: FrozenModel, iso: list[int], b3: list[int],
          b3_repetitions: int | None = None) -> list[Row]:
    p = Plan(c, model)
    p.add("HOME", "state_enter", note="E0 warm-up is external; re-home both ends and centre")
    p.add("INIT", "state_enter", note="print model/table hashes, B3 seed and home MSCNT")
    p.reference_block(1)
    order = c["b2_mode_order"]
    for run in range(1, 6):
        p.add("B2A", "state_enter", block="B2", run=run, mode=order[run - 1])
        p.bidirectional("B2A", "B2", run, order[run - 1], iso)
    p.reference_block(2)
    for run in range(6, 11):
        p.add("B2B", "state_enter", block="B2", run=run, mode=order[run - 1])
        p.bidirectional("B2B", "B2", run, order[run - 1], iso)
    p.reference_block(3)
    repetitions = c["b3_repetitions"] if b3_repetitions is None else b3_repetitions
    for rep in range(1, repetitions + 1):
        p.add("B3", "state_enter", block="B3", run=rep, mode="C")
        p.bidirectional("B3", "B3", rep, "C", b3)
    p.reference_block(4)
    p.add("PARK", "state_enter", note="park and hold after sequence")
    p.move("PARK", round(c["park_mm"] * c["ustep_per_mm"]), "park", "")
    return p.rows


def schedule(rows: list[Row], c: dict) -> float:
    elapsed, position = 0.0, c["park_mm"]
    for row in rows:
        if row.u_cmd is not None:
            target = row.u_cmd / c["ustep_per_mm"]
            elapsed += move_time_s(target - position, c)
            position = target
        if row.action == "stop":
            elapsed += c["settle_s"] + c["record_s"]
    return elapsed


def iso_plan_stops(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [(1 if r["direction"] == "up" else -1, int(r["target_usteps"]))
                for r in csv.DictReader(fh)
                if r["action"] == "stop" and r["state"] == "ISO" and int(r["run"]) == 1]


def validate(c: dict, model: FrozenModel, iso: list[int], b3: list[int],
             offsets: list[int], rows: list[Row], exclude: np.ndarray) -> dict:
    upm = c["ustep_per_mm"]
    stops = [r for r in rows if r.action == "stop"]
    b2 = [r for r in stops if r.block == "B2"]
    b2c = [r for r in b2 if r.mode == "C"]
    b3s = [r for r in stops if r.block == "B3"]
    refs = [r for r in stops if r.block.startswith("R")]
    compensated = b2c + b3s
    duration = schedule(rows, c) / 60.0
    one_rep = schedule(build(c, model, iso, b3, b3_repetitions=1), c) / 60.0
    order = c["b2_mode_order"]
    e0 = iso_plan_stops((HERE / c["iso_plan_source"]).resolve())
    u_runs = sorted({r.run for r in b2 if r.mode == "U"})
    u_match = all([(r.dir, r.u_cmd) for r in b2 if r.run == run] == e0 for run in u_runs)
    diffs = []
    for block_rows in (b2c, b3s):
        for run in sorted({r.run for r in block_rows}):
            up = {r.target_id: r.u_cmd for r in block_rows if r.run == run and r.dir == 1}
            dn = {r.target_id: r.u_cmd for r in block_rows if r.run == run and r.dir == -1}
            diffs.extend((up[k] - dn[k]) / upm * 1000.0 for k in up)
    expected = -model.beta0_um
    commands = [r.u_cmd for r in rows if r.u_cmd is not None]
    lo_u, hi_u = round(c["soft_min_mm"] * upm), round(c["soft_max_mm"] * upm)
    nominal_offsets = [signed_el_off(u) for u in b3]
    commanded_offsets = sorted({r.el_off for r in b3s})
    edge_iso = sorted({r.target_id for r in b2 if r.edge})
    checks = {
        "model_is_E2_only_H2": model.H == 2 and model.data["experiment"] == "E2_global_map",
        "five_U_five_C_ABBA_runs": order == "UCCUUCCUUC" and order.count("U") == order.count("C") == 5,
        "B2_has_480_stops": len(b2) == 480,
        "four_bidirectional_reference_blocks": len(refs) == 16,
        "reference_positions_on_grid_and_uncompensated": all(
            r.mode == "U" and r.u_cmd == r.u_nom and r.u_cmd % 64 == 0 for r in refs),
        "B3_has_requested_stops": len(b3s) == 2 * len(b3) * c["b3_repetitions"],
        "U_runs_match_E0_v4": u_match,
        "max_abs_correction_below_limit": max(abs(r.corr_um) for r in compensated) < c["max_abs_correction_um"],
        "quantisation_within_half_ustep": max(abs(r.q_um) for r in compensated) <= 1000 / upm / 2 + 1e-9,
        "up_minus_down_sign_and_size": all(abs(v - expected) <= 1000 / upm + 1e-9 for v in diffs),
        "B3_nominal_offsets_cover_full_cycle": sorted(nominal_offsets) == list(range(-32, 32, 2)),
        "B3_targets_inside_5_to_110_mm": min(b3) / upm >= 5 and max(b3) / upm <= 110,
        "B3_targets_clear_E0_E2_E3": all(np.min(np.abs(exclude - u)) > 64 for u in b3),
        "all_motion_inside_soft_limits": min(commands) >= lo_u and max(commands) <= hi_u,
        "full_plan_below_58_min_gate": duration < c["watchdog_gate_min"],
    }
    warnings = {
        "eight_minute_margin_to_58_min_gate": duration <= c["watchdog_gate_min"] - c["recording_margin_min"],
        "all_U_rows_on_64_ustep_grid": all(r.u_cmd % 64 == 0 for r in b2 if r.mode == "U"),
        "B3_commanded_offsets_form_exact_even_cycle": commanded_offsets == list(range(-32, 32, 2)),
    }
    corr = np.array([r.corr_um for r in compensated])
    q = np.array([r.q_um for r in compensated])
    return {
        "checks": checks,
        "procedure_warnings": warnings,
        "warning_explanations": {
            "all_U_rows_on_64_ustep_grid": "Impossible together with using the unchanged E0 v4 targets; E0 v4 contains off-grid targets.",
            "eight_minute_margin_to_58_min_gate": "The full two-repeat B3 plan uses measured E0 motion timing; use the documented one-repeat B3 fallback if margin is mandatory.",
            "B3_commanded_offsets_form_exact_even_cycle": "The nominal B3 targets cover the even offsets exactly; direction-dependent command compensation shifts final u_cmd offsets.",
        },
        "duration_min": duration,
        "duration_one_b3_repetition_min": one_rep,
        "stops": {"total": len(stops), "B2": len(b2), "B2_U": len(b2) - len(b2c),
                  "B2_C": len(b2c), "B3": len(b3s), "reference": len(refs)},
        "correction_um": {"min": float(corr.min()), "max": float(corr.max()),
                          "max_abs": float(np.abs(corr).max())},
        "quantisation_um": {"max_abs": float(np.abs(q).max()),
                            "rms": float(np.sqrt(np.mean(q ** 2)))},
        "up_minus_down_um": {"expected": expected, "min": min(diffs), "max": max(diffs)},
        "B3_nominal_offsets": nominal_offsets,
        "B3_commanded_offsets_observed": commanded_offsets,
        "edge_ISO_target_ids": edge_iso,
        "edge_ISO_target_count": len(edge_iso),
    }


def write_table(path: Path, rows: list[Row]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(TABLE_COLUMNS)
        for row in rows:
            if row.action != "stop":
                continue
            data = asdict(row)
            writer.writerow(["" if data[k] is None else
                             f"{data[k]:.6f}" if isinstance(data[k], float) else data[k]
                             for k in TABLE_COLUMNS])


def write_plan(path: Path, rows: list[Row]) -> None:
    names = [f.name for f in fields(Row)]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        writer.writeheader()
        for row in rows:
            data = asdict(row)
            writer.writerow({k: "" if data[k] is None else
                             f"{data[k]:.6f}" if isinstance(data[k], float) else data[k]
                             for k in names})


def write_header(path: Path, c: dict, config_hash: str, model_hash: str,
                 table_hash: str, b3_seed: int, rows: list[Row]) -> None:
    actions = []
    for row in rows:
        actions.append("  {%d, %d, %d, %d, %d, %d, %d, %d, %ld, %ld}" % (
            ACTION_CODES[row.action], STATE_CODES[row.state], BLOCK_CODES[row.block],
            row.run or 0, MODE_CODES[row.mode], row.dir or 0, row.target_id or 0,
            row.row or 0, row.u_nom or 0, row.u_cmd or 0))
    state_names = ", ".join(['"IDLE"'] +
                            [f'"{name}"' for name, _ in sorted(STATE_CODES.items(), key=lambda kv: kv[1])])
    block_names = ", ".join(f'"{name or "NONE"}"' for name, _ in sorted(BLOCK_CODES.items(), key=lambda kv: kv[1]))
    lines = [
        "/* GENERATED by e5_make_table.py from the frozen v2 command table. */",
        "#pragma once", "#include <stdint.h>", "",
        f'constexpr char CFG_HASH[] = "{config_hash}";',
        f'constexpr char CFG_MODEL_HASH[] = "{model_hash}";',
        f'constexpr char CFG_TABLE_HASH[] = "{table_hash}";',
        f"constexpr uint32_t CFG_B3_SEED = {b3_seed}UL;",
        f"constexpr double CFG_LEAD_MM = {c['lead_mm']:.6f};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {c['steps_per_rev']};",
        f"constexpr uint16_t CFG_MRES = {c['mres']};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {c['run_current_ma']};",
        f"constexpr uint16_t CFG_HOLD_CURRENT_MA = {c['hold_current_ma']};",
        f"constexpr uint8_t CFG_TPOWERDOWN = {c['tpowerdown']};",
        f"constexpr uint8_t CFG_IHOLDDELAY = {c['iholddelay']};",
        f"constexpr double CFG_V_MM_S = {c['velocity_mm_s']:.6f};",
        f"constexpr double CFG_A_MM_S2 = {c['acceleration_mm_s2']:.6f};",
        f"constexpr uint32_t CFG_SETTLE_MS = {round(c['settle_s'] * 1000)}UL;",
        f"constexpr uint32_t CFG_RECORD_MS = {round(c['record_s'] * 1000)}UL;",
        f"constexpr int32_t CFG_PARK_U = {round(c['park_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {round(c['soft_min_mm'] * c['ustep_per_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {round(c['soft_max_mm'] * c['ustep_per_mm'])};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {c['motor_cycle_usteps']};",
        f"constexpr uint16_t CFG_SCREW_REV_U = {c['screw_rev_usteps']};",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(c['time_limit_min'] * 60000)}UL;",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['watchdog_gate_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {c['aquadb_clock_ns']};",
        f"constexpr uint32_t CFG_RMT_HZ = {c['rmt_resolution_hz']}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {c['rmt_chunk_pulses']};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {c['rmt_queue_depth']};",
        f"constexpr uint8_t CFG_RMT_RING = {c['rmt_ring_buffers']};",
        f"constexpr int8_t CFG_SYNC_GPIO = {c['sync_gpio']};", "",
        "enum E5ActionType : uint8_t { E5_STATE = 0, E5_MOVE = 1, E5_STOP = 2 };",
        "struct E5Action { uint8_t type, state, block, run, mode; int8_t direction; "
        "uint8_t target; uint16_t row; int32_t u_nom, u_cmd; };", "",
        f"constexpr uint16_t CFG_ACTION_COUNT = {len(rows)};",
        f"constexpr uint16_t CFG_STOP_COUNT = {sum(r.action == 'stop' for r in rows)};",
        f"constexpr const char *CFG_STATE_NAMES[] = {{{state_names}}};",
        f"constexpr const char *CFG_BLOCK_NAMES[] = {{{block_names}}};", "",
        "constexpr E5Action CFG_ACTIONS[CFG_ACTION_COUNT] = {",
        ",\n".join(actions), "};", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    c, config_hash = load_config(CONFIG)
    model = FrozenModel((HERE / c["model_file"]).resolve())
    iso = read_usteps((HERE / c["iso_target_source"]).resolve())
    exclude = np.array(sorted({u for source in c["b3_exclusion_sources"]
                               for u in read_usteps((HERE / source).resolve())}))
    b3, offsets, seed = b3_targets(c, exclude)
    rows = build(c, model, iso, b3)
    summary = validate(c, model, iso, b3, offsets, rows, exclude)
    summary.update({
        "experiment": c["experiment"], "procedure_version": c["procedure_version"],
        "config_sha256": config_hash, "model_file": c["model_file"],
        "model_sha256": model.sha256, "model": {"selection": model.data["selection"],
        "selection_basis": model.data["selection_basis"], "variant": model.data["variant"],
        "a0_um": model.a0, "beta0_um": model.beta0_um,
        "gain_ppm": model.data["gain_ppm"], "domain_mm": list(model.domain)},
        "b3_seed": seed, "b3_targets_usteps": b3, "b3_offsets_assigned": offsets,
        "config": c,
    })
    failed = [name for name, passed in summary["checks"].items() if not passed]
    if not args.check_only:
        SEQ.mkdir(parents=True, exist_ok=True)
        table_path = SEQ / "e5_command_table.csv"
        write_table(table_path, rows)
        write_plan(SEQ / "e5_measurement_plan.csv", rows)
        with (SEQ / "e5_b3_targets.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["index", "target_usteps", "target_mm", "nominal_offset_usteps"])
            for i, (u, offset) in enumerate(zip(b3, offsets), 1):
                writer.writerow([i, u, f"{u / c['ustep_per_mm']:.6f}", offset])
        summary["table_sha256"] = sha256(table_path)
        (SEQ / "e5_plan_summary.json").write_text(
            json.dumps(summary, indent=2, default=float) + "\n", encoding="utf-8")
        write_header(FIRMWARE_HEADER, c, config_hash, model.sha256,
                     summary["table_sha256"], seed, rows)
    warnings = [name for name, passed in summary["procedure_warnings"].items() if not passed]
    print(f"E5 v2 {'VALID' if not failed else 'INVALID'}: {summary['stops']['total']} stops; "
          f"B2={summary['stops']['B2']}, refs={summary['stops']['reference']}, "
          f"B3={summary['stops']['B3']}; {summary['duration_min']:.1f} min "
          f"(one B3 rep {summary['duration_one_b3_repetition_min']:.1f} min)")
    print(f"model H={model.H} Ks={model.Ks} Kc={model.Kc} {model.sha256[:12]}; "
          f"B3 seed={seed}; correction {summary['correction_um']['min']:+.2f}.."
          f"{summary['correction_um']['max']:+.2f} um; |q|<="
          f"{summary['quantisation_um']['max_abs']:.3f} um")
    if warnings:
        print("PROCEDURE WARNINGS: " + ", ".join(warnings))
    if failed:
        print("FAILED: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
