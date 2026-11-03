#!/usr/bin/env python3
"""Generate and validate the E5 compensated-validation command table.

Every stop of the session, in execution order, with the step count the
firmware will command. No model is evaluated on the ESP32: the host applies

    x_c = x* - [G(x*) - a0] - (d/2) beta0,     u_c = round(x_c * U)

with G the frozen E2 model (model/geom_model_E2_frozen.json, hash recorded),
beta0 its D2 reversal, d = +1 up / -1 down, and records the quantisation
residual q = u_c/U - x_c. Uncompensated (U) rows command u_nom unchanged.

Session (procedure §3), with E0 v4's ISO cycle exactly:

    B1a  3 x reference                       (DATUM)
    B2   10 ISO runs, order U C C U U C C U U C; each run opens with E0 v4's
         x_ref visit, then 24 up / 24 down stops between 3 mm overruns
    B1b  3 x reference after run 5            (MID_REF)
    B1c  3 x reference                        (END_REF)
    --   recording split (operator gate; stage held at x_ref, no re-home)
    B1d  3 x reference                        (B3_REF)
    B3   32 off-grid targets, compensated, up/down, 2 repetitions
    B1e  3 x reference                        (B3_END_REF)

The split exists because E0 v4's measured cycle is 5.05 min per run, not the
~2.5 min the procedure's 3 s/stop assumed: B2 alone fills one 52 min
recording, so B3 gets a second recording in the same homing rather than
being cut. Set split_recording_before_b3 false to see the single-window plan
fail its time check.

    python3 e5_make_table.py --check-only   # validate, write nothing
    python3 e5_make_table.py                # write sequence/
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
sys.path.insert(0, str(HERE))
from e5_geometry import FrozenModel, theta_of  # noqa: E402


@dataclass
class Row:
    """One firmware action. Stops carry the full command-table columns."""
    order: int
    recording: int
    state: str
    action: str
    block: str = ""
    run: int | None = None
    mode: str = ""
    dir: int | None = None
    target_id: int | None = None
    row: int | None = None          # stop index, logged by the firmware
    u_nom: int | None = None
    u_cmd: int | None = None
    x_nom_mm: float | None = None
    G_rel_um: float | None = None
    corr_um: float | None = None
    q_um: float | None = None
    theta_deg: float | None = None
    el_off: int | None = None
    note: str = ""


TABLE_COLUMNS = ["row", "recording", "block", "run", "mode", "dir", "target_id",
                 "u_nom", "u_cmd", "corr_um", "q_um", "theta_deg", "el_off",
                 "x_nom_mm", "G_rel_um"]


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


# ---- B3 off-grid targets ---------------------------------------------------
def b3_targets(c: dict, exclude: np.ndarray) -> tuple[list[int], list[int], int]:
    """32 targets, one per equal stratum of the span, each at its own offset.

    The offsets -32, -30, ..., +30 usteps from the 64-ustep grid cover the full
    electrical cycle once. They are assigned to strata by a seeded permutation,
    not in x order, so electrical phase is not confounded with position (and
    with the slow error s(x)). The base grid point in each stratum is drawn
    from those whose offset target keeps > 1 grid point from every E0/E2/E3
    target. Seeds increment until every stratum has a candidate.
    """
    upm, n = c["ustep_per_mm"], c["b3_target_count"]
    lo, hi = (round(v * upm) for v in c["b3_span_mm"])
    edges = np.linspace(lo, hi, n + 1)
    offsets = [-32 + c["b3_offset_step_usteps"] * k for k in range(n)]
    guard = c["b3_exclusion_grid_points"] * c["motor_cycle_usteps"]
    for seed in range(c["b3_seed_start"], c["b3_seed_start"] + 1000):
        rng = random.Random(seed)
        perm = offsets[:]; rng.shuffle(perm)
        chosen = []
        for k in range(n):
            a, b = edges[k], edges[k + 1]
            grid = np.arange(math.ceil(a / 64) * 64, math.floor(b / 64) * 64 + 1, 64)
            cand = grid + perm[k]
            cand = cand[(cand >= a) & (cand < b)]
            cand = [int(u) for u in cand if np.min(np.abs(exclude - u)) > guard]
            if not cand:
                break
            chosen.append(rng.choice(cand))
        if len(chosen) == n:
            return chosen, perm, seed
    raise SystemExit("no B3 seed found within 1000 tries")


# ---- plan ------------------------------------------------------------------
class Plan:
    def __init__(self, c: dict, model: FrozenModel):
        self.c, self.m, self.rows, self.recording, self.stops = c, model, [], 1, 0

    def add(self, state: str, action: str, **kw) -> Row:
        r = Row(len(self.rows) + 1, self.recording, state, action, **kw)
        self.rows.append(r)
        return r

    def move(self, state: str, u: int, note: str, block: str, run=None) -> None:
        self.add(state, "move", block=block, run=run, u_nom=u, u_cmd=u,
                 x_nom_mm=u / self.c["ustep_per_mm"], note=note)

    def stop(self, state: str, block: str, run, mode: str, d: int,
             target_id: int, u_nom: int, note: str = "") -> None:
        upm = self.c["ustep_per_mm"]
        x = u_nom / upm
        if mode == "C":
            g = float(self.m.G_rel(np.array([x]), theta_of(np.array([u_nom])))[0])
            x_c = x - g / 1000.0 - d * self.m.beta0_um / 2000.0
            u_cmd = int(round(x_c * upm))
            q = (u_cmd / upm - x_c) * 1000.0
        else:
            g, u_cmd, q = None, u_nom, 0.0
        self.stops += 1
        self.add(state, "stop", block=block, run=run, mode=mode, dir=d,
                 target_id=target_id, row=self.stops, u_nom=u_nom, u_cmd=u_cmd,
                 x_nom_mm=x, G_rel_um=g, corr_um=(u_cmd - u_nom) / upm * 1000.0,
                 q_um=q, theta_deg=float(np.degrees(theta_of(np.array([u_cmd]))[0])),
                 el_off=signed_el_off(u_cmd), note=note)

    def reference(self, state: str, block: str, visits: int, run=None) -> None:
        upm = self.c["ustep_per_mm"]
        self.add(state, "state_enter", block=block, run=run)
        for v in range(1, visits + 1):
            self.move(state, round(self.c["ref_preapproach_mm"] * upm),
                      "fixed-side pre-approach", block, run)
            self.stop(state, block, run, "U", +1, v, round(self.c["x_ref_mm"] * upm),
                      "reference visit")

    def bidirectional(self, state: str, block: str, run, mode: str,
                      targets: list[int]) -> None:
        """E0 v4's pass: overruns from the NOMINAL end targets, up then down."""
        ov = round(self.c["overrun_mm"] * self.c["ustep_per_mm"])
        self.move(state, targets[0] - ov, "lower overrun", block, run)
        for i, u in enumerate(targets, 1):
            self.stop(state, block, run, mode, +1, i, u)
        self.move(state, targets[-1] + ov, "upper overrun", block, run)
        for i in range(len(targets), 0, -1):
            self.stop(state, block, run, mode, -1, i, targets[i - 1])
        self.move(state, targets[0] - ov, "lower overrun", block, run)


def build(c: dict, model: FrozenModel, iso: list[int], b3: list[int]) -> list[Row]:
    p = Plan(c, model)
    p.add("HOME", "state_enter", note="home both ends (E2 routine), centre at x_ref, arm")
    p.add("INIT", "state_enter", note="driver, switches, header: config, model and B3 hashes")
    p.reference("DATUM", "B1a", c["b1_visits"])
    for run, mode in enumerate(c["b2_mode_order"], 1):
        p.reference("ISO_REF", "B2", 1, run)            # E0 v4's per-run x_ref visit
        p.add("ISO", "state_enter", block="B2", run=run, mode=mode)
        p.bidirectional("ISO", "B2", run, mode, iso)
        if run == c["b1b_after_run"]:
            p.reference("MID_REF", "B1b", c["b1_visits"], run)
    p.reference("END_REF", "B1c", c["b1_visits"])
    if c["split_recording_before_b3"]:
        p.add("SPLIT", "operator_gate", note="hold at x_ref; save IDS file 1, start file 2, press c")
        p.recording = 2
    p.reference("B3_REF", "B1d", c["b1_visits"])
    for rep in range(1, c["b3_repetitions"] + 1):
        p.add("B3", "state_enter", block="B3", run=rep, mode="C")
        p.bidirectional("B3", "B3", rep, "C", b3)
    p.reference("B3_END_REF", "B1e", c["b1_visits"])
    p.add("PARK", "state_enter", note="remain at x_ref, motor held; SEQ_DONE then RAM dump")
    return p.rows


def schedule(rows: list[Row], c: dict) -> dict[int, float]:
    """Nominal duration of each recording, in seconds (moves + settle + record)."""
    t, pos, per = {}, c["x_ref_mm"], {}
    for r in rows:
        if r.u_cmd is not None:
            x = r.u_cmd / c["ustep_per_mm"]
            per[r.recording] = per.get(r.recording, 0.0) + move_time_s(x - pos, c)
            pos = x
        if r.action == "stop":
            per[r.recording] = per.get(r.recording, 0.0) + c["settle_s"] + c["record_s"]
    return per


def iso_plan_stops(path: Path) -> list[tuple[int, int, int]]:
    """(run, direction, target_usteps) of every ISO stop in the E0 v4 plan."""
    out = []
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["action"] == "stop" and r["state"] == "ISO":
                out.append((int(r["run"]), 1 if r["direction"] == "up" else -1,
                            int(r["target_usteps"])))
    return out


def validate(c: dict, model: FrozenModel, iso: list[int], b3: list[int],
             perm: list[int], rows: list[Row], exclude: np.ndarray) -> dict:
    upm = c["ustep_per_mm"]
    stops = [r for r in rows if r.action == "stop"]
    b2 = [r for r in stops if r.block == "B2" and r.state == "ISO"]
    b2c = [r for r in b2 if r.mode == "C"]
    b3s = [r for r in stops if r.block == "B3"]
    comp = b2c + b3s
    lo_u, hi_u = round(c["soft_min_mm"] * upm), round(c["soft_max_mm"] * upm)
    commanded = [r.u_cmd for r in rows if r.u_cmd is not None]

    # Up-minus-down command difference per compensated target and run.
    diffs = []
    for block, group in (("B2", b2c), ("B3", b3s)):
        for run in sorted({r.run for r in group}):
            ups = {r.target_id: r.u_cmd for r in group if r.run == run and r.dir == 1}
            dns = {r.target_id: r.u_cmd for r in group if r.run == run and r.dir == -1}
            diffs += [(ups[k] - dns[k]) / upm * 1000.0 for k in ups]
    expected = -model.beta0_um

    # U runs against the E0 v4 plan, run by run.
    e0 = iso_plan_stops((HERE / c["iso_plan_source"]).resolve())
    e0_run = [(d, u) for run, d, u in e0 if run == 1]
    u_runs = sorted({r.run for r in b2 if r.mode == "U"})
    u_match = all([(r.dir, r.u_cmd) for r in b2 if r.run == run] == e0_run for run in u_runs)

    per = schedule(rows, c)
    budget = (c["time_limit_min"] - c["recording_margin_min"]) * 60
    order = c["b2_mode_order"]
    # ABBA: a linear drift over run index cancels when the mean run index of
    # U equals that of C.
    mean_u = np.mean([i for i, m in enumerate(order, 1) if m == "U"])
    mean_c = np.mean([i for i, m in enumerate(order, 1) if m == "C"])
    in_span = model.in_domain(np.array([r.x_nom_mm for r in comp]))
    checks = {
        "iso_targets_are_e0_v4": len(iso) == 24,
        "u_runs_identical_to_e0_v4_commands": u_match,
        "five_u_five_c_runs": order.count("U") == 5 and order.count("C") == 5,
        "abba_drift_balance": bool(abs(mean_u - mean_c) < 0.6),
        "max_abs_correction_below_limit":
            max(abs(r.corr_um) for r in comp) < c["max_abs_correction_um"],
        "quantisation_within_half_ustep":
            max(abs(r.q_um) for r in comp) <= 1000 / upm / 2 + 1e-9,
        "up_minus_down_equals_minus_beta0":
            all(abs(v - expected) <= 1000 / upm + 1e-9 for v in diffs),
        "compensated_targets_inside_model_domain": bool(in_span.all()),
        "approach_direction_never_reversed":
            max(abs(r.corr_um) for r in comp) / 1000 < c["overrun_mm"],
        "all_motion_inside_proven_soft_limits":
            min(commanded) >= lo_u and max(commanded) <= hi_u,
        "b3_offsets_cover_cycle":
            sorted(perm) == list(range(-32, 32, c["b3_offset_step_usteps"])),
        "b3_no_target_within_one_grid_point_of_e0_e2_e3": bool(all(
            np.min(np.abs(exclude - u)) > c["b3_exclusion_grid_points"] * 64 for u in b3)),
        "each_recording_within_budget": all(v <= budget for v in per.values()),
        "aquadb_velocity_limit_respected":
            c["velocity_mm_s"] <= 1000 / c["aquadb_clock_ns"],
    }
    corr_c = np.array([r.corr_um for r in comp]); q_c = np.array([r.q_um for r in comp])
    b3_el = [signed_el_off(u) for u in b3]
    return {
        "checks": checks,
        "recording_min": {str(k): v / 60 for k, v in per.items()},
        "recording_budget_min": budget / 60,
        "single_window_min": sum(per.values()) / 60,
        "stops": {"total": len(stops), "B2": len(b2), "B2_U": len(b2) - len(b2c),
                  "B2_C": len(b2c), "B3": len(b3s),
                  "reference": len([r for r in stops if r.block.startswith("B1")
                                    or r.state == "ISO_REF"])},
        "correction_um": {"min": float(corr_c.min()), "max": float(corr_c.max()),
                          "max_abs": float(np.abs(corr_c).max())},
        "quantisation_um": {"max_abs": float(np.abs(q_c).max()),
                            "rms": float(np.sqrt(np.mean(q_c**2)))},
        "up_minus_down_um": {"expected": expected, "min": min(diffs), "max": max(diffs)},
        "abba_mean_run_index": {"U": float(mean_u), "C": float(mean_c)},
        "b3_commanded_el_off_range": [min(b3_el), max(b3_el)],
    }


def spot_checks(model: FrozenModel, iso: list[int], c: dict) -> list[dict]:
    """Three targets, one within 5 mm of each end, re-derived step by step."""
    upm = c["ustep_per_mm"]
    out = []
    for u in (iso[0], iso[len(iso) // 2], iso[-1]):
        x = u / upm
        theta = float(theta_of(np.array([u]))[0])
        g = float(model.G_rel(np.array([x]), np.array([theta]))[0])
        for d in (+1, -1):
            x_c = x - g / 1000 - d * model.beta0_um / 2000
            out.append({"u_nom": u, "x_mm": x, "theta_deg": math.degrees(theta),
                        "G_rel_um": g, "dir": d, "x_c_mm": x_c,
                        "u_cmd": int(round(x_c * upm)),
                        "corr_um": (round(x_c * upm) - u) / upm * 1000})
    return out


def write_table(path: Path, rows: list[Row]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(TABLE_COLUMNS)
        for r in rows:
            if r.action != "stop":
                continue
            v = asdict(r)
            w.writerow([("" if v[k] is None else
                         f"{v[k]:.6f}" if isinstance(v[k], float) else v[k])
                        for k in TABLE_COLUMNS])


def write_plan(path: Path, rows: list[Row]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        names = [f.name for f in fields(Row)]
        w = csv.DictWriter(fh, fieldnames=names)
        w.writeheader()
        for r in rows:
            v = asdict(r)
            w.writerow({k: ("" if v[k] is None else
                            f"{v[k]:.6f}" if isinstance(v[k], float) else v[k])
                        for k in names})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()
    c, digest = load_config(CONFIG)
    model = FrozenModel((HERE / c["model_file"]).resolve())
    iso = read_usteps((HERE / c["iso_target_source"]).resolve())
    exclude = np.array(sorted({u for src in c["b3_exclusion_sources"]
                               for u in read_usteps((HERE / src).resolve())}))
    b3, perm, seed = b3_targets(c, exclude)
    rows = build(c, model, iso, b3)
    summary = validate(c, model, iso, b3, perm, rows, exclude)
    summary.update({
        "experiment": c["experiment"], "config_sha256": digest,
        "model_file": c["model_file"], "model_sha256": model.sha256,
        "model": {"selection": model.data["selection"], "variant": model.data["variant"],
                  "a0_um": model.a0, "beta0_um": model.beta0_um,
                  "gain_ppm": model.data["gain_ppm"], "domain_mm": list(model.domain)},
        "b3_seed": seed, "b3_targets_usteps": b3, "b3_offsets_assigned": perm,
        "spot_checks": spot_checks(model, iso, c), "config": c,
    })
    failed = [k for k, v in summary["checks"].items() if not v]
    if not args.check_only:
        SEQ.mkdir(parents=True, exist_ok=True)
        write_table(SEQ / "e5_command_table.csv", rows)
        write_plan(SEQ / "e5_measurement_plan.csv", rows)
        with (SEQ / "e5_b3_targets.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["index", "target_usteps", "target_mm", "assigned_offset_usteps",
                        "screw_phase_deg"])
            for i, (u, o) in enumerate(zip(b3, perm), 1):
                w.writerow([i, u, f"{u / c['ustep_per_mm']:.6f}", o,
                            f"{math.degrees(theta_of(np.array([u]))[0]):.2f}"])
        (SEQ / "e5_plan_summary.json").write_text(
            json.dumps(summary, indent=2, default=float) + "\n", encoding="utf-8")
    rec = ", ".join(f"R{k} {v:.1f} min" for k, v in summary["recording_min"].items())
    print(f"E5 {'valid' if not failed else 'INVALID'}: {summary['stops']['total']} stops "
          f"(B2 {summary['stops']['B2']}, B3 {summary['stops']['B3']}); {rec} "
          f"(budget {summary['recording_budget_min']:.0f} min each; single window "
          f"{summary['single_window_min']:.1f} min)")
    print(f"model {model.sha256[:12]}  B3 seed {seed}  corrections "
          f"{summary['correction_um']['min']:+.2f}..{summary['correction_um']['max']:+.2f} um  "
          f"|q| <= {summary['quantisation_um']['max_abs']:.3f} um  up-down "
          f"{summary['up_minus_down_um']['min']:.3f}..{summary['up_minus_down_um']['max']:.3f} um")
    if failed:
        print("FAILED: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
