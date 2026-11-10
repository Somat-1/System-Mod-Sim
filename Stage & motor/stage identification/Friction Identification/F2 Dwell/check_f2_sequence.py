#!/usr/bin/env python3
"""Check the written F2 artefacts against the pre-run list in procedure 5.

This reads only what will be taken to the bench -- the command table, the
trial order, the plan and table summaries and the generated firmware header
-- and re-derives each check from them. It does not import either builder,
so a mistake in a builder cannot pass its own check.

    python3 check_f2_sequence.py

Exit status is 0 only if every check passes.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f2_config.json"
SEQ = HERE / "sequence"
TABLE = SEQ / "f2_command_table.csv"
ORDER = SEQ / "f2_trial_order.csv"
TABLE_SUMMARY = SEQ / "f2_table_summary.json"
PLAN_SUMMARY = SEQ / "f2_plan_summary.json"
HEADER = HERE / "firmware" / "f2_dwell" / "f2_config.h"

USTEP_PER_MM = 25600
MOTOR_CYCLE_U = 1024
SOFT_MIN_U = round(2.0 * USTEP_PER_MM)
SOFT_MAX_U = round(120.695 * USTEP_PER_MM)
RMT_MIN_RATE_HZ = 1e6 / 30000
TYPES = {"move": 0, "ramp": 1, "step": 2, "dwell": 3, "mark": 4}
KINDS = ["transit", "setup", "settle", "dwell", "measure", "hold", "ref"]
BLOCKS = ["R1", "B", "A", "L", "R3"]
SUBS = ["ref", "R", "C", "setup"]


def read_table(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["row"] = int(r["row"])
        r["trial"] = int(r["trial"])
        r["rep"] = int(r["rep"])
        r["u_target"] = int(r["u_target"])
        r["td_s"] = float(r["td_s"])
        r["rate_hz"] = float(r["rate_hz"])
        r["dwell_s"] = float(r["dwell_s"])
        r["t0_s"] = float(r["t0_s"])
        r["t1_s"] = float(r["t1_s"])
    return rows


def header_segments(text: str) -> list[tuple]:
    body = text.split("CFG_SEGS[CFG_SEG_COUNT] = {", 1)[1]
    out = []
    for line in body.splitlines():
        m = re.match(r"\s*\{(.+?)\},?\s*$", line)
        if not m:
            continue
        out.append(tuple(int(v.strip().rstrip("UL")) for v in m.group(1).split(",")))
    return out


def header_value(text: str, name: str) -> str:
    m = re.search(rf"\b{name}(?:\[\])?\s*=\s*([^;]+);", text)
    return m.group(1).strip().strip('"') if m else ""


def pearson(xs: list, ys: list) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def main() -> int:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    rows = read_table(TABLE)
    summary = json.loads(TABLE_SUMMARY.read_text(encoding="utf-8"))
    variant = summary["variant"]
    blocks = c["reduced_blocks"] if variant == "reduced" else c["blocks"]
    text = HEADER.read_text(encoding="utf-8")
    segs = header_segments(text)
    with ORDER.open(newline="", encoding="utf-8") as fh:
        order_rows = list(csv.DictReader(fh))

    results: dict[str, bool] = {}
    notes: list[str] = []

    def check(name: str, value: bool, note: str = "") -> None:
        results[name] = bool(value)
        if note:
            notes.append(f"  {name}: {note}")

    # -- 5.1 seed logged, order printed, coverage as planned ---------------
    check("seed_in_config", isinstance(c["shuffle_seed"], int))
    check("seed_in_header",
          header_value(text, "CFG_SHUFFLE_SEED") == f"{c['shuffle_seed']}ULL")
    check("seed_in_summary", summary["shuffle_seed"] == c["shuffle_seed"])
    check("trial_order_file_complete", len(order_rows) == summary["trials"])

    planned = {f"{b['name']}.{d['s']:g}s.{m}": d["reps"]
               for b in blocks for d in b["dwells"] for m in b["modes"]}
    seen: dict[str, int] = {}
    for o in order_rows:
        key = f"{o['block']}.{float(o['td_s']):g}s.{o['mode']}"
        seen[key] = seen.get(key, 0) + 1
    check("coverage_matches_grid", seen == planned,
          f"{len(seen)} (T_d, mode) cells, {sum(seen.values())} trials")

    per_rep: dict[tuple, list] = {}
    for o in order_rows:
        per_rep.setdefault((o["block"], o["rep"]), []).append(
            (float(o["td_s"]), o["mode"]))
    check("each_cell_once_per_repetition",
          all(len(v) == len(set(v)) for v in per_rep.values()))

    shuffled_blocks = {b["name"] for b in blocks if b.get("shuffle", True)
                       and len(b["dwells"]) > 1}
    actually_shuffled = any(
        v != sorted(v) for (block, _rep), v in per_rep.items()
        if block in shuffled_blocks)
    check("order_is_shuffled", actually_shuffled or not shuffled_blocks)

    starts = {}
    for r in rows:
        if r["trial"] and r["trial"] not in starts:
            starts[r["trial"]] = r
    corr = {}
    for block in dict.fromkeys(r["block"] for r in starts.values()):
        rb = [r for r in starts.values() if r["block"] == block]
        if len(rb) >= 3:
            corr[block] = pearson([math.log10(r["td_s"]) for r in rb],
                                  [r["t0_s"] for r in rb])
    check("dwell_time_decorrelated",
          all(abs(v) <= c["max_dwell_time_correlation"] for v in corr.values()),
          ", ".join(f"{b} r={v:+.3f}" for b, v in corr.items()))

    # -- 5.2 the dwell itself ---------------------------------------------
    dwells = [r for r in rows if r["kind"] == "dwell"]
    check("one_dwell_per_trial", len(dwells) == summary["trials"])
    check("dwell_durations_match_grid",
          all(abs(r["dwell_s"] - r["td_s"]) < 1e-9 for r in dwells))
    check("dwells_whole_milliseconds",
          all(abs(r["dwell_s"] * 1000 - round(r["dwell_s"] * 1000)) < 1e-6
              for r in rows if r["seg"] in ("dwell", "mark")))
    check("dwells_motionless",
          all(r["u_target"] == rows[r["row"] - 1]["u_target"] for r in dwells))
    check("dwell_tolerance_in_header",
          header_value(text, "CFG_DWELL_TOL_MS")
          == f"{round(c['dwell_timing_tolerance_ms'])}UL",
          f"{c['dwell_timing_tolerance_ms']:g} ms")
    check("longest_dwell_in_header",
          header_value(text, "CFG_LONGEST_DWELL_MS")
          == f"{round(summary['longest_dwell_s'] * 1000)}UL",
          f"{summary['longest_dwell_s']:g} s")

    # -- 5.3 wipe from below; R down, C up --------------------------------
    last, wiped, ramps = 0, 0, 0
    prev_u = rows[0]["u_target"]
    direction_ok = True
    for r in rows:
        d = (r["u_target"] > prev_u) - (r["u_target"] < prev_u)
        if r["kind"] == "dwell":
            if last != 1:
                direction_ok = False
            wiped += 1
        if r["kind"] == "measure":
            want = -last if r["sub"] == "R" else last
            if d != want:
                direction_ok = False
            ramps += 1
        if d:
            last = d
        prev_u = r["u_target"]
    check("all_dwells_follow_a_wipe_from_below", direction_ok and wiped == len(dwells),
          f"{wiped} dwells, {ramps} post-dwell ramps")
    check("post_dwell_ramp_directions", direction_ok)
    check("ramp_count_matches_trials", ramps == summary["trials"])
    check("ramp_lengths_equal",
          len({abs(r["u_target"] - rows[r["row"] - 1]["u_target"])
               for r in rows if r["kind"] == "measure"}) == 1,
          f"{summary['ramp_usteps']} usteps = {summary['ramp_um']:.4f} um")
    check("ramp_rate_on_timed_path",
          all(r["rate_hz"] < RMT_MIN_RATE_HZ for r in rows
              if r["kind"] == "measure"),
          f"{summary['ramp_rate_hz']:.2f} Hz < {RMT_MIN_RATE_HZ:.2f} Hz")

    # -- 5.4 duration and recording ---------------------------------------
    duration_min = rows[-1]["t1_s"] / 60.0
    check("inside_gate", duration_min <= c["gate_min"],
          f"{duration_min:.2f} min of {c['gate_min']:.0f} min")
    check("inside_scope_record", duration_min <= c["scope_record_min"])
    check("gate_plus_longest_dwell_inside_record",
          c["gate_min"] + summary["longest_dwell_s"] / 60.0
          <= c["scope_record_min"])
    check("table_times_contiguous",
          all(abs(b["t0_s"] - a["t1_s"]) < 2e-3 for a, b in zip(rows, rows[1:])))

    # -- bench safety and geometry ----------------------------------------
    check("inside_soft_limits",
          all(SOFT_MIN_U <= r["u_target"] <= SOFT_MAX_U for r in rows))
    check("positions_on_electrical_grid",
          all(round(c["positions_mm"][b["position"]] * USTEP_PER_MM)
              % MOTOR_CYCLE_U == 0 for b in blocks))
    check("references_present_at_both_ends",
          (rows[0]["block"], rows[-1]["block"]) == ("R1", "R3")
          or variant == "reduced" and not c["reduced_includes_references"])
    check("irun_equals_ihold_in_header",
          header_value(text, "CFG_RUN_CURRENT_MA")
          == header_value(text, "CFG_HOLD_CURRENT_MA")
          and header_value(text, "CFG_TPOWERDOWN") == "0"
          and header_value(text, "CFG_IHOLDDELAY") == "0")
    check("interpolation_and_mres_fixed",
          header_value(text, "CFG_MRES") == "256")

    # -- the header is the table -------------------------------------------
    check("header_segment_count",
          header_value(text, "CFG_SEG_COUNT") == str(len(rows))
          and len(segs) == len(rows))
    mismatch = next(
        (r["row"] for r, s in zip(rows, segs)
         if s != (TYPES[r["seg"]], KINDS.index(r["kind"]),
                  BLOCKS.index(r["block"]), SUBS.index(r["sub"]), r["row"],
                  r["trial"], r["u_target"], round(r["rate_hz"] * 1000),
                  round(r["dwell_s"] * 1000))),
        None)
    check("header_matches_table", mismatch is None,
          "row for row" if mismatch is None else f"first mismatch at row {mismatch}")
    check("header_trial_count",
          header_value(text, "CFG_TRIAL_COUNT") == str(summary["trials"]))
    check("header_dwell_windows",
          header_value(text, "CFG_DWELL_WINDOWS") == str(len(dwells)))
    check("header_variant_matches", header_value(text, "CFG_VARIANT") == variant)
    check("header_marks_match",
          header_value(text, "CFG_MARK_COUNT")
          == str(sum(r["seg"] == "mark" for r in rows)))

    # -- the two builders agreed -------------------------------------------
    if PLAN_SUMMARY.exists():
        plan = json.loads(PLAN_SUMMARY.read_text(encoding="utf-8"))
        check("both_builders_same_variant", plan["variant"] == variant)
        check("both_builders_same_segments", plan["segments"] == len(rows))
        check("both_builders_same_duration",
              abs(plan["duration_min"] - duration_min) <= 0.01,
              f"plan {plan['duration_min']:.2f} min, table {duration_min:.2f} min")
        check("both_builders_same_trial_order",
              [[o["block"], int(o["rep"]), float(o["td_s"]), o["mode"]]
               for o in order_rows]
              == [[b, r, td, m] for b, r, td, m in plan["trial_order"]])
    else:
        check("plan_summary_present", False, "run f2_sequence.py first")

    # -- provisional and conditional ---------------------------------------
    check("provisional_flagged_in_header",
          header_value(text, "CFG_PROVISIONAL") == "1"
          or c["f0b_dependent"]["dwell_probe_delta_um"] is not None)
    probe = c["f0b_dependent"]["dwell_probe_delta_um"]
    threshold = 3.0 * c["f0b_dependent"]["d_min_um"]
    if probe is None:
        notes.append("  variant_justified_by_f0b: F0b dwell probe not processed "
                     f"yet; running '{variant}' on the procedure's defaults")
        results["variant_justified_by_f0b"] = True
    else:
        want = "full" if probe > threshold else "reduced"
        check("variant_justified_by_f0b", variant == want,
              f"probe {probe:g} um vs 3*d_min {threshold:g} um -> {want}")

    bad = [k for k, v in results.items() if not v]
    print(f"F2 pre-run check ({variant}): {len(results) - len(bad)}/{len(results)} "
          f"checks pass")
    print(f"  {summary['segments']} segments, {summary['trials']} trials, "
          f"{duration_min:.2f} min, {summary['dwelling_min']:.2f} min at rest")
    for line in notes:
        print(line)
    if bad:
        print("FAILED: " + ", ".join(bad))
    else:
        print("  procedure section 5: seed and order logged, dwell grid complete, "
              "dwell timing bounded, directions correct, inside the gate")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
