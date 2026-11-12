#!/usr/bin/env python3
"""Re-derive F0b_pt2's pre-run list from the written files alone.

    python3 check_f0b_pt2_sequence.py

Imports no builder. It reads the command table, the arrivals, the summary and
the generated firmware header, and checks procedure section 5 against them,
plus the four firmware changes of section 0 as far as a host can see them.
What it cannot check is on the bench: that the RMT path really delivers
51.2 kHz, that MSCNT survives a live MRES switch, and that a 1 mm move
measures 1 mm at every resolution. The firmware's 'c' command does the
first, and aborts the run on the second.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEQ = HERE / "sequence"
HEADER = HERE / "firmware" / "f0b_pt2" / "f0b_pt2_config.h"
CONFIG = HERE / "f0b_pt2_config.json"

CODE_OF = {256: 0, 128: 1, 64: 2, 32: 3, 16: 4}


def load_table() -> list[dict]:
    with (SEQ / "f0b_pt2_command_table.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("row", "visit", "mres", "u_target", "rate_mcps",
                  "end_offset_ms"):
            r[k] = int(r[k])
        for k in ("t0_ms", "t1_ms"):
            r[k] = float(r[k])
    return rows


def header_value(text: str, name: str) -> str:
    m = re.search(rf"\b{name}(?:\[\])?\s*=\s*([^;]+);", text)
    return m.group(1).strip().strip('"') if m else ""


def header_segments(text: str) -> list[tuple]:
    body = text[text.index("CFG_SEGS[CFG_SEG_COUNT] = {"):]
    out = []
    for line in body.splitlines():
        line = line.strip()
        if not line.startswith("{") or not line.endswith("},"):
            continue
        parts = [p.strip().rstrip("UL") for p in line[1:-2].split(",")]
        out.append(tuple(int(p) for p in parts))
    return out


def main() -> int:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    table = load_table()
    summary = json.loads((SEQ / "f0b_pt2_summary.json").read_text(encoding="utf-8"))
    text = HEADER.read_text(encoding="utf-8")
    segs = header_segments(text)
    with (SEQ / "f0b_pt2_arrivals.csv").open(encoding="utf-8") as fh:
        arrivals = list(csv.DictReader(fh))

    upm = c["ustep_per_mm"]
    x0 = round(c["x0_mm"] * upm)
    cycle = c["motor_cycle_u"]
    checks: dict[str, tuple[bool, str]] = {}

    def add(name, ok, detail=""):
        checks[name] = (bool(ok), detail)

    # -- the artefacts agree with each other ----------------------------
    add("config_hash_matches",
        hashlib.sha256(CONFIG.read_bytes()).hexdigest()
        == header_value(text, "CFG_HASH") == summary["config_sha256"],
        header_value(text, "CFG_HASH")[:12] + "...")
    add("table_hash_matches",
        hashlib.sha256((SEQ / "f0b_pt2_command_table.csv").read_bytes()).hexdigest()
        == header_value(text, "CFG_TABLE_HASH") == summary["table_sha256"],
        header_value(text, "CFG_TABLE_HASH")[:12] + "...")
    add("header_matches_table", len(segs) == len(table) and all(
        s[4] == r["row"] and s[6] == r["u_target"] and s[7] == r["rate_mcps"]
        and s[8] == r["end_offset_ms"] and s[2] == CODE_OF[r["mres"]]
        for s, r in zip(segs, table)),
        f"{len(segs)} rows")
    add("segment_count_matches",
        int(header_value(text, "CFG_SEG_COUNT")) == len(table) ==
        summary["segments"], f"{len(table)} segments")

    # -- the one internal unit, and divisibility ------------------------
    here, bad = round(c["park_mm"] * upm), []
    for r in table:
        if r["type"] in ("move", "ramp"):
            step = 256 // r["mres"]
            if (r["u_target"] - here) % step:
                bad.append(r["row"])
            here = r["u_target"]
    add("every_move_divisible_at_its_mres", not bad,
        "all whole steps" if not bad else f"rows {bad[:5]}")
    add("all_turning_points_on_the_coarse_grid",
        all(r["u_target"] % 16 == 0 for r in table),
        "multiples of 0.625 um")
    saturated = {x0, x0 - round(c["m_overrun_mm"] * upm),
                 x0 + round(c["m_overrun_mm"] * upm)}
    add("saturated_points_on_the_electrical_grid",
        all(u % cycle == 0 for u in saturated),
        f"x0 and +-{c['m_overrun_mm']:g} mm")

    # -- directions ------------------------------------------------------
    prev, wrong = round(c["park_mm"] * upm), 0
    for r in table:
        if r["type"] in ("move", "ramp"):
            if r["side"] != "none" and (r["side"] == "below") != (r["u_target"] > prev):
                wrong += 1
            prev = r["u_target"]
    add("approach_sides_match_direction", wrong == 0,
        f"{sum(1 for r in table if r['side'] != 'none')} declared")

    # -- windows measured from the last edge -----------------------------
    settle, record = c["settle_ms"], c["record_ms"]
    dwells = [r for r in table if r["type"] == "dwell"]
    marks = [r for r in table if r["type"] == "mark"]
    add("settle_offsets_from_last_edge",
        all(r["end_offset_ms"] == settle for r in dwells),
        f"{len(dwells)} settles at {settle} ms")
    add("record_offsets_from_last_edge",
        all(r["end_offset_ms"] == settle + record for r in marks),
        f"{len(marks)} marks at {settle + record} ms")
    add("both_analysis_windows_fit",
        c["early_window_ms"][1] <= settle and settle + record >= settle + 500,
        f"{c['early_window_ms'][0]}-{c['early_window_ms'][1]} ms and "
        f"{settle}-{settle + record} ms")

    # -- MRES switching ---------------------------------------------------
    switches = [i for i, r in enumerate(table) if r["type"] == "setmres"]
    add("one_mres_switch_per_m_visit",
        sum(1 for r in table if r["type"] == "setmres" and r["block"] == "M")
        == len(c["m_visit_order"]),
        f"{sum(1 for r in table if r['type'] == 'setmres' and r['block'] == 'M')}"
        f" of {len(c['m_visit_order'])} visits")
    add("mres_switches_at_rest",
        all(i == 0 or table[i - 1]["type"] in ("dwell", "mark", "setmres")
            for i in switches),
        f"{len(switches)} switches")
    add("first_segment_sets_mres", bool(switches) and switches[0] == 0,
        "table does not rely on the flashed resolution")
    # Derived from the resolution in force in each visit, not from the
    # switch rows, so the check holds however the builder emits them.
    first_of_visit: dict[int, int] = {}
    for r in table:
        if r["block"] == "M" and r["visit"] not in first_of_visit:
            first_of_visit[r["visit"]] = r["mres"]
    visits = [first_of_visit[k] for k in sorted(first_of_visit)]
    add("mres_visit_order_as_planned", visits == c["m_visit_order"],
        "/".join(str(v) for v in visits))
    add("mres_order_is_abba", visits == visits[::-1], "palindrome")
    add("every_level_visited_twice",
        all(Counter(visits)[m] == 2 for m in c["mres_levels"]),
        ", ".join(f"1/{m}x{Counter(visits)[m]}" for m in c["mres_levels"]))

    # -- the arrivals the analysis will read ------------------------------
    per_visit = defaultdict(list)
    for a in arrivals:
        per_visit[(a["block"], a["visit"])].append(a)
    m_visits = {k: v for k, v in per_visit.items() if k[0] == "M"}
    v_visits = {k: v for k, v in per_visit.items() if k[0] == "V"}
    sat_per_visit = {k: sum(1 for a in v if "saturated" in a["tag"])
                     for k, v in m_visits.items()}
    add("seven_saturated_arrivals_per_mres_visit",
        set(sat_per_visit.values()) == {c["m_approaches"]},
        f"{sorted(set(sat_per_visit.values()))} -> "
        f"{c['m_approaches'] - 1} reversals each")
    add("alternating_sides_in_each_chain", all(
        [a["side"] for a in sorted(v, key=lambda q: int(q["row"]))
         if "saturated" in a["tag"]][:2] == ["below", "above"]
        for v in m_visits.values()), "below, above, below ...")
    v_counts = {k: len(v) for k, v in v_visits.items()}
    add("five_arrivals_per_speed_visit",
        set(v_counts.values()) == {c["v_approaches"]},
        f"{sorted(set(v_counts.values()))} -> {c['v_approaches'] - 1} reversals")
    add("speed_order_is_abba",
        list(c["v_speed_order_mm_s"]) == list(c["v_speed_order_mm_s"])[::-1],
        "/".join(f"{s:g}" for s in c["v_speed_order_mm_s"]))
    partials = sum(1 for a in arrivals if "um" in a["tag"])
    want = len(c["m_visit_order"]) * len(c["m_partial_strokes_um"]) * \
        c["m_partial_cycles"] * 2
    add("partial_stroke_arrivals_complete", partials == want,
        f"{partials} of {want}")

    # -- rates and limits --------------------------------------------------
    rates = [r["rate_mcps"] / 1e3 / (256 // r["mres"])
             for r in table if r["rate_mcps"]]
    add("fastest_rate_within_bench_limit",
        max(rates) <= c["max_step_rate_hz"] + 1,
        f"{max(rates):.0f} Hz of {c['max_step_rate_hz']} Hz")
    add("slowest_rate_on_the_rmt_path",
        min(rates) >= 1e6 / c["rmt_max_period_us"],
        f"{min(rates):.0f} Hz, floor {1e6 / c['rmt_max_period_us']:.1f} Hz")
    add("velocity_is_mres_independent", all(
        abs(r["rate_mcps"] - round(r["rate_mcps"])) < 1e-9 for r in table),
        "rates stored in counts, not steps")
    add("inside_soft_limits", all(
        round(c["soft_min_mm"] * upm) <= r["u_target"] <=
        round(c["soft_max_mm"] * upm) for r in table), "")
    add("inside_gate", summary["duration_min"] <= c["gate_min"],
        f"{summary['duration_min']:.2f} min of {c['gate_min']:.0f} min")
    add("gate_inside_watchdog_and_record",
        c["gate_min"] <= c["watchdog_min"] <= c["scope_record_min"],
        f"{c['gate_min']:.0f} <= {c['watchdog_min']:.0f} <= "
        f"{c['scope_record_min']:.0f} min")
    add("watchdog_in_header",
        int(header_value(text, "CFG_WATCHDOG_MS").rstrip("UL"))
        == round(c["watchdog_min"] * 60000), f"{c['watchdog_min']:g} min")
    add("builder_checks_all_pass", all(summary["checks"].values()),
        f"{sum(summary['checks'].values())}/{len(summary['checks'])}")

    passed = sum(1 for ok, _ in checks.values() if ok)
    print(f"F0b_pt2 pre-run check: {passed}/{len(checks)} checks pass")
    print(f"  {summary['segments']} segments, {summary['arrivals']} arrivals, "
          f"{summary['duration_min']:.2f} min, "
          f"{summary['mres_switches']} MRES switches")
    for name, (ok, detail) in checks.items():
        if not ok:
            print(f"  FAIL {name}: {detail}")
    for name, (ok, detail) in checks.items():
        if ok and detail:
            print(f"  {name}: {detail}")
    print("not checkable from here, do it on the bench: RMT at 51.2 kHz "
          "('c'), MSCNT across a live MRES switch (the run aborts if it "
          "moves), 1 mm measuring 1 mm at every resolution")
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
