#!/usr/bin/env python3
"""Build the F2 command table and the firmware header from f2_config.json.

f2_sequence.py plans the session in time, for the overview figure. This
builds the table the board actually replays: absolute microstep targets,
step rates and dwells, in the same five segment types as F0b and F1.

    python3 f2_make_table.py                  # procedure 4.1, the full session
    python3 f2_make_table.py --variant reduced   # procedure 4.3

writes

    sequence/f2_command_table.csv          the table
    sequence/f2_trial_order.csv            the shuffled (T_d, mode) order
    sequence/f2_table_summary.json         durations, counts, checks, hashes
    firmware/f2_dwell/f2_config.h          the generated firmware header

Everything here is derived from f2_config.json by its own route, including a
second, separately written copy of the seeded shuffle. The result is
cross-checked against sequence/f2_plan_summary.json when that file exists:
segment count, duration, per-block times and the trial order must agree.

THE POST-DWELL RAMP IS PROVISIONAL. Its length is 2*s90, and s90 comes from
F0b, which is not processed yet; f2_config.json carries the procedure's
placeholder defaults (s90 = 3 um, d_min = 5 nm) and the header records
provisional=1. F2 is also conditional: the full version is justified only
if F0b's dwell probe (1 s vs 60 s) exceeded 3*d_min. Re-run this script
after F0b and reflash.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f2_config.json"
SEQ = HERE / "sequence"
HEADER = HERE / "firmware" / "f2_dwell" / "f2_config.h"
PLAN = SEQ / "f2_plan_summary.json"

TYPES = {"move": 0, "ramp": 1, "step": 2, "dwell": 3, "mark": 4}
KINDS = ["transit", "setup", "settle", "dwell", "measure", "hold", "ref"]
BLOCKS = ["R1", "B", "A", "L", "R3"]
SUBS = ["ref", "R", "C", "setup"]

# Fixed bench values, identical to F0b and F1: the three sessions are tied to
# the same driver configuration and the same transit profile.
USTEP_PER_MM = 25600
LEAD_MM = 2.0
STEPS_PER_REV = 200
MRES = 256
MOTOR_CYCLE_U = 1024
SCREW_REV_U = 51200
RUN_CURRENT_MA = 400
HOLD_CURRENT_MA = 400          # IRUN = IHOLD: the dwell is held at run current
TPOWERDOWN = 0
IHOLDDELAY = 0
SOFT_MIN_MM = 2.0
SOFT_MAX_MM = 120.695
MAX_STEP_RATE_HZ = 25600
AQUADB_CLOCK_NS = 120
RMT_HZ = 1_000_000
RMT_CHUNK = 128
RMT_QUEUE_DEPTH = 4
RMT_RING = 6
RMT_MAX_PERIOD_US = 30000
SYNC_GPIO = -1
TIME_LIMIT_MIN = 60.0
MASK64 = (1 << 64) - 1


# -- deterministic shuffle ------------------------------------------------
# Second implementation of the same SplitMix64 + Fisher-Yates used by
# f2_sequence.py, written separately on purpose: if the two disagree the
# cross-check on the trial order fails, which is what it is there for.
def _mix(state: int) -> tuple[int, int]:
    state = (state + 0x9E3779B97F4A7C15) & MASK64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return state, (z ^ (z >> 31)) & MASK64


def shuffle_with_seed(items: list, seed: int) -> list:
    out, state = list(items), seed & MASK64
    for i in reversed(range(1, len(out))):
        state, value = _mix(state)
        j = value % (i + 1)
        out[i], out[j] = out[j], out[i]
    return out


def sub_seed(seed: int, block_index: int, rep: int) -> int:
    return (seed + 1000003 * block_index + 101 * rep) & MASK64


def trial_order(c: dict, blocks: list) -> list:
    order = []
    for index, b in enumerate(blocks):
        for rep in range(1, max(d["reps"] for d in b["dwells"]) + 1):
            batch = []
            for d in b["dwells"]:
                if rep <= d["reps"]:
                    batch += [(d["s"], mode) for mode in b["modes"]]
            if b.get("shuffle", True):
                batch = shuffle_with_seed(batch, sub_seed(c["shuffle_seed"],
                                                          index, rep))
            order += [(b["name"], rep, td, mode) for td, mode in batch]
    return order


# -- session --------------------------------------------------------------
@dataclass
class Row:
    row: int
    block: str
    sub: str
    kind: str
    seg: str
    trial: int
    rep: int
    td_s: float
    u_target: int
    rate_hz: float
    dwell_s: float
    tag: str
    t0: float
    t1: float


class Session:
    def __init__(self, c: dict) -> None:
        self.c = c
        self.rows: list[Row] = []
        self.t = 0.0
        self.block = "R1"
        self.sub = "setup"
        self.kind = "setup"
        self.trial = 0
        self.rep = 0
        self.td = 0.0
        self.u = self.u_of_mm(c["positions_mm"]["B"])
        self.last_dir = 0

    # -- helpers ----------------------------------------------------------
    def u_of_mm(self, mm: float) -> int:
        return round(mm * USTEP_PER_MM)

    def u_of_um(self, um: float) -> int:
        return round(um * USTEP_PER_MM / 1000.0)

    def rate_of_um_s(self, um_s: float) -> float:
        return um_s * USTEP_PER_MM / 1000.0

    def profile_time(self, n_usteps: int, rate_hz: float) -> float:
        """Trapezoidal time for n microsteps peaking at rate_hz (as F0b/F1)."""
        n = abs(n_usteps)
        if n == 0 or rate_hz <= 0.0:
            return 0.0
        a = self.c["accel_mm_s2"] * USTEP_PER_MM
        ramp = rate_hz * rate_hz / a
        if n <= ramp:
            return 2.0 * math.sqrt(n / a)
        return 2.0 * rate_hz / a + (n - ramp) / rate_hz

    def transit_rate(self) -> float:
        return self.c["transit_mm_s"] * USTEP_PER_MM

    # -- emitters ---------------------------------------------------------
    def _add(self, seg: str, u_target: int, rate_hz: float, dwell_s: float,
             tag: str, duration: float) -> None:
        self.rows.append(Row(len(self.rows), self.block, self.sub, self.kind,
                             seg, self.trial, self.rep, self.td, u_target,
                             rate_hz, dwell_s, tag, self.t, self.t + duration))
        self.t += duration
        if u_target != self.u:
            self.last_dir = 1 if u_target > self.u else -1
        self.u = u_target

    def move(self, u_target: int, tag: str) -> None:
        self._add("move", u_target, 0.0, 0.0, tag,
                  self.profile_time(u_target - self.u, self.transit_rate()))

    def ramp(self, u_target: int, rate_hz: float, tag: str) -> None:
        self._add("ramp", u_target, rate_hz, 0.0, tag,
                  self.profile_time(u_target - self.u, rate_hz))

    def dwell(self, seconds: float, tag: str) -> None:
        self._add("dwell", self.u, 0.0, seconds, tag, seconds)

    def mark(self, seconds: float, tag: str) -> None:
        self._add("mark", self.u, 0.0, seconds, tag, seconds)

    def wipe(self, u_target: int, tag: str) -> None:
        """The one preparation every trial shares: 3 mm from below at 1 mm/s."""
        off = round(self.c["wipe_mm"] * USTEP_PER_MM)
        far = abs(u_target - off - self.u) > self.u_of_mm(10.0)
        self.kind = "transit" if far else "setup"
        self.move(u_target - off, tag + (" transit to run-up" if far
                                         else " run-up"))
        self.kind = "setup"
        self.move(u_target, tag + " wipe from below")
        assert self.last_dir == 1, tag


# -- blocks ---------------------------------------------------------------
def ramp_usteps(c: dict) -> int:
    """Post-dwell ramp length in whole microsteps: 2*s90 from F0b."""
    raw = c["ramp_factor_of_s90"] * c["f0b_dependent"]["s90_um"]
    n = max(1, round(raw * USTEP_PER_MM / 1000.0))
    assert n * 1000.0 / USTEP_PER_MM >= 5.0 * c["f0b_dependent"]["d_min_um"], raw
    return n


def reference_block(s: Session, name: str) -> None:
    """15 / 105 mm, both directions, nearer first: the drift and gain pair."""
    c = s.c
    s.block, s.sub, s.kind, s.trial, s.rep, s.td = name, "ref", "ref", 0, 0, 0.0
    pre = round(c["wipe_mm"] * USTEP_PER_MM)
    for u in sorted((s.u_of_mm(mm) for mm in c["reference_mm"]),
                    key=lambda v: abs(v - s.u)):
        mm = u / USTEP_PER_MM
        for side in (-1, 1):
            s.move(u + side * pre, f"{mm:.0f} mm approach")
            s.move(u, f"{mm:.0f} mm " + ("from below" if side < 0 else "from above"))
            s.dwell(c["reference_settle_s"],
                    f"{mm:.0f} mm " + ("up" if side < 0 else "down") + " settle")
            s.mark(c["reference_record_s"],
                   f"{mm:.0f} mm " + ("up" if side < 0 else "down"))


def trial(s: Session, xp: int, td: float, mode: str, number: int, rep: int) -> None:
    """wipe from below -> 0.5 s settle -> dwell T_d -> 6 um ramp -> 2 s hold."""
    c = s.c
    s.sub, s.trial, s.rep, s.td = mode, number, rep, td
    label = f"t{number} r{rep} Td={td:g}s {mode}"
    s.wipe(xp, label)
    s.kind = "settle"
    s.dwell(c["pre_dwell_settle_s"], label + " settle (dwell t=0 at its end)")
    s.kind = "dwell"
    s.dwell(td, label + " DWELL")
    s.kind = "measure"
    direction = -1 if mode == "R" else 1
    s.ramp(xp + direction * ramp_usteps(c), s.rate_of_um_s(c["ramp_rate_um_s"]),
           label + (" post-dwell ramp down (reversal)" if direction < 0
                    else " post-dwell ramp up (continuation)"))
    s.kind = "hold"
    s.dwell(c["post_ramp_hold_s"], label + " hold")


def blocks_of(c: dict, variant: str) -> tuple[list, bool]:
    if variant == "reduced":
        return c["reduced_blocks"], c["reduced_includes_references"]
    return c["blocks"], True


def build(c: dict, blocks: list, references: bool) -> Session:
    s = Session(c)
    order = trial_order(c, blocks)
    position = {b["name"]: c["positions_mm"][b["position"]] for b in blocks}
    if references:
        reference_block(s, "R1")
    for number, (name, rep, td, mode) in enumerate(order, start=1):
        s.block = name
        trial(s, s.u_of_mm(position[name]), td, mode, number, rep)
    if references:
        reference_block(s, "R3")
    return s


# -- checks ---------------------------------------------------------------
def recheck_directions(rows: list[Row]) -> tuple[int, int]:
    """Recompute, from the finished table alone, that every trial is wiped
    from below, that its dwell is motionless, and that the post-dwell ramp
    reverses the wipe in mode R and continues it in mode C."""
    last, ramps, dwells = 0, 0, 0
    prev_u = rows[0].u_target
    for r in rows:
        d = (r.u_target > prev_u) - (r.u_target < prev_u)
        if r.kind == "dwell":
            assert d == 0, f"row {r.row}: dwell moves"
            assert last == 1, f"row {r.row}: dwell not preceded by a wipe up"
            dwells += 1
        if r.kind == "measure":
            want = -last if r.sub == "R" else last
            assert d == want, (f"row {r.row} trial {r.trial}: mode {r.sub} ramp "
                               f"{d:+d} after last motion {last:+d}")
            ramps += 1
        if d:
            last = d
        prev_u = r.u_target
    return ramps, dwells


def coverage(order: list) -> dict:
    out: dict[str, int] = {}
    for name, _rep, td, mode in order:
        key = f"{name}.{td:g}s.{mode}"
        out[key] = out.get(key, 0) + 1
    return out


def planned_coverage(blocks: list) -> dict:
    return {f"{b['name']}.{d['s']:g}s.{m}": d["reps"]
            for b in blocks for d in b["dwells"] for m in b["modes"]}


def pearson(xs: list, ys: list) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def dwell_time_correlation(rows: list[Row]) -> dict:
    """log10(T_d) against when the trial starts, per block. Shuffling is only
    worth doing if this stays small; a residual trend in time would otherwise
    be indistinguishable from an ageing effect."""
    first: dict[int, Row] = {}
    for r in rows:
        if r.trial and r.trial not in first:
            first[r.trial] = r
    out = {}
    for block in dict.fromkeys(r.block for r in first.values()):
        rows_b = [r for r in first.values() if r.block == block]
        if len(rows_b) >= 3:
            out[block] = pearson([math.log10(r.td_s) for r in rows_b],
                                 [r.t0 for r in rows_b])
    return out


def cross_check_plan(s: Session, order: list) -> dict:
    """Compare with f2_sequence.py's independent plan, if it has been run."""
    if not PLAN.exists():
        return {"plan_summary_present": False}
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    by_block: dict[str, float] = {}
    for r in s.rows:
        by_block[r.block] = by_block.get(r.block, 0.0) + (r.t1 - r.t0)
    same_blocks = all(
        abs(by_block.get(b, 0.0) / 60.0 - m) <= 0.01
        for b, m in plan["blocks_min"].items())
    return {
        "plan_summary_present": True,
        "plan_variant_matches": plan["variant"] == s.c["_resolved_variant"],
        "plan_segments_match": plan["segments"] == len(s.rows),
        "plan_trials_match": plan["trials"] == len({r.trial for r in s.rows if r.trial}),
        "plan_duration_matches_0p01min": abs(plan["duration_min"] - s.t / 60.0) <= 0.01,
        "plan_blocks_match_0p01min": same_blocks,
        "plan_ramp_usteps_match": plan["ramp_usteps"] == ramp_usteps(s.c),
        "plan_trial_order_matches": [list(x) for x in order] == [
            [b, r, td, m] for b, r, td, m in plan["trial_order"]],
    }


def checks(s: Session, blocks: list, order: list) -> dict:
    c = s.c
    u_min, u_max = s.u_of_mm(SOFT_MIN_MM), s.u_of_mm(SOFT_MAX_MM)
    rates = [r.rate_hz for r in s.rows if r.seg == "ramp"]
    positions = [s.u_of_mm(c["positions_mm"][b["position"]]) for b in blocks]
    lo_mm, hi_mm = c["constrained_span_mm"]
    dwells = [r for r in s.rows if r.kind == "dwell"]
    longest_s = max(r.dwell_s for r in s.rows if r.seg in ("dwell", "mark"))
    corr = dwell_time_correlation(s.rows)
    out = {
        "positions_on_electrical_grid": all(u % MOTOR_CYCLE_U == 0 for u in positions),
        "positions_inside_constrained_span": all(
            lo_mm <= c["positions_mm"][b["position"]] <= hi_mm for b in blocks),
        "inside_soft_limits": all(u_min <= r.u_target <= u_max for r in s.rows),
        "under_driver_limit": all(r <= MAX_STEP_RATE_HZ + 1e-9 for r in rates),
        "ramp_rate_on_timed_path": all(
            r.rate_hz < 1e6 / RMT_MAX_PERIOD_US + 1e-9
            for r in s.rows if r.kind == "measure"),
        "dwells_are_whole_ms": all(
            abs(r.dwell_s * 1000.0 - round(r.dwell_s * 1000.0)) < 1e-6
            for r in s.rows if r.seg in ("dwell", "mark")),
        "dwells_motionless": all(r.u_target == s.rows[r.row - 1].u_target
                                 for r in dwells),
        "every_trial_has_one_dwell": len(dwells) == len(order),
        "coverage_matches_grid": coverage(order) == planned_coverage(blocks),
        "baseline_dwell_present": any(td == c["baseline_dwell_s"]
                                      for _b, _r, td, _m in order),
        "dwell_time_decorrelated": all(abs(v) <= c["max_dwell_time_correlation"]
                                       for v in corr.values()),
        "irun_equals_ihold": RUN_CURRENT_MA == HOLD_CURRENT_MA and TPOWERDOWN == 0,
        "inside_gate": s.t / 60.0 <= c["gate_min"],
        "gate_plus_longest_dwell_inside_record": (
            c["gate_min"] + longest_s / 60.0 <= c["scope_record_min"]),
        "inside_scope_record": s.t / 60.0 <= c["scope_record_min"],
    }
    out.update(cross_check_plan(s, order))
    return out


# -- summary and writers --------------------------------------------------
def summarise(s: Session, variant: str, order: list, ok: dict) -> dict:
    by_block, by_kind, by_seg = {}, {}, {}
    for r in s.rows:
        by_block[r.block] = by_block.get(r.block, 0.0) + (r.t1 - r.t0)
        by_kind[r.kind] = by_kind.get(r.kind, 0.0) + (r.t1 - r.t0)
        by_seg[r.seg] = by_seg.get(r.seg, 0.0) + (r.t1 - r.t0)
    dwelling = sum(r.dwell_s for r in s.rows if r.kind == "dwell")
    overhead = sum(r.t1 - r.t0 for r in s.rows if r.trial and r.kind != "dwell")
    n = ramp_usteps(s.c)
    return {
        "variant": variant,
        "segments": len(s.rows),
        "trials": len(order),
        "duration_min": s.t / 60.0,
        "gate_min": s.c["gate_min"],
        "margin_to_gate_min": s.c["gate_min"] - s.t / 60.0,
        "provisional": True,
        "conditional_on_f0b_dwell_probe": {
            "delta_um": s.c["f0b_dependent"]["dwell_probe_delta_um"],
            "threshold_um": 3.0 * s.c["f0b_dependent"]["d_min_um"],
            "significant": s.c["f0b_dependent"]["dwell_probe_significant"],
        },
        "s90_um": s.c["f0b_dependent"]["s90_um"],
        "d_min_um": s.c["f0b_dependent"]["d_min_um"],
        "ramp_usteps": n,
        "ramp_um": n * 1000.0 / USTEP_PER_MM,
        "ramp_rate_hz": s.rate_of_um_s(s.c["ramp_rate_um_s"]),
        "ramp_seconds": s.profile_time(n, s.rate_of_um_s(s.c["ramp_rate_um_s"])),
        "shuffle_seed": s.c["shuffle_seed"],
        "dwelling_min": dwelling / 60.0,
        "trial_overhead_s": overhead / len(order) if order else 0.0,
        "longest_dwell_s": max(td for _b, _r, td, _m in order),
        "dwell_times_s": sorted({td for _b, _r, td, _m in order}),
        "blocks_min": {b: by_block[b] / 60.0 for b in BLOCKS if b in by_block},
        "kinds_min": {k: by_kind[k] / 60.0 for k in KINDS if k in by_kind},
        "segment_types_min": {k: v / 60.0 for k, v in sorted(by_seg.items())},
        "segment_type_counts": {k: sum(r.seg == k for r in s.rows)
                                for k in sorted(by_seg)},
        "coverage": coverage(order),
        "dwell_time_correlation": dwell_time_correlation(s.rows),
        "checks": ok,
    }


def write_table(path: Path, rows: list[Row]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["row", "block", "sub", "kind", "seg", "trial", "rep",
                    "td_s", "u_target", "rate_hz", "dwell_s", "tag",
                    "t0_s", "t1_s"])
        for r in rows:
            w.writerow([r.row, r.block, r.sub, r.kind, r.seg, r.trial, r.rep,
                        f"{r.td_s:g}", r.u_target, f"{r.rate_hz:.6f}",
                        f"{r.dwell_s:.3f}", r.tag, f"{r.t0:.3f}", f"{r.t1:.3f}"])
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_trial_order(path: Path, s: Session) -> str:
    """The shuffled order, one row per trial: the pre-run check in section 5
    is read off this file, and the analysis joins to it by trial number."""
    first: dict[int, dict] = {}
    for r in s.rows:
        if not r.trial:
            continue
        slot = first.setdefault(r.trial, {"trial": r.trial, "block": r.block,
                                          "sub": r.sub, "rep": r.rep,
                                          "td_s": r.td_s, "t_start_s": r.t0})
        if r.kind == "dwell":
            slot["t_dwell_start_s"] = r.t0
            slot["t_dwell_end_s"] = r.t1
        if r.kind == "measure":
            slot["t_ramp_start_s"] = r.t0
            slot["u_ramp_end"] = r.u_target
            slot["ramp_dir"] = "down" if r.u_target < s.rows[r.row - 1].u_target else "up"
    fields = ["trial", "block", "rep", "order_in_rep", "td_s", "mode", "ramp_dir",
              "t_start_s", "t_dwell_start_s", "t_dwell_end_s", "t_ramp_start_s",
              "u_ramp_end"]
    seen: dict[tuple, int] = {}
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for trial_no in sorted(first):
            v = first[trial_no]
            key = (v["block"], v["rep"])
            seen[key] = seen.get(key, 0) + 1
            w.writerow({
                "trial": v["trial"], "block": v["block"], "rep": v["rep"],
                "order_in_rep": seen[key], "td_s": f"{v['td_s']:g}",
                "mode": v["sub"], "ramp_dir": v["ramp_dir"],
                "t_start_s": f"{v['t_start_s']:.3f}",
                "t_dwell_start_s": f"{v['t_dwell_start_s']:.3f}",
                "t_dwell_end_s": f"{v['t_dwell_end_s']:.3f}",
                "t_ramp_start_s": f"{v['t_ramp_start_s']:.3f}",
                "u_ramp_end": v["u_ramp_end"]})
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_header(path: Path, c: dict, s: Session, variant: str, order: list,
                 cfg_hash: str, table_hash: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    marks = sum(r.seg == "mark" for r in s.rows)
    longest_ms = round(max(r.dwell_s for r in s.rows if r.seg == "dwell") * 1000)
    out = [
        "/* GENERATED by f2_make_table.py. Do not edit by hand.",
        " *",
        " * PROVISIONAL: the post-dwell ramp is 2*s90 long, built from the",
        f" * placeholder s90 = {c['f0b_dependent']['s90_um']} um and "
        f"d_min = {c['f0b_dependent']['d_min_um']} um in f2_config.json.",
        " * F2 is also CONDITIONAL: the full grid is justified only if F0b's",
        " * dwell probe (1 s vs 60 s) exceeded 3*d_min; otherwise run the",
        " * reduced variant. Re-run f2_make_table.py after F0b and reflash. */",
        "#pragma once",
        "#include <stdint.h>",
        "",
        f'constexpr char CFG_HASH[] = "{cfg_hash}";',
        f'constexpr char CFG_TABLE_HASH[] = "{table_hash}";',
        'constexpr char CFG_EXPERIMENT[] = "F2_dwell";',
        f'constexpr char CFG_VARIANT[] = "{variant}";',
        "constexpr uint8_t CFG_PROVISIONAL = 1;",
        f"constexpr double CFG_S90_UM = {c['f0b_dependent']['s90_um']:.6f};",
        f"constexpr double CFG_D_MIN_UM = {c['f0b_dependent']['d_min_um']:.6f};",
        f"constexpr double CFG_LEAD_MM = {LEAD_MM:.6f};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {STEPS_PER_REV};",
        f"constexpr uint16_t CFG_MRES = {MRES};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {RUN_CURRENT_MA};",
        f"constexpr uint16_t CFG_HOLD_CURRENT_MA = {HOLD_CURRENT_MA};",
        f"constexpr uint8_t CFG_TPOWERDOWN = {TPOWERDOWN};",
        f"constexpr uint8_t CFG_IHOLDDELAY = {IHOLDDELAY};",
        f"constexpr double CFG_V_MM_S = {c['transit_mm_s']:.6f};",
        f"constexpr double CFG_A_MM_S2 = {c['accel_mm_s2']:.6f};",
        f"constexpr uint32_t CFG_MAX_STEP_RATE_HZ = {MAX_STEP_RATE_HZ}UL;",
        f"constexpr int32_t CFG_PARK_U = {s.u_of_mm(c['positions_mm']['B'])};",
        f"constexpr int32_t CFG_X_A_U = {s.u_of_mm(c['positions_mm']['A'])};",
        f"constexpr int32_t CFG_X_B_U = {s.u_of_mm(c['positions_mm']['B'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {s.u_of_mm(SOFT_MIN_MM)};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {s.u_of_mm(SOFT_MAX_MM)};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {MOTOR_CYCLE_U};",
        f"constexpr uint32_t CFG_SCREW_REV_U = {SCREW_REV_U}UL;",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(TIME_LIMIT_MIN * 60000)}UL;",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['gate_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {AQUADB_CLOCK_NS};",
        f"constexpr uint32_t CFG_RMT_HZ = {RMT_HZ}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {RMT_CHUNK};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {RMT_QUEUE_DEPTH};",
        f"constexpr uint8_t CFG_RMT_RING = {RMT_RING};",
        f"constexpr uint32_t CFG_RMT_MAX_PERIOD_US = {RMT_MAX_PERIOD_US}UL;",
        f"constexpr int8_t CFG_SYNC_GPIO = {SYNC_GPIO};",
        f"constexpr uint32_t CFG_PLAN_DURATION_MS = {round(s.t * 1000)}UL;",
        "",
        "/* F2-specific. The dwell is the measurement, so its timing is",
        "   checked against the table at run time and reported per segment. */",
        f"constexpr uint64_t CFG_SHUFFLE_SEED = {c['shuffle_seed']}ULL;",
        f"constexpr uint16_t CFG_TRIAL_COUNT = {len(order)};",
        f"constexpr uint32_t CFG_LONGEST_DWELL_MS = {longest_ms}UL;",
        f"constexpr uint32_t CFG_DWELL_TOL_MS = "
        f"{round(c['dwell_timing_tolerance_ms'])}UL;",
        f"constexpr uint32_t CFG_RAMP_USTEPS = {ramp_usteps(c)}UL;",
        f"constexpr double CFG_RAMP_RATE_HZ = "
        f"{s.rate_of_um_s(c['ramp_rate_um_s']):.6f};",
        f"constexpr uint32_t CFG_PRE_DWELL_SETTLE_MS = "
        f"{round(c['pre_dwell_settle_s'] * 1000)}UL;",
        f"constexpr uint32_t CFG_POST_RAMP_HOLD_MS = "
        f"{round(c['post_ramp_hold_s'] * 1000)}UL;",
        "",
        "enum F2SegType : uint8_t { F2_MOVE = 0, F2_RAMP = 1, F2_STEP = 2,",
        "                           F2_DWELL = 3, F2_MARK = 4 };",
        "",
        "/* rate_mhz is the step rate in milli-hertz: 12800 means 12.8 steps/s.",
        "   dwell_ms is the duration of F2_DWELL and F2_MARK. kind separates",
        "   the trial's timed dwell from the settle, the hold and the",
        "   reference stops, which are the same segment type. */",
        "struct F2Seg {",
        "  uint8_t type;",
        "  uint8_t kind;",
        "  uint8_t block;",
        "  uint8_t sub;",
        "  uint16_t row;",
        "  uint16_t trial;",
        "  int32_t u_target;",
        "  uint32_t rate_mhz;",
        "  uint32_t dwell_ms;",
        "};",
        "",
        f"constexpr uint16_t CFG_SEG_COUNT = {len(s.rows)};",
        f"constexpr uint16_t CFG_MARK_COUNT = {marks};",
        f"constexpr uint16_t CFG_DWELL_WINDOWS = "
        f"{sum(r.kind == 'dwell' for r in s.rows)};",
        "constexpr const char *CFG_BLOCK_NAMES[] = {"
        + ", ".join(f'"{b}"' for b in BLOCKS) + "};",
        "constexpr const char *CFG_SUB_NAMES[] = {"
        + ", ".join(f'"{x}"' for x in SUBS) + "};",
        "constexpr const char *CFG_KIND_NAMES[] = {"
        + ", ".join(f'"{k}"' for k in KINDS) + "};",
        "",
        "constexpr F2Seg CFG_SEGS[CFG_SEG_COUNT] = {",
    ]
    for r in s.rows:
        out.append("  {%d, %d, %d, %d, %d, %d, %d, %dUL, %dUL}," % (
            TYPES[r.seg], KINDS.index(r.kind), BLOCKS.index(r.block),
            SUBS.index(r.sub), r.row, r.trial, r.u_target,
            round(r.rate_hz * 1000.0), round(r.dwell_s * 1000.0)))
    out[-1] = out[-1].rstrip(",")
    out += ["};", ""]
    path.write_text("\n".join(out), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variant", choices=("full", "reduced"), default=None,
                    help="procedure 4.1 (full) or 4.3 (reduced); "
                         "default: the variant in f2_config.json")
    args = ap.parse_args()

    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    cfg_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    variant = args.variant or c["variant"]
    c["_resolved_variant"] = variant
    blocks, references = blocks_of(c, variant)
    order = trial_order(c, blocks)

    s = build(c, blocks, references)
    ramps, dwells = recheck_directions(s.rows)
    ok = checks(s, blocks, order)
    summary = summarise(s, variant, order, ok)
    summary["ramps_checked"] = ramps
    summary["dwells_checked"] = dwells
    summary["references_included"] = references

    SEQ.mkdir(exist_ok=True)
    table_hash = write_table(SEQ / "f2_command_table.csv", s.rows)
    order_hash = write_trial_order(SEQ / "f2_trial_order.csv", s)
    summary["config_sha256"] = cfg_hash
    summary["table_sha256"] = table_hash
    summary["trial_order_sha256"] = order_hash
    (SEQ / "f2_table_summary.json").write_text(json.dumps(summary, indent=1) + "\n",
                                               encoding="utf-8")
    write_header(HEADER, c, s, variant, order, cfg_hash, table_hash)

    probe = c["f0b_dependent"]["dwell_probe_delta_um"]
    print(f"F2 table ({variant}): {summary['segments']} segments, "
          f"{summary['trials']} trials, {summary['duration_min']:.2f} min "
          f"(gate {c['gate_min']:.0f} min, "
          f"{summary['margin_to_gate_min']:.1f} min of margin)")
    for b, m in summary["blocks_min"].items():
        print(f"  {b:3s} {m:6.2f} min")
    print(f"  dwelling {summary['dwelling_min']:.2f} min over "
          f"{summary['trials']} trials, {summary['trial_overhead_s']:.1f} s of "
          f"overhead each, longest dwell {summary['longest_dwell_s']:g} s")
    print(f"  post-dwell ramp {summary['ramp_um']:.4f} um "
          f"({summary['ramp_usteps']} usteps) at {summary['ramp_rate_hz']:.2f} Hz "
          f"= {summary['ramp_seconds']:.2f} s, on the timed pulse path")
    print(f"  seed {c['shuffle_seed']}; log10(T_d) vs session time "
          + ", ".join(f"{b} r={v:+.3f}"
                      for b, v in summary["dwell_time_correlation"].items()))
    print(f"  {ramps} post-dwell ramps and {dwells} dwells re-checked from the "
          f"table alone")
    print(f"  PROVISIONAL from s90={summary['s90_um']:g} um; F0b dwell probe "
          + ("not processed yet" if probe is None else f"delta={probe:g} um")
          + f", threshold 3*d_min="
            f"{summary['conditional_on_f0b_dwell_probe']['threshold_um']:g} um")
    bad = [k for k, v in ok.items() if not v]
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    print(f"wrote {HEADER.relative_to(HERE)}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
