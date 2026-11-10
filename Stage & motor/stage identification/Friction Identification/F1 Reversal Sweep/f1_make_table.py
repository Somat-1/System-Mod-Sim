#!/usr/bin/env python3
"""Build the F1 command table and the firmware header from f1_config.json.

f1_sequence.py plans the session in time, for the overview figure. This
builds the table the board actually replays: absolute microstep targets,
step rates and dwells, in the same five segment types as F0b.

    python3 f1_make_table.py

writes

    sequence/f1_command_table.csv              the table
    sequence/f1_table_summary.json             durations, counts and the checks
    firmware/f1_reversal_sweep/f1_config.h     the generated firmware header

THE STROKE GRID IS PROVISIONAL. P, M and Q are built from s90 and d_min,
which F0b measures. Until F0b is processed, f1_config.json carries the
procedure's placeholder defaults (s90 = 3 um, d_min = 5 nm) and the header
records provisional=1. Re-run this script after F0b and reflash.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f1_config.json"
SEQ = HERE / "sequence"
HEADER = HERE / "firmware" / "f1_reversal_sweep" / "f1_config.h"

TYPES = {"move": 0, "ramp": 1, "step": 2, "dwell": 3, "mark": 4}
GROUPS = ["R1", "A", "B", "R2", "C", "R3"]
SUBS = ["R", "V", "P", "M", "Q", "setup"]

# Fixed bench values shared with F0b; the procedure ties the two sessions to
# the same driver configuration and the same transit profile.
USTEP_PER_MM = 25600
LEAD_MM = 2.0
STEPS_PER_REV = 200
MRES = 256
MOTOR_CYCLE_U = 1024
SCREW_REV_U = 51200
RUN_CURRENT_MA = 400
HOLD_CURRENT_MA = 400
TPOWERDOWN = 0
IHOLDDELAY = 0
SOFT_MIN_MM = 2.0
SOFT_MAX_MM = 120.695
SETTLE_S = 2.0
RECORD_S = 0.5
MAX_STEP_RATE_HZ = 25600
AQUADB_CLOCK_NS = 120
RMT_HZ = 1_000_000
RMT_CHUNK = 128
RMT_QUEUE_DEPTH = 4
RMT_RING = 6
RMT_MAX_PERIOD_US = 30000
SYNC_GPIO = -1
TIME_LIMIT_MIN = 60.0


@dataclass
class Row:
    row: int
    group: str
    sub: str
    seg: str
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
        self.group = "R1"
        self.sub = "setup"
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
        """Trapezoidal time for n microsteps peaking at rate_hz (see F0b)."""
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
        self.rows.append(Row(len(self.rows), self.group, self.sub, seg,
                             u_target, rate_hz, dwell_s, tag,
                             self.t, self.t + duration))
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

    def steps(self, u_target: int, dwell_s: float, tag: str) -> None:
        n = abs(u_target - self.u)
        self._add("step", u_target, 0.0, dwell_s, tag, n * dwell_s)

    def dwell(self, seconds: float, tag: str) -> None:
        self._add("dwell", self.u, 0.0, seconds, tag, seconds)

    def mark(self, tag: str) -> None:
        self.dwell(SETTLE_S, tag + " settle")
        self._add("mark", self.u, 0.0, RECORD_S, tag, RECORD_S)

    def wipe(self, u_target: int, from_below: bool, tag: str) -> None:
        off = round(self.c["wipe_mm"] * USTEP_PER_MM)
        self.move(u_target - off if from_below else u_target + off,
                  tag + " run-up")
        self.move(u_target, tag + " wipe " +
                  ("from below" if from_below else "from above"))
        assert self.last_dir == (1 if from_below else -1), tag
        self.dwell(self.c["settle_s"], tag + " settle")

    def half(self, u_target: int, rate_hz: float, tag: str,
             step_mode: bool) -> None:
        """One half-cycle followed by its turning pause."""
        if step_mode:
            self.steps(u_target, self.c["step_dwell_s"], tag + " steps")
        else:
            self.ramp(u_target, rate_hz, tag)
        self.dwell(self.c["turn_pause_s"], tag + " turn")


# -- stroke grid ----------------------------------------------------------
def strokes_usteps(c: dict) -> list[int]:
    """P stroke grid in whole microsteps, symmetric about the test position."""
    ustep_um = 1000.0 / USTEP_PER_MM
    s90 = c["f0b_dependent"]["s90_um"]
    d_min = c["f0b_dependent"]["d_min_um"]
    raw = [s90 * f for f in c["p_stroke_factors"]] + c["p_anchor_strokes_um"]
    out = []
    for s in raw:
        half = max(1, round(s / 2.0 / ustep_um))
        stroke = 2 * half
        assert stroke * ustep_um >= max(ustep_um, 5.0 * d_min), s
        out.append(stroke)
    return out


def centred_loops(s: Session, xp: int, stroke: int, cycles: int,
                  half_s: float, tag: str) -> None:
    c = s.c
    rate = min(s.rate_of_um_s(stroke * 1000.0 / USTEP_PER_MM / half_s),
               s.rate_of_um_s(c["p_rate_cap_um_s"]))
    step_mode = stroke <= c["step_segment_max_usteps"]
    lo, hi = xp - stroke // 2, xp + stroke // 2
    s.half(lo, rate, f"{tag} to start", step_mode)
    for k in range(cycles):
        s.half(hi, rate, f"{tag} c{k + 1} up", step_mode)
        s.half(lo, rate, f"{tag} c{k + 1} down", step_mode)
    s.half(xp, rate, f"{tag} return", step_mode)


def reference_block(s: Session, name: str) -> None:
    s.group, s.sub = name, "R"
    c = s.c
    pre = round(c["wipe_mm"] * USTEP_PER_MM)
    targets = sorted((s.u_of_mm(mm) for mm in c["reference_mm"]),
                     key=lambda u: abs(u - s.u))
    for u in targets:
        for side in (-1, 1):
            mm = u / USTEP_PER_MM
            s.move(u + side * pre, f"{mm:.0f} mm approach")
            s.move(u, f"{mm:.0f} mm " +
                   ("from below" if side < 0 else "from above"))
            s.mark(f"{mm:.0f} mm " + ("up" if side < 0 else "down"))


def position_block(s: Session, group: str, xp_mm: float,
                   decreasing: bool = False, rate_check: bool = False) -> None:
    c = s.c
    s.group = group
    xp = s.u_of_mm(xp_mm)
    s90 = c["f0b_dependent"]["s90_um"]

    # V: virgin curves. Reverse off a wipe, slowly for the first 2*s90, then
    # faster out to 100 um.
    s.sub = "V"
    slow = s.u_of_um(c["v_slow_factor_of_s90"] * s90)
    far = s.u_of_um(c["v_ramp_um"])
    for side in (-1, 1):
        s.wipe(xp, side < 0, f"V {'up' if side < 0 else 'down'} prep")
        d = side
        s.ramp(xp + d * slow, s.rate_of_um_s(c["v_slow_um_s"]),
               "V virgin slow (reversal)")
        s.ramp(xp + d * far, s.rate_of_um_s(c["v_fast_um_s"]),
               "V virgin fast (continuation)")

    # P: centred loops, increasing, and decreasing at the middle position.
    s.sub = "P"
    s.wipe(xp, True, "P prep")
    grid = strokes_usteps(c)
    order = grid + (grid[::-1] if decreasing else [])
    for i, stroke in enumerate(order):
        label = "inc" if i < len(grid) else "dec"
        centred_loops(s, xp, stroke, c["p_cycles"], c["p_half_cycle_s"],
                      f"P s={stroke}u {label}")

    # M: nested loops. Even microstep counts keep +/- S_o/2 symmetric.
    s.sub = "M"
    so = 2 * max(1, round(c["m_outer_factor"] * s90 * USTEP_PER_MM / 2000.0))
    si = max(1, s.u_of_um(c["m_inner_factor"] * s90))
    ro = s.rate_of_um_s(so * 1000.0 / USTEP_PER_MM / c["p_half_cycle_s"])
    ri = s.rate_of_um_s(si * 1000.0 / USTEP_PER_MM / c["p_half_cycle_s"])
    outer_steps = so <= c["step_segment_max_usteps"]
    inner_steps = si <= c["step_segment_max_usteps"]
    s.wipe(xp - so // 2, True, "M prep")
    for k in range(c["m_repeats"]):
        s.half(xp + so // 2, ro, f"M{k + 1} outer up", outer_steps)
        s.half(xp, ro, f"M{k + 1} outer down to xp", outer_steps)
        s.half(xp + si, ri, f"M{k + 1} inner up", inner_steps)
        s.half(xp, ri, f"M{k + 1} inner down", inner_steps)
        s.half(xp - so // 2, ro, f"M{k + 1} outer down rest", outer_steps)

    # Q: rate check at s90, middle position only.
    if rate_check:
        s.sub = "Q"
        q = strokes_usteps({**c, "p_stroke_factors": [1.0],
                            "p_anchor_strokes_um": []})[0]
        s.wipe(xp, True, "Q prep")
        for half_s in c["q_half_cycle_s"]:
            centred_loops(s, xp, q, c["q_cycles"], half_s, f"Q {half_s:g}s")


def build(c: dict, include_r2: bool) -> Session:
    s = Session(c)
    reference_block(s, "R1")
    position_block(s, "A", c["positions_mm"]["A"])
    position_block(s, "B", c["positions_mm"]["B"], decreasing=True,
                   rate_check=True)
    if include_r2:
        reference_block(s, "R2")
    position_block(s, "C", c["positions_mm"]["C"])
    reference_block(s, "R3")
    return s


def plan(c: dict) -> tuple[Session, bool]:
    """R2 is dropped if the session would exceed the gate (procedure section 4)."""
    s = build(c, True)
    if s.t / 60.0 > c["gate_min"]:
        return build(c, False), False
    return s, True


# -- checks ---------------------------------------------------------------
def checks(s: Session, r2: bool) -> dict[str, bool]:
    c = s.c
    u_min, u_max = s.u_of_mm(SOFT_MIN_MM), s.u_of_mm(SOFT_MAX_MM)
    grid = strokes_usteps(c)
    rates = [r.rate_hz for r in s.rows if r.seg == "ramp"]
    cap = s.rate_of_um_s(c["p_rate_cap_um_s"])
    loops = [r for r in s.rows if r.sub in ("P", "Q", "M") and r.seg == "ramp"]
    positions = [s.u_of_mm(v) for v in c["positions_mm"].values()]
    lo, hi = c["constrained_span_mm"]
    return {
        "positions_on_electrical_grid": all(u % MOTOR_CYCLE_U == 0
                                            for u in positions),
        "positions_inside_constrained_span": all(
            lo <= v <= hi for v in c["positions_mm"].values()),
        "inside_soft_limits": all(u_min <= r.u_target <= u_max
                                  for r in s.rows),
        "under_driver_limit": all(r <= MAX_STEP_RATE_HZ + 1e-9 for r in rates),
        "loops_under_rate_cap": all(r.rate_hz <= cap + 1e-9 for r in loops),
        "strokes_even_usteps": all(v % 2 == 0 for v in grid),
        "strokes_monotonic": all(b > a for a, b in zip(grid, grid[1:])),
        "strokes_above_5_dmin": all(
            v * 1000.0 / USTEP_PER_MM >= 5.0 * c["f0b_dependent"]["d_min_um"]
            for v in grid),
        "inside_gate": s.t / 60.0 <= c["gate_min"],
        "inside_scope_record": s.t / 60.0 <= c["scope_record_min"],
    }


def summarise(s: Session, r2: bool, ok: dict) -> dict:
    by_group, by_sub, by_seg = {}, {}, {}
    for r in s.rows:
        by_group[r.group] = by_group.get(r.group, 0.0) + (r.t1 - r.t0)
        by_sub[f"{r.group}.{r.sub}"] = by_sub.get(f"{r.group}.{r.sub}", 0.0) \
            + (r.t1 - r.t0)
        by_seg[r.seg] = by_seg.get(r.seg, 0.0) + (r.t1 - r.t0)
    return {
        "segments": len(s.rows),
        "duration_min": s.t / 60.0,
        "r2_included": r2,
        "gate_min": s.c["gate_min"],
        "margin_to_gate_min": s.c["gate_min"] - s.t / 60.0,
        "provisional": True,
        "s90_um": s.c["f0b_dependent"]["s90_um"],
        "d_min_um": s.c["f0b_dependent"]["d_min_um"],
        "strokes_usteps": strokes_usteps(s.c),
        "groups_min": {g: by_group[g] / 60.0 for g in GROUPS if g in by_group},
        "sub_blocks_min": {k: v / 60.0 for k, v in by_sub.items()},
        "segment_types_min": {k: v / 60.0 for k, v in sorted(by_seg.items())},
        "segment_type_counts": {k: sum(r.seg == k for r in s.rows)
                                for k in sorted(by_seg)},
        "slowest_ramp_hz": min((r.rate_hz for r in s.rows if r.seg == "ramp"),
                               default=0.0),
        "fastest_ramp_hz": max((r.rate_hz for r in s.rows if r.seg == "ramp"),
                               default=0.0),
        "checks": ok,
    }


# -- writers --------------------------------------------------------------
def write_table(path: Path, rows: list[Row]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["row", "group", "sub", "seg", "u_target", "rate_hz",
                    "dwell_s", "tag", "t0_s", "t1_s"])
        for r in rows:
            w.writerow([r.row, r.group, r.sub, r.seg, r.u_target,
                        f"{r.rate_hz:.6f}", f"{r.dwell_s:.3f}", r.tag,
                        f"{r.t0:.3f}", f"{r.t1:.3f}"])
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_header(path: Path, c: dict, s: Session, r2: bool, cfg_hash: str,
                 table_hash: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    marks = sum(r.seg == "mark" for r in s.rows)
    park = s.u_of_mm(c["positions_mm"]["B"])
    out = [
        "/* GENERATED by f1_make_table.py. Do not edit by hand.",
        " *",
        " * PROVISIONAL: the P, M and Q strokes are built from the placeholder",
        f" * s90 = {c['f0b_dependent']['s90_um']} um and "
        f"d_min = {c['f0b_dependent']['d_min_um']} um in f1_config.json.",
        " * Re-run f1_make_table.py after F0b is processed and reflash. */",
        "#pragma once",
        "#include <stdint.h>",
        "",
        f'constexpr char CFG_HASH[] = "{cfg_hash}";',
        f'constexpr char CFG_TABLE_HASH[] = "{table_hash}";',
        f'constexpr char CFG_EXPERIMENT[] = "F1_reversal_sweep";',
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
        f"constexpr uint32_t CFG_SETTLE_MS = {round(SETTLE_S * 1000)}UL;",
        f"constexpr uint32_t CFG_RECORD_MS = {round(RECORD_S * 1000)}UL;",
        f"constexpr int32_t CFG_PARK_U = {park};",
        f"constexpr int32_t CFG_X_A_U = {s.u_of_mm(c['positions_mm']['A'])};",
        f"constexpr int32_t CFG_X_B_U = {s.u_of_mm(c['positions_mm']['B'])};",
        f"constexpr int32_t CFG_X_C_U = {s.u_of_mm(c['positions_mm']['C'])};",
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
        f"constexpr uint8_t CFG_R2_INCLUDED = {1 if r2 else 0};",
        f"constexpr uint32_t CFG_PLAN_DURATION_MS = {round(s.t * 1000)}UL;",
        "",
        "enum F1SegType : uint8_t { F1_MOVE = 0, F1_RAMP = 1, F1_STEP = 2,",
        "                           F1_DWELL = 3, F1_MARK = 4 };",
        "",
        "/* rate_mhz is the step rate in milli-hertz: 6400 means 6.4 steps/s.",
        "   dwell_ms is the dwell for F1_DWELL and F1_MARK, and the dwell",
        "   after each single microstep for F1_STEP. */",
        "struct F1Seg {",
        "  uint8_t type;",
        "  uint8_t group;",
        "  uint8_t sub;",
        "  uint16_t row;",
        "  int32_t u_target;",
        "  uint32_t rate_mhz;",
        "  uint32_t dwell_ms;",
        "};",
        "",
        f"constexpr uint16_t CFG_SEG_COUNT = {len(s.rows)};",
        f"constexpr uint16_t CFG_MARK_COUNT = {marks};",
        "constexpr const char *CFG_GROUP_NAMES[] = {"
        + ", ".join(f'"{g}"' for g in GROUPS) + "};",
        "constexpr const char *CFG_SUB_NAMES[] = {"
        + ", ".join(f'"{x}"' for x in SUBS) + "};",
        "",
        "constexpr F1Seg CFG_SEGS[CFG_SEG_COUNT] = {",
    ]
    for r in s.rows:
        out.append("  {%d, %d, %d, %d, %d, %dUL, %dUL}," % (
            TYPES[r.seg], GROUPS.index(r.group), SUBS.index(r.sub), r.row,
            r.u_target, round(r.rate_hz * 1000.0), round(r.dwell_s * 1000.0)))
    out[-1] = out[-1].rstrip(",")
    out += ["};", ""]
    path.write_text("\n".join(out), encoding="utf-8")


def main() -> int:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    cfg_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    s, r2 = plan(c)
    ok = checks(s, r2)
    summary = summarise(s, r2, ok)

    SEQ.mkdir(exist_ok=True)
    table_hash = write_table(SEQ / "f1_command_table.csv", s.rows)
    summary["config_sha256"] = cfg_hash
    summary["table_sha256"] = table_hash
    (SEQ / "f1_table_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    write_header(HEADER, c, s, r2, cfg_hash, table_hash)

    print(f"F1 table: {summary['segments']} segments, "
          f"{summary['duration_min']:.1f} min "
          f"(R2 {'included' if r2 else 'dropped'}; gate {c['gate_min']:.0f} min), "
          f"{summary['margin_to_gate_min']:.1f} min of margin")
    for g, m in summary["groups_min"].items():
        print(f"  {g:3s} {m:5.2f} min")
    print(f"  ramp rates {summary['slowest_ramp_hz']:.2f} .. "
          f"{summary['fastest_ramp_hz']:.0f} Hz "
          f"(driver limit {MAX_STEP_RATE_HZ} Hz)")
    print(f"  PROVISIONAL grid from s90={summary['s90_um']} um, "
          f"d_min={summary['d_min_um']} um: "
          f"{len(summary['strokes_usteps'])} strokes "
          f"{summary['strokes_usteps'][0]}..{summary['strokes_usteps'][-1]} usteps")
    bad = [k for k, v in ok.items() if not v]
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    print(f"wrote {HEADER.relative_to(HERE)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
