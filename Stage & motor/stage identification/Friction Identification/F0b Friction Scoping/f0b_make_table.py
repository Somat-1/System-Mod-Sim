#!/usr/bin/env python3
"""Build the F0b command table and the firmware header from f0b_sequence_config.json.

Every row is one segment the board replays: a trapezoidal transit, a
constant-rate ramp, a run of single microsteps with a dwell after each, a
dwell in place, or a recorded reference mark. The host owns the arithmetic;
the firmware owns only the timing.

    python3 f0b_make_table.py

writes

    sequence/f0b_command_table.csv          the table, procedure section 2.1 format
    sequence/f0b_plan_summary.json          durations, counts and the checks
    firmware/f0b_friction_scoping/f0b_config.h   the generated firmware header

and exits non-zero if any check fails, so it can gate the session.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f0b_sequence_config.json"
RENDER_CONFIG = HERE / "f0b_config.json"
SEQ = HERE / "sequence"
HEADER = HERE / "firmware" / "f0b_friction_scoping" / "f0b_config.h"

# Segment type codes, mirrored by the enum in the firmware header.
TYPES = {"move": 0, "ramp": 1, "step": 2, "dwell": 3, "mark": 4}
BLOCKS = ["R1", "N", "L", "R2", "S", "D", "V", "R3"]


@dataclass
class Row:
    row: int
    block: str
    seg: str
    u_target: int
    rate_hz: float
    dwell_s: float
    tag: str
    t0: float
    t1: float


class Session:
    """Accumulates the segment list and the clock as the sequence is built."""

    def __init__(self, c: dict) -> None:
        self.c = c
        self.rows: list[Row] = []
        self.t = 0.0
        self.u = round(c["park_mm"] * c["ustep_per_mm"])
        self.block = "R1"

    # -- geometry helpers -------------------------------------------------
    def u_of(self, mm: float) -> int:
        return round(mm * self.c["ustep_per_mm"])

    def mm_of(self, u: int) -> float:
        return u / self.c["ustep_per_mm"]

    def profile_time(self, n_usteps: int, rate_hz: float) -> float:
        """Trapezoidal time for n microsteps peaking at rate_hz.

        The firmware runs every motion segment through one generator: it
        accelerates at CFG_A_MM_S2 towards the segment's rate and decelerates
        at the far end. Below about 1 kHz the acceleration phase is shorter
        than one microstep, so a slow measurement ramp is constant-rate in
        practice; at 19.2 and 25.6 kHz it is a real trapezoid, which is what
        keeps the stage from being asked to start at 0.75 or 1 mm/s from rest.
        """
        n = abs(n_usteps)
        if n == 0 or rate_hz <= 0.0:
            return 0.0
        a = self.c["acceleration_mm_s2"] * self.c["ustep_per_mm"]
        ramp = rate_hz * rate_hz / a
        if n <= ramp:
            return 2.0 * math.sqrt(n / a)
        return 2.0 * rate_hz / a + (n - ramp) / rate_hz

    def transit_rate(self) -> float:
        return self.c["transit_mm_s"] * self.c["ustep_per_mm"]

    def transit_time(self, d_usteps: int) -> float:
        return self.profile_time(d_usteps, self.transit_rate())

    # -- segment emitters -------------------------------------------------
    def _add(self, seg: str, u_target: int, rate_hz: float, dwell_s: float,
             tag: str, duration: float) -> None:
        self.rows.append(Row(len(self.rows), self.block, seg, u_target,
                             rate_hz, dwell_s, tag, self.t, self.t + duration))
        self.t += duration
        self.u = u_target

    def move(self, u_target: int, tag: str) -> None:
        self._add("move", u_target, 0.0, 0.0, tag,
                  self.transit_time(u_target - self.u))

    def move_mm(self, mm: float, tag: str) -> None:
        self.move(self.u_of(mm), tag)

    def ramp(self, u_target: int, rate_hz: float, tag: str) -> None:
        self._add("ramp", u_target, rate_hz, 0.0, tag,
                  self.profile_time(u_target - self.u, rate_hz))

    def steps(self, u_target: int, dwell_s: float, tag: str) -> None:
        n = abs(u_target - self.u)
        self._add("step", u_target, 0.0, dwell_s, tag, n * dwell_s)

    def dwell(self, seconds: float, tag: str) -> None:
        self._add("dwell", self.u, 0.0, seconds, tag, seconds)

    def mark(self, tag: str) -> None:
        """Settle, then the recorded window that defines a reference stop."""
        self.dwell(self.c["settle_s"], tag + " settle")
        self._add("mark", self.u, 0.0, self.c["record_s"], tag,
                  self.c["record_s"])

    def wipe(self, u_target: int, from_below: bool, tag: str) -> None:
        """Standard memory wipe: 3 mm run-up in one direction, then settle."""
        off = round(self.c["wipe_mm"] * self.c["ustep_per_mm"])
        self.move(u_target - off if from_below else u_target + off,
                  tag + " run-up")
        self.move(u_target, tag + " wipe " +
                  ("from below" if from_below else "from above"))
        self.dwell(self.c["settle_s"], tag + " settle")


# -- blocks ---------------------------------------------------------------
def reference_block(s: Session, name: str) -> None:
    """15 / 105 mm, both directions, nearest position first (as E5)."""
    s.block = name
    c = s.c
    pre = round(c["reference_preapproach_mm"] * c["ustep_per_mm"])
    targets = sorted((s.u_of(mm) for mm in c["reference_positions_mm"]),
                     key=lambda u: abs(u - s.u))
    for u in targets:
        for side in (-1, 1):
            s.move(u + side * pre, f"{s.mm_of(u):.0f} mm approach")
            s.move(u, f"{s.mm_of(u):.0f} mm " +
                   ("from below" if side < 0 else "from above"))
            s.mark(f"{s.mm_of(u):.0f} mm " +
                   ("up" if side < 0 else "down"))


def block_n(s: Session) -> None:
    """Noise floor at three positions spanning the travel (procedure section 3)."""
    s.block = "N"
    for hold in s.c["n_holds"]:
        u = s.u_of(hold["position_mm"])
        s.wipe(u, True, hold["tag"])
        s.dwell(hold["hold_s"], hold["tag"] + " hold")


def block_l(s: Session) -> None:
    """Coarse reversal ladder: anchored loops above x0, increasing then decreasing."""
    s.block = "L"
    c = s.c
    x0 = s.u_of(c["x0_mm"])
    strokes = list(c["l_strokes_usteps"])
    s.wipe(x0, True, "L")
    for pass_name, order in (("inc", strokes), ("dec", strokes[::-1])):
        for stroke in order:
            tag = f"L s={stroke} {pass_name}"
            step_mode = stroke <= c["l_step_segment_max_usteps"]
            rate = min(stroke / c["l_half_cycle_s"], c["l_rate_cap_hz"])
            for half in range(c["l_halves_per_stroke"]):
                up = half % 2 == 0
                target = x0 + stroke if up else x0
                label = f"{tag} h{half + 1} {'up' if up else 'down'}"
                if step_mode:
                    s.steps(target, c["l_step_dwell_s"], label)
                else:
                    s.ramp(target, rate, label)
                s.dwell(c["turn_pause_s"], f"{tag} turn")


def block_s(s: Session) -> None:
    """Slow ramps: 20 um at 0.25 um/s, four preparations, twice."""
    s.block = "S"
    c = s.c
    x0 = s.u_of(c["x0_mm"])
    span = round(c["s_ramp_um"] * 1e-3 * c["ustep_per_mm"])
    rate = c["s_rate_um_s"] * 1e-3 * c["ustep_per_mm"]
    for rep in range(c["s_repetitions"]):
        for spec in c["s_ramps"]:
            tag = f"{spec['tag']} r{rep + 1}"
            s.wipe(x0, spec["wipe_from"] == "below", tag)
            s.ramp(x0 + spec["direction"] * span, rate,
                   f"{tag} {spec['type']} ramp")


def block_d(s: Session) -> None:
    """Dwell probe: wipe, dwell T_d, then a 5 um reversal ramp. ABBA order."""
    s.block = "D"
    c = s.c
    x0 = s.u_of(c["x0_mm"])
    span = round(c["d_probe_um"] * 1e-3 * c["ustep_per_mm"])
    rate = c["d_rate_um_s"] * 1e-3 * c["ustep_per_mm"]
    for i, td in enumerate(c["d_dwells_s"]):
        tag = f"D{i + 1} Td={td:g}s"
        s.wipe(x0, True, tag)
        s.dwell(td, tag + " dwell")
        s.ramp(x0 - span, rate, tag + " reversal ramp")


def block_v(s: Session) -> None:
    """Speed probe: 1 mm up and down, fast-fast-slow-fast-fast."""
    s.block = "V"
    c = s.c
    x0 = s.u_of(c["x0_mm"])
    top = x0 + round(c["v_span_mm"] * c["ustep_per_mm"])
    s.wipe(x0, True, "V")
    for i, speed in enumerate(c["v_pass_order"]):
        mm_s = c["v_fast_mm_s"] if speed == "fast" else c["v_slow_mm_s"]
        rate = mm_s * c["ustep_per_mm"]
        s.ramp(top, rate, f"V{i + 1} {speed} up")
        s.dwell(c["turn_pause_s"], f"V{i + 1} turn")
        s.ramp(x0, rate, f"V{i + 1} {speed} down")
        s.dwell(c["turn_pause_s"], f"V{i + 1} turn")


def build(c: dict) -> Session:
    s = Session(c)
    reference_block(s, "R1")
    block_n(s)
    block_l(s)
    reference_block(s, "R2")
    block_s(s)
    block_d(s)
    block_v(s)
    reference_block(s, "R3")
    return s


# -- checks ---------------------------------------------------------------
def checks(s: Session) -> dict[str, bool]:
    c = s.c
    x0 = s.u_of(c["x0_mm"])
    u_min, u_max = s.u_of(c["soft_min_mm"]), s.u_of(c["soft_max_mm"])
    strokes = c["l_strokes_usteps"]
    ustep_um = 1e3 / c["ustep_per_mm"]
    rates = [r.rate_hz for r in s.rows if r.seg == "ramp"]
    ladder = [r for r in s.rows if r.block == "L" and r.seg == "ramp"]
    return {
        "x0_on_electrical_grid": x0 % c["grid_usteps"] == 0,
        "x0_matches_60mm": x0 == 1_536_000,
        "inside_soft_limits": all(u_min <= r.u_target <= u_max for r in s.rows),
        "under_driver_limit": all(r <= c["max_step_rate_hz"] + 1e-9 for r in rates),
        "ladder_under_cap": all(r.rate_hz <= c["l_rate_cap_hz"] + 1e-9
                                for r in ladder),
        "strokes_monotonic": all(b > a for a, b in zip(strokes, strokes[1:])),
        "ustep_above_floor": ustep_um > 5.0 * c["interferometer_floor_um"],
        "ladder_brackets_beta0": (strokes[0] * ustep_um < c["beta_0_um"]
                                  < strokes[-1] * ustep_um),
        "inside_watchdog_gate": s.t / 60.0 <= c["watchdog_gate_min"],
        "inside_scope_record": s.t / 60.0 <= c["scope_record_min"],
    }


def summarise(s: Session, ok: dict[str, bool]) -> dict:
    by_block, by_seg = {}, {}
    for r in s.rows:
        by_block[r.block] = by_block.get(r.block, 0.0) + (r.t1 - r.t0)
        by_seg[r.seg] = by_seg.get(r.seg, 0.0) + (r.t1 - r.t0)
    return {
        "segments": len(s.rows),
        "duration_min": s.t / 60.0,
        "gate_min": s.c["watchdog_gate_min"],
        "margin_to_gate_min": s.c["watchdog_gate_min"] - s.t / 60.0,
        "blocks_min": {b: by_block[b] / 60.0 for b in BLOCKS if b in by_block},
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
        w.writerow(["row", "block", "seg", "u_target", "rate_hz", "dwell_s",
                    "tag", "t0_s", "t1_s"])
        for r in rows:
            w.writerow([r.row, r.block, r.seg, r.u_target, f"{r.rate_hz:.6f}",
                        f"{r.dwell_s:.3f}", r.tag, f"{r.t0:.3f}", f"{r.t1:.3f}"])
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_header(path: Path, c: dict, s: Session, cfg_hash: str,
                 table_hash: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    marks = sum(r.seg == "mark" for r in s.rows)
    out = [
        "/* GENERATED by f0b_make_table.py. Do not edit by hand. */",
        "#pragma once",
        "#include <stdint.h>",
        "",
        f'constexpr char CFG_HASH[] = "{cfg_hash}";',
        f'constexpr char CFG_TABLE_HASH[] = "{table_hash}";',
        f'constexpr char CFG_EXPERIMENT[] = "{c["experiment"]}";',
        f"constexpr double CFG_LEAD_MM = {c['lead_mm']:.6f};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {c['steps_per_rev']};",
        f"constexpr uint16_t CFG_MRES = {c['mres']};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {c['run_current_ma']};",
        f"constexpr uint16_t CFG_HOLD_CURRENT_MA = {c['hold_current_ma']};",
        f"constexpr uint8_t CFG_TPOWERDOWN = {c['tpowerdown']};",
        f"constexpr uint8_t CFG_IHOLDDELAY = {c['iholddelay']};",
        f"constexpr double CFG_V_MM_S = {c['transit_mm_s']:.6f};",
        f"constexpr double CFG_A_MM_S2 = {c['acceleration_mm_s2']:.6f};",
        f"constexpr uint32_t CFG_MAX_STEP_RATE_HZ = {c['max_step_rate_hz']}UL;",
        f"constexpr uint32_t CFG_SETTLE_MS = {round(c['settle_s'] * 1000)}UL;",
        f"constexpr uint32_t CFG_RECORD_MS = {round(c['record_s'] * 1000)}UL;",
        f"constexpr int32_t CFG_PARK_U = {s.u_of(c['park_mm'])};",
        f"constexpr int32_t CFG_X0_U = {s.u_of(c['x0_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {s.u_of(c['soft_min_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {s.u_of(c['soft_max_mm'])};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {c['motor_cycle_usteps']};",
        f"constexpr uint32_t CFG_SCREW_REV_U = {c['screw_rev_usteps']}UL;",
        f"constexpr uint32_t CFG_TIME_LIMIT_MS = {round(c['time_limit_min'] * 60000)}UL;",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['watchdog_gate_min'] * 60000)}UL;",
        f"constexpr uint16_t CFG_AQUADB_CLOCK_NS = {c['aquadb_clock_ns']};",
        f"constexpr uint32_t CFG_RMT_HZ = {c['rmt_resolution_hz']}UL;",
        f"constexpr uint16_t CFG_RMT_CHUNK = {c['rmt_chunk_pulses']};",
        f"constexpr uint8_t CFG_RMT_QUEUE_DEPTH = {c['rmt_queue_depth']};",
        f"constexpr uint8_t CFG_RMT_RING = {c['rmt_ring_buffers']};",
        f"constexpr uint32_t CFG_RMT_MAX_PERIOD_US = {c['rmt_max_period_us']}UL;",
        f"constexpr int8_t CFG_SYNC_GPIO = {c['sync_gpio']};",
        f"constexpr uint32_t CFG_PLAN_DURATION_MS = {round(s.t * 1000)}UL;",
        "",
        "enum F0bSegType : uint8_t { F0B_MOVE = 0, F0B_RAMP = 1, F0B_STEP = 2,",
        "                            F0B_DWELL = 3, F0B_MARK = 4 };",
        "",
        "/* rate_mhz is the step rate in milli-hertz: 6400 means 6.4 steps/s.",
        "   dwell_ms is the dwell for F0B_DWELL and F0B_MARK, and the dwell",
        "   after each single microstep for F0B_STEP. */",
        "struct F0bSeg {",
        "  uint8_t type;",
        "  uint8_t block;",
        "  uint16_t row;",
        "  int32_t u_target;",
        "  uint32_t rate_mhz;",
        "  uint32_t dwell_ms;",
        "};",
        "",
        f"constexpr uint16_t CFG_SEG_COUNT = {len(s.rows)};",
        f"constexpr uint16_t CFG_MARK_COUNT = {marks};",
        "constexpr const char *CFG_BLOCK_NAMES[] = {"
        + ", ".join(f'"{b}"' for b in BLOCKS) + "};",
        "",
        "constexpr F0bSeg CFG_SEGS[CFG_SEG_COUNT] = {",
    ]
    for r in s.rows:
        out.append(
            "  {%d, %d, %d, %d, %dUL, %dUL}," % (
                TYPES[r.seg], BLOCKS.index(r.block), r.row, r.u_target,
                round(r.rate_hz * 1000.0), round(r.dwell_s * 1000.0)))
    out[-1] = out[-1].rstrip(",")
    out += ["};", ""]
    path.write_text("\n".join(out), encoding="utf-8")


def cross_check_render_config(s: Session, summary: dict) -> list[str]:
    """Warn when the renderers' presentation values drift from the real plan."""
    if not RENDER_CONFIG.exists():
        return []
    r = json.loads(RENDER_CONFIG.read_text(encoding="utf-8"))
    warn = []
    if abs(r.get("duration_min", 0) - summary["duration_min"]) > 0.3:
        warn.append(f"f0b_config.json duration_min={r.get('duration_min')} "
                    f"but the plan is {summary['duration_min']:.1f}")
    if r.get("segments") not in (None, summary["segments"]):
        warn.append(f"f0b_config.json segments={r.get('segments')} "
                    f"but the plan has {summary['segments']}")
    for b in r.get("blocks", []):
        got = summary["blocks_min"].get(b["name"])
        if got is not None and abs(got - b["duration_min"]) > 0.3:
            warn.append(f"block {b['name']}: f0b_config.json says "
                        f"{b['duration_min']:.1f} min, plan says {got:.1f}")
    return warn


def main() -> int:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    cfg_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    s = build(c)
    ok = checks(s)
    summary = summarise(s, ok)

    SEQ.mkdir(exist_ok=True)
    table_hash = write_table(SEQ / "f0b_command_table.csv", s.rows)
    summary["config_sha256"] = cfg_hash
    summary["table_sha256"] = table_hash
    (SEQ / "f0b_plan_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    write_header(HEADER, c, s, cfg_hash, table_hash)

    print(f"F0b plan: {summary['segments']} segments, "
          f"{summary['duration_min']:.1f} min, "
          f"{summary['margin_to_gate_min']:.1f} min inside the "
          f"{c['watchdog_gate_min']:.0f} min gate")
    for b, m in summary["blocks_min"].items():
        print(f"  {b:3s} {m:5.2f} min")
    print(f"  ramp rates {summary['slowest_ramp_hz']:.2f} .. "
          f"{summary['fastest_ramp_hz']:.0f} Hz "
          f"(driver limit {c['max_step_rate_hz']} Hz)")
    for w in cross_check_render_config(s, summary):
        print(f"  WARNING: {w}")
    bad = [k for k, v in ok.items() if not v]
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    print(f"wrote {HEADER.relative_to(HERE)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
