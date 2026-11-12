#!/usr/bin/env python3
"""Build the F0b_pt2 command table and the firmware header it is flashed as.

    python3 f0b_pt2_make_table.py

writes

    sequence/f0b_pt2_command_table.csv   one row per segment, in execution order
    sequence/f0b_pt2_arrivals.csv        the arrivals the analysis reads
    sequence/f0b_pt2_summary.json        counts, durations, hashes, checks
    firmware/f0b_pt2/f0b_pt2_config.h    generated; do not edit

THE ONE INTERNAL UNIT IS 1/256 COUNTS (procedure section 0, firmware change 1).
Every position in the table is in those counts whatever MRES is in force, and
every rate is in milli-counts per second, so a velocity means the same thing
at every MRES and only the step rate changes. A segment is legal at MRES
1/m only if its travel is divisible by 256/m counts; this builder asserts
that for every segment against the MRES actually in force at that point in
the sequence, and the firmware asserts it again at run time.

WINDOWS ARE OFFSETS FROM THE LAST STEP EDGE (firmware change 3). dwell_ms in
the table is not a duration from when the segment was issued but the offset
from t_last_edge at which the segment ends, so the settle ends 2.000 s after
the final pulse and the mark 2.500 s after it, whatever the UART and the
scheduler did in between.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f0b_pt2_config.json"
SEQ = HERE / "sequence"
FIRMWARE = HERE / "firmware" / "f0b_pt2"

MOVE, RAMP, DWELL, MARK, SETMRES = range(5)
TYPE_NAMES = ["move", "ramp", "dwell", "mark", "setmres"]
SIDE_NAMES = ["none", "below", "above"]


def mres_code(mres: int) -> int:
    """CHOPCONF.MRES field: 0 = 1/256, 1 = 1/128 ... 4 = 1/16."""
    code = {256: 0, 128: 1, 64: 2, 32: 3, 16: 4}.get(mres)
    if code is None:
        raise ValueError(f"unsupported MRES 1/{mres}")
    return code


def counts_per_step(mres: int) -> int:
    return 256 // mres


@dataclass
class Seg:
    row: int
    block: str
    visit: int
    kind: int
    mres: int
    side: int
    u_target: int          # absolute, 1/256 counts
    rate_mcps: int         # milli-counts per second; 0 for dwell/mark/setmres
    end_offset_ms: int     # for dwell/mark: offset from t_last_edge
    t0_ms: float = 0.0
    t1_ms: float = 0.0
    tag: str = ""


@dataclass
class Session:
    c: dict
    segs: list = field(default_factory=list)
    t: float = 0.0                 # ms
    u: int = 0                     # 1/256 counts
    mres: int = 256
    mres_set: bool = False
    block: str = "R1"
    visit: int = 0
    arrivals: list = field(default_factory=list)

    # -- helpers ---------------------------------------------------------
    def u_of_mm(self, mm: float) -> int:
        counts = mm * self.c["ustep_per_mm"]
        out = int(round(counts))
        assert abs(counts - out) < 1e-6, f"{mm} mm is not a whole count"
        return out

    def move_time_ms(self, du: int, mm_s: float) -> float:
        """Trapezoid at mm_s with accel_mm_s2, in counts."""
        d = abs(du)
        if d == 0:
            return 0.0
        a = self.c["accel_mm_s2"] * self.c["ustep_per_mm"]
        v = mm_s * self.c["ustep_per_mm"]
        if d <= v * v / a:
            return 2000.0 * math.sqrt(d / a)
        return 1000.0 * (d / v + v / a)

    def add(self, kind: int, target: int, rate_mm_s: float = 0.0,
            end_offset_ms: int = 0, side: int = 0, tag: str = "",
            duration_ms: float | None = None) -> Seg:
        if kind in (MOVE, RAMP):
            step = counts_per_step(self.mres)
            delta = target - self.u
            assert delta % step == 0, (
                f"row {len(self.segs)}: {delta} counts is not a whole number "
                f"of 1/{self.mres} steps ({step} counts each)")
        rate_mcps = int(round(rate_mm_s * self.c["ustep_per_mm"] * 1000))
        if duration_ms is None:
            duration_ms = (self.move_time_ms(target - self.u, rate_mm_s)
                           if kind in (MOVE, RAMP) else float(end_offset_ms))
        seg = Seg(len(self.segs), self.block, self.visit, kind, self.mres,
                  side, target, rate_mcps if kind in (MOVE, RAMP) else 0,
                  end_offset_ms, self.t, self.t + duration_ms, tag)
        self.segs.append(seg)
        self.t += duration_ms
        self.u = target
        return seg

    def go(self, target_u: int, mm_s: float, tag: str, side: int = 0,
           kind: int = RAMP) -> None:
        self.add(kind, target_u, rate_mm_s=mm_s, side=side, tag=tag)

    def stop(self, side: int, tag: str) -> None:
        """A settle and a record window, both measured from the last edge."""
        c = self.c
        self.add(DWELL, self.u, end_offset_ms=c["settle_ms"], side=side,
                 tag=tag + " settle", duration_ms=c["settle_ms"])
        self.add(MARK, self.u, end_offset_ms=c["settle_ms"] + c["record_ms"],
                 side=side, tag=tag + " record", duration_ms=c["record_ms"])
        self.arrivals.append({
            "row": len(self.segs) - 1, "block": self.block,
            "visit": self.visit, "mres": self.mres,
            "side": SIDE_NAMES[side], "u": self.u,
            "mm": self.u / self.c["ustep_per_mm"], "tag": tag,
            "t_ms": round(self.t, 1)})

    def set_mres(self, mres: int, tag: str, force: bool = False) -> None:
        """Always emitted the first time, so the table never relies on the
        resolution the board happened to be flashed at, and always at the
        head of an M visit even when the level does not change: the switch
        carries the MSCNT check, and procedure section 3 asks for it once
        per visit, including the repeat at the same level.
        """
        if mres == self.mres and self.mres_set and not force:
            return
        self.mres, self.mres_set = mres, True
        self.add(SETMRES, self.u, end_offset_ms=0, tag=tag, duration_ms=60.0)


# -- blocks ---------------------------------------------------------------
def reference_pair(s: Session, name: str) -> None:
    """15 / 105 mm, both directions, at the reference MRES and E5 speed."""
    c = s.c
    s.block, s.visit = name, 0
    s.set_mres(c["reference_mres"], f"{name} to 1/{c['reference_mres']}")
    over = s.u_of_mm(c["reference_overrun_mm"])
    speed = c["reference_speed_mm_s"]
    for mm in sorted(c["reference_mm"], key=lambda r: abs(s.u_of_mm(r) - s.u)):
        target = s.u_of_mm(mm)
        for side in (1, 2):                      # from below, then from above
            launch = target - over if side == 1 else target + over
            s.go(launch, c["transit_mm_s"], f"{name} {mm:g} mm launch",
                 kind=MOVE)
            s.go(target, speed, f"{name} {mm:g} mm approach", side=side)
            s.stop(side, f"{name} {mm:g} mm from {SIDE_NAMES[side]}")


def mres_visit(s: Session, mres: int, visit: int) -> None:
    """One MRES level: saturated reversals, then two partial-stroke loops."""
    c = s.c
    s.block, s.visit = "M", visit
    x0 = s.u_of_mm(c["x0_mm"])
    speed = c["m_speed_mm_s"]
    over = s.u_of_mm(c["m_overrun_mm"])

    s.set_mres(mres, f"M visit {visit} to 1/{mres}", force=True)

    # ISO-style chain: 7 arrivals at x0, alternating below / above
    for k in range(c["m_approaches"]):
        side = 1 if k % 2 == 0 else 2
        launch = x0 - over if side == 1 else x0 + over
        s.go(launch, speed, f"M{visit} overrun {k + 1}")
        s.go(x0, speed, f"M{visit} approach {k + 1}", side=side)
        s.stop(side, f"M{visit} saturated {k + 1}")

    # partial strokes, anchored at x0 approached from below
    s.go(x0 - over, speed, f"M{visit} wipe launch")
    s.go(x0, speed, f"M{visit} wipe", side=1)
    s.stop(1, f"M{visit} anchor")
    for stroke_um in c["m_partial_strokes_um"]:
        top = x0 + s.u_of_mm(stroke_um / 1000.0)
        for cycle in range(c["m_partial_cycles"]):
            s.go(top, speed, f"M{visit} {stroke_um:g} um up c{cycle + 1}", side=1)
            s.stop(1, f"M{visit} {stroke_um:g} um top c{cycle + 1}")
            s.go(x0, speed, f"M{visit} {stroke_um:g} um down c{cycle + 1}", side=2)
            s.stop(2, f"M{visit} {stroke_um:g} um bottom c{cycle + 1}")


def speed_visit(s: Session, mm_s: float, mres: int, visit: int) -> None:
    """One approach speed: a wipe, then a chain of saturated arrivals."""
    c = s.c
    s.block, s.visit = "V", visit
    x0 = s.u_of_mm(c["x0_mm"])
    s.set_mres(mres, f"V visit {visit} to 1/{mres}")
    wipe = s.u_of_mm(c["v_wipe_mm"])
    over = s.u_of_mm(c["v_overrun_um"] / 1000.0)
    s.go(x0 - wipe, c["transit_mm_s"], f"V{visit} wipe launch", kind=MOVE)
    s.go(x0, c["transit_mm_s"], f"V{visit} wipe", side=1, kind=MOVE)
    for k in range(c["v_approaches"]):
        side = 1 if k % 2 == 0 else 2
        launch = x0 - over if side == 1 else x0 + over
        s.go(launch, mm_s, f"V{visit} overrun {k + 1}")
        s.go(x0, mm_s, f"V{visit} approach {k + 1}", side=side)
        s.stop(side, f"V{visit} {mm_s:g} mm/s {k + 1}")


def build(c: dict) -> Session:
    s = Session(c)
    s.u = s.u_of_mm(c["park_mm"])
    s.mres = c["reference_mres"]           # the board is flashed at 1/16
    reference_pair(s, "R1")
    for visit, mres in enumerate(c["m_visit_order"], start=1):
        mres_visit(s, mres, visit)
    speeds = list(c["v_speed_order_mm_s"])
    for visit, mm_s in enumerate(speeds, start=1):
        speed_visit(s, mm_s, c["v_mres"], visit)
    if c["v_extra_coarse"]:
        speed_visit(s, 2.0, c["reference_mres"], len(speeds) + 1)
    reference_pair(s, "R3")
    return s


# -- checks ---------------------------------------------------------------
def checks(s: Session) -> dict:
    c = s.c
    x0 = s.u_of_mm(c["x0_mm"])
    cycle = c["motor_cycle_u"]
    coarse = counts_per_step(min(c["mres_levels"]))

    turning = {seg.u_target for seg in s.segs if seg.kind in (MOVE, RAMP)}
    saturated = {x0 - s.u_of_mm(c["m_overrun_mm"]),
                 x0 + s.u_of_mm(c["m_overrun_mm"]), x0}

    # every segment divisible at the MRES in force, re-derived from the table
    divisible, here = True, s.u_of_mm(c["park_mm"])
    for seg in s.segs:
        if seg.kind in (MOVE, RAMP):
            if (seg.u_target - here) % counts_per_step(seg.mres):
                divisible = False
            here = seg.u_target

    # every approach's direction matches its declared side
    sides_ok, prev = True, s.u_of_mm(c["park_mm"])
    for seg in s.segs:
        if seg.kind in (MOVE, RAMP) and seg.side:
            rising = seg.u_target > prev
            if (seg.side == 1) != rising:
                sides_ok = False
        if seg.kind in (MOVE, RAMP):
            prev = seg.u_target

    rates = [seg.rate_mcps / 1e3 / counts_per_step(seg.mres)
             for seg in s.segs if seg.rate_mcps]
    lo, hi = c["constrained_span_mm"]
    return {
        "all_targets_divisible_at_their_mres": divisible,
        "all_turning_points_on_the_coarse_grid": all(
            u % coarse == 0 for u in turning),
        "saturated_points_on_the_electrical_grid": all(
            u % cycle == 0 for u in saturated),
        "approach_sides_match_direction": sides_ok,
        "inside_soft_limits": all(
            s.u_of_mm(c["soft_min_mm"]) <= u <= s.u_of_mm(c["soft_max_mm"])
            for u in turning),
        "test_positions_inside_constrained_span": all(
            lo <= u / c["ustep_per_mm"] <= hi for u in saturated),
        "max_step_rate_within_bench": max(rates) <= c["max_step_rate_hz"] + 1,
        "mres_switches_at_rest": all(
            seg.kind != SETMRES or i == 0 or s.segs[i - 1].kind in (DWELL, MARK, SETMRES)
            for i, seg in enumerate(s.segs)),
        "inside_gate": s.t / 60000.0 <= c["gate_min"],
        "gate_inside_watchdog": c["gate_min"] <= c["watchdog_min"] <= c["scope_record_min"],
        "every_mres_level_visited_twice": all(
            c["m_visit_order"].count(m) == 2 for m in c["mres_levels"]),
        "mres_order_is_abba": c["m_visit_order"] == c["m_visit_order"][::-1],
        "speed_order_is_abba": (list(c["v_speed_order_mm_s"])
                                == list(c["v_speed_order_mm_s"])[::-1]),
    }


# -- outputs ---------------------------------------------------------------
def write_table(path: Path, s: Session) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["row", "block", "visit", "type", "mres", "side",
                    "u_target", "rate_mcps", "end_offset_ms", "t0_ms", "t1_ms",
                    "tag"])
        for seg in s.segs:
            w.writerow([seg.row, seg.block, seg.visit, TYPE_NAMES[seg.kind],
                        seg.mres, SIDE_NAMES[seg.side], seg.u_target,
                        seg.rate_mcps, seg.end_offset_ms,
                        f"{seg.t0_ms:.1f}", f"{seg.t1_ms:.1f}", seg.tag])
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_arrivals(path: Path, s: Session) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(s.arrivals[0]))
        w.writeheader()
        w.writerows(s.arrivals)


def write_header(path: Path, c: dict, s: Session, cfg_hash: str,
                 table_hash: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    blocks = list(dict.fromkeys(seg.block for seg in s.segs))
    out = [
        "/* GENERATED by f0b_pt2_make_table.py. Do not edit by hand. */",
        "#pragma once",
        "#include <stdint.h>",
        "",
        f'constexpr char CFG_HASH[] = "{cfg_hash}";',
        f'constexpr char CFG_TABLE_HASH[] = "{table_hash}";',
        f'constexpr char CFG_EXPERIMENT[] = "{c["experiment"]}";',
        f"constexpr double CFG_LEAD_MM = {c['lead_mm']:.6f};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {c['steps_per_rev']};",
        f"constexpr uint16_t CFG_START_MRES = {c['reference_mres']};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {c['run_current_ma']};",
        f"constexpr uint16_t CFG_HOLD_CURRENT_MA = {c['hold_current_ma']};",
        f"constexpr uint8_t CFG_TPOWERDOWN = {c['tpowerdown']};",
        f"constexpr uint8_t CFG_IHOLDDELAY = {c['iholddelay']};",
        f"constexpr double CFG_V_MM_S = {c['transit_mm_s']:.6f};",
        f"constexpr double CFG_A_MM_S2 = {c['accel_mm_s2']:.6f};",
        f"constexpr uint32_t CFG_MAX_STEP_RATE_HZ = {c['max_step_rate_hz']}UL;",
        f"constexpr uint32_t CFG_SETTLE_MS = {c['settle_ms']}UL;",
        f"constexpr uint32_t CFG_RECORD_MS = {c['record_ms']}UL;",
        f"constexpr int32_t CFG_PARK_U = {s.u_of_mm(c['park_mm'])};",
        f"constexpr int32_t CFG_X0_U = {s.u_of_mm(c['x0_mm'])};",
        f"constexpr int32_t CFG_SOFT_MIN_U = {s.u_of_mm(c['soft_min_mm'])};",
        f"constexpr int32_t CFG_SOFT_MAX_U = {s.u_of_mm(c['soft_max_mm'])};",
        f"constexpr uint16_t CFG_MOTOR_CYCLE_U = {c['motor_cycle_u']};",
        f"constexpr uint32_t CFG_SCREW_REV_U = {c['screw_rev_u']}UL;",
        f"constexpr uint32_t CFG_WATCHDOG_MS = {round(c['watchdog_min'] * 60000)}UL;",
        f"constexpr uint32_t CFG_RMT_MAX_PERIOD_US = {c['rmt_max_period_us']}UL;",
        "constexpr uint32_t CFG_RMT_HZ = 1000000UL;",
        "constexpr uint16_t CFG_RMT_CHUNK = 128;",
        "constexpr uint8_t CFG_RMT_RING = 6;",
        "constexpr int8_t CFG_SYNC_GPIO = -1;",
        f"constexpr uint32_t CFG_PLAN_DURATION_MS = {round(s.t)}UL;",
        "",
        "enum Pt2SegType : uint8_t { PT2_MOVE = 0, PT2_RAMP = 1,",
        "                           PT2_DWELL = 2, PT2_MARK = 3,",
        "                           PT2_SETMRES = 4 };",
        "",
        "/* u_target is in 1/256 counts, the one internal unit, whatever MRES",
        "   is in force; the firmware converts to driver steps by dividing by",
        "   (1 << mres_code) and refuses a segment that does not divide.",
        "   rate_mcps is milli-counts per second, so a velocity is the same",
        "   number at every MRES.  end_offset_ms is measured from the last",
        "   step edge, not from when the segment was issued: a DWELL ends at",
        "   t_last_edge + settle, a MARK at t_last_edge + settle + record.",
        "   side is 0 none, 1 arriving from below, 2 from above. */",
        "struct Pt2Seg {",
        "  uint8_t type;",
        "  uint8_t block;",
        "  uint8_t mres_code;",
        "  uint8_t side;",
        "  uint16_t row;",
        "  uint16_t visit;",
        "  int32_t u_target;",
        "  uint32_t rate_mcps;",
        "  uint32_t end_offset_ms;",
        "};",
        "",
        f"constexpr uint16_t CFG_SEG_COUNT = {len(s.segs)};",
        f"constexpr uint16_t CFG_MARK_COUNT = "
        f"{sum(1 for q in s.segs if q.kind == MARK)};",
        "constexpr const char *CFG_BLOCK_NAMES[] = {"
        + ", ".join(f'"{b}"' for b in blocks) + "};",
        "constexpr const char *CFG_SIDE_NAMES[] = {\"none\", \"below\", \"above\"};",
        "",
        "constexpr Pt2Seg CFG_SEGS[CFG_SEG_COUNT] = {",
    ]
    index = {b: i for i, b in enumerate(blocks)}
    for seg in s.segs:
        out.append(
            f"  {{{seg.kind}, {index[seg.block]}, {mres_code(seg.mres)}, "
            f"{seg.side}, {seg.row}, {seg.visit}, {seg.u_target}, "
            f"{seg.rate_mcps}UL, {seg.end_offset_ms}UL}},")
    out += ["};", ""]
    path.write_text("\n".join(out), encoding="utf-8")


def main() -> int:
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    cfg_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    s = build(c)
    ok = checks(s)

    SEQ.mkdir(exist_ok=True)
    table_hash = write_table(SEQ / "f0b_pt2_command_table.csv", s)
    write_arrivals(SEQ / "f0b_pt2_arrivals.csv", s)
    write_header(FIRMWARE / "f0b_pt2_config.h", c, s, cfg_hash, table_hash)

    by_block: dict[str, float] = {}
    for seg in s.segs:
        by_block[seg.block] = by_block.get(seg.block, 0.0) + (seg.t1_ms - seg.t0_ms)
    rates = [seg.rate_mcps / 1e3 / counts_per_step(seg.mres)
             for seg in s.segs if seg.rate_mcps]
    summary = {
        "experiment": c["experiment"],
        "segments": len(s.segs),
        "arrivals": len(s.arrivals),
        "marks": sum(1 for q in s.segs if q.kind == MARK),
        "mres_switches": sum(1 for q in s.segs if q.kind == SETMRES),
        "duration_min": s.t / 60000.0,
        "gate_min": c["gate_min"],
        "margin_to_gate_min": c["gate_min"] - s.t / 60000.0,
        "blocks_min": {b: v / 60000.0 for b, v in by_block.items()},
        "fastest_step_rate_hz": max(rates),
        "slowest_step_rate_hz": min(rates),
        "m_speed_mm_s": c["m_speed_mm_s"],
        "config_sha256": cfg_hash,
        "table_sha256": table_hash,
        "checks": ok,
    }
    (SEQ / "f0b_pt2_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8")

    print(f"F0b_pt2: {len(s.segs)} segments, {len(s.arrivals)} arrivals, "
          f"{s.t / 60000:.2f} min (gate {c['gate_min']:.0f} min, "
          f"{summary['margin_to_gate_min']:.1f} min of margin)")
    for b, v in by_block.items():
        print(f"  {b:3s} {v / 60000:6.2f} min")
    print(f"  step rate {min(rates):.2f} Hz to {max(rates):.0f} Hz "
          f"(bench limit {c['max_step_rate_hz']} Hz)")
    print(f"  {summary['mres_switches']} MRES switches, "
          f"block M at {c['m_speed_mm_s']:g} mm/s")
    bad = [k for k, v in ok.items() if not v]
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
