#!/usr/bin/env python3
"""Build the timed F1 segment plan from f1_config.json.

Every move, settle, turning pause and reference stop is one segment with its
start time, absolute start/end position (rounded to 1/256 microsteps) and
tags. The plan drives the overview figure and gives the segment count and
duration; the firmware command table is generated once F0b fixes the
angle-bracket values (s90, d_min).

    python3 f1_sequence.py          # writes sequence/f1_segments.csv + summary
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f1_config.json"
SEQ = HERE / "sequence"


@dataclass
class Segment:
    group: str          # R1, A, B, R2, C, R3
    sub: str            # V, P, M, Q, R, or '' for transits
    kind: str           # measure, pause, setup, transit, ref
    t0: float           # s
    t1: float
    u0: float           # absolute commanded position, um
    u1: float
    ref_um: float       # test position the trajectory panel is drawn against
    tag: str = ""
    last_dir: int = 0   # direction of the last motion before this segment (+1 up, -1 down)
    ramp_type: str = "" # 'reversal' / 'continuation' for declared measurement ramps


@dataclass
class Plan:
    c: dict
    segs: list = field(default_factory=list)
    now: float = 0.0
    u: float = 0.0
    last_dir: int = 0

    def snap(self, x_um):
        q = self.c["ustep_um"]
        return round(x_um / q) * q

    def move_time(self, d_um, v_um_s):
        a = self.c["accel_mm_s2"] * 1000
        d = abs(d_um)
        if d == 0:
            return 0.0
        return 2 * math.sqrt(d / a) if d <= v_um_s ** 2 / a else d / v_um_s + v_um_s / a

    def add(self, group, sub, kind, target, duration, ref, tag="", ramp_type=""):
        target = self.snap(target)
        self.segs.append(Segment(group, sub, kind, self.now, self.now + duration,
                                 self.u, target, ref, tag, self.last_dir, ramp_type))
        self.now += duration
        if target != self.u:
            self.last_dir = 1 if target > self.u else -1
        self.u = target

    def move(self, group, sub, kind, target, v, ref, tag=""):
        self.add(group, sub, kind, target, self.move_time(target - self.u, v), ref, tag)

    def hold(self, group, sub, kind, seconds, ref, tag=""):
        self.add(group, sub, kind, self.u, seconds, ref, tag)

    def wipe(self, group, sub, target, side, ref):
        """3 mm approach from below (side=-1) or above (+1) to target, then settle."""
        c, vt = self.c, self.c["transit_mm_s"] * 1000
        start = target + side * c["wipe_mm"] * 1000
        kind = "transit" if abs(start - self.u) > 10000 else "setup"
        self.move(group, sub, kind, start, vt, ref, "to wipe start")
        self.move(group, sub, "setup", target, vt, ref, "wipe")
        # Wipe from below (side -1) ends moving up; from above (+1) ends moving down.
        assert self.last_dir == -side, (group, sub, side, self.last_dir)
        self.hold(group, sub, "setup", c["settle_s"], ref, "settle")

    def ramp(self, group, sub, target, v, ref, tag, ramp_type):
        """Measurement ramp whose direction is checked against the last motion."""
        d = 1 if self.snap(target) > self.u else -1
        expected = {"reversal": -self.last_dir, "continuation": self.last_dir}[ramp_type]
        assert d == expected, f"{group}.{sub} {tag}: {ramp_type} ramp {d:+d} after {self.last_dir:+d}"
        self.add(group, sub, "measure", target, self.move_time(target - self.u, v), ref, tag,
                 ramp_type)

    def half(self, group, sub, target, rate, ref, tag, step_mode=False):
        if step_mode:
            n = round(abs(target - self.u) / self.c["ustep_um"])
            self.add(group, sub, "measure", target, n * self.c["step_dwell_s"], ref, tag + " steps")
        else:
            self.move(group, sub, "measure", target, rate, ref, tag)
        self.hold(group, sub, "pause", self.c["turn_pause_s"], ref, "turn")


def strokes(c):
    """P stroke grid: s90 factors then anchors, symmetric in whole microsteps."""
    q, s90, dmin = c["ustep_um"], c["f0b_dependent"]["s90_um"], c["f0b_dependent"]["d_min_um"]
    raw = [s90 * f for f in c["p_stroke_factors"]] + c["p_anchor_strokes_um"]
    out = []
    for s in raw:
        half = max(1, round(s / 2 / q))
        s_eff = 2 * half * q
        assert s_eff >= max(q, 5 * dmin), s
        out.append(s_eff)
    return out


def centred_loops(p, group, sub, xp, stroke, cycles, half_s, tag):
    c = p.c
    rate = min(stroke / half_s, c["p_rate_cap_um_s"])
    step_mode = round(stroke / c["ustep_um"]) <= c["step_segment_max_usteps"]
    lo, hi = xp - stroke / 2, xp + stroke / 2
    p.half(group, sub, lo, rate, xp, f"{tag} to start", step_mode)
    for k in range(cycles):
        p.half(group, sub, hi, rate, xp, f"{tag} c{k+1} up", step_mode)
        p.half(group, sub, lo, rate, xp, f"{tag} c{k+1} down", step_mode)
    p.half(group, sub, xp, rate, xp, f"{tag} return", step_mode)


def reference_block(p, group, ref_um):
    c, vt = p.c, p.c["transit_mm_s"] * 1000
    refs = sorted((r * 1000 for r in c["reference_mm"]), key=lambda r: abs(r - p.u))
    for r in refs:
        for side in (-1, 1):
            p.move(group, "R", "ref", r + side * c["wipe_mm"] * 1000, vt, ref_um, "approach")
            p.move(group, "R", "ref", r, vt, ref_um, "from below" if side < 0 else "from above")
            p.hold(group, "R", "ref", c["reference_stop_s"], ref_um, "stop")


def position_block(p, group, xp_mm, decreasing=False, rate_check=False):
    c = p.c
    xp, s90 = xp_mm * 1000, c["f0b_dependent"]["s90_um"]
    # V: virgin curves, reversal down after a wipe from below, then up after one from above.
    slow = c["v_slow_factor_of_s90"] * s90
    for side in (-1, 1):
        p.wipe(group, "V", xp, side, xp)
        d = side                         # reverse: back towards the side the wipe came from
        p.ramp(group, "V", xp + d * slow, c["v_slow_um_s"], xp, "virgin slow", "reversal")
        p.ramp(group, "V", xp + d * c["v_ramp_um"], c["v_fast_um_s"], xp, "virgin fast",
               "continuation")
    # P: centred loops, increasing (and decreasing at B).
    p.wipe(group, "P", xp, -1, xp)
    grid = strokes(c)
    order = grid + (grid[::-1] if decreasing else [])
    for i, s in enumerate(order):
        tag = f"P s={s:.4g} {'inc' if i < len(grid) else 'dec'}"
        centred_loops(p, group, "P", xp, s, c["p_cycles"], c["p_half_cycle_s"], tag)
    # M: nested loops; each branch follows the half-cycle rule of its own loop.
    so, si = p.snap(c["m_outer_factor"] * s90), p.snap(c["m_inner_factor"] * s90)
    ro, ri = so / c["p_half_cycle_s"], si / c["p_half_cycle_s"]
    p.wipe(group, "M", xp - so / 2, -1, xp)
    for k in range(c["m_repeats"]):
        p.half(group, "M", xp + so / 2, ro, xp, f"M{k+1} outer up")
        p.half(group, "M", xp, ro, xp, f"M{k+1} outer down to x_p")
        p.half(group, "M", xp + si, ri, xp, f"M{k+1} inner up")
        p.half(group, "M", xp, ri, xp, f"M{k+1} inner down")
        p.half(group, "M", xp - so / 2, ro, xp, f"M{k+1} outer down rest")
    # Q: rate check at s90.
    if rate_check:
        s = strokes({**c, "p_stroke_factors": [1.0], "p_anchor_strokes_um": []})[0]
        p.wipe(group, "Q", xp, -1, xp)
        for T in c["q_half_cycle_s"]:
            centred_loops(p, group, "Q", xp, s, c["q_cycles"], T, f"Q {T:g} s")


def build(c, include_r2=True):
    p = Plan(c)
    pos = {k: v * 1000 for k, v in c["positions_mm"].items()}
    p.u = pos["B"]                       # session starts homed and centred at x_B
    reference_block(p, "R1", pos["A"])
    position_block(p, "A", c["positions_mm"]["A"])
    position_block(p, "B", c["positions_mm"]["B"], decreasing=True, rate_check=True)
    if include_r2:
        reference_block(p, "R2", pos["B"])
    position_block(p, "C", c["positions_mm"]["C"])
    reference_block(p, "R3", pos["C"])
    return p


def plan(c):
    """Full plan; R2 is dropped if the session would exceed the gate (procedure §4)."""
    p = build(c, include_r2=True)
    if p.now / 60 > c["gate_min"]:
        return build(c, include_r2=False), False
    return p, True


def check_ramp_types(segs):
    """Recompute each declared ramp's type from the table itself (independent of build)."""
    last, checked = 0, 0
    for s in segs:
        d = (s.u1 > s.u0) - (s.u1 < s.u0)
        if s.ramp_type:
            expected = -last if s.ramp_type == "reversal" else last
            assert d == expected, f"{s.group}.{s.sub} {s.tag}: declared {s.ramp_type}, " \
                                  f"ramp {d:+d} after last motion {last:+d}"
            checked += 1
        if d:
            last = d
    return checked


def summarise(p, r2):
    groups, subs = {}, {}
    for s in p.segs:
        groups.setdefault(s.group, [s.t0, s.t1])[1] = s.t1
        key = (s.group, s.sub)
        subs.setdefault(key, [s.t0, s.t1])[1] = s.t1
    return {
        "segments": len(p.segs),
        "duration_min": p.now / 60,
        "r2_included": r2,
        "strokes_um": strokes(p.c),
        "groups_min": {g: (b - a) / 60 for g, (a, b) in groups.items()},
        "sub_blocks_min": {f"{g}.{s}": (b - a) / 60 for (g, s), (a, b) in subs.items()},
        "time_by_kind_min": {k: sum(s.t1 - s.t0 for s in p.segs if s.kind == k) / 60
                             for k in ("measure", "pause", "setup", "transit", "ref")},
    }


def main():
    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    p, r2 = plan(c)
    checked = check_ramp_types(p.segs)
    SEQ.mkdir(exist_ok=True)
    with (SEQ / "f1_segments.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(p.segs[0])))
        w.writeheader()
        for s in p.segs:
            row = asdict(s)
            row.update({k: f"{row[k]:.4f}" for k in ("t0", "t1", "u0", "u1", "ref_um")})
            w.writerow(row)
    summary = summarise(p, r2)
    (SEQ / "f1_plan_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"F1 plan: {summary['segments']} segments, {summary['duration_min']:.1f} min "
          f"(R2 {'included' if r2 else 'dropped'}; gate {c['gate_min']:.0f} min); "
          f"{checked} declared ramps match their reversal/continuation type")
    for g, m in summary["groups_min"].items():
        print(f"  {g:3s} {m:5.1f} min")


if __name__ == "__main__":
    main()
