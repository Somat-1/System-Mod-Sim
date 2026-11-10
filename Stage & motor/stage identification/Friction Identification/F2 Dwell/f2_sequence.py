#!/usr/bin/env python3
"""Build the timed F2 segment plan from f2_config.json.

Every transit, run-up, wipe, settle, dwell, post-dwell ramp, hold and
reference stop is one segment with its start time, absolute start/end
position (rounded to 1/256 microsteps) and tags. The plan drives the
overview figure and gives the segment count, the trial order and the
duration.

    python3 f2_sequence.py                # procedure 4.1, the full session
    python3 f2_sequence.py --variant reduced   # procedure 4.3

writes

    sequence/f2_segments.csv          the plan, one row per segment
    sequence/f2_plan_summary.json     counts, durations, coverage, trial order

This is the first of two independent builders: f2_make_table.py builds the
table the board replays from the same configuration, by its own route, and
cross-checks itself against the summary written here. Two separately
written builders agreeing is the check that the grid in f2_config.json was
read the same way twice -- including the shuffled trial order, which each
builder derives from the logged seed with its own copy of the generator.

THE POST-DWELL RAMP IS PROVISIONAL. Its length is 2*s90, and s90 comes from
F0b. Until F0b is processed, f2_config.json carries the procedure's
placeholder defaults (s90 = 3 um, d_min = 5 nm). F2 is also conditional:
the full version is justified only if F0b's dwell probe exceeded 3*d_min.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "f2_config.json"
SEQ = HERE / "sequence"

MASK64 = (1 << 64) - 1


# -- deterministic shuffle -------------------------------------------------
# SplitMix64 feeding a downward Fisher-Yates. Written out rather than taken
# from the random module so that the shuffled order is reproducible from the
# logged seed alone, on any interpreter, and can be re-derived independently
# by f2_make_table.py and by the analysis.
def split_mix_64(seed: int):
    state = seed & MASK64

    def nxt() -> int:
        nonlocal state
        state = (state + 0x9E3779B97F4A7C15) & MASK64
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        return (z ^ (z >> 31)) & MASK64

    return nxt


def shuffled(items: list, seed: int) -> list:
    nxt = split_mix_64(seed)
    out = list(items)
    for i in range(len(out) - 1, 0, -1):
        j = nxt() % (i + 1)
        out[i], out[j] = out[j], out[i]
    return out


def block_seed(seed: int, block_index: int, rep: int) -> int:
    """One sub-seed per (block, repetition), so order within a repetition is
    independent of the other repetitions and of the block order."""
    return (seed + 1000003 * block_index + 101 * rep) & MASK64


# -- plan -----------------------------------------------------------------
@dataclass
class Segment:
    block: str          # R1, B, A, L, R3
    sub: str            # R (reference), R/C (trial mode), setup
    kind: str           # transit, setup, settle, dwell, measure, hold, ref
    trial: int          # trial number in execution order, 0 outside a trial
    rep: int            # repetition inside the block, 0 outside a trial
    td_s: float         # dwell time of the trial this segment belongs to
    t0: float           # s
    t1: float
    u0: float           # absolute commanded position, um
    u1: float
    ref_um: float       # test position the trajectory panel is drawn against
    tag: str = ""
    last_dir: int = 0   # direction of the last motion before this segment
    ramp_type: str = "" # 'reversal' / 'continuation' for the post-dwell ramp


@dataclass
class Plan:
    c: dict
    segs: list = field(default_factory=list)
    now: float = 0.0
    u: float = 0.0
    last_dir: int = 0
    block: str = "R1"
    sub: str = "setup"
    trial: int = 0
    rep: int = 0
    td: float = 0.0

    def snap(self, x_um: float) -> float:
        q = self.c["ustep_um"]
        return round(x_um / q) * q

    def move_time(self, d_um: float, v_um_s: float) -> float:
        a = self.c["accel_mm_s2"] * 1000.0
        d = abs(d_um)
        if d == 0.0 or v_um_s <= 0.0:
            return 0.0
        return 2 * math.sqrt(d / a) if d <= v_um_s ** 2 / a else d / v_um_s + v_um_s / a

    def add(self, kind, target, duration, ref, tag="", ramp_type=""):
        target = self.snap(target)
        self.segs.append(Segment(self.block, self.sub, kind, self.trial, self.rep,
                                 self.td, self.now, self.now + duration, self.u,
                                 target, ref, tag, self.last_dir, ramp_type))
        self.now += duration
        if target != self.u:
            self.last_dir = 1 if target > self.u else -1
        self.u = target

    def move(self, kind, target, v, ref, tag=""):
        self.add(kind, target, self.move_time(target - self.u, v), ref, tag)

    def hold(self, kind, seconds, ref, tag=""):
        self.add(kind, self.u, seconds, ref, tag)

    def wipe(self, target, ref, tag="wipe"):
        """3 mm approach from below at the transit speed. Every F2 trial uses
        this one preparation, so the state before the dwell is always the
        same: the last motion was upward."""
        c, vt = self.c, self.c["transit_mm_s"] * 1000.0
        start = target - c["wipe_mm"] * 1000.0
        far = abs(start - self.u) > 10000.0
        self.move("transit" if far else "setup", start, vt, ref,
                  "transit to wipe start" if far else "to wipe start")
        self.move("setup", target, vt, ref, tag + " from below")
        assert self.last_dir == 1, (self.block, self.trial, "wipe must end moving up")

    def ramp(self, target, v, ref, tag, ramp_type):
        """The post-dwell ramp. Its direction is checked against the last
        motion: mode R reverses the wipe, mode C continues it. The dwell in
        between is motionless, so last_dir still describes the wipe."""
        d = 1 if self.snap(target) > self.u else -1
        expected = {"reversal": -self.last_dir, "continuation": self.last_dir}[ramp_type]
        assert d == expected, (f"trial {self.trial} {tag}: {ramp_type} ramp "
                               f"{d:+d} after {self.last_dir:+d}")
        self.add("measure", target, self.move_time(target - self.u, v), ref, tag,
                 ramp_type)


def reference_block(p: Plan, name: str, ref_um: float) -> None:
    """15 / 105 mm, both directions, nearer reference first (as in F1)."""
    c, vt = p.c, p.c["transit_mm_s"] * 1000.0
    p.block, p.sub, p.trial, p.rep, p.td = name, "R", 0, 0, 0.0
    refs = sorted((r * 1000.0 for r in c["reference_mm"]), key=lambda r: abs(r - p.u))
    for r in refs:
        for side in (-1, 1):
            p.move("ref", r + side * c["wipe_mm"] * 1000.0, vt, ref_um, "approach")
            p.move("ref", r, vt, ref_um, "from below" if side < 0 else "from above")
            p.hold("ref", c["reference_settle_s"], ref_um, "stop settle")
            p.hold("ref", c["reference_record_s"], ref_um, "stop record")


def ramp_um(c: dict) -> float:
    """Post-dwell ramp length, 2*s90 rounded to a whole number of microsteps."""
    raw = c["ramp_factor_of_s90"] * c["f0b_dependent"]["s90_um"]
    n = max(1, round(raw / c["ustep_um"]))
    out = n * c["ustep_um"]
    assert out >= 5.0 * c["f0b_dependent"]["d_min_um"], raw
    return out


def trial_order(c: dict, blocks: list) -> list:
    """(block, rep, T_d, mode) in execution order, shuffled per repetition."""
    out = []
    for index, b in enumerate(blocks):
        reps = max(d["reps"] for d in b["dwells"])
        for rep in range(1, reps + 1):
            batch = [(d["s"], m) for d in b["dwells"] if rep <= d["reps"]
                     for m in b["modes"]]
            if b.get("shuffle", True):
                batch = shuffled(batch, block_seed(c["shuffle_seed"], index, rep))
            out += [(b["name"], rep, td, mode) for td, mode in batch]
    return out


def trial(p: Plan, xp_um: float, td: float, mode: str, number: int, rep: int) -> None:
    """One F2 trial: wipe from below, settle, dwell T_d, post-dwell ramp, hold."""
    c = p.c
    p.sub, p.trial, p.rep, p.td = mode, number, rep, td
    span = ramp_um(c)
    p.wipe(xp_um, xp_um, f"t{number} Td={td:g}s {mode}")
    p.hold("settle", c["pre_dwell_settle_s"], xp_um,
           f"t{number} settle (t=0 of the dwell at its end)")
    p.hold("dwell", td, xp_um, f"t{number} dwell Td={td:g}s")
    direction = -1 if mode == "R" else 1
    p.ramp(xp_um + direction * span, c["ramp_rate_um_s"], xp_um,
           f"t{number} post-dwell ramp {'down' if direction < 0 else 'up'}",
           "reversal" if mode == "R" else "continuation")
    p.hold("hold", c["post_ramp_hold_s"], xp_um, f"t{number} hold")


def build(c: dict, blocks: list, references: bool) -> Plan:
    p = Plan(c)
    pos = {k: v * 1000.0 for k, v in c["positions_mm"].items()}
    p.u = pos["B"]                       # the session starts homed and parked at x_B
    if references:
        reference_block(p, "R1", pos[blocks[0]["position"]])
    by_name = {b["name"]: b for b in blocks}
    number = 0
    for name, rep, td, mode in trial_order(c, blocks):
        b = by_name[name]
        p.block = name
        number += 1
        trial(p, pos[b["position"]], td, mode, number, rep)
    if references:
        reference_block(p, "R3", pos[blocks[-1]["position"]])
    return p


def blocks_of(c: dict, variant: str) -> tuple[list, bool]:
    if variant == "reduced":
        return c["reduced_blocks"], c["reduced_includes_references"]
    return c["blocks"], True


# -- checks ---------------------------------------------------------------
def check_ramp_types(segs: list) -> int:
    """Recompute every declared ramp type from the plan itself, independent
    of how it was built: mode R must reverse the last motion, mode C must
    continue it, and the dwell in between must not move."""
    last, checked = 0, 0
    for s in segs:
        d = (s.u1 > s.u0) - (s.u1 < s.u0)
        if s.kind == "dwell":
            assert d == 0, f"trial {s.trial}: dwell is not motionless"
        if s.ramp_type:
            expected = -last if s.ramp_type == "reversal" else last
            assert d == expected, (f"trial {s.trial} {s.tag}: declared "
                                   f"{s.ramp_type}, ramp {d:+d} after {last:+d}")
            checked += 1
        if d:
            last = d
    return checked


def check_wipes(segs: list) -> int:
    """Every trial's last motion before its dwell must be the upward wipe."""
    pre, checked = {}, 0
    for s in segs:
        if not s.trial:
            continue
        d = (s.u1 > s.u0) - (s.u1 < s.u0)
        if s.kind == "dwell":
            assert pre.get(s.trial) == 1, f"trial {s.trial}: dwell not after a wipe up"
            checked += 1
        if d:
            pre[s.trial] = d
    return checked


def coverage(order: list) -> dict:
    out: dict[str, int] = {}
    for name, _rep, td, mode in order:
        out[f"{name}.{td:g}s.{mode}"] = out.get(f"{name}.{td:g}s.{mode}", 0) + 1
    return out


def planned_coverage(blocks: list) -> dict:
    out: dict[str, int] = {}
    for b in blocks:
        for d in b["dwells"]:
            for m in b["modes"]:
                out[f"{b['name']}.{d['s']:g}s.{m}"] = d["reps"]
    return out


def correlation(xs: list, ys: list) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return 0.0 if sxx == 0 or syy == 0 else sxy / math.sqrt(sxx * syy)


def dwell_time_correlation(p: Plan) -> dict:
    """Randomisation check: log10(T_d) against the trial's start time, per
    block. The point of shuffling is that this stays near zero, so a trend
    in the data cannot be read as a dwell effect."""
    starts: dict[int, tuple] = {}
    for s in p.segs:
        if s.trial and s.trial not in starts:
            starts[s.trial] = (s.block, s.t0, s.td_s)
    out = {}
    for block in dict.fromkeys(b for b, _t, _d in starts.values()):
        rows = [(t, td) for b, t, td in starts.values() if b == block]
        if len(rows) >= 3:
            out[block] = correlation([math.log10(td) for _t, td in rows],
                                     [t for t, _td in rows])
    return out


def checks(p: Plan, blocks: list, order: list) -> dict:
    c = p.c
    q = c["ustep_um"]
    span = ramp_um(c)
    lo_mm, hi_mm = c["constrained_span_mm"]
    positions = [c["positions_mm"][b["position"]] for b in blocks]
    longest = max(td for _b, _r, td, _m in order)
    corr = dwell_time_correlation(p)
    return {
        "positions_on_electrical_grid": all(
            round(x * 1000.0 / q) % round(c["electrical_cycle_um"] / q) == 0
            for x in positions),
        "positions_inside_constrained_span": all(lo_mm <= x <= hi_mm for x in positions),
        "ramp_whole_usteps": abs(span / q - round(span / q)) < 1e-9,
        "ramp_above_5_dmin": span >= 5.0 * c["f0b_dependent"]["d_min_um"],
        "dwells_motionless": all(s.u0 == s.u1 for s in p.segs if s.kind == "dwell"),
        "coverage_matches_grid": coverage(order) == planned_coverage(blocks),
        "baseline_dwell_present": any(td == c["baseline_dwell_s"]
                                      for _b, _r, td, _m in order),
        "dwell_time_decorrelated": all(abs(v) <= c["max_dwell_time_correlation"]
                                       for v in corr.values()),
        "inside_gate": p.now / 60.0 <= c["gate_min"],
        "gate_plus_longest_dwell_inside_record": (
            c["gate_min"] + longest / 60.0 <= c["scope_record_min"]),
        "inside_scope_record": p.now / 60.0 <= c["scope_record_min"],
    }


# -- summary --------------------------------------------------------------
def summarise(p: Plan, variant: str, blocks: list, order: list, ok: dict) -> dict:
    by_block, by_kind = {}, {}
    for s in p.segs:
        by_block[s.block] = by_block.get(s.block, 0.0) + (s.t1 - s.t0)
        by_kind[s.kind] = by_kind.get(s.kind, 0.0) + (s.t1 - s.t0)
    trials = {s.trial for s in p.segs if s.trial}
    dwell_s = sum(s.t1 - s.t0 for s in p.segs if s.kind == "dwell")
    overhead = sum(s.t1 - s.t0 for s in p.segs
                   if s.trial and s.kind != "dwell")
    return {
        "variant": variant,
        "segments": len(p.segs),
        "trials": len(trials),
        "duration_min": p.now / 60.0,
        "gate_min": p.c["gate_min"],
        "margin_to_gate_min": p.c["gate_min"] - p.now / 60.0,
        "provisional": True,
        "s90_um": p.c["f0b_dependent"]["s90_um"],
        "d_min_um": p.c["f0b_dependent"]["d_min_um"],
        "ramp_um": ramp_um(p.c),
        "ramp_usteps": round(ramp_um(p.c) / p.c["ustep_um"]),
        "shuffle_seed": p.c["shuffle_seed"],
        "dwelling_min": dwell_s / 60.0,
        "trial_overhead_s": overhead / len(trials) if trials else 0.0,
        "longest_dwell_s": max(td for _b, _r, td, _m in order),
        "blocks_min": {b: v / 60.0 for b, v in by_block.items()},
        "time_by_kind_min": {k: v / 60.0 for k, v in sorted(by_kind.items())},
        "coverage": coverage(order),
        "dwell_time_correlation": dwell_time_correlation(p),
        "trial_order": [[b, r, td, m] for b, r, td, m in order],
        "checks": ok,
    }


def write_segments(path: Path, segs: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(segs[0])))
        w.writeheader()
        for s in segs:
            row = asdict(s)
            row.update({k: f"{row[k]:.4f}" for k in
                        ("t0", "t1", "u0", "u1", "ref_um", "td_s")})
            w.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variant", choices=("full", "reduced"), default=None,
                    help="procedure 4.1 (full) or 4.3 (reduced); "
                         "default: the variant in f2_config.json")
    args = ap.parse_args()

    c = json.loads(CONFIG.read_text(encoding="utf-8"))
    variant = args.variant or c["variant"]
    blocks, references = blocks_of(c, variant)
    order = trial_order(c, blocks)
    p = build(c, blocks, references)

    ramps = check_ramp_types(p.segs)
    wipes = check_wipes(p.segs)
    ok = checks(p, blocks, order)
    summary = summarise(p, variant, blocks, order, ok)
    summary["ramps_checked"] = ramps
    summary["wipes_checked"] = wipes

    SEQ.mkdir(exist_ok=True)
    write_segments(SEQ / "f2_segments.csv", p.segs)
    (SEQ / "f2_plan_summary.json").write_text(json.dumps(summary, indent=1) + "\n",
                                              encoding="utf-8")

    print(f"F2 plan ({variant}): {summary['segments']} segments, "
          f"{summary['trials']} trials, {summary['duration_min']:.2f} min "
          f"(gate {c['gate_min']:.0f} min, "
          f"{summary['margin_to_gate_min']:.1f} min of margin)")
    for b, m in summary["blocks_min"].items():
        print(f"  {b:3s} {m:6.2f} min")
    print(f"  dwelling {summary['dwelling_min']:.2f} min, "
          f"{summary['trial_overhead_s']:.1f} s of overhead per trial, "
          f"longest dwell {summary['longest_dwell_s']:g} s")
    print(f"  post-dwell ramp {summary['ramp_um']:.4f} um "
          f"({summary['ramp_usteps']} usteps) at {c['ramp_rate_um_s']:g} um/s, "
          f"PROVISIONAL from s90={summary['s90_um']:g} um")
    print(f"  seed {c['shuffle_seed']}, log10(T_d) vs session time: "
          + ", ".join(f"{b} r={v:+.3f}"
                      for b, v in summary["dwell_time_correlation"].items()))
    print(f"  {ramps} post-dwell ramps match their declared type, "
          f"{wipes} dwells follow a wipe from below")
    bad = [k for k, v in ok.items() if not v]
    print("checks: " + ("all pass" if not bad else "FAILED " + ", ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
