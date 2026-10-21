#!/usr/bin/env python3
"""ISO 230-2 session design: check the config, place the targets, write the header.

    python scripts/iso_design.py                 # ROM = rom_expected_mm from the config
    python scripts/iso_design.py --rom 122.69    # ROM as measured by homing

Reads `../config/iso_config.json` and writes:

  firmware/iso230_2/iso_config.h        the config, compiled into the sketch
  sequence/iso_targets_seed<S>_rom<R>.csv   the targets for that ROM
  sequence/iso_plan_summary.json        derived values, checks, duration

The board places its own targets after homing, from the travel it measures.
`place_targets()` here is the same algorithm - same xorshift32 generator, same
order of floating-point operations, same rounding - so for a given (ROM, seed)
it gives the same list, to the microstep. Pass the ROM the board printed in
`# HOMED` to reproduce its list offline.

Coordinates are mm from the home switch's slow trigger, increasing away from
it. All fixed-side approaches come from below, i.e. from the home side.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "iso_config.json"
HEADER = ROOT / "firmware" / "iso230_2" / "iso_config.h"
OUT = ROOT / "sequence"

U256 = 256
MAX_PULSE_RATE_HZ = 20_000        # firmware step-generator ceiling
MAX_TARGETS = 64                  # must match the sketch
TCLK_HZ = 12e6                    # TMC2209 internal clock


# --------------------------------------------------------------- derived

def upmm(c) -> float:
    """u256 per mm."""
    return c["steps_per_rev"] * U256 / c["lead_mm"]


def fs_per_mm(c) -> float:
    return c["steps_per_rev"] / c["lead_mm"]


def t_hold_s(c) -> float:
    """TPOWERDOWN delay plus the IHOLDDELAY ramp, upper bound (31 steps)."""
    unit = 2**18 / TCLK_HZ
    reduces = c["ihold_multiplier"] < 1.0
    return c["tpowerdown"] * unit + (c["iholddelay"] * 31 * unit if reduces else 0.0)


def t_settle_s(c) -> float:
    return t_hold_s(c) + c["settle_margin_s"]


# --------------------------------------------------------------- targets

class XorShift32:
    """Identical to the sketch's generator, so the offsets match bit for bit."""

    def __init__(self, seed: int):
        self.s = (seed & 0xFFFFFFFF) or 1

    def uniform(self) -> float:
        x = self.s
        x ^= (x << 13) & 0xFFFFFFFF
        x ^= x >> 17
        x ^= (x << 5) & 0xFFFFFFFF
        self.s = x
        return (x >> 8) / 16777216.0


def snap_u256(x_mm: float, c) -> int:
    """mm -> u256 on the microstep grid of the configured MRES."""
    step = U256 // c["mres"]
    u = math.floor(x_mm * upmm(c) + 0.5)
    return ((u + step // 2) // step) * step


def place_targets(rom_mm: float, c):
    """Distribute ISO targets over the measured travel. Mirrors placeTargets()."""
    soft_min = c["end_margin_mm"]
    soft_max = rom_mm - c["end_margin_mm"]
    x_start = soft_min + c["overrun_mm"]
    x_end_raw = soft_max - c["overrun_mm"]
    q = c["target_end_quantum_mm"]
    x_end = math.floor(x_end_raw / q) * q
    if x_end - x_start < c["p_nom_mm"]:
        raise ValueError(f"travel {rom_mm} mm leaves no room for two targets")
    m = int(math.floor((x_end - x_start) / c["p_nom_mm"] + 0.5)) + 1
    if m > MAX_TARGETS:
        raise ValueError(f"{m} targets exceeds the sketch's {MAX_TARGETS}")
    p = (x_end - x_start) / (m - 1)
    rng = XorShift32(c["seed"])
    mm = []
    for k in range(m):
        u = rng.uniform()
        t = x_start + k * p + (u * 2.0 - 1.0) * p / 4.0
        if k == 0:
            t = x_start
        if k == m - 1:
            t = x_end
        mm.append(math.floor(t * 1000.0 + 0.5) / 1000.0)
    u256 = [snap_u256(t, c) for t in mm]
    if c["x_ref_mm"] is None:
        x_ref_mm = (soft_min + soft_max) / 2.0
    else:
        x_ref_mm = float(c["x_ref_mm"])
    x_ref = math.floor(x_ref_mm * upmm(c) / U256 + 0.5) * U256    # on a full step
    return dict(soft_min=soft_min, soft_max=soft_max, x_start=x_start,
                x_end=x_end, m=m, p=p, targets_mm=mm, targets_u256=u256,
                x_ref_u256=x_ref)


# -------------------------------------------------------------- duration

def move_s(dist_mm: float, c) -> float:
    d = abs(dist_mm) * fs_per_mm(c)
    v = c["v_mm_s"] * fs_per_mm(c)
    a = c["a_mm_s2"] * fs_per_mm(c)
    if d == 0:
        return 0.0
    if d >= v * v / a:
        return v / a + d / v
    return 2.0 * math.sqrt(d / a)


def run_s(plan, c) -> float:
    """One ISO run (= one warm-up cycle), including its x_ref visit."""
    dwell = t_settle_s(c) + c["t_record_s"]
    to_mm = lambda u: u / upmm(c)
    xs = [to_mm(u) for u in plan["targets_u256"]]
    xr = to_mm(plan["x_ref_u256"])
    pre = c["ref_preapproach_mm"]
    ov = c["overrun_mm"]
    path = [xr - pre, xr, xs[0] - ov] + xs + [xs[-1] + ov] + xs[::-1] + [xs[0] - ov]
    t, pos = 0.0, xs[0] - ov
    for x in path:
        t += move_s(x - pos, c)
        pos = x
    return t + dwell * (1 + 2 * len(xs))


def ref_visit_s(plan, c) -> float:
    dwell = t_settle_s(c) + c["t_record_s"]
    return 2 * move_s(c["ref_preapproach_mm"], c) + dwell


# ---------------------------------------------------------------- checks

def check(c, plan):
    issues = []
    for key in ("v_mm_s", "home_seek_v_mm_s", "home_slow_v_mm_s"):
        rate = c[key] * fs_per_mm(c) * c["mres"]
        if rate > MAX_PULSE_RATE_HZ:
            issues.append(f"{key}={c[key]} mm/s needs {rate:.0f} Hz at 1/{c['mres']}, "
                          f"over the {MAX_PULSE_RATE_HZ} Hz step generator")
    if c["mres"] < 1 or c["mres"] > 256 or c["mres"] & (c["mres"] - 1):
        issues.append(f"mres {c['mres']} is not a power of two in 1..256")
    if c["home_side"] not in ("min", "max"):
        issues.append("home_side must be 'min' or 'max'")
    if c["end_motor_state"] not in ("hold", "disabled"):
        issues.append("end_motor_state must be 'hold' or 'disabled'")
    lo = plan["soft_min"]
    hi = plan["soft_max"]
    xr = plan["x_ref_u256"] / upmm(c)
    for name, x in (("first target - overrun", plan["targets_mm"][0] - c["overrun_mm"]),
                    ("last target + overrun", plan["targets_mm"][-1] + c["overrun_mm"]),
                    ("x_ref - pre-approach", xr - c["ref_preapproach_mm"]),
                    ("x_ref", xr)):
        if not lo - 1e-9 <= x <= hi + 1e-9:
            issues.append(f"{name} = {x:.3f} mm is outside the soft limits {lo}..{hi}")
    if c["home_backoff_mm"] * 2 > hi - lo:
        issues.append("home_backoff_mm is too large for this travel")
    return issues


# ---------------------------------------------------------------- header

def write_header(c, digest: str) -> None:
    def f(x):
        return f"{float(x):.6f}"
    x_ref = "-1.0" if c["x_ref_mm"] is None else f(c["x_ref_mm"])
    lines = [
        "/* GENERATED by scripts/iso_design.py from config/iso_config.json.",
        " * Do not edit - change the JSON and regenerate, then reflash.",
        f" * config sha256 {digest}",
        " */",
        "#pragma once",
        "#include <stdint.h>",
        "",
        f'constexpr char     CFG_HASH[] = "{digest[:12]}";',
        f"constexpr double   CFG_LEAD_MM = {f(c['lead_mm'])};",
        f"constexpr uint16_t CFG_STEPS_PER_REV = {int(c['steps_per_rev'])};",
        f"constexpr uint16_t CFG_MRES = {int(c['mres'])};",
        f"constexpr uint16_t CFG_RUN_CURRENT_MA = {int(c['run_current_ma'])};",
        f"constexpr float    CFG_IHOLD_MULT = {f(c['ihold_multiplier'])}F;",
        f"constexpr uint8_t  CFG_TPOWERDOWN = {int(c['tpowerdown'])};",
        f"constexpr uint8_t  CFG_IHOLDDELAY = {int(c['iholddelay'])};",
        f"constexpr bool     CFG_INTPOL = {'true' if c['intpol'] else 'false'};",
        f"constexpr double   CFG_V_MM_S = {f(c['v_mm_s'])};",
        f"constexpr double   CFG_A_MM_S2 = {f(c['a_mm_s2'])};",
        f"constexpr bool     CFG_HOME_MIN = {'true' if c['home_side'] == 'min' else 'false'};",
        f"constexpr double   CFG_HOME_SEEK_MM_S = {f(c['home_seek_v_mm_s'])};",
        f"constexpr double   CFG_HOME_SLOW_MM_S = {f(c['home_slow_v_mm_s'])};",
        f"constexpr double   CFG_HOME_BACKOFF_MM = {f(c['home_backoff_mm'])};",
        f"constexpr bool     CFG_MEASURE_FAR_END = {'true' if c['measure_far_end'] else 'false'};",
        f"constexpr double   CFG_ROM_EXPECTED_MM = {f(c['rom_expected_mm'])};",
        f"constexpr double   CFG_END_MARGIN_MM = {f(c['end_margin_mm'])};",
        f"constexpr double   CFG_OVERRUN_MM = {f(c['overrun_mm'])};",
        f"constexpr double   CFG_P_NOM_MM = {f(c['p_nom_mm'])};",
        f"constexpr double   CFG_TARGET_END_QUANTUM_MM = {f(c['target_end_quantum_mm'])};",
        f"constexpr uint32_t CFG_SEED = {int(c['seed'])}UL;",
        f"constexpr double   CFG_X_REF_MM = {x_ref};          // < 0: centre of the soft limits",
        f"constexpr double   CFG_REF_PREAPPROACH_MM = {f(c['ref_preapproach_mm'])};",
        f"constexpr uint32_t CFG_T_SETTLE_MS = {round(t_settle_s(c) * 1000)};   // t_hold + margin",
        f"constexpr uint32_t CFG_T_RECORD_MS = {round(c['t_record_s'] * 1000)};",
        f"constexpr uint8_t  CFG_N_RUNS = {int(c['n_runs'])};",
        f"constexpr uint8_t  CFG_DATUM_VISITS = {int(c['datum_visits'])};",
        f"constexpr uint8_t  CFG_END_REF_VISITS = {int(c['end_ref_visits'])};",
        f"constexpr uint8_t  CFG_WARMUP_REF_EVERY = {int(c['warmup_ref_every_cycles'])};",
        f"constexpr uint16_t CFG_WARMUP_T_MAX_MIN = {int(c['warmup_t_max_min'])};",
        f"constexpr bool     CFG_END_HOLD = {'true' if c['end_motor_state'] == 'hold' else 'false'};",
        f"constexpr bool     CFG_AUTO_RUN = {'true' if c.get('auto_run_on_boot') else 'false'};",
        f"constexpr uint16_t CFG_AUTO_COUNTDOWN_S = {int(c.get('auto_countdown_s', 10))};",
        f"constexpr uint16_t CFG_AUTO_CENTRE_WAIT_S = {int(c.get('auto_centre_wait_s', 15))};",
        "",
    ]
    HEADER.parent.mkdir(parents=True, exist_ok=True)
    HEADER.write_text("\n".join(lines), encoding="utf-8")


# ------------------------------------------------------------------ main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rom", type=float, help="travel between the switch triggers [mm]")
    args = ap.parse_args()

    raw = CONFIG.read_bytes()
    c = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    rom = args.rom if args.rom else c["rom_expected_mm"]
    plan = place_targets(rom, c)
    issues = check(c, plan)

    run = run_s(plan, c)
    ref = ref_visit_s(plan, c)
    iso_min = (c["datum_visits"] * ref + c["n_runs"] * run
               + c["end_ref_visits"] * ref) / 60.0
    stops = c["n_runs"] * 2 * plan["m"]

    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"seed{c['seed']}_rom{rom:g}"
    with (OUT / f"iso_targets_{tag}.csv").open("w", encoding="utf-8") as fh:
        fh.write("i,target_mm,cmd_u256,cmd_mm\n")
        for i, (t, u) in enumerate(zip(plan["targets_mm"], plan["targets_u256"]), 1):
            fh.write(f"{i},{t:.3f},{u},{u / upmm(c):.6f}\n")

    phase = sorted((t % c["lead_mm"]) / c["lead_mm"] for t in plan["targets_mm"])
    summary = {
        "config_sha256": digest,
        "rom_mm": rom,
        "u256_per_mm": upmm(c),
        "microstep_um": 1000.0 / (fs_per_mm(c) * c["mres"]),
        "soft_limits_mm": [plan["soft_min"], plan["soft_max"]],
        "span_mm": [plan["x_start"], plan["x_end"]],
        "targets": plan["m"],
        "interval_mm": plan["p"],
        "offset_max_mm": plan["p"] / 4,
        "x_ref_mm": plan["x_ref_u256"] / upmm(c),
        "t_hold_s": t_hold_s(c),
        "t_settle_s": t_settle_s(c),
        "dwell_s": t_settle_s(c) + c["t_record_s"],
        "run_minutes": run / 60.0,
        "iso_minutes": iso_min,
        "measured_stops": stops,
        "pulse_rate_hz_at_v": c["v_mm_s"] * fs_per_mm(c) * c["mres"],
        "lead_phase_min_gap": min(b - a for a, b in zip(phase, phase[1:] + [phase[0] + 1])),
        "lead_phase_max_gap": max(b - a for a, b in zip(phase, phase[1:] + [phase[0] + 1])),
        "checks": issues,
    }
    (OUT / "iso_plan_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                              encoding="utf-8")
    write_header(c, digest)

    print(f"ROM {rom:g} mm -> {plan['m']} targets over {plan['x_start']:g}..{plan['x_end']:g} mm, "
          f"interval {plan['p']:.4f} mm (offsets +/-{plan['p'] / 4:.3f} mm)")
    print(f"x_ref {summary['x_ref_mm']:.4f} mm, soft limits {plan['soft_min']:g}..{plan['soft_max']:g} mm")
    print(f"dwell {summary['dwell_s']:.2f} s (settle {summary['t_settle_s']:.2f} + record "
          f"{c['t_record_s']:g}); one run {run / 60:.1f} min; ISO sequence "
          f"{iso_min:.1f} min for {c['n_runs']} runs, {stops} measured stops")
    print(f"lead-phase coverage: largest gap {summary['lead_phase_max_gap']:.2f} of a "
          f"lead period over {plan['m']} targets")
    for i in issues:
        print(f"CHECK: {i}")
    print(f"wrote {HEADER.relative_to(ROOT)} and sequence/iso_targets_{tag}.csv")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
