#!/usr/bin/env python3
"""Build the F0b command table from the procedure, and time it.

The host generates the whole session as rows and the board replays them, as
in E5. Nothing here is measured: this is the plan, and its only job is to be
the single place the session's numbers come from - the command table that
goes to the firmware, the duration that has to clear the 58 min gate, and the
figure that explains the session, all read the same file.

Writes, into `sequence/`:

    f0b_command_table.csv   one row per segment, the format of procedure 2.1
    f0b_schedule.csv        the same rows replayed on the timing model
    f0b_plan_summary.json   per-block durations, counts, and the checks

    python f0b_plan.py
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEQ = HERE / "sequence"


# --------------------------------------------------------------- timing

def move_time(distance_mm: float, v: float, a: float) -> float:
    """Trapezoidal transit: triangular below the ramp distance, else flat-top."""
    d = abs(distance_mm)
    if d == 0.0:
        return 0.0
    ramp = v * v / a
    if d <= ramp:
        return 2.0 * math.sqrt(d / a)
    return 2.0 * v / a + (d - ramp) / v


class Session:
    """Replays the plan as it is built, so every row carries its own clock."""

    def __init__(self, cfg: dict):
        self.cfg = cfg
        s = cfg["session"]
        self.u_per_mm = cfg["drive"]["usteps_per_mm"]
        self.v = s["transit_velocity_mm_s"]
        self.a = s["transit_acceleration_mm_s2"]
        self.t = 0.0
        self.x = s["x0_mm"]              # recording starts after centring
        self.block = ""
        self.rows: list[dict] = []
        self.sched: list[dict] = []
        self.rates: list[float] = []

    # -- bookkeeping ----------------------------------------------------
    def _emit(self, seg: str, x_end: float, duration: float, tag: str,
              rate_hz: float | None = None, dwell_s: float = 0.0) -> None:
        self.rows.append({
            "row": len(self.rows) + 1,
            "block": self.block,
            "seg": seg,
            "u_target": round(x_end * self.u_per_mm),
            "rate_hz": "" if rate_hz is None else f"{rate_hz:.4f}",
            "dwell_s": f"{dwell_s:g}",
            "tag": tag,
        })
        self.sched.append({
            "row": len(self.rows),
            "block": self.block,
            "seg": seg,
            "t_start_s": f"{self.t:.4f}",
            "t_end_s": f"{self.t + duration:.4f}",
            "x_start_mm": f"{self.x:.9f}",
            "x_end_mm": f"{x_end:.9f}",
            "tag": tag,
        })
        self.t += duration
        self.x = x_end

    # -- segment types --------------------------------------------------
    def move(self, x: float, tag: str) -> None:
        """Trapezoidal transit. Not analysed; it only has to cost its time."""
        self._emit("move", x, move_time(x - self.x, self.v, self.a), tag)

    def ramp(self, x: float, rate_hz: float, tag: str) -> None:
        """Constant step rate to a target - the segments that are analysed."""
        usteps = abs(x - self.x) * self.u_per_mm
        self.rates.append(rate_hz)
        self._emit("ramp", x, usteps / rate_hz, tag, rate_hz=rate_hz)

    def steps(self, n: int, direction: int, dwell: float, tag: str) -> None:
        """n single microsteps, each followed by `dwell`."""
        x = self.x + direction * n / self.u_per_mm
        self._emit("step", x, n * dwell, tag, dwell_s=dwell)

    def hold(self, seconds: float, tag: str) -> None:
        self._emit("dwell", self.x, seconds, tag, dwell_s=seconds)

    # -- composites -----------------------------------------------------
    def wipe(self, x: float, from_below: bool, tag: str) -> None:
        """Procedure 2: a 3 mm approach in the stated direction, then settle.

        The run-up is what wipes the memory, so it is a move in its own
        right and is timed as one.
        """
        s = self.cfg["session"]
        sign = -1.0 if from_below else 1.0
        self.move(x + sign * s["wipe_mm"], f"{tag}:run-up")
        self.move(x, f"{tag}:approach-{'up' if from_below else 'down'}")
        self.hold(s["wipe_settle_s"], f"{tag}:settle")


# ----------------------------------------------------------------- blocks

def reference_block(S: Session, name: str) -> None:
    """15 and 105 mm, each from both directions, in the procedure's order."""
    r, s = S.cfg["reference"], S.cfg["session"]
    S.block = name
    for x in r["positions_mm"]:
        for direction, word in ((1, "up"), (-1, "down")):
            S.move(x - direction * r["overrun_mm"], f"{name}:run-up:{x:g}")
            S.move(x, f"{name}:stop:{x:g}:{word}")
            S.hold(s["stop_settle_s"] + s["stop_record_s"],
                   f"{name}:settle+record:{x:g}:{word}")
    S.move(s["x0_mm"], f"{name}:return")


def block_n(S: Session) -> None:
    """Noise floor: two long holds, at x0 and 20 mm above it."""
    n, x0 = S.cfg["N"], S.cfg["session"]["x0_mm"]
    S.block = "N"
    S.wipe(x0, True, "N1")
    S.hold(n["n1_hold_s"], "N1:hold")
    S.wipe(x0 + n["n2_offset_mm"], True, "N2")
    S.hold(n["n2_hold_s"], "N2:hold")
    S.move(x0, "N:return")


def block_l(S: Session) -> None:
    """Coarse reversal ladder: anchored loops above x0, increasing then back."""
    cfg, x0 = S.cfg["L"], S.cfg["session"]["x0_mm"]
    pause = S.cfg["session"]["turn_pause_s"]
    u_per_mm = S.u_per_mm
    cap_hz = cfg["rate_cap_mm_s"] * u_per_mm
    S.block = "L"
    S.wipe(x0, True, "L")
    order = {"increasing": cfg["strokes_usteps"],
             "decreasing": list(reversed(cfg["strokes_usteps"]))}
    for pass_name in cfg["passes"]:
        for s_u in order[pass_name]:
            s_mm = s_u / u_per_mm
            for half in range(4):
                up = half % 2 == 0
                target = x0 + s_mm if up else x0
                # half 0 is a continuation after the wipe; 1-3 are reversals
                kind = "continuation" if half == 0 else "reversal"
                tag = f"L:{pass_name}:s{s_u}:half{half + 1}:{kind}"
                if s_u <= cfg["step_mode_max_usteps"]:
                    S.steps(s_u, 1 if up else -1, cfg["step_dwell_s"], tag)
                else:
                    S.ramp(target, min(s_u / cfg["half_cycle_s"], cap_hz), tag)
                S.hold(pause, f"{tag}:turn")


def block_s(S: Session) -> None:
    """Slow ramps: 20 um at 0.25 um/s, reversal and continuation, both ways."""
    cfg, x0 = S.cfg["S"], S.cfg["session"]["x0_mm"]
    rate = cfg["ramp_rate_um_s"] * S.u_per_mm / 1000.0
    span = cfg["ramp_um"] / 1000.0
    S.block = "S"
    for rep in range(1, cfg["repeats"] + 1):
        for r in cfg["ramps"]:
            tag = f"{r['name']}:rep{rep}:{r['kind']}"
            S.wipe(x0, r["wipe_from_below"], tag)
            S.ramp(x0 + r["direction"] * span, rate, f"{tag}:ramp")
            S.move(x0, f"{tag}:return")


def block_d(S: Session) -> None:
    """Dwell probe: 1 s and 60 s at rest, ABBA, each followed by a reversal."""
    cfg, x0 = S.cfg["D"], S.cfg["session"]["x0_mm"]
    rate = cfg["probe_rate_um_s"] * S.u_per_mm / 1000.0
    span = cfg["probe_um"] / 1000.0
    S.block = "D"
    for i, td in enumerate(cfg["dwells_s"], start=1):
        tag = f"D{i}:Td{td:g}s"
        S.wipe(x0, True, tag)
        S.hold(td, f"{tag}:dwell")
        S.ramp(x0 - span, rate, f"{tag}:reverse-probe")
        S.move(x0, f"{tag}:return")


def block_v(S: Session) -> None:
    """Speed probe: the same 1 mm segment fast and slow, ABA."""
    cfg, x0 = S.cfg["V"], S.cfg["session"]["x0_mm"]
    pause = S.cfg["session"]["turn_pause_s"]
    S.block = "V"
    for i, kind in enumerate(cfg["order"], start=1):
        speed = cfg["fast_mm_s"] if kind == "fast" else cfg["slow_mm_s"]
        rate = speed * S.u_per_mm
        S.ramp(x0 + cfg["segment_mm"], rate, f"V{i}:{kind}:up")
        S.hold(pause, f"V{i}:{kind}:turn")
        S.ramp(x0, rate, f"V{i}:{kind}:down")
        S.hold(pause, f"V{i}:{kind}:turn")


# ----------------------------------------------------------------- checks

def build(cfg: dict) -> Session:
    S = Session(cfg)
    for name in cfg["session"]["order"]:
        if name.startswith("R"):
            reference_block(S, name)
        else:
            {"N": block_n, "L": block_l, "S": block_s,
             "D": block_d, "V": block_v}[name](S)
    return S


def block_times(S: Session) -> dict:
    out: dict[str, float] = {}
    for row in S.sched:
        span = float(row["t_end_s"]) - float(row["t_start_s"])
        out[row["block"]] = out.get(row["block"], 0.0) + span
    return out


def checks(S: Session, cfg: dict) -> dict:
    xs = [float(r["x_start_mm"]) for r in S.sched]
    xs += [float(r["x_end_mm"]) for r in S.sched]
    lo, hi = cfg["session"]["soft_limits_mm"]
    strokes = cfg["L"]["strokes_usteps"]
    grid = cfg["drive"]["electrical_cycle_usteps"]
    ustep_um = cfg["drive"]["ustep_um"]
    floor = cfg["reference_values"]["interferometer_floor_um"]
    beta0 = cfg["reference_values"]["beta0_um"]
    minutes = S.t / 60.0
    return {
        "x0_on_electrical_cycle_grid": cfg["session"]["u0"] % grid == 0,
        "u0_matches_x0": (cfg["session"]["u0"]
                          == round(cfg["session"]["x0_mm"] * S.u_per_mm)),
        "all_motion_inside_soft_limits": lo <= min(xs) and max(xs) <= hi,
        "no_rate_above_driver_limit": (
            max(S.rates) <= cfg["drive"]["max_step_rate_hz"] + 1e-6),
        "ladder_rate_within_cap": all(
            float(r["rate_hz"]) <= cfg["L"]["rate_cap_mm_s"] * S.u_per_mm + 1e-6
            for r in S.rows if r["block"] == "L" and r["rate_hz"]),
        "strokes_match_procedure_table": strokes == sorted(strokes),
        "smallest_stroke_above_interferometer_floor": ustep_um > 3 * floor,
        "ladder_brackets_beta0": (strokes[0] * ustep_um < beta0
                                  < strokes[-1] * ustep_um),
        "full_plan_below_58_min_gate": minutes < cfg["session"]["watchdog_gate_min"],
        "full_plan_inside_scope_record": minutes < cfg["session"]["scope_record_min"],
    }


def main() -> int:
    cfg = json.loads((HERE / "f0b_config.json").read_text(encoding="utf-8"))
    S = build(cfg)
    SEQ.mkdir(parents=True, exist_ok=True)

    with (SEQ / "f0b_command_table.csv").open("w", newline="",
                                              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(S.rows[0]))
        w.writeheader()
        w.writerows(S.rows)
    with (SEQ / "f0b_schedule.csv").open("w", newline="",
                                         encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(S.sched[0]))
        w.writeheader()
        w.writerows(S.sched)

    times = block_times(S)
    ladder = [r for r in S.rows if r["block"] == "L"
              and r["seg"] in ("ramp", "step")]
    summary = {
        "total_s": round(S.t, 2),
        "total_min": round(S.t / 60.0, 2),
        "block_minutes": {k: round(v / 60.0, 2) for k, v in times.items()},
        "segments": len(S.rows),
        "counts": {
            "ladder_half_cycles": len(ladder),
            "ladder_strokes": len(cfg["L"]["strokes_usteps"]),
            "slow_ramps": sum(1 for r in S.rows
                              if r["block"] == "S" and r["tag"].endswith(":ramp")),
            "dwell_trials": len(cfg["D"]["dwells_s"]),
            "speed_passes": 2 * len(cfg["V"]["order"]),
            "reference_stops": sum(1 for r in S.rows
                                   if ":stop:" in r["tag"]),
        },
        "amplitudes_um": {
            "ustep": cfg["drive"]["ustep_um"],
            "ladder_min": cfg["L"]["strokes_usteps"][0] * cfg["drive"]["ustep_um"],
            "ladder_max": cfg["L"]["strokes_usteps"][-1] * cfg["drive"]["ustep_um"],
            "slow_ramp": cfg["S"]["ramp_um"],
            "dwell_probe": cfg["D"]["probe_um"],
            "speed_segment": cfg["V"]["segment_mm"] * 1000.0,
        },
        "margin_to_gate_min": round(cfg["session"]["watchdog_gate_min"]
                                    - S.t / 60.0, 2),
        "checks": checks(S, cfg),
    }
    (SEQ / "f0b_plan_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"  {len(S.rows)} segments, {S.t / 60.0:.1f} min, "
          f"{summary['margin_to_gate_min']:.1f} min inside the gate")
    for name in cfg["session"]["order"]:
        print(f"    {name:<3} {times[name] / 60.0:6.2f} min")
    failed = [k for k, v in summary["checks"].items() if not v]
    print("  checks: all pass" if not failed
          else "  CHECKS FAILED: " + ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
