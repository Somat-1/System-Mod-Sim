"""v7 settling-error campaign: the commanded sequence, as data.

This module builds the proposed stepping/displacement sequence and nothing
else -- no plotting, no hardware. `render_sequence.py` imports it to draw it,
and a firmware author can import it (or read the exported CSV) to emit the
move list.

The experiment's question: **for a commanded position, where does the stage
actually come to rest?** Everything below is in service of separating that
error from the things that would otherwise masquerade as it -- backlash,
thermal drift, and the motor's own cyclic error.

Every number here is grounded in this repository's own measurements rather
than in generic guidance; the constants carry their source.

## Phases

| Phase | What it isolates |
|---|---|
| Warm-up | thermal transient -- excluded from measurement entirely |
| Datum | a single, direction-consistent origin |
| Unidirectional sweep | absolute accuracy, free of backlash by construction |
| Bidirectional pass | reversal error (backlash), by deliberate contrast |
| Cyclic comb | the 4-full-step electrical-cycle error, sampled at 1/4 step |
| Settling repeats | repeatability at one target, same move every time |
| Drift witness | the thermal drift that would otherwise be read as error |

## Where the numbers come from

* **1 full step = 10 um** -- `FULL_STEPS_PER_MM = 100` in the v6 firmware.
* **Electrical cycle = 4 full steps = 40 um** -- the cyclic error's period.
* **Travel = 44 mm** -- measured in `../../v4/scripts/tmp_v4_manual_range_sweep`.
* **Move rate 70 full-steps/s = 700 um/s** -- v6's MODERATE. One rate is used
  for *every* measured move: changing the rate changes the rotor dynamics and
  so changes the settling, which would confound the thing being measured.
* **Settling discard 0.20 s** -- v2.3 measured the free decay at 343.8-350.0 Hz
  with tau = 15-21 ms. 200 ms is ~10 time constants, so the ring-down is gone
  by a wide margin.
* **Dwell 1.5 s** -- long enough to average 1300 samples after the discard,
  short enough that v6.1's slow wander (+0.49 nm/s, and larger excursions over
  tens of seconds) does not accumulate inside one measurement.
* **Back-off 500 um** -- v6 measured reversal errors up to 239 um at MRES 1.
  The back-off before any re-approach must clear the worst backlash with
  margin, or the "unidirectional" approach is not actually unidirectional.
* **Targets snapped to whole full steps** -- the motor is accurate at full-step
  positions; a target between them confounds placement error with microstep
  non-linearity.
* **Targets jittered off the 2 mm lead pitch** -- a uniform 2 mm grid would
  alias exactly with the lead screw, and a once-per-revolution error would
  then look like a constant offset instead of a periodic one.
"""
import json
from pathlib import Path

import numpy as np

V7_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = V7_DIR / "sequence"

# --- Rig constants (see module docstring for provenance) --------------------
FULL_STEP_UM = 10.0
ELECTRICAL_CYCLE_STEPS = 4
ELECTRICAL_CYCLE_UM = ELECTRICAL_CYCLE_STEPS * FULL_STEP_UM   # 40 um
LEAD_PITCH_UM = 2000.0
TRAVEL_TOTAL_UM = 44_000.0
GUARD_UM = 2_000.0

# --- Motion ----------------------------------------------------------------
RATE_FULLSTEPS_S = {"SLOW": 27.5, "MODERATE": 70.0, "FAST": 200.0}
MEASURE_RATE = "MODERATE"            # the one rate every measured move uses
TRAVERSE_RATE = "FAST"               # unmeasured repositioning only
WARMUP_RATE = "MODERATE"

# --- Timing ----------------------------------------------------------------
DWELL_MEASURE_S = 1.5
SETTLE_DISCARD_S = 0.20
DATUM_DWELL_S = 3.0
WARMUP_MIN = 20.0

# --- Geometry --------------------------------------------------------------
SWEEP_HALF_UM = 18_000.0
N_SWEEP_TARGETS = 19
TARGET_JITTER_UM = 700.0
JITTER_SEED = 7
BACKOFF_UM = 500.0

COMB_CENTRES_UM = (-12_000.0, 0.0, 12_000.0)
COMB_STEP_UM = FULL_STEP_UM / 4.0            # 2.5 um = 1/4 full step
COMB_SPAN_CYCLES = 3                          # 120 um per comb
WITNESS_EVERY = 5
N_SETTLE_REPEATS = 20
SETTLE_TARGET_UM = 0.0
SETTLE_APPROACH_UM = 2_000.0

MRES_VALUES = (1, 4, 16, 32)

PHASES = ("warmup", "datum", "unidirectional", "bidirectional", "cyclic",
          "settling", "witness")


def snap_full_step(x_um):
    """Nearest whole full step. Targets must land on full-step boundaries."""
    return np.round(np.asarray(x_um, dtype=float) / FULL_STEP_UM) * FULL_STEP_UM


def sweep_targets():
    """Global targets: full-step aligned, deliberately off the lead pitch.

    A plain linspace over this span lands on an exact 2000 um spacing -- the
    lead screw's pitch -- so every target would sample the same screw phase and
    a once-per-revolution error would read as a constant offset. A bounded
    pseudo-random jitter (fixed seed, so the plan is reproducible) spreads the
    sampled phase across the revolution instead, then each target is snapped
    back onto a whole full step.
    """
    base = np.linspace(-SWEEP_HALF_UM, SWEEP_HALF_UM, N_SWEEP_TARGETS)
    rng = np.random.default_rng(JITTER_SEED)
    jitter = rng.uniform(-TARGET_JITTER_UM, TARGET_JITTER_UM, base.size)
    jitter[0] = jitter[-1] = 0.0                  # keep the span's ends exact
    t = snap_full_step(base + jitter)
    if np.any(np.diff(t) <= 0):
        raise RuntimeError("jitter reordered the targets; reduce it")
    if np.max(np.abs(t)) > SWEEP_HALF_UM + TARGET_JITTER_UM:
        raise RuntimeError("a target escaped the sweep span")
    return t


def comb_step_for(mres):
    """The comb's spacing at a given MRES.

    The target is 1/4 of a full step, but a driver cannot be asked for less
    than one microstep: at MRES 1 that *is* a full step, so a 2.5 um comb is
    not commandable there and silently rounding it would produce a plan the
    hardware cannot execute. The comb therefore widens to the finest step the
    MRES allows. MRES 1 still samples the 4-full-step cycle, just at one point
    per full step (4 per cycle) instead of 16.
    """
    return max(COMB_STEP_UM, FULL_STEP_UM / mres)


def comb_targets(centre_um, mres=4):
    """One electrical-cycle comb, at the finest spacing this MRES allows."""
    step = comb_step_for(mres)
    span = COMB_SPAN_CYCLES * ELECTRICAL_CYCLE_UM
    n = int(round(span / step)) + 1
    return centre_um - span / 2.0 + np.arange(n) * step


class Sequence:
    """Accumulates (time, position) segments with a phase label on each."""

    def __init__(self, x0_um=0.0):
        self.t = 0.0
        self.x = float(x0_um)
        self.segments = []
        self.measurements = []

    def _add(self, kind, phase, dt, x1, rate=None, label=None, measured=False):
        seg = {
            "kind": kind, "phase": phase, "label": label,
            "t_start_s": self.t, "t_end_s": self.t + dt, "duration_s": dt,
            "x_start_um": self.x, "x_end_um": x1,
            "rate": rate, "measured": measured,
        }
        self.segments.append(seg)
        self.t += dt
        self.x = x1
        return seg

    def move_to(self, x_um, phase, rate=MEASURE_RATE, label=None):
        v = RATE_FULLSTEPS_S[rate] * FULL_STEP_UM
        dt = abs(x_um - self.x) / v
        return self._add("move", phase, dt, float(x_um), rate=rate, label=label)

    def dwell(self, dt, phase, label=None, measured=False):
        return self._add("dwell", phase, dt, self.x, label=label,
                         measured=measured)

    def measure_at(self, x_um, phase, label, approach="+"):
        """Move to a target and hold it long enough to measure where it landed.

        `approach` is the direction the target is finally reached from. For "+"
        the stage must already be below it; if it is not, it backs off past the
        target by BACKOFF_UM first, so the final approach is always the same
        direction and backlash cannot enter.
        """
        if approach == "+" and self.x > x_um:
            self.move_to(x_um - BACKOFF_UM, phase, rate=TRAVERSE_RATE,
                         label="back-off")
        elif approach == "-" and self.x < x_um:
            self.move_to(x_um + BACKOFF_UM, phase, rate=TRAVERSE_RATE,
                         label="back-off")
        self.move_to(x_um, phase, rate=MEASURE_RATE, label=label)
        seg = self.dwell(DWELL_MEASURE_S, phase, label=label, measured=True)
        self.measurements.append({
            "phase": phase, "label": label, "target_um": float(x_um),
            "approach": approach,
            "window_start_s": seg["t_start_s"] + SETTLE_DISCARD_S,
            "window_end_s": seg["t_end_s"],
        })
        return seg

    def datum(self, phase="datum", label="datum"):
        """The origin, always entered from below so it means one thing."""
        self.move_to(-BACKOFF_UM, phase, rate=TRAVERSE_RATE, label="back-off")
        self.move_to(0.0, phase, rate=MEASURE_RATE, label=label)
        seg = self.dwell(DATUM_DWELL_S, phase, label=label, measured=True)
        self.measurements.append({
            "phase": phase, "label": label, "target_um": 0.0, "approach": "+",
            "window_start_s": seg["t_start_s"] + SETTLE_DISCARD_S,
            "window_end_s": seg["t_end_s"],
        })
        return seg


def build_warmup(seq):
    """Full-range sweeps until the stage is thermally settled. Not measured.

    v6.1 measures a residual drift of +0.49 nm/s and slow wander of 50-150 nm
    over tens of seconds in a stage that had already been running 19 minutes.
    Starting cold would put a far larger thermal ramp straight through the
    measurement, where it is indistinguishable from a position-dependent error.
    """
    reach = SWEEP_HALF_UM
    v = RATE_FULLSTEPS_S[WARMUP_RATE] * FULL_STEP_UM
    leg_s = 2 * reach / v
    n_legs = max(2, int(round(WARMUP_MIN * 60.0 / leg_s)))
    for i in range(n_legs):
        seq.move_to(reach if i % 2 == 0 else -reach, "warmup",
                    rate=WARMUP_RATE, label="warm-up sweep")
    return n_legs


def build_block(seq, mres):
    """One MRES block: every measured phase, in order."""
    tag = f"MRES {mres}"
    seq.datum(label=f"{tag} datum")

    # --- Unidirectional: absolute accuracy, backlash excluded by construction.
    targets = sweep_targets()
    for i, x in enumerate(targets):
        seq.measure_at(x, "unidirectional", f"{tag} uni #{i + 1}", approach="+")
        if (i + 1) % WITNESS_EVERY == 0 and i + 1 < targets.size:
            seq.datum(phase="witness", label=f"{tag} drift witness")

    # --- Bidirectional: the same targets, now reached from the other side.
    # Forward then backward over one target list; the difference at each target
    # is the reversal error. This is the only phase that deliberately lets
    # backlash in -- it is the quantity being measured.
    for i, x in enumerate(targets):
        seq.measure_at(x, "bidirectional", f"{tag} bi+ #{i + 1}", approach="+")
    for i, x in enumerate(targets[::-1]):
        seq.measure_at(x, "bidirectional",
                       f"{tag} bi- #{targets.size - i}", approach="-")

    # --- Cyclic: map the 4-full-step error at 1/4-step spacing.
    for c, centre in enumerate(COMB_CENTRES_UM):
        comb = comb_targets(centre, mres)
        for j, x in enumerate(comb):
            seq.measure_at(x, "cyclic", f"{tag} comb{c + 1} #{j + 1}",
                           approach="+")
        seq.datum(phase="witness", label=f"{tag} drift witness")

    # --- Settling: the identical move, repeated, for the repeatability figure.
    for r in range(N_SETTLE_REPEATS):
        seq.move_to(SETTLE_TARGET_UM - SETTLE_APPROACH_UM, "settling",
                    rate=TRAVERSE_RATE, label="approach")
        seq.measure_at(SETTLE_TARGET_UM, "settling", f"{tag} repeat #{r + 1}",
                       approach="+")

    seq.datum(phase="witness", label=f"{tag} closing witness")


def build_campaign(mres_values=MRES_VALUES, warmup=True):
    seq = Sequence()
    n_legs = build_warmup(seq) if warmup else 0
    block_bounds = []
    for mres in mres_values:
        t0 = seq.t
        build_block(seq, mres)
        block_bounds.append({"mres": mres, "t_start_s": t0, "t_end_s": seq.t})
    return seq, block_bounds, n_legs


def summarize(seq, block_bounds, n_legs, mres_values):
    targets = sweep_targets()
    phase_time = {p: 0.0 for p in PHASES}
    for s in seq.segments:
        phase_time[s["phase"]] += s["duration_s"]
    # Blocks differ in length: a coarser MRES cannot run as many comb points.
    per_block = {str(b["mres"]): b["t_end_s"] - b["t_start_s"]
                 for b in block_bounds}
    phase_mod = np.sort(targets % LEAD_PITCH_UM)
    return {
        "rig": {
            "full_step_um": FULL_STEP_UM,
            "electrical_cycle_um": ELECTRICAL_CYCLE_UM,
            "lead_pitch_um": LEAD_PITCH_UM,
            "travel_total_um": TRAVEL_TOTAL_UM,
            "guard_um": GUARD_UM,
            "sweep_span_um": 2 * SWEEP_HALF_UM,
        },
        "motion": {
            "measure_rate": MEASURE_RATE,
            "measure_velocity_um_s": RATE_FULLSTEPS_S[MEASURE_RATE] * FULL_STEP_UM,
            "traverse_rate": TRAVERSE_RATE,
            "backoff_um": BACKOFF_UM,
            "note": "one rate for every measured move; the traverse rate is "
                    "used only for unmeasured repositioning",
        },
        "timing": {
            "dwell_measure_s": DWELL_MEASURE_S,
            "settle_discard_s": SETTLE_DISCARD_S,
            "averaging_window_s": DWELL_MEASURE_S - SETTLE_DISCARD_S,
            "samples_per_measurement_at_1khz":
                int(round((DWELL_MEASURE_S - SETTLE_DISCARD_S) * 1000)),
            "datum_dwell_s": DATUM_DWELL_S,
            "warmup_min": WARMUP_MIN,
            "warmup_legs": n_legs,
        },
        "targets": {
            "n_sweep_targets": int(targets.size),
            "sweep_targets_um": [float(v) for v in targets],
            "all_on_full_step": bool(np.allclose(targets % FULL_STEP_UM, 0.0)),
            "min_spacing_um": float(np.min(np.diff(targets))),
            "max_spacing_um": float(np.max(np.diff(targets))),
            "lead_phase_coverage": {
                "note": "targets modulo the 2 mm lead pitch; a uniform grid "
                        "would collapse these to a single value",
                "min_gap_um": float(np.min(np.diff(phase_mod))),
                "max_gap_um": float(np.max(np.diff(phase_mod))),
                "spread_um": float(phase_mod[-1] - phase_mod[0]),
            },
            "comb": {
                "centres_um": list(COMB_CENTRES_UM),
                "target_step_um": COMB_STEP_UM,
                "span_um": COMB_SPAN_CYCLES * ELECTRICAL_CYCLE_UM,
                "mres_required_for_target_step": int(FULL_STEP_UM / COMB_STEP_UM),
                "per_mres": {
                    str(m): {
                        "step_um": comb_step_for(m),
                        "step_fraction_of_full_step":
                            comb_step_for(m) / FULL_STEP_UM,
                        "points_per_comb": int(comb_targets(0.0, m).size),
                        "points_per_electrical_cycle":
                            ELECTRICAL_CYCLE_UM / comb_step_for(m),
                        "limited_by_mres": comb_step_for(m) > COMB_STEP_UM,
                    } for m in mres_values
                },
            },
            "settling_repeats": N_SETTLE_REPEATS,
        },
        "budget": {
            "mres_values": list(mres_values),
            "n_segments": len(seq.segments),
            "n_measurements": len(seq.measurements),
            "per_block_s": per_block,
            "per_block_min": {k: v / 60.0 for k, v in per_block.items()},
            "total_s": seq.t,
            "total_min": seq.t / 60.0,
            "phase_time_s": phase_time,
        },
    }


def export(seq, block_bounds, n_legs, mres_values, out_dir=OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "sequence_plan.csv", "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("t_start_s,t_end_s,duration_s,kind,phase,label,"
                 "x_start_um,x_end_um,rate,measured\n")
        for s in seq.segments:
            fh.write(f"{s['t_start_s']:.4f},{s['t_end_s']:.4f},"
                     f"{s['duration_s']:.4f},{s['kind']},{s['phase']},"
                     f"{s['label'] or ''},{s['x_start_um']:.3f},"
                     f"{s['x_end_um']:.3f},{s['rate'] or ''},"
                     f"{int(s['measured'])}\n")

    with open(out_dir / "sequence_targets.csv", "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("phase,label,target_um,approach,window_start_s,window_end_s\n")
        for m in seq.measurements:
            fh.write(f"{m['phase']},{m['label']},{m['target_um']:.3f},"
                     f"{m['approach']},{m['window_start_s']:.4f},"
                     f"{m['window_end_s']:.4f}\n")

    summary = summarize(seq, block_bounds, n_legs, mres_values)
    summary["blocks"] = block_bounds
    with open(out_dir / "sequence_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def main():
    seq, blocks, n_legs = build_campaign()
    summary = export(seq, blocks, n_legs, MRES_VALUES)
    b = summary["budget"]
    t = summary["targets"]
    print(f"targets          {t['n_sweep_targets']} global, "
          f"spacing {t['min_spacing_um']:.0f}-{t['max_spacing_um']:.0f} um, "
          f"all on full step: {t['all_on_full_step']}")
    print(f"lead phase       spread {t['lead_phase_coverage']['spread_um']:.0f} um "
          f"of {LEAD_PITCH_UM:.0f} um, largest gap "
          f"{t['lead_phase_coverage']['max_gap_um']:.0f} um")
    for m, c in t["comb"]["per_mres"].items():
        flag = "  [limited by MRES]" if c["limited_by_mres"] else ""
        print(f"comb MRES {m:<3}    {c['points_per_comb']:>3} points at "
              f"{c['step_um']:.3f} um "
              f"({c['points_per_electrical_cycle']:.0f}/electrical cycle)"
              f"{flag}")
    print(f"measurements     {b['n_measurements']} across "
          f"{b['n_segments']} segments")
    for m, v in b["per_block_min"].items():
        print(f"block MRES {m:<3}   {v:.1f} min")
    print(f"total campaign   {b['total_min']:.1f} min "
          f"({len(b['mres_values'])} MRES + {WARMUP_MIN:.0f} min warm-up)")
    for f in sorted(OUT_DIR.iterdir()):
        print(f"  wrote {f.name}")


if __name__ == "__main__":
    main()
