#!/usr/bin/env python3
"""The two chirp sweep profiles, defined as data.

A mechanical frequency sweep of the stage: the commanded position oscillates at
a frequency that walks from low to high and back, so one run covers the band and
the response can be read against frequency.

Two profiles, because of something that happened on this rig:

  FULL  1 Hz -> 1 kHz and back, **with a hard 150-250 Hz gap**. The band is not
        visited at all. See `SAFETY` below - this is not a tuning choice.
  HIGH  250 Hz -> 1 kHz and back, continuous. Entirely above the gap, so it
        needs no gap of its own, and it starts with a dwell at 250 Hz instead
        of arriving there by a jump.

Both keep v3's shared sweep rate of 7.8125 Hz/s and its 750 Hz sync marker
bursts. Amplitude (2026-10-13): a constant +-1.875 full steps (30 microsteps at
MRES 16) at 500 mA across the whole band (first tried at +-1.5 / 400 mA, which
gave a small signal). It is under the 2-full-step lag at which a
stepper slips (v3's +-3 full steps was not), and under the step-rate ceiling at
1 kHz with 6% margin, so it is never clamped. A constant-acceleration law
(+-1 full step rolling off as 1/f^2 above 100 Hz) was tried first and dropped:
at 250-1000 Hz it commanded 0.01-0.16 full steps, inside the motor's friction
dead band, and the accelerometer showed nothing.

SAFETY
------
Two hardware incidents on 2026-09-15 produced uncontrolled motor motion,
reported live at around 200 Hz. A smooth amplitude notch over that band did not
prevent them - not at 90% suppression, and not when widened to 65-285 Hz at 97%.
Reducing the amplitude there was *not sufficient, however aggressively*. That is
why FULL has a genuine frequency gap rather than a notch: between 150 and 250 Hz
it commands nothing at all.

Do not narrow the gap to "get the missing data". The missing data is the point.
If the 150-250 Hz band has to be characterised, it needs a different method than
driving the motor through it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

# --- carried over from the hardware-validated v3 chirp ---------------------
MRES = 8
CRUISE_FULL_STEPS = 2.5                               # above the 2-fs slip limit
CRUISE_MICROSTEPS = CRUISE_FULL_STEPS * MRES          # 20
F_CORNER_HZ = float("inf")                            # no roll-off: constant
MAX_STEP_RATE = 200_000.0                             # firmware ceiling, microsteps/s
CURRENT_RMS_MA = 1000                                 # run 8; firmware default 500
SWEEP_RATE_HZ_PER_S = 7.8125                          # (1000-250)/96

# The band that caused the incidents. FULL jumps across it; HIGH starts above it.
RESONANCE_GAP = (150.0, 250.0)



def amplitude(f_hz: float) -> float:
    """Commanded amplitude in microsteps at MRES; must match the firmware."""
    if f_hz <= F_CORNER_HZ:
        return CRUISE_MICROSTEPS
    return CRUISE_MICROSTEPS * (F_CORNER_HZ / f_hz) ** 2


SYNC_FREQ_HZ = 750.0
SYNC_BURST_S = 0.1
SYNC_GAP_S = 0.1
SYNC_BURSTS = 3
SYNC_MARKER_S = SYNC_BURSTS * SYNC_BURST_S + (SYNC_BURSTS - 1) * SYNC_GAP_S  # 0.5

PRE_ROLL_S = 30.0          # long, so the recorder is definitely running
TAIL_S = 5.0
MID_DWELL_S = 2.0
SETTLE_S = 3.0             # HIGH only: dwell at the start frequency


@dataclass
class Segment:
    kind: str              # silent | marker | sweep | dwell
    label: str
    duration_s: float
    f_start: float = 0.0
    f_end: float = 0.0
    amplitude_microsteps: float = 0.0
    note: str = ""


@dataclass
class Profile:
    key: str
    title: str
    purpose: str
    segments: list = field(default_factory=list)

    @property
    def duration_s(self) -> float:
        return sum(s.duration_s for s in self.segments)

    def timeline(self):
        t, out = 0.0, []
        for s in self.segments:
            out.append((t, t + s.duration_s, s))
            t += s.duration_s
        return out


def _sweep_time(f0: float, f1: float) -> float:
    return abs(f1 - f0) / SWEEP_RATE_HZ_PER_S


def _marker(label: str) -> Segment:
    return Segment("marker", label, SYNC_MARKER_S, SYNC_FREQ_HZ, SYNC_FREQ_HZ,
                   amplitude(SYNC_FREQ_HZ),
                   note=f"{SYNC_BURSTS} bursts at {SYNC_FREQ_HZ:g} Hz for alignment")


def _top_marker() -> list:
    """The 1 kHz turn-around: silence, the sync marker, silence."""
    return [Segment("silent", "top gap", 1.0, note="up sweep done; 1 kHz reached"),
            _marker("top marker (1 kHz turn-around)"),
            Segment("silent", "top gap", 1.0, note="down sweep follows")]


def full_profile() -> Profile:
    """1 Hz to 1 kHz and back, skipping 150-250 Hz entirely.

    The jump across the gap is an instantaneous change of target frequency.
    Bounded by the cruise amplitude it is not a new class of event for the step
    generator - but it *is* a discontinuity in the excitation, and anything odd
    seen right at 250 Hz should be suspected of being this before it is
    attributed to the stage. That is what the HIGH profile is for.
    """
    lo, hi = RESONANCE_GAP
    p = Profile("full", "FULL — 1 Hz to 1 kHz, gap 150–250 Hz",
                "cover the whole band except the one that is unsafe to drive")
    p.segments += [
        Segment("silent", "pre-roll", PRE_ROLL_S, note="no excitation; start the recorder"),
        _marker("sync marker 1"),
        Segment("sweep", "up: low leg", _sweep_time(1.0, lo), 1.0, lo, CRUISE_MICROSTEPS),
        Segment("dwell", "gap jump", 0.0, lo, hi, 0.0,
                note=f"instantaneous {lo:g}->{hi:g} Hz; nothing commanded in between"),
        Segment("sweep", "up: high leg", _sweep_time(hi, 1000.0), hi, 1000.0, CRUISE_MICROSTEPS),
        *_top_marker(),
        Segment("sweep", "down: high leg", _sweep_time(1000.0, hi), 1000.0, hi, CRUISE_MICROSTEPS),
        Segment("dwell", "gap jump", 0.0, hi, lo, 0.0,
                note=f"instantaneous {hi:g}->{lo:g} Hz"),
        Segment("sweep", "down: low leg", _sweep_time(lo, 1.0), lo, 1.0, CRUISE_MICROSTEPS),
        _marker("sync marker 2"),
        Segment("silent", "tail", TAIL_S),
    ]
    return p


def high_profile() -> Profile:
    """250 Hz to 1 kHz and back, continuous — and started cleanly.

    Two differences from FULL's high leg, both deliberate:

    * It never goes near 150-250 Hz, so there is no gap and no jump.
    * It **dwells at 250 Hz before sweeping**, instead of arriving there
      discontinuously. In FULL, 250 Hz is exactly the frequency the sweep jumps
      back in at, and an instantaneous frequency change injects a transient
      across a wide band. If something odd appears at 250 Hz in FULL and not
      here, it was the jump; if it appears in both, it is the stage.
    """
    _, start = RESONANCE_GAP
    p = Profile("high", f"HIGH — {start:g} Hz to 1 kHz, continuous",
                "resolve whether the 250 Hz feature is the stage or the gap edge")
    p.segments += [
        Segment("silent", "pre-roll", PRE_ROLL_S, note="no excitation; start the recorder"),
        _marker("sync marker 1"),
        Segment("dwell", f"settle at {start:g} Hz", SETTLE_S, start, start,
                CRUISE_MICROSTEPS,
                note="steady excitation before the sweep, so the start is not a jump"),
        Segment("sweep", "up", _sweep_time(start, 1000.0), start, 1000.0, CRUISE_MICROSTEPS),
        Segment("dwell", "mid dwell", MID_DWELL_S, 1000.0, 1000.0, CRUISE_MICROSTEPS),
        *_top_marker(),
        Segment("sweep", "down", _sweep_time(1000.0, start), 1000.0, start, CRUISE_MICROSTEPS),
        Segment("dwell", f"settle at {start:g} Hz", SETTLE_S, start, start,
                CRUISE_MICROSTEPS),
        _marker("sync marker 2"),
        Segment("silent", "tail", TAIL_S),
    ]
    return p


def profiles():
    return [full_profile(), high_profile()]


def summarise(ps):
    out = {
        "sweep_rate_hz_per_s": SWEEP_RATE_HZ_PER_S,
        "amplitude_law": "constant A0 at every frequency",
        "a0_full_steps": CRUISE_FULL_STEPS,
        "a0_microsteps": CRUISE_MICROSTEPS,
        "peak_step_rate_at_1khz": round(2 * 3.141592653589793 * 1000 * CRUISE_MICROSTEPS),
        "amplitude_microsteps_at_hz": {str(f): round(amplitude(f), 3)
                                       for f in (1, 10, 50, 100, 150, 250, 500, 750, 1000)},
        "mres": MRES,
        "current_rms_ma": CURRENT_RMS_MA,
        "resonance_gap_hz": list(RESONANCE_GAP),
        "safety": "150-250 Hz caused uncontrolled motion on 2026-09-15; "
                  "FULL omits the band entirely rather than attenuating it",
        "profiles": {},
    }
    for p in ps:
        sweeps = [s for s in p.segments if s.kind == "sweep"]
        covered = sorted({(min(s.f_start, s.f_end), max(s.f_start, s.f_end))
                          for s in sweeps})
        out["profiles"][p.key] = {
            "title": p.title,
            "purpose": p.purpose,
            "duration_s": round(p.duration_s, 2),
            "duration_min": round(p.duration_s / 60.0, 2),
            "segments": len(p.segments),
            "bands_hz": [[a, b] for a, b in covered],
            "sweep_seconds": round(sum(s.duration_s for s in sweeps), 2),
        }
    return out


def main():
    import argparse
    import csv
    from pathlib import Path

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = Path(args.out) if args.out else Path(__file__).resolve().parent.parent / "sequence"
    out.mkdir(parents=True, exist_ok=True)

    ps = profiles()
    rows = []
    for p in ps:
        for t0, t1, s in p.timeline():
            rows.append(dict(profile=p.key, kind=s.kind, label=s.label,
                             t_start=round(t0, 4), t_end=round(t1, 4),
                             duration_s=round(s.duration_s, 4),
                             f_start=s.f_start, f_end=s.f_end,
                             amplitude_microsteps=round(amplitude(s.f_start), 3)
                             if s.amplitude_microsteps else 0.0,
                             amplitude_end=round(amplitude(s.f_end), 3)
                             if s.amplitude_microsteps else 0.0,
                             note=s.note))
    with open(out / "chirp_plan.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    summary = summarise(ps)
    (out / "chirp_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                            encoding="utf-8")
    for p in ps:
        info = summary["profiles"][p.key]
        bands = " + ".join(f"{a:g}-{b:g} Hz" for a, b in info["bands_hz"])
        print(f"  {p.key.upper():<5} {info['duration_min']:5.2f} min  "
              f"{info['segments']:2d} segments   {bands}")
    print(f"wrote {out/'chirp_plan.csv'} and chirp_summary.json")


if __name__ == "__main__":
    main()
