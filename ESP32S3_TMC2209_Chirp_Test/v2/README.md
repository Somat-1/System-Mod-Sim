# v2 — resonance-band chirp, hard notch replaced by a smooth taper

## Status

**Design preview only.** This folder currently holds the analytic schedule
and a pre-hardware render of what the excitation looks like
(`scripts/plot_planned_sequence.py` -> `rendered_assets/`), plus a
machine-generated reference of the exact configuration
(`scripts/generate_sweep_config_log.py` -> `SWEEP_CONFIG_LOG.md`). It does
**not** yet contain ESP32/TMC2209 firmware implementing the procedure below
— v1 (`../ESP32S3_TMC2209_Chirp_Test.ino`) still hard-notches 120-230 Hz and
runs a fixed ±1 microstep amplitude at native 1/256 microstepping.

This is revision 2 of the v2 design. Revision 1 derived command amplitude
from an inverted following-error budget (`e_max`) against an assumed SDOF
plant model at MRES=256. It is kept only in git history: it asked for
sub-microstep amplitude above ~816 Hz (not just small — literally
uncommandable), and `e_max` itself was never a satisfying way to reason
about the excitation. Revision 2 (below) picks amplitude directly and drops
to MRES=16.

## Why the notch comes out, and why MRES changed

v1 hard-excludes STEP commands in 120-230 Hz around the modeled resonance,
at native 1/256 microstepping. Two problems with keeping that as-is:

- The exclusion band was set from the current *model's* 176.7 Hz, and a 10%
  error in f_n at Q≈20 is enough to put the notch in the wrong place and
  let the pole slip through unexcited but also unprotected. v2 measures f_n
  first (Section 1) rather than trusting the model.
- **1/256 microstepping is not real resolution near standstill.** The
  driver will happily accept a 1/256 command, but detent torque means the
  rotor cannot be trusted to actually resolve individual microsteps at that
  fineness — the true positioning floor is set by the motor's magnetic
  cogging, not by what MRES the driver electrically supports. v2 uses
  **MRES=16** instead: coarser, but each commanded microstep (1/16 = 0.0625
  full steps = 0.625 µm) is a real, resolvable mechanical event rather than
  a sub-detent-torque fiction.

Given that, v2 does not try to hold a tiny following-error budget through
resonance — at MRES=16 you cannot resolve one anyway. Instead it commands a
directly-chosen "cruise" amplitude almost everywhere, and only winds it
down — smoothly, not a hard cliff — in the band around resonance where the
motor actually struggles. Higher amplitude away from resonance is fine; it
does not need to be justified by a model, only kept inside the two hard
physical clamps in Section 3.

## Does tapering the amplitude bias the frequency response?

Not for a linear plant — the FRF is `H(f) = response(f) / command(f)`, so
whatever amplitude was actually commanded at each frequency cancels out of
the ratio. The real risk is that this stepper/leadscrew system is not
linear (detent torque, friction, backlash), so the measured response can
depend on drive amplitude as well as on frequency — a softening backbone
can shift the apparent peak with level, and a single sweep can't tell you
whether a change in response came from the frequency changing or the
amplitude changing, since the taper moves both together.

The procedure does not try to hide that; it is set up to measure it:

- Three independent cruise levels (Section 3) are run as complete separate
  sweeps, not blended into one variable-amplitude pass, so their FRFs can
  be compared level-to-level.
- Up-sweep vs. down-sweep separation near resonance is the explicit
  quasi-linearity check (the classic jump/bifurcation signature of a
  nonlinear resonance) and is one of the two stop conditions.
- Where possible, the FRF should be built from the *measured* step signal
  (a spare DAQ channel logging the actual command), not from the analytic
  amplitude law, so the taper's shape is never assumed in the result.

## 1. Measure the plant before exciting it

Before every campaign, not just once:

1. Energize at the run current, let it settle, command a single microstep,
   and log the ringdown on the interferometer.
2. f_n comes from the ringing period; ζ from the log decrement over 10-20
   cycles. Average 20-30 ringdowns from alternating directions.
3. Repeat at every current setting used. f_n scales as √I, so 200 mA and
   400 mA runs should sit ~1.41x apart and each needs its own campaign
   (`chirp_v2_schedule.scale_f_n_for_current` is a *planning* aid for this
   ratio only — it is not a substitute for measuring each setting).
4. Redo a short ringdown between runs; if f_n has drifted more than ~2%, the
   motor's thermal state changed and the notch needs re-centering.

`chirp_v2_schedule.F_N_HZ_PLACEHOLDER` (176.7 Hz) is the prior analytical-
model value, kept only so the render below has a number to center the notch
on and so Section 4's dwell-time argument (which also needs `Q_PLACEHOLDER`)
has something to compute with. Replace both with measured values before any
real run — the notch is centered directly on f_n, so an f_n error puts the
notch in the wrong place exactly as it would have for v1's hard cutoff.

## 2. Driver registers

Set explicitly and dump the full register set before and after each run;
any of these changing mid-campaign invalidates the comparison.

- **MRES = 4** (1/16 step) — see "Why the notch comes out, and why MRES
  changed" above. Not native 1/256: detent torque, not the driver, is what
  actually limits resolvable microstepping near standstill.
- **Force SpreadCycle** for the whole run; set the StealthChop velocity
  threshold so no mode switch can occur anywhere in the sweep.
- **IHOLD = IRUN**, hold delay zero, power-down delay at maximum — the
  excitation is centered on standstill, so any standstill current reduction
  changes electromagnetic stiffness and moves f_n mid-run.
- **Fixed chopper timing** — pin TOFF, TBL, HSTRT/HEND once for this motor
  and do not touch them again; chopper frequency should land above 1 kHz.
- StallGuard may stay logged but must not gate anything — it is unreliable
  near standstill and useless as a slip detector here.

## 3. Amplitude: a direct cruise level plus a smooth notch

No following-error model, no inversion. `chirp_v2_schedule.amplitude_profile()`
is:

```
tapered(f) = cruise_full_steps * notch(f)
notch(f)   = 1 - (1 - depth_ratio) * 0.5 * (1 + cos(pi * clip((f - f_n) / halfwidth, -1, 1)))
A_cmd(f)   = min(tapered(f), stroke_clamp, step_rate_clamp / (2*pi*f*MRES))
```

`notch(f)` is a raised-cosine (Hann-shaped) window: exactly 1.0 outside
`f_n ± halfwidth`, exactly `depth_ratio` at `f_n`, and C1-continuous
everywhere in between — zero slope at both edges and at the center, so
there is no kink where the flat cruise region meets the taper. This is a
deliberately simpler thing to reason about than a curve inverted from a
resonance transmissibility (revision 1): "flat, then a smooth cosine dip to
10% right around resonance, then flat again" is the whole model.

| Parameter | Meaning | Value |
|---|---|---:|
| `CRUISE_LEVELS_FULL_STEPS` | Commanded amplitude away from resonance | 3.0 full steps (single level; no escalation) |
| `NOTCH_HALFWIDTH_HZ` | Notch half-width, set to bracket v1's old hard-notch band (120-230 Hz) around f_n so the protected range is not shrinking | 56.7 Hz (-> 120.0-233.4 Hz around the 176.7 Hz placeholder) |
| `NOTCH_DEPTH_RATIO` | Amplitude at f_n, as a fraction of cruise | 0.10 |
| `STROKE_CLAMP_MM` | Hard ceiling from interferometer window / stage travel | 0.5 mm (inactive at 3.0 full steps; would only bind at much higher amplitude) |
| `STEP_RATE_CLAMP_HZ` | Hard ceiling from driver STEP frequency / generator jitter floor | 200,000 Hz (bites for roughly the top third of the sweep, above ~663 Hz — an honest hardware limit, not a modeling artifact, though see the timing note below: this ceiling itself is not yet benchmarked) |

At the notch center, amplitude is `cruise * 0.10` = 0.3 full steps, i.e. 4.8
microsteps at MRES=16 — comfortably above the 1-microstep floor (see
below), unlike revision 1.

**Timing budget, unverified**: an *unclamped* 3.0 full-step cruise would
need a peak microstep pulse rate of ~301.6 kHz at 1000 Hz (`2*pi*f*A*MRES`)
— i.e. the step-rate clamp above is doing real work, cutting the top of the
sweep to fit inside the 200 kHz placeholder. TMC2209's own STEP input
timing has enormous headroom at either number (datasheet minimum STEP
high/low pulse width is well under 1 µs). The open question is the ESP32-S3
firmware's pulse generation: v1's actual `.ino` (the only implementation
that exists) bit-bangs STEP via `digitalWrite()` in a busy-wait loop, and
its own hardcoded minimum inter-edge interval (`DIR_SETUP_US + STEP_HIGH_US
+ 2` = 9 µs) implies only a ~111 kHz ceiling for that specific approach —
below both the 200 kHz placeholder and the 301.6 kHz an unclamped sweep
would ask for. A hardware-timer- or RMT-driven generator would have far
more headroom, but neither that firmware nor a preflight benchmark exists
yet. Until one does, treat 200 kHz as optimistic and prefer bench-testing
the real achievable rate (the same pattern the v4 campaign already uses —
`benchmark_individual_rate()` in
`../../Microstepping Test Data/v4/scripts/run_mres_trajectory_campaign.py`)
before trusting this number on real hardware.

A velocity-based presliding clamp is **not applied by default** in this
revision (it was in revision 1, and is why the excitation there looked much
smaller everywhere). It is fine to run at higher amplitude generally,
provided the band where the motor struggles is wound down — that is the
notch's job, not a global velocity limit. The Stribeck-velocity numbers
that clamp was based on are kept in the module for reference, in case a
*specific* run is deliberately meant to characterize presliding behavior:
source is the Rev 4 LuGre model
(`../../Rev 4 Analytical Model Derivation/index.html`), Stribeck velocities
of 2.0e-4 m/s (nut) and 2.5e-4 m/s (guideway), both flagged there as
pre-emptive estimates reused from Rev 3, not independently identified.

### Microstep quantization floor

One microstep at MRES=16 is `1/16 = 0.0625` full steps. `plot_planned_sequence.py`
prints the minimum amplitude actually reached, in microsteps, for every
cruise level on every run; as of this design that minimum is 4.8 microsteps
(right at f_n) — never below 1, unlike revision 1's MRES=256 design, which
dropped below 1 microstep above ~816 Hz. Re-check
this print after changing any of `CRUISE_LEVELS_FULL_STEPS`,
`NOTCH_DEPTH_RATIO`, or `MRES`; a few microsteps is still coarse
quantization (a commanded "sine" is a visible staircase, not a clean tone),
so do not read "above the floor" as "clean" — only as "not physically
impossible."

## 4. Sweep timing

Two-segment law, per direction:

| Segment | Range | Duration | Notes |
|---|---|---:|---|
| Log | 1 -> 60 Hz | 120 s | ~20 s/octave |
| Linear | 60 -> 1000 Hz | 120 s | 7.83 Hz/s |

The 1 Hz start is not load-bearing — anywhere in 1-5 Hz is fine; the band
below 1 Hz is not a priority for this test. Kept at 1 Hz here only because
that is what the original spec stated; change `LOG_START_HZ` if a higher
floor is preferred.

At that rate the dwell across the half-power bandwidth at resonance is
~1.1 s (~30 time constants at Q=20) — deliberate full steady state; the
notch carries the safety margin, not the sweep rate. Re-verify this once
f_n and Q are actually measured (`Q_PLACEHOLDER` is still used here even
though it no longer feeds the amplitude profile).

Per amplitude level: 30 s idle pre-roll, sync marker, up-sweep, 5 s dwell,
sync marker, down-sweep, 30 s tail —
`chirp_v2_schedule.LEVEL_DURATION_S` = 546.0 s (~9.1 min). Three levels is
roughly half an hour with the ringdown characterization.

The sync marker is three 0.1 s bursts at a fixed frequency clear of
resonance and ambient tones (placeholder 750 Hz), separated by 0.1 s gaps.
**Make each sweep net-zero in displacement**: compare interferometer
position at the start and end of every segment; any residual offset near a
multiple of the full-step pitch (lead ÷ 200 = 0.010 mm) means a pole
slipped and the run is void.

## 5. DAQ

- Sample at ≥10 kHz with anti-alias engaged (chopper harmonics otherwise
  fold into the band).
- Accelerometers and interferometer share one time base.
- Treat anything within 3x of the accelerometer's high-pass corner
  (typically 0.5-1 Hz for IEPE units) as unusable; rely on the interferometer
  there.
- Two pre-roll noise floors: driver energized/not stepping, then driver
  disabled, to separate self-inflicted floor from ambient.
- **The existing Bode pipeline needs a normalization step before it can be
  reused on v2 data.** `../analysis/scripts/process_chirp_bode.py` (built
  for v1) plots raw synchronous-detection response magnitude directly,
  because v1 commands the same fixed amplitude at every frequency, so raw
  magnitude is already proportional to the true FRF there. v2's amplitude
  is not constant (cruise, then a 10x notch dip, then cruise again), so
  running v2 data through that script unmodified would produce a curve
  shaped partly by the excitation's own taper, not purely by the plant.
  Divide the measured response at each frequency by the actual commanded
  amplitude at that frequency — ideally from a logged command/step channel,
  not the analytic profile — before treating the result as `H(f)`.
  `SWEEP_CONFIG_LOG.md` (below) exists specifically to make that
  normalization step reproducible without re-deriving it from this README.

## 6. Guardrails and run order

Warm up at IRUN for 10-15 minutes before the first run. Soft travel limits
well inside hard stops. Confirm screw backdrivability before relying on
driver-disable as an abort. Run order: ringdowns at each current, then a
scout up-sweep at 3.0 full steps with the down-sweep suppressed, check for
net offset, then enable the down-sweep. Always run up before down — the
softening backbone makes the down direction the one that jumps onto the
high-amplitude branch. Either stop condition (a net displacement offset, or
up/down curves separating near resonance) means stop and reduce
`CRUISE_LEVELS_FULL_STEPS` before trying again — with only one level
configured now there is no automatic lower rung to fall back to, so a
failure here is a real stop, not a step down to a pre-planned lower level.

## Files

- `scripts/chirp_v2_schedule.py` — the analytic design: plant placeholders,
  the cruise/notch amplitude profile and its two hard clamps, and the
  two-segment sweep law. No serial or hardware I/O.
- `scripts/plot_planned_sequence.py` — renders `rendered_assets/planned_chirp_v2_excitation.png`,
  four stacked panels:
  1. **Full-run envelope** — peak commanded displacement vs. time across the
     *entire* ~9.1 min level (idle, both sync markers, both sweep
     directions, mid-dwell, tail), log-scaled, with a right-hand axis in
     microsteps and the 1-microstep floor marked. The carrier itself is not
     drawn here — at this time scale it would just be a dense fill — only
     its peak envelope. Flat at 3.0 full steps away from resonance, with an
     honest step-rate-clamp roll-off near the top of each sweep (above
     ~663 Hz) and the smooth notch dip at each resonance crossing. Plots
     every level in `CRUISE_LEVELS_FULL_STEPS` if more than one is
     configured again later.
  2. Beginning-of-run window (idle tail, sync marker, first ~20 s of the
     log sweep) with the actual oscillating command, not just its envelope.
  3. Mid-run window centered on the placeholder resonance crossing inside
     the linear segment, showing the smooth notch shape directly.
  4. Turnaround window at 1000 Hz (last second of the up-sweep, the 5 s
     mid-dwell, the second sync marker, first second of the down-sweep
     starting back down from 1000 Hz) — at MRES=16 cruise amplitude and the
     sync marker are now comparable in size, so both are visible on one
     scale without the clipping revision 1 needed.
  Panels 2-4 show the real carrier because that is only tractable over a
  few seconds at a time; panel 1 is the only view of the complete run, and
  it deliberately only shows the envelope for that reason.
- `scripts/generate_sweep_config_log.py` — introspects `chirp_v2_schedule.py`
  and (re)writes `SWEEP_CONFIG_LOG.md`. Run it after changing any constant
  in the schedule module so the log never drifts from the code that
  actually generated a given run's commanded signal.
- `SWEEP_CONFIG_LOG.md` — generated, not hand-edited. The reference to pull
  up when writing or reviewing the processing pipeline for a real recording:
  exact amplitude/notch/clamp/timing parameters, without re-reading this
  narrative README or the module source.
- `data/hardware_runs/` — not created yet; will hold ringdown and sweep logs
  once there is real hardware to log from.

## Rendering the preview

```powershell
python scripts\plot_planned_sequence.py
python scripts\generate_sweep_config_log.py
```

The first prints the placeholder resonance crossing time and the minimum
microstep count reached per cruise level, and saves the PNG above. The
second regenerates `SWEEP_CONFIG_LOG.md`. Edit the placeholders in
`chirp_v2_schedule.py` (f_n, the cruise levels, the notch and clamp values)
once real numbers are measured, and re-run both.
