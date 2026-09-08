#!/usr/bin/env python3
'''Analytic design for the v2 resonance-band chirp (hard notch replaced).

Revision 2 of this module. Revision 1 derived the command amplitude from an
inverted following-error budget (e_max) against an assumed SDOF plant model
-- physically motivated, but it produced two problems in practice: (1) it
was not intuitive ("what does e_max actually mean" was never a fully
satisfying answer), and (2) at MRES=256 it asked for sub-microstep amplitude
above ~816 Hz, which is not just "small" but literally uncommandable, and
detent torque makes the real quantization floor higher than 1/256 anyway --
native 1/256 microstepping is not something this motor can be trusted to
resolve mechanically near standstill regardless of what the driver accepts
electrically.

This revision drops both the inverted-model amplitude and the fine MRES:

- **MRES = 16** (1/16 step). Coarser than native 1/256, but each microstep
  is large enough to be a real, resolvable mechanical event rather than a
  sub-detent-torque fiction. 1 microstep = 1/16 = 0.0625 full steps.
- **Amplitude is chosen directly**, not inverted from a target following
  error. A fixed "cruise" amplitude (in full steps) is commanded across
  most of the band; a smooth, symmetric window multiplies it down to a
  small fraction only within a defined band around resonance, then back up.
  This is deliberately simple: cruise amplitude and notch depth are numbers
  you pick and can directly picture, not the output of inverting a plant
  model you have not measured yet.
- **The notch is smooth by construction**, not because a resonance curve
  happened to have that shape. A raised-cosine (Hann-shaped) window is
  C1-continuous at both edges and at the center, so there is no kink where
  the flat cruise region meets the taper -- unlike a following-error curve
  inverted from a resonance transmissibility, whose slope genuinely does
  change character where the plant model's asymptotic regimes hand off to
  each other. The window's *width* is set to bracket the v1 hard-notch band
  (120-230 Hz) so the protected frequency range is unchanged from v1; only
  its edges are smoothed.
- **Two hard physical clamps remain**: a stroke ceiling (interferometer
  window / stage travel) and a step-rate ceiling (driver STEP frequency /
  generator jitter floor). Both are real "cannot exceed" limits, unlike the
  old velocity/presliding clamp, which is no longer applied by default --
  it is fine to run at a higher amplitude generally, provided the resonance
  band itself is wound down (the notch's job). The Stribeck-velocity numbers
  it was based on are kept below for reference, in case a *specific* run is
  deliberately meant to characterize presliding behavior.

No serial/hardware I/O lives here. This is the same "plan first, render
before flashing" split used by the v4 MRES/trajectory campaign
(../../Microstepping Test Data/v4/scripts/run_mres_trajectory_campaign.py).

f_n below is the prior analytical-model placeholder (176.7 Hz), carried
over ONLY so this module can render a representative plot before hardware
exists, and only to place and size the notch -- it is no longer used to
compute amplitude magnitude anywhere. Every real campaign must replace it
with a same-day ringdown measurement per current setting (Section 1 of the
README); nothing here should be read back into firmware unchanged.
'''

from __future__ import annotations

import math

import numpy as np

# --- Mechanics / drive -------------------------------------------------

LEAD_MM_PER_REV = 2.0
FULL_STEPS_PER_REV = 200
FULL_STEP_MM = LEAD_MM_PER_REV / FULL_STEPS_PER_REV  # 0.010 mm

# 1/16, not native 1/256: each microstep needs to be a real mechanical
# event. At 1/256, detent torque means the rotor cannot be trusted to
# resolve individual commanded microsteps near standstill regardless of
# what the driver accepts electrically -- the effective floor is set by the
# motor, not by MRES. 1/16 is coarser but the motion it commands is real.
MRES = 16
ONE_MICROSTEP_FULL_STEPS = 1.0 / MRES

# --- Plant placeholders (replace with same-day ringdown measurements) --

CURRENT_RMS_MA = 400.0
F_N_HZ_PLACEHOLDER = 176.7  # Prior model value; DO NOT trust for run safety.
Q_PLACEHOLDER = 20.0  # Only used for the Section-4 dwell-time argument now.
ZETA_PLACEHOLDER = 1.0 / (2.0 * Q_PLACEHOLDER)

# f_n scales ~ sqrt(I); a 200 mA run needs its own measured f_n, not this one.
CURRENT_SCALE_REFERENCE_MA = 400.0


def scale_f_n_for_current(f_n_hz: float, current_ma: float) -> float:
    '''Project a measured f_n to another current via the sqrt(I) rule.

    Only a planning aid: every real current setting still needs its own
    ringdown per the memo. Do not substitute this for a measurement.
    '''
    return f_n_hz * math.sqrt(current_ma / CURRENT_SCALE_REFERENCE_MA)


# --- Amplitude: directly-chosen cruise level + a smooth notch -----------

# Directly the peak command amplitude in full steps away from resonance --
# no inversion, no model. Single level: escalating 1.0 -> 2.0 -> 3.0 full
# steps was the revision-1/e_max-era safety ramp; 3.0 full steps is used
# directly now. Re-introduce lower levels here (a tuple) if a cautious
# first pass at lower amplitude is wanted again for a specific motor/setup.
CRUISE_LEVELS_FULL_STEPS = (3.0,)

# The notch is centered on f_n and its half-width is set to bracket v1's
# hard-notch band (120-230 Hz) on both sides, so the protected range is not
# shrinking relative to what was already validated -- only its edges are
# now smooth instead of a cliff.
_V1_NOTCH_LOW_HZ = 120.0
_V1_NOTCH_HIGH_HZ = 230.0
NOTCH_HALFWIDTH_HZ = max(
    F_N_HZ_PLACEHOLDER - _V1_NOTCH_LOW_HZ, _V1_NOTCH_HIGH_HZ - F_N_HZ_PLACEHOLDER
)
NOTCH_DEPTH_RATIO = 0.10  # Amplitude at f_n = 10% of cruise, not nearly zero.

# Hard physical clamps (Section 3). Stroke is inactive at these cruise
# levels (50 full steps vs. 3) but is kept as a real ceiling; step-rate
# does bite near the top of the sweep -- that roll-off is meant to be an
# honest driver/generator limit, not a modeling artifact, but the 200 kHz
# figure below is NOT yet benchmarked against real firmware. TMC2209's own
# STEP input timing has enormous headroom here (datasheet minimum STEP
# high/low well under 1 us), so it is not the bottleneck. The ESP32-S3
# firmware's pulse-generation approach is: v1's actual .ino (the only
# implementation that exists so far) bit-bangs STEP via digitalWrite() in a
# busy-wait loop, and its own hardcoded minimum inter-edge interval
# (DIR_SETUP_US + STEP_HIGH_US + 2 = 9 us) implies a ~111 kHz ceiling for
# that approach specifically -- well below both this placeholder and the
# ~301.6 kHz an unclamped cruise=3.0 sweep would ask for at 1000 Hz. A
# hardware-timer- or RMT-driven generator (not yet written) would have much
# more headroom, but until either that exists or this number is measured
# with a preflight benchmark (the same pattern the v4 campaign already uses
# -- see benchmark_individual_rate() in
# ../../Microstepping Test Data/v4/scripts/run_mres_trajectory_campaign.py),
# treat 200 kHz as unverified and prefer the more conservative ~111 kHz if
# in doubt.
STROKE_CLAMP_MM = 0.5
STEP_RATE_CLAMP_HZ = 200_000.0

# Not applied by default (see module docstring) -- kept for a run that
# specifically wants to stay in presliding. Source: Rev 4 LuGre model
# (../../Rev 4 Analytical Model Derivation/index.html), Stribeck velocities
# of 2.0e-4 m/s (nut) and 2.5e-4 m/s (guideway), both flagged there as
# pre-emptive estimates reused from Rev 3, not independently identified.
STRIBECK_VELOCITY_NUT_MM_S = 0.20
STRIBECK_VELOCITY_GUIDEWAY_MM_S = 0.25


def resonance_notch_taper(
    f_hz: np.ndarray,
    *,
    f_n_hz: float = F_N_HZ_PLACEHOLDER,
    halfwidth_hz: float = NOTCH_HALFWIDTH_HZ,
    depth_ratio: float = NOTCH_DEPTH_RATIO,
) -> np.ndarray:
    '''Raised-cosine multiplier: 1.0 outside the notch, depth_ratio at f_n.

    C1-continuous everywhere (zero slope at both edges and at the center),
    unlike a curve inverted from a resonance transmissibility -- there is no
    frequency at which this window's character changes abruptly.
    '''
    f_hz = np.asarray(f_hz, dtype=float)
    x = np.clip((f_hz - f_n_hz) / halfwidth_hz, -1.0, 1.0)
    window = 0.5 * (1.0 + np.cos(np.pi * x))  # 1 at center (x=0), 0 at edges
    return 1.0 - (1.0 - depth_ratio) * window


def amplitude_profile(
    f_hz: np.ndarray,
    cruise_full_steps: float,
    *,
    f_n_hz: float = F_N_HZ_PLACEHOLDER,
    halfwidth_hz: float = NOTCH_HALFWIDTH_HZ,
    depth_ratio: float = NOTCH_DEPTH_RATIO,
    stroke_clamp_mm: float = STROKE_CLAMP_MM,
    step_rate_clamp_hz: float = STEP_RATE_CLAMP_HZ,
) -> np.ndarray:
    '''Peak command amplitude (full steps): cruise * notch, clamped twice.'''
    f_hz = np.asarray(f_hz, dtype=float)
    omega = 2.0 * np.pi * f_hz
    tapered = cruise_full_steps * resonance_notch_taper(
        f_hz, f_n_hz=f_n_hz, halfwidth_hz=halfwidth_hz, depth_ratio=depth_ratio,
    )
    a_stroke = stroke_clamp_mm / FULL_STEP_MM
    a_step_rate = step_rate_clamp_hz / (omega * MRES)
    return np.minimum(np.minimum(tapered, a_stroke), a_step_rate)


# --- Sweep timing (Section 4) --------------------------------------------

LOG_START_HZ = 1.0
LOG_END_HZ = 60.0
LOG_DURATION_S = 120.0

LINEAR_START_HZ = LOG_END_HZ
LINEAR_END_HZ = 1000.0
LINEAR_DURATION_S = 120.0

SWEEP_DURATION_S = LOG_DURATION_S + LINEAR_DURATION_S  # 240 s, one direction

PRE_ROLL_S = 30.0
SYNC_BURST_S = 0.1
SYNC_GAP_S = 0.1
SYNC_BURST_COUNT = 3
SYNC_MARKER_DURATION_S = (
    SYNC_BURST_COUNT * SYNC_BURST_S + (SYNC_BURST_COUNT - 1) * SYNC_GAP_S
)
SYNC_MARKER_FREQ_HZ = 750.0  # Clear of f_n and of 50/60 Hz mains harmonics.
SYNC_MARKER_AMPLITUDE_FULL_STEPS = 2.0

MID_DWELL_S = 5.0
TAIL_S = 30.0

LEVEL_DURATION_S = (
    PRE_ROLL_S + SYNC_MARKER_DURATION_S + SWEEP_DURATION_S + MID_DWELL_S
    + SYNC_MARKER_DURATION_S + SWEEP_DURATION_S + TAIL_S
)


# --- Frequency-vs-time law (Section 4) ------------------------------------

def up_sweep_frequency(t_s: np.ndarray) -> np.ndarray:
    '''Instantaneous frequency for the ascending sweep, t in [0, SWEEP_DURATION_S].'''
    t_s = np.asarray(t_s, dtype=float)
    in_log = t_s < LOG_DURATION_S
    log_ratio = LOG_END_HZ / LOG_START_HZ
    f_log = LOG_START_HZ * log_ratio ** (t_s / LOG_DURATION_S)
    t_lin = t_s - LOG_DURATION_S
    f_lin = LINEAR_START_HZ + (LINEAR_END_HZ - LINEAR_START_HZ) * (
        t_lin / LINEAR_DURATION_S
    )
    return np.where(in_log, f_log, f_lin)


def down_sweep_frequency(t_s: np.ndarray) -> np.ndarray:
    '''Instantaneous frequency for the descending sweep (mirror of the up-sweep).'''
    t_s = np.asarray(t_s, dtype=float)
    in_linear = t_s < LINEAR_DURATION_S
    f_lin = LINEAR_END_HZ - (LINEAR_END_HZ - LINEAR_START_HZ) * (
        t_s / LINEAR_DURATION_S
    )
    t_log = t_s - LINEAR_DURATION_S
    log_ratio = LOG_START_HZ / LOG_END_HZ
    f_log = LOG_END_HZ * log_ratio ** (t_log / LOG_DURATION_S)
    return np.where(in_linear, f_lin, f_log)


def integrate_phase_rad(t_s: np.ndarray, f_hz: np.ndarray) -> np.ndarray:
    '''Cumulative phase (rad) from instantaneous frequency via trapezoidal integration.'''
    from scipy.integrate import cumulative_trapezoid
    return cumulative_trapezoid(2.0 * np.pi * f_hz, t_s, initial=0.0)


def _render_sweep_command(
    t_s: np.ndarray,
    frequency_law,
    *,
    cruise_full_steps: float,
    f_n_hz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Commanded position (full steps) at times t_s (from that sweep's own start).

    Returns (frequency_hz, amplitude_full_steps, position_full_steps). Phase
    is integrated fresh from t_s[0], so each direction's phase reference is
    independent -- there is no continuity requirement across the mid-sweep
    dwell and sync marker that separates them.
    '''
    f_hz = frequency_law(t_s)
    amplitude = amplitude_profile(f_hz, cruise_full_steps, f_n_hz=f_n_hz)
    phase_rad = integrate_phase_rad(t_s, f_hz)
    position = amplitude * np.sin(phase_rad)
    return f_hz, amplitude, position


def render_up_sweep_command(
    t_s: np.ndarray,
    *,
    cruise_full_steps: float,
    f_n_hz: float = F_N_HZ_PLACEHOLDER,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Commanded position (full steps) across the up-sweep at times t_s (from sweep start).'''
    return _render_sweep_command(
        t_s, up_sweep_frequency, cruise_full_steps=cruise_full_steps,
        f_n_hz=f_n_hz,
    )


def render_down_sweep_command(
    t_s: np.ndarray,
    *,
    cruise_full_steps: float,
    f_n_hz: float = F_N_HZ_PLACEHOLDER,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Commanded position (full steps) across the down-sweep at times t_s (from sweep start).'''
    return _render_sweep_command(
        t_s, down_sweep_frequency, cruise_full_steps=cruise_full_steps,
        f_n_hz=f_n_hz,
    )


def sync_marker_command(t_s: np.ndarray) -> np.ndarray:
    '''Three short bursts at SYNC_MARKER_FREQ_HZ, separated by silent gaps.'''
    t_s = np.asarray(t_s, dtype=float)
    position = np.zeros_like(t_s)
    period = SYNC_BURST_S + SYNC_GAP_S
    for burst_index in range(SYNC_BURST_COUNT):
        burst_start = burst_index * period
        burst_end = burst_start + SYNC_BURST_S
        in_burst = (t_s >= burst_start) & (t_s < burst_end)
        local_t = t_s[in_burst] - burst_start
        position[in_burst] = SYNC_MARKER_AMPLITUDE_FULL_STEPS * np.sin(
            2.0 * np.pi * SYNC_MARKER_FREQ_HZ * local_t
        )
    return position
