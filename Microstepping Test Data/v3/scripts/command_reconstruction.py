#!/usr/bin/env python3
"""Shared "commanded trajectory" reconstruction for Block D plateaus.

Software-paced rates (0.125, 0.375, 1.25 full-steps/s) genuinely execute as
one MA/MR command per single microstep -- the logged MOVE_ACK stream already
*is* the real commanded trajectory (a staircase of individual pulses) and
needs no further reconstruction; a zero-order hold between consecutive
MOVE_ACK events is exact.

Controller-paced rates (3.5, 9.5, 27.5, 70, 200 full-steps/s) execute as a
*single* MA/MR command per direction; the controller then runs its own
ACCEL/CONST/DECEL ramp in hardware (SS ACCEL=DECEL=628, RAMPTYPE=1). Only the
issue instant and the final target are logged, so a zero-order hold between
those two points draws an instantaneous step where the stage actually took
several seconds to ramp up, hold speed, and ramp down. See README.md,
"Controller-paced D-rate command reconstruction".

This module reconstructs that trapezoid analytically from the same
ACCEL/rate parameters the controller itself used
(run_identification_dedicated_controller.py: speed_code, plateau_duration_s,
PLATEAU_ACCEL_CODE), timed so it starts and ends exactly on the two real
logged MOVE_ACK/window-boundary instants -- so the small number of
(MRES, rate) combinations already known to complete faster than the
theoretical ramp (README.md's D-rate execution-timing anomaly) get a
correctly time-scaled trapezoid rather than an overlapping/out-of-order one.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_identification_dedicated_controller import (  # noqa: E402
    CONDITIONING_FULL_STEPS_S,
    MOTOR_FULL_STEPS_PER_REV,
    PLATEAU_ACCEL_CODE,
    SLOW_PLATEAU_RATES_FS_S,
    plateau_duration_s,
    speed_code,
)

# Blocks whose moves never stamp their own rate_full_steps_s (only Block D's
# MOVE_ACK rows do) but which are known, from the driver script, to run at a
# single fixed constant-velocity rate with configure_speed(..., constant_start
# =True) -- i.e. SS min=max, no acceleration phase at all. Without this, these
# moves fall back to an instantaneous zero-order-hold jump, which is a much
# larger and faster step than the controller ever actually commanded (e.g.
# the C block's approach/return moves execute over ~27 ms at 150 full-steps/s,
# not instantly) -- see README.md, "C-block command reconstruction".
BURST_RATE_BY_BLOCK_FS_S = {
    'C': CONDITIONING_FULL_STEPS_S,
}


def is_controller_paced(rate_full_steps_s_str: str) -> bool:
    """True for MOVE_ACK rows belonging to a controller-paced D-rate plateau
    (as opposed to a software-paced one, or a non-D-rate move that never
    stamps rate_full_steps_s at all)."""
    if not rate_full_steps_s_str:
        return False
    return float(rate_full_steps_s_str) not in SLOW_PLATEAU_RATES_FS_S


def trapezoid_profile(rate_full_steps_s: float):
    """(t_accel_s, duration_s, t_total_s) for one direction of a
    controller-paced plateau move at this rate, from the same formula as
    run_supported_plateau_direction (ACCEL=DECEL=628, RAMPTYPE=1)."""
    maximum_code = speed_code(rate_full_steps_s)
    actual_omega = maximum_code * 0.01
    actual_rate_fs_s = actual_omega * MOTOR_FULL_STEPS_PER_REV / (2.0 * np.pi)
    actual_accel = PLATEAU_ACCEL_CODE * 0.01
    actual_accel_fs_s2 = actual_accel * MOTOR_FULL_STEPS_PER_REV / (2.0 * np.pi)
    t_accel_s = actual_rate_fs_s / actual_accel_fs_s2
    duration_s = plateau_duration_s(rate_full_steps_s)
    return t_accel_s, duration_s, duration_s + 2.0 * t_accel_s


def trapezoid_fraction(t_rel: float, t_accel_s: float, duration_s: float) -> float:
    """Fraction (0..1) of total travel completed t_rel seconds into a
    trapezoid move with ramp time t_accel_s and constant-speed time
    duration_s (unit velocity/acceleration -- only the ratio between the two
    times matters, so no rate/accel-code scaling is needed here)."""
    t_total = duration_s + 2.0 * t_accel_s
    t_rel = min(max(t_rel, 0.0), t_total)
    total_dist = duration_s + t_accel_s
    if total_dist <= 0.0:
        return 1.0
    if t_rel <= t_accel_s:
        dist = 0.5 * t_rel * t_rel / t_accel_s if t_accel_s > 0 else 0.0
    elif t_rel <= t_accel_s + duration_s:
        dist = 0.5 * t_accel_s + (t_rel - t_accel_s)
    else:
        t2 = t_rel - t_accel_s - duration_s
        dist = 0.5 * t_accel_s + duration_s + t2 - 0.5 * t2 * t2 / t_accel_s
    return dist / total_dist


def trapezoid_fraction_array(t_rel: np.ndarray, t_accel_s: float, duration_s: float) -> np.ndarray:
    """Vectorized form of trapezoid_fraction, for sampling a whole solve_ivp
    output array at once (e.g. to derive torque time series post hoc)."""
    t_total = duration_s + 2.0 * t_accel_s
    t_rel = np.clip(t_rel, 0.0, t_total)
    total_dist = duration_s + t_accel_s
    if total_dist <= 0.0:
        return np.ones_like(t_rel)
    accel_end = t_accel_s
    cruise_end = t_accel_s + duration_s
    dist = np.empty_like(t_rel)
    in_accel = t_rel <= accel_end
    in_cruise = (~in_accel) & (t_rel <= cruise_end)
    in_decel = ~(in_accel | in_cruise)
    dist[in_accel] = (
        0.5 * t_rel[in_accel] ** 2 / t_accel_s if t_accel_s > 0 else 0.0
    )
    dist[in_cruise] = 0.5 * t_accel_s + (t_rel[in_cruise] - t_accel_s)
    t2 = t_rel[in_decel] - t_accel_s - duration_s
    dist[in_decel] = (
        0.5 * t_accel_s + duration_s + t2 - 0.5 * t2 * t2 / t_accel_s
    )
    return dist / total_dist


def reconstruct_segments(rows, run_index, block_start_s, block_end_s, block_name=None):
    """List of segments describing the reconstructed commanded trajectory
    across one block, in ideal_position_rev units (revolutions), each
    dict either:
      {'kind': 'hold', 't0', 't1', 'value_start', 'value_end'} (flat,
        value_start == value_end, exact for software-paced/instant moves)
      {'kind': 'trapezoid', 't0', 't1', 'value_start', 'value_end',
        't_accel_s', 'duration_s'} (continuous ramp, real-time-scaled --
        t_accel_s=0 degenerates to a pure constant-velocity linear ramp,
        used for fixed-rate burst moves such as the C block's approach/
        return, which have no accel phase at all)
    Times are seconds since block_start_s. `block_name` enables the
    fixed-rate burst-move reconstruction for blocks in
    BURST_RATE_BY_BLOCK_FS_S; omit it (or pass an unlisted name) to fall
    back to a zero-order hold for moves with no rate_full_steps_s, as
    before."""
    moves = [
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['ideal_position_rev']
        and block_start_s <= r['ids_time_s'] <= block_end_s
    ]
    block_len = block_end_s - block_start_s
    fallback_rate = BURST_RATE_BY_BLOCK_FS_S.get(block_name)
    segments = []
    prev_t = 0.0
    prev_val = 0.0
    for index, row in enumerate(moves):
        row_t = row['ids_time_s'] - block_start_s
        row_val = float(row['ideal_position_rev'])
        rate_str = row.get('rate_full_steps_s', '')
        if row_t > prev_t:
            segments.append({
                'kind': 'hold', 't0': prev_t, 't1': row_t,
                'value_start': prev_val, 'value_end': prev_val,
            })
        next_t = (
            moves[index + 1]['ids_time_s'] - block_start_s
            if index + 1 < len(moves) else block_len
        )
        if is_controller_paced(rate_str):
            t_accel_theory, duration_theory, t_total_theory = trapezoid_profile(
                float(rate_str)
            )
            real_avail = max(next_t - row_t, 1.0e-6)
            scale = real_avail / t_total_theory if t_total_theory > 0 else 1.0
            segments.append({
                'kind': 'trapezoid', 't0': row_t, 't1': row_t + real_avail,
                'value_start': prev_val, 'value_end': row_val,
                't_accel_s': t_accel_theory * scale,
                'duration_s': duration_theory * scale,
            })
            prev_t, prev_val = row_t + real_avail, row_val
        elif not rate_str and fallback_rate and row_val != prev_val:
            # Unlike the D-rate case above, do NOT scale this to the gap
            # until the next logged event: burst moves here are routinely
            # followed by a deliberate dwell() of many seconds, so that gap
            # is "move time + dwell time", not just "how long this move
            # took". Use the theoretical constant-velocity duration as-is;
            # the bridging 'hold' segment built on the next loop iteration
            # (or the block-end hold below) correctly absorbs any dwell
            # that follows.
            delta_full_steps = abs(row_val - prev_val) * MOTOR_FULL_STEPS_PER_REV
            duration_theory = delta_full_steps / fallback_rate
            duration_theory = min(duration_theory, max(next_t - row_t, 1.0e-6))
            segments.append({
                'kind': 'trapezoid', 't0': row_t, 't1': row_t + duration_theory,
                'value_start': prev_val, 'value_end': row_val,
                't_accel_s': 0.0,
                'duration_s': duration_theory,
            })
            prev_t, prev_val = row_t + duration_theory, row_val
        else:
            segments.append({
                'kind': 'hold', 't0': row_t, 't1': row_t,
                'value_start': row_val, 'value_end': row_val,
            })
            prev_t, prev_val = row_t, row_val
    if block_len > prev_t:
        segments.append({
            'kind': 'hold', 't0': prev_t, 't1': block_len,
            'value_start': prev_val, 'value_end': prev_val,
        })
    return segments


def sample_segments(segments, scale=1.0, n_trapezoid=80, truncate=False):
    """Flatten segments into (times, values) point arrays for step/line
    plotting, with `value`s multiplied by `scale` (e.g. LEAD_PITCH_M*1e6 for
    µm, or 2*pi for radians). truncate=True drops the final block-end hold
    tail (see command_series's docstring in the plotting scripts)."""
    times, values = [], []
    usable = segments
    if truncate and usable and usable[-1]['kind'] == 'hold' and usable[-1]['t1'] > usable[-1]['t0']:
        usable = usable[:-1]
    for seg in usable:
        if seg['kind'] == 'hold':
            times.append(seg['t0'])
            values.append(seg['value_start'] * scale)
            if seg['t1'] > seg['t0']:
                times.append(seg['t1'])
                values.append(seg['value_end'] * scale)
        else:
            t0, t1 = seg['t0'], seg['t1']
            ts = np.linspace(t0, t1, n_trapezoid)
            fracs = np.array([
                trapezoid_fraction(t - t0, seg['t_accel_s'], seg['duration_s'])
                for t in ts
            ])
            vals = (
                seg['value_start']
                + (seg['value_end'] - seg['value_start']) * fracs
            ) * scale
            times.extend(ts.tolist())
            values.extend(vals.tolist())
    return np.asarray(times), np.asarray(values)


def commanded_value_at(segments, t_query, scale=1.0):
    """Evaluate the reconstructed commanded trajectory at arbitrary query
    times (e.g. a measured or simulated trace's own time grid), so it can be
    subtracted pointwise for a tracking-error series. `t_query` is an
    ndarray of seconds since block start, `scale` is the same unit
    multiplier as sample_segments (e.g. LEAD_PITCH_M*1e6 for µm)."""
    t_query = np.asarray(t_query, dtype=float)
    values = np.empty_like(t_query)
    block_end = segments[-1]['t1'] if segments else 0.0
    for seg in segments:
        if seg is segments[-1]:
            mask = (t_query >= seg['t0']) & (t_query <= seg['t1'])
        else:
            mask = (t_query >= seg['t0']) & (t_query < seg['t1'])
        if not np.any(mask):
            continue
        if seg['kind'] == 'hold':
            values[mask] = seg['value_start']
        else:
            frac = trapezoid_fraction_array(
                t_query[mask] - seg['t0'], seg['t_accel_s'], seg['duration_s']
            )
            values[mask] = (
                seg['value_start'] + (seg['value_end'] - seg['value_start']) * frac
            )
    values[t_query < (segments[0]['t0'] if segments else 0.0)] = (
        segments[0]['value_start'] if segments else 0.0
    )
    values[t_query > block_end] = segments[-1]['value_end'] if segments else 0.0
    return values * scale
