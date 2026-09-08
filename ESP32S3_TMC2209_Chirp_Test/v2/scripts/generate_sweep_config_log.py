#!/usr/bin/env python3
'''Regenerate ../SWEEP_CONFIG_LOG.md from chirp_v2_schedule.py.

This is the reference to hand to a processing script: exact amplitude,
notch, clamp, and timing parameters for the current design, without having
to re-derive them from the README or read the module source. It is
generated, not hand-edited -- run this after any constant in
chirp_v2_schedule.py changes, so the log can never silently drift from the
code that actually produced a given run's commanded signal.

Run: `python generate_sweep_config_log.py`
'''

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import chirp_v2_schedule as sched  # noqa: E402


def timing_breakdown() -> list[tuple[str, float, float]]:
    '''(label, start_s, end_s) for every segment in one amplitude level.'''
    t = 0.0
    rows: list[tuple[str, float, float]] = []

    def add(label: str, duration: float) -> None:
        nonlocal t
        rows.append((label, t, t + duration))
        t += duration

    add('idle pre-roll', sched.PRE_ROLL_S)
    add('sync marker 1', sched.SYNC_MARKER_DURATION_S)
    add('up-sweep: log segment', sched.LOG_DURATION_S)
    add('up-sweep: linear segment', sched.LINEAR_DURATION_S)
    add('mid-dwell', sched.MID_DWELL_S)
    add('sync marker 2', sched.SYNC_MARKER_DURATION_S)
    add('down-sweep: linear segment', sched.LINEAR_DURATION_S)
    add('down-sweep: log segment', sched.LOG_DURATION_S)
    add('tail (idle)', sched.TAIL_S)
    return rows


def resonance_crossings(rows: list[tuple[str, float, float]]) -> tuple[float, float]:
    up_lin_start = next(t0 for label, t0, _ in rows if label == 'up-sweep: linear segment')
    down_lin_start = next(t0 for label, t0, _ in rows if label == 'down-sweep: linear segment')
    rate = (sched.LINEAR_END_HZ - sched.LINEAR_START_HZ) / sched.LINEAR_DURATION_S
    t_up = up_lin_start + (sched.F_N_HZ_PLACEHOLDER - sched.LINEAR_START_HZ) / rate
    t_down = down_lin_start + (sched.LINEAR_END_HZ - sched.F_N_HZ_PLACEHOLDER) / rate
    return t_up, t_down


def microstep_floor_summary() -> list[tuple[float, float, float]]:
    '''(cruise_full_steps, min_microsteps, frequency_of_min_hz) per level.'''
    f_grid = np.logspace(0, 3, 4000)
    out = []
    for cruise in sched.CRUISE_LEVELS_FULL_STEPS:
        microsteps = sched.amplitude_profile(f_grid, cruise) * sched.MRES
        idx = int(microsteps.argmin())
        out.append((cruise, float(microsteps[idx]), float(f_grid[idx])))
    return out


def render_markdown() -> str:
    rows = timing_breakdown()
    t_res_up, t_res_down = resonance_crossings(rows)
    floor_rows = microstep_floor_summary()
    notch_low = sched.F_N_HZ_PLACEHOLDER - sched.NOTCH_HALFWIDTH_HZ
    notch_high = sched.F_N_HZ_PLACEHOLDER + sched.NOTCH_HALFWIDTH_HZ
    generated = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')

    lines: list[str] = []
    w = lines.append

    w('# v2 chirp -- sweep configuration log')
    w('')
    w(f'Generated {generated} by `generate_sweep_config_log.py` from '
      '`chirp_v2_schedule.py`. **Do not hand-edit** -- change the schedule '
      'module and regenerate instead, or this stops being trustworthy as a '
      'processing reference.')
    w('')
    w('## Mechanics / drive')
    w('')
    w('| Parameter | Value |')
    w('|---|---:|')
    w(f'| Lead | {sched.LEAD_MM_PER_REV:g} mm/rev |')
    w(f'| Full steps/rev | {sched.FULL_STEPS_PER_REV:g} |')
    w(f'| Full step size | {sched.FULL_STEP_MM * 1000:g} um ({sched.FULL_STEP_MM:g} mm) |')
    w(f'| MRES | {sched.MRES:g} (1/{sched.MRES:g} step) |')
    w(f'| 1 microstep | {sched.ONE_MICROSTEP_FULL_STEPS:.6g} full steps '
      f'({sched.FULL_STEP_MM * 1000 / sched.MRES:.4g} um) |')
    w(f'| Run current | {sched.CURRENT_RMS_MA:g} mA RMS |')
    w('')
    w('## Plant placeholders -- REPLACE with same-day ringdown measurements')
    w('')
    w('| Parameter | Placeholder value | Source |')
    w('|---|---:|---|')
    w(f'| f_n | {sched.F_N_HZ_PLACEHOLDER:g} Hz | Prior analytical model; not measured |')
    w(f'| Q | {sched.Q_PLACEHOLDER:g} | Prior analytical model; not measured |')
    w(f'| zeta | {sched.ZETA_PLACEHOLDER:.4g} | = 1/(2*Q) |')
    w('')
    w('## Amplitude profile')
    w('')
    w('`A_cmd(f) = min(cruise * notch(f), stroke_clamp, step_rate_clamp / '
      '(2*pi*f*MRES))`, `notch(f)` a raised-cosine window centered on f_n.')
    w('')
    w('| Parameter | Value |')
    w('|---|---:|')
    w('| Cruise levels | ' + ' / '.join(f'{c:g}' for c in sched.CRUISE_LEVELS_FULL_STEPS) + ' full steps |')
    w(f'| Notch half-width | {sched.NOTCH_HALFWIDTH_HZ:.4g} Hz |')
    w(f'| Notch band (at f_n placeholder) | {notch_low:.1f} - {notch_high:.1f} Hz |')
    w(f'| Notch depth ratio | {sched.NOTCH_DEPTH_RATIO:g} (amplitude at f_n = {sched.NOTCH_DEPTH_RATIO * 100:g}% of cruise) |')
    w(f'| Stroke clamp | {sched.STROKE_CLAMP_MM:g} mm ({sched.STROKE_CLAMP_MM / sched.FULL_STEP_MM:g} full steps) |')
    w(f'| Step-rate clamp | {sched.STEP_RATE_CLAMP_HZ:,.0f} Hz |')
    w('| Velocity/presliding clamp | not applied by default (see README Section 3) |')
    w('')
    w('Minimum amplitude actually reached per cruise level (from '
      '`amplitude_profile`, swept 1-1000 Hz):')
    w('')
    w('| Cruise (full steps) | Minimum amplitude | At frequency | Above 1 microstep? |')
    w('|---:|---:|---:|:---:|')
    for cruise, min_microsteps, f_at_min in floor_rows:
        ok = 'yes' if min_microsteps >= 1.0 else '**NO**'
        w(f'| {cruise:g} | {min_microsteps:.3g} microsteps | {f_at_min:.1f} Hz | {ok} |')
    w('')
    w('Stribeck-velocity reference (not applied as a default clamp -- see '
      'README Section 3):')
    w('')
    w(f'- Nut: {sched.STRIBECK_VELOCITY_NUT_MM_S:g} mm/s')
    w(f'- Guideway: {sched.STRIBECK_VELOCITY_GUIDEWAY_MM_S:g} mm/s')
    w('')
    w('## Sweep frequency law')
    w('')
    w('| Segment | Range | Duration |')
    w('|---|---|---:|')
    w(f'| Log | {sched.LOG_START_HZ:g} -> {sched.LOG_END_HZ:g} Hz | {sched.LOG_DURATION_S:g} s |')
    w(f'| Linear | {sched.LINEAR_START_HZ:g} -> {sched.LINEAR_END_HZ:g} Hz | '
      f'{sched.LINEAR_DURATION_S:g} s ({(sched.LINEAR_END_HZ - sched.LINEAR_START_HZ) / sched.LINEAR_DURATION_S:.4g} Hz/s) |')
    w('')
    w('Down-sweep is the exact mirror (linear 1000->60 Hz, then log 60->1 Hz).')
    w('')
    w('## Sync marker')
    w('')
    w('| Parameter | Value |')
    w('|---|---:|')
    w(f'| Frequency | {sched.SYNC_MARKER_FREQ_HZ:g} Hz |')
    w(f'| Amplitude | {sched.SYNC_MARKER_AMPLITUDE_FULL_STEPS:g} full steps |')
    w(f'| Burst count | {sched.SYNC_BURST_COUNT:g} |')
    w(f'| Burst duration | {sched.SYNC_BURST_S:g} s |')
    w(f'| Gap duration | {sched.SYNC_GAP_S:g} s |')
    w(f'| Total marker duration | {sched.SYNC_MARKER_DURATION_S:g} s |')
    w('')
    w('## Level timing breakdown (one amplitude level, one direction pass)')
    w('')
    w('| Segment | Start (s) | End (s) | Duration (s) |')
    w('|---|---:|---:|---:|')
    for label, t0, t1 in rows:
        w(f'| {label} | {t0:.2f} | {t1:.2f} | {t1 - t0:.2f} |')
    w(f'| **Total** | | | **{sched.LEVEL_DURATION_S:.2f} s ({sched.LEVEL_DURATION_S / 60.0:.2f} min)** |')
    w('')
    w(f'Resonance crossing (f_n placeholder = {sched.F_N_HZ_PLACEHOLDER:g} Hz): '
      f'up-sweep at t={t_res_up:.2f} s, down-sweep at t={t_res_down:.2f} s '
      '(both relative to level start).')
    w('')
    w('## Full campaign (all cruise levels)')
    w('')
    n_levels = len(sched.CRUISE_LEVELS_FULL_STEPS)
    level_word = 'level' if n_levels == 1 else 'levels'
    w(f'{n_levels} cruise {level_word} x {sched.LEVEL_DURATION_S:.1f} s = '
      f'{n_levels * sched.LEVEL_DURATION_S / 60.0:.1f} min of sweeping, plus '
      'ringdown characterization time (Section 1) not counted here.')
    w('')
    return '\n'.join(lines) + '\n'


def main() -> None:
    out_path = HERE.parent / 'SWEEP_CONFIG_LOG.md'
    out_path.write_text(render_markdown(), encoding='utf-8')
    print(f'Wrote: {out_path}')


if __name__ == '__main__':
    main()
