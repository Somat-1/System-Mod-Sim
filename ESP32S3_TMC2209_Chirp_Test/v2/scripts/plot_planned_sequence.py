#!/usr/bin/env python3
'''Render the planned v2 chirp excitation before any hardware exists.

The full per-level sequence is ~9 minutes and the interesting content near
resonance lasts a couple of seconds, so plotting the whole thing at once
would either hide the fine structure or be unreadably wide. Instead this
renders four views:

1. The peak-displacement envelope across the *entire* ~9 minute level: idle,
   both sync markers, both sweep directions, the mid-dwell, and the tail --
   traced (not the dense oscillating carrier itself) for all three cruise
   levels, against a marked 1-microstep floor (see MRES note below).
2. A ~24 s window at the start of a run: idle pre-roll tail, sync marker,
   and the first slow cycles of the log sweep.
3. A ~3 s window mid-run, centered on the placeholder resonance crossing
   inside the linear segment, showing the smooth notch.
4. A ~7.5 s window at the turnaround: the last second of the up-sweep at
   1000 Hz, the 5 s mid-dwell, the second sync marker, and the first second
   of the down-sweep starting back down from 1000 Hz.

MRES is 16 (per chirp_v2_schedule.MRES), not native 1/256 -- see the module
docstring there for why. 1 microstep is 1/16 = 0.0625 full steps; the
resonance notch bottoms out at 10% of cruise, which stays a few microsteps
even at the lowest cruise level, so nothing here asks for less than 1
microstep (unlike the revision-1 e_max design, which did above ~816 Hz).

Run: `python plot_planned_sequence.py`
'''

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chirp_v2_schedule import (  # noqa: E402
    CRUISE_LEVELS_FULL_STEPS,
    F_N_HZ_PLACEHOLDER,
    LEVEL_DURATION_S,
    LINEAR_DURATION_S,
    LINEAR_END_HZ,
    LINEAR_START_HZ,
    LOG_DURATION_S,
    MID_DWELL_S,
    MRES,
    NOTCH_DEPTH_RATIO,
    NOTCH_HALFWIDTH_HZ,
    ONE_MICROSTEP_FULL_STEPS,
    PRE_ROLL_S,
    SWEEP_DURATION_S,
    SYNC_BURST_COUNT,
    SYNC_BURST_S,
    SYNC_GAP_S,
    SYNC_MARKER_AMPLITUDE_FULL_STEPS,
    SYNC_MARKER_DURATION_S,
    TAIL_S,
    amplitude_profile,
    down_sweep_frequency,
    render_down_sweep_command,
    render_up_sweep_command,
    sync_marker_command,
    up_sweep_frequency,
)

SWEEP_START_ABS_S = PRE_ROLL_S + SYNC_MARKER_DURATION_S
UP_SWEEP_END_ABS_S = SWEEP_START_ABS_S + SWEEP_DURATION_S
MARKER2_START_ABS_S = UP_SWEEP_END_ABS_S + MID_DWELL_S
DOWN_SWEEP_START_ABS_S = MARKER2_START_ABS_S + SYNC_MARKER_DURATION_S
DOWN_SWEEP_END_ABS_S = DOWN_SWEEP_START_ABS_S + SWEEP_DURATION_S

GRID_DT_S = 1.0 / 4000.0
GRID_END_ABS_S = 170.0

START_WINDOW_ABS_S = (PRE_ROLL_S - 3.0, SWEEP_START_ABS_S + 20.5)

_RESONANCE_LINEAR_FRACTION_S = (F_N_HZ_PLACEHOLDER - LINEAR_START_HZ) / (
    (LINEAR_END_HZ - LINEAR_START_HZ) / LINEAR_DURATION_S
)
RESONANCE_ABS_S = (
    SWEEP_START_ABS_S + LOG_DURATION_S + _RESONANCE_LINEAR_FRACTION_S
)
MIDRUN_WINDOW_ABS_S = (RESONANCE_ABS_S - 1.5, RESONANCE_ABS_S + 1.5)

_RESONANCE_LINEAR_FRACTION_DOWN_S = (LINEAR_END_HZ - F_N_HZ_PLACEHOLDER) / (
    (LINEAR_END_HZ - LINEAR_START_HZ) / LINEAR_DURATION_S
)
RESONANCE_DOWN_ABS_S = (
    DOWN_SWEEP_START_ABS_S + _RESONANCE_LINEAR_FRACTION_DOWN_S
)

TURNAROUND_TAIL_S = 1.0  # How much of each sweep's edge to show.
TURNAROUND_WINDOW_ABS_S = (
    UP_SWEEP_END_ABS_S - TURNAROUND_TAIL_S,
    DOWN_SWEEP_START_ABS_S + TURNAROUND_TAIL_S,
)

MID_LEVEL_INDEX = 0  # Index into CRUISE_LEVELS_FULL_STEPS for the zoom panels.


def build_level_timeline() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Idle -> sync marker -> start of the up-sweep, at MID_LEVEL_INDEX.

    Returns (t_abs_s, position_full_steps, frequency_hz). frequency_hz is
    NaN outside the sweep, where there is no defined instantaneous tone.
    '''
    t_abs = np.arange(0.0, GRID_END_ABS_S, GRID_DT_S)
    position = np.zeros_like(t_abs)
    frequency = np.full_like(t_abs, np.nan)

    marker_mask = (t_abs >= PRE_ROLL_S) & (t_abs < SWEEP_START_ABS_S)
    position[marker_mask] = sync_marker_command(t_abs[marker_mask] - PRE_ROLL_S)

    sweep_mask = t_abs >= SWEEP_START_ABS_S
    local_t = t_abs[sweep_mask] - SWEEP_START_ABS_S
    cruise = CRUISE_LEVELS_FULL_STEPS[MID_LEVEL_INDEX]
    f_hz, _amplitude, sweep_position = render_up_sweep_command(
        local_t, cruise_full_steps=cruise,
    )
    position[sweep_mask] = sweep_position
    frequency[sweep_mask] = f_hz
    return t_abs, position, frequency


ENVELOPE_DT_S = 0.02  # Coarse: an envelope only needs the instantaneous
                      # frequency, not enough samples to resolve the carrier.


def build_full_run_envelope(
    cruise_full_steps: float,
) -> tuple[np.ndarray, np.ndarray]:
    '''Peak commanded displacement (full steps) across the whole ~9 min level.

    No carrier is rendered here (that is what the three zoom panels are
    for) -- just A_cmd(f(t)) traced through idle/marker/sweep/dwell/tail,
    including both sync-marker bursts at their own fixed amplitude.
    '''
    t_abs = np.arange(0.0, LEVEL_DURATION_S, ENVELOPE_DT_S)
    envelope = np.zeros_like(t_abs)

    period = SYNC_BURST_S + SYNC_GAP_S
    for marker_start in (PRE_ROLL_S, MARKER2_START_ABS_S):
        local = t_abs - marker_start
        for burst_index in range(SYNC_BURST_COUNT):
            burst_lo = burst_index * period
            burst_hi = burst_lo + SYNC_BURST_S
            in_burst = (local >= burst_lo) & (local < burst_hi)
            envelope[in_burst] = SYNC_MARKER_AMPLITUDE_FULL_STEPS

    up_mask = (t_abs >= SWEEP_START_ABS_S) & (t_abs < UP_SWEEP_END_ABS_S)
    up_f_hz = up_sweep_frequency(t_abs[up_mask] - SWEEP_START_ABS_S)
    envelope[up_mask] = amplitude_profile(up_f_hz, cruise_full_steps)

    down_mask = (
        (t_abs >= DOWN_SWEEP_START_ABS_S) & (t_abs < DOWN_SWEEP_END_ABS_S)
    )
    down_f_hz = down_sweep_frequency(t_abs[down_mask] - DOWN_SWEEP_START_ABS_S)
    envelope[down_mask] = amplitude_profile(down_f_hz, cruise_full_steps)

    return t_abs, envelope


def build_turnaround_timeline(
    cruise_full_steps: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Last second of the up-sweep -> mid-dwell -> sync marker -> first second of the down-sweep.

    Returns (t_abs_s, position_full_steps, frequency_hz), frequency_hz NaN
    outside the two sweep slivers. The up-sweep tail is computed from its
    own full 240 s history so its phase (and therefore instantaneous
    amplitude, from the notch taper) is correct at the edge; the down-sweep
    only needs its own first second, integrated fresh from its own start.
    '''
    up_local_full = np.arange(0.0, SWEEP_DURATION_S, GRID_DT_S)
    up_f_hz, _up_amplitude, up_position = render_up_sweep_command(
        up_local_full, cruise_full_steps=cruise_full_steps,
    )
    up_tail_mask = up_local_full >= SWEEP_DURATION_S - TURNAROUND_TAIL_S
    t_before = UP_SWEEP_END_ABS_S + (
        up_local_full[up_tail_mask] - SWEEP_DURATION_S
    )
    position_before = up_position[up_tail_mask]
    frequency_before = up_f_hz[up_tail_mask]

    gap_local = np.arange(0.0, MID_DWELL_S + SYNC_MARKER_DURATION_S, GRID_DT_S)
    t_gap = UP_SWEEP_END_ABS_S + gap_local
    position_gap = np.zeros_like(gap_local)
    marker_mask = gap_local >= MID_DWELL_S
    position_gap[marker_mask] = sync_marker_command(
        gap_local[marker_mask] - MID_DWELL_S
    )
    frequency_gap = np.full_like(gap_local, np.nan)

    down_local = np.arange(0.0, TURNAROUND_TAIL_S, GRID_DT_S)
    down_f_hz, _down_amplitude, down_position = render_down_sweep_command(
        down_local, cruise_full_steps=cruise_full_steps,
    )
    t_after = DOWN_SWEEP_START_ABS_S + down_local

    t_abs = np.concatenate([t_before, t_gap, t_after])
    position = np.concatenate([position_before, position_gap, down_position])
    frequency = np.concatenate([frequency_before, frequency_gap, down_f_hz])
    return t_abs, position, frequency


def window_slice(
    t_abs: np.ndarray, window: tuple[float, float]
) -> np.ndarray:
    return (t_abs >= window[0]) & (t_abs <= window[1])


def plot_time_window(
    axis: plt.Axes,
    t_abs: np.ndarray,
    position: np.ndarray,
    window: tuple[float, float],
    *,
    envelope: tuple[np.ndarray, np.ndarray] | None = None,
) -> None:
    mask = window_slice(t_abs, window)
    axis.plot(t_abs[mask], position[mask], color='#2474b5', lw=0.8)
    if envelope is not None:
        env_t, env_a = envelope
        axis.plot(env_t, env_a, color='#c0472c', lw=1.1, ls='--',
                   label='+/-A_cmd(f(t)) envelope')
        axis.plot(env_t, -env_a, color='#c0472c', lw=1.1, ls='--')
    axis.axhline(0.0, color='#333333', lw=0.5, alpha=0.6)
    axis.set_xlim(*window)
    axis.set_ylabel('Commanded position (full steps)')
    axis.grid(True, alpha=0.25)


def shade_pre_sweep_context(axis: plt.Axes) -> None:
    axis.axvspan(START_WINDOW_ABS_S[0], PRE_ROLL_S, color='#c7ccd1',
                 alpha=0.5, lw=0, label='idle pre-roll (tail)')
    axis.axvspan(PRE_ROLL_S, SWEEP_START_ABS_S, color='#8e5bb7',
                 alpha=0.35, lw=0, label='sync marker (3x0.1 s bursts)')
    axis.axvspan(SWEEP_START_ABS_S, START_WINDOW_ABS_S[1], color='#e8eef4',
                 alpha=0.5, lw=0, label='log sweep, 1 Hz start')


def shade_turnaround_context(axis: plt.Axes) -> None:
    axis.axvspan(TURNAROUND_WINDOW_ABS_S[0], UP_SWEEP_END_ABS_S,
                 color='#e8eef4', alpha=0.5, lw=0, label='up-sweep, 1000 Hz end')
    axis.axvspan(UP_SWEEP_END_ABS_S, MARKER2_START_ABS_S, color='#c7ccd1',
                 alpha=0.5, lw=0, label=f'{MID_DWELL_S:g} s mid-dwell')
    axis.axvspan(MARKER2_START_ABS_S, DOWN_SWEEP_START_ABS_S, color='#8e5bb7',
                 alpha=0.35, lw=0, label='sync marker (3x0.1 s bursts)')
    axis.axvspan(DOWN_SWEEP_START_ABS_S, TURNAROUND_WINDOW_ABS_S[1],
                 color='#fbeee6', alpha=0.6, lw=0,
                 label='down-sweep, 1000 Hz start')


def shade_full_run_context(axis: plt.Axes) -> None:
    bounds = [
        (0.0, PRE_ROLL_S, '#c7ccd1', 'idle (pre-roll/dwell/tail)'),
        (PRE_ROLL_S, SWEEP_START_ABS_S, '#8e5bb7', 'sync marker'),
        (SWEEP_START_ABS_S, UP_SWEEP_END_ABS_S, '#e8eef4', 'up-sweep'),
        (UP_SWEEP_END_ABS_S, MARKER2_START_ABS_S, '#c7ccd1', None),
        (MARKER2_START_ABS_S, DOWN_SWEEP_START_ABS_S, '#8e5bb7', None),
        (DOWN_SWEEP_START_ABS_S, DOWN_SWEEP_END_ABS_S, '#fbeee6', 'down-sweep'),
        (DOWN_SWEEP_END_ABS_S, LEVEL_DURATION_S, '#c7ccd1', None),
    ]
    for t0, t1, color, label in bounds:
        axis.axvspan(t0 / 60.0, t1 / 60.0, color=color, alpha=0.5, lw=0,
                     label=label)


_LEVEL_STYLES = (
    ('#0b3d68', 2.2, '-'),
    ('#8fb7d9', 1.0, ':'),
    ('#1f8f6b', 1.6, '--'),
    ('#c0472c', 1.6, '-.'),
)


def plot_full_run_envelope(axis: plt.Axes) -> None:
    shade_full_run_context(axis)
    for cruise, (color, lw, ls) in zip(CRUISE_LEVELS_FULL_STEPS, _LEVEL_STYLES):
        t_abs, envelope = build_full_run_envelope(cruise)
        axis.plot(t_abs / 60.0, envelope, ls, color=color, lw=lw,
                   label=f'cruise = {cruise:g} full steps')

    axis.axhline(ONE_MICROSTEP_FULL_STEPS, color='#c0472c', lw=1.3, ls='--')
    axis.annotate(
        f'1 microstep floor (1/{MRES:g} full step)',
        xy=(LEVEL_DURATION_S / 60.0 - 0.1, ONE_MICROSTEP_FULL_STEPS),
        xytext=(0, 6), textcoords='offset points', ha='right',
        fontsize=7.5, color='#c0472c',
    )
    for t_res in (RESONANCE_ABS_S, RESONANCE_DOWN_ABS_S):
        axis.axvline(t_res / 60.0, color='#555555', lw=0.7, ls=':')

    axis.set_yscale('log')
    axis.set_xlim(0.0, LEVEL_DURATION_S / 60.0)
    axis.set_ylabel('Peak commanded displacement (full steps)')
    axis.set_xlabel('Time since level start (min)')
    axis.set_title(
        'Full run envelope (cruise = '
        f'{" / ".join(f"{c:g}" for c in CRUISE_LEVELS_FULL_STEPS)} full steps): '
        f'{LEVEL_DURATION_S / 60.0:.1f} min total, both sweep directions, '
        'carrier not shown -- see the zoom panels below for that'
    )
    axis.grid(True, which='both', alpha=0.2)
    axis.legend(loc='lower left', fontsize=7.5, frameon=False, ncols=1,
                labelspacing=0.3)

    microstep_axis = axis.secondary_yaxis(
        'right', functions=(lambda a: a * MRES, lambda m: m / MRES)
    )
    microstep_axis.set_ylabel('Peak commanded displacement (microsteps)')


def main() -> None:
    t_abs, position, frequency = build_level_timeline()
    cruise_mid = CRUISE_LEVELS_FULL_STEPS[MID_LEVEL_INDEX]
    turn_t_abs, turn_position, turn_frequency = build_turnaround_timeline(
        cruise_mid
    )

    fig = plt.figure(figsize=(13, 12.5), constrained_layout=True)
    grid = fig.add_gridspec(4, 1, height_ratios=(3.0, 2.0, 2.0, 2.0))
    overview_ax = fig.add_subplot(grid[0, 0])
    start_ax = fig.add_subplot(grid[1, 0])
    midrun_ax = fig.add_subplot(grid[2, 0])
    turnaround_ax = fig.add_subplot(grid[3, 0])

    plot_full_run_envelope(overview_ax)

    shade_pre_sweep_context(start_ax)
    plot_time_window(start_ax, t_abs, position, START_WINDOW_ABS_S)
    start_ax.set_xlabel('Time since level start (s)')
    start_ax.set_title(
        f'Beginning of a run (cruise = {cruise_mid:g} full steps): idle '
        'tail, sync marker, first ~20 s of the 1->60 Hz log sweep'
    )
    start_ax.legend(loc='upper left', fontsize=7.5, frameon=False, ncols=3)

    mid_mask = window_slice(t_abs, MIDRUN_WINDOW_ABS_S)
    envelope = (
        t_abs[mid_mask],
        amplitude_profile(frequency[mid_mask], cruise_mid),
    )
    plot_time_window(midrun_ax, t_abs, position, MIDRUN_WINDOW_ABS_S,
                      envelope=envelope)
    midrun_ax.axvline(RESONANCE_ABS_S, color='#555555', lw=0.8, ls=':')
    midrun_ax.set_xlabel('Time since level start (s)')
    midrun_ax.set_title(
        'Mid-run, linear segment crossing the placeholder resonance '
        f'({F_N_HZ_PLACEHOLDER:g} Hz at t={RESONANCE_ABS_S:.1f} s): smooth '
        f'notch, {NOTCH_HALFWIDTH_HZ:.0f} Hz half-width, bottoms out at '
        f'{NOTCH_DEPTH_RATIO * 100:.0f}% of cruise'
    )
    midrun_ax.legend(loc='upper right', fontsize=7.5, frameon=False)

    shade_turnaround_context(turnaround_ax)
    turn_envelope_amplitude = np.full_like(turn_frequency, np.nan)
    turn_has_tone = ~np.isnan(turn_frequency)
    turn_envelope_amplitude[turn_has_tone] = amplitude_profile(
        turn_frequency[turn_has_tone], cruise_mid,
    )
    plot_time_window(turnaround_ax, turn_t_abs, turn_position,
                      TURNAROUND_WINDOW_ABS_S,
                      envelope=(turn_t_abs, turn_envelope_amplitude))
    turnaround_ax.axvline(UP_SWEEP_END_ABS_S, color='#555555', lw=0.8, ls=':')
    turnaround_ax.axvline(DOWN_SWEEP_START_ABS_S, color='#555555', lw=0.8,
                           ls=':')
    turnaround_ax.set_xlabel('Time since level start (s)')
    turnaround_ax.set_title(
        f'Turnaround at {LINEAR_END_HZ:g} Hz (cruise = {cruise_mid:g} full '
        f'steps): {TURNAROUND_TAIL_S:g} s of up-sweep, {MID_DWELL_S:g} s '
        'dwell, sync marker, then the down-sweep starts back down from '
        f'{LINEAR_END_HZ:g} Hz'
    )
    turnaround_ax.legend(loc='upper left', fontsize=7.5, frameon=False, ncols=2)

    fig.suptitle(
        'v2 chirp excitation (design preview, pre-hardware) - hard notch '
        'replaced by a smooth taper, cruise amplitude chosen directly',
        fontsize=13,
    )

    out_dir = HERE.parent / 'rendered_assets'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'planned_chirp_v2_excitation.png'
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    print(f'Saved: {out_path}')
    print(
        f'Placeholder resonance crossing inside the up-sweep: '
        f'{RESONANCE_ABS_S:.2f} s since level start '
        f'({RESONANCE_ABS_S - SWEEP_START_ABS_S:.2f} s into the sweep).'
    )

    f_grid = np.logspace(0, 3, 4000)
    for cruise in CRUISE_LEVELS_FULL_STEPS:
        microsteps = amplitude_profile(f_grid, cruise) * MRES
        print(
            f'cruise={cruise:g} full steps: minimum amplitude anywhere in '
            f'the sweep is {microsteps.min():.2f} microsteps '
            f'(at f={f_grid[microsteps.argmin()]:.1f} Hz)'
        )


if __name__ == '__main__':
    main()
