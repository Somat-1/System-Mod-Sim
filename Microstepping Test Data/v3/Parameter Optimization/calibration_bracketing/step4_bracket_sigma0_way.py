#!/usr/bin/env python3
"""Step 4: bracket sigma0_way and Fs_way from two independent measured
signatures, before any optimizer sees the data:

(a) D_0.125 at MRES 1/1 (runs 5, 6): at mres=1 each MOVE_ACK is one full
    10 um step, and pulses are ~8 s apart (pulse_rate = 0.125*1 Hz), so
    each step's transient is fully isolated in time -- ideal for reading
    off how much of the 10 um happens as a slow compliant ramp
    (presliding, ~Fs/sigma0) versus a sudden jump.
(b) The C block's 40 um approach move (fast burst, not isolated the same
    way, but the pre-slip compliant buildup at the very start of the move
    is directly visible at 1 ms sample resolution).

If presliding really cost Fs/sigma0 = 4.5/7.6e5 = 5.9 um (the Rev 3
carry-over placeholder), a 10 um full step should show roughly half its
travel as a slow ramp before release. This script measures how much it
actually shows.
"""
import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3\scripts")
sys.path.insert(0, str(SCRIPTS_DIR))
from plot_block_montage import parse_ids, load_log_rows, find_block, IDS_PATH, LOG_PATH

OUT_DIR = Path(__file__).resolve().parent


def load_common():
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)
    return time_s, position_nm, sample_period_s, rows


def measure_d0125_treads(run_index, time_s, position_nm, sample_period_s, rows):
    """For each full-step MOVE_ACK in D_0.125 (run_index at MRES 1/1),
    measure the pre-jump creep distance and the jump itself."""
    start_s, end_s = find_block(rows, run_index, 'D_0.125')
    moves = [
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['ideal_position_rev']
        and start_s <= float(r['ids_time_s']) <= end_s
    ]
    results = []
    for move in moves:
        move_abs_t = float(move['ids_time_s'])
        move_t = move_abs_t - start_s
        idx0 = int(round((move_abs_t - 2.0) / sample_period_s))
        idx1 = int(round((move_abs_t + 2.0) / sample_period_s))
        idx0, idx1 = max(0, idx0), min(position_nm.size, idx1)
        window = position_nm[idx0:idx1] / 1000.0  # um
        t_window = time_s[idx0:idx1] - move_abs_t
        if window.size < 10:
            continue
        baseline = np.median(window[: max(1, int(0.5 / sample_period_s))])
        final = np.median(window[-max(1, int(0.5 / sample_period_s)):])
        step_height = final - baseline
        if abs(step_height) < 1.0:
            continue  # not a real step (noise/edge of block)
        # jump instant = sample of steepest slope
        dy = np.diff(window)
        jump_idx = int(np.argmax(np.abs(dy)))
        pre_jump_disp = window[jump_idx] - baseline
        creep_frac = abs(pre_jump_disp / step_height)
        results.append({
            'move_t': move_t, 'move_abs_t': move_abs_t,
            'step_height_um': step_height,
            'pre_jump_disp_um': pre_jump_disp, 'creep_frac': creep_frac,
        })
    return results


def measure_c_block_approach(run_index, time_s, position_nm, sample_period_s, rows):
    """Fine-resolution look at the very start of C's first approach move."""
    start_s, end_s = find_block(rows, run_index, 'C')
    approach = next(
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['label'] == 'positive_approach'
        and start_s <= float(r['ids_time_s']) <= end_s
    )
    move_abs_t = float(approach['ids_time_s'])
    idx0 = int(round((move_abs_t - 0.05) / sample_period_s))
    idx1 = int(round((move_abs_t + 0.3) / sample_period_s))
    idx0, idx1 = max(0, idx0), min(position_nm.size, idx1)
    window = position_nm[idx0:idx1] / 1000.0
    t_window = (time_s[idx0:idx1] - move_abs_t) * 1000.0  # ms
    baseline = np.median(window[: max(1, int(0.02 / sample_period_s))])
    return t_window, window - baseline


def main():
    time_s, position_nm, sample_period_s, rows = load_common()

    print('=== D_0.125 tread analysis, MRES 1/1 (runs 5, 6) ===')
    all_creep_fracs = []
    for run_index in (5, 6):
        results = measure_d0125_treads(run_index, time_s, position_nm, sample_period_s, rows)
        print(f'\nRun {run_index}: {len(results)} isolated full-step treads')
        for r in results:
            print(f"  t={r['move_t']:6.2f}s  step={r['step_height_um']:+7.3f} um  "
                  f"pre-jump creep={r['pre_jump_disp_um']:+7.3f} um  "
                  f"creep_frac={r['creep_frac']:.3f}")
            all_creep_fracs.append(r['creep_frac'])

    median_creep_frac = np.median(all_creep_fracs) if all_creep_fracs else float('nan')
    print(f'\nMedian creep fraction (pre-jump displacement / total step) = '
          f'{median_creep_frac:.3f}')
    print('Rev-3-placeholder prediction (Fs/sigma0=5.9um out of 10um step) '
          'would give creep_frac ~ 0.59')

    # sigma0_way bracket: Fs_way / (creep_frac * 10um), using the measured
    # creep fraction and the current Fs_way=4.5N as a fixed reference point.
    fs_way_placeholder = 4.5
    if median_creep_frac > 1.0e-3:
        implied_sigma0 = fs_way_placeholder / (median_creep_frac * 10.0e-6)
        print(f'\nImplied sigma0_way from measured creep fraction '
              f'(holding Fs_way={fs_way_placeholder} N fixed): '
              f'{implied_sigma0:.3g} N/m')

    print('\n=== C-block first approach move: fine-resolution presliding look ===')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0))
    for run_index, ax in zip((5, 6), axes):
        t_ms, y_um = measure_c_block_approach(run_index, time_s, position_nm, sample_period_s, rows)
        ax.plot(t_ms, y_um, color='#1f77b4', lw=1.0, marker='.', markersize=2)
        ax.axhline(0.0, color='#9a9a9a', lw=0.7)
        ax.axvline(0.0, color='#d62728', lw=0.8, linestyle='--', label='commanded move issued')
        ax.set_xlabel('Time since MOVE_ACK (ms)')
        ax.set_ylabel('Position (µm, baselined)')
        ax.set_title(f'Run {run_index}: C-block first approach move (-40 µm commanded)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle('Step 4b -- presliding look at C-block approach onset', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_path = OUT_DIR / 'step4b_c_approach_presliding.png'
    fig.savefig(out_path, dpi=150)
    print(f'Saved {out_path}')

    fig2, axes2 = plt.subplots(1, 2, figsize=(13.0, 5.0))
    for run_index, ax in zip((5, 6), axes2):
        results = measure_d0125_treads(run_index, time_s, position_nm, sample_period_s, rows)
        for r in results:
            idx0 = int(round((r['move_abs_t'] - 1.0) / sample_period_s))
            idx1 = int(round((r['move_abs_t'] + 1.0) / sample_period_s))
            idx0, idx1 = max(0, idx0), min(position_nm.size, idx1)
            window = position_nm[idx0:idx1] / 1000.0
            t_window = (time_s[idx0:idx1] - r['move_abs_t']) * 1000.0
            baseline = np.median(window[:50])
            ax.plot(t_window, window - baseline, lw=0.9, alpha=0.8,
                    label=f"t={r['move_t']:.1f}s")
        ax.axvline(0.0, color='#d62728', lw=0.8, linestyle='--')
        ax.axhline(0.0, color='#9a9a9a', lw=0.7)
        ax.set_xlabel('Time since MOVE_ACK (ms)')
        ax.set_ylabel('Position (µm, baselined)')
        ax.set_title(f'Run {run_index}: D_0.125 isolated full-step treads (MRES 1/1)')
        ax.legend(fontsize=6.5, loc='lower right')
        ax.grid(True, alpha=0.3)
    fig2.suptitle('Step 4a -- D_0.125 tread shape at MRES 1/1', fontsize=11)
    fig2.tight_layout(rect=(0, 0, 1, 0.94))
    out_path2 = OUT_DIR / 'step4a_d0125_treads.png'
    fig2.savefig(out_path2, dpi=150)
    print(f'Saved {out_path2}')


if __name__ == '__main__':
    main()
