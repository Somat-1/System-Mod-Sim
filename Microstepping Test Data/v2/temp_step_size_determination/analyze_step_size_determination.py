#!/usr/bin/env python3
'''Overlay the integer MR probe commands on the IDS encoder measurement.'''

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'scripts' / 'StepSizeDepermination.csv'
OUTPUT_DIR = ROOT / 'temp_step_size_determination'
PLOT = OUTPUT_DIR / 'step_size_command_overlay.png'
SUMMARY = OUTPUT_DIR / 'step_size_summary.json'

# IDS linear encoder is nm-native (1 count = 1 nm); everything downstream of
# load_ids() is normalized to micrometres.
ENCODER_NM_PER_COUNT = 1.0
NM_PER_UM = 1.0e3
EXPECTED_COMMAND_GAP_S = 3.336

# Probe configuration from probe_dedicated_controller_integer_unit.py:
# `SM X 4 200 0` -> MRES 1/4, 200 full steps/rev. Lead screw pitch is 2 mm/rev.
LEAD_PITCH_MM = 2.0
FULL_STEPS_PER_REV = 200
MRES = 4
THEORETICAL_MICROSTEP_UM = (LEAD_PITCH_MM * 1000.0) / (FULL_STEPS_PER_REV * MRES)


def load_ids(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[int, int]] = []
    pattern = re.compile(r'^\s*(\d+)\t(\d+)\s*$')
    for line in path.read_text(
        encoding='utf-8-sig', errors='replace'
    ).splitlines():
        match = pattern.match(line)
        if match:
            rows.append((int(match.group(1)), int(match.group(2))))
    if not rows:
        raise RuntimeError(f'No numeric IDS samples found in {path}')
    array = np.asarray(rows, dtype=np.int64)
    time_s = array[:, 0].astype(float) / 1000.0
    raw = array[:, 1]
    increments = np.diff(raw)
    increments = np.where(
        increments > 2**31, increments - 2**32, increments
    )
    increments = np.where(
        increments < -(2**31), increments + 2**32, increments
    )
    relative_counts = np.r_[0, np.cumsum(increments)].astype(float)
    position_um = relative_counts * ENCODER_NM_PER_COUNT / NM_PER_UM
    return time_s, position_um


def robust_level(
    time_s: np.ndarray, position_um: np.ndarray, start: float, end: float,
) -> tuple[float, float]:
    values = position_um[(time_s >= start) & (time_s < end)]
    if values.size < 10:
        raise RuntimeError('Plateau window contains too few samples')
    median = float(np.median(values))
    mad_sigma = float(1.4826 * np.median(np.abs(values - median)))
    return median, mad_sigma


def infer_alignment(
    time_s: np.ndarray, position_um: np.ndarray,
) -> dict[str, float]:
    best: dict[str, float] | None = None
    for first_s in np.arange(1.0, time_s[-1] - 5.0, 0.001):
        second_s = first_s + EXPECTED_COMMAND_GAP_S
        pre, noise_pre = robust_level(
            time_s, position_um, first_s - 0.70, first_s - 0.15
        )
        plateau, noise_plateau = robust_level(
            time_s, position_um, first_s + 0.25, second_s - 0.25
        )
        post, noise_post = robust_level(
            time_s, position_um, second_s + 0.20, second_s + 0.75
        )
        step_um = ((plateau - pre) + (plateau - post)) / 2.0
        return_error = abs(post - pre)
        noise = max(0.002, (noise_pre + noise_plateau + noise_post) / 3.0)
        score = abs(step_um) / (noise + return_error)
        candidate = {
            'score': score,
            'first_command_s': first_s,
            'second_command_s': second_s,
            'pre_level_um': pre,
            'plateau_level_um': plateau,
            'post_level_um': post,
            'step_um': step_um,
            'return_error_um': return_error,
            'noise_pre_um': noise_pre,
            'noise_plateau_um': noise_plateau,
            'noise_post_um': noise_post,
        }
        if best is None or candidate['score'] > best['score']:
            best = candidate
    if best is None:
        raise RuntimeError('Could not infer command alignment')
    return best


def make_plot(
    time_s: np.ndarray, position_um: np.ndarray,
    result: dict[str, float],
) -> None:
    first_s = result['first_command_s']
    second_s = result['second_command_s']
    pre = result['pre_level_um']
    plateau = result['plateau_level_um']
    post = result['post_level_um']
    smooth = median_filter(position_um, size=101, mode='nearest')
    command_t = np.array([
        time_s[0], first_s, first_s, second_s, second_s, time_s[-1]
    ])
    command_y = np.array([pre, pre, plateau, plateau, post, post])

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 7.5), constrained_layout=True,
        gridspec_kw={'height_ratios': [1.0, 1.25]},
    )
    ax = axes[0]
    ax.plot(time_s, position_um, color='0.72', linewidth=0.45,
            label='IDS, 1 nm/count')
    ax.plot(time_s, smooth, color='#18794e', linewidth=1.3,
            label='101 ms median')
    ax.plot(command_t, command_y, color='#d97706', linewidth=2.0,
            label='Aligned MR X 1 / MR X -1 command')
    ax.axvspan(first_s, second_s, color='#f59e0b', alpha=0.10)
    ax.set(xlabel='Time from IDS export start (s)',
           ylabel='Relative displacement (µm)',
           title='Integer-position calibration: full IDS record')
    ax.grid(True, alpha=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right')

    ax = axes[1]
    zoom = (time_s >= first_s - 1.5) & (time_s <= second_s + 1.5)
    ax.plot(time_s[zoom], position_um[zoom], color='0.72', linewidth=0.55,
            label='IDS samples')
    ax.plot(time_s[zoom], smooth[zoom], color='#18794e', linewidth=1.6,
            label='101 ms median')
    ax.plot(command_t, command_y, color='#d97706', linewidth=2.2,
            label='Command overlay')
    ax.axvline(first_s, color='#b45309', linestyle='--', linewidth=1.1)
    ax.axvline(second_s, color='#b45309', linestyle='--', linewidth=1.1)
    ax.text(first_s, ax.get_ylim()[1], ' MR X 1', va='top', ha='left',
            color='#92400e')
    ax.text(second_s, ax.get_ylim()[1], ' MR X -1', va='top', ha='left',
            color='#92400e')
    uncertainty = max(
        0.003,
        (result['noise_pre_um'] + result['noise_plateau_um']
         + result['noise_post_um']) / 3.0,
    )
    step_um = result['step_um']
    return_error_um = result['return_error_um']
    deviation_pct = result['deviation_from_theoretical_pct']
    note = (
        f'Plateau medians: {pre:.3f} → {plateau:.3f} → '
        f'{post:.3f} µm\n'
        f'Measured one-count step: {step_um:.3f} '
        f'± {uncertainty:.3f} µm\n'
        f'Theoretical microstep (2 mm lead, 1/{MRES}, '
        f'{FULL_STEPS_PER_REV} steps/rev): {THEORETICAL_MICROSTEP_UM:.3f} '
        f'µm\n'
        f'Deviation from theoretical: {deviation_pct:+.1f}%\n'
        f'Return residual: {return_error_um:.3f} µm'
    )
    ax.text(0.02, 0.96, note, transform=ax.transAxes, va='top', ha='left',
            bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.88})
    ax.set_xlim(first_s - 1.5, second_s + 1.5)
    ax.set(xlabel='Time from IDS export start (s)',
           ylabel='Relative displacement (µm)',
           title='Inferred command alignment and robust plateau levels')
    ax.grid(True, alpha=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right')
    fig.suptitle(
        'Dedicated-controller integer unit: inferred physical displacement',
        fontsize=14, fontweight='semibold',
    )
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    time_s, position_um = load_ids(INPUT)
    result = infer_alignment(time_s, position_um)
    uncertainty = max(
        0.003,
        (result['noise_pre_um'] + result['noise_plateau_um']
         + result['noise_post_um']) / 3.0,
    )
    deviation_pct = (
        (result['step_um'] - THEORETICAL_MICROSTEP_UM)
        / THEORETICAL_MICROSTEP_UM * 100.0
    )
    result['deviation_from_theoretical_pct'] = deviation_pct
    result.update({
        'input': str(INPUT),
        'sample_count': int(time_s.size),
        'sample_period_ms': 1.0,
        'encoder_scale_nm_per_count': ENCODER_NM_PER_COUNT,
        'expected_command_gap_s': EXPECTED_COMMAND_GAP_S,
        'estimated_uncertainty_um': uncertainty,
        'lead_pitch_mm_per_rev': LEAD_PITCH_MM,
        'full_steps_per_rev': FULL_STEPS_PER_REV,
        'mres': MRES,
        'theoretical_microstep_um': THEORETICAL_MICROSTEP_UM,
        'interpretation': (
            'One controller integer MR unit produced approximately the '
            'reported step_um displacement. Command times are inferred '
            'because the standalone probe did not timestamp commands. '
            'The theoretical microstep size follows from the 2 mm lead '
            'screw pitch, 200 full steps/rev, and MRES 1/4 used in the '
            'probe (SM X 4 200 0); the deviation quantifies how well the '
            'commanded single-microstep move is actually followed by the '
            'stage, versus the ideal kinematic prediction.'
        ),
    })
    make_plot(time_s, position_um, result)
    SUMMARY.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    print(f'Plot: {PLOT}')
    print(f'Summary: {SUMMARY}')


if __name__ == '__main__':
    main()
