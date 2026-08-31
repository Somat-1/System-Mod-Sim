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
COUNTS_TO_NM = 1.0
EXPECTED_COMMAND_GAP_S = 3.336


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
    return time_s, relative_counts * COUNTS_TO_NM


def robust_level(
    time_s: np.ndarray, position_nm: np.ndarray, start: float, end: float,
) -> tuple[float, float]:
    values = position_nm[(time_s >= start) & (time_s < end)]
    if values.size < 10:
        raise RuntimeError('Plateau window contains too few samples')
    median = float(np.median(values))
    mad_sigma = float(1.4826 * np.median(np.abs(values - median)))
    return median, mad_sigma


def infer_alignment(
    time_s: np.ndarray, position_nm: np.ndarray,
) -> dict[str, float]:
    best: dict[str, float] | None = None
    for first_s in np.arange(1.0, time_s[-1] - 5.0, 0.001):
        second_s = first_s + EXPECTED_COMMAND_GAP_S
        pre, noise_pre = robust_level(
            time_s, position_nm, first_s - 0.70, first_s - 0.15
        )
        plateau, noise_plateau = robust_level(
            time_s, position_nm, first_s + 0.25, second_s - 0.25
        )
        post, noise_post = robust_level(
            time_s, position_nm, second_s + 0.20, second_s + 0.75
        )
        step_nm = ((plateau - pre) + (plateau - post)) / 2.0
        return_error = abs(post - pre)
        noise = max(2.0, (noise_pre + noise_plateau + noise_post) / 3.0)
        score = abs(step_nm) / (noise + return_error)
        candidate = {
            'score': score,
            'first_command_s': first_s,
            'second_command_s': second_s,
            'pre_level_nm': pre,
            'plateau_level_nm': plateau,
            'post_level_nm': post,
            'step_nm': step_nm,
            'return_error_nm': return_error,
            'noise_pre_nm': noise_pre,
            'noise_plateau_nm': noise_plateau,
            'noise_post_nm': noise_post,
        }
        if best is None or candidate['score'] > best['score']:
            best = candidate
    if best is None:
        raise RuntimeError('Could not infer command alignment')
    return best


def make_plot(
    time_s: np.ndarray, position_nm: np.ndarray,
    result: dict[str, float],
) -> None:
    first_s = result['first_command_s']
    second_s = result['second_command_s']
    pre = result['pre_level_nm']
    plateau = result['plateau_level_nm']
    post = result['post_level_nm']
    smooth = median_filter(position_nm, size=101, mode='nearest')
    command_t = np.array([
        time_s[0], first_s, first_s, second_s, second_s, time_s[-1]
    ])
    command_y = np.array([pre, pre, plateau, plateau, post, post])

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 7.5), constrained_layout=True,
        gridspec_kw={'height_ratios': [1.0, 1.25]},
    )
    ax = axes[0]
    ax.plot(time_s, position_nm, color='0.72', linewidth=0.45,
            label='IDS, 1 nm/count')
    ax.plot(time_s, smooth, color='#18794e', linewidth=1.3,
            label='101 ms median')
    ax.plot(command_t, command_y, color='#d97706', linewidth=2.0,
            label='Aligned MR X 1 / MR X -1 command')
    ax.axvspan(first_s, second_s, color='#f59e0b', alpha=0.10)
    ax.set(xlabel='Time from IDS export start (s)',
           ylabel='Relative displacement (nm)',
           title='Integer-position calibration: full IDS record')
    ax.grid(True, alpha=0.25)
    ax.legend(loc='upper right')

    ax = axes[1]
    zoom = (time_s >= first_s - 1.5) & (time_s <= second_s + 1.5)
    ax.plot(time_s[zoom], position_nm[zoom], color='0.72', linewidth=0.55,
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
        3.0,
        (result['noise_pre_nm'] + result['noise_plateau_nm']
         + result['noise_post_nm']) / 3.0,
    )
    step_nm = result['step_nm']
    return_error_nm = result['return_error_nm']
    note = (
        f'Plateau medians: {pre:.0f} -> {plateau:.0f} -> {post:.0f} nm\n'
        f'Measured one-count step: {step_nm:.0f} '
        f'+/- {uncertainty:.0f} nm\n'
        f'Return residual: {return_error_nm:.0f} nm'
    )
    ax.text(0.02, 0.96, note, transform=ax.transAxes, va='top', ha='left',
            bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.88})
    ax.set_xlim(first_s - 1.5, second_s + 1.5)
    ax.set(xlabel='Time from IDS export start (s)',
           ylabel='Relative displacement (nm)',
           title='Inferred command alignment and robust plateau levels')
    ax.grid(True, alpha=0.25)
    ax.legend(loc='lower right')
    fig.suptitle(
        'Dedicated-controller integer unit: inferred physical displacement',
        fontsize=14,
    )
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    time_s, position_nm = load_ids(INPUT)
    result = infer_alignment(time_s, position_nm)
    uncertainty = max(
        3.0,
        (result['noise_pre_nm'] + result['noise_plateau_nm']
         + result['noise_post_nm']) / 3.0,
    )
    result.update({
        'input': str(INPUT),
        'sample_count': int(time_s.size),
        'sample_period_ms': 1.0,
        'encoder_scale_nm_per_count': COUNTS_TO_NM,
        'expected_command_gap_s': EXPECTED_COMMAND_GAP_S,
        'estimated_uncertainty_nm': uncertainty,
        'interpretation': (
            'One controller integer MR unit produced approximately the '
            'reported step_nm displacement. Command times are inferred '
            'because the standalone probe did not timestamp commands.'
        ),
    })
    make_plot(time_s, position_nm, result)
    SUMMARY.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    print(f'Plot: {PLOT}')
    print(f'Summary: {SUMMARY}')


if __name__ == '__main__':
    main()
