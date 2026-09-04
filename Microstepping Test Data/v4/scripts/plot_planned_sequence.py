#!/usr/bin/env python3
"""Preview the v4 settling-sweep trajectory from a --dry-run log -- no
hardware involved. Reconstructs the ideal commanded position over time
for ONE representative configuration (one MRES/current pair; the full
campaign repeats this same pattern once per configuration -- see
run_settling_dedicated_controller.py's run_campaign) and shades every
block so the segmentation structure (marker -> outbound settle ->
return settle, repeated per test distance) is visible before running
anything on real hardware.

Reconstruction note: DryRunTransport.command() does not advance the
virtual clock for MA/MR (only dwell()'s explicit clock.sleep() does), so
move ACKs in the dry-run log are logged at the same instant as whatever
preceded them. Move durations are therefore computed analytically here
(full_steps / configured rate) rather than read off monotonic_ns; dwell
durations, by contrast, ARE correctly captured by the virtual clock and
are read directly from each DWELL_START/DWELL_END pair.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_settling_dedicated_controller import (  # noqa: E402
    MOTOR_FULL_STEPS_PER_REV, TEST_DISTANCES_FULL_STEPS, SETTLE_DWELL_S,
)

LEAD_MM_PER_REV = 2.0  # matches Rev 4.2 model_parameters.json's L=2e-3 m
BLOCK_COLORS = {
    'reference': '#7f7f7f', 'marker': '#9467bd', 'settle_outbound': '#1f77b4',
    'settle_return': '#d62728', 'other': '#2ca02c',
}


def classify_block(block_name: str) -> str | None:
    if block_name.startswith('BLOCK_0'):
        return 'reference'
    if block_name.startswith('MARKER_'):
        return 'marker'
    return None  # settling-test blocks: left unshaded on purpose, see BLOCK_COLORS


def load_rows(path: Path, run_index: int) -> list[dict]:
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        rows = [r for r in csv.DictReader(handle) if r['run_index'] == str(run_index)]
    if not rows:
        raise SystemExit(f'No rows found for run_index={run_index} in {path}')
    return rows


def reconstruct_timeline(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Returns (segments, blocks). segments: dicts with t0,t1,pos0_rev,pos1_rev,
    kind,label. blocks: dicts with t0,t1,name,color, for shading."""
    segments: list[dict] = []
    blocks: list[dict] = []
    open_blocks: dict[str, float] = {}
    current_rate_full_steps_s = None
    t = 0.0
    pos_rev = 0.0

    i = 0
    n = len(rows)
    while i < n:
        row = rows[i]
        event = row['event']
        if event == 'DRY_COMMAND' and row['command'].startswith('SS '):
            tokens = row['command'].split()
            maximum_code = int(tokens[3])
            current_rate_full_steps_s = (
                maximum_code * 0.01 * MOTOR_FULL_STEPS_PER_REV / (2.0 * math.pi)
            )
        elif event == 'BLOCK_START':
            open_blocks[row['block']] = t
        elif event == 'BLOCK_END':
            start_t = open_blocks.pop(row['block'], t)
            blocks.append({'t0': start_t, 't1': t, 'name': row['block']})
        elif event == 'MOVE_ACK':
            new_pos_rev = float(row['ideal_position_rev'])
            rate = current_rate_full_steps_s
            if row['rate_full_steps_s']:
                rate = float(row['rate_full_steps_s'])
            delta_rev = new_pos_rev - pos_rev
            full_steps = abs(delta_rev) * MOTOR_FULL_STEPS_PER_REV
            duration = full_steps / rate if rate else 0.0
            segments.append({
                't0': t, 't1': t + duration, 'pos0_rev': pos_rev,
                'pos1_rev': new_pos_rev, 'kind': 'move', 'label': row['label'],
                'block': row['block'],
            })
            t += duration
            pos_rev = new_pos_rev
        elif event == 'DWELL_START':
            start_t, start_mono = t, int(row['monotonic_ns'])
            end_row = rows[i + 1]
            assert end_row['event'] == 'DWELL_END' and end_row['label'] == row['label'], (
                f'Expected paired DWELL_END for {row["label"]!r}, got {end_row["event"]!r}'
            )
            duration = (int(end_row['monotonic_ns']) - start_mono) / 1.0e9
            segments.append({
                't0': start_t, 't1': start_t + duration, 'pos0_rev': pos_rev,
                'pos1_rev': pos_rev, 'kind': 'dwell', 'label': row['label'],
                'block': row['block'],
            })
            t += duration
            i += 1  # consume the paired DWELL_END too
        i += 1
    return segments, blocks


def main() -> None:
    log_candidates = sorted(
        (HERE.parents[0] / 'data' / 'hardware_runs').glob('settling_dry_run_*.csv')
    )
    if not log_candidates:
        raise SystemExit(
            'No dry-run log found; run '
            '`python run_settling_dedicated_controller.py --dry-run` first.'
        )
    log_path = log_candidates[-1]
    run_index = 1
    print(f'Reading {log_path} (run_index={run_index})', flush=True)

    rows = load_rows(log_path, run_index)
    segments, blocks = reconstruct_timeline(rows)
    total_s = segments[-1]['t1']
    print(f'Reconstructed {len(segments)} segments, {len(blocks)} blocks, '
          f'total duration {total_s:.1f} s ({total_s / 60.0:.1f} min)', flush=True)

    fig, ax = plt.subplots(figsize=(16.0, 6.0))

    for block in blocks:
        kind = classify_block(block['name'])
        if kind is None:
            continue
        color = BLOCK_COLORS[kind]
        alpha = 0.35 if kind == 'marker' else 0.12
        ax.axvspan(block['t0'], block['t1'], color=color, alpha=alpha, lw=0)

    for seg in segments:
        pos0_mm = seg['pos0_rev'] * LEAD_MM_PER_REV
        pos1_mm = seg['pos1_rev'] * LEAD_MM_PER_REV
        color = '#1f77b4' if seg['kind'] == 'move' else '#9a9a9a'
        lw = 1.6 if seg['kind'] == 'move' else 1.0
        ax.plot([seg['t0'], seg['t1']], [pos0_mm, pos1_mm], color=color, lw=lw)

    top_mm = max(s['pos1_rev'] for s in segments) * LEAD_MM_PER_REV
    ax.set_ylim(top=top_mm * 1.30 + 0.1)

    # Label each settling test's block span along the top.
    label_y = top_mm * 1.06 + 0.05
    for block in blocks:
        if block['name'].startswith('SETTLE_'):
            mid_t = 0.5 * (block['t0'] + block['t1'])
            ax.annotate(
                block['name'].replace('SETTLE_', ''), xy=(mid_t, 0),
                xytext=(mid_t, label_y), fontsize=7, ha='center', va='bottom',
                rotation=90, color='#1f77b4',
                arrowprops=dict(arrowstyle='-', color='#c7c7c7', lw=0.6),
            )

    ax.set_xlabel('Time within this configuration (s)')
    ax.set_ylabel('Ideal commanded stage position (mm)')
    ax.grid(True, alpha=0.3)
    ax.set_title(
        'v4 planned settling sweep -- ONE representative configuration '
        '(MRES/current pair; full campaign repeats this 6x)\n'
        f'{len(TEST_DISTANCES_FULL_STEPS)} test distances '
        f'({TEST_DISTANCES_FULL_STEPS[0]}..{TEST_DISTANCES_FULL_STEPS[-1]} full steps), '
        f'{SETTLE_DWELL_S:g} s settle dwell after every move, '
        f'duration = {total_s / 60.0:.1f} min\n'
        'Shading: grey = reference block, purple = segmentation marker '
        '(blue trace = commanded move, grey trace = dwell)',
        fontsize=10.5,
    )
    fig.tight_layout()

    out_dir = HERE.parents[0] / 'rendered_assets'
    out_path = out_dir / 'planned_settling_sequence_preview.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}', flush=True)


if __name__ == '__main__':
    main()
