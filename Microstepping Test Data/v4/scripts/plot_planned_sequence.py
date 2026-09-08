#!/usr/bin/env python3
"""Render the complete intended v4 MRES/trajectory measurement sequence.

The large panel shows all four MRES runs and every 25 mm trajectory. The
four zoom panels show the 15 one-microstep oscillations that are too small to
see on the 25 mm scale. A recent dry-run log is checked before rendering so
the preview cannot silently drift away from the executable campaign.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_mres_trajectory_campaign import (  # noqa: E402
    CAMPAIGN_LEAD_IN_S,
    CAMPAIGN_TAIL_S,
    MRES_VALUES,
    OSCILLATION_CYCLES,
    OSCILLATION_HALF_DWELL_S,
    OSCILLATION_MOVE_RATE_FULL_STEPS_S,
    TRAJECTORY_DISTANCE_MM,
    TRAJECTORY_ENDPOINT_DWELL_S,
    TRAJECTORY_FULL_STEPS,
    TRAJECTORY_RATE_NAMES,
    TRAJECTORY_RATES_FULL_STEPS_S,
    marker_amplitude,
    trajectory_block_name,
)
from run_settling_dedicated_controller import (  # noqa: E402
    MARKER_RATE_FULL_STEPS_S,
    MARKER_REVERSE_DWELL_S,
    MARKER_SETTLE_S,
)


FULL_STEP_MM = 0.010
COLORS = {
    'lead': '#8a8a8a',
    'marker': '#8e5bb7',
    'oscillation': '#2a9d67',
    'direct': '#2474b5',
    'individual': '#dd7f24',
}


@dataclass(frozen=True)
class Segment:
    t0: float
    t1: float
    p0_mm: float
    p1_mm: float
    kind: str
    label: str
    mres: int | None = None


@dataclass(frozen=True)
class Block:
    t0: float
    t1: float
    kind: str
    label: str
    mres: int | None = None


class Plan:
    def __init__(self) -> None:
        self.t = 0.0
        self.position_mm = 0.0
        self.segments: list[Segment] = []
        self.blocks: list[Block] = []
        self.oscillation_blocks: dict[int, Block] = {}
        self.config_spans: dict[int, tuple[float, float]] = {}
        self.marker_index = 0

    def move(
        self, target_mm: float, duration_s: float, kind: str, label: str,
        mres: int | None = None,
    ) -> None:
        self.segments.append(
            Segment(
                self.t, self.t + duration_s, self.position_mm, target_mm,
                kind, label, mres,
            )
        )
        self.t += duration_s
        self.position_mm = target_mm

    def dwell(
        self, duration_s: float, kind: str, label: str,
        mres: int | None = None,
    ) -> None:
        self.move(self.position_mm, duration_s, kind, label, mres)

    def marker(self, label: str, mres: int) -> None:
        self.marker_index += 1
        amplitude_steps = marker_amplitude(self.marker_index)
        amplitude_mm = amplitude_steps * FULL_STEP_MM
        t0 = self.t
        duration = amplitude_steps / MARKER_RATE_FULL_STEPS_S
        self.move(-amplitude_mm, duration, 'marker', label, mres)
        self.dwell(MARKER_REVERSE_DWELL_S, 'marker', label, mres)
        self.move(0.0, duration, 'marker', label, mres)
        self.dwell(MARKER_SETTLE_S, 'marker', label, mres)
        self.blocks.append(Block(t0, self.t, 'marker', label, mres))

    def oscillation(self, mres: int) -> None:
        label = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        t0 = self.t
        displacement_mm = FULL_STEP_MM / mres
        move_s = 1.0 / (mres * OSCILLATION_MOVE_RATE_FULL_STEPS_S)
        for cycle in range(1, OSCILLATION_CYCLES + 1):
            self.move(
                displacement_mm, move_s, 'oscillation',
                f'{label}_cycle_{cycle:02d}_forward', mres,
            )
            self.dwell(
                OSCILLATION_HALF_DWELL_S, 'oscillation',
                f'{label}_cycle_{cycle:02d}_forward_dwell', mres,
            )
            self.move(
                0.0, move_s, 'oscillation',
                f'{label}_cycle_{cycle:02d}_return', mres,
            )
            self.dwell(
                OSCILLATION_HALF_DWELL_S, 'oscillation',
                f'{label}_cycle_{cycle:02d}_return_dwell', mres,
            )
        block = Block(t0, self.t, 'oscillation', label, mres)
        self.blocks.append(block)
        self.oscillation_blocks[mres] = block

    def trajectory(
        self, mres: int, mode: str, rate_name: str, rate: float
    ) -> None:
        label = trajectory_block_name(mres, mode, rate_name, rate)
        t0 = self.t
        leg_s = TRAJECTORY_FULL_STEPS / rate
        self.move(float(TRAJECTORY_DISTANCE_MM), leg_s, mode, label, mres)
        self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, mode, label, mres)
        self.move(0.0, leg_s, mode, label, mres)
        self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, mode, label, mres)
        self.blocks.append(Block(t0, self.t, mode, label, mres))


def build_plan() -> Plan:
    plan = Plan()
    plan.dwell(CAMPAIGN_LEAD_IN_S, 'lead', 'campaign_lead_in')
    for run_index, mres in enumerate(MRES_VALUES, start=1):
        config_t0 = plan.t
        plan.marker(f'CONFIG_{run_index:02d}_MRES_{mres}', mres)
        oscillation_label = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        plan.marker(oscillation_label, mres)
        plan.oscillation(mres)
        for mode in ('direct', 'individual'):
            for rate_name, rate in zip(
                TRAJECTORY_RATE_NAMES, TRAJECTORY_RATES_FULL_STEPS_S
            ):
                label = trajectory_block_name(mres, mode, rate_name, rate)
                plan.marker(label, mres)
                plan.trajectory(mres, mode, rate_name, rate)
        plan.config_spans[mres] = (config_t0, plan.t)
    plan.dwell(CAMPAIGN_TAIL_S, 'lead', 'campaign_tail')
    return plan


def latest_dry_run() -> Path | None:
    log_dir = HERE.parent / 'data' / 'hardware_runs'
    candidates = sorted(log_dir.glob('mres_trajectory_dry_run_*.csv'))
    return candidates[-1] if candidates else None


def validate_dry_run(path: Path) -> None:
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not any(row['event'] == 'CAMPAIGN_COMPLETE' for row in rows):
        raise SystemExit(f'Dry-run did not complete: {path}')
    run_mres = [
        int(row['mres']) for row in rows if row['event'] == 'RUN_CONFIG'
    ]
    if tuple(run_mres) != MRES_VALUES:
        raise SystemExit(
            f'Dry-run MRES order {tuple(run_mres)} != planned {MRES_VALUES}'
        )
    expected_trajectories = len(MRES_VALUES) * 2 * len(
        TRAJECTORY_RATES_FULL_STEPS_S
    )
    trajectory_starts = sum(
        row['event'] == 'TRAJECTORY_START' for row in rows
    )
    if trajectory_starts != expected_trajectories:
        raise SystemExit(
            f'Dry-run has {trajectory_starts} trajectories; '
            f'expected {expected_trajectories}'
        )
    oscillation_starts = sum(
        row['event'] == 'OSCILLATION_START' for row in rows
    )
    if oscillation_starts != len(MRES_VALUES):
        raise SystemExit(
            f'Dry-run has {oscillation_starts} oscillations; '
            f'expected {len(MRES_VALUES)}'
        )


def plot_full_plan(plan: Plan, out_path: Path) -> None:
    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    grid = fig.add_gridspec(4, 4, height_ratios=(3.8, 0.75, 1.65, 1.65))
    overview = fig.add_subplot(grid[0, :])
    schedule = fig.add_subplot(grid[1, :], sharex=overview)
    oscillation_axes = [
        fig.add_subplot(
            grid[2 + index // 2, (index % 2) * 2:(index % 2) * 2 + 2]
        )
        for index in range(len(MRES_VALUES))
    ]

    for config_index, (mres, (t0, t1)) in enumerate(plan.config_spans.items()):
        if config_index % 2 == 0:
            overview.axvspan(
                t0 / 60.0, t1 / 60.0, color='#e8eef4', alpha=0.55
            )
        overview.axvline(t0 / 60.0, color='#5d6872', lw=0.7, alpha=0.55)
        overview.text(
            (t0 + t1) / 120.0, 26.2, f'MRES {mres}', ha='center',
            va='bottom', fontsize=10, weight='bold', color='#34424e',
        )

    for segment in plan.segments:
        overview.plot(
            [segment.t0 / 60.0, segment.t1 / 60.0],
            [segment.p0_mm, segment.p1_mm],
            color=COLORS[segment.kind],
            lw=1.45 if segment.kind in {'direct', 'individual'} else 1.0,
            alpha=0.95,
        )

    overview.axhline(0.0, color='#3a3a3a', lw=0.6)
    overview.set_xlim(0.0, plan.t / 60.0)
    overview.set_ylim(-1.65, 27.4)
    overview.set_ylabel('Ideal stage position (mm)')
    overview.grid(True, alpha=0.22)
    overview.tick_params(labelbottom=False)

    lane_y = {
        'marker': 0.05,
        'oscillation': 0.26,
        'direct': 0.47,
        'individual': 0.68,
    }
    lane_h = 0.16
    for block in plan.blocks:
        schedule.add_patch(
            Rectangle(
                (block.t0 / 60.0, lane_y[block.kind]),
                (block.t1 - block.t0) / 60.0,
                lane_h,
                facecolor=COLORS[block.kind], edgecolor='none', alpha=0.88,
            )
        )
        if block.kind in {'direct', 'individual'}:
            short = block.label.split('_')[-2][0]
            schedule.text(
                (block.t0 + block.t1) / 120.0,
                lane_y[block.kind] + lane_h / 2,
                short, ha='center', va='center', color='white',
                fontsize=7, weight='bold',
            )
    schedule.set_yticks(
        [
            lane_y[k] + lane_h / 2
            for k in ('marker', 'oscillation', 'direct', 'individual')
        ],
        ['marker', '30 s oscillation', 'direct 25 mm', 'individual 25 mm'],
    )
    schedule.set_ylim(0.0, 0.90)
    schedule.set_xlabel('Measurement time (min)')
    schedule.grid(axis='x', alpha=0.2)

    for axis, mres in zip(oscillation_axes, MRES_VALUES):
        block = plan.oscillation_blocks[mres]
        for segment in plan.segments:
            if segment.kind != 'oscillation' or segment.mres != mres:
                continue
            axis.plot(
                [segment.t0 - block.t0, segment.t1 - block.t0],
                [segment.p0_mm * 1000.0, segment.p1_mm * 1000.0],
                color=COLORS['oscillation'], lw=1.35,
            )
        amplitude_um = FULL_STEP_MM * 1000.0 / mres
        axis.set_xlim(0.0, block.t1 - block.t0)
        axis.set_ylim(-0.05 * amplitude_um, 1.20 * amplitude_um)
        axis.set_title(
            f'MRES {mres}: 15 complete one-microstep cycles in 30 s '
            f'({amplitude_um:g} µm amplitude)',
            fontsize=9.5,
        )
        axis.set_xlabel('Time within oscillation block (s)')
        axis.set_ylabel('Displacement (µm)')
        axis.grid(True, alpha=0.25)

    total_min = plan.t / 60.0
    legend = [
        Patch(facecolor=COLORS['marker'], label='unique separation marker'),
        Patch(facecolor=COLORS['oscillation'], label='one-microstep oscillation'),
        Patch(facecolor=COLORS['direct'], label='direct distance command'),
        Patch(
            facecolor=COLORS['individual'],
            label='one full-step command at a time',
        ),
    ]
    overview.legend(handles=legend, loc='lower right', ncols=2, fontsize=8.5)
    fig.suptitle(
        'v4 intended complete measurement trajectory — '
        f'{total_min:.2f} min planned (55 min recording limit)\n'
        'MRES 1, 4, 16, 32; 15-cycle/30 s microstep oscillation; '
        '25 mm out-and-back at 27.5, 70, and 200 full steps/s; '
        'direct and individually commanded approaches\n'
        'Each colored block is preceded by a unique negative marker; '
        'S/M/F in trajectory strips denote slow/moderate/fast.',
        fontsize=12,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Render the complete planned v4 MRES trajectory campaign.'
    )
    parser.add_argument(
        '--log-file', type=Path,
        help='Dry-run CSV to validate; defaults to the newest campaign dry run.',
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    log_path = args.log_file or latest_dry_run()
    if log_path is None:
        raise SystemExit(
            'No campaign dry-run log found. Run '
            '`python run_mres_trajectory_campaign.py --dry-run` first.'
        )
    validate_dry_run(log_path)
    plan = build_plan()
    if plan.t >= 55.0 * 60.0:
        raise SystemExit(
            f'Planned sequence is {plan.t / 60.0:.2f} min, exceeding 55 min.'
        )
    print(f'Validated dry-run log: {log_path}')
    print(
        f'Complete planned duration: {plan.t:.1f} s '
        f'({plan.t / 60.0:.2f} min)'
    )

    out_dir = HERE.parent / 'rendered_assets'
    campaign_path = out_dir / 'planned_mres_trajectory_campaign.png'
    compatibility_path = out_dir / 'planned_settling_sequence_preview.png'
    plot_full_plan(plan, campaign_path)
    plot_full_plan(plan, compatibility_path)
    print(f'Saved: {campaign_path}')
    print(f'Updated existing preview: {compatibility_path}')


if __name__ == '__main__':
    main()
