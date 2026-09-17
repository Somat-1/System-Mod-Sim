#!/usr/bin/env python3
'''Render the v4 SpreadCycle + MicroPlyer MRES/trajectory campaign plan.

Clone of plot_planned_sequence.py for the MRES 1/4/16 SpreadCycle +
MicroPlyer variant (see SPREADCYCLE_MICROPLYER_VARIANT_MEMO.md and
scripts/esp32_v4_mres134_spreadcycle_microplyer_campaign/). MRES 32 is
dropped from this variant to make room for MicroPlyer (step interpolation
on) within the 55-minute recording limit; every rate, dwell, marker, and
separation for MRES 1/4/16 is otherwise identical to the baseline plan.
This script models timing/position only -- it does not encode driver
chopper mode -- and renders the complete plan: three oscillation panels,
the full 18-trajectory plot, and the overall schedule.
'''

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

from dedicated_controller_support import (  # noqa: E402
    MARKER_RATE_FULL_STEPS_S,
    MARKER_REVERSE_DWELL_S,
    MARKER_SETTLE_S,
)
from run_mres_trajectory_campaign import (  # noqa: E402
    CAMPAIGN_LEAD_IN_S,
    CAMPAIGN_TAIL_S,
    EXPERIMENT_SEPARATION_S,
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

# This variant drops MRES 32 (vs. the baseline's (1, 4, 16, 32)) to make
# room for MicroPlyer within the 55-minute recording limit.
MRES_VALUES = (1, 4, 16)

FULL_STEP_MM = 0.010
COLORS = {
    'lead': '#8a8a8a',
    'separation': '#c7ccd1',
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
        self,
        target_mm: float,
        duration_s: float,
        kind: str,
        label: str,
        mres: int | None = None,
    ) -> None:
        if duration_s < 0.0:
            raise ValueError('Segment duration cannot be negative')
        self.segments.append(
            Segment(
                self.t,
                self.t + duration_s,
                self.position_mm,
                target_mm,
                kind,
                label,
                mres,
            )
        )
        self.t += duration_s
        self.position_mm = target_mm

    def dwell(
        self,
        duration_s: float,
        kind: str,
        label: str,
        mres: int | None = None,
    ) -> None:
        self.move(self.position_mm, duration_s, kind, label, mres)

    def separation(self, next_label: str, mres: int) -> None:
        t0 = self.t
        self.dwell(
            EXPERIMENT_SEPARATION_S,
            'separation',
            f'before_{next_label}',
            mres,
        )
        self.blocks.append(
            Block(t0, self.t, 'separation', next_label, mres)
        )

    def marker(self, label: str, mres: int) -> None:
        self.marker_index += 1
        start_position = self.position_mm
        amplitude_steps = marker_amplitude(self.marker_index)
        amplitude_mm = amplitude_steps * FULL_STEP_MM
        t0 = self.t
        duration = amplitude_steps / MARKER_RATE_FULL_STEPS_S
        self.move(
            start_position - amplitude_mm,
            duration,
            'marker',
            label,
            mres,
        )
        self.dwell(MARKER_REVERSE_DWELL_S, 'marker', label, mres)
        self.move(start_position, duration, 'marker', label, mres)
        self.dwell(MARKER_SETTLE_S, 'marker', label, mres)
        if self.position_mm != start_position:
            raise AssertionError(f'Marker {label} did not return to its start')
        self.blocks.append(Block(t0, self.t, 'marker', label, mres))

    def oscillation(self, mres: int) -> None:
        label = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        t0 = self.t
        start_position = self.position_mm
        displacement_mm = FULL_STEP_MM / mres
        move_s = 1.0 / (mres * OSCILLATION_MOVE_RATE_FULL_STEPS_S)
        for cycle in range(1, OSCILLATION_CYCLES + 1):
            self.move(
                start_position + displacement_mm,
                move_s,
                'oscillation',
                f'{label}_cycle_{cycle:02d}_forward',
                mres,
            )
            self.dwell(
                OSCILLATION_HALF_DWELL_S,
                'oscillation',
                f'{label}_cycle_{cycle:02d}_forward_dwell',
                mres,
            )
            self.move(
                start_position,
                move_s,
                'oscillation',
                f'{label}_cycle_{cycle:02d}_return',
                mres,
            )
            self.dwell(
                OSCILLATION_HALF_DWELL_S,
                'oscillation',
                f'{label}_cycle_{cycle:02d}_return_dwell',
                mres,
            )
        if self.position_mm != start_position:
            raise AssertionError(
                f'{label} ended at {self.position_mm:g} mm, '
                f'not its {start_position:g} mm start'
            )
        block = Block(t0, self.t, 'oscillation', label, mres)
        self.blocks.append(block)
        self.oscillation_blocks[mres] = block

    def trajectory(
        self,
        mres: int,
        mode: str,
        rate_name: str,
        rate: float,
    ) -> None:
        label = trajectory_block_name(mres, mode, rate_name, rate)
        t0 = self.t
        start_position = self.position_mm
        leg_s = TRAJECTORY_FULL_STEPS / rate
        self.move(
            start_position + float(TRAJECTORY_DISTANCE_MM),
            leg_s,
            mode,
            label,
            mres,
        )
        self.dwell(
            TRAJECTORY_ENDPOINT_DWELL_S,
            mode,
            f'{label}_outbound_endpoint',
            mres,
        )
        self.move(start_position, leg_s, mode, label, mres)
        self.dwell(
            TRAJECTORY_ENDPOINT_DWELL_S,
            mode,
            f'{label}_return_endpoint',
            mres,
        )
        if self.position_mm != start_position:
            raise AssertionError(
                f'{label} did not return to its start position'
            )
        self.blocks.append(Block(t0, self.t, mode, label, mres))


def build_plan() -> Plan:
    plan = Plan()
    plan.dwell(CAMPAIGN_LEAD_IN_S, 'lead', 'campaign_lead_in')
    for run_index, mres in enumerate(MRES_VALUES, start=1):
        if run_index > 1:
            plan.separation(f'CONFIG_{run_index:02d}_MRES_{mres}', mres)
        config_t0 = plan.t
        plan.marker(f'CONFIG_{run_index:02d}_MRES_{mres}', mres)
        oscillation_label = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        plan.marker(oscillation_label, mres)
        plan.oscillation(mres)
        for mode in ('direct', 'individual'):
            for rate_name, rate in zip(
                TRAJECTORY_RATE_NAMES,
                TRAJECTORY_RATES_FULL_STEPS_S,
            ):
                label = trajectory_block_name(
                    mres,
                    mode,
                    rate_name,
                    rate,
                )
                plan.separation(label, mres)
                plan.marker(label, mres)
                plan.trajectory(mres, mode, rate_name, rate)
        plan.config_spans[mres] = (config_t0, plan.t)
    plan.dwell(CAMPAIGN_TAIL_S, 'lead', 'campaign_tail')
    if plan.position_mm != 0.0:
        raise AssertionError(
            f'Complete plan ended at {plan.position_mm:g} mm, not origin'
        )
    return plan


def latest_dry_run() -> Path | None:
    log_dir = HERE.parent / 'data' / 'hardware_runs'
    candidates = sorted(
        log_dir.glob('mres134_spreadcycle_microplyer_dry_run_*.csv')
    )
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

    expected_trajectories = (
        len(MRES_VALUES) * 2 * len(TRAJECTORY_RATES_FULL_STEPS_S)
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

    expected_separations = (
        expected_trajectories + len(MRES_VALUES) - 1
    )
    separation_starts = sum(
        row['event'] == 'SEPARATION_START' for row in rows
    )
    if separation_starts != expected_separations:
        raise SystemExit(
            f'Dry-run has {separation_starts} separations; '
            f'expected {expected_separations}'
        )

    endpoint_dwells = [
        row for row in rows
        if row['event'] == 'DWELL_START'
        and row['label'] in {'outbound_endpoint', 'return_endpoint'}
    ]
    if len(endpoint_dwells) != 2 * expected_trajectories:
        raise SystemExit(
            f'Dry-run has {len(endpoint_dwells)} endpoint dwells; '
            f'expected {2 * expected_trajectories}'
        )
    expected_detail = f'{TRAJECTORY_ENDPOINT_DWELL_S:.9g} s'
    if any(row['detail'] != expected_detail for row in endpoint_dwells):
        raise SystemExit(
            f'Not every endpoint dwell is {expected_detail}'
        )

    for mres in MRES_VALUES:
        block = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        moves = [
            row for row in rows
            if row['event'] == 'MOVE_ACK' and row['block'] == block
        ]
        if len(moves) != 2 * OSCILLATION_CYCLES:
            raise SystemExit(
                f'{block} has {len(moves)} moves; '
                f'expected {2 * OSCILLATION_CYCLES}'
            )
        final_position_rev = float(moves[-1]['ideal_position_rev'])
        if abs(final_position_rev) > 1e-12:
            raise SystemExit(
                f'{block} ends at {final_position_rev:g} rev'
            )


def plot_oscillation_axis(
    axis: plt.Axes,
    plan: Plan,
    mres: int,
) -> None:
    block = plan.oscillation_blocks[mres]
    relevant = [
        segment for segment in plan.segments
        if segment.kind == 'oscillation' and segment.mres == mres
    ]
    for segment in relevant:
        axis.plot(
            [segment.t0 - block.t0, segment.t1 - block.t0],
            [segment.p0_mm * 1000.0, segment.p1_mm * 1000.0],
            color=COLORS['oscillation'],
            lw=1.35,
        )
    start_um = relevant[0].p0_mm * 1000.0
    end_um = relevant[-1].p1_mm * 1000.0
    axis.scatter(
        [0.0, block.t1 - block.t0],
        [start_um, end_um],
        s=28,
        color='#222222',
        zorder=5,
        label='start = end',
    )
    amplitude_um = FULL_STEP_MM * 1000.0 / mres
    duration_s = block.t1 - block.t0
    # A small left margin keeps the near-instantaneous first step-up (right
    # at t=0) visible instead of flush against the axis spine.
    axis.set_xlim(-0.03 * duration_s, duration_s)
    axis.set_ylim(
        start_um - 0.10 * amplitude_um,
        start_um + 1.22 * amplitude_um,
    )
    axis.set_title(
        f'MRES {mres}: 15 complete one-microstep cycles '
        f'({amplitude_um:g} um; net displacement = 0)',
        fontsize=9.5,
    )
    axis.set_xlabel('Time within oscillation block (s)')
    axis.set_ylabel('Displacement (um)')
    axis.grid(True, alpha=0.25)
    axis.legend(loc='upper right', fontsize=7.2, frameon=False)


def shade_separations(
    axis: plt.Axes,
    blocks: list[Block],
    *,
    origin_s: float = 0.0,
) -> None:
    for block in blocks:
        if block.kind != 'separation':
            continue
        axis.axvspan(
            (block.t0 - origin_s) / 60.0,
            (block.t1 - origin_s) / 60.0,
            color=COLORS['separation'],
            alpha=0.50,
            lw=0,
        )


def plot_complete_trajectory(axis: plt.Axes, plan: Plan) -> None:
    for config_index, (mres, span) in enumerate(plan.config_spans.items()):
        t0, t1 = span
        if config_index % 2 == 0:
            axis.axvspan(
                t0 / 60.0,
                t1 / 60.0,
                color='#e8eef4',
                alpha=0.55,
                lw=0,
            )
        axis.axvline(t0 / 60.0, color='#5d6872', lw=0.7, alpha=0.55)
        axis.text(
            (t0 + t1) / 120.0,
            26.2,
            f'MRES {mres}',
            ha='center',
            va='bottom',
            fontsize=10,
            weight='bold',
            color='#34424e',
        )
    shade_separations(axis, plan.blocks)
    for segment in plan.segments:
        axis.plot(
            [segment.t0 / 60.0, segment.t1 / 60.0],
            [segment.p0_mm, segment.p1_mm],
            color=COLORS[segment.kind],
            lw=1.45 if segment.kind in {'direct', 'individual'} else 1.0,
            alpha=0.95,
        )
    axis.axhline(0.0, color='#3a3a3a', lw=0.6)
    axis.set_xlim(0.0, plan.t / 60.0)
    axis.set_ylim(-1.65, 27.4)
    axis.set_xlabel('Measurement time (min)')
    axis.set_ylabel('Ideal stage position (mm)')
    axis.set_title('Complete experimental trajectory')
    axis.grid(True, alpha=0.22)


LANE_ORDER = (
    'separation',
    'marker',
    'oscillation',
    'direct',
    'individual',
)
LANE_LABELS = (
    '10 s separation',
    'marker jump',
    '30 s oscillation',
    'direct 25 mm',
    'individual 25 mm',
)


def plot_schedule(
    axis: plt.Axes,
    blocks: list[Block],
    *,
    origin_s: float,
    duration_s: float,
    xlabel: str,
) -> None:
    lane_y = {kind: 0.04 + index * 0.18 for index, kind in enumerate(LANE_ORDER)}
    lane_h = 0.135
    for block in blocks:
        x0 = (block.t0 - origin_s) / 60.0
        width = (block.t1 - block.t0) / 60.0
        axis.add_patch(
            Rectangle(
                (x0, lane_y[block.kind]),
                width,
                lane_h,
                facecolor=COLORS[block.kind],
                edgecolor='none',
                alpha=0.90,
            )
        )
        if block.kind in {'direct', 'individual'}:
            short = block.label.split('_')[-2][0]
            axis.text(
                x0 + width / 2.0,
                lane_y[block.kind] + lane_h / 2.0,
                short,
                ha='center',
                va='center',
                color='white',
                fontsize=7,
                weight='bold',
            )
    axis.set_yticks(
        [lane_y[kind] + lane_h / 2.0 for kind in LANE_ORDER],
        LANE_LABELS,
    )
    axis.set_xlim(0.0, duration_s / 60.0)
    axis.set_ylim(0.0, 0.97)
    axis.set_xlabel(xlabel)
    axis.grid(axis='x', alpha=0.2)


def plot_full_plan(plan: Plan, out_path: Path) -> None:
    fig = plt.figure(figsize=(18, 11.5), constrained_layout=True)
    # Only 3 MRES values here (vs. 4 in the baseline), so the oscillation
    # panels fit in a single row of 6 columns (2 columns per panel) instead
    # of the baseline's 2x2 layout.
    grid = fig.add_gridspec(
        3,
        6,
        height_ratios=(1.55, 3.8, 0.95),
    )
    oscillation_axes = [
        fig.add_subplot(grid[0, index * 2:index * 2 + 2])
        for index in range(len(MRES_VALUES))
    ]
    overview = fig.add_subplot(grid[1, :])
    schedule = fig.add_subplot(grid[2, :], sharex=overview)

    for axis, mres in zip(oscillation_axes, MRES_VALUES):
        plot_oscillation_axis(axis, plan, mres)
    plot_complete_trajectory(overview, plan)
    plot_schedule(
        schedule,
        plan.blocks,
        origin_s=0.0,
        duration_s=plan.t,
        xlabel='Overall measurement time (min)',
    )

    legend = [
        Patch(
            facecolor=COLORS['separation'],
            label=f'{EXPERIMENT_SEPARATION_S:g} s experiment separation',
        ),
        Patch(
            facecolor=COLORS['marker'],
            label='unique separation marker',
        ),
        Patch(
            facecolor=COLORS['oscillation'],
            label='one-microstep oscillation',
        ),
        Patch(
            facecolor=COLORS['direct'],
            label='direct distance command',
        ),
        Patch(
            facecolor=COLORS['individual'],
            label='one full-step command at a time',
        ),
    ]
    overview.legend(
        handles=legend,
        loc='lower right',
        ncols=3,
        fontsize=8.2,
    )
    total_min = plan.t / 60.0
    fig.suptitle(
        'v4 SpreadCycle + MicroPlyer MRES measurement plan - '
        f'{total_min:.2f} min (55 min recording limit)\n'
        'MRES 1, 4, 16 (32 dropped to fit MicroPlyer); every microstep '
        'oscillation returns to its start; 25 mm out-and-back at 27.5, 70, '
        'and 200 full steps/s\n'
        'SpreadCycle enabled with MicroPlyer step interpolation ON for the '
        f'entire sequence; every velocity run dwells '
        f'{TRAJECTORY_ENDPOINT_DWELL_S:g} s at +25 mm and '
        f'{TRAJECTORY_ENDPOINT_DWELL_S:g} s after returning; '
        f'{EXPERIMENT_SEPARATION_S:g} s gaps retain unique marker jumps',
        fontsize=12,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            'Render the planned v4 SpreadCycle + MicroPlyer MRES '
            'trajectory campaign (MRES 1/4/16).'
        )
    )
    parser.add_argument(
        '--log-file',
        type=Path,
        help=(
            'Dry-run CSV to validate; defaults to the newest matching '
            'dry run. Validation is skipped with a warning if none exists '
            'yet -- this sketch has not been compiled/flashed/run.'
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    log_path = args.log_file or latest_dry_run()
    if log_path is None:
        print(
            'No SpreadCycle + MicroPlyer dry-run log found; rendering the '
            'plan model only (unvalidated against a captured run).'
        )
    else:
        validate_dry_run(log_path)
        print(f'Validated dry-run log: {log_path}')

    plan = build_plan()
    if plan.t >= 55.0 * 60.0:
        raise SystemExit(
            f'Planned sequence is {plan.t / 60.0:.2f} min, '
            'exceeding 55 min.'
        )
    print(
        f'Complete planned duration: {plan.t:.1f} s '
        f'({plan.t / 60.0:.2f} min); '
        f'{55.0 - plan.t / 60.0:.2f} min under the 55 min recording limit.'
    )

    out_dir = HERE.parent / 'rendered_assets'
    campaign_path = out_dir / 'planned_mres134_spreadcycle_microplyer_campaign.png'
    plot_full_plan(plan, campaign_path)
    print(f'Saved: {campaign_path}')


if __name__ == '__main__':
    main()
