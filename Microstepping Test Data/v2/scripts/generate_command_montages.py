#!/usr/bin/env python3
"""Generate command-only motion-sequence montages for the v2 stage test."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / 'rendered_assets'
PLOT_DIR = ASSET_DIR / 'trajectory_visualization_plots'
DATA_DIR = ROOT / 'data'
SUMMARY = ASSET_DIR / 'command_montage_summary.json'
CONFIG = DATA_DIR / 'motion_sequence_config.json'

FULL_STEP_TRAVEL_M = 5.0e-6
MRES = (('1/16', 16), ('1/4', 4), ('1/2', 2), ('1/1', 1))
BURST_FULL_STEPS_S = 250.0
CONDITIONING_FULL_STEPS_S = 150.0
DWELL_LADDER = 0.4
EXEC_LADDER_REPEATS = 25
PREVIEW_LADDER_REPEATS = 3
EXEC_LOOP_REPEATS = 10
PREVIEW_LOOP_REPEATS = 2
EXEC_DOUBLET_REPEATS = 20
PREVIEW_DOUBLET_REPEATS = 3
LOOP_DWELL_S = 0.30
DOUBLET_DWELL_S = 1.0
PLATEAU_RAMP_S = 0.5
PLATEAU_RATES_FS_S = (0.125, 0.375, 1.25, 3.5, 9.5, 27.5, 70.0, 200.0)

REFERENCE = (16, -16, 4, -4, 1, -1, -16, 16, -4, 4, -1, 1)
N_VALUES = (1, 2, 4, 8, 16, 32)
NEST_DESC = (32, -16, 8, -4, 2, -1, 1, -2, 4, -8, 16, -32)
NEST_ASYM = (8, -3, 2, -5, 6, -8, 4, -4)
NEST_MINOR = (64, -16, 2, -2, -16, 2, -2, -16, 2, -2, -16)


@dataclass
class Trace:
    mres: int
    time_s: list[float] = field(default_factory=lambda: [0.0])
    position_m: list[float] = field(default_factory=lambda: [0.0])

    @property
    def t(self) -> float:
        return self.time_s[-1]

    @property
    def x(self) -> float:
        return self.position_m[-1]

    @property
    def pulse_travel_m(self) -> float:
        return FULL_STEP_TRAVEL_M / self.mres

    def event(self, time_s: float, delta_pulses: int) -> None:
        self.time_s.append(float(time_s))
        self.position_m.append(
            self.position_m[-1] + delta_pulses * self.pulse_travel_m
        )

    def hold(self, duration_s: float) -> None:
        if duration_s < 0.0:
            raise ValueError('Hold duration must be non-negative')
        self.time_s.append(self.t + duration_s)
        self.position_m.append(self.x)

    def driver_move(self, pulses: int, pulse_rate_hz: float) -> None:
        if pulses == 0:
            return
        direction = 1 if pulses > 0 else -1
        start = self.t
        interval = 1.0 / pulse_rate_hz
        for index in range(abs(pulses)):
            self.event(start + (index + 1) * interval, direction)

    def burst(self, pulses: int) -> None:
        self.driver_move(pulses, BURST_FULL_STEPS_S * self.mres)

    def physical_move(self, full_steps: int, full_step_rate_hz: float) -> None:
        self.driver_move(
            full_steps * self.mres, full_step_rate_hz * self.mres
        )

    def linear_rate_ramp(
        self, direction: int, target_full_steps_s: float, ramp_in: bool,
    ) -> None:
        pulse_rate = target_full_steps_s * self.mres
        pulse_count = max(int(round(0.5 * pulse_rate * PLATEAU_RAMP_S)), 1)
        start = self.t
        for index in range(1, pulse_count + 1):
            fraction = index / pulse_count
            normalized_time = (
                np.sqrt(fraction) if ramp_in
                else 1.0 - np.sqrt(1.0 - fraction)
            )
            self.event(
                start + PLATEAU_RAMP_S * float(normalized_time), direction
            )

    def plateau_direction(
        self, direction: int, full_step_rate_hz: float,
    ) -> None:
        duration = plateau_duration(full_step_rate_hz)
        use_ramp = full_step_rate_hz >= 10.0
        if use_ramp:
            self.linear_rate_ramp(direction, full_step_rate_hz, True)
        pulse_rate = full_step_rate_hz * self.mres
        pulse_count = int(np.floor(pulse_rate * duration + 1.0e-12))
        start = self.t
        for index in range(pulse_count):
            self.event(start + (index + 1) / pulse_rate, direction)
        remainder = duration - pulse_count / pulse_rate
        if remainder > 0.0:
            self.hold(remainder)
        if use_ramp:
            self.linear_rate_ramp(direction, full_step_rate_hz, False)

    def arrays(self) -> tuple[np.ndarray, np.ndarray]:
        return np.asarray(self.time_s), np.asarray(self.position_m)


def plateau_duration(full_step_rate_hz: float) -> float:
    return float(np.clip(200.0 / full_step_rate_hz, 5.0, 20.0))


def reference_trace(mres: int) -> Trace:
    trace = Trace(mres)
    trace.hold(2.0)
    for pulses in REFERENCE:
        trace.burst(pulses)
        trace.hold(1.0)
    trace.hold(2.0)
    return trace


def conditioning_trace(mres: int) -> Trace:
    trace = Trace(mres)
    trace.physical_move(4, CONDITIONING_FULL_STEPS_S)
    trace.physical_move(-4, CONDITIONING_FULL_STEPS_S)
    trace.hold(2.0)
    return trace


def a1_trace(mres: int, n_value: int) -> Trace:
    trace = Trace(mres)
    for _ in range(PREVIEW_LADDER_REPEATS):
        trace.burst(n_value)
        trace.hold(DWELL_LADDER)
    for _ in range(PREVIEW_LADDER_REPEATS):
        trace.burst(-n_value)
        trace.hold(DWELL_LADDER)
    return trace


def a2_trace(mres: int, n_value: int) -> Trace:
    trace = Trace(mres)
    for _ in range(PREVIEW_LADDER_REPEATS):
        trace.burst(n_value)
        trace.hold(DWELL_LADDER)
        trace.burst(-n_value)
        trace.hold(DWELL_LADDER)
    return trace


def loop_trace(mres: int, pattern: tuple[int, ...]) -> Trace:
    trace = Trace(mres)
    for _ in range(PREVIEW_LOOP_REPEATS):
        for pulses in pattern:
            trace.burst(pulses)
            trace.hold(LOOP_DWELL_S)
    return trace


def creep_trace(mres: int, approach: str) -> Trace:
    trace = Trace(mres)
    sign = 1 if approach == 'positive' else -1
    trace.physical_move(-4 * sign, CONDITIONING_FULL_STEPS_S)
    trace.hold(1.0)
    trace.physical_move(4 * sign, CONDITIONING_FULL_STEPS_S)
    trace.hold(60.0)
    return trace


def plateau_trace(mres: int, full_step_rate_hz: float) -> Trace:
    trace = Trace(mres)
    trace.plateau_direction(1, full_step_rate_hz)
    trace.plateau_direction(-1, full_step_rate_hz)
    return trace


def doublet_trace(mres: int, n_value: int) -> Trace:
    trace = Trace(mres)
    for _ in range(PREVIEW_DOUBLET_REPEATS):
        trace.burst(n_value)
        trace.burst(-n_value)
        trace.hold(DOUBLET_DWELL_S)
    return trace


def validate_trace(trace: Trace, name: str) -> None:
    if not np.all(np.diff(np.asarray(trace.time_s)) >= 0.0):
        raise AssertionError(f'{name}: time is not monotonic')
    if not np.isclose(trace.x, 0.0, atol=1.0e-15):
        raise AssertionError(f'{name}: final position is {trace.x:.6e} m')


def trace_metadata(trace: Trace) -> dict[str, float | int]:
    _, position = trace.arrays()
    return {
        'duration_s': trace.t,
        'samples': len(trace.time_s),
        'minimum_command_m': float(np.min(position)),
        'maximum_command_m': float(np.max(position)),
        'final_command_m': trace.x,
    }


def render_montage(
    slug: str,
    group: str,
    title: str,
    row_specs: list[tuple[str, Callable[[int], Trace]]],
    execution_note: str,
) -> tuple[Path, dict[str, dict[str, dict[str, float | int]]]]:
    row_count = len(row_specs)
    fig, axes = plt.subplots(
        row_count, len(MRES),
        figsize=(15.5, 2.35 * row_count + 1.5),
        squeeze=False,
    )
    colors = ('#1565c0', '#00897b', '#ef6c00', '#c62828')
    metadata: dict[str, dict[str, dict[str, float | int]]] = {}

    for row_index, (row_label, factory) in enumerate(row_specs):
        traces = [factory(microsteps) for _, microsteps in MRES]
        for (mres_label, _), trace in zip(MRES, traces):
            validate_trace(trace, f'{slug}/{row_label}/{mres_label}')
        row_positions = [
            trace.arrays()[1] * 1.0e6 for trace in traces
        ]
        y_min = min(float(np.min(values)) for values in row_positions)
        y_max = max(float(np.max(values)) for values in row_positions)
        y_span = max(y_max - y_min, 0.1)
        y_limits = (y_min - 0.08 * y_span, y_max + 0.08 * y_span)
        x_max = max(trace.t for trace in traces)
        metadata[row_label] = {}

        for column, ((mres_label, _), trace, color) in enumerate(
            zip(MRES, traces, colors)
        ):
            axis = axes[row_index, column]
            time_values, position_values = trace.arrays()
            axis.step(
                time_values, position_values * 1.0e6,
                where='post', color=color, linewidth=1.0,
            )
            axis.set_xlim(0.0, x_max)
            axis.set_ylim(*y_limits)
            axis.axhline(0.0, color='#555555', linestyle=':', linewidth=0.55)
            axis.grid(True, color='#d0d0d0', linewidth=0.4)
            axis.set_xlabel('Local time (s)')
            if column == 0:
                axis.set_ylabel(f'{row_label}\nCommand travel (µm)')
            if row_index == 0:
                pulse_um = FULL_STEP_TRAVEL_M / MRES[column][1] * 1.0e6
                axis.set_title(
                    f'MRES {mres_label}\n1 pulse = {pulse_um:g} µm'
                )
            axis.text(
                0.98, 0.96, execution_note,
                transform=axis.transAxes, ha='right', va='top',
                fontsize=7.3, color='#37474f',
                bbox={
                    'boxstyle': 'round,pad=0.22',
                    'facecolor': 'white',
                    'edgecolor': '#b0bec5',
                    'alpha': 0.85,
                },
            )
            metadata[row_label][mres_label] = trace_metadata(trace)

    fig.suptitle(
        title + '\n'
        'Command-only preview; shared y-scale across each row; '
        'full-step travel = 5.0 µm',
        fontsize=13,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    output_dir = PLOT_DIR / group
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f'{slug}_command_montage.png'
    temporary = output.with_name(f'{output.stem}.tmp{output.suffix}')
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    replace_with_retry(temporary, output)
    return output, metadata


def replace_with_retry(
    temporary: Path, output: Path, attempts: int = 8,
) -> None:
    for attempt in range(attempts):
        try:
            temporary.replace(output)
            return
        except PermissionError:
            if attempt + 1 == attempts:
                raise
            time.sleep(0.15 * (attempt + 1))


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f'{path.stem}.tmp{path.suffix}')
    temporary.write_text(
        json.dumps(payload, indent=2) + '\n', encoding='utf-8'
    )
    replace_with_retry(temporary, path)


def motion_config() -> dict:
    return {
        'status': 'command preview; not executable controller instructions',
        'full_step_travel_m': FULL_STEP_TRAVEL_M,
        'mres': [{
            'label': label,
            'microsteps_per_full_step': microsteps,
            'u_per_driver_pulse': 16 // microsteps,
            'burst_pulse_rate_hz': BURST_FULL_STEPS_S * microsteps,
        } for label, microsteps in MRES],
        'outer_loop': {
            'holding_currents': ['I_lo', 'I_mid', 'I_hi'],
            'executions': 12,
            'executions_per_runner': 12,
            'hardware_acquisition_segments_total': 24,
            'stage_start': 'same fixed starting position for every execution',
        },
        'active_dedicated_controller_campaign': {
            'mres': [4, 2, 1],
            'current_levels': [
                {
                    'name': 'I_50pct',
                    'sc_peak_ma': 200,
                    'relative_percent': 50,
                },
                {
                    'name': 'I_100pct',
                    'sc_peak_ma': 400,
                    'relative_percent': 100,
                },
            ],
            'executions': 6,
            'planned_motion_duration_s': 2552.067,
            'configuration_marker_amplitude_full_steps': (
                '64 + 4 * one_based_run_index'
            ),
            'test_marker_amplitude_full_steps': {
                'COND_C': 12, 'COND_D': 16, 'C': 20,
                'D_0.125': 24, 'D_0.375': 28,
                'D_1.25': 32, 'D_3.5': 36, 'D_9.5': 40,
                'D_27.5': 44, 'D_70': 48, 'D_200': 52,
                'BLOCK_0_END': 56,
            },
            'marker_pattern': (
                'negative leap, 1.0 s dwell, equal positive return, '
                '0.5 s settle'
            ),
        },
        'current_levels': [
            {
                'name': 'I_lo',
                'tmc_set_rms_ma': 360,
                'tmc_measured_rms_ma': 355,
                'dedicated_controller_sc_peak_ma': 502,
                'hold_equals_run': True,
            },
            {
                'name': 'I_mid',
                'tmc_set_rms_ma': 600,
                'tmc_measured_rms_ma': 556,
                'dedicated_controller_sc_peak_ma': 786,
                'hold_equals_run': True,
            },
            {
                'name': 'I_hi',
                'tmc_set_rms_ma': 750,
                'tmc_measured_rms_ma': 715,
                'dedicated_controller_sc_peak_ma': 1011,
                'hold_equals_run': True,
            },
        ],
        'hardware_runners': {
            'esp32_s3_tmc2209': [
                'block_0', 'conditioning', 'block_A1', 'block_A2',
                'block_B', 'block_E', 'block_0_end',
            ],
            'dedicated_controller_axis_X': [
                'block_0', 'conditioning', 'block_C', 'block_D',
                'block_0_end',
            ],
        },
        'execution': {
            'dwell_ladder_s': DWELL_LADDER,
            'ladder_repeats': EXEC_LADDER_REPEATS,
            'loop_repeats': EXEC_LOOP_REPEATS,
            'doublet_repeats': EXEC_DOUBLET_REPEATS,
            'loop_dwell_s': LOOP_DWELL_S,
            'doublet_dwell_s': DOUBLET_DWELL_S,
            'conditioning_full_steps': 4,
            'conditioning_rate_full_steps_s': (
                CONDITIONING_FULL_STEPS_S
            ),
            'plateau_rates_full_steps_s': list(PLATEAU_RATES_FS_S),
            'plateau_duration': 'clip(200/f_fs, 5, 20) seconds',
            'plateau_ramp_s': PLATEAU_RAMP_S,
            'plateau_ramp_rule': 'enabled only for f_fs >= 10',
            'dedicated_controller_plateau_execution': {
                'software_paced_rates_full_steps_s': [0.125, 0.375, 1.25],
                'software_pacing_command': 'absolute MA microsteps',
                'timestamp': 'every command acknowledgement',
                'controller_paced_rates_full_steps_s': [
                    3.5, 9.5, 27.5, 70.0, 200.0,
                ],
                'accel_code': 628,
                'decel_code': 628,
                'ramp_type': 1,
                'identification_head_discard_s': 0.5,
            },
        },
        'preview': {
            'ladder_repeats': PREVIEW_LADDER_REPEATS,
            'loop_repeats': PREVIEW_LOOP_REPEATS,
            'doublet_repeats': PREVIEW_DOUBLET_REPEATS,
            'current_traces': 'not plotted; command is current-independent',
        },
        'reference_moves_driver_pulses': list(REFERENCE),
        'nested_patterns_driver_pulses': {
            'descending': list(NEST_DESC),
            'asymmetric': list(NEST_ASYM),
            'minor': list(NEST_MINOR),
        },
        'unresolved_for_hardware': [
            'pilot confirmation of dwell_ladder',
            'pilot confirmation of execution repeat count',
        ],
    }


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    details: dict[str, dict] = {}

    def create(
        key: str, group: str, title: str,
        rows: list[tuple[str, Callable[[int], Trace]]],
        note: str,
    ) -> None:
        path, metadata = render_montage(key, group, title, rows, note)
        outputs[key] = str(path.relative_to(ROOT))
        details[key] = metadata

    create(
        'block_0_reference', '00_reference_and_conditioning',
        'Block 0 — reference fingerprint',
        [('12-move fingerprint', reference_trace)],
        'start + end in execution',
    )
    create(
        'conditioning', '00_reference_and_conditioning',
        'Conditioning wipe before each major block',
        [('±4 full steps', conditioning_trace)],
        'execution: +4 fs, −4 fs, 2 s settle',
    )
    create(
        'block_A1_unidirectional', '01_step_ladders',
        'Block A1 — unidirectional step-and-settle',
        [(f'N={n}', lambda m, n=n: a1_trace(m, n)) for n in N_VALUES],
        'preview ×3; execution ×25',
    )
    create(
        'block_A2_alternating', '01_step_ladders',
        'Block A2 — alternating ±N',
        [(f'N={n}', lambda m, n=n: a2_trace(m, n)) for n in N_VALUES],
        'preview ×3; execution ×25',
    )
    create(
        'block_B_nested_loops', '02_reversal_and_creep',
        'Block B — nested reversal loops',
        [
            ('descending', lambda m: loop_trace(m, NEST_DESC)),
            ('asymmetric', lambda m: loop_trace(m, NEST_ASYM)),
            ('minor loop', lambda m: loop_trace(m, NEST_MINOR)),
        ],
        'preview ×2; execution ×10',
    )
    create(
        'block_C_creep', '02_reversal_and_creep',
        'Block C — creep approach and 60 s record',
        [
            ('arrive from +', lambda m: creep_trace(m, 'positive')),
            ('arrive from −', lambda m: creep_trace(m, 'negative')),
        ],
        'physical ±4-full-step approach',
    )
    create(
        'block_D_velocity_plateaus', '03_velocity_and_doublets',
        'Block D — velocity plateaus',
        [(
            f'{rate:g} fs/s; ramp={PLATEAU_RAMP_S if rate >= 10 else 0:g}s',
            lambda m, rate=rate: plateau_trace(m, rate),
        ) for rate in PLATEAU_RATES_FS_S],
        'both directions; duration clip(200/f, 5, 20)',
    )
    create(
        'block_E_doublets', '03_velocity_and_doublets',
        'Block E — immediate-return doublets',
        [(
            f'N={n}', lambda m, n=n: doublet_trace(m, n)
        ) for n in (1, 2, 4, 8, 16)],
        'preview ×3; execution ×20',
    )

    config = motion_config()
    write_json(CONFIG, config)
    write_json(SUMMARY, {
        'configuration': str(CONFIG.relative_to(ROOT)),
        'figures': outputs,
        'validation': {
            'all_preview_cells_return_to_zero': True,
            'row_y_scales_shared_across_mres': True,
            'current_dimension_omitted_from_command_preview': True,
        },
        'cell_metadata': details,
    })
    print(f'Wrote {len(outputs)} command montage figures')
    for output in outputs.values():
        print(f'  {output}')
    print(f'Wrote {CONFIG.relative_to(ROOT)}')
    print(f'Wrote {SUMMARY.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
