#!/usr/bin/env python3
"""Run the v4 MRES/trajectory campaign on the EVO dedicated controller.

The measured campaign contains, for each MRES in (1, 4, 16, 32):

* 15 one-microstep out/return oscillations in 30 seconds;
* a 25 mm out/return trajectory at slow, moderate, and fast rates;
* the same trajectory generated once as a direct distance command and once
  as 2,500 individually issued full-step commands per leg.

The individual-command preflight runs before IDS acquisition.  It measures
the closed-loop host/controller command throughput and refuses to begin the
recorded campaign if the fastest requested rate cannot be sustained.  No
StealthChop, SpreadCycle, or StallGuard setting is changed by this script.
"""

from __future__ import annotations

import argparse
import signal
import sys
from dataclasses import dataclass
from datetime import datetime
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_settling_dedicated_controller import (  # noqa: E402
    AXIS,
    Clock,
    CsvEventLog,
    Cancelled,
    DEFAULT_LOG_DIR,
    DryRunTransport,
    MARKER_RATE_FULL_STEPS_S,
    MARKER_REVERSE_DWELL_S,
    MARKER_SETTLE_S,
    RunContext,
    SerialTransport,
    SettlingRunner,
    Transport,
    build_parser as build_base_parser,
    validate_args,
    wait_for_acquisition,
)


MRES_VALUES = (1, 4, 16, 32)
CURRENT_NAME = 'I_100pct'
CURRENT_PEAK_MA = 400

LEAD_MM_PER_REV = Decimal('2')
FULL_STEPS_PER_MM = Decimal('100')
TRAJECTORY_DISTANCE_MM = Decimal('25')
TRAJECTORY_FULL_STEPS = int(TRAJECTORY_DISTANCE_MM * FULL_STEPS_PER_MM)

# Selected from the v3 plateau rates.  The slower v3 rates cannot fit a
# complete 25 mm out-and-back campaign under the 55 minute recording limit.
TRAJECTORY_RATES_FULL_STEPS_S = (27.5, 70.0, 200.0)
TRAJECTORY_RATE_NAMES = ('slow', 'moderate', 'fast')

OSCILLATION_CYCLES = 15
OSCILLATION_HALF_DWELL_S = 1.0
OSCILLATION_DURATION_S = (
    OSCILLATION_CYCLES * 2 * OSCILLATION_HALF_DWELL_S
)
OSCILLATION_MOVE_RATE_FULL_STEPS_S = 250.0

TRAJECTORY_ENDPOINT_DWELL_S = 1.0
CAMPAIGN_LEAD_IN_S = 2.0
CAMPAIGN_TAIL_S = 2.0

THROUGHPUT_CANDIDATES_HZ = (
    27.5, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0,
    750.0, 1000.0,
)
THROUGHPUT_COMMANDS_PER_RATE = 20
THROUGHPUT_MIN_RATE_RATIO = 0.95
THROUGHPUT_MAX_LATE_PERIODS = 0.50

RECORDING_LIMIT_S = 55.0 * 60.0
PLANNED_MARGIN_S = 2.0 * 60.0


def rate_token(rate: float) -> str:
    return f'{rate:g}'.replace('.', 'p')


def trajectory_block_name(mres: int, mode: str, rate_name: str, rate: float) -> str:
    return (
        f'TRAJECTORY_MRES_{mres}_{mode.upper()}_'
        f'{rate_name.upper()}_{rate_token(rate)}FSPS'
    )


def marker_amplitude(marker_index: int) -> int:
    """Return a small, globally unique negative marker amplitude."""
    return 8 + 4 * marker_index


def planned_campaign_duration_s() -> float:
    """Ideal motion time, including all measured markers and dwells."""
    markers_per_mres = 2 + 2 * len(TRAJECTORY_RATES_FULL_STEPS_S)
    marker_count = len(MRES_VALUES) * markers_per_mres
    marker_s = sum(
        2.0 * marker_amplitude(index) / MARKER_RATE_FULL_STEPS_S
        + MARKER_REVERSE_DWELL_S
        + MARKER_SETTLE_S
        for index in range(1, marker_count + 1)
    )
    oscillation_s = sum(
        OSCILLATION_DURATION_S
        + 2.0 * OSCILLATION_CYCLES
        / (mres * OSCILLATION_MOVE_RATE_FULL_STEPS_S)
        for mres in MRES_VALUES
    )
    trajectories_s = len(MRES_VALUES) * 2.0 * sum(
        2.0 * TRAJECTORY_FULL_STEPS / rate
        + 2.0 * TRAJECTORY_ENDPOINT_DWELL_S
        for rate in TRAJECTORY_RATES_FULL_STEPS_S
    )
    return (
        CAMPAIGN_LEAD_IN_S + marker_s + oscillation_s
        + trajectories_s + CAMPAIGN_TAIL_S
    )


@dataclass(frozen=True)
class ThroughputResult:
    requested_hz: float
    achieved_hz: float
    max_lateness_us: float
    supported: bool


class CampaignRunner(SettlingRunner):
    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, **kwargs)
        self.marker_index = 0
        self.preflight_results: list[ThroughputResult] = []

    def next_marker(self, label: str) -> None:
        self.marker_index += 1
        self.run_marker(label, marker_amplitude(self.marker_index))

    def _simulate_move_time(self, full_steps: int, rate: float) -> None:
        # Live wait_ready() already accounts for physical travel.  The dry-run
        # transport is instantaneous, so explicitly advance its virtual clock.
        if self.clock.dry_run:
            self.clock.sleep(abs(full_steps) / rate)

    def run_oscillation(self, mres: int) -> None:
        label = f'OSCILLATION_MRES_{mres}_15CYCLES_30S'
        self.next_marker(label)
        with self.block(label):
            self.configure_speed(
                OSCILLATION_MOVE_RATE_FULL_STEPS_S, constant_start=True
            )
            self.log.log(
                'OSCILLATION_START', **self.context.fields(), label=label,
                detail=(
                    f'cycles={OSCILLATION_CYCLES}; driver_pulses_per_leg=1; '
                    f'half_dwell_s={OSCILLATION_HALF_DWELL_S:g}; '
                    f'planned_duration_s={OSCILLATION_DURATION_S:g}'
                ),
            )
            for cycle in range(1, OSCILLATION_CYCLES + 1):
                self.move_driver_pulses(
                    1, f'cycle_{cycle:02d}_forward',
                    rate_full_steps_s=OSCILLATION_MOVE_RATE_FULL_STEPS_S,
                )
                self.dwell(
                    OSCILLATION_HALF_DWELL_S,
                    f'cycle_{cycle:02d}_forward_dwell',
                )
                self.move_driver_pulses(
                    -1, f'cycle_{cycle:02d}_return',
                    rate_full_steps_s=OSCILLATION_MOVE_RATE_FULL_STEPS_S,
                )
                self.dwell(
                    OSCILLATION_HALF_DWELL_S,
                    f'cycle_{cycle:02d}_return_dwell',
                )
            self.log.log(
                'OSCILLATION_END', **self.context.fields(), label=label
            )
        self.assert_origin(label)

    def _run_direct_leg(self, direction: int, rate: float, label: str) -> None:
        full_steps = direction * TRAJECTORY_FULL_STEPS
        self.move_full_steps(
            Decimal(full_steps), label, rate_full_steps_s=rate
        )
        self._simulate_move_time(full_steps, rate)

    def run_direct_trajectory(
        self, mres: int, rate_name: str, rate: float
    ) -> None:
        block_name = trajectory_block_name(mres, 'direct', rate_name, rate)
        self.next_marker(block_name)
        with self.block(block_name):
            self.configure_speed(rate, constant_start=True)
            self.log.log(
                'TRAJECTORY_START', **self.context.fields(), label='direct',
                rate_full_steps_s=rate,
                detail=(
                    f'distance_mm={TRAJECTORY_DISTANCE_MM}; '
                    f'full_steps_per_leg={TRAJECTORY_FULL_STEPS}; '
                    'commands_per_leg=1'
                ),
            )
            self._run_direct_leg(1, rate, 'direct_outbound_25mm')
            self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, 'outbound_endpoint')
            self._run_direct_leg(-1, rate, 'direct_return_25mm')
            self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, 'return_endpoint')
            self.log.log(
                'TRAJECTORY_END', **self.context.fields(), label='direct',
                rate_full_steps_s=rate,
            )
        self.assert_origin(block_name)

    def _run_individual_leg(
        self, direction: int, rate: float, label: str
    ) -> tuple[float, float]:
        period_ns = int(round(1e9 / rate))
        start_ns = self.clock.now_ns()
        max_lateness_us = 0.0
        for pulse_index in range(1, TRAJECTORY_FULL_STEPS + 1):
            deadline_ns = start_ns + pulse_index * period_ns
            self.clock.sleep_until_ns(deadline_ns)
            ack_ns = self.move_full_steps(
                Decimal(direction), label,
                rate_full_steps_s=rate, deadline_ns=deadline_ns,
            )
            lateness_us = max(0.0, (ack_ns - deadline_ns) / 1000.0)
            max_lateness_us = max(max_lateness_us, lateness_us)
            if pulse_index in (1, TRAJECTORY_FULL_STEPS):
                self.log.log(
                    'INDIVIDUAL_LEG_BOUNDARY', **self.context.fields(),
                    label=label, rate_full_steps_s=rate,
                    pulse_index=pulse_index, lateness_us=lateness_us,
                )
        elapsed_s = (self.clock.now_ns() - start_ns) / 1e9
        achieved_hz = TRAJECTORY_FULL_STEPS / elapsed_s
        return achieved_hz, max_lateness_us

    def run_individual_trajectory(
        self, mres: int, rate_name: str, rate: float
    ) -> None:
        block_name = trajectory_block_name(mres, 'individual', rate_name, rate)
        self.next_marker(block_name)
        with self.block(block_name):
            self.configure_speed(rate, constant_start=True)
            self.log.log(
                'TRAJECTORY_START', **self.context.fields(), label='individual',
                rate_full_steps_s=rate,
                detail=(
                    f'distance_mm={TRAJECTORY_DISTANCE_MM}; '
                    f'full_steps_per_leg={TRAJECTORY_FULL_STEPS}; '
                    f'commands_per_leg={TRAJECTORY_FULL_STEPS}; '
                    f'driver_pulses_per_command={mres}'
                ),
            )
            outbound_hz, outbound_late_us = self._run_individual_leg(
                1, rate, 'individual_outbound_25mm'
            )
            self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, 'outbound_endpoint')
            return_hz, return_late_us = self._run_individual_leg(
                -1, rate, 'individual_return_25mm'
            )
            self.dwell(TRAJECTORY_ENDPOINT_DWELL_S, 'return_endpoint')
            self.log.log(
                'INDIVIDUAL_RATE_SUMMARY', **self.context.fields(),
                label=block_name, rate_full_steps_s=rate,
                lateness_us=max(outbound_late_us, return_late_us),
                detail=(
                    f'outbound_achieved_hz={outbound_hz:.6g}; '
                    f'return_achieved_hz={return_hz:.6g}; '
                    f'commands_per_leg={TRAJECTORY_FULL_STEPS}'
                ),
            )
            self.log.log(
                'TRAJECTORY_END', **self.context.fields(), label='individual',
                rate_full_steps_s=rate,
            )
        self.assert_origin(block_name)

    def benchmark_individual_rate(self, rate: float) -> ThroughputResult:
        self.configure_speed(rate, constant_start=True)
        period_ns = int(round(1e9 / rate))
        start_ns = self.clock.now_ns()
        max_lateness_us = 0.0
        for command_index in range(1, THROUGHPUT_COMMANDS_PER_RATE + 1):
            deadline_ns = start_ns + command_index * period_ns
            self.clock.sleep_until_ns(deadline_ns)
            direction = 1 if command_index % 2 else -1
            self.move_full_steps(
                Decimal(direction), f'preflight_{rate_token(rate)}',
                rate_full_steps_s=rate, deadline_ns=deadline_ns,
            )
            lateness_us = max(
                0.0, (self.clock.now_ns() - deadline_ns) / 1000.0
            )
            max_lateness_us = max(max_lateness_us, lateness_us)
        elapsed_s = (self.clock.now_ns() - start_ns) / 1e9
        achieved_hz = THROUGHPUT_COMMANDS_PER_RATE / elapsed_s
        allowed_lateness_us = (
            THROUGHPUT_MAX_LATE_PERIODS * period_ns / 1000.0
        )
        supported = (
            achieved_hz >= rate * THROUGHPUT_MIN_RATE_RATIO
            and max_lateness_us <= allowed_lateness_us
        )
        result = ThroughputResult(
            requested_hz=rate,
            achieved_hz=achieved_hz,
            max_lateness_us=max_lateness_us,
            supported=supported,
        )
        self.log.log(
            'THROUGHPUT_RESULT', **self.context.fields(),
            label='individual_full_step_commands',
            rate_full_steps_s=rate, lateness_us=max_lateness_us,
            detail=(
                f'achieved_hz={achieved_hz:.6g}; '
                f'supported={int(supported)}; '
                f'samples={THROUGHPUT_COMMANDS_PER_RATE}; '
                f'criterion_rate_ratio={THROUGHPUT_MIN_RATE_RATIO:g}; '
                f'criterion_max_late_periods={THROUGHPUT_MAX_LATE_PERIODS:g}'
            ),
        )
        self.assert_origin(f'throughput {rate:g} Hz')
        return result

    def run_throughput_preflight(self) -> list[ThroughputResult]:
        self.context.block = 'THROUGHPUT_PREFLIGHT_UNRECORDED'
        self.context.run_index = 0
        self.context.mres = MRES_VALUES[-1]
        self.context.current = CURRENT_NAME
        self.configure_mechanics(self.context.mres)
        self.configure_current(CURRENT_PEAK_MA)
        self.command(f'ME {AXIS}')
        self.log.log(
            'THROUGHPUT_PREFLIGHT_START', **self.context.fields(),
            detail=(
                'unrecorded; alternating one-full-step MR commands; '
                f'driver_pulses_per_command={self.context.mres}'
            ),
        )
        self.preflight_results = []
        for rate in THROUGHPUT_CANDIDATES_HZ:
            result = self.benchmark_individual_rate(rate)
            self.preflight_results.append(result)
            if not result.supported:
                break
        supported = [r.requested_hz for r in self.preflight_results if r.supported]
        max_supported = max(supported, default=0.0)
        self.log.log(
            'THROUGHPUT_PREFLIGHT_END', **self.context.fields(),
            rate_full_steps_s=max_supported,
            detail=f'maximum_tested_supported_hz={max_supported:g}',
        )
        return self.preflight_results

    def require_campaign_throughput(self) -> None:
        if not self.preflight_results:
            return
        supported = [r.requested_hz for r in self.preflight_results if r.supported]
        max_supported = max(supported, default=0.0)
        required = max(TRAJECTORY_RATES_FULL_STEPS_S)
        if max_supported < required:
            raise RuntimeError(
                f'Individual-command preflight sustained only {max_supported:g} '
                f'full-step commands/s; campaign requires {required:g}. '
                'The recorded run was not started.'
            )

    def run_campaign(self) -> None:
        planned_s = planned_campaign_duration_s()
        self.log.log(
            'CAMPAIGN_START', **self.context.fields(),
            detail=(
                'v4 MRES campaign; baseline motion only; '
                'no chopper-mode or StallGuard reconfiguration; '
                f'ideal_planned_duration_s={planned_s:.6g}'
            ),
        )
        self.dwell(CAMPAIGN_LEAD_IN_S, 'campaign_lead_in')
        for run_index, mres in enumerate(MRES_VALUES, start=1):
            self.check_cancelled()
            self.assert_origin(f'before MRES {mres}')
            self.context.run_index = run_index
            self.context.mres = mres
            self.context.current = CURRENT_NAME
            self.configure_mechanics(mres)
            self.configure_current(CURRENT_PEAK_MA)
            self.command(f'ME {AXIS}')
            self.log.log(
                'RUN_CONFIG', **self.context.fields(),
                detail=(
                    f'SC_peak_mA={CURRENT_PEAK_MA}; '
                    'modes=controller_default_unchanged'
                ),
            )

            config_label = f'CONFIG_{run_index:02d}_MRES_{mres}'
            self.next_marker(config_label)
            self.run_oscillation(mres)
            for mode in ('direct', 'individual'):
                for rate_name, rate in zip(
                    TRAJECTORY_RATE_NAMES, TRAJECTORY_RATES_FULL_STEPS_S
                ):
                    if mode == 'direct':
                        self.run_direct_trajectory(mres, rate_name, rate)
                    else:
                        self.run_individual_trajectory(mres, rate_name, rate)
            self.assert_origin(f'after MRES {mres}')
            self.log.log('RUN_COMPLETE', **self.context.fields())
        self.dwell(CAMPAIGN_TAIL_S, 'campaign_tail')
        self.log.log('CAMPAIGN_COMPLETE', **self.context.fields())


def build_parser() -> argparse.ArgumentParser:
    parser = build_base_parser()
    parser.description = (
        'Execute the v4 MRES oscillation and 25 mm trajectory campaign.'
    )
    parser.add_argument(
        '--throughput-only', action='store_true',
        help='Run the individual-command throughput preflight, then stop.',
    )
    parser.add_argument(
        '--skip-throughput-preflight', action='store_true',
        help=(
            'Skip the live rate check. Not recommended: the campaign may exceed '
            'the recording limit if UART commands cannot sustain 200 Hz.'
        ),
    )
    return parser


def default_log_path(dry_run: bool) -> Path:
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    suffix = 'dry_run' if dry_run else 'live'
    return DEFAULT_LOG_DIR / f'mres_trajectory_{suffix}_{stamp}.csv'


def main() -> int:
    args = build_parser().parse_args()
    validate_args(args)
    clock = Clock(args.dry_run)
    log_path = args.log_file or default_log_path(args.dry_run)
    log = CsvEventLog(log_path, clock)
    context = RunContext()
    transport: Transport
    if args.dry_run:
        transport = DryRunTransport(log, context)
    else:
        transport = SerialTransport(
            args.port, args.command_timeout_s, log, context
        )

    runner = CampaignRunner(
        transport, log, clock, context, direction=args.direction,
        positive_limit_rev=Decimal(args.positive_limit_rev),
        negative_limit_rev=Decimal(args.negative_limit_rev),
        status_timeout_s=args.status_timeout_s,
        status_poll_s=args.status_poll_s,
    )

    def request_cancel(_signum: int, _frame: object) -> None:
        runner.cancel()

    signal.signal(signal.SIGINT, request_cancel)
    if hasattr(signal, 'SIGTERM'):
        signal.signal(signal.SIGTERM, request_cancel)

    exit_code = 0
    reason = 'NORMAL'
    try:
        planned_s = planned_campaign_duration_s()
        if planned_s > RECORDING_LIMIT_S - PLANNED_MARGIN_S:
            raise RuntimeError(
                f'Ideal campaign duration {planned_s / 60.0:.2f} min exceeds '
                f'the {RECORDING_LIMIT_S / 60.0:.0f} min recording limit '
                f'after the {PLANNED_MARGIN_S / 60.0:.0f} min safety margin.'
            )
        runner.initialise_session(args)
        if not args.skip_throughput_preflight:
            runner.run_throughput_preflight()
            runner.require_campaign_throughput()
        if not args.throughput_only:
            wait_for_acquisition(args, log, context)
            runner.run_campaign()
    except Cancelled:
        reason = 'CANCELLED'
        exit_code = 130
    except Exception as exc:
        reason = f'FAILED: {type(exc).__name__}: {exc}'
        print(reason, file=sys.stderr)
        log.log('ERROR', **context.fields(), detail=reason)
        exit_code = 1
    finally:
        runner.safe_shutdown(reason)
        transport.close()
        log.close()
        print(f'Log: {log_path}')
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
