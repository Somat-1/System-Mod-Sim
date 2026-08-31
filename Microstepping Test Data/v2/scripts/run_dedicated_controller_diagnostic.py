#!/usr/bin/env python3
'''Run an approximately 92-second preflight of the dedicated-controller moves.'''

from __future__ import annotations

import signal
import sys
from datetime import datetime
from decimal import Decimal
from pathlib import Path

import run_identification_dedicated_controller as controller


DIAGNOSTIC_MRES = 4
DIAGNOSTIC_CURRENT_PEAK_MA = 400
DIAGNOSTIC_CREEP_RECORD_S = 5.0
DIAGNOSTIC_SLOW_RATE_FS_S = 1.25
DIAGNOSTIC_SLOW_DURATION_S = 8.0
DIAGNOSTIC_FAST_RATE_FS_S = 200.0
DIAGNOSTIC_FAST_DURATION_S = 5.0


class DiagnosticRunner(controller.IdentificationRunner):
    def run_short_creep(self) -> None:
        with self.block('DIAG_C'):
            self.configure_speed(
                controller.CONDITIONING_FULL_STEPS_S, constant_start=True
            )
            for approach, sign in (('positive', 1), ('negative', -1)):
                self.move_full_steps(
                    Decimal(-4 * sign), f'{approach}_approach'
                )
                self.dwell(1.0, f'{approach}_pre_return')
                self.move_full_steps(
                    Decimal(4 * sign), f'{approach}_return'
                )
                self.dwell(
                    DIAGNOSTIC_CREEP_RECORD_S,
                    f'{approach}_short_creep_record',
                )
        self.assert_origin('diagnostic C')

    def run_diagnostic(self) -> None:
        self.context.run_index = 1
        self.context.mres = DIAGNOSTIC_MRES
        self.context.current = 'I_100pct'
        self.configure_mechanics(DIAGNOSTIC_MRES)
        self.configure_current(DIAGNOSTIC_CURRENT_PEAK_MA)
        self.command(f'ME {controller.AXIS}')
        self.log.log(
            'DIAGNOSTIC_START', **self.context.fields(),
            detail=(
                f'SC_peak_mA={DIAGNOSTIC_CURRENT_PEAK_MA}; '
                f'planned_motion_s=approximately_92'
            ),
        )

        self.run_marker('DIAG_CONFIG_MRES_4_I_100pct', 68)
        self.run_reference('DIAG_BLOCK_0_START')
        self.run_marker(
            'COND_C',
            controller.TEST_MARKER_AMPLITUDES_FULL_STEPS['COND_C'],
        )
        self.run_conditioning('DIAG_C')
        self.run_marker(
            'C', controller.TEST_MARKER_AMPLITUDES_FULL_STEPS['C']
        )
        self.run_short_creep()
        self.run_marker(
            'COND_D',
            controller.TEST_MARKER_AMPLITUDES_FULL_STEPS['COND_D'],
        )
        self.run_conditioning('DIAG_D')

        slow_label = f'D_{DIAGNOSTIC_SLOW_RATE_FS_S:g}'
        self.run_marker(
            slow_label,
            controller.TEST_MARKER_AMPLITUDES_FULL_STEPS[slow_label],
        )
        with self.block(f'DIAG_{slow_label}'):
            self.command(
                f'SS {controller.AXIS} 10 10 '
                f'{controller.PLATEAU_ACCEL_CODE} '
                f'{controller.PLATEAU_ACCEL_CODE} '
                f'{controller.PLATEAU_RAMP_TYPE}'
            )
            self.run_slow_plateau_direction(
                DIAGNOSTIC_SLOW_RATE_FS_S, 1,
                duration_s=DIAGNOSTIC_SLOW_DURATION_S,
            )
            self.run_slow_plateau_direction(
                DIAGNOSTIC_SLOW_RATE_FS_S, -1,
                duration_s=DIAGNOSTIC_SLOW_DURATION_S,
            )
        self.assert_origin('diagnostic slow D')

        fast_label = f'D_{DIAGNOSTIC_FAST_RATE_FS_S:g}'
        self.run_marker(
            fast_label,
            controller.TEST_MARKER_AMPLITUDES_FULL_STEPS[fast_label],
        )
        with self.block(f'DIAG_{fast_label}'):
            self.run_supported_plateau_direction(
                DIAGNOSTIC_FAST_RATE_FS_S, 1,
                duration_s=DIAGNOSTIC_FAST_DURATION_S,
            )
            self.run_supported_plateau_direction(
                DIAGNOSTIC_FAST_RATE_FS_S, -1,
                duration_s=DIAGNOSTIC_FAST_DURATION_S,
            )
        self.assert_origin('diagnostic fast D')

        self.run_marker(
            'BLOCK_0_END',
            controller.TEST_MARKER_AMPLITUDES_FULL_STEPS['BLOCK_0_END'],
        )
        self.run_reference('DIAG_BLOCK_0_END')
        self.assert_origin('diagnostic end')
        self.log.log('DIAGNOSTIC_COMPLETE', **self.context.fields())


def default_log_path(dry_run: bool) -> Path:
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    suffix = 'dry_run' if dry_run else 'live'
    return (
        controller.DEFAULT_LOG_DIR
        / f'dedicated_controller_diagnostic_{suffix}_{stamp}.csv'
    )


def apply_rig_defaults(args: object) -> None:
    args.direction = 0
    args.reuse_working_origin = False
    args.use_current_position_as_origin = True
    args.working_position_rev = None
    args.positive_limit_rev = '30'
    args.negative_limit_rev = '30'
    args.skip_home = True
    if args.execute:
        args.confirm_position_units = 'REVOLUTIONS'
        args.wait_for_acquisition = True


def main() -> int:
    parser = controller.build_parser()
    parser.description = (
        'Execute the short motion diagnostic with the EVO controller.'
    )
    args = parser.parse_args()
    apply_rig_defaults(args)
    controller.validate_args(args)
    clock = controller.Clock(args.dry_run)
    log_path = args.log_file or default_log_path(args.dry_run)
    log = controller.CsvEventLog(log_path, clock)
    context = controller.RunContext()
    transport: controller.Transport
    if args.dry_run:
        transport = controller.DryRunTransport(log, context)
    else:
        transport = controller.SerialTransport(
            args.port, args.command_timeout_s, log, context
        )

    runner = DiagnosticRunner(
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
        runner.initialise_session(args)
        controller.wait_for_acquisition(args, log, context)
        runner.run_diagnostic()
    except controller.Cancelled:
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
