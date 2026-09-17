#!/usr/bin/env python3
"""Reusable EVO dedicated-controller infrastructure for the v4 campaign.

Protocol: ../../v2/docs/Stepper Motor Controller Command list.pdf.
This module owns transport, logging, guarded motion primitives, markers,
origin checks, CLI safety validation, and shutdown handling. Campaign ordering
and experiment-specific constants live in `run_mres_trajectory_campaign.py`.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, ROUND_HALF_UP, getcontext
from pathlib import Path
from typing import Iterator, Optional, Union

try:
    from typing import Protocol
except ImportError:  # Python 3.7
    class Protocol:
        pass


getcontext().prec = 28
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOG_DIR = ROOT / 'data' / 'hardware_runs'

AXIS = 'X'
DEFAULT_SERIAL_PORT = 'COM5'
EXPECTED_MODULE_TYPE = 0
BAUD_RATE = 115200
MOTOR_FULL_STEPS_PER_REV = 200
POSITION_QUANTUM_REV = Decimal('0.000001')
TRIGGER_MASK = 32

# Generic motion defaults and marker signature used by the active campaign.
DEFAULT_ACCEL_CODE = 628
DEFAULT_RAMP_TYPE = 1
MARKER_RATE_FULL_STEPS_S = 150.0
MARKER_REVERSE_DWELL_S = 1.0
MARKER_SETTLE_S = 0.5


class Cancelled(RuntimeError):
    """Raised after a local interrupt or explicit safety cancellation."""


class Clock:
    def __init__(self, dry_run: bool) -> None:
        self.dry_run = dry_run
        self._virtual_ns = 0

    def now_ns(self) -> int:
        return self._virtual_ns if self.dry_run else time.perf_counter_ns()

    def sleep(self, duration_s: float) -> None:
        if duration_s < 0:
            raise ValueError('Sleep duration cannot be negative')
        if self.dry_run:
            self._virtual_ns += int(round(duration_s * 1e9))
        else:
            time.sleep(duration_s)

    def sleep_until_ns(self, deadline_ns: int) -> None:
        if self.dry_run:
            self._virtual_ns = max(self._virtual_ns, deadline_ns)
            return
        while True:
            remaining_ns = deadline_ns - time.perf_counter_ns()
            if remaining_ns <= 0:
                return
            if remaining_ns > 2_000_000:
                time.sleep((remaining_ns - 1_000_000) / 1e9)


class CsvEventLog:
    FIELDNAMES = (
        'utc', 'monotonic_ns', 'event', 'run_index', 'current', 'mres',
        'block', 'label', 'command', 'response', 'ideal_position_rev',
        'commanded_position_rev', 'rate_full_steps_s', 'pulse_index',
        'lateness_us', 'detail',
    )

    def __init__(self, path: Path, clock: Clock) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.clock = clock
        self._handle = path.open('w', newline='', encoding='utf-8')
        self._writer = csv.DictWriter(self._handle, fieldnames=self.FIELDNAMES)
        self._writer.writeheader()

    def close(self) -> None:
        self._handle.flush()
        self._handle.close()

    def log(self, event: str, **fields: object) -> None:
        row = {name: '' for name in self.FIELDNAMES}
        row.update(
            utc=datetime.now(timezone.utc).isoformat(timespec='microseconds'),
            monotonic_ns=self.clock.now_ns(),
            event=event,
        )
        row.update({key: value for key, value in fields.items() if key in row})
        self._writer.writerow(row)
        self._handle.flush()


@dataclass
class RunContext:
    run_index: int = 0
    current: str = ''
    mres: int = 0
    block: str = 'SESSION'

    def fields(self) -> dict[str, object]:
        return {
            'run_index': self.run_index,
            'current': self.current,
            'mres': self.mres,
            'block': self.block,
        }


class Transport(Protocol):
    def command(self, text: str, expected_prefix: Optional[str] = None) -> str:
        ...

    def close(self) -> None:
        ...


class SerialTransport:
    def __init__(
        self, port: str, timeout_s: float, log: CsvEventLog,
        context: RunContext,
    ) -> None:
        try:
            import serial
        except ImportError as exc:
            raise RuntimeError(
                'pyserial is required: python -m pip install pyserial'
            ) from exc
        self.serial = serial.Serial(
            port=port, baudrate=BAUD_RATE, bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE,
            timeout=0.1, write_timeout=timeout_s,
        )
        self.timeout_s = timeout_s
        self.log = log
        self.context = context
        self.serial.reset_input_buffer()

    def command(self, text: str, expected_prefix: Optional[str] = None) -> str:
        if '\r' in text or '\n' in text:
            raise ValueError('Controller command must not contain CR/LF')
        expected = expected_prefix if expected_prefix is not None else text
        self.serial.write((text.upper() + '\r').encode('ascii'))
        self.serial.flush()
        deadline = time.monotonic() + self.timeout_s
        received: list[str] = []
        while time.monotonic() < deadline:
            raw = self.serial.readline()
            if not raw:
                continue
            line = raw.decode('ascii', errors='replace').strip()
            if not line:
                continue
            received.append(line)
            if line.upper().startswith(expected.upper()):
                self.log.log(
                    'COMMAND_ACK', **self.context.fields(),
                    command=text, response=line,
                )
                return line
        joined = ' | '.join(received)
        raise TimeoutError(
            f'No {expected!r} response for {text!r}; '
            f'received {joined or "<nothing>"}'
        )

    def close(self) -> None:
        self.serial.close()


class DryRunTransport:
    def __init__(self, log: CsvEventLog, context: RunContext) -> None:
        self.log = log
        self.context = context
        self.position = Decimal('0')
        self.output_byte = 0

    def command(self, text: str, expected_prefix: Optional[str] = None) -> str:
        parts = text.upper().split()
        if parts[0] == 'DM':
            response = 'DM 0'
        elif parts[0] == 'DS':
            response = f'DS {AXIS} 1'
        elif parts[0] == 'DP':
            response = f'DP {AXIS} {format_position(self.position)}'
        elif parts[0] == 'SO':
            self.output_byte |= int(parts[-1])
            response = f'SO {self.output_byte}'
        elif parts[0] == 'CO':
            self.output_byte &= ~int(parts[-1])
            response = f'CO {self.output_byte}'
        else:
            response = text.upper()
            if parts[0] in {'MA', 'SP'}:
                self.position = Decimal(parts[-1])
            elif parts[0] == 'MR':
                self.position += Decimal(parts[-1])
        self.log.log(
            'DRY_COMMAND', **self.context.fields(),
            command=text, response=response,
        )
        return response

    def close(self) -> None:
        return


def format_position(value: Decimal) -> str:
    quantized = value.quantize(POSITION_QUANTUM_REV, rounding=ROUND_HALF_UP)
    if quantized == Decimal('-0.000000'):
        quantized = Decimal('0.000000')
    return f'{quantized:.6f}'


def speed_code(full_steps_s: float) -> int:
    omega_rad_s = 2.0 * math.pi * full_steps_s / MOTOR_FULL_STEPS_PER_REV
    return int(round(omega_rad_s * 100.0))


class DedicatedControllerRunner:
    def __init__(
        self, transport: Transport, log: CsvEventLog, clock: Clock,
        context: RunContext, *, direction: int,
        positive_limit_rev: Decimal, negative_limit_rev: Decimal,
        status_timeout_s: float, status_poll_s: float,
    ) -> None:
        self.transport = transport
        self.log = log
        self.clock = clock
        self.context = context
        self.direction = direction
        self.positive_limit_rev = positive_limit_rev
        self.negative_limit_rev = negative_limit_rev
        self.status_timeout_s = status_timeout_s
        self.status_poll_s = status_poll_s
        self.ideal_position_rev = Decimal('0')
        self.controller_origin_rev = Decimal('0')
        self.controller_position_pulses = 0
        self.relative_commands = False
        self.cancelled = False

    def cancel(self) -> None:
        self.cancelled = True

    def check_cancelled(self) -> None:
        if self.cancelled:
            raise Cancelled('Execution cancelled')

    def command(self, text: str, expected_prefix: Optional[str] = None) -> str:
        self.check_cancelled()
        return self.transport.command(text, expected_prefix)

    def wait_ready(self) -> None:
        deadline = time.monotonic() + self.status_timeout_s
        while True:
            self.check_cancelled()
            response = self.command(f'DS {AXIS}', f'DS {AXIS} ')
            try:
                status = int(response.split()[-1])
            except (ValueError, IndexError) as exc:
                raise RuntimeError(f'Invalid DS response: {response!r}') from exc
            if status == 1:
                return
            if status != 2:
                raise RuntimeError(f'Unexpected motor status: {response!r}')
            if not self.clock.dry_run and time.monotonic() >= deadline:
                raise TimeoutError('Motor did not become ready before timeout')
            self.clock.sleep(self.status_poll_s)

    def configure_mechanics(self, mres: int) -> None:
        expected = (
            f'SM {AXIS} {mres} {MOTOR_FULL_STEPS_PER_REV} {self.direction}'
        )
        response = self.command(expected)
        if response.upper() != expected:
            raise RuntimeError(f'SM verification failed: {response!r}')

    def configure_current(self, peak_ma: int) -> None:
        expected = f'SC {AXIS} {peak_ma} {peak_ma}'
        response = self.command(expected)
        if response.upper() != expected:
            raise RuntimeError(f'SC verification failed: {response!r}')

    def configure_speed(
        self, rate_full_steps_s: float, *, constant_start: bool,
        accel_code: int = DEFAULT_ACCEL_CODE,
        ramp_type: int = DEFAULT_RAMP_TYPE,
    ) -> int:
        maximum = speed_code(rate_full_steps_s)
        if not 10 <= maximum <= 12800:
            raise ValueError(
                f'{rate_full_steps_s:g} full steps/s maps to unsupported '
                f'MAXSPEED code {maximum}'
            )
        minimum = maximum if constant_start else 10
        self.command(
            f'SS {AXIS} {minimum} {maximum} '
            f'{accel_code} {accel_code} {ramp_type}'
        )
        return maximum

    def verify_target(self, target: Decimal) -> None:
        if target > self.positive_limit_rev:
            raise RuntimeError(
                f'Positive travel guard: {target} rev exceeds '
                f'{self.positive_limit_rev} rev'
            )
        if target < -self.negative_limit_rev:
            raise RuntimeError(
                f'Negative travel guard: {target} rev exceeds '
                f'-{self.negative_limit_rev} rev'
            )
        if self.relative_commands:
            controller_target = self.controller_origin_rev + target
            if abs(controller_target) > Decimal('30'):
                raise RuntimeError(
                    f'Relative move would take displayed controller position '
                    f'to {controller_target} rev, outside +/-30 rev'
                )

    def move_to_ideal(
        self, target: Decimal, label: str, *,
        rate_full_steps_s: Optional[float] = None,
        deadline_ns: Optional[int] = None,
    ) -> int:
        self.verify_target(target)
        commanded = format_position(target)
        if self.relative_commands:
            target_pulses = int(
                (target * MOTOR_FULL_STEPS_PER_REV * self.context.mres)
                .to_integral_value(rounding=ROUND_HALF_UP)
            )
            delta_pulses = target_pulses - self.controller_position_pulses
            command_text = f'MR {AXIS} {delta_pulses}'
            response = self.command(command_text, f'MR {AXIS} ')
            self.controller_position_pulses = target_pulses
        else:
            command_text = f'MA {AXIS} {commanded}'
            response = self.command(command_text)
        ack_ns = self.clock.now_ns()
        lateness_us: Union[str, float] = ''
        if deadline_ns is not None:
            lateness_us = (ack_ns - deadline_ns) / 1000.0
        self.ideal_position_rev = target
        self.log.log(
            'MOVE_ACK', **self.context.fields(), label=label,
            command=command_text, response=response,
            ideal_position_rev=str(target),
            commanded_position_rev=commanded,
            rate_full_steps_s=(
                '' if rate_full_steps_s is None else rate_full_steps_s
            ),
            lateness_us=lateness_us,
        )
        self.wait_ready()
        return ack_ns

    def move_delta_rev(
        self, delta_rev: Decimal, label: str, **kwargs: object,
    ) -> int:
        return self.move_to_ideal(
            self.ideal_position_rev + delta_rev, label, **kwargs
        )

    def move_driver_pulses(
        self, pulses: int, label: str, **kwargs: object,
    ) -> int:
        denominator = Decimal(MOTOR_FULL_STEPS_PER_REV * self.context.mres)
        return self.move_delta_rev(Decimal(pulses) / denominator, label, **kwargs)

    def move_full_steps(self, full_steps: Decimal, label: str, **kwargs: object) -> int:
        return self.move_delta_rev(
            full_steps / Decimal(MOTOR_FULL_STEPS_PER_REV), label, **kwargs
        )

    def dwell(self, duration_s: float, label: str) -> None:
        self.log.log(
            'DWELL_START', **self.context.fields(), label=label,
            detail=f'{duration_s:.9g} s',
        )
        self.clock.sleep(duration_s)
        self.log.log('DWELL_END', **self.context.fields(), label=label)

    def set_trigger(self, active: bool, *, safety: bool = False) -> None:
        opcode = 'SO' if active else 'CO'
        command = f'{opcode} {TRIGGER_MASK}'
        send = self.transport.command if safety else self.command
        response = send(command, f'{opcode} ')
        try:
            output_byte = int(response.split()[-1])
        except (ValueError, IndexError) as exc:
            raise RuntimeError(
                f'Invalid {opcode} response: {response!r}'
            ) from exc
        bit_is_set = bool(output_byte & TRIGGER_MASK)
        if bit_is_set != active:
            expected = 'set' if active else 'cleared'
            raise RuntimeError(
                f'Trigger bit was not {expected}: {response!r}'
            )

    @contextmanager
    def block(self, name: str) -> Iterator[None]:
        previous = self.context.block
        self.context.block = name
        self.set_trigger(True)
        self.log.log('BLOCK_START', **self.context.fields())
        body_failed = False
        try:
            yield
        except BaseException:
            body_failed = True
            raise
        finally:
            try:
                self.set_trigger(False, safety=True)
            except Exception as exc:
                self.log.log(
                    'TRIGGER_CLEAR_WARNING', **self.context.fields(),
                    detail=f'{type(exc).__name__}: {exc}',
                )
                if not body_failed:
                    raise
            finally:
                self.log.log('BLOCK_END', **self.context.fields())
                self.context.block = previous

    def assert_origin(self, where: str) -> None:
        if self.ideal_position_rev != 0:
            raise RuntimeError(
                f'{where} did not return to the ideal origin: '
                f'{self.ideal_position_rev} rev'
            )
        response = self.command(f'DP {AXIS}', f'DP {AXIS} ')
        try:
            displayed = Decimal(response.split()[-1])
        except (IndexError, ArithmeticError) as exc:
            raise RuntimeError(f'Invalid DP response: {response!r}') from exc
        expected_displayed = (
            self.controller_origin_rev
            if self.relative_commands else Decimal('0')
        )
        tolerance = Decimal('0') if self.relative_commands else POSITION_QUANTUM_REV
        if abs(displayed - expected_displayed) > tolerance:
            raise RuntimeError(
                f'{where} did not return to the captured controller '
                f'position: expected {expected_displayed}, '
                f'received {displayed} rev'
            )
        self.log.log(
            'ORIGIN_CHECK_OK', **self.context.fields(),
            label=where, response=response,
        )

    def run_marker(self, label: str, amplitude_full_steps: int) -> None:
        if amplitude_full_steps <= 0:
            raise ValueError('Marker amplitude must be positive')
        with self.block(f'MARKER_{label}'):
            self.configure_speed(
                MARKER_RATE_FULL_STEPS_S, constant_start=True
            )
            self.log.log(
                'MARKER_SIGNATURE', **self.context.fields(), label=label,
                detail=(
                    f'negative_then_positive; '
                    f'amplitude_full_steps={amplitude_full_steps}; '
                    f'reverse_dwell_s={MARKER_REVERSE_DWELL_S:g}'
                ),
            )
            self.move_full_steps(
                Decimal(-amplitude_full_steps), f'{label}_negative'
            )
            self.dwell(MARKER_REVERSE_DWELL_S, f'{label}_reverse_dwell')
            self.move_full_steps(
                Decimal(amplitude_full_steps), f'{label}_return'
            )
            self.dwell(MARKER_SETTLE_S, f'{label}_settle')
        self.assert_origin(f'marker {label}')

    def initialise_session(
        self,
        args: argparse.Namespace,
        *,
        initial_mres: int,
        current_name: str,
        current_peak_ma: int,
    ) -> None:
        module = self.command('DM', 'DM ')
        try:
            module_type = int(module.split()[-1])
        except (ValueError, IndexError) as exc:
            raise RuntimeError(f'Invalid DM response: {module!r}') from exc
        if module_type != EXPECTED_MODULE_TYPE:
            raise RuntimeError(
                f'Expected XY test box DM 0, received {module!r}'
            )

        self.context.mres = initial_mres
        self.context.current = current_name
        self.configure_mechanics(initial_mres)
        self.configure_current(current_peak_ma)
        self.command(f'ME {AXIS}')

        if not args.skip_home:
            self.command(
                f'SS {AXIS} {args.home_min_speed} {args.home_max_speed} '
                f'{args.home_accel} {args.home_decel} {args.home_ramp_type}'
            )
            self.command(
                f'MH {AXIS} {args.home_max_steps} {args.home_sensor_mask} '
                f'{args.home_overtravel}'
            )
            self.wait_ready()
            self.log.log('HOME_COMPLETE', **self.context.fields())

        # Same operating convention as v3: capture the current position
        # with DP and issue every move as a relative MR from there.
        self.wait_ready()
        response = self.command(f'DP {AXIS}', f'DP {AXIS} ')
        try:
            self.controller_origin_rev = Decimal(response.split()[-1])
        except (IndexError, ArithmeticError) as exc:
            raise RuntimeError(f'Invalid DP response: {response!r}') from exc
        self.relative_commands = True
        self.ideal_position_rev = Decimal('0')
        self.controller_position_pulses = 0
        self.log.log(
            'CURRENT_POSITION_CAPTURED_AS_ORIGIN',
            **self.context.fields(), response=response,
            commanded_position_rev=str(self.controller_origin_rev),
        )

    def safe_shutdown(self, reason: str) -> None:
        try:
            self.set_trigger(False, safety=True)
        except Exception as exc:
            self.log.log(
                'SHUTDOWN_WARNING', **self.context.fields(),
                detail=f'CO {TRIGGER_MASK} failed: {exc}',
            )
        try:
            self.transport.command(f'MO {AXIS}')
        except Exception as exc:
            self.log.log(
                'SHUTDOWN_WARNING', **self.context.fields(),
                detail=f'MO {AXIS} failed: {exc}',
            )
        self.log.log(
            'SAFE_SHUTDOWN', **self.context.fields(), detail=reason
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Execute a guarded v4 campaign on the EVO dedicated controller.'
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        '--dry-run', action='store_true',
        help='Validate all commands without opening a serial port.',
    )
    mode.add_argument(
        '--execute', action='store_true',
        help='Enable live communication and physical motor motion.',
    )
    parser.add_argument(
        '--port', default=DEFAULT_SERIAL_PORT,
        help=f'Serial port; defaults to {DEFAULT_SERIAL_PORT}.',
    )
    parser.add_argument(
        '--direction', type=int, choices=(0, 1), default=0,
        help='SM direction flag verified with the small direction test.',
    )
    parser.add_argument(
        '--positive-limit-rev', default='30',
        help='Available positive travel from the working origin.',
    )
    parser.add_argument(
        '--negative-limit-rev', default='30',
        help='Available negative travel magnitude from the working origin.',
    )
    parser.add_argument(
        '--skip-home', action='store_true',
        help='Only when this controller session has already been homed.',
    )
    parser.add_argument('--home-max-steps', type=int)
    parser.add_argument('--home-sensor-mask', type=int)
    parser.add_argument('--home-overtravel', type=int)
    parser.add_argument('--home-min-speed', type=int)
    parser.add_argument('--home-max-speed', type=int)
    parser.add_argument('--home-accel', type=int)
    parser.add_argument('--home-decel', type=int)
    parser.add_argument('--home-ramp-type', type=int, choices=(0, 1))
    parser.add_argument('--command-timeout-s', type=float, default=2.0)
    parser.add_argument('--status-timeout-s', type=float, default=120.0)
    parser.add_argument('--status-poll-s', type=float, default=0.01)
    parser.add_argument(
        '--wait-for-acquisition', action='store_true',
        help='After initialization, wait for Enter before starting the first marker.',
    )
    parser.add_argument(
        '--log-file', type=Path,
        help='CSV output path; defaults to data/hardware_runs.',
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.dry_run:
        if args.direction is None:
            args.direction = 0
        args.skip_home = True
        return

    args.skip_home = True  # this campaign always captures the current position as origin
    required = {
        '--port': args.port,
        '--direction': args.direction,
        '--positive-limit-rev': args.positive_limit_rev,
        '--negative-limit-rev': args.negative_limit_rev,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise SystemExit('Live execution requires: ' + ', '.join(missing))
    for label in ('positive_limit_rev', 'negative_limit_rev'):
        value = Decimal(getattr(args, label))
        if value <= 0 or value > Decimal('30'):
            option = '--' + label.replace('_', '-')
            raise SystemExit(f'{option} must be in (0, 30] rev')


def wait_for_acquisition(
    args: argparse.Namespace, log: CsvEventLog, context: RunContext,
) -> None:
    if args.dry_run or not args.wait_for_acquisition:
        return
    print(
        'Controller initialized at the working origin. Start IDS acquisition, '
        'then press Enter to begin motion.'
    )
    try:
        input()
    except EOFError as exc:
        raise RuntimeError(
            'Acquisition confirmation requires an interactive terminal'
        ) from exc
    log.log('ACQUISITION_CONFIRMED', **context.fields())
