#!/usr/bin/env python3
"""Upload a generated motor-identification plan and record board telemetry."""

from __future__ import annotations
import argparse
import csv
import datetime as dt
import re
import sys
import time
from pathlib import Path

ACTIONS = {"hold": 0, "move": 1, "stair": 2,
           "step_response": 3, "velocity": 4, "reversal": 5}
MODE_REQUIRED = ("spreadcycle=1", "intpol=0", "stealth=0")


def load_plan(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("plan is empty")
    for i, row in enumerate(rows):
        if row["action"] not in ACTIONS:
            raise ValueError(f"row {i}: unknown action {row['action']}")
        mres = int(row["mres"])
        if not 1 <= mres <= 256 or mres & (mres - 1):
            raise ValueError(f"row {i}: invalid MRES {mres}")
        if not 10 <= int(row["sample_hz"]) <= 5000:
            raise ValueError(f"row {i}: sample rate outside 10..5000 Hz")
    if int(rows[-1]["end_u256"]):
        raise ValueError("plan does not end at command zero")
    return rows


def tag_for(index: int, row: dict) -> str:
    stem = re.sub(r"[^A-Za-z0-9_-]", "_", row["block"])
    return f"{index:02d}_{stem}"[:23]


def plan_command(index: int, row: dict) -> str:
    speed_milli = round(float(row["speed_fs_s"]) * 1000)
    return (f"PLAN ADD {ACTIONS[row['action']]} {int(row['microsteps'])} "
            f"{speed_milli} {int(row['dwell_ms'])} {int(row['repeats'])} "
            f"{int(row['current_ma'])} {int(row['mres'])} "
            f"{int(row['sample_hz'])} {tag_for(index, row)}")


def choose_port(explicit: str | None) -> str:
    if explicit:
        return explicit
    from serial.tools import list_ports
    ports = list(list_ports.comports())
    if not ports:
        raise RuntimeError("no serial ports found; pass --port")
    preferred = [p for p in ports if any(x in (p.description or "").lower()
                 for x in ("esp32", "usb jtag", "usb serial", "wch"))]
    if len(preferred) == 1:
        return preferred[0].device
    shown = ", ".join(p.device for p in (preferred or ports))
    raise RuntimeError(f"port is ambiguous ({shown}); pass --port")


def read_line(board, deadline: float) -> str:
    while time.monotonic() < deadline:
        raw = board.readline()
        if raw:
            return raw.decode("utf-8", errors="replace").rstrip("\r\n")
    raise TimeoutError("board response timed out")


def send_expect(board, command: str, prefix: str, timeout: float = 3.0) -> list[str]:
    board.write((command + "\n").encode())
    board.flush()
    seen = []
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        line = read_line(board, deadline)
        seen.append(line)
        if line.startswith("# ERR"):
            raise RuntimeError(f"board rejected {command!r}: {line}")
        if line.startswith(prefix):
            return seen
    raise TimeoutError(f"no {prefix!r} after {command!r}")


def verify_live_board(board) -> list[str]:
    """Query live state; do not depend on opening the port resetting USB CDC."""
    pong = send_expect(board, "PING", "# PONG", timeout=10)
    if "ready=1" not in pong[-1].lower():
        raise RuntimeError(f"firmware is not ready: {pong[-1]}")
    mode = send_expect(board, "DRIVER", "# DRIVER,", timeout=5)
    driver_line = mode[-1].lower()
    if not all(item in driver_line for item in MODE_REQUIRED):
        raise RuntimeError("driver readback did not confirm SpreadCycle=1, "
                           "INTPOL=0, and StealthChop inactive")
    return pong + mode


def run_board(port: str, rows: list[dict], output: Path,
              confirm: bool) -> None:
    try:
        import serial
    except ImportError as exc:
        raise RuntimeError("pyserial is required: python -m pip install pyserial") from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    with serial.Serial(port, 115200, timeout=0.2, write_timeout=2) as board:
        # Native USB may or may not reset on open. Wait for either case, then
        # query the live firmware and driver registers explicitly.
        time.sleep(1.5)
        board.reset_input_buffer()
        boot = verify_live_board(board)
        for line in boot:
            print(line)
        send_expect(board, "PLAN CLEAR", "# PLAN_CLEARED")
        for index, row in enumerate(rows):
            replies = send_expect(board, plan_command(index, row), "# PLAN_ADD")
            if index % 10 == 0 or index == len(rows) - 1:
                print(replies[-1])

        if confirm:
            answer = input("Motor decoupled/known-load, encoder on shaft, safe to rotate? [y/N] ")
            if answer.strip().lower() not in {"y", "yes"}:
                print("Not started; plan remains loaded.")
                return

        with output.open("w", encoding="utf-8", newline="\n") as log:
            log.write(f"# host_start_utc={dt.datetime.now(dt.timezone.utc).isoformat()}\n")
            log.write(f"# port={port}\n")
            log.write("record_type,seq,t_us,action,event,command_u256,encoder,"
                      "index_count,current_ma,mres\n")
            board.write(b"PLAN RUN\n")
            board.flush()
            print(f"Recording to {output}; Ctrl+C sends ABORT.")
            try:
                while True:
                    line = read_line(board, time.monotonic() + 5)
                    log.write(line + "\n")
                    if line.startswith("# ACTION") or line.startswith("# RUN_"):
                        print(line)
                    if line.startswith("# RUN_COMPLETE"):
                        break
                    if line.startswith("# ERR"):
                        raise RuntimeError(line)
            except KeyboardInterrupt:
                board.write(b"ABORT\n")
                board.flush()
                print("\nABORT sent; waiting for acknowledgement.")
                try:
                    line = read_line(board, time.monotonic() + 3)
                    log.write(line + "\n")
                    print(line)
                except TimeoutError:
                    print("No abort acknowledgement received.")


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path,
                        default=root / "sequence" / "motor_identification_plan.csv")
    parser.add_argument("--port")
    parser.add_argument("--out", type=Path,
                        default=root / "recordings" / f"motor_identification_{stamp}.csv")
    parser.add_argument("--yes", action="store_true",
                        help="skip the final physical-safety prompt")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate and print board commands without opening a port")
    args = parser.parse_args()
    rows = load_plan(args.plan)
    print(f"Validated {len(rows)} actions; estimated "
          f"{float(rows[-1]['t_end_s']) / 60:.2f} min.")
    if args.dry_run:
        for index, row in enumerate(rows):
            print(plan_command(index, row))
        return
    run_board(choose_port(args.port), rows, args.out, not args.yes)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, TimeoutError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
