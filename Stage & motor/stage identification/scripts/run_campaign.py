#!/usr/bin/env python3
"""Flash, validate, upload, run, and record the stage-identification campaign.

The default action never flashes and never moves: it connects to an already
flashed stage-identification board, performs read-back checks, uploads the
generated plan with an acknowledgement for every row, and asks for a physical
safety confirmation before ZERO and PLAN RUN.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from collections import deque
from pathlib import Path

BAUD = 115200
FQBN = "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc"
MAX_MOVES = 512
DEFAULT_LIMIT = 1100
MAX_PULSE_RATE_HZ = 20_000
REQUIRED_MODE = ("spreadcycle=1", "intpol=0", "stealth=0")
REQUIRED_COLUMNS = {
    "phase", "label", "t_start", "t_end", "start", "target", "speed",
    "dwell_s", "mres", "micro_steps", "micro_dwell_s", "group", "note",
}

ROOT = Path(__file__).resolve().parent.parent
SKETCH_DIR = ROOT / "firmware" / "stage_identification"
DEFAULT_PLAN = ROOT / "sequence" / "sequence_plan.csv"


def as_int(text: str, name: str, row_number: int) -> int:
    value = float(text)
    if not value.is_integer():
        raise ValueError(f"row {row_number}: {name} must be an integer, got {text!r}")
    return int(value)


def load_plan(path: Path, limit: int = DEFAULT_LIMIT) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"plan is missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError("plan is empty")
    if len(rows) > MAX_MOVES:
        raise ValueError(f"plan has {len(rows)} rows; firmware limit is {MAX_MOVES}")

    previous = 0
    previous_end = 0.0
    for index, row in enumerate(rows):
        number = index + 2  # CSV header is line 1
        start = as_int(row["start"], "start", number)
        target = as_int(row["target"], "target", number)
        speed = as_int(row["speed"], "speed", number)
        mres = as_int(row["mres"], "mres", number)
        micro_steps = as_int(row["micro_steps"], "micro_steps", number)
        dwell_ms = round(float(row["dwell_s"]) * 1000)
        micro_dwell_ms = round(float(row["micro_dwell_s"]) * 1000)
        t_start = float(row["t_start"])
        t_end = float(row["t_end"])

        if start != previous:
            raise ValueError(
                f"row {number}: start {start} does not match previous target {previous}")
        if abs(target) > limit:
            raise ValueError(f"row {number}: target {target} exceeds +/-{limit} soft limit")
        if speed <= 0:
            raise ValueError(f"row {number}: speed must be positive")
        if not 1 <= mres <= 256 or mres & (mres - 1):
            raise ValueError(f"row {number}: MRES {mres} is not a power of two in 1..256")
        if speed * mres > MAX_PULSE_RATE_HZ:
            raise ValueError(f"row {number}: pulse rate exceeds {MAX_PULSE_RATE_HZ} Hz")
        if min(dwell_ms, micro_dwell_ms, micro_steps) < 0:
            raise ValueError(f"row {number}: dwell and microstep fields must be non-negative")
        if micro_steps:
            expected = abs(target - start) * mres
            if micro_steps != expected:
                raise ValueError(
                    f"row {number}: staircase has {micro_steps} microsteps, expected {expected}")
            if micro_dwell_ms == 0:
                raise ValueError(f"row {number}: staircase dwell must be positive")
        elif micro_dwell_ms:
            raise ValueError(f"row {number}: micro_dwell_s set on a normal move")
        if abs(t_start - previous_end) > 0.001:
            raise ValueError(f"row {number}: timeline is not contiguous")
        if t_end <= t_start:
            raise ValueError(f"row {number}: t_end must be after t_start")

        row["_target"] = target
        row["_speed"] = speed
        row["_dwell_ms"] = dwell_ms
        row["_mres"] = mres
        row["_micro_steps"] = micro_steps
        row["_micro_dwell_ms"] = micro_dwell_ms
        previous = target
        previous_end = t_end

    if previous != 0:
        raise ValueError(f"plan ends at {previous}, not command zero")
    return rows


def tag_for(index: int, row: dict) -> str:
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{row['phase']}_{row['label']}").strip("_")
    return f"{index:03d}_{stem}"[:19]


def plan_command(index: int, row: dict) -> str:
    return (f"PLAN ADD {row['_target']} {row['_speed']} {row['_dwell_ms']} "
            f"{row['_mres']} {row['_micro_steps']} {row['_micro_dwell_ms']} "
            f"{tag_for(index, row)}")


def plan_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def choose_port(explicit: str | None) -> str:
    if explicit:
        return explicit
    try:
        from serial.tools import list_ports
    except ImportError as exc:
        raise RuntimeError("pyserial is required: python3 -m pip install pyserial") from exc
    ports = list(list_ports.comports())
    if not ports:
        raise RuntimeError("no serial ports found; pass --port")
    preferred = [p for p in ports if any(
        token in f"{p.description} {p.manufacturer}".lower()
        for token in ("esp32", "usb jtag", "usb serial", "wch", "cp210"))]
    candidates = preferred or ports
    if len(candidates) != 1:
        shown = ", ".join(f"{p.device} ({p.description})" for p in candidates)
        raise RuntimeError(f"serial port is ambiguous: {shown}; pass --port")
    return candidates[0].device


def read_line(board, deadline: float) -> str:
    while time.monotonic() < deadline:
        raw = board.readline()
        if raw:
            return raw.decode("utf-8", errors="replace").rstrip("\r\n")
    raise TimeoutError("board response timed out")


def send_expect(board, command: str, prefix: str,
                timeout: float = 5.0) -> list[str]:
    board.write((command + "\n").encode())
    board.flush()
    seen = []
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        line = read_line(board, deadline)
        seen.append(line)
        if line.startswith("# ERR") or line.startswith("# CHECK_FAILED"):
            raise RuntimeError(f"board rejected {command!r}: {line}")
        if line.startswith(prefix):
            return seen
    raise TimeoutError(f"no {prefix!r} response after {command!r}")


def verify_board(board) -> list[str]:
    lines = send_expect(board, "PING", "# PONG", timeout=10)
    pong = lines[-1].lower()
    if "fw=stage_identification" not in pong:
        raise RuntimeError(f"wrong firmware: {lines[-1]}")
    if "ready=1" not in pong or "driver=1" not in pong or "encoder=1" not in pong:
        raise RuntimeError(f"board is not ready: {lines[-1]}")
    mode = send_expect(board, "DRIVER", "# DRIVER", timeout=5)
    driver = mode[-1].lower()
    if not all(item in driver for item in REQUIRED_MODE):
        raise RuntimeError(
            "driver readback did not confirm SpreadCycle=1, INTPOL=0, and StealthChop off")
    status = send_expect(board, "STATUS", "# STATUS", timeout=3)
    return lines + mode + status


def compile_and_flash(port: str) -> None:
    cli = shutil.which("arduino-cli")
    if not cli:
        raise RuntimeError("arduino-cli is required for --flash")
    if subprocess.run([cli, "core", "list"], capture_output=True, text=True).stdout.find(
            "esp32:esp32") < 0:
        raise RuntimeError("ESP32 Arduino core is not installed")
    libs = subprocess.run([cli, "lib", "list"], capture_output=True, text=True)
    if libs.returncode or "TMCStepper" not in libs.stdout:
        raise RuntimeError("TMCStepper is not installed: arduino-cli lib install TMCStepper")
    print("Compiling stage-identification firmware from source ...")
    subprocess.run([cli, "compile", "--fqbn", FQBN, str(SKETCH_DIR)], check=True)
    print(f"Flashing {port} ...")
    subprocess.run(
        [cli, "upload", "--fqbn", FQBN, "--port", port, str(SKETCH_DIR)], check=True)
    print("Flash complete; waiting for USB to re-enumerate ...")
    time.sleep(3)


class SimulatedBoard:
    """Fast protocol model used by --simulate; it never touches hardware."""

    def __init__(self):
        self.responses = deque()
        self.moves: list[str] = []
        self.limit = DEFAULT_LIMIT

    def write(self, raw: bytes) -> int:
        for command in raw.decode().splitlines():
            self._command(command.strip())
        return len(raw)

    def _command(self, command: str) -> None:
        if command == "PING":
            self.responses.append(
                f"# PONG,fw=stage_identification,ready=1,driver=1,encoder=1,moves={len(self.moves)}")
        elif command == "DRIVER":
            self.responses.append("# DRIVER,spreadcycle=1,intpol=0,stealth=0")
        elif command == "STATUS":
            self.responses.append(
                f"# STATUS,pos=0,enc=0,moves={len(self.moves)},running=0,limit={self.limit}")
        elif command == "PLAN CLEAR":
            self.moves.clear()
            self.responses.append("# PLAN_CLEARED")
        elif command.startswith("LIMIT "):
            self.limit = int(command.split()[1])
            self.responses.append(f"# LIMIT_OK,limit={self.limit}")
        elif command.startswith("PLAN ADD "):
            self.moves.append(command)
            tag = command.split()[-1]
            index = len(self.moves) - 1
            self.responses.append(f"# PLAN_ADD,index={index},count={len(self.moves)},tag={tag}")
        elif command == "ZERO":
            self.responses.append("# ZERO_OK,pos=0,enc=0")
        elif command == "CHECK":
            self.responses.append(
                "# CHECK,command_steps=10,forward_delta=1000,return_delta=5,direction=1")
        elif command == "PLAN RUN":
            self.responses.append(f"# RUN_START,moves={len(self.moves)},pos=0,enc=0")
            for index, move in enumerate(self.moves):
                target = move.split()[2]
                tag = move.split()[-1]
                self.responses.append(f"# MOVE,{index},{tag},target={target}")
                self.responses.append(f"# ARRIVED,{index},pos={target},enc={int(target) * 200}")
            self.responses.append("# RUN_COMPLETE,pos=0,enc=0")
        elif command == "ABORT":
            self.responses.append("# ABORT")
            self.responses.append("# RUN_ABORTED,pos=0,enc=0")
        else:
            self.responses.append(f"# ERR,reason=unknown_command,command={command}")

    def readline(self) -> bytes:
        return ((self.responses.popleft() + "\n").encode()
                if self.responses else b"")

    def flush(self) -> None:
        pass

    def reset_input_buffer(self) -> None:
        self.responses.clear()

    def close(self) -> None:
        pass


def open_board(port: str):
    try:
        import serial
    except ImportError as exc:
        raise RuntimeError("pyserial is required: python3 -m pip install pyserial") from exc
    return serial.Serial(port, BAUD, timeout=0.2, write_timeout=2)


def upload_plan(board, rows: list[dict], limit: int) -> None:
    print(send_expect(board, f"LIMIT {limit}", "# LIMIT_OK")[-1])
    print(send_expect(board, "PLAN CLEAR", "# PLAN_CLEARED")[-1])
    for index, row in enumerate(rows):
        reply = send_expect(board, plan_command(index, row), "# PLAN_ADD")[-1]
        expected = f"index={index},count={index + 1},tag={tag_for(index, row)}"
        if expected not in reply:
            raise RuntimeError(f"upload acknowledgement mismatch at row {index}: {reply}")
        if index % 25 == 0 or index == len(rows) - 1:
            print(f"Uploaded {index + 1}/{len(rows)} moves")
    status = send_expect(board, "STATUS", "# STATUS")[-1]
    if f"moves={len(rows)}" not in status:
        raise RuntimeError(f"board did not retain the complete plan: {status}")


def record_run(board, rows: list[dict], plan_path: Path, output: Path,
               blind: bool = False) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    arrived = 0
    with output.open("w", encoding="utf-8", newline="\n", buffering=1) as log:
        metadata = {
            "host_start_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "plan": str(plan_path.resolve()),
            "plan_sha256": plan_sha256(plan_path),
            "moves": len(rows),
            "blind": blind,
        }
        log.write("# " + json.dumps(metadata, sort_keys=True) + "\n")
        for line in send_expect(board, "ZERO", "# ZERO_OK"):
            log.write(line + "\n")
        if blind:
            # No encoder attached: the preflight would see zero counts and
            # refuse, and the run records the command only. Say so in the log.
            print("Blind run: encoder preflight skipped; the log holds commanded "
                  "position only (encoder columns will read 0).")
        else:
            check = send_expect(board, "CHECK", "# CHECK,", timeout=10)
            for line in check:
                log.write(line + "\n")
            match = re.search(r"forward_delta=(-?\d+),return_delta=(-?\d+)", check[-1])
            if not match:
                raise RuntimeError(f"could not parse preflight encoder result: {check[-1]}")
            forward_delta, return_delta = map(int, match.groups())
            if not 200 <= abs(forward_delta) <= 5000:
                raise RuntimeError(
                    f"preflight encoder response is implausible: forward_delta={forward_delta}"
                    " (no encoder attached? use --blind)")
            if abs(return_delta) > max(300, abs(forward_delta) // 2):
                raise RuntimeError(
                    f"preflight did not return near its start: return_delta={return_delta}")
            print(f"Preflight passed: encoder {forward_delta:+d} counts over +0.1 mm; "
                  f"return error {return_delta:+d} counts.")
        # CHECK returns the commanded position to zero. Re-zero measured state
        # so the campaign starts from a fresh, unambiguous encoder origin.
        for line in send_expect(board, "ZERO", "# ZERO_OK"):
            log.write(line + "\n")
        start = send_expect(board, "PLAN RUN", "# RUN_START")
        for line in start:
            log.write(line + "\n")
        print(f"Running and recording to {output}. Ctrl+C sends ABORT.")
        try:
            while True:
                line = read_line(board, time.monotonic() + 10)
                log.write(line + "\n")
                if line.startswith("# ARRIVED"):
                    arrived += 1
                if line.startswith("# MOVE") or line.startswith("# RUN_"):
                    print(line)
                if line.startswith("# ERR") or line.startswith("# RUN_FAILED"):
                    raise RuntimeError(line)
                if line.startswith("# RUN_ABORTED"):
                    raise RuntimeError("campaign aborted")
                if line.startswith("# RUN_COMPLETE"):
                    if arrived != len(rows):
                        raise RuntimeError(
                            f"run completed with {arrived} arrivals; expected {len(rows)}")
                    if "pos=0" not in line:
                        raise RuntimeError(f"run did not finish at command zero: {line}")
                    print(f"Completed all {arrived} moves; log closed cleanly.")
                    return
        except KeyboardInterrupt:
            board.write(b"ABORT\n")
            board.flush()
            print("\nABORT sent; waiting for the board to stop ...")
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                line = read_line(board, deadline)
                log.write(line + "\n")
                if line.startswith("# RUN_ABORTED") or line.startswith("# ABORT"):
                    print(line)
                    return


def safety_confirmation(port: str, rows: list[dict], limit: int,
                        blind: bool = False) -> bool:
    print("\nThe plan is uploaded but motion has NOT started.")
    print(f"It will ZERO the present position, then use {len(rows)} moves within +/-{limit/100:.1f} mm.")
    print("Confirm: warm-up stopped; stage centred; STEP/DIR/EN on GPIO7/8/9;")
    print("BLIND - no encoder, nothing will be measured;" if blind else
          "encoder A/B/Z on GPIO4/5/6;")
    print("travel clear; motor power on; and emergency power removal within reach.")
    answer = input("Type RUN to zero and start: ")
    return answer.strip().upper() in ("RUN", f"RUN {port}".upper())


def main() -> int:
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--port", help="serial port, e.g. /dev/cu.usbmodem1101 or COM11")
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT,
                        help="soft limit in full steps about ZERO (default 1100 = 11 mm)")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "recordings" / f"stage_identification_{stamp}.log")
    parser.add_argument("--flash", action="store_true",
                        help="compile and flash the identification firmware first")
    parser.add_argument("--yes", action="store_true",
                        help="skip the typed physical-safety confirmation")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate the plan and commands without opening a serial port")
    parser.add_argument("--simulate", action="store_true",
                        help="exercise upload, run, and logging against a fast protocol model")
    parser.add_argument("--blind", action="store_true",
                        help="no encoder attached: skip the encoder preflight and run "
                             "the motion only")
    args = parser.parse_args()

    if not 1050 <= args.limit <= 2200:
        raise ValueError("limit must be 1050..2200 full steps")
    rows = load_plan(args.plan, args.limit)
    print(f"Validated {len(rows)} moves; planned duration "
          f"{float(rows[-1]['t_end']) / 60:.2f} min; max target "
          f"{max(abs(row['_target']) for row in rows) / 100:.2f} mm.")
    if args.dry_run:
        for index, row in enumerate(rows):
            print(plan_command(index, row))
        return 0

    if args.simulate:
        port = "SIMULATED"
        board = SimulatedBoard()
    else:
        port = choose_port(args.port)
        if args.flash:
            compile_and_flash(port)
        board = open_board(port)

    try:
        if not args.simulate:
            time.sleep(1.5)
            board.reset_input_buffer()
        for line in verify_board(board):
            print(line)
        upload_plan(board, rows, args.limit)
        if not args.yes and not args.simulate and not safety_confirmation(
                port, rows, args.limit, args.blind):
            print("Not started; the uploaded plan remains motionless on the board.")
            return 1
        record_run(board, rows, args.plan, args.out, blind=args.blind)
    finally:
        board.close()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TimeoutError, ValueError,
            subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
