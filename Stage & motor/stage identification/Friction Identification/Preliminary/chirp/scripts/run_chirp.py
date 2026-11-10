#!/usr/bin/env python3
"""Run a chirp profile on the board and record it.

    python scripts/run_chirp.py --port /dev/cu.usbmodem3101 --profile high
    python scripts/run_chirp.py --port ... --profile full --out <folder>

Checks that the board runs the chirp firmware and that its driver and encoder
are ready, prints CHECK, then sends RUN <profile> (the firmware counts down 7 s
before moving; Ctrl+C during the countdown or the run sends ABORT). Writes two
files into the output folder, both line-buffered so they can be followed live:

  chirp_<profile>_run.log       every line the board sent
  chirp_<profile>_windows.csv   the demodulated windows (# W lines), one row each

Columns of the windows file: t_s, f_hz, cycles, f_issued_hz, amp_target,
amp_issued, amp_enc, phase_deg, lag_max, dc_err. Amplitudes are microsteps at
MRES 256 (the encoder is converted at 2.56 microsteps/count).
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SILENCE_S = 45.0       # longest gap between board lines during a run
COLUMNS = ["t_s", "f_hz", "cycles", "f_issued_hz", "amp_target", "amp_issued",
           "amp_enc", "phase_deg", "lag_max", "dc_err"]


def read_line(board, deadline):
    while time.monotonic() < deadline:
        raw = board.readline()
        if raw:
            return raw.decode("utf-8", errors="replace").rstrip("\r\n")
    raise TimeoutError("board response timed out")


def ask(board, command, prefix, timeout=5.0):
    board.write((command + "\n").encode())
    board.flush()
    deadline = time.monotonic() + timeout
    while True:
        line = read_line(board, deadline)
        if line.startswith("# ERR"):
            raise RuntimeError(f"board rejected {command!r}: {line}")
        if line.startswith(prefix):
            return line


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", required=True)
    ap.add_argument("--profile", choices=("high", "full"), default="high")
    ap.add_argument("--current-ma", type=int,
                    help="coil current for this run, RMS mA (default: firmware's 500)")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "recordings" / dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
    args = ap.parse_args()
    import serial

    args.out.mkdir(parents=True, exist_ok=True)
    log_path = args.out / f"chirp_{args.profile}_run.log"
    win_path = args.out / f"chirp_{args.profile}_windows.csv"

    with serial.Serial(args.port, 115200, timeout=0.2) as board:
        time.sleep(1.5)
        board.reset_input_buffer()
        pong = ask(board, "PING", "# PONG", timeout=10)
        print(pong)
        if "fw=stage_chirp_sweep" not in pong:
            raise RuntimeError("the board is not running the chirp firmware; flash "
                               "firmware/stage_chirp_sweep first")
        if "ready=1" not in pong:
            raise RuntimeError("driver or encoder not ready (see # READY at boot)")
        if args.current_ma:
            print(ask(board, f"CURRENT {args.current_ma}", "# CURRENT_OK"))
        check = ask(board, "CHECK", "# CHECK")
        print(check)

        with log_path.open("w", encoding="utf-8", buffering=1) as log, \
                win_path.open("w", encoding="utf-8", buffering=1) as win:
            log.write(f"# host_start_utc={dt.datetime.now(dt.timezone.utc).isoformat()}\n")
            log.write(f"# profile={args.profile}\n{pong}\n{check}\n")
            win.write(",".join(COLUMNS) + "\n")
            board.write(f"RUN {args.profile.upper()}\n".encode())
            board.flush()
            print(f"RUN {args.profile.upper()} sent; 7 s countdown, then motion. "
                  f"Ctrl+C sends ABORT.\nRecording to {args.out}")
            try:
                while True:
                    # The pre-roll is 30 s with nothing to report; allow more.
                    try:
                        line = read_line(board, time.monotonic() + SILENCE_S)
                    except TimeoutError:
                        # Never leave the board running with nobody listening.
                        board.write(b"ABORT\n")
                        board.flush()
                        raise
                    log.write(line + "\n")
                    if line.startswith("# W,"):
                        win.write(line[4:] + "\n")
                    elif line.startswith(("# SEG", "# COUNTDOWN", "# GAP_JUMP",
                                          "# TOP_1KHZ", "# DRV_STATUS", "# WIGGLE_DONE",
                                          "# RUN_", "# ERR", "# WARNING")):
                        print(line)
                    # RUN_START also begins "# RUN_", so name the endings exactly.
                    if line.startswith(("# RUN_COMPLETE", "# RUN_FAILED", "# RUN_ABORTED", "# ABORTED")):
                        return 0 if line.startswith("# RUN_COMPLETE") else 1
            except KeyboardInterrupt:
                board.write(b"ABORT\n")
                board.flush()
                print("\nABORT sent; waiting for the board ...")
                deadline = time.monotonic() + 5
                while time.monotonic() < deadline:
                    try:
                        line = read_line(board, deadline)
                    except TimeoutError:
                        break
                    log.write(line + "\n")
                    if line.startswith(("# RUN_COMPLETE", "# RUN_FAILED", "# RUN_ABORTED", "# ABORTED")):
                        print(line)
                        break
                return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, TimeoutError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
