#!/usr/bin/env python3
"""Run the MRES test from a terminal menu.

    python scripts/run_mres.py                      # menu
    python scripts/run_mres.py --mode stepping      # straight to a mode
    python scripts/run_mres.py --mode stepping --simulate --yes

Modes:

  1  warm-up        part 1: homing to the end switches + warm-up (not defined yet)
  2  stepping       part 2: x_ref baseline, six MRES blocks, end reference
  3  full session   part 1 then part 2
  4  single block   B1-B3 at one MRES, for bring-up or a repeat
  5  preview        the plan and the exact commands, no hardware
  6  simulate       mode 2 against a protocol model, no hardware
  7  regenerate     rebuild the plan CSV, the summary and the figure

The plan is built in memory from `mres_design.py`, never read back from a
CSV, so what runs is what the design says. Each segment (reference, each
block, end) is uploaded and run as its own plan; every one starts and ends at
x_start. Nothing moves before the typed RUN confirmation, and Ctrl+C sends
ABORT and stops the session.

The board protocol is the one specified in `../firmware/README.md`. No sketch
implements it yet, so modes 2-4 will stop at the firmware check until one does;
mode 6 exercises everything else.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import subprocess
import sys
import time
from collections import deque
from dataclasses import asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mres_design as md  # noqa: E402

ROOT = HERE.parent
DATA_DIR = ROOT / "data"
BAUD = 115200
FIRMWARE_ID = "fw=mres_test"
MAX_MOVES = 512
LINE_TIMEOUT_S = 10.0             # the board streams samples at 50 Hz while running

MODES = [
    ("1", "warmup", "Warm-up", "part 1: homing to the end switches + warm-up"),
    ("2", "stepping", "Stepping sequence", "part 2: x_ref baseline, six MRES blocks, end reference"),
    ("3", "full", "Full session", "part 1 then part 2"),
    ("4", "block", "Single block", "B1-B3 at one MRES, for bring-up or a repeat"),
    ("5", "preview", "Preview", "the plan and the exact commands, no hardware"),
    ("6", "simulate", "Simulate", "the stepping sequence against a protocol model"),
    ("7", "regenerate", "Regenerate", "rebuild the plan CSV, summary and figure"),
]


# ------------------------------------------------------------------ plan

def segments(rows):
    """Split a plan into its upload units, in order."""
    out = []
    for r in rows:
        if not out or out[-1][0] != r.segment:
            out.append((r.segment, []))
        out[-1][1].append(r)
    return out


def segment_file(name: str, rows) -> str:
    """S1_B2_MRES16-style names from the procedure, with one session."""
    if name.startswith("block"):
        return f"B{rows[0].block}_MRES{rows[0].mres}.log"
    return f"{name.upper()}.log"


def tag_for(index: int, r) -> str:
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{r.phase}_{r.tier or r.kind}").strip("_")
    return f"{index:03d}_{stem}"[:19]


def plan_command(index: int, r) -> str:
    return (f"PLAN ADD {r.kind} {r.target} {r.speed_fs} {round(r.dwell_s * 1000)} "
            f"{r.mres} {r.n_inc} {r.inc_usteps} {round(r.inc_dwell_s * 1000)} "
            f"{tag_for(index, r)}")


def validate(rows) -> None:
    """The checks the board will make, made first so a bad plan never uploads."""
    if not rows:
        raise ValueError("plan is empty")
    pos = rows[0].start
    for i, r in enumerate(rows):
        if r.start != pos:
            raise ValueError(f"row {i}: starts at {r.start}, previous ended at {pos}")
        step = md.U // r.mres
        if (r.target - r.start) % step:
            raise ValueError(f"row {i}: target off the 1/{r.mres} grid")
        if r.kind == "set_mres" and r.start % md.U:
            raise ValueError(f"row {i}: MRES change off a full step")
        if r.kind == "stair" and abs(r.target - r.start) != r.n_inc * r.inc_usteps * step:
            raise ValueError(f"row {i}: staircase distance does not match its increments")
        if r.speed_fs * r.mres > md.MAX_PULSE_RATE_HZ:
            raise ValueError(f"row {i}: step rate above {md.MAX_PULSE_RATE_HZ} Hz")
        pos = r.target
    for name, seg in segments(rows):
        if len(seg) > MAX_MOVES:
            raise ValueError(f"{name}: {len(seg)} rows, firmware holds {MAX_MOVES}")
        if seg[0].start != md.X_START or seg[-1].target != md.X_START:
            raise ValueError(f"{name}: does not start and end at x_start")


# ---------------------------------------------------------------- serial

def choose_port() -> str:
    try:
        from serial.tools import list_ports
    except ImportError as exc:
        raise RuntimeError("pyserial is required: python3 -m pip install pyserial") from exc
    ports = list(list_ports.comports())
    if not ports:
        raise RuntimeError("no serial ports found")
    print("\nSerial ports:")
    for i, p in enumerate(ports, 1):
        print(f"  {i}  {p.device}  ({p.description})")
    choice = input("Port number: ").strip()
    if not choice.isdigit() or not 1 <= int(choice) <= len(ports):
        raise RuntimeError("no port selected")
    return ports[int(choice) - 1].device


def open_board(port: str):
    import serial
    board = serial.Serial(port, BAUD, timeout=0.2, write_timeout=2)
    time.sleep(1.5)
    board.reset_input_buffer()
    return board


def read_line(board, deadline: float) -> str:
    while time.monotonic() < deadline:
        raw = board.readline()
        if raw:
            return raw.decode("utf-8", errors="replace").rstrip("\r\n")
    raise TimeoutError("board response timed out")


def send_expect(board, command: str, prefix: str, timeout: float = 5.0) -> list[str]:
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
    raise TimeoutError(f"no {prefix!r} response after {command!r}")


class SimulatedBoard:
    """Instant model of the protocol in ../firmware/README.md; touches nothing."""

    def __init__(self):
        self.responses = deque()
        self.moves: list[str] = []
        self.pos = 0

    def write(self, raw: bytes) -> int:
        for command in raw.decode().splitlines():
            self._command(command.strip())
        return len(raw)

    def _command(self, c: str) -> None:
        d = md.DRIVER
        if c == "PING":
            self.responses.append(f"# PONG,{FIRMWARE_ID},ready=1,driver=1,moves={len(self.moves)}")
        elif c == "DRIVER":
            self.responses.append(
                f"# DRIVER,spreadcycle=1,intpol={d['intpol']},stealth=0,"
                f"ihold_mult={d['ihold_multiplier']:.2f},tpowerdown={d['tpowerdown']},"
                f"iholddelay={d['iholddelay']},mres=16,mscnt=0")
        elif c.startswith("HOLD "):
            self.responses.append("# HOLD_OK," + c[5:].replace(" ", ","))
        elif c.startswith("INTPOL "):
            self.responses.append(f"# INTPOL_OK,intpol={c.split()[1]}")
        elif c == "ZERO":
            self.pos = 0
            self.responses.append("# ZERO_OK,pos=0,enc=0")
        elif c == "STATUS":
            self.responses.append(f"# STATUS,pos={self.pos},moves={len(self.moves)},running=0")
        elif c == "PLAN CLEAR":
            self.moves.clear()
            self.responses.append("# PLAN_CLEARED")
        elif c.startswith("PLAN ADD "):
            self.moves.append(c)
            self.responses.append(f"# PLAN_ADD,index={len(self.moves) - 1},"
                                  f"count={len(self.moves)},tag={c.split()[-1]}")
        elif c == "PLAN RUN":
            self.responses.append(f"# RUN_START,moves={len(self.moves)},pos={self.pos}")
            for i, m in enumerate(self.moves):
                f = m.split()
                kind, target, mres, tag = f[2], int(f[3]), f[6], f[-1]
                self.responses.append(f"# MOVE,{i},{kind},{tag},target={target},mres={mres}")
                if kind == "set_mres":
                    self.responses.append(f"# MRES,{mres},mscnt=0")
                self.responses.append(f"# S,{i * 20},{target},0,0,{i}")
                self.responses.append(f"# ARRIVED,{i},pos={target},enc=0,mscnt=0")
                self.pos = target
            self.responses.append(f"# RUN_COMPLETE,pos={self.pos}")
        elif c == "ABORT":
            self.responses.append("# ABORT")
        else:
            self.responses.append(f"# ERR,reason=unknown_command,command={c}")

    def readline(self) -> bytes:
        return (self.responses.popleft() + "\n").encode() if self.responses else b""

    def flush(self) -> None:
        pass

    def close(self) -> None:
        pass


def prepare_board(board) -> list[str]:
    """Firmware identity, driver mode, and the section 2 settings, read back."""
    lines = send_expect(board, "PING", "# PONG", timeout=10)
    pong = lines[-1].lower()
    if FIRMWARE_ID not in pong:
        raise RuntimeError(f"board is not running the MRES-test firmware: {lines[-1]}\n"
                           "  (the protocol is specified in firmware/README.md; "
                           "no sketch implements it yet - use Simulate)")
    if "ready=1" not in pong:
        raise RuntimeError(f"board is not ready: {lines[-1]}")
    d = md.DRIVER
    lines += send_expect(board, f"HOLD {d['ihold_multiplier']} {d['tpowerdown']} "
                                f"{d['iholddelay']}", "# HOLD_OK")
    lines += send_expect(board, f"INTPOL {d['intpol']}", "# INTPOL_OK")
    lines += send_expect(board, "DRIVER", "# DRIVER")
    drv = lines[-1].lower()
    if "spreadcycle=1" not in drv or "stealth=0" not in drv or f"intpol={d['intpol']}" not in drv:
        raise RuntimeError(f"driver readback does not match section 2: {lines[-1]}")
    return lines


# ------------------------------------------------------------------- run

def run_segment(board, name: str, rows, path: Path, header: dict) -> None:
    """Upload one segment, run it, and log every line the board sends."""
    send_expect(board, "PLAN CLEAR", "# PLAN_CLEARED")
    for i, r in enumerate(rows):
        reply = send_expect(board, plan_command(i, r), "# PLAN_ADD")[-1]
        if f"index={i},count={i + 1}" not in reply:
            raise RuntimeError(f"{name}: upload acknowledgement mismatch at row {i}: {reply}")
    status = send_expect(board, "STATUS", "# STATUS")[-1]
    if f"moves={len(rows)}" not in status:
        raise RuntimeError(f"{name}: board did not keep the whole plan: {status}")

    arrived = 0
    with path.open("w", encoding="utf-8", buffering=1) as log:
        log.write("# " + json.dumps({**header, "segment": name, "rows": len(rows),
                                     "host_start_utc": dt.datetime.now(dt.timezone.utc).isoformat()},
                                    sort_keys=True) + "\n")
        for line in send_expect(board, "PLAN RUN", "# RUN_START"):
            log.write(f"{time.time():.3f} {line}\n")
        while True:
            line = read_line(board, time.monotonic() + LINE_TIMEOUT_S)
            log.write(f"{time.time():.3f} {line}\n")
            if line.startswith("# ARRIVED"):
                arrived += 1
            if line.startswith("# ERR") or line.startswith("# RUN_FAILED"):
                raise RuntimeError(f"{name}: {line}")
            if line.startswith("# RUN_ABORTED"):
                raise RuntimeError(f"{name}: aborted")
            if line.startswith("# RUN_COMPLETE"):
                if arrived != len(rows):
                    raise RuntimeError(f"{name}: {arrived} arrivals, expected {len(rows)}")
                if f"pos={md.X_START}" not in line:
                    raise RuntimeError(f"{name}: did not finish at x_start: {line}")
                return


def confirm(rows, simulate: bool, assume_yes: bool) -> bool:
    if simulate or assume_yes:
        return True
    reach = max(abs(r.target) for r in rows) / md.U / md.FULL_STEPS_PER_MM
    print("\nNothing has moved. Before starting, confirm:")
    print("  - the stage is at x_start (centre of travel). Part 1 is not defined, so")
    print("    there is no homing: position it by hand or with the jog tool;")
    print("  - IDS 3010 acquisition is running (>= 100 Hz), environment compensation on;")
    print(f"  - travel is clear for +/-{reach:.1f} mm; motor power on; power removal in reach.")
    print(f"It will ZERO here, then run {len(rows)} rows, "
          f"{rows[-1].t_end / 60 - rows[0].t_start / 60:.1f} min.")
    return input("Type RUN to start: ").strip().upper() == "RUN"


def execute(rows, label: str, simulate: bool, assume_yes: bool, seed: int) -> int:
    validate(rows)
    board = SimulatedBoard() if simulate else open_board(choose_port())
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = DATA_DIR / (f"{stamp}_{label}" + ("_SIMULATED" if simulate else ""))
    try:
        board_lines = prepare_board(board)
        for line in board_lines:
            print(line)
        if not confirm(rows, simulate, assume_yes):
            print("Not started.")
            return 1
        out.mkdir(parents=True, exist_ok=True)
        with (out / "plan.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()))
            w.writeheader()
            w.writerows(asdict(r) for r in rows)
        header = {"mode": label, "seed": seed, "simulated": simulate,
                  "driver": md.DRIVER, "board": board_lines,
                  "parameters": md.summarise(seed, *md.build(seed))["parameters"]}
        (out / "session.json").write_text(json.dumps(header, indent=2) + "\n")
        send_expect(board, "ZERO", "# ZERO_OK")
        t0 = time.monotonic()
        segs = segments(rows)
        for k, (name, seg) in enumerate(segs, 1):
            mres = f"1/{seg[0].mres}"
            print(f"[{k}/{len(segs)}] {name:<9} {mres:>6}  "
                  f"{(seg[-1].t_end - seg[0].t_start) / 60:4.1f} min planned ...",
                  end="", flush=True)
            run_segment(board, name, seg, out / segment_file(name, seg), header)
            print(f" done ({time.monotonic() - t0:.0f} s elapsed)")
        print(f"Complete. Logs in {out}")
        return 0
    except KeyboardInterrupt:
        board.write(b"ABORT\n")
        board.flush()
        print("\nABORT sent. The session is stopped; logs so far are in", out)
        return 1
    finally:
        board.close()


# ----------------------------------------------------------------- modes

def mode_warmup(args) -> int:
    rows = md.build_warmup()
    if not rows:
        print("\nPart 1 (homing + warm-up) is not defined yet - it waits for the two\n"
              "end switches. Edit warmup() in scripts/mres_design.py when they are fitted.")
        return 1
    return execute(rows, "warmup", args.simulate, args.yes, args.seed)


def mode_stepping(args) -> int:
    order, rows = md.build_stepping(args.seed)
    print(f"\nStepping sequence: {' -> '.join(f'1/{m}' for m in order)}, "
          f"{rows[-1].t_end / 60:.1f} min")
    return execute(rows, "stepping", args.simulate, args.yes, args.seed)


def mode_full(args) -> int:
    if not md.build_warmup():
        print("\nA full session needs part 1, which is not defined yet.\n"
              "Run the stepping sequence (2) with the stage centred by hand instead.")
        return 1
    _, stepping = md.build_stepping(args.seed)
    warm = md.build_warmup()
    shift = warm[-1].t_end
    for r in stepping:
        r.t_start += shift
        r.t_end += shift
    return execute(warm + stepping, "full", args.simulate, args.yes, args.seed)


def mode_block(args) -> int:
    mres = args.mres
    if mres is None:
        choices = "/".join(str(m) for m in md.MRES_TESTED)
        text = input(f"MRES ({choices}): ").strip().removeprefix("1/")
        mres = int(text) if text.isdigit() else 0
    if mres not in md.MRES_TESTED:
        print(f"MRES must be one of {md.MRES_TESTED}")
        return 1
    rows = md.build_single_block(mres)
    print(f"\nSingle block at 1/{mres}: {rows[-1].t_end / 60:.1f} min")
    return execute(rows, f"block_MRES{mres}", args.simulate, args.yes, args.seed)


def mode_preview(args) -> int:
    order, rows = md.build_stepping(args.seed)
    validate(rows)
    print(f"\nPart 1  warm-up: not defined yet")
    print(f"Part 2  stepping: {' -> '.join(f'1/{m}' for m in order)}, "
          f"{len(rows)} rows, {rows[-1].t_end / 60:.1f} min\n")
    print(f"  {'segment':<9} {'MRES':>6} {'rows':>5} {'min':>6}  file")
    for name, seg in segments(rows):
        print(f"  {name:<9} {'1/' + str(seg[0].mres):>6} {len(seg):5d} "
              f"{(seg[-1].t_end - seg[0].t_start) / 60:6.1f}  {segment_file(name, seg)}")
    if args.commands or input("\nShow every command? [y/N] ").strip().lower() == "y":
        for name, seg in segments(rows):
            print(f"\n# --- {name}")
            for i, r in enumerate(seg):
                print(plan_command(i, r))
    return 0


def mode_simulate(args) -> int:
    args.simulate = True
    return mode_stepping(args)


def mode_regenerate(args) -> int:
    for script in ("mres_design.py", "render_mres.py"):
        result = subprocess.run([sys.executable, str(HERE / script), "--seed", str(args.seed)])
        if result.returncode:
            print(f"{script} failed ({result.returncode})")
            return result.returncode
    return 0


HANDLERS = {"warmup": mode_warmup, "stepping": mode_stepping, "full": mode_full,
            "block": mode_block, "preview": mode_preview, "simulate": mode_simulate,
            "regenerate": mode_regenerate}


def menu(args) -> int:
    order, rows = md.build_stepping(args.seed)
    warm = "not defined yet" if not md.build_warmup() else "defined"
    while True:
        print("\n" + "=" * 64)
        print(" MRES test - microstep-resolution selection")
        print(f" part 1 warm-up: {warm}   |   part 2: {len(order)} blocks, "
              f"{rows[-1].t_end / 60:.1f} min")
        print("=" * 64)
        for key, _, title, desc in MODES:
            print(f"  {key}  {title:<18} {desc}")
        print("  q  Quit")
        choice = input("\nSelect a mode: ").strip().lower()
        if choice in ("q", "quit", "exit"):
            return 0
        mode = next((m for k, m, *_ in MODES if choice in (k, m)), None)
        if mode is None:
            print("Unknown choice.")
            continue
        try:
            HANDLERS[mode](args)
        except (OSError, RuntimeError, TimeoutError, ValueError) as exc:
            print(f"error: {exc}")
        args.simulate = args.simulate_flag
        input("\nPress Enter to return to the menu ...")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=list(HANDLERS), help="skip the menu")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--mres", type=int, help="MRES for --mode block")
    ap.add_argument("--simulate", action="store_true", help="use the protocol model")
    ap.add_argument("--yes", action="store_true", help="skip the typed RUN confirmation")
    ap.add_argument("--commands", action="store_true", help="preview: print every command")
    args = ap.parse_args()
    args.simulate_flag = args.simulate
    if args.mode:
        return HANDLERS[args.mode](args)
    return menu(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TimeoutError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
    except (KeyboardInterrupt, EOFError):
        print()
        raise SystemExit(130)
