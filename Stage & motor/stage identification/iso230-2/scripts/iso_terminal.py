#!/usr/bin/env python3
"""Control terminal for the ISO 230-2 firmware.

    /usr/bin/python3 scripts/iso_terminal.py [port]

Every key goes to the board the moment it is pressed - no Enter needed. The
firmware's menus act on single keys; a key stops homing, latches the end of
the warm-up, and shows progress during the ISO sequence. Prompts for numbers
echo what you type and take Enter.

  Ctrl+C       ABORT - sent as the single byte 0x03, which the firmware always
               treats as ABORT (a typed "ABORT" would be read as keys while a
               state is running)
  Ctrl+]       close the terminal

The 50 Hz `# S,...` position samples are not printed - they would scroll over
everything else - but everything the board sends, samples included, is logged
to data/terminal_<stamp>.log.

The port is opened with DTR and RTS low. On the ESP32-S3's built-in USB serial
those lines drive reset and the boot strap, and the usual default (both high)
restarts the chip into the ROM bootloader, where it waits silently.
"""

from __future__ import annotations

import datetime as dt
import glob
import os
import sys
import termios
import threading
import time
import tty
from pathlib import Path

import serial

BAUD = 115200
LOG_DIR = Path(__file__).resolve().parent.parent / "data"
QUIT = b"\x1d"                     # Ctrl+]
INTERRUPT = b"\x03"                # Ctrl+C


def find_port() -> str:
    ports = sorted(glob.glob("/dev/cu.usbmodem*"))
    if not ports:
        raise SystemExit("no /dev/cu.usbmodem* port - is the board plugged in?")
    return ports[0]


def open_port(port: str) -> serial.Serial:
    s = serial.Serial()
    s.port, s.baudrate, s.timeout = port, BAUD, 0.05
    s.dtr = s.rts = False
    s.open()
    return s


def reader(s: serial.Serial, stop: threading.Event, log) -> None:
    pending = ""
    while not stop.is_set():
        try:
            data = s.read(4096)
        except serial.SerialException as exc:
            sys.stdout.write(f"\r\n[serial error: {exc}]\r\n")
            stop.set()
            return
        if not data:
            continue
        pending += data.decode("utf-8", errors="replace").replace("\r", "")
        *lines, pending = pending.split("\n")
        for line in lines:
            log.write(f"{time.time():.3f} {line}\n")
            if not line.startswith("# S,"):
                sys.stdout.write(line + "\n")
        # A prompt has no newline after it; show it rather than hold it back.
        # Protocol lines all start with '#' and always end in one, so wait.
        if pending and not pending.startswith("#"):
            sys.stdout.write(pending)
            log.write(f"{time.time():.3f} {pending}\n")
            pending = ""
        sys.stdout.flush()


def main() -> int:
    port = sys.argv[1] if len(sys.argv) > 1 else find_port()
    s = open_port(port)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / dt.datetime.now().strftime("terminal_%Y%m%d_%H%M%S.log")
    log = log_path.open("w", encoding="utf-8", buffering=1)
    print(f"MRES test control terminal - {port}")
    print(f"Logging everything (position samples included) to {log_path.name}")
    print("Keys act at once. ? = menu, Ctrl+C = ABORT, Ctrl+] = close.\n")

    stop = threading.Event()
    threading.Thread(target=reader, args=(s, stop, log), daemon=True).start()
    s.write(b"?")

    fd = sys.stdin.fileno()
    interactive = os.isatty(fd)
    saved = termios.tcgetattr(fd) if interactive else None
    try:
        if interactive:
            # Raw keys in, but keep output processing so "\n" still starts a
            # new line. Ctrl+C arrives as a byte, not a signal, so it can
            # become ABORT.
            tty.setraw(fd, termios.TCSANOW)
            attrs = termios.tcgetattr(fd)
            attrs[1] |= termios.OPOST | termios.ONLCR
            termios.tcsetattr(fd, termios.TCSANOW, attrs)
        while not stop.is_set():
            key = os.read(fd, 1)
            if not key or key == QUIT:
                break
            if key == INTERRUPT:
                s.write(INTERRUPT)
                sys.stdout.write("\n[ABORT sent]\n")
                sys.stdout.flush()
                continue
            if key == b"\r":
                key = b"\n"
            s.write(key)
            if not interactive:
                time.sleep(0.05)
    finally:
        if saved is not None:
            termios.tcsetattr(fd, termios.TCSADRAIN, saved)
        stop.set()
        time.sleep(0.1)
        s.close()
        log.close()
        print("\n[closed]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
