#!/usr/bin/env python3
"""Raw-key terminal and lossless logger for the E3 ESP32 firmware."""

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
DATA_DIR = Path(__file__).resolve().parent / "data"
QUIT = b"\x1d"       # Ctrl+]
ABORT = b"\x03"      # Ctrl+C


def find_port() -> str:
    ports = sorted(glob.glob("/dev/cu.usbmodem*"))
    if not ports:
        raise SystemExit("no /dev/cu.usbmodem* port found")
    return ports[0]


def open_port(port: str) -> serial.Serial:
    connection = serial.Serial()
    connection.port = port
    connection.baudrate = BAUD
    connection.timeout = 0.05
    connection.dtr = connection.rts = False
    connection.open()
    return connection


def reader(connection: serial.Serial, stop: threading.Event, log) -> None:
    pending = ""
    while not stop.is_set():
        try:
            payload = connection.read(4096)
        except serial.SerialException as exc:
            sys.stdout.write(f"\r\n[serial error: {exc}]\r\n")
            stop.set()
            return
        if not payload:
            continue
        pending += payload.decode("utf-8", errors="replace").replace("\r", "")
        *lines, pending = pending.split("\n")
        for line in lines:
            log.write(f"{time.time():.3f} {line}\n")
            sys.stdout.write(line + "\n")
        if pending and not pending.startswith("#"):
            sys.stdout.write(pending)
            log.write(f"{time.time():.3f} {pending}\n")
            pending = ""
        sys.stdout.flush()


def main() -> int:
    port = sys.argv[1] if len(sys.argv) > 1 else find_port()
    connection = open_port(port)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    log_path = DATA_DIR / dt.datetime.now().strftime("e3_terminal_%Y%m%d_%H%M%S.log")
    log = log_path.open("w", encoding="utf-8", buffering=1)
    print(f"E3 terminal — {port}")
    print(f"Logging to {log_path}")
    print("Keys act immediately. Ctrl+C aborts; Ctrl+] closes the terminal.\n")

    stop = threading.Event()
    thread = threading.Thread(target=reader, args=(connection, stop, log), daemon=True)
    thread.start()
    fd = sys.stdin.fileno()
    interactive = os.isatty(fd)
    saved = termios.tcgetattr(fd) if interactive else None
    try:
        if interactive:
            tty.setraw(fd, termios.TCSANOW)
            attrs = termios.tcgetattr(fd)
            attrs[1] |= termios.OPOST | termios.ONLCR
            termios.tcsetattr(fd, termios.TCSANOW, attrs)
        while not stop.is_set():
            key = os.read(fd, 1)
            if not key or key == QUIT:
                break
            connection.write(ABORT if key == ABORT else key)
            if key == ABORT:
                sys.stdout.write("\n[ABORT sent]\n")
                sys.stdout.flush()
    finally:
        if saved is not None:
            termios.tcsetattr(fd, termios.TCSADRAIN, saved)
        stop.set()
        time.sleep(0.1)
        connection.close()
        log.close()
        print("\n[closed]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
