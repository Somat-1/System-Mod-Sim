#!/usr/bin/env python3
"""Drive the ESP32-S3 / TMC2209 stage: pick two positions, shuttle between them.

Standalone: one file, one dependency (pyserial). Run it and it finds the ESP32,
asks for two positions as whole numbers, and then moves the stage back and
forth between them until you stop it.

    python jog_stepper.py                 # find the port, then ask
    python jog_stepper.py --port COM11    # skip the search
    python jog_stepper.py --simulate      # no hardware; try the interface out
    python jog_stepper.py --a -200 --b 200 --go    # straight into motion

It needs esp32s3_tmc2209_jog.ino running on the board - a host script cannot
pulse STEP and DIR by itself, those are ESP32 pins. `--simulate` stands in a
model of that firmware so the interface can be exercised with nothing attached.

Positions are whole **full steps**; 100 full steps = 1 mm on this rig, so one
unit is 10 um. Every prompt shows the millimetre equivalent. Speed is a
full-step rate in steps/s, matching the campaign sketches: 70 steps/s is the
"MODERATE" rate they use, i.e. 0.7 mm/s.

Safety, because this moves a real stage with no homing and no end stops:

* Position 0 means "wherever the stage was when it was last zeroed", nothing
  more. There is no absolute reference.
* A soft limit of +/-500 full steps (+/-5 mm) about zero is applied by the
  firmware and refuses targets outside it. The stage's measured travel is
  44 mm, so widen it with --limit once you know where you are within that.
* Nothing moves until you ask. Ctrl-C, or any key while jogging, stops it.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

FULL_STEPS_PER_MM = 100
DEFAULT_RATE_STEPS_S = 70          # 0.7 mm/s, the rigs "MODERATE"
DEFAULT_LIMIT = 500                # +/-5 mm
BAUD = 115200

# ESP32-S3 native USB CDC, plus the usual bridge chips, so an older board on
# the same rig is still found.
USB_VIDS = {
    0x303A: "Espressif (ESP32-S3 native USB)",
    0x10C4: "Silicon Labs CP210x",
    0x1A86: "WCH CH34x",
    0x0403: "FTDI",
}


def mm(steps: float) -> str:
    return f"{steps / FULL_STEPS_PER_MM:+.2f} mm"


# --------------------------------------------------------------------------
# transport
# --------------------------------------------------------------------------

@dataclass
class Link:
    """A line-oriented connection to the board, read on a background thread."""

    port: str
    ser: object
    lines: list = field(default_factory=list)
    lock: threading.Lock = field(default_factory=threading.Lock)
    alive: bool = True
    verbose: bool = False

    def start(self):
        self._t = threading.Thread(target=self._reader, daemon=True)
        self._t.start()

    def _reader(self):
        buf = b""
        while self.alive:
            try:
                chunk = self.ser.read(256)
            except Exception:
                break
            if not chunk:
                continue
            buf += chunk
            while b"\n" in buf:
                raw, buf = buf.split(b"\n", 1)
                text = raw.decode("utf-8", "replace").strip()
                if not text:
                    continue
                with self.lock:
                    self.lines.append(text)
                if self.verbose:
                    print(f"  <- {text}")

    def send(self, cmd: str):
        if self.verbose:
            print(f"  -> {cmd}")
        self.ser.write((cmd + "\n").encode())
        try:
            self.ser.flush()
        except Exception:
            pass

    def drain(self) -> list:
        with self.lock:
            out, self.lines = self.lines, []
        return out

    def wait_for(self, prefixes, timeout=3.0) -> str | None:
        """Collect lines until one starts with any of `prefixes`."""
        if isinstance(prefixes, str):
            prefixes = (prefixes,)
        deadline = time.monotonic() + timeout
        seen = []
        while time.monotonic() < deadline:
            for line in self.drain():
                seen.append(line)
                for p in prefixes:
                    if line.startswith(p):
                        self._replay(seen, line)
                        return line
            time.sleep(0.01)
        self._replay(seen, None)
        return None

    def _replay(self, seen, matched):
        """Put back everything that was not what we were waiting for."""
        keep = [s for s in seen if s is not matched]
        if keep:
            with self.lock:
                self.lines = keep + self.lines

    def close(self):
        self.alive = False
        try:
            self.ser.close()
        except Exception:
            pass


def parse_kv(line: str) -> dict:
    """'# STATUS,pos=12,mode=IDLE' -> {'_tag': 'STATUS', 'pos': '12', ...}"""
    body = line.lstrip("#").strip()
    parts = [p for p in body.split(",") if p]
    out = {}
    if parts and "=" not in parts[0]:
        out["_tag"] = parts[0]
        parts = parts[1:]
    for p in parts:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    return out


# --------------------------------------------------------------------------
# port discovery
# --------------------------------------------------------------------------

def list_ports():
    from serial.tools import list_ports as lp
    return list(lp.comports())


def describe(p) -> str:
    vid = getattr(p, "vid", None)
    tag = USB_VIDS.get(vid, "")
    bits = [p.device]
    if p.description and p.description != "n/a":
        bits.append(p.description)
    if vid is not None:
        bits.append(f"VID:PID {vid:04X}:{getattr(p, 'pid', 0) or 0:04X}")
    if tag:
        bits.append(tag)
    return "  ".join(bits)


def probe(port: str, verbose=False, settle=1.2) -> Link | None:
    """Open a port and see whether the jog firmware answers."""
    import serial
    try:
        ser = serial.Serial(port, BAUD, timeout=0.1)
    except Exception as exc:
        print(f"  {port}: cannot open ({exc})")
        return None
    link = Link(port=port, ser=ser, verbose=verbose)
    link.start()
    # The board may reset when the port opens; give its banner time to arrive.
    time.sleep(settle)
    for _ in range(3):
        link.send("PING")
        reply = link.wait_for("# PONG", timeout=1.2)
        if reply:
            return link
    link.close()
    return None


def find_board(preferred: str | None, verbose=False) -> Link | None:
    """Try the likeliest ports first, confirming each with a PING."""
    ports = list_ports()
    if not ports:
        print("No serial ports found at all. Is the board plugged in?")
        return None

    print("Serial ports:")
    for p in ports:
        print("  " + describe(p))

    names = [p.device for p in ports]
    order: list[str] = []

    def add(name):
        if name in names and name not in order:
            order.append(name)

    if preferred:
        add(preferred)
        if preferred not in names:
            print(f"\n{preferred} was asked for but is not present; trying the rest.")
    add("COM11")                                     # the rig's usual port
    for p in ports:                                  # then known USB-serial VIDs
        if getattr(p, "vid", None) in USB_VIDS:
            add(p.device)
    for name in names:
        add(name)

    print()
    for name in order:
        print(f"Probing {name} ...")
        link = probe(name, verbose=verbose)
        if link:
            print(f"Found the jog firmware on {name}.")
            return link
        print(f"  {name}: no response to PING")
    print("\nNothing answered. Either the board is on another port, or it is not")
    print("running esp32s3_tmc2209_jog.ino (see firmware/ in this folder).")
    return None


# --------------------------------------------------------------------------
# flashing
# --------------------------------------------------------------------------
# Two ways to get the firmware onto the board, because they suit very different
# machines.
#
#   prebuilt  a committed flash image written with esptool. esptool is pure
#             Python and pip-installs in seconds, so this needs no compiler
#             and no toolchain - the right choice for a locked-down rig PC.
#   source    arduino-cli compiles the sketch and uploads it. Always current
#             with the .ino, but the first run pulls the ESP32 core and its
#             toolchains - measured at over 5 GB, so not something to start on
#             a metered or slow connection.
#
# Prebuilt is tried first and source is the fallback, so the common case is
# fast and the uncommon case still works.

HERE = Path(__file__).resolve().parent
SKETCH_DIR = HERE / "firmware" / "esp32s3_tmc2209_jog"
PREBUILT_BIN = HERE / "firmware" / "prebuilt" / "esp32s3_tmc2209_jog_merged.bin"
TOOLS_DIR = HERE / "tools"

# Both settings matter and both are in the FQBN: without CDCOnBoot=cdc the
# sketch runs but Serial never appears on USB, which looks exactly like a board
# that will not answer PING.
FQBN = "esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc"
CORE = "esp32:esp32"
LIBS = ("TMCStepper",)
FLASH_OFFSET = "0x0"
CHIP = "esp32s3"

ACLI_URL = {
    "win32": "https://downloads.arduino.cc/arduino-cli/arduino-cli_latest_Windows_64bit.zip",
    "linux": "https://downloads.arduino.cc/arduino-cli/arduino-cli_latest_Linux_64bit.tar.gz",
    "darwin": "https://downloads.arduino.cc/arduino-cli/arduino-cli_latest_macOS_64bit.tar.gz",
}


def run(cmd, **kw):
    print("    $ " + " ".join(str(c) for c in cmd))
    kw.setdefault("check", False)
    return subprocess.run([str(c) for c in cmd], **kw)


def prebuilt_is_current() -> bool | None:
    """Does the committed image still match the sketch on disk?

    A prebuilt binary that silently lags an edited .ino is the obvious trap
    here: you fix the firmware, flash, and get the old build with no hint that
    anything was skipped. The manifest records the sketch's hash at build time,
    so the mismatch can be reported instead. None = cannot tell.
    """
    import hashlib
    import json
    manifest = PREBUILT_BIN.parent / "manifest.json"
    ino = SKETCH_DIR / "esp32s3_tmc2209_jog.ino"
    if not manifest.exists() or not ino.exists():
        return None
    try:
        recorded = json.loads(manifest.read_text(encoding="utf-8")).get("sketch_sha256")
    except (OSError, ValueError):
        return None
    if not recorded:
        return None
    return hashlib.sha256(ino.read_bytes()).hexdigest() == recorded


def flash_prebuilt(port: str) -> bool:
    """Write the committed flash image with esptool. No compiler involved."""
    import importlib
    if not PREBUILT_BIN.exists():
        print(f"  no prebuilt image at {PREBUILT_BIN.relative_to(HERE)}")
        return False
    current = prebuilt_is_current()
    if current is False:
        print("  NOTE: the sketch has been edited since this image was built.")
        print("  The prebuilt image is what will be flashed, which will NOT")
        print("  include your changes. Use --flash-mode source to compile them.")
    elif current is None:
        print("  (cannot confirm the image matches the sketch - no manifest)")
    try:
        importlib.import_module("esptool")
    except ImportError:
        print("  esptool is not installed. It is pure Python and small:")
        print(f"    {sys.executable} -m pip install esptool")
        return False
    size_mb = PREBUILT_BIN.stat().st_size / 1e6
    print(f"  flashing {PREBUILT_BIN.name} ({size_mb:.2f} MB) to {port} ...")
    # esptool 5 renamed the commands and option values to hyphens and warns
    # loudly about the old spellings; esptool 4 only knows the old ones. Try
    # the modern spelling first so current installs are quiet, and fall back
    # rather than requiring a particular version.
    variants = (
        ("write-flash", "default-reset", "hard-reset"),
        ("write_flash", "default_reset", "hard_reset"),
    )
    for verb, before, after in variants:
        cmd = [sys.executable, "-m", "esptool", "--chip", CHIP, "--port", port,
               "--baud", "921600", "--before", before, "--after", after,
               verb, FLASH_OFFSET, str(PREBUILT_BIN)]
        res = run(cmd, capture_output=True, text=True)
        out = (res.stdout or "") + (res.stderr or "")
        if res.returncode == 0:
            print("  flashed.")
            return True
        # Only a spelling the installed esptool does not know is worth retrying;
        # anything else is a real failure and retrying would just hide it.
        if not any(s in out for s in
                   ("invalid choice", "unrecognized arguments", "No such command")):
            print(out.strip()[-1500:])
            return False
    print("  esptool rejected both command spellings; is it a supported version?")
    return False


def find_arduino_cli(allow_download=False):
    """arduino-cli from PATH, from tools/, or downloaded on request."""
    exe = "arduino-cli.exe" if sys.platform == "win32" else "arduino-cli"
    found = shutil.which("arduino-cli") or (
        str(TOOLS_DIR / exe) if (TOOLS_DIR / exe).exists() else None)
    if found:
        return found
    if not allow_download:
        print("  arduino-cli was not found on PATH or in tools/.")
        print("  re-run with --allow-download to fetch it (~18 MB) into tools/,")
        print("  or install it yourself: https://arduino.github.io/arduino-cli/")
        return None
    url = ACLI_URL.get(sys.platform)
    if not url:
        print(f"  no arduino-cli download known for platform {sys.platform}")
        return None
    import io
    import tarfile
    import urllib.request
    import zipfile
    print(f"  downloading arduino-cli from {url}")
    TOOLS_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=180) as r:
        blob = r.read()
    if url.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            z.extractall(TOOLS_DIR)
    else:
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as t:
            t.extractall(TOOLS_DIR)
    target = TOOLS_DIR / exe
    if not target.exists():
        print("  download did not produce arduino-cli")
        return None
    if sys.platform != "win32":
        target.chmod(0o755)
    print(f"  arduino-cli in {target}")
    return str(target)


def flash_from_source(port: str, allow_download=False) -> bool:
    """Compile the sketch with arduino-cli and upload it."""
    cli = find_arduino_cli(allow_download)
    if cli is None:
        return False
    if not SKETCH_DIR.exists():
        print(f"  sketch folder missing: {SKETCH_DIR}")
        return False
    print("  checking the ESP32 core and libraries "
          "(the first run downloads >5 GB) ...")
    env = dict(os.environ)
    if run([cli, "core", "list"], capture_output=True, text=True,
           env=env).stdout.find("esp32:esp32") < 0:
        run([cli, "core", "update-index"], env=env)
        if run([cli, "core", "install", CORE], env=env).returncode != 0:
            print("  could not install the ESP32 core")
            return False
    for lib in LIBS:
        run([cli, "lib", "install", lib], env=env, capture_output=True, text=True)
    print("  compiling ...")
    if run([cli, "compile", "--fqbn", FQBN, SKETCH_DIR], env=env).returncode != 0:
        print("  compile failed")
        return False
    print(f"  uploading to {port} ...")
    if run([cli, "upload", "--fqbn", FQBN, "--port", port, SKETCH_DIR],
           env=env).returncode != 0:
        print("  upload failed")
        return False
    print("  flashed.")
    return True


def do_flash(port: str, mode="auto", allow_download=False) -> bool:
    print(f"\nFlashing the jog firmware to {port}")
    if mode in ("auto", "prebuilt"):
        if flash_prebuilt(port):
            return True
        if mode == "prebuilt":
            return False
        print("  falling back to compiling from source")
    return flash_from_source(port, allow_download)


# --------------------------------------------------------------------------
# simulator - stands in for the firmware so the interface can be tried out
# --------------------------------------------------------------------------

class FakeFirmware:
    """A model of the sketch, accurate enough to exercise the host script.

    Deliberately not a physics model: it reproduces the protocol, the soft
    limit, the timing and the mode transitions, which is what the host code
    actually depends on.
    """

    def __init__(self):
        self.pos = 0.0           # full steps, fractional while moving
        self.a = 0
        self.b = 0
        self.target = 0
        self.mode = "IDLE"
        self.toward_b = True
        self.rate = DEFAULT_RATE_STEPS_S * 1000
        self.mres = 16
        self.current = 360
        self.limit = DEFAULT_LIMIT
        self.enc = 0
        self.calib_end = 0.0
        self.calib_dwell = 2.0
        self.calib_next = 0.0
        self.calib_cycles = 0
        self.calib_out = True
        self.out: list[str] = []
        self.last = time.monotonic()
        self.last_tel = 0.0
        self.lock = threading.Lock()
        self.emit("# BOOT,fw=esp32s3_tmc2209_jog (SIMULATED)")
        self.emit("# PINS,step=5,dir=6,uart_rx=18,uart_tx=17,en=grounded")
        self.emit(f"# READY,ready=1,mres={self.mres},current_ma={self.current},"
                  f"rate_millihz={self.rate},limit={self.limit},"
                  f"steps_per_mm={FULL_STEPS_PER_MM},enc={self.enc},index=0")
        self.emit("# IDLE,nothing moves until commanded")

    def emit(self, line):
        with self.lock:
            self.out.append(line)

    def read_lines(self):
        with self.lock:
            out, self.out = self.out, []
        return out

    def status(self):
        return (f"# STATUS,pos={round(self.pos)},target={self.target},"
                f"a={self.a},b={self.b},mode={self.mode},"
                f"rate_millihz={self.rate},mres={self.mres},"
                f"current_ma={self.current},limit={self.limit},"
                f"steps_per_mm={FULL_STEPS_PER_MM},enc={self.enc},index=0")

    def within(self, v):
        return self.limit <= 0 or -self.limit <= v <= self.limit

    @staticmethod
    def _num(text):
        try:
            return int(text)
        except ValueError:
            return None

    def command(self, line):
        line = line.strip()
        if not line:
            return
        verb, _, rest = line.partition(" ")
        verb, rest = verb.upper(), rest.strip()
        if verb == "PING":
            self.emit("# PONG,fw=esp32s3_tmc2209_jog,ready=1")
        elif verb == "STATUS":
            self.emit(self.status())
        elif verb == "POS?":
            self.emit(f"# POS,pos={round(self.pos)}")
        elif verb == "STOP":
            self.mode = "IDLE"
            self.emit(f"# STOP,pos={round(self.pos)},reason=command")
        elif verb == "ZERO":
            if self.mode != "IDLE":
                self.emit("# ERR,reason=busy")
            else:
                self.pos = 0.0
                self.emit("# ZERO_OK,pos=0")
        elif verb == "LIMIT":
            v = self._num(rest)
            if v is None or v < 0:
                self.emit("# ERR,reason=bad_limit")
            else:
                self.limit = v
                self.emit(f"# LIMIT_OK,limit={self.limit}")
        elif verb == "RATE":
            v = self._num(rest)
            # mirror the firmware's pulse-rate ceiling
            if v is None or v <= 0 or v * self.mres > 20000 * 1000:
                self.emit("# ERR,reason=bad_rate")
            else:
                self.rate = v
                self.emit(f"# RATE_OK,rate_millihz={self.rate},"
                          f"pulse_hz={self.rate*self.mres/1000:.1f}")
        elif verb == "CURRENT":
            v = self._num(rest)
            if v is None or not 50 <= v <= 1200:
                self.emit("# ERR,reason=bad_current")
            else:
                self.current = v
                self.emit(f"# CURRENT_OK,current_ma={self.current}")
        elif verb == "MRES":
            v = self._num(rest)
            if self.mode != "IDLE":
                self.emit("# ERR,reason=busy")
            elif v is None or not (1 <= v <= 256 and v & (v - 1) == 0):
                self.emit("# ERR,reason=bad_mres")
            else:
                self.mres = v
                self.emit(f"# MRES_OK,mres={self.mres},"
                          f"pulse_hz={self.rate*self.mres/1000:.1f}")
        elif verb == "MOVETO":
            v = int(rest)
            if not self.within(v):
                self.emit("# ERR,reason=outside_limit")
            else:
                self.target, self.mode = v, "MOVING"
                self.emit(f"# MOVETO,target={v}")
        elif verb == "JOG":
            try:
                a, b = (int(x) for x in rest.split())
            except ValueError:
                self.emit("# ERR,reason=jog_needs_two")
                return
            if not (self.within(a) and self.within(b)):
                self.emit("# ERR,reason=outside_limit")
            elif a == b:
                self.emit("# ERR,reason=equal_positions")
            else:
                self.a, self.b = a, b
                self.emit(f"# JOG_SET,a={a},b={b}")
        elif verb == "START":
            if self.a == self.b:
                self.emit("# ERR,reason=no_jog_set")
            else:
                self.toward_b = abs(self.b - self.pos) >= abs(self.a - self.pos)
                self.mode = "JOGGING"
                self.emit(f"# START,a={self.a},b={self.b},"
                          f"first={'B' if self.toward_b else 'A'}")
        elif verb == "ENC?":
            self.emit(f"# ENC,enc={self.enc},index=0")
        elif verb == "CALIB":
            parts = rest.split()
            secs = int(parts[0]) if parts else 30
            dwell = int(parts[1]) if len(parts) > 1 else 2000
            if not (0 < secs <= 3600 and 50 <= dwell <= 60000):
                self.emit("# ERR,reason=bad_calib_args")
                return
            self.mode = "CALIB"
            self.calib_end = time.monotonic() + secs
            self.calib_dwell = dwell / 1000.0
            self.calib_next = time.monotonic()
            self.calib_cycles = 0
            self.calib_out = True
            self.emit(f"# CALIB_START,seconds={secs},dwell_ms={dwell},"
                      f"mres={self.mres},microstep_um={10.0/self.mres:.4f}")
        else:
            self.emit(f"# ERR,reason=unknown_command,verb={verb}")

    def tick_calib(self, now):
        if now >= self.calib_end:
            self.emit(f"# CALIB_DONE,cycles={self.calib_cycles},"
                      f"pos={round(self.pos)},enc={self.enc}")
            self.mode = "IDLE"
            return
        if now < self.calib_next:
            return
        micro = 1.0 / self.mres                      # in full steps
        self.pos += micro if self.calib_out else -micro
        # A plausible stand-in for the real thing: the encoder follows the
        # command with a little lost motion on reversal. It is invented, and
        # exists only so the host's display can be exercised.
        self.enc += round(micro * 200 * (0.94 if not self.calib_out else 1.0))
        self.calib_next = now + self.calib_dwell
        if not self.calib_out:
            self.calib_cycles += 1
        self.emit(f"# CALIB,{'+' if self.calib_out else '-'},"
                  f"cycle={self.calib_cycles},pos_u256={round(self.pos*256)},"
                  f"enc={self.enc},t_ms={int(now*1000)}")
        self.calib_out = not self.calib_out

    def tick(self):
        now = time.monotonic()
        dt, self.last = now - self.last, now
        if self.mode == "CALIB":
            self.tick_calib(now)
            return
        if self.mode != "IDLE":
            tgt = (self.b if self.toward_b else self.a) if self.mode == "JOGGING" \
                else self.target
            span = self.rate / 1000.0 * dt
            if abs(tgt - self.pos) <= span:
                self.pos = float(tgt)
                if self.mode == "JOGGING":
                    self.emit(f"# REACHED,pos={round(self.pos)},"
                              f"next={'A' if self.toward_b else 'B'}")
                    self.toward_b = not self.toward_b
                else:
                    self.emit(f"# ARRIVED,pos={round(self.pos)}")
                    self.mode = "IDLE"
            else:
                self.pos += span if tgt > self.pos else -span
            if now - self.last_tel >= 0.1:
                self.last_tel = now
                self.emit(f"# POS,pos={round(self.pos)},mode={self.mode}")


class FakeSerial:
    """Just enough of a pyserial object for Link to drive the simulator."""

    def __init__(self):
        self.fw = FakeFirmware()
        self._buf = b""
        self._stop = False
        threading.Thread(target=self._run, daemon=True).start()

    def _run(self):
        while not self._stop:
            self.fw.tick()
            time.sleep(0.005)

    def read(self, _n=1):
        for line in self.fw.read_lines():
            self._buf += (line + "\n").encode()
        out, self._buf = self._buf, b""
        if not out:
            time.sleep(0.02)
        return out

    def write(self, data):
        for line in data.decode().splitlines():
            self.fw.command(line)

    def flush(self):
        pass

    def close(self):
        self._stop = True


# --------------------------------------------------------------------------
# keyboard: "any key stops it" without blocking the status display
# --------------------------------------------------------------------------

class KeyWatcher:
    """Non-blocking 'was a key pressed' check, so jogging can be interrupted.

    A plain input() would block the status line, and a reader thread printing
    over a half-typed command makes a mess of the terminal. Checking for a
    keypress instead keeps the display in one place and still lets a single key
    stop the stage.
    """

    def __init__(self):
        self.kind = None
        if sys.platform == "win32":
            try:
                import msvcrt  # noqa: F401
                self.kind = "msvcrt"
            except ImportError:
                pass
        elif sys.stdin.isatty():
            try:
                import termios, tty  # noqa: F401
                self.kind = "termios"
            except ImportError:
                pass

    def __enter__(self):
        if self.kind == "termios":
            import termios, tty
            self._fd = sys.stdin.fileno()
            self._saved = termios.tcgetattr(self._fd)
            tty.setcbreak(self._fd)
        return self

    def __exit__(self, *exc):
        if self.kind == "termios":
            import termios
            termios.tcsetattr(self._fd, termios.TCSADRAIN, self._saved)

    def pressed(self) -> bool:
        if self.kind == "msvcrt":
            import msvcrt
            hit = False
            while msvcrt.kbhit():
                msvcrt.getwch()
                hit = True
            return hit
        if self.kind == "termios":
            import select
            if select.select([sys.stdin], [], [], 0)[0]:
                sys.stdin.read(1)
                return True
        return False


# --------------------------------------------------------------------------
# the interface
# --------------------------------------------------------------------------

HELP = """
Commands
  <number> <number>   set both positions at once, e.g.  -200 200
  a <number>          set position A            b <number>   set position B
  go                  shuttle between A and B until a key is pressed
  go <n>              ... or for exactly n legs, then stop
  calib [s] [dwell]   microstep calibration: +1 / -1 microstep with long
                      dwells, default 30 s and 2000 ms. Set the step size
                      first with 'mres' - at mres 16 one microstep is
                      0.625 um.
  moveto <number>     one move to an absolute position
  speed <steps/s>     full-step rate      (now: {rate} steps/s, {mmps:.2f} mm/s)
  mres <1..256>       microstep resolution
  current <mA>        coil current
  limit <steps>       soft limit, +/- about zero; 0 disables
  zero                call the present position 0
  status              ask the board where it is
  help                this list            quit / q      leave

Positions are whole full steps; {spm} full steps = 1 mm, so 1 unit = 10 um.
"""


class Console:
    def __init__(self, link: Link, limit: int, rate: int):
        self.link = link
        self.state = {}
        self.a = None
        self.b = None
        self.limit = limit
        self.rate = rate
        self.refresh()

    # -- board interaction ------------------------------------------------
    def refresh(self) -> dict:
        self.link.drain()
        self.link.send("STATUS")
        line = self.link.wait_for("# STATUS", timeout=2.0)
        if line:
            self.state = parse_kv(line)
        return self.state

    def ask(self, cmd: str, expect="# ", timeout=2.0):
        self.link.send(cmd)
        line = self.link.wait_for(("# ERR", expect), timeout=timeout)
        if line is None:
            print("  (no reply from the board)")
        elif line.startswith("# ERR"):
            kv = parse_kv(line)
            print(f"  refused: {kv.get('reason', line)}")
            return None
        return line

    def pos(self) -> int:
        try:
            return int(self.state.get("pos", 0))
        except ValueError:
            return 0

    # -- the jog loop -----------------------------------------------------
    def run_jog(self, max_legs=None):
        if self.a is None or self.b is None:
            print("  set both positions first (e.g.  -200 200 )")
            return
        if self.a == self.b:
            print("  the two positions are the same")
            return
        # Unbounded jogging is stopped by a keypress, which needs a real
        # terminal. Driven from a pipe or a script there is nobody to press a
        # key, so insist on a leg count rather than running until killed.
        if max_legs is None and not sys.stdin.isatty():
            print("  stdin is not a terminal, so there is no key to stop on.")
            print("  give a number of legs instead:  go 6   (or --max-legs 6)")
            return
        if self.ask(f"JOG {self.a} {self.b}", "# JOG_SET") is None:
            return
        if self.ask("START", "# START") is None:
            return

        span = abs(self.b - self.a)
        secs = span / max(self.rate, 1)
        print(f"\n  Jogging {self.a} <-> {self.b}  ({mm(self.a)} to {mm(self.b)}, "
              f"{span} steps = {span/FULL_STEPS_PER_MM:.2f} mm, ~{secs:.1f} s per leg)")
        print(f"  stopping after {max_legs} legs" if max_legs
              else "  press any key to stop")
        print()

        legs = 0
        pos = self.pos()
        mode = "JOGGING"
        t0 = time.monotonic()
        try:
            with KeyWatcher() as keys:
                while True:
                    for line in self.link.drain():
                        kv = parse_kv(line)
                        tag = kv.get("_tag", "")
                        if "pos" in kv:
                            try:
                                pos = int(kv["pos"])
                            except ValueError:
                                pass
                        if tag == "REACHED":
                            legs += 1
                        elif tag in ("STOP", "ARRIVED"):
                            mode = "IDLE"
                        elif tag == "ERR":
                            print(f"\n  board reported: {kv.get('reason', line)}")
                            mode = "IDLE"
                    bar = self._bar(pos)
                    sys.stdout.write(
                        f"\r  {pos:+6d} steps  {mm(pos)}  {bar}  "
                        f"legs {legs:3d}  {time.monotonic()-t0:6.1f} s   ")
                    sys.stdout.flush()
                    if mode == "IDLE":
                        break
                    if max_legs is not None and legs >= max_legs:
                        break
                    if keys.pressed():
                        break
                    time.sleep(0.05)
        except KeyboardInterrupt:
            pass
        finally:
            self.link.send("STOP")
            time.sleep(0.15)
            self.refresh()
        print(f"\n  stopped at {self.pos()} steps ({mm(self.pos())}), {legs} legs\n")

    def run_calibration(self, seconds=30, dwell_ms=2000):
        """One microstep forward, hold, one back, hold - for `seconds`.

        This is the smallest commanded motion the rig can make, held still long
        enough to be measured. A +1/-1 pair returns to the same commanded
        position, so the encoder reading at the end of each pair is the stage's
        error rather than the command's: lost motion, stiction and the
        dead-band around a direction reversal all show up here and nowhere
        else.

        The microstep size is whatever `mres` is set to - that is the point of
        being able to change it. At mres 16 one microstep is 0.625 um; at
        mres 1 it is 10 um.
        """
        mres = int(self.state.get("mres", 16) or 16)
        micro_um = 10.0 / mres
        print(f"\n  Microstep calibration")
        print(f"    microstep   1/{mres} full step = {micro_um:.4f} um")
        print(f"    pattern     +1, hold {dwell_ms} ms, -1, hold {dwell_ms} ms")
        print(f"    duration    {seconds} s  "
              f"(~{int(seconds * 1000 / (2 * dwell_ms))} full cycles)")
        if self.ask(f"CALIB {seconds} {dwell_ms}", "# CALIB_START") is None:
            return
        print("    press any key to stop early\n")

        rows, cycles, enc0 = [], 0, None
        t0 = time.monotonic()
        try:
            with KeyWatcher() as keys:
                while True:
                    done = False
                    for line in self.link.drain():
                        kv = parse_kv(line)
                        tag = kv.get("_tag", "")
                        if tag == "CALIB":
                            enc = int(kv.get("enc", 0))
                            if enc0 is None:
                                enc0 = enc
                            cycles = int(kv.get("cycle", 0))
                            rows.append((kv.get("", ""), enc - enc0))
                            sys.stdout.write(
                                f"\r    cycle {cycles:4d}   encoder "
                                f"{enc - enc0:+7d} counts   "
                                f"{time.monotonic()-t0:6.1f} s    ")
                            sys.stdout.flush()
                        elif tag in ("CALIB_DONE", "STOP"):
                            done = True
                        elif tag == "ERR":
                            print(f"\n  board reported: {kv.get('reason', line)}")
                            done = True
                    if done or keys.pressed():
                        break
                    time.sleep(0.05)
        except KeyboardInterrupt:
            pass
        finally:
            self.link.send("STOP")
            time.sleep(0.2)
            self.refresh()
        print(f"\n  stopped after {cycles} cycles\n")

    def _bar(self, pos, width=30):
        lo, hi = min(self.a, self.b), max(self.a, self.b)
        if hi == lo:
            return ""
        frac = min(max((pos - lo) / (hi - lo), 0.0), 1.0)
        i = int(frac * (width - 1))
        return "[" + "-" * i + "o" + "-" * (width - 1 - i) + "]"

    # -- prompt -----------------------------------------------------------
    def banner(self):
        s = self.state
        print("\n" + "=" * 72)
        print(f"  ESP32-S3 / TMC2209 stage jog     port {self.link.port}")
        print("=" * 72)
        print(f"  position   {self.pos():+d} steps  ({mm(self.pos())})")
        print(f"  speed      {self.rate} full steps/s  "
              f"({self.rate/FULL_STEPS_PER_MM:.2f} mm/s)")
        print(f"  mres       1/{s.get('mres', '?')}        "
              f"current {s.get('current_ma', '?')} mA")
        print(f"  soft limit +/-{self.limit} steps  "
              f"(+/-{self.limit/FULL_STEPS_PER_MM:.2f} mm about zero)")
        print(f"  A {self.a if self.a is not None else '-'}   "
              f"B {self.b if self.b is not None else '-'}")
        print("=" * 72)
        print("  'help' for commands. Nothing moves until you say 'go'.")

    def set_pos(self, which, value):
        if self.limit and abs(value) > self.limit:
            print(f"  {value} is outside the soft limit of +/-{self.limit} "
                  f"({mm(value)}). Widen it with 'limit <steps>' if that is "
                  f"really where you want to go.")
            return
        setattr(self, which, value)
        print(f"  {which.upper()} = {value} steps  ({mm(value)})")

    def loop(self):
        self.banner()
        while True:
            try:
                raw = input("\njog> ").strip()
            except (EOFError, KeyboardInterrupt):
                print()
                return
            if not raw:
                continue
            parts = raw.split()
            verb = parts[0].lower()
            args = parts[1:]

            # two bare numbers: set both positions
            nums = self._as_ints(parts)
            if nums is not None and len(nums) == 2:
                self.set_pos("a", nums[0])
                self.set_pos("b", nums[1])
                continue
            if nums is not None and len(nums) == 1:
                print("  one number is ambiguous - use 'a <n>', 'b <n>' or "
                      "'moveto <n>', or give both positions at once")
                continue

            if verb in ("q", "quit", "exit"):
                return
            if verb in ("help", "?", "h"):
                print(HELP.format(rate=self.rate,
                                  mmps=self.rate / FULL_STEPS_PER_MM,
                                  spm=FULL_STEPS_PER_MM))
                continue
            if verb == "status":
                self.refresh()
                self.banner()
                continue
            if verb == "calib":
                secs = self._int(args[0]) if len(args) >= 1 else 30
                dwell = self._int(args[1]) if len(args) >= 2 else 2000
                if secs is None or dwell is None:
                    continue
                if not sys.stdin.isatty():
                    print("  'calib' needs a terminal to stop on; it runs for "
                          "the full duration otherwise")
                self.run_calibration(secs, dwell)
                continue
            if verb == "go":
                legs = self._int(args[0]) if len(args) == 1 else None
                if len(args) == 1 and legs is None:
                    continue
                self.run_jog(max_legs=legs)
                continue
            if verb == "zero":
                if self.ask("ZERO", "# ZERO_OK"):
                    self.refresh()
                    print(f"  zeroed; position is now {self.pos()} steps")
                continue
            if verb in ("a", "b") and len(args) == 1:
                v = self._int(args[0])
                if v is not None:
                    self.set_pos(verb, v)
                continue
            if verb == "moveto" and len(args) == 1:
                v = self._int(args[0])
                if v is None:
                    continue
                if self.ask(f"MOVETO {v}", "# MOVETO") is None:
                    continue
                self._wait_arrive(v)
                continue
            if verb == "speed" and len(args) == 1:
                v = self._int(args[0])
                if v is None or v <= 0:
                    print("  speed must be a positive number of full steps/s")
                    continue
                if self.ask(f"RATE {v * 1000}", "# RATE_OK"):
                    self.rate = v
                    print(f"  speed = {v} full steps/s "
                          f"({v/FULL_STEPS_PER_MM:.2f} mm/s)")
                continue
            if verb == "limit" and len(args) == 1:
                v = self._int(args[0])
                if v is None or v < 0:
                    continue
                if self.ask(f"LIMIT {v}", "# LIMIT_OK"):
                    self.limit = v
                    print(f"  soft limit = +/-{v} steps "
                          f"(+/-{v/FULL_STEPS_PER_MM:.2f} mm)")
                continue
            if verb == "mres" and len(args) == 1:
                v = self._int(args[0])
                if v is None:
                    continue
                if not (1 <= v <= 256 and v & (v - 1) == 0):
                    print("  mres must be a power of two from 1 to 256 "
                          "(1, 2, 4, 8, 16, 32, 64, 128, 256)")
                    continue
                if self.ask(f"MRES {v}", "# MRES_OK"):
                    self.refresh()
                    print(f"  mres = 1/{self.state.get('mres', v)}")
                continue
            if verb == "current" and len(args) == 1:
                v = self._int(args[0])
                if v is None:
                    continue
                if not 50 <= v <= 1200:
                    print("  current must be 50-1200 mA RMS")
                    continue
                if self.ask(f"CURRENT {v}", "# CURRENT_OK"):
                    print(f"  current = {v} mA RMS")
                continue
            print(f"  don't know '{raw}'. 'help' lists the commands.")

    def _wait_arrive(self, target):
        print(f"  moving to {target} ({mm(target)}) ...")
        t0 = time.monotonic()
        while time.monotonic() - t0 < 300:
            done = False
            for line in self.link.drain():
                kv = parse_kv(line)
                if kv.get("_tag") in ("ARRIVED", "STOP"):
                    done = True
                if kv.get("_tag") == "ERR":
                    print(f"  board reported: {kv.get('reason', line)}")
                    done = True
                if "pos" in kv:
                    sys.stdout.write(f"\r    at {kv['pos']:>7} steps   ")
                    sys.stdout.flush()
            if done:
                break
            time.sleep(0.05)
        self.refresh()
        print(f"\r  at {self.pos()} steps ({mm(self.pos())})      ")

    @staticmethod
    def _int(text):
        try:
            return int(text)
        except ValueError:
            print(f"  '{text}' is not a whole number")
            return None

    @staticmethod
    def _as_ints(parts):
        try:
            return [int(p) for p in parts]
        except ValueError:
            return None


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", help="serial port, e.g. COM11; otherwise searched for")
    ap.add_argument("--simulate", action="store_true",
                    help="no hardware: run against a model of the firmware")
    ap.add_argument("--a", type=int, help="position A in full steps")
    ap.add_argument("--b", type=int, help="position B in full steps")
    ap.add_argument("--speed", type=int, default=DEFAULT_RATE_STEPS_S,
                    help=f"full-step rate [steps/s] (default {DEFAULT_RATE_STEPS_S})")
    ap.add_argument("--limit", type=int, default=DEFAULT_LIMIT,
                    help=f"soft limit +/- about zero [full steps] (default {DEFAULT_LIMIT})")
    ap.add_argument("--go", action="store_true", help="start jogging immediately")
    ap.add_argument("--calib", nargs="?", const=30, type=int, metavar="SECONDS",
                    help="run the microstep calibration oscillation and exit")
    ap.add_argument("--calib-dwell", type=int, default=2000, metavar="MS",
                    help="dwell at each microstep during --calib (default 2000)")
    ap.add_argument("--mres", type=int,
                    help="set microstep resolution before anything else")
    ap.add_argument("--no-prompt", action="store_true",
                    help="exit after --go instead of opening the prompt")
    ap.add_argument("--max-legs", type=int,
                    help="stop after this many legs instead of on a keypress")
    ap.add_argument("--list", action="store_true", help="list serial ports and exit")
    ap.add_argument("--flash", action="store_true",
                    help="flash the firmware before connecting")
    ap.add_argument("--flash-only", action="store_true",
                    help="flash and exit, without opening the prompt")
    ap.add_argument("--flash-mode", choices=("auto", "prebuilt", "source"),
                    default="auto",
                    help="auto tries the prebuilt image then compiles; "
                         "source always compiles")
    ap.add_argument("--allow-download", action="store_true",
                    help="permit fetching arduino-cli (~18 MB) into tools/")
    ap.add_argument("--verbose", action="store_true", help="echo the serial traffic")
    args = ap.parse_args()

    # pyserial is only needed to touch real hardware; --simulate must work
    # without it so the interface can be tried on any machine.
    if not args.simulate:
        try:
            import serial  # noqa: F401
        except ImportError:
            print("pyserial is missing.  pip install pyserial")
            print("(or run with --simulate to try the interface without it)")
            return 2

    if args.list:
        ports = list_ports()
        if not ports:
            print("No serial ports found.")
        for p in ports:
            print(describe(p))
        return 0

    if (args.flash or args.flash_only) and not args.simulate:
        port = args.port
        if not port:
            # Flashing needs a specific port; pick the one that looks right
            # rather than writing firmware to a guess.
            cands = [p.device for p in list_ports()
                     if getattr(p, "vid", None) in USB_VIDS]
            if "COM11" in [p.device for p in list_ports()]:
                port = "COM11"
            elif len(cands) == 1:
                port = cands[0]
            else:
                print("Which port is the board on? Pass --port; "
                      "--list shows what is attached.")
                for p in list_ports():
                    print("  " + describe(p))
                return 2
            print(f"No --port given; flashing {port}")
        if not do_flash(port, args.flash_mode, args.allow_download):
            print("\nFlashing did not succeed - see the log above.")
            return 1
        if args.flash_only:
            return 0
        # The board reboots into the new firmware; let USB re-enumerate.
        print("  waiting for the board to come back ...")
        time.sleep(3.0)
        if not args.port:
            args.port = port

    if args.simulate:
        print("SIMULATION - no hardware. The stage is a model; nothing will move.")
        link = Link(port="simulated", ser=FakeSerial(), verbose=args.verbose)
        link.start()
        time.sleep(0.3)
        link.send("PING")
        link.wait_for("# PONG", timeout=2.0)
    else:
        link = find_board(args.port, verbose=args.verbose)
        if link is None:
            print("\nTry --list to see the ports, or --simulate to try the "
                  "interface without hardware.")
            return 1

    try:
        ready = link.wait_for("# READY", timeout=0.4)
        if ready and parse_kv(ready).get("ready") == "0":
            print("\nThe board says its driver did not come up (TMC2209 UART).")
            print("Check the UART wiring on GPIO17/18 and the motor supply;")
            print("moves will be refused until it does.")

        console = Console(link, limit=args.limit, rate=args.speed)
        if args.limit != DEFAULT_LIMIT:
            console.ask(f"LIMIT {args.limit}", "# LIMIT_OK")
        if args.speed != DEFAULT_RATE_STEPS_S:
            console.ask(f"RATE {args.speed * 1000}", "# RATE_OK")
        if args.mres:
            if not (1 <= args.mres <= 256 and args.mres & (args.mres - 1) == 0):
                print("  mres must be a power of two from 1 to 256")
                return 2
            console.ask(f"MRES {args.mres}", "# MRES_OK")
        console.refresh()
        if args.calib:
            console.run_calibration(args.calib, args.calib_dwell)
            if args.no_prompt or not sys.stdin.isatty():
                return 0
        if args.a is not None:
            console.set_pos("a", args.a)
        if args.b is not None:
            console.set_pos("b", args.b)
        if args.go:
            console.run_jog(max_legs=args.max_legs)
            # A scripted invocation has nobody to answer the prompt, so falling
            # into it would just block on input() forever.
            if args.no_prompt or not sys.stdin.isatty():
                return 0
        console.loop()
    finally:
        try:
            link.send("STOP")
            time.sleep(0.1)
        except Exception:
            pass
        link.close()
    print("closed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
