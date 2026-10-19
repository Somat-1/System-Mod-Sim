#!/usr/bin/env python3
"""Live monitor / recorder for the ESP32-S3 quadrature encoder stream.

Pairs with firmware/esp32s3_encoder_stream. Opens a window with live plots of
position, velocity and the raw A/B/Z logic levels; lets you record a segment
and export it to CSV (plus a PNG of the current view).

Usage:
    python encoder_monitor.py                       # pick the port in the GUI
    python encoder_monitor.py --port /dev/cu.usbmodem3101 --cpr 4000
    python encoder_monitor.py --simulate            # no hardware, synthetic encoder
    python encoder_monitor.py --headless --port ... --duration 10 --out run.csv
"""

import argparse
import csv
import datetime as dt
import math
import os
import sys
import threading
import time
import traceback
from pathlib import Path


def _bootstrap():
    """Re-launch under the project venv / Homebrew expat if this interpreter can't run the GUI.

    Homebrew's Python on this Mac links pyexpat against an older system libexpat
    (missing _XML_SetAllocTrackerActivationThreshold), which breaks matplotlib.
    Homebrew's own expat has the symbol, so point the dynamic loader at it.
    pyserial lives in "Stage & motor/motor identification/.venv".
    """
    if os.environ.get("ENCODER_MONITOR_BOOTSTRAPPED"):
        return
    env = dict(os.environ, ENCODER_MONITOR_BOOTSTRAPPED="1")
    exe = sys.executable
    relaunch = False
    try:
        import pyexpat  # noqa: F401
    except ImportError:
        brew_expat = "/opt/homebrew/opt/expat/lib"
        if os.path.isdir(brew_expat):
            env["DYLD_LIBRARY_PATH"] = brew_expat + os.pathsep + env.get("DYLD_LIBRARY_PATH", "")
            relaunch = True
    venv_py = Path(__file__).resolve().parent.parent / ".venv" / "bin" / "python"
    if venv_py.exists() and Path(sys.prefix).resolve() != venv_py.parent.parent.resolve():
        exe = str(venv_py)
        relaunch = True
    if relaunch:
        os.execve(exe, [exe, os.path.abspath(__file__), *sys.argv[1:]], env)


if __name__ == "__main__":
    _bootstrap()

import numpy as np

try:
    import serial
    import serial.tools.list_ports
except ImportError:  # allow --simulate without pyserial
    serial = None

HERE = Path(__file__).resolve().parent
RECORDINGS_DIR = HERE / "recordings"

CSV_COLUMNS = ["host_time_s", "seq", "t_us", "count", "A", "B", "Z", "Z_edge", "z_count", "z_pos"]


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

def parse_line(line, host_t):
    """'D,seq,t_us,count,levels,z_count,z_pos' -> sample tuple (CSV_COLUMNS order)."""
    parts = line.split(",")
    if len(parts) != 7 or parts[0] != "D":
        return None
    try:
        seq, t_us, count, levels, zc, zp = (int(p) for p in parts[1:])
    except ValueError:
        return None
    return (host_t, seq, t_us, count,
            levels & 1, (levels >> 1) & 1, (levels >> 2) & 1, (levels >> 3) & 1, zc, zp)


class SerialSource:
    def __init__(self, port, on_sample, on_info):
        if serial is None:
            raise RuntimeError("pyserial is not installed (pip install pyserial)")
        self.ser = serial.Serial(port, 921600, timeout=0.1)
        self.port = port
        self.on_sample = on_sample
        self.on_info = on_info
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self):
        self.ser.reset_input_buffer()
        self.send("S1")
        self.send("I")
        self._thread.start()

    def send(self, cmd):
        self.ser.write((cmd.strip() + "\n").encode())

    def stop(self):
        self._stop.set()
        self._thread.join(timeout=1.0)
        try:
            self.ser.close()
        except Exception:
            pass

    def _run(self):
        buf = b""
        while not self._stop.is_set():
            try:
                chunk = self.ser.read(self.ser.in_waiting or 1)
            except Exception as exc:
                self.on_info(f"#ERR,serial: {exc}")
                return
            if not chunk:
                continue
            host_t = time.time()
            buf += chunk
            *lines, buf = buf.split(b"\n")
            for raw in lines:
                line = raw.decode(errors="replace").strip()
                if line.startswith("D,"):
                    s = parse_line(line, host_t)
                    if s:
                        self.on_sample(s)
                elif line:
                    self.on_info(line)


class SimSource:
    """Synthetic encoder: oscillating motion with index pulses, for testing the GUI."""

    def __init__(self, on_sample, on_info, cpr=4000, rate=1000):
        self.on_sample, self.on_info = on_sample, on_info
        self.cpr, self.rate = cpr, rate
        self.port = "simulated"
        self.offset = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self):
        self._thread.start()
        self.on_info(f"#INFO,fw=simulator,rate_hz={self.rate},cpr={self.cpr}")

    def send(self, cmd):
        cmd = cmd.strip().upper()
        if cmd.startswith("R"):
            self.rate = max(10, min(5000, int(cmd[1:] or 1000)))
            self.on_info(f"#ACK,rate_hz={self.rate}")
        elif cmd == "Z":
            self.offset = self._pos(time.time() - self.t0)
            self.zc, self.zp = 0, 0
            self.on_info("#ACK,zero")

    def stop(self):
        self._stop.set()
        self._thread.join(timeout=1.0)

    def _pos(self, t):
        return int(round(self.cpr * (1.5 * math.sin(2 * math.pi * 0.25 * t) + 0.3 * t)))

    def _run(self):
        self.t0 = time.time()
        seq, prev, self.zc, self.zp = 0, 0, 0, 0
        next_t = self.t0
        while not self._stop.is_set():
            now = time.time()
            while next_t <= now:
                t = next_t - self.t0
                raw = self._pos(t)
                phase = raw % 4
                a = 1 if phase in (1, 2) else 0
                b = 1 if phase in (2, 3) else 0
                edge = 0
                if raw // self.cpr != prev // self.cpr:
                    self.zc += 1
                    self.zp = raw - self.offset
                    edge = 1
                prev = raw
                z = 1 if raw % self.cpr == 0 else 0
                self.on_sample((now, seq, int(t * 1e6), raw - self.offset, a, b, z, edge, self.zc, self.zp))
                seq += 1
                next_t += 1.0 / self.rate
            time.sleep(0.005)


# --------------------------------------------------------------------------
# Buffering
# --------------------------------------------------------------------------

class Ring:
    """Fixed-capacity numpy ring buffer of sample rows."""

    def __init__(self, capacity, ncols=len(CSV_COLUMNS)):
        self.data = np.zeros((capacity, ncols), dtype=np.float64)
        self.cap, self.n, self.i = capacity, 0, 0

    def extend(self, rows):
        rows = np.asarray(rows, dtype=np.float64)
        if len(rows) >= self.cap:
            rows = rows[-self.cap:]
        k = len(rows)
        end = self.i + k
        if end <= self.cap:
            self.data[self.i:end] = rows
        else:
            first = self.cap - self.i
            self.data[self.i:] = rows[:first]
            self.data[:k - first] = rows[first:]
        self.i = end % self.cap
        self.n = min(self.cap, self.n + k)

    def latest(self, m=None):
        m = self.n if m is None else min(m, self.n)
        idx = (np.arange(self.i - m, self.i)) % self.cap
        return self.data[idx]

    def clear(self):
        self.n, self.i = 0, 0


class Collector:
    """Thread-safe hand-off between the reader thread and the consumer."""

    def __init__(self):
        self.lock = threading.Lock()
        self.pending, self.info = [], []

    def add_sample(self, s):
        with self.lock:
            self.pending.append(s)

    def add_info(self, line):
        with self.lock:
            self.info.append(line)

    def drain(self):
        with self.lock:
            p, i = self.pending, self.info
            self.pending, self.info = [], []
        return p, i


def write_csv(path, rows, meta):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        for k, v in meta.items():
            f.write(f"# {k}: {v}\n")
        w = csv.writer(f)
        w.writerow(CSV_COLUMNS)
        for r in rows:
            w.writerow([f"{r[0]:.6f}", *r[1:]])
    return path


def default_filename(prefix="encoder"):
    return RECORDINGS_DIR / f"{prefix}_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"


def measured_cpr(z_positions):
    """Counts between consecutive index pulses (median), or None."""
    if len(z_positions) < 2:
        return None
    d = np.abs(np.diff(z_positions))
    d = d[d > 0]
    return float(np.median(d)) if len(d) else None


# --------------------------------------------------------------------------
# Headless recording
# --------------------------------------------------------------------------

def run_headless(args):
    col = Collector()
    src = SimSource(col.add_sample, col.add_info, rate=args.rate) if args.simulate \
        else SerialSource(args.port, col.add_sample, col.add_info)
    src.start()
    src.send(f"R{args.rate}")
    time.sleep(0.3)
    col.drain()
    if args.zero:
        src.send("Z")
        time.sleep(0.1)
        col.drain()
    rows, t_end = [], time.time() + args.duration
    while time.time() < t_end:
        time.sleep(0.05)
        p, info = col.drain()
        rows.extend(p)
        for line in info:
            print(line)
    src.stop()
    out = Path(args.out) if args.out else default_filename()
    meta = {"source": src.port, "rate_hz": args.rate, "cpr": args.cpr or "unknown",
            "start": dt.datetime.fromtimestamp(rows[0][0]).isoformat() if rows else "",
            "samples": len(rows)}
    write_csv(out, rows, meta)
    if rows:
        seq = np.array([r[1] for r in rows])
        drops = int(np.sum(np.diff(seq) - 1)) if len(seq) > 1 else 0
        print(f"{len(rows)} samples, final count {rows[-1][3]}, index pulses {rows[-1][8]}, "
              f"dropped {drops}")
    print(f"saved {out}")


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------

def run_gui(args):
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    MAX_RATE = 5000
    MAX_WINDOW_S = 60
    PLOT_POINTS = 4000  # decimate plotted traces to roughly this many points
    REDRAW_MS = 100     # 10 FPS is ample for a live instrument display

    class App:
        def __init__(self, root):
            self.root = root
            root.title("Encoder Monitor — Motor Identification v1")
            self.col = Collector()
            self.src = None
            self.ring = Ring(MAX_RATE * MAX_WINDOW_S)
            self.recording = False
            self.rec_rows = []
            self.rec_saved = True
            self.rec_start = None
            self.last_seq = None
            self.drops = 0
            self.z_hist = []
            self.rate_est = 0.0
            self.frozen = False
            self.data_version = 0
            self.drawn_version = -1

            self._build_controls()
            self._build_plot()
            root.protocol("WM_DELETE_WINDOW", self.on_close)
            if args.simulate:
                self.port_var.set("SIMULATE")
                self.connect()
            elif args.port:
                self.port_var.set(args.port)
                self.connect()
            self.root.after(50, self.tick)

        # ---- layout ------------------------------------------------------
        def _build_controls(self):
            bar = ttk.Frame(self.root, padding=6)
            bar.pack(side=tk.TOP, fill=tk.X)

            ttk.Label(bar, text="Port").pack(side=tk.LEFT)
            self.port_var = tk.StringVar()
            self.port_box = ttk.Combobox(bar, textvariable=self.port_var, width=28)
            self.port_box.pack(side=tk.LEFT, padx=4)
            ttk.Button(bar, text="↻", width=3, command=self.refresh_ports).pack(side=tk.LEFT)
            self.conn_btn = ttk.Button(bar, text="Connect", command=self.toggle_connect)
            self.conn_btn.pack(side=tk.LEFT, padx=(4, 12))

            ttk.Label(bar, text="Rate [Hz]").pack(side=tk.LEFT)
            self.rate_var = tk.StringVar(value=str(args.rate))
            ttk.Combobox(bar, textvariable=self.rate_var, width=6,
                         values=["100", "250", "500", "1000", "2000", "5000"]).pack(side=tk.LEFT, padx=4)
            ttk.Button(bar, text="Set", command=self.apply_rate).pack(side=tk.LEFT, padx=(0, 12))

            ttk.Label(bar, text="CPR (x4)").pack(side=tk.LEFT)
            self.cpr_var = tk.StringVar(value=str(args.cpr) if args.cpr else "")
            ttk.Entry(bar, textvariable=self.cpr_var, width=7).pack(side=tk.LEFT, padx=4)
            ttk.Button(bar, text="Use Z", command=self.use_measured_cpr).pack(side=tk.LEFT, padx=(0, 12))

            ttk.Label(bar, text="Window [s]").pack(side=tk.LEFT)
            self.win_var = tk.StringVar(value=str(args.window))
            ttk.Combobox(bar, textvariable=self.win_var, width=4,
                         values=["1", "2", "5", "10", "20", "60"]).pack(side=tk.LEFT, padx=4)

            bar2 = ttk.Frame(self.root, padding=(6, 0, 6, 6))
            bar2.pack(side=tk.TOP, fill=tk.X)
            ttk.Button(bar2, text="Zero", command=self.zero).pack(side=tk.LEFT)
            self.freeze_btn = ttk.Button(bar2, text="Freeze view", command=self.toggle_freeze)
            self.freeze_btn.pack(side=tk.LEFT, padx=4)
            self.rec_btn = ttk.Button(bar2, text="● Record", command=self.toggle_record)
            self.rec_btn.pack(side=tk.LEFT, padx=(12, 4))
            ttk.Button(bar2, text="Export CSV…", command=self.export_csv).pack(side=tk.LEFT, padx=4)
            ttk.Button(bar2, text="Save plot PNG…", command=self.export_png).pack(side=tk.LEFT, padx=4)
            self.rec_label = ttk.Label(bar2, text="not recording")
            self.rec_label.pack(side=tk.LEFT, padx=12)

            stat = ttk.Frame(self.root, padding=(6, 0))
            stat.pack(side=tk.TOP, fill=tk.X)
            self.stat_var = tk.StringVar(value="disconnected")
            ttk.Label(stat, textvariable=self.stat_var, font=("Menlo", 12)).pack(side=tk.LEFT)

            self.info_var = tk.StringVar(value="")
            ttk.Label(self.root, textvariable=self.info_var, foreground="#666",
                      padding=(6, 0)).pack(side=tk.BOTTOM, fill=tk.X)
            self.refresh_ports()

        def _build_plot(self):
            # Recomputing constrained layout on every live draw is expensive.
            # Lay the fixed three-axis arrangement out once instead.
            self.fig = Figure(figsize=(11, 7), dpi=100)
            self.ax_pos, self.ax_vel, self.ax_dig = self.fig.subplots(
                3, 1, sharex=True, gridspec_kw={"height_ratios": [3, 2, 2]})
            self.l_pos, = self.ax_pos.plot([], [], lw=1.2, color="#2a6fdb")
            self.l_vel, = self.ax_vel.plot([], [], lw=1.0, color="#d9822b")
            self.l_a, = self.ax_dig.step([], [], where="post", lw=1, color="#2a6fdb", label="A")
            self.l_b, = self.ax_dig.step([], [], where="post", lw=1, color="#3fa34d", label="B")
            self.l_z, = self.ax_dig.step([], [], where="post", lw=1, color="#c0392b", label="Z (edge-latched)")
            self.ax_dig.set_yticks([0.5, 2.0, 3.5], ["Z", "B", "A"])
            self.ax_dig.set_ylim(-0.3, 4.3)
            self.ax_dig.set_xlabel("time [s]")
            for ax in (self.ax_pos, self.ax_vel, self.ax_dig):
                ax.grid(alpha=0.3)
            self.fig.tight_layout()
            canvas = FigureCanvasTkAgg(self.fig, master=self.root)
            NavigationToolbar2Tk(canvas, self.root).update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            self.canvas = canvas

        # ---- connection --------------------------------------------------
        def refresh_ports(self):
            ports = []
            if serial is not None:
                ports = [p.device for p in serial.tools.list_ports.comports()
                         if "Bluetooth" not in p.device and "debug-console" not in p.device]
            self.port_box["values"] = ports + ["SIMULATE"]
            if not self.port_var.get() and ports:
                usb = [p for p in ports if "usbmodem" in p or "usbserial" in p]
                self.port_var.set((usb or ports)[0])

        def toggle_connect(self):
            self.disconnect() if self.src else self.connect()

        def connect(self):
            port = self.port_var.get().strip()
            try:
                if port == "SIMULATE":
                    self.src = SimSource(self.col.add_sample, self.col.add_info,
                                         cpr=args.cpr or 4000, rate=int(self.rate_var.get()))
                else:
                    self.src = SerialSource(port, self.col.add_sample, self.col.add_info)
                self.src.start()
            except Exception as exc:
                self.src = None
                messagebox.showerror("Connect failed", str(exc))
                return
            self.reset_view()
            self.apply_rate()
            self.conn_btn.config(text="Disconnect")

        def disconnect(self):
            if self.recording:
                self.toggle_record()
            if self.src:
                self.src.stop()
            self.src = None
            self.conn_btn.config(text="Connect")
            self.stat_var.set("disconnected")

        def apply_rate(self):
            if self.src:
                try:
                    self.src.send(f"R{int(float(self.rate_var.get()))}")
                except ValueError:
                    pass

        def reset_view(self):
            self.ring.clear()
            self.last_seq = None
            self.drops = 0
            self.z_hist = []

        def zero(self):
            if self.src:
                self.src.send("Z")
                self.z_hist = []

        # ---- record / export ---------------------------------------------
        def toggle_record(self):
            if not self.recording:
                if not self.rec_saved and self.rec_rows and not messagebox.askyesno(
                        "Discard recording?", f"{len(self.rec_rows)} unsaved samples will be discarded."):
                    return
                self.rec_rows, self.rec_saved = [], False
                self.rec_start = dt.datetime.now()
                self.recording = True
                self.rec_btn.config(text="■ Stop")
            else:
                self.recording = False
                self.rec_btn.config(text="● Record")

        def _meta(self):
            return {"source": self.src.port if self.src else "",
                    "rate_hz": self.rate_var.get(),
                    "cpr": self.cpr() or "unknown",
                    "start": self.rec_start.isoformat() if self.rec_start else "",
                    "samples": len(self.rec_rows),
                    "dropped_samples_session": self.drops,
                    "pins": "A=GPIO4 B=GPIO5 Z=GPIO6 (ESP32-S3, PCNT x4)"}

        def export_csv(self):
            if not self.rec_rows:
                messagebox.showinfo("Nothing to export", "Record something first.")
                return
            RECORDINGS_DIR.mkdir(exist_ok=True)
            d = default_filename()
            path = filedialog.asksaveasfilename(initialdir=d.parent, initialfile=d.name,
                                                defaultextension=".csv",
                                                filetypes=[("CSV", "*.csv")])
            if path:
                write_csv(path, self.rec_rows, self._meta())
                self.rec_saved = True
                self.info_var.set(f"saved {path}")

        def export_png(self):
            RECORDINGS_DIR.mkdir(exist_ok=True)
            d = default_filename("encoder_plot").with_suffix(".png")
            path = filedialog.asksaveasfilename(initialdir=d.parent, initialfile=d.name,
                                                defaultextension=".png",
                                                filetypes=[("PNG", "*.png")])
            if path:
                self.fig.savefig(path, dpi=150)
                self.info_var.set(f"saved {path}")

        # ---- helpers -----------------------------------------------------
        def cpr(self):
            try:
                v = float(self.cpr_var.get())
                return v if v > 0 else None
            except ValueError:
                return None

        def use_measured_cpr(self):
            m = measured_cpr(self.z_hist)
            if m:
                self.cpr_var.set(f"{m:.0f}")
            else:
                messagebox.showinfo("No index data", "Need at least two Z pulses (turn the shaft > 1 rev).")

        def toggle_freeze(self):
            self.frozen = not self.frozen
            self.freeze_btn.config(text="Live view" if self.frozen else "Freeze view")

        # ---- periodic update --------------------------------------------
        def tick(self):
            try:
                self._ingest()
                if not self.frozen and self.data_version != self.drawn_version:
                    self._redraw()
                    self.drawn_version = self.data_version
            except Exception as exc:
                # Tk otherwise prints callback failures only to the launching
                # terminal, which looks like an unexplained GUI freeze/crash.
                traceback.print_exc()
                self.info_var.set(f"#ERR,monitor: {type(exc).__name__}: {exc}")
            finally:
                self.root.after(REDRAW_MS, self.tick)

        def _ingest(self):
            rows, info = self.col.drain()
            for line in info:
                self.info_var.set(line)
                if line.startswith("#ACK,zero"):
                    self.reset_view()
            if not rows:
                return
            seqs = np.fromiter((r[1] for r in rows), dtype=np.int64, count=len(rows))
            if self.last_seq is not None:
                gaps = np.diff(np.concatenate(([self.last_seq], seqs))) - 1
                self.drops += int(np.sum(gaps[gaps > 0]))
            self.last_seq = int(seqs[-1])
            prev_zc = self.ring.latest(1)[0, 8] if self.ring.n else None
            self.ring.extend(rows)
            self.data_version += 1
            for r in rows:
                if prev_zc is not None and r[8] != prev_zc:
                    self.z_hist.append(r[9])
                    self.z_hist = self.z_hist[-50:]
                prev_zc = r[8]
            if self.recording:
                self.rec_rows.extend(rows)

        def _redraw(self):
            if self.ring.n < 2:
                return
            try:
                win = float(self.win_var.get())
            except ValueError:
                win = 5.0
            # Copy only enough recent rows to cover the requested time window.
            # The old implementation copied the entire 300k-row ring at every
            # redraw and then discarded almost all of it. After a few minutes
            # that caused ~1 GB RSS, 100% CPU, reader backlog, and GUI failure.
            try:
                configured_rate = max(1.0, float(self.rate_var.get()))
            except ValueError:
                configured_rate = 1000.0
            display_rate = min(MAX_RATE, max(configured_rate, self.rate_est, 10.0))
            rows_needed = min(
                self.ring.n,
                max(2, int(math.ceil(win * display_rate * 1.25)) + 2),
            )
            data = self.ring.latest(rows_needed)
            t = data[:, 2] * 1e-6
            data = data[t >= t[-1] - win]
            t = data[:, 2] * 1e-6
            if len(t) < 2:
                return
            dt_med = float(np.median(np.diff(t[-200:])))
            self.rate_est = 1.0 / dt_med if dt_med > 0 else 0.0
            t_rel = t - t[-1]
            count = data[:, 3]
            cpr = self.cpr()

            # velocity over ~10 ms (at least one sample)
            k = max(1, int(round(0.01 * self.rate_est)))
            if len(count) > k:
                v_cps = (count[k:] - count[:-k]) / np.maximum(t[k:] - t[:-k], 1e-9)
                v_t = t_rel[k:]
            else:
                v_cps, v_t = np.zeros(1), t_rel[-1:]

            step = max(1, len(t_rel) // PLOT_POINTS)
            if cpr:
                self.l_pos.set_data(t_rel[::step], count[::step] / cpr * 360.0)
                self.ax_pos.set_ylabel("angle [deg]")
                self.l_vel.set_data(v_t[::step], v_cps[::step] / cpr * 60.0)
                self.ax_vel.set_ylabel("speed [rpm]")
            else:
                self.l_pos.set_data(t_rel[::step], count[::step])
                self.ax_pos.set_ylabel("position [counts]")
                self.l_vel.set_data(v_t[::step], v_cps[::step])
                self.ax_vel.set_ylabel("speed [counts/s]")

            # digital traces: no decimation of the last ~PLOT_POINTS samples would
            # hide edges, so show them full-resolution over the visible window.
            dstep = 1 if len(t_rel) <= 4 * PLOT_POINTS else step
            td = t_rel[::dstep]
            self.l_a.set_data(td, data[::dstep, 4] + 3.0)
            self.l_b.set_data(td, data[::dstep, 5] + 1.5)
            z_vis = np.maximum(data[:, 6], data[:, 7])
            if dstep > 1:  # keep short Z pulses visible when decimating
                n = len(z_vis) // dstep * dstep
                z_vis = z_vis[:n].reshape(-1, dstep).max(axis=1) if n else z_vis[:1]
                td = td[:len(z_vis)]
            self.l_z.set_data(td, z_vis)

            self.ax_dig.set_xlim(-win, 0)
            for ax in (self.ax_pos, self.ax_vel):
                ax.relim()
                ax.autoscale_view(scalex=False)
            self.canvas.draw_idle()

            last = data[-1]
            pos_txt = f"{int(last[3]):>9d} cnt"
            if cpr:
                pos_txt += f"  {last[3] / cpr:+9.3f} rev  {(last[3] % cpr) / cpr * 360:7.2f}°"
            vel = v_cps[-1]
            vel_txt = f"{vel / cpr * 60:+8.1f} rpm" if cpr else f"{vel:+9.0f} cnt/s"
            m = measured_cpr(self.z_hist)
            z_txt = f"Z pulses {int(last[8])}" + (f" (CPR≈{m:.0f})" if m else "")
            self.stat_var.set(
                f"pos {pos_txt} | vel {vel_txt} | A{int(last[4])} B{int(last[5])} Z{int(last[6])} | "
                f"{z_txt} | {self.rate_est:6.0f} Hz | dropped {self.drops}")
            if self.recording:
                secs = (dt.datetime.now() - self.rec_start).total_seconds()
                self.rec_label.config(text=f"REC {secs:6.1f} s  {len(self.rec_rows)} samples")
            elif self.rec_rows:
                self.rec_label.config(text=f"{len(self.rec_rows)} samples recorded"
                                           f"{'' if self.rec_saved else ' (unsaved)'}")

        def on_close(self):
            if self.rec_rows and not self.rec_saved and not messagebox.askyesno(
                    "Quit?", f"{len(self.rec_rows)} recorded samples are not exported. Quit anyway?"):
                return
            self.disconnect()
            self.root.destroy()

    root = tk.Tk()
    App(root)
    root.mainloop()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", help="serial port, e.g. /dev/cu.usbmodem3101")
    ap.add_argument("--simulate", action="store_true", help="use a synthetic encoder")
    ap.add_argument("--rate", type=int, default=1000, help="sample rate requested from the ESP32 [Hz]")
    ap.add_argument("--cpr", type=float, default=0, help="counts per revolution after x4 decoding (= 4 x PPR)")
    ap.add_argument("--window", type=float, default=5, help="plot window [s]")
    ap.add_argument("--headless", action="store_true", help="record without GUI and exit")
    ap.add_argument("--duration", type=float, default=10, help="headless record duration [s]")
    ap.add_argument("--zero", action="store_true", help="headless: zero the count before recording")
    ap.add_argument("--out", help="headless: output CSV path (default recordings/encoder_<timestamp>.csv)")
    args = ap.parse_args()
    if args.headless:
        if not (args.port or args.simulate):
            ap.error("--headless needs --port or --simulate")
        run_headless(args)
    else:
        run_gui(args)


if __name__ == "__main__":
    main()
