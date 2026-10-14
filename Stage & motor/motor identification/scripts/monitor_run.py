#!/usr/bin/env python3
"""Live progress window for a motor-identification recording.

Follows the CSV that run_campaign.py is writing and shows, refreshed twice a
second:

* overall progress through the plan, the current block and action, and the
  time elapsed and remaining;
* commanded and measured shaft angle over the last few seconds;
* the tracking error between them, in degrees of shaft rotation.

It only reads the recording, never the serial port, so it cannot disturb the
run and can be closed and reopened at any time.

    python scripts/monitor_run.py                  # newest recording
    python scripts/monitor_run.py recordings/motor_identification_<stamp>.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import deque
from pathlib import Path

import tkinter as tk
from tkinter import ttk

ROOT = Path(__file__).resolve().parent.parent
PLAN = ROOT / "sequence" / "motor_identification_plan.csv"
FULL_STEPS_PER_REV = 200
U256 = 256
WINDOW_S = 8.0                 # seconds of trace shown


class Recording:
    """Incrementally parses the recording as run_campaign.py appends to it."""

    def __init__(self, path: Path, cpr: int, plan: list[dict]):
        self.path = path
        self.cpr = cpr
        self.plan = plan
        self.total_s = float(plan[-1]["t_end_s"]) if plan else 0.0
        self.offset = 0
        self.partial = ""
        self.t0_us = None
        self.samples = deque(maxlen=40000)     # (t_s, cmd_deg, enc_deg)
        self.action = -1
        self.tag = ""
        self.state = "waiting for RUN_START"
        self.last_t = 0.0
        self.max_abs_err = 0.0

    def cmd_deg(self, u256: int) -> float:
        return u256 / U256 / FULL_STEPS_PER_REV * 360.0

    def enc_deg(self, counts: int) -> float:
        return counts / self.cpr * 360.0

    def poll(self):
        try:
            with self.path.open("r", encoding="utf-8", errors="replace") as f:
                f.seek(self.offset)
                chunk = f.read()
                self.offset = f.tell()
        except FileNotFoundError:
            return
        if not chunk:
            return
        text = self.partial + chunk
        lines = text.split("\n")
        self.partial = lines.pop()             # keep an incomplete last line
        for line in lines:
            self._line(line.strip())

    def _line(self, line: str):
        if line.startswith("D,"):
            p = line.split(",")
            if len(p) != 10:
                return
            try:
                t_us, cmd, enc = int(p[2]), int(p[5]), int(p[6])
            except ValueError:
                return
            if self.t0_us is None:
                self.t0_us = t_us
            t = (t_us - self.t0_us) / 1e6
            self.last_t = t
            c, e = self.cmd_deg(cmd), self.enc_deg(enc)
            self.samples.append((t, c, e))
            self.max_abs_err = max(self.max_abs_err, abs(e - c))
        elif line.startswith("# ACTION_BEGIN"):
            m = re.search(r"index=(\d+),tag=([^,]+)", line)
            if m:
                self.action, self.tag = int(m.group(1)), m.group(2)
            self.state = "running"
        elif line.startswith("# RUN_START"):
            self.state = "running"
        elif line.startswith("# RUN_"):
            self.state = line[2:].split(",")[0]       # RUN_COMPLETE / ABORTED / FAILED
        elif line.startswith("# ERR"):
            self.state = line[2:]

    def block(self) -> str:
        if 0 <= self.action < len(self.plan):
            row = self.plan[self.action]
            return f"{row['block']} - {row['label']}"
        return ""


def newest_recording() -> Path | None:
    files = sorted((p for p in (ROOT / "recordings").glob("*.csv")
                    if not p.stem.endswith(("_steps", "_summary"))),
                   key=lambda p: p.stat().st_mtime)
    return files[-1] if files else None


PANEL_BG = "#ffffff"
GRID = "#e4e4e4"
INK = "#333333"
C_CMD = "#999999"
C_ENC = "#2f6db5"
C_ERR = "#c0392b"


class Plot:
    """A minimal line plot on a tkinter Canvas: axes, grid, labelled ranges."""

    def __init__(self, parent, height, ylabel):
        self.c = tk.Canvas(parent, height=height, bg=PANEL_BG, highlightthickness=0)
        self.c.pack(fill="both", expand=True, padx=8, pady=4)
        self.ylabel = ylabel

    def draw(self, xs, series, xlim, ylim):
        c = self.c
        c.delete("all")
        w, h = max(c.winfo_width(), 200), max(c.winfo_height(), 80)
        l, r, t, b = 70, 12, 10, 24
        x0, x1 = xlim
        y0, y1 = ylim

        def px(x):
            return l + (x - x0) / (x1 - x0) * (w - l - r)

        def py(y):
            return h - b - (y - y0) / (y1 - y0) * (h - t - b)

        for i in range(5):
            y = y0 + (y1 - y0) * i / 4
            c.create_line(l, py(y), w - r, py(y), fill=GRID)
            c.create_text(l - 6, py(y), text=f"{y:.3g}", anchor="e", fill=INK,
                          font=("Helvetica", 10))
        for i in range(5):
            x = x0 + (x1 - x0) * i / 4
            c.create_text(px(x), h - b + 12, text=f"{x:.1f} s", fill=INK,
                          font=("Helvetica", 10))
        c.create_text(8, t, text=self.ylabel, anchor="nw", fill=INK,
                      font=("Helvetica", 10, "bold"))
        for ys, colour in series:
            pts = []
            for x, y in zip(xs, ys):
                pts += [px(x), py(min(max(y, y0), y1))]
            if len(pts) >= 4:
                c.create_line(*pts, fill=colour, width=1.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("recording", nargs="?", type=Path)
    ap.add_argument("--cpr", type=int, default=20000,
                    help="encoder counts per shaft revolution (bring-up measured 20000)")
    ap.add_argument("--plan", type=Path, default=PLAN,
                    help="the plan being run (default: the full campaign)")
    args = ap.parse_args()

    path = args.recording or newest_recording()
    if path is None:
        print("no recording found in recordings/; start run_campaign.py first")
        return 2
    with args.plan.open(newline="", encoding="utf-8") as f:
        plan = list(csv.DictReader(f))
    rec = Recording(path, args.cpr, plan)

    root = tk.Tk()
    root.title(f"Motor identification - {path.name}")
    root.geometry("1000x720")
    root.configure(bg=PANEL_BG)
    head = tk.Label(root, text="", anchor="w", bg=PANEL_BG, fg=INK,
                    font=("Helvetica", 14, "bold"))
    head.pack(fill="x", padx=10, pady=(10, 2))
    bar = ttk.Progressbar(root, maximum=1000)
    bar.pack(fill="x", padx=10, pady=4)
    info = tk.Label(root, text="", anchor="w", bg=PANEL_BG, fg=INK,
                    font=("Menlo", 12))
    info.pack(fill="x", padx=10, pady=(0, 6))
    legend = tk.Label(root, text="grey: commanded    blue: encoder    "
                      f"(last {WINDOW_S:.0f} s)", anchor="w", bg=PANEL_BG,
                      fg=INK, font=("Helvetica", 11))
    legend.pack(fill="x", padx=10)
    angle = Plot(root, 330, "shaft angle [deg]")
    error = Plot(root, 200, "encoder - command [deg]")

    def update():
        rec.poll()
        frac = min(rec.last_t / rec.total_s, 1.0) if rec.total_s else 0.0
        if rec.state == "RUN_COMPLETE":
            frac = 1.0
        bar["value"] = frac * 1000
        remaining = max(rec.total_s - rec.last_t, 0.0)
        head.config(text=f"{rec.state}    action {max(rec.action, 0) + 1}/{len(plan)}"
                         f"    {rec.block()}",
                    fg=C_ERR if rec.state not in ("running", "RUN_COMPLETE",
                                                  "waiting for RUN_START") else INK)
        info.config(text=f"{frac * 100:5.1f} %   elapsed {rec.last_t / 60:5.1f} min   "
                         f"remaining ~{remaining / 60:4.1f} min   "
                         f"max |error| {rec.max_abs_err:.3f} deg")
        if rec.samples:
            t_end = rec.samples[-1][0]
            pts = [p for p in rec.samples if p[0] >= t_end - WINDOW_S]
            pts = pts[::max(1, len(pts) // 1500)]
            t = [p[0] for p in pts]
            cmd = [p[1] for p in pts]
            enc = [p[2] for p in pts]
            err = [e - c for c, e in zip(cmd, enc)]
            xlim = (max(t_end - WINDOW_S, 0.0), max(t_end, WINDOW_S))
            lo, hi = min(cmd + enc), max(cmd + enc)
            pad = max((hi - lo) * 0.1, 0.5)
            angle.draw(t, [(cmd, C_CMD), (enc, C_ENC)], xlim, (lo - pad, hi + pad))
            m = max(max(abs(x) for x in err), 0.05) * 1.2
            error.draw(t, [(err, C_ERR)], xlim, (-m, m))
        root.after(500, update)

    root.after(200, update)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
