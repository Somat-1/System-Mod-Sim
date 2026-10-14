#!/usr/bin/env python3
"""Live window for a chirp run: progress, and the response building up.

Follows chirp_<profile>_windows.csv as run_chirp.py writes it (it never touches
the serial port, so it can be closed and reopened freely) and shows:

  * progress through the profile, the current frequency, and any fault;
  * generator check - issued amplitude / target and issued frequency / target,
    against frequency. Both should sit on 1.0: that is the proof the step
    generator kept up;
  * response - encoder amplitude / issued amplitude, against frequency
    (1.0 = the rotor follows the command), with the phase underneath.

    /opt/homebrew/bin/python3 scripts/monitor_chirp.py <folder or windows.csv>
"""

from __future__ import annotations

import csv
import json
import sys
import tkinter as tk
from pathlib import Path
from tkinter import ttk

ROOT = Path(__file__).resolve().parent.parent
BG, GRID, INK = "#ffffff", "#e4e4e4", "#333333"
C_UP, C_DOWN, C_REF = "#2f6db5", "#d4762c", "#999999"


class Plot:
    def __init__(self, parent, height, ylabel):
        self.c = tk.Canvas(parent, height=height, bg=BG, highlightthickness=0)
        self.c.pack(fill="both", expand=True, padx=8, pady=3)
        self.ylabel = ylabel

    def draw(self, series, xlim, ylim, ref=None):
        c = self.c
        c.delete("all")
        w, h = max(c.winfo_width(), 200), max(c.winfo_height(), 80)
        l, r, t, b = 64, 12, 10, 22
        (x0, x1), (y0, y1) = xlim, ylim
        px = lambda x: l + (x - x0) / (x1 - x0) * (w - l - r)
        py = lambda y: h - b - (min(max(y, y0), y1) - y0) / (y1 - y0) * (h - t - b)
        for i in range(5):
            y = y0 + (y1 - y0) * i / 4
            c.create_line(l, py(y), w - r, py(y), fill=GRID)
            c.create_text(l - 5, py(y), text=f"{y:.3g}", anchor="e", fill=INK,
                          font=("Helvetica", 10))
        for x in range(0, int(x1) + 1, 100):
            c.create_text(px(x), h - b + 11, text=f"{x}", fill=INK, font=("Helvetica", 9))
        c.create_text(8, t, text=self.ylabel, anchor="nw", fill=INK,
                      font=("Helvetica", 10, "bold"))
        if ref is not None:
            c.create_line(l, py(ref), w - r, py(ref), fill=C_REF, dash=(4, 3))
        for pts, colour in series:
            for x, y in pts:
                X, Y = px(x), py(y)
                c.create_oval(X - 1.5, Y - 1.5, X + 1.5, Y + 1.5, fill=colour, outline="")


def main() -> int:
    arg = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if arg is None:
        runs = sorted((ROOT / "recordings").glob("*/chirp_*_windows.csv"),
                      key=lambda p: p.stat().st_mtime)
        if not runs:
            print("no chirp recording found; pass the windows.csv or its folder")
            return 2
        arg = runs[-1]
    path = arg if arg.suffix == ".csv" else next(arg.glob("chirp_*_windows.csv"))
    log = path.with_name(path.name.replace("_windows.csv", "_run.log"))
    profile = "full" if "_full_" in path.name else "high"
    summary = json.loads((ROOT / "sequence" / "chirp_summary.json").read_text())
    total = summary["profiles"][profile]["duration_s"] + 7    # + countdown

    root = tk.Tk()
    root.title(f"Chirp {profile.upper()} - {path.parent.name}")
    root.geometry("1000x780")
    root.configure(bg=BG)
    head = tk.Label(root, anchor="w", bg=BG, fg=INK, font=("Helvetica", 14, "bold"))
    head.pack(fill="x", padx=10, pady=(10, 2))
    bar = ttk.Progressbar(root, maximum=1000)
    bar.pack(fill="x", padx=10, pady=3)
    info = tk.Label(root, anchor="w", bg=BG, fg=INK, font=("Menlo", 12))
    info.pack(fill="x", padx=10)
    tk.Label(root, text="blue: up sweep   orange: down sweep   x axis: frequency [Hz]",
             anchor="w", bg=BG, fg=INK, font=("Helvetica", 11)).pack(fill="x", padx=10)
    gen = Plot(root, 170, "generator: issued / target (amplitude)")
    resp = Plot(root, 230, "response: encoder / issued amplitude")
    phase = Plot(root, 170, "phase, encoder vs issued [deg]")

    state = {"offset": 0, "rows": [], "status": "waiting", "partial": ""}

    def poll():
        try:
            with path.open(encoding="utf-8") as f:
                f.seek(state["offset"])
                chunk = f.read()
                state["offset"] = f.tell()
        except FileNotFoundError:
            chunk = ""
        text = state["partial"] + chunk
        lines = text.split("\n")
        state["partial"] = lines.pop()
        for ln in lines:
            if not ln or ln.startswith("t_s"):
                continue
            try:
                state["rows"].append([float(x) for x in ln.split(",")])
            except ValueError:
                pass
        try:
            for ln in log.read_text(encoding="utf-8").splitlines()[-40:]:
                if ln.startswith(("# RUN_", "# ERR", "# ABORT", "# COUNTDOWN")):
                    state["status"] = ln[2:]
                elif ln.startswith("# SEG") and not state["status"].startswith(("RUN_", "ERR")):
                    state["status"] = "running  " + ln.split(",")[1]
        except FileNotFoundError:
            pass

    def update():
        poll()
        rows = state["rows"]
        t = rows[-1][0] if rows else 0.0
        bar["value"] = 1000 if state["status"].startswith("RUN_COMPLETE") else min(t / total, 1) * 1000
        bad = state["status"].startswith(("ERR", "RUN_FAILED", "RUN_ABORTED"))
        head.config(text=state["status"], fg="#c0392b" if bad else INK)
        if rows:
            r = rows[-1]
            gen_ok = [x[5] / x[4] for x in rows if x[4] > 0]
            info.config(text=f"t {t:6.1f} s   f {r[1]:7.1f} Hz   target {r[4]:7.2f}   "
                             f"issued {r[5]:7.2f}   encoder {r[6]:7.2f} us   "
                             f"worst issued/target {min(gen_ok):.3f}..{max(gen_ok):.3f}")
            # up vs down: frequency rising or falling relative to the previous window
            up, down, gup, gdn, pup, pdn = [], [], [], [], [], []
            for prev, x in zip(rows, rows[1:]):
                if x[4] <= 0 or x[5] <= 0:
                    continue
                g = x[6] / x[5]
                tgt = (up, gup, pup) if x[1] >= prev[1] else (down, gdn, pdn)
                tgt[0].append((x[1], x[5] / x[4]))
                tgt[1].append((x[1], g))
                tgt[2].append((x[1], x[7]))
            xlim = (0, 1000)
            gen.draw([(up, C_UP), (down, C_DOWN)], xlim, (0.8, 1.2), ref=1.0)
            gmax = max([g for _, g in gup + gdn] + [1.2])
            resp.draw([(gup, C_UP), (gdn, C_DOWN)], xlim, (0, min(gmax * 1.1, 5)), ref=1.0)
            phase.draw([(pup, C_UP), (pdn, C_DOWN)], xlim, (-180, 180), ref=0)
        root.after(500, update)

    root.after(200, update)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
