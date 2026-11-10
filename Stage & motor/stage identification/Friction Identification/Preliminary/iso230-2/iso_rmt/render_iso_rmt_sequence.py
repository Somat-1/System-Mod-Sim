#!/usr/bin/env python3
"""Render the RMT ISO 230-2 sequence: full timeline and one-run detail."""
from __future__ import annotations

import csv, json, math
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "rendered"
OVERVIEW = OUT / "iso_rmt_sequence.png"
COL = {"HOME": "#6b7280", "INIT": "#6b7280", "DATUM": "#7c5a9b",
       "ISO_REF": "#a56a2a", "ISO": "#176d8c", "END_REF": "#7c5a9b",
       "PARK": "#4c5560", "record": "#6d4c91"}

def load():
    config = json.loads((ROOT / "iso_rmt_config.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "sequence/iso_rmt_plan_summary.json").read_text(encoding="utf-8"))
    with (ROOT / "sequence/iso_rmt_measurement_plan.csv").open(newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh))
    return config, summary, plan

def move_time(distance, c):
    d, v, a = abs(distance), c["velocity_mm_s"], c["acceleration_mm_s2"]
    return 2 * math.sqrt(d / a) if d <= v * v / a else 2 * v / a + (d - v * v / a) / v

def schedule(plan, c):
    """Commanded position against time, from the plan rows."""
    now, pos, state = 0.0, c["x_ref_mm"], "INIT"
    segments, records, runs = [], [], {}
    for row in plan:
        if row["action"] == "state_enter":
            state = row["state"]
            if state == "ISO_REF":
                runs.setdefault(int(row["run"]), [now, now])
            continue
        if row["target_mm"]:
            target = float(row["target_mm"])
            end = now + move_time(target - pos, c)
            segments.append((now, end, pos, target, state))
            now, pos = end, target
        settle, record = float(row["settle_s"] or 0), float(row["record_s"] or 0)
        if settle or record:
            segments.append((now, now + settle + record, pos, pos, state))
            records.append(((now + settle) / 60, pos, state))
            now += settle + record
        if row["run"] and int(row["run"]) in runs:
            runs[int(row["run"])][1] = now
    return segments, records, runs, now

def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

def draw(ax, segments, records, t0=0.0, t1=None):
    for a, b, x0, x1, state in segments:
        if t1 is not None and (b / 60 < t0 or a / 60 > t1):
            continue
        ax.plot([a / 60, b / 60], [x0, x1], color=COL[state], lw=1.0,
                solid_capstyle="butt")
    for t, x, state in records:
        if t1 is None or t0 <= t <= t1:
            ax.plot(t, x, "o", ms=2.4 if t1 is None else 3.6, color=COL["record"],
                    mec="white", mew=0.4, zorder=3)
    ax.grid(color="#e3e6e8", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)

def render(c, summary, segments, records, runs, total):
    fig = plt.figure(figsize=(7.45, 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.6], left=0.105, right=0.97,
                          top=0.815, bottom=0.08, hspace=0.42)
    fig.text(0.105, 0.975, "ISO 230-2 v4 — RMT sequence", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.105, 0.928,
             f"{summary['targets']} v3 targets (seed {c['seed']}) · {c['n_runs']} bidirectional runs · "
             f"{summary['analysis_stops']} stops + {summary['reference_stops']} x_ref visits · "
             f"~{total / 60:.1f} min after homing",
             fontsize=8.5, color="#4a5157", va="top")
    fig.text(0.105, 0.902,
             f"RMT motion, silent while moving · {c['velocity_mm_s']:g} mm/s at "
             f"{c['acceleration_mm_s2']:g} mm/s² · settle {c['settle_s']:g} s + record {c['record_s']:g} s · "
             f"{c['run_current_ma']} mA, hold = run",
             fontsize=8.5, color="#4a5157", va="top")

    ax = fig.add_subplot(gs[0])
    draw(ax, segments, records)
    for run, (a, b) in runs.items():
        ax.text((a + b) / 120, c["soft_max_mm"] + 3, f"run {run}", ha="center",
                fontsize=7.5, color="#41474c")
    ax.axhline(c["x_ref_mm"], color="#a56a2a", lw=0.6, ls=(0, (4, 3)))
    ax.text(0.2, c["x_ref_mm"] + 2, f"x_ref {c['x_ref_mm']:g} mm", fontsize=7, color="#7a4f1f")
    ax.set_xlim(0, total / 60)
    ax.set_ylim(0, c["soft_max_mm"] + 8)
    ax.set_xlabel("time after centring [min]")
    ax.set_ylabel("commanded position [mm]")
    ax.set_title("Full sequence: DATUM ×3, 5 × (x_ref, up, down), END_REF ×3, PARK",
                 loc="left", fontsize=9.5)
    handles = [Line2D([], [], color=COL[k], lw=1.6, label=label)
               for k, label in (("DATUM", "DATUM / END REF"), ("ISO_REF", "ISO REF"),
                                ("ISO", "ISO"))]
    handles.append(Line2D([], [], marker="o", ls="", color=COL["record"], ms=4,
                          label="record window"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.095, 0.892), ncol=4,
               frameon=False, fontsize=7.5, handlelength=1.8, columnspacing=1.6)

    a, b = runs[1]
    ax2 = fig.add_subplot(gs[1])
    draw(ax2, segments, records, a / 60 - 0.1, b / 60 + 0.1)
    ax2.set_xlim(a / 60 - 0.1, b / 60 + 0.1)
    ax2.set_ylim(0, c["soft_max_mm"] + 4)
    ax2.set_xlabel("time after centring [min]")
    ax2.set_ylabel("position [mm]")
    ax2.set_title(f"Run 1 ({(b - a) / 60:.1f} min): x_ref from {c['ref_preapproach_mm']:g} mm, "
                  f"overrun {c['overrun_mm']:g} mm at both ends, {summary['targets']} stops each way",
                  loc="left", fontsize=9.5)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OVERVIEW, dpi=220)
    plt.close(fig)

def main() -> int:
    c, summary, plan = load()
    style()
    segments, records, runs, total = schedule(plan, c)
    render(c, summary, segments, records, runs, total)
    print(f"wrote {OVERVIEW.relative_to(ROOT)} ({total / 60:.1f} min)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
