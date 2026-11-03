#!/usr/bin/env python3
"""Render publication-style E2 sequence and target-design figures."""
from __future__ import annotations

import csv, json, math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "rendered"
OVERVIEW = OUT / "e2_global_map_sequence.png"
DETAIL = OUT / "e2_target_design_and_rmt.png"
COL = {"INIT": "#6b7280", "DATUM": "#7c5a9b", "UP": "#176d8c",
       "TURN_REF": "#a56a2a", "DOWN": "#b4473f", "END_REF": "#7c5a9b",
       "PARK": "#4c5560", "record": "#6d4c91", "settle": "#c49a54",
       "small": "#b4473f", "grid": "#dde1e4"}
# The fine zoom sits on a single-grid-point gap flanked by typical ones, so the
# smallest and the median step appear at the same scale.
FINE_FIRST, FINE_LAST = 167, 175

def load():
    config = json.loads((ROOT / "e2_config.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "sequence/e2_plan_summary.json").read_text(encoding="utf-8"))
    with (ROOT / "sequence/e2_measurement_plan.csv").open(newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh))
    with (ROOT / "sequence/e2_targets_seeded_grid64.csv").open(newline="", encoding="utf-8") as fh:
        targets = list(csv.DictReader(fh))
    iso_path = (ROOT / config["iso_target_source"]).resolve()
    with iso_path.open(newline="", encoding="utf-8") as fh:
        iso = np.array([float(row["cmd_mm"]) for row in csv.DictReader(fh)])
    return config, summary, plan, targets, iso

def move_time(distance, c):
    d, v, a = abs(distance), c["velocity_mm_s"], c["acceleration_mm_s2"]
    return 2 * math.sqrt(d / a) if d <= v * v / a else 2 * v / a + (d - v * v / a) / v

def schedule(plan, c):
    now, pos, state = 0.0, c["x_ref_mm"], "INIT"
    segments, blocks, records, stops = [], [], [], []
    block_start = 0.0
    for row in plan:
        if row["action"] == "state_enter":
            if now > block_start:
                blocks.append((block_start, now, state))
            state, block_start = row["state"], now
            continue
        if row["target_mm"]:
            target = float(row["target_mm"])
            end = now + move_time(target - pos, c)
            segments.append((now, end, pos, target, state, "move"))
            now, pos = end, target
        settle, record = float(row["settle_s"] or 0), float(row["record_s"] or 0)
        hold_start = now
        if settle:
            segments.append((now, now + settle, pos, pos, state, "settle")); now += settle
        if record:
            start = now
            segments.append((now, now + record, pos, pos, state, "record")); now += record
            if state in ("UP", "DOWN"):
                records.append((int(row["index"]), state, (start + now) / 2))
        if row["record"] == "True":
            stops.append((hold_start, now, pos, state,
                          int(row["index"]) if row["index"] else 0, row["direction"]))
    if now > block_start:
        blocks.append((block_start, now, state))
    return segments, blocks, records, stops, now


def staircase(stops, lo_index, hi_index, direction):
    """The commanded position against session time over a run of stops."""
    chosen = [s for s in stops if s[5] == direction and lo_index <= s[4] <= hi_index]
    time, position = [], []
    for start, end, pos, _state, _index, _dir in chosen:
        time += [start, end]
        position += [pos, pos]
    return np.array(time), np.array(position), chosen

def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

def render_overview(c, summary, segments, blocks, stops, total):
    total_min = total / 60
    fig = plt.figure(figsize=(7.45, 8.7))
    gs = fig.add_gridspec(4, 6, height_ratios=[0.62, 2.45, 1.75, 0.80],
                          left=0.095, right=0.975, top=0.895, bottom=0.068,
                          hspace=0.52, wspace=1.05)
    fig.text(0.095, 0.978, "Experiment E2 — global geometric map", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.095, 0.942,
             "400 stratified targets · mirrored directional passes · seven drift references · buffered RMT motion",
             fontsize=8.5, color="#566069", va="top")

    ax0 = fig.add_subplot(gs[0, :])
    for start, end, state in blocks:
        width = (end - start) / 60
        ax0.barh(0, width, left=start / 60, height=0.38, color=COL.get(state, "#777"),
                 edgecolor="white", linewidth=0.6)
        label = state.replace("TURN_REF", "TURN").replace("END_REF", "E")
        if width > 0.40:
            ax0.text((start + end) / 120, 0, label, ha="center", va="center",
                     fontsize=6.3, color="white", weight="bold")
    ax0.set(xlim=(0, total_min), ylim=(-0.35, 0.35), yticks=[])
    ax0.set_title("(a) Time allocation by controller state", loc="left", pad=6)
    ax0.tick_params(axis="x", bottom=False, labelbottom=False)
    ax0.spines[:].set_visible(False)

    ax = fig.add_subplot(gs[1, :], sharex=ax0)
    for state in ("DATUM", "UP", "TURN_REF", "DOWN", "END_REF", "PARK"):
        lines = [[(a / 60, x0), (b / 60, x1)] for a, b, x0, x1, st, _ in segments if st == state]
        if lines:
            ax.add_collection(LineCollection(lines, colors=COL[state], linewidths=1.05))
    ax.axhline(c["x_ref_mm"], color=COL["DATUM"], lw=0.7, ls=(0, (3, 2)))
    ax.text(total_min * 0.015, c["x_ref_mm"] + 2.1, r"$x_{ref}=61.35$ mm",
            color=COL["DATUM"], fontsize=7)
    ax.axhspan(c["target_min_mm"], c["target_max_mm"], color="#176d8c", alpha=0.035)
    ax.set(xlim=(0, total_min), ylim=(0, 123), ylabel="commanded position [mm]")
    ax.set_title("(b) Complete commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color="#dde1e4", lw=0.55); ax.spines[["top", "right"]].set_visible(False)
    coarse_t, coarse_x, _ = staircase(stops, FINE_FIRST - 25, FINE_LAST + 25, "up")
    ax.axvspan(coarse_t[0] / 60, coarse_t[-1] / 60, color=COL["UP"], alpha=0.14, lw=0)
    ax.annotate("(c)", xy=(coarse_t[0] / 60, coarse_x.max() + 8), fontsize=7.4,
                color=COL["UP"], weight="bold", ha="center")

    step_mm = np.diff(np.array([float(row) for row in
                                [s[2] for s in stops if s[3] == "UP"]]))
    box = (f"seed  {summary['seed_used']}\n"
           f"targets  {summary['targets']} × 2 directions\n"
           f"reference records  {summary['reference_stops']}\n"
           f"step  {step_mm.min():.2f}–{step_mm.max():.2f} mm "
           f"(mean {step_mm.mean():.2f})\n"
           f"nominal duration  {total_min:.1f} min\n"
           f"hard limit  {c['time_limit_min']:.0f} min")
    # Top left: the up pass climbs away from this corner, so it stays clear.
    ax.text(0.015, 0.965, box, transform=ax.transAxes, ha="left", va="top",
            fontsize=7.2, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.48", fc="white", ec="#aeb4b9", lw=0.7))

    # (c) a coarse zoom: enough stops to show the staircase as a recording sees it.
    ax = fig.add_subplot(gs[2, 0:3])
    ax.plot(coarse_t - coarse_t[0], coarse_x, color=COL["UP"], lw=1.0,
            solid_joinstyle="miter")
    fine_t, fine_x, fine_stops = staircase(stops, FINE_FIRST, FINE_LAST, "up")
    ax.axvspan(fine_t[0] - coarse_t[0], fine_t[-1] - coarse_t[0],
               color=COL["small"], alpha=0.13, lw=0)
    ax.annotate("(d)", xy=((fine_t[0] + fine_t[-1]) / 2 - coarse_t[0], fine_x.max() + 1.6),
                fontsize=7.4, color=COL["small"], weight="bold", ha="center")
    ax.set(xlabel="time within the up pass [s]", ylabel="position [mm]")
    ax.set_title(f"(c) Up pass, targets {FINE_FIRST - 25}–{FINE_LAST + 25}",
                 loc="left", pad=6)
    ax.grid(color=COL["grid"], lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.035, 0.955,
            f"each tread = {c['settle_s']:.1f} s settle + {c['record_s']:.1f} s record",
            transform=ax.transAxes, va="top", fontsize=6.8, color="#566069")

    # (d) every riser labelled, the single-grid-point one included: at 40 µm
    # against a 2.4 mm span it is nearly invisible, which is the point.
    ax = fig.add_subplot(gs[2, 3:6])
    base = fine_x[0]
    ax.plot(fine_t - fine_t[0], (fine_x - base) * 1000, color=COL["UP"], lw=1.25,
            solid_joinstyle="miter")
    for j in range(len(fine_stops) - 1):
        x_lo = (fine_stops[j][2] - base) * 1000
        x_hi = (fine_stops[j + 1][2] - base) * 1000
        grid_points = int(round((x_hi - x_lo) / 40.0))
        when = fine_stops[j + 1][0] - fine_t[0]
        smallest = grid_points == 1
        # A 1 gp riser is thinner than the line, so its label clears the trace.
        ax.text(when + (-0.4 if smallest else 0.5),
                (x_lo + x_hi) / 2 + (190.0 if smallest else 0.0),
                f"{grid_points * 0.04:.2f} mm · {grid_points} gp",
                fontsize=6.1, va="bottom" if smallest else "center",
                ha="right" if smallest else "left",
                color=COL["small"] if smallest else "#4b5257",
                weight="bold" if smallest else "normal",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8)
                if smallest else None)
    ax.set(xlabel="time within the up pass [s]", ylabel=f"position − {base:.2f} mm [µm]")
    ax.set_title(f"(d) Nine stops, targets {FINE_FIRST}–{FINE_LAST}", loc="left", pad=6)
    span = fine_t[-1] - fine_t[0]
    ax.set_xlim(-span * 0.03, span * 1.34)
    ax.margins(y=0.12)
    ax.grid(color=COL["grid"], lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    single = int((np.round(step_mm / 0.04).astype(int) == 1).sum())
    ax.text(0.035, 0.965,
            "gp = grid point = 64 µsteps = 0.040 mm\n"
            "     = one motor electrical cycle\n"
            f"the 1 gp step is one of {single} in the 400-target list",
            transform=ax.transAxes, va="top", fontsize=6.4, color="#566069",
            linespacing=1.45)

    ax2 = fig.add_subplot(gs[3, 0:6], sharex=ax0)
    dwell = c["settle_s"] + c["record_s"]
    components = [("motion + transits", total - 807 * dwell, "#6d7d88"),
                  ("settle windows", 807 * c["settle_s"], COL["settle"]),
                  ("record windows", 807 * c["record_s"], COL["record"])]
    left = 0.0
    for label, seconds, color in components:
        width = seconds / 60
        ax2.barh(0, width, left=left, height=0.48, color=color, edgecolor="white")
        if width > 2:
            ax2.text(left + width / 2, 0, f"{label}\n{width:.1f} min", ha="center",
                     va="center", fontsize=6.7, color="white", weight="bold")
        left += width
    ax2.set(yticks=[], xlabel="elapsed recording time [min]", ylim=(-0.42, 0.42))
    ax2.set_title("(e) Recording-time budget", loc="left", pad=6)
    ax2.spines[["left", "right", "top"]].set_visible(False)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OVERVIEW, dpi=320); plt.close(fig)

def render_detail(c, summary, targets, iso, records, total):
    x = np.array([float(row["target_mm"]) for row in targets])
    phase = np.array(summary["phase_counts"])
    fig, axes = plt.subplots(2, 2, figsize=(8.27, 7.15))
    fig.subplots_adjust(left=0.09, right=0.965, top=0.88, bottom=0.10, hspace=0.42, wspace=0.29)
    fig.text(0.09, 0.975, "E2 sampling design and command fidelity", fontsize=16,
             weight="bold", color="#20262b", va="top")
    fig.text(0.09, 0.94,
             "Restricted stratified random design on the fixed electrical phase grid",
             fontsize=8.6, color="#566069", va="top")

    ax = axes[0, 0]
    ax.vlines(x, 0.0, 0.70, color=COL["UP"], lw=0.55, alpha=0.82)
    ax.scatter(iso, np.full_like(iso, 0.90), marker="x", s=18, color="#b4473f",
               linewidths=0.8, label="E0 validation targets")
    ax.hlines(0.90, 2, 120, color="#c4c9cd", lw=0.5)
    ax.text(3, 0.60, "400 E2 targets — one per 0.295 mm stratum", color=COL["UP"], fontsize=7.2,
            bbox=dict(fc="white", ec="none", alpha=0.82, pad=1.5))
    ax.set(xlim=(0, 122), ylim=(-0.03, 1.08), yticks=[], xlabel="commanded position [mm]")
    ax.set_title("(a) Full-travel coverage and independent E0 set", loc="left", pad=7)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    ax.spines[["left", "right", "top"]].set_visible(False)

    ax = axes[0, 1]
    bins = np.arange(50)
    ax.bar(bins, phase, width=0.82, color="#2f7d74", edgecolor="white", linewidth=0.25)
    ax.axhline(8, color="#545b61", lw=0.8, ls=(0, (3, 2)), label="uniform mean = 8")
    ax.axhline(c["phase_min_count"], color="#b4473f", lw=0.7, ls=":", label="minimum = 5")
    ax.set(xlim=(-1, 50), ylim=(0, max(phase) + 1.5), xlabel="screw-phase bin (64 µsteps)",
           ylabel="target count")
    ax.set_title("(b) Phase balance over one 2 mm screw revolution", loc="left", pad=7)
    ax.legend(frameon=False, fontsize=7); ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e1e3e5", lw=0.5)

    ax = axes[1, 0]
    up = {i: t / 60 for i, state, t in records if state == "UP"}
    down = {i: t / 60 for i, state, t in records if state == "DOWN"}
    indices = sorted(set(up) & set(down))
    xx = np.array([x[i - 1] for i in indices])
    tu = np.array([up[i] for i in indices]); td = np.array([down[i] for i in indices])
    ax.plot(xx, tu, color=COL["UP"], lw=1.1, label="ascending visit")
    ax.plot(xx, td, color=COL["DOWN"], lw=1.1, label="descending visit")
    ax.plot(xx, 0.5 * (tu + td), color="#30363b", lw=0.9, ls="--", label="paired mean time")
    ax.fill_between(xx, tu, td, color="#788691", alpha=0.08)
    ax.set(xlabel="target position [mm]", ylabel="record-window midpoint [min]")
    ax.set_title("(c) Mirrored timing suppresses linear drift bias", loc="left", pad=7)
    ax.legend(frameon=False, fontsize=7); ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e1e3e5", lw=0.5)

    ax = axes[1, 1]
    pulses = int(np.median(np.diff(np.array([int(row["target_usteps"]) for row in targets]))))
    distance = pulses / c["ustep_per_mm"]
    v, a = c["velocity_mm_s"], c["acceleration_mm_s2"]
    if distance <= v * v / a:
        t_acc, t_cruise, peak = math.sqrt(distance / a), 0.0, math.sqrt(distance * a)
    else:
        t_acc, peak = v / a, v
        t_cruise = (distance - peak * peak / a) / peak
    t_end = 2 * t_acc + t_cruise
    time = np.linspace(0, t_end, 500)
    vel = np.where(time <= t_acc, a * time,
          np.where(time <= t_acc + t_cruise, peak, a * (t_end - time)))
    pos = np.zeros_like(time)
    for j in range(1, len(time)):
        pos[j] = pos[j - 1] + 0.5 * (vel[j - 1] + vel[j]) * (time[j] - time[j - 1])
    ax.fill_between(time * 1000, vel, color=COL["UP"], alpha=0.15)
    ax.plot(time * 1000, vel, color=COL["UP"], lw=1.6)
    ax2 = ax.twinx(); ax2.plot(time * 1000, pos, color="#66517b", lw=1.0, ls=":")
    ax2.set_ylabel("displacement [mm]", color="#66517b"); ax2.set_ylim(0, distance * 1.12)
    chunks = [min(c["rmt_chunk_pulses"], pulses - i)
              for i in range(0, pulses, c["rmt_chunk_pulses"])]
    ax.text(0.03, 0.96,
            f"median move  {distance:.3f} mm = {pulses} STEP pulses\n"
            f"RMT transactions  {' + '.join(map(str, chunks))} symbols\n"
            f"then SETTLE {c['settle_s']:.2f} s + RECORD {c['record_s']:.2f} s",
            transform=ax.transAxes, va="top", fontsize=7.1, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.38", fc="white", ec="#b8bdc1", lw=0.6))
    ax.set(xlabel="time from move command [ms]", ylabel="commanded velocity [mm/s]")
    ax.set_title("(d) Representative buffered RMT move", loc="left", pad=7)
    ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False)
    ax.grid(color="#e1e3e5", lw=0.5)

    fig.text(0.09, 0.035,
             "Safety adaptation: nominal −1/123 mm approach points lie outside the proven E1 soft envelope; "
             "the executable plan uses 2.000/120.695 mm and records this deviation in every session header.",
             fontsize=7.1, color="#5b6268")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(DETAIL, dpi=320); plt.close(fig)

def main() -> int:
    c, summary, plan, targets, iso = load()
    style()
    segments, blocks, records, stops, total = schedule(plan, c)
    render_overview(c, summary, segments, blocks, stops, total)
    render_detail(c, summary, targets, iso, records, total)
    print(f"wrote {OVERVIEW.relative_to(ROOT)}")
    print(f"wrote {DETAIL.relative_to(ROOT)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
