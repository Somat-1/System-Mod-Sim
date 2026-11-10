#!/usr/bin/env python3
"""One figure that explains what the F0b session does.

Same form as the E1-E5 sequence renders: block allocation, the complete
commanded trajectory, and where the session goes - with one panel added that
only F0b needs.

F0b is a scoping session, so what matters about it is not where the stage
goes but how large the moves are. They span five and a half decades, from a
single 1/256 microstep to a 3 mm stroke, which is why the trajectory panel
carries a symmetric-log axis and why panel (c) exists at all: a linear axis
would draw the whole session as one flat line at x0.

Everything plotted is read from `sequence/`, which `f0b_plan.py` writes from
the procedure. Nothing is measured yet.

    python render_f0b_sequence.py
"""
from __future__ import annotations

import csv
import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

HERE = Path(__file__).resolve().parent
SEQ = HERE / "sequence"
OUT = HERE / "rendered"
FIG = OUT / "f0b_sequence_overview.png"

C_INK, C_MUTED, C_GRID = "#20262b", "#566069", "#dde1e4"
COL = {"R1": "#73578f", "R2": "#73578f", "R3": "#73578f",
       "N": "#7b858c", "L": "#11536b", "S": "#c08a3e",
       "D": "#4f7d6a", "V": "#b4473f"}
BUDGET = {"move": "#697b86", "ramp": "#11536b", "step": "#4f7d6a",
          "dwell": "#c59b50"}
CHAR_INCH = 0.0455


def caption(fig, left, text, bottom=0.012, size=6.9) -> None:
    usable = fig.get_figwidth() * (1.0 - left) - 0.12
    width = max(40, int(usable / CHAR_INCH))
    lines = textwrap.wrap(" ".join(text.split()), width=width)
    fig.text(left, bottom, "\n".join(lines), fontsize=size, color=C_MUTED,
             va="bottom", linespacing=1.5)


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def load():
    cfg = json.loads((HERE / "f0b_config.json").read_text(encoding="utf-8"))
    summary = json.loads((SEQ / "f0b_plan_summary.json")
                         .read_text(encoding="utf-8"))
    with (SEQ / "f0b_schedule.csv").open(newline="", encoding="utf-8") as fh:
        sched = list(csv.DictReader(fh))
    return cfg, summary, sched


def blocks_of(sched) -> list[tuple[str, float, float]]:
    """Contiguous runs of one block name, as (name, t_start, t_end) minutes."""
    out = []
    for row in sched:
        name = row["block"]
        t0, t1 = float(row["t_start_s"]) / 60.0, float(row["t_end_s"]) / 60.0
        if out and out[-1][0] == name:
            out[-1][2] = t1
        else:
            out.append([name, t0, t1])
    return [(n, a, b) for n, a, b in out]


# ------------------------------------------------------------------ panels

def panel_allocation(ax, blocks, total_min) -> None:
    for name, t0, t1 in blocks:
        width = t1 - t0
        ax.barh(0, width, left=t0, height=0.42, color=COL[name],
                edgecolor="white", linewidth=0.6)
        if width > 1.4:
            ax.text((t0 + t1) / 2.0, 0.055, name, ha="center", va="center",
                    fontsize=7.2, color="white", weight="bold")
            ax.text((t0 + t1) / 2.0, -0.10, f"{width:.1f}", ha="center",
                    va="center", fontsize=5.8, color="white")
    ax.set(xlim=(0, total_min), ylim=(-0.36, 0.36), yticks=[])
    ax.set_title("(a) Block allocation, with minutes", loc="left", pad=6)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines[:].set_visible(False)


def panel_trajectory(ax, cfg, sched, total_min) -> None:
    """The whole session, as commanded displacement from x0.

    Symmetric log, because this one axis has to hold a 45 mm transit to the
    reference positions and a 39 nm microstep.
    """
    x0 = cfg["session"]["x0_mm"]
    grouped: dict[str, list] = {}
    for row in sched:
        seg = [(float(row["t_start_s"]) / 60.0,
                (float(row["x_start_mm"]) - x0) * 1000.0),
               (float(row["t_end_s"]) / 60.0,
                (float(row["x_end_mm"]) - x0) * 1000.0)]
        grouped.setdefault(COL[row["block"]], []).append(seg)
    for colour, lines in grouped.items():
        ax.add_collection(LineCollection(lines, colors=colour, linewidths=1.0))

    ustep = cfg["drive"]["ustep_um"]
    ax.set_yscale("symlog", linthresh=ustep, linscale=0.45)
    ax.set(xlim=(0, total_min), ylim=(-9e4, 9e4),
           ylabel="commanded position − x$_0$ [µm]")
    ticks = [-1e4, -1e2, -1, 0, 1, 1e2, 1e4]
    ax.set_yticks(ticks)
    ax.set_yticklabels(["−10 000", "−100", "−1", "0",
                        "1", "100", "10 000"])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    for y, text, side in ((45000, "reference pair, 105 mm", "right"),
                          (-45000, "reference pair, 15 mm", "right"),
                          (20000, "N2, 20 mm above x$_0$", "left")):
        ax.axhline(y, color="#9aa2a8", lw=0.5, ls=(0, (3, 2)))
        ax.text(total_min * (0.008 if side == "left" else 0.992), y * 1.10,
                text, fontsize=6.0, color="#6a737a", ha=side, va="bottom")
    ax.set_title("(b) Complete commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color=C_GRID, lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], color=COL[k], lw=1.8, label=lab)
               for k, lab in (("R1", "reference pair"), ("N", "N  noise floor"),
                              ("L", "L  reversal ladder"),
                              ("S", "S  slow ramps"), ("D", "D  dwell probe"),
                              ("V", "V  speed probe"))]
    ax.legend(handles=handles, frameon=False, fontsize=6.8, ncol=6,
              loc="upper center", bbox_to_anchor=(0.5, -0.095),
              handlelength=1.5, columnspacing=1.2)


def panel_amplitudes(ax, cfg, summary) -> None:
    """What each block interrogates, on one logarithmic amplitude axis."""
    ustep = cfg["drive"]["ustep_um"]
    strokes = [s * ustep for s in cfg["L"]["strokes_usteps"]]
    floor = cfg["reference_values"]["interferometer_floor_um"]
    beta0 = cfg["reference_values"]["beta0_um"]
    biggest = f"{cfg['L']['strokes_usteps'][-1]:,}".replace(",", " ")

    rows = [
        ("N", [floor], floor, floor,
         "noise floor and drift at rest; d$_{min}$ follows from it"),
        ("L", strokes, strokes[0], strokes[-1],
         f"{len(strokes)} strokes, 1 to {biggest} µsteps, "
         "four half-cycles each, twice"),
        ("S", [ustep, cfg["S"]["ramp_um"]], ustep, cfg["S"]["ramp_um"],
         "20 µm at 0.25 µm/s, resolved microstep by microstep"),
        ("D", [cfg["D"]["probe_um"]], cfg["D"]["probe_um"],
         cfg["D"]["probe_um"], "5 µm after 1 s and after 60 s at rest"),
        ("V", [1000.0], 1000.0, 1000.0,
         "1 mm traversed at 0.01 and at 1 mm/s"),
    ]
    # A note starting at the right-hand end of the axis has nowhere to run,
    # so the rows near the top of the range hang their note to the left.
    side = {"V": "right"}
    for i, (key, points, lo, hi, note) in enumerate(rows):
        y = len(rows) - 1 - i
        if hi > lo:
            ax.plot([lo, hi], [y, y], color=COL[key], lw=2.6, alpha=0.30,
                    solid_capstyle="round", zorder=2)
        ax.plot(points, [y] * len(points), marker="|", ms=7, mew=1.3,
                ls="none", color=COL[key], zorder=3)
        # The note sits above its own row and starts where the row starts,
        # so it crosses neither another row nor the reference labels below.
        ha = side.get(key, "left")
        ax.text(lo * (1.16 if ha == "right" else 0.86), y + 0.30, note,
                fontsize=6.3, color=C_MUTED, va="center", ha=ha)

    marks = (
        (floor, "interferometer floor\n≈2 nm", "#6a737a", "left"),
        (ustep, "1 µstep\n39 nm", "#6a737a", "center"),
        (beta0, f"β$_0$ = {beta0:.2f} µm\n(expected ℓ$_{{sat}}$)",
         "#8d3934", "center"),
        (cfg["drive"]["electrical_cycle_um"], "electrical cycle\n40 µm",
         "#6a737a", "center"),
    )
    for x, label, colour, ha in marks:
        ax.axvline(x, color=colour, lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.text(x, -1.26, label, fontsize=6.0, color=colour, ha=ha,
                va="bottom", linespacing=1.3)

    ax.set_xscale("log")
    ax.set(xlim=(1.1e-3, 6.5e3), ylim=(-1.40, len(rows) - 0.28),
           yticks=list(range(len(rows))),
           xlabel="commanded amplitude [µm]")
    ax.set_yticklabels([r[0] for r in reversed(rows)], fontsize=8.5,
                       weight="bold")
    for tick, row in zip(reversed(ax.get_yticklabels()), rows):
        tick.set_color(COL[row[0]])
    ax.tick_params(axis="y", length=0)
    ax.set_title("(c) What each block interrogates", loc="left", pad=6)
    ax.grid(axis="x", color=C_GRID, lw=0.55)
    ax.spines[["left", "top", "right"]].set_visible(False)


def panel_budget(ax, cfg, sched, total_min) -> None:
    spend: dict[str, float] = {}
    for row in sched:
        span = float(row["t_end_s"]) - float(row["t_start_s"])
        spend[row["seg"]] = spend.get(row["seg"], 0.0) + span
    parts = [("transits and approaches", spend.get("move", 0.0)),
             ("analysed ramps", spend.get("ramp", 0.0)),
             ("single microsteps", spend.get("step", 0.0)),
             ("holds, dwells and settles", spend.get("dwell", 0.0))]
    key = {"transits and approaches": "move", "analysed ramps": "ramp",
           "single microsteps": "step", "holds, dwells and settles": "dwell"}
    left = 0.0
    for label, seconds in parts:
        width = seconds / 60.0
        ax.barh(0, width, left=left, height=0.50, color=BUDGET[key[label]],
                edgecolor="white")
        if width > 3.0:
            ax.text(left + width / 2.0, 0, f"{label}\n{width:.1f} min",
                    ha="center", va="center", fontsize=6.7, color="white",
                    weight="bold")
        left += width
    gate = cfg["session"]["watchdog_gate_min"]
    ax.axvline(gate, color="#b4473f", lw=1.0, ls="--")
    ax.text(gate - 0.4, -0.40, f"{gate:.0f} min gate", ha="right", va="top",
            fontsize=6.5, color="#8d3934")
    ax.set(yticks=[], xlabel="elapsed session time [min]", ylim=(-0.45, 0.45),
           xlim=(0, max(total_min, gate) * 1.06))
    ax.set_title("(d) Where the session goes", loc="left", pad=6)
    ax.spines[["left", "right", "top"]].set_visible(False)


# ------------------------------------------------------------------ figure

def figure(cfg, summary, sched) -> None:
    total = summary["total_min"]
    blocks = blocks_of(sched)

    fig = plt.figure(figsize=(7.45, 9.2))
    gs = fig.add_gridspec(4, 1, height_ratios=[0.55, 2.30, 1.95, 0.66],
                          left=0.105, right=0.975, top=0.893, bottom=0.168,
                          hspace=0.70)
    fig.text(0.105, 0.978, "Experiment F0b — friction scoping",
             fontsize=16, weight="bold", color=C_INK, va="top")
    fig.text(0.105, 0.940,
             "One session that sizes every friction effect the model will "
             "need, from a single 39 nm microstep to a 3 mm stroke",
             fontsize=8.5, color=C_MUTED, va="top")

    ax0 = fig.add_subplot(gs[0])
    panel_allocation(ax0, blocks, total)

    spend: dict[str, float] = {}
    for row in sched:
        span = float(row["t_end_s"]) - float(row["t_start_s"])
        spend[row["seg"]] = spend.get(row["seg"], 0.0) + span
    at_rest = spend.get("dwell", 0.0) / (total * 60.0)

    facts = [("segments", f"{summary['segments']}"),
             ("ladder strokes", f"{summary['counts']['ladder_strokes']}"
                                f" × 4 half-cycles × 2 passes"),
             ("slow ramps", f"{summary['counts']['slow_ramps']}"),
             ("smallest / largest stroke",
              f"{summary['amplitudes_um']['ladder_min']:.3f} / "
              f"{summary['amplitudes_um']['ladder_max']:.0f} µm"),
             ("predicted duration", f"{total:.1f} min"),
             ("margin to gate",
              f"{summary['margin_to_gate_min']:.1f} min")]
    columns = (0.105, 0.420, 0.735)
    for i, (k, v) in enumerate(facts):
        x = columns[i % 3]
        y = 0.828 - 0.0215 * (i // 3)
        fig.text(x, y, k, fontsize=7.0, color=C_MUTED, va="top")
        fig.text(x + 0.148, y, v, fontsize=7.0, color=C_INK, va="top")

    panel_trajectory(fig.add_subplot(gs[1], sharex=ax0), cfg, sched, total)
    panel_amplitudes(fig.add_subplot(gs[2]), cfg, summary)
    panel_budget(fig.add_subplot(gs[3]), cfg, sched, total)

    caption(fig, 0.105,
            f"The F0b plan as the host will generate it: {summary['segments']} segments, "
            f"{total:.1f} min, {summary['margin_to_gate_min']:.1f} min inside the 58 min gate. "
            f"(a) Eight blocks in the procedure's order, bracketed and split by three "
            f"reference pairs that carry the drift correction. (b) The session commands "
            f"displacements from 39 nm to 45 mm, so the axis is symmetric-log; the three "
            f"excursions to ±45 mm are the reference pairs and the plateau at +20 mm is the "
            f"second noise hold. The flat stretch after each reference pair is not inactivity "
            f"- it is N, L, S and D working within a few microns of x$_0$. (c) The same "
            f"session by amplitude rather than by time, which is the axis F0b is actually "
            f"built on: the ladder spans 1 to 76 800 microsteps and brackets the frozen "
            f"reversal constant β$_0$ by about a decade and a half either side, so the "
            f"saturation level it measures can be compared with it. (d) {at_rest * 100:.0f}% of "
            f"the session is spent at rest, by design: the noise floor and the dwell probe are "
            f"measurements made standing still, and every one of the "
            f"{summary['counts']['ladder_strokes']} ladder strokes and "
            f"{summary['counts']['slow_ramps']} slow ramps is preceded by a 3 mm memory wipe "
            f"and a settle.")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=320)
    plt.close(fig)


def luminance(hexstr: str) -> float:
    r, g, b = (int(hexstr[i:i + 2], 16) for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def main() -> int:
    style()
    cfg, summary, sched = load()
    figure(cfg, summary, sched)
    print(f"  {FIG.name}")
    seen = {}
    for k, v in COL.items():
        seen.setdefault(v, k if not k.startswith("R") else "R")
    print("  block colour luminance: " + ", ".join(
        f"{name} {luminance(hexstr):.0f}" for hexstr, name in seen.items()))
    print("  colour names the block; every block is also labelled in (a) "
          "and (c), so grey separation is not load-bearing")
    failed = [k for k, v in summary["checks"].items() if not v]
    if failed:
        print("  PLAN CHECKS FAILED: " + ", ".join(failed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
