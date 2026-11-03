#!/usr/bin/env python3
"""Two figures that explain what E5 does.

`e5_sequence_overview.png` answers "how is the session executed" in the same
form as the E1-E4 sequence renders: controller-state allocation, the complete
commanded trajectory, and where the hour goes.

`e5_compensation_chain.png` answers "how does the model improve the stage".
The correction is feed-forward and open loop - it is computed on the host
before the run and baked into a command table, and the interferometer never
enters the loop - so the figure is a signal path with a worked target on it,
not a control loop with feedback.

    python render_e5_sequence.py
"""
from __future__ import annotations

import csv
import json
import math
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered"
SEQUENCE = OUT / "e5_sequence_overview.png"
CHAIN = OUT / "e5_compensation_chain.png"

# Chapter convention for E5: uncompensated grey, compensated blue.
C_U, C_C = "#7b858c", "#176d8c"
C_REF, C_B3 = "#73578f", "#a16a2d"
C_INK, C_MUTED, C_GRID = "#20262b", "#566069", "#dde1e4"
COL = {"HOME": "#4c5560", "INIT": "#69747c", "R1": C_REF, "R2": C_REF,
       "R3": C_REF, "R4": C_REF, "B2A": C_C, "B2B": C_C, "B3": C_B3,
       "PARK": "#4c5560", "settle": "#c59b50", "record": "#67527c",
       "motion": "#697b86"}
CHAR_INCH = 0.0455


def caption(fig, left, text, bottom=0.013, size=6.9) -> None:
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


def tidy(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=C_GRID, lw=0.5)


def load():
    config = json.loads((HERE / "e5_config.json").read_text(encoding="utf-8"))
    summary = json.loads((HERE / "sequence/e5_plan_summary.json")
                         .read_text(encoding="utf-8"))
    with (HERE / "sequence/e5_measurement_plan.csv").open(
            newline="", encoding="utf-8") as fh:
        plan = list(csv.DictReader(fh))
    with (HERE / "sequence/e5_command_table.csv").open(
            newline="", encoding="utf-8") as fh:
        table = list(csv.DictReader(fh))
    return config, summary, plan, table


def move_time(distance_mm: float, c: dict) -> float:
    distance = abs(distance_mm)
    v, a = c["velocity_mm_s"], c["acceleration_mm_s2"]
    ramp = v * v / a
    return 2.0 * math.sqrt(distance / a) if distance <= ramp else 2.0 * v / a + (distance - ramp) / v


def fill_modes(plan: list[dict]) -> list[str]:
    """Mode for every row, carried back onto the moves that lead into a stop.

    Moves have no mode of their own, but a move belongs to the run it is
    approaching: without this the whole of B2 draws in one colour and the
    U/C alternation - the reason the session is built this way - disappears.
    """
    modes = [row["mode"] for row in plan]
    pending = ""
    for i in range(len(plan) - 1, -1, -1):
        if plan[i]["action"] == "state_enter":
            pending = ""
        elif modes[i]:
            pending = modes[i]
        else:
            modes[i] = pending
    return modes


def schedule(plan: list[dict], c: dict):
    """Replay the plan on the generator's own timing model."""
    modes = fill_modes(plan)
    now, position, state = 0.0, c["park_mm"], "HOME"
    segments, blocks, stops = [], [], []
    block_start, block_mode, block_run = 0.0, "", ""
    for row, mode in zip(plan, modes):
        if row["action"] == "state_enter":
            if now > block_start:
                blocks.append((block_start, now, state, block_mode, block_run))
            state, block_start = row["state"], now
            block_mode, block_run = "", ""
            continue
        if not block_mode:
            block_mode, block_run = mode, row["run"]
        if row["x_nom_mm"]:
            target = float(row["x_nom_mm"])
            end = now + move_time(target - position, c)
            segments.append((now, end, position, target, state, mode))
            now, position = end, target
        if row["action"] == "stop":
            hold = c["settle_s"] + c["record_s"]
            segments.append((now, now + hold, position, position, state, mode))
            stops.append((now, now + hold, position, state, mode))
            now += hold
    if now > block_start:
        blocks.append((block_start, now, state, block_mode, block_run))
    return segments, blocks, stops, now


def segment_colour(state: str, mode: str) -> str:
    """ISO runs are coloured by what they command, everything else by block."""
    if state.startswith("B2"):
        return C_C if mode == "C" else C_U
    return COL.get(state, "#7a7f84")


# ------------------------------------------------------------- figure 1

def figure_sequence(c, summary, segments, blocks, stops, total) -> None:
    minutes = total / 60.0
    fig = plt.figure(figsize=(7.45, 6.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.74, 3.05, 0.92],
                          left=0.100, right=0.975, top=0.862, bottom=0.150,
                          hspace=0.52)
    fig.text(0.100, 0.975, "Experiment E5 — compensated validation", fontsize=16,
             weight="bold", color=C_INK, va="top")
    fig.text(0.100, 0.932,
             "Ten ISO 230-2 runs, five uncompensated and five compensated, "
             "interleaved ABBA in one session",
             fontsize=8.5, color=C_MUTED, va="top")

    ax0 = fig.add_subplot(gs[0])
    for start, end, state, mode, run in blocks:
        width = (end - start) / 60.0
        ax0.barh(0, width, left=start / 60.0, height=0.42,
                 color=segment_colour(state, mode), edgecolor="white", linewidth=0.6)
        if width > 1.1:
            # Name the ISO runs by what they command, so the ABBA order reads
            # straight off the bar rather than having to be taken on trust.
            label = f"{mode}{run}" if state.startswith("B2") else state
            ax0.text((start + end) / 120.0, 0, label, ha="center", va="center",
                     fontsize=6.4, color="white", weight="bold")
    ax0.set(xlim=(0, minutes), ylim=(-0.36, 0.36), yticks=[])
    ax0.set_title("(a) Controller-state allocation", loc="left", pad=6)
    ax0.tick_params(axis="x", bottom=False, labelbottom=False)
    ax0.spines[:].set_visible(False)

    ax = fig.add_subplot(gs[1], sharex=ax0)
    grouped: dict[str, list] = {}
    for a, b, x0, x1, st, mode in segments:
        grouped.setdefault(segment_colour(st, mode), []).append(
            [(a / 60.0, x0), (b / 60.0, x1)])
    for colour, lines in grouped.items():
        ax.add_collection(LineCollection(lines, colors=colour, linewidths=0.9))
    for x in c["reference_positions_mm"]:
        ax.axhline(x, color=C_REF, lw=0.6, ls=(0, (3, 2)), alpha=0.8)
    ax.text(0.4, c["reference_positions_mm"][1] + 2.5,
            "reference pair, 15 and 105 mm", fontsize=6.6, color=C_REF)
    ax.set(xlim=(0, minutes), ylim=(0, 124), ylabel="commanded position [mm]")
    ax.set_title("(b) Complete commanded trajectory", loc="left", pad=6)
    ax.grid(axis="x", color=C_GRID, lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], color=C_U, lw=1.6, label="uncompensated run"),
               plt.Line2D([], [], color=C_C, lw=1.6, label="compensated run"),
               plt.Line2D([], [], color=C_REF, lw=1.6, label="reference block"),
               plt.Line2D([], [], color=C_B3, lw=1.6, label="B3 off-grid")]
    ax.legend(handles=handles, frameon=False, fontsize=6.8, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.055),
              handlelength=1.6, columnspacing=1.6)
    box = (f"ISO targets  24 (E0 v4)\n"
           f"runs  5 U + 5 C, order {c['b2_mode_order']}\n"
           f"reference stops  {summary['stops']['reference']} in 4 blocks\n"
           f"B3 off-grid stops  {summary['stops']['B3']}\n"
           f"recorded stops  {summary['stops']['total']}\n"
           f"predicted duration  {minutes:.1f} min")
    ax.text(0.015, 0.975, box, transform=ax.transAxes, ha="left", va="top",
            fontsize=7.2, linespacing=1.35,
            bbox=dict(boxstyle="square,pad=0.48", fc="white", ec="#aeb4b9", lw=0.7))

    ax2 = fig.add_subplot(gs[2], sharex=ax0)
    n = summary["stops"]["total"]
    parts = [("motion + transits", total - n * (c["settle_s"] + c["record_s"]),
              COL["motion"]),
             ("settle windows", n * c["settle_s"], COL["settle"]),
             ("record", n * c["record_s"], COL["record"])]
    left = 0.0
    for label, seconds, colour in parts:
        width = seconds / 60.0
        ax2.barh(0, width, left=left, height=0.50, color=colour, edgecolor="white")
        if width > 2.0:
            ax2.text(left + width / 2.0, 0, f"{label}\n{width:.1f} min",
                     ha="center", va="center", fontsize=6.7, color="white",
                     weight="bold")
        left += width
    ax2.axvline(c["watchdog_gate_min"], color="#b4473f", lw=1.0, ls="--")
    ax2.text(c["watchdog_gate_min"] - 0.4, -0.40,
             f"{c['watchdog_gate_min']:.0f} min gate", ha="right", va="top",
             fontsize=6.5, color="#8d3934")
    ax2.set(yticks=[], xlabel="elapsed session time [min]", ylim=(-0.45, 0.45),
            xlim=(0, max(minutes, c["watchdog_gate_min"]) * 1.06))
    ax2.set_title("(c) Where the session goes", loc="left", pad=6)
    ax2.spines[["left", "right", "top"]].set_visible(False)

    caption(fig, 0.100,
            f"One session, {minutes:.1f} min, {summary['stops']['total']} recorded stops. The ten "
            f"ISO runs alternate {c['b2_mode_order']} so that any drift over the session falls "
            f"equally on both modes and the U/C comparison stays paired. Four reference blocks "
            f"bracket the ISO blocks; each visits 15 and 105 mm from both directions, giving drift "
            f"anchors, a session gain estimate, and reversal at fixed positions against time. B3 "
            f"then sweeps the electrical cycle off the 64-µstep grid. Only the C runs and B3 carry "
            f"a correction; U and the references are commanded at their nominal targets.")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(SEQUENCE, dpi=320)
    plt.close(fig)


# ------------------------------------------------------------- figure 2

def box(ax, x, y, w, h, text, fc, ec, fontsize=6.9, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=fc, ec=ec, lw=0.9, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, color=C_INK, zorder=3, linespacing=1.35,
            weight=weight)


def arrow(ax, x0, y0, x1, y1, style="-|>", colour=C_INK, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=9, lw=0.9, color=colour,
                                 linestyle=ls, zorder=2,
                                 shrinkA=0, shrinkB=0))


def figure_chain(c, summary, table) -> None:
    rows = [r for r in table if r["mode"] == "C" and r["block"] == "B2"]
    # An interior target: the edge-flagged ones sit where the spline is weakly
    # constrained, the worst place to take a worked example from.
    interior = [r for r in rows if r.get("edge") != "1"] or rows
    pick = max(interior, key=lambda r: abs(float(r["G_rel_um"])))
    target = pick["target_id"]
    pair = {r["dir"]: r for r in rows if r["target_id"] == target
            and r["run"] == pick["run"]}
    beta0 = summary["up_minus_down_um"]["expected"]
    x_nom = float(pair["1"]["x_nom_mm"])
    g = float(pair["1"]["G_rel_um"])

    fig = plt.figure(figsize=(7.8, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.52, 1.0],
                          left=0.085, right=0.975, top=0.885, bottom=0.215,
                          hspace=0.40, wspace=0.34)
    fig.text(0.085, 0.975, "E5 — how the model is applied to the stage",
             fontsize=16, weight="bold", color=C_INK, va="top")
    fig.text(0.085, 0.932,
             "Feed-forward pre-correction: computed on the host before the run, "
             "baked into a command table",
             fontsize=8.5, color=C_MUTED, va="top")

    ax = fig.add_subplot(gs[0, :])
    ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[])
    ax.spines[:].set_visible(False)

    y, h = 0.745, 0.205
    box(ax, 0.005, y, 0.135, h, "nominal\nISO target\n$x^*$", "#f2f4f5", "#9aa3a9")
    box(ax, 0.175, y, 0.175, h,
        "frozen model\n$\\hat{G}(x^*,\\theta)$\nH = 2, E2 only", "#e8f0f3", C_C)
    box(ax, 0.385, y, 0.165, h,
        "subtract\ngeometry\n$-[\\hat{G}-a_0]$", "#e8f0f3", C_C)
    box(ax, 0.585, y, 0.165, h,
        "subtract\nreversal\n$-\\frac{d}{2}\\beta_0$", "#e8f0f3", C_C)
    box(ax, 0.785, y, 0.21, h,
        "round to µstep\n$u_{cmd}=\\mathrm{round}(x_c\\!\\cdot\\!1600)$",
        "#f6efe3", C_B3)
    for x0, x1 in ((0.140, 0.175), (0.350, 0.385), (0.550, 0.585), (0.750, 0.785)):
        arrow(ax, x0, y + h / 2, x1, y + h / 2)

    y2, hh = 0.435, 0.175
    box(ax, 0.785, y2, 0.21, hh, "ESP32 RMT\nstep pulses", "#f2f4f5", "#9aa3a9")
    arrow(ax, 0.890, y, 0.890, y2 + hh)
    box(ax, 0.515, y2, 0.21, hh, "stage\n(screw, nut, motor)", "#f2f4f5", "#9aa3a9")
    arrow(ax, 0.785, y2 + hh / 2, 0.725, y2 + hh / 2)
    box(ax, 0.235, y2, 0.22, hh,
        "interferometer\nmeasures the result", "white", "#9aa3a9")
    arrow(ax, 0.515, y2 + hh / 2, 0.455, y2 + hh / 2)

    box(ax, 0.005, 0.015, 0.465, 0.265,
        f"worked target {target}: $x^*$ = {x_nom:.3f} mm\n"
        f"$\\hat{{G}}$ = {g:+.3f} µm,  $\\beta_0$ = −{beta0:.3f} µm\n"
        f"up    $u$ {pair['1']['u_nom']} → {pair['1']['u_cmd']}"
        f"   ({float(pair['1']['corr_um']):+.3f} µm)\n"
        f"down  $u$ {pair['-1']['u_nom']} → {pair['-1']['u_cmd']}"
        f"   ({float(pair['-1']['corr_um']):+.3f} µm)",
        "#fbfbfc", "#aeb4b9", fontsize=6.9)

    arrow(ax, 0.235, y2 + hh / 2, 0.152, y2 + hh / 2, style="-|>",
          colour="#b4473f", ls=(0, (3, 2)))
    ax.text(0.148, y2 + hh / 2, "validation only ", ha="right", va="center",
            fontsize=6.6, color="#8d3934")
    ax.text(0.515, 0.255,
            "No feedback. The model is evaluated once per stop on the\n"
            "host, before the session. During the run the board replays\n"
            "step counts and nothing more; the interferometer records\n"
            "the outcome without ever influencing a command.",
            ha="left", va="top", fontsize=6.9, color="#8d3934",
            linespacing=1.55, style="italic")
    ax.set_title("(a) The signal path for one commanded stop", loc="left", pad=4)

    # (b) what gets commanded along the travel
    ax = fig.add_subplot(gs[1, 0])
    up = sorted(((float(r["x_nom_mm"]), float(r["corr_um"]))
                 for r in rows if r["dir"] == "1"), key=lambda t: t[0])
    dn = sorted(((float(r["x_nom_mm"]), float(r["corr_um"]))
                 for r in rows if r["dir"] == "-1"), key=lambda t: t[0])
    geo = sorted(((float(r["x_nom_mm"]), -float(r["G_rel_um"]))
                  for r in rows if r["dir"] == "1"), key=lambda t: t[0])
    ax.plot([a for a, _ in geo], [b for _, b in geo], lw=1.1, color=C_INK,
            label=r"$-\hat{G}(x)$")
    ax.plot([a for a, _ in up], [b for _, b in up], ls="none", marker="^", ms=3.2,
            color=C_C, label="commanded, up")
    ax.plot([a for a, _ in dn], [b for _, b in dn], ls="none", marker="v", ms=3.2,
            mfc="none", mec=C_C, mew=0.9, label="commanded, down")
    ax.set(xlabel="target position [mm]", ylabel="correction [µm]")
    ax.set_title("(b) What each run commands", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=6.4, loc="upper left", handlelength=1.5,
              labelspacing=0.25)
    tidy(ax)

    # (c) quantisation
    ax = fig.add_subplot(gs[1, 1])
    q = [float(r["q_um"]) for r in table if r["mode"] == "C"]
    ax.hist(q, bins=24, color=C_B3, edgecolor="white", linewidth=0.5)
    limit = 1000.0 / 1600.0 / 2.0
    for sign in (-1, 1):
        ax.axvline(sign * limit, color=C_INK, lw=0.9, ls=(0, (3, 2)))
    ax.text(0.5, 0.96, f"half a µstep = ±{limit:.3f} µm", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.5, color=C_INK)
    ax.set(xlabel="quantisation residual q [µm]", ylabel="stops")
    ax.set_title("(c) The price of whole µsteps", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=C_GRID, lw=0.5)

    # (d) what that should leave
    ax = fig.add_subplot(gs[1, 2])
    q_rms = summary["quantisation_um"]["rms"]
    e3_transfer, e3_floor = 0.522, 0.369
    expected = math.hypot(e3_transfer, q_rms / math.sqrt(2))
    bars = [("uncompensated\n(E3 measured)", 1.489, C_U),
            ("expected after\ncompensation", expected, C_C),
            ("E3 transfer\n(offline)", e3_transfer, "#9aa3a9"),
            ("same-day limit", e3_floor, "#c3c9cd")]
    ypos = np.arange(len(bars))[::-1]
    ax.barh(ypos, [v for _, v, _ in bars], height=0.6,
            color=[col for _, _, col in bars])
    for yy, (_, v, _) in zip(ypos, bars):
        ax.text(v + 0.03, yy, f"{v:.2f}", va="center", fontsize=6.4, color=C_INK)
    ax.set_yticks(ypos)
    ax.set_yticklabels([n for n, _, _ in bars], fontsize=6.3)
    ax.set(xlabel="µm rms", xlim=(0, 1.75))
    ax.set_title("(d) What it should leave", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=C_GRID, lw=0.5)

    caption(fig, 0.085,
            f"(a) The correction is computed once per stop on the host and stored as a step count; "
            f"the board replays the table and never evaluates the model, so nothing measured during "
            f"the run changes what is commanded. The two correction terms are separate: "
            f"$-[\\hat{{G}}-a_0]$ depends on position and screw angle, $-\\frac{{d}}{{2}}\\beta_0$ "
            f"only on approach direction, which is why the up and down commands for one target "
            f"differ by {beta0:.3f} µm. (b) Across the travel the correction spans "
            f"{summary['correction_um']['min']:+.2f} to {summary['correction_um']['max']:+.2f} µm; "
            f"the up and down markers straddle the geometry curve by half the reversal term each. "
            f"(c) Commands are whole µsteps, so a residual of up to ±{limit:.3f} µm is unavoidable "
            f"({q_rms:.3f} µm rms here). (d) Adding that to E3's measured cross-day transfer error "
            f"predicts roughly {expected:.2f} µm rms — above the same-day limit, and the number E5 "
            f"is really testing.")
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(CHAIN, dpi=320)
    plt.close(fig)


def main() -> int:
    config, summary, plan, table = load()
    style()
    segments, blocks, stops, total = schedule(plan, config)
    if abs(total / 60.0 - summary["duration_min"]) > 1e-6:
        raise RuntimeError("render schedule disagrees with the generated plan")
    figure_sequence(config, summary, segments, blocks, stops, total)
    figure_chain(config, summary, table)
    print(f"wrote {SEQUENCE.name} and {CHAIN.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
