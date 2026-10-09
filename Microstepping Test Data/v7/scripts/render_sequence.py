"""Render the v7 settling-error sequence defined in `sequence_design.py`.

Four panels, answering four different questions about the plan:

1. **Campaign** -- what the whole run looks like end to end, and how much of it
   is warm-up rather than measurement.
2. **One MRES block** -- the phase structure, at a scale where the individual
   moves are visible.
3. **One measurement** -- where the settling error actually gets read off, and
   why the first 200 ms of every dwell is thrown away.
4. **Target design** -- that the targets avoid aliasing with the lead screw,
   which is the one design choice that is invisible in a time plot.

Panel 3's ring-down is a **model**, drawn from v2.3's measured frequency and
damping; it is there to show where the measurement window sits relative to the
transient, not to predict an amplitude. Every other panel is the actual
commanded plan.

Run from this directory:

    python render_sequence.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sequence_design import (
    build_campaign, summarize, sweep_targets, comb_targets, comb_step_for,
    MRES_VALUES, COMB_CENTRES_UM, LEAD_PITCH_UM, FULL_STEP_UM,
    SWEEP_HALF_UM, DWELL_MEASURE_S, SETTLE_DISCARD_S, RATE_FULLSTEPS_S,
    MEASURE_RATE, N_SWEEP_TARGETS, WARMUP_MIN, OUT_DIR,
)

# Phase colours: the dataviz reference categorical order, slots 1/2/3/5/7.
# Yellow (slot 4) is deliberately skipped -- it is the slot that fails the
# all-pairs gate beside orange. Warm-up takes neutral grey rather than a
# categorical hue because it is the one phase that is NOT measured, and grey
# says that. Every band is additionally direct-labelled, so identity never
# rests on colour (and that is also the documented relief for magenta, which
# sits below 3:1 on white).
C_PHASE = {
    "unidirectional": "#2a78d6",
    "bidirectional": "#eb6834",
    "cyclic": "#1baf7a",
    "settling": "#e87ba4",
    "datum": "#4a3aa7",
    "witness": "#4a3aa7",
    "warmup": "#9aa3ab",
}
PHASE_LABEL = {
    "warmup": "Warm-up (not measured)",
    "datum": "Datum",
    "witness": "Drift witness",
    "unidirectional": "Unidirectional sweep",
    "bidirectional": "Bidirectional (reversal)",
    "cyclic": "Cyclic comb (¼ full step)",
    "settling": "Settling repeats",
}
C_INK = "#17202a"
C_INK2 = "#46515c"
C_MUTED = "#77818a"
C_GRID = "#dde3e8"
C_BASELINE = "#aeb7bf"
C_TRACE = "#2a3440"
LABEL_BBOX = dict(facecolor="white", edgecolor=C_GRID, linewidth=0.6,
                  alpha=0.95, pad=2.6)

RINGDOWN_HZ = 347.0        # v2.3: 343.8-350.0 Hz across four MRES values
RINGDOWN_TAU_S = 0.018     # v2.3: 15-21 ms
RINGDOWN_A_NM = 400.0      # illustrative amplitude -- NOT a measured value
SETTLING_ERROR_NM = 180.0  # illustrative offset, to show what gets read off


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": C_BASELINE,
        "axes.labelcolor": C_INK2,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlepad": 8,
        "axes.labelsize": 9.5,
        "axes.grid": False,
        "grid.color": C_GRID,
        "grid.linewidth": 0.7,
        "xtick.color": C_MUTED,
        "ytick.color": C_MUTED,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "font.size": 9.5,
        "text.color": C_INK,
        "figure.dpi": 150,
        "savefig.facecolor": "white",
    })


def style(ax, xlabel=None, ylabel=None, title=None, subtitle=None):
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(C_BASELINE)
    ax.grid(axis="y")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        # Padded for however many lines the subtitle runs to; without this the
        # two render on top of each other.
        lines = subtitle.count("\n") + 1 if subtitle else 0
        ax.set_title(title, loc="left", pad=14 + 11 * lines if subtitle else 8)
    if subtitle:
        ax.text(0.0, 1.015, subtitle, transform=ax.transAxes, fontsize=8.5,
                color=C_MUTED, ha="left", va="bottom", linespacing=1.4)


def polyline(segments):
    """Segments -> a single (t, x) polyline in minutes and mm."""
    t = [segments[0]["t_start_s"]]
    x = [segments[0]["x_start_um"]]
    for s in segments:
        t.append(s["t_end_s"])
        x.append(s["x_end_um"])
    return np.asarray(t) / 60.0, np.asarray(x) / 1000.0


def phase_runs(segments):
    """Contiguous runs of one phase, as (phase, t0_s, t1_s)."""
    runs = []
    for s in segments:
        if runs and runs[-1][0] == s["phase"]:
            runs[-1][2] = s["t_end_s"]
        else:
            runs.append([s["phase"], s["t_start_s"], s["t_end_s"]])
    return runs


def shade_phases(ax, segments, label=(), fontsize=8.0, min_label_s=45.0):
    """Shade each phase run; direct-label the longest run of named phases.

    Only phases in `label` get text, and only if their longest run is wide
    enough to hold it -- the datum and witness runs are seconds long and their
    labels would collide with whatever follows. They stay in the legend.
    """
    runs = phase_runs(segments)
    for phase, t0, t1 in runs:
        ax.axvspan(t0 / 60.0, t1 / 60.0, color=C_PHASE[phase], alpha=0.16,
                   lw=0, zorder=0)
    best = {}
    for phase, t0, t1 in runs:
        if phase not in label:
            continue
        if phase not in best or (t1 - t0) > (best[phase][1] - best[phase][0]):
            best[phase] = (t0, t1)
    for phase, (t0, t1) in best.items():
        if (t1 - t0) < min_label_s:
            continue
        ax.text((t0 + t1) / 120.0, 0.965, PHASE_LABEL[phase],
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=fontsize, color=C_INK2, bbox=LABEL_BBOX, zorder=6)


def panel_campaign(ax, seq, blocks, summary):
    t, x = polyline(seq.segments)
    shade_phases(ax, seq.segments)
    ax.plot(t, x, color=C_TRACE, lw=0.5, zorder=3)
    for b in blocks:
        mid = (b["t_start_s"] + b["t_end_s"]) / 120.0
        ax.axvline(b["t_start_s"] / 60.0, color=C_BASELINE, lw=0.9,
                   ls=(0, (4, 3)), zorder=2)
        ax.text(mid, 0.965, f"MRES {b['mres']}",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8.5, color=C_INK2, fontweight="semibold",
                bbox=LABEL_BBOX, zorder=6)
    wu = summary["timing"]["warmup_min"]
    ax.text(wu / 2.0, 0.965, "Warm-up\n(not measured)",
            transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=8.0, color=C_INK2, bbox=LABEL_BBOX, zorder=6)
    ax.set_xlim(0, t[-1])
    b = summary["budget"]
    style(ax, "Time [min]", "Commanded position [mm]",
          "The whole campaign, end to end",
          f"{b['total_min']:.0f} min total · {WARMUP_MIN:.0f} min warm-up + "
          f"{len(b['mres_values'])} MRES blocks · "
          f"{b['n_measurements']} measured positions")


def mark_zoom_windows(ax, windows, t_offset_s=0.0):
    """Show on a wide panel where each narrower panel below zooms into."""
    for i, (t0, t1, name) in enumerate(windows):
        ax.axvspan(t0 / 60.0, t1 / 60.0, color=C_INK, alpha=0.10, lw=0,
                   zorder=1)
        ax.annotate(name, xy=((t0 + t1) / 120.0, 0.012),
                    xycoords=ax.get_xaxis_transform(),
                    xytext=(0, 16 if i % 2 == 0 else 34),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=7.5, color=C_INK2, bbox=LABEL_BBOX, zorder=7,
                    arrowprops=dict(arrowstyle="-", color=C_MUTED, lw=0.8))


def panel_block(ax, seq, blocks, mres=16):
    b = next(x for x in blocks if x["mres"] == mres)
    segs = [s for s in seq.segments
            if s["t_start_s"] >= b["t_start_s"] - 1e-9
            and s["t_end_s"] <= b["t_end_s"] + 1e-9]
    t, x = polyline(segs)
    t0 = b["t_start_s"] / 60.0
    shade_phases(ax, segs, label=("unidirectional", "bidirectional", "cyclic",
                                  "settling"))
    ax.plot(t, x, color=C_TRACE, lw=0.7, zorder=3)
    ax.set_xlim(t0, b["t_end_s"] / 60.0)
    style(ax, "Time [min]", "Commanded position [mm]",
          f"One MRES block in detail (MRES {mres})",
          f"{(b['t_end_s'] - b['t_start_s']) / 60.0:.1f} min · every measured "
          "move at one rate (700 µm/s) · each flat top is a 1.5 s dwell · "
          "narrow violet bands are datum / drift witnesses")


def find_window(seq, blocks, mres, phase, nth=0, count=5, pad_s=2.0,
                approach=None):
    """Time window spanning `count` measured dwells of a phase, with padding.

    Picked from the plan rather than hard-coded seconds, so the zooms stay
    pointed at the right transitions if any timing constant changes.
    """
    b = next(x for x in blocks if x["mres"] == mres)
    dwells = [s for s in seq.segments
              if s["phase"] == phase and s["measured"]
              and b["t_start_s"] <= s["t_start_s"] <= b["t_end_s"]]
    if approach is not None:
        marks = [m for m in seq.measurements if m["phase"] == phase]
        keep = {round(m["window_end_s"], 4) for m in marks
                if m["approach"] == approach}
        dwells = [d for d in dwells if round(d["t_end_s"], 4) in keep]
    sel = dwells[nth:nth + count]
    if not sel:
        raise RuntimeError(f"no {phase} dwells found for MRES {mres}")
    return sel[0]["t_start_s"] - pad_s, sel[-1]["t_end_s"] + pad_s


def find_reversal_window(seq, blocks, mres, pad_s=3.0):
    """The instant the bidirectional pass turns around -- where backlash acts."""
    b = next(x for x in blocks if x["mres"] == mres)
    marks = [m for m in seq.measurements
             if m["phase"] == "bidirectional"
             and b["t_start_s"] <= m["window_start_s"] <= b["t_end_s"]]
    for i in range(1, len(marks)):
        if marks[i]["approach"] == "-" and marks[i - 1]["approach"] == "+":
            return (marks[i - 2]["window_start_s"] - pad_s,
                    marks[i + 1]["window_end_s"] + pad_s)
    raise RuntimeError("no +/- reversal found in the bidirectional pass")


def panel_zoom(ax, seq, t0, t1, title, subtitle, unit="s"):
    """A time window of the plan, with every measured dwell marked.

    This is the scale the campaign and block panels cannot show: one move, its
    arrival, and the hold that the position is actually read from.
    """
    segs = [s for s in seq.segments
            if s["t_end_s"] >= t0 and s["t_start_s"] <= t1]
    t = [segs[0]["t_start_s"]] + [s["t_end_s"] for s in segs]
    x = [segs[0]["x_start_um"]] + [s["x_end_um"] for s in segs]
    t = np.asarray(t) - t0
    x = np.asarray(x)
    span = t1 - t0

    for s in segs:
        if not s["measured"]:
            continue
        ax.axvspan(s["t_start_s"] - t0, s["t_end_s"] - t0,
                   color=C_PHASE[s["phase"]], alpha=0.22, lw=0, zorder=0)
        ax.axvspan(s["t_start_s"] - t0, s["t_start_s"] - t0 + SETTLE_DISCARD_S,
                   color=C_PHASE["settling"], alpha=0.30, lw=0, zorder=1)
    ax.plot(t, x, color=C_TRACE, lw=1.3, zorder=3,
            marker="o", ms=2.6, mew=0, mfc=C_TRACE)
    ax.set_xlim(0, span)
    style(ax, f"Time within window [{unit}]", "Commanded position [µm]",
          title, subtitle)


def panel_measurement(ax):
    """Anatomy of a single measured position -- where the error is read off."""
    v = RATE_FULLSTEPS_S[MEASURE_RATE] * FULL_STEP_UM      # um/s
    dt = 0.0001
    t_app = np.arange(-0.05, 0.0, dt)
    x_app = t_app * v * 1000.0                              # nm, target at 0

    t_dw = np.arange(0.0, DWELL_MEASURE_S, dt)
    ring = (RINGDOWN_A_NM * np.exp(-t_dw / RINGDOWN_TAU_S)
            * np.cos(2 * np.pi * RINGDOWN_HZ * t_dw))
    x_dw = SETTLING_ERROR_NM + ring
    envelope = RINGDOWN_A_NM * np.exp(-t_dw / RINGDOWN_TAU_S)

    ax.axvspan(0.0, SETTLE_DISCARD_S, color=C_PHASE["settling"], alpha=0.20,
               lw=0, zorder=0)
    ax.axvspan(SETTLE_DISCARD_S, DWELL_MEASURE_S,
               color=C_PHASE["unidirectional"], alpha=0.13, lw=0, zorder=0)
    ax.axhline(0.0, color=C_BASELINE, lw=1.1, zorder=2)
    ax.plot(t_app, x_app, color=C_TRACE, lw=1.4, zorder=3)
    ax.plot(t_dw, x_dw, color=C_TRACE, lw=0.6, zorder=3)
    # The envelope, so the decay is legible where the oscillation itself packs
    # 347 cycles a second into a few millimetres of axis.
    ax.plot(t_dw, SETTLING_ERROR_NM + envelope, color=C_PHASE["settling"],
            lw=1.2, zorder=4)
    ax.plot(t_dw, SETTLING_ERROR_NM - envelope, color=C_PHASE["settling"],
            lw=1.2, zorder=4)
    ax.axhline(SETTLING_ERROR_NM, color=C_PHASE["bidirectional"], lw=1.3,
               ls=(0, (5, 3)), zorder=5)

    ax.annotate("", xy=(1.08, SETTLING_ERROR_NM), xytext=(1.08, 0.0),
                arrowprops=dict(arrowstyle="<->",
                                color=C_PHASE["bidirectional"], lw=1.4),
                zorder=6)
    ax.text(1.12, SETTLING_ERROR_NM / 2 - 40,
            "settling error\n= mean(window) − commanded",
            fontsize=8.0, color=C_PHASE["bidirectional"],
            fontweight="semibold", va="top", zorder=7, bbox=LABEL_BBOX)
    ax.text(SETTLE_DISCARD_S / 2, 0.045, "discard\n0.20 s",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=8.0, color=C_INK2, bbox=LABEL_BBOX, zorder=7)
    ax.text(0.72, 0.045, "averaging window — 1.30 s, 1300 samples at 1 kHz",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=8.0, color=C_INK2, bbox=LABEL_BBOX, zorder=7)
    ax.text(0.015, 0.5, "commanded target", transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=8.0, color=C_INK2,
            bbox=LABEL_BBOX, zorder=7)

    ax.set_xlim(-0.05, DWELL_MEASURE_S)
    ax.set_ylim(-620, 900)
    style(ax, "Time from arrival at target [s]",
          "Position relative to target [nm]",
          "One measurement: where the settling error is read off",
          f"ring-down modelled at {RINGDOWN_HZ:.0f} Hz, τ = "
          f"{RINGDOWN_TAU_S * 1000:.0f} ms (v2.3's measured decay); "
          "amplitude and error illustrative")

    # Inset: the first 60 ms at a scale where the oscillation is visible.
    axi = ax.inset_axes([0.33, 0.56, 0.42, 0.40])
    m = t_dw <= 0.06
    axi.axvspan(0.0, 0.06, color=C_PHASE["settling"], alpha=0.20, lw=0)
    axi.plot(t_dw[m] * 1000, x_dw[m], color=C_TRACE, lw=0.8)
    axi.plot(t_dw[m] * 1000, SETTLING_ERROR_NM + envelope[m],
             color=C_PHASE["settling"], lw=1.1)
    axi.plot(t_dw[m] * 1000, SETTLING_ERROR_NM - envelope[m],
             color=C_PHASE["settling"], lw=1.1)
    axi.axhline(SETTLING_ERROR_NM, color=C_PHASE["bidirectional"], lw=0.9,
                ls=(0, (4, 3)))
    axi.set_xlim(0, 60)
    axi.set_title("first 60 ms of the discard", fontsize=7.5, loc="left",
                  pad=3, color=C_INK2)
    axi.tick_params(labelsize=7, colors=C_MUTED)
    axi.set_xlabel("ms", fontsize=7.5, labelpad=1)
    for s in ("top", "right"):
        axi.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        axi.spines[s].set_color(C_BASELINE)


def panel_targets(ax):
    """Do the targets avoid aliasing with the 2 mm lead screw pitch?"""
    targets = sweep_targets()
    uniform = np.linspace(-SWEEP_HALF_UM, SWEEP_HALF_UM, N_SWEEP_TARGETS)

    for c in COMB_CENTRES_UM:
        ax.axvspan((c - 400) / 1000.0, (c + 400) / 1000.0,
                   color=C_PHASE["cyclic"], alpha=0.18, lw=0, zorder=0)
    ax.scatter(uniform / 1000.0, uniform % LEAD_PITCH_UM, s=42,
               facecolor="white", edgecolor=C_MUTED, linewidth=1.1, zorder=3,
               label="A uniform 2 mm grid — every target at the same screw phase")
    ax.scatter(targets / 1000.0, targets % LEAD_PITCH_UM, s=46,
               color=C_PHASE["unidirectional"], edgecolor="white",
               linewidth=0.9, zorder=4,
               label="This plan — jittered, then snapped to whole full steps")

    ax.set_xlim(-SWEEP_HALF_UM / 1000.0 - 2, SWEEP_HALF_UM / 1000.0 + 2)
    ax.set_ylim(-150, LEAD_PITCH_UM + 1050)
    style(ax, "Target position along travel [mm]",
          "Position within one lead-screw turn [µm]",
          "Target design: spreading the lead-screw phase",
          f"{targets.size} global targets · green bands are the three "
          "¼-full-step combs · a uniform grid would sample one screw phase "
          "only")
    ax.legend(loc="upper left", ncol=1)


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    seq, blocks, n_legs = build_campaign()
    summary = summarize(seq, blocks, n_legs, MRES_VALUES)

    ZOOM_MRES = 16
    w_uni = find_window(seq, blocks, ZOOM_MRES, "unidirectional", nth=2,
                        count=4)
    w_rev = find_reversal_window(seq, blocks, ZOOM_MRES)
    w_comb = find_window(seq, blocks, ZOOM_MRES, "cyclic", nth=6, count=6,
                         pad_s=1.0)

    fig = plt.figure(figsize=(14.2, 16.4))
    gs = fig.add_gridspec(4, 3, height_ratios=[0.95, 0.95, 0.92, 1.05],
                          hspace=0.60, wspace=0.27,
                          left=0.062, right=0.98, top=0.898, bottom=0.072)

    panel_campaign(fig.add_subplot(gs[0, :]), seq, blocks, summary)
    ax_block = fig.add_subplot(gs[1, :])
    panel_block(ax_block, seq, blocks, mres=ZOOM_MRES)
    mark_zoom_windows(ax_block, [
        (*w_uni, "zoom A"), (*w_rev, "zoom B"), (*w_comb, "zoom C")])

    panel_zoom(fig.add_subplot(gs[2, 0]), seq, *w_uni,
               "A — unidirectional: move, arrive, hold",
               "four consecutive targets · each shaded band is a\n"
               "1.5 s dwell, its darker head the 0.20 s discard")
    panel_zoom(fig.add_subplot(gs[2, 1]), seq, *w_rev,
               "B — the reversal, where backlash acts",
               "the pass turns around · the 500 µm back-off puts\n"
               "the final approach on the other side")
    panel_zoom(fig.add_subplot(gs[2, 2]), seq, *w_comb,
               "C — the ¼-full-step comb",
               "2.5 µm commanded increments — the transitions\n"
               "the two panels above cannot resolve")

    panel_measurement(fig.add_subplot(gs[3, 0:2]))
    panel_targets(fig.add_subplot(gs[3, 2]))

    handles = [Line2D([0], [0], color=C_PHASE[p], lw=7, alpha=0.5,
                      label=PHASE_LABEL[p])
               for p in ("warmup", "datum", "unidirectional", "bidirectional",
                         "cyclic", "settling")]
    fig.legend(handles=handles, loc="lower center", ncol=6,
               bbox_to_anchor=(0.5, 0.012), frameon=False, fontsize=9)

    fig.suptitle("v7 — proposed settling-error measurement sequence",
                 x=0.062, ha="left", y=0.982, fontsize=17,
                 fontweight="semibold")
    fig.text(0.062, 0.952,
             "How well does the stage reach a commanded position? "
             "Every phase below exists to keep one confounder "
             "(backlash, thermal drift, cyclic error) out of that number.",
             ha="left", fontsize=10, color=C_MUTED)

    for ext in ("png", "svg"):
        fig.savefig(OUT_DIR / f"experimental_sequence.{ext}", dpi=150,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"wrote experimental_sequence.png/.svg -> {OUT_DIR}")


if __name__ == "__main__":
    main()
