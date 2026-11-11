#!/usr/bin/env python3
"""The F0b results figures, to the rendering instructions.

    python3 render_f0b_results.py

writes, as vector PDF and 300 dpi PNG, into rendered/f0b/

    f0b_1_noise_floor      what the measurement can resolve
    f0b_2_lost_motion      how much a reversal costs, and over what distance
    f0b_3_slow_ramps       what happens microstep by microstep
    f0b_4_dwell_speed      whether rest time or speed change the response
    f0b_A1_validity        drift anchors, windows, MSCNT, checks

and the scales card as f0b/f0b_scales_card.md.

Findings live in the captions, which are printed to stdout and written to
rendered/f0b/captions.md; the figures carry no titles and no statistics
boxes.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import f0b_process as P
from thesis_style import (C, FULL_WIDTH, HALF_WIDTH, LABEL, apply, finish,
                          panel_letter, ref_line, save)

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "f0b"
DODGE = 1.015            # 1.5 % of the axis, for up/down at the same x


def dminband(ax, value_um: float, label: str = "below detection") -> None:
    ax.axhspan(-value_um, value_um, color=C["ref"], alpha=0.18, lw=0, zorder=0)
    ax.text(0.012, value_um, label, transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=8, color=C["ref"])


# -- F0b-1 ----------------------------------------------------------------
def figure_noise(d: dict) -> tuple[list[Path], str]:
    nz = d["noise"]
    fig, ax = plt.subplots(figsize=(HALF_WIDTH, 2.6))
    styles = {"N-mid": ("-", 1.6), "N-end1": ((0, (4, 2)), 1.0),
              "N-end2": ((0, (1, 2)), 1.0)}
    for hold in nz["holds"]:
        tau = np.array(hold["taus_s"])
        dev = np.array([np.nan if v is None else v for v in hold["allan_nm"]])
        ls, lw = styles[hold["name"]]
        ax.plot(tau, dev, ls=ls, lw=lw, color=C["N"],
                label=f"{hold['name']} ({hold['position_mm']:.0f} mm)")
        k = hold["taus_s"].index(0.5)
        ax.plot(0.5, dev[k], "o", mfc="white", mec=C["N"], mew=0.9, ms=5)

    first = np.array(nz["holds"][0]["allan_nm"][:3], float)
    t0 = np.array(nz["holds"][0]["taus_s"][:3])
    scale = float(np.mean(first * np.sqrt(t0)))
    tt = np.array([1e-3, 0.3])
    ax.plot(tt, scale / np.sqrt(tt), ls=(0, (3, 2)), lw=0.7, color=C["ref"])
    ax.text(1.4e-3, scale / np.sqrt(1.4e-3) * 0.55, "white noise", fontsize=8,
            color=C["ref"], rotation=-28)

    ax.axhspan(0, nz["d_min_nm"], color=C["ref"], alpha=0.18, lw=0, zorder=0)
    ax.text(0.985, nz["d_min_nm"], "minimum detectable step",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=8, color=C["ref"])
    ref_line(ax, P.USTEP_NM, "one microstep", side="right")

    ax.set(xscale="log", yscale="log", xlim=(1e-3, 10), ylim=(0.2, 90),
           xlabel=LABEL["tau"], ylabel=LABEL["precision"])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncols=3,
              handlelength=1.8, columnspacing=1.0)
    finish(ax)
    panel_letter(ax, "a", dx=-0.17)

    mid = next(h for h in nz["holds"] if h["name"] == "N-mid")
    caption = (
        "**F0b-1.** Measurement floor at rest. (a) Precision of a window mean "
        "against averaging time at the middle and both ends of the travel, as "
        "the overlapping Allan deviation after removing a linear trend from "
        f"each hold; the dashed line is the white-noise slope. At the 0.5 s "
        f"record window (open markers) the precision is {nz['sigma_500ms_nm']:.2f} nm, "
        f"giving a minimum detectable step of {nz['d_min_nm']:.1f} nm, "
        f"{P.USTEP_NM / nz['d_min_nm']:.0f} times below one microstep. The floor "
        f"does not depend on position: the three holds agree within "
        f"{nz['position_spread_nm']:.2f} nm at 0.5 s. Averaging stops helping "
        f"beyond about 0.2 s, where drift takes over at "
        f"{mid['drift_nm_per_s'] * 60:+.0f} nm/min at the middle of the travel.")
    return save(fig, OUT, "f0b_1_noise_floor"), caption


# -- F0b-2 ----------------------------------------------------------------
def first_halves(halves, cap_um):
    """The first half-cycle of each stroke: it starts from whatever the
    previous stroke left behind, not from a closed loop, so it is shown
    separately rather than averaged into the reversal points."""
    by: dict[float, list] = {}
    for h in halves:
        if h.index != 1 or h.pass_name != "increasing":
            continue
        by.setdefault(round(h.stroke_um, 4), []).append(h)
    s = np.array(sorted(k for k in by if k <= cap_um))
    lost = np.array([np.mean([q.lost_travel_um for q in by[k]]) for k in s])
    err = np.array([2 * np.std([q.lost_travel_um for q in by[k]]) for k in s])
    frac = np.array([np.mean([q.tau for q in by[k]]) for k in s])
    ferr = np.array([2 * np.std([q.tau for q in by[k]]) for k in s])
    return s, lost, err, frac, ferr


def figure_lost_motion(d: dict) -> tuple[list[Path], str]:
    lad, halves = d["ladder"], d["halves"]
    d_min_um = d["noise"]["d_min_nm"] / 1000
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH, 3.0),
                             gridspec_kw={"width_ratios": [1.3, 1]})
    ax, bx = axes

    def series(pass_name, kind="reversal"):
        by: dict[float, list] = {}
        for h in halves:
            if h.kind != kind or h.pass_name != pass_name or h.index < 2:
                continue
            by.setdefault(round(h.stroke_um, 4), []).append(h)
        s = np.array(sorted(by))
        keep = s <= lad["usable_above_um"]
        s = s[keep]
        lost = np.array([np.mean([q.lost_travel_um for q in by[k]]) for k in s])
        err = np.array([2 * np.std([q.lost_travel_um for q in by[k]]) for k in s])
        frac = np.array([np.mean([q.tau for q in by[k]]) for k in s])
        ferr = np.array([2 * np.std([q.tau for q in by[k]]) for k in s])
        return s, lost, err, frac, ferr

    s_i, l_i, e_i, f_i, fe_i = series("increasing")
    s_d, l_d, e_d, f_d, fe_d = series("decreasing")
    s_c, l_c, e_c, f_c, fe_c = first_halves(halves, lad["usable_above_um"])

    ax.plot([0.02, 5], [0.02, 5], ls=(0, (4, 2)), lw=0.7, color=C["ref"])
    ax.text(0.055, 0.030, "nothing arrives", fontsize=8, color=C["ref"],
            rotation=34)
    ax.errorbar(s_i, l_i, yerr=e_i, fmt="o", ms=4, color=C["L"], capsize=2,
                elinewidth=0.8, label="increasing pass")
    ax.errorbar(s_d * DODGE, l_d, yerr=e_d, fmt="o", ms=4, mfc="white",
                mec=C["L"], ecolor=C["L"], capsize=2, elinewidth=0.8,
                label="decreasing pass")
    ax.plot(s_c, l_c, "D", ms=3.2, color=C["ref"], ls="none",
            label="first half-cycle of each stroke")
    ref_line(ax, lad["saturation_um"],
             f"$\\ell_{{sat}}$ = {lad['saturation_um']:.2f} µm", side="left")
    ref_line(ax, P.BETA_E5_UM, "$\\beta_0$ from E2–E5", side="right",
             style=(0, (4, 2)))
    ref_line(ax, lad["s90_um"], f"$s_{{90}}$ = {lad['s90_um']:.1f} µm", axis="x",
             at=0.05)
    dminband(ax, d_min_um)
    ax.set(xscale="log", yscale="log", xlim=(0.025, 160), ylim=(0.0015, 6),
           xlabel=LABEL["stroke"], ylabel=LABEL["lost"])
    ax.legend(loc="lower right", handlelength=1.6)
    finish(ax)
    panel_letter(ax, "a")

    bx.plot(s_i, f_i, lw=0.8, color=C["L"], alpha=0.45, zorder=1)
    bx.errorbar(s_i, f_i, yerr=fe_i, fmt="o", ms=4, color=C["L"], capsize=2,
                elinewidth=0.8, zorder=3)
    bx.errorbar(s_d * DODGE, f_d, yerr=fe_d, fmt="o", ms=4, mfc="white",
                mec=C["L"], ecolor=C["L"], capsize=2, elinewidth=0.8, zorder=3)
    bx.plot(s_c, f_c, "D", ms=3.2, color=C["ref"], ls="none", zorder=2)
    for value, text in ((0.0, "nothing"), (1.0, "all")):
        bx.axhline(value, color=C["ref"], lw=0.7, ls=(0, (1, 2)))
        bx.text(0.985, value, text, transform=bx.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8, color=C["ref"])
    bx.set(xscale="log", xlim=(0.025, 160), ylim=(-0.08, 1.12),
           xlabel=LABEL["stroke"], ylabel=LABEL["fraction"])
    finish(bx)
    panel_letter(bx, "b")

    caption = (
        "**F0b-2.** Lost motion after a reversal against stroke, from the "
        "reversal ladder at $x_0$ = 60 mm. (a) Lost motion as commanded "
        "stroke minus measured travel between the settled window means at "
        "the two turning points; filled: increasing pass, open: decreasing "
        "pass (dodged right), grey diamonds: the first half-cycle of each "
        "stroke, which starts from whatever the previous stroke left rather "
        "than from a closed loop. "
        "Error bars are ±2σ over half-cycles. Small strokes lose almost all "
        "their motion (dashed diagonal); lost motion saturates at "
        f"$\\ell_{{sat}}$ = {lad['saturation_um']:.2f} µm, about 30 % above "
        f"$\\beta_0$ = {P.BETA_E5_UM:g} µm from E2–E5, with a pre-sliding "
        f"length $s_{{90}}$ = {lad['s90_um']:.1f} µm — two and a half times "
        "the 3 µm that F1 and F2 were planned on. The two passes agree to "
        f"{lad['pass_difference_um'] * 1000:.0f} nm, so the ladder carries no "
        "history beyond the stroke itself and F1 needs no extra wipes between "
        "strokes. (b) The same data as the fraction of the commanded stroke "
        "that reaches the carriage; the light line is a reading aid. Strokes "
        "above 100 µm are not shown: their rungs run at over 25 µm/s and do "
        "not settle within the 0.5 s turning-point pause.")
    return save(fig, OUT, "f0b_2_lost_motion"), caption


# -- F0b-3 ----------------------------------------------------------------
def figure_slow_ramps(d: dict) -> tuple[list[Path], str]:
    run, ramps, grid = d["run"], d["ramps"], d["grid"]
    s = d["summary"]
    d_min_um = d["noise"]["d_min_nm"] / 1000
    fig = plt.figure(figsize=(FULL_WIDTH, 4.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.62, wspace=0.26)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    cx = fig.add_subplot(gs[1, :])

    for axis, (rev, cont), delta, dead, letter in (
            (ax, ("S1", "S4"), d["delta_down"], s["deadband_down_um"], "a"),
            (bx, ("S2", "S3"), d["delta_up"], s["deadband_up_um"], "b")):
        axis.plot([0, 20], [0, -20], ls=(0, (4, 2)), lw=0.7, color=C["ref"])
        axis.annotate("carriage stuck", xy=(3.5, -3.5), xytext=(4.7, -3.3),
                      fontsize=8, color=C["ref"], ha="left", va="center",
                      annotation_clip=True, zorder=6,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=C["ref"]))
        for r in ramps:
            if r.name not in (rev, cont):
                continue
            t, e = P.ramp_curve(run, r)
            axis.plot(t, e, lw=0.9, color=C["S"], alpha=0.85,
                      ls="-" if r.name == rev else (0, (3, 2)))
        axis.plot(grid, delta, lw=1.8, color=C["delta"])
        axis.plot([], [], lw=0.9, color=C["S"], label=f"reversal ({rev})")
        axis.plot([], [], lw=0.9, color=C["S"], ls=(0, (3, 2)),
                  label=f"continuation ({cont})")
        axis.plot([], [], lw=1.8, color=C["delta"], label="Δ")
        if np.isfinite(dead):
            axis.plot([dead, dead], [-0.55, 0.35], lw=0.9, color=C["ref"])
            axis.annotate(f"deadband {dead * 1000:.0f} nm", xy=(dead, 0.35),
                          xytext=(2.2, 1.45), fontsize=8, color=C["ref"],
                          arrowprops=dict(arrowstyle="-", lw=0.6,
                                          color=C["ref"]))
        axis.set(xlim=(0, 20), ylim=(-4.2, 2.0), xlabel=LABEL["travel"],
                 ylabel=LABEL["lag"])
        axis.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
                    handlelength=1.6, columnspacing=1.1)
        finish(axis)
        panel_letter(axis, letter)

    st = d["steps"]
    r1 = next(r for r in ramps if r.name == "S1")
    r4 = next(r for r in ramps if r.name == "S4")
    n = 30
    for r, ls in ((r1, "-"), (r4, (0, (4, 2)))):
        t = np.arange(0, int(n * st["dt_ms"])) / 1000.0
        x = run.x_nm[r.t0_ms:r.t0_ms + len(t)]
        x = (x - np.mean(x[:60])) * r.direction
        cx.plot(t, x, lw=0.5, color=C["S"], alpha=0.35, ls=ls)
        edges = np.arange(n + 1) * st["dt_ms"]
        means = [np.mean(x[int(a) + int(0.3 * st["dt_ms"]):int(b)])
                 for a, b in zip(edges[:-1], edges[1:])]
        cx.plot((edges[:-1] + st["dt_ms"] / 2) / 1000, means, "o", ms=3,
                mfc="white" if ls != "-" else C["S"], mec=C["S"], mew=0.8,
                label=f"{r.name} {'reversal' if r.reversal else 'continuation'}")
    stair_t = np.repeat(np.arange(n + 1) * st["dt_ms"], 2)[1:-1] / 1000
    stair_x = np.repeat(np.arange(n) * P.USTEP_NM, 2)
    cx.plot(stair_t, stair_x, lw=0.8, color=C["ref"], label="commanded, 39 nm steps")
    cx.set(xlim=(0, n * st["dt_ms"] / 1000), xlabel="time since the ramp start (s)",
           ylabel="carriage position (nm)")
    cx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
              handlelength=1.8)
    finish(cx)
    panel_letter(cx, "c", dx=-0.05)

    caption = (
        "**F0b-3.** Slow ramps at 0.25 µm/s after a reversal and in "
        "continuation over the same commanded positions, both repetitions of "
        "each. (a) Downward: reversal (S1, solid) and continuation (S4, "
        "dashed); their difference Δ (thick) isolates the friction transient, "
        "since the geometric and microstep error are identical on both. The "
        f"carriage leaves the stuck line after {s['deadband_down_um'] * 1000:.0f} nm, "
        "about one microstep, but only half of the commanded motion is "
        f"arriving by {s['half_arrival_down_um']:.1f} µm and Δ has not "
        "flattened within 20 µm — consistent with the 7.5 µm pre-sliding "
        f"length of F0b-2. (b) The same upward, deadband "
        f"{s['deadband_up_um'] * 1000:.0f} nm. (c) The first 30 microsteps "
        "against time: the commanded staircase advances 39 nm every 156 ms, "
        f"and the carriage follows it smoothly at a mean "
        f"{s['microstep_mean_increment_nm']:.0f} nm per step after a reversal, "
        f"with {s['microstep_stuck_fraction'] * 100:.0f} % of steps delivering "
        f"less than d_min and {s['microstep_jump_fraction'] * 100:.0f} % more "
        "than twice the mean. There is "
        f"{'a sawtooth' if s['stick_slip'] else 'no sawtooth'}: the carriage "
        "does not stick and jump, it simply arrives short.")
    return save(fig, OUT, "f0b_3_slow_ramps"), caption


# -- F0b-4 ----------------------------------------------------------------
def figure_dwell_speed(d: dict) -> tuple[list[Path], str]:
    dw, sp = d["dwell"], d["speed"]
    d_min_um = d["noise"]["d_min_nm"] / 1000
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(FULL_WIDTH, 2.7),
                                 gridspec_kw={"wspace": 0.34})

    short = [t for t in dw["trials"] if t["dwell_s"] < 10]
    long = [t for t in dw["trials"] if t["dwell_s"] >= 10]
    for group, ls, name in ((short, "-", f"{dw['short_s']:g} s at rest"),
                            (long, (0, (4, 2)), f"{dw['long_s']:g} s at rest")):
        for k, t in enumerate(group):
            ax.plot(t["travel_um"], t["lost_um"], ls=ls, lw=1.1, color=C["D"],
                    label=name if k == 0 else None)
    ref = np.mean([np.interp(np.linspace(0, 5, 200), t["travel_um"], t["lost_um"])
                   for t in short], axis=0)
    gx = np.linspace(0, 5, 200)
    ax.fill_between(gx, ref - d_min_um, ref + d_min_um, color=C["ref"],
                    alpha=0.25, lw=0)
    ax.set(xlim=(0, 5), ylim=(-0.2, 2.2), xlabel="commanded travel since the reversal (µm)",
           ylabel=LABEL["lost"])
    ax.legend(loc="lower right", handlelength=1.8)
    finish(ax)
    panel_letter(ax, "a")

    creep = next(t["creep"] for t in long if t["creep"])
    ins = ax.inset_axes([0.10, 0.56, 0.38, 0.38])
    tt = np.array(creep["t_s"])
    ins.plot(tt[tt > 0.05], np.array(creep["x_nm"])[tt > 0.05], lw=0.7,
             color=C["D"])
    ins.set(xscale="log", xlabel="time at rest (s)", ylabel="creep (nm)")
    ins.tick_params(labelsize=7)
    ins.xaxis.label.set_size(7)
    ins.yaxis.label.set_size(7)
    finish(ins)

    slow = sp["by_rate"][f"{sp['slow_um_s']:g}"]
    fast = sp["by_rate"][f"{sp['fast_um_s']:g}"]
    pos = np.array(slow["position_um"]) / 1000
    bx.fill_between(pos, -fast["timing_uncertainty_um"],
                    fast["timing_uncertainty_um"], color=C["ref"], alpha=0.22,
                    lw=0, label=f"{sp['fast_um_s'] / 1000:g} mm/s: timing "
                                "uncertainty alone")
    bx.plot(pos, slow["lag_um"], ls=(0, (4, 2)), lw=1.4, color=C["V"],
            label=f"{sp['slow_um_s'] / 1000:g} mm/s")
    ref_line(bx, P.BETA_E5_UM / 2, "$\\beta_0/2$", side="right")
    bx.set(xlabel="position within the segment (mm)",
           ylabel="friction lag (µm)", ylim=(-30, 30))
    bx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), handlelength=1.8)
    finish(bx)
    panel_letter(bx, "b")

    caption = (
        "**F0b-4.** (a) Lost motion after a reversal following "
        f"{dw['short_s']:g} s and {dw['long_s']:g} s at rest, two trials each "
        "in ABBA order; the band is ±d_min around the mean of the short-dwell "
        f"trials. The difference at the end of the ramp is "
        f"{dw['delta_nm']:+.0f} nm against a trial-to-trial scatter of "
        f"{dw['within_sd_nm']:.0f} nm, so although it exceeds the "
        f"{dw['d_min_nm']:.1f} nm detection limit it is not separable from the "
        "scatter on two repetitions: the probe is underpowered rather than "
        "null, and F2 is not decided by it. Inset: creep during the 60 s rest. "
        "(b) Friction lag during steady motion at "
        f"{sp['slow_um_s'] / 1000:g} mm/s, taken as half the difference of the "
        "up and down passes at the same commanded position, which cancels "
        f"the geometric error exactly: {sp['slow_lag_um']:.2f} µm, flat along "
        f"the segment to ±{sp['slow_lag_sd_um']:.2f} µm and close to "
        f"$\beta_0/2$. The {sp['fast_um_s'] / 1000:g} mm/s passes are not "
        "comparable with it and are not drawn: reading a lag part-way "
        "through a move needs the commanded position at that instant, the "
        "segment log gives only start and end times, and the end timestamp "
        "carries up to 26 ms of overhead beyond the step train. The grey "
        "band is what that uncertainty alone is worth at 1 mm/s — tens of "
        "microns, larger than any friction effect. F0b therefore does not "
        "settle the speed question; F3 needs the step-train end logged, or "
        "the two speeds compared between settled endpoints rather than "
        "mid-move.")
    return save(fig, OUT, "f0b_4_dwell_speed"), caption


# -- F0b-A1 ---------------------------------------------------------------
def figure_validity(d: dict) -> tuple[list[Path], str]:
    val, anchors, tracks = d["validity"], d["anchors"], d["tracks"]
    fig, axes = plt.subplots(2, 2, figsize=(FULL_WIDTH, 4.3),
                             gridspec_kw={"hspace": 0.62, "wspace": 0.30})
    (ax, bx), (cx, dx) = axes

    for mm, (t, y) in sorted(tracks.items()):
        ax.plot(t, y / 1000, marker="o", ms=4, lw=1.0, color=C["R"],
                mfc=C["R"] if mm < 60 else "white", mec=C["R"],
                ls="-" if mm < 60 else (0, (4, 2)), label=f"{mm:.0f} mm")
    ax.set(xlabel="session time (min)", ylabel="drift of the pair mean (µm)")
    ax.legend(loc="upper left", handlelength=1.8)
    finish(ax)
    panel_letter(ax, "a")

    sds = np.array(val["window_sd_nm"])
    bx.hist(sds, bins=np.linspace(0, 60, 25), color=C["N"], alpha=0.8)
    bx.axvline(val["threshold_nm"], color=C["ref"], lw=0.8, ls=(0, (4, 2)))
    bx.text(val["threshold_nm"] - 1.5, 0.95, "rejection threshold",
            transform=bx.get_xaxis_transform(), ha="right", va="top",
            fontsize=8, color=C["ref"], rotation=90)
    bx.set(xlabel="standard deviation of a record window (nm)",
           ylabel="windows")
    finish(bx)
    panel_letter(bx, "b")

    t = [v[0] for v in val["mscnt_returns"]]
    m = [v[1] for v in val["mscnt_returns"]]
    cx.plot(t, m, "o", ms=3, color=C["L"])
    cx.set(xlabel="session time (min)", ylabel="MSCNT at $x_0$",
           ylim=(min(m) - 8, max(m) + 8))
    finish(cx)
    panel_letter(cx, "c")

    dx.axis("off")
    lines = [f"{'check':<34}{'value':<26}{'':<4}"]
    for c in val["checks"]:
        mark = "pass" if c["pass"] else "FAIL"
        lines.append(f"{c['name']:<34}{c['value']:<26}{mark:<4}")
    dx.text(0.0, 1.0, "\n".join(lines), transform=dx.transAxes, va="top",
            ha="left", family="monospace", fontsize=7, color=C["ink"],
            linespacing=1.55)
    panel_letter(dx, "d", dx=-0.02)

    caption = (
        "**F0b-A1.** Session validity. (a) Drift of the reference-pair mean "
        "at 15 mm (filled, solid) and 105 mm (open, dashed), each referred to "
        "its own first visit; averaging the up and down approaches removes "
        "the reversal error and leaves the slow change. The two ends move "
        "differently — 15 mm walks "
        f"{max(abs(y[-1]) for _t, y in [tracks[min(tracks)]]) / 1000:.1f} µm "
        "over the session while 105 mm stays put — so the drift is a change "
        "in the stage, not in the instrument, and the correction is "
        "interpolated in position as well as in time. (b) Standard deviation "
        "of each 0.5 s record window; all are below the threshold once the "
        "50 ms trailing guard is applied, without which three of the twelve "
        "reference windows reach 2.1–2.9 µm because the two clocks have "
        "slipped far enough by the end of the session for the next move to "
        "start inside the window. (c) MSCNT at every return to $x_0$: one "
        "value throughout, so no steps were lost. (d) Checks.")
    return save(fig, OUT, "f0b_A1_validity"), caption


# -- T1 -------------------------------------------------------------------
def scales_card(d: dict) -> Path:
    s, lad = d["summary"], d["ladder"]
    rows = [
        ("Precision of a 0.5 s window mean", "σ", f"{s['sigma_500ms_nm']:.2f} nm", "F0b-1"),
        ("Minimum detectable step", "d_min", f"{s['d_min_nm']:.1f} nm", "F0b-1"),
        ("Saturation lost motion", "ℓ_sat",
         f"{s['saturation_um']:.2f} µm (β₀ = {P.BETA_E5_UM:g} µm)", "F0b-2"),
        ("Pre-sliding length", "s₉₀", f"{s['s90_um']:.1f} µm", "F0b-2"),
        ("Deadband after reversal, down / up", "—",
         f"{s['deadband_down_um'] * 1000:.0f} / {s['deadband_up_um'] * 1000:.0f} nm",
         "F0b-3"),
        ("Smallest stroke that moves after reversal", "—",
         f"{s['smallest_moving_stroke_usteps']:.0f} µstep "
         f"({s['smallest_moving_stroke_um'] * 1000:.0f} nm)", "F0b-2/3"),
        ("Stick-slip in continuation", "—",
         f"{'yes' if s['stick_slip'] else 'no'} "
         f"({s['microstep_stuck_fraction'] * 100:.0f} % of steps below d_min, "
         f"{s['microstep_jump_fraction'] * 100:.0f} % over twice the mean)",
         "F0b-3c"),
        ("Dwell effect (60 s vs 1 s)", "—",
         f"{s['dwell_delta_nm']:+.0f} ± {s['dwell_within_sd_nm']:.0f} nm "
         f"({'significant' if s['dwell_significant'] else 'not separable'})",
         "F0b-4a"),
        ("Sliding lag at 0.01 mm/s", "—",
         f"{s['speed_slow_lag_um']:.2f} µm (β₀/2 = {P.BETA_E5_UM / 2:.2f} µm)",
         "F0b-4b"),
        ("Speed effect (0.01 vs 1 mm/s)", "—",
         "not resolvable from this logging", "F0b-4b"),
    ]
    text = ["# T1 — F0b scales card", "",
            "The numbers that set F1–F3, generated from "
            "`f0b/f0b_summary.json`.", "",
            "| Quantity | Symbol | Value | From |", "|---|---|---|---|"]
    text += [f"| {a} | {b} | {c} | {e} |" for a, b, c, e in rows]
    path = P.OUT / "f0b_scales_card.md"
    path.write_text("\n".join(text) + "\n", encoding="utf-8")
    return path


def main() -> int:
    apply()
    d = P.process()
    captions = []
    for build in (figure_noise, figure_lost_motion, figure_slow_ramps,
                  figure_dwell_speed, figure_validity):
        paths, caption = build(d)
        captions.append(caption)
        print(f"  {paths[0].stem}: " + ", ".join(p.suffix[1:] for p in paths))
    card = scales_card(d)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "captions.md").write_text(
        "# F0b results figures — captions\n\n" + "\n\n".join(captions) + "\n",
        encoding="utf-8")
    print(f"  {card.relative_to(HERE)}")
    print(f"  {(OUT / 'captions.md').relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
