#!/usr/bin/env python3
"""The F0b results figures.

Two figures carry the findings and one appendix figure carries the session
validity, in the shared chapter style of the rendering instructions.

    python3 render_f0b_results.py

writes, as vector PDF and 300 dpi PNG, into rendered/f0b/

    f0b_1_reversal     what a reversal costs, and over what distance
    f0b_2_decisions    what the session resolves, and what it changes
    f0b_A1_validity    drift, windows, MSCNT, checks (appendix)

with the captions in rendered/f0b/captions.md and the scales card in
f0b/f0b_scales_card.md.

This replaces an earlier set of five figures, one per block. The two
slow-ramp panels became one, drawn as lost motion rather than lag so the
ladder points can be laid straight over them and the two independent
measurements of the same curve can be compared by eye; the noise floor,
the microstep detail, the dwell probe and the comparison with what F1 and
F2 were planned on now sit together as the decisions figure.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import f0b_process as P
from thesis_style import (C, FULL_WIDTH, LABEL, apply, finish, panel_letter,
                          ref_line, save)

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "f0b"
DODGE = 1.015            # 1.5 % of the axis, so error bars do not overlap


def ladder_series(halves, cap_um, pass_name=None, first=False):
    """Mean and 2 sigma of the reversal half-cycles at each stroke.

    first=True selects the first half-cycle of each stroke instead, which
    starts from whatever the previous stroke left rather than from a closed
    loop and is therefore drawn apart from the rest.
    """
    by: dict[float, list] = {}
    for h in halves:
        if first:
            if h.index != 1 or h.pass_name != "increasing":
                continue
        elif h.kind != "reversal" or h.index < 2 or h.pass_name != pass_name:
            continue
        by.setdefault(round(h.stroke_um, 4), []).append(h)
    s = np.array(sorted(k for k in by if k <= cap_um))
    lost = np.array([np.mean([q.lost_travel_um for q in by[k]]) for k in s])
    lost_e = np.array([2 * np.std([q.lost_travel_um for q in by[k]]) for k in s])
    frac = np.array([np.mean([q.tau for q in by[k]]) for k in s])
    frac_e = np.array([2 * np.std([q.tau for q in by[k]]) for k in s])
    return s, lost, lost_e, frac, frac_e


# -- figure 1: what a reversal costs --------------------------------------
def figure_reversal(d: dict) -> tuple[list[Path], str]:
    lad, halves, s = d["ladder"], d["halves"], d["summary"]
    d_min_um = d["noise"]["d_min_nm"] / 1000
    cap = lad["usable_above_um"]

    fig = plt.figure(figsize=(FULL_WIDTH, 4.6))
    gs = fig.add_gridspec(2, 2, hspace=0.78, wspace=0.27)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    cx = fig.add_subplot(gs[1, :])

    s_i, l_i, e_i, f_i, fe_i = ladder_series(halves, cap, "increasing")
    s_d, l_d, e_d, f_d, fe_d = ladder_series(halves, cap, "decreasing")
    s_c, l_c, _, f_c, _ = ladder_series(halves, cap, first=True)

    # (a) lost motion against stroke
    ax.plot([0.02, 5], [0.02, 5], ls=(0, (4, 2)), lw=0.7, color=C["ref"])
    ax.text(0.052, 0.028, "nothing arrives", fontsize=8, color=C["ref"],
            rotation=34)
    ax.plot(s_c, l_c, "D", ms=3.2, ls="none", color=C["ref"],
            label="first half-cycle")
    ax.errorbar(s_i, l_i, yerr=e_i, fmt="o", ms=4, color=C["L"], capsize=2,
                elinewidth=0.8, label="increasing")
    ax.errorbar(s_d * DODGE, l_d, yerr=e_d, fmt="o", ms=4, mfc="white",
                mec=C["L"], ecolor=C["L"], capsize=2, elinewidth=0.8,
                label="decreasing")
    ref_line(ax, lad["saturation_um"],
             f"$\\ell_{{sat}}$ = {lad['saturation_um']:.2f} µm", side="left")
    ref_line(ax, P.BETA_E5_UM, "$\\beta_0$ (E2–E5)", side="right",
             style=(0, (4, 2)))
    ref_line(ax, lad["s90_um"], f"$s_{{90}}$ = {lad['s90_um']:.1f} µm",
             axis="x", at=0.05)
    ax.axhspan(-d_min_um, d_min_um, color=C["ref"], alpha=0.18, lw=0, zorder=0)
    ax.text(0.012, d_min_um, "below detection",
            transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=8, color=C["ref"])
    ax.set(xscale="log", yscale="log", xlim=(0.025, 160), ylim=(0.0015, 6),
           xlabel=LABEL["stroke"], ylabel=LABEL["lost"])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
              handlelength=1.3, columnspacing=0.9, handletextpad=0.4)
    finish(ax)
    panel_letter(ax, "a")

    # (b) the same, as the question a reader actually asks
    bx.plot(s_i, f_i, lw=0.8, color=C["L"], alpha=0.45, zorder=1)
    bx.errorbar(s_i, f_i, yerr=fe_i, fmt="o", ms=4, color=C["L"], capsize=2,
                elinewidth=0.8, zorder=3)
    bx.errorbar(s_d * DODGE, f_d, yerr=fe_d, fmt="o", ms=4, mfc="white",
                mec=C["L"], ecolor=C["L"], capsize=2, elinewidth=0.8, zorder=3)
    bx.plot(s_c, f_c, "D", ms=3.2, ls="none", color=C["ref"], zorder=2)
    for value, text in ((0.0, "nothing"), (1.0, "all")):
        bx.axhline(value, color=C["ref"], lw=0.7, ls=(0, (1, 2)))
        bx.text(0.985, value, text, transform=bx.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8, color=C["ref"])
    bx.axvline(lad["s90_um"], color=C["ref"], lw=0.7, ls=(0, (1, 2)))
    bx.text(lad["s90_um"] * 1.15, 0.04, "$s_{90}$", fontsize=8, color=C["ref"])
    bx.set(xscale="log", xlim=(0.025, 160), ylim=(-0.08, 1.12),
           xlabel=LABEL["stroke"], ylabel=LABEL["fraction"])
    finish(bx)
    panel_letter(bx, "b")

    # (c) the same quantity again, continuous, from the slow ramps
    grid = d["grid"]
    cx.plot([0, 20], [0, 20], ls=(0, (4, 2)), lw=0.7, color=C["ref"])
    cx.annotate("carriage stuck", xy=(2.6, 2.6), xytext=(0.9, 3.5),
                fontsize=8, color=C["ref"], ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C["ref"]))
    cx.plot(grid, -d["rev_down"], lw=1.6, color=C["S"],
            label="after a reversal, down")
    cx.plot(grid, -d["rev_up"], lw=1.3, ls=(0, (4, 2)), color=C["S"],
            label="after a reversal, up")
    cx.plot(grid, -d["rev_down"] + d["delta_down"], lw=1.0, color=C["S"],
            alpha=0.45)
    cx.plot(grid, -d["rev_up"] + d["delta_up"], lw=1.0, ls=(0, (4, 2)),
            color=C["S"], alpha=0.45)
    cx.text(13.2, 0.26, "continuing", fontsize=8, color=C["S"], alpha=0.9)
    inside = (s_i <= 20) & (s_i >= 1)
    cx.errorbar(s_i[inside], l_i[inside], yerr=e_i[inside], fmt="o", ms=4.5,
                color=C["L"], capsize=2, elinewidth=0.8, zorder=4,
                label="ladder, reversal half-cycles")
    dead = s["deadband_down_um"]
    cx.plot([dead, dead], [0, 0.45], lw=0.9, color=C["ref"])
    cx.annotate(f"deadband {dead * 1000:.0f} nm, one microstep",
                xy=(dead, 0.45), xytext=(0.7, 1.35), fontsize=8,
                color=C["ref"],
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C["ref"]))
    ref_line(cx, lad["saturation_um"], "$\\ell_{sat}$", side="right")
    cx.set(xlim=(0, 20), ylim=(0, 3.9), xlabel=LABEL["travel"],
           ylabel=LABEL["lost"])
    cx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
              handlelength=1.8, columnspacing=1.4)
    finish(cx)
    panel_letter(cx, "c", dx=-0.05)

    caption = (
        f"**F0b-1. A reversal costs {lad['saturation_um']:.2f} µm, and takes "
        f"about {lad['s90_um']:.0f} µm of travel to pay off — two and a half "
        "times the pre-sliding length F1 and F2 were planned on.** "
        "(a) Lost motion against stroke from the reversal ladder at "
        "$x_0$ = 60 mm, as commanded stroke minus measured travel between the "
        "settled window means at the two turning points. Filled: increasing "
        "pass; open: decreasing pass, dodged right; grey diamonds: the first "
        "half-cycle of each stroke, which starts from whatever the previous "
        "stroke left rather than from a closed loop. Error bars are ±2σ over "
        "half-cycles. Small strokes lose almost all their motion (dashed "
        f"diagonal); the loss saturates at {lad['saturation_um']:.2f} µm, "
        f"about 30 % above $\\beta_0$ = {P.BETA_E5_UM:g} µm from E2–E5. The "
        f"two passes agree to {lad['pass_difference_um'] * 1000:.0f} nm, so "
        "the ladder carries no history beyond the stroke itself and F1 needs "
        "no extra wipes between strokes. Strokes above 100 µm are not shown: "
        "their rungs run at over 25 µm/s and do not settle within the 0.5 s "
        "turning-point pause. (b) The same data as the fraction of the "
        "commanded stroke that reaches the carriage; the light line is a "
        "reading aid. (c) The same quantity measured continuously instead, "
        "from the slow ramps at 0.25 µm/s after a 3 mm wipe (heavy) and "
        "continuing in the same direction (light); their difference is the "
        "part that a reversal adds, and it is taken over the same commanded "
        "positions so the geometric and microstep error cancel. The carriage "
        f"leaves the stuck line after {dead * 1000:.0f} nm, one microstep, "
        "and has still not saturated after 20 µm. The ladder points fall "
        "within the band the two ramp directions span over the whole "
        "1–20 µm range, two independent measurements of the same quantity "
        "from different blocks of the session; the difference between the "
        "up and down ramps, 0.3–0.5 µm, is as large as any difference "
        "between the methods, and is itself the direction asymmetry F1 "
        "checks at three positions.")
    return save(fig, OUT, "f0b_1_reversal"), caption


# -- figure 2: what the session resolves, and what it changes -------------
def figure_decisions(d: dict) -> tuple[list[Path], str]:
    nz, dw, sp, s = d["noise"], d["dwell"], d["speed"], d["summary"]
    d_min_um = nz["d_min_nm"] / 1000

    fig = plt.figure(figsize=(FULL_WIDTH, 4.8))
    gs = fig.add_gridspec(2, 2, hspace=0.85, wspace=0.32)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    cx, dx = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    # (a) the measurement floor
    styles = {"N-mid": ("-", 1.6), "N-end1": ((0, (4, 2)), 1.0),
              "N-end2": ((0, (1, 2)), 1.0)}
    for hold in nz["holds"]:
        tau = np.array(hold["taus_s"])
        dev = np.array([np.nan if v is None else v for v in hold["allan_nm"]])
        ls, lw = styles[hold["name"]]
        ax.plot(tau, dev, ls=ls, lw=lw, color=C["N"],
                label=f"{hold['position_mm']:.0f} mm")
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
    ax.text(0.985, nz["d_min_nm"], "$d_{\\min}$",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=8, color=C["ref"])
    ref_line(ax, P.USTEP_NM, "one microstep", side="right")
    ax.set(xscale="log", yscale="log", xlim=(1e-3, 10), ylim=(0.2, 90),
           xlabel=LABEL["tau"], ylabel=LABEL["precision"])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
              handlelength=1.6, columnspacing=1.0)
    finish(ax)
    panel_letter(ax, "a")

    # (b) the first 30 microsteps
    st, ramps, run = d["steps"], d["ramps"], d["run"]
    n = 30
    for name, ls, fill, label in (("S1", "-", C["S"], "after a reversal"),
                                  ("S4", (0, (4, 2)), "white", "continuing")):
        r = next(q for q in ramps if q.name == name)
        t = np.arange(0, int(n * st["dt_ms"])) / 1000.0
        x = run.x_nm[r.t0_ms:r.t0_ms + len(t)]
        x = (x - np.mean(x[:60])) * r.direction
        bx.plot(t, x, lw=0.5, color=C["S"], alpha=0.30, ls=ls)
        edges = np.arange(n + 1) * st["dt_ms"]
        means = [np.mean(x[int(a) + int(0.3 * st["dt_ms"]):int(b)])
                 for a, b in zip(edges[:-1], edges[1:])]
        bx.plot((edges[:-1] + st["dt_ms"] / 2) / 1000, means, "o", ms=3,
                mfc=fill, mec=C["S"], mew=0.8, ls="none", label=label)
    stair_t = np.repeat(np.arange(n + 1) * st["dt_ms"], 2)[1:-1] / 1000
    stair_x = np.repeat(np.arange(n) * P.USTEP_NM, 2)
    bx.plot(stair_t, stair_x, lw=0.8, color=C["ref"], label="commanded")
    bx.set(xlim=(0, n * st["dt_ms"] / 1000),
           xlabel="time since the ramp start (s)",
           ylabel="carriage position (nm)")
    bx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=3,
              handlelength=1.6, columnspacing=1.0)
    finish(bx)
    panel_letter(bx, "b")

    # (c) the dwell probe that was meant to decide F2
    short = [t for t in dw["trials"] if t["dwell_s"] < 10]
    long = [t for t in dw["trials"] if t["dwell_s"] >= 10]
    for group, ls, name in ((short, "-", f"{dw['short_s']:g} s at rest"),
                            (long, (0, (4, 2)), f"{dw['long_s']:g} s at rest")):
        for k, t in enumerate(group):
            cx.plot(t["travel_um"], t["lost_um"], ls=ls, lw=1.1, color=C["D"],
                    label=name if k == 0 else None)
    gx = np.linspace(0, 5, 200)
    ref = np.mean([np.interp(gx, t["travel_um"], t["lost_um"]) for t in short],
                  axis=0)
    cx.fill_between(gx, ref - d_min_um, ref + d_min_um, color=C["ref"],
                    alpha=0.3, lw=0)
    ins = cx.inset_axes([0.14, 0.60, 0.34, 0.34])
    creep = next(t["creep"] for t in long if t["creep"])
    tt = np.array(creep["t_s"])
    ins.plot(tt[tt > 0.05], np.array(creep["x_nm"])[tt > 0.05], lw=0.7,
             color=C["D"])
    ins.set(xscale="log", xlabel="time at rest (s)", ylabel="creep (nm)")
    ins.tick_params(labelsize=7)
    ins.xaxis.label.set_size(7)
    ins.yaxis.label.set_size(7)
    finish(ins)
    cx.set(xlim=(0, 5), ylim=(-0.2, 2.3),
           xlabel="commanded travel since the reversal (µm)",
           ylabel=LABEL["lost"])
    cx.legend(loc="lower right", handlelength=1.8)
    finish(cx)
    panel_letter(cx, "c")

    # (d) what the session changes about the plan
    rows = [
        ("$d_{\\min}$", 5.0, s["d_min_nm"], "nm", C["N"]),
        ("$s_{90}$", 3.0, s["s90_um"], "µm", C["L"]),
        ("$\\ell_{sat}$ vs $\\beta_0$", P.BETA_E5_UM, s["saturation_um"],
         "µm", C["L"]),
    ]
    dx.axvline(1.0, color=C["ref"], lw=0.8, ls=(0, (4, 2)))
    dx.text(1.0, 1.02, "as planned", transform=dx.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=8, color=C["ref"])
    for i, (name, planned, measured, unit, colour) in enumerate(rows):
        y = len(rows) - 1 - i
        ratio = measured / planned
        dx.plot(1.0, y, "o", ms=5, mfc="white", mec=colour, mew=1.0, zorder=3)
        dx.annotate("", xy=(ratio, y), xytext=(1.0, y),
                    arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.3,
                                    mutation_scale=10), zorder=2)
        dx.plot(ratio, y, "o", ms=6, color=colour, zorder=3)
        mid = (1.0 + ratio) / 2
        dx.text(mid, y + 0.17, f"{planned:.3g} → {measured:.3g} {unit}",
                ha="center", va="bottom", fontsize=8, color=colour)
        dx.text(mid, y - 0.17, f"×{ratio:.2f}", ha="center", va="top",
                fontsize=8, color=C["muted"])
    dx.set(xlim=(0, 3.2), ylim=(-0.65, len(rows) - 0.35),
           yticks=list(range(len(rows))),
           yticklabels=[r[0] for r in reversed(rows)],
           xlabel="measured / planned")
    finish(dx)
    dx.grid(axis="y", visible=False)
    panel_letter(dx, "d")

    caption = (
        f"**F0b-2. The session resolves {nz['d_min_nm']:.1f} nm, "
        f"{P.USTEP_NM / nz['d_min_nm']:.0f} times finer than a microstep, and "
        "moves two of the three numbers F1 and F2 were planned on.** "
        "(a) Precision of a window mean against averaging time at the middle "
        "and both ends of the travel, as the overlapping Allan deviation "
        "after removing a linear trend from each hold; the dashed line is the "
        "white-noise slope. At the 0.5 s record window (open markers) the "
        f"precision is {nz['sigma_500ms_nm']:.2f} nm, and the three positions "
        f"agree within {nz['position_spread_nm']:.2f} nm, so the floor does "
        "not depend on position. (b) The first 30 microsteps of a slow ramp: "
        "the commanded staircase advances 39 nm every 156 ms and the carriage "
        f"follows it smoothly, at a mean {s['microstep_mean_increment_nm']:.0f} "
        "nm per step after a reversal. No step delivers less than $d_{\\min}$ "
        "and none delivers more than twice the mean, so there is no "
        "stick-slip: the carriage simply arrives short, and single 1/256 "
        "steps are usable. (c) Lost motion after a reversal following "
        f"{dw['short_s']:g} s and {dw['long_s']:g} s at rest, two trials each "
        "in ABBA order; the band is ±$d_{\\min}$ about the mean of the "
        f"short-dwell trials. The difference is {dw['delta_nm']:+.0f} nm "
        f"against a trial-to-trial scatter of {dw['within_sd_nm']:.0f} nm, so "
        "although it exceeds the detection limit it is not separable from the "
        "scatter on two repetitions: the probe is underpowered rather than "
        "null, and F2 is not decided by it. Inset: creep during the 60 s "
        "rest. (d) Each F0b-dependent number against what F1 and F2 were "
        "planned on; open circle planned, filled measured. The speed probe "
        "is not drawn: reading a lag part-way through a move needs the "
        "commanded position at that instant, and the segment log gives only "
        "start and end times with up to 26 ms of overhead, worth tens of "
        f"microns at 1 mm/s. At {sp['slow_um_s'] / 1000:g} mm/s, where that "
        f"is sub-micron, the sliding lag is {sp['slow_lag_um']:.2f} µm "
        f"against $\\beta_0/2$ = {P.BETA_E5_UM / 2:.2f} µm.")
    return save(fig, OUT, "f0b_2_decisions"), caption


# -- appendix: session validity -------------------------------------------
def figure_validity(d: dict) -> tuple[list[Path], str]:
    val, tracks = d["validity"], d["tracks"]
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
    cx.plot(t, m, "o", ms=3, ls="none", color=C["L"])
    cx.set(xlabel="session time (min)", ylabel="MSCNT at $x_0$",
           ylim=(min(m) - 8, max(m) + 8))
    finish(cx)
    panel_letter(cx, "c")

    dx.axis("off")
    lines = [f"{'check':<24}{'value':<21}{'':<4}"]
    for c in val["checks"]:
        mark = "pass" if c["pass"] else "FAIL"
        lines.append(f"{c['name'][:23]:<24}{c['value'][:20]:<21}{mark:<4}")
    dx.text(0.0, 1.0, "\n".join(lines), transform=dx.transAxes, va="top",
            ha="left", family="monospace", fontsize=7, color=C["ink"],
            linespacing=1.7)
    panel_letter(dx, "d", dx=-0.02)

    caption = (
        "**F0b-A1. The session is valid: 346 of 346 segments, no lost steps, "
        "no driver faults.** (a) Drift of the reference-pair mean at 15 mm "
        "(filled, solid) and 105 mm (open, dashed), each referred to its own "
        "first visit; averaging the up and down approaches removes the "
        "reversal error and leaves the slow change. The two ends move "
        "differently, so the drift is a change in the stage rather than in "
        "the instrument, and the correction is interpolated in position as "
        "well as in time. (b) Standard deviation of each 0.5 s record window; "
        "all are below the threshold once the 50 ms trailing guard is "
        "applied, without which three of the twelve reference windows reach "
        "2.1–2.9 µm because the two free-running clocks have slipped far "
        "enough by the end of the session for the next move to start inside "
        "the window. (c) MSCNT at every return to $x_0$: one value "
        "throughout. (d) Checks.")
    return save(fig, OUT, "f0b_A1_validity"), caption


# -- T1 -------------------------------------------------------------------
def scales_card(d: dict) -> Path:
    s = d["summary"]
    rows = [
        ("Precision of a 0.5 s window mean", "σ",
         f"{s['sigma_500ms_nm']:.2f} nm", "F0b-2a"),
        ("Minimum detectable step", "d_min", f"{s['d_min_nm']:.1f} nm", "F0b-2a"),
        ("Saturation lost motion", "ℓ_sat",
         f"{s['saturation_um']:.2f} µm (β₀ = {P.BETA_E5_UM:g} µm)", "F0b-1a"),
        ("Pre-sliding length", "s₉₀", f"{s['s90_um']:.1f} µm", "F0b-1a"),
        ("Deadband after reversal, down / up", "—",
         f"{s['deadband_down_um'] * 1000:.0f} / "
         f"{s['deadband_up_um'] * 1000:.0f} nm", "F0b-1c"),
        ("Smallest stroke that moves after reversal", "—",
         f"{s['smallest_moving_stroke_usteps']:.0f} µstep "
         f"({s['smallest_moving_stroke_um'] * 1000:.0f} nm)", "F0b-1a"),
        ("Stick-slip in continuation", "—",
         f"{'yes' if s['stick_slip'] else 'no'} "
         f"({s['microstep_stuck_fraction'] * 100:.0f} % of steps below d_min, "
         f"{s['microstep_jump_fraction'] * 100:.0f} % over twice the mean)",
         "F0b-2b"),
        ("Dwell effect (60 s vs 1 s)", "—",
         f"{s['dwell_delta_nm']:+.0f} ± {s['dwell_within_sd_nm']:.0f} nm "
         f"({'significant' if s['dwell_significant'] else 'not separable'})",
         "F0b-2c"),
        ("Sliding lag at 0.01 mm/s", "—",
         f"{s['speed_slow_lag_um']:.2f} µm (β₀/2 = {P.BETA_E5_UM / 2:.2f} µm)",
         "F0b-2 caption"),
        ("Speed effect (0.01 vs 1 mm/s)", "—",
         "not resolvable from this logging", "F0b-2 caption"),
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
    for build in (figure_reversal, figure_decisions, figure_validity):
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
