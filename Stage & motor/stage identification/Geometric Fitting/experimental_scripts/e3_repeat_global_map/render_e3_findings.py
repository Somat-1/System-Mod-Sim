#!/usr/bin/env python3
"""Render the E3 figure set: F1-F6 and the appendix, per section 11.

Style and marker language follow the final E2 set, which section 11 asks for
and which matters more than any generic palette: these figures sit beside the
E1/E2/E4 ones. E2 is filled circles on a solid line, E3 open circles on a
dashed line, the pooled model a bold line with a band, and E4/E1/E0 keep their
dots, squares and diamonds.

The script refits what it draws rather than deserialising it. Every fit here is
a single least-squares solve - the expensive part of process_e3 is the cross
validation, and none of that is needed to draw a curve - so the figures cannot
drift out of step with the tables.
"""
from __future__ import annotations

import json
import math
import textwrap
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
E2_DIR = HERE.parent / "e2_global_map"
for directory in (E2_DIR, HERE):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

import e2_run                                   # noqa: E402
import e3_run                                   # noqa: E402
import render_e2_findings as e2f                 # noqa: E402
import process_e3 as p                           # noqa: E402

ANALYSIS = HERE / "analysis"
OUT = HERE / "rendered" / "experimental_findings"
PRESENTATION = HERE / "rendered" / "presentation"

E2_STYLE = dict(color=e2f.INK, marker="o", ms=2.6, ls="-", lw=1.0)
E3_STYLE = dict(color=e2f.BLUE, marker="o", ms=3.0, ls="--", lw=1.0,
                mfc="none", mew=.7)
POOLED = dict(color=e2f.GREEN, lw=1.8)
FLOOR = 0.24


CHAR_INCH = 0.0455          # mean advance of the caption face at 6.9 pt


def caption(fig: plt.Figure, left: float, text: str, width: int | None = None,
            bottom: float = .014, size: float = 6.9) -> int:
    """A wrapped footnote, anchored to the bottom left of the figure.

    matplotlib does not wrap fig.text, so a long caption runs off the canvas
    and takes its last clause with it - which is worse than truncating, since
    the figure still looks finished. The wrap width is derived from the
    figure's own width rather than fixed, because these figures run from 5.2
    to 7.5 inches and a count that suits one spills off another.

    Returns the line count, so a caller can check its bottom margin.
    """
    if width is None:
        usable = fig.get_figwidth() * (1.0 - left) - 0.12
        width = max(40, int(usable / CHAR_INCH))
    lines = textwrap.wrap(" ".join(text.split()), width=width)
    fig.text(left, bottom, "\n".join(lines), fontsize=size, color=e2f.INK2,
             va="bottom", linespacing=1.5)
    return len(lines)


def band(ax, x, centre, sigma, colour, alpha=.16, **kw):
    ax.fill_between(x, centre - 2 * sigma, centre + 2 * sigma, color=colour,
                    alpha=alpha, lw=0, **kw)


def hide_outside(ax, lo=6.0, hi=114.0):
    """Section 11: the fit is not drawn where it is not supported."""
    for x in (lo, hi):
        ax.axvline(x, color=e2f.LIGHT, lw=.6, ls=(0, (1, 2)), zorder=0)


def load_everything() -> dict:
    summary = json.loads((ANALYSIS / "e3_findings_summary.json").read_text())
    e2f.LOG = E2_DIR / "data" / "e2_terminal_20261027_151332.log"
    s2 = p.reduce_session(e2_run, "E2")
    e2f.LOG = e3_run.LOG
    s3 = p.reduce_session(e3_run, "E3")
    config = json.loads(p.CONFIG.read_text())
    flags = {"E2": p.window_quality(s2, config), "E3": p.window_quality(s3, config)}
    keep = flags["E2"]["keep"] & flags["E3"]["keep"]
    d2 = p.subset_data(s2["data"], keep)
    d3 = p.subset_data(s3["data"], keep)
    domain = p.basis_domain(d2)
    published = json.loads(p.E2_SUMMARY.read_text())
    row_b = published["model_selection"]["M_B_one_se"]
    row_c = published["model_selection"]["M_C_one_se"]
    models = {label: {"M_B": p.fit_with_settings(data, row_b, domain),
                      "M_C": p.fit_with_settings(data, row_c, domain)}
              for label, data in (("E2", d2), ("E3", d3))}
    chosen = summary["pooled"]["selected_model"]
    kc = None if chosen["K_C"] == "stationary" else int(chosen["K_C"])
    lam = float(chosen["lambda"]) if chosen["lambda"] not in ("", None) else None
    state = p.make_state(d2, int(chosen["H"]), int(chosen["K_s"]), kc, domain)
    pooled = p.fit_pooled([d2, d3], state, lam)
    state_b = p.make_state(d2, int(summary["pooled"]["selection"]["M_B_one_se"]["H"]),
                           int(summary["pooled"]["selection"]["M_B_one_se"]["K_s"]),
                           None, domain)
    pooled_b = p.fit_pooled([d2, d3], state_b, None)
    counts = d2["x"].size
    return {"summary": summary, "d2": d2, "d3": d3, "domain": domain,
            "models": models, "pooled": pooled, "pooled_b": pooled_b,
            "single": p.as_single_session(pooled),
            "single_b": p.as_single_session(pooled_b),
            "counts": counts, "sessions": (s2, s3), "config": config,
            "row_b": row_b, "row_c": row_c}


# ------------------------------------------------------------------ F1

def figure_f1(ctx: dict) -> None:
    d3, models, summary = ctx["d3"], ctx["models"], ctx["summary"]
    fig = plt.figure(figsize=(7.0, 5.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], left=.095,
                          right=.975, top=.915, bottom=.135, hspace=.44)
    ax = fig.add_subplot(gs[0])
    cross = p.cross_residual(models["E2"]["M_C"], d3)
    own = models["E3"]["M_C"]["residual"] - models["E3"]["M_C"]["residual"].mean()
    ax.axhline(0, color=e2f.GREY, lw=.55)
    ax.plot(d3["x"], cross, "o", ms=2.6, mfc="none", mew=.6, color=e2f.RUST,
            label="E3 after the E2 model (cross)")
    ax.plot(d3["x"], own, "o", ms=2.0, color=e2f.INK, alpha=.55,
            label="E3 after its own model")
    for values, colour, style in ((cross, e2f.RUST, "-"), (own, e2f.INK, "--")):
        profile = e2f.running_profile(d3["x"], values, 4.0)
        ax.plot(profile["x"], profile["mean"], style, color=colour, lw=1.4)
    e2f.endpoint_bands(ax)
    hide_outside(ax)
    ax.set(xlabel="position x [mm]", ylabel="residual [µm]")
    e2f.panel_title(ax, "a", "E3 residual, the E2 map applied without refitting",
                    raised=True)
    ax.legend(loc="lower left", frameon=False, ncols=2)
    e2f.finish_axis(ax)

    ax = fig.add_subplot(gs[1])
    rows = [row for row in _read_csv(ANALYSIS / "e3_t2_cross_prediction.csv")
            if row["direction"] == "E2 -> E3" and row["block"] != "overall"]
    labels = [row["block"].replace(" mm", "") for row in rows]
    series = [("rms_raw_um", "raw", e2f.LIGHT),
              ("rms_cross_M_B_um", "E2 $\\to$ E3, $M_B$", e2f.GREY),
              ("rms_cross_M_C_um", "E2 $\\to$ E3, $M_C$", e2f.BLUE),
              ("rms_own_M_C_um", "E3 own $M_C$", e2f.GREEN)]
    width = .8 / len(series)
    base = np.arange(len(rows))
    for j, (key, label, colour) in enumerate(series):
        ax.bar(base + j * width - .4 + width / 2,
               [float(row[key]) for row in rows], width * .92,
               color=colour, label=label, lw=0)
    ax.axhline(FLOOR, color=e2f.RUST, lw=1.0, ls=(0, (3, 2)),
               label=f"E4 local floor, {FLOOR} µm")
    ax.set_xticks(base, labels, fontsize=6.6)
    ax.set(xlabel="position block [mm]", ylabel="residual rms [µm]")
    e2f.panel_title(ax, "b", "Per 10 mm block", raised=True)
    ax.legend(loc="upper center", frameon=False, ncols=5, fontsize=6.4)
    e2f.finish_axis(ax)
    overall = summary["cross_prediction_overall"]["E2 -> E3"]
    caption(fig, .095,
             "F1. Same targets as E2, a different day. Cross "
             f"{overall['rms_cross_M_C_um']:.3f} µm against "
             f"{overall['rms_own_M_C_um']:.3f} µm for E3's own fit and "
             f"{overall['rms_raw_um']:.3f} µm uncorrected: "
             f"{overall['gap_closed_M_C']:.0%} of the gap to the noise floor. "
             "New positions are not tested by this session.")
    e2f.save(fig, "e3_f1_day_to_day", PRESENTATION)


# ------------------------------------------------------------------ F2

def figure_f2(ctx: dict) -> None:
    d2, d3 = ctx["d2"], ctx["d3"]
    grid = np.arange(6.0, 114.01, 1.0)
    sliding = {"E2": e2f.sliding_profile(d2), "E3": e2f.sliding_profile(d3)}
    windows1, windows4 = e2f.external_windows()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    fig.subplots_adjust(left=.085, right=.98, top=.9, bottom=.145,
                        hspace=.42, wspace=.26)
    pooled_profiles = {k: e2f.harmonic_profile(ctx["single"], grid, k)
                       for k in (1, 2)}
    own = {label: {k: e2f.harmonic_profile(ctx["models"][label]["M_C"], grid, k)
                   for k in (1, 2)} for label in ("E2", "E3")}
    for column, k in enumerate((1, 2)):
        ax = axes[0, column]
        profile = pooled_profiles[k]
        band(ax, grid, profile["amplitude"], profile["amplitude_sigma"],
             e2f.GREEN)
        ax.plot(grid, profile["amplitude"], label="pooled model", **POOLED)
        for label, style in (("E2", E2_STYLE), ("E3", E3_STYLE)):
            rows = sliding[label]
            ax.plot([row["x"] for row in rows], [row[f"A{k}"] for row in rows],
                    label=f"{label} sliding fit",
                    **{**style, "ms": style["ms"], "lw": .9})
        for source, windows, marker in (("E1", windows1, "s"),
                                        ("E4", windows4, ".")):
            centres = [w["centre"] for w in windows.values()]
            values = [abs(w["spectrum"][k - 1]["complex"]) for w in windows.values()]
            ax.plot(centres, values, marker, ms=5 if marker == "s" else 7,
                    mfc="none" if marker == "s" else e2f.RUST,
                    color=e2f.RUST, mew=.8, ls="none", label=source)
        e2f.endpoint_bands(ax, labels=column == 0)
        ax.set(xlabel="position x [mm]", ylabel=f"$A_{k}$ [µm]")
        e2f.panel_title(ax, "ab"[column], f"$k = {k}$ amplitude", raised=True)
        e2f.finish_axis(ax)
        if column == 0:
            ax.legend(loc="lower left", frameon=False, fontsize=6.2, ncols=2)

        ax = axes[1, column]
        for label, style in (("E2", E2_STYLE), ("E3", E3_STYLE)):
            profile = own[label][k]
            reference = e2f.phase_reference(profile)
            shift = -np.degrees(profile["phase"] - reference) / k
            sigma = np.degrees(profile["phase_sigma"]) / k
            band(ax, grid, shift, sigma,
                 e2f.INK if label == "E2" else e2f.BLUE, alpha=.1)
            ax.plot(grid, shift, label=f"{label} $M_C$",
                    **{**style, "marker": "none"})
        ax.axhline(0, color=e2f.GREY, lw=.55)
        e2f.endpoint_bands(ax, labels=False)
        ax.set(xlabel="position x [mm]",
               ylabel=rf"$\Delta\theta_{k}$ [deg]")
        e2f.panel_title(ax, "cd"[column], f"$k = {k}$ pattern shift", raised=True)
        # Where A_k is small the phase is barely determined, so the k = 2
        # curve runs to -65 deg at the support end and flattens everything
        # that matters. Clip to the range the band supports.
        inner = (grid >= 10.0) & (grid <= 110.0)
        limit = 1.9 * max(np.nanmax(np.abs(
            -np.degrees(own[label][k]["phase"][inner]
                        - e2f.phase_reference(own[label][k])) / k))
            for label in ("E2", "E3"))
        ax.set_ylim(-limit, limit)
        e2f.finish_axis(ax)
        if column == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=6.4)
    summary = ctx["summary"]["coefficient_comparison"]
    caption(fig, .085,
             "F2. The harmonic profile is a stage property, not a session "
             f"effect: max |ΔC₁| = {summary['C1']['max_complex_difference_um']:.3f} µm "
             f"and max |ΔC₂| = {summary['C2']['max_complex_difference_um']:.3f} µm "
             "between days, and the drift variants move the amplitudes by less "
             "than 0.3 nm.")
    e2f.save(fig, "e3_f2_harmonic_profile", PRESENTATION)


# ------------------------------------------------------------------ F3

def figure_f3(ctx: dict) -> None:
    grid = np.arange(6.0, 114.01, 1.0)
    fig = plt.figure(figsize=(6.9, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 1.0], left=.105,
                          right=.975, top=.915, bottom=.225, hspace=.3)
    ax = fig.add_subplot(gs[0])
    parts = {}
    for label, style in (("E2", E2_STYLE), ("E3", E3_STYLE)):
        model = ctx["models"][label]["M_C"]
        component = e2f.components(model, grid)
        parts[label] = component
        band(ax, grid, component["slow"], component["slow_sigma"],
             e2f.INK if label == "E2" else e2f.BLUE, alpha=.12)
        gain = model["beta"][1] * 1000.0
        ax.plot(grid, component["slow"],
                label=f"{label}  $s(x)$, gain {gain:+.1f} ppm",
                **{**style, "marker": "none"})
    pooled_component = e2f.components(ctx["single"], grid)
    ax.plot(grid, pooled_component["slow"],
            label=f"pooled, gain {ctx['single']['beta'][1] * 1000:+.1f} ppm",
            **POOLED)
    band(ax, grid, pooled_component["slow"], pooled_component["slow_sigma"],
         e2f.GREEN)
    ax.axhline(0, color=e2f.GREY, lw=.55)
    e2f.endpoint_bands(ax)
    ax.set(ylabel="slow error $s(x)$ [µm]")
    e2f.panel_title(ax, "a", "Slow error and gain, both sessions", raised=True)
    ax.legend(loc="lower left", frameon=False, fontsize=6.6)
    e2f.finish_axis(ax)

    ax = fig.add_subplot(gs[1])
    difference = parts["E3"]["slow"] - parts["E2"]["slow"]
    sigma = np.hypot(parts["E2"]["slow_sigma"], parts["E3"]["slow_sigma"])
    band(ax, grid, difference, sigma, e2f.RUST, alpha=.16)
    ax.plot(grid, difference, color=e2f.RUST, lw=1.3, label="E3 $-$ E2")
    # The correction D2 applied, in the same subspace s(x) lives in: its
    # linear part goes to the gain by construction and is removed here, so
    # the two curves are comparable rather than merely adjacent.
    correction = ctx["d3"]["D1"]["ybar"] - ctx["d3"]["D2"]["ybar"]
    x3 = ctx["d3"]["x"]
    q = np.column_stack((np.ones(x3.size), x3 - x3.mean()))
    coefficients, *_ = np.linalg.lstsq(q, correction, rcond=None)
    profile = p.running_mean(x3, correction - q @ coefficients, grid, 4.0)
    ax.plot(grid, profile, color=e2f.GREY, lw=1.1, ls=(0, (4, 2)),
            label="curvature of the D2 correction")
    ax.axhline(0, color=e2f.GREY, lw=.55)
    e2f.endpoint_bands(ax, labels=False)
    ax.set(xlabel="position x [mm]", ylabel="difference [µm]")
    e2f.panel_title(ax, "b", "Session difference against the drift it could be",
                    raised=True)
    ax.legend(loc="upper left", frameon=False, fontsize=6.4, ncols=2)
    e2f.finish_axis(ax)
    verdict = ctx["summary"]["coefficient_comparison"]["slow_vs_drift"]
    caption(fig, .105,
             f"F3. The D2 drift correction this session needed is "
             f"{verdict['correction_linear_ppm']:+.0f} ppm of gain and "
             f"{verdict['correction_curvature_range_um']:.3f} µm of everything "
             "else — linear in x exactly, because the up leg runs with time "
             "and the down leg against it. So the drift reaches the gain and "
             f"nothing else. The {verdict['gain_difference_ppm']:+.1f} ppm gain "
             f"difference is {verdict['gain_difference_over_correction']:.0%} "
             "of that correction and is therefore not evidence about the "
             f"screw; the {verdict['slow_difference_range_um']:.2f} µm s(x) "
             "difference in (b) is untouched by it, and is either a real "
             "day-to-day change or drift curvature within a segment, which "
             "three reference visits cannot separate.")
    e2f.save(fig, "e3_f3_slow_error", PRESENTATION)


# ------------------------------------------------------------------ F4

def figure_f4(ctx: dict) -> None:
    d2, d3, summary = ctx["d2"], ctx["d3"], ctx["summary"]
    structure = summary["residual_structure"]
    r2 = ctx["models"]["E2"]["M_C"]["residual"]
    r3 = ctx["models"]["E3"]["M_C"]["residual"]
    r2, r3 = r2 - r2.mean(), r3 - r3.mean()
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.3))
    fig.subplots_adjust(left=.075, right=.985, top=.855, bottom=.275, wspace=.33)

    ax = axes[0]
    grid = np.arange(6.0, 114.01, .5)
    for values, label, style in ((r2, "E2", E2_STYLE), (r3, "E3", E3_STYLE)):
        profile = p.running_mean(d2["x"], values, grid, 2.0)
        ax.plot(grid, profile, label=label, **{**style, "marker": "none"})
    ax.axhline(0, color=e2f.GREY, lw=.55)
    e2f.endpoint_bands(ax, labels=False)
    ax.set(xlabel="position x [mm]", ylabel="2 mm mean residual [µm]")
    e2f.panel_title(ax, "a", "Residual under each session's own fit",
                    raised=True)
    ax.legend(loc="upper left", frameon=False, fontsize=6.6)
    e2f.finish_axis(ax)

    ax = axes[1]
    ax.plot(r2, r3, "o", ms=2.8, mfc="none", mew=.6, color=e2f.BLUE)
    limit = 1.06 * max(np.abs(r2).max(), np.abs(r3).max())
    ax.plot([-limit, limit], [-limit, limit], color=e2f.GREY, lw=.7,
            ls=(0, (3, 2)))
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit),
           xlabel="E2 residual [µm]", ylabel="E3 residual [µm]")
    paired = structure["paired_own_models"]
    e2f.panel_title(ax, "b", f"Paired, ρ = {paired['rho']:+.2f}", raised=True)
    e2f.finish_axis(ax)

    ax = axes[2]
    for label, data, style in (("E2", d2, E2_STYLE), ("E3", d3, E3_STYLE)):
        up = data["D2"]["e_up"] - data["D2"]["e_up"].mean()
        down = data["D2"]["e_down"] - data["D2"]["e_down"].mean()
        ax.plot(up, down, style["marker"], ms=2.6,
                mfc="none" if label == "E3" else style["color"],
                mew=.6, color=style["color"], ls="none",
                label=f"{label}, ρ = {structure[f'up_down_{label}']['rho']:+.2f}")
    ax.set(xlabel="e↑ [µm]", ylabel="e↓ [µm]")
    e2f.panel_title(ax, "c", "Within session, up against down", raised=True)
    ax.legend(loc="upper left", frameon=False, fontsize=6.4)
    e2f.finish_axis(ax)

    running = structure["running_mean_2mm"]
    caption(fig, .075,
             "F4. The residual repeats between days at the same targets: "
             f"paired ρ = {paired['rho']:+.2f} over {paired['n']} targets "
             f"({running['rho']:+.2f} on the prescribed 2 mm running mean, "
             f"p = {running['p_circular_shift']:.3f}, circular-shift null). "
             f"That splits the residual into {paired['repeatable_rms_um']:.2f} µm "
             f"repeatable and {paired['independent_rms_um']:.2f} µm independent, "
             f"the latter near the {paired['e4_local_floor_um']} µm E4 floor. "
             "The repeatable half is geometric and compensable.")
    e2f.save(fig, "e3_f4_residual_structure", PRESENTATION)


# ------------------------------------------------------------------ F5

def figure_f5(ctx: dict) -> None:
    d2, d3, summary = ctx["d2"], ctx["d3"], ctx["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    fig.subplots_adjust(left=.085, right=.98, top=.855, bottom=.275, wspace=.26)

    ax = axes[0]
    for label, data, colour in (("E2", d2, e2f.INK), ("E3", d3, e2f.BLUE)):
        for variant in ("D0", "D1", "D2"):
            profile = e2f.running_profile(data["x"], data[variant]["B"], 6.0)
            ax.plot(profile["x"], profile["mean"],
                    color=e2f.DRIFT_COLOURS[variant],
                    lw=1.5 if label == "E2" else 1.1,
                    ls="-" if label == "E2" else "--")
    ax.axhline(0, color=e2f.GREY, lw=.55)
    e2f.endpoint_bands(ax, labels=False)
    ax.set(xlabel="position x [mm]", ylabel="reversal $B(x)$ [µm]")
    e2f.panel_title(ax, "a", "6 mm mean reversal, every drift variant",
                    raised=True)
    handles = [Line2D([], [], color=e2f.DRIFT_COLOURS[v], lw=1.4, label=v)
               for v in ("D0", "D1", "D2")]
    handles += [Line2D([], [], color=e2f.GREY, lw=1.4, ls="-", label="E2"),
                Line2D([], [], color=e2f.GREY, lw=1.1, ls="--", label="E3")]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=6.2,
              ncols=5)
    e2f.finish_axis(ax)

    ax = axes[1]
    leakage = summary["leakage"]
    for name, colour, style in (("E2", e2f.INK, "-"), ("E3", e2f.BLUE, "--")):
        block = leakage[name]
        x = np.asarray(block["difference_profile"]["x_mm"])
        difference = np.asarray(block["difference_profile"]["difference_um"])
        order = np.argsort(x)
        statistic = block["F"]
        ax.plot(x[order], difference[order], color=colour, lw=1.3, ls=style,
                label=f"{name}: $s_↑-s_↓$, F = {statistic['statistic']:.1f}")
    ax.axhline(0, color=e2f.GREY, lw=1.0, ls=(0, (3, 2)),
               label="constant difference (H$_0$)")
    e2f.endpoint_bands(ax, labels=False)
    ax.set(xlabel="position x [mm]", ylabel="direction difference [µm]")
    e2f.panel_title(ax, "b", "Split-direction leakage test", raised=True)
    ax.legend(loc="upper right", frameon=False, fontsize=6.0)
    e2f.finish_axis(ax)

    spread = summary["reversal"]["session_spread"]
    caption(fig, .085,
             "F5. β₀ over the sessions that can measure it: "
             + ", ".join(f"{name} {spread['values_um'][name]:+.2f}"
                         for name in spread["applied_sessions"])
             + f" µm, range {spread['range_um']:.2f} µm. E3 is excluded: its "
             f"three drift variants span {spread['E3_variant_spread_um']:.2f} µm, "
             "because B differences the directions and so carries the whole "
             "drift over the up-to-down gap.")
    e2f.save(fig, "e3_f5_reversal_leakage", PRESENTATION)


# ------------------------------------------------------------------ F6

def figure_f6(ctx: dict) -> None:
    summary = ctx["summary"]
    v2 = summary["validation"]["V2"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    fig.subplots_adjust(left=.085, right=.98, top=.855, bottom=.275, wspace=.28)

    ax = axes[0]
    centres = np.arange(len(v2))
    series = [("rms_raw_um", "raw", e2f.LIGHT),
              ("rms_M_B_um", "pooled $M_B$", e2f.GREY),
              ("rms_M_C_um", "pooled $M_C$", e2f.GREEN),
              ("rms_local_um", "local floor", e2f.RUST)]
    width = .8 / len(series)
    for j, (key, label, colour) in enumerate(series):
        ax.bar(centres + j * width - .4 + width / 2,
               [row[key] for row in v2], width * .9, color=colour, label=label,
               lw=0)
    ax.set_xticks(centres,
                  [f"{row['source']} {row['window']}\n{row['x_centre_mm']:.0f}"
                   for row in v2], fontsize=6.0)
    ax.set(xlabel="E1 / E4 window, centre [mm]", ylabel="residual rms [µm]")
    e2f.panel_title(ax, "a", "V2, independent windows", raised=True)
    ax.legend(loc="upper right", frameon=False, fontsize=6.2, ncols=2)
    e2f.finish_axis(ax)

    ax = axes[1]
    metrics = summary["validation"]["V3"]
    prediction = summary["validation"]["E5_prediction"]
    keys = ("E", "M", "A")
    base = np.arange(len(keys))
    width = .19
    for j, (session, colour, hatch) in enumerate((("v3", e2f.BLUE, None),
                                                  ("v4", e2f.GREEN, "///"))):
        ax.bar(base + (2 * j - 1.5) * width,
               [metrics[session]["raw"][key] for key in keys], width * .92,
               color=colour, alpha=.35, lw=0,
               label=f"E0 {session}, measured")
        values = [prediction[session][key]["mean_um"] for key in keys]
        errors = [2 * prediction[session][key]["sd_um"] for key in keys]
        ax.bar(base + (2 * j - .5) * width, values, width * .92, color=colour,
               lw=0, hatch=hatch, edgecolor="white", yerr=errors,
               error_kw=dict(ecolor=e2f.INK, lw=.8, capsize=2),
               label=f"E0 {session}, after G + B")
    ax.set_xticks(base, [f"{key} [µm]" for key in keys])
    ax.set(ylabel="ISO 230-2 value [µm]")
    e2f.panel_title(ax, "b", "V3, the E5 prediction", raised=True)
    ax.legend(loc="upper left", frameon=False, fontsize=6.2)
    e2f.finish_axis(ax)

    mean = summary["validation"]["V2_mean_rms_um"]
    caption(fig, .085,
             f"F6. Pooled model on E1/E4: {mean['raw']:.2f} µm raw to "
             f"{mean['M_C']:.2f} µm, against a {mean['local']:.2f} µm local "
             "floor. On E0 the predicted compensated values carry ±2σ from "
             f"{p.N_BOOTSTRAP} posterior draws of the model and of β₀; the "
             "spread between the v3 and v4 sessions is the larger term.")
    e2f.save(fig, "e3_f6_final_validation", PRESENTATION)


# ------------------------------------------------------------- appendix

def appendix_drift(ctx: dict) -> None:
    """The session's dominant defect, and why it spares the geometry."""
    d2, d3, summary = ctx["d2"], ctx["d3"], ctx["summary"]
    drift = summary["drift"]["E3"]
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.3))
    fig.subplots_adjust(left=.075, right=.985, top=.85, bottom=.30, wspace=.36)

    ax = axes[0]
    for label, data, style in (("E2", d2, E2_STYLE), ("E3", d3, E3_STYLE)):
        times = np.asarray(summary["drift"][label]["anchor_time_s"]) / 60.0
        errors = np.asarray(summary["drift"][label]["anchor_error_um"])
        ax.plot(times, errors, label=label, **{**style, "ms": 5.5})
    for limit in (-.5, .5):
        ax.axhline(limit, color=e2f.RUST, lw=.8, ls=(0, (3, 2)))
    ax.set(xlabel="time [min]", ylabel="reference error [µm]")
    e2f.panel_title(ax, "a", "Datum at $x_{ref}$, ±0.5 µm limit", raised=True)
    ax.legend(loc="lower left", frameon=False, fontsize=6.6)
    e2f.finish_axis(ax)

    # How long the drift has to act before it reaches each estimator. B sees
    # the whole up-to-down gap; ybar sees only the departure of the mean visit
    # time from its own mean, which the mirrored design holds at zero. Both
    # are durations, so they belong on one axis, and the contrast is the point.
    ax = axes[1]
    gap = (d3["down_time_s"] - d3["up_time_s"]) / 60.0
    mean_time = 0.5 * (d3["up_time_s"] + d3["down_time_s"]) / 60.0
    ax.plot(d3["x"], gap, color=e2f.RUST, lw=1.4,
            label="$B$: up-to-down gap")
    ax.plot(d3["x"], np.abs(mean_time - mean_time.mean()), color=e2f.INK,
            lw=1.4, ls=(0, (4, 2)),
            label="$\\bar y$: departure of the mean visit time")
    ax.annotate(f"{np.ptp(mean_time) * 60:.3f} s", (d3["x"][len(gap) // 2], 0),
                textcoords="offset points", xytext=(0, 7), fontsize=6.6,
                color=e2f.INK, ha="center")
    ax.set_ylim(bottom=-1.5)
    ax.set(xlabel="position x [mm]", ylabel="time the drift acts over [min]")
    e2f.panel_title(ax, "b", "Why $\\bar y$ survives and $B$ does not",
                    raised=True)
    ax.legend(loc="upper right", frameon=False, fontsize=6.2)
    e2f.finish_axis(ax)

    ax = axes[2]
    rows = [("$A_1$", "A1_mean_um"), ("$A_2$", "A2_mean_um"),
            ("$s(x)$ range", "slow_range_um"), ("gain, ppm/10", "gain_ppm"),
            ("$β_0$", "beta0_um")]
    base = np.arange(len(rows))
    floor = 1e-5
    for j, (label, session) in enumerate((("E2", "E2"), ("E3", "E3"))):
        spread = summary["drift_variant_robustness"][session]["spread"]
        values = [max(spread[key] / (10.0 if key == "gain_ppm" else 1.0), floor)
                  for _, key in rows]
        ax.barh(base + (j - .5) * .38, values, .34, left=floor,
                color=e2f.INK if j == 0 else e2f.BLUE, lw=0, label=label)
    ax.set_yticks(base, [label for label, _ in rows], fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(floor, 10.0)
    ax.axvline(0.01, color=e2f.RUST, lw=.9, ls=(0, (3, 2)))
    ax.text(0.012, -.46, "usable", fontsize=6.4, color=e2f.RUST, va="bottom")
    ax.set(xlabel="D0 / D1 / D2 spread [µm]")
    e2f.panel_title(ax, "c", "Drift-variant spread", raised=True)
    ax.legend(loc="lower right", frameon=False, fontsize=6.4)
    e2f.finish_axis(ax)

    caption(fig, .075,
             f"Appendix. The reference moved {drift['total_change_um']:+.2f} µm, "
             f"at {drift['segments'][-1]['rate_um_per_h']:+.0f} µm/h over the "
             "second half. The mirrored design holds each target's mean visit "
             f"time constant to {drift['mean_visit_time_spread_s']:.3f} s, so "
             "linear drift cancels from ȳ and only its curvature survives; B "
             "gets the full excursion. MSCNT is identical at all seven "
             "references, so no steps were lost.")
    e2f.save(fig, "e3_appendix_drift", OUT)


def appendix_validity(ctx: dict) -> None:
    summary = ctx["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    fig.subplots_adjust(left=.08, right=.98, top=.85, bottom=.265, wspace=.26)
    ax = axes[0]
    for label, data, style in (("E2", ctx["d2"], E2_STYLE),
                               ("E3", ctx["d3"], E3_STYLE)):
        counts = np.bincount(data["phase_bin"], minlength=50)[:50]
        ax.step(np.arange(50), counts, where="mid", label=label,
                color=style["color"], lw=1.1,
                ls="-" if label == "E2" else "--")
    ax.axhline(5, color=e2f.RUST, lw=.9, ls=(0, (3, 2)), label="minimum 5")
    ax.set(xlabel="screw-phase bin", ylabel="targets")
    e2f.panel_title(ax, "a", "Phase coverage", raised=True)
    ax.legend(loc="lower center", frameon=False, fontsize=6.4, ncols=3)
    e2f.finish_axis(ax)

    ax = axes[1]
    rows = summary["pooled"]["selection"]
    table = _read_csv(ANALYSIS / "e3_t4_pooled_model_selection.csv")
    for kind, colour, marker in (("M_B", e2f.GREY, "s"), ("M_C", e2f.BLUE, "o")):
        subset = [row for row in table if row["model"] == kind]
        ax.plot([float(row["effective_parameters"]) for row in subset],
                [float(row["CV_RMSE_um"]) for row in subset], marker, ms=2.4,
                mfc="none", mew=.5, color=colour, ls="none", label=kind)
        chosen = rows[f"{kind}_one_se"]
        ax.plot([float(chosen["effective_parameters"])],
                [float(chosen["CV_RMSE_um"])], "*", ms=10, color=e2f.GREEN,
                ls="none")
    ax.plot([], [], "*", ms=9, color=e2f.GREEN, ls="none",
            label="one-SE choice")
    ax.set(xlabel="effective parameters", ylabel="blocked CV rmse [µm]")
    e2f.panel_title(ax, "b", "Pooled model selection", raised=True)
    ax.legend(loc="upper right", frameon=False, fontsize=6.4)
    e2f.finish_axis(ax)
    caption(fig, .08,
             "Appendix. Phase coverage and the pooled selection surface. "
             "Leave-one-session-out on two sessions at identical targets is a "
             "replication test, not a test of new positions.")
    e2f.save(fig, "e3_appendix_validity_selection", OUT)


def appendix_homing(ctx: dict) -> None:
    summary = ctx["summary"]
    homing = summary["homing_check"]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    fig.subplots_adjust(left=.14, right=.97, top=.86, bottom=.30)
    base = np.arange(2)
    values = [homing["delta1_deg"], homing["delta2_deg"]]
    errors = [2 * homing["delta1_sigma_deg"], 2 * homing["delta2_sigma_deg"]]
    ax.bar(base, values, .30, color=e2f.BLUE, lw=0, yerr=errors,
           error_kw=dict(ecolor=e2f.INK, lw=.9, capsize=3))
    # The sign convention between MSCNT and the harmonic phase is not fixed by
    # the log, so the prediction is drawn as a magnitude, both ways.
    predicted = abs(homing["mscnt"]["screw_deg"])
    for sign in (+1, -1):
        ax.axhline(sign * predicted, color=e2f.RUST, lw=1.1, ls=(0, (3, 2)))
    for sign in (+1, -1):
        ax.axhline(sign * homing["requirement_deg"], color=e2f.GREY, lw=.8,
                   ls=(0, (1, 2)))
    ax.text(.27, -predicted - .14,
            f"driver MSCNT: ±{predicted:.2f}°", color=e2f.RUST, fontsize=6.6,
            va="top")
    ax.text(-.42, homing["requirement_deg"] - .12,
            f"requirement ±{homing['requirement_deg']:.0f}°", color=e2f.INK2,
            fontsize=6.6, va="top")
    ax.set_xlim(-.55, 1.55)
    ax.set_ylim(-4.3, 3.6)
    ax.axhline(0, color=e2f.GREY, lw=.6)
    ax.set_xticks(base, [r"$\delta_1$ from $k=1$", r"$\delta_2$ from $k=2$"])
    ax.set(ylabel="homing shift δ [deg]")
    ax.set_title("Screw-angle reference between sessions", loc="left", pad=4)
    e2f.finish_axis(ax)
    caption(fig, .14,
             f"Appendix. δ₁ = {homing['delta1_deg']:+.2f}° ± "
             f"{homing['delta1_sigma_deg']:.2f}° and δ₂ = "
             f"{homing['delta2_deg']:+.2f}° ± {homing['delta2_sigma_deg']:.2f}° "
             "agree, so the rotation is a reference offset and not a change in "
             f"the screw. The driver's MSCNT shifted {homing['mscnt']['difference_counts']} "
             f"counts between sessions, {homing['mscnt']['usteps']:.0f} µsteps "
             f"= {abs(predicted):.2f}°, an independent prediction of the same "
             "quantity. Both are inside the 3° requirement.")
    e2f.save(fig, "e3_appendix_homing", OUT)


def appendix_periodogram(ctx: dict) -> None:
    counts = ctx["counts"]
    pooled = ctx["pooled"]["residual"]
    r2, r3 = pooled[:counts], pooled[counts:]
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    fig.subplots_adjust(left=.11, right=.97, top=.86, bottom=.28)
    x = np.r_[ctx["d2"]["x"], ctx["d3"]["x"]]
    frequency, power = e2f.residual_periodogram(
        x, np.r_[r2 - r2.mean(), r3 - r3.mean()], 0.02, 0.5)
    ax.plot(frequency, power, color=e2f.BLUE, lw=.9)
    peak = ctx["summary"]["residual_structure"]["pooled_periodogram"]
    ax.plot([peak["peak_cycles_per_mm"]], [peak["peak_power"]], "v", ms=6,
            color=e2f.RUST)
    ax.annotate(f"{peak['peak_period_mm']:.2f} mm",
                (peak["peak_cycles_per_mm"], peak["peak_power"]),
                textcoords="offset points", xytext=(6, -2), fontsize=6.8,
                color=e2f.RUST)
    ax.set(xlabel="spatial frequency [cycles/mm]", ylabel="normalised power")
    ax.set_title("Pooled residual periodogram", loc="left", pad=4)
    e2f.finish_axis(ax)
    caption(fig, .11,
             "Appendix. 0.02-0.5 cycles/mm. A peak at the rail bolt pitch "
             "would point at the guide; a smooth red spectrum at drift or air.")
    e2f.save(fig, "e3_appendix_periodogram", OUT)


def _read_csv(path: Path) -> list[dict]:
    import csv
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in row.items():
            if key in ("direction", "block", "model", "K_C", "source",
                       "window", "session", "check", "criterion"):
                continue
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def main() -> int:
    e2f.style()
    OUT.mkdir(parents=True, exist_ok=True)
    PRESENTATION.mkdir(parents=True, exist_ok=True)
    ctx = load_everything()
    figure_f1(ctx)
    figure_f2(ctx)
    figure_f3(ctx)
    figure_f4(ctx)
    figure_f5(ctx)
    figure_f6(ctx)
    appendix_drift(ctx)
    appendix_validity(ctx)
    appendix_homing(ctx)
    appendix_periodogram(ctx)
    for directory in (PRESENTATION, OUT):
        for path in sorted(directory.glob("e3_*.png")):
            print("wrote", path.relative_to(HERE))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
