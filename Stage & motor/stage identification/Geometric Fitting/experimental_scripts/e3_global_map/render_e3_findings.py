#!/usr/bin/env python3
"""Figures for the 3 Nov E3 session.

The session's job is one question - does E2's model survive at positions it
never saw - so F2 is the figure the rest support. F1 establishes that the
capture is worth interpreting, F3 says which parts of the model transfer and
which do not, F4 does the same for reversal, and the appendix carries the
drift and validity evidence.
"""
from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
E2 = HERE.parent / "e2_global_map"
E5 = HERE.parent / "e5_compensated_validation"
for path in (E2, E5, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import render_e2_findings as e2f                       # noqa: E402
import e2_run                                          # noqa: E402
import e3_run                                          # noqa: E402
from e5_geometry import FrozenModel                    # noqa: E402
from process_e3 import reduce_session, M_B, M_C, own_fit  # noqa: E402

OUT = HERE / "rendered" / "findings"
SUMMARY = HERE / "analysis" / "e3_findings_summary.json"

C_E2, C_E3 = "#176d8c", "#b4473f"
C_MUTED, C_GRID = "#566069", "#dde1e4"
CHAR_INCH = 0.0455


def caption(fig, left, text, width=None, bottom=0.014, size=6.9) -> None:
    if width is None:
        usable = fig.get_figwidth() * (1.0 - left) - 0.12
        width = max(40, int(usable / CHAR_INCH))
    lines = textwrap.wrap(" ".join(text.split()), width=width)
    fig.text(left, bottom, "\n".join(lines), fontsize=size, color=C_MUTED,
             va="bottom", linespacing=1.5)


def header(fig, left, title, subtitle) -> None:
    """Title and subtitle, spaced in inches so short figures do not collide."""
    height = fig.get_figheight()
    top = 1.0 - 0.22 / height
    fig.text(left, top, title, fontsize=15, weight="bold", color="#20262b",
             va="top")
    fig.text(left, top - 0.30 / height, subtitle, fontsize=8.5, color=C_MUTED,
             va="top")


def tidy(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=C_GRID, lw=0.5)


# ---------------------------------------------------------------- figure 1

def figure_f1(d2: dict, d3: dict, summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 6.1))
    gs = fig.add_gridspec(2, 2, left=0.095, right=0.975, top=0.845,
                          bottom=0.165, hspace=0.46, wspace=0.26)
    header(fig, 0.095, "E3 — the session, against E2",
           "400 positions from a new seed, reduced by the same code as E2")

    ax = fig.add_subplot(gs[0, :])
    ax.plot(d2["x"], d2["D2"]["ybar"], lw=0.9, color=C_E2, label="E2 (seed 5315)")
    ax.plot(d3["x"], d3["D2"]["ybar"], lw=0.9, color=C_E3,
            label="E3 (seed 9101, new positions)")
    ax.set(xlabel="position [mm]", ylabel=r"$\bar{y}$  [µm]")
    ax.set_title("(a) Mean bidirectional error along the travel", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    tidy(ax)

    ax = fig.add_subplot(gs[1, 0])
    for data, colour, label in ((d2, C_E2, "E2"), (d3, C_E3, "E3")):
        t = data["reference_time_s"] / 60.0
        ax.plot(t, data["reference_error_um"] - data["reference_error_um"][0],
                marker="o", ms=3.4, lw=0.9, color=colour, label=label)
    ax.set(xlabel="session time [min]", ylabel="reference error [µm]")
    ax.set_title("(b) Drift at the reference", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)

    ax = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0, 0.6, 7) * 0.04 / 0.04
    for data, colour, label in ((d2, C_E2, "E2"), (d3, C_E3, "E3")):
        gaps = np.diff(data["usteps"]) / 64.0
        ax.hist(gaps * 0.04, bins=np.arange(0, 0.62, 0.04), histtype="step",
                lw=1.1, color=colour, label=label)
    ax.set(xlabel="step between consecutive targets [mm]", ylabel="count")
    ax.set_title("(c) The two target sets are independent", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)

    drift = summary["drift"]
    caption(fig, 0.095,
            f"The session drifted {summary['validity'][-1]['value']} µm from datum to end, "
            f"against 5.41 µm on 30 Oct, and the reference spread over the whole run was "
            f"{summary['validity'][-2]['value']} µm. Mirrored timing holds the mean visit time to "
            f"{drift['mean_visit_time_spread_s']:.3f} s across all 400 targets, so D0 and D1 agree "
            f"to {drift['D0_vs_D1_max_ybar_um'] * 1000:.3f} nm; the piecewise correction D2 moves "
            f"ȳ by at most {drift['D1_vs_D2_max_ybar_um']:.3f} µm. Every result below uses D2.")
    e2f.save(fig, "e3_f1_session", OUT)


# ---------------------------------------------------------------- figure 2

def figure_f2(frozen: dict, transfer: list[dict], summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 6.3))
    gs = fig.add_gridspec(2, 2, left=0.095, right=0.975, top=0.855,
                          bottom=0.170, hspace=0.48, wspace=0.28)
    header(fig, 0.095, "E3 — does E2's model predict positions it never saw?",
           "The model frozen for E5, applied unchanged to 400 new targets")

    x, measured, predicted = frozen["_x"], frozen["_measured"], frozen["_predicted"]
    residual = frozen["_residual"]

    ax = fig.add_subplot(gs[0, :])
    ax.plot(x, measured, lw=0.9, color=C_E3, label="E3 measured")
    ax.plot(x, predicted, lw=0.9, color=C_E2, ls="--", label="E2 model, predicted")
    ax.set(xlabel="position [mm]", ylabel=r"$\bar{y}$  [µm]")
    ax.set_title("(a) Prediction against measurement", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    tidy(ax)

    ax = fig.add_subplot(gs[1, 0])
    ax.scatter(predicted, measured, s=6, color=C_E3, alpha=0.65, linewidths=0)
    lim = 1.08 * max(np.abs(predicted).max(), np.abs(measured).max())
    ax.plot([-lim, lim], [-lim, lim], lw=0.8, color="#4b5257", ls=(0, (3, 2)))
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="predicted [µm]",
           ylabel="measured [µm]")
    ax.set_title("(b) One point per new target", loc="left", pad=6)
    ax.text(0.04, 0.955, f"r = {frozen['correlation']:.3f}\n"
                         f"$R^2$ = {frozen['variance_explained']:.3f}",
            transform=ax.transAxes, va="top", fontsize=7.4, linespacing=1.4)
    tidy(ax)

    ax = fig.add_subplot(gs[1, 1])
    names = [row["name"] for row in transfer]
    values = [row["e3_residual_um"] for row in transfer]
    colours = [C_E2 if "selected" in n else "#8d979e" for n in names]
    bars = ax.barh(range(len(names)), values, color=colours, height=0.62)
    floor = summary["fits"]["E3_M_C"]["rms_um"]
    ax.axvline(floor, color=C_E3, lw=1.0, ls=(0, (3, 2)))
    ax.text(floor + 0.014, -0.62, f"E3's own fit  {floor:.3f} µm",
            fontsize=6.6, color=C_E3, va="center")
    for i, value in enumerate(values):
        ax.text(value - 0.012, i, f"{value:.3f}", ha="right", va="center",
                fontsize=6.4, color="white", weight="bold")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=6.6)
    ax.set(xlabel="E3 residual rms [µm]", xlim=(0, max(values) * 1.18))
    ax.set_title("(c) What transfers, and what the extra terms buy",
                 loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=C_GRID, lw=0.5)

    caption(fig, 0.095,
            f"E2's model, frozen before this capture existed, takes E3's error from "
            f"{frozen['measured_rms_um']:.3f} µm rms to {frozen['residual_rms_um']:.3f} µm at 400 "
            f"positions drawn from a different seed — a factor of {frozen['reduction_factor']:.2f}, "
            f"explaining {frozen['variance_explained'] * 100:.1f}% of the variance. The stationary "
            f"model M_B leaves {transfer[0]['e3_residual_um']:.3f} µm on the same targets, so the "
            f"varying coefficients are carrying real structure rather than fitting E2's noise. What "
            f"is left, {frozen['residual_rms_um']:.3f} µm, sits above E3's own-fit floor of "
            f"{floor:.3f} µm: in quadrature the day-to-day term is "
            f"{np.sqrt(frozen['residual_rms_um']**2 - floor**2):.3f} µm.")
    e2f.save(fig, "e3_f2_out_of_sample", OUT)


# ---------------------------------------------------------------- figure 3

def figure_f3(model2: dict, model3: dict, grid: np.ndarray, rows: list[dict],
              floor: float) -> None:
    fig = plt.figure(figsize=(7.45, 6.2))
    gs = fig.add_gridspec(2, 2, left=0.095, right=0.975, top=0.850,
                          bottom=0.168, hspace=0.48, wspace=0.26)
    header(fig, 0.095, "E3 — which harmonics reproduce",
           "Amplitude and phase of each screw harmonic along the travel, both sessions")

    ax = fig.add_subplot(gs[0, 0])
    for k, style in zip((1, 2), ("-", "--")):
        a = e2f.harmonic_profile(model2, grid, k)
        b = e2f.harmonic_profile(model3, grid, k)
        ax.plot(grid, a["amplitude"], ls=style, lw=1.1, color=C_E2)
        ax.plot(grid, b["amplitude"], ls=style, lw=1.1, color=C_E3)
        ax.text(grid[-1] + 1.5, b["amplitude"][-1], f"k={k}", fontsize=7,
                color="#4b5257", va="center")
    ax.set(xlabel="position [mm]", ylabel="amplitude [µm]")
    ax.set_title("(a) k = 1 and 2: reproduced", loc="left", pad=6)
    ax.plot([], [], color=C_E2, lw=1.1, label="E2")
    ax.plot([], [], color=C_E3, lw=1.1, label="E3")
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower left")
    tidy(ax)

    ax = fig.add_subplot(gs[0, 1])
    for k, style in zip((3, 4), ("-", "--")):
        a = e2f.harmonic_profile(model2, grid, k)
        b = e2f.harmonic_profile(model3, grid, k)
        ax.plot(grid, a["amplitude"], ls=style, lw=1.1, color=C_E2)
        ax.plot(grid, b["amplitude"], ls=style, lw=1.1, color=C_E3)
        ax.text(grid[-1] + 1.5, b["amplitude"][-1], f"k={k}", fontsize=7,
                color="#4b5257", va="center")
    ax.axhline(floor, color="#4b5257", lw=0.9, ls=(0, (3, 2)))
    ax.text(grid[0] + 1, floor + 0.012, f"own-fit residual {floor:.2f} µm",
            fontsize=6.5, color="#4b5257")
    ax.set(xlabel="position [mm]", ylabel="amplitude [µm]")
    ax.set_title("(b) k = 3 and 4: not reproduced", loc="left", pad=6)
    ax.plot([], [], color=C_E2, lw=1.1, label="E2")
    ax.plot([], [], color=C_E3, lw=1.1, label="E3")
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="center left")
    tidy(ax)

    ax = fig.add_subplot(gs[1, 0])
    k_values = [row["harmonic"] for row in rows]
    width = 0.38
    ax.bar([k - width / 2 for k in k_values],
           [row["E2_mean_amplitude_um"] for row in rows], width, color=C_E2,
           label="E2")
    ax.bar([k + width / 2 for k in k_values],
           [row["E3_mean_amplitude_um"] for row in rows], width, color=C_E3,
           label="E3")
    ax.axhline(floor, color="#4b5257", lw=0.9, ls=(0, (3, 2)))
    ax.set(xticks=k_values, xlabel="harmonic k", ylabel="mean amplitude [µm]")
    ax.set_title("(c) Mean amplitude", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=C_GRID, lw=0.5)

    ax = fig.add_subplot(gs[1, 1])
    shift = [abs(row["mean_phase_shift_deg"]) for row in rows]
    colours = [C_E2 if s < 10 else C_E3 for s in shift]
    ax.bar(k_values, shift, 0.6, color=colours)
    for k, s in zip(k_values, shift):
        ax.text(k, s + 1.6, f"{s:.0f}°", ha="center", fontsize=6.6,
                color="#4b5257")
    ax.set(xticks=k_values, xlabel="harmonic k",
           ylabel="mean phase shift E3 − E2 [deg]", ylim=(0, max(shift) * 1.25))
    ax.set_title("(d) Phase agreement", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=C_GRID, lw=0.5)

    caption(fig, 0.095,
            f"The first harmonic repeats to {abs(rows[0]['amplitude_ratio'] - 1) * 100:.1f}% in "
            f"amplitude and {abs(rows[0]['mean_phase_shift_deg']):.1f}° in phase, and the second to "
            f"{abs(rows[1]['amplitude_ratio'] - 1) * 100:.1f}% and "
            f"{abs(rows[1]['mean_phase_shift_deg']):.0f}°. Harmonics 3 and 4 do not repeat: their "
            f"amplitudes fall by roughly {(1 - rows[2]['amplitude_ratio']) * 100:.0f}% and "
            f"{(1 - rows[3]['amplitude_ratio']) * 100:.0f}% and their phases move "
            f"{abs(rows[2]['mean_phase_shift_deg']):.0f}° and "
            f"{abs(rows[3]['mean_phase_shift_deg']):.0f}°. Both sit at or under the "
            f"{floor:.2f} µm own-fit residual, so they are session-specific structure that the "
            f"within-session cross-validation could not distinguish from physics. E3's own model "
            f"selection, run independently on this capture alone, reaches the same verdict: it "
            f"chooses H = 2 where E2 chose H = 4.")
    e2f.save(fig, "e3_f3_harmonics", OUT)


# ---------------------------------------------------------------- figure 4

def figure_f4(d2: dict, d3: dict, summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 4.4))
    gs = fig.add_gridspec(1, 2, left=0.095, right=0.975, top=0.800,
                          bottom=0.255, wspace=0.26)
    header(fig, 0.095, "E3 — reversal",
           "The up-minus-down difference B, and how much of it is a constant")

    ax = fig.add_subplot(gs[0, 0])
    for data, colour, label in ((d2, C_E2, "E2"), (d3, C_E3, "E3")):
        ax.plot(data["x"], data["D2"]["B"], lw=0.8, color=colour, alpha=0.9,
                label=label)
    for key, colour in (("E2", C_E2), ("E3", C_E3)):
        ax.axhline(summary["reversal"][key]["beta0_um"], color=colour, lw=0.9,
                   ls=(0, (3, 2)))
    ax.set(xlabel="position [mm]", ylabel="B = e↑ − e↓  [µm]")
    ax.set_title("(a) B along the travel", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    tidy(ax)

    ax = fig.add_subplot(gs[0, 1])
    labels = [r"$\beta_0$ [µm]", "slope [µm/100 mm]", "1st harmonic [µm]"]
    e2v = [summary["reversal"]["E2"]["beta0_um"],
           summary["reversal"]["E2"]["slope_um_per_100mm"],
           summary["reversal"]["E2"]["harmonic_amplitude_um"]]
    e3v = [summary["reversal"]["E3"]["beta0_um"],
           summary["reversal"]["E3"]["slope_um_per_100mm"],
           summary["reversal"]["E3"]["harmonic_amplitude_um"]]
    err2 = [summary["reversal"]["E2"]["beta0_sigma_um"],
            summary["reversal"]["E2"]["slope_sigma_um_per_100mm"], 0.0]
    err3 = [summary["reversal"]["E3"]["beta0_sigma_um"],
            summary["reversal"]["E3"]["slope_sigma_um_per_100mm"], 0.0]
    y = np.arange(len(labels))
    ax.barh(y - 0.19, e2v, 0.34, xerr=err2, color=C_E2, label="E2",
            error_kw=dict(lw=0.8, ecolor="#30363b"))
    ax.barh(y + 0.19, e3v, 0.34, xerr=err3, color=C_E3, label="E3",
            error_kw=dict(lw=0.8, ecolor="#30363b"))
    ax.axvline(0, color="#8d979e", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("value (bars show ±1σ where fitted)")
    ax.set_title("(b) Fitted structure of B", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=C_GRID, lw=0.5)

    r2, r3 = summary["reversal"]["E2"], summary["reversal"]["E3"]
    delta = r3["beta0_um"] - r2["beta0_um"]
    sigma = np.hypot(r2["beta0_sigma_um"], r3["beta0_sigma_um"])
    caption(fig, 0.095,
            f"The constant part of reversal moved from {r2['beta0_um']:.3f} µm to "
            f"{r3['beta0_um']:.3f} µm, a change of {abs(delta):.3f} µm or {abs(delta) / sigma:.1f}σ: "
            f"small in size, far too large to be noise. The slope changed sign, from "
            f"{r2['slope_um_per_100mm']:+.3f} to {r3['slope_um_per_100mm']:+.3f} µm per 100 mm, "
            f"which is a real change in how reversal varies along the screw rather than a shift in "
            f"its level. A compensator that applies E2's β₀ unchanged would leave about "
            f"{abs(delta):.2f} µm of uncorrected reversal.")
    e2f.save(fig, "e3_f4_reversal", OUT)


# ------------------------------------------------------------- appendix

def appendix_validity(d3: dict, summary: dict) -> None:
    fig = plt.figure(figsize=(7.45, 5.8))
    gs = fig.add_gridspec(2, 2, left=0.095, right=0.975, top=0.850,
                          bottom=0.095, hspace=0.50, wspace=0.26)
    header(fig, 0.095, "E3 — validity evidence",
           "The checks that decide whether the session may be interpreted")

    ax = fig.add_subplot(gs[0, 0])
    t = d3["anchor_time_s"] / 60.0
    ax.plot(t, d3["anchor_error_um"], marker="o", ms=4.2, lw=1.0, color=C_E3)
    ax.set(xlabel="session time [min]", ylabel="anchor error [µm]")
    ax.set_title("(a) The three drift anchors", loc="left", pad=6)
    for name, time, value in zip(("DATUM", "TURN", "END"), t, d3["anchor_error_um"]):
        ax.annotate(name, (time, value), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=6.6, color=C_MUTED)
    tidy(ax)

    ax = fig.add_subplot(gs[0, 1])
    ax.hist(d3["quality"]["std_um"] * 1000.0, bins=40, color=C_E3)
    ax.set(xlabel="record-window std [nm]", ylabel="stops")
    ax.set_title("(b) Record-window quiet", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=C_GRID, lw=0.5)

    ax = fig.add_subplot(gs[1, 0])
    for name, variant, colour in (("D0 none", "D0", "#8d979e"),
                                  ("D1 linear", "D1", "#a16a2d"),
                                  ("D2 piecewise", "D2", C_E3)):
        ax.plot(d3["x"], d3[variant]["ybar"] - d3["D2"]["ybar"], lw=0.9,
                color=colour, label=name)
    ax.set(xlabel="position [mm]", ylabel=r"$\bar{y}$ − $\bar{y}_{D2}$  [µm]")
    ax.set_title("(c) What the drift variant changes", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)

    ax = fig.add_subplot(gs[1, 1])
    rows = summary["validity"]
    ax.axis("off")
    lines = []
    for row in rows:
        mark = "pass" if row["pass"] else "FAIL"
        value = str(row["value"])
        if len(value) > 26:
            value = value[:23] + "..."
        lines.append(f"{mark:4s}  {row['check'][:38]:38s}  {value}")
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=5.4, family="monospace", linespacing=1.55, color="#2d3338")
    ax.set_title("(d) All checks", loc="left", pad=6)

    e2f.save(fig, "e3_appendix_validity", OUT)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    e2f.style()
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))

    d2 = reduce_session(e2_run, "E2")
    d3 = reduce_session(e3_run, "E3")

    model = FrozenModel()
    x3, th3 = d3["x"], d3["theta"]
    inside = model.in_domain(x3)
    predicted = model.G_rel(x3[inside], th3[inside])
    measured = d3["D2"]["ybar"][inside]
    frozen = dict(summary["frozen_E2_prediction_on_E3"])
    frozen["_x"] = x3[inside]
    frozen["_measured"] = measured - measured.mean()
    frozen["_predicted"] = predicted - predicted.mean()
    frozen["_residual"] = frozen["_measured"] - frozen["_predicted"]

    dom2 = (float(d2["x"].min()), float(d2["x"].max()))
    xi, ti = x3[inside], th3[inside]
    yi = frozen["_measured"]
    transfer = []
    for name, H, Ks, Kc, lam in [("M_B  H2, stationary", 2, 2, None, None),
                                 ("H2, varying", 2, 2, 8, 1.0),
                                 ("H3, varying", 3, 2, 8, 1.0),
                                 ("M_C  H4, varying (selected)", 4, 2, 8, 1.0)]:
        m = e2f.fit_linear(d2["x"], d2["theta"], d2["D2"]["ybar"], H, Ks, Kc,
                           dom2, robust=True, penalty=lam)
        pred = e2f.predict(m, xi, ti)
        pred = pred - pred.mean()
        transfer.append({"name": name,
                         "e3_residual_um": float(np.sqrt(np.mean((yi - pred)**2))),
                         "edf": float(m["edf"])})

    model2_c = own_fit(d2, M_C)
    model3_c = own_fit(d3, M_C)
    grid = np.linspace(max(d2["x"].min(), d3["x"].min()),
                       min(d2["x"].max(), d3["x"].max()), 400)
    floor = summary["fits"]["E3_M_C"]["rms_um"]

    figure_f1(d2, d3, summary)
    figure_f2(frozen, transfer, summary)
    figure_f3(model2_c, model3_c, grid, summary["harmonic_comparison"], floor)
    figure_f4(d2, d3, summary)
    appendix_validity(d3, summary)

    (HERE / "analysis" / "e3_transfer_table.csv").write_text(
        "name,edf,e3_residual_um\n" + "".join(
            f"{row['name']},{row['edf']:.2f},{row['e3_residual_um']:.4f}\n"
            for row in transfer), encoding="utf-8")
    print(f"wrote 5 figures to {OUT.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
