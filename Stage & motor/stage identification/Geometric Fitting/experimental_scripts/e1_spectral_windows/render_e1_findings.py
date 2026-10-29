#!/usr/bin/env python3
"""Analyse the completed E1 run and render the publication findings set.

Outputs three main figures, one harmonics table, one appendix validity figure,
and machine-readable harmonic results. Observations are always markers; lines
are reserved for fitted models and acceptance/reference levels.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import f as f_distribution, norm as normal_distribution

import e1_run
from e1_run import WINDOWS, X_REF_MM


HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "experimental_findings"
ANALYSIS = HERE / "analysis"
STEP_A = HERE.parents[1] / "results" / "step_a_results.json"
CONFIG = HERE / "e1_config.json"

COLOURS = {"W1": "#4c78a8", "W2": "#598a76", "W3": "#9d6770"}
MARKERS = {"W1": "o", "W2": "s", "W3": "^"}
INK = "#252b30"
INK2 = "#59636b"
GRID = "#dfe3e6"
LIMIT = "#a44f45"
COMBINED = "#202428"
HMAX = 8
N_NULL = 10_000


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.edgecolor": "#7d858c",
        "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c",
        "ytick.color": "#41474c",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "legend.fontsize": 7.5,
    })


def harmonic_columns(theta: np.ndarray, harmonics: int) -> np.ndarray:
    cols = [np.ones(theta.size)]
    for k in range(1, harmonics + 1):
        cols.extend((np.cos(k * theta), np.sin(k * theta)))
    return np.column_stack(cols)


def least_squares(x: np.ndarray, y: np.ndarray) -> dict:
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    dof = y.size - x.shape[1]
    variance = float(residual @ residual / dof)
    covariance = variance * np.linalg.pinv(x.T @ x)
    return {"beta": beta, "cov": covariance, "residual": residual,
            "sigma": math.sqrt(variance), "dof": dof}


def spectrum_from_fit(fit: dict, coefficient_start: int = 1) -> dict:
    amplitudes, amplitude_sigma, phases, phase_sigma = [], [], [], []
    beta, covariance = fit["beta"], fit["cov"]
    for k in range(HMAX):
        i = coefficient_start + 2 * k
        c, s = beta[i:i + 2]
        block = covariance[i:i + 2, i:i + 2]
        amplitude = float(np.hypot(c, s))
        if amplitude > 1e-15:
            g_amp = np.array([c, s]) / amplitude
            sigma_amp = float(np.sqrt(max(g_amp @ block @ g_amp, 0.0)))
            # e = A cos(k theta + phi), hence phi = atan2(-s, c).
            g_phi = np.array([s, -c]) / (amplitude * amplitude)
            sigma_phi = math.degrees(float(np.sqrt(max(g_phi @ block @ g_phi, 0.0))))
        else:
            sigma_amp, sigma_phi = float("nan"), float("nan")
        amplitudes.append(amplitude)
        amplitude_sigma.append(sigma_amp)
        phases.append(float(np.degrees(np.arctan2(-s, c)) % 360.0))
        phase_sigma.append(sigma_phi)
    return {"amplitude": np.array(amplitudes),
            "amplitude_sigma": np.array(amplitude_sigma),
            "phase_deg": np.array(phases),
            "phase_sigma_deg": np.array(phase_sigma)}


def stop_quality(run: e1_run.Run) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means, stds, slopes = [], [], []
    for start, end in zip(run.stops.record_start_s, run.stops.record_end_s):
        i, j = run.index_at([start, end])
        values = run.pos[i:j] * 1000.0
        seconds = run.t[i:j]
        if values.size < 3:
            raise RuntimeError("record window contains fewer than three samples")
        means.append(float(values.mean() / 1000.0))
        stds.append(float(values.std(ddof=1)))
        centred = seconds - seconds.mean()
        slopes.append(float(np.dot(centred, values - values.mean()) /
                            np.dot(centred, centred)))
    return np.array(means), np.array(stds), np.array(slopes)


def prepare_windows(run: e1_run.Run, means_mm: np.ndarray,
                    gain_slope: float, gain_centre: float) -> dict:
    st = run.stops
    error_um = (means_mm - st.target_mm) * 1000.0
    windows = {}
    for name, low, high in WINDOWS:
        mask = st.window == name
        up = np.flatnonzero(mask & (st.direction == 1))
        down = np.flatnonzero(mask & (st.direction == -1))
        up = up[np.argsort(st.target_mm[up])]
        down = down[np.argsort(st.target_mm[down])]
        x = st.target_mm[up]
        if not np.allclose(x, st.target_mm[down], atol=1e-12):
            raise RuntimeError(f"{name}: upward/downward target mismatch")
        gain = gain_slope * (x - gain_centre)
        e_up = error_um[up] - gain
        e_down = error_um[down] - gain
        ybar = 0.5 * (e_up + e_down)
        offset = float(ybar.mean())
        e_up -= offset
        e_down -= offset
        ybar -= offset
        reversal = e_up - e_down
        target_usteps = np.rint(x * 1600.0).astype(int)
        phase_bin = np.mod(target_usteps, 3200)
        theta = 2.0 * np.pi * phase_bin / 3200.0
        windows[name] = {
            "x": x, "low": low, "high": high, "theta": theta,
            "phase_deg": np.degrees(theta), "e_up": e_up,
            "e_down": e_down, "ybar": ybar, "B": reversal,
            "offset_um": offset,
        }
    return windows


def combined_design(windows: dict, key: str, harmonics: int,
                    separate_harmonics: bool = False) -> tuple[np.ndarray, np.ndarray]:
    names = list(COLOURS)
    n = sum(windows[name][key].size for name in names)
    if separate_harmonics:
        design = np.zeros((n, len(names) * (1 + 2 * harmonics)))
    else:
        design = np.zeros((n, len(names) + 2 * harmonics))
    values = np.empty(n)
    row = 0
    for w, name in enumerate(names):
        y = windows[name][key]
        theta = windows[name]["theta"]
        rows = slice(row, row + y.size)
        values[rows] = y
        if separate_harmonics:
            base = w * (1 + 2 * harmonics)
            design[rows, base] = 1.0
            for k in range(1, harmonics + 1):
                design[rows, base + 2 * k - 1] = np.cos(k * theta)
                design[rows, base + 2 * k] = np.sin(k * theta)
        else:
            design[rows, w] = 1.0
            for k in range(1, harmonics + 1):
                design[rows, len(names) + 2 * k - 2] = np.cos(k * theta)
                design[rows, len(names) + 2 * k - 1] = np.sin(k * theta)
        row += y.size
    return design, values


def phase_repeat_sigma(windows: dict) -> float:
    sum_sq, degrees = 0.0, 0
    for data in windows.values():
        bins = np.rint(data["phase_deg"] / 360.0 * 3200.0).astype(int)
        for phase in np.unique(bins):
            values = data["ybar"][bins == phase]
            sum_sq += float(np.sum((values - values.mean()) ** 2))
            degrees += values.size - 1
    return math.sqrt(sum_sq / degrees)


def false_alarm_amplitude(windows: dict, sigma: float) -> float:
    """1% family-wise amplitude threshold over W1-W3, combined, and k=1..8."""
    rng = np.random.default_rng(20261027)
    maximum = np.zeros(N_NULL)
    simulated = []
    for name in COLOURS:
        theta = windows[name]["theta"]
        design = harmonic_columns(theta, HMAX)
        noise = rng.normal(0.0, sigma, (theta.size, N_NULL))
        beta = np.linalg.pinv(design) @ noise
        amp = np.hypot(beta[1::2], beta[2::2])
        maximum = np.maximum(maximum, amp.max(axis=0))
        simulated.append(noise)
    combined_x, _ = combined_design(windows, "ybar", HMAX)
    combined_noise = np.vstack(simulated)
    beta = np.linalg.pinv(combined_x) @ combined_noise
    start = len(COLOURS)
    amp = np.hypot(beta[start::2], beta[start + 1::2])
    maximum = np.maximum(maximum, amp.max(axis=0))
    return float(np.quantile(maximum, 0.99))


def circular_difference(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def stationarity_verdict(spectra: dict, significant: np.ndarray) -> tuple[bool, list]:
    checks = []
    names = list(COLOURS)
    for k in np.flatnonzero(significant):
        harmonic_ok = True
        comparisons = []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                left, right = spectra[names[i]], spectra[names[j]]
                da = abs(left["amplitude"][k] - right["amplitude"][k])
                la = 2.0 * math.hypot(left["amplitude_sigma"][k],
                                      right["amplitude_sigma"][k])
                dp = circular_difference(left["phase_deg"][k], right["phase_deg"][k])
                lp = 2.0 * math.hypot(left["phase_sigma_deg"][k],
                                      right["phase_sigma_deg"][k])
                passed = da <= la and dp <= lp
                harmonic_ok &= passed
                comparisons.append({"windows": [names[i], names[j]],
                                    "amplitude_difference_um": da,
                                    "amplitude_2sigma_um": la,
                                    "phase_difference_deg": dp,
                                    "phase_2sigma_deg": lp,
                                    "pass": bool(passed)})
        checks.append({"harmonic": int(k + 1), "pass": bool(harmonic_ok),
                       "comparisons": comparisons})
    return all(item["pass"] for item in checks), checks


def reversal_test(windows: dict, harmonics: int) -> dict:
    null_x, values = combined_design(windows, "B", 0)
    full_x, _ = combined_design(windows, "B", harmonics,
                                separate_harmonics=True)
    null = least_squares(null_x, values)
    full = least_squares(full_x, values)
    df_num = full_x.shape[1] - null_x.shape[1]
    df_den = full["dof"]
    rss_null = float(null["residual"] @ null["residual"])
    rss_full = float(full["residual"] @ full["residual"])
    statistic = ((rss_null - rss_full) / df_num) / (rss_full / df_den)
    p_value = float(f_distribution.sf(statistic, df_num, df_den))
    return {"f": float(statistic), "df_num": int(df_num),
            "df_den": int(df_den), "p": p_value,
            "varies": bool(p_value < 0.01)}


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=320)
    plt.close(fig)


def figure_raw(windows: dict, gain_ppm: float) -> None:
    maximum = max(np.max(np.abs(data[key])) for data in windows.values()
                  for key in ("e_up", "e_down", "ybar"))
    limit = math.ceil(maximum * 1.12 * 2.0) / 2.0
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.2), sharey=True)
    fig.subplots_adjust(left=0.12, right=0.975, top=0.965, bottom=0.155,
                        hspace=0.30)
    for panel, (ax, (name, data)) in enumerate(zip(axes, windows.items())):
        colour = COLOURS[name]
        ax.scatter(data["x"], data["e_up"], s=8, marker="^", color=colour,
                   alpha=0.22, linewidths=0, label=r"$e_{\uparrow}$" if panel == 0 else None)
        ax.scatter(data["x"], data["e_down"], s=8, marker="v", color=colour,
                   alpha=0.22, linewidths=0, label=r"$e_{\downarrow}$" if panel == 0 else None)
        ax.scatter(data["x"], data["ybar"], s=15, marker="o", facecolor="white",
                   edgecolor=colour, linewidth=0.75,
                   label=r"$\bar{y}$" if panel == 0 else None, zorder=3)
        ax.axhline(0, color="#8c949a", lw=0.6)
        ax.set_xlim(data["low"] - 0.25, data["high"] + 0.25)
        ax.set_ylim(-limit, limit)
        ax.set_ylabel("corrected error [µm]")
        ax.set_title(f"({chr(97 + panel)}) {name}  ·  {data['low']:.0f}–{data['high']:.0f} mm",
                     loc="left", pad=5)
        ax.grid(color=GRID, lw=0.48)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(loc="upper right", frameon=False, ncol=3)
    axes[-1].set_xlabel("commanded position [mm]")
    fig.text(0.12, 0.030,
             f"Markers are measured stops; no fitted curve is shown. Step A gain ({gain_ppm:+.2f} ppm) and the "
             "direction-averaged offset of each window were subtracted.\nEach panel spans six screw revolutions; "
             "the common y-scale permits direct amplitude comparison.",
             fontsize=7.2, color=INK2)
    save(fig, "e1_findings_figure1_raw_windows")


def fit_curve(data: dict, key: str, harmonics: int,
              degrees: np.ndarray) -> np.ndarray:
    fit = least_squares(harmonic_columns(data["theta"], harmonics), data[key])
    return harmonic_columns(np.radians(degrees), harmonics) @ fit["beta"]


def figure_folded(windows: dict, harmonics: int) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.975, top=0.895, bottom=0.17,
                        hspace=0.24)
    degrees = np.linspace(0, 360, 1000)
    for name, data in windows.items():
        colour, marker = COLOURS[name], MARKERS[name]
        order = np.argsort(data["phase_deg"])
        axes[0].scatter(data["phase_deg"][order], data["ybar"][order],
                        s=12, marker=marker, facecolor="white", edgecolor=colour,
                        linewidth=0.65, alpha=0.72, label=name)
        axes[0].plot(degrees, fit_curve(data, "ybar", harmonics, degrees),
                     color=colour, lw=1.45)
        axes[1].scatter(data["phase_deg"][order], data["B"][order],
                        s=11, marker=marker, facecolor="white", edgecolor=colour,
                        linewidth=0.62, alpha=0.68, label=name)
        axes[1].plot(degrees, fit_curve(data, "B", harmonics, degrees),
                     color=colour, lw=1.25)
    mean_b = float(np.mean(np.concatenate([data["B"] for data in windows.values()])))
    for ax in axes:
        ax.axhline(0, color="#8c949a", lw=0.6)
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 60))
        ax.grid(color=GRID, lw=0.48)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"direction average $\bar{y}$ [µm]")
    axes[0].set_title(f"(a) Direction-averaged error; order-{harmonics} harmonic fits", loc="left", pad=5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.54, 0.982), columnspacing=1.8, handletextpad=0.5)
    axes[1].axhline(mean_b, color=INK2, lw=1.05, ls=(0, (5, 3)), zorder=2)
    axes[1].set_ylabel(r"reversal $B=e_{\uparrow}-e_{\downarrow}$ [µm]")
    axes[1].set_xlabel("absolute screw angle θ [deg]")
    axes[1].set_title("(b) Reversal on the same absolute screw reference", loc="left", pad=5)
    fig.text(0.115, 0.030,
             "Absolute angle is computed from commanded microsteps modulo 3200. Markers are measurements; solid lines "
             "are least-squares harmonic fits.\nStep A gain and per-window offsets are removed only from ȳ; reversal B "
             "is unchanged by those subtractions.\n"
             f"The dashed line in (b) is the across-window mean, B̄ = {mean_b:.3f} µm.",
             fontsize=7.2, color=INK2)
    save(fig, "e1_findings_figure2_folded_revolution")


def figure_spectrum(spectra: dict, threshold: float, significant: np.ndarray,
                    repeat_sigma: float) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.35))
    fig.subplots_adjust(left=0.115, right=0.975, top=0.95, bottom=0.235)
    harmonics = np.arange(1, HMAX + 1)
    for name in COLOURS:
        ax.scatter(harmonics, spectra[name]["amplitude"], marker=MARKERS[name],
                   s=27, facecolor="white", edgecolor=COLOURS[name], linewidth=0.95,
                   label=name, zorder=3)
    ax.scatter(harmonics, spectra["combined"]["amplitude"], marker="D", s=24,
               color=COMBINED, label="combined", zorder=4)
    ax.axhline(threshold, color=LIMIT, lw=1.0, ls=(0, (4, 2)),
               label="1% family-wise false-alarm level")
    for k in harmonics[significant]:
        ax.axvspan(k - 0.16, k + 0.16, color="#dce8df", alpha=0.65, lw=0)
    ax.set(xlim=(0.7, 8.3), xticks=harmonics, xlabel="harmonic number k [cycles per revolution]",
           ylabel="least-squares amplitude [µm]")
    ax.set_title("Amplitude spectrum on the absolute screw-angle basis", loc="left", pad=6)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.text(0.115, 0.035,
             "The 1% threshold controls the maximum across W1–W3, the combined spectrum, and k=1–8.\n"
             f"Its null scale (σ={repeat_sigma:.3f} µm) comes from repeated absolute phases and includes unmodelled "
             "slow structure\n(Appendix validity panel d), making the threshold conservative. Green bands mark "
             "harmonics above threshold in at least two windows.",
             fontsize=7.15, color=INK2)
    save(fig, "e1_findings_figure3_harmonic_spectrum")


def render_table(spectra: dict, threshold: float, significant: np.ndarray,
                 chosen_h: int, stationary: bool, reversal: dict) -> tuple[list[dict], dict]:
    rows, csv_rows = [], []
    for k in range(HMAX):
        amplitudes = [spectra[name]["amplitude"][k] for name in COLOURS]
        sigmas = [spectra[name]["amplitude_sigma"][k] for name in COLOURS]
        phases = [spectra[name]["phase_deg"][k] for name in COLOURS]
        count = sum(value > threshold for value in amplitudes)
        label = f"Yes ({count}/3)" if significant[k] else f"No ({count}/3)"
        rows.append([str(k + 1)] +
                    [f"{a:.3f} ± {s:.3f}" for a, s in zip(amplitudes, sigmas)] +
                    [" / ".join(f"{p:.0f}°" for p in phases), label])
        csv_rows.append({
            "harmonic": k + 1,
            **{f"{name}_amplitude_um": spectra[name]["amplitude"][k] for name in COLOURS},
            **{f"{name}_amplitude_sigma_um": spectra[name]["amplitude_sigma"][k] for name in COLOURS},
            **{f"{name}_phase_deg": spectra[name]["phase_deg"][k] for name in COLOURS},
            "significant_windows": count,
            "significant_at_1pct_fwer": bool(significant[k]),
        })
    a1 = float(spectra["W1"]["amplitude"][0])
    a2 = float(spectra["W2"]["amplitude"][0])
    s1 = float(spectra["W1"]["amplitude_sigma"][0])
    s2 = float(spectra["W2"]["amplitude_sigma"][0])
    ratio = a1 / a2
    ratio_sigma = ratio * math.hypot(s1 / a1, s2 / a2)
    z_score = abs(a1 - a2) / math.hypot(s1, s2)
    p_value = float(2.0 * normal_distribution.sf(z_score))
    contrast = {
        "harmonic": 1,
        "contrast": "W1/W2 amplitude ratio",
        "ratio": ratio,
        "ratio_sigma": ratio_sigma,
        "difference_z": z_score,
        "two_sided_p": p_value,
        "significant_at_0p05": bool(p_value < 0.05),
    }

    columns = ["Harmonic k", "W1 amplitude ± σ", "W2 amplitude ± σ",
               "W3 amplitude ± σ", "Phase W1 / W2 / W3", "Significant?"]

    fig, ax = plt.subplots(figsize=(8.25, 4.85))
    fig.subplots_adjust(left=0.035, right=0.965, top=0.90, bottom=0.20)
    ax.axis("off")
    ax.text(0.0, 1.055, "Table 1. Harmonic coefficients by E1 window (k = 1–4)",
            transform=ax.transAxes, fontsize=10.2, weight="bold", color=INK, va="bottom")
    table = ax.table(cellText=rows[:4], colLabels=columns, cellLoc="center", colLoc="center",
                     colWidths=[0.09, 0.16, 0.16, 0.16, 0.25, 0.14],
                     bbox=[0.0, 0.40, 1.0, 0.59])
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#c5c9cc")
        cell.set_linewidth(0.55)
        if row == 0:
            cell.set_facecolor("#e9ecee")
            cell.set_text_props(weight="bold", color=INK)
        elif significant[row - 1]:
            cell.set_facecolor("#edf4ee")
        else:
            cell.set_facecolor("white")

    contrast_table = ax.table(
        cellText=[["k = 1 amplitude: W1/W2", f"{ratio:.3f} ± {ratio_sigma:.3f}",
                   f"z = {z_score:.2f}; p = {p_value:.2e}",
                   "Yes" if p_value < 0.05 else "No"]],
        colLabels=["Non-stationarity contrast", "Ratio ± σ", "Difference test", "Significant?"],
        cellLoc="center", colLoc="center", colWidths=[0.29, 0.21, 0.30, 0.20],
        bbox=[0.0, 0.10, 1.0, 0.21])
    contrast_table.auto_set_font_size(False)
    contrast_table.set_fontsize(8.0)
    for (row, col), cell in contrast_table.get_celld().items():
        cell.set_edgecolor("#c5c9cc")
        cell.set_linewidth(0.55)
        if row == 0:
            cell.set_facecolor("#e9ecee")
            cell.set_text_props(weight="bold", color=INK)
        elif col == 3 and p_value < 0.05:
            cell.set_facecolor("#f4e8e6")
        else:
            cell.set_facecolor("white")
    verdict = "yes" if stationary else "no"
    reversal_text = "yes" if reversal["varies"] else "no"
    fig.text(0.035, 0.145,
             f"Chosen H: {chosen_h}    ·    Stationary within 2σ: {verdict}    ·    "
             f"Reversal varies per revolution: {reversal_text}",
             fontsize=8.6, weight="bold", color=INK)
    fig.text(0.035, 0.045,
             f"Harmonic significance requires the 1% family-wise threshold ({threshold:.3f} µm) in at least two "
             "windows.\nThe ratio uncertainty uses first-order independent-fit propagation; its p-value is a "
             "two-sided normal difference test.\nPhase uses e = A cos(kθ + φ) on the absolute screw reference.",
             fontsize=7.2, color=INK2)
    save(fig, "e1_findings_table_harmonics")

    fig, ax = plt.subplots(figsize=(8.25, 4.75))
    fig.subplots_adjust(left=0.035, right=0.965, top=0.91, bottom=0.15)
    ax.axis("off")
    ax.text(0.0, 1.065, "Appendix Table A1. Full E1 harmonic coefficients",
            transform=ax.transAxes, fontsize=10.2, weight="bold", color=INK, va="bottom")
    table = ax.table(cellText=rows, colLabels=columns, cellLoc="center", colLoc="center",
                     colWidths=[0.09, 0.16, 0.16, 0.16, 0.25, 0.14],
                     bbox=[0.0, 0.08, 1.0, 0.90])
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#c5c9cc")
        cell.set_linewidth(0.55)
        if row == 0:
            cell.set_facecolor("#e9ecee")
            cell.set_text_props(weight="bold", color=INK)
        elif significant[row - 1]:
            cell.set_facecolor("#edf4ee")
        else:
            cell.set_facecolor("white")
    fig.text(0.035, 0.035,
             f"Amplitude significance requires the 1% family-wise threshold ({threshold:.3f} µm) in at least two "
             "windows. Phase uses e = A cos(kθ + φ) on the absolute screw reference.",
             fontsize=7.2, color=INK2)
    save(fig, "e1_findings_appendix_table_harmonics")
    return csv_rows, contrast


def figure_appendix(run: e1_run.Run, windows: dict, stds: np.ndarray,
                    slopes: np.ndarray, config: dict, chosen_h: int) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.25, 7.1))
    fig.subplots_adjust(left=0.095, right=0.965, top=0.96, bottom=0.13,
                        hspace=0.38, wspace=0.30)

    ax = axes[0, 0]
    ref = run.stops.window == "REF"
    means = run.hold()[ref] * 1000.0
    times = run.stops.record_start_s[ref] / 60.0
    datum = means[run.stops.state[ref] == "DATUM"].mean()
    drift = means - datum
    states = run.stops.state[ref]
    groups = [
        ("datum", states == "DATUM"),
        ("R1", (states == "REF") & (run.stops.index[ref] == 1)),
        ("R2", (states == "REF") & (run.stops.index[ref] == 2)),
        ("end", states == "END_REF"),
    ]
    for label, mask in groups:
        x, y = float(times[mask].mean()), float(drift[mask].mean())
        spread = float(np.ptp(drift[mask]) / 2.0) if mask.sum() > 1 else 0.0
        ax.errorbar(x, y, yerr=spread, fmt="o", ms=5.0, mfc="white",
                    mec="#73578f", mew=0.9, ecolor="#73578f", elinewidth=0.8,
                    capsize=2.0)
        ax.annotate(label, (x, y), xytext=(0, 7), textcoords="offset points",
                    ha="center", fontsize=6.9, color=INK2)
    warning = config["reference_drift_warning_um"]
    ax.axhline(0, color="#8c949a", lw=0.6)
    ax.axhline(warning, color=LIMIT, lw=0.9, ls="--", label="±0.5 µm warning")
    ax.axhline(-warning, color=LIMIT, lw=0.9, ls="--")
    ax.set(xlabel="ESP experiment time [min]", ylabel="reference drift from datum mean [µm]")
    ax.set_title("(a) Fixed-position reference stability", loc="left", pad=5)
    ax.legend(frameon=False, loc="lower left")
    ax.grid(color=GRID, lw=0.48)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[0, 1]
    std_limit = config["expected_noise_um"]
    std_lo = max(float(stds.min()) * 0.8, 1e-5)
    std_bins = np.geomspace(std_lo, std_limit * 1.12, 32)
    ax.hist(stds, bins=std_bins, color="#6f8794", alpha=0.90,
            edgecolor="white", lw=0.35)
    ax.axvline(std_limit, color=LIMIT, lw=1.0, ls="--", label=f"limit {std_limit:.2f} µm")
    ax.set_xscale("log")
    ax.set(xlabel="record-window standard deviation [µm]", ylabel="stop count")
    ax.set_title("(b) In-window scatter", loc="left", pad=5)
    ax.legend(frameon=False)
    ax.grid(axis="y", color=GRID, lw=0.48)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 0]
    slope_effect = np.abs(slopes) * config["record_s"]
    effect_limit = 0.5 * config["expected_noise_um"]
    positive_effect = slope_effect[slope_effect > 0]
    effect_lo = max(float(positive_effect.min()) * 0.8, 1e-6)
    effect_bins = np.geomspace(effect_lo, effect_limit * 1.12, 32)
    ax.hist(positive_effect, bins=effect_bins, color="#9a7b4f", alpha=0.90,
            edgecolor="white", lw=0.35)
    ax.axvline(effect_limit, color=LIMIT, lw=1.0, ls="--",
               label=f"limit {effect_limit:.2f} µm")
    ax.set_xscale("log")
    ax.set(xlabel=r"$|$slope$|\,t_{record}$ [µm]", ylabel="stop count")
    ax.set_title("(c) Drift during each record window", loc="left", pad=5)
    ax.legend(frameon=False)
    ax.grid(axis="y", color=GRID, lw=0.48)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 1]
    maximum = 0.0
    for name, data in windows.items():
        design = harmonic_columns(data["theta"], chosen_h)
        fit = least_squares(design, data["ybar"])
        residual = fit["residual"]
        maximum = max(maximum, float(np.max(np.abs(residual))))
        ax.scatter(data["x"], residual, s=12, marker=MARKERS[name],
                   facecolor="white", edgecolor=COLOURS[name], linewidth=0.65,
                   alpha=0.78, label=name)
    ax.axhline(0, color="#8c949a", lw=0.6)
    ax.set(xlabel="commanded position [mm]", ylabel="harmonic-fit residual [µm]")
    ax.set_title(f"(d) Residual after the order-{chosen_h} fit", loc="left", pad=5)
    ax.legend(frameon=False, ncol=3)
    ax.grid(color=GRID, lw=0.48)
    ax.spines[["top", "right"]].set_visible(False)

    ref_drift = float(run.hold()[ref][-3:].mean() - run.hold()[ref][:3].mean()) * 1000.0
    std_pass = float(np.mean(stds <= std_limit) * 100.0)
    slope_pass = float(np.mean(slope_effect <= effect_limit) * 100.0)
    fig.text(0.095, 0.035,
             f"All 654 logged stops are included in the quality histograms. Reference datum-to-end drift: {ref_drift:+.3f} µm; "
             f"std pass: {std_pass:.1f}%; slope-effect pass: {slope_pass:.1f}%; largest harmonic residual: {maximum:.3f} µm.\n"
             f"Panel (d) retains coherent position-dependent structure after the order-{chosen_h} fit; the residual is not white noise.",
             fontsize=7.2, color=INK2)
    save(fig, "e1_findings_appendix_validity")


def main() -> int:
    style()
    run = e1_run.load()
    means, stds, slopes = stop_quality(run)
    step_a = json.loads(STEP_A.read_text(encoding="utf-8"))["fit_uncorrected"]
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    gain_slope = float(step_a["slope_um_per_mm"])
    gain_centre = float(step_a["x_center_mm"])
    windows = prepare_windows(run, means, gain_slope, gain_centre)

    fits, spectra = {}, {}
    for name, data in windows.items():
        fits[name] = least_squares(harmonic_columns(data["theta"], HMAX), data["ybar"])
        spectra[name] = spectrum_from_fit(fits[name])
    combined_x, combined_y = combined_design(windows, "ybar", HMAX)
    combined_fit = least_squares(combined_x, combined_y)
    spectra["combined"] = spectrum_from_fit(combined_fit, len(COLOURS))

    repeat_sigma = phase_repeat_sigma(windows)
    threshold = false_alarm_amplitude(windows, repeat_sigma)
    exceed = np.vstack([spectra[name]["amplitude"] > threshold for name in COLOURS])
    significant = exceed.sum(axis=0) >= 2
    chosen_h = int(np.flatnonzero(significant)[-1] + 1) if significant.any() else 0
    if chosen_h == 0:
        raise RuntimeError("no screw harmonic meets the predeclared significance rule")
    stationary, stationarity_checks = stationarity_verdict(spectra, significant)
    reversal = reversal_test(windows, chosen_h)

    figure_raw(windows, float(step_a["slope_ppm"]))
    figure_folded(windows, chosen_h)
    figure_spectrum(spectra, threshold, significant, repeat_sigma)
    table_rows, nonstationarity = render_table(spectra, threshold, significant, chosen_h,
                                               stationary, reversal)
    figure_appendix(run, windows, stds, slopes, config, chosen_h)

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    with (ANALYSIS / "e1_harmonics_table.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(table_rows[0]))
        writer.writeheader(); writer.writerows(table_rows)
    serial_spectra = {
        name: {key: value.tolist() for key, value in spectrum.items()}
        for name, spectrum in spectra.items()
    }
    result = {
        "experiment": "E1_spectral_windows",
        "source_stops": int(run.stops.state.size),
        "analysis_stops": int(sum(data["ybar"].size * 2 for data in windows.values())),
        "step_a_gain_removed": {
            "slope_um_per_mm": gain_slope,
            "slope_ppm": float(step_a["slope_ppm"]),
            "x_center_mm": gain_centre,
            "source": str(STEP_A),
        },
        "harmonic_basis": "e = a0 + sum[c_k cos(k theta) + s_k sin(k theta)]",
        "phase_convention": "phi_deg = atan2(-s_k, c_k); e_k = A_k cos(k theta + phi_k)",
        "absolute_phase": "theta = 2 pi (commanded_usteps mod 3200) / 3200",
        "hmax_tested": HMAX,
        "repeat_noise_sigma_um": repeat_sigma,
        "false_alarm": {
            "probability": 0.01,
            "family": "maximum over W1, W2, W3, combined and harmonics 1..8",
            "simulations": N_NULL,
            "amplitude_threshold_um": threshold,
        },
        "significant_harmonics": [int(k + 1) for k in np.flatnonzero(significant)],
        "chosen_H": chosen_h,
        "stationary_within_2sigma": bool(stationary),
        "stationarity_checks": stationarity_checks,
        "nonstationarity_k1_w1_w2": nonstationarity,
        "reversal_periodic_test": reversal,
        "quality": {
            "std_limit_um": float(config["expected_noise_um"]),
            "std_pass_fraction": float(np.mean(stds <= config["expected_noise_um"])),
            "slope_effect_limit_um": float(0.5 * config["expected_noise_um"]),
            "slope_pass_fraction": float(np.mean(np.abs(slopes) * config["record_s"] <=
                                                 0.5 * config["expected_noise_um"])),
        },
        "spectra": serial_spectra,
    }
    (ANALYSIS / "e1_findings_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(f"chosen H={chosen_h}; significant harmonics="
          f"{','.join(map(str, result['significant_harmonics']))}")
    print(f"1% family-wise amplitude threshold={threshold:.4f} µm; "
          f"phase-repeat sigma={repeat_sigma:.4f} µm")
    print(f"stationary within 2 sigma={stationary}; "
          f"reversal varies per revolution={reversal['varies']} (p={reversal['p']:.3g})")
    for path in sorted(OUT.glob("e1_findings_*")):
        print("wrote", path.relative_to(HERE))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
