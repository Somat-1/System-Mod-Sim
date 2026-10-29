#!/usr/bin/env python3
"""Process E4, refit E1, and render the stationarity-profile findings set."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import block_diag
from scipy.optimize import least_squares, minimize_scalar
from scipy.stats import pearsonr

import e4_run

HERE = Path(__file__).resolve().parent
E1_DIR = HERE.parent / "e1_spectral_windows"
if str(E1_DIR) not in sys.path:
    sys.path.insert(0, str(E1_DIR))
import e1_run  # noqa: E402

ANALYSIS = HERE / "analysis"
OUT = HERE / "rendered" / "experimental_findings"
PRESENTATION = HERE / "rendered" / "presentation"
STEP_A = HERE.parents[1] / "results" / "step_a_results.json"
CONFIG = HERE / "e4_config.json"

# Physical endpoint assumptions from the 122.770 mm homed travel. These are
# deliberately centralized because the processing brief requires measurement.
MOTOR_END_MM = 0.0
SUPPORT_BEARING_MM = 122.770
ENDPOINT_STATUS = "inferred from homed travel; mechanically confirm"
STD_LIMIT_UM = 0.5
SLOPE_EFFECT_LIMIT_UM = 0.25
SATURATION_INCREMENT_NM = 8333.0
K_FIT = 4
K_HEATMAP = 8
INK = "#252b30"
INK2 = "#59636b"
GRID = "#dfe3e6"
LIMIT = "#a44f45"
HARMONIC_COLOURS = {1: "#315b7d", 2: "#a3553d"}
REVERSAL_COLOUR = "#2f6b4f"
POSITION_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "e4_position", ("#243b62", "#24677a", "#2f8a7a", "#9a7b3f"))
POSITION_NORM = plt.Normalize(MOTOR_END_MM, SUPPORT_BEARING_MM)


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8,
        "axes.titlesize": 9.5, "axes.labelsize": 9,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "legend.fontsize": 7.2, "axes.unicode_minus": True,
    })


def save(fig: plt.Figure, stem: str, directory: Path = OUT) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{stem}.png", dpi=320)
    plt.close(fig)


def endpoint_lines(ax: plt.Axes) -> None:
    ax.axvline(MOTOR_END_MM, color="#725247", lw=0.8, ls=":", zorder=0)
    ax.axvline(SUPPORT_BEARING_MM, color="#526d72", lw=0.8, ls=":", zorder=0)


def endpoint_bands(ax: plt.Axes) -> None:
    """Show physical ends without implying fitted values there."""
    ax.axvspan(MOTOR_END_MM - 1.4, MOTOR_END_MM + 1.4,
               color="#725247", alpha=.09, lw=0, zorder=0)
    ax.axvspan(SUPPORT_BEARING_MM - 1.4, SUPPORT_BEARING_MM + 1.4,
               color="#526d72", alpha=.09, lw=0, zorder=0)
    transform = ax.get_xaxis_transform()
    ax.text(MOTOR_END_MM + 1.8, 1.025, "motor end", transform=transform,
            ha="left", va="bottom", fontsize=6.8, color=INK2, clip_on=False)
    ax.text(SUPPORT_BEARING_MM - 1.8, 1.025, "support bearing", transform=transform,
            ha="right", va="bottom", fontsize=6.8, color=INK2, clip_on=False)


def panel_title(ax: plt.Axes, label: str, title: str) -> None:
    # Leave a separate line immediately above the axes for endpoint labels.
    ax.set_title(rf"$\bf{{({label})}}$ {title}", loc="left", y=1.065, pad=4)


def quality(run) -> dict:
    means, stds, slopes, saturated = [], [], [], []
    for start, end in zip(run.stops.record_start_s, run.stops.record_end_s):
        i, j = run.index_at([start, end])
        values = run.pos[i:j] * 1000.0
        seconds = run.t[i:j]
        means.append(float(values.mean() / 1000.0))
        stds.append(float(values.std(ddof=1)))
        centred = seconds - seconds.mean()
        slopes.append(float(np.dot(centred, values - values.mean()) /
                            np.dot(centred, centred)))
        saturated.append(int(np.count_nonzero(np.abs(np.diff(run.pos[i:j])) * 1e6
                                               >= SATURATION_INCREMENT_NM)))
    return {"means_mm": np.asarray(means), "std_um": np.asarray(stds),
            "slope_um_s": np.asarray(slopes), "saturated": np.asarray(saturated)}


def drift_corrected_error(run, means_mm: np.ndarray) -> tuple[np.ndarray, dict]:
    raw = (means_mm - run.stops.target_mm) * 1000.0
    ref = run.stops.window == "REF"
    ref_t = 0.5 * (run.stops.record_start_s[ref] + run.stops.record_end_s[ref])
    ref_error = raw[ref]
    all_t = 0.5 * (run.stops.record_start_s + run.stops.record_end_s)
    interpolated = np.interp(all_t, ref_t, ref_error)
    datum = float(ref_error.mean())
    corrected = raw - interpolated
    return corrected, {
        "times_min": ref_t / 60.0, "raw_error_um": ref_error,
        "drift_from_datum_um": ref_error - datum, "datum_um": datum,
        "datum_to_end_um": float(ref_error[-3:].mean() - ref_error[:3].mean()),
        "spread_um": float(np.ptp(ref_error)),
    }


def linear_fit(design: np.ndarray, values: np.ndarray) -> dict:
    beta, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ beta
    dof = values.size - design.shape[1]
    variance = float(residual @ residual / dof)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return {"beta": beta, "cov": covariance, "residual": residual,
            "dof": dof, "sigma": math.sqrt(variance),
            "condition": float(np.linalg.cond(design))}


def geometry_design(x: np.ndarray, theta: np.ndarray, harmonics: int,
                    with_linear: bool = True) -> np.ndarray:
    columns = [np.ones(x.size)]
    if with_linear:
        columns.append(x - x.mean())
    for k in range(1, harmonics + 1):
        columns.extend((np.cos(k * theta), np.sin(k * theta)))
    return np.column_stack(columns)


def spectrum(fit: dict, harmonics: int, with_linear: bool = True) -> list[dict]:
    start = 2 if with_linear else 1
    result = []
    for k in range(1, harmonics + 1):
        index = start + 2 * (k - 1)
        a, b = fit["beta"][index:index + 2]
        cov_ab = fit["cov"][index:index + 2, index:index + 2]
        amplitude = float(np.hypot(a, b))
        if amplitude:
            grad_a = np.asarray([a, b]) / amplitude
            sigma_a = float(np.sqrt(max(grad_a @ cov_ab @ grad_a, 0.0)))
            grad_phi = np.asarray([b, -a]) / amplitude**2
            sigma_phi = float(np.sqrt(max(grad_phi @ cov_ab @ grad_phi, 0.0)))
        else:
            sigma_a = sigma_phi = float("nan")
        transform = np.diag([1.0, -1.0])
        cov_complex = transform @ cov_ab @ transform
        result.append({
            "k": k, "a_um": float(a), "b_um": float(b),
            "complex": complex(float(a), float(-b)), "cov_complex": cov_complex,
            "amplitude_um": amplitude, "amplitude_sigma_um": sigma_a,
            "phase_rad": float(np.arctan2(-b, a)),
            "phase_deg": float(np.degrees(np.arctan2(-b, a))),
            "phase_sigma_rad": sigma_phi,
            "phase_sigma_deg": float(np.degrees(sigma_phi)),
            "significance": amplitude / sigma_a if sigma_a else float("inf"),
        })
    return result


def make_windows(run, error_um: np.ndarray, labels: list[tuple[str, float, float]],
                 source: str) -> dict:
    windows = {}
    for name, low, high in labels:
        mask = run.stops.window == name
        target_values = np.unique(run.stops.target_mm[mask])
        rows = []
        missing = []
        for target in target_values:
            here = mask & np.isclose(run.stops.target_mm, target, atol=1e-12)
            up = error_um[here & (run.stops.direction == 1)]
            down = error_um[here & (run.stops.direction == -1)]
            if not up.size or not down.size:
                missing.append(float(target)); continue
            rows.append((target, float(up.mean()), float(down.mean())))
        if not rows:
            raise RuntimeError(f"{source} {name}: no paired targets")
        x, e_up, e_down = np.asarray(rows).T
        ybar = 0.5 * (e_up + e_down)
        reversal = e_up - e_down
        usteps = np.rint(x * 1600.0).astype(int)
        theta = 2 * np.pi * np.mod(usteps, 3200) / 3200.0
        geom_fit = linear_fit(geometry_design(x, theta, K_FIT, True), ybar)
        geom_nolinear = linear_fit(geometry_design(x, theta, K_FIT, False), ybar)
        heat_fit = linear_fit(geometry_design(x, theta, K_HEATMAP, True), ybar)
        reversal_fit = linear_fit(geometry_design(x, theta, 2, True), reversal)
        windows[name] = {
            "source": source, "name": name, "low": low, "high": high,
            "centre": float(x.mean()), "x": x, "theta": theta,
            "phase_deg": np.degrees(theta), "e_up": e_up, "e_down": e_down,
            "ybar": ybar, "B": reversal, "missing": missing,
            "fit": geom_fit, "fit_without_linear": geom_nolinear,
            "spectrum": spectrum(geom_fit, K_FIT, True),
            "spectrum_without_linear": spectrum(geom_nolinear, K_FIT, False),
            "heat_spectrum": spectrum(heat_fit, K_HEATMAP, True),
            "reversal_fit": reversal_fit,
            "reversal_spectrum": spectrum(reversal_fit, 2, True),
        }
    return windows


def complex_linear_fit(basis: np.ndarray, values: np.ndarray,
                       covariances: list[np.ndarray], name: str) -> dict:
    n, m = basis.shape
    design = np.zeros((2 * n, 2 * m))
    observed = np.empty(2 * n)
    for i in range(n):
        design[2*i, 0::2] = basis[i]
        design[2*i + 1, 1::2] = basis[i]
        observed[2*i:2*i + 2] = (values[i].real, values[i].imag)
    weight = block_diag(*[np.linalg.pinv(cov) for cov in covariances])
    normal = design.T @ weight @ design
    beta = np.linalg.pinv(normal) @ design.T @ weight @ observed
    residual = observed - design @ beta
    chi2 = float(residual @ weight @ residual)
    covariance = np.linalg.pinv(normal)
    parameters = np.asarray([complex(beta[2*j], beta[2*j + 1]) for j in range(m)])
    return {"name": name, "basis": basis, "parameters": parameters,
            "beta": beta, "cov": covariance, "chi2": chi2,
            "bic": chi2 + beta.size * math.log(2 * n), "parameter_count": beta.size}


def profile_models(x: np.ndarray, values: np.ndarray,
                   covariances: list[np.ndarray]) -> dict:
    xbar = float(x.mean())
    models = {}
    models["P0 constant"] = complex_linear_fit(np.ones((x.size, 1)), values,
                                                covariances, "P0 constant")
    models["P1 linear"] = complex_linear_fit(
        np.column_stack((np.ones(x.size), x - xbar)), values, covariances, "P1 linear")

    def p2_at(decay: float) -> dict:
        basis = np.column_stack((np.ones(x.size),
                                 np.exp(-np.abs(x - MOTOR_END_MM) / decay)))
        return complex_linear_fit(basis, values, covariances, "P2 decay")

    optimization = minimize_scalar(lambda log_l: p2_at(math.exp(log_l))["chi2"],
                                   bounds=(math.log(1.0), math.log(1000.0)),
                                   method="bounded")
    decay = float(math.exp(optimization.x))
    initial = p2_at(decay)
    whiteners = [np.linalg.cholesky(np.linalg.pinv(cov)) for cov in covariances]

    def nonlinear_residual(params: np.ndarray) -> np.ndarray:
        cinf = complex(params[0], params[1]); c0 = complex(params[2], params[3])
        lam = math.exp(params[4])
        prediction = cinf + c0 * np.exp(-np.abs(x - MOTOR_END_MM) / lam)
        blocks = []
        for i, delta in enumerate(values - prediction):
            blocks.append(whiteners[i].T @ np.asarray([delta.real, delta.imag]))
        return np.concatenate(blocks)

    fit_nl = least_squares(nonlinear_residual,
                           np.r_[initial["beta"], math.log(decay)], method="trf",
                           bounds=(np.r_[[-np.inf] * 4, math.log(1.0)],
                                   np.r_[[np.inf] * 4, math.log(1000.0)]))
    cov_nl = np.linalg.pinv(fit_nl.jac.T @ fit_nl.jac)
    p = fit_nl.x
    p2 = initial | {
        "beta": p[:4], "parameters": np.asarray([complex(p[0], p[1]),
                                                   complex(p[2], p[3])]),
        "chi2": float(fit_nl.fun @ fit_nl.fun),
        "bic": float(fit_nl.fun @ fit_nl.fun + 5 * math.log(2 * x.size)),
        "parameter_count": 5, "decay_mm": float(math.exp(p[4])),
        "decay_sigma_mm": float(math.exp(p[4]) * math.sqrt(max(cov_nl[4, 4], 0))),
        "nonlinear_cov": cov_nl,
    }
    models["P2 decay"] = p2

    knot = float(np.median(x))
    knots = np.r_[[x.min()] * 3, knot, [x.max()] * 3]
    spline_basis = BSpline.design_matrix(x, knots, 2).toarray()
    spline = complex_linear_fit(spline_basis, values, covariances,
                                "P3 quadratic spline")
    spline["knots"] = knots
    models["P3 quadratic spline"] = spline
    best_bic = min(model["bic"] for model in models.values())
    if models["P2 decay"]["bic"] <= models["P3 quadratic spline"]["bic"] and \
            models["P2 decay"]["bic"] <= best_bic + 2.0:
        chosen = "P2 decay"
    else:
        chosen = min(models, key=lambda key: models[key]["bic"])
    models["chosen"] = chosen
    models["xbar"] = xbar
    return models


def profile_predict(model: dict, x: np.ndarray, xbar: float) -> np.ndarray:
    if model["name"] == "P0 constant":
        return np.full(x.size, model["parameters"][0], dtype=complex)
    if model["name"] == "P1 linear":
        return model["parameters"][0] + model["parameters"][1] * (x - xbar)
    if model["name"] == "P2 decay":
        return model["parameters"][0] + model["parameters"][1] * np.exp(
            -np.abs(x - MOTOR_END_MM) / model["decay_mm"])
    basis = BSpline.design_matrix(x, model["knots"], 2, extrapolate=True).toarray()
    return basis @ model["parameters"]


def window_colours(windows: dict) -> dict:
    # Dark at the low-coordinate motor end, light at the support end.
    return {name: POSITION_CMAP(POSITION_NORM(window["centre"]))
            for name, window in windows.items()}


def figure_f1(windows: dict, gain_slope: float, gain_centre: float,
              colours: dict) -> None:
    corrected = {}
    maximum = 0.0
    for name, w in windows.items():
        gain = gain_slope * (w["x"] - gain_centre)
        offset = float(np.mean(w["ybar"] - gain))
        corrected[name] = {key: w[key] - gain - offset for key in ("e_up", "e_down", "ybar")}
        maximum = max(maximum, *(float(np.max(np.abs(v))) for v in corrected[name].values()))
    limit = math.ceil(maximum * 1.08 * 2) / 2
    fig, axes = plt.subplots(4, 2, figsize=(8.25, 9.4), sharey=True)
    fig.subplots_adjust(left=0.10, right=0.97, top=0.97, bottom=0.095,
                        hspace=0.36, wspace=0.17)
    for panel, (ax, (name, w)) in enumerate(zip(axes.flat, windows.items())):
        data = corrected[name]; colour = colours[name]
        ax.scatter(w["x"], data["e_up"], s=7, marker="^", color=colour,
                   alpha=0.22, linewidths=0)
        ax.scatter(w["x"], data["e_down"], s=7, marker="v", color=colour,
                   alpha=0.22, linewidths=0)
        ax.scatter(w["x"], data["ybar"], s=13, facecolor="white", edgecolor=colour,
                   linewidth=0.65, zorder=3)
        ax.axhline(0, color="#8c949a", lw=0.55)
        ax.set(xlim=(w["low"] - .15, w["high"] + .15), ylim=(-limit, limit))
        ax.set_title(f"({chr(97 + panel)}) {name} · {w['low']:.0f}–{w['high']:.0f} mm",
                     loc="left", pad=4)
        ax.grid(color=GRID, lw=.45); ax.spines[["top", "right"]].set_visible(False)
        if panel % 2 == 0: ax.set_ylabel("corrected deviation [µm]")
        if panel >= 6: ax.set_xlabel("commanded position [mm]")
    axes.flat[-1].axis("off")
    handles = [Line2D([], [], marker="^", ls="", color="#666", alpha=.35,
                      label=r"$e_\uparrow$"),
               Line2D([], [], marker="v", ls="", color="#666", alpha=.35,
                      label=r"$e_\downarrow$"),
               Line2D([], [], marker="o", ls="", mfc="white", mec="#444",
                      label=r"$\bar y$")]
    axes.flat[-1].legend(handles=handles, loc="center", frameon=False, ncol=1)
    fig.text(.10, .025,
             f"Figure F1. Directional observations and paired geometry after removal of the Step A gain "
             f"({gain_slope:+.5f} µm/mm) and each window offset.\nCommon ordinate; markers are measured data.",
             fontsize=7.4, color=INK2)
    save(fig, "e4_f1_raw_windows")


def fitted_periodic(w: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0, 2*np.pi, 721)
    if key == "geometry":
        beta = w["fit"]["beta"]; start = 2; harmonics = 2
    else:
        beta = w["reversal_fit"]["beta"]; start = 2; harmonics = 2
    values = np.full(theta.size, beta[0] if key == "reversal" else 0.0)
    for k in range(1, harmonics + 1):
        values += beta[start + 2*(k-1)] * np.cos(k*theta)
        values += beta[start + 2*(k-1) + 1] * np.sin(k*theta)
    return np.degrees(theta), values


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(color=GRID, lw=.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=.6)


def endpoint_shading(ax: plt.Axes, label: bool = False) -> None:
    """Light bands at the physical ends; text only where asked."""
    for x in (MOTOR_END_MM, SUPPORT_BEARING_MM):
        ax.axvspan(x - 1.4, x + 1.4, color="#000000", alpha=.06, lw=0, zorder=0)
    if label:
        transform = ax.get_xaxis_transform()
        ax.text(MOTOR_END_MM + 1.8, 1.02, "motor end", transform=transform,
                ha="left", va="bottom", fontsize=6.6, color=INK2, clip_on=False)
        ax.text(SUPPORT_BEARING_MM - 1.8, 1.02, "support bearing", transform=transform,
                ha="right", va="bottom", fontsize=6.6, color=INK2, clip_on=False)


def spread_labels(values: list[float], gap: float) -> list[float]:
    """Nudge label heights apart by at least `gap`, keeping their order."""
    order = np.argsort(values)
    placed = np.asarray(values, dtype=float)[order]
    for i in range(1, placed.size):
        placed[i] = max(placed[i], placed[i - 1] + gap)
    shift = np.mean(placed - np.asarray(values, dtype=float)[order])
    placed -= shift
    out = np.empty_like(placed)
    out[order] = placed
    return out.tolist()


# Line styles carry window identity, so the figure survives greyscale print.
FOLD_STYLES = {"W1": dict(color=INK, lw=1.9, ls="-"),
               "W4": dict(color=INK, lw=1.5, ls=(0, (5.5, 2.2))),
               "W7": dict(color=INK, lw=1.3, ls=(0, (1.1, 1.7)))}


def figure_f2(windows: dict, colours: dict) -> None:
    """Presentation figure 1: geometry folded onto one screw revolution."""
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    fig.subplots_adjust(left=.10, right=.86, top=.90, bottom=.20)
    curves = {name: fitted_periodic(w, "geometry") for name, w in windows.items()}
    others = [name for name in windows if name not in FOLD_STYLES]
    for name in others:
        phase, curve = curves[name]
        ax.plot(phase, curve, color="#c9cdd0", lw=.9, zorder=1,
                solid_capstyle="round")
    ends, texts = [], []
    for name, style_kw in FOLD_STYLES.items():
        phase, curve = curves[name]
        ax.plot(phase, curve, zorder=3, dash_capstyle="round",
                solid_capstyle="round", **style_kw)
        ends.append(float(curve[-1]))
        texts.append(f"{windows[name]['centre']:.0f} mm")
    low = min(windows[n]["centre"] for n in others)
    high = max(windows[n]["centre"] for n in others)
    placed = spread_labels(ends, gap=.32)
    for y_end, y_text, text in zip(ends, placed, texts):
        ax.annotate(text, (360, y_end), xytext=(372, y_text), textcoords="data",
                    va="center", ha="left", fontsize=7.6,
                    color=INK, annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=.5, color="#9aa1a6",
                                    shrinkA=0, shrinkB=1))
    ax.axhline(0, color="#8c949a", lw=.6, zorder=0)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_xlabel("screw angle θ [°]")
    ax.set_ylabel("geometric error [µm]")
    ax.set_title("Geometric error over one screw revolution, by window position",
                 loc="left", pad=6)
    clean_axes(ax)
    fig.text(.10, .045,
             "Figure 1. Direction-averaged k = 1–2 fit of each E4 window on the absolute screw angle. "
             "Line style identifies\nthe window centre; the four intermediate windows "
             f"({low:.0f}–{high:.0f} mm) are drawn in light grey.",
             fontsize=7.3, color=INK2, linespacing=1.3)
    save(fig, "e4_figure1_geometry_fold", PRESENTATION)


def figure_f2_animation(windows: dict, fps: int = 12, hold_s: float = 1.5,
                        fade_s: float = 0.9) -> None:
    """Figure 1 as a GIF: the highlight walks along the travel, window by window."""
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.colors import to_rgb

    order = sorted(windows, key=lambda name: windows[name]["centre"])
    curves = {name: fitted_periodic(windows[name], "geometry") for name in order}
    grey, ink = np.asarray(to_rgb("#c9cdd0")), np.asarray(to_rgb(INK))
    hold, fade = round(hold_s * fps), round(fade_s * fps)
    per_window = hold + fade

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    fig.subplots_adjust(left=.10, right=.97, top=.90, bottom=.20)
    lines = {}
    for name in order:
        phase, curve = curves[name]
        lines[name], = ax.plot(phase, curve, color=grey, lw=.9,
                               solid_capstyle="round")
    lows = min(float(c.min()) for _, c in curves.values())
    highs = max(float(c.max()) for _, c in curves.values())
    ax.set_ylim(lows - .25, highs + .35)
    ax.axhline(0, color="#8c949a", lw=.6, zorder=0)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_xlabel("screw angle θ [°]")
    ax.set_ylabel("geometric error [µm]")
    ax.set_title("Geometric error over one screw revolution, by window position",
                 loc="left", pad=6)
    clean_axes(ax)
    fig.text(.10, .045,
             "Direction-averaged k = 1–2 fit of each E4 window on the absolute screw angle. "
             "The highlighted window\nsteps along the travel from the motor end to the support bearing; "
             "the others stay in light grey.",
             fontsize=7.3, color=INK2, linespacing=1.3)

    # One label: the old window fades out, then the new one fades in.
    big = ax.text(.025, .955, "", transform=ax.transAxes, ha="left", va="top",
                  fontsize=15, weight="bold", color=INK)
    small = ax.text(.025, .835, "", transform=ax.transAxes, ha="left", va="top",
                    fontsize=7.6, color=INK2)

    # Travel strip: where the highlighted window sits between the two ends.
    strip = ax.inset_axes([.025, .69, .30, .07])
    strip.set_xlim(MOTOR_END_MM - 3, SUPPORT_BEARING_MM + 3)
    strip.set_ylim(-1, 1)
    strip.axis("off")
    strip.plot([MOTOR_END_MM, SUPPORT_BEARING_MM], [0, 0], color="#9aa1a6", lw=.8)
    for x in (MOTOR_END_MM, SUPPORT_BEARING_MM):
        strip.plot([x, x], [-.45, .45], color="#9aa1a6", lw=.8)
    strip.text(MOTOR_END_MM, -1.0, "motor end", ha="left", va="top",
               fontsize=6.0, color=INK2)
    strip.text(SUPPORT_BEARING_MM, -1.0, "support bearing", ha="right", va="top",
               fontsize=6.0, color=INK2)
    dots = {name: strip.plot(windows[name]["centre"], 0, "o", ms=4.2, mfc="white",
                             mec="#9aa1a6", mew=.8, zorder=3)[0] for name in order}

    def weights(frame: int) -> dict:
        index, step = divmod(frame, per_window)
        current, upcoming = order[index], order[(index + 1) % len(order)]
        level = {name: 0.0 for name in order}
        if step < hold:
            level[current] = 1.0
        else:
            t = (step - hold + 1) / fade
            t = t * t * (3 - 2 * t)              # smoothstep
            level[current], level[upcoming] = 1.0 - t, t
        return level

    def update(frame: int):
        level = weights(frame)
        for name in order:
            w = level[name]
            lines[name].set_color(grey + (ink - grey) * w)
            lines[name].set_linewidth(.9 + 1.0 * w)
            lines[name].set_zorder(1 + 2 * w)
            dots[name].set_markerfacecolor(tuple(np.asarray(to_rgb("white"))
                                                 + (ink - 1.0) * w))
            dots[name].set_markeredgecolor(tuple(np.asarray(to_rgb("#9aa1a6"))
                                                 + (ink - np.asarray(to_rgb("#9aa1a6"))) * w))
        shown = [n for n in order if level[n] > 0]
        name = max(shown, key=lambda n: level[n])
        alpha = 1.0 if len(shown) == 1 else max(0.0, 2 * level[name] - 1)
        w = windows[name]
        big.set_text(f"{w['centre']:.0f} mm")
        small.set_text(f"window {name} · {w['low']:.0f}–{w['high']:.0f} mm")
        big.set_alpha(alpha); small.set_alpha(alpha)
        return []

    PRESENTATION.mkdir(parents=True, exist_ok=True)
    animation = FuncAnimation(fig, update, frames=per_window * len(order),
                              blit=False)
    animation.save(PRESENTATION / "e4_figure1_geometry_fold.gif",
                   writer=PillowWriter(fps=fps), dpi=150)
    plt.close(fig)


def collect_profile(e1_windows: dict, e4_windows: dict, k: int):
    entries = [(w["centre"], w["spectrum"][k-1], "E1", name) for name, w in e1_windows.items()]
    entries += [(w["centre"], w["spectrum"][k-1], "E4", name) for name, w in e4_windows.items()]
    entries.sort(key=lambda item: item[0])
    x = np.asarray([item[0] for item in entries])
    values = np.asarray([item[1]["complex"] for item in entries])
    covs = [item[1]["cov_complex"] for item in entries]
    return entries, x, values, covs


E1_NUDGE_MM = 1.4   # E1 markers sit just right of coincident E4 windows


def profile_marker(ax, x, y, yerr, source, weak=False):
    ink = "#a3a9ad" if weak else INK
    ax.errorbar(x + (E1_NUDGE_MM if source == "E1" else 0.0), y, yerr=yerr,
                fmt="o" if source == "E4" else "s", ms=4.6 if source == "E4" else 4.3,
                mfc=("white" if source == "E1" or weak else ink), mec=ink, mew=.85,
                ecolor=ink, elinewidth=.75, capsize=1.8, capthick=.75, zorder=3,
                ls="none")


def figure_f3(e1_windows: dict, e4_windows: dict, profiles: dict) -> None:
    """Presentation figure 2: k = 1/2 amplitude and phase profiles, one per panel."""
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.9), sharex=True)
    fig.subplots_adjust(left=.095, right=.98, top=.87, bottom=.145,
                        hspace=.22, wspace=.24)
    amplitude_top = 0.0
    weak_any = False
    for col, k in enumerate((1, 2)):
        entries, x, _, _ = collect_profile(e1_windows, e4_windows, k)
        chosen = profiles[k][profiles[k]["chosen"]]
        constant = profiles[k]["P0 constant"]
        grid = np.linspace(float(x.min()), float(x.max()), 600)
        ax = axes[0, col]
        ax.plot(grid, np.abs(profile_predict(constant, grid, profiles[k]["xbar"])),
                color="#8c949a", lw=.9, ls=(0, (4, 2.5)), zorder=1)
        ax.plot(grid, np.abs(profile_predict(chosen, grid, profiles[k]["xbar"])),
                color=INK, lw=1.3, zorder=2)
        for xpos, spec, source, _ in entries:
            profile_marker(ax, xpos, spec["amplitude_um"],
                           2*spec["amplitude_sigma_um"], source)
            amplitude_top = max(amplitude_top,
                                spec["amplitude_um"] + 2*spec["amplitude_sigma_um"])

        raw_phase = np.asarray([item[1]["phase_rad"] for item in entries])
        unwrapped_deg = np.degrees(np.unwrap(raw_phase))
        reliable = np.asarray([(item[1]["significance"] >= 3.0 and
                                (k == 1 or item[0] < 105.0)) for item in entries])
        phase_sigma = np.asarray([item[1]["phase_sigma_deg"] for item in entries])
        weights = 1.0 / np.maximum(phase_sigma[reliable], 1e-9)**2
        mean_phase = float(np.average(unwrapped_deg[reliable], weights=weights))
        screw_shift = -(unwrapped_deg - mean_phase) / k
        ax = axes[1, col]
        ax.axhline(0, color="#8c949a", lw=.6, zorder=0)
        for index, (xpos, spec, source, _) in enumerate(entries):
            weak = not bool(reliable[index])
            weak_any |= weak
            profile_marker(ax, xpos, screw_shift[index],
                           2*spec["phase_sigma_deg"] / k, source, weak=weak)
            if weak and source == "E4":
                ax.annotate("weakly\nidentified", (xpos, screw_shift[index]),
                            xytext=(-30, -4), textcoords="offset points",
                            ha="right", va="center", fontsize=6.4, color="#7d858c",
                            arrowprops=dict(arrowstyle="-", lw=.5, color="#a3a9ad"))
        phase_limit = float(x.max()) if k == 1 else float(x[reliable].max())
        phase_grid = np.linspace(float(x.min()), phase_limit, 500)
        model_phase = np.unwrap(np.angle(
            profile_predict(chosen, phase_grid, profiles[k]["xbar"])))
        shift = 360 * round((np.median(unwrapped_deg[reliable])
                             - np.median(np.degrees(model_phase))) / 360)
        ax.plot(phase_grid, -(np.degrees(model_phase) + shift - mean_phase) / k,
                color=INK, lw=1.3, zorder=2)

        axes[0, col].set_title(
            f"k = {k}  ·  {'once' if k == 1 else 'twice'} per revolution",
            loc="left", pad=15, fontsize=9.2)
    for ax in axes[0]:
        ax.set_ylim(0, math.ceil(amplitude_top * 1.08 * 4) / 4)
    for row, ax in enumerate(axes.flat):
        endpoint_shading(ax, label=row < 2)
        ax.set_xlim(MOTOR_END_MM - 2.0, SUPPORT_BEARING_MM + 2.0)
        clean_axes(ax)
    axes[0, 0].set_ylabel(r"amplitude $A_k$ [µm]")
    axes[1, 0].set_ylabel(r"pattern shift $\Delta\theta_k$ [°]")
    for ax in axes[1]:
        ax.set_xlabel("window centre [mm]")
    handles = [
        Line2D([], [], marker="o", ls="", ms=4.6, mfc=INK, mec=INK, label="E4 window"),
        Line2D([], [], marker="s", ls="", ms=4.3, mfc="white", mec=INK, label="E1 window (refit)"),
        Line2D([], [], color=INK, lw=1.3, label="smooth profile fit"),
        Line2D([], [], color="#8c949a", lw=.9, ls=(0, (4, 2.5)), label="constant (stationary) model"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.085, .975),
               ncol=4, frameon=False, fontsize=7.4, handlelength=2.2,
               columnspacing=1.6)
    for letter, ax in zip("abcd", axes.flat):
        ax.text(.015, .97, f"({letter})", transform=ax.transAxes, ha="left",
                va="top", fontsize=8, weight="bold", color=INK)
    fig.text(.095, .035,
             "Figure 2. Harmonic amplitude (top) and pattern shift around the screw (bottom) along the travel; "
             r"error bars ±2σ, $\Delta\theta_k=-(\phi_k-\bar\phi_k)/k$." "\n"
             "E1 markers are offset 1.4 mm to the right where they coincide with E4. "
             "Shaded bands mark the physical ends, inferred from homing.",
             fontsize=7.1, color=INK2, linespacing=1.3)
    save(fig, "e4_figure2_harmonic_profile", PRESENTATION)


def reversal_metrics(w: dict) -> dict:
    fit = w["reversal_fit"]
    mean = float(fit["beta"][0]); mean_sigma = float(np.sqrt(fit["cov"][0, 0]))
    harmonic = w["reversal_spectrum"][0]
    theta = np.linspace(0, 2*np.pi, 1000)
    g = w["spectrum"]
    geom_derivative = (-g[0]["a_um"]*np.sin(theta) + g[0]["b_um"]*np.cos(theta)
                       - 2*g[1]["a_um"]*np.sin(2*theta) + 2*g[1]["b_um"]*np.cos(2*theta))
    rb = w["reversal_spectrum"]
    reversal = (rb[0]["a_um"]*np.cos(theta) + rb[0]["b_um"]*np.sin(theta)
                + rb[1]["a_um"]*np.cos(2*theta) + rb[1]["b_um"]*np.sin(2*theta))
    correlation = float(np.corrcoef(geom_derivative, reversal)[0, 1])
    return {"mean_um": mean, "mean_sigma_um": mean_sigma,
            "amplitude_um": harmonic["amplitude_um"],
            "amplitude_sigma_um": harmonic["amplitude_sigma_um"],
            "phase_deg": harmonic["phase_deg"],
            "derivative_correlation_r": correlation}


def figure_f4(windows: dict, metrics: dict,
              e1_windows: dict, e1_metrics: dict) -> None:
    """Presentation figure 3: reversal profile and direct co-variation test."""
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.9),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    fig.subplots_adjust(left=.085, right=.975, top=.82, bottom=.24, wspace=.30)
    centres = np.asarray([w["centre"] for w in windows.values()])
    mean_b = np.asarray([metrics[name]["mean_um"] for name in windows])
    mean_sigma = np.asarray([metrics[name]["mean_sigma_um"] for name in windows])
    geometry = np.asarray([windows[name]["spectrum"][0]["amplitude_um"] for name in windows])
    geometry_sigma = np.asarray([windows[name]["spectrum"][0]["amplitude_sigma_um"]
                                 for name in windows])
    reversal = np.asarray([metrics[name]["amplitude_um"] for name in windows])
    reversal_sigma = np.asarray([metrics[name]["amplitude_sigma_um"] for name in windows])
    e1_centres = np.asarray([w["centre"] for w in e1_windows.values()])
    e1_mean = np.asarray([e1_metrics[name]["mean_um"] for name in e1_windows])
    e1_mean_sigma = np.asarray([e1_metrics[name]["mean_sigma_um"] for name in e1_windows])

    ax = axes[0]
    all_window_mean = float(np.mean(mean_b))
    ax.axhline(all_window_mean, color="#8c949a", lw=.9, ls=(0, (4, 2.5)), zorder=1)
    ax.text(66, all_window_mean - .012,
            f"E4 mean {all_window_mean:.2f} µm".replace("-", "\u2212"),
            ha="left", va="top", fontsize=6.8, color=INK2)
    for xpos, y, s in zip(centres, mean_b, mean_sigma):
        profile_marker(ax, xpos, y, 2*s, "E4")
    for xpos, y, s in zip(e1_centres, e1_mean, e1_mean_sigma):
        profile_marker(ax, xpos, y, 2*s, "E1")
    endpoint_shading(ax, label=True)
    ax.set_xlim(MOTOR_END_MM - 2.0, SUPPORT_BEARING_MM + 2.0)
    span = np.concatenate([mean_b - 2*mean_sigma, mean_b + 2*mean_sigma,
                           e1_mean - 2*e1_mean_sigma, e1_mean + 2*e1_mean_sigma])
    pad = .15 * (span.max() - span.min())
    ax.set_ylim(span.min() - pad, span.max() + pad)
    ax.set_xlabel("window centre [mm]")
    ax.set_ylabel(r"mean reversal $\beta_0$ [µm]")
    ax.set_title(r"$\bf{(a)}$ Mean reversal along the travel", loc="left", pad=15)

    ax = axes[1]
    for g, gs, r, rs in zip(geometry, geometry_sigma, reversal, reversal_sigma):
        ax.errorbar(g, r, xerr=2*gs, yerr=2*rs, fmt="o", ms=4.6, mfc=INK, mec=INK,
                    ecolor="#7d858c", elinewidth=.7, capsize=1.8, capthick=.7,
                    zorder=3, ls="none")
    # Direct position labels; nudge the crowded low-amplitude cluster apart.
    label_offsets = {}
    for i in np.argsort(geometry):
        label_offsets[i] = (6, 5)
    crowded = [i for i in range(geometry.size)
               if any(j != i and abs(geometry[j] - geometry[i]) < .12 and
                      abs(reversal[j] - reversal[i]) < .1 for j in range(geometry.size))]
    for n, i in enumerate(sorted(crowded, key=lambda i: reversal[i])):
        label_offsets[i] = (6, -9) if n % 2 == 0 else (6, 6)
    for i, centre in enumerate(centres):
        ax.annotate(f"{centre:.0f} mm", (geometry[i], reversal[i]),
                    xytext=label_offsets[i], textcoords="offset points",
                    fontsize=6.8, color=INK2, ha="left", va="center")
    correlation, pvalue = pearsonr(geometry, reversal)
    ax.text(.5, .04,
            rf"Pearson $r$ = {correlation:.2f},  $n$ = {geometry.size},  $p \approx {pvalue:.2f}$",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.4, color=INK)
    ax.set_xlabel(r"geometric amplitude $A_1$ [µm]")
    ax.set_ylabel(r"reversal amplitude $A_{B,1}$ [µm]")
    ax.set_title(r"$\bf{(b)}$ Reversal against geometry, k = 1", loc="left", pad=15)
    ax.set_ylim(bottom=min(0.0, float(np.min(reversal - 2*reversal_sigma))) - .02)
    for ax in axes:
        clean_axes(ax)
    handles = [
        Line2D([], [], marker="o", ls="", ms=4.6, mfc=INK, mec=INK, label="E4 window"),
        Line2D([], [], marker="s", ls="", ms=4.3, mfc="white", mec=INK, label="E1 window (refit)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.075, .985),
               ncol=2, frameon=False, fontsize=7.4, columnspacing=1.6)
    fig.text(.085, .045,
             "Figure 3. (a) Mean directional reversal per window. (b) Once-per-revolution reversal against geometric "
             "amplitude for the\nE4 windows, labelled by window centre; no clear co-variation. "
             "Error bars ±2σ; E1 markers offset 1.4 mm where they coincide with E4.",
             fontsize=7.1, color=INK2, linespacing=1.3)
    save(fig, "e4_figure3_reversal", PRESENTATION)


def figure_validity(run, q: dict, drift: dict, windows: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.25, 7.0))
    fig.subplots_adjust(left=.10, right=.97, top=.96, bottom=.12, hspace=.38, wspace=.30)
    ax = axes[0, 0]
    ax.plot(drift["times_min"], drift["drift_from_datum_um"], color="#73578f", lw=.8)
    ax.scatter(drift["times_min"], drift["drift_from_datum_um"], s=20, facecolor="white",
               edgecolor="#73578f", linewidth=.8)
    ax.axhline(.5, color=LIMIT, ls="--", lw=.8); ax.axhline(-.5, color=LIMIT, ls="--", lw=.8)
    ax.set(xlabel="ESP experiment time [min]", ylabel="reference drift [µm]")
    ax.set_title("(a) Reference stability", loc="left", pad=4)
    ax = axes[0, 1]
    ax.hist(q["std_um"], bins=32, color="#6f8794", edgecolor="white", lw=.3)
    ax.axvline(STD_LIMIT_UM, color=LIMIT, ls="--", lw=.9)
    ax.set(xlabel="record-window standard deviation [µm]", ylabel="stop count")
    ax.set_title("(b) In-window scatter", loc="left", pad=4)
    ax = axes[1, 0]
    effect = np.abs(q["slope_um_s"]) * .5
    ax.hist(effect, bins=32, color="#9a7b4f", edgecolor="white", lw=.3)
    ax.axvline(SLOPE_EFFECT_LIMIT_UM, color=LIMIT, ls="--", lw=.9)
    ax.set(xlabel=r"$|$slope$|t_{record}$ [µm]", ylabel="stop count")
    ax.set_title("(c) Drift within records", loc="left", pad=4)
    ax = axes[1, 1]
    values = [w["fit"]["condition"] for w in windows.values()]
    ax.bar(list(windows), values,
           color=[POSITION_CMAP(POSITION_NORM(w["centre"])) for w in windows.values()])
    ax.set(xlabel="window", ylabel="design condition number")
    ax.set_title("(d) Harmonic-fit conditioning", loc="left", pad=4)
    for ax in axes.flat:
        ax.grid(axis="y", color=GRID, lw=.45); ax.spines[["top", "right"]].set_visible(False)
    fig.text(.10, .03,
             f"Appendix A1. All {run.stops.state.size} recorded stops. Std pass "
             f"{100*np.mean(q['std_um'] <= STD_LIMIT_UM):.1f}%; slope-effect pass "
             f"{100*np.mean(np.abs(q['slope_um_s'])*.5 <= SLOPE_EFFECT_LIMIT_UM):.1f}%; "
             f"saturated record increments {int(q['saturated'].sum())}; reference datum-to-end "
             f"{drift['datum_to_end_um']:+.3f} µm.", fontsize=7.3, color=INK2)
    save(fig, "e4_appendix_validity")


def figure_residuals(windows: dict, colours: dict) -> None:
    maximum = max(float(np.max(np.abs(w["fit"]["residual"]))) for w in windows.values())
    limit = math.ceil(maximum * 1.1 * 10) / 10
    fig, axes = plt.subplots(4, 2, figsize=(8.25, 8.5), sharey=True)
    fig.subplots_adjust(left=.10, right=.97, top=.96, bottom=.10, hspace=.36, wspace=.17)
    for panel, (ax, (name, w)) in enumerate(zip(axes.flat, windows.items())):
        ax.scatter(w["x"], w["fit"]["residual"], s=12, facecolor="white",
                   edgecolor=colours[name], linewidth=.65)
        ax.axhline(0, color="#8c949a", lw=.55); ax.set_ylim(-limit, limit)
        ax.set_title(f"({chr(97+panel)}) {name} · {w['low']:.0f}–{w['high']:.0f} mm",
                     loc="left", pad=4)
        ax.grid(color=GRID, lw=.45); ax.spines[["top", "right"]].set_visible(False)
        if panel % 2 == 0: ax.set_ylabel("residual [µm]")
        if panel >= 6: ax.set_xlabel("commanded position [mm]")
    axes.flat[-1].axis("off")
    fig.text(.10, .025, "Appendix A2. Residuals after the intercept + local trend + k=1…4 OLS fit. "
             "All panels share the ordinate.", fontsize=7.4, color=INK2)
    save(fig, "e4_appendix_fit_residuals")


def figure_heatmap(windows: dict) -> None:
    matrix = np.asarray([[h["amplitude_um"] for h in w["heat_spectrum"]]
                         for w in windows.values()])
    fig, ax = plt.subplots(figsize=(7.2, 3.7))
    fig.subplots_adjust(left=.12, right=.90, top=.91, bottom=.18)
    image = ax.imshow(matrix, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_xticks(np.arange(K_HEATMAP), np.arange(1, K_HEATMAP+1))
    ax.set_yticks(np.arange(len(windows)),
                  [f"{name} · {w['centre']:.0f} mm" for name, w in windows.items()])
    ax.set(xlabel="screw harmonic k", ylabel="window (motor to support)")
    ax.set_title("Appendix A3. E4 harmonic-amplitude heatmap", loc="left", pad=7)
    bar = fig.colorbar(image, ax=ax, pad=.025); bar.set_label("amplitude [µm]")
    fig.text(.12, .04, "Each row is an independent intercept + local trend + k=1…8 fit.",
             fontsize=7.4, color=INK2)
    save(fig, "e4_appendix_harmonic_heatmap")


def figure_e1_comparison(windows: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.3), sharex=True)
    fig.subplots_adjust(left=.10, right=.97, top=.95, bottom=.12, hspace=.29, wspace=.26)
    x = np.asarray([w["centre"] for w in windows.values()])
    for k, colour in HARMONIC_COLOURS.items():
        before = [w["spectrum_without_linear"][k-1] for w in windows.values()]
        after = [w["spectrum"][k-1] for w in windows.values()]
        axes[0, k-1].errorbar(x, [s["amplitude_um"] for s in before],
                              yerr=[2*s["amplitude_sigma_um"] for s in before],
                              fmt="o", mfc="white", mec=colour, ecolor=colour,
                              capsize=2, label="without local trend")
        axes[0, k-1].errorbar(x, [s["amplitude_um"] for s in after],
                              yerr=[2*s["amplitude_sigma_um"] for s in after],
                              fmt="s", mfc=colour, mec=colour, ecolor=colour,
                              capsize=2, label="with local trend")
        p0 = np.unwrap([s["phase_rad"] for s in before]); p1 = np.unwrap([s["phase_rad"] for s in after])
        axes[1, k-1].plot(x, np.degrees(p0), "o--", color=colour, mfc="white", lw=.8)
        axes[1, k-1].plot(x, np.degrees(p1), "s-", color=colour, lw=1.0)
        axes[0, k-1].set_title(f"({'a' if k==1 else 'b'}) k={k} amplitude", loc="left")
        axes[1, k-1].set_title(f"({'c' if k==1 else 'd'}) k={k} phase", loc="left")
        axes[0, k-1].set_ylabel("amplitude [µm]"); axes[1, k-1].set_ylabel("phase [deg]")
        axes[1, k-1].set_xlabel("E1 window centre [mm]")
    axes[0, 0].legend(frameon=False)
    for ax in axes.flat:
        ax.grid(color=GRID, lw=.45); ax.spines[["top", "right"]].set_visible(False)
    fig.text(.10, .03, "Appendix A4. E1 refit comparison; error bars are ±2σ. The refit uses the same "
             "intercept + local trend + k=1…4 model as E4.", fontsize=7.4, color=INK2)
    save(fig, "e4_appendix_e1_refit_comparison")


def serial_complex(value: complex) -> dict:
    return {"real_um": float(value.real), "imag_um": float(value.imag)}


def main() -> int:
    style(); ANALYSIS.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    step_a = json.loads(STEP_A.read_text(encoding="utf-8"))["fit_uncorrected"]

    run4 = e4_run.load(); q4 = quality(run4)
    error4, drift4 = drift_corrected_error(run4, q4["means_mm"])
    windows4 = make_windows(run4, error4, list(e4_run.WINDOWS), "E4")
    colours = window_colours(windows4)

    run1 = e1_run.load(); q1 = quality(run1)
    error1, _ = drift_corrected_error(run1, q1["means_mm"])
    windows1 = make_windows(run1, error1, list(e1_run.WINDOWS), "E1")

    profiles = {}
    for k in range(1, K_FIT + 1):
        _, x, values, covariances = collect_profile(windows1, windows4, k)
        profiles[k] = profile_models(x, values, covariances)

    reversal = {name: reversal_metrics(window) for name, window in windows4.items()}
    e1_reversal = {name: reversal_metrics(window) for name, window in windows1.items()}
    figure_f1(windows4, float(step_a["slope_um_per_mm"]), float(step_a["x_center_mm"]), colours)
    figure_f2(windows4, colours)
    figure_f2_animation(windows4)
    figure_f3(windows1, windows4, profiles)
    figure_f4(windows4, reversal, windows1, e1_reversal)
    figure_validity(run4, q4, drift4, windows4)
    figure_residuals(windows4, colours)
    figure_heatmap(windows4)
    figure_e1_comparison(windows1)

    rows = []
    for source, windows in (("E1", windows1), ("E4", windows4)):
        for name, w in windows.items():
            row = {"source": source, "window": name, "x_low_mm": w["low"],
                   "x_high_mm": w["high"], "x_centre_mm": w["centre"],
                   "paired_targets": len(w["x"]), "missing_targets": len(w["missing"]),
                   "slope_um_per_mm": float(w["fit"]["beta"][1]),
                   "residual_rms_um": float(np.sqrt(np.mean(w["fit"]["residual"]**2))),
                   "condition_number": w["fit"]["condition"]}
            for k in range(1, 5):
                h = w["spectrum"][k-1]
                row |= {f"A{k}_um": h["amplitude_um"], f"A{k}_sigma_um": h["amplitude_sigma_um"],
                        f"phi{k}_deg": h["phase_deg"], f"phi{k}_sigma_deg": h["phase_sigma_deg"],
                        f"A{k}_over_sigma": h["significance"]}
            rb = reversal.get(name, {}) if source == "E4" else reversal_metrics(w)
            row |= {"mean_B_um": rb["mean_um"], "mean_B_sigma_um": rb["mean_sigma_um"],
                    "A_B1_um": rb["amplitude_um"], "A_B1_sigma_um": rb["amplitude_sigma_um"]}
            rows.append(row)
    with (ANALYSIS / "e4_t1_window_harmonics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    model_rows, model_json = [], {}
    for k in range(1, K_FIT + 1):
        model_json[str(k)] = {"chosen": profiles[k]["chosen"], "models": {}}
        for name, model in profiles[k].items():
            if not isinstance(model, dict): continue
            parameters = [serial_complex(value) for value in model["parameters"]]
            record = {"harmonic": k, "model": name, "parameter_count": model["parameter_count"],
                      "chi2": model["chi2"], "bic": model["bic"],
                      "chosen": name == profiles[k]["chosen"],
                      "parameters_json": json.dumps(parameters)}
            if name == "P2 decay":
                cinf, c0 = model["parameters"]
                record |= {"lambda_mm": model["decay_mm"],
                           "lambda_sigma_mm": model["decay_sigma_mm"],
                           "C0_amplitude_um": abs(c0), "Cinf_amplitude_um": abs(cinf),
                           "phase_offset_deg": float(np.degrees(np.angle(c0/cinf)))}
            model_rows.append(record)
            model_json[str(k)]["models"][name] = record
    fields = sorted({key for row in model_rows for key in row})
    with (ANALYSIS / "e4_t2_profile_models.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(model_rows)

    ref = run4.stops.window == "REF"
    reference_mscnt = run4.stops.mscnt[ref]
    validity = [
        {"check": "analysis stop count", "value": int(np.count_nonzero(run4.stops.window != "REF")),
         "criterion": 714, "pass": int(np.count_nonzero(run4.stops.window != "REF")) == 714},
        {"check": "paired targets", "value": sum(len(w["x"]) for w in windows4.values()),
         "criterion": 357, "pass": all(not w["missing"] for w in windows4.values())},
        {"check": "record std", "value": float(np.mean(q4["std_um"] <= STD_LIMIT_UM)),
         "criterion": "report pass fraction; limit 0.5 um", "pass": bool(np.all(q4["std_um"] <= STD_LIMIT_UM))},
        {"check": "record slope effect", "value": float(np.mean(np.abs(q4["slope_um_s"])*config["record_s"] <= SLOPE_EFFECT_LIMIT_UM)),
         "criterion": "report pass fraction; limit 0.25 um", "pass": bool(np.all(np.abs(q4["slope_um_s"])*config["record_s"] <= SLOPE_EFFECT_LIMIT_UM))},
        {"check": "saturated record samples", "value": int(q4["saturated"].sum()),
         "criterion": 0, "pass": int(q4["saturated"].sum()) == 0},
        {"check": "reference MSCNT consistency", "value": reference_mscnt.tolist(),
         "criterion": "all reference visits identical", "pass": bool(np.all(reference_mscnt == reference_mscnt[0]))},
        {"check": "reference drift warning", "value": drift4["datum_to_end_um"],
         "criterion": "abs <= 0.5 um", "pass": abs(drift4["datum_to_end_um"]) <= .5},
        {"check": "environment", "value": "snapshot requested; no values logged",
         "criterion": "report availability", "pass": None},
    ]
    with (ANALYSIS / "e4_validity_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "value", "criterion", "pass"])
        writer.writeheader(); writer.writerows(validity)

    phase_statement = {}
    for k in (1, 2):
        entries, _, _, _ = collect_profile(windows1, windows4, k)
        phase = np.unwrap([entry[1]["phase_rad"] for entry in entries])
        jumps = np.abs(np.diff(np.degrees(phase)))
        phase_statement[str(k)] = {"maximum_adjacent_jump_deg": float(jumps.max()),
                                   "smooth_under_90_deg": bool(jumps.max() < 90)}
    reversal_amplitudes = np.asarray([reversal[name]["amplitude_um"] for name in windows4])
    geometry_amplitudes = np.asarray([windows4[name]["spectrum"][0]["amplitude_um"]
                                      for name in windows4])
    amplitude_correlation = float(np.corrcoef(reversal_amplitudes,
                                               geometry_amplitudes)[0, 1])
    higher_harmonics = {}
    for k in (3, 4):
        detections = []
        for source, windows in (("E1", windows1), ("E4", windows4)):
            for name, window in windows.items():
                h = window["spectrum"][k - 1]
                if h["significance"] >= 3.0:
                    detections.append({"source": source, "window": name,
                                       "x_centre_mm": window["centre"],
                                       "significance": h["significance"]})
        higher_harmonics[str(k)] = {"detections_at_or_above_3sigma": detections,
                                    "include_in_model_c": bool(detections)}
    summary = {
        "experiment": "E4_stationarity_windows", "source_stops": int(run4.stops.state.size),
        "analysis_stops": int(np.count_nonzero(run4.stops.window != "REF")),
        "scope_samples": int(run4.t.size), "alignment_fit_mad_um": run4.fit_mad_um,
        "esp_clock_rate_ppm": run4.rate_ppm,
        "endpoint_assumptions": {"motor_end_mm": MOTOR_END_MM,
                                 "support_bearing_mm": SUPPORT_BEARING_MM,
                                 "status": ENDPOINT_STATUS},
        "drift": {key: value for key, value in drift4.items() if not isinstance(value, np.ndarray)},
        "validity": validity, "profile_models": model_json,
        "phase_continuity": phase_statement,
        "higher_harmonic_screen": higher_harmonics,
        "reversal_geometry_exploratory": {
            "per_window": reversal,
            "A_B1_vs_A1_pearson_r": amplitude_correlation,
            "interpretation": "exploratory; require a consistent sign across windows",
        },
        "model_c_recommendation": {
            f"k{k}": {"model": profiles[k]["chosen"],
                       "term": f"Re[C_{k}(x) exp(i {k} theta)]",
                       "status": ("core" if k <= 2 else
                                  "retain: at least one E1/E4 window is >=3 sigma")}
            for k in range(1, K_FIT + 1)},
        "convention": {"deviation_sign": "positive means measured travel exceeds command",
                       "theta": "2 pi (commanded microsteps mod 3200) / 3200",
                       "phase": "e = A cos(k theta + phi)",
                       "drift": "linear interpolation through fixed-position reference visits"},
    }
    (ANALYSIS / "e4_findings_summary.json").write_text(
        json.dumps(summary, indent=2, default=lambda value: bool(value) if isinstance(value, np.bool_) else float(value)) + "\n",
        encoding="utf-8")

    print(f"E4: {sum(len(w['x']) for w in windows4.values())} paired targets; "
          f"std pass {100*np.mean(q4['std_um'] <= STD_LIMIT_UM):.1f}%; "
          f"slope pass {100*np.mean(np.abs(q4['slope_um_s'])*config['record_s'] <= SLOPE_EFFECT_LIMIT_UM):.1f}%")
    print(f"reference drift {drift4['datum_to_end_um']:+.3f} um; "
          f"MSCNT values {np.unique(reference_mscnt).tolist()}")
    for k in range(1, K_FIT + 1):
        chosen = profiles[k][profiles[k]["chosen"]]
        print(f"k={k}: chosen {profiles[k]['chosen']}, BIC={chosen['bic']:.2f}")
    for path in sorted(OUT.glob("e4_*")):
        print("wrote", path.relative_to(HERE))
    for path in sorted(PRESENTATION.glob("e4_*")):
        print("wrote", path.relative_to(HERE))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
