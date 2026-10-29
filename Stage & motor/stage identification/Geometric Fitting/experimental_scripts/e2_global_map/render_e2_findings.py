#!/usr/bin/env python3
"""Process E2 and render the global geometric-map findings.

The existing E2 figures remain acquisition/design diagnostics.  This script is
the inferential analysis requested after E1/E4: three drift variants, blocked
cross-validation of stationary and varying-coefficient models, independent
E1/E4/E0 validation, reversal analysis, and publication-style PNG figures.
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.ticker
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import least_squares
from scipy.signal import lombscargle

import e2_run

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"
OUT = HERE / "rendered" / "experimental_findings"
PRESENTATION = HERE / "rendered" / "presentation"
CONFIG = HERE / "e2_config.json"
LOG = HERE / "data" / "e2_terminal_20261027_151332.log"
# The v3 session is read through its corrected evaluation. The original
# iso230_2_v3_results.json placed each averaging window with a fitted clock
# whose stop detector ran late, so every window closed ~38 ms into the next
# move; that biases xbar_up and xbar_dn by -+0.62 um and nothing else. The
# direction average is untouched, but the reversal, R, E and A that this
# script validates against are not. Both paths stay overridable so the old
# inputs can still be reproduced for comparison.
STEP_A = Path(os.environ.get(
    "E2_STEP_A", HERE.parents[1] / "results_v3_corrected" / "step_a_results.json"))
ISO_RESULTS = Path(os.environ.get(
    "E2_ISO_RESULTS",
    HERE.parents[2] / "iso230-2" / "analysis" / "iso230_2_v3_results_blind.json"))
E1_DIR = HERE.parent / "e1_spectral_windows"
E4_DIR = HERE.parent / "e4_stationarity_windows"
for directory in (E1_DIR, E4_DIR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
import e1_run  # noqa: E402
import e4_run  # noqa: E402
import render_e4_findings as e4f  # noqa: E402

CV_OFFSETS = (0.0, 0.2, 0.4, 0.6, 0.8)   # block-edge shifts, in block widths
CV_MIN_TEST = 5                          # a shifted end block can be narrow
PEN_KNOTS = 8                            # deliberately rich; lambda does the work
PEN_ORDER = 2                            # second differences between neighbours
PEN_LAMBDAS = (1e-2, 3.2e-2, 1e-1, 3.2e-1, 1.0, 3.2, 1e1, 3.2e1, 1e2)
CV_END_FACTOR = 5.0                      # held-out end block, over interior CV
MOTOR_END_MM = 0.0
SUPPORT_BEARING_MM = 122.770
STD_LIMIT_UM = 0.5
SLOPE_EFFECT_LIMIT_UM = 0.25
SATURATION_INCREMENT_NM = 8333.0
INK = "#252b30"
INK2 = "#59636b"
GRID = "#dfe3e6"
BLUE = "#315b7d"
RUST = "#a3553d"
GREEN = "#2f6b4f"
GREY = "#8b949a"
LIGHT = "#c7ccd0"
DRIFT_COLOURS = {"D0": "#9a6b3f", "D1": "#657a8a", "D2": "#2f6b4f"}
POSITION_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "e2_position", ("#243b62", "#24677a", "#2f8a7a", "#9a7b3f"))
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
        "legend.fontsize": 7.1, "axes.unicode_minus": True,
    })


def save(fig: plt.Figure, stem: str, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{stem}.png", dpi=320)
    plt.close(fig)


def finish_axis(ax: plt.Axes) -> None:
    ax.grid(color=GRID, lw=.45)
    ax.spines[["top", "right"]].set_visible(False)


def endpoint_bands(ax: plt.Axes, labels: bool = True) -> None:
    ax.axvspan(MOTOR_END_MM - 1.4, MOTOR_END_MM + 1.4,
               color="#725247", alpha=.09, lw=0, zorder=0)
    ax.axvspan(SUPPORT_BEARING_MM - 1.4, SUPPORT_BEARING_MM + 1.4,
               color="#526d72", alpha=.09, lw=0, zorder=0)
    if labels:
        transform = ax.get_xaxis_transform()
        ax.text(MOTOR_END_MM + 1.8, 1.02, "motor end", transform=transform,
                ha="left", va="bottom", fontsize=6.8, color=INK2, clip_on=False)
        ax.text(SUPPORT_BEARING_MM - 1.8, 1.02, "support bearing", transform=transform,
                ha="right", va="bottom", fontsize=6.8, color=INK2, clip_on=False)


def panel_title(ax: plt.Axes, label: str, title: str, raised: bool = False) -> None:
    ax.set_title(rf"$\bf{{({label})}}$ {title}", loc="left",
                 y=1.055 if raised else 1.0, pad=4)


def record_quality(run: e2_run.Run) -> dict:
    means = run.hold()
    std, slope, saturated, samples, time = [], [], [], [], []
    for i in range(run.plan.target_mm.size):
        lo, hi = run.record_window(i)
        segment = run.pos[lo:hi] * 1000.0
        seconds = run.t[lo:hi]
        centred = seconds - seconds.mean()
        std.append(float(segment.std(ddof=1)))
        slope.append(float(np.dot(centred, segment - segment.mean()) /
                           np.dot(centred, centred)))
        saturated.append(int(np.count_nonzero(
            np.abs(np.diff(run.pos[lo:hi])) * 1e6 >= SATURATION_INCREMENT_NM)))
        samples.append(int(hi - lo))
        time.append(float(0.5 * (seconds[0] + seconds[-1])))
    return {"means_mm": means, "std_um": np.asarray(std),
            "slope_um_s": np.asarray(slope), "saturated": np.asarray(saturated),
            "samples": np.asarray(samples), "time_s": np.asarray(time)}


def parse_serial_checks() -> dict:
    stop_re = re.compile(
        r"# STOP,([A-Z_0-9]+),([A-Z_0-9]+),(\d+),(-?\d+),(\d+),(\d+),"
        r"([\d.]+),(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)")
    states, kinds, usteps, mscnt = [], [], [], []
    summary = None
    for line in LOG.read_text(encoding="utf-8", errors="replace").splitlines():
        match = stop_re.search(line)
        if match:
            states.append(match.group(1)); kinds.append(match.group(2))
            usteps.append(int(match.group(6))); mscnt.append(int(match.group(13)))
        if "# SUMMARY," in line:
            summary = line.split("# SUMMARY,", 1)[1]
    states = np.asarray(states); kinds = np.asarray(kinds)
    usteps = np.asarray(usteps); mscnt = np.asarray(mscnt)
    ref = kinds == "REF"
    expected = np.mod(usteps, 3200)
    # MSCNT is the driver's 0..1023 phase; reference equality is the decisive
    # lost-step check because every reference uses the same absolute command.
    return {"stop_count": int(states.size), "reference_mscnt": mscnt[ref],
            "map_mscnt": mscnt[~ref], "map_usteps_mod_rev": expected[~ref],
            "summary": summary or "missing"}


def preprocess(run: e2_run.Run, quality: dict) -> dict:
    plan = run.plan
    raw = (quality["means_mm"] - plan.target_mm) * 1000.0
    ref = np.isin(plan.state, ("DATUM", "TURN_REF", "END_REF"))
    anchors = []
    for state in ("DATUM", "TURN_REF", "END_REF"):
        mask = plan.state == state
        anchors.append((float(quality["time_s"][mask].mean()), float(raw[mask].mean())))
    anchor_t = np.asarray([item[0] for item in anchors])
    anchor_y = np.asarray([item[1] for item in anchors])
    corrections = {
        "D0": np.zeros_like(raw),
        "D1": np.interp(quality["time_s"], anchor_t[[0, 2]], anchor_y[[0, 2]]) - anchor_y[0],
        "D2": np.interp(quality["time_s"], anchor_t, anchor_y) - anchor_y[0],
    }
    up = (plan.state == "UP") & run.matched
    down = (plan.state == "DOWN") & run.matched
    common, iu, idn = np.intersect1d(plan.index[up], plan.index[down], return_indices=True)
    up_rows = np.flatnonzero(up)[iu]; down_rows = np.flatnonzero(down)[idn]
    x = plan.target_mm[up_rows]
    usteps = np.rint(x * 1600.0).astype(int)
    theta = 2 * np.pi * np.mod(usteps, 3200) / 3200.0
    order = np.argsort(x)
    out = {"x": x[order], "usteps": usteps[order], "theta": theta[order],
           "phase_bin": plan.phase_bin[up_rows][order],
           "up_time_s": quality["time_s"][up_rows][order],
           "down_time_s": quality["time_s"][down_rows][order],
           "missing": sorted(set(range(1, 401)) - set(common.tolist())),
           "anchor_time_s": anchor_t, "anchor_error_um": anchor_y,
           "reference_time_s": quality["time_s"][ref],
           "reference_error_um": raw[ref]}
    for name, correction in corrections.items():
        e_up = raw[up_rows] - correction[up_rows]
        e_down = raw[down_rows] - correction[down_rows]
        ybar = 0.5 * (e_up + e_down)
        ybar -= ybar.mean()
        out[name] = {"e_up": e_up[order], "e_down": e_down[order],
                     "ybar": ybar[order], "B": (e_up - e_down)[order],
                     "time_s": (0.5 * (quality["time_s"][up_rows] +
                                        quality["time_s"][down_rows]))[order]}
    return out


def raw_bspline(x: np.ndarray, interior: int, domain: tuple[float, float]) -> np.ndarray:
    degree = 3
    lo, hi = domain
    knots_i = np.linspace(lo, hi, interior + 2)[1:-1] if interior else np.empty(0)
    knots = np.r_[[lo] * (degree + 1), knots_i, [hi] * (degree + 1)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, degree,
                                 extrapolate=True).toarray()


def prepare_design(x: np.ndarray, theta: np.ndarray, harmonics: int,
                   slow_knots: int, coefficient_knots: int | None,
                   domain: tuple[float, float], state: dict | None = None
                   ) -> tuple[np.ndarray, dict]:
    if state is None:
        state = {"xbar": float(x.mean()), "domain": domain, "H": harmonics,
                 "Ks": slow_knots, "Kc": coefficient_knots}
        if slow_knots:
            basis = raw_bspline(x, slow_knots, domain)
            q = np.column_stack((np.ones(x.size), x - state["xbar"]))
            projection = np.linalg.pinv(q) @ basis
            projected = basis - q @ projection
            _, singular, vt = np.linalg.svd(projected, full_matrices=False)
            rank = int(np.count_nonzero(singular > singular[0] * 1e-10))
            state["slow_projection"] = projection
            state["slow_transform"] = vt[:rank].T
    columns = [np.ones(x.size), x - state["xbar"]]
    names = ["offset", "gain"]
    groups = {}
    if state["Ks"]:
        basis = raw_bspline(x, state["Ks"], state["domain"])
        q = np.column_stack((np.ones(x.size), x - state["xbar"]))
        slow = ((basis - q @ state["slow_projection"]) @ state["slow_transform"])
        start = sum(column.shape[1] if column.ndim == 2 else 1 for column in columns)
        columns.append(slow); groups["slow"] = slice(start, start + slow.shape[1])
        names.extend([f"slow_{j}" for j in range(slow.shape[1])])
    if state["Kc"] is None:
        for k in range(1, state["H"] + 1):
            start = sum(column.shape[1] if column.ndim == 2 else 1 for column in columns)
            columns.extend((np.cos(k * theta), np.sin(k * theta)))
            groups[f"C{k}"] = (slice(start, start + 1), slice(start + 1, start + 2))
            names.extend((f"a{k}", f"b{k}"))
    else:
        coefficient_basis = raw_bspline(x, state["Kc"], state["domain"])
        state["coefficient_count"] = coefficient_basis.shape[1]
        for k in (1, 2):
            start = sum(column.shape[1] if column.ndim == 2 else 1 for column in columns)
            columns.append(coefficient_basis * np.cos(k * theta)[:, None])
            middle = start + coefficient_basis.shape[1]
            columns.append(coefficient_basis * np.sin(k * theta)[:, None])
            end = middle + coefficient_basis.shape[1]
            groups[f"C{k}"] = (slice(start, middle), slice(middle, end))
            names.extend([f"a{k}_{j}" for j in range(coefficient_basis.shape[1])])
            names.extend([f"b{k}_{j}" for j in range(coefficient_basis.shape[1])])
        for k in range(3, state["H"] + 1):
            start = sum(column.shape[1] if column.ndim == 2 else 1 for column in columns)
            columns.extend((np.cos(k * theta), np.sin(k * theta)))
            groups[f"C{k}"] = (slice(start, start + 1), slice(start + 1, start + 2))
            names.extend((f"a{k}", f"b{k}"))
    state["groups"] = groups; state["names"] = names
    return np.column_stack(columns), state


def coefficient_penalty(state: dict, width: int) -> np.ndarray:
    """Second-difference penalty on each varying-coefficient spline block.

    Only the C1 and C2 coefficient curves are penalised. The offset, gain,
    slow trend and the stationary higher harmonics are left alone, so the
    penalty controls how fast the screw harmonics are allowed to change along
    the travel and nothing else.
    """
    blocks = []
    for k in (1, 2):
        if f"C{k}" not in state["groups"]:
            continue
        for part in state["groups"][f"C{k}"]:
            count = part.stop - part.start
            if count <= PEN_ORDER:
                continue
            differences = np.eye(count)
            for _ in range(PEN_ORDER):
                differences = np.diff(differences, axis=0)
            block = np.zeros((differences.shape[0], width))
            block[:, part] = differences
            blocks.append(block)
    return np.vstack(blocks) if blocks else np.zeros((0, width))


def fit_linear(x: np.ndarray, theta: np.ndarray, y: np.ndarray, harmonics: int,
               slow_knots: int, coefficient_knots: int | None,
               domain: tuple[float, float], robust: bool = False,
               penalty: float | None = None) -> dict:
    design, state = prepare_design(x, theta, harmonics, slow_knots,
                                   coefficient_knots, domain)
    if penalty is None:
        beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        normal_inverse = np.linalg.pinv(design.T @ design)
        edf = float(design.shape[1])
    else:
        roughness = coefficient_penalty(state, design.shape[1])
        augmented = np.vstack([design, math.sqrt(penalty) * roughness])
        beta, _, _, _ = np.linalg.lstsq(
            augmented, np.r_[y, np.zeros(roughness.shape[0])], rcond=None)
        normal_inverse = np.linalg.inv(design.T @ design
                                       + penalty * (roughness.T @ roughness))
        edf = float(np.trace(design @ normal_inverse @ design.T))
    residual = y - design @ beta
    dof = max(y.size - edf, 1.0)
    variance = float(residual @ residual / dof)
    covariance = variance * normal_inverse
    result = {"beta": beta, "cov": covariance, "residual": residual,
              "sigma": math.sqrt(variance), "rms": float(np.sqrt(np.mean(residual**2))),
              "bic": float(y.size * np.log(np.mean(residual**2)) +
                           edf * np.log(y.size)),
              "p": design.shape[1], "edf": edf, "lambda": penalty,
              "state": state, "design": design}
    # Conditioning of the system actually solved: with a penalty that is the
    # augmented matrix, which is the whole point of adding one.
    solved = design if penalty is None else augmented
    scaled = solved / np.maximum(np.linalg.norm(solved, axis=0), 1e-12)
    result["condition_scaled"] = float(np.linalg.cond(scaled))
    if robust:
        rows = (design if penalty is None
                else np.vstack([design, math.sqrt(penalty) * roughness]))
        target = (y if penalty is None
                  else np.r_[y, np.zeros(roughness.shape[0])])
        huber = least_squares(lambda values: rows @ values - target, beta,
                              loss="huber", f_scale=2 * result["sigma"])
        movement = np.abs(huber.x - beta) / np.maximum(np.sqrt(np.diag(covariance)), 1e-12)
        result["huber_beta"] = huber.x
        result["huber_over_sigma_max"] = float(movement.max())
        result["huber_over_sigma_count"] = int(np.count_nonzero(movement > 1))
    return result


def predict(model: dict, x: np.ndarray, theta: np.ndarray) -> np.ndarray:
    design, _ = prepare_design(x, theta, model["state"]["H"],
                               model["state"]["Ks"], model["state"]["Kc"],
                               model["state"]["domain"], state=model["state"])
    return design @ model["beta"]


def components(model: dict, x: np.ndarray) -> dict:
    zero = np.zeros_like(x)
    design, _ = prepare_design(x, zero, model["state"]["H"],
                               model["state"]["Ks"], model["state"]["Kc"],
                               model["state"]["domain"], state=model["state"])
    gain = model["beta"][1] * (x - model["state"]["xbar"])
    slow = np.zeros_like(x)
    slow_sigma = np.zeros_like(x)
    if "slow" in model["state"]["groups"]:
        group = model["state"]["groups"]["slow"]
        slow_design = design[:, group]
        slow = slow_design @ model["beta"][group]
        cov = model["cov"][group, group]
        slow_sigma = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", slow_design, cov,
                                                  slow_design), 0))
    return {"gain": gain, "slow": slow, "slow_sigma": slow_sigma}


def harmonic_profile(model: dict, grid: np.ndarray, k: int) -> dict:
    group = model["state"]["groups"][f"C{k}"]
    alpha, beta_group = group
    if model["state"]["Kc"] is not None and k <= 2:
        basis = raw_bspline(grid, model["state"]["Kc"], model["state"]["domain"])
    else:
        basis = np.ones((grid.size, 1))
    a = basis @ model["beta"][alpha]
    b = basis @ model["beta"][beta_group]
    amplitude = np.hypot(a, b)
    phase = np.unwrap(np.arctan2(-b, a))
    amp_sigma, phase_sigma = [], []
    for row, av, bv, magnitude in zip(basis, a, b, amplitude):
        ga = np.zeros(model["p"]); gp = np.zeros(model["p"])
        if magnitude > 1e-12:
            ga[alpha] = row * av / magnitude
            ga[beta_group] = row * bv / magnitude
            gp[alpha] = row * bv / magnitude**2
            gp[beta_group] = -row * av / magnitude**2
        amp_sigma.append(math.sqrt(max(float(ga @ model["cov"] @ ga), 0)))
        phase_sigma.append(math.sqrt(max(float(gp @ model["cov"] @ gp), 0)))
    return {"complex": a - 1j*b, "amplitude": amplitude,
            "amplitude_sigma": np.asarray(amp_sigma), "phase": phase,
            "phase_sigma": np.asarray(phase_sigma)}


def cv_assignments(x: np.ndarray, domain: tuple[float, float]
                   ) -> list[np.ndarray]:
    """Blocked 10-fold assignments, repeated with the block edges shifted.

    A single set of contiguous blocks gives only eight interior scores, and
    their spread sets the standard error the one-SE rule compares against.
    With eight scores that error is large enough (~0.044 um) to admit models
    well away from the CV optimum - it was what let K_C = 1 through, a model
    that is worst of its family on every other measure. Sliding the edges and
    pooling the folds scores the same data under five different partitions, so
    a target that sat at a block boundary in one partition sits inside a block
    in another. The estimate is the same quantity, measured more precisely:
    the standard error drops to ~0.013 um and the rule can discriminate.
    """
    width = (domain[1] - domain[0]) / 10.0
    out = []
    for offset in CV_OFFSETS:
        edges = np.linspace(domain[0], domain[1], 11) + offset * width
        out.append(np.clip(np.digitize(x, edges[1:-1]), 0, 9))
    return out


def blocked_cv(data: dict) -> tuple[list[dict], dict, dict]:
    x, theta, y = data["x"], data["theta"], data["D2"]["ybar"]
    domain = (float(x.min()), float(x.max()))
    folds = cv_assignments(x, domain)
    rows = []
    # M_B: stationary harmonics. M_C_fixed: varying coefficients on a fixed
    # knot grid, kept as the diagnostic that shows what K_C costs and buys.
    # M_C: the selectable varying-coefficient family - one deliberately rich
    # basis, with lambda rather than the knot count doing the regularising.
    candidates = [("M_B", h, ks, None, None) for h in (2, 3, 4) for ks in range(9)]
    candidates += [("M_C_fixed", h, ks, kc, None) for h in (2, 3, 4)
                   for ks in range(9) for kc in range(1, 6)]
    candidates += [("M_C", h, ks, PEN_KNOTS, lam) for h in (2, 3, 4)
                   for ks in range(9) for lam in PEN_LAMBDAS]
    for kind, h, ks, kc, lam in candidates:
        full = fit_linear(x, theta, y, h, ks, kc, domain, penalty=lam)
        scored, interior_rmse, low, high = [], [], [], []
        for fold in folds:
            for held in range(10):
                train = fold != held; test = ~train
                if test.sum() < CV_MIN_TEST:
                    continue
                trial = fit_linear(x[train], theta[train], y[train], h, ks, kc,
                                   domain, penalty=lam)
                value = float(np.sqrt(np.mean(
                    (y[test] - predict(trial, x[test], theta[test]))**2)))
                scored.append(value)
                if held == 0:
                    low.append(value)
                elif held == 9:
                    high.append(value)
                else:
                    interior_rmse.append(value)
        interior = np.asarray(interior_rmse)
        rows.append({"model": kind, "H": h, "K_s": ks,
                     "K_C": "stationary" if kc is None else kc,
                     "lambda": "" if lam is None else lam,
                     "parameters": full["p"],
                     "effective_parameters": round(full["edf"], 2),
                     "residual_rms_um": full["rms"],
                     "BIC": full["bic"], "CV_RMSE_um": float(interior.mean()),
                     "CV_SE_um": float(interior.std(ddof=1) / np.sqrt(interior.size)),
                     "CV_folds_scored": int(interior.size),
                     "end_low_RMSE_um": float(np.mean(low)),
                     "end_high_RMSE_um": float(np.mean(high)),
                     "fold_RMSE_json": json.dumps(scored)})

    def one_se(kind: str) -> dict:
        """Cheapest stable model within one standard error of the best CV.

        Cost is the effective parameter count, so a penalised fit is compared
        with a fixed-knot one on the flexibility it actually uses rather than
        on how many columns its basis happens to have.

        The stability guard is not decoration. Interior folds never ask a
        model to predict beyond its training support, so they cannot see a
        basis whose edge coefficients are only weakly determined - and the
        penalty here covers the C1/C2 coefficient curves, not the slow trend,
        so a large K_s buys flexibility that the interior score does not
        charge for. The two held-out end blocks do ask exactly that question.
        Models that fail it miss by orders of magnitude, not by a little:
        among the candidates here the passes sit near 3x the interior score
        and the failures run from 60x to 20000x, so the cut is not delicate.
        """
        subset = [row for row in rows if row["model"] == kind]
        best = min(subset, key=lambda row: row["CV_RMSE_um"])
        eligible = [row for row in subset
                    if row["CV_RMSE_um"] <= best["CV_RMSE_um"] + best["CV_SE_um"]]
        stable = [row for row in eligible
                  if max(row["end_low_RMSE_um"], row["end_high_RMSE_um"])
                  <= CV_END_FACTOR * row["CV_RMSE_um"]]
        return min(stable or eligible,
                   key=lambda row: (row["effective_parameters"],
                                    row["CV_RMSE_um"]))

    return rows, one_se("M_B"), one_se("M_C")


def fit_selected(data: dict, row: dict, variant: str = "D2", robust: bool = True) -> dict:
    kc = None if row["K_C"] == "stationary" else int(row["K_C"])
    lam = float(row["lambda"]) if row.get("lambda", "") != "" else None
    return fit_linear(data["x"], data["theta"], data[variant]["ybar"],
                      int(row["H"]), int(row["K_s"]), kc,
                      (float(data["x"].min()), float(data["x"].max())),
                      robust=robust, penalty=lam)


def slow_coefficient_correlation(model: dict) -> float:
    """Largest absolute sample correlation between s(x) and C1/C2 columns."""
    groups = model["state"]["groups"]
    if "slow" not in groups or model["state"]["Kc"] is None:
        return 0.0
    slow = model["design"][:, groups["slow"]]
    coefficient_columns = []
    for key in ("C1", "C2"):
        coefficient_columns.extend(range(groups[key][0].start, groups[key][0].stop))
        coefficient_columns.extend(range(groups[key][1].start, groups[key][1].stop))
    varying = model["design"][:, coefficient_columns]
    a = slow - slow.mean(axis=0); b = varying - varying.mean(axis=0)
    a /= np.maximum(np.linalg.norm(a, axis=0), 1e-12)
    b /= np.maximum(np.linalg.norm(b, axis=0), 1e-12)
    return float(np.max(np.abs(a.T @ b)))


def sliding_profile(data: dict, field: str = "ybar") -> list[dict]:
    x, theta = data["x"], data["theta"]
    values = data["D2"][field]
    centres = np.arange(math.ceil(x.min() + 6), math.floor(x.max() - 6) + .1, 2.0)
    out = []
    for centre in centres:
        mask = np.abs(x - centre) <= 6.0
        bins = np.floor(np.mod(theta[mask], 2*np.pi) / (2*np.pi) * 50).astype(int)
        if mask.sum() < 25 or np.unique(bins).size < 25:
            continue
        local_x = x[mask] - centre
        design = np.column_stack((np.ones(mask.sum()), local_x,
                                  np.cos(theta[mask]), np.sin(theta[mask]),
                                  np.cos(2*theta[mask]), np.sin(2*theta[mask])))
        beta, _, _, _ = np.linalg.lstsq(design, values[mask], rcond=None)
        residual = values[mask] - design @ beta
        variance = float(residual @ residual / (mask.sum() - design.shape[1]))
        cov = variance * np.linalg.pinv(design.T @ design)
        record = {"x": float(centre), "n": int(mask.sum()),
                  "phase_bins": int(np.unique(bins).size), "rms": float(np.sqrt(np.mean(residual**2)))}
        for k, start in ((1, 2), (2, 4)):
            a, b = beta[start:start+2]
            cab = cov[start:start+2, start:start+2]
            amp = float(np.hypot(a, b))
            ga = np.asarray([a, b]) / max(amp, 1e-12)
            gp = np.asarray([b, -a]) / max(amp**2, 1e-12)
            record[f"C{k}"] = complex(float(a), float(-b))
            record[f"A{k}"] = amp
            record[f"A{k}_sigma"] = math.sqrt(max(float(ga @ cab @ ga), 0))
            record[f"phi{k}"] = float(np.arctan2(-b, a))
            record[f"phi{k}_sigma"] = math.sqrt(max(float(gp @ cab @ gp), 0))
        out.append(record)
    return out


def running_profile(x: np.ndarray, values: np.ndarray, width: float,
                    step: float = 1.0) -> dict:
    grid = np.arange(math.ceil(x.min()), math.floor(x.max()) + .1, step)
    mean, sigma = [], []
    for point in grid:
        here = np.abs(x - point) <= width / 2
        sample = values[here]
        mean.append(float(sample.mean()))
        sigma.append(float(sample.std(ddof=1) / np.sqrt(sample.size)))
    return {"x": grid, "mean": np.asarray(mean), "sigma": np.asarray(sigma)}


def reversal_fit(data: dict, variant: str) -> dict:
    x, theta, values = data["x"], data["theta"], data[variant]["B"]
    xc = x - x.mean()
    design = np.column_stack((np.ones(x.size), xc,
                              np.cos(theta), np.sin(theta),
                              np.cos(2*theta), np.sin(2*theta)))
    beta, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ beta
    variance = float(residual @ residual / (values.size - design.shape[1]))
    cov = variance * np.linalg.pinv(design.T @ design)
    return {"beta": beta, "cov": cov, "residual": residual,
            "mean_um": float(beta[0]), "mean_sigma_um": float(np.sqrt(cov[0, 0])),
            "slope_um_per_100mm": float(beta[1] * 100),
            "slope_sigma_um_per_100mm": float(np.sqrt(cov[1, 1]) * 100)}


def external_windows() -> tuple[dict, dict]:
    run1 = e1_run.load(); q1 = e4f.quality(run1)
    error1, _ = e4f.drift_corrected_error(run1, q1["means_mm"])
    windows1 = e4f.make_windows(run1, error1, list(e1_run.WINDOWS), "E1")
    run4 = e4_run.load(); q4 = e4f.quality(run4)
    error4, _ = e4f.drift_corrected_error(run4, q4["means_mm"])
    windows4 = e4f.make_windows(run4, error4, list(e4_run.WINDOWS), "E4")
    return windows1, windows4


def validate_windows(model_b: dict, model_c: dict, windows1: dict,
                     windows4: dict) -> list[dict]:
    rows = []
    for source, windows in (("E1", windows1), ("E4", windows4)):
        for name, window in windows.items():
            x, theta, observed = window["x"], window["theta"], window["ybar"]
            raw = observed - observed.mean()
            rb = observed - predict(model_b, x, theta); rb -= rb.mean()
            rc = observed - predict(model_c, x, theta); rc -= rc.mean()
            local = window["fit"]["residual"]
            metrics = {"raw": float(np.sqrt(np.mean(raw**2))),
                       "M_B": float(np.sqrt(np.mean(rb**2))),
                       "M_C": float(np.sqrt(np.mean(rc**2))),
                       "local": float(np.sqrt(np.mean(local**2)))}
            gap = metrics["raw"] - metrics["local"]
            record = {"source": source, "window": name,
                      "x_centre_mm": window["centre"],
                      **{f"rms_{key}_um": value for key, value in metrics.items()},
                      "gap_closed_MB": (metrics["raw"] - metrics["M_B"]) / gap if gap else np.nan,
                      "gap_closed_MC": (metrics["raw"] - metrics["M_C"]) / gap if gap else np.nan}
            for k in (1, 2):
                actual = window["spectrum"][k-1]["complex"]
                mb = harmonic_profile(model_b, np.asarray([window["centre"]]), k)["complex"][0]
                mc = harmonic_profile(model_c, np.asarray([window["centre"]]), k)["complex"][0]
                record[f"C{k}_difference_MB_um"] = float(abs(actual - mb))
                record[f"C{k}_difference_MC_um"] = float(abs(actual - mc))
            rows.append(record)
    return sorted(rows, key=lambda row: (row["x_centre_mm"], row["source"]))


def iso_metrics(up: np.ndarray, down: np.ndarray, s_up: np.ndarray,
                s_down: np.ndarray) -> dict:
    mid = 0.5 * (up + down)
    high = np.maximum(up + 2*s_up, down + 2*s_down)
    low = np.minimum(up - 2*s_up, down - 2*s_down)
    return {"E_up": float(np.ptp(up)), "E_down": float(np.ptp(down)),
            "E": float(max(up.max(), down.max()) - min(up.min(), down.min())),
            "M": float(np.ptp(mid)), "A": float(high.max() - low.min()),
            "R_up": float(4*s_up.max()), "R_down": float(4*s_down.max()),
            "R": float(max(np.max(2*s_up + 2*s_down + np.abs(up-down)),
                           4*s_up.max(), 4*s_down.max()))}


def validate_iso(model: dict, reversal: dict) -> tuple[list[dict], dict]:
    source = json.loads(ISO_RESULTS.read_text(encoding="utf-8"))["per_target"]
    rows0 = list(source.values())
    x = np.asarray([row["cmd_mm"] for row in rows0])
    usteps = np.rint(x * 1600).astype(int)
    theta = 2*np.pi*np.mod(usteps, 3200)/3200
    geometry = predict(model, x, theta)
    up = np.asarray([row["xbar_up"] for row in rows0])
    down = np.asarray([row["xbar_dn"] for row in rows0])
    s_up = np.asarray([row["s_up"] for row in rows0])
    s_down = np.asarray([row["s_dn"] for row in rows0])
    offset = float(np.mean(0.5*(up + down) - geometry))
    g_up = up - geometry - offset; g_down = down - geometry - offset
    beta0 = reversal["mean_um"]
    gb_up = g_up - beta0/2; gb_down = g_down + beta0/2
    rows = []
    for i in range(x.size):
        rows.append({"target_mm": x[i], "electrical_phase": int(usteps[i] % 64),
                     "raw_mean_um": float(0.5*(up[i]+down[i])),
                     "G_mean_um": float(0.5*(g_up[i]+g_down[i])),
                     "G_B_mean_um": float(0.5*(gb_up[i]+gb_down[i])),
                     "raw_up_um": up[i], "raw_down_um": down[i],
                     "G_up_um": g_up[i], "G_down_um": g_down[i],
                     "G_B_up_um": gb_up[i], "G_B_down_um": gb_down[i]})
    metrics = {"raw": iso_metrics(up, down, s_up, s_down),
               "G": iso_metrics(g_up, g_down, s_up, s_down),
               "G_plus_reversal": iso_metrics(gb_up, gb_down, s_up, s_down),
               "note": "R and repeatability standard deviations are unchanged by deterministic compensation"}
    # Per the E5 hand-over convention, compensation changes the predicted
    # systematic terms E/M/A; repeatability R and the measured scatters s do
    # not change.  Preserve the measured R values rather than reinterpreting a
    # session-to-session reversal mismatch as repeatability.
    for level in ("G", "G_plus_reversal"):
        for key in ("R_up", "R_down", "R"):
            metrics[level][key] = metrics["raw"][key]
    return rows, metrics


def residual_periodogram(x: np.ndarray, values: np.ndarray,
                         low: float, high: float, count: int = 3000) -> tuple[np.ndarray, np.ndarray]:
    frequency = np.linspace(low, high, count)
    centred = values - values.mean()
    power = lombscargle(x, centred, 2*np.pi*frequency, normalize=True)
    return frequency, power


def phase_reference(profile: dict) -> float:
    weight = 1 / np.maximum(profile["phase_sigma"], np.radians(.2))**2
    return float(np.average(profile["phase"], weights=weight))


def nearest_unwrap(phase: float, reference: float) -> float:
    return float(phase + 2*np.pi*round((reference-phase)/(2*np.pi)))


# ---- presentation figures: greyscale-safe, identity by line style and marker
BAND = "#e2e5e8"          # light fill for the E2 / D2 ±2σ bands
HATCH_EDGE = "#aab1b6"    # hatched fill for the model ±2σ bands
MID = "#8b949a"


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(color=GRID, lw=.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=.6)


def endpoint_shading(ax: plt.Axes, label: bool = False) -> None:
    """Neutral bands at the physical ends; text only where asked."""
    for x in (MOTOR_END_MM, SUPPORT_BEARING_MM):
        ax.axvspan(x - 1.4, x + 1.4, color="#000000", alpha=.06, lw=0, zorder=0)
    if label:
        transform = ax.get_xaxis_transform()
        ax.text(MOTOR_END_MM + 1.8, 1.02, "motor end", transform=transform,
                ha="left", va="bottom", fontsize=6.6, color=INK2, clip_on=False)
        ax.text(SUPPORT_BEARING_MM - 1.8, 1.02, "support bearing", transform=transform,
                ha="right", va="bottom", fontsize=6.6, color=INK2, clip_on=False)


def letter(ax: plt.Axes, text: str) -> None:
    ax.text(.012, .975, f"({text})", transform=ax.transAxes, ha="left", va="top",
            fontsize=8.2, weight="bold", color=INK, zorder=6,
            bbox=dict(boxstyle="square,pad=.15", fc="white", ec="none", alpha=.85))


def screw_theta(x: np.ndarray) -> np.ndarray:
    """Commanded screw angle, exactly as preprocess() derives it per target."""
    return 2 * np.pi * np.mod(np.rint(x * 1600.0).astype(int), 3200) / 3200.0


ZOOMS_MM = ((6.0, 14.0), (70.0, 78.0))   # four screw turns each


ZOOM_SPAN_UM = 5.0   # both zooms share this y-span so their shapes compare fairly


def _f1_title(ax, label, text, raised=False):
    ax.set_title(rf"$\bf{{{label}}}$ {text}", loc="left", y=1.06 if raised else 1.0, pad=4)


def _f1_legend(ax, handles, raised=False):
    ax.legend(handles=handles, frameon=False, loc="lower right",
              bbox_to_anchor=(1.0, 1.055 if raised else 1.0), ncol=len(handles),
              fontsize=6.8, handlelength=1.8, columnspacing=1.1, borderaxespad=.2)


def figure_f1(data: dict, model: dict) -> None:
    """Figure 1a: what the error looks like, and whether the model captures it."""
    x, y = data["x"], data["D2"]["ybar"]
    grid = np.linspace(x.min(), x.max(), 800)
    comp = components(model, grid)
    slow_part = comp["gain"] + comp["slow"]
    envelope = (harmonic_profile(model, grid, 1)["amplitude"]
                + harmonic_profile(model, grid, 2)["amplitude"])

    fig = plt.figure(figsize=(7.45, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.4, .85], left=.10, right=.975,
                          top=.92, bottom=.15, hspace=.45, wspace=.16)
    ax_a = fig.add_subplot(gs[0, :])
    ax_z = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

    for lo, hi in ZOOMS_MM:
        ax_a.axvspan(lo, hi, color="#000000", alpha=.07, lw=0, zorder=0)
        ax_a.text((lo + hi) / 2, 1.0, "a′" if lo < 40 else "a″",
                  transform=ax_a.get_xaxis_transform(), ha="center", va="top",
                  fontsize=7, color=INK2)
    ax_a.fill_between(grid, slow_part - envelope, slow_part + envelope, facecolor="#c9ced2",
                      edgecolor="#aab1b6", lw=.5, alpha=.85, zorder=1)
    ax_a.plot(x, y, "o", ms=1.7, mfc="#5d656b", mec="none", zorder=2)
    ax_a.plot(grid, slow_part, color=INK, lw=1.7, zorder=3)
    ax_a.set_ylabel(r"direction-averaged deviation $\bar y$ [µm]")
    ax_a.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0))
    ax_a.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(.5))
    ax_a.grid(which="minor", axis="y", color=GRID, lw=.3)
    ax_a.tick_params(which="minor", length=0)
    ax_a.set_xlabel("commanded position x [mm]")
    _f1_title(ax_a, "(a)", "Full travel", raised=True)
    _f1_legend(ax_a, [matplotlib.patches.Patch(facecolor="#c9ced2", edgecolor="#aab1b6", alpha=.85,
                                               label="Cumulative error (slow part) ± 1st & 2nd "
                                                     "Harmonics Envelope"),
                      Line2D([], [], marker="o", ls="", ms=2.6, mfc="#5d656b", mec="none",
                             label=f"{x.size} targets"),
                      Line2D([], [], color=INK, lw=1.7, label="slow part (gain + trend)")],
               raised=True)
    endpoint_shading(ax_a, label=True)
    clean_axes(ax_a)
    ax_a.set_xlim(MOTOR_END_MM - 2, SUPPORT_BEARING_MM + 2)

    for ax, (lo, hi), label in zip(ax_z, ZOOMS_MM, ("(a′)", "(a″)")):
        dense = np.linspace(lo, hi, 2000)
        curve = predict(model, dense, screw_theta(dense))
        here = (x >= lo) & (x <= hi)
        ax.plot(dense, curve, color=MID, lw=.9, zorder=1)
        ax.plot(x[here], y[here], "o", ms=3.0, mfc=INK, mec="white", mew=.4, zorder=2)
        values = np.concatenate([curve, y[here]])
        centre = .5 * (values.min() + values.max())
        ax.set_ylim(centre - ZOOM_SPAN_UM / 2, centre + ZOOM_SPAN_UM / 2)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1.0))
        ax.set_xlim(lo, hi)
        ax.set_xticks(np.arange(lo, hi + .1, 2))
        ax.tick_params(labelsize=7)
        ax.set_xlabel("commanded position x [mm]", fontsize=8)
        _f1_title(ax, label, f"Zoom {lo:.0f}–{hi:.0f} mm")
        clean_axes(ax)
    ax_z[0].set_ylabel(r"$\bar y$ [µm]")
    _f1_legend(ax_z[1], [Line2D([], [], color=MID, lw=.9, label="full model"),
                         Line2D([], [], marker="o", ls="", ms=3, mfc=INK, mec="white",
                                label="target")])
    fig.text(.10, .03,
             f"Figure 1a. E2 deviation over the full travel ({x.size} targets, D2 drift correction): the bold line "
             "is the slow part and the band the\nonce- and twice-per-turn wobble around it. Zooms span four screw "
             f"turns each, on a common {ZOOM_SPAN_UM:g} µm scale.",
             fontsize=7.0, color=INK2, linespacing=1.3)
    save(fig, "e2_figure1_full_travel_map", PRESENTATION)


def figure_f1b(data: dict, model: dict, v2: list[dict]) -> None:
    """Figure 1b: the slow error, and what the model leaves unexplained."""
    x, theta, y = data["x"], data["theta"], data["D2"]["ybar"]
    residual = y - predict(model, x, theta)
    grid = np.linspace(x.min(), x.max(), 800)
    comp = components(model, grid)
    slow_part = comp["gain"] + comp["slow"]
    running = running_profile(x, residual, 2.0, .5)   # ±1 mm: one 2 mm screw turn
    gain_ppm = float(model["beta"][1] * 1000)
    rms = float(np.sqrt(np.mean(residual**2)))
    floor = float(np.median([row["rms_local_um"] for row in v2 if row["source"] == "E4"]))

    fig = plt.figure(figsize=(7.45, 5.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, .8], left=.10, right=.975,
                          top=.90, bottom=.16, hspace=.32)
    ax_b = fig.add_subplot(gs[0])
    ax_c = fig.add_subplot(gs[1], sharex=ax_b)

    ax_b.plot(grid, comp["gain"], color=INK2, lw=1.0, ls=(0, (5, 2.5)), zorder=2)
    ax_b.plot(grid, slow_part, color=INK, lw=1.7, zorder=3)
    ax_b.axhline(0, color=MID, lw=.55, zorder=0)
    ax_b.set_ylabel("slow error [µm]")
    ax_b.tick_params(labelbottom=False)
    _f1_title(ax_b, "(b)", "Cumulative slow error", raised=True)
    _f1_legend(ax_b, [Line2D([], [], color=INK2, lw=1.0, ls=(0, (5, 2.5)),
                             label=f"gain {gain_ppm:+.1f} ppm".replace("-", "−")),
                      Line2D([], [], color=INK, lw=1.7,
                             label=f"total slow error, {np.ptp(slow_part):.1f} µm peak-to-peak")],
               raised=True)

    ax_c.axhspan(-floor, floor, color=BAND, lw=0, zorder=0)
    ax_c.plot(x, residual, "o", ms=1.6, mfc="#a9afb3", mec="none", zorder=1)
    ax_c.plot(running["x"], running["mean"], color=INK, lw=1.15, zorder=2)
    ax_c.axhline(0, color=MID, lw=.55, zorder=1)
    ax_c.set_ylabel("residual [µm]")
    ax_c.set_xlabel("commanded position x [mm]")
    _f1_title(ax_c, "(c)", "Residual after the model")
    _f1_legend(ax_c, [matplotlib.patches.Patch(facecolor=BAND, label=f"E4 local floor ±{floor:.2f} µm"),
                      Line2D([], [], marker="o", ls="", ms=2.6, mfc="#a9afb3", mec="none",
                             label="target residual"),
                      Line2D([], [], color=INK, lw=1.15, label="2 mm running mean"),
                      Line2D([], [], ls="", label=f"rms = {rms:.2f} µm")])
    for i, ax in enumerate((ax_b, ax_c)):
        endpoint_shading(ax, label=i == 0)
        clean_axes(ax)
        ax.set_xlim(MOTOR_END_MM - 2, SUPPORT_BEARING_MM + 2)
    fig.text(.10, .03,
             "Figure 1b. (b) The slow error with the screw wobble averaged out, against its straight-line gain. "
             "(c) What the model leaves,\nagainst the E4 local-fit floor; the 2 mm running mean spans exactly one "
             "screw turn, so it holds no leftover wobble.",
             fontsize=7.0, color=INK2, linespacing=1.3)
    save(fig, "e2_figure1b_slow_error_and_residual", PRESENTATION)


E1_NUDGE_MM = 1.4   # E1 markers sit just right of coincident E4 windows


def window_marker(ax, x, y, yerr, source):
    ax.errorbar(x + (E1_NUDGE_MM if source == "E1" else 0.0), y, yerr=yerr,
                fmt="o" if source == "E4" else "s", ms=4.4 if source == "E4" else 4.1,
                mfc=INK if source == "E4" else "white", mec=INK, mew=.85,
                ecolor=INK, elinewidth=.7, capsize=1.6, capthick=.7, zorder=5,
                ls="none")


def figure_f2(sliding: list[dict], model: dict, windows1: dict, windows4: dict) -> None:
    sx = np.asarray([row["x"] for row in sliding])
    grid = np.linspace(model["state"]["domain"][0], model["state"]["domain"][1], 600)
    # The model profile is drawn only where the sliding fit also has centres, so
    # the two curves are compared over one span. The full grid is still carried
    # for the phase reference and for unwrapping windows that sit outside it.
    drawn = np.linspace(sx.min(), sx.max(), 600)
    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.1), sharex=True)
    fig.subplots_adjust(left=.085, right=.98, top=.865, bottom=.14,
                        hspace=.20, wspace=.24)
    for column, k in enumerate((1, 2)):
        mc = harmonic_profile(model, grid, k)
        mcd = harmonic_profile(model, drawn, k)
        slide_amp = np.asarray([row[f"A{k}"] for row in sliding])
        slide_amp_sigma = np.asarray([row[f"A{k}_sigma"] for row in sliding])
        slide_phase = np.unwrap(np.asarray([row[f"phi{k}"] for row in sliding]))
        slide_phase_sigma = np.asarray([row[f"phi{k}_sigma"] for row in sliding])
        reference = phase_reference(mc)
        model_shift = -np.degrees(mcd["phase"] - reference) / k
        slide_phase = np.asarray([nearest_unwrap(p, np.interp(x, grid, mc["phase"]))
                                  for p, x in zip(slide_phase, sx)])
        slide_shift = -np.degrees(slide_phase - reference) / k
        rows = ((slide_amp, 2*slide_amp_sigma, mcd["amplitude"], 2*mcd["amplitude_sigma"]),
                (slide_shift, 2*np.degrees(slide_phase_sigma)/k,
                 model_shift, 2*np.degrees(mcd["phase_sigma"])/k))
        for row, (s_val, s_err, m_val, m_err) in enumerate(rows):
            ax = axes[row, column]
            ax.fill_between(sx, s_val - s_err, s_val + s_err, facecolor=BAND,
                            lw=0, zorder=1)
            ax.fill_between(drawn, m_val - m_err, m_val + m_err, facecolor="none",
                            edgecolor=HATCH_EDGE, hatch="//////", lw=0, zorder=2)
            ax.plot(drawn, m_val, color=INK2, lw=1.0, ls=(0, (5, 2.5)), zorder=3)
            ax.plot(sx, s_val, color=INK, lw=1.35, zorder=4)
        axes[1, column].axhline(0, color=MID, lw=.55, zorder=0)
        for source, windows in (("E4", windows4), ("E1", windows1)):
            for window in windows.values():
                spec = window["spectrum"][k-1]
                xpos = window["centre"]
                phase = nearest_unwrap(spec["phase_rad"], np.interp(xpos, grid, mc["phase"]))
                window_marker(axes[0, column], xpos, spec["amplitude_um"],
                              2*spec["amplitude_sigma_um"], source)
                window_marker(axes[1, column], xpos, -np.degrees(phase - reference)/k,
                              2*spec["phase_sigma_deg"]/k, source)
        axes[0, column].set_title(
            f"k = {k}  ·  {'once' if k == 1 else 'twice'} per revolution",
            loc="left", pad=15, fontsize=9.2)
        axes[1, column].set_xlabel("position along the travel [mm]")
    axes[0, 0].set_ylabel(r"amplitude $A_k$ [µm]")
    axes[1, 0].set_ylabel(r"pattern shift $\Delta\theta_k$ [°]")
    for i, ax in enumerate(axes.flat):
        endpoint_shading(ax, label=i < 2)
        clean_axes(ax)
        ax.set_xlim(MOTOR_END_MM - 2, SUPPORT_BEARING_MM + 2)
        letter(ax, "abcd"[i])
    handles = [
        Line2D([], [], color=INK, lw=1.35, label="E2 sliding fit"),
        matplotlib.patches.Patch(facecolor=BAND, label="±2σ"),
        Line2D([], [], color=INK2, lw=1.0, ls=(0, (5, 2.5)), label="Fitted Global Model"),
        matplotlib.patches.Patch(facecolor="white", edgecolor=HATCH_EDGE, hatch="//////",
                                 lw=0, label="±2σ"),
        Line2D([], [], marker="o", ls="", ms=4.4, mfc=INK, mec=INK, label="E4 window"),
        Line2D([], [], marker="s", ls="", ms=4.1, mfc="white", mec=INK, label="E1 window"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.075, .985),
               ncol=6, frameon=False, fontsize=7.2, handlelength=1.9,
               columnspacing=1.3)
    fig.text(.085, .035,
             "Figure 2. Continuous E2 harmonic profile (sliding fit and varying-coefficient model $M_C$) "
             "against the independent E1/E4\nwindow estimates. All bands and error bars are ±2σ; E1 markers "
             "are offset 1.4 mm to the right where they coincide with E4.",
             fontsize=7.1, color=INK2, linespacing=1.3)
    save(fig, "e2_figure2_continuous_harmonic_profile", PRESENTATION)


VALIDATION_LEVELS = (  # (key, label, marker, filled, size)
    ("raw", "raw (no correction)", "o", False, 4.6),
    ("M_B", r"$M_B$ (constant harmonics)", "^", False, 4.8),
    ("M_C", r"$M_C$ (varying harmonics)", "o", True, 4.6),
    ("local", "local floor (window's own fit)", "_", True, 8.5),
)


def figure_f3(v2: list[dict], v3: list[dict], metrics: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.25),
                             gridspec_kw={"width_ratios": [1.05, 1]})
    fig.subplots_adjust(left=.075, right=.975, top=.76, bottom=.23, wspace=.25)

    ax = axes[0]
    xs = np.arange(len(v2))
    values = {key: np.asarray([row[f"rms_{key}_um"] for row in v2])
              for key, *_ in VALIDATION_LEVELS}
    for i in xs:
        column = [values[key][i] for key, *_ in VALIDATION_LEVELS]
        ax.plot([i, i], [min(column), max(column)], color="#c3c8cc", lw=1.0, zorder=1)
    for key, label, marker, filled, size in VALIDATION_LEVELS:
        ax.plot(xs, values[key], marker, ls="none", ms=size,
                mfc=INK if filled else "white", mec=INK, mew=1.1 if marker == "_" else .85,
                zorder=3, label=label)
    ax.set_xticks(xs, [f"{row['x_centre_mm']:.0f}\n{row['source']}" for row in v2])
    ax.tick_params(axis="x", length=0, labelsize=7.2)
    ax.set_xlim(-.6, len(v2) - .4)
    ax.set_ylim(0, math.ceil(max(values["raw"]) * 1.1 * 5) / 5)
    ax.set_xlabel("validation window: centre [mm] and experiment")
    ax.set_ylabel("residual rms [µm]")
    ax.set_title(r"$\bf{(a)}$ Independent E1/E4 window validation", loc="left", pad=6)
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(-.01, 1.10), ncol=2,
              fontsize=6.9, handlelength=1.4, columnspacing=1.2, borderaxespad=0)

    ax = axes[1]
    xr = np.asarray([row["target_mm"] for row in v3])
    raw = np.asarray([row["raw_mean_um"] for row in v3])
    after_g = np.asarray([row["G_mean_um"] for row in v3])
    for x, a, b in zip(xr, raw, after_g):
        ax.plot([x, x], [a, b], color="#c3c8cc", lw=.9, zorder=1)
    ax.plot(xr, raw, "o", ms=4.2, mfc="white", mec=MID, mew=.85, zorder=2,
            label=rf"before: $E={metrics['raw']['E']:.2f}$, $M={metrics['raw']['M']:.2f}$ µm")
    ax.plot(xr, after_g, "o", ms=4.2, mfc=INK, mec=INK, zorder=3,
            label=rf"after $G$: $E={metrics['G']['E']:.2f}$, $M={metrics['G']['M']:.2f}$ µm")
    ax.text(.12, .965,
            rf"$G+B$ leaves the mean deviation unchanged;" "\n"
            rf"predicted $E={metrics['G_plus_reversal']['E']:.2f}$ µm",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.9, color=INK2)
    ax.axhline(0, color=MID, lw=.55, zorder=0)
    endpoint_shading(ax, label=True)
    ax.set_xlim(MOTOR_END_MM - 2, SUPPORT_BEARING_MM + 2)
    ax.set_xlabel("commanded position [mm]")
    ax.set_ylabel("ISO mean deviation [µm]")
    ax.set_title(r"$\bf{(b)}$ Predicted E0 ISO improvement", loc="left", pad=15)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(.53, .02), ncol=1,
              fontsize=6.9, handlelength=1.2, borderaxespad=0)
    for ax in axes:
        clean_axes(ax)
    axes[0].grid(axis="x", visible=False)
    fig.text(.075, .045,
             "Figure 3. (a) Residual rms in windows not used for fitting; the grey line spans each window's "
             "range from raw to its own local floor.\n(b) E0 ISO mean deviation before and after the E2 geometry "
             "correction G; repeatability R is unchanged.",
             fontsize=7.1, color=INK2, linespacing=1.3)
    save(fig, "e2_figure3_validation", PRESENTATION)


def figure_f4(data: dict, reversal: dict, windows4: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.35, 4.3))
    fig.subplots_adjust(left=.10, right=.975, top=.84, bottom=.19)
    profiles = {name: running_profile(data["x"], data[name]["B"], 6.0)
                for name in ("D0", "D1", "D2")}
    d2 = profiles["D2"]
    ax.fill_between(d2["x"], d2["mean"] - 2*d2["sigma"], d2["mean"] + 2*d2["sigma"],
                    facecolor=BAND, lw=0, zorder=1)
    styles = {"D0": dict(color=MID, lw=.95, ls=(0, (5, 2.5))),
              "D1": dict(color=MID, lw=1.0, ls=(0, (1.2, 1.6)))}
    for name, kw in styles.items():
        ax.plot(profiles[name]["x"], profiles[name]["mean"], zorder=2, **kw)
    ax.plot(d2["x"], d2["mean"], color=INK, lw=1.45, zorder=3)
    ax.axhline(reversal["mean_um"], color=INK2, lw=.8, ls=(0, (6, 2, 1, 2)), zorder=2)
    ax.text(46.5, reversal["mean_um"] - .012,
            rf"$\beta_0$ = {reversal['mean_um']:.2f} µm".replace("-", "−"),
            ha="left", va="top", fontsize=6.9, color=INK2)
    for window in windows4.values():
        metric = e4f.reversal_metrics(window)
        ax.errorbar(window["centre"], metric["mean_um"], yerr=2*metric["mean_sigma_um"],
                    fmt="o", ms=4.4, mfc=INK, mec="white", mew=.6, ecolor=INK,
                    elinewidth=.75, capsize=1.6, capthick=.75, zorder=5)
    endpoint_shading(ax, label=True)
    clean_axes(ax)
    ax.set_xlim(MOTOR_END_MM - 2, SUPPORT_BEARING_MM + 2)
    ax.set_xlabel("position along the travel [mm]")
    ax.set_ylabel(r"reversal $B$ [µm]")
    ax.set_title("Directional reversal along the travel", loc="left", pad=15)
    handles = [
        Line2D([], [], color=INK, lw=1.45, label="E2, D2 drift correction (default)"),
        matplotlib.patches.Patch(facecolor=BAND, label="±2σ"),
        Line2D([], [], color=MID, lw=.95, ls=(0, (5, 2.5)), label="D0"),
        Line2D([], [], color=MID, lw=1.0, ls=(0, (1.2, 1.6)), label="D1"),
        Line2D([], [], marker="o", ls="", ms=4.4, mfc=INK, mec="white", label="E4 window"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.09, .99), ncol=5,
               frameon=False, fontsize=7.1, handlelength=1.9, columnspacing=1.3)
    fig.text(.10, .045,
             "Figure 4. Six-millimetre running mean of the E2 reversal. D0 and D1 are the alternative drift "
             "corrections, shown for sensitivity;\nD1 lies under D2 (they differ by < 0.001 µm). "
             "The band and the E4 error bars are ±2σ.",
             fontsize=7.1, color=INK2, linespacing=1.3)
    save(fig, "e2_figure4_reversal_profile", PRESENTATION)


def appendix_validity(data: dict, quality: dict, validity: list[dict]) -> None:
    fig = plt.figure(figsize=(8.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.05, 1), hspace=.36, wspace=.28,
                          left=.08, right=.97, top=.95, bottom=.10)
    ax0 = fig.add_subplot(gs[0, 0]); ax1 = fig.add_subplot(gs[0, 1]); ax2 = fig.add_subplot(gs[1, :])
    counts = np.bincount(data["phase_bin"], minlength=50)[:50]
    ax0.bar(np.arange(50), counts, color=BLUE, width=.85)
    ax0.axhline(5, color=RUST, ls="--", lw=.8, label="minimum 5")
    ax0.set(xlabel="screw-phase bin", ylabel="targets per bin")
    panel_title(ax0, "a", "Phase coverage"); ax0.legend(frameon=False)
    time = data["reference_time_s"]/60
    error = data["reference_error_um"]
    ax1.plot(time, error, "o-", color=GREEN, mfc="white", mew=.8, lw=1)
    for variant, ls in (("D1", "--"), ("D2", "-")):
        anchors = data["anchor_error_um"]
        at = data["anchor_time_s"]/60
        if variant == "D1":
            ax1.plot(at[[0, 2]], anchors[[0, 2]], ls, color=DRIFT_COLOURS[variant], lw=.9)
        else:
            ax1.plot(at, anchors, ls, color=DRIFT_COLOURS[variant], lw=.9)
    ax1.set(xlabel="capture time [min]", ylabel="reference error [µm]")
    panel_title(ax1, "b", "Reference drift and correction anchors")
    def concise(row: dict) -> list[str]:
        value = row["value"]
        if row["check"] == "firmware completion":
            value = "807 stops; complete; no abort"
        elif isinstance(value, float):
            value = f"{value:.3f}"
        status = "yes" if row["pass"] is True else "no" if row["pass"] is False else "n/a"
        return [row["check"], str(value), status]

    table_rows = [concise(row) for row in validity]
    ax2.axis("off")
    table = ax2.table(cellText=table_rows, colLabels=("check", "value", "pass"),
                      colWidths=(.45, .39, .12), loc="center", cellLoc="left")
    table.auto_set_font_size(False); table.set_fontsize(6.8); table.scale(1, 1.25)
    for key, cell in table.get_celld().items():
        cell.set_edgecolor("#d7dbde"); cell.set_linewidth(.45)
        if key[0] == 0: cell.set_text_props(weight="bold"); cell.set_facecolor("#eef0f1")
    ax2.set_title(r"$\bf{(c)}$ Validity table", loc="left", pad=4)
    for ax in (ax0, ax1): finish_axis(ax)
    fig.text(.08, .025, "Appendix A1. E2 validity checks. Environment data were not recorded.",
             fontsize=7.3, color=INK2)
    save(fig, "e2_appendix_validity", OUT)


def appendix_cv(rows: list[dict], chosen_b: dict, chosen_c: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.55))
    fig.subplots_adjust(left=.07, right=.98, top=.88, bottom=.31, wspace=.31)
    kc_values = ["stationary", 1, 2, 3, 4, 5]
    best_kc = [min((r for r in rows if r["K_C"] == kc and r["model"] != "M_C"),
                   key=lambda r: r["CV_RMSE_um"]) for kc in kc_values]
    axes[0].errorbar(range(6), [r["CV_RMSE_um"] for r in best_kc],
                     yerr=[r["CV_SE_um"] for r in best_kc], fmt="o-", color=BLUE,
                     lw=1, ms=4, capsize=2, label="fixed knots")
    axes[0].errorbar([6], [chosen_c["CV_RMSE_um"]], yerr=[chosen_c["CV_SE_um"]],
                     fmt="D", color=GREEN, ms=5, capsize=2,
                     label=rf"penalised, $\lambda$={float(chosen_c['lambda']):.2g}")
    axes[0].set_xticks(range(7),
                       ["stationary", "1", "2", "3", "4", "5", "penalised"],
                       rotation=35)
    axes[0].set(xlabel=r"coefficient knots $K_C$", ylabel="interior CV-RMSE [µm]")
    panel_title(axes[0], "a", "Varying-coefficient flexibility")
    axes[0].legend(frameon=False, fontsize=6.6, loc="upper right")
    for kind, colour in (("M_B", GREY), ("M_C", BLUE)):
        values = []
        for ks in range(9):
            subset = [r for r in rows if r["model"] == kind and r["K_s"] == ks]
            values.append(min(subset, key=lambda r: r["CV_RMSE_um"])["CV_RMSE_um"])
        axes[1].plot(range(9), values, "o-", ms=3.5, lw=1, color=colour, label=kind)
    axes[1].set(xlabel=r"slow-trend knots $K_s$", ylabel="interior CV-RMSE [µm]")
    panel_title(axes[1], "b", "Slow-trend flexibility"); axes[1].legend(frameon=False)
    for kind, colour in (("M_B", GREY), ("M_C", BLUE)):
        values = []
        for h in (2, 3, 4):
            subset = [r for r in rows if r["model"] == kind and r["H"] == h]
            values.append(min(subset, key=lambda r: r["CV_RMSE_um"])["CV_RMSE_um"])
        axes[2].plot((2, 3, 4), values, "o-", ms=3.5, lw=1, color=colour, label=kind)
    axes[2].set_xticks((2, 3, 4)); axes[2].set(xlabel="maximum harmonic H", ylabel="interior CV-RMSE [µm]")
    panel_title(axes[2], "c", "Harmonic order")
    for ax in axes: finish_axis(ax)
    fig.text(.07, .045,
             "Appendix A2. Blocked 10-fold CV repeated over five shifted sets of "
             "block edges; only interpolation folds are scored.\n"
             rf"M_C carries a rich $K_C$={PEN_KNOTS} basis with a second-difference "
             r"penalty, so $\lambda$ sets the flexibility; (b) and (c) are that "
             "penalised family.\n"
             f"One-SE choices: M_B H={chosen_b['H']}, Ks={chosen_b['K_s']}; "
             f"M_C H={chosen_c['H']}, Ks={chosen_c['K_s']}, "
             rf"$\lambda$={float(chosen_c['lambda']):.2g} "
             f"({chosen_c['effective_parameters']:.0f} effective parameters).",
             fontsize=7.1, color=INK2, linespacing=1.25)
    save(fig, "e2_appendix_model_selection", OUT)


def appendix_residuals(data: dict, model: dict) -> None:
    x, theta = data["x"], data["theta"]
    residual = data["D2"]["ybar"] - predict(model, x, theta)
    time = data["D2"]["time_s"] - np.mean(data["D2"]["time_s"])
    frequency, power = residual_periodogram(x, residual, .05, 2.2)
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.1))
    fig.subplots_adjust(left=.09, right=.97, top=.94, bottom=.11, hspace=.34, wspace=.28)
    axes[0, 0].plot(x, residual, ".", ms=2.3, color=GREY)
    running = running_profile(x, residual, 2.0, .5)
    axes[0, 0].plot(running["x"], running["mean"], color=RUST, lw=1.2)
    axes[0, 0].set(xlabel="position [mm]", ylabel="residual [µm]")
    panel_title(axes[0, 0], "a", "Residual against position")
    phase = np.degrees(theta)
    axes[0, 1].plot(phase, residual, ".", ms=2.3, color=GREY)
    axes[0, 1].set(xlabel="absolute screw angle [°]", ylabel="residual [µm]")
    panel_title(axes[0, 1], "b", "Residual folded on one screw turn")
    axes[1, 0].plot(time, residual, ".", ms=2.3, color=GREY)
    axes[1, 0].set(xlabel="paired-visit mean time − mean [s]", ylabel="residual [µm]")
    panel_title(axes[1, 0], "c", "Residual against visit time")
    axes[1, 1].plot(frequency, power, color=BLUE, lw=1)
    for k in range(1, 5): axes[1, 1].axvline(k/2, color=LIGHT, lw=.6, ls="--")
    axes[1, 1].set(xlabel="spatial frequency [cycles/mm]", ylabel="normalised power")
    panel_title(axes[1, 1], "d", "Residual Lomb–Scargle periodogram")
    for ax in axes.ravel():
        if ax is not axes[1, 1]:
            ax.axhline(0, color=GREY, lw=.45)
        finish_axis(ax)
    fig.text(.09, .03, "Appendix A3. Residual checks for the selected E2 model under D2.",
             fontsize=7.3, color=INK2)
    save(fig, "e2_appendix_residual_diagnostics", OUT)


def appendix_exploratory(data: dict, sliding: list[dict], model: dict,
                         windows4: dict) -> dict:
    x = data["x"]
    comp = components(model, x)
    base = model["beta"][0] + comp["gain"] + comp["slow"]
    remainder = data["D2"]["ybar"] - base
    f1, p1 = residual_periodogram(x, remainder, .40, .60, 4000)
    f2, p2 = residual_periodogram(x, remainder, .90, 1.10, 4000)
    sx = np.asarray([row["x"] for row in sliding])
    fig, axes = plt.subplots(1, 2, figsize=(8.1, 3.65))
    fig.subplots_adjust(left=.09, right=.96, top=.89, bottom=.20, wspace=.30)
    axes[0].plot(f1, p1, color=BLUE, lw=1.1, label="near k = 1")
    axes[0].plot(f2, p2, color=RUST, lw=1.1, label="near k = 2")
    axes[0].axvline(.5, color=GREY, ls="--", lw=.7)
    axes[0].axvline(1.0, color=GREY, ls="--", lw=.7)
    axes[0].set(xlabel="spatial frequency [cycles/mm]", ylabel="normalised power")
    panel_title(axes[0], "a", "Fine screw-frequency search")
    axes[0].legend(frameon=False)
    for k, colour in ((1, BLUE), (2, RUST)):
        values = np.asarray([row[f"C{k}"] for row in sliding])
        axes[1].plot(values.real, values.imag, "-", color=colour, lw=1.1,
                     label=rf"E2 $C_{k}(x)$")
        points = axes[1].scatter(values.real, values.imag, c=sx, cmap=POSITION_CMAP,
                                 norm=POSITION_NORM, s=8, zorder=3)
        for window in windows4.values():
            c = window["spectrum"][k-1]["complex"]
            axes[1].plot(c.real, c.imag, "o", ms=3.8, color=colour)
    axes[1].axhline(0, color=GREY, lw=.45); axes[1].axvline(0, color=GREY, lw=.45)
    axes[1].set(xlabel="real coefficient [µm]", ylabel="imaginary coefficient [µm]")
    panel_title(axes[1], "b", "Sliding harmonic phasor paths")
    axes[1].legend(frameon=False, loc="best")
    cbar = fig.colorbar(points, ax=axes[1], pad=.025, fraction=.055)
    cbar.set_label("window centre [mm]"); cbar.outline.set_linewidth(.5)
    for ax in axes: finish_axis(ax)
    fig.text(.09, .055,
             "Appendix A4. Exploratory beat test. Filled circles in the phasor panel are E4 windows.",
             fontsize=7.2, color=INK2)
    save(fig, "e2_appendix_beat_test", OUT)
    peaks = []
    for label, frequency, power in (("k1", f1, p1), ("k2", f2, p2)):
        index = int(np.argmax(power)); peaks.append({"band": label,
                                                     "frequency_per_mm": float(frequency[index]),
                                                     "power": float(power[index])})
    return {"peaks": peaks,
            "interpretation": ("The dominant k=1 peak is unresolved from exactly "
                               "0.5 cycles/mm at the approximately 0.009 cycles/mm "
                               "span resolution; no resolved beat sideband is supported.")}


def appendix_drift(data: dict, selected_row: dict) -> dict:
    grid = np.linspace(data["x"].min(), data["x"].max(), 500)
    slow, bprofiles, models = {}, {}, {}
    for name in ("D0", "D1", "D2"):
        models[name] = fit_selected(data, selected_row, name, robust=False)
        slow[name] = components(models[name], grid)["slow"]
        bprofiles[name] = running_profile(data["x"], data[name]["B"], 6.0)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.5))
    fig.subplots_adjust(left=.09, right=.97, top=.89, bottom=.20, wspace=.28)
    for name in ("D0", "D1", "D2"):
        axes[0].plot(grid, slow[name], color=DRIFT_COLOURS[name], lw=1.1, label=name)
        axes[1].plot(bprofiles[name]["x"], bprofiles[name]["mean"],
                     color=DRIFT_COLOURS[name], lw=1.1, label=name)
    axes[0].set(xlabel="position [mm]", ylabel="slow trend [µm]")
    axes[1].set(xlabel="position [mm]", ylabel="6 mm mean reversal [µm]")
    panel_title(axes[0], "a", "Slow geometric trend")
    panel_title(axes[1], "b", "Reversal profile")
    axes[0].legend(frameon=False); axes[1].legend(frameon=False)
    for ax in axes: endpoint_bands(ax, labels=False); finish_axis(ax)
    fig.text(.09, .055, "Appendix A5. Drift-correction sensitivity; D2 is the default.",
             fontsize=7.2, color=INK2)
    save(fig, "e2_appendix_drift_sensitivity", OUT)
    slow_stack = np.vstack([slow[name] for name in ("D0", "D1", "D2")])
    b_stack = np.vstack([bprofiles[name]["mean"] for name in ("D0", "D1", "D2")])
    return {"max_slow_difference_um": float(np.max(np.ptp(slow_stack, axis=0))),
            "max_reversal_difference_um": float(np.max(np.ptp(b_stack, axis=0)))}


def appendix_iso_phase(v3: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    fig.subplots_adjust(left=.12, right=.97, top=.89, bottom=.22)
    phase = np.asarray([row["electrical_phase"] for row in v3])
    residual = np.asarray([row["G_B_mean_um"] for row in v3])
    ax.plot(phase, residual, "D", ms=4, color=BLUE)
    ax.axhline(0, color=GREY, lw=.55)
    ax.set(xlabel="electrical phase [microsteps modulo 64]",
           ylabel="E0 residual after G + B [µm]")
    ax.set_title("E0 residual against motor electrical phase", loc="left", pad=4)
    finish_axis(ax)
    fig.text(.12, .055,
             "Appendix A6. E2 cannot identify this motor-phase term because all E2 targets lie on the 64-microstep grid.",
             fontsize=7.2, color=INK2)
    save(fig, "e2_appendix_e0_electrical_phase", OUT)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows: return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    style(); ANALYSIS.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    run = e2_run.load(); quality = record_quality(run); data = preprocess(run, quality)
    serial = parse_serial_checks()
    assert np.all(data["usteps"] % 64 == 0), "E2 target escaped the 64-microstep grid"
    phase_counts = np.bincount(data["phase_bin"], minlength=50)[:50]
    cv_rows, chosen_b_row, chosen_c_row = blocked_cv(data)
    model_b = fit_selected(data, chosen_b_row)
    model_c = fit_selected(data, chosen_c_row)
    windows1, windows4 = external_windows()
    v2 = validate_windows(model_b, model_c, windows1, windows4)
    v2_pass = float(np.mean([row["rms_M_C_um"] for row in v2])) < float(
        np.mean([row["rms_M_B_um"] for row in v2]))
    cv_win = chosen_b_row["CV_RMSE_um"] - chosen_c_row["CV_RMSE_um"]
    cv_required = max(chosen_b_row["CV_SE_um"], chosen_c_row["CV_SE_um"])
    choose_mc = cv_win > cv_required and v2_pass
    selected_row = chosen_c_row if choose_mc else chosen_b_row
    selected = model_c if choose_mc else model_b
    reversal = {name: reversal_fit(data, name) for name in ("D0", "D1", "D2")}
    sliding = sliding_profile(data)
    reversal_sliding = sliding_profile(data, field="B")
    v3, iso = validate_iso(selected, reversal["D2"])

    ref_mscnt = serial["reference_mscnt"]
    validity = [
        {"check": "planned/matched stops", "value": f"{run.matched.sum()}/{run.plan.target_mm.size}",
         "criterion": "807/807", "pass": bool(run.matched.all())},
        {"check": "paired map targets", "value": len(data["x"]), "criterion": 400,
         "pass": len(data["x"]) == 400 and not data["missing"]},
        {"check": "record std pass fraction", "value": float(np.mean(quality["std_um"] <= STD_LIMIT_UM)),
         "criterion": "limit 0.5 um", "pass": bool(np.all(quality["std_um"] <= STD_LIMIT_UM))},
        {"check": "record slope-effect pass fraction",
         "value": float(np.mean(np.abs(quality["slope_um_s"])*config["record_s"] <= SLOPE_EFFECT_LIMIT_UM)),
         "criterion": "limit 0.25 um", "pass": bool(np.all(np.abs(quality["slope_um_s"])*config["record_s"] <= SLOPE_EFFECT_LIMIT_UM))},
        {"check": "saturated record samples", "value": int(quality["saturated"].sum()),
         "criterion": 0, "pass": int(quality["saturated"].sum()) == 0},
        {"check": "reference MSCNT", "value": np.unique(ref_mscnt).tolist(),
         "criterion": "identical at all 7 references", "pass": bool(np.all(ref_mscnt == ref_mscnt[0]))},
        {"check": "64-microstep grid", "value": int(np.count_nonzero(data["usteps"] % 64)),
         "criterion": 0, "pass": bool(np.all(data["usteps"] % 64 == 0))},
        {"check": "phase-bin minimum", "value": int(phase_counts.min()),
         "criterion": ">=5", "pass": bool(phase_counts.min() >= 5)},
        {"check": "DATUM to END_REF drift", "value": float(data["anchor_error_um"][-1]-data["anchor_error_um"][0]),
         "criterion": "abs <=0.5 um", "pass": bool(abs(data["anchor_error_um"][-1]-data["anchor_error_um"][0]) <= .5)},
        {"check": "firmware completion", "value": serial["summary"],
         "criterion": "complete=1, aborted=0", "pass": "complete=1" in serial["summary"] and "aborted=0" in serial["summary"]},
        {"check": "environment", "value": "not recorded", "criterion": "report availability", "pass": None},
    ]

    figure_f1(data, selected)
    figure_f1b(data, selected, v2)
    figure_f2(sliding, model_c, windows1, windows4)
    figure_f3(v2, v3, iso)
    figure_f4(data, reversal["D2"], windows4)
    appendix_validity(data, quality, validity)
    appendix_cv(cv_rows, chosen_b_row, chosen_c_row)
    appendix_residuals(data, selected)
    exploratory = appendix_exploratory(data, sliding, model_c, windows4)
    sensitivity = appendix_drift(data, selected_row)
    appendix_iso_phase(v3)

    target_rows = []
    for i, x in enumerate(data["x"]):
        row = {"target_mm": x, "target_usteps": int(data["usteps"][i]),
               "screw_angle_deg": float(np.degrees(data["theta"][i])),
               "phase_bin": int(data["phase_bin"][i]),
               "up_time_s": data["up_time_s"][i], "down_time_s": data["down_time_s"][i]}
        for variant in ("D0", "D1", "D2"):
            row |= {f"e_up_{variant}_um": data[variant]["e_up"][i],
                    f"e_down_{variant}_um": data[variant]["e_down"][i],
                    f"ybar_{variant}_um": data[variant]["ybar"][i],
                    f"B_{variant}_um": data[variant]["B"][i]}
        target_rows.append(row)
    write_csv(ANALYSIS / "e2_preprocessed_targets.csv", target_rows)
    write_csv(ANALYSIS / "e2_model_selection.csv", cv_rows)
    write_csv(ANALYSIS / "e2_v2_window_validation.csv", v2)
    write_csv(ANALYSIS / "e2_v3_iso_prediction.csv", v3)
    write_csv(ANALYSIS / "e2_validity_table.csv", validity)
    reversal_sliding_rows = []
    for row in reversal_sliding:
        reversal_sliding_rows.append({
            key: (value.real if key.endswith("_real") else value)
            for key, value in row.items() if key not in ("C1", "C2")})
        reversal_sliding_rows[-1] |= {
            "C1_real_um": row["C1"].real, "C1_imag_um": row["C1"].imag,
            "C2_real_um": row["C2"].real, "C2_imag_um": row["C2"].imag,
        }
    write_csv(ANALYSIS / "e2_reversal_sliding_profile.csv", reversal_sliding_rows)

    step_a = json.loads(STEP_A.read_text(encoding="utf-8"))["fit_uncorrected"]
    profile_summary = {}
    grid = np.arange(data["x"].min(), data["x"].max()+.001, .5)
    for k in (1, 2):
        profile = harmonic_profile(model_c, grid, k)
        profile_summary[str(k)] = {"amplitude_min_um": float(profile["amplitude"].min()),
                                   "amplitude_max_um": float(profile["amplitude"].max())}
    reversal_summary = {}
    for name, result in reversal.items():
        reversal_summary[name] = {
            "mean_um": result["mean_um"], "mean_sigma_um": result["mean_sigma_um"],
            "slope_um_per_100mm": result["slope_um_per_100mm"],
            "slope_sigma_um_per_100mm": result["slope_sigma_um_per_100mm"],
        }
    slope_robust = (len({np.sign(result["slope_um_per_100mm"])
                         for result in reversal.values()}) == 1 and
                    all(abs(result["slope_um_per_100mm"]) >
                        2*result["slope_sigma_um_per_100mm"]
                        for result in reversal.values()))
    region95 = np.abs(data["x"] - 95.0) <= 3.0
    reversal_summary["robust_trend_across_D0_D1_D2"] = slope_robust
    reversal_summary["region_95mm"] = {
        "mean_um": float(data["D2"]["B"][region95].mean()),
        "sigma_of_mean_um": float(data["D2"]["B"][region95].std(ddof=1) /
                                  np.sqrt(region95.sum())),
        "target_count": int(region95.sum()),
        "interpretation": "within the continuous E2 reversal profile; no isolated anomaly",
    }
    reversal_amp = np.asarray([row["A1"] for row in reversal_sliding])
    reversal_summary["sliding_A_B1_um"] = {
        "minimum": float(reversal_amp.min()), "maximum": float(reversal_amp.max()),
        "near_motor": float(np.mean(reversal_amp[:3])),
        "mid_travel": float(np.mean(reversal_amp[len(reversal_amp)//2-1:
                                                  len(reversal_amp)//2+2])),
        "near_support": float(np.mean(reversal_amp[-3:])),
    }
    map_stop = np.isin(run.plan.state, ("UP", "DOWN"))
    record_mean_noise = quality["std_um"][map_stop] / np.sqrt(quality["samples"][map_stop])
    summary = {
        "experiment": "E2_global_map", "scope_samples": int(run.t.size),
        "validity": validity,
        "drift": {"anchors_time_s": data["anchor_time_s"].tolist(),
                  "anchors_error_um": data["anchor_error_um"].tolist(),
                  **sensitivity},
        "model_selection": {"M_B_one_se": chosen_b_row, "M_C_one_se": chosen_c_row,
                            "M_C_CV_improvement_um": cv_win,
                            "required_one_se_um": cv_required,
                            "V2_M_C_better": v2_pass,
                            "selected": "M_C" if choose_mc else "M_B"},
        "selected_model": {"H": selected_row["H"], "K_s": selected_row["K_s"],
                           "K_C": selected_row["K_C"],
                           "lambda": selected_row.get("lambda", ""),
                           "effective_parameters": selected_row["effective_parameters"],
                           "parameters": selected_row["parameters"],
                           "residual_rms_um": selected["rms"],
                           "condition_number_scaled": selected["condition_scaled"],
                           "max_abs_slow_vs_C_column_correlation": slow_coefficient_correlation(selected),
                           "huber_coefficients_over_1sigma": selected["huber_over_sigma_count"],
                           "huber_max_shift_sigma": selected["huber_over_sigma_max"],
                           "median_record_mean_noise_um": float(np.median(record_mean_noise)),
                           "mean_external_local_floor_um": float(np.mean(
                               [row["rms_local_um"] for row in v2]))},
        "gain": {"E2_ppm": float(selected["beta"][1]*1000),
                 "E2_sigma_ppm": float(np.sqrt(selected["cov"][1, 1])*1000),
                 "Step_A_ppm": step_a["slope_ppm"],
                 "Step_A_ci95_ppm": step_a["slope_ci95_ppm"]},
        "harmonic_profile": profile_summary,
        "reversal": reversal_summary,
        "validation": {"V2": v2, "V3_ISO_metrics_um": iso},
        "exploratory_beat_test": exploratory,
        "limitations": {"environment": "not recorded",
                        "endpoint_locations": "inferred from homed travel",
                        "E0_microstep": "E0 residual includes motor electrical-phase error not identifiable by E2"},
        "convention": {"deviation": "positive means measured travel exceeds command",
                       "theta": "2 pi (commanded microsteps mod 3200) / 3200",
                       "phase": "e = A cos(k theta + phi)",
                       "pattern_shift": "-(phi_k - mean(phi_k))/k"},
    }
    (ANALYSIS / "e2_findings_summary.json").write_text(
        json.dumps(summary, indent=2, default=lambda value: value.tolist() if isinstance(value, np.ndarray) else float(value)) + "\n",
        encoding="utf-8")

    print(f"E2: {len(data['x'])} paired targets; reference drift "
          f"{data['anchor_error_um'][-1]-data['anchor_error_um'][0]:+.3f} um")
    print(f"M_B one-SE: H={chosen_b_row['H']}, Ks={chosen_b_row['K_s']}, "
          f"CV={chosen_b_row['CV_RMSE_um']:.3f} um")
    print(f"M_C one-SE: H={chosen_c_row['H']}, Ks={chosen_c_row['K_s']}, "
          f"Kc={chosen_c_row['K_C']}, lambda={chosen_c_row['lambda']}, "
          f"edf={chosen_c_row['effective_parameters']:.1f}, "
          f"CV={chosen_c_row['CV_RMSE_um']:.3f} um")
    print(f"selected {'M_C' if choose_mc else 'M_B'}; residual rms {selected['rms']:.3f} um; "
          f"gain {selected['beta'][1]*1000:+.2f} ppm")
    print(f"ISO prediction E {iso['raw']['E']:.2f} -> {iso['G_plus_reversal']['E']:.2f} um; "
          f"M {iso['raw']['M']:.2f} -> {iso['G_plus_reversal']['M']:.2f} um")
    for directory in (PRESENTATION, OUT):
        for path in sorted(directory.glob("e2_*.png")):
            print("wrote", path.relative_to(HERE))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
