#!/usr/bin/env python3
"""Process the 30 Oct repeat session against E2 and identify the pooled map.

This follows E3_processing_instructions.md section by section, with one
substantive departure forced by the data: the session that landed is a repeat
of E2 at the *same* targets, not the new-seed session the procedure specifies.
The firmware was re-run unmodified, so it re-derived seed 5315 and emitted the
same 400 positions (sequence/PROVENANCE.md). Every question that asks whether a
stage property reproduces between days is answerable, and answerable more
sharply than planned, because the two sessions pair target by target and no
interpolation stands between them. The one question that is not answerable is
Q1 as written - whether the map predicts *new positions* - and the summary says
so rather than quietly answering the easier question instead.

All of the modelling machinery is imported from render_e2_findings: the design
builder, the penalised fit, the harmonic profile, the blocked CV, the E1/E4
window validation and the ISO 230-2 projection. Reducing the two sessions with
one implementation is the whole basis of a day-to-day comparison.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.stats import f as f_dist

HERE = Path(__file__).resolve().parent
E2_DIR = HERE.parent / "e2_global_map"
ISO_DIR = HERE.parents[2] / "iso230-2" / "analysis"
for directory in (E2_DIR, ISO_DIR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

import e2_run                                   # noqa: E402
import e3_run                                   # noqa: E402
import render_e2_findings as e2f                 # noqa: E402

ANALYSIS = HERE / "analysis"
CONFIG = E2_DIR / "e2_config.json"
E2_SUMMARY = E2_DIR / "analysis" / "e2_findings_summary.json"
V3_ISO = ISO_DIR / "iso230_2_v3_results_blind.json"
V4_ISO = ISO_DIR / "iso230_2_v4_results_blind.json"

WELL_SUPPORTED = (6.0, 114.0)   # where the varying coefficients are determined
GRID_STEP = 1.0                 # the 1 mm comparison grid of section 5.3
E4_FLOOR_UM = 0.24              # E4 median local noise floor
BLOCK_MM = 10.0                 # per-block reporting width
SHIFT_MIN_MM = 10.0             # circular-shift null, minimum shift
N_PERM = 2000
STEEL_PPM_PER_K = 11.5
RNG_SEED = 0
N_BOOTSTRAP = 400


# ---------------------------------------------------------------- sessions

def reduce_session(module, label: str) -> dict:
    """Alignment, record quality, per-target preprocessing and serial checks."""
    e2f.LOG = module.LOG if hasattr(module, "LOG") else e2f.LOG
    run = module.load()
    quality = e2f.record_quality(run)
    data = e2f.preprocess(run, quality)
    serial = e2f.parse_serial_checks()
    data["label"] = label
    return {"run": run, "quality": quality, "data": data, "serial": serial}


def basis_domain(data: dict) -> tuple[float, float]:
    return float(data["x"].min()), float(data["x"].max())


# ---------------------------------------------------------------- validity

def validity_rows(session: dict, config: dict, domain: tuple[float, float],
                  iso_usteps: np.ndarray, reference: dict | None) -> list[dict]:
    run, quality, data = session["run"], session["quality"], session["data"]
    serial = session["serial"]
    phase_counts = np.bincount(data["phase_bin"], minlength=50)[:50]
    ref_mscnt = serial["reference_mscnt"]
    slope_effect = np.abs(quality["slope_um_s"]) * config["record_s"]
    drift = float(data["anchor_error_um"][-1] - data["anchor_error_um"][0])
    rows = [
        {"check": "planned/matched stops",
         "value": f"{int(run.matched.sum())}/{run.plan.target_mm.size}",
         "criterion": "807/807", "pass": bool(run.matched.all())},
        {"check": "paired map targets", "value": len(data["x"]),
         "criterion": 400,
         "pass": len(data["x"]) == 400 and not data["missing"]},
        {"check": "record std pass fraction",
         "value": float(np.mean(quality["std_um"] <= e2f.STD_LIMIT_UM)),
         "criterion": "limit 0.5 um",
         "pass": bool(np.all(quality["std_um"] <= e2f.STD_LIMIT_UM))},
        {"check": "record slope-effect pass fraction",
         "value": float(np.mean(slope_effect <= e2f.SLOPE_EFFECT_LIMIT_UM)),
         "criterion": "limit 0.25 um",
         "pass": bool(np.all(slope_effect <= e2f.SLOPE_EFFECT_LIMIT_UM))},
        {"check": "saturated record samples",
         "value": int(quality["saturated"].sum()), "criterion": 0,
         "pass": int(quality["saturated"].sum()) == 0},
        {"check": "reference MSCNT", "value": np.unique(ref_mscnt).tolist(),
         "criterion": "identical at all 7 references",
         "pass": bool(np.all(ref_mscnt == ref_mscnt[0]))},
        {"check": "64-microstep grid",
         "value": int(np.count_nonzero(data["usteps"] % 64)), "criterion": 0,
         "pass": bool(np.all(data["usteps"] % 64 == 0))},
        {"check": "phase-bin minimum", "value": int(phase_counts.min()),
         "criterion": ">=5", "pass": bool(phase_counts.min() >= 5)},
        {"check": "DATUM to END_REF drift", "value": drift,
         "criterion": "abs <=0.5 um", "pass": bool(abs(drift) <= .5)},
        {"check": "firmware completion", "value": serial["summary"],
         "criterion": "complete=1, aborted=0",
         "pass": "complete=1" in serial["summary"]
                 and "aborted=0" in serial["summary"]},
        {"check": "recording span [min]",
         "value": float(run.t[-1] / 60.0), "criterion": "report", "pass": None},
        {"check": "environment", "value": "not recorded",
         "criterion": "report availability", "pass": None},
    ]
    # Section 2 additions. Independence from E0 is a property of the target
    # list, so it is re-checked rather than inherited: the exclusion was
    # applied by e2_sequence, and a check that only trusts the generator would
    # not notice a plan that had been edited by hand afterwards.
    near_e0 = int(np.count_nonzero(
        np.min(np.abs(data["usteps"][:, None] - iso_usteps[None, :]), axis=1) <= 64))
    rows.append({"check": "independence from E0 (within +-64 usteps)",
                 "value": near_e0, "criterion": 0, "pass": near_e0 == 0})
    rows.append({"check": "basis span used [mm]",
                 "value": f"{domain[0]:.3f}..{domain[1]:.3f}",
                 "criterion": "identical in both sessions", "pass": None})
    if reference is not None:
        overlap = int(np.count_nonzero(
            np.min(np.abs(data["usteps"][:, None]
                          - reference["usteps"][None, :]), axis=1) <= 64))
        rows.append({"check": "overlap with E2 (within +-64 usteps)",
                     "value": f"{overlap}/{len(data['x'])}",
                     "criterion": "a few expected from independent seeds",
                     "pass": None})
        identical = bool(np.array_equal(data["usteps"], reference["usteps"]))
        rows.append({"check": "target list identical to E2",
                     "value": identical,
                     "criterion": "False for a new-seed session",
                     "pass": not identical})
    return rows


def iso_target_usteps(config: dict) -> np.ndarray:
    """E0 target positions in this experiment's microsteps.

    The ISO sequence stores its commands as cmd_u256, sixteenths of the mres=16
    microstep the map is gridded on, so the comparison has to come down to the
    coarser unit rather than the other way round.
    """
    path = (E2_DIR / config["iso_target_source"]).resolve()
    rows = list(csv.DictReader(path.open()))
    return np.asarray([int(round(int(row["cmd_u256"]) / 16)) for row in rows])


# ---------------------------------------------------------------- section 3

def drift_sensitivity(data: dict) -> dict:
    """How much the three drift variants move ybar and B, per section 3.

    The mirrored pass design makes the mean visit time of a target constant
    along the travel: going up, target i is reached after the gaps below it;
    coming down, after the gaps above it, and the two sum to the same total
    for every i. A drift that is *linear* in time therefore adds the same
    constant to every ybar and cancels when the mean is removed - which is why
    D0 and D1 return identical ybar shapes. Only the drift's curvature reaches
    ybar, so D1 against D2 measures exactly that and nothing else.

    B is the opposite case. It differs the two directions, so it carries the
    full drift accumulated over the up-to-down gap, and that gap runs from one
    minute at the turnaround to the whole session at the ends.
    """
    ybar = np.vstack([data[name]["ybar"] for name in ("D0", "D1", "D2")])
    reversal = np.vstack([data[name]["B"] for name in ("D0", "D1", "D2")])
    curvature = data["D1"]["ybar"] - data["D2"]["ybar"]
    times = data["anchor_time_s"]
    errors = data["anchor_error_um"]
    segments = [{"from_s": float(times[i]), "to_s": float(times[i + 1]),
                 "change_um": float(errors[i + 1] - errors[i]),
                 "rate_um_per_h": float((errors[i + 1] - errors[i])
                                        / (times[i + 1] - times[i]) * 3600.0)}
                for i in range(times.size - 1)]
    gap = data["down_time_s"] - data["up_time_s"]
    return {"anchor_time_s": times.tolist(), "anchor_error_um": errors.tolist(),
            "segments": segments,
            "total_change_um": float(errors[-1] - errors[0]),
            "mean_visit_time_spread_s": float(np.ptp(
                0.5 * (data["up_time_s"] + data["down_time_s"]))),
            "up_to_down_gap_s": [float(gap.min()), float(gap.max())],
            "max_ybar_spread_um": float(np.max(np.ptp(ybar, axis=0))),
            "max_reversal_spread_um": float(np.max(np.ptp(reversal, axis=0))),
            "D0_vs_D1_max_ybar_um": float(np.max(np.abs(
                data["D0"]["ybar"] - data["D1"]["ybar"]))),
            "D1_vs_D2_max_ybar_um": float(np.max(np.abs(curvature))),
            "D1_vs_D2_peak_at_x_mm": float(data["x"][np.argmax(np.abs(curvature))]),
            "D3_available": False,
            "D3_note": "the firmware emits DATUM, TURN_REF and END_REF only, "
                       "so D2 already runs through every reference visit"}


def window_quality(session: dict, config: dict) -> dict:
    """Record-window std and slope for each paired target, both directions.

    e2f.preprocess reduces the plan to paired targets but does not carry the
    per-window quality across, so the row mapping is rebuilt here the same way
    it builds it. One failing window in 807 is the difference between a
    reported flag and a 1.8 um outlier sitting in the fit.
    """
    run, quality = session["run"], session["quality"]
    plan = run.plan
    up = (plan.state == "UP") & run.matched
    down = (plan.state == "DOWN") & run.matched
    _, iu, idn = np.intersect1d(plan.index[up], plan.index[down],
                                return_indices=True)
    up_rows = np.flatnonzero(up)[iu]
    down_rows = np.flatnonzero(down)[idn]
    order = np.argsort(plan.target_mm[up_rows])
    up_rows, down_rows = up_rows[order], down_rows[order]
    out = {}
    for name, rows in (("up", up_rows), ("down", down_rows)):
        out[f"std_{name}_um"] = quality["std_um"][rows]
        out[f"slope_effect_{name}_um"] = (np.abs(quality["slope_um_s"][rows])
                                          * config["record_s"])
        out[f"samples_{name}"] = quality["samples"][rows]
    bad = ((out["std_up_um"] > e2f.STD_LIMIT_UM)
           | (out["std_down_um"] > e2f.STD_LIMIT_UM)
           | (out["slope_effect_up_um"] > e2f.SLOPE_EFFECT_LIMIT_UM)
           | (out["slope_effect_down_um"] > e2f.SLOPE_EFFECT_LIMIT_UM))
    out["flagged"] = bad
    out["keep"] = ~bad
    return out


# ---------------------------------------------------------------- section 5

def fit_with_settings(data: dict, row: dict, domain: tuple[float, float],
                      variant: str = "D2", robust: bool = True) -> dict:
    kc = None if row["K_C"] == "stationary" else int(row["K_C"])
    lam = float(row["lambda"]) if row.get("lambda", "") != "" else None
    return e2f.fit_linear(data["x"], data["theta"], data[variant]["ybar"],
                          int(row["H"]), int(row["K_s"]), kc, domain,
                          robust=robust, penalty=lam)


def cross_residual(model: dict, data: dict, variant: str = "D2") -> np.ndarray:
    """Apply a model fitted elsewhere: no refit, one constant offset only."""
    residual = (data[variant]["ybar"]
                - e2f.predict(model, data["x"], data["theta"]))
    return residual - residual.mean()


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values) ** 2)))


def block_edges(domain: tuple[float, float]) -> np.ndarray:
    lo = math.floor(domain[0] / BLOCK_MM) * BLOCK_MM
    hi = math.ceil(domain[1] / BLOCK_MM) * BLOCK_MM
    return np.arange(lo, hi + .1, BLOCK_MM)


def cross_prediction_table(target: dict, own_b: dict, own_c: dict,
                           other_b: dict, other_c: dict,
                           direction: str, domain: tuple[float, float]
                           ) -> list[dict]:
    """T2: residual rms per level, overall and per 10 mm block."""
    x = target["x"]
    raw = target["D2"]["ybar"] - target["D2"]["ybar"].mean()
    levels = {
        "raw": raw,
        "cross_M_B": cross_residual(other_b, target),
        "cross_M_C": cross_residual(other_c, target),
        "own_M_B": own_b["residual"] - own_b["residual"].mean(),
        "own_M_C": own_c["residual"] - own_c["residual"].mean(),
    }
    edges = block_edges(domain)
    rows = []
    for name, lo, hi in ([("overall", domain[0] - 1, domain[1] + 1)]
                         + [(f"{a:.0f}-{b:.0f} mm", a, b)
                            for a, b in zip(edges[:-1], edges[1:])]):
        mask = (x > lo) & (x <= hi) if name != "overall" else np.ones(x.size, bool)
        if mask.sum() < 3:
            continue
        row = {"direction": direction, "block": name, "n": int(mask.sum()),
               **{f"rms_{key}_um": rms(value[mask])
                  for key, value in levels.items()},
               "floor_um": E4_FLOOR_UM}
        for kind in ("M_B", "M_C"):
            gap = row["rms_raw_um"] - row[f"rms_own_{kind}_um"]
            row[f"gap_closed_{kind}"] = (
                (row["rms_raw_um"] - row[f"rms_cross_{kind}_um"]) / gap
                if gap > 1e-9 else float("nan"))
        rows.append(row)
    return rows


def coefficient_comparison(model_a: dict, model_b: dict, grid: np.ndarray
                           ) -> tuple[list[dict], dict]:
    """Section 5.3: gain, slow error and the varying harmonics, side by side."""
    gain_a = float(model_a["beta"][1] * 1000.0)
    gain_b = float(model_b["beta"][1] * 1000.0)
    sigma_a = float(np.sqrt(model_a["cov"][1, 1]) * 1000.0)
    sigma_b = float(np.sqrt(model_b["cov"][1, 1]) * 1000.0)
    delta = gain_b - gain_a
    gain = {"E2_ppm": gain_a, "E2_sigma_ppm": sigma_a,
            "E3_ppm": gain_b, "E3_sigma_ppm": sigma_b,
            "difference_ppm": delta,
            "z": delta / math.hypot(sigma_a, sigma_b),
            "equivalent_steel_dT_K": delta / STEEL_PPM_PER_K}
    comp_a = e2f.components(model_a, grid)
    comp_b = e2f.components(model_b, grid)
    rows = []
    profiles = {}
    for k in (1, 2):
        profiles[k] = (e2f.harmonic_profile(model_a, grid, k),
                       e2f.harmonic_profile(model_b, grid, k))
    shift_reference = {}
    for k, (pa, pb) in profiles.items():
        shift_reference[k] = (e2f.phase_reference(pa), e2f.phase_reference(pb))
    for i, x in enumerate(grid):
        row = {"x_mm": float(x),
               "slow_E2_um": float(comp_a["slow"][i]),
               "slow_E3_um": float(comp_b["slow"][i]),
               "slow_difference_um": float(comp_b["slow"][i] - comp_a["slow"][i]),
               "slow_difference_sigma_um": float(math.hypot(
                   comp_a["slow_sigma"][i], comp_b["slow_sigma"][i]))}
        for k, (pa, pb) in profiles.items():
            ref_a, ref_b = shift_reference[k]
            # Pattern shift, as the convention defines it: -(phi_k - mean)/k,
            # with each session's own mean phase as its reference.
            shift_a = -np.degrees(pa["phase"][i] - ref_a) / k
            shift_b = -np.degrees(pb["phase"][i] - ref_b) / k
            row |= {
                f"A{k}_E2_um": float(pa["amplitude"][i]),
                f"A{k}_E3_um": float(pb["amplitude"][i]),
                f"A{k}_difference_um": float(pb["amplitude"][i] - pa["amplitude"][i]),
                f"A{k}_difference_sigma_um": float(math.hypot(
                    pa["amplitude_sigma"][i], pb["amplitude_sigma"][i])),
                f"shift{k}_E2_deg": float(shift_a),
                f"shift{k}_E3_deg": float(shift_b),
                f"shift{k}_difference_deg": float(shift_b - shift_a),
                f"shift{k}_difference_sigma_deg": float(np.degrees(math.hypot(
                    pa["phase_sigma"][i], pb["phase_sigma"][i])) / k),
                f"C{k}_complex_difference_um": float(abs(
                    pb["complex"][i] - pa["complex"][i])),
            }
        rows.append(row)
    summary = {"gain": gain}
    for k in (1, 2):
        summary[f"C{k}"] = {
            "max_abs_amplitude_difference_um": max(
                abs(row[f"A{k}_difference_um"]) for row in rows),
            "max_abs_shift_difference_deg": max(
                abs(row[f"shift{k}_difference_deg"]) for row in rows),
            "max_complex_difference_um": max(
                row[f"C{k}_complex_difference_um"] for row in rows),
            "rms_complex_difference_um": rms(np.asarray(
                [row[f"C{k}_complex_difference_um"] for row in rows])),
        }
    summary["slow"] = {
        "max_abs_difference_um": max(abs(row["slow_difference_um"]) for row in rows),
        "max_difference_over_sigma": max(
            abs(row["slow_difference_um"]) / max(row["slow_difference_sigma_um"], 1e-9)
            for row in rows)}
    return rows, summary


def mscnt_prediction(serial_a: dict, serial_b: dict, mres: int = 16) -> dict:
    """The homing shift the driver's own phase counter implies.

    MSCNT indexes one electrical period - four full steps - in 1024 counts, so
    at mres = 16 one commanded microstep advances it by 16. Both sessions take
    every reference at the same absolute command, so a difference in MSCNT
    there is a difference in where homing put the microstep origin, in units
    the serial log reports directly and the harmonic phases know nothing about.
    Comparing it with delta1 and delta2 turns section 5.4 from a consistency
    check between two harmonics into a test against an independent measurement.
    """
    a = int(np.unique(serial_a["reference_mscnt"])[0])
    b = int(np.unique(serial_b["reference_mscnt"])[0])
    counts_per_ustep = 1024 // (4 * mres)
    difference = b - a
    # The counter wraps, so take the shorter way round.
    difference -= 1024 * round(difference / 1024)
    usteps = difference / counts_per_ustep
    return {"E2_mscnt": a, "E3_mscnt": b, "difference_counts": difference,
            "counts_per_ustep": counts_per_ustep,
            "usteps": usteps,
            "screw_deg": usteps / 3200.0 * 360.0,
            "travel_um": usteps / 1600.0 * 1000.0}


def homing_delta(model_a: dict, model_b: dict, grid: np.ndarray) -> dict:
    """Section 5.4: a homing shift rotates phi_k by k*delta for every k.

    Estimating delta from k = 1 and again from k = 2 is the whole test. A
    reference offset gives the same delta twice; a change in the screw's own
    harmonic content does not.
    """
    out = {}
    for k in (1, 2):
        pa = e2f.harmonic_profile(model_a, grid, k)
        pb = e2f.harmonic_profile(model_b, grid, k)
        weight = 1.0 / np.maximum(np.hypot(pa["phase_sigma"], pb["phase_sigma"]),
                                  np.radians(.2)) ** 2
        difference = pb["phase"] - pa["phase"]
        difference = difference - 2 * np.pi * np.round(difference / (2 * np.pi))
        mean = float(np.average(difference, weights=weight))
        sigma = float(np.sqrt(1.0 / weight.sum()))
        out[f"delta{k}_deg"] = np.degrees(mean) / k
        out[f"delta{k}_sigma_deg"] = np.degrees(sigma) / k
    out["delta_difference_deg"] = out["delta1_deg"] - out["delta2_deg"]
    out["consistent"] = bool(abs(out["delta_difference_deg"]) <= 2.0 * math.hypot(
        out["delta1_sigma_deg"], out["delta2_sigma_deg"]))
    out["delta1_travel_um"] = out["delta1_deg"] / 360.0 * 2000.0
    out["requirement_deg"] = 3.0
    out["within_requirement"] = bool(abs(out["delta1_deg"]) < 3.0
                                     and abs(out["delta2_deg"]) < 3.0)
    return out


# ---------------------------------------------------------------- section 6

def running_mean(x: np.ndarray, values: np.ndarray, grid: np.ndarray,
                 width: float) -> np.ndarray:
    out = np.full(grid.shape, np.nan)
    for j, point in enumerate(grid):
        mask = np.abs(x - point) <= width / 2
        if mask.sum() >= 3:
            out[j] = values[mask].mean()
    return out


def shift_null_correlation(a: np.ndarray, b: np.ndarray, step: float,
                           seed: int = RNG_SEED) -> tuple[float, float]:
    """Correlation with a circular-shift null, because neighbours correlate."""
    rho = float(np.corrcoef(a, b)[0, 1])
    rng = np.random.default_rng(seed)
    k_min = int(SHIFT_MIN_MM / step)
    if b.size <= 2 * k_min + 1:
        return rho, float("nan")
    null = np.asarray([
        np.corrcoef(a, np.roll(b, int(rng.integers(k_min, b.size - k_min))))[0, 1]
        for _ in range(N_PERM)])
    return rho, float(np.mean(np.abs(null) >= abs(rho)))


def residual_structure(d2: dict, d3: dict, residual2: np.ndarray,
                       residual3: np.ndarray, pooled2: np.ndarray,
                       pooled3: np.ndarray, pooled_edf: float) -> dict:
    """Section 6, with the paired test the repeat design makes available.

    The procedure prescribes a running-mean correlation on a common grid
    because two seeds put their targets in different places. Here they are the
    same places, so the residuals pair one to one and the correlation needs no
    smoothing at all - a stronger test, since smoothing both series is what
    forces the awkward circular-shift null in the first place. Both are
    reported: the paired one because it is the real answer, the smoothed one
    because it is what the procedure asks for and what a future new-seed
    session will have to use.
    """
    assert np.array_equal(d2["usteps"], d3["usteps"])
    out = {}

    # The model-free statement first. Two measurements of the same target on
    # two days differ by the session effects and the noise and by nothing
    # else, so their difference needs no model at all.
    y2, y3 = d2["D2"]["ybar"], d3["D2"]["ybar"]
    difference = y2 - y3
    out["model_free"] = {
        "rho_ybar": float(np.corrcoef(y2, y3)[0, 1]),
        "rms_ybar_E2_um": rms(y2), "rms_ybar_E3_um": rms(y3),
        "rms_difference_um": rms(difference),
        "day_to_day_per_session_um": rms(difference) / math.sqrt(2.0),
        "note": "the deviation at a target reproduces to this, before any "
                "model and including the sessions' smooth terms",
    }

    rho = float(np.corrcoef(residual2, residual3)[0, 1])
    n = residual2.size
    independent = rms(residual2 - residual3) / math.sqrt(2.0)
    shared = math.sqrt(max(rms(residual2) * rms(residual3) - independent ** 2,
                           0.0))
    out["paired_own_models"] = {
        "rho": rho, "n": int(n),
        # Fisher z: the targets are shared, the measurements are not.
        "fisher_z": 0.5 * math.log((1 + rho) / (1 - rho)) * math.sqrt(n - 3),
        "shared_variance_fraction": rho ** 2,
        "repeatable_rms_um": shared,
        "independent_rms_um": independent,
        "e4_local_floor_um": E4_FLOOR_UM,
        "smooth_terms_cost_um": math.sqrt(max(
            out["model_free"]["day_to_day_per_session_um"] ** 2
            - independent ** 2, 0.0)),
        "note": "each session fitted on its own, so the residuals share only "
                "physics; this is the test Q3 asks for",
    }

    # The same correlation under the pooled model is reported because it is
    # prescribed, and immediately discounted because it cannot answer the
    # question. Sharing one beta between the sessions makes the two residual
    # vectors structurally dependent: in every fitted direction each session's
    # residual carries minus half the other's noise, which biases the
    # correlation towards -p_eff/(2n) however geometric the truth is. The
    # measured value sits on that bias, so it says nothing further.
    rho_pooled = float(np.corrcoef(pooled2, pooled3)[0, 1])
    bias = -pooled_edf / (2.0 * n)
    out["paired_pooled_model"] = {
        "rho": rho_pooled, "n": int(n),
        "structural_bias": bias,
        "effective_parameters": pooled_edf,
        "consistent_with_bias_alone": bool(abs(rho_pooled - bias)
                                           < 2.0 / math.sqrt(n)),
        "note": "not a test of the physics: shared parameters induce this "
                "correlation. Use paired_own_models.",
    }
    grid = np.arange(WELL_SUPPORTED[0], WELL_SUPPORTED[1] + .01, 0.5)
    m2 = running_mean(d2["x"], residual2, grid, 2.0)
    m3 = running_mean(d3["x"], residual3, grid, 2.0)
    ok = ~np.isnan(m2) & ~np.isnan(m3)
    rho, p = shift_null_correlation(m2[ok], m3[ok], 0.5)
    out["running_mean_2mm"] = {"rho": rho, "p_circular_shift": p,
                               "grid_points": int(ok.sum())}
    for label, data, residual in (("E2", d2, residual2), ("E3", d3, residual3)):
        up = data["D2"]["e_up"] - data["D2"]["e_up"].mean()
        down = data["D2"]["e_down"] - data["D2"]["e_down"].mean()
        out[f"up_down_{label}"] = {
            "rho": float(np.corrcoef(up, down)[0, 1]),
            "note": "correlated means geometry, uncorrelated means time or air"}
    frequency, power = e2f.residual_periodogram(
        np.r_[d2["x"], d3["x"]], np.r_[pooled2, pooled3], 0.02, 0.5)
    peak = int(np.argmax(power))
    out["pooled_periodogram"] = {
        "peak_cycles_per_mm": float(frequency[peak]),
        "peak_period_mm": float(1.0 / frequency[peak]),
        "peak_power": float(power[peak]),
        "band": "0.02-0.5 cycles/mm"}
    return out


# ---------------------------------------------------------------- section 7

def session_design(data: dict, state: dict) -> np.ndarray:
    design, _ = e2f.prepare_design(data["x"], data["theta"], state["H"],
                                   state["Ks"], state["Kc"], state["domain"],
                                   state=state)
    return design


def pooled_system(datas: list[dict], state: dict, variant: str = "D2",
                  field: str = "ybar") -> tuple[np.ndarray, np.ndarray]:
    """Stack sessions with one offset column each; everything else shared."""
    designs = [session_design(data, state) for data in datas]
    counts = [data["x"].size for data in datas]
    offsets = np.zeros((sum(counts), len(datas)))
    start = 0
    for j, count in enumerate(counts):
        offsets[start:start + count, j] = 1.0
        start += count
    shared = np.vstack([design[:, 1:] for design in designs])
    y = np.concatenate([data[variant][field] for data in datas])
    return np.hstack([offsets, shared]), y


def pooled_penalty(state: dict, width: int, sessions: int) -> np.ndarray:
    """e2f.coefficient_penalty, remapped onto the pooled column layout."""
    single = e2f.coefficient_penalty(state, width - sessions + 1)
    if single.size == 0:
        return np.zeros((0, width))
    out = np.zeros((single.shape[0], width))
    out[:, sessions:] = single[:, 1:]
    return out


def fit_pooled(datas: list[dict], state: dict, penalty: float | None,
               variant: str = "D2", field: str = "ybar") -> dict:
    design, y = pooled_system(datas, state, variant, field)
    sessions = len(datas)
    if penalty is None:
        beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        normal_inverse = np.linalg.pinv(design.T @ design)
        edf = float(design.shape[1])
        solved = design
    else:
        roughness = pooled_penalty(state, design.shape[1], sessions)
        solved = np.vstack([design, math.sqrt(penalty) * roughness])
        beta, _, _, _ = np.linalg.lstsq(
            solved, np.r_[y, np.zeros(roughness.shape[0])], rcond=None)
        normal_inverse = np.linalg.inv(design.T @ design
                                       + penalty * (roughness.T @ roughness))
        edf = float(np.trace(design @ normal_inverse @ design.T))
    residual = y - design @ beta
    variance = float(residual @ residual / max(y.size - edf, 1.0))
    scaled = solved / np.maximum(np.linalg.norm(solved, axis=0), 1e-12)
    return {"beta": beta, "cov": variance * normal_inverse, "residual": residual,
            "sigma": math.sqrt(variance), "rms": rms(residual),
            "bic": float(y.size * np.log(np.mean(residual ** 2))
                         + edf * np.log(y.size)),
            "p": design.shape[1], "edf": edf, "lambda": penalty,
            "sessions": sessions, "state": state, "design": design, "y": y,
            "condition_scaled": float(np.linalg.cond(scaled))}


def as_single_session(pooled: dict, index: int = 0) -> dict:
    """Pooled beta and covariance in single-session shape, for reuse.

    e2f.predict, components and harmonic_profile all expect the design that
    prepare_design builds, whose first column is one offset. Taking that one
    session's offset column out of the pooled fit gives exactly that object,
    so the pooled model can be pushed through the same validation code as a
    single-session one with nothing reimplemented.
    """
    sessions = pooled["sessions"]
    keep = np.r_[index, np.arange(sessions, pooled["p"])]
    return {"beta": pooled["beta"][keep],
            "cov": pooled["cov"][np.ix_(keep, keep)],
            "state": pooled["state"], "p": keep.size,
            "sigma": pooled["sigma"], "rms": pooled["rms"],
            "edf": pooled["edf"], "lambda": pooled["lambda"]}


def make_state(data: dict, harmonics: int, slow_knots: int,
               coefficient_knots: int | None, domain: tuple[float, float]
               ) -> dict:
    _, state = e2f.prepare_design(data["x"], data["theta"], harmonics,
                                  slow_knots, coefficient_knots, domain)
    return state


def pooled_selection(datas: list[dict], domain: tuple[float, float]
                     ) -> tuple[list[dict], dict, dict]:
    """Section 7.1: re-run selection on the pooled data, two CV schemes."""
    x_all = np.concatenate([data["x"] for data in datas])
    folds = e2f.cv_assignments(x_all, domain)
    candidates = [("M_B", h, ks, None, None) for h in (2, 3, 4) for ks in range(9)]
    candidates += [("M_C", h, ks, e2f.PEN_KNOTS, lam) for h in (2, 3, 4)
                   for ks in range(9) for lam in e2f.PEN_LAMBDAS]
    rows = []
    for kind, h, ks, kc, lam in candidates:
        state = make_state(datas[0], h, ks, kc, domain)
        full = fit_pooled(datas, state, lam)
        interior, low, high = [], [], []
        for fold in folds:
            for held in range(10):
                train = fold != held
                if (~train).sum() < e2f.CV_MIN_TEST:
                    continue
                parts, masks = [], []
                start = 0
                for data in datas:
                    count = data["x"].size
                    masks.append(train[start:start + count])
                    start += count
                trimmed = [subset_data(data, mask)
                           for data, mask in zip(datas, masks)]
                trial = fit_pooled(trimmed, state, lam)
                single = as_single_session(trial)
                errors = []
                for data, mask in zip(datas, masks):
                    test = ~mask
                    if not test.any():
                        continue
                    prediction = e2f.predict(single, data["x"][test],
                                             data["theta"][test])
                    residual = data["D2"]["ybar"][test] - prediction
                    errors.append(residual - residual.mean()
                                  if residual.size > 1 else residual)
                value = rms(np.concatenate(errors))
                (low if held == 0 else high if held == 9 else interior).append(value)
        interior = np.asarray(interior)
        loso = leave_one_session_out(datas, state, lam)
        rows.append({"model": kind, "H": h, "K_s": ks,
                     "K_C": "stationary" if kc is None else kc,
                     "lambda": "" if lam is None else lam,
                     "parameters": full["p"],
                     "effective_parameters": round(full["edf"], 2),
                     "residual_rms_um": full["rms"], "BIC": full["bic"],
                     "CV_RMSE_um": float(interior.mean()),
                     "CV_SE_um": float(interior.std(ddof=1) / np.sqrt(interior.size)),
                     "CV_folds_scored": int(interior.size),
                     "end_low_RMSE_um": float(np.mean(low)),
                     "end_high_RMSE_um": float(np.mean(high)),
                     "LOSO_RMSE_um": loso["mean"],
                     "LOSO_worst_um": loso["worst"]})

    def one_se(kind: str) -> dict:
        subset = [row for row in rows if row["model"] == kind]
        best = min(subset, key=lambda row: row["CV_RMSE_um"])
        eligible = [row for row in subset
                    if row["CV_RMSE_um"] <= best["CV_RMSE_um"] + best["CV_SE_um"]]
        stable = [row for row in eligible
                  if max(row["end_low_RMSE_um"], row["end_high_RMSE_um"])
                  <= e2f.CV_END_FACTOR * row["CV_RMSE_um"]]
        return min(stable or eligible,
                   key=lambda row: (row["effective_parameters"], row["CV_RMSE_um"]))

    return rows, one_se("M_B"), one_se("M_C")


def subset_data(data: dict, mask: np.ndarray) -> dict:
    """The same preprocessing dict, restricted to the kept targets.

    Per-target arrays are masked; session-level entries (the drift anchors,
    the label) are carried through untouched, so a masked dict is still a
    valid input to drift_sensitivity and the reversal fits.
    """
    count = data["x"].size
    out = {}
    for key, value in data.items():
        if key in ("D0", "D1", "D2"):
            out[key] = {name: array[mask] for name, array in value.items()}
        elif isinstance(value, np.ndarray) and value.size == count:
            out[key] = value[mask]
        else:
            out[key] = value
    return out


def leave_one_session_out(datas: list[dict], state: dict,
                          penalty: float | None) -> dict:
    """Fit on one day, predict the other: the strictest test available."""
    scores = []
    for held in range(len(datas)):
        train = [data for j, data in enumerate(datas) if j != held]
        model = as_single_session(fit_pooled(train, state, penalty))
        residual = cross_residual(model, datas[held])
        scores.append(rms(residual))
    return {"mean": float(np.mean(scores)), "worst": float(np.max(scores)),
            "scores": scores}


# ------------------------------------------------- section 7.2 empirical Bayes

def prior_blocks(state: dict, width: int, sessions: int) -> dict:
    """Column groups the prior treats differently, in pooled layout."""
    groups = state["groups"]
    shift = sessions - 1
    blocks = {"slow": None, "harmonic": {}}
    if "slow" in groups:
        group = groups["slow"]
        blocks["slow"] = np.arange(group.start + shift, group.stop + shift)
    for k in range(1, state["H"] + 1):
        key = f"C{k}"
        if key not in groups:
            continue
        parts = groups[key]
        blocks["harmonic"][k] = [
            np.arange(part.start + shift, part.stop + shift) for part in parts]
    return blocks


def difference_matrix(count: int, order: int = e2f.PEN_ORDER) -> np.ndarray:
    out = np.eye(count)
    for _ in range(order):
        out = np.diff(out, axis=0)
    return out


def prior_precision(state: dict, width: int, sessions: int, log_tau2: float,
                    gamma: float, log_lambda_slow: float,
                    log_lambda_coefficient: float) -> np.ndarray:
    """Gaussian prior: harmonics shrink with order, curves get P-spline roughness.

    The offsets and the gain are left flat. They are the datum and the lead
    scale - the two things the experiment is there to measure - and shrinking
    them towards zero would be an assumption, not a regularisation.
    """
    precision = np.zeros((width, width))
    blocks = prior_blocks(state, width, sessions)
    tau2 = math.exp(log_tau2)
    if blocks["slow"] is not None and blocks["slow"].size > e2f.PEN_ORDER:
        index = blocks["slow"]
        d = difference_matrix(index.size)
        precision[np.ix_(index, index)] += math.exp(log_lambda_slow) * (d.T @ d)
    for k, parts in blocks["harmonic"].items():
        scale = (1.0 / tau2) * float(k) ** (2.0 * gamma)
        for index in parts:
            precision[np.ix_(index, index)] += scale * np.eye(index.size)
            if index.size > e2f.PEN_ORDER:
                d = difference_matrix(index.size)
                precision[np.ix_(index, index)] += (
                    math.exp(log_lambda_coefficient) * (d.T @ d))
    return precision


def marginal_likelihood(design: np.ndarray, y: np.ndarray,
                        precision: np.ndarray) -> tuple[float, np.ndarray, float]:
    """Profiled log marginal likelihood of the Gaussian linear model.

    beta ~ N(0, sigma^2 P^+) on the penalised block and flat elsewhere, so the
    flat directions cancel between log|P| and log|XtX + P| and only the
    penalised block contributes a pseudo-determinant.
    """
    gram = design.T @ design + precision
    try:
        chol = np.linalg.cholesky(gram)
    except np.linalg.LinAlgError:
        return -np.inf, np.zeros(design.shape[1]), np.inf
    rhs = design.T @ y
    beta = np.linalg.solve(gram, rhs)
    residual = y - design @ beta
    penalised = float(residual @ residual + beta @ precision @ beta)
    n = y.size
    sigma2 = penalised / n
    log_det_gram = 2.0 * float(np.sum(np.log(np.diag(chol))))
    eigenvalues = np.linalg.eigvalsh(precision)
    positive = eigenvalues[eigenvalues > eigenvalues.max() * 1e-10]
    log_det_prior = float(np.sum(np.log(positive)))
    value = 0.5 * (log_det_prior - log_det_gram) - 0.5 * n * (1.0 + math.log(sigma2))
    return value, beta, sigma2


def empirical_bayes(datas: list[dict], state: dict) -> dict:
    """Section 7.2: tau^2, gamma and the two smoothnesses by marginal likelihood."""
    design, y = pooled_system(datas, state)
    sessions = len(datas)
    width = design.shape[1]

    def objective(theta):
        precision = prior_precision(state, width, sessions, *theta)
        value, _, _ = marginal_likelihood(design, y, precision)
        return -value

    best, start = None, None
    for log_tau2 in (-2.0, 0.0, 2.0):
        for gamma in (0.0, 1.0, 2.0):
            theta = np.array([log_tau2, gamma, 0.0, 0.0])
            value = objective(theta)
            if best is None or value < best:
                best, start = value, theta
    result = minimize(objective, start, method="Nelder-Mead",
                      options={"maxiter": 4000, "xatol": 1e-3, "fatol": 1e-4})
    precision = prior_precision(state, width, sessions, *result.x)
    value, beta, sigma2 = marginal_likelihood(design, y, precision)
    gram = design.T @ design + precision
    inverse = np.linalg.inv(gram)
    residual = y - design @ beta
    posterior = {
        "beta": beta, "cov": sigma2 * inverse, "residual": residual,
        "sigma": math.sqrt(sigma2), "rms": rms(residual),
        "p": width, "edf": float(np.trace(design @ inverse @ design.T)),
        "lambda": None, "sessions": sessions, "state": state,
        "design": design, "y": y,
        "condition_scaled": float(np.linalg.cond(
            design / np.maximum(np.linalg.norm(design, axis=0), 1e-12))),
    }
    hyper = {"tau2": math.exp(result.x[0]), "gamma": float(result.x[1]),
             "lambda_slow": math.exp(result.x[2]),
             "lambda_coefficient": math.exp(result.x[3]),
             "sigma_um": math.sqrt(sigma2),
             "log_marginal_likelihood": value,
             "effective_parameters": posterior["edf"],
             "residual_rms_um": posterior["rms"],
             "converged": bool(result.success), "iterations": int(result.nit)}
    return {"model": posterior, "hyperparameters": hyper}


# ---------------------------------------------------------------- section 8

def split_direction_leakage(datas: list[dict], state: dict,
                            penalty: float | None) -> dict:
    """Section 8.3: is the up map the down map plus a constant?

    H0 gives the two directions one shared geometry and their own offsets.
    H1 adds a direction-specific slow spline and stationary C1/C2 on top. If
    friction is translation-invariant the added terms buy nothing, and the F
    test says so.
    """
    designs, ys, offsets = [], [], []
    rows = 0
    pieces = []
    for session_index, data in enumerate(datas):
        for direction, field in ((+1, "e_up"), (-1, "e_down")):
            design = session_design(data, state)
            pieces.append((design, data, direction, session_index))
            ys.append(data["D2"][field])
            rows += data["x"].size
    groups = len(pieces)
    offset = np.zeros((rows, groups))
    start = 0
    shared, flags, xs, thetas = [], [], [], []
    for j, (design, data, direction, _) in enumerate(pieces):
        count = data["x"].size
        offset[start:start + count, j] = 1.0
        start += count
        shared.append(design[:, 1:])
        flags.append(np.full(count, direction, float))
        xs.append(data["x"]); thetas.append(data["theta"])
    shared = np.vstack(shared)
    flag = np.concatenate(flags)
    x = np.concatenate(xs); theta = np.concatenate(thetas)
    y = np.concatenate(ys)
    design0 = np.hstack([offset, shared])
    # The difference basis is deliberately coarse: a slow spline plus
    # stationary C1/C2. The test is whether the difference is constant, and a
    # rich basis would spend its power on fitting noise.
    extra, _ = e2f.prepare_design(x, theta, 2, 2, None, state["domain"])
    design1 = np.hstack([design0, flag[:, None] * extra[:, 1:]])
    out = {}
    for name, design in (("H0_constant_difference", design0),
                         ("H1_structured_difference", design1)):
        beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ beta
        out[name] = {"rss": float(residual @ residual), "p": design.shape[1],
                     "rms_um": rms(residual)}
    n = y.size
    p0, p1 = out["H0_constant_difference"]["p"], out["H1_structured_difference"]["p"]
    rss0, rss1 = (out["H0_constant_difference"]["rss"],
                  out["H1_structured_difference"]["rss"])
    statistic = ((rss0 - rss1) / (p1 - p0)) / (rss1 / (n - p1))
    out["F"] = {"statistic": float(statistic), "df1": p1 - p0, "df2": n - p1,
                "p_value": float(f_dist.sf(statistic, p1 - p0, n - p1))}
    # Size of the leakage in microns, whatever the test says about it.
    beta1, _, _, _ = np.linalg.lstsq(design1, y, rcond=None)
    delta = extra[:, 1:] @ beta1[design0.shape[1]:]
    order = np.argsort(x)
    out["difference_profile"] = {
        "x_mm": x[order][::8].tolist(),
        "difference_um": (2.0 * delta)[order][::8].tolist(),
        "max_abs_um": float(np.max(np.abs(2.0 * delta))),
        "rms_um": rms(2.0 * delta),
        "note": "s_up(x) - s_down(x), twice the fitted direction term"}
    return out


def reversal_ripple(data: dict, variant: str = "D2") -> dict:
    """Section 8.4: B = beta0 + Re[D1(x) exp(i theta)], coarse spline in x."""
    x, theta = data["x"], data["theta"]
    values = data[variant]["B"]
    basis = e2f.raw_bspline(x, 2, (float(x.min()), float(x.max())))
    design = np.hstack([np.ones((x.size, 1)),
                        basis * np.cos(theta)[:, None],
                        basis * np.sin(theta)[:, None]])
    beta, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ beta
    width = basis.shape[1]
    grid = np.arange(math.ceil(x.min()), math.floor(x.max()) + .1, 1.0)
    grid_basis = e2f.raw_bspline(grid, 2, (float(x.min()), float(x.max())))
    a = grid_basis @ beta[1:1 + width]
    b = grid_basis @ beta[1 + width:]
    amplitude = np.hypot(a, b)
    return {"beta0_um": float(beta[0]), "rms_um": rms(residual),
            "x_mm": grid.tolist(), "amplitude_um": amplitude.tolist(),
            "amplitude_min_um": float(amplitude.min()),
            "amplitude_max_um": float(amplitude.max()),
            "amplitude_near_motor_um": float(amplitude[:5].mean()),
            "amplitude_mid_um": float(amplitude[amplitude.size // 2 - 2:
                                                amplitude.size // 2 + 3].mean()),
            "amplitude_near_support_um": float(amplitude[-5:].mean())}


# ---------------------------------------------------------------- section 9

def variant_robustness(data: dict, row: dict, domain: tuple[float, float]
                       ) -> dict:
    """Every headline number under D0, D1 and D2.

    Section 3 says a feature is real only if it survives all three variants.
    With E2 the three agreed to a nanometre and the rule was a formality; with
    this session it decides which results may be quoted at all, so the spread
    is carried into the summary as a number rather than a reassurance.
    """
    out = {}
    for variant in ("D0", "D1", "D2"):
        model = fit_with_settings(data, row, domain, variant=variant,
                                 robust=False)
        reversal = e2f.reversal_fit(data, variant)
        grid = np.arange(WELL_SUPPORTED[0], WELL_SUPPORTED[1] + .01, GRID_STEP)
        slow = e2f.components(model, grid)["slow"]
        profile = {k: e2f.harmonic_profile(model, grid, k) for k in (1, 2)}
        out[variant] = {
            "gain_ppm": float(model["beta"][1] * 1000.0),
            "residual_rms_um": model["rms"],
            "slow_range_um": float(np.ptp(slow)),
            "beta0_um": reversal["mean_um"],
            "A1_mean_um": float(profile[1]["amplitude"].mean()),
            "A2_mean_um": float(profile[2]["amplitude"].mean()),
        }
    keys = ("gain_ppm", "slow_range_um", "beta0_um", "A1_mean_um", "A2_mean_um")
    out["spread"] = {key: float(np.ptp([out[v][key] for v in ("D0", "D1", "D2")]))
                     for key in keys}
    # The spread is not an error bar, and treating it as one would be wrong in
    # both directions. D0 applies no correction and D1, given the mirrored
    # design, provably removes only a constant from ybar - so with a drift
    # this large neither is a candidate answer, and the D0-to-D2 distance
    # measures how much contamination D2 took out rather than how uncertain
    # the result is. What is genuinely uncertain is whether a straight line
    # between three references describes the real drift, and three references
    # cannot answer that. Quantities whose spread is negligible are safe
    # whichever view one takes, which is the useful half of the rule.
    out["interpretation"] = {
        "admissible_variant": "D2",
        "D0_D1_inadmissible_because": "no correction, and a linear-in-time "
                                      "correction cancels against the mirrored "
                                      "pass design",
        "spread_meaning": "contamination removed by D2, not an uncertainty",
        "unbounded_residual_risk": "the adequacy of a piecewise-linear drift "
                                   "model through three references",
        "drift_insensitive": sorted(key for key in keys
                                    if out["spread"][key] < 0.01),
    }
    return out


def slow_difference_is_drift(data_b: dict, model_a: dict, model_b: dict,
                             grid: np.ndarray, gain_difference_ppm: float
                             ) -> dict:
    """How much of the gain and s(x) difference could be drift mismodelling?

    The drift correction D2 applied to this session is known exactly - it is
    ybar_D1 minus ybar_D2, built from the reference visits alone - and it has
    to be decomposed before it can be compared with anything. s(x) is
    orthogonalised against {1, x} by construction, so every linear part of the
    correction lands in the gain and none of it in s(x). Correlating the raw
    correction against s(x) therefore compares quantities living in different
    subspaces and returns nothing whatever the truth is; the comparison has to
    be made within each subspace separately.

    This is a comparison of scales, not a correlation test. The correction is
    what a piecewise-linear drift model took out; what cannot be known from
    three references is how wrong that model was. So the useful statement is
    the fraction of the correction that the between-session difference amounts
    to: if the difference is a small part of a correction this large, a modest
    error in the drift model accounts for it and the difference is not
    evidence about the screw.
    """
    correction = data_b["D1"]["ybar"] - data_b["D2"]["ybar"]
    x = data_b["x"]
    q = np.column_stack((np.ones(x.size), x - x.mean()))
    coefficients, *_ = np.linalg.lstsq(q, correction, rcond=None)
    linear_ppm = float(coefficients[1] * 1000.0)
    curvature = correction - q @ coefficients
    profile = running_mean(x, curvature, grid, 4.0)
    slow_difference = (e2f.components(model_b, grid)["slow"]
                       - e2f.components(model_a, grid)["slow"])
    ok = ~np.isnan(profile)
    rho = float(np.corrcoef(slow_difference[ok] - slow_difference[ok].mean(),
                            profile[ok])[0, 1])
    # With a piecewise-linear drift and two mirrored monotone passes the
    # correction to ybar is linear in x exactly: on the up leg the time runs
    # with x, on the down leg against it, so the two segment slopes combine
    # into a single slope and the curvature cancels. The measured curvature
    # range confirms it at the nanometre. The useful consequence is that the
    # drift this session suffered reaches the gain and reaches nothing else -
    # so the gain difference between the days is explicable as drift
    # mismodelling, and the s(x) difference is not.
    curvature_range = float(np.ptp(profile[ok])) if ok.any() else float("nan")
    degenerate = not (curvature_range > 0.01)
    return {
        "correction_linear_ppm": linear_ppm,
        "correction_curvature_range_um": curvature_range,
        "correction_is_linear_in_x": bool(degenerate),
        "gain_difference_ppm": gain_difference_ppm,
        "gain_difference_over_correction": abs(gain_difference_ppm
                                               / linear_ppm) if linear_ppm else float("nan"),
        "slow_difference_range_um": float(np.ptp(slow_difference)),
        "slow_difference_over_correction": (
            float("nan") if degenerate
            else float(np.ptp(slow_difference) / curvature_range)),
        "curvature_rho": float("nan") if degenerate else rho,
        "gain_verdict": "the correction applied is far larger than the "
                        "between-day difference, so the gain difference is "
                        "within reach of an error in the drift model and is "
                        "not evidence about the screw",
        "slow_verdict": "the correction D2 applied moves s(x) by nothing, so "
                        "the s(x) difference is either a real day-to-day "
                        "change or drift curvature *within* a segment, which "
                        "three reference visits cannot see. Not separable "
                        "here; more reference visits would separate it.",
    }


def iso_validation(model: dict, beta0_um: float, path: Path) -> tuple[list, dict]:
    e2f.ISO_RESULTS = path
    return e2f.validate_iso(model, {"mean_um": beta0_um})


def iso_bootstrap(model: dict, beta0_um: float, beta0_sigma_um: float,
                  path: Path, draws: int = N_BOOTSTRAP) -> dict:
    """E5 prediction with uncertainty: posterior draws of the model and beta0."""
    rng = np.random.default_rng(RNG_SEED)
    mean = np.asarray(model["beta"])
    covariance = np.asarray(model["cov"])
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, vectors = np.linalg.eigh(covariance)
    root = vectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    collected = {"E": [], "M": [], "A": []}
    for _ in range(draws):
        draw = dict(model)
        draw["beta"] = mean + root @ rng.standard_normal(mean.size)
        _, metrics = iso_validation(draw, beta0_um + beta0_sigma_um
                                    * rng.standard_normal(), path)
        for key in collected:
            collected[key].append(metrics["G_plus_reversal"][key])
    return {key: {"mean_um": float(np.mean(values)),
                  "sd_um": float(np.std(values, ddof=1)),
                  "p5_um": float(np.percentile(values, 5)),
                  "p95_um": float(np.percentile(values, 95))}
            for key, values in collected.items()}


# ---------------------------------------------------------------- output

def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return float(value)


def main() -> int:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    published = json.loads(E2_SUMMARY.read_text(encoding="utf-8"))

    e2f.LOG = E2_DIR / "data" / "e2_terminal_20261027_151332.log"
    s2 = reduce_session(e2_run, "E2")
    e2f.LOG = e3_run.LOG
    s3 = reduce_session(e3_run, "E3")
    d2, d3 = s2["data"], s3["data"]

    # The basis is E2's realised span and reference position, used unchanged
    # for E3 and for the pooled fit. With identical targets this is automatic,
    # and asserting it is cheaper than trusting it.
    domain = basis_domain(d2)
    assert basis_domain(d3) == domain, "the two sessions span different travel"
    assert np.array_equal(d2["usteps"], d3["usteps"])
    assert np.all(d2["usteps"] % 64 == 0), "a target escaped the 64-ustep grid"
    state_reference = make_state(d2, 4, 2, e2f.PEN_KNOTS, domain)
    xbar = state_reference["xbar"]

    iso_usteps = iso_target_usteps(config)
    validity = {"E2": validity_rows(s2, config, domain, iso_usteps, None),
                "E3": validity_rows(s3, config, domain, iso_usteps, d2)}

    # One E3 window in 807 fails the std and slope gates: an UP stop whose
    # hold was short enough that the averaging window overlapped settling. It
    # is the session's largest single residual. A flagged window is excluded
    # rather than down-weighted, and the same mask is applied to E2 - where it
    # removes nothing - so both sessions still go through one code path.
    flags = {"E2": window_quality(s2, config), "E3": window_quality(s3, config)}
    keep = flags["E2"]["keep"] & flags["E3"]["keep"]
    excluded = [{"session": label,
                 "target_mm": float(d3["x"][i]),
                 "target_usteps": int(d3["usteps"][i]),
                 "std_up_um": float(flags[label]["std_up_um"][i]),
                 "std_down_um": float(flags[label]["std_down_um"][i]),
                 "slope_effect_up_um": float(flags[label]["slope_effect_up_um"][i]),
                 "slope_effect_down_um": float(flags[label]["slope_effect_down_um"][i]),
                 "samples_up": int(flags[label]["samples_up"][i]),
                 "samples_down": int(flags[label]["samples_down"][i])}
                for label in ("E2", "E3")
                for i in np.flatnonzero(flags[label]["flagged"])]
    write_csv(ANALYSIS / "e3_excluded_windows.csv", excluded)
    full2, full3 = d2, d3
    d2, d3 = subset_data(d2, keep), subset_data(d3, keep)
    for name, rows in validity.items():
        rows.append({"check": "record-window exclusions",
                     "value": int(np.count_nonzero(flags[name]["flagged"])),
                     "criterion": "0 preferred, logged if not",
                     "pass": bool(not flags[name]["flagged"].any())})
        rows.append({"check": "targets carried into the fit",
                     "value": int(keep.sum()), "criterion": 400,
                     "pass": bool(keep.all())})
    write_csv(ANALYSIS / "e3_t1_validity_table.csv",
              [{"session": name, **row} for name, rows in validity.items()
               for row in rows])

    drift = {"E2": drift_sensitivity(d2), "E3": drift_sensitivity(d3)}

    # Section 5.1 - E3 under the published E2 settings, and E2 refitted here
    # so the two sides of every comparison come from this run of this script.
    row_b = published["model_selection"]["M_B_one_se"]
    row_c = published["model_selection"]["M_C_one_se"]
    models = {}
    for label, data in (("E2", d2), ("E3", d3)):
        models[label] = {
            "M_B": fit_with_settings(data, row_b, domain),
            "M_C": fit_with_settings(data, row_c, domain)}

    # Section 5.2 - cross prediction, both directions
    cross_rows = cross_prediction_table(
        d3, models["E3"]["M_B"], models["E3"]["M_C"],
        models["E2"]["M_B"], models["E2"]["M_C"], "E2 -> E3", domain)
    cross_rows += cross_prediction_table(
        d2, models["E2"]["M_B"], models["E2"]["M_C"],
        models["E3"]["M_B"], models["E3"]["M_C"], "E3 -> E2", domain)
    write_csv(ANALYSIS / "e3_t2_cross_prediction.csv", cross_rows)

    # Section 5.3 / 5.4
    grid = np.arange(WELL_SUPPORTED[0], WELL_SUPPORTED[1] + .01, GRID_STEP)
    coefficient_rows, coefficient_summary = coefficient_comparison(
        models["E2"]["M_C"], models["E3"]["M_C"], grid)
    write_csv(ANALYSIS / "e3_coefficient_comparison.csv", coefficient_rows)
    homing = homing_delta(models["E2"]["M_C"], models["E3"]["M_C"], grid)
    homing["mscnt"] = mscnt_prediction(s2["serial"], s3["serial"],
                                       config["mres"])
    # Magnitudes, not signed values: the serial log fixes neither the
    # direction homing moved the microstep origin nor how MSCNT's count
    # direction relates to the harmonic phase convention, so a sign agreement
    # here would be luck rather than evidence. The magnitudes are a real test
    # and they are what gets compared.
    for k in (1, 2):
        residual = (abs(homing[f"delta{k}_deg"])
                    - abs(homing["mscnt"]["screw_deg"]))
        homing[f"delta{k}_vs_mscnt_deg"] = residual
        homing[f"delta{k}_matches_mscnt"] = bool(
            abs(residual) <= 2.0 * homing[f"delta{k}_sigma_deg"])
    homing["sign_convention"] = (
        "compared on magnitude; MSCNT's direction relative to the harmonic "
        "phase convention is not fixed by the log")

    # Section 3's rule, applied before anything is quoted.
    robustness = {"E2": variant_robustness(d2, row_c, domain),
                  "E3": variant_robustness(d3, row_c, domain)}
    coefficient_summary["slow_vs_drift"] = slow_difference_is_drift(
        d3, models["E2"]["M_C"], models["E3"]["M_C"], grid,
        coefficient_summary["gain"]["difference_ppm"])

    sliding = {label: e2f.sliding_profile(data)
               for label, data in (("E2", d2), ("E3", d3))}
    reversal_sliding = {label: e2f.sliding_profile(data, field="B")
                        for label, data in (("E2", d2), ("E3", d3))}

    # Section 7 - pooled identification
    datas = [d2, d3]
    pooled_rows, pooled_b_row, pooled_c_row = pooled_selection(datas, domain)
    write_csv(ANALYSIS / "e3_t4_pooled_model_selection.csv", pooled_rows)
    state_b = make_state(d2, int(pooled_b_row["H"]), int(pooled_b_row["K_s"]),
                         None, domain)
    kc_c = (None if pooled_c_row["K_C"] == "stationary"
            else int(pooled_c_row["K_C"]))
    state_c = make_state(d2, int(pooled_c_row["H"]), int(pooled_c_row["K_s"]),
                         kc_c, domain)
    lam_c = (float(pooled_c_row["lambda"])
             if pooled_c_row["lambda"] != "" else None)
    pooled_b = fit_pooled(datas, state_b, None)
    pooled_c = fit_pooled(datas, state_c, lam_c)
    choose_c = (pooled_b_row["CV_RMSE_um"] - pooled_c_row["CV_RMSE_um"]
                > max(pooled_b_row["CV_SE_um"], pooled_c_row["CV_SE_um"]))
    pooled_selected = pooled_c if choose_c else pooled_b
    selected_row = pooled_c_row if choose_c else pooled_b_row
    single_b = as_single_session(pooled_b)
    single_c = as_single_session(pooled_c)
    single_selected = single_c if choose_c else single_b

    bayes = empirical_bayes(datas, state_c)
    single_bayes = as_single_session(bayes["model"])

    # Section 6 - residual structure, under own models and under the pooled one
    counts = [data["x"].size for data in datas]
    pooled_residual = pooled_selected["residual"]
    pooled2 = pooled_residual[:counts[0]]
    pooled3 = pooled_residual[counts[0]:]
    structure = residual_structure(
        d2, d3,
        models["E2"]["M_C"]["residual"] - models["E2"]["M_C"]["residual"].mean(),
        models["E3"]["M_C"]["residual"] - models["E3"]["M_C"]["residual"].mean(),
        pooled2 - pooled2.mean(), pooled3 - pooled3.mean(),
        pooled_selected["edf"])

    # Section 8 - reversal and leakage
    reversal = {label: {name: e2f.reversal_fit(data, name)
                        for name in ("D0", "D1", "D2")}
                for label, data in (("E2", d2), ("E3", d3))}
    reversal_summary = {
        label: {name: {key: result[key] for key in
                       ("mean_um", "mean_sigma_um", "slope_um_per_100mm",
                        "slope_sigma_um_per_100mm")}
                for name, result in variants.items()}
        for label, variants in reversal.items()}
    # Section 3's rule decides whether E3 may contribute a reversal at all.
    # B differences the two directions, so it carries the whole drift over the
    # up-to-down gap; this session's three variants put beta0 between -1.1 and
    # -3.8 um. That is not a measurement of the reversal, it is a measurement
    # of the drift, and it is reported as such rather than averaged in.
    e3_beta0_spread = robustness["E3"]["spread"]["beta0_um"]
    e3_beta0_usable = bool(e3_beta0_spread <= 0.5)
    beta0_all = {"E2": reversal["E2"]["D2"]["mean_um"],
                 "E3": reversal["E3"]["D2"]["mean_um"],
                 "E4": -1.66, "E0": -1.60}
    applied = {name: value for name, value in beta0_all.items()
               if name != "E3" or e3_beta0_usable}
    reversal_summary["session_spread"] = {
        "values_um": beta0_all,
        "E3_variant_spread_um": e3_beta0_spread,
        "E3_usable": e3_beta0_usable,
        "applied_sessions": sorted(applied),
        "mean_um": float(np.mean(list(applied.values()))),
        "range_um": float(np.ptp(list(applied.values()))),
        "sd_um": float(np.std(list(applied.values()), ddof=1)),
        "range_including_E3_um": float(np.ptp(list(beta0_all.values())))}
    reversal_summary["ripple"] = {label: reversal_ripple(data)
                                  for label, data in (("E2", d2), ("E3", d3))}
    # Section 8.3, run per session as well as pooled. s_up(x) - s_down(x) is
    # the smooth part of B(x), and B is the drift-sensitive combination, so a
    # pooled test that includes this session measures its drift and calls it
    # friction. E2 alone is the admissible answer; the other two are reported
    # so the contamination is visible rather than hidden.
    leakage = {name: split_direction_leakage(subset, state_c, lam_c)
               for name, subset in (("E2", [d2]), ("E3", [d3]),
                                    ("pooled", datas))}
    leakage["admissible"] = "E2"
    leakage["note"] = ("E3 and the pooled test are contaminated: the "
                       "direction difference is the smooth part of B, which "
                       "carries the session drift over the up-to-down gap")

    # Section 9 - validation
    windows1, windows4 = e2f.external_windows()
    v2 = e2f.validate_windows(single_b, single_c, windows1, windows4)
    write_csv(ANALYSIS / "e3_t5_v2_window_validation.csv", v2)
    beta0 = reversal_summary["session_spread"]["mean_um"]
    beta0_sigma = reversal_summary["session_spread"]["sd_um"]
    v3_tables, iso_metrics_all = {}, {}
    for name, path in (("v3", V3_ISO), ("v4", V4_ISO)):
        rows, metrics = iso_validation(single_selected, beta0, path)
        v3_tables[name] = rows
        iso_metrics_all[name] = metrics
        write_csv(ANALYSIS / f"e3_t5_v3_iso_{name}.csv", rows)
    e5 = {name: iso_bootstrap(single_selected, beta0, beta0_sigma, path)
          for name, path in (("v3", V3_ISO), ("v4", V4_ISO))}

    target_rows = []
    for i, x in enumerate(d3["x"]):
        row = {"target_mm": float(x), "target_usteps": int(d3["usteps"][i]),
               "screw_angle_deg": float(np.degrees(d3["theta"][i])),
               "phase_bin": int(d3["phase_bin"][i]),
               "up_time_s": float(d3["up_time_s"][i]),
               "down_time_s": float(d3["down_time_s"][i])}
        for variant in ("D0", "D1", "D2"):
            row |= {f"e_up_{variant}_um": float(d3[variant]["e_up"][i]),
                    f"e_down_{variant}_um": float(d3[variant]["e_down"][i]),
                    f"ybar_{variant}_um": float(d3[variant]["ybar"][i]),
                    f"B_{variant}_um": float(d3[variant]["B"][i])}
        row |= {"ybar_E2_D2_um": float(d2["D2"]["ybar"][i]),
                "residual_own_E2_um": float(models["E2"]["M_C"]["residual"][i]),
                "residual_own_E3_um": float(models["E3"]["M_C"]["residual"][i]),
                "residual_pooled_E2_um": float(pooled2[i] - pooled2.mean()),
                "residual_pooled_E3_um": float(pooled3[i] - pooled3.mean())}
        target_rows.append(row)
    write_csv(ANALYSIS / "e3_preprocessed_targets.csv", target_rows)

    overall = {row["direction"]: row for row in cross_rows
               if row["block"] == "overall"}
    summary = {
        "experiment": "E3_repeat_global_map",
        "session": {
            "label": e3_run.SESSION_LABEL,
            "capture": e3_run.SCOPE.name,
            "log": e3_run.LOG.name,
            "scope_samples": int(s3["run"].t.size),
            "recording_min": float(s3["run"].t[-1] / 60.0),
            "firmware": "e2_global_map, config 96579c689d4d, unmodified",
            "seed": 5315,
            "targets_paired": int(full3["x"].size),
            "targets_fitted": int(d3["x"].size),
        },
        "design_departure": {
            "specified": "E3 = same design, different seed, new positions",
            "delivered": "same seed 5315, identical 400 targets, different day",
            "evidence": "the 400 # TARGET lines of the two logs are identical; "
                        "e2_config.json seed_start was not changed",
            "Q1_new_positions": "NOT ANSWERABLE from this session; the map is "
                                "tested only where it was fitted",
            "Q1_new_day": "answerable, and reported below",
            "gain_for_Q2_Q3_Q5": "residuals pair target by target, so the "
                                 "between-day comparison needs no smoothing "
                                 "and no interpolation",
        },
        "basis": {"domain_mm": list(domain), "xbar_mm": float(xbar),
                  "source": "E2 realised span and reference position, "
                            "used unchanged for E3 and the pooled fit"},
        "validity": validity,
        "excluded_windows": excluded,
        "drift": drift,
        "drift_variant_robustness": robustness,
        "single_session_models": {
            label: {kind: {"H": int(row["H"]), "K_s": int(row["K_s"]),
                           "K_C": row["K_C"], "lambda": row.get("lambda", ""),
                           "residual_rms_um": models[label][kind]["rms"],
                           "sigma_um": models[label][kind]["sigma"],
                           "effective_parameters": models[label][kind]["edf"],
                           "gain_ppm": float(models[label][kind]["beta"][1] * 1000),
                           "gain_sigma_ppm": float(np.sqrt(
                               models[label][kind]["cov"][1, 1]) * 1000)}
                    for kind, row in (("M_B", row_b), ("M_C", row_c))}
            for label in ("E2", "E3")},
        "cross_prediction_overall": overall,
        "coefficient_comparison": coefficient_summary,
        "homing_check": homing,
        "residual_structure": structure,
        "pooled": {
            "selection": {"M_B_one_se": pooled_b_row, "M_C_one_se": pooled_c_row,
                          "selected": "M_C" if choose_c else "M_B"},
            "selected_model": {
                "H": int(selected_row["H"]), "K_s": int(selected_row["K_s"]),
                "K_C": selected_row["K_C"], "lambda": selected_row["lambda"],
                "parameters": pooled_selected["p"],
                "effective_parameters": pooled_selected["edf"],
                "residual_rms_um": pooled_selected["rms"],
                "sigma_um": pooled_selected["sigma"],
                "condition_number_scaled": pooled_selected["condition_scaled"],
                "session_offsets_um": pooled_selected["beta"][:2].tolist(),
                "gain_ppm": float(single_selected["beta"][1] * 1000),
                "gain_sigma_ppm": float(np.sqrt(single_selected["cov"][1, 1]) * 1000),
                "LOSO_RMSE_um": selected_row["LOSO_RMSE_um"]},
            "empirical_bayes": bayes["hyperparameters"],
        },
        "reversal": reversal_summary,
        "leakage": leakage,
        "validation": {
            "V2": v2,
            "V2_mean_rms_um": {
                kind: float(np.mean([row[f"rms_{kind}_um"] for row in v2]))
                for kind in ("raw", "M_B", "M_C", "local")},
            "V3": {name: metrics for name, metrics in iso_metrics_all.items()},
            "E5_prediction": e5,
        },
        "limitations": {
            "new_positions": "not tested; see design_departure",
            "E3_drift": f"the reference moved "
                        f"{drift['E3']['total_change_um']:+.2f} um over the "
                        f"session against a +-0.5 um criterion, and not "
                        f"linearly: see drift_variant_robustness for which "
                        f"results survive D0/D1/D2 and which do not",
            "E3_reversal": "beta0 is not identifiable in this session; B "
                           "carries the full drift over the up-to-down gap",
            "environment": "neither session recorded air or stage temperature, "
                           "so a gain difference cannot be attributed thermally "
                           "and the E3 drift cannot be confirmed as thermal",
            "model_D_gp": "section 7.3 structured GP not run; it is marked "
                          "optional and the P-spline posterior bands already "
                          "carry the smoothness uncertainty",
            "endpoint_locations": "inferred from homed travel",
        },
        "convention": published["convention"],
    }
    (ANALYSIS / "e3_findings_summary.json").write_text(
        json.dumps(summary, indent=2, default=jsonable) + "\n", encoding="utf-8")

    np.savez_compressed(
        ANALYSIS / "e3_models.npz",
        pooled_beta=pooled_selected["beta"], pooled_cov=pooled_selected["cov"],
        bayes_beta=bayes["model"]["beta"], bayes_cov=bayes["model"]["cov"],
        single_beta=single_selected["beta"], single_cov=single_selected["cov"],
        e2_beta=models["E2"]["M_C"]["beta"], e3_beta=models["E3"]["M_C"]["beta"],
        e2_cov=models["E2"]["M_C"]["cov"], e3_cov=models["E3"]["M_C"]["cov"],
        grid=grid, domain=np.asarray(domain), xbar=xbar,
        pooled_H=int(selected_row["H"]), pooled_Ks=int(selected_row["K_s"]),
        pooled_Kc=-1 if selected_row["K_C"] == "stationary" else int(selected_row["K_C"]),
        pooled_lambda=float(selected_row["lambda"] or "nan"),
        residual_e2=pooled2 - pooled2.mean(), residual_e3=pooled3 - pooled3.mean(),
        own_residual_e2=models["E2"]["M_C"]["residual"],
        own_residual_e3=models["E3"]["M_C"]["residual"],
    )
    for label in ("E2", "E3"):
        write_csv(ANALYSIS / f"e3_sliding_{label.lower()}_ybar.csv",
                  [{k: (v.real if isinstance(v, complex) else v)
                    for k, v in row.items() if not isinstance(v, complex)}
                   | {f"{k}_real_um": row[k].real for k in ("C1", "C2")}
                   | {f"{k}_imag_um": row[k].imag for k in ("C1", "C2")}
                   for row in sliding[label]])
        write_csv(ANALYSIS / f"e3_sliding_{label.lower()}_reversal.csv",
                  [{k: v for k, v in row.items() if not isinstance(v, complex)}
                   | {f"{k}_real_um": row[k].real for k in ("C1", "C2")}
                   | {f"{k}_imag_um": row[k].imag for k in ("C1", "C2")}
                   for row in reversal_sliding[label]])

    print(f"E2 {d2['x'].size} paired targets; E3 {d3['x'].size} paired targets")
    print(f"basis: {domain[0]:.3f}..{domain[1]:.3f} mm, xbar {xbar:.4f} mm")
    print(f"E3 own M_C rms {models['E3']['M_C']['rms']:.3f} um; "
          f"E2 own M_C rms {models['E2']['M_C']['rms']:.3f} um")
    for row in cross_rows:
        if row["block"] == "overall":
            print(f"{row['direction']}: raw {row['rms_raw_um']:.3f} -> "
                  f"M_B {row['rms_cross_M_B_um']:.3f} -> "
                  f"M_C {row['rms_cross_M_C_um']:.3f} "
                  f"(own {row['rms_own_M_C_um']:.3f}, floor {E4_FLOOR_UM}); "
                  f"gap closed {row['gap_closed_M_C']:.2f}")
    print(f"gain E2 {coefficient_summary['gain']['E2_ppm']:+.2f} vs E3 "
          f"{coefficient_summary['gain']['E3_ppm']:+.2f} ppm, "
          f"z {coefficient_summary['gain']['z']:+.2f}")
    for k in (1, 2):
        block = coefficient_summary[f"C{k}"]
        print(f"  C{k}: max |dA| {block['max_abs_amplitude_difference_um']:.3f} um, "
              f"max |dC| {block['max_complex_difference_um']:.3f} um, "
              f"rms |dC| {block['rms_complex_difference_um']:.3f} um")
    sd = coefficient_summary["slow_vs_drift"]
    print(f"  drift correction D2 applied: {sd['correction_linear_ppm']:+.1f} ppm "
          f"linear + {sd['correction_curvature_range_um']:.4f} um curvature "
          f"(linear in x: {sd['correction_is_linear_in_x']})")
    print(f"  gain difference is {sd['gain_difference_over_correction']:.1%} of "
          f"the correction applied -> explicable as drift mismodelling")
    print(f"  s(x) difference {sd['slow_difference_range_um']:.3f} um is not "
          f"moved by that correction -> real change or within-segment drift "
          f"curvature, not separable with three references")
    print(f"E3 drift: {drift['E3']['total_change_um']:+.2f} um over the run, "
          f"D0 vs D1 on ybar {drift['E3']['D0_vs_D1_max_ybar_um']:.2e} um, "
          f"D1-D2 curvature {drift['E3']['D1_vs_D2_max_ybar_um']:.3f} um")
    for label in ("E2", "E3"):
        spread = robustness[label]["spread"]
        print(f"  {label} variant spread: gain {spread['gain_ppm']:.2f} ppm, "
              f"s(x) range {spread['slow_range_um']:.3f} um, beta0 "
              f"{spread['beta0_um']:.3f} um, A1 {spread['A1_mean_um']:.3f} um")
    print(f"homing delta1 {homing['delta1_deg']:+.3f} deg, delta2 "
          f"{homing['delta2_deg']:+.3f} deg, consistent {homing['consistent']}; "
          f"MSCNT predicts {homing['mscnt']['screw_deg']:+.4f} deg "
          f"({homing['mscnt']['usteps']:+.1f} usteps, "
          f"{homing['mscnt']['travel_um']:+.2f} um)")
    free = structure["model_free"]
    paired = structure["paired_own_models"]
    print(f"day-to-day, model free: corr(ybar) {free['rho_ybar']:.3f}, "
          f"per-session difference {free['day_to_day_per_session_um']:.3f} um")
    print(f"paired residual rho {paired['rho']:+.3f} (n={paired['n']}): "
          f"{paired['repeatable_rms_um']:.3f} um repeatable + "
          f"{paired['independent_rms_um']:.3f} um independent "
          f"(E4 floor {paired['e4_local_floor_um']}); running-mean rho "
          f"{structure['running_mean_2mm']['rho']:+.3f} "
          f"p={structure['running_mean_2mm']['p_circular_shift']:.4f}")
    pooled_note = structure["paired_pooled_model"]
    print(f"  pooled-residual rho {pooled_note['rho']:+.3f} against a "
          f"structural bias of {pooled_note['structural_bias']:+.3f}: "
          f"uninformative, as expected")
    print(f"pooled {'M_C' if choose_c else 'M_B'}: H={selected_row['H']} "
          f"Ks={selected_row['K_s']} Kc={selected_row['K_C']} "
          f"lambda={selected_row['lambda']} rms {pooled_selected['rms']:.3f} um, "
          f"LOSO {selected_row['LOSO_RMSE_um']:.3f} um")
    hyper = bayes["hyperparameters"]
    print(f"empirical Bayes: tau2 {hyper['tau2']:.4g}, gamma {hyper['gamma']:.2f}, "
          f"lambda_s {hyper['lambda_slow']:.4g}, lambda_C "
          f"{hyper['lambda_coefficient']:.4g}, edf {hyper['effective_parameters']:.1f}")
    print(f"beta0: E2 {beta0_all['E2']:+.3f}, E3 {beta0_all['E3']:+.3f} um; "
          f"four-session range {reversal_summary['session_spread']['range_um']:.3f} um")
    for name in ("E2", "E3", "pooled"):
        block = leakage[name]
        print(f"leakage {name}: F({block['F']['df1']},{block['F']['df2']}) = "
              f"{block['F']['statistic']:.2f}, p = {block['F']['p_value']:.3g}, "
              f"max |s_up - s_down| "
              f"{block['difference_profile']['max_abs_um']:.3f} um"
              + ("   <- admissible" if name == "E2" else ""))
    for name in ("v3", "v4"):
        raw = iso_metrics_all[name]["raw"]
        after = iso_metrics_all[name]["G_plus_reversal"]
        print(f"ISO {name}: E {raw['E']:.2f} -> {after['E']:.2f}, "
              f"M {raw['M']:.2f} -> {after['M']:.2f}, "
              f"A {raw['A']:.2f} -> {after['A']:.2f} um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
