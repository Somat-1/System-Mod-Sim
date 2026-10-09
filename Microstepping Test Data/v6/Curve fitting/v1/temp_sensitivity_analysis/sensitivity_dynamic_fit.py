"""Initial-testing sensitivity analysis for the v6 curve-fitting "dynamic" pass.

Same type of analysis as `v3/Parameter Optimization/calibration_bracketing/
step6_sensitivity_table.py` (see its SENSITIVITY_ANALYSIS.md), adapted from a
LuGre friction-model identification (parameters perturbed in a nonlinear
time-domain simulation) to the much cheaper "dynamic" smoothing filter used
in `fit_v6_curves.py` (median despike + centered Gaussian, applied directly
to already-baselined `measured` arrays already cached in each MRES group's
`fitted_blocks.npz` -- no need to reparse the raw EL5101 capture).

Parameters perturbed independently, +-10% each, holding the other fixed:
  - median_ms  (despike window)
  - sigma_ms   (Gaussian sigma)

Cost metric, per (segment, tag) job:
    cost = RMS( measured_um(t) - dynamic_fit_um(t) )
over the full segment window (matches v3's per-block RMS definition).

Sensitivity, per (segment, parameter) cell, from the +-10% perturbation:
    sensitivity = 0.5 * ( |cost(+10%) - cost(base)| + |cost(-10%) - cost(base)| )
                  / cost(base)
a normalized, symmetric relative sensitivity -- identical formula to v3.
NOISE_FLOOR = 0.01 (1%): at or below this, the fit's RMS-to-measured is
"insensitive" to that +-10% change; above it, "sensitive".

Scope (this is initial testing, not the full sequence -- picked to contrast
the two ends of the ramp-dynamic clip range plus one MRES 1 vs MRES 32
oscillation comparison):
  - osc_MRES1     oscillation dynamic fit, MRES 1  (fixed params: 3.0/2.0 ms)
  - osc_MRES32    oscillation dynamic fit, MRES 32 (same fixed params, but
                  the commanded step is ~30x smaller relative to sensor noise)
  - leg_MRES1_Dmod   ramp dynamic fit, MRES 1 D-mod  (mid-range, unclipped:
                      3.57/2.14 ms -- ramp_dynamic_params picked this freely)
  - leg_MRES32_Dfast ramp dynamic fit, MRES 32 D-fast (floored at the hard
                      lower bound: 1.0/0.5 ms -- see CURVE_FITTING_PROCESS.md
                      Sec on MRES 32 flooring out)

For the floored segment, the -10% perturbation is evaluated unclipped (i.e.
below the production RAMP_DYNAMIC_*_BOUNDS_MS floor) specifically to probe
the filter's *inherent* sensitivity there, separate from the clipping that
the real pipeline would apply. This is flagged explicitly in the output.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CURVE_FIT_V1 = HERE.parent
sys.path.insert(0, str(CURVE_FIT_V1 / "scripts"))
from fit_v6_curves import smooth, OSC_DYNAMIC_MEDIAN_MS, OSC_DYNAMIC_SIGMA_MS  # noqa: E402

PERTURB_FRAC = 0.10
NOISE_FLOOR = 0.01
DT_S = 1.0e-3  # 1 kHz EL5101 capture, confirmed via fit_summary.json for both groups


def load_leg_params(mres, label):
    summary = json.loads((CURVE_FIT_V1 / f"MRES_{mres}" / "fit_summary.json").read_text())
    p = summary["trajectory_legs"]["dynamic"]["params_by_leg_ms"][label]
    return float(p["median_ms"]), float(p["sigma_ms"])


def load_segment(mres, kind, label=None):
    npz = np.load(CURVE_FIT_V1 / f"MRES_{mres}" / "fitted_blocks.npz")
    if kind == "oscillation":
        t = npz["oscillation_time_s"]
        measured = npz["oscillation_measured_um"]
        median_ms, sigma_ms = OSC_DYNAMIC_MEDIAN_MS, OSC_DYNAMIC_SIGMA_MS
    else:
        t = npz[f"leg_{label}_time_s"]
        measured = npz[f"leg_{label}_measured_um"]
        median_ms, sigma_ms = load_leg_params(mres, label)
    return t, measured, median_ms, sigma_ms


SEGMENTS = [
    {"name": "osc_MRES1", "mres": 1, "kind": "oscillation", "label": None,
     "note": "oscillation dynamic fit, fixed 3.0/2.0 ms filter"},
    {"name": "osc_MRES32", "mres": 32, "kind": "oscillation", "label": None,
     "note": "oscillation dynamic fit, same fixed filter, ~30x finer commanded step"},
    {"name": "leg_MRES1_Dmod", "mres": 1, "kind": "leg", "label": "D-mod",
     "note": "ramp dynamic fit, unclipped mid-range filter (3.57/2.14 ms)"},
    {"name": "leg_MRES32_Dfast", "mres": 32, "kind": "leg", "label": "D-fast",
     "note": "ramp dynamic fit, floored at the hard lower bound (1.0/0.5 ms)"},
]


def rms_cost(measured, fitted):
    return float(np.sqrt(np.mean((measured - fitted) ** 2)))


def run_segment(seg):
    t, measured, base_median_ms, base_sigma_ms = load_segment(
        seg["mres"], seg["kind"], seg["label"])
    base_fit = smooth(measured, DT_S, median_ms=base_median_ms, sigma_ms=base_sigma_ms)
    base_cost = rms_cost(measured, base_fit)

    jobs = {"base": {"median_ms": base_median_ms, "sigma_ms": base_sigma_ms,
                     "cost_um": base_cost}}
    cells = {}
    for param, base_val in (("median_ms", base_median_ms), ("sigma_ms", base_sigma_ms)):
        costs = {}
        for tag, sign in (("plus", +1.0), ("minus", -1.0)):
            perturbed_val = base_val * (1.0 + sign * PERTURB_FRAC)
            kwargs = {"median_ms": base_median_ms, "sigma_ms": base_sigma_ms}
            kwargs[param] = perturbed_val
            fitted = smooth(measured, DT_S, **kwargs)
            cost = rms_cost(measured, fitted)
            costs[tag] = cost
            jobs[f"{tag}|{param}"] = {param: perturbed_val, "cost_um": cost}
        sensitivity = 0.5 * (abs(costs["plus"] - base_cost)
                             + abs(costs["minus"] - base_cost)) / base_cost
        cells[param] = {
            "base_value_ms": base_val,
            "cost_plus10_um": costs["plus"],
            "cost_minus10_um": costs["minus"],
            "sensitivity": sensitivity,
            "verdict": "sensitive" if sensitivity > NOISE_FLOOR else "insensitive",
        }
    return {
        "segment": seg["name"], "mres": seg["mres"], "kind": seg["kind"],
        "label": seg["label"], "note": seg["note"],
        "n_samples": int(measured.size), "duration_s": float(t[-1] - t[0]),
        "base_cost_um": base_cost, "cells": cells, "jobs": jobs,
        "_t": t, "_measured": measured, "_base_fit": base_fit,
        "_base_median_ms": base_median_ms, "_base_sigma_ms": base_sigma_ms,
    }


def make_overview_plot(results, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    names = [r["segment"] for r in results]
    x = np.arange(len(names))
    width = 0.35
    for ax, param, title in zip(axes, ("median_ms", "sigma_ms"),
                                ("Despike window (median_ms)", "Gaussian sigma (sigma_ms)")):
        vals = [r["cells"][param]["sensitivity"] for r in results]
        colors = ["#c0392b" if v > NOISE_FLOOR else "#2a7f62" for v in vals]
        ax.bar(x, vals, color=colors, width=width)
        ax.axhline(NOISE_FLOOR, color="#555555", linewidth=0.8, linestyle="--",
                  label=f"noise floor ({NOISE_FLOOR:.0%})")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("relative sensitivity")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(True, axis="y", alpha=0.3)
    fig.suptitle(f"Dynamic-fit sensitivity to ±{PERTURB_FRAC:.0%} filter-parameter "
                f"perturbation (RMS-to-measured cost)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def make_floor_case_plot(result, out_path):
    """The floored segment's absolute RMS cost is ~0.06 um -- already near
    the encoder noise floor, so a raw measured-vs-fit overlay makes any
    perturbation invisible regardless of its relative sensitivity. Plot the
    fit-vs-fit residual (perturbed minus base) instead, which is where the
    effect actually lives."""
    t, measured, base_fit = result["_t"], result["_measured"], result["_base_fit"]
    base_median_ms, base_sigma_ms = result["_base_median_ms"], result["_base_sigma_ms"]
    plus_median = smooth(measured, DT_S, median_ms=base_median_ms * 1.10, sigma_ms=base_sigma_ms)
    minus_median = smooth(measured, DT_S, median_ms=base_median_ms * 0.90, sigma_ms=base_sigma_ms)
    plus_sigma = smooth(measured, DT_S, median_ms=base_median_ms, sigma_ms=base_sigma_ms * 1.10)
    minus_sigma = smooth(measured, DT_S, median_ms=base_median_ms, sigma_ms=base_sigma_ms * 0.90)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    axes[0].plot(t, plus_median - base_fit, color="#c0392b", lw=1.0, label="+10% median - base")
    axes[0].plot(t, minus_median - base_fit, color="#2a7f62", lw=1.0, label="-10% median - base")
    axes[0].set_title(f"median_ms perturbation (base {base_median_ms:.2f} ms; "
                      f"{max(1, round(base_median_ms))}-sample window either way)", fontsize=9.5)
    axes[1].plot(t, plus_sigma - base_fit, color="#c0392b", lw=1.0, label="+10% sigma - base")
    axes[1].plot(t, minus_sigma - base_fit, color="#2a7f62", lw=1.0, label="-10% sigma - base")
    axes[1].set_title(f"sigma_ms perturbation (base {base_sigma_ms:.2f} ms "
                      f"= {base_sigma_ms * 1e-3 / DT_S:.2f} samples)", fontsize=9.5)
    for ax in axes:
        ax.axhline(0.0, color="#333333", lw=0.6)
        ax.set_xlabel("Time (s)", fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7.5)
    axes[0].set_ylabel("Fit change vs base (µm)", fontsize=8)
    fig.suptitle(f"{result['segment']} — floored filter (1.0/0.5 ms): fit-vs-fit "
                f"residual under ±10% perturbation (base RMS-to-measured = "
                f"{result['base_cost_um']:.4f} µm)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main():
    results = [run_segment(seg) for seg in SEGMENTS]

    print(f"{'segment':<20}{'param':<12}{'base_ms':>9}{'cost_um':>10}"
          f"{'sens':>9}  verdict")
    for r in results:
        for param in ("median_ms", "sigma_ms"):
            c = r["cells"][param]
            print(f"{r['segment']:<20}{param:<12}{c['base_value_ms']:>9.3f}"
                  f"{r['base_cost_um']:>10.4f}{c['sensitivity']:>9.4f}  {c['verdict']}")

    make_overview_plot(results, HERE / "plots" / "sensitivity_overview.png")
    floored = next(r for r in results if r["segment"] == "leg_MRES32_Dfast")
    make_floor_case_plot(floored, HERE / "plots" / "floor_case_perturbation.png")

    out = {
        "schema": "v6-curve-fit-dynamic-sensitivity-v1",
        "perturb_frac": PERTURB_FRAC,
        "noise_floor": NOISE_FLOOR,
        "dt_s": DT_S,
        "segments": [
            {k: v for k, v in r.items() if not k.startswith("_")}
            for r in results
        ],
    }
    (HERE / "sensitivity_results.json").write_text(
        json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote sensitivity_results.json, plots/sensitivity_overview.png, "
          f"plots/floor_case_perturbation.png -> {HERE}")


if __name__ == "__main__":
    main()
