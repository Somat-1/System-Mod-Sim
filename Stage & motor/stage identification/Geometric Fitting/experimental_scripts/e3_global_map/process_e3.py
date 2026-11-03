#!/usr/bin/env python3
"""Reduce the 3 Nov E3 session and test E2's model on it.

E3 is the first session whose targets E2 never saw. That makes one test
available here that no earlier session could run: take the E2 model exactly as
it was frozen for E5 - before this capture existed - and ask it to predict 400
positions drawn from a different seed. Everything else in this file supports
reading that number: whether the session is valid, whether the drift
correction changes the answer, and what E3's own fit says when it is allowed
to choose its own coefficients.

Both sessions are reduced by the same code (`render_e2_findings`), so a
difference between them is a difference in the stage, not in the pipeline.

    python process_e3.py             # core analysis
    python process_e3.py --select    # also re-run E3's own model selection
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

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

# The model E3 is tested against is the one frozen BEFORE E3 was analysed.
# e5_geometry.FrozenModel defaults to whichever model E5 currently uses, and
# E5 v2 switched that default to an H=2 model whose H was chosen FROM E3, so
# leaving it implicit would quietly make this test circular. Pin it.
E2_FROZEN = E5 / "model" / "geom_model_E2_frozen.json"

ANALYSIS = HERE / "analysis"
SUMMARY = ANALYSIS / "e3_findings_summary.json"

USTEP_PER_MM = 1600.0
SCREW_REV_USTEPS = 3200

# E2's published selection, one-SE rule, from e2_findings_summary.json.
M_B = {"model": "M_B", "H": 2, "K_s": 2, "K_C": "stationary", "lambda": ""}
M_C = {"model": "M_C", "H": 4, "K_s": 2, "K_C": 8, "lambda": 1.0}


# ---------------------------------------------------------------- reduction

def reduce_session(module, label: str) -> dict:
    """Run one session through E2's reduction, with its own log attached."""
    saved = e2f.LOG
    e2f.LOG = module.log_path() if hasattr(module, "log_path") else saved
    try:
        run = module.load(use_cache=True)
        quality = e2f.record_quality(run)
        data = e2f.preprocess(run, quality)
        serial = e2f.parse_serial_checks()
    finally:
        e2f.LOG = saved
    data["label"] = label
    data["quality"] = quality
    data["serial"] = serial
    data["run"] = run
    return data


def domain_of(data: dict) -> tuple[float, float]:
    return float(data["x"].min()), float(data["x"].max())


# ---------------------------------------------------------------- validity

def validity_rows(data: dict) -> list[dict]:
    """The checks that decide whether the session may be interpreted at all."""
    quality, serial, run = data["quality"], data["serial"], data["run"]
    rows: list[dict] = []

    def add(check, value, criterion, ok):
        rows.append({"check": check, "value": value, "criterion": criterion,
                     "pass": bool(ok)})

    matched = int(run.matched.sum())
    add("stops matched to the plan", matched, run.plan.target_mm.size,
        matched == run.plan.target_mm.size)
    add("targets with both directions", data["x"].size, 400, data["x"].size == 400)
    add("missing targets", len(data["missing"]), 0, not data["missing"])

    # MSCNT indexes one electrical period in 1024 counts, so one ustep is 16
    # counts and the 64-ustep grid is exactly one full wrap. Every map target
    # therefore lands on the same count - that is the whole point of the grid -
    # and the useful checks are that it is constant, not what value it takes.
    reference_mscnt = set(serial["reference_mscnt"].tolist())
    add("reference MSCNT identical (no lost steps)", sorted(reference_mscnt),
        "one value", len(reference_mscnt) == 1)
    map_mscnt = set(serial["map_mscnt"].tolist())
    add("map MSCNT identical (all targets on one electrical phase)",
        sorted(map_mscnt), "one value", len(map_mscnt) == 1)
    if len(reference_mscnt) == 1 and len(map_mscnt) == 1:
        # x_ref is 98160 usteps, 48 past a grid point, so the reference must sit
        # a fixed 48 x 16 = 768 counts from the map phase. Finding that exact
        # offset ties the references to the map through the command chain.
        offset = (next(iter(reference_mscnt)) - next(iter(map_mscnt))) % 1024
        add("reference-to-map MSCNT offset", offset,
            f"{(-48 * 16) % 1024} from x_ref's 48-ustep grid offset",
            offset == (-48 * 16) % 1024)

    add("serial stop count", serial["stop_count"], run.plan.target_mm.size,
        serial["stop_count"] == run.plan.target_mm.size)
    add("firmware summary", serial["summary"], "complete=1, aborted=0",
        "complete=1" in serial["summary"] and "aborted=0" in serial["summary"])

    worst_std = float(np.max(quality["std_um"]))
    add("worst record-window std [um]", round(worst_std, 4), "< 0.20",
        worst_std < 0.20)
    worst_slope = float(np.max(np.abs(quality["slope_um_s"])))
    add("worst record-window slope [um/s]", round(worst_slope, 4), "< 0.50",
        worst_slope < 0.50)
    saturated = int(np.sum(quality["saturated"]))
    add("saturated increments", saturated, 0, saturated == 0)

    reference = data["reference_error_um"]
    spread = float(np.ptp(reference))
    add("reference spread over the session [um]", round(spread, 3),
        "reported, not a gate", True)
    end = data["anchor_error_um"]
    add("drift datum to end [um]", round(float(end[-1] - end[0]), 3),
        "reported, not a gate", True)
    return rows


def drift_sensitivity(data: dict) -> dict:
    """How much the drift variant moves the quantity each result rests on."""
    d0, d1, d2 = data["D0"], data["D1"], data["D2"]
    visit = 0.5 * (data["up_time_s"] + data["down_time_s"])
    return {
        "anchor_time_s": [round(v, 2) for v in data["anchor_time_s"].tolist()],
        "anchor_error_um": [round(v, 4) for v in data["anchor_error_um"].tolist()],
        "D0_vs_D1_max_ybar_um": float(np.max(np.abs(d0["ybar"] - d1["ybar"]))),
        "D1_vs_D2_max_ybar_um": float(np.max(np.abs(d1["ybar"] - d2["ybar"]))),
        "D0_vs_D2_max_ybar_um": float(np.max(np.abs(d0["ybar"] - d2["ybar"]))),
        "D0_vs_D2_max_B_um": float(np.max(np.abs(d0["B"] - d2["B"]))),
        "mean_visit_time_spread_s": float(np.ptp(visit)),
        "up_to_down_gap_s": [float(np.min(data["down_time_s"] - data["up_time_s"])),
                             float(np.max(data["down_time_s"] - data["up_time_s"]))],
    }


# ------------------------------------------------- the out-of-sample test

def frozen_prediction(data: dict) -> dict:
    """E2's frozen model, evaluated on E3's own positions.

    The model file was written for E5 and hashed before this session ran, so
    nothing here can have tuned it. `G_rel` is the model's error relative to
    its own intercept, which is what `ybar` is too once centred.
    """
    model = FrozenModel(E2_FROZEN)
    x, theta = data["x"], data["theta"]
    inside = model.in_domain(x)
    predicted = model.G_rel(x[inside], theta[inside])
    measured = data["D2"]["ybar"][inside]
    predicted = predicted - predicted.mean()
    measured = measured - measured.mean()
    residual = measured - predicted
    out = {
        "model_sha256_prefix": model.sha256[:12] if hasattr(model, "sha256") else "",
        "targets_in_domain": int(inside.sum()),
        "targets_outside_domain": int((~inside).sum()),
        "measured_rms_um": float(np.sqrt(np.mean(measured**2))),
        "predicted_rms_um": float(np.sqrt(np.mean(predicted**2))),
        "residual_rms_um": float(np.sqrt(np.mean(residual**2))),
        "residual_max_um": float(np.max(np.abs(residual))),
        "variance_explained": float(1.0 - np.var(residual) / np.var(measured)),
        "correlation": float(np.corrcoef(measured, predicted)[0, 1]),
    }
    out["reduction_factor"] = out["measured_rms_um"] / out["residual_rms_um"]
    out["_x"] = x[inside]
    out["_theta"] = theta[inside]
    out["_measured"] = measured
    out["_predicted"] = predicted
    out["_residual"] = residual
    return out


def own_fit(data: dict, row: dict, variant: str = "D2") -> dict:
    return e2f.fit_selected(data, row, variant=variant)


def harmonic_comparison(model_a: dict, model_b: dict, grid: np.ndarray,
                        harmonics: int) -> list[dict]:
    """Amplitude and phase of each harmonic along the travel, both sessions."""
    rows = []
    for k in range(1, harmonics + 1):
        a = e2f.harmonic_profile(model_a, grid, k)
        b = e2f.harmonic_profile(model_b, grid, k)
        reference = e2f.phase_reference(a)
        phase_a = np.asarray([e2f.nearest_unwrap(p, reference) for p in a["phase"]])
        phase_b = np.asarray([e2f.nearest_unwrap(p, reference) for p in b["phase"]])
        rows.append({
            "harmonic": k,
            "E2_mean_amplitude_um": float(np.mean(a["amplitude"])),
            "E3_mean_amplitude_um": float(np.mean(b["amplitude"])),
            "amplitude_ratio": float(np.mean(b["amplitude"]) / np.mean(a["amplitude"])),
            "max_amplitude_difference_um": float(np.max(np.abs(b["amplitude"] - a["amplitude"]))),
            "mean_phase_shift_deg": float(np.degrees(np.mean(phase_b - phase_a))),
            "max_phase_shift_deg": float(np.max(np.abs(np.degrees(phase_b - phase_a)))),
        })
    return rows


def transfer_table(d2: dict, d3: dict) -> list[dict]:
    """E2 fitted on E2, predicting E3, at a few harmonic orders.

    This is what decides whether the extra harmonics E2's cross-validation
    selected are worth carrying: within a session they lower the residual, and
    the question is whether they lower it on a session they have not seen.
    """
    model = FrozenModel(E2_FROZEN)
    inside = model.in_domain(d3["x"])
    xi, ti = d3["x"][inside], d3["theta"][inside]
    yi = d3["D2"]["ybar"][inside]
    yi = yi - yi.mean()
    dom2 = domain_of(d2)
    rows = []
    for name, H, Ks, Kc, lam in [("M_B  H2, stationary", 2, 2, None, None),
                                 ("H2, varying", 2, 2, 8, 1.0),
                                 ("H3, varying", 3, 2, 8, 1.0),
                                 ("M_C  H4, varying (selected)", 4, 2, 8, 1.0)]:
        fit = e2f.fit_linear(d2["x"], d2["theta"], d2["D2"]["ybar"], H, Ks, Kc,
                             dom2, robust=True, penalty=lam)
        predicted = e2f.predict(fit, xi, ti)
        predicted = predicted - predicted.mean()
        rows.append({"name": name, "edf": round(float(fit["edf"]), 2),
                     "e2_rms_um": round(float(fit["rms"]), 4),
                     "e3_residual_um": round(
                         float(np.sqrt(np.mean((yi - predicted)**2))), 4)})
    return rows


def reversal(data: dict, variant: str = "D2") -> dict:
    fit = e2f.reversal_fit(data, variant)
    B = data[variant]["B"]
    return {"beta0_um": fit["mean_um"], "beta0_sigma_um": fit["mean_sigma_um"],
            "slope_um_per_100mm": fit["slope_um_per_100mm"],
            "slope_sigma_um_per_100mm": fit["slope_sigma_um_per_100mm"],
            "harmonic_amplitude_um": float(np.hypot(fit["beta"][2], fit["beta"][3])),
            "B_mean_um": float(np.mean(B)), "B_sd_um": float(np.std(B, ddof=1)),
            "B_min_um": float(np.min(B)), "B_max_um": float(np.max(B)),
            "fit": fit}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def jsonable(value):
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items() if not k.startswith("_")}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def main() -> int:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    select = "--select" in sys.argv

    print("reducing E2 ...", flush=True)
    d2 = reduce_session(e2_run, "E2 (seed 5315)")
    print("reducing E3 ...", flush=True)
    d3 = reduce_session(e3_run, "E3 (seed 9101, new positions)")

    validity = validity_rows(d3)
    write_csv(ANALYSIS / "e3_validity_table.csv", validity)
    failed = [row["check"] for row in validity if not row["pass"]]

    drift = drift_sensitivity(d3)
    frozen = frozen_prediction(d3)

    grid = np.linspace(max(d2["x"].min(), d3["x"].min()),
                       min(d2["x"].max(), d3["x"].max()), 400)
    model2_c = own_fit(d2, M_C)
    model3_c = own_fit(d3, M_C)
    model3_b = own_fit(d3, M_B)

    harmonics = harmonic_comparison(model2_c, model3_c, grid, int(M_C["H"]))
    write_csv(ANALYSIS / "e3_harmonic_comparison.csv", harmonics)

    transfer = transfer_table(d2, d3)
    write_csv(ANALYSIS / "e3_transfer_table.csv", transfer)

    rev2, rev3 = reversal(d2), reversal(d3)

    if select:
        print("running E3's own model selection (slow) ...", flush=True)
        rows, best_b, best_c = e2f.blocked_cv(d3)
        write_csv(ANALYSIS / "e3_model_selection.csv", rows)
        selection = {"M_B": best_b, "M_C": best_c}
    else:
        # The selection takes minutes; a plain run keeps the stored result
        # rather than silently dropping it from the summary.
        previous = (json.loads(SUMMARY.read_text(encoding="utf-8"))
                    if SUMMARY.exists() else {})
        selection = previous.get("own_selection", {})

    summary = {
        "session": "E3_global_map",
        "capture": "E3_Experimental_Run.csv",
        "log": e3_run.log_path().name,
        "axis_sign": float(np.load(e3_run.CACHE)["axis_sign"]),
        "samples": int(d3["run"].t.size),
        "duration_min": float(d3["run"].t[-1] / 60.0),
        "datum_mm": float(d3["run"].datum_mm),
        "validity": validity,
        "validity_failures": failed,
        "drift": drift,
        "frozen_E2_prediction_on_E3": frozen,
        "harmonic_comparison": harmonics,
        "transfer_table": transfer,
        "fits": {
            "E2_M_C": {"rms_um": model2_c["rms"], "edf": model2_c["edf"],
                       "p": model2_c["p"], "gain_ppm": model2_c.get("gain_ppm")},
            "E3_M_C": {"rms_um": model3_c["rms"], "edf": model3_c["edf"],
                       "p": model3_c["p"], "gain_ppm": model3_c.get("gain_ppm")},
            "E3_M_B": {"rms_um": model3_b["rms"], "edf": model3_b["edf"],
                       "p": model3_b["p"]},
        },
        "reversal": {"E2": {k: v for k, v in rev2.items() if k != "fit"},
                     "E3": {k: v for k, v in rev3.items() if k != "fit"}},
        "electrical_phase": {
            "E2_map_mscnt": sorted(set(d2["serial"]["map_mscnt"].tolist())),
            "E3_map_mscnt": sorted(set(d3["serial"]["map_mscnt"].tolist())),
            "E2_reference_mscnt": sorted(set(d2["serial"]["reference_mscnt"].tolist())),
            "E3_reference_mscnt": sorted(set(d3["serial"]["reference_mscnt"].tolist())),
            "session_shift_counts": int(
                (next(iter(set(d3["serial"]["map_mscnt"].tolist())))
                 - next(iter(set(d2["serial"]["map_mscnt"].tolist())))) % 1024),
        },
        "own_selection": selection,
    }
    SUMMARY.write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")

    print(f"\nvalidity: {len(validity) - len(failed)}/{len(validity)} pass"
          + (f"  FAILED: {failed}" if failed else ""))
    print(f"frozen E2 model on E3's new positions:")
    print(f"  measured rms {frozen['measured_rms_um']:.3f} um -> "
          f"residual rms {frozen['residual_rms_um']:.3f} um "
          f"({frozen['reduction_factor']:.2f}x, R2 {frozen['variance_explained']:.3f})")
    print(f"E3 own fit M_C rms {model3_c['rms']:.4f} um (E2 {model2_c['rms']:.4f})")
    print(f"reversal beta0: E2 {rev2['beta0_um']:+.3f}  E3 {rev3['beta0_um']:+.3f} um")
    print(f"wrote {SUMMARY.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
