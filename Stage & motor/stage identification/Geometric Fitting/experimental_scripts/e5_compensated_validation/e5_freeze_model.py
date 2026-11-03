#!/usr/bin/env python3
"""Freeze the E2-only geometric model for E5.

Refits exactly the model E2 selected (render_e2_findings: M_C, H = 4, K_s = 2,
K_C = 8, lambda = 1, drift variant D2) on the E2 capture, together with the
D2 reversal fit, and writes every number the standalone evaluator needs to
model/geom_model_E2_frozen.json. The file is written only after
e5_geometry.FrozenModel reproduces the E2 predict() on a dense (x, theta)
grid to 1e-9 um, so E5 never depends on the E2 pipeline again.

The model is E2-only by design: E5 must stay an independent test, so nothing
from E3 enters it. Rerunning this script reproduces the same file; its hash
is what e5_make_table.py records in the command table and the firmware log.

    DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib \
      /opt/homebrew/bin/python3 e5_freeze_model.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
E2_DIR = HERE.parent / "e2_global_map"
sys.path.insert(0, str(E2_DIR))
import e2_run  # noqa: E402
import render_e2_findings as e2f  # noqa: E402

sys.path.insert(0, str(HERE))
from e5_geometry import MODEL_FILE, FrozenModel, theta_of  # noqa: E402

VARIANT = "D2"


def selected_row() -> dict:
    """The M_C row E2 selected, read back from its own model-selection table."""
    summary = json.loads((E2_DIR / "analysis" / "e2_findings_summary.json").read_text())
    chosen = summary["model_selection"]["M_C_one_se"]
    assert summary["model_selection"]["selected"] == "M_C", "E2 did not select M_C"
    return {"H": chosen["H"], "K_s": chosen["K_s"], "K_C": chosen["K_C"],
            "lambda": chosen["lambda"]}, summary


def main() -> int:
    row, summary = selected_row()
    run = e2_run.load()
    data = e2f.preprocess(run, e2f.record_quality(run))
    model = e2f.fit_selected(data, row, variant=VARIANT)
    reversal = e2f.reversal_fit(data, VARIANT)
    s = model["state"]
    assert abs(model["rms"] - summary["selected_model"]["residual_rms_um"]) < 1e-9, \
        "refit does not reproduce the E2 selected-model residual"
    assert abs(reversal["mean_um"] - summary["reversal"][VARIANT]["mean_um"]) < 1e-9

    frozen = {
        "experiment": "E2_global_map",
        "purpose": "E5 command pre-correction; frozen before E3 is analysed",
        "source_log": e2f.LOG.name,
        "source_capture": "E2_Experimental_Run.csv",
        "variant": VARIANT,
        "selection": row,
        "convention": summary["convention"],
        "state": {"xbar": s["xbar"], "domain": list(s["domain"]), "H": s["H"],
                  "Ks": s["Ks"], "Kc": s["Kc"],
                  "slow_projection": s["slow_projection"].tolist(),
                  "slow_transform": s["slow_transform"].tolist(),
                  "names": s["names"]},
        "beta": model["beta"].tolist(),
        "beta_sigma": np.sqrt(np.diag(model["cov"])).tolist(),
        "cov": model["cov"].tolist(),
        "fit": {"residual_rms_um": model["rms"], "edf": model["edf"],
                "sigma_um": model["sigma"],
                "cv_rmse_um": summary["model_selection"]["M_C_one_se"]["CV_RMSE_um"]},
        "gain_ppm": float(model["beta"][1]) * 1000.0,
        "reversal": {"beta0_um": reversal["mean_um"],
                     "beta0_sigma_um": reversal["mean_sigma_um"],
                     "slope_um_per_100mm": reversal["slope_um_per_100mm"]},
    }
    MODEL_FILE.parent.mkdir(parents=True, exist_ok=True)
    candidate = MODEL_FILE.with_suffix(".candidate.json")
    candidate.write_text(json.dumps(frozen, indent=1) + "\n", encoding="utf-8")

    check = FrozenModel(candidate)
    x = np.linspace(s["domain"][0] - 3.0, s["domain"][1] + 3.0, 4001)
    theta = theta_of(np.rint(x * 1600).astype(int))
    error = float(np.max(np.abs(check.G(x, theta) - e2f.predict(model, x, theta))))
    if error > 1e-9:
        candidate.unlink()
        raise SystemExit(f"standalone evaluator disagrees with E2 predict(): {error:.3e} um")
    candidate.replace(MODEL_FILE)
    final = FrozenModel(MODEL_FILE)

    # The paired E2 targets and the fitted values, for the validation figures.
    with (HERE / "model" / "e2_fit_points.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["x_mm", "usteps", "ybar_um", "fit_um", "B_um"])
        fit = e2f.predict(model, data["x"], data["theta"])
        for i in range(data["x"].size):
            w.writerow([f"{data['x'][i]:.6f}", int(data["usteps"][i]),
                        f"{data[VARIANT]['ybar'][i]:.6f}", f"{fit[i]:.6f}",
                        f"{data[VARIANT]['B'][i]:.6f}"])
    print(f"frozen {MODEL_FILE.name}: sha256 {final.sha256[:12]}  "
          f"max |standalone - E2| = {error:.1e} um  "
          f"rms {model['rms']:.3f} um  beta0 {reversal['mean_um']:+.3f} um  "
          f"gain {frozen['gain_ppm']:+.2f} ppm  domain {s['domain'][0]:.3f}-{s['domain'][1]:.3f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
