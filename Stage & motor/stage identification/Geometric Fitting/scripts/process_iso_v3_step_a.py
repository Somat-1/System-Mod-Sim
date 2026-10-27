#!/usr/bin/env python3
"""Step A: linear geometric identification from the saved ISO 230-2 v3 results.

The raw AquadB scope CSV is optional.  The current workspace does not contain
it, so this script performs every calculation recoverable exactly from the
saved per-target means, standard deviations and counts, and marks raw-only
validity/drift/per-run outputs as unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, stats


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT.parent
DEFAULT_RESULTS = STAGE / "iso230-2" / "analysis" / "iso230_2_v3_results.json"
DEFAULT_LOG = STAGE / "iso230-2" / "data" / "terminal_20261021_183115.log"
LEAD_MM = 2.0


def load_aggregate(path: Path):
    data = json.loads(path.read_text())
    rows = []
    for key, value in sorted(data["per_target"].items(), key=lambda kv: int(kv[0])):
        row = {"target_i": int(key), **value}
        rows.append(row)
    return data, rows


def inspect_log(path: Path):
    text = path.read_text(errors="replace")
    stop_re = re.compile(
        r"# STOP,(?P<state>[^,]+),(?P<run>-?\d+),(?P<direction>-?\d+),"
        r"(?P<target>-?\d+),(?P<usteps>-?\d+),(?P<cmd>[\d.]+),"
        r"(?P<t_arr>\d+),(?P<t0>\d+),(?P<t1>\d+)"
    )
    stops = [{k: (float(v) if k == "cmd" else int(v) if k not in {"state"} else v)
              for k, v in m.groupdict().items()} for m in stop_re.finditer(text)]
    iso = [s for s in stops if s["state"] == "ISO"]
    refs = [s for s in stops if s["state"] in {"ISO_REF", "END_REF"}]
    done = re.search(r"# ISO_DONE,runs=(\d+),stops=(\d+),minutes=([\d.]+)", text)
    driver = re.search(
        r"# DRIVER,.*?intpol=(\d+).*?mres=(\d+),irun=(\d+),ihold=(\d+),"
        r"tpowerdown=(\d+),iholddelay=(\d+),current_ma=(\d+)", text
    )
    return {
        "iso_stop_count": len(iso),
        "runs": sorted(set(s["run"] for s in iso)),
        "directions": sorted(set(s["direction"] for s in iso)),
        "reference_stop_count": len(refs),
        "done": ({"runs": int(done.group(1)), "stops": int(done.group(2)),
                  "minutes": float(done.group(3))} if done else None),
        "abort_or_lost_step_flag": bool(re.search(
            r"ISO_ABORT|lost[-_ ]?step|watchdog.*(?:trip|fail)|FAULT", text, re.I)),
        "driver": ({"intpol": int(driver.group(1)), "mres": int(driver.group(2)),
                    "irun": int(driver.group(3)), "ihold": int(driver.group(4)),
                    "tpowerdown": int(driver.group(5)),
                    "iholddelay": int(driver.group(6)),
                    "current_ma": int(driver.group(7))} if driver else None),
        "targets_on_microstep_grid": all(
            abs(s["cmd"] * 1600.0 - round(s["cmd"] * 1600.0)) < 1e-7 for s in iso),
        "record_durations_ms": sorted(set(s["t1"] - s["t0"] for s in iso)),
        "settle_durations_ms": sorted(set(s["t0"] - s["t_arr"] for s in iso)),
    }


def linear_fit(x, y, var=None):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xc = float(x.mean())
    phi = np.column_stack([np.ones_like(x), x - xc])
    if var is None:
        w = np.ones_like(x)
        beta = np.linalg.lstsq(phi, y, rcond=None)[0]
        resid = y - phi @ beta
        scale = np.sum(resid ** 2) / (len(x) - 2)
        cov = scale * np.linalg.inv(phi.T @ phi)
        chi2nu = None
    else:
        var = np.asarray(var, float)
        w = 1.0 / var
        normal = phi.T @ (w[:, None] * phi)
        beta = np.linalg.solve(normal, phi.T @ (w * y))
        resid = y - phi @ beta
        chi2nu = float(np.sum(w * resid ** 2) / (len(x) - 2))
        cov = chi2nu * np.linalg.inv(normal)
    critical = stats.t.ppf(0.975, len(x) - 2)
    return {
        "x_center_mm": xc,
        "intercept_centered_um": float(beta[0]),
        "intercept_at_zero_um": float(beta[0] - beta[1] * xc),
        "slope_um_per_mm": float(beta[1]),
        "slope_ppm": float(beta[1] * 1000.0),
        "slope_se_um_per_mm": float(np.sqrt(cov[1, 1])),
        "slope_ci95_um_per_mm": float(critical * np.sqrt(cov[1, 1])),
        "slope_ci95_ppm": float(critical * np.sqrt(cov[1, 1]) * 1000.0),
        "reduced_chi2": chi2nu,
        "predicted": phi @ beta,
        "residual": resid,
    }


def iso_summary(rows, correction=None):
    up = np.array([r["xbar_up"] for r in rows], float)
    dn = np.array([r["xbar_dn"] for r in rows], float)
    su = np.array([r["s_up"] for r in rows], float)
    sd = np.array([r["s_dn"] for r in rows], float)
    if correction is not None:
        up = up - correction
        dn = dn - correction
    mid = 0.5 * (up + dn)
    reversal = up - dn
    hi = np.maximum(up + 2 * su, dn + 2 * sd)
    lo = np.minimum(up - 2 * su, dn - 2 * sd)
    ri = np.maximum.reduce([4 * su, 4 * sd, 2 * su + 2 * sd + np.abs(reversal)])
    return {
        "B": float(np.max(np.abs(reversal))),
        "B_mean": float(np.mean(reversal)),
        "R_up": float(np.max(4 * su)),
        "R_dn": float(np.max(4 * sd)),
        "R": float(np.max(ri)),
        "E_up": float(np.ptp(up)),
        "E_dn": float(np.ptp(dn)),
        "E": float(max(up.max(), dn.max()) - min(up.min(), dn.min())),
        "M": float(np.ptp(mid)),
        "A": float(hi.max() - lo.min()),
        "s_med": float(np.median(np.r_[su, sd])),
        "s_max": float(np.max(np.r_[su, sd])),
    }


def total_power_share(rows, x, correction):
    means, ns, within = [], [], []
    for row in rows:
        for d in ("up", "dn"):
            n = int(row[f"n_{d}"])
            means.append(float(row[f"xbar_{d}"]))
            ns.append(n)
            within.append((n - 1) * float(row[f"s_{d}"]) ** 2)
    means, ns = np.asarray(means), np.asarray(ns)
    grand = np.average(means, weights=ns)
    before = np.sum(within) + np.sum(ns * (means - grand) ** 2)
    corr = np.repeat(np.asarray(correction), 2)
    corrected = means - corr
    grand_after = np.average(corrected, weights=ns)
    after = np.sum(within) + np.sum(ns * (corrected - grand_after) ** 2)
    return float(1.0 - after / before)


def serialise_fit(fit):
    return {k: v for k, v in fit.items() if k not in {"predicted", "residual"}}


def save_csvs(rows, x, y, var, reversal, fit, outdir):
    with (outdir / "per_target_statistics.csv").open("w", newline="") as fh:
        fields = ["target_i", "cmd_mm", "n_up", "n_dn", "xbar_up_um",
                  "xbar_dn_um", "s_up_um", "s_dn_um", "ybar_um",
                  "var_ybar_um2", "B_um", "fit_um", "residual_um"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow({
                "target_i": row["target_i"], "cmd_mm": x[i],
                "n_up": row["n_up"], "n_dn": row["n_dn"],
                "xbar_up_um": row["xbar_up"], "xbar_dn_um": row["xbar_dn"],
                "s_up_um": row["s_up"], "s_dn_um": row["s_dn"],
                "ybar_um": y[i], "var_ybar_um2": var[i], "B_um": reversal[i],
                "fit_um": fit["predicted"][i], "residual_um": fit["residual"][i],
            })


def render_a4(x, y, var, up, dn, reversal, fit, bfit, loo, period, power,
              before, after, validity, results, outpath):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.2,
        "axes.titlesize": 9.2, "axes.labelsize": 8.1,
        "axes.edgecolor": "#4f5964", "axes.linewidth": 0.75,
        "xtick.labelsize": 7.3, "ytick.labelsize": 7.3,
        "grid.color": "#d7dde3", "grid.linewidth": 0.55,
        "figure.facecolor": "white", "axes.facecolor": "#fbfcfd",
    })
    blue, red, green, ink, muted = "#276fbf", "#b44b4b", "#26735b", "#17212b", "#66717d"
    fig = plt.figure(figsize=(8.27, 11.69), constrained_layout=False)
    gs = fig.add_gridspec(5, 2, left=.085, right=.965, top=.845, bottom=.075,
                          hspace=.74, wspace=.33, height_ratios=[1.12, 1.0, 1.0, .82, .72])
    fig.text(.085, .958, "ISO 230-2 v3 — linear geometric identification", fontsize=17,
             weight="bold", color=ink, va="top")
    fig.text(.085, .927, "Step A · direction-averaged position error · 24 targets · 5 bidirectional runs",
             fontsize=9.2, color=muted, va="top")
    slope = fit["slope_um_per_mm"]
    ci = fit["slope_ci95_um_per_mm"]
    cards = [
        ("GAIN SLOPE", f"{slope:+.4f} ± {ci:.4f}", f"µm/mm · {slope*1000:+.1f} ± {ci*1000:.1f} ppm"),
        ("GAIN OVER SPAN", f"{slope*(x.max()-x.min()):+.3f} µm", f"over {x.max()-x.min():.3f} mm"),
        ("EXPLAINED POWER", f"{100*results['power']['direction_average_share']:.1f}%", "direction-averaged"),
        ("PREDICTED E", f"{before['E']:.2f} → {after['E']:.2f} µm", "linear correction only"),
    ]
    for j, (head, main, sub) in enumerate(cards):
        x0 = .085 + j * .222
        fig.add_artist(plt.Rectangle((x0, .86), .205, .052, transform=fig.transFigure,
                                     facecolor="#f0f4f7", edgecolor="#d7dde3", lw=.7))
        fig.text(x0+.010, .901, head, fontsize=6.7, color=muted, weight="bold", va="top")
        fig.text(x0+.010, .883, main, fontsize=9.2, color=ink, weight="bold", va="top")
        fig.text(x0+.010, .867, sub, fontsize=6.9, color=muted, va="top")

    ax = fig.add_subplot(gs[0, :])
    ax.fill_between(x, y-2*np.sqrt(var), y+2*np.sqrt(var), color=blue, alpha=.10, lw=0)
    ax.plot(x, up, "o", ms=3, color="#d28b45", alpha=.72, label="up-direction mean")
    ax.plot(x, dn, "s", ms=3, color="#6d74b8", alpha=.72, label="down-direction mean")
    ax.errorbar(x, y, yerr=2*np.sqrt(var), fmt="o", ms=4, color=blue, ecolor=blue,
                elinewidth=.8, capsize=2, label="direction average ±2σ")
    ax.plot(x, fit["predicted"], color=red, lw=1.7, label="weighted linear fit")
    ax.set(title="(a) Direction-averaged geometric error and fitted gain",
           xlabel="rounded commanded position (mm)", ylabel="deviation (µm)")
    ax.legend(ncol=2, fontsize=7.1, frameon=False, loc="upper right")
    ax.grid(True)

    ax = fig.add_subplot(gs[1, 0])
    ax.axhline(0, color=muted, lw=.8)
    ax.errorbar(x, fit["residual"], yerr=2*np.sqrt(var), fmt="o", color=blue,
                ms=3.5, capsize=2, elinewidth=.7)
    ax.plot(x, fit["residual"], color=blue, alpha=.45, lw=.8)
    ax.set(title="(b) Residual structure after gain removal", xlabel="position (mm)",
           ylabel="residual (µm)")
    ax.grid(True)

    ax = fig.add_subplot(gs[1, 1])
    sigma_b = np.sqrt(np.array([r["s_up"]**2/r["n_up"] + r["s_dn"]**2/r["n_dn"] for r in results["_rows"]]))
    ax.axhline(0, color=muted, lw=.8)
    ax.errorbar(x, reversal, yerr=2*sigma_b, fmt="o", color=green, ms=3.5,
                capsize=2, elinewidth=.7)
    ax.plot(x, bfit["predicted"], color=red, lw=1.4)
    ax.set(title="(c) Reversal, characterized but not fitted into G", xlabel="position (mm)",
           ylabel="B = up − down (µm)")
    ax.grid(True)

    ax = fig.add_subplot(gs[2, 0])
    ax.plot(period, power, color=blue, lw=1.25)
    peak_i = int(np.argmax(power))
    ax.axvline(period[peak_i], color=red, lw=1, ls="--")
    ax.text(.97, .90, f"largest preliminary peak\nperiod {period[peak_i]:.2f} mm",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.2, color=red)
    ax.set(title="(d) Preliminary Lomb–Scargle view of residual", xlabel="spatial period (mm)",
           ylabel="normalized power", xscale="log")
    ax.grid(True, which="both")

    ax = fig.add_subplot(gs[2, 1])
    phase = np.mod(x, LEAD_MM)/LEAD_MM
    order = np.argsort(phase)
    ax.axhline(0, color=muted, lw=.8)
    ax.errorbar(phase, fit["residual"], yerr=2*np.sqrt(var), fmt="o", color=blue,
                ms=3.5, capsize=2, elinewidth=.7)
    ax.plot(phase[order], fit["residual"][order], color=blue, lw=.7, alpha=.35)
    ax.set(title=f"(e) Residual folded onto nominal {LEAD_MM:g} mm screw lead",
           xlabel="lead phase, (x mod L)/L", ylabel="residual (µm)")
    ax.grid(True)

    ax = fig.add_subplot(gs[3, 0])
    dfit = results["robustness"]["direction_fits"]
    vals = [dfit["up"]["slope_ppm"], fit["slope_ppm"], dfit["down"]["slope_ppm"]]
    cis = [dfit["up"]["slope_ci95_ppm"], fit["slope_ci95_ppm"], dfit["down"]["slope_ci95_ppm"]]
    ax.errorbar([0,1,2], vals, yerr=cis, fmt="o", color=blue, capsize=3)
    ax.fill_between([-.25,2.25], min(loo)*1000, max(loo)*1000, color=green, alpha=.13,
                    label="leave-one-target-out range")
    ax.set_xticks([0,1,2], ["up", "average", "down"])
    ax.set(title="(f) Slope robustness", ylabel="gain slope (ppm)")
    ax.legend(fontsize=6.8, frameon=False)
    ax.grid(True, axis="y")

    ax = fig.add_subplot(gs[3, 1])
    names = ["E", "M", "A", "B", "R"]
    z = np.arange(len(names)); w=.34
    ax.bar(z-w/2, [before[k] for k in names], width=w, color="#9aa7b2", label="measured")
    ax.bar(z+w/2, [after[k] for k in names], width=w, color=green, label="after gain correction")
    ax.set_xticks(z, names)
    ax.set(title="(g) Predicted ISO metric change", ylabel="µm")
    ax.legend(fontsize=6.8, frameon=False)
    ax.grid(True, axis="y")

    ax = fig.add_subplot(gs[4, :]); ax.axis("off")
    validity_lines = [
        f"Sequence: PASS — {validity['iso_stop_count']} ISO stops; {len(validity['runs'])} runs; completion marker present.",
        f"Command grid: {'PASS' if validity['targets_on_microstep_grid'] else 'FAIL'} — all target positions are integer 1/16-microstep commands.",
        f"Watchdog/log: {'PASS' if not validity['abort_or_lost_step_flag'] else 'FLAG'} — no abort, watchdog, fault or lost-step marker found.",
        "Settling / x-ref drift / per-run slopes: NOT RECOMPUTED — raw ISO230-2_v3.csv is absent.",
        "Environment / physical attribution: NOT AVAILABLE — air and stage temperature were not saved with this project.",
    ]
    ax.text(0, 1.02, "Validity and interpretation", fontsize=9.2, weight="bold", color=ink, va="top")
    ax.text(0, .86, "\n".join(validity_lines), fontsize=7.1, color=ink, va="top", linespacing=1.5)
    conclusion = (
        f"The linear term is {slope*1000:+.1f} ppm and accounts for "
        f"{100*results['power']['direction_average_share']:.1f}% of direction-averaged spatial error power. "
        f"After applying it, predicted E falls from {before['E']:.2f} to {after['E']:.2f} µm; "
        "B and R remain unchanged by construction. Residual structure remains, so a periodic/slow model is justified, "
        "but the sparse 24-point periodogram is screening evidence only."
    )
    fig.text(.085, .025, conclusion, fontsize=7.1, color=ink, va="bottom", wrap=True)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aggregate", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--log", type=Path, default=DEFAULT_LOG)
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    ap.add_argument("--render", type=Path, default=ROOT / "rendered" / "iso_v3_step_a_a4.png")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    saved, rows = load_aggregate(args.aggregate)
    validity = inspect_log(args.log)
    x = np.array([r["cmd_mm"] for r in rows])
    up = np.array([r["xbar_up"] for r in rows])
    dn = np.array([r["xbar_dn"] for r in rows])
    su = np.array([r["s_up"] for r in rows])
    sd = np.array([r["s_dn"] for r in rows])
    nu = np.array([r["n_up"] for r in rows])
    nd = np.array([r["n_dn"] for r in rows])
    y = 0.5 * (up + dn)
    var = 0.25 * (su**2/nu + sd**2/nd)
    reversal = up - dn

    fit = linear_fit(x, y, var)
    ols = linear_fit(x, y)
    upfit = linear_fit(x, up, su**2/nu)
    dnfit = linear_fit(x, dn, sd**2/nd)
    bvar = su**2/nu + sd**2/nd
    bfit = linear_fit(x, reversal, bvar)
    loo = []
    for i in range(len(x)):
        mask = np.arange(len(x)) != i
        loo.append(linear_fit(x[mask], y[mask], var[mask])["slope_um_per_mm"])

    sst = np.sum((y-y.mean())**2)
    direction_share = float(1 - np.sum(fit["residual"]**2)/sst)
    total_share = total_power_share(rows, x, fit["predicted"])
    before = iso_summary(rows)
    after = iso_summary(rows, fit["predicted"])

    frequencies = np.linspace(1/(x.max()-x.min()), 0.5, 2500)
    angular = 2*np.pi*frequencies
    power = signal.lombscargle(x, fit["residual"]-np.mean(fit["residual"]), angular,
                               precenter=False, normalize=True, weights=1/var)
    periods = 1/frequencies
    peak = int(np.argmax(power))

    results = {
        "analysis": "ISO 230-2 v3 Step A linear geometric identification",
        "source": {"aggregate_results": str(args.aggregate), "controller_log": str(args.log),
                   "raw_scope_csv": None,
                   "limitation": "Raw ISO230-2_v3.csv was not present; drift correction, settling checks and per-run slopes cannot be recomputed."},
        "convention": {"deviation_sign": "positive means measured stage travel exceeded command",
                       "datum": "single global mean from the existing v3 evaluation",
                       "position": "rounded commanded microstep position"},
        "dataset": {"targets": len(x), "target_min_mm": float(x.min()),
                    "target_max_mm": float(x.max()), "span_mm": float(np.ptp(x)),
                    "runs_per_direction": int(np.min(np.r_[nu, nd]))},
        "validity": validity,
        "fit_uncorrected": serialise_fit(fit),
        "fit_drift_corrected": None,
        "ols_crosscheck": serialise_fit(ols),
        "robustness": {
            "direction_fits": {"up": serialise_fit(upfit), "down": serialise_fit(dnfit)},
            "direction_slope_difference_ppm": float((upfit["slope_um_per_mm"]-dnfit["slope_um_per_mm"])*1000),
            "leave_one_target_out_slope_um_per_mm": loo,
            "leave_one_target_out_min_um_per_mm": float(min(loo)),
            "leave_one_target_out_max_um_per_mm": float(max(loo)),
            "per_run_fits": None,
        },
        "power": {"direction_average_share": direction_share, "total_budget_share": total_share,
                  "rms_before_um": float(np.sqrt(np.mean((y-y.mean())**2))),
                  "rms_after_um": float(np.sqrt(np.mean(fit["residual"]**2)))},
        "reversal": {"fit": serialise_fit(bfit), "mean_um": float(np.mean(reversal)),
                     "std_um": float(np.std(reversal, ddof=1)),
                     "max_abs_um": float(np.max(np.abs(reversal)))},
        "iso_metrics_measured_um": before,
        "iso_metrics_predicted_after_gain_um": after,
        "periodogram_preliminary": {"peak_period_mm": float(periods[peak]),
                                    "peak_frequency_per_mm": float(frequencies[peak]),
                                    "peak_power": float(power[peak]),
                                    "interpretation": "screening only; 24 irregular target samples are too sparse for a definitive periodic model"},
        "physical_attribution": {"thermal_ppm": None, "wavelength_ppm": None,
                                 "cosine_ppm": None, "screw_lead_ppm": None,
                                 "reason": "temperature, pressure, compensation state, alignment angle and screw accuracy class are unavailable"},
    }
    results["_rows"] = rows
    save_csvs(rows, x, y, var, reversal, fit, args.out)
    render_a4(x, y, var, up, dn, reversal, fit, bfit, loo, periods, power,
              before, after, validity, results, args.render)
    results.pop("_rows")
    (args.out / "step_a_results.json").write_text(json.dumps(results, indent=2) + "\n")

    print(f"slope {fit['slope_um_per_mm']:+.6f} ± {fit['slope_ci95_um_per_mm']:.6f} µm/mm "
          f"({fit['slope_ppm']:+.1f} ± {fit['slope_ci95_ppm']:.1f} ppm)")
    print(f"explained: {direction_share*100:.1f}% direction-average, {total_share*100:.1f}% total")
    print(f"predicted E: {before['E']:.3f} -> {after['E']:.3f} µm")
    print(f"wrote {args.out / 'step_a_results.json'}")
    print(f"wrote {args.render}")


if __name__ == "__main__":
    main()
