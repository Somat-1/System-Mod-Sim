#!/usr/bin/env python3
"""Block 6 - velocity plateaus: friction against speed, and ripple magnification.

Actions 46-67 of the 2026-10-13 identification run. Eleven constant-rate
plateaus from 0.25 to 400 full-step/s, each immediately repeated in reverse
over the same angles, all at MRES=16 and sampled at 2 kHz. Every plateau starts
at its final rate with no ramp.

Two quantities come out of each plateau:

  * mean lag (theta_cmd - theta_enc) over a whole number of electrical cycles.
    The 50/100/200 static errors average to zero over whole cycles, so what is
    left is load torque over stiffness. Half the forward-reverse difference is
    friction; their mean should be about zero and is a check.
  * ripple, the error left after removing the mean lag, decomposed into the
    50/100/200 per-rev harmonics and divided by the static map's amplitudes to
    give the measured magnification.

Excitation frequencies at command speed v full-step/s: 50/rev at v/4, 100/rev
at v/2, 200/rev at v. Those are compared against the second-order magnification
|H(r)| = 1 / sqrt((1-r^2)^2 + (2 zeta r)^2) predicted by block 5, with no new
fitting.

Figure 5a  friction lag against speed, with the Stribeck + viscous fit.
Figure 5b  ripple magnification against frequency ratio, against |H(r)|.
"""
from pathlib import Path
import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from staircase_figures import (                                    # noqa: E402
    CACHE, ROOT, CPR, DPC, C2C, FULL_STEP_DEG, ELEC_U256, ELEC_DEG,
    REV_U256, INK, INK2, INK3, SURFACE, GRID, BAND, style, harmonic_fit,
)
from current_sweep_figures import (                                # noqa: E402
    geometry_model, titleblock3, AXES_RECT,
)

OUT = Path(os.environ.get("VELOCITY_OUT", ROOT / "analysis" / "velocity"))
WORK = OUT / "working"

# speed in full-step/s -> (forward action, reverse action)
PLATEAUS = {0.25: (46, 47), 0.5: (48, 49), 1: (50, 51), 2: (52, 53),
            5: (54, 55), 10: (56, 57), 20: (58, 59), 50: (60, 61),
            100: (62, 63), 200: (64, 65), 400: (66, 67)}

# Static map harmonic amplitudes, in counts (block 2).
STATIC_AMP_COUNTS = {50: 0.7526186425930725, 100: 1.328467058890878,
                     200: 0.5459225721629228}
F_N_HZ = 153.04              # block 5, settled full-step ensemble
ZETA = 0.0528                # block 5, from the decay envelope
HYST_DEG = 0.04418874        # block 2 static-map hysteresis half-width
# K_theta from the nameplate holding torque read as an RMS rating, scaled to
# the realised 331.5 mA. The peak reading would be sqrt(2) larger; see README.
K_THETA = 2.4862             # N m / rad
# A stepper cannot hold a lag beyond one full step: that is 90 deg electrical,
# the peak of the torque curve. Any larger apparent lag is lost position, not
# friction, and the plateau that produced it is unusable for the Stribeck curve.
MAX_SUSTAINABLE_LAG_DEG = FULL_STEP_DEG
DRIFT_LIMIT_COUNTS = 40      # lag drift across a plateau interior
ELEC_CYCLE_COUNTS = CPR / 50.0   # 400 counts = 7.2 deg = one electrical cycle

ORDER_COLOR = {50: "#2a78d6", 100: "#eb6834", 200: "#1baf7a"}
ORDER_MARKER = {50: "o", 100: "^", 200: "D"}
FWD = "#2a78d6"
REV = "#eb6834"


def load_velocity(path):
    a = np.loadtxt(path, delimiter=",", dtype=np.int64)
    return {"action": a[:, 0], "t": a[:, 1], "event": a[:, 2],
            "cmd": a[:, 3], "enc": a[:, 4]}


def plateau(d, action, geo_eval):
    """Trim one plateau to a whole number of electrical cycles of motion."""
    m = d["action"] == action
    t, cmd, enc = d["t"][m], d["cmd"][m], d["enc"][m]
    corrected = enc - geo_eval(enc / CPR * 360.0)
    lag = cmd * C2C - corrected          # counts; positive = rotor behind

    moving = np.flatnonzero(np.diff(cmd) != 0)
    if len(moving) < 4:
        return None
    i0, i1 = moving[0], moving[-1] + 1

    # Discard the first and last electrical cycle of commanded travel. The two
    # slowest plateaus span only two electrical cycles, so a full cycle at each
    # end would leave nothing; those trim a fifth of the span instead.
    c0, c1 = cmd[i0], cmd[i1]
    step = np.sign(c1 - c0)
    span_all = abs(c1 - c0)
    trim = ELEC_U256 if span_all > 3 * ELEC_U256 else span_all // 5
    lo, hi = c0 + step * trim, c1 - step * trim
    inner = (cmd - lo) * step >= 0
    inner &= (hi - cmd) * step >= 0
    inner[:i0] = False
    inner[i1:] = False
    if inner.sum() < 20:
        return None

    # Then keep a whole number of electrical cycles.
    ci = cmd[inner]
    span = abs(ci[-1] - ci[0])
    whole = int(span // ELEC_U256)
    if whole >= 1:
        end = ci[0] + step * whole * ELEC_U256
        inner &= (end - cmd) * step >= 0

    li = lag[inner]
    k = max(20, len(li) // 12)
    drift = float(li[-k:].mean() - li[:k].mean())
    return {"t": t[inner], "cmd": cmd[inner], "enc": enc[inner],
            "lag": li, "direction": int(step), "cycles": whole,
            "trim_u256": int(trim),
            "drift_counts": drift,
            "drifted": bool(abs(drift) > DRIFT_LIMIT_COUNTS)}


def ripple_harmonics(cmd, resid):
    """Fold the residual on the electrical cycle and fit 50/100/200 per rev."""
    phase = np.mod(cmd, ELEC_U256) / ELEC_U256 * ELEC_DEG
    bins = np.linspace(0, ELEC_DEG, 65)
    which = np.clip(np.digitize(phase, bins) - 1, 0, 63)
    folded = np.full(64, np.nan)
    for k in range(64):
        s = which == k
        if s.any():
            folded[k] = resid[s].mean()
    good = ~np.isnan(folded)
    if good.sum() < 16:
        return None
    centres = (bins[:-1] + bins[1:]) / 2
    _, names, _, out = harmonic_fit(centres[good], folded[good], (1, 2, 4),
                                    linear=False, period_deg=ELEC_DEG)
    return {50: out[1]["amp"], 100: out[2]["amp"], 200: out[4]["amp"],
            "folded": folded, "centres": centres,
            "coverage": int(good.sum())}


def magnification(r, zeta):
    return 1.0 / np.sqrt((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2)


def stribeck(v_deg_s, fc, fs, vs, sigma2):
    return fc + (fs - fc) * np.exp(-(v_deg_s / vs) ** 2) + sigma2 * v_deg_s


def fit_stribeck(v, f):
    def resid(p):
        return stribeck(v, *p) - f

    p0 = [np.median(f), f.max(), max(v.min(), 1.0), 0.0]
    r = least_squares(resid, p0,
                      bounds=([0, 0, 1e-3, 0], [1, 1, 1e4, 1e-2]),
                      xtol=1e-12, ftol=1e-12)
    return {"F_c_deg": float(r.x[0]), "F_s_deg": float(r.x[1]),
            "v_s_deg_s": float(r.x[2]), "sigma2_deg_per_deg_s": float(r.x[3]),
            "rms_residual_deg": float(np.sqrt(np.mean(r.fun ** 2)))}


# ------------------------------------------------------------------- figures

def figure5a(res, fit, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes([AXES_RECT[0], AXES_RECT[1], 0.78, AXES_RECT[3] - 0.06])

    speeds = np.array(sorted(res))
    v_deg = speeds * FULL_STEP_DEG
    fric = np.array([abs(res[s]["friction_deg"]) for s in speeds])
    lost = np.array([res[s]["step_loss"] for s in speeds])
    neg = np.array([res[s]["negative_friction"] and not res[s]["step_loss"]
                    for s in speeds])
    good = ~(lost | neg)

    v_step = metrics["stepwise_below_fs_s"] * FULL_STEP_DEG
    ax.axvspan(v_deg.min() * 0.6, v_step, color=BAND, alpha=0.18, lw=0,
               zorder=1, label="stepwise motion (microstep rate < $f_n$)")
    ax.axhspan(MAX_SUSTAINABLE_LAG_DEG, 1e3, color="#b83a3a", alpha=0.10,
               lw=0, zorder=1)
    ax.axhline(MAX_SUSTAINABLE_LAG_DEG, color="#b83a3a", lw=1.6,
               ls=(0, (5, 3)), zorder=3,
               label=f"one full step ({MAX_SUSTAINABLE_LAG_DEG}°) = "
                     "90° electrical: no larger lag can be held")
    ax.axhline(HYST_DEG, color=INK3, lw=1.6, ls=(0, (5, 3)), zorder=3,
               label=f"static-map hysteresis {HYST_DEG:.3f}°")

    ax.plot(v_deg[good], fric[good], "o", color=INK, ms=9, mec=SURFACE,
            mew=1.3, zorder=6, label="usable plateau")
    ax.plot(v_deg[neg], fric[neg], "s", color=SURFACE, mec="#eda100", mew=2.0,
            ms=9, zorder=6, label="lag has the wrong sign - not usable")
    ax.plot(v_deg[lost], fric[lost], "o", color=SURFACE, mec="#b83a3a",
            mew=2.0, ms=9, zorder=6, label="position lost - not usable")
    for x, y, b in zip(v_deg, fric, lost):
        ax.plot([x, x], [max(y, 1e-4), 1e-4], color=INK3, lw=0.6, alpha=0.4,
                zorder=2)

    if fit is not None:
        g = np.logspace(np.log10(v_deg.min()), np.log10(v_deg.max()), 200)
        ax.plot(g, stribeck(g, fit["F_c_deg"], fit["F_s_deg"],
                            fit["v_s_deg_s"], fit["sigma2_deg_per_deg_s"]),
                color=INK2, lw=1.6, ls=(0, (6, 3)), zorder=5,
                label="Stribeck + viscous fit")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rotor speed (degrees/s)")
    ax.set_ylabel("apparent friction lag (degrees)")
    ax.legend(loc="lower right", fontsize=7.8, ncol=1, framealpha=0.0)
    ax.set_ylim(0.01, 12)

    top = ax.secondary_xaxis("top", functions=(lambda v: v / FULL_STEP_DEG,
                                               lambda v: v * FULL_STEP_DEG))
    top.set_xlabel("command speed (full-step/s)", color=INK2)
    top.spines["top"].set_color(GRID)
    right = ax.secondary_yaxis(
        "right", functions=(lambda a: np.radians(a) * K_THETA * 1000,
                            lambda T: np.degrees(T / 1000 / K_THETA)))
    right.set_ylabel("implied torque (mN·m)", color=INK2)
    right.spines["right"].set_color(GRID)

    titleblock3(
        fig,
        "Figure 5a  Apparent friction lag against speed",
        "Half the forward-reverse lag difference per plateau, over a whole "
        "number of electrical cycles. Torque axis uses "
        f"K = {K_THETA:.2f} N·m/rad from block 5.",
        metrics["headline_5a"], width=124)
    return fig


def figure5b(points, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes([AXES_RECT[0], AXES_RECT[1], 0.80, AXES_RECT[3] - 0.04])

    r_grid = np.logspace(-2.6, 0.75, 400)
    ax.plot(r_grid, magnification(r_grid, ZETA), color=INK, lw=2.2,
            zorder=5, label=f"|H(r)| from block 5  ($\\zeta$ = {ZETA:.3f})")
    ax.axvline(1.0, color=INK3, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(1.0, 0.985, " resonance", transform=ax.get_xaxis_transform(),
            color=INK2, size=8.6, va="top")
    ax.axhline(1.0, color=INK3, lw=0.9, zorder=2)

    for order in (50, 100, 200):
        p = points[order]
        if not len(p["r"]):
            continue
        ok = ~np.asarray(p["step_loss"], dtype=bool)
        ax.plot(np.asarray(p["r"])[ok], np.asarray(p["mag"])[ok],
                ORDER_MARKER[order], color=ORDER_COLOR[order], ms=8,
                mec=SURFACE, mew=1.2, zorder=6, label=f"{order}/rev")
        if (~ok).any():
            ax.plot(np.asarray(p["r"])[~ok], np.asarray(p["mag"])[~ok],
                    ORDER_MARKER[order], color=SURFACE,
                    mec=ORDER_COLOR[order], mew=1.6, ms=8, zorder=6)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("frequency ratio  $r = f_{excitation} / f_n$")
    ax.set_ylabel("ripple amplitude / static-map amplitude")
    ax.legend(loc="lower left", fontsize=8.4, framealpha=0.0)

    titleblock3(
        fig,
        "Figure 5b  Ripple magnification against frequency ratio",
        "Each harmonic of each plateau, divided by its static-map amplitude. "
        "The line is the second-order magnification predicted by block 5, with "
        "nothing refitted. Hollow marks are plateaus with step loss.",
        metrics["headline_5b"], width=124)
    return fig


def working_traces(samples, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (speed, s) in zip(axes, samples.items()):
        ang = (s["cmd"] - s["cmd"][0]) / REV_U256 * 360.0
        ax.plot(ang, (s["lag"] - s["lag"].mean()) * DPC, color=FWD, lw=0.8)
        ax.set_title(f"{speed:g} full-step/s"
                     + (f"  - drift {s['drift_counts']:+.0f} counts"
                        if s["drifted"] else "  - clean"), size=10.5)
        ax.set_xlabel("shaft angle travelled (deg)")
        ax.set_ylabel("error about the mean (deg)")
    fig.suptitle("Working plot: raw error against angle for three plateaus "
                 "(transients and step loss)", size=12, weight="bold",
                 x=0.055, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------- main

def main():
    src = CACHE / "velocity.csv"
    if not src.exists():
        raise SystemExit(
            f"missing {src}\n\n"
            "  awk -F, '$1==\"D\" && $4>=46 && $4<=67"
            "{print $4\",\"$3\",\"$5\",\"$6\",\"$7}' \\\n"
            "    data/2026-10-13_motor_identification/"
            "motor_identification_run.csv \\\n"
            "    > \"$STAIRCASE_CACHE/velocity.csv\"\n")
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    style()

    geo_eval, _, _, _ = geometry_model()
    d = load_velocity(src)

    res, points = {}, {o: {"r": [], "mag": [], "speed": [], "step_loss": []}
                       for o in (50, 100, 200)}
    samples = {}
    for speed, (fa, ra) in PLATEAUS.items():
        pf, pr = plateau(d, fa, geo_eval), plateau(d, ra, geo_eval)
        if pf is None or pr is None:
            continue
        lag_f, lag_r = float(pf["lag"].mean()), float(pr["lag"].mean())
        entry = {
            "speed_fs_s": speed, "speed_deg_s": speed * FULL_STEP_DEG,
            "lag_fwd_deg": lag_f * DPC, "lag_rev_deg": lag_r * DPC,
            "friction_deg": (lag_f - lag_r) / 2 * DPC,
            "offset_deg": (lag_f + lag_r) / 2 * DPC,
            "cycles_fwd": pf["cycles"], "cycles_rev": pr["cycles"],
            "drift_fwd_counts": pf["drift_counts"],
            "drift_rev_counts": pr["drift_counts"],
            "drifted": bool(pf["drifted"] or pr["drifted"]),
            "ripple_deg": {},
        }
        entry["exceeds_physical_lag"] = bool(
            abs(entry["friction_deg"]) > MAX_SUSTAINABLE_LAG_DEG)
        # Friction always opposes motion, so the rotor lags the direction of
        # travel and this difference must be positive. A negative value is
        # proof that the plateau pair carries an offset larger than the
        # friction it is meant to measure.
        entry["negative_friction"] = bool(entry["friction_deg"] <= 0)
        entry["step_loss"] = bool(entry["drifted"]
                                  or entry["exceeds_physical_lag"])
        entry["stepwise"] = bool(speed < F_N_HZ / 16.0)
        entry["valid_for_friction"] = not (entry["step_loss"]
                                           or entry["negative_friction"])
        entry["valid_for_ripple"] = not (entry["step_loss"] or entry["stepwise"])
        for p in (pf, pr):
            h = ripple_harmonics(p["cmd"], p["lag"] - p["lag"].mean())
            if h is None:
                continue
            for order in (50, 100, 200):
                entry["ripple_deg"].setdefault(str(order), []).append(
                    h[order] * DPC)
        for order in (50, 100, 200):
            vals = entry["ripple_deg"].get(str(order))
            if not vals:
                continue
            amp = float(np.mean(vals))
            entry["ripple_deg"][str(order)] = amp
            f_exc = order * speed / 200.0
            points[order]["r"].append(f_exc / F_N_HZ)
            points[order]["mag"].append(amp / (STATIC_AMP_COUNTS[order] * DPC))
            points[order]["speed"].append(speed)
            points[order]["step_loss"].append(not entry["valid_for_ripple"])
        res[speed] = entry
        if speed in (1, 50, 400):
            samples[speed] = pf

    speeds = np.array(sorted(res))
    v_deg = speeds * FULL_STEP_DEG
    fric = np.array([res[s]["friction_deg"] for s in speeds])
    ok = np.array([res[s]["valid_for_friction"] for s in speeds])
    # Only fit the Stribeck curve if enough plateaus survive both tests;
    # with four free parameters, fewer than five points is not a fit.
    # Four free parameters need more than a handful of points, and the
    # surviving set is both small and does not span the Stribeck knee.
    fit = fit_stribeck(v_deg[ok], fric[ok]) if ok.sum() >= 6 else None

    stepwise_below = F_N_HZ / 16.0        # microstep rate = 16 v
    peak_order = max(
        ((o, points[o]["speed"][int(np.argmax(points[o]["mag"]))],
          max(points[o]["mag"])) for o in (50, 100, 200) if points[o]["mag"]),
        key=lambda x: x[2])

    metrics = {
        "source": {"run": "data/2026-10-13_motor_identification/"
                          "motor_identification_run.csv",
                   "actions": list(range(46, 68)),
                   "trim": "first and last electrical cycle discarded, then a "
                           "whole number of electrical cycles kept"},
        "carried_over": {"f_n_hz": F_N_HZ, "zeta": ZETA,
                         "K_theta_N_m_rad": K_THETA,
                         "static_hysteresis_deg": HYST_DEG,
                         "static_amplitudes_counts": STATIC_AMP_COUNTS},
        "stepwise_below_fs_s": stepwise_below,
        "plateaus": {str(s): res[s] for s in speeds},
        "stribeck_fit_deg": fit,
        "stribeck_fit_torque_mNm": None if fit is None else {
            "F_c": float(np.radians(fit["F_c_deg"]) * K_THETA * 1000),
            "F_s": float(np.radians(fit["F_s_deg"]) * K_THETA * 1000),
            "v_s_deg_s": fit["v_s_deg_s"],
            "sigma2_mNm_per_deg_s": float(
                np.radians(fit["sigma2_deg_per_deg_s"]) * K_THETA * 1000)},
        "max_sustainable_lag_deg": MAX_SUSTAINABLE_LAG_DEG,
        "drift_limit_counts": DRIFT_LIMIT_COUNTS,
        "usable_speeds_fs_s": [float(x) for x in speeds
                               if res[x]["valid_for_friction"]],
        "step_loss_speeds_fs_s": [float(x) for x in speeds
                                  if res[x]["step_loss"]],
        "max_offset_deg": float(np.max(np.abs(
            [res[s]["offset_deg"] for s in speeds]))),
        "peak_magnification": {"order": int(peak_order[0]),
                               "speed_fs_s": float(peak_order[1]),
                               "value": float(peak_order[2])},
    }
    n_ok = int(ok.sum())
    metrics["headline_5a"] = (
        f"Only {n_ok} of {len(speeds)} plateaus survive: the rest either drift "
        f"by more than {DRIFT_LIMIT_COUNTS} counts across the plateau or show "
        f"an apparent lag beyond one full step, which no stepper can hold. The "
        f"Stribeck curve therefore cannot be fitted from this run - the "
        f"velocity block needs repeating with closed-loop position checks.")
    solid = [(o, r, m) for o in (50, 100, 200)
             for r, m, bad in zip(points[o]["r"], points[o]["mag"],
                                  points[o]["step_loss"]) if not bad]
    metrics["usable_ripple_points"] = len(solid)
    metrics["headline_5b"] = (
        f"{len(solid)} of {sum(len(points[o]['r']) for o in (50, 100, 200))} "
        f"points survive. Below {stepwise_below:.0f} full-step/s the rotor "
        f"moves in discrete steps rather than following the command, so its "
        f"ripple is the stepping itself and not the static error being "
        f"magnified; those plateaus are shown hollow with the step-loss ones.")

    f5a = figure5a(res, fit, metrics)
    f5a.savefig(OUT / "fig5a_friction_vs_speed.png")
    plt.close(f5a)
    f5b = figure5b(points, metrics)
    f5b.savefig(OUT / "fig5b_ripple_magnification.png")
    plt.close(f5b)
    if samples:
        working_traces(samples, WORK / "w8_plateau_traces.png")

    (OUT / "velocity_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    hdr = "speed_fs_s,speed_deg_s,lag_fwd_deg,lag_rev_deg,friction_deg," \
          "offset_deg,ripple50_deg,ripple100_deg,ripple200_deg,step_loss"
    rows = [[s, res[s]["speed_deg_s"], res[s]["lag_fwd_deg"],
             res[s]["lag_rev_deg"], res[s]["friction_deg"], res[s]["offset_deg"],
             res[s]["ripple_deg"].get("50", np.nan),
             res[s]["ripple_deg"].get("100", np.nan),
             res[s]["ripple_deg"].get("200", np.nan),
             int(res[s]["step_loss"])] for s in speeds]
    np.savetxt(OUT / "fig5_plateaus.csv", np.asarray(rows, dtype=float),
               delimiter=",", comments="", header=hdr, fmt="%.6f")

    print(f"{'v fs/s':>8} {'fwd deg':>9} {'rev deg':>9} {'fric deg':>9} "
          f"{'offset':>8} {'r50':>8} {'r100':>8} {'r200':>8} {'cyc':>5} loss")
    for s in speeds:
        e = res[s]
        rp = e["ripple_deg"]
        print(f"{s:8g} {e['lag_fwd_deg']:9.4f} {e['lag_rev_deg']:9.4f} "
              f"{e['friction_deg']:9.4f} {e['offset_deg']:8.4f} "
              f"{rp.get('50', float('nan')):8.4f} "
              f"{rp.get('100', float('nan')):8.4f} "
              f"{rp.get('200', float('nan')):8.4f} "
              f"{e['cycles_fwd']:5d} {'YES' if e['step_loss'] else ''}")
    print(f"\nstepwise below      {stepwise_below:.1f} fs/s")
    print(f"usable for friction {metrics['usable_speeds_fs_s']}")
    print(f"rejected            {metrics['step_loss_speeds_fs_s']} (drift or "
          f"lag beyond one full step)")
    print(f"negative friction   "
          f"{[float(x) for x in speeds if res[x]['negative_friction']]}")
    if fit is None:
        print("Stribeck            NOT FITTED - too few usable plateaus")
    else:
        t = metrics["stribeck_fit_torque_mNm"]
        print(f"Stribeck (torque)   F_c {t['F_c']:.3f} mN.m, "
              f"F_s {t['F_s']:.3f} mN.m")
    print(f"max |offset|        {metrics['max_offset_deg']:.4f} deg "
          f"(should be ~0)")
    print(f"step loss at        {metrics['step_loss_speeds_fs_s']} fs/s")
    print(f"peak magnification  {peak_order[2]:.2f}x on {peak_order[0]}/rev "
          f"at {peak_order[1]:g} fs/s")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
