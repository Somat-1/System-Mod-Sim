#!/usr/bin/env python3
"""Block 5 - dynamics: natural frequency, damping type, and absolute torque.

Actions 26-45 of the 2026-10-13 identification run. At each of four rotor
angles (0, 50, 100, 150 full steps) the firmware issues twelve out-and-back
single-step excitations at MRES=1 (one full step, 1.8 deg) and twelve at
MRES=16 (one microstep, 0.1125 deg), in both directions, each followed by a
350 ms dwell. Every step edge carries event=+-1, which is t = 0 for its
response. Both the excitation and the reset edge are usable single-step
responses, so each action yields 24.

THE FULL-STEP RESPONSES ARE MOSTLY UNUSABLE, and that is the block's main
result. A full step is 90 deg electrical and the measured overshoot is 78%,
which puts the rotor's peak at about 160 deg electrical - just short of the
180 deg pull-out angle. Only the first twelve responses settle one step away;
every later one slips into the adjacent equilibrium and travels roughly two
full steps for a one-step command. See `classify_full_step`.

Model, in full-step units (x = 1 is one full step = 90 deg electrical):

    x'' = -(wn^2 / (pi/2)) sin((pi/2)(x - x_cmd)) - 2 zeta wn x' - f_c sgn(x')

The sine is the magnetic restoring torque; for small (x - x_cmd) it reduces to
-wn^2 (x - x_cmd), so wn is the small-signal natural frequency. f_c = wn^2 *
delta, with delta the Coulomb dead band in full steps. Fitted on the clean
full-step ensemble only, then used to *predict* the microstep responses, which
are validation and are never fitted.

Figure 4a  full-step responses, normalised, clean and slipped separated.
Figure 4b  microstep responses in counts, with the model prediction.
"""
from pathlib import Path
import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.signal import find_peaks

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from staircase_figures import (                                    # noqa: E402
    CACHE, ROOT, CPR, DPC, FULL_STEP_DEG, INK, INK2, INK3, SURFACE, GRID,
    BAND, style,
)
from current_sweep_figures import (                                # noqa: E402
    geometry_model, titleblock3, AXES_RECT,
)

OUT = Path(os.environ.get("DYNAMICS_OUT", ROOT / "analysis" / "dynamics"))
WORK = OUT / "working"

N_R = 50                       # rotor teeth / electrical cycles per revolution
I_M_DATASHEET = 9.0e-7         # kg m^2, motor rotor alone (repo model params)
T_D_DATASHEET = 0.005          # N m, detent torque amplitude (spec sheet)
T_HOLD_NAMEPLATE = 0.060       # N m at the 400 mA nameplate current
NAMEPLATE_MA = 400.0
REALISED_MA_RMS = 331.5        # CS = 9, confirmed by register readback
HYST_DEG = 0.04418874          # half-width from the static map (block 2)

COUNTS_PER_FULL_STEP = CPR / 200.0                    # 100
COUNTS_PER_MICROSTEP = COUNTS_PER_FULL_STEP / 16.0    # 6.25
SETTLE_FROM_MS = 250.0
CLEAN_TOL = 0.15               # |final - 1| below this counts as settled

# (full-step +, full-step -, microstep +, microstep -) per rotor angle.
ANGLE_ACTIONS = {0: (26, 27, 28, 29), 50: (31, 32, 33, 34),
                 100: (36, 37, 38, 39), 150: (41, 42, 43, 44)}
ANGLE_COLOR = {0: "#6da7ec", 50: "#2a78d6", 100: "#1c5cab", 150: "#0d366b"}

FWD = "#2a78d6"
ACCENT = "#eb6834"
SLIP = "#eda100"


def load_dynamics(path):
    a = np.loadtxt(path, delimiter=",", dtype=np.int64)
    return {"action": a[:, 0], "t": a[:, 1], "event": a[:, 2],
            "cmd": a[:, 3], "enc": a[:, 4], "mres": a[:, 5]}


def extract(d, action, target_counts, geo_eval, grid_ms):
    """One normalised response per marked edge, resampled onto a common grid.

    The geometry correction is applied at the shaft angle the *encoder*
    reports, not the commanded angle: the dynamics block's absolute command
    bookkeeping is disturbed by its own MRES changes, while the encoder is not.
    """
    m = d["action"] == action
    t, enc, ev = d["t"][m], d["enc"][m], d["event"][m]
    corrected = enc - geo_eval(enc / CPR * 360.0)
    out, sign = [], []
    for i in np.flatnonzero(ev != 0):
        rel = (t - t[i]) / 1000.0
        sel = (rel >= grid_ms[0] - 2) & (rel <= grid_ms[-1] + 2)
        if sel.sum() < 50 or rel[sel].max() < grid_ms[-1]:
            continue
        base = corrected[sel][rel[sel] < -0.25]
        if len(base) < 2:
            continue
        out.append(np.interp(grid_ms,
                             rel[sel],
                             (corrected[sel] - base.mean()) * float(ev[i])
                             / target_counts))
        sign.append(int(ev[i]))
    return np.asarray(out), np.asarray(sign)


def finals(grid_ms, y):
    return y[:, grid_ms >= SETTLE_FROM_MS].mean(axis=1)


def classify_full_step(grid_ms, y):
    """Split responses that settled one step away from those that slipped."""
    f = finals(grid_ms, y)
    clean = np.abs(f - 1.0) < CLEAN_TOL
    return clean, f


def ringdown(grid_ms, y):
    """Extrema, half-period frequencies and the settled position."""
    final = float(np.mean(y[grid_ms >= SETTLE_FROM_MS]))
    dev = y - final
    hi, _ = find_peaks(dev, prominence=0.03)
    lo, _ = find_peaks(-dev, prominence=0.03)
    ext = np.sort(np.r_[hi, lo])
    ext = ext[grid_ms[ext] > 0]
    out = {"final": final, "peak_t": grid_ms[ext], "peak_a": np.abs(dev[ext]),
           "overshoot": float(dev[grid_ms > 0].max())}
    if len(ext) >= 3:
        half = np.diff(grid_ms[ext]) / 1000.0
        out["f_half"] = 1.0 / (2 * half)
        out["f_amp"] = (out["peak_a"][:-1] + out["peak_a"][1:]) / 2
    else:
        out["f_half"] = np.array([])
        out["f_amp"] = np.array([])
    return out


def first_peak_fn(grid_ms, y):
    """Small-signal f_n and zeta from the first overshoot of a step response.

    For a second-order step response the first peak sits at t = pi / w_d and
    its height gives zeta through the standard overshoot relation. This is the
    only route to f_n that the microstep responses support: they are too
    coarsely quantised for a multi-cycle ringdown fit, but their first peak is
    resolved.
    """
    final = float(np.mean(y[grid_ms >= SETTLE_FROM_MS]))
    w = (grid_ms > 0) & (grid_ms < 40)
    g, v = grid_ms[w], y[w]
    # The FIRST local maximum, not the global one: these traces are quantised
    # to one encoder count and a later ringing bump can exceed the first
    # overshoot, which would read as a much lower natural frequency.
    sm = np.convolve(v, np.ones(9) / 9, mode="same")
    pk_idx, _ = find_peaks(sm)
    if len(pk_idx) == 0:
        return None
    k = int(pk_idx[0])
    t_pk, pk = g[k], sm[k]
    if final <= 0 or pk <= final:
        return None
    os_ratio = (pk - final) / final
    ln = np.log(max(os_ratio, 1e-6))
    zeta = float(-ln / np.sqrt(np.pi ** 2 + ln ** 2))
    w_d = np.pi / (t_pk / 1000.0)
    w_n = w_d / max(np.sqrt(1 - zeta ** 2), 1e-3)
    return {"t_peak_ms": float(t_pk), "overshoot_ratio": float(os_ratio),
            "zeta": zeta, "f_n_hz": float(w_n / (2 * np.pi)),
            "final": final}


def model_response(t_s, wn, zeta, delta, x_cmd=1.0, eps=1e-4):
    """Integrate the nonlinear single-step model in full-step units."""
    k = wn ** 2 / (np.pi / 2)
    fc = wn ** 2 * delta

    def rhs(_t, s):
        x, v = s
        return [v,
                -k * np.sin((np.pi / 2) * (x - x_cmd))
                - 2 * zeta * wn * v - fc * np.tanh(v / eps)]

    sol = solve_ivp(rhs, (0.0, t_s[-1]), [0.0, 0.0], t_eval=t_s,
                    method="LSODA", rtol=1e-8, atol=1e-10, max_step=2e-4)
    return sol.y[0]


def fit_greybox(grid_ms, mean_y, p0=(150.0, 0.05, 0.0246), t_max_ms=200.0):
    """Fit f_n, zeta and the Coulomb band to the measured ensemble.

    The window has to run well past the first few cycles. Over 0-80 ms a large
    Coulomb band and a near-zero viscous term reproduce the early decay just as
    well as the reverse, and the two only separate once the late, low-amplitude
    ringing is included: Coulomb decay is linear and would have stopped the
    rotor by then, viscous decay is exponential and does not.
    """
    t_s = grid_ms / 1000.0
    use = (grid_ms >= 0) & (grid_ms <= t_max_ms)

    def resid(p):
        return model_response(t_s[use], p[0] * 2 * np.pi, p[1], p[2]) - mean_y[use]

    r = least_squares(resid, p0, bounds=([50, 0.0, 0.0], [400, 0.6, 0.2]),
                      xtol=1e-12, ftol=1e-12)
    return {"f_n_hz": float(r.x[0]), "zeta": float(r.x[1]),
            "delta_full_steps": float(r.x[2]),
            "delta_deg": float(r.x[2] * FULL_STEP_DEG),
            "rms_residual": float(np.sqrt(np.mean(r.fun ** 2))),
            "fit_window_ms": t_max_ms}


def fit_greybox_fixed_zeta(grid_ms, mean_y, zeta, p0=(150.0, 0.012),
                           t_max_ms=200.0):
    """Fit f_n and the Coulomb band with zeta taken from the decay envelope.

    The free three-parameter fit is poorly conditioned: viscous and Coulomb
    damping trade off against each other and it settles in a Coulomb-heavy
    minimum that cannot reproduce the late ringing. Taking zeta from the
    envelope's decay ratio, which is a model-free measurement, leaves a
    better-conditioned two-parameter problem and a lower residual.
    """
    t_s = grid_ms / 1000.0
    use = (grid_ms >= 0) & (grid_ms <= t_max_ms)

    def resid(p):
        return model_response(t_s[use], p[0] * 2 * np.pi, zeta, p[1]) - mean_y[use]

    r = least_squares(resid, p0, bounds=([50, 0.0], [400, 0.2]),
                      xtol=1e-12, ftol=1e-12)
    return {"f_n_hz": float(r.x[0]), "zeta": float(zeta),
            "zeta_source": "decay envelope",
            "delta_full_steps": float(r.x[1]),
            "delta_deg": float(r.x[1] * FULL_STEP_DEG),
            "rms_residual": float(np.sqrt(np.mean(r.fun ** 2))),
            "fit_window_ms": t_max_ms}


def envelope_type(cycle, amp):
    """A linear envelope means Coulomb damping, log-linear means viscous."""
    good = amp > 0
    c, a = np.asarray(cycle)[good], np.asarray(amp)[good]
    lin = np.polyfit(c, a, 1)
    r2_lin = 1 - np.var(a - np.polyval(lin, c)) / np.var(a)
    ex = np.polyfit(c, np.log(a), 1)
    pred = np.exp(np.polyval(ex, c))
    r2_exp = 1 - np.var(a - pred) / np.var(a)
    return {"r2_linear": float(r2_lin), "r2_exponential": float(r2_exp),
            "linear_slope_per_half_cycle": float(lin[0]),
            "log_slope_per_half_cycle": float(ex[0]),
            "decay_ratio_per_half_cycle": float(np.exp(ex[0])),
            # amplitude ratio per half cycle = exp(-pi zeta / sqrt(1-zeta^2))
            "zeta_from_envelope": float(-ex[0] / np.pi),
            "verdict": "viscous (exponential)" if r2_exp > r2_lin
                       else "Coulomb (linear)"}


# ------------------------------------------------------------------- figures

def figure4a(grid_ms, clean_y, slip_y, mean_clean, fit, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes(AXES_RECT)
    ax.grid(False)
    ax.grid(True, axis="y")

    band = HYST_DEG / FULL_STEP_DEG
    ax.axhspan(1 - band, 1 + band, color=BAND, alpha=0.25, lw=0, zorder=2)
    ax.axhline(1.0, color=INK3, lw=0.9, zorder=3)
    ax.axhline(2.0, color=SLIP, lw=0.9, ls=(0, (4, 3)), zorder=3)

    for y in slip_y[::4]:
        ax.plot(grid_ms, y, color=SLIP, lw=0.5, alpha=0.16, zorder=1)
    for y in clean_y:
        ax.plot(grid_ms, y, color=FWD, lw=0.7, alpha=0.40, zorder=4)
    ax.plot(grid_ms, mean_clean, color=FWD, lw=1.8, zorder=6)
    t_s = grid_ms[grid_ms >= 0] / 1000.0
    pred = model_response(t_s, fit["f_n_hz"] * 2 * np.pi, fit["zeta"],
                          fit["delta_full_steps"])
    ax.plot(grid_ms[grid_ms >= 0], pred, color=INK, lw=2.4, zorder=7)

    h = [mpl.lines.Line2D([], [], color=FWD, lw=0.9,
                          label=f"settled at one step (n={len(clean_y)})"),
         mpl.lines.Line2D([], [], color=SLIP, lw=0.9,
                          label=f"slipped (n={len(slip_y)})"),
         mpl.lines.Line2D([], [], color=FWD, lw=1.8, label="mean of settled"),
         mpl.lines.Line2D([], [], color=INK, lw=2.4, label="fitted model"),
         mpl.patches.Patch(facecolor=BAND, alpha=0.25,
                           label=f"±{HYST_DEG:.3f}° friction band")]
    ax.legend(handles=h, loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=5,
              handlelength=1.6, columnspacing=1.0, borderaxespad=0.0,
              fontsize=8.0)

    ax.set_xlim(-4, 90)
    ax.set_ylim(-0.45, 2.75)
    ax.set_xlabel("time from step edge (ms)")
    ax.set_ylabel("position  (0 = start, 1 = one full step)")
    sec = ax.secondary_yaxis("right", functions=(lambda v: v * FULL_STEP_DEG,
                                                 lambda v: v / FULL_STEP_DEG))
    sec.set_ylabel("degrees", color=INK2)
    sec.spines["right"].set_color(GRID)

    env = metrics["envelope"]
    titleblock3(
        fig,
        "Figure 4a  Full-step responses",
        f"All {metrics['full_step']['n_total']} single-full-step responses, "
        "normalised and sign-flipped to overlay. A full step is 90° "
        "electrical, so this is a large-signal test, not a small-signal one.",
        f"Only {metrics['full_step']['n_clean']} settle one step away: "
        f"overshoot {metrics['full_step']['overshoot_pct']:.0f}% puts the peak "
        f"near the 180° pull-out angle, and the other "
        f"{metrics['full_step']['n_slipped']} slip into the next equilibrium. "
        f"Fitted on the settled ones: f_n = {fit['f_n_hz']:.1f} Hz, "
        f"zeta = {fit['zeta']:.3f}. Envelope is {env['verdict']}.", width=124)
    return fig


def figure4b(grid_ms, micro, mean_micro, fit, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes(AXES_RECT)
    ax.grid(False)
    ax.grid(True, axis="y")

    band_counts = HYST_DEG / DPC
    ax.axhspan(COUNTS_PER_MICROSTEP - band_counts,
               COUNTS_PER_MICROSTEP + band_counts,
               color=BAND, alpha=0.25, lw=0, zorder=2)
    ax.axhline(COUNTS_PER_MICROSTEP, color=INK3, lw=0.9, zorder=3)

    for y in micro[::3]:
        ax.plot(grid_ms, y * COUNTS_PER_MICROSTEP, color=ACCENT, lw=0.5,
                alpha=0.13, zorder=1)
    ax.plot(grid_ms, mean_micro * COUNTS_PER_MICROSTEP, color=ACCENT, lw=2.0,
            zorder=5)
    t_s = grid_ms[grid_ms >= 0] / 1000.0
    pred = model_response(t_s, fit["f_n_hz"] * 2 * np.pi, fit["zeta"],
                          fit["delta_full_steps"], x_cmd=1 / 16)
    ax.plot(grid_ms[grid_ms >= 0], pred * COUNTS_PER_FULL_STEP, color=INK,
            lw=2.4, zorder=6)

    h = [mpl.lines.Line2D([], [], color=ACCENT, lw=0.9, alpha=0.5,
                          label=f"{len(micro)} responses"),
         mpl.lines.Line2D([], [], color=ACCENT, lw=2.0, label="ensemble mean"),
         mpl.lines.Line2D([], [], color=INK, lw=2.4,
                          label="model prediction (not fitted here)"),
         mpl.patches.Patch(facecolor=BAND, alpha=0.25,
                           label=f"±{band_counts:.1f} count friction band")]
    ax.legend(handles=h, loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=4,
              handlelength=1.7, columnspacing=1.2, borderaxespad=0.0,
              fontsize=8.4)

    ax.set_xlim(-4, 50)
    ax.set_xlabel("time from step edge (ms)")
    ax.set_ylabel("encoder counts  (commanded = 6.25)")
    sec = ax.secondary_yaxis("right", functions=(lambda v: v * DPC,
                                                 lambda v: v / DPC))
    sec.set_ylabel("degrees", color=INK2)
    sec.spines["right"].set_color(GRID)

    mm = metrics["microstep"]
    titleblock3(
        fig,
        "Figure 4b  Microstep responses",
        f"All {len(micro)} single-microstep responses in encoder counts, "
        "against the model carried over unchanged from the full-step fit.",
        f"One small overshoot then a stop, as a friction-dominated step should "
        f"give: first-peak zeta = {mm['zeta_first_peak']:.2f} against "
        f"{metrics['fit']['zeta']:.3f} for the full step. Settled travel is "
        f"{mm['final_counts_mean']:.2f} counts against 6.25 commanded "
        f"({mm['final_pct']:.0f}%), varying {mm['ratio_min']:.2f}-"
        f"{mm['ratio_max']:.2f} with rotor angle - the detent modulating "
        f"stiffness.", width=124)
    return fig


def working_envelope(metrics, cyc, amp, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    env = metrics["envelope"]
    for ax, logy in zip(axes, (False, True)):
        ax.plot(cyc, amp, "o", color=FWD, ms=8, mec=SURFACE, mew=1.3,
                label="measured peak amplitude")
        g = np.linspace(cyc.min(), cyc.max(), 50)
        ax.plot(g, env["linear_slope_per_half_cycle"] * (g - cyc[0]) + amp[0],
                color=INK3, lw=1.6, ls=(0, (5, 3)), label="linear (Coulomb)")
        ax.plot(g, amp[0] * np.exp(env["log_slope_per_half_cycle"]
                                   * (g - cyc[0])),
                color=ACCENT, lw=1.6, ls=(0, (5, 3)),
                label="exponential (viscous)")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel("extremum number")
        ax.set_ylabel("amplitude about the settled value (full steps)")
        ax.set_title("log scale" if logy else "linear scale", size=10.5)
        ax.legend(fontsize=8.4)
    fig.suptitle(f"Working plot: decay envelope - {env['verdict']} "
                 f"(R2 linear {env['r2_linear']:.3f}, "
                 f"exponential {env['r2_exponential']:.3f})",
                 size=12, weight="bold", x=0.055, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path)
    plt.close(fig)


def working_slip(grid_ms, per_action, path):
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.grid(False)
    ax.grid(True, axis="y")
    ax.axhline(1.0, color=INK3, lw=1.0)
    ax.axhspan(1 - CLEAN_TOL, 1 + CLEAN_TOL, color=BAND, alpha=0.22, lw=0,
               label="counted as settled")
    ax.axhline(2.0, color=SLIP, lw=1.0, ls=(0, (4, 3)),
               label="adjacent equilibrium")
    for (angle, action), f in per_action.items():
        ax.plot(np.arange(len(f)), f, "o-", ms=4, lw=1.0,
                color=ANGLE_COLOR[angle], alpha=0.85,
                label=f"{angle} fs (action {action})")
    ax.set_xlabel("response number within the action")
    ax.set_ylabel("settled position (full steps)")
    ax.legend(fontsize=7.8, ncol=5, loc="upper center")
    ax.margins(y=0.22)
    ax.set_title("Working plot: where each full-step response settled - "
                 "only the start of action 26 stays at one step",
                 size=12, weight="bold", loc="left", pad=12)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------- main

def main():
    src = CACHE / "dynamics.csv"
    if not src.exists():
        raise SystemExit(
            f"missing {src}\n\n"
            "  awk -F, '$1==\"D\" && $4>=26 && $4<=45"
            "{print $4\",\"$3\",\"$5\",\"$6\",\"$7\",\"$10}' \\\n"
            "    data/2026-10-13_motor_identification/"
            "motor_identification_run.csv \\\n"
            "    > \"$STAIRCASE_CACHE/dynamics.csv\"\n")
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    style()

    geo_eval, _, _, _ = geometry_model()
    d = load_dynamics(src)
    grid_ms = np.arange(-4.0, 340.0, 0.25)

    full, slip, micro, per_action, micro_ratio, micro_fn = [], [], [], {}, {}, {}
    for angle, (fp, fm, mp, mm) in ANGLE_ACTIONS.items():
        for a in (fp, fm):
            y, _ = extract(d, a, COUNTS_PER_FULL_STEP, geo_eval, grid_ms)
            clean, f = classify_full_step(grid_ms, y)
            per_action[(angle, a)] = f
            full.append(y[clean])
            slip.append(y[~clean])
        for a in (mp, mm):
            y, _ = extract(d, a, COUNTS_PER_MICROSTEP, geo_eval, grid_ms)
            micro.append(y)
            micro_ratio[(angle, a)] = float(np.mean(finals(grid_ms, y)))
            fpk = first_peak_fn(grid_ms, y.mean(axis=0))
            if fpk:
                micro_fn[(angle, a)] = fpk
    clean_y = np.vstack(full)
    slip_y = np.vstack(slip)
    micro = np.vstack(micro)

    mean_clean = clean_y.mean(axis=0)
    mean_micro = micro.mean(axis=0)
    ens = ringdown(grid_ms, mean_clean)
    cyc = np.arange(1, len(ens["peak_a"]) + 1, dtype=float)
    env = envelope_type(cyc, ens["peak_a"])
    free_fit = fit_greybox(grid_ms, mean_clean)
    fit = fit_greybox_fixed_zeta(grid_ms, mean_clean, env["zeta_from_envelope"])

    fn_by_angle = {}
    for angle in ANGLE_ACTIONS:
        vals = [v["f_n_hz"] for (a, _), v in micro_fn.items() if a == angle]
        if vals:
            fn_by_angle[angle] = float(np.mean(vals))
    fn_arr = np.array(list(fn_by_angle.values()))
    fn_spread = float((fn_arr.max() - fn_arr.min()) / fn_arr.mean() * 100) \
        if len(fn_arr) > 1 else float("nan")

    finals_clean = finals(grid_ms, clean_y)
    micro_fin = finals(grid_ms, micro) * COUNTS_PER_MICROSTEP
    micro_pk = first_peak_fn(grid_ms, mean_micro)
    ratios = np.array(list(micro_ratio.values()))

    wn = fit["f_n_hz"] * 2 * np.pi
    hyst_rad = np.radians(HYST_DEG)
    k_rotor = I_M_DATASHEET * wn ** 2
    t_hold_peak = T_HOLD_NAMEPLATE * (REALISED_MA_RMS * np.sqrt(2)) / NAMEPLATE_MA
    t_hold_rms = T_HOLD_NAMEPLATE * REALISED_MA_RMS / NAMEPLATE_MA
    scale = {
        "f_n_hz": fit["f_n_hz"], "omega_n_rad_s": float(wn),
        "I_m_datasheet_kg_m2": I_M_DATASHEET,
        "K_theta_if_J_is_rotor_only": float(k_rotor),
        "T_hold_if_J_is_rotor_only": float(k_rotor / N_R),
        "T_f_if_J_is_rotor_only": float(k_rotor * hyst_rad),
        "T_hold_nameplate_scaled_rms_reading": float(t_hold_rms),
        "T_hold_nameplate_scaled_peak_reading": float(t_hold_peak),
        "K_theta_nameplate_rms_reading": float(N_R * t_hold_rms),
        "K_theta_nameplate_peak_reading": float(N_R * t_hold_peak),
        "J_implied_rms_reading": float(N_R * t_hold_rms / wn ** 2),
        "J_implied_peak_reading": float(N_R * t_hold_peak / wn ** 2),
        "T_f_nameplate_rms_reading": float(N_R * t_hold_rms * hyst_rad),
        "T_f_nameplate_peak_reading": float(N_R * t_hold_peak * hyst_rad),
        "friction_pct_of_holding_torque": float(np.sin(hyst_rad * N_R) * 100),
        "detent_datasheet_N_m": T_D_DATASHEET,
    }

    metrics = {
        "source": {"run": "data/2026-10-13_motor_identification/"
                          "motor_identification_run.csv",
                   "actions": list(range(26, 46)),
                   "sample_interval_us_median": float(np.median(np.diff(
                       d["t"][d["action"] == 26]))),
                   "dwell_ms": 350, "settle_window_from_ms": SETTLE_FROM_MS},
        "full_step": {
            "n_total": int(len(clean_y) + len(slip_y)),
            "n_clean": int(len(clean_y)), "n_slipped": int(len(slip_y)),
            "clean_tolerance_full_steps": CLEAN_TOL,
            "overshoot_pct": float(ens["overshoot"] * 100),
            "settled_mean_full_steps": float(np.mean(finals_clean)),
            "settled_scatter_deg": float(np.std(finals_clean) * FULL_STEP_DEG),
            "settled_offset_deg": float((np.mean(finals_clean) - 1)
                                        * FULL_STEP_DEG),
            "slipped_mean_full_steps": float(np.mean(finals(grid_ms, slip_y))),
            "peak_times_ms": ens["peak_t"].tolist(),
            "peak_amplitudes_full_steps": ens["peak_a"].tolist(),
            "half_period_frequencies_hz": ens["f_half"].tolist()},
        "fit": fit,
        "free_fit_unconstrained": free_fit,
        "envelope": env,
        "microstep": {
            "n": int(len(micro)),
            "final_counts_mean": float(np.mean(micro_fin)),
            "final_counts_std": float(np.std(micro_fin)),
            "final_pct": float(np.mean(micro_fin) / COUNTS_PER_MICROSTEP * 100),
            "ratio_min": float(ratios.min()), "ratio_max": float(ratios.max()),
            "ratio_by_action": {f"{k[0]}fs_a{k[1]}": v
                                for k, v in micro_ratio.items()},
            "zeta_first_peak": micro_pk["zeta"] if micro_pk else None,
            "f_n_first_peak_hz": micro_pk["f_n_hz"] if micro_pk else None,
            "t_peak_ms": micro_pk["t_peak_ms"] if micro_pk else None},
        "f_n_by_angle_hz_from_microstep": fn_by_angle,
        "f_n_spread_pct": fn_spread,
        "static_map_hysteresis_deg": HYST_DEG,
        "absolute_scale": scale,
    }

    f4a = figure4a(grid_ms, clean_y, slip_y, mean_clean, fit, metrics)
    f4a.savefig(OUT / "fig4a_full_step_responses.png")
    plt.close(f4a)
    f4b = figure4b(grid_ms, micro, mean_micro, fit, metrics)
    f4b.savefig(OUT / "fig4b_microstep_responses.png")
    plt.close(f4b)
    working_envelope(metrics, cyc, ens["peak_a"], WORK / "w6_decay_envelope.png")
    working_slip(grid_ms, per_action, WORK / "w7_settled_position.png")

    (OUT / "dynamics_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    np.savetxt(OUT / "fig4_ensembles.csv",
               np.column_stack([grid_ms, mean_clean,
                                mean_micro * COUNTS_PER_MICROSTEP]),
               delimiter=",", comments="",
               header="t_ms,full_step_normalised,microstep_counts", fmt="%.6f")

    fs = metrics["full_step"]
    print(f"full-step responses  {fs['n_total']} total, {fs['n_clean']} settled, "
          f"{fs['n_slipped']} slipped (mean {fs['slipped_mean_full_steps']:.2f} steps)")
    print(f"overshoot            {fs['overshoot_pct']:.1f}%  "
          f"-> peak at {(1 + ens['overshoot']) * 90:.0f} deg electrical")
    print(f"grey-box fit         f_n {fit['f_n_hz']:.2f} Hz, zeta {fit['zeta']:.4f} "
          f"(from envelope), delta {fit['delta_deg']:.4f} deg "
          f"(rms resid {fit['rms_residual']:.4f})")
    print(f"  free 3-param fit   f_n {free_fit['f_n_hz']:.2f} Hz, "
          f"zeta {free_fit['zeta']:.4f}, delta {free_fit['delta_deg']:.4f} deg "
          f"(rms {free_fit['rms_residual']:.4f}) - ill-conditioned, not used")
    print(f"half-period f        {np.round(ens['f_half'], 1).tolist()} Hz")
    print(f"envelope             {env['verdict']} (R2 lin {env['r2_linear']:.3f}, "
          f"exp {env['r2_exponential']:.3f}, ratio {env['decay_ratio_per_half_cycle']:.3f})")
    print(f"settled              {fs['settled_mean_full_steps']:.4f} steps, "
          f"offset {fs['settled_offset_deg']:+.4f} deg, "
          f"scatter {fs['settled_scatter_deg']:.4f} deg "
          f"(hysteresis {HYST_DEG:.4f})")
    mm = metrics["microstep"]
    print(f"microstep            {mm['final_counts_mean']:.2f} counts "
          f"({mm['final_pct']:.0f}%), ratio {mm['ratio_min']:.2f}-{mm['ratio_max']:.2f}, "
          f"zeta {mm['zeta_first_peak']:.2f}, f_n {mm['f_n_first_peak_hz']:.0f} Hz")
    print(f"f_n by angle         "
          f"{ {k: round(v, 1) for k, v in fn_by_angle.items()} } spread {fn_spread:.1f}%")
    print(f"K (J = rotor only)   {scale['K_theta_if_J_is_rotor_only']:.4f} N.m/rad "
          f"-> T_hold {scale['T_hold_if_J_is_rotor_only'] * 1000:.1f} mN.m")
    print(f"nameplate T_hold     {t_hold_rms * 1000:.1f} (RMS) / "
          f"{t_hold_peak * 1000:.1f} (peak) mN.m")
    print(f"J implied            {scale['J_implied_rms_reading']:.2e} / "
          f"{scale['J_implied_peak_reading']:.2e} kg.m2 "
          f"(rotor alone {I_M_DATASHEET:.2e})")
    print(f"T_f implied          {scale['T_f_nameplate_rms_reading'] * 1000:.2f} / "
          f"{scale['T_f_nameplate_peak_reading'] * 1000:.2f} mN.m")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
