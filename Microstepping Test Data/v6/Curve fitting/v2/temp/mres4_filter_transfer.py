"""Apply the tail-tuned v2 filter chain to real MRES 4 commanded motion.

The v2 parameters (median w=4, Gaussian W=10 / sigma=W/6, cubic B-spline at
~750-sample knots) were derived from the *idle* tail, where the only content is
noise. This script applies exactly those operators to MRES 4 stepping and ramp
blocks -- where the structure the simulation model has to match actually lives --
and overlays v1's purpose-built "dynamic" fit for comparison.

The point is the frequency collision. At MRES 4 the microstep is 2.5 um, so the
commanded step rate is:

    slow (27.5 fsps,  275 um/s)  ->  9.09 ms/step  ->  110 Hz
    mod  (70.0 fsps,  700 um/s)  ->  3.57 ms/step  ->  280 Hz
    fast (200.0 fsps, 2000 um/s) ->  1.25 ms/step  ->  800 Hz

against a filter whose -3 dB point is 112 Hz and a capture whose Nyquist is
500 Hz. Filters are applied over each whole block (as on the tail), then a
window is zoomed for display.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.interpolate import make_lsq_spline

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent.parent / "analysis" / "scripts"))
from ids_common import apply_style, C_MEASURED, C_GRID

NPZ = HERE.parent.parent / "v1" / "MRES_4" / "fitted_blocks.npz"

# v2 tail parameters, applied verbatim.
MEDIAN_W = 4
GAUSS_W = 10
GAUSS_SIGMA = 10 / 6
SPLINE_KNOTS = 750
FS = 1000.0

FULL_STEP_UM, MRES = 10.0, 4
STEP_UM = FULL_STEP_UM / MRES
RATE_FSPS = {"slow": 27.5, "mod": 70.0, "fast": 200.0}

C_DESPIKE, C_SMOOTH, C_SPLINE, C_V1 = "#8a4fd6", "#d94801", "#0f9d58", "#111111"


def despike(y):
    return median_filter(y, size=MEDIAN_W, mode="nearest")


def smooth(y):
    return gaussian_filter1d(despike(y), sigma=GAUSS_SIGMA, truncate=3.0, mode="nearest")


def bspline(t, y, spacing=SPLINE_KNOTS, k=3):
    idx = np.arange(spacing, t.size - spacing, spacing)
    if idx.size < 1:
        return None
    knots = np.r_[(t[0],) * (k + 1), t[idx], (t[-1],) * (k + 1)]
    return make_lsq_spline(t, y, knots, k=k)(t)


def style(ax):
    ax.grid(True, color=C_GRID, alpha=.7, linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.set_xlabel("Time [ms from window start]", fontsize=9)
    ax.set_ylabel("Position [µm]", fontsize=9)


def main():
    d = np.load(NPZ)

    # --- pick the snippets -------------------------------------------------
    snippets = []

    # 1. Stepping: one commanded step transition from the oscillation block,
    #    centred on a seed time so the overshoot/ring-down is in frame.
    ot, om = d["oscillation_time_s"], d["oscillation_measured_um"]
    seeds = d["oscillation_seed_times_s"]
    seed = seeds[len(seeds) // 2]
    snippets.append(("stepping (oscillation, one commanded step)",
                     ot, om, d["oscillation_dynamic_um"],
                     seed - 0.05, seed + 0.25, None))

    # 2-4. Ramps: a mid-ramp window on each rate, sized to show ~8 microsteps
    #      (or 40 ms when the steps are too fast to resolve).
    for label in ("D-slow", "D-mod", "D-fast"):
        lt = d[f"leg_{label}_time_s"]
        lm = d[f"leg_{label}_measured_um"]
        ld = d[f"leg_{label}_dynamic_um"]
        rate = RATE_FSPS[label.split("-")[1]]
        step_ms = 1e3 * STEP_UM / (rate * FULL_STEP_UM)
        width_s = max(8 * step_ms, 40.0) / 1e3
        # Centre on the ramp itself, not the leg's temporal midpoint -- the
        # stored legs carry dwell on either side, and the midpoint can easily
        # land in one of them.
        speed = np.abs(np.gradient(gaussian_filter1d(lm, sigma=20.0, mode="nearest"), lt))
        moving = np.flatnonzero(speed > 0.5 * speed.max())
        mid = float(lt[moving[moving.size // 2]])
        snippets.append((f"ramp {label}  ({rate:g} fsps, {step_ms:.2f} ms/step, "
                         f"{1e3/step_ms:.0f} Hz)",
                         lt, lm, ld, mid - width_s / 2, mid + width_s / 2, step_ms))

    # --- filter whole blocks, then zoom ------------------------------------
    apply_style()
    fig, axes = plt.subplots(len(snippets), 1, figsize=(13, 4.6 * len(snippets)))

    for ax, (title, t, meas, v1dyn, lo, hi, step_ms) in zip(axes, snippets):
        des = despike(meas)
        smo = smooth(meas)
        spl = bspline(t, meas)

        m = (t >= lo) & (t <= hi)
        tt = (t[m] - t[m][0]) * 1e3          # ms from window start
        base = meas[m][0]

        ax.plot(tt, meas[m] - base, color=C_MEASURED, lw=1.0, marker=".", ms=3,
                alpha=.55, label="raw", zorder=1)
        ax.plot(tt, des[m] - base, color=C_DESPIKE, lw=1.3, alpha=.9,
                label=f"median despike (w={MEDIAN_W})", zorder=2)
        ax.plot(tt, smo[m] - base, color=C_SMOOTH, lw=2.0, alpha=.95,
                label=f"+ Gaussian (W={GAUSS_W}, $\\sigma$={GAUSS_SIGMA:.2f})  "
                      f"[tail params]", zorder=4)
        if spl is not None:
            ax.plot(tt, spl[m] - base, color=C_SPLINE, lw=2.0, alpha=.95,
                    label=f"cubic B-spline ({SPLINE_KNOTS} smp knots)", zorder=3)
        ax.plot(tt, v1dyn[m] - base, color=C_V1, lw=1.4, ls="--", alpha=.85,
                label="v1 dynamic fit (per-leg adaptive)", zorder=5)

        # Mark the commanded step cadence so the smearing is legible.
        if step_ms is not None and step_ms >= 2.0:
            for k in range(1, int(tt[-1] / step_ms) + 1):
                ax.axvline(k * step_ms, color="#999999", lw=.6, ls=":", alpha=.55)

        style(ax)
        ax.set_xlim(tt[0], tt[-1])
        # Scale to the measurement, not to whatever the spline does -- at these
        # knot spacings the spline can leave the window entirely, and that
        # shouldn't squash the traces being compared.
        yy = meas[m] - base
        pad = 0.15 * max(np.ptp(yy), 1e-3)
        ax.set_ylim(yy.min() - pad, yy.max() + pad)
        ax.set_title(title, loc="left", fontsize=9.5)
        ax.legend(loc="best", fontsize=7.5, framealpha=.92)

    fig.tight_layout()
    out = HERE / "mres4_filter_transfer.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

    # --- how much step structure survives each filter ----------------------
    print("\nresidual step structure (std of first difference, nm) over the "
          "whole block:")
    print(f"{'block':12s} {'raw':>9s} {'despike':>9s} {'+gauss':>9s} {'v1 dyn':>9s}")
    for title, t, meas, v1dyn, lo, hi, step_ms in snippets:
        name = title.split()[0] + ("" if step_ms is None else " " + title.split()[1])
        row = [np.diff(x).std(ddof=1) * 1e3 for x in
               (meas, despike(meas), smooth(meas), v1dyn)]
        print(f"{name:12s} " + " ".join(f"{v:9.1f}" for v in row))


if __name__ == "__main__":
    main()
