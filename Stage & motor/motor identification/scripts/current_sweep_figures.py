#!/usr/bin/env python3
"""Current sweep: separate torque-driven error from electrical error.

Actions 17-25 of the 2026-10-13 identification run. At each of three current
levels: a 4 s settle, then four electrical cycles (256 native 1/16 microsteps,
28.8 degrees of shaft) out and back as a staircase with 120 ms dwells.

Method, per the analysis brief:
  * segment and settle exactly as the static map (final 28% of each dwell);
  * apply the static map's geometry fit at absolute shaft angle - never refit
    it, because 28.8 degrees is not a whole revolution;
  * per level, mean of the two directions is the static error and half their
    difference is the hysteresis; fold both onto one electrical cycle;
  * fit 50/100/200 per-rev as complex amplitudes (cos and sin separately);
  * separate parts with A_k(I) = a_k + b_k / I and h(I) = h0 / I.

Figure 3a  static error over one electrical cycle at the three currents, each
           with its hysteresis band shaded.
Figure 3b  amplitudes against normalised inverse current I100 / I; intercept is
           the electrical part, slope the torque-driven part.

Shares its loading, settling and harmonic helpers with staircase_figures.py,
which sits beside it, and reuses that module's geometry fit.
"""
from pathlib import Path
import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from staircase_figures import (                                    # noqa: E402
    CACHE, ROOT, DPC, C2C, FULL_STEP_DEG, ELEC_U256, ELEC_DEG,
    REV_U256, TAIL, INK, INK2, INK3, SURFACE, GRID, BAND,
    design, harmonic_fit, settled, paired, load_static, style, titleblock,
    wrap,
)

OUT = Path(os.environ.get("CURRENT_OUT", ROOT / "analysis" / "current_sweep"))
WORK = OUT / "working"

# Axes box leaving room for a wrapped three-part title block above it.
AXES_RECT = [0.075, 0.105, 0.845, 0.595]


def titleblock3(fig, title, subtitle, stats, width=116):
    """Title block that wraps, so long lines cannot widen the saved canvas.

    savefig(bbox_inches="tight") grows the output to fit any text that
    overflows the figure, which silently stretches the aspect ratio.
    """
    import textwrap
    y = 0.975
    fig.text(0.055, y, title, size=13.5, weight="bold", color=INK, va="top")
    y -= 0.062
    for i, line in enumerate(textwrap.wrap(subtitle, width)):
        fig.text(0.055, y, line, size=9.6, color=INK2, va="top")
        y -= 0.034
    y -= 0.016
    for line in textwrap.wrap(stats, width):
        fig.text(0.055, y, line, size=9.6, color=INK, va="top")
        y -= 0.034

# Ordinal ramp, light to dark, for the three current levels (sequential blue,
# steps 250/450/650 - the ordinal floor on a light surface is step 250).
LEVEL_COLOR = {180: "#6da7ec", 270: "#256abf", 360: "#0d366b"}
# Categorical slots 1-3 for the three harmonic orders; the hysteresis is a
# different kind of quantity and wears ink, which needs no hue validation.
ORDER_COLOR = {50: "#2a78d6", 100: "#eb6834", 200: "#1baf7a"}
ORDER_MARKER = {50: "o", 100: "^", 200: "D"}

# (forward action, reverse action, preceding settle hold) per requested current.
LEVELS = {180: (18, 19, 17), 270: (21, 22, 20), 360: (24, 25, 23)}

# IRUN codes. CS=9 at 360 mA is confirmed by the run's own register readback;
# the other two follow from the same TMCStepper conversion. Only the absolute
# milliamp scale depends on that assumption - every ratio below depends solely
# on the integer CS codes.
CS_CODE = {180: 4, 270: 7, 360: 9}
VSENSE_FS_V = 0.180
R_TOTAL_OHM = 0.12


def realised_ma(requested):
    return ((CS_CODE[requested] + 1) / 32 * VSENSE_FS_V / R_TOTAL_OHM
            / np.sqrt(2) * 1000.0)


def load_current(path):
    a = np.loadtxt(path, delimiter=",", dtype=np.int64)
    return {"action": a[:, 0], "t": a[:, 1], "event": a[:, 2],
            "cmd": a[:, 3], "enc": a[:, 4], "ma": a[:, 5]}


def geometry_model():
    """Refit nothing: return the static map's fit so it can be evaluated."""
    d = load_static(CACHE / "static.csv")
    p = paired(d, 15, 16)
    ideal = p["cmd"] * C2C
    mean = ((p["ef"] - ideal) + (p["er"] - ideal)) / 2.0
    angle = p["cmd"] / REV_U256 * 360.0
    coef, names, _, geo = harmonic_fit(angle, mean, (1, 2), linear=False)

    def evaluate(angle_deg):
        A, _ = design(angle_deg, (1, 2), linear=False)
        return A @ coef

    return evaluate, geo, p, mean


def fold_exact(cmd, signal):
    """Fold one direction's own window, which is exactly four cycles.

    Each direction covers 256 commanded microsteps, so every one of the 64
    phases gets exactly 4 samples and the fold is perfectly balanced. Pairing
    the two directions first would leave 255 points and unbalance phase 0.
    """
    phase = np.mod(cmd, ELEC_U256).astype(int)
    keys = np.unique(phase)
    mean = np.array([signal[phase == k].mean() for k in keys])
    sem = np.array([signal[phase == k].std(ddof=1) / np.sqrt((phase == k).sum())
                    for k in keys])
    counts = np.array([(phase == k).sum() for k in keys])
    return keys, mean, sem, counts


def level_result(d, requested, geo_eval):
    fwd_a, rev_a, hold_a = LEVELS[requested]
    out = {}
    folded = {}
    for name, action in (("f", fwd_a), ("r", rev_a)):
        cmd, enc, std, n = settled(d, action)
        err = enc - cmd * C2C
        corrected = err - geo_eval(cmd / REV_U256 * 360.0)
        keys, m, sem, cnt = fold_exact(cmd, corrected)
        folded[name] = {"keys": keys, "mean": m, "sem": sem, "n": cnt}
    keys = folded["f"]["keys"]
    static = (folded["f"]["mean"] + folded["r"]["mean"]) / 2.0
    hyst = (folded["r"]["mean"] - folded["f"]["mean"]) / 2.0

    # The hold immediately before the pair gives the starting rotor position,
    # so the pair's return error is measured against it.
    hcmd, henc, hstd, hn = settled(d, hold_a)
    start_enc = float(henc[-1])
    rev_cmd, rev_enc, _, _ = settled(d, rev_a)
    end_enc = float(rev_enc[rev_cmd == 0][0])

    out.update({
        "requested_ma": requested, "realised_ma": realised_ma(requested),
        "cs": CS_CODE[requested],
        "keys": keys, "x_deg": keys / ELEC_U256 * ELEC_DEG,
        "static": static - static.mean(), "static_offset": float(static.mean()),
        "hyst": hyst, "sem": (folded["f"]["sem"] + folded["r"]["sem"]) / 2.0,
        "fold_n": folded["f"]["n"],
        "return_error": end_enc - start_enc,
        "hold_std": float(hstd[-1]),
    })
    return out


def complex_amplitudes(x_deg, signal):
    """cos and sin parts at electrical harmonics 1, 2, 4 = 50/100/200 per rev.

    Fitted on the folded cycle, which is exactly periodic, so the three orders
    are mutually orthogonal and cannot leak into one another.

    The standard error treats everything the three-harmonic model does not
    explain as noise. That is deliberately conservative: the residual contains
    real structure at other orders as well as measurement noise, so a term
    that clears this bar is genuinely resolved.
    """
    A, names = design(x_deg, (1, 2, 4), linear=False, period_deg=ELEC_DEG)
    coef, *_ = np.linalg.lstsq(A, signal, rcond=None)
    resid = signal - A @ coef
    dof = max(1, len(signal) - A.shape[1])
    sigma = float(np.sqrt(np.sum(resid ** 2) / dof))
    se = sigma * np.sqrt(2.0 / len(signal))   # per cos/sin coefficient
    out = {}
    for elec, mech in ((1, 50), (2, 100), (4, 200)):
        c = float(coef[names.index(f"cos{elec}")])
        s = float(coef[names.index(f"sin{elec}")])
        amp = float(np.hypot(c, s))
        out[mech] = {"cos": c, "sin": s, "amp": amp,
                     "phase_deg": float(np.degrees(np.arctan2(s, c)) % 360),
                     "se": float(se), "snr": float(amp / se) if se else np.inf,
                     "phase_se_deg": float(np.degrees(se / amp)) if amp else 180.0}
    out["_sigma"] = sigma
    return out


def fit_vs_inverse_current(x, y):
    """y = a + b*x with x = I100/I; returns intercept, slope and the residual.

    Three levels and two parameters leave one degree of freedom, so the
    residual is a real test of the model rather than a consequence of the fit.
    """
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    return {"intercept": float(coef[0]), "slope": float(coef[1]),
            "residual": [float(v) for v in (y - pred)],
            "max_abs_residual": float(np.max(np.abs(y - pred)))}


def fit_through_origin(x, y):
    """y = s*x, the friction model h = h0 / I. Two degrees of freedom."""
    s = float(np.sum(x * y) / np.sum(x * x))
    pred = s * x
    return {"slope": s, "residual": [float(v) for v in (y - pred)],
            "max_abs_residual": float(np.max(np.abs(y - pred)))}


# ------------------------------------------------------------------- figures

def figure3a(levels, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes(AXES_RECT)
    ax.grid(False)
    ax.grid(True, axis="y")
    for k in range(4):
        ax.axvline(k * FULL_STEP_DEG, color=INK3, lw=1.2, alpha=0.7, zorder=1)
    step_proxy = mpl.lines.Line2D([], [], color=INK3, lw=1.2, alpha=0.7,
                                  label="full-step positions")
    ax.axhline(0, color=INK3, lw=0.8, zorder=2)

    handles = []
    for req in (180, 270, 360):
        L = levels[req]
        col = LEVEL_COLOR[req]
        xs, ys = wrap(L["x_deg"], L["static"] * DPC)
        _, lo = wrap(L["x_deg"], (L["static"] - np.abs(L["hyst"])) * DPC)
        _, hi = wrap(L["x_deg"], (L["static"] + np.abs(L["hyst"])) * DPC)
        ax.fill_between(xs, lo, hi, color=col, alpha=0.15, lw=0, zorder=3)
        line, = ax.plot(xs, ys, color=col, lw=2.1, zorder=5,
                        label=f"{L['realised_ma']:.0f} mA (req. {req}, CS {L['cs']})")
        handles.append(line)
    band_proxy = mpl.patches.Patch(facecolor=INK3, alpha=0.20,
                                   label="hysteresis band")

    ax.set_xlim(0, ELEC_DEG)
    ax.set_xticks(np.arange(0, 7.21, 0.9))
    ax.margins(y=0.12)
    ax.set_xlabel("position within one electrical cycle (degrees of shaft "
                  "angle; 7.2° = 4 full steps)")
    ax.set_ylabel("static error, geometry removed (degrees)")
    ax.legend(handles=handles + [band_proxy, step_proxy], loc="lower right",
              bbox_to_anchor=(1.0, 1.01), ncol=5, handlelength=1.8,
              columnspacing=1.1, borderaxespad=0.0, fontsize=8.3)

    sec = ax.secondary_yaxis("right", functions=(lambda v: v / DPC,
                                                 lambda v: v * DPC))
    sec.set_ylabel("encoder counts", color=INK2)
    sec.spines["right"].set_color(GRID)

    pp = {r: np.ptp(levels[r]["static"]) * DPC for r in (180, 270, 360)}
    titleblock3(
        fig,
        "Figure 3a  Static error over one electrical cycle, at three currents",
        "Geometry fit from the static map applied at absolute shaft angle, then "
        "folded over the four electrical cycles of each level.",
        f"Peak-to-peak grows from {pp[360]:.3f}° at 331 mA to "
        f"{pp[180]:.3f}° at 166 mA ({pp[180] / pp[360]:.2f}×) and the "
        f"hysteresis band widens with it. Features that grow as current falls "
        f"are torque-driven; those that stay put are electrical.", width=124)
    return fig


def figure3b(series, metrics):
    fig = plt.figure(figsize=(11.0, 6.3))
    ax = fig.add_axes([AXES_RECT[0], AXES_RECT[1], 0.78, AXES_RECT[3]])
    x = np.array(metrics["x_inverse_current"])
    grid = np.linspace(0, 2.25, 50)

    ax.axvline(0, color=INK3, lw=1.0, zorder=1)
    ax.text(0.025, 0.97,
            "open marks at x = 0:\ncurrent-independent\n(electrical) part",
            transform=ax.transAxes, color=INK2, size=8.6, va="top")

    handles = []
    for mech in (50, 100, 200):
        s = series[f"order_{mech}"]
        y = np.array(s["amp_deg"])
        sep = metrics["separation"][f"order_{mech}"]
        col, mk = ORDER_COLOR[mech], ORDER_MARKER[mech]
        # Model the cos and sin parts separately, then take the magnitude.
        # When the electrical and torque parts sit at different phases the
        # magnitude is not linear in 1/I, so a straight line would misread
        # the intercept as the electrical part.
        solid = sep["resolved"]
        # An unresolved order must not be extrapolated back to x = 0: that
        # intercept is exactly the quantity the data cannot constrain. Draw
        # its curve only across the measured range.
        g = grid if solid else np.linspace(0.9, 2.1, 25)
        cm = np.hypot(sep["cos"]["intercept"] + sep["cos"]["slope"] * g,
                      sep["sin"]["intercept"] + sep["sin"]["slope"] * g)
        ax.plot(g, cm, color=col, lw=1.5, ls=(0, (6, 3)),
                alpha=1.0 if solid else 0.45, zorder=3)
        ax.errorbar(x, y, yerr=sep["se_deg"], fmt="none", ecolor=col,
                    elinewidth=1.3, capsize=3, zorder=4)
        h, = ax.plot(x, y, mk, color=col, ms=9, mec=SURFACE, mew=1.4, zorder=5,
                     label=f"{mech}/rev" + ("" if solid else " (unresolved)"))
        if solid:
            ax.plot([0], [sep["electrical_amp_deg"]], mk, color=SURFACE,
                    mec=col, mew=1.8, ms=8, zorder=5)
        handles.append(h)
        ax.annotate(f"{mech}/rev", xy=(x[0], y[0]),
                    xytext=(12, {50: -9, 100: 0, 200: 9}[mech]),
                    textcoords="offset points", color=col, size=9.2,
                    va="center", weight="bold",
                    alpha=1.0 if solid else 0.55)

    y = np.array(series["hysteresis"]["amp_deg"])
    fit = metrics["separation"]["hysteresis"]["through_origin"]
    ax.plot(grid, fit["slope"] * grid, color=INK, lw=1.5, ls=(0, (6, 3)),
            zorder=3)
    h, = ax.plot(x, y, "s", color=INK, ms=9, mec=SURFACE, mew=1.4, zorder=5,
                 label="hysteresis half-width")
    handles.append(h)
    ax.annotate("hysteresis", xy=(x[0], y[0]), xytext=(12, 0),
                textcoords="offset points", color=INK, size=9.2, va="center",
                weight="bold")

    ax.set_xlim(-0.09, 2.3)
    ax.set_ylim(bottom=0)
    ax.set_xticks([0, 0.5, 1.0, 1.25, 1.5, 2.0])
    ax.set_xlabel("normalised inverse current  $I_{100}/I$   "
                  "(1.00 = 331 mA, 1.25 = 265 mA, 2.00 = 166 mA)")
    ax.set_ylabel("amplitude (degrees)")
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.01),
              ncol=4, handlelength=1.3, columnspacing=1.5, borderaxespad=0.0,
              fontsize=8.8)

    hy = metrics["separation"]["hysteresis"]
    o100 = metrics["separation"]["order_100"]
    titleblock3(
        fig,
        "Figure 3b  Error amplitude against normalised inverse current",
        "Points measured (error bars 1$\\sigma$); dashed curves are "
        "$A(I) = a + b/I$ fitted to the cos and sin parts separately, then "
        "combined as vectors.",
        f"Hysteresis runs through the origin to within "
        f"{hy['through_origin']['max_abs_residual'] / DPC:.2f} counts, as friction "
        f"should. Only 100/rev is resolved well enough to split: "
        f"{100 - o100['electrical_share_at_nominal_pct']:.0f}% torque-driven at "
        f"nominal current.", width=124)
    return fig


def working_settling(d, path):
    """At 50% current omega_n is lower and friction relatively larger, so the
    120 ms dwell has to be re-checked rather than assumed from the static map."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    ring_hz = {}
    for ax, req in zip(axes, (180, 270, 360)):
        fwd_a = LEVELS[req][0]
        m = d["action"] == fwd_a
        cmd, enc, t = d["cmd"][m], d["enc"][m], d["t"][m]
        edges = np.r_[0, np.flatnonzero(np.diff(cmd) != 0) + 1, len(cmd)]
        grid = np.linspace(0, 120, 61)
        acc = []
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 10:
                continue
            acc.append(np.interp(grid, (t[a:b] - t[a]) / 1000.0,
                                 enc[a:b] - enc[b - 1],
                                 left=np.nan, right=np.nan))
        acc = np.asarray(acc)
        col = LEVEL_COLOR[req]
        ax.axhline(0, color=INK3, lw=0.8)
        ax.axvspan(120 * (1 - TAIL), 120, color=BAND, alpha=0.14, lw=0,
                   label="settled tail")
        ax.fill_between(grid, np.nanpercentile(acc, 5, axis=0),
                        np.nanpercentile(acc, 95, axis=0), color=col,
                        alpha=0.22, lw=0, label="5-95%")
        ax.plot(grid, np.nanmedian(acc, axis=0), color=col, lw=2.0,
                label="median")
        tail = acc[:, grid >= 120 * (1 - TAIL)]
        spread = np.nanpercentile(tail, 95) - np.nanpercentile(tail, 5)
        # Ring frequency of the ensemble median. Holding stiffness rises with
        # current, so omega_n should scale as sqrt(I); this measures whether
        # it does. Sampling is 500 Hz, so the FFT bin spacing is coarse.
        med = np.nanmedian(acc, axis=0)[:30]
        yy = (med - med[-10:].mean()) * np.hanning(len(med))
        ff = np.fft.rfftfreq(len(yy), d=float(np.diff(grid)[0]) / 1000.0)
        ring = float(ff[np.argmax(np.abs(np.fft.rfft(yy))[1:]) + 1])
        ring_hz[req] = ring
        ax.set_title(f"{realised_ma(req):.0f} mA  (requested {req})", size=10.5)
        ax.set_xlabel("time since step edge (ms)")
        ax.annotate(f"ring {ring:.0f} Hz\ntail spread (p95-p5) {spread:.2f} counts",
                    xy=(0.97, 0.07), xycoords="axes fraction", ha="right",
                    size=8.6, color=INK2)
        ax.legend(loc="upper right", fontsize=8.2)
    axes[0].set_ylabel("counts relative to dwell end")
    pred = {r: ring_hz[360] * np.sqrt(realised_ma(r) / realised_ma(360))
            for r in (180, 270)}
    fig.suptitle("Working plot: settling within each 120 ms dwell, by current "
                 "level", size=12, weight="bold", x=0.055, ha="left")
    fig.text(0.055, 0.905,
             "Ring frequency follows sqrt(I) as expected: "
             f"{ring_hz[180]:.0f} and {ring_hz[270]:.0f} Hz measured, "
             f"{pred[180]:.0f} and {pred[270]:.0f} Hz predicted from the 331 mA "
             "result. Every tail is settled, so the 120 ms dwell holds at 50% "
             "current too.",
             size=9.2, color=INK2, va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(path)
    plt.close(fig)
    return ring_hz


def working_crosscheck(cross, path):
    """Level differences cancel geometry, so they test the correction."""
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.grid(False)
    ax.grid(True, axis="y")
    ax.axhline(0, color=INK3, lw=0.8)
    for pair, col in ((("180", "360"), LEVEL_COLOR[180]),
                      (("270", "360"), LEVEL_COLOR[270])):
        c = cross[f"{pair[0]}_minus_{pair[1]}"]
        ax.plot(c["x_deg"], np.array(c["direct"]) * DPC, color=col, lw=2.2,
                label=f"{pair[0]}−{pair[1]} mA, raw (geometry cancels)")
        ax.plot(c["x_deg"], np.array(c["corrected"]) * DPC, color=INK, lw=1.0,
                ls=(0, (4, 3)),
                label=f"{pair[0]}−{pair[1]} mA, via geometry correction")
    ax.set_xlim(0, ELEC_DEG)
    ax.set_xlabel("position within one electrical cycle (degrees)")
    ax.set_ylabel("difference in static error (degrees)")
    ax.legend(fontsize=8.4, ncol=2, loc="upper center")
    ax.margins(y=0.3)
    ax.set_title("Working plot: level differences, with and without the "
                 "geometry correction", size=12, weight="bold", loc="left",
                 pad=12)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------- main

def main():
    src = CACHE / "current.csv"
    if not src.exists():
        raise SystemExit(
            f"missing {src}\n\n"
            "  awk -F, '$1==\"D\" && $4>=17 && $4<=25"
            "{print $4\",\"$3\",\"$5\",\"$6\",\"$7\",\"$9}' \\\n"
            "    data/2026-10-13_motor_identification/"
            "motor_identification_run.csv \\\n"
            "    > \"$STAIRCASE_CACHE/current.csv\"\n")
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    style()

    geo_eval, geo, static_pair, static_mean = geometry_model()
    d = load_current(src)
    levels = {req: level_result(d, req, geo_eval) for req in LEVELS}

    x = np.array([realised_ma(360) / levels[r]["realised_ma"]
                  for r in (180, 270, 360)])

    series = {f"order_{m}": {"cos_deg": [], "sin_deg": [], "amp_deg": [],
                             "phase_deg": []} for m in (50, 100, 200)}
    series["hysteresis"] = {"amp_deg": []}
    for req in (180, 270, 360):
        L = levels[req]
        amps = complex_amplitudes(L["x_deg"], L["static"])
        L["amplitudes"] = amps
        for mech in (50, 100, 200):
            series[f"order_{mech}"]["cos_deg"].append(amps[mech]["cos"] * DPC)
            series[f"order_{mech}"]["sin_deg"].append(amps[mech]["sin"] * DPC)
            series[f"order_{mech}"]["amp_deg"].append(amps[mech]["amp"] * DPC)
            series[f"order_{mech}"]["phase_deg"].append(amps[mech]["phase_deg"])
        series["hysteresis"]["amp_deg"].append(float(np.mean(L["hyst"])) * DPC)

    separation = {}
    for mech in (50, 100, 200):
        s = series[f"order_{mech}"]
        fc = fit_vs_inverse_current(x, np.array(s["cos_deg"]))
        fs = fit_vs_inverse_current(x, np.array(s["sin_deg"]))
        # The electrical and torque-driven contributions are vectors, not
        # magnitudes. Combining them as vectors is what keeps the separation
        # honest when the two parts sit at different phases: the magnitude of
        # the total is then not linear in 1/I at all.
        elec = np.hypot(fc["intercept"], fs["intercept"])
        torque = np.hypot(fc["slope"], fs["slope"])
        d_phase = (np.degrees(np.arctan2(fs["slope"], fc["slope"]))
                   - np.degrees(np.arctan2(fs["intercept"], fc["intercept"])))
        separation[f"order_{mech}"] = {
            "amp": fit_vs_inverse_current(x, np.array(s["amp_deg"])),
            "cos": fc, "sin": fs,
            "electrical_amp_deg": float(elec),
            "torque_amp_at_nominal_deg": float(torque),
            "phase_between_parts_deg": float((d_phase + 180) % 360 - 180),
            "electrical_share_at_nominal_pct": float(
                elec / np.hypot(fc["intercept"] + fc["slope"],
                                fs["intercept"] + fs["slope"]) * 100),
            "phase_deg": s["phase_deg"],
            "phase_swing_deg": float(np.ptp(
                np.degrees(np.unwrap(np.radians(s["phase_deg"]))))),
            # The fit runs in encoder counts; the figures and table are in
            # degrees, so the standard error has to be converted with them.
            "se_deg": [levels[r]["amplitudes"][mech]["se"] * DPC
                       for r in (180, 270, 360)],
            "snr": [levels[r]["amplitudes"][mech]["snr"]
                    for r in (180, 270, 360)],
            "phase_se_deg": [levels[r]["amplitudes"][mech]["phase_se_deg"]
                             for r in (180, 270, 360)],
            # A two-parameter split needs the order to be resolved at every
            # level. Four electrical cycles is a short window: only the
            # strongest order clears this bar, and the others' amplitudes are
            # not even monotonic in current, which is the giveaway.
            "resolved": bool(min(levels[r]["amplitudes"][mech]["snr"]
                                 for r in (180, 270, 360)) >= 8.0),
            "monotonic_in_current": bool(
                np.all(np.diff(np.array(s["amp_deg"])[::-1]) >= 0)),
        }
    separation["static_pp"] = fit_vs_inverse_current(
        x, np.array([np.ptp(levels[r]["static"]) * DPC
                     for r in (180, 270, 360)]))
    hy = np.array(series["hysteresis"]["amp_deg"])
    separation["hysteresis"] = {
        "free": fit_vs_inverse_current(x, hy),
        "through_origin": fit_through_origin(x, hy),
        "h_times_I_deg_mA": [float(hy[i] * levels[r]["realised_ma"])
                             for i, r in enumerate((180, 270, 360))],
    }

    # Cross-check: at identical commanded positions the geometry term is the
    # same, so a difference between two levels removes it without relying on
    # the fit at all.
    def direct_mean(req):
        fwd_a, rev_a, _ = LEVELS[req]
        cf, ef, _, _ = settled(d, fwd_a)
        cr, er, _, _ = settled(d, rev_a)
        cmd, i, j = np.intersect1d(cf, cr, return_indices=True)
        return cmd, ((ef[i] - cmd * C2C) + (er[j] - cmd * C2C)) / 2.0

    cross = {}
    for a, b in (("180", "360"), ("270", "360")):
        ca, ma_ = direct_mean(int(a))
        cb, mb = direct_mean(int(b))
        cmd, i, j = np.intersect1d(ca, cb, return_indices=True)
        keys, folded_direct, _, _ = fold_exact(cmd, ma_[i] - mb[j])
        corrected = levels[int(a)]["static"] - levels[int(b)]["static"]
        corrected = corrected - corrected.mean()
        fd = folded_direct - folded_direct.mean()
        cross[f"{a}_minus_{b}"] = {
            "x_deg": (keys / ELEC_U256 * ELEC_DEG).tolist(),
            "direct": fd.tolist(), "corrected": corrected.tolist(),
            "rms_disagreement_counts": float(np.sqrt(np.mean(
                (fd - corrected) ** 2))),
        }

    # Repeatability: the 100% level should reproduce the matching section of
    # the static map, at the same commanded positions, to about one count.
    cs_cmd, cs_mean = direct_mean(360)
    common, a_i, b_i = np.intersect1d(static_pair["cmd"], cs_cmd,
                                      return_indices=True)
    rep = static_mean[a_i] - cs_mean[b_i]
    offset = float(np.mean(rep))
    rep = rep - offset
    repeat = {"positions": int(len(common)),
              "rms_counts": float(np.sqrt(np.mean(rep ** 2))),
              "max_abs_counts": float(np.max(np.abs(rep))),
              "offset_counts": offset}

    metrics = {
        "source": {"run": "data/2026-10-13_motor_identification/"
                          "motor_identification_run.csv",
                   "actions": list(range(17, 26)), "dwell_ms": 120,
                   "settled_tail_fraction": TAIL,
                   "electrical_cycles_per_level": 4,
                   "shaft_span_deg": 4 * ELEC_DEG},
        "geometry_source": {
            "note": "fit taken from the static map (actions 15/16) and "
                    "evaluated at absolute shaft angle; not refitted here",
            "one_per_rev_counts": geo[1]["amp"],
            "one_per_rev_deg": geo[1]["amp"] * DPC,
            "two_per_rev_counts": geo[2]["amp"]},
        "x_inverse_current": x.tolist(),
        "levels": {
            str(r): {
                "requested_ma": r, "cs": levels[r]["cs"],
                "realised_ma": levels[r]["realised_ma"],
                "static_pp_deg": float(np.ptp(levels[r]["static"]) * DPC),
                "static_pp_counts": float(np.ptp(levels[r]["static"])),
                "hysteresis_half_deg": float(np.mean(levels[r]["hyst"]) * DPC),
                "hysteresis_half_counts": float(np.mean(levels[r]["hyst"])),
                "return_error_counts": levels[r]["return_error"],
                "fold_sem_median_counts": float(np.median(levels[r]["sem"])),
                "hold_std_counts": levels[r]["hold_std"],
                "amplitudes_deg": {
                    str(m): {"amp": levels[r]["amplitudes"][m]["amp"] * DPC,
                             "cos": levels[r]["amplitudes"][m]["cos"] * DPC,
                             "sin": levels[r]["amplitudes"][m]["sin"] * DPC,
                             "phase_deg": levels[r]["amplitudes"][m]["phase_deg"]}
                    for m in (50, 100, 200)}}
            for r in (180, 270, 360)},
        "separation": separation,
        "cross_check": {k: {"rms_disagreement_counts":
                            v["rms_disagreement_counts"]}
                        for k, v in cross.items()},
        "repeatability_vs_static_map": repeat,
    }

    f3a = figure3a(levels, metrics)
    f3a.savefig(OUT / "fig3a_static_error_by_current.png")
    plt.close(f3a)
    f3b = figure3b(series, metrics)
    f3b.savefig(OUT / "fig3b_amplitude_vs_inverse_current.png")
    plt.close(f3b)

    ring_hz = working_settling(d, WORK / "w4_settling_by_current.png")
    metrics["ring_hz"] = {str(k): v for k, v in ring_hz.items()}
    working_crosscheck(cross, WORK / "w5_level_difference_crosscheck.png")

    (OUT / "current_sweep_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    rows = np.column_stack(
        [levels[180]["x_deg"]]
        + [levels[r][k] * DPC for r in (180, 270, 360)
           for k in ("static", "hyst")])
    np.savetxt(OUT / "fig3a_curves.csv", rows, delimiter=",", comments="",
               header="shaft_deg_in_cycle,static_180_deg,hyst_180_deg,"
                      "static_270_deg,hyst_270_deg,static_360_deg,hyst_360_deg",
               fmt="%.6f")

    print(f"x = I100/I            {np.round(x, 4).tolist()}")
    for r in (180, 270, 360):
        L = metrics["levels"][str(r)]
        print(f"{r} mA (CS {L['cs']}, {L['realised_ma']:.0f} mA real): "
              f"pp {L['static_pp_deg']:.4f} deg, "
              f"hyst {L['hysteresis_half_deg']:.4f} deg, "
              f"return {L['return_error_counts']:+.0f} counts, "
              f"SEM {L['fold_sem_median_counts']:.3f}")
    for mech in (50, 100, 200):
        sp = separation[f"order_{mech}"]
        print(f"order {mech:<4}: electrical {sp['electrical_amp_deg']:.5f} deg, "
              f"torque@nominal {sp['torque_amp_at_nominal_deg']:.5f} deg, "
              f"elec share {sp['electrical_share_at_nominal_pct']:.0f}%, "
              f"parts {sp['phase_between_parts_deg']:+.0f} deg apart")
        print(f"            SNR {[round(v, 1) for v in sp['snr']]}, "
              f"phase se {[round(v, 1) for v in sp['phase_se_deg']]} deg, "
              f"phase swing {sp['phase_swing_deg']:.1f} deg")
    hyf = separation["hysteresis"]
    print(f"hysteresis free fit : intercept {hyf['free']['intercept']:.4f} deg, "
          f"slope {hyf['free']['slope']:.4f} deg")
    print(f"hysteresis origin   : slope {hyf['through_origin']['slope']:.4f} deg, "
          f"max resid {hyf['through_origin']['max_abs_residual'] / DPC:.2f} counts")
    print(f"h x I (deg.mA)      {[round(v, 4) for v in hyf['h_times_I_deg_mA']]}")
    for k, v in metrics["cross_check"].items():
        print(f"cross-check {k}: rms disagreement "
              f"{v['rms_disagreement_counts']:.3f} counts")
    print(f"repeatability vs static map: rms {repeat['rms_counts']:.2f} counts, "
          f"max {repeat['max_abs_counts']:.2f}, offset "
          f"{repeat['offset_counts']:+.2f}, over {repeat['positions']} positions")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
