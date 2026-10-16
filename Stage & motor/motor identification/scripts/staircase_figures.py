#!/usr/bin/env python3
"""Two presentation figures and one metrics table from the 2026-10-13 staircase.

The staircase is actions 15 and 16 of the identification run: every native 1/16
microstep over one shaft revolution forward, then the identical commanded
positions in reverse, with a 150 ms settling dwell at each.

Figure 1  theta_enc - theta_cmd over one revolution, both directions, with the
          fitted 1/rev + 2/rev geometry term overlaid.
Figure 2  the same error minus that fit, folded onto one electrical cycle and
          averaged over all 50 cycles; forward and reverse with the gap shaded.

Working plots (settling, spectrum, microstep size) go to working/ and are not
part of the presentation set.

Colours are slots 1 and 2 of the validated categorical palette, used unchanged.
"""
from pathlib import Path
import json
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
MAIN = ROOT / "data" / "2026-10-13_motor_identification"

# Both are overridable so the figures can be rendered somewhere with free space
# when the repository volume is full.
CACHE = Path(os.environ.get("STAIRCASE_CACHE",
                            Path(__file__).resolve().parent / "_staircase_cache"))
OUT = Path(os.environ.get("STAIRCASE_OUT", ROOT / "analysis" / "staircase"))
WORK = OUT / "working"

CPR = 20000.0                 # encoder counts per shaft revolution
DPC = 360.0 / CPR             # degrees per count = 0.018
APC = DPC * 3600.0            # arcsec per count = 64.8
C2C = CPR / (200 * 256)       # counts per 1/256 full step = 0.390625
FULL_STEP_DEG = 360.0 / 200   # 1.8
ELEC_U256 = 4 * 256           # one electrical cycle = 4 full steps = 1024 u256
ELEC_DEG = 4 * FULL_STEP_DEG  # 7.2
REV_U256 = 200 * 256          # 51200
TAIL = 0.28                   # fraction of each dwell treated as settled

INK = "#0b0b0b"
INK2 = "#52514e"
INK3 = "#8a8880"
SURFACE = "#fcfcfb"
GRID = "#e4e3dd"
FWD = "#2a78d6"               # categorical slot 1
REV = "#eb6834"               # categorical slot 2
BAND = "#6b6a63"


# ---------------------------------------------------------------- data loading

def load_static(path):
    """action, t_us, event, command_u256, encoder for the two staircase actions."""
    a = np.loadtxt(path, delimiter=",", dtype=np.int64)
    return {"action": a[:, 0], "t": a[:, 1], "event": a[:, 2],
            "cmd": a[:, 3], "enc": a[:, 4]}


def load_warmup(path):
    """action, t_us, encoder, index_count for the warm-up revolutions."""
    a = np.loadtxt(path, delimiter=",", dtype=np.int64)
    return {"action": a[:, 0], "t": a[:, 1], "enc": a[:, 2], "index": a[:, 3]}


def settled(d, action):
    """Mean and spread of the settled tail of every dwell, sorted by command."""
    m = d["action"] == action
    cmd, enc = d["cmd"][m], d["enc"][m]
    edges = np.r_[0, np.flatnonzero(np.diff(cmd) != 0) + 1, len(cmd)]
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        n = max(3, int(np.ceil((b - a) * TAIL)))
        z = enc[max(a, b - n):b]
        rows.append((cmd[a], z.mean(), z.std(), b - a))
    rows = np.asarray(rows, dtype=float)
    return rows[np.argsort(rows[:, 0])].T


def paired(d, fwd_action=15, rev_action=16):
    """Forward and reverse settled positions on the common command grid."""
    cf, ef, sf, nf = settled(d, fwd_action)
    cr, er, sr, nr = settled(d, rev_action)
    cmd, i, j = np.intersect1d(cf, cr, return_indices=True)
    return {"cmd": cmd, "ef": ef[i], "er": er[j], "sf": sf[i], "sr": sr[j],
            "nf": nf[i], "nr": nr[j]}


# ------------------------------------------------------------------- harmonics

def design(angle_deg, orders, linear=True, period_deg=360.0):
    """Design matrix for a harmonic fit, and the name of each column.

    Split out so a fit made on one window can be *evaluated* on another -
    the current sweep reuses the static map's geometry fit at absolute shaft
    angle rather than refitting it over its own 28.8 degree window.
    """
    th = 2 * np.pi * np.asarray(angle_deg, dtype=float) / period_deg
    cols = [np.ones_like(th)]
    names = ["const"]
    if linear:
        cols.append(np.asarray(angle_deg, dtype=float) / period_deg)
        names.append("linear")
    for k in orders:
        cols += [np.cos(k * th), np.sin(k * th)]
        names += [f"cos{k}", f"sin{k}"]
    return np.column_stack(cols), names


def harmonic_fit(angle_deg, signal, orders, linear=True, period_deg=360.0):
    """Simultaneous least squares at exact mechanical orders.

    Projecting at the exact orders avoids the spectral leakage an FFT would
    introduce, since the common command grid is one microstep short of a
    whole revolution.
    """
    A, names = design(angle_deg, orders, linear, period_deg)
    coef, *_ = np.linalg.lstsq(A, signal, rcond=None)
    out = {}
    for k in orders:
        c = coef[names.index(f"cos{k}")]
        s = coef[names.index(f"sin{k}")]
        out[k] = {"amp": float(np.hypot(c, s)),
                  "phase_deg": float(np.degrees(np.arctan2(s, c)) % 360)}
    return coef, names, A, out


def fold(cmd, signal, period=ELEC_U256):
    """Average a signal over every electrical cycle in the revolution."""
    phase = np.mod(cmd, period).astype(int)
    keys = np.unique(phase)
    mean = np.array([signal[phase == k].mean() for k in keys])
    sem = np.array([signal[phase == k].std(ddof=1) / np.sqrt((phase == k).sum())
                    for k in keys])
    n = np.array([(phase == k).sum() for k in keys])
    return keys, mean, sem, n


def wrap(x_deg, y):
    """Close a folded cycle so the curve joins up at the period boundary."""
    return np.r_[x_deg, ELEC_DEG], np.r_[y, y[0]]


# ------------------------------------------------------------------ scale check

def counts_per_rev(w):
    """Encoder counts between consecutive index pulses in the warm-up block."""
    deltas = []
    for action in np.unique(w["action"]):
        m = w["action"] == action
        idx, enc = w["index"][m], w["enc"][m]
        hits = np.flatnonzero(np.diff(idx) != 0) + 1
        if len(hits) < 2:
            continue
        e = enc[hits].astype(float)
        deltas.append(np.abs(np.diff(e)))
    d = np.concatenate(deltas)
    # Sampling at 500 Hz during 8000 counts/s motion puts each index crossing
    # within about 16 counts, so the spread is the sampling grid, not the scale.
    return d


# ------------------------------------------------------------------- rendering

def style():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans", "font.size": 9.5,
        "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": GRID,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.labelsize": 10, "axes.titlesize": 12.5,
        "axes.linewidth": 0.8, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.7, "grid.alpha": 1.0,
        "legend.frameon": False, "legend.fontsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
    })


def titleblock(fig, title, subtitle, stats):
    fig.text(0.055, 0.975, title, size=13.5, weight="bold", color=INK,
             va="top")
    fig.text(0.055, 0.922, subtitle, size=9.6, color=INK2, va="top")
    fig.text(0.055, 0.858, stats, size=9.6, color=INK, va="top",
             weight="medium")


def outside_legend(ax, ncol):
    """Legend on its own row above the axes, so it never covers data."""
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=ncol,
              handlelength=2.2, columnspacing=1.8, borderaxespad=0.0)


def figure1(st, fit_curve, metrics):
    fig = plt.figure(figsize=(11.0, 5.8))
    ax = fig.add_axes([0.075, 0.115, 0.845, 0.645])
    ang = st["angle"]
    ax.axhline(0, color=INK3, lw=0.8, zorder=1)
    ax.plot(ang, st["f_deg"], color=FWD, lw=0.9, zorder=3,
            label="forward pass")
    ax.plot(ang, st["r_deg"], color=REV, lw=0.9, zorder=3,
            label="reverse pass")
    ax.plot(ang, fit_curve * DPC, color=INK, lw=2.1, ls=(0, (6, 3)), zorder=4,
            label="fitted 1/rev + 2/rev (encoder geometry)")
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.margins(y=0.10)
    ax.set_xlabel("shaft angle (degrees)")
    ax.set_ylabel(r"$\theta_{enc} - \theta_{cmd}$  (degrees)")
    outside_legend(ax, 3)

    sec = ax.secondary_yaxis("right", functions=(lambda v: v / DPC,
                                                 lambda v: v * DPC))
    sec.set_ylabel("encoder counts", color=INK2)
    sec.spines["right"].set_color(GRID)

    e = metrics["geometry"]
    titleblock(
        fig,
        "Figure 1  Position error over one shaft revolution",
        "Settled tail of every 150 ms dwell: all 3,199 commanded 1/16 microsteps "
        "shared by the forward (increasing command) and reverse passes.",
        f"The slow wave is encoder geometry — 1/rev amplitude "
        f"{e['one_per_rev_deg']:.3f}°, 2/rev {e['two_per_rev_deg']:.3f}°, "
        f"together {e['variance_explained_pct']:.0f}% of the variance. "
        f"The fine ripple riding on it is the motor.")
    return fig


def figure2(fd, metrics):
    fig = plt.figure(figsize=(11.0, 5.8))
    ax = fig.add_axes([0.075, 0.115, 0.845, 0.645])

    # Horizontal grid only: a vertical grid at the same tick spacing would be
    # indistinguishable from the full-step markers.
    ax.grid(False)
    ax.grid(True, axis="y")
    for k in range(4):
        ax.axvline(k * FULL_STEP_DEG, color=INK3, lw=1.2, alpha=0.7, zorder=1)
    step_proxy = mpl.lines.Line2D([], [], color=INK3, lw=1.2, alpha=0.7,
                                  label="full-step positions")

    xf, yf = wrap(fd["x"], fd["f"])
    xr, yr = wrap(fd["x"], fd["r"])
    ax.axhline(0, color=INK3, lw=0.8, zorder=2)
    band = ax.fill_between(xf, yf, yr, color=BAND, alpha=0.20, lw=0, zorder=3,
                           label="direction gap = friction hysteresis")
    lf, = ax.plot(xf, yf, color=FWD, lw=2.1, zorder=5, label="forward pass")
    lr, = ax.plot(xr, yr, color=REV, lw=2.1, zorder=5, label="reverse pass")

    ax.set_xlim(0, ELEC_DEG)
    ax.set_xticks(np.arange(0, 7.21, 0.9))
    ax.margins(y=0.12)
    ax.set_xlabel("position within one electrical cycle (degrees of shaft "
                  "angle; 7.2° = 4 full steps)")
    ax.set_ylabel("error after removing geometry fit (degrees)")
    ax.legend(handles=[lf, lr, band, step_proxy], loc="lower right",
              bbox_to_anchor=(1.0, 1.01), ncol=4, handlelength=2.2,
              columnspacing=1.8, borderaxespad=0.0)

    sec = ax.secondary_yaxis("right", functions=(lambda v: v / DPC,
                                                 lambda v: v * DPC))
    sec.set_ylabel("encoder counts", color=INK2)
    sec.spines["right"].set_color(GRID)

    m = metrics["motor"]
    titleblock(
        fig,
        "Figure 2  The motor's own error over one electrical cycle",
        "Geometry fit subtracted, then folded onto one electrical cycle and "
        "averaged over all 50 cycles.",
        f"Curve shape is the repeatable microstep error: "
        f"{m['static_pp_deg']:.3f}° peak-to-peak "
        f"({m['static_pp_pct_full_step']:.1f}% of a full step). "
        f"Band width is friction hysteresis: {m['hysteresis_half_deg']:.3f}° "
        f"mean half-width ({m['hysteresis_half_arcsec']:.0f} arcsec).")
    return fig


# -------------------------------------------------------------- working plots

def working_settling(d, path):
    """Did 150 ms actually settle? Ensemble of dwells, aligned to the step edge."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for ax, action, name, col in ((axes[0], 15, "forward", FWD),
                                  (axes[1], 16, "reverse", REV)):
        m = d["action"] == action
        cmd, enc, t = d["cmd"][m], d["enc"][m], d["t"][m]
        edges = np.r_[0, np.flatnonzero(np.diff(cmd) != 0) + 1, len(cmd)]
        grid = np.linspace(0, 150, 76)
        acc = []
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 10:
                continue
            tt = (t[a:b] - t[a]) / 1000.0
            yy = enc[a:b] - enc[b - 1]
            acc.append(np.interp(grid, tt, yy, left=np.nan, right=np.nan))
        acc = np.asarray(acc[:600])
        ax.axhline(0, color=INK3, lw=0.8)
        ax.axvspan(150 * (1 - TAIL), 150, color=BAND, alpha=0.14, lw=0,
                   label=f"settled tail ({TAIL:.0%})")
        ax.plot(grid, np.nanmedian(acc, axis=0), color=col, lw=2.0,
                label="median")
        ax.fill_between(grid, np.nanpercentile(acc, 5, axis=0),
                        np.nanpercentile(acc, 95, axis=0), color=col,
                        alpha=0.18, lw=0, label="5-95%")
        ax.set_title(f"{name} dwells", size=11)
        ax.set_xlabel("time since step edge (ms)")
        ax.set_ylabel("counts relative to dwell end")
        ax.legend(loc="upper right")
    fig.suptitle("Working plot: settling within each 150 ms dwell", size=12,
                 weight="bold", x=0.08, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


def working_spectrum(st, path):
    y = st["mean"][:-1].copy()
    x = np.linspace(-1, 1, len(y), endpoint=False)
    y = y - np.polyval(np.polyfit(x, y, 1), x)
    spec = 2 * np.abs(np.fft.rfft(y)) / len(y)
    order = np.arange(len(spec))
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.semilogy(order[1:], np.maximum(spec[1:], 1e-4), color=INK2, lw=0.9)
    for k, lab in ((1, "1/rev\ngeometry"), (2, "2/rev"), (50, "50/rev\nelectrical"),
                   (100, "100/rev"), (200, "200/rev\nfull step")):
        if k < len(spec):
            ax.plot([k], [spec[k]], "o", ms=7, color=FWD, zorder=5)
            ax.annotate(lab, xy=(k, spec[k]), xytext=(0, 11),
                        textcoords="offset points", ha="center", size=8.2,
                        color=INK2)
    ax.set_xlim(0, 420)
    ax.set_ylim(top=spec[1] * 8)          # headroom for the order-1 label
    ax.set_xlabel("mechanical order (cycles per revolution)")
    ax.set_ylabel("amplitude (counts)")
    ax.set_title("Working plot: spatial harmonics of the direction-pair mean",
                 size=12, weight="bold", loc="left", pad=16)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def working_microsteps(st, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for ax, key, name, col in ((axes[0], "ef", "forward", FWD),
                               (axes[1], "er", "reverse", REV)):
        step = np.abs(np.diff(st[key]))
        axes_ideal = 6.25
        ax.hist(step, bins=np.arange(0, 14.5, 0.5), color=col, alpha=0.75,
                edgecolor=SURFACE, lw=1.2)
        ax.axvline(axes_ideal, color=INK, lw=1.8, ls=(0, (5, 3)),
                   label=f"ideal {axes_ideal} counts")
        ax.axvline(step.mean(), color=INK2, lw=1.4,
                   label=f"mean {step.mean():.2f}")
        ax.set_xlabel("settled increment per commanded 1/16 microstep (counts)")
        ax.set_ylabel("dwells")
        ax.set_title(f"{name} pass", size=11)
        ax.legend(loc="upper right")
    fig.suptitle("Working plot: realised microstep size", size=12,
                 weight="bold", x=0.08, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------------------ main

def main():
    src_static = CACHE / "static.csv"
    src_warm = CACHE / "warmup.csv"
    if not (src_static.exists() and src_warm.exists()):
        run = MAIN / "motor_identification_run.csv"
        raise SystemExit(
            f"missing cache in {CACHE}\n\n"
            "Pre-filtering the 86 MB recording keeps this script fast. Both\n"
            "paths are overridable, which matters while the repository volume\n"
            "is full:\n\n"
            "  export STAIRCASE_CACHE=/somewhere/with/space\n"
            "  export STAIRCASE_OUT=$STAIRCASE_CACHE/staircase\n"
            "  mkdir -p \"$STAIRCASE_CACHE\"\n\n"
            "  # staircase: action, t_us, event, command_u256, encoder\n"
            "  awk -F, '$1==\"D\" && ($4==15||$4==16)"
            "{print $4\",\"$3\",\"$5\",\"$6\",\"$7}' \\\n"
            f"    '{run}' \\\n"
            "    > \"$STAIRCASE_CACHE/static.csv\"\n\n"
            "  # warm-up revolutions for the counts/rev check:\n"
            "  #   action, t_us, encoder, index_count\n"
            "  awk -F, '$1==\"D\" && $4>=1 && $4<=12"
            "{print $4\",\"$3\",\"$7\",\"$8}' \\\n"
            f"    '{run}' \\\n"
            "    > \"$STAIRCASE_CACHE/warmup.csv\"\n")

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    style()

    d = load_static(src_static)
    p = paired(d)

    ideal = p["cmd"] * C2C
    f = p["ef"] - ideal
    r = p["er"] - ideal
    mean = (f + r) / 2.0
    hys = (f - r) / 2.0
    angle = p["cmd"] / REV_U256 * 360.0

    # Geometry: 1/rev + 2/rev fitted to the direction-pair mean, so that
    # hysteresis cannot bias it. No linear term, so nothing but true geometry
    # orders is removed.
    coef, names, A, geo = harmonic_fit(angle, mean, (1, 2), linear=False)
    fit_curve = A @ coef
    resid = mean - fit_curve
    var_expl = 100.0 * (1 - np.var(resid) / np.var(mean))

    # Harmonic amplitudes for the table, projected at exact orders.
    _, _, _, orders = harmonic_fit(angle, mean, (1, 2, 50, 100, 200),
                                   linear=True)

    # Fold the geometry-corrected passes onto one electrical cycle. A single
    # common offset is removed so the band width survives the centring.
    cf = f - fit_curve
    cr = r - fit_curve
    keys, fold_f, sem_f, nf = fold(p["cmd"], cf)
    _, fold_r, sem_r, _ = fold(p["cmd"], cr)
    _, fold_m, sem_m, _ = fold(p["cmd"], (cf + cr) / 2.0)
    offset = fold_m.mean()
    fold_f -= offset
    fold_r -= offset
    fold_m -= offset
    x_deg = keys / ELEC_U256 * ELEC_DEG

    half = (fold_r - fold_f) / 2.0          # positive = rotor lags travel

    w = load_warmup(src_warm)
    cpr = counts_per_rev(w)

    static_closure = float(r[-1] - f[0])

    metrics = {
        "source": {
            "run": "data/2026-10-13_motor_identification/motor_identification_run.csv",
            "actions": [15, 16],
            "commanded_positions_paired": int(len(p["cmd"])),
            "dwell_ms": 150, "settled_tail_fraction": TAIL,
            "samples_per_dwell_median": float(np.median(p["nf"])),
        },
        "geometry": {
            "one_per_rev_counts": geo[1]["amp"],
            "one_per_rev_deg": geo[1]["amp"] * DPC,
            "one_per_rev_arcsec": geo[1]["amp"] * APC,
            "one_per_rev_phase_deg": geo[1]["phase_deg"],
            "two_per_rev_counts": geo[2]["amp"],
            "two_per_rev_deg": geo[2]["amp"] * DPC,
            "variance_explained_pct": float(var_expl),
            "implied_scale_radius_mm": 159.52 / (geo[1]["amp"] * DPC
                                                 * np.pi / 180.0) / 1000.0,
        },
        "motor": {
            "static_pp_counts": float(np.ptp(fold_m)),
            "static_pp_deg": float(np.ptp(fold_m) * DPC),
            "static_pp_arcsec": float(np.ptp(fold_m) * APC),
            "static_pp_pct_full_step": float(np.ptp(fold_m) * DPC
                                             / FULL_STEP_DEG * 100),
            "hysteresis_half_counts": float(half.mean()),
            "hysteresis_half_deg": float(half.mean() * DPC),
            "hysteresis_half_arcsec": float(half.mean() * APC),
            "hysteresis_half_pct_full_step": float(half.mean() * DPC
                                                   / FULL_STEP_DEG * 100),
            "hysteresis_half_range_counts": [float(half.min()),
                                             float(half.max())],
            "fold_sem_median_counts": float(np.median(sem_m)),
            "full_rev_corrected_pp_deg": float(np.ptp(resid) * DPC),
        },
        "harmonics_of_mean": {
            f"order_{k}": {"counts": v["amp"], "arcsec": v["amp"] * APC,
                           "deg": v["amp"] * DPC}
            for k, v in orders.items()
        },
        "checks": {
            "counts_per_rev_mean": float(cpr.mean()),
            "counts_per_rev_std": float(cpr.std()),
            "counts_per_rev_n_intervals": int(len(cpr)),
            "counts_per_rev_nominal": 20000,
            "settled_noise_median_counts": float(np.median(np.r_[p["sf"],
                                                                 p["sr"]])),
            "staircase_closure_counts": static_closure,
            "staircase_closure_arcsec": static_closure * APC,
            "whole_run_drift_counts": 267.0,
        },
    }

    st = {"angle": angle, "mean": mean, "hys": hys,
          "f_deg": f * DPC, "r_deg": r * DPC,
          "ef": p["ef"], "er": p["er"]}
    fd = {"x": x_deg, "f": fold_f * DPC, "r": fold_r * DPC,
          "m": fold_m * DPC}

    f1 = figure1(st, fit_curve, metrics)
    f1.savefig(OUT / "fig1_error_over_revolution.png")
    plt.close(f1)

    f2 = figure2(fd, metrics)
    f2.savefig(OUT / "fig2_electrical_cycle.png")
    plt.close(f2)

    working_settling(d, WORK / "w1_settling.png")
    working_spectrum(st, WORK / "w2_spectrum.png")
    working_microsteps(st, WORK / "w3_microstep_size.png")

    (OUT / "staircase_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    np.savetxt(OUT / "fig2_electrical_cycle.csv",
               np.column_stack([x_deg, fold_f * DPC, fold_r * DPC,
                                fold_m * DPC, half * DPC]),
               delimiter=",", comments="",
               header="shaft_deg_in_cycle,forward_deg,reverse_deg,mean_deg,"
                      "hysteresis_half_deg", fmt="%.6f")

    g, mo, ch = metrics["geometry"], metrics["motor"], metrics["checks"]
    print(f"paired positions        {len(p['cmd'])}")
    print(f"1/rev amplitude         {g['one_per_rev_counts']:.2f} counts "
          f"= {g['one_per_rev_deg']:.4f} deg")
    print(f"2/rev amplitude         {g['two_per_rev_counts']:.2f} counts "
          f"= {g['two_per_rev_deg']:.4f} deg")
    print(f"fit explains            {g['variance_explained_pct']:.1f}% of variance")
    print(f"implied scale radius    {g['implied_scale_radius_mm']:.1f} mm")
    print(f"static error pp         {mo['static_pp_deg']:.4f} deg "
          f"= {mo['static_pp_pct_full_step']:.2f}% of a full step")
    print(f"hysteresis half-width   {mo['hysteresis_half_deg']:.4f} deg "
          f"= {mo['hysteresis_half_arcsec']:.1f} arcsec")
    print(f"fold SEM (median)       {mo['fold_sem_median_counts']:.3f} counts")
    for k in (1, 2, 50, 100, 200):
        v = metrics["harmonics_of_mean"][f"order_{k}"]
        print(f"order {k:<4}             {v['counts']:7.3f} counts "
              f"= {v['arcsec']:8.2f} arcsec")
    print(f"counts/rev              {ch['counts_per_rev_mean']:.1f} "
          f"+- {ch['counts_per_rev_std']:.1f} "
          f"({ch['counts_per_rev_n_intervals']} index intervals)")
    print(f"settled noise median    {ch['settled_noise_median_counts']:.3f} counts")
    print(f"staircase closure       {ch['staircase_closure_counts']:+.1f} counts")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
