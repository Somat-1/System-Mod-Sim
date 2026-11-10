#!/usr/bin/env python3
"""Companion figures for F1 (V1, P1, P2, M1, Q1); illustrative model curves only.

All curves come from the quasi-static Maxwell-slip model used for the F0b S
animation (placeholder slip limits, saturation beta0/2 per side), driven along
the F1 command paths. Lag e = x - u: an up branch lies below, a down branch
above. X1 (position comparison and hold-out) is data-only and is not drawn.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

from f1_sequence import CONFIG, strokes
from render_f1_sequence import GRID, INK, MUTED, style

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "F0b Friction Scoping"))
from render_f0b_s_animation import GAMMA0, LIMITS, advance  # noqa: E402

OUT = HERE / "rendered" / "companion"
C = json.loads(CONFIG.read_text(encoding="utf-8"))
COL = C["colors"]
BETA = C["beta_0_um"]
S90_CFG = C["f0b_dependent"]["s90_um"]
DMIN = C["f0b_dependent"]["d_min_um"]
G = GAMMA0 * (BETA / 2) / (GAMMA0 @ LIMITS)
GUIDE = "#9aa2a8"


# --------------------------------------------------------------------- model

def run(z0, u_points):
    """Drive the model through commanded positions u_points; return u, e, final z."""
    z, e = np.array(z0, float), [-G @ z0]
    for a, b in zip(u_points[:-1], u_points[1:]):
        z, _ = advance(z, b - a, G)
        e.append(-G @ z)
    return np.asarray(u_points, float), np.asarray(e), z


def path(waypoints, n=200):
    """Piecewise-linear command path, n samples per leg (midpoints included)."""
    pts = [waypoints[0]]
    for a, b in zip(waypoints[:-1], waypoints[1:]):
        pts.extend(np.linspace(a, b, n + 1)[1:])
    return np.asarray(pts)


def virgin(direction, travel):
    """Reversal curve after a saturating wipe: direction -1 = down after wipe from below."""
    z0 = -direction * LIMITS            # wipe ends moving opposite to the ramp
    _, e, _ = run(z0, np.r_[0.0, direction * travel])
    return np.r_[0.0, travel], e


def loop(stroke, cycles=3, n=200):
    """Centred loop after a wipe from below to x_p; returns per-cycle up/down branches."""
    h = stroke / 2
    wp = [0.0, -h] + [h, -h] * cycles + [0.0]
    u, e, _ = run(LIMITS, path(wp, n))
    branches = []
    for k in range(cycles):
        i0 = n + 2 * k * n                  # start of cycle k (at -h)
        up = slice(i0, i0 + n + 1)
        dn = slice(i0 + n, i0 + 2 * n + 1)
        branches.append(((u[up], e[up]), (u[dn], e[dn])))
    return branches


def width_at_centre(stroke):
    (uu, eu), (ud, ed) = loop(stroke, cycles=3)[-1]
    return np.interp(0, ud[::-1], ed[::-1]) - np.interp(0, uu, eu)


def s90_model():
    s = np.geomspace(1e-3, 100, 4000)
    _, e = virgin(-1, s)
    flip = np.abs(e[1:] - e[0])
    return s[np.argmax(flip >= .9 * BETA)]


# ------------------------------------------------------------------ helpers

def finish(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(color=GRID, lw=.5)


def illustrative(ax, x=.99, y=.02, ha="right"):
    ax.text(x, y, "illustrative model", transform=ax.transAxes, ha=ha, va="bottom",
            fontsize=7, color=MUTED, style="italic")


def arrow_on(ax, u, e, frac, color, size=8):
    i = int(frac * (len(u) - 1))
    ax.annotate("", xy=(u[i + 2], e[i + 2]), xytext=(u[i - 2], e[i - 2]),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0, mutation_scale=size))


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    path_ = OUT / name
    fig.savefig(path_, dpi=300, bbox_inches="tight", pad_inches=.03)
    plt.close(fig)
    return path_


# ----------------------------------------------------------------------- V1

def figure_v1(s90):
    col = COL["V"]
    travel = np.geomspace(.01, C["v_ramp_um"], 1500)
    sd, ed = virgin(-1, travel)
    su, eu = virgin(+1, travel)
    fig, ax = plt.subplots(figsize=(5.2, 3.3))
    ax.plot(sd[1:], ed[1:], color=col, lw=1.6, label="reversal down (after wipe from below)")
    ax.plot(su[1:], eu[1:], color=col, lw=1.6, ls=(0, (4, 2)),
            label="reversal up (after wipe from above)")
    for v in (BETA / 2, -BETA / 2):
        ax.axhline(v, color=GUIDE, lw=.7, ls=":")
    ax.axvline(s90, color=INK, lw=.7, ls=(0, (3, 2)))
    ax.text(s90 / 1.1, -1.2, f"$s_{{90}}$ = {s90:.1f} µm\npre-sliding length", fontsize=7,
            color=INK, va="bottom", ha="right", linespacing=.95)
    # Regimes along the down curve.
    ax.annotate("carriage nearly stuck\n($|de/ds| \\approx 1$)", xy=(.06, eu[np.searchsorted(su, .06)]),
                xytext=(.012, .55), fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6))
    ax.annotate("pre-sliding: contacts\ndeform and slip one by one",
                xy=(1.2, ed[np.searchsorted(sd, 1.2)]), xytext=(.2, 1.02), fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6))
    ax.annotate("", xy=(55, BETA / 2), xytext=(55, -BETA / 2),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=.7, mutation_scale=7))
    ax.text(48, 0, f"sliding: lag flipped\nby $\\beta_0$ = {BETA:.2f} µm", fontsize=7, color=INK,
            ha="right", va="center", linespacing=.95)
    # Slow / fast ramp strip along the bottom.
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    slow = C["v_slow_factor_of_s90"] * S90_CFG
    ax.fill_between([.01, slow], 0, .035, transform=tr, color=col, alpha=.55, lw=0)
    ax.fill_between([slow, C["v_ramp_um"]], 0, .035, transform=tr, color=col, alpha=.18, lw=0)
    ax.text(np.sqrt(.01 * slow), .045, f"{C['v_slow_um_s']:g} µm/s", transform=tr, fontsize=6.5,
            color=MUTED, ha="center")
    ax.text(np.sqrt(slow * C["v_ramp_um"]), .045, f"{C['v_fast_um_s']:g} µm/s", transform=tr,
            fontsize=6.5, color=MUTED, ha="center")
    ax.set_xscale("log")
    ax.set(xlim=(.01, C["v_ramp_um"]), ylim=(-1.35, 1.35), xlabel="travel since reversal [µm]",
           ylabel=r"lag $e = x - u$ [µm]")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -.2), ncol=2, frameon=False, fontsize=7)
    illustrative(ax, x=.01, y=.93, ha="left")
    finish(ax)
    # Inset: the commanded V path (wipe from below, then reversal down).
    ins = ax.inset_axes([.03, .33, .2, .24])
    ins.plot([0, 1, 2], [0, -1, 0], color=COL["setup"], lw=1.2)
    ins.plot([2, 2.4], [0, 0], color=COL["setup"], lw=1.2)
    ins.annotate("", xy=(4.4, -.55), xytext=(2.4, 0),
                 arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3, mutation_scale=8))
    ins.text(.5, -1.0, "wipe", fontsize=6.5, color=MUTED, ha="center", va="top")
    ins.text(3.5, -.05, "reversal", fontsize=6.5, color=col, ha="center", va="bottom")
    ins.set(xlim=(-.2, 4.7), ylim=(-1.35, .45), xticks=[], yticks=[])
    ins.set_ylabel("$u$", fontsize=7, labelpad=1)
    ins.spines[["top", "right"]].set_visible(False)
    return save(fig, "f1_V1_virgin_curve.png")


# ----------------------------------------------------------------------- P1

def figure_p1():
    col = COL["P"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(width_ratios=[1, 1.15],
                                                                        wspace=.3))
    picks = [(.3 * S90_CFG, .45, (0, (1.2, 1.2)), "0.3 $s_{90}$"),
             (S90_CFG, .7, (0, (4, 2)), "$s_{90}$"),
             (100.0, 1.0, "-", "100 µm")]
    for s, alpha, ls, lab in picks:
        (uu, eu), (ud, ed) = loop(s)[-1]
        h = s / 2
        a.plot(uu / h, eu, color=col, alpha=alpha, lw=1.4, ls=ls, label=lab)
        a.plot(ud / h, ed, color=col, alpha=alpha, lw=1.4, ls=ls)
        if s == 100.0:
            arrow_on(a, uu / h, eu, .5, col)
            arrow_on(a, ud / h, ed, .5, col)
    for v in (BETA / 2, -BETA / 2):
        a.axhline(v, color=GUIDE, lw=.7, ls=":")
    a.set(xlim=(-1.08, 1.08), ylim=(-1.15, 1.15), xticks=[-1, 0, 1],
          xlabel=r"$(u - x_p)\,/\,(s/2)$", ylabel=r"lag $e$ [µm]")
    a.legend(fontsize=7, frameon=False, loc="lower right", bbox_to_anchor=(1, .99), ncol=3,
             handlelength=1.8, columnspacing=1.0)
    a.set_title("(a)", loc="left", pad=3)
    illustrative(a)
    finish(a)

    s_line = np.geomspace(.05, 5000, 70)
    w_line = np.array([width_at_centre(s) for s in s_line])
    grid = np.array(strokes(C))
    w_grid = np.array([width_at_centre(s) for s in grid])
    b.axvspan(S90_CFG / 2, 2 * S90_CFG, color="#9aa7ad", alpha=.16, lw=0)
    b.text(S90_CFG, 1.97, "$s_{90}$ band\n(F0b)", ha="center", va="top", fontsize=7,
           color="#5d6b72", linespacing=.95)
    b.axhline(BETA, color=GUIDE, lw=.7, ls=":")
    b.text(.06, BETA + .05, f"$\\beta_0$ = {BETA:.2f} µm", fontsize=7, color=MUTED)
    b.plot(s_line, w_line, color=col, lw=1.1, alpha=.6)
    b.scatter(grid, w_grid, s=16, facecolor="white", edgecolor=col, lw=.9, zorder=3,
              label="F1 strokes")
    b.set_xscale("log")
    b.set(xlim=(.05, 5000), ylim=(0, 2.05), xlabel="stroke $s$ [µm]",
          ylabel=r"loop width $w(x_p) = e_\downarrow - e_\uparrow$ [µm]")
    b.text(.97, .3, "same commanded position on both\nbranches: geometric and microstep\n"
           "error cancel", transform=b.transAxes, ha="right", va="bottom", fontsize=7,
           color=INK, linespacing=1.0)
    b.set_title("(b)", loc="left", pad=3)
    illustrative(b)
    finish(b)
    return save(fig, "f1_P1_loops_and_width.png")


# ----------------------------------------------------------------------- P2

def figure_p2():
    col = COL["P"]
    s = 2 * S90_CFG
    h = s / 2
    (uu, eu), (ud, ed) = loop(s, n=400)[-1]
    # Masing: each branch is twice the initial-loading curve f, stretched by two,
    # evaluated on the branch's own samples (f has kinks; no interpolation).
    d = uu - uu[0]
    _, e0, _ = run(np.zeros_like(LIMITS), d / 2)            # initial loading from rest, up
    f = -e0                                                 # magnitude of the loading curve
    up_m = eu[0] - 2 * f
    dn_m = ed[0] + 2 * f
    assert np.max(np.abs(eu - up_m)) < 1e-9
    assert np.max(np.abs(ed - dn_m)) < 1e-9

    fig = plt.figure(figsize=(4.4, 3.8))
    a = fig.add_axes([.15, .36, .8, .58])
    r = fig.add_axes([.15, .12, .8, .14], sharex=a)
    a.plot(uu, eu, color=col, lw=1.8, label="loop branches")
    a.plot(ud, ed, color=col, lw=1.8)
    a.plot(-h + d, up_m, color=INK, lw=1.0, ls=(0, (3, 2)),
           label=r"$2\,f(\Delta u/2)$ from each turning point")
    a.plot(h - d, dn_m, color=INK, lw=1.0, ls=(0, (3, 2)))
    a.plot([-h, h], [eu[0], ed[0]], "o", ms=4, color=INK, zorder=4)
    arrow_on(a, uu, eu, .55, col)
    arrow_on(a, ud, ed, .55, col)
    a.set(ylabel=r"lag $e$ [µm]", ylim=(-1.15, 1.15))
    plt.setp(a.get_xticklabels(), visible=False)
    a.legend(loc="lower left", fontsize=7, frameon=False)
    a.text(.98, .97, f"stroke {s:g} µm $\\approx$ 2 $s_{{90}}$\n$f$: loading curve from rest",
           transform=a.transAxes, ha="right", va="top", fontsize=7, color=MUTED,
           linespacing=1.0)
    a.set_title("(a)", loc="left", pad=3)
    finish(a)
    r.axhspan(-DMIN * 1000, DMIN * 1000, color="#9aa7ad", alpha=.25, lw=0)
    r.axhline(0, color=INK, lw=.9)
    r.text(.01, .9, r"with data: branch $-\ 2f$ ; band $\pm d_{\min}$", transform=r.transAxes,
           fontsize=6.5, color=MUTED, va="top")
    r.set(xlim=(-h * 1.08, h * 1.08), ylim=(-12, 12), yticks=[-10, 0, 10],
          xlabel=r"$u - x_p$ [µm]", ylabel="[nm]")
    r.set_title("(b)", loc="left", pad=2)
    finish(r)
    illustrative(a, y=.78)
    return save(fig, "f1_P2_masing_overlay.png")


# ----------------------------------------------------------------------- M1

def figure_m1():
    col = COL["M"]
    so = C["m_outer_factor"] * S90_CFG
    si = C["m_inner_factor"] * S90_CFG
    n = 300
    wp = [-so / 2, so / 2, 0.0, si, 0.0, -so / 2]
    u, e, _ = run(LIMITS, path(wp, n))
    ref_u, ref_e, _ = run(LIMITS, path([-so / 2, so / 2, -so / 2], n))
    i_off, i_on = 2 * n, 4 * n                      # arrival at x_p, end of inner loop
    rp_err = e[i_on] - e[i_off]
    assert abs(rp_err) < 1e-9

    fig, ax = plt.subplots(figsize=(4.9, 3.4))
    ax.plot(ref_u[n:], ref_e[n:], color=GUIDE, lw=3.2, alpha=.45,
            label="outer down-branch, uninterrupted")
    ax.plot(u[:i_off + 1], e[:i_off + 1], color=col, lw=1.5, label="outer loop")
    ax.plot(u[i_on:], e[i_on:], color=col, lw=1.5)
    ax.plot(u[i_off:i_on + 1], e[i_off:i_on + 1], color=INK, lw=1.1, label="inner loop")
    for frac, (lo, hi) in ((.5, (0, n)), (.5, (n, 2 * n))):
        arrow_on(ax, u[lo:hi + 1], e[lo:hi + 1], frac, col)
    ax.plot(u[i_off], e[i_off], "o", ms=5, mfc="white", mec=INK, mew=1, zorder=5)
    ax.set(xlim=(-so / 2 - .6, so / 2 + .6), ylim=(-1.15, 1.15),
           xlabel=r"$u - x_p$ [µm]", ylabel=r"lag $e$ [µm]")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=7, frameon=False,
              handlelength=1.8, columnspacing=1.0)
    ax.text(-so / 2 + .1, .66, "GMS: after the inner loop the\nstage resumes the outer branch\n"
            "where it left it (return-point\nmemory); a single-state model\nwould not",
            fontsize=7, color=INK, va="top", linespacing=1.0)
    illustrative(ax, x=.01, ha="left")
    finish(ax)
    # Zoom on the branch-off / rejoin point.
    ins = ax.inset_axes([.07, .1, .34, .36])
    sl = slice(i_off - 60, i_on + 61)
    ins.plot(ref_u[n:], ref_e[n:], color=GUIDE, lw=3.2, alpha=.45)
    ins.plot(u[sl], e[sl], color=INK, lw=1)
    ins.plot(u[i_off], e[i_off], "o", ms=5, mfc="white", mec=INK, mew=1)
    ins.set(xlim=(-.5, si + .25), ylim=(e[i_off] - .55, e[i_off] + .25), xticks=[0, si],
            yticks=[])
    ins.set_xticklabels(["$x_p$", f"+{si:g}"], fontsize=6.5)
    ins.set_title("branch-off = rejoin;\nreturn-point error $\\approx$ 0", fontsize=6.5, pad=2,
                  color=INK, linespacing=.95)
    ins.tick_params(length=2, pad=1)
    ax.indicate_inset_zoom(ins, edgecolor=MUTED, lw=.6)
    return save(fig, "f1_M1_nested_loop.png")


# ----------------------------------------------------------------------- Q1

def figure_q1():
    col = COL["Q"]
    s = strokes({**C, "p_stroke_factors": [1.0], "p_anchor_strokes_um": []})[0]
    (uu, eu), (ud, ed) = loop(s)[-1]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw=dict(width_ratios=[1.25, 1],
                                                                        wspace=.32))
    styles = [("1.3 s", "-", 1.6), ("4 s", (0, (4, 2)), 1.4), ("12 s", (0, (1, 1.3)), 1.4)]
    for lab, ls, lw in styles:
        a.plot(uu, eu, color=col, lw=lw, ls=ls, label=f"half-cycle {lab}")
        a.plot(ud, ed, color=col, lw=lw, ls=ls)
    a.plot(uu, eu, color=INK, lw=.5, label="4 s repeat (drift)")
    a.plot(ud, ed, color=INK, lw=.5)
    arrow_on(a, uu, eu, .5, col)
    arrow_on(a, ud, ed, .5, col)
    a.set(xlabel=r"$u - x_p$ [µm]", ylabel=r"lag $e$ [µm]", ylim=(-1.1, 1.1))
    a.legend(fontsize=7, frameon=False, loc="lower left")
    a.text(.98, .96, f"stroke $s_{{90}}$ = {s:.2f} µm", transform=a.transAxes, ha="right",
           va="top", fontsize=7, color=MUTED)
    a.set_title("(a)", loc="left", pad=3)
    finish(a)
    w = width_at_centre(s)
    order = [(4, "4"), (1.3, "1.3"), (12, "12"), (4, "4 (repeat)")]
    b.axhspan(w - DMIN, w + DMIN, color="#9aa7ad", alpha=.25, lw=0)
    b.text(1.05, w + DMIN * 1.15, r"$\pm d_{\min}$", fontsize=7, color=MUTED, va="bottom")
    for T, lab in order:
        rep = "repeat" in lab
        b.plot(T * (1.07 if rep else 1), w, "o", ms=6, mfc="white" if rep else col, mec=col,
               mew=1.1)
    b.set_xscale("log")
    b.set(xlim=(1, 16), ylim=(w - 4 * DMIN, w + 4 * DMIN), xticks=[1.3, 4, 12],
          xlabel="half-cycle time [s]", ylabel=r"loop width $w(x_p)$ [µm]")
    b.set_xticklabels(["1.3", "4", "12"])
    b.text(.5, .1, "rate-independent model:\nidentical widths", transform=b.transAxes,
           ha="center", fontsize=7, color=INK, linespacing=1.0)
    b.set_title("(b)", loc="left", pad=3)
    finish(b)
    illustrative(b, y=.88)
    return save(fig, "f1_Q1_rate_check.png")


def main():
    style()
    s90 = s90_model()
    for p in (figure_v1(s90), figure_p1(), figure_p2(), figure_m1(), figure_q1()):
        print(p)
    print(f"model s90 = {s90:.2f} um (config placeholder {S90_CFG:g} um)")


if __name__ == "__main__":
    main()
