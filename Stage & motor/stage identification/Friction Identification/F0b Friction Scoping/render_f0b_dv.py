#!/usr/bin/env python3
"""Standalone figures for the F0b D (dwell) and V (speed) blocks.

Each figure pairs the commanded segments, drawn on an event axis in the S
render's visual grammar (grey positioning, block colour measured, arrows for
direction), with an illustrative sketch of the expected observable. Expected
signals are placeholders, not F0b measurements.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch
import numpy as np

from render_f0b_sequence_details import GRID, INK, MUTED, SETUP, colours, finish, style
from render_f0b_s_animation import GAMMA0, LIMITS, ramp

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "sequence_details"
DARK = "#2b3640"
SETTLE_S, D_RAMP_UM, D_RATE = 2.0, 5.0, .25     # s, um, um/s
DWELLS = (1, 60, 60, 1)                         # ABBA
STATIC_RISE_UM = .1                             # placeholder dwell effect
CREEP_NM = 30.0                                 # order of the E5 hold creep
DELTA_UM = {"fast": .88, "slow": .95}           # placeholder speed effect
ARROW = dict(arrowstyle="-|>", mutation_scale=9, shrinkA=0, shrinkB=0)


def bracket(ax, x0, x1, y, text, color=MUTED, h=.06, fs=7.5, transform=None):
    """Horizontal bracket with end ticks and a centred label above it."""
    tr = transform or ax.get_xaxis_transform()
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], color=color, lw=.8,
            transform=tr, clip_on=False)
    ax.text((x0 + x1) / 2, y + .02, text, ha="center", va="bottom", fontsize=fs,
            color=color, transform=tr)


# --------------------------------------------------------------------------- D

def d_lag_curves(beta):
    """1 s: S-model reversal curve. 60 s: stuck line held for an extra static rise."""
    gamma = GAMMA0 * (beta / 2) / (GAMMA0 @ LIMITS)
    r = ramp(+1, gamma)                        # wipe from below, then ramp down
    keep = r["s"] <= D_RAMP_UM + 1e-9
    s, e1 = r["s"][keep], r["e"][keep]
    k0 = r["slope"][0]                         # initial (stuck) slope
    stuck = k0 * s
    s_star = s[np.argmax(stuck < e1 - STATIC_RISE_UM)]
    e60 = np.where(s < s_star, stuck,
                   e1 - STATIC_RISE_UM * np.exp(-(s - s_star) / .06))
    return s, e1, e60, s_star


def figure_d(config, beta):
    c = colours(config)["D"]
    fig = plt.figure(figsize=(7.0, 6.4))
    ev = fig.add_axes([.1, .76, .88, .18])
    ruler = fig.add_axes([.1, .56, .88, .085])
    lag = fig.add_axes([.1, .08, .88, .34])

    # (a) four trials on an event axis; wipe and amplitudes are symbolic.
    widths = {1: .14, 60: .5}
    x, centres = 0.0, []
    for k, dwell in enumerate(DWELLS):
        w0 = x
        ev.annotate("", xy=(x + .16, 0), xytext=(x + .04, -1.7),
                    arrowprops=dict(**ARROW, color=SETUP, lw=1.6))
        x += .2
        ev.plot([x, x + widths[dwell]], [0, 0], color=c, lw=7, solid_capstyle="butt")
        label_y, va, tcol = (0, "center", "white") if dwell == 60 else (.32, "bottom", c)
        ev.text(x + widths[dwell] / 2, label_y, rf"$T_d$ = {dwell} s", ha="center",
                va=va, fontsize=7.5, color=tcol)
        x += widths[dwell]
        ev.annotate("", xy=(x + .3, -1), xytext=(x, 0),
                    arrowprops=dict(**ARROW, color=INK, lw=1.6))
        x += .3
        centres.append((w0 + x) / 2)
        if k == 0:
            ev.text(w0 + .02, -1.85, "3 mm wipe\nfrom below", ha="left", va="top",
                    fontsize=7, color=MUTED, linespacing=.95)
            ev.text(x + .02, -.62, "5 µm\nreversal", ha="left", va="center", fontsize=7,
                    color=INK, linespacing=.95)
        x += .22
    ev.set_xlim(-.05, x - .1)
    ev.set_ylim(-2.4, .75)
    ev.set_xticks(centres, [f"D{k}" for k in range(1, 5)])
    ev.set_yticks([0, -1], [r"$x_0$", r"$x_0 - 5$ µm"])
    ev.tick_params(length=0)
    ev.spines[["top", "right", "left"]].set_visible(False)
    bracket(ev, centres[0] - .3, centres[-1] + .3, 1.12, "ABBA order cancels drift")
    ev.set_title("(a)", loc="left", pad=16, fontsize=9, color=INK)

    # (b) time ruler: one trial of each kind, aligned at the reversal.
    T = D_RAMP_UM / D_RATE
    for yy, dwell in ((1, 60), (0, 1)):
        ruler.barh(yy, SETTLE_S, left=-dwell - SETTLE_S, height=.56, color=SETUP, lw=0)
        ruler.barh(yy, dwell, left=-dwell, height=.56, color=c, lw=0)
        ruler.barh(yy, T, left=0, height=.56, color=INK, lw=0)
    ruler.axvline(0, color="#7d858c", lw=.7, ls=(0, (3, 2)), zorder=0)
    ruler.set_yticks([1, 0], [r"$T_d$ = 60 s", r"$T_d$ = 1 s"])
    ruler.set(ylim=(-.6, 1.6), xlim=(-66, 22))
    ruler.xaxis.tick_top()
    ruler.xaxis.set_label_position("top")
    ruler.set_xlabel("time from reversal [s]", labelpad=3)
    ruler.spines[["bottom", "right"]].set_visible(False)
    ruler.set_axisbelow(True)
    ruler.grid(axis="x", color=GRID, lw=.55)
    ruler.tick_params(axis="y", length=0)
    ruler.legend(handles=[Patch(color=SETUP, label="settle"), Patch(color=c, label="dwell"),
                          Patch(color=INK, label="reversal ramp")],
                 loc="upper left", bbox_to_anchor=(0, -.08), ncol=3, frameon=False,
                 fontsize=7.5, handlelength=1.2, columnspacing=1.2, borderaxespad=0)
    fig.text(.02, .67, "(b)", fontsize=9, color=INK)

    # (c) expected observable over the ramp; connectors make (b) its time ruler.
    s, e1, e60, s_star = d_lag_curves(beta)
    lag.plot(s, e1, color=c, lw=1.6, label=r"$T_d$ = 1 s")
    lag.plot(s, e60, color=c, lw=1.6, ls=(0, (4, 2)), label=r"$T_d$ = 60 s")
    lag.plot([0, 1.9], [0, -1.9], color=GREY_GUIDE, ls=":", lw=.9)
    lag.text(1.95, -1.9, "slope −1", fontsize=7, color=MUTED, va="center")
    i = np.argmin(np.abs(s - s_star))
    lag.annotate("", xy=(s_star + .06, e60[i]), xytext=(s_star + .06, e1[i]),
                 arrowprops=dict(arrowstyle="<->", color=INK, lw=.8, shrinkA=0, shrinkB=0,
                                 mutation_scale=7))
    lag.annotate("dwell effect,\nmeasured here", xy=(s_star + .1, (e1[i] + e60[i]) / 2),
                 xytext=(1.25, -.42), fontsize=7.5, color=INK, va="center",
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6))
    lag.set(xlim=(0, D_RAMP_UM), ylim=(-2.0, .15), xlabel="travel after reversal [µm]",
            ylabel="lag $e$ [µm]")
    lag.legend(loc="lower left", fontsize=7.5, frameon=False)
    lag.text(.995, .02, "illustrative", transform=lag.transAxes, ha="right", va="bottom",
             fontsize=7, color=MUTED, style="italic")
    top = lag.secondary_xaxis("top", functions=(lambda u: u / D_RATE, lambda t: t * D_RATE))
    top.set_xticks([0, 5, 10, 15, 20], ["0", "5", "10", "15", "20 s"])
    top.tick_params(labelsize=7, colors=MUTED, pad=1)
    fig.text(.02, .44, "(c)", fontsize=9, color=INK)
    finish(lag)
    for xr, xl in ((0, 0), (T, 1)):
        fig.add_artist(ConnectionPatch(xyA=(xr, -.6), coordsA=ruler.transData,
                                       xyB=(xl, 1.09), coordsB=lag.transAxes,
                                       color="#9aa2a8", lw=.6, ls=(0, (2, 2))))

    # Creep during the 60 s dwell: the other D observable.
    ins = lag.inset_axes([.58, .5, .38, .36])
    t = np.linspace(0, 60, 400)
    rng = np.random.default_rng(3)
    creep = CREEP_NM * (1 - np.exp(-t / 18)) + rng.normal(0, 1.5, t.size)
    ins.plot(t, creep, color=c, lw=.8)
    ins.set(xlim=(0, 60), ylim=(-6, 42), xticks=[0, 30, 60], yticks=[0, 20, 40])
    ins.tick_params(labelsize=6.5, pad=1)
    ins.set_xlabel("time in 60 s dwell [s]", fontsize=7, labelpad=1)
    ins.set_ylabel(r"$x - x(0)$ [nm]", fontsize=7, labelpad=1)
    ins.set_title("creep during dwell", fontsize=7.5, loc="left", pad=2, color=INK)
    finish(ins)

    path = OUT / "f0b_D_dwell_probe.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


GREY_GUIDE = "#9aa2a8"


# --------------------------------------------------------------------------- V

def figure_v(config):
    c = colours(config)["V"]
    fig = plt.figure(figsize=(7.0, 5.4))
    ev = fig.add_axes([.09, .67, .89, .22])
    fast = fig.add_axes([.09, .1, .4, .4])
    slow = fig.add_axes([.58, .1, .4, .4], sharey=fast)

    # (a) eight equal slots: speed by line weight, direction by arrow.
    slots = [("fast", 1)] * 0 + [("fast", d) for d in (1, -1, 1, -1)] + \
            [("slow", d) for d in (1, -1)] + [("fast", d) for d in (1, -1)]
    for i, (speed, d) in enumerate(slots):
        y0, y1 = (0, 1) if d > 0 else (1, 0)
        ev.annotate("", xy=(i + .85, y1), xytext=(i + .15, y0),
                    arrowprops=dict(arrowstyle="-|>", mutation_scale=10 if speed == "fast" else 8,
                                    shrinkA=0, shrinkB=0, color=c,
                                    lw=2.4 if speed == "fast" else .9))
        ev.text(i + .5, -.36, "1 s" if speed == "fast" else "100 s", ha="center",
                va="top", fontsize=7, color=MUTED)
    for x in range(1, 8):
        ev.axvline(x, color=GRID, lw=.5, zorder=0)
    ev.set_xticks(np.arange(8) + .5, [rf"{sp} $\uparrow$" if d > 0 else rf"{sp} $\downarrow$"
                                      for sp, d in slots])
    ev.tick_params(axis="x", length=0, pad=3)
    ev.set(xlim=(0, 8), ylim=(-.12, 1.12), yticks=[0, 1],
           yticklabels=[r"$x_0$", r"$x_0$ + 1 mm"])
    ev.spines[["top", "right"]].set_visible(False)
    tr = ev.get_xaxis_transform()
    for x0, x1, txt in ((.1, 3.9, "fast"), (4.1, 5.9, "slow"), (6.1, 7.9, "fast")):
        bracket(ev, x0, x1, 1.06, txt, color=c, h=.04)
    bracket(ev, .1, 7.9, 1.27, "fast – slow – fast brackets drift", h=.04)
    ev.set_title("(a)", loc="left", pad=30, fontsize=9, color=INK)

    # (b, c) expected observable: lag over the segment, up vs down, per speed.
    x = np.linspace(0, 1, 400)
    G = 1.25 * np.sin(np.pi * x + .7) + .35 * np.sin(2 * np.pi * x + 1.9)   # 2 mm screw wobble
    for ax, key, title in ((fast, "fast", "(b)  1 mm/s"), (slow, "slow", "(c)  0.01 mm/s")):
        d = DELTA_UM[key]
        ax.plot(x, G - d, color=c, lw=1.4, label=r"up, $e_\uparrow$")
        ax.plot(x, G + d, color=c, lw=1.4, ls=(0, (4, 2)), label=r"down, $e_\downarrow$")
        ax.plot(x, G, color=GREY_GUIDE, lw=.9, ls=":", label="geometry (cancels)")
        ax.axhline(d, color=DARK, lw=2.2, label=r"$\frac{1}{2}(e_\downarrow - e_\uparrow)$")
        for xa, y, sgn in ((.75, G[299] - d, 1), (.3, G[120] + d, -1)):
            ax.annotate("", xy=(xa + .07 * sgn, y + .0), xytext=(xa - .07 * sgn, y),
                        arrowprops=dict(arrowstyle="->", color=c, lw=.9), zorder=4)
        j = 220
        ax.annotate("", xy=(x[j], G[j] + d), xytext=(x[j], G[j] - d),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=.7, shrinkA=0,
                                    shrinkB=0, mutation_scale=7))
        ax.text(x[j] + .02, G[j] - .1, r"$2\delta$", fontsize=8, color=INK, va="center")
        ax.text(.03, d - .12, rf"$\delta$ = {d:.2f} µm", ha="left", va="top",
                fontsize=7.5, color=DARK)
        ax.set(xlim=(0, 1), ylim=(-3.1, 3.1), xticks=[0, .5, 1],
               xlabel=r"$u - x_0$ [mm]")
        ax.set_title(title, loc="left", pad=4, fontsize=9, color=INK)
        ax.text(.02, .03, "illustrative", transform=ax.transAxes, fontsize=7,
                color=MUTED, style="italic")
        finish(ax)
    fast.set_ylabel("lag $e$ [µm]")
    plt.setp(slow.get_yticklabels(), visible=False)
    fast.legend(loc="lower center", bbox_to_anchor=(1.12, 1.1), ncol=4, frameon=False,
                fontsize=7.5, handlelength=2.0, columnspacing=1.3)
    # Explicit comparison between the two flat lines.
    y_fig = fast.transData.transform((0, DELTA_UM["fast"]))
    y_fig = fig.transFigure.inverted().transform(y_fig)[1]
    fig.text(.535, y_fig + .045, "compare\n" + r"$\delta$", ha="center", va="bottom",
             fontsize=7.5, color=INK, linespacing=1.1)
    for ax, key, xe, xs in ((fast, "fast", 1.0, .525), (slow, "slow", 0.0, .545)):
        fig.add_artist(ConnectionPatch(xyA=(xs, y_fig + .035), coordsA=fig.transFigure,
                                       xyB=(xe, DELTA_UM[key]), coordsB=ax.transData,
                                       arrowstyle="->", color=INK, lw=.6, mutation_scale=8,
                                       shrinkB=1))
    path = OUT / "f0b_V_speed_probe.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def main():
    config = json.loads((HERE / "f0b_config.json").read_text(encoding="utf-8"))
    beta = config["beta_0_um"]
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    print(figure_d(config, beta))
    print(figure_v(config))
    s, e1, e60, s_star = d_lag_curves(beta)
    print(f"D: 60 s breakaway at {s_star:.3f} um travel; 1 s lag at 5 um {e1[-1]:+.3f} um")


if __name__ == "__main__":
    main()
