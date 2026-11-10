#!/usr/bin/env python3
"""F1 M block as a plain-language GIF: does the stage remember where it was?

The brief is F1_M_Q_gifs_simple.md. It replaces the technical M1 figure
(nested loop, return-point detail) with one animation aimed at a reader with
no background: the stage makes a big back-and-forth, halfway back it takes a
small detour, and afterwards it picks up exactly where it left off.

    python3 render_f1_M_gif.py     # writes rendered/companion/f1_M_memory.gif

No measured data. The response comes from the same 30-element quasi-static
Maxwell-slip solver the other companion figures use, imported from
render_f1_companions so the three GIFs and the remaining figure cannot drift
apart.

THE GREY ALTERNATIVE. "A stage that forgets" is the same solver with the
outer excursion erased: after the detour's turning point it responds as a
stage whose whole history is that detour, so it comes back onto the big path
0.09 um low and closes the difference only over the next few microns. It is
a worked alternative rather than a drawn squiggle, but it is still an
illustration of one way to lack memory, not a claim about every model.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_f1_companions import LIMITS, G, run, path as fine_path, BETA, S90_CFG
from f1_gif_common import (DOTTED, FPS, GHOST, INK, MUTED, ORANGE, TEAL,
                           capture, footer, new_figure, report, save_gif)

HERE = Path(__file__).resolve().parent
GIF = HERE / "rendered" / "companion" / "f1_M_memory.gif"

# Geometry, brief section "Curves": big move 10 um, detour 3 um from halfway.
OUTER, DETOUR = 10.0, 3.0
HALF = OUTER / 2

# Pace, brief section "Timing": 20 fps, about 20 s.
INTRO, F_UP, F_DOWN, F_DET, F_CONT, WHATIF, HOLD = 40, 60, 40, 80, 50, 60, 60
K_UP, K_DOWN = F_UP, F_UP + F_DOWN
K_DET, K_CONT = K_DOWN + F_DET, K_DOWN + F_DET + F_CONT

MESSAGES = {
    0: "The stage makes a big move back and forth, with a small detour on the way back.",
    1: "Going up: the stage follows, a little behind.",
    2: "Coming back down…",
    3: "A small detour: up a little, then back.",
    4: "After the detour, it continues exactly where it left off.",
    5: "A stage that forgets would take a different path.",
    6: "The stage remembers where it was.",
}

# The zoom window in panel (C), in microns of commanded position and of lag.
# Panel (B) spans 2.1 um of lag; this window is the magnification of it.
ZOOM_U = (-5.4, 3.0)
ZOOM_E = (0.44, 0.90)
ZOOM_FACTOR = 2.1 / (ZOOM_E[1] - ZOOM_E[0])


def to_x(u):
    """Panel (B) and (C) share one horizontal coordinate: 0 at the bottom of
    the big move, 1 at the top."""
    return (np.asarray(u) + HALF) / OUTER


# -- the motion -----------------------------------------------------------
def build() -> dict:
    """The drawn path, the path it would have taken without the detour, and
    the path a stage without memory would take after the detour.

    SETTLING, as in the P GIF: one full big move is run through the solver
    before the drawn one, so the drawn loop is the repeatable one and the
    reversal at the bottom shows its bend rather than a flat saturated line.
    """
    state = LIMITS.copy()                       # wiped: the last motion was up
    _, _, state = run(state, fine_path([-HALF, HALF, -HALF], 240))

    u = np.concatenate([
        [-HALF],
        np.linspace(-HALF, HALF, F_UP + 1)[1:],
        np.linspace(HALF, 0.0, F_DOWN + 1)[1:],
        np.linspace(0.0, DETOUR, F_DET // 2 + 1)[1:],
        np.linspace(DETOUR, 0.0, F_DET // 2 + 1)[1:],
        np.linspace(0.0, -HALF, F_CONT + 1)[1:],
    ])
    assert len(u) == K_CONT + 1
    _, lag, _ = run(state, u)

    # What it would have done without the detour, from the same settled state.
    ref_u, ref_e, _ = run(state, fine_path([-HALF, HALF, -HALF], 600))
    ref_u, ref_e = ref_u[600:], ref_e[600:]          # the down branch only

    # A stage that forgets: the same solver with the outer excursion erased,
    # driven down from the detour's turning point with no earlier history.
    f_u, f_e, _ = run(LIMITS.copy(), fine_path([DETOUR, -HALF], 600))
    keep = f_u <= 0.0
    return {"u": u, "lag": lag, "ref_u": ref_u, "ref_e": ref_e,
            "forget_u": f_u[keep], "forget_e": f_e[keep],
            "branch": float(lag[K_DOWN])}


def phase_of(k: int) -> int:
    if k <= 0:
        return 0
    return 1 if k <= K_UP else 2 if k <= K_DOWN else 3 if k <= K_DET else 4


# -- figure ---------------------------------------------------------------
def build_figure(d: dict) -> dict:
    fig = new_figure()
    motor = fig.add_axes([0.075, 0.555, 0.385, 0.300])
    stage = fig.add_axes([0.580, 0.555, 0.385, 0.300])
    zoom = fig.add_axes([0.075, 0.215, 0.890, 0.170])
    art = {"fig": fig, "motor": motor, "stage": stage, "zoom": zoom}

    # (A) what the motor is told to do
    motor.set(xlim=(-0.03, 1.03), ylim=(-1.45, 1.45), xticks=[], yticks=[])
    motor.set_xlabel("time", labelpad=2)
    motor.set_ylabel("motor position", labelpad=6)
    motor.set_title("What the motor is told to do", loc="left", pad=20)
    motor.spines[["left", "bottom"]].set_color(GHOST)
    motor.axvspan(K_DOWN / K_CONT, K_DET / K_CONT, color="#eceff1", lw=0, zorder=0)
    motor.text((K_DOWN + K_DET) / 2 / K_CONT, 1.42, "small detour", ha="center",
               va="top", fontsize=10, color=MUTED)
    art["plan"], = motor.plot([], [], color=GHOST, lw=2.5, solid_capstyle="round")
    art["trail"], = motor.plot([], [], color=INK, lw=2.5, solid_capstyle="round")
    art["dot"], = motor.plot([], [], "o", color=INK, ms=8)
    art["up"], = motor.plot([], [], "^", color=INK, ms=7)
    art["down"], = motor.plot([], [], "v", color=INK, ms=7)

    # (B) what the stage does
    stage.set(xlim=(-0.05, 1.05), ylim=(-1.05, 1.05), xticks=[0, 1],
              yticks=[-BETA / 2, BETA / 2])
    stage.set_xticklabels(["start", "end"])
    stage.set_yticklabels(["behind\n(going up)", "behind\n(going down)"],
                          fontsize=10)
    stage.get_yticklabels()[0].set_color(TEAL)
    stage.get_yticklabels()[1].set_color(ORANGE)
    stage.set_xlabel("where the motor is", labelpad=2)
    stage.set_title("What the stage does", loc="left", pad=20)
    stage.spines[["left", "bottom"]].set_color(GHOST)
    stage.text(0.00, 1.03, "going up", transform=stage.transAxes, color=TEAL,
               fontsize=11, va="bottom")
    stage.text(0.30, 1.03, "coming back", transform=stage.transAxes,
               color=ORANGE, fontsize=11, va="bottom")
    stage.text(0.72, 1.03, "detour", transform=stage.transAxes, color=INK,
               fontsize=11, va="bottom")

    art["heading"], = stage.plot([], [], color=GHOST, lw=1.8, zorder=1,
                                 solid_capstyle="round")
    # The faint line runs along the top of the panel, so its label sits in the
    # empty middle-left with a leader rather than on top of the line.
    art["heading_tag"] = stage.annotate(
        "where it was heading", xy=(0.25, float(np.interp(
            -2.5, d["ref_u"][::-1], d["ref_e"][::-1]))),
        xytext=(0.09, 0.40), fontsize=10, color=MUTED, ha="left", va="center",
        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8), visible=False)
    art["seg"] = {
        "up": stage.plot([], [], color=TEAL, lw=3, solid_capstyle="round",
                         zorder=3)[0],
        "down": stage.plot([], [], color=ORANGE, lw=3, solid_capstyle="round",
                           zorder=3)[0],
        "det": stage.plot([], [], color=INK, lw=3.6, solid_capstyle="round",
                          zorder=5)[0],
        "cont": stage.plot([], [], color=ORANGE, lw=3, solid_capstyle="round",
                           zorder=3)[0],
    }
    art["mark"], = stage.plot([], [], "o", mfc="white", mec=INK, mew=1.6, ms=9,
                              zorder=6)
    # Beside the circle, on a white patch: the faint line runs right past it.
    art["mark_tag"] = stage.text(0, 0, "", fontsize=10, color=INK, ha="right",
                                 va="center", visible=False, zorder=7,
                                 bbox=dict(fc="white", ec="none", pad=1.0))
    art["live"], = stage.plot([], [], "o", color=INK, ms=7, zorder=8)
    art["arrow"] = None

    # (C) the zoom, hidden until the detour starts
    zoom.set(xlim=to_x(ZOOM_U), ylim=ZOOM_E, xticks=[to_x(0.0)], yticks=[])
    zoom.set_xticklabels(["where the detour starts"])
    zoom.set_xlabel("where the motor is", labelpad=11)
    zoom.set_ylabel("how far behind", labelpad=6)
    zoom.set_title("Zoom on the moment it rejoins  ·  "
                   f"vertical scale ×{ZOOM_FACTOR:.0f}", loc="left", pad=8)
    zoom.spines[["left", "bottom"]].set_color(GHOST)
    zoom.plot(to_x(d["ref_u"]), d["ref_e"], color=GHOST, lw=2.0, zorder=1)
    art["z_before"], = zoom.plot([], [], color=ORANGE, lw=3, zorder=3,
                                 solid_capstyle="round")
    art["z_det"], = zoom.plot([], [], color=INK, lw=3.6, zorder=4,
                              solid_capstyle="round")
    art["z_cont"], = zoom.plot([], [], color=ORANGE, lw=3, zorder=3,
                               solid_capstyle="round")
    art["z_forget"], = zoom.plot([], [], color=DOTTED, lw=2.2, ls=":", zorder=2)
    art["z_forget_tag"] = zoom.annotate(
        "a stage that forgets would take this path",
        xy=(to_x(-1.1), float(np.interp(-1.1, d["forget_u"][::-1],
                                        d["forget_e"][::-1]))),
        xytext=(to_x(-3.0), ZOOM_E[0] + 0.07), fontsize=10, color=DOTTED,
        ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color=DOTTED, lw=0.8), visible=False)
    art["z_mark"], = zoom.plot([], [], "o", mfc="white", mec=INK, mew=1.6,
                               ms=10, zorder=6)
    art["z_tag"] = zoom.text(to_x(0.35), d["branch"] + 0.055, "", fontsize=10,
                             color=INK, ha="left", va="bottom", visible=False)
    zoom.set_visible(False)

    art["message"] = fig.text(0.075, 0.080, "", fontsize=13, color=INK,
                              va="center")
    footer(fig, "the big move is repeated until the loop settles; "
                "the detour is drawn larger than the real test")
    return art


def set_arrow(art: dict, xs, ys, colour: str) -> None:
    """An arrowhead riding the branch being drawn, pointing the way it goes."""
    if art["arrow"] is not None:
        art["arrow"].remove()
        art["arrow"] = None
    if len(xs) < 9:
        return
    mid = len(xs) // 2
    art["arrow"] = art["stage"].annotate(
        "", xy=(xs[mid + 3], ys[mid + 3]), xytext=(xs[mid - 3], ys[mid - 3]),
        arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.0,
                        mutation_scale=17), zorder=7)


# -- frames ---------------------------------------------------------------
def render(d: dict) -> list:
    art = build_figure(d)
    fig = art["fig"]
    u, lag = d["u"], d["lag"]
    x = to_x(u)
    plan_t = np.linspace(0, 1, len(u))
    plan_y = u / HALF
    frames = []

    # phase 0: the opening line, nothing drawn yet
    art["message"].set_text(MESSAGES[0])
    art["plan"].set_data(plan_t, plan_y)
    frames += [capture(fig)] * INTRO

    spans = {"up": (0, K_UP), "down": (K_UP, K_DOWN),
             "det": (K_DOWN, K_DET), "cont": (K_DET, K_CONT)}
    colours = {"up": TEAL, "down": ORANGE, "det": INK, "cont": ORANGE}
    phase = 0
    for k in range(1, K_CONT + 1):
        new_phase = phase_of(k)
        if new_phase != phase:
            phase = new_phase
            art["message"].set_text(MESSAGES[phase])
            if phase == 2:
                art["heading"].set_data(to_x(d["ref_u"]), d["ref_e"])
                art["heading_tag"].set_visible(True)
            if phase == 3:
                art["zoom"].set_visible(True)
                art["mark"].set_data([x[K_DOWN]], [lag[K_DOWN]])
                art["z_mark"].set_data([x[K_DOWN]], [lag[K_DOWN]])
                art["mark_tag"].set_position((x[K_DOWN] - 0.03, lag[K_DOWN]))
                art["mark_tag"].set_text("leaves here")
                art["mark_tag"].set_visible(True)
                art["z_tag"].set_text("leaves here")
                art["z_tag"].set_visible(True)
            if phase == 4:
                art["mark_tag"].set_text("back on track")
                art["z_tag"].set_text("back on track")

        for name, (a, b) in spans.items():
            stop = min(k, b)
            if stop <= a:
                art["seg"][name].set_data([], [])
            else:
                art["seg"][name].set_data(x[a:stop + 1], lag[a:stop + 1])
        art["trail"].set_data(plan_t[:k + 1], plan_y[:k + 1])
        art["dot"].set_data([plan_t[k]], [plan_y[k]])
        art["live"].set_data([x[k]], [lag[k]])
        rising = u[k] > u[k - 1]
        art["up"].set_data([plan_t[k]] if rising else [],
                           [plan_y[k] + 0.30] if rising else [])
        art["down"].set_data([] if rising else [plan_t[k]],
                             [] if rising else [plan_y[k] - 0.30])
        a, b = spans["up" if phase == 1 else "down" if phase == 2 else
                     "det" if phase == 3 else "cont"]
        set_arrow(art, x[a:min(k, b) + 1], lag[a:min(k, b) + 1],
                  colours["up" if phase == 1 else "down" if phase == 2 else
                          "det" if phase == 3 else "cont"])

        # the zoom follows the same samples, once it is showing
        if phase >= 3:
            art["z_before"].set_data(x[K_UP:K_DOWN + 1], lag[K_UP:K_DOWN + 1])
            art["z_det"].set_data(x[K_DOWN:min(k, K_DET) + 1],
                                  lag[K_DOWN:min(k, K_DET) + 1])
            if k > K_DET:
                art["z_cont"].set_data(x[K_DET:k + 1], lag[K_DET:k + 1])
        frames.append(capture(fig))

    # phase 5: the alternative
    art["message"].set_text(MESSAGES[5])
    art["live"].set_data([], [])
    if art["arrow"] is not None:
        art["arrow"].remove()
        art["arrow"] = None
    art["z_forget"].set_data(to_x(d["forget_u"]), d["forget_e"])
    art["z_forget_tag"].set_visible(True)
    frames += [capture(fig)] * WHATIF

    # phase 6: hold on the finished picture
    art["message"].set_text(MESSAGES[6])
    frames += [capture(fig)] * HOLD

    import matplotlib.pyplot as plt
    plt.close(fig)
    return frames


def main() -> int:
    d = build()
    u, lag = d["u"], d["lag"]

    # The claim the GIF makes, checked rather than drawn: the detour leaves
    # and rejoins at the same point, and what follows is the branch the big
    # move was on, not merely something close to it.
    assert abs(u[K_DOWN]) < 1e-12 and abs(u[K_DET]) < 1e-12
    closes = abs(lag[K_DET] - lag[K_DOWN])
    assert closes < 1e-9, closes
    after = np.interp(u[K_DET:][::-1], d["ref_u"][::-1], d["ref_e"][::-1])[::-1]
    rejoin = float(np.abs(lag[K_DET:] - after).max())
    assert rejoin < 1e-4, rejoin

    # The grey alternative has to miss by enough to be worth drawing, and
    # then close the difference rather than stay off.
    miss = float(lag[K_DOWN] - np.interp(0.0, d["forget_u"][::-1],
                                         d["forget_e"][::-1]))
    end = float(np.interp(-HALF, d["forget_u"][::-1], d["forget_e"][::-1]))
    assert 0.02 < miss < 0.5, miss
    assert abs(lag[-1] - end) < miss / 2, (lag[-1], end)

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)   # catch missing glyphs
        frames = render(d)
    save_gif(frames, GIF)

    report(GIF, frames)
    print(f"  big move {OUTER:g} µm, detour {DETOUR:g} µm from halfway down")
    print(f"  leaves at {lag[K_DOWN]:+.4f} µm behind, returns to "
          f"{lag[K_DET]:+.4f} µm: closes to {closes:.1e} µm")
    print(f"  after the detour it follows the uninterrupted branch to "
          f"{rejoin:.1e} µm")
    print(f"  a stage that forgets comes back {miss:.3f} µm low and is still "
          f"{abs(lag[-1] - end):.3f} µm off at the end of the move")
    print(f"  maximum gap (saturation) {BETA:g} µm, illustrative s90 "
          f"{S90_CFG:g} µm, {len(LIMITS)} elements, stiffness {G.sum():.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
