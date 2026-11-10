#!/usr/bin/env python3
"""F1 Q block as a plain-language GIF: does speed matter?

The brief is F1_M_Q_gifs_simple.md. It replaces the technical Q1 figure
(coincident loops, flat width) with one animation aimed at a reader with no
background: the same small back-and-forth is done fast, medium and slow, and
the stage draws the same loop every time.

    python3 render_f1_Q_gif.py      # writes rendered/companion/f1_Q_speed.gif

No measured data. The response comes from the same 30-element quasi-static
Maxwell-slip solver the other companion figures use, imported from
render_f1_companions so the three GIFs and the remaining figure cannot drift
apart.

The solver has no rate in it at all, so the three loops are not merely close,
they are the same numbers; main() asserts that before the GIF is written.
That is the honest form of the claim: F1's Q block is the measurement that
would show a rate effect if there were one, and this is what the model being
tested predicts.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_f1_companions import (LIMITS, G, run, path as fine_path, BETA,
                                  S90_CFG, width_at_centre)
from f1_gif_common import (DOTTED, FPS, GHOST, INK, MUTED, ORANGE, TEAL,
                           capture, footer, lighten, new_figure, report,
                           save_gif)

HERE = Path(__file__).resolve().parent
GIF = HERE / "rendered" / "companion" / "f1_Q_speed.gif"

STROKE = 3.0                      # brief: around the size where the loop opens
HALF = STROKE / 2
CYCLES = 2                        # back-and-forths per pass

# name, the real half-cycle time, frames per half, line width, how pale.
# The drawing stretches time 1 : 2 : 4; the labels carry the real 1 : 3 : 9.
PASSES = [
    ("fast", "1 s", 10, 7.0, 0.58),
    ("medium", "4 s", 20, 3.6, 0.30),
    ("slow", "12 s", 40, 1.7, 0.00),
]
RATE_EFFECT = 0.10                # the grey alternative: +10% from fast to slow

INTRO, PAUSE, WHATIF, FINAL = 40, 20, 50, 60

MESSAGES = [
    "The same small move, done at three speeds.",
    "Fast: this is the loop.",
    "Medium speed: the same loop.",
    "Slow: still the same loop.",
    "If speed mattered, the gap would change like this.",
    "Speed doesn't change what the stage does.",
]


def to_x(u):
    return (np.asarray(u) + HALF) / STROKE


# -- the motion -----------------------------------------------------------
def build() -> dict:
    """Three passes of the same back-and-forth, run one after another through
    one solver, after one settling cycle so every drawn loop is the
    repeatable one."""
    state = LIMITS.copy()                       # wiped: the last motion was up
    _, _, state = run(state, fine_path([-HALF, HALF, -HALF], 240))

    u, bounds = [-HALF], []
    for _name, _time, fph, _lw, _pale in PASSES:
        start = len(u) - 1
        for _ in range(CYCLES):
            u += list(np.linspace(-HALF, HALF, fph + 1)[1:])
            u += list(np.linspace(HALF, -HALF, fph + 1)[1:])
        bounds.append((start, len(u) - 1))
    u = np.asarray(u)
    _, lag, _ = run(state, u)
    return {"u": u, "lag": lag, "x": to_x(u), "bounds": bounds}


def half_slice(d: dict, p: int, index: int) -> slice:
    """Samples of one half-cycle of pass p, turning points at both ends."""
    fph = PASSES[p][2]
    start = d["bounds"][p][0] + index * fph
    return slice(start, start + fph + 1)


def gap(d: dict, p: int, up: int = 2, down: int = 3) -> float:
    """How far apart the two branches are halfway through the move; measured
    on the second back-and-forth, as in the P GIF."""
    a, b = half_slice(d, p, up), half_slice(d, p, down)
    out = np.interp(0.5, d["x"][a], d["lag"][a])
    back = np.interp(0.5, d["x"][b][::-1], d["lag"][b][::-1])
    return float(back - out)


def branches(d: dict, p: int, k: int) -> tuple[list, list, list, list]:
    """The up and down polylines of pass p as far as global sample k."""
    fph = PASSES[p][2]
    start = d["bounds"][p][0]
    up_x, up_y, down_x, down_y = [], [], [], []
    for index in range(2 * CYCLES):
        a = start + index * fph
        if k <= a:
            break
        b = min(k, a + fph)
        target = (up_x, up_y) if index % 2 == 0 else (down_x, down_y)
        target[0].extend([*d["x"][a:b + 1], np.nan])
        target[1].extend([*d["lag"][a:b + 1], np.nan])
    return up_x, up_y, down_x, down_y


# -- figure ---------------------------------------------------------------
def build_figure(d: dict) -> dict:
    fig = new_figure()
    motor = fig.add_axes([0.075, 0.555, 0.385, 0.300])
    stage = fig.add_axes([0.580, 0.555, 0.385, 0.300])
    panel = fig.add_axes([0.075, 0.205, 0.890, 0.190])
    art = {"fig": fig, "motor": motor, "stage": stage, "panel": panel}
    last = d["bounds"][-1][1]

    # (A) what the motor is told to do: the same move three times over
    motor.set(xlim=(-0.03, 1.03), ylim=(-1.50, 1.50), xticks=[], yticks=[])
    motor.set_xlabel("time", labelpad=2)
    motor.set_ylabel("motor position", labelpad=6)
    motor.set_title("What the motor is told to do", loc="left", pad=20)
    motor.spines[["left", "bottom"]].set_color(GHOST)
    for p, (name, seconds, _fph, _lw, _pale) in enumerate(PASSES):
        a, b = (v / last for v in d["bounds"][p])
        if p:
            motor.axvline(a, color=GHOST, lw=0.8, zorder=0)
        motor.text((a + b) / 2, 1.47, f"{name}\n{seconds}", ha="center",
                   va="top", fontsize=10, color=MUTED, linespacing=1.15)
    art["plan"], = motor.plot([], [], color=GHOST, lw=2.2, solid_capstyle="round")
    art["trail"], = motor.plot([], [], color=INK, lw=2.2, solid_capstyle="round")
    art["dot"], = motor.plot([], [], "o", color=INK, ms=8)

    # (B) what the stage does: each pass traces over the one before it
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
    stage.text(0.33, 1.03, "coming back", transform=stage.transAxes,
               color=ORANGE, fontsize=11, va="bottom")
    art["loops"] = []
    for p, (_name, _s, _fph, lw, pale) in enumerate(PASSES):
        art["loops"].append((
            stage.plot([], [], color=lighten(TEAL, pale), lw=lw, zorder=3 + p,
                       solid_capstyle="round")[0],
            stage.plot([], [], color=lighten(ORANGE, pale), lw=lw, zorder=3 + p,
                       solid_capstyle="round")[0]))
    # Outside the lens the loops make, in the empty upper right.
    art["same"] = stage.text(0.80, 0.56, "", fontsize=11, color=MUTED,
                             ha="center", va="center", zorder=9)
    art["gap_arrow"] = None
    art["gap_text"] = stage.text(0.545, 0.0, "", fontsize=11, color=INK,
                                 ha="left", va="center", zorder=9)
    art["live"], = stage.plot([], [], "o", color=INK, ms=7, zorder=10)

    # (C) the gap each pass produced
    g = gap(d, 0)
    # Zero-based, so a 10% rate effect looks like the 10% it is.
    panel.set(xlim=(-0.45, 2.75), ylim=(0, 1.25), yticks=[],
              xticks=range(len(PASSES)))
    panel.set_xticklabels([f"{n}\n{s}" for n, s, *_ in PASSES])
    panel.tick_params(axis="x", length=0, pad=6)
    panel.set_ylabel("gap", labelpad=6)
    panel.set_title("The gap is the same at every speed", loc="left", pad=8)
    panel.spines[["left", "bottom"]].set_color(GHOST)
    art["flat"], = panel.plot([], [], color=TEAL, lw=2.5, zorder=2)
    art["dots"] = [panel.plot([], [], "o", color=INK, ms=9, zorder=3)[0]
                   for _ in PASSES]
    art["slope"], = panel.plot([], [], color=DOTTED, lw=2.0, ls=":", zorder=1)
    art["slope_tag"] = panel.text(2.08, g * (1 + RATE_EFFECT),
                                  "  if speed mattered", fontsize=10,
                                  color=DOTTED, ha="left", va="center",
                                  visible=False)

    art["message"] = fig.text(0.075, 0.080, "", fontsize=13, color=INK,
                              va="center")
    footer(fig, "each speed is repeated until the loop settles; the drawing "
                "stretches time 1 : 2 : 4, the real passes are 1 : 3 : 9")
    return art


def draw_gap(art: dict, d: dict) -> None:
    """The one arrow that says what 'gap' means, kept for the whole GIF."""
    a, b = half_slice(d, 0, 2), half_slice(d, 0, 3)
    low = float(np.interp(0.5, d["x"][a], d["lag"][a]))
    high = float(np.interp(0.5, d["x"][b][::-1], d["lag"][b][::-1]))
    art["gap_arrow"] = art["stage"].annotate(
        "", xy=(0.5, low), xytext=(0.5, high),
        arrowprops=dict(arrowstyle="<->", color=INK, lw=1.4,
                        mutation_scale=16), zorder=8)
    art["gap_text"].set_text("gap")


# -- frames ---------------------------------------------------------------
def render(d: dict) -> list:
    art = build_figure(d)
    fig = art["fig"]
    last = d["bounds"][-1][1]
    plan_t = np.linspace(0, 1, last + 1)
    plan_y = d["u"] / HALF
    frames = []

    art["message"].set_text(MESSAGES[0])
    art["plan"].set_data(plan_t, plan_y)
    frames += [capture(fig)] * INTRO

    for p, (_name, _s, _fph, _lw, _pale) in enumerate(PASSES):
        art["message"].set_text(MESSAGES[p + 1])
        start, stop = d["bounds"][p]
        for k in range(start + 1, stop + 1):
            up_x, up_y, down_x, down_y = branches(d, p, k)
            art["loops"][p][0].set_data(up_x, up_y)
            art["loops"][p][1].set_data(down_x, down_y)
            art["trail"].set_data(plan_t[:k + 1], plan_y[:k + 1])
            art["dot"].set_data([plan_t[k]], [plan_y[k]])
            art["live"].set_data([d["x"][k]], [d["lag"][k]])
            frames.append(capture(fig))

        if p == 0:
            draw_gap(art, d)
        if p:
            art["same"].set_text("same loop")
        art["dots"][p].set_data([p], [gap(d, p)])
        if p:
            art["flat"].set_data(range(p + 1), [gap(d, i) for i in range(p + 1)])
        frames += [capture(fig)] * PAUSE

    # the grey alternative
    art["message"].set_text(MESSAGES[4])
    art["live"].set_data([], [])
    g = gap(d, 0)
    art["slope"].set_data([0, len(PASSES) - 1], [g, g * (1 + RATE_EFFECT)])
    art["slope_tag"].set_visible(True)
    frames += [capture(fig)] * WHATIF

    art["message"].set_text(MESSAGES[5])
    frames += [capture(fig)] * FINAL

    import matplotlib.pyplot as plt
    plt.close(fig)
    return frames


def main() -> int:
    d = build()

    # The claim the GIF makes: not "close", the same. Every pass samples the
    # stroke at a divisor of the slow pass's spacing, so the three can be
    # compared at coincident commanded positions with no interpolation.
    gaps = [gap(d, p) for p in range(len(PASSES))]
    for p in range(1, len(PASSES)):
        assert abs(gaps[p] - gaps[0]) < 1e-12, (p, gaps)
    assert abs(gaps[0] - width_at_centre(STROKE)) < 1e-9, gaps[0]
    slow = len(PASSES) - 1
    exact = 0.0
    for p in range(slow):
        step = PASSES[slow][2] // PASSES[p][2]
        assert PASSES[slow][2] == step * PASSES[p][2]
        for index in range(2 * CYCLES):
            mine = d["lag"][half_slice(d, p, index)]
            theirs = d["lag"][half_slice(d, slow, index)][::step]
            assert np.allclose(d["u"][half_slice(d, p, index)],
                               d["u"][half_slice(d, slow, index)][::step])
            exact = max(exact, float(np.abs(mine - theirs).max()))
    assert exact < 1e-12, exact

    # What is left is only how coarsely each pass is drawn: the fast pass
    # joins 10 points per half-cycle, so it cuts the corner at a reversal by
    # well under a pixel. Reported, not hidden.
    chord = 0.0
    ref = half_slice(d, slow, 2)
    for p in range(slow):
        a = half_slice(d, p, 2)
        here = np.interp(d["x"][ref], d["x"][a], d["lag"][a])
        chord = max(chord, float(np.abs(here - d["lag"][ref]).max()))
    assert chord < 0.01, chord

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)   # catch missing glyphs
        frames = render(d)
    save_gif(frames, GIF)

    report(GIF, frames)
    print(f"  stroke {STROKE:g} µm, {CYCLES} back-and-forths per speed")
    for (name, seconds, fph, *_), g in zip(PASSES, gaps):
        print(f"  {name:7s} {seconds:>4s}  {fph:3d} frames per half  "
              f"gap {g:.6f} µm")
    print(f"  the three gaps agree to {max(gaps) - min(gaps):.1e} µm and the "
          f"loops to {exact:.1e} µm at shared positions")
    print(f"  the coarse fast pass cuts the reversal corner by {chord * 1000:.1f} nm, "
          f"about {chord / 2.1 * 260:.1f} px")
    print(f"  maximum gap (saturation) {BETA:g} µm, illustrative s90 "
          f"{S90_CFG:g} µm, {len(LIMITS)} elements, stiffness {G.sum():.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
