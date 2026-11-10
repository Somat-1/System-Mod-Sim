#!/usr/bin/env python3
"""F1 P block as a plain-language GIF: back and forth, and what the stage does.

The brief is F1_P_loops_gif_simple.md. It replaces the two technical P
figures (P1 loops-and-width, P2 Masing overlay) with one animation aimed at
a reader with no background: the motor is moved back and forth by three
different amounts, the stage follows a little behind, and drawing how far
behind against where the motor is gives a loop that opens up with the size
of the move and then stops opening.

    python3 render_f1_P_gif.py      # writes rendered/companion/f1_P_loops.gif

No measured data. The response comes from the same 30-element quasi-static
Maxwell-slip solver the other companion figures use, imported from
render_f1_companions so that the GIF and the remaining figures cannot drift
apart: saturation beta_0 = 1.757 um, initial stiffness 20, and the
placeholder s90 = 3 um from f1_config.json.

The brief sketches its own small solver in section 5 and calls it a
placeholder; the repository's calibrated one is used instead. Both produce
the same picture, and the calibrated one is the one already checked by
check_f1_companions.py.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_f1_companions import LIMITS, G, run, width_at_centre, BETA, S90_CFG
from f1_gif_common import (FPS, GHOST, INK, MUTED, ORANGE, SIZE_PX, TEAL,
                           capture, footer, new_figure, report, save_gif)

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered" / "companion"
GIF = OUT / "f1_P_loops.gif"

# Pace from the brief, section 4; colours and canvas from f1_gif_common.
FPH = 30            # frames per half-cycle (up, or back down)
FPA = 10            # frames to get into position before each move
CYCLES = 2          # back-and-forths per move
HOLD = 20           # frames held after each move, while its point drops in
INTRO = 40
SUMMARY = 60
TAIL = 20           # the 1 s pause before the loop restarts
STROKES = [1.0, 3.0, 100.0]

MESSAGES = [
    "The motor moves the stage back and forth. Watch how far the stage falls behind.",
    "Small move: the stage hardly moves, so the loop stays thin.",
    "Medium move: the stage starts to follow, so the loop opens.",
    "Large move: the stage follows, but always a bit behind. "
    "The gap is now as big as it gets.",
    "The gap grows with the size of the move, then stops growing.",
]
LABELS = ["small move", "medium move", "large move"]


# -- the three moves ------------------------------------------------------
def build_moves() -> list[dict]:
    """Three back-and-forth moves, run one after another through one solver.

    The state carries across the moves, exactly as it does on the bench: the
    P block wipes once, then walks up the stroke grid. Each move starts by
    travelling to its own lower turning point, which is drawn but is not
    part of the loop.

    SETTLING. The first back-and-forth of a move inherits the state the
    previous one left behind, so its loop is not yet the repeatable one: for
    the 1 um move the first cycle is more than ten times as open as the
    second. The bench handles this by running three cycles per stroke and
    reading the later ones. The GIF does the same thing: one cycle is run
    through the solver before the two that are drawn, so both shown cycles
    are the settled loop and the gap does not jump halfway through. The
    footer says so, and main() checks that the two drawn cycles agree.
    """
    state = LIMITS.copy()            # wiped: the last motion was upward
    position, moves = 0.0, []
    for stroke in STROKES:
        half = stroke / 2
        cycle = (list(np.linspace(-half, half, FPH + 1)[1:])
                 + list(np.linspace(half, -half, FPH + 1)[1:]))
        lead = ([position] + list(np.linspace(position, -half, FPA + 1)[1:]))
        _, lead_lag, state = run(state, np.asarray(lead + cycle))
        shown = np.asarray([-half] + cycle * CYCLES)
        _, shown_lag, state = run(state, shown)
        u = np.asarray(lead + list(shown[1:]))
        lag = np.r_[np.full(FPA, np.nan), shown_lag]   # the lead is not drawn
        assert len(u) == len(lag) == 1 + FPA + 2 * CYCLES * FPH
        moves.append({
            "s": stroke,
            "u": u,
            "lag": lag,
            "x": (u + half) / stroke,          # 0 at the start, 1 at the end
            "n": u / half,                     # normalised command, for panel A
            "lead_lag": lead_lag[:FPA + 1],
        })
        position = -half
    for move in moves:
        move["gap"] = gap_between(move, 2, 3)  # the second back-and-forth
    return moves


def half_slice(index: int) -> slice:
    """Samples of one half-cycle, turning points included at both ends."""
    start = FPA + index * FPH
    return slice(start, start + FPH + 1)


def gap_between(move: dict, up: int, down: int) -> float:
    """How far apart the two branches are halfway through the move."""
    a, b = half_slice(up), half_slice(down)
    out = np.interp(0.5, move["x"][a], move["lag"][a])
    back = np.interp(0.5, move["x"][b][::-1], move["lag"][b][::-1])
    return float(back - out)


def branches(move: dict, k: int) -> tuple[list, list, list, list]:
    """The teal and orange polylines as far as sample k, NaN-separated."""
    up_x, up_y, down_x, down_y = [], [], [], []
    for index in range(2 * CYCLES):
        start = FPA + index * FPH
        if k <= start:
            break
        stop = min(k, start + FPH)
        xs, ys = move["x"][start:stop + 1], move["lag"][start:stop + 1]
        target = (up_x, up_y) if index % 2 == 0 else (down_x, down_y)
        target[0].extend([*xs, np.nan])
        target[1].extend([*ys, np.nan])
    return up_x, up_y, down_x, down_y


# -- figure ---------------------------------------------------------------
def build_figure(moves: list[dict]):
    fig = new_figure()
    motor = fig.add_axes([0.075, 0.520, 0.385, 0.325])
    stage = fig.add_axes([0.580, 0.520, 0.385, 0.325])
    summary = fig.add_axes([0.075, 0.175, 0.890, 0.205])
    art: dict = {"fig": fig, "motor": motor, "stage": stage, "summary": summary}

    # (A) what the motor is told to do
    motor.set(xlim=(-0.03, 1.03), ylim=(-1.45, 1.35), xticks=[], yticks=[])
    motor.set_xlabel("time", labelpad=2)
    motor.set_ylabel("motor position", labelpad=6)
    motor.set_title("What the motor is told to do", loc="left", pad=22)
    motor.spines[["left", "bottom"]].set_color(GHOST)
    art["plan"], = motor.plot([], [], color=GHOST, lw=2.5, solid_capstyle="round")
    art["trail"], = motor.plot([], [], color=INK, lw=2.5, solid_capstyle="round")
    art["dot"], = motor.plot([], [], "o", color=INK, ms=8)
    art["up"], = motor.plot([], [], "^", color=INK, ms=7)
    art["down"], = motor.plot([], [], "v", color=INK, ms=7)
    art["size_tag"] = motor.text(0.0, 1.03, "", transform=motor.transAxes,
                                 ha="left", va="bottom", fontsize=11,
                                 color=MUTED)
    art["faint_plans"] = []
    art["faint_tags"] = []

    # (B) what the stage does. The two y labels sit where the stage is
    # furthest behind in each direction, so no rotated axis text is needed.
    stage.set(xlim=(-0.05, 1.05), ylim=(-1.05, 1.05), xticks=[0, 1],
              yticks=[-BETA / 2, BETA / 2])
    stage.set_xticklabels(["start", "end"])
    stage.set_yticklabels(["behind\n(going up)", "behind\n(going down)"],
                          fontsize=10)
    stage.get_yticklabels()[0].set_color(TEAL)
    stage.get_yticklabels()[1].set_color(ORANGE)
    stage.set_xlabel("where the motor is", labelpad=2)
    stage.set_title("What the stage does", loc="left", pad=22)
    stage.spines[["left", "bottom"]].set_color(GHOST)
    stage.text(0.0, 1.03, "going up", transform=stage.transAxes, color=TEAL,
               fontsize=11, va="bottom")
    stage.text(0.33, 1.03, "coming back", transform=stage.transAxes,
               color=ORANGE, fontsize=11, va="bottom")
    art["old_loops"] = []
    art["branch_up"], = stage.plot([], [], color=TEAL, lw=3,
                                   solid_capstyle="round")
    art["branch_down"], = stage.plot([], [], color=ORANGE, lw=3,
                                     solid_capstyle="round")
    art["live"], = stage.plot([], [], "o", color=INK, ms=7)
    art["arrow_up"] = None
    art["arrow_down"] = None
    art["gap_arrow"] = None
    art["gap_text"] = stage.text(0, 0, "", fontsize=11, color=INK,
                                 ha="center", va="bottom", visible=False,
                                 zorder=8)

    # (C) loop size for each move size
    biggest = max(STROKES)
    summary.set(xscale="log", xlim=(0.45, biggest * 2.6), ylim=(0, 2.25),
                yticks=[])
    summary.set_xticks(STROKES)
    summary.set_xticklabels([f"{s:g} µm" for s in STROKES])
    summary.minorticks_off()
    summary.set_xlabel("size of the move", labelpad=2)
    summary.set_ylabel("gap", labelpad=6)
    summary.set_title("Bigger moves, bigger gap, up to a limit", loc="left",
                      pad=8)
    summary.spines[["left", "bottom"]].set_color(GHOST)
    summary.axhline(BETA, color=INK, ls=":", lw=1.4)
    summary.text(0.48, BETA + 0.12, "maximum gap", fontsize=10, color=MUTED,
                 va="bottom")
    summary.text(biggest * 2.4, BETA + 0.12, f"≈ {BETA:.1f} µm", fontsize=10,
                 color=MUTED, ha="right", va="bottom")
    art["curve"], = summary.plot([], [], color=TEAL, lw=2.5, zorder=2)
    art["points"] = [summary.plot([], [], "o", color=INK, ms=9, zorder=3)[0]
                     for _ in moves]

    art["message"] = fig.text(0.075, 0.072, "", fontsize=13, color=INK,
                              va="center")
    footer(fig, "each move is repeated until the loop settles")
    return art


# -- frames ---------------------------------------------------------------
def clear_move(art: dict) -> None:
    art["branch_up"].set_data([], [])
    art["branch_down"].set_data([], [])
    art["live"].set_data([], [])
    art["gap_text"].set_visible(False)
    for key in ("gap_arrow", "arrow_up", "arrow_down"):
        if art[key] is not None:
            art[key].remove()
            art[key] = None


def set_branch_arrow(art: dict, key: str, move: dict, half: int | None,
                     k: int, colour: str) -> None:
    """An arrowhead riding along the branch, pointing the way it is drawn."""
    if art[key] is not None:
        art[key].remove()
        art[key] = None
    if half is None or half < 0:
        return
    start = FPA + half * FPH
    stop = min(k, start + FPH)
    if stop - start < 8:
        return
    mid = (start + stop) // 2
    a, b = mid - 3, mid + 3
    art[key] = art["stage"].annotate(
        "", xy=(move["x"][b], move["lag"][b]),
        xytext=(move["x"][a], move["lag"][a]),
        arrowprops=dict(arrowstyle="-|>", color=colour, lw=2.0,
                        mutation_scale=17), zorder=7)


def fade_current(art: dict, move: dict) -> None:
    """Leave the finished loop behind in light grey, for comparison."""
    cycling = slice(FPA, None)
    line, = art["stage"].plot(move["x"][cycling], move["lag"][cycling],
                              color=GHOST, lw=1.6, zorder=1,
                              solid_capstyle="round")
    art["old_loops"].append(line)


def draw_gap(art: dict, move: dict, up: int, down: int) -> None:
    a, b = half_slice(up), half_slice(down)
    low = float(np.interp(0.5, move["x"][a], move["lag"][a]))
    high = float(np.interp(0.5, move["x"][b][::-1], move["lag"][b][::-1]))
    if art["gap_arrow"] is not None:
        art["gap_arrow"].remove()
    # A thin gap needs small heads, or the heads alone would make it look
    # several times wider than it is.
    ends = art["stage"].transData.transform([(0.5, low), (0.5, high)])
    scale = float(np.clip(abs(ends[1][1] - ends[0][1]) * 0.8, 5, 16))
    art["gap_arrow"] = art["stage"].annotate(
        "", xy=(0.5, low), xytext=(0.5, high),
        arrowprops=dict(arrowstyle="<->", color=INK, lw=1.4,
                        mutation_scale=scale), zorder=6)
    # Beside the arrow when there is room, above its top end when the gap is
    # too thin for a label to sit inside it.
    if abs(ends[1][1] - ends[0][1]) > 40:
        art["gap_text"].set_position((0.54, (low + high) / 2))
        art["gap_text"].set_ha("left")
        art["gap_text"].set_va("center")
    else:
        art["gap_text"].set_position((0.5, max(low, high) + 0.06))
        art["gap_text"].set_ha("center")
        art["gap_text"].set_va("bottom")
    art["gap_text"].set_text("gap")
    art["gap_text"].set_visible(True)


def render(moves: list[dict]) -> list[Image.Image]:
    art = build_figure(moves)
    fig = art["fig"]
    frames: list[Image.Image] = []

    # phase 0: the empty panels and the opening line
    art["message"].set_text(MESSAGES[0])
    frames += [capture(fig)] * INTRO

    for index, move in enumerate(moves):
        art["message"].set_text(MESSAGES[index + 1])
        clear_move(art)
        steps = len(move["u"]) - 1
        time = np.linspace(0, 1, steps + 1)
        art["plan"].set_data(time, move["n"])
        art["size_tag"].set_text(
            f"{LABELS[index]} · {move['s']:g} µm" if index == 0 else
            f"{LABELS[index]} · {move['s']:g} µm  "
            f"(×{move['s'] / moves[0]['s']:g} bigger)")
        for k in range(steps + 1):
            art["trail"].set_data(time[:k + 1], move["n"][:k + 1])
            art["dot"].set_data([time[k]], [move["n"][k]])
            rising = k > 0 and move["u"][k] > move["u"][k - 1]
            art["up"].set_data([time[k]] if rising else [],
                               [move["n"][k] + 0.33] if rising else [])
            art["down"].set_data([] if rising or k == 0 else [time[k]],
                                 [] if rising or k == 0 else [move["n"][k] - 0.33])
            up_x, up_y, down_x, down_y = branches(move, k)
            art["branch_up"].set_data(up_x, up_y)
            art["branch_down"].set_data(down_x, down_y)
            if k >= FPA:
                art["live"].set_data([move["x"][k]], [move["lag"][k]])
            current = min(2 * CYCLES - 1, max(0, (k - FPA - 1) // FPH))
            latest_up = current if current % 2 == 0 else current - 1
            latest_down = current if current % 2 else current - 1
            set_branch_arrow(art, "arrow_up", move, latest_up, k, TEAL)
            set_branch_arrow(art, "arrow_down", move, latest_down, k, ORANGE)
            if k == FPA + 2 * FPH:
                draw_gap(art, move, 0, 1)
            if k == FPA + 4 * FPH:
                draw_gap(art, move, 2, 3)
            frames.append(capture(fig))

        # the move is over: its point drops into the bottom panel
        art["points"][index].set_data([move["s"]], [move["gap"]])
        frames += [capture(fig)] * HOLD
        fade_current(art, move)

    # phase 4: everything faint, the three points joined
    clear_move(art)
    art["trail"].set_data([], [])
    art["dot"].set_data([], [])
    art["up"].set_data([], [])
    art["down"].set_data([], [])
    art["plan"].set_data([], [])
    art["size_tag"].set_text("")
    starts = np.linspace(0, 1, len(moves) + 1)
    for index, move in enumerate(moves):
        span = np.linspace(starts[index], starts[index + 1] - 0.02,
                           len(move["u"]))
        line, = art["motor"].plot(span, move["n"], color=GHOST, lw=2)
        art["faint_plans"].append(line)
        art["faint_tags"].append(art["motor"].text(
            (starts[index] + starts[index + 1]) / 2 - 0.01, -1.36,
            f"{move['s']:g} µm", ha="center", va="center", fontsize=10,
            color=MUTED))
    sizes = np.geomspace(STROKES[0] * 0.5, STROKES[-1] * 2.0, 90)
    art["curve"].set_data(sizes, [width_at_centre(float(v)) for v in sizes])
    art["message"].set_text(MESSAGES[4])
    frames += [capture(fig)] * SUMMARY
    frames += [frames[-1]] * TAIL
    plt.close(fig)
    return frames


def main() -> int:
    moves = build_moves()

    # Both drawn cycles must be the same loop, or the settling run before
    # them did not do its job and the gap would change on screen.
    for move in moves:
        first, second = gap_between(move, 0, 1), gap_between(move, 2, 3)
        assert abs(first - second) < 1e-9, (move["s"], first, second)

    # The three animated gaps must land on the curve drawn in panel (C);
    # the curve is the settled width for a fresh loop of the same size.
    for move in moves:
        settled = width_at_centre(move["s"])
        assert abs(move["gap"] - settled) < 2e-3, (move["s"], move["gap"], settled)
    assert abs(moves[-1]["gap"] - BETA) < 1e-6, "the largest move must saturate"
    assert moves[0]["gap"] < moves[1]["gap"] < moves[2]["gap"], "gap must grow"

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)   # catch missing glyphs
        frames = render(moves)
    save_gif(frames, GIF)

    report(GIF, frames)
    for move, label in zip(moves, LABELS):
        print(f"  {label:12s} {move['s']:6g} µm move -> gap "
              f"{move['gap']:.3f} µm ({move['gap'] / BETA * 100:.0f}% of the "
              f"maximum)")
    print(f"  maximum gap (saturation) {BETA:g} µm, "
          f"illustrative s90 {S90_CFG:g} µm, {len(LIMITS)} elements, "
          f"stiffness {G.sum():.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
