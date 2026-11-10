#!/usr/bin/env python3
"""The friction identification roadmap: what gets measured, in what order.

Written to be read by someone who has never seen the rig. Three rows -
the lab sequence, the one test that runs alongside it, and the desk work -
read left to right, one plain sentence per step, and no symbol that is not
spelled out on the page.

An earlier draft followed ROADMAP_SPEC.md to the letter: five lines in every
box, three decision diamonds with their questions, a feed-forward arc over
the whole lane, a bypass branch and a route-B footnote. All of it was true
and none of it was readable - a roadmap that needs the thesis open beside it
is not doing its job. What stayed is the structure; what went is everything
that only made sense to someone already inside the project. The detail now
lives in README.md and in the caption, where a reader can take it at their
own pace. ROADMAP_SPEC.md records what was dropped.

It is the companion to `Stage & motor/paper_figures/fig_identification.pdf`
and is built the same way: millimetres at final size, orthogonal arrows, one
small filled head, a serif face, and a check after every run that no run of
text is wider than the box it sits in.

    python render_friction_roadmap.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered"

MM = 1 / 25.4
PT = 25.4 / 72.0                        # one point, in millimetres
W = 180.0

# ------------------------------------------------------------------ style
# The spec's colour table. Its three pale fills sit within five grey levels
# of each other, which would normally be a fault - but the three rows
# are separated by vertical position and the flow types by solid against
# dashed, which is what the spec's own black-and-white clause asks for. The
# fills are decoration, not an encoding. The one fill that has to survive a
# grey print is F5's, and at luminance 71 against 239 it does.
LOOK = {
    "lab": {"face": "#E8F0F5", "edge": "#1F4E79", "title": "#1F4E79",
            "body": "#222222"},
    "goal": {"face": "#1F4E79", "edge": "#1F4E79", "title": "#FFFFFF",
             "body": "#FFFFFF"},
    "parallel": {"face": "#F5ECE0", "edge": "#A0692A", "title": "#A0692A",
                 "body": "#222222"},
    "desk": {"face": "#F1F2F4", "edge": "#7A8088", "title": "#222222",
             "body": "#222222"},
}
BODY = "#222222"
TAG = "#6B7075"
ARROW = "#404040"
PT_TITLE, PT_BODY, PT_EDGE = 7.6, 7.0, 6.4

LW = 0.6
LW_CENTRAL = 1.3                        # F1, the central experiment
HEAD = 1.4                              # arrow head length, mm
CLEAR = 0.7                             # mm of air at a head
INSET = 1.5                             # mm from a box edge to its text
DASH = (0, (2.0 / LW, 1.2 / LW))        # [2.0, 1.2] pt at 0.6 pt width

# ----------------------------------------------------------------- layout
# One column grid for the lab row: six positions, equal widths, with the
# gap before F5 widened so the desk work has somewhere to sit underneath.
X0, BOX_W, GAP, WIDE_GAP = 4.5, 23.0, 5.6, 10.0
LAB_X = [X0 + i * (BOX_W + GAP) for i in range(5)]
LAB_X.append(LAB_X[-1] + BOX_W + WIDE_GAP)
RIGHT = LAB_X[-1] + BOX_W               # 174.9 mm

PAR_X, PAR_W = 63.5, 48.0               # F4, spanning under F1 to F3
F0A_X, F0A_W = 4.5, 38.0
ID_X, ID_W = 101.5, 36.0
COMP_X, COMP_W = 143.5, 32.0

DATA_IN = 119.5                         # where the measurements enter, in ID
FORCE_DOWN = 87.5                       # F4 hands off down its own column
UP_TO_F0B = 38.5                        # the F0a riser, inside both boxes

LAB_Y = 10.5

# Row geometry, frozen by `set_layout`.
LAB_H = PAR_Y = PAR_H = DESK_Y = DESK_H = BUS_Y = 0.0
LAB_MID = PAR_MID = DESK_MID = 0.0
HEIGHT = LEGEND_Y = 0.0

OVERFLOW: list[tuple] = []


def set_layout() -> None:
    """Freeze the row geometry."""
    global LAB_H, PAR_Y, PAR_H, DESK_Y, DESK_H, BUS_Y
    global LAB_MID, PAR_MID, DESK_MID, HEIGHT, LEGEND_Y
    LAB_H = 12.5
    BUS_Y = 27.5
    # F4 carries one line, not two, so its row is shallower; it stays
    # centred on where a full-height row would have put it.
    PAR_Y, PAR_H = 32.8, 10.0
    DESK_Y, DESK_H = 51.5, 12.5
    LAB_MID = LAB_Y + LAB_H / 2
    PAR_MID = PAR_Y + PAR_H / 2
    DESK_MID = DESK_Y + DESK_H / 2
    HEIGHT = DESK_Y + DESK_H + 11.5
    LEGEND_Y = HEIGHT - 5.0


def style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": PT_BODY,
        "axes.linewidth": LW,
        "lines.linewidth": LW,
    })


# ------------------------------------------------------------- primitives

def watch(artist, limit_mm, where):
    OVERFLOW.append((artist, limit_mm, where))
    return artist


def box(ax, x, y, w, h, title, body=(), look="lab", dashed=False,
        central=False):
    """One step: what it is called, and in one sentence what it asks."""
    c = LOOK[look]
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
        linewidth=LW_CENTRAL if central else LW,
        edgecolor=c["edge"], facecolor=c["face"],
        linestyle=DASH if dashed else "solid", zorder=2))
    limit = w - 2 * INSET
    watch(ax.text(x + INSET, y + 1.5, title, fontsize=PT_TITLE, weight="bold",
                  color=c["title"], va="top", ha="left", zorder=3),
          limit, f"{title} / title")
    for i, text in enumerate(body):
        watch(ax.text(x + INSET, y + 5.3 + i * 2.9, text, fontsize=PT_BODY,
                      color=c["body"], va="top", ha="left", zorder=3),
              limit, f"{title} / line {i + 1}")


def run(ax, points, colour=ARROW, dashed=False, zorder=4, lw=LW):
    """A plain orthogonal run of segments, no head."""
    ls = DASH if dashed else "solid"
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        ax.add_line(Line2D([x0, x1], [y0, y1], color=colour, lw=lw,
                           linestyle=ls, zorder=zorder,
                           solid_capstyle="butt"))


def arrow(ax, points, colour=ARROW, dashed=False, zorder=4, head=HEAD):
    """An orthogonal run with one small filled head on the last segment."""
    if len(points) > 2:
        run(ax, points[:-1], colour, dashed, zorder)
    (x0, y0), (x1, y1) = points[-2], points[-1]
    dx, dy = x1 - x0, y1 - y0
    length = (dx * dx + dy * dy) ** 0.5 or 1.0
    x1 -= CLEAR * dx / length
    y1 -= CLEAR * dy / length
    # mutation_scale is pinned to 1 so the head is the size asked for;
    # left alone, annotate scales head_length by the font size.
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=zorder,
                arrowprops=dict(
                    arrowstyle=f"-|>,head_length={head / PT:.3f},"
                               f"head_width={0.45 * head / PT:.3f}",
                    mutation_scale=1.0, linewidth=LW, color=colour,
                    linestyle=DASH if dashed else "solid",
                    shrinkA=0, shrinkB=0))


def width_mm(ax, artist) -> float:
    renderer = ax.figure.canvas.get_renderer()
    bbox = artist.get_window_extent(renderer=renderer)
    return bbox.transformed(ax.transData.inverted()).width


# ---------------------------------------------------------------- content
# One sentence per box, in words that need no glossary. What each step
# actually produces - pre-sliding elements, the Stribeck curve, the GMS
# parameters - is in README.md and in the caption.

LAB = [
    dict(title="Setup", look="lab",
         body=["get the rig ready:", "fine steps, logging"]),
    dict(title="F0b · scoping", look="lab",
         body=["how big is friction", "here at all?"]),
    dict(title="F1 · reversal", look="lab", central=True,
         body=["what changes when", "the stage reverses?"]),
    dict(title="F3 · velocity", look="lab",
         body=["how does it change", "with speed?"]),
    dict(title="F2 · dwell", look="lab", dashed=True,
         body=["does standing still", "matter?"]),
    dict(title="F5 · validation", look="goal",
         body=["does the correction", "actually work?"]),
]

LEGEND = [("line", False, "order of work"),
          ("line", True, "results feed in"),
          ("box", True, "run if earlier results call for it")]



# ----------------------------------------------------------------- figure

def draw() -> plt.Figure:
    OVERFLOW.clear()
    set_layout()
    fig = plt.figure(figsize=(W * MM, HEIGHT * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(HEIGHT, 0)              # millimetres, y downwards
    ax.axis("off")

    # --- the lab row ---------------------------------------------------
    for x, spec in zip(LAB_X, LAB):
        box(ax, x, LAB_Y, BOX_W, LAB_H, **spec)
    for i in range(4):
        arrow(ax, [(LAB_X[i] + BOX_W, LAB_MID), (LAB_X[i + 1], LAB_MID)])

    # Two notes above the row, where nothing else sits. F1's heavier
    # border and F2's dashed one both mean something, and a border is not
    # a thing a first-time reader can be expected to decode on its own.
    for i, note in ((2, "the main experiment"),
                    (4, "run if F0b finds an effect")):
        ax.text(LAB_X[i] + BOX_W / 2, LAB_Y - 1.4, note, fontsize=PT_EDGE,
                style="italic", color=TAG, va="bottom", ha="center")

    # --- the test that runs alongside ----------------------------------
    box(ax, PAR_X, PAR_Y, PAR_W, PAR_H, "F4 · guide friction",
        body=["how much comes from the guides alone?"],
        look="parallel")

    # --- the desk work -------------------------------------------------
    box(ax, F0A_X, DESK_Y, F0A_W, DESK_H, "F0a · desk review",
        body=["what the measurements we", "already have can tell us"],
        look="desk")
    box(ax, ID_X, DESK_Y, ID_W, DESK_H, "Identification",
        body=["turn the measurements", "into a friction model"],
        look="desk")
    box(ax, COMP_X, DESK_Y, COMP_W, DESK_H, "Compensation",
        body=["use the model to", "correct the command"], look="desk")

    # --- what feeds what -----------------------------------------------
    drops = [x + BOX_W / 2 for x in LAB_X[2:5]]
    for x in drops:
        run(ax, [(x, LAB_Y + LAB_H), (x, BUS_Y)], dashed=True)
    run(ax, [(drops[0], BUS_Y), (drops[-1], BUS_Y)], dashed=True)
    arrow(ax, [(DATA_IN, BUS_Y), (DATA_IN, DESK_Y)], dashed=True)
    ax.text(drops[0] + 2.2, BUS_Y - 1.0, "measurements", fontsize=PT_EDGE,
            color=BODY, va="bottom", ha="left")

    arrow(ax, [(UP_TO_F0B, DESK_Y), (UP_TO_F0B, LAB_Y + LAB_H)])
    ax.text(UP_TO_F0B - 1.8, PAR_MID, "what we already know",
            fontsize=PT_EDGE, color=BODY, va="center", ha="right")

    # F4 hands off down its own column and into the side of Identification,
    # not into its top: the top already carries the measurements.
    ochre = LOOK["parallel"]["edge"]
    arrow(ax, [(FORCE_DOWN, PAR_Y + PAR_H), (FORCE_DOWN, DESK_MID),
               (ID_X, DESK_MID)], ochre, dashed=True)
    ax.text(FORCE_DOWN + 2.2, DESK_Y - 1.6, "the guides' share",
            fontsize=PT_EDGE, color=ochre, va="bottom", ha="left")

    arrow(ax, [(ID_X + ID_W, DESK_MID), (COMP_X, DESK_MID)])

    up = LAB_X[-1] + BOX_W / 2
    arrow(ax, [(up, DESK_Y), (up, LAB_Y + LAB_H)])
    ax.text(up - 1.8, PAR_MID, "corrected commands", fontsize=PT_EDGE,
            color=BODY, va="center", ha="right")

    # --- legend --------------------------------------------------------
    fig.canvas.draw()
    x = X0
    for kind, dashed, label in LEGEND:
        if kind == "line":
            run(ax, [(x, LEGEND_Y), (x + 4.0, LEGEND_Y)], ARROW,
                dashed=dashed, zorder=3)
        else:
            ax.add_patch(Rectangle((x, LEGEND_Y - 1.1), 4.0, 2.2,
                                   facecolor=LOOK["lab"]["face"],
                                   edgecolor=LOOK["lab"]["edge"], lw=LW,
                                   linestyle=DASH))
        text = ax.text(x + 5.2, LEGEND_Y, label, fontsize=PT_EDGE,
                       color=BODY, va="center", ha="left")
        x += 5.2 + width_mm(ax, text) + 7.0
    if x - 7.0 > RIGHT:
        print(f"  legend runs to {x - 7.0:.1f} mm, past {RIGHT:.1f} mm")
    return fig


def luminance(hexstr: str) -> float:
    r, g, b = (int(hexstr[i:i + 2], 16) for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def grayscale_report() -> None:
    """What the fills are worth once the colour is gone."""
    vals = {k: luminance(v["face"]) for k, v in LOOK.items()}
    pale = {k: v for k, v in vals.items() if k != "goal"}
    print("  fill luminance: " + ", ".join(f"{k} {v:.0f}"
                                           for k, v in vals.items()))
    print(f"  the three pale fills span "
          f"{max(pale.values()) - min(pale.values()):.1f} grey levels - "
          "decoration only; row position and dash style carry meaning")
    gap = min(pale.values()) - vals["goal"]
    print(f"  F5 against the pale fills: {gap:.0f} levels - "
          f"{'ok' if gap >= 60 else 'TOO CLOSE'}")


def report(ax) -> None:
    """Any run of text wider than the box it sits in, at final size."""
    bad = [(where, width_mm(ax, art), limit)
           for art, limit, where in OVERFLOW
           if width_mm(ax, art) > limit + 0.05]
    if not bad:
        print("  text fits: no block overflows")
        return
    print(f"  {len(bad)} overflowing text runs:")
    for where, width, limit in bad:
        print(f"    {where:<34} {width:5.1f} mm in {limit:5.1f} mm")


def main() -> int:
    style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = draw()
    fig.canvas.draw()
    report(fig.axes[0])
    grayscale_report()
    name = "fig_friction_roadmap.png"
    fig.savefig(OUT / name, dpi=600)
    print(f"  {name}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
