#!/usr/bin/env python3
"""The 2.500 mm/s limit, one idea per page.

Four pages, each making a single point and nothing else:

    1  what the experiment was
    2  the limit, and where the number comes from
    3  what it does - the record runs late, the stage is not slow
    4  what it costs the measurement

Writes plots/temp/05_recorded_velocity_limit_explained.pdf (4 pages) and the
same four pages as separate PNGs.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

_HERE = Path(__file__).resolve().parent
_MOTOR = _HERE.parent.parent.parent / "motor identification" / "scripts"
if str(_MOTOR) not in sys.path:
    sys.path.insert(0, str(_MOTOR))
from staircase_figures import INK, INK2, INK3, SURFACE, GRID, style   # noqa

SRC = Path(os.environ.get("ISO_TEST_CSV", _HERE.parent.parent / "IsoTest.csv"))
OUT = Path(os.environ.get("ISO_TEMP_OUT", _HERE / "plots" / "temp"))

HEADER_ROWS = 23
WRAP = 2 ** 32
NM_PER_COUNT = 1.0          # Resolution AquadB 1000 pm
CLOCK_NS = 400.0            # AquadB clock
CEILING = 1e9 / CLOCK_NS * NM_PER_COUNT / 1e6     # 2.500 mm/s
CEILING_NM = CEILING * 1e3                        # per 1 ms sample
L_MM = 59.338

PROFILES = [("A", 5.0, 40.0), ("B", 5.0, 10.0),
            ("C", 2.0, 40.0), ("D", 1.0, 5.0)]
COLOUR = {"A": "#b3261e", "B": "#e08214", "C": "#2a78d6", "D": "#1f7a4d"}
RED, GREEN = "#b3261e", "#1f7a4d"


def load(path=SRC):
    idx, val = [], []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i < HEADER_ROWS:
                continue
            line = line.strip()
            if not line or line.startswith("EOF"):
                continue
            p = line.split("\t")
            if len(p) == 2:
                idx.append(int(p[0]))
                val.append(int(p[1]))
    val = np.array(val, np.int64)
    dv = np.diff(val)
    dv[dv > WRAP // 2] -= WRAP
    dv[dv < -WRAP // 2] += WRAP
    pos = np.empty_like(val)
    pos[0] = val[0]
    pos[1:] = val[0] + np.cumsum(dv)
    return np.array(idx, np.int64), pos


def transfers(pos, min_mm=50.0):
    d = np.diff(pos)
    sp = np.convolve(np.abs(d).astype(float), np.ones(101) / 101,
                     mode="same") / 1e3
    moving = sp > 0.03
    out, i = [], 0
    while i < len(moving):
        if moving[i]:
            j = i
            while j < len(moving) and moving[j]:
                j += 1
            if abs(pos[j] - pos[i]) / 1e6 > min_mm:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def measure(pos, runs):
    """Four profiles, four runs each; keep the three plain transfers."""
    d = np.diff(pos).astype(float)
    per = len(runs) // len(PROFILES)
    rows = []
    for k, (name, v, acc) in enumerate(PROFILES):
        block = runs[k * per:(k + 1) * per]
        clean = sorted(block, key=lambda ab: ab[1] - ab[0])[:3]
        a, b = clean[0]
        seg = np.abs(d[a + 2000:b - 2000])
        sm = np.convolve(np.abs(d[b - 1500:b + 100]), np.ones(15) / 15,
                         mode="valid")
        pk = float(np.median(sm[:400]))
        hi = np.where(sm >= 0.9 * pk)[0]
        lo = np.where(sm <= 0.10 * pk)[0]
        lo = lo[lo > hi[-1]] if len(hi) else lo
        tail = np.convolve(np.abs(d[b:b + 12000]), np.ones(51) / 51,
                           mode="same")
        n = 0
        for s in tail < 20.0:
            if s:
                n += 1
            elif n > 200:
                break
            else:
                n = 0
        rows.append(dict(
            name=name, v=v, acc=acc, slice=(a, b), block=(block[0][0],
                                                          block[-1][1]),
            dist=np.mean([abs(pos[y] - pos[x]) / 1e6 for x, y in clean]),
            recorded=np.mean([(y - x) / 1000.0 for x, y in clean]),
            physics=L_MM / v + v / acc,
            drain=L_MM / CEILING,
            clip=100.0 * float(np.mean(seg == CEILING_NM)),
            decel=int(lo[0] - hi[-1]) if len(lo) and len(hi) else -1,
            decel_needed=1000.0 * min(v, CEILING) / acc,
            still=n / 1000.0))
    return rows


def head(fig, title, line):
    fig.text(0.055, 0.955, title, size=19, weight="bold", color=INK, va="top")
    fig.text(0.055, 0.895, line, size=11.5, color=INK2, va="top")


# --------------------------------------------------------------- page 1
def page_experiment(pos, rows):
    fig = plt.figure(figsize=(12.6, 7.0))
    head(fig, "1.  What the experiment was",
         "The same 59.338 mm move, driven four different ways. Nothing else "
         "changes between them.")
    ax = fig.add_axes([0.075, 0.300, 0.895, 0.520])
    t = np.arange(len(pos)) / 1000.0
    ax.plot(t, (pos - pos.min()) / 1e6, color=INK3, lw=0.9)
    for r in rows:
        a, b = r["block"]
        ax.axvspan(a / 1000.0, b / 1000.0, color=COLOUR[r["name"]],
                   alpha=0.16, lw=0)
        ax.text((a + b) / 2000.0, 128, r["name"], ha="center", size=15,
                weight="bold", color=COLOUR[r["name"]])
        ax.text((a + b) / 2000.0, 118, f"{r['v']:g} mm/s\n{r['acc']:g} mm/s²",
                ha="center", size=9.2, color=INK2)
    ax.text(340, 128, "switch baseline", ha="center", size=10.5, color=INK2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("stage position (mm)")
    ax.set_xlim(0, t[-1])
    ax.set_ylim(-5, 140)
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")

    rowtxt = [("", "commanded speed", "commanded acceleration", "")]
    for r in rows:
        rowtxt.append((r["name"], f"{r['v']:g} mm/s", f"{r['acc']:g} mm/s²",
                       ""))
    y = 0.215
    fig.text(0.075, y, "profile", size=10, weight="bold", color=INK)
    fig.text(0.175, y, "speed", size=10, weight="bold", color=INK)
    fig.text(0.300, y, "acceleration", size=10, weight="bold", color=INK)
    fig.text(0.460, y, "what it isolates", size=10, weight="bold", color=INK)
    notes = ["the setting v2 ran",
             "same speed, acceleration cut 4x",
             "same acceleration, speed cut 2.5x",
             "both well below the limit"]
    for i, (r, note) in enumerate(zip(rows, notes)):
        yy = y - 0.040 * (i + 1)
        fig.text(0.075, yy, r["name"], size=11, weight="bold",
                 color=COLOUR[r["name"]])
        fig.text(0.175, yy, f"{r['v']:g} mm/s", size=10.5, color=INK)
        fig.text(0.300, yy, f"{r['acc']:g} mm/s²", size=10.5, color=INK)
        fig.text(0.460, yy, note, size=10.5, color=INK2)
    fig.text(0.075, 0.035,
             "Each profile runs the move twice, holding 10 s at each end, "
             "then re-homes. A against B changes only acceleration; "
             "A against C changes only speed.",
             size=10, color=INK2)
    return fig


# --------------------------------------------------------------- page 2
def page_limit(pos, rows):
    fig = plt.figure(figsize=(12.6, 7.0))
    head(fig, "2.  The limit",
         "The recorded speed can never exceed 2.500 mm/s, and that number "
         "comes from two interface settings.")
    ax = fig.add_axes([0.055, 0.300, 0.390, 0.520])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    items = [(0.93, "Resolution AquadB", "1000 pm", "one increment = 1 nm"),
             (0.76, "AquadB clock", "400 ns",
              "at best one increment per 400 ns")]
    for yy, k, v, note in items:
        ax.text(0.0, yy, k, size=12, weight="bold", color=INK, va="center")
        ax.text(0.0, yy - 0.075, v, size=15, color=RED, va="center")
        ax.text(0.30, yy - 0.075, note, size=10.5, color=INK2, va="center")
    ax.plot([0.0, 1.0], [0.60, 0.60], color=GRID, lw=1.2)
    ax.text(0.0, 0.50, "1 ÷ 400 ns  =  2 500 000 increments/s", size=12.5,
            color=INK, va="center")
    ax.text(0.0, 0.40, "× 1 nm each", size=12.5, color=INK, va="center")
    ax.plot([0.0, 1.0], [0.33, 0.33], color=GRID, lw=1.2)
    ax.text(0.0, 0.21, "2.500 mm/s", size=26, weight="bold", color=INK,
            va="center")
    ax.text(0.0, 0.07, "= 2500 counts in every 1 ms sample", size=11.5,
            color=INK2, va="center")

    bx = fig.add_axes([0.560, 0.300, 0.410, 0.520])
    d = np.diff(pos).astype(float)
    steady = np.concatenate([np.abs(d[r["slice"][0] + 2000:
                                     r["slice"][1] - 2000]) for r in rows])
    below = steady[steady < CEILING_NM]
    bx.hist(below, bins=np.arange(0, CEILING_NM + 1, 50), color=INK3, lw=0)
    bx.bar([CEILING_NM], [int(np.sum(steady == CEILING_NM))], width=46,
           align="edge", color=RED, lw=0)
    bx.set_yscale("log")
    bx.set_xlim(0, CEILING_NM + 320)
    bx.set_xlabel("distance moved in one 1 ms sample (nm)")
    bx.set_ylabel("number of samples")
    bx.grid(True, axis="y")
    bx.grid(False, axis="x")
    bx.annotate(f"{int(np.sum(steady == CEILING_NM)):,} samples sit exactly\n"
                f"on 2500 nm, and not one goes past it",
                xy=(CEILING_NM, np.sum(steady == CEILING_NM)),
                xytext=(-30, -70), textcoords="offset points", ha="right",
                size=10.5, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    fig.text(0.055, 0.195, "The measured wall and the arithmetic agree "
                           "exactly.", size=12, weight="bold", color=INK)
    fig.text(0.055, 0.140,
             "This is a limit of the link that carries position from the "
             "interferometer to the counter — not of the stage, the motor, "
             "or the interferometer itself.",
             size=10.5, color=INK2)
    return fig


# --------------------------------------------------------------- page 3
def page_late(pos, rows):
    fig = plt.figure(figsize=(12.6, 7.0))
    head(fig, "3.  What the limit does",
         "Above 2.500 mm/s the record falls behind the stage. The stage is "
         "not slow — the recording is late.")
    r = rows[0]
    a, b = r["slice"]
    d = np.diff(pos).astype(float)
    ax = fig.add_axes([0.075, 0.300, 0.895, 0.470])
    sm = np.convolve(np.abs(d[a:a + 34000]), np.ones(101) / 101,
                     mode="same") / 1e3
    t = np.arange(len(sm)) / 1000.0
    arrive = r["physics"]
    ax.axvspan(arrive, arrive + 10.0, color=GREEN, alpha=0.12, lw=0)
    ax.plot(t, sm, color=COLOUR["A"], lw=1.8)
    ax.axhline(CEILING, color=INK, lw=1.2, ls=(0, (5, 3)))
    ax.axvline(arrive, color=INK, lw=1.4)
    ax.axvline(r["recorded"], color=RED, lw=1.4, ls=(0, (4, 2)))
    ax.annotate("", xy=(arrive, 2.80), xytext=(r["recorded"], 2.80),
                arrowprops=dict(arrowstyle="<->", color=RED, lw=1.6))
    ax.text((arrive + r["recorded"]) / 2, 2.86,
            f"{r['recorded'] - arrive:.1f} s of backlog still arriving",
            ha="center", size=12, color=RED, weight="bold")
    ax.text(arrive - 0.5, 1.30, "the stage\nstops here", ha="right", size=11,
            color=INK)
    ax.text(r["recorded"] + 0.6, 1.85, "the record\nstops here", ha="left",
            size=11, color=RED)
    ax.text(arrive + 5.0, 0.16, "the intended 10 s hold", ha="center",
            size=11, color=GREEN)
    ax.set_xlabel("time from the start of the move (s)")
    ax.set_ylabel("recorded speed (mm/s)")
    ax.set_xlim(0, 32)
    ax.set_ylim(0, 3.15)
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")
    ax.set_title("profile A — commanded 5 mm/s", size=12, loc="left",
                 color=INK2)

    fig.text(0.055, 0.205, "How we know the stage was not simply slow",
             size=12, weight="bold", color=INK)
    fig.text(0.055, 0.150,
             f"Profile B is commanded 5 mm/s at 10 mm/s², which physically "
             f"needs {rows[1]['decel_needed']:.0f} ms to come to a stop. Its "
             f"recorded move ends in {rows[1]['decel']} ms — far too fast to "
             f"be real motion.",
             size=10.8, color=INK2)
    fig.text(0.055, 0.100,
             f"Profiles C and D stay under the limit and show exactly the "
             f"stopping time their acceleration implies "
             f"({rows[2]['decel']} ms against {rows[2]['decel_needed']:.0f}, "
             f"{rows[3]['decel']} ms against {rows[3]['decel_needed']:.0f}). "
             f"A sudden stop is a queue running dry, not a stage braking.",
             size=10.8, color=INK2)
    return fig


# --------------------------------------------------------------- page 4
def page_effect(pos, rows):
    fig = plt.figure(figsize=(12.6, 7.0))
    head(fig, "4.  What it costs the measurement",
         "Positions and distances are unharmed. Timing, and the stationary "
         "window an ISO stop needs, are not.")
    names = [r["name"] for r in rows]
    xs = np.arange(len(rows))

    ax = fig.add_axes([0.075, 0.330, 0.370, 0.470])
    dists = [r["dist"] for r in rows]
    ax.bar(xs, [dd - 59.330 for dd in dists], bottom=59.330, width=0.6,
           color=[COLOUR[n] for n in names], lw=0)
    ax.set_ylim(59.330, 59.345)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r['name']}\n{r['clip']:.0f}% at the limit"
                        for r in rows], size=9.6)
    ax.set_ylabel("distance actually moved (mm)")
    ax.set_title("unharmed:  every move lands in the same place",
                 size=12, loc="left", color=GREEN)
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")
    ax.text(0.5, 0.06,
            f"spread {1000*(max(dists)-min(dists)):.1f} µm across profiles "
            f"that\nsaturate anywhere from {min(r['clip'] for r in rows):.0f}%"
            f" to {max(r['clip'] for r in rows):.0f}% of the time",
            transform=ax.transAxes, ha="center", size=10, color=INK2)

    bx = fig.add_axes([0.600, 0.330, 0.370, 0.470])
    bx.bar(xs, [r["still"] for r in rows], width=0.6,
           color=[COLOUR[n] for n in names], lw=0)
    bx.axhline(10.0, color=INK, lw=1.3, ls=(0, (5, 3)))
    bx.text(3.45, 10.25, "the 10 s hold that was asked for", ha="right",
            size=10, color=INK)
    for i, r in enumerate(rows):
        bx.text(i, r["still"] + 0.3, f"{r['still']:.1f} s", ha="center",
                size=11, color=INK)
    bx.set_xticks(xs)
    bx.set_xticklabels([f"{r['name']}\n{r['v']:g} mm/s" for r in rows],
                       size=9.6)
    bx.set_ylabel("hold time recorded as stationary (s)")
    bx.set_ylim(0, 12)
    bx.set_title("harmed:  the backlog eats the stop window",
                 size=12, loc="left", color=RED)
    bx.grid(True, axis="y")
    bx.grid(False, axis="x")

    fig.text(0.055, 0.225, "What this means in practice", size=12,
             weight="bold", color=INK)
    for i, line in enumerate([
            "Keep the commanded speed at or below 2.5 mm/s and the recording "
            "follows the stage honestly.",
            "Above it, positions are still correct, but the trace is "
            "stretched in time and a 10 s hold can shrink to under 2 s of "
            "usable stationary data.",
            "ISO 230-2 reads only settled positions, so a run made under the "
            "limit — like v3 at 2 mm/s — is unaffected."]):
        fig.text(0.055, 0.170 - 0.047 * i, "•  " + line, size=10.8,
                 color=INK2)
    return fig


def main() -> int:
    if not SRC.is_file():
        print(f"no such file: {SRC}")
        return 2
    style()
    idx, pos = load()
    rows = measure(pos, transfers(pos))

    print(f"{'':2} {'cmd':>11} {'dist mm':>8} {'physics':>8} {'drain':>7} "
          f"{'recorded':>9} {'clip%':>6} {'decel':>7} {'needed':>7} "
          f"{'still s':>8}")
    for r in rows:
        print(f"{r['name']:2} {r['v']:5.1f}/{r['acc']:<5.0f} {r['dist']:8.3f} "
              f"{r['physics']:8.2f} {r['drain']:7.2f} {r['recorded']:9.2f} "
              f"{r['clip']:6.1f} {r['decel']:7d} {r['decel_needed']:7.0f} "
              f"{r['still']:8.1f}")

    OUT.mkdir(parents=True, exist_ok=True)
    pages = [("05a_experiment", page_experiment),
             ("05b_the_limit", page_limit),
             ("05c_record_runs_late", page_late),
             ("05d_effect_on_measurement", page_effect)]
    pdf_path = OUT / "05_recorded_velocity_limit_explained.pdf"
    with PdfPages(pdf_path) as pdf:
        for stem, fn in pages:
            fig = fn(pos, rows)
            fig.savefig(OUT / f"{stem}.png")
            pdf.savefig(fig)
            plt.close(fig)
            print(f"  {stem}.png")
    print(f"\nwrote {pdf_path} (4 pages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
