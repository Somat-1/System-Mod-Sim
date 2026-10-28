#!/usr/bin/env python3
"""The E1 / ISO 230-2 reversal disagreement, and where it came from.

It was not the stage. The last few tens of milliseconds of every ISO 230-2
record window land on the start of the next move, because the board's record
end maps about 37 ms early onto the scope timeline. The next move always runs
in the direction of the pass in progress, so the contamination pushes the
forward means up and the reverse means down - straight into B, which is their
difference. Dropping the moving tail brings the published -0.365 um to
-1.600 um, against E1's -1.680 um.

The direction average is almost untouched (the two biases cancel), so the gain
fit built on it survives: -24.17 ppm becomes -24.10 ppm. B, R, E and A do not.

Needs the raw ISO capture at stage identification/ISO230-2_v3.csv.

Writes rendered/reversal_reconciliation.{png,pdf}.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import e1_run
from e1_run import WINDOWS

_HERE = Path(__file__).resolve().parent
_STAGE = _HERE.parents[2]
for extra in (_HERE.parents[3] / "motor identification" / "scripts",
              _STAGE / "iso230-2" / "analysis"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))
from staircase_figures import INK, INK2, style                    # noqa: E402
from iso230_2_evaluate import load_scope, load_stops, fit_clock   # noqa: E402

ISO_SCOPE = _STAGE / "ISO230-2_v3.csv"
ISO_LOG = _STAGE / "iso230-2" / "data" / "terminal_20261021_183115.log"
OUT = Path(os.environ.get("E1_OUT", _HERE / "rendered"))
COLOURS = {"W1": "#4c78a8", "W2": "#598a76", "W3": "#9d6770"}
ISO_COLOUR = "#161616"
BAD, GOOD = "#99483e", "#356449"
MINUS = "−"
MOVING_NM = 60.0          # 60 nm in 1 ms: the stage is no longer holding


def num(value, places=2):
    return f"{value:.{places}f}".replace("-", MINUS)


def thesis_style():
    style()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "figure.facecolor": "#ffffff",
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": "#ffffff",
    })


def iso_windows():
    """Each ISO record window as the evaluation maps it, 500 samples long."""
    pos = load_scope(ISO_SCOPE)
    stops, _ = load_stops(ISO_LOG)
    rate, off, _ = fit_clock(stops, pos)
    iso = [s for s in stops if s["state"] == "ISO"]
    seg, meta = [], []
    for s in iso:
        a = int(round(s["t0"] * rate + off))
        b = int(round(s["t1"] * rate + off))
        if b - a != 500 or b > pos.size:
            continue
        seg.append(pos[a:b])
        meta.append((s["i"], s["dir"], s["cmd_mm"]))
    return np.array(seg), np.array(meta)


def reversal(seg, meta, trim):
    """B_i with the last `trim` samples of each window dropped."""
    meas = seg[:, :500 - trim].mean(axis=1)
    idx, dirs, cmd = meta[:, 0], meta[:, 1], meta[:, 2]
    dev = (meas - np.mean(meas - cmd) - cmd) * 1000.0
    return np.array([dev[(idx == i) & (dirs == 1)].mean()
                     - dev[(idx == i) & (dirs == -1)].mean()
                     for i in np.unique(idx)])


def e1_windows(run):
    st = run.stops
    i0 = run.index_at(st.record_start_s)
    return np.array([run.pos[i:i + 500] for i in i0]), st.direction


def main():
    thesis_style()
    seg, meta = iso_windows()
    run = e1_run.load()
    e1seg, e1dir = e1_windows(run)

    # windows referenced to their own settled part, signed along the travel
    iso_curve = ((seg - seg[:, :400].mean(axis=1, keepdims=True)) * 1000.0
                 * meta[:, 1][:, None])
    e1_curve = ((e1seg - e1seg[:, :400].mean(axis=1, keepdims=True)) * 1000.0
                * e1dir[:, None])

    trims = np.arange(0, 320, 5)
    Bt = np.array([reversal(seg, meta, t).mean() for t in trims])
    B0 = reversal(seg, meta, 0)
    Bc = reversal(seg, meta, 100)
    st = run.stops
    resid = (run.hold() - st.target_mm) * 1000.0
    e1B = []
    for name, _, _ in WINDOWS:
        m = st.window == name
        up, dn = m & (st.direction == 1), m & (st.direction == -1)
        ou, od = np.argsort(st.target_mm[up]), np.argsort(st.target_mm[dn])
        e1B.append(resid[up][ou] - resid[dn][od])
    e1B = np.concatenate(e1B)

    fig = plt.figure(figsize=(11.0, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.95), hspace=0.45,
                          wspace=0.26)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # (a) the contaminated tail
    t = np.arange(500)
    for curve, colour, label in ((iso_curve, ISO_COLOUR,
                                  "ISO 230-2 v3, 240 windows"),
                                 (e1_curve, COLOURS["W2"],
                                  "E1, 654 windows")):
        lo, hi = np.percentile(curve, [5, 95], axis=0)
        ax_a.fill_between(t, lo, hi, color=colour, alpha=0.18, lw=0)
        ax_a.plot(t, curve.mean(axis=0), lw=1.6, color=colour, label=label)
    first = np.array([np.flatnonzero(np.abs(np.diff(s)) * 1e6 > MOVING_NM)
                      for s in seg], dtype=object)
    onset = np.array([f[0] if len(f) else 500 for f in first])
    ax_a.axvline(np.median(onset), color=BAD, lw=1.0, ls="--")
    ax_a.annotate(f"the next move starts here,\n{int(500 - np.median(onset))} ms "
                  "before the window ends",
                  (np.median(onset), 2.4), xytext=(-12, 0),
                  textcoords="offset points", ha="right", fontsize=9,
                  color=BAD)
    ax_a.set_xlabel("sample within the 500 ms record window [ms]")
    ax_a.set_ylabel("position, signed along\nthe travel [µm]")
    ax_a.set_ylim(-0.6, 4.0)
    ax_a.text(498, 3.85, "(continues to +23 µm)", ha="right", va="top",
              fontsize=8.5, color=ISO_COLOUR)
    ax_a.set_xlim(0, 500)
    ax_a.legend(loc="upper left", fontsize=9, framealpha=1.0)
    ax_a.set_title("(a) every record window, referenced to its own settled "
                   "part — median and 5–95 %", fontsize=11, loc="left")

    # (b) B against how much of the tail is dropped
    ax_b.plot(trims, Bt, "-", lw=1.8, color=ISO_COLOUR, zorder=4)
    ax_b.plot([0], [Bt[0]], "o", ms=7, mfc="#ffffff", mec=BAD, mew=1.6,
              zorder=5)
    ax_b.annotate(f"published\n{num(Bt[0])} µm", (0, Bt[0]),
                  xytext=(14, -4), textcoords="offset points", fontsize=9,
                  color=BAD, va="top")
    ax_b.axhline(e1B.mean(), color=COLOURS["W2"], lw=1.2, ls="--", zorder=3)
    ax_b.set_ylim(e1B.mean() - 0.12, Bt[0] + 0.18)
    ax_b.annotate(f"E1, {num(e1B.mean())} µm", (310, e1B.mean()),
                  xytext=(0, 5), textcoords="offset points", ha="right",
                  fontsize=9, color=COLOURS["W2"])
    ax_b.annotate(f"plateau, {num(Bt[-1])} µm", (200, Bt[-1]),
                  xytext=(0, 12), textcoords="offset points", ha="center",
                  fontsize=9, color=GOOD)
    ax_b.set_xlabel("samples dropped from the end of each window [ms]")
    ax_b.set_ylabel("mean reversal $\\bar{B}$ [µm]")
    ax_b.set_title("(b) drop the moving tail and the mean settles at once",
                   fontsize=11, loc="left")

    # (c) before and after, against E1
    groups = [("ISO\npublished", B0, BAD),
              ("ISO\ntail dropped", Bc, GOOD),
              ("E1", e1B, COLOURS["W2"])]
    rng = np.random.default_rng(3)
    for k, (label, B, colour) in enumerate(groups):
        ax_c.plot(k + rng.uniform(-0.16, 0.16, B.size), B, ".", ms=3.0,
                  color=colour, alpha=0.6, zorder=2)
        ax_c.plot([k - 0.3, k + 0.3], [B.mean()] * 2, "-", lw=2.2,
                  color=colour, zorder=4)
        ax_c.annotate(num(B.mean()), (k, B.mean()), xytext=(22, -3),
                      textcoords="offset points", fontsize=9, color=colour)
    ax_c.axhline(0.0, color=INK2, lw=0.8)
    ax_c.set_xticks(range(len(groups)))
    ax_c.set_xticklabels([g[0] for g in groups], fontsize=9)
    ax_c.set_xlim(-0.6, len(groups) - 0.35)
    ax_c.set_ylabel("reversal $B_i$ [µm]")
    ax_c.set_title("(c) the two runs agree to "
                   f"{num(abs(Bc.mean() - e1B.mean()))} µm", fontsize=11,
                   loc="left")

    for ax in (ax_a, ax_b, ax_c):
        ax.set_yticklabels([num(v, 1) for v in ax.get_yticks()])

    fig.text(0.012, 0.013,
             "The direction average is barely touched by this - the forward "
             "and reverse biases cancel - so the gain fit built on it holds: "
             f"{MINUS}24.17 ppm becomes {MINUS}24.10 ppm. "
             "B, R, E and A do not.", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.028, 1, 1))
    OUT.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"reversal_reconciliation.{ext}", dpi=200)
    print("wrote", OUT / "reversal_reconciliation.png")
    print(f"published {B0.mean():+.3f}  tail dropped {Bc.mean():+.3f}  "
          f"E1 {e1B.mean():+.3f} um")
    print(f"median onset of the next move: sample {int(np.median(onset))} of 500")


if __name__ == "__main__":
    main()
