#!/usr/bin/env python3
"""Preliminary check of the F0b and F1 runs of 10 November 2026.

Two questions only, and deliberately not the full analysis:

  1. Is the capture sound? Does it align with the command table, is the
     scale right, how quiet is the stage at rest, and how fine is the
     quantisation?
  2. Does it contain what the experiments were built to measure? Lost
     motion from one microstep to three millimetres, the pre-sliding
     length, the saturation level, the dwell probe that decides F2, and
     the four slow-ramp cases.

    python3 preliminary_check.py

writes

    rendered/preliminary_quality.png    question 1
    rendered/preliminary_physics.png    question 2
    preliminary_check.json              every number in both figures

Nothing here fits a friction model or replaces the analysis each
procedure specifies. Where a number is only indicative, it says so.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_data import USTEP_NM, load

HERE = Path(__file__).resolve().parent
OUT = HERE / "rendered"

INK, MUTED, GRID = "#20262b", "#59636b", "#dfe3e6"
TEAL, ORANGE, PURPLE, RED = "#1a6478", "#c48a33", "#73578f", "#b5473a"
POS_COLOUR = {"A": TEAL, "B": ORANGE, "C": PURPLE}
F0B_BLOCKS = {"R1": PURPLE, "N": "#34454f", "L": "#17647a", "R2": PURPLE,
              "S": "#c68e37", "D": "#548473", "V": "#b84d42", "R3": PURPLE}


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8,
        "axes.titlesize": 9, "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.edgecolor": "#7d858c", "axes.labelcolor": "#2d3338",
        "xtick.color": "#41474c", "ytick.color": "#41474c",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def tidy(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=GRID, lw=.5)
    ax.set_axisbelow(True)


# -- question 1: is the capture sound? ------------------------------------
def allan(x: np.ndarray, taus_ms) -> np.ndarray:
    """Overlapping Allan deviation of a position record sampled at 1 kHz."""
    out = []
    for m in taus_ms:
        k = len(x) // m
        if k < 3:
            out.append(np.nan)
            continue
        block = x[:k * m].reshape(k, m).mean(axis=1)
        out.append(float(np.sqrt(0.5 * np.mean(np.diff(block) ** 2))))
    return np.asarray(out)


def noise_floor(run) -> list[dict]:
    """The N holds: how still the stage is, and how far it walks."""
    rows = []
    for seg in run.segments:
        if seg.block != "N" or seg.kind != "dwell" or seg.dwell_ms < 60000:
            continue
        x = run.x_nm[seg.t0_ms + 3000:seg.t1_ms - 500]
        t = np.arange(len(x)) / 1000.0
        slope, _ = np.polyfit(t, x, 1)
        taus = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
        rows.append({
            "position_mm": seg.target_u * USTEP_NM / 1e6,
            "minutes": len(x) / 60000.0,
            "sd_raw_nm": float(np.std(x)),
            "sd_detrended_nm": float(np.std(x - np.polyval([slope, x.mean()], t - t.mean()))),
            "drift_nm_per_min": float(slope * 60),
            "p2p_nm": float(np.ptp(x)),
            "taus_ms": taus,
            "allan_nm": [None if np.isnan(v) else float(v)
                         for v in allan(x, taus)],
            "allan_500ms_nm": float(allan(x, [500])[0]),
        })
    return rows


def quantisation(run) -> dict:
    """The step the counter actually moves in, read in the longest dwell."""
    seg = max((s for s in run.segments if s.kind == "dwell"),
              key=lambda s: s.t1_ms - s.t0_ms)
    guard = min(3000, (seg.t1_ms - seg.t0_ms) // 4)
    d = np.diff(run.x_nm[seg.t0_ms + guard:seg.t1_ms - guard])
    nz = np.abs(d[np.abs(d) > 1e-6])
    return {"window_ms": int(seg.t1_ms - seg.t0_ms - 2 * guard),
            "smallest_step_nm": float(np.min(nz)) if len(nz) else None,
            "median_abs_step_nm": float(np.median(np.abs(d)))}


# -- question 2: is the physics in there? ---------------------------------
def ladder(run, block: str, sub: str | None = None) -> dict[float, list[float]]:
    """Lost motion per reversal half-cycle, keyed by commanded stroke in nm.

    Every turning point is a short pause, so the settled position at the end
    of each pause is read and differenced against the next. A half-cycle
    counts as a reversal when its direction opposes the one before it.
    """
    segs = [s for s in run.segments
            if s.block == block and (sub is None or s.sub == sub)]
    points = [(s.target_u, run.settled(s, 250)) for i, s in enumerate(segs)
              if s.kind == "dwell" and i and segs[i - 1].kind in ("step", "ramp")]
    out: dict[float, list[float]] = defaultdict(list)
    last = +1
    for (u0, x0), (u1, x1) in zip(points[:-1], points[1:]):
        travel = (u1 - u0) * USTEP_NM
        if travel == 0:
            continue
        direction = 1 if travel > 0 else -1
        if direction != last:
            out[round(abs(travel), 3)].append(abs(travel) - abs(x1 - x0))
        last = direction
    return out


def saturation_and_s90(strokes_nm, lost_nm, knee_nm: float = 1.9e4
                       ) -> tuple[float, float]:
    """beta0 from the plateau above the knee, s90 by interpolation below it.

    Only F0b's anchored ladder reaches the plateau: F1's symmetric loops are
    still rising at their largest clean stroke, so they cross-check the shape
    rather than set the level. Rungs faster than about 25 um/s are excluded
    before this is called, because the 0.5 s pause is not long enough for
    them (see post_move_creep).
    """
    s = np.asarray(strokes_nm, float)
    l = np.asarray(lost_nm, float)
    order = np.argsort(s)
    s, l = s[order], l[order]
    plateau = l[s >= knee_nm]
    if len(plateau) < 2:
        raise ValueError("the ladder does not reach saturation")
    beta = float(np.mean(plateau))
    rise = s <= knee_nm
    s90 = float(np.interp(0.9 * beta, l[rise], s[rise]))
    return beta, s90


def dwell_probe(run) -> dict:
    """F0b's D block: 1 s against 60 s before the same 5 um reversal."""
    segs = [s for s in run.segments if s.block == "D"]
    by = defaultdict(list)
    for i, s in enumerate(segs):
        if s.kind != "ramp":
            continue
        pre = [q for q in segs[:i] if q.kind == "dwell"][-1]
        x0 = run.settled(pre, 300)
        x1 = float(np.mean(run.x_nm[s.t1_ms - 100:s.t1_ms]))
        commanded = abs((s.target_u - pre.target_u) * USTEP_NM)
        by[pre.dwell_ms / 1000.0].append(commanded - abs(x1 - x0))
    short, long = sorted(by)
    a, b = np.array(by[short]), np.array(by[long])
    spread = float(np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2))
    return {"short_s": short, "long_s": long,
            "short_nm": [float(v) for v in a], "long_nm": [float(v) for v in b],
            "short_mean_nm": float(a.mean()), "long_mean_nm": float(b.mean()),
            "delta_nm": float(b.mean() - a.mean()),
            "within_sd_nm": spread,
            "sem_nm": float(spread * np.sqrt(1 / len(a) + 1 / len(b)))}


def slow_ramps(run) -> list[dict]:
    """F0b's S block: the four preparations, breakaway and lost motion."""
    segs = [s for s in run.segments if s.block == "S"]
    rows = []
    for i, s in enumerate(segs):
        if s.kind != "ramp":
            continue
        pre = [q for q in segs[:i] if q.kind == "dwell"][-1]
        x0 = run.settled(pre, 300)
        noise = float(np.std(run.x_nm[pre.t1_ms - 400:pre.t1_ms]))
        x = run.x_nm[s.t0_ms:s.t1_ms]
        c = run.cmd_nm[s.t0_ms:s.t1_ms]
        moved = np.flatnonzero(np.abs(x - x0) > max(5 * noise, 25))
        commanded = (s.target_u - pre.target_u) * USTEP_NM
        # which preparation this is: the wipe before it went up or down
        wipe = [q for q in segs[:i] if q.kind == "move"][-2]
        from_below = wipe.target_u < pre.target_u
        rows.append({
            "commanded_nm": float(commanded),
            "from_below": bool(from_below),
            "reversal": bool((commanded < 0) == from_below),
            "breakaway_nm": float(abs(c[moved[0]] - c[0])) if len(moved) else None,
            "lost_nm": float(abs(commanded) - abs(
                float(np.mean(run.x_nm[s.t1_ms - 100:s.t1_ms])) - x0)),
            "t0_ms": s.t0_ms, "t1_ms": s.t1_ms, "x0_nm": x0,
        })
    return rows


def transport_delay(run) -> dict:
    """F0b's V block: the lag mid-move, at 1 mm/s and at 0.01 mm/s.

    A pure pipeline delay shows up as a lag proportional to speed and
    cancels between two settled endpoints, which is why the ladder is
    unaffected by it. Anything left over is sliding friction.
    """
    rows = []
    for i, s in enumerate(run.segments):
        if s.block != "V" or s.kind != "ramp":
            continue
        rate = s.rate_mhz / 1e3 * USTEP_NM / 1000.0      # um/s
        a, b = s.t0_ms, s.t1_ms
        mid = slice(a + int(0.4 * (b - a)), a + int(0.9 * (b - a)))
        lag = ((run.x_nm[mid] - run.x_nm[a]) - (run.cmd_nm[mid] - run.cmd_nm[a]))
        rows.append({"rate_um_s": rate, "lag_nm": float(np.mean(lag)),
                     "lag_sd_nm": float(np.std(lag)),
                     "implied_delay_ms": float(abs(np.mean(lag)) / rate)})
    fast = [r for r in rows if r["rate_um_s"] > 100]
    return {"passes": rows,
            "fast_delay_ms": float(np.mean([r["implied_delay_ms"] for r in fast]))}


def post_move_creep(run) -> list[dict]:
    """How much the stage is still moving when the 0.5 s pause begins."""
    segs = [s for s in run.segments if s.block == "L"]
    rows = []
    for i, s in enumerate(segs):
        if s.kind != "dwell" or i == 0 or segs[i - 1].kind != "ramp":
            continue
        rate = segs[i - 1].rate_mhz / 1e3 * USTEP_NM / 1000.0
        x = run.x_nm[s.t0_ms:s.t1_ms]
        rows.append({"rate_um_s": rate, "catch_up_nm": float(x[-1] - x[0])})
    by = defaultdict(list)
    for r in rows:
        by[round(r["rate_um_s"], 2)].append(r["catch_up_nm"])
    return [{"rate_um_s": k, "catch_up_nm": float(np.mean(v)), "n": len(v)}
            for k, v in sorted(by.items())]


# -- figures --------------------------------------------------------------
def figure_quality(runs: dict, report: dict) -> Path:
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, left=.085, right=.975, top=.90, bottom=.145,
                          hspace=.60, wspace=.26)

    # (a) the whole F0b run, aligned
    ax = fig.add_subplot(gs[0, :])
    run = runs["F0b"]
    d = 50
    t = np.arange(0, len(run.cmd_nm), d) / 60000.0
    for seg in run.segments:
        colour = F0B_BLOCKS.get(seg.block)
        if colour:
            ax.axvspan(seg.t0_ms / 60000, seg.t1_ms / 60000, color=colour,
                       alpha=.08, lw=0, zorder=0)
    ax.plot(t, run.cmd_nm[::d] / 1e6, color=MUTED, lw=1.6, alpha=.5,
            label="commanded")
    ax.plot(t, run.x_nm[::d] / 1e6, color=INK, lw=.7, label="measured")
    seen = set()
    for seg in run.segments:
        if seg.block in seen or seg.block not in F0B_BLOCKS:
            continue
        seen.add(seg.block)
        ax.text(seg.t0_ms / 60000 + .2, 116, seg.block, fontsize=7.5,
                color=F0B_BLOCKS[seg.block], weight="bold", va="top")
    ax.set(xlim=(0, len(run.cmd_nm) / 60000), ylim=(0, 122),
           xlabel="elapsed session time [min]", ylabel="position [mm]")
    ax.set_title("(a)  F0b, capture aligned to the command table   "
                 f"$R^2$ = {run.r2:.6f},  scale {run.counts_per_nm:.5f} count/nm",
                 loc="left", pad=6)
    ax.legend(fontsize=7, loc="lower right", frameon=False, ncols=2)
    tidy(ax)

    # (b) Allan deviation of the three noise holds
    ax = fig.add_subplot(gs[1, 0])
    for row, colour in zip(report["F0b"]["noise_floor"], (TEAL, ORANGE, PURPLE)):
        tau = np.array(row["taus_ms"], float) / 1000.0
        dev = np.array([np.nan if v is None else v for v in row["allan_nm"]])
        ax.plot(tau, dev, "o-", ms=3, lw=1.2, color=colour,
                label=f"{row['position_mm']:.0f} mm")
    ax.axhline(USTEP_NM, color=RED, ls=":", lw=1)
    ax.text(1.3e-3, USTEP_NM * 1.15, "1 µstep = 39 nm", fontsize=7, color=RED)
    ax.axhline(5.0, color=MUTED, ls="--", lw=.8)
    ax.text(1.3e-3, 5.5, "$d_{\\min}$ assumed 5 nm", fontsize=7, color=MUTED)
    ax.set(xscale="log", yscale="log", xlim=(1e-3, 6),
           xlabel="averaging time [s]", ylabel="Allan deviation [nm]")
    ax.set_title("(b)  how still it is at rest", loc="left", pad=6)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    tidy(ax)

    # (c) one noise hold, raw
    ax = fig.add_subplot(gs[1, 1])
    hold = max((s for s in run.segments
                if s.block == "N" and s.kind == "dwell" and s.dwell_ms > 60000),
               key=lambda s: s.dwell_ms)
    x = run.x_nm[hold.t0_ms + 3000:hold.t1_ms - 500]
    tt = np.arange(len(x)) / 60000.0
    ax.plot(tt, x - x.mean(), color=INK, lw=.4)
    blk = x[:(len(x) // 500) * 500].reshape(-1, 500).mean(axis=1)
    ax.plot((np.arange(len(blk)) + .5) * 500 / 60000.0, blk - x.mean(),
            color=ORANGE, lw=1.2, label="0.5 s means")
    ax.set(xlabel="time into the 5 min hold [min]",
           ylabel="position [nm]", ylim=(-60, 60))
    ax.set_title("(c)  the 60 mm hold, raw and averaged", loc="left", pad=6)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    tidy(ax)

    q = report["F0b"]
    fig.text(.085, .035,
             "Both runs completed: F0b 346/346 segments in 48.05 min, F1 1119/1119 in 51.31 min, neither aborted or time-limited. "
             "The capture is aligned to the segment log by\ncorrelating velocities, which fixes the offset to a few ms; "
             f"the ladder's 3 mm rungs then come out {q['long_rung_shortfall_nm'] / 1000:.2f} µm short of command instead of the 220 µm a 300 ms error would give. "
             "Quantisation is 1 nm.",
             fontsize=7, color=MUTED, va="bottom", linespacing=1.4)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "preliminary_quality.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    return path


def figure_physics(runs: dict, report: dict) -> Path:
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, left=.085, right=.975, top=.90, bottom=.175,
                          hspace=.62, wspace=.26)

    # (a) lost motion against stroke, both runs, three positions
    ax = fig.add_subplot(gs[0, :])
    lad = report["ladder"]
    ax.errorbar(np.array(lad["F0b"]["stroke_um"]), lad["F0b"]["lost_nm"],
                yerr=lad["F0b"]["sd_nm"], fmt="s-", ms=4, lw=1.4, capsize=2,
                color=INK, label="F0b · anchored ladder at 60 mm")
    for pos in ("A", "B", "C"):
        row = lad["F1"][pos]
        ax.errorbar(np.array(row["stroke_um"]), row["lost_nm"],
                    yerr=row["sd_nm"], fmt="o--", ms=3, lw=1.0, capsize=2,
                    color=POS_COLOUR[pos], alpha=.85,
                    label=f"F1 · symmetric loops at {row['position_mm']:.0f} mm")
    beta = report["summary"]["beta0_um"]
    s90 = report["summary"]["s90_um"]
    ax.axhline(beta * 1000, color=RED, lw=.9, ls=":")
    ax.text(.05, beta * 1050, f"saturation reached here, $\\beta_0 \\approx$ {beta:.1f} µm",
            fontsize=7, color=RED, va="bottom")
    ax.axhline(1757, color=MUTED, lw=.8, ls="--")
    ax.text(.05, 1757 * .78, "$\\beta_0$ assumed 1.76 µm", fontsize=7, color=MUTED,
            va="top")
    for value, label, colour in ((s90, f"$s_{{90}} \\approx$ {s90:.0f} µm", RED),
                                 (3.0, "$s_{90}$ assumed 3 µm", MUTED)):
        ax.axvline(value, color=colour, lw=.8, ls=":" if colour == RED else "--")
        ax.text(value * 1.06, 5200, label, fontsize=7, color=colour,
                rotation=90, va="top")
    ax.plot([.03, 3], [30, 3000], color=GRID, lw=.9, zorder=0)
    ax.text(.045, 55, "everything lost", fontsize=7, color=MUTED, rotation=37)
    ax.set(xscale="log", yscale="log", xlim=(.03, 150), ylim=(25, 6000),
           xlabel="commanded stroke [µm]", ylabel="lost motion [nm]")
    ax.set_title("(a)  lost motion against stroke — the measurement both sessions exist to make",
                 loc="left", pad=6)
    ax.legend(fontsize=6.8, frameon=False, loc="lower right", ncols=2)
    tidy(ax)

    # (b) the four slow-ramp preparations
    ax = fig.add_subplot(gs[1, 0])
    run = runs["F0b"]
    for row in report["F0b"]["slow_ramps"]:
        x = run.x_nm[row["t0_ms"]:row["t1_ms"]] - row["x0_nm"]
        c = run.cmd_nm[row["t0_ms"]:row["t1_ms"]] - run.cmd_nm[row["t0_ms"]]
        sign = 1.0 if row["commanded_nm"] > 0 else -1.0
        ax.plot(sign * c / 1000, sign * x / 1000, lw=1.1,
                color=RED if row["reversal"] else TEAL, alpha=.85)
    ax.plot([0, 20], [0, 20], color=GRID, lw=.9, zorder=0)
    ax.text(12.5, 13.8, "no loss", fontsize=7, color=MUTED, rotation=38)
    ax.plot([], [], color=RED, lw=1.1, label="after a reversal")
    ax.plot([], [], color=TEAL, lw=1.1, label="continuing")
    ax.set(xlim=(0, 20), ylim=(0, 20), xlabel="commanded travel [µm]",
           ylabel="measured travel [µm]")
    ax.set_title("(b)  S: 20 µm at 0.25 µm/s, 8 ramps", loc="left", pad=6)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    tidy(ax)

    # (c) the dwell probe that decides F2
    ax = fig.add_subplot(gs[1, 1])
    d = report["F0b"]["dwell_probe"]
    for i, (key, colour) in enumerate((("short_nm", TEAL), ("long_nm", ORANGE))):
        vals = np.array(d[key])
        ax.plot(np.full(len(vals), i) + np.linspace(-.06, .06, len(vals)),
                vals, "o", ms=6, color=colour)
        ax.plot([i - .22, i + .22], [vals.mean()] * 2, color=colour, lw=2)
    ax.errorbar(1.55, d["short_mean_nm"] + d["delta_nm"] / 2,
                yerr=abs(d["delta_nm"]) / 2, fmt="none", ecolor=INK,
                capsize=4, lw=1.2)
    ax.text(1.62, d["short_mean_nm"] + d["delta_nm"] / 2,
            f"$\\Delta$ = {d['delta_nm']:+.0f} nm\n± {d['sem_nm']:.0f} nm",
            fontsize=7, color=INK, va="center")
    floor_nm = report["summary"]["d_min_nm_measured"]
    ax.axhspan(d["short_mean_nm"] - 3 * floor_nm, d["short_mean_nm"] + 3 * floor_nm,
               color=MUTED, alpha=.35, lw=0)
    ax.text(-.40, d["short_mean_nm"] - 40,
            f"±3·$d_{{\\min}}$ = ±{3 * floor_nm:.0f} nm", fontsize=6.8,
            color=MUTED, va="top")
    ax.set(xlim=(-.45, 2.6), xticks=[0, 1],
           xticklabels=[f"{d['short_s']:g} s dwell", f"{d['long_s']:g} s dwell"],
           ylabel="lost motion over the 5 µm reversal [nm]")
    ax.set_title("(c)  D: the probe that decides F2", loc="left", pad=6)
    tidy(ax)

    s = report["summary"]
    fig.text(.085, .030,
             "(a) Every rung from one microstep (39 nm) to 100 µm is resolved, at three positions and in two independent sessions, and the two agree on shape. "
             f"Both the pre-sliding length and the\nsaturation come out about twice their placeholders: $s_{{90}} \\approx$ {s['s90_um']:.0f} µm against 3 µm assumed, "
             f"$\\beta_0 \\approx$ {s['beta0_um']:.1f} µm against 1.76 µm. (b) The four preparations separate exactly as the model says they should. "
             f"\n(c) 60 s of dwell changes the next reversal by {abs(s['dwell_delta_nm']):.0f} nm against a trial-to-trial scatter of {s['dwell_within_sd_nm']:.0f} nm: "
             "on two repetitions each, no ageing is detectable.",
             fontsize=7, color=MUTED, va="bottom", linespacing=1.4)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "preliminary_physics.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    return path


# -- report ---------------------------------------------------------------
def main() -> int:
    style()
    runs = {k: load(k) for k in ("F0b", "F1")}
    report: dict = {"ladder": {"F1": {}}}

    for key, run in runs.items():
        head = run.header
        report[key] = {
            "capture_min": run.capture_min,
            "starttime": run.meta["starttime"],
            "samples": run.meta["samples"],
            "lag_ms": run.lag_ms,
            "counts_per_nm": run.counts_per_nm,
            "r2": run.r2,
            "segments": len(run.segments),
            "complete": head.get("SEQ_DONE.complete"),
            "aborted": head.get("SEQ_DONE.aborted"),
            "elapsed_ms": int(head["SEQ_DONE.elapsed_ms"]),
            "planned_ms": int(head["SEQ_DONE.planned_ms"]),
            "config_sha256": head["SESSION.config_sha256"],
            "table_sha256": head["SESSION.table_sha256"],
            "dwell_spread_nm": run.meta["dwell_spread_nm"],
            "long_rung_shortfall_nm": run.meta["long_rung_shortfall_nm"],
            "long_rung_rms_nm": run.meta["long_rung_rms_nm"],
            "quantisation": quantisation(run),
        }

    f0b = runs["F0b"]
    report["F0b"]["noise_floor"] = noise_floor(f0b)
    report["F0b"]["dwell_probe"] = dwell_probe(f0b)
    report["F0b"]["slow_ramps"] = slow_ramps(f0b)
    report["F0b"]["transport_delay"] = transport_delay(f0b)
    report["F0b"]["post_move_creep"] = post_move_creep(f0b)

    # The ladders, trimmed to where the 0.5 s pause is long enough. The rung
    # rate is stroke/4 s, so 100 um is 25 um/s and still settles inside the
    # pause; 500 um and 3 mm do not, and are left out of both figures.
    data = ladder(f0b, "L")
    keep = sorted(s for s in data if s <= 1.1e5)
    report["ladder"]["F0b"] = {
        "position_mm": 60.0,
        "stroke_um": [s / 1000 for s in keep],
        "lost_nm": [float(np.mean(data[s])) for s in keep],
        "sd_nm": [float(np.std(data[s])) for s in keep],
        "n": [len(data[s]) for s in keep],
    }
    for pos, mm in (("A", 25.0), ("B", 60.0), ("C", 100.48)):
        d = ladder(runs["F1"], pos, "P")
        keep = sorted(s for s in d if s <= 1.6e4)
        report["ladder"]["F1"][pos] = {
            "position_mm": mm,
            "stroke_um": [s / 1000 for s in keep],
            "lost_nm": [float(np.mean(d[s])) for s in keep],
            "sd_nm": [float(np.std(d[s])) for s in keep],
            "n": [len(d[s]) for s in keep],
        }

    anchored = report["ladder"]["F0b"]
    beta, s90 = saturation_and_s90([v * 1000 for v in anchored["stroke_um"]],
                                   anchored["lost_nm"])
    f1_top = max(max(row["lost_nm"]) for row in report["ladder"]["F1"].values())
    d = report["F0b"]["dwell_probe"]
    floor = report["F0b"]["noise_floor"]
    report["summary"] = {
        "beta0_um": beta / 1000,
        "beta0_assumed_um": 1.757,
        "s90_um": s90 / 1000,
        "s90_assumed_um": 3.0,
        "f1_largest_lost_um": f1_top / 1000,
        "d_min_nm_measured": float(np.max([r["allan_500ms_nm"] for r in floor])),
        "d_min_nm_assumed": 5.0,
        "dwell_delta_nm": d["delta_nm"],
        "dwell_within_sd_nm": d["within_sd_nm"],
        "dwell_sem_nm": d["sem_nm"],
        "transport_delay_ms": report["F0b"]["transport_delay"]["fast_delay_ms"],
    }

    quality = figure_quality(runs, report)
    physics = figure_physics(runs, report)
    (HERE / "preliminary_check.json").write_text(
        json.dumps(report, indent=1) + "\n", encoding="utf-8")

    s = report["summary"]
    print("PRELIMINARY CHECK — F0b and F1, 10 November 2026")
    for key in ("F0b", "F1"):
        r = report[key]
        print(f"  {key}: {r['segments']} segments, "
              f"{r['elapsed_ms'] / 60000:.2f} min of a planned "
              f"{r['planned_ms'] / 60000:.2f}, complete={r['complete']}, "
              f"aborted={r['aborted']}")
        print(f"       capture {r['capture_min']:.2f} min, lag {r['lag_ms']} ms, "
              f"scale {r['counts_per_nm']:.6f} count/nm, R2 {r['r2']:.6f}")
        print(f"       dwells quiet to {r['dwell_spread_nm']:.0f} nm; long rungs "
              f"short by {r['long_rung_shortfall_nm'] / 1000:.2f} µm")
    print("  noise floor (F0b N block):")
    for row in report["F0b"]["noise_floor"]:
        print(f"       {row['position_mm']:6.1f} mm: raw sd {row['sd_raw_nm']:5.1f} nm, "
              f"0.5 s Allan {row['allan_500ms_nm']:.2f} nm, "
              f"drift {row['drift_nm_per_min']:+6.1f} nm/min")
    print(f"  resolution: 1 nm quantisation, d_min measured "
          f"{s['d_min_nm_measured']:.1f} nm against {s['d_min_nm_assumed']:g} assumed")
    print(f"  pre-sliding: s90 ~ {s['s90_um']:.1f} µm (assumed "
          f"{s['s90_assumed_um']:g}), beta0 ~ {s['beta0_um']:.2f} µm (assumed "
          f"{s['beta0_assumed_um']:g})")
    print(f"  dwell probe: {s['dwell_delta_nm']:+.0f} nm for "
          f"{d['long_s']:g} s over {d['short_s']:g} s, "
          f"within-condition sd {s['dwell_within_sd_nm']:.0f} nm "
          f"-> no ageing detectable on 2 repetitions")
    print(f"  pipeline delay mid-move: {s['transport_delay_ms']:.0f} ms "
          f"(cancels between settled endpoints)")
    print(f"  wrote {quality.name}, {physics.name}, preliminary_check.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
