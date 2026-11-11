#!/usr/bin/env python3
"""Process the F0b run into the numbers the results figures and the scales
card are drawn from.

Follows procedure section 5. The one deliberate departure is how lost motion
is defined; see `ladder()`.

    python3 f0b_process.py

writes

    f0b/f0b_noise.json      sigma(tau) per hold, d_min, drift rate
    f0b/f0b_ladder.csv      one row per half-cycle: s, pass, kind, dx, tau, l
    f0b/f0b_summary.json    the scales card, plus every number in the captions

Every window mean is read with a 50 ms trailing guard. The ESP and the scope
free-run on their own crystals, and by the end of a 48 min session they have
slipped far enough that the last samples of three of the twelve reference
record windows had already caught the start of the next move: their standard
deviation was 2.1-2.9 um without the guard and 4-8 nm with it.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from run_data import USTEP_NM, load, Run, Segment

HERE = Path(__file__).resolve().parent
OUT = HERE / "f0b"

GUARD_MS = 50           # dropped from the end of every window mean
WINDOW_MS = 300         # length of a window mean
X0_MM = 60.0
BETA_E5_UM = 1.757      # saturation lost motion from E2-E5, for comparison


# -- windows --------------------------------------------------------------
def window(run: Run, seg: Segment, length: int = WINDOW_MS) -> slice:
    stop = min(seg.t1_ms - GUARD_MS, len(run.x_nm))
    return slice(max(seg.t0_ms, stop - length), stop)


def mean_at(run: Run, seg: Segment, length: int = WINDOW_MS) -> float:
    w = window(run, seg, length)
    return float(np.mean(run.x_nm[w])) if w.stop > w.start else np.nan


def sd_at(run: Run, seg: Segment, length: int = WINDOW_MS) -> float:
    w = window(run, seg, length)
    return float(np.std(run.x_nm[w])) if w.stop > w.start else np.nan


# -- 5.1 drift ------------------------------------------------------------
@dataclass
class Anchor:
    block: str
    position_mm: float
    t_min: float
    up_nm: float
    down_nm: float
    mean_nm: float
    reversal_nm: float
    window_sd_nm: float


def anchors(run: Run) -> list[Anchor]:
    """The reference pairs: the up/down mean at 15 and 105 mm, three times.

    Averaging the two approach directions cancels the reversal error and
    leaves the slow change in the stage, which is what the drift correction
    is for. The difference between them is that reversal error, and it is a
    free cross-check on the ladder's saturation level.
    """
    marks = [s for s in run.segments if s.kind == "mark"]
    rows, bucket = [], {}
    for i, seg in enumerate(marks):
        before = [q for q in run.segments if q.kind == "move" and q.t1_ms <= seg.t0_ms]
        rising = before[-2].target_u < seg.target_u
        key = (seg.block, round(seg.target_u * USTEP_NM / 1e6, 3))
        entry = bucket.setdefault(key, {"t": [], "sd": []})
        entry["up" if rising else "down"] = mean_at(run, seg)
        entry["t"].append(seg.t0_ms / 60000.0)
        entry["sd"].append(sd_at(run, seg))
    for (block, mm), e in bucket.items():
        if "up" not in e or "down" not in e:
            continue
        rows.append(Anchor(block, mm, float(np.mean(e["t"])), e["up"], e["down"],
                           (e["up"] + e["down"]) / 2, e["down"] - e["up"],
                           float(np.max(e["sd"]))))
    return sorted(rows, key=lambda a: a.t_min)


def drift_model(rows: list[Anchor]):
    """Piecewise-linear drift in time, linear in position between 15 and 105.

    Each reference is referred to its own first visit, so the model returns
    how much the stage has moved since the start of the session at a given
    place and time.
    """
    by_pos: dict[float, list[Anchor]] = {}
    for a in rows:
        by_pos.setdefault(a.position_mm, []).append(a)
    tracks = {}
    for mm, items in by_pos.items():
        items = sorted(items, key=lambda a: a.t_min)
        tracks[mm] = (np.array([a.t_min for a in items]),
                      np.array([a.mean_nm - items[0].mean_nm for a in items]))
    places = sorted(tracks)

    def drift(t_min, x_mm: float = X0_MM) -> np.ndarray:
        values = [np.interp(t_min, *tracks[p]) for p in places]
        if len(places) == 1:
            return np.asarray(values[0])
        return np.interp(x_mm, places, [float(np.atleast_1d(v)[0])
                                        for v in values]) * 0 + \
            np.asarray([np.interp(x_mm, places, col)
                        for col in np.array(values).T]) if np.ndim(t_min) else \
            float(np.interp(x_mm, places, [float(v) for v in values]))

    return drift, tracks


# -- 5.2 noise floor ------------------------------------------------------
def overlapping_allan(x: np.ndarray, taus_ms) -> np.ndarray:
    """Allan deviation of the window means, after a linear detrend."""
    t = np.arange(len(x)) / 1000.0
    detrended = x - np.polyval(np.polyfit(t, x, 1), t)
    out = []
    for m in taus_ms:
        k = len(detrended) // m
        if k < 3:
            out.append(np.nan)
            continue
        block = detrended[:k * m].reshape(k, m).mean(axis=1)
        out.append(float(np.sqrt(0.5 * np.mean(np.diff(block) ** 2))))
    return np.asarray(out)


TAUS_MS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]


def noise(run: Run) -> dict:
    holds = []
    for seg in run.segments:
        if seg.block != "N" or seg.kind != "dwell" or seg.dwell_ms < 60000:
            continue
        x = run.x_nm[seg.t0_ms + 3000:seg.t1_ms - GUARD_MS]
        t = np.arange(len(x)) / 1000.0
        slope = float(np.polyfit(t, x, 1)[0])
        dev = overlapping_allan(x, TAUS_MS)
        holds.append({
            "name": {15.0: "N-end1", 60.0: "N-mid", 105.0: "N-end2"}.get(
                round(seg.target_u * USTEP_NM / 1e6, 1), "N"),
            "position_mm": seg.target_u * USTEP_NM / 1e6,
            "minutes": len(x) / 60000.0,
            "taus_s": [m / 1000 for m in TAUS_MS],
            "allan_nm": [None if np.isnan(v) else float(v) for v in dev],
            "sigma_500ms_nm": float(dev[TAUS_MS.index(500)]),
            "drift_nm_per_s": slope,
            "sd_raw_nm": float(np.std(x)),
        })
    mid = next(h for h in holds if h["name"] == "N-mid")
    sigma = mid["sigma_500ms_nm"]
    return {"holds": holds, "sigma_500ms_nm": sigma,
            "d_min_nm": 3 * np.sqrt(2) * sigma,
            "primary": "N-mid",
            "position_spread_nm": float(
                max(h["sigma_500ms_nm"] for h in holds)
                - min(h["sigma_500ms_nm"] for h in holds))}


def spectrum(run: Run, nperseg: int = 10000) -> dict:
    """Welch displacement spectrum of the long hold, Hann, 50 % overlap."""
    seg = max((s for s in run.segments
               if s.block == "N" and s.kind == "dwell" and s.dwell_ms >= 300000),
              key=lambda s: s.dwell_ms)
    x = run.x_nm[seg.t0_ms + 3000:seg.t1_ms - GUARD_MS]
    t = np.arange(len(x)) / 1000.0
    x = x - np.polyval(np.polyfit(t, x, 1), t)
    win = np.hanning(nperseg)
    step = nperseg // 2
    pieces = [x[i:i + nperseg] * win
              for i in range(0, len(x) - nperseg + 1, step)]
    spec = np.mean([np.abs(np.fft.rfft(p)) ** 2 for p in pieces], axis=0)
    scale = 2.0 / (1000.0 * np.sum(win ** 2))      # one-sided PSD, nm^2/Hz
    psd = spec * scale
    freq = np.fft.rfftfreq(nperseg, 1 / 1000.0)
    asd = np.sqrt(psd)
    band = (freq > 1) & (freq < 400)
    order = np.argsort(asd[band])[::-1][:4]
    return {"freq_hz": freq.tolist(), "asd_nm_per_rthz": asd.tolist(),
            "lines_hz": sorted(float(f) for f in freq[band][order])}


# -- 5.3 the ladder -------------------------------------------------------
@dataclass
class Half:
    stroke_um: float
    pass_name: str          # increasing / decreasing
    index: int              # 1..4 within the stroke
    kind: str               # continuation / reversal
    direction: int
    dx_um: float            # measured travel
    tau: float              # fraction arriving
    lost_travel_um: float   # s - dx, the procedure's definition
    loop_width_um: float    # width at the stroke midpoint, geometry cancelled
    t_min: float


def branch(run: Run, seg: Segment) -> tuple[np.ndarray, np.ndarray]:
    sl = slice(seg.t0_ms, seg.t1_ms)
    return run.cmd_nm[sl], run.x_nm[sl]


def ladder(run: Run, drift=None) -> list[Half]:
    """One row per half-cycle of the L block.

    LOST MOTION IS TAKEN TWO WAYS, and which one the figures use is a
    decision the data forced.

    The results-plot instructions ask for the width of the hysteresis loop
    read at one commanded position, because the screw's geometric error is
    then the same on both branches and cancels exactly. That is the better
    definition in general, and it is computed here as `loop_width_um`, at
    the midpoint of the stroke.

    It cannot be used for this ladder. Each rung runs at stroke / 4 s, so
    every stroke is measured at a different speed, from 10 nm/s to
    750 um/s. The carriage lags the screw in proportion to speed -- a lag
    that relaxes within about 15 ms of the command stopping, so it is
    mechanical and not a measurement delay -- and that lag enters the loop
    width twice, once per branch. The loop width therefore keeps growing
    with the rung rate instead of saturating: 2.0 um at the 20 um rung,
    3.2 um at 100 um and 41 um at 3 mm, where ell = s - dx stays near 2 um
    throughout. It cannot give ell_sat or s90.

    The procedure's own definition, ell = s - dx between the settled
    window means at the two turning points, is free of that: the rate lag
    has relaxed before either window is read. It does carry the geometric
    error, but over a stroke below 100 um that is under 0.6 um and below
    10 um it is under 0.06 um (procedure 5.1), well under ell itself. The
    figures and the scales card use it, and the strokes above 100 um are
    excluded for the separate reason that the 0.5 s pause is too short for
    them to settle.
    """
    segs = [s for s in run.segments if s.block == "L"]
    moves = [(i, s) for i, s in enumerate(segs) if s.kind in ("step", "ramp")]
    rows: list[Half] = []
    last_dir, seen = +1, {}
    for n, (i, seg) in enumerate(moves):
        pause = segs[i + 1] if i + 1 < len(segs) else None
        if pause is None or pause.kind != "dwell":
            continue
        start = segs[i - 1]
        x0 = mean_at(run, start)
        x1 = mean_at(run, pause)
        travel = (seg.target_u - start.target_u) * USTEP_NM
        if travel == 0:
            continue
        direction = 1 if travel > 0 else -1
        stroke = abs(travel)
        if drift is not None:
            x0 -= drift(start.t1_ms / 60000.0)
            x1 -= drift(pause.t1_ms / 60000.0)
        dx = abs(x1 - x0)
        key = round(stroke, 3)
        seen[key] = seen.get(key, 0) + 1
        index = (seen[key] - 1) % 4 + 1
        pass_name = "increasing" if seen[key] <= 4 else "decreasing"

        # the loop width at the midpoint, against the opposite branch of the
        # same cycle; only defined once both branches of a cycle exist
        width = np.nan
        if index >= 2:
            prev_seg = moves[n - 1][1]
            u_mid = (min(seg.target_u, prev_seg.target_u)
                     + max(seg.target_u, prev_seg.target_u)) / 2 * USTEP_NM
            ua, xa = branch(run, prev_seg)
            ub, xb = branch(run, seg)
            if len(ua) > 2 and len(ub) > 2:
                fa = np.interp(u_mid, ua, xa) if ua[-1] > ua[0] else \
                    np.interp(u_mid, ua[::-1], xa[::-1])
                fb = np.interp(u_mid, ub, xb) if ub[-1] > ub[0] else \
                    np.interp(u_mid, ub[::-1], xb[::-1])
                width = abs(fb - fa)
        rows.append(Half(stroke / 1000, pass_name, index,
                         "continuation" if direction == last_dir else "reversal",
                         direction, dx / 1000, dx / stroke,
                         (stroke - dx) / 1000, width / 1000,
                         seg.t0_ms / 60000.0))
        last_dir = direction
    return rows


def ladder_summary(rows: list[Half], d_min_um: float) -> dict:
    """Saturation, pre-sliding length and the smallest stroke that moves."""
    rev = [r for r in rows if r.kind == "reversal"]
    by: dict[float, list[Half]] = {}
    for r in rev:
        by.setdefault(round(r.stroke_um, 4), []).append(r)
    strokes = np.array(sorted(by))
    width = np.array([float(np.mean([q.lost_travel_um for q in by[s]]))
                      for s in strokes])
    width_sd = np.array([float(np.std([q.lost_travel_um for q in by[s]]))
                         for s in strokes])
    loop = np.array([float(np.nanmean([q.loop_width_um for q in by[s]]))
                     for s in strokes])
    # rungs faster than ~25 um/s do not settle inside the 0.5 s pause
    usable = strokes <= 100.0
    s_u, w_u = strokes[usable], width[usable]
    knee = s_u >= 19.0
    sat = float(np.mean(w_u[knee])) if knee.any() else float(np.max(w_u))
    rise = s_u <= 19.0
    s90 = float(np.interp(0.9 * sat, w_u[rise], s_u[rise]))

    moved = [r for r in rev if r.dx_um > d_min_um]
    smallest = min((r.stroke_um for r in moved), default=np.nan)
    inc = {round(r.stroke_um, 4): r.lost_travel_um for r in rev
           if r.pass_name == "increasing"}
    dec = {round(r.stroke_um, 4): r.lost_travel_um for r in rev
           if r.pass_name == "decreasing"}
    shared = sorted(set(inc) & set(dec))
    hysteresis = float(np.mean([abs(inc[s] - dec[s]) for s in shared])) if shared else np.nan
    return {
        "stroke_um": strokes.tolist(), "lost_um": width.tolist(),
        "lost_sd_um": width_sd.tolist(), "loop_width_um": loop.tolist(),
        "saturation_um": sat, "s90_um": s90,
        "beta_e5_um": BETA_E5_UM,
        "smallest_moving_stroke_um": smallest,
        "smallest_moving_stroke_usteps": smallest * 1000 / USTEP_NM,
        "pass_difference_um": hysteresis,
        "usable_above_um": 100.0,
    }


# -- 5.4 slow ramps -------------------------------------------------------
@dataclass
class Ramp:
    name: str               # S1..S4
    repeat: int
    direction: int
    reversal: bool
    t0_ms: int
    t1_ms: int
    rate_um_s: float


def slow_ramps(run: Run) -> list[Ramp]:
    """S1 reverse down, S2 reverse up, S3 continue up, S4 continue down."""
    segs = [s for s in run.segments if s.block == "S"]
    out, order = [], 0
    for i, seg in enumerate(segs):
        if seg.kind != "ramp":
            continue
        pre = [q for q in segs[:i] if q.kind == "dwell"][-1]
        wipe = [q for q in segs[:i] if q.kind == "move"][-2]
        from_below = wipe.target_u < pre.target_u
        direction = 1 if seg.target_u > pre.target_u else -1
        reversal = (direction < 0) == from_below
        name = {(-1, True): "S1", (1, True): "S2",
                (1, False): "S3", (-1, False): "S4"}[(direction, reversal)]
        out.append(Ramp(name, order // 4 + 1, direction, reversal,
                        seg.t0_ms, seg.t1_ms,
                        seg.rate_mhz / 1e3 * USTEP_NM / 1000))
        order += 1
    return out


def ramp_curve(run: Run, r: Ramp) -> tuple[np.ndarray, np.ndarray]:
    """Lag against commanded travel since the ramp start, both in um."""
    sl = slice(r.t0_ms, r.t1_ms)
    u = run.cmd_nm[sl] - run.cmd_nm[r.t0_ms]
    x = run.x_nm[sl] - run.x_nm[r.t0_ms]
    return np.abs(u) / 1000, (x - u) * r.direction / 1000


def ramp_mean(run: Run, ramps: list[Ramp], name: str,
              grid: np.ndarray) -> np.ndarray:
    """Mean lag of the repetitions of one ramp case, on a common grid."""
    curves = []
    for r in ramps:
        if r.name != name:
            continue
        s, e = ramp_curve(run, r)
        curves.append(np.interp(grid, s, e))
    return np.mean(curves, axis=0)


def ramp_pair(run: Run, ramps: list[Ramp], a: str, b: str,
              grid: np.ndarray) -> np.ndarray:
    """Mean of a minus mean of b on a common travel grid; isolates friction.

    The two ramps cover the same commanded positions in the same direction,
    so the geometric and microstep error are identical on both and cancel.
    """
    def mean_of(name):
        curves = []
        for r in ramps:
            if r.name != name:
                continue
            s, e = ramp_curve(run, r)
            curves.append(np.interp(grid, s, e))
        return np.mean(curves, axis=0)
    return mean_of(a) - mean_of(b)


def deadband(grid: np.ndarray, lag: np.ndarray, d_min_um: float) -> float:
    """Commanded travel before the carriage has moved by more than d_min.

    Read on the mean reversal ramp, where "stuck" is the slope -1 line
    lag = -s, so the arrived displacement is s + lag. The first crossing
    has to stick: at these travels a single sample is +-10 nm, so the test
    is the first travel beyond which the arrival stays above d_min.
    """
    arrived = grid + lag
    above = arrived > d_min_um
    for i in range(len(grid)):
        if above[i] and above[i:i + 50].all():
            return float(grid[i])
    return np.nan


def half_arrival(grid: np.ndarray, lag: np.ndarray) -> float:
    """Travel at which half of the commanded motion has arrived.

    The deadband at d_min is below the sampling here, so this is the number
    worth quoting next to it: the carriage starts moving at once, but only
    reaches half transmission after this much travel.
    """
    arrived = grid + lag
    above = arrived > 0.5 * grid
    for i in range(len(grid)):
        if grid[i] > 0.05 and above[i] and above[i:i + 50].all():
            return float(grid[i])
    return np.nan


def microstep_steps(run: Run, r: Ramp, n: int = 30) -> dict:
    """The first n microsteps: does each one arrive?"""
    dt = 1000.0 / (r.rate_um_s * 1000 / USTEP_NM)      # ms per microstep
    edges = r.t0_ms + np.arange(n + 1) * dt
    means, scatter = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        w = run.x_nm[int(a) + int(0.3 * dt):int(b)]
        means.append(float(np.mean(w)) if len(w) else np.nan)
        scatter.append(float(np.std(w)) if len(w) else np.nan)
    means = np.asarray(means) * r.direction      # along the direction of travel
    inc = np.diff(means)
    return {"dt_ms": dt, "edges_ms": edges.tolist(),
            "step_means_nm": means.tolist(),
            "per_step_increment_nm": inc.tolist(),
            "within_step_sd_nm": scatter,
            "mean_increment_nm": float(np.nanmean(inc))}


def rate_lag(run: Run) -> dict:
    """How far the carriage is behind when a ramp's command stops, by rate.

    Read at the first sample after the last commanded step and again 400 ms
    later, once it has caught up. A measurement delay would not catch up;
    this does, within about 15 ms, so it is the transmission lagging under
    load and not the instrument.
    """
    rows = []
    for seg in run.segments:
        if seg.kind != "ramp" or seg.block not in ("L", "S", "V"):
            continue
        rate = seg.rate_mhz / 1e3 * USTEP_NM / 1000
        after = run.x_nm[seg.t1_ms:seg.t1_ms + 400]
        if len(after) < 300 or rate <= 0:
            continue
        deficit = float(after[-1] - after[0])
        target = 0.63 * abs(deficit)
        hit = np.flatnonzero(np.abs(after - after[0]) > target)
        rows.append({"rate_um_s": rate, "deficit_um": abs(deficit) / 1000,
                     "equivalent_ms": abs(deficit) / rate,
                     "tau_ms": int(hit[0]) if len(hit) else None,
                     "block": seg.block})
    fast = [r for r in rows if r["rate_um_s"] > 1.0]
    if fast:
        v = np.array([r["rate_um_s"] for r in fast])
        d = np.array([r["deficit_um"] for r in fast])
        slope = float(np.polyfit(v, d, 1)[0])
    else:
        slope = np.nan
    taus = [r["tau_ms"] for r in rows if r["tau_ms"]]
    return {"rows": rows, "deficit_per_speed_s": slope,
            "relaxation_ms": float(np.median(taus)) if taus else np.nan}


# -- 5.5 dwell, 5.6 speed -------------------------------------------------
def dwell(run: Run, d_min_um: float) -> dict:
    segs = [s for s in run.segments if s.block == "D"]
    trials = []
    for i, seg in enumerate(segs):
        if seg.kind != "ramp":
            continue
        pre = [q for q in segs[:i] if q.kind == "dwell"][-1]
        x0 = mean_at(run, pre)
        sl = slice(seg.t0_ms, seg.t1_ms - GUARD_MS)
        travel = np.abs(run.cmd_nm[sl] - run.cmd_nm[seg.t0_ms]) / 1000
        lost = travel - np.abs(run.x_nm[sl] - x0) / 1000
        creep = None
        if pre.dwell_ms >= 10000:
            w = run.x_nm[pre.t0_ms:pre.t1_ms - GUARD_MS]
            creep = {"t_s": (np.arange(len(w)) / 1000).tolist(),
                     "x_nm": (w - np.mean(w[:200])).tolist()}
        trials.append({"dwell_s": pre.dwell_ms / 1000.0,
                       "travel_um": travel.tolist(),
                       "lost_um": lost.tolist(),
                       "final_lost_um": float(lost[-1]),
                       "creep": creep})
    short = [t for t in trials if t["dwell_s"] < 10]
    long = [t for t in trials if t["dwell_s"] >= 10]
    a = np.array([t["final_lost_um"] for t in short])
    b = np.array([t["final_lost_um"] for t in long])
    sd = float(np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2))
    return {"trials": trials,
            "short_s": short[0]["dwell_s"], "long_s": long[0]["dwell_s"],
            "short_mean_um": float(a.mean()), "long_mean_um": float(b.mean()),
            "delta_nm": float((b.mean() - a.mean()) * 1000),
            "within_sd_nm": sd * 1000,
            "sem_nm": float(sd * np.sqrt(1 / len(a) + 1 / len(b)) * 1000),
            "d_min_nm": d_min_um * 1000,
            "significant": bool(abs(b.mean() - a.mean()) > 3 * d_min_um
                                and abs(b.mean() - a.mean())
                                > 2 * sd * np.sqrt(1 / len(a) + 1 / len(b)))}


def speed(run: Run) -> dict:
    """Friction lag as half the difference of the up and down passes.

    Taken at the same commanded position, so the geometric error is common
    to both directions and cancels exactly.

    ONLY THE SLOW PASS IS USABLE. Reading a lag part-way through a move
    needs the commanded position at that instant, and the segment log gives
    only the start and end times. The end timestamp carries an overhead of
    up to 26 ms beyond the step train itself, and the global alignment of
    the two clocks is itself flat to tens of ms, so the commanded position
    mid-ramp is uncertain by a few tens of milliseconds of travel. At
    10 um/s that is well under a micron; at 1 mm/s it is tens of microns,
    larger than any friction effect. The fast passes are therefore
    computed and reported but marked unusable, and the speed question is
    left to F3 with the step-train end logged.
    """
    segs = [s for s in run.segments if s.block == "V"]
    passes = []
    for i, seg in enumerate(segs):
        if seg.kind != "ramp":
            continue
        pre = [q for q in segs[:i] if q.kind == "dwell"][-1]
        sl = slice(seg.t0_ms, seg.t1_ms)
        u = run.cmd_nm[sl]
        x = run.x_nm[sl]
        direction = 1 if seg.target_u > pre.target_u else -1
        passes.append({"rate_um_s": seg.rate_mhz / 1e3 * USTEP_NM / 1000,
                       "direction": direction, "u": u, "x": x,
                       "t_min": seg.t0_ms / 60000.0})
    # One absolute window of commanded position, trimmed clear of the ends
    # so neither pass is still accelerating inside it.
    lo = max(min(p["u"].min() for p in passes) + 100e3, 0)
    hi = max(p["u"].max() for p in passes) - 100e3
    grid = np.linspace(lo, hi, 200)
    out = {}
    for rate in sorted({round(p["rate_um_s"], 3) for p in passes}):
        ups, downs = [], []
        for p in passes:
            if round(p["rate_um_s"], 3) != rate:
                continue
            order = np.argsort(p["u"])
            # x - u at the same commanded position: the geometric error is
            # common to both directions and cancels in the difference below.
            resampled = np.interp(grid, p["u"][order],
                                  (p["x"] - p["u"])[order])
            (ups if p["direction"] > 0 else downs).append(resampled)
        if not ups or not downs:
            continue
        lag = (np.mean(downs, axis=0) - np.mean(ups, axis=0)) / 2 / 1000
        out[f"{rate:g}"] = {"position_um": ((grid - lo) / 1000).tolist(),
                            "lag_um": lag.tolist(),
                            "mean_lag_um": float(np.mean(lag)),
                            "passes": len(ups) + len(downs)}
    rates = sorted(out, key=float)
    slow = out[rates[0]]
    # how much commanded position a 26 ms timing uncertainty is worth
    for rate in rates:
        out[rate]["timing_uncertainty_um"] = 0.026 * float(rate)
        out[rate]["usable"] = out[rate]["timing_uncertainty_um"] < 0.5
    return {"by_rate": out, "slow_um_s": float(rates[0]),
            "fast_um_s": float(rates[-1]),
            "slow_lag_um": slow["mean_lag_um"],
            "slow_lag_sd_um": float(np.std(slow["lag_um"])),
            "comparable": False}


# -- A1 validity ----------------------------------------------------------
def validity(run: Run, d_min_nm: float) -> dict:
    windows = [(s.row, s.block, sd_at(run, s))
               for s in run.segments
               if s.kind in ("dwell", "mark") and s.t1_ms - s.t0_ms >= 450]
    sds = np.array([w[2] for w in windows if np.isfinite(w[2])])
    threshold = 50.0
    returns = [(s.t0_ms / 60000.0, s.mscnt1) for s in run.segments
               if s.target_u == round(X0_MM * 1000 * 1000 / USTEP_NM)
               and s.kind in ("dwell", "mark")]
    head = run.header
    mscnt = sorted({m for _t, m in returns})
    checks = [
        ("segments executed", f"{head['SEQ_DONE.segments']} of "
                              f"{head['PLAN.segments']}",
         head["SEQ_DONE.segments"] == head["PLAN.segments"]),
        ("run completed", "complete, not aborted"
         if head["SEQ_DONE.complete"] == "1" and head["SEQ_DONE.aborted"] == "0"
         else "incomplete",
         head["SEQ_DONE.complete"] == "1" and head["SEQ_DONE.aborted"] == "0"),
        ("MSCNT constant at x0", f"{mscnt[0]} throughout" if len(mscnt) == 1
         else f"{len(mscnt)} values", len(mscnt) == 1),
        ("driver faults", "none in DRV_STATUS",
         all((s.status & 0x0FFF0000) in (0x000B0000, 0x000B0000) or
             (s.status & ~0x80000000) == 0x000B0000 for s in run.segments)),
        ("record windows in spec",
         f"{int(np.sum(sds < threshold))}/{len(sds)} under {threshold:g} nm",
         bool(np.sum(sds < threshold) == len(sds))),
        ("scope covers the run",
         f"{run.capture_min:.1f} min / "
         f"{int(head['SEQ_DONE.elapsed_ms']) / 60000:.1f} min",
         run.capture_min > int(head["SEQ_DONE.elapsed_ms"]) / 60000),
    ]
    return {"window_sd_nm": sds.tolist(), "threshold_nm": threshold,
            "mscnt_returns": returns,
            "checks": [{"name": n, "value": v, "pass": bool(p)}
                       for n, v, p in checks]}


# -- assemble -------------------------------------------------------------
def process() -> dict:
    run = load("F0b")
    rows = anchors(run)
    drift, tracks = drift_model(rows)
    nz = noise(run)
    d_min_um = nz["d_min_nm"] / 1000

    halves = ladder(run, drift)
    lad = ladder_summary(halves, d_min_um)
    ramps = slow_ramps(run)
    grid = np.linspace(0, 20, 4001)
    down = ramp_pair(run, ramps, "S1", "S4", grid)
    up = ramp_pair(run, ramps, "S2", "S3", grid)
    rev_down = ramp_mean(run, ramps, "S1", grid)
    rev_up = ramp_mean(run, ramps, "S2", grid)
    dw = dwell(run, d_min_um)
    sp = speed(run)
    val = validity(run, nz["d_min_nm"])

    rl = rate_lag(run)
    steps = microstep_steps(run, next(r for r in ramps if r.name == "S1"))
    increments = np.array(steps["per_step_increment_nm"])
    arrives = float(np.mean(increments > nz["d_min_nm"]))
    # Stick-slip is a sawtooth: steps that deliver nothing, then one that
    # delivers several steps' worth. Scatter alone is not evidence of it.
    stuck = float(np.mean(increments < nz["d_min_nm"]))
    jumps = float(np.mean(increments > 2 * steps["mean_increment_nm"]))

    summary = {
        "run": {"segments": len(run.segments), "r2": run.r2,
                "lag_ms": run.lag_ms, "counts_per_nm": run.counts_per_nm,
                "config_sha256": run.header["SESSION.config_sha256"],
                "table_sha256": run.header["SESSION.table_sha256"]},
        "sigma_500ms_nm": nz["sigma_500ms_nm"],
        "d_min_nm": nz["d_min_nm"],
        "floor_position_spread_nm": nz["position_spread_nm"],
        "saturation_um": lad["saturation_um"],
        "beta_e5_um": BETA_E5_UM,
        "s90_um": lad["s90_um"],
        "deadband_down_um": deadband(grid, rev_down, d_min_um),
        "deadband_up_um": deadband(grid, rev_up, d_min_um),
        "half_arrival_down_um": half_arrival(grid, rev_down),
        "half_arrival_up_um": half_arrival(grid, rev_up),
        "smallest_moving_stroke_um": lad["smallest_moving_stroke_um"],
        "smallest_moving_stroke_usteps": lad["smallest_moving_stroke_usteps"],
        "pass_difference_um": lad["pass_difference_um"],
        "microstep_arrival_fraction": arrives,
        "microstep_mean_increment_nm": steps["mean_increment_nm"],
        "microstep_stuck_fraction": stuck,
        "microstep_jump_fraction": jumps,
        "stick_slip": bool(stuck > 0.10 and jumps > 0.05),
        "stick_slip_spread_nm": float(np.nanmax(increments)
                                      - np.nanmin(increments)),
        "dwell_delta_nm": dw["delta_nm"],
        "dwell_within_sd_nm": dw["within_sd_nm"],
        "dwell_significant": dw["significant"],
        "speed_slow_lag_um": sp["slow_lag_um"],
        "speed_slow_um_s": sp["slow_um_s"],
        "speed_comparable": sp["comparable"],
        "rate_lag_deficit_per_speed_s": rl["deficit_per_speed_s"],
        "rate_lag_relaxation_ms": rl["relaxation_ms"],
        "reference_reversal_um": float(np.mean([abs(a.reversal_nm) for a in rows]) / 1000),
        "drift_at_x0_um": float(abs(drift(rows[-1].t_min)) / 1000),
        "checks_passed": sum(c["pass"] for c in val["checks"]),
        "checks_total": len(val["checks"]),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f0b_noise.json").write_text(json.dumps(nz, indent=1) + "\n",
                                        encoding="utf-8")
    with (OUT / "f0b_ladder.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(halves[0])))
        w.writeheader()
        for h in halves:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                        for k, v in asdict(h).items()})
    (OUT / "f0b_summary.json").write_text(json.dumps(summary, indent=1) + "\n",
                                          encoding="utf-8")
    return {"run": run, "anchors": rows, "drift": drift, "tracks": tracks,
            "noise": nz, "halves": halves, "ladder": lad, "ramps": ramps,
            "grid": grid, "delta_down": down, "delta_up": up,
            "rev_down": rev_down, "rev_up": rev_up,
            "steps": steps, "dwell": dw, "speed": sp, "validity": val,
            "rate_lag": rl,
            "summary": summary}


if __name__ == "__main__":
    d = process()
    s = d["summary"]
    print(f"F0b processed: {s['checks_passed']}/{s['checks_total']} validity checks pass")
    print(f"  sigma(0.5 s) {s['sigma_500ms_nm']:.2f} nm -> d_min {s['d_min_nm']:.1f} nm")
    print(f"  saturation {s['saturation_um']:.3f} um (beta_0 from E5 {s['beta_e5_um']:g}), "
          f"s90 {s['s90_um']:.2f} um")
    print(f"  deadband down {s['deadband_down_um'] * 1000:.0f} nm, "
          f"up {s['deadband_up_um'] * 1000:.0f} nm; half the motion arrives after "
          f"{s['half_arrival_down_um']:.2f} / {s['half_arrival_up_um']:.2f} um")
    print(f"  smallest moving stroke {s['smallest_moving_stroke_um'] * 1000:.0f} nm "
          f"({s['smallest_moving_stroke_usteps']:.0f} usteps)")
    print(f"  dwell {s['dwell_delta_nm']:+.0f} nm (sd {s['dwell_within_sd_nm']:.0f}), "
          f"sliding lag at {s['speed_slow_um_s']:g} um/s "
          f"{s['speed_slow_lag_um']:+.2f} um")
    print(f"  rate lag relaxes in {s['rate_lag_relaxation_ms']:.0f} ms after the "
          f"command stops")
    print(f"  wrote f0b/f0b_noise.json, f0b_ladder.csv, f0b_summary.json")
