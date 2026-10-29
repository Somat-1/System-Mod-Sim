#!/usr/bin/env python3
"""Load and align the completed E4 stationarity-windows experiment."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ISO_ANALYSIS = HERE.parents[2] / "iso230-2" / "analysis"
if str(ISO_ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ISO_ANALYSIS))
from render_raw_measurement import read_scope  # noqa: E402

DATA = HERE / "data"
SCOPE = DATA / "E4_Experimental_Run.csv"
LOG = DATA / "e4_terminal_20261027_155838.log"
CACHE = DATA / "e4_run_aligned.npz"

WINDOWS = tuple((f"W{i + 1}", c - 3.0, c + 3.0)
                for i, c in enumerate((10., 27., 44., 61., 78., 95., 112.)))
X_REF_MM = 61.35
MOVING_MM_S = 0.05
STOP_RE = re.compile(
    r"# STOP,([A-Z_0-9]+),([A-Z_0-9]+),(\d+),(-?\d+),(\d+),(\d+),([\d.]+),"
    r"(\d+),(\d+),(\d+),(\d+),(\d+),(\d+)")
EVENT_RE = re.compile(
    r"# EVENT,(\d+),([A-Z_0-9]+),(\w+)(?:,label=(\w+))?(?:,target_usteps=(-?\d+))?")


@dataclass
class Stops:
    state: np.ndarray
    window: np.ndarray
    direction: np.ndarray
    index: np.ndarray
    target_usteps: np.ndarray
    target_mm: np.ndarray
    phase_bin: np.ndarray
    move_start_s: np.ndarray
    arrive_s: np.ndarray
    record_start_s: np.ndarray
    record_end_s: np.ndarray
    mscnt: np.ndarray


@dataclass
class Events:
    ms: np.ndarray
    state: np.ndarray
    kind: np.ndarray
    label: np.ndarray


@dataclass
class Run:
    t: np.ndarray
    pos: np.ndarray
    scope_index: np.ndarray
    stops: Stops
    events: Events
    offset_mm: float
    rate_ppm: float
    fit_mad_um: float
    rollover_count: int

    def index_at(self, seconds):
        return np.searchsorted(self.t, seconds)

    def hold(self):
        values = []
        for start, end in zip(self.stops.record_start_s, self.stops.record_end_s):
            i, j = self.index_at([start, end])
            values.append(float(self.pos[i:j].mean()))
        return np.asarray(values)


def parse_log(path: Path = LOG) -> tuple[Stops, Events]:
    stops, events = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = STOP_RE.search(line)
        if match:
            stops.append(match.groups())
            continue
        match = EVENT_RE.search(line)
        if match:
            events.append(match.groups())
    if not stops or not events:
        raise RuntimeError(f"no E4 stops/events found in {path}")
    g = list(zip(*stops))
    s = Stops(
        state=np.asarray(g[0]), window=np.asarray(g[1]),
        direction=np.asarray(g[3], dtype=int), index=np.asarray(g[4], dtype=int),
        target_usteps=np.asarray(g[5], dtype=int), target_mm=np.asarray(g[6], dtype=float),
        phase_bin=np.asarray(g[7], dtype=int),
        move_start_s=np.asarray(g[8], dtype=float) / 1000.0,
        arrive_s=np.asarray(g[9], dtype=float) / 1000.0,
        record_start_s=np.asarray(g[10], dtype=float) / 1000.0,
        record_end_s=np.asarray(g[11], dtype=float) / 1000.0,
        mscnt=np.asarray(g[12], dtype=int),
    )
    h = list(zip(*events))
    e = Events(ms=np.asarray(h[0], dtype=int), state=np.asarray(h[1]),
               kind=np.asarray(h[2]), label=np.asarray([x or "" for x in h[3]]))
    return s, e


def _coarse_shift(t: np.ndarray, disp: np.ndarray, stops: Stops) -> float:
    middle = 0.5 * (stops.record_start_s + stops.record_end_s)
    best = None
    for shift in np.arange(-400.0, 400.0, 0.25):
        indices = np.searchsorted(t, middle + shift)
        if indices.max() >= t.size:
            continue
        residual = stops.target_mm - disp[indices]
        cost = float(np.median(np.abs(residual - np.median(residual))))
        if best is None or cost < best[0]:
            best = (cost, shift)
    if best is None:
        raise RuntimeError("could not seed E4 scope/log alignment")
    return float(best[1])


def _clock_fit(t: np.ndarray, disp: np.ndarray, events: Events,
               shift0: float) -> tuple[float, float]:
    velocity = np.abs(np.diff(disp)) * 1000.0
    moving = np.convolve(velocity, np.ones(5) / 5, mode="same") > MOVING_MM_S
    cumulative = np.r_[0, np.cumsum(moving)]
    starts = events.ms[events.kind == "move_start"].astype(float) / 1000.0
    arrivals = events.ms[events.kind == "arrive"].astype(float) / 1000.0
    ends = arrivals[np.clip(np.searchsorted(arrivals, starts, side="right"),
                            0, arrivals.size - 1)]
    t0, n = t[0], moving.size

    def covered(epsilon: float, shift: float) -> float:
        i = np.clip((((1 + epsilon) * starts + shift - t0) * 1000).astype(int), 0, n)
        j = np.clip((((1 + epsilon) * ends + shift - t0) * 1000).astype(int), 0, n)
        j = np.maximum(j, i)
        return float((cumulative[j] - cumulative[i]).sum() / max((j - i).sum(), 1))

    epsilon_grid = np.arange(-120e-6, 120e-6, 1e-6)
    shift_grid = np.arange(shift0 - 1.5, shift0 + 1.5, 0.002)
    good = np.asarray([[covered(e, s) >= 1.0 - 1e-5 for s in shift_grid]
                       for e in epsilon_grid])
    if not good.any():
        raise RuntimeError("no clock solution puts every logged E4 move inside measured motion")
    ei, si = np.nonzero(good)
    return float(epsilon_grid[ei].mean()), float(shift_grid[si].mean())


def load(use_cache: bool = True) -> Run:
    stops, events = parse_log()
    if use_cache and CACHE.exists():
        z = np.load(CACHE)
        return Run(t=z["t"], pos=z["pos"], scope_index=z["scope_index"],
                   stops=stops, events=events, offset_mm=float(z["offset_mm"]),
                   rate_ppm=float(z["rate_ppm"]), fit_mad_um=float(z["fit_mad_um"]),
                   rollover_count=int(z["rollover_count"]))
    scope = read_scope(SCOPE)
    shift0 = _coarse_shift(scope.time_s, scope.displacement_mm, stops)
    epsilon, shift = _clock_fit(scope.time_s, scope.displacement_mm, events, shift0)
    t = (scope.time_s - shift) / (1 + epsilon)
    middle = 0.5 * (stops.record_start_s + stops.record_end_s)
    indices = np.searchsorted(t, middle)
    residual = stops.target_mm - scope.displacement_mm[indices]
    offset = float(np.median(residual))
    mad = float(np.median(np.abs(residual - offset))) * 1000.0
    pos = scope.displacement_mm + offset
    np.savez_compressed(CACHE, t=t, pos=pos, scope_index=scope.index,
                        offset_mm=offset, rate_ppm=epsilon * 1e6,
                        fit_mad_um=mad, rollover_count=scope.rollover_count)
    return Run(t=t, pos=pos, scope_index=scope.index, stops=stops, events=events,
               offset_mm=offset, rate_ppm=epsilon * 1e6, fit_mad_um=mad,
               rollover_count=scope.rollover_count)


if __name__ == "__main__":
    run = load(use_cache=False)
    print(f"{run.t.size} samples, {run.t[-1] - run.t[0]:.1f} s, "
          f"fit MAD {run.fit_mad_um:.3f} um, offset {run.offset_mm:.6f} mm, "
          f"ESP clock {run.rate_ppm:+.1f} ppm")
    print(f"{run.stops.state.size} stops; position {run.pos.min():.4f} .. "
          f"{run.pos.max():.4f} mm; counter rollovers {run.rollover_count}")
