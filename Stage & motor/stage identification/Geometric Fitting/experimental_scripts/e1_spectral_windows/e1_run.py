#!/usr/bin/env python3
"""Load the E1 run: the 1 kHz interferometer capture and the ESP/TMC log.

The scope export carries no absolute position - only the EL5101 counter - and
its clock is not the ESP's. Both are recovered from the log, which records every
stop with its commanded target and the millisecond window over which the stage
was held. Fitting one time shift and one position offset to those 654 stops puts
the capture into stage coordinates; the fit residual is then the stage's own
position error, which is what every analysis here is about.

The capture starts mid-homing, so it is ~4 min longer than the measurement.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np

_HERE = Path(__file__).resolve().parent
_ISO = _HERE.parents[2] / "iso230-2" / "analysis"
if str(_ISO) not in sys.path:
    sys.path.insert(0, str(_ISO))
from render_raw_measurement import read_scope   # noqa: E402

DATA = _HERE / "data"
SCOPE = DATA / "E1_Experimental_Run.csv"
LOG = DATA / "e1_terminal_20261027_141545.log"
CACHE = DATA / "e1_run_aligned.npz"

WINDOWS = (("W1", 6.0, 18.0), ("W2", 55.0, 67.0), ("W3", 104.0, 116.0))
OVERRUN_MM = 3.0
X_REF_MM = 61.35
MOVING_MM_S = 0.05          # the threshold that calls the stage "in motion"
STOP_RE = re.compile(
    r"# STOP,([A-Z_0-9]+),([A-Z_0-9]+),(\d+),(-?\d+),(\d+),(\d+),([\d.]+),"
    r"(\d+),(\d+),(\d+),(\d+),(\d+)")
EVENT_RE = re.compile(
    r"# EVENT,(\d+),([A-Z_0-9]+),(\w+)(?:,label=(\w+))?(?:,target_usteps=(-?\d+))?")


@dataclass
class Stops:
    state: np.ndarray
    window: np.ndarray
    direction: np.ndarray
    index: np.ndarray
    target_mm: np.ndarray
    move_start_s: np.ndarray
    arrive_s: np.ndarray
    record_start_s: np.ndarray
    record_end_s: np.ndarray


@dataclass
class Events:
    ms: np.ndarray
    state: np.ndarray
    kind: np.ndarray
    label: np.ndarray

    def when(self, state, kind, label=None):
        m = (self.state == state) & (self.kind == kind)
        if label is not None:
            m &= self.label == label
        return self.ms[m] / 1000.0


@dataclass
class Run:
    t: np.ndarray            # scope seconds, shifted onto the ESP clock
    pos: np.ndarray          # stage position, mm
    stops: Stops
    events: Events
    offset_mm: float
    rate_ppm: float
    fit_mad_um: float

    def index_at(self, seconds):
        return np.searchsorted(self.t, seconds)

    def hold(self, half_ms=250):
        """Mean position over each stop's record window."""
        k = self.index_at(0.5 * (self.stops.record_start_s +
                                 self.stops.record_end_s))
        return np.array([self.pos[i - half_ms:i + half_ms].mean() for i in k])


def _parse_log(path: Path):
    stops, events = [], []
    for line in path.read_text(errors="replace").splitlines():
        m = STOP_RE.search(line)
        if m:
            stops.append(m.groups())
            continue
        m = EVENT_RE.search(line)
        if m:
            events.append(m.groups())
    g = list(zip(*stops))
    s = Stops(state=np.array(g[0]), window=np.array(g[1]),
              direction=np.array(g[3], dtype=int),
              index=np.array(g[4], dtype=int),
              target_mm=np.array(g[6], dtype=float),
              move_start_s=np.array(g[7], dtype=float) / 1000.0,
              arrive_s=np.array(g[8], dtype=float) / 1000.0,
              record_start_s=np.array(g[9], dtype=float) / 1000.0,
              record_end_s=np.array(g[10], dtype=float) / 1000.0)
    h = list(zip(*events))
    e = Events(ms=np.array(h[0], dtype=int), state=np.array(h[1]),
               kind=np.array(h[2]),
               label=np.array([x or "" for x in h[3]]))
    return s, e


def _coarse_shift(t, disp, stops):
    """Where the capture sits on the ESP clock, from the commanded holds.

    Position alone fixes this only to about a second - the stage is standing
    still for 2.5 s around every stop - but that is close enough to seed the
    motion fit below.
    """
    mid = 0.5 * (stops.record_start_s + stops.record_end_s)
    best = None
    for shift in np.arange(-400.0, 400.0, 0.25):
        k = np.searchsorted(t, mid + shift)
        if k.max() >= t.size:
            continue
        r = stops.target_mm - disp[k]
        cost = np.median(np.abs(r - np.median(r)))
        if best is None or cost < best[0]:
            best = (cost, shift)
    return best[1]


def _clock_fit(t, disp, events, shift0):
    """Rate and offset of the ESP clock against the scope's, from the moves.

    The ESP's millisecond counter is a free-running crystal and drifts tens of
    ppm against the TwinCAT 1 kHz clock - over half an hour that is a tenth of
    a second, which is most of a 0.12 mm move. Holds cannot see it; motion can.
    Every logged move interval is required to fall inside measured motion, and
    the solution is the centre of the region where all 676 of them do.
    """
    v = np.abs(np.diff(disp)) * 1000.0
    moving = np.convolve(v, np.ones(5) / 5, mode="same") > MOVING_MM_S
    cumulative = np.r_[0, np.cumsum(moving)]
    starts = events.ms[events.kind == "move_start"].astype(float) / 1000.0
    arrivals = events.ms[events.kind == "arrive"].astype(float) / 1000.0
    ends = arrivals[np.clip(np.searchsorted(arrivals, starts, side="right"),
                            0, arrivals.size - 1)]
    t0, n = t[0], moving.size

    def covered(eps, shift):
        i = np.clip((((1 + eps) * starts + shift - t0) * 1000).astype(int), 0, n)
        j = np.clip((((1 + eps) * ends + shift - t0) * 1000).astype(int), 0, n)
        j = np.maximum(j, i)
        return (cumulative[j] - cumulative[i]).sum() / max((j - i).sum(), 1)

    eps_grid = np.arange(-120e-6, 120e-6, 1e-6)
    shift_grid = np.arange(shift0 - 1.5, shift0 + 1.5, 0.002)
    good = np.array([[covered(e, s) >= 1.0 - 1e-5 for s in shift_grid]
                     for e in eps_grid])
    if not good.any():
        raise RuntimeError("no (rate, shift) puts every logged move inside motion")
    ei, si = np.nonzero(good)
    return float(eps_grid[ei].mean()), float(shift_grid[si].mean())


def load(use_cache=True) -> Run:
    stops, events = _parse_log(LOG)
    if use_cache and CACHE.exists():
        z = np.load(CACHE)
        return Run(t=z["t"], pos=z["pos"], stops=stops, events=events,
                   offset_mm=float(z["offset_mm"]), rate_ppm=float(z["rate_ppm"]),
                   fit_mad_um=float(z["fit_mad_um"]))
    scope = read_scope(SCOPE)
    shift0 = _coarse_shift(scope.time_s, scope.displacement_mm, stops)
    eps, shift = _clock_fit(scope.time_s, scope.displacement_mm, events, shift0)
    t = (scope.time_s - shift) / (1 + eps)       # now the ESP's own seconds
    mid = 0.5 * (stops.record_start_s + stops.record_end_s)
    k = np.searchsorted(t, mid)
    r = stops.target_mm - scope.displacement_mm[k]
    offset = float(np.median(r))
    mad = float(np.median(np.abs(r - offset))) * 1000.0
    pos = scope.displacement_mm + offset
    np.savez_compressed(CACHE, t=t, pos=pos, offset_mm=offset,
                        rate_ppm=eps * 1e6, fit_mad_um=mad)
    return Run(t=t, pos=pos, stops=stops, events=events, offset_mm=offset,
               rate_ppm=eps * 1e6, fit_mad_um=mad)


if __name__ == "__main__":
    run = load(use_cache=False)
    print(f"{run.t.size} samples, {run.t[-1] - run.t[0]:.1f} s, "
          f"fit MAD {run.fit_mad_um:.3f} um, offset {run.offset_mm:.6f} mm, "
          f"ESP clock {run.rate_ppm:+.1f} ppm")
    print(f"position {run.pos.min():.4f} .. {run.pos.max():.4f} mm")
