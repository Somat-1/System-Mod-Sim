#!/usr/bin/env python3
"""Load an F-series run: the scope capture, the ESP segment log, and the
alignment between them.

Both experiments record the same two things. The TwinCAT YT-Scope exports the
EL5101 counter reading the IDS at 1 kHz, as a tab-separated index/counter
table behind a free-form header; the ESP prints one `# SEG` line per segment
carrying its commanded target and its start and end times in millis() since
boot. Neither carries the other's clock, so the capture is aligned to the
sequence by correlating the measured trace against the commanded one.

    from run_data import load
    run = load(RUNS["F0b"])
    run.t, run.x_nm, run.segments, run.lag_ms, run.r2

Nothing here is specific to F0b or F1 beyond the two entries in RUNS.
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FRICTION = HERE.parent
CACHE = HERE / "cache"

USTEP_NM = 2e6 / (200 * 256)          # 39.0625 nm at 1/256 on a 2 mm lead
DATA_ROW = re.compile(rb"^\d+\t-?\d+")

RUNS = {
    "F0b": {
        "title": "F0b friction scoping",
        "csv": FRICTION / "F0b Friction Scoping/data/F0b_Experimental_Run.csv",
        "log": FRICTION / "F0b Friction Scoping/data/f0b_terminal_20261110_124042.log",
    },
    "F1": {
        "title": "F1 reversal sweep",
        "csv": FRICTION / "F1 Reversal Sweep/data/F1_Experimental_Run.csv",
        "log": FRICTION / "F1 Reversal Sweep/data/f1_terminal_20261110_154639.log",
    },
}


# -- the capture ----------------------------------------------------------
def parse_scope(csv: Path, cache: Path) -> Path:
    """YT-Scope export -> npz of the raw counter. Done once per capture."""
    meta: dict[str, str] = {}
    t0 = time.time()
    with csv.open("rb") as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if not line:
                raise ValueError(f"no sample rows in {csv}")
            if DATA_ROW.match(line):
                fh.seek(pos)
                break
            key, _, value = line.decode("latin-1").strip().partition("\t")
            if key:
                meta.setdefault(key, value.strip())
        body = fh.read().replace(b"EOF", b"")
    flat = np.array(body.split(), dtype=np.int64).reshape(-1, 2)
    index, counter = flat[:, 0], flat[:, 1]
    if not np.array_equal(index, np.arange(index.size)):
        raise ValueError("sample index is not contiguous; the export has gaps")
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache, counter=counter,
        sample_ms=float(meta.get("SampleTime[ms]", 1.0)),
        starttime=meta.get("Starttime of export", ""),
        endtime=meta.get("Endtime of export", ""),
        channel=meta.get("SymbolName", ""))
    print(f"  parsed {csv.name}: {counter.size} samples "
          f"({time.time() - t0:.0f} s) -> {cache.name}")
    return cache


def unwrap(raw: np.ndarray) -> np.ndarray:
    """The EL5101 counter is UINT32 and wraps; return a continuous count."""
    d = np.diff(raw.astype(np.int64))
    d -= (d > 2**31) * 2**32
    d += (d < -2**31) * 2**32
    return np.concatenate([[0], np.cumsum(d)]).astype(np.float64)


# -- the segment log ------------------------------------------------------
@dataclass
class Segment:
    row: int
    block: str
    sub: str
    kind: str          # move, ramp, step, dwell, mark
    target_u: int      # commanded position, microsteps
    rate_mhz: int
    dwell_ms: int
    t0_ms: int         # rebased on the first segment
    t1_ms: int
    mscnt0: int = 0    # microstep counter at the segment start
    mscnt1: int = 0    # and at its end
    status: int = 0    # TMC2209 DRV_STATUS


def load_segments(log: Path) -> list[Segment]:
    out = []
    for line in log.read_text(errors="replace").splitlines():
        if not line.startswith("# SEG,"):
            continue
        f = line[6:].split(",")
        if len(f) == 13:                       # F1 carries group and sub
            row, block, sub, kind, rest = f[0], f[1], f[2], f[3], f[4:]
        elif len(f) == 12:                     # F0b carries the block only
            row, block, sub, kind, rest = f[0], f[1], "", f[2], f[3:]
        else:
            continue
        out.append(Segment(int(row), block, sub, kind, int(rest[0]),
                           int(rest[1]), int(rest[2]), int(rest[6]),
                           int(rest[7]), int(rest[4]), int(rest[5]),
                           int(rest[8], 16)))
    if not out:
        raise ValueError(f"no SEG lines in {log}")
    base = out[0].t0_ms
    for s in out:
        s.t0_ms -= base
        s.t1_ms -= base
    return out


def header_values(log: Path) -> dict[str, str]:
    """The key=value pairs from the SESSION, SEQ_DONE and CENTRED lines."""
    out: dict[str, str] = {}
    for line in log.read_text(errors="replace").splitlines():
        if not line.startswith("#"):
            continue
        tag, _, rest = line[2:].partition(",")
        if tag not in ("SESSION", "SEQ_DONE", "CENTRED", "HOMED", "SUMMARY",
                       "PLAN", "SETTINGS", "BOOT"):
            continue
        for pair in rest.split(","):
            key, _, value = pair.partition("=")
            out.setdefault(f"{tag}.{key.strip()}", value.strip())
    return out


def commanded(segments: list[Segment], start_u: int) -> np.ndarray:
    """Commanded position in nm at 1 kHz, t = 0 at the first segment.

    Linear between segment endpoints: inside a move that is the constant-rate
    part of a trapezoid, which is what all but the first and last 25 ms of
    each move is. Good enough to align on, and exact wherever the stage is
    standing still, which is where everything is measured.
    """
    cuts = np.array([0.0] + [s.t1_ms for s in segments], float)
    vals = np.array([start_u] + [s.target_u for s in segments], float)
    grid = np.arange(0.0, cuts[-1] + 1.0)
    return np.interp(grid, cuts, vals) * USTEP_NM


# -- alignment ------------------------------------------------------------
def best_lag(measured: np.ndarray, command: np.ndarray, velocity: bool = True) -> int:
    """Where the sequence sits inside the capture, in ms.

    On position the correlation is dominated by the 45 mm transits, which are
    45 s long: a second of error costs almost nothing there, and the peak
    lands up to 1.4 s out. Differencing first weights the starts and stops
    instead, which are the sharp features, and the peak then sits within a
    few ms. The difference matters: at the 0.75 mm/s ladder rate, 300 ms of
    misalignment turns into 225 µm of apparent lost motion.
    """
    a, b = (np.diff(measured), np.diff(command)) if velocity else (measured, command)
    a = a - a.mean()
    b = b - b.mean()
    n = 1 << int(np.ceil(np.log2(len(a) + len(b))))
    corr = np.fft.irfft(np.fft.rfft(a, n) * np.conj(np.fft.rfft(b, n)), n)
    return int(np.argmax(corr[:len(a) - len(b) + 1]))


def stillness(measured: np.ndarray, segments: list[Segment], lag: int,
              guard: int = 300) -> tuple[float, int]:
    """Mean 95th-percentile spread inside the dwell windows, in nm.

    The correlation above is dominated by the 45 mm transits, where a second
    of error costs almost nothing, so on its own it lands about a second out.
    Every dwell is a window where the stage is standing still by construction,
    so the lag that makes the dwells quietest is the lag that is right.
    """
    total, n = 0.0, 0
    for s in segments:
        if s.kind not in ("dwell", "mark") or s.t1_ms - s.t0_ms < 1000:
            continue
        a, b = lag + s.t0_ms + guard, lag + s.t1_ms - guard
        if b <= a or b > len(measured):
            continue
        w = measured[a:b]
        total += float(np.percentile(np.abs(w - np.median(w)), 95))
        n += 1
    return (total / n if n else np.inf), n


def excursion_error(measured: np.ndarray, segments: list[Segment], lag: int,
                    min_nm: float = 4e5) -> np.ndarray:
    """Commanded minus measured travel over the ladder's long rungs, in nm.

    The sharpest check that the alignment is right: these rungs run at
    0.75 mm/s, so every millisecond of error is 750 nm of apparent lost
    motion. What is left once the alignment is right is the real thing, and
    it is about a micron.
    """
    reads = []
    for seg, nxt in zip(segments[:-1], segments[1:]):
        if seg.kind != "ramp" or nxt.kind != "dwell":
            continue
        a, b = lag + nxt.t1_ms - 300, lag + nxt.t1_ms
        if b > len(measured):
            continue
        reads.append((seg.target_u, float(np.mean(measured[a:b]))))
    out = []
    for (u0, m0), (u1, m1) in zip(reads[:-1], reads[1:]):
        travel = abs(u1 - u0) * USTEP_NM
        if travel >= min_nm:
            out.append(travel - abs(m1 - m0))
    return np.asarray(out)


@dataclass
class Run:
    name: str
    title: str
    segments: list[Segment]
    header: dict
    x_nm: np.ndarray              # measured, aligned so index = ms of sequence
    cmd_nm: np.ndarray            # commanded, same clock
    lag_ms: int
    counts_per_nm: float
    r2: float
    capture_min: float
    meta: dict = field(default_factory=dict)

    def window(self, seg: Segment, last_ms: int = 0) -> slice:
        """Samples of a segment; last_ms > 0 takes only its final stretch."""
        a, b = seg.t0_ms, min(seg.t1_ms, len(self.x_nm))
        if last_ms:
            a = max(a, b - last_ms)
        return slice(a, b)

    def settled(self, seg: Segment, last_ms: int = 300) -> float:
        """Mean position over the end of a standing-still segment, in nm."""
        w = self.window(seg, last_ms)
        return float(np.mean(self.x_nm[w])) if w.stop > w.start else np.nan


def load(name: str, force: bool = False) -> Run:
    spec = RUNS[name]
    csv, log = Path(spec["csv"]), Path(spec["log"])
    cache = CACHE / f"{name}.npz"
    if force or not cache.exists():
        parse_scope(csv, cache)
    z = np.load(cache, allow_pickle=True)

    # The counter runs against +x: both captures correlate with scale -1.
    measured = -unwrap(z["counter"])
    segments = load_segments(log)
    head = header_values(log)
    start_u = int(head["CENTRED.x_usteps"])
    command = commanded(segments, start_u)

    lag = best_lag(measured, command)
    quiet, windows = stillness(measured, segments, lag, guard=150)
    piece = measured[lag:lag + len(command)]
    a, b = piece - piece.mean(), command - command.mean()
    counts_per_nm = float(np.dot(a, b) / np.dot(b, b))
    r2 = float(1.0 - np.var(a - counts_per_nm * b) / np.var(a))

    err = excursion_error(measured, segments, lag)
    x_nm = (measured[lag:lag + len(command)] - measured[lag]) / counts_per_nm
    x_nm += command[0]
    return Run(name=name, title=spec["title"], segments=segments, header=head,
               x_nm=x_nm, cmd_nm=command, lag_ms=lag,
               counts_per_nm=counts_per_nm, r2=r2,
               capture_min=len(measured) / 60000.0,
               meta={"starttime": str(z["starttime"]),
                     "endtime": str(z["endtime"]),
                     "channel": str(z["channel"]),
                     "samples": int(z["counter"].size),
                     "dwell_spread_nm": quiet,
                     "dwell_windows": windows,
                     "long_rung_shortfall_nm": float(np.mean(err)) if len(err) else np.nan,
                     "long_rung_rms_nm": float(np.std(err)) if len(err) else np.nan,
                     "long_rungs": int(len(err))})


if __name__ == "__main__":
    for key in RUNS:
        run = load(key)
        m = run.meta
        print(f"{key}: lag {run.lag_ms} ms, {run.counts_per_nm:.6f} count/nm, "
              f"R2 {run.r2:.6f}, {len(run.segments)} segments, "
              f"capture {run.capture_min:.2f} min")
        print(f"     dwells quiet to {m['dwell_spread_nm']:.1f} nm over "
              f"{m['dwell_windows']} windows; long ladder rungs short by "
              f"{m['long_rung_shortfall_nm'] / 1000:.3f} +- "
              f"{m['long_rung_rms_nm'] / 1000:.3f} um ({m['long_rungs']} rungs)")
