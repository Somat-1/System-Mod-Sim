"""Shared loading and styling for the v4 EL5101 IDS captures.

Both exports are Beckhoff EL5101-0011 encoder traces written by a TWinCAT "YT
Scope Project": a tab-separated header block, then `index<TAB>counter` rows at
the rate given by `SampleTime[ms]`. The counter is a raw UINT32 count; the
repository's established scale is 1 nm per count (see
`v2/temp_step_size_determination/analyze_step_size_determination.py`), so a
count is a nanometre and 1 full step of the stage is 10 um.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V4_DIR = Path(__file__).resolve().parents[2]
RUNS_DIR = V4_DIR / "data" / "hardware_runs"
ANALYSIS_DIR = V4_DIR / "analysis"

NM_PER_COUNT = 1.0
FULL_STEP_UM = 10.0          # 2 mm lead / 200 full steps per revolution.

# Validated pair (OKLab dE: 32 normal vision, 24 worst-case CVD, contrast
# >= 3:1 on both surfaces). Identity never rests on colour alone here - the
# commanded trace is also dashed and both series are named in the legend.
C_MEASURED = "#1769aa"
C_COMMANDED = "#d94801"
C_INK = "#1a1a1a"
C_MUTED = "#6b7280"
C_GRID = "#d7dbe0"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": C_GRID,
        "axes.labelcolor": C_INK,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.labelsize": 9,
        "axes.grid": True,
        "grid.color": C_GRID,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.7,
        "xtick.color": C_MUTED,
        "ytick.color": C_MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "font.size": 9,
        "text.color": C_INK,
        "figure.dpi": 150,
    })


def load_ids(path):
    """Return (t_seconds, position_um, meta) from a YT Scope export."""
    counts, dt_ms, meta = [], None, {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and not parts[0].isdigit():
                if parts[0] in ("SampleTime[ms]", "Name", "Unit", "Data-Type",
                                "Starttime of export", "Endtime of export"):
                    meta.setdefault(parts[0], parts[1])
                if parts[0] == "SampleTime[ms]":
                    dt_ms = float(parts[1])
            elif len(parts) == 2 and parts[0].isdigit() and parts[1].strip().lstrip("-").isdigit():
                counts.append(int(parts[1]))
    if dt_ms is None:
        raise RuntimeError(f"{path.name}: no SampleTime[ms] in header")
    c = np.asarray(counts, dtype=np.int64)
    if c.size == 0:
        raise RuntimeError(f"{path.name}: no data rows parsed")
    t = np.arange(c.size) * (dt_ms / 1000.0)
    pos_um = (c - c[0]) * NM_PER_COUNT / 1000.0
    meta["dt_s"] = dt_ms / 1000.0
    meta["fs_hz"] = 1000.0 / dt_ms
    meta["n"] = int(c.size)
    meta["duration_s"] = float(t[-1])
    return t, pos_um, meta


def decimate_minmax(t, y, target=4000):
    """Envelope-preserving decimation: keeps min and max per bucket, so a
    2.5 M-sample trace renders honestly instead of dropping its extremes."""
    n = t.size
    if n <= target * 2:
        return t, y
    bucket = int(np.ceil(n / target))
    usable = (n // bucket) * bucket
    tb = t[:usable].reshape(-1, bucket)
    yb = y[:usable].reshape(-1, bucket)
    lo_i = yb.argmin(axis=1)
    hi_i = yb.argmax(axis=1)
    rows = np.arange(tb.shape[0])
    pts_t = np.stack([tb[rows, np.minimum(lo_i, hi_i)],
                      tb[rows, np.maximum(lo_i, hi_i)]], axis=1).ravel()
    pts_y = np.stack([yb[rows, np.minimum(lo_i, hi_i)],
                      yb[rows, np.maximum(lo_i, hi_i)]], axis=1).ravel()
    if usable < n:
        pts_t = np.append(pts_t, t[usable:])
        pts_y = np.append(pts_y, y[usable:])
    return pts_t, pts_y


def find_marker_onsets(t, pos_um, threshold_um=50.0, min_gap_s=5.0):
    """Marker moves are 20 full steps (200 um) negative-then-return, orders of
    magnitude larger than any microstep, so a plain threshold finds them."""
    below = pos_um < (np.median(pos_um) - threshold_um)
    onsets = np.flatnonzero(np.diff(below.astype(np.int8)) == 1) + 1
    if onsets.size == 0:
        return np.array([])
    keep = [onsets[0]]
    for i in onsets[1:]:
        if t[i] - t[keep[-1]] >= min_gap_s:
            keep.append(i)
    return np.asarray(keep)


def plot_full_capture(t, pos_um, meta, out_path, title, subtitle=None,
                      marks=None):
    """The whole raw capture, always rendered, at full extent."""
    apply_style()
    fig, ax = plt.subplots(figsize=(13, 4.2))
    td, yd = decimate_minmax(t, pos_um)
    ax.plot(td, yd, color=C_MEASURED, lw=0.7, label="Measured (EL5101)")
    if marks is not None and len(marks):
        for i, m in enumerate(marks):
            ax.axvline(t[m], color=C_COMMANDED, lw=0.9, ls=":", alpha=0.8,
                       label="Marker onset" if i == 0 else None)
    ax.set_xlabel("Time from start of capture [s]")
    ax.set_ylabel("Position [µm]")
    ax.set_xlim(t[0], t[-1])
    ax.set_title(title, loc="left", pad=24 if subtitle else 8)
    if subtitle:
        ax.text(0.0, 1.012, subtitle, transform=ax.transAxes,
                fontsize=8.5, color=C_MUTED, ha="left", va="bottom")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
