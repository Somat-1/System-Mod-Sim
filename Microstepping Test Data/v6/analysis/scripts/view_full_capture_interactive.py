"""Interactive matplotlib window over the full EXPrun capture.

A companion to `make_interactive_capture.py`, which builds a Plotly HTML page.
This opens a native matplotlib window instead -- faster to drive, and more
faithful at depth, for one specific reason:

**It re-decimates on every zoom.** The HTML viewer decimates the 1.3M-sample
capture once, down to 100k points, and what you see at the top level is what
you see at every zoom depth -- zooming in never reveals a sample that was
dropped on the way in. Here the decimation is recomputed from the raw arrays
for whatever window is currently on screen: ~4000 buckets wide when zoomed
out, and once the window holds fewer samples than that, the raw samples
themselves, drawn with visible markers. So the far view stays cheap and the
close view is exact, which is the combination the HTML page cannot give.

Decimation is the same envelope-preserving min/max used everywhere else in v6
(`ids_common.decimate_minmax`): each bucket keeps its min and its max, so a
one-sample spike survives instead of being averaged out of existence.

Controls
--------
  toolbar / scroll   pan and box-zoom as usual
  y                  toggle auto-fit of the y-axis to the visible window
                     (on by default -- this is what makes zooming into a flat
                     stretch actually show anything)
  o                  toggle the MRES/tail guide bands
  r                  reset to the full capture
  q                  close

Run from this directory:

    python view_full_capture_interactive.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ids_common import (load_ids, decimate_minmax, still_mask, apply_style,
                        C_MEASURED, C_MUTED, C_GRID, C_INK)
from plot_v6_exp_run import CSV, find_trajectories, MRES_VALUES

# AFTER those imports, deliberately: ids_common calls matplotlib.use("Agg") at
# import time so that every other script in v6 renders headlessly to a file.
# Selecting the GUI backend before them is silently undone by that call, and
# plt.show() then does nothing but warn. force=True switches an already-loaded
# pyplot over. tkinter is the only GUI toolkit present here (no Qt bindings).
_GUI_BACKEND = "TkAgg"
matplotlib.use(_GUI_BACKEND, force=True)

# ~4000 buckets -> ~8000 plotted points. Redrawn on every zoom, so this is a
# per-frame cost, not a one-off: it has to stay small enough to feel instant.
DECIMATE_BUCKETS = 4000
# Below this many samples in view, skip decimation and draw the real samples.
RAW_THRESHOLD = 2 * DECIMATE_BUCKETS

BAND_COLORS = ("#eef2f7", "#f7f2ec")   # alternating, very low contrast
TAIL_COLOR = "#eaf3ee"


class CaptureViewer:
    def __init__(self, t, pos, meta):
        self.t, self.pos, self.meta = t, pos, meta
        self.autoscale_y = True
        self.show_bands = True
        self._updating = False

        apply_style()
        self.fig, self.ax = plt.subplots(figsize=(15, 7.5))
        self.fig.canvas.manager.set_window_title(
            f"EXPrun full capture — {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz")

        (self.line,) = self.ax.plot([], [], color=C_MEASURED, linewidth=0.7,
                                    label="Measured (EL5101)")
        self.ax.set_xlabel("Time from start of capture [s]")
        self.ax.set_ylabel("Position [µm]")
        self.ax.set_axisbelow(True)
        self.ax.spines[["top", "right"]].set_visible(False)
        self.ax.grid(True, color=C_GRID, linewidth=0.6, alpha=0.9)

        self.bands = []
        self._build_bands()

        self.status = self.ax.text(
            0.004, 1.015, "", transform=self.ax.transAxes, fontsize=8.5,
            color=C_MUTED, ha="left", va="bottom")
        self.ax.set_title(
            "Full EXPrun capture — zoom in and the view re-decimates from the "
            "raw samples   ·   keys: y auto-y   o bands   r reset   q quit",
            loc="left", fontsize=10.5, pad=22)

        self.full_xlim = (float(t[0]), float(t[-1]))
        self.ax.set_xlim(*self.full_xlim)
        self.ax.callbacks.connect("xlim_changed", self._on_xlim_changed)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._redraw()
        self.fig.tight_layout()

    # -- guide bands ------------------------------------------------------
    def _build_bands(self):
        """Alternating bands over each MRES group's trajectory legs, plus the
        idle tail, purely as navigation landmarks for a 22-minute capture."""
        spans = find_trajectories(self.t, self.pos, self.meta["dt_s"])
        if len(spans) != 6 * len(MRES_VALUES):
            return
        top = float(np.max(self.pos))
        for gi, mres in enumerate(MRES_VALUES):
            group = spans[gi * 6:(gi + 1) * 6]
            t0, t1 = self.t[group[0][0]], self.t[group[-1][2]]
            self.bands.append(self.ax.axvspan(
                t0, t1, color=BAND_COLORS[gi % 2], zorder=0, linewidth=0))
            self.bands.append(self.ax.text(
                0.5 * (t0 + t1), top, f"MRES {mres}", ha="center", va="top",
                fontsize=9, color=C_INK, alpha=0.65, zorder=1))
        # The post-recording idle tail (what v2 / v2.2 analyse).
        still = still_mask(self.pos, self.meta["dt_s"])
        last_end = spans[-1][2]
        rev = still[last_end:][::-1]
        not_still = np.flatnonzero(~rev)
        if not_still.size:
            tail0 = last_end + (still.size - last_end - not_still[0])
            self.bands.append(self.ax.axvspan(
                self.t[tail0], self.t[-1], color=TAIL_COLOR, zorder=0,
                linewidth=0))
            self.bands.append(self.ax.text(
                0.5 * (self.t[tail0] + self.t[-1]), top, "idle tail",
                ha="center", va="top", fontsize=9, color=C_INK, alpha=0.65,
                zorder=1))

    # -- the adaptive redraw ----------------------------------------------
    def _redraw(self):
        lo_t, hi_t = self.ax.get_xlim()
        i0, i1 = np.searchsorted(self.t, (lo_t, hi_t))
        i0 = max(0, int(i0) - 1)
        i1 = min(self.t.size, int(i1) + 1)
        if i1 - i0 < 2:
            return
        tv, yv = self.t[i0:i1], self.pos[i0:i1]
        n_raw = i1 - i0

        if n_raw <= RAW_THRESHOLD:
            td, yd = tv, yv
            # Individual samples become meaningful at this depth, so show them.
            self.line.set_marker(".")
            self.line.set_markersize(3.5)
            self.line.set_markeredgewidth(0)
            mode = f"raw samples ({n_raw:,})"
        else:
            td, yd = decimate_minmax(tv, yv, target=DECIMATE_BUCKETS)
            self.line.set_marker("None")
            mode = f"{n_raw:,} samples → {td.size:,} min/max points"

        self.line.set_data(td, yd)
        if self.autoscale_y and yv.size:
            lo, hi = float(np.min(yv)), float(np.max(yv))
            pad = max(0.04 * (hi - lo), 1e-4)
            self.ax.set_ylim(lo - pad, hi + pad)
        span = hi_t - lo_t
        self.status.set_text(
            f"window {span:.4g} s   ·   {mode}   ·   "
            f"auto-y {'on' if self.autoscale_y else 'off'}")

    def _on_xlim_changed(self, _ax):
        if self._updating:
            return
        self._updating = True          # set_ylim below re-enters the callbacks
        try:
            self._redraw()
        finally:
            self._updating = False

    def _on_key(self, event):
        if event.key == "y":
            self.autoscale_y = not self.autoscale_y
        elif event.key == "o":
            self.show_bands = not self.show_bands
            for artist in self.bands:
                artist.set_visible(self.show_bands)
        elif event.key == "r":
            self.ax.set_xlim(*self.full_xlim)
        else:
            return
        self._redraw()
        self.fig.canvas.draw_idle()


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    backend = matplotlib.get_backend()
    if backend.lower() == "agg":
        raise SystemExit(
            f"matplotlib fell back to {backend!r}, which cannot open a window. "
            f"Check that tkinter is importable, or set MPLBACKEND={_GUI_BACKEND}.")
    t, pos, meta = load_ids(CSV)
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz "
          f"= {meta['duration_s']:.1f} s ({meta['duration_s'] / 60:.2f} min)")
    print(f"backend: {backend}")
    CaptureViewer(t, pos, meta)
    print("window open — close it (or press q) to return")
    plt.show()


if __name__ == "__main__":
    main()
