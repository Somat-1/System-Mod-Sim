"""Interactive, zoomable HTML view of the full EXPrun capture.

The static `overview_full_capture.png` is fixed at print resolution -- fine
for a report, useless for digging into one specific leg or marker by hand.
This renders the same trace as a self-contained Plotly HTML page (WebGL,
pan/zoom/box-zoom, a range slider, hover readout) that opens directly in a
browser with no server and no internet connection needed.

The raw capture is 1.3M samples. Handed to Plotly unchanged, that is slow to
pan and zoom even in WebGL, and every sample past the first million or so on
some browsers silently fails to render at all. It is decimated with the same
envelope-preserving min/max approach as the static plots (`decimate_minmax` in
ids_common.py): each time bucket keeps only its min and max sample, so spikes
and dips survive instead of being averaged away, and re-zooming the page does
not re-decimate -- what you see at the top level is what is actually in the
file at every zoom depth.
"""
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ids_common import (ANALYSIS_DIR, load_ids, decimate_minmax, still_mask,
                        C_MEASURED, C_COMMANDED)
from plot_v6_exp_run import CSV, find_trajectories, leg_rest_levels, MRES_VALUES

OUT_DIR = ANALYSIS_DIR / "exp_run"
C_ORIGIN = "#178a5f"   # distinct from measured/commanded; matches the "true
                       # rest" colour in rest_window_diagnostic.png.

# ~200k points after min/max pairing (2 per bucket) -- smooth pan/zoom on
# ordinary hardware in WebGL, versus 1.3M raw. Far more than the "half" first
# guessed, because that guess was made without knowing the true sample count;
# this is chosen for browser performance instead, not as a fraction of it.
DECIMATE_TARGET = 100_000


def main():
    t, pos, meta = load_ids(CSV)
    n_raw = meta["n"]
    print(f"loaded {CSV.name}: {n_raw:,} samples @ {meta['fs_hz']:.0f} Hz "
          f"= {meta['duration_s']:.1f} s ({meta['duration_s']/60:.2f} min)")

    td, yd = decimate_minmax(t, pos, target=DECIMATE_TARGET)
    print(f"decimated to {td.size:,} points ({td.size / n_raw:.1%} of raw) "
          f"via min/max envelope, {DECIMATE_TARGET:,} buckets")

    spans = find_trajectories(t, pos, meta["dt_s"])

    # Origin trace: the actual resting position immediately before and after
    # each leg's ramp (see leg_rest_levels in plot_v6_exp_run.py -- bounded to
    # that leg's own settled dwell, not the gap to the previous leg, which can
    # otherwise land on an unrelated marker or oscillation dwell; see
    # plot_v6_rest_window_diagnostic.py). At full zoom this sits almost
    # exactly on the raw trace's own resting plateaus and is easy to miss;
    # zoom in vertically on any dwell to see the gap between this line and
    # where the axis started (the dashed reference) open up over the
    # campaign -- that gap is the cumulative drift.
    still = still_mask(pos, meta["dt_s"])
    idx_b, rest_b, idx_a, rest_a = leg_rest_levels(t, pos, spans, still, meta["dt_s"])
    order = np.argsort(np.concatenate([idx_b, idx_a]))
    origin_t = t[np.concatenate([idx_b, idx_a])[order]]
    origin_y = np.concatenate([rest_b, rest_a])[order]

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=td, y=yd, mode="lines",
        line=dict(color=C_MEASURED, width=1),
        name="Measured (EL5101)",
        hovertemplate="t=%{x:.3f} s<br>pos=%{y:.2f} µm<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[t[0], t[-1]], y=[origin_y[0], origin_y[0]], mode="lines",
        line=dict(color=C_ORIGIN, width=1, dash="dash"),
        name="Starting origin (leg 1's rest, held constant)",
        hovertemplate="reference: %{y:.2f} µm<extra></extra>",
    ))
    fig.add_trace(go.Scattergl(
        x=origin_t, y=origin_y, mode="lines+markers",
        line=dict(color=C_ORIGIN, width=1.5),
        marker=dict(size=4),
        name="Measured origin (rest before/after each leg)",
        hovertemplate="t=%{x:.3f} s<br>rest=%{y:.2f} µm<extra></extra>",
    ))

    for i, (s, _apex, e) in enumerate(spans):
        fig.add_vrect(
            x0=t[s], x1=t[e], fillcolor=C_COMMANDED, opacity=0.10, line_width=0,
            annotation_text=f"#{i+1}" if i % 2 == 0 else None,
            annotation_position="top left", annotation_font_size=9,
            annotation_font_color=C_COMMANDED,
        )
    for i in range(0, len(spans), 6):
        mres = MRES_VALUES[i // 6] if i // 6 < len(MRES_VALUES) else None
        if mres is not None:
            fig.add_vline(x=t[spans[i][0]], line_width=1, line_dash="dot",
                          line_color="#6b7280")

    fig.update_layout(
        title=dict(text=(f"EXPrun — full v6 campaign, interactive"
                         f"<br><sub>{n_raw:,} raw samples decimated to "
                         f"{td.size:,} (min/max envelope, {meta['fs_hz']:.0f} Hz "
                         f"native)  ·  1 nm/count  ·  shaded bands = detected "
                         f"10 mm trajectory legs  ·  green = measured resting "
                         f"origin per leg vs. its starting value  ·  drag on "
                         f"the chart to box-zoom both axes, scroll to zoom, "
                         f"double-click to reset, drag the strip below to "
                         f"pan</sub>"),
                  x=0.01, xanchor="left"),
        # dragmode="zoom" is the default, but set it explicitly: a rectangle
        # dragged on the main plot zooms x AND y to the box's extent. The
        # range slider strip below the axis is x-only by design (it is a
        # pan/overview control, not a substitute for the box-zoom above).
        dragmode="zoom",
        xaxis=dict(title="Time from start of capture [s]",
                  rangeslider=dict(visible=True), showspikes=True),
        yaxis=dict(title="Position [µm]", showspikes=True, fixedrange=False,
                  autorange=True),
        template="plotly_white",
        hovermode="x unified",
        margin=dict(t=100),
    )

    out_path = OUT_DIR / "overview_full_capture_interactive.html"
    # scrollZoom: mouse-wheel zooms both axes around the cursor, in addition
    # to the drag-box zoom above. doubleClick "reset" restores the full view
    # on both axes (default is "reset+autosize", which is equivalent here).
    fig.write_html(out_path, include_plotlyjs=True, full_html=True,
                   config=dict(scrollZoom=True, doubleClick="reset",
                              displaylogo=False))
    size_mb = out_path.stat().st_size / 1e6
    print(f"wrote {out_path.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
