# ISO 230-2 raw measurement view

The raw scope export is rendered with:

```powershell
python .\render_raw_measurement.py
```

Open the full-resolution interactive Matplotlib viewer with:

```powershell
python .\render_raw_measurement.py --interactive
```

The viewer plots all recorded samples. Use the Matplotlib toolbar to pan or
box-zoom, or place the cursor over either panel and use the scroll wheel to
zoom horizontally. Visible-window vertical autoscaling is enabled by default,
so small settled variations become visible after horizontal zooming.

Keyboard controls:

- `A`: toggle visible-window vertical autoscaling
- `Y`: rescale both vertical axes to the current time window
- `R`: reset to the full recording
- `C`: show/hide the controller command reference
- `G`: show/hide the grid

Data treatment is intentionally minimal: the UINT32 hardware counter is
unwrapped at exact modulo rollovers, the first sample is subtracted, and the
nominal scale of 1 count = 1 nm is used. No filtering, averaging, detrending,
calibration, or eccentricity compensation is applied.
