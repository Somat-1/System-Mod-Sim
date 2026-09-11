# v3 — chirp repeat run, 2026-09-10

Second hardware execution of the v2 chirp sequence. Same firmware, same
schedule, same excitation parameters as the v2 campaign — this is a repeat
capture, not a redesign.

## Status

**Recorded and processed.** Bode plots are rendered in `analysis/plots/`.

## The run

Executed 2026-09-10, 18:29:44 → 18:38:50, driven by the committed
`v2/scripts/esp32_v2_chirp_test/esp32_v2_chirp_test.ino` (unmodified `HEAD`
version, not the working-tree edit).

Driver configuration read back off the TMC2209 before motion:

```
GCONF      = 0x000001C4   bit 2  en_spreadCycle = 1   -> SpreadCycle
CHOPCONF   = 0x04030155   bit 28 intpol         = 0   -> MicroPlyer off
                          bits 27:24            = 4   -> MRES 1/16
TPWMTHRS   = 0x00000000                               -> no StealthChop switchover
IHOLD_IRUN = 0x00000404                               -> IHOLD = IRUN
```

Excitation: cruise 3.0 full steps (30 µm peak), notch centred 176.7 Hz,
half-width 56.7 Hz, depth ratio 0.10, 400 mA RMS.

### Measured segment timing

Every segment landed on its design duration with no drift:

| Segment | Start (s) | Duration (s) | Design (s) |
|---|---:|---:|---:|
| PRE_ROLL | 0.00 | 30.00 | 30 |
| MARKER_1 | 30.00 | 0.50 | 0.5 |
| UP_SWEEP | 30.50 | 240.00 | 240 |
| MID_DWELL | 270.50 | 5.00 | 5 |
| MARKER_2 | 275.50 | 0.50 | 0.5 |
| DOWN_SWEEP | 276.00 | 240.00 | 240 |
| TAIL | 516.00 | 30.00 | 30 |
| **Total** | | **546.00** | **546.0** |

Up-sweep started at 1.000000 Hz, down-sweep at 999.999939 Hz, and
`CHIRP_COMPLETE` logged position 0 — the axis returned to its origin, so no
steps were lost across 480 s of continuous excitation.

Because these timings match the design exactly, the analysis constants in
`process_chirp_bode_v3.py` (marker times 30.0 s and 275.5 s, sweep law,
sequence length) carry over from v2 unchanged.

## Layout

```
v3/
├── Ftest2.dxd                  raw Dewesoft capture (146 MB, local-only)
├── analysis/
│   ├── cache/                  .npy cache built from the .dxd (gitignored)
│   ├── plots/                  rendered output
│   ├── scripts/
│   │   ├── cache_dxd_v3.py         .dxd -> .npy cache (reads the binary directly)
│   │   └── process_chirp_bode_v3.py  cache -> Bode plots + bode_data.npz
│   └── chirp_serial_*.csv      ESP32 serial event logs (see below)
```

This mirrors `v2/` so the two revisions are directly comparable, minus
`rendered_assets/` and `scripts/`: those hold v2's pre-hardware design
preview and the schedule generators that produced it. v3 reuses that
design unchanged, so it has nothing of its own to render there.

### Serial event logs

Three logs are kept, because two of them document a hardware failure worth
remembering:

| File | Outcome |
|---|---|
| `chirp_serial_20260910_182203.csv` | **Failed.** Died 65 s into the up-sweep, `rst:0x15 USB_UART_CHIP_RESET` |
| `chirp_serial_20260910_182833.csv` | **Failed.** Reset in the same second `RUN` was sent |
| `chirp_serial_20260910_182936.csv` | **Good run.** Matches `Ftest2.dxd` |

The first failure was caused by a host-side capture script reopening the serial
port mid-run; reopening asserts DTR/RTS, which resets an ESP32-S3 over native
USB. That drove the USB-Serial-JTAG peripheral into a state where every
subsequent `RUN` tripped another reset, including an attempt that opened the
port exactly once. **Re-flashing cleared it**, and the identical firmware then
ran the full 546 s without a hiccup.

Rig note: if this board starts resetting on connect, re-flash it before
suspecting a firmware fault. Never reopen the serial port during a run.

## Rendering the plots

**No CSV export is needed.** v2 went via a 710 MB `Ftest.csv`; v3 reads the
native `.dxd` directly through `dwdatareader`, which wraps Dewesoft's own
`DWDataReaderLib`. The 146 MB binary caches in about two seconds.

```powershell
pip install dwdatareader pandas          # once
python "v3/analysis/scripts/cache_dxd_v3.py"
python "v3/analysis/scripts/process_chirp_bode_v3.py"
```

Output lands in `v3/analysis/plots/`: per-axis Bode magnitude plots (raw linear
and log y, f²-normalised, transmissibility), an all-axis overview,
`recording_overview.png`, and `bode_data.npz` with the curves — the same
deliverable set as v2.

### Capture contents

| | |
|---|---|
| Duration | 574.398 s (546 s sequence + idle margin) |
| Sample rate | 20 000 Hz, verified from the time base |
| Channels | `X`, `Y`, `Z` acceleration in m/s², 11 487 957 samples each |

v2 named its axes `AI 1/2/3`; v3 names them `X/Y/Z`. Column order is unchanged
(X→0, Y→1, Z→2), so only the display labels differ.

### Results of the render

Both sync markers were located and the sequence origin recovered directly,
with no fitting:

| | |
|---|---|
| MARKER_1 | 19.744 s (3 bursts, ~106 ms each) |
| MARKER_2 | 265.237 s |
| Measured separation | 245.4925 s vs 245.500 s design |
| Clock error | **−7.5 ms = −30.5 ppm** DAQ vs ESP32 |
| Sequence t=0 | recording t = −10.255 s |
| Sweep coverage | complete; 38.7 s of idle tail truncated |

Median up-sweep SNR against the idle noise floor: **X 2231×, Y 250×, Z 155×**.
The detector is valid at 999 of 1000 grid points; below 2.00 Hz the sweep does
not dwell for the three cycles the estimator requires.

> ⚠️ These plots use the **committed** excitation parameters — cruise 3.0 full
> steps, notch centred 176.7 Hz, half-width 56.7 Hz, depth 0.10. The
> working-tree edit to `v2/.../esp32_v2_chirp_test.ino` measures the pole at
> 175.0 Hz with Q≈10.7 rather than the assumed 20, which implies this run
> over-attenuates in-band by roughly 2×. Treat v3 as a faithful repeat of the
> v2 configuration, not an improved one.

## Large files

`Ftest2.dxd` (146 MB) follows the same local-only policy as `v2/Ftest.csv`
and is gitignored. No CSV export is produced. **The `.dxd` is a
direct hardware capture and is not regenerable — keep a backup outside git.**
The `.npy` cache is derived data and rebuilds in about two seconds.

## Note on the analysis scripts

`cache_dxd_v3.py` and `process_chirp_bode_v3.py` are adapted from the v2
scripts, following this repository's convention that each revision is
self-contained and resolves paths relative to its own directory. Every analysis constant is identical, because the firmware and
schedule were identical; only the paths, the reader, and the channel labels
differ.

Two changes were made in copying:

1. The v2 `cache_ftest.py` resolves its root with `parents[1]`, which points at
   `v2/analysis/` rather than `v2/`. It must have been run before it was moved
   into `analysis/scripts/`, and is broken in place. The v3 copy uses
   `parents[2]`.
2. The CSV parser was replaced entirely by `cache_dxd_v3.py`, which reads the
   `.dxd` through `dwdatareader`. Note that `DWChannel.scaled()` returns a
   `(2, N)` array — row 0 the time base, row 1 the data — so taking it flat
   silently yields 2N samples.
