"""Dewesoft .dxd -> .npy cache for the v3 chirp capture.

Replaces v2's CSV route entirely. v2 exported Ftest.csv (710 MB of text) and
parsed it in chunks; here `dwdatareader` reads the native .dxd through
Dewesoft's own DWDataReaderLib, so no manual export step is needed and the
140 MB binary is read directly.

Channel mapping: the v3 capture names its accelerometer axes X, Y, Z where v2
called them AI 1, AI 2, AI 3. Column order is preserved (X->0, Y->1, Z->2), so
the downstream Bode processing is unchanged apart from the display labels.

Writes analysis/cache/{accel.npy, time.npy} in the same layout
process_chirp_bode_v3.py expects.

Requires: pip install dwdatareader pandas
"""
import time
from pathlib import Path

import numpy as np
import dwdatareader as dw

V3 = Path(__file__).resolve().parents[2]
DXD = V3 / "Ftest2.dxd"
CACHE = V3 / "analysis" / "cache"
AXES = ("X", "Y", "Z")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log(f"opening {DXD.name} ({DXD.stat().st_size/1e6:.1f} MB)")

    with dw.open_file(str(DXD)) as f:
        fs = float(f.info.sample_rate)
        duration = float(f.info.duration)
        log(f"duration={duration:.5f} s  fs={fs:.1f} Hz")

        # DWChannel.scaled() returns a (2, N) array: row 0 is the time base,
        # row 1 the scaled data. Taking it flat silently yields 2N samples.
        cols, n_ref, t = [], None, None
        for name in AXES:
            ch = f[name]
            raw = np.asarray(ch.scaled())
            if raw.ndim != 2 or raw.shape[0] != 2:
                raise RuntimeError(f"{name}: unexpected scaled() shape {raw.shape}")
            ts, data = raw[0], raw[1]
            log(f"  {name}: {data.size:,} samples, unit={ch.unit}, "
                f"range {data.min():.2f}..{data.max():.2f}")
            if n_ref is None:
                n_ref, t = data.size, ts.astype(np.float64)
            elif data.size != n_ref:
                raise RuntimeError(
                    f"channel {name} has {data.size} samples, expected {n_ref}")
            cols.append(data.astype(np.float32))

    accel = np.stack(cols, axis=1)
    dt = float(np.median(np.diff(t[:20000])))
    log(f"dt={dt:.9f} s  fs={1/dt:.4f} Hz (header said {fs:.1f} Hz)")

    log(f"accel shape={accel.shape}  span={t[-1]:.5f} s")
    np.save(CACHE / "accel.npy", accel)
    np.save(CACHE / "time.npy", t)
    log(f"saved cache in {time.time()-t0:.1f}s total")


if __name__ == "__main__":
    main()
