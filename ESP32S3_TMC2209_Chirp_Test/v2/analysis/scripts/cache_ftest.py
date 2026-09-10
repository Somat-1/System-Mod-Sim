"""One-time CSV -> .npy cache for Ftest.csv.

np.loadtxt on a 744 MB / 10.9 M-row text file is far too slow to iterate
against. This parses it in 64 MB binary chunks with numpy's C-accelerated
text-to-float conversion (newlines rewritten to the field delimiter so a
whole chunk becomes one flat parse), then stores the three acceleration
channels as float32 and the time base as float64. Later analysis passes
mmap the .npy and start instantly.
"""
import time
from pathlib import Path
import numpy as np

V2 = Path(__file__).resolve().parents[1]
CSV = V2 / "Ftest.csv"
CACHE = V2 / "analysis" / "cache"
NCOL = 4
CHUNK = 64 << 20


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    size = CSV.stat().st_size
    log(f"parsing {CSV.name} ({size/1e9:.3f} GB)")
    t0 = time.time()

    blocks, rem, nrow, done = [], b"", 0, 0
    with open(CSV, "rb") as fh:
        fh.readline()                                  # header
        done += fh.tell()
        while True:
            raw = fh.read(CHUNK)
            if not raw:
                break
            done += len(raw)
            buf = rem + raw
            cut = buf.rfind(b"\n")                     # keep partial line for next pass
            if cut == -1:
                rem = buf
                continue
            rem = buf[cut + 1:]
            vals = np.fromstring(buf[:cut].replace(b"\n", b","), dtype=np.float64, sep=",")
            if vals.size % NCOL:
                raise RuntimeError(f"chunk not a multiple of {NCOL}: {vals.size}")
            blocks.append(vals.reshape(-1, NCOL))
            nrow += blocks[-1].shape[0]
            log(f"  {done/size*100:5.1f}%  rows={nrow:,}  {time.time()-t0:6.1f}s")

    if rem.strip():
        vals = np.fromstring(rem.strip().replace(b"\n", b","), dtype=np.float64, sep=",")
        blocks.append(vals.reshape(-1, NCOL))
        nrow += blocks[-1].shape[0]

    log(f"stacking {len(blocks)} blocks, {nrow:,} rows ...")
    arr = np.vstack(blocks)
    del blocks

    t = arr[:, 0].copy()
    accel = arr[:, 1:].astype(np.float32)
    del arr

    dt = np.median(np.diff(t[:20000]))
    jitter = np.max(np.abs(np.diff(t) - dt))
    log(f"dt={dt:.9f}s  fs={1/dt:.4f} Hz  max sample-interval deviation={jitter:.3e}s")
    log(f"duration={t[-1]-t[0]:.5f}s  rows={len(t):,}")

    np.save(CACHE / "accel.npy", accel)
    np.save(CACHE / "time.npy", t)
    log(f"saved cache in {time.time()-t0:.1f}s total")


if __name__ == "__main__":
    main()
