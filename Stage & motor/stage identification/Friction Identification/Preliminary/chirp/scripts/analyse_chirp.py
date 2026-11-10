#!/usr/bin/env python3
"""Check a chirp recording and reduce it to a frequency response.

    python scripts/analyse_chirp.py <folder or chirp_<profile>_windows.csv>

1. Generator check - did the board issue what it was asked to?
     issued/target amplitude   should be 1.0 everywhere
     issued/target frequency   (zero-crossing count; windows of >= 5 cycles)
     worst lag of issued steps behind the target, in microsteps
2. Frequency response, binned (5 Hz to 150 Hz, 25 Hz above), up and down
   sweeps separately: gain = encoder amplitude / issued amplitude
   (1.0 = the rotor follows the command), and phase of the encoder relative
   to the issued steps.

Writes <stem>_frf.csv beside the input.
"""

from __future__ import annotations

import csv
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def bin_of(f):
    return 5 * math.floor(f / 5) if f < 150 else 25 * math.floor(f / 25)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    p = Path(sys.argv[1])
    path = p if p.suffix == ".csv" else next(p.glob("chirp_*_windows.csv"))
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(path.open())]
    sw = [r for r in rows if r["amp_target"] > 0 and r["amp_issued"] > 0]
    if not sw:
        print("no excited windows in", path)
        return 1

    ratio = [r["amp_issued"] / r["amp_target"] for r in sw]
    fr = [r["f_issued_hz"] / r["f_hz"] for r in sw if r["cycles"] >= 5]
    print(f"{path.name}: {len(rows)} windows, {len(sw)} excited, "
          f"{sw[0]['f_hz']:.1f}..{max(r['f_hz'] for r in sw):.1f} Hz")
    print("\nGenerator check (should all be ~1.0)")
    print(f"  issued/target amplitude   min {min(ratio):.3f}  median {st.median(ratio):.3f}  "
          f"max {max(ratio):.3f}")
    if fr:
        print(f"  issued/target frequency   min {min(fr):.3f}  median {st.median(fr):.3f}  "
              f"max {max(fr):.3f}")
    worst = max(sw, key=lambda r: r["lag_max"])
    print(f"  worst issued-step lag     {worst['lag_max']:.0f} microsteps "
          f"(at {worst['f_hz']:.0f} Hz, amplitude {worst['amp_target']:.1f})")
    off = [r for r, q in zip(sw, ratio) if abs(q - 1) > 0.05]
    if off:
        print(f"  {len(off)} windows off by >5%, e.g. at "
              + ", ".join(f"{r['f_hz']:.0f} Hz" for r in off[:6]))

    if max(r["amp_enc"] for r in sw) < 0.5:
        print("\nNo encoder response in this recording (no encoder fitted, or not "
              "counting): the generator check above is the whole result. Take the "
              "response from the external sensor, aligned on the 750 Hz markers.")
        return 0

    # Direction of sweep: frequency rising or falling between windows.
    bins = defaultdict(list)
    for prev, r in zip(sw, sw[1:]):
        d = "up" if r["f_hz"] >= prev["f_hz"] else "down"
        bins[(d, bin_of(r["f_hz"]))].append(r)
    out = []
    for (d, b), rs in sorted(bins.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        gain = st.median(r["amp_enc"] / r["amp_issued"] for r in rs)
        ph = math.degrees(math.atan2(
            sum(math.sin(math.radians(r["phase_deg"])) for r in rs),
            sum(math.cos(math.radians(r["phase_deg"])) for r in rs)))
        out.append(dict(direction=d, f_bin_hz=b, windows=len(rs),
                        amp_issued=round(st.median(r["amp_issued"] for r in rs), 3),
                        amp_enc=round(st.median(r["amp_enc"] for r in rs), 3),
                        gain=round(gain, 4), phase_deg=round(ph, 1)))
    print("\nFrequency response (median per bin; amplitudes in microsteps at MRES 256)")
    print(f"  {'f [Hz]':>7} {'dir':>5} {'issued':>8} {'encoder':>8} {'gain':>7} {'phase':>7}")
    for o in out:
        print(f"  {o['f_bin_hz']:>7} {o['direction']:>5} {o['amp_issued']:8.2f} "
              f"{o['amp_enc']:8.2f} {o['gain']:7.3f} {o['phase_deg']:7.1f}")
    frf = path.with_name(path.stem.replace("_windows", "_frf") + ".csv")
    with frf.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {frf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
