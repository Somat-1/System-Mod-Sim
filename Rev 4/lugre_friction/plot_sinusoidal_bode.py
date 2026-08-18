#!/usr/bin/env python3
"""Standalone Bode-style plot (magnitude + phase, linear frequency axis in
Hz -- matching the dense linear 0-2000 Hz grid the sweep now runs on, and
the same axis convention as the linear model's ../rendered_assets/bode_rev4.svg)
from the sinusoidal describing-function sweep. Reads the data saved by
run_sinusoidal_sweep.py rather than re-running the (long) simulation.
Companion to sinusoidal_sweep_montage.svg, which additionally shows two
example time-domain traces.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
DATA_FILE = ASSET_DIR / "sinusoidal_sweep_data.npz"


def main() -> None:
    d = np.load(DATA_FILE)
    freq_hz = d["freq_hz"]
    mags_small, phases_small = d["mags_small"], d["phases_small"]
    mags_large, phases_large = d["mags_large"], d["phases_large"]
    A_small, A_large = float(d["A_small"]), float(d["A_large"])

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)

    ax_mag.plot(freq_hz, mags_small, color="#2b6cb0", linewidth=1.0,
                label=f"A = {A_small:.0e} rad (small, stiction-trapped)")
    ax_mag.plot(freq_hz, mags_large, color="#c05621", linewidth=1.0,
                label=f"A = {A_large:.0e} rad (large, sliding)")
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(r"Rev 4 LuGre sub-branch: empirical sinusoidal sweep, $x_n(t)/\theta_{cmd}(t)$")
    ax_mag.grid(True, linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=9)

    ax_phase.plot(freq_hz, phases_small, color="#2b6cb0", linewidth=1.0)
    ax_phase.plot(freq_hz, phases_large, color="#c05621", linewidth=1.0)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_xlim(0.0, 2000.0)
    ax_phase.grid(True, linewidth=0.4, color="#cccccc")

    fig.tight_layout()
    out_path = ASSET_DIR / "sinusoidal_bode.svg"
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=110)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
