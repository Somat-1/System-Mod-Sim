"""Overlay v2 AI1 and v3 X primary-axis f^2-normalized Bode curves."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
V3_DIR = HERE.parents[1]
CHIRP_DIR = V3_DIR.parent
V2_DATA = CHIRP_DIR / "v2" / "analysis" / "plots" / "bode_data.npz"
V3_DATA = V3_DIR / "analysis" / "plots" / "bode_data.npz"
OUT = V3_DIR / "analysis" / "plots" / "bode_X_v2_v3_f2norm_logY_overlay.png"

NOTCH_LOW_HZ = 120.0
NOTCH_HIGH_HZ = 233.4
C_V2 = "#0072B2"
C_V3 = "#D55E00"


def main():
    v2 = np.load(V2_DATA, allow_pickle=False)
    v3 = np.load(V3_DATA, allow_pickle=False)

    freq = np.asarray(v2["freq_grid"], dtype=float)
    v3_freq = np.asarray(v3["freq_grid"], dtype=float)
    if freq.shape != v3_freq.shape or not np.allclose(freq, v3_freq):
        raise ValueError("v2 and v3 frequency grids differ")

    f2 = freq ** 2
    v2_up = np.asarray(v2["up_mag_AI1"], dtype=float) / f2
    v2_down = np.asarray(v2["down_mag_AI1"], dtype=float) / f2
    v3_up = np.asarray(v3["up_mag_AI1"], dtype=float) / f2
    v3_down = np.asarray(v3["down_mag_AI1"], dtype=float) / f2
    v2.close()
    v3.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(freq, v2_up, color=C_V2, lw=1.4,
            label="v2 AI1 (X) -- up sweep")
    ax.plot(freq, v2_down, color=C_V2, lw=1.2, ls="--", alpha=0.82,
            label="v2 AI1 (X) -- down sweep")
    ax.plot(freq, v3_up, color=C_V3, lw=1.4,
            label="v3 X -- up sweep")
    ax.plot(freq, v3_down, color=C_V3, lw=1.2, ls="--", alpha=0.82,
            label="v3 X -- down sweep")
    ax.axvspan(
        NOTCH_LOW_HZ,
        NOTCH_HIGH_HZ,
        color="gray",
        alpha=0.15,
        label="resonance-notch band",
    )

    ax.set_xlim(1.0, 1000.0)
    ax.set_yscale("log")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("acceleration / f$^2$  (m/s$^2$ per Hz$^2$)")
    ax.set_title("v2 vs v3 -- X-axis response, f$^2$-normalized, log y-scale")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    plt.close(fig)
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
