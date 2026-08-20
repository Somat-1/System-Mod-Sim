#!/usr/bin/env python3
"""x_n RMS amplitude vs. frequency, dB re 1 m, one curve per RMS command
amplitude -- pure arithmetic on the already-saved broadband_id_data.npz
(no re-simulation). Drive-integrity check per the 2026-08-20 request: plot
|G(f)|*A directly (not the per-unit-command gain |G(f)| alone) against a
dashed ideal-no-slip-tracking reference (L/2*pi)*A per amplitude. The
vertical gap to that dashed line is the tracking error; if curves cross or
invert relative to EACH OTHER, that's the excitation failing to drive the
plant as commanded (motor pull-out), not a friction effect -- friction can
flatten/soften the curves, it cannot invert them.

Also prints the DERIVED torque-limited frequency cap for each amplitude
(not a flat 100 Hz -- see module docstring in run_broadband_id.py history /
conversation): the peak angular acceleration of a sinusoidal command
A*sin(wt) is A*w^2, so I_eff*A*w^2 < T_hold gives

    f_max(A) = (1/2pi) * sqrt(T_hold / (I_eff * A)),   I_eff = I_m+I_c+I_s+I_sb
               + (M_screw+M_s)*(L/2pi)^2

marked on the plot per amplitude -- pull-out should start there, not at a
fixed frequency. Also prints the amplitude-independent breakaway lag angle
Delta_theta = arcsin(T_fric_total/T_hold) (T_fric_total = Fs_way and Fs_nut
reflected through L/2pi, plus Ts_sb): if this is already a large fraction
of the tested amplitudes, no frequency cap fixes that amplitude regime.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from broadband_signals import F_HI_HZ, F_LO_HZ
from lugre_model import load_parameters

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"


def torque_limit_f_max(p: dict[str, float], A: np.ndarray) -> np.ndarray:
    lead_ratio = p["L"] / (2.0 * np.pi)
    I_eff = p["I_m"] + p["I_c"] + p["I_s"] + p["I_sb"] + (p["M_screw"] + p["M_s"]) * lead_ratio ** 2
    return (1.0 / (2.0 * np.pi)) * np.sqrt(p["T_hold"] / (I_eff * A))


def breakaway_lag_angle(p: dict[str, float]) -> float:
    lead_ratio = p["L"] / (2.0 * np.pi)
    T_fric_total = p["Fs_way"] * lead_ratio + p["Fs_nut"] * lead_ratio + p["Ts_sb"]
    return float(np.arcsin(T_fric_total / p["T_hold"])), T_fric_total


def main() -> None:
    p = load_parameters()
    lead_ratio = p["L"] / (2.0 * np.pi)

    dtheta_breakaway, T_fric_total = breakaway_lag_angle(p)
    print(f"Breakaway lag angle Delta_theta = arcsin(T_fric_total/T_hold) = "
          f"{dtheta_breakaway*1e3:.3f} mrad  (T_fric_total={T_fric_total*1e3:.4f} mN*m, "
          f"T_hold={p['T_hold']*1e3:.1f} mN*m)")

    d = np.load(ASSET_DIR / "npz" / "broadband_id_data.npz")
    amplitudes = d["amplitudes"]
    n_amp = len(amplitudes)
    f_max = torque_limit_f_max(p, amplitudes)

    print("\nTorque-limited f_max(A) = (1/2pi)*sqrt(T_hold/(I_eff*A)) per amplitude:")
    for A, fm in zip(amplitudes, f_max):
        print(f"  A={A*1e3:9.4f} mrad RMS -> f_max={fm:8.2f} Hz  "
              f"(2000 Hz sweep exceeds it by {2000.0/fm:.1f}x)")

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, n_amp - 1)) for i in range(n_amp)]

    fig, ax = plt.subplots(figsize=(10.5, 7.0))

    f_prbs = d["prbs_f_amp0"]
    mask_p = (f_prbs >= F_LO_HZ) & (f_prbs <= F_HI_HZ)
    fp = f_prbs[mask_p]
    f_chirp = d["chirp_f_amp0"]
    mask_c = (f_chirp >= F_LO_HZ) & (f_chirp <= F_HI_HZ)
    fc = f_chirp[mask_c]

    for i, A in enumerate(amplitudes):
        color = colors[i]
        G_p = d[f"prbs_G_amp{i}"][mask_p]
        xn_db_p = 20.0 * np.log10(np.abs(G_p)) + 20.0 * np.log10(A)
        ax.plot(fp, xn_db_p, color=color, linewidth=1.0,
                label=f"A={A*1e3:.3f} mrad RMS (PRBS)")

        mag_db_c = d[f"chirp_mag_db_amp{i}"][mask_c]
        xn_db_c = mag_db_c + 20.0 * np.log10(A)
        ax.plot(fc, xn_db_c, color=color, linewidth=1.0, linestyle="--")

        ideal_db = 20.0 * np.log10(lead_ratio * A)
        ax.axhline(ideal_db, color=color, linestyle=":", linewidth=0.8, alpha=0.6)

        ax.axvline(f_max[i], color=color, linestyle="-", linewidth=0.6, alpha=0.4)

    ax.set_xscale("log")
    ax.set_xlim(F_LO_HZ, F_HI_HZ)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("x_n RMS amplitude (dB re 1 m) = 20*log10(|G|*A)")
    ax.set_title("Driven stage amplitude vs. frequency -- solid=PRBS, dashed=chirp, "
                 "dotted=ideal no-slip tracking\nthin vertical lines: torque-limited f_max(A) per amplitude")
    ax.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=7, ncol=2)

    fig.tight_layout()
    out_path = ASSET_DIR / "broadband_driven_amplitude.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
