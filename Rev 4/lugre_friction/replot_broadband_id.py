#!/usr/bin/env python3
"""One-off: regenerate broadband_id_montage.png from the already-saved
broadband_id_data.npz, using the corrected plotting logic (excitation-band
frequency slicing, wrapped phase) without re-running the ~2.4 hour sweep.
Mirrors run_broadband_id.main()'s plotting block exactly -- if that block
changes again, re-sync this file or delete it."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from broadband_signals import F_HI_HZ, F_LO_HZ

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"
COH_GATE = 0.9


def main() -> None:
    d = np.load(ASSET_DIR / "npz" / "broadband_id_data.npz")
    amplitudes = d["amplitudes"]
    n_amp = len(amplitudes)

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, n_amp - 1)) for i in range(n_amp)]

    fig, (ax_mag, ax_phase, ax_coh) = plt.subplots(3, 1, figsize=(10.0, 11.0), sharex=True)
    for i, A in enumerate(amplitudes):
        color = colors[i]
        f_p, G_p, gamma2_p = d[f"prbs_f_amp{i}"], d[f"prbs_G_amp{i}"], d[f"prbs_gamma2_amp{i}"]
        f_c, mag_c, phase_c = d[f"chirp_f_amp{i}"], d[f"chirp_mag_db_amp{i}"], d[f"chirp_phase_deg_amp{i}"]

        mp = (f_p >= F_LO_HZ) & (f_p <= F_HI_HZ)
        mc = (f_c >= F_LO_HZ) & (f_c <= F_HI_HZ)
        fp, fc_ = f_p[mp], f_c[mc]

        prbs_mag_db = 20.0 * np.log10(np.maximum(np.abs(G_p[mp]), 1e-300))
        prbs_phase_deg = ((np.angle(G_p[mp]) * 180.0 / np.pi + 180.0) % 360.0) - 180.0
        chirp_phase_deg = ((phase_c[mc] + 180.0) % 360.0) - 180.0
        low_coh = gamma2_p[mp] < COH_GATE

        ax_mag.plot(fp, np.ma.masked_where(low_coh, prbs_mag_db), color=color,
                    linewidth=1.0, label=f"A={A*1e3:.3f} mrad RMS (PRBS)")
        ax_mag.plot(fp, np.ma.masked_where(~low_coh, prbs_mag_db), color="#bbbbbb",
                    linewidth=0.8, zorder=0)
        ax_mag.plot(fc_, mag_c[mc], color=color, linewidth=1.0, linestyle="--")

        ax_phase.plot(fp, np.ma.masked_where(low_coh, prbs_phase_deg), color=color,
                      linewidth=0.0, marker=".", markersize=1.5)
        ax_phase.plot(fp, np.ma.masked_where(~low_coh, prbs_phase_deg), color="#bbbbbb",
                      linewidth=0.0, marker=".", markersize=1.2, zorder=0)
        ax_phase.plot(fc_, chirp_phase_deg, color=color, linewidth=0.0,
                      marker=".", markersize=1.5)

        ax_coh.plot(fp, gamma2_p[mp], color=color, linewidth=1.0)

    ax_coh.axhline(COH_GATE, color="#333333", linestyle=":", linewidth=1.0)

    ax_mag.set_xscale("log")
    ax_mag.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("x_n / theta_cmd magnitude -- solid=PRBS (grey where coherence<0.9), dashed=chirp")
    ax_mag.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=7, ncol=2)

    ax_phase.set_xscale("log")
    ax_phase.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_phase.set_ylabel("Phase (deg, wrapped)")
    ax_phase.set_title("x_n / theta_cmd phase -- wrapped to (-180, 180] (dots: PRBS + chirp, "
                        "not a continuous unwrap -- see caption)")
    ax_phase.grid(True, which="both", linewidth=0.4, color="#cccccc")

    ax_coh.set_xscale("log")
    ax_coh.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_coh.set_xlabel("Frequency (Hz)")
    ax_coh.set_ylabel("Coherence gamma^2 (PRBS only)")
    ax_coh.set_title("Coherence -- <1 is nonlinear distortion, not measurement noise (noiseless sim)")
    ax_coh.set_ylim(0.0, 1.05)
    ax_coh.grid(True, which="both", linewidth=0.4, color="#cccccc")

    fig.suptitle("Rev 4 LuGre nonlinear-drive broadband ID: chirp + PRBS, "
                  f"{n_amp} RMS amplitudes")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.97])
    montage_path = ASSET_DIR / "broadband_id_montage.png"
    fig.savefig(montage_path, dpi=110)
    plt.close(fig)
    print(f"Wrote {montage_path}")


if __name__ == "__main__":
    main()
