#!/usr/bin/env python3
"""Modal superposition of the frictionless baseline's command-to-stage
transfer function x_n(s)/theta_cmd(s), requested 2026-08-20; revised
2026-08-20 to run 10 Hz-10 kHz (nothing below 10 Hz, and this now shows
modes 5/6 which the old 0.1-2000 Hz range cut off), drop the plotted
analytic overlay (kept internally only to drive the mismatch shading),
drop antiresonance markers, and compress the phase panel.

x_n(s)/theta_cmd(s) = C_y (M s^2 + C s + K)^-1 B_u[:,0] is expanded as a
sum of six independent second-order (mass-normalized-mode) systems:

    G(s) ~= sum_j  R_j / (s^2 + 2*zeta_j*omega_j*s + omega_j^2)
    R_j = (C_y @ phi_j) * (phi_j.T @ B_u[:,0])

phi from scipy.linalg.eigh(K, M) is mass-normalized (phi.T @ M @ phi = I)
by construction. zeta_j = diag(phi.T @ C @ phi)_j / (2*omega_j) is the
PROPORTIONAL-damping approximation -- it only uses the diagonal of
phi.T @ C @ phi, discarding any off-diagonal (non-proportional/modal-
coupling) terms. The analytic G(jw) (direct (sI-A)^-1 solve, no
proportional-damping assumption, same laplace_transfer_function as
build_bode_rev4.py, evaluated on this script's own extended grid rather
than the saved 0-2000 Hz bode_rev4_data.npz) is computed to check that
approximation -- shaded red wherever the two disagree by more than
MAG_SHADE_THRESHOLD_DB -- but is not itself drawn, per review. The shaded
region (2026-08-20 revision) marks everywhere beyond the first frequency
where the summed magnitude drops below MAG_SHADE_THRESHOLD_DB, not just the
narrow bands where the 6-mode sum disagrees with the analytic curve --
still computed and printed (max/median diff) as a numeric check, just no
longer the shading criterion.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from build_bode_rev4 import (
    STATE_LABELS,
    build_matrices,
    build_state_space,
    laplace_transfer_function,
    load_parameters,
)

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "rendered_assets" / "temp"
F_LO_HZ = 10.0
F_HI_HZ = 10000.0
N_FREQ = 4000
MAG_SHADE_THRESHOLD_DB = -110.0   # shade everything beyond the first crossing below this
N_DIMMED = 3   # modes beyond this index (0-based) are drawn dimmed/dashed


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)
    n = M.shape[0]

    C_y = np.zeros((1, n))
    C_y[0, STATE_LABELS.index("x_n")] = 1.0
    b_col = B_u[:, 0]   # theta_cmd column

    lam, phi = eigh(K, M)   # phi.T @ M @ phi = I (mass-normalized) by construction
    omega = np.sqrt(lam)
    freq_hz = omega / (2.0 * np.pi)

    zeta = np.diag(phi.T @ C @ phi) / (2.0 * omega)

    R = (C_y @ phi).flatten() * (phi.T @ b_col)
    R_norm = np.abs(R) / np.abs(R[0])

    print(f"{'mode':>4s} {'f (Hz)':>10s} {'zeta':>10s} {'|R|':>12s} {'|R|/|R1|':>10s}")
    for j in range(n):
        print(f"{j+1:4d} {freq_hz[j]:10.2f} {zeta[j]:10.5f} {abs(R[j]):12.4e} {R_norm[j]:10.4f}")

    frequencies_hz = np.logspace(np.log10(F_LO_HZ), np.log10(F_HI_HZ), N_FREQ)
    s_grid = 1j * 2.0 * np.pi * frequencies_hz

    mode_terms = np.array([
        R[j] / (s_grid ** 2 + 2.0 * zeta[j] * omega[j] * s_grid + omega[j] ** 2)
        for j in range(n)
    ])
    G_modal = mode_terms.sum(axis=0)
    mag_modal_db = 20.0 * np.log10(np.maximum(np.abs(G_modal), 1e-300))
    phase_modal_deg = np.unwrap(np.angle(G_modal)) * 180.0 / np.pi

    # Analytic G(jw): same machinery as build_bode_rev4.py, own grid -- kept
    # only to drive the mismatch shading below, not plotted.
    A_fo, B_fo, C_y_fo = build_state_space(M, C, K, B_u)
    b_col_fo = B_fo[:, 0]
    G_analytic = np.array([laplace_transfer_function(A_fo, b_col_fo, C_y_fo, s) for s in s_grid])
    mag_analytic_db = 20.0 * np.log10(np.maximum(np.abs(G_analytic), 1e-15))
    mag_diff = np.abs(mag_modal_db - mag_analytic_db)
    print(f"\nModal-sum vs. analytic G(jw): max |magnitude diff| = {mag_diff.max():.4f} dB, "
          f"median = {np.median(mag_diff):.4f} dB")

    peak_db_main = max(np.max(20.0 * np.log10(np.maximum(np.abs(mode_terms[j]), 1e-300)))
                        for j in range(N_DIMMED))
    peak_db_dimmed = max(np.max(20.0 * np.log10(np.maximum(np.abs(mode_terms[j]), 1e-300)))
                          for j in range(N_DIMMED, n))
    down_db = peak_db_main - peak_db_dimmed
    print(f"Modes {N_DIMMED+1}-{n} peak {down_db:.1f} dB below the tallest of modes 1-{N_DIMMED}.")

    fig, (ax_mag, ax_phase) = plt.subplots(
        2, 1, figsize=(10.5, 8.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]},
    )
    cmap = plt.get_cmap("tab10")
    for j in range(n):
        mode_mag_db = 20.0 * np.log10(np.maximum(np.abs(mode_terms[j]), 1e-300))
        dimmed = j >= N_DIMMED
        ax_mag.plot(frequencies_hz, mode_mag_db,
                    color="#bbbbbb" if dimmed else cmap(j),
                    linewidth=0.7 if dimmed else 0.9,
                    linestyle="--" if dimmed else "-",
                    alpha=0.7,
                    label=f"mode {j+1} ({freq_hz[j]:.1f} Hz, ζ={zeta[j]:.4f}, "
                          f"|R|/|R1|={R_norm[j]:.3f})")

    # Shade everything beyond the first point where the summed magnitude
    # drops below MAG_SHADE_THRESHOLD_DB -- one continuous region from that
    # crossing to the top of the grid, not perforated by the few points
    # (mode 3/4 peaks) that briefly poke back above it.
    below_thresh = mag_modal_db < MAG_SHADE_THRESHOLD_DB
    if np.any(below_thresh):
        f_shade_start = frequencies_hz[np.argmax(below_thresh)]
        ax_mag.axvspan(f_shade_start, frequencies_hz[-1], color="#e53e3e", alpha=0.10, zorder=0)
        ax_phase.axvspan(f_shade_start, frequencies_hz[-1], color="#e53e3e", alpha=0.10, zorder=0)
        print(f"Shading from {f_shade_start:.1f} Hz (first crossing below {MAG_SHADE_THRESHOLD_DB:.0f} dB) "
              f"to {frequencies_hz[-1]:.0f} Hz.")

    ax_mag.plot(frequencies_hz, mag_modal_db, color="#000000", linewidth=1.8, label="sum of 6 modes")
    ax_mag.set_xscale("log")
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(r"Modal superposition of $x_n(s)/\theta_{cmd}(s)$ (frictionless baseline)" + "\n"
                     f"modes {N_DIMMED+1}-{n} dimmed (>{down_db:.0f} dB below modes 1-{N_DIMMED})",
                     fontsize=10)
    ax_mag.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=6.5, ncol=2, loc="lower left")

    ax_phase.plot(frequencies_hz, phase_modal_deg, color="#000000", linewidth=1.4)
    ax_phase.set_xscale("log")
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_phase.grid(True, which="both", linewidth=0.4, color="#cccccc")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_dir = OUT_DIR / "npz"
    npz_dir.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    out_path = OUT_DIR / "modal_superposition.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_superposition_data.npz",
              omega=omega, freq_hz=freq_hz, zeta=zeta, R=R, R_norm=R_norm,
              frequencies_hz=frequencies_hz, G_modal=G_modal,
              mag_modal_db=mag_modal_db, phase_modal_deg=phase_modal_deg)

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_superposition_data.npz'}")


if __name__ == "__main__":
    main()
