#!/usr/bin/env python3
"""Modal superposition of the frozen-linearized (LuGre local-linearization)
system's command-to-stage transfer function x_n(s)/theta_cmd(s), requested
2026-08-21 -- mirror of ../../scripts/plot_modal_superposition.py, applied
to run_local_linearization_bode.py's single-input K/C/M/B_em (frozen bristle
equivalent stiffness/damping at V_STAGE=5 mm/s) instead of the frictionless
baseline's three-input B_u.

x_n(s)/theta_cmd(s) = C_y (M s^2 + C s + K)^-1 B_em is expanded as a sum of
six independent second-order (mass-normalized-mode) systems, same
proportional-damping approximation (zeta_j = diag(phi.T @ C @ phi)_j /
(2*omega_j)) as the frictionless version. The analytic G(jw) (direct
(sI-A)^-1 solve, laplace_transfer_function cross-imported from
../../scripts/build_bode_rev4.py -- that function is fully generic, no
frictionless-specific assumption in it) drives the mismatch shading but is
not itself plotted, same convention as the frictionless version.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from lugre_model import load_parameters
from run_local_linearization_bode import STATE_LABELS, build_linearized_matrices, build_state_space

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parents[1]
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import laplace_transfer_function  # noqa: E402

F_LO_HZ = 10.0
F_HI_HZ = 10000.0
N_FREQ = 4000
# The frictionless baseline hardcodes -110 dB, tuned to ITS OWN DC gain
# (~-77 dB, comfortably above -110). This system's low-frequency gain sits
# at ~-114 dB -- already below a flat -110 dB threshold before the sweep
# even reaches mode 1's resonance -- so the shading is computed as a fixed
# headroom BELOW this system's own peak magnitude instead of a transplanted
# absolute number, keeping the same INTENT (shade only the genuinely-deep,
# near-noise-floor antiresonances) rather than the same literal constant.
SHADE_HEADROOM_DB = 100.0
N_DIMMED = 3


def main() -> None:
    params = load_parameters()
    M, K, C, B_em = build_linearized_matrices(params)
    n = M.shape[0]

    C_y = np.zeros((1, n))
    C_y[0, STATE_LABELS.index("x_n")] = 1.0
    b_col = B_em   # single input: theta_cmd only

    lam, phi = eigh(K, M)
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

    A_fo, B_fo, C_y_fo = build_state_space(M, K, C, B_em)
    G_analytic = np.array([laplace_transfer_function(A_fo, B_fo, C_y_fo, s) for s in s_grid])
    mag_analytic_db = 20.0 * np.log10(np.maximum(np.abs(G_analytic), 1e-15))
    mag_diff = np.abs(mag_modal_db - mag_analytic_db)
    i_max = int(np.argmax(mag_diff))
    print(f"\nModal-sum vs. analytic G(jw): max |magnitude diff| = {mag_diff.max():.4f} dB "
          f"(at {frequencies_hz[i_max]:.1f} Hz: modal={mag_modal_db[i_max]:.1f} dB, "
          f"analytic={mag_analytic_db[i_max]:.1f} dB), median = {np.median(mag_diff):.4f} dB")
    print("  NOTE: max diff is much larger than the frictionless baseline's (~4.5 dB) -- the "
          "proportional-damping approximation (zeta_j from the DIAGONAL of phi.T@C@phi only) is "
          "worse here because two modes are heavily overdamped/non-proportional (sb-dominated "
          "zeta up to ~21, see plot_pole_map-style analysis), and the peak error concentrates "
          "exactly where the true analytic curve dives into a deep antiresonance the 6-mode sum "
          "doesn't reproduce -- both curves are near/below the numerical noise floor there, so a "
          "large dB gap at a near-zero magnitude is expected, not alarming on its own.")

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

    mag_shade_threshold_db = mag_modal_db.max() - SHADE_HEADROOM_DB
    below_thresh = mag_modal_db < mag_shade_threshold_db
    # Find the LAST contiguous below-threshold run (the true high-frequency
    # noise-floor tail), not just the first True -- unlike the frictionless
    # system, this one's low-frequency plateau can itself dip transiently
    # near the threshold before climbing into mode 1's resonance, and that
    # is not the region this shading is meant to flag.
    if np.any(below_thresh):
        run_start = len(below_thresh)
        for k in range(len(below_thresh) - 1, -1, -1):
            if below_thresh[k]:
                run_start = k
            else:
                break
        f_shade_start = frequencies_hz[run_start]
        ax_mag.axvspan(f_shade_start, frequencies_hz[-1], color="#e53e3e", alpha=0.10, zorder=0)
        ax_phase.axvspan(f_shade_start, frequencies_hz[-1], color="#e53e3e", alpha=0.10, zorder=0)
        print(f"Shading from {f_shade_start:.1f} Hz (start of the trailing run below "
              f"{mag_shade_threshold_db:.0f} dB = peak - {SHADE_HEADROOM_DB:.0f} dB) "
              f"to {frequencies_hz[-1]:.0f} Hz.")

    ax_mag.plot(frequencies_hz, mag_modal_db, color="#000000", linewidth=1.8, label="sum of 6 modes")
    ax_mag.set_xscale("log")
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(r"Modal superposition of $x_n(s)/\theta_{cmd}(s)$ (frozen-linearized LuGre system)" + "\n"
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
    out_path = OUT_DIR / "modal_superposition_linearized.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)

    np.savez(npz_dir / "modal_superposition_linearized_data.npz",
              omega=omega, freq_hz=freq_hz, zeta=zeta, R=R, R_norm=R_norm,
              frequencies_hz=frequencies_hz, G_modal=G_modal,
              mag_modal_db=mag_modal_db, phase_modal_deg=phase_modal_deg)

    print(f"\nWrote {out_path}")
    print(f"Wrote {npz_dir / 'modal_superposition_linearized_data.npz'}")


if __name__ == "__main__":
    main()
