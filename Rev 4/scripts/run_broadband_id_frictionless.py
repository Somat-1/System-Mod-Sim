#!/usr/bin/env python3
"""Sanity check requested 2026-08-20: push the FRICTIONLESS linear baseline
(build_bode_rev4.py's (A, B, C_y), theta_cmd -> x_n, friction-port inputs
held at zero) through the exact same chirp + PRBS + Welch identification
pipeline used for the nonlinear LuGre model (lugre_friction/run_broadband_id.py),
at the SAME 5 command amplitudes, and compare.

Reuses lugre_friction/broadband_signals.py (ChirpSignal, PRBSSignal -- same
sweep/clock/period design) and lugre_friction/broadband_estimators.py
(prbs_estimate, chirp_estimate, fractional_octave_smooth -- same Welch
segment-per-period / full-record-FFT-ratio / 1/24-octave smoothing) rather
than reimplementing them, so this is genuinely the same pipeline, not a
lookalike. Only the plant changes: scipy.signal.lsim on the linear (A, B,
C_y) system instead of solve_ivp on the nonlinear LuGre RHS -- correct and
sufficient for an LTI system, and far cheaper (no stiffness, no friction
kink, no atol tuning).

Expected result for a genuinely linear, time-invariant, noiseless plant:
G(f) = x_n(f)/theta_cmd(f) must be IDENTICAL across all 5 amplitudes
(superposition), and coherence must sit at ~1.0 everywhere. Any amplitude-
dependence or coherence loss here would mean the pipeline itself (signal
generation, Welch/FFT estimation, smoothing) has a bug -- since the LuGre
run showed both amplitude-dependence and structure, this is the cross-
check for whether that behavior is real friction physics or a pipeline
artifact. Also overlaid against the existing analytic Bode
(rendered_assets/bode_rev4_data.npz, built directly from the Laplace-domain
transfer function) as an independent ground truth for the estimator itself.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import StateSpace, lsim

ROOT = Path(__file__).resolve().parent.parent
LUGRE_DIR = ROOT / "lugre_friction"
sys.path.insert(0, str(LUGRE_DIR / "scripts"))

from broadband_estimators import chirp_estimate, fractional_octave_smooth, prbs_estimate  # noqa: E402
from broadband_signals import (  # noqa: E402
    ChirpSignal,
    F_HI_HZ,
    F_LO_HZ,
    FS_HZ,
    PRBSSignal,
    PRBS_DISCARD_PERIODS,
)

from build_bode_rev4 import build_matrices, build_state_space, load_parameters  # noqa: E402

OUT_DIR = ROOT / "rendered_assets" / "temp"
NPZ_DIR = OUT_DIR / "npz"
LUGRE_DATA = LUGRE_DIR / "rendered_assets" / "npz" / "broadband_id_data.npz"
ANALYTIC_BODE_DATA = ROOT / "rendered_assets" / "npz" / "bode_rev4_data.npz"

COH_GATE = 0.9


def simulate(sys_lin: StateSpace, u_full: np.ndarray, dt: float) -> np.ndarray:
    n = len(u_full)
    t = np.arange(n) * dt
    U = np.zeros((n, 3))
    U[:, 0] = u_full   # T_fric_sb, F_fric_way stay zero -- frictionless test
    _tout, yout, _xout = lsim(sys_lin, U, t)
    return yout   # x_n(t) directly -- C_y has a single row


def run_prbs(sys_lin: StateSpace, A: float):
    sig = PRBSSignal(A)
    n = int(round(sig.total_duration_s * FS_HZ)) + 1
    t = np.arange(n) / FS_HZ
    u_full = sig.command_at(t)
    t0 = time.time()
    y_full = simulate(sys_lin, u_full, 1.0 / FS_HZ)
    wall_s = time.time() - t0
    f, G, gamma2, n_segments = prbs_estimate(u_full, y_full, FS_HZ, sig.samples_per_period, PRBS_DISCARD_PERIODS)
    return dict(f=f, G=G, gamma2=gamma2, n_segments=n_segments, wall_s=wall_s)


def run_chirp(sys_lin: StateSpace, A: float):
    sig = ChirpSignal(A)
    n = int(round(sig.total_duration_s * FS_HZ)) + 1
    t = np.arange(n) / FS_HZ
    u_full = sig.command_at(t)
    t0 = time.time()
    y_full = simulate(sys_lin, u_full, 1.0 / FS_HZ)
    wall_s = time.time() - t0
    prehold_samples = int(round(1.0 * FS_HZ))
    u_trim, y_trim = u_full[prehold_samples:], y_full[prehold_samples:]
    f, G = chirp_estimate(u_trim, y_trim, FS_HZ)
    f_s, mag_s, phase_s = fractional_octave_smooth(f, G, frac=24.0)
    return dict(f=f_s, mag_db=mag_s, phase_deg=phase_s, wall_s=wall_s)


def main() -> None:
    p = load_parameters()
    M, C, K, B_u = build_matrices(p)
    A_mat, B_mat, C_y = build_state_space(M, C, K, B_u)
    sys_lin = StateSpace(A_mat, B_mat, C_y, np.zeros((1, 3)))

    amplitudes = np.load(LUGRE_DATA)["amplitudes"]   # same 5 amplitudes as the LuGre run
    print(f"Amplitudes (same 5 as the LuGre run): {amplitudes*1e3} mrad RMS")

    prbs_results, chirp_results = [], []
    for A in amplitudes:
        rp = run_prbs(sys_lin, float(A))
        rc = run_chirp(sys_lin, float(A))
        print(f"  A={A*1e3:9.4f} mrad: PRBS {rp['wall_s']:.2f}s ({rp['n_segments']} segments), "
              f"chirp {rc['wall_s']:.2f}s")
        prbs_results.append(rp)
        chirp_results.append(rc)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    npz_payload = {"amplitudes": amplitudes}
    for i, A in enumerate(amplitudes):
        npz_payload[f"prbs_f_amp{i}"] = prbs_results[i]["f"]
        npz_payload[f"prbs_G_amp{i}"] = prbs_results[i]["G"]
        npz_payload[f"prbs_gamma2_amp{i}"] = prbs_results[i]["gamma2"]
        npz_payload[f"chirp_f_amp{i}"] = chirp_results[i]["f"]
        npz_payload[f"chirp_mag_db_amp{i}"] = chirp_results[i]["mag_db"]
        npz_payload[f"chirp_phase_deg_amp{i}"] = chirp_results[i]["phase_deg"]
    np.savez(NPZ_DIR / "broadband_id_frictionless_data.npz", **npz_payload)

    d_analytic = np.load(ANALYTIC_BODE_DATA)
    f_analytic, mag_analytic, phase_analytic = (
        d_analytic["frequencies_hz"], d_analytic["magnitude_db"], d_analytic["phase_deg"],
    )
    mask_a = (f_analytic >= F_LO_HZ) & (f_analytic <= F_HI_HZ)
    phase_analytic_wrapped = ((phase_analytic[mask_a] + 180.0) % 360.0) - 180.0

    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, len(amplitudes) - 1)) for i in range(len(amplitudes))]

    fig, (ax_mag, ax_phase, ax_coh) = plt.subplots(3, 1, figsize=(10.0, 11.0), sharex=True)
    for i, A in enumerate(amplitudes):
        color = colors[i]
        rp, rc = prbs_results[i], chirp_results[i]
        mp = (rp["f"] >= F_LO_HZ) & (rp["f"] <= F_HI_HZ)
        mc = (rc["f"] >= F_LO_HZ) & (rc["f"] <= F_HI_HZ)
        fp, fc_ = rp["f"][mp], rc["f"][mc]

        prbs_mag_db = 20.0 * np.log10(np.maximum(np.abs(rp["G"][mp]), 1e-300))
        prbs_phase_deg = ((np.angle(rp["G"][mp]) * 180.0 / np.pi + 180.0) % 360.0) - 180.0
        chirp_phase_deg = ((rc["phase_deg"][mc] + 180.0) % 360.0) - 180.0
        low_coh = rp["gamma2"][mp] < COH_GATE

        ax_mag.plot(fp, np.ma.masked_where(low_coh, prbs_mag_db), color=color,
                    linewidth=1.0, label=f"A={A*1e3:.3f} mrad RMS (PRBS)")
        ax_mag.plot(fp, np.ma.masked_where(~low_coh, prbs_mag_db), color="#bbbbbb",
                    linewidth=0.8, zorder=0)
        ax_mag.plot(fc_, rc["mag_db"][mc], color=color, linewidth=1.0, linestyle="--")

        ax_phase.plot(fp, np.ma.masked_where(low_coh, prbs_phase_deg), color=color,
                      linewidth=0.0, marker=".", markersize=1.5)
        ax_phase.plot(fp, np.ma.masked_where(~low_coh, prbs_phase_deg), color="#bbbbbb",
                      linewidth=0.0, marker=".", markersize=1.2, zorder=0)
        ax_phase.plot(fc_, chirp_phase_deg, color=color, linewidth=0.0, marker=".", markersize=1.5)

        ax_coh.plot(fp, rp["gamma2"][mp], color=color, linewidth=1.0)

    ax_mag.plot(f_analytic[mask_a], mag_analytic[mask_a], color="#333333", linewidth=1.2,
                linestyle="-.", label="analytic (Laplace-domain, bode_rev4.py)", zorder=5)
    ax_phase.plot(f_analytic[mask_a], phase_analytic_wrapped, color="#333333", linewidth=0.0,
                  marker=".", markersize=1.0, zorder=5)

    ax_coh.axhline(COH_GATE, color="#333333", linestyle=":", linewidth=1.0)

    ax_mag.set_xscale("log")
    ax_mag.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("Frictionless x_n / theta_cmd magnitude -- solid=PRBS (grey where coherence<0.9), "
                      "dashed=chirp, dash-dot=analytic\nsame chirp+PRBS+Welch pipeline as the LuGre run, "
                      "all 5 amplitude curves should overlap exactly (linear system)")
    ax_mag.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=7, ncol=2)

    ax_phase.set_xscale("log")
    ax_phase.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_phase.set_ylim(-180.0, 180.0)
    ax_phase.set_yticks([-180, -90, 0, 90, 180])
    ax_phase.set_ylabel("Phase (deg, wrapped)")
    ax_phase.set_title("Frictionless x_n / theta_cmd phase")
    ax_phase.grid(True, which="both", linewidth=0.4, color="#cccccc")

    ax_coh.set_xscale("log")
    ax_coh.set_xlim(F_LO_HZ, F_HI_HZ)
    ax_coh.set_xlabel("Frequency (Hz)")
    ax_coh.set_ylabel("Coherence gamma^2 (PRBS only)")
    ax_coh.set_title("Coherence -- should sit at ~1.0 everywhere (linear, noiseless)")
    ax_coh.set_ylim(0.0, 1.05)
    ax_coh.grid(True, which="both", linewidth=0.4, color="#cccccc")

    fig.suptitle("Rev 4 frictionless baseline: same chirp+PRBS+Welch pipeline as the LuGre broadband ID")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    montage_path = OUT_DIR / "broadband_id_frictionless_montage.png"
    fig.savefig(montage_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {montage_path}")
    print(f"Wrote {NPZ_DIR / 'broadband_id_frictionless_data.npz'}")

    # ---- Amplitude-independence check: max spread across amplitudes at each freq ----
    mag_stack = np.array([
        20.0 * np.log10(np.maximum(np.abs(prbs_results[i]["G"]), 1e-300)) for i in range(len(amplitudes))
    ])
    spread_db = mag_stack.max(axis=0) - mag_stack.min(axis=0)
    print(f"\nPRBS magnitude spread across the 5 amplitudes: max {spread_db.max():.4f} dB, "
          f"median {np.median(spread_db):.4f} dB (should be ~0 for a linear system)")
    coh_stack = np.array([prbs_results[i]["gamma2"] for i in range(len(amplitudes))])
    print(f"Coherence: min {coh_stack.min():.5f}, median {np.median(coh_stack):.5f} "
          f"(should be ~1.0 everywhere)")


if __name__ == "__main__":
    main()
