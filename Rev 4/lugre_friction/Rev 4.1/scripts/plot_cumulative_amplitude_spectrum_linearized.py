#!/usr/bin/env python3
"""Cumulative amplitude spectrum (CAS) of the frozen-linearized (LuGre
local-linearization) system's stepping tracking error, requested 2026-08-21
-- mirror of ../../scripts/plot_cumulative_amplitude_spectrum.py, applied to
run_local_linearization_bode.py's single-input K/C/M/B_em instead of the
frictionless baseline's three-input system.

Reuses ../../scripts/generate_stepping_trajectory.py's build_error_system
(already generic in B's column count -- it only ever touches column 0,
theta_cmd) and MOVES/FULL_STEP/edge_list/DT/DWELL_SETTLE (the commanded
step sequence itself is a property of the trajectory generator, not the
plant). run_case_segments is NOT reused as-is: that function hardcodes a
3-column input array (u = np.zeros((n_hold, 3))) to match the frictionless
baseline's 3-input B_u (theta_cmd, T_fric_sb, F_fric_way -- both held at
zero). The frozen-linearized system has already absorbed T_fric_sb and
F_fric_way into K/C as static feedback (see run_local_linearization_bode.py
docstring), so B_em is a single column; run_case_segments_single_input below
is the same lsim-per-dwell logic with u = np.zeros((n_hold, 1)) instead.

Same PSD/CAS/per-mode-band methodology as the frictionless version:
cas(f) = sqrt(reverse-cumsum(S*df)) from scipy.signal.welch, per-mode
quadrature contribution cas(f_lo)^2 - cas(f_hi)^2 at geometric-mean band
edges between this system's OWN modal frequencies (from
modal_superposition_linearized_data.npz, not the frictionless baseline's).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.signal import StateSpace, lsim, welch

from lugre_model import load_parameters
from run_local_linearization_bode import build_linearized_matrices, build_state_space

LUGRE_DIR = Path(__file__).resolve().parent.parent
REV4_DIR = LUGRE_DIR.parents[1]
OUT_DIR = LUGRE_DIR / "rendered_assets" / "temp"
NPZ_DIR = OUT_DIR / "npz"

sys.path.insert(0, str(REV4_DIR / "scripts"))
from generate_stepping_trajectory import DT, DWELL_SETTLE, build_error_system, edge_list  # noqa: E402

NPERSEG = 4096


def run_case_segments_single_input(sys_ss: StateSpace, micro: int, t_fire: float, dt: float, decim: int = 20):
    """Same per-dwell lsim-with-carried-state logic as
    generate_stepping_trajectory.run_case_segments, generalized to a
    single-column B (u = np.zeros((n_hold, 1))) instead of that function's
    hardcoded 3 columns."""
    n_hold = round(t_fire / dt)
    assert abs(t_fire / dt - n_hold) < 1e-9, "edges must land on samples"

    edges = edge_list(micro)
    t_seg = np.arange(n_hold) * dt
    z = np.zeros(sys_ss.A.shape[0])
    S, E, CMD, ends = [], [], [], []

    for j, theta_cmd in enumerate(edges):
        u = np.zeros((n_hold, 1))
        u[:, 0] = theta_cmd
        _, y, x = lsim(sys_ss, u, t_seg, X0=z)
        z = x[-1]
        s = (j + np.arange(n_hold) / n_hold) / micro
        S.append(s[::decim])
        E.append(y[::decim])
        CMD.append(np.full(n_hold, theta_cmd)[::decim])
        ends.append(y[-1])

    return np.concatenate(S), np.concatenate(E), np.concatenate(CMD), np.array(ends)


def main() -> None:
    params = load_parameters()
    M, K, C, B_em = build_linearized_matrices(params)
    A, B, _Cy = build_state_space(M, K, C, B_em)
    lead_ratio = params["L"] / (2.0 * np.pi)
    B_col = B.reshape(-1, 1)   # build_error_system expects a 2-D B (uses B.shape[1])
    sys_err = build_error_system(A, B_col, lead_ratio)

    decim = 20
    fs = 1.0 / (DT * decim)

    print(f"Running full-step, settled (t_fire={DWELL_SETTLE*1e3:.0f} ms) case, fs={fs:.1f} Hz...")
    _S, e, _CMD, _ends = run_case_segments_single_input(sys_err, micro=1, t_fire=DWELL_SETTLE, dt=DT, decim=decim)
    print(f"  {len(e)} samples, {len(e)/fs:.3f} s span")

    f, S_psd = welch(e, fs=fs, nperseg=NPERSEG)
    df = f[1] - f[0]
    cas = np.sqrt(np.cumsum(S_psd[::-1] * df)[::-1])

    static_droop = float(np.mean(e))
    total_rms_cas = cas[0]
    total_rms_direct = float(np.std(e))
    gap_pct = abs(total_rms_cas - total_rms_direct) / total_rms_direct * 100.0
    print(f"\nSanity check: cas(f->0) = {total_rms_cas*1e9:.4f} nm  vs.  "
          f"np.std(e) = {total_rms_direct*1e9:.4f} nm  (gap {gap_pct:.2f}%)")
    print(f"Static droop (mean e, excluded from PSD by detrend='constant'): {static_droop*1e9:+.3f} nm")

    # The frictionless baseline's equivalent gap is ~0.6% (Parseval's
    # theorem, as expected for a stationary-ish signal). A gap this large
    # means the "settled" DWELL_SETTLE=250 ms dwell -- a constant tuned to
    # the FRICTIONLESS system's damping, reused unchanged here -- does NOT
    # actually settle this system: mode 1's damping collapsed from
    # zeta=0.0200 (frictionless) to zeta=0.0008 here (a ~25x drop), pushing
    # its decay time constant from ~45 ms to ~500 ms and its ~4-tau settling
    # time from ~180 ms to ~2.0 s -- 8x longer than the 250 ms dwell. The
    # signal is therefore still ringing from one step when the next one
    # lands, which is non-stationary, exactly what Welch's PSD (and the
    # reverse-cumsum CAS built from it) assumes away. The huge low-frequency
    # (<10 Hz) CAS content below is very likely this artifact, not a real
    # sub-10 Hz structural mode -- there isn't one.
    GAP_WARN_THRESHOLD_PCT = 5.0
    if gap_pct > GAP_WARN_THRESHOLD_PCT:
        print(f"  WARNING: gap {gap_pct:.1f}% >> frictionless baseline's ~0.6% -- the signal is very "
              "likely non-stationary within the 250 ms dwell (mode 1's zeta dropped ~25x vs. the "
              "frictionless case, so its settling time, ~2.0 s, is ~8x the dwell). Treat this CAS, "
              "especially its low-frequency content, as unreliable rather than a genuine steady-"
              "periodic error spectrum -- see script docstring/comments for the full diagnosis.")

    # ---- Modal frequencies, from plot_modal_superposition_linearized.py's saved eigenanalysis ----
    modal_data_path = NPZ_DIR / "modal_superposition_linearized_data.npz"
    if modal_data_path.exists():
        freq_hz_modes = np.load(modal_data_path)["freq_hz"]
    else:
        lam, _phi = eigh(K, M)
        freq_hz_modes = np.sqrt(lam) / (2.0 * np.pi)
    nyquist = fs / 2.0
    modes_in_band = freq_hz_modes[freq_hz_modes < nyquist]
    modes_above_nyquist = freq_hz_modes[freq_hz_modes >= nyquist]

    f_lo_edge = f[1]
    edges = [f_lo_edge]
    for j in range(len(modes_in_band) - 1):
        edges.append(np.sqrt(modes_in_band[j] * modes_in_band[j + 1]))
    edges.append(nyquist)

    cas_at_edges = np.interp(edges, f, cas)
    band_power = cas_at_edges[:-1] ** 2 - cas_at_edges[1:] ** 2
    band_nm = np.sqrt(np.maximum(band_power, 0.0)) * 1e9
    band_pct = band_power / (total_rms_cas ** 2) * 100.0

    print("\nPer-mode quadrature contribution:")
    for j, fm in enumerate(modes_in_band):
        print(f"  mode {j+1} ({fm:7.1f} Hz): {band_nm[j]:8.2f} nm  ({band_pct[j]:5.1f} % of total)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(NPZ_DIR / "cumulative_amplitude_spectrum_linearized_data.npz",
              f=f, S_psd=S_psd, cas=cas, e=e, fs=fs, static_droop=static_droop,
              band_edges=edges, band_nm=band_nm, band_pct=band_pct, modes_in_band=modes_in_band)

    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    ax.plot(f, cas * 1e9, color="#2b6cb0", linewidth=1.4,
            label=f"settled, {DWELL_SETTLE*1e3:.0f} ms/step (df={df:.2f} Hz)")

    y_top = cas[0] * 1e9 * 1.15
    trans = ax.get_xaxis_transform()
    for j, fm in enumerate(modes_in_band):
        ax.axvline(fm, color="#555555", linestyle="--", linewidth=0.8, zorder=0)
        ax.text(fm, 0.97, f"mode {j+1}, {fm:.1f} Hz", transform=trans, fontsize=7,
                color="#555555", ha="right", va="top", rotation=90, clip_on=False)
        mid_x = np.sqrt(edges[j] * edges[j + 1])
        ax.annotate(f"{band_pct[j]:.1f}%",
                    xy=(mid_x, np.interp(mid_x, f, cas) * 1e9),
                    xytext=(0, 10), textcoords="offset points", fontsize=7,
                    color="#2b6cb0", ha="center", va="bottom")

    ax.set_xscale("log")
    ax.set_xlim(f_lo_edge, nyquist)
    ax.set_ylim(0.0, y_top)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Cumulative RMS error above f (nm)")
    title = "Frozen-linearized LuGre system, stepping tracking error -- cumulative amplitude spectrum"
    if gap_pct > GAP_WARN_THRESHOLD_PCT:
        title += (f"\nCAUTION: Parseval gap {gap_pct:.0f}% (frictionless ~0.6%) -- mode 1 settling\n"
                  "time ~2.0 s >> 250 ms dwell; low-f content likely a non-stationarity artifact")
    ax.set_title(title, fontsize=9.5)
    ax.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, loc="center left")

    fig.tight_layout()

    out_path = OUT_DIR / "cumulative_amplitude_spectrum_linearized.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {out_path}")
    print(f"Wrote {NPZ_DIR / 'cumulative_amplitude_spectrum_linearized_data.npz'}")


if __name__ == "__main__":
    main()
