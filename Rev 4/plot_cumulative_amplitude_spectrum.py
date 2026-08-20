#!/usr/bin/env python3
"""Cumulative amplitude spectrum (CAS) of the frictionless model's stepping
tracking error, requested 2026-08-20; revised 2026-08-20 per review: settled
case only, per-mode annotations show percentage only (no nm value), title
trimmed to one line, and the bottom footnote removed.

Signal: e(t) = (L/2pi)*theta_cmd(t) - x_n(t) from the full-step, settled
(250 ms/step) stepping run -- generate_stepping_trajectory.py's panel (b),
the baseline/reference case Figure B in that script is also built around.
Reused via import (run_case_segments, build_error_system, DWELL_SETTLE, DT),
not recomputed by hand, since that machinery is already the validated
source of this signal; generate_stepping_trajectory.py itself never saves
the raw array, only the rendered figures, so this script reruns that one
case rather than modifying that file.

f, S = welch(e, fs, nperseg=...) -> power spectral density (per-segment
mean removed -- scipy's default detrend='constant' -- so this is the AC/
fluctuating error only).
cas(f) = sqrt(reverse-cumsum(S*df)) -- reverse, not forward, because the
mechatronics-convention question is "how much RMS error remains above this
frequency" (directly the residual you'd have left after an ideal controller
with unity gain up to that frequency), not "how much has accumulated below
it."

Per-mode quadrature contribution: since cas(f)^2 is the integral of S from
f to Nyquist, the power inside a band [f_lo, f_hi] is exactly
cas(f_lo)^2 - cas(f_hi)^2. Bands are split at the geometric mean between
consecutive modal frequencies (from plot_modal_superposition.py's saved
eigenanalysis) -- geometric, not arithmetic, because the axis is log.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

from build_bode_rev4 import build_matrices, build_state_space, load_parameters
from generate_stepping_trajectory import DT, DWELL_SETTLE, build_error_system, run_case_segments

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "rendered_assets" / "temp"
NPZ_DIR = OUT_DIR / "npz"
NPERSEG = 4096


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)
    A, B, _Cy = build_state_space(M, C, K, B_u)
    lead_ratio = params["L"] / (2.0 * np.pi)
    sys_err = build_error_system(A, B, lead_ratio)

    decim = 20
    fs = 1.0 / (DT * decim)

    print(f"Running full-step, settled (t_fire={DWELL_SETTLE*1e3:.0f} ms) case, fs={fs:.1f} Hz...")
    _S, e, _CMD, _ends = run_case_segments(sys_err, micro=1, t_fire=DWELL_SETTLE, dt=DT, decim=decim)
    print(f"  {len(e)} samples, {len(e)/fs:.3f} s span")

    f, S_psd = welch(e, fs=fs, nperseg=NPERSEG)
    df = f[1] - f[0]
    cas = np.sqrt(np.cumsum(S_psd[::-1] * df)[::-1])

    static_droop = float(np.mean(e))
    total_rms_cas = cas[0]
    total_rms_direct = float(np.std(e))
    print(f"\nSanity check: cas(f->0) = {total_rms_cas*1e9:.4f} nm  vs.  "
          f"np.std(e) = {total_rms_direct*1e9:.4f} nm  "
          f"(gap {abs(total_rms_cas-total_rms_direct)/total_rms_direct*100:.2f}%)")
    print(f"Static droop (mean e, excluded from PSD by detrend='constant'): {static_droop*1e9:+.3f} nm")

    # ---- Modal frequencies, from plot_modal_superposition.py's saved eigenanalysis ----
    modal_data_path = NPZ_DIR / "modal_superposition_data.npz"
    if modal_data_path.exists():
        freq_hz_modes = np.load(modal_data_path)["freq_hz"]
    else:
        from scipy.linalg import eigh
        lam, _phi = eigh(K, M)
        freq_hz_modes = np.sqrt(lam) / (2.0 * np.pi)
    nyquist = fs / 2.0
    modes_in_band = freq_hz_modes[freq_hz_modes < nyquist]
    modes_above_nyquist = freq_hz_modes[freq_hz_modes >= nyquist]

    # ---- Per-mode quadrature contribution: band edges at the geometric mean
    # between consecutive in-band modal frequencies. ----
    f_lo_edge = f[1]   # first non-zero bin -- the record's actual frequency resolution
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
    np.savez(NPZ_DIR / "cumulative_amplitude_spectrum_data.npz",
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
    ax.set_title("Frictionless model, stepping tracking error -- cumulative amplitude spectrum")
    ax.grid(True, which="both", linewidth=0.4, color="#cccccc")
    ax.legend(fontsize=8, loc="center left")

    fig.tight_layout()

    out_path = OUT_DIR / "cumulative_amplitude_spectrum.png"
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"\nWrote {out_path}")
    print(f"Wrote {NPZ_DIR / 'cumulative_amplitude_spectrum_data.npz'}")


if __name__ == "__main__":
    main()
