#!/usr/bin/env python3
"""Quantify how the multisine BLA's resonance peaks settle with amplitude.

Reads the aggregated multisine_rev42_results.npz (already produced by
run_multisine_rev42.py) and, for each of the two visible resonance peaks,
tracks peak frequency and peak magnitude across the 30 RMS amplitude levels.
Convergence is defined relative to the highest-amplitude (1000% Fs) case:
the amplitude beyond which every subsequent point stays within
MAGNITUDE_TOLERANCE_DB of that asymptotic magnitude.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULT_NPZ = ROOT / 'rendered_assets' / 'npz' / 'multisine_rev42_results.npz'
OUTPUT_PLOT = ROOT / 'rendered_assets' / 'multisine_resonance_convergence.png'
OUTPUT_JSON = ROOT / 'rendered_assets' / 'multisine_resonance_convergence.json'

SOLVER = 'Radau'
BANDS = {
    'primary (~160 Hz)': (50.0, 300.0),
    'secondary (~700-730 Hz)': (600.0, 900.0),
}
MAGNITUDE_TOLERANCE_DB = 0.2


def magnitude_db(response: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(np.abs(response), 1.0e-300))


def peak_in_band(freq: np.ndarray, mag_db: np.ndarray, lo: float, hi: float):
    mask = (freq >= lo) & (freq <= hi)
    f = freq[mask]
    peak_f, peak_m = [], []
    for row in mag_db:
        m = row[mask]
        i = int(np.argmax(m))
        peak_f.append(float(f[i]))
        peak_m.append(float(m[i]))
    return np.asarray(peak_f), np.asarray(peak_m)


def first_converged_amplitude(amplitude, delta_db, tolerance_db):
    """Smallest amplitude beyond which |delta_db| never exceeds tolerance again."""
    ok = np.abs(delta_db) <= tolerance_db
    # Find the last index where it is NOT converged; the answer is the next one.
    bad = np.flatnonzero(~ok)
    if bad.size == 0:
        return float(amplitude[0])
    last_bad = bad[-1]
    if last_bad + 1 >= amplitude.size:
        return None  # never converges within the tested range
    return float(amplitude[last_bad + 1])


def main() -> None:
    data = np.load(RESULT_NPZ)
    amplitude = data['amplitude_percent_Fs']
    freq = data['frequency_hz']
    mag_db = magnitude_db(data[f'{SOLVER}_G_bla'])

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex=True)
    results = {}

    for col, (band_name, (lo, hi)) in enumerate(BANDS.items()):
        peak_f, peak_m = peak_in_band(freq, mag_db, lo, hi)
        final_f, final_m = peak_f[-1], peak_m[-1]
        delta_m = peak_m - final_m
        converge_amp = first_converged_amplitude(
            amplitude, delta_m, MAGNITUDE_TOLERANCE_DB
        )

        ax_f, ax_m = axes[0, col], axes[1, col]
        ax_f.plot(amplitude, peak_f, 'o-', color='#136f63', markersize=3.5)
        ax_f.set_title(f'{band_name}')
        ax_f.set_ylabel('Peak frequency (Hz)')
        ax_f.grid(True, which='both', alpha=0.3)

        ax_m.plot(amplitude, peak_m, 'o-', color='#d1495b', markersize=3.5)
        ax_m.axhline(final_m, color='0.5', linestyle=':', linewidth=1.0)
        ax_m.fill_between(
            amplitude, final_m - MAGNITUDE_TOLERANCE_DB,
            final_m + MAGNITUDE_TOLERANCE_DB, color='0.5', alpha=0.15,
            label=f'±{MAGNITUDE_TOLERANCE_DB:g} dB of 1000% Fs value',
        )
        if converge_amp is not None:
            ax_m.axvline(
                converge_amp, color='#136f63', linestyle='--', linewidth=1.1,
                label=f'converged beyond {converge_amp:g}% Fs',
            )
        ax_m.set_ylabel('Peak magnitude (dB re 1 m/rad)')
        ax_m.set_xlabel('RMS excitation amplitude (% of static friction Fs)')
        ax_m.set_xscale('log')
        ax_m.grid(True, which='both', alpha=0.3)
        ax_m.legend(loc='lower right', fontsize=8)

        results[band_name] = {
            'frequency_hz': peak_f.tolist(),
            'magnitude_db': peak_m.tolist(),
            'final_frequency_hz': final_f,
            'final_magnitude_db': final_m,
            'converged_beyond_amplitude_pct_fs': converge_amp,
            'tolerance_db': MAGNITUDE_TOLERANCE_DB,
        }
        print(f'{band_name}:')
        print(f'  frequency range: {peak_f.min():.2f} - {peak_f.max():.2f} Hz '
              f'(final {final_f:.2f} Hz)')
        print(f'  magnitude range: {peak_m.min():.2f} - {peak_m.max():.2f} dB '
              f'(final {final_m:.2f} dB)')
        if converge_amp is not None:
            print(f'  within +/-{MAGNITUDE_TOLERANCE_DB:g} dB of the 1000% Fs '
                  f'value from {converge_amp:g}% Fs onward')
        else:
            print('  never settles within tolerance across the tested range')
        print()

    fig.suptitle(
        'Rev 4.2 multisine BLA: resonance peak convergence vs. excitation '
        'amplitude', fontsize=13.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTPUT_PLOT, dpi=160)
    plt.close(fig)

    OUTPUT_JSON.write_text(
        json.dumps({
            'amplitude_percent_Fs': amplitude.tolist(),
            'tolerance_db': MAGNITUDE_TOLERANCE_DB,
            'bands': results,
        }, indent=2) + '\n', encoding='utf-8',
    )
    print(f'Plot: {OUTPUT_PLOT}')
    print(f'Summary: {OUTPUT_JSON}')


if __name__ == '__main__':
    main()
