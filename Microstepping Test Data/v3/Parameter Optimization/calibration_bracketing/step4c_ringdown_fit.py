#!/usr/bin/env python3
"""Step 4c (from the reviewer's reframing of 4b): fit the post-move
ring-down visible right after the C block's first approach move, across
all six runs, to get a direct measurement of the dominant drivetrain
mode's frequency and damping ratio -- currently a 2%-critical placeholder
(`sigma1_target_zeta = 0.7` is the LuGre bristle-damping target, NOT this
structural mode) everywhere in the parameter table, and the number that
governs how hard detent forcing rings (Figure B's ~140 Hz feature).

Model: y(t) = y_inf + A*exp(-zeta*wn*t)*cos(wd*t - phi), wd = wn*sqrt(1-zeta^2)
fit by nonlinear least squares to the measured position residual after the
approach move's main jump.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

SCRIPTS_DIR = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Microstepping Test Data\v3\scripts")
sys.path.insert(0, str(SCRIPTS_DIR))
from plot_block_montage import parse_ids, load_log_rows, find_block, IDS_PATH, LOG_PATH

OUT_DIR = Path(__file__).resolve().parent
WINDOW_S = 0.30  # look 300 ms after the jump


def damped_sinusoid(t, y_inf, A, zeta, wn, phi):
    wd = wn * np.sqrt(max(1.0 - zeta ** 2, 1.0e-6))
    return y_inf + A * np.exp(-zeta * wn * t) * np.cos(wd * t - phi)


RAMP_DURATION_S = 4.0 / 150.0  # 4 full steps at CONDITIONING_FULL_STEPS_S=150, no accel phase
FREQ_BAND_HZ = (50.0, 500.0)  # restrict the search to the plausible structural-mode band


def extract_ringdown(run_index, time_s, position_nm, sample_period_s, rows):
    start_s, end_s = find_block(rows, run_index, 'C')
    approach = next(
        r for r in rows if r['event'] == 'MOVE_ACK'
        and r['run_index'] == str(run_index) and r['label'] == 'positive_approach'
        and start_s <= float(r['ids_time_s']) <= end_s
    )
    move_abs_t = float(approach['ids_time_s'])
    idx0 = int(round((move_abs_t - 0.05) / sample_period_s))
    idx1 = int(round((move_abs_t + WINDOW_S) / sample_period_s))
    window = position_nm[idx0:idx1] / 1000.0
    t_window = np.arange(idx0, idx1) * sample_period_s - move_abs_t

    # find the jump (steepest slope) as the ramp's start
    dy = np.diff(window)
    jump_idx = int(np.argmax(np.abs(dy)))
    # the free ring-down starts once the commanded ramp itself has finished
    # (RAMP_DURATION_S later) -- fitting from the jump itself would conflate
    # the ramp's own shape with the drivetrain's free response after it.
    ramp_end_idx = jump_idx + int(round(RAMP_DURATION_S / sample_period_s))
    baseline = window[max(0, jump_idx - 20):jump_idx].mean()
    t_ring = t_window[ramp_end_idx:] - t_window[ramp_end_idx]
    y_ring = window[ramp_end_idx:] - baseline
    return t_ring, y_ring


def fit_ringdown(t_ring, y_ring):
    y_inf0 = np.median(y_ring[-max(1, len(y_ring) // 10):])
    resid0 = y_ring - y_inf0
    # initial frequency guess from a quick FFT of the residual, restricted
    # to a plausible structural-mode band so a slow settling trend (left
    # over from the ramp) can't be mistaken for the oscillation itself.
    dt = np.median(np.diff(t_ring))
    spectrum = np.abs(np.fft.rfft(resid0 * np.hanning(len(resid0))))
    freqs = np.fft.rfftfreq(len(resid0), d=dt)
    valid = (freqs >= FREQ_BAND_HZ[0]) & (freqs <= FREQ_BAND_HZ[1])
    f0_guess = freqs[valid][np.argmax(spectrum[valid])] if np.any(valid) else 150.0
    wn0 = 2.0 * np.pi * f0_guess
    A0 = resid0[0] if abs(resid0[0]) > 1.0e-6 else np.max(np.abs(resid0))

    p0 = [y_inf0, A0, 0.05, wn0, 0.0]
    bounds = ([-np.inf, -np.inf, 0.0, 2 * np.pi * FREQ_BAND_HZ[0], -np.pi],
              [np.inf, np.inf, 0.9, 2 * np.pi * FREQ_BAND_HZ[1], np.pi])
    popt, _ = curve_fit(damped_sinusoid, t_ring, y_ring, p0=p0, bounds=bounds, maxfev=20000)
    y_inf, A, zeta, wn, phi = popt
    wd = wn * np.sqrt(max(1.0 - zeta ** 2, 1.0e-6))
    return {
        'y_inf': y_inf, 'A': A, 'zeta': zeta, 'wn': wn, 'wd': wd, 'phi': phi,
        'f_n_hz': wn / (2.0 * np.pi), 'f_d_hz': wd / (2.0 * np.pi),
    }


def main():
    time_s, position_nm, start_epoch_s, sample_period_s = parse_ids(IDS_PATH)
    rows = load_log_rows(LOG_PATH, start_epoch_s)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.0), sharex=True)
    results = {}
    print(f"{'run':>4s} {'f_n (Hz)':>10s} {'f_d (Hz)':>10s} {'zeta':>8s} {'A (um)':>8s} {'y_inf (um)':>10s}")
    for run_index, ax in zip(range(1, 7), axes.flat):
        t_ring, y_ring = extract_ringdown(run_index, time_s, position_nm, sample_period_s, rows)
        try:
            fit = fit_ringdown(t_ring, y_ring)
        except RuntimeError as exc:
            print(f'{run_index:>4d}  FIT FAILED: {exc}')
            ax.plot(t_ring * 1000, y_ring, color='#1f77b4', lw=1.0)
            ax.set_title(f'Run {run_index}: fit failed')
            continue
        results[run_index] = fit
        print(f"{run_index:>4d} {fit['f_n_hz']:>10.2f} {fit['f_d_hz']:>10.2f} "
              f"{fit['zeta']:>8.4f} {fit['A']:>8.3f} {fit['y_inf']:>10.3f}")

        t_fine = np.linspace(t_ring[0], t_ring[-1], 2000)
        y_fit = damped_sinusoid(t_fine, fit['y_inf'], fit['A'], fit['zeta'], fit['wn'], fit['phi'])
        ax.plot(t_ring * 1000, y_ring, color='#1f77b4', lw=0.9, label='measured')
        ax.plot(t_fine * 1000, y_fit, color='#d62728', lw=1.3, linestyle='--', label='damped-sinusoid fit')
        ax.set_title(f"Run {run_index}: f_n={fit['f_n_hz']:.1f} Hz, zeta={fit['zeta']:.3f}", fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    for ax in axes[-1]:
        ax.set_xlabel('Time since commanded ramp end (ms)')
    for ax in axes[:, 0]:
        ax.set_ylabel('Position (µm, baselined)')

    if results:
        f_ns = [r['f_n_hz'] for r in results.values()]
        zetas = [r['zeta'] for r in results.values()]
        print(f"\nAcross {len(results)} runs: f_n = {np.mean(f_ns):.1f} +/- {np.std(f_ns):.1f} Hz "
              f"(median {np.median(f_ns):.1f}), zeta = {np.mean(zetas):.4f} +/- {np.std(zetas):.4f} "
              f"(median {np.median(zetas):.4f})")
        print(f"For comparison: sigma1_target_zeta placeholder = 0.7 "
              f"(this is a different quantity -- LuGre bristle damping ratio, not this structural mode)")

    fig.suptitle('Step 4c -- C-block approach ring-down: damped-sinusoid fit, all 6 runs', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_path = OUT_DIR / 'step4c_ringdown_fit.png'
    fig.savefig(out_path, dpi=150)
    print(f'\nSaved {out_path}')


if __name__ == '__main__':
    main()
