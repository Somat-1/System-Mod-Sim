#!/usr/bin/env python3
"""Frequency-response estimators for run_broadband_id.py.

PRBS -> Welch, segment length = exactly one PRBS period, noverlap=0,
window='boxcar'. Because the excitation is periodic and segments align
exactly with its period, there is no spectral leakage and no windowing
penalty -- a Hann window here would only throw away information. Coherence
falls out of the same Welch cross/auto spectra:
    G(f)      = Pxy(f) / Pxx(f)
    gamma^2(f) = |Pxy(f)|^2 / (Pxx(f) * Pyy(f))
With a noiseless simulation, gamma^2 < 1 can only come from the response
not repeating identically period-to-period -- i.e. genuine nonlinear
distortion (harmonic/subharmonic generation, incomplete settling), not
measurement noise. It is therefore a direct nonlinearity detector, not
just a data-quality flag.

Chirp -> full-record FFT ratio, NOT Welch. Segmenting a chirp the way PRBS
is segmented would put only a narrow instantaneous band in each segment,
so Pxx(f) ~= 0 across most of the spectrum in every segment and Pxy/Pxx
becomes division by noise almost everywhere. G(f) = Y(f)/U(f) is computed
once over the entire record instead, then smoothed with a fractional-
octave (1/24 octave) moving average -- constant-Hz smoothing would
over-smooth the low end (where the sweep lingers) and under-smooth the
high end (where it's already sparse in time). Repeating the identical
chirp run and averaging (as one would with real measurement noise) buys
nothing here: this simulation is deterministic, so repeats are bit-for-bit
identical and there is no noise to average out. No chirp-side coherence is
produced for the same reason -- coherence needs independent realizations
to be anything other than trivially 1.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import csd, welch


def prbs_estimate(u: np.ndarray, y: np.ndarray, fs: float, samples_per_period: int,
                   n_discard_periods: int):
    """u, y: full (undiscarded) PRBS input/output records at fs. Returns
    (f, G, gamma2, n_segments_used)."""
    start = n_discard_periods * samples_per_period
    n_segments = (len(u) - start) // samples_per_period
    end = start + n_segments * samples_per_period
    u_use, y_use = u[start:end], y[start:end]

    f, Pxx = welch(u_use, fs=fs, window="boxcar", nperseg=samples_per_period, noverlap=0)
    _, Pyy = welch(y_use, fs=fs, window="boxcar", nperseg=samples_per_period, noverlap=0)
    _, Pxy = csd(u_use, y_use, fs=fs, window="boxcar", nperseg=samples_per_period, noverlap=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        G = Pxy / Pxx
        gamma2 = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    return f, G, gamma2, n_segments


def chirp_estimate(u: np.ndarray, y: np.ndarray, fs: float):
    """u, y: full chirp input/output records (transient pre-hold already
    stripped) at fs. Returns (f, G) from a single full-record FFT ratio."""
    n = len(u)
    U = np.fft.rfft(u)
    Y = np.fft.rfft(y)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    with np.errstate(divide="ignore", invalid="ignore"):
        G = Y / U
    return f, G


def fractional_octave_smooth(f: np.ndarray, G: np.ndarray, frac: float = 24.0,
                              f_out: np.ndarray | None = None, n_out: int = 500):
    """Smooth complex G(f) (magnitude in dB, phase unwrapped, smoothed
    separately -- averaging raw complex values would let phase variation
    within the band cancel the magnitude) over a 1/frac-octave window, at
    a log-spaced set of output frequencies (standard for fractional-octave
    reporting; also far cheaper than smoothing at the original FFT bin
    resolution). f must be sorted ascending, as rfftfreq's output is."""
    mag_db = 20.0 * np.log10(np.maximum(np.abs(G), 1e-300))
    phase_deg = np.unwrap(np.angle(G)) * 180.0 / np.pi

    if f_out is None:
        f_lo = f[1] if f[0] <= 0.0 else f[0]   # skip DC
        f_out = np.logspace(np.log10(f_lo), np.log10(f[-1]), n_out)

    log2f = np.log2(np.maximum(f, 1e-300))
    half_width = 1.0 / (2.0 * frac)

    mag_out = np.empty_like(f_out)
    phase_out = np.empty_like(f_out)
    for i, fc in enumerate(f_out):
        log2fc = np.log2(fc)
        lo = np.searchsorted(log2f, log2fc - half_width, side="left")
        hi = np.searchsorted(log2f, log2fc + half_width, side="right")
        if hi <= lo:
            idx = min(max(np.searchsorted(log2f, log2fc), 0), len(log2f) - 1)
            mag_out[i], phase_out[i] = mag_db[idx], phase_deg[idx]
        else:
            mag_out[i] = np.mean(mag_db[lo:hi])
            phase_out[i] = np.mean(phase_deg[lo:hi])

    return f_out, mag_out, phase_out
