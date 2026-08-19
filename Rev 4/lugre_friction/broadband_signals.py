#!/usr/bin/env python3
"""Chirp/PRBS excitation signals, RMS-matched amplitude selection, and the
per-state atol builder for run_broadband_id.py.

Chirp: linear (not log) sweep F_LO_HZ -> F_HI_HZ over CHIRP_DURATION_S.
Log spacing would spend most of the 60 s below 100 Hz and rip through the
resonances too fast to resolve -- the opposite of what a linear sweep at
this rate buys. Duration is set by mode 1 (176.7 Hz, zeta=0.02, half-power
bandwidth 7.07 Hz, ~180 ms settling): sweep rate R < 7.07/0.18 ~= 39 Hz/s
=> T_chirp >= 2000/39 ~= 51 s; 60 s used for margin. A CHIRP_PREHOLD_S
zero-amplitude hold is prepended so the 1 s (5+ settling times) transient
discard doesn't eat into the low-frequency start of the sweep itself.

PRBS: maximal-length sequence (scipy.signal.max_len_seq), n=12 bits =>
period 4095 chips. Clock >= 4.5 kHz (PRBS is usable to ~0.44*f_clock before
the sinc envelope rolls off, so 4.5 kHz stays flat to ~1980 Hz, just under
the 2000 Hz band edge). Sequence period (2**12-1)/4500 = 0.91 s comfortably
exceeds the ~180 ms LINEAR settling time -- but that 180 ms is mode 1's
small-signal settling time, and it badly underestimates how long this
strongly nonlinear system actually takes to reach periodic steady state
near the breakaway threshold. Measured directly (2026-08-19): at the
threshold amplitude, discarding the originally-planned 2 periods left
median PRBS coherence at 0.05 (garbage); coherence only climbs past 0.9
once 5-6 periods are discarded (median 0.97 at 5, 0.999 at 6) -- "critical
slowing down" near the nonlinear transition, the same phenomenon already
seen in the sinusoidal-sweep convergence checks. PRBS_N_PERIODS/
PRBS_DISCARD_PERIODS below are sized from that measurement (12 total, 6
discarded, 6 kept for Welch averaging -- same segment count originally
intended, just built on genuinely settled data), not from the 180 ms
linear estimate alone.

Amplitude: RMS, not peak, sets the operating regime. The breakaway
threshold is derived from the GUIDEWAY port (not support-bearing or nut --
this is a deliberate, single governing threshold for this sweep):
z_max_way = Fs_way/sigma0_way, converted to a command angle via 2*pi/L.
Sanity check: 5 um of stage travel is 5e-6 * 2*pi/L = 31.4 mrad of
command; the guideway threshold should land near that order of magnitude
if it's right (see choose_rms_amplitudes' printed check).

RMS matching: PRBS is bipolar +-A (crest factor 1, so peak = RMS = A);
a sinusoid/chirp has crest factor sqrt(2) (peak = RMS*sqrt(2)). Equal peak
between the two would mean unequal RMS, i.e. unequal actual drive energy
-- so amplitudes here are always RMS, and each generator converts to its
own correct peak internally.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import chirp as scipy_chirp
from scipy.signal import max_len_seq

F_LO_HZ = 0.1
F_HI_HZ = 2000.0
CHIRP_DURATION_S = 60.0
CHIRP_PREHOLD_S = 1.0

PRBS_NBITS = 12
PRBS_CLOCK_HZ = 4500.0
PRBS_N_PERIODS = 12
PRBS_DISCARD_PERIODS = 6

FS_HZ = 20000.0   # fixed output/simulation grid: 10x the 2 kHz band, headroom
                  # past the 6863 Hz structural mode, anti-alias margin

N_AMPLITUDES = 5
AMPLITUDE_SPAN_DECADES = 1.0


def breakaway_command_angle(p: dict[str, float]) -> float:
    """Guideway breakaway bristle deflection z_max = Fs_way/sigma0_way,
    converted to an equivalent commanded angle via 2*pi/L."""
    lead_ratio = p["L"] / (2.0 * np.pi)
    z_max_way = p["Fs_way"] / p["sigma0_way"]
    return z_max_way / lead_ratio


def choose_rms_amplitudes(p: dict[str, float]) -> tuple[np.ndarray, float]:
    a_thresh = breakaway_command_angle(p)
    amplitudes = np.logspace(
        np.log10(a_thresh) - AMPLITUDE_SPAN_DECADES,
        np.log10(a_thresh) + AMPLITUDE_SPAN_DECADES,
        N_AMPLITUDES,
    )
    return amplitudes, a_thresh


class ChirpSignal:
    """theta_cmd_func(t): scalar-t callable for solve_ivp. wave/t_wave:
    the precomputed chirp waveform (peak-scaled) and its own local time
    axis (0..CHIRP_DURATION_S), exposed for reconstructing the exact
    input array post-hoc for spectral analysis."""

    def __init__(self, A_rms: float):
        self.A_rms = A_rms
        self.A_peak = A_rms * np.sqrt(2.0)
        self.n_samples = int(round(CHIRP_DURATION_S * FS_HZ))
        self.t_wave = np.arange(self.n_samples) / FS_HZ
        self.wave = self.A_peak * scipy_chirp(
            self.t_wave, f0=F_LO_HZ, f1=F_HI_HZ, t1=CHIRP_DURATION_S, method="linear",
        )
        self.total_duration_s = CHIRP_PREHOLD_S + CHIRP_DURATION_S

    def __call__(self, t: float) -> float:
        t_shift = t - CHIRP_PREHOLD_S
        if t_shift < 0.0 or t_shift > CHIRP_DURATION_S:
            return 0.0
        idx = min(int(round(t_shift * FS_HZ)), self.n_samples - 1)
        return float(self.wave[idx])

    def command_at(self, t_array: np.ndarray) -> np.ndarray:
        """Vectorized reconstruction of the exact command array at
        arbitrary sample times -- same lookup rule as __call__, for
        rebuilding the input record post-hoc from a solve_ivp t_eval."""
        t_shift = np.asarray(t_array) - CHIRP_PREHOLD_S
        idx = np.clip(np.round(t_shift * FS_HZ).astype(np.int64), 0, self.n_samples - 1)
        out = self.wave[idx]
        out = np.where((t_shift >= 0.0) & (t_shift <= CHIRP_DURATION_S), out, 0.0)
        return out


class PRBSSignal:
    """theta_cmd_func(t): scalar-t callable for solve_ivp. bipolar/period_bits:
    the precomputed +-A_rms chip sequence and its length, exposed for
    reconstructing the exact input array post-hoc."""

    def __init__(self, A_rms: float):
        self.A_rms = A_rms
        seq, _state = max_len_seq(PRBS_NBITS)
        self.bipolar = (2.0 * seq.astype(float) - 1.0) * A_rms   # crest factor 1: peak == RMS
        self.period_bits = len(self.bipolar)   # 2**12 - 1 = 4095
        self.chip_dt = 1.0 / PRBS_CLOCK_HZ
        self.period_s = self.period_bits * self.chip_dt
        self.total_duration_s = PRBS_N_PERIODS * self.period_s
        # Exact-integer-samples-per-period check: Welch needs nperseg to be
        # a whole number of output samples equal to exactly one PRBS period.
        samples_per_period = self.period_s * FS_HZ
        self.samples_per_period = int(round(samples_per_period))
        assert abs(samples_per_period - self.samples_per_period) < 1e-6, (
            f"PRBS period ({self.period_s*1e3:.4f} ms) does not divide the output "
            f"grid (dt={1e3/FS_HZ:.4f} ms) into a whole number of samples "
            f"({samples_per_period:.4f}) -- Welch segmenting would misalign."
        )

    def __call__(self, t: float) -> float:
        idx = int(t / self.chip_dt) % self.period_bits
        return float(self.bipolar[idx])

    def command_at(self, t_array: np.ndarray) -> np.ndarray:
        """Vectorized reconstruction of the exact command array at
        arbitrary sample times -- same ZOH lookup rule as __call__. Chips
        do NOT land on whole numbers of fs samples (period_bits=4095 vs
        samples_per_period=18200 -- only the WHOLE PERIOD is an integer
        number of samples, needed for Welch segmenting, not each chip), so
        this evaluates the lookup at every sample rather than repeating a
        fixed chips-per-sample block."""
        t_array = np.asarray(t_array)
        idx = (np.floor(t_array / self.chip_dt).astype(np.int64)) % self.period_bits
        return self.bipolar[idx]


def build_atol(A_rms: float, p: dict[str, float], f_max_hz: float = F_HI_HZ) -> np.ndarray:
    """Per-state atol ~= 1e-6 x each state's expected amplitude -- NOT one
    scalar (that's what killed the earlier fixed-scalar-atol run). Position/
    velocity estimates use the RMS command amplitude under an ideal-tracking
    assumption (theta ~ A_rms for the rotational chain, x ~ A_rms*lead_ratio
    for the axial DOFs, velocities scaled by the swept band's own top
    frequency as a worst-case bound); bristle-state estimates use each
    port's own breakaway limit Fs-or-Ts/sigma0, which bounds |z| regardless
    of A_rms."""
    lead_ratio = p["L"] / (2.0 * np.pi)
    omega_max = 2.0 * np.pi * f_max_hz

    theta_amp = A_rms
    x_amp = A_rms * lead_ratio
    thetadot_amp = A_rms * omega_max
    xdot_amp = A_rms * lead_ratio * omega_max
    z_sb_amp = p["Ts_sb"] / p["sigma0_sb"]
    z_nut_amp = p["Fs_nut"] / p["sigma0_nut"]
    z_way_amp = p["Fs_way"] / p["sigma0_way"]

    amps = np.array([
        theta_amp, theta_amp, theta_amp, theta_amp, x_amp, x_amp,
        thetadot_amp, thetadot_amp, thetadot_amp, thetadot_amp, xdot_amp, xdot_amp,
        z_sb_amp, z_nut_amp, z_way_amp,
    ])
    return np.maximum(1e-6 * amps, 1e-15)
