#!/usr/bin/env python3
"""Empirical sinusoidal-sweep describing-function analysis of the 15-state
LuGre model (lugre_model.py), per the requested protocol:

  theta_cmd(t) = A*sin(omega_i * t), with omega_i drawn from a DENSE LINEAR
  grid in Hertz, 0-2000 Hz (FREQ_HZ), scaled to rad/s only for the ODE loop
  (OMEGA = FREQ_HZ * 2*pi). The very first grid point is set to 0.1 Hz
  rather than exactly 0 Hz to avoid a divide-by-zero on the period
  T = 2*pi/omega. All reporting/plotting is in Hz.
  Integrate the full nonlinear ODE (scipy.integrate.solve_ivp, stiff
  method) for N_PERIODS cycles, then measure the steady-state fundamental
  response of x_n(t) per three post-processing rules (2026-08-18 revision,
  replacing the original zero-crossing-based extraction):
    Rule 1: keep only the last 50% of the simulated window (discard
            transients from the first half, not just N_TRANSIENT_PERIODS).
    Rule 2: if the steady-state peak-to-peak of x_n is below a realistic
            sensor floor, phase is forced to exactly 0 deg.
    Rule 3: magnitude and phase both come from a single-frequency DFT at
            the exact known driving frequency (correlate the steady window
            against exp(-j*omega*t)), not zero-crossing timing:
              U = sum(theta_cmd(t)*exp(-j*omega*t))
              Y = sum(x_n(t)*exp(-j*omega*t))
              G = Y/U
              Magnitude (dB) = 20*log10(abs(G))
              Phase (deg)    = unwrap(angle(G)) * 180/pi  (unwrapped across
                                the whole frequency sweep, then Rule 2 is
                                applied on top)

Run twice, at a small amplitude (stays trapped in presliding/stiction
across most of the sweep) and a large amplitude (breaks free into gross
sliding across most of the sweep) -- the nonlinearity means these are NOT
the same curve.

Solver note (2026-08-18): Radau and BDF both diverge to NaN on this model
regardless of rtol/atol/max_step tuning, even after fixing an unrelated
parameter bug (sigma0_sb was 3000x too stiff for I_sb, confirmed via these
same solver tests -- see model_parameters.json parameter_notes). The
remaining failure is consistent with LuGre's non-smooth kink at v=0
defeating their Newton-iteration/finite-difference-Jacobian machinery.
LSODA, which switches to cheap explicit steps whenever the local problem
isn't actually stiff, integrates the same cases successfully. LSODA is
therefore the working default here; method= is still a plain solve_ivp
argument if Radau/BDF are revisited with an analytical Jacobian later.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from lugre_model import LuGreModel, N_STATES

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"

METHOD = "LSODA"           # see module docstring: Radau/BDF diverge on this problem

# Dense LINEAR grid in Hz, 0-2000 Hz. N_FREQ=1001 gives exactly 2 Hz steps
# (the coarser of the two options given); set N_FREQ=2001 for 1 Hz steps.
N_FREQ = 1001
FREQ_HZ = np.linspace(0.0, 2000.0, N_FREQ)
FREQ_HZ[0] = 0.1   # protect the 0 Hz boundary: T = 2*pi/omega would divide by zero
OMEGA = FREQ_HZ * 2.0 * np.pi   # rad/s, ODE loop only -- all reporting stays in Hz

N_PERIODS = 6           # Rule 1 keeps the last 50% of this -- final 3 cycles
POINTS_PER_PERIOD = 100

# Per-state absolute tolerances: the 15 states span meters (~1e-7-1e-4),
# radians (~1e-6-1e-1), their time derivatives, and bristle deflections on
# the same scales as their host DOF. A single scalar atol is either too
# loose for the position states or absurdly tight for the velocity states;
# both caused solver instability during testing.
ATOL = np.array([
    1e-10, 1e-10, 1e-10, 1e-10, 1e-13, 1e-13,      # q: theta x4, x_s, x_n
    1e-8, 1e-8, 1e-8, 1e-8, 1e-10, 1e-10,           # qdot
    1e-10, 1e-13, 1e-13,                             # z_sb, z_nut, z_way
])
RTOL = 1e-6

# Amplitudes (rad). Sized against the nut interface's own presliding range
# (Fs_nut/sigma0_nut = 0.75/2e6 = 3.75e-7 m) and Stribeck velocity
# (vs_nut = 2e-4 m/s): A_SMALL keeps the peak nut velocity (L/2pi)*omega*A
# below vs_nut across the whole sweep (stays in presliding); A_LARGE pushes
# it well past vs_nut for most of the sweep (breaks into sliding). See
# README.md for the arithmetic.
A_SMALL = 5.0e-6   # rad
A_LARGE = 5.0e-3   # rad

# Post-processing rules (2026-08-18 revision -- replaces the zero-crossing
# based extraction, which was noise-sensitive and produced a spurious
# low-frequency magnitude spike):
#   Rule 1: only the LAST 50% of the simulated window (which, at
#           N_PERIODS=6, leaves the final 3 cycles) is passed to the
#           magnitude/phase measurement -- not just N_TRANSIENT_PERIODS/
#           N_PERIODS as before.
#   Rule 2: if the steady-state peak-to-peak of x_n falls below a realistic
#           sensor floor, phase is forced to exactly 0 deg rather than
#           reporting a meaningless number extracted from numerical noise.
#   Rule 3: magnitude AND phase both come from a single-frequency DFT
#           (correlate the steady-state window against exp(-j*omega*t)) at
#           the exact known driving frequency, not zero-crossing timing.
TRIM_FRACTION = 0.5
AMPLITUDE_THRESHOLD = 1.0e-9  # m, peak-to-peak


def run_single(model: LuGreModel, A: float, omega: float, method: str = METHOD):
    """Integrate one (A, omega) case; return (mag_db, phase_deg, t, theta_cmd, x_n)
    where the trace arrays cover the full run (transient + steady state)."""
    T = 2.0 * np.pi / omega

    def theta_cmd_func(t):
        return A * np.sin(omega * t)

    t_span = (0.0, N_PERIODS * T)
    t_eval = np.linspace(0.0, N_PERIODS * T, N_PERIODS * POINTS_PER_PERIOD)
    x0 = np.zeros(N_STATES)

    sol = solve_ivp(
        model.rhs, t_span, x0, method=method, args=(theta_cmd_func,),
        t_eval=t_eval, rtol=RTOL, atol=ATOL,
    )
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed at A={A}, omega={omega}: {sol.message}")

    # Rule 1: discard the first 50% of the total simulated window.
    t_trim_start = TRIM_FRACTION * (N_PERIODS * T)
    mask = sol.t >= t_trim_start
    t_ss = sol.t[mask]
    theta_cmd_ss = A * np.sin(omega * t_ss)
    x_n_ss = sol.y[5, mask]

    peak_to_peak = float(x_n_ss.max() - x_n_ss.min())

    # Rule 3: single-frequency DFT at the exact driving frequency.
    kernel = np.exp(-1j * omega * t_ss)
    U = np.sum(theta_cmd_ss * kernel)
    Y = np.sum(x_n_ss * kernel)
    G = Y / U
    mag_db = 20.0 * np.log10(max(abs(G), 1e-300))
    phase_raw_deg = float(np.angle(G)) * (180.0 / np.pi)  # NOT unwrapped yet --
    # run_sweep unwraps across the whole frequency sweep first, THEN applies
    # the Rule 2 override below it, so a handful of forced-0 points can't
    # corrupt np.unwrap's continuity for their genuinely-computed neighbors.

    return mag_db, phase_raw_deg, peak_to_peak, sol.t, A * np.sin(omega * sol.t), sol.y[5, :]


def run_sweep(model: LuGreModel, A: float, label: str):
    mags = np.empty(N_FREQ)
    phases_raw = np.empty(N_FREQ)
    below_threshold = np.zeros(N_FREQ, dtype=bool)
    example = None
    example_target_hz = 15.9  # ~100 rad/s -- shared example frequency for the montage
    example_idx = int(np.argmin(np.abs(FREQ_HZ - example_target_hz)))
    for i, omega in enumerate(OMEGA):
        t0 = time.time()
        mag_db, phase_raw_deg, peak_to_peak, t, theta_cmd, x_n = run_single(model, A, omega)
        mags[i] = mag_db
        phases_raw[i] = phase_raw_deg
        below_threshold[i] = peak_to_peak < AMPLITUDE_THRESHOLD
        if i == example_idx:
            example = (FREQ_HZ[i], t, theta_cmd, x_n)
        print(f"  [{label}] {i+1:4d}/{N_FREQ} f={FREQ_HZ[i]:8.2f} Hz  mag={mag_db:8.2f} dB  "
              f"phase_raw={phase_raw_deg:8.2f} deg  p2p={peak_to_peak:.3e} m  "
              f"({time.time()-t0:.2f}s)", flush=True)

    # Unwrap across the whole sweep first (Rule 3), THEN pin the low-amplitude
    # points to exactly 0 deg (Rule 2) -- doing it in this order keeps the
    # forced zeros from disturbing np.unwrap's continuity for their neighbors.
    phases = np.unwrap(phases_raw, period=360.0)
    phases[below_threshold] = 0.0
    n_pinned = int(below_threshold.sum())
    if n_pinned:
        print(f"  [{label}] {n_pinned}/{N_FREQ} points pinned to 0 deg "
              f"(peak-to-peak x_n < {AMPLITUDE_THRESHOLD:.0e} m)")
    return mags, phases, example


def main():
    model = LuGreModel()
    print(f"Sweep: {N_FREQ} frequencies, {FREQ_HZ[0]:.2f}-{FREQ_HZ[-1]:.2f} Hz, method={METHOD}")

    print("Running A_SMALL sweep...")
    mags_small, phases_small, example_small = run_sweep(model, A_SMALL, "small")

    print("Running A_LARGE sweep...")
    mags_large, phases_large, example_large = run_sweep(model, A_LARGE, "large")

    ASSET_DIR.mkdir(exist_ok=True)
    np.savez(
        ASSET_DIR / "sinusoidal_sweep_data.npz",
        freq_hz=FREQ_HZ, omega=OMEGA, mags_small=mags_small, phases_small=phases_small,
        mags_large=mags_large, phases_large=phases_large,
        A_small=A_SMALL, A_large=A_LARGE,
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.5))
    ax_mag, ax_phase, ax_ex_small, ax_ex_large = axes.flat

    ax_mag.plot(FREQ_HZ, mags_small, color="#2b6cb0", linewidth=1.0,
                label=f"A = {A_SMALL:.0e} rad (small)")
    ax_mag.plot(FREQ_HZ, mags_large, color="#c05621", linewidth=1.0,
                label=f"A = {A_LARGE:.0e} rad (large)")
    ax_mag.set_xlabel("Frequency (Hz)")
    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title("x_n / theta_cmd magnitude")
    ax_mag.grid(True, linewidth=0.4, color="#cccccc")
    ax_mag.legend(fontsize=8)

    ax_phase.plot(FREQ_HZ, phases_small, color="#2b6cb0", linewidth=1.0)
    ax_phase.plot(FREQ_HZ, phases_large, color="#c05621", linewidth=1.0)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.set_title("x_n / theta_cmd phase")
    ax_phase.grid(True, linewidth=0.4, color="#cccccc")

    for ax, example, label, color in [
        (ax_ex_small, example_small, f"A={A_SMALL:.0e} rad (stiction-trapped)", "#2b6cb0"),
        (ax_ex_large, example_large, f"A={A_LARGE:.0e} rad (sliding)", "#c05621"),
    ]:
        f_ex, t, theta_cmd, x_n = example
        lead_ratio = model.lead_ratio
        ax.plot(t * 1e3, theta_cmd * lead_ratio * 1e6, color="#999999", linewidth=1.0,
                linestyle="--", label="theta_cmd (ideal, scaled)")
        ax.plot(t * 1e3, x_n * 1e6, color=color, linewidth=1.1, label="x_n (actual)")
        ax.set_title(f"Example @ f={f_ex:.1f} Hz, {label}")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Position (um)")
        ax.grid(True, linewidth=0.4, color="#cccccc")
        ax.legend(fontsize=7)

    fig.suptitle("Rev 4 LuGre sub-branch: sinusoidal describing-function sweep")
    fig.tight_layout()
    montage_path = ASSET_DIR / "sinusoidal_sweep_montage.svg"
    fig.savefig(montage_path)
    fig.savefig(montage_path.with_suffix(".png"), dpi=110)
    plt.close(fig)
    print(f"Wrote {montage_path}")


if __name__ == "__main__":
    main()
