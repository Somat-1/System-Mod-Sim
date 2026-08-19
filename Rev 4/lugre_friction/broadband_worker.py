#!/usr/bin/env python3
"""Per-run worker for the parallel broadband ID sweep (run_broadband_id.py).

Each of the (amplitude, signal_type) combinations -- 5 amplitudes x
{chirp, prbs} -- is a fully independent simulation, so they're dispatched
across a ProcessPoolExecutor. Same pattern as the earlier sinusoidal
sweep's worker: deliberately no numpy/scipy import at module level, so the
pool's initializer can set single-threaded BLAS env vars before each
worker's first numpy import, rather than after (oversubscribed multi-
threaded BLAS inside every worker would fight the pool for cores instead
of adding throughput).
"""

from __future__ import annotations

_model = None


def init_worker() -> None:
    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    from lugre_model_broadband import LuGreModelNonlinearDrive

    global _model
    _model = LuGreModelNonlinearDrive()


def run_one(payload: dict) -> dict:
    """payload: {'A': float, 'signal_type': 'prbs' or 'chirp'}."""
    import time

    import numpy as np
    from scipy.integrate import solve_ivp

    from broadband_estimators import chirp_estimate, fractional_octave_smooth, prbs_estimate
    from broadband_signals import (
        ChirpSignal,
        FS_HZ,
        PRBSSignal,
        PRBS_DISCARD_PERIODS,
        build_atol,
    )
    from lugre_model import N_STATES

    global _model
    model = _model
    p = model.p
    A, signal_type = payload["A"], payload["signal_type"]
    atol = build_atol(A, p)

    def simulate(theta_cmd_func, total_duration_s):
        n_samples = int(round(total_duration_s * FS_HZ)) + 1
        t_eval = np.arange(n_samples) / FS_HZ
        x0 = np.zeros(N_STATES)
        sol = solve_ivp(
            model.rhs, (0.0, t_eval[-1]), x0, method="LSODA",
            args=(theta_cmd_func,), t_eval=t_eval, rtol=1e-6, atol=atol,
        )
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed (A={A}, {signal_type}): {sol.message}")
        return sol.t, sol.y

    def saturation(y):
        z_sb, z_nut, z_way = y[12], y[13], y[14]
        z_max_sb = p["Ts_sb"] / p["sigma0_sb"]
        z_max_nut = p["Fs_nut"] / p["sigma0_nut"]
        z_max_way = p["Fs_way"] / p["sigma0_way"]
        return {
            "sb": float(np.max(np.abs(z_sb)) / z_max_sb),
            "nut": float(np.max(np.abs(z_nut)) / z_max_nut),
            "way": float(np.max(np.abs(z_way)) / z_max_way),
        }

    t0 = time.time()
    if signal_type == "prbs":
        sig = PRBSSignal(A)
        t, y = simulate(sig, sig.total_duration_s)
        u_full = sig.command_at(t)
        y_full = y[5]
        f, G, gamma2, n_segments = prbs_estimate(
            u_full, y_full, FS_HZ, sig.samples_per_period, PRBS_DISCARD_PERIODS,
        )
        result = dict(f=f, G=G, gamma2=gamma2, n_segments=n_segments)
    elif signal_type == "chirp":
        sig = ChirpSignal(A)
        t, y = simulate(sig, sig.total_duration_s)
        prehold_samples = int(round(1.0 * FS_HZ))
        u_full = sig.command_at(t)
        y_full = y[5]
        u_trim, y_trim = u_full[prehold_samples:], y_full[prehold_samples:]
        f, G = chirp_estimate(u_trim, y_trim, FS_HZ)
        f_s, mag_s, phase_s = fractional_octave_smooth(f, G, frac=24.0)
        result = dict(f=f_s, mag_db=mag_s, phase_deg=phase_s)
    else:
        raise ValueError(f"unknown signal_type {signal_type!r}")

    result.update(
        A=A, signal_type=signal_type, saturation=saturation(y),
        wall_s=time.time() - t0, n_samples=len(t),
    )
    return result
