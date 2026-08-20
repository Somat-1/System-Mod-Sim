#!/usr/bin/env python3
"""Direct sliding-vs-presliding check for the smallest broadband-ID
amplitude (A_rms=3.7203 mrad), requested 2026-08-20. The earlier saturation
ratio (peak |z|/(Fs-or-Ts/sigma0)) is a bristle-deflection proxy; the
textbook LuGre criterion for gross sliding is velocity-based: |v| > vs
(Stribeck velocity) means the port has left the Stribeck-blended regime and
is sliding at a rate where g(v) ~= Fc, not presliding near v=0 where
g(v) ~= Fs. This reruns the PRBS case (same settings as
run_broadband_id.py: 12 periods, atol from build_atol, rest start) because
the production sweep never saved the raw velocity time series, only
summary stats.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.integrate import solve_ivp

from broadband_signals import FS_HZ, PRBSSignal, PRBS_DISCARD_PERIODS, build_atol
from lugre_model_broadband import LuGreModelNonlinearDrive
from lugre_model import N_STATES

A_RMS = 3.72030709e-3   # smallest of the 5 broadband-ID amplitudes


def main() -> None:
    model = LuGreModelNonlinearDrive()
    p = model.p
    lead_ratio = model.lead_ratio

    sig = PRBSSignal(A_RMS)
    atol = build_atol(A_RMS, p)
    n = int(round(sig.total_duration_s * FS_HZ)) + 1
    t_eval = np.arange(n) / FS_HZ
    x0 = np.zeros(N_STATES)

    print(f"A_rms = {A_RMS*1e3:.4f} mrad, PRBS total duration = {sig.total_duration_s:.2f} s "
          f"({sig.total_duration_s*FS_HZ:.0f} samples)")
    t0 = time.time()
    sol = solve_ivp(model.rhs, (0.0, t_eval[-1]), x0, method="LSODA",
                     args=(sig,), t_eval=t_eval, rtol=1e-6, atol=atol)
    print(f"solve_ivp done in {time.time()-t0:.1f}s, success={sol.success}")
    if not sol.success:
        raise RuntimeError(sol.message)
    y = sol.y

    start = PRBS_DISCARD_PERIODS * sig.samples_per_period
    theta_s_dot, theta_sb_dot = y[8, start:], y[9, start:]
    x_s_dot, x_n_dot = y[10, start:], y[11, start:]
    v_sb = theta_sb_dot
    v_way = x_n_dot
    v_nut = lead_ratio * theta_s_dot - (x_n_dot - x_s_dot)

    print(f"\nSteady-state window: {len(v_sb)} samples "
          f"({PRBS_DISCARD_PERIODS} of {sig.samples_per_period*12} periods discarded)")
    print(f"\n{'port':6s} {'vs (Stribeck)':>16s} {'peak|v|':>14s} {'RMS|v|':>14s} "
          f"{'peak|v|/vs':>12s} {'time sliding (|v|>vs)':>24s}")
    for name, v, vs in [("sb", v_sb, p["vs_sb"]), ("nut", v_nut, p["vs_nut"]), ("way", v_way, p["vs_way"])]:
        peak_v = float(np.max(np.abs(v)))
        rms_v = float(np.sqrt(np.mean(v ** 2)))
        frac_sliding = float(np.mean(np.abs(v) > vs))
        unit = "rad/s" if name == "sb" else "m/s"
        print(f"{name:6s} {vs:16.4e} {peak_v:14.4e} {rms_v:14.4e} "
              f"{peak_v/vs:12.4f} {frac_sliding*100:22.3f} %   [{unit}]")
        verdict = "SLIDING at times" if frac_sliding > 0.0 else "never exceeds vs -- presliding throughout"
        print(f"       -> {verdict}")


if __name__ == "__main__":
    main()
