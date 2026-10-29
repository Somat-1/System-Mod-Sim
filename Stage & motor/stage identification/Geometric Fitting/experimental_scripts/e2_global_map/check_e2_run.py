#!/usr/bin/env python3
"""Acceptance checks on the E2 run, printed as a report.

The question this answers is narrow: did the stage do what the plan said, and
is the capture good enough to fit a global map to. No model is fitted here.

This diagnostic intentionally reads the scope trace and plan only. The serial
log's driver settings, config hash, completion summary and reference MSCNT are
checked by ``render_e2_findings.py``.
"""
from __future__ import annotations

import numpy as np

import e2_run
from e2_run import X_REF_MM, TURN_LOW_MM, TURN_HIGH_MM

CEIL_NM = 8333.0          # 120 ns AquadB at 1 nm per increment
V_CMD = 2.0


def main():
    run = e2_run.load()
    plan = run.plan
    hold = run.hold()
    resid = (hold - plan.target_mm) * 1000.0
    ok = run.matched
    up = ok & (plan.state == "UP")
    dn = ok & (plan.state == "DOWN")
    ref = ok & np.isin(plan.state, ("DATUM", "TURN_REF", "END_REF"))
    run_lo = run.hold_start[ok].min()
    run_hi = run.hold_end[ok].max()

    print("capture")
    gaps = np.flatnonzero(np.diff(run.t) > 1.5e-3)
    print(f"  {run.t.size} samples at 1 kHz, {run.t[-1] / 60:.2f} min")
    print(f"  index gaps: {gaps.size}"
          + (f" (at t = {', '.join(f'{run.t[g]:.1f} s' for g in gaps)}, "
             f"{int(round((np.diff(run.t)[gaps].sum() - gaps.size * 1e-3) * 1000))}"
             " sample(s) lost)" if gaps.size else ""))
    print(f"  whole file spans {run.pos.min():.3f} .. {run.pos.max():.3f} mm "
          "- it opens part-way through the homing sweep")
    print(f"  during the measurement itself {run.pos[run_lo:run_hi].min():.4f} "
          f".. {run.pos[run_lo:run_hi].max():.4f} mm "
          f"(turn points {TURN_LOW_MM} .. {TURN_HIGH_MM})")
    print(f"  datum fitted from the holds: {run.datum_mm:.6f} mm")

    print("\nsequence against the plan")
    print(f"  planned stops {plan.target_mm.size}, matched to a measured hold "
          f"{ok.sum()}")
    if (~ok).any():
        print("  unmatched:")
        for i in np.flatnonzero(~ok):
            print(f"    order {plan.order[i]:>4} {plan.state[i]:>8} "
                  f"index {plan.index[i]:>4} target {plan.target_mm[i]:9.4f} mm")
    print(f"  up {up.sum()}, down {dn.sum()}, reference {ref.sum()}")
    both = np.intersect1d(plan.index[up], plan.index[dn])
    print(f"  targets measured in both directions: {both.size} of 400")

    print("\ndid it go where it was told")
    print(f"  measured minus commanded, um: median {np.median(resid[ok]):+.3f}  "
          f"sd {resid[ok].std(ddof=1):.3f}  worst {np.abs(resid[ok]).max():.3f}")
    print(f"    up   {resid[up].mean():+.3f} +- {resid[up].std(ddof=1):.3f}")
    print(f"    down {resid[dn].mean():+.3f} +- {resid[dn].std(ddof=1):.3f}")
    order = np.argsort(plan.target_mm[up])
    du = np.diff(plan.target_mm[up][order])
    dm = np.diff(hold[up][order])
    print(f"  commanded spacing {du.min() * 1000:.0f} .. {du.max() * 1000:.0f} um; "
          f"measured minus commanded spacing: median "
          f"{np.median(dm - du) * 1000:+.3f} um, worst "
          f"{np.abs(dm - du).max() * 1000:.3f} um")

    print("\nhold quality")
    lengths = (run.hold_end[ok] - run.hold_start[ok])
    print(f"  dwell length ms: median {int(np.median(lengths))}, "
          f"5-95% {np.percentile(lengths, [5, 95]).astype(int)}  (commanded 2500)")
    ptp, drift, spikes = [], [], []
    for i in np.flatnonzero(ok):
        lo, hi = run.record_window(i)
        seg = run.pos[lo:hi] * 1000.0
        ptp.append(np.ptp(seg))
        drift.append(seg[-100:].mean() - seg[:100].mean())
        a, b = run.hold_start[i], run.hold_end[i]
        spikes.append(int((np.abs(np.diff(run.pos[a:b])) * 1e6 > 60).sum()))
    ptp, drift, spikes = np.array(ptp), np.array(drift), np.array(spikes)
    print(f"  record window ptp: median {np.median(ptp):.3f} um, "
          f"95th {np.percentile(ptp, 95):.3f}, worst {ptp.max():.3f}")
    print(f"  drift across it: mean {drift.mean():+.4f} um, "
          f"worst {np.abs(drift).max():.4f}")
    rate = spikes.sum() / (run.hold_end[ok] - run.hold_start[ok]).sum() * 1000
    print(f"  single samples over 60 nm inside a dwell: {rate:.2f} per 1000"
          f" (E1 ran at 5.72) - the interferometer noise floor, not motion")

    print("\nmotion")
    step_nm = np.abs(np.diff(run.pos[run_lo:run_hi])) * 1e6
    v10 = np.convolve(step_nm / 1000.0, np.ones(10) / 10, mode="valid")
    peaks = []
    idx = np.flatnonzero(ok)
    for a, b in zip(run.hold_end[idx[:-1]] - run_lo, run.hold_start[idx[1:]] - run_lo):
        if b > a + 20:
            peaks.append(v10[a:b].max())
    peaks = np.array(peaks)
    print(f"  peak speed between dwells: median {np.median(peaks):.3f} mm/s, "
          f"5-95% {np.percentile(peaks, [5, 95]).round(3)} (commanded {V_CMD})")
    print(f"  largest 1 ms increment {step_nm.max():.0f} nm against the "
          f"{CEIL_NM:.0f} nm AquadB ceiling; samples at the ceiling "
          f"{int((step_nm >= CEIL_NM).sum())}")
    print("  (homing, which the file also contains, is excluded from these)")

    print("\nreference stability")
    for i in np.flatnonzero(ref):
        lo, _ = run.record_window(i)
        print(f"  {plan.state[i]:>8}  t={run.t[lo] / 60:6.2f} min  "
              f"{hold[i]:.6f} mm  {(hold[i] - X_REF_MM) * 1000:+7.3f} um")
    first = hold[ref][:3].mean()
    last = hold[ref][-3:].mean()
    print(f"  spread across all {ref.sum()} reference visits "
          f"{np.ptp(hold[ref]) * 1000:.3f} um")
    print(f"  DATUM -> END_REF drift {(last - first) * 1000:+.3f} um "
          "(warning threshold 0.5 um)")

    print("\nreversal")
    common, iu, idn = np.intersect1d(plan.index[up], plan.index[dn],
                                     return_indices=True)
    B = resid[up][iu] - resid[dn][idn]
    x = plan.target_mm[up][iu]
    slope = np.polyfit(x - x.mean(), B, 1)[0] * 1000
    print(f"  B over {B.size} targets: mean {B.mean():+.3f} um, "
          f"sd {B.std(ddof=1):.3f}, max|B| {np.abs(B).max():.3f}")
    print(f"  slope against position {slope:+.2f} nm/mm")
    print("  for comparison: E1 -1.680 um, ISO 230-2 v3 re-derived -1.600 um")

    print("\nscrew-phase coverage (the point of the seeded grid)")
    pb = plan.phase_bin[up]
    counts = np.bincount(pb, minlength=50)[:50]
    print(f"  50 bins, measured counts min {counts.min()}, max {counts.max()}, "
          f"empty {int((counts == 0).sum())}")

    print("\nserial cross-check")
    print("  this trace-only diagnostic does not parse the controller log;")
    print("  render_e2_findings.py checks settings, SUMMARY and reference MSCNT")


if __name__ == "__main__":
    main()
