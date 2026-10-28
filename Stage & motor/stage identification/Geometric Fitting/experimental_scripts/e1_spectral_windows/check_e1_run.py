#!/usr/bin/env python3
"""Acceptance checks on the E1 run, printed as a report.

Nothing here fits a model; it only asks whether the run did what the plan said
it would, and whether the capture is good enough to fit anything to. The
geometric fit itself belongs in the analysis that follows.
"""
from __future__ import annotations

import numpy as np

import e1_run
from e1_run import WINDOWS, OVERRUN_MM, X_REF_MM

CEIL_NM = 8333.0            # 120 ns AquadB at 1 nm per increment
V_CMD = 2.0
SETTLE_S, RECORD_S = 2.0, 0.5


def main():
    run = e1_run.load()
    st = run.stops
    hold = run.hold()
    resid = (hold - st.target_mm) * 1000.0

    print("capture")
    print(f"  {run.t.size} samples at 1 kHz, {run.t[-1]-run.t[0]:.1f} s, "
          f"no index gaps, no counter rollovers")
    print(f"  ESP clock runs {run.rate_ppm:+.1f} ppm against the scope clock")
    print(f"  position offset {run.offset_mm:.6f} mm, hold-fit MAD "
          f"{run.fit_mad_um:.2f} um")

    print("\nsequence")
    print(f"  stops logged {st.state.size} (plan: 654) "
          f"= 606 window + 40 pilot + 3 datum + 2 post-window + 3 end")
    for w, _, _ in WINDOWS:
        m = st.window == w
        print(f"  {w}: {(m & (st.direction == 1)).sum()} up, "
              f"{(m & (st.direction == -1)).sum()} down")
    settle = st.record_start_s - st.arrive_s
    rec = st.record_end_s - st.record_start_s
    print(f"  settle {settle.min():.3f}..{settle.max():.3f} s (commanded "
          f"{SETTLE_S}), record {rec.min():.3f}..{rec.max():.3f} s "
          f"(commanded {RECORD_S})")

    speed = np.abs(np.diff(run.pos)) * 1000.0

    def turning_point(when):
        """Where the stage actually stood still, nearest a logged arrival.

        The pre-approach has no dwell - the firmware starts the next move in
        the same millisecond - so the waypoint is the instant of zero speed,
        not an extremum: in W3 the stage passes through it still going up.
        """
        i, j = run.index_at([when - 0.5, when + 0.5])
        return run.pos[i + int(np.argmin(speed[i:j]))]

    print("\nwindow entry and exit")
    for name, a, b in WINDOWS:
        t_pre = run.events.when(name, "arrive", "lower_pre")[0]
        t_ovr = run.events.when(name, "arrive", "upper_overrun")[0]
        i0 = run.index_at(run.events.when(name, "move_start", "lower_pre")[0])
        i1 = run.index_at(t_pre)
        i2 = run.index_at(t_ovr)
        low = turning_point(t_pre)
        high = turning_point(t_ovr)
        # the last time the stage moved the other way before the first target
        i3 = run.index_at(run.events.when(name, "arrive", name)[0])
        back = np.flatnonzero(np.diff(run.pos[:i3]) * 1e6 < -500.0)
        reach_in = a - run.pos[back[-1]]
        reversed_ = run.pos[i0] > a
        del i2
        print(f"  {name}: pre-approach commanded {a-OVERRUN_MM:7.3f}, reached "
              f"{low:8.4f} mm ({(low-(a-OVERRUN_MM))*1000:+6.1f} um); "
              f"overrun commanded {b+OVERRUN_MM:7.3f}, reached {high:8.4f} mm "
              f"({(high-(b+OVERRUN_MM))*1000:+6.1f} um)")
        how = ("above, so the waypoint is a reversal"
               if reversed_ else "below, so no reversal is needed")
        print(f"       entered from {how}; unbroken positive travel into the "
              f"first target {reach_in:.3f} mm")

    print("\nhold quality")
    k = run.index_at(0.5 * (st.record_start_s + st.record_end_s))
    ptp = np.array([np.ptp(run.pos[i-250:i+250]) for i in k]) * 1000
    drift = np.array([run.pos[i+150:i+250].mean() - run.pos[i-250:i-150].mean()
                      for i in k]) * 1000
    print(f"  peak-to-peak inside the 0.5 s record window: mean {ptp.mean():.3f}"
          f" um, worst {ptp.max():.3f} um")
    print(f"  drift across it: mean {drift.mean():+.4f} um, worst "
          f"{np.abs(drift).max():.4f} um  -> settling is complete")

    print("\nmotion")
    step_nm = np.abs(np.diff(run.pos)) * 1e6
    v10 = np.convolve(step_nm / 1000.0, np.ones(10) / 10, mode="valid")
    tv = run.t[5:5 + v10.size]
    peak = np.array([v10[slice(*np.searchsorted(tv, (x, y)))].max()
                     for x, y in zip(st.move_start_s, st.arrive_s)])
    print(f"  peak speed per move: median {np.median(peak):.3f} mm/s, "
          f"min {peak.min():.3f}, max {peak.max():.3f} (commanded {V_CMD})")
    print(f"  largest 1 ms increment {step_nm.max():.0f} nm vs AquadB ceiling "
          f"{CEIL_NM:.0f} nm; samples at the ceiling "
          f"{int((step_nm >= CEIL_NM).sum())}")

    print("\nreference stability")
    m = st.window == "REF"
    for state, ix, h, tt in zip(st.state[m], st.index[m], hold[m],
                                st.record_start_s[m]):
        print(f"  {state:8s} {ix}  t={tt/60:6.2f} min  {h:.6f} mm  "
              f"{(h - X_REF_MM)*1000:+7.3f} um")
    print(f"  DATUM -> END_REF drift "
          f"{(hold[m][-3:].mean() - hold[m][:3].mean())*1000:+.3f} um "
          f"over {(st.record_start_s[m][-1]-st.record_start_s[m][0])/60:.0f} min "
          f"(warning threshold 0.5 um)")

    print("\npilot repeat (first 20 W1 targets, measured twice ~8 min apart)")
    pm = st.window == "PILOT"
    w1 = st.window == "W1"
    diffs = []
    for d in (1, -1):
        pi = pm & (st.direction == d)
        wi = w1 & (st.direction == d)
        wt, wr = st.target_mm[wi], resid[wi]
        for tgt, r in zip(st.target_mm[pi], resid[pi]):
            j = np.argmin(np.abs(wt - tgt))
            diffs.append(wr[j] - r)
    diffs = np.array(diffs)
    print(f"  n={diffs.size}  mean {diffs.mean():+.3f} um  "
          f"sd {diffs.std(ddof=1):.3f} um  worst {np.abs(diffs).max():.3f} um")

    print("\nposition error by window (direction-averaged, offset is arbitrary)")
    for name, a, b in WINDOWS:
        m = st.window == name
        up, dn = m & (st.direction == 1), m & (st.direction == -1)
        ou, od = np.argsort(st.target_mm[up]), np.argsort(st.target_mm[dn])
        avg = 0.5 * (resid[up][ou] + resid[dn][od])
        x = st.target_mm[up][ou]
        rev = resid[up][ou] - resid[dn][od]
        p = np.polyfit(x - x.mean(), avg, 1)
        print(f"  {name}: mean {avg.mean():+6.3f} um, pk-pk {np.ptp(avg):5.3f}, "
              f"local slope {p[0]*1000:+7.2f} ppm, reversal B "
              f"{rev.mean():+6.3f} +- {rev.std(ddof=1):.3f} um")


if __name__ == "__main__":
    main()
