#!/usr/bin/env python3
"""Inspect the open-loop trajectory diagnostic from the raw IFM counter."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks

from render_raw_measurement import read_scope


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE.parent.parent / "IsoTest.csv"
DEFAULT_OUTPUT = HERE / "plots" / "iso_open_loop_diagnostic_overview.png"
PROFILE_NAMES = ("current", "low_accel", "low_speed", "conservative")
PROFILE_LABELS = ("current: 5 mm/s, 40 mm/s²",
                  "low acceleration: 5 mm/s, 10 mm/s²",
                  "low speed: 2 mm/s, 40 mm/s²",
                  "conservative: 1 mm/s, 5 mm/s²")
EXPECTED_SPAN_MM = 59.35


def block_mean(a: np.ndarray, n: int) -> np.ndarray:
    end = len(a) - len(a) % n
    return a[:end].reshape(-1, n).mean(axis=1)


def motion_segments(t: np.ndarray, x: np.ndarray):
    dt = float(np.median(np.diff(t)))
    v = np.gradient(x, dt)
    moving = np.abs(v) > 0.05
    moving = ndimage.binary_closing(moving, structure=np.ones(11, dtype=bool))
    moving = ndimage.binary_opening(moving, structure=np.ones(5, dtype=bool))
    edges = np.diff(np.r_[False, moving, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1)
    rows = []
    for a, b in zip(starts, stops):
        if b - a < 5:
            continue
        dx = x[min(b, len(x) - 1)] - x[a]
        rows.append((t[a], t[min(b, len(t) - 1)], x[a],
                     x[min(b, len(x) - 1)], dx, dx / (t[b - 1] - t[a])))
    return v, rows


def window_stats(t: np.ndarray, x: np.ndarray, a: float, b: float):
    use = (t >= a) & (t <= b)
    xx = x[use]
    tt = t[use]
    slope = np.polyfit(tt - tt.mean(), xx, 1)[0] if len(xx) > 1 else np.nan
    return {
        "median": float(np.median(xx)),
        "pp_um": float(np.ptp(xx) * 1000),
        "slope_um_s": float(slope * 1000),
    }


def infer_profiles(t: np.ndarray, x: np.ndarray, raw_step_nm: np.ndarray, rows):
    """Infer profile stops from the known four-profile sequence.

    A terminal log was not supplied. The 16 full 59.35 mm transfers after the
    switch baseline are unambiguous in the IFM trace. Every group of four is:
    return 1, outward 2, return 2, then closure + outward 1 of the next profile.
    """
    full = [r for r in rows if r[0] > 650 and abs(r[4]) > 55]
    if len(full) != 16:
        raise RuntimeError(f"expected 16 diagnostic transfers, found {len(full)}")
    result = []
    for p, name in enumerate(PROFILE_NAMES):
        k = 4 * p
        # Stable gaps around the three pure moves. Trim 0.2 s from each edge.
        if p == 0:
            high1_window = (full[k][0] - 1.4, full[k][0] - 0.2)
        else:
            high1_window = (full[k - 1][1] + 0.2, full[k][0] - 0.2)
        low1_window = (full[k][1] + 0.2, full[k + 1][0] - 0.2)
        high2_window = (full[k + 1][1] + 0.2, full[k + 2][0] - 0.2)
        low2_window = (full[k + 2][1] + 0.2, full[k + 3][0] - 0.2)
        stops = [window_stats(t, x, *w) for w in
                 (high1_window, low1_window, high2_window, low2_window)]
        spans = [stops[0]["median"] - stops[1]["median"],
                 stops[2]["median"] - stops[3]["median"]]
        move = full[k]  # first return is never merged with a switch closure
        use = (t[:-1] >= move[0]) & (t[:-1] <= move[1])
        step = np.abs(raw_step_nm[use])
        result.append({
            "name": name,
            "stops": stops,
            "spans": spans,
            "return_s": move[1] - move[0],
            "return_speed": abs(move[4]) / (move[1] - move[0]),
            "clip_pct": 100.0 * float(np.mean(step == 2500)),
            "visible_holds_s": [high1_window[1] - high1_window[0],
                                low1_window[1] - low1_window[0],
                                high2_window[1] - high2_window[0],
                                low2_window[1] - low2_window[0]],
            "start_s": high1_window[0],
            "end_s": full[k + 3][0],
        })
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    scope = read_scope(args.source)
    x = scope.displacement_mm
    t = scope.time_s
    block = max(1, round(100.0 / scope.sample_ms))  # 100 ms
    tb = block_mean(t, block)
    xb = block_mean(x, block)
    vb, rows = motion_segments(tb, xb)
    raw_step_nm = np.diff(scope.counter_unwrapped)
    profiles = infer_profiles(t, x, raw_step_nm, rows)
    home_i, _ = find_peaks(-xb, prominence=0.8,
                           distance=max(1, round(8.0 / np.median(np.diff(tb)))))
    home_i = home_i[(tb[home_i] > 350) & (xb[home_i] < -60.5)]

    print(f"samples={len(t)} duration_s={t[-1]:.3f} span_mm={np.ptp(x):.6f} "
          f"rollovers={scope.rollover_count}")
    print("start_s,end_s,duration_s,x0_mm,x1_mm,dx_mm,mean_v_mm_s")
    for row in rows:
        print(f"{row[0]:.2f},{row[1]:.2f},{row[1]-row[0]:.2f},"
              f"{row[2]:.5f},{row[3]:.5f},{row[4]:.5f},{row[5]:.5f}")
    print("\nprofile,cycle,span_mm,error_um,high_pp_um,low_pp_um")
    for p in profiles:
        for cycle in range(2):
            hi = p["stops"][2 * cycle]
            lo = p["stops"][2 * cycle + 1]
            span = p["spans"][cycle]
            print(f"{p['name']},{cycle + 1},{span:.6f},"
                  f"{(span - EXPECTED_SPAN_MM) * 1000:.2f},"
                  f"{hi['pp_um']:.2f},{lo['pp_um']:.2f}")
        print(f"profile_summary,{p['name']},return_s={p['return_s']:.2f},"
              f"ifm_speed_mm_s={p['return_speed']:.4f},clip_pct={p['clip_pct']:.1f},"
              f"visible_holds_s={','.join(f'{v:.1f}' for v in p['visible_holds_s'])}")
    print("\nhome_switch_like_minima_s,position_mm")
    for i in home_i:
        print(f"{tb[i]:.2f},{xb[i]:.6f}")
    if len(home_i):
        print(f"home_minima_range_um={np.ptp(xb[home_i]) * 1000:.2f}")

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(15, 8), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(tb / 60, xb, lw=0.9, color="#1f5a94")
    ax.set_ylabel("IFM displacement from record start (mm)")
    ax.grid(True, alpha=0.25)
    bx.plot(tb / 60, vb, lw=0.65, color="#b34c27")
    bx.axhline(0, color="0.3", lw=0.6)
    bx.set_ylabel("100 ms velocity\n(mm/s)")
    bx.set_xlabel("time from recording start (min)")
    bx.grid(True, alpha=0.25)
    fig.suptitle("Open-loop diagnostic: complete raw interferometer trajectory")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(f"wrote={args.output}")

    # Result-focused view: diagnostic profiles only, plus the numerical verdict.
    result_output = args.output.with_name("iso_open_loop_diagnostic_profiles.png")
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_axes([0.07, 0.49, 0.90, 0.34])
    bx = fig.add_axes([0.07, 0.30, 0.90, 0.15], sharex=ax)
    tx = fig.add_axes([0.07, 0.055, 0.90, 0.17])
    colours = ("#a33f30", "#d08927", "#3b78a8", "#33805b")
    for p, label, colour in zip(profiles, PROFILE_LABELS, colours):
        ax.axvspan(p["start_s"] / 60, p["end_s"] / 60, color=colour,
                   alpha=0.08, lw=0)
        ax.text((p["start_s"] + p["end_s"]) / 120, 1.015, label,
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=9, color=colour)
    show = (tb >= profiles[0]["start_s"] - 4) & (tb <= profiles[-1]["end_s"] + 5)
    ax.plot(tb[show] / 60, xb[show], color="#194f80", lw=1.0)
    ax.set_ylabel("IFM displacement\nfrom record start (mm)")
    ax.grid(True, alpha=0.22)
    bx.plot(tb[show] / 60, vb[show], color="#b34c27", lw=0.7)
    bx.axhline(2.5, color="#8f1d14", lw=0.9, ls=(0, (5, 3)))
    bx.axhline(-2.5, color="#8f1d14", lw=0.9, ls=(0, (5, 3)))
    bx.text(0.995, 0.92, "recorded derivative ceiling: ±2.500 mm/s",
            transform=bx.transAxes, ha="right", va="top", fontsize=9,
            color="#8f1d14")
    bx.set_ylabel("100 ms velocity\n(mm/s)")
    bx.set_xlabel("time from recording start (min)")
    bx.grid(True, alpha=0.22)
    tx.axis("off")
    table_rows = []
    for p, label in zip(profiles, PROFILE_LABELS):
        span = np.mean(p["spans"])
        repeat = abs(np.diff(p["spans"])[0]) * 1000
        hold = np.mean(p["visible_holds_s"])
        table_rows.append([label, f"{span:.6f}",
                           f"{(span - EXPECTED_SPAN_MM) * 1000:+.1f}",
                           f"{repeat:.2f}", f"{p['return_speed']:.3f}",
                           f"{hold:.1f}"])
    table = tx.table(cellText=table_rows,
                     colLabels=["profile", "mean span (mm)", "error (µm)",
                                "cycle difference (µm)", "recorded ramp (mm/s)",
                                "visible stable hold (s)"],
                     cellLoc="center", colLoc="center", loc="center",
                     colWidths=[0.28, 0.13, 0.11, 0.16, 0.15, 0.17])
    table.auto_set_font_size(False)
    table.set_fontsize(8.7)
    table.scale(1, 1.45)
    fig.suptitle("Open-loop trajectory diagnostic — lower speed resolves the apparent shortfall",
                 x=0.07, y=0.965, ha="left", fontsize=16, fontweight="bold")
    fig.text(0.07, 0.905,
             "All four profiles reach the same endpoints. At both 5 mm/s profiles the recorded "
             "counter remains pinned near 2.5 mm/s through most of the 10 s hold.\n"
             "At 2 and 1 mm/s the full plateau is visible.", fontsize=10.5,
             linespacing=1.35)
    fig.text(0.07, 0.018,
             "Profile boundaries inferred from the known sequence because no matching terminal "
             "log was supplied. Endpoint displacement is physical IFM evidence; requested/emitted "
             "STEP counts cannot be audited from this CSV alone.", fontsize=8.5, color="0.35")
    fig.savefig(result_output, dpi=180)
    print(f"wrote={result_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
