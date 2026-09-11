"""Render the TMC2209ExpRun IDS capture: the full v4 MRES/trajectory campaign.

Source: `data/hardware_runs/TMC2209ExpRun.csv`, an EL5101 capture at 1 kHz
covering 17:20:06 - 18:03:00 on 2026-09-10, i.e. the complete 41 min campaign
(4 MRES values x [oscillation + 6 x 25 mm out-and-back trajectories]).

PROVENANCE WARNING. This run predates the SpreadCycle fix. The firmware on the
board at the time had MicroPlyer interpolation disabled but left the chopper at
its power-on default, so the campaign ran under **StealthChop**, not
SpreadCycle. See TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md.

No serial event log exists for this execution - it was started from the board's
own RUN command, not from a host capture - so blocks are segmented from the
measurement itself rather than read from timestamps. Every segment boundary
here is therefore inferred, and is labelled as such.
"""
import numpy as np
import matplotlib.pyplot as plt

from ids_common import (ANALYSIS_DIR, RUNS_DIR, C_MEASURED, C_COMMANDED,
                        C_INK, C_MUTED, apply_style, load_ids,
                        plot_full_capture, decimate_minmax)

OUT_DIR = ANALYSIS_DIR / "tmc2209_exp_run"
CSV = RUNS_DIR / "TMC2209ExpRun.csv"

MRES_VALUES = [1, 4, 16, 32]
RATE_NAMES = ["slow", "moderate", "fast"]
RATES_FSPS = [27.5, 70.0, 200.0]
TRAJECTORY_MM = 25.0
TRAJ_THRESHOLD_UM = 5000.0     # a 25 mm leg is 25 000 um; markers are <= 140 um
OSC_DURATION_S = 30.0


def find_trajectories(t, pos, dt_s):
    """Return (start, apex, end) index triples for each 25 mm out-and-back.

    The profile is triangular and the resting baseline drifts downward across
    the run, so the threshold is taken against a low percentile (the resting
    level), not the median, which on a triangle sits mid-ramp. The span is the
    above-threshold run itself: expanding each one down to the resting band
    makes neighbouring legs collide, because consecutive trajectories are
    separated only by a 1 s endpoint dwell and a small marker.
    """
    baseline = np.percentile(pos, 5.0)
    above = pos > (baseline + 0.32 * TRAJECTORY_MM * 1000.0)
    d = np.diff(above.astype(np.int8))
    starts = list(np.flatnonzero(d == 1) + 1)
    ends = list(np.flatnonzero(d == -1) + 1)
    if above[0]:
        starts.insert(0, 0)
    if above[-1]:
        ends.append(above.size - 1)

    spans = []
    for s0, e0 in zip(starts, ends):
        if (e0 - s0) * dt_s < 2.0:
            continue
        apex = s0 + int(np.argmax(pos[s0:e0 + 1]))
        spans.append((s0, apex, e0))
    return spans


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t, pos, meta = load_ids(CSV)
    pos_mm = pos / 1000.0
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz "
          f"= {meta['duration_s']:.1f} s ({meta['duration_s']/60:.2f} min)")

    # 1. The whole raw capture, always.
    plot_full_capture(
        t, pos, meta, OUT_DIR / "overview_full_capture.png",
        "TMC2209ExpRun — complete EL5101 capture (full v4 campaign)",
        subtitle=(f"{meta['n']:,} samples at {meta['fs_hz']:.0f} Hz  ·  "
                  f"{meta['duration_s']/60:.1f} min  ·  1 nm/count  ·  "
                  "MicroPlyer off, StealthChop (predates the SpreadCycle fix)"))
    print("wrote overview_full_capture.png")

    spans = find_trajectories(t, pos, meta["dt_s"])
    print(f"detected {len(spans)} trajectory excursions (expected 24)")

    # 2. Full capture in mm with the detected trajectories shaded.
    fig, ax = plt.subplots(figsize=(13.5, 4.8))
    td, yd = decimate_minmax(t, pos_mm, target=6000)
    ax.plot(td, yd, color=C_MEASURED, lw=0.7, label="Measured (EL5101)")
    for i, (s, _apex, e) in enumerate(spans):
        ax.axvspan(t[s], t[e], color=C_COMMANDED, alpha=0.13, lw=0,
                   label="Detected 25 mm trajectory" if i == 0 else None)
    ax.set_xlabel("Time from start of capture [s]")
    ax.set_ylabel("Position [mm]")
    ax.set_xlim(t[0], t[-1])
    ax.set_title("TMC2209ExpRun — campaign structure, 25 mm trajectories highlighted",
                 loc="left", pad=24)
    ax.text(0.0, 1.012, f"{len(spans)} excursions detected from the measurement; "
                       "no serial event log exists for this run",
            transform=ax.transAxes, fontsize=8.5, color=C_MUTED,
            va="bottom")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "trajectories_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote trajectories_overview.png")

    # 3. Per-trajectory montage: 24 out-and-back legs.
    if spans:
        n = len(spans)
        ncol = 6
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.3 * nrow),
                                 squeeze=False)
        for i, ax in enumerate(axes.ravel()):
            if i >= n:
                ax.axis("off")
                continue
            s, _apex, e = spans[i]
            pad = int(2.0 / meta["dt_s"])
            a, b = max(0, s - pad), min(pos.size - 1, e + pad)
            ax.plot(t[a:b] - t[s], pos_mm[a:b], color=C_MEASURED, lw=0.7)
            ax.axhline(0.0, color=C_MUTED, lw=0.6, ls=":")
            peak = pos_mm[s:e].max() - pos_mm[s:e].min()
            dur = (e - s) * meta["dt_s"]
            ax.set_title(f"#{i+1}  ·  {peak:.2f} mm  ·  {dur:.0f} s",
                         loc="left", fontsize=9)
            ax.tick_params(labelsize=7)
            if i % ncol == 0:
                ax.set_ylabel("Position [mm]", fontsize=8)
            if i >= n - ncol:
                ax.set_xlabel("Time in block [s]", fontsize=8)
        fig.suptitle("TMC2209ExpRun — every detected 25 mm out-and-back leg",
                     x=0.005, ha="left", fontsize=12, fontweight="semibold")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(OUT_DIR / "trajectory_montage.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print("wrote trajectory_montage.png")

    # 4. Oscillation blocks: the quiet 30 s stretches before each MRES group's
    #    first trajectory. With 6 trajectories per MRES value, group boundaries
    #    fall every 6 detected excursions.
    if len(spans) >= 6:
        group_firsts = [spans[i][0] for i in range(0, len(spans), 6)]
        fig, axes = plt.subplots(len(group_firsts), 1,
                                 figsize=(12, 2.6 * len(group_firsts)),
                                 squeeze=False)
        for ax, gf, mres in zip(axes.ravel(), group_firsts,
                                MRES_VALUES + [None] * 4):
            # Search back for the last long quiet stretch before this group.
            win = int(OSC_DURATION_S / meta["dt_s"])
            lo = max(0, gf - int(220.0 / meta["dt_s"]))
            seg = pos[lo:gf]
            if seg.size < win:
                ax.axis("off")
                continue
            # Oscillation = small-amplitude region: pick the window with the
            # smallest peak-to-peak that still contains motion.
            k = max(1, seg.size // 400)
            ptp = np.array([np.ptp(seg[i:i + win])
                            for i in range(0, seg.size - win, k)])
            best = int(np.argmin(np.where(ptp > 0.02, ptp, np.inf))) * k
            o0 = lo + best
            sl = slice(o0, o0 + win)
            y = pos[sl] - np.median(pos[sl])
            ax.plot(t[sl] - t[o0], y, color=C_MEASURED, lw=0.8)
            title = f"MRES {mres}" if mres else "block"
            ax.set_title(f"{title}  ·  inferred oscillation window  ·  "
                         f"peak-to-peak {np.ptp(y):.3f} µm", loc="left")
            ax.set_ylabel("Position [µm]")
            ax.set_xlabel("Time in window [s]")
        fig.suptitle("TMC2209ExpRun — one-microstep oscillation blocks "
                     "(windows inferred from the measurement)",
                     x=0.005, ha="left", fontsize=12, fontweight="semibold")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(OUT_DIR / "oscillation_blocks.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print("wrote oscillation_blocks.png")

    # Resting level between legs. Sampling a fixed time before the threshold
    # crossing lands mid-ramp and reports the leg's speed rather than the axis
    # position, so isolate genuine dwells: samples whose local slope is small.
    slope = np.abs(np.diff(pos, prepend=pos[0])) / meta["dt_s"]   # um/s
    still = slope < 20.0

    rests = []
    for i, (s0, _apex, _e0) in enumerate(spans):
        g0 = spans[i - 1][2] if i else 0
        seg = pos_mm[g0:s0]
        seg_still = still[g0:s0]
        rests.append(float(np.median(seg[seg_still]))
                     if seg_still.sum() >= 50 else float(np.median(seg)))
    rests = np.asarray(rests)
    apexes = np.array([pos_mm[a] for _s, a, _e in spans])
    travel = apexes - rests

    print("\n  #  start[s]   apex[s]    end[s]   dur[s]  apex[mm]  rest[mm]  travel[mm]")
    for i, (s0, apex, e0) in enumerate(spans):
        print(f"{i+1:3d} {t[s0]:9.1f} {t[apex]:9.1f} {t[e0]:9.1f} "
              f"{(e0-s0)*meta['dt_s']:7.1f} {apexes[i]:9.3f} {rests[i]:9.3f} "
              f"{travel[i]:10.3f}")

    print(f"\ntravel per leg: mean {travel.mean():.3f} mm, "
          f"sd {travel.std(ddof=1):.3f} mm  (commanded {TRAJECTORY_MM:.1f} mm)")
    print(f"resting level: {rests[0]:+.3f} -> {rests[-1]:+.3f} mm over "
          f"{len(spans)} legs = {(rests[-1]-rests[0])*1000:+.0f} um total, "
          f"{(rests[-1]-rests[0])*1000/len(spans):+.1f} um per leg")

    # 5. Travel and origin drift per leg - the headline result for this run.
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True)
    legs = np.arange(1, len(spans) + 1)
    a1.plot(legs, travel, color=C_MEASURED, marker="o", ms=6, lw=1.6,
            label="Measured travel")
    a1.axhline(TRAJECTORY_MM, color=C_COMMANDED, ls="--", lw=1.4,
               label=f"Commanded ({TRAJECTORY_MM:.0f} mm)")
    a1.set_ylabel("Travel per leg [mm]")
    a1.set_title("Travel per leg", loc="left")
    a1.legend(loc="lower left")

    a2.plot(legs, rests, color=C_MEASURED, marker="o", ms=6, lw=1.6,
            label="Measured resting level")
    a2.axhline(rests[0], color=C_COMMANDED, ls="--", lw=1.4,
               label="Commanded origin (constant)")
    a2.set_xlabel("Trajectory leg number (chronological)")
    a2.set_ylabel("Resting level [mm]")
    a2.set_title(f"Cumulative origin drift: {(rests[-1]-rests[0])*1000:+.0f} µm "
                 f"over the campaign", loc="left")
    a2.set_xticks(legs[::2])
    a2.legend(loc="lower left")
    for a in (a1, a2):
        for g in range(6, len(spans), 6):
            a.axvline(g + 0.5, color=C_MUTED, lw=0.8, ls=":")
    fig.suptitle("TMC2209ExpRun — per-leg travel and origin drift  "
                 "(dotted lines mark MRES group boundaries)",
                 x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT_DIR / "travel_and_drift.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote travel_and_drift.png")


if __name__ == "__main__":
    main()
