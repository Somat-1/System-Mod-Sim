"""Render the EXPrun IDS capture: the full v6 MRES/trajectory campaign.

Source: `../../EXPrun.csv`, an EL5101 capture at 1 kHz covering
13:48:58 - 14:10:40 on 2026-09-15, i.e. the complete ~21.7 min v6 campaign
(4 MRES values x [oscillation + 6 x 10 mm out-and-back trajectories]).

v6 is v4's campaign firmware (`esp32_v4_mres_trajectory_campaign`) forked with
two corrections: R_SENSE_OHM = 0.11 (the resistor actually fitted on this
board, not v4's 0.03, which under-delivered current ~3x) and a 10 mm leg
instead of v4's 25 mm (the stage's real travel is ~44 mm total; 25 mm one-way
does not fit). `configureDriver()` refuses to start the campaign unless
GCONF/CHOPCONF/DRV_STATUS all confirm SpreadCycle active and interpolation
off, so — unlike the pre-fix v4 TMC2209ExpRun capture — this run is not
suspected of having run under StealthChop. See `../../README.md`.

As with TMC2209ExpRun in v4, no serial event log accompanies this export, so
blocks are segmented from the measurement itself rather than read from
timestamps. Every segment boundary here is therefore inferred, and is
labelled as such.
"""
import numpy as np
import matplotlib.pyplot as plt

from ids_common import (ANALYSIS_DIR, V6_DIR, C_MEASURED, C_COMMANDED,
                        C_INK, C_MUTED, apply_style, load_ids,
                        plot_full_capture, decimate_minmax,
                        still_mask, ramp_extent)

OUT_DIR = ANALYSIS_DIR / "exp_run"
CSV = V6_DIR / "EXPrun.csv"

MRES_VALUES = [1, 4, 16, 32]
RATE_NAMES = ["slow", "moderate", "fast"]
RATES_FSPS = [27.5, 70.0, 200.0]
TRAJECTORY_MM = 10.0
OSC_DURATION_S = 30.0


def find_trajectories(t, pos, dt_s):
    """Return (start, apex, end) index triples for each 10 mm out-and-back.

    The profile is triangular and the resting baseline can drift across the
    run, so the threshold is taken against a low percentile (the resting
    level), not the median, which on a triangle sits mid-ramp. The span is the
    above-threshold run itself: expanding each one down to the resting band
    makes neighbouring legs collide, because consecutive trajectories are
    separated only by a short endpoint dwell and a small marker.
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


def leg_rest_levels(t, pos, spans, still, dt_s):
    """Resting position (um) immediately before and after each leg's ramp,
    each with the sample index it was measured at.

    Bounded by that leg's own `ramp_extent`, not the gap to the previous leg:
    the gap also contains that leg's marker and, at the start of each MRES
    group, the full 30 s oscillation block, so a median over the whole gap can
    land on an unrelated dwell instead of the genuine rest. See
    `plot_v6_rest_window_diagnostic.py`.

    Returns four arrays of length len(spans): idx_before, rest_before_um,
    idx_after, rest_after_um.
    """
    idx_before, rest_before, idx_after, rest_after = [], [], [], []
    for s0, _apex, e0 in spans:
        a, b = ramp_extent(still, s0, e0, dt_s, pad_s=1.0)
        bs, as_ = still[a:s0], still[e0:b]
        # Representative time for each rest point: the middle of the window
        # it was measured over. Only used to place the point on a time axis,
        # not part of the median itself.
        idx_before.append((a + s0) // 2)
        idx_after.append((e0 + b) // 2)
        rest_before.append(float(np.median(pos[a:s0][bs])) if bs.sum() > 20
                           else float(np.median(pos[a:s0])))
        rest_after.append(float(np.median(pos[e0:b][as_])) if as_.sum() > 20
                          else float(np.median(pos[e0:b])))
    return (np.asarray(idx_before), np.asarray(rest_before),
            np.asarray(idx_after), np.asarray(rest_after))


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
        "EXPrun — complete EL5101 capture (full v6 campaign)",
        subtitle=(f"{meta['n']:,} samples at {meta['fs_hz']:.0f} Hz  ·  "
                  f"{meta['duration_s']/60:.1f} min  ·  1 nm/count  ·  "
                  "MicroPlyer off, SpreadCycle, 0.11 ohm R_SENSE"))
    print("wrote overview_full_capture.png")

    spans = find_trajectories(t, pos, meta["dt_s"])
    print(f"detected {len(spans)} trajectory excursions (expected 24)")

    # 2. Full capture in mm with the detected trajectories shaded.
    fig, ax = plt.subplots(figsize=(13.5, 4.8))
    td, yd = decimate_minmax(t, pos_mm, target=6000)
    ax.plot(td, yd, color=C_MEASURED, lw=0.7, label="Measured (EL5101)")
    for i, (s, _apex, e) in enumerate(spans):
        ax.axvspan(t[s], t[e], color=C_COMMANDED, alpha=0.13, lw=0,
                   label="Detected 10 mm trajectory" if i == 0 else None)
    ax.set_xlabel("Time from start of capture [s]")
    ax.set_ylabel("Position [mm]")
    ax.set_xlim(t[0], t[-1])
    ax.set_title("EXPrun — campaign structure, 10 mm trajectories highlighted",
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

    # 3. Per-trajectory montage: 24 out-and-back legs, full extent.
    #
    # find_trajectories returns only the part of each ramp above the detection
    # threshold, so plotting that span alone cuts both ends off mid-ramp. Walk
    # out to where motion actually started and stopped, then show a second of
    # the dwell either side, and report the real travel (apex minus the resting
    # level) rather than the range inside the clipped window.
    still = still_mask(pos, meta["dt_s"])
    if spans:
        n = len(spans)
        ncol = 6
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.5 * nrow),
                                 squeeze=False)
        for i, ax in enumerate(axes.ravel()):
            if i >= n:
                ax.axis("off")
                continue
            s0, apex, e0 = spans[i]
            a, b = ramp_extent(still, s0, e0, meta["dt_s"], pad_s=1.0)
            seg_still = still[a:b]
            rest = (float(np.median(pos_mm[a:b][seg_still]))
                    if seg_still.sum() > 50 else float(np.min(pos_mm[a:b])))
            travel = pos_mm[apex] - rest

            ax.plot(t[a:b] - t[a], pos_mm[a:b], color=C_MEASURED, lw=0.8)
            ax.axhline(rest, color=C_MUTED, lw=0.7, ls=":")
            ax.axhline(rest + TRAJECTORY_MM, color=C_COMMANDED, lw=0.9,
                       ls="--", alpha=0.8)
            mode = "D" if (i % 6) < 3 else "I"
            rate = ["slow", "mod", "fast"][i % 3]
            ax.set_title(f"#{i+1}  MRES {MRES_VALUES[i//6]}  {mode}·{rate}"
                         f"\n{travel:.3f} mm  ·  {(b-a)*meta['dt_s']:.0f} s",
                         loc="left", fontsize=8.5)
            ax.tick_params(labelsize=7)
            ax.set_ylim(rest - 2.0, rest + TRAJECTORY_MM + 2.0)
            if i % ncol == 0:
                ax.set_ylabel("Position [mm]", fontsize=8)
            if i >= n - ncol:
                ax.set_xlabel("Time in block [s]", fontsize=8)
        fig.suptitle("EXPrun — every 10 mm out-and-back leg, full ramp "
                     "plus 1 s of dwell either side"
                     "\ndotted = resting level, dashed = commanded 10 mm; "
                     "D = DIRECT, I = INDIVIDUAL",
                     x=0.005, ha="left", fontsize=12, fontweight="semibold")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
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
        fig.suptitle("EXPrun — one-microstep oscillation blocks "
                     "(windows inferred from the measurement)",
                     x=0.005, ha="left", fontsize=12, fontweight="semibold")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(OUT_DIR / "oscillation_blocks.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print("wrote oscillation_blocks.png")

    # Resting level before and after each leg, from genuine dwells only
    # (still_mask). This is NOT the same as "the gap between leg i-1 and
    # leg i": that gap also contains the marker for leg i (a rapid move to
    # a large, leg-specific offset and back, held for its own ~1 s dwell),
    # and an earlier version of this script medianed the whole gap as one
    # window. Whichever cluster of still samples was larger -- the marker's
    # own dwell or the true pre-ramp rest -- won the median, and which one
    # that was varied leg to leg almost at random. That produced a spurious
    # ~0.07 mm sd in travel and made the MRES 1 group look uniquely
    # unrepeatable (rest swinging up to ~240 um leg-to-leg there, versus
    # ~2 um for MRES >= 4). It was a windowing bug, not a hardware effect:
    # see plot_v6_rest_window_diagnostic.py, which renders the raw signal in
    # both the old (whole-gap) and new (ramp-local) windows side by side for
    # the legs where the bug was most visible.
    #
    # The fix is to reuse each leg's own ramp_extent(a, b) -- already bounded
    # to "motion stopped for settle_s" on both sides of the ramp -- and take
    # the still-sample median separately from its leading edge (a to s0) and
    # trailing edge (e0 to b), rather than pooling a whole inter-leg gap that
    # may contain unrelated motion. Doing this for all 24 legs drops travel's
    # sd from 0.067 mm to 0.001 mm.
    _idx_b, rest_before_um, _idx_a, rest_after_um = leg_rest_levels(
        t, pos, spans, still, meta["dt_s"])
    rest_before, rest_after = rest_before_um / 1000.0, rest_after_um / 1000.0
    apexes = np.array([pos_mm[a] for _s, a, _e in spans])
    travel_out = apexes - rest_before      # outbound: rest -> apex
    travel_back = apexes - rest_after      # return: apex -> rest
    travel = 0.5 * (travel_out + travel_back)
    # Net drift contributed by this one leg (rest after minus rest before);
    # cumulative sum gives the same "origin drift" quantity the old plot
    # showed, now built leg-by-leg instead of read off noisy absolute levels.
    leg_drift = rest_after - rest_before
    cum_drift = np.cumsum(leg_drift)

    print("\n  #  start[s]   apex[s]    end[s]   dur[s]  apex[mm]  "
          "out[mm]  back[mm]  travel[mm]  leg_drift[um]")
    for i, (s0, apex, e0) in enumerate(spans):
        print(f"{i+1:3d} {t[s0]:9.1f} {t[apex]:9.1f} {t[e0]:9.1f} "
              f"{(e0-s0)*meta['dt_s']:7.1f} {apexes[i]:9.3f} "
              f"{travel_out[i]:8.4f} {travel_back[i]:9.4f} {travel[i]:11.4f} "
              f"{leg_drift[i]*1000:+9.1f}")

    print(f"\ntravel per leg: mean {travel.mean():.4f} mm, "
          f"sd {travel.std(ddof=1):.4f} mm  (commanded {TRAJECTORY_MM:.1f} mm)")
    print(f"cumulative origin drift: {cum_drift[-1]*1000:+.0f} um over "
          f"{len(spans)} legs = {cum_drift[-1]*1000/len(spans):+.2f} um/leg")

    # 5. Travel and origin drift per leg - the headline result for this run.
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True)
    legs = np.arange(1, len(spans) + 1)
    a1.plot(legs, travel_out, color=C_MEASURED, marker="o", ms=5, lw=1.4,
            alpha=0.55, label="Outbound (rest -> apex)")
    a1.plot(legs, travel_back, color=C_MEASURED, marker="s", ms=5, lw=1.4,
            alpha=0.55, ls="--", label="Return (apex -> rest)")
    a1.axhline(TRAJECTORY_MM, color=C_COMMANDED, ls="--", lw=1.4,
               label=f"Commanded ({TRAJECTORY_MM:.0f} mm)")
    a1.set_ylabel("Travel per leg [mm]")
    a1.set_title("Travel per leg, outbound vs return", loc="left")
    a1.legend(loc="lower left", fontsize=8)

    a2.plot(legs, cum_drift, color=C_MEASURED, marker="o", ms=6, lw=1.6,
            label="Cumulative origin drift")
    a2.axhline(0.0, color=C_COMMANDED, ls="--", lw=1.4,
               label="No drift (constant origin)")
    a2.set_xlabel("Trajectory leg number (chronological)")
    a2.set_ylabel("Cumulative drift [mm]")
    a2.set_title(f"Cumulative origin drift: {cum_drift[-1]*1000:+.0f} µm "
                 f"over the campaign", loc="left")
    a2.set_xticks(legs[::2])
    a2.legend(loc="lower left")
    for a in (a1, a2):
        for g in range(6, len(spans), 6):
            a.axvline(g + 0.5, color=C_MUTED, lw=0.8, ls=":")
    fig.suptitle("EXPrun — per-leg travel and origin drift  "
                 "(dotted lines mark MRES group boundaries; rest levels from "
                 "each leg's own ramp-local dwell, not the inter-leg gap)",
                 x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT_DIR / "travel_and_drift.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote travel_and_drift.png")


if __name__ == "__main__":
    main()
