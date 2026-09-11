"""Render the TMCStepping IDS capture: the MRES oscillation sweep.

Source: `data/hardware_runs/TMCStepping.csv`, an EL5101 capture at 1 kHz of the
`tmp_v4_oscillation_mres_sweep` firmware (SpreadCycle, MicroPlyer disabled,
MRES 1-2-4-8-16-32, 20 s per size, 2 s dwell at each end of the cycle).

The commanded sequence is known exactly from the firmware, so it is overlaid on
the measurement rather than merely described. Both traces are positions in µm,
so they share one axis - never a second y-scale.

Timing alignment: each block is anchored on its own measured marker RETURN
edge — the sharp crossing back up from -200 um — plus the firmware's 0.5 s
settle. Anchoring on the marker *onset* instead means adding a nominal 1.7667 s
spanning two rate-dependent moves and two dwells, which measured a consistent
-67 ms too long across all six blocks and biased every downstream window.

Within a block the oscillation cycle period is taken from the measured block
spacing rather than the nominal 4.000 s, because the firmware spends ~0.3 s per
block on UART register reads and logging that the nominal schedule does not
model.

The first cycle of each block is excluded from the statistics: it sits on the
settling transient from the marker's 200 um move, which at fine MRES is larger
than the commanded step. It is still drawn, shaded, so the exclusion is visible.
"""
import numpy as np
import matplotlib.pyplot as plt

from ids_common import (ANALYSIS_DIR, RUNS_DIR, FULL_STEP_UM, C_MEASURED,
                        C_COMMANDED, C_INK, C_MUTED, apply_style, load_ids,
                        find_marker_onsets, plot_full_capture, decimate_minmax)

OUT_DIR = ANALYSIS_DIR / "tmc_stepping"
CSV = RUNS_DIR / "TMCStepping.csv"

# Firmware constants (tmp_v4_oscillation_mres_sweep.ino).
MRES_LADDER = [1, 2, 4, 8, 16, 32]
MARKER_FULL_STEPS = 20
MARKER_RATE_FSPS = 150.0
MARKER_REVERSE_DWELL_S = 1.0
MARKER_SETTLE_S = 0.5
CYCLES = 5
HALF_DWELL_S = 2.0
MARKER_MOVE_S = MARKER_FULL_STEPS / MARKER_RATE_FSPS          # 0.1333 s
MARKER_TOTAL_S = 2 * MARKER_MOVE_S + MARKER_REVERSE_DWELL_S + MARKER_SETTLE_S

# The first cycle sits on the stage's settling transient from the marker's
# 200 um move, which at fine MRES is larger than the commanded step itself.
# Exclude it from the statistics; it is still drawn, shaded, so the
# exclusion is visible rather than silent.
SKIP_CYCLES = 1

# Two-sided 99% critical values of Student's t, indexed by degrees of
# freedom. Dropping the first cycle leaves n=4, so df=3 and the criterion is
# 5.841 - a fixed '3 standard errors' rule is far too loose at this sample
# size and manufactures detections out of scatter.
T_CRIT_99 = {1: 63.657, 2: 9.925, 3: 5.841, 4: 4.604, 5: 4.032,
             6: 3.707, 7: 3.499, 8: 3.355, 9: 3.250, 10: 3.169}


def marker_return_index(pos, onset_idx, level):
    """First sample after the marker onset where the axis has come back up.

    Anchoring the oscillation on the marker *onset* means adding a nominal
    1.767 s that contains two rate-dependent moves and two dwells, and every
    downstream window inherits that error. The return edge is a sharp 200 um
    crossing only MARKER_SETTLE_S away from the oscillation, so it is a much
    tighter anchor."""
    j = int(onset_idx)
    n = pos.size
    while j < n and pos[j] >= level:      # still descending into the marker
        j += 1
    while j < n and pos[j] < level:       # down at the marker floor
        j += 1
    return min(j, n - 1)


def commanded_block(t0, cycle_period_s, mres):
    """Commanded position (µm) for one MRES block, as (times, positions)."""
    step_um = FULL_STEP_UM / mres
    m = MARKER_FULL_STEPS * FULL_STEP_UM
    ts = [t0, t0 + MARKER_MOVE_S,
          t0 + MARKER_MOVE_S + MARKER_REVERSE_DWELL_S,
          t0 + 2 * MARKER_MOVE_S + MARKER_REVERSE_DWELL_S,
          t0 + MARKER_TOTAL_S]
    ys = [0.0, -m, -m, 0.0, 0.0]
    osc0 = t0 + MARKER_TOTAL_S
    for k in range(CYCLES):
        up = osc0 + k * cycle_period_s
        dn = up + cycle_period_s / 2.0
        ts += [up, up, dn, dn]
        ys += [0.0, step_um, step_um, 0.0]
    end = osc0 + CYCLES * cycle_period_s
    ts += [end]
    ys += [0.0]
    return np.asarray(ts), np.asarray(ys)


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t, pos, meta = load_ids(CSV)
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz "
          f"= {meta['duration_s']:.2f} s")

    onsets = find_marker_onsets(t, pos)
    print(f"marker onsets: {np.round(t[onsets], 3)}")
    if len(onsets) != len(MRES_LADDER):
        print(f"WARNING: found {len(onsets)} markers, expected {len(MRES_LADDER)}")

    # 1. The whole raw capture, always.
    plot_full_capture(
        t, pos, meta, OUT_DIR / "overview_full_capture.png",
        "TMCStepping — complete EL5101 capture",
        subtitle=(f"{meta['n']:,} samples at {meta['fs_hz']:.0f} Hz  ·  "
                  f"{meta['duration_s']:.1f} s  ·  1 nm/count  ·  "
                  "MRES sweep 1-2-4-8-16-32, SpreadCycle, MicroPlyer off"),
        marks=onsets)
    print("wrote overview_full_capture.png")

    # Per-block anchors, computed once and shared by every plot below so
    # the overlay and the small multiples cannot disagree about where a
    # block starts.
    level = float(np.median(pos)) - 50.0
    returns = [marker_return_index(pos, o, level) for o in onsets]
    osc0_list = [t[r] + MARKER_SETTLE_S for r in returns]
    print(f"\nalignment: anchoring on the marker return edge "
          f"(was onset + {MARKER_TOTAL_S:.4f} s nominal)")
    for o, r, a, mres in zip(onsets, returns, osc0_list, MRES_LADDER):
        print(f"  MRES {mres:2d}: return edge {t[r]:8.3f} s -> osc0 {a:8.3f} s "
              f"(onset-anchored {t[o] + MARKER_TOTAL_S:8.3f} s, "
              f"shift {1000*(a - t[o] - MARKER_TOTAL_S):+7.1f} ms)")

    # Cycle period measured from block spacing (see module docstring).
    if len(onsets) >= 2:
        block_period = float(np.median(np.diff(t[onsets])))
    else:
        block_period = MARKER_TOTAL_S + CYCLES * 2 * HALF_DWELL_S
    cycle_period = (block_period - MARKER_TOTAL_S) / CYCLES
    print(f"block period {block_period:.3f} s -> cycle period "
          f"{cycle_period:.3f} s (nominal {2*HALF_DWELL_S:.3f} s)")

    # 2. Full sweep with the commanded sequence overlaid.
    fig, ax = plt.subplots(figsize=(13, 4.6))
    td, yd = decimate_minmax(t, pos, target=6000)
    ax.plot(td, yd, color=C_MEASURED, lw=0.8, label="Measured (EL5101)", zorder=3)
    for i, (a, mres) in enumerate(zip(osc0_list, MRES_LADDER)):
        ct, cy = commanded_block(a - MARKER_TOTAL_S, cycle_period, mres)
        ax.plot(ct, cy, color=C_COMMANDED, lw=1.3, ls="--", zorder=4,
                label="Commanded" if i == 0 else None)
        ax.text(a + block_period / 2 - MARKER_TOTAL_S, 40, f"MRES {mres}",
                ha="center", fontsize=8.5, color=C_INK, zorder=5)
    ax.set_xlabel("Time from start of capture [s]")
    ax.set_ylabel("Position [µm]")
    ax.set_xlim(t[onsets[0]] - 4, t[onsets[-1]] + block_period + 4)
    ax.set_title("TMCStepping — measured vs commanded sequence",
                 loc="left", pad=24)
    ax.text(0.0, 1.012, "Markers are 20 full steps (200 µm); oscillations are "
                       "one microstep, 2 s at each end",
            transform=ax.transAxes, fontsize=8.5, color=C_MUTED,
            va="bottom")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "commanded_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote commanded_overlay.png")

    # 3. Per-MRES small multiples of the oscillation plateaus.
    #
    # Each block is re-zeroed on its own pre-oscillation baseline. The axis
    # does not return exactly to the same absolute position between blocks, so
    # a global zero would push the measured trace off-scale at fine MRES and
    # invite a false reading. The oscillation is a relative measurement.
    fig, axes = plt.subplots(3, 2, figsize=(13, 9.5), sharex=False)
    summary = []
    # Encoder noise floor: the quietest 1 s window anywhere in the capture.
    # Taking a window "just before the first marker" is not safe - that region
    # is not necessarily quiet, and using it inflated the floor 600-fold.
    w = int(round(1.0 / meta["dt_s"]))
    stds = np.array([pos[i:i + w].std() for i in range(0, pos.size - w, w)])
    noise_um = float(stds.min())
    level = float(np.median(pos)) - 50.0
    print(f"\nalignment: anchoring on the marker return edge "
          f"(was onset + {MARKER_TOTAL_S:.4f} s nominal)")
    for ax, o, mres in zip(axes.ravel(), onsets, MRES_LADDER):
        step_um = FULL_STEP_UM / mres
        ret = marker_return_index(pos, o, level)
        osc0 = t[ret] + MARKER_SETTLE_S
        osc0_old = t[o] + MARKER_TOTAL_S
        print(f"  MRES {mres:2d}: return edge {t[ret]:8.3f} s -> osc0 {osc0:8.3f} s "
              f"(onset-anchored {osc0_old:8.3f} s, shift {1000*(osc0-osc0_old):+7.1f} ms)")
        base_m = (t >= osc0 - 0.40) & (t <= osc0 - 0.05)
        baseline = float(np.mean(pos[base_m])) if base_m.any() else 0.0

        osc_end = osc0 + CYCLES * cycle_period
        lo, hi = osc0 - 0.5, osc_end - 0.10
        sel = (t >= lo) & (t <= hi)
        y = pos[sel] - baseline
        ax.plot(t[sel] - osc0, y, color=C_MEASURED, lw=0.8, label="Measured")
        if SKIP_CYCLES:
            ax.axvspan(-0.5, SKIP_CYCLES * cycle_period, color=C_MUTED,
                       alpha=0.13, lw=0, zorder=0,
                       label="excluded (marker settling)")
        ct, cy = commanded_block(osc0 - MARKER_TOTAL_S, cycle_period, mres)
        cm = (ct >= lo) & (ct <= hi)
        ax.plot(ct[cm] - osc0, cy[cm], color=C_COMMANDED, lw=1.4, ls="--",
                label="Commanded")

        # Plateau means over the settled middle of each half-dwell.
        ups, up_t, dns, dn_t = [], [], [], []
        for k in range(SKIP_CYCLES, CYCLES):
            for phase, acc, acc_t in ((0.0, ups, up_t), (0.5, dns, dn_t)):
                a = osc0 + (k + phase) * cycle_period + 0.20 * cycle_period
                b = osc0 + (k + phase) * cycle_period + 0.45 * cycle_period
                m = (t >= a) & (t <= b)
                if m.any():
                    acc.append(float(np.mean(pos[m] - baseline)))
                    acc_t.append(0.5 * (a + b))

        # Detrend against the return plateaus: the stage creeps across a block,
        # and an un-detrended up-minus-down difference absorbs that creep.
        if len(dns) >= 2:
            fit = np.polyfit(dn_t, dns, 1)
            steps = np.asarray(ups) - np.polyval(fit, up_t)
        else:
            steps = np.asarray(ups) - (np.mean(dns) if dns else 0.0)
        measured_step = float(np.mean(steps))
        step_sd = float(np.std(steps, ddof=1)) if steps.size > 1 else 0.0
        step_se = step_sd / np.sqrt(steps.size) if steps.size else np.inf
        # Resolved when the mean step differs from zero at 99% on a
        # two-sided t-test over the cycles kept, not a fixed multiple of
        # the standard error.
        dof = max(1, steps.size - 1)
        t_stat = measured_step / step_se if step_se > 0 else 0.0
        t_crit = T_CRIT_99.get(dof, 3.0)
        resolved = abs(t_stat) > t_crit
        summary.append((mres, step_um, measured_step, step_sd, step_se,
                        t_stat, dof, resolved))
        ax.set_title(f"MRES {mres}  ·  commanded {step_um:.4g} µm  ·  "
                     f"measured {measured_step:.3g} µm", loc="left")
        if not resolved:
            ax.text(0.985, 0.06, f"not resolved  (t={t_stat:.1f}, crit {t_crit:.1f})",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8.5, color=C_MUTED)
        ax.set_xlabel("Time from oscillation start [s]")
        ax.set_ylabel("Position [µm]")
        ax.set_xlim(-0.5, osc_end - osc0 - 0.10)
        # Scale on the oscillation itself, excluding the marker ramps.
        core = y[(t[sel] >= osc0) & (t[sel] <= osc_end - 0.10)]
        span = max(np.ptp(core) if core.size else step_um,
                   step_um * 1.15, 8.0 * noise_um)
        mid = 0.5 * (min(core.min() if core.size else 0.0, 0.0)
                     + max(core.max() if core.size else step_um, step_um))
        ax.set_ylim(mid - 0.75 * span, mid + 0.75 * span)
        if mres == MRES_LADDER[0]:
            ax.legend(loc="upper right")
    fig.suptitle("TMCStepping — one-microstep oscillation, measured vs commanded"
                 "\neach block re-zeroed on its own pre-oscillation baseline",
                 x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_DIR / "oscillation_by_mres.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote oscillation_by_mres.png")

    # 4. Commanded vs realised step size across the ladder.
    mres_v = np.array([s[0] for s in summary], dtype=float)
    cmd_v = np.array([s[1] for s in summary])
    meas_v = np.array([s[2] for s in summary])
    sd_v = np.array([s[3] for s in summary])
    se_v = np.array([s[4] for s in summary])
    resolved_v = np.array([r[7] for r in summary])
    floor = 3.0 * noise_um

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    axL.plot(mres_v, cmd_v, color=C_COMMANDED, ls="--", marker="s", ms=7,
             lw=1.6, label="Commanded", zorder=3)
    resolved = resolved_v & (meas_v > 0)
    axL.errorbar(mres_v[resolved], meas_v[resolved], yerr=se_v[resolved],
                 color=C_MEASURED, marker="o", ms=8, lw=1.8, capsize=3,
                 label="Measured", zorder=4)
    # A log axis cannot show a null or negative result, and silently dropping
    # those points would overstate the measurement. Mark them at the noise
    # floor as upper bounds instead.
    if (~resolved).any():
        axL.scatter(mres_v[~resolved], np.full((~resolved).sum(), floor),
                    facecolors="white", edgecolors=C_MEASURED, s=70,
                    linewidths=1.6, zorder=5,
                    label="Not resolved (t-test, 99%)")
    axL.axhspan(1e-4, floor, color=C_MUTED, alpha=0.12, zorder=1)
    axL.text(mres_v[-1], floor * 0.75, "encoder noise floor", ha="right",
             va="top", fontsize=8, color=C_MUTED)
    axL.set_xscale("log", base=2)
    axL.set_yscale("log")
    axL.set_xticks(mres_v)
    axL.set_xticklabels([f"{int(m)}" for m in mres_v])
    axL.set_ylim(floor * 0.35, cmd_v.max() * 2.2)
    axL.set_xlabel("MRES (microsteps per full step)")
    axL.set_ylabel("Step size [µm]")
    axL.set_title("Commanded vs realised microstep size", loc="left")
    axL.legend(loc="lower left")

    frac = np.where(cmd_v > 0, meas_v / cmd_v, np.nan)
    axR.axhline(1.0, color=C_COMMANDED, ls="--", lw=1.4, label="Ideal (1.0)")
    axR.plot(mres_v, frac, color=C_MEASURED, marker="o", ms=8, lw=1.8,
             label="Realised fraction")
    for m, f in zip(mres_v, frac):
        axR.annotate(f"{f:.2f}", (m, f), textcoords="offset points",
                     xytext=(0, 10), ha="center", fontsize=8, color=C_INK)
    axR.set_xscale("log", base=2)
    axR.set_xticks(mres_v)
    axR.set_xticklabels([f"{int(m)}" for m in mres_v])
    axR.set_ylim(-0.2, 1.35)
    axR.set_xlabel("MRES (microsteps per full step)")
    axR.set_ylabel("Measured / commanded")
    axR.set_title("Fraction of the commanded step actually realised", loc="left")
    axR.legend(loc="lower left")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "step_size_vs_mres.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote step_size_vs_mres.png")

    print(f"\nencoder noise (1 sigma, quiet window) = {noise_um*1000:.1f} nm  "
          f"-> 3 sigma = {floor*1000:.2f} nm")
    print(f"cycles used: {CYCLES - SKIP_CYCLES} of {CYCLES} (first excluded: marker settling)")
    print("MRES  commanded[µm]  measured[µm]   sd[µm]   se[µm]  ratio      t   crit  resolved")
    for (m, c, v, sd, se, ts, dof, res) in summary:
        print(f"{m:5d}  {c:12.4f}  {v:12.4f}  {sd:7.4f}  {se:7.4f}  {v/c:6.3f} "
              f"{ts:6.1f} {T_CRIT_99.get(dof, 3.0):6.2f}  {'yes' if res else 'NO'}")


if __name__ == "__main__":
    main()
