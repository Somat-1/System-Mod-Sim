"""Render the v5 microstepping-only IDS capture: commanded vs realised step
size at MRES 1, 4, 16, 32.

Source: `../../MREStesting.csv`, an EL5101 capture at 1 kHz of
`esp32_v5_microstepping_only` (SpreadCycle, MicroPlyer disabled, 0.11 ohm
sense resistor, ~336 mA RMS). Each MRES block is: a small CONFIG marker, a
slightly larger OSCILLATION marker (amplitude grows by 4 full steps every
marker across the whole campaign), then 15 cycles of one-microstep
forward/return, 1 s dwell each side.

Anchor detection deliberately does NOT use a single global threshold (an
earlier version anchored on `median(position) - 50 um`, which broke for the
first block: the two markers there do not return the stage all the way to 0
-- there is a large, real residual offset -- so the global median is itself
dragged well negative by time spent off-baseline, and the "below threshold"
region spanning both markers plus the residual crept past the threshold at
the wrong point). Instead:

  1. Find large (>30 um over a 60 ms window) transitions -- these are the
     marker moves, unambiguous at any MRES since marker amplitudes are always
     >=120 um. Group them by a >15 s gap into one group per MRES block; the
     block's oscillation follows its last large transition.
  2. Within a short window after that, look for small transitions sized to
     that block's own commanded step (a threshold scaled to the commanded
     step, not a fixed value, since the step shrinks from 10 um at MRES 1 to
     0.3 um at MRES 32). The first run of several such transitions spaced at
     the nominal ~1.0 s half-dwell IS the oscillation; its first edge is the
     true osc0. A block or two show one extra "creep" transition right after
     the marker, before the regular pattern starts -- this discards it rather
     than anchoring on it.
  3. MRES 32 has no such run: the commanded 0.3 um step produced no
     detectable motion anywhere in its ~30 s block (see
     `oscillation_by_mres.png` -- it is visibly flat). When no regular run is
     found, the anchor falls back to last-transition + a fixed settle; it
     does not much matter where within a flat block it lands.

The first oscillation cycle of every block is excluded from the statistics:
it sits on the settling transient from the marker's move. It is still drawn,
shaded.
"""
from pathlib import Path

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from ids_common import (ANALYSIS_DIR, FULL_STEP_UM, C_MEASURED, C_COMMANDED,
                        C_INK, C_MUTED, apply_style, load_ids,
                        plot_full_capture, decimate_minmax)

OUT_DIR = ANALYSIS_DIR / "microstepping_only"
CSV = Path(__file__).resolve().parents[2] / "MREStesting.csv"

MRES_LADDER = [1, 4, 16, 32]
CYCLES = 15
HALF_DWELL_S = 1.0
NOMINAL_CYCLE_S = 2.0 * HALF_DWELL_S
SKIP_CYCLES = 1  # First cycle sits on the marker's settling transient.

BLOCK_GAP_S = 15.0     # Separates one MRES block's markers from the next's.
LARGE_JUMP_UM = 30.0   # Always well below a marker (>=120 um), above noise.
EDGE_WINDOW_S = 0.060  # ~60 ms: shorter than any dwell, longer than a pulse.
REGULAR_LO_S, REGULAR_HI_S = 0.85, 1.25  # Tolerance band around HALF_DWELL_S.
FALLBACK_SETTLE_S = 0.6


def windowed_edges(t, pos, threshold_um, start_t, end_t, window_s=EDGE_WINDOW_S):
    """Indices where |pos[i+w] - pos[i]| exceeds threshold_um, collapsed to
    one index per contiguous crossing (>=150 ms apart)."""
    dt = t[1] - t[0]
    w = max(1, int(round(window_s / dt)))
    lo = np.searchsorted(t, start_t)
    hi = min(np.searchsorted(t, end_t), pos.size - w)
    if hi <= lo:
        return np.array([], dtype=int)
    d = pos[lo + w:hi + w] - pos[lo:hi]
    idx = np.flatnonzero(np.abs(d) > threshold_um) + lo + w
    out = []
    for e in idx:
        if not out or t[e] - t[out[-1]] > 0.15:
            out.append(e)
    return np.array(out, dtype=int)


def group_by_gap(t, idx, gap_s):
    if idx.size == 0:
        return []
    groups = [[idx[0]]]
    for e in idx[1:]:
        if t[e] - t[groups[-1][-1]] > gap_s:
            groups.append([e])
        else:
            groups[-1].append(e)
    return groups


def find_osc0(t, pos, last_marker_edge, mres):
    """First edge of the longest run of consecutive small transitions spaced
    within [REGULAR_LO_S, REGULAR_HI_S] -- i.e. the real oscillation cadence,
    as opposed to marker settling. Falls back to a fixed offset if the
    commanded step never produced a detectable run (MRES 32)."""
    step_um = FULL_STEP_UM / mres
    small_edges = windowed_edges(t, pos, step_um * 0.3,
                                 t[last_marker_edge] + 0.05,
                                 t[last_marker_edge] + 5.0)
    if small_edges.size >= 2:
        gaps = np.diff(t[small_edges])
        good = (gaps >= REGULAR_LO_S) & (gaps <= REGULAR_HI_S)
        best_start, best_len, cur_start, cur_len = None, 0, None, 0
        for i, g in enumerate(good):
            if g:
                if cur_len == 0:
                    cur_start = i
                cur_len += 1
                if cur_len > best_len:
                    best_len, best_start = cur_len, cur_start
            else:
                cur_len = 0
        if best_start is not None:
            return t[small_edges[best_start]], True
    return t[last_marker_edge] + FALLBACK_SETTLE_S, False


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t, pos, meta = load_ids(CSV)
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {meta['fs_hz']:.0f} Hz "
          f"= {meta['duration_s']:.2f} s")

    large = windowed_edges(t, pos, LARGE_JUMP_UM, t[0], t[-1])
    groups = group_by_gap(t, large, BLOCK_GAP_S)
    anchors = [g[-1] for g in groups]
    print(f"marker-block anchors ({len(anchors)}): {np.round(t[anchors], 3)}")
    if len(anchors) != len(MRES_LADDER):
        print(f"WARNING: found {len(anchors)} block anchors, expected "
              f"{len(MRES_LADDER)}")

    osc0_list, resolved_run = [], []
    for a, mres in zip(anchors, MRES_LADDER):
        osc0, found_run = find_osc0(t, pos, a, mres)
        osc0_list.append(osc0)
        resolved_run.append(found_run)
        print(f"  MRES {mres:2d}: last marker edge {t[a]:8.3f} s -> osc0 "
              f"{osc0:8.3f} s ({'regular run found' if found_run else 'FALLBACK: no regular run (flat block)'})")

    plot_full_capture(
        t, pos, meta, OUT_DIR / "overview_full_capture.png",
        "MREStesting — complete EL5101 capture (v5 microstepping-only)",
        subtitle=(f"{meta['n']:,} samples at {meta['fs_hz']:.0f} Hz  ·  "
                  f"{meta['duration_s']:.1f} s  ·  1 nm/count  ·  "
                  "MRES 1, 4, 16, 32, SpreadCycle, MicroPlyer off, 0.11 ohm/336 mA"),
        marks=anchors)
    print("wrote overview_full_capture.png")

    cycle_period = NOMINAL_CYCLE_S  # Per-cycle timing has no UART calls; see header.

    # Per-MRES small multiples, normalised to fraction of the commanded step
    # so all four panels share one comparable scale regardless of MRES.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.6), sharex=False, sharey=True)
    summary = []
    w = int(round(1.0 / meta["dt_s"]))
    stds = np.array([pos[i:i + w].std() for i in range(0, pos.size - w, w)])
    noise_um = float(stds.min())

    for ax, osc0, mres in zip(axes.ravel(), osc0_list, MRES_LADDER):
        step_um = FULL_STEP_UM / mres
        base_m = (t >= osc0 - 0.40) & (t <= osc0 - 0.05)
        baseline = float(np.mean(pos[base_m])) if base_m.any() else 0.0

        osc_end = osc0 + CYCLES * cycle_period
        lo, hi = osc0 - 0.5, osc_end - 0.10
        sel = (t >= lo) & (t <= hi)
        y_norm = (pos[sel] - baseline) / step_um
        ax.plot(t[sel] - osc0, y_norm, color=C_MEASURED, lw=0.8, label="Measured")
        if SKIP_CYCLES:
            ax.axvspan(-0.5, SKIP_CYCLES * cycle_period, color=C_MUTED,
                       alpha=0.13, lw=0, zorder=0, label="excluded (settling)")
        ax.axhline(1.0, color=C_COMMANDED, lw=1.2, ls="--", zorder=2,
                   label="Commanded (1.0)")
        ax.axhline(0.0, color=C_COMMANDED, lw=1.2, ls="--", zorder=2)

        ups, up_t, up_span, dns, dn_t, dn_span = [], [], [], [], [], []
        for k in range(SKIP_CYCLES, CYCLES):
            for phase, acc, acc_t, acc_span in (
                    (0.0, ups, up_t, up_span), (0.5, dns, dn_t, dn_span)):
                a = osc0 + (k + phase) * cycle_period + 0.20 * cycle_period
                b = osc0 + (k + phase) * cycle_period + 0.45 * cycle_period
                m = (t >= a) & (t <= b)
                if m.any():
                    acc.append(float(np.mean(pos[m] - baseline)))
                    acc_t.append(0.5 * (a + b))
                    acc_span.append((a - osc0, b - osc0))

        # Average-value overlay: a short horizontal mark at the mean level
        # actually used for the statistics, spanning the window it was
        # averaged over -- makes the plateau-mean method visible, not just
        # asserted in the title.
        for span, val in zip(up_span, ups):
            ax.hlines(val / step_um, span[0], span[1], color=C_INK, lw=2.4,
                      zorder=6, alpha=0.85)
        for span, val in zip(dn_span, dns):
            ax.hlines(val / step_um, span[0], span[1], color=C_INK, lw=2.4,
                      zorder=6, alpha=0.85)

        if len(dns) >= 2:
            fit = np.polyfit(dn_t, dns, 1)
            steps = np.asarray(ups) - np.polyval(fit, up_t)
        else:
            steps = np.asarray(ups) - (np.mean(dns) if dns else 0.0)
        measured_step = float(np.mean(steps))
        step_sd = float(np.std(steps, ddof=1)) if steps.size > 1 else 0.0
        step_se = step_sd / np.sqrt(steps.size) if steps.size else np.inf
        dof = max(1, steps.size - 1)
        t_stat = measured_step / step_se if step_se > 0 else 0.0
        t_crit = float(scipy.stats.t.ppf(0.995, dof))
        resolved = abs(t_stat) > t_crit
        summary.append((mres, step_um, measured_step, step_sd, step_se,
                        t_stat, dof, resolved))
        ax.set_title(f"MRES {mres}  ·  commanded {step_um:.4g} µm  ·  "
                     f"measured {measured_step:.3g} µm ({measured_step/step_um:.2f}x)",
                     loc="left")
        if not resolved:
            ax.text(0.985, 0.06, f"not resolved  (t={t_stat:.1f}, crit {t_crit:.1f})",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8.5, color=C_MUTED)
        ax.set_xlabel("Time from oscillation start [s]")
        ax.set_ylabel("Position ÷ commanded step  [dimensionless]")
        ax.set_xlim(-0.5, osc_end - osc0 - 0.10)
        ax.set_ylim(-1.6, 1.6)
        if mres == MRES_LADDER[0]:
            ax.legend(loc="upper right", fontsize=7.5)
    fig.suptitle("MREStesting — one-microstep oscillation, measured vs commanded"
                 "\nnormalised to fraction of the commanded step (shared 1.0 = target scale); "
                 "black marks are the plateau means actually used for the statistics",
                 x=0.005, ha="left", fontsize=11.5, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT_DIR / "oscillation_by_mres.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote oscillation_by_mres.png")

    # Commanded vs realised step size across the ladder.
    mres_v = np.array([s[0] for s in summary], dtype=float)
    cmd_v = np.array([s[1] for s in summary])
    meas_v = np.array([s[2] for s in summary])
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
        t_crit = float(scipy.stats.t.ppf(0.995, dof))
        print(f"{m:5d}  {c:12.4f}  {v:12.4f}  {sd:7.4f}  {se:7.4f}  {v/c:6.3f} "
              f"{ts:6.1f} {t_crit:6.2f}  {'yes' if res else 'NO'}")


if __name__ == "__main__":
    main()
