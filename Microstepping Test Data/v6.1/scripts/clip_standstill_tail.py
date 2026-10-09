"""Clip the standstill tail out of v6's EXPrun capture, DC-removed, for
frequency analysis.

`../../v6/EXPrun.csv` runs 1302.3 s. The campaign's last commanded motion (leg
#24 of 24, MRES 32 I-fast) decelerates into its final rest at t=1148.67 s, and
the scope keeps recording for another 153.65 s with nothing commanded at all.
That trailing segment is what this script extracts.

Two things happen to it before export:

1. **A settling guard.** The detected boundary's first ~4 ms carry a single
   ~1.15 um transient (see `../../v6/Curve fitting/v2/noise_summary.json`,
   `large_outliers_gt_300nm`). Left in, that impulse smears broadband energy
   across the whole spectrum and swamps a 40 nm signal. `GUARD_S` seconds are
   dropped from the front; the render still shows the dropped span, so it is
   visible rather than quietly gone.
2. **DC removal.** The raw counter sits ~194 um below the capture's first
   sample, which is a parked position, not signal. The exported `pos_um_ac`
   column is the clip with its own mean subtracted, so bin 0 of a DFT does not
   dominate everything else.

Boundary detection reuses `ids_common.still_mask` -- the same windowed-slope
test the rest of v6 uses -- scanned backward from the end of the capture for
the earliest sample after which the axis is continuously still. Not a fixed
offset from the last leg; an actual detected boundary. It reproduces the
t=1148.667 s that `Curve fitting/v2` found independently.

Outputs land in `../standstill_tail/`.
"""
import json
import sys
from pathlib import Path

import numpy as np

V61_DIR = Path(__file__).resolve().parents[1]
MTD_DIR = V61_DIR.parent
V6_DIR = MTD_DIR / "v6"
OUT_DIR = V61_DIR / "standstill_tail"

sys.path.insert(0, str(V6_DIR / "analysis" / "scripts"))
from ids_common import (load_ids, still_mask, decimate_minmax, apply_style,
                        C_MEASURED, C_COMMANDED, C_INK, C_MUTED)

import matplotlib.pyplot as plt

GUARD_S = 1.0        # dropped from the front of the detected tail (transient)
ZOOM_S = 1.0         # width of the native-resolution excerpt panel


def find_tail_start(t, pos_um, dt_s):
    """Earliest sample after which the axis is still all the way to the end."""
    still = still_mask(pos_um, dt_s)
    moving = np.flatnonzero(~still)
    if moving.size == 0:
        raise RuntimeError("no motion detected anywhere in the capture")
    return int(moving[-1]) + 1


def summarize(t_rel, ac_nm):
    """Metrics that matter for judging a segment destined for a DFT."""
    d = np.diff(ac_nm)
    slope, intercept = np.polyfit(t_rel, ac_nm, 1)
    return {
        "n_samples": int(ac_nm.size),
        "duration_s": float(t_rel[-1] - t_rel[0]),
        "std_nm": float(ac_nm.std(ddof=1)),
        "ptp_nm": float(np.ptp(ac_nm)),
        "min_nm": float(ac_nm.min()),
        "max_nm": float(ac_nm.max()),
        "mean_abs_diff_nm": float(np.abs(d).mean()),
        "median_abs_diff_nm": float(np.median(np.abs(d))),
        "residual_linear_slope_nm_per_s": float(slope),
        "residual_linear_total_nm": float(slope * (t_rel[-1] - t_rel[0])),
        "std_after_linear_detrend_nm": float(
            (ac_nm - (slope * t_rel + intercept)).std(ddof=1)),
    }


def headroom(ax, y, frac=0.38, pad_frac=0.06):
    """Reserve empty space at the top of an axis for the legend to sit in.

    These traces fill their own extent, so every in-axes legend placement
    lands on data unless the axis is given somewhere for it to go.
    """
    lo, hi = float(np.min(y)), float(np.max(y))
    span = max(hi - lo, 1e-9)
    ax.set_ylim(lo - span * pad_frac, hi + span * frac)


def render(t, pos_um, i0_detected, i0, t_rel, ac_nm, stats, meta, out_stem):
    apply_style()
    fig, axes = plt.subplots(3, 1, figsize=(13, 9.6))
    ax_ctx, ax_tail, ax_zoom = axes

    # --- Panel 1: provenance -- where in the capture this clip comes from.
    td, yd = decimate_minmax(t, pos_um)
    ax_ctx.plot(td, yd, color=C_MEASURED, lw=0.7, label="Measured (EL5101)")
    ax_ctx.axvspan(t[i0], t[-1], color=C_COMMANDED, alpha=0.13,
                   label=f"Exported clip ({stats['duration_s']:.2f} s)")
    ax_ctx.axvline(t[i0_detected], color=C_COMMANDED, lw=1.0, ls=":",
                   label=f"Detected standstill start (t={t[i0_detected]:.3f} s)")
    ax_ctx.set_xlim(t[0], t[-1])
    ax_ctx.set_xlabel("Time from start of capture [s]")
    ax_ctx.set_ylabel("Position [µm]")
    ax_ctx.set_title("Full v6 campaign capture — the standstill tail clipped "
                     "from its end", loc="left", pad=22)
    ax_ctx.text(0.0, 1.012,
                f"EXPrun.csv · {meta['duration_s']:.1f} s at "
                f"{meta['fs_hz']:.0f} Hz · 24 trajectory legs, then nothing "
                "commanded", transform=ax_ctx.transAxes, fontsize=8.5,
                color=C_MUTED, ha="left", va="bottom")
    headroom(ax_ctx, yd, frac=0.30)
    ax_ctx.legend(loc="upper left", ncol=3)

    # --- Panel 2: the clip itself, DC removed. nm, because +/-50 nm plotted in
    # um is an axis of zeroes.
    td2, yd2 = decimate_minmax(t_rel, ac_nm)
    ax_tail.plot(td2, yd2, color=C_MEASURED, lw=0.6,
                 label="Measured, DC removed")
    ax_tail.axhline(0.0, color=C_INK, lw=0.8, alpha=0.55,
                    label="Mean of clip (subtracted DC)")
    z0 = t_rel[-1] / 2.0 - ZOOM_S / 2.0
    ax_tail.axvspan(z0, z0 + ZOOM_S, color=C_COMMANDED, alpha=0.16,
                    label=f"Panel 3 zoom ({ZOOM_S:.0f} s)")
    ax_tail.set_xlim(t_rel[0], t_rel[-1])
    ax_tail.set_xlabel("Time from start of clip [s]")
    ax_tail.set_ylabel("Position, DC removed [nm]")
    ax_tail.set_title("Standstill tail, DC offset subtracted — full clip",
                      loc="left", pad=22)
    ax_tail.text(0.0, 1.012,
                 f"σ = {stats['std_nm']:.2f} nm · peak-to-peak = "
                 f"{stats['ptp_nm']:.1f} nm · mean |Δsample| = "
                 f"{stats['mean_abs_diff_nm']:.2f} nm · residual trend = "
                 f"{stats['residual_linear_slope_nm_per_s']:+.3f} nm/s",
                 transform=ax_tail.transAxes, fontsize=8.5, color=C_MUTED,
                 ha="left", va="bottom")
    headroom(ax_tail, yd2, frac=0.34)
    ax_tail.legend(loc="upper left", ncol=3)

    # --- Panel 3: native 1 kHz resolution, so the sample-level character and
    # the 1 nm counter quantisation are both visible.
    a = int(round(z0 / meta["dt_s"]))
    b = a + int(round(ZOOM_S / meta["dt_s"]))
    ax_zoom.plot(t_rel[a:b], ac_nm[a:b], color=C_MEASURED, lw=0.8,
                 marker="o", ms=1.8, mew=0, label="Measured, DC removed")
    ax_zoom.axhline(ac_nm[a:b].mean(), color=C_INK, lw=0.8, ls="--",
                    alpha=0.55, label="Local mean of this window")
    ax_zoom.set_xlim(t_rel[a], t_rel[b - 1])
    ax_zoom.set_xlabel("Time from start of clip [s]")
    ax_zoom.set_ylabel("Position, DC removed [nm]")
    ax_zoom.set_title(f"Same signal, {ZOOM_S:.0f} s at full 1 kHz resolution "
                      "(no decimation)", loc="left", pad=22)
    ax_zoom.text(0.0, 1.012,
                 f"window σ = {ac_nm[a:b].std(ddof=1):.2f} nm · "
                 f"peak-to-peak = {np.ptp(ac_nm[a:b]):.1f} nm · "
                 "1 count = 1 nm, so the vertical steps are the encoder's own "
                 "quantisation", transform=ax_zoom.transAxes,
                 fontsize=8.5, color=C_MUTED, ha="left", va="bottom")
    headroom(ax_zoom, ac_nm[a:b], frac=0.30)
    ax_zoom.legend(loc="upper left", ncol=2)

    fig.tight_layout(h_pad=2.4)
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    src = V6_DIR / "EXPrun.csv"
    t, pos_um, meta = load_ids(src)

    i0_detected = find_tail_start(t, pos_um, meta["dt_s"])
    i0 = i0_detected + int(round(GUARD_S / meta["dt_s"]))

    t_abs = t[i0:]
    t_rel = t_abs - t_abs[0]
    raw_um = pos_um[i0:]
    dc_um = float(raw_um.mean())
    ac_um = raw_um - dc_um
    ac_nm = ac_um * 1000.0

    stats = summarize(t_rel, ac_nm)

    np.savez_compressed(
        OUT_DIR / "standstill_tail.npz",
        t_s=t_abs, t_rel_s=t_rel, pos_um=raw_um, pos_um_ac=ac_um,
        pos_nm_ac=ac_nm,
        fs_hz=meta["fs_hz"], dt_s=meta["dt_s"], dc_offset_um=dc_um,
        clip_start_s=t_abs[0], clip_end_s=t_abs[-1],
        detected_start_s=t[i0_detected], guard_s=GUARD_S,
    )

    hdr = ("# v6 EXPrun.csv standstill tail, DC removed\n"
           f"# clip {t_abs[0]:.3f}-{t_abs[-1]:.3f} s of capture, "
           f"{meta['fs_hz']:.0f} Hz, uniformly sampled, no gaps\n"
           f"# subtracted DC offset = {dc_um:.6f} um\n"
           "# t_s = time from start of capture; pos_um = raw (capture frame); "
           "pos_um_ac = DC removed\n")
    with open(OUT_DIR / "standstill_tail.csv", "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write(hdr)
        fh.write("t_s,pos_um,pos_um_ac\n")
        np.savetxt(fh, np.column_stack([t_abs, raw_um, ac_um]),
                   fmt=("%.3f", "%.3f", "%.6f"), delimiter=",")

    summary = {
        "source": str(src),
        "sample_rate_hz": meta["fs_hz"],
        "uniformly_sampled_no_gaps": True,
        "capture_duration_s": meta["duration_s"],
        "detected_standstill_start_s": float(t[i0_detected]),
        "guard_s": GUARD_S,
        "guard_rationale": ("drops the ~1.15 um settling transient in the "
                            "first ~4 ms of the detected tail, which would "
                            "otherwise smear broadband energy across the "
                            "spectrum"),
        "clip_start_s": float(t_abs[0]),
        "clip_end_s": float(t_abs[-1]),
        "dc_offset_subtracted_um": dc_um,
        "spectral_properties": {
            "nyquist_hz": meta["fs_hz"] / 2.0,
            "frequency_resolution_hz": 1.0 / stats["duration_s"],
            "counter_quantisation_nm": 1.0,
        },
        "clip_stats_dc_removed": stats,
    }
    with open(OUT_DIR / "clip_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    render(t, pos_um, i0_detected, i0, t_rel, ac_nm, stats, meta,
           str(OUT_DIR / "standstill_tail"))

    print(f"source                 {src}")
    print(f"detected standstill at t={t[i0_detected]:.3f} s")
    print(f"clip (guard {GUARD_S:.1f} s)    t={t_abs[0]:.3f}-{t_abs[-1]:.3f} s  "
          f"{stats['duration_s']:.3f} s  {stats['n_samples']} samples")
    print(f"DC offset subtracted   {dc_um:.4f} um")
    print(f"sigma                  {stats['std_nm']:.2f} nm")
    print(f"peak-to-peak           {stats['ptp_nm']:.1f} nm")
    print(f"freq resolution        {1.0 / stats['duration_s']:.5f} Hz   "
          f"Nyquist {meta['fs_hz'] / 2:.0f} Hz")
    for f in sorted(OUT_DIR.iterdir()):
        print(f"  wrote {f.name}")


if __name__ == "__main__":
    main()
