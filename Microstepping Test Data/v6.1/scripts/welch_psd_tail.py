"""Welch PSD (Hann) of the full standstill tail exported by
`clip_standstill_tail.py`.

Input is `../standstill_tail/standstill_tail.npz` -- the whole 152.65 s clip,
DC already removed. Nothing is re-clipped or re-detrended here beyond Welch's
own per-segment mean removal, so what this transforms is exactly the series in
that file.

## The configuration, and why

| | |
|---|---|
| Window | Hann |
| `nperseg` | 16384 samples = 16.384 s |
| Overlap | 50 % (8192 samples) |
| Segments averaged | 17 |
| Resolution | 0.061 Hz |
| `detrend` | `"constant"` (per segment) |

`nperseg` is held at **16384 to match
`../../v6/Curve fitting/v2/v2.2`**, which ran the same estimator over the
first 60 s of this tail. Same window, same resolution, same per-segment
detrend means the two spectra are directly comparable rather than differing
for reasons of configuration. What the full tail buys at that setting is
averaging: 17 segments against v2.2's 6, so the variance of the estimate falls
by roughly sqrt(17/6) ~ 1.7x. That is the point of running it over the whole
tail rather than a slice.

The trade `nperseg` controls is resolution against variance -- doubling it
halves Δf but halves the number of averages. 16.4 s segments resolve the
sub-Hz wander (the dominant feature at the low end) while still leaving enough
averages for the high-frequency floor to be smooth. `sweep_nperseg` prints the
trade at 2^12..2^15 so the choice can be seen rather than taken on faith.

**`detrend="constant"` removes each segment's mean, not its slope.** The tail
carries a slow wander (see `../README.md`), so the lowest bins measure that
wander, and it leaks upward through the Hann sidelobes. That is a property of
the signal, reported rather than fitted away -- v2.2 is the version that
subtracts a polynomial first, and it is the one to compare against for a
residual-noise number.

Run from this directory:

    python welch_psd_tail.py

Outputs land in `../psd/`.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

V61_DIR = Path(__file__).resolve().parents[1]
V6_DIR = V61_DIR.parent / "v6"
CLIP = V61_DIR / "standstill_tail" / "standstill_tail.npz"
OUT_DIR = V61_DIR / "psd"

sys.path.insert(0, str(V6_DIR / "analysis" / "scripts"))
from ids_common import apply_style, C_MEASURED, C_COMMANDED, C_INK, C_MUTED

import matplotlib.pyplot as plt

WINDOW = "hann"
NPERSEG = 16384
OVERLAP_FRAC = 0.5
DETREND = "constant"

FLOOR_BAND = (10.0, 250.0)      # where the spectrum is flat enough to call a floor
MODE_BAND = (330.0, 370.0)      # the peak v2/v2.2/v2.3 locate
BANDS = ((0, 1), (1, 10), (10, 100), (100, 300), (300, 400), (400, 500))

LABEL_BBOX = dict(facecolor="white", edgecolor="#d7dbe0", linewidth=0.6,
                  alpha=0.94, pad=3.0)


def n_segments(n, nperseg, noverlap):
    return 1 + (n - nperseg) // (nperseg - noverlap)


def sweep_nperseg(x, fs):
    """The resolution-vs-averaging trade, so NPERSEG is a visible choice."""
    rows = []
    for nps in (4096, 8192, 16384, 32768):
        if nps > x.size:
            continue
        rows.append({
            "nperseg": nps,
            "segment_length_s": nps / fs,
            "frequency_resolution_hz": fs / nps,
            "n_segments_averaged": int(n_segments(x.size, nps, nps // 2)),
        })
    return rows


def render(f, pxx, x, fs, nseg, stats, out_stem):
    apply_style()
    fig, (ax_psd, ax_cum) = plt.subplots(2, 1, figsize=(13, 8.6))

    cumulative = stats["_cumulative"]
    q = stats["cumulative_rms_quantiles"]
    floor_asd = stats["white_floor"]["asd_nm_per_sqrt_hz"]

    # --- Panel 1: the PSD itself.
    ax_psd.axvspan(*MODE_BAND, color=C_COMMANDED, alpha=0.12, zorder=0,
                   label=f"{MODE_BAND[0]:.0f}–{MODE_BAND[1]:.0f} Hz mode band")
    ax_psd.loglog(f[1:], pxx[1:], color=C_MEASURED, lw=0.85, zorder=3,
                  label="Welch PSD (Hann, 50 % overlap)")
    ax_psd.axhline(floor_asd ** 2, color=C_INK, lw=1.0, ls=(0, (5, 3)),
                   alpha=0.7, zorder=4,
                   label=f"Floor {FLOOR_BAND[0]:.0f}–{FLOOR_BAND[1]:.0f} Hz: "
                         f"{floor_asd:.3f} nm/√Hz")
    pk = stats["peak_above_10hz"]
    ax_psd.plot([pk["frequency_hz"]], [pk["psd_nm2_per_hz"]], marker="o",
                ms=6, color=C_COMMANDED, mec="white", mew=1.2, zorder=5)
    ax_psd.annotate(f"{pk['frequency_hz']:.1f} Hz",
                    (pk["frequency_hz"], pk["psd_nm2_per_hz"]),
                    textcoords="offset points", xytext=(-14, 12),
                    fontsize=9, color=C_COMMANDED, fontweight="semibold",
                    ha="right", bbox=LABEL_BBOX, zorder=6)
    ax_psd.set_xlim(f[1], fs / 2)
    ax_psd.set_ylabel("PSD [nm²/Hz]")
    ax_psd.set_xlabel("Frequency [Hz]")
    ax_psd.set_title("Welch power spectral density — full v6 standstill tail",
                     loc="left", pad=22)
    ax_psd.text(0.0, 1.012,
                f"Hann · nperseg {NPERSEG:,} ({NPERSEG / fs:.3f} s) · "
                f"{OVERLAP_FRAC:.0%} overlap · {nseg} segments averaged · "
                f"Δf {fs / NPERSEG:.4f} Hz · detrend '{DETREND}' per segment",
                transform=ax_psd.transAxes, fontsize=8.5, color=C_MUTED,
                ha="left", va="bottom")
    ax_psd.grid(True, which="major")
    ax_psd.grid(True, which="minor", lw=0.45, alpha=0.4)
    ax_psd.legend(loc="lower left", ncol=3)

    # --- Panel 2: where that power actually accumulates.
    sd = float(x.std(ddof=1))
    ax_cum.semilogx(f[1:], cumulative[1:], color=C_MEASURED, lw=1.8, zorder=3,
                    label="Cumulative RMS from DC")
    ax_cum.axhline(sd, color=C_INK, lw=1.0, ls=(0, (5, 3)), alpha=0.7,
                   zorder=2, label=f"Time-domain σ: {sd:.2f} nm")
    # The 50 % point lands in the lowest resolvable bin, so it is reported as a
    # bound rather than marked at a frequency the estimate cannot support; 90 %
    # and 99 % are resolvable and are what get marked.
    for frac in ("0.90", "0.99"):
        fq = q[frac]["frequency_hz"]
        ax_cum.axvline(fq, color=C_COMMANDED, lw=1.0, ls=(0, (3, 3)), zorder=2)
        ax_cum.plot([fq], [float(frac) * cumulative[-1]], marker="o", ms=6,
                    color=C_COMMANDED, mec="white", mew=1.2, zorder=4)
        ax_cum.annotate(f"{float(frac):.0%} of RMS by {fq:.3g} Hz",
                        (fq, float(frac) * cumulative[-1]),
                        textcoords="offset points", xytext=(12, -22),
                        fontsize=9, color=C_COMMANDED, fontweight="semibold",
                        bbox=LABEL_BBOX, zorder=5)
    ax_cum.annotate(f"{q['0.50']['power_fraction_in_first_two_bins']:.0%} of "
                    f"the power is already in the two lowest bins\n"
                    f"(≤ {f[1]:.3f} Hz) — the wander is below this estimate's "
                    f"resolution",
                    (f[1], cumulative[1]), textcoords="offset points",
                    xytext=(14, 62), fontsize=8.5, color=C_MUTED,
                    bbox=LABEL_BBOX, zorder=5)
    ax_cum.set_xlim(f[1], fs / 2)
    ax_cum.set_ylim(0, max(cumulative[-1], sd) * 1.32)
    ax_cum.set_xlabel("Frequency [Hz]")
    ax_cum.set_ylabel("Cumulative RMS [nm]")
    ax_cum.set_title("Cumulative RMS — which frequencies carry the motion",
                     loc="left", pad=22)
    ax_cum.text(0.0, 1.012,
                f"integrated {stats['integrated_rms_nm']:.2f} nm against a "
                f"time-domain σ of {sd:.2f} nm "
                f"(Parseval ratio {stats['parseval_ratio']:.3f})",
                transform=ax_cum.transAxes, fontsize=8.5, color=C_MUTED,
                ha="left", va="bottom")
    ax_cum.grid(True, which="major")
    ax_cum.grid(True, which="minor", axis="x", lw=0.45, alpha=0.4)
    ax_cum.legend(loc="upper left", ncol=2)

    fig.tight_layout(h_pad=2.6)
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    d = np.load(CLIP)
    x = d["pos_nm_ac"]                  # DC-removed, nm
    fs = float(d["fs_hz"])
    noverlap = int(NPERSEG * OVERLAP_FRAC)
    nseg = int(n_segments(x.size, NPERSEG, noverlap))

    f, pxx = signal.welch(x, fs=fs, window=WINDOW, nperseg=NPERSEG,
                          noverlap=noverlap, detrend=DETREND,
                          return_onesided=True, scaling="density")
    df = float(f[1] - f[0])

    cumulative = np.sqrt(np.cumsum(pxx) * df)
    integrated = float(cumulative[-1])

    # Where the RMS budget fills up. Reported as the first bin at or below
    # which the fraction is reached -- NOT interpolated between bins: the 50 %
    # point falls inside the lowest bins, where an interpolated "frequency"
    # would be a number the estimate cannot resolve rather than a measurement.
    quantiles = {}
    for frac in (0.5, 0.9, 0.95, 0.99):
        i = int(np.searchsorted(cumulative, frac * integrated))
        quantiles[f"{frac:.2f}"] = {
            "frequency_hz": float(f[i]),
            "bin_index": i,
            "resolved": bool(i > 1),
            "power_fraction_in_first_two_bins":
                float(np.sum(pxx[:2]) * df / integrated ** 2),
        }

    floor_m = (f >= FLOOR_BAND[0]) & (f < FLOOR_BAND[1])
    hf = f > 10.0                       # above the wander, where mechanics live
    pk = int(np.flatnonzero(hf)[int(np.argmax(pxx[hf]))])
    wander_pk = int(np.argmax(pxx[1:])) + 1
    mode_m = (f >= MODE_BAND[0]) & (f < MODE_BAND[1])

    stats = {
        "integrated_rms_nm": integrated,
        "parseval_ratio": integrated / float(x.std(ddof=1)),
        "parseval_note": ("below 1 because detrend='constant' removes each "
                          "16.384 s segment's own mean, which discards the "
                          "wander between segments -- i.e. content below "
                          "~1/nperseg. The missing variance is that wander, "
                          "not an estimator error"),
        "cumulative_rms_quantiles": quantiles,
        "peak_above_10hz": {
            "frequency_hz": float(f[pk]),
            "psd_nm2_per_hz": float(pxx[pk]),
            "times_above_floor": float(pxx[pk] / np.median(pxx[floor_m])),
        },
        "wander_peak": {
            "frequency_hz": float(f[wander_pk]),
            "psd_nm2_per_hz": float(pxx[wander_pk]),
        },
        "white_floor": {
            "band_hz": list(FLOOR_BAND),
            "median_psd_nm2_per_hz": float(np.median(pxx[floor_m])),
            "asd_nm_per_sqrt_hz": float(np.sqrt(np.median(pxx[floor_m]))),
        },
        "mode_band_rms_nm": float(np.sqrt(np.sum(pxx[mode_m]) * df)),
        "band_rms_nm": {
            f"{lo:g}-{hi:g}_Hz":
                float(np.sqrt(np.sum(pxx[(f >= lo) & (f < hi)]) * df))
            for lo, hi in BANDS
        },
        "_cumulative": cumulative,
    }

    render(f, pxx, x, fs, nseg, stats, str(OUT_DIR / "welch_psd_tail"))

    np.savez_compressed(OUT_DIR / "welch_psd_tail.npz", f_hz=f,
                        psd_nm2_per_hz=pxx, cumulative_rms_nm=cumulative,
                        fs_hz=fs, nperseg=NPERSEG, noverlap=noverlap,
                        n_segments=nseg, window=WINDOW, detrend=DETREND)

    stats.pop("_cumulative")
    summary = {
        "source_clip": str(CLIP),
        "input": {
            "series": "pos_nm_ac (DC already removed by clip_standstill_tail.py)",
            "n_samples": int(x.size),
            "duration_s": float(x.size / fs),
            "sample_rate_hz": fs,
            "time_domain_std_nm": float(x.std(ddof=1)),
        },
        "welch": {
            "window": WINDOW,
            "nperseg": NPERSEG,
            "segment_length_s": NPERSEG / fs,
            "overlap_fraction": OVERLAP_FRAC,
            "noverlap": noverlap,
            "n_segments_averaged": nseg,
            "frequency_resolution_hz": fs / NPERSEG,
            "nyquist_hz": fs / 2,
            "detrend_per_segment": DETREND,
            "scaling": "density",
            "nperseg_rationale": ("matched to v6/Curve fitting/v2/v2.2 so the "
                                  "two spectra are directly comparable; the "
                                  "full tail buys 17 averages against its 6"),
            "nperseg_trade": sweep_nperseg(x, fs),
        },
        "results": stats,
        "caveat_aliasing": ("the spectrum rises into Nyquist, the signature of "
                            "folding: a feature at f is indistinguishable from "
                            "1000-f at this 1 kHz rate. See ../README.md"),
        "outputs": {
            "plot": "psd/welch_psd_tail.png",
            "plot_vector": "psd/welch_psd_tail.svg",
            "spectrum": "psd/welch_psd_tail.npz",
        },
    }
    with open(OUT_DIR / "welch_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"input            {x.size:,} samples, {x.size / fs:.3f} s, "
          f"σ = {x.std(ddof=1):.2f} nm")
    print(f"welch            {WINDOW}, nperseg {NPERSEG:,} "
          f"({NPERSEG / fs:.3f} s), {OVERLAP_FRAC:.0%} overlap")
    print(f"                 {nseg} segments averaged, Δf = {fs / NPERSEG:.4f} Hz, "
          f"Nyquist {fs / 2:.0f} Hz")
    print(f"integrated RMS   {integrated:.2f} nm "
          f"(Parseval ratio {stats['parseval_ratio']:.3f})")
    q50 = stats["cumulative_rms_quantiles"]["0.50"]
    print(f"50% of RMS       by bin {q50['bin_index']} "
          f"({q50['frequency_hz']:.4f} Hz)"
          f"{'' if q50['resolved'] else '  [unresolved: lowest bins]'}")
    print(f"90% / 99% of RMS "
          f"{stats['cumulative_rms_quantiles']['0.90']['frequency_hz']:.3f} Hz"
          f" / {stats['cumulative_rms_quantiles']['0.99']['frequency_hz']:.1f} Hz")
    print(f"wander peak      {stats['wander_peak']['frequency_hz']:.4f} Hz")
    print(f"peak above 10 Hz {stats['peak_above_10hz']['frequency_hz']:.2f} Hz "
          f"({stats['peak_above_10hz']['times_above_floor']:.0f}x the floor)")
    print(f"floor {FLOOR_BAND[0]:.0f}-{FLOOR_BAND[1]:.0f} Hz "
          f"{stats['white_floor']['asd_nm_per_sqrt_hz']:.4f} nm/√Hz")
    for f_ in sorted(OUT_DIR.iterdir()):
        print(f"  wrote {f_.name}")


if __name__ == "__main__":
    main()
