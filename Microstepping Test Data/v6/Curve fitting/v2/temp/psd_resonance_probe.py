"""Indicative PSD probe: find the mechanical resonance in the idle tail and
derive the Gaussian filter parameters from it.

    [ raw 10 s block ] -> [ mean + linear detrend ] -> [ Welch PSD ] -> [ peak fc ]

Then, per the specified relations at fs = 1000 Hz:

    T       = 1 / fc                 noise wave period
    W       = 2 * T * fs             window spans exactly two full cycles
    sigma   = W / 6

Two renders:
  1. psd_detrended.png   -- the PSD itself, peak annotated, noise floor visible
  2. gaussian_overlay.png -- 0.1 s zoom, raw vs Gaussian filtered with the
                             parameters the PSD peak implies

This is exploratory: it exists to check whether the ~140-150 Hz resonance
assumed when the filter was first parameterised is actually in the data.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch, detrend
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
sys.path.insert(0, str(HERE.parent.parent.parent / "analysis" / "scripts"))
from ids_common import load_ids, apply_style, C_MEASURED, C_COMMANDED, C_INK, C_MUTED, C_GRID
from plot_v6_exp_run import CSV
from analyze_post_recording_tail import find_post_recording_tail, window_start_level

BLOCK_S = 10.0          # length of the steady uncommanded block to analyse
NPERSEG = 2048          # Welch segment length -> fs/NPERSEG ~= 0.49 Hz resolution
PEAK_SEARCH_LO_HZ = 5.0  # ignore the drift/wander end when hunting for the peak
ZOOM_S = 0.1            # width of the time-domain verification window

# Window width in samples. None = take it from the PSD peak (2 cycles of fc,
# rounded to the nearest odd integer). Set explicitly to try a wider kernel
# than the ripple period alone would justify.
WINDOW_SAMPLES_OVERRIDE = 10


def style(ax):
    ax.grid(True, color=C_GRID, alpha=.7, linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def main():
    t, pos, meta = load_ids(CSV)
    dt_s, fs = meta["dt_s"], meta["fs_hz"]
    start_idx, end_idx, _ = find_post_recording_tail(t, pos, dt_s)

    # A steady block from the middle of the tail: no commanded motion, and far
    # enough from the boundary transient at the tail's own start.
    seg_t, seg = t[start_idx:end_idx], pos[start_idx:end_idx]
    center = 0.5 * (seg_t[0] + seg_t[-1])
    m = (seg_t >= center - BLOCK_S / 2) & (seg_t <= center + BLOCK_S / 2)
    bt, block_um = seg_t[m], seg[m]
    block_nm = (block_um - window_start_level(block_um)) * 1000.0

    # Step 1: kill the mean and the linear (thermal) trend so the macro wander
    # can't leak into the high-frequency part of the spectrum.
    block_dt = detrend(block_nm, type="linear")

    print(f"block: {bt[0]:.3f}-{bt[-1]:.3f} s, n={block_nm.size}, fs={fs:.0f} Hz")
    print(f"  raw std {block_nm.std(ddof=1):.2f} nm -> detrended std {block_dt.std(ddof=1):.2f} nm")

    # Step 2: Welch PSD -- averaging many windowed FFTs flattens the random
    # background so a genuine mechanical line stands out from it.
    f, pxx = welch(block_dt, fs=fs, nperseg=NPERSEG, noverlap=NPERSEG // 2,
                   window="hann", detrend="linear", scaling="density")

    # Step 3: locate the peak above the drift-dominated low end.
    band = f >= PEAK_SEARCH_LO_HZ
    fb, pb = f[band], pxx[band]
    fc = float(fb[np.argmax(pb)])
    peak_power = float(pb.max())

    # Characterise the floor so we can say how prominent the peak really is.
    floor = float(np.median(pb))
    prominence = peak_power / floor

    order = np.argsort(pb)[::-1][:8]
    print(f"\n  peak fc = {fc:.2f} Hz   power {peak_power:.4g} nm^2/Hz")
    print(f"  median floor {floor:.4g} nm^2/Hz  -> peak/floor = {prominence:.2f}x")
    print("  top bins:")
    for i in order:
        print(f"    {fb[i]:8.2f} Hz   {pb[i]:.4g}   ({pb[i]/floor:.2f}x floor)")

    # Is the peak a property of the machine or of this one block? Sweep
    # 10 s blocks across the whole tail and see whether fc holds still.
    print("\n  stability sweep (10 s blocks across the tail):")
    fcs = []
    for t0 in np.arange(seg_t[0] + 1, seg_t[-1] - BLOCK_S, 15.0):
        bm = (seg_t >= t0) & (seg_t < t0 + BLOCK_S)
        yb = detrend((seg[bm] - seg[bm][0]) * 1000.0, type="linear")
        fb2, pb2 = welch(yb, fs=fs, nperseg=NPERSEG, noverlap=NPERSEG // 2,
                         window="hann", detrend="linear")
        k = fb2 >= PEAK_SEARCH_LO_HZ
        fpk = float(fb2[k][np.argmax(pb2[k])])
        fcs.append(fpk)
        print(f"    t0={t0:8.1f} s   fc={fpk:7.1f} Hz   "
              f"peak/floor={pb2[k].max()/np.median(pb2[k]):6.1f}x")
    fcs = np.array(fcs)
    core = fcs[fcs < 450]           # drop blocks that latch onto the Nyquist rise
    print(f"    -> {core.size}/{fcs.size} blocks agree: "
          f"fc = {core.mean():.1f} +/- {core.std(ddof=1):.1f} Hz")

    # Where does the power actually sit? Contrast the band the peak occupies
    # with the 120-160 Hz band the filter was originally parameterised for.
    print("\n  band power [nm^2]:")
    for lo, hi in [(5, 100), (100, 200), (120, 160), (200, 300),
                   (300, 400), (400, 480), (480, 500)]:
        bb = (f >= lo) & (f < hi)
        print(f"    {lo:3d}-{hi:3d} Hz : {np.trapezoid(pxx[bb], f[bb]):9.4f}")

    # Step 4: frequency -> Gaussian kernel geometry. The window spans two full
    # cycles of the ripple, rounded to an odd integer so the kernel is centred
    # on a sample; sigma then follows the 3-sigma rule (+/-3 sigma covers the
    # window, so the bell tapers to ~zero at its bounds rather than being
    # chopped off mid-slope).
    T = 1.0 / fc
    W_exact = 2.0 * T * fs
    if WINDOW_SAMPLES_OVERRIDE is None:
        W = float(int(round((W_exact - 1) / 2.0) * 2 + 1))   # nearest odd integer
        w_note = "nearest odd integer of 2*T*fs"
    else:
        W = float(WINDOW_SAMPLES_OVERRIDE)
        w_note = f"set explicitly ({W/W_exact:.2f}x the 2-cycle width)"
    sigma = W / 6.0
    print(f"\n  T = 1/fc            = {T*1000:.3f} ms")
    print(f"  2*T*fs              = {W_exact:.2f} samples -> W = {W:g} ({w_note})")
    print(f"  sigma = W/6         = {sigma:.3f} samples  (3-sigma rule)")
    print(f"  kernel spans        = {W/ (T*fs):.2f} cycles of {fc:.1f} Hz "
          f"({W/fs*1000:.2f} ms)")

    # What that kernel actually does to the ripple it was built for.
    sigma_s = sigma / fs
    atten = np.exp(-2 * np.pi**2 * sigma_s**2 * fc**2)
    f3db = np.sqrt(np.log(2) / (2 * np.pi**2 * sigma_s**2))
    print(f"  -> -3dB cutoff      = {f3db:.1f} Hz")
    print(f"  -> gain at fc       = {atten:.3f}  ({100*(1-atten):.1f}% of the "
          f"ripple removed)")

    # truncate=3.0 IS the 3-sigma rule: the kernel is cut at +/-3 sigma = +/-W/2,
    # i.e. it spans exactly the target window and tapers to ~zero at its bounds.
    filtered = gaussian_filter1d(block_dt, sigma=sigma, truncate=3.0, mode="nearest")
    radius = int(3.0 * sigma + 0.5)
    print(f"  realised kernel     = {2*radius+1} taps (+/-{radius} samples)")
    # Did it actually take the 300-400 Hz band out?
    _, pxx_f = welch(filtered, fs=fs, nperseg=NPERSEG, noverlap=NPERSEG // 2,
                     window="hann", detrend="linear")
    bb = (f >= 300) & (f < 400)
    print(f"  -> 300-400 Hz power {np.trapezoid(pxx[bb], f[bb]):.3f} -> "
          f"{np.trapezoid(pxx_f[bb], f[bb]):.3f} nm^2")
    print(f"  filtered std        = {filtered.std(ddof=1):.2f} nm "
          f"(from {block_dt.std(ddof=1):.2f} nm detrended)")
    print(f"  samples per cycle   = {fs/fc:.2f}   "
          f"(Nyquist limit is 2.00 -- see the caveat below)")
    print("\n  CAVEAT: fs/2 = 500 Hz, so anything above 500 Hz in the real "
          "machine folds\n  back into this band. A measured 365 Hz is "
          "indistinguishable from 635 Hz (or\n  1365 Hz, ...) at this sample "
          "rate, and the PSD's steep rise toward 500 Hz\n  is the usual "
          "signature of that folding. 365 Hz is the right number for "
          "designing\n  a filter that runs on THIS 1 kHz series; it is not "
          "proof of a 365 Hz mode.")

    apply_style()

    # --- render 1: the PSD -------------------------------------------------
    fig, (ax_log, ax_lin) = plt.subplots(2, 1, figsize=(12, 8.5))

    for ax, logy in [(ax_log, True), (ax_lin, False)]:
        ax.plot(f, pxx, color=C_MEASURED, lw=.9)
        ax.axvline(fc, color=C_COMMANDED, lw=1.2, ls="--", alpha=.9,
                   label=f"peak  $f_c$ = {fc:.1f} Hz")
        ax.axhline(floor, color=C_MUTED, lw=.9, ls=":", alpha=.8,
                   label=f"median floor ({floor:.3g} nm²/Hz)")
        ax.axvspan(0, PEAK_SEARCH_LO_HZ, color=C_INK, alpha=.06, lw=0)
        style(ax)
        ax.set_xlabel("Frequency [Hz]", fontsize=9)
        ax.set_ylabel("PSD [nm²/Hz]", fontsize=9)
        ax.set_xlim(0, fs / 2)
        if logy:
            ax.set_yscale("log")
        ax.legend(loc="upper right", fontsize=8, framealpha=.92)

    ax_log.set_title(f"Detrended PSD (Welch, nperseg={NPERSEG}, Hann, 50% overlap) — "
                     f"log scale  ·  block {bt[0]:.1f}-{bt[-1]:.1f} s  ·  "
                     f"peak/floor = {prominence:.1f}×", loc="left", fontsize=9.5)
    ax_lin.set_title("same, linear scale (shaded = below the peak-search cutoff, "
                     "drift-dominated)", loc="left", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "psd_detrended.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\nwrote", HERE / "psd_detrended.png")

    # --- render 2: time-domain verification --------------------------------
    zc = 0.5 * (bt[0] + bt[-1])
    zm = (bt >= zc - ZOOM_S / 2) & (bt <= zc + ZOOM_S / 2)

    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.plot(bt[zm], block_dt[zm], color=C_MEASURED, lw=1.0, marker=".", ms=3.5,
            alpha=.75, label=f"raw (detrended)  ·  std = {block_dt[zm].std(ddof=1):.1f} nm")
    ax.plot(bt[zm], filtered[zm], color=C_COMMANDED, lw=2.0,
            label=f"Gaussian W={W:.1f}, $\\sigma$={sigma:.2f} smp  ·  "
                  f"std = {filtered[zm].std(ddof=1):.1f} nm")
    style(ax)
    ax.set_xlabel("Time from start of capture [s]", fontsize=9)
    ax.set_ylabel("Detrended position [nm]", fontsize=9)
    ax.set_xlim(bt[zm][0], bt[zm][-1])
    ax.set_title(f"{ZOOM_S:g} s verification window  ·  kernel spans "
                 f"{W/(T*fs):.2f} cycles of {fc:.1f} Hz ({W/fs*1000:.1f} ms)  ·  "
                 f"{100*(1-atten):.1f}% of $f_c$ removed", loc="left", fontsize=9.5)
    ax.legend(loc="upper right", fontsize=8, framealpha=.92)
    fig.tight_layout()
    fig.savefig(HERE / "gaussian_overlay.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", HERE / "gaussian_overlay.png")


if __name__ == "__main__":
    main()
