"""
Process a Dewesoft chirp-test recording (tempTVAl.csv, time + 3 accel
channels, no embedded trigger) into Bode magnitude plots.

The firmware (ESP32S3_TMC2209_Chirp_Test.ino) commands a bounded
alternating-microstep excitation whose frequency is a known linear law:
  2 s lead-in (stationary) -> up-sweep 1..1000 Hz over 30 s
  -> 0.5 s turnaround (stationary) -> down-sweep 1000..1 Hz over 30 s
  -> 2 s tail (stationary)
with STEP commands suppressed for 120-230 Hz (resonance exclusion notch).

The CSV has no trigger channel and the recording brackets the sequence with
extra idle padding, so the exact start time of each sweep is recovered from
the data itself:
  1. Locate the two ~30 s high-energy blocks via a short-time RMS envelope.
  2. Within each block, estimate the instantaneous excitation frequency via
     ridge-picking on a windowed FFT, then fit t -> f as a line (excluding
     the notch band) to recover the precise start time/rate.
  3. Using that fitted law, evaluate the response magnitude at each
     command frequency via a synchronous (single-frequency, Hann-windowed)
     correlation -- i.e. a swept-sine analysis -- for both the up and the
     down sweep.

No pandas/scipy dependency: numpy's C-accelerated loadtxt is fast enough
for this file size, and the STFT/synchronous-detection is plain FFT/numpy.
"""

import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_PATH = "/Users/tomasvalentinas/Documents/System-Mod-Sim/ESP32S3_TMC2209_Chirp_Test/tempTVAl.csv"
PLOT_DIR = "/Users/tomasvalentinas/Documents/System-Mod-Sim/ESP32S3_TMC2209_Chirp_Test/analysis/plots"

# Column order assumed [time, AI0, AI1, AI2]; AI1 (index 2) is the primary/x axis.
COL_TIME = 0
COL_AI0 = 1
COL_AI1 = 2
COL_AI2 = 3
CHANNEL_NAMES = {COL_AI0: "AI0", COL_AI1: "AI1 (x, primary)", COL_AI2: "AI2"}

NOTCH_LOW_HZ, NOTCH_HIGH_HZ = 120.0, 230.0
F_LO, F_HI = 1.0, 1000.0


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_csv(path):
    log(f"loading {path} ...")
    t0 = time.time()
    arr = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float32)
    log(f"loaded shape={arr.shape} in {time.time()-t0:.1f}s")
    return arr


def sample_rate(t):
    dt = np.median(np.diff(t[:20000]).astype(np.float64))
    fs = 1.0 / dt
    log(f"dt={dt:.9f}s -> fs={fs:.3f} Hz (checked first {20000} samples)")
    return fs, dt


def ridge_frequencies(x, fs, t_lo, t_hi, win_s=0.5, hop_s=0.1, f_search=(0.5, 1050.0)):
    n0 = max(0, int(t_lo * fs))
    n1 = min(len(x), int(t_hi * fs))
    win = int(round(win_s * fs))
    hop = int(round(hop_s * fs))
    hann = np.hanning(win)
    freqs = np.fft.rfftfreq(win, d=1.0 / fs)
    band = (freqs >= f_search[0]) & (freqs <= f_search[1])
    band_freqs = freqs[band]

    times, ridge = [], []
    for start in range(n0, n1 - win, hop):
        seg = x[start:start + win].astype(np.float64)
        seg = seg - seg.mean()
        spec = np.abs(np.fft.rfft(seg * hann))[band]
        if spec.size == 0:
            continue
        peak = band_freqs[np.argmax(spec)]
        times.append((start + win / 2.0) / fs)
        ridge.append(peak)
    return np.array(times), np.array(ridge)


def robust_line_fit(times, ridge, t_lo, t_hi, f_lo, f_hi, iters=10, min_points=20):
    """Iterative sigma-clipped least-squares fit of ridge ~ rate*t + b,
    seeded from a (t_lo,t_hi) x (f_lo,f_hi) window and refined by rejecting
    outliers each pass (ambient tones, notch-band dropout, harmonics)."""
    mask = (times >= t_lo) & (times <= t_hi) & (ridge >= f_lo) & (ridge <= f_hi)
    rate = b = None
    for i in range(iters):
        tt, ff = times[mask], ridge[mask]
        if tt.size < min_points:
            raise RuntimeError(f"robust_line_fit: too few points ({tt.size}) to continue")
        A = np.vstack([tt, np.ones_like(tt)]).T
        rate, b = np.linalg.lstsq(A, ff, rcond=None)[0]
        pred = rate * times + b
        resid = ridge - pred
        sigma = np.std(resid[mask])
        newmask = (times >= t_lo) & (times <= t_hi) & (np.abs(resid) < max(3 * sigma, 6.0))
        log(f"    iter{i}: rate={rate:+.4f} Hz/s  b={b:+.2f}  sigma={sigma:.2f}  n={newmask.sum()}")
        if newmask.sum() == mask.sum():
            mask = newmask
            break
        mask = newmask
    return rate, b, mask.sum()


def locate_and_fit_sweeps(x, fs, total_s):
    """Find the up/down sweep linear frequency laws f(t)=rate*t+b directly
    from a full-record spectral ridge trace (argmax per STFT frame). No
    trigger channel is available, and raw-energy (RMS) localization is not
    reliable here because the excitation is weak relative to ambient/motor
    noise -- but the swept ridge is still the dominant spectral peak once it
    rises above roughly the notch band, which is enough to calibrate the
    exact timing law; low-frequency points are then recovered later via
    synchronous detection at the extrapolated known frequency, not by
    peak-picking."""
    log("computing full-record ridge trace for sweep-law calibration ...")
    times, ridge = ridge_frequencies(x, fs, 0.0, total_s)
    log(f"  {len(times)} ridge frames computed")

    near_top = ridge > 0.95 * ridge.max()
    t_peak = float(np.median(times[near_top]))
    log(f"  approximate up/down transition (near-max ridge) at t={t_peak:.2f}s")

    # Fit DOWN first: it is reliably dominant across a wide window on its own.
    log("  fitting DOWN sweep:")
    rate_dn, b_dn, n_dn = robust_line_fit(
        times, ridge, t_lo=max(0.0, t_peak - 3.0), t_hi=total_s, f_lo=1.0, f_hi=1050.0)

    # Use DOWN's rate (sweep is symmetric in duration) to predict when the UP
    # sweep's instantaneous frequency first reaches the search band's floor,
    # so the UP fit window skips the early stretch dominated by a persistent
    # ambient/motor tone below that frequency (the raw ridge is otherwise
    # pulled toward that flat tone instead of the true rising chirp).
    f_lo_up = NOTCH_HIGH_HZ - 10
    expected_duration = (F_HI - F_LO) / abs(rate_dn)
    t_start_up_est = t_peak - expected_duration
    t_lo_up = t_start_up_est + (f_lo_up - F_LO) / abs(rate_dn)
    log(f"  UP search window informed by DOWN: expected_duration={expected_duration:.2f}s, "
        f"predicted UP start={t_start_up_est:.2f}s, t_lo_up={t_lo_up:.2f}s")

    log("  fitting UP sweep (restricted to above-notch band to dodge ambient tones):")
    rate_up, b_up, n_up = robust_line_fit(
        times, ridge, t_lo=t_lo_up, t_hi=t_peak + 3.0, f_lo=f_lo_up, f_hi=1050.0)

    for label, rate, b in [("UP", rate_up, b_up), ("DOWN", rate_dn, b_dn)]:
        t_at_1 = (1.0 - b) / rate
        t_at_1000 = (1000.0 - b) / rate
        log(f"  {label}: rate={rate:+.4f} Hz/s, t(f=1)={t_at_1:.2f}s, "
            f"t(f=1000)={t_at_1000:.2f}s, implied duration={abs(t_at_1000-t_at_1):.2f}s")

    return (rate_up, b_up), (rate_dn, b_dn)


def synchronous_magnitude(x, fs, t_center, freq, cycles=10, win_min_s=0.01, win_max_s=0.6):
    win_s = np.clip(cycles / freq, win_min_s, win_max_s)
    n = int(round(win_s * fs))
    if n < 8:
        n = 8
    start = int(round(t_center * fs)) - n // 2
    stop = start + n
    if start < 0 or stop > len(x):
        return np.nan
    seg = x[start:stop].astype(np.float64)
    seg = seg - seg.mean()
    w = np.hanning(n)
    tt = np.arange(n) / fs
    c = np.cos(2 * np.pi * freq * tt)
    s = np.sin(2 * np.pi * freq * tt)
    re = np.sum(seg * w * c)
    im = np.sum(seg * w * s)
    coherent_gain = w.sum()
    mag = 2.0 * np.hypot(re, im) / coherent_gain
    return mag


def sweep_magnitude_curve(x, fs, rate, b, freq_grid):
    """f(t) = rate*t + b  =>  t(f) = (f - b) / rate."""
    mags = np.empty_like(freq_grid)
    for i, f in enumerate(freq_grid):
        t_center = (f - b) / rate
        mags[i] = synchronous_magnitude(x, fs, t_center, f)
    return mags


def mask_notch(freq_grid, values):
    out = values.copy()
    out[(freq_grid >= NOTCH_LOW_HZ) & (freq_grid <= NOTCH_HIGH_HZ)] = np.nan
    return out


def make_bode_plots(freq_grid, up_mag, down_mag, out_prefix, title_prefix):
    for norm, norm_label, norm_suffix in [
        (False, "raw", "raw"),
        (True, "normalized by f^2", "f2norm"),
    ]:
        if norm:
            f2 = freq_grid ** 2
            up_y = up_mag / f2
            down_y = down_mag / f2
            ylabel = "acceleration magnitude / f^2  (m/s^2 per Hz^2)"
        else:
            up_y = up_mag
            down_y = down_mag
            ylabel = "acceleration magnitude (m/s^2)"

        for yscale, scale_suffix in [("linear", "linY"), ("log", "logY")]:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(freq_grid, mask_notch(freq_grid, up_y), label="up sweep (1->1000 Hz)",
                    lw=1.2)
            ax.plot(freq_grid, mask_notch(freq_grid, down_y), label="down sweep (1000->1 Hz)",
                    lw=1.2, alpha=0.8)
            ax.axvspan(NOTCH_LOW_HZ, NOTCH_HIGH_HZ, color="gray", alpha=0.15,
                       label="resonance notch (unexcited)")
            ax.set_xlabel("frequency (Hz)")
            ax.set_ylabel(ylabel)
            ax.set_yscale(yscale)
            ax.set_xlim(F_LO, F_HI)
            ax.set_title(f"{title_prefix} -- {norm_label}, {yscale} y-scale")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend()
            fig.tight_layout()
            out_path = f"{PLOT_DIR}/{out_prefix}_{norm_suffix}_{scale_suffix}.png"
            fig.savefig(out_path, dpi=150)
            plt.close(fig)
            log(f"  saved {out_path}")


def make_offaxis_overview(freq_grid, mags_by_channel, out_path):
    fig, ax = plt.subplots(figsize=(10, 6))
    for col, mag in mags_by_channel.items():
        ax.plot(freq_grid, mask_notch(freq_grid, mag), label=CHANNEL_NAMES[col], lw=1.1)
    ax.axvspan(NOTCH_LOW_HZ, NOTCH_HIGH_HZ, color="gray", alpha=0.15,
               label="resonance notch (unexcited)")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("acceleration magnitude (m/s^2)")
    ax.set_yscale("log")
    ax.set_xlim(F_LO, F_HI)
    ax.set_title("All-axis overview (up sweep, raw, log y-scale)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log(f"  saved {out_path}")


def main():
    arr = load_csv(CSV_PATH)
    t = arr[:, COL_TIME].astype(np.float64)
    fs, dt = sample_rate(t)
    total_s = t[-1] - t[0]
    log(f"total recorded duration: {total_s:.2f}s ({len(t)} samples)")

    ai1 = arr[:, COL_AI1]
    (rate_up, b_up), (rate_dn, b_dn) = locate_and_fit_sweeps(ai1, fs, total_s)

    freq_grid = np.linspace(F_LO, F_HI, 1000)

    log("computing synchronous-detection magnitude curves ...")
    results = {}
    for col in (COL_AI1, COL_AI0, COL_AI2):
        ch = arr[:, col]
        t0s = time.time()
        up_mag = sweep_magnitude_curve(ch, fs, rate_up, b_up, freq_grid)
        down_mag = sweep_magnitude_curve(ch, fs, rate_dn, b_dn, freq_grid)
        results[col] = (up_mag, down_mag)
        log(f"  {CHANNEL_NAMES[col]} done in {time.time()-t0s:.1f}s "
            f"(up: {np.isfinite(up_mag).sum()}/{len(freq_grid)} valid, "
            f"down: {np.isfinite(down_mag).sum()}/{len(freq_grid)} valid)")

    up_mag1, down_mag1 = results[COL_AI1]
    make_bode_plots(freq_grid, up_mag1, down_mag1, "bode_AI1", "Bode: AI1 (x, primary)")

    make_offaxis_overview(
        freq_grid,
        {col: results[col][0] for col in (COL_AI0, COL_AI1, COL_AI2)},
        f"{PLOT_DIR}/bode_all_axes_overview.png",
    )

    np.savez(
        f"{PLOT_DIR}/bode_data.npz",
        freq_grid=freq_grid,
        up_mag_AI0=results[COL_AI0][0], down_mag_AI0=results[COL_AI0][1],
        up_mag_AI1=results[COL_AI1][0], down_mag_AI1=results[COL_AI1][1],
        up_mag_AI2=results[COL_AI2][0], down_mag_AI2=results[COL_AI2][1],
        rate_up=rate_up, b_up=b_up, rate_down=rate_dn, b_down=b_dn,
        fs=fs,
    )
    log("saved bode_data.npz")
    log("done.")


if __name__ == "__main__":
    main()
