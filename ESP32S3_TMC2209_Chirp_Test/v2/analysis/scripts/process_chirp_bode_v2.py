"""Process the v2 chirp recording (Ftest.csv) into Bode magnitude plots.

Same deliverable as the v1 analysis (analysis/scripts/process_chirp_bode.py):
per-axis Bode magnitude plots, an all-axis overview (f^2-normalized, log y)
and a .npz of the curves. Raw magnitude is rendered on both linear and log y;
the f^2-normalized and transmissibility views span several decades and are
produced on log y only. Three things about the v2 design make the v1 *method*
inapplicable, though, and each is handled differently here:

1. **Timing comes from the sync markers, not a ridge fit.** v1 had no trigger
   channel, so it recovered the sweep law by peak-picking a spectral ridge and
   fitting a line to it. v2 emits two 750 Hz marker bursts (3 x 100 ms, 100 ms
   gaps) at known sequence times 30.0 s and 275.5 s. Locating those pins the
   sequence origin directly, and their measured separation cross-checks the
   DAQ-vs-ESP32 clock ratio. This is exact where v1 was a fit.

2. **The frequency law is piecewise, not linear.** v1 swept 1->1000 Hz
   linearly. v2 sweeps log 1->60 Hz over 120 s, then linear 60->1000 Hz over
   120 s. Both branches are inverted analytically here (t(f)), so no fitting is
   involved at all.

3. **The notch is a taper, not an exclusion.** v1 suppressed STEP commands over
   120-230 Hz, so that band had no excitation and had to be masked to NaN. v2
   tapers amplitude to 10% at f_n but never stops commanding, so the band is
   real data and is NOT masked -- it is shaded to mark reduced drive.

That third point generalizes into the main analysis upgrade over v1. Because
the commanded amplitude is known in closed form at every frequency, the
response can be divided by it to give a genuine transfer function rather than a
bare response curve. This matters more than it sounds: the step-rate clamp
(200 kHz / (2*pi*f*MRES)) starts binding at 663 Hz and rolls the commanded
amplitude off as 1/f above there (3.00 -> 1.99 full steps by 1 kHz), so a raw
plot shows a ~34% high-frequency droop that belongs to the *excitation*, not to
the mechanism. The `transmissibility` variant divides it out; `raw` and
`f2norm` are kept for direct visual comparison against the v1 plots.

Units: the DAQ channels are acceleration (m/s^2) and the command is a
displacement (full steps). Dividing acceleration by (2*pi*f)^2 converts to
displacement, so `transmissibility` is dimensionless
(displacement out / displacement commanded).

Reads the .npy cache written by cache_ftest.py, not the 744 MB CSV directly.
"""

import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Layout mirrors the v1 campaign (ESP32S3_TMC2209_Chirp_Test/analysis/): this
# file lives in analysis/scripts/, outputs go to analysis/plots/, and only the
# primary channel sits at the plots root -- the others go to <name>_archive/
# subfolders, as v1 does with AI1_archive/.
V2_DIR = Path(__file__).resolve().parents[2]
CACHE_DIR = V2_DIR / "analysis" / "cache"
PLOT_DIR = V2_DIR / "analysis" / "plots"

FS = 20000.0

# --- Sequence constants, mirroring chirp_v2_schedule.py -------------------
MRES = 16
FULL_STEP_MM = 0.010
CRUISE_FULL_STEPS = 3.0
F_N_HZ = 176.7
NOTCH_HALFWIDTH_HZ = 56.7
NOTCH_DEPTH_RATIO = 0.10
STROKE_CLAMP_MM = 0.5
STEP_RATE_CLAMP_HZ = 200_000.0

LOG_START_HZ, LOG_END_HZ, LOG_DURATION_S = 1.0, 60.0, 120.0
LINEAR_START_HZ, LINEAR_END_HZ, LINEAR_DURATION_S = 60.0, 1000.0, 120.0
SWEEP_DURATION_S = LOG_DURATION_S + LINEAR_DURATION_S

PRE_ROLL_S = 30.0
SYNC_MARKER_DURATION_S = 0.5
MID_DWELL_S = 5.0
TAIL_S = 30.0

SEQ_MARKER_1_S = PRE_ROLL_S                                    # 30.0
SEQ_UP_START_S = SEQ_MARKER_1_S + SYNC_MARKER_DURATION_S       # 30.5
SEQ_MID_DWELL_S = SEQ_UP_START_S + SWEEP_DURATION_S            # 270.5
SEQ_MARKER_2_S = SEQ_MID_DWELL_S + MID_DWELL_S                 # 275.5
SEQ_DOWN_START_S = SEQ_MARKER_2_S + SYNC_MARKER_DURATION_S     # 276.0
SEQ_TAIL_S = SEQ_DOWN_START_S + SWEEP_DURATION_S               # 516.0
SEQ_END_S = SEQ_TAIL_S + TAIL_S                                # 546.0

MARKER_FREQ_HZ = 750.0

# Grid and plot conventions follow the v1 analysis
# (analysis/scripts/process_chirp_bode.py) so the two sets of figures read the
# same way: linear frequency axis, 1 Hz linear grid, matplotlib default series
# colours, shaded notch band.
F_LO, F_HI = 1.0, 1000.0
N_FREQ = 1000

CHANNELS = {0: "AI1", 1: "AI2", 2: "AI3"}
# AI1 carries by far the largest response and is the counterpart of v1's
# "AI0 (x, primary)"; the other two are archived alongside, as in v1.
PRIMARY_COL = 0
NOTCH_LOW_HZ = F_N_HZ - NOTCH_HALFWIDTH_HZ
NOTCH_HIGH_HZ = F_N_HZ + NOTCH_HALFWIDTH_HZ


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --- Sync-marker detection -----------------------------------------------

def marker_envelope(x, t_lo, t_hi, win_s=0.02, hop_s=0.0005):
    """Synchronous 750 Hz magnitude envelope over a time window."""
    win = int(win_s * FS)
    hop = max(1, int(hop_s * FS))
    n0, n1 = int(t_lo * FS), min(len(x), int(t_hi * FS))
    seg_all = np.asarray(x[n0:n1], dtype=np.float64)
    w = np.hanning(win)
    ph = 2 * np.pi * MARKER_FREQ_HZ * np.arange(win) / FS
    c, s = np.cos(ph) * w, np.sin(ph) * w
    starts = np.arange(0, len(seg_all) - win, hop)
    env = np.empty(len(starts))
    for i, st in enumerate(starts):
        sg = seg_all[st:st + win]
        sg = sg - sg.mean()
        env[i] = np.hypot(sg @ c, sg @ s) * 2.0 / w.sum()
    return (n0 + starts + win / 2.0) / FS, env


def find_marker_start(x, t_lo, t_hi, label):
    """Return the leading edge of the first of the three marker bursts."""
    t, env = marker_envelope(x, t_lo, t_hi)
    above = env > 0.25 * env.max()
    d = np.diff(above.astype(int))
    rises, falls = np.where(d == 1)[0] + 1, np.where(d == -1)[0] + 1
    n = min(len(rises), len(falls))
    bursts = [(t[rises[i]], t[falls[i]]) for i in range(n)]
    log(f"  {label}: peak={env.max():.4f}, {len(bursts)} bursts")
    for i, (b0, b1) in enumerate(bursts):
        log(f"      burst{i+1}: {b0:.4f} -> {b1:.4f} ({(b1-b0)*1000:.1f} ms)")
    if len(bursts) != 3:
        raise RuntimeError(f"{label}: expected 3 bursts, found {len(bursts)}")
    return bursts[0][0]


def establish_timebase(x):
    """Pin sequence t=0 in recording time and measure the DAQ/ESP32 clock ratio.

    Marker 1 is searched in the first 60 s, which is unambiguous: the up-sweep
    does not reach 750 Hz until ~238 s of sequence time, so nothing else in
    that window puts energy at the marker frequency.
    """
    log("locating sync markers ...")
    m1 = find_marker_start(x, 0.0, 60.0, "MARKER_1")
    m2_nominal = m1 + (SEQ_MARKER_2_S - SEQ_MARKER_1_S)
    m2 = find_marker_start(x, m2_nominal - 1.0, m2_nominal + 1.5, "MARKER_2")

    measured = m2 - m1
    design = SEQ_MARKER_2_S - SEQ_MARKER_1_S
    clock = measured / design
    t0 = m1 - SEQ_MARKER_1_S * clock
    log(f"  marker separation: measured={measured:.4f}s design={design:.4f}s "
        f"err={(measured-design)*1000:+.1f}ms ({(clock-1)*1e6:+.1f} ppm)")
    log(f"  sequence t=0 at recording t={t0:.4f}s, clock ratio={clock:.9f}")
    return t0, clock


# --- Frequency law, inverted analytically --------------------------------

def t_of_f_up(f):
    """Sequence time within the up-sweep at which instantaneous freq == f."""
    f = np.asarray(f, dtype=float)
    t_log = LOG_DURATION_S * np.log(f / LOG_START_HZ) / np.log(LOG_END_HZ / LOG_START_HZ)
    t_lin = LOG_DURATION_S + LINEAR_DURATION_S * (f - LINEAR_START_HZ) / (
        LINEAR_END_HZ - LINEAR_START_HZ)
    return np.where(f < LOG_END_HZ, t_log, t_lin)


def t_of_f_down(f):
    """Sequence time within the down-sweep at which instantaneous freq == f."""
    f = np.asarray(f, dtype=float)
    t_lin = LINEAR_DURATION_S * (LINEAR_END_HZ - f) / (LINEAR_END_HZ - LINEAR_START_HZ)
    t_log = LINEAR_DURATION_S + LOG_DURATION_S * np.log(LOG_END_HZ / f) / np.log(
        LOG_END_HZ / LOG_START_HZ)
    return np.where(f >= LOG_END_HZ, t_lin, t_log)


# --- Commanded amplitude -------------------------------------------------

def notch_taper(f):
    x = np.clip((np.asarray(f, float) - F_N_HZ) / NOTCH_HALFWIDTH_HZ, -1.0, 1.0)
    window = 0.5 * (1.0 + np.cos(np.pi * x))
    return 1.0 - (1.0 - NOTCH_DEPTH_RATIO) * window


def commanded_amplitude_full_steps(f):
    """cruise * notch, clamped by stroke and by peak microstep pulse rate."""
    f = np.asarray(f, float)
    tapered = CRUISE_FULL_STEPS * notch_taper(f)
    a_stroke = STROKE_CLAMP_MM / FULL_STEP_MM
    a_rate = STEP_RATE_CLAMP_HZ / (2 * np.pi * f * MRES)
    return np.minimum(np.minimum(tapered, a_stroke), a_rate)


# --- Synchronous (swept-sine) detection ----------------------------------

MIN_VALID_CYCLES = 3.0
SMEAR_TOLERANCE = 0.06


def detector_window_s(freq, cycles=20.0):
    """Window length for the synchronous detector at `freq`.

    Two competing constraints, unlike v1's flat 0.6 s cap:

    * enough cycles to resolve the tone at all, and
    * short enough that the sweep does not move the frequency appreciably
      *during* the window.

    The second is what a flat cap gets wrong at the bottom of the band: v1's
    ceiling is shorter than a single cycle below ~1.7 Hz, so the "magnitude"
    there is a fragment of a cycle -- which, once divided by (2*pi*f)^2 for
    transmissibility, inflates into a spurious low-frequency peak. The smear
    limit is derived per segment from the actual law: df/dt = f*ln(60)/120 in
    the log branch (so the limit is a constant 1.76 s), and 940/120 = 7.833
    Hz/s in the linear branch.
    """
    if freq < LOG_END_HZ:
        smear_s = SMEAR_TOLERANCE / (np.log(LOG_END_HZ / LOG_START_HZ) / LOG_DURATION_S)
    else:
        smear_s = SMEAR_TOLERANCE * freq / (
            (LINEAR_END_HZ - LINEAR_START_HZ) / LINEAR_DURATION_S)
    return max(min(cycles / freq, smear_s), 0.02)


def synchronous_magnitude(x, t_center, freq, cycles=20):
    """Hann-windowed single-frequency correlation.

    Returns NaN where the sweep simply does not dwell long enough at `freq`
    to support an estimate (fewer than MIN_VALID_CYCLES fit inside the
    smear-limited window), rather than returning a number that looks real.
    """
    win_s = detector_window_s(freq, cycles)
    if win_s * freq < MIN_VALID_CYCLES:
        return np.nan
    n = max(8, int(round(win_s * FS)))
    start = int(round(t_center * FS)) - n // 2
    stop = start + n
    if start < 0 or stop > len(x):
        return np.nan
    seg = np.asarray(x[start:stop], dtype=np.float64)
    seg = seg - seg.mean()
    w = np.hanning(n)
    tt = np.arange(n) / FS
    re = np.sum(seg * w * np.cos(2 * np.pi * freq * tt))
    im = np.sum(seg * w * np.sin(2 * np.pi * freq * tt))
    return 2.0 * np.hypot(re, im) / w.sum()


def sweep_curve(x, freq_grid, seq_times, t0, clock):
    mags = np.empty_like(freq_grid)
    for i, (f, t_seq) in enumerate(zip(freq_grid, seq_times)):
        mags[i] = synchronous_magnitude(x, t0 + t_seq * clock, f)
    return mags


def noise_floor_curve(x, freq_grid, t0, clock):
    """Idle noise floor: same detector run in the pre-roll, where the drive is
    commanding nothing. Anything within this is not a measurement."""
    centres = t0 + np.linspace(2.0, PRE_ROLL_S - 2.0, 8) * clock
    out = np.empty_like(freq_grid)
    for i, f in enumerate(freq_grid):
        vals = [synchronous_magnitude(x, tc, f) for tc in centres]
        out[i] = np.nanmedian(vals)
    return out


# --- Plotting ------------------------------------------------------------

def _decorate(ax, ylabel, title, yscale):
    """v1's axis treatment: linear frequency axis, shaded notch, framed legend."""
    ax.axvspan(NOTCH_LOW_HZ, NOTCH_HIGH_HZ, color="gray", alpha=0.15,
               label=f"resonance notch (drive tapered to {NOTCH_DEPTH_RATIO:.0%})")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel(ylabel)
    ax.set_yscale(yscale)
    ax.set_xlim(F_LO, F_HI)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()


# Only the raw magnitude is rendered on both y-scales. The f^2-normalised and
# transmissibility views span several decades and are only legible on log y,
# so their linear-y counterparts are not produced.
VARIANTS = {
    "raw": ("linear", "log"),
    "f2norm": ("log",),
    "transmissibility": ("log",),
}


def transform(variant, freq_grid, mag):
    """Return (y, ylabel) for one of the three magnitude conventions."""
    if variant == "raw":
        return mag, "acceleration magnitude (m/s$^2$)"
    if variant == "f2norm":
        return mag / freq_grid ** 2, "acceleration / f$^2$  (m/s$^2$ per Hz$^2$)"
    disp_mm = mag / (2 * np.pi * freq_grid) ** 2 * 1000.0
    cmd_mm = commanded_amplitude_full_steps(freq_grid) * FULL_STEP_MM
    return disp_mm / cmd_mm, "displacement transmissibility (out/commanded)"


def make_bode_plots(freq_grid, up_mag, down_mag, floor, out_prefix, title_prefix,
                    out_dir=None):
    for variant, yscales in VARIANTS.items():
        up_y, ylabel = transform(variant, freq_grid, up_mag)
        down_y, _ = transform(variant, freq_grid, down_mag)
        for yscale in yscales:
            suffix = "linY" if yscale == "linear" else "logY"
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(freq_grid, up_y, label="up sweep (1->1000 Hz)", lw=1.2)
            ax.plot(freq_grid, down_y, label="down sweep (1000->1 Hz)", lw=1.2,
                    alpha=0.8)
            _decorate(ax, ylabel, f"{title_prefix} -- {variant}, {yscale} y-scale", yscale)
            if variant == "transmissibility":
                # Measured on this record: the response dips only ~1.3x across
                # the notch where the command dips 10x, so dividing by the
                # commanded amplitude over-corrects and manufactures a peak.
                # See the module docstring / README note on the presliding regime.
                ax.annotate(
                    "normalisation unreliable in the notch band:\n"
                    "response does not scale with commanded amplitude here",
                    xy=(F_N_HZ, 0.97), xycoords=("data", "axes fraction"),
                    xytext=(0, -28), textcoords="offset points",
                    ha="center", va="top", fontsize=7.5, color="#B00020")
            fig.tight_layout()
            target = out_dir if out_dir is not None else PLOT_DIR
            target.mkdir(parents=True, exist_ok=True)
            out = target / f"{out_prefix}_{variant}_{suffix}.png"
            fig.savefig(out, dpi=150)
            plt.close(fig)
            log(f"  saved {out.relative_to(PLOT_DIR)}")


def make_overview(freq_grid, up_by_ch, floor_by_ch, out_path):
    """Single-panel all-axis overview: up sweep, f^2-normalised, log y."""
    fig, ax = plt.subplots(figsize=(10, 6))
    for col, mag in up_by_ch.items():
        y, ylabel = transform("f2norm", freq_grid, mag)
        ax.plot(freq_grid, y, label=CHANNELS[col], lw=1.1)
    _decorate(ax, ylabel, "All-axis overview (up sweep, f^2-normalized, log y-scale)",
              "log")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log(f"  saved {out_path.name}")


def make_general_visualization(accel, t0, clock, out_path):
    """Time-domain overview: what the recording actually contains, per axis,
    with the commanded sequence segments marked."""
    segs = [
        ("pre-roll", 0.0, SEQ_MARKER_1_S),
        ("marker 1", SEQ_MARKER_1_S, SEQ_UP_START_S),
        ("up-sweep", SEQ_UP_START_S, SEQ_MID_DWELL_S),
        ("mid-dwell", SEQ_MID_DWELL_S, SEQ_MARKER_2_S),
        ("marker 2", SEQ_MARKER_2_S, SEQ_DOWN_START_S),
        ("down-sweep", SEQ_DOWN_START_S, SEQ_TAIL_S),
        ("tail", SEQ_TAIL_S, SEQ_END_S),
    ]
    fig, axes = plt.subplots(4, 1, figsize=(13, 7.5), sharex=True,
                             gridspec_kw={"height_ratios": [0.32, 1, 1, 1]})

    # Dedicated segment ruler, so labels never collide with the traces.
    ruler = axes[0]
    for k, (name, s0, s1) in enumerate(segs):
        r0, r1 = t0 + s0 * clock, t0 + s1 * clock
        ruler.axvspan(r0, r1, color=("#DDDDDD" if k % 2 == 0 else "#EFEFEF"), lw=0)
        if (s1 - s0) >= 20:
            ruler.text((r0 + r1) / 2, 0.5, name, ha="center", va="center", fontsize=8.5)
    # Short segments (0.5 s markers, 5 s dwell) get leader labels above the
    # ruler. marker 2 and mid-dwell are 0.5 s apart on a 543 s axis, so the
    # labels are staggered in height and nudged apart horizontally -- at this
    # scale they would otherwise print on top of each other.
    short = [(n, s0, s1) for n, s0, s1 in segs if (s1 - s0) < 20]
    span = accel.shape[0] / FS
    for k, (name, s0, s1) in enumerate(short):
        rc = t0 + (s0 + s1) / 2 * clock
        y = 1.55 + 0.75 * (k % 2)
        dx = (-0.035 if k % 2 == 0 else 0.035) * span
        ruler.annotate(name, xy=(rc, 1.0), xytext=(rc + dx, y),
                       textcoords="data", ha="center", fontsize=7.5,
                       arrowprops=dict(arrowstyle="-", lw=0.7, color="#777777",
                                       shrinkA=0, shrinkB=1))
    ruler.set_ylim(0, 3.2)
    ruler.set_yticks([])
    for side in ("top", "right", "left"):
        ruler.spines[side].set_visible(False)

    # Envelope per axis (decimated peak-hold, so bursts stay visible).
    dec = 200
    n = (accel.shape[0] // dec) * dec
    t_dec = (np.arange(n // dec) * dec + dec / 2) / FS
    for col, ax in zip((0, 1, 2), axes[1:4]):
        block = np.asarray(accel[:n, col], dtype=np.float32).reshape(-1, dec)
        # Same default-cycle colour each axis gets in the all-axis overview, so
        # AI1/AI2/AI3 read consistently across the two figures.
        ax.fill_between(t_dec, block.min(axis=1), block.max(axis=1),
                        lw=0, color=f"C{col}", label=CHANNELS[col])
        ax.set_ylabel(f"{CHANNELS[col]}\n(m/s$^2$)")
        ax.grid(True, alpha=0.25, lw=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(frameon=False, loc="upper right", fontsize=8.5)

    axes[3].set_xlabel("recording time (s)")

    for ax in axes[1:]:
        for _, s0, _ in segs:
            ax.axvline(t0 + s0 * clock, color="#555555", lw=0.6, alpha=0.45)
    ruler.set_title("v2 chirp recording (Ftest.csv) -- full record, all axes, "
                    "with commanded sequence overlaid")
    axes[3].set_xlim(0, accel.shape[0] / FS)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log(f"  saved {out_path.name}")


def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    accel = np.load(CACHE_DIR / "accel.npy", mmap_mode="r")
    log(f"cache: {accel.shape[0]:,} samples x {accel.shape[1]} ch, "
        f"{accel.shape[0]/FS:.3f}s at {FS:g} Hz")

    primary = np.asarray(accel[:, 0], dtype=np.float64)
    t0, clock = establish_timebase(primary)

    rec_end = accel.shape[0] / FS
    down_end = t0 + SEQ_TAIL_S * clock
    if down_end > rec_end:
        raise RuntimeError(f"recording ends at {rec_end:.1f}s, before the down-sweep "
                           f"finishes at {down_end:.1f}s -- sweep data is missing")
    log(f"  down-sweep ends at {down_end:.2f}s, recording ends {rec_end:.2f}s: "
        f"all sweep data present ({t0 + SEQ_END_S*clock - rec_end:.1f}s of idle tail truncated)")

    freq_grid = np.linspace(F_LO, F_HI, N_FREQ)
    up_times = SEQ_UP_START_S + t_of_f_up(freq_grid)
    down_times = SEQ_DOWN_START_S + t_of_f_down(freq_grid)

    log("running synchronous detection ...")
    results, floors = {}, {}
    for col, name in CHANNELS.items():
        ch = np.asarray(accel[:, col], dtype=np.float64)
        ts = time.time()
        up = sweep_curve(ch, freq_grid, up_times, t0, clock)
        dn = sweep_curve(ch, freq_grid, down_times, t0, clock)
        fl = noise_floor_curve(ch, freq_grid, t0, clock)
        results[col] = (up, dn)
        floors[col] = fl
        snr = np.nanmedian(up / fl)
        log(f"  {name}: {time.time()-ts:.1f}s, median up-sweep SNR vs idle floor = {snr:.1f}x")

    for col, name in CHANNELS.items():
        up, dn = results[col]
        primary = col == PRIMARY_COL
        make_bode_plots(
            freq_grid, up, dn, floors[col], f"bode_{name}",
            f"Bode: {name}{' (primary)' if primary else ''}",
            out_dir=None if primary else PLOT_DIR / f"{name}_archive")

    make_overview(freq_grid, {c: results[c][0] for c in CHANNELS}, floors,
                  PLOT_DIR / "bode_all_axes_overview.png")
    make_general_visualization(accel, t0, clock, PLOT_DIR / "recording_overview.png")
    valid = np.isfinite(results[0][0])
    log(f"detector validity: {valid.sum()}/{len(freq_grid)} grid points; "
        f"lowest resolvable frequency = {freq_grid[valid][0]:.2f} Hz "
        f"(below this the sweep does not dwell for {MIN_VALID_CYCLES:g} cycles)")

    np.savez(
        PLOT_DIR / "bode_data.npz",
        freq_grid=freq_grid,
        commanded_full_steps=commanded_amplitude_full_steps(freq_grid),
        up_mag_AI1=results[0][0], down_mag_AI1=results[0][1], floor_AI1=floors[0],
        up_mag_AI2=results[1][0], down_mag_AI2=results[1][1], floor_AI2=floors[1],
        up_mag_AI3=results[2][0], down_mag_AI3=results[2][1], floor_AI3=floors[2],
        t0_recording_s=t0, clock_ratio=clock, fs=FS,
    )
    log("saved bode_data.npz")
    log("done.")


if __name__ == "__main__":
    main()
