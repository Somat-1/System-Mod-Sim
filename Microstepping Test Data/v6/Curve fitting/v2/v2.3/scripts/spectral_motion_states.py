"""v2.3: is the ~360 Hz peak a real mode? Spectra of moving segments.

v2 found a peak at 362.8 +/- 3.0 Hz in the post-recording idle tail, and
[v2.2](../../v2.2/README.md) reproduced it at 357.9 Hz on a different segment
with a different detrend. Both measured it in the **same state**: nothing
commanded, axis at rest. A peak in the noise floor of one state cannot tell
you whether it is a mechanical mode of the stage, an artifact of the driver
holding position, or an alias of something else entirely.

This version tests it against the states the tail cannot reach, using the same
Welch machinery on short windows cut out of the campaign itself:

  1. **Idle tail** -- the control, carried over from v2.2.
  2. **Oscillation dwell** -- the settled second between two commanded
     microsteps. Axis at rest, driver energised and holding, but in the middle
     of the campaign rather than after it.
  3. **Mid-ramp** -- a constant-velocity stretch of a trajectory leg. Axis
     moving.
  4. **Ring-down** -- the free decay right after a commanded microstep. This
     is the decisive one: a free decay rings at the system's own damped
     natural frequency, so a peak that shows up *here* is a mechanical
     property, not a noise-floor curiosity.

The three quasi-stationary states (1-3) share `psd_by_motion_state.png` and
are computed identically -- same `NPERSEG`, same Hann window, same overlap --
so their levels are directly comparable. The ring-down is deliberately NOT on
that figure: it is a short decaying burst, not a stationary process, so its
spectrum is an amplitude spectrum of a transient and does not share units with
a PSD. It gets `ringdown_fit.png` to itself, where a damped-sinusoid fit
reports f_d, zeta and Q per MRES.

### Reading the ramp spectra

`psd_ramp_grid.png` exists because the moving spectra are dominated by
commanded stepping, not by the mode. Two families of line are marked on every
panel:

  * **Harmonics of the FULL-STEP rate** (`rate_fsps`, i.e. 27.5 / 70 / 200 Hz)
    -- not the microstep rate. The stage advances one full step at a time
    regardless of MRES; microstepping subdivides how it gets there. These
    harmonics are what actually dominate: at MRES 1 and 4 they run six orders
    of magnitude above the broadband floor.
  * **The aliased microstep pulse rate** `rate_fsps x mres`, folded into
    [0, f_s/2]. At MRES 4 fast that is 800 Hz -> 200 Hz; at MRES 32 fast,
    6400 Hz -> 400 Hz. Above Nyquist these do not disappear, they fold.

The 27.5 Hz fundamental puts a harmonic within ~14 Hz of 350 Hz no matter
what, and 70 x 5 = 350 exactly, so the slow and moderate legs cannot separate
the mode from stepping content. **The fast legs are the clean test**: their
harmonics land at 200 and 400 Hz, leaving 330-370 Hz free of any commanded
line.

Run from this directory:

    python spectral_motion_states.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy import signal, optimize
import matplotlib
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
V23_DIR = HERE.parent
V2_DIR = V23_DIR.parent
CF_DIR = V2_DIR.parent
V6_DIR = CF_DIR.parent

sys.path.insert(0, str(V6_DIR / "analysis" / "scripts"))
sys.path.insert(0, str(CF_DIR / "v1" / "scripts"))
sys.path.insert(0, str(V2_DIR / "scripts"))

from ids_common import load_ids, still_mask, ramp_extent, FULL_STEP_UM
from plot_v6_exp_run import CSV, find_trajectories, MRES_VALUES, RATES_FSPS
from fit_v6_curves import find_oscillation_seed_times, fit_steps
from analyze_post_recording_tail import find_post_recording_tail

matplotlib.use("Agg", force=True)
OUT_DIR = V23_DIR / "rendered_assets"

# One Welch configuration for every quasi-stationary state, so the three are
# directly comparable. 512 is set by the shortest window available -- a single
# oscillation dwell yields 650 usable samples -- and every other state is cut
# with the same nperseg rather than a longer one it could afford.
NPERSEG = 512
WINDOW = "hann"
OVERLAP = 0.5
DETREND_ORDER = 3          # matches v2.2's selected order

DWELL_START_S, DWELL_END_S = 0.30, 0.95   # settled part of a 1 s dwell
RING_S = 0.160                            # ring-down capture after an edge
RING_FIT_FROM_S, RING_FIT_TO_S = 0.005, 0.150
MODE_BAND = (330.0, 370.0)
FLOOR_BAND = (200.0, 300.0)
RAMP_CORE = (0.20, 0.80)   # fraction of the out-ramp used, clear of both ends

LABELS = ["D-slow", "D-mod", "D-fast", "I-slow", "I-mod", "I-fast"]
RATE_BY_LABEL = {f"{m}-{n}": r for m in "DI"
                 for n, r in zip(("slow", "mod", "fast"), RATES_FSPS)}

# Validated trio (dataviz reference slots 1-3, aqua darkened to clear 3:1 on
# white); all pairs pass normal/protan/deutan separation. One state, one hue.
C_TAIL = "#2563a6"
C_DWELL = "#d55e00"
C_RAMP = "#008b6d"
C_INK, C_INK2, C_MUTED = "#17202a", "#46515c", "#77818a"
C_GRID, C_BASELINE, C_BAND = "#dde3e8", "#aeb7bf", "#fff1e8"
LABEL_BBOX = dict(facecolor="white", edgecolor=C_GRID, linewidth=0.6,
                  alpha=0.94, pad=3.0)


def apply_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": C_BASELINE, "axes.labelcolor": C_INK2,
        "axes.titlesize": 11, "axes.titleweight": "semibold",
        "axes.titlepad": 9, "axes.labelsize": 10, "axes.grid": False,
        "grid.color": C_GRID, "grid.linewidth": 0.7, "grid.alpha": 0.8,
        "xtick.color": C_MUTED, "ytick.color": C_MUTED,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "legend.frameon": False, "legend.edgecolor": C_GRID,
        "legend.fontsize": 9, "font.size": 10, "text.color": C_INK,
        "figure.dpi": 150, "savefig.dpi": 200,
        "savefig.facecolor": "white", "savefig.bbox": "tight",
    })


def style(ax, xlabel=None, ylabel=None):
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(C_BASELINE)
    ax.grid(axis="y")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


# ------------------------------------------------------------- primitives ---

def detrend_poly(y_nm, dt_s, order=DETREND_ORDER):
    x = np.arange(y_nm.size) * dt_s
    return y_nm - Polynomial.fit(x, y_nm, order)(x)


def psd_of(segments, fs):
    """Welch PSD averaged over a list of equal-treatment segments (nm -> nm^2/Hz).

    Averaging across many short windows rather than concatenating them: the
    dwells are separated by commanded steps, so joining them end to end would
    manufacture step discontinuities that smear across the whole spectrum.
    """
    acc, f = None, None
    for seg in segments:
        if seg.size < NPERSEG:
            continue
        f, p = signal.welch(seg, fs=fs, window=WINDOW, nperseg=NPERSEG,
                            noverlap=int(NPERSEG * OVERLAP), detrend="constant")
        acc = p if acc is None else acc + p
    if acc is None:
        return None, None, 0
    n = sum(1 for s in segments if s.size >= NPERSEG)
    return f, acc / n, n


def alias_into_band(f_hz, fs):
    """Where a frequency above Nyquist actually appears after folding."""
    r = f_hz % fs
    return min(r, fs - r)


def band_stat(f, p, band, floor=FLOOR_BAND):
    inb = (f >= band[0]) & (f <= band[1])
    inf = (f >= floor[0]) & (f <= floor[1])
    if not inb.any() or not inf.any():
        return float("nan"), float("nan"), float("nan")
    pk = float(np.max(p[inb]))
    fl = float(np.median(p[inf]))
    return pk, fl, (pk / fl if fl > 0 else float("nan"))


# ------------------------------------------------------ segment extraction ---

def oscillation_segments(t, pos, dt_s, leg1_idx):
    """(dwell segments [nm], ring-down stack [nm], edge count) for one block."""
    seeds, method = find_oscillation_seed_times(t, pos, dt_s, leg1_idx)
    lo = max(0, int(np.searchsorted(t, seeds[0] - 0.15)))
    hi = min(pos.size, int(np.searchsorted(t, seeds[-1] + 0.85)))
    local_t = t[lo:hi] - t[lo]
    y_um = pos[lo:hi] - np.median(pos[lo:hi])
    _, edges = fit_steps(local_t, y_um, seeds - t[lo], dt_s)
    y_nm = y_um * 1000.0

    a0, a1 = int(round(DWELL_START_S / dt_s)), int(round(DWELL_END_S / dt_s))
    dwells = [y_nm[int(e) + a0:int(e) + a1] - np.mean(y_nm[int(e) + a0:int(e) + a1])
              for e in edges if int(e) + a1 <= y_nm.size]

    n_ring = int(round(RING_S / dt_s))
    # Forward transitions only (even index): the two directions are not the
    # same excitation -- lost motion makes the return differ -- and averaging
    # them together would blur the decay this figure is trying to measure.
    rings = [y_nm[int(e):int(e) + n_ring]
             for k, e in enumerate(edges) if k % 2 == 0 and int(e) + n_ring <= y_nm.size]
    return dwells, np.asarray(rings), method


def midramp_segment(t, pos, dt_s, still, span):
    """Constant-velocity core of a leg's out-ramp, in nm, cubic-detrended."""
    s0, _apex, e0 = span
    a, b = ramp_extent(still, s0, e0, dt_s, pad_s=1.0)
    seg_still = still[a:b]
    moving = ~seg_still
    d = np.diff(moving.astype(np.int8))
    st = list(np.flatnonzero(d == 1) + 1)
    en = list(np.flatnonzero(d == -1) + 1)
    if moving[0]:
        st.insert(0, 0)
    if moving[-1]:
        en.append(moving.size - 1)
    runs = [(x, y) for x, y in zip(st, en) if (y - x) * dt_s > 0.5]
    if len(runs) != 2:
        return None
    o0, o1 = runs[0]
    n = o1 - o0
    y = pos[a + o0 + int(RAMP_CORE[0] * n):a + o0 + int(RAMP_CORE[1] * n)] * 1000.0
    if y.size < NPERSEG:
        return None
    return detrend_poly(y, dt_s)


# ------------------------------------------------------------ ring-down fit ---

def damped_cosine(x, amp, zeta, wn, phase, offset):
    wd = wn * np.sqrt(max(1.0 - zeta * zeta, 1e-9))
    return amp * np.exp(-zeta * wn * x) * np.cos(wd * x + phase) + offset


def ringdown_spectrum(ring_median, dt_s):
    """Amplitude spectrum (nm) of the decaying burst."""
    y = ring_median - ring_median.mean()
    w = np.hanning(y.size)
    return (np.fft.rfftfreq(y.size, dt_s),
            np.abs(np.fft.rfft(y * w)) * 2 / w.sum())


def spectral_peak(fr, amp, band):
    m = (fr >= band[0]) & (fr <= band[1])
    if not m.any():
        return float("nan"), float("nan")
    i = int(np.argmax(amp[m]))
    return float(fr[m][i]), float(amp[m][i])


def fit_ringdown(ring_median, dt_s, f_guess=None):
    x = np.arange(ring_median.size) * dt_s
    m = (x >= RING_FIT_FROM_S) & (x <= RING_FIT_TO_S)
    xf, yf = x[m] - RING_FIT_FROM_S, ring_median[m]
    # Seed from the MEASURED spectrum, not a fixed number. At MRES 1 the
    # commanded step is a full 10 um -- 32x the MRES 32 microstep -- and the
    # response carries a second, lower, more heavily damped component beside
    # the mode the finer steps excite. Started at a fixed 350 Hz the optimiser
    # converged onto that ~173 Hz component instead, so MRES 1 reported a
    # different mode from the other three for the same measurement. Seeding at
    # the dominant peak inside MODE_BAND makes all four fits measure the same
    # feature; the lower component is reported separately rather than being
    # picked silently by an optimiser.
    if f_guess is None:
        _fr, _amp = ringdown_spectrum(ring_median, dt_s)
        f_guess, _ = spectral_peak(_fr, _amp, MODE_BAND)
        if not np.isfinite(f_guess):
            f_guess = 350.0
    p0 = [np.ptp(yf) / 2, 0.05, 2 * np.pi * f_guess, 0.0, float(np.mean(yf[-30:]))]
    try:
        popt, _ = optimize.curve_fit(damped_cosine, xf, yf, p0=p0, maxfev=40000)
    except Exception:
        return None
    amp, zeta, wn, phase, offset = popt
    zeta = abs(zeta)
    fn = abs(wn) / (2 * np.pi)
    resid = yf - damped_cosine(xf, *popt)
    r2 = float(1 - np.var(resid) / np.var(yf))
    _fr2, _sp2 = ringdown_spectrum(ring_median, dt_s)
    f_band, a_band = spectral_peak(_fr2, _sp2, MODE_BAND)
    f_low, a_low = spectral_peak(_fr2, _sp2, (80.0, MODE_BAND[0]))
    return dict(amp_nm=float(abs(amp)), zeta=float(zeta), fn_hz=float(fn),
                fd_hz=float(fn * np.sqrt(max(1 - zeta ** 2, 0.0))),
                q=float(1 / (2 * zeta)) if zeta > 0 else float("inf"),
                tau_ms=float(1e3 / (zeta * abs(wn))) if zeta > 0 else float("inf"),
                r2=r2,
                # A low R2 means the single-mode model did not describe this
                # decay; zeta and Q are then curve-fit outputs, not measured
                # damping, and are labelled as such rather than quoted.
                damping_reliable=bool(r2 >= 0.60),
                spectral_peak_in_band_hz=f_band,
                spectral_peak_in_band_nm=a_band,
                secondary_peak_80_330hz_hz=f_low,
                secondary_peak_80_330hz_nm=a_low,
                popt=[float(v) for v in popt], x=xf, y=yf)


# ------------------------------------------------------------------ plots ----

def plot_states(states, fs, out_dir):
    """Compare states on one shared scale so level differences stay honest."""
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8), sharex=True,
                             sharey=True)
    for panel, (ax, mres) in enumerate(zip(axes.ravel(), MRES_VALUES)):
        st = states[mres]
        ax.axvspan(*MODE_BAND, color=C_BAND, zorder=0, linewidth=0,
                   label="Mode band" if panel == 0 else None)
        for key, colour, label in (
                ("tail", C_TAIL, "Idle tail"),
                ("dwell", C_DWELL, "At rest, holding"),
                ("ramp", C_RAMP, "Moving, fast leg")):
            f, p, _n = st[key]
            if f is None:
                continue
            ax.loglog(f[1:], p[1:], color=colour, linewidth=1.05, zorder=3,
                      label=label if panel == 0 else None)
        ax.set_title(f"MRES {mres}", loc="left")
        ax.set_xlim(fs / NPERSEG, fs / 2)
        ax.grid(True, which="major", axis="both")
        ax.grid(True, which="minor", axis="both", linewidth=0.4, alpha=0.35)
        style(ax)
        ax.text(0.98, 0.95,
                f"mode / floor\nhold  {st['dwell_band_ratio']:.0f}×   ·   move  {st['ramp_band_ratio']:.0f}×",
                transform=ax.transAxes, fontsize=8.5, color=C_INK2,
                ha="right", va="top", linespacing=1.35, bbox=LABEL_BBOX)
        ax.label_outer()

    handles, labels_ = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper left", bbox_to_anchor=(0.065, 0.905),
               ncol=4, frameon=False, columnspacing=1.7, handlelength=2.3)
    fig.suptitle("The same mode band across operating states", x=0.065,
                 ha="left", y=0.98, fontsize=16, fontweight="semibold")
    fig.text(0.065, 0.93,
             "Shared axes and identical Welch treatment · moving spectra are 3–4 orders higher than at-rest spectra",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.supxlabel("Frequency [Hz]", x=0.53, y=0.035, fontsize=10,
                  color=C_INK2)
    fig.supylabel("PSD [nm²/Hz]", x=0.018, fontsize=10, color=C_INK2)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.84, bottom=0.09,
                        hspace=0.18, wspace=0.10)
    fig.savefig(out_dir / "psd_by_motion_state.png")
    plt.close(fig)


def plot_ramp_grid(ramps, fs, out_dir):
    """Small multiples with one legend and row-consistent PSD scales."""
    rates = ("slow", "mod", "fast")
    headings = ("Slow", "Moderate", "Fast")
    fig, axes = plt.subplots(4, 3, figsize=(14.2, 12.4), sharex=True,
                             sharey="row")
    for row, mres in enumerate(MRES_VALUES):
        for col, (rate_name, heading) in enumerate(zip(rates, headings)):
            ax = axes[row, col]
            entry = ramps[mres].get(rate_name)
            if entry is None:
                ax.axis("off")
                continue
            f, p, fsps = entry["f"], entry["p"], entry["fsps"]
            first_panel = row == 0 and col == 0
            ax.axvspan(*MODE_BAND, color=C_BAND, zorder=0, linewidth=0,
                       label="330–370 Hz mode band" if first_panel else None)
            ax.loglog(f[1:], p[1:], color=C_RAMP, linewidth=0.95, zorder=3,
                      label="Moving-stage PSD" if first_panel else None)
            for k in range(1, 15):
                fh = k * fsps
                if fh >= fs / 2:
                    break
                ax.axvline(fh, color=C_MUTED, linewidth=0.65, linestyle=":",
                           alpha=0.65, zorder=2,
                           label="Full-step harmonics"
                                 if first_panel and k == 1 else None)
            pulse = fsps * mres
            pulse_visible = alias_into_band(pulse, fs)
            ax.axvline(pulse_visible, color=C_DWELL, linewidth=1.25,
                       linestyle=(0, (4, 3)), zorder=4,
                       label="Microstep rate (folded)" if first_panel else None)
            ax.set_xlim(fs / NPERSEG, fs / 2)
            ax.grid(True, which="major", axis="both")
            ax.grid(True, which="minor", axis="both", linewidth=0.4,
                    alpha=0.32)
            style(ax)
            if row == 0:
                ax.set_title(f"{heading}\n{fsps:g} full-step/s", loc="center")
            if col == 0:
                ax.set_ylabel(f"MRES {mres}", fontsize=10,
                              fontweight="semibold", labelpad=17)
            folded = " folded" if pulse > fs / 2 else ""
            ax.text(0.02, 0.05, f"pulse → {pulse_visible:.0f} Hz{folded}",
                    transform=ax.transAxes, fontsize=7.8, color=C_DWELL,
                    ha="left", va="bottom", bbox=LABEL_BBOX)
            ax.text(0.98, 0.95, f"mode / floor  {entry['band_ratio']:.0f}×",
                    transform=ax.transAxes, fontsize=8.0, color=C_INK2,
                    ha="right", va="top", bbox=LABEL_BBOX)
            ax.label_outer()

    handles, labels_ = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper left", bbox_to_anchor=(0.06, 0.91),
               ncol=4, frameon=False, columnspacing=1.5, handlelength=2.2)
    fig.suptitle("Commanded stepping dominates the ramp spectra", x=0.06,
                 ha="left", y=0.985, fontsize=16, fontweight="semibold")
    fig.text(0.06, 0.94,
             "Columns compare speed; rows compare microstep resolution · fast legs leave the mode band clear of step harmonics",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.supxlabel("Frequency [Hz]", x=0.53, y=0.025, fontsize=10,
                  color=C_INK2)
    fig.supylabel("PSD [nm²/Hz]", x=0.012, fontsize=10, color=C_INK2)
    fig.subplots_adjust(left=0.08, right=0.985, top=0.84, bottom=0.07,
                        hspace=0.18, wspace=0.10)
    fig.savefig(out_dir / "psd_ramp_grid.png")
    plt.close(fig)


def plot_ringdown(rings, fits, dt_s, out_dir):
    """Give every MRES a readable row: decay on the left, spectrum on the right."""
    fig, axes = plt.subplots(4, 2, figsize=(14.2, 12.2),
                             gridspec_kw={"width_ratios": [1.75, 1.0]},
                             sharex="col")
    for row, mres in enumerate(MRES_VALUES):
        stack = rings[mres]
        fit = fits[mres]
        x_ms = np.arange(stack.shape[1]) * dt_s * 1e3
        med = np.median(stack, axis=0)

        ax = axes[row, 0]
        ax.axhline(0, color=C_BASELINE, linewidth=0.9, zorder=1)
        for cycle_i, cycle in enumerate(stack):
            ax.plot(x_ms, cycle, color=C_RAMP, alpha=0.16, linewidth=0.65,
                    zorder=2, label="Individual cycles" if row == 0 and cycle_i == 0 else None)
        ax.plot(x_ms, med, color=C_TAIL, linewidth=1.8, zorder=4,
                label="Median response" if row == 0 else None)
        if fit:
            ax.plot((fit["x"] + RING_FIT_FROM_S) * 1e3,
                    damped_cosine(fit["x"], *fit["popt"]), color=C_DWELL,
                    linewidth=1.65, linestyle=(0, (5, 2)), zorder=5,
                    label="Single-mode fit" if row == 0 else None)
        ax.set_title(f"MRES {mres}  ·  {FULL_STEP_UM / mres:g} µm command",
                     loc="left", fontsize=10.5)
        ax.set_xlim(0, x_ms[-1])
        style(ax)
        if fit:
            if fit["damping_reliable"]:
                result = (f"peak {fit['spectral_peak_in_band_hz']:.1f} Hz   ·   "
                          f"τ {fit['tau_ms']:.1f} ms   ·   Q {fit['q']:.0f}   ·   "
                          f"R² {fit['r2']:.2f}")
            else:
                result = (f"peak {fit['spectral_peak_in_band_hz']:.1f} Hz   ·   "
                          f"R² {fit['r2']:.2f}   ·   single-mode damping unreliable")
            ax.text(0.985, 0.08, result, transform=ax.transAxes, fontsize=8.3,
                    color=C_INK2, ha="right", va="bottom", bbox=LABEL_BBOX)

        ax = axes[row, 1]
        fr, amp = ringdown_spectrum(med, dt_s)
        ax.axvspan(*MODE_BAND, color=C_BAND, zorder=0, linewidth=0,
                   label="330–370 Hz mode band" if row == 0 else None)
        ax.semilogy(fr[1:], amp[1:], color=C_RAMP, linewidth=1.2, zorder=3,
                    label="Transient spectrum" if row == 0 else None)
        if fit:
            peak = fit["spectral_peak_in_band_hz"]
            ax.axvline(peak, color=C_DWELL, linewidth=1.2,
                       linestyle=(0, (4, 3)), zorder=4,
                       label="Measured band peak" if row == 0 else None)
            ax.text(peak, 0.92, f"{peak:.0f} Hz",
                    transform=ax.get_xaxis_transform(), fontsize=8.2,
                    color=C_DWELL, ha="center", va="top", bbox=LABEL_BBOX)
        ax.set_xlim(0, 1 / dt_s / 2)
        style(ax)
        ax.grid(True, which="major", axis="both")
        ax.grid(True, which="minor", axis="y", linewidth=0.4, alpha=0.32)
        if row < len(MRES_VALUES) - 1:
            axes[row, 0].tick_params(labelbottom=False)
            axes[row, 1].tick_params(labelbottom=False)

    left_h, left_l = axes[0, 0].get_legend_handles_labels()
    right_h, right_l = axes[0, 1].get_legend_handles_labels()
    fig.legend(left_h + right_h, left_l + right_l, loc="upper left",
               bbox_to_anchor=(0.06, 0.90), ncol=6, frameon=False,
               columnspacing=1.35, handlelength=2.1, fontsize=8.7)
    axes[-1, 0].set_xlabel("Time after command [ms]")
    axes[-1, 1].set_xlabel("Frequency [Hz]")
    fig.text(0.018, 0.47, "Position [nm]", rotation=90, va="center",
             fontsize=10, color=C_INK2)
    fig.text(0.655, 0.47, "Amplitude [nm]", rotation=90, va="center",
             fontsize=10, color=C_INK2)
    fig.suptitle("Free decay after one commanded microstep", x=0.06,
                 ha="left", y=0.985, fontsize=16, fontweight="semibold")
    fig.text(0.06, 0.94,
             "Each row pairs the time-domain ring-down with its transient spectrum · damping is shown only when the single-mode fit is credible",
             ha="left", fontsize=9.5, color=C_MUTED)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.83, bottom=0.07,
                        hspace=0.34, wspace=0.22)
    fig.savefig(out_dir / "ringdown_fit.png")
    plt.close(fig)

def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t, pos, meta = load_ids(CSV)
    dt_s, fs = meta["dt_s"], meta["fs_hz"]
    print(f"loaded {CSV.name}: {meta['n']:,} samples @ {fs:.0f} Hz")

    spans = find_trajectories(t, pos, dt_s)
    still = still_mask(pos, dt_s)

    # Control: the same 60 s window v2.2 analyses.
    tstart, tend, _ = find_post_recording_tail(t, pos, dt_s)
    tl = tstart + int(round(0.050 / dt_s))
    tail_nm = detrend_poly((pos[tl:tl + 60000] -
                            np.median(pos[tl:tl + 60000])) * 1000.0, dt_s)
    f_tail, p_tail, n_tail = psd_of([tail_nm], fs)
    tail_pk, tail_fl, tail_ratio = band_stat(f_tail, p_tail, MODE_BAND)
    tail_peak_hz = float(f_tail[(f_tail >= MODE_BAND[0]) & (f_tail <= MODE_BAND[1])][
        np.argmax(p_tail[(f_tail >= MODE_BAND[0]) & (f_tail <= MODE_BAND[1])])])
    print(f"idle tail: {MODE_BAND[0]:.0f}-{MODE_BAND[1]:.0f} Hz peak at "
          f"{tail_peak_hz:.1f} Hz, {tail_ratio:.0f}x the 200-300 Hz floor")

    states, ramps, rings, fits, summary_mres = {}, {}, {}, {}, {}
    for gi, mres in enumerate(MRES_VALUES):
        group = spans[gi * 6:(gi + 1) * 6]
        dwells, ring_stack, seed_method = oscillation_segments(
            t, pos, dt_s, group[0][0])
        f_dw, p_dw, n_dw = psd_of(dwells, fs)
        dw_pk, dw_fl, dw_ratio = band_stat(f_dw, p_dw, MODE_BAND)
        dwell_peak_hz = float(f_dw[(f_dw >= MODE_BAND[0]) & (f_dw <= MODE_BAND[1])][
            np.argmax(p_dw[(f_dw >= MODE_BAND[0]) & (f_dw <= MODE_BAND[1])])])

        ramps[mres] = {}
        for (span, label) in zip(group, LABELS):
            if not label.startswith("D-"):
                continue
            rate_name = label.split("-")[1]
            seg = midramp_segment(t, pos, dt_s, still, span)
            if seg is None:
                continue
            f_r, p_r, _ = psd_of([seg], fs)
            _, _, r_ratio = band_stat(f_r, p_r, MODE_BAND)
            ramps[mres][rate_name] = dict(f=f_r, p=p_r, band_ratio=r_ratio,
                                          fsps=RATE_BY_LABEL[label],
                                          n_samples=int(seg.size))

        fast = ramps[mres].get("fast")
        states[mres] = {
            "tail": (f_tail, p_tail, n_tail),
            "dwell": (f_dw, p_dw, n_dw),
            "ramp": (fast["f"], fast["p"], 1) if fast else (None, None, 0),
            "dwell_band_ratio": dw_ratio,
            "ramp_band_ratio": fast["band_ratio"] if fast else float("nan"),
        }

        rings[mres] = ring_stack
        fit = fit_ringdown(np.median(ring_stack, axis=0), dt_s)
        fits[mres] = fit
        print(f"  MRES {mres:2d}: dwell peak {dwell_peak_hz:6.1f} Hz "
              f"({dw_ratio:6.0f}x floor) | moving {states[mres]['ramp_band_ratio']:6.0f}x"
              f" | ring-down peak {fit['spectral_peak_in_band_hz']:6.1f} Hz, "
              f"fit f_d {fit['fd_hz']:6.1f} Hz, zeta {fit['zeta']:.4f}"
              f"{'' if fit['damping_reliable'] else ' (UNRELIABLE)'}, "
              f"R2 {fit['r2']:.2f}, 2nd {fit['secondary_peak_80_330hz_hz']:.0f} Hz")

        summary_mres[f"MRES_{mres}"] = {
            "oscillation_seed_method": seed_method,
            "dwell": {"n_windows": n_dw, "peak_hz": dwell_peak_hz,
                      "peak_psd_nm2_per_hz": dw_pk,
                      "band_vs_floor_ratio": dw_ratio},
            "mid_ramp_band_vs_floor_ratio": {
                rn: e["band_ratio"] for rn, e in ramps[mres].items()},
            "ringdown": {k: fit[k] for k in
                         ("spectral_peak_in_band_hz", "spectral_peak_in_band_nm",
                          "secondary_peak_80_330hz_hz", "secondary_peak_80_330hz_nm",
                          "fd_hz", "fn_hz", "zeta", "q", "tau_ms", "r2",
                          "damping_reliable", "amp_nm")}
            if fit else None,
        }

    plot_states(states, fs, OUT_DIR)
    plot_ramp_grid(ramps, fs, OUT_DIR)
    plot_ringdown(rings, fits, dt_s, OUT_DIR)

    summary = {
        "source": str(CSV),
        "question": "is the ~360 Hz peak found in the idle tail a mechanical "
                    "mode, and is it dominant?",
        "sample_rate_hz": fs, "nyquist_hz": fs / 2,
        "welch": {"window": WINDOW, "nperseg": NPERSEG,
                  "overlap_fraction": OVERLAP,
                  "frequency_resolution_hz": fs / NPERSEG,
                  "detrend_polynomial_order": DETREND_ORDER},
        "mode_band_hz": list(MODE_BAND), "floor_band_hz": list(FLOOR_BAND),
        "idle_tail_control": {"peak_hz": tail_peak_hz,
                              "band_vs_floor_ratio": tail_ratio},
        "by_mres": summary_mres,
        "aliasing_caveat": (
            f"Nyquist is {fs/2:.0f} Hz. A measured {tail_peak_hz:.0f} Hz is "
            f"indistinguishable from {fs - tail_peak_hz:.0f} Hz or "
            f"{fs + tail_peak_hz:.0f} Hz at this sample rate; the ring-down "
            "establishes that it is a real free-decay frequency of the "
            "mechanics, NOT which absolute frequency that is."),
        "outputs": {
            "states": "rendered_assets/psd_by_motion_state.png",
            "ramp_grid": "rendered_assets/psd_ramp_grid.png",
            "ringdown": "rendered_assets/ringdown_fit.png",
        },
    }
    (V23_DIR / "spectral_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"wrote 3 figures -> {OUT_DIR}")
    print(f"wrote spectral_summary.json -> {V23_DIR}")


if __name__ == "__main__":
    main()
