#!/usr/bin/env python3
"""Stepping-sequence time-domain simulation for the Rev 4 6-DOF model.

Pipeline: build theta_cmd(t) as a full-step-equivalent staircase at a given
microstep divisor -> u(t) = [theta_cmd, T_fric_sb=0, F_fric_way=0] (frictionless
test) -> scipy.signal.lsim against (A, B, C_out, D) from build_bode_rev4's
state-space assembly, with

    C_out row 1 -> theta_m   (Motor Rotor Position)
    C_out row 2 -> x_s       (Screw Axial Compression)
    C_out row 3 -> x_n       (Stage Carriage Position)
    D = 0(3x3)

Four runs: {full step, 16x microstep} x {fast firing, settled firing}. The
move-count sequence is fired in full-step-equivalent units at every divisor
(a "2 steps forth" move is 2 elementary pulses at divisor=1, or 32 finer
pulses at divisor=16, same net travel and same total move duration either
way -- only the graduation gets finer).

Outputs two figures in rendered_assets/:
  stepping_montage.svg     -- 4-panel macro overview (all four runs)
  stepping_diagnostics.svg -- 3-panel zoom (rotor snap / structural lag /
                               axial bounce), from the full-step, settled run
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import StateSpace, find_peaks, lsim

from build_bode_rev4 import STATE_LABELS, build_matrices, build_state_space, load_parameters

ROOT = Path(__file__).resolve().parent
ASSET_DIR = ROOT / "rendered_assets"

# Full-step-equivalent move sequence: 2 forth, 1 back, 1 forth, 1 back, 1 forth,
# 4 back, 1 forth, 1 back, 1 forth, 1 back, 2 forth. Net zero (ends on the same line).
MOVE_SEQUENCE = [2, -1, 1, -1, 1, -4, 1, -1, 1, -1, 2]
assert sum(MOVE_SEQUENCE) == 0

DWELL_FAST = 4.0e-3     # s, per full-step-equivalent move -- short vs. the ~176 ms
                        # settling time of the slowest mode, so ringing overlaps
DWELL_SETTLE = 250.0e-3  # s, per full-step-equivalent move -- comfortably longer
                        # than the ~176 ms settling time of the slowest mode

MICROSTEP_DIVISORS = {"full": 1, "micro16": 16}

DT = 5.0e-6  # 145.7 us / 5 us ~= 29 samples/cycle of the fastest mode (6863 Hz).
             # Also divides both elementary dwells exactly (250 ms full-step ->
             # 50000 samples, 15.625 ms 16x -> 3125 samples), so every step edge
             # lands on a sample instead of being smeared into a ramp of duration
             # dt by lsim's linear interpolation.

NAMES = ["theta_m", "x_s", "x_n"]  # C_out row order, see build_C_out_D


def highest_mode_period(A: np.ndarray) -> float:
    eig = np.linalg.eigvals(A)
    osc = [e for e in eig if e.imag > 1e-6]
    f_max = max(abs(e) / (2.0 * np.pi) for e in osc)
    return 1.0 / f_max


def all_mode_properties(A: np.ndarray) -> list[tuple[float, float]]:
    """Distinct oscillatory modes of A as (natural frequency Hz, damping
    ratio), sorted ascending by frequency -- index 0 is mode 1. Frequency is
    |lambda|/(2*pi) (matches all_mode_frequencies); damping ratio is
    -Re(lambda)/|lambda|."""
    eig = np.linalg.eigvals(A)
    osc = [e for e in eig if e.imag > 1e-6]
    props = {}
    for e in osc:
        f = round(abs(e) / (2.0 * np.pi), 3)
        props[f] = -e.real / abs(e)
    return sorted(props.items())


def all_mode_frequencies(A: np.ndarray) -> np.ndarray:
    """Distinct oscillatory eigenfrequencies of A, in Hz. Any FFT peak of the
    time-domain output that does not land on one of these is aliasing, not
    physics -- there is no other source of oscillation in a linear model."""
    return np.array([f for f, _ in all_mode_properties(A)])


def grid_halving(A, B, C_out, D, divisor, dwell, step_full, dt) -> np.ndarray:
    """Run at dt and dt/2, compare on the dt grid. Returns per-channel max
    relative error (NAMES order). This is the decisive convergence check --
    if it's not near the solver noise floor, the plot isn't trustworthy
    regardless of how it looks."""
    t1, theta_cmd1, _ = build_theta_cmd(divisor, dwell, dt, step_full)
    _, y1 = run_case(A, B, C_out, D, dt, theta_cmd1)
    t2, theta_cmd2, _ = build_theta_cmd(divisor, dwell, dt / 2.0, step_full)
    _, y2 = run_case(A, B, C_out, D, dt / 2.0, theta_cmd2)
    n_common = min(len(y1), len(y2[::2]))
    y1c = y1[:n_common]
    y2c = y2[::2][:n_common]
    return np.max(np.abs(y1c - y2c), axis=0) / np.max(np.abs(y1c), axis=0)


def compute_spectrum(x_s_full, edge_idx, dt, window_s=0.05):
    """Hanning-windowed FFT of a window_s segment of x_s starting at edge_idx.
    Returns (f, mag) -- the full spectrum, not just samples at known
    eigenfrequencies, so it can be plotted as evidence rather than only
    quoted as a table of ratios."""
    n = int(round(window_s / dt))
    seg = x_s_full[edge_idx:edge_idx + n]
    f = np.fft.rfftfreq(len(seg), dt)
    mag = np.abs(np.fft.rfft(seg * np.hanning(len(seg))))
    return f, mag


def fft_against_eigenfrequencies(A, f, mag):
    """Sample a precomputed spectrum (f, mag) at each known structural
    eigenfrequency, normalized by the spectrum's peak magnitude. Confirms
    nothing aliased (every real peak lands on a known eigenfrequency) and
    identifies which mode actually dominates x_s."""
    freqs_hz = all_mode_frequencies(A)
    mag_max = mag.max()
    return [(fm, mag[np.argmin(np.abs(f - fm))] / mag_max) for fm in freqs_hz]


def full_step_angle(p: dict[str, float]) -> float:
    """One full motor step, per Sec. 3.2: 4*N_r detent cycles per revolution."""
    return 2.0 * np.pi / (4.0 * p["N_r"])


def build_theta_cmd(divisor: int, dwell_full_step: float, dt: float, step_full: float) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Piecewise-constant theta_cmd(t) staircase for MOVE_SEQUENCE at the given
    microstep divisor. One dwell_full_step of pre-roll at zero and one of
    tail hold are added. Returns (t, theta_cmd, elementary_edge_times)."""
    step_elem = step_full / divisor
    dwell_elem = dwell_full_step / divisor

    edge_times: list[float] = []
    seg_starts = [0.0]
    seg_values = [0.0]

    t_cursor = dwell_full_step  # pre-roll at zero
    theta_cursor = 0.0
    for move in MOVE_SEQUENCE:
        direction = 1.0 if move > 0 else -1.0
        for _ in range(abs(move) * divisor):
            edge_times.append(t_cursor)
            theta_cursor += direction * step_elem
            seg_starts.append(t_cursor)
            seg_values.append(theta_cursor)
            t_cursor += dwell_elem

    t_total = t_cursor + dwell_full_step  # tail hold
    t = np.arange(0.0, t_total, dt)
    seg_starts_arr = np.array(seg_starts)
    seg_values_arr = np.array(seg_values)
    idx = np.searchsorted(seg_starts_arr, t, side="right") - 1
    idx = np.clip(idx, 0, len(seg_values_arr) - 1)
    theta_cmd = seg_values_arr[idx]
    return t, theta_cmd, edge_times


def build_C_out_D() -> tuple[np.ndarray, np.ndarray]:
    C_out = np.zeros((3, 12))
    C_out[0, STATE_LABELS.index("theta_m")] = 1.0
    C_out[1, STATE_LABELS.index("x_s")] = 1.0
    C_out[2, STATE_LABELS.index("x_n")] = 1.0
    D = np.zeros((3, 3))
    return C_out, D


def run_case(A, B, C_out, D, dt, theta_cmd):
    n = len(theta_cmd)
    t = np.arange(n) * dt
    U = np.zeros((n, 3))
    U[:, 0] = theta_cmd  # T_fric_sb, F_fric_way stay zero: frictionless test
    sys = StateSpace(A, B, C_out, D)
    tout, yout, _xout = lsim(sys, U=U, T=t)
    return tout, yout


def main() -> None:
    params = load_parameters()
    M, C, K, B_u = build_matrices(params)
    A, B, _Cy = build_state_space(M, C, K, B_u)
    C_out, D = build_C_out_D()

    period_min = highest_mode_period(A)
    dt = DT
    print(f"highest resonant mode period = {period_min*1e3:.4f} ms "
          f"({1.0/period_min:.1f} Hz) -> dt = {dt*1e6:.2f} us "
          f"({period_min/dt:.2f} samples/cycle)")

    step_full = full_step_angle(params)
    lead_ratio = params["L"] / (2.0 * np.pi)

    cases = {}  # (divisor_name, speed_name) -> dict(t, theta_cmd, y, edges)
    for divisor_name, divisor in MICROSTEP_DIVISORS.items():
        for speed_name, dwell in [("fast", DWELL_FAST), ("settled", DWELL_SETTLE)]:
            # Edge-alignment check: the elementary step dwell must divide dt
            # exactly, or lsim's linear interpolation smears that edge into a
            # ramp of duration dt at a phase-dependent offset, so steps stop
            # being identical excitations.
            dwell_elem = dwell / divisor
            n_edge = dwell_elem / dt
            assert abs(n_edge - round(n_edge)) < 1e-9, (
                f"{divisor_name}/{speed_name}: elementary dwell {dwell_elem*1e6:.4f} us "
                f"does not divide dt={dt*1e6:.4f} us evenly (lands at sample {n_edge:.4f})"
            )

            t, theta_cmd, edges = build_theta_cmd(divisor, dwell, dt, step_full)
            print(f"{divisor_name:8s} {speed_name:8s}: {len(t):8d} samples, "
                  f"{t[-1]*1e3:9.2f} ms span, divisor={divisor}")
            tout, yout = run_case(A, B, C_out, D, dt, theta_cmd)
            cases[(divisor_name, speed_name)] = dict(
                t=tout, theta_cmd=theta_cmd, y=yout, edges=edges, dwell=dwell,
            )
            print(f"  -> lsim done")

            # Grid halving: converged when halving dt moves nothing (rel error
            # near the solver noise floor). This is the decisive check -- a
            # plot can look fine and still not be resolved.
            rel = grid_halving(A, B, C_out, D, divisor, dwell, step_full, dt)
            print(f"  -> dt vs dt/2 grid halving:")
            for name, e in zip(NAMES, rel):
                print(f"       {name:10s} {e:.2e}")

    # Supplementary: was the original coarse default (50 us) already adequate?
    # If its grid-halving error is also at the noise floor, the honest finding
    # is "the coarse grid was fine, the original plot was a rendering
    # artifact" rather than a claim that DT=5us fixed an integration error.
    print("\n50 us vs 25 us grid halving, full/settled case (was the original default adequate?):")
    rel_50 = grid_halving(A, B, C_out, D, MICROSTEP_DIVISORS["full"], DWELL_SETTLE, step_full, 50e-6)
    for name, e in zip(NAMES, rel_50):
        print(f"  {name:10s} {e:.2e}")

    # FFT check: is the beat visible in x_s real (all energy on known
    # eigenfrequencies) or an artifact of a mismatched sampling/beat period?
    # Also identifies which mode actually dominates x_s, which is what the
    # diagnostics-figure caption below reports samples/period for. Reuses the
    # full/settled case already computed above rather than rerunning lsim.
    ref = cases[("full", "settled")]
    ref_edge_idx = int(round(ref["edges"][0] / dt))
    f_spec, mag_spec = compute_spectrum(ref["y"][:, 1], ref_edge_idx, dt)
    fft_results = fft_against_eigenfrequencies(A, f_spec, mag_spec)
    print("\nFFT vs eigenfrequencies, 50 ms segment after first step edge, x_s, full/settled case:")
    for fm, ratio in fft_results:
        print(f"  {fm:8.2f} Hz -> {ratio:.3f}")
    # "Dominant" by peak ratio isn't what sets the sampling requirement -- the
    # fastest mode that's actually present (above a noise threshold) does.
    # Modes with negligible ratio here (e.g. 5/6, the ones the original DT=50us
    # concern was about) don't participate in x_s and impose no requirement.
    SIGNIFICANCE_THRESHOLD = 0.05
    significant = [(fm, r) for fm, r in fft_results if r > SIGNIFICANCE_THRESHOLD]
    print(f"  significant modes in x_s (ratio > {SIGNIFICANCE_THRESHOLD}): "
          + ", ".join(f"{fm:.1f} Hz ({r:.3f})" for fm, r in significant))
    fastest_significant_freq, fastest_significant_ratio = max(significant, key=lambda r: r[0])
    print(f"  fastest significant mode: {fastest_significant_freq:.2f} Hz "
          f"(ratio {fastest_significant_ratio:.3f}) -- sets the sampling requirement for x_s")

    ASSET_DIR.mkdir(exist_ok=True)

    # ---- Figure A: 4-panel macro montage ----
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0))
    panel_order = [("full", "fast"), ("full", "settled"), ("micro16", "fast"), ("micro16", "settled")]
    titles = {
        ("full", "fast"): "Full step -- fast firing",
        ("full", "settled"): "Full step -- settled",
        ("micro16", "fast"): "16x microstep -- fast firing",
        ("micro16", "settled"): "16x microstep -- settled",
    }
    RENDER_POINTS_TARGET = 5000  # integrate at DT, decimate only at render time
    for ax, key in zip(axes.flat, panel_order):
        d = cases[key]
        stride = max(1, len(d["t"]) // RENDER_POINTS_TARGET)
        t_ms = d["t"][::stride] * 1e3
        x_cmd_um = (lead_ratio * d["theta_cmd"] * 1e6)[::stride]
        x_n_um = (d["y"][:, 2] * 1e6)[::stride]
        ax.plot(t_ms, x_cmd_um, color="#999999", linewidth=1.0, linestyle="--", label="commanded (ideal)")
        ax.plot(t_ms, x_n_um, color="#2b6cb0", linewidth=1.1, label="x_n (actual)")
        ax.set_title(titles[key])
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Position (µm)")
        ax.grid(True, linewidth=0.4, color="#cccccc")
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Rev 4 stepping sequence: 2f,1b,1f,1b,1f,4b,1f,1b,1f,1b,2f (net zero)")
    fig.tight_layout()
    montage_path = ASSET_DIR / "stepping_montage.svg"
    fig.savefig(montage_path)
    fig.savefig(montage_path.with_suffix(".png"), dpi=110)
    plt.close(fig)
    print(f"Wrote {montage_path}")

    # ---- Figure B: 3-panel diagnostics, full-step + settled run, first step edge ----
    d = cases[("full", "settled")]
    t0 = d["edges"][0]
    t1 = d["edges"][1]  # next edge -- isolation boundary so "a single full step" is literal
    dwell = d["dwell"]

    # Single-step isolation window: a little pre-roll context, then strictly
    # before the next edge. The previous window (t0 + 1.05*dwell) ran 12.5 ms
    # past t1, so "settled"/"peak"/"target" silently mixed in the onset of
    # the second step -- e.g. target read as 3.6 deg (two steps) instead of
    # 1.8 deg. Confirmed against edges[1] before trusting any number below.
    lo, hi = t0 - 0.15 * dwell, t1
    mask = (d["t"] >= lo) & (d["t"] < hi)
    tt = (d["t"][mask] - t0) * 1e3  # ms relative to the step edge

    theta_cmd_deg = np.rad2deg(d["theta_cmd"][mask])
    theta_m_deg = np.rad2deg(d["y"][mask, 0])
    x_n_zoom = d["y"][mask, 2]
    e_um = (lead_ratio * d["theta_cmd"][mask] - x_n_zoom) * 1e6

    target_deg = theta_cmd_deg[-1]
    settled_deg = theta_m_deg[-1]
    peak_idx = int(np.argmax(theta_m_deg))
    peak_deg = theta_m_deg[peak_idx]
    peak_t_ms = tt[peak_idx]
    overshoot_pct = (peak_deg - settled_deg) / settled_deg * 100.0

    mode_props = all_mode_properties(A)
    freq1_hz, zeta1 = mode_props[0]  # mode 1 -- lowest frequency, dominates this ringing

    settle_tol = 0.02 * abs(target_deg)
    outside_band = np.where(np.abs(theta_m_deg - target_deg) > settle_tol)[0]
    settle_idx = outside_band[-1] + 1 if len(outside_band) else 0
    settle_t_ms = tt[settle_idx] if settle_idx < len(tt) else tt[-1]

    peak_idxs, _ = find_peaks(theta_m_deg)
    p1_t = p2_t = period_ms = None
    if len(peak_idxs) >= 3:
        # peaks[0] is the initial big-overshoot peak; use two later, more
        # regular cycles for a clean period read.
        p1_t, p2_t = tt[peak_idxs[1]], tt[peak_idxs[2]]
        period_ms = p2_t - p1_t

    # Fixed, not autoscaled: sized for the current (post-detent-fix) peaks
    # with headroom, so a before/after B_u comparison stays visible instead
    # of each run silently rescaling to fill the frame.
    ROTOR_XLIM = ((lo - t0) * 1e3, (hi - t0) * 1e3)
    ROTOR_YLIM = (-0.3, 4.4)  # extra headroom above the peak for the stacked annotations
    ERROR_YLIM = (-6.5, 6.5)
    BOUNCE_XLIM = (-1.0, 20.0)
    BOUNCE_YLIM = (-0.7, 0.7)

    fig2, (ax_rotor, ax_err, ax_bounce) = plt.subplots(3, 1, figsize=(9.5, 13.0), sharex=False)

    # --- Panel 1: rotor response, annotated with the numbers actually quoted ---
    # Everything below the peak (y <= peak_deg) is data-occupied at some x in
    # this window (peak_deg is the global max), so every label here is parked
    # in the strip above it (y > peak_deg), stacked to avoid stacking on
    # each other: overshoot label lowest, period bracket above that, the
    # panel-3-window marker highest (it spans the widest x-range).
    ax_rotor.plot(tt, theta_cmd_deg, color="#999999", linewidth=1.2, linestyle="--", label="theta_cmd")
    ax_rotor.plot(tt, theta_m_deg, color="#2b6cb0", linewidth=1.2, label="theta_m")

    ax_rotor.axhline(target_deg, color="#333333", linestyle=":", linewidth=1.0, zorder=1)
    ax_rotor.text(ROTOR_XLIM[1], target_deg, f" target = {target_deg:.2f} deg", fontsize=7.5,
                  color="#333333", va="bottom", ha="right", clip_on=False)

    overshoot_y = peak_deg + 0.18
    ax_rotor.annotate("", xy=(peak_t_ms, peak_deg), xytext=(peak_t_ms, target_deg),
                       arrowprops=dict(arrowstyle="<->", color="#c05621", linewidth=1.1))
    ax_rotor.text(peak_t_ms, overshoot_y,
                  f"{overshoot_pct:.0f}% overshoot (ζ₁ = {zeta1:.3f})",
                  fontsize=7.5, color="#c05621", ha="center", va="bottom")

    if period_ms is not None:
        y_bracket = peak_deg + 0.45
        ax_rotor.annotate("", xy=(p2_t, y_bracket), xytext=(p1_t, y_bracket),
                           arrowprops=dict(arrowstyle="<->", color="#555555", linewidth=1.0))
        ax_rotor.text((p1_t + p2_t) / 2.0, y_bracket + 0.06,
                      f"{period_ms:.2f} ms = {freq1_hz:.1f} Hz (mode 1)",
                      fontsize=7.5, color="#555555", ha="center", va="bottom")

    ax_rotor.axhspan(target_deg * (1.0 - 0.02), target_deg * (1.0 + 0.02),
                      color="#2b6cb0", alpha=0.08, zorder=0)
    ax_rotor.axvline(settle_t_ms, color="#555555", linestyle=":", linewidth=1.0)
    ax_rotor.text(settle_t_ms + 3.0, ROTOR_YLIM[0] + 0.15,
                  f"settled (±2%)\nat {settle_t_ms:.0f} ms", fontsize=7.5, color="#555555")

    # Marks where panel 3's own (much tighter) time axis sits within this
    # one -- placed above the axis frame itself, not inside the data area.
    ax_rotor.axvspan(BOUNCE_XLIM[0], BOUNCE_XLIM[1], color="#c05621", alpha=0.08, zorder=0)
    win_trans = ax_rotor.get_xaxis_transform()  # x in data coords, y in axes fraction
    ax_rotor.annotate("", xy=(BOUNCE_XLIM[0], 1.025), xytext=(BOUNCE_XLIM[1], 1.025),
                       xycoords=win_trans, textcoords=win_trans,
                       arrowprops=dict(arrowstyle="<->", color="#c05621", linewidth=1.0),
                       annotation_clip=False)
    ax_rotor.text((BOUNCE_XLIM[0] + BOUNCE_XLIM[1]) / 2.0, 1.045, "panel 3 window",
                  transform=win_trans, ha="center", va="bottom", fontsize=7.5,
                  color="#c05621", clip_on=False)

    ax_rotor.set_title("Rotor response to a single full step", pad=30)
    ax_rotor.set_xlabel("Time since step edge (ms)")
    ax_rotor.set_ylabel("Angle (deg)")
    ax_rotor.set_xlim(*ROTOR_XLIM)
    ax_rotor.set_ylim(*ROTOR_YLIM)
    ax_rotor.grid(True, linewidth=0.4, color="#cccccc")
    ax_rotor.legend(fontsize=8, loc="upper right")

    # --- Panel 2: stage tracking error e(t) = (L/2pi)*theta_cmd - x_n. The
    # residual, not the nearly-overlapping position traces it replaces --
    # decays to ~0 here (frictionless); with friction live this is where a
    # non-zero dead-band would show up. Shares panel 1's x-axis exactly. ---
    ax_err.plot(tt, e_um, color="#2b6cb0", linewidth=1.1)
    ax_err.axhline(0.0, color="#333333", linestyle=":", linewidth=1.0)
    ax_err.text(ROTOR_XLIM[1], 0.0, f" e -> {e_um[-1]:.3f} µm", fontsize=7.5, color="#333333",
                va="bottom", ha="right", clip_on=False)
    ax_err.set_title(r"Stage tracking error, $e(t) = (L/2\pi)\,\theta_{cmd}(t) - x_n(t)$")
    ax_err.set_xlabel("Time since step edge (ms)")
    ax_err.set_ylabel("Error (µm)")
    ax_err.set_xlim(*ROTOR_XLIM)
    ax_err.set_ylim(*ERROR_YLIM)
    ax_err.grid(True, linewidth=0.4, color="#cccccc")

    # Inset: +-6.5 um hides the residual entirely (it's ~500x smaller than
    # the ring amplitude). Zoom to +-50 nm over the last 100 ms -- the actual
    # number this panel will get compared against once friction is live and
    # this stops decaying to zero.
    INSET_YLIM_UM = (-0.05, 0.05)  # +-50 nm
    inset_lo_ms = tt[-1] - 100.0
    inset_mask = tt >= inset_lo_ms
    ax_err_inset = ax_err.inset_axes([0.60, 0.60, 0.36, 0.34])
    ax_err_inset.plot(tt[inset_mask], e_um[inset_mask], color="#2b6cb0", linewidth=1.0)
    ax_err_inset.axhline(0.0, color="#333333", linestyle=":", linewidth=0.7)
    ax_err_inset.set_xlim(inset_lo_ms, tt[-1])
    ax_err_inset.set_ylim(*INSET_YLIM_UM)
    ax_err_inset.set_title(f"last 100 ms, ±50 nm (settles to {e_um[-1]*1000:.0f} nm)", fontsize=6.5)
    ax_err_inset.tick_params(labelsize=6)
    ax_err_inset.grid(True, linewidth=0.3, color="#cccccc")
    ax_err.indicate_inset_zoom(ax_err_inset, edgecolor="#888888", alpha=0.6, linewidth=0.8)

    # --- Panel 3: screw axial deflection at the thrust bearing. x_s is the
    # screw body's translation against the grounded k_brg preload spring, not
    # compression of the screw shaft -- different physical quantity. ---
    x_s_um_full = d["y"][:, 1] * 1e6
    peak_x_s_um = np.max(np.abs(x_s_um_full))
    # ~20 ms post-edge: at 2.5 ms the trace is still growing (Fig. montage);
    # 20 ms gives ~15 cycles of the fastest significant mode found by the FFT
    # check above and covers its decay envelope, so the panel shows what's
    # actually there instead of a truncated transient.
    bounce_lo, bounce_hi = t0 + BOUNCE_XLIM[0] * 1e-3, t0 + BOUNCE_XLIM[1] * 1e-3
    bounce_mask = (d["t"] >= bounce_lo) & (d["t"] <= bounce_hi)
    samples_per_period_fast = (1.0 / fastest_significant_freq) / dt
    ax_bounce.plot((d["t"][bounce_mask] - t0) * 1e3, x_s_um_full[bounce_mask], color="#c05621", linewidth=0.9)
    ax_bounce.axvline(0.0, color="#dddddd", linewidth=0.5, zorder=0)
    ax_bounce.set_title(
        f"Screw axial deflection at the thrust bearing (x_s), ~20 ms window\n"
        f"dt = {dt*1e6:.2f} us, {samples_per_period_fast:.0f} samples/period @ {fastest_significant_freq:.0f} Hz "
        f"| peak = {peak_x_s_um:.3f} µm"
    )
    ax_bounce.set_xlabel("Time since step edge (ms)")
    ax_bounce.set_ylabel("Deflection (µm)")
    ax_bounce.set_xlim(*BOUNCE_XLIM)
    ax_bounce.set_ylim(*BOUNCE_YLIM)
    ax_bounce.grid(True, linewidth=0.4, color="#cccccc")

    # k_brg is a single value in model_parameters.json but the source spec
    # (Barden duplex, light preload, Sec. 8 of state_space_6dof.md) gives a
    # range, not a point value -- the N-axis inherits that same uncertainty.
    K_BRG_RANGE_N_PER_M = (7.5e6, 15.3e6)
    k_brg = params["k_brg"]
    ax_bounce_n = ax_bounce.twinx()
    ax_bounce_n.set_ylim(BOUNCE_YLIM[0] * k_brg * 1e-6, BOUNCE_YLIM[1] * k_brg * 1e-6)
    ax_bounce_n.set_ylabel(f"F = k_brg·x_s (N), k_brg = {k_brg/1e6:.1f} MN/m", fontsize=9, labelpad=8)

    peak_force_lo_n = peak_x_s_um * K_BRG_RANGE_N_PER_M[0] * 1e-6
    peak_force_hi_n = peak_x_s_um * K_BRG_RANGE_N_PER_M[1] * 1e-6
    fig2.tight_layout(rect=[0.0, 0.05, 1.0, 1.0])
    fig2.subplots_adjust(right=0.90, top=0.93, hspace=0.55)
    fig2.text(
        0.5, 0.026,
        f"Panel 3 N-axis uses k_brg = {k_brg/1e6:.1f} MN/m (model_parameters.json), a midpoint of the "
        f"Barden duplex light-preload spec's {K_BRG_RANGE_N_PER_M[0]/1e6:.1f}–{K_BRG_RANGE_N_PER_M[1]/1e6:.1f} MN/m range.",
        ha="center", va="bottom", fontsize=7.5, color="#555555",
    )
    fig2.text(
        0.5, 0.008,
        f"The force reading inherits that same uncertainty -- the peak deflection shown "
        f"({peak_x_s_um:.3f} µm) is {peak_force_lo_n:.1f}–{peak_force_hi_n:.1f} N across the range.",
        ha="center", va="bottom", fontsize=7.5, color="#555555",
    )
    diag_path = ASSET_DIR / "stepping_diagnostics.svg"
    fig2.savefig(diag_path)
    fig2.savefig(diag_path.with_suffix(".png"), dpi=110)
    plt.close(fig2)
    print(f"Wrote {diag_path}")

    # ---- Figure C: x_s spectrum, evidence for the "fastest significant mode"
    # claim in Figure B's caption -- the actual FFT, not just the sampled
    # ratio table it's derived from.
    fig3, ax_spec = plt.subplots(figsize=(9.5, 5.5))
    mag_db = 20.0 * np.log10(mag_spec / mag_spec.max() + 1e-12)
    ax_spec.plot(f_spec, mag_db, color="#2b6cb0", linewidth=0.8)
    f_max_plot = 1.15 * max(fm for fm, _ in fft_results)
    ax_spec.set_ylim(-80.0, 20.0)
    for fm, ratio in fft_results:
        is_fastest_sig = abs(fm - fastest_significant_freq) < 1e-6
        is_significant = ratio > SIGNIFICANCE_THRESHOLD
        color = "#c05621" if is_fastest_sig else ("#2b6cb0" if is_significant else "#999999")
        ax_spec.axvline(fm, color=color, linestyle="--", linewidth=1.2 if is_significant else 0.8)
        ax_spec.text(fm, 3.0, f"{fm:.0f} Hz", rotation=90, va="bottom", ha="center",
                     fontsize=8, color=color, clip_on=True)
    ax_spec.axhline(20.0 * np.log10(SIGNIFICANCE_THRESHOLD), color="#888888", linestyle=":",
                     linewidth=1.0, label=f"significance threshold (ratio={SIGNIFICANCE_THRESHOLD})")
    ax_spec.set_xlim(0.0, f_max_plot)
    ax_spec.set_xlabel("Frequency (Hz)")
    ax_spec.set_ylabel("Magnitude (dB, normalized to spectrum peak)")
    fig3.suptitle(
        f"x_s spectrum, 50 ms segment after first step edge (full/settled case)\n"
        f"Hanning-windowed FFT, dt = {dt*1e6:.2f} us -- dashed lines: known eigenfrequencies "
        f"(orange = fastest significant mode, {fastest_significant_freq:.0f} Hz)",
        fontsize=10,
    )
    ax_spec.grid(True, linewidth=0.4, color="#cccccc")
    ax_spec.legend(fontsize=8, loc="upper right")
    fig3.tight_layout(rect=[0.0, 0.0, 1.0, 0.88])
    spec_path = ASSET_DIR / "stepping_xs_spectrum.svg"
    fig3.savefig(spec_path)
    fig3.savefig(spec_path.with_suffix(".png"), dpi=110)
    plt.close(fig3)
    print(f"Wrote {spec_path}")


if __name__ == "__main__":
    main()
