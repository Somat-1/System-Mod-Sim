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
from scipy.signal import StateSpace, lsim

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


def highest_mode_period(A: np.ndarray) -> float:
    eig = np.linalg.eigvals(A)
    osc = [e for e in eig if e.imag > 1e-6]
    f_max = max(abs(e) / (2.0 * np.pi) for e in osc)
    return 1.0 / f_max


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
    dt = period_min / 10.0
    print(f"highest resonant mode period = {period_min*1e3:.4f} ms -> dt = {dt*1e6:.2f} us")

    step_full = full_step_angle(params)
    lead_ratio = params["L"] / (2.0 * np.pi)

    cases = {}  # (divisor_name, speed_name) -> dict(t, theta_cmd, y, edges)
    for divisor_name, divisor in MICROSTEP_DIVISORS.items():
        for speed_name, dwell in [("fast", DWELL_FAST), ("settled", DWELL_SETTLE)]:
            t, theta_cmd, edges = build_theta_cmd(divisor, dwell, dt, step_full)
            print(f"{divisor_name:8s} {speed_name:8s}: {len(t):8d} samples, "
                  f"{t[-1]*1e3:9.2f} ms span, divisor={divisor}")
            tout, yout = run_case(A, B, C_out, D, dt, theta_cmd)
            cases[(divisor_name, speed_name)] = dict(
                t=tout, theta_cmd=theta_cmd, y=yout, edges=edges, dwell=dwell,
            )
            print(f"  -> lsim done")

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
    for ax, key in zip(axes.flat, panel_order):
        d = cases[key]
        x_cmd_um = lead_ratio * d["theta_cmd"] * 1e6
        x_n_um = d["y"][:, 2] * 1e6
        ax.plot(d["t"] * 1e3, x_cmd_um, color="#999999", linewidth=1.0, linestyle="--", label="commanded (ideal)")
        ax.plot(d["t"] * 1e3, x_n_um, color="#2b6cb0", linewidth=1.1, label="x_n (actual)")
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
    dwell = d["dwell"]
    lo, hi = t0 - 0.15 * dwell, t0 + 1.05 * dwell
    mask = (d["t"] >= lo) & (d["t"] <= hi)
    tt = (d["t"][mask] - t0) * 1e3  # ms relative to the step edge

    theta_cmd_deg = np.rad2deg(d["theta_cmd"][mask])
    theta_m_deg = np.rad2deg(d["y"][mask, 0])
    x_s_um_zoom = d["y"][mask, 1] * 1e6
    x_n_um_zoom = d["y"][mask, 2] * 1e6

    fig2, (ax_rotor, ax_lag, ax_bounce) = plt.subplots(3, 1, figsize=(9.0, 10.0), sharex=False)

    ax_rotor.plot(tt, theta_cmd_deg, color="#999999", linewidth=1.2, linestyle="--", label="theta_cmd")
    ax_rotor.plot(tt, theta_m_deg, color="#2b6cb0", linewidth=1.2, label="theta_m")
    ax_rotor.set_title("Rotor Snap: command vs. rotor position")
    ax_rotor.set_xlabel("Time since step edge (ms)")
    ax_rotor.set_ylabel("Angle (deg)")
    ax_rotor.grid(True, linewidth=0.4, color="#cccccc")
    ax_rotor.legend(fontsize=8)

    ax_lag.plot(tt, x_n_um_zoom, color="#2b6cb0", linewidth=1.2)
    ax_lag.set_title("Structural Lag: stage carriage position (x_n), single step")
    ax_lag.set_xlabel("Time since step edge (ms)")
    ax_lag.set_ylabel("Position (µm)")
    ax_lag.grid(True, linewidth=0.4, color="#cccccc")

    ax_bounce.plot((d["t"]) * 1e3, d["y"][:, 1] * 1e6, color="#c05621", linewidth=0.9)
    for e in d["edges"]:
        ax_bounce.axvline(e * 1e3, color="#dddddd", linewidth=0.5, zorder=0)
    ax_bounce.set_title("Axial Bounce: screw axial compression (x_s), full sequence")
    ax_bounce.set_xlabel("Time (ms)")
    ax_bounce.set_ylabel("Position (µm)")
    ax_bounce.grid(True, linewidth=0.4, color="#cccccc")

    fig2.tight_layout()
    diag_path = ASSET_DIR / "stepping_diagnostics.svg"
    fig2.savefig(diag_path)
    fig2.savefig(diag_path.with_suffix(".png"), dpi=110)
    plt.close(fig2)
    print(f"Wrote {diag_path}")


if __name__ == "__main__":
    main()
