#!/usr/bin/env python3
"""Rev 4.2 nonlinear stepping suite and frictionless comparison overlay.

Mirrors the Rev 4 command sequence and four montage cases:
  {full step, 16x microstep} x {4 ms, 250 ms firing interval}.
The nonlinear segments are integrated with Radau and the analytical Rev 4.2
Jacobian.  The frictionless baseline uses its exact LTI state-space model on
the same output grids.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp
from scipy.signal import StateSpace, lsim

from lugre_model_rev42 import (
    N_Q,
    N_STATES,
    PORTS,
    LuGreModelRev42,
    lugre_terms,
    lugre_terms_exact,
    _port_values,
)


ROOT = Path(__file__).resolve().parent.parent
REV4_DIR = ROOT.parents[1]
ASSET_DIR = ROOT / "rendered_assets"
NPZ_DIR = ASSET_DIR / "npz"
sys.path.insert(0, str(REV4_DIR / "scripts"))
from build_bode_rev4 import (  # noqa: E402
    INPUT_LABELS,
    build_matrices as build_baseline_matrices,
    build_state_space as build_baseline_state_space,
    load_parameters as load_baseline_parameters,
)


MOVES = [2, -1, 1, -1, 1, -4, 1, -1, 1, -1, 2]
assert sum(MOVES) == 0
FULL_STEP = np.deg2rad(1.8)
CASES = (
    ("full", 1, "settled", 250.0e-3),
    ("micro16", 16, "settled", 250.0e-3),
    ("full", 1, "fast", 4.0e-3),
    ("micro16", 16, "fast", 4.0e-3),
)
PLOT_DT = 100.0e-6
DIAGNOSTIC_DT = 5.0e-6
RTOL = 1.0e-6
ATOL = np.array([
    1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-12, 1.0e-12,
    1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9,
    1.0e-11, 1.0e-11, 1.0e-9,
])


def edge_list(micro: int) -> np.ndarray:
    theta = 0.0
    edges: list[float] = []
    for move in MOVES:
        for _ in range(abs(move) * micro):
            theta += np.sign(move) * FULL_STEP / micro
            edges.append(theta)
    return np.asarray(edges)


def _segment_grid(duration: float, dt: float) -> np.ndarray:
    count = int(round(duration / dt))
    if not np.isclose(count * dt, duration, rtol=0.0, atol=1.0e-14):
        raise ValueError(f"Segment duration {duration} is not divisible by dt={dt}")
    return np.linspace(0.0, duration, count + 1)


def sampled_interface_power(model: LuGreModelRev42, states: np.ndarray) -> np.ndarray:
    velocity = states[N_Q:2 * N_Q].T
    power = np.zeros(states.shape[1])
    for index, port in enumerate(PORTS):
        v = velocity @ model.jacobians[port]
        z = states[2 * N_Q + index]
        term_function = (
            lugre_terms_exact if model.regularization == 'exact' else lugre_terms
        )
        force = term_function(v, z, *_port_values(model.p, port))[0]
        power += force * v
    return power


def simulate_rev42(
    model: LuGreModelRev42,
    micro: int,
    firing_interval: float,
    output_dt: float,
    rtol: float = RTOL,
    progress: bool = True,
) -> dict[str, np.ndarray | float | int]:
    commands = edge_list(micro)
    local_time = _segment_grid(firing_interval, output_dt)
    state = np.zeros(N_STATES)
    s_parts: list[np.ndarray] = []
    t_parts: list[np.ndarray] = []
    command_parts: list[np.ndarray] = []
    state_parts: list[np.ndarray] = []
    nfev = njev = nlu = 0
    start = time.perf_counter()
    jacobian = (
        None if model.regularization == 'exact'
        else lambda _t, y: model.analytical_linearization(y)[0]
    )

    for segment, command in enumerate(commands):
        solution = solve_ivp(
            lambda t, y, cmd=command: model.rhs(t, y, cmd),
            (0.0, firing_interval),
            state,
            method="Radau",
            jac=jacobian,
            t_eval=local_time,
            rtol=rtol,
            atol=ATOL,
        )
        if not solution.success:
            raise RuntimeError(
                f"Rev 4.2 integration failed at segment {segment}: {solution.message}"
            )
        state = solution.y[:, -1]
        nfev += solution.nfev
        njev += solution.njev
        nlu += solution.nlu
        kept_time = local_time[:-1]
        s_parts.append((segment + kept_time / firing_interval) / micro)
        t_parts.append(segment * firing_interval + kept_time)
        command_parts.append(np.full(kept_time.size, command))
        state_parts.append(solution.y[:, :-1])
        if progress and ((segment + 1) % max(1, len(commands) // 8) == 0):
            print(
                f"    Rev 4.2: {segment + 1:3d}/{len(commands)} segments, "
                f"elapsed={time.perf_counter() - start:.1f}s",
                flush=True,
            )

    states = np.concatenate(state_parts, axis=1)
    command = np.concatenate(command_parts)
    time_values = np.concatenate(t_parts)
    lead = model.p["L"] / (2.0 * np.pi)
    power = sampled_interface_power(model, states)
    return {
        "s": np.concatenate(s_parts),
        "time_s": time_values,
        "theta_cmd": command,
        "theta_m": states[0],
        "x_s": states[4],
        "x_n": states[5],
        "error": lead * command - states[5],
        "interface_power_W": power,
        "interface_work_J": float(np.trapezoid(power, time_values)),
        "minimum_interface_power_W": float(np.min(power)),
        "final_state": state,
        "nfev": nfev,
        "njev": njev,
        "nlu": nlu,
        "elapsed_s": time.perf_counter() - start,
    }


def baseline_system() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    parameters = load_baseline_parameters()
    mass, damping, stiffness, input_matrix = build_baseline_matrices(parameters)
    A, B, _ = build_baseline_state_space(mass, damping, stiffness, input_matrix)
    output = np.zeros((3, A.shape[0]))
    output[0, 0] = 1.0
    output[1, 4] = 1.0
    output[2, 5] = 1.0
    feedthrough = np.zeros((3, B.shape[1]))
    lead = parameters["L"] / (2.0 * np.pi)
    return A, B, output, feedthrough, lead


def simulate_baseline(
    system_parts: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float],
    micro: int,
    firing_interval: float,
    output_dt: float,
) -> dict[str, np.ndarray]:
    A, B, output, feedthrough, lead = system_parts
    system = StateSpace(A, B, output, feedthrough)
    commands = edge_list(micro)
    local_time = _segment_grid(firing_interval, output_dt)
    state = np.zeros(A.shape[0])
    s_parts: list[np.ndarray] = []
    t_parts: list[np.ndarray] = []
    command_parts: list[np.ndarray] = []
    output_parts: list[np.ndarray] = []
    command_column = INPUT_LABELS.index("theta_cmd")

    for segment, command in enumerate(commands):
        inputs = np.zeros((local_time.size, B.shape[1]))
        inputs[:, command_column] = command
        _, values, states = lsim(system, U=inputs, T=local_time, X0=state)
        state = states[-1]
        kept_time = local_time[:-1]
        s_parts.append((segment + kept_time / firing_interval) / micro)
        t_parts.append(segment * firing_interval + kept_time)
        command_parts.append(np.full(kept_time.size, command))
        output_parts.append(values[:-1].T)

    values = np.concatenate(output_parts, axis=1)
    command = np.concatenate(command_parts)
    return {
        "s": np.concatenate(s_parts),
        "time_s": np.concatenate(t_parts),
        "theta_cmd": command,
        "theta_m": values[0],
        "x_s": values[1],
        "x_n": values[2],
        "error": lead * command - values[2],
    }


def single_step_rev42(
    model: LuGreModelRev42, output_dt: float, rtol: float
) -> dict[str, np.ndarray | int]:
    local_time = _segment_grid(250.0e-3, output_dt)
    command = FULL_STEP
    solution = solve_ivp(
        lambda t, y: model.rhs(t, y, command),
        (0.0, local_time[-1]),
        np.zeros(N_STATES),
        method="Radau",
        jac=lambda _t, y: model.analytical_linearization(y)[0],
        t_eval=local_time,
        rtol=rtol,
        atol=ATOL,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    lead = model.p["L"] / (2.0 * np.pi)
    power = sampled_interface_power(model, solution.y)
    return {
        "time_s": solution.t,
        "theta_cmd": np.full(solution.t.size, command),
        "theta_m": solution.y[0],
        "x_s": solution.y[4],
        "x_n": solution.y[5],
        "error": lead * command - solution.y[5],
        "interface_power_W": power,
        "nfev": solution.nfev,
    }


def single_step_baseline(
    system_parts: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float],
    output_dt: float,
) -> dict[str, np.ndarray]:
    A, B, output, feedthrough, lead = system_parts
    local_time = _segment_grid(250.0e-3, output_dt)
    inputs = np.zeros((local_time.size, B.shape[1]))
    inputs[:, INPUT_LABELS.index("theta_cmd")] = FULL_STEP
    _, values, _ = lsim(StateSpace(A, B, output, feedthrough), U=inputs, T=local_time)
    return {
        "time_s": local_time,
        "theta_cmd": np.full(local_time.size, FULL_STEP),
        "theta_m": values[:, 0],
        "x_s": values[:, 1],
        "x_n": values[:, 2],
        "error": lead * FULL_STEP - values[:, 2],
    }


def render_montage(cases: dict[tuple[str, str], dict], lead: float) -> Path:
    layout = (("settled", 250.0e-3), ("fast", 4.0e-3))
    columns = (("full", 1, "full step"), ("micro16", 16, "16x microstep"))
    fig, axes = plt.subplots(4, 2, figsize=(11.0, 13.0), sharex="col")
    for column, (name, micro, label) in enumerate(columns):
        for group, (speed, firing_interval) in enumerate(layout):
            case = cases[(name, speed)]
            s = case["s"]
            command_um = lead * case["theta_cmd"] * 1.0e6
            position_um = case["x_n"] * 1.0e6
            error_um = case["error"] * 1.0e6
            ax_position = axes[2 * group, column]
            ax_error = axes[2 * group + 1, column]
            linewidth = 1.1 if speed == "settled" else 0.65
            ax_position.plot(s, command_um, "--", color="#888888", linewidth=0.8)
            ax_position.plot(s, position_um, color="#c0392b", linewidth=linewidth)
            ax_error.plot(s, error_um, color="#c0392b", linewidth=linewidth)
            ax_error.axhline(0.0, color="#333333", linestyle=":", linewidth=0.7)
            for axis in (ax_position, ax_error):
                axis.grid(True, linewidth=0.4, color="#cccccc")
                axis.set_xlim(0.0, 16.0)
            ax_position.set_title(
                f"{label} — {speed} ({firing_interval * 1e3:.0f} ms)", fontsize=9
            )
            peak = np.max(np.abs(error_um))
            ax_error.text(
                0.98, 0.90, f"peak |e| = {peak:.3f} µm",
                transform=ax_error.transAxes, ha="right", va="top", fontsize=7.5,
            )
            if column == 0:
                ax_position.set_ylabel("Position (µm)")
                ax_error.set_ylabel("Error (µm)")
            if group == 1:
                ax_error.set_xlabel("Steps fired (full-step units)")
    fig.legend(
        handles=[
            Line2D([0], [0], color="#888888", linestyle="--", label="commanded"),
            Line2D([0], [0], color="#c0392b", label="Rev 4.2 nonlinear x_n"),
        ],
        loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Rev 4.2 nonlinear stepping sequence: 2f,1b,1f,1b,1f,4b,1f,1b,1f,1b,2f",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    path = ASSET_DIR / "stepping_montage_rev42.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def render_diagnostics(case: dict, lead: float) -> Path:
    time_ms = case["time_s"] * 1.0e3
    command_deg = np.rad2deg(case["theta_cmd"])
    rotor_deg = np.rad2deg(case["theta_m"])
    error_um = case["error"] * 1.0e6
    screw_um = case["x_s"] * 1.0e6
    power_mw = case["interface_power_W"] * 1.0e3
    fig, axes = plt.subplots(4, 1, figsize=(10.0, 13.0), sharex=True)
    axes[0].plot(time_ms, command_deg, "--", color="#888888", label="theta_cmd")
    axes[0].plot(time_ms, rotor_deg, color="#c0392b", label="theta_m")
    axes[0].set_ylabel("Angle (deg)")
    axes[0].set_title("Rev 4.2 response to one full step")
    axes[0].legend(fontsize=8)
    axes[1].plot(time_ms, error_um, color="#c0392b")
    axes[1].axhline(0.0, color="#333333", linestyle=":", linewidth=0.8)
    axes[1].set_ylabel("Error (µm)")
    axes[1].set_title(r"Stage tracking error, $(L/2\pi)\theta_{cmd}-x_n$")
    axes[2].plot(time_ms, screw_um, color="#c05621")
    axes[2].set_ylabel("x_s (µm)")
    axes[2].set_title("Screw axial deflection at the thrust bearing")
    axes[3].plot(time_ms, power_mw, color="#6b46c1")
    axes[3].axhline(0.0, color="#333333", linestyle=":", linewidth=0.8)
    axes[3].set_ylabel("Σ F_p v_p (mW)")
    axes[3].set_xlabel("Time since step edge (ms)")
    axes[3].set_title("Instantaneous interface power (negative = bristle energy return)")
    for axis in axes:
        axis.grid(True, linewidth=0.4, color="#cccccc")
        axis.set_xlim(0.0, 250.0)
    fig.tight_layout()
    path = ASSET_DIR / "stepping_diagnostics_rev42.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def render_overlay(
    nonlinear_cases: dict[tuple[str, str], dict],
    baseline_cases: dict[tuple[str, str], dict],
) -> Path:
    panels = (
        ("full", "settled", "Full step — 250 ms"),
        ("micro16", "settled", "16x microstep — 250 ms"),
        ("full", "fast", "Full step — 4 ms"),
        ("micro16", "fast", "16x microstep — 4 ms"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.5), sharex=True)
    for axis, (step_name, speed, title) in zip(axes.flat, panels):
        nonlinear = nonlinear_cases[(step_name, speed)]
        baseline = baseline_cases[(step_name, speed)]
        axis.plot(
            baseline["s"], baseline["error"] * 1.0e6,
            color="#2b6cb0", linewidth=1.0, label="frictionless baseline",
        )
        axis.plot(
            nonlinear["s"], nonlinear["error"] * 1.0e6,
            color="#c0392b", linewidth=0.9, label="Rev 4.2 nonlinear LuGre",
        )
        axis.axhline(0.0, color="#333333", linestyle=":", linewidth=0.7)
        axis.set_title(title)
        axis.set_xlim(0.0, 16.0)
        axis.set_xlabel("Steps fired (full-step units)")
        axis.set_ylabel("Tracking error (µm)")
        axis.grid(True, linewidth=0.4, color="#cccccc")
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("Stepping tracking-error overlay: Rev 4.2 vs frictionless baseline")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    path = ASSET_DIR / "stepping_overlay_frictionless.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def render_spectrum(case: dict, model: LuGreModelRev42) -> tuple[Path, list[float]]:
    window = case["time_s"] <= 50.0e-3
    signal = case["x_s"][window]
    frequency = np.fft.rfftfreq(signal.size, DIAGNOSTIC_DT)
    magnitude = np.abs(np.fft.rfft(signal * np.hanning(signal.size)))
    magnitude_db = 20.0 * np.log10(magnitude / np.max(magnitude) + 1.0e-14)
    A, _, _ = model.analytical_linearization(model.cruise_state(5.0e-3))
    eigenvalues = np.linalg.eigvals(A)
    modes = sorted(
        float(abs(value) / (2.0 * np.pi))
        for value in eigenvalues if np.imag(value) > 1.0e-5
    )
    fig, axis = plt.subplots(figsize=(10.0, 5.8))
    axis.plot(frequency, magnitude_db, color="#c0392b", linewidth=0.8)
    for mode in modes:
        axis.axvline(mode, color="#777777", linestyle="--", linewidth=0.8)
        axis.text(mode, 2.0, f"{mode:.0f}", rotation=90, ha="center", va="bottom", fontsize=7)
    axis.set_xlim(0.0, 8000.0)
    axis.set_ylim(-100.0, 10.0)
    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel("Magnitude (dB, normalized)")
    axis.set_title("Rev 4.2 x_s spectrum, first 50 ms after a full step")
    axis.grid(True, linewidth=0.4, color="#cccccc")
    fig.tight_layout()
    path = ASSET_DIR / "stepping_xs_spectrum_rev42.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path, modes


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    model = LuGreModelRev42(enforce_interface_power=False)
    exact_model = LuGreModelRev42(
        enforce_interface_power=False, regularization='exact'
    )
    baseline_parts = baseline_system()
    lead = model.p["L"] / (2.0 * np.pi)
    nonlinear_cases: dict[tuple[str, str], dict] = {}
    exact_cases: dict[tuple[str, str], dict] = {}
    baseline_cases: dict[tuple[str, str], dict] = {}
    total_start = time.perf_counter()

    for step_name, micro, speed, firing_interval in CASES:
        print(
            f"Running {step_name}/{speed}: micro={micro}, "
            f"firing={firing_interval * 1e3:.1f} ms",
            flush=True,
        )
        nonlinear = simulate_rev42(model, micro, firing_interval, PLOT_DT)
        exact = simulate_rev42(exact_model, micro, firing_interval, PLOT_DT)
        baseline = simulate_baseline(baseline_parts, micro, firing_interval, PLOT_DT)
        nonlinear_cases[(step_name, speed)] = nonlinear
        exact_cases[(step_name, speed)] = exact
        baseline_cases[(step_name, speed)] = baseline
        print(
            f"  done in {nonlinear['elapsed_s']:.1f}s; "
            f"peak |e|={np.max(np.abs(nonlinear['error'])) * 1e6:.4f} µm; "
            f"net interface work={nonlinear['interface_work_J']:.4e} J",
            flush=True,
        )

    print("Running 5 µs single-step diagnostic and convergence check...", flush=True)
    diagnostic = single_step_rev42(model, DIAGNOSTIC_DT, RTOL)
    diagnostic_strict = single_step_rev42(model, DIAGNOSTIC_DT, 3.0e-7)
    diagnostic_baseline = single_step_baseline(baseline_parts, DIAGNOSTIC_DT)
    convergence = {}
    for key in ("theta_m", "x_s", "x_n"):
        scale = max(float(np.max(np.abs(diagnostic_strict[key]))), 1.0e-30)
        convergence[key] = float(np.max(np.abs(diagnostic[key] - diagnostic_strict[key])) / scale)

    montage_path = render_montage(nonlinear_cases, lead)
    diagnostic_path = render_diagnostics(diagnostic, lead)
    overlay_path = render_overlay(nonlinear_cases, baseline_cases)
    spectrum_path, modes = render_spectrum(diagnostic, model)

    data: dict[str, np.ndarray] = {}
    for prefix, collection in (("rev42", nonlinear_cases), ("baseline", baseline_cases)):
        for (step_name, speed), case in collection.items():
            key_prefix = f"{prefix}_{step_name}_{speed}"
            for field in ("s", "time_s", "theta_cmd", "theta_m", "x_s", "x_n", "error"):
                data[f"{key_prefix}_{field}"] = np.asarray(case[field])
    for prefix, case in (("diagnostic_rev42", diagnostic), ("diagnostic_baseline", diagnostic_baseline)):
        for field in ("time_s", "theta_cmd", "theta_m", "x_s", "x_n", "error"):
            data[f"{prefix}_{field}"] = np.asarray(case[field])
    data_path = NPZ_DIR / "stepping_rev42_and_frictionless.npz"
    for (step_name, speed), case in exact_cases.items():
        for field in ('s', 'time_s', 'theta_cmd', 'theta_m', 'x_s', 'x_n', 'error'):
            data[f'exact_{step_name}_{speed}_{field}'] = np.asarray(case[field])
    np.savez_compressed(data_path, **data)

    case_summary = {}
    for key, nonlinear in nonlinear_cases.items():
        baseline = baseline_cases[key]
        name = f"{key[0]}_{key[1]}"
        case_summary[name] = {
            "duration_s": float(nonlinear["time_s"][-1] + PLOT_DT),
            "peak_absolute_error_rev42_um": float(np.max(np.abs(nonlinear["error"])) * 1e6),
            "peak_absolute_error_baseline_um": float(np.max(np.abs(baseline["error"])) * 1e6),
            "final_error_rev42_nm": float(nonlinear["error"][-1] * 1e9),
            "final_error_baseline_nm": float(baseline["error"][-1] * 1e9),
            "minimum_instantaneous_interface_power_W": float(nonlinear["minimum_interface_power_W"]),
            "net_sampled_interface_work_J": float(nonlinear["interface_work_J"]),
            "solver_elapsed_s": float(nonlinear["elapsed_s"]),
            "solver_nfev": int(nonlinear["nfev"]),
        }
    for key, exact in exact_cases.items():
        name = f'{key[0]}_{key[1]}'
        smoothed = nonlinear_cases[key]
        case_summary[name].update({
            'peak_absolute_error_exact_um': float(np.max(np.abs(exact['error'])) * 1e6),
            'final_error_exact_nm': float(exact['error'][-1] * 1e9),
            'maximum_exact_vs_smoothed_difference_nm': float(
                np.max(np.abs(exact['error'] - smoothed['error'])) * 1e9
            ),
            'exact_solver_elapsed_s': float(exact['elapsed_s']),
            'exact_solver_nfev': int(exact['nfev']),
        })
    summary = {
        'exact_method': 'piecewise Radau, exact abs(v), numerical Jacobian',
        "method": "piecewise Radau with analytical Rev 4.2 Jacobian",
        "plot_output_dt_s": PLOT_DT,
        "diagnostic_output_dt_s": DIAGNOSTIC_DT,
        "solver_rtol": RTOL,
        "instantaneous_power_note": "F*v may be negative while a LuGre bristle returns stored energy; sign validation uses the exact -J.T*F virtual-power identity.",
        "cases": case_summary,
        "single_step_relative_convergence": convergence,
        "rev42_tangent_modes_hz": modes,
        "total_wall_time_s": time.perf_counter() - total_start,
        "figures": [
            montage_path.name,
            diagnostic_path.name,
            overlay_path.name,
            spectrum_path.name,
        ],
        "data": str(data_path.relative_to(ROOT)),
    }
    summary_path = ASSET_DIR / "stepping_rev42_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    print(f"Wrote {montage_path}")
    print(f"Wrote {diagnostic_path}")
    print(f"Wrote {overlay_path}")
    print(f"Wrote {spectrum_path}")
    print(f"Wrote {data_path}")


if __name__ == "__main__":
    main()
