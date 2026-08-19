#!/usr/bin/env python3
"""Numerical Jacobian and eigenvalue check for the nonlinear 15-state LuGre
RHS (lugre_model.LuGreModel.rhs).

J = d f(x, u) / dx, computed by central finite differences of
model.rhs(t, x, theta_cmd_func) with theta_cmd held constant, evaluated at
rest: x = 0 (every state, including all three bristle deflections) and
theta_cmd = 0.

Rest is also exactly the LuGre non-smooth kink: v = 0 at every one of the
three ports simultaneously, where the |v|/g(v) term in z_dot = v -
sigma0*|v|/g(v)*z is not differentiable (its one-sided derivatives from
v>0 and v<0 disagree). backlog.md and lugre_friction/README.md flag this
kink as what defeats Radau/BDF's Newton-iteration/finite-difference-
Jacobian machinery (both diverge to NaN on this model; LSODA, which
tolerates it via explicit/implicit switching, is used everywhere else in
this sub-branch instead). This script computes that exact Jacobian
directly, at the exact point where it's ill-defined, rather than only
inferring the problem from solver failure logs: central differencing
through a kink doesn't raise an error, it just silently averages the two
disagreeing one-sided slopes into a number that isn't the true derivative
in any direction -- the eigenvalues below should be read with that caveat,
not taken as a clean small-signal model at this operating point.

Per-state finite-difference step sizes are NOT one scalar: positions
(theta_*, x_s, x_n), velocities (theta_*_dot, x_s_dot, x_n_dot), and
bristle deflections (z_sb rotational; z_nut, z_way translational) span
very different physical scales, and the two translational bristles are
already ~1e-7 m at their own breakaway limit (Fs/sigma0) -- a step sized
for a radian would either swamp them or, sized for them, be numerical
noise for a radian. Rotational/translational pairs are scaled by
lead_ratio = L/(2*pi), the same ratio used for the atol vector in the
sinusoidal sweep and flagged in state_space_6dof.md Sec. 9 item 3.
"""

from __future__ import annotations

import numpy as np

from lugre_model import LuGreModel, N_STATES, STATE_LABELS

# Base step sizes for a rotational-position, rotational-velocity, and
# rotational-bristle state; translational counterparts are this times
# lead_ratio (~1.59e-4 for L=1.0e-3 m), so the perturbation represents an
# "equally sized" nudge in the coupled rotational/translational sense.
STEP_THETA = 1.0e-6      # rad
STEP_THETA_DOT = 1.0e-4  # rad/s
STEP_Z_ROT = 1.0e-8      # rad (z_sb)


def build_step_vector(lead_ratio: float) -> np.ndarray:
    step_x = STEP_THETA * lead_ratio
    step_x_dot = STEP_THETA_DOT * lead_ratio
    step_z_trans = STEP_Z_ROT * lead_ratio
    return np.array([
        STEP_THETA, STEP_THETA, STEP_THETA, STEP_THETA, step_x, step_x,               # q
        STEP_THETA_DOT, STEP_THETA_DOT, STEP_THETA_DOT, STEP_THETA_DOT, step_x_dot, step_x_dot,  # qdot
        STEP_Z_ROT, step_z_trans, step_z_trans,                                       # z_sb, z_nut, z_way
    ])


def numerical_jacobian(model: LuGreModel, x0: np.ndarray, theta_cmd: float,
                        step: np.ndarray) -> np.ndarray:
    def theta_cmd_func(_t):
        return theta_cmd

    def f(x):
        return model.rhs(0.0, x, theta_cmd_func)

    n = len(x0)
    J = np.zeros((n, n))
    for j in range(n):
        dx = np.zeros(n)
        dx[j] = step[j]
        f_plus = f(x0 + dx)
        f_minus = f(x0 - dx)
        J[:, j] = (f_plus - f_minus) / (2.0 * step[j])
    return J


def main() -> None:
    model = LuGreModel()
    x0 = np.zeros(N_STATES)
    theta_cmd0 = 0.0
    step = build_step_vector(model.lead_ratio)

    print("=" * 78)
    print("Nonlinear LuGre model: numerical Jacobian + eigenvalues at rest")
    print("=" * 78)
    print(f"States ({N_STATES}): {STATE_LABELS}")
    print(f"Operating point: x = 0 (all states), theta_cmd = {theta_cmd0}")
    print("  -> v_sb = v_way = v_nut = 0 at this point: the LuGre non-smooth")
    print("     kink. See this script's module docstring before trusting the")
    print("     velocity-state columns of J below at face value.")
    print(f"lead_ratio = L/(2*pi) = {model.lead_ratio:.6e}")
    print()

    J = numerical_jacobian(model, x0, theta_cmd0, step)

    np.set_printoptions(precision=4, suppress=False, linewidth=200)
    print("Jacobian J = df/dx (15x15):")
    print(J)
    print()

    eigvals = np.linalg.eigvals(J)
    order = np.argsort(-eigvals.real)  # least stable (largest real part) first

    print("Eigenvalues of J, sorted by real part descending (least stable first):")
    print(f"  {'Re(lambda)':>14s} {'Im(lambda)':>14s} {'|lambda|':>12s} {'f (Hz)':>10s} {'zeta':>9s}")
    for k in order:
        lam = eigvals[k]
        mag = abs(lam)
        if abs(lam.imag) > 1e-9:
            freq_hz = mag / (2.0 * np.pi)
            zeta = -lam.real / mag if mag > 0 else float("nan")
            print(f"  {lam.real:+14.4e} {lam.imag:+14.4e} {mag:12.4e} {freq_hz:10.2f} {zeta:+9.4f}")
        else:
            print(f"  {lam.real:+14.4e} {0.0:+14.4e} {mag:12.4e} {'--':>10s} {'--':>9s}")

    n_unstable = int(np.sum(eigvals.real > 1e-9))
    n_marginal = int(np.sum(np.abs(eigvals.real) <= 1e-9))
    print()
    print(f"{n_unstable}/{N_STATES} eigenvalues with positive real part (unstable),")
    print(f"{n_marginal}/{N_STATES} with |Re(lambda)| <= 1e-9 (marginal/rigid-body-like).")


if __name__ == "__main__":
    main()
