#!/usr/bin/env python3
"""Two-master Guyan model with parallel Rev 4.2 LuGre friction.

Reduced coordinates are q_r = [theta_m, x_n].  The baseline nut stiffness
and damping remain in K_r/C_r.  Way, nut and support-bearing friction act
through the projected full-model port Jacobians.  The detent is the exact
periodic torque T_d*sin(4*N_r*theta_m), not a linear spring in K_r.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
REV4 = PARENT.parent
LUGRE_PARAMETERS = REV4 / "lugre_friction" / "Rev 4.2" / "model_parameters.json"
PORTS = ("way", "nut", "sb")
N_Q = 2
N_STATES = 2 * N_Q + len(PORTS)


def load_parameters() -> dict[str, float]:
    return json.loads(LUGRE_PARAMETERS.read_text(encoding="utf-8"))["parameters"]


def build_reduced_structure(p: dict[str, float]) -> dict[str, np.ndarray]:
    """Project M/C/K with the linear detent removed before projection."""
    ell = p["L"] / (2.0 * np.pi)
    k_em = p["N_r"] * p["T_hold"]
    mass = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])
    stiffness = np.array([
        [p["k_c"] + k_em, -p["k_c"], 0, 0, 0, 0],
        [-p["k_c"], p["k_c"] + p["k_s1"], -p["k_s1"], 0, 0, 0],
        [0, -p["k_s1"], p["k_s1"] + p["k_s2"] + ell**2*p["k_nut"], -p["k_s2"], ell*p["k_nut"], -ell*p["k_nut"]],
        [0, 0, -p["k_s2"], p["k_s2"], 0, 0],
        [0, 0, ell*p["k_nut"], 0, p["k_brg"] + p["k_nut"], -p["k_nut"]],
        [0, 0, -ell*p["k_nut"], 0, -p["k_nut"], p["k_nut"]],
    ], dtype=float)
    damping = np.array([
        [p["c_c"] + p["c_EM"], -p["c_c"], 0, 0, 0, 0],
        [-p["c_c"], p["c_c"] + p["c_s1"], -p["c_s1"], 0, 0, 0],
        [0, -p["c_s1"], p["c_s1"] + p["c_s2"] + ell**2*p["c_nut"], -p["c_s2"], ell*p["c_nut"], -ell*p["c_nut"]],
        [0, 0, -p["c_s2"], p["c_s2"], 0, 0],
        [0, 0, ell*p["c_nut"], 0, p["c_brg"] + p["c_nut"], -p["c_nut"]],
        [0, 0, -ell*p["c_nut"], 0, -p["c_nut"], p["c_nut"]],
    ], dtype=float)
    command = np.array([k_em, 0, 0, 0, 0, 0], dtype=float)
    beta = p["k_nut"] / (p["k_nut"] + p["k_brg"])
    kappa = 1.0 / (1.0/p["k_nut"] + 1.0/p["k_brg"])
    k_ch = 1.0 / (1.0/p["k_c"] + 1.0/p["k_s1"])
    nu = ell**2*kappa / (k_ch + ell**2*kappa)
    mu = p["k_s1"] / (p["k_c"] + p["k_s1"])
    transform = np.array([
        [1, 0], [1-mu*nu, mu*nu/ell], [1-nu, nu/ell],
        [1-nu, nu/ell], [-ell*beta*(1-nu), beta*(1-nu)], [0, 1],
    ], dtype=float)
    full_ports = {
        "way": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        "nut": np.array([0.0, 0.0, ell, 0.0, 1.0, -1.0]),
        "sb": np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
    }
    result = {
        "T": transform,
        "M": transform.T @ mass @ transform,
        "C": transform.T @ damping @ transform,
        "K": transform.T @ stiffness @ transform,
        "b": transform.T @ command,
    }
    for port, jacobian in full_ports.items():
        result[f"J_{port}"] = jacobian @ transform
    return result


def _port_parameters(p: dict[str, float], port: str) -> tuple[float, ...]:
    force_prefix = "T" if port == "sb" else "F"
    return (
        p[f"sigma0_{port}"], p[f"sigma1_{port}"], p[f"sigma2_{port}"],
        p[f"{force_prefix}c_{port}"], p[f"{force_prefix}s_{port}"],
        p[f"vs_{port}"], p["smooth_velocity_epsilon"],
    )


def lugre_terms(v: float, z: float, values: tuple[float, ...]):
    sigma0, sigma1, sigma2, fc, fs, vs, epsilon = values
    exponential = np.exp(-(v / vs) ** 2)
    g = fc + (fs - fc) * exponential
    g_prime = (fs - fc) * exponential * (-2.0 * v / vs**2)
    speed = np.sqrt(v**2 + epsilon**2)
    speed_prime = v / speed
    decay = sigma0 * speed / g
    decay_prime = sigma0 * (speed_prime / g - speed * g_prime / g**2)
    z_dot = v - decay * z
    dzdot_dv = 1.0 - decay_prime * z
    dzdot_dz = -decay
    force = sigma0 * z + sigma1 * z_dot + sigma2 * v
    dforce_dv = sigma1 * dzdot_dv + sigma2
    dforce_dz = sigma0 + sigma1 * dzdot_dz
    return force, z_dot, dforce_dv, dforce_dz, dzdot_dv, dzdot_dz


class GuyanFrictionModel:
    def __init__(self, parameters: dict[str, float] | None = None):
        self.p = dict(parameters) if parameters is not None else load_parameters()
        reduced = build_reduced_structure(self.p)
        self.T = reduced["T"]
        self.mass = reduced["M"]
        self.damping = reduced["C"]
        self.stiffness = reduced["K"]
        self.command = reduced["b"]
        self.mass_inverse = np.linalg.inv(self.mass)
        self.jacobians = {port: reduced[f"J_{port}"] for port in PORTS}

    def rhs(self, _time: float, state: np.ndarray, theta_command: float) -> np.ndarray:
        q = state[:N_Q]
        velocity = state[N_Q:2 * N_Q]
        reaction = np.zeros(N_Q)
        z_dot = np.empty(len(PORTS))
        for index, port in enumerate(PORTS):
            jacobian = self.jacobians[port]
            force, zdot, *_ = lugre_terms(
                float(jacobian @ velocity), state[2 * N_Q + index],
                _port_parameters(self.p, port),
            )
            reaction -= jacobian * force
            z_dot[index] = zdot
        detent = np.array([
            self.p["T_d"] * np.sin(4.0 * self.p["N_r"] * q[0]), 0.0
        ])
        acceleration = self.mass_inverse @ (
            self.command * theta_command + reaction - self.damping @ velocity
            - self.stiffness @ q - detent
        )
        return np.concatenate((velocity, acceleration, z_dot))

    def analytical_linearization(self, state: np.ndarray | None = None):
        """Return A/B/C/D at a state; default is the zero rest equilibrium."""
        if state is None:
            state = np.zeros(N_STATES)
        velocity = state[N_Q:2 * N_Q]
        A = np.zeros((N_STATES, N_STATES))
        A[:N_Q, N_Q:2 * N_Q] = np.eye(N_Q)
        position_block = -self.stiffness.copy()
        position_block[0, 0] -= (
            4.0 * self.p["N_r"] * self.p["T_d"]
            * np.cos(4.0 * self.p["N_r"] * state[0])
        )
        velocity_block = -self.damping.copy()
        for index, port in enumerate(PORTS):
            jacobian = self.jacobians[port]
            terms = lugre_terms(
                float(jacobian @ velocity), state[2 * N_Q + index],
                _port_parameters(self.p, port),
            )
            _, _, dforce_dv, dforce_dz, dzdot_dv, dzdot_dz = terms
            velocity_block -= dforce_dv * np.outer(jacobian, jacobian)
            A[N_Q:2 * N_Q, 2 * N_Q + index] = (
                self.mass_inverse @ (-jacobian * dforce_dz)
            )
            A[2 * N_Q + index, N_Q:2 * N_Q] = dzdot_dv * jacobian
            A[2 * N_Q + index, 2 * N_Q + index] = dzdot_dz
        A[N_Q:2 * N_Q, :N_Q] = self.mass_inverse @ position_block
        A[N_Q:2 * N_Q, N_Q:2 * N_Q] = self.mass_inverse @ velocity_block
        B = np.zeros((N_STATES, 1))
        B[N_Q:2 * N_Q, 0] = self.mass_inverse @ self.command
        C = np.zeros((1, N_STATES)); C[0, 1] = 1.0
        D = np.zeros((1, 1))
        return A, B, C, D
