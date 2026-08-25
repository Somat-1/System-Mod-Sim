#!/usr/bin/env python3
"""Rev 4.2 15-state LuGre model with Jacobian-defined friction ports.

State: x = [q(6), qdot(6), z_way, z_nut, z_sb].  Each port is defined by
v_p = J_p qdot and contributes -J_p.T F_p to the generalized force.  The
baseline k_nut/c_nut path remains in K/C; the nut LuGre force is parallel
pre-rolling drag, not the load-transmission element.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
PARAMETER_FILE = ROOT / "model_parameters.json"
N_Q = 6
PORTS = ("way", "nut", "sb")
N_STATES = 2 * N_Q + len(PORTS)
STATE_LABELS = (
    "theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n",
    "theta_m_dot", "theta_c_dot", "theta_s_dot", "theta_sb_dot",
    "x_s_dot", "x_n_dot", "z_way", "z_nut", "z_sb",
)


def load_parameters() -> dict[str, float]:
    return json.loads(PARAMETER_FILE.read_text(encoding="utf-8"))["parameters"]


def build_structural_matrices(
    p: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build M/C/K/B with k_nut/c_nut intact and nonlinear detent omitted.

    K includes electromagnetic stiffness but not k_d.  B contains only
    k_EM, as required for the single theta_cmd input in the nonlinear model.
    """
    k_em = p["N_r"] * p["T_hold"]
    lead = p["L"] / (2.0 * np.pi)
    mass = np.diag([
        p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"],
    ])

    k_c, k_s1, k_s2 = p["k_c"], p["k_s1"], p["k_s2"]
    k_nut, k_brg = p["k_nut"], p["k_brg"]
    stiffness = np.array([
        [k_c + k_em, -k_c, 0.0, 0.0, 0.0, 0.0],
        [-k_c, k_c + k_s1, -k_s1, 0.0, 0.0, 0.0],
        [0.0, -k_s1, k_s1 + k_s2 + lead**2 * k_nut, -k_s2,
         lead * k_nut, -lead * k_nut],
        [0.0, 0.0, -k_s2, k_s2, 0.0, 0.0],
        [0.0, 0.0, lead * k_nut, 0.0, k_brg + k_nut, -k_nut],
        [0.0, 0.0, -lead * k_nut, 0.0, -k_nut, k_nut],
    ])

    c_c, c_s1, c_s2 = p["c_c"], p["c_s1"], p["c_s2"]
    c_nut, c_brg, c_em = p["c_nut"], p["c_brg"], p["c_EM"]
    damping = np.array([
        [c_c + c_em, -c_c, 0.0, 0.0, 0.0, 0.0],
        [-c_c, c_c + c_s1, -c_s1, 0.0, 0.0, 0.0],
        [0.0, -c_s1, c_s1 + c_s2 + lead**2 * c_nut, -c_s2,
         lead * c_nut, -lead * c_nut],
        [0.0, 0.0, -c_s2, c_s2, 0.0, 0.0],
        [0.0, 0.0, lead * c_nut, 0.0, c_brg + c_nut, -c_nut],
        [0.0, 0.0, -lead * c_nut, 0.0, -c_nut, c_nut],
    ])
    command = np.array([k_em, 0.0, 0.0, 0.0, 0.0, 0.0])
    return mass, damping, stiffness, command


def port_jacobians(p: dict[str, float]) -> dict[str, np.ndarray]:
    lead = p["L"] / (2.0 * np.pi)
    return {
        "way": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        "nut": np.array([0.0, 0.0, lead, 0.0, 1.0, -1.0]),
        "sb": np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
    }


def _port_values(p: dict[str, float], port: str) -> tuple[float, ...]:
    force_prefix = "T" if port == "sb" else "F"
    return (
        p[f"sigma0_{port}"], p[f"sigma1_{port}"], p[f"sigma2_{port}"],
        p[f"{force_prefix}c_{port}"], p[f"{force_prefix}s_{port}"],
        p[f"vs_{port}"], p["smooth_velocity_epsilon"],
    )


def lugre_terms(
    v: float,
    z: float,
    sigma0: float,
    sigma1: float,
    sigma2: float,
    fc: float,
    fs: float,
    vs: float,
    epsilon: float,
) -> tuple[float, float, float, float, float, float]:
    """Return F, zdot and their analytical derivatives with respect to v,z."""
    exp_term = np.exp(-(v / vs) ** 2)
    g = fc + (fs - fc) * exp_term
    g_prime = (fs - fc) * exp_term * (-2.0 * v / vs**2)
    smooth_speed = np.sqrt(v**2 + epsilon**2)
    smooth_prime = v / smooth_speed
    decay = sigma0 * smooth_speed / g
    decay_prime = sigma0 * (
        smooth_prime / g - smooth_speed * g_prime / g**2
    )
    z_dot = v - decay * z
    dzdot_dv = 1.0 - decay_prime * z
    dzdot_dz = -decay
    force = sigma0 * z + sigma1 * z_dot + sigma2 * v
    dforce_dv = sigma1 * dzdot_dv + sigma2
    dforce_dz = sigma0 + sigma1 * dzdot_dz
    return force, z_dot, dforce_dv, dforce_dz, dzdot_dv, dzdot_dz


class LuGreModelRev42:
    def __init__(
        self,
        parameters: dict[str, float] | None = None,
        enforce_interface_power: bool = True,
        power_tolerance: float = 1.0e-12,
    ) -> None:
        self.p = dict(parameters) if parameters is not None else load_parameters()
        self.mass, self.damping, self.stiffness, self.command = (
            build_structural_matrices(self.p)
        )
        self.mass_inverse = np.diag(1.0 / np.diag(self.mass))
        self.jacobians = port_jacobians(self.p)
        self.enforce_interface_power = enforce_interface_power
        self.power_tolerance = power_tolerance

    def port_observables(
        self, state: np.ndarray
    ) -> dict[str, dict[str, float]]:
        velocity = state[N_Q:2 * N_Q]
        observations: dict[str, dict[str, float]] = {}
        for index, port in enumerate(PORTS):
            v = self.jacobians[port] @ velocity
            z = state[2 * N_Q + index]
            terms = lugre_terms(v, z, *_port_values(self.p, port))
            observations[port] = {
                "velocity": v,
                "z": z,
                "force": terms[0],
                "z_dot": terms[1],
                "dforce_dv": terms[2],
                "dforce_dz": terms[3],
                "dzdot_dv": terms[4],
                "dzdot_dz": terms[5],
            }
        return observations

    def rhs(
        self,
        _time: float,
        state: np.ndarray,
        theta_command: float | Callable[[float], float],
    ) -> np.ndarray:
        q = state[:N_Q]
        velocity = state[N_Q:2 * N_Q]
        observations = self.port_observables(state)
        theta_cmd = theta_command(_time) if callable(theta_command) else theta_command

        friction_reaction = np.zeros(N_Q, dtype=np.result_type(state, float))
        interface_power = 0.0
        for port in PORTS:
            observation = observations[port]
            friction_reaction -= self.jacobians[port] * observation["force"]
            interface_power += observation["force"] * observation["velocity"]

        if (
            self.enforce_interface_power
            and np.isrealobj(state)
            and interface_power < -self.power_tolerance
        ):
            raise AssertionError(
                f"Negative total interface power {interface_power:.6e} W"
            )

        detent = np.zeros(N_Q, dtype=np.result_type(state, float))
        detent[0] = self.p["T_d"] * np.sin(4.0 * self.p["N_r"] * q[0])
        acceleration = self.mass_inverse @ (
            self.command * theta_cmd
            + friction_reaction
            - self.damping @ velocity
            - self.stiffness @ q
            - detent
        )

        derivative = np.empty(N_STATES, dtype=np.result_type(state, float))
        derivative[:N_Q] = velocity
        derivative[N_Q:2 * N_Q] = acceleration
        for index, port in enumerate(PORTS):
            derivative[2 * N_Q + index] = observations[port]["z_dot"]
        return derivative

    def analytical_linearization(
        self, state: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the exact tangent A/B/C matrices at a real 15-state point."""
        observations = self.port_observables(state)
        system = np.zeros((N_STATES, N_STATES))
        system[:N_Q, N_Q:2 * N_Q] = np.eye(N_Q)

        position_block = -self.stiffness.copy()
        detent_slope = (
            4.0 * self.p["N_r"] * self.p["T_d"]
            * np.cos(4.0 * self.p["N_r"] * state[0])
        )
        position_block[0, 0] -= detent_slope
        velocity_block = -self.damping.copy()

        for index, port in enumerate(PORTS):
            jacobian = self.jacobians[port]
            observation = observations[port]
            velocity_block -= observation["dforce_dv"] * np.outer(jacobian, jacobian)
            system[N_Q:2 * N_Q, 2 * N_Q + index] = (
                self.mass_inverse
                @ (-jacobian * observation["dforce_dz"])
            )
            system[2 * N_Q + index, N_Q:2 * N_Q] = (
                observation["dzdot_dv"] * jacobian
            )
            system[2 * N_Q + index, 2 * N_Q + index] = observation["dzdot_dz"]

        system[N_Q:2 * N_Q, :N_Q] = self.mass_inverse @ position_block
        system[N_Q:2 * N_Q, N_Q:2 * N_Q] = self.mass_inverse @ velocity_block
        input_vector = np.zeros(N_STATES)
        input_vector[N_Q:2 * N_Q] = self.mass_inverse @ self.command
        output = np.zeros((1, N_STATES))
        output[0, 5] = 1.0
        return system, input_vector, output

    def cruise_state(self, stage_velocity: float) -> np.ndarray:
        """Frozen tangent point with ideal screw tracking at stage_velocity."""
        state = np.zeros(N_STATES)
        omega = stage_velocity * 2.0 * np.pi / self.p["L"]
        state[N_Q:2 * N_Q] = [omega, omega, omega, omega, 0.0, stage_velocity]
        for index, port in enumerate(PORTS):
            v = self.jacobians[port] @ state[N_Q:2 * N_Q]
            sigma0, _, _, fc, fs, vs, epsilon = _port_values(self.p, port)
            g = fc + (fs - fc) * np.exp(-(v / vs) ** 2)
            state[2 * N_Q + index] = (
                v * g / (sigma0 * np.sqrt(v**2 + epsilon**2))
                if sigma0 > 0.0 else 0.0
            )
        return state


if __name__ == "__main__":
    model = LuGreModelRev42(enforce_interface_power=True)
    operating_state = model.cruise_state(5.0e-3)
    model.rhs(0.0, operating_state, 0.0)
    print("Rev 4.2 model constructed; interface-power assertion passed at cruise.")
