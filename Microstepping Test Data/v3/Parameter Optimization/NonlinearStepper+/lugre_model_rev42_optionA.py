#!/usr/bin/env python3
"""Rev 4.2 model with Option A: the electromagnetic drive replaced by the
current-projection torque law, instead of its small-angle linearization.

Reference: MathWorks Simscape Stepper Motor block electrical torque
  Te = -Km*(iA - eA/Rm)*sin(Nr*theta) + Km*(iB - eB/Rm)*cos(Nr*theta)
       - Td*sin(4*Nr*theta)
  (https://www.mathworks.com/help/sps/ref/steppermotor.html)

The -Td*sin(4*Nr*theta) detent term is already exactly what this repo's
detent implementation is (see lugre_model_rev42.py line 216) -- unchanged
here, and independently validated by that reference.

The electromagnetic part in that reference is a current projection, not a
spring: assuming an ideal chopper (iA, iB track the commanded current
vector exactly, eA/eB back-EMF neglected -- Option A, no new states), the
two-phase sum collapses to
  T_motor(theta_err) = T_hold * sin(N_r * theta_err),   theta_err = theta_cmd - theta_m
which saturates at +-T_hold when N_r*theta_err = +-pi/2 (theta_err = one
full step, pi/(2*N_r) rad = 1.8 degrees for N_r=50) and rolls back down
past that -- the actual pull-out/stall mechanism. This is exactly the
first-order term's parent nonlinearity: for small theta_err,
T_hold*sin(N_r*theta_err) ~ T_hold*N_r*theta_err = k_em*theta_err, i.e. it
reduces to the existing lugre_model_rev42.py linear drive for small lag,
by construction.

Nothing else changes: same 15 states, same mass/damping/K_c,K_s1,K_s2,
K_nut,K_brg structural chain, same three LuGre friction ports, same exact
sin() detent. Only the (0,0) stiffness entry (which held the linear k_em)
and the command vector (which also held k_em) are replaced by the
nonlinear torque evaluated directly in rhs()/analytical_linearization().

T_hold retains this repo's existing current-dependence convention
(T_hold(I) = sqrt(2)*K_t*I, K_t from model_parameters.json) -- deliberately
UNCHANGED from lugre_model_rev42.py for this comparison, to isolate the
one variable under test (linear vs. sin drive law) from the separate,
still-open question of which datasheet Km (0674A: 0.090 N*m/A, 0956A:
0.045 N*m/A) actually applies to this motor -- see README.md.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

REV42_SCRIPTS = Path(r"\\mult-fp01.hitdom.lan\project2\Internships\Tomas Valentinas\Sytem Mod & Sim\Rev 4\lugre_friction\Rev 4.2\scripts")
sys.path.insert(0, str(REV42_SCRIPTS))
from lugre_model_rev42 import (  # noqa: E402
    N_Q, N_STATES, PORTS, STATE_LABELS, load_parameters, port_jacobians,
    _port_values, lugre_terms, lugre_terms_exact,
)


def build_structural_matrices_optionA(p: dict[str, float]):
    """Same as lugre_model_rev42.build_structural_matrices, EXCEPT the
    (0,0) stiffness entry has no k_em (the electromagnetic coupling is now
    computed nonlinearly in rhs(), not folded into a constant spring) and
    there is no command vector (theta_cmd enters only through the
    nonlinear torque, not a linear b*theta_cmd forcing term)."""
    lead = p["L"] / (2.0 * np.pi)
    mass = np.diag([
        p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"],
    ])

    k_c, k_s1, k_s2 = p["k_c"], p["k_s1"], p["k_s2"]
    k_nut, k_brg = p["k_nut"], p["k_brg"]
    stiffness = np.array([
        [k_c, -k_c, 0.0, 0.0, 0.0, 0.0],          # <- no + k_em here (Option A)
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
    return mass, damping, stiffness


class LuGreModelRev42OptionA:
    def __init__(
        self,
        parameters: dict[str, float] | None = None,
        enforce_interface_power: bool = True,
        power_tolerance: float = 1.0e-12,
        regularization: str = 'smooth',
    ) -> None:
        self.p = dict(parameters) if parameters is not None else load_parameters()
        if regularization not in {'smooth', 'exact'}:
            raise ValueError('regularization must be smooth or exact')
        self.regularization = regularization
        self.mass, self.damping, self.stiffness = build_structural_matrices_optionA(self.p)
        self.mass_inverse = np.diag(1.0 / np.diag(self.mass))
        self.jacobians = port_jacobians(self.p)
        self.enforce_interface_power = enforce_interface_power
        self.power_tolerance = power_tolerance

    def motor_torque(self, theta_err: float) -> float:
        """T_hold*sin(N_r*theta_err) -- saturates at +-T_hold when
        theta_err reaches one full step (pi/(2*N_r) rad), then rolls off:
        the pull-out mechanism absent from the linear k_em*theta_err law."""
        return self.p["T_hold"] * np.sin(self.p["N_r"] * theta_err)

    def motor_torque_slope(self, theta_err: float) -> float:
        """d(motor_torque)/d(theta_err) = N_r*T_hold*cos(N_r*theta_err) --
        the local (operating-point-dependent) tangent stiffness, replacing
        the constant k_em used by the linear model's Jacobian."""
        return self.p["N_r"] * self.p["T_hold"] * np.cos(self.p["N_r"] * theta_err)

    def port_observables(self, state: np.ndarray) -> dict[str, dict[str, float]]:
        velocity = state[N_Q:2 * N_Q]
        observations: dict[str, dict[str, float]] = {}
        for index, port in enumerate(PORTS):
            v = self.jacobians[port] @ velocity
            z = state[2 * N_Q + index]
            term_function = (
                lugre_terms_exact if self.regularization == 'exact' else lugre_terms
            )
            terms = term_function(v, z, *_port_values(self.p, port))
            observations[port] = {
                "velocity": v, "z": z, "force": terms[0], "z_dot": terms[1],
                "dforce_dv": terms[2], "dforce_dz": terms[3],
                "dzdot_dv": terms[4], "dzdot_dz": terms[5],
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
        theta_err = theta_cmd - q[0]

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

        drive = np.zeros(N_Q, dtype=np.result_type(state, float))
        drive[0] = self.motor_torque(theta_err)  # Option A: nonlinear, not command*theta_cmd

        detent = np.zeros(N_Q, dtype=np.result_type(state, float))
        detent[0] = self.p["T_d"] * np.sin(4.0 * self.p["N_r"] * q[0])

        acceleration = self.mass_inverse @ (
            drive + friction_reaction - self.damping @ velocity
            - self.stiffness @ q - detent
        )

        derivative = np.empty(N_STATES, dtype=np.result_type(state, float))
        derivative[:N_Q] = velocity
        derivative[N_Q:2 * N_Q] = acceleration
        for index, port in enumerate(PORTS):
            derivative[2 * N_Q + index] = observations[port]["z_dot"]
        return derivative

    def analytical_linearization(self, state: np.ndarray, theta_cmd: float):
        """Same structure as lugre_model_rev42.LuGreModelRev42's, except
        the drive term's slope is the LOCAL tangent N_r*T_hold*cos(...) at
        the current (theta_cmd, theta_m), not a constant k_em -- so, unlike
        the linear model, this genuinely needs theta_cmd as an argument."""
        observations = self.port_observables(state)
        theta_err = theta_cmd - state[0]
        drive_slope = self.motor_torque_slope(theta_err)

        system = np.zeros((N_STATES, N_STATES))
        system[:N_Q, N_Q:2 * N_Q] = np.eye(N_Q)

        position_block = -self.stiffness.copy()
        position_block[0, 0] -= drive_slope  # d(drive)/d(theta_m) = -drive_slope
        detent_slope = (
            4.0 * self.p["N_r"] * self.p["T_d"] * np.cos(4.0 * self.p["N_r"] * state[0])
        )
        position_block[0, 0] -= detent_slope
        velocity_block = -self.damping.copy()

        for index, port in enumerate(PORTS):
            jacobian = self.jacobians[port]
            observation = observations[port]
            velocity_block -= observation["dforce_dv"] * np.outer(jacobian, jacobian)
            system[N_Q:2 * N_Q, 2 * N_Q + index] = (
                self.mass_inverse @ (-jacobian * observation["dforce_dz"])
            )
            system[2 * N_Q + index, N_Q:2 * N_Q] = observation["dzdot_dv"] * jacobian
            system[2 * N_Q + index, 2 * N_Q + index] = observation["dzdot_dz"]

        system[N_Q:2 * N_Q, :N_Q] = self.mass_inverse @ position_block
        system[N_Q:2 * N_Q, N_Q:2 * N_Q] = self.mass_inverse @ velocity_block
        input_vector = np.zeros(N_STATES)
        input_vector[N_Q] = drive_slope / self.p["I_m"]  # d(drive)/d(theta_cmd) = +drive_slope
        output = np.zeros((1, N_STATES))
        output[0, 5] = 1.0
        return system, input_vector, output
