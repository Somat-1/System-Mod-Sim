#!/usr/bin/env python3
"""15-state nonlinear LuGre-friction model for the Rev 4 ball-screw stage.

Sub-branch of Rev 4: builds on the frictionless recoupled baseline
(../model_parameters.json, ../build_bode_rev4.py) but does NOT modify or
import from it. This model:

  - drops k_nut, c_nut from K/C (the baseline's rigid, non-slipping nut
    coupling) and replaces the screw-nut interface with an explicit LuGre
    friction port, alongside LuGre ports at the support bearing and the
    guideway -- three interfaces total, matching Sec. 7 of
    state_space_6dof.md before the interfaces were assumed rigid/frozen.
  - keeps every other linear coupling from the baseline unchanged: the
    rotational chain (k_c/c_c, k_s1/c_s1, k_s2/c_s2, k_EM/k_d/c_EM on
    theta_m), and the axial bearing preload spring (k_brg/c_brg on x_s,
    which is NOT a friction interface -- it is the thrust-bearing pair's
    own elastic/damping behavior, distinct from T_fric,sb).

State vector (15,): [q(6), qdot(6), z_sb, z_nut, z_way]
  q    = [theta_m, theta_c, theta_s, theta_sb, x_s, x_n]
  qdot = [theta_m_dot, ..., x_n_dot]
  z_*  = LuGre internal bristle deflection at each of the 3 interfaces

Relative (sliding) velocities feeding each LuGre port:
  v_sb  = theta_sb_dot
  v_way = x_n_dot
  v_nut = (L/2pi)*theta_s_dot - (x_n_dot - x_s_dot)

IMPORTANT SIGN NOTE on v_nut and F_fric,nut injection: v_nut as defined
above is the *negative* of the "closing velocity" u_dot = x_n_dot - x_s_dot
- (L/2pi)*theta_s_dot that the baseline's corrected K/C sign convention
(state_space_6dof.md Sec. 3.4/5.3, 2026-08-18 fix) was built around --
that fix established x_n - x_s - (L/2pi)*theta_s as the physically correct
coupling variable (positive theta_s drives positive x_n). Because a LuGre
force is odd in its velocity argument, F_fric,nut computed from v_nut here
is the negative of what the baseline's F_nut meant, so the generalized
force injection signs on the x_s and x_n rows are flipped relative to
Sec. 4 of state_space_6dof.md to stay physically self-consistent. Derived
by virtual work: for v_nut = c_ts*theta_s_dot + c_xs*x_s_dot + c_xn*x_n_dot
(here c_ts=+L/2pi, c_xs=+1, c_xn=-1), the generalized force on q_i from an
interface force F opposing v_nut is Q_i = -F*c_i:
  theta_s row: -(L/2pi)*F_nut   (same sign as Sec. 4 -- coefficient itself
                                  did not flip, since both v_nut's sign and
                                  F_nut's meaning flipped, cancelling)
  x_s row:     -F_nut            (flipped vs. Sec. 4's "+F_nut")
  x_n row:     +F_nut            (flipped vs. Sec. 4's "-F_nut")
See lugre_friction/README.md for the full derivation and backlog.md for the
K/C sign fix this is built on top of.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PARAMETER_FILE = ROOT / "model_parameters.json"

STATE_LABELS = [
    "theta_m", "theta_c", "theta_s", "theta_sb", "x_s", "x_n",
    "theta_m_dot", "theta_c_dot", "theta_s_dot", "theta_sb_dot", "x_s_dot", "x_n_dot",
    "z_sb", "z_nut", "z_way",
]
N_STATES = 15


def load_parameters() -> dict[str, float]:
    payload = json.loads(PARAMETER_FILE.read_text(encoding="utf-8"))
    return payload["parameters"]


def build_linear_baseline(p: dict[str, float]):
    """M, K_lin, C_lin, B_em -- everything EXCEPT the three friction ports.

    K_lin/C_lin retain the rotational chain and the (non-friction) axial
    bearing preload spring k_brg/c_brg; they have NO k_nut/c_nut terms and
    NO cross-coupling between theta_s and the axial DOFs -- that coupling
    is now carried entirely by the nonlinear F_fric,nut term computed at
    each RHS evaluation.
    """
    k_EM = p["N_r"] * p["T_hold"]
    k_d = 4.0 * p["N_r"] * p["T_d"]

    M = np.diag([p["I_m"], p["I_c"], p["I_s"], p["I_sb"], p["M_screw"], p["M_s"]])

    k_c, k_s1, k_s2, k_brg = p["k_c"], p["k_s1"], p["k_s2"], p["k_brg"]
    K_lin = np.array([
        [k_c + k_EM + k_d, -k_c, 0.0, 0.0, 0.0, 0.0],
        [-k_c, k_c + k_s1, -k_s1, 0.0, 0.0, 0.0],
        [0.0, -k_s1, k_s1 + k_s2, -k_s2, 0.0, 0.0],
        [0.0, 0.0, -k_s2, k_s2, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, k_brg, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ])

    c_c, c_s1, c_s2, c_brg, c_EM = p["c_c"], p["c_s1"], p["c_s2"], p["c_brg"], p["c_EM"]
    C_lin = np.array([
        [c_c + c_EM, -c_c, 0.0, 0.0, 0.0, 0.0],
        [-c_c, c_c + c_s1, -c_s1, 0.0, 0.0, 0.0],
        [0.0, -c_s1, c_s1 + c_s2, -c_s2, 0.0, 0.0],
        [0.0, 0.0, -c_s2, c_s2, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, c_brg, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ])

    B_em = np.array([k_EM, 0.0, 0.0, 0.0, 0.0, 0.0])

    return M, K_lin, C_lin, B_em


def lugre_force(v: float, z: float, sigma0: float, sigma1: float, sigma2: float,
                 Fc: float, Fs: float, vs: float) -> tuple[float, float]:
    """One LuGre bristle: returns (friction force/torque, z_dot)."""
    g = Fc + (Fs - Fc) * np.exp(-(v / vs) ** 2)
    z_dot = v - sigma0 * abs(v) / g * z
    F = sigma0 * z + sigma1 * z_dot + sigma2 * v
    return F, z_dot


class LuGreModel:
    """Bundles the linear baseline + friction parameters; callable as an RHS
    for scipy.integrate.solve_ivp: model(t, state, theta_cmd_func)."""

    def __init__(self, p: dict[str, float] | None = None):
        self.p = p if p is not None else load_parameters()
        self.M, self.K_lin, self.C_lin, self.B_em = build_linear_baseline(self.p)
        self.M_inv = np.diag(1.0 / np.diag(self.M))
        self.lead_ratio = self.p["L"] / (2.0 * np.pi)

    def rhs(self, t: float, state: np.ndarray, theta_cmd_func) -> np.ndarray:
        p = self.p
        q = state[0:6]
        qdot = state[6:12]
        z_sb, z_nut, z_way = state[12], state[13], state[14]

        theta_sb_dot = qdot[3]
        x_s_dot = qdot[4]
        x_n_dot = qdot[5]
        theta_s_dot = qdot[2]

        v_sb = theta_sb_dot
        v_way = x_n_dot
        v_nut = self.lead_ratio * theta_s_dot - (x_n_dot - x_s_dot)

        T_fric_sb, z_sb_dot = lugre_force(
            v_sb, z_sb, p["sigma0_sb"], p["sigma1_sb"], p["sigma2_sb"],
            p["Tc_sb"], p["Ts_sb"], p["vs_sb"],
        )
        F_fric_way, z_way_dot = lugre_force(
            v_way, z_way, p["sigma0_way"], p["sigma1_way"], p["sigma2_way"],
            p["Fc_way"], p["Fs_way"], p["vs_way"],
        )
        F_fric_nut, z_nut_dot = lugre_force(
            v_nut, z_nut, p["sigma0_nut"], p["sigma1_nut"], p["sigma2_nut"],
            p["Fc_nut"], p["Fs_nut"], p["vs_nut"],
        )

        F_nl = np.zeros(6)
        F_nl[2] = -self.lead_ratio * F_fric_nut   # theta_s
        F_nl[3] = -T_fric_sb                       # theta_sb
        F_nl[4] = -F_fric_nut                       # x_s   (flipped vs. Sec. 4, see module docstring)
        F_nl[5] = F_fric_nut - F_fric_way           # x_n   (flipped vs. Sec. 4, see module docstring)

        theta_cmd = theta_cmd_func(t)
        accel = self.M_inv @ (F_nl + self.B_em * theta_cmd - self.K_lin @ q - self.C_lin @ qdot)

        d_state = np.empty(N_STATES)
        d_state[0:6] = qdot
        d_state[6:12] = accel
        d_state[12] = z_sb_dot
        d_state[13] = z_nut_dot
        d_state[14] = z_way_dot
        return d_state

    def __call__(self, t: float, state: np.ndarray, theta_cmd_func) -> np.ndarray:
        return self.rhs(t, state, theta_cmd_func)
