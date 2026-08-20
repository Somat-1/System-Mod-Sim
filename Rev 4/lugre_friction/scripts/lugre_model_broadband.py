#!/usr/bin/env python3
"""Nonlinear-drive variant of LuGreModel: the electromagnetic and detent
torques on theta_m are kept in their true nonlinear (sin) form instead of
the small-signal k_EM/k_d linearization the rest of the sub-branch uses.

state_space_6dof.md Sec 3.1/3.2 give the pre-linearization forms:
    T_EM  = T_hold * sin(N_r * (theta_cmd - theta_m))
    T_det = T_d    * sin(4*N_r * theta_m)
Both linearizations (k_EM = N_r*T_hold, k_d = 4*N_r*T_d) are first-order
Taylor terms of these sin()s at zero argument, valid only for
theta_m << 1/(4*N_r) ~= 5 mrad (backlog.md item 5) and, by the same
reasoning, |theta_cmd - theta_m| << 1/N_r = 20 mrad for T_EM. The broadband
identification amplitudes in run_broadband_id.py are tens of mrad RMS,
comparable to or past both thresholds, so neither linearization holds and
both must be replaced by their nonlinear form here -- the same reasoning
already applied to the three LuGre friction ports elsewhere in this
sub-branch.

Everything else (friction ports, rotational chain k_c/k_s1/k_s2, axial
k_brg/c_brg, viscous c_EM) is unchanged from lugre_model.LuGreModel; this
subclass only overrides K_lin[0,0]/B_em (drops k_EM, k_d) and rhs() (adds
the nonlinear T_EM/T_det torques directly instead of the linear terms they
replace).
"""

from __future__ import annotations

import numpy as np

from lugre_model import N_STATES, LuGreModel, lugre_force


class LuGreModelNonlinearDrive(LuGreModel):
    def __init__(self, p: dict[str, float] | None = None):
        super().__init__(p)
        k_EM = self.p["N_r"] * self.p["T_hold"]
        k_d = 4.0 * self.p["N_r"] * self.p["T_d"]
        # Undo the small-signal embedding: K_lin[0,0] was k_c + k_EM + k_d
        # (build_linear_baseline); the nonlinear T_EM/T_det below replace
        # both, leaving only k_c (theta_m<->theta_c coupling). c_EM is
        # viscous damping against ground, not part of either sin() torque,
        # and stays in C_lin[0,0] unchanged.
        self.K_lin = self.K_lin.copy()
        self.K_lin[0, 0] -= (k_EM + k_d)
        # B_em fed k_EM*theta_cmd as a linear feedforward; T_EM below
        # carries theta_cmd directly, so B_em no longer contributes.
        self.B_em = np.zeros_like(self.B_em)

    def rhs(self, t: float, state: np.ndarray, theta_cmd_func) -> np.ndarray:
        p = self.p
        q = state[0:6]
        qdot = state[6:12]
        z_sb, z_nut, z_way = state[12], state[13], state[14]

        theta_m = q[0]
        theta_s_dot, theta_sb_dot, x_s_dot, x_n_dot = qdot[2], qdot[3], qdot[4], qdot[5]

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

        theta_cmd = theta_cmd_func(t)
        T_EM = p["T_hold"] * np.sin(p["N_r"] * (theta_cmd - theta_m))
        T_det = p["T_d"] * np.sin(4.0 * p["N_r"] * theta_m)

        F_nl = np.zeros(6)
        F_nl[0] = T_EM - T_det                       # theta_m -- Sec 4 Eq (1) sign convention
        F_nl[2] = -self.lead_ratio * F_fric_nut       # theta_s
        F_nl[3] = -T_fric_sb                          # theta_sb
        F_nl[4] = -F_fric_nut                         # x_s
        F_nl[5] = F_fric_nut - F_fric_way             # x_n

        accel = self.M_inv @ (F_nl - self.K_lin @ q - self.C_lin @ qdot)

        d_state = np.empty(N_STATES)
        d_state[0:6] = qdot
        d_state[6:12] = accel
        d_state[12] = z_sb_dot
        d_state[13] = z_nut_dot
        d_state[14] = z_way_dot
        return d_state

    def __call__(self, t: float, state: np.ndarray, theta_cmd_func) -> np.ndarray:
        return self.rhs(t, state, theta_cmd_func)
