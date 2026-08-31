# Guyan with nonlinear friction

This folder extends the existing two-master Guyan reduction with the Rev 4.2
parallel LuGre friction ports and the exact nonlinear detent torque. The retained
coordinates remain

\[
q_r=[\theta_m,\;x_n]^T.
\]

The reduced structural equation is

\[
M_r\ddot q_r+C_r\dot q_r+K_rq_r
=b_r\theta_{cmd}-\sum_p J_{p,r}^T F_p
-[T_d\sin(4N_r\theta_m),\;0]^T,
\]

where \(p\in\{way,nut,sb\}\). The original \(k_{nut}\) and \(c_{nut}\)
load-transmission path remains in \(K_r\) and \(C_r\); nut LuGre friction is
added in parallel. The linear detent stiffness is removed from \(K_r\), so the
nonlinear detent torque is not double counted.

The nonlinear state has seven entries:

\[
x=[q_r,\dot q_r,z_{way},z_{nut},z_{sb}]^T.
\]

## Bode definition

A nonlinear model has no unique amplitude-independent Bode response. The supplied
figure is therefore the analytical small-signal linearization at the zero-velocity,
zero-bristle rest equilibrium. The exact periodic detent contributes its local
slope \(4N_rT_d\) only inside that linearization; the time-domain right-hand side
continues to use the exact sine torque.

Run `python3 generate_bode.py` from this folder to regenerate the single Bode
figure and its numerical data.
