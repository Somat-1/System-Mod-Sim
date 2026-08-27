# Rev 4.2 — Jacobian-injected LuGre friction

Rev 4.2 retains the Rev 4 structural plant and adds three LuGre friction
ports through Jacobian-transposed generalized forces. Its main change from
Rev 4.1 is at the screw–nut interface: the load-bearing structural contact
`k_nut`/`c_nut` remains in the structural matrices, while the LuGre element
is an additional pre-rolling drag path in **parallel** across the same
relative coordinate. LuGre friction therefore does not replace the contact
stiffness or damping and cannot cap the transmitted thrust at its breakaway
force.

## 1. Coordinates and retained structural plant

Let

$$
\mathbf q=
\begin{bmatrix}
\theta_m&\theta_c&\theta_s&\theta_{sb}&x_s&x_n
\end{bmatrix}^{T},
\qquad
r=\frac{L}{2\pi},
$$

where the first four coordinates are rotations and the last two are axial
translations. The structural equation is

$$
\mathbf M\ddot{\mathbf q}
+\mathbf C\dot{\mathbf q}
+\mathbf K\mathbf q
=\mathbf G\theta_{cmd}
+\mathbf Q_{fric}
-\mathbf e_m T_d\sin(4N_r\theta_m),
$$

with

$$
\mathbf M=\operatorname{diag}
\left(I_m,I_c,I_s,I_{sb},M_{screw},M_s\right),
\qquad
\mathbf G=\begin{bmatrix}k_{EM}&0&0&0&0&0\end{bmatrix}^{T},
$$

and $\mathbf e_m=[1,0,0,0,0,0]^T$. The electromagnetic spring is already
included in $\mathbf K$ through
$T_{EM}=k_{EM}(\theta_{cmd}-\theta_m)$. The electromagnetic damping
$c_{EM}$ acts from the motor coordinate to ground. Detent is evaluated as
the nonlinear torque shown above; its tangent $k_d$ is not embedded in the
structural stiffness matrix.

The retained stiffness matrix is

$$
\mathbf K=
\begin{bmatrix}
k_c+k_{EM} &-k_c&0&0&0&0\\
-k_c&k_c+k_{s1}&-k_{s1}&0&0&0\\
0&-k_{s1}&k_{s1}+k_{s2}+r^2k_{nut}&-k_{s2}&rk_{nut}&-rk_{nut}\\
0&0&-k_{s2}&k_{s2}&0&0\\
0&0&rk_{nut}&0&k_{brg}+k_{nut}&-k_{nut}\\
0&0&-rk_{nut}&0&-k_{nut}&k_{nut}
\end{bmatrix}.
$$

The damping matrix has the same element topology:

$$
\mathbf C=
\begin{bmatrix}
c_c+c_{EM} &-c_c&0&0&0&0\\
-c_c&c_c+c_{s1}&-c_{s1}&0&0&0\\
0&-c_{s1}&c_{s1}+c_{s2}+r^2c_{nut}&-c_{s2}&rc_{nut}&-rc_{nut}\\
0&0&-c_{s2}&c_{s2}&0&0\\
0&0&rc_{nut}&0&c_{brg}+c_{nut}&-c_{nut}\\
0&0&-rc_{nut}&0&-c_{nut}&c_{nut}
\end{bmatrix}.
$$

Thus `k_nut` and `c_nut` are present before any LuGre force is added.

## 2. Friction ports and reaction signs

For port $p$, define its relative velocity by

$$
v_p=\mathbf J_p\dot{\mathbf q}.
$$

The port Jacobians used in Rev 4.2 are

| Port | $\mathbf J_p$ | Relative velocity |
|---|---|---|
| Guideway | $[0,0,0,0,0,1]$ | $v_{way}=\dot x_n$ |
| Screw–nut | $[0,0,r,0,1,-1]$ | $v_{nut}=r\dot\theta_s+\dot x_s-\dot x_n$ |
| Support bearing | $[0,0,0,1,0,0]$ | $v_{sb}=\dot\theta_{sb}$ |

The sign of the generalized reaction follows from virtual power. Here
$F_p$ is positive when it opposes positive $v_p$, so the power delivered to
the mechanical coordinates is

$$
\dot{\mathbf q}^{T}\mathbf Q_{fric,p}=-v_pF_p.
$$

Using $v_p=\mathbf J_p\dot{\mathbf q}$ gives

$$
\dot{\mathbf q}^{T}\mathbf Q_{fric,p}
=-\dot{\mathbf q}^{T}\mathbf J_p^{T}F_p
\quad\Longrightarrow\quad
\boxed{\mathbf Q_{fric,p}=-\mathbf J_p^{T}F_p}.
$$

Therefore

$$
\mathbf Q_{fric}
=-\mathbf J_{way}^{T}F_{way}
-\mathbf J_{nut}^{T}F_{nut}
-\mathbf J_{sb}^{T}T_{sb},
$$

or, expanded,

$$
\mathbf Q_{fric}=
\begin{bmatrix}
0\\0\\-rF_{nut}\\-T_{sb}\\-F_{nut}\\F_{nut}-F_{way}
\end{bmatrix}.
$$

This convention removes the need to assign reaction signs separately in
each mechanical equation. It also resolves the sign ambiguity in the source
notes: the expanded vector above is already $-\sum_p\mathbf J_p^TF_p$ and
must enter the structural equation as $+\mathbf Q_{fric}$.

## 3. Why the nut branches are parallel

Define the oriented screw–nut relative coordinate

$$
\rho_{nut}=\mathbf J_{nut}\mathbf q
=r\theta_s+x_s-x_n,
\qquad
\dot\rho_{nut}=v_{nut}.
$$

The opposite orientation,
$\delta_{nut}=x_n-x_s-r\theta_s=-\rho_{nut}$, describes the same physical
contact deformation. The structural nut potential and Rayleigh dissipation
are

$$
V_{nut}=\frac12 k_{nut}\rho_{nut}^{2},
\qquad
\mathcal R_{nut}=\frac12 c_{nut}\dot\rho_{nut}^{2}.
$$

They produce

$$
\mathbf Q_{nut,struct}
=-\mathbf J_{nut}^{T}
\left(k_{nut}\rho_{nut}+c_{nut}v_{nut}\right).
$$

The LuGre element sees exactly the same velocity $v_{nut}$ and produces

$$
\mathbf Q_{nut,LuGre}=-\mathbf J_{nut}^{T}F_{nut}.
$$

The complete screw–nut contribution is consequently

$$
\boxed{
\mathbf Q_{nut,total}
=-\mathbf J_{nut}^{T}
\left(k_{nut}\rho_{nut}+c_{nut}v_{nut}+F_{nut}\right)
}.
$$

All three constitutive terms experience the same relative motion and their
forces add. This is the analytical statement that the elastic contact,
viscous contact damping and LuGre pre-rolling drag are parallel branches.
In the implementation, the first two terms are assembled in $\mathbf K$ and
$\mathbf C$; the third is evaluated as nonlinear feedback.

## 4. LuGre law at each port

Each port has one bristle state $z_p$. To avoid the nondifferentiability of
$|v_p|$ at zero velocity, Rev 4.2 uses

$$
s_p=\sqrt{v_p^2+\varepsilon^2}.
$$

The Stribeck curve and decay rate are

$$
g_p(v_p)=F_{c,p}+\left(F_{s,p}-F_{c,p}\right)
\exp\!\left[-\left(\frac{v_p}{v_{s,p}}\right)^2\right],
\qquad
\lambda_p=\frac{\sigma_{0,p}s_p}{g_p(v_p)}.
$$

The bristle evolution and friction force are

$$
\boxed{
\dot z_p=v_p-\lambda_pz_p,
\qquad
F_p=\sigma_{0,p}z_p+\sigma_{1,p}\dot z_p+\sigma_{2,p}v_p
}.
$$

For the rotational support-bearing port, $F_p$ is interpreted as torque and
the corresponding parameters carry rotational units.

The computation sequence is therefore:

1. Evaluate the command and the current mechanical state.
2. Extract each port velocity with $v_p=\mathbf J_p\dot{\mathbf q}$.
3. Evaluate $g_p$, $\lambda_p$, $\dot z_p$ and $F_p$.
4. Map each scalar friction output back with
   $-\mathbf J_p^TF_p$.
5. Solve the mechanical acceleration using the retained structural
   $\mathbf M$, $\mathbf C$ and $\mathbf K$.

## 5. Nonlinear state-space model

With

$$
\mathbf v=\dot{\mathbf q},
\qquad
\mathbf z=
\begin{bmatrix}z_{way}&z_{nut}&z_{sb}\end{bmatrix}^{T},
\qquad
\mathbf x=
\begin{bmatrix}\mathbf q^T&\mathbf v^T&\mathbf z^T\end{bmatrix}^{T},
$$

the model has 15 states and can be written as

$$
\dot{\mathbf q}=\mathbf v,
$$

$$
\dot{\mathbf v}=\mathbf M^{-1}
\left[
\mathbf G\theta_{cmd}
-\mathbf C\mathbf v
-\mathbf K\mathbf q
-\mathbf e_mT_d\sin(4N_r\theta_m)
-\sum_p\mathbf J_p^TF_p
\right],
$$

$$
\dot z_p
=\mathbf J_p\mathbf v
-\frac{\sigma_{0,p}
\sqrt{(\mathbf J_p\mathbf v)^2+\varepsilon^2}}
{g_p(\mathbf J_p\mathbf v)}z_p.
$$

The code supplies the analytical Jacobian of this complete 15-state vector
field. The analytical result is checked against a complex-step numerical
Jacobian before the Bode response is accepted.

## 6. Bode comparison

`scripts/build_bode_rev42.py` linearizes the complete 15-state model at a
frozen 5 mm/s cruise point and compares $x_n/\theta_{cmd}$ with the Rev 4
frictionless baseline over 0–8 kHz. It saves the overlay, response arrays and
a JSON summary under `rendered_assets`.

The support-bearing bristle stiffness is
`sigma0_sb = 0.076 N m/rad`, not the earlier trial value of
`500 N m/rad`. Because $\partial F/\partial z$ contains $\sigma_0$, an
incorrectly large value behaves as a dominant tangent stiffness in the
linearized model and moves the associated resonance substantially.

## 7. Parameter caveats

The friction levels are reused from Rev 4.1 rather than newly identified.
The values of `sigma1` are initialized using

$$
\sigma_1=2\zeta\sqrt{\sigma_0m_{eff}},
\qquad \zeta=0.7,
$$

with the appropriate translational mass or rotational inertia at each port.
This is an initialization rule, not an identification result.

The baseline `c_nut = 101 N s/m` is retained even though LuGre `sigma1` and
`sigma2` add parallel dissipation. The split between structural contact
damping and frictional damping must be identified or repartitioned before
the combined damping is treated as calibrated. Retaining both here is
intentional for the topology comparison, but it can otherwise double-count
dissipation.

## 8. Nonlinear stepping

`scripts/generate_stepping_rev42.py` mirrors the Rev 4 full-step/16× and
fast/settled command cases. It uses piecewise Radau integration with the
analytical Jacobian, generates the montage, single-step diagnostic and axial
spectrum, and adds a four-case tracking-error overlay against the
frictionless baseline. The accepted trajectory's interface power and
cumulative work are saved with convergence and solver statistics in
`stepping_rev42_summary.json`.

Instantaneous $F_pv_p$ can become negative while a LuGre bristle returns
stored energy. The reaction-sign audit therefore checks the exact
$\dot{\mathbf q}^{T}\mathbf Q_{fric,p}=-v_pF_p$ mapping; it does not assume
that each port must be dissipative at every instant.
