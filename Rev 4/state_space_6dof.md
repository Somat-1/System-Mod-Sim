# Simplified 6-Coordinate Drivetrain Model: State-Space Derivation

Reduced lumped model of the open-loop stepper-driven positioning axis.
Retained degrees of freedom: motor rotation, coupling rotation, screw rotation, support-bearing rotation, screw axial translation, nut/stage axial translation.

---

## 1. Conventions

Positive rotation `θ+` about the screw axis drives positive translation `x+`.
The lead `L` is the axial advance per screw revolution, so the ideal kinematic ratio is

$$
x = \frac{L}{2\pi}\,\theta , \qquad
\frac{\partial \theta}{\partial x} = \frac{2\pi}{L}
$$

All stiffness and damping elements are linear and act on relative coordinates.
All friction elements are treated as **ports**. Their outputs enter the model as generalized forces, not as constitutive matrix entries.

---

## 2. Generalized coordinates

$$
\mathbf{q} =
\begin{bmatrix}
\theta_m & \theta_c & \theta_s & \theta_{sb} & x_s & x_n
\end{bmatrix}^{\mathsf{T}}
$$

| Index | Symbol | Body | Inertia | Unit |
|---|---|---|---|---|
| 1 | $\theta_m$ | Motor rotor | $I_m$ | rad |
| 2 | $\theta_c$ | Bellows coupling | $I_c$ | rad |
| 3 | $\theta_s$ | Screw shaft (rotation) | $I_s$ | rad |
| 4 | $\theta_{sb}$ | Support bearing inner ring | $I_{sb}$ | rad |
| 5 | $x_s$ | Screw shaft (axial) | $M_{screw}$ | m |
| 6 | $x_n$ | Nut and stage assembly | $M_s$ | m |

> Note on the sketch: the hand derivation labels the screw axial mass $m_s$ and the stage mass $M_s$. In this document the screw axial mass is written $M_{screw}$ to remove the ambiguity with $M_s$.

---

## 3. Element constitutive relations

### 3.1 Electromagnetic torque (microstep detent well)

$$
T_{EM} = k_{EM}\,(\theta_{cmd} - \theta_m)
\;\big[\; + \; c_{EM}\,(\dot{\theta}_{cmd} - \dot{\theta}_m)\;\big]
$$

Linearized about the commanded equilibrium, with $N_r$ rotor teeth and holding torque $T_{hold}$:

$$
k_{EM} = N_r\,T_{hold}
$$

### 3.2 Detent (cogging) torque

$$
T_{det} = T_d \sin(4 N_r \theta_m)
\;\xrightarrow[\;\sin\phi \approx \phi\;]{}\;
k_d\,\theta_m ,
\qquad
k_d = 4 N_r T_d
$$

The factor $4N_r$ gives 200 detent cycles per revolution for a 1.8 deg hybrid motor ($N_r = 50$), one period per full step.

### 3.3 Torsional couplings

$$
\begin{aligned}
T_{c} &= k_{c}\,(\theta_m - \theta_c) \;\big[\; + \; c_{c}\,(\dot{\theta}_m - \dot{\theta}_c)\;\big] \\
T_{s1} &= k_{s1}\,(\theta_c - \theta_s) \;\big[\; + \; c_{s1}\,(\dot{\theta}_c - \dot{\theta}_s)\;\big] \\
T_{s2} &= k_{s2}\,(\theta_s - \theta_{sb}) \;\big[\; + \; c_{s2}\,(\dot{\theta}_s - \dot{\theta}_{sb})\;\big]
\end{aligned}
$$

Each appears once as a drive torque on the downstream body and once as a reaction torque on the upstream body.

### 3.4 Axial couplings

$$
\begin{aligned}
F_{nut} &= k_{nut}\,\Big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\Big) \;\big[\; + \; c_{nut}\,\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\;\big] \\
F_{ax,\text{react}} &= k_{brg}\,x_s \;\big[\; + \; c_{brg}\,\dot{x}_s\;\big]
\end{aligned}
$$

$k_{nut}$ is the screw-nut axial contact stiffness. As of this revision it is coupled directly to screw rotation through the lead ratio $L/2\pi$, rather than injected as a separate friction port (Sec. 7, Sec. 10 items 1-3). The sign inside $F_{nut}$ is fixed by the Sec. 1 convention: $x_n - x_s - \tfrac{L}{2\pi}\theta_s$ is the nut's actual position minus the ideal no-slip position $x_s + \tfrac{L}{2\pi}\theta_s$, so positive $\theta_s$ drives positive $x_n$ (2026-08-18 correction; an earlier draft of this revision had the opposite sign on the $\theta_s$ term, which inverted the sign of the command-to-stage DC gain). $k_{brg}$ is the axial rigidity of the support bearing pair, grounded.

---

## 4. Equations of motion

### 4.1 Compact form

$$
\begin{aligned}
\text{(1) Motor rotor:} && T_{EM} - T_{det} - T_{c,\text{react}} &= I_m \ddot{\theta}_m \\
\text{(2) Coupling:} && T_{c,\text{drive}} - T_{s1,\text{react}} &= I_c \ddot{\theta}_c \\
\text{(3) Screw rotation:} && T_{s1,\text{drive}} - \tfrac{L}{2\pi}F_{nut} - T_{s2,\text{react}} &= I_s \ddot{\theta}_s \\
\text{(4) Support bearing:} && T_{s2,\text{drive}} - T_{fric,sb} &= I_{sb} \ddot{\theta}_{sb} \\
\text{(5) Screw axial:} && F_{nut} - F_{ax,\text{react}} &= M_{screw} \ddot{x}_s \\
\text{(6) Nut and stage:} && -F_{nut} - F_{fric,way} &= M_s \ddot{x}_n
\end{aligned}
$$

### 4.2 Fully expanded

$$
\begin{aligned}
\text{(1)}\quad
& k_{EM}(\theta_{cmd} - \theta_m) - T_d \sin(4N_r\theta_m)
- \Big[k_c(\theta_m - \theta_c) + c_c(\dot{\theta}_m - \dot{\theta}_c)\Big]
= I_m \ddot{\theta}_m \\[1.5ex]
\text{(2)}\quad
& \Big[k_c(\theta_m - \theta_c) + c_c(\dot{\theta}_m - \dot{\theta}_c)\Big]
- \Big[k_{s1}(\theta_c - \theta_s) + c_{s1}(\dot{\theta}_c - \dot{\theta}_s)\Big]
= I_c \ddot{\theta}_c \\[1.5ex]
\text{(3)}\quad
& \Big[k_{s1}(\theta_c - \theta_s) + c_{s1}(\dot{\theta}_c - \dot{\theta}_s)\Big]
- \frac{L}{2\pi}\Big[k_{nut}\big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\big) + c_{nut}\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\Big]
- \Big[k_{s2}(\theta_s - \theta_{sb}) + c_{s2}(\dot{\theta}_s - \dot{\theta}_{sb})\Big]
= I_s \ddot{\theta}_s \\[1.5ex]
\text{(4)}\quad
& \Big[k_{s2}(\theta_s - \theta_{sb}) + c_{s2}(\dot{\theta}_s - \dot{\theta}_{sb})\Big]
- T_{fric,sb}
= I_{sb} \ddot{\theta}_{sb} \\[1.5ex]
\text{(5)}\quad
& \Big[k_{nut}\big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\big) + c_{nut}\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\Big]
- \Big[k_{brg}x_s + c_{brg}\dot{x}_s\Big]
= M_{screw} \ddot{x}_s \\[1.5ex]
\text{(6)}\quad
& -\Big[k_{nut}\big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\big) + c_{nut}\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\Big]
- F_{fric,way}
= M_s \ddot{x}_n
\end{aligned}
$$

---

## 5. Second-order matrix form

$$
\mathbf{M}\ddot{\mathbf{q}} + \mathbf{C}\dot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{B}_u \mathbf{u}
$$

### 5.1 Input vector

$$
\mathbf{u} =
\begin{bmatrix}
\theta_{cmd} & T_{fric,sb} & F_{fric,way}
\end{bmatrix}^{\mathsf{T}}
$$

As of this revision the screw-nut interface is no longer an independent input; it is embedded directly in $\mathbf{K}$ and $\mathbf{C}$ (Sec. 3.4, Sec. 5.3-5.4). See Sec. 7 and Sec. 10 items 1-3.

### 5.2 Mass matrix

$$
\mathbf{M} = \operatorname{diag}\big(I_m,\; I_c,\; I_s,\; I_{sb},\; M_{screw},\; M_s\big)
$$

### 5.3 Stiffness matrix

$$
\mathbf{K} =
\begin{bmatrix}
k_c + k_{EM} + k_d & -k_c & 0 & 0 & 0 & 0 \\
-k_c & k_c + k_{s1} & -k_{s1} & 0 & 0 & 0 \\
0 & -k_{s1} & k_{s1} + k_{s2} + \left(\dfrac{L}{2\pi}\right)^2 k_{nut} & -k_{s2} & \dfrac{L}{2\pi}k_{nut} & -\dfrac{L}{2\pi}k_{nut} \\
0 & 0 & -k_{s2} & k_{s2} & 0 & 0 \\
0 & 0 & \dfrac{L}{2\pi}k_{nut} & 0 & k_{brg} + k_{nut} & -k_{nut} \\
0 & 0 & -\dfrac{L}{2\pi}k_{nut} & 0 & -k_{nut} & k_{nut}
\end{bmatrix}
$$

The screw-nut interface (Sec. 3.4) is now reflected through the lead ratio $L/2\pi$ directly onto $\theta_s$, coupling the rotational and axial blocks that were fully decoupled in the prior revision (Sec. 10 items 1-3). The sign of the four off-diagonal $L/2\pi \, k_{nut}$ terms was corrected on 2026-08-18: an earlier draft had the opposite sign here, which made the command-to-stage DC gain negative (a 180 deg phase offset at low frequency, opposite to the Sec. 1 convention that positive $\theta_s$ drives positive $x_n$). This form matches $F_{nut}$ in Sec. 3.4.

### 5.4 Damping matrix

$$
\mathbf{C} =
\begin{bmatrix}
c_c + c_{EM} & -c_c & 0 & 0 & 0 & 0 \\
-c_c & c_c + c_{s1} & -c_{s1} & 0 & 0 & 0 \\
0 & -c_{s1} & c_{s1} + c_{s2} + \left(\dfrac{L}{2\pi}\right)^2 c_{nut} & -c_{s2} & \dfrac{L}{2\pi}c_{nut} & -\dfrac{L}{2\pi}c_{nut} \\
0 & 0 & -c_{s2} & c_{s2} & 0 & 0 \\
0 & 0 & \dfrac{L}{2\pi}c_{nut} & 0 & c_{brg} + c_{nut} & -c_{nut} \\
0 & 0 & -\dfrac{L}{2\pi}c_{nut} & 0 & -c_{nut} & c_{nut}
\end{bmatrix}
$$

As of 2026-08-18, $\mathbf{C}_{11} = c_c + c_{EM}$ is applied by default rather than treated as an optional variant. Reason: mode 1 (176.7 Hz, the whole rotational chain moving almost rigidly against the $k_{EM}$ spring) barely engages the relative-velocity dampers $c_c, c_{s1}, c_{s2}$, so without $c_{EM}$ it comes out at $\zeta_1 \approx 5\times10^{-5}$ (a ~74 s settling time) regardless of how those placeholders are set. $c_{EM}$ acts directly on $\dot\theta_m$ and is the efficient lever on this mode; the value in Sec. 8/`model_parameters.json` is chosen to put $\zeta_1$ at 2%. Unlike the literal Sec. 3.1 relation ($T_{EM}$ depends on $\dot\theta_{cmd}-\dot\theta_m$), only the $\dot\theta_m$ half is implemented here, i.e. plain viscous damping of the rotor against a fixed frame, not against the moving command. The paired $\dot\theta_{cmd}$ feedforward column in $\mathbf{B}_u$ (below) is intentionally omitted: $\theta_{cmd}$ is a step staircase in the time-domain simulations this model feeds, and its derivative is a train of Dirac impulses, which is not representable on a finite time grid.

### 5.5 Input matrix

$$
\mathbf{B}_u =
\begin{bmatrix}
k_{EM} & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & -1 & 0 \\
0 & 0 & 0 \\
0 & 0 & -1
\end{bmatrix}
$$

---

## 6. First-order state space

State vector, 12 states:

$$
\mathbf{z} =
\begin{bmatrix}\mathbf{q} \\ \dot{\mathbf{q}}\end{bmatrix}
=
\begin{bmatrix}
\theta_m & \theta_c & \theta_s & \theta_{sb} & x_s & x_n &
\dot{\theta}_m & \dot{\theta}_c & \dot{\theta}_s & \dot{\theta}_{sb} & \dot{x}_s & \dot{x}_n
\end{bmatrix}^{\mathsf{T}}
$$

$$
\dot{\mathbf{z}} = \mathbf{A}\mathbf{z} + \mathbf{B}\mathbf{u}
$$

$$
\mathbf{A} =
\begin{bmatrix}
\mathbf{0}_{6\times6} & \mathbf{I}_{6\times6} \\[0.5ex]
-\mathbf{M}^{-1}\mathbf{K} & -\mathbf{M}^{-1}\mathbf{C}
\end{bmatrix},
\qquad
\mathbf{B} =
\begin{bmatrix}
\mathbf{0}_{6\times4} \\[0.5ex]
\mathbf{M}^{-1}\mathbf{B}_u
\end{bmatrix}
$$

### 6.1 Output

Stage position and the open-loop tracking error against the ideal kinematic command:

$$
y = x_n = \mathbf{C}_y \mathbf{z}, \qquad
\mathbf{C}_y = \begin{bmatrix} 0&0&0&0&0&1&0&0&0&0&0&0 \end{bmatrix}
$$

$$
e(t) = \frac{L}{2\pi}\theta_{cmd}(t) - x_n(t)
$$

---

## 7. Friction ports

Two of the three inputs are friction outputs (support bearing and guideway). As of this revision, the screw-nut interface is no longer a frozen-friction port: it is embedded directly in $\mathbf{K}$ and $\mathbf{C}$ as a linear reflected coupling between $\theta_s$ and the axial coordinates, scaled by the lead ratio $L/2\pi$ (Sec. 3.4, Sec. 5.3-5.4). This removes $T_{fric,nut}$ from $\mathbf{u}$ and resolves Sec. 10 items 1-3, at the cost of treating the nut contact as always-stuck (no slip). The remaining two friction outputs are functions of the states, so the system is linear only when they are frozen. In simulation the plant is evaluated as an LTI core wrapped by nonlinear feedback from these two ports.

| Port | Output | Relative velocity argument | Element |
|---|---|---|---|
| Support bearing | $T_{fric,sb}$ | $\omega_{sb} = \dot{\theta}_{sb}$ (grounded ring) | GMS or LuGre, rotational |
| Guideway | $F_{fric,way}$ | $v_{way} = \dot{x}_n$ | GMS or LuGre, translational |

Corresponding relative displacements, needed for presliding state propagation:

$$
\delta_{sb} = \theta_{sb}, \qquad
\delta_{way} = x_n
$$

The former screw-nut relative displacement, $x_n - x_s + \tfrac{L}{2\pi}\theta_s$, is retained only as the argument of the linear elastic/damping coupling in $\mathbf{K}$ and $\mathbf{C}$ (Sec. 3.4); it is no longer a friction-port state and does not need presliding propagation.

### 7.1 Simulation form

$$
\dot{\mathbf{z}} = \mathbf{A}\mathbf{z} + \mathbf{B}\,\mathbf{u}\big(\theta_{cmd}(t),\,\mathbf{z},\,\mathbf{w}\big),
\qquad
\dot{\mathbf{w}} = f_{GMS}(\mathbf{z}, \mathbf{w})
$$

where $\mathbf{w}$ collects the internal friction states, that is the $N$ Maxwell-slip element deflections $z_i$ per port together with their stick/slip flags. Integration requires a switched or stiff solver because each element toggles between stick and slip.

---

## 8. Parameter table

| Symbol | Description | Value | Unit | Source |
|---|---|---|---|---|
| $L$ | Screw lead | 1e-3 | m | BOM, KGT-F1-08-01 |
| $N_r$ | Rotor teeth | 50 | - | 1.8 deg step angle |
| $T_{hold}$ | Holding torque | TBD | N·m | Motor datasheet |
| $T_d$ | Detent torque amplitude | ~3e-3 | N·m | Datasheet / measurement |
| $k_{EM}$ | Electromagnetic stiffness | $N_r T_{hold}$ | N·m/rad | Derived |
| $k_d$ | Linearized detent stiffness | $4 N_r T_d$ | N·m/rad | Derived |
| $I_m$ | Rotor inertia | TBD | kg·m² | Motor datasheet |
| $I_c$ | Coupling inertia | TBD | kg·m² | Coupling datasheet |
| $I_s$ | Screw rotational inertia | TBD | kg·m² | Geometry |
| $I_{sb}$ | Support bearing inner ring inertia | TBD | kg·m² | Geometry |
| $M_{screw}$ | Screw axial mass | TBD | kg | Geometry |
| $M_s$ | Nut and stage mass | TBD | kg | CAD |
| $k_c$ | Coupling torsional stiffness | TBD | N·m/rad | Datasheet |
| $k_{s1}$ | Screw torsional stiffness, motor side | TBD | N·m/rad | Geometry |
| $k_{s2}$ | Screw torsional stiffness, bearing side | TBD | N·m/rad | Geometry |
| $k_{nut}$ | Screw-nut axial contact stiffness | TBD | N/m | Datasheet / Hertz |
| $k_{brg}$ | Support bearing axial stiffness | 7.5e6 to 15.3e6 | N/m | Barden duplex, light preload |
| $c_\bullet$ | Damping terms | TBD | - | Condensed from full model |

---

## 9. Implementation notes

1. Build $\mathbf{M}$, $\mathbf{C}$, $\mathbf{K}$, $\mathbf{B}_u$ as separate functions of a single parameter dictionary. Assemble $\mathbf{A}$ and $\mathbf{B}$ from them rather than hard-coding entries.
2. Keep $\mathbf{M}^{-1}$ as an explicit diagonal inverse. The mass matrix is diagonal by construction, so no factorization is needed.
3. Scale the state vector before eigenanalysis. Rotational and translational entries differ by roughly $2\pi/L \approx 6.28\times10^3$, and the raw $\mathbf{A}$ is badly conditioned.
4. Verify the rigid-body check first. With all friction ports zeroed, $k_d = 0$ and $k_{brg} = 0$, the axial block should show a free rigid-body mode.
5. Validate against the linearized detent case before enabling friction. The motor-only mode should appear near $\sqrt{(k_{EM}+k_d+k_c)/I_m}$.
6. Log the friction port outputs and the relative velocities alongside the states. The port signals are the quantities that later get compared to the identification data.

---

## 10. Open items in the derivation

These are unresolved points in the model as written, not transcription errors.

1. ~~Axial reaction pair at the screw-nut interface.~~ **Resolved, this revision.** The screw-nut traction is no longer an independently-injected port force; embedding it in $\mathbf{K}$ and $\mathbf{C}$ (Sec. 3.4) makes the reaction on $\theta_s$ and the reaction on $x_s$ symmetric by construction (both matrices are symmetric), so the force/torque balance closes automatically instead of needing a manual $\mathbf{B}_u$ patch.
2. ~~Interface double path.~~ **Resolved, this revision.** There is now exactly one screw-nut element: the linear coupling in $\mathbf{K}$/$\mathbf{C}$ between $\theta_s$, $x_s$, and $x_n$. The traction no longer exists as a separate path in parallel with $k_{nut}$.
3. ~~No kinematic constraint between $\theta_s$ and $x_n$.~~ **Resolved, with a new tradeoff.** $\theta_s$ and the axial block are now coupled directly through $\mathbf{K}$ and $\mathbf{C}$. This is only a stiff reflected spring/damper ($k_{nut}$, $c_{nut}$ scaled by $L/2\pi$), not a rigid constraint, but it assumes the nut interface never slips. The microslip/gross-slip behavior the original friction-port formulation was meant to capture at that interface is no longer represented anywhere in the model. If nut slip is significant, the port formulation should be restored for this interface specifically, and the direct $\mathbf{K}$/$\mathbf{C}$ coupling removed to avoid double-counting the compliance.
4. **Electromagnetic damping.** $c_{EM}$ appears in the constitutive relation but not in $\mathbf{C}$ or $\mathbf{B}_u$. Either drop it from Section 3.1 or add it per the variant in Section 5.4.
5. **Detent linearization validity.** $\sin(4N_r\theta_m) \approx 4N_r\theta_m$ holds only for $\theta_m \ll 1/(4N_r) \approx 5$ mrad, roughly a quarter step. For any motion larger than a single step the nonlinear form must be retained and $k_d$ removed from $\mathbf{K}$.
6. **Lead value.** The derivation notes reference a 2 mm lead. The BOM part number gives 1 mm. All numerical evaluation should use $L = 1$ mm.
