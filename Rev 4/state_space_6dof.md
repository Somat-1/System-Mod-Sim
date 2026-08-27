# Rev 4 Frictionless Six-Coordinate Drivetrain: State-Space Derivation

Reduced, frictionless lumped model of the open-loop stepper-driven positioning axis.
Retained degrees of freedom: motor rotation, coupling rotation, screw rotation, support-bearing rotation, screw axial translation, nut/stage axial translation.

An independent energy-based reconstruction and its frictionless Bode result
are provided in `Lagrange Derivation/README.md`.

This document derives the structural Rev 4 reference plant only. Nonlinear
friction laws, friction states, port Jacobians and alternative screw–nut
topologies belong to the LuGre sub-revisions:

- `lugre_friction/Rev 4.1/README.md`: LuGre replaces the structural
  `k_nut`/`c_nut` branch.
- `lugre_friction/Rev 4.2/README.md`: LuGre is added in parallel while the
  structural `k_nut`/`c_nut` branch is retained.

---

## 1. Conventions

Positive rotation `θ+` about the screw axis drives positive translation `x+`.
The lead `L` is the axial advance per screw revolution, so the ideal kinematic ratio is

$$
x = \frac{L}{2\pi}\,\theta , \qquad
\frac{\partial \theta}{\partial x} = \frac{2\pi}{L}
$$

All stiffness and damping elements are linear and act on relative coordinates.
Friction forces are zero in this reference model.

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

$k_{nut}$ and $c_{nut}$ form the linear structural screw–nut contact used by
the frictionless reference plant. They are coupled to screw rotation through
the lead ratio $L/2\pi$. The sign inside $F_{nut}$ is fixed by the Sec. 1
convention: $x_n-x_s-\tfrac{L}{2\pi}\theta_s$ is the nut's actual position
minus the ideal no-slip position $x_s+\tfrac{L}{2\pi}\theta_s$, so positive
$\theta_s$ drives positive $x_n$. The 2026-08-18 correction established this
sign after an earlier draft produced a negative command-to-stage DC gain.
$k_{brg}$ is the grounded axial rigidity of the support-bearing pair.

---

## 4. Equations of motion

### 4.1 Compact form

$$
\begin{aligned}
\text{(1) Motor rotor:} && T_{EM} - T_{det} - T_{c,\text{react}} &= I_m \ddot{\theta}_m \\
\text{(2) Coupling:} && T_{c,\text{drive}} - T_{s1,\text{react}} &= I_c \ddot{\theta}_c \\
\text{(3) Screw rotation:} && T_{s1,\text{drive}} - \tfrac{L}{2\pi}F_{nut} - T_{s2,\text{react}} &= I_s \ddot{\theta}_s \\
\text{(4) Support bearing:} && T_{s2,\text{drive}} &= I_{sb} \ddot{\theta}_{sb} \\
\text{(5) Screw axial:} && F_{nut} - F_{ax,\text{react}} &= M_{screw} \ddot{x}_s \\
\text{(6) Nut and stage:} && -F_{nut} &= M_s \ddot{x}_n
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
= I_{sb} \ddot{\theta}_{sb} \\[1.5ex]
\text{(5)}\quad
& \Big[k_{nut}\big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\big) + c_{nut}\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\Big]
- \Big[k_{brg}x_s + c_{brg}\dot{x}_s\Big]
= M_{screw} \ddot{x}_s \\[1.5ex]
\text{(6)}\quad
& -\Big[k_{nut}\big(x_n - x_s - \tfrac{L}{2\pi}\theta_s\big) + c_{nut}\big(\dot{x}_n - \dot{x}_s - \tfrac{L}{2\pi}\dot{\theta}_s\big)\Big]
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
\mathbf{u}=\theta_{cmd}
$$

The baseline has one command input. The structural screw–nut contact is
embedded directly in $\mathbf K$ and $\mathbf C$; it is not an input.

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

The screw–nut structural contact is reflected through $L/2\pi$ onto
$\theta_s$, coupling the rotational and axial blocks. The four off-diagonal
$L/2\pi\,k_{nut}$ signs enforce the convention that positive $\theta_s$
drives positive $x_n$. This form matches $F_{nut}$ in Sec. 3.4.

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
k_{EM}\\
0\\
0\\
0\\
0\\
0
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
\mathbf{0}_{6\times1} \\[0.5ex]
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

## 7. Scope boundary for friction extensions

No friction state or friction generalized force belongs to the Rev 4
reference equations above. The reference transfer function is evaluated with
all friction forces equal to zero. The implementation retains two zeroed
legacy input columns for support-bearing torque and guideway force so older
analysis scripts remain compatible; those columns are not part of the
one-input frictionless derivation.

The two nonlinear alternatives are documented independently:

1. **Rev 4.1 — replacement experiment.** The nut LuGre port replaces
   $k_{nut}$ and $c_{nut}$. This tests the interpretation that LuGre carries
   the complete screw–nut contact compliance and friction.
2. **Rev 4.2 — parallel experiment.** The structural $k_{nut}/c_{nut}$ path
   remains load-bearing and the nut LuGre port supplies an additional
   pre-rolling friction force across the same relative coordinate.

Neither alternative changes the definition of the frictionless reference
plant in this document.

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
4. Verify the rigid-body check first. With $k_d=0$ and $k_{brg}=0$, the axial block should show a free rigid-body mode.
5. Validate the optional linearized-detent case separately. The motor-only mode should appear near $\sqrt{(k_{EM}+k_d+k_c)/I_m}$.

---

## 10. Open items in the derivation

These are unresolved points in the model as written, not transcription errors.

1. **Structural nut-contact identification.** $k_{nut}$ and $c_{nut}$ remain
   placeholders until identified from contact data or an assembled-axis test.
2. **Detent linearization validity.** $\sin(4N_r\theta_m) \approx
   4N_r\theta_m$ holds only for $\theta_m \ll 1/(4N_r) \approx 5$ mrad,
   roughly a quarter step. For larger motion the nonlinear form must be
   retained and $k_d$ removed from $\mathbf K$.
3. **Lead value.** The derivation notes reference a 2 mm lead. The BOM part
   number gives 1 mm. All numerical evaluation uses $L=1$ mm.
