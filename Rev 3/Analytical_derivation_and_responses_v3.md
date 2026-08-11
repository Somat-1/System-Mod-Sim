# Revision 3: Analytical Derivation and Executable Responses

This is the executable companion to [the Revision 3 model specification](ball_screw_stage_dynamic_derivation_v3.html). It derives the ten-coordinate plant, audits the two-DOF reduction, and compares the friction cases.

> **Reproducibility boundary.** Browser edits update derived values and the live transfer panel. Rebuild to refresh static figures and simulated metrics.

Amber cells still require identification.

<details>
<summary>How to reproduce both HTML documents and every figure</summary>

From the `Rev 3` folder, run:

```text
python build_model_documentation.py
```

The build writes both HTML documents and every figure, including the response overlay, reduction audit, and position sweep. There is no second simulation or rendering script.

</details>

## 1. Model hierarchy and case map

There are two structural models and seven executed friction cases.

| Layer | Coordinates | Purpose |
|---|---:|---|
| Full Revision 3 plant | 10 mechanical DOFs | Trace every inertia, compliance, and nut-interface sign; expose discarded modes |
| Reduced plant | 2 mechanical DOFs | Execute the ≤900 Hz response and keep nut/guideway friction identifiable |
| Friction states | 0, 1/site, or 4/site | Internal constitutive memory; these are states, not additional mechanical DOFs |

| Case | Active sites | Law | Role |
|---|---|---|---|
| 0 | none | none | frictionless modal baseline |
| A | drivetrain + guideway | LuGre | guideway hypothesis with common drivetrain loss |
| A2 | drivetrain + guideway | GMS | topology-matched alternative to A |
| B | drivetrain + nut rolling + nut microslip | LuGre | nut hypothesis |
| B2 | drivetrain + nut rolling + nut microslip | GMS | topology-matched alternative to B |
| C | all four ports | LuGre | combined hypothesis |
| C2 | all four ports | GMS | topology-matched alternative to C |

$F_{f,d}$ is active in every friction case. Case 0 remains the only frictionless run.

## 2. Entry parameters

Open only the parameter group you need. Browser edits update derived values and the live transfer panel. Rebuild to refresh static figures.

<details class="parameter-group">
<summary>Geometry, reduced plant, and excitation</summary>

| Symbol | Meaning | Executed value | Unit |
|---|---|---:|---|
| $L$ | screw lead | [[input:lead=1.000e-3]] | m/rev |
| $N_r$ | rotor teeth | [[input:rotor_teeth=50]] | – |
| $r=L/(2\pi)$ | transmission ratio, derived | [[derived:transmission_ratio=1.59155e-4]] | m/rad |
| $m_d=J_\Sigma/r^2$ | reflected drivetrain mass, derived | [[derived:reduced_drive_mass=121.994]] | kg |
| $m_s$ | reduced nut + stage effective mass | [[input:reduced_stage_mass=0.600]] | kg |
| $T_{max}$ | rated-current holding torque, 0674A | [[input:holding_torque=0.060]] | N·m |
| $\hat T_{det}$ | published detent torque, enabled | [[input:detent_torque=0.005]] | N·m |
| $\phi_{det}$ | detent phase at the stable report origin | [[assumed:detent_phase=0.0]] | rad |
| $K_m=N_rT_{max}/r^2$ | commutation tangent, derived | [[derived:magnetic_stiffness=1.18435e8]] | N/m |
| $K_{det}=4N_r\hat T_{det}\cos\phi_{det}/r^2$ | detent tangent, derived | [[derived:detent_stiffness=3.94784e7]] | N/m |
| $k_{ax}$ | measured reduced axial-path stiffness | [[input:reduced_axial_stiffness=1.140e7]] | N/m |
| $c_{ax}$ | retained structural damping | [[assumed:axial_damping=55.0]] | N·s/m |
| $\zeta_m$ | provisional open-loop drive damping ratio | [[assumed:electromagnetic_zeta=0.05]] | – |
| $n_\mu$ | external STEP/DIR microstep divisor | [[assumed:microstep_divisor=64]] | – |
| $p_{step}$ | 1.8° full-step linear pitch, derived | [[derived:full_step_pitch=5.000e-6]] | m |
| $p_{step}/4$ | maximum command increment, derived | [[derived:quarter_step_bound=1.250e-6]] | m |
| $p_{step}/n_\mu$ | executed STEP/DIR quantum, derived | [[derived:command_step=7.81250e-8]] | m |
| $p_{step}/256$ | optional interpolated quantum, derived | [[derived:interpolated_step=1.95313e-8]] | m |
| axial play | accuracy grade O | 0.0 | m |

</details>

<details class="parameter-group">
<summary>Ten-DOF inertias and masses</summary>

| Symbol | Meaning | Executed value | Unit |
|---|---|---:|---|
| $J_m$ | 0674A motor rotor inertia, datasheet | [[input:J_m=9.000e-7]] | kg·m² |
| $J_c$ | coupling inertia estimate from the 23.8 g annulus | [[assumed:J_c=1.180e-6]] | kg·m² |
| $m_c$ | coupling mass, datasheet | 0.0238 | kg |
| $L_s$ | screw length L2 | [[input:screw_length=0.320]] | m |
| $d_s$ | nominal screw diameter | [[input:screw_diameter=8.000e-3]] | m |
| $\rho_s$ | steel density | [[assumed:screw_density=7850]] | kg/m³ |
| $J_s$ | complete screw polar inertia, derived | [[derived:screw_inertia=1.01014e-6]] | kg·m² |
| $J_{s1}=J_{s2}=J_{s3}$ | one-third screw inertia | [[derived:screw_segment_inertia=3.36712e-7]] | kg·m² |
| $m_{screw}$ | complete screw mass, derived | [[derived:screw_mass=0.126267]] | kg |
| $m_b=m_e=m_f$ | one-third axial screw mass | [[derived:screw_segment_mass=0.0420890]] | kg |
| $m_n$ | nut body mass | [[assumed:m_n=0.050]] | kg |
| $m_{stage}$ | stage body portion used with $m_n$ | [[assumed:m_stage=0.550]] | kg |

</details>

<details class="parameter-group">
<summary>Full-model stiffnesses and damping</summary>

| Symbol | Meaning | Executed value | Unit |
|---|---|---:|---|
| $k_c$ | coupling series torsion, datasheet | [[input:k_c_series=68.7549]] | N·m/rad |
| $k_{c1}=k_{c2}=2k_c$ | two equal half-springs, derived | [[derived:k_c_half=137.510]] | N·m/rad |
| $k_{\theta a}$ | screw torsion before nut | [[assumed:k_theta_a=211.0]] | N·m/rad |
| $k_{\theta b}$ | screw torsion beyond nut | [[assumed:k_theta_b=211.0]] | N·m/rad |
| $k_{brg}$ | support-bearing axial stiffness | [[assumed:k_brg=2.500e7]] | N/m |
| $k_{sha}$ | screw axial stiffness before nut | [[assumed:k_sha=6.700e7]] | N/m |
| $k_{shb}$ | screw axial stiffness beyond nut | [[assumed:k_shb=2.000e8]] | N/m |
| $k_{ball}$ | ball-contact stiffness from closure | [[assumed:k_ball=4.387e7]] | N/m |
| $k_{mnt}$ | nut-mount stiffness | [[assumed:k_mnt=1.000e8]] | N/m |
| $\zeta_{int}$ | proportional element damping ratio | [[assumed:zeta_internal=0.010]] | – |

</details>

The screw uses $m=\rho\pi d_s^2L_s/4$ and $J_s=md_s^2/8$. The three screw coordinates receive equal thirds. The coupling value remains an estimate because its datasheet publishes 23.8 g mass but not polar inertia. No target value is imposed on $m_d$.

## 3. Kinematic diagram and degrees of freedom

![Ten-DOF topology, retained two-DOF reduction, and rejected one-DOF collapse](rendered_assets/kinematic_diagram.svg)

The main axial load path is ground, $k_{brg}$, $u_b$, $k_{sha}(x_s)$, $u_e$, the screw transformer and $k_{ball}$, $u_n$, $k_{mnt}$, then $x_s$. The $u_f$ and $\theta_{s3}$ coordinates are overhang stubs. They do not carry stage load.

The figure uses two independent encodings. Box fill shows whether a mass migrates to $m_d$, migrates to $m_s$, or is dropped. Spring stroke shows whether its compliance is retained or discarded. Thus a dropped axial mass can sit between retained springs without implying that its inertia survives. The one-DOF panel is rejected because it sets the nut-port velocity to zero and merges the friction sites.

The honest full enumeration is ten mechanical DOFs:

| Index | Coordinate | Physical motion | Native unit |
|---:|---|---|---|
| 1 | $\theta_m$ | motor rotor rotation | rad |
| 2 | $\theta_c$ | coupling body rotation | rad |
| 3 | $\theta_{s1}$ | screw rotation at drive end | rad |
| 4 | $\theta_{s2}$ | screw rotation at nut plane | rad |
| 5 | $\theta_{s3}$ | screw rotation beyond nut | rad |
| 6 | $u_b$ | axial screw motion at bearing | m |
| 7 | $u_e$ | axial screw motion at nut plane | m |
| 8 | $u_f$ | axial screw motion beyond nut | m |
| 9 | $u_n$ | nut-body axial motion | m |
| 10 | $x_s$ | stage axial motion | m |

$x_{cmd}$ is an imposed electromagnetic-field input. It has no independently solved inertia and is not a DOF. LuGre bristle deflections and GMS element forces are constitutive internal states, not mechanical generalized coordinates.

## 4. Full ten-DOF derivation

Define

$$\mathbf q=[\theta_m,\theta_c,\theta_{s1},\theta_{s2},\theta_{s3},u_b,u_e,u_f,u_n,x_s]^T.$$

The full linear mechanical system has the form

$$\mathbf M\ddot{\mathbf q}+\mathbf C\dot{\mathbf q}+\mathbf K\mathbf q=\mathbf b\,x_{cmd}+\mathbf Q_f.$$

<details>
<summary>Step 1: kinetic energy and diagonal mass matrix</summary>

Because every coordinate is defined at a physical inertia or lumped axial mass, no kinematic substitution is made before the energy is written:

$$
\begin{aligned}
\mathcal T={}&\tfrac12J_m\dot\theta_m^2+\tfrac12J_c\dot\theta_c^2
+\tfrac12J_{s1}\dot\theta_{s1}^2+\tfrac12J_{s2}\dot\theta_{s2}^2
+\tfrac12J_{s3}\dot\theta_{s3}^2\\
&+\tfrac12m_b\dot u_b^2+\tfrac12m_e\dot u_e^2+\tfrac12m_f\dot u_f^2
+\tfrac12m_n\dot u_n^2+\tfrac12m_{stage}\dot x_s^2.
\end{aligned}
$$

Therefore

$$\mathbf M=\operatorname{diag}(J_m,J_c,J_{s1},J_{s2},J_{s3},m_b,m_e,m_f,m_n,m_{stage}).$$

The diagonal form is a consequence of coordinate choice, not an assumption that the branches are uncoupled. Their coupling enters through potential energy.

</details>

<details>
<summary>Step 2: every elastic deflection and the complete potential energy</summary>

The elementary deformations are

$$
\begin{array}{lll}
d_m=\theta_m-\theta_{cmd}, & d_{c1}=\theta_m-\theta_c, & d_{c2}=\theta_c-\theta_{s1},\\
d_{\theta a}=\theta_{s1}-\theta_{s2}, & d_{\theta b}=\theta_{s2}-\theta_{s3}, & d_{brg}=u_b,\\
d_{sha}=u_e-u_b, & d_{shb}=u_f-u_e, & d_{mnt}=u_n-x_s.
\end{array}
$$

The thread/ball contact is the only mixed rotational–axial deformation:

$$\boxed{\delta_n=u_n-u_e-r\theta_{s2}}.$$

The complete linearized potential is

$$
\begin{aligned}
\mathcal V={}&\tfrac12 k_m d_m^2+\tfrac12k_{c1}d_{c1}^2+\tfrac12k_{c2}d_{c2}^2
+\tfrac12k_{\theta a}d_{\theta a}^2+\tfrac12k_{\theta b}d_{\theta b}^2\\
&+\tfrac12k_{brg}d_{brg}^2+\tfrac12k_{sha}d_{sha}^2+\tfrac12k_{shb}d_{shb}^2
+\tfrac12k_{ball}\delta_n^2+\tfrac12k_{mnt}d_{mnt}^2.
\end{aligned}
$$

For any scalar element $\tfrac12k(\mathbf h^T\mathbf q)^2$, its matrix contribution is $k\mathbf h\mathbf h^T$. This outer-product rule is the direct bridge from the physical diagram to code.

</details>

<details>
<summary>Step 3: nut-interface virtual work and sign audit</summary>

Let the ball-contact force be positive in extension:

$$F_n=k_{ball}\delta_n+c_{ball}\dot\delta_n.$$

The gradient of the deformation is

$$\mathbf h_n=\frac{\partial\delta_n}{\partial\mathbf q}
=[0,0,0,-r,0,0,-1,0,+1,0]^T.$$

The contact’s generalized force on the coordinates is $-F_n\mathbf h_n$. Hence:

- on $\theta_{s2}$: $+rF_n$;
- on $u_e$: $+F_n$;
- on $u_n$: $-F_n$.

The corresponding stiffness and damping contributions are

$$\mathbf K_n=k_{ball}\mathbf h_n\mathbf h_n^T,\qquad
\mathbf C_n=c_{ball}\mathbf h_n\mathbf h_n^T.$$

This construction guarantees equal-and-opposite internal work and symmetric passive matrices. It also prevents the common sign mistake of applying the same axial force direction to $u_e$ and $u_n$.

</details>

<details>
<summary>Step 4: Rayleigh dissipation and the electromagnetic damping repair</summary>

For each structural element with deformation $d_j=\mathbf h_j^T\mathbf q$, use

$$\mathcal R_j=\tfrac12c_j\dot d_j^2,\qquad \mathbf C_j=c_j\mathbf h_j\mathbf h_j^T.$$

The executable provisional values use $c_j=2\zeta_{int}\sqrt{k_jm_{rel,j}}$ with $\zeta_{int}=0.01$. These are damping assumptions, not identified loss factors.

Revision 2 exposed a separate missing term: the ideal position-source stepper model supplied restoring stiffness but no effective drive damping. That lossless oscillator rang around every command. The same phenomenological term is retained here:

$$c_{\theta m}=2\zeta_m\sqrt{k_mJ_\Sigma},\qquad
c_m=\frac{c_{\theta m}}{r^2}=2\zeta_m\sqrt{K_mm_d}.$$

It enters as $-c_{\theta m}\dot\theta_m$ in the full rotor equation and $-c_m\dot x_d$ in the reduced drive equation. The executed $\zeta_m=0.05$ is not identified. The sensitivity plot spans 0.02 to 0.50. The TMC2209 can use StealthChop2 or SpreadCycle, so driver mode and current settings are required before assigning a physical damping ratio. See the [official TMC2209 datasheet](https://www.analog.com/media/en/technical-documentation/data-sheets/TMC2209_datasheet_rev1.09.pdf).

</details>

<details>
<summary>Step 5: scalar equations recovered from Lagrange’s equation</summary>

Using $\frac{d}{dt}(\partial\mathcal T/\partial\dot q_i)-\partial\mathcal T/\partial q_i+\partial\mathcal V/\partial q_i+\partial\mathcal R/\partial\dot q_i=Q_{f,i}$ gives:

$$J_m\ddot\theta_m=T_{mag}+T_{det}-c_{\theta m}\dot\theta_m-k_{c1}(\theta_m-\theta_c)-c_{c1}(\dot\theta_m-\dot\theta_c)-T_{h1}-T_{mb},$$

$$J_c\ddot\theta_c=k_{c1}(\theta_m-\theta_c)+c_{c1}(\dot\theta_m-\dot\theta_c)+T_{h1}-k_{c2}(\theta_c-\theta_{s1})-c_{c2}(\dot\theta_c-\dot\theta_{s1})-T_{h2},$$

$$J_{s1}\ddot\theta_{s1}=k_{c2}(\theta_c-\theta_{s1})+c_{c2}(\dot\theta_c-\dot\theta_{s1})+T_{h2}-k_{\theta a}(\theta_{s1}-\theta_{s2})-T_{brg},$$

$$J_{s2}\ddot\theta_{s2}=k_{\theta a}(\theta_{s1}-\theta_{s2})-k_{\theta b}(\theta_{s2}-\theta_{s3})+rF_n-T_{f,n}-T_{f,r},$$

$$J_{s3}\ddot\theta_{s3}=k_{\theta b}(\theta_{s2}-\theta_{s3}),$$

$$m_b\ddot u_b=-k_{brg}u_b-c_{brg}\dot u_b+k_{sha}(u_e-u_b)+c_{sha}(\dot u_e-\dot u_b),$$

$$m_e\ddot u_e=-k_{sha}(u_e-u_b)-c_{sha}(\dot u_e-\dot u_b)+k_{shb}(u_f-u_e)+c_{shb}(\dot u_f-\dot u_e)+F_n-F_{f,n},$$

$$m_f\ddot u_f=-k_{shb}(u_f-u_e)-c_{shb}(\dot u_f-\dot u_e),$$

$$m_n\ddot u_n=-F_n+F_{f,n}-k_{mnt}(u_n-x_s)-c_{mnt}(\dot u_n-\dot x_s),$$

$$m_{stage}\ddot x_s=k_{mnt}(u_n-x_s)+c_{mnt}(\dot u_n-\dot x_s)-F_{f,g}.$$

The nut contact microslip port uses

$$v_n=r\dot\theta_{s2}+\dot u_e-\dot u_n=-\dot\delta_n,\qquad T_{f,n}=rF_{f,n}.$$

It applies $-T_{f,n}$ to $\theta_{s2}$, $-F_{f,n}$ to $u_e$, and $+F_{f,n}$ to $u_n$. The port is internal and power consistent before reduction.

Gross nut rolling drag is separate. It uses $v_r=\dot x_d$ in the reduced model and $T_{f,r}=rF_{f,r}$ in the full-model bookkeeping. It remains active during common drive and stage motion. This prevents the gross 5 N level from being assigned to an elastic deformation rate that vanishes at steady motion.

</details>

## 5. Stepper input: nonlinear law, linearization, and bound

The commanded linear position maps to field angle through $\theta_{cmd}=x_{cmd}/r$. With $N_r$ rotor teeth,

The 0674A motor is assumed to run at rated current, so its 0.060 N·m holding torque is used as $T_{max}$.

$$T_{mag}=T_{max}\sin\!\left[N_r(\theta_{cmd}-\theta_m)\right].$$

Under the reduced coordinate $x_d=r\theta_m$,

$$F_{mag}=\frac{T_{max}}r\sin\!\left[\frac{N_r}{r}(x_{cmd}-x_d)\right]
=F_{max}\sin[\kappa(x_{cmd}-x_d)].$$

The small-signal stiffness is

$$K_m=\left.\frac{\partial F_{mag}}{\partial(x_{cmd}-x_d)}\right|_0
=\frac{N_rT_{max}}{r^2}=F_{max}\kappa.$$

For the current inputs, one full step is [[derived:full_step_pitch=5.000e-6]] m. The nonlinear audit bound remains [[derived:quarter_step_bound=1.250e-6]] m. The executed response step is [[derived:command_step=7.81250e-8]] m at 64 microsteps per full step. The optional 256-interpolated quantum is [[derived:interpolated_step=1.95313e-8]] m. The TMC2209 supports 8, 16, 32, or 64 STEP/DIR settings and interpolation to 256, according to the [official datasheet](https://www.analog.com/media/en/technical-documentation/data-sheets/TMC2209_datasheet_rev1.09.pdf).

At the quarter-step audit bound, $\kappa e=\pi/8$ and

$$\frac{\sin(\pi/8)}{\pi/8}=0.9745.$$

Thus the sine law is only 2.55% below its tangent at the largest commanded increment. It can shift amplitude and frequency slightly, but it was not the source of the old sustained oscillation; missing damping was.

### 5.1 Why the 150–250 Hz stepper feature is difficult to see

![Low-frequency stepper-mode visibility versus damping and selected output](rendered_assets/stepper_resonance_visibility.svg)

The enabled detent and component-derived inertia give a low pole of [[derived:mode_1_hz=180.02]] Hz. It is the common drive mode, approximately

$$f_m\approx\frac{1}{2\pi}\sqrt{\frac{K_m+K_{det}}{m_d}}.$$

The executed $\zeta_m=0.05$ exposes this pole. The plotted output still matters: stage motion $X_s/X_{cmd}$ has weaker low-mode participation than the internal drive coordinate. The 0.02, 0.05, 0.10, and 0.50 curves are sensitivity cases, not identified driver properties.

### 5.2 Rotor-equivalent drive and stage transfer functions

![Command-to-rotor, command-to-stage, and rotor-to-stage Bode functions](rendered_assets/rotor_stage_transfer_functions.svg)

For the frictionless two-DOF linear model define

$$
a(s)=m_ds^2+(c_m+c_{ax})s+K_m+K_{det}+k_{ax},\qquad
b(s)=-(c_{ax}s+k_{ax}),\qquad
d(s)=m_ss^2+c_{ax}s+k_{ax},
$$

and $\Delta(s)=a(s)d(s)-b(s)^2$. The three relevant transfer functions are

$$\boxed{\frac{X_d}{X_{cmd}}=\frac{K_md(s)}{\Delta(s)}},$$

$$\boxed{\frac{X_s}{X_{cmd}}=\frac{-K_mb(s)}{\Delta(s)}},$$

$$\boxed{\frac{X_s}{X_d}=\frac{-b(s)}{d(s)}
=\frac{c_{ax}s+k_{ax}}{m_ss^2+c_{ax}s+k_{ax}}}.$$

$X_d/X_{cmd}$ is the clearest view of the low rotor/drive pole. $X_s/X_{cmd}$ contains both plant poles but can show weaker low-mode participation. $X_s/X_d$ treats rotor-equivalent drive motion as a prescribed input, so the common motor/magnetic pole cancels and only the stage-following dynamics remain. Therefore rotor-to-stage alone cannot be used to prove that the lower mode is absent.

The panel below is computed in the browser. Changing lead, component inertias, screw geometry, holding torque, detent torque, or plant parameters updates the derived values and these curves. Static SVGs still require a Python rebuild.

<div id="live-transfer-panel" class="live-transfer-panel" data-live-transfer-plots></div>

<details>
<summary>Enabled detent torque and its linearization</summary>

The published $\hat T_{det}=0.005$ N·m is enabled in the nonlinear simulation:

$$T_{det}(\theta_m)=-\hat T_{det}\sin(4N_r\theta_m+\phi_{det}).$$

Around an equilibrium $\theta_0$, it contributes the position-dependent tangent stiffness

$$k_{det,lin}=4N_r\hat T_{det}\cos(4N_r\theta_0+\phi_{det}).$$

The report origin uses $\phi_{det}=0$, a stable detent equilibrium. Its current tangent is $K_{det}=$ [[derived:detent_stiffness=3.94784e7]] N/m.

Balancing commutation and detent torque gives the worst-case drive-only bound

$$|x_{err}|_{max}=\frac{r}{N_r}\sin^{-1}\!\left(\frac{\hat T_{det}}{T_{max}}\right)=266\ \mathrm{nm}.$$

Its spatial period is $p_{step}/4=1.25\ \mu$m. The earlier 371 nm estimate used the 0.043 N·m motor variant. With the executed 0.060 N·m motor, 266 nm is the consistent value. Compliance and friction modify the realized stage error, but the periodic term must not be absorbed into friction identification.

Detent alone is not the whole “stepper resonance” mechanism. Current-loop dynamics, back-EMF, driver delay, current quantization, and mechanical damping set the response amplitude and stability. A higher-fidelity model still needs detent phase, an identified damping ratio, and current-controller states.

</details>

<details>
<summary>Cautious comparison with the physical modal-test folder</summary>

The motor-excited chirp results in `TempScripts/Modal Comparison` show a broad 155–190 Hz feature in all payload configurations. Only the +1 kg up/down pair, about 159–160 Hz, passed the report’s 3× local-floor criterion. Its weak payload dependence is compatible with a motor/base/grounded mode, but does not uniquely prove detent resonance.

Important limitations prevent treating all plotted features as identified modes:

- the chirp normalization divides acceleration by $f^2$, amplifying the low-SNR region below about 250 Hz;
- the second-pass impact analysis starts at 200 Hz, so it cannot validate most of the 155–190 Hz band;
- its +1 kg PLA result retained only 4 of 15 impacts and reported a weak 256.3 Hz candidate;
- DitherV2 applies a high-pass centered at 250 Hz and therefore intentionally removes the frequency range in question;
- the fixed approximately 345 Hz feature does not shift with payload and was already flagged as probable chopper/measurement contamination;
- the 0 kg 685–700 Hz chirp feature is direction-dependent and did not clear the formal prominence threshold.

The physical data support comparison with the 178 Hz drive mode, but do not identify detent phase or damping.

</details>

### 5.3 Can the present model reproduce every measured feature?

**No - not every feature.** The updated two-DOF modes are [[derived:mode_1_hz=180.02]] Hz and [[derived:mode_2_hz=695.56]] Hz. They now align with the broad 155–190 Hz response and the 685–700 Hz band. The transfer matrix still has only two mechanical modal pairs:

$$\det\!\left(\mathbf M_rs^2+\mathbf C_rs+\mathbf K_r\right)=0,$$

which is fourth order in $s$. The rebuilt ten-DOF model predicts 179, 676, 1717, and 2839 Hz below 3 kHz. It still has no independent mode near 256, 345, or 1007–1044 Hz. Friction can move or damp existing poles, but cannot create missing coordinates.

<details>
<summary>Evidence audit from the referenced Modal Comparison folder</summary>

| Observed feature | Local-test evidence and caution | Present-model interpretation |
|---|---|---|
| Broad 155-190 Hz response | Seen across chirp runs, but only the +1 kg up/down pair at about 159-160 Hz clears the 3x local-floor rule. It changes little with payload. | The updated 178 Hz drive mode lies in this band. Matching still requires the same measured input/output definition. |
| 256.3 Hz hammer candidate | Appears in the +1 kg PLA second-pass result, whose accepted-tap accounting is weak and internally cautioned. | Insufficient evidence for a new plant pole; repeat with a reliable mount before adding a coordinate. |
| Approximately 345 Hz notch | Fixed across payload and sweep direction; the local report identifies chopper/ambient contamination. | Must not be fitted as a mechanical mode. |
| Approximately 592-614 Hz hammer and 685-700 Hz chirp features | The tests use different inputs, outputs, boundaries, and normalization. | The updated 676–696 Hz relative mode is a plausible match. Co-located FRFs are still needed. |
| Approximately 1007-1044 Hz PLA-mounted impact peaks | Coherent hammer peaks appear with added PLA-mounted payloads. | Likely requires an explicit payload/bracket or sensor-mount compliance; that coordinate is absent from both the reduced model and the present rigid-stage output. |
| Dither candidates | No configuration passes both the prominence and visible-fade criteria in the DitherV2 report. | Do not use these candidates for parameter fitting yet. |

The audited sources are Results Summary/modal_testing_summary.pdf, Results Summary/chirp_results_summary_v2.pdf, MotorExcitation/Chirp Tests/chirp_results/summary.md, V2 Testing/v2_modal_second_pass/summary.md, and MotorExcitation/DitherV2/ditherv2_coherent_results/summary.md under TempScripts/Modal Comparison.

</details>

<details>
<summary>Minimum defensible extensions, in order</summary>

1. **Match the measured transfer function first.** Generate impact inertance $A_s/F_{impact}$ at the accelerometer point and motor-excited $A_s/X_{cmd}$ or $A_s/I_{cmd}$, rather than comparing both tests directly with displacement ratio $X_s/X_{cmd}$. Poles are shared only when the same structure and boundary are excited; modal residues and antiresonances depend strongly on input/output location.
2. **Add a grounded base/support coordinate** $x_b$ with identified $m_b,k_b,c_b$, and express guideway/bearing forces relative to $x_b$ instead of an immovable ground. This is the leading candidate for a broad, weakly payload-dependent 155-190 Hz mode.
3. **Add payload/bracket compliance** between the stage body and accelerometer/payload coordinate. This is the leading candidate for the approximately 1 kHz PLA-mounted impact peaks.
4. **Identify the motor dynamics.** Detent amplitude is now sourced and enabled. Its phase, current-controller dynamics, and low-mode damping still require dedicated tests.
5. **Add transverse/rocking coordinates only if cross-axis measurements support them.** Rail bending, stage pitch/yaw, screw lateral bending, and bearing-housing motion are excluded by the present one-axis topology but can appear in an accelerometer FRF through sensor misalignment or structural coupling.

This is a representability issue, not a reason to add arbitrary fitted resonances. Each new coordinate needs a mode-shape or mass/stiffness experiment that distinguishes it from measurement artifacts.

</details>

## 6. Reduction from ten DOFs to two

<details>
<summary>Why two model DOFs when the end effector moves along only one axis?</summary>

The stage has one measured output direction, but output dimension is not the same as system DOF count. $x_s$ is the end-effector translation; $x_d$ is an internal reflected rotor/screw coordinate on the same axis. Finite axial compliance permits $x_d-x_s\ne0$, so two independent initial positions and velocities are required.

A two-mass system connected by a spring has two DOFs even when both masses move along the same line and only the second mass is measured. Here the relative coordinate $x_d-x_s$ stores axial elastic energy and produces the [[derived:mode_2_hz=695.56]] Hz mode.

If the complete drivetrain is imposed rigidly, then $x_d=x_s=x$ and a legitimate one-DOF model results:

$$
(m_d+m_s)\ddot x+c_m\dot x+(K_m+K_{det})x=K_mx_{cmd}-F_{f,aggregate}.
$$

That model retains the [[derived:mode_1_hz=180.02]] Hz common-motion pole,

$$f_1\approx\frac{1}{2\pi}\sqrt{\frac{K_m+K_{det}}{m_d+m_s}},$$

but it removes the relative [[derived:mode_2_hz=695.56]] Hz mode, makes the modeled nut-port velocity $\dot x_d-\dot x_s$ identically zero, and merges the remaining friction sites. It is suitable only below the axial mode.

</details>

<details>
<summary>Step 1: retain compliance even when an internal mass is collapsed</summary>

Collapsing a coordinate does not mean deleting its spring. The load-path stiffness is retained as a series compliance:

$$\boxed{\frac1{k_{ax}}=\frac1{k_{brg}}+\frac1{k_{sha}}+\frac1{k_{ball}}+\frac1{k_{mnt}}}.$$

At the executable 150 mm position,

$$
\frac1{25.0\times10^6}+\frac1{67.0\times10^6}
+\frac1{43.87\times10^6}+\frac1{100\times10^6}
\approx\frac1{11.4\times10^6}.
$$

$k_{ball}$ is derived from the remainder after selecting the highlighted $k_{brg}=25$ MN/m. The published 15 MN/m upper value alone consumes more compliance than the measured chain permits; this is the original Revision 3 closure warning.

</details>

<details>
<summary>Step 2: reflect rotational inertia and define retained coordinates</summary>

With $x_d=r\theta$ and equal kinetic energy,

$$\tfrac12J_\Sigma\dot\theta^2=\tfrac12m_d\dot x_d^2
\quad\Rightarrow\quad
\boxed{m_d=J_\Sigma/r^2}.$$

The stage-side retained mass is $m_s=m_n+m_{stage}=0.60$ kg in this implementation. The nut mass therefore migrates to the stage coordinate. The axial screw masses $m_b$, $m_e$, and $m_f$ are dropped, while the four series load-path compliances are retained. The retained coordinates are

$$\mathbf x=[x_d,x_s]^T.$$

$x_d$ represents the collapsed motor/coupling/screw drive side; $x_s$ remains the measured stage coordinate.

</details>

<details>
<summary>Step 3: assemble the reduced equations and force ports</summary>

Let $F_{ax}=k_{ax}(x_d-x_s)+c_{ax}(\dot x_d-\dot x_s)$. Then

$$m_d\ddot x_d=F_{mag}+F_{det}-c_m\dot x_d-F_{ax}-F_{f,n}-F_{f,r}-F_{f,d},$$

$$m_s\ddot x_s=F_{ax}+F_{f,n}-F_{f,g}.$$

In matrix form for the frictionless linear baseline,

$$
\underbrace{\begin{bmatrix}m_d&0\\0&m_s\end{bmatrix}}_{\mathbf M_r}\ddot{\mathbf x}
+\underbrace{\left[c_{ax}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}
+c_m\begin{bmatrix}1&0\\0&0\end{bmatrix}\right]}_{\mathbf C_r}\dot{\mathbf x}
+\underbrace{\begin{bmatrix}K_m+K_{det}+k_{ax}&-k_{ax}\\-k_{ax}&k_{ax}\end{bmatrix}}_{\mathbf K_r}\mathbf x
=\underbrace{\begin{bmatrix}K_m\\0\end{bmatrix}}_{\mathbf b_r}x_{cmd}.
$$

The nut microslip port is internal and equal/opposite. Guideway, gross nut rolling, and drivetrain ports act against ground in the reduced coordinates. That distinction is preserved through the reduction.

</details>

<details>
<summary>Step 4: analytical modal polynomial and transfer function</summary>

Ignoring damping for the closed-form modal calculation,

$$\det(\mathbf K_r-\omega^2\mathbf M_r)=0,$$

which expands to

$$m_dm_s\omega^4-\left[m_dk_{ax}+m_s(K_m+K_{det}+k_{ax})\right]\omega^2+(K_m+K_{det})k_{ax}=0.$$

The command-to-stage transfer function used for the Bode plot is

$$G(s)=\frac{X_s(s)}{X_{cmd}(s)}
=\mathbf e_2^T\left(\mathbf M_rs^2+\mathbf C_rs+\mathbf K_r\right)^{-1}\mathbf b_r.$$

The full ten-DOF response uses the identical dynamic-stiffness expression with $\mathbf e_{10}$ and the full matrices. No fitted transfer-function numerator is introduced.

</details>

## 7. Full-versus-reduced verification

![Full versus reduced Bode, bounded stepping, and reduction residual](rendered_assets/full_vs_reduced_verification.svg)

The comparison is deliberately linear and frictionless so that it tests structural reduction rather than confounding it with different friction memories. The same zero-order-held closed command sequence drives both models: 0 → +78.125 nm → 0 → −78.125 nm → 0. The sequence uses the conservative 64-microstep STEP/DIR quantum and ends at its starting level.

The full model includes the discarded internal resonances. Agreement is expected only up to 900 Hz. This is a reduction check, not independent modal validation. The same measured 690 Hz feature set $k_{ax}$ and is then reproduced by the reduced model. The effective mass and mode identity remain alternative explanations if the compliance budget does not close. See [Appendix A](#appendix-a-position-dependent-axial-stiffness) for the carriage-position stiffness sweep.

## 8. Friction constitutive laws

At each active site the Stribeck level is

$$s(v)=F_c+(F_s-F_c)\exp\left[-\left|v/v_s\right|^\delta\right].$$

<details>
<summary>LuGre derivation used in A, B, and C</summary>

The average bristle displacement $z$ evolves as

$$\dot z=v-\sigma_0\frac{|v|}{s(v)}z,$$

and the friction output is

$$F_f=\sigma_0z+\sigma_1\dot z+\sigma_2v.$$

Near rest, $\dot z\approx v$ and $F_f\approx\sigma_0z+(\sigma_1+\sigma_2)v$, so the frequency-response linearization adds $\sigma_0$ stiffness and $\sigma_1+\sigma_2$ damping along the site’s velocity vector.

</details>

<details>
<summary>GMS derivation and corrected slip-attractor sign used in A2, B2, and C2</summary>

Four parallel elements carry force states $F_i$. While an element sticks,

$$\dot F_i=k_iv,\qquad |F_i|<\nu_i s(v).$$

Once it yields in the current direction, the stable slip branch is

$$\boxed{\dot F_i=C\left[\operatorname{sgn}(v)-\frac{F_i}{\nu_i s(v)}\right]}.$$

Its equilibrium is $F_i=\operatorname{sgn}(v)\nu_i s(v)$ and the derivative points back toward that equilibrium on either velocity sign. Writing an extra $\operatorname{sgn}(v)$ outside the entire bracket makes the negative-velocity branch repelling; that was the sign defect corrected in the Revision 2 work and is not repeated here.

The output is

$$F_f=\sum_{i=1}^{4}F_i+\sigma_2v.$$

The element force fractions and stiffnesses are normalized so that

$$\boxed{\sum_{i=1}^{4}\nu_i=1},\qquad
\boxed{\sum_{i=1}^{4}k_i=\sigma_0}.$$

The builder verifies both identities before any simulation or rendering starts. It aborts rather than silently renormalizing an inconsistent parameter set. Therefore the element breakaway forces sum to the site Stribeck force and the stuck elements sum to the specified aggregate presliding stiffness.

On reversal the affected element returns to the stuck branch. Different $k_i$ and $\nu_i$ create different yield distances and preserve non-local presliding memory.

<details>
<summary>Exact re-stick test and derivative-evaluation ordering</summary>

Each GMS call is evaluated from the **current Runge–Kutta trial state** in this order:

1. Read the current site velocity $v$ and element forces $F_i$; compute $s(v)$, the thresholds $\nu_i s(v)$, and $k_i$.
2. If $|v|\le10^{-14}$ m/s, hold every element state with $\dot F_i=0$. No branch transition is inferred from a zero-velocity sign.
3. Otherwise evaluate the reversal/re-stick predicate **before assigning a derivative**: $vF_i\le0$. If true, select the stuck derivative $\dot F_i=k_iv$.
4. If it is not a reversal, test the current-state yield condition $|F_i|<\nu_i s(v)$. A sub-threshold element also receives $\dot F_i=k_iv$.
5. Only when neither test is true is the stable slip-attractor derivative evaluated.
6. Compute the friction output from the unadvanced trial-state forces, $F_f=\sum_iF_i+\sigma_2v$, and return all derivatives to RK4. RK4 then forms its next trial state and repeats every test.

Thus a derivative never selects its own branch during the same evaluation. There is no event localization at the threshold crossing. Branch switching uses the RK trial grid, so Section 13 includes a time-step check.

</details>

</details>

### 8.1 Executed provisional friction values

<details class="parameter-group">
<summary>Friction and GMS entry parameters</summary>

| Site | $\sigma_0$ | $\sigma_1$ | $\sigma_2$ | $F_s$ | $F_c$ | $v_s$ | GMS $C$ (N/s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Guideway | [[input:g_sigma0=7.600e5]] | [[assumed:g_sigma1=3.0]] | [[assumed:g_sigma2=0.40]] | [[assumed:g_Fs=3.0]] | [[assumed:g_Fc=2.4]] | [[assumed:g_vs=2.5e-4]] | [[assumed:g_C=5.000e3]] |
| Nut microslip $n$ | [[assumed:n_sigma0=2.000e6]] | [[assumed:n_sigma1=5.0]] | [[assumed:n_sigma2=0.25]] | [[assumed:n_Fs=1.6]] | [[assumed:n_Fc=1.2]] | [[assumed:n_vs=2.0e-4]] | [[assumed:n_C=5.000e3]] |
| Nut rolling $r$ | [[assumed:r_sigma0=2.000e6]] | [[assumed:r_sigma1=5.0]] | [[assumed:r_sigma2=0.25]] | [[assumed:r_Fs=5.0]] | [[assumed:r_Fc=4.0]] | [[assumed:r_vs=2.0e-4]] | [[assumed:r_C=5.000e3]] |
| Drivetrain $d$ | [[assumed:d_sigma0=1.000e6]] | [[assumed:d_sigma1=4.0]] | [[assumed:d_sigma2=0.20]] | [[assumed:d_Fs=2.0]] | [[assumed:d_Fc=1.5]] | [[assumed:d_vs=3.0e-4]] | [[assumed:d_C=5.000e3]] |

The four executed GMS elements use shared force fractions $\nu_i$ and site-scaled stiffnesses $k_i$:

| Element $i$ | Force fraction $\nu_i$ | Guideway $k_{i,g}$ | Nut $k_{i,n}=k_{i,r}$ | Drivetrain $k_{i,d}$ |
|---:|---:|---:|---:|---:|
| 1 | [[assumed:gms_nu1=0.10]] | [[assumed:g_k1=3.040e5]] | [[assumed:n_k1=8.000e5]] | [[assumed:d_k1=4.000e5]] |
| 2 | [[assumed:gms_nu2=0.20]] | [[assumed:g_k2=2.280e5]] | [[assumed:n_k2=6.000e5]] | [[assumed:d_k2=3.000e5]] |
| 3 | [[assumed:gms_nu3=0.30]] | [[assumed:g_k3=1.520e5]] | [[assumed:n_k3=4.000e5]] | [[assumed:d_k3=2.000e5]] |
| 4 | [[assumed:gms_nu4=0.40]] | [[assumed:g_k4=7.600e4]] | [[assumed:n_k4=2.000e5]] | [[assumed:d_k4=1.000e5]] |
| **Executed sum** | **1.00** | **$7.600\times10^5$** | **$2.000\times10^6$** | **$1.000\times10^6$** |

These values support comparison, but still require identification.

</details>

### 8.2 Force locations

Define the reduced velocity vector $\dot{\mathbf x}=[\dot x_d,\dot x_s]^T$. Each friction site is a power-conjugate port:

| Site | Velocity row $\mathbf H_\alpha$ | Driving velocity $v_\alpha=\mathbf H_\alpha\dot{\mathbf x}$ | Applied generalized force $-\mathbf H_\alpha^TF_{f,\alpha}$ |
|---|---|---|---|
| Guideway $g$ | $[0,1]$ | $\dot x_s$ | $[0,-F_{f,g}]^T$ |
| Nut microslip $n$ | $[1,-1]$ | $\dot x_d-\dot x_s$ | $[-F_{f,n},+F_{f,n}]^T$ |
| Nut rolling $r$ | $[1,0]$ | $\dot x_d$ | $[-F_{f,r},0]^T$ |
| Drivetrain $d$ | $[1,0]$ | $\dot x_d$ | $[-F_{f,d},0]^T$ |

The minus-transpose rule guarantees dissipated power $\dot{\mathbf x}^T(-\mathbf H^TF_f)=-vF_f\le0$ when the constitutive force opposes motion.

The presliding tangent alone cannot separate $k_{ax}$ from $\sigma_{0,n}$ because both enter the differential stiffness as a sum. Separation requires finite-amplitude reversal data. Microslip yields and dissipates; $k_{ax}$ remains conservative.

[Appendix B](#appendix-b-reduced-model-bond-graph) draws the same rows as power bonds. The paired nut bonds make the equal-opposite force application visible.

<details>
<summary>Nonlinear time-domain implementation</summary>

At every Runge–Kutta evaluation the model performs the following operations.

1. Compute $v_g=\dot x_s$, $v_n=\dot x_d-\dot x_s$, and $v_r=v_d=\dot x_d$ from the current mechanical velocities.
2. For each active site, advance either one LuGre state $z_\alpha$ or four GMS force states $F_{i,\alpha}$.
3. Evaluate the site force from that state and velocity.
4. Apply guideway friction only to the stage, nut microslip equal-and-opposite across the two bodies, and both rolling and drivetrain drag to the drive body.
5. Integrate the friction states together with $x_d,x_s,\dot x_d,\dot x_s$; the memory is not evaluated afterward as a plotting correction.

Cases A/A2 activate $d,g$. Cases B/B2 activate $d,r,n$. Cases C/C2 activate all four sites. Case 0 has no friction state. All friction parameters remain provisional until identified, but no defined gross-drag port is silently disabled.

</details>

<details>
<summary>Linear Bode implementation versus nonlinear stepping implementation</summary>

A nonlinear hysteretic law has no single amplitude-independent Bode response. The displayed Bode curves use the zero-velocity presliding tangent. For each active site,

$$\Delta\mathbf K=\sigma_0\mathbf H^T\mathbf H.$$

LuGre adds tangent damping $(\sigma_1+\sigma_2)\mathbf H^T\mathbf H$; the present GMS tangent adds $\sigma_2\mathbf H^T\mathbf H$, while its four elastic states supply the presliding stiffness. The time-domain plots use the complete nonlinear state equations, including Stribeck variation, yielding, and reversal memory.

</details>

<details>
<summary>What changes in the simulated outcome?</summary>

- Guideway friction is ground-referenced. Its presliding stiffness changes the local tangent and adds reversal hysteresis.
- Nut microslip is an internal equal-and-opposite port. It changes differential deformation and the relative mode.
- Gross nut rolling drag $F_{f,r}$ acts on $v_d$ and remains active when $v_n$ approaches zero.
- Aggregated drivetrain drag $F_{f,d}$ is active in every friction case.
- LuGre and GMS can share the same small-signal stiffness and therefore similar modal frequencies, while producing different reversal loops and final errors. GMS retains non-local memory through separately yielding elements; LuGre has one local bristle state.
- Friction can add damping at some amplitudes but can also produce stick–slip, lost motion, amplitude-dependent apparent stiffness, and nonzero final error. These effects cannot be inferred from the linear Bode curve alone.

</details>

## 9. Force-instrumented partial-slip memory experiment

The main 64-microstep sequence does not revisit nested reversal points. This separate test uses cases A and A2. Force is the primary discriminator. Interferometer displacement is secondary because the earlier modeled LuGre/GMS difference was comparable to the 4.6 nm project ADEV floor.

![Presliding nested-reversal motion, modeled command-stage deviation, memory loops, and comparison metrics](rendered_assets/presliding_memory_comparison.svg)

<details>
<summary>9.1 Exact 64-microstep command</summary>

One external STEP/DIR quantum is

$$q_\mu=\frac{5\ \mu\mathrm m}{64}=78.125\ \mathrm{nm}.$$

After a 5 ms zero dwell, each listed level is held for 10 ms:

| Plateau | Command (microsteps) | Command (nm) | Purpose |
|---:|---:|---:|---|
| 1 | 0 | 0.00 | origin |
| 2 | +48 | +3750.00 | positive outer reversal |
| 3 | +12 | +937.50 | first inner return level |
| 4 | +42 | +3281.25 | nested reversal |
| 5 | +12 | +937.50 | revisit +12 |
| 6 | +48 | +3750.00 | revisit +48 |
| 7 | 0 | 0.00 | close positive branch |
| 8 | -46 | -3593.75 | negative outer reversal |
| 9 | -12 | -937.50 | second inner return level |
| 10 | -40 | -3125.00 | nested reversal |
| 11 | -12 | -937.50 | revisit -12 |
| 12 | -46 | -3593.75 | revisit -46 |
| 13 | 0 | 0.00 | final positive step back to the origin |

The largest increment is 48 microsteps = 3.7500 µm. This dedicated identification test exceeds the quarter-step linearization bound, so it uses the nonlinear magnetic law. It remains below one full step and ends at its starting command.

</details>

<details>
<summary>9.2 Why this remains presliding while still activating GMS memory</summary>

For the first guideway GMS element, the zero-speed yield displacement predicted by the provisional parameters is

$$z_{y,1}=\frac{\nu_1F_s}{k_1}
=\frac{0.10(3.0)}{0.40(7.60\times10^5)}
=0.987\ \mu\mathrm m.$$

The second guideway element yields at

$$z_{y,2}=\frac{0.20(3.0)}{0.30(7.60\times10^5)}=2.63\ \mu\mathrm m.$$

The 3.7500 µm outer command crosses the nominal yield distances of two elements. This gives more distributed memory than the earlier 1.094 µm test, which reached only the first element. The generated force audit checks whether the aggregate interface remains below gross breakaway.

This distinction matters. If every element stayed perfectly elastic, both laws would reduce almost to a spring and their loops would be indistinguishable. If every element entered gross sliding, the nested presliding memory would be erased. The chosen amplitude lies between those two uninformative limits for the current provisional parameters.

</details>

<details>
<summary>9.3 Nonlocal-memory mechanism: one LuGre state versus four GMS states</summary>

LuGre compresses the guideway interface into one average bristle state $z_g$. At a given current $z_g$ and velocity it has no independent record of several earlier reversal points. It produces a local hysteresis loop, but nested minor-loop closure is not an independently stored property.

GMS carries four element-force states $F_{1,g},\ldots,F_{4,g}$ with different stiffnesses and yield thresholds. A reversal can unload one element while another remains on a different branch. The vector of retained states therefore depends on more than the latest displacement and preserves the order of prior extrema. This is the nonlocal memory being exercised when +12, +48, -12, and -46 microsteps are revisited.

The plotted force-position loops use the friction forces produced inside the time integration. They are not reconstructed from position afterward. Faint lines show the full dynamic trace; markers show the mean over the final 2 ms of each plateau. The experiment requires a force measurement across the compliant path. A displacement-only result is insufficient unless the measured LuGre/GMS difference clears the metrology floor with margin.

</details>

<details>
<summary>9.4 Metrics, equations, and interpretation</summary>

Let $\bar d_j$ and $\bar F_j$ be the mean modeled command-stage deviation and guideway friction force over the final 2 ms of plateau $j$. For the repeated-level pair set

$$\mathcal P=\{(2,6),(3,5),(8,12),(9,11)\},$$

where the plateau numbers are one-based as in the table, define

$$E_{ret}=\frac{1}{|\mathcal P|}\sum_{(i,j)\in\mathcal P}|\bar d_i-\bar d_j|,$$

$$F_{ret}=\frac{1}{|\mathcal P|}\sum_{(i,j)\in\mathcal P}|\bar F_i-\bar F_j|.$$

$E_{ret}$ measures how closely the modeled plant response returns to the same result at a repeated command level; $F_{ret}$ directly measures constitutive return-point closure. The final-origin metric is $|\bar d_{13}|$. Whole-sequence RMS is also reported, but it is dominated by the commanded jumps and is less sensitive to hysteretic memory.

<!-- BEGIN GENERATED PRESLIDING SUMMARY -->
| Executed metric | LuGre A | GMS A2 | GMS minus LuGre |
|---|---:|---:|---:|
| Whole-sequence RMS command-stage deviation | 2636.84 nm | 2677.65 nm | +40.81 nm |
| Peak absolute command-stage deviation | 5975.60 nm | 5953.19 nm | -22.41 nm |
| Mean repeated-return deviation mismatch | 1540.21 nm | 1621.63 nm | +81.42 nm |
| Mean repeated-return friction-force mismatch | 0.8157 N | 0.7157 N | -0.1001 N |
| Absolute mean error after final return to zero | 2223.86 nm | 2306.85 nm | +82.98 nm |

The maximum executed guideway friction magnitude is **2.270 N**, or **75.7%** of the provisional 3.0 N macro breakaway level. The sequence therefore probes partial slip rather than gross sliding.

The whole-sequence RMS includes the unavoidable error at every instantaneous command edge. The repeated-return and final-origin measures isolate the history dependence that this experiment is intended to distinguish. The provisional parameters do not guarantee that GMS closes more tightly than LuGre; measured force loops must decide that question.
<!-- END GENERATED PRESLIDING SUMMARY -->

The comparison does not assume that GMS is better. With the current provisional parameters, GMS reduces the force return mismatch by about 0.10 N while its displacement return and final-origin metrics are worse. This mixed result is the correct outcome of an identification design: measured force loops must select and fit the constitutive law.

</details>

## 10. Response comparison across friction cases

All seven cases use the same mechanical plant. The table shows the active force placement.

| Cases | Active port | Generalized force |
|---|---|---|
| 0 | none | $[F_{mag}+F_{det}-c_m\dot x_d,\ 0]^T$ |
| A, A2 | drivetrain + guideway | $[F_{mag}+F_{det}-c_m\dot x_d-F_{f,d},\ -F_{f,g}]^T$ |
| B, B2 | drivetrain + nut rolling + microslip | $[F_{mag}+F_{det}-c_m\dot x_d-F_{f,d}-F_{f,r}-F_{f,n},\ +F_{f,n}]^T$ |
| C, C2 | all four | $[F_{mag}+F_{det}-c_m\dot x_d-F_{f,d}-F_{f,r}-F_{f,n},\ F_{f,n}-F_{f,g}]^T$ |

![All-case Bode overlay, resonance zoom, and matched-law magnitude differences](rendered_assets/lugre_gms_pairwise_comparison.svg)

<!-- BEGIN GENERATED BODE COMPARISON -->
| Topology | Local peak | Shift from Case 0 | Largest GMS/LuGre gap | Cause |
|---|---:|---:|---:|---|
| Case 0 | 695.5 Hz | reference | not applicable | No friction tangent |
| A/A2 | 718.1 Hz | +22.6 Hz, +3.3% | 0.45 dB at 718 Hz | Guideway presliding stiffness acts against ground |
| B/B2 | 753.4 Hz | +57.9 Hz, +8.3% | 0.74 dB at 753 Hz | Nut microslip shifts the relative mode; rolling and drivetrain tangents act on the drive |
| C/C2 | 774.8 Hz | +79.3 Hz, +11.4% | 1.16 dB at 775 Hz | All four friction tangents are active |
<!-- END GENERATED BODE COMPARISON -->

Each matched LuGre/GMS pair has the same presliding stiffness, so its resonance frequency is the same. LuGre adds $\sigma_1+\sigma_2$ tangent damping. The current GMS tangent adds only $\sigma_2$. This damping difference changes peak height, most visibly in C/C2.

## 11. Generated numerical summary

<!-- BEGIN GENERATED RESPONSE SUMMARY -->
| Case | Friction law | Presliding modes (Hz) | Presliding tangent gain $X_s/X_{cmd}$ | Smallest first-yield travel | First-step overshoot | Full-sequence RMS deviation | Peak absolute deviation | Final-window RMS deviation |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | none | 180.6, 695.6 | 0.75000 | not applicable | 48.3% | 32.9 nm | 78.1 nm | 20.5 nm |
| A | LuGre | 181.6, 718.1 | 0.69558 | 0.500 µm | 37.4% | 32.3 nm | 78.1 nm | 19.8 nm |
| A2 | GMS | 181.6, 718.1 | 0.69558 | 0.500 µm | 37.3% | 32.4 nm | 78.1 nm | 19.9 nm |
| B | LuGre | 182.3, 754.1 | 0.73602 | 0.200 µm | 44.4% | 32.6 nm | 78.1 nm | 21.1 nm |
| B2 | GMS | 182.3, 754.1 | 0.73602 | 0.200 µm | 44.4% | 32.7 nm | 78.1 nm | 21.1 nm |
| C | LuGre | 182.8, 775.0 | 0.69341 | 0.200 µm | 35.4% | 32.2 nm | 78.1 nm | 19.9 nm |
| C2 | GMS | 182.8, 775.0 | 0.69341 | 0.200 µm | 35.4% | 32.2 nm | 78.1 nm | 20.0 nm |

The tangent gain is a local presliding linearization. It is valid only below the listed first-yield travel and is not a full-range tracking gain. Sustained travel produces bounded friction offsets in the nonlinear model. The three deviation columns use $d(t)=x_{cmd}(t)-x_s(t)$. They describe the open-loop modeled plant response under each friction law, not closed-loop servo tracking performance. The final column summarizes the last 2 ms of the nonlinear run and is not an identified settling specification. All cases include rated-current commutation, enabled detent torque, and the highlighted electromagnetic damping assumption. Case 0 remains frictionless.

### Generated reduction audit

| Quantity | Executed value |
|---|---:|
| Closure-derived $k_{ball}$ | 43.871 MN/m |
| Motor rotor inertia | 9.000e-07 kg m² |
| Coupling inertia | 1.180e-06 kg m² |
| 0.320 m screw inertia | 1.010e-06 kg m² |
| 0.320 m screw mass | 0.1263 kg |
| Full-model reflected drivetrain mass | 121.994 kg |
| Rated-current holding torque | 0.060 N m |
| Enabled detent torque | 0.005 N m |
| Full/reduced sequence RMS residual | 5.666 nm |
| Full/reduced sequence peak residual | 17.304 nm |

The reduced drive mass is derived from the listed component inertias and the current lead. It is not an independent input.
<!-- END GENERATED RESPONSE SUMMARY -->

## 12. Interpreting commanded and actual motion

The plotted difference is defined as

$$d_{model}(t)=x_{cmd}(t)-x_s(t).$$

The reported metrics are

$$d_{RMS}=\sqrt{\frac{1}{T}\int_0^T d_{model}^2(t)\,dt},\qquad
d_{max}=\max_{0\le t\le T}|d_{model}(t)|.$$

| Metric | Window | Interpretation |
|---|---|---|
| Whole-sequence RMS | full 85 ms | Includes every command edge |
| Peak absolute deviation | full 85 ms | Usually occurs at a command edge |
| Final-window RMS | final 2 ms | Describes the last zero-command dwell |

These are open-loop model descriptors, not tracking specifications. The model has no position controller, estimator, sensor dynamics, or shaped trajectory. The damping term $c_m$ removes the earlier unphysical sustained ringing. Remaining values are provisional until $\zeta_m$, $c_{ax}$, and the friction parameters are identified.

### 12.1 Lead accuracy and offline compensation

Lead accuracy can dominate the full-range error even when local friction is well compensated. Karl Hipp lists the following tolerance for useful travel up to 315 mm:

| Accuracy class | Lead tolerance | Share of a 10 µm budget |
|---|---:|---:|
| IT1 | 6 µm | 60% |
| IT3 | 12 µm | 120% |

The source table assigns 6 µm to IT1, not IT3. See [Karl Hipp ballscrew technology](https://www.karl-hipp.de/en/technology). The installed screw class is not recorded in this project, so neither value is executed as an error waveform. A mapped lead-error curve should be a primary offline pre-distortion input. The class tolerance is a bound, not the measured correction table.

## 13. Verification checks and limitations

<details>
<summary>Checks performed by construction</summary>

1. $\mathbf M$ is diagonal and positive for all executed parameters.
2. Every passive spring and damper is added by a positive-semidefinite outer product.
3. The nut virtual-work vector applies $+rF_n$, $+F_n$, and $-F_n$ with consistent power.
4. The GMS negative-velocity slip equilibrium is attracting.
5. The nonlinear command is held constant over all four RK4 stages at a discontinuity.
6. The main response uses 78.125 nm increments. The separate memory-identification test intentionally reaches 3.7500 µm and uses the nonlinear magnetic law.
7. Full and reduced verification use the same command, sample grid, and damping repair.
8. The generated metrics table is rewritten by the builder, tying numbers to executed code.
9. The builder asserts $\sum_i\nu_i=1$ and $\sum_i k_i=\sigma_0$ for every defined GMS site before simulation.

</details>

### 13.1 GMS step-halving convergence

The production nonlinear plots use fixed-step RK4 with $h=5$ µs. To test sensitivity of the requested final-window RMS result, the builder reruns A2, B2, and C2 using $h=10$, 5, and 2.5 µs. All command transitions fall exactly on all three grids, and the command remains one zero-order-held value across the four RK stages of each step.

<!-- BEGIN GENERATED STEP HALVING SUMMARY -->
| Case | 10.0 us | 5.0 us | 2.5 us | $\Delta R_{10\to5}$ | $\Delta R_{5\to2.5}$ | Difference ratio |
|---|---:|---:|---:|---:|---:|---:|
| A2 | 19.86722 nm | 19.88457 nm | 19.89325 nm | 0.01735 nm | 0.00868 nm | 2.00 |
| B2 | 21.13011 nm | 21.14869 nm | 21.15798 nm | 0.01858 nm | 0.00930 nm | 2.00 |
| C2 | 20.00969 nm | 20.02846 nm | 20.03786 nm | 0.01877 nm | 0.00939 nm | 2.00 |

The successive change decreases for all three GMS cases, which is consistent with time-step convergence for this reported metric. The largest 5.0-to-2.5 us relative change is **0.0469%**.

These values use the identical 85 ms zero-order-held command and the identical final 2 ms RMS definition. Since GMS branch switching is evaluated at RK trial states without event localization, the difference ratio is a sensitivity indicator, not a claimed fourth-order convergence rate for the hybrid trajectory.
<!-- END GENERATED STEP HALVING SUMMARY -->

<details>
<summary>Known limitations and measurements that would remove assumptions</summary>

- Coupling inertia and torsional stiffness require CAD or datasheet values.
- Bearing stiffness/contact angle and preload require BOM confirmation or static loading.
- $k_{ball}$ is a closure-derived remainder, not a direct Hertzian calculation or measurement.
- Driver mode and effective damping must be identified; $\zeta_m=0.05$ is provisional and is accompanied by a 0.02 to 0.50 sensitivity sweep.
- LuGre and GMS values require velocity sweeps and nested reversal tests.
- The installed screw accuracy class and measured lead-error map are missing. This is a first-order full-range uncertainty, not a minor residual.
- Yaw, pitch, roll, rail bending, cyclic error, runout, temperature, and load-dependent nut friction are omitted.
- The electrical winding/current-controller dynamics are represented only by effective stiffness and damping.
- Editing inputs in the rendered HTML does not recompute the static plots; the browser cannot safely execute the local Python model.

</details>

## Appendix A. Position-dependent axial stiffness

![Axial stiffness and stage-mode prediction versus nut position](rendered_assets/position_dependence.svg)

For the screw segment before the nut, $k_{sha}=EA/L_{free}$. A longer free length reduces both $k_{sha}$ and the series stiffness $k_{ax}$. Test this prediction at carriage positions of 50, 150, and 250 mm.

## Appendix B. Reduced-model bond graph

![Reduced-model bond graph and power-port audit](rendered_assets/reduced_bond_graph.svg)

The two 1-junctions carry $\dot x_d$ and $\dot x_s$. The central 0-junction carries the common internal force. Structural compliance, damping, and nut microslip are distinct parallel constitutive elements. Gross nut rolling and drivetrain drag connect to the drive junction. The bond directions reproduce $\mathbf Q_f=-\mathbf H^TF_f$ and $P_f=-v_fF_f\le0$.

This graph is the visual form of the Section 8.2 incidence rows. It adds no model elements.

## Appendix C. Critical-error disposition

| Item | Evaluation | Implemented disposition |
|---:|---|---|
| 1 | Confirmed | Split gross nut rolling drag onto $v_d$. Reduced differential $F_s$ to 1.6 N, giving a 0.20 µm first GMS yield. The tangent correlation with $k_{ax}$ is stated. |
| 2 | Confirmed | Execute $F_{f,d}$ in A/A2, B/B2, and C/C2. |
| 3 | Stale in the supplied review | Detent torque was already enabled. Added the 0.060 N·m motor bound of 266 nm with 1.25 µm period. |
| 4 | Confirmed | Made force the primary memory-test metric and increased the excursion to cross two nominal guideway yield distances. |
| 5 | Confirmed | Main runs now use 64 microsteps, 78.125 nm. The optional 256-interpolated quantum is 19.531 nm. Hardware configuration remains open. |
| 6 | Partly confirmed | The old 0.50 baseline was unsupported. Execute 0.05 and show 0.02 to 0.50 sensitivity. The driver mode must be measured. |
| 7 | Confirmed | Keep the failed compliance budget prominent. State that reproducing 694 Hz is calibration, not validation. |
| 8 | Confirmed | Rename DC gain to presliding tangent gain and report the first-yield validity travel. |
| 9 | Source claim corrected | The manufacturer gives 6 µm for IT1 and 12 µm for IT3 up to 315 mm. Lead mapping is a primary pre-distortion input. |
| 10 | Stale in the supplied review | Full equations already apply equal and opposite microslip reactions. The bond graph audits the signs. |
| 11 | Partly stale | $J_m$ and the 0.320 m screw inertia were already corrected. Set $J_c=1.18\times10^{-6}$ kg·m² and rebuild $m_d$. |

## Appendix D. Variable and parameter glossary

<details>
<summary>Expand all symbols, states, ports, and units</summary>

| Symbol | Definition | Unit |
|---|---|---|
| $L$ | screw lead | m/rev |
| $r=L/(2\pi)$ | screw transmission ratio | m/rad |
| $N_r$ | rotor tooth count | – |
| $\theta_m,\theta_c,\theta_{s1..s3}$ | full-model torsional coordinates | rad |
| $u_b,u_e,u_f,u_n,x_s$ | full-model axial coordinates | m |
| $x_d$ | reduced effective drive coordinate | m |
| $x_{cmd}$ | commanded field-equivalent linear position | m |
| $d_{model}=x_{cmd}-x_s$ | modeled open-loop command-stage deviation; not a servo tracking specification | m |
| $J_m,J_c,J_{s1..s3}$ | rotational inertias | kg·m² |
| $m_b,m_e,m_f,m_n,m_{stage}$ | axial lumped masses | kg |
| $m_d,m_s$ | reduced drive and stage masses | kg |
| $k_m$ | rotational magnetic stiffness | N·m/rad |
| $K_m=k_m/r^2$ | reflected magnetic stiffness | N/m |
| $T_{max},F_{max}$ | peak magnetic torque and reflected force | N·m, N |
| $\kappa=N_r/r$ | magnetic spatial wavenumber | rad/m |
| $k_{c1},k_{c2}$ | coupling torsional stiffnesses | N·m/rad |
| $k_{\theta a},k_{\theta b}$ | screw torsional stiffnesses | N·m/rad |
| $k_{brg}$ | support-bearing axial stiffness | N/m |
| $k_{sha},k_{shb}$ | screw axial stiffnesses before/after nut | N/m |
| $k_{ball}$ | normal ball-contact axial stiffness | N/m |
| $k_{mnt}$ | nut-mount axial stiffness | N/m |
| $k_{ax}$ | reduced series axial stiffness | N/m |
| $c_{ax}$ | reduced axial structural damping | N·s/m |
| $c_m,c_{\theta m}$ | reduced/rotational electromagnetic damping | N·s/m, N·m·s/rad |
| $u_t=u_e+r\theta_{s2}$ | full-model screw-transformer output | m |
| $\delta_n$ | ball-contact deformation $u_n-u_e-r\theta_{s2}$ | m |
| $F_n$ | normal axial ball-contact force | N |
| $F_{f,g}$ | guideway friction force | N |
| $F_{f,n}$ | axial-equivalent nut differential microslip force | N |
| $F_{f,r}$ | gross nut rolling-drag force on the drive coordinate | N |
| $F_{f,d}$ | aggregated drivetrain-to-ground friction | N |
| $s(v)$ | velocity-dependent Stribeck force level | N |
| $F_s,F_c$ | static and Coulomb force levels | N |
| $v_s,\delta$ | Stribeck velocity and exponent | m/s, – |
| $z$ | LuGre bristle deflection | m |
| $\sigma_0,\sigma_1,\sigma_2$ | LuGre stiffness, microdamping, viscous coefficients | N/m, N·s/m, N·s/m |
| $F_i,k_i,\nu_i$ | GMS element force, stiffness, threshold weight | N, N/m, – |
| $C$ | GMS slip-attraction rate coefficient | N/s |
| $\mathbf M,\mathbf C,\mathbf K$ | mass, damping, stiffness matrices | mixed native units |
| $\mathbf h_n$ | nut-deformation incidence/virtual-work vector | mixed |

</details>
