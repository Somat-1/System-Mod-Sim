# Ball-screw stage system modelling workspace

This workspace documents and executes a physics-based model of the stepper-driven ball-screw positioning stage. Revision 3 is the active model. Revision 2 is retained as the previous baseline and should not be edited when changing Revision 3.

## Start here

- [Revision 3 model specification](Rev%203/ball_screw_stage_dynamic_derivation_v3.html) — assumptions, coordinates, topology, parameters, reduction, and friction variants.
- [Revision 3 analytical derivation and responses](Rev%203/Analytical_derivation_and_responses_v3.html) — comprehensive equations, expandable derivation steps, focused Bode plots, step responses, modeled command-stage deviations, kinematic diagram, presliding memory test, and verification results.
- [Revision 3 Markdown derivation](Rev%203/Analytical_derivation_and_responses_v3.md) — editable source for the analytical document.
- [Single Revision 3 builder](Rev%203/build_model_documentation.py) — simulations, plots, generated metric tables, and Markdown-to-HTML rendering.
- `Rev 2/` — preserved earlier revision for comparison.

## Model derivation in brief

1. **Choose independent coordinates.** The full model retains five rotational coordinates, $	heta_m,\theta_c,\theta_{s1},\theta_{s2},\theta_{s3}$, and five axial coordinates, $u_b,u_e,u_f,u_n,x_s$. The command is an input and friction memories are internal constitutive states; neither changes the mechanical degree-of-freedom count.
2. **Define elastic deflections.** Each coupling, screw segment, bearing, ball contact, and nut mount is represented by a relative deformation. These deflections define the total potential energy; the inertias and masses define the kinetic energy.
3. **Assemble the equations.** Lagrange's equation produces the full system in the form

   $$\mathbf M\ddot{\mathbf q}+\mathbf C\dot{\mathbf q}+\mathbf K\mathbf q
   =\mathbf f_{mag}+\mathbf f_{friction}.$$

   Passive two-node elements are assembled as positive-semidefinite outer products, which preserves sign and energy consistency.
4. **Apply ball-screw kinematics.** With $r=L/(2\pi)$, the torsional/axial nut-contact deformation couples screw rotation to translation. Virtual work fixes the equal-and-opposite force signs at the interface.
5. **Reduce only after auditing the full model.** Internal coordinates whose modes lie outside the intended bandwidth are condensed into a two-coordinate model: effective drive displacement $x_d$ and stage displacement $x_s$. The stage still has one externally observable translation; the second retained coordinate represents compliant internal drive motion.
6. **Retain the nonlinear stepper law.** The nonlinear simulations use the bounded magnetic force $F_{max}\sin[\kappa(x_{cmd}-x_d)]$ plus provisional electromagnetic damping. Linear Bode plots use its small-signal tangent.
7. **Apply friction through power-conjugate ports.** Guideway friction uses $v_g=\dot x_s$; nut friction uses $v_n=\dot x_d-\dot x_s$. Nut friction is therefore equal-and-opposite internally, while guideway friction acts against ground.

## Executed cases

| Case | Active friction site(s) | Constitutive law |
|---|---|---|
| 0 | none | frictionless modal baseline |
| A / A2 | guideway | LuGre / GMS |
| B / B2 | nut differential | LuGre / GMS |
| C / C2 | guideway and nut | LuGre / GMS |

LuGre uses one average bristle state per active site. GMS uses four force states with different thresholds, allowing nested return-point memory. The builder verifies $\sum_i\nu_i=1$ and $\sum_i k_i=\sigma_0$ before running any case.

## Simulations and numerical checks

- Frequency responses use the linearized presliding tangent and are shown per case, with topology-matched LuGre/GMS comparisons.
- The standard nonlinear command is bounded to one quarter of the 5 µm full-step pitch and ends at its starting level.
- A separate 1/32-full-step nested-reversal sequence exercises partial-slip memory and reports both ordinary tracking RMS and return-point closure metrics.
- Fixed-step RK4 holds a discontinuous command constant over all four stages. GMS re-stick and yield tests use the current RK trial state before the derivative is assigned.
- A generated step-halving study compares the final-window RMS error for A2, B2, and C2 at 10, 5, and 2.5 µs.
- Full-versus-reduced residuals, parameter-closure checks, modal comparisons, and command-bound checks are documented in the analytical HTML.

## Rebuilding Revision 3

From the workspace root, run:

```powershell
python ".\Rev 3\build_model_documentation.py"
```

The builder requires NumPy and Matplotlib. It rewrites the generated result blocks in the analytical Markdown, regenerates SVG files in `Rev 3/rendered_assets/`, and renders both Revision 3 Markdown documents to HTML. Edit the Markdown and the single builder; do not manually edit generated HTML or SVG output.

## Editable HTML parameters

Amber fields are provisional or assumed values. Browser edits persist in browser storage and the page URL, and **Save HTML copy** embeds the chosen values in a saved HTML file. The live transfer-function panel recalculates in the browser. Publication SVG plots and generated numerical tables remain build-time results and require rerunning the Python builder after changing source parameters.

## Interpretation boundary

The present numerical values combine measured inputs, closure-derived quantities, and highlighted assumptions. The model is suitable for transparent hypothesis testing and experiment design, but friction coefficients, electromagnetic damping, contact stiffnesses, and several inertias still require identification before the simulated amplitudes should be treated as validated hardware predictions.
