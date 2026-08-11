# Implementation issue backlog

This backlog tracks defects, rendering problems, build failures, and behavior that differs from the intended implementation. It does not track model features, scientific scope, or new script functionality.

## Status key

| Status | Meaning |
|---|---|
| Open | Reproduced and still needs a fix |
| Investigating | Cause or fix is being tested |
| Resolved | Fixed and verified |
| Deferred | Valid issue, but intentionally postponed |

## Active issues

| ID | Status | Area | Observed behavior | Expected behavior | Next check |
|---|---|---|---|---|---|
| IMP-005 | Open | Offline HTML | Equations depend on MathJax from a CDN. They do not typeset without network access unless the library is already cached. | A saved report should state this dependency or bundle the renderer. | Decide whether to bundle MathJax or add an explicit offline notice. |
| IMP-006 | Open | Generated artifacts | A rebuild changes SVG metadata and internal IDs even when a figure's data and layout are unchanged. | Unchanged figures should produce stable files and small review diffs. | Set deterministic SVG metadata and a fixed Matplotlib SVG hash salt. |
| IMP-021 | Open | Hardware configuration | The installed TMC2209 mode and effective external/interpolated microstep configuration are not recorded in the workspace. | Simulations should read the executed driver configuration rather than rely on the conservative 64-step assumption. | Capture UART registers or the standalone pin state and record the motion-controller STEP setting. |

## Resolved issues

| ID | Status | Area | Observed behavior | Resolution | Verified in |
|---|---|---|---|---|---|
| IMP-001 | Resolved | Report layout | Main sections blended together in the rendered HTML. | Added stronger section blocks, spacing, and distinct appendix styling. | Revision 3 rebuild |
| IMP-002 | Resolved | Parameter UI | Entry tables occupied the page even when parameters were not being edited. | Moved parameter groups into collapsed toggle panels. | Revision 3 rebuild |
| IMP-003 | Resolved | Response comparison | Seven similar case figures obscured the useful Bode differences. | Put each case response in an HTML toggle and retained one labeled pairwise comparison at the end. | Revision 3 settled-response rebuild |
| IMP-004 | Resolved | Technical writing | Long passages and em dashes made the analytical report harder to scan. | Shortened key explanations, replaced em dashes, and moved reference material to appendices. | Revision 3 rebuild |
| IMP-007 | Resolved | Parameter dependencies | Derived values such as $r$, $m_d$, and $K_m$ were editable independent inputs and could become inconsistent. | Made component values the inputs and recalculated dependent outputs and inline reduced-model results in the browser. | Revision 3 component-value rebuild |
| IMP-008 | Resolved | Drive model | Detent was either disabled or later promoted to a global spring, producing a false 0.75 Case 0 DC gain. | Retained the 0.005 N·m periodic nonlinear torque, removed detent from global stiffness, and report only a declared local 137–194 Hz tangent band with a 5 µm period. | Revision 3 settled-response rebuild |
| IMP-009 | Resolved | Generated SVGs | Matplotlib emitted trailing spaces that caused repository whitespace checks to fail. | Normalized SVG lines in the shared writer after every figure save. | Revision 3 component-value rebuild |
| IMP-010 | Resolved | Kinematic diagram | The axial load path incorrectly placed $u_f$ between $u_e$ and $u_n$, and $\theta_{s3}$ appeared inline. | Redrew both coordinates as dead-end overhang stubs and aligned the remaining coordinates by physical station. | Revision 3 topology-diagram rebuild |
| IMP-011 | Resolved | Nut friction port | The full equations applied $T_{f,n}$ only to $\theta_{s2}$, while the reduced equations used an internal equal-opposite port. | Defined $v_n=r\dot\theta_{s2}+\dot u_e-\dot u_n$, added the paired axial forces, and documented $T_{f,n}=rF_{f,n}$. | Revision 3 topology-diagram rebuild |
| IMP-012 | Resolved | Reduction diagram | Mass destination and compliance retention shared one fill color, which incorrectly mapped $u_n$ and implied that axial screw masses survived. | Made box fill encode mass aggregation and spring stroke encode compliance retention. | Revision 3 reduction-map rebuild |
| IMP-013 | Resolved | Diagram registration | Station headers, the nut transformer, and axial coordinates did not share a column grid; several numeric labels collided with elements. | Added station guides, aligned the nut column, and moved reduction evidence into a dedicated annotation column. | Revision 3 reduction-map rebuild |
| IMP-014 | Resolved | Friction topology | Gross rolling and drivetrain laws were both placed on the identical $[1,0]$ incidence row and could not be separated. | Combined their provisional force/tangent budget into one identifiable drive-side law; retained the lower-force differential microslip port with a 0.20 µm first yield. | Revision 3 settled-response rebuild |
| IMP-015 | Resolved | Executed cases | $F_{f,d}$ was defined but disabled in every case, so modeled drivetrain loss was silently zero. | Activated $F_{f,d}$ in A/A2, B/B2, and C/C2 and updated the diagrams and case map. | Revision 3 critical-error rebuild |
| IMP-016 | Resolved | Command generation | The later 78.125 nm main run never crossed yield and made all friction cases nearly identical. | Kept 64x quanta but span 0.234–2.031 µm absolute levels, cap adjacent increments at 1.016 µm, and finish with a positive return to zero. | Revision 3 settled-response rebuild |
| IMP-017 | Resolved | Response table | The zero-speed linearization was labeled "DC gain", which could be mistaken for full-range tracking gain. | Renamed it presliding tangent gain, added first-yield validity travel, and stated that nonlinear travel produces bounded offsets. | Revision 3 critical-error rebuild |
| IMP-018 | Resolved | Memory experiment | A/A2 alone did not test the exact $k_{ax}$/$\sigma_{0,n}$ correlation; a free-stage B/B2 run did not traverse nut yield; and 10 ms plateaus sampled ringing. | Added a dedicated blocked-stage, force-instrumented B/B2 identification loop, retained the normal free-stage plant for response plots, and derive at least 100 ms dwell from damping; endpoint metrics use settled 20 ms windows. | Revision 3 settled-response rebuild |
| IMP-019 | Resolved | Source interpretation | The installed screw class and travel geometry were not recorded correctly. | Record IT1, 192 mm complete screw length, approximately 170 mm usable distance, and 150 mm stage travel; remove the unnecessary Section 12.1 comparison. | Revision 3 settled-response rebuild |
| IMP-020 | Resolved | Parameter provenance | The coupling estimate remained 1.20×10⁻⁶ kg·m² after the annulus estimate had been revised. | Set $J_c=1.18×10^{-6}$ kg·m² and regenerated reflected mass and modal results. | Revision 3 critical-error rebuild |
| IMP-022 | Resolved | Reduction verification | Nanometre residuals appeared to improve when command amplitude was reduced, even though the residual ratio was unchanged. | Report RMS and peak residuals both absolutely and normalized by command amplitude. | Revision 3 settled-response rebuild |
| IMP-023 | Resolved | Browser model | Editable inputs saved, but the live transfer equations still added the local detent tangent as a global spring. | Browser and Python now use the same global commutation model and expose the local detent low-pole band separately. | Revision 3 settled-response rebuild |
| IMP-024 | Resolved | Kinematic diagram | The topology had only a column grid, several dangling ground-referenced ports, floating compliance percentages, scattered case tags, and approximate rather than literal reduction alignment. | Rebuilt it on six row bands and named columns; added one hatched datum per panel, an orthogonal transformer summing node, a stacked compliance bar, a four-by-four source/port matrix, opacity-based ghosting, and vertical full-to-reduced correspondence lines. | Revision 3 diagram-grid rebuild |

## New issue template

Copy one row into **Active issues** and assign the next ID.

| ID | Status | Area | Observed behavior | Expected behavior | Next check |
|---|---|---|---|---|---|
| IMP-XXX | Open | Component or file | What went wrong, including the shortest reproducible condition. | What should have happened. | The next concrete diagnostic or fix. |
