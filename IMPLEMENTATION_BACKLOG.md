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
| IMP-003 | Resolved | Response comparison | Seven similar case figures obscured the useful Bode differences. | Grouped the cases into one overlay with a resonance zoom and labeled pairwise differences. | Revision 3 rebuild |
| IMP-004 | Resolved | Technical writing | Long passages and em dashes made the analytical report harder to scan. | Shortened key explanations, replaced em dashes, and moved reference material to appendices. | Revision 3 rebuild |
| IMP-007 | Resolved | Parameter dependencies | Derived values such as $r$, $m_d$, and $K_m$ were editable independent inputs and could become inconsistent. | Made component values the inputs and recalculated dependent outputs and inline reduced-model results in the browser. | Revision 3 component-value rebuild |
| IMP-008 | Resolved | Drive model | Detent torque was disabled even though a published 0.005 N·m value was available. | Enabled the nonlinear detent force and its equilibrium tangent in the linear models. | Revision 3 component-value rebuild |
| IMP-009 | Resolved | Generated SVGs | Matplotlib emitted trailing spaces that caused repository whitespace checks to fail. | Normalized SVG lines in the shared writer after every figure save. | Revision 3 component-value rebuild |
| IMP-010 | Resolved | Kinematic diagram | The axial load path incorrectly placed $u_f$ between $u_e$ and $u_n$, and $\theta_{s3}$ appeared inline. | Redrew both coordinates as dead-end overhang stubs and aligned the remaining coordinates by physical station. | Revision 3 topology-diagram rebuild |
| IMP-011 | Resolved | Nut friction port | The full equations applied $T_{f,n}$ only to $\theta_{s2}$, while the reduced equations used an internal equal-opposite port. | Defined $v_n=r\dot\theta_{s2}+\dot u_e-\dot u_n$, added the paired axial forces, and documented $T_{f,n}=rF_{f,n}$. | Revision 3 topology-diagram rebuild |
| IMP-012 | Resolved | Reduction diagram | Mass destination and compliance retention shared one fill color, which incorrectly mapped $u_n$ and implied that axial screw masses survived. | Made box fill encode mass aggregation and spring stroke encode compliance retention. | Revision 3 reduction-map rebuild |
| IMP-013 | Resolved | Diagram registration | Station headers, the nut transformer, and axial coordinates did not share a column grid; several numeric labels collided with elements. | Added station guides, aligned the nut column, and moved reduction evidence into a dedicated annotation column. | Revision 3 reduction-map rebuild |
| IMP-014 | Resolved | Friction topology | A 5 N gross nut-drag level was assigned to the differential elastic rate $\dot x_d-\dot x_s$, which vanishes during common motion. | Split gross nut rolling drag onto $v_d$ and retained a lower-force differential microslip port with a 0.20 µm first-yield distance. | Revision 3 critical-error rebuild |
| IMP-015 | Resolved | Executed cases | $F_{f,d}$ was defined but disabled in every case, so modeled drivetrain loss was silently zero. | Activated $F_{f,d}$ in A/A2, B/B2, and C/C2 and updated the diagrams and case map. | Revision 3 critical-error rebuild |
| IMP-016 | Resolved | Command generation | The main and memory tests used 4x and 32x steps despite the stated TMC2209 hardware resolution. | Changed the main sequence to a 64x STEP/DIR quantum and exposed both 64x and 256x derived values in the browser. | Revision 3 critical-error rebuild |
| IMP-017 | Resolved | Response table | The zero-speed linearization was labeled "DC gain", which could be mistaken for full-range tracking gain. | Renamed it presliding tangent gain, added first-yield validity travel, and stated that nonlinear travel produces bounded offsets. | Revision 3 critical-error rebuild |
| IMP-018 | Resolved | Memory experiment | The original displacement discriminator was close to the interferometer floor and yielded only one GMS element. | Made force the primary metric and increased the 64x-quantized excursion to cross two nominal guideway yield distances. | Revision 3 critical-error rebuild |
| IMP-019 | Resolved | Source interpretation | The report request attributed 6 µm over 315 mm to IT3. The manufacturer table assigns 6 µm to IT1 and 12 µm to IT3. | Corrected the class mapping, elevated lead error, and left the waveform unexecuted until the installed class or measured map is known. | Revision 3 critical-error rebuild |
| IMP-020 | Resolved | Parameter provenance | The coupling estimate remained 1.20×10⁻⁶ kg·m² after the annulus estimate had been revised. | Set $J_c=1.18×10^{-6}$ kg·m² and regenerated reflected mass and modal results. | Revision 3 critical-error rebuild |

## New issue template

Copy one row into **Active issues** and assign the next ID.

| ID | Status | Area | Observed behavior | Expected behavior | Next check |
|---|---|---|---|---|---|
| IMP-XXX | Open | Component or file | What went wrong, including the shortest reproducible condition. | What should have happened. | The next concrete diagnostic or fix. |
