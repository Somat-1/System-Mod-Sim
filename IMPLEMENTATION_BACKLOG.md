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
| IMP-029 | Open | GMS integrator | A2 final/settled-window RMS diverges under step halving (observed order about -1.4) because branch switching is evaluated only at RK trial states. | Hybrid branch events should be localized so reported A2 metrics have a bounded discretization error. | Add event localization or an equivalent hybrid integrator, then repeat the 50/25/12.5/6.25 µs A2 sequence. |

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
| IMP-019 | Resolved | Source interpretation | The installed screw class and travel geometry were not recorded correctly. | Record IT3 from the BOM, 192 mm complete screw length, approximately 170 mm usable distance, and 150 mm stage travel; remove the unnecessary Section 12.1 comparison. | Revision 3 1/16-microstep rebuild |
| IMP-020 | Resolved | Parameter provenance | The coupling estimate remained 1.20×10⁻⁶ kg·m² after the annulus estimate had been revised. | Set $J_c=1.18×10^{-6}$ kg·m² and regenerated reflected mass and modal results. | Revision 3 critical-error rebuild |
| IMP-022 | Resolved | Reduction verification | Nanometre residuals appeared to improve when command amplitude was reduced, even though the residual ratio was unchanged. | Report RMS and peak residuals both absolutely and normalized by command amplitude. | Revision 3 settled-response rebuild |
| IMP-023 | Resolved | Browser model | Editable inputs saved, but the live transfer equations still added the local detent tangent as a global spring. | Browser and Python now use the same global commutation model and expose the local detent low-pole band separately. | Revision 3 settled-response rebuild |
| IMP-024 | Resolved | Kinematic diagram | The topology had only a column grid, several dangling ground-referenced ports, floating compliance percentages, scattered case tags, and approximate rather than literal reduction alignment. | Rebuilt it on six row bands and named columns; added one hatched datum per panel, an orthogonal transformer summing node, a stacked compliance bar, a four-by-four source/port matrix, opacity-based ghosting, and vertical full-to-reduced correspondence lines. | Revision 3 diagram-grid rebuild |
| IMP-025 | Resolved | Source/render consistency | The Markdown fallback for the derived STEP/DIR quantum disagreed with the builder-rendered HTML. | Corrected both Rev 3 Markdown sources to the production 312.5 nm value and retained `command_step` as a derived output sourced from the builder. | Revision 3 1/16-microstep rebuild |
| IMP-026 | Resolved | Build performance | Independent nonlinear case and step-halving RK4 runs executed sequentially and left available CPU cores idle. | Added a bounded process pool and a `--jobs` option; simulations remain deterministic and artifact writes stay in the parent process. | Revision 3 1/16-microstep rebuild |
| IMP-021 | Resolved | Hardware configuration | The executed external microstep setting was not recorded, so the model used an idealized finer grid. | Record the production Stepper-Board limit as 1/16 and rebuild every command-dependent result on a 312.5 nm grid. | Revision 3 1/16-microstep rebuild |
| IMP-027 | Resolved | Render-only consistency | Rendering HTML in a fresh process replaced the live GMS census sentence with a “not rebuilt” notice even though the generated Markdown table contained current results. | Recover the live sentence from the generated Markdown census when no in-memory simulation result is present. | Revision 3 1/16-microstep rebuild |
| IMP-006 | Resolved | Generated artifacts | A rebuild changed SVG metadata and internal IDs even when a figure's data and layout were unchanged. | Suppressed SVG date metadata and fixed Matplotlib's SVG hash salt; verified with two consecutive final builds. | Revision 3 Section 9 repair rebuild |
| IMP-028 | Resolved | Generator ownership | Section 9 and the rendered HTML regressed because corrections were not represented in every generator-owned block and case map. | Restored G/G2 and all requested prose/metric changes in the Python generator and Markdown sources, added a report build ID, and verified consecutive generated outputs. | Revision 3 Section 9 repair rebuild |
| IMP-030 | Resolved | HTML navigation | The fixed left outline consumed reading width and could not be hidden. | Added a persistent Hide/Show outline control that collapses the left column without changing document content. | Revision 3 Section 9 final rebuild |
| IMP-031 | Resolved | Figure readability | Dense SVG diagrams and response plots could not be inspected beyond their inline page width. | Made report images keyboard-accessible and added an expanded lightbox with zoom, mouse-wheel zoom, drag-to-pan, reset, and close controls. | Revision 3 Section 9 final rebuild |
| IMP-032 | Resolved | Parameter persistence | Browser-edited values were not a reliable input to the next Python plot-generation run. | Added a single Save parameters action that writes `model_parameters.json`; the builder loads that file before constructing executable defaults and creates a complete handoff file on first run. | Revision 3 Section 9 final rebuild |
| IMP-033 | Resolved | Retired generated content | Retiring the IDS report subsection left generator code and generated SVG artifacts able to preserve or recreate stale content. | Removed the IDS analysis/render path and its two generated SVGs, while leaving the unprocessed measurement records untouched. | Revision 3 Section 9 final rebuild |

## New issue template

Copy one row into **Active issues** and assign the next ID.

| ID | Status | Area | Observed behavior | Expected behavior | Next check |
|---|---|---|---|---|---|
| IMP-XXX | Open | Component or file | What went wrong, including the shortest reproducible condition. | What should have happened. | The next concrete diagnostic or fix. |
