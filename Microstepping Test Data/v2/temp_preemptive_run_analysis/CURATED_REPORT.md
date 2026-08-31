# Curated preemptive-run assessment

## Best synchronized evidence (not an accepted quantitative fit)

Block A2 at N = 1, 2, 4, 8, 16, and 32 pulses is the clearest synchronized
family. All conditions show the alternating timing and reproduce across the two
completed current settings. However, the raw plateau excursion remains roughly
similar while commanded travel changes by 32x. Therefore A2 is retained only as
timing/shape evidence, not as a valid quantitative model fit.

## Model overlay convention

The ESP commands drive the Rev 4 frictionless tangent, exact periodic-detent
frictionless, and Rev 4.2 parallel-LuGre models. Model output is stage coordinate
`x_n`. The EL5101 CSV contains counter values but no documented counts-to-length
calibration, so each measured panel is offset and amplitude-normalized to its own
two plateau medians. The overlays test timing and normalized response shape, not
absolute gain. A physical gain fit would be misleading until counter calibration
is supplied. The command-inferred ratios span 6.02 to 257
counts/µm, confirming that one global inferred scale is not defensible here.
`04_A2_gain_consistency_rejection.png` exposes this rather than normalizing it away.
`05_A2_single_scale_model_overlay_rejection_montage.png` uses one scale anchored
at N=16 and shows the resulting amplitude mismatch at every other step size.

## Bottom line

No segment in this preemptive capture supports an absolute measured-versus-model
fit. The A2 montage is the best timing-aligned comparison available; its normalized
overlay must not be interpreted as parameter validation or gain agreement.

## Excluded from fit claims

- A1: dominated by drift/creep and inconsistent monotonic response.
- E: spike-dominated; doublets are too fast relative to this acquisition for a
  trustworthy displacement comparison.
- B descending/minor: retained only as secondary qualitative evidence.
- Block 0: useful as a repeatability fingerprint, not a clean model-gain test.
