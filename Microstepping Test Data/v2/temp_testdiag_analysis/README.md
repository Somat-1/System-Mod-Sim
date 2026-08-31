# QUICK12 testDiag assessment

Source: `../data/testDiag.csv`, 130,123 EL5101 samples at 1 kHz.

The reconstructed 107.5 s QUICK12 command was aligned at approximately 2.45 s
after recording start. The exact absolute count-to-length calibration is not
documented in the export, so amplitudes are assessed in encoder counts.

## Result: rejected

This diagnostic is **not suitable for identification**, and the long campaign
must not be started from this result. The earlier provisional acceptance was
incorrect.

- Baseline robust noise sigma: 5.93 counts.
- Commanded-interval 1--99% span: 195 counts.
- Global command/alignment score: only 0.218.
- Section correlations and inferred gains vary strongly between experiments.
- Native half-step A2: correlation magnitude 0.589; span 104 counts.
- Native half-step B: correlation magnitude 0.404; span 56 counts.
- Whole-step endpoints A2: correlation magnitude 0.635; span 47 counts.
- Whole-step endpoints B: correlation magnitude 0.514; span 56 counts.

A2 and B contain visually structured changes, but that alone does not establish
correct stage tracking. A1 does not reproduce the commanded staircase, E is
dominated by transient spikes, and the inferred count-per-step gains are neither
consistent across blocks nor proportional between half-step and whole-step
endpoints. These failures invalidate quantitative use of the recording.

Large counter discontinuities occur at both ends of the aligned interval and are
far larger than the purported motion signal. Their source is unresolved. Merely
excluding them from the metrics does not make the intervening data trustworthy.

## Required before another campaign

Do not execute the existing four-resolution `RUN`. Native TMC2209
full-step MRES code 8 did not persist on the attached module, while half-step
MRES code 7 was repeatedly verified. The long ESP campaign must therefore:

1. exclude 1/16 and 1/4 configurations;
2. use verified native half-step for the half-step data;
3. label whole-step endpoints explicitly and generate them as paired contiguous
   half-step pulses unless native full-step register control is repaired;
4. retain hard MRES readback checks and abort on mismatch;
5. retain the large configuration/phase markers and at least a few seconds of
   IDS lead-in and tail recording;
6. run one simple IDS validation trajectory first: large one-direction steps
   with multi-second settled plateaus, followed by an explicit return;
7. require consistent measured step sign, monotonic plateaus, return-to-origin,
   and approximately proportional displacement before testing A1/A2/B/E.

Plots:

- `plots/00_testdiag_alignment.png`
- `plots/01_testdiag_sections.png`
