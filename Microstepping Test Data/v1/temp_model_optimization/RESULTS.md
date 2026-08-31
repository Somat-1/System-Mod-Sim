# Optimization results

## Selected refined result

The corrected three-start search selected start 2 at a weighted squared
physical residual of 7.8860615. A final single-start linear least-squares polish
with explicit monotonic acceptance reduced this to **7.8856681**. All four
parameter-group updates were accepted, and no locked structural parameter was
changed.

Final physical trajectory RMSE:

- StepSize 1: 1.03522 um
- StepSize 2: 0.56340 um
- StepSize 16: 0.12880 um

The small 0.005% objective reduction during polishing shows that the selected
restart was already on a shallow local minimum. These friction and detent values
should therefore be treated as an effective fit to these three stepping records,
not as a uniquely identified set of physical LuGre parameters.

The final trajectory comparison is in `plots/optimized_stepping_montage.png`,
the corresponding linearized response in `plots/optimized_bode.png`, and the
full corrected restart history in `plots/refined_multistart_selection.png`.

## Earlier pilot

The configured pilot completed three starts, one alternating pass and two
measured cycles from each of StepSize 1, 2 and 16.

The best result came from start 1:

- initial least-squares cost: 61.045615
- optimized least-squares cost: 61.043960
- relative reduction from that start: approximately 0.0027%
- numerical runtime for the best start: 714.8 s

The optimization machinery is functioning, checkpoints correctly and preserves
all locked structural stiffness and damping values. However, this pilot does not
produce a strongly identified physical parameter set. The StepSize 1 and 2
plateaus are reproduced reasonably, but the optimized model substantially
underpredicts StepSize 16. The cost surface is shallow and the three starts end
at similar costs.

The pilot is preserved as `pilot_optimization_*`; the non-prefixed result files
contain the selected and polished refined result above.
