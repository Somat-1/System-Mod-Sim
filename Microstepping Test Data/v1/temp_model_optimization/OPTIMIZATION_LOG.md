# Optimization setup and refinement log

## Fixed model structure

The optimizer uses the two-master Guyan model with retained coordinates
`theta_m` and `x_n`. Every mechanical structural stiffness and damping value is
copied from Rev 4.2 and asserted unchanged at every objective evaluation:

`k_c`, `k_s1`, `k_s2`, `k_nut`, `k_brg`, `c_c`, `c_s1`, `c_s2`, `c_nut`,
`c_brg`, and `c_EM`.

The nut structural path remains in parallel with nut LuGre friction. The time
model uses the exact `T_d*sin(4*N_r*theta_m)` detent torque.

## Free variables

Nineteen bounded log multipliers are varied: `T_d` plus `sigma0`, `sigma1`,
`sigma2`, Coulomb level, static-minus-Coulomb gap and Stribeck velocity for each
of the way, nut and support-bearing ports. Gap parameterization guarantees
`Fs > Fc` and `Ts > Tc`.

## Physical data and objective

StepSize 1, 2 and 16 are fitted jointly using their detected physical command
edges and IDS displacement. StepSize 8 remains excluded because its encoder
record was previously rejected. The residual before weighting is always

`IDS measured position [um] - simulated x_n [um]`.

The refined run divides each trajectory's residual by its own robust 5–95%
range so StepSize1 does not overwhelm the smaller trajectories solely through
amplitude. The initial refined search used soft-L1; after the objective
inconsistency described below, final selection and polishing use linear least
squares. A weak log-parameter regularizer limits drift along non-identifiable
directions.

## Optimization schedule

Each start alternates through four groups: detent, way friction, nut friction
and support-bearing friction. The refined run starts at the preserved pilot
optimum and uses two independently perturbed restarts around it. Exact state
Jacobians are supplied to the stiff Radau integrator. `optimization_run_log.jsonl`
records the configuration and every completed group with timestamp, cost,
function-evaluation count and parameter vector.

## Run history

- Pilot: 3 starts, 2 cycles of data, 1 alternating pass, 5 evaluations/group.
  Preserved as `pilot_optimization_*`. It reduced its best-start raw cost by
  only about 0.0027% and underpredicted StepSize16.
- Refined: 3 starts, 3 cycles of data, 2 alternating passes, 8 evaluations/group,
  balanced trajectory weighting and soft-L1 loss. Results overwrite only the
  non-prefixed `optimization_*` files.

## Refined-run objective correction

The completed soft-L1 run revealed that its robust internal objective could
accept a terminal point whose recorded weighted squared physical residual was
worse than its start. The terminal files are preserved as
`refined_terminal_optimization_*`. The selected result was recovered from every
logged start and group endpoint using one consistent weighted squared residual;
start 2 at cost 7.8860615 was the lowest logged candidate. Future executions use
linear least squares plus an explicit monotonic physical-cost acceptance guard,
so a group update cannot worsen that reported objective.

The recovered three-start selection is preserved as `selected_restart_*`, and
the complete historical JSONL stream as `refined_softl1_run_log.jsonl`. A final
single-start polish then begins from that selected restart, uses linear least
squares, and applies the monotonic acceptance guard. This is a polishing stage
after the documented multi-start search, not a replacement for it.

## Final monotonic polish

The polish completed one alternating pass (detent, way, nut and support-bearing
groups), with five function evaluations permitted per group. Every update was
accepted and the consistent objective decreased monotonically:

- selected restart: 7.8860615
- after detent: 7.8859052
- after way friction: 7.8859041
- after nut friction: 7.8857728
- after support-bearing friction: 7.8856681

Runtime was 4083 s. The final physical RMSE values are 1.03522 um, 0.56340 um
and 0.12880 um for StepSize 1, 2 and 16 respectively. The small improvement
confirms a shallow optimum and does not establish unique identifiability of the
nineteen nonlinear parameters.
