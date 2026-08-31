# Temporary reduced-model optimization

This workspace calibrates the nonlinear two-master Guyan model against the v1
IDS stepping records for StepSize 1, 2 and 16. StepSize 8 remains excluded
because its command/encoder file was previously identified as faulty.

The objective residual is measured IDS stage position minus simulated reduced
model stage position, in micrometres. The command trajectory uses the edge times
and directions detected from each matching EL5101 file. The configured number of
early cycles is fitted jointly across all three step sizes.

All mechanical structural stiffnesses and dampings are locked and asserted
unchanged during every evaluation. Free positive variables are represented by
bounded log multipliers: detent torque and each LuGre port's sigma0, sigma1,
sigma2, Coulomb level, static-minus-Coulomb gap and Stribeck velocity. This
enforces positive values and Fs > Fc (or Ts > Tc).

Optimization alternates through detent, guideway, nut and support-bearing groups
for every starting point. A weak log-parameter regularizer limits otherwise
severe non-identifiability from fitting three friction ports to one output.

Run:

```text
/usr/bin/python3 optimize_guyan_lugre.py
python3 plot_optimization_results.py
```

The first command uses the system Python installation that contains SciPy on
this machine. The second generates the result figures after the numerical run.
