# ESP32 v6 MRES trajectory firmware

This folder is an Arduino sketch folder, forked directly from
`../../../v4/scripts/esp32_v4_mres_trajectory_campaign`. See the v4 root
[README](../../../v4/README.md) for the full sequence, hardware pin map,
measured command-rate preflight, marker scheme, and build/upload commands —
all unchanged here. See the [v6 root README](../../README.md) for what's
different.

In short: `R_SENSE_OHM` is corrected to 0.11 (the resistor actually fitted on
this board, not v4's 0.03), and the trajectory legs are 10 mm each way, not
v4's 25 mm — the stage's actual travel was since measured at 44 mm total, so
25 mm one-way never fit.

The sketch starts once automatically after boot. It leaves StallGuard and
CoolStep untouched, but explicitly forces **SpreadCycle** and **disables
MicroPlyer interpolation**, verifying both by readback before it will run.

Neither may be left at its default: `en_spreadCycle` powers up 0 (StealthChop,
lower microstep positional fidelity) and `intpol` powers up 1 (each commanded
microstep smeared across the interval to the next one). See the v4 root
[README](../../../v4/README.md) for the measurements behind this.

It configures MRES over the verified TMC2209 UART connection and drives motion
through STEP GPIO 5 and DIR GPIO 6 while EN/ENN is externally grounded.

**This sketch has no range-finding or centering step and starts moving
automatically at boot.** A 10 mm one-way leg only fits safely if the stage
already has at least 10 mm of clearance on both sides of wherever it's
sitting when flashed — confirm that first (e.g. with
`../../../v4/scripts/tmp_v4_manual_range_sweep`) rather than assuming it.
