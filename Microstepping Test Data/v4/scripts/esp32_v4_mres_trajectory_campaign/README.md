# ESP32 v4 MRES trajectory firmware

This folder is an Arduino sketch folder. See the v4 root
[README](../../README.md) for the full sequence, hardware pin map, measured
command-rate preflight, duration, build/upload commands, and safety notes.

The sketch starts once automatically after boot. It leaves StallGuard and
CoolStep untouched, but explicitly forces **SpreadCycle** and **disables
MicroPlyer interpolation**, verifying both by readback before it will run.

Neither may be left at its default: `en_spreadCycle` powers up 0 (StealthChop,
lower microstep positional fidelity) and `intpol` powers up 1 (each commanded
microstep smeared across the interval to the next one). See the v4 root
[README](../../README.md) for the measurements behind this.

It configures MRES over the verified TMC2209 UART connection and drives motion
through STEP GPIO 5 and DIR GPIO 6 while EN/ENN is externally grounded.
