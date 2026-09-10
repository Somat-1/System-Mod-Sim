# ESP32 v4 MRES trajectory firmware

This folder is an Arduino sketch folder. See the v4 root
[README](../../README.md) for the full sequence, hardware pin map, measured
command-rate preflight, duration, build/upload commands, and safety notes.

The sketch starts once automatically after boot. It deliberately leaves
StealthChop, SpreadCycle, StallGuard, and CoolStep untouched, but explicitly
disables MicroPlyer interpolation, which powers up enabled and would otherwise
smear every commanded microstep across the interval to the next one. It
configures MRES over the verified TMC2209 UART connection and drives motion
through STEP GPIO 5 and DIR GPIO 6 while EN/ENN is externally grounded.
