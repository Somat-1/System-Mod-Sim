# Oscillating Stepping Test

`FB_OscillatingStepping.st` contains only the small equilibrium-step oscillation from the ESP32 v5/v6 campaign. All marker moves, 10–25 mm travel legs, throughput tests, and long back-and-forth ramps are omitted.

## Sequence

For requested MRES values 1, 2, 4, 8, 16, 32, and 64:

1. Settle at the captured origin for 5 s before the first level, and add 5 s of settling between levels.
2. Command one requested-MRES microstep forward.
3. Dwell/settle for 5 s.
4. Return to the exact captured origin.
5. Dwell/settle for 5 s.
6. Repeat for 15 cycles.

The final command is the same captured origin, so position commands do not accumulate.

## Important interpretation

The source v5/v6 test does **not** sweep oscillation frequency. It sweeps microstep resolution while using a fixed, very slow forward/return equilibrium test. The Beckhoff terminal remains on its fine 1/64 grid and emulates the requested levels with these target increments:

| Requested level | Counts on 1/64 grid | Mechanical increment |
|---:|---:|---:|
| MRES 1 | 64 | 1 full step / 10 µm |
| MRES 2 | 32 | 1/2 step / 5 µm |
| MRES 4 | 16 | 1/4 step / 2.5 µm |
| MRES 8 | 8 | 1/8 step / 1.25 µm |
| MRES 16 | 4 | 1/16 step / 0.625 µm |
| MRES 32 | 2 | 1/32 step / 0.3125 µm |
| MRES 64 | 1 | 1/64 step / 0.15625 µm |

This makes the mechanical equilibrium positions comparable, but it is not electrically identical to changing a TMC2209's MRES register. The Beckhoff current controller and its 1/64 internal grid remain active throughout.

## Default duration

With 15 cycles, 5 s at each endpoint, and 5 s additional settling between resolution levels, the run takes 1,090 seconds, or approximately 18 minutes 10 seconds.

Run this test before the chirp. Watch the external position sensor and terminal current/temperature, and verify that each return endpoint coincides with the original position.

## Rendered preview

![Complete stepping sequence](Rendered/stepping_sequence_overview.png)

The source is `render_stepping_sequence.py`; the generated file and regeneration instructions are in `Rendered/README.md`.