# Connections — one ESP32-S3 driving the TMC2209 *and* reading the encoder

Everything in `Stage & motor/` assumes a single ESP32-S3 that both **commands**
the stage (STEP/DIR to a TMC2209) and **measures** it (quadrature encoder on the
screw). That pairing is what makes identification possible at all: commanded and
measured position share one clock, so an error is a real error and not two
instruments disagreeing about time.

---

## 1. Read this first: the two existing sketches collide

The encoder firmware and the TMC firmware in this repository were written for
the same board at different times, and **they claim the same two pins**:

| GPIO | `esp32s3_encoder_stream.ino` | `esp32s3_tmc2209_jog.ino` | |
|---|---|---|---|
| 4 | A | — | ok |
| **5** | **B** | **STEP** | **collision** |
| **6** | **Z (index)** | **DIR** | **collision** |
| 17 | — | UART TX → PDN | ok |
| 18 | — | UART RX ← PDN | ok |
| 48 | — | onboard RGB LED | ok |

Flashing either sketch unchanged onto a board wired for both gives nonsense: the
PCNT unit would count DIR transitions as encoder edges, and the index interrupt
would fire on every STEP pulse.

**Resolution used here: the encoder moves, the driver stays.** The TMC pins
(5, 6, 17, 18) appear in every campaign sketch from v2 to v6 and in the prebuilt
binary under `stage identification/jog/firmware/prebuilt/`. Moving them would
invalidate all of it. The encoder appears in one recent sketch, so it is the
cheaper side to change — two jumper wires.

---

## 2. Pin map used by everything in this folder

| Signal | GPIO | Direction | Note |
|---|---|---|---|
| **TMC2209 STEP** | **5** | out | rising edge = one microstep |
| **TMC2209 DIR** | **6** | out | set ≥20 µs before the edge |
| **TMC2209 PDN_UART** | **17** TX / **18** RX | both | single wire via 1 kΩ, driver address 0 |
| TMC2209 EN / ENN | — | — | **grounded externally**; no enable line in firmware |
| **Encoder A** | **4** | in | PCNT edge input — *unchanged* |
| **Encoder B** | **7** | in | **moved from GPIO5** |
| **Encoder Z (index)** | **15** | in | **moved from GPIO6** |
| Onboard RGB LED | 48 | out | status |
| GND | GND | — | common ground is required, not optional |

```
        TMC2209                     ESP32-S3                    Encoder
      ┌──────────┐              ┌───────────────┐            ┌──────────┐
      │     STEP ├──────────────┤ GPIO5   GPIO4 ├────────────┤ A        │
      │      DIR ├──────────────┤ GPIO6   GPIO7 ├────────────┤ B        │
      │ PDN_UART ├──[1k]──┬─────┤ GPIO17  GPIO15├────────────┤ Z        │
      │          │        └─────┤ GPIO18        │            │          │
      │  EN/ENN  ├── GND        │           GND ├────────────┤ GND      │
      │      GND ├── GND        │           3V3 ├────────────┤ V+       │
      └──────────┘              │           USB ├── host     └──────────┘
                                └───────────────┘
```

### Why these two replacements

GPIO7 and GPIO15 are ordinary pins on the ESP32-S3. They are not strapping pins
(0, 3, 45, 46), not the USB data pins (19, 20), and not used for flash or PSRAM.
GPIO1 and GPIO2 were left free on purpose — earlier sketches use them as a
Dewesoft trigger out/echo pair, and that may come back.

PCNT reaches both A and B through the GPIO matrix, so neither needs to be a
particular pin. Only the pair matters: A is the edge input with B as level, and
the second channel is wired the other way round, which is what makes the count
x4 and direction-aware.

---

## 3. Levels, power and the things that quietly break

**The ESP32-S3 is 3.3 V and is not 5 V tolerant.** Check the encoder's output
before connecting it. A 5 V push-pull encoder needs a divider or a level shifter
on A, B and Z; a 5 V open-collector encoder can be pulled up to 3.3 V instead
and wired directly.

**Common ground is required.** The motor supply, the driver and the ESP32 must
share a ground reference or the STEP/DIR levels are undefined.

**`R_SENSE` is 0.11 Ω** — the resistor actually fitted on this board. v4 used
0.03 (a different board's value) and under-delivered current by about 3×, which
is the kind of error that looks like a mechanical problem for a week.

**Motor supply before logic.** Powering the ESP32 while the driver is unpowered
leaves STEP/DIR driving an unpowered input.

---

## 4. Driver configuration, and why it is checked rather than assumed

Both sketches here force **SpreadCycle** and disable **MicroPlyer interpolation**,
then verify both by reading `GCONF`, `CHOPCONF` and `DRV_STATUS` back. A run that
starts is therefore evidence of the configuration rather than an assumption
about it.

This is not defensive habit, it is a specific scar. On the TMC2209 both relevant
bits reset to the wrong value for this work:

| Register bit | Resets to | What "leave it alone" gives you |
|---|---|---|
| `GCONF` bit 2 `en_spreadCycle` | 0 | StealthChop |
| `CHOPCONF` bit 28 `intpol` | 1 | interpolation **on** |

Interpolation on means every commanded microstep is smeared into the interval to
the next one — which is exactly the thing a microstep staircase sets out to
measure. See
[`../TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md`](../TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md)
for the full account, including the 41-minute run lost to an unvalidated UART
read.

---

## 5. Scale

| Quantity | Value | Source |
|---|---|---|
| Screw pitch | 2 mm/rev | rig |
| Full steps per revolution | 200 | motor |
| **Full steps per mm** | **100** | ⇒ 1 full step = 10 µm |
| Electrical period | 4 full steps = 40 µm | 2-phase stepper |
| Encoder counts/rev (nominal) | 20 000 (5000 PPR × 4) | assumed |
| Encoder counts/rev (measured) | 20 052 over one revolution | `Stepper Isolated` recording |
| Stage travel | 44 mm | `tmp_v4_manual_range_sweep` |

The 0.26 % gap between nominal and measured counts/rev is **not yet resolved** —
it survives filtering, so it is systematic rather than noise, and it is one of
the things step 1 of the identification sequence exists to settle. Until then,
treat 20 000 counts/rev as provisional.
