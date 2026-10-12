# ESP32S3_TMC2209_Jog — drive the stage to two positions and shuttle between them

Interactive positioning for the ESP32-S3 / TMC2209 rig: type two numbers, and
the stage moves back and forth between them until you stop it.

```
ESP32S3_TMC2209_Jog/
├── README.md                         ← this document
├── jog_stepper.py                    ← the host script: find the port, prompt, drive
└── firmware/esp32s3_tmc2209_jog/
    └── esp32s3_tmc2209_jog.ino       ← must be flashed to the board first
```

---

## 1. Why there are two pieces

The request was for a standalone Python script, and `jog_stepper.py` is one —
single file, single dependency. But a host script cannot pulse STEP and DIR: on
this rig those are **ESP32 pins**, four metres of USB away from Python. Something
on the board has to turn a position request into a pulse train.

None of the existing sketches can do it. They are campaign firmware: they run a
fixed, pre-programmed trajectory from boot and accept only `RUN`, `STATUS` and
`ABORT`. There is no "go to position *n*" among them. So this folder adds the
missing half.

Flash the sketch once; after that the Python script is the only thing you touch.

If you want to see the interface before wiring anything up:

```bash
python jog_stepper.py --simulate
```

which runs it against a model of the firmware. Nothing moves, and no hardware
or pyserial is needed.

---

## 2. Wiring

Unchanged from the rig's established wiring, taken from
[`esp32_v6_mres_trajectory_campaign`](../Microstepping%20Test%20Data/v6/scripts/esp32_v6_mres_trajectory_campaign)
so the same board runs both without rewiring.

| TMC2209 | ESP32-S3 | Note |
|---|---|---|
| STEP | **GPIO5** | |
| DIR | **GPIO6** | |
| PDN_UART | **GPIO17 (ESP TX)** via 1 kΩ, **GPIO18 (ESP RX)** | single-wire UART, address 0 |
| EN / ENN | **grounded externally** | always enabled; there is no enable line to toggle |
| — | GPIO48 | onboard RGB LED |

`R_SENSE` is **0.11 Ω**, the resistor actually fitted on this board. This matters:
v4 used 0.03 (a different board's value) and under-delivered current by about 3×.

---

## 3. Quickstart

```bash
pip install pyserial

python jog_stepper.py                      # find the board, then prompt
python jog_stepper.py --port COM11         # skip the search
python jog_stepper.py --simulate           # no hardware
python jog_stepper.py --a -200 --b 200 --go --max-legs 10
```

A session looks like:

```
jog> -200 200        set both positions
jog> speed 70        full steps per second
jog> go              shuttle until a key is pressed
```

### Finding the port

With no `--port`, the script lists every serial port and probes them in order of
likelihood — your `--port` if given, then `COM11`, then anything with a known
USB-serial vendor ID (Espressif `303A`, CP210x, CH34x, FTDI), then the rest. Each
candidate gets a `PING`, and only a board that answers is used, so it cannot
latch onto the wrong device.

Worth knowing: **probing opens each candidate port, and opening a port resets an
ESP32.** On a bench with more than one board that will restart the others
harmlessly, but if something is mid-run, pass `--port` so nothing else is
touched. `--list` shows the ports without probing any of them.

---

## 4. Commands

| Command | Effect |
|---|---|
| `-200 200` | set both positions at once |
| `a <n>` / `b <n>` | set one of them |
| `go` | shuttle between A and B until a key is pressed |
| `go <n>` | ... or for exactly *n* legs |
| `moveto <n>` | a single move to an absolute position |
| `speed <steps/s>` | full-step rate; 70 = 0.7 mm/s |
| `mres <1..256>` | microstep resolution, powers of two, when idle |
| `current <mA>` | coil current, RMS |
| `limit <steps>` | soft limit, ± about zero; `0` disables |
| `zero` | call the present position 0 |
| `status` | ask the board where it is |
| `help`, `quit` | |

While jogging you get a live position readout, a travel bar and a leg count.
**Any key stops it**, as does Ctrl-C.

---

## 5. Units

Positions are whole **full steps**. On this rig 100 full steps = 1 mm, so

| 1 full step | 1 mm | stage travel |
|---|---|---|
| 10 µm | 100 full steps | 44 mm measured |

Every prompt shows the millimetre equivalent, so `-200 200` is announced as
−2.00 mm to +2.00 mm. Speed is a full-step rate in steps/s, matching the
campaign sketches' convention — their `MODERATE` is 70 steps/s, which is the
default here.

Microstepping affects smoothness, not scale: at `mres 16` one full step is 16
pulses, so a given full-step rate means 16× the pulse rate but the same stage
speed.

---

## 6. Safety

This moves a real stage that has **no homing and no end stops**, so:

- **Position 0 means "wherever the stage was when it was last zeroed"**, nothing
  more. There is no absolute reference, and power-cycling the board loses it.
- **A soft limit of ±500 full steps (±5 mm) about zero is on by default.** The
  firmware refuses targets outside it, and the host refuses to send them. The
  measured travel is 44 mm, so there is room to widen it with `limit <steps>`
  once you know where the stage is sitting within that — but the default is
  deliberately tight, because a mistyped target is otherwise a full-speed run
  into an end stop.
- **Nothing moves at boot or at connect.** The campaign sketches start their
  trajectory automatically, which is only safe when clearance is already known;
  this one always waits to be told.
- `STOP` acts within one pulse period. Step generation is non-blocking — one
  `loop()` pass emits at most one pulse — unlike the campaign sketches, which
  busy-wait inside a move and can only check for an abort every 32 pulses.

---

## 7. Driver configuration

SpreadCycle is forced, MicroPlyer interpolation disabled, and both are
**verified by readback** before any move is accepted — `GCONF`, `CHOPCONF` and
`DRV_STATUS` must all agree, or the firmware reports
`ERR,reason=driver_mode_not_confirmed` and refuses to move.

This follows
[`TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md`](../TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md):
on this chip `en_spreadCycle` resets to 0 and `intpol` resets to 1, so leaving
the registers alone does not mean "neutral", it means StealthChop with every
microstep smeared into the next. Register reads are taken twice and required to
agree, because a failed TMC2209 UART read returns garbage rather than an error —
that is what cost v4 a 41-minute run.

---

## 8. Protocol

Newline-terminated ASCII in, `#`-prefixed lines out. Positions in whole full
steps, rates as a full-step rate in milliHz.

| Host → board | Board → host |
|---|---|
| `PING` | `# PONG,fw=…,ready=1` |
| `STATUS` | `# STATUS,pos=…,a=…,b=…,mode=…,…` |
| `POS?` | `# POS,pos=…` |
| `ZERO` | `# ZERO_OK,pos=0` |
| `LIMIT <steps>` | `# LIMIT_OK,limit=…` |
| `RATE <milliHz>` | `# RATE_OK,rate_millihz=…,pulse_hz=…` |
| `MRES <n>` | `# MRES_OK,mres=…` |
| `CURRENT <mA>` | `# CURRENT_OK,current_ma=…` |
| `MOVETO <steps>` | `# MOVETO,target=…` then `# ARRIVED,pos=…` |
| `JOG <a> <b>` | `# JOG_SET,a=…,b=…` |
| `START` | `# START,…`, then `# REACHED,…` each leg |
| `STOP` | `# STOP,pos=…,reason=…` |

Unsolicited: `# POS,pos=…,mode=…` at 10 Hz while moving, and `# ERR,reason=…`
for anything refused. Position is held internally in 1/256 full-step units so
every supported MRES divides it exactly and changing MRES cannot accumulate a
rounding drift.

---

## 9. What has and has not been tested

Worth being explicit, because this was written on a machine with no board attached.

**Tested.** The host script end to end against the built-in simulator: port
discovery and ranking, the PING handshake, the prompt, setting positions,
jogging with leg counting and reversal, `moveto`, `zero`, `speed`, `mres`,
`current`, `limit`, the soft-limit refusals, and the malformed-input paths. The
simulator enforces the same validation as the firmware — same soft limit, same
power-of-two MRES, same current range, same pulse-rate ceiling — so that it
cannot accept something the real board would reject.

**Not tested.** The firmware has never been compiled, let alone run: there is no
Arduino toolchain on this machine and no board on the bench. It needs
**TMCStepper** installed and the **ESP32** board package. Expect to fix something
on the first build.

**Worth checking on the first run**, before trusting it with a long move:

1. `# READY,ready=1` appears at connect. `ready=0` means the TMC2209 UART did
   not answer — check GPIO17/18 and the motor supply.
2. `moveto 10` moves 0.1 mm in the direction you expect. If it goes the wrong
   way, swap the motor pair or invert `DIR_PIN` in the sketch.
3. The distance is right: `moveto 100` should be exactly 1 mm. If it is out by a
   constant factor, `FULL_STEPS_PER_MM` does not match the fitted screw.
4. Only then widen `limit` and run a long jog.
