# Firmware for the MRES test

`mres_test/` holds the sketch. It is built on the stage-identification sketch
([`../../firmware/stage_identification`](../../firmware/stage_identification/stage_identification.ino))
— same pins, same PCNT encoder — and adds what that one could not do: u256
targets, k-microstep staircase increments, MSCNT in telemetry, and settable
IHOLD / TPOWERDOWN / intpol.

```
mres_test/
├── mres_test.ino     the runner
└── mres_plan.h       GENERATED from ../sequence/mres_plan.csv
```

> **Not compiled, never run on hardware.** No toolchain was available where it
> was written. Compile it, set the two switch pins, and bring it up with the
> motor uncoupled before trusting it with the stage.

## Before the first run

1. **Confirm the end-switch pins.** `SW_NEG_PIN` is GPIO10 and `SW_POS_PIN`
   is GPIO11, as wired on the bench. The switches are 10T85 microswitches with COM on common
   ground and NO to the GPIO, so `SWITCH_ACTIVE_LEVEL` is `LOW` and the pins use
   the internal pull-up. Menu option **8** reads both switches so you can press
   each one by hand and confirm which is which before anything moves.
2. **Check the direction sense.** `SW_NEG_PIN` must be the switch reached by
   moving in the DIR-low direction. If homing seeks the wrong way, swap the two
   pin numbers.
3. **Confirm the MSCNT full-step convention.** `MSCNT_FULL_STEP_OFFSET = 128`
   assumes a full step sits at MSCNT 128, 384, 640, 896. Check it against the
   datasheet for your table; `set_mres` refuses to write anywhere else, so a
   wrong constant shows up as every MRES change failing, not as silent drift.

The switch wiring as fitted is **not fail-safe**: a cut wire reads as "clear"
and the stage would drive into the hard stop. Moving each wire to the NC
contact and setting `SWITCH_ACTIVE_LEVEL` to `HIGH` would make a broken cable
stop the stage instead.

## Connecting a terminal

Use the control terminal in this folder:

```bash
/usr/bin/python3 scripts/mres_terminal.py
```

It opens the port with **DTR and RTS off**, sends each key as it is pressed,
hides the 50 Hz position samples (logging everything to
`data/terminal_<stamp>.log`), turns Ctrl+C into `ABORT`, and closes on Ctrl+].
`arduino-cli monitor` and miniterm both open with DTR on by default; even with
it off, `arduino-cli monitor` showed nothing back in an interactive window.

On the ESP32-S3's built-in USB serial, DTR and RTS drive the chip's reset and
boot strap. Opening with them on, which is the default for `arduino-cli
monitor`, pyserial and miniterm, restarts the chip into the ROM bootloader. It then sits
silently at `waiting for download`, and the terminal looks connected but
nothing answers. If that has happened, `esptool --before no-reset --after
watchdog-reset chip-id` boots the sketch again. `run_mres.py` already opens the
port this way.

## Regenerating the plan

Part 2 runs standalone, so the plan is compiled in. The header is generated
from the same CSV `mres_design.py` writes, which keeps one definition of the
sequence:

```bash
python scripts/mres_design.py        # -> sequence/mres_plan.csv
python scripts/gen_plan_header.py    # -> firmware/mres_test/mres_plan.h
# then reflash
```

Editing `mres_plan.h` by hand defeats the point — it is overwritten.

## The menu

On boot the sketch prints a menu and waits. Nothing moves until a key is
pressed. Use `scripts/mres_terminal.py`, which sends each key the moment it is
pressed; keys act at once, without Enter. Prompts that ask for a number echo
what you type and take Enter.

```
  1  Home              find both end switches, then wait here
  2  Warm-up           shuttle between the limits; any key pauses
  3  MRES sequence     centre, then 133 rows (~30 min); any key stops
  4  Settings          current, speed, MRES, amplitude
  5  Go to centre      move to x_start on a full step
  6  Host plan mode    PLAN ADD / PLAN RUN for run_mres.py
  7  Driver readback   the # DRIVER line
  8  Switch state      read both end switches now
  0  Release motor     EN high, outputs off
```

**Any key stops homing, centring and the sequence.** Stopping the sequence
part-way leaves the block in progress incomplete, so that block has to be rerun.

**Homing waits when it is done.** From there, set the speed, current and
amplitude with 4, then start the warm-up with 2.

**Any key pauses the warm-up** wherever the stage is, and asks what next:

```
  c continue   s settings   r run MRES sequence   m main menu
```

`r` asks you to type `RUN`, then centres the stage and starts the sequence.

**Before the first row, the sequence aligns 0 to an electrical full step.** The
commanded position starts at 0 at boot, but the driver's MSCNT keeps whatever
phase earlier moves left it at, so a commanded full step need not be an
electrical one. When it isn't, every MRES change in the plan is refused
(`mscnt_off_full_step`). After centring, the firmware switches to 1/256, steps
at most half a full step (5 µm) until MSCNT is on a full-step value, and makes
that point the new 0. It reports the shift in `# ALIGNED,...`. The switch
positions and soft limits shift with it. An
MRES change is refused unless the stage is on a full step; the warm-up
endpoints and the centre are both on one.

## Homing

One fast seek finds each switch, then **three slow touches** measure it. The
three must agree within `HOME_TOLERANCE_U256` (20 µm); a wider spread means the
contact is bouncing or soft, and homing retries once and then fails rather than
putting a soft limit in the wrong place. The most conservative touch — the one
that stopped soonest — is the one kept.

```
fast approach  ──────────────▶ |TRIG|
back off       ◀── 0.5 mm
slow approach  ───▶ |TRIG|  ×3, must agree within 20 µm
limit = most conservative trigger, 0.5 mm inside, snapped to a full step
```

Both ends done, `0` is the centre of the measured travel snapped to a full
step, and the soft limits sit a margin inside each switch. `# HOMED` reports
the travel, both trigger positions and both limits. Releasing the motor (`0`)
invalidates the homing — the stage can be pushed by hand once the outputs are
off, so the travel has to be measured again.

Tuning constants are grouped under `// ---- homing` in the sketch:
`HOME_FAST_FS`, `HOME_SLOW_FS`, `HOME_TOUCHES`, `HOME_BACKOFF_U256`,
`HOME_TOLERANCE_U256`, `HOME_MARGIN_U256`, `HOME_MAX_TRAVEL_U256`.

## Protocol

Unchanged from what `run_mres.py` already speaks, and available at the menu
prompt at any time — the host does not need to pick mode 3 first.

```
PING                 -> # PONG,fw=mres_test,ready=1,driver=1,encoder=1,homed=..,moves=..
HOLD <m> <tpd> <hd>  -> # HOLD_OK,ihold_mult=..,tpowerdown=..,iholddelay=..,ihold=..
INTPOL <0|1>         -> # INTPOL_OK,intpol=<0|1>
DRIVER               -> # DRIVER,spreadcycle=..,intpol=..,stealth=..,ihold_mult=..,
                                 tpowerdown=..,iholddelay=..,mres=..,mscnt=..
ZERO                 -> # ZERO_OK,pos=0,enc=..
STATUS               -> # STATUS,pos=..,enc=..,mres=..,homed=..,lim_neg=..,lim_pos=..
PLAN CLEAR / ADD     -> # PLAN_CLEARED / # PLAN_ADD,index=<i>,count=<n>,tag=<tag>
PLAN RUN             -> # RUN_START ... # RUN_COMPLETE,pos=<u256>
HOME                 -> # HOME_OK,lim_neg=..,lim_pos=..
ABORT                -> stops the current motion
```

Plan rows, positions in u256:

```
PLAN ADD <kind> <target_u256> <speed_fs> <dwell_ms> <mres> <n_inc> <inc_usteps> <inc_dwell_ms> <tag>
```

- `kind` is `move`, `stair` or `set_mres`.
- `set_mres` refuses unless the position is a whole full step **and** MSCNT
  reads back on the full-step grid, then emits `# MRES,<mres>,mscnt=<n>`.
- `stair` checks `n_inc * inc_usteps * (256 / mres) == |target - start|`.

Telemetry, at 50 Hz throughout:

```
# S,<t_ms>,<pos_u256>,<encoder>,<mscnt>,<row>
# MOVE,<row>,<kind>,<tag>,target=<u256>,mres=<m>
# ARRIVED,<row>,pos=<u256>,enc=<counts>,mscnt=<n>
```

## Safety

Both switches are read, debounced, before **every** pulse. Outside homing a
trigger stops the motion immediately and reports `# RUN_FAILED` or
`# WARMUP_STOPPED`. Once homed, the soft limits are enforced the same way, so a
move that would cross one stops before it does. `ABORT` is checked in the same
place and takes effect within one pulse.

## Still open

- **Part 1 lives here now, not in `mres_design.py`.** `warmup()` there is still
  an empty stub and `run_mres.py` modes 1 and 3 still refuse. Either point them
  at the firmware's menu or delete them; they currently describe the old plan.
- **Homing is not in the plan CSV**, so a part-2 session driven from the host
  still needs the operator to home from the menu first.
