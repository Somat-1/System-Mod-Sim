# Bench procedure — settling matrix

**Written for: the operator at the bench.** Read §1 before powering anything.
The firmware has never been compiled or run; §2 is bring-up, not a formality.

Total time once bring-up is done: **32 min** of automatic sequence, plus about
3 min of homing.

---

## 1. Before power

**Set the end-switch pins.** `firmware/settling_test/settling_test.ino` ships
with placeholders:

```c
constexpr uint8_t SW_NEG_PIN = 10;
constexpr uint8_t SW_POS_PIN = 11;
constexpr int SWITCH_ACTIVE_LEVEL = LOW;
```

Wiring assumed: 10T85 microswitch, **COM to common ground, NO to the GPIO**,
internal pull-up, so a pressed switch reads LOW. This is **not fail-safe** — a
cut wire reads as clear and the stage will run into the end. Moving both
switches to their NC contacts and setting `SWITCH_ACTIVE_LEVEL` to `HIGH` makes
a broken cable stop the stage instead. Worth doing before an unattended 32 min
run.

**Check the travel.** The plan uses **±20 mm about the centre** of the measured
travel, and homing holds back a 2 mm margin at each end, so the stage needs
**at least 44 mm of usable travel**. The firmware checks this before it moves
and refuses with `# ERR,reason=plan_exceeds_travel` if the plan does not fit.

**Uncouple the motor** for the first bring-up (§2). Homing drives toward the
switches at 2 mm/s and the switch logic has never been exercised.

---

## 2. Bring-up — first time only, motor uncoupled

1. **Compile.** `arduino-cli compile --fqbn esp32:esp32:esp32s3
   firmware/settling_test`. The sketch needs `TMCStepper`. It has been checked
   only statically — brace balance, `printf` arity, 64-bit casts — so expect to
   fix something here.
2. **Flash and open the terminal:**
   `python3 scripts/settling_terminal.py` (it finds the port; pass one to
   override). It logs everything to `data/terminal_<stamp>.log`.
3. **Read the boot lines.** You want:
   - `# DRIVER_ID,...,version=0x21,expected=0x21` — UART talks to the TMC2209.
     Anything else means the UART is wrong, and nothing below is meaningful.
   - `# SWITCHES,neg=0,pos=0,active_level=low` — both switches read clear with
     nothing pressed. If either reads 1, the pin or the contact is wrong.
   - `# READY,driver=1,rows=344,measured=240`.
4. **Press each switch by hand** and press `6` to read them back. Confirm the
   right switch reads 1, and that it returns to 0 when released. Getting these
   swapped sends homing the wrong way at the wrong speed.
5. **Press `5`** and confirm `spreadcycle=1,intpol=0,stealth=0`. Settling under
   StealthChop is a different measurement; if `stealth=1`, stop.
6. **Press `7`.** This is the check that the whole campaign rests on: that the
   microstep size follows the register. It pulses once and five times at each
   size and checks `MSCNT` moved by exactly `256/MRES` per pulse. It needs
   almost no travel (under 4 µm), so it is safe uncoupled.

   ```
     MRES   per pulse   1 pulse   5 pulses   closes   verdict
     1/16          16        16         80      yes   OK
     ...
   # MRES_CHECK,pass
   ```

   A `MISMATCH` means the comparison would be meaningless: check
   `GCONF.mstep_reg_select` is set and `intpol` is off. The matrix runs this
   automatically and refuses to start if it fails, so do not skip it here —
   find out now, not at the start of a 32 min run.

Only once all six pass, couple the motor.

---

## 3. The run

The menu acts on **single keys, no Enter**. `?` reprints it. **Ctrl+C aborts
from anywhere** — it is sent as a raw byte and the firmware honours it mid-move.
Ctrl+] closes the terminal without touching the board.

```
  1  Home + run matrix    homing, then 344 rows (~32 min)
  2  Home only
  3  Run matrix           uses the travel already measured
  4  Go to centre         the start position, on a full step
  5  Driver readback
  6  Switch state
  7  Verify MRES switching    MSCNT per pulse, no travel
  0  Release motor        EN high, outputs off
```

### 3.1 Start the encoder recording first

The board does not see the stage. Start the TwinCAT YT scope recording
**before** pressing `1`, and stop it after `# RUN_COMPLETE`. A recording that
starts mid-run loses the block-1 marker and with it the clock alignment.

### 3.2 Home

Press `2` (or `1`, which homes then runs). Homing seeks each end at 2 mm/s,
backs off 2 mm, and touches the switch **three times** at 0.15 mm/s, requiring
all three triggers within 20 µm:

```
# HOMING,neg,touch=1,pos=-...
# HOMING,neg,touch=2,pos=-...
# HOMING,neg,touch=3,pos=-...
# HOMING,neg,spread=...,tolerance=512
```

A spread over tolerance is retried once, then fails. If it fails repeatedly the
switch is bouncing or the actuator is loose — fix it rather than widening the
tolerance, because the tolerance is what makes the measured travel trustworthy.

Success prints the travel:

```
# HOMED,travel_u256=...,travel_mm=...,lim_neg=...,lim_pos=...
```

**Check `travel_mm` against the stage's spec.** If it is short, a switch
triggered early and every limit below it is wrong.

### 3.3 Run

Press `3`. The firmware prints the extent the plan will use and compares it
against the measured limits:

```
# PLAN_EXTENT,min_mm=-20.000,max_mm=20.000,limits_mm=-...
```

Then it asks for confirmation — type `RUN` and Enter. Anything else cancels.

It then re-runs the MRES check from §2 step 6 and prints `# MRES_CHECK,pass`
before the first row. If that fails the run does not start, and nothing is
lost — fix the driver configuration and press `3` again.

During the run, expect per row:

```
# MOVE,<row>,<tag>,target=...,d=...,dir=...,mres=...
# STOP,<row>,<tag>,d=...,dir=...,measured=1,pos=...,mscnt=...
```

and `# S,` samples at 50 Hz (logged, not printed — they would scroll over
everything). Block boundaries are visible as `# MRES,<size>,mscnt=...` followed
by a marker excursion.

The run ends with one of:

| | |
|---|---|
| `# RUN_COMPLETE,pos=0` | all 344 rows done, back at the start |
| `# RUN_ABORTED,pos=...` | Ctrl+C or `ABORT` |
| `# RUN_FAILED,pos=...,reason=<n>` | see below |

`reason`: 1 abort, 2 switch hit, 3 soft limit, 4 error (off-grid target, pulse
rate, or an MRES write that would not confirm).

### 3.4 Afterwards

1. Stop the scope recording; export the CSV into `data/`.
2. Keep the terminal log — the analysis needs both halves, and the clock
   alignment needs the `# S,` samples that only the log has.
3. Press `0` to release the motor.

---

## 4. If it fails mid-run

Do not restart from row 0 on top of a partial recording. Note the row number
from the last `# STOP,`, keep the partial data with its log, and start a fresh
recording for the retry. The blocks are independent — a run that failed in
block 3 still gives usable 1/16 and 1/32 data, because every measurement is a
difference within a block.

Known failure modes and what they mean:

| line | cause |
|---|---|
| `# ERR,reason=mscnt_off_full_step` | an MRES change was attempted off the full-step grid; the plan should prevent this, so this is a bug or a lost step |
| `# ERR,reason=plan_exceeds_travel` | the stage is shorter than the plan needs; shorten the long distances in `settling_design.py` and regenerate |
| `# ERR,reason=pulse_rate` | a row asks for more than 20 kHz; the design should prevent it |
| `# ERR,reason=driver_mode_not_confirmed` | SpreadCycle/intpol/StealthChop not as intended; do not run |
| `reason=3` (soft limit) | a lost-step accumulation walked the stage out of bounds; the encoder record will show it |

A switch hit (`reason=2`) during the matrix should never happen — the plan stays
2 mm inside the soft limits, which are themselves 2 mm inside the switches. If
it does, the stage has lost a significant number of steps and the block's data
is suspect.
