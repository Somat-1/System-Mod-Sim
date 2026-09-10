# TMC2209 driver configuration: interpolation, register writes, and what they cost us

Background note written 2026-09-10, after the v4 MRES campaign was found to be
measuring something other than what it intended to measure.

This note exists because the defect was invisible in the firmware source. The
v4 campaign sketch contained no bug in the usual sense: no wrong constant, no
inverted logic, no bad arithmetic. It was wrong because of something it *did
not* write.

---

## 1. Summary

Two independent problems were found in the v4 campaign firmware.

| # | Problem | Effect | Status |
|---|---|---|---|
| 1 | MicroPlyer interpolation left at its power-on default (enabled) | Every commanded microstep was smeared over the interval to the next step instead of moving discretely | Fixed: `intpol(false)` forced and verified |
| 2 | `setMres()` could not distinguish a failed UART read from a real register value | Aborted a 41-minute run; could also have silently disabled the driver | Fixed: reads validated, write/verify retried |

**Only the v4 campaign was affected by problem 1.** Every other sketch in this
repository already disables interpolation explicitly. The v2 chirp test is
clean — see [section 5](#5-which-experiments-were-affected).

**Problem 2 is a shared pattern** and the same latent hazard exists in the
other sketches, including the chirp tests. See
[section 7](#7-the-cross-cutting-hazard-the-chopconf-read-modify-write).

---

## 2. What MicroPlyer interpolation does

The TMC2209 contains a step interpolator ("MicroPlyer"). When enabled, it takes
each incoming STEP pulse and expands it into 256 internal microsteps, which it
plays out gradually over an interval derived from the time between the previous
step pulses.

The intent is noise and smoothness: a printer commanding coarse steps gets
256-microstep smoothness for free. For that use it is a good feature.

For an experiment that measures *the physical response to one commanded
microstep*, it is fatal. The thing being measured — a discrete step and the
mechanical settling that follows it — is precisely what the interpolator
removes.

### The default is ON

This is the part that caught us. Per the TMC2209 datasheet, `CHOPCONF` resets
to `0x10000053`. Bit 28 (`intpol`) is **set** in that reset value.

So "we do not configure interpolation" does not mean interpolation is off. It
means interpolation is **on**. The v4 firmware's design intent — configure as
little as possible to keep the run a clean baseline — is exactly what caused
the problem. Leaving a register alone is a choice with a value attached, and
here the value was the wrong one.

---

## 3. How it was found

The firmware's own position accounting could not detect this. It counted the
pulses it emitted, and it emitted them correctly. Both the commanded position
and the block timings were right. From the serial log the run looked perfect.

The diagnosis came from asking the driver what *it* thought its position was,
using two registers the campaign never read:

- **`MSCNT`** (0x6A) — the driver's own index into its 1024-entry sine table.
  That table spans four full steps, so one commanded microstep should advance
  `MSCNT` by exactly `256 / MRES`.
- **`TSTEP`** (0x12) — the measured interval between step pulses, which
  saturates at `0xFFFFF` when no steps arrive.

`MSCNT` is the ground truth. It is incremented by the driver's own step edge
detection, so comparing it against what the firmware believes it commanded
separates "the pulse never arrived" from "the pulse arrived and the driver did
something other than what was expected".

### The measurement

At MRES=1, a commanded full step should move `MSCNT` by 256. Reading `MSCNT`
immediately before and immediately after the step:

```
OSC_FORWARD  mres=1  mscnt 640 -> 627   delta = -13   expected 256
```

Only 13 counts. But reading again one second later, before the return step, it
had reached 384 — a total of 256 counts from 640. The step was not lost. It was
*arriving late*, dribbling out over the following second.

That signature — correct total, wrong distribution in time — is interpolation.

### Confirmation

Setting `driver.intpol(false)` and changing nothing else:

| MRES | With `intpol=1` (default) | With `intpol=0` | Expected |
|---|---|---|---|
| 1 | 640 → 627, delta **−13** | 640 → 384, delta **−256** | 256 |
| 4 | 736 → 734, delta **−2** | delta **−64** | 64 |
| 16 | 760 → 759, delta **−1** | delta **−16** | 16 |
| 32 | 756 → 755, delta **−1** | delta **−8** | 8 |

With interpolation off, every commanded microstep moves the driver's position
by exactly the right amount, immediately, at every resolution tested.

> The sign is negative because a commanded "positive" direction drives DIR high,
> which decreases `MSCNT` on this rig. Magnitudes are exact and every sequence
> returns to its origin, so this is a sign convention, not an error.

---

## 4. What this did to the v4 data

The oscillation blocks are the clearest casualty. Their whole purpose is fifteen
discrete one-microstep out-and-back cycles. With interpolation on, each
"microstep" was a slow creep spread across the one-second dwell that was
supposed to be the *settling observation window*. The motion and the
measurement window overlapped, so what was recorded was the interpolator's
output profile rather than the mechanical response to a step input.

**This was not confined to the oscillation blocks.** Interpolation applies to
every step pulse the driver receives, so the 25 mm trajectories ran with it too.
At the trajectory rates the effect is far less visible — steps arrive fast
enough that the interpolation interval is short and the smearing is small
relative to the motion — but it is not zero, and it is not characterised.

**Practical consequence:** v4 data recorded before 2026-09-10 is not directly
comparable with data recorded after. The oscillation blocks should be
considered invalid for step-response purposes. The trajectory blocks are more
salvageable but carry an uncharacterised low-pass effect at the step level.

---

## 5. Which experiments were affected

Every sketch in the repository was audited for the interpolation call.

| Sketch | Disables `intpol`? | Affected? |
|---|---|---|
| `ESP32S3_TMC2209_Chirp_Test.ino` | Yes — `intpol(false)`, line 113 | **No** |
| `v2/scripts/esp32_v2_chirp_test.ino` | Yes — `intpol(false)`, line 472 | **No** |
| `Microstepping Test Data/v2/.../run_identification_esp32_tmc2209.ino` | Yes — lines 239 and 283 | **No** |
| `v4/ESP32_TMC2209_StepSize_Sweep/.../esp32_tmc2209_stepsize_sweep.ino` | Yes — line 237 | **No** |
| `v4/scripts/esp32_v4_mres_trajectory_campaign.ino` | **No** (until 2026-09-10) | **Yes** |

### The chirp tests are clean

Both chirp sketches configure the chopper explicitly and deliberately, and
disable interpolation as part of that. From `esp32_v2_chirp_test.ino`:

```cpp
driver.en_spreadCycle(true);   // force SpreadCycle for the whole run
driver.TPWMTHRS(0);
driver.ihold(driver.irun());   // standstill current must not drop
driver.iholddelay(0);
driver.TPOWERDOWN(0);
driver.intpol(false);          // <-- interpolation explicitly off
```

The ordering is also correct: `intpol(false)` is issued *before* the
`CHOPCONF` read-modify-write that sets MRES, and that read-modify-write
preserves bit 28, so nothing re-enables it afterwards.

**Conclusion: the v2 chirp test measurements are not contaminated by this
defect.** The frequency-response data from those runs stands.

### Why v4 was the odd one out

The chirp and identification sketches were written to pin the driver into one
known chopper configuration, because a mode switch mid-sweep would corrupt a
frequency response. Disabling interpolation fell naturally out of that intent.

The v4 campaign was written with the opposite intent — "baseline motion only,
change nothing" — to avoid biasing a comparison between microstep resolutions.
That intent silently inherited a default that was wrong for the measurement.

The lesson generalises: **for a measurement rig, "leave it at default" is not a
neutral choice and should never be stated as one.** Every register that affects
the quantity being measured needs to be set deliberately and verified by
readback, including the ones being set to their default value.

---

## 6. The fix

In `configureDriver()`:

```cpp
driver.intpol(false);  // Powers up enabled; see header note.
const uint32_t chopconf = driver.CHOPCONF();
if (((chopconf >> 28) & 1U) != 0U) {
  Serial.printf("# INTPOL_DISABLE_FAILED,chopconf=0x%08lX\n", ...);
  return false;   // refuses to start the campaign
}
```

The verification matters as much as the write. Returning `false` here makes
`configureDriver()` fail, which drives the firmware into its failure-blink state
instead of starting a campaign. A run that reaches `CAMPAIGN_START` is therefore
*proof* that interpolation was off, rather than an assumption that it was.

---

## 7. The cross-cutting hazard: the CHOPCONF read-modify-write

A second, independent problem surfaced when a 41-minute run aborted ten minutes
in with:

```
# MRES_READBACK_FAILED,wrote=6,read=0
```

### Root cause

`TMCStepper`'s UART read retries twice (`max_retries = 2`) and then **returns 0
while setting a public `CRCerror` flag**. A caller that ignores that flag cannot
distinguish a dead transaction from a register that genuinely reads zero.

Here the readback of `CHOPCONF` failed, returned 0, and 0 was interpreted as
"MRES code 0" — a mismatch against the requested code 6, which aborted the run.
Nothing was wrong with the driver; one serial datagram was corrupt.

### The more dangerous variant

The same ambiguity on the *initial* read of a read-modify-write is worse:

```cpp
uint32_t chopconf = driver.CHOPCONF();          // returns 0 if the read fails
driver.CHOPCONF((chopconf & ~MRES_MASK) | (code << 24));
```

If that read fails, every other field in `CHOPCONF` is written as zero —
including `toff`. **`toff = 0` disables the driver entirely.** The subsequent
MRES readback would still return the requested code, so the write appears to
succeed. The result is a campaign that runs to completion, logs correct
positions, and steps nothing at all.

That failure mode is worth internalising, because from the serial log it is
indistinguishable from a successful run.

### The fix in v4

```cpp
bool readChopconfChecked(uint32_t &value) {
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    const uint32_t candidate = driver.CHOPCONF();
    // CRC clean, and toff non-zero -- configureDriver always sets toff=5,
    // so a zero there means the read failed, not that the field is zero.
    if (!driver.CRCerror && (candidate & 0x0FUL) != 0UL) {
      value = candidate;
      return true;
    }
    delay(5);
  }
  return false;
}
```

`setMres()` now retries the whole write/verify cycle up to five times and logs
`MRES_READ_RETRY` / `MRES_MISMATCH_RETRY`, so a marginal UART link becomes
visible in the data instead of killing a long run.

### This pattern is still present elsewhere

**The other sketches use the same unguarded read-modify-write** and carry the
same silent-disable hazard. In `esp32_v2_chirp_test.ino`:

```cpp
uint32_t chopconf = driver.CHOPCONF();          // unchecked
chopconf = (chopconf & ~MRES_MASK) | (MRES_REGISTER_CODE << 24);
driver.CHOPCONF(chopconf);
```

Interpolation is not at risk there (a zeroed write leaves `intpol` at 0, which
is what those sketches want anyway), but `toff` is. A single corrupt read during
startup would disable the driver while the MRES readback still passed.

This has not been observed on the chirp rig, and a failure at *startup* is far
more likely to be noticed than one mid-run, since no motion at all would result.
It is recorded here as a known risk rather than a fixed defect.

---

## 8. Open question: the sense resistor

Every sketch in the repository declares:

```cpp
constexpr float R_SENSE_OHM = 0.03F;
```

`TMCStepper`'s own default is `0.11`, which is the value used by SilentStepStick,
BTT and Fysetc TMC2209 modules. 0.03 Ω is unusual.

This constant is a *declaration of a hardware fact*, not a tuning parameter. It
only feeds the current-scale calculation:

```cpp
CS = 32 * 1.41421 * mA/1000 * (Rsense + 0.02) / 0.325 - 1;
if (CS < 16) { vsense(true); CS = ... / 0.180 - 1; }
```

With `Rsense = 0.03` and a 360 mA request this yields `vsense = true, CS = 3`,
which matches the `cs_actual = 3` read back from the driver — so the model is
confirmed against the hardware.

The consequence depends on the physical resistors:

| Physical Rsense | Actual RMS current | Against 360 mA intent |
|---|---|---|
| 0.03 Ω | ~318 mA | Correct (quantisation only) |
| 0.11 Ω | **~122 mA** | Roughly one third |

If the modules are actually 0.11 Ω, every run in this repository has been
executing at about a third of its intended current, with correspondingly reduced
holding torque and margin against skipped steps.

Because the constant is identical across all sketches, any error is *systematic*
rather than differential, so comparisons between experiments remain internally
valid. Absolute current figures would not be.

**To resolve:** read the markings on the sense resistors — typically `R030` for
0.03 Ω, `R110` for 0.11 Ω. Until then, treat stated current values as nominal.

---

## 9. Verification recipe for future runs

Any change to driver configuration on this rig should be checked against
`MSCNT` rather than against the firmware's own position accounting, because the
firmware cannot detect a driver that ignores or reshapes its pulses.

1. Set the MRES under test.
2. Read `MSCNT`.
3. Command **ten** microsteps in one direction, reading `MSCNT` after each.
   Ten, not one — a single step can be masked by quantisation, and a consistent
   arithmetic sequence is much stronger evidence than one delta.
4. Each reading should differ from the previous by exactly `256 / MRES`.
5. Command ten in the reverse direction and confirm the original value returns.
6. Check `TSTEP` is not `0xFFFFF`, which would mean no step edges arrived.

`tmp_v4_oscillation_check.ino` implements exactly this as a bench tool, along
with a static STEP/DIR hold mode for probing the lines with a multimeter and an
`IOIN`-based check of what the driver sees on its own input pins.

> One caveat on that tool, recorded honestly: during diagnosis its `IOIN`,
> `MSCNT` and `TSTEP` reads all indicated that no step pulses were reaching the
> driver, which led to a wiring investigation that turned out to be a dead end —
> the lines measured correct at both ends and the real cause was interpolation.
> Those registers later read correctly on the same hardware with no wiring
> change. The reason for the earlier readings is not understood. Treat a
> wholesale "no steps arriving" indication with suspicion and confirm it with a
> meter before pulling wiring apart.

---

## 10. References

- TMC2209 datasheet, `CHOPCONF` (0x6C) — reset value `0x10000053`, bit 28
  `intpol`, bits 27:24 `MRES`, bits 3:0 `toff`
- TMC2209 datasheet, `MSCNT` (0x6A) — microstep table counter, 0..1023
- TMC2209 datasheet, `TSTEP` (0x12) — measured step interval, saturates at
  `0xFFFFF`
- `TMCStepper/src/TMCStepper.h` — `max_retries = 2`, `replyDelay = 2`,
  public `CRCerror`, default `Rsense = 0.11`
- `TMCStepper/src/TMCStepper.cpp` — `rms_current()` and the `vsense` selection
