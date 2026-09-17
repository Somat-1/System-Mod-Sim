# TwinCAT Commissioning and Test Procedure

## Scope and hardware assumption

These files target the Beckhoff EL70x1 family described by `el70x1_en.pdf`, using direct **Position control** rather than an NC axis or the Positioning Interface. Confirm the exact terminal type and revision in the live EtherCAT configuration before running. If the terminal does not expose `STM Position` (PDO `0x1603`) and `STM Internal position` (`0x1A07`), do not run these files unchanged.

The code was prepared from the supplied documentation but cannot be compiled or hardware-tested from this workspace. Compile it in the actual TwinCAT project and resolve any project-version-specific type or library differences before enabling power.

## Required TwinCAT objects

Create/import these items in this order:

1. DUT `DUT_BeckhoffTests.st`.
2. GVL named exactly `GVL_BeckhoffIO` from `GVL_BeckhoffIO.st`.
3. Function block `Stepping/FB_OscillatingStepping.st`.
4. Function block `Chirp/FB_Chirp250To1000.st`.
5. Program `MAIN_Tester.st`.
6. Call `MAIN_Tester` from the task that exchanges the EL70x1 PDOs.

The function blocks use `R_TRIG` and `TON` from the standard TwinCAT PLC library.

## Terminal configuration

Set and verify the motor parameters against the physical motor nameplate and wiring. Do not copy current or voltage values blindly.

Important CoE/PDO settings from the supplied EL70x1 manual:

- `0x8010:01` — maximum motor current.
- `0x8010:02` — reduced current.
- `0x8010:03` — nominal motor voltage.
- `0x8010:04` — motor coil resistance.
- `0x8010:06` — motor full steps per revolution (normally 200 here).
- `0x8012:01` — operating mode; select the mode required by the terminal for position control.
- `0x8012:05` — speed range. The chirp calculations assume setting 5, or 32,000 full steps/s.
- EL7041-1000 only: `0x8012:45` — microstepping. The code assumes setting 6, or 1/64 step.

The manual warns against changing CoE settings while the axis is active. Disable the output stage before changing them and retain the settings in the TwinCAT startup list.

Select the predefined **Position control** PDO assignment, or make the equivalent assignment manually:

- SM2: `0x1602 STM Control` and `0x1603 STM Position`.
- SM3: `0x1A03 STM Status` and `0x1A07 STM Internal position`.

Link every variable in `GVL_BeckhoffIO` to the matching PDO subindex stated in its comment. Do not rely on declaration order or automatic structure packing.

## Task timing

- Stepping test: a 1 ms task is sufficient because the endpoints dwell for 5 s.
- Chirp test: use a Distributed Clocks synchronized 250 us EtherCAT/PLC task and set `udiTaskCycleUs := 250` in `MAIN_Tester`.
- Verify that the task has no sustained overruns and that the terminal reports no sync error.

The EL70x1 manual states that its internal cycle is 250 us. At a 250 us command cycle, the 1 kHz endpoint receives only four position samples per sine period. This is enough to produce a coarse sampled command but not enough for a high-fidelity 1 kHz sinusoidal excitation. The chirp block therefore requires `bAcknowledgeMarginalChirpSampling := TRUE` before it will start.

For ten command samples per period, the practical maximum at 250 us is 400 Hz. Use faster/specialized waveform hardware if the full 1 kHz band must be quantitatively accurate.

## Dry commissioning

1. Mechanically decouple the stage or put it in a safe test condition.
2. Keep `bArmDrive := FALSE`.
3. Activate the configuration and run the PLC.
4. Confirm all PDO links and confirm `udiStmInternalPosition` updates.
5. Pulse `bResetDrive` for one scan if a known, investigated terminal error needs resetting.
6. Set `bArmDrive := TRUE`; confirm `bStmReady` and no `bStmError`/`bStmSyncError`.
7. Run the stepping test first and verify direction, one-count response, current, temperature, and return to origin.
8. Test a reduced chirp range/amplitude before the full requested sequence.
9. Record the external position sensor, `diSelectedCommand`, `diActualPosition`, `lrCommandFrequencyHz`, `uiActiveSegment`, and terminal status bits in TwinCAT Scope.

## Operation

Expose the operator variables at the top of `MAIN_Tester` in a visualization or watch list:

- `eSelectedTest`: `TEST_IDLE`, `TEST_STEPPING`, or `TEST_CHIRP`.
- `bArmDrive`: enables the terminal output stage.
- `bStartTest`: pulse for one scan to start the selected test.
- `bAbortAll`: freezes the target at the reported internal position.
- `bResetDrive`: pulse only after investigating a terminal error.
- `bAcknowledgeMarginalChirpSampling`: required for the 1 kHz chirp.

Normal completion returns the command to the captured origin. Abort deliberately freezes at the actual reported internal position instead of commanding an automatic recovery move.

## Chirp feasibility

At the requested 3-full-step amplitude and 1 kHz endpoint, the ideal sine requires approximately 18,850 full steps/s peak velocity. This is below the documented 32,000 full steps/s terminal range assumed by the code. That establishes only a command-rate margin; it does not prove that the motor/load can supply the required acceleration or remain synchronized.

The limiting issue is command resolution: four samples per period at 1 kHz. Therefore:

- Oscillating stepping: suitable for this EL70x1 approach.
- 250–1000 Hz chirp: nominally executable at 250 us and within the selected speed range, but marginal and not sufficient for a clean measurement-grade 1 kHz waveform.
- A 1 ms task is not usable for this chirp; 1 kHz would be sampled once per period and can alias to little or no commanded motion.

Because this is open-loop stepper control, neither the terminal's internal position nor its target position proves that the rotor followed. Use the independent position sensor and abort on the first sign of pole slip or origin drift.