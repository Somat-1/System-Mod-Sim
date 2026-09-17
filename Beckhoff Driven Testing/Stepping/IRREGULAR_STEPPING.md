# Irregular 20 mm Stepping Test

## Purpose

`FB_IrregularStepping20mm.st` commands reproducible irregular target positions within a **20 mm total centered span**, interpreted as −10 mm to +10 mm from the position captured when the test starts.

The test alternates positive and negative targets so every random move crosses back and forth through the centered range. It executes 16 seeded irregular moves by default, waits for the EL70x1 internal position to reach each target, dwells for 5 seconds, and finishes with a logged return to the captured origin.

## Running from `MAIN_Tester`

1. Run the normal commissioning and confirm the physical stage has at least 10 mm safe travel on both sides of its current position.
2. Set `udiIrregularSeed` to the desired nonzero seed. Reusing a seed reproduces the same sequence.
3. Set `uiIrregularRandomMoves` from 1 to 30; the default is 16.
4. Set `eSelectedTest := TEST_IRREGULAR_STEPPING`.
5. Arm the drive and pulse `bStartTest` for one PLC scan.
6. Monitor `fbIrregular.bBusy`, `fbIrregular.bError`, and the external position sensor.

A normal run performs `uiIrregularRandomMoves + 1` commands because the final command returns to origin.

## Command record

The function block records the quantized commands actually sent to the terminal:

- `fbIrregular.aCommandedTargetOffsetMm` — target relative to the captured origin.
- `fbIrregular.aCommandedMoveDistanceMm` — signed distance from the previous target.
- `fbIrregular.aCommandTimeS` — elapsed test time when the command was issued.
- `fbIrregular.aCommandedPositionCounts` — absolute EL70x1 position command.
- `fbIrregular.uiLogCount` — number of valid records.

The arrays remain populated after completion until the next run starts. Record them with TwinCAT Scope/ADS or copy/export them from the watch window. The external sensor remains authoritative because the terminal position is open-loop.

## Preview and CSV

`render_irregular_stepping.py` reproduces the same xorshift32 sequence and writes both a plot and a CSV of the planned commands. With the default seed:

```powershell
$env:MPLBACKEND='Agg'
python .\render_irregular_stepping.py
```

For a different run:

```powershell
python .\render_irregular_stepping.py --seed 0x12345678 --moves 20
```

The preview command times assume instantaneous arrival followed by the 5 s dwell. Actual PLC timestamps include motion-to-target time and are therefore the authoritative timing record.

## Safety interpretation

The default 20 mm range means **20 mm total**, not ±20 mm. If ±20 mm was intended, do not simply increase the value until the physical limits, centered origin, cabling, and stopping margin have been verified.