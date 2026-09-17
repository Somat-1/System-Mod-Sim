# Beckhoff Driven Testing

## Beckhoff Test Implementation

- `COMMISSIONING.md` — TwinCAT import order, EL70x1 PDO/CoE setup, timing limits, dry commissioning, and operating procedure.
- `DUT_BeckhoffTests.st` — shared TwinCAT enumeration types.
- `GVL_BeckhoffIO.st` — variables to link to the EL70x1 Position-control PDOs.
- `MAIN_Tester.st` — top-level selector, safety interlocks, output ownership, and Scope channels.
- `Stepping/` — v5/v6-inspired one-microstep oscillation, reproducible renderer, and `Rendered/` previews.
- `Chirp/` — unnotched 250–1000–250 Hz position chirp, reproducible renderer, and `Rendered/` previews.

The supplied EL70x1 documentation is the configuration reference. Confirm the exact installed terminal revision and compile in the real TwinCAT project before enabling motion.
## Purpose

This folder contains the files, notes, configurations, test procedures, scripts, data, and results associated with Beckhoff-driven testing.

## General Guidance

- Keep all Beckhoff-driven test material within this folder or an appropriately named subfolder.
- Use clear, descriptive filenames and include dates or revision identifiers where they improve traceability.
- Document the test objective, setup, hardware and software configuration, procedure, and acceptance criteria before testing.
- Preserve raw test data. Store processed data and generated plots separately from their source data.
- Record relevant Beckhoff details, including the controller, I/O or motion hardware, TwinCAT version, project/configuration revision, task cycle times, and communication settings.
- Note dependencies, assumptions, known limitations, safety considerations, and required recovery or reset steps.
- Do not overwrite significant results without retaining the earlier version or documenting why it was replaced.
- Update this document whenever the folder structure, workflow, test setup, or important implementation details change.

## Suggested Structure

Create subfolders as needed, for example:

- `Configurations/` — Beckhoff and TwinCAT configuration exports or snapshots.
- `Scripts/` — automation, control, processing, and analysis scripts.
- `Test Procedures/` — repeatable test instructions and checklists.
- `Data/Raw/` — original, unmodified measurements.
- `Data/Processed/` — cleaned or transformed measurements.
- `Results/` — plots, reports, and summarized findings.
- `References/` — manuals, interface definitions, and supporting notes.

## Development Log

Maintain a chronological log of developments in this section. Add an entry whenever a meaningful change is made, including changes to hardware, software, configuration, test methods, scripts, data-processing methods, or conclusions. Do not remove earlier entries; append new entries at the top so the latest development is easiest to find.

Each entry should include:

- Date and author.
- Summary of the change or activity.
- Files or systems affected.
- Reason for the change.
- Verification performed and outcome.
- Open issues or next steps.

### Entry Template

#### YYYY-MM-DD — Author — Short title

- **Summary:**
- **Affected files/systems:**
- **Reason:**
- **Verification/outcome:**
- **Open issues/next steps:**

### 2026-09-17 — Expanded stepping ladder and simplified preview

- **Summary:** Expanded the oscillating stepping test to MRES 1, 2, 4, 8, 16, 32, and 64; removed the per-resolution plot; regenerated the single overview; and latched the chirp microstep grid at run start.
- **Affected files/systems:** `Stepping/FB_OscillatingStepping.st`, `Stepping/render_stepping_sequence.py`, `Stepping/Rendered/`, `Stepping/README.md`, `Chirp/FB_Chirp250To1000.st`, and `Chirp/README.md`.
- **Reason:** Cover the complete power-of-two microstep ladder and keep only the useful full-run visualization.
- **Verification/outcome:** Confirmed a 1,090 s stepping run, increments of 64/32/16/8/4/2/1 terminal counts, seven labeled overview blocks, and a fixed 1/64-step chirp grid.
- **Open issues/next steps:** Compile the revised ST in TwinCAT and confirm the installed terminal is configured for 1/64 microstepping.
### 2026-09-17 — Added rendered run previews

- **Summary:** Added reproducible overview and detail plots showing the configured Beckhoff stepping and chirp runs.
- **Affected files/systems:** `Stepping/render_stepping_sequence.py`, `Stepping/Rendered/`, `Chirp/render_chirp_sequence.py`, `Chirp/Rendered/`, and the related README files.
- **Reason:** Make the expected timing, displacement levels, chirp schedule, and 250 µs cyclic sampling visually reviewable before hardware execution.
- **Verification/outcome:** Generated and visually checked four PNGs. Confirmed the 625 s stepping duration, 207 s chirp duration, MRES block amplitudes, 250–1000–250 Hz schedule, and four samples per 1 kHz period.
- **Open issues/next steps:** Regenerate the plots whenever the ST timing, amplitude, frequency, or MRES parameters change.
### 2026-09-17 — Beckhoff stepping and chirp ST implementation

- **Summary:** Added EL70x1 Position-control Structured Text for the oscillating MRES step test and an unnotched 250–1000–250 Hz chirp, plus a shared top-level tester.
- **Affected files/systems:** `COMMISSIONING.md`, `DUT_BeckhoffTests.st`, `GVL_BeckhoffIO.st`, `MAIN_Tester.st`, `Stepping/`, and `Chirp/`.
- **Reason:** Run the existing ESP32-inspired measurements directly from TwinCAT/Beckhoff without the unwanted long trajectory ramps.
- **Verification/outcome:** Source sequences and the supplied EL70x1 manual were reviewed; file structure, constants, phase equations, speed demand, PDO assignments, and timing guards were checked locally. TwinCAT compilation and hardware validation remain required on the target system.
- **Open issues/next steps:** Confirm the exact EL70x1 terminal revision, link PDOs, compile in TwinCAT, dry-commission the stepping test, and validate whether four samples per period at 1 kHz are acceptable for the intended chirp measurement.
### 2026-09-17 — Initial folder setup

- **Summary:** Created the Beckhoff Driven Testing folder and its general guidance document.
- **Affected files/systems:** `README.md`.
- **Reason:** Establish a dedicated, traceable workspace for Beckhoff-driven testing.
- **Verification/outcome:** Guidance and development-log conventions documented.
- **Open issues/next steps:** Add project-specific setup details and subfolders as testing work begins.