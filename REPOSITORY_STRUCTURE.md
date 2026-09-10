# Repository structure and revision guide

Map of this workspace: what lives where, how the revisions relate, which ones
are current, and which can be safely ignored.

Written 2026-09-10 from a survey of the tree, the per-directory READMEs, and
`git log` recency per area. Where a directory's own README disagrees with the
evidence, that is flagged rather than silently resolved.

---

## 1. The two workstreams

The repository holds two parallel efforts that feed each other. Almost every
navigation mistake comes from not knowing which one a directory belongs to.

| Workstream | Question it answers | Directories |
|---|---|---|
| **Modelling** | What *should* the stage do? | `Rev 2/`, `Rev 3/`, `Rev 4/`, `Rev 4 Analytical Model Derivation/` |
| **Hardware measurement** | What does the stage *actually* do? | `Microstepping Test Data/`, `ESP32S3_TMC2209_Chirp_Test/` |

The physical system throughout is a stepper-driven ball-screw positioning
stage: 2 mm lead, 200 full steps/rev, ~150 mm travel, TMC2209 driver on an
ESP32-S3, measured with a Dewesoft IDS linear encoder and accelerometers.

---

## 2. Quick orientation — where to start

| If you want… | Go to |
|---|---|
| The current model | `Rev 4/state_space_6dof.md` |
| The current friction model | `Rev 4/lugre_friction/Rev 4.2/` |
| The current hardware campaign | `Microstepping Test Data/v4/` |
| Frequency-response measurements | `ESP32S3_TMC2209_Chirp_Test/v2/` |
| Known driver-configuration pitfalls | `TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md` |
| Build/rendering defects | `IMPLEMENTATION_BACKLOG.md` |

> ⚠️ **The root `README.md` is out of date.** It states "Revision 3 is the
> active model" and does not mention Rev 4 at all, but Rev 4 exists, is far
> larger, and was last touched 2026-09-03 against Rev 3's 2026-08-17. Treat
> the root README's *derivation narrative* as still broadly valid and its
> *"which revision is active"* claim as stale.

---

## 3. Modelling lineage

Each revision supersedes the last. They are kept side by side deliberately so
results stay comparable, not because all are live.

### `Rev 2/` — superseded baseline · **archive, do not edit**

First documented revision. The root README is explicit that it is retained as
the previous baseline and *should not be edited* when changing later
revisions. Last touched 2026-08-10.

Contains `Analytical_derivation_and_responses.{md,html}`,
`Simulation_description.{md,html}`, one builder script, and 8 response SVGs.

**Relevance: low.** Historical comparison only. 1.7 MB.

### `Rev 3/` — previous active model · **reference**

The revision the root README still describes. Adds the kinematic diagrams,
flowcharts, `model_parameters.json`, and a much larger rendered asset set.
Last touched 2026-08-17.

Its build entry point is `Rev 3/build_model_documentation.py`, which
regenerates the SVGs and renders both Markdown documents to HTML.

**Relevance: medium.** Superseded by Rev 4 in substance, but it remains the
most *complete* single narrative document set (the Rev 4 material is spread
across sub-approaches). Useful for understanding the derivation as a whole.
5.6 MB.

### `Rev 4/` — current model · **active**

The live modelling work, and the most structurally complex part of the
repository. 107 MB, 223 files. It is not one model but a **frictionless
reference plant plus several competing treatments layered on it**.

```
Rev 4/
├── state_space_6dof.md          ← the frictionless 6-DOF reference plant
├── model_parameters.json        ← parameters for the root model only
├── backlog.md                   ← issues for the frictionless plant only
├── scripts/                     ← build_bode_rev4.py, stepping trajectory, plots
├── Lagrange Derivation/         ← independent frictionless re-derivation
├── Guyan Model Reduction/       ← classical Guyan reduction
│   └── Guyan w Friction/        ← the same reduction with friction added
└── lugre_friction/
    ├── Rev 4.1/                 ← LuGre REPLACEMENT approach
    └── Rev 4.2/                 ← LuGre PARALLEL approach  ← most recent
```

**The root of `Rev 4/` is the common comparison baseline.** Its `backlog.md`
says so directly: it tracks only the structural frictionless plant, and it
"does not select between" the friction approaches. `Rev 4.1` and `Rev 4.2` are
each self-contained — `Rev 4.1/model_parameters.json` is a standalone copy,
not an import — so do not assume a change at the root propagates downward.

#### `Rev 4.1` vs `Rev 4.2` — the distinction that matters most

Both add LuGre friction, and they differ in exactly one structural decision at
the screw–nut interface:

| | Rev 4.1 — *replacement* | Rev 4.2 — *parallel* |
|---|---|---|
| Structural `k_nut` / `c_nut` | **Removed**, LuGre replaces them | **Retained** in the structural matrices |
| LuGre element | Is the load path | An additional pre-rolling drag branch |
| Consequence | LuGre can cap transmitted thrust at breakaway force | It cannot |

`Rev 4.2` is the more recent and carries the latest calibration work (last
commit 2026-09-03, "Rev 4.2 friction recalibration, C-block reconstruction
fix"). `Rev 4.1` also holds a large modal-analysis asset set (Co-MAC, mode
tracking, pole maps) that has no Rev 4.2 equivalent.

**Relevance: Rev 4.2 high (current), Rev 4.1 medium (alternative hypothesis
plus unique modal analysis), Lagrange/Guyan medium (independent cross-checks
of the same frictionless plant).**

### `Rev 4 Analytical Model Derivation/` — **stub, low relevance**

Four files: an `index.html` and three assets, two of which are named
`*_placeholder.*`. Appears to be a started-but-unfinished presentation layer
for the Rev 4 derivation. Not referenced by the root README.

**Relevance: low.** Check whether it is still wanted before investing in it.

---

## 4. Hardware measurement lineage

### `Microstepping Test Data/` — step-level positional accuracy

117 MB. Its own README describes only v1 and v2 and is therefore also stale;
v3 and v4 exist and are more recent.

| Version | What it is | Last touched | Relevance |
|---|---|---|---|
| `v1/` | Migrated legacy measurements — `StepSize{1,2,16}.csv`, `IDSdata.txt`, plus a `temp_model_optimization/` scratch area. No README. | 2026-08-31 | **Low** — archival |
| `v2/` | Motion-sequence *specification* and command-only preview figures, plus the ESP32 identification runner for Blocks A/B/E. `BACKLOG.md` is the authoritative sequence rationale. | 2026-08-31 | **Medium** — the sequence design rationale is still cited |
| `v3/` | The EVO dedicated-controller reliability campaign: six MRES/current configurations through reference, conditioning, creep/settling (C) and velocity-plateau (D) blocks. Real hardware data. | 2026-09-04 | **High** — source of the plateau rates v4 reuses |
| `v4/` | **Current campaign.** MRES ∈ {1,4,16,32} oscillation + 25 mm trajectory, ESP32-S3 firmware, direct vs individual command generation. | 2026-09-10 | **Highest — active work** |

Note `v3/data/raw_local/SteppingSequenceID.csv` is `.gitignore`d at >50 MB and
is **not regenerable** — it is a direct Dewesoft export. Keep a backup outside
git.

Inside `v4/`, the important distinction:

- `scripts/esp32_v4_mres_trajectory_campaign/` — **the live firmware**, the one
  that actually runs the campaign.
- `scripts/run_mres_trajectory_campaign.py` — a host-side *dry-run mirror* used
  to validate block count, order, timing and CSV schema. It does not drive the
  ESP32.
- `scripts/run_settling_dedicated_controller.py` — retained legacy v3-era
  settling experiment. **Not** the current campaign.
- `ESP32_TMC2209_StepSize_Sweep/` — a separate 7-value MRES sweep sketch with
  StallGuard/DRV_STATUS logging. Its README states it was never compiled or
  flashed. Independent of the main campaign.
- `scripts/tmp_v4_oscillation_check/` — bench diagnostic added 2026-09-10 while
  investigating the interpolation defect. Not part of the campaign.
- `STEALTHCHOP_STALLGUARD_REPEAT_MEMO.md` — specification for a *future* repeat
  with StealthChop and StallGuard enabled. Not yet executed.

### `ESP32S3_TMC2209_Chirp_Test/` — frequency-response measurement

921 MB, the largest area by far, almost entirely raw accelerometer captures.

| Version | What it is | Relevance |
|---|---|---|
| root (`v1`) | Bounded chirp, 1 Hz → 1 kHz → 1 Hz over 64.5 s, alternating ±1 microstep so commanded position never leaves the origin. | **Medium** — superseded by v2 |
| `v2/` | Resonance-band chirp with the hard notch replaced by a smooth taper. Bode analysis added 2026-09-10. | **High — most recent** |

> ⚠️ **`v2/README.md` is stale.** It states "nothing here has run on real
> hardware" and describes the firmware as an unflashed draft. But
> `.gitignore` documents `v2/Ftest.csv` as a 744 MB, 10.9-million-row capture
> of 3 accelerometer axes at 20 kHz, and recent commits add v2 chirp Bode
> analysis. The v2 campaign has clearly run. The status section needs updating.

Both chirp raw captures are local-only and **not regenerable**:
`tempTVAl.csv` (~2.9 GB, v1) and `v2/Ftest.csv` (744 MB). Back them up outside
git. The `v2/analysis/cache/` `.npy` parse is derived and rebuilds in ~60 s.

---

## 5. Root-level files

| File | Purpose | Relevance |
|---|---|---|
| `README.md` | Workspace overview and derivation narrative | High, but see the staleness warning in §2 |
| `IMPLEMENTATION_BACKLOG.md` | Rendering defects, build failures, unexpected behaviour | Medium |
| `TMC2209_DRIVER_CONFIGURATION_BACKGROUND.md` | Driver register pitfalls: interpolation default, UART read hazards, sense-resistor question | High for any hardware run |
| `REPOSITORY_STRUCTURE.md` | This file | — |

---

## 6. What to disregard

### Scratch and superseded areas

| Path | Why |
|---|---|
| `temp/` (root) | Three files — one script, two PNGs, on a notch comparison. No README, not referenced anywhere. Scratch. |
| `Microstepping Test Data/v1/temp_model_optimization/` | Superseded optimisation scratch |
| `Microstepping Test Data/v2/temp_preemptive_run_analysis/` | Pre-run analysis, superseded by the actual run |
| `Microstepping Test Data/v2/temp_step_size_determination/` | Folded into the v2 sequence spec |
| `Microstepping Test Data/v2/temp_testdiag_analysis/` | One-off diagnostic |
| `Rev 4/rendered_assets/temp/` | Intermediate render output |
| `Rev 4/lugre_friction/Rev 4.1/rendered_assets/temp/` | Intermediate Co-MAC renders |

The `temp_`-prefix convention is used consistently for throwaway analysis, and
that convention is worth keeping.

### Committed build artefacts and OS cruft

Roughly **70 files that should not be tracked** are currently in git:

- **`__pycache__/*.pyc`** — ~60 files across `Rev 2`, `Rev 3`, `Rev 4` and its
  sub-approaches, and `Microstepping Test Data/v1`–`v2`. Several are duplicated
  across Python versions (`cpython-312` *and* `cpython-314`).
- **`.DS_Store`** — `Rev 4/`, `Rev 4/lugre_friction/`
- **`Thumbs.db`** — 6 copies under `Rev 4/rendered_assets/` and Rev 4.1/4.2
- **`desktop.ini`** — `ESP32S3_TMC2209_Chirp_Test/v2/`
- **`esp32_tmc2209_stepsize_sweep.ino.cpp`** — an Arduino build artefact,
  currently untracked but sitting in the working tree

Suggested `.gitignore` additions:

```gitignore
__pycache__/
*.pyc
.DS_Store
Thumbs.db
desktop.ini
*.ino.cpp
```

Removing the already-tracked ones needs `git rm --cached` as a separate,
deliberate commit — listed here as a recommendation, not something done.

---

## 7. Conventions that hold across the tree

- **`rendered_assets/`** — generated output. Never hand-edit; regenerate with
  the owning script.
- **`scripts/`** — the generators. Edit these, not their output.
- **Markdown is source, HTML is generated.** Edit the `.md`; the builder
  renders the `.html`.
- **Each revision resolves paths relative to its own directory**, so scripts
  are not portable between revisions without adjustment.
- **Git LFS** tracks `*.npz`, `*.csv`, `*.png`, `*.mat`.
- **Large raw hardware captures are deliberately local-only.** Three files are
  not regenerable and have no backup inside git — see §4.

---

## 8. Maintenance notes

Three documentation drifts were found during this survey and are recorded here
rather than fixed, since correcting them changes documents owned by other
workstreams:

1. Root `README.md` names Rev 3 as active; Rev 4 is active and unmentioned.
2. `Microstepping Test Data/README.md` describes only v1 and v2; v3 and v4
   exist and v4 is current.
3. `ESP32S3_TMC2209_Chirp_Test/v2/README.md` says nothing has run on hardware;
   a 744 MB capture and a Bode analysis say otherwise.
