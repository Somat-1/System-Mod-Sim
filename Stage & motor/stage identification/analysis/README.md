# Analysis of the 13 October 2026 recordings

Processing for the three hardware exports taken on 13 October 2026, and the
report built from them.

```
analysis/
├── cache_recordings.py   raw CSV -> npz caches (run this first)
├── stage_analysis.py     IDdata.csv      -> stage_results.json + 12 figures
├── chirp_analysis.py     tempTvalm*.csv  -> chirp_results.json + 5 figures
├── build_report.py       the two json files + the figures -> the PDF
├── vizstyle.py           shared figure style
├── plots/                one PNG per figure
├── cache/                npz caches, gitignored, rebuilt in ~45 s
└── stage_identification_report.pdf
```

Run in order; the whole chain takes about four minutes, most of it in the
chirp response extraction:

```powershell
python cache_recordings.py
python stage_analysis.py
python chirp_analysis.py
python build_report.py
```

Needs numpy, scipy and matplotlib. Nothing else, and no network.

## The recordings

| File | What it is | Length |
|---|---|---|
| `../IDdata.csv` | TwinCAT YT Scope export, one channel: the EL5101 counter at 1 kHz | 44.9 min, 45 MB |
| `../chirp/tempTvalm.csv` | Dewesoft export, 3 accelerometer axes at 20 kHz | 4.6 min, 375 MB |
| `../chirp/tempTvalm2.csv` | same format | 10.7 min, 865 MB |

All three are raw hardware exports and are **not regenerable**. They are
local-only by the same policy as the other captures in this repository — keep
backups outside git.

## What the processing has to do first

Neither dataset carries the information that would normally tie it to the
plan, so both pipelines start by recovering it from the signal.

**The stage recording has no command channel and no sync marker.** It is the
encoder counter and nothing else. A fixed time offset does not align it: the
plan's timeline models a move as distance over speed with no acceleration, so
the rig falls behind, and after 230 moves the error is minutes. Instead every
interval in which the stage moved is found by thresholding a smoothed
`|dcounts|`, and that list is aligned to the planned move list by dynamic
programming over commanded displacement, with gaps allowed on both sides. The
alignment is judged by what it skips:

- 222 of 230 planned moves matched, displacement agreeing to 0.26 full steps rms;
- the 8 unmatched planned moves are exactly the 6 quasi-static staircase legs,
  whose microsteps are too small to trip the motion threshold, plus 2
  zero-displacement rows;
- the 8 unmatched measured segments are the operator's 2 mm jog before the run,
  plus 7 noise blips during dwells.

The counts-per-full-step scale is then fitted on phase 1 rather than assumed,
because phase 1 exists to measure it. The staircase legs are read back on the
recovered 0.751 s microstep cadence rather than by detecting each step, so that
microsteps during which the stage did not visibly move are still counted —
those are the interesting ones.

**The chirp recordings do not say which profile they hold.** Peak-picking the
spectrogram does not find out either: the rig has a loud band near 685 Hz that
outweighs the drive wherever the sweep is not passing through it, so a tracker
follows the resonance. Each planned sweep leg is instead slid along the
recording as a line of known start frequency, end frequency and duration, and
scored against the spectrogram under that line after each bin has had its own
spectral neighbourhood subtracted. A profile counts as present only when every
leg it specifies is found — which is what separates FULL from HIGH, since the
two share their 250–1000 Hz leg.

## What was found in the files

Both accelerometer files hold the **HIGH** profile only: three complete runs,
and no sweep below 250 Hz anywhere. **The FULL profile was not recorded.**

| File | Runs | Windows |
|---|---|---|
| `tempTvalm.csv` | 1 × HIGH | 79–273 s |
| `tempTvalm2.csv` | 2 × HIGH | 58–254 s and 394–590 s |

The stage recording holds the full 230-move sequence, starting 115 s in.

## Headline results

| Quantity | Result |
|---|---|
| Step to distance gain | 1 full step = 9.9973 µm (−0.027 % of nominal) |
| Unidirectional repeatability | R↑ 1.21 µm, R↓ 1.25 µm |
| Bidirectional accuracy | A = 8.69 µm |
| Reversal value (backlash) | −1.71 µm, near-constant along the travel |
| Static microstep error | up to 1.41 µm peak-to-valley |
| Post-move settling | none resolvable above 0.2 µm |
| Dominant resonance | 685 Hz, all three axes, both sweep directions |

The 0.26 % counts-per-revolution anomaly recorded in
[`../../CONNECTIONS.md`](../../CONNECTIONS.md) **is not reproduced**: this
measurement puts the gain 0.027 % the other way.

Full numbers, figures and the caveats that go with them are in
[`stage_identification_report.pdf`](stage_identification_report.pdf). The two
json files hold everything the report draws on, so a number in the PDF can
always be traced back to the run that produced it.

## Artefacts

Three were found and corrected in processing; they do not affect the numbers
above. The report page *Artefacts found, and what was done about each* has the
detail.

| Artefact | Correction |
|---|---|
| A 252 Hz peak on one down sweep | the 3 s dwell at 250 Hz was leaking into the first bin of the leg; half a second is now trimmed from each leg end |
| Drift check 1 reading +31.6 µm | its dwell was being measured through the staircase leg that follows; the dwell window is now capped at the planned length |
| Inflated staircase peak-to-valley | the apex slot belongs to both passes at once; both now stop short of it |

Three are properties of how the runs were taken and cannot be processed away:

- **The −5 mm staircase leg is entered against its own direction.** The
  positioning move arrives travelling negative and the leg then runs positive,
  so that leg alone spends its first microsteps taking up lost motion. Its
  hysteresis figure (+0.41 µm) is not comparable with the other two (−2.34,
  −1.89 µm). The staircase needs a run-up, as the repeatability phase has.
- **Seven impulsive disturbances of 0.21–0.36 µm** appear during settling
  dwells, in two clusters near 1509–1512 s and 1678–1698 s. They return to
  where they started and are ~20× the 15 nm encoder noise, so they are
  environmental rather than lost steps or sensor noise. Settled positions are
  read as a median over the dwell, so they do not bias these results — but they
  are a quarter of the repeatability figure.
- **The two chirp files came from different configurations.** The 750 Hz sync
  markers are specified identically in every run and read 0.0110 m/s² in
  `tempTvalm` against 0.0918 and 0.0894 in `tempTvalm2`; the turnaround dwell
  is 2.00 s in the first and 4.56 s in both runs of the second, against a
  planned 2.00 s. Amplitudes are comparable within `tempTvalm2` and not across
  the files, and `tempTvalm2` did not come from the profile `chirp_design.py`
  currently generates.

## Other caveats before quoting any of this

- **The FULL chirp profile is missing**, so nothing here describes the stage
  below 250 Hz, and the HIGH-against-FULL comparison the campaign was designed
  around cannot be made from these files.
- **The encoder scale is one session's measurement**, not a calibration, and
  the scope channel's own scaling has not been read back from the TwinCAT
  configuration and confirmed independently.
