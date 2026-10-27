# Data inventory

No source files are copied into this folder; the analysis uses the canonical
files already held by the stage-identification project.

| Source | Status | Role |
|---|---|---|
| `../../iso230-2/analysis/iso230_2_v3_results.json` | available | 24-target directional means, sample standard deviations, counts and ISO metrics |
| `../../iso230-2/data/terminal_20261021_183115.log` | available | commanded targets, microstep values, timing, run completion and driver settings |
| `../../ISO230-2_v3.csv` | missing | raw 1 kHz AquadB position capture required for drift correction, per-stop settling and per-run reconstruction |
| environment log | missing | air/stage temperature, pressure and interferometer compensation state |

The missing raw capture was also searched for in the workspace, Downloads,
Desktop, mounted volumes and the local Spotlight index. Restore it as
`../../ISO230-2_v3.csv` before extending this project with raw-window analyses.

