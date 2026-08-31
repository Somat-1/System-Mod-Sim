# Microstepping Test Data v2

This revision contains the physical-stage motion-sequence specification and its
command-only verification figures. It is intentionally separate from the
archived v1 analysis workflow.

- `BACKLOG.md`: authoritative sequence setup, rationale, open hardware inputs,
  and the later measured-response plotting plan.
- `scripts/generate_command_montages.py`: deterministic command-preview
  generator. It does not issue hardware commands.
- `scripts/run_identification_esp32_tmc2209/`: Arduino-ESP32/TMCStepper/RMT
  runner for Blocks A, B, and E.
- `scripts/run_identification_dedicated_controller.py`: host-side serial
  runner for Blocks C and D.
- `data/motion_sequence_config.json`: machine-readable execution and preview
  parameters.
- `rendered_assets/trajectory_visualization_plots/`: command montages grouped
  by test purpose.
- `rendered_assets/command_montage_summary.json`: cell-level validation metadata
  and paths to every trajectory figure.
- `docs/`: reserved for derivations, conventions, and run notes.
- `docs/HARDWARE_RUNNERS.md`: wiring, safety gates, backend responsibilities,
  dry-run checks, and live invocation.

Regenerate from this directory with:

```powershell
python .\scripts\generate_command_montages.py
```

All executions start from the same fixed stage position. Hardware execution
remains gated on pilot confirmation of the ladder dwell/repeat count and the
live safety arguments documented in `docs/HARDWARE_RUNNERS.md`.
