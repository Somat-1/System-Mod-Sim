# Microstepping Test Data v2

This revision contains the physical-stage motion-sequence specification and its
command-only verification figures. It is intentionally separate from the
archived v1 analysis workflow.

- `BACKLOG.md`: authoritative sequence setup, rationale, open hardware inputs,
  and the later measured-response plotting plan.
- `scripts/generate_command_montages.py`: deterministic command-preview
  generator. It does not issue hardware commands.
- `data/motion_sequence_config.json`: machine-readable execution and preview
  parameters.
- `rendered_assets/trajectory_visualization_plots/`: command montages grouped
  by test purpose.
- `rendered_assets/command_montage_summary.json`: cell-level validation metadata
  and paths to every trajectory figure.
- `docs/`: reserved for derivations, conventions, and run notes.

Regenerate from this directory with:

```powershell
python .\scripts\generate_command_montages.py
```

All 12 executions start from the same fixed stage position. Hardware execution
remains gated on the holding-current values and pilot confirmation of the
ladder dwell and execution repeat count. See `BACKLOG.md` before converting
this preview into a controller-specific motion program.
