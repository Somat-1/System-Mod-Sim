# Sequence provenance

This session has **no sequence files of its own**, and that is the finding, not an
omission.

The 30 Oct run was produced by the unmodified `e2_global_map` firmware
(`config=96579c689d4d`), which derives its target list from `seed_start: 5200` in
`e2_config.json`. That value was not changed between sessions, so the firmware
settled on the same accepted seed, `5315`, and emitted the same 400 targets. The
400 `# TARGET` lines of the two terminal logs diff clean.

The plan is therefore read directly from

    ../../e2_global_map/sequence/e2_measurement_plan.csv

rather than copied here. Reading the one file both sessions were driven by makes
the identity of the two target sets structural: it cannot drift out of date, and
no later reader can mistake a copy for an independently generated plan.

**Consequence for the procedure.** `E3_processing_instructions.md` specifies E3 as
the same design with a *different* seed, so that the map is tested at *new*
positions. This session is a different day at the *same* positions. See
`../README.md` for which questions that still answers and which it does not.
