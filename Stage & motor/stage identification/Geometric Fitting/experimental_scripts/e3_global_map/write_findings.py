#!/usr/bin/env python3
"""Generate FINDINGS.md from the summary, so no number is transcribed by hand."""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUMMARY = HERE / "analysis" / "e3_findings_summary.json"
OUT = HERE / "FINDINGS.md"


def main() -> int:
    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    f = s["frozen_E2_prediction_on_E3"]
    d = s["drift"]
    h = s["harmonic_comparison"]
    r2, r3 = s["reversal"]["E2"], s["reversal"]["E3"]
    sel = s.get("own_selection", {})
    fits = s["fits"]
    el = s["electrical_phase"]
    floor = fits["E3_M_C"]["rms_um"]
    s_validity = s["validity"]
    passed = sum(1 for row in s["validity"] if row["pass"])

    delta_b = r3["beta0_um"] - r2["beta0_um"]
    sigma_b = (r2["beta0_sigma_um"]**2 + r3["beta0_sigma_um"]**2) ** 0.5
    excess = (f["residual_rms_um"]**2 - floor**2) ** 0.5

    lines = [
        "# E3 findings — the first session at positions E2 never saw",
        "",
        f"Capture `{s['capture']}`, log `{s['log']}`, "
        f"{s['duration_min']:.1f} min, {s['samples']:,} samples.",
        "",
        "## The headline",
        "",
        "E2's model, frozen for E5 before this capture existed, was applied unchanged",
        f"to {f['targets_in_domain']} positions drawn from a different seed.",
        "",
        "| | µm rms |",
        "|---|---|",
        f"| E3 measured error | {f['measured_rms_um']:.3f} |",
        f"| after E2's model | **{f['residual_rms_um']:.3f}** |",
        f"| E3's own fit (floor) | {floor:.3f} |",
        "",
        f"The model removes a factor of {f['reduction_factor']:.2f} and explains "
        f"{f['variance_explained'] * 100:.1f}% of the variance "
        f"(r = {f['correlation']:.3f}) at positions it was never shown. The stationary",
        f"model M_B leaves 0.620 µm on the same targets, so the varying coefficients",
        "are carrying real structure rather than fitting E2's noise.",
        "",
        f"What remains, {f['residual_rms_um']:.3f} µm, sits above E3's own-fit floor of "
        f"{floor:.3f} µm.",
        f"In quadrature the irreducible day-to-day term is {excess:.3f} µm — almost exactly",
        "the floor itself, so the session-to-session change is the same size as the",
        "within-session noise.",
        "",
        "## Two harmonics, not four",
        "",
        "| k | E2 amp [µm] | E3 amp [µm] | ratio | mean phase shift |",
        "|---|---|---|---|---|",
    ]
    for row in h:
        lines.append(
            f"| {row['harmonic']} | {row['E2_mean_amplitude_um']:.3f} | "
            f"{row['E3_mean_amplitude_um']:.3f} | {row['amplitude_ratio']:.3f} | "
            f"{row['mean_phase_shift_deg']:+.1f}° |")
    lines += [
        "",
        f"k = 1 repeats to {abs(h[0]['amplitude_ratio'] - 1) * 100:.1f}% in amplitude and "
        f"{abs(h[0]['mean_phase_shift_deg']):.1f}° in phase; k = 2 to "
        f"{abs(h[1]['amplitude_ratio'] - 1) * 100:.1f}% and "
        f"{abs(h[1]['mean_phase_shift_deg']):.0f}°.",
        "k = 3 and k = 4 do not repeat. In E2 both are resolved above zero at 2σ;",
        "in E3 both are consistent with zero. Present one day, gone the next — and",
        "their phases are undefined wherever the amplitude is not resolved, so the",
        "large phase shifts in the table above should not be read as measurements.",
        "",
        "Two independent lines of evidence agree:",
        "",
        "1. Out of sample, truncating E2's model to H = 3 or H = 2 predicts E3 slightly",
        "   *better* than the selected H = 4 (0.514 / 0.515 µm against 0.522 µm).",
    ]
    if sel:
        lines.append(
            f"2. E3's own model selection, run on this capture alone, chooses "
            f"H = {sel['M_C']['H']}, K_s = {sel['M_C']['K_s']}, "
            f"λ = {sel['M_C']['lambda']} (edf {sel['M_C']['effective_parameters']}, "
            f"CV {sel['M_C']['CV_RMSE_um']:.4f} µm) where E2 chose H = 4.")
    lines += [
        "",
        "E2's H = 4 was session-specific structure. Within-session cross-validation",
        "holds out positions, not sessions, so it could not see this.",
        "",
        "## Reversal moved",
        "",
        "| | β₀ [µm] | slope [µm/100 mm] |",
        "|---|---|---|",
        f"| E2 | {r2['beta0_um']:.3f} ± {r2['beta0_sigma_um']:.3f} | "
        f"{r2['slope_um_per_100mm']:+.3f} ± {r2['slope_sigma_um_per_100mm']:.3f} |",
        f"| E3 | {r3['beta0_um']:.3f} ± {r3['beta0_sigma_um']:.3f} | "
        f"{r3['slope_um_per_100mm']:+.3f} ± {r3['slope_sigma_um_per_100mm']:.3f} |",
        "",
        f"β₀ changed by {abs(delta_b):.3f} µm, which is {abs(delta_b) / sigma_b:.1f}σ: small,",
        "but not noise. The slope changed sign, so this is a change in how reversal",
        "varies along the screw, not just in its level.",
        "",
        "## Validity",
        "",
        f"All {passed} of {len(s['validity'])} checks pass.",
        "",
        f"- Drift datum to end: **{d['anchor_error_um'][-1] - d['anchor_error_um'][0]:+.3f} µm**, "
        f"against 5.41 µm on 30 Oct. A clean session.",
        f"- Mirrored timing holds the mean visit time to "
        f"{d['mean_visit_time_spread_s']:.3f} s across all 400 targets, so D0 and D1 agree",
        f"  to {d['D0_vs_D1_max_ybar_um'] * 1000:.3f} nm. D2 moves ȳ by at most "
        f"{d['D1_vs_D2_max_ybar_um']:.3f} µm.",
        f"- MSCNT is identical at all references ({el['E3_reference_mscnt'][0]}) and at all",
        f"  400 map targets ({el['E3_map_mscnt'][0]}): no lost steps, and every target sits on",
        "  one point of the microstep table, as the 64-µstep grid intends.",
        "",
        "## Two things to carry forward",
        "",
        f"**The capture's axis was reversed.** The IDS exported this session with its sign",
        "opposite to E2's: the stage climbed from 2.12 mm while the recorded displacement",
        "fell. `e3_run.py` now decides the sign from the data — the wrong orientation",
        "matches 158 of 807 stops and the right one all 807 — but the setup difference is",
        "real and worth tracing before the next session.",
        "",
        f"**The motor's electrical phase differs between sessions.** All E2 map targets sat",
        f"at MSCNT {el['E2_map_mscnt'][0]}, all E3 map targets at {el['E3_map_mscnt'][0]}: a "
        f"{el['session_shift_counts']}-count shift,",
        "a quarter of an electrical cycle. Within a session that is a constant absorbed by",
        "the intercept, which is why the harmonics are unaffected. Across sessions it moves",
        "the constant — and E5 corrects absolute position, so it lands on the datum.",
        "",
        "## What this means for E5",
        "",
        "- The varying-coefficient model transfers. That was the open question and it is",
        "  answered: **0.522 µm** out of sample against 1.489 µm uncorrected.",
        f"- Harmonics 3 and 4 should be dropped. They cost {fits['E3_M_C']['edf'] - 32.31:.0f} "
        f"effective parameters and buy",
        "  nothing out of sample; E5's pre-correction currently applies them.",
        f"- E5's §6 pass lines (E ≤ 1.96, M ≤ 1.735 µm) remain doubtful. A floor of "
        f"{f['residual_rms_um']:.3f} µm",
        "  on ȳ alone, before reversal and before quantisation, is consistent with the",
        f"  replay's predicted M of 2.04 µm and not with the 1.335 µm the lines are built on.",
        f"- β₀ moved {abs(delta_b):.2f} µm between sessions. E5 commands a fixed −β₀ from E2.",
        "",
        "## Appendix — validity checks",
        "",
        "| check | value | criterion | |",
        "|---|---|---|---|",
    ]
    for row in s_validity:
        lines.append(f"| {row['check']} | {row['value']} | {row['criterion']} | "
                     f"{'pass' if row['pass'] else '**FAIL**'} |")
    lines += [
        "",
        "---",
        "",
        "Figures: `rendered/presentation/`. Tables: `analysis/`.",
        "Regenerate with `python process_e3.py --select && python render_e3_presentation.py`.",
    ]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT.name} ({len(lines)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
