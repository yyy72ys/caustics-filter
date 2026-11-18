# BRIEF — 第4段 訂正再計算（Snell 精確 f＋体制地図＋50mm、√(2CD) 式）

- Job: `caustics-v4-4-recalc`
- Source of this brief: `RALLY-V4-4-1-CONFIRM-R2.md` §6, mechanical items 1/3/4 (the integration ruling of the stage-4 confirmation rally, round 2, 2026-11-16). This job is PURE COMPUTATION. It makes no design decisions and changes no design text. Its job is to produce corrected numbers that the R3 draft will later adopt.
- Inputs to read on the PC: `D:\yohey2026\dsh\muse\caustics-filter\DESIGN-V4-4-1-R2-DRAFT.md` (N1–N6, especially N1/N2/N5), `D:\yohey2026\dsh\muse\caustics-filter\offline\out\v4-4-map\` (REPORT + numbers + tables), `D:\yohey2026\dsh\muse\caustics-filter\offline\out\v4-4-map50\` (REPORT + numbers + tables). If any of these is missing on the PC, say so in the report and work from the input numbers below — do not invent replacements.

## Background (why this job exists)

Two errors were found in the round-2 rally (Gemini seat, independently verified by the integrator):

1. **The focal lengths in the R2 draft were computed with geometric distance and dropped Snell refraction.** The draft computed f from s = h_a + d_w: main 38.11 mm (s = 1.6 m, W = 1.4755 m), sub 37.22 mm, contrast 23.65 mm. The Snell-exact values are main ≈ 34.79 mm, sub ≈ 28.24 mm, contrast ≈ 20.86 mm. With f = 38.11 the real bed window is 1.350 m (8.5% narrower than the claimed 1.4755 m); for the sub frame the error is −22.2%. Every regime-map and 50 mm number built on those f values inherits the error.
2. **The b_eff minima in the first regime-map job were overestimated by a factor of √2.** The formula used was 2·sqrt(C·D); the correct minimum is sqrt(2·C·D). The corrected minima for the OLD f are main 2.06 / sub 3.12 / contrast 1.26 px (the A/B conclusions did not change with the old f). The 50 mm job used the correct formula.

## Input numbers (given; the job computes from them, it does not re-derive them)

- Water refractive index n = 1.333 (air n = 1). Calibration relation (draft N2): apparent distance s_app = h_a + d_w/n.
- Standard camera (all three frames, h_a = 1.0 m): main — window W = 1.4755 m, depth d_w = 0.6 m; sub — W = 4.8 m, d_w = 4.0 m; contrast — W = 2.4 m, d_w = 0.6 m. Sensor width W_s = 36 mm (draft N1). The draft's window formula: W = W_s·(s/f − 1) with s = h_a + d_w (geometric — this is the WRONG relation to be replaced by the Snell-exact one).
- 50 mm same-window variants (windows unchanged): main h_a = 1.4993 m, sub h_a = 2.7167 m, contrast h_a = 2.7833 m.
- 50 mm same-height variants (h_a = 1.0 m): depth 0.6 m frames window 1.1160 m (main and contrast share it), sub window 3.5640 m.
- The N sweep for the regime maps is the same as in the v4-4-map / v4-4-map50 jobs (read it from their tables; do not change the sweep).

## STAGE 1 — Snell-exact focal lengths (the job's own derivation)

Derive, from first principles, the exact thin-lens focal length f that makes the bed window W fill the sensor width W_s when the lens sits at height h_a above a flat water surface of depth d_w. The ray path refracts at the surface (Snell, n = 1.333). Do NOT hard-code the expected answers. The expected values below are a VERIFICATION GATE only: your computed f must land within 0.05 mm of main 34.79 / sub 28.24 / contrast 20.86 mm. If it does not, stop and report the discrepancy — do not adjust the derivation to fit. State the derivation (which ray, which apparent-distance relation, paraxial or exact) in the report so it can be checked. Also recompute the bed window that the OLD f values actually give (main should come out ≈ 1.350 m) as a cross-check of the same derivation.

## STAGE 2 — Regime maps with corrected f and corrected minima formula

Recompute the full regime maps (same N sweep, same threshold: the other face's effective blur b_eff ≤ 2 px is regime A) for:
(a) the standard camera, all three frames, with the Snell-exact f;
(b) the 50 mm variants (same-window main/sub/contrast, same-height main/contrast at 0.6 m and sub), with the Snell-exact f for each variant (derive each variant's f the same way as Stage 1).
Use b_eff minima via sqrt(2·C·D) throughout. For each frame/variant report: the A interval in N (if any), the minimum b_eff and the N at which it occurs, and which N cells (if any) are diffraction-dominated (regime D). State explicitly whether the main frame now has an A interval (the integrator's rough estimate said its minimum falls below 2 px, ≈ 1.9 px — verify, do not assume).

## STAGE 3 — Kappa critical table and Coddington split with both denominators

(a) κ_crit table: κ_crit = n/((n−1)·d_w) for d_w = 0.6 m and 4.0 m (expected 6.67 /m and 1.00 /m). List, for each frame, the declared κ_max from the draft (main |κ1| p99 = 53 /m; sub unmeasured, provisional) against κ_crit, and state in numbers how far beyond the single-image-branch limit each frame sits. No new κ measurements — use the draft's numbers.
(b) Coddington astigmatism split: recompute the split ratio at the draft's angles with BOTH denominators: (d_s − d_t)/d_w (the draft's 10% figure at 31.82°) and (d_s − d_t)/d_s (expected ≈ 14.4% at the same angle). For each denominator, give the angle θ at which the split reaches 10% (θ_max under that definition). Report both θ_max values and state which one the draft should adopt.

## Deliverables (new files only, under out/v4-4-recalc/)

- `REPORT-V4-4-RECALC.md` — derivation (Stage 1, checkable), all tables, the A/B verdict per frame/variant, the κ_crit table, both θ_max values, and a caveats section (what was NOT recomputed: the sensor-geometry definition fix is a draft-text matter for the R3 draft stage; the glitter D_glit applicability and the pupil-sub-ray question are design decisions, not this job).
- `numbers-v4-4-recalc.json` — every number in the report, machine-readable, with the Stage-1 gate values (computed f, expected f, difference) as explicit fields.
- `tables-v4-4-recalc.md` — the regime-map tables in the same format as the v4-4-map tables, so they can be diffed directly.
- Checkpoint: write `out/v4-4-recalc/checkpoint.json` with progress units as the stages complete (the room's watcher reads it).

## Rules

- Do not modify any existing *.py or any file outside out/v4-4-recalc/. New files only.
- Do not hard-code the expected f values as outputs. They are a gate, not an input.
- All numbers in the report must also appear in the JSON. If a number cannot be computed, say so and mark it missing — never fill it with a guess.
- The report's verdict section must state, for each of the standard frames and each 50 mm variant, A or B and the interval. No hedging: if the main frame's minimum is below 2 px, say it has an A interval and give it.
