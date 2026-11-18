# RALLY-V4-2-4 — v4-2.4 patch design (texture patch + paving)

Topic: v4-2.4 bed multi-scale organisation design (DESIGN-V4-2-4.md, BRIEF-V4-2-4.md).
New topic → new threads (one per AI, never fork).

- GPT: https://chatgpt.com/c/6af6375e-b860-83ec-a2fa-4fbf653bd9d8
- Grok: https://grok.com/c/139fd63b-4ebb-4df3-944c-8e7cc6d821b5?rid=f3614d8f-d00a-4a5e-8858-3bf32c19db3c (account むせすぱ; leg took 42 min — grok.com root was bot-blocked, recovered via grok.com/chat direct, email-login Turnstile failed repeatedly, Google sign-in succeeded, anonymous first send lost and re-sent signed-in)
- Gemini (AI Studio): https://aistudio.google.com/prompts/1SP8mvzHbg2MrthZmrVVzXBpFVQudTDD6 (model Gemini 4 Hydrogen, thinking High)

Design under attack: texture-patch field P(x,y) tying local grain size and height envelope (coarse patches higher), anisotropic Gaussian correlation 24 cm x 12 cm (2:1), "paving" (intra-patch grain height scatter a <= 2.5 mm), local grading narrowed by sigma_ln,local = sqrt(sigma_ln,global^2 - alpha^2) with alpha = 0.64, substrate sigma lowered to re-close sigma_z; calibration stage first (H >= 0.72 at sigma_z = 0.40*D50 or stop). Gate and band NOT moved.

Established evidence given to all three: 2.3 audit — composite H 0.219–0.232 vs band 0.70–0.85; substrate-only H 0.796–0.809; grain relief H 0.103 dominating S (95% @4 cm, 83% @10 cm, 66% @20 cm); pure fBm control 0.58 (0.6 m) → 0.74 (3.0 m).

## Verdicts

| Claim | GPT (07:03) | Grok | Gemini (07:03) |
|---|---|---|---|
| D1 patch layer is physical (one field P + paving) | WEAKENED | WEAKENED | WEAKENED |
| D2 H gate legitimate as written (vs mis-specified lag band) | WRONG | WRONG | WRONG |
| D3 variance (S-ratio) budget correct / corner forced | WEAKENED | WEAKENED | WRONG |
| D4 other 7 gates protected | WRONG | WEAKENED | WRONG |

### GPT detail (thread above)
- D1 WEAKENED: coarse clusters rising / fines sheltering shows a size–height link, but does not support generating both from one stationary anisotropic Gaussian field (2:1), nor flattening intra-patch relief to 2.5 mm. Counter-example: measured clusters with cluster-scale humps AND large grain-scale relief.
- D2 WRONG (as a claim about natural-bed H): Nikora 0.79 is a small-lag scaling-region value; 6.6–20 cm (≈1.8–5.6 D50) is transition/saturation, so recovering Nikora's exponent there is unjustified. Gate should move to the scaling region; if 6.6–20 cm stays, rename it as an independent requirement (e.g. H_meso) with its own reference population. The patch itself may remain as a meso-organisation device.
- D3 WEAKENED: composite variance is not a sum of three independent variances — covariance terms (esp. P × grain relief coupling) exist, so the corner (a ≤ 2.5 mm, b ≈ 14 mm, sigma_form ≈ 2–4 mm) is not forced. Anisotropy effects on the estimator, correlation length not fixing the contribution, double-counting risk of cluster + patch envelope. Calibration sweeps should measure the actual composite field, not marginal-variance arithmetic.
- D4 WRONG (asserted, not shown): even preserving global GSD and sigma_z, the effect of local narrowing + flattening on packing, porosity and imbrication dip is unmeasured; the 7 gates are "expected to survive", not protected. Lowering substrate sigma to recover sigma_z may shift skewness.
- Most likely failure mode: the patch field numerically hits target H but fabricates a meso-scale covariance pattern that does not represent a water-worked gravel bed — main gate passes while physical interpretation and the other 7 gates are damaged.
- Fastest killer measurement: an independent high-resolution DEM of a water-worked gravel surface with comparable D50 and armour ratio; measure directional/2-D structure functions over the same 4–20 cm band and compare directional covariance / spectral distribution with the generated surface. H matching alone with mismatched directional structure kills the central physical claim.

### Gemini detail (thread above)
- D1 WEAKENED: qualitative support from Bertin et al. (2017) for the common coupling to P, but the 24 cm × 12 cm Gaussian field and a ≤ 2.5 mm paving are unphysical for a water-worked gravel bed. At D50 = 36 mm, 24×12 cm is only 6.7 D50 × 3.3 D50; a ≤ 2.5 mm is < 0.07 D50, inconsistent with natural gravel-crest curvature (≈0.2–0.5 D50, 7–18 mm).
- D2 WRONG: the H gate (0.70–0.85, lags 6.6–20 cm) is mis-specified; the patch layer is an artificial fix to an invalid target. Nikora H ≈ 0.79 is the scaling-region value (lags ≲ D50 ≈ 3.6 cm); 6.6–20 cm (1.8–5.6 D50) is transition/saturation. Move the lag band below the break (≈0.4–3.6 cm), or gate the 10–50 cm band by real meso-scale correlation length / variance indices instead of H.
- D3 WRONG: the budget is not physically forced to a ≤ 2.5 mm and misses coupling and discretisation terms. Escape routes: Matérn/fBm-type patch kernel, expanding Lx/Ly to ~40–50 cm, streamwise-only S(d) evaluation. Missed: positive Cov(z_patch, z_grain); heterogeneous variance proportional to local grain size (mean variance inflated by exp(2α²) ≈ 2.27); in 3-D grains, even with coplanar centres, sphere/ellipsoid caps and interstitial saddles push a past 3 mm.
- D4 WRONG: the other 7 gates are actively damaged by alpha = 0.64 and a ≤ 2.5 mm. A 15–30°-tilted 36 mm grain spans 9–18 mm vertically, so a ≤ 2.5 mm forces grains nearly flat, breaking imbrication dip. Local narrowing: at D84/D50 = 1.85 (sigma_ln,global = 0.615) sigma_ln,local becomes imaginary; at D84/D50 = 2.1, sigma_ln,local = 0.37, nearly monodisperse — pore filling by small grains is lost and porosity / armour ratio / negative skewness destabilise.
- Most likely failure mode: calibration-to-realization collapse. Heightfield-only calibration passes H ≥ 0.72, but DEM-ising discretised imbricated 3-D grains pushes a_eff past the 3.0 mm cliff and H falls back below 0.65 (or the bed becomes a flat cobbled floor and dip/porosity/skewness fail).
- Fastest killer measurement: on a flat substrate (b = 0, sigma_form = 0), build ONE 0.5 m × 0.5 m 3-D grain fill satisfying the imbrication-dip and GSD gates, and measure the DEM structure function S_grain(6.6 cm) (or in-patch RMS a_eff). If a_eff > 3.0 mm (S_grain(6.6 cm) > 18 mm²), the variance budget and the whole construction pipeline fail before any full-bed build.

### Grok detail (thread above)
- D1 WEAKENED: Bertin's observation (coarse clusters rise, fines in lows) is real and the single anisotropic Gaussian field P tying size and height has qualitative physical basis. But paving (flattening intra-patch relief to ≤ 2.5 mm) differs from measured water-worked beds: grain-scale relief inside coarse patches also contributes to the structure function, so the implementation is not equivalent; a flow-deposited version would diverge, especially in the 10–30 cm band.
- D2 WRONG: the gate itself is mis-specified. Nikora H ≈ 0.79 is the scaling-region value (lags ≲ D50 ≈ 3.6 cm); the gate's 6.6–20 cm is transition/saturation, where measured structure functions flatten. The lag band should be lowered to the true scaling region; the patch layer is an artifact of a misplaced gate, not the answer to a correctly placed one.
- D3 WEAKENED: the corner (a ≤ 2.5 mm, b ≈ 14 mm, sigma_form ≈ 2–4 mm) is a direct consequence of the variance budget, and the fragility at a = 3 mm is conceded. But escape routes are not fully closed: the 24 cm × 12 cm anisotropy changes directional structure functions and the effective lag sampling of the 2-D estimator, and covariance between the grain-size field and local packing density can move variance between the grain term and the shape term without changing the marginal amplitudes of a and b.
- D4 WEAKENED: the protection of the 7 gates by local narrowing + lowered substrate sigma is algebraically consistent. But paving lowers porosity and narrows the long-axis dip distribution relative to a natural armour, so even holding mean dip, hydraulic roughness and the entrainment threshold differ — the protection is superficial.
- Most likely failure mode: a surface numerically passes the H gate at 6.6–20 cm while its grain-scale packing, local porosity and directional anisotropy remain distinguishable from water-worked gravel.
- Fastest killer measurement: direct comparison of a true water-worked armour (lab/field DEM, same D50 and armour ratio) with the generated surface, in the true scaling region (lags ≲ D50) and cross-stream, via 2nd-order structure functions (or 2-D power spectra).

## Synthesis (complete, 2026-11-13 07:43)
- Unanimous: D2 WRONG (3/3) — the H gate's lag band is mis-specified. Nikora's 0.79 belongs to lags ≲ D50; demanding it at 6.6–20 cm (transition/saturation) is unjustified. The patch is an artifact of a misplaced gate.
- D1 WEAKENED (3/3): the single-field coupling has qualitative support, but paving (a ≤ 2.5 mm) is not equivalent to measured beds.
- D3 WEAKENED ×2 + WRONG ×1: the budget misses covariance (P × grain relief, size × packing) and discretisation; the corner is not forced; escape routes exist (Matérn/fBm patch kernel, Lx/Ly ~40–50 cm, streamwise-only S(d)).
- D4 WRONG ×2 + WEAKENED ×1: the 7 gates are asserted, not shown protected; paving forces grains nearly flat (tilted 36 mm grains span 9–18 mm) and local narrowing can go imaginary (D84/D50 = 1.85) or near-monomodal.
- Convergent killer test (all three, in different forms): compare the generated surface against a real water-worked armour DEM at the same D50/armour ratio, directional 2-D structure functions in the true scaling region — H alone matching proves nothing.
- FINAL: do NOT implement v4-2.4 as written. The decision is the user's: (a) move the H gate's lag band to the scaling region (≈0.4–3.6 cm) and re-audit existing beds — the patch becomes unnecessary; (b) keep 6.6–20 cm as a renamed independent meso requirement (H_meso) with its own reference population — needs research on what real beds do there; (c) rework the patch design against the D1/D3/D4 attacks (bigger scales, Matérn/fBm kernel, streamwise-only S(d), real 3-D grain killer test first) and keep the current gate. Presented to the user as a choice 2026-11-13 07:43.
