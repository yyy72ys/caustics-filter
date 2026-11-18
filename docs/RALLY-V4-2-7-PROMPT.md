# RALLY-V4-2-7-PROMPT — 床の次の一手：設計の攻撃検証

You are an ADVERSARIAL REVIEWER, not a friendly assistant. Your job is to attack the design below and find where it is wrong, weak, or self-deceived. Do not praise it. Do not agree by default. For every claim, give a verdict: STANDS / WEAKENED / WRONG, with the concrete reason (a number, a mechanism, a counterexample, or a named paper). Then list: counterarguments the design missed, failure conditions under which it breaks, mechanisms it ignores, and your TOP 3 fixes ranked by impact.

## Context (what this project is)

We build a water-surface + caustics filter by reproducing Earth's environment in code. The bed (riverbed) must be a physically real water-worked gravel bed, judged by measured statistics, not by looks. The standard is REPRODUCING nature. A method defended only on "it looks natural" is rejected by definition.

## Measured evidence (all numbers are measured, not assumed)

- **v4-2.6 STAGE 1 STOP**: 0 of 22 mechanism combinations held local exponent >= 0.70 across 0.4–3.6 cm with sigma_z(matrix) in [0.32, 0.46]*D50 (D50 = 0.0362 m). Best worst-direction local-min 0.36. Support/deposition packing raises H by at most +0.084; imbrication 25–35 deg is nearly inert (delta H_x = +0.019); GSD width does NOT raise H (sigma_g 1.5 -> 2.7: H_x 0.315 -> 0.270); raising sigma_z within the band does not move H.
- **Ruler ceiling (v4-2.6 §4.3, extra check)**: even an ideal H_gen=0.79 self-affine field returns band local-min 0.65 in a 0.6 m box (finite-box bias ~0.1 at the band's upper end). Passing local >= 0.70 on ALL seeds needs H_gen >= 0.92. The H_fit band 0.70–0.85 IS reachable by ideal fields (0.719 @ 0.6 m / 0.750 @ 2.4 m).
- **Decomposition (v4-2.5, v4-2.6)**: in-band variance is ~100% grain relief (its own H ~0.336–0.34; grain-only H = 0.287); the fBm substrate (H_gen=0.79) contributes only ~8% in-band. The composite H IS the grain layer's H.
- **Rally settlement (v4-2.4, unanimous 3/3)**: Nikora H=0.79 is a SCALING-REGION fit value at lags ≲ D50; demanding it at 6.6–20 cm was a misplaced gate. Paving (flatten intra-patch scatter ≤ 2.5 mm) was weakened as unphysical by all three reviewers.
- **Hole-filling sensitivity**: the choice of hole-fill moves H_x by up to 0.134 (form-fill 0.473) but never reaches 0.70.
- **Meso band real-bed anchors (Bertin/Groom/Friedrich 2017; Powell et al. 2016)**: sigma_B/sigma_total median 0.64; sigma_B/D50 median 0.26; variogram ranges ~200 mm streamwise / ~100 mm cross-stream; apparent meso exponent 0.13–0.19 (span 0.09–0.26); anisotropy S_x/S_y ~0.83.
- **Real-bed organisation literature**: textural patches (Nelson 2010/2015, Buffington & Montgomery 1999); clusters cover 5–34% (Heays et al. 2014: only ~1/5 of potential anchors form clusters; Little & Venditti 2014: up to 8x flow resistance); hiding/exposure is a TRANSPORT concept (Wilcock & Crowe 2003, Parker 1990).

## The design under attack (DESIGN-V4-2-7-BED-NEXT.md)

**(A) Gate-form redesign — claims A1–A8:**
- A1: The old STAGE 1 gate ("local exponent >= 0.70 everywhere") is ABOLISHED — it fails ideal H_gen=0.79 fields and passes near-smooth H_gen>=0.92 fields, i.e. it gets both ends wrong.
- A2: New gate = profile-vs-profile comparison: the generated bed's directional H_loc(lag) profile against the reference profile of an ideal H_gen=0.79 fBm field measured with the SAME ruler and SAME box (per-box references; never cross-box).
- A3: Reference profile construction is a fixed procedure (spectral synthesis, PSD ~ k^-(2H+2), periodic boundary, R=8 for 0.6 m patch / R=4 for 3.0 m strip, measured verbatim through the same estimator functions).
- A4: Comparison metric = per-lag, per-direction z-scores; pass rule = three conditions ANDed (C1 mean |z| <= 1.0, C2 max |z| <= 3.0, C3 profile-drop match within 0.10) — ALL values explicitly provisional, to be re-decided by calibration, NOT derivations.
- A5: Estimator-bias calibration suite FIRST (synthetic fields H_gen in {0.70, 0.79, 0.85, 0.92}; acceptance: 0.79 fields pass >=7/8, 0.70/0.92 fields pass <=1/8, existing grain beds 0/3; hole-fill rule fixed before measurement; fill-flipping beds = "unresolved").
- A6: H_fit band 0.70–0.85 demoted to a coarse screen (profile is primary; profile-PASS + H_fit-out = unresolved).
- A7–A8: PASS/FAIL redefined as two-sided (too smooth also fails); PASS means "matches the reference profile", NOT "is a real bed"; fixed procedure with pre-registered thresholds.

**(B) Patch-scale organisation — claims B1–B7:**
- B1: Inherits the v4-2.4 patch design but DROPS paving, stretches correlation to 40–50 cm, switches kernel to Matérn nu=1.5.
- B2–B3: Concrete parameters (correlation length x 45 cm sweep 40–50, y 22 cm sweep 20–25, alpha 0.64 sweep 0.40/0.55/0.64, armour modulation R_a = 2.0*m^gamma, hiding/exposure as a secondary wake-sheltering variant only) with a feasibility guard (alpha_max ~ 0.70 at sigma_g=2.1).
- B4: bed_gen3.py composition sketch (patch field P -> per-grain modulation -> support/deposition scaled by local D -> imbrication kept for the dip gate, not H -> cluster anchors biased to coarse patches -> sigma_form re-iteration).
- B5: HONEST expected movement — meso band should close (Phi(6.6 cm) 0.84 -> 0.56–0.60, H_meso 0.08 -> 0.13–0.19) but scaling-band H_fit is expected to reach only 0.35–0.50, NOT 0.79; includes a cheap pre-check (mixture-of-cutoffs calculation) before any build.
- B6: Six falsifiers (F-B0–F-B5, including a 3-D grain killer test with effective grain size > 3.0 mm).
- B7: The two expected effects live in different boxes — scaling in the 0.6 m patch, meso only in the 3.0 m strip.

**Real-bed inventory — claims B8–B12:** textural patches, clusters, meso bedforms, hiding/exposure correlation, and a scale-hierarchy table — with citations.

**(C) M1–M6 meso-gate validation — claims C1–C3:** calibration suite first, re-centre bands on calibration readings, measure existing v4-2.2/2.3 strip beds (predicted: M3-in-band-for-the-wrong-reason, M4/M5 likely fail, no retroactive pass), then apply to bed_gen3 strips; A4 and M1–M6 are complementary non-overlapping bands, never merged into one verdict.

## Questions to answer (in addition to per-claim verdicts)

- Q1: Does the reference-profile construction (pure fBm vs knee-broken) change the verdict? How should "ideal field" be defined?
- Q2: Are the A4 provisional thresholds calibrated or merely placed? Are the calibration-suite acceptance criteria strict enough?
- Q3: Is two-sided testing (failing too-smooth fields) correct, and consistent with "reproduce nature"?
- Q4: Is the 40–50 cm correlation length, scaled by D50, defensible? Powell is one measurement; Bertin's recurrence is a wavelength, not a range.
- Q5: Does the "mixture effect" of varying local D50 really raise the scaling-band slope? What condition must the area-weighted average of local curves satisfy? Is the cheap F-B0 check enough, or should the 3-D grain killer test (F-B4) come first?
- Q6: Is translating hiding/exposure (a transport concept) into a placement rule ("fine grains in coarse grains' wakes") legitimate, or is it "looks-like" tuning without a formation process?
- Q7: Should the killer test (directional 2-D structure functions vs a real measured armour-bed DEM) come before building bed_gen3? Which comes first — obtaining Bertin's DEM or the calibration suite?
- Q8: M1–M6 bands — is the M3 lower bound 0.08 too close to a grains-only bed (0.079)? Is M6's D50 scaling defensible? How sensitive is M5's separation (2.5*D90) to D90 uncertainty (75–95 mm)?

## Answer format (fixed)

For EACH claim A1–A8, B1–B12, C1–C3: verdict STANDS / WEAKENED / WRONG + one paragraph of concrete reasoning. Then: (1) counterarguments the design missed, (2) failure conditions, (3) mechanisms it ignores, (4) your TOP 3 fixes ranked by impact. Then answer Q1–Q8 briefly. No vague agreement. If a claim survives, say exactly what evidence saves it.
