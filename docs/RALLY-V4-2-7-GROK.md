# RALLY-V4-2-7-GROK (verbatim, Grok, device browser)

Thread: https://grok.com/c/ea1ac1d2-9d66-4c0b-b6d0-3e3f63b71bdf (bed topic, continued; via device browser — cloud browser was blocked)
Prompt: RALLY-V4-2-7-PROMPT.md (verbatim in the task)

---

A1
Verdict: STANDS
The old local-exponent gate is falsified by its own numbers: ideal H_gen=0.79 fields return band local-min ≈0.65 under the same ruler and box size that the gate demands ≥0.70, while fields with H_gen≥0.92 can pass. That double failure is decisive; the gate cannot be retained.

A2
Verdict: WEAKENED
Profile-vs-profile comparison is a necessary correction to a broken scalar gate, but it still treats an ideal self-affine field as the sole reference. Real armour beds are not pure fBm; they carry knee breaks, textural-patch anisotropy and meso-scale saturation that a pure H=0.79 spectrum does not contain. Matching only the ideal profile therefore risks accepting beds that are still optically wrong.

A3
Verdict: STANDS
A fixed, pre-registered spectral synthesis procedure (PSD ~k^-(2H+2), periodic boundaries, stated R values, identical estimator code) is the only way to keep the reference free of post-hoc tuning. The claim is saved by procedural transparency, not by any claim that the resulting field is "natural."

A4
Verdict: WRONG
The three conditions and the numerical thresholds (mean |z|≤1.0, max |z|≤3.0, drop match ≤0.10) are explicitly labelled provisional yet are written into the gate as if they already possess statistical meaning. No calibration suite has yet measured the sampling distribution of those z-scores on pure fBm realisations of the same size; the numbers are therefore arbitrary and cannot be used as a pass rule.

A5
Verdict: WEAKENED
Running a calibration suite first is correct in principle, but the acceptance criteria (≥7/8 for H=0.79, ≤1/8 for the extremes) are themselves uncalibrated guesses about the estimator's Type-I/II error rates. Until those rates are measured, the suite cannot certify the gate.

A6
Verdict: STANDS
Demoting the scalar H_fit band to a coarse screen is required once the profile is declared primary; the finite-box bias already shows that H_fit alone is an unreliable discriminator. The claim is saved by consistency with the ruler-ceiling evidence.

A7
Verdict: STANDS
Two-sided failure (too-smooth also fails) follows directly from the definition of the reference profile. A bed that is smoother than the ideal field at the lags that matter for caustics is not a match, regardless of any "more natural" intuition.

A8
Verdict: STANDS
PASS is correctly redefined as "matches the reference profile under the fixed ruler," not "is a real bed." That semantic restriction prevents the classic circular claim that a bed is natural because it passed a gate that was itself tuned to look natural.

B1
Verdict: WEAKENED
Dropping paving is required by the earlier unanimous review, but stretching the correlation length to 40–50 cm and switching to Matérn ν=1.5 are free parameters whose only justification is "to close the meso band." No formation process or measured covariance function is cited that would produce those lengths from grain-scale mechanics.

B2
Verdict: WRONG
The concrete numbers (45 cm / 22 cm, α=0.64, R_a=2.0 m^γ) are taken from single-study medians or from unconstrained sweeps. Powell's variogram range is one measurement under one flow regime; Bertin's recurrence wavelength is not a correlation length. Treating them as design targets converts measured statistics into free knobs.

B3
Verdict: WEAKENED
The feasibility guard (α_max≈0.70 at σ_g=2.1) is a useful numerical check, but it does not address whether the resulting local D50 modulation still preserves the measured armour-ratio distribution or the observed cluster-size statistics. The guard is necessary but far from sufficient.

B4
Verdict: WEAKENED
The composition sketch (patch field → per-grain modulation → support scaled by local D → imbrication kept only for the dip gate) is a coherent engineering pipeline, yet it still embeds the same grain-layer dominance that produced H≈0.3 in every prior experiment. Nothing in the sketch alters the fact that in-band variance remains ~100 % grain relief.

B5
Verdict: STANDS
The honest forecast that scaling-band H_fit will reach only 0.35–0.50 while the meso band may close is the single most valuable statement in the design. It correctly refuses the earlier demand that grain-based methods deliver Nikora's H=0.79 at the wrong lags.

B6
Verdict: WEAKENED
Six falsifiers are better than none, but F-B4 (3-D grain killer with effective grain size >3 mm) is still downstream of the build. A falsifier that can be applied only after bed_gen3 exists cannot protect the design decisions that precede the build.

B7
Verdict: STANDS
Separating the two effects into different box sizes (scaling in 0.6 m, meso in 3.0 m) is required by the measured scale hierarchy and by the finite-box bias already quantified. The claim is saved by that scale separation.

B8
Verdict: STANDS
Textural patches are a documented real-bed organisation (Nelson, Buffington & Montgomery). Their inclusion in the inventory is correct; the claim does not yet assert that the Matérn construction reproduces them.

B9
Verdict: STANDS
Cluster coverage 5–34 % and the observation that only a minority of potential anchors form clusters are measured facts. Listing them is not a design claim and therefore survives.

B10
Verdict: STANDS
Meso-bedform statistics (σ_B/σ_total median 0.64, ranges ~200 mm / 100 mm, apparent exponents 0.13–0.19) are taken from the cited field studies. The inventory is accurate.

B11
Verdict: WRONG
Hiding/exposure is a transport and entrainment concept (Wilcock & Crowe, Parker). Converting it into a static placement rule ("fine grains in the wakes of coarse grains") has no measured formation process behind it; it is a visual heuristic.

B12
Verdict: WEAKENED
A scale-hierarchy table is useful, but the table as written still treats the grain layer as the carrier of the scaling-band H that the earlier decomposition showed it cannot carry. The hierarchy is therefore incomplete until the substrate contribution is restored or the grain-layer H problem is solved.

C1
Verdict: STANDS
Calibration suite first, then re-centring of the M-bands, is the only statistically defensible order. The claim is saved by that procedural insistence.

C2
Verdict: WEAKENED
Measuring the existing v4-2.2/2.3 strips and predicting that M3 will pass "for the wrong reason" is useful diagnostic work, but the prediction itself is not yet a measurement. Until the strips are actually run through the new estimators, the claim remains an expectation.

C3
Verdict: STANDS
Keeping the profile gate (A4) and the meso gates (M1–M6) in non-overlapping lag bands prevents the classic error of averaging incompatible scales into a single verdict. The separation is required by the measured break points.

(1) Counterarguments the design missed

Finite-box bias is not the only bias; the structure-function estimator itself has lag-dependent variance that grows with the number of patch boundaries inside the 0.6 m box.
Real armour beds exhibit long-range imbrication coherence and subsurface packing that no surface-only Matérn modulation can reproduce.
The 8 % substrate contribution measured in earlier decompositions is an upper bound only under the particular hole-fill rule then in use; different fills move that fraction by >0.1 in H.
Optical relevance is never quantified: a bed that matches the profile gate can still produce caustic flash-rate errors >30 % if the high-frequency grain roughness is wrong.

(2) Failure conditions

Any generated strip whose directional structure-function profile falls outside the calibrated z-score envelope of pure H=0.79 fBm in the 0.6 m box.
Meso-band σ_B/σ_total or variogram range that remains statistically indistinguishable from a pure grain bed after the patch modulation is applied.
Cluster coverage or armour-ratio histograms that deviate from the measured 5–34 % / 1–4.7 ranges once the wake-sheltering rule is active.
Caustic intensity histograms under the new bed that differ from those measured on a real armour DEM by more than the already-documented morphology-induced 30 % reflectance bias.

(3) Mechanisms it ignores

Force-chain locking and true kinematic imbrication that control contact-angle distributions.
Subsurface filtration and vertical grading that set the porosity profile seen by refracted rays.
Flow-history dependence of patch wavelength (the same D50 can produce different correlation lengths under different formative Shields numbers).
The fact that Nikora's H≈0.79 is a scaling-region exponent at lags ≲ D50; demanding its recovery at 6–20 cm remains a category error even after the gate redesign.

(4) TOP 3 fixes ranked by impact

Obtain at least one published high-resolution armour-bed DEM (Bertin or equivalent) and run the identical structure-function estimator on it before any further parameter is frozen; use the real profile, not pure fBm, as the primary reference.
Execute the 3-D grain killer test (F-B4) on a minimal grain-only realisation first; if effective grain size cannot be driven below ~3 mm without destroying the measured σ_z/D50, the entire grain-primary architecture is falsified and must be abandoned.
Replace the provisional A4 thresholds with thresholds measured on an ensemble of pure fBm realisations of identical size and boundary conditions; only then re-open the calibration suite.

Q1
Yes. A pure fBm reference systematically under-states the knee and the meso-scale saturation present in real beds. "Ideal field" should be redefined as the ensemble-averaged structure-function profile of measured armour DEMs processed with the identical ruler, not a synthetic spectrum.

Q2
They are merely placed. The acceptance criteria (≥7/8, ≤1/8) have no measured Type-I/II rates behind them and are therefore not strict enough to certify the gate.

Q3
Two-sided testing is correct for a profile-matching gate. It is consistent with "reproduce nature" only if the reference itself is taken from nature, not from an idealised spectrum that real beds do not match at every lag.

Q4
Not defensible as stated. Powell's range is one data point; Bertin's recurrence is a spectral peak, not a Matérn correlation length. Scaling by D50 without a formative-flow argument is free-parameter tuning.

Q5
No evidence is given that an area-weighted mixture of local curves raises the global scaling-band slope above the grain-layer floor (~0.34). The mixture can at best produce a piecewise spectrum; the cheap F-B0 check is insufficient and the 3-D grain killer must come first.

Q6
Illegitimate. Hiding/exposure is a transport and entrainment concept. Mapping it onto a static placement rule without a measured depositional process is precisely the "looks-like" tuning the project claims to reject.

Q7
Yes. The killer test against a real DEM must precede bed_gen3. Obtaining Bertin's (or equivalent) DEM is higher priority than the synthetic calibration suite; the suite can only characterise estimator bias, not biological or sedimentary reality.

Q8
M3's lower bound 0.08 sits inside the grains-only noise floor (0.079) and is therefore too low to be discriminative. M6's D50 scaling inherits all the uncertainty of the surface D50 measurement. M5's 2.5×D90 separation moves by >20 mm when D90 ranges 75–95 mm; that sensitivity must be quantified before the band is locked.

---

(End of Grok's reply body. The follow-up suggestions shown on the Grok page — "Explore force-chain locking effects", "Investigate subsurface filtration roles" — are UI suggestions, not part of the reply.)
