# RALLY-V4-2-7-GPT (verbatim, ChatGPT)

Thread: https://chatgpt.com/c/6af5c4f8-c6ac-83ee-b6a0-e7dfb275f5a2 (bed topic, continued; thread title "Review grain architecture")
Prompt: RALLY-V4-2-7-PROMPT.md (verbatim in the task)

---

A. Gate-form redesign
A1 — Abolish the old STAGE 1 gate

Verdict: STANDS

The old gate is mathematically miscalibrated. The supplied ruler-ceiling experiment establishes that an ideal H_gen = 0.79 field can produce a measured local minimum of 0.65 in the 0.6 m box. Requiring every local estimate to exceed 0.70 therefore rejects the very reference population it is supposed to accept. Conversely, requiring an artificially high generating exponent (≥ 0.92) rewards the wrong physical spectrum. This is a demonstrated false-rejection/false-acceptance problem, not merely a preference about thresholds. Abolishing this gate is justified. However, abolishing it does not establish that the replacement gate is valid.

A2 — Compare directional H_loc profiles against matched reference profiles

Verdict: STANDS

This is the correct direction of correction because it compares like with like: same box, same directional estimator, same lag range, and same measurement procedure. It avoids treating a finite-box estimator as an unbiased measurement of the generating exponent. The strongest limitation is that a reference profile is only meaningful if the reference field represents the correct physical scaling model. Matching a synthetic reference does not establish matching a riverbed. In particular, a pure fBm reference can produce a misleading acceptance region if the real bed has a knee, multiple scaling regimes, or anisotropic spectral structure. The reference must be justified independently of the generator being tested.

A3 — Fixed fBm reference construction

Verdict: WEAKENED

The reproducibility provisions are useful, but the choice of pure fBm, periodic boundaries, and replication factors R = 8 and R = 4 is not yet physically justified. A PSD proportional to k^-(2H+2) describes an ideal two-dimensional self-affine surface over a specified scaling range. It does not describe a surface with a measured scaling break near D50, meso-scale saturation, and potentially different directional correlation lengths. Periodicity also introduces artificial long-range correlation and can contaminate large-lag estimates. Increasing replication reduces some sampling variability but does not remove model misspecification. A reference should represent a defined physical hypothesis, not simply the easiest synthetic field to generate.

A4 — Provisional z-score thresholds and three-part AND gate

Verdict: WEAKENED

The design correctly labels its thresholds provisional, but the numbers 1.0, 3.0, and 0.10 currently have no demonstrated statistical interpretation. A mean absolute z-score threshold of 1.0 can conceal systematic regional deviations; a maximum threshold of 3.0 depends heavily on the number of lag-direction observations; and the profile-drop condition can be redundant with the z-score conditions. More seriously, lag estimates are correlated. Treating the number of lag bins as the number of independent observations would exaggerate the evidence. An AND gate can also become arbitrarily strict as the number of tested lags increases. The thresholds need empirical operating-characteristic calibration, including false acceptance and false rejection rates, rather than being selected because they look conservative.

A5 — Estimator-bias calibration suite

Verdict: WEAKENED

Calibration before application to generated beds is essential, but the acceptance criteria are insufficient. Requiring 0.79 fields to pass at least 7/8 and 0.70/0.92 fields to pass at most 1/8 is a small-sample demonstration, not a reliable error-rate estimate. With eight trials, one field changes the observed pass fraction by 12.5 percentage points. Furthermore, the calibration only tests three generating exponents against the proposed gate. It does not test different random seeds, spectral cutoffs, anisotropy, knees, non-Gaussian elevation distributions, or nonstationarity. Existing grain beds at 0/3 are useful negative controls, but three failures cannot establish specificity.

A6 — Demote H_fit = 0.70–0.85 to a coarse screen

Verdict: STANDS

The supplied decomposition establishes why this is necessary: approximately 100% of the in-band variance comes from grain relief, whose measured exponent is approximately 0.34. The fBm substrate contributes only approximately 8%. Therefore, expecting the composite bed to recover the substrate's generating exponent is a category error. The profile is the more informative measurement. However, the coarse screen must not become a mechanism for dismissing contradictory evidence: profile-PASS with H_fit outside the band should remain unresolved, as proposed, until the cause is independently established.

A7 — Two-sided PASS/FAIL

Verdict: STANDS

A physical reproduction gate must reject both excessive roughness and excessive smoothness. The old one-sided threshold could reward a surface with an unrealistically high exponent, even though the ideal reference is 0.79. Two-sided testing is therefore consistent with the stated objective. The qualification is important: too-smooth and too-rough failures must be defined relative to a physically justified reference distribution, not merely relative to an arbitrarily selected fBm. A two-sided gate with a wrong reference is still a wrong gate.

A8 — PASS means reference-profile match, not real-bed validation

Verdict: STANDS

This distinction is scientifically necessary. A profile match demonstrates compatibility with a specified statistical reference, not reproduction of nature. The supplied measurements do not establish that H_loc profiles uniquely identify real water-worked gravel beds. Passing this gate cannot establish grain packing, transport history, cluster topology, pore connectivity, or optical equivalence. The wording must remain narrow: "PASS: compatible with the specified reference profile under the registered measurement procedure." Calling it a natural bed would exceed the evidence.

B. Patch-scale organisation
B1 — Matérn ν=1.5, 40–50 cm correlation, no paving

Verdict: WEAKENED

Dropping paving is correct because it artificially suppresses intra-patch variability rather than reproducing a formation process. Matérn ν=1.5 is a reasonable covariance choice for a differentiable correlated field, but the choice is not established by the cited riverbed evidence. The 40–50 cm streamwise scale is also not independently established as a universal patch scale. The measured approximately 200 mm variogram range and the proposed 450 mm correlation length are not necessarily contradictory: a Matérn correlation parameter and an empirical variogram range are different quantities. But the conversion between them must be calculated, not assumed. The design also needs to establish that the covariance model does not introduce artificial smoothness at grain-to-meso transitions.

B2 — Concrete patch parameters and alpha modulation

Verdict: WEAKENED

The parameters are testable, which is preferable to unrestricted artistic tuning. However, α = 0.64 is a measured variance-partition target, not proof that multiplying grain sizes or relief by a patch field produces the correct bedform component. The feasibility guard α_max ≈ 0.70 at σ_g = 2.1 is conditional on the assumed generator and should not be treated as a physical upper bound. Also, R_a = 2.0 m^γ needs a precise definition of m, its normalization, and the mapping from armour modulation to actual surface grain-size distributions. A parameter can match the variance ratio while generating the wrong spatial organization.

B3 — Hiding/exposure as a secondary wake-sheltering variant

Verdict: WEAKENED

Correctly demoting hiding/exposure is important. Wilcock and Crowe (2003) formulate hiding as a transport-threshold effect relative to surface grain-size distribution, not a universal static spatial-placement rule. The proposed wake-sheltering variant risks turning a hydraulic interaction into a visual placement heuristic. Fine sediment deposition downstream of coarse particles is physically plausible, but the position, probability, and size distribution depend on local flow and sediment supply. A wake-only rule is not enough to reproduce transport sorting.

B4 — bed_gen3.py composition architecture

Verdict: WEAKENED

The composition sketch contains a serious coupling problem. Patch field P modifies grain properties; support/deposition is then scaled by local D; clusters are biased toward coarse patches; and σ_form is re-iterated. These operations are not independent. Changing grain sizes changes packing, exposed area, support geometry, pore volume, and local elevation statistics. Changing deposition changes cluster probability and local grain-size distributions. Re-iteration can therefore converge numerically while producing a statistically inconsistent bed. The design needs a defined causal order and a conservation or packing constraint, not just an iterative adjustment until target statistics are met.

B5 — Honest expected movement: meso closes, scaling remains 0.35–0.50

Verdict: STANDS

This is one of the most defensible claims because it refuses to promise an impossible result. The measured decomposition says that the grain layer dominates the scaling-band variance. The supplied experiments show that changing support, imbrication, GSD width, and substrate amplitude has not solved the grain-scale exponent problem. Expecting meso-scale correlation to improve the meso band without magically raising the grain-dominated scaling exponent is consistent with that evidence. The numerical prediction of 0.35–0.50 is nevertheless a hypothesis, not a result. The cheap mixture calculation must be checked against the actual spatially varying grain field and its measured structure functions.

B6 — Six falsifiers, including the 3-D grain killer test

Verdict: STANDS

Explicit falsifiers are essential because the previous 22 combinations demonstrate how easily a generator can accumulate mechanisms without solving the measured discrepancy. The 3-D grain killer test is particularly important: a 2-D elevation construction can conceal unrealistic particle overlap, burial, and support geometry. However, the test's effective grain-size threshold of 3.0 mm needs justification relative to D50 = 36.2 mm. If 3.0 mm refers to voxel resolution or a minimum resolvable grain dimension, it is a numerical resolution requirement, not evidence of physical adequacy. The falsifiers must also be independent: using the same data to tune and kill a mechanism weakens their evidential value.

B7 — Scaling and meso effects live in different boxes

Verdict: STANDS

This is methodologically sound. A 0.6 m patch and a 3.0 m strip have different finite-size effects and support different lag ranges. Combining their results into one statistic would confound scaling behavior with meso-scale correlation and finite-box bias. The supplied ruler-ceiling measurements directly support box-specific reference profiles. The remaining weakness is boundary treatment: a 3.0 m strip must have enough independent meso structures to estimate its range and saturation reliably. A large nominal domain does not guarantee adequate independent spatial samples.

B8 — Textural patches as a required component of natural-bed organisation

Verdict: WEAKENED

Textural patches are supported by the gravel-bed literature, but the existence of patches does not establish that a Gaussian or Matérn random field is an adequate generator of them. Natural patches are related to hydraulic sorting, sediment supply, local topography, and transport history. Their boundaries need not be smooth, isotropic, or stationary. The Bertin, Groom, and Friedrich work supports the importance of separating grain and bedform roughness, but also demonstrates that their contributions are not simply additive. A field that produces spatially correlated grain-size variation is not automatically a model of textural patches.

B9 — Clusters as a required component

Verdict: WEAKENED

The evidence strongly supports cluster inclusion, but the design risks treating coverage as the principal validation target. Heays et al. (2014) found maximum cluster coverage between 5% and 34%, with only approximately 22% of potential anchor stones forming clusters. Their work also identifies obstacle, stoss, and wake regions with different particle-size populations. A cluster generator that simply biases cluster anchors toward coarse patches can match coverage while failing the internal organization. The cluster is a structured sedimentary arrangement, not a collection of coarse particles grouped by proximity.

B10 — Meso bedforms as a distinct roughness scale

Verdict: STANDS

The distinction between grain roughness and larger bedform roughness is supported by direct DEM analysis. Bertin et al. (2017) separated grain and bedform roughness using moving-window detrending and found that coarse sediment controls bedform development through surface humps, with finer sediment sheltered in hollows. This supports the separate meso-scale gate. The qualification is that the decomposition itself is scale-dependent: moving-window detrending can leak energy between bands, and the resulting variance fractions depend on the chosen window and filtering method.

B11 — Hiding/exposure correlation as an inventory item

Verdict: WEAKENED

Hiding and exposure are real physical phenomena, but the inventory risks conflating three different things: geometric exposure, hydraulic entrainment thresholds, and post-transport spatial sorting. Wilcock and Crowe (2003) explicitly formulate hiding in relation to fractional transport rates and the surface grain-size distribution. A spatial correlation between grain size and local elevation or sheltering is not by itself evidence of a correct hiding mechanism. The model needs a transport state or an empirically validated static equivalent.

B12 — Scale-hierarchy table

Verdict: WEAKENED

A hierarchy table is useful for preventing grain, cluster, and bedform scales from being conflated. But a table of nominal scales is not a demonstrated spectral decomposition. The same physical structure can contribute to multiple bands, particularly when grain-size distributions are broad and cluster dimensions overlap the grain and meso scales. Powell et al. (2016) explicitly describe scaling, transition, and saturation regions, and warn that spectral overlap complicates the definition of variogram range. Therefore, the hierarchy must specify operational filters, overlap, and uncertainty rather than presenting each scale as a cleanly isolated physical object.

C. Meso-gate validation
C1 — Calibration suite first, then recenter bands

Verdict: WEAKENED

Calibration must precede acceptance testing, but recentering the bands on calibration readings introduces a risk of moving the target to accommodate the instrument. Calibration should estimate measurement bias and uncertainty, not redefine the physical reference population. If the calibration result differs substantially from the literature-derived target, the correct response is to investigate the discrepancy, not automatically move the acceptance interval. The design needs a distinction between a fixed physical target, an estimator correction, and a calibrated acceptance tolerance.

C2 — Measure existing v4-2.2/2.3 strips, predict M3 false-positive and M4/M5 failures

Verdict: STANDS

The proposed retroactive test is scientifically useful precisely because it does not award historical beds a pass merely because a new gate has been designed around them. M3-in-band-for-the-wrong-reason is a particularly important diagnostic: a scalar metric can be correct for the wrong physical composition. M4/M5 failures would also establish whether the old generators reproduce the meso structure or merely happen to land inside a numerical interval. The qualification is that these beds are not independent ground-truth examples. They are useful negative and diagnostic controls, not a substitute for measured reference DEMs.

C3 — A4 and M1–M6 are complementary, non-overlapping bands

Verdict: WEAKENED

Separating the gates is correct if they measure distinct physical scales and use distinct lag domains. But "non-overlapping" must be demonstrated in the actual measurement operators, not declared in the design document. Grain, meso, and profile statistics can share variance through filtering, detrending, and spectral leakage. If A4 and M1–M6 are calculated from the same elevations using filters with overlapping transfer functions, their verdicts are statistically dependent. They must not be presented as independent confirmations merely because their names and nominal lag bands differ.

(1) Counterarguments the design missed

1. Identifiability is still unresolved. Multiple physically different beds can satisfy the same profile, variance ratio, and grain-size statistics. The design is improving the gate but has not demonstrated that the inverse problem has a unique or sufficiently constrained solution.

2. A static DEM cannot validate transport history. Armour formation, cluster survival, and hiding/exposure are consequences of sediment transport. Matching the final spatial distribution does not prove that the proposed formation mechanism is correct.

3. The covariance model may be doing too much work. Matérn ν=1.5 imposes a particular smoothness structure. Real bedforms can be nonstationary, directional, and multimodal. A good fit to a variogram does not establish the correct physical mechanism.

4. The variance decomposition is method-dependent. Bertin et al. show that patch-scale vertical roughness cannot simply be calculated as the sum of grain and bedform roughness. The weighting and decomposition procedure matter.

5. A DEM is not a complete bed. It measures the exposed surface, not buried grain geometry, pore connectivity, or the three-dimensional support structure.

6. Statistical agreement does not imply optical agreement. Two beds with similar elevation structure functions may have different facet-normal distributions, shadowing, and caustic energy distributions.

7. The reference population is not yet defined. A single H = 0.79 fBm reference may not represent the variability of water-worked gravel surfaces across different hydraulic conditions, grain distributions, and patch sizes.

(2) Failure conditions

- A2–A5: Synthetic H = 0.79 fields pass, but real water-worked DEMs systematically fail
- A3: Periodic-boundary artifacts alter the largest tested lags
- A4: Pass probability changes substantially when the number of lag bins changes
- A5: Calibration works for isotropic fBm but fails for anisotropic or knee-broken fields
- B1: Matérn fields reproduce the variogram but not the directional 2-D structure function
- B2: Target variance fraction is achieved while local grain-size distribution is wrong
- B3: Fine-grain wake placement does not match observed stoss/wake populations
- B4: Re-iteration converges but violates porosity, packing, or elevation constraints
- B5: Actual mixture-of-cutoffs calculation predicts no meaningful increase in meso-band slope
- B6: 3-D grain test shows overlapping or unsupported grains despite passing 2-D metrics
- B7: Strip contains too few independent meso structures for reliable range estimation
- B8–B12: Inventory statistics match individually but joint spatial statistics fail
- C1: Calibration recentering shifts acceptance bands away from independent field measurements
- C2: Existing beds pass because their measurement errors resemble calibration errors
- C3: Filtering leakage makes supposedly independent gates strongly correlated

(3) Mechanisms the design ignores

1. Hydraulic forcing and sediment supply — The current patch model is principally a spatial stochastic generator. It lacks a defined flow field, shear-stress distribution, transport rate, and sediment supply. Those variables control where particles move and deposit.

2. Particle contact mechanics — Support/deposition scaling is not equivalent to mechanically stable packing. A grain must have a plausible support configuration, contact network, and resistance to displacement.

3. Cluster formation and disintegration — Heays et al. (2014) demonstrate that only a minority of potential anchor grains form clusters and that cluster coverage changes with flow and time. A coarse-patch-biased point process does not reproduce this selection mechanism.

4. Transport-induced anisotropy — The measured anisotropy may differ across scaling, transition, and meso regimes. A single directional covariance parameterization may conceal this.

5. Surface/subsurface exchange — Armouring is not just a larger surface D50. It is a dynamically selected surface distribution related to subsurface composition and selective transport.

6. Measurement resolution and filtering — DEM resolution, interpolation, detrending, and hole filling all affect the estimated structure functions. The observed hole-fill sensitivity of 0.134 is already large enough to change a verdict.

7. Hydraulic feedback from bedforms — Bedforms change flow resistance and local shear stress, which in turn influence transport and bedform development. A one-way spatial modulation does not capture this feedback.

(4) TOP 3 fixes ranked by impact
1 — Highest impact: Obtain an actual water-worked gravel-bed DEM and make it the external reference

The most serious remaining weakness is that the design validates itself against synthetic reference fields without demonstrating that those fields represent measured riverbeds.

Required work:

Obtain a high-resolution, water-worked gravel-bed DEM with known grain-size distribution.

Preserve the original DEM and its measurement resolution.

Calculate directional 1-D and 2-D structure functions.

Calculate variograms, grain/bedform variance decomposition, and surface-normal distributions.

Apply exactly the same estimator and preprocessing pipeline to synthetic and measured data.

Keep the measured DEM out of generator calibration.

Without this, the design can demonstrate synthetic statistical consistency but not reproduction of nature.

2 — High impact: Run the 3-D grain killer before bed_gen3

Do not spend time building the full generator before testing whether the grain geometry itself can satisfy the measured constraints.

Test:

Realistic non-spherical grain geometry.

Grain contact and support.

Packing and porosity.

Exposed-surface elevation statistics.

Directional 2-D structure functions.

Surface grain-size recovery.

Sensitivity to DEM resolution and hole filling.

The critical question is whether the grain layer intrinsically imposes an exponent around 0.34 regardless of the meso substrate.

If so, the patch generator cannot repair the fundamental scaling problem without introducing an artificial mechanism.

3 — Medium-high impact: Replace the arbitrary patch field with a constrained physical or empirical model

Keep Matérn as a candidate, not a presumed truth.

Compare:

Matérn ν=1.5.

An empirically measured covariance model.

A knee-broken spectral model.

A nonstationary patch model conditioned on coarse-grain anchors and measured bedform geometry.

Evaluate the models against the same independent DEM.

Do not optimize solely for α = 0.64. Require agreement in directional covariance, spatial grain-size distribution, bedform wavelength, and grain/bedform coupling.

(5) Q1–Q8: Direct technical answers
Q1 — Pure fBm or knee-broken reference?

The verdict changes if the physical reference changes.

A pure fBm is acceptable as a null/reference model for testing estimator behavior. It is not automatically the correct physical reference for a gravel bed.

The supplied evidence already identifies a scaling break near D50 and saturation near D90. Therefore, the strongest physical reference is a knee-broken model with independently specified regimes.

I would define two references:

- Pure fBm, H=0.79 — Estimator calibration and finite-box bias
- Knee-broken reference — Physical comparison with measured gravel-bed scaling

The ideal field should be defined by its generating spectrum, directional properties, boundary conditions, resolution, and finite domain—not merely by its nominal H.

Crucially, do not tune the knee-broken reference to make bed_gen3 pass. Derive its parameters from independent measurements.

Q2 — Are A4 thresholds calibrated?

No. They are placed provisionally.

The 1.0, 3.0, and 0.10 thresholds have not yet earned statistical meaning.

The 8-seed acceptance criteria are also weak. Seven successes out of eight do not establish reliable acceptance probability.

Improvements:

Increase independent seeds substantially.

Include intermediate exponents, not just 0.70, 0.79, 0.85, and 0.92.

Include anisotropic and knee-broken fields.

Include different spectral cutoffs and boundary conditions.

Measure false acceptance and false rejection.

Test sensitivity to lag-bin count.

Use seed-level bootstrap confidence intervals.

One particularly important correction: calibrate the maximum absolute z-score against the distribution of the maximum over the entire profile. Do not interpret each lag's z-score independently.

Q3 — Is two-sided testing correct?

Yes.

The old one-sided gate rewards excessively smooth fields. That contradicts the objective of reproducing measured gravel-bed structure.

However, two-sided testing is only meaningful when the reference distribution is physically justified.

A bed can be too smooth, too rough, or have the wrong scale-dependent transition even when its overall variance is correct.

The appropriate test is a two-sided profile comparison, not simply a two-sided H threshold.

Q4 — Is 40–50 cm defensible?

WEAKLY defensible as an experimental sweep; not defensible as an established physical correlation length.

The distinction between correlation length and variogram range is critical.

Powell et al. (2016) distinguish the scaling region, transition region, and saturation region. They explicitly caution that the upper transition/range estimates can be highly variable and that the definition is affected by scale overlap.

Bertin et al. (2017) show that grain and bedform roughness are coupled and that their contributions cannot simply be added.

Therefore:

200 mm empirical range does not automatically imply 200 mm Matérn length.

A 450 mm Matérn parameter cannot be justified merely by multiplying D50.

The conversion depends on the exact covariance parameterization and range definition.

My recommendation: sweep 20–50 cm, but make the measured DEM determine the accepted range.

Q5 — Does the mixture effect actually raise the scaling-band slope?

Not automatically. This is a major unresolved mathematical assumption.

Suppose local patches have structure functions:

S_i(ℓ) = A_i ℓ^(2 H_i) f_i(ℓ / L_i)

where L_i is the local cutoff scale.

Ignoring boundary and mean-offset effects temporarily, the area-weighted structure function is:

S_mix(ℓ) = Σ_i w_i S_i(ℓ)

The effective exponent is:

H_eff(ℓ) = (1/2) d log S_mix(ℓ) / d log ℓ

For a target exponent H*, the weighted mixture must satisfy:

[Σ_i w_i S_i(ℓ) H_i(ℓ)] / [Σ_i w_i S_i(ℓ)] ≈ H*

throughout the entire required lag interval.

This is the key condition. It is not enough for local D50 values to vary. The area-weighted contributions of the local structure functions must produce the required effective slope.

There is an additional complication: a spatially nonstationary mixture contains cross-patch pairs. The actual structure function is not necessarily the simple weighted average above. Patch boundaries and differences in local mean elevation can introduce extra terms.

Which test first? F-B0, then F-B4.

F-B0 is cheap and should first determine whether the proposed mixture of cutoff scales can mathematically produce the target slope.

F-B4 then tests whether actual 3-D grains, packing, and spatial organization preserve that result.

But F-B0 cannot authorize bed_gen3. It is only a necessary preliminary plausibility test.

If the mixture calculation predicts no increase beyond 0.35–0.50, there is no reason to build the full generator around that mechanism.

Q6 — Is wake-based hiding/exposure placement legitimate?

Not as currently described.

There is a physically plausible mechanism, but the placement rule is not equivalent to the transport process.

Heays et al. (2014) provide particularly relevant evidence: typical clusters have an anchor, stoss, and wake, with different grain-size populations. They also found that only about 22% of potential anchor stones actually formed clusters.

That supports a conditional deposition mechanism, not a universal rule that fine grains belong behind coarse grains.

Wilcock and Crowe (2003) model hiding through size-dependent transport thresholds referenced to the surface distribution.

The missing causal chain is:

Coarse grain creates local flow disturbance.

Local hydraulic conditions modify entrainment and transport.

Incoming sediment has a size-dependent probability of deposition.

Deposited particles modify subsequent flow and sheltering.

Without at least a calibrated representation of that chain, "fine grains in coarse-grain wakes" is a spatial prior, not a demonstrated formation mechanism.

It can be used as a secondary hypothesis, but must not be counted as evidence of physical reproduction.

Q7 — Should the 3-D killer precede bed_gen3? DEM or calibration first?

Yes, the killer should precede the full bed_gen3 build.

My ordering:

Acquire the Bertin/Powell DEM

Secure the actual data, resolution, coordinate system, and grain-size metadata. A paper's plotted variogram is not an adequate substitute for the original elevation data.

Run synthetic estimator calibration

This can proceed immediately, without waiting for the DEM. Validate finite-box bias, periodicity, lag dependence, and z-score behavior.

Run the 3-D grain killer

Use the acquired DEM as an external benchmark. Test whether realistic grain geometry can reproduce the observed grain-scale signatures.

Run the mixture-of-cutoffs test

Determine whether meso heterogeneity can mathematically produce the expected effective slope.

Authorize or reject bed_gen3

Only proceed if the grain layer and meso mechanism each pass their own necessary tests.

Important distinction: obtaining the DEM should happen before the killer, but synthetic calibration does not need to wait for DEM acquisition.

The source literature identifies relevant DEM datasets, but I cannot establish from the available publication records that the original raw DEM files are publicly downloadable.

Q8 — M3, M6, and M5 sensitivity

M3: 0.08 is too close to 0.079.

The difference is:

0.080 − 0.079 = 0.001

That is approximately 1.27% relative to the grains-only value.

This is not a defensible separation unless measurement and estimator uncertainty are substantially below 0.001 in the relevant metric.

A grains-only bed could therefore pass M3 for the wrong reason. The proposed prediction of M3-in-band-for-the-wrong-reason is credible, but the lower bound needs recalibration.

M5: highly sensitive to D90 uncertainty.

Given D90 = 75–95 mm:

2.5 D90 = 187.5–237.5 mm

Total threshold uncertainty: 50 mm, approximately 24% of the midpoint.

That is large enough to change the classification of a bed whose transition lies near the boundary.

Use a measured D90 per bed, propagate its uncertainty, and avoid treating 2.5D90 as a universal physical separation.

M6: D50 scaling is conditionally defensible.

Normalizing a grain-related scale by D50 is useful for comparing beds with different characteristic grain sizes. But it is not automatically valid for meso structures.

Meso bedforms may scale with D84, D90, hydraulic forcing, channel geometry, or sediment supply rather than D50 alone.

The decisive test is whether normalized M6 distributions collapse across independent beds with substantially different D50 values.

If they do not, D50 normalization is an imposed convenience rather than a demonstrated physical law.

I cannot assign a stronger verdict to M6's exact numerical boundary without its full definition and formula; the design summary only identifies it as a D50-scaled gate.

---

Verdict tally from the reply: STANDS — A1, A2, A6, A7, A8, B5, B6, B7, B10, C2 (10 claims); WEAKENED — A3, A4, A5, B1, B2, B3, B4, B8, B9, B11, B12, C1, C3 (13 claims); WRONG — none.
