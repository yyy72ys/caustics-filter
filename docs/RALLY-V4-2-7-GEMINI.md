# RALLY-V4-2-7-GEMINI (verbatim, AI Studio, Gemini 3.8 Flash, thinking Medium)

Thread: https://aistudio.google.com/prompts/17kvdsfgdvfZxATYk3CzMa2s1XNmNo-eQ (bed topic, continued)
Prompt: RALLY-V4-2-7-PROMPT.md (verbatim in the task)

---

Claim-by-Claim Verdicts
(A) Gate-Form Redesign (Claims A1–A8)

Claim A1 (Abolish old STAGE 1 gate: local H ≥ 0.70 everywhere)

Verdict: STANDS

Reasoning: The old gate exhibits an intolerable 100% Type I error rate against its own analytical target. Because second-order structure functions (D2(Δx) = ⟨[z(x+Δx) − z(x)]²⟩) evaluated inside a finite window of length L suffer from boundary truncation, unrepresented long wavelengths, and sample variance, an ideal synthetic fractional Brownian motion (fBm) surface with input H_gen = 0.79 systematically yields an apparent local exponent H_loc(Δx) ≈ 0.65 at the upper lag boundary in a 0.6 m box (finite-box attenuation ΔH ∼ 0.10–0.14). Demanding H_loc ≥ 0.70 across the entire 0.4–3.6 cm band forces the acceptance of unphysical surfaces (H_gen ≥ 0.92) that violate the Nikora et al. (1998) scaling laws.

Claim A2 (New gate = profile-vs-profile comparison against ideal H_gen = 0.79 fBm)

Verdict: WEAKENED

Reasoning: An ideal fBm surface is mono-affine and scale-invariant across all spatial frequencies (S(k) ∝ k^−(2H+1)), whereas natural water-worked gravel beds possess two distinct physical scale breaks: a micro-scale grain-geometry regime (H ≈ 0.75–0.85 at lags Δx < D50), a structural transition regime (D50 < Δx < D90), and a saturation regime (H → 0 at Δx > D90) (Nikora et al., 1998; Aberle & Nikora, 2006). Comparing a natural bed's profile against an un-broken, scale-free fBm reference across the entire lag domain forces the generator to match a mathematical abstraction rather than fluvial reality, penalizing the bed for exhibiting natural clast-scale saturation.

Claim A3 (Fixed reference profile construction: spectral synthesis, PSD ∼ k^−(2H+2), identical estimator)

Verdict: STANDS

Reasoning: This is a methodologically necessary baseline protocol. Enforcing identical grid resolutions, boundary conditions (periodic wrapping vs. mirror padding), and estimator execution paths between synthetic reference fields and candidate beds eliminates algorithmic leakage, discretization bias, and windowing artifacts from the statistical comparison.

Claim A4 (Comparison metric = per-lag z-scores with C1 mean |z| ≤ 1.0, C2 max |z| ≤ 3.0, C3 drop ≤ 0.10)

Verdict: WEAKENED

Reasoning: Structure function estimates across consecutive lags Δx_i and Δx_{i+1} are non-independent; their errors are bound by a dense, off-diagonal covariance matrix governed by the spatial autocorrelation of the surface. Treating per-lag deviations as independent univariate z-scores via mean and max thresholds ignores this covariance structure, allowing systematic low-frequency profile drifts or oscillatory distortions to pass C1 and C2 as long as their amplitude does not breach arbitrary point thresholds.

Claim A5 (Estimator-bias calibration suite executed prior to candidate evaluation)

Verdict: STANDS

Reasoning: Pre-registering a statistical discriminator on known synthetic benchmarks (H ∈ {0.70, 0.79, 0.85, 0.92}) and requiring a ≥ 7/8 hit rate for H = 0.79 alongside ≤ 1/8 false-positive rate for adjacent steps enforces basic statistical power. It guarantees that the testing apparatus can actually falsify candidate beds before compute is squandered on generator sweeps.

Claim A6 (H_fit band 0.70–0.85 demoted to coarse screen; profile is primary)

Verdict: STANDS

Reasoning: Linear regression across an arbitrary lag interval on a log-log variogram collapses multi-scale curvature into a single scalar, obscuring inflection points, localized flattening, and non-monotonic behavior. Profile matching preserves local slope derivatives (∂ log D2 / ∂ log Δx), making it strictly superior to aggregate line fitting.

Claim A7–A8 (Two-sided testing; "PASS" denotes reference match, not ontological truth)

Verdict: STANDS

Reasoning: Fluvial beds are bounded on both sides: an unworked bed is too rough (H ≈ 0.45–0.55, screeded/uncorrelated), while an over-smoothed synthetic surface (H → 1.0) violates turbulent shear packing. Restricting the semantic scope of "PASS" to mean statistical indistinguishability from the reference profile under a pre-registered protocol prevents epistemic inflation and confirmation bias.

(B) Patch-Scale Organisation (Claims B1–B7)

Claim B1 (Drop paving, stretch correlation to 40–50 cm, switch to Matérn ν = 1.5)

Verdict: WEAKENED

Reasoning: Matérn ν = 1.5 (C(r) ∝ (1 + 3r/ℓ) exp(−3r/ℓ)) correctly provides once-differentiable surfaces, avoiding the unphysical infinite differentiability of squared-exponential kernels. However, stretching the correlation length ℓ to 40–50 cm within a 60 cm evaluation box creates an extreme domain starvation issue (L_box/ℓ ≈ 1.2–1.5). At this ratio, a single realization contains barely one correlation cell, causing the patch field to manifest as a trivial planar tilt or low-order dome across the entire box, destroying statistical stationarity and causing massive seed-to-seed variance.

Claim B2–B3 (Parameterization: L_x = 45 cm, L_y = 22 cm, α = 0.64, armour modulation, wake sheltering)

Verdict: WEAKENED

Reasoning: Patch-scale sediment sorting in rivers (Nelson et al., 2010; Buffington & Montgomery, 1999) is an emergent instability of bedload transport over mobile beds, governed by boundary shear stress ratios τ*/τ_c*. Imposing an ad-hoc continuous scalar field to deterministically scale discrete grain dimensions (D50, R_a) without hydrodynamic coupling forces an arbitrary spatial segregation pattern. Furthermore, parameterizing wake sheltering without solving the recirculating eddy length behind clasts (L_w ≈ 2–4 D_clast) produces artificial coarse-fine spatial correlations that do not match measured fluvial cluster tails (Heays et al., 2014).

Claim B4 (bed_gen3.py multi-stage serial composition)

Verdict: WRONG

Reasoning: The proposed pipeline strings together six serial, uncoupled heuristic stages: patch field → grain modulation → kinematic support/deposition → forced imbrication → biased cluster anchoring → form-scale variance iteration. Because kinematic settling and cluster insertion displace surrounding particles, each subsequent step destructively overwrites the geometric moments established by prior steps. In granular physics, packing, imbrication, and cluster stability are simultaneously coupled via contact force networks; treating them as an uncoupled open-loop sequence guarantees perpetual parameter drift and destabilizes both the target porosity (0.30–0.45) and the measured imbrication dip (25°–35°).

Claim B5 (Honest movement: meso band closes, scaling-band H_fit stays low at 0.35–0.50)

Verdict: STANDS

Reasoning: This is a mathematically sound concession. The v4-2.5 variance decomposition showed that in-band variance is ≈ 100% dominated by grain-scale surface relief (which possesses an intrinsic individual clast exponent H ≈ 0.30–0.34). Modulating grain size across patches shifts the local break lag D50(x, y), but as shown by superposition theorems for self-affine processes, mixing non-stationary cutoffs does not steepen the sub-break power-law slope; it only broadens the transition knee. Admitting that H_fit cannot reach 0.79 via patch modulation alone matches analytical reality.

Claim B6 (Six falsifiers F-B0–F-B5 pre-registered, including 3D grain killer test)

Verdict: STANDS

Reasoning: Pre-registering analytical and numerical barriers—specifically F-B0 (proving whether a mixture of cutoffs can mathematically raise H) and F-B4 (the 3D grain killer test)—prevents sunk-cost fallacies. If a simple mixture model cannot produce an apparent exponent ≥ 0.70, terminating bed_gen3.py before building complex generators is strictly compliant with empirical rigor.

Claim B7 (Scaling effects live in 0.6 m patch; meso effects live in 3.0 m strip)

Verdict: STANDS

Reasoning: Resolves domain contamination. Measuring meso-scale spatial structures (λ ≈ 100–200 mm) inside a 0.6 m box violates the Nyquist-like requirement for spatial statistics (L_box ≥ 5λ). Allocating the meso-band variogram to the 3.0 m strip (L/ℓ ≈ 15) and grain-scale micro-topography to the high-resolution 0.6 m patch correctly isolates the respective spatial support windows.

Real-Bed Inventory (Claims B8–B12)

Claim B8–B12 (Literature inventory: patches, clusters, bedforms, transport-based hiding/exposure, scale hierarchy)

Verdict: STANDS

Reasoning: The inventory accurately synthesizes foundational fluvial literature: textural patches (Buffington & Montgomery, 1999; Nelson et al., 2010), obstacle clusters and anchor clast ratios (Heays et al., 2014; Little & Venditti, 2014), and the physical distinction between geometric packing and transport hiding/exposure (Parker, 1990; Wilcock & Crowe, 2003). Framing these phenomena as distinct physical regimes organized by the scale-hierarchy table is geomorphologically accurate.

(C) Meso-Gate Validation (Claims C1–C3)

Claim C1 (Run calibration suite first; re-centre bands on calibration readings)

Verdict: STANDS

Reasoning: Calibrating empirical bounds against an explicit numerical baseline before testing candidate models prevents shifting goalposts during active development. It provides the only objective mechanism to separate estimator bias from actual physical bed properties.

Claim C2 (Audit existing v4-2.2/2.3 strips; predict M3-in-band-for-wrong-reason, M4/M5 fail; no retroactive pass)

Verdict: STANDS

Reasoning: Applying the new diagnostic suite to known failed baselines serves as negative control validation. If an old, discarded generator accidentally passes the new gates, the gates are degenerate; predicting specific failure modes (M4/M5 failure due to lack of spatial patchiness) pre-tests the discriminative power of the metrics.

Claim C3 (A4 and M1–M6 are complementary non-overlapping bands, never merged into one verdict)

Verdict: STANDS

Reasoning: Merging orthogonal spatial regimes into a single composite score or scalar metric enables compensation errors, where excellent performance in the meso-band masks total physical failure in the grain-roughness band. Evaluating them as strictly decoupled gates prevents dimensional cross-contamination.

Counterarguments the Design Missed

The Superposition Dilemma (Lowest H Dominates Small Lags): If a surface z(x) is constructed as the sum of a grain field z_g(x) (with H_g ≈ 0.33) and a patch field z_p(x) (with H_p ≈ 0.80), the second-order structure function at small lags Δx → 0 satisfies:

D2(Δx) = D2,g(Δx) + D2,p(Δx) ∼ A Δx^(2H_g) + B Δx^(2H_p)

Because H_g < H_p, as Δx → 0, the Δx^(2H_g) term asymptotically dominates the sum. You cannot increase the scaling exponent of a rough surface at small lags by superimposing a smoother, higher-H long-wavelength field; the rougher component invariably dictates the small-lag slope.

Ellipsoidal Rasterization Distortion: Representing clasts as analytical ellipsoids and rasterizing them onto a regular grid generates artificial, singular geometric cusps along intersecting clast boundaries. Natural water-worked beds do not have sharp geometrical intersections; they exhibit smooth fine-sediment fillets, sand infilling, and abraded grain edges. These artificial raster cusps introduce excessive high-frequency, uncorrelated step discontinuities that systematically depress H_loc at millimetric lags.

Hydrodynamic Clast-Top Shearing vs. Deposition Geometry: In real gravel beds, the anomalously high Hurst exponent (H ≈ 0.79) reported by Nikora et al. is not produced by clast shapes or spatial clustering; it is produced by hydrodynamic fluid shear acting on the bed boundary. Water flow "planes off" protruding clast summits through selective entrainment, while depositing fine material in local pockets, flattening high-elevation variance and creating positive elevation skewness (+0.5–+0.7). Pure kinematic settling without hydrodynamic surface truncation cannot reproduce this mechanism.

Failure Conditions Under Which the Design Breaks

Analytical Breakdown of F-B0: The F-B0 mixture-of-cutoffs will analytically fail to raise H_fit from 0.35 to ≥ 0.70. Convolving or area-weighting multi-scale structure functions with locally varying turnover lengths D50(x, y) broadens the knee between the scaling and saturation regimes, flattening the local slope rather than steepening it.

Domain Starvation in the 0.6 m Evaluation Box: With Matérn correlation lengths of L_x = 45 cm within a 60 cm window, individual seeds will be dominated by random single-patch realizations. When evaluating C1–C3 across a batch of 8 seeds, the sample variance will cause max |z| > 3.0 on at least 3 of 8 seeds, failing the calibration acceptance criterion (≥ 7/8).

Inter-Clast Voids Destabilizing Porosity: Because bed_gen3.py relies on kinematic deposition of coarse grains first, followed by ad-hoc cluster insertion, the void spaces between coarse particles in un-armoured patches will remain completely unfilled by matrix fines, causing local porosity to exceed the physical upper limit of 0.45, directly violating the bulk porosity gate.

Mechanisms It Ignores

Fluvial Abrasion and Grain Angularity: Natural fluvial gravels are not smooth mathematical ellipsoids; they are fractured, abraded clasts possessing flat cleavage facets, angular edges, and variable roundness (Wadell roundness 0.2–0.6). Clast faceting introduces localized planar orientation correlations that elevate H at sub-clast scales.

Pore-Throat Matrix Infiltration (Hydrodynamic Sieving): In natural streams, fine gravel and sand infiltrate downward through the coarse skeletal framework until reaching pore-throat bridging limits (D_pore ≈ 0.15–0.25 D_framework). The design treats matrix and framework as a combined single-layer kinematic deposit, completely ignoring gravitational and fluid-driven sieving.

Hydrodynamic Boundary-Layer Truncation (Topographic Planation): Turbulent boundary shear stress imposes a physical ceiling on clast protrusion. Protruding grains experience drag forces scaling with u²(z), mobilizing them until the bed surface forms a hydrodynamically stable, sheared armoring plane.

TOP 3 Fixes Ranked by Impact

Replace Smooth Ellipsoids with Faceted, Sheared Clasts & Implement Bed-Top Elevation Clipping (Impact: HIGHEST):

Action: Introduce angular/faceted mesh geometries (e.g., Voronoi-fractured particles) and apply a non-linear hydrodynamic shearing operator that numerically truncates/depresses the highest 5–10% of protruding clast summits while preserving positive skewness (+0.5–+0.7). This directly alters the sub-clast height variance ⟨Δz²⟩, bridging the gap between H ≈ 0.35 and H ≈ 0.79.

Replace Pure Synthetic fBm Reference in A2/A3 with Empirical Gravel-Bed Reference Curves (Impact: HIGH):

Action: Discard the scale-free H_gen = 0.79 fBm reference profile. Extract the baseline reference profile directly from open, calibrated laboratory and field laser-scan DEMs of water-worked gravel beds (e.g., Bertin et al., 2017 open data repository). Match the candidate generator against real physical structure functions that already contain the correct physical scaling breaks at D50 and D90.

Execute the F-B0 Analytical Proof as an Absolute Blocking Gate (Impact: MEDIUM):

Action: Halt all code development on bed_gen3.py until the F-B0 mixture-of-cutoffs derivation is complete. If the mathematical integration of spatially distributed structure functions proves that the aggregate slope cannot exceed max(H_local), immediately abandon patch-scale modulation as a mechanism for raising scaling-band H and pivot resources to clast-scale physics.

Answers to Questions Q1–Q8

Q1: Reference profile construction (fBm vs. knee-broken)?
Pure fBm changes the verdict of A2 to WEAKENED. A pure fBm field lacks physical clast turnover. The "ideal field" must be defined as a knee-broken profile with three distinct analytical or empirical zones: a scaling sub-range (H ≈ 0.79 for lags < D50), a transitional roll-off (D50 ≤ lag ≤ D90), and a flat saturation plateau (H → 0 for lags > D90).

Q2: A4 provisional thresholds: calibrated or merely placed? Calibration suite rigor?
They are merely placed. Thresholds like |z| ≤ 1.0 and |z| ≤ 3.0 are textbook Gaussian heuristics chosen without calculating the empirical covariance of D2(Δx). The calibration suite is not strict enough: requiring ≥ 7/8 on H = 0.79 while allowing 1/8 on 0.70 and 0.92 gives an unacceptably wide margin of error (≈ 12.5%), which is too coarse to distinguish H = 0.75 from 0.79.

Q3: Two-sided testing consistent with "reproduce nature"?
Yes. Natural gravel beds are neither uncorrelated white noise (H → 0) nor differentiable, frictionless slopes (H → 1). Water working organizes the topography into a specific self-affine window (H ∈ [0.70, 0.85]). Failing surfaces that are too smooth (H ≥ 0.90) is mandatory to prevent convergence toward non-physical, plastic-looking surfaces.

Q4: Is the 40–50 cm correlation length defensible?
No. Scaling it directly from Powell et al. (2016) or Bertin et al. (2017) ignores the ratio of domain size to correlation length. A 45 cm correlation length in a 60 cm box yields L/ℓ ≈ 1.3, rendering spatial statistics completely non-ergodic. It is only defensible within the 3.0 m strip (L/ℓ ≈ 6.7). For the 0.6 m patch, the maximum defensible correlation length is 10–12 cm (≈ 3 D50).

Q5: Mixture effect on scaling slope? F-B0 vs. F-B4?
It does not raise the slope. The area-weighted average D̄2(Δx) = ∫ D2(Δx; D50) p(D50) dD50 of curves that each have an asymptotic power-law slope of 2H ≈ 0.68 (H ≈ 0.34) cannot yield a composite slope of 2(0.79) = 1.58. The envelope of lower bounds remains strictly bounded by the dominant local power. The cheap mathematical check F-B0 must come first; running F-B4 before solving F-B0 is an inefficient use of engineering effort.

Q6: Translating hiding/exposure into wake-sheltering placement?
It is "looks-like" tuning without a physical formation process. In sediment transport, hiding/exposure describes the modulation of critical Shields stress (τ*_ci = τ*_c50 (D_i/D50)^−β) due to relative pocket depth and flow protrusion. Hardcoding fine grains into downstream geometric shadow zones without simulating fluid drag forces mimics the visual outcome while bypassing the hydrodynamics, making it fragile under statistical auditing.

Q7: Killer test vs. calibration suite vs. bed_gen3?
Obtaining Bertin's real-bed DEM comes first. You cannot calibrate an estimator or validate a generator against an ungrounded mathematical abstraction. The absolute operational sequence must be: (1) Obtain real-bed DEM → (2) Extract true empirical reference profile and covariance matrix → (3) Run calibration suite on empirical profile → (4) Execute F-B0 proof → (5) Only then consider building bed_gen3.

Q8: M1–M6 bands: M3 lower bound, M6 scaling, M5 separation sensitivity?
The M3 lower bound (0.08) is dangerously close to the unorganized grains-only baseline (0.079), rendering it incapable of discriminating between true meso-bedform organization and numerical noise. M6's D50 scaling is defensible because cluster spacing scales directly with boundary-layer thickness and grain obstacle size (3–5 D50). M5's separation threshold (2.5 D90) is hyper-sensitive to D90 uncertainty: varying D90 from 75 to 95 mm shifts the cutoff from 187.5 to 237.5 mm—a 50 mm swing that spans nearly half of the entire cross-stream variogram range.

---

Verdict tally in the reply: STANDS — A1, A3, A5, A6, A7–A8, B5, B6, B7, B8–B12, C1, C2, C3; WEAKENED — A2, A4, B1, B2–B3; WRONG — B4. The thread remains open in the retained tab.
