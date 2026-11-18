# RALLY-V4-3-1-GEMINI — 第3段設計ラリー Gemini 脚（verbatim）

- 脚: Gemini（Google AI Studio、モデル Gemini 4 Hydrogen / gemini-4-hydrogen、Thinking level: High — 選択を確認してから投稿）
- 対象: RALLY-V4-3-1-PROMPT.md（DESIGN-V4-3-1-DRAFT の主張 C1–C10 への敵対的査読）
- 日時: 2026-11-14 16:39–16:44 JST
- プロンプトURL: https://aistudio.google.com/prompts/1beUsA77lJffmNs5_1QIldixqmasoo9Qk
- 注記: 以下は返答本文の verbatim 保存。改変なし。

---

Verdicts on Individual Claims (C1–C10)

C1. GRIN bending is first-order on both the sun side and view side

Verdict: WEAKENED

Counterargument: Two-way GRIN is required in principle: sun-side and view-side rays traverse different regions of n(x) when the sun is oblique, and at near-horizontal interfaces (stratification, haloclines) vertical top-down view rays suffer near-zero Snell deflection while oblique sun rays refract. Moreover, sun-side GRIN acts as a refractive lens modulating caustic irradiance via ∇²n (shadowgraph focusing), not merely a shift. However, claiming sun-side GRIN is generally first-order for Stage 3's chosen first field (vertical bed-heated plumes) fails on scaling. A sun ray refracted into water at θw ≈ 20°–45° crosses a vertical plume of width w transversely over path w/sinθw rather than grazing along its full 0.5 m vertical height. For ΔT = 3 K (Δn ≈ 3×10⁻⁴), transverse crossing yields δθ ∼ (Δn/n) cotθw ≈ 0.2–0.6 mrad, shifting the caustic at 0.3 m below the plume by ∼0.06–0.18 mm — sub-pixel, and dwarfed by the ∼3.5–4 mm solar-disc blur and centimetre-scale surface-wave caustics.

Failure condition: On the main verification window (0.6 m depth) with solar zenith angle ≥30° and bed-heated plumes, toggling sun-side GRIN on/off produces <0.2 mm caustic shift and <1% change in deposit band energy after solar-disc quadrature.

C2. Bed-heated plumes produce 1–4 mm apparent displacement (several pixels), while cool-skin convection is sub-pixel

Verdict: WEAKENED

Counterargument: Using Quan & Fry (|∂n/∂T| ≈ 0.9×10⁻⁴ K⁻¹ at 20°C, n = 1.333), achieving apparent displacement δ = (1/2)(∇⊥n/n)L² = 1–4 mm over L = 0.5 m requires a sustained horizontal gradient ⟨∇⊥n⟩ = 0.011–0.043 m⁻¹, i.e., ∇⊥T ≈ 110–430 K/m (1–4 K/cm) maintained coherently along the entire 0.5 m vertical path. Real turbulent buoyant plumes from solar bed heating (ΔT ∼ 0.5–2 K at the bed) dilute rapidly with height (ΔT(z) ∝ z⁻⁵ᐟ³), widen by entrainment (b ∼ 0.1z), and meander laterally so a vertical view ray rarely stays inside a steep edge wall for 0.5 m. A realistic smooth plume (ΔT = 2 K, Gaussian σ = 0.08 m) yields δmax ≈ 0.17 mm — strictly sub-pixel, just like cool-skin plumes. The 1–4 mm figure only holds for unmixed, razor-straight vertical thermal sheets of ∼1 cm wall thickness or strong localized springs (ΔT ≳ 5 K).

Failure condition: Computing δ from physically entraining, meandering thermal plumes driven by realistic pond surface energy balance (≤500 W/m² bed heating) yields δ95% < 0.4 mm (<0.5 px) across the main verification window.

C3. A parametrised plume field from energy balance and bed inputs is admissible as a "labelled input"

Verdict: WEAKENED

Counterargument: External boundary conditions (bed albedo map, spring inflow locations, heat flux) are legitimate labelled inputs; synthesizing a 3-D volumetric temperature field T(x,t) from them without buoyancy/entrainment conservation laws is not. Unlike optical depth (an integral), GRIN schlieren depends on ∇T (a spatial derivative), which is directly controlled by the assumed plume profile shape and wall thickness b(z). Assigning "spacing of order the depth" (which is near-threshold laminar Rayleigh–Bénard spacing, not high-Ra ∼10⁸–10¹⁰ turbulent pond plumes that erupt from thin boundary layers) and choosing a profile width without Morton–Taylor–Turner plume scaling (b ∝ αe z, ΔT ∝ B²ᐟ³ z⁻⁵ᐟ³) leaves the entire schlieren signal amplitude governed by an unphysical free shape parameter — violating principle (b).

Failure condition: Two geometrically plausible profile shapes (e.g., σ = 0.01 m sharp sheet vs σ = 0.10 m Gaussian plume) satisfying the identical surface energy balance and bed inputs produce schlieren deflections differing by >5×, demonstrating that the result is dictated by arbitrary profile parameterization rather than physical inputs.

C4. Generalising to inhomogeneous RTE needs no new kernel shape — only evaluating v16 FF-CDF at local (τ, b/c)

Verdict: WRONG

Counterargument: In inhomogeneous media, the spatial halo PSF is not a function of path-integrated (τ, b/c) alone; it depends critically on the geometric lever arm (L − s) between the scattering depth s and the target plane L. By small-angle RTE (Wells / Dolin beam spread), the radial PSF scales with (L − s) tanθ. Compare two water columns with identical total (τ, b/c):

A turbid river plume in the top 20% over clear bottom water (block 3.3's priority scene!): scattering near the top has a ∼0.5 m lever arm to the bed, producing a broad, wide halo and strong bed blur.

A benthic nepheloid layer in the bottom 20% hugging the bed: identical (τ, b/c), but lever arm ≈ 0, producing almost zero spatial blur on the bed.

Evaluating the homogeneous v16 FF-CDF kernel at (τ, b/c) collapses both cases to the same halo, and also forces radial symmetry where transverse turbidity edges (plume fronts) produce asymmetric halos. This directly violates the governing principle by substituting a homogeneous-medium look-alike for depth-weighted scattering physics.

Failure condition: Running 3.2's own delta-tracking MC reference on two columns with equal (τ = 0.6, b/c = 0.8) — turbid top 20% vs turbid bottom 20% — will show a 3–5× difference in halo half-energy radius while the FF-CDF (τ, b/c) model predicts identical containment curves.

C5. Path radiance (veil) can be built now against the fixed top-down view, with Stage 4 only re-pointing the eye end

Verdict: WEAKENED

Counterargument: Formulating veil as a ray-marched integral along arbitrary view rays is sound in principle, but building and validating it only against a fixed top-down view hides three failures that break at Stage 4:

Phase-function regime mismatch: Top-down view looks roughly opposite the sun (θs > 110°, near the minimum backscatter plateau of the Fournier–Forand phase function), whereas tilted Stage 4 cameras look toward forward-scattering angles where βFF(θs) rises by 10²–10³×.

Missing 3-D volumetric caustic irradiance: H4 explicitly cites "visible light shafts in the water". Light shafts require evaluating the 3-D focused caustic irradiance field E(x) inside the water volume (driven by Stage 1 waves + 3.1 GRIN), not 1-D attenuated sunlight. Top-down views foreshorten shafts; Stage 4 oblique views expose them immediately.

Stage 1 view-ray refraction: For an above-water Stage 4 camera, view rays refract across the wavy Stage 1 surface — changing both path geometry and θs(x) pixel-by-pixel, not merely "re-pointing" at the camera.

Failure condition: Validating veil top-down to <5% error, then tilting the Stage 4 camera to 45° toward the sun and observing >30% radiance error or absent volumetric light shafts against MC.

C6. Delta-tracking Monte Carlo is the right verification reference (HydroLight rejected)

Verdict: WEAKENED

Counterargument: Rejecting 1-D plane-parallel HydroLight for 3-D finite patches + bed is physically correct. However, invoking a generic "delta-tracking (Woodcock) Monte Carlo tracer" as the reference overlooks two fatal numerical couplings in Stage 3:

Curved rays vs Woodcock stepping: Standard delta-tracking assumes straight free-flight segments in a homogeneous majorant cmaj. In 3.1's GRIN medium (and 3.3's sharp interfaces), rays curve continuously — requiring Woodcock null-collision sampling integrated directly inside the RK4 ray-marching step (and boundary intersection at 3.3 steps).

Catastrophic variance / broken NEE for the veil & sun-bed-eye paths: You cannot connect a scattering vertex inside water to a 0.5° sun disc using a straight shadow ray when light must cross a wavy refractive Stage 1 surface AND a curved 3-D GRIN field. Without bidirectional/photon-beam estimators or manifold specular-chain sampling, backward MC has near-zero hit rate on the sun and forward MC has near-zero hit rate on a pinhole/orthographic pixel — rendering the "reference" unconverged or biased if straight shadow rays are substituted.

Failure condition: The MC reference either drops GRIN curvature during null-collision tracking, uses unrefracted/straight shadow rays for path radiance NEE (systematic bias >10%), or fails to achieve <2% Monte Carlo relative standard error on the separated veil image and halo tail.

C7. Single scattering is adequate at cL ≈ 0.2, multiple scattering becomes first-order at cL → 1+, decided by order-capped MC

Verdict: WEAKENED

Counterargument: Using controlled order-capped MC comparison is the right methodology, and multiple scattering is certainly first-order by cL ∼ 1. However, the test design misses three physical controls:

cL is the wrong control variable: Multiple scattering depends on scattering optical depth bL = ω₀cL, not cL. Pure CDOM-dominated water (ω₀ = 0.15) at cL = 2 has bL = 0.3 (weak multiple scattering), whereas clay/calcite turbid water (ω₀ = 0.9) at cL = 0.8 and solar zenith 45° has two-way scattering depth bL 2-way ≈ 1.6, where multiple scattering already dominates.

The 4.0 m sub-verification scene already exceeds cL = 1 even in fairly clear water (c = 0.3 m⁻¹ ⟹ cL = 1.2).

Observable blindness: Measuring only deposit band energies and halo containment misses path radiance (veil). Because FF backscatter (θs > 110°) is very weak, two-step scattering (one small-angle forward scatter + one moderate scatter, or bed-reflected diffuse scatter) becomes first-order in the veil image at much lower bL (∼0.3–0.5) than in the bed deposit.

Failure condition: Sweeping ω₀ ∈ [0.2, 0.9] at fixed cL = 1.0 changes multiple-scattering error from <3% to >40%, or order-1 vs uncapped MC disagrees by >25% in the separated veil image at cL = 0.5 (ω₀ = 0.85) while deposit metrics still agree.

C8. Halocline as a sharp step vs thermal plume as continuous GRIN

Verdict: WEAKENED

Counterargument: The second half stands: staircasing a continuous thermal plume crushes ∇⊥T or introduces faceted voxel-boundary artifacts. However, treating a dynamic brackish river-plume front / salt wedge as a discontinuous step ignores shear-layer physics and creates an optical trap. Even though salt diffusivity is 100× smaller than heat (Le ≈ 100), a turbulent sheared front with ΔU ≈ 0.2 m/s and ΔS = 10 psu (g′ ≈ 0.075 m/s²) maintains an interfacial mixing thickness limited by critical gradient Richardson number (Rig ≈ 0.25): hi ∼ 0.25ΔU²/g′ ≈ 0.13 m (≫ λ). While Snell's invariant makes exit angle insensitive to hi on a flat horizontal interface, at a steeply curved river-plume front a discontinuous step triggers spurious Fresnel reflection and total internal reflection (TIR) (since Δn ≈ 0.0028 gives TIR within ∼3.7° of tangent) plus sharp caustic discontinuities, whereas a continuous 3–15 cm transition turns rays smoothly with zero Fresnel reflection. Furthermore, 3.3(iii) runs a freshwater thermal overflow through this same sharp-step lane — directly contradicting C8's physical distinction between salt and temperature interfaces.

Failure condition: At a tilted river-plume front (30°–75° slope), comparing a discontinuous step (with Fresnel/TIR) against a continuous GRIN interface of realistic shear thickness hi = 0.25ΔU²/g′ produces spurious bright TIR/reflection bands or torn caustic discontinuities absent in the continuous transition.

C9. Ellison & Turner E(Ri) transfers to 0.6 m creeks without modification, like Dolcetti's data in Stage 1

Verdict: WEAKENED

Counterargument: Depth-wise, 0.6 m creeks (Re ∼10⁴) sit near laboratory flume scale, so Re is not the main flaw. The analogy to Dolcetti's data fails because Dolcetti measured the exact target system (rough-bed open-channel surface waves), whereas Ellison & Turner (1959) / Hallberg curves describe steady, uniform, wall-bounded dense overflows or unbounded jets — missing three first-order mechanisms in a 0.6 m brackish river plume:

No bottom-wall turbulence at a surface plume interface, but active Stage-1 wind/wave mixing: Surface plumes lack bottom drag at the pycnocline, while wind shear and breaking/orbital waves drive entrainment even at Ri > 0.8 where Ellison–Turner shuts off (E → 0).

Front-head dynamics vs steady body: 3.3 prioritizes river-plume fronts. Front speed and intense head mixing/downwelling are governed by Benjamin (1968) / Simpson front shock dynamics and lobe-and-cleft instability, not steady body E(Ri).

Finite depth ratio h₁/H ∼ 0.3–0.7: Verification 3.3(i) cites front speed ∼ √(g′h₁), which is the deep-ambient limit (h₁/H → 0). In a 0.6 m creek where h₁ ∼ 0.2–0.4 m, lower-layer confinement and bed drag slow the front to ∼ √(g′h₁h₂/H) — up to 40–50% slower.

Failure condition: Tested against shallow-creek/flume river-plume fronts at h₁/H ∈ [0.3, 0.7], √(g′h₁) overpredicts front speed by >30%, or wind-ruffled plumes continue entraining at Ri > 1 where the unmodified E(Ri) predicts zero mixing.

C10. Reduced models producing (n, a, b) plus exact transport can satisfy "reproducing nature" without a 3-D fluid PDE solver

Verdict: WEAKENED

Counterargument: Two-layer shallow-water equations are a physically grounded asymptotic PDE for H/L ≪ 1, so C10 holds for depth-integrated macro-quantities (front propagation, total optical depth τ, bulk layer thickness). It fails wherever exact light transport is sensitive to spatial derivatives or sub-layer vertical profiles that the reduced model discards: (1) ∇n and ∇²n in thermal schlieren (where reduced plume models invent the gradient width), (2) vertical b(z) profiles inside a layer (which control the halo lever arm (L − s), per C4), and (3) shared advection with Stage 1 (wave orbital velocities in a 0.6 m pond penetrate to the bed, oscillating subsurface turbidity/plumes in phase with surface ripples). Running an "exact" ray tracer on ungrounded derivatives does not reproduce nature; it reproduces the reduced model's shape assumptions.

Failure condition: Two reduced-model parametrisations conserving identical bulk heat/sediment budgets yield >2× discrepancy in GRIN deflection variance or halo half-energy radius because they assume different unmodeled vertical profiles or gradient thicknesses.

What This Design Misses Entirely

1. 3-D Volumetric Caustic Irradiance Field E(x,t) (Required for Light Shafts) — Block 3.2 & Stage 1 → 3 Coupling

Why it matters physically: The current stack projects caustics only onto the 2-D bed mesh. Without tracking/depositing focused irradiance inside the 3-D water column (e.g. photon beams / ray-marching against surface-refracted sun rays), H4's "visible light shafts in the water" cannot exist — path radiance reduces to a flat, structureless haze.

2. Pure-Water Molecular Scattering bw(λ) vs Particulate Scattering bp(λ) — Block 3.2

Why it matters physically: Natural water scattering is b(λ) = bw(λ) + bp(λ). Pure-water molecular scattering (∝ λ⁻⁴·³²) has a Rayleigh-like phase function (∼50% backscatter), whereas particulates follow Fournier–Forand (≪1% backscatter). In clear ponds (cL ≈ 0.2 or the 4.0 m clear sub-scene), bw can provide 20–60% of total backscattering at blue/green wavelengths. Applying the FF phase function to total b severely underpredicts the upward path radiance (veil).

3. Inelastic Scattering (Raman Scattering by H₂O and CDOM/Chlorophyll Fluorescence) — Block 3.2

Why it matters physically: In clear to moderately CDOM-rich waters (especially at 4.0 m depth), Raman scattering and CDOM fluorescence shift abundant short-wavelength solar photons into 500–650 nm where pure-water absorption drops, contributing 10–30% of upwelling path radiance in natural waters. Neither 3.2 nor the MC reference includes cross-band energy transfer.

4. Internal Total Internal Reflection (TIR) at the Wavy Stage-1 Surface & Fresnel Coupling at Haloclines — Blocks 3.2, 3.3 & Stage 1 ↔ 3 Coupling

Why it matters physically: Upwelling light from the bed and in-water scattering striking the wavy air-water interface beyond θc ≈ 48.6° undergoes 100% TIR back into the volume, raising bed illumination and multiple-scattered veil in turbid/shallow scenes by 10–30%. Additionally, view rays from an above-water camera must refract through the wavy Stage-1 surface before traversing Stage 3. Neither surface TIR feedback nor view-side surface refraction is accounted for in 3.2/3.4.

5. Polarization in RTE (Vector vs Scalar Transport) — Block 3.2

Why it matters physically: Skylight and oblique Fresnel-transmitted sunlight are strongly polarized, and Rayleigh/small-particle water scattering is polarization-sensitive. Scalar RTE introduces systematic 5–15% errors in upwelling veil radiance — growing larger once Stage 4 introduces tilted camera angles.

6. Chromatic Dispersion (n(λ)) in Ray Integration & Interpolation Regularity — Block 3.1

Why it matters physically: Water's chromatic dispersion across 400–700 nm is Δnλ ≈ 0.013 — roughly 4× larger than a ΔS = 20 psu halocline and 40× larger than a 3 K thermal plume. The draft cites Quan & Fry n(S,T,λ) but does not state whether GRIN rays are traced separately per wavelength band (and whether C¹-smooth interpolation is used on the n(x) grid; standard trilinear interpolation creates discontinuous ∇n across voxel faces, generating artificial caustic streaks).

7. Stage-1 Wave-Orbital Advection & Stage-2 Bathymetry/Wakes — Blocks 3.1, 3.3, 3.4 & Stage 1/2 ↔ 3 Coupling

Why it matters physically: Surface wave orbital motion penetrates to depth ∼ λwave/2 (spanning the entire 0.6 m window for ordinary wind chop), oscillating subsurface plumes, generating internal waves on the 3.3 halocline, and driving wave-turbulence entrainment. Likewise, 3.3 treats two-layer shallow water without Stage-2 bed bathymetry (which steers salt wedges and pools dense water) or bed roughness / submerged object wakes (which drive bottom drag and local sediment resuspension b(x)).

8. Distinct Optical Model for Stage-2 Bubbles — Block 3.2 & Stage 2 → 3 Interface

Why it matters physically: Stage 3 claims to "consume a concentration field if Stage 2 supplies one" for bubbles, but air-in-water microbubbles (nrel = 0.75 < 1, critical-angle scattering, ω₀ ≈ 1, high backscatter) have a fundamentally different phase function from hydrosols and cannot simply be summed into 3.2's FF particulate b(x).

9. Critical Validation Gaps — Blocks 3.1–3.4

- Periodic loop-closure vs non-periodic fronts (3.4): Claiming loop closure via "periodic boundaries" does not bend field physics is false for 3.3 river-plume fronts and salt wedges, which are non-periodic inflow-to-ambient propagating shocks that will collide with their own periodic tail inside an 8-second loop.
- Zero empirical external validation: Unlike the settled homogeneous stack (validated against published phase-function measurements), Stage 3 only checks analytic toy limits and self-consistency against its own MC solver — never testing against published real-water BOS/schlieren deflection measurements, lab-tank inhomogeneous MTF/halo data, or shallow estuarine front mixing profiles.
- Missing multi-component energy-closure audit: No verification check ensures that spectral solar input equals water-absorbed + bed-absorbed + upwelling-transmitted energy across the composite (deposit + halo + veil).

Final Verdict on the Claim Under Test

FALLS — The four-block decomposition cannot stand as currently shaped and must change in five concrete ways:

1. Replace Block 3.2's homogeneous (τ, b/c) FF-CDF halo evaluation with depth-weighted beam-spread / lever-arm (L − s) integration, split scattering into Rayleigh pure-water bw(λ) + particulate FF bp(λ) (plus bubble phase function if supplied by Stage 2), and add a 3-D volumetric caustic irradiance field E(x,t) so path radiance actually reproduces light shafts.
2. Re-parameterise Block 3.2's multiple-scattering sweep and MC reference around scattering optical depth bL = ω₀cL and separated veil accuracy, and specify curved-ray RK4-Woodcock stepping plus wavy-surface/GRIN photon-beam or bidirectional sampling so the MC reference does not fail on NEE.
3. Ground Block 3.1's thermal plume gradients in physical scaling (Morton–Taylor–Turner buoyancy/entrainment and boundary-layer thickness, plus C¹ grid interpolation and step-size convergence) rather than unpinned shape assumptions.
4. Fix Block 3.3's halocline and front physics: replace the discontinuous Fresnel/TIR step at tilted fronts with a shear-limited interfacial thickness hi ∼ 0.25ΔU²/g′, include finite-depth-ratio h₁/H front speed, Stage-2 bathymetry, and Stage-1 wave/wind entrainment above Ri > 0.8.
5. Fix Block 3.4's validation and loop design: remove the contradiction of periodic boundary conditions for propagating river-plume fronts, and add empirical external benchmarks (lab/field schlieren and inhomogeneous water optics) alongside internal MC checks.
