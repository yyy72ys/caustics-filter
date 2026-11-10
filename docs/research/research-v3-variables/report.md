# v3 Variable Research — Mathematical Grounding for New Nodes

Chunk 0 research for the caustics/water filter v3. Purpose: every proposed new variable gets (i) a definition, (ii) the established formula, (iii) a real-world range, (iv) its coupling to variables already in the pipeline (wind speed/direction, wave fineness, depth, transparency, turbidity, pigment, sun position), and (v) an implementation note naming the stage it plugs into — **state simulation** (the Phillips-FFT water surface), **light transport** (deterministic Jacobian density + sun-disc quadrature), **composite** (the frozen E_bottom = E_amb + E_base,dir × C form), or a new **observer post-stage**.

Conventions: formulas are stated as published; numbers I evaluated myself from published formulas are marked *(computed)*. URLs are given inline per claim. Nothing here is an arbitrary slider: where a quantity has no exact definition, §7 says so and names the closest defensible proxy.

---

## 1. Water surface

### 1.1 Significant wave height H_s vs wind speed (Pierson–Moskowitz)

**Definition.** Significant wave height is defined spectrally as H_m0 = 4√m₀, where m₀ = ∫ S(ω) dω is the total variance of surface elevation (the modern form of the older H₁⁄₃, the mean height of the highest third of waves). (WMO, *Guide to Wave Analysis and Forecasting*, ch. 1: https://www.yumpu.com/en/document/view/47477002/guide-wave-analysis-and-forecasting-wmo/25)

**Formula.** The Pierson–Moskowitz (1964) spectrum for a *fully developed sea* (unlimited fetch and duration) is a one-parameter family:

  S(ω) = α g² ω⁻⁵ · exp[ −β (ω₀/ω)⁴ ],  α = 0.0081, β = 0.74, ω₀ = g / U₁₉.₅

with peak frequency f_p = 0.877 g / (2π U). Integrating gives the closed-form wind relation

  **H_s = 0.0246 · U₁₀²** (H_s in m, U₁₀ = wind speed at 10 m in m/s),

equivalently H₁⁄₃ = 0.21 U₁₉.₅² / g ≈ 0.22 U₁₀² / g. (WMO guide, eqs. 1.28–1.31, same URL as above; Stewart, *Introduction to Physical Oceanography*, ch. 16, eqs. 16.30–16.35: https://nnty.fun/downloads/books/cdn.preterhuman.net/texts/science_and_technology/physics/Introduction_to_Physical_Oceanography/Chap16.pdf) The PM spectrum is single-parameter: specifying any one of U, ω_p, or H_s fixes the whole spectrum (spectrum notes: https://api.bitbucket.org/2.0/repositories/brian_bingham/ocean_notes/patch/0f477d1ac23a9798e6917189c5fecea8a336167b).

**Real-world range** *(computed from H_s = 0.0246 U₁₀²)*:

| U₁₀ (m/s) | 2 | 5 | 8 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|
| H_s (m), fully developed | 0.10 | 0.62 | 1.57 | 2.46 | 5.5 | 9.8 |

These are open-ocean values. A pond, river reach, or the ~1.5 m patch this pipeline simulates is **never** fully developed: fetch and duration cap the growth long before PM equilibrium. PM therefore must not be applied raw at scene scale (see Honesty notes).

**Coupling.** Wind speed already drives the pipeline's Phillips spectrum amplitude and peak placement. H_s adds nothing new physically if it is *derived* from wind via PM — its real use is the reverse direction: H_s is the quantity a user can set from observation (a photo of a lake), and then m₀ = (H_s/4)² fixes the spectrum's total variance directly.

**Implementation note (state simulation).** Two defensible modes:
- *Wind-driven (existing path):* keep the Phillips spectrum as is; *report* the implied H_s by numerically integrating the realised spectrum, H_s = 4√m₀, and display it as a read-only derived node value.
- *Height-driven (new input):* rescale the spectrum so its integral equals m₀ = (H_s/4)² while keeping the peak wavenumber where wind/fetch put it. This is exactly the two-parameter (Bretschneider/modified-PM) idea — amplitude and peak frequency specified independently — and it stays inside the FFT framework: one scalar multiplier on the spectrum before the transform. A physical cap applies: individual-wave steepness cannot exceed the Stokes limiting steepness H/λ ≈ 1/7 ≈ 0.14 before breaking (classical result, cf. Stewart ch. 16 above); in-patch, cap the realised crest slope accordingly rather than letting H_s grow without bound on a 1.5 m patch.

### 1.2 "Calmness" — the derivable quantity is mean square slope (Cox–Munk)

**Definition.** There is no physical unit of "calmness". The measured quantity that behaves like it — and that directly controls sun glitter and caustic texture — is the **mean square slope** (mss, σ²) of the surface, the variance of the surface gradient, estimated by Cox & Munk (1954) from aerial photographs of sun glitter. The slope distribution is approximately bivariate Gaussian in along-wind (a) and cross-wind (c) components.

**Formula** (U = wind speed in m/s at 12.5 m):

  σ_a² = 3.16×10⁻³ · U
  σ_c² = 0.003 + 1.92×10⁻³ · U
  **σ² = σ_a² + σ_c² = 0.003 + 5.12×10⁻³ · U**   (± 0.004)

(Cox & Munk 1954, as tabulated in Mobley, *Ocean Optics Web Book*, "Cox–Munk Sea Surface Slope Statistics": https://www.oceanopticsbook.info/view/surfaces/cox-munk-sea-surface-slope-statistics)

**Real-world range** *(computed)*: σ² goes from 0.003 at U = 0 (rms slope ≈ 3.1° — a residual not caused by local wind) to ≈ 0.029 at U = 5 (rms ≈ 9.7°), ≈ 0.054 at U = 10 (rms ≈ 13.5°), ≈ 0.075 at U = 14 (rms ≈ 15.9°; Cox & Munk's own summary: mss reaches (tan 16°)² at 14 m/s — Cox & Munk 1954, JOSA 44:838, abstract: https://misclab.umeoce.maine.edu/ftp/classes/OO2017/readings/Cox_Munk_JOSA1954.pdf). Two honest modifiers from the same literature: natural surface slicks (surface films) strongly suppress mss at given wind — this is the physical basis of "oily calm" — and an air–sea temperature difference shifts roughness at fixed wind (water warmer than air → rougher than Cox–Munk predicts) (NOAA PSL, "Glittering Light on Water": https://psl.noaa.gov/outreach/education/science/glitter/).

**Coupling.** σ² is *the* coupling hub for v3:
- It is an integral of the pipeline's own spectrum: σ² = ∫ k² S(k) dk over the simulated band. So "calmness" needs no new input at all — it is a readout of the existing state, and it is also the natural target the state should hit when wind is the input (calibrate the spectrum's high-k tail so realised σ² matches Cox–Munk at the chosen U, over the resolved band).
- It drives glitter width in the composite (glitter is literally the image of the slope pdf).
- It drives wave shadowing (§1.3) and the rain-ripple visibility (§3: rain rings are only visible because the background σ² is near zero).

**Implementation note.**
- *State simulation:* compute realised σ² per frame from the FFT slopes (free — the transport already needs gradients) and expose it as the node **"calmness" = σ²** (label it with units rad² and the Cox–Munk reference value at the current wind). If a user-facing 0–1 calmness display is wanted, define it once as C = σ²_ref(U)/σ² or a fixed mapping against σ² ∈ [0.003, 0.075], and document it as a display transform, not physics.
- *Composite:* the sun-glitter term should sample the *anisotropic* (σ_a², σ_c²) Gaussian, which automatically narrows the glitter path as calmness rises. The current isotropic microfacet glitter is the σ_a = σ_c special case.

### 1.3 Wave shadowing — Smith/Beckmann masking-shadowing, and shadows on the bottom

**Definition.** On a rough surface, some facets are hidden from the light (shadowing) or from the viewer (masking) by neighbouring surface. Smith (1967) derived the visible fraction for a random rough surface from its slope distribution alone. In modern microfacet form, the one-direction function is

  **G₁(ω) = 1 / (1 + Λ(ω))**,

where Λ(ω) is the mean invisible (masked) facet area per unit visible area, a functional of the slope pdf (derivation and definitions: GAMES-101 lecture notes, "Microfacet Material", appendix: https://iskxcr.github.io/CG-Note-Misc/GAMES-101/readable-html/Lecture17.html; Smith 1967: B. G. Smith, "Geometrical shadowing of a random rough surface", IEEE Trans. Antennas Propag. 15(5):668–671, cited in https://yusuketokuyoshi.com/papers/2019/ImprovedGeometricSpecularAA(supplemental).pdf).

**Formulas.** For a Beckmann (Gaussian-slope) surface, Λ has an exact form in terms of the complementary error function; the closed-form rational approximation in standard use (Walter et al. 2007) is, with a = 1/(α tan θ), α = rms slope:

  Λ(a) ≈ (1 − 1.259a + 0.396a²) / (3.535a + 2.181a²) for a < 1.6, and Λ = 0 for a ≥ 1.6.

For a GGX slope distribution the Smith term has the simple closed form

  G₁(θ) = 2 cos θ / ( cos θ + √( α² + (1 − α²) cos²θ ) )

(Smith-GGX implementation reference: https://deepwiki.com/hunterloftis/pbr/3.7-bsdf). The two-direction term should use the height-correlated form G₂ = χ⁺/(1 + Λ(ω_i) + Λ(ω_o)) (Heitz 2014, as given in the Square Enix technical report: http://www.jp.square-enix.com/tech/library/pdf/Error%20Reduction%20and%20Simplification%20for%20Shading%20Anti-Aliasing.pdf).

**Where it applies in this pipeline — two different places, do not conflate them:**
1. **Surface reflection (glitter) at grazing view/sun.** The specular term must be multiplied by G₂(sun, view) with α set by the *resolved* slope statistics (Cox–Munk σ_a, σ_c for anisotropy). Without it, grazing glitter is systematically too bright. Plugs into: **composite** (glitter term) — cheap, analytic.
2. **Shadowing of the bottom ("wave shadows").** Physically this is not a BRDF term: it is the fraction of the solar disc's rays, refracted at the surface, that never reach a bottom point because the ray's entry facet is oriented away or the straight in-water path is interrupted by another wave. Two hard physical bounds keep it honest in this scene: refraction compresses all in-water sun directions into a cone of half-angle ≤ 48.6° from vertical (the critical angle, sin⁻¹(1/1.333); Mobley, Fresnel equations for polarization: https://www.oceanopticsbook.info/view/surfaces/level-2/fresnel-equations-for-polarization/pdf), and wave slopes are modest (rms ≲ 16° even at 14 m/s wind, §1.2). True bottom shadowing is therefore a *low-sun, steep-short-wave* effect, not an all-day one.

**Implementation note (light transport).** The pipeline's transport already shoots a sun-disc angular quadrature per bottom sample. The exact treatment costs nothing structurally: for each quadrature ray, test whether the segment from the surface entry point, traced back along the incident direction, re-intersects the heightfield within the patch (a 1-D march along the ray's surface footprint); multiply that sample's weight by V ∈ {0,1} (soft: by the unoccluded fraction of the sample's footprint). The cheap statistical alternative is a single multiplier G₁(sun) with α = realised rms slope from §1.2, applied to the direct term in the **composite** — acceptable as a v3-first step, labelled as the statistical (Smith) approximation, not ray-exact. Do not apply both (double counting).

---

## 2. Bottom: surface reflectance vs bottom albedo, and depth

### 2.1 Two different "reflectances" — keep them as separate nodes

**Surface reflectance (Fresnel).** This is a property of the air–water interface alone, fixed by the refractive index (n ≈ 1.333, with the small per-channel dispersion the pipeline already uses). For unpolarised light the Fresnel reflectance at incidence angle θᵢ is R(θᵢ) = ½(R_s + R_p) with the standard amplitude coefficients; at normal incidence it reduces to

  **R₀ = ((n₂ − n₁)/(n₂ + n₁))² = ((1.333 − 1)/(1.333 + 1))² ≈ 0.020 (≈ 2 %)**

(worked example with R + T = 1 check: https://www.physics2u.com/vault/fresnel-equations-boundary.html; Mobley's Fresnel treatment gives 0.021 for n = 1.34 and the 48.6° critical angle for the reverse direction: https://www.oceanopticsbook.info/view/surfaces/level-2/fresnel-equations-for-polarization/pdf). R rises slowly at first and then steeply: ≈ 2.0 % at 0°, 2.1 % at 30°, 2.8 % at 45°, 6.0 % at 60°, 21 % at 75°, → 100 % at grazing (angle table: https://www.illuminationpros.net/articles/calculating_the_refraction_index_of_light_penetrating_moving_water/). Per-channel with the pipeline's IOR set (1.33132 / 1.33434 / 1.33925), R₀ = 0.0202 / 0.0205 / 0.0210 *(computed)* — a real but tiny dispersion effect.

Consequence for the node design: **"surface reflectance" is not a free variable.** The only honest inputs are n(λ) (fixed for water) and the angle (already in the sun/view geometry). What a slider could legitimately control is nothing — the Fresnel term should stay computed. Its role in the energy split is large at grazing (the surface mirror wins) and small near noon (98 % of the light enters the water and can reach the bottom).

**Bottom reflectance (albedo).** The bottom is (to first order) a Lambertian reflector with spectral albedo ρ_b(λ): irradiance arriving at the bottom E(0⁻, λ) leaves as radiance L = ρ_b(λ)/π · E. This *is* a genuine material variable, set by substrate type. Retrieved broadband bottom albedos from shallow-water remote sensing (Roatán, IKONOS; Mishra et al. 2005, IEEE TGRS 43(7):1592): **seagrass ≈ 0.12, coral ≈ 0.12–0.24, sand ≈ 0.24**, summarised by the authors as "from very dark (coral ≈ 0.12) to very bright (sand ≈ 0.3)" (https://core.ac.uk/download/pdf/17238883.pdf). Spectrally resolved in-situ values spread wider: Posidonia seagrass stays **below 0.10 across 500–600 nm**, while carbonate sand peaks as high as ≈ 0.63 at its brightest wavelengths in spectroradiometer measurements (Crete study: https://hal.science/hal-02425707v1/file/Lemenkova-Rostock.pdf) — peak spectral values and broadband albedo are different quantities; do not mix them in one node. Fogarty et al. (2018, JGR Oceans 123:53–73) model exactly the three cases this filter cares about — bright sand, seagrass canopy, turbid water — and find the *combined* surface albedo over shallow bright bottoms reaching ≈ 0.20 vs the 0.05 open-ocean parameterisation (abstract: https://www.mendeley.com/catalogue/887da399-a21b-3d06-a80f-10cb96f7d849/).

**Real-world range for the node:** ρ_b ≈ 0.05 (dark mud/organic sediment — darker than seagrass; no single canonical table, see Honesty notes) to ≈ 0.30 (bright carbonate sand, broadband); spectral shape matters more than the scalar: sand rises toward red, seagrass has the vegetation green peak and red absorption, rock is spectrally flatter and typically 0.1–0.25 dry, less when wet/submerged.

**Coupling.** ρ_b multiplies only the direct-bottom term (E_base,dir × C) in the composite — the ambient term already arrived via a different path. It couples to depth and turbidity through §2.2: what the camera receives is ρ_b × (light that survived the trip down) × (survival of the trip back up). In the web app the "bottom" is the user's photo itself, so ρ_b maps to a *spectral gain on the photo's linear values before the caustic multiplier*, not to a colour tint after it.

**Implementation note (composite).** Add ρ_b as three per-channel gains (or a small substrate preset table: sand / rock / vegetation / mud with cited spectra) applied to E_base,dir inside the existing separated-pass composite. The surface Fresnel needs no new node — it is already implied by sun angle and IOR in transport; if the reflection/glitter view is ever rendered, use R(θ) with the Smith G₂ of §1.3.

### 2.2 Bottom range (depth): Beer–Lambert with real absorption spectra

**Definition/formula.** Along any straight in-water path of length s, radiance is attenuated exponentially:

  **L(s, λ) = L(0, λ) · exp( −c(λ) · s ),  c(λ) = a(λ) + b(λ)**

(Beer–Lambert / beam attenuation; a = absorption, b = scattering). For the direct (unscattered) sun beam that makes caustics, this is exact; scattering removes light from the beam (and returns some as diffuse ambient, which is what the pipeline's E_amb term represents).

**Pure-water absorption a_w(λ)** — Pope & Fry (1997), Appl. Opt. 36:8710–8723, the accepted reference values (as tabulated in the NASA Ocean Optics Protocols, Table 1.1: https://oceancolor.gsfc.nasa.gov/files/resources/docs/technical/volivch1err1.pdf; same values in Kirk's *Light and Photosynthesis in Aquatic Ecosystems*, Table 3.1: https://1library.co/article/recursos-personales-dise%C3%B1o-programaci%C3%B3n-geograf%C3%ADa-escolar-partir-adquisici.lzg01w2q):

| λ (nm) | 400 | 420 | 440 | 500 | 550 | 600 | 650 | 700 |
|---|---|---|---|---|---|---|---|
| a_w (m⁻¹) | 0.0066 | 0.0045 | 0.0064 | 0.0204 | 0.0565 | 0.2224 | 0.340 | 0.624 |

The minimum is in the blue (~420–440 nm, a ≈ 0.0045/m → 1 % loss per 2.2 m); the red end is brutal: at 650 nm a 1 m path already transmits only e^(−0.34) ≈ 0.71, and at 700 nm only ≈ 0.54. This asymmetry *is* the physical origin of the water-colour tint. Pure-water scattering b_w is small (order 0.001–0.01 /m, Rayleigh-like, ∝ λ⁻⁴·³²; same NASA table) — in clear water, beam attenuation is absorption-dominated; the pipeline's **turbidity** node is precisely the addition Δa(λ) + Δb(λ) on top of these pure-water floors (particulate scattering is spectrally much flatter than a_w, which is why turbid water goes green-grey rather than deep blue).

**Path geometry (the "bottom range" coupling).** Depth z enters twice with different path lengths: down-path s↓ = z / cos θ_w (θ_w = in-water sun angle, ≤ 48.6° after refraction) for the direct term, and — for a viewer above the surface — an up-path of the same form for the reflected bottom radiance. So the bottom term carries exp(−c(λ)·z·(1/cos θ_w↓ + 1/cos θ_w↑)) ≈ exp(−2 c(λ) z) near vertical sun. Worked example *(computed)*: clear water, z = 2 m, vertical sun, two-way: at 550 nm, exp(−2·0.0565·2) ≈ 0.80; at 650 nm, exp(−2·0.34·2) ≈ 0.26. Red dies with depth even in perfectly clear water — a depth node that only dims uniformly is physically wrong; it must be spectral.

**Implementation note (light transport + composite).** The transport already produces C as a *ratio* against the flat-surface baseline, which correctly cancels the common attenuation; the spectral Beer–Lambert factor belongs on the *absolute* terms in the composite: E_base,dir(λ) ← E_base,dir(λ) · exp(−c(λ) s↓) and the returned bottom light × exp(−c(λ) s↑). Implementation is three exponentials per frame (RGB), with a(λ) interpolated from the Pope & Fry table at the pipeline's channel wavelengths (use the R/G/B effective wavelengths implied by its dispersion set, ≈ 650/550/470 nm as an approximation — and say so). Depth z and turbidity are already nodes; this section just replaces any ad-hoc "depth darkening" with c(λ).

### 2.3 How caustic contrast changes with depth

There is **no universal contrast-vs-depth law**, and the report does not invent one. What the literature supports:

- The caustic pattern is a flux-density map: total flux is conserved (mean of C stays 1 in the pipeline's normalisation), so depth redistributes contrast rather than destroying it — bright filaments sharpen where rays focus and the pattern's *spatial scale grows with depth* (the lensing geometry magnifies).
- Empirically, the pattern's spatial frequency is highest at a depth of about **five times the waves' crest-to-crest distance**, and the pattern **lessens with increasing depth, turbidity, and diffuseness of the incident light** (Venables et al. 2022, Sci. Adv. 8:eabq2770, reviewing the underwater light-field literature: https://www.researchgate.net/publication/363424652_Polarization_vision_mitigates_visual_noise_from_flickering_light_underwater).
- The depth-decay of *contrast* in real water is dominated by scattering, not geometry: each scattering event randomises ray directions, low-passing the pattern. Tessendorf's ocean notes model the underwater light field with a diffuse attenuation coefficient K in the range 0.03–0.1 /m for ocean water and still show structured caustics at 100 m in clear-water simulation — i.e., in the zero-scattering limit the pattern does not vanish with depth (https://3dmodelizm.ru/ext/graphics.ucsd.edu/courses/rendering/2005/jdewall/tessendorf.pdf).

**Defensible model for the pipeline (implementation note, light transport).** Keep the geometric transport contrast as computed (it already has the correct z-dependence through ray geometry), and let depth act on contrast only through the two physical channels already derived: (i) spectral Beer–Lambert on the direct term (§2.2), and (ii) a turbidity-driven scattering blur: convolve C with a kernel whose width grows with the scattering optical depth b·z (e.g., σ_blur ∝ √(b z) · (pattern scale)), representing forward-scatter diffusion. Both knobs are existing nodes (depth, turbidity) — no new "contrast" slider is needed or honest.

---

## 3. Rain on windless (calm) water — the new scene

Scene inputs: rain rate R (mm/h) is the one new physical input; wind = 0 (or near it), so the background surface is the σ² ≈ 0.003 Cox–Munk floor of §1.2. Rain rate categories for reference: light 0.2–2, moderate 2–10, heavy 10–30, very heavy 30–50, intense 50–100 mm/h (radar-based categories tabulated in https://refubium.fu-berlin.de/bitstream/handle/fub188/5927/02_andi_walther_chapter2.pdf?sequence=3&isAllowed=y).

### 3.1 Drop size distribution — Marshall–Palmer

**Definition/formula.** The number concentration of drops of diameter D (mm) per unit volume of air:

  **N(D) = N₀ · exp(−ΛD),  N₀ = 8000 m⁻³ mm⁻¹,  Λ = 4.1 · R^(−0.21) mm⁻¹** (R in mm/h)

(Marshall & Palmer 1948, for stratiform rain; https://en.wikipedia.org/wiki/Raindrop_size_distribution, and the same formula with the Z–R context in the FU Berlin chapter above). Consequences *(computed)*: mean diameter 1/Λ rises from 0.24 mm at R = 1 to 0.40 mm at R = 10 and 0.56 mm at R = 50; the mass-median diameter D₀ = 3.67/Λ is ≈ 1.45 mm at R = 10. Real drops span ≈ 0.3–6 mm; drops larger than ~6 mm break up in fall (Gunn & Kinzer data discussion: https://www.madsci.org/posts/archives/2000-07/962626446.Ph.r.html).

**Caveat (feeds Honesty notes).** Marshall–Palmer is an *average over many stratiform events* with a fixed N₀; measured DSDs vary event-to-event and are better fit by normalised gamma distributions. For a renderer's statistics (ring sizes, impact rate) it is the defensible default; it is not a law of nature for any given minute of rain (review discussion of exactly this limitation: https://amt.copernicus.org/preprints/amt-2020-174/amt-2020-174-AR1.pdf).

### 3.2 Fall speed — terminal velocity

**Formula.** Atlas & Ulbrich (1977), fit to the Gunn & Kinzer (1949) measurements (sea level, still air), D in mm, v in m/s:

  **v(D) = 9.65 − 10.3 · exp(−0.6 D)**, valid for D ≈ 0.5–5 mm

(the formula and its validity band 0.05–0.50 cm are stated in https://weather.missouri.edu/rains/cooleyresearch.pdf). Values *(computed)*: v(1 mm) ≈ 4.0 m/s, v(2) ≈ 6.5, v(3) ≈ 8.0, v(5) ≈ 9.1 m/s — matching the Gunn & Kinzer table (1.2 mm → 4.64 m/s; 4.6 mm → 9.03 m/s; same madsci.org source as above). Impact kinetic energy per drop, ½ m v² with m = (π/6)ρD³, is what sets ring-wave amplitude scaling in §3.4. Altitude/air-density corrections exist (thinner air aloft → faster fall) but are irrelevant at a pond surface.

### 3.3 Impact frequency per unit area

**Definition/formula.** The number of impacts per m² per second is the flux of the DSD through a horizontal plane:

  **F = ∫₀^∞ N(D) · v(D) dD**  [drops m⁻² s⁻¹]

(consistency: the same integral weighted by drop volume, (π/6)∫N(D)D³v(D)dD, returns R by construction). Evaluating the MP + Atlas–Ulbrich pair *(computed, integral taken over the validity band and slightly beyond; treat as ±factor ~2 because of the unmodelled sub-0.5 mm population)*:

| R (mm/h) | 1 | 10 | 25 | 50 |
|---|---|---|---|---|
| F (impacts m⁻² s⁻¹) | ≈ 1.3×10³ | ≈ 4.2×10³ | ≈ 6.4×10³ | ≈ 8.6×10³ |

Note the honest shape of this: impact *count* grows slowly with R (most extra rain volume arrives as bigger drops, not more drops), while mean drop size and per-impact energy grow steadily. In the pipeline's ≈ 1.5 m × 1.1 m patch at R = 10 mm/h that is ≈ 7,000 impacts per second over the patch — far too many to simulate individually, which dictates the statistical treatment in §3.5.

### 3.4 Ring waves from a single impact — wavelengths, lifetime, decay

**Phenomenology and timescales** (high-speed imaging + theory, reviewed in Chang et al. 2016, J. Geophys. Res. Oceans, "A study of radar backscattering from water surface in response to rainfall": https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JC010975?cookieSet=1): an impact makes a crater/cavity and a crown (growth + decay ≈ 32 ms), then a central stalk (Worthington jet, up to ≈ 110 ms; stalk heights ~2–3 cm), and radiating **ring waves**; the whole impact process lasts ≈ 150 ms.

**Ring-wave wavelengths.** The ring waves are gravity–capillary waves with **wavelengths of 1–2 cm** (theory: Le Méhauté 1988; experiments by Bliven and others), remaining visible only within a circular region of **diameter ≈ 25–50 cm** around the impact (Wetzel 1990) — both as reviewed in the JGR paper above; the TOPEX rain study gives the same numbers (rings expanding to 25–50 cm diameter before disappearing): https://www.academia.edu/28145116/The_effects_of_rain_on_TOPEX_Poseidon_altimeter_data. The spectral peak frequency of the ring waves decreases as drop diameter increases and is nearly independent of rain rate (JGR paper above) — i.e., *bigger drops → longer rings*, which is the coupling that makes R affect ring scale via Λ(R), not via impact count.

**Dispersion (why 1–2 cm).** Gravity–capillary waves obey ω² = gk + (σ/ρ)k³ (σ = surface tension 0.072 N/m, ρ = 1000 kg/m³). The two restoring forces cross at λ_c = 2π√(σ/(ρg)) ≈ 1.7 cm, where the phase speed is minimal (≈ 0.23 m/s); rain rings sit right at this crossover (overview of the crossover and damping behaviour: https://vectree.io/pdf/c/ripple). So the ring-wave band to add is k ≈ 2π/(0.01–0.02 m) ≈ 310–630 rad/m.

**Decay.** Two mechanisms, in order of importance at this scale: (i) *geometric spreading* — a ring's energy is distributed over a growing circumference, so amplitude falls ≈ r^(−1/2) even without dissipation (https://vectree.io/pdf/c/ripple); (ii) *viscous damping*, which grows as k² and is dramatically increased by surface films/surfactants (same source) — clean-water rings survive to the observed 25–50 cm; a filmed surface kills them much sooner. A defensible per-ring envelope is A(r) = A₀ (r₀/r)^(1/2) · exp(−β k² r / c_g) with β an effective viscosity factor; beyond the functional form, β must be fit, not derived (Honesty notes).

**The second rain effect (do not omit).** Rain also generates near-surface turbulence that *damps* pre-existing waves longer than a few cm while *enhancing* the capillary band — net surface roughness rises at low wind but wind-wave energy falls (JGR paper above; also the GNSS-R rain characterisation: https://websites.umich.edu/~cruf/pubs/Balasubramaniam-Ruf_GNSS-R_Rain-Impact_RSE_2020.pdf). In the windless scene this second effect is nearly a no-op (there is nothing to damp), which is exactly why the scene is clean; but if the user later mixes rain + wind, the coupling sign is: ring band ↑, gravity band ↓.

### 3.5 How rain ripples are modelled — and the pipeline fit

Established approaches, in increasing physical fidelity:
1. **Procedural ring textures / per-drop analytic packets** (games): spawn a damped, expanding ring packet per drop (analytic in the shader, tiled or particle-based). Looks right, statistics wrong (drop sizes/counts decoupled from R). Not acceptable here given the project's "no fake patterns" rule, except as a *visualisation of* approach 2's statistics.
2. **Spectral (the remote-sensing standard):** represent the rain-roughened surface by adding a ring-wave *spectrum* to the surface spectrum. Bliven et al. (1997) published a ring-wave spectrum for simulated rain that "has been widely used" for exactly this (https://websites.umich.edu/~cruf/pubs/Balasubramaniam-Ruf_GNSS-R_Rain-Impact_RSE_2020.pdf). This is the approach that fits this pipeline natively.
3. **Direct wave-equation heightfield simulation** with drop impulses: physically appealing but at 7,000 impacts/s/patch and cm-scale rings it needs ~mm grid spacing over the patch — heavier than the current 2048² FFT state and statistically *less* controlled than approach 2.

**Implementation note (state simulation).** Add a rain band to the existing Phillips spectrum before the FFT: a spectral bump centred at k_ring ≈ 2π/(0.015 m) ≈ 420 rad/m whose *total variance* is set by drop statistics — variance per unit area ∝ F · (mean ring energy per impact) · (mean ring lifetime), with ring energy ∝ mean impact kinetic energy from §3.1–3.2 and lifetime from §3.4. The random-phase FFT machinery then produces a statistically correct, seamlessly looping ring field with the right wavelength, isotropy (rings are isotropic → the bump is direction-uniform, unlike the wind spectrum), and density — and the Jacobian transport turns it into caustics with no further changes. The one number that cannot be derived from first principles is the impact→ring energy conversion efficiency; expose it as a single calibrated constant, documented as such (Honesty notes). Rain rate R is the node; Λ(R), v(D), F(R), k_ring band, and ring variance are all derived.

---

---

## 4. The observer — a human-retina post-stage

Everything in §§1–3 produces a radiance field. What a viewer *perceives* of it is filtered by the eye. This section covers only the parts of visual physiology that can honestly become image operations, each driven by one physical input: the **adapting luminance** L_a (cd/m²) — the time-averaged luminance the eye has adapted to. That input has a hard prerequisite: a scene→display luminance scale (e.g., "this render is displayed with peak = 200 cd/m², mean = L_a"). Without that scale, pupil diameters and Weber fractions are numbers without a referent (Honesty notes).

### 4.1 Pupil diameter vs luminance

**Formula (primary).** Watson & Yellott (2012, J. Vision 12(10):12), a unified formula fitted to the major datasets, built on the Stanley & Davies core, with L = adapting luminance (cd/m²) and a = adapting field area (deg²):

  **D_SD(L, a) = 7.75 − 5.75 · [ (La/846)^0.41 / ( (La/846)^0.41 + 2 ) ]**  (mm)

The full unified form adds two corrections: field area enters through the effective corneal flux density F = L·a·M(e) (M = an eccentricity/duty factor), and age shifts the result linearly, D_U = D_SD(F, 1) + (y − 28.58)·[0.02132 − 0.009562·D_SD(F, 1)]. (Full text: https://escholarship.org/content/qt4fh0726w/qt4fh0726w.pdf?v=lg)

**Simpler form** (Moon & Spencer 1944, still widely used): D = 4.9 − 3·tanh(0.4·log₁₀L − 0.00114). Values *(computed)*: L = 0.01 cd/m² → ≈ 6.9 mm; L = 1 → ≈ 4.9 mm; L = 100 → ≈ 2.9 mm; L = 10⁴ → ≈ 2.1 mm. Overall physiological range ≈ **2–8 mm**; these are population means — individual pupils vary by roughly ±1 mm at the same luminance.

**Couplings.** Pupil area sets retinal illuminance (trolands: Td = L × pupil area in mm²), which is the quantity adaptation (§4.2) actually responds to. Pupil diameter also sets the blur scale (§4.3) and the diffraction/aberration balance. In this pipeline the pupil's image effect is therefore *indirect but real*: it modulates the defocus blur and the effective adaptation state, not the rendered radiance itself.

**Implementation note (observer post-stage).** Node: L_a (from the display scale) + optional age (default 28.6 yr, the formula's reference) → D via the unified formula. Cost: one scalar per frame (animate it slowly — real pupils lag luminance changes by ~0.2–1 s).

### 4.2 Light/dark adaptation and the Weber fraction

**Established behaviour.** Over the photopic range the just-noticeable luminance difference follows Weber's law with ΔL/L ≈ **0.01–0.02** for large uniform fields; at low luminance it crosses to the de Vries–Rose regime where ΔL ∝ √L (photon-noise limited); the transition sits around 0.1–1 cd/m². Mechanistically, adaptation models (e.g., Snippe, Poot & van Hateren 2000, Vis. Neurosci. 17) decompose into a fast divisive gain control (with the √-law regime), slower log-like gain loops, a subtractive (background-removal) stage, and contrast gain control (model structure reviewed in https://pure.rug.nl/ws/files/3118722/2000VisNeurosciSnippe.pdf; Weber-fraction values as summarised in https://library.imaging.org/admin/apis/public/api/ist/website/downloadArticle/cic/30/1/6).

**What this means as an image operation.** The honest core is a **divisive normalisation by the adaptation level**: perceived contrast ≈ (L − L_a)/L_a near L_a (Weber), compressing to a √-law at low light. This is exactly the family of operations a tone curve approximates — but here the curve is *derived* from L_a rather than chosen aesthetically.

**What cannot be included honestly.** The *time course* of adaptation (dark adaptation takes ~30 min, with separate cone/rod phases) is a property of a viewer watching changing content over minutes; a still or an 8 s loop has a single adaptation state. The node is a steady-state gain, not a simulation of adaptation dynamics.

**Implementation note (observer post-stage).** Apply after composite: L_out = f(L_in; L_a) with a Weber-divisive form and the de Vries–Rose blend below ~0.3 cd/m². Expose L_a (shared with §4.1) — no separate "adaptation" slider.

### 4.3 Retinal blur — accommodation and the eye's depth of focus

**Formula (geometric).** For a pupil of diameter p (m) and a defocus of ΔD diopters, the angular blur circle is ≈ p·ΔD radians (small-angle geometry; explainer: https://explorehall.com/en/how/human-eye/). Example *(computed)*: p = 5 mm, ΔD = 1 D → 0.005 rad ≈ 17 arcmin — far above the ~1 arcmin acuity limit, which is why defocus dominates retinal blur whenever it is present. Accommodation range (Hofstetter's formulas): amplitude ≈ 18.5 − 0.3·age D (average), 15 − 0.25·age (minimum) — a 20-year-old accommodates ~12 D, a 50-year-old ~3 D (https://www.visiontherapy.co.in/2025/09/Unit-1-accommodation-and-presbyopia-visual-optics-4th-semester-optometry.html). The eye's depth of focus is of order ±0.25–0.5 D, depending on pupil size and the blur criterion adopted (derive it as ΔD_max = θ_tol/p with θ_tol a stated criterion — it is a criterion choice, not a constant of nature).

**The flat-image problem (state plainly).** The pipeline's output is a flat image of a scene whose *depicted* depths vary (surface to bottom), but the viewer's eye accommodates to the *display*, one distance. There is no honest per-pixel defocus from the depicted water depth, because the eye cannot accommodate to a picture's contents. Two defensible options only: (i) **uniform defocus blur** — focus error at the display, blur = p·ΔD with ΔD a small stated value; (ii) a **labelled artistic pseudo-depth blur** using the water-depth map as a fake focal plane — permitted as art, forbidden as physics, and must be named as such in the node.

**Implementation note (observer post-stage).** Option (i): a single Gaussian PSF whose width comes from p (§4.1) × ΔD. Option (ii) if wanted: depth-keyed blur with the node named "pseudo-focus (artistic)".

### 4.4 Contrast sensitivity (CSF)

**Formula.** The Mannos & Sakrison (1974) CSF, f in cycles/degree:

  **H(f) = 2.6 · (0.0192 + 0.114f) · exp( −(0.114f)^1.1 )**,

used normalised; it peaks near f ≈ 8 cpd and falls to ~0 by ~60 cpd (as reproduced in Pappas et al., SPM 2003: https://www.eecs.northwestern.edu/~pappas/papers/pappas_spm03.pdf). Caveat: this is a fit at one (high) luminance; the real CSF shifts down and left as adaptation luminance falls (Barten's model is the luminance-aware form). At the pond scene's typical display luminance the fixed fit is a reasonable first order.

**Honest uses.** (i) *Analysis/display:* weight the caustic pattern's spatial-frequency content by H(f) to predict which detail a viewer actually sees — e.g., to check whether the rain-ring band (§3.4, seen at typical viewing distances) survives perception, or whether glitter sparkle is above threshold. (ii) *Perceptual shaping:* a soft, labelled contrast reshaping in the frequency domain. What the CSF must not become is a silent bake-in that "improves" the artwork — it describes a viewer's sensitivity, not the water.

**Implementation note (observer post-stage, optional).** Frequency-domain weighting on the final image, parameterised by viewing geometry (display size + distance → pixels-per-degree, which converts image frequencies to cpd). Without the viewing geometry the CSF has no frequency axis — same class of prerequisite as §4's luminance scale.

---

## 5. Variable map — every v3 variable, its stage, and what it derives from

| Variable | Stage | Input or derived | Derives from / couples to |
|---|---|---|---|
| Wind speed/direction U | state sim | input (existing) | spectrum level & peak; σ² target; glitter width |
| Significant wave height H_s | state sim | derived readout, or input in height-driven mode | m₀ of the spectrum (§1.1); capped by Stokes 1/7 |
| "Calmness" = mean square slope σ² | state sim readout + composite | derived | ∫k²S(k)dk of realised field; Cox–Munk reference (§1.2); drives glitter + shadowing α |
| Wave shadowing G₁/G₂ (Smith) | composite (statistical) *or* transport (ray-exact) — choose one | derived | realised rms slope (§1.3); sun/view geometry |
| Surface Fresnel R(θ) | transport (already implied) | **computed, no node** | n(λ) fixed; sun/view angles (§2.1) |
| Bottom albedo ρ_b(λ) | composite | input (substrate preset) | multiplies E_base,dir only (§2.1) |
| Depth z ("bottom range") | transport + composite | input (existing) | Beer–Lambert paths with c(λ) (§2.2); contrast via scattering blur (b·z) |
| Turbidity Δa, Δb | composite | input (existing) | added to pure-water a_w, b_w floors (§2.2) |
| Rain rate R | state sim | input (new scene) | Λ(R), v(D), F(R), ring band & variance (§3) |
| Gravity g | state sim | input (scene preset; §6.1) | dispersion, √(gh), growth laws, seiche periods |
| Fetch F; basin length L | state sim (scene) | input (new; §6.3) | JONSWAP α, ω_p; SMB caps; Merian modes |
| Current U_c | state sim | input (new; §6.5) | Doppler-shifted phase evolution; blocking |
| Seismic event (t₀, ground motion) | state sim (event) | event input (new; §6.6) | basin-mode ringing; tilt; slow depth modulation |
| Magnetic field B | — | **rejected** (§6.2) | no coupling exists at environmental strengths |
| Pupil D; adaptation gain; defocus blur; CSF weighting | observer post-stage | derived from L_a (+ age, viewing geometry) | display luminance scale is a prerequisite (§4) |


---

## 6. Environmental coupling — the wave field as a function of everything at once

Sections 1–4 treated variables one at a time. This section treats the couplings themselves: the places where one environmental quantity reaches inside the formulas of another. The organising fact is the full gravity–capillary dispersion relation, because nearly every coupling below is a term in it or a boundary condition on it:

  **ω² = ( gk + (σ/ρ) k³ ) · tanh(kh)**

with σ = surface tension (0.072 N/m), ρ = water density — gravity restores the long waves, surface tension the short ones, and tanh(kh) is the bottom reaching up into the dynamics (classical form as stated e.g. in https://arxiv.org/html/2602.17539v1, eq. 1, and with a current added in §6.5's source).

### 6.1 Gravity — what changing g actually does

**Where g sits.** Not in one place but in five: (i) the gravity restoring term gk above; (ii) the shallow-water speed c = √(gh) (the tanh(kh) → kh limit); (iii) every wave-growth law, through the dimensionless groups gF/U² and gH_s/U² (§6.3) and H_s ∝ U²/g (§1.1); (iv) basin-mode periods T ∝ 1/√(gh) (§6.3); (v) drop fall speeds (§3.2 is an Earth calibration).

**The crossover it controls.** Gravity and surface tension trade dominance where gk = (σ/ρ)k³, i.e., at the **capillary length**

  **l_c = √(σ / (ρg)) ≈ 2.7 mm** for water on Earth, k_c = 1/l_c ≈ 370 rad/m, crossover wavelength λ_c = 2π l_c ≈ 1.7 cm

(https://arxiv.org/html/2602.17539v1). This is the same 1.7 cm crossover the rain rings of §3.4 live on. The minimum possible phase speed of a water wave also sits there: c_min = (4gσ/ρ)^{1/4} ≈ 0.23 m/s on Earth.

**Gravity as a variable — the Moon worked example** *(all computed from the formulas above, g_Moon = 1.62 m/s²)*:
- Wave clock: ω ∝ √g at fixed k → all periods ×√(9.81/1.62) ≈ **×2.46 slower**.
- Capillary length: l_c ∝ g^(−1/2) → **6.7 mm**; crossover wavelength **4.2 cm** — "capillary" behaviour reaches up to visibly larger ripples, and the rain-ring band would sit deeper inside the capillary regime.
- Slowest wave: c_min ∝ g^{1/4} → **0.147 m/s**.
- Fully developed sea: H_s = 0.21 U²/g → **×6.05 taller** for the same wind.
- Seiche periods: ×2.46 longer.

**Honesty inside the physics.** The formulas accept any g, and inside them the change is perfectly consistent — g is a parameter of the dynamics, not of the planet's name. But two caveats are part of the variable, not footnotes: the Moon has essentially no atmosphere, so lunar *wind* waves are a counterfactual (the worked example is "same air and water, weaker gravity"); and several calibrations used elsewhere in this report are Earth measurements that do not travel with g — the Cox–Munk mss-vs-wind relation (§1.2) and the Atlas–Ulbrich fall speeds (§3.2) were measured at 1 g. The spectrum-derived σ² *readout* (§1.2) stays self-consistent at any g because it is computed from the realised field; the Cox–Munk *target* does not.

**Implementation note (state simulation).** g is already a constant baked into the spectrum construction and phase animation. Exposing it as a scene node (presets: Earth 9.81, Mars 3.71, Moon 1.62 m/s²) is a one-scalar change that propagates correctly through dispersion, √(gh) factors, growth caps and seiche modes — the single cheapest "environmental coupling" in the whole report, and the cleanest demonstration of the user's point: one variable, everywhere at once.

### 6.2 Magnetic field — the honest zero

**Question:** does a magnetic field affect waves on ordinary water? **Answer: no — not at any field strength found in nature on Earth, and the deficit is not marginal.**

**The mechanism that would have to work.** Water is diamagnetic, with volume susceptibility **χ_v = −9.05×10⁻⁶** (SI) (https://www.scienceabc.com/pure-sciences/can-a-magnetic-field-deflect-a-stream-of-water). A *uniform* field exerts no net force on a diamagnetic body at all; force requires a gradient of B², with force density f = (χ_v/2μ₀)·∇B². For that force density to merely equal water's weight density ρg requires

  B·|∇B| = μ₀ρg/|χ_v| ≈ **1.4×10³ T²/m**  *(computed)*.

**What nature provides.** Earth's field is **25–65 µT** (https://www.sollau.com/glossary-of-terms?signal=openCookiesDetailWindow), varying over planetary distances — its vertical gradient is at most ~0.03 nT/m (lecture values: https://www.karlancer.com/api/file/1747829992-QdgP.pdf). So B·∇B ≈ 5×10⁻⁵ × 3×10⁻¹¹ ≈ **1.5×10⁻¹⁵ T²/m** *(computed)* — some **18 orders of magnitude** below the level where the force would matter. Even the crudest pressure comparison fails: magnetic pressure B²/2μ₀ ≈ 10⁻³ Pa at 50 µT, vs ≈ 9.8 Pa for a single millimetre of water head *(computed)* — and a uniform field's pressure does no net work on the surface shape anyway. The most favourable alternative channel — MHD tension, since seawater conducts — scales with the Alfvén speed v_A = B/√(μ₀ρ) ≈ 1.4 mm/s at Earth's field *(computed)*, against wave phase speeds that never fall below 0.23 m/s (§3.4): two orders too slow to organise anything.

**Scale anchors.** Making water's diamagnetism *visible* is a famous laboratory stunt precisely because of the field required: levitating a frog (mostly water) took ≈ **16 T** with extreme gradients (Berry & Geim 1997, Nijmegen; an idealised textbook estimate gives B ≈ 11.7 T under the most favourable gradient assumption: http://fab.cba.mit.edu/classes/MAS.862/people/Jiaming/pdf/Week%2010%20Magnetic%20Materials.pdf). Clinical MRI fields (1.5–3 T, near-uniform) do nothing visible to water surfaces.

**Verdict for the node graph: no magnetic node for water.** Any visible effect from a "magnetic field (µT)" slider on water would be fabricated physics. If the artwork wants the *aesthetic* anyway, there is exactly one honest route: relabel the material. **Ferrofluids** have χ ~ +0.1–1 (four to five orders of magnitude larger) and genuinely grow magnetic surface spikes (the Rosensweig instability) in millitesla fields — a "magnetic field" node is real physics *for that material*, and connects naturally to the parked Blender "water held in a shape by a field" idea, which is also, physically, a ferrofluid/levitation aesthetic rather than water. Ship it, if at all, as a material switch labelled "ferrofluid (not water)".

### 6.3 Fetch, depth, and basin size — "the amount of water", translated

**The translation, plainly.** The user's intuition that "the fluctuation is governed by how much water there is" has no *volume* term anywhere in wave physics. What actually enters is three separate lengths, each governing a different phenomenon — **fetch F** (the distance wind has blown over the water), **depth h**, and **basin length L** (for standing modes) — plus wind duration t. Two basins of identical volume but different shape behave completely differently; volume is the wrong variable, and these three are the right ones.

**Fetch-limited growth (why a pond only ripples).** The SMB (Sverdrup–Munk–Bretschneider) relations, in deep-water fetch-limited form:

  gH_s/U² = 0.283 · tanh[ 0.0125 · (gF/U²)^0.42 ],   gT_s/U = 7.54 · tanh[ 0.077 · (gF/U²)^0.25 ]

(the general SPM form multiplies in a depth factor tanh[0.53·(gh/U²)^0.75]; the dimensionless equations as published in https://www.thefreelibrary.com/Novel+Ocean+Wave+Height+and+Energy+Spectrum+Forecasting+Approaches%3a...-a0771803391). At pond fetches the dimensionless fetch gF/U² is O(1–100), deep inside the tanh's linear region, so growth is strangled no matter how long the wind blows *(all computed)*:

| F | U | H_s | T_s | λ (deep-water equiv.) |
|---|---|---|---|---|
| 10 m | 5 m/s | ≈ 1.6 cm | ≈ 0.41 s | ≈ 0.27 m |
| 100 m | 5 m/s | ≈ 4.2 cm | ≈ 0.73 s | ≈ 0.84 m |
| 100 m | 10 m/s | ≈ 9.4 cm | ≈ 1.04 s | ≈ 1.7 m |

Even reaching the 100 m-fetch state takes time: the minimum duration t_min (gt/U = 68.8·(gF/U²)^0.67) is ≈ 7 min at U = 5 m/s *(computed)*. The spectrum-level version of the same physics is JONSWAP: α = 0.076·(gF/U²)^(−0.22), ω_p = 22·(g²/(UF))^(1/3), peak enhancement γ = 3.3 (Hasselmann et al. 1973; https://wikiwaves.org/index.php?title=Ocean-Wave_Spectra&diff=next&oldid=4012) — at F = 100 m, U = 5 m/s this gives T_p ≈ 0.5 s *(computed)*, vs the SMB significant period 0.73 s above; the two differ because peak period and significant period are different statistics, not because either is wrong.

**Depth limits.** Depth caps waves three ways: (i) through tanh(kh) in dispersion — long waves feel the bottom and slow toward √(gh); (ii) through the SMB depth factor above; (iii) as a hard breaking ceiling, the Miche/McCowan criterion H_b = (0.88/k)·tanh(0.78·kh/0.88), whose shallow-water limit is **H_b ≈ 0.78 h** (USACE coastal-engineering reference: https://coastalinlets.tpub.com/TR-02-22-RPT_2/TR-02-22-RPT_20024.htm). In a 0.3 m-deep pond margin, no wave taller than ≈ 23 cm can exist without breaking, at any wind speed.

**Basin modes (seiche).** A bounded basin owns a discrete set of standing-wave modes; for a rectangular basin the fundamental (Merian 1828) is

  **T_n = 2L / ( n √(gh) )**,  n = 1, 2, 3, …

(https://en.wikipedia.org/wiki/Seiche). Worked values *(computed)*: a pond L = 30 m, h = 1 m → T₁ ≈ 19 s; a garden pond L = 3 m, h = 0.4 m → T₁ ≈ 3.0 s; at the other extreme, Lake Geneva's observed fundamental is ≈ 73 min (Forel's measurements, same source). Anything that pushes the basin near these periods — gusts, pressure jumps, or ground motion (§6.6) — rings them, and they decay slowly.

**Implementation note (state simulation).** Fetch F is a *scene* node (the depicted water body's wind run, in metres; default ∞ recovers §1.1's PM spectrum). With finite F, the spectrum should be built in JONSWAP form (α(F, U), ω_p(F, U) above) rather than PM — a parameter change inside the existing FFT construction, not a new mechanism. Basin L and h jointly define the scene's Merian modes, which enter only through §6.6's event response and as a slow tilt/depth modulation: the 1.5 m FFT patch has periodic boundaries and cannot itself hold a standing mode, so the seiche must be composited as a scene-level component, not simulated in-patch.


### 6.4 Wave–wave interaction — "waves make waves", literally

**The mechanism.** In a random sea the spectral components are not independent: pairs of components exchange energy with other pairs through **resonant four-wave interactions**, allowed only where the quartet satisfies both resonance conditions

  **k₁ + k₂ = k₃ + k₄,   ω₁ + ω₂ = ω₃ + ω₄**   (ω(k) = the dispersion relation of §6).

Phillips (1960) first described the transfer; Hasselmann (1962) and Zakharov (1968) cast it as a kinetic equation for the wave-action density whose source term is a sixfold integral with two δ-functions enforcing conservation of energy, wave action, and momentum (formulation and computational survey: https://www.ecmwf.int/sites/default/files/elibrary/2012/12895-efficient-algorithms-non-linear-four-wave-interactions.pdf). Its observed consequences are the sea's self-organisation: the spectral peak **downshifts** to lower frequencies as the sea develops (an inverse cascade of wave action), the forward face of the spectrum grows, and the spectral shape stabilises. The JONSWAP analysis (§6.3) attributed the downshift and forward-face growth to exactly this transfer — and self-similar solutions of the *conservative* kinetic equation (wind input and dissipation dropped) reproduce the JONSWAP fetch-growth laws quantitatively (Zakharov, fetch-limited theory: http://zakharov75.itp.ac.ru/static/local/zve75/zakharov/2002/2002-03-1BanffFetch02.pdf). In other words, the empirical fetch laws of §6.3 *are* the fingerprint of wave–wave interaction. (In shallow water the dispersion curve straightens toward ω ∝ k, and three-wave (triad) resonances switch on as well — *reasoned from the geometry of the dispersion relation*: a linear ω(k) makes collinear triads satisfy both resonance conditions identically.)

**What this means for the pipeline — stated honestly.** The FFT state is a linear superposition with randomised phases: by construction its components *never* exchange energy. That is not a hidden defect. PM and JONSWAP spectra are **post-interaction shapes** — the interactions have already happened, in the ocean data the spectra were fitted to. Using an evolved spectrum is the physically correct shortcut, and it is the only feasible one: evaluating the kinetic integral is a sixfold integral that operational wave models themselves replace with an approximation (the DIA) for cost (ECMWF source above). Three honest in-pipeline traces of the nonlinearity remain available:
- *Use evolved spectra* (JONSWAP with §6.3's fetch laws) rather than young-sea shapes — done by construction once F is a node.
- *Bound-wave sharpening:* the horizontal (choppy/Gerstner-type) displacement already used in the offline pipeline sharpens crests and flattens troughs — that asymmetry is the leading bound-harmonic signature of Stokes waves, i.e., the visible face of weak nonlinearity, and it is already in the transport's input.
- *Steepness caps* (§1.1's H/λ ≈ 1/7) keep the linear field from reaching slopes the real, interacting field would have broken away.

Not representable: modulational (Benjamin–Feir) instability and genuine spectral evolution over time — both need long evolution of a narrow-banded field, which a static-spectrum filter does not have. Noted here so the absence is a documented decision, not an oversight.

### 6.5 River current — the same water, in a moving frame

**Dispersion with a uniform current U.** Frequencies measured from the bank (fixed frame) obey

  **(ω − k·U)² = ( gk + (σ/ρ) k³ ) · tanh(kh)**

(Rousseaux et al., "Some aspects of dispersive horizons": https://germain-rousseaux.cnrs.fr/wp-content/uploads/2020/07/AG7C2013.pdf). The combination σ_int = ω − k·U is the intrinsic frequency — the only one the water itself knows about. Everything else is bookkeeping of the observer's frame, but the bookkeeping has visible teeth:

- **Doppler and reshaping at fixed ω:** a following current stretches waves (k falls) and flattens them; an opposing current shortens them and *raises* them — "when a gravity wave meets a countercurrent, the incident wavelength diminishes and the wave height increases" (same source). Wave-action conservation is the mechanism: energy flux piles up where the current slows the energy's escape.
- **Blocking:** when the opposing current reaches the point where the absolute group velocity U + c_g passes through zero, wave energy can no longer propagate upstream at all — the waves are blocked, steepen, and break or reflect. Evaluating the blocking point from the dispersion relation for an incident deep-water wave *(computed)* gives |U_block| = g/(4ω) = ½·c_g0: a T = 1 s wave is blocked by an opposing current of only ≈ **0.39 m/s**. Worse for small texture: the capillary crossover band (§6.1) has c_g ≈ 0.23 m/s at its minimum, so currents of a few tens of cm/s — ordinary stream flow — already block and reorganise centimetre ripples. A river's visible surface texture is governed by its current field first and its local wind second; the §3 rain rings, conversely, are a *still-water* phenomenon and wash out under even gentle flow.

**Implementation note (state simulation).** A uniform current is almost free and exactly this cheap: advance each spectral component's phase with the fixed-frame frequency ω(k) = σ_int(k) + k·U_c instead of σ_int(k), and advect the whole pattern at U_c. The caustics inherit the consequences automatically, because the heightfield statistics genuinely changed (compressed, steepened ripples against the flow). What uniform-U does *not* cover: spatially varying current (shear, eddies, confluences), which refracts the wave field toward low-speed regions — that requires a ray/WKB step on top of the spectrum and is out of scope for v3; the node is one uniform vector, labelled as such.

### 6.6 Earthquakes and plates — the event-type coupling

An earthquake is not a continuous environmental variable like wind; it is an **impulse** delivered to the whole system at once. It couples to water through two completely different doors, and keeping them apart is most of the physics.

**6.6.1 The ocean door — tsunamis (the bottom moves).** At a subduction zone the plates are locked while convergence (cm/yr) accumulates strain for decades to centuries; rupture releases it in seconds to minutes, and the seafloor jumps vertically by metres over hundreds of kilometres. The rupture is much faster than water can flow aside, so the sea surface initially *copies the deformed seafloor* and the entire water column above it is lifted or dropped — plate motion is the ultimate source, the earthquake merely its release (mechanism steps: https://thesea.org/how-do-tsunamis-form/). Because the wavelength (tens to hundreds of km) vastly exceeds the depth, a tsunami is a shallow-water wave even mid-ocean:

  c = √(gh) ≈ 198 m/s ≈ 713 km/h at h = 4000 m; deep-ocean amplitude typically < 1 m; period (10–60 min) unchanged in propagation

It grows only where the depth gives it back: shoaling follows **Green's law, a ∝ h^(−1/4)** — 4000 m → 10 m multiplies amplitude ×≈4.5 — with further amplification where a bay's own seiche period (§6.3) matches the incoming period (worked numbers: https://bpsc.careercraftlabs.org/read/gate-xe7-ocean-physical-ch03; Green's-law derivation from energy-flux conservation: https://vectree.io/pdf/c/tsunami).

**6.6.2 The lake door — seismic seiches (the ground moves).** In an enclosed basin nothing arrives *through the water* — the basin's floor and shores themselves oscillate as seismic waves pass (body waves, then surface waves), and horizontal ground motion pumps the basin's standing modes: exactly the Merian modes of §6.3, T_n = 2L/(n√(gh)). The response is resonant and lightly damped: modes whose period matches the ground motion's energetic periods ring up far beyond the ground's own amplitude and then decay over hours. Tank "sloshing" in engineering is the same phenomenon at appliance scale.

**6.6.3 The documented far-field case — the user's butterfly effect, with numbers.** On 11 March 2011, about 30 minutes after the M9.0–9.1 Tōhoku earthquake, fjords in western Norway — more than **8,300 km** away — began to seiche: trough-to-peak amplitudes **1.0–1.5 m**, periods **67–100 s**, oscillations lasting ~3 hours, *only* in fjords whose axes pointed SW–NE, toward the epicentre. The trigger was not the (later) surface waves but the passage of horizontal **S waves**, with ground displacement of only ~**1 cm** — reproduced by the authors with a shallow-water model of the fjords (Bondevik, Gjevik & Sørensen 2013, *Geophys. Res. Lett.* 40, doi:10.1002/grl.50639; abstract and reporting: http://ktwop.com/tag/fjords/ ; https://www.earthmagazine.org/article/giant-quake-sloshed-fjords-half-world-away/ ; https://www.livescience.com/38960-japan-earthquake-norway-seiche-wave.html). The amplification — centimetre ground motion → metre-scale water oscillation — is a factor of ~**100** *(computed from the reported values)*: resonance, not force. It was also not a first: the very term "seismic seiche" was coined by Kvale (1955) for lake oscillations in Norway and England caused by the 1950 Assam (M8.6) earthquake (same sources). An earthquake on the other side of the planet really does shake a pond — provided the pond's own period is listening at the right frequency.

**6.6.4 Implementation note — earthquakes are events, not sliders.** The defensible formulation stops at *local ground motion*: the chain magnitude → ground motion at the scene → basin response has two empirical links (ground-motion prediction, basin damping), so the event node takes the local input directly:

  event = { t₀, peak ground displacement (or acceleration), dominant period T_g, duration }

In the basin's frame, horizontal ground acceleration a_h acts as a tilting of effective gravity — the equilibrium surface tilts by ≈ a_h/g (the same g-vector that §6.1 exposes) — and each Merian mode responds as a driven, lightly damped oscillator, ringing at T_n with resonant build-up when T_g ≈ T_n. In-patch, a ~1.5 m periodic FFT tile has no shores and no true modes, so the honest in-frame signatures are scene-level: (i) a slow whole-frame **tilt oscillation** at the scene basin's T₁ (scene L, h from §6.3); (ii) a slow **depth modulation** h(t) at the same period, feeding the Beer–Lambert paths of §2.2; (iii) the tilt entering the surface-slope statistics. The Norway factor-~100 is the calibration anchor for the resonant amplitude — order-of-magnitude, flagged as calibrated, since basin damping (seiche Q) is empirical at every scale. And a tsunami inside a frame is simpler still, and must not be faked as waves: at frame timescales it is a slow water-level ramp (the depth node moving over minutes) plus, near shore, a strong current (§6.5). No ripple-field change.

### 6.7 Synthesis — where the sensitivity actually lives

Reading §§6.1–6.6 as one system, each environmental lever moves the field through a distinct first domino:

| Lever | First domino | What falls next |
|---|---|---|
| g | dispersion clock (√g); crossover l_c (g^(−1/2)) | all periods; which ripples are "capillary"; growth H_s ∝ 1/g; seiche T |
| Wind U | spectrum level (H_s ∝ U²), peak (ω_p ∝ 1/U) | σ² (∝ U) → glitter width, shadowing; steepness caps engage |
| Fetch F | growth ceiling (tanh of (gF/U²)^0.42) | pond vs ocean regime; JONSWAP α, ω_p |
| Depth h | tanh(kh); breaking cap 0.78 h | wave speeds; Beer–Lambert paths; seiche T (h^(−1/2)) |
| Current U_c | Doppler term k·U_c | wavelength/height reshaping; blocking threshold |
| Rain R | isotropic capillary band added | glitter texture on calm water; gravity band damped |
| Seismic event | basin modes rung at T_n | tilt + depth oscillation for hours from a 1-minute cause |
| Magnetic B | — | **nothing: the system's one true zero** (§6.2) |

Two properties deserve to be said plainly, because they are the accurate version of the user's "environmental butterfly effect". First, the couplings are mostly smooth power laws — double the wind and the field changes a lot, but *predictably*. Second, the exceptions are concentrated at **thresholds and resonances**, and that is where small causes flip regimes: blocking (U_c ≈ c_g, §6.5), seiche resonance (T_g ≈ T_n, §6.6), breaking caps (steepness 1/7; H ≈ 0.78 h), and the gravity–capillary crossover where rain rings live (§§3.4, 6.1). A centimetre of ground motion moving a fjord by a metre is not chaos; it is a resonance doing exactly what resonances do. And the linear wave field itself is not chaotic in the Lorenz sense — the butterfly here is metaphor for *coupling density plus a few cliff edges*, and the artwork is better served by the true version: every knob in §§1–6 either does nothing (magnetism), moves everything smoothly (gravity, wind), or matters enormously at one specific coincidence (fetch caps, blocking, resonance).


---

## 7. Honesty notes — what resists exact definition, and the closest defensible proxy

- **"Calmness" has no physical unit.** Proxy adopted: Cox–Munk mean square slope σ² (§1.2) — measurable, derived from the pipeline's own spectrum, and the quantity that actually controls glitter and caustic texture. Any 0–1 "calmness" display is a labelled transform of σ², not physics.
- **PM, SMB, and JONSWAP are empirical fits, not laws.** PM assumes a fully developed open ocean (§1.1); SMB/JONSWAP were fitted to ocean and North Sea data. Applying their dimensionless forms at pond scale (§6.3) is defensible *scaling*, and the pond numbers are computed from those forms — they are not measurements of ponds.
- **"The amount of water" is not a variable.** No volume term exists in wave physics. The user's intuition resolves into fetch F, depth h, and basin length L (§6.3); equal volumes with different shapes behave differently.
- **Magnetic field is not a proxy problem — the effect is zero.** At environmental field strengths the diamagnetic coupling sits ~18 orders of magnitude below relevance (§6.2). There is no "approximately right" magnetic slider for water; the only honest magnetic node is a ferrofluid material relabel.
- **Gravity travels only part-way.** The dispersion and growth *formulas* accept any g (§6.1), but Cox–Munk mss targets, PM constants, and raindrop fall speeds are Earth measurements. A Moon-g render is self-consistent scaling, not validated lunar oceanography — and the real Moon has no atmosphere to make wind waves with.
- **Wave–wave interaction is inherited, not simulated.** PM/JONSWAP spectra are post-interaction equilibria; the Hasselmann kinetic integral is not evaluated (§6.4). Spectral evolution over time, and modulational instability, are absent by documented decision.
- **Marshall–Palmer is an average, not a law.** Fixed N₀, stratiform rain; real drop-size distributions vary event to event (§3.1). It grounds the *statistics* of the rain scene, not any given minute of it.
- **Two rain numbers must be calibrated, not derived:** the ring-wave decay factor β and the impact→ring energy conversion efficiency (§§3.4–3.5). Each is a single documented constant; everything else in the rain chain (Λ(R), v(D), F(R), the k-band) is derived.
- **Bottom albedo mixes two different quantities in the literature.** Broadband retrieved albedos (sand ≈ 0.24–0.30) and peak *spectral* values (sand up to ≈ 0.63) must not share a node (§2.1). Mud/organic bottoms (≈ 0.05) have no canonical table; substrate presets are literature ranges.
- **There is no caustic-contrast-vs-depth law.** The model is restricted to the two channels that exist: spectral Beer–Lambert attenuation and turbidity-driven scattering blur (§2.3). A "contrast" slider would be fiction.
- **Bottom shadowing has two non-combinable treatments.** Smith G₁ is a statistical proxy; ray-exact visibility costs a heightfield intersection test per quadrature sample (§1.3). Choose one; applying both double-counts.
- **The seismic chain has empirical links.** Magnitude → local ground motion → basin amplitude cannot be one formula; the event node honestly takes *local ground motion* as input, and seiche damping (Q) plus the resonant amplification (Norway anchor: ~×100) are calibrations, not derivations (§6.6).
- **Uniform current only.** The Doppler treatment (§6.5) is exact for uniform U_c; shear refraction by spatially varying currents is not modelled.
- **The retina formulas are population means on a calibrated stage.** Individual pupils vary ~±1 mm; Watson–Yellott and Mannos–Sakrison are fits (the CSF at a single luminance). The whole observer stage is uncalibrated without a scene→display luminance scale and viewing geometry (§4). A flat image admits only *uniform* defocus honestly; depth-keyed blur is a labelled artistic device. Adaptation is a steady-state gain here — its real time dynamics (minutes) do not apply to a still or a short loop.
