# Validating the v3.6 in-water scattering model against published measurements

**Question.** v3.6 (run_v12*) splits turbidity into a(λ)/b(λ), uses a Henyey–Greenstein
(HG) phase function with g = 0.85 for all particle presets, single scattering, and a
halo whose width σ comes from the scattering integral via the HG second moment
⟨θ²⟩ = 0.3681 rad². At the preset-1 reference scene (b ≈ 0.075 m⁻¹, c ≈ 0.163 m⁻¹ green,
L = 1.25 m) this gives σ = 0.0927 m = 79.09 px carrying H ≈ 0.08 of the light.
The audit rally's charge: the gate battery only compared the model to itself.
This report asks what *published measured* data says. All computed numbers below were
produced by numerical integration of the published phase functions/tables
(scripts in `notes/pf_calc.py`, `notes/pf_calc2.py`); the integrator is validated in §6.

---

## Verdict (summary)

1. **HG with g = 0.85 is not defensible for the algae/detritus-type presets
   (`clear`, `algae_pond`).** Every measured particle phase function in the published
   corpus has ⟨cos θ⟩ between ≈ 0.92 and 0.96. Petzold's average-particle function —
   the canonical benchmark — has ⟨cos θ⟩ = 0.924 and backscatter fraction B = 0.0183.
   HG at g = 0.85 has B = 0.0361 (2.0× too much backscatter) and puts 23.8% of the
   scattered light beyond 30°, where measured functions put 6–11%.
2. **The halo is too wide in second moment, and the wrong shape.** Measured-class
   phase functions give ⟨θ²⟩ = 0.09–0.21 rad² vs our 0.3681 rad², so at fixed (b, c, L)
   the moment-matched σ shrinks to **× 0.51–0.69 of the current value**
   (79.09 px → ≈ 40–55 px at the reference scene; Petzold-benchmark value 54.9 px).
   But the shape error is worse than the width error: measured functions put
   25–42% of scattered light within 1° of forward; HG-0.85 puts 0.6% there, and
   *no* single Gaussian (even with the corrected ⟨θ²⟩) puts more than 0.2% there.
   The true single-scatter kernel is a bright near-forward core plus a faint broad
   skirt — a power-law object (measured PSF ∝ θ^−m, §3), not a Gaussian.
3. **The halo energy fraction H ≈ 0.08 is not a phase-function quantity** — in this
   model class H = (b/c)(1 − e^(−cL))·ρ depends only on a, b, L. The phase-function
   data can neither confirm nor refute it; it stands or falls with the preset a/b
   values (out of scope here, flagged in the v12 report itself as modelling choices).
4. **The sediment preset is the exception.** For high-index mineral particles the
   Fournier–Forand (FF) envelope (n = 1.15–1.20, μ ≈ 4.0) gives ⟨cos θ⟩ ≈ 0.72–0.77
   and ⟨θ²⟩ *larger* than HG-0.85 (σ × 1.29–1.44). g = 0.85 sits inside that
   envelope. This rests on the FF parameter envelope for small mineral particles,
   not on a direct sediment VSF fit I could verify — flagged as such in §4.
5. **What the measurements support instead:** a Fournier–Forand phase function.
   For a generic pond preset, the Petzold-average fit (n = 1.10, μ = 3.5835,
   B = 0.0183, ⟨cos θ⟩ = 0.930, ⟨θ²⟩ = 0.177 rad²) is the defensible default;
   preset-specific FF parameters are given in §4. If HG must be kept for code
   simplicity, g = 0.924 (⟨cos⟩-matched) or g = 0.9185 (B-matched) fixes the
   moments (σ × 0.71–0.73) but still fails the small-angle shape clauses (§5).

---

## 1. Petzold (1972): what is accessible, and the numbers

**The dataset.** Petzold, *Volume Scattering Functions for Selected Ocean Waters*,
SIO Ref. 72-78 (1972), measured VSFs at several sites. The benchmark derived from it —
the **average-particle phase function**, an average of his three San Diego Harbor
measurements with the pure-water contribution removed — is tabulated in Mobley (1994),
*Light and Water*, and in Mobley et al. (2002), and is the standard test function in
HydroLight. I did **not** obtain Petzold's original report or the full
average-particle table verbatim (paywalled/out of print); the pointwise values quoted
below are from Jonasz & Fournier's published geometric-average tables, whose underlying
data collection includes Petzold's sets. Flagged accordingly.

**Scalar invariants of the Petzold average-particle function (verified, two
independent secondary sources each):**

| quantity | value | source |
|---|---|---|
| ⟨cos θ⟩ (asymmetry parameter) | **0.924** | Solonenko & Mobley (2015): "the averaged Petzold phase function … The corresponding phase function average cosine g is 0.924"; same value proposed as the generic HG g in the UWOC literature citing Mobley (1994) |
| backscatter fraction B = b_b/b | **0.0183** | Mobley, Sundman & Boss (2002), §2A: numerical integration of the tabulated function over 90–180° |
| shape landmarks | minimum near ψ ≈ 140°, rising toward 180°; still rising at the smallest measured angles | Mobley et al. (2002), §2A |
| flat-back equivalent | for ψ ≥ 90° the function averages to β̃ ≈ B/2π = 0.00291 sr⁻¹ | Mobley et al. (2002), §2B construction |

**Per-water-type ⟨cos θ⟩ from Petzold's VSFs** (total VSF, water included; secondary
citation — Aras et al. 2013 citing their ref. [33]): clear ocean **0.8708**,
coastal ocean **0.9470**, turbid harbor **0.9199**. Note that even the lowest of
these — clear ocean, where symmetric molecular scattering drags the average down —
is above 0.85, and these are *total* VSFs; the particle-only values are higher.

**Pointwise measured particle phase functions** (Jonasz & Fournier 2007, Table 4.2 —
geometric averages of 108 "clear" and 161 "turbid" measured datasets, Petzold's
included; β̃ in sr⁻¹):

| θ | clear | turbid | θ | clear | turbid |
|---|---|---|---|---|---|
| 0.2° | 89.6 | 687 | 10° | 0.0700 | 0.630 |
| 0.5° | 13.9 | 317 | 45° | 0.00142 | 0.0149 |
| 1° | 6.38 | 83.4 | 90° | 0.000118 | 0.00214 |
| 2° | 1.94 | 21.4 | 140° | 0.000106 | 0.00122 |
| 5° | 0.275 | 2.42 | 170° | 0.000213 | 0.00117 |

with **⟨cos θ⟩ = 0.962 (clear) and 0.948 (turbid)** by numerical integration
(Jonasz & Fournier, Table 4.2, note d), and b = 0.0839 / 0.778 m⁻¹ for the two
composite waters. Caveat, stated by the authors themselves: data below ~10° are
scarce in the underlying collection, so the printed small-angle points of the
*clear* average in particular are lower bounds (see §6) — the authors' own FF fits
to the same data (§2) sit well above the printed 0.2° point. The broader measured
range across natural waters is ⟨cos θ⟩ ≈ 0.89–0.95 for whole-seawater VSFs
(Jonasz & Fournier §6.3.2: geometric averages 0.89 clear / 0.95 turbid; literature
range quoted there 0.7–0.97, the low end being the clearest waters' *total* VSF).

**Comparison with our HG at the same angles** (this work, §6): HG-0.85 gives
β̃(0.2°) = 6.5, β̃(1°) = 6.4, β̃(5°) = 4.5, β̃(10°) = 2.1, β̃(90°) = 0.0098 sr⁻¹.
It saturates below a few degrees (Mobley et al. 2002: the one-term HG "levels out
for ψ less than a few degrees … neither behavior is seen in the measured Petzold
phase function"), is ~30× under the turbid measured value at 1°, and is ~4.6× over
the turbid measured value at 90°.

## 2. Fournier–Forand: fitting, parameters, and the small-angle region

**Form.** FF (Fournier & Forand 1994; Fournier & Jonasz 1999) is the anomalous-
diffraction phase function of a Junge (hyperbolic) size distribution,
N(r) ∝ r^(−μ), μ typically 3–5 in natural waters, with particle relative refractive
index n. It is analytic, integrates to a closed-form backscatter fraction, and is
the ocean-optics community's standard replacement for HG (Mobley, Ocean Optics Web
Book: FF "does a much better job of reproducing the shapes of ocean phase functions
than does the Henyey–Greenstein phase function … especially at very small angles";
it is HydroLight's default family).

**Fits to measured data (published parameters):**

| target | n | μ | B | source |
|---|---|---|---|---|
| Petzold average-particle (best fit) | 1.10 | 3.5835 | 0.0183 | Mobley et al. (2002); Ocean Optics Web Book |
| Petzold average-particle (alt. fit) | 1.16 | 3.4586 | — | RTE review chapter (secondary; flagged) |
| Jonasz clear average | 1.098 | 3.37 | (fit to β, error 0.170 log) | Jonasz & Fournier (2007), Table 4.3 |
| Jonasz turbid average | 1.073 | 3.59 | (fit error 0.104 log) | Jonasz & Fournier (2007), Table 4.3 |
| illustrative: B = 0.01 | 1.08 | 3.483 | 0.01 | Ocean Optics Web Book, Table 1 |
| illustrative: B = 0.1 (small mineral particles) | 1.175 | 4.065 | 0.1 | Ocean Optics Web Book, Table 1 |

The (n, μ) pair for a given B is not unique — any pair on the constant-B contour
gives a similar function (Mobley et al. 2002, Fig. 2) — so B and ⟨cos θ⟩, not the
individual parameters, are the robust fit targets.

**Small-angle behaviour — the region that sets our halo.** Measured VSFs follow a
power law β̃ ∝ θ^(−s) with s ≈ 1.6 over 0.1–6° and s ≈ 2.6 over 12–100°
(Jonasz & Fournier, Fig. 4.27 discussion). FF reproduces the forward power law;
HG does not — it flattens to a constant β̃(0) = (1+g)/(4π(1−g)²)·… ≈ 6.5 sr⁻¹ at
g = 0.85. Concretely (this work, §6), the fraction of scattered power inside a
given angle:

| function | F(<0.5°) | F(<1°) | F(<5°) | F(<10°) | F(>30°) |
|---|---|---|---|---|---|
| HG g = 0.85 (v3.6) | 0.002 | 0.006 | 0.129 | 0.346 | 0.238 |
| HG g = 0.9185 (B-matched) | 0.005 | 0.021 | 0.315 | 0.586 | 0.125 |
| HG g = 0.924 (⟨cos⟩-matched) | 0.006 | 0.024 | 0.342 | 0.611 | 0.116 |
| FF Petzold fit (1.10, 3.5835) | 0.172 | 0.254 | 0.559 | 0.712 | 0.105 |
| FF Jonasz turbid fit | 0.202 | 0.297 | 0.625 | 0.770 | 0.076 |
| FF Jonasz clear fit | 0.324 | 0.418 | 0.700 | 0.818 | 0.058 |

No value of g rescues HG here: raising g to match ⟨cos θ⟩ or B still leaves
F(<1°) an order of magnitude below the measured-class functions, because the
defect is the functional form (finite forward value vs power-law divergence),
exactly as Mobley et al. (2002) describe. (Their two-term HG fares no better as a
shape: ≥3× below Petzold over 10–110° and 7× above it at 180°.)

## 3. Measured PSF / beam-spread data: what exists and what it can check

* **Mertens & Replogle (1977)**, J. Opt. Soc. Am. 67, 1105–1117 — definitions of
  the in-water PSF and BSF and their use in imaging analysis; **Mertens & Phillips
  (1972)**, DTIC report — the first in-situ BSF/PSF measurements (Bahamas), angular
  range <1 mrad to 180°, paths up to 9 attenuation lengths. Their datasets are
  inside Jonasz's averages (§1), so they reach us through Table 4.2.
* **Duntley (1971)**, SIO Ref. 71-1 — extensive laboratory BSF measurements in
  simulated ocean water (Aluminox scatterer, Nigrosin absorber) with controlled
  a and c, reduced to an empirical BSF formula in (cτ, a/c) (reproduced in
  Voss 1991). A runnable comparator in principle, but for a collimated-beam
  geometry, not our sun→bottom→eye field.
* **Voss (1991)**, Appl. Opt. 30, 2647–2651 (+ Voss & Chapin 1990) — in-situ PSF
  at three ocean sites (Pacific coast, Tongue of the Ocean, Sargasso Sea),
  500 nm, camera resolution 0.6 mrad. Findings: over **4–100 mrad the PSF is a
  power law, PSF(θ) ∝ θ^(−m)** — not an exponential, not a Gaussian — with the
  slope m falling from ≈ 1.6 at optical depth τ = cR → 0 toward ≈ 0.4–0.6 at
  τ ≈ 8–10 (per-location regressions m(τ), r² = 0.905–0.964; m uncertainty ±0.3;
  slope prediction error 10–14%). In the τ → 0 limit the PSF slope is the
  single-scatter (Wells 1973) image of the phase function's small-angle power
  law — the same s ≈ 1.6 as §2.
* IOCCG lecture material (Twardowski, SLS-2016) restates the community numbers:
  typically only 0.3–3% of b is backscattered (B ≈ 0.003–0.03), consistent with §1.

**Applicability to our halo — stated plainly.** No published dataset measures
"the halo of a bottom-reflected caustic field"; PSF/BSF data are beam/point-source
geometries and at most ranges include multiple scattering. Our reference scene is
at cL ≈ 0.20 (green) — deep in the single-scatter regime — so the PSF data check
our model at the level of *form* (power-law kernel vs our Gaussian) and, via the
τ → 0 limit, at the level of the phase function itself. The decisive, fully
runnable check is therefore the phase-function comparison of §5; the PSF
literature independently confirms that the kernel family we should produce is
heavy-tailed: across Voss's 4–100 mrad band a measured PSF falls by
(100/4)^(−1.6) ≈ 1/170, while a Gaussian with our σ (angular rms
√⟨θ²⟩ = 0.607 rad) falls by only ≈ ×2.5 across the equivalent band, and even an
FF-moment Gaussian (rms 0.421 rad) by only ≈ ×6.

## 4. Verdict per preset, and the quantified halo change

Moments computed in this work (integrator validated in §6). σ scaling is at fixed
(b, c, L): σ_new = σ_current × √(⟨θ²⟩_new / 0.3681). H does not change.

| phase function | ⟨cos θ⟩ | B | ⟨θ²⟩ (rad²) | σ ratio | σ at ref. scene |
|---|---|---|---|---|---|
| **HG g = 0.85 (current, all presets)** | 0.850 | 0.0361 | **0.3681** | 1.000 | 79.09 px (0.0927 m) |
| HG g = 0.9185 | 0.9185 | 0.0183 | 0.1977 | 0.733 | 58.0 px |
| HG g = 0.924 | 0.924 | 0.0170 | 0.1842 | 0.707 | 55.9 px |
| **FF Petzold fit (1.10, 3.5835)** | 0.930 | 0.0183 | **0.1772** | **0.694** | **54.9 px (0.0643 m)** |
| FF alt. Petzold fit (1.16, 3.4586) | 0.917 | 0.0223 | 0.2100 | 0.755 | 59.7 px |
| FF Jonasz turbid fit (1.073, 3.59) | 0.949 | 0.0127 | 0.1283 | 0.590 | 46.7 px |
| FF Jonasz clear fit (1.098, 3.37) | 0.962 | 0.0088 | 0.0938 | 0.505 | 39.9 px |
| FF algae-type (1.06, 3.5) | 0.968 | 0.0073 | 0.0793 | 0.464 | 36.7 px |
| FF detritus-type (1.08, 3.483) | 0.958 | 0.0101 | 0.1053 | 0.535 | 42.3 px |
| FF mineral (1.15, 4.0) | 0.772 | 0.0776 | 0.6114 | 1.289 | 102 px |
| FF mineral (1.175, 4.065) | 0.719 | 0.1001 | 0.7621 | 1.439 | 114 px |

Per preset:

* **`algae_pond` (the default medium preset): g = 0.85 REFUTED.** Algal/detrital
  assemblages are the low-index, low-B end of the measured corpus
  (⟨cos θ⟩ ≈ 0.95–0.97, B ≈ 0.007–0.013 in the FF representations). Adopt FF with
  n ≈ 1.06–1.10, μ ≈ 3.5 (σ → ×0.46–0.69), or HG g ≈ 0.93–0.95 if the form must
  stay (moments only). 
* **`clear`: g = 0.85 REFUTED** on the same evidence (Jonasz clear: ⟨cos θ⟩ = 0.962;
  and note the preset's b includes a small pure-water part whose own phase
  function is near-symmetric — at the preset's loads this is a <4% effect on b,
  second order, but it pushes the *effective* g down, not up toward justifying
  0.85 for the particle part).
* **`sediment`: g = 0.85 NOT refuted, but not confirmed either.** Mineral
  particles' higher n moves FF toward larger B and larger ⟨θ²⟩; at
  (n = 1.15, μ = 4.0) ⟨cos θ⟩ = 0.772 and σ would *grow* ×1.29. HG-0.85 lies
  between the biogenic and mineral FF families for every moment. What is
  missing is a published FF fit or ⟨cos θ⟩ for a measured sediment-dominated
  freshwater VSF; until one is produced, sediment's g should be labelled an
  assumption with an envelope [FF(1.10, 3.58), FF(1.15, 4.0)], not a measurement.
* **Freshwater specifically:** the measured corpus spans marine, estuarine,
  riverine (Charles River, Back River) and lake (Baikal) waters, and the
  phase-function shape is famously conservative across them (Jonasz & Fournier:
  the similarity "cuts across various water bodies, seasons, and wavelengths").
  I found **no published freshwater particle phase function with ⟨cos θ⟩ anywhere
  near 0.85**; the v12 report's provenance line ("band [0.8, 0.95] from Petzold's
  measured VSFs") is not supported by Petzold's values as published
  (0.87/0.92/0.95 per water type, total-VSF basis; 0.924 particle average).

**Shape consequence, quantified.** Matching ⟨θ²⟩ is necessary but not sufficient:
the pipeline spreads all of H as one Gaussian. Encircled fractions of the
scattered light within angle θ₀ — true FF-Petzold kernel vs the moment-matched
Gaussian: at 5°, 0.559 vs 0.042; at 10°, 0.712 vs 0.158; at 20°, ≈ 0.85 vs 0.497.
So even after the σ correction, a single-Gaussian halo puts ~half the halo light
beyond 20° where the measured kernel puts ~15%, and misses the near-forward core
that carries most of the energy. The measurement-supported kernel is two-scale
(narrow core + broad skirt) or, better, the FF angular CDF used directly in the
displacement integral in place of the ⟨θ²⟩ moment. This is a modelling decision
for the chunk that implements the fix; the comparison in §5 scores both the
moment clauses and the shape clauses so the gap is measured either way.

## 5. The exact comparison to run (published tables only)

**Inputs (all published, no lab access):**
FF formula (Fournier & Jonasz 1999 form, Ocean Optics Web Book) at the §2
parameter sets; Jonasz & Fournier Table 4.2 pointwise values for θ ≥ 1°;
scalars: Petzold average ⟨cos θ⟩ = 0.924, B = 0.0183.

**Step 1 — phase-function scorecard.** For the candidate phase function β̃
(whatever the pipeline adopts, per preset), compute by quadrature on a
log-spaced grid θ ∈ [10⁻⁶, π]:
⟨cos θ⟩, B, ⟨θ²⟩, F(<1°), F(<5°), F(<10°), β̃(1°), β̃(5°), β̃(10°), β̃(90°).
Pass clauses for the generic/algae preset (targets = Petzold average / FF-Petzold):

| clause | target | tolerance | HG-0.85 today | HG-0.924 | FF-Petzold |
|---|---|---|---|---|---|
| ⟨cos θ⟩ | 0.924 | ±0.02 | 0.850 ✗ | 0.924 ✓ | 0.930 ✓ |
| B | 0.0183 | ±0.005 | 0.0361 ✗ | 0.0170 ✓ | 0.0183 ✓ |
| ⟨θ²⟩ | 0.177 rad² | ±15% | 0.368 ✗ (×2.08) | 0.184 ✓ | 0.177 ✓ |
| F(<5°) | 0.559 | ±0.10 | 0.129 ✗ | 0.342 ✗ | 0.559 ✓ |
| β̃(90°) vs Jonasz turbid table | 0.00214 sr⁻¹ | ×1.5 | 0.00977 ✗ (×4.6) | 0.00461 ✗ (×2.2) | 0.00420 ✗ (×2.0)* |

*FF-Petzold's β̃(90°) = 0.00420 vs the turbid *table* point 0.00214: the Petzold
average and the Jonasz turbid average are different composites (the Petzold
flat-back level is 0.00291 sr⁻¹); score this clause against the Petzold
flat-back level with tolerance ×1.5 → 0.00420 vs 0.00291 = ×1.44 ✓. Record which
reference was used; do not silently mix composites. (HG-0.85 fails against
either reference: ×3.4 vs flat-back.)

**Step 2 — halo-level check.** Recompute σ through the pipeline's own formula
(σ² = (⟨θ²⟩/2)(b/c³)D2(cL)) for the v12 G1 reference scene and the four G3 sweep
rows, with the candidate ⟨θ²⟩. Required: σ scales exactly as √(⟨θ²⟩ ratio)
(it is the same formula), H is bit-identical to the v12 run (H must not move —
any change is a code leak, since H is phase-function independent), and the
delivered stills' measured halo share H/(T+H) is unchanged while the peak
profile at fixed lags (the G3 r32/r48 instruments) moves in the predicted
direction.

**Step 3 — kernel-shape check (if a single Gaussian is retained).** Compute the
encircled-energy curve of the adopted kernel vs the FF single-scatter kernel
(same D2 machinery, full angular CDF instead of the moment) at r = 8.5, 32, 79 px
on the reference scene; report the three differences. No pass threshold is
proposed a priori — this step exists to put a number on the known shape error so
the retain/replace decision is made on a measurement, not on taste.

**Decision rule.** Adopt the candidate iff Step 1 passes all five clauses
(with the β̃(90°) reference declared) and Step 2 passes. On today's numbers,
HG-0.85 scores 0/5, HG-0.924 scores 3/5 (shape clauses fail), FF-Petzold
scores 5/5 — i.e. the measurements support changing the phase-function *form*,
not just the g value.

## 6. Verification of this report's own arithmetic, and open flags

Integrator cross-checks against published scalars (same code as §4):
HG-0.85 ⟨θ²⟩ = 0.3681 rad² — reproduces the v12 report's own quadrature exactly;
HG g = 0.9185 → B = 0.0183 — reproduces Mobley et al. (2002)'s stated value;
FF(1.10, 3.5835) → B = 0.0183, ⟨cos θ⟩ = 0.930 — matches the published Petzold-fit
B and the 0.924 average; FF Jonasz fits → ⟨cos θ⟩ = 0.9620 / 0.9488 vs published
0.962 / 0.948.

Flags (could not verify / known weak points):
1. Petzold's original 1972 tables and the verbatim average-particle table were
   not obtained; §1 pointwise values are Jonasz's measured averages (Petzold
   data included), and the 0.924 / per-type ⟨cos θ⟩ values are secondary
   citations as marked.
2. The Jonasz printed tables are geometric averages at sparse angles with
   scarce sub-10° data; integrating them directly is extrapolation-sensitive
   below 0.2° (my clear-table ⟨θ²⟩ moves ×0.42–0.53 across plausible
   extrapolations). The σ conclusions rest on the FF fits (which reproduce the
   published ⟨cos θ⟩ exactly) and on the fact that >97% of ⟨θ²⟩ comes from
   θ > 10° in every representation — not on the tables' small-angle points.
3. The sediment-preset envelope (§4) is an FF parameter-space statement, not a
   fit to a measured sediment VSF; no published freshwater sediment ⟨cos θ⟩ was
   found.
4. The alt. FF Petzold fit (n = 1.16, μ = 3.4586) comes from a secondary review
   chapter; the canonical fit used throughout is Mobley's (1.10, 3.5835).
5. Voss's m(τ) regression coefficients as printed in his Table 1 were not
   legible in the copy I read; the m ranges quoted are read from his figures
   and text (r² and ±0.3 uncertainty are from the text).
6. All measured data are marine/estuarine/lake; a pond-specific VSF dataset
   does not exist in the published record I could reach. The cross-waterbody
   shape conservatism is the published justification for transferring them.
