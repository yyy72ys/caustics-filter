# Physically-Grounded Water / Caustics Renderer — Environment Stack Research
## v5 (Part A: environment factors) + v8 (Part B: cosmic factors)

**Fixed scale for every verdict:** a shallow pond / river scene, depth 0.3–5 m, metres
across (order 1–30 m), scene timescales seconds to a day unless a factor only acts
over weeks–years, in which case that timescale is stated. The renderer already has:
wind-driven wave spectra, depth, turbidity, sun position, rain. "Real" below means
*adds information or visual change not already implied by those inputs, at this scale*.

**Verdict key (as assigned):** (a) real independent variable · (b) derived from other
variables (derivation given) · (c) honestly negligible / no effect at this scale.

**Confidence basis:** numbers below are from the web index (search-result full text of
the cited pages/papers), read 2026-11-10, cross-checked where load-bearing. No factor
required a live interactive page check — all are published physical quantities, not
prices/availability — so flags are `index` throughout, with primary papers preferred.

---

## Summary — the short answer for the implementer

**Add (Part A):** water temperature as a prognostic state (it moves viscosity by
~×1.8 per 20 °C, surface tension, and n by −10⁻⁴/K); humidity + air temperature
*only as the pair that drives evaporation / surface heat balance over hours*;
surface-film coverage (the single largest missing surface factor — a monomolecular
slick can suppress capillary waves by an order of magnitude); salinity *if* the scene
is brackish; and CDOM absorption as a spectral, not grey, turbidity component.

**Do not add (Part A):** water pressure / hydrostatic pressure as an optical or wave
variable (Δn ≤ ~10⁻⁵ over the full 5 m), atmospheric pressure as a surface-loading
or wave variable (uniform load, no slope; its fluctuating part *is* the wind you
already have), air-temperature mirage optics at metre scale.

**Add (Part B):** nothing as a separate knob. Every cosmic factor that is real at all
reaches the pond *through* quantities already in the stack (surface irradiance
spectrum, weather) and is 2–6 orders of magnitude below ordinary weather variation
at scene scale. Feed TSI / Kp / flare class in only as *annotation / provenance*
metadata if the artwork wants it, never as a physical forcing. The one genuinely
real cosmic-scale variable — Earth's orbital eccentricity (±3.4 % irradiance,
6.9 % peak-to-trough over a year) — is already inside "sun position + date".

---

# PART A — v5: Environment factors

## A1. Water pressure (hydrostatic)

**Physics.** Hydrostatic pressure is not a free input at all:

$$P(z)=P_{atm}+\int_0^z \rho g\,dz' \approx P_{atm}+\rho g z$$

At z = 5 m: ρgz = 1000 × 9.81 × 5 ≈ **49 kPa ≈ 0.48 atm gauge**, total ≈ 1.48 atm.
Everything pressure could do runs through compression. Water's bulk modulus is
K ≈ 2.2 GPa, compressibility β = 1/K ≈ 4.6–5.1×10⁻¹⁰ Pa⁻¹; even at 100 m depth the
fractional volume change is only 5.2×10⁻⁴ (0.05 %) [Compadre/Nexus bulk-modulus page].

- **Density:** Δρ/ρ = β·ρgz = 4.6×10⁻¹⁰ × 4.9×10⁴ ≈ **2.3×10⁻⁵ (0.0023 %)** at 5 m.
- **Refractive index:** via Lorentz–Lorenz, (n²−1)/(n²+2) ∝ ρ. With n = 1.333,
  d/dn of that function ≈ 0.56, so Δn ≈ 0.206 × (Δρ/ρ) / 0.56 ≈ **8×10⁻⁶ at 5 m**
  (order 10⁻⁵; consistent with the classic Waxler/NIST tables in which n only moves
  in the third decimal over *hundreds* of bar). Compare: 0.1 °C of temperature moves
  n by 10⁻⁵ (A2). Pressure over the whole pond is worth ≈ 0.08 °C.
- **Surface (i):** no path — surface pressure is P_atm, not ρgz; hydrostatic pressure
  does not enter the gravity-capillary dispersion ω² = gk + (σ/ρ)k³ at all.
- **Underwater light (ii):** Δn ~ 10⁻⁵ changes ray angles by ~10⁻⁵ rad and caustic
  positions by ~10⁻⁵ × path ≈ 0.05 mm over 5 m. Invisible.
- **Bottom (iii):** none optical. The only real bottom effect is mechanical load,
  which does not render.
- **Dissolved gas / bubbles:** this is the one place pressure is honest physics.
  Henry solubility rises with total pressure, and a bubble formed at 5 m expands by
  Boyle's law by a factor 1.48 on rising to the surface, and is *less* likely to
  nucleate at depth for the same gas tension. But in a 0.3–5 m pond there is no
  decompression regime (divers need ~atm-scale swings and supersaturated tissue);
  spontaneous bubble formation is controlled by photosynthetic O₂ supersaturation
  and nucleation sites, not by the 0–0.5 atm hydrostatic gradient. If the renderer
  ever grows a bubble-plume feature, bubble *size vs depth* (r ∝ P^(−1/3), +14 % in
  radius from 5 m to surface) is the correct, and only, pressure term to implement.

**Implementation:** none for optics/waves. Optionally r(z) for rising bubbles as above.
**Live data:** none needed — P(z) is computed from depth, which is already an input.

> **VERDICT: (b) derived — P(z) = P_atm + ρgz exactly — and (c) negligible in every
> optical/wave consequence at 0.3–5 m.** Reason: β ≈ 5×10⁻¹⁰ Pa⁻¹ caps Δρ/ρ at
> 2×10⁻⁵ and Δn at ~10⁻⁵ over the deepest scene; two orders below the temperature
> effect the renderer should already carry.

## A2. Water temperature — the biggest real factor in Part A

**Physics / numbers.**

- **Refractive index:** dn/dT ≈ **−1.0×10⁻⁴ K⁻¹** at ~20 °C (multiple sources agree
  on −0.0001/°C; n_D = 1.33434 just below 0 °C (its maximum), 1.3330 at 20 °C,
  1.31854 at 100 °C). So a realistic pond swing 5→25 °C moves n by **−0.002** —
  200× the full-depth pressure effect, ~0.15 % of n. That shifts critical angle
  (48.6°) by ~0.1° and Snell refraction measurably; a *correct* renderer should use
  n(T,λ,S), not a constant 1.333.
- **Density:** ρ = 999.84 kg/m³ at 0 °C, maximum **999.97 at 3.98 °C**, 998.20 at
  20 °C, 992.2 at 40 °C. Thermal expansion α ≈ 0 at 4 °C, ≈ 2.1×10⁻⁴ K⁻¹ at 20 °C.
  The 4 °C anomaly matters for stratification sign (below), not for optics directly.
- **Viscosity — the dominant wave path:** μ = **1.79 mPa·s at 0 °C, 1.00 at 20 °C,
  0.65 at 40 °C** (Kestin et al. correlation); ν ≈ μ/ρ ≈ 1.79→1.00→0.66 ×10⁻⁶ m²/s.
  Viscous amplitude damping of a free wave is γ = 2νk². At λ = 1 cm (k = 628 m⁻¹):
  γ ≈ 0.79 s⁻¹ at 20 °C vs **1.41 s⁻¹ at 0 °C** — capillary ripples die almost
  twice as fast in cold water. At λ = 10 cm γ ≈ 0.008 s⁻¹ (negligible either way).
  So temperature selectively reshapes the *capillary tail* of the spectrum — exactly
  the part that makes glitter and fine caustic texture.
- **Surface tension:** σ = 75.6 mN/m at 0 °C, **72.8 at 20 °C**, dσ/dT ≈ **−0.152
  mN m⁻¹ K⁻¹** (IAPWS correlation). In ω² = gk + σk³/ρ the capillary crossover
  λ_c = 2π√(σ/ρg) ≈ 1.7 cm and c_min ≈ 0.23 m/s shift by ~2 % per 10 °C — small,
  but free once T is carried.
- **Stratification / schlieren (ii):** *if* T varies with depth, n varies with depth:
  dn/dz = (dn/dT)(dT/dz), and a ray curves with radius R ≈ n/|dn/dz|. Well-mixed
  pond (dT/dz ≲ 0.1 K/m): R ≳ 10⁵ m — nothing. Sharp thermocline, 5 K over 0.1 m:
  dn/dz ≈ 5×10⁻³ m⁻¹, R ≈ 270 m; over a 5 m ray the bend is ~1° and displacement
  ~5 cm — borderline. A *turbulent mixing interface* (inflow plume, diver-visible
  "wrinkled glass" thermocline, where a few °C change over ~cm gives R ~ 10–30 m)
  is genuinely visible as shimmering schlieren — this is a documented scuba
  observation caused by "the altered refractive index of the cold or warm water
  column" [LibreTexts thermocline]. Implement only as a *localised* gradient /
  animated n-perturbation at a mixing interface, never as a bulk effect.
- **Bottom (iii):** only via the two paths above (n, stratification shimmer in front
  of the bottom), plus the slow thermal effect on CDOM/biology (out of scope).
- **Evaporation coupling:** T_water sets saturation vapour pressure at the surface,
  e_s(T) by Clausius–Clapeyron / Tetens (≈ +6–7 %/K at 20 °C) — the driver in A4.

**Is T independent?** In nature T_water is *prognostic* — the integral of the surface
heat balance (A3/A4). For a renderer it should be a **state variable**: either a user
input, or integrated from net radiation + sensible + latent flux. Do not derive it
instantaneously from air temperature; a 0.3 m layer has heat capacity
ρc_p h ≈ 1.3×10⁶ J m⁻² K⁻¹, so 200 W/m² changes it by only ~0.55 K/h — the lag *is*
the physics (ponds stay warm into the evening; that is visible in evaporation "steam"
and in n).

**Implementation:** carry T (ideally T(z) in 2–3 layers). n(λ,T,S) from a standard
formula (e.g. Quan & Fry / IAPWS-style polynomial); ν(T), σ(T) from 1-D lookups;
γ_damp(k) = 2ν(T)k² added to the spectral dissipation; optional schlieren term only
where a prescribed thermocline/plume interface exists.
**Live data:** no public feed knows *your* pond's T. In-situ probe, or model it from
Open-Meteo / ERA5 forcing (radiation, T_air, RH, wind) via the heat-balance ODE above.

> **VERDICT: (a) real independent (prognostic) variable — the strongest in Part A.**
> Reason: ν changes ×1.8 from 20→0 °C (capillary damping ×~1.8), n changes
> 2×10⁻³ over a normal seasonal swing, σ changes 0.2 %/K. Stratification optics are
> (c) negligible *except* at sharp (≳ K/cm-scale) mixing interfaces, where schlieren
> are real and documented.

## A3. Air temperature

**Paths in.** T_air itself touches nothing optical in the water. It enters through:

1. **Sensible heat flux** (surface heat balance): Q_H = ρ_a c_p C_H U (T_s − T_a),
   ρ_a c_p ≈ 1200 J m⁻³ K⁻¹, C_H ≈ 1.3×10⁻³. At U = 5 m/s and ΔT = 5 K:
   Q_H ≈ **39 W/m²** — same order as the latent flux (A4), ~5–20 % of midday solar.
   Sign flips when air is warmer than water (common on spring afternoons over a cold
   pond — then air *heats* the water).
2. **Evaporation**, jointly with RH, via q_a = RH · q_s(T_a) (A4).
3. **Near-surface air refraction / mirage — the question asked.** Air:
   n_a − 1 ≈ 2.9×10⁻⁴ (P/1013 hPa)(273 K/T); dn_a/dT ≈ −1×10⁻⁶ K⁻¹. Ray curvature
   needs an extreme gradient: an inferior mirage requires dT/dz < **−3.5 °C per
   10 cm** (≈ −35 K/m), a condition met over sun-baked asphalt, and over water only
   in a mm–cm skin when warm water heats cold air from below [Stull, *Practical
   Meteorology*, LibreTexts]. Superior mirage / looming over cold water is real and
   photographed — at **6 km** object distance, air 27–30 °C over water 12–14 °C
   [Galicia mirage photographs, ResearchGate]. Scaling is the killer: angular
   deviation accumulates ∝ path length through the gradient layer. At 6 km it is
   arcminutes-to-degrees; across a **5–30 m pond it is ~10⁻³ of that** — sub-arcsecond
   to a few arcseconds, far below a pixel and below the surface-wave slope noise
   (wave slopes are ~0.1 rad). What *can* appear at pond scale is generic heat
   shimmer (turbulent C_n² fluctuations) in a long, grazing *reflection* path on a
   calm, strongly ΔT day — an aesthetic noise term, not mirage geometry.
4. **Air density** ρ_a = P/(R_d T_v) enters wind stress τ = ρ_a C_d U² — a 10 K
   swing changes ρ_a, hence stress at fixed U, by ~3.5 %. Second-order vs the
   U² and C_d uncertainties already in the wind model; fold in only if T_a is
   already being fed for (1).

- **Surface (i):** only via (4) and via T_s evolution. **Underwater (ii):** none
  direct; via T_s(t) → n(T). **Bottom (iii):** none.

**Implementation:** input to the heat-balance ODE for T_s (A2) and to evaporation
(A4); do *not* build an above-water graded-index ray marcher for a pond scene.
Optional: shimmer noise amplitude ∝ |T_s − T_a| / U on grazing reflected rays.
**Live data:** Open-Meteo forecast/archive API (temperature_2m, hourly, no key) or
any met station / ERA5.

> **VERDICT: (b) derived / indirect — a real forcing, but only through the surface
> heat and moisture fluxes (formulas above); it is not an optical variable.**
> Mirage/looming over water at scene scale: **(c) negligible** — requires
> −35 K/m gradients and/or kilometre paths; at metres the deviation is ≲10⁻³ of the
> documented sea cases.

## A4. Humidity

**Physics.** Humidity's one real job is evaporation. Bulk aerodynamic (Dalton) form,
as used in lake models (DYRESM):

$$E=\rho_a C_E U\,[\,q_s(T_s)-q_a\,],\qquad Q_E=L_vE,\qquad q_a = RH\cdot q_s(T_a)$$

with C_E ≈ 1.3×10⁻³, L_v ≈ 2.45×10⁶ J/kg, q_s(20 °C) ≈ 0.0147 kg/kg.
Worked example: U = 5 m/s, T_s = T_a = 20 °C, RH = 50 % → Δq ≈ 0.0074 →
E ≈ 5.8×10⁻⁵ kg m⁻² s⁻¹ ≈ **5 mm/day**, Q_E ≈ **142 W/m²**. That is a first-order
term in the pond's energy budget — comparable to the daily-mean solar input.
Cooling rate of a 0.3 m surface layer: Q_E/(ρc_p h) ≈ **0.4 K/h** (0.04 K/h if mixed
over 3 m). At RH = 100 % with T_a = T_s, E → 0 exactly; with warm humid air over
cold water q_a > q_s(T_s) and the flux reverses — condensation, i.e. the physical
origin of "steam fog" / mist over a warm pond on a cold morning (the reverse case).

- **Surface (i):** no direct mechanical/optical effect. Indirectly, evaporative
  cooling creates a cool skin (~0.1–0.5 K cooler, mm thick) which very slightly
  raises σ and ν *in the skin* — a <1 % correction to capillary damping; note it,
  don't model it separately from T_skin if you model a skin at all.
- **Underwater (ii):** none. Water vapour changes air's index by only a few ×10⁻⁶
  between dry and saturated at fixed P,T (and in the "wrong" direction — moist air
  has *lower* n) — utterly below wave-slope noise at the interface.
- **Bottom (iii):** none except through multi-day water level / T evolution.
- **Timescale honesty:** over a rendered *minute*, evaporation removes ~3.5 mg/m²
  and cools by ~10⁻⁴ K. Humidity matters to the artwork only if it simulates
  hours–days, or if it wants the *visible* by-product: evaporation fog when
  e_s(T_s) ≫ e_a (cold air over warm water) — that is a particle/fog effect, and it
  is dramatic and real at pond scale in autumn mornings.

**Implementation:** prognostic only: E and Q_E in the T_s / level ODEs; a fog switch
when T_s − T_dewpoint,a ≳ ~5–10 K with low wind (empirical steam-fog regime).
**Live data:** same Open-Meteo fields (relative_humidity_2m, dew_point_2m) / ERA5.

> **VERDICT: (a) real independent variable — but slow.** Its effect is exactly the
> bulk formula above: up to ~140 W/m² latent cooling at RH 50 %, zero at saturation.
> Direct optical role: **(c) negligible** (Δn_air ~ 10⁻⁶). On sub-hour scene time
> it can be frozen as a constant; on day time it must be carried or T_s will drift
> by degrees.

## A5. Atmospheric pressure

Separate the *mean* pressure from its *fluctuations* — they have opposite verdicts.

**(i) Mean / synoptic pressure (±30 hPa weather swings).**

- **Direct loading:** a uniform ΔP over a pond produces **no surface slope and no
  current** — Pascal: the same ΔP acts on every point of the surface and (through
  the water) on the bottom. The famous inverse-barometer response,
  Δη = −ΔP/(ρg) ≈ **−1.02 cm per hPa** (≈ ±30 cm for ±30 hPa) [Padman et al.,
  J. Glaciology, after Gill 1982], is an *ocean* equilibrium against a *spatially
  varying* pressure field with room to move water horizontally. Over a metres-wide
  enclosed pond the pressure gradient across the scene is ~ (30 hPa / 1000 km) ×
  10 m ≈ 3×10⁻⁴ hPa — a tilt of ~3 µm. If the pond is a reach of a river/lake
  system connected to a large water body, the level change arrives as a *boundary
  water-level input*, not as local physics.
- **Air index:** n_a − 1 ∝ P/T, so ±30 hPa (±3 %) moves n_a by ±9×10⁻⁶ — refraction
  at the surface changes by ~10⁻⁵ relative. Negligible (cf. A1's Δn).
- **Water compression:** ΔP = 3 kPa → Δρ/ρ = βΔP ≈ 1.4×10⁻⁶, Δn ≈ 5×10⁻⁷. Nothing.
- **Henry's law / bubbles / foam:** equilibrium dissolved-gas concentration
  c = k_H · p_gas scales with partial pressure, so a ±30 hPa swing is a **±3 %
  change in gas solubility** [dissolved-oxygen/Henry references]. Bubble formation
  needs supersaturation *ratios* of order 1.5–several (or a pressure *drop* of a
  substantial fraction of an atmosphere, as in a waterfall plunge or a propeller) —
  a 3 % solubility drift crosses no nucleation threshold, creates no foam, and
  changes nothing visible. Foam in ponds is surfactant-stabilised (A6), not
  pressure-made.

**(ii) Fluctuating pressure — wave generation.** This is real, classical, and
*already in the renderer*. Phillips (1957): turbulent pressure fluctuations advected
over a flat surface resonantly force the (k,ω) components matching the free-wave
dispersion relation, growing wave amplitude *linearly* in t from rest. Miles (1957):
once waves exist, the critical layer (where U(z) = c) feeds energy exponentially.
Together the Miles–Phillips mechanism [Wikipedia/encyclopedia summary of
Phillips 1957, Miles 1957]. Crucial implementation point: the pressure spectrum
Π(k,ω) here is a property of the **wind's turbulence**, parameterised in practice
by U₁₀ / u* in every modern wave-growth source term (which the existing wind-driven
spectrum already uses). Adding a separate term driven by *barometric* pressure
would double-count, and synoptic ΔP itself contributes ~nothing to Π at wave
scales (its energy is at 100–1000 km, days).

- **Surface (i):** mean: none; fluctuations: the wave source itself (already had).
  **Underwater (ii):** none beyond waves. **Bottom (iii):** none.

**Implementation:** do nothing with P_atm except — if a connected large water body
is ever simulated — Δη_IB as a boundary level offset. Keep Phillips/Miles growth
inside the wind source term, parameterised by wind, as now.
**Live data:** Open-Meteo surface_pressure / pressure_msl if the annotation or the
boundary-level case is wanted.

> **VERDICT: mean pressure — (c) honestly negligible at pond scale** (uniform load,
> no gradient; ±3 % gas-solubility drift; Δn ~ 10⁻⁶). **Fluctuating pressure —
> (b) derived from the wind field** (Phillips resonance + Miles critical layer);
> it is a real wave-generation mechanism but must not be added twice.

## A6. Further factors a specialist would insist on

### A6.1 Surface films / slicks — the largest omitted surface factor. REAL.
A monomolecular layer (10⁻⁹–10⁻⁸ m) of natural surfactant — from phytoplankton
exudates, decaying organics, pollen, or a trace of oil — changes the surface
boundary condition, not just σ. Stretching a film-covered surface locally dilutes
the surfactant, raising local σ; the resulting Marangoni (Gibbs elasticity) stress
resists surface divergence and dissipates wave energy in the viscous boundary
layer. Measured consequences: films "dampen small-scale waves… waves with lengths
of centimetres to decimetres" strongly enough to be visible in SAR imagery
[Alpers & Hühnerfuss]; in a tank with a fish-oil film "the wave field does not grow
from the rest condition" under wind alone, and with paddle waves the standard
deviation of surface elevation "shrinks by one order of magnitude" [Copernicus
Ocean Science discussion paper os-2018-111]. Damping is *resonant* in k
(Marangoni resonance, peaking in the cm range) and grows ~k^(3/2) outside it —
precisely the capillary range that produces sun glitter and fine caustics.
Films also lower mean σ (clean 72.8 → commonly 50–65 mN/m under a slick) and cut
wind input and gas transfer. A pond is *more* prone to films than the sea: low
fetch, low wind, high biological productivity, and films persist until wind
≳ ~5–7 m/s breaks them up (empirical, variable — see Could-not-verify).

**Implementation:** a coverage fraction F ∈ [0,1] (user control or a slow state
advected by wind, accumulating at downwind shores). Multiply the spectral
dissipation by a Marangoni factor y_M(k; |E|) with film elasticity |E| ~ 10–50
mN/m as the tuning parameter, and scale the wind-input source by (1 − aF) in the
resonant band. Even a crude version — extra damping β_film(k) = F · β₀ ·
exp[−(ln(k/k_res))²/2s²], k_res ≈ 2π/0.03 m — will read correctly on screen:
glass-smooth patches inside a rippled pond.
**Live data:** none public for a specific pond; F is an artistic/site parameter.

> **VERDICT: (a) real independent variable.** Reason: order-of-magnitude suppression
> of cm-wave amplitude at monolayer coverage — larger than any other Part A surface
> effect except the wind itself.

### A6.2 Salinity / dissolved solids. CONDITIONALLY REAL.
For seawater S = 35: n = **1.33940 vs 1.33300 fresh** (Δn = +0.0064, i.e.
dn/dS ≈ **+1.8×10⁻⁴ per g/kg**), density ×1.0264 (∂ρ/∂S ≈ +0.8 kg m⁻³ per g/kg)
[Holmes-Farley, Reefkeeping, standard seawater values]. Scaling down: a freshwater
pond at S ≲ 0.5 g/kg differs from pure water by Δn ≲ 10⁻⁴, Δρ ≲ 0.4 kg/m³ —
borderline invisible. A brackish estuary reach, salt-marsh pond, road-salt or
evaporation-concentrated pool at S = 5–35 is a different scene optically
(critical angle shifts ~0.3° at S = 35) and dynamically (haline stratification is
much stronger per unit than thermal: 1 g/kg ≈ the density effect of ~5 K).
Also note seawater σ is ~+1 mN/m higher and ν ~5 % higher — second-order.

**Implementation:** if the scene can be brackish, add S as a constant (or 2-layer)
parameter feeding n(T,λ,S) and ρ(S,T); skip entirely for a declared freshwater pond.

> **VERDICT: (a) real independent variable in brackish/saline scenes (Δn up to
> 6×10⁻³, the largest index change in this whole report); (c) negligible for a
> freshwater pond with S < ~0.5 g/kg.** The verdict is scene-dependent, and the
> renderer should ask, not assume.

### A6.3 CDOM (coloured dissolved organic matter). REAL — and it is *not* turbidity.
The existing "turbidity" is (presumably) particulate scattering. CDOM is *dissolved*
absorption, strongly blue-weighted: a_CDOM(λ) = a_CDOM(440)·exp[−S_slope(λ−440)],
S_slope ≈ 0.014–0.019 nm⁻¹. Measured inland/coastal range a_CDOM(440) ≈
**0.1–8.5 m⁻¹, up to ~20 m⁻¹** in humic freshwaters [Remote Sensing of Environment
2014 assessment] — against pure water's own absorption of only ~0.01–0.02 m⁻¹ at
440 nm (but 0.3–0.6 m⁻¹ at 650 nm, which is why clear water is blue). So in a
peaty pond CDOM exceeds water absorption in the blue by **2–3 orders of magnitude**,
turns transmitted light and caustics yellow→tea-brown, and shortens blue
penetration to ~1/a ≈ 0.1–1 m while red still dies by water itself. No grey
turbidity knob can imitate that spectrum. Chlorophyll absorption (peaks 440/675 nm)
is the biological sibling — add only if algal blooms are in scope.

**Implementation:** split attenuation into absorption a(λ) = a_water(λ) +
a_CDOM(440)e^(−S(λ−440)) + a_chl(λ) and scattering b(λ) (existing turbidity);
Beer–Lambert per spectral band in the ray transport — the renderer is already
spectral, so this is one extra exponential per band.
**Live data:** none for a pond; Sentinel-2/3 ocean-colour retrievals exist only for
large lakes/coasts. Site parameter / water sample.

> **VERDICT: (a) real independent variable.** Reason: a(440) spanning 0.1–20 m⁻¹
> vs 0.015 m⁻¹ for pure water is the dominant control on underwater colour in
> most real ponds and rivers.

### A6.4 Shorter verdicts in this class
- **Ice cover / freezing:** not a perturbation — a phase change that replaces the
  surface boundary condition entirely (no waves, TIR from below, σ undefined).
  (a) real *if* winter scenes are in scope; otherwise out of scope by declaration.
  Note the 4 °C density maximum (A2) is what makes ice form on top.
- **Dissolved O₂ / gases (non-bubble):** no refractive or mechanical effect at
  mg/L concentrations beyond the Henry discussion in A5. **(c).**
- **pH / hardness / nutrients per se:** no direct optics; they act only through
  biology → CDOM/chlorophyll/films above. **(c) direct, (b) via A6.1/A6.3.**
- **Suspended sediment pulses / bottom resuspension:** real, but already inside
  the existing turbidity input; a specialist would ask that turbidity be allowed
  to vary with rain/runoff and current, which the renderer already models as rain.
- **Groundwater inflow / spring plume:** real locally as a T/S/CDOM *inhomogeneity*
  — visible exactly by the schlieren path of A2. Model as a plume boundary
  condition on A2/A6.2/A6.3, not as a new variable.

---

# PART B — v8: Cosmic-level factors

Framing number for all of Part B: ordinary weather at the pond surface swings
shortwave irradiance by **hundreds of W/m² in seconds** (a cloud edge), air
temperature by ~10 K/day, pressure by ±30 hPa. A cosmic factor earns a knob only
if it can beat ~1 % of that, through a stated path, on a stated timescale.
None does.

## B1. Solar activity — TSI and spectral (UV) variation

**Numbers.** TSI ≈ **1361 W/m²** at 1 AU; over the 11-year cycle it varies by
**≈ 0.1 %, ≈ 1.3 W/m² peak-to-trough**, higher at solar maximum (faculae beat
sunspots) [Frontiers review of TSI composites; solar-constant references].
In global-mean forcing terms that is 1.3/4 ≈ 0.33 W/m² at TOA, ~0.24 W/m²
absorbed — and the IPCC AR6 net solar forcing change *since 1750* is only
**+0.05 W/m² vs +2.72 W/m² for greenhouse gases** (secondary synthesis of AR6;
flagged in Could-not-verify). The observed 11-year global surface-temperature
imprint is ~0.07 °C. Spectrally the variation is far larger but in the wrong
place: **a few % at 200–300 nm, 10–20 % at 150–200 nm, 50–100 % at Lyman-α
(121.6 nm)** [Krivova et al., A&A 2006]; the UV share of TSI swings 6–8 %.
But λ < 200 nm is absorbed in the thermosphere/mesosphere and 200–300 nm
almost entirely by stratospheric ozone — *none of the strongly varying part
reaches the pond*. Surface UV (315–400 nm) varies by of order 1 % over a cycle,
visible light by ~0.1 %. Compare the variation the renderer *must* already
handle: Earth's orbital eccentricity makes TOA irradiance vary **6.9 %
peak-to-trough over one year** (perihelion ~1408 → aphelion ~1317 W/m²) — 70×
the solar-cycle signal, and it is already inside sun position + date.

- **Paths:** (i) surface — only via total irradiance heating (already an input if
  irradiance is fed, or via sun geometry). (ii) underwater — the varying UV does
  not arrive; a real *slow* path exists: UV drives CDOM photobleaching and
  phytoplankton biology over weeks–seasons, so solar-cycle UV could in principle
  modulate A6.3 on decadal time — unmeasured at pond scale, and swamped by
  catchment/seasonal CDOM changes. (iii) bottom — none.
**Implementation:** none. If a multi-year artwork mode exists, scale the *input
irradiance spectrum* by a TSI/SSI dataset — i.e. change the input data, not the model.
**Live data (annotation-grade):** LASP LISIRD (TSI + spectral SSI from TSIS-1),
NOAA/NASA SORCE heritage data, SILSO sunspot number.

> **VERDICT: (b) derived — real at the top of the atmosphere (0.1 % TSI, %–100 %
> in UV bands that never reach the surface), but at the pond it is entirely
> contained in the irradiance input the renderer already takes; (c) negligible
> as a separate variable.** Reason: 1.3 W/m² TOA vs hundreds of W/m² of weather
> and 91 W/m² of annual orbital variation.

## B2. Solar flares and coronal mass ejections

**What actually happens.** A flare is minutes of enhanced short-λ emission:
averaged over X-class flares, irradiance at peak rises **+680 % at 0–14 nm,
+65 % at 14–25 nm, +120 % at 25–105 nm, +35 % at 122–190 nm** [FISM modelling,
Huang et al., JASTP 2014] — huge percentages of bands that carry a negligible
share of TSI and are absorbed at 100–400 km altitude, heating the thermosphere
by tens of K *up there* and enhancing ionospheric electron density (sudden
ionospheric disturbance, HF radio blackout on the sunlit side). The *total*
irradiance excursion even in the largest recorded flares is at/below the
~0.1 % level for minutes (see Could-not-verify for the exact X28 figure).
A CME is plasma, arriving 1–3 days later *if* Earth-directed: it drives a
geomagnetic storm (B3), aurora, and geomagnetically induced currents (below).
Energetic-particle (SEP) events raise ionisation down to stratospheric levels
at polar latitudes (ozone chemistry over weeks) — a polar stratosphere effect,
not a troposphere/pond effect.

**What reaches the surface, itemised:**
- Extra visible/near-IR light or heat at the pond: **none measurable.**
- Ionisation change in the *lower* troposphere: flare X-rays are absorbed
  ≳ 60–90 km; no quantified surface-level ionisation change.
- Geomagnetically induced currents: real and sometimes destructive — in
  **long grounded conductors**: GIC of tens to hundreds of amperes recorded in
  power grids/pipelines; the 13 March 1989 storm collapsed the Hydro-Québec grid
  [GIC references]. The geoelectric field is ~volts per *kilometre*; across a
  10 m pond that is ~10–30 mV distributed through water/earth, driving no
  chemistry, no bubbles, no visible effect. A pond is not a transmission line.
- Indirect human path (grid outage changes artificial light at night) is real
  sociology, not scene physics; noted only to be dismissed.

**Implementation:** none. **Live data:** NOAA SWPC GOES X-ray flux / flare list,
R/S/G scales, NASA DONKI (CME) — legitimate as *display metadata* only.

> **VERDICT: (c) honestly negligible for the pond, on every timescale.**
> Reason: 100 % of the enhanced radiation is deposited above ~60 km; the only
> surface manifestation (GIC) scales with conductor length in kilometres and is
> ~10⁻² V across the whole scene. No quantified path to surface weather or
> water exists in the flare/CME literature.

## B3. Geomagnetic activity (Kp, storms)

**The claim to assess:** that geomagnetic activity couples to tropospheric
weather — classically a lagged correlation between the aa/Kp index and the
winter North Atlantic Oscillation / Arctic Oscillation, and the related
"top-down" solar-UV → stratospheric ozone → polar vortex → NAO mechanism.

**State of evidence, honestly:**
- The *top-down UV mechanism* (Gray et al. and successors) is physically
  constructed and appears in some model studies: solar-max UV heats the upper
  stratosphere, alters the polar vortex, and a small NAO-like surface signal
  (~fraction of an NAO standard deviation, winter only, high-latitude North
  Atlantic) emerges in composites. It is a **UV/irradiance** mechanism (B1),
  not a Kp mechanism — Kp is at best a correlated proxy.
- Direct **aa/Kp–NAO correlations** have been reported repeatedly (Bucha &
  Bucha; Boberg & Lundstedt; Palamara & Bryant; Maliniemi et al.) and reviewed
  in the geomagnetic-stratosphere literature [ScienceDirect, *Impact of
  geomagnetic activity on stratosphere and upper troposphere*], but the
  relationship is documented as **non-stationary** — e.g. "weakly negative" in
  some epochs, positive in others — and day-to-day autoregressive studies find
  only small, conditional (QBO-phase-dependent) effects. A correlation that
  changes sign between epochs, selected post hoc over seasons/regions/lags,
  is not an established forcing; IPCC AR6 does not attribute any robust
  tropospheric weather effect to geomagnetic activity.
- Energetics bound it independently (B6): even extreme magnetospheric input
  is ≲ 0.03 W/m² global-mean, deposited at ≳ 90 km, mostly over the polar caps.

**At pond scale:** even taking the most favourable published NAO claim at face
value, it is a *statistical shift in the odds of winter circulation regimes
over an ocean basin over weeks* — i.e. it would enter, if at all, as a tiny
change in the *weather inputs* (wind, pressure, temperature), which the
renderer takes as given. There is no local Kp → wave/water path: the pond
does not respond to a 100 nT field change (Lorentz force on water ions at
thermal/flow speeds is ~10⁻¹⁵ N scale; magnetohydrodynamic effects need
laboratory tesla fields and conducting melts, not 5×10⁻⁵ T ± 0.2 %).

**Implementation:** none. **Live data:** NOAA SWPC planetary K-index JSON
(3-hourly, plus 1-minute estimates), GFZ Potsdam definitive Kp — annotation only.

> **VERDICT: (c) negligible / no defensible effect at scene scale.** Reason:
> the proposed tropospheric link is (i) contested and non-stationary in the
> literature, (ii) at best a small statistical modulation of the *weather
> inputs already supplied*, and (iii) energetically capped ~4 orders below
> tropospheric energy fluxes (B6). Kp is a magnetosphere index, not a
> water variable.

## B4. Galactic cosmic rays and cloud nucleation (Svensmark / CERN CLOUD)

**The hypothesis:** GCR ionisation → aerosol nucleation → more cloud
condensation nuclei (CCN) → more/brighter low cloud → cooling; solar activity
modulates GCRs, closing a Sun–climate amplifier.

**What is measured, link by link:**
1. **GCR modulation is real and large:** neutron-monitor count rates vary
   ~20 % over a solar cycle (factor ~2 at ~100 MeV energies); Forbush decreases
   during CME passages are typically < 3 %, rarely > 10 %, with 3–20 %
   depressions reported for major events. Surface-level ionisation from GCRs
   is a few ion pairs cm⁻³ s⁻¹ (order 2–10, altitude/latitude dependent).
2. **Ion-induced nucleation is real in the chamber:** CERN CLOUD (Kirkby et al.,
   *Nature* 2011) found ionisation enhances nucleation of H₂SO₄–H₂O clusters
   by a factor **2–10 or more**, mainly in the cold mid/upper troposphere;
   ammonia raises nucleation 100–1000× by comparison [CERN/Physics Today
   summaries]. The nucleated particles are ~1–2 nm; CCN need ~50–100 nm,
   so survival/growth is the bottleneck.
3. **The amplifier fails at the global step:** Dunne et al. (*Science* 2016),
   building a global aerosol model directly from a decade of CLOUD rates,
   concluded variations in cosmic-ray intensity "do not appreciably affect
   climate through nucleation in the present-day atmosphere" [White Rose /
   *Science* abstract]; the CLOUD team's own summary quantifies it: ions are
   involved in forming ~⅓ of new particles, yet **CCN concentration changes by
   only 0.1 % over the solar cycle** [Goethe-Universität Frankfurt summary of
   Dunne et al. 2016]. Satellite GCR–cloud-cover correlations from the original
   Svensmark claims have not survived reanalysis robustly (data artefacts,
   El Niño confounding), and IPCC assessments assign no significant
   GCR–cloud forcing.
4. **At the pond:** the hypothesised chain ends in *cloud cover* — which, in
   this renderer, is not even a variable; irradiance/sun is an input. A 0.1 %
   CCN change → ≪0.1 % cloud change → ≪0.1 W/m² at the surface, against
   hundreds of W/m² of actual cloud weather. Direct GCR ionisation of pond
   water: dose rates ~nGy/h scale; no optical, chemical or bubble consequence.

**Implementation:** none. **Live data:** Neutron Monitor DataBase (NMDB) count
rates, Oulu cosmic-ray station — annotation only, if ever.

> **VERDICT: (c) negligible — a real microphysical effect with a measured,
> broken amplifier.** Reason: CLOUD confirms step 1 (×2–10 nucleation) and
> simultaneously measures the full chain to be a **0.1 % CCN change per solar
> cycle**, i.e. ~3 orders below anything visible, before it even reaches the
> renderer's already-given sky.

## B5. Lunar / tidal and Earth-rotation factors at pond scale

- **Direct lunar gravity on the pond (ocean-tide mechanism):** ocean tides need
  a basin comparable to the tidal wavelength and a connection to the ocean.
  For an enclosed pond, the lunar tidal acceleration is a_t ≈ 1.1×10⁻⁶ m/s²
  (horizontal, max). It tilts the equipotential — and with it the water surface,
  which simply follows the equipotential — by a_t/g ≈ **1.1×10⁻⁷ rad**, i.e. a
  **~1 µm height difference across 10 m**, and, critically, the *ground and
  the observer tilt with it*: there is no relative signal to render. The
  associated change in effective g is Δg/g ~ 10⁻⁷, shifting wave dispersion by
  5×10⁻⁸ relative. Nothing.
- **Solid Earth tide:** the crust itself rises/falls — M2 vertical amplitude
  can reach **~55 cm at the equator**, typically 20–30 cm twice daily
  [Earth tide references, after Agnew 2007]. Spectacular geodesy, zero scene
  content: pond basin, water, banks and camera move *together* over a
  ~10⁴ km wavelength deformation. Only a surveyor comparing against a
  satellite-fixed frame can see it.
- **Atmospheric lunar tide:** real and measured — surface-pressure amplitude
  **~0.06 hPa in the tropics, ~0.02 hPa in mid-latitudes** (AMS Glossary of
  Meteorology; other compilations: ~5–10 Pa max at low latitudes, equatorial
  wind amplitude ~0.03 m/s). That is 300–1500× smaller than a weather swing
  and, per A5, uniform pressure does nothing to a pond anyway; the tidal wind
  (3 cm/s, at altitude-weighted levels) is 2 orders below the lightest air
  movement the wave spectrum can feel, and it is *in* any measured wind input.
- **The one real exception — connection:** if the "pond" is a tidal river reach
  or an estuary pond connected to the sea, lunar influence arrives *massively*
  — as metres of water-level and current variation. But that is a **boundary
  condition** (downstream ocean tide propagated up the channel), correctly
  modelled as prescribed level/current forcing, not as a local lunar force.
  The renderer should ask whether the scene is tidally connected; if yes, feed
  a tide-gauge / harmonic level series as an input.
- **Earth rotation / Coriolis:** f = 2Ω sinφ ≈ **1.0×10⁻⁴ s⁻¹** at 45°. For pond
  flow U ~ 0.1 m/s over L ~ 10 m, Rossby number Ro = U/(fL) ≈ **100**
  (a draining bath: Ro ~ 3200). Coriolis acceleration fU ≈ 10⁻⁵ m/s² ≈ 10⁻⁶ g;
  the inertial radius U/f ≈ 1 km ≫ the pond. Rotation organises weather and
  ocean currents (Ro ≲ 1) — it does not touch metre-scale water motion, wave
  dispersion (Ω/ω ~ 10⁻⁴ for a 1 Hz wave), or light. Length-of-day variations
  (ms) and polar motion are smaller still by orders.

**Implementation:** none for an enclosed pond; harmonic tide level/current
input *iff* tidally connected. **Live data:** tide predictions/gauges (e.g.
national hydrographic services) for the connected case only.

> **VERDICT: (c) honestly negligible for an enclosed pond** — differential lunar
> effect ~1 µm across the scene, atmospheric lunar tide ≤0.06 hPa, Coriolis
> Ro ~ 100. **Exception: (b) derived/real as a boundary level input if, and
> only if, the water body is tidally connected to the ocean.**

## B6. Anything else at this level with a defensible measured influence

### B6.1 Solar wind / magnetosphere energy budget — the master bound for B2–B3
At 1 AU the solar-wind kinetic and thermal energy fluxes are each
~**5×10⁻⁴ W/m²**, and the interplanetary magnetic field's ~5×10⁻⁶ W/m²,
against **1365 W/m²** of radiation — a factor ~10⁶ smaller [Haigh et al. /
arXiv solar-forcing review]. Only 1–10 % of even that couples into the
magnetosphere. In absolute terms: magnetospheric input is ~0.1–0.2 TW quiet,
up to ~2 TW in active times; for the extreme October 2003 storms, Cluster-based
estimates put global input at **~17 TW** (empirical proxies said 37 TW), of
which ~30 % dissipated as ionospheric Joule heating [Rosenqvist et al. 2006,
ESA/Cluster]. Spread over Earth, 17 TW is **0.033 W/m²**, deposited at
≳ 90 km altitude, concentrated over the polar caps — vs ~240 W/m² absorbed
solar and ~80 W/m² of global latent-heat flux in the troposphere. There is
simply no energy available for space weather to *drive* tropospheric weather
directly; any influence would have to be a trigger/amplifier, which is exactly
what B3/B4 assessed and found unproven/too small.

> **VERDICT: (c).** This budget is the quantitative reason B2 and B3 cannot be
> rescued by "unknown mechanisms" of any direct kind.

### B6.2 Stratospheric aerosol events (volcanic, and meteoric smoke) — real, but derived
A major eruption (Pinatubo 1991) raised stratospheric aerosol optical depth to
~0.1–0.15, cut surface direct sunlight by ~10 %+ while increasing diffuse light,
and imposed ~ −3 W/m² global forcing for ~2 years. At a pond this *is* visible:
a milkier sky, softer shadows, redder sunsets, slightly dimmer/cooler caustics.
But note the path: it changes the **irradiance spectrum and direct/diffuse
ratio at the surface** — inputs the renderer already takes (or should: a single
aerosol optical depth / turbidity-of-air parameter fed by data). Ordinary
tropospheric aerosol (dust, smoke, pollution) does the same thing more often
and more strongly.

> **VERDICT: (b) derived — implement as surface irradiance spectrum + direct/diffuse
> split (one air-optical-depth parameter), fed by e.g. Copernicus CAMS aerosol
> data or a local sun-photometer; do not implement "volcanoes" in the renderer.**

### B6.3 Checked and dismissed, with the number that dismisses each
- **Orbital eccentricity / precession / Milankovitch:** real, and the largest
  "cosmic" signal of all (6.9 % annual irradiance swing, B1) — **(b)**, already
  inside sun position + date. Millennial terms are irrelevant to a scene.
- **Geomagnetic *field itself* (main field, secular variation):** ~25–65 µT,
  changing ~0.01 %/yr. No coupling to neutral water (see B3). **(c).**
- **Schumann resonances / global electric circuit:** fair-weather field
  ~100–150 V/m vertical; the energy density and any force on the water surface
  (electrostatic pressure ε₀E²/2 ~ 10⁻⁷ Pa) are ~9 orders below capillary
  pressures (σk ~ 10–100 Pa). Claims linking the circuit to cloud microphysics
  remain unproven and, per B6.1, energy-starved. **(c).**
- **Geothermal heat flux:** global mean ~0.09 W/m² (≈47 TW / Earth's area;
  continental ~0.06 W/m²) vs 100s of W/m² solar. **(c)** for an open pond —
  with one honest footnote: in *ice-covered* winter lakes, geothermal + stored
  heat helps maintain the 4 °C bottom layer; that is a boundary condition on
  A2 in a winter model, not a cosmic factor.
- **Neutrinos, dark matter wind, planetary tides, "planetary alignments":**
  solar neutrino energy flux is ~2–3 % of photon flux *in number-energy terms
  at the Sun's surface accounting* but the interaction cross-section is so
  small that the deposited fraction in a pond is ~10⁻¹⁵ of it; planetary tidal
  accelerations are ~10⁻⁵ of the Moon's (which was already 1 µm, B5). **(c)**,
  recorded here so the negative is on the record as instructed.

---

## Cross-cutting verdict table

| # | Factor | Verdict | One-line reason (quantified) |
|---|--------|---------|------------------------------|
| A1 | Water (hydrostatic) pressure | (b) derived, (c) in effect | P=P_atm+ρgz; Δρ/ρ ≤2×10⁻⁵, Δn ~10⁻⁵ over 5 m (β≈5×10⁻¹⁰/Pa) |
| A2 | Water temperature | **(a) real** | ν ×1.8 per 20 °C (capillary damping), dn/dT −10⁻⁴/K, dσ/dT −0.15 mN/m/K |
| A2′ | Thermal stratification optics | (c), except sharp interfaces | R ≳10⁵ m well-mixed; schlieren real only at ~K/cm mixing interfaces |
| A3 | Air temperature | (b) derived | Acts only via sensible flux (~40 W/m² at ΔT 5 K, U 5 m/s) + evaporation |
| A3′ | Mirage over water | (c) | Needs −35 K/m / km paths; ≲10⁻³ of sea cases at metres |
| A4 | Humidity | **(a) real, slow** | Bulk evaporation: 0 at RH 100 % → ~140 W/m² latent at RH 50 %; direct optics (c) |
| A5 | Atmospheric pressure (mean) | (c) | Uniform load, no slope; IB 1 cm/hPa needs an ocean; ±30 hPa = ±3 % gas solubility, Δn ~10⁻⁶ |
| A5′ | Pressure fluctuations (Phillips/Miles) | (b) derived from wind | The wave-growth mechanism itself — already in the wind source term; adding barometric P double-counts |
| A6.1 | Surface films / slicks | **(a) real — biggest omission** | Monolayer suppresses cm-wave amplitude up to ~×10 (Marangoni damping) |
| A6.2 | Salinity | (a) if brackish, (c) if fresh | dn/dS +1.8×10⁻⁴/g/kg → Δn 0.0064 at S=35; <10⁻⁴ below S≈0.5 |
| A6.3 | CDOM | **(a) real** | a(440) 0.1–20 m⁻¹ vs 0.015 pure water; spectral, not grey turbidity |
| A6.4 | Ice cover | (a) if winter in scope | Phase change replaces the surface entirely |
| B1 | Solar cycle TSI / UV | (b) derived, (c) as a knob | 0.1 % TSI = 1.3 W/m² TOA; varying UV absorbed above stratosphere; orbit gives 6.9 % already |
| B2 | Flares / CMEs | (c) | Enhancement deposited >60 km; GIC needs km conductors (~10⁻² V across scene) |
| B3 | Geomagnetic activity / Kp | (c) | NAO links contested & non-stationary; no local path; energy-capped (B6.1) |
| B4 | GCR → clouds (Svensmark) | (c) | CLOUD: nucleation ×2–10 real, but CCN change over cycle = 0.1 % (Dunne 2016) |
| B5 | Lunar / Earth tide / Coriolis | (c); (b) if tidally connected | Differential lunar tilt ~1 µm/scene; atm. lunar tide ≤0.06 hPa; Ro≈100. Connected tide = boundary input |
| B6.1 | Solar-wind/magnetosphere energy | (c) | ~5×10⁻⁴ W/m² vs 1365 radiative; extreme storms 0.03 W/m² at >90 km |
| B6.2 | Volcanic/stratospheric aerosol | (b) derived | Real (~−3 W/m², Pinatubo) but enters only as irradiance spectrum + direct/diffuse input |
| B6.3 | Geothermal, Schumann/GEC, main B-field, neutrinos, planets | (c) | Each 5–15 orders below the relevant scene energy/force scale |

**Net build list for v5:** T_water (prognostic) → n, ν, σ; RH + T_air → evaporation/
heat-balance (+ steam-fog switch); film coverage F → Marangoni damping; S (if
brackish); CDOM a(440) + slope → spectral absorption. **Net build list for v8:**
nothing — optionally log TSI (LISIRD), Kp/flare class (NOAA SWPC) as scene metadata.

---

## Could not verify

1. **Exact piezo-optic coefficient dn/dP for water at pond conditions.** I derived
   Δn ≈ 8×10⁻⁶ at 5 m from Lorentz–Lorenz + measured β, consistent in order with the
   NIST/Waxler high-pressure tables (which resolve n only to ~10⁻⁵ per ~10 bar in the
   indexed excerpts). The verdict is insensitive to a factor of several here.
2. **TSI excursion during the largest single flare (X28, 4 Nov 2003).** Reported in
   the literature as a small (~0.1 %-class, minutes-long) TSI spike; I did not open a
   primary source stating the exact figure, so the report says "at/below ~0.1 %" rather
   than a precise number. Verdict unaffected.
3. **A single Kp/aa–NAO correlation coefficient.** The literature itself reports a
   non-stationary, sign-changing relationship; no single defensible number exists to
   quote, which is part of the B3 finding. I relied on review-level characterisations
   rather than one primary correlation study.
4. **IPCC AR6 solar forcing +0.05 W/m² and ~0.07 °C cycle imprint** were taken from a
   secondary synthesis explicitly citing AR6 WG1, not from the AR6 chapter PDF itself.
   Both are consistent with the primary TSI numbers verified independently.
5. **Natural film elasticity |E| and the wind speed that destroys pond films** vary
   widely by film chemistry (published |E| spans roughly 5–100 mN/m; breakup winds
   ~3–10 m/s in different studies). The report gives ranges and a tunable parameter
   rather than a single value; the *existence and order* of the damping (×~10 in tank
   experiments) is solidly sourced.
6. **No live interactive verification was performed** for the data feeds named
   (Open-Meteo, NOAA SWPC JSON, LISIRD, NMDB): their existence, fields and access
   pattern are documented in the indexed sources, but I did not fetch a live payload
   in the research browser during this run. Treat endpoint details as `index`-grade.
7. **Steam-fog threshold** (ΔT ≳ 5–10 K) is an empirical rule of thumb from
   meteorological practice, not a sharp physical constant.

*No page encountered instructed this research to do anything; no such instructions
were followed or needed flagging beyond this note.*

---

## Sources (read 2026-11-10, all `index` unless noted)

- Water compressibility / bulk modulus 2.2 GPa; 0.05 % at 100 m — https://www.compadre.org/nexusph/course/Bulk_modulus_--_liquids
- Refractive index vs T (max 1.33434 <0 °C) and piezo-optic behaviour — https://vitroid.github.io/water-science/water/density_anomalies.html ; dn/dT ≈ −0.0001/°C — https://patents.justia.com/patent/8239144
- Viscosity vs T (1.792 / 1.002 mPa·s at 0 / 20 °C, Kestin correlation) — https://statgarden.com/physics/water-viscosity-calculator.html
- Surface tension 72.8 mN/m at 20 °C, IAPWS — https://fluidprops.com/fluids/water/surface-tension/
- Evaporation bulk/Dalton formulae; DYRESM QE = ρa λ CE U (qa−qs), CE≈1.3×10⁻³ — https://storage.googleapis.com/nintione-resources/DKCRC-Working-19-Evaporation.pdf
- Mirages: −3.5 °C/10 cm inferior threshold; superior over cold seas (Stull) — https://geo.libretexts.org/Bookshelves/Meteorology_and_Climate_Science/Practical_Meteorology_(Stull)/22%3A_Atmospheric_Optics/22.6%3A_Mirages
- Thermocline schlieren / "wrinkled glass" from altered refractive index — https://geo.libretexts.org/Courses/Lumen_Learning/Earth_Science_(Lumen)/11%3A_Ocean_Water_and_Life/11.02%3A_Thermocline
- Inverse barometer ~1 cm/hPa — https://www.cambridge.org/core/journals/journal-of-glaciology/article/iceshelf-elevation-changes-due-to-atmospheric-pressure-variations/B582D7579755103D9565CFAAAA7B0DE8
- Miles–Phillips mechanism (Phillips resonance, Miles critical layer) — https://en.wikipedia.org/wiki/Miles-Phillips_mechanism
- Henry's law / dissolved O₂ ∝ atmospheric pressure — https://www.consort.be/Support/Technical_documentation/APP_Dissolved_Oxygen.pdf
- Surface films: cm–dm wave damping, Marangoni resonance — https://www.researchgate.net/publication/262522788_Marine_Surface_Films_in_the_Coastal_Zone_What_can_we_Learn_Using_Remote_Sensing_Data ; fish-oil tank, σ(z) shrinks ~×10, no wind growth from rest — https://os.copernicus.org/preprints/os-2018-111/os-2018-111-manuscript-version6.pdf
- Seawater vs freshwater n (1.33940 vs 1.33300) and density (×1.0264) at S=35 — https://reefkeeping.com/issues/2005-11/rhf/
- CDOM a(440) 0.1–8.5 m⁻¹, up to 20 m⁻¹ inland — https://www.sciencedirect.com/science/article/abs/pii/S0034425713003830
- TSI 1361 W/m², ~0.1 % cycle variation, UV share 6–8 % — https://www.frontiersin.org/articles/10.3389/feart.2023.1164636/full
- Solar-cycle UV variation (few % 200–300 nm … 50–100 % Ly-α) — https://round-lake.dustinice.workers.dev:443/https/www.aanda.org/articles/aa/pdf/2006/23/aa4809-06.pdf
- X-class flare band enhancements (FISM: +680 % 0–14 nm … +35 % 122–190 nm) — https://ia803207.us.archive.org/13/items/NASA_NTRS_Archive_20150007959/NASA_NTRS_Archive_20150007959.pdf
- Magnetospheric energy budget, Oct 2003: 17 TW input, ~30 % Joule heating — https://sci.esa.int/web/cluster/-/40261-rosenqvist-l-et-al-2006 ; solar-wind flux ~5×10⁻⁴ W/m² vs radiative — https://arxiv.org/pdf/0905.1812
- GIC: tens–hundreds of A, Hydro-Québec 1989 collapse — https://en.wikipedia.org/wiki/Geomagnetically_induced_currents
- Geomagnetic activity / stratosphere, NAO-correlation literature — https://www.sciencedirect.com/science/article/abs/pii/S1364682624001159
- CLOUD 2011: ions enhance nucleation up to ×10+ — https://physicstoday.aip.org/news/cern-clears-the-sky-on-cloud-formation ; Dunne et al. 2016 (Science): GCR variations do not appreciably affect climate via nucleation — https://eprints.whiterose.ac.uk/id/eprint/107269/1/global%2520atmospheric%2520particle%2520formation.pdf ; CCN change 0.1 % over cycle (CLOUD team summary) — https://bibbase.org/network/publication/dunne-gordon-krten-almeida-duplissy-williamson-ortega-pringle-etal-globalatmosphericparticleformationfromcerncloudmeasurements-2016
- Lunar atmospheric tide 0.06 mb tropics / 0.02 mb midlatitudes — https://glossary.ametsoc.org/wiki/lunar-atmospheric-tide/
- Solid Earth tide, M2 up to ~55 cm — https://en.wikipedia.org/wiki/Earth_tide
- Coriolis scaling / Rossby number (bath Ro ~3200) — https://www.fluid-flow-aerodynamics.com/essays/the-bath-does-not-know-the-hemisphere/
- Open-Meteo fields (temperature_2m, relative_humidity_2m, surface_pressure, radiation; no key) — https://github.com/lekhanpro/suryagrid/blob/HEAD/docs/DATA_SOURCE_CATALOG.md
- NOAA SWPC Kp / GOES flare products as published feeds — https://voyaggora.com/travel-live/today/space-weather-today/
