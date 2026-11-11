# PIPELINE-V9 — chunk 7 = version 3.7: two-regime spectrum determination

Status: DESIGN, synthesised 2026-11-11 from (1) the measured diagnosis
(delivery/spectrum-diag/), (2) deep research
(~/workspace/research_notes/unified-spectrum-pond-slope-height-20261110-2200/),
(3) the three-way rally in the v3-gaps threads (digests
rally-v37-{gpt,grok,gemini}.md). USER PRE-APPROVED 2026-11-11 07:05:
implement immediately after design, report results.
Code family: run_v10* (sequence v5=3.1 … v9=3.5, v10=3.7).
Scope: the ENVIRONMENTAL STATE GENERATOR only (spectrum → surface
fields). Transport, frozen composite, sky, retina: untouched, guarded
by regression gates.

## 1. What was wrong (measured, 2026-11-11 diagnosis)

Two INDEPENDENT defects in the v3.4 state generator:
- D-A: realised H_s = 4 × waveHeightM(env), a wind-only amplitude
  curve (no fetch, no depth). The spectrum never sets the height.
  Ratio vs the fetch law 2.9–14.9× in all 8 sweep cells; ×10 fetch
  moves realised H_s ≤ 1.2%.
- D-B: no capillary tail (<1 cm band ≈ 1e-8 of σ²; missing capillary
  σ² 0.011–0.034 rad²/cell), no resolved gravity peak (Phillips peak
  at 2.22× patch; λ_p = patch size in every cell), and the σ² anchor
  is a taper that reaches its target by DELETING k > k_c. At U=10 the
  anchor silently no-ops (σ² misses by 59–67%).
GPT's framing is adopted: v3.4 is rejected as a physically valid
environmental wave generator; its transport machinery stands.

## 2. Conventions (pinned — GPT's trap, closed by gate G0)

One convention in the new module: omnidirectional wavenumber spectrum
S(k) [m² per rad/m] with m0 = ∫₀^∞ S(k) dk and σ² = ∫₀^∞ k² S(k) dk.
Directional form Ψ(k,θ) = S(k)·Φ(k,θ)/(2πk) on the (kx,ky) grid with
∫₀^{2π} Φ dθ = 1. Wind convention frozen: U10 everywhere (Open-Meteo
real mode already supplies wind_speed_10m); u* derived from U10 only.
Gate G0 (Parseval): the built fields' numeric m0 and σ² must match
the analytic integrals of the assembled S(k), per band, within the
model's own realisation-scatter band. Convention errors cannot hide.

## 3. Spectrum construction — additive, NO splice

S(k) = S_g(k) + S_h(k). No blend weights, no splice wavenumber, no
post-hoc renormalisation of either part (GPT's forbidden list).

- S_g (gravity): JONSWAP in frequency form S_J(ω; H_s, f_p) from the
  fetch laws, converted by the exact change of variables
  S_g(k) = S_J(ω(k)) · dω/dk with ω(k) the existing finite-depth
  gravity–capillary dispersion. This plays Elfouhaily B_l's role
  without extrapolating B_l's wave-age fit beyond its range
  (Ω_c ≈ 6.3 > 5 at F=100).
- S_h (capillary): Elfouhaily B_h(k; u*) converted into the same
  S(k) convention and ADDED. Its level is set by u* alone; it is
  never scaled to hit any σ² target.
- Overlap honesty: the 50–200 rad/m band contributions of S_g and
  S_h are reported separately in every gate run (double counting is
  bounded by JONSWAP's steep tail; it is measured, not assumed away).
- Φ(k,θ): ECKV's unified azimuthal function applied to the total
  (Δ(k) interpolating gravity ↔ capillary regimes). Validation:
  σ²_along/σ²_cross vs the Cox–Munk anisotropy table per cell (G5).
- Low-wind domain (GPT N2): α_m = 0.01·[1 + ln(u*/c_m)] < 0 for
  u* ≲ 0.085 m/s. Treatment is EXPLICIT: α_m floored at 0.0, the
  numbers file carries `eckv_domain: low_wind_floor` for such cells,
  and gate G4 checks B(k) ≥ 0 for all k in all cells. No silent
  clamping anywhere: unattainable conditions RAISE.
- High-k termination: B_h's own roll-off plus the established
  viscous envelope; the integration cutoff is decided by the G6
  sensitivity sweep (k_cut ∈ {800, 1000, 1200, 1500}), not fixed
  by fiat.

## 4. Amplitude (D1 revised)

- H_s comes from the fetch laws (the pipeline's existing law family)
  — waveHeightM is RETIRED as a source.
- Depth: a Young–Verhagen-style depth-limited GROWTH correction is
  applied only when the regime binds (k_p·h < 1.0); its magnitude is
  logged per cell either way. Finite-depth dispersion was already
  correct and stays.
- Every state is labelled fetch-equilibrium; the minimum-duration
  estimate t_min(F, U) is reported in the numbers (GPT). Duration-
  limited growth is recorded as future work, not silently assumed.
- Recipe explicit H_s override: honoured (it is an environment
  statement, like depth), and logged as `hs_source: law|override`
  in every numbers file — never silent (Gemini's escape-hatch point).
- The breaking clamp min(0.142λ, 0.78h) stays as a post-check.
- OPEN REGRESSION, stated plainly: preset realised values WILL
  change (preset-1 σ² target 0.0210 retired with the taper; H_s
  drops to law values). The physics-strict baseline is re-anchored
  to the laws by this chunk; the comparison stills are deliverables.

## 5. Patch and cascades (D4 revised)

- Patch is a SCENE property (recipe/scene scale), wind-independent.
  The wind-coupled patch growth of v3.4 is retired.
- The ONE spectrum of §3 is realised on two cascades when one grid
  cannot serve both ends: gravity cascade (tile S_lo ≥ 4λ_p,
  convergence-checked by tile doubling) + capillary cascade (tile
  S_hi with k_Nyq,hi ≥ 2·k_cut). The split is a partition of unity
  in ln k (W_lo + W_hi = 1) applied to the single spectrum — no
  double counting by construction; both cascades quantise ω to the
  SAME 2π/8 s grid and share the seed family (Gemini's beat trap).
- Single cascade only when one grid provably satisfies both ends;
  which mode ran is recorded in the numbers.
- Gates are CONVERGENCE tests, not the bare λ_p ≤ S/4 rule: doubling
  S_lo and N must move peak location, MSS, slope covariance and
  C-statistics within stated tolerances.

## 6. Time (D5 revised)

Dispersion unchanged. The 8 s loop quantisation is a frozen product
interface and stays, but its error is MEASURED per band (max
|Δω|/ω, worst at the longest modes; v3.4 reference: shallow phase
speed −4.16%) and printed in every gate report (G7).

## 7. Acceptance battery (replaces the retired scalar anchor)

- G0 Parseval per §2.
- G1 H_s: realised (spectral AND spatial readouts) vs law ≤ ±15%
  across the 8-cell sweep (U ∈ {1.5,3,6,10} × F ∈ {100,1000}, z=1.0).
- G2 σ² honesty: realised σ² vs the assembled spectrum's own
  integral; band shares (>1 m / 0.1–1 m / 1–10 cm / <1 cm) reported;
  fractions above k=100 and k=370 compared with the research
  envelopes (≈50–70% / ≈25–40%) as VALIDATION, never forced.
- G3 MSS references: realised σ² printed beside Cox–Munk clean,
  the retired 0.0548 target, and the full-ECKV integration. None is
  a target (GPT/Gemini consensus: MSS is an output).
- G4 positivity + domain log; zero silent no-ops.
- G5 anisotropy per §3.
- G6 cutoff sensitivity convergence per §3.
- G7 quantisation error report per §6.
- G8 REGRESSION (frozen paths): (a) frozen composite on identical
  buffers bit-identical; (b) Beer–Lambert outputs for identical
  optical inputs unchanged (wave depth and optical path stay
  separate dependencies — no double counting); (c) η=1 direct route
  golden test on the reference scene; (d) rain R=0 equals the new
  base state exactly (3.2 gate form adapted to the new generator).
- G9 transport sanity on the new states: photon cross-check +
  C percentiles on the reference scenes (C mean = 1.0 invariant;
  p99.9/mean vs the 6–12 band; deliverable-still MS-SSIM). This is
  where Gemini's Hessian/quadrature concern gets its verdict. If G9
  fails by shredding, STOP and report — the sub-footprint
  micro-normal route (research §, Bruneton-style) is a separate
  design step, never improvised inside this chunk.
- G10 windy-pond re-render: the 2026-11-11 windy recipe (wind 6,
  fetch 100, z=1.0, moth photo) under v3.7; warp yardsticks
  (block-matching + the pipeline's own shift field, in px and mm)
  compared against the v3.4 still and yesterday's production frame;
  picture delivered. Expectation stated in advance: with the true
  H_s the broad warp weakens and the fine web strengthens — that is
  the physics, and the numbers must show which band moved.

## 8. Recorded, out of scope

- Directional fetch F(θ_wind, shoreline) (GPT N1): fetch stays a
  recipe input in v3.7; directional effective fetch is a noted
  recipe extension for later.
- Surface-film damping of B_h: this is where v5's film factor F
  will couple in; recorded at the seam, not built here.
- Rain re-derivation (user's standing rain revisit): the rain bump
  at k≈420 rad/m now sits inside a real capillary band; the revisit
  will be cleaner on top of v3.7. Separate decision, unchanged.
- Duration-limited growth: labelled only (§4).

## 9. Implementation staging (jobkit, foreground, deepseek-v4.1-flash)

1. v10_spectrum.py: assembly + analytic integrals + G0/G4 checks.
   No state or render changes.
2. State builder v10 (single + two-cascade) + 8-cell sweep runner:
   G1/G2/G3/G5/G6/G7 tables.
3. Transport hookup on the new states + G8 regression battery +
   G9 photon/C gates on the reference scenes.
4. G10 windy re-render + comparison + REPORT-V10-GATES.md.
Deliverables: PC offline/out/v10/, VM delivery/offline-v10/.
