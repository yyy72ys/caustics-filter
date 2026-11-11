# PIPELINE-V10 — chunk 8 = version 3.8: gravity-tail termination

Design basis: research `research_notes/gravity-tail-termination-spectrum-20261111-0528/report.md`
(Rank 1) + three-way rally D1–D5 (digests `rally-v38-{gpt,grok,gemini}.md`, round in
`RALLY-V3-GAPS.md`). User decision 2026-11-11: proceed with the fix — "correcting the
counting, not dressing the look".

## 1. What changes, and what does not

v3.7's measured defect, restated precisely: the total slope variance is right (within
1–2% of the full-ECKV total) but its band distribution is wrong. In 50–200 rad/m the
converted JONSWAP tail sits at the plateau B_g ≈ α_J/2 ≈ 0.006–0.023 and carries 87%
of the band's slope variance; the sub-centimetre band is correspondingly starved.
The dispersion conversion already kills the tail beyond ~370 rad/m (factor 0.25 at
370, 0.02 at 1000) — so this is a LEVEL defect in 50–200, not a cutoff defect, and no
rule acting at k_t can fix it (rally D1: survives at Gemini/Grok; GPT's convention
check is adopted below as gate A1).

**Change**: the gravity part of the omnidirectional spectrum is replaced wholesale:

```
S(k) = k^-3 * ( B_l(k) + B_h(k) )
```

B_h is v3.7's Elfouhaily high-frequency part, unchanged (imported from the v10
modules). B_l is Elfouhaily's low-frequency part (§3). There is no join, no blend,
no splice, and no renormalisation of either part — the disqualifications recorded
in the research (§8) and the v3.7 rally (D2) stand.

**Honest label** (rally D3, all three AIs): the result is NOT "the Elfouhaily
model". It is an **Elfouhaily-form spectrum with fetch-derived sea-state
parameters** — ECKV's shape functions driven by this pipeline's pond fetch laws.
Code, numbers files and reports use that label. Claiming more would be claiming a
validation that does not exist.

**Unchanged**: transport, frozen composite (A2.6 + Blow + E2.2), sky (chunk 3.5),
retina, directional distribution Φ (ECKV unified function on the total spectrum),
patch/cascade realisation, 8 s quantisation interface, recipe schema. The v3.7
state generator remains in code as the legacy source (see §5, challenger rule).

## 2. Conventions

Inherited verbatim from PIPELINE-V9 §2: omnidirectional S(k) [m² per rad/m],
m0 = ∫S dk, σ² = ∫k²S dk, Ψ = S·Φ/k with ∫Φ dθ = 1, wind = U10 everywhere,
u* from U10 only. One addition: the curvature spectrum B(k) := k³·S(k) is a
first-class reported quantity in this chunk (the defect and its fix both live
in B, not in S).

## 3. The gravity part B_l (Elfouhaily et al. 1997, verified formulas)

```
B_l(k) = (1/2) * α_p * (c_p / c(k)) * F_p(k)
F_p(k) = L_PM * J_p * exp( -(Ω_c/√10) * ( √(k/k_p) - 1 ) )
L_PM   = exp( -1.25 * (k_p/k)^2 )
J_p    = γ^Γ,   Γ = exp( -( √(k/k_p) - 1 )^2 / (2σ^2) ),   σ = 0.08 * (1 + 4*Ω_c^-3)
α_p    = 6e-3 * Ω_c^0.55
γ      = 1.7                                  (Ω_c <= 1)
       = 1.7 + 6 * log10(Ω_c)                 (Ω_c > 1)
```

Parameter sources — every input is already in the pipeline, nothing new is fitted:

- k_p: from the existing fetch laws exactly as v3.7 derived the JONSWAP peak
  (f_p from the fetch law, k_p through the existing finite-depth gravity–capillary
  dispersion). c_p = ω(k_p)/k_p from that same dispersion. Ω_c = U10 / c_p.
- c(k): phase speed from the same dispersion. Recorded deviation: ECKV-as-published
  uses deep-water c(k); at pond depths kh ≫ 1 throughout the tail (kh ≥ 15 in the
  defect band), and the peak region carries the difference — this is a stated
  deviation, not a free parameter (research §1).
- u* does not enter B_l (it sets B_h only, as in v3.7).
- ECKV's own fetch→Ω_c relation is NOT used (ocean-fetch fit; research §1).

Domain labels, per cell, in every numbers file:

- `eckv_domain: published` for Ω_c ≤ 5.
- `eckv_domain: wave_age_extrapolation` for Ω_c > 5 — the γ log-formula is continued
  (γ = 6.49 at Ω_c = 6.3), α_p likewise; nothing is clamped, nothing diverges, and
  the label says plainly that these cells are uncalibrated extrapolation (rally D3).
- The v3.7 low-wind B_h floor (α_m at U = 1.5) is inherited unchanged and keeps its
  `low_wind_floor` label.

## 4. H_s semantics under the new source

B_l is α_p-normalised, not H_s-normalised. Under the new source:

- H_s is a PREDICTION: H_s_pred = 4√m0 of the assembled spectrum, computed and
  reported per cell beside the fetch-law value (hs_law) — never rescaled, never
  tuned toward (research §1: a per-cell rescaling factor is a free parameter and
  would rescale the tail back up by the same factor).
- A recipe H_s override is NOT applied to the spectrum under the new source. It is
  carried into the numbers as `hs_override` (a statement the scene made) beside
  `hs_predicted` (what the spectrum says). Under the legacy source the v3.7
  override semantics are unchanged.
- Depth-limited growth (v3.7's regime-gated correction) does not transfer: it was
  a property of the fetch-law H_s route. Under the new source, depth enters via
  the dispersion (c(k), c_p) and is reported via kh_p per cell. The breaking clamp
  min(0.142λ, 0.78h) stays as a post-check on built fields, as in v3.7.

## 5. Challenger rule and merge criteria

The new spectrum is implemented as a selectable state-generator source
(`spectrum_source: jonswap | elfouhaily-form`) alongside v3.7's, sharing everything
downstream. The default flips to `elfouhaily-form` only when the merge criteria
below all hold in the implementation's own gate report; otherwise the default
stays `jonswap`, the challenger remains available behind the switch, and the
report states why. No third option may be improvised by the implementer — in
particular no blend of the two sources and no rescaling (rally D5 / GPT's
compatibility-gate disposition).

### A — algebra gates
- **A1 (D1 verification)**: compute B_g(k) numerically from the CURRENT v3.7 code
  path (implemented S_ω → exact transform → B = k³S) under the pipeline's pinned
  convention, and compare pointwise with the research's closed form
  B_g ≈ (α_J/2)·(1+3κ²)/(1+κ²)³·exp[−1.25(ω_p/ω)⁴], κ = k/k_m: deep-water ratio
  must be 1 to numerical precision across 50–200 rad/m; finite-depth (h = 0.3,
  0.5, 1.0 m) ratios reported. If the finite-depth plateau moves ≥ 30% from the
  deep-water value, the D1 diagnosis is weakened — STOP and report before
  proceeding (Grok's failure condition).
- **A2**: Parseval on a built field from the new spectrum (v3.7 G0 analogue):
  numeric m0/σ² vs mode sums vs continuum integrals, |z| ≤ 3 convention test.

### C — conformance (internal consistency; NOT physical validation)
- **C1**: per-cell band table [k_p,10] / [10,50] / [50,200] / [200,370] / [370,1000]
  / [1000,∞), split by part. Gate: gravity-part share of the 50–200 slope variance
  ≤ 50% in all 8 cells (research expectation ≲ 20%); max B_gravity in 50–200
  ≤ 0.010.
- **C2**: total σ² per cell within ±10% of the v3.7 totals (the totals already sit
  within 1–2% of full-ECKV; this anchors continuity of the total while the
  distribution changes). Band ratios vs the full-ECKV reference at matched inputs
  reported per band WITH absolute contributions (a ratio alone can hide
  compensating errors — GPT).
- **C3**: fraction of σ² above k = 370 vs ECKV's own fraction at matched
  (k_p, Ω_c, u*), target ±0.03. Labelled conformance in the report, in those words.

### V — validation (nature-facing, uncertainty-bearing; no pass/fail theatre)
- **V1**: total σ² vs clean Cox–Munk (0.003 + 0.00512·U) per cell, factor reported
  (v3.7 observed 1.14–2.34×).
- **V2**: band shares beside the Laxague et al. partition (gravity–capillary bulk;
  pure capillary k>1173 ≈ 5%) as a comparison with stated uncertainty, including
  the caveat that Laxague's seas were Ω_c ≤ 2.
- **V3**: B(k) across 50–200 beside the measured saturation level B̂ ≈ 0.007–0.008:
  the gravity part must no longer sit above it across the band (report the max
  ratio per cell).
- **Recorded known gap** (Gemini, D2): ECKV has no parasitic-capillary source;
  at Ω_c ≈ 4–6.3 a real pond's steep crests may carry MORE sub-cm slope than this
  spectrum. This chunk does not close that gap; it is the stated reason the
  Kudryavtsev challenger remains on the roadmap (D4: deferred at all three AIs).

### K — compatibility (the renderer must not be silently re-calibrated)
- **K1**: ΔH_s = H_s_pred − H_s_law per cell, reported. |ΔH_s| > 35% of H_s_law
  FLAGS the cell: a flag blocks the merge pending review (it does not declare
  the spectrum physically wrong — the fetch law is not ground truth in young,
  duration-limited cells; GPT + Gemini agree the v3.8 draft's 15% rejection was
  a poison pill because the two H_s scalings have incompatible exponents).
- **K2**: height quantiles of built fields (η p50/p90/p99), old vs new source,
  same seed family, per reference scene.
- **K3**: the v3.7 G10 windy scene (wind 6, fetch 100, z = 1.0, moth photo,
  η = 1.0, sun 2026-06-10 09:30, frozen composite) rendered with BOTH sources
  through the unchanged transport + composite: C percentiles, warp yardsticks
  (own geometric shift field in mm + block matching), fine-texture statistic
  (mean |dL/dy|), crush/blow. Expectation stated in advance: fine web holds or
  strengthens, broad warp follows the new H_s; any change must be attributable
  to the band redistribution or to ΔH_s, and the report must say which.
- **K4**: the v3.7 G8 regression battery verbatim — frozen composite bit-identical
  ON IDENTICAL INPUT BUFFERS (this is a code-path regression test; renders under
  the new state are expected to differ and are judged by K3, never by bit
  identity); Beer–Lambert unchanged for identical optical inputs; rain R = 0
  equals the new base state exactly.
- **K5**: 8 s quantisation per-band report for the new spectrum (v3.7 G7 analogue).
- **K6**: continuity audit — B(k) and d ln B / d ln k sampled on a log grid across
  1–3000 rad/m: no jumps; B ≥ 0 everywhere; domain labels present in every cell.

### Merge rule
Default flips to `elfouhaily-form` iff: A1, A2, C1, C2 pass; C3, V1–V3 reported;
K4, K5, K6 pass; K2, K3 reported; and no K1 flag is open. Any failure or flag:
default stays, report says why, and the fix — if one exists within this design —
is a separate decision, never an in-session improvisation.

## 6. Out of scope (recorded, not silently dropped)
- Kudryavtsev B_p + B_eq — future challenger (needs its α(k/k_g), n(k/k_g),
  φ_filt forms verified from print sources first).
- Directional fetch F(θ); surface-film damping of B_h (the v5 seam); rain
  re-derivation (standing user decision, before v4); duration-limited growth
  (label only, as in v3.7).

## 7. Code family
run_v11* (sequence: v5=3.1 … v10=3.7, v11=3.8). New: v11_spectrum.py (B_l +
assembled spectrum + algebra-gate helpers), v11_state.py (source switch over
the v10 state builder). v5–v10 modules are imported, never edited.
