# PIPELINE-V11 — chunk 9 = version 3.9: the low-wind capillary floor, honestly labelled

Origin: v3.8's merge was withheld on C1/C2. Research (research_notes/
low-wind-capillary-floor-20261111) reframed C1's failure as TWO defects: an
implementation-invented hard zero in ECKV alpha_m, and a gate that is
unsatisfiable at U = 1.5. The three-way rally (RALLY-V3-GAPS.md, v3.9 round;
digests rally-v39-{gpt,grok,gemini}.md) then REJECTED the research's proposed
fix (a)+(d) at all three AIs. This design implements what the rally converged
on, not what the research proposed.

## Rally consensus (the constraints on this design)
- R1 SURVIVES (3/3): the alpha_m zero is a fit artefact; the hard zero must go
  as a physical claim.
- R3 WEAKENED (3/3): a floor at sigma2_h = 0.003 is NOT a measured pond
  capillary level. Admissible only as a labelled provisional prior, with its
  uncertainty swept and its provenance printed wherever it acts.
- Decision (a)+(d) REJECTED (3/3): C1 is not re-scoped into green. Gates become
  THREE-STATE (PASS / FAIL / UNRESOLVED); low-wind C1 = UNRESOLVED.
- Grok's condition: one wind-independent minimum-capillary test for ALL cells,
  so the floor must prove it supplies short-wave variance.
- Gemini's condition: report the absolute saturation B in 50–200 against the
  ECKV/Phillips envelope (informational), and note the viscous-transition
  physics (a real capillary collapse exists near u* ~ 0.04–0.08, far below the
  fit crossing at 0.085 — the floor prior must not be mistaken for it).
- GPT's condition: matched-definition care for R2 (recorded reservation: C1's
  unsatisfiability is stated under the pipeline's own MSS conventions, which
  are the conventions the gate itself uses).
- All three name the same decisive measurement (lab/field slope spectrum at
  U = 1.5–3 resolving k > 200 rad/m). We cannot make it; the design lives with
  UNRESOLVED honestly instead of papering over it.

## The change (spectrum level only; transport and water untouched)
Base: v11 code family (v11_spectrum.py / v11_state.py), imported, never edited.
New family: v13_spectrum.py / v13_state.py, out/v13/.

1. **The hard zero is removed.** alpha_m_raw is used as ECKV gives it inside
   the fit domain (u* >= c_m/e).
2. **Below the domain, a labelled prior.** alpha_m_used = max(alpha_m_raw,
   alpha_m_floor(cell)), where alpha_m_floor(cell) solves sigma2_h,total =
   0.003 for that cell (research: alpha_m ≈ 0.0022–0.0027). Every state built
   with the prior active carries `low_wind_prior: {sigma2_h_target: 0.003,
   source: "Cox-Munk crosswind intercept — ocean residual, mechanism NOT
   claimed", rally: "v3.9 R3 weakened at 3/3"}`. The old label
   `eckv_domain: low_wind_floor` is retired.
3. **Sensitivity sweep (mandatory deliverable).** sigma2_h target swept over
   {0, 0.0015, 0.003, 0.006} on the preset-1 scene and the windy scene: C
   percentiles, own-field warp, fine texture |dL/dy|, total sigma2 — so the
   prior's influence on the picture is measured, not asserted.

## Gates (three-state; this chunk's battery)
- A1/A2 (algebra, Parseval): as v11, re-run on the new source. PASS/FAIL.
- C1 (50–200 gravity share ≤ 50%): three-state. u* >= c_m/e cells: PASS/FAIL as
  before. Below: UNRESOLVED — reported with the numbers, never green, never
  waived.
- C1b (Grok's floor-proof, NEW, all cells, PASS/FAIL): capillary share of the
  50–200 band ≥ 15% OR sigma2_h,total ≥ 0.002. The prior must do measurable
  work everywhere, or the chunk fails.
- C2: RETIRED with its reason recorded (premise false, v11 §2.2). Replaced by
  reporting: total sigma2 beside Cox–Munk clean and the Hu law per cell, and
  the full band partition table. Reporting is not a gate.
- V1–V3, K1, K5, K6: as v11, re-run. K1 flag rule unchanged (>35% blocks).
- K3 (windy compatibility): re-rendered with the floored source vs v11 legacy;
  attribution table as in v11 §3.2.
- K4: v11 regression battery verbatim (bit-identity on identical buffers holds
  for all cells where the prior is inactive — i.e. u* >= c_m/e cells and the
  whole legacy path).

## Merge rule (re-run of the v3.8 evaluation)
Default spectrum_source flips to `elfouhaily-form` iff: A1/A2 PASS; C1 PASS on
every in-domain cell; C1b PASS on all cells; V reported; K1 no flag; K4–K6
PASS; the sensitivity sweep is delivered; and the low-wind cells are listed in
the report as UNRESOLVED (rally consensus) with the prior labelled. UNRESOLVED
is disclosed, not counted as PASS. If any FAIL: default stays `jonswap`, the
report says which and why. No blends, no rescaling, no third option.
