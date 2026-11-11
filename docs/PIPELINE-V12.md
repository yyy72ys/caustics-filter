# PIPELINE-V12 — chunk 10: the 3.9 redo (boundary treatment, three-arm measurement)

Written 2026-11-11 after the three-way audit round of 3.1–3.9 (digests rally-v40-{gpt,grok,gemini}.md; round appended to RALLY-V3-GAPS.md). User order 2026-11-11 18:14: audit 3.1–3.9 in the rally, then do 3.9 on the results.

## Audit standings adopted (rally consensus)
- 3.1 STANDS (scope-limited baseline). 3.3 STANDS. 3.4 STANDS as transport machinery (follow-up recorded: enforce the generator/transport separation in code+tests — GPT).
- 3.2: flagged **REJECTED_INVALID_DERIVATION** in bookkeeping (Gemini's label, matching the user's own suspicion); the revisit-before-v4 promise stands. Not deleted — flagged.
- 3.5 WEAKENED: η=1 is a scope decision, not a repair; the diffuse branch stays known-broken and unused.
- 3.6 WEAKENED: merged as best-available, NOT validated (its battery compares the model to itself; the halo was never tested against a measured in-water PSF). Known limit, recorded; no rework in this chunk.
- 3.7: **"validated physical default" status REVOKED (3/3)**. It remains the operational fallback WITH KNOWN DEFECTS DISCLOSED (H_s not from spectrum; ECKV agreement circular, totals 0.81–3.32× true ECKV). Report/board language updated accordingly.
- 3.8 STANDS as challenger. 3.9 WEAKENED: the fit crossing is a fit artefact, not a physical boundary (3/3) — the boundary treatment is this chunk's subject.

## The fork, and how it is settled
- GPT + Grok rank (b) [extend the labelled prior until raw α_m reaches the floor level on its own] first; both reject (d) [redesign the share gate] UNLESS independently justified.
- Gemini ranks (d) first, with a physical argument the other two did not weigh: capillary dissipation ∝ 4νk² is massive at k>100 while gravity dissipation is negligible, so capillary dominance of the 50–200 band at U≤3 may be physically unattainable and the 50% share gate (a high-wind intuition) may itself be wrong at low wind. Its legitimacy test: the ABSOLUTE saturation B(k) must stay inside the equilibrium envelope (B ≤ 0.012 everywhere; B_l ≤ β̂ = 0.008 in the band) and totals inside [Cox–Munk, Hu].
- **Settlement = measure both.** Three arms, same 8-cell matrix, same code path family (run_v14 on the v13 base):
  - **Arm RAW** — v13 as shipped (prior below the crossing only). Control.
  - **Arm B** — prior zone extended to every cell where the raw spectrum delivers σ²_h < 0.003: applied σ²_h = max(raw, 0.003), the α_m solve continuous by construction; the hand-off u* reported per cell; provenance label `low_wind_prior` wherever the top-up acts (including U=3 when it acts).
  - **Arm D** — raw spectrum untouched; at U≤3 the band-share gates are REPORTED but replaced as gates by the absolute criteria (B-envelope + totals envelope + C¹ smoothness). Provenance label `low_wind_absolute_gate`.
- Decision rule (fixed in advance): the arm that satisfies the INDEPENDENT constraints (totals envelope, B-envelope, construction gates, transport sanity) with the fewest invented inputs wins merge-as-challenger. Tie → Arm B (majority rank). Arm D winning requires its B-envelope checks to actually pass at U=3 — assertion is not evidence.

## Gates (synthesis of the three batteries)
Construction (all arms, all cells; any failure = FAIL): A1 positivity; A2 continuity — B(k) continuous and dlnB/dlnk without steps across every join and the prior hand-off (Grok: ≤1% jumps); A3 prior accounting — raw α_m, applied α_m, prior contribution, σ²_h printed separately in every state; A4 Parseval + independent band integrals agree; A5 resolution convergence (grid + k-cutoff).
Physical (three-state, strict meanings — UNRESOLVED is not a soft PASS):
- High-wind cells (U=6, 10): C1 (gravity share ≤50%), C1b (capillary share ≥0.15 AND σ²_h ≥0.002), totals envelope, B-envelope — must PASS.
- Low-wind cells (U=1.5, 3): at best UNRESOLVED (3/3 consensus). Reported: C1, C1b, totals vs [Cox–Munk clean, Hu], B(k) max, B_l band max, σ²_h delivered, prior state. FAIL if: totals outside [0.5×, 2.0×] envelope, any B > 0.012, negative energy, a step in S(k), or |ΔH_s| > 35% (v3.8's K1 line).
- Sensitivity: Arm B floor sweep {0.0015, 0.003, 0.006} (carried from v13).
- Transport: windy cell (must stay in the v3.8/v3.9 observed band: warp 20–45 mm, p99.9 ≥ 3.0, crush ≤5%, blow 0) + one U=3/F=1000 render per arm (RAW vs winner) for the optical consequence of the boundary choice.
- Regression: v13 legacy lane byte-identical; v12 water default path unchanged; frozen composite unchanged.

## Merge rule
Winner arm merges as the default spectrum (elfouhaily-form + boundary treatment, provenance labels on) ONLY IF: construction gates all pass, zero FAIL, high-wind cells PASS, regression green. Low-wind cells stay UNRESOLVED and labelled — no claim of validation at U≤3, ever, from this chunk. 3.7 (jonswap) remains selectable, labelled legacy fallback with defects disclosed. If neither arm clears: no merge; the report says why; default unchanged.

## Explicitly NOT adopted
- Gemini's runtime fallback (flat surface for U≤3): the project's scenes are low-wind ponds; physics-first here means honestly labelled uncertainty in the state, not deleting the surface. Recorded disagreement.
- Grok's merge rule requiring all U≥3 cells PASS: incompatible with the 3/3 consensus that low-wind cannot be PASS without measurement; replaced by the UNRESOLVED-with-labels rule above.
- Gate (d) as a standalone act (changing thresholds because a cell failed): only Arm D's pre-registered absolute criteria are tested, and they must pass on the raw spectrum.

## Outcome addendum 2026-11-11 (chunk 10 measured + B-envelope adjudication)

Chunk 10 (v14) measured all arms: construction gates all PASS; low wind UNRESOLVED in every arm, no FAIL; Arm B buys nothing at U=3 (optics identical to RAW); winner reading = Arm D (fewest invented inputs). MERGE was withheld on ONE clause: the B-envelope cap failed 3 of 4 high-wind cells in ALL arms (0.018324 @ Omega 7.30; 0.025004 @ Omega 10.26; 0.012485 @ Omega 4.87, capillary peak).
Adjudication (literature research + rally round, see RALLY-V3-GAPS.md): the clause is RETIRED as a category error (Phillips beta mean misapplied as a total-spectrum cap; the failing values are Elfouhaily's own published peaks). Replacement = the region-wise two-column gate specified in BRIEF-V15-MERGE.md. Merge approved 3/3 with conditions (labels travel; low-wind UNRESOLVED; two-column reporting; JONSWAP anchor printed; forward 15% H_s tripwire). Chunk 11 (run_v15) performs the merge; on its PASS, SPECTRUM default = elfouhaily-form (Arm-D rules) and v3.9's line closes as v3.10.
