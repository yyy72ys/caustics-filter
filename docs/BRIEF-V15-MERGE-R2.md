# STAGE 0 — LOCATION PREFLIGHT (do this FIRST, before anything else)

Your working directory MUST be D:\yohey2026\dsh\muse\caustics-filter\offline
(the runner log's cwd= line proves it). The base this brief builds on lives
THERE — verified by the dispatcher today: run_v14_spectrum.py and
out\v14\REPORT-V14-GATES.md both exist at that path.

Before stage 1, confirm ALL of these exist (absolute paths):
- D:\yohey2026\dsh\muse\caustics-filter\offline\run_v14_spectrum.py
- D:\yohey2026\dsh\muse\caustics-filter\offline\out\v14\REPORT-V14-GATES.md
- D:\yohey2026\dsh\muse\caustics-filter\offline\out\v14\numbers-v14-spectrum.json
- D:\yohey2026\dsh\muse\caustics-filter\offline\out\v10\ (jonswap reference outputs)
- D:\yohey2026\dsh\muse\caustics-filter\offline\run_v12.py and v12_water.py

If any is missing, STOP and report exactly which — do NOT search other
trees and do NOT substitute another directory. In particular,
D:\mutaro-work\caustics-pub is the PUBLISH copy (v5 generation only); it is
NOT the pipeline and must not be touched.

---

# BRIEF-V15-MERGE — chunk 11: the merge (v3.10)

## What this is

The B-envelope clause (max B(k) <= 0.012 over [10,1500]) has been RETIRED by
literature + a three-way rally (2026-11-11; research report:
research_notes/benvelope-cap-literature-20261111-1014; digests rally-v41-*).
The clause was a category error: 0.012 is Phillips' equilibrium-range
coefficient beta (a MEAN, measured range 0.006-0.024), never a cap on the
total saturation spectrum. Elfouhaily's own B_l(k_p) = 0.1433*alpha_p*gamma
reproduces the three "failing" v14 values exactly (0.0183 @ Omega 7.3 flank,
0.0250 @ Omega 10.3 flank, 0.01249 = (1/2)*alpha_m at U=10).

All three reviewers APPROVED merging the elfouhaily-form family as the
default spectrum, under conditions. This chunk performs the merge. The
winner reading from chunk 10 stands: Arm D (v13/v14 RAW numbers; low wind
judged by absolute criteria; low-wind cells UNRESOLVED and labelled).

Code family: NEW files run_v15* in the same offline directory, built on
run_v14 (which carries v13's spectrum). DO NOT modify any v11/v12/v13/v14
file. Jobkit checkpoint per stage (out/v15/checkpoint.json); run in the
FOREGROUND until every deliverable below is written. If any gate FAILs,
write the report with the failure and STOP — no tuning, no threshold edits.

## Stage 1 — default switch + region-wise gate (run_v15_spectrum.py)

1. SPECTRUM DEFAULT: expose SPECTRUM_SOURCE_DEFAULT = 'elfouhaily-form' in
   the v15 module, with the Arm-D rule set: alpha_m = v13 as shipped (prior
   below the ECKV fit crossing solving sigma2_h,total = 0.003, raw inside
   the domain; provenance low_wind_prior / eckv_fit_domain_raw as in v14).
   The 'jonswap' source remains selectable (rollback).
2. RETIRE the universal B <= 0.012 clause everywhere in the v15 gate code.
3. REGION-WISE B GATE, reported in TWO columns per cell (never collapsed):
   CONFORMANCE (model vs its own published form):
   (a) gravity: B_l flank max over [k_p/2, 2k_p] vs
       0.1433 * alpha_p * gamma(Omega_c) * 1.04; PASS if ratio <= 1.5.
   (c) capillary: B_h max in [k_m/2, 2k_m] vs (1/2)*alpha_m(u*);
       PASS if the ratio is within [0.8, 1.25].
   OBSERVATION (independent anchors):
   (b) equilibrium beta: reconstruct E(f) from the assembled omnidirectional
       spectrum via deep-water dispersion (omega^2 = g k, f = omega/2pi);
       least-squares fit in log space of E(f) = beta * I(p) * g * u* * f^-4
       over f in [1.3 f_p, 3 f_p] with I(p) = 2.5 (state the convention);
       PASS if fitted beta in [0.006, 0.024]. If inside that band B_h exceeds
       B_l anywhere (capillary contamination), mark this sub-check
       UNRESOLVED (band-contaminated), never a silent PASS/FAIL.
   (d) JONSWAP anchor (context, not a gate): print
       B_J(k_p) = 0.5 * alpha_J * exp(-1.25) * 3.3,
       alpha_J = 0.076 * (g F / U^2)^-0.22, beside the model's B_l(k_p) for
       every cell; print the sub-km extrapolation caveat for F=100 cells.
4. C1/C1b band shares: computed and reported for all cells; GATING only for
   U >= 6 cells (as in v13/v14); at U <= 3 reported-not-gated.
5. Totals envelope [Cox-Munk clean, Hu] with the 0.5x-2.0x clause: unchanged,
   all cells, gating.
6. K1 |dHs| (spectrum-integrated H_s vs fetch law): FAIL if > 35%
   (pre-registered threshold, unchanged). Additionally set hs_flag = true
   when |dHs| > 15% (a forward tripwire adopted from the rally; on THIS
   matrix it annotates, it does not re-gate — say so in the report).
7. PROVENANCE in every per-cell record AND in the numbers JSON: omega_c,
   warn_extrapolated_wave_age = (Omega_c > 5), provenance label, low-wind
   status. Low-wind cells (U <= 3): verdict UNRESOLVED unless a hard clause
   fires (totals outside 0.5x-2.0x, negative energy, a step in S(k),
   |dHs| > 35%, conformance ratio (a) > 1.5) — then FAIL.
   High-wind cells: PASS requires construction gates + totals + C1/C1b +
   K1 + conformance (a),(c) PASS + observation (b) PASS or flagged
   UNRESOLVED-with-reason. Zero FAILs required for the merge to stand.
8. Print the full 8-cell matrix (U in {1.5,3,6,10} x F in {100,1000}) with
   both columns, all flags, and verdicts into numbers-v15-merge.json.

## Stage 2 — regression (all must hold; report each)

- D_STRICT branch == v11.assemble, byte-identical spectrum (as in v13/v14).
- v15 RAW arm numbers == v14 RAW numbers per cell (alpha_m, sigma2_g,
  sigma2_h, B max, H_s pred): identical code path, expect exact equality;
  report any delta.
- v12 water model untouched: WATER_MODEL_DEFAULT still 'ab'; its G1 legacy
  lane still byte-identical to delivery still-A-v37 (sha cfb086b22a0e...).
- Frozen composite modules (the v12 hash list): unchanged.
- 'jonswap' still selectable in v15 and reproduces the v10 scene-A spectrum
  statistics (compare realised sigma2 and H_s to numbers-v10 within 1%).

## Stage 3 — transport through the MERGED DEFAULT stack

Stack: v15 spectrum (default elfouhaily-form) + v12 water (ab) + frozen
composite (A2.6 + Blow + E2.2), eta = 1.0, retina OFF, rain 0. Three stills:
- T1: U=6, F=1000, depth 1.0 (the standing windy comparison scene)
- T2: U=3, F=1000, depth 0.6 (low-wind; UNRESOLVED label must appear in its
  numbers JSON)
- T3: U=6, F=100, depth 1.0 (Omega_c = 7.30; warn_extrapolated_wave_age must
  appear as true in its numbers JSON — the warning travelling with an
  actual output is a merge condition)
Each numbers JSON carries: spectrum_source, omega_c,
warn_extrapolated_wave_age, cell status, C mean / p99 / p99.9 / max on the
unattenuated field, crush %, blow %, own geometric shift (mm).
Gates: C mean 1 +/- 0.02; crush < 5%; blow == 0. Any FAIL: report and STOP.

## Stage 4 — report (out/v15/REPORT-V15-MERGE.md)

- The two-column gate table for all 8 cells + verdicts.
- The merge statement: default spectrum is now elfouhaily-form (Arm-D rules);
  jonswap retained as rollback; the 0.012 clause retired with its reason.
- Flags inventory: cells carrying warn_extrapolated_wave_age; UNRESOLVED
  cells; hs_flag cells (with their dHs values).
- Regression results (stage 2, item by item).
- Transport results (stage 3) incl. proof the labels travelled into T2/T3
  numbers JSONs.
- Hash list of new v15 files; explicit statement that no v11-v14 file was
  modified (hash-compare the v14 report's recorded list where available).
- Label honesty: call the new battery a "region-wise conformance +
  observation battery", never "physical truth" (rally condition).

## Deliverables (VM copies will be pulled from out/v15)

out/v15/REPORT-V15-MERGE.md, out/v15/numbers-v15-merge.json,
out/v15/numbers-v15-transport.json, out/v15/still-T1-windy.png,
out/v15/still-T2-lowwind.png, out/v15/still-T3-warned.png,
out/v15/checkpoint.json (jobkit, one unit per stage).
