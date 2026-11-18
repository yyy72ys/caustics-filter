# BRIEF-V14-39REDO — chunk 10: the 3.9 redo, three-arm boundary measurement

Read FIRST (in D:\mutaro-work\rally-in\): PIPELINE-V12.md — the design; it governs everything below. Also on the PC already: run_v13.py + v13_spectrum.py/v13_state.py (the base), run_v11.py/run_v12.py (references).

## Job
New code family run_v14* writing ONLY to offline/out/v14 (your own folder; other sessions run concurrently — never write elsewhere, never kill any python.exe/blender.exe process). Foreground execution; jobkit checkpoint after each stage; do not end the session until REPORT-V14-GATES.md is written.

Base: the v13 tree (Elfouhaily-form spectrum + labelled provisional prior below the fit crossing, spectrum_source switch, frozen composite, WATER_MODEL_DEFAULT=ab untouched).

Implement THREE arms over the same 8-cell matrix (U ∈ {1.5, 3, 6, 10} × F ∈ {100, 1000}, depth from the v13 cells):
- Arm RAW: exactly v13 behaviour (prior only where u* < c_m/e).
- Arm B: prior zone extended — wherever the raw spectrum delivers σ²_h,total < 0.003, top up α_m by the same solve v13 uses so σ²_h,total = 0.003; applied level = max(raw, floor) so the hand-off is continuous by construction. Report per cell: raw α_m, applied α_m, prior contribution, whether the top-up acted. Provenance string `low_wind_prior` wherever it acts.
- Arm D: raw spectrum untouched; it changes NO numbers — it is an evaluation arm. At U≤3 cells, evaluate Gemini's absolute criteria instead of the band-share gates: max B(k) over k∈[10,1500] ≤ 0.012; B_l in 50–200 ≤ 0.008; totals inside [Cox–Munk clean, Hu] (matched definitions, as in v13's reporting); C¹ smoothness of dlnB/dlnk (no steps). Provenance string `low_wind_absolute_gate`.

## Gates (from PIPELINE-V12.md — apply exactly)
Construction, all arms/cells, any failure = FAIL: positivity; continuity (B continuous, no step in dlnB/dlnk across joins and the prior hand-off; jumps ≤1%); prior accounting printed per cell (arm B); Parseval + independent band integrals agree; resolution convergence (one cell re-run at finer grid + higher k-cutoff, totals stable).
Physical, three-state (PASS/FAIL/UNRESOLVED; UNRESOLVED is not a soft pass):
- U=6/10 cells: C1 gravity share ≤0.50, C1b (capillary share ≥0.15 AND σ²_h ≥0.002), totals envelope, B-envelope — must PASS (arms RAW/B identical here; arm D same numbers).
- U=1.5/3 cells: at best UNRESOLVED. FAIL if totals outside [0.5×, 2.0×] of the [Cox–Munk, Hu] envelope, any B(k) > 0.012, negative energy, a step in S(k), or |ΔH_s| > 35%.
- Arm B sensitivity: floor sweep {0.0015, 0.003, 0.006} at the U=3 cells (σ²_h delivered + shares per floor).
- Decision data: for each arm, the count of invented inputs acting per cell and every independent-constraint result, tabulated side by side. DO NOT decide the winner yourself — report the table; the merge decision follows PIPELINE-V12's pre-registered rule and is taken outside this session. Set NO default change in this job: spectrum_source default stays as v13 left it.
Transport: windy cell (U=6/F=1000/depth 1.0) full stats per arm (expect identical for RAW/B; report anyway); PLUS U=3/F=1000 stills: arm RAW vs arm B rendered through the frozen pipeline (same seed/frame) → comp-u3-raw-vs-priorB.png side by side, with C stats (mean must be 1.0, p99, p99.9, max) and own-shift for each.
Regression: v13 legacy lane byte-identical re-check; v12 water path untouched (import check only).

## Deliverables (PC offline/out/v14/)
- REPORT-V14-GATES.md — per-arm per-cell tables (U, u*, u*/u*_fit, raw/applied α_m, σ²_g, σ²_h, total MSS, gravity share 50–200, capillary share 50–200, >370 capillary fraction, H_s predicted vs fetch law, B max, B_l band max, three-state verdict per arm with the arm's rule), construction-gate results, sensitivity table, transport stats, the decision table (constraints × arms), regression results. Deviations section mandatory.
- numbers-v14-{spectrum,sweep,transport}.json — every number in the report, machine-readable, exported spectra sufficient to reproduce every band integral independently.
- comp-u3-raw-vs-priorB.png + the windy still.

## Rules
Numerical image statistics only (no vision available). ASCII output. If a number cannot be computed, report it as missing with the reason — never estimate silently. H_s is always spectrum-integrated, never rescaled.
