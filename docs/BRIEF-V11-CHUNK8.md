# BRIEF-V11-CHUNK8 — implement version 3.8 (gravity-tail termination)

Design: `PIPELINE-V10.md` (read it in full first; it is the specification).
Machine: PC, `D:\yohey2026\dsh\muse\caustics-filter\offline`. Model: deepseek-v4.1-flash.
New code family: run_v11* (`v11_spectrum.py`, `v11_state.py`, `run_v11.py`,
`v11_report.py`). v5–v10 modules are imported, NEVER edited.

Run all four stages FOREGROUND in this one session, jobkit-checkpointed
(`out/v11/checkpoint.json`, one unit per delivered artefact, atomic commits).
Do not end the session until `out/v11/REPORT-V11-GATES.md` is written.

## Stage 1 — spectrum module + algebra gates (design §3, §5-A)
- Implement B_l exactly as designed (formulas in PIPELINE-V10 §3; parameter sources:
  k_p/c_p/Ω_c from the existing fetch laws through the existing finite-depth
  gravity–capillary dispersion; NO ECKV fetch relation; γ log-formula continued
  past Ω_c = 5 with the `wave_age_extrapolation` domain label).
- A1: numerically compute B_g from the CURRENT v10 code path and compare pointwise
  with the closed form, deep water and h = 0.3/0.5/1.0. If the finite-depth plateau
  in 50–200 moves >= 30% from the deep-water value: STOP, write the finding into
  the stage-1 numbers file, and end the session with the report stating D1 is
  weakened. Do not proceed to stage 2 on a failed A1.
- A2: Parseval on a built test field from the new spectrum (v10 G0 analogue,
  same test configuration as v10's G0 where applicable).

## Stage 2 — 8-cell sweep (design §5-C, §5-V, K1, K5, K6)
Same 8 cells as v10 (U in {1.5, 3, 6, 10} x F in {100, 1000}, h = 1.0, production
seed family). Per cell, for the new source: band table split by part; gravity
share of 50–200; max B in 50–200 split by part; total sigma^2 vs the v10 totals,
vs full-ECKV at matched inputs, vs clean Cox–Munk; >370 fraction vs ECKV's own
fraction at matched inputs; B_l-predicted H_s vs fetch-law H_s (K1 delta + flag
at >35%); kh_p; domain labels; K5 quantisation table; K6 continuity audit on the
assembled B(k). Legacy-source numbers come from the v10 sweep JSON (do not
recompute the legacy sweep).

## Stage 3 — compatibility (design §5-K2, K3, K4)
- K2: height quantiles of built fields, old vs new source, two reference scenes
  (v10 stage-3 scenes A and B if identifiable from the v10 outputs; otherwise
  U=1.5/F=100/h=0.6 and U=6/F=100/h=1.0), same seed family.
- K3: the v10 G10 windy scene rendered with BOTH sources through the unchanged
  transport + frozen composite (Q=25, master 2048^2, eta=1.0, moth photo native
  735x539): C percentiles, warp yardsticks (own geometric shift in px and mm +
  block matching, same method as v10 G10), fine-texture mean |dL/dy|, crush/blow.
  Deliver both stills.
- K4: regression battery exactly as v10 G8: frozen composite bit-identical on
  identical buffers; Beer–Lambert identical for identical optical inputs;
  rain R=0 equals the NEW base state exactly.

## Stage 4 — merge evaluation + report
- Evaluate the merge rule (PIPELINE-V10 §5) item by item. If and only if it
  passes with no K1 flag open: flip the default `spectrum_source` to
  `elfouhaily-form` in v11_state.py and record the flip in the report. Otherwise
  leave the default at `jonswap` and state the blocking items.
- Write `out/v11/REPORT-V11-GATES.md` in the v10 report's style: conventions and
  any deviation stated openly, per-gate tables with realised numbers, the preset
  regression stated plainly, recorded limitations (including the parasitic-
  capillary gap and the wave-age extrapolation), deliverables list, reproduce
  commands.

## Session rules (standing, from previous chunks)
- Never `taskkill` python.exe or blender.exe. Numerical image statistics only
  (no image inspection). ASCII outputs. Atomic writes via jobkit.
- Physical gate failures are recorded openly with their numbers (chunk-3.5
  gate-3 precedent) — never tune a constant to make a gate pass. There are no
  free parameters in this design; if you find yourself adding one, stop and
  write it up instead.
- Deliverables on the PC in `offline/out/v11/`: REPORT-V11-GATES.md,
  numbers-v11-sweep.json, numbers-v11-transport.json, comp-windy-v38-new.png,
  comp-windy-v38-legacy.png (K3 stills), spectrum-curves-v11.json.
