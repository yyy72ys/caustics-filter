# BRIEF-V16-SCATTER — chunk 12: correct the v3.6 phase function against published measurements

## STAGE 0 — absolute-path location preflight (do this FIRST)
Your working directory must be `D:\yohey2026\dsh\muse\caustics-filter\offline`.
Verify these exact files exist: `v12_water.py`, `run_v12.py`, `run_v15_spectrum.py`,
`run_v15.py`, `out\v12\REPORT-V12-GATES.md` (on the VM it is
delivery/offline-v12/REPORT-V12-GATES.md; on the PC use your out\v12 records),
`out\v15\numbers-v15-merge.json`.
If ANY is missing: STOP, write a BLOCKED report to `out\v16\REPORT-V16-SCATTER.md`
stating exactly which path is missing, and end. Do NOT search other trees.
`D:\mutaro-work\caustics-pub` is the PUBLISH copy — it is NOT the pipeline; never work there.
Write ONLY inside `out\v16\` and new `run_v16*` files. Never modify v11–v15 files.
Use jobkit checkpointing (stage per stage) in `out\v16\checkpoint.json`.
Execute in the foreground until all deliverables below are written; do not end early.

## Why (the measurement verdict, 2026-11-12 research)
The full research report is `research_notes/scattering-validation-20261112/report.md`
on the VM; its load-bearing numbers are reproduced here so this brief is self-contained.
- v3.6's Henyey–Greenstein phase function with g = 0.85 (all particle presets) is
  REFUTED by published measurements for the algae/detritus/clear presets:
  Petzold's average-particle phase function has ⟨cos θ⟩ = 0.924 (Solonenko & Mobley 2015)
  and backscatter fraction B = 0.0183 (Mobley et al. 2002). HG-0.85 has B = 0.0361 (2× too
  much), ⟨θ²⟩ = 0.3681 rad² vs the measured-class 0.09–0.21, and puts 23.8% of scattered
  light beyond 30° vs 6–11% measured. The v12 provenance line ("band [0.8, 0.95] from
  Petzold") is not supported by Petzold's published values.
- The Fournier–Forand (FF) fit to Petzold (n = 1.10, μ = 3.5835) gives ⟨cos θ⟩ = 0.930,
  B = 0.0183, ⟨θ²⟩ = 0.1772 rad². On the research's five-clause scorecard (below):
  HG-0.85 scores 0/5, HG-0.924 scores 3/5 (shape clauses fail), FF-Petzold scores 5/5.
  The measurements support changing the phase-function FORM, not just g.
- Shape: measured kernels put 25–42% of scattered light within 1° of forward
  (FF-Petzold: F(<5°) = 0.559); a moment-matched Gaussian puts 4.2% within 5°.
  No g fixes HG's form, and a single Gaussian misplaces the halo energy even with
  the right ⟨θ²⟩. The measurement-supported kernel is the FF angular CDF used
  directly in the displacement integral.
- Sediment preset: NOT resolved either way. HG-0.85 sits between the biogenic and
  mineral FF families on every moment. No published FF fit for a measured
  sediment-dominated freshwater VSF was found. Sediment's phase function must be
  labelled an ASSUMPTION with its envelope printed, not a measurement.
- H (halo energy fraction, ≈0.08 at the reference scene) is phase-function
  independent in this model class: it must NOT move. Any movement = code leak.

## What to build (new family run_v16*, on the v12 water + v15 spectrum base)
1. **FF phase function** β̃(θ; n, μ) in the Fournier & Jonasz (1999) analytic form
   (Ocean Optics Web Book formulation), normalised by quadrature on a log-spaced
   grid θ ∈ [1e-6, π]. Per-preset parameters:
   - `algae_pond` (default medium): FF n = 1.10, μ = 3.5835 (Petzold fit).
   - `clear`: FF n = 1.098, μ = 3.37 (Jonasz clear fit).
   - `sediment`: FF n = 1.15, μ = 4.0 nominal, labelled ASSUMPTION; every sediment
     state prints the envelope [FF(1.10, 3.5835), FF(1.15, 4.0)] and provenance
     flag `phase_assumption`.
   - HG remains selectable as `phase: hg` (legacy lane, g per preset as in v12).
2. **Kernel**: implement the halo kernel from the FF angular CDF used directly in
   the existing displacement integral (replacing the ⟨θ²⟩-moment Gaussian) as the
   candidate default. Keep the Gaussian kernel selectable (`kernel: gaussian`).
3. **Stage A — scorecard** (research §5 Step 1, verbatim targets): by quadrature
   compute ⟨cos θ⟩, B, ⟨θ²⟩, F(<1°), F(<5°), F(<10°), β̃(1°), β̃(5°), β̃(10°), β̃(90°)
   for HG-0.85, HG-0.924, and each preset's FF parameters. Pass clauses for the
   algae/generic preset: ⟨cos θ⟩ = 0.924 ± 0.02; B = 0.0183 ± 0.005;
   ⟨θ²⟩ = 0.177 ± 15%; F(<5°) = 0.559 ± 0.10; β̃(90°) within ×1.5 of the DECLARED
   reference — declare it in the report: Petzold flat-back level 0.00291 sr⁻¹
   (FF-Petzold gives 0.00420 = ×1.44, PASS). Do not mix composites silently;
   also print the Jonasz-turbid-table comparison for the record.
4. **Stage B — halo-level check** (research §5 Step 2): recompute the v12 G1
   reference scene and the four G3 sweep rows through the v15 default stack with
   the FF water. Required: the moment-diagnostic σ scales as √(⟨θ²⟩ ratio) per
   preset; **H bit-identical to the v12 run** (report the comparison exactly);
   delivered stills' halo share H/(T+H) unchanged; peak profile at fixed lags
   (G3 r32/r48 instruments) moves in the predicted direction.
5. **Stage C — kernel-shape check** (research §5 Step 3): encircled-energy of
   the Gaussian kernel vs the FF-CDF kernel at r = 8.5, 32, 79 px on the
   reference scene; report the three differences for both kernels.
6. **Stage D — regression**: (i) `phase: hg` + `kernel: gaussian` reproduces the
   v12 numbers (G1 legacy lane sha cfb086b22a0e where applicable; per-number
   comparison otherwise); (ii) v15 spectrum numbers unchanged (recompute one
   cell and compare to numbers-v15-merge.json); (iii) frozen composite module
   hashes unchanged; (iv) no v11–v15 file modified (hash statement).

## Merge rule (pre-registered; apply it, do not improvise)
The FF-CDF kernel becomes the water default (`PHASE_DEFAULT = ff-cdf`, per-preset
parameters above) IF AND ONLY IF: Stage A passes all five clauses for the
algae/generic preset (β̃(90°) reference declared), Stage B passes (H bit-identical,
σ scaling exact), Stage C is reported, and Stage D is fully green.
Pre-registered expectation for Stage C: the Gaussian's F(<5°) error is ~0.5, far
beyond any plausible tolerance — if Stage C instead showed the Gaussian within
0.05 encircled-energy of the FF kernel at all three radii, retain the Gaussian
and say so. If any other gate fails: write the report with the failure, keep the
v12 default, STOP. No tuning of n/μ to force a pass — the parameters above are
the published fits, not knobs.

## Deliverables (out\v16\)
- REPORT-V16-SCATTER.md — scorecard table (all candidates, all clauses, the
  declared β̃(90°) reference), Stage B/C results, regression item-by-item,
  merge statement (fired or not, and why), sediment assumption labelling proof
  (a sediment state's printed provenance), hash list, explicit "no v11–v15 file
  modified".
- numbers-v16-scorecard.json, numbers-v16-halo.json
- pair-halo-old-vs-new.png — the G3 mid scene (b = 0.25) rendered with the v12
  default (HG-0.85 + Gaussian) beside the merged candidate (FF-CDF), same
  exposure, labelled.
- checkpoint.json (jobkit, one unit per stage A–D + report).
