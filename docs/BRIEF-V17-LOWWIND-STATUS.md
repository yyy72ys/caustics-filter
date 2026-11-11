# BRIEF-V17-LOWWIND-STATUS — settle the low-wind status in code (GPT wording, user decision 2026-11-12)

## Background (one paragraph)

Rally v42 (Gemini + GPT; Grok skipped by the user's one-time exception) ruled on adopting the low-wind capillary prior (σ²_h = 0.003, provenance `low_wind_prior`) as "ADOPTED_WITH_BASIS". GPT REJECTED the adoption claim; Gemini allowed it only under a rename whose headline leads with UNRESOLVED. The user chose GPT's wording (2026-11-12). This job installs that settled status in code as a NEW module plus an automated provenance test. **No existing file may be modified** (the v11–v16 families are frozen and hash-verified).

## STAGE 0 — absolute-path preflight (STOP if any is missing; do NOT search other trees)

- Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
- Must exist: `run_v15_spectrum.py` in that workdir (read it to find the cell-evaluation API and how provenance/`low_wind_prior` appears in cell records). Also confirm `out\v15\numbers-v15-merge.json` exists (the low-wind cell records the test compares against).
- If `run_v15_spectrum.py` has no importable cell-evaluation entry point, STOP and report what its public functions are instead of guessing.

## STAGE 1 — create `lowwind_status.py` (new file in the workdir)

Contents:

1. `STATUS_ID = "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED"`
2. `STATUS_HEADING = "UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a gate pass)"`
3. `STATUS_TEXT` — exactly this wording (the rally-settled text):

> Status: UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a gate pass). For winds below the Elfouhaily et al. (1997) α_m fit domain (u* < c_m/e ≈ 0.085 m/s; U ≲ 2.7 m/s in this pipeline), the published logarithmic fit crosses zero and does not prescribe a capillary-spectrum level. The renderer retains its existing operational convention: α_m is floored per cell to give capillary slope variance σ²_h = 0.003, retaining the ECKV B_h shape. This is a modelling prior, not a measured low-wind capillary spectrum; results using affected cells must report the sensitivity interval σ²_h ∈ [0.0015, 0.006].
>
> Basis and limitation: 0.003 is borrowed as a magnitude anchor from the Cox & Munk (1954) crosswind zero-wind MSS intercept. The mechanism behind that residual is not claimed to operate on a pond, and the intercept does not directly measure the capillary-band component used here. The measured low-wind total-MSS ranges from Cox–Munk and Hu et al. provide broad context for total surface roughness; they do not independently bracket σ²_h. No published measurement currently resolves the capillary-band slope spectrum of a pond at U = 1–3 m/s.
>
> Interpretation of the zero alternative: Measurements of low-wind wave inception and nonzero surface roughness argue against treating the entire water surface as perfectly smooth. They do not establish a nonzero sustained capillary-band variance of 0.003, or rule out a very small capillary contribution under sheltered, film-covered, or below-threshold conditions.
>
> Known limits: The floor is wind- and fetch-independent where active, while real capillary levels may depend on wind, fetch, shelter, surface films and growth thresholds. The small-scale MSS measured on a sheltered lake demonstrates that substantially smoother conditions occur. Band-share gate C1 at these winds remains not passed and is printed, not waived. Every state built with the floor active carries the provenance label low_wind_prior.
>
> Supersession: This operational convention must be reconsidered when a published slope spectrum at U = 1–3 m/s resolves k > 200 rad/m and distinguishes the capillary contribution from total MSS. A faithful implementation of Kudryavtsev et al. (2013) may trigger a model-comparison review, but does not by itself count as independent observational validation.

4. `STATUS_DICT` — machine-readable:
```python
{
  "status_id": "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED",
  "sigma2_h": 0.003,
  "sensitivity_interval": [0.0015, 0.006],
  "interval_kind": "declared sensitivity interval (not an evidence bracket)",
  "provenance_label": "low_wind_prior",
  "gate_C1_low_wind": "NOT PASSED (printed, not waived)",
  "supersession": "published slope spectrum at U = 1-3 m/s resolving k > 200 rad/m that distinguishes the capillary contribution from total MSS",
}
```
5. `def status_block() -> dict` returning STATUS_DICT plus `"heading": STATUS_HEADING` and `"text": STATUS_TEXT`.
6. The string `ADOPTED_WITH_BASIS` must NOT appear anywhere in the module.

## STAGE 2 — create and run `test_lowwind_provenance.py` (new file in the workdir)

The test must:
1. Import `lowwind_status` and assert: STATUS_TEXT starts with "Status: UNRESOLVED"; STATUS_TEXT contains "not a gate pass"; "ADOPTED" does not appear in STATUS_TEXT or STATUS_HEADING.
2. Evaluate the four low-wind cells (U ∈ {1.5, 3.0} × F ∈ {100, 1000}) through the v15 spectrum API (the same entry point whose numbers are stored in `out\v15\numbers-v15-merge.json`; compare against that file where the API shape allows) and assert for every low-wind cell record:
   - the record's provenance shows the low-wind prior was active (the `low_wind_prior` provenance, however the v15 family spells it — find it from the code/records, do not invent a field name; if the stored records carry it under a different key, use that key and say so in the report);
   - the cell's state is NOT an unqualified PASS (it must be UNRESOLVED or carry the unresolved label, exactly as the v15 report recorded).
3. Print PASS/FAIL per assertion and a final verdict; exit non-zero on any failure.
4. Write the results to `out\v17\numbers-v17-lowwind.json` (per-cell: U, F, state, provenance keys found, assertions).

## STAGE 3 — no-modification guard

- At Stage 0, snapshot sha256 of every `*.py` directly in the workdir. At the end, re-hash and assert: the ONLY new/changed files are `lowwind_status.py`, `test_lowwind_provenance.py`, and files under `out\v17\`. Report the comparison.

## Deliverables (all under `out\v17\` except the two new code files)

- `lowwind_status.py`, `test_lowwind_provenance.py` (workdir root)
- `out\v17\REPORT-V17-LOWWIND-STATUS.md` — what was created, the test output verbatim, the guard result, and any place where the v15 family's actual field names differed from this brief's expectations.
- `out\v17\numbers-v17-lowwind.json`
- `out\v17\checkpoint.json` (jobkit: one unit per stage, committed atomically per stage)

## Rules

- Execute in the foreground inside the job; do not end until all deliverables are written.
- No tuning, no extra features, no edits to existing files. If the v15 API cannot be called as expected, STOP after Stage 0/1 and report precisely what blocked Stage 2 (the module + docs still count as delivered; the test outcome is then reported as BLOCKED with the reason, not silently passed).
