# BRIEF — det 堆積の局所化（安定 counting sort でレイをビン順化）

- Job: `caustics-det-binsort`
- Source of this brief: `delivery/det-speedup/REPORT-DET-SPEEDUP.md` §7 item 5 (the recorded next lever of the det_frame speedup job, cf-0051, left unissued until now). This job is a SPEED job with a hard equivalence gate. It does NOT wire anything into production — wiring is a separate follow-up after a PASS.
- Inputs to read on the PC: `D:\yohey2026\dsh\muse\caustics-filter\offline\det_fast.py` (the current production det, wired in by cf-0056), `D:\yohey2026\dsh\muse\caustics-filter\offline\det.py` (the original production det), `D:\yohey2026\dsh\muse\caustics-filter\offline\out\det-speedup\` (the speedup job's report, numbers, gate script `gate_det.py`, and bench scripts). If any of these is missing on the PC, say so in the report and stop — do not invent replacements.

## Background (measured, from the speedup report)

- Production det (2048², Q=25 sun directions) costs 10.96 s/frame effective (3 workers × 4 threads) after the 4× speedup; the remaining dominant cost is the **deposit stage: 5.5 s/frame, 32% of one frame's 17.2 s single-process fast time**.
- The deposit is 4 passes × (scatter + commit) with random reads/writes into a 32 MiB buffer. The report's next lever: **a stable counting sort that orders rays by destination bin** localizes the accesses while preserving the addition order — the relative order of rays within one bin is preserved by the stable sort, so the per-bin accumulation order (and hence the bit-exact result) is unchanged.
- Hard constraints carried over from the speedup job (do not "improve" them away):
  - `row_block = 256` is fixed — it defines the production addition grouping.
  - `fmod` stays for exactness — do NOT replace it with `P - floor(P/S)·S` (that can shift by 1 ulp and break bit-exactness).
  - Run with `python -B`, Numba `cache=False` (as the speedup job did). Do not leave fresh `.pyc` files.

## STAGE 1 — Implement the binned deposit as a NEW module

- Create `det_fast2.py` (new file, offline root). It must expose the same `det_frame_fast` entry point signature as `det_fast.py` so the existing gate and bench scripts can call it unchanged.
- Mechanism: before the deposit passes, stably sort the rays by destination bin (counting sort on the bin index); then run the deposit passes in bin order, keeping the per-bin accumulation order identical to the current implementation. The geometry stage is out of scope — do not touch it.
- If the stable sort itself costs more than the locality it buys, say so with measurements and stop at a negative result — a measured "this lever does not pay" is a valid outcome, not a failure to hide.

## STAGE 2 — Equivalence gate (mandatory, PASS required)

- Reuse the speedup job's gate approach (`out/det-speedup/gate_det.py`): compare the candidate's C against the **current production det** (`det_fast` via the `run_prod_v4.py` path) on **frames 0 and 24**.
- Required: **max abs diff = 0.0 (bit-exact)** on both frames, and the final PNG pixels identical. Any nonzero diff = FAIL — report it with the frame, the max diff, and where it first appears. Do not "fix" a diff by loosening the comparison.

## STAGE 3 — Measure

- Per-frame wall time of the deposit stage alone and of the whole det_frame for frames 0 and 24, single process, same machine-load conditions as the speedup job's measurements where possible. Report the deposit time before/after and the whole-frame time before/after, with the same-worker/thread configuration the speedup job used for its headline numbers (3 workers × 4 threads if the machine allows; state the configuration used).
- Also measure the sort cost separately (time spent in the counting sort per frame) so the net gain is visible.

## Deliverables (new files only, under out/det-binsort/)

- `REPORT-DET-BINSORT.md` — mechanism (how the bin order preserves addition order), the equivalence-gate result (with the diff numbers), the time measurements (deposit and whole-frame, before/after, with configuration), the sort cost, and a verdict: PASS (bit-exact + faster) / FAIL (diff nonzero) / NO-GAIN (bit-exact but not faster, with the numbers). Include a caveats section (what was NOT done: production wiring, 192-frame full run — the same scope limits as the speedup job).
- `numbers-det-binsort.json` — every number in the report, machine-readable.
- `det_fast2.py` — the candidate module (offline root, new file).
- Checkpoint: write `out/det-binsort/checkpoint.json` with progress units as the stages complete.

## Rules

- Do not modify any existing *.py file. `det_fast.py`, `det.py`, `run_prod_v4.py`, `make_video_v4.py` are read-only for this job. Verify with a before/after SHA256 hash list of offline/*.py (same method as the speedup job's §6) and report `n_changed: 0`.
- Do not write to `out/prod/` (read-only, 583 files must stay untouched).
- Do NOT wire the candidate into `run_prod_v4.py` or change `DET_IMPL` defaults. Production wiring is a separate follow-up after a PASS, with its own gate (the det-wire job pattern, cf-0056).
- All numbers in the report must also appear in the JSON. If a number cannot be measured, say so and mark it missing — never fill it with a guess.
