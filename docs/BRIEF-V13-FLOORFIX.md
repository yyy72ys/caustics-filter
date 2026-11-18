# BRIEF-V13-FLOORFIX: chunk 9 = v3.9, the low-wind capillary floor

Design: PIPELINE-V11.md (read it first, in D:\mutaro-work\rally-in\). Research:
research_notes/low-wind-capillary-floor-20261111 (VM-side; the key numbers are
in the design and in out/v11/numbers-v11-spectrum.json). Rally digests:
rally-v39-{gpt,grok,gemini}.md (copied alongside this brief).

## Task
Implement the design exactly:
- New family v13_spectrum.py / v13_state.py / run_v13.py, importing the v11
  family unchanged (never edit v5–v12 files). Write ONLY in offline/out/v13/
  and the v13 files. Another job is writing out/v12 — do not touch it.
- Remove the hard zero: inside the ECKV fit domain (u* >= c_m/e) alpha_m is
  used raw. Below it, alpha_m_used = max(alpha_m_raw, alpha_m_floor(cell))
  where alpha_m_floor solves sigma2_h,total = 0.003 for the cell. States built
  with the prior active carry the design's `low_wind_prior` label verbatim
  (target, source, rally note). The `low_wind_floor` label is retired.
- Sensitivity sweep (mandatory): sigma2_h target in {0, 0.0015, 0.003, 0.006}
  on the preset-1 scene and the G10 windy scene — C percentiles, own-field
  warp, fine texture |dL/dy|, total sigma2 per level.

## Stages (jobkit checkpoints in offline/out/v13/, one unit per stage)
1. v13 spectrum/state + A1/A2 re-run (algebra vs v11 conventions, Parseval).
2. Eight-cell sweep: C1 three-state (PASS/FAIL in-domain, UNRESOLVED below,
   reported with numbers), C1b floor-proof (all cells), V1–V3 reports, total
   sigma2 beside Cox–Munk clean + Hu per cell, full band partitions, K1/K5/K6.
3. K3 windy compatibility (floored source vs v11 legacy, attribution table as
   in REPORT-V11 §3.2) + K4 regression battery (bit-identity where the prior
   is inactive; legacy path untouched) + the sensitivity sweep renders/numbers.
4. Merge evaluation per PIPELINE-V11's rule + REPORT-V13-GATES.md +
   numbers-v13-*.json. Flip the default ONLY if the rule says so; if not, the
   report states which item blocked it. No blends, no rescaling, no third
   option, no gate redefinition beyond the design.

## Rules
Foreground execution; do not end the session until stage-4 deliverables are
written. NEVER taskkill python.exe/blender.exe. Numerical image statistics
only (no vision). ASCII output, atomic writes. Failures, deviations and
UNRESOLVED items recorded openly — an UNRESOLVED cell is a result, not an
embarrassment to hide.
