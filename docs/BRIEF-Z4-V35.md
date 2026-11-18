# BRIEF: deep z=4.0 condition through the v3.5-state pipeline

## Goal
The user saw the v3.1 depth-pair right panel (offline/out/v5/pair2b-depth-z1.5-vs-z4.png,
right half: "depth z=4.00 m, sand, wind mode, two-way beer 0.02/0.26/0.36, blur 22.0 px")
and liked its light. Reproduce THAT CONDITION with the pipeline as it stood at v3.5:
run_v9.py (accumulates chunk 3.1 water/bottom + 3.2 rain + 3.4 spectrum + 3.5 sky).
Chunk 3.3 retina post-stage: OFF. Chunk 3.6 is design-only everywhere; do not invent it.

## Condition (match the v3.1 pair2b right panel's environment)
- Base photo: D:\mutaro-work\moth-on-bark.jpg (735x539 native).
- Base recipe: the v3.1 preset-1 baseline (see offline/out/v5/numbers-v5.json, its
  "recipe" block; seed 3255045276382562994).
- Overrides: depthM = 4.0; substrate = sand; wind mode = preset wind (as in the baseline);
  eta = 1.0 (direct sun only, matching the pre-3.5 light); rain R = 0; retina OFF.
- Use run_v9.py and its normal recipe/state/transport/composite path (frozen composite
  A2.6 + Blow + E2.2). Do NOT pixel-match the old v5 render; the environment is the target.
  Record the exact recipe JSON you used.

## Rules
- Foreground execution; do not end the session until the deliverables are written.
- jobkit checkpoints in offline/out/z4-v35/ (one per stage).
- NEVER taskkill python.exe/blender.exe. Another job is writing offline/out/v11 — do not
  touch that folder; write ONLY in offline/out/z4-v35/.
- Numerical image statistics only (no vision). ASCII output. Atomic writes.
- If run_v9.py cannot express depth 4.0 / sand / eta through its recipe, say so in the
  numbers file and use the closest supported path; record any deviation openly.

## Deliverables (offline/out/z4-v35/)
- still-z4-v35.png (735x539 composite still)
- numbers-z4.json: recipe used, realised state (Hs, sigma2, patch), C percentiles,
  crush/blow rates, deviations if any.
