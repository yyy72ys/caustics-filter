# BRIEF-V20-PROD8S-V2COND — the 8 s piece at the v2 loop's place & weather, honest waves

## Why this run exists (owner decision, 2026-11-12)

The 2026-11-10 production loop (the "v2" piece) moved the owner — he could almost hear the water. Its environment is reproduced here as faithfully as the current honest stack allows: same photo, same place and weather (preset 1: a shallow green stream in the morning, depth 0.6 m), same seed. What is NOT reproduced is the old wave law: the v2 waves came from the pre-v3.7 amplitude-curve machinery whose H_s was not spectrum-derived — that defect is why the machinery was replaced, and physics-first forbids reviving it. This run is the controlled comparison: same place, same weather, honest waves.

## STAGE 0 — absolute-path preflight (STOP and report if any item is missing; do NOT search other trees)

- Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
- Base photo: `D:\mutaro-work\moth-on-bark.jpg` — must decode at 735×539.
- Code (all in the workdir): `run_v19_prod.py` (the production pattern to follow), `run_v15_spectrum.py`, `run_v16_water.py` (assert `PHASE_DEFAULT = "ff"`, `KERNEL_DEFAULT = "ff-cdf"`), `run_v16.py`, `v12_water.py`, `v5_composite.py`, `transport.py`, `det.py`, `image.py`.

## Scene (fixed; do not tune)

- Spectrum: v3.10 default (elfouhaily-form, Arm-D rules), **U = 3 m/s, F = 100 m** (light wind over a small stream — the honest mapping of preset 1's place).
- **Depth z = 0.6 m** (the v2 recipe's depthM, exactly).
- Water: **algae_pond** medium + phase (the green-stream character of preset 1), v16 defaults.
- Sky: direct sun only, η = 1.0, sun handling identical to v19. Retina OFF. Rain 0.
- Seed: **3255045276382562994** (same as the v2 loop, v18 and v19).
- C_mesh 2048, det Q = 25, frozen composite A2.6 + Blow (σ 3 px) + E = 2.2 (literal, as in v18/v19).
- Base handling: photo 735×539 → right+bottom 1 px edge-replicated pad to **736×540** (the v2 loop's own convention). Output 736×540, 24 fps, 192 frames, t = n/24 s, frequencies on the 2π/8 s grid (inherited from the state build).

## Gates (same as v18/v19; a failed gate does NOT license tuning — report and stop)

- G1 seam: quantisation proof as in v18 (branch (b) acceptable, branch (a) result printed either way).
- G2 levels: per-frame C mean within 1 ± 0.02 for ≥ 95% of frames; crush < 5% per frame; blowout = 0 per frame.
- G3 encode: ffprobe h264, 736×540, 24 fps, 192 frames, 8.000 s.
- G4 integrity: no v10–v19 file modified (hash the `*.py` set at Stage 0 and at the end; only `run_v20_prod.py` and `out\v20\*` may be new).

## Also report (measured, NOT gated)

- Realised patch size (m), realised σ², spectrum H_s, any warn flags.
- The geometric warp of the delivered piece: mean / p50 / p99 of |displacement| in output pixels at frame 96 — for direct comparison with the v2 loop's yardstick (≈ 23 px mean at its scale, block-matched) and with v19's reported value.

## Deliverables (under `out\v20\`)

- `piece-8s-v20-v2cond.mp4`
- `frame-0096.png`
- `numbers-v20-prod.json`
- `REPORT-V20-PROD8S.md`
- `checkpoint.json` (jobkit; chunks of 24 frames, idempotent per chunk)

## Rules

- Execute in the foreground inside the job; do not end until the MP4 and report are written.
- New code only in `run_v20_prod.py`. Follow the run_v19_prod.py pattern. No edits to any existing file, no parameter invention beyond this brief.
