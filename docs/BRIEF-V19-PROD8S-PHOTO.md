# BRIEF-V19-PROD8S-PHOTO — the 8 s piece, photo base, full v3 stack, impact-real conditions

## Why this run exists (owner feedback, 2026-11-12)

The v18 piece (Blender base, U=6/F=1000/depth 1.0, algae_pond) verified perfectly on data but reads visually weak. The owner orders one more run: **the original photo as base** (the still-image moth version, like the 2026-11-10 production loop), the **full current v3 stack**, and **conditions changed before the run** to a real environment chosen for visual impact — the project's own established route: pick a real condition where the fine-warp band is physically large; never break the physics to chase the look.

## STAGE 0 — absolute-path preflight (STOP and report if any item is missing; do NOT search other trees)

- Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
- Base photo: `D:\mutaro-work\moth-on-bark.jpg` — must decode at 735×539 (the project test input).
- Code (all in the workdir): `run_v18_prod.py` (the production pattern to follow), `run_v15_spectrum.py`, `run_v16_water.py` (assert `PHASE_DEFAULT = "ff"`, `KERNEL_DEFAULT = "ff-cdf"`), `run_v16.py`, `v12_water.py`, `v5_composite.py`, `transport.py`, `det.py`, `image.py`.
- Water clarity: the v16/v12 water code must expose a CLEAR medium preset and the clear phase fit (Fournier–Forand Jonasz clear, n = 1.098, μ = 3.37). Find their exact identifiers from the code and print them in the report. If no clear medium preset exists, STOP and report what presets do exist — do not substitute algae_pond silently.

## Scene (fixed; do not tune)

- Spectrum: v3.10 default (elfouhaily-form, Arm-D rules), **U = 6 m/s, F = 100 m** (a small pond in a breeze). This cell carries `warn_extrapolated_wave_age` — the report must print the flag, as always.
- **Depth z = 4.0 m** (deep clear pond — the owner's proven favourite depth for large irregular warp).
- Water: **clear preset** for both the a/b medium and the phase function.
- Sky: direct sun only, η = 1.0. Retina OFF. Rain 0.
- Seed: **3255045276382562994** (same as the v2-era production loop and v18, for comparability).
- C_mesh 2048, det Q = 25, frozen composite A2.6 + Blow (σ 3 px) + E = 2.2 (literal, as in v18).
- Base handling: photo 735×539 → right+bottom 1 px edge-replicated pad to **736×540** (the v2 production loop's own convention). Output 736×540, 24 fps, 192 frames, t = n/24 s, frequencies on the 2π/8 s grid (inherited from the state build, as in v18).

## Gates (same as v18; a failed gate does NOT license tuning — report and stop)

- G1 seam: quantisation proof as in v18 (branch (b) acceptable, branch (a) result printed either way).
- G2 levels: per-frame C mean within 1 ± 0.02 for ≥ 95% of frames; crush < 5% per frame; blowout = 0 per frame.
- G3 encode: ffprobe h264, 736×540, 24 fps, 192 frames, 8.000 s.
- G4 integrity: no v10–v18 file modified (hash the `*.py` set at Stage 0 and at the end; only `run_v19_prod.py` and `out\v19\*` may be new).

## Also report (measured, NOT gated)

- Realised patch size (m), realised σ², spectrum H_s, the cell's warn flag.
- The geometric warp of the delivered piece: mean / p50 / p99 of |displacement| in output pixels, taken from the transport field of a mid frame (frame 96) — the project's warp yardstick, printed for comparison with the historical values (v2-era loop ≈ 23 px mean at its scale; v18 was not measured this way).

## Deliverables (under `out\v19\`)

- `piece-8s-v19-photo.mp4`
- `frame-0096.png`
- `numbers-v19-prod.json` (parameters as actually read/used, per-frame C stats, warp stats, gates, chunk hashes)
- `REPORT-V19-PROD8S.md`
- `checkpoint.json` (jobkit; chunks of 24 frames, idempotent per chunk as in v18)

## Rules

- Execute in the foreground inside the job; do not end until the MP4 and report are written.
- New code only in `run_v19_prod.py`. Follow the run_v18_prod.py pattern. No edits to any existing file, no parameter invention beyond this brief.
