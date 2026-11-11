# BRIEF-V18-PROD8S — the 8-second piece on the v3.10 + v16 stack (v3 ALPHA close)

## What this is

The closing piece of v3: the shipped Blender moth animation (base v2) composited under the current standard water stack — v3.10 spectrum (elfouhaily-form, Arm-D rules), v16 water (FF phase function, FF-CDF kernel), sky = direct sun only (η = 1.0), frozen composite A2.6 + Blow + E2.2. Output: a seamless 8 s loop MP4 at 832×480, 24 fps, 192 frames.

## STAGE 0 — absolute-path preflight (STOP and report if any item is missing; do NOT search other trees)

- Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
- Base frames: `D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\renders_v2\` — must contain exactly the 192 shipped frames at 832×480 (confirm count and size from the files themselves; the hyphenated `renders-v2` folder is a SUPERSEDED different animation and must NOT be used).
- Code (all in the workdir): `run_v15_spectrum.py`, `run_v15.py` (the merged-stack still path used for the T1–T3 transport stills), `v12_water.py`, `run_v16_water.py` (must carry `PHASE_DEFAULT = "ff"`, `KERNEL_DEFAULT = "ff-cdf"` — assert both), `run_v16.py`, `v5_composite.py` (frozen composite), `transport.py`, `det.py`, `image.py`.
- Scene definition: the v15 T1 transport configuration exactly — U = 6 m/s, F = 1000 m, depth = 1.0 m — as recorded in `out\v15\numbers-v15-transport.json` (read it; if the T1 parameters differ from this line, the FILE wins and the difference is reported). Water preset: the v16 default lane (algae_pond phase parameters). Retina OFF, rain 0, η = 1.0.

## Production design (follow the pipeline's own conventions)

- Each output frame n (0–191) uses base frame n as the submerged scene. The wave field is time-evolved at t = n/24 s with all component frequencies quantised to multiples of 2π/8 s (the v2 production convention), so the water layer is exactly 8 s-periodic and the loop seam carries no water discontinuity (the base animation is already pixel-exact at its own seam).
- Per frame: build the surface from the v15 spectrum at time t (same seed/parameters as the T1 still), transport through the v16 water (a/b partition + FF-CDF halo), composite with the frozen v5 composite (A2.6 + Blow + E2.2). Reuse the run_v15/run_v16 still code paths — do not re-derive physics in this job.
- Work in chunks of 24 frames. Each chunk is committed atomically (frames to a temp name, then rename) and recorded in `out\v18\checkpoint.json` (jobkit). Re-running the job must skip chunks whose 24 frames already exist with the recorded hashes (idempotent per chunk).

## Gates (report pass/fail per gate; a failed gate does NOT license tuning — report and stop)

- G1 seam: the caustics/transport field at t = 8 s equals the field at t = 0 (max abs difference 0 on the field arrays, or the quantisation proof printed from the frequency table).
- G2 levels: per-frame C mean within 1 ± 0.02 for at least 95% of frames (report min/max over frames); crush (< 0.02 linear) < 5% per frame; blowout = 0 per frame.
- G3 encode: ffprobe on the final MP4 shows h264, 832×480, 24 fps, 192 frames, duration 8.000 s (± 1 frame period).
- G4 integrity: no file of the v10–v17 families modified (hash the `*.py` set at Stage 0 and at the end, as in the v17 guard; only `run_v18_prod.py` and `out\v18\*` may be new).

## Deliverables (under `out\v18\`)

- `piece-8s-v3alpha.mp4` — the finished piece.
- `frame-0096.png` — the mid-piece frame, final composite.
- `numbers-v18-prod.json` — scene parameters as actually read, per-frame C mean/crush/blow summary, gate results, chunk hashes.
- `REPORT-V18-PROD8S.md` — parameters, gate results, wall time, anything that deviated from this brief and why.
- `checkpoint.json` (jobkit).

## Rules

- Execute in the foreground inside the job; do not end until the MP4 and report are written.
- New code only in `run_v18_prod.py`. No edits to any existing file. No parameter invention: anything this brief leaves open is taken from the v15 T1 record or the v16 defaults, and the choice is printed in the report.
