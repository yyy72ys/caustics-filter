# BRIEF — 1-second composite samples: eta=1.0 vs eta=0.8 (user-ordered comparison)

PC repo: D:\yohey2026\dsh\muse\caustics-filter\ . Use jobkit checkpoints
(D:\mutaro-work\jobkit\checkpoint.py) per stage. Foreground execution: do NOT end
the session until all deliverable files are written and verified. NEVER
taskkill python.exe (parallel work runs on this PC).

## Purpose
The user chooses the composite light mix for the 8 s main piece by comparing two
1-second video samples over the SAME base frames: (A) eta=1.0 (direct sun only,
the v3.4 route) and (B) eta=0.8 (direct + Perez diffuse, the v3.5 route).

## Base frames
The moth base render lives in blender-piece8s\ (frames f0000..f0191 PNG). Use
frames f0062..f0085 (24 frames: wing-raise -> liftoff -> exit start) as the
sample window. If the PNG frames do not exist, extract the same frames from the
delivery mp4 and record that fallback in the report.

## Per-frame computation (both versions share everything except the light mix)
- Water state: preset 1 exactly as in production (run_v8 / run_v9 code paths,
  same seed, same 8 s loop quantisation), evaluated at each frame's time.
- Version A: C = direct transport only (eta=1) via the v8 route.
- Version B: C_total at eta=0.8 via the v9 route (Perez sky, 96 dome directions,
  diffuse master at 2048^2). FIRST measure the per-frame wall time on frame
  f0062; if the projected total for the 24 B-frames exceeds ~40 min, reduce the
  dome direction count stepwise (96 -> 64 -> 48) to fit ~40 min and RECORD the
  count actually used in the report (the eta comparison is robust within the
  measured convergence band; hiding a reduction is not acceptable, recording it
  is).
- Composite: the FROZEN composite, unchanged (ambient x2.6 + low-pass sigma 3 px
  blow + exposure E=2.2 on the direct term; Beer-Lambert on the direct path as
  in production/chunk 5). Resize/cover to the base frame size 832x480 exactly
  as the production composite does.

## Deliverables (copy to the VM folder delivery\offline-v9\samples-1s\ via the
project's usual transfer, AND keep PC copies under offline\out\v9samples\)
- comp-eta100-1s.mp4 and comp-eta080-1s.mp4 (24 fps, 24 frames, 1.000 s,
  832x480, h264_nvenc is fine)
- still-eta100.png and still-eta080.png (the same mid-window frame, f0074)
- SAMPLE-REPORT.md: per-version wall times, dome direction count used for B,
  mean crush/blow fractions per version, any fallback used.

## Gates (verify and record)
- Both mp4s: ffprobe = 24 frames, 1.000 s, 832x480.
- The two stills differ (they must not be byte-identical); report mean abs diff.
- A-version frame f0062's C statistics match the v8 preset-1 reference within
  tolerance (sanity that the production route was used).
