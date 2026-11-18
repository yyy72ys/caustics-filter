# BRIEF — one still: v3.4, windy pond (physics-first maximum fine warp)

One deliverable still. PC, `D:\yohey2026\dsh\muse\caustics-filter\offline\`
(v8/v3.4 stack). Model: deepseek-v4.1-flash. Foreground, jobkit
checkpoints, numerical verification only.

## Scene (all numbers are the environment; nothing is tuned for looks)
- Base image: the moth photo `D:\mutaro-work\moth-on-bark.jpg`
  (736x540) — the SAME photo as yesterday's production-loop frame 1,
  so the outputs are directly comparable.
- Pipeline: v3.4 exactly as chunk 4 (run_v8), eta = 1.0, master
  2048^2, Q=25, composite frozen A2.6-Blow + E2.2.
- Environment: preset-1 pond, with ONLY these changes:
  - wind = 6.0 m/s (a breezy day; fetch stays the pond's 100 m)
  - depth z = 1.0 m
  - rain = 0 (rain is explicitly OUT of scope this round)
  - surface film: none
- H_s from the pipeline's own SMB computation; all caps (Stokes,
  depth) must hold and be reported.

## Why (context, not instructions to tune)
The target is yesterday-level refraction warp reached honestly:
fine-scale slope variance (wind chop) x depth. Expected realised
sigma^2 ~ 0.04-0.05. If the realised sigma^2 lands below 0.035,
re-check the scene computation for errors and report what you find —
do NOT raise parameters beyond the recipe above.

## Outputs (leave in `offline\out\windy-still\` on the PC)
1. `comp-windy-still.png` (736x540)
2. `numbers-windy.json`: recipe, realised H_s, sigma^2, C stats,
   crush/blow, AND the warp measurement: block-matching displacement
   between the base photo and the composite (48 px blocks, search
   +/-24 px at 2 px steps, luminance): mean shift px, shift-field
   std px. (Yesterday's yardstick, same method: mean 23.2 px,
   std 17.1 px. Today's stream still: 9.2 / 9.5.)
3. `REPORT-WINDY.md`: the numbers + one line each on what the depth
   cost (colour attenuation at z=1.0) did to the image stats.
Coordinator fetches; do not push to the VM.
