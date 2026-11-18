# BRIEF — one still: v3.4 pipeline, gentle-stream environment numbers

One deliverable still image. Work on the PC in
`D:\yohey2026\dsh\muse\caustics-filter\offline\` (the v8/v3.4 stack:
run_v8.py and friends). Model: deepseek-v4.1-flash. Foreground
execution, jobkit checkpoints per step, numerical verification only.

## What to make
The piece's shipped base frame `blender-piece8s\renders_v2\frame_0000.png`
(832x480, moth on bark) composited with the v3.4 caustics pipeline at
**eta = 1.0** (direct sun only — the user's decision today; do NOT use
the v9/Perez diffuse route), composite frozen exactly as production:
A2.6-Blow + E2.2, master 2048^2, Q=25.

The environment is NOT preset 1's pond. It is a **gentle stream edge in
the morning** — same sun/sky family as preset 1, but the water is a
stream's. Set the recipe's water parameters to honest stream numbers:
- depth z = 0.30 m (shallow stream edge)
- wind = 2.5 m/s, blowing along the stream corridor
- fetch = 1000 m along the flow (a stream corridor, not a pond's 100 m)
- H_s must come out of the pipeline's own SMB/fetch computation, never
  hand-set; report the realised value
- directional distribution: narrow it toward the flow axis as far as
  the spectrum code honestly supports; report the parameter used
- rain = 0, surface film = none

## Honesty rules
- If a "stream" property has no slot in the current code (e.g. current
  advection / Doppler steepening of the spectrum), do NOT hack it in
  silently. Use the closest honest slots (the numbers above) and record
  the gap in the report in one line: what a still cannot show anyway
  (advection is motion) vs what is genuinely missing (steepening).
- Record every number you actually used and every realised value
  (H_s target vs realised, slope variance sigma^2, C stats).

## Outputs (leave them in `offline\out\stream-still\` on the PC)
1. `comp-stream-still.png` — the composited still, 832x480
2. `numbers-stream.json` — recipe used + realised H_s, sigma^2, C
   p50/p99/p99.9/max, crush/blow fractions
3. `REPORT-STREAM.md` — short: the numbers table; a side-by-side of
   what changed vs preset 1 pond (fetch 100 m / wind 1.5 / z 0.6) in
   realised values; the honest-gap lines
The coordinator fetches the files; do not try to push them to the VM.

## Gates before you finish
- Output is 832x480, opens, and is NOT bit-identical to the base frame
- C mean = 1.0 (within tolerance), crush < 2%, blow = 0
- The pond-vs-stream comparison in the report uses REALISED values
  from this run, not preset labels
