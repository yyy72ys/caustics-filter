# BRIEF proto4-B2: COMPLETION ONLY for the proto4 detail polish stage

You are running headless on the user's PC. A previous session executing `D:\mutaro-work\rally-in\brief-proto4-details.md` in `D:\yohey2026\dsh\muse\caustics-filter\blender-moth-proto4\` died on a transport socket error near the end. Do NOT redo completed work. Finish only what is missing, exactly as specified below, then stop.

## Known state (verify, do not assume)
- `checkpoint.json` in blender-moth-proto4 shows units done: copy, colour, thorax, shadow. Units silhouette / stills / report are NOT marked done.
- The previous session left try renders: base-f0/ns/sil.png (proto3-state baselines), c1-try1/2/3/3b/3c.png + c1-final.png (colour), c2-thorax.png (thorax), c3-f0/ns.png, c3b-f0.png, c3c-f0.png (shadow), s4a/s4b/s4c/s4d-sil.png (silhouette tries), a `tools\` folder with its measurement scripts (incl. shape.py), and an updated build_moth.py (mtime 20:35).
- The last silhouette measurement visible from that session (tryD): A1 1.457, A3 7 designed / 5 clear peaks, A2 fan_best 0.777 / render_best 0.776, photo_mirror 0.901 by its method.

## Tasks (in order)
1. Re-read `D:\mutaro-work\rally-in\brief-proto4-details.md` for the required report contents and measurement protocols. Work only in blender-moth-proto4.
2. SILHOUETTE VERIFY: Render a fresh silhouette still from the CURRENT build_moth.py (`python run_blender.py build_moth.py -- still s4verify-sil.png --silhouette`) and measure it with tools\shape.py. Compare with the tryD numbers above. If they match (within measurement noise), the current code is the final silhouette state — say so. If they do not match, check whether the try definitions survive in tools\ or any notes/scripts in the folder; if the best try's parameter values are recoverable, apply them and re-verify; if NOT recoverable, keep the current code state and record that honestly in the report. Do NOT start new silhouette exploration tries in this session.
3. FINAL STILLS: Render `stageB-f0.png` (normal), `stageB-obl.png` (--oblique 35), `stageB-ns.png` (--noshadow), `stageB-sil.png` (--silhouette) from the current build_moth.py.
4. MEASUREMENTS FOR THE REPORT: Gather the before/after numbers required by the original brief. Reuse the previous session's measurement outputs wherever they exist in the folder (tools\ outputs, json/txt logs, _dbg\). If a specific number is missing, re-measure it from the EXISTING try renders (c1-final.png for colour, c2-thorax.png for thorax, stageB-f0/stageB-ns for the E3 shadow bands, stageB-sil.png for A1/A2/A3/A4) rather than re-rendering tries. Required set:
   - Colour: reference-photo wing stats (median + p25/p75) if already computed, else compute from D:\mutaro-work\moth-on-bark.jpg wing crops now; rendered wing stats for proto3 baseline (stage3b-f0.png equivalent = base-f0.png) and for proto4 final (stageB-f0.png); D1 micro-contrast before/after (proto3 value 4.14 by its method; measure proto4 the same way); final material parameter values (from build_moth.py).
   - Thorax: silhouette bounding aspect of the thorax blob before (base-f0.png) / after (stageB-f0.png), parameters changed (from the build_moth.py diff vs proto3's build_moth.py — a plain text diff summary is fine), and the C1 CHECK: crop-zoom the centre seam area of stageB-f0.png and state whether any dark abdomen triangle shows between the fan inner margins (must be NO).
   - Shadow: E3 band table (0-3 / 3-6 / 6-12 / 12-24 px mean luminance drop) for proto3 (from base-f0.png vs base-ns.png) and proto4 (stageB-f0.png vs stageB-ns.png), same differencing method both times; parameters changed.
   - Silhouette: A1 / A2 (fan + render, raw + best) / A3 / A4 for proto4 final, next to proto3's values (A1 1.494; A2 fan 0.771, render raw 0.739 / best 0.794; A3 7 designed / 5 clear; A4 0.494).
5. REPORT: Write `stageB-results.md` in blender-moth-proto4 with one section per change (colour / thorax / shadow / silhouette): what changed (parameters), before/after measurement tables, checks (C1, band readability), and honest remaining issues. Append a "proto4 stage B (details)" section to `STAGE-LOG.md` in blender-moth-proto4 summarising renders + results + the fact that this stage was finished by a completion session after a transport failure. Mark checkpoint units silhouette, stills, report done via the jobkit checkpoint (D:\mutaro-work\jobkit\checkpoint.py).

## Rules
- Completion only: do not re-tune colour/thorax/shadow, do not touch choreography/keyframes, do not render the 48-frame video.
- If a required earlier measurement is unrecoverable AND cannot be re-measured from existing renders, write that plainly in the report instead of inventing it.
- Run to completion; never wait for input.
