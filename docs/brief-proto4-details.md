# BRIEF proto4-B: proto4 detail polish stage (colour / thorax / contact shadow / silhouette)

You are running headless on the user's PC. Execute this brief exactly, to completion, without stopping to ask questions. If a step fails, log it and continue. Change ONE element at a time, in the order given, rendering and measuring after each change. All work is code (Blender Python); NO mesh-generation AI, NO photo cutout sprites.

## Base and working copy
- Base: `D:\yohey2026\dsh\muse\caustics-filter\blender-moth-proto3\` (build_moth.py, run_blender.py, assemble.py, STAGE-LOG.md, RESULT-PROTO3.md). Do NOT modify proto3.
- Create a full copy as `D:\yohey2026\dsh\muse\caustics-filter\blender-moth-proto4\` and work ONLY there. Render route (from proto3): `python run_blender.py build_moth.py -- still <name.png>` with optional flags `--silhouette`, `--oblique 35`, `--noshadow`. Blender 4.2.3 LTS, Eevee 64 samples, motion blur on, 832x480, production camera fixed.
- Use jobkit `D:\mutaro-work\jobkit\checkpoint.py` with `checkpoint.json` in blender-moth-proto4. Units: copy, colour, thorax, shadow, silhouette, stills, report. unit_done each; on start skip units already done.
- KEEP UNCHANGED in this stage: choreography/keyframes, background plate, camera, lights, band design (lunulate modulation), seam patches, frenulum single-fan structure. The motion stage comes later and is NOT part of this brief. Do NOT render the full 48-frame video in this stage; stills only.
- Reference photo (for colour sampling): `D:\mutaro-work\moth-on-bark.jpg`.

## proto3 baseline measurements (your before-values; reuse proto3's own measurement scripts/methods if present in the folder, otherwise re-implement the same method and say so)
- A1 forewing L/W 1.494 (PASS band 1.35-1.55). A2 silhouette IoU: fan-only 0.771, render raw 0.739 / best-aligned 0.794 (gate 0.90, FAIL — the reference photo's own mirror IoU is 0.93, i.e. pose asymmetry caps it). A3 scallops: 7 designed, 5 clear peaks. A4 body/wingspan 0.494.
- D1 grain micro-contrast (flat patch high-frequency sigma): proto3 4.14 (proto2 was 3.50). Implementation: wing UV noise mapping scale (110, 200), grain colour multiply (0.94 + 0.11*g), bump strength 0.15 / distance 0.0015.
- D4 fade mottle: low-frequency noise amplitude max 0.16 (hindwing x1.5 toward #a8c4b2), base 0.90 -> tip 1.00 gradient. Overall render still reads PALE and UNIFORM vs the photo (user + reviewer verdict).
- C3 thorax: radii (0.062, 0.080, 0.050), centre z 0.036; fuzz = bump 0.95 / distance 0.009, sheen 0.75, alpha-eroded x1.13 shell. Still reads as a smooth egg-like blob at full view. C1 (abdomen invisible from above, covered by sunken abdomen + centre seam patches on each fan) PASSES and must not regress.
- E3 contact shadow (differencing vs --noshadow render, mean luminance drop by distance band from wing outline): 0-3px -2.4%, 3-6px -3.4%, 6-12px -6.1%, 12-24px -6.9%, beyond ~0. Broad soft component dominates; the fine edge line is WEAK (proto2's localised outline effect was stronger).

## Change 1 — COLOUR (deepen + mottle toward the reference photo)
1. Sample the reference photo: crop 2-3 wing regions of the moth in moth-on-bark.jpg, compute per-channel median and p25/p75 of the wing pixels. Also measure proto3's rendered wing pixels the same way (from stage3b-f0.png in the proto3 folder, wing mask by the silhouette method or manual crops — state your method). Write both sets of numbers down; they are your target and your before.
2. Adjust mat_wing() so the RENDERED wing colour closes the gap toward the photo: expect to deepen the base green (material base colours, currently sampled around base #6e9c90 / dark #5f826c / pale #83afa4) AND increase low-frequency mottle amplitude (0.16 -> up to ~0.25) AND widen the grain colour multiply (e.g. 0.90 + 0.20*g). These are starting points, not orders: iterate at most 3 tries, measuring the rendered wing median after each try, and keep the try whose rendered median is closest to the photo median without (a) breaking the lunulate bands' readability, (b) adding any specular/iridescent look (the green is matte pigment — keep roughness/sheen regime), (c) darkening so far that the wing reads as a different species (photo p25 is your floor guide).
3. Re-measure D1 with proto3's method (flat patch high-frequency sigma on the rendered f0 wing). Target: at least +30% over proto3's 4.14. If colour changes alone do not reach it, strengthen the grain multiply/bump within the same try budget and say what you changed.
4. Record: photo stats, before/after rendered stats, final parameter values, D1 before/after.

## Change 2 — THORAX (kill the egg-blob read at full view)
Directions (combine as needed, keep it one change-set): lower and broaden the thorax (reduce z radius and/or centre height slightly, widen x/y slightly), add a small SEEDED displacement noise to the thorax mesh surface so its outline is not a smooth ellipse, and strengthen the fuzz read at distance (shell scale and/or sheen up a step, alpha erosion kept). The dossier picture: a fluffy piliform-scale mass between the wing bases, wider than the abdomen, not a polished oval.
HARD CHECK after rendering f0: C1 must still pass — inspect the centre seam area closely (crop and zoom): NO dark abdomen triangle may show between the fan inner margins, and the seam patches must still cover. If C1 regresses, revert this change-set and log it.
Record: thorax silhouette bounding aspect (width/height of the thorax blob measured on the f0 render) before/after, parameters changed, C1 check result.

## Change 3 — CONTACT SHADOW (strengthen the fine edge line)
The shadow is a rasterised silhouette texture (make_shadow_image) plus GTAO. Add or strengthen a TIGHT dark core along the outline in the shadow texture (e.g. a high-alpha stroke/erode pass covering roughly the first 2 px inside the outline) while keeping the existing blurred soft layer.
Verify with proto3's E3 protocol exactly: render f0 and f0 --noshadow, compute mean luminance drop in the same distance bands (0-3, 3-6, 6-12, 12-24 px). Target: 0-3px band clearly stronger than proto3's -2.4% (aim for -6% or stronger) with the soft falloff preserved. Then eyeball a 200% crop along a wing margin in the normal f0: if the edge line reads as a drawn cartoon outline, back the core alpha off one step and re-measure.
Record: band table proto3 vs proto4, parameters changed, crop verdict.

## Change 4 — SILHOUETTE (symmetric tweaks only)
The fan outline lists (FAN_LE / FAN_TE) are shared by both sides by construction — keep it that way. NEVER warp asymmetrically; the 0.90 IoU gate is NOT to be chased (documented cap ~0.93 from the photo's own asymmetry).
Allowed tweaks, guided by the dossier's reading of the photo: forewing front margin nearly straight to a gentle arc; forewing apex slightly pointed (NOT hooked/falcate); margin curve smoothness. Make small adjustments, re-measuring after each: A1 (must stay within 1.35-1.55), A2 fan IoU and render IoU (proto3's exact method), A3 scallop peak count, A4 ratio. Stop when fan IoU stops improving or A1/A3 start to degrade; keep the best state. Report the numbers honestly even if A2 remains below 0.90.

## Final stills and report
Render final stills into blender-moth-proto4: `stageB-f0.png` (normal), `stageB-obl.png` (--oblique 35), `stageB-ns.png` (--noshadow), `stageB-sil.png` (--silhouette).
Write `stageB-results.md` in blender-moth-proto4 with, per change: what code/parameters changed, before/after measurements (photo colour stats, rendered colour stats, D1, D4 values, thorax aspect + C1 check, E3 band table, A1/A2/A3/A4), and honest remaining issues. Append a proto4 section to `STAGE-LOG.md` in blender-moth-proto4 (copy proto3's log as the starting point if the copy did not already include it) summarising this stage's renders and results. Mark all checkpoint units done.

## Rules
- Work only inside blender-moth-proto4. ASCII log text preferred (the console mangles Japanese); file contents may be English.
- If a change regresses a currently-PASS item (C1, A1, B2 band readability), revert that change and log the revert.
- Run to completion; never wait for input.
