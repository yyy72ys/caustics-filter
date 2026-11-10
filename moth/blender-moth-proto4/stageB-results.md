# proto4 stage B — detail polish results (colour / thorax / contact shadow / silhouette)

Base: `blender-moth-proto3` copied verbatim to `blender-moth-proto4` (robocopy /E, 26 files).
proto3 itself was never written to (verified: same 26 files, all timestamps <= 19:52).
All "before" numbers come from renders of the **byte-identical copy** in proto4
(`base-f0.png`, `base-ns.png`, `base-sil.png`, produced before any edit).

Render route for every still:
`python run_blender.py build_moth.py -- still <name>.png` (+ `--silhouette`, `--oblique 35`, `--noshadow`).
Blender 4.2.3 LTS, Eevee NEXT 64 samples, motion blur on, 832x480, production camera.
Stills only — the 48-frame video is out of scope for this stage and was NOT rendered.

---

## 0. Measurement tooling (important caveat)

proto3 shipped **no measurement scripts** (only build/render/assemble code), so every
metric was re-implemented from the brief's definition and **calibrated against the
brief's own proto3 before-values**. Parity:

| metric | brief (proto3) | re-implementation on the proto4 copy of proto3 | parity |
|---|---|---|---|
| A1 forewing L/W | 1.494 | **1.483** | -0.7 % |
| A2 photo own-mirror IoU | 0.93 | **0.901** | -2.7 % |
| A2 fan IoU (best aligned) | 0.771 | **0.775** | +0.5 % |
| A2 render IoU (best aligned) | 0.794 | **0.775** | -2.4 % |
| A3 scallops designed / clear | 7 / 5 | **7 / 5** | exact |
| D1 grain sigma (flat patch) | 4.14 | **3.96** | -4.4 % |
| A4 body/wingspan | 0.494 | **0.345** | **not reproducible** |
| E3 band drops | -2.4 / -3.4 / -6.1 / -6.9 / ~0 % | **-4.3 / -18.7 / -26.6 / -15.4 / -1.4 %** | **not reproducible** |

Two methods could not be recovered:

- **A4**: no reading of "body" / "wingspan" gives 0.494 from this build (body
  0.369 / span 1.07 = 0.345; photo- and render-derived variants all give
  0.33-0.36). Reported with an explicit definition instead: *body length
  (head front y=0.134 .. abdomen tip y=-0.235) / full x span incl. wing-root
  offsets (2*(0.505+0.030))*. Definition unchanged, so before == after.
- **E3**: the brief's absolute magnitudes cannot be reproduced by any distance
  reference I tried (segmented mask, fixed mask, exact projected outline) — the
  curve *shape* matches (peak away from the edge, ~0 past 24 px) but is 3-6x
  larger. proto4 is judged against **my own proto3 baseline measured the same
  way**; the brief's table is quoted alongside.

A2 "raw" (no rotation, area-matched, centroid-aligned) = 0.51 for both render
and fan; the brief's 0.739 raw could not be reproduced. The **best-aligned**
numbers match within 2.4 %, which is what the gate discussion uses.

Tools (kept in `tools/`): `measure.py`, `geomask.py`, `e3.py`, `e3seg.py`,
`final_metrics.py`, `shape.py`, `quick.py`, `thorax.py`, `crop.py`, `sbs.py`,
`gridcrop.py`, `ckpt.py`. Completion-session verification scripts added alongside
them: `e3geo.py` (independent E3/D1 with a fixed projected-outline reference),
`bandcheck.py` (band-readability proxy), `darkaudit.py`, `barkaudit.py`,
`c1audit.py`, `c1shape.py` (centre-seam / thorax-footprint checks),
`basecheck.py` (proto3 copy parity), `colourhistory.py` (colour try table).

---

## Change 1 — COLOUR

### 1.1 Sampling method
- **Photo** (`D:\mutaro-work\moth-on-bark.jpg`, 735x539): moth mask = smoothed
  greenness `G-(R+B)/2 > 8` AND local luminance std (9x9) < 10, 11x11 closing,
  hole fill, largest component. 3 crops of 80x80 on the mask's principal axes
  (all three 100 % pure wing interior, n = 19 200 px).
- **Render**: same segmentation on the perched moth, wing interior = mask eroded
  7 px (n = 32 000 px). Colour-based, i.e. "wing pixels" not "background".

### 1.2 Colour statistics (sRGB 0-255, R,G,B)

| set | median | p25 | p75 | std | n |
|---|---|---|---|---|---|
| **photo (target)** | **112, 157, 148** | 104, 144, 134 | 121, 167, 160 | 14.8, 15.9, 19.9 | 19 200 |
| proto3 render (before) | 90, 122, 108 | 87, 117, 103 | 93, 125, 112 | 6.8, 6.5, 8.1 | 31 971 |
| **proto4 render (after)** | **111, 156, 148** | 106, 149, 140 | 114, 160, 153 | 7.6, 10.7, 16.4 | 32 132 |

Honest reading: against the photo the proto3 wing was **not brighter, it was
~23 % darker and far flatter** (p25-p75 spread 6 levels vs 14; std 6.5 vs 15.9).
The brief says "pale" and "deepen", but its own numeric rule is *"keep the try
whose rendered median is closest to the photo median"* with *"photo p25 as your
floor"*. I followed the numeric rule (lift + slight cool), and attacked the
flatness with mottle + grain. Logged as a deviation from the word "deepen".

### 1.3 Parameter changes (all in `mat_wing()`)
| item | before | after |
|---|---|---|
| ramp dark / base / pale | `#55755f` / `#6e9c90` / `#83afa4` | `#6e9b87` / `#8ecccc` / `#aae8e9` |
| margin darkening | `#5f826c` | `#7cac99` |
| hindwing tint | `#7cab9f` | `#a1e3e2` |
| vein colour | `#59755f` | `#739b87` |
| discal marks (m1, m2) | `#4c6e52` | `#629274` |
| fade-mottle target (D4) | `#a8c4b2` | `#ccded2` |
| **D4 mottle amplitude** | **0.16** (hind x1.5) | **0.25** (hind x1.5 kept) |
| grain colour multiply (D1 source) | `0.94 + 0.11*g` | `0.85 + 0.30*g` |
| **new fine grain octave** | — | Mapping (420,760,1) -> Noise (scale 1, detail 2, rough 0.7) -> `0.935 + 0.13*g2` |
| bump | strength 0.15 / dist 0.0015 | strength 0.26 / dist 0.0018, height = 0.55*g + 0.45*g2 |
| roughness / sheen / specular | 0.74-0.93 / 0.10 / 0.05 | **unchanged** (matte pigment regime kept) |

Iteration log (rendered median; goal 112/157/148):
1. palette x(1.24,1.29,1.37), mottle 0.25, grain 0.90+0.20g -> 107,153,143
2. palette x(1.047,1.026,1.035), grain 0.88+0.24g -> 111,156,148  (**kept**)
3. grain 0.85+0.30g, bump 0.26 -> 111,156,148 (median stable, grain up)

### 1.4 D1 grain micro-contrast
Method: Laplace(4-neighbour) std of Rec.709 luminance, median over 24x24 blocks
fully inside the wing and in the bottom 30 % of low-frequency (sigma 6) std.

| sampling mask | proto3 | proto4 | delta |
|---|---|---|---|
| segmented wing (eroded 7 px) | 3.96 | **4.96** | **+25.3 %** |
| geometric fan outline (no body) | 3.78 | **4.62** | **+22.2 %** |

**Target (+30 % over 4.14 = 5.38) NOT reached.** Evidence that this is metric
saturation at 832x480, not a weak material:
- grain multiply 0.24->0.30 and bump 0.18->0.26: 5.00 -> 5.10;
- bump 0.42 / dist 0.0025: **5.12** (no gain);
- speckle screen 0.08->0.16: **4.58** (loss).
Residual high-frequency content is dominated by Eevee AA/TAA. For reference the
**photo's own D1** (same metric) is **12.1**, so the render is still much
smoother than the reference: "grain +25 % by the metric, visibly less flat,
still smoother than the photo".

### 1.5 Honest remaining issues (colour)
- Photo has ~2x the tonal spread of the render (std 15-20 vs 7-16).
- The body/thorax was deliberately **not** re-coloured (brief lists the wing
  palette only). With the wings now brighter, the thorax is a darker olive
  island; the photo's thorax is brown/tan. This is the most visible colour
  mismatch left.
- The bark plate is fixed and measured darker than the photo's bark
  (render overall mean 70,80,71 vs photo 76,90,80), so matching the photo's
  wing median raises moth-vs-bark contrast above the photo's: at 1x the moths
  read slightly more "cut out" than in proto3. That is a plate/exposure
  interaction, not a wing-colour error, and the plate must stay unchanged.
- Lunulate bands (B2) stay readable at 1x and 3x (checked visually) but are
  lower contrast than in proto3 because the base moved up toward them.

---

## Change 2 — THORAX

### 2.1 Code changes (`build_moth()` thorax block, `mat_fur_shell()`, mats table)
| item | before | after |
|---|---|---|
| thorax radii (x, y, z) | 0.062, 0.080, 0.050 | **0.072, 0.086, 0.042** (wider, flatter) |
| thorax centre z | 0.036 | **0.030** (lower) |
| **new seeded displacement** | — | `r *= 1 + 0.13*lump(u)` on the unit sphere; 3-octave sum-of-sines, seed 4.17, 5-9 lumps across the body |
| fur shell scale | 1.13 | **1.18** |
| fur-shell noise / alpha cut | scale 40 / > 0.45 | **scale 90 / > 0.36** |
| fur-shell sheen | 0 | **0.60** (sheen roughness 0.40) |
| body material sheen | 0.75 | **0.85** |
| alpha erosion of the shell | present | **kept** |

### 2.2 Measurements on the f0 render
Body blob = pixels with `G-B > 14` inside a 130x120 box at the body centre.

| item | proto3 | proto4 |
|---|---|---|
| visible blob area | 1 871 px | **2 276 px** (+21.6 %) |
| bounding box | 57 x 45 px | 59 x 53 px |
| bounding aspect (w/h, image axes) | 1.267 | 1.113 |
| **blob PCA axis ratio (orientation-free)** | **1.505** | **1.358** (10 % squatter) |
| mean body colour | 92,118,86 | 94,121,87 |

The image-space bbox is axis-aligned while the moth is rotated -60 deg, so the
bbox mixes the thorax axes; the PCA ratio is the meaningful width-vs-depth
number and it moves the intended way. The bbox got taller because lowering the
thorax exposes slightly more of the head side.

### 2.3 C1 hard check (centre seam)
- Numeric proxy: dark-pixel fraction in a 36x36 box at the body centre:
  `dark_frac = 0.0` before and after.
- Visual: 4x / 5x centre crops (`_dbg/cmp_c2_centre.png`, `_dbg/cmp_th_final.png`):
  "no dark abdomen triangle visible; C1 PASSES", and the stronger lump
  amplitude (0.085 -> 0.13) was re-checked after rendering.

**C1: PASS, no regression. No revert.** The +16 % x radius helps cover the seam.

### 2.4 Honest remaining issues (thorax)
- At 1x the thorax is still a compact mass; the fuzz reads as a ragged
  silhouette plus sheen, not as individually visible hair.
- In the **oblique** view (`stageB-obl.png`, 35 deg) the body still reads as a
  smooth, slightly glossy capsule: the displacement/alpha fuzz works on the
  silhouette at full view but does not survive at grazing angles. Raising body
  roughness (0.88 -> ~0.95) and adding a second, finer hair shell would be the
  next step; not done here because the brief's Change 2 list stops at shell
  scale / sheen / displacement.
- Because the wings got brighter, the thorax separates by *value* more than by
  texture — still a blob, but no longer a polished egg.

---

## Change 3 — CONTACT SHADOW (with the "back off" iteration the brief asks for)

### 3.1 Code changes (`make_shadow_image()`, `main()`)
Final alpha stack:

```
a = min(0.90, 0.62*sharp + 0.45*wide + 0.85*stroke + 0.22*core
              + 0.22*contact_edge + 0.08*contact_fill)
```
- `stroke`: 4 px / blur r=2 / weight 0.85 -> **5 px / blur r=2 / weight 0.85**.
- **new `core`**: `(box_blur(mask,2) - 0.50)/0.32` clamped = the first ~2 px
  inside the outline (the brief's "high-alpha stroke/erode pass").
- **new `contact_edge` / `contact_fill`**: a second copy of the silhouette drawn
  inside the texture, shifted so that after the decal's rotate/scale/+SH_OFF it
  lands exactly under the wing — this is what puts a line *on* the wing margin
  (the decal itself is offset by SH_OFF = 0.075,-0.075 ~ 24 px, so without it
  the "fine line" sits 24 px from the wing edge).
- SH_OFF is applied **after** the decal's rotation/scale, so the shift is per
  decal: `shift = -R(rot)^-1 * SH_OFF / scale`. `make_shadow_image()` now takes
  `(rot_deg, scale, name)` and `main()` builds **two textures**:
  L rot 50 scale 1.059 -> shift (0.0087, 0.0997); R rot -60 scale 1.157 ->
  shift (-0.0885, -0.0237).
- Soft layers, decal offset, opacity animation, sun and GTAO untouched.

### 3.2 Iteration log (this is where the brief's "back off if it looks drawn" applied)
| build | core / contact / cap | 0-3 px mean (geo) | 1x verdict |
|---|---|---|---|
| proto3 | — / — / 0.88 | -15.26 % | broad dark rim |
| c3 core only (no contact) | 0.45 / 0 / 0.96 | -15.03 % | rim unchanged |
| c3c strong contact | 0.45 / 0.60 / 0.96 | -15.08 % | **hard black halo / "cutout-like", doubling suspected** |
| c3d back off step 1 | 0.30 / 0.32 / 0.92 | — | halo reduced, "still fairly graphic", no doubling |
| **c3f final (kept)** | **0.22 / 0.22 / 0.90** | **-12.87 %** | **"no longer a graphic dark halo… behaves like natural contact darkening; proto4 looks better overall"** |

The strong-contact build was **rejected** exactly as the brief instructs ("if the
edge line reads as a drawn cartoon outline, back the core alpha off one step and
re-measure"). Two steps were needed for the 1x view to read naturally.

### 3.3 E3 band tables (f0 vs f0 --noshadow, frame 0, perched moth)

**(a) segmented mask, held fixed to the proto3 shadow-free render** (proto3's method)
| band | brief (proto3) | proto3 (mine) | **proto4 final** |
|---|---|---|---|
| 0-3 px | -2.4 % | -4.28 % | **-3.43 %** |
| 3-6 px | -3.4 % | -18.69 % | **-16.53 %** |
| 6-12 px | -6.1 % | -26.60 % | **-27.49 %** |
| 12-24 px | -6.9 % | -15.41 % | **-17.53 %** |
| 24-60 px | ~0 | -1.43 % | **-1.74 %** |
| 0-3 px darkest 5 % | — | -50.94 % | **-55.75 %** |

**(b) exact projected fan outline** (`geomask.py`, stable reference)
| band | proto3 | **proto4 final** | change |
|---|---|---|---|
| 0-3 px | -15.26 % | **-12.87 %** | -2.4 pt |
| 3-6 px | -22.83 % | **-23.48 %** | +0.7 pt |
| 6-12 px | -22.49 % | **-24.23 %** | +1.7 pt |
| 12-24 px | -12.42 % | **-14.59 %** | +2.2 pt (+17 %) |
| 24-60 px | -1.01 % | **-1.19 %** | +0.2 pt |
| 0-3 px darkest 5 % | -55.08 % | **-59.26 %** | +4 pt |

### 3.4 Honest reading of E3
- The **mean of the 0-3 px annulus got weaker** (both references), not stronger.
  The only strong version (contact 0.60) reached parity with proto3 (-15.08 vs
  -15.26, i.e. still not stronger) while producing a hard halo — so the brief's
  "+6 % or stronger with a natural look" could not be had with this decal
  geometry.
- What did improve: the 6-24 px bands (+2..+17 %), and the darkest pixels of the
  0-3 px band (-55.1 -> -59.3 %), i.e. the line is tighter, not broader.
- Root cause (worth recording for the motion stage): the shadow is a *cast*
  decal offset 24 px; a genuine contact line needs to live at the wing outline,
  which required the per-decal contact copy. That copy is what creates the halo
  as soon as it is strong enough to move the 0-3 px mean. A real fix would be
  a true AO/contact term (e.g. a zero-offset decal on its own layer with a
  narrow outline), not more alpha in the cast texture.
- The contact copy is clipped at the canvas border for the right decal
  (shift -0.0885 in x takes the mirrored fan to x = -0.62). No visible effect in
  these stills, but a latent limitation.

---

## Change 4 — SILHOUETTE

Fan outline stays **shared and symmetric**; no asymmetric warp. FAN_LE unchanged;
FAN_TE rows t=0.454 / 0.580 / 0.706 deepened by 0.016:

| row | before | after |
|---|---|---|
| (0.454, 0.132) | -0.310 | **-0.326** |
| (0.580, 0.196) | -0.335 | **-0.352** |
| (0.706, 0.235) | -0.270 | **-0.284** |

Tries (each rendered `--silhouette`, measured):

| try | what | A1 | A2 render | A2 fan | A3 clear | verdict |
|---|---|---|---|---|---|---|
| base | proto3 outline | 1.483 | 0.775 | 0.775 | 5 | — |
| A | straighter LE (bulge x0.6) + pointed apex | 1.526 | 0.768 | 0.769 | 5 | **rejected** (A2 worse) |
| B | A's LE + deeper TE | 1.470 | 0.770 | 0.777 | 5 | rejected |
| C | deeper TE only | **1.434** | **0.780** | **0.777** | 5 | **kept** |
| D | deeper TE, half amplitude | 1.457 | 0.776 | 0.777 | 5 | rejected (no gain) |

Final shape metrics:

| metric | proto3 | **proto4** | gate |
|---|---|---|---|
| A1 forewing L/W | 1.483 | **1.434** | PASS 1.35-1.55 |
| A2 fan IoU (best aligned) | 0.775 | **0.777** | FAIL (< 0.90, documented cap) |
| A2 render IoU (best aligned) | 0.775 | **0.780** | FAIL (< 0.90) |
| A2 photo own-mirror IoU | 0.901 | 0.901 | reference / cap |
| A3 designed / peaks / clear | 7 / 7 / 5 | **7 / 7 / 5** | unchanged |
| A4 body/wingspan | 0.345 | **0.345** | definition unchanged |

Honest reading: +0.005 render / +0.002 fan IoU — essentially nothing — while A1
moves away from the photo's 1.494 (still inside the PASS band). The 0.90 gate
stays unreachable: the photo's own mirror IoU is 0.901, and the residual error is
structural (a symmetric two-lobe fan vs the photo's asymmetric wedge with the
body at the tip). Kept because A1/A3 stay in band and A2 is marginally better,
but this is a cosmetic tweak, not a fix.

---

## Final numbers summary

| item | proto3 | proto4 | verdict |
|---|---|---|---|
| wing median (photo 112,157,148) | 90,122,108 | **111,156,148** | matched |
| wing std | 6.8, 6.5, 8.1 | 7.6, 10.7, 16.4 | livelier (photo 15-20) |
| D1 grain | 3.96 (3.78 geo) | **4.96 (4.62 geo)** | +25 % / +22 %, target +30 % MISSED |
| D4 mottle amplitude | 0.16 | **0.25** (hind x1.5) | as specified |
| thorax blob aspect (PCA) | 1.505 | **1.358** | squatter |
| thorax blob area | 1 871 px | **2 276 px** | bigger mass |
| C1 (no abdomen triangle) | PASS | **PASS** | no regression |
| E3 3-6 / 6-12 / 12-24 px | -22.8 / -22.5 / -12.4 % | **-23.5 / -24.2 / -14.6 %** | slightly stronger |
| E3 0-3 px mean | -15.3 % | -12.9 % | **weaker** (trade-off, see 3.4) |
| E3 0-3 px darkest 5 % | -55.1 % | **-59.3 %** | tighter/blacker line |
| A1 | 1.483 | 1.434 | PASS band |
| A2 fan / render IoU | 0.775 / 0.775 | 0.777 / 0.780 | still < 0.90 |
| A3 | 7 / 5 | 7 / 5 | unchanged |
| A4 | 0.345 | 0.345 | unchanged |

Untouched as required: choreography/keyframes, background plate, camera, lights,
lunulate band design, seam patches, frenulum single-fan structure. No video render.

## Remaining issues (all changes)
1. **D1 below target and saturated** — the metric cannot be pushed at 832x480
   with this material; the wing is still smoother than the photo (12.1 vs 4.96).
2. **Body colour** is olive, not the photo's tan, and now much darker than the
   wings: the most visible remaining mismatch.
3. **A2 stays far below 0.90** (+0.005 gained); structural, capped by the
   photo's own asymmetry (0.901 mirror IoU).
4. **E3 0-3 px mean is weaker than proto3** — the natural-looking contact line
   and a stronger 0-3 px mean could not be achieved together with this cast-decal
   geometry; a real contact/AO term is the proper next step.
5. **Full-frame perception**: at 1x the brighter wings separate from the fixed
   (darker-than-photo) bark plate more than in proto3; a re-graded plate would
   fix it, but the plate is frozen in this stage.

---

## Completion-session verification (2nd session)

This stage was finished by a completion session after the first session died on a
transport/socket error near the end. That session re-tuned nothing; it only verified and
finished. What it checked, independently of the first session's own logs:

1. **Checkpoint honoured, not assumed.** On entry `checkpoint.json` had copy/colour/
   thorax/shadow done; silhouette/stills/report were still open and were completed.
2. **Silhouette state verified from a fresh render.** `s4verify-sil.png` was rendered from
   the current `build_moth.py` (exit 0) and measured with `tools/shape.py`: A1 1.434,
   A2 fan 0.777 / render 0.780, A3 7 designed / 5 clear, A4 0.345, photo own-mirror 0.901.
   These agree with the kept try C. The last value the first session had printed was try D
   (A1 1.457, fan 0.777, render 0.776); try D is the half-amplitude variant from the
   iteration log, reconstructible in principle, but it loses on A2 (0.776 vs 0.780) while
   A1 stays in band either way, so the brief's stop rule keeps try C. Recorded plainly:
   the current code is **not** byte-identical to try D, and try D was not re-applied.
   No new silhouette exploration was started.
3. **Stills re-rendered from the final code.** `build_moth.py` was edited again after the
   stills were first written (speck fix / shadow back-off), so all four were re-rendered:
   stageB-f0 689453 B, stageB-obl 223076 B, stageB-ns 693601 B, stageB-sil 156875 B, all
   exit 0, all newer than the shipped `build_moth.py`. Re-renders reproduce every metric
   within TAA noise (Eevee is not bit-deterministic: mean |diff| 0.001, 49 px > 2 levels).
4. **Numbers re-measured with an independent tool.** `tools/e3geo.py` (a separate
   implementation using the projected wing outline as a fixed distance reference, applied
   to both states) gives E3 0-3 -14.22 -> -12.86 %, 3-6 -22.29 -> -23.48, 6-12 -22.38 ->
   -24.23, 12-24 -12.84 -> -14.59, 24-60 -1.02 -> -1.19, 0-3 px darkest 5 % -54.96 ->
   -59.26, and D1 on the same fixed mask 3.845 -> 4.752 (+23.6 %). `tools/e3seg.py`
   (segmented reference fixed to `base-ns.png`) gives -4.28 -> -3.43, -18.69 -> -16.49,
   -26.60 -> -27.49, -15.41 -> -17.53, -1.43 -> -1.74, p5 -50.94 -> -55.75. Colour median
   111,156,148 (n=32177), thorax 59x53 / aspect 1.113 / PCA 1.358 / 2276 px, A1 1.434,
   A2 0.777 / 0.780, A3 7/5, A4 0.345 - all agree with the tables above.
5. **C1 re-checked on the final still.** No abdomen wedge. The dark central shape that can
   read as an abdomen triangle at 1x is the **thorax**: the dark patch (bbox 564..617 x
   206..256, 2034 px) overlaps the thorax blob footprint (564..622 x 203..255, 2276 px) by
   **99.7 %**, and nothing inside the moth silhouette falls below 0.65x the wing median
   (min in the centre box 93.7, i.e. brighter than proto3's 89.7 in absolute terms). The
   "dark triangle" impression is the body-vs-wing value gap opened by brightening the wings
   (issue 2 above), not a seam-coverage failure. C1 stays PASS.
6. **Nothing invented.** A4 (brief 0.494) and the E3 absolute normalisation could not be
   reproduced from this build and remain reported with explicit definitions only.
