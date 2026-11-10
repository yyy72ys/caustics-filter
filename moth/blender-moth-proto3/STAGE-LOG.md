# proto3 STAGE 1 log — shape only

Source copied: blender-moth-proto2 -> blender-moth-proto3 (full folder, work only in proto3; proto2 untouched, verified Nt,Nu still 44,18 there).
Preserved (not listed, untouched): background plate, choreography (left exits by f24, right 3->6.5 deg breathing), Eevee 64 samples, motion blur, 832x480, 48 frames @24fps, fixed production camera, sun/ambient rig. BLEND_OUT left as proto2.blend (not listed); proto number in build_params.json left as 2.

## Changes
- Change 1 (new wing outline FAN_LE/FAN_TE + grid Nt 44->64, Nu 18->22): applied. Replaced both constant lists with the 11-row LE3 / 14-row TE3 values from the brief; make_fan_wing now Nt,Nu = 64,22.
- Change 2 (scalloped hindwing helper fan_outline, shared mesh+shadow): applied. Added fan_outline(n) sampling via sample_curve/smooth_profile with TE normal displacement for t in (0.17,0.60), w = sin(pi*(t-0.17)/0.43)**0.7, d = 0.0070*sin(2*pi*7*(t-0.17)/0.43+0.5)*w, normal oriented toward LE. make_fan_wing and make_shadow_image both call it. Old inline outer-margin scallop (sstep 0.42-0.80, 0.0075*sin(2*pi*6t+0.7)) removed so mesh and shadow agree.
- Change 3 (wing root z 0.026 -> 0.015): applied (wing object location line only; sole occurrence).
- Change 4 (fan surface z formula): applied. Note: code did not match brief verbatim — actual old line was 0.010*sin(pi*u)*(0.35+0.65t) + 0.0022*sin(2*pi*2.5t+3u)*t - 0.012t^2, not the brief's 0.008*sin... - 0.012t^2. Replaced whole computation with the brief's 5-term formula (camber 0.008 + cover dome 0.024*exp(-(xr/0.045)^2)*(1-sstep(0.02,0.16,t)) + bark relief 0.0035*sin(2*pi*(3.1lx+2.3ly)+1.7)*(0.3+0.7u) - margin droop 0.010*u^6*sstep(0.25,0.75,t)*(0.55+0.45*sin(2*pi*2.6t+0.8)) - tip 0.006t^2), lx/ly pre-mirror, xr = side*(lx+0.030), _sstep as smoothstep.
- Change 5 (abdomen sunk/flattened): applied. Radius profile x0.85 overall; zb 0.042-0.024t^2 -> 0.013-0.006t^2; cross-section x factor 0.92->0.80, z factor 1.0->0.40.
- Change 6 (thorax/head): applied. Thorax radii (0.058,0.082,0.055)->(0.062,0.080,0.050), centre z 0.040->0.036; head centre z 0.046->0.042.
- Change 7 (legs tucked): applied. Hip z 0.038->0.020, mid z 0.018->0.010, foot 0.001 unchanged.
- Change 8 (antennae filiform): applied. Deleted entire pecten barb loop (6 barbs/side). Main tube radii x0.8. Notes: brief cited base 0.0016->0.0013 / tip 0.0006->0.0005 but code had base 0.0036 / tip ~0.00162; applied x0.8 to actual values (base 0.00288 / tip ~0.00130). Tip y 0.244 matched code (0.124+0.150-0.030); changed linear coeff 0.150->0.164 so tip = 0.124+0.164-0.030 = 0.258.
- Change 9 (material t remap): applied. Bands 0.50->0.58, 0.72->0.765, 0.875->0.895; crease ul 0.42+0.58t -> 0.309+0.690t (values + comment); vein fade sstep(0.55,1.0)->sstep(0.62,1.0); margin sstep(0.84,1.0)->sstep(0.87,1.0). Nothing else touched.
- Change 10 (debug flags --silhouette / --oblique DEG / --noshadow): applied. Single-moth setup hides Moth_L + Shadow_L, zeroes Moth_R rotation/scale/location, clears wing animation to flat 0. Silhouette: white emission on all right-moth meshes, bark + both shadows hidden, world black (color 0,0,0 strength 1.0), new ortho top camera (identity rotation, ortho_scale 3.4 at moth x,y,z=5). Oblique: normal materials, world grey (0.35,0.35,0.35 strength 1.0), persp 50mm cam at elevation DEG / azimuth -30 deg / dist 3.2 with track-quat look_at moth centre (x,y,0.05). Noshadow: both decals hide_render. Extra: brief runs use positional `still <file.png>` (not --still/--out); added positional-still support (STILL=0, OUT=<file> under HERE) so the literal commands work.
- Change 11 (assemble.py output): applied. OUT proto.mp4 -> proto3.mp4 (note: code had proto.mp4, not proto2.mp4; docstring updated too). Checkpoint verify->encode->ffprobe->atomic flow unchanged.

## Renders (in blender-moth-proto3, Blender 4.2.3 LTS)
1. python run_blender.py build_moth.py -- still stage1-f0.png -> exit 0, stage1-f0.png 690921 bytes
2. python run_blender.py build_moth.py -- still stage1-sil.png --silhouette -> exit 0, stage1-sil.png 159288 bytes
3. python run_blender.py build_moth.py -- still stage1-obl.png --oblique 35 -> exit 0, stage1-obl.png 575427 bytes

## STAGE 1b — trailing-edge correction only
- Change applied: FAN_TE (file's name for the brief's FAN_TE3) first six rows replaced with the continuous convex rear arc values: (0.00,-0.030,-0.312), (0.055,-0.029,-0.305), (0.11,-0.024,-0.298), (0.16,0.000,-0.295), (0.244,0.030,-0.300), (0.345,0.078,-0.318). All other rows and everything else untouched.
- Clarification fix: --oblique now hides the bark plate and the shadow decal (plate.hide_render/viewport + sh_r.hide_render/viewport, same as --silhouette; sh_l was already hidden by the single-moth block), keeping normal wing/body materials and the grey world.

## Renders 1b (in blender-moth-proto3, Blender 4.2.3 LTS)
1. python run_blender.py build_moth.py -- still stage1-f0.png -> exit 0, stage1-f0.png 688572 bytes
2. python run_blender.py build_moth.py -- still stage1-sil.png --silhouette -> exit 0, stage1-sil.png 156994 bytes
3. python run_blender.py build_moth.py -- still stage1-obl.png --oblique 35 -> exit 0, stage1-obl.png 219045 bytes

## STAGE 2 — materials / texture (brief proto3-stage2-brief.md)
Conventions note: brief header says fan UV.x = u / UV.y = t and cites
make_wing_material / make_body_material, but the code uses mesh UV = (t, u)
(UV.x = t span, UV.y = u chord), material t_ = SeparateXYZ X / u_ = Y, and
functions mat_wing() / mat_body(). Code mapping preserved everywhere below;
fringe UV kept as (t, u). No silent guessing beyond this.

- Change 2-0 (cover dome 0.024->0.016, width 0.045->0.055): applied.
  make_fan_wing z dome term only.
- Change 2-1 (lunulate bands): applied to all three bands. Shared
  L(u) = 0.25+0.75*|sin(pi*u*6.5)|^0.8 and W(u) = 0.85+0.47*exp(-((u-0.48)/0.30)^2)
  built once from u_; inside build_band the distance |t-cen| is divided by W
  before the MapRange smoothstep pulse, and the band strength term is
  multiplied by L (amt4). Band colour stays #edf2e9; no black introduced.
- Change 2-2 (discal marks): applied. Two gaussian mixes toward #4c6e52
  inserted after the vein mix and before the bands: m1 = 0.16*exp(-(((t-0.50)/0.05)^2
  + ((u-0.38)/0.10)^2)), m2 = 0.13*exp(-(((t-0.34)/0.045)^2 + ((u-0.78)/0.09)^2))
  via SUBTRACT/DIVIDE/POWER(2)/ADD/MULTIPLY(-1)/EXPONENT/MULTIPLY(amp).
- Change 2-3 (scale grain g): applied. New Noise on wing UV through Mapping
  scale (110, 200, 1) exactly as listed (Noise Scale left at default 5.0 —
  brief gives none; Detail 2, Roughness 0.7), Fac = g. Colour: col multiplied
  by (0.94+0.11*g) via MULTIPLY_ADD + CombineColor + Mix MULTIPLY just before
  BSDF Base Color. Roughness: (g-0.5)*0.10 added to MapRange result (see 2-5).
  Bump: g feeds Bump (strength 0.15, distance 0.0015), replacing the old n2
  Height. Note: no separate vein-bump node existed (single n2 bump only), so
  there was nothing to add to — old broad-noise contribution dropped, n2 kept
  for its colour-ramp and roughness-base roles.
- Change 2-4 (fade mottle + gradient): applied after rim, before grain
  multiply. Fade noise on UV through Mapping (3, 2, 1) exactly (Noise
  Scale/Detail/Roughness not listed — left at defaults 5.0/2/0.5), p = Fac,
  smoothstep(0.55, 0.80, p)*0.16, hind factor 1.0+0.5*hindF0 reusing the
  existing crease MapRange (dsub 0..0.07 -> 0..1: 1.0 forewing, 1.5 hindwing),
  mix toward #a8c4b2. Radial: col multiplied by (0.90+0.10*smoothstep(0.16, 0.90, t)).
- Change 2-5 (matte): applied. Specular IOR Level 0.12->0.05; Sheen Weight
  0.15->0.10; roughness MapRange To Min/Max 0.66/0.85 -> 0.74/0.93; rim mix
  factor 0.08->0.035. Note: grain offset (g-0.5)*0.10 is added after the
  MapRange, so the final roughness spans ~0.69-0.98 around the 0.74-0.93 base.
- Change 2-6 (fringe): applied and KEPT. Grid extended with u = 1.012, 1.024,
  1.035 (NCOL = 23+3 = 26); x/y use the same linear LE->TE extrapolation, which
  equals TE + normalised(TE-LE)*(u-1)*chord_length; z lowered 0.0012 per step.
  UV kept (t, u). Alpha = 1-smoothstep(1.005, 1.035, u)*(0.5+0.5*g) to BSDF
  Alpha; blend_method BLEND + use_transparent_shadow True. make_shadow_image
  untouched (fan_outline base outline only), so the shadow texture excludes the
  fringe with no code change. No Blender error (both renders exit 0); stage2-f0
  shows no ribbon artifact; stage2-obl shows a serrated grain-gated trailing
  edge judged to be the intended cutout, not a detached ribbon — kept.
- Change 2-7 (thorax fur): applied to BodyMat only (thorax+head). mat_body()
  extended with sheen / bump_dist params (defaults = old values) so Leg/Abdomen/Eye
  are untouched; BodyMat now fur bump 0.55->0.95, distance 0.006->0.009,
  Sheen 0.4->0.75. New mat_fur_shell("FurShellMat", "#6f8f66"): base = body
  colour x1.1, Alpha = (Object-noise scale 40 detail 2 Fac > 0.45), BLEND,
  visible_shadow False where available. build_moth duplicates ONLY thorax and
  head meshes, scaled 1.13 about each mesh's own vertex centroid, parented to
  the moth root (same static parenting as the body). Shells appended to
  body_parts so the existing single-moth hide + silhouette white-override apply.

## Renders 2 (in blender-moth-proto3, Blender 4.2.3 LTS)
1. python run_blender.py build_moth.py -- still stage2-f0.png -> exit 0, stage2-f0.png 687653 bytes
2. python run_blender.py build_moth.py -- still stage2-obl.png --oblique 35 -> exit 0, stage2-obl.png 220563 bytes

## STAGE 3 — light / shadow + full render (brief proto3-stage3-brief.md)
- Change 3-0 (flatten centreline tent): applied. Abdomen zb 0.013-0.006t^2 -> 0.006-0.006t^2 (top now ~=0.018); cross-section z factor 0.40 -> 0.35. make_fan_wing dome amplitude 0.016 -> 0.006 (width 0.055 kept; sheet over centreline ~=0.024 > 0.018).
- Change 3-0b (grain scale): applied. Scale-grain Noise (_gn) Scale input 5.0 (default) -> 1.0; Mapping (110,200,1), Detail 2, Roughness 0.7 unchanged.
- Change 3-1 (shadow in ROOT coords): applied. make_shadow_image rasterises in moth-root coords: right fan poly points += (+0.030,+0.012), left (mirrored) += (-0.030,+0.012); body ellipse untouched (already root coords). Canvas SH_X0..X1 = -0.62..0.62, SH_Y0..Y1 = -0.42..0.30, texture 512x296. build_shadow_decal (SH_* rect) now spans exactly that root rect: width 1.24, height 0.72, centred (0,-0.06), so texture and plane align. Deviation: brief says "the decal stays a child of the moth root as before", but the code parents nothing — both decals live in the scene collection at REST_LOC+SH_OFF with the moth's rotation/scale (plus the left shadow's fade/drift keyframes). Kept that world-placed scheme (only the rect changed) to preserve the sun offset and existing shadow animation; no reparenting.
- Change 3-2 (margin stroke): applied. New _stroke_mask helper rasterises each fan's closed outline polyline (root coords, scalloped TE from fan_outline) as a 4.0px-wide line (pixel-centre distance <= 2.0), blurred r=2 via _box_blur. Final alpha = min(0.88, 0.62*sharp + 0.45*wide + 0.85*stroke). Body ellipse excluded (fans only, per brief).
- Change 3-3 (AO): applied. Eevee GTAO factor 0.6 -> 0.9, distance 0.10 kept; sun/ambient unchanged. Dropped: nothing.
- Render-command note: brief lists the full animation as `python run_blender.py build_moth.py -- render renders`, but build_moth.py only triggers on the `--render` flag (bare `render` is ignored and renders nothing). Ran `python run_blender.py build_moth.py -- --render` instead, which writes the same default `renders/` output the brief names.

## Renders 3 (in blender-moth-proto3, Blender 4.2.3 LTS)
1. python run_blender.py build_moth.py -- still stage3-f0.png -> exit 0, stage3-f0.png 685806 bytes
2. python run_blender.py build_moth.py -- still stage3-ns.png --noshadow -> exit 0, stage3-ns.png 690168 bytes
3. python run_blender.py build_moth.py -- --render -> exit 0, 48 frames into renders/ (frame_0000..0047.png)
4. python assemble.py -> exit 0, frames ok=48 bad=0; proto3.mp4 257847 bytes (251.8 KB); ffprobe: width=832 height=480 fps=24/1 frames(nb_read_frames)=48 duration=2.000000 codec=h264. Raw line: {"codec_name": "h264", "width": 832, "height": 480, "r_frame_rate": "24/1", "duration": "2.000000", "nb_read_frames": "48"}

## STAGE 3b — centreline seam patch (brief proto3-stage3b-brief.md)
- Change 3b-1 (half-patch appended to each fan mesh in make_fan_wing): applied. After the fan grid verts/faces, 11 rows x 3 pts (A=P(s), B=midpoint, C=centreline at same y; s=k/10, P(s)=LE0+s*(TE0-LE0) with LE0=(0.002,0.052), TE0=(-0.030,-0.312)) appended to the SAME mesh/UV/material; z = fan height formula at t=0.02, u=s with the point's own lx,ly + 0.0015 lip; UV=(0.03, 0.05+0.9*s); 20 quads (A-B, B-C strips) with the fan's side-based winding (side>0: [A,B,Bn,An]/[B,C,Cn,Bn], else reversed). Flap/perch animation, silhouette override, shadow behaviour unchanged; shadow texture untouched (covered-area change negligible, as brief allows).
- Deviation: brief lists centreline local x_c = -side*0.030 (wing root offset side*0.030). The object offset is actually sign*0.030 with sign=mir*side, so the true mesh-local centreline is -sign*0.030 = -mir*side*0.030. Implemented in grid coords (pre-mirror) as x_c=-0.030 then mesh x=sign*lx, which equals the brief for mir=1 (Moth_R: both stills, right-moth closeups) and additionally puts Moth_L's patches on the correct side (the literal formula would mirror Moth_L's patches away from the centreline). z uses xr=side*(lx+0.030) exactly as the fan formula (sign irrelevant: squared in the dome term). Nothing else touched.

## Renders 3b (in blender-moth-proto3, Blender 4.2.3 LTS)
1. python run_blender.py build_moth.py -- still stage3-f0.png -> exit 0, stage3-f0.png 685091 bytes
2. python run_blender.py build_moth.py -- still stage3-obl.png --oblique 35 -> exit 0, stage3-obl.png 219894 bytes (new file; stage 3 had stage3-ns.png instead of an oblique)
3. python run_blender.py build_moth.py -- --render -> exit 0, 48 frames into renders/ (frame_0000..0047.png)
4. python assemble.py -> exit 0, frames ok=48 bad=0; proto3.mp4 259253 bytes (253.2 KB); ffprobe: width=832 height=480 fps=24/1 frames(nb_read_frames)=48 duration=2.000000 codec=h264. Raw line: {"codec_name": "h264", "width": 832, "height": 480, "r_frame_rate": "24/1", "duration": "2.000000", "nb_read_frames": "48"}

## STAGE 3b — follow-up run (no code changes)
- Code state: the centreline seam patch (3b-1 above) was applied in the PREVIOUS session. Nothing in build_moth.py / assemble.py was edited in this run.
- f0: stage3-f0.png was re-rendered in the previous session (685,091 bytes) and was not re-rendered here.
- Full render: python run_blender.py build_moth.py -- --render -> exit 0, 48 frames (renders/frame_0000..0047.png), Blender 4.2.3 LTS, 145.7s wall.
- Assemble: python assemble.py -> exit 0, frames ok=48 bad=0; encode 0.32s.
- proto3.mp4: 259366 bytes (253.3 KB). Prior 3b encode was 259253 bytes; same code, re-encode only (+113 bytes, not a content change signal).
- ffprobe (this run): width=832 height=480 fps=24/1 (avg 24/1) frames(nb_read_frames)=48 duration=2.000000 codec=h264. Raw line: {"codec_name": "h264", "width": 832, "height": 480, "r_frame_rate": "24/1", "avg_frame_rate": "24/1", "duration": "2.000000", "nb_read_frames": "48", "size": 259366}
