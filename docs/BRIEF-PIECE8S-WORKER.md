# BRIEF — piece8s base render (worker session, PC)

Working folder (PC): `D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\` (create it this session).
Base: `..\blender-moth-proto5\` (build_moth.py v5, md5 `8858284bd1d3df4ffb6f5289aacc3fe0`) — NEVER modify the proto5 folder.
Spec source: `..\BRIEF-PIECE8S.md` (VM-originated master brief; summary of its hard requirements is folded into this brief — where they differ, BRIEF-PIECE8S.md wins).
Deliverables stay in the piece folder on the PC; the coordinator fetches them to the VM. Do NOT copy anything to the VM yourself.
Model: deepseek-v4.1-flash only (already fixed by the caller).

## What this is

The 8-second main piece's BASE render: the proto5 moth pair, same scene, same asset,
but a 192-frame (8.000 s @ 24 fps) choreography with a loop seam at f0/f191.
No caustics/water composite — that is a later stage. Eevee 64 samples, motion blur,
832×480, camera fully locked, background plate unchanged — all as proto5.

## Step 1 — working copy

Copy `blender-moth-proto5` to `blender-piece8s`, excluding: `renders\` contents,
`*.mp4`, `checkpoint.json`, `encode-report.json`, `render-timing.json`, `__pycache__`.
Keep: `build_moth.py`, `run_blender.py`, `assemble.py`, `tools\`, `background-bark-832x480.png`,
`f00.png` (proto5's delivered f0 — needed for the parity gate), `_dbg\` may be re-created fresh.
Use jobkit (`D:\mutaro-work\jobkit\checkpoint.py`, see its README) with a Checkpoint at
`blender-piece8s\checkpoint.json`; record units as you complete them:
`copy, build_edit, f0_parity, render_start, render_done, assemble, seam_verify, report`.
Log the session in `PIECE-LOG.md` in the piece folder (append as you go, incl. incidents).

## Step 2 — build_piece.py

Create `build_piece.py` as a copy of `build_moth.py` with ONLY these changes:
- `F_END = 191` (frames 0..191), `BLEND_OUT = piece.blend`.
- The ANIMATION section (MOTION_TABLE + the per-frame functions + shadow keys) replaced
  by the piece choreography in Step 3 below.
- `params` dict: `proto: "piece8s"`, a new `choreography` string summarising the piece
  timeline, `frames: [0, 191]`.
- The motion log (`_dbg/motion-log.json`) extended with the left shadow's per-frame
  opacity + location + scale so the seam check can include it.
Geometry, materials, hair, scale layer, camera, sun, background, moth construction:
byte-identical construction to proto5. Do not "improve" anything else.

### Parity gate (before the full render)
Render frame 0 from build_piece.py (`--still 0`) and diff against proto5's `f00.png`
in the folder: PASS = mean |Δ| ≤ 0.5 (0–255 scale) and no structural difference
(report % of pixels with |Δ| > 8; expect TAA-level, proto5 measured 0.107 mean on a
double build). If it fails, the copy/edit broke the asset — fix before continuing.

## Step 3 — choreography (all values are exact; functions of frame number only)

Constants: rest wing angle 3.0°, raise target 28.0°, stroke top 32.0°, twist 7.0°,
BASE_HZ 8.0 (period 3.0 f), beat factors [1.00, 0.90, 1.12, 0.94] cycling in order
(period = BASE_T × factor[i mod 4]), seeded literals as proto5. Liftoff LIFT = 72.
Beat grid starts at 72 and runs to ≥196. The lag wing (proto5's `left_up_lag` logic,
rebased to LIFT=72) differs only during the first beat after liftoff; after that both
wings share `stroke_phase`. Twist ramps in with clamp01((f−71)/3) after liftoff and
ramps out ×(1 − sstep((f−160)/8)), so twist = 0 from f168 on.

### LEFT moth (flyer) — root rest = REST_LOC_L, base yaw 50°
- f0–f61: at rest. Wing 3.0°, root at rest, rotation offsets 0.
  (proto5's clip-start f3 twitch and f6 settle dip are NOT carried over — f0 is a loop
  seam, the left moth starts from clean rest. Note this in the RESULT.)
- One antenna event: `_piecewise([[37,0],[40,2.0],[42,0]], f)` deg, both antennae, 0 elsewhere.
- Takeoff (MEASURED structure, re-timed by back-calculation from liftoff f72;
  proto's span was raise-start→liftoff = 10 f with hold in the last 2 f):
  - Wing raise f62–f69: 3.0 + 25.0 × sstep((f−62)/7.0) → 28.0 at f69.
  - Hold f70–f71: 28.0.
  - Liftoff f72 (first downstroke simultaneous; 1-frame lead on the leading wing as proto5).
- Flight wing (f72–f167), lead wing: flight(f) = 32.0 × cos(2π·stroke_phase(f)).
  Braking blend from f156: brake(f) = 18.0 + 14.0 × cos(2π·stroke_phase(f));
  b(f) = sstep((f−156)/8) clamped; wing = (1−b)·flight + b·brake for 72 ≤ f < 168.
- Settle (f ≥ 168): amp(f) = 14 × (1 − sstep((f−168)/10)) → 0 at f178;
  mean(f) = 18 − 15 × sstep((f−168)/12) → 3.0 at f180;
  wing(f) = mean(f) + amp(f)·cos(2π·stroke_phase(f)). From f180: exactly 3.0.
- Root path — `_pchip` monotone cubics through these knots (use the script's REST_LOC_L
  values, not literals, where a knot says "rest"):
  - X: (0,rest)(71,rest)(72,rest)(76,−0.95)(80,−1.75)(84,−2.60)(92,−3.40)(132,−3.40)
       (138,−2.90)(142,−2.25)(148,−1.80)(156,−1.05)(162,−0.72)(166,−0.63)(168,rest), flat_end=True
  - Y: (0,rest)(71,rest)(72,rest)(76,0.30)(80,0.62)(84,0.95)(92,1.25)(132,1.25)(140,1.15)
       (148,0.95)(156,0.62)(162,0.28)(166,0.05)(168,rest), flat_end=True
  - Z: (0,0)(71,0)(72,0.015)(76,0.55)(80,0.85)(84,1.05)(92,1.15)(132,1.15)(148,1.02)
       (156,0.85)(162,0.55)(166,0.18)(168,0.015)(170,0.0), plus a settle dip
       dip(f) = _piecewise([[169,0],[171,−0.012],[175,0]], f) added on top (0 after f175).
  - PITCH deg: (0,0)(70,0)(71,0.5)(72,2.5)(76,18)(80,21)(84,20)(92,18)(148,18)(156,20)
       (162,26)(166,14)(168,2.5)(172,0.5)(176,0) — clamped 0 afterwards.
  - ROLL deg: (0,0)(71,0)(72,0)(76,−22)(80,−30)(84,−33)(92,−34)(148,−30)(156,−20)
       (162,−8)(166,−2)(168,0) — clamped 0 afterwards.
  - YAW deg: (0,0)(71,0)(72,0)(76,11)(80,17)(84,20)(92,22)(148,20)(156,14)(162,6)
       (166,1)(168,0) — clamped 0 afterwards.
  Touchdown is f168: root location == REST_LOC_L exactly at f168 (path knots) and for
  every f ≥ 170 (Z dip returns to 0 at f175; X/Y clamped at rest from f168).
- Pixel gates (measure with the same moth-pixel method the proto rounds used for
  "fully out = 0 px" — see tools/exitcap.py / motion_verify.py; state the method):
  px(f80) > 0; px(f84) == 0; px(f120) == 0; first frame with px > 0 during the return
  ∈ [139,143]; root at f168 == rest exactly (animation data). You may steepen/soften
  path knots slightly to hit these gates, but keep every curve monotone and log any change.

### RIGHT moth (never translates — root Δ = 0 for all 192 frames)
Reuse proto5's functions/shapes verbatim; only the tables extend:
- breath_keys: [[0,3.0],[6,6.5],[14,3.0],[23,3.0],[28,6.5],[37,3.0],[47,3.0],[56,6.5],
  [66,3.0],[75,3.0],[84,6.5],[94,3.0],[103,3.0],[112,6.5],[122,3.0],[131,3.0],[140,6.5],
  [150,3.0],[160,3.0],[170,6.5],[180,3.0],[191,3.0]]
- antenna_both: [[4,0],[5,2.5],[8,0]] (proto5, unchanged)
- wing flick: [[12,0],[15,2.5],[17,0]] (proto5 code, unchanged)
- antenna sweeps (both antennae, additive): [[28,0],[29,3],[30,0]] and [[104,0],[105,3],[106,0]]
- partial raises (additive): [[31,0],[37,9],[39,0]] and [[118,0],[124,9],[126,0]]
- No reaction to the other's departure/return (no data — deliberate, as proto5).
- Right shadow decal: untouched (static, as proto5).

### LEFT shadow decal
- Opacity keys: (0,0.95)(71,0.95)(75,0.62)(80,0.30)(86,0.10)(94,0.03)(100,0.0)(150,0.0)
  (160,0.25)(166,0.60)(168,0.95)(191,0.95).
- Location/scale offset keys (dx, dy, scale-mult) relative to the proto5 base
  (base_loc = REST_LOC_L + SH_OFF, base scale SCALE_L), same construction as proto5:
  {0:(0,0,1.0), 71:(0,0,1.0), 76:(−0.09,0.05,1.03), 82:(−0.20,0.11,1.08),
   90:(−0.28,0.15,1.12), 100:(−0.34,0.18,1.15), 140:(−0.20,0.11,1.08),
   155:(−0.09,0.05,1.03), 164:(0,0,1.0), 168:(0,0,1.0), 191:(0,0,1.0)}.
  At f191 the decal must be exactly at its f0 state (it is, by the keys above).

## Step 4 — seam verification (tools/seam_verify.py, new; run BEFORE assembly)
1. Animation data: iterate every action/fcurve in the built scene; compare each
   fcurve's value at frame 0 vs frame 191 (keyframe values where keyed per-frame;
   fcurve.evaluate() for sparse ones like the shadow opacity node). PASS = max abs
   difference over ALL fcurves == 0.0 exactly. Also assert the motion-log rows
   f0 == f191 field-by-field. Write `seam-verify.json` with the max diff, the fcurve
   count, and any offender. (Hair particle systems carry no fcurves — excluded by
   nature; state this in the RESULT.)
2. Right root: assert location identical at every frame (Δ = 0 for all 192 rows of
   the motion log).
3. Left root: assert loc(f168..f191) == REST_LOC_L for f ≥ 175 and loc(f191) == loc(f0)
   exactly; assert loc(f0..f71) == REST_LOC_L except the Z dip is 0 there (root fully
   planted until liftoff: X/Y at rest f0–f72, Z = 0 until f71, 0.015 at f72).

## Step 5 — render (DETACHED) + assemble
- Launch the full render DETACHED (WMI Win32_Process.Create) so it survives this
  session: `python run_blender.py build_piece.py --render`, stdout/stderr to
  `render.log` in the piece folder. Poll `renders\` frame count + render.log every
  ~60 s. Record render wall time. If this session dies, the render continues; the
  coordinator will dispatch a completion session — so keep checkpoint.json and
  PIECE-LOG.md current at every stage boundary.
- Assemble with an adapted `assemble_piece.py` (from assemble.py): expect 192 frames,
  ffmpeg libx264 crf 16 yuv420p 24 fps → `piece-base.mp4`, ffprobe-verify BEFORE the
  atomic commit (832×480, 24 fps, nb_read_frames == 192, duration 8.0 s, codec h264),
  write `encode-report.json`, checkpoint `assemble`.
- Render-diff seam check (after assembly): PIL/numpy diff of renders/frame_0000.png
  vs renders/frame_0191.png — overall mean |Δ| (report on BOTH 0–255 and 0–1 scales);
  moth mask = |frame_0000 − background-bark-832x480.png| > 12 on any channel, dilated
  3 px; report mean |Δ| inside mask / outside mask and % pixels with |Δ| > 8 in/out.
  GATE: outside-mask mean ≤ 0.002 (0–1 scale). Inside: report only (expect TAA/hair
  shimmer magnitude). Append results to seam-verify.json.
- Stills: copy renders/frame_0000.png→f0.png, 0072→f72.png, 0084→f84.png,
  0120→f120.png, 0168→f168.png, 0191→f191.png (piece folder root).

## Step 6 — RESULT-PIECE.md (piece folder)
Contents:
- Timeline as implemented, with MEASURED actuals: liftoff frame, 0-px interval
  (first/last), return entry frame, touchdown frame, and the pixel-gate table.
- Seam verification: fcurve max diff, motion-log f0/f191 equality, render-diff stats
  (inside/outside mask, both scales), gate verdicts.
- Provenance label per segment (use the dossiers' vocabulary):
  rest life + right-moth proto events = proto5 labels unchanged (TRANSFERRED /
  MEASURED shape / INVENTION as in RESULT-PROTO5); right-moth added events (sweep
  f104–106, raise f118–126) = same shape/order labels, placement INVENTION;
  breathing extension = INVENTION; left antenna event = INVENTION; takeoff structure
  f62–72 = MEASURED (re-timed by back-calculation); exit path = proto3 measured
  envelope in shape/direction, timing redesigned to the user's ≈0.5 s rule (DESIGN);
  return path / braking / touchdown / settle dip = INVENTION (design, no footage);
  wingbeat 8.0 Hz + factors = TRANSFERRED; first-stroke 1 f lead = INVENTION;
  left clip-start micro-events dropped = editorial (seam), stated.
- ffprobe block, render wall time, mp4 bytes.
- Honest remains: hair/TAA shimmer across the seam (~0.9 % of pixels, params exact);
  landing is invented; anything that failed a gate or was tuned, with numbers.

Final message: summary with the key numbers (ffprobe line, seam diffs, pixel gates,
wall time). Do NOT copy deliverables to the VM.
