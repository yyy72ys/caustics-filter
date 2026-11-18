# BRIEF — 8 s main piece, base render v2 (user review rework)

PC work dir: D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\ (v1 build:
build_piece.py; moth = proto5 = moth v1, DO NOT remodel the moth). Use jobkit
checkpoints per stage. If the render may exceed ~10 min, launch it DETACHED via
WMI (Win32_Process.Create) and verify detachment from a NEW session (process
alive AND output advancing, twice). NEVER taskkill python.exe / blender.exe
belonging to other jobs. Foreground otherwise; do not end until deliverables
are written and verified.

## The user's review of v1 (fix ALL five points)
1. "The shadow is weird." — diagnose first (below), then fix.
2. "The opening stillness is long and I can't tell what's happening." — v1
   rested 2.5 s. Cut the opening rest to ~1 s and put the antenna event early.
3. "The wings should be fluttering more." — in v1 the flight was visible for
   only ~8 frames. In v2 the wingbeat must be LARGE and clearly visible
   throughout the flight (wide sweep angle, body bob coupled to the beats).
4. "It should fly around inside the frame." — NO exit from the frame at all.
5. "Going out and coming back feels lame." — replace the out-and-back with a
   circuit INSIDE the frame (below).

## Step 0 — shadow diagnosis (log before fixing)
Inspect v1 frames f0068..f0084 (takeoff/flight) and the perched frames: state
in PIECE-V2-LOG.md what the shadow does wrong (candidate defects: hard dark
blob glued under the airborne moth; shadow not offsetting along the sun
direction with height; shadow not softening/lightening with height; contact
shadow offset from the body). Then implement the fix:
- Perched: contact shadow as in v1 (it passed proto gates).
- Airborne: shadow offset along the sun azimuth proportional to height,
  blur radius growing with height, opacity falling with height (monotone in
  all three; record the three curves' sample values in the log).

## New choreography (192 frames, 24 fps, seam f191 == f0 PIXEL-IDENTICAL,
right/staying moth untouched, root delta 0 over all frames)
- f0..~f24: rest (~1 s). Antenna event early (~f12..f18, the proto5 sweep).
- ~f30..f40: wing raise (measured structure, proto4/5) then liftoff.
- ~f40..f150: FLIGHT INSIDE THE FRAME. A fluttery circuit: rise and cross to
  the left across the upper frame, descend along the left side, sweep back
  rightwards low (passing below/in front of the staying moth's area WITHOUT
  overlapping it), one small loop or hover-pause with heavy flutter, then the
  approach. The moth stays FULLY inside the frame at all times (margin >=
  ~10 px; verify per frame from the alpha bbox). Path smooth at the large
  scale, with deterministic SEEDED irregularity at the small scale (bobbing
  loosely phase-locked to the wingbeats, small heading wobbles) — the proto4/5
  method: measured structure + seeded irregularity, no random drift.
- Wingbeat in flight: 8 Hz (proto5 transferred value), LARGE amplitude (wings
  sweeping a wide angle — visibly flapping in every flight frame pair), body
  pitch/bob per beat. Braking flare before touchdown as in v1 (invention,
  labelled).
- ~f160..f168: touchdown at the EXACT original perch. Settle to the f0 pose by
  ~f180. Final hold only as long as the pixel-exact seam needs.
- Takeoff/landing structure: measured (proto4/5 dossiers) as before; the
  circuit path itself is choreographic invention — label it INVENTION in the
  RESULT (as v1 did for landing).

## Gates (verify and record in RESULT-PIECE-V2.md)
- 192 frames, 8.000 s, 832x480 (ffprobe).
- Seam: f191 vs f0 pixel-identical (0 differing pixels); all fcurves periodic.
- Right moth root delta = 0 over all 192 frames.
- Moth fully in frame on every flight frame (alpha-bbox margin log, min
  margin reported).
- Shadow: offset/blur/opacity monotone in height (sampled values in the log).
- Wingbeat: beat onsets during flight at ~8 Hz with the large amplitude
  (report the wing-angle range in flight vs v1's).

## Deliverables -> VM delivery\piece-8s-base-v2\ (keep PC copies)
piece-base-v2.mp4; stills f0, f60, f96, f144, f168, f191; build_piece_v2.py;
PIECE-V2-LOG.md (incl. shadow diagnosis); RESULT-PIECE-V2.md; seam-verify.json.
Render idempotent per frame (skip existing), as v1.
