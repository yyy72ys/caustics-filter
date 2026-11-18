# BRIEF proto4-C: choreography v2 — measured takeoff sequence + seeded irregular micro-motion

You are running headless on the user's PC. Execute this brief exactly, to completion, without stopping to ask questions. If a step fails, log it and continue.

## Context
Working folder: `D:\yohey2026\dsh\muse\caustics-filter\blender-moth-proto4\`. Its build_moth.py is proto3's model with stage-B detail polish (colour/thorax/shadow/silhouette) already applied — do NOT change any geometry, materials, shadow, camera, lights, or the background plate in this stage. You are changing ONLY the animation (choreography) of the two moths, then rendering the final video.

The full motion research is in `D:\mutaro-work\rally-in\MOTION-DOSSIER.md` — READ IT FIRST. The measured facts that drive this choreography:
- A real geometrid takeoff starts with the WINGS: a distinct wing-raise phase (wings lift from the resting posture into a raised V) happens while the feet are still planted; NO body translation and NO body pitch before the raise is substantially complete.
- The raise occupies ~79% of the interval from first movement to liftoff (measured on *Nemoria lixaria*, 6000 fps). Then a short planted hold with a small forward pitch, then the first downstroke and liftoff in quick succession.
- Both wings sweep together (no alternating stroke). Fore/hind wings are one fan mesh here, so phase lock is structural — keep it.
- No periodic perched rhythm was measured anywhere (candidate intervals were irregular, CV 0.5-1.16, and the events themselves unconfirmed). Perched micro-motion in this choreography is therefore small, brief, and IRREGULARLY timed — and the items the dossier marks as pure invention (antenna twitches, reaction to the other moth, a 1-frame left/right lead) are kept small and labelled as artistic in the results file.

## Hard constraints (unchanged from proto3)
Blender 4.2.3, Eevee 64 samples, motion blur on, 832x480, 48 frames @24fps (2.000 s), production camera truly locked, background plate unchanged, left moth out of frame by f24 (proto3 state: ~85% out at f24, fully out by f47 — preserve that exit envelope), right moth perched for the whole clip, abdomen covered while perched and exposed only during takeoff, NO mesh-gen AI, NO photo sprites. Deterministic: the same code must render bit-identical frames on re-run (no wall-clock or unseeded randomness anywhere).

## Choreography v2 spec (frames are 0-based, 48 frames total)
Implement the animation so the following timeline holds. Where the existing code computes wing/body animation procedurally, restructure its parameters/keyframes to match; keep the left moth's exit PATH (positions from f18 to f47) as in proto3 — you are re-timing the wing/body behaviour inside the same envelope.

SEED = 20261110. Every irregular value below is a fixed literal from this brief (they ARE the seeded values); if you need any additional jitter, derive it from a literal table you write into the code and into the results file — never from runtime randomness.

### Left moth (the one that leaves)
- f0-f7 perched: rest pose as proto3 (wings flat, ~3 deg). Two micro events only:
  - f3: antenna twitch — both antennae rotate ~2 deg and back over 2 frames (artistic, dossier-labelled).
  - f6: wing settle — both wings dip 1.5 deg over 1 frame and return over the next (barely visible).
- f8: wing raise BEGINS. From here to liftoff the moth's root position must not translate and the body must not pitch (wings only).
- f8-f15: RAISE phase (8 frames): wing dihedral rises from rest (~3 deg) to raised (~+28 deg, mirroring the measured +26 deg span-inclination change), smooth ease-in-out, both wings symmetric.
- f16-f17: HOLD at the raised posture (angle within +/-1 deg). At f17 a small forward body pitch (2-3 deg) begins while still perched.
- f18: LIFTOFF: body leaves the bark and the FIRST DOWNSTROKE starts in the same frame. On this first downstroke only, the left wing leads the right wing by exactly 1 frame (artistic asymmetry, dossier-labelled).
- f18-f24: powered exit: 3-4 wingbeats with SEEDED irregular periods — beat period factors [1.00, 0.90, 1.12, 0.94] applied in order to the base flap period (proto3 flap_hz 9.5). Body pitch increases into a climb attitude; translation follows proto3's exit path so the moth is ~85% out of frame at f24 and fully out by f47. After f24 keep the existing proto3 behaviour for the remaining frames (it is off-screen anyway).
- Abdomen: stays covered by the seam/hindwing structure while perched; becomes visible only once the wings raise/cup during f15+ and in flight (as in proto3's f12 takeoff still).

### Right moth (stays perched the whole clip)
- Never translates; feet planted; rotations only. The existing slow breathing (wing angle drifting 3 -> 6.5 deg across the clip) CONTINUES underneath everything below.
- f5: antenna adjustment — both antennae rotate ~2-3 deg and settle over 3 frames.
- f13: wing micro-flick — both wings lift +2.5 deg over 2 frames and settle back over 2 frames (4 frames total). (Timing/amp are artistic; brief, small, non-periodic — dossier-consistent.)
- f21: single-antenna adjustment — RIGHT antenna only, ~3 deg over 4 frames, out and back (asymmetric, artistic).
- f31-f39: partial wing-raise cycle (after the left moth has gone; cf. the perched wing-sweep cycling observed in the Ant Lab *Hypagyrtis* segment — the behaviour class is observed, its use here as post-departure staging is artistic): wings raise +9 deg over f31-f34, hold f35-f36, settle back over f37-f39.
- Nothing else moves. In particular NO reaction flick at the moment of the left moth's liftoff (a reaction was NOT measured; do not stage one at f18).

### f0 identity check
Frame 0 must be identical to stage B's still: after your changes, render frame 0 and compare with `stageB-f0.png` in the same folder (mean abs pixel diff). It must be ~0 (the micro events start at f3/f5, and the breathing phase at f0 must be unchanged). If it differs, fix the phase/offsets until it matches, and record the diff value in the results file.

## Verification (implement, run, and report the numbers — do not just assert)
Write a small measurement script (reuse tools\ from stage B where useful) and produce `motion-results.md` in blender-moth-proto4 with:
1. Wing-angle and root-position series for the LEFT moth, frames 0-24 (from the animation parameters your code sets — log them per frame from the build script itself, e.g. by evaluating the same functions the keyframing uses): root position constant until f18 (report max |delta| in world units for f0-f17), wing angle at f8/f15/f16/f17/f18, pitch start frame.
2. Raise proportion: (liftoff frame - first-move frame) and (raise-end frame - first-move frame); ratio must land in 0.6-0.85 (spec gives 8/10 = 0.80).
3. Exit wingbeat intervals (frames between successive downstroke starts, f18-f24) as a list + their CV (must be > 0.05 — i.e. not metronomic).
4. Left/right wing phase: confirm identical angles every frame EXCEPT the first downstroke's 1-frame lead (report the max L/R angle difference outside f18-f20 and the lead at the first downstroke).
5. Right moth event log: frames and peak amplitudes of the f5/f13/f21/f31-39 events as rendered (from the same parameter log), and confirmation its root never moves (max |delta| over all 48 frames).
6. Determinism: render frame 24 twice in two separate Blender invocations; report both PNG md5 hashes — they must be equal.
7. f0 identity: mean abs diff vs stageB-f0.png (target ~0; report the number).

## Render and deliverables (in blender-moth-proto4)
- Full render + encode via the existing assemble.py flow, output `proto4.mp4`. Verify with ffprobe: 832x480, 24fps, 48 frames, 2.000 s, h264. Update the assemble output name/params if it still says proto3 (check build_params/assemble.py; proto number fields may still say 2 or 3 — set them to 4 where they feed file naming, and note it).
- Extract stills `f00.png`, `f12.png`, `f24.png`, `f47.png` from the final render (ffmpeg).
- `motion-results.md` as above, including: the full seeded event table as implemented, which items are measured vs artistic (copy the distinction from the dossier), verification numbers, and honest remaining issues.
- Append a "proto4 stage C (choreography v2)" section to STAGE-LOG.md. Use the jobkit checkpoint (D:\mutaro-work\jobkit\checkpoint.py) in this folder with units: choreo, verify, render, report — unit_done each.

## Rules
- Animation only. If you find you must touch geometry/material code to achieve the timeline, STOP that change, log it in motion-results.md as a blocker, and achieve the timeline within the animation system instead.
- If a verification target fails, adjust the animation ONCE and re-verify; if it still fails, report the failing numbers honestly rather than tuning silently.
- Run to completion; never wait for input.
