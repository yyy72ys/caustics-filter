# motion-results.md — proto4 STAGE C: choreography v2

Blender 4.2.3 LTS, Eevee NEXT 64 samples, motion blur on, 832x480, 48 frames @ 24 fps
(2.000 s), production camera locked, bark plate unchanged. Blender render time for the
full 48-frame clip: **165.6 s** (3.45 s/frame mean).

This stage changes **animation only**. Evidence that geometry / materials / contact
shadow / camera / lights / background plate are untouched: frame 0 of the shipped code
differs from stage B's `stageB-f0.png` by **mean |diff| 0.00058** (39 px of 399 360 above
2 levels) — and two renders of the *same* frame with the *same* code differ by
**0.00051** (27 px), i.e. the f0 difference is the renderer's own noise floor (§7).

---

## 1. What changed

| file | change |
|---|---|
| `build_moth.py` | the whole `ANIMATION` block was replaced (one block, lines "PROTO4 STAGE C"); `params` now reports `proto: 4`, `flap_hz 9.5`, the beat grid and the motion-log path. **Nothing else was edited** — no geometry, material, shadow-texture, decal, camera, light or plate line was touched. |
| `assemble.py` | output `proto3.mp4` → **`proto4.mp4`** (docstring + `OUT`); everything else unchanged. |
| `_dbg/build_params.json` | proto number field `proto: 2` → `4` (it is a parameter record, not a filename). |
| new tools | `tools/exitcap.py` (reads proto3's animated values out of the stage-B scene), `tools/motion_verify.py` (all §4 numbers), `tools/f0diff.py`, `tools/stills.py`, `tools/abdomen.py`, `tools/sheet.py`, `tools/cropzoom.py`, `tools/showlog*.py`. |

Deliberately left alone and noted: `BLEND_OUT` is still `proto2.blend` (the working scene
file, historically named; renaming it would orphan `tools/exitcap.py` and the stage-B
artefacts). It feeds no deliverable name.

### How the animation is now built

The old code used sparse keyframes with Bezier interpolation and a `sin()` flap whose
phase came from a fixed 9.5 Hz rate. v2 replaces it with **one pure function of the frame
number per channel**, evaluated per frame, keyed on every frame, and switched to LINEAR:

* `left_up(f)`, `left_up_lag(f)` — left-moth wing dihedral; `right_wing_up(f)` — right moth.
* `stroke_phase(f)` — the irregular beat grid (cumulative sums of the listed period factors).
* `_pchip()` (monotone cubic Hermite) — the left-moth exit path, so the position travels
  smoothly through proto3's envelope points without overshooting.
* `_piecewise()` (smoothstep between explicit frame/value keys) — every micro event.
* antenna rotation about the **base point** (pivot compensation on the object location) —
  a pure transform, computed from the tube's own first vertex ring. No mesh edit.

Consequence: **no wall clock, no unseeded randomness, no sub-frame sampling of a
frame-rate-dependent formula anywhere.** Two separate build runs emit a
**byte-identical** `_dbg/motion-log.json` (22 866 B, verified with `fc /b`).

The `_dbg/motion-log.json` written by the build is produced inside the same per-frame
loop that sets the keyframes, so every number in §4 is the value that was keyed, not a
re-derivation.

---

## 2. Seeded event table as implemented (SEED = 20261110)

`SEED` names the table; every "irregular" value is a **literal in `MOTION_TABLE`** in
`build_moth.py`. No additional jitter was needed, so no extra literal table was written.

### LEFT moth (leaves)

| frame | event | value as implemented | kind |
|---|---|---|---|
| 0–7 | perched rest | 3.00 deg, no translation, no pitch | measured (rest pose) |
| 3 | antenna twitch (both) | 0 → **2.00** (f3) → 1.00 (f4) → 0 (f5), 2 frames | **artistic** (dossier: 触角の動き = 純粋な創作) |
| 6 | wing settle | 3.00 → **1.50** (f6) → 3.00 (f7), 1 frame + return | artistic (below dossier resolution; 実測なし) |
| 8 | wing raise begins | 3.00 deg, root stationary | measured (raise phase) |
| 9–14 | raise | 4.39, 7.96, 12.84, 18.16, 23.04, 26.62 deg | measured |
| 15 | raise complete | **28.00** deg (+25.0 from rest) | measured (dossier +26 deg span-inclination → this is the 3-D dihedral) |
| 16–17 | planted hold | 28.00 / 28.00 deg (Δ = 0.000) | measured (short hold) |
| 17 | forward pitch begins | pitch **0.50** deg, still planted | measured (small fwd pitch while planted) |
| 18 | **liftoff + first downstroke** | wing to **32.00** deg (top of stroke), root z +0.015, pitch **2.50** deg | measured (liftoff, first downstroke) |
| 18 | first downstroke: left wing leads | leading wing starts at f18.0, lagging wing at f19.0 → **exactly 1 frame** | **artistic** (dossier: 最初の打ち下ろしの左右1フレーム差) |
| 18–28 | 4 beats, irregular | period factors **1.00, 0.90, 1.12, 0.94** on the 2.52632-frame base period → starts at **18.000, 20.526, 22.800, 25.629, 28.004** | measured mechanism (no periodic rhythm) + seeded factors |
| 18–24 | exit | proto3 path, `pitch → 20.45`, `roll → −26.22`, `yaw → +13.26` at f24 | measured (climb attitude) |
| 24–47 | off-screen | proto3 path continues (f28/f36/f47 = proto3 values exactly), flap continues | proto3 behaviour kept |

### RIGHT moth (perched throughout)

| frame | event | value as implemented | kind |
|---|---|---|---|
| 0–47 | breathing | 3.00 → 6.50 → 3.00 → 3.00 → 6.50 → 3.00 → 3.00 deg (proto3's own key list, smoothstepped; f0 = 3.00) | measured-ish (proto3 baseline, kept) |
| 5 | antenna adjustment (both) | 0 → **2.50** (f5) → 1.85 (f6) → 0.65 (f7) → 0 (f8), 3-frame settle | **artistic** |
| 13–17 | wing micro-flick | event 0 (f12) → **+2.50** (f15) → 0 (f17); 4 frames total; total wing angle peaks 5.50 (= 3.00 breath + 2.50) | **artistic** (dossier-consistent: short, non-periodic) |
| 21–24 | **right** antenna only | 0 (f20) → **3.00** (f22) → 0 (f24), 4 frames | **artistic** |
| 31–39 | partial wing-raise cycle | event 0 (f30) → **+9.00** (f34), hold (f36), 0 (f39); total peaks **12.91** = 3.91 breath + 9.00 | observed behaviour class (Hypagyrtis wing-sweeping); **its staging here is artistic** |
| 18 | — | **no reaction**: wing f17 3.000 / f18 3.000 / f19 3.000 deg | correct — a reaction was NOT measured, none staged |
| all | root | max |Δ| over 48 frames = **0.000000** world units | measured constraint |

### Measured vs artistic (distinction copied from MOTION-DOSSIER.md)

| rule | dossier basis | type |
|---|---|---|
| takeoff starts with the wings; no body translation/pitch until the raise is nearly done | 2 species measured | **measured** |
| raise occupies 6–8.5 tenths of first-move → liftoff (Nemoria 79 %) | Nemoria lixaria, 6000 fps | **measured** |
| after the raise: short planted hold + small forward pitch → downstroke → liftoff | measured | **measured** |
| first downstroke immediately after liftoff, both wings together | measured | **measured** |
| perched micro-motion: short, small, irregular — never periodic | periodicity NOT measured (candidate CV 0.505, unconfirmed; the 1.164 set rejected as codec artefact) | **contradicts nothing measured / artistic** |
| perched "small repeated wing lift" is allowed | Hypagyrtis observation (function uninterpreted) | **observed behaviour class** |
| antenna twitches, reaction to the neighbour, 1-frame L/R lead on the first downstroke | measurably absent (antennae move < 1 px; no such scene) | **pure invention, kept small** (as instructed) |
| all irregularity seeded/reproducible | production rule | **production rule** |

---

## 3. Frame-0 identity check

| render | mean abs diff vs `stageB-f0.png` | max | px > 2 levels |
|---|---|---|
| `c4-f0.png` (still) | **0.000508** | 11 | 33 (0.0083 %) |
| `renders/frame_0000.png` (from the full render) | **0.000579** | 11 | 39 (0.0098 %) |
| `det-f0b.png` vs `c4-f0.png` — **same code, second run (noise floor)** | 0.000511 | 12 | 27 (0.0068 %) |

Target ~0 → **met**: the f0 difference is smaller than the difference between two renders
of the identical frame, i.e. it is Eevee sampling noise, not a phase/offset change. The
parameters themselves are identical to proto3's f0 (wing 3.000 deg both sides, root at
(−0.60000, −0.03461, 0.00000), pitch/roll/yaw 0, right moth angle 3.000).

---

## 4. Verification numbers (all from `tools/motion_verify.py`)

### 4.1 LEFT moth — wing angle + root position, f0–f24

```
f   loc (x, y, z)                    pitch    roll     yaw   wing lead  wing lag
00  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
01  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
02  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
03  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
04  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
05  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
06  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +1.500    +1.500
07  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
08  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +3.000    +3.000
09  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +4.385    +4.385
10  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00     +7.956    +7.956
11  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +12.840   +12.840
12  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +18.160   +18.160
13  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +23.044   +23.044
14  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +26.615   +26.615
15  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +28.000   +28.000
16  -0.60000 -0.03461 +0.00000        0.00    0.00    0.00    +28.000   +28.000
17  -0.60000 -0.03461 +0.00000        0.50    0.00    0.00    +28.000   +28.000
18  -0.60000 -0.03461 +0.01500        2.50    0.00    0.00    +32.000   +28.000
19  -0.62612 -0.01749 +0.06707        5.33   -1.74    0.84    -25.387   +28.000
20  -0.69452 +0.02715 +0.16528        8.90   -6.14    2.99     +8.282   -15.713
21  -0.79023 +0.08920 +0.28913       12.70  -11.99    5.87     +8.282    +8.282
22  -0.89828 +0.15857 +0.41811       16.23  -18.10    8.92    -19.109   -19.109
23  -1.00372 +0.22514 +0.53171       18.99  -23.24   11.57    +28.896   +28.896
24  -1.09158 +0.27882 +0.60945       20.45  -26.22   13.26    -28.430   -28.430
```

* **Root constant until f18: max |Δ| over f0–f17 = 0.000000 world units.** First frame
  with any root delta = **f18** (+0.015 in z; x/y still at rest, so the moth visibly leaves
  the bark from f19).
* **Wing angle f8 = 3.00, f15 = 28.00, f16 = 28.00, f17 = 28.00, f18 = 32.00 deg**
  (rest 3.00). Hold f16→f17 delta = **0.000 deg** (spec ±1). f6 settle dip = **1.50 deg**.
* **Pitch start frame = f17** (0.50 deg planted, 2.50 at liftoff — inside the 2–3 deg the
  brief asks for at liftoff — then 20.45 deg by f24).
* f24 position/rotation are **numerically identical to proto3's** evaluated values
  (`(−1.09158, 0.27882, 0.60945)`, pitch 20.45, roll −26.22, yaw 63.26).

### 4.2 Raise proportion

| definition | value |
|---|---|
| first move / raise end / liftoff | f8 / f15 / f18 |
| (raise-end − first-move) / (liftoff − first-move) = 7/10 | **0.700** |
| counting the 8-frame raise as f8…f15 and its completion at the first hold frame f16 | **0.800** (= the brief's 8/10) |
| dossier value / spec band | 0.79 / **0.6–0.85** → **PASS** |

Both readings land in the band; nothing was tuned to hit the number.

### 4.3 Exit wingbeat intervals (f18–f24)

* base period = 24/9.5 = **2.52632 frames**
* downstroke starts inside f18–f24: **18.000, 20.526, 22.800**
* intervals: **[2.5263, 2.2737] frames**
* **CV = 0.0526 (population) / 0.0744 (sample, n−1)** → target **> 0.05 PASS**
  (only two intervals fit inside the window — the 3rd/4th beats start at 25.629 and
  28.004, i.e. after f24 but still in the takeoff — so the CV is computed on a 2-sample
  list and sits close to the gate by construction of the brief's own factors).

### 4.4 Left/right wing phase

* max |leading − lagging| **outside f18–f20 = 0.000000 deg** (identical every frame).
* inside the first downstroke: f18 **4.000**, f19 **53.387**, f20 **23.995** deg.
* the lead is exactly **1 frame** (leading stroke starts 18.0, lagging 19.0), and the two
  wings are **identical again from f21 to f47** (verified: all deltas < 1e-9).
* implementation note: the lagging wing's first stroke is compressed into the remainder of
  the same beat (18.526→20.526) so the pair re-syncs at the second beat instead of
  drifting a permanent frame apart.

### 4.5 RIGHT moth event log + root

| event | as rendered (total) | event component | spec |
|---|---|---|---|
| breathing | f0 3.00, min 3.00, max 6.50, f47 3.00 | — | 3 → 6.5 drift ✔ |
| f5 antenna both | peak **2.50** @ f5 (1.85 / 0.65 / 0 over 3 frames) | 2.50 | ~2–3 deg ✔ |
| f13 wing flick | peak **5.50** @ f15 = 3.00 breath + **2.50** | 2.50 | 2.5 deg ✔ (4 frames f13–f17) |
| f21 right antenna only | peak **3.00** @ f22 (0 @ f20, 1.50 @ f21/f23, 0 @ f24) | 3.00 | ~3 deg, 4 frames ✔ |
| f31–f39 wing raise | peak **12.91** @ f34 = 3.91 breath + **9.00**; hold f34–f36; 0 @ f39 | 9.00 | +9 raise, hold, settle ✔ |

**Right root max |Δ| over all 48 frames = 0.000000 world units** (never translates).
**No reaction at f18** (wing 3.000 / 3.000 / 3.000 for f17/f18/f19).

### 4.6 Determinism — **FAIL at the pixel level, pass at the parameter level**

Frame 24 rendered three times in three separate Blender invocations:

| file | md5 |
|---|---|
| `det-a.png` | `9718e7ff8fd70f73c5f626b1a1f7780b` |
| `det-b.png` | `df79d1177d47d8e2d4d64dcdaa8ebd78` |
| `det-c.png` | `7b51b39388bbd890025d1afaa6a127cb` |

* **md5 NOT equal.** Pairwise pixel differences: mean |d| 0.0017 / 0.0010 / 0.0010,
  max 28 / 28 / 10 levels, 74 / 34 / 44 px of 399 360 above 2 levels (0.019 % / 0.009 % /
  0.011 %).
* **The animation is bit-deterministic**: two full build runs write a byte-identical
  `motion-log.json` (`fc /b` → no differences, 22 866 B). Every animated value is a
  literal-table-driven pure function of the frame number; there is no `random`, no
  `time`, no frame-rate-dependent sampling.
* Cause of the pixel difference: Eevee's stochastic sampling (TAA + shadow jitter), not
  our code. Stage B recorded the same effect independently (`verify-f0` vs `stageB-f0`:
  0.001 / 49 px). **This target cannot be met with Eevee at any sample count**; it is
  reported as failed rather than worked around. (Rendering in Cycles with fixed seed
  would change the look — out of scope for an animation-only stage.)

### 4.7 Exit envelope

| frame | left-moth visible silhouette | bbox (x0,y0,x1,y1) |
|---|---|---|
| f0 | 29 211 px | 148,131 – 364,366 |
| f16 | 27 352 px | 150,140 – 359,358 |
| f24 | 25 922 px | **15,3 – 210,234** (clipped at the top/left edges) |
| f47 | **0 px** | — (fully out) |
| proto3 f24 (from `proto3.mp4`) | 31 514 px | 11,7 – 205,288 |

* Fully out by f47: **verified (0 px)**.
* Exit path preserved: f24/f28/f36/f47 positions and rotations are **exactly proto3's**.
* **The brief's "~85 % out of frame at f24" does not reproduce** — see honest issues.

### 4.8 Abdomen cover (stage-B C1 proxy, `tools/abdomen.py`)

| frame | centre-box median | dark fraction (< 0.65 × wing median) |
|---|---|---|
| f0 perched | 131.0 | **0.002** (covered, = stage B's C1 PASS) |
| f12 raising | 128.8 | 0.000 |
| f16 raised | 125.3 | 0.029 |
| f18 liftoff | 120.3 | **0.069** |
| proto3 f12 "takeoff" still (same method) | 96.0 | 0.016 |

The dark-centre fraction rises 30× once the wings leave the flat pose and is *larger* than
proto3's own reference still, and the f18 close-up shows the segmented abdomen between the
raised wings. The mechanism is purely the wing angle (the wings sit above the abdomen when
flat, swing clear when raised/cupped) — no geometry, seam or opacity change was made.

---

## 5. Render + deliverables

| item | value |
|---|---|
| frames | `renders/frame_0000.png … frame_0047.png`, 48/48 verified (832×480, > 100 kB) |
| encode | `assemble.py` → **`proto4.mp4`**, h264 crf 16, yuv420p, **252 750 B** |
| ffprobe | `{"codec_name":"h264","width":832,"height":480,"pix_fmt":"yuv420p","r_frame_rate":"24/1","duration":"2.000000","nb_read_frames":"48"}` |
| stills | `f00.png`, `f12.png`, `f24.png`, `f47.png` (extracted from `proto4.mp4` with ffmpeg, `-vf select=eq(n\,N)`) |

---

## 6. Honest remaining issues

1. **Determinism (brief §6) fails as specified.** PNG md5s differ between invocations
   (§4.6). The parameter series is byte-identical; the pixel difference is Eevee's
   stochastic sampler (~0.01–0.02 % of pixels, ≤ 28/255 on a handful). Not fixable inside
   an animation-only stage without changing the renderer.
2. **"~85 % out of frame at f24" is not what the inherited proto3 state does.** The f24
   root position/rotation are bit-identical to proto3's, so the envelope was preserved —
   but measured by silhouette area only **11.3 %** of the moth is outside the frame at f24
   (it is clipped at the top/left corner). Measured on **proto3's own f24 frame** the
   moth is likewise still mostly in frame (31 514 px visible vs 25 922 px now; the
   difference is the wing pose — proto3's f24 had wings at −4.4 deg spread, v2 has them at
   −28.4 deg mid-stroke — plus motion blur). I did **not** speed the exit up to chase the
   85 % figure, because the brief's operative instruction is to keep proto3's exit path.
   Flagged as a discrepancy in the brief/inherited state, not silently tuned.
3. **The 1-frame witness asymmetry is visible for ~2 frames.** At 9.5 Hz / 24 fps a
   one-frame lead is 53 deg of wing-angle difference at f19 (leading wing already at the
   bottom of the stroke, lagging wing still at the top). That is the literal consequence of
   the instructed 1-frame lead and it is labelled artistic; it is not smoothed away.
4. **Exit-interval CV is close to the gate** (0.0526 population / 0.0744 sample) because
   only two intervals fall inside f18–f24 at 9.5 Hz. The intervals are exactly the brief's
   factors; no extra irregularity was invented to raise the CV.
5. **The raised hold (f15–f17) shows slight motion-blur edge doubling** on the wings —
   the wings are still moving ~2.5 deg/frame at f15 and the 0.5-frame shutter catches it.
   Cosmetic; removing it would mean changing the shutter (a render parameter, not animation).
6. **Abdomen exposure is subtle from this top-down 85 mm camera.** The proxy moves
   0.002 → 0.069 and the f18 close-up shows the abdomen, but at 1× the read is a darker
   centre rather than a dramatic reveal — the same as proto3's f12 still (0.016). No
   geometry was touched to "improve" it.
7. **`BLEND_OUT` still writes `proto2.blend`** and `stage*-f0/png` naming is unchanged;
   only the deliverable name (`proto4.mp4`) and the `proto` parameter field were set to 4.
