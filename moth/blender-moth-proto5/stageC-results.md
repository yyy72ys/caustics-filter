# STAGE C — proto5 v3 (choreography v3) vs proto4

Session 3 / S3b completion, 2026-11-11. Working folder `blender-moth-proto5/`,
normative spec `../MOTION-DOSSIER-2.md` §1/§4. Single writer; the folder was idle.

**Build provenance.** Session-start file `build_moth.py` md5 **248545747cd8**
(90 709 B on disk — the brief quotes 90 698 B; the file is byte-identical to
`_dbg/build_moth.s2-accepted.py`, same md5, so the byte count in the brief is a
typo, not a difference) → after the choreography edit **8858284bd1d3df4ffb6f5289aacc3fe0**
(93 959 B), snapshot `_dbg/build_moth.s3b-v3.py`. The v3 md5 was re-checked
immediately before the 48-frame render and again after it: unchanged.

---

## 1. What changed (choreography only)

| # | change | where |
|---|---|---|
| 1 | `BASE_HZ` 9.5 → **8.0** (period 2.52632 → **3.0** frames at 24 fps); comment now carries the dossier §1 transfer argument | `build_moth.py` ANIMATION |
| 2 | **f21 right-antenna-only event DELETED**; **both-antennae sweep** keys (28, 0.0) (29, 3.0) (30, 0.0) added and applied to **both** antennae of the perched moth | `right_ant_sweep()` |
| 3 | partial wing raise (31, 0) → (37, +9.0) → **(39, 0)**, rest from f39 (6 f up / 2 f down = Xanthotype 75:25 shape) | `wing_raise_cycle` |
| 4 | unchanged byte-for-byte: f5 antenna event, f13 flick keys, breathing keys, f8-15 raise / f16-17 hold / f18 liftoff, exit path, factors [1.00, 0.90, 1.12, 0.94], 1-frame left lead | verified by log diff below |
| 5 | description string → `proto5 v3 (8.0 Hz transferred; Xanthotype-shaped perched events; dossier-2)` | `_dbg/build_params.json` |

Geometry, materials, the shadow stack and the S2 detail layers were **not**
touched. Frame-0 identity confirms it: the S2 accepted `final-f0.png` vs this
build's `renders/frame_0000.png` differ by mean |d| **0.0068** (285 px >2, all in
the pile, `tools/detcheck.py`) — inside the hair-noise band measured below.

---

## 2. Motion verification (v3) — **PASS**

From the build's own `_dbg/motion-log.json` (`_dbg/s3b_motion_check.py`,
log `_dbg/s3b-motion-check.log`). Note on tooling: `tools/motion_verify.py` is
proto4-shaped — its block 5 hard-codes the v2 event list and reads
`R_ant_right_deg`, which no longer exists — so it cannot run on a v3 log. The
checker re-implements its block 1-5 logic on the same log with the v3 field
names; the silhouette block (8) and the proto4 comparisons are re-run separately
in §3.1. The proto4 tool itself is left unmodified.

```
base_hz = 8.0      base_period = 3.0 frames @24 fps
beat_starts[:5] = [18.0, 21.0, 23.7, 27.06, 29.88]      (brief: same, max|d| 0.00000)
PASS BASE_HZ == 8.0 / period == 3.0 / factors [1.0, 0.9, 1.12, 0.94] unchanged
PASS f5 both-antenna event kept (+2.5 deg at f5)
PASS both-antennae sweep keys (28,0)(29,3)(30,0) present
PASS no right-antenna-only event at f19-f27 (f21 event deleted; max 0.0000 deg)
PASS f13 flick kept (+2.5 deg over breathing)
PASS raise peak +9.0 deg at f37; keys (31,0)/(39,0) -> 6 f up / 2 f down
PASS at rest from f39 onwards
PASS right moth root never translates
PASS rest 3.0 / f6 dip 1.5 / raise 28.0 @f15 / hold f16-f17 / liftoff f18
PASS exit path starts at f18, root fixed f0-f17
PASS 1-frame left lead at the first downstroke, resynced after
RESULT: ALL CHECKS PASS
```

Per-object check of the antenna change (loaded `proto2.blend`, `_dbg/s3b_ant_probe.py`):
`Moth_R_antenna_L` and `Moth_R_antenna_R` are **identical** at every sampled frame
and both read **+2.50° at f5**, **+3.00° at f29**, **0.00° at f19/f21/f23/f28/f30/f39** —
i.e. both antennae sweep, and the old f21 event is gone. `Moth_L_*` stays 0 at those frames.

Log diff vs the v2 baseline (`_dbg/motion-log.v2-baseline.json`):
frames **f0-f18 are identical**; differences start at f19 and are confined to
`L_wing_lead/lag`, `L_twist_lead/lag` (the new beat grid) and, on the perched moth,
`R_antennas` at f21-23 (deleted event) / f29 (new sweep) and `R_wing` at f31-38
(new raise shape). Nothing else moved.

---

## 3. Stage C parity table — proto4 vs proto5 v3

Frame-0 stills: `s3b-f0.png` / `s3b-ns.png` (--noshadow) / `s3b-sil.png`
(--silhouette), same route as every previous session. Tools: `measure.py`,
`final_metrics.py`, `bodymean.py`, `antenna.py`, `c1parity_s2.py`.

| metric | proto4 | **proto5 v3** | delta | note |
|---|---|---|---|---|
| A1 forewing L/W | 1.434 | **1.434** | 0 | exact |
| A2 fan IoU (best aligned) | 0.777 | **0.777** | 0 | exact |
| A2 render IoU (best aligned) | 0.780 | **0.779** | -0.001 | inside hair noise (see §6.3) |
| A2 photo own-mirror (reference) | 0.901 | 0.901 | 0 | the structural cap, see §6.1 |
| A3 designed / peaks / clear | 7 / 7 / 5 | **7 / 7 / 5** | 0 | exact |
| A4 body / wingspan | 0.345 | **0.345** | 0 | exact |
| D1 grain, segmented mask | 4.963 | **5.782** | +16.5 % | S2 Change 3 (kept), not this session |
| D1 grain, geometric outline | 4.620 | **4.948** | +7.1 % | S2 layer; S2 measured 4.962/4.992, so this is TAA spread |
| E3 geom 0-3 px | -12.87 % | **-11.76 %** | +1.11 pt | S2 Change 2 attribute, see §6.2 |
| E3 geom 3-6 px | -23.48 % | -21.80 % | +1.68 pt | ditto |
| E3 geom 6-12 px | -24.23 % | -23.02 % | +1.21 pt | ditto |
| E3 geom 12-24 px | -14.59 % | -14.20 % | +0.39 pt | ditto |
| E3 geom 24-60 px | -1.19 % | -1.18 % | +0.01 pt | unchanged |
| E3 geom 0-3 px darkest 5 % | -59.26 % | **-58.52 %** | +0.74 pt | guard -68 % never in danger |
| wing median (n≈32 k) | (111,156,148) | **(111,156,148)** | exact | wing materials untouched |
| body footprint mean (2 273 px) | (93.9,121.0,86.6) | **(108.9,151.2,138.7)** | S2 colour | IN box R105-115/G148-158/B138-150 |
| antenna cores R / L | (67,67,60) / (66,65,55) | **(157,149,127) / (153,144,122)** | S2 colour | photo median (157,146,131) |
| C1 dark-centre rate (0.55 / 0.65 × wing median) | 0.002 | **0.000 / 0.000** | stricter | also 0.000 at 0.75 |
| C1 dark patch (<0.80 × wing median) | 4 px (99.7 % overlap) | **0 px (criterion vacuous)** | S2 pile/colour | documented S2 gap 2 |
| 48-frame render wall | 165.6 s | **182.9 s** (3.81 s/frame, 1.10×) | +17.3 s | same as S2's 182.7 s alone |

The body-colour and antenna rows are the S2 Change 2 result, not a v3 regression;
they are in the table because the brief asks for proto4 → proto5 side by side.
Every frame-0 metric **reproduces the merged-build baseline within TAA noise**:
A1/A2-fan/A3/A4 exact, wing median exact, C1 exact, body mean -0.1/-0.1/0.0,
D1-geometric -0.9 %, E3 bands ≤0.16 pt (0-3 px -11.76 vs -11.73, 3-6 -21.80 vs
-21.88). Nothing moved materially.

### 3.1 Exit check (rendered frames)

Same silhouette detector as proto4 block 8, run on **both** builds
(`_dbg/s3b-exitcheck.py`, `_dbg/s3b-exit-series.py`):

| frame | proto4 left silhouette | v3 left silhouette |
|---|---|---|
| f0 | 29 212 px | 29 393 px |
| f22 | 27 754 px (74,48-273,294) | 27 469 px (75,48-273,293) |
| **f24** | 25 883 px (15,3-209,234) | **27 976 px (14,63-208,312)** |
| f30 | 4 027 px | 2 756 px |
| f36 / f40 / f47 | 0 / 0 / 0 px | **0 / 0 / 0 px** |

* **Full exit at f47: verified 0 px**, three independent ways: silhouette
  detector 0 px (also 0 px from f36 on); `frame_0047.png` vs `frame_0040.png`
  (both moth-free in the left half) differ by max 3 levels, mean 0.296, 102 px
  >2 and **0 px >8** — i.e. TAA/hair noise plus the last shadow fade, no moth;
  and the bark-plate diff in the left half matches the never-moth control patches
  (median 4.0 vs 3.0/3.0), so the plate's global tone offset, not a moth, makes
  that number.
* **f24 out-of-frame fraction = 4.8 %** by the area-ratio method (proto4: 11.4 %).
  **Investigated, not silent:** the root position *and* pitch/roll/yaw at
  f24/f28/f36/f47 are **identical** between v2 and v3 logs (bit-equal), and no
  part of the f24 bbox is off-frame in either build (off-frame extent 0 px in
  both). The metric is pose-sensitive: at f24 the v3 beat grid (18 / 21 / 23.7)
  puts the wings at **+27.1°** (top of the stroke) where v2 (18 / 20.53 / 22.8)
  had them at **-28.4°** (bottom) — the same 20.45° pitch / -26.22° roll shows
  more dorsal area, which raises the visible fraction from 88.6 % to 95.2 % and
  moves the bbox from y3-234 to y63-312. The exit **path** is unchanged; the
  number is a wing-phase artefact, and it does **not** indicate a slower exit
  (both builds are fully out by f36).

---

## 4. Render + assembly

| item | value |
|---|---|
| render | `python run_blender.py build_moth.py --render`, Eevee 64 samples, motion blur, 832x480, 48 frames |
| wall time | **182.9 s** (3.81 s/frame; proto4 165.6 s = 1.10×; S2 accepted build 182.7 s alone) |
| frames | `renders/frame_0000…0047.png`, 48/48 verified (832x480, >100 kB) |
| encode | `assemble.py` → **`proto5.mp4`**, h264 crf 16, yuv420p, 250 366 B |
| stills | `f00.png`, `f12.png`, `f24.png`, `f47.png` copied from `renders/` |

ffprobe block (`ffprobe -v error -select_streams v:0 -count_frames -show_entries
stream=codec_name,width,height,r_frame_rate,avg_frame_rate,duration,nb_read_frames,pix_fmt,bit_rate
-show_entries format=duration,size,format_name -of json proto5.mp4`):

```json
{"streams": [{"codec_name": "h264", "width": 832, "height": 480,
              "pix_fmt": "yuv420p", "r_frame_rate": "24/1", "avg_frame_rate": "24/1",
              "duration": "2.000000", "bit_rate": "995664", "nb_read_frames": "48"}],
 "format": {"format_name": "mov,mp4,m4a,3gp,3g2,mj2", "duration": "2.000000",
            "size": "250366"}}
```

832x480 ✔ 24 fps ✔ 48 frames ✔ 2.000 s ✔

---

## 5. What the 8 s clip now does

* f0-f17 — perched pair, unchanged pose; perched moth breathes (3 → 6.5°), f5 both
  antennae +2.5°, f13 flick +2.5°.
* f18 — left moth liftoff with the proto4 measured raise/hold/liftoff sequence and
  the 1-frame left lead; unreliable-but-seeded beat grid **8.0 Hz**
  (3.0 f base, factors 1.00/0.90/1.12/0.94) → beats 18, 21, 23.7, 27.06, 29.88…
* f28-f30 — perched moth: both antennae sweep +3° (reflex order, before the wing).
* f31-f39 — perched moth: partial wing raise 0 → +9° → 0, slow open (6 f) / fast
  close (2 f); rest from f39. No f21 antenna event any more.
* f36-f47 — left moth fully out of frame; perched moth alone to the end.

---

## 6. Honest remains

1. **A2 is still 0.777 fan / 0.779 render against a 0.90 gate, and cannot get
   there by tuning.** The photo's own left/right mirror under best alignment is
   **0.901**, i.e. the gate is at the reference's own asymmetry ceiling
   (proto4 recorded the same cap: "~0.90-0.93 from the photo's own asymmetry").
   Reaching 0.90 against the photo would require reproducing that asymmetry
   rather than a symmetric fan, which contradicts the A2 construction itself.
   Carried forward unchanged; not a v3 regression (exact match to proto4).
2. **The contact-shadow layer stays reverted (S2 Change 4 negative result).**
   A zero-offset per-moth contact decal cannot darken the E3 0-3 px band because
   the rendered moth covers 100 % of the geometric outline — the layer's only
   exposed part is its outer tail (measured: 380 of 3 042 band px changed,
   +0.37 pt at alpha 0.40, +1.05 pt even at alpha 1.0, gate ≤ -15 %). The v3
   0-3 px band is **-11.76 %**, and ~1.1 pt of the proto4 gap is the brighter
   body/pile (S2 Change 2) raising the `--noshadow` reference, not the shadow:
   the absolute drop is unchanged (S2 A/B: -13.37 levels vs -13.36 at S1).
3. **The hair pile is not bit-deterministic.** Two renders of the *same* build
   differ by mean |d| **0.107**, max 55, 3 778 px >2 (0.95 % of the frame), bbox
   x210-627 y198-278 — all inside the pile. The parameter log is byte-identical
   and the motion log is reproducible, so this is Eevee TAA on sub-pixel strands,
   inherited from proto4. Adjacent frames will shimmer very slightly in the pile;
   judge determinism on the log + mean |d|, not on PNG md5. (The S2-vs-v3 f0 pair
   happened to land at 0.0068 mean, i.e. the noise varies run to run within the
   band 0.007-0.11.)
4. **Perched antenna behaviour is still unmeasured in footage** (3 data sets).
   Both the f5 event and the new f28-f30 sweep are **TRANSFERRED** in *order*
   (Manjunath et al. 2025: antennae → legs → wings) with invented amplitude;
   the sweep's timing is chosen to sit immediately before the wing raise, which
   is the only part the dossier supports. If footage ever resolves an antennae
   displacement, these two events must be re-measured, not defended.
5. **The perched events are a behaviour-class transfer, not this species.**
   The 6 f/2 f open:close asymmetry is the *shape* of Xanthotype's perched
   open/close cycle (75:25) with the amplitude shrunk to the partial range
   (+9°, proto4's number); the cycle period itself (8.3 s playback, speed-ramped
   → real time unknown) is deliberately not transferred. The dossier's §4 table
   labels each event; the same labels are in the build's choreography comment.
6. **The f24 out-of-frame number is pose-sensitive** (4.8 % now vs 11.4 % in
   proto4) even though the exit path is bit-identical — see §3.1. Anyone
   comparing this metric across builds must compare root position, not this ratio.
7. **8.0 Hz is still a band transfer, not a measurement.** Morpho 5.9 ± 1.7 Hz,
   winter geometrids ~4 Hz and f ∝ √(m/A_f) place a large low-wing-loading
   geometrid in the 5-9 Hz band; 8.0 Hz was chosen inside it. No footage of
   *this* species at a verifiable frame rate exists (dossier §5.2), so the
   number must not be quoted as measured.
8. **Cosmetic leftover:** `_dbg/build_params.json` still carries `"proto": 4`.
   It is an unrelated lineage field, left untouched because the brief restricted
   the session to the choreography parts; the `choreography` string is correct.
9. **E3's segmented variant is still unusable for comparisons** (S2 gap 3): its
   mask now includes the brighter thorax/pile. Geometric only, as before.

---

## 7. Artefacts

| kind | files |
|---|---|
| build | `build_moth.py` (md5 `8858284bd1d3df4ffb6f5289aacc3fe0`), snapshots `_dbg/build_moth.s3-baseline.py` (= accepted S2), `_dbg/build_moth.s3b-v3.py` |
| verify scripts | `_dbg/s3b_motion_check.py`, `_dbg/s3b_ant_probe.py`, `_dbg/s3b_exitcheck.py`, `_dbg/s3b_exit_series.py`, `_dbg/s3b_f24pose.py` |
| logs | `_dbg/s3b-motion-check.log`, `_dbg/s3b-ant-probe.log`, `_dbg/s3b-exitcheck.log`, `_dbg/s3b-exit-series.log`, `_dbg/s3b-f24pose.log`, `_dbg/s3b-measure.log`, `_dbg/s3b-final.log`, `_dbg/s3b-body.log`, `_dbg/s3b-ant.log`, `_dbg/s3b-c1.log`, `_dbg/s3b-det-still-vs-render.log`, `_dbg/s3b-det-s2-vs-v3.log`, `_dbg/s3b-render48.log`, `_dbg/s3b-assemble.log`, `_dbg/s3b-build-f0.log` |
| motion | `_dbg/motion-log.json` (v3), `_dbg/motion-log.v2-baseline.json` (v2, for the diff) |
| measures | `_dbg/measure-s3b.json`, `_dbg/build_params.json` |
| media | `proto5.mp4`, `f00/f12/f24/f47.png`, `renders/frame_0000…0047.png`, `s3b-f0/ns/sil.png`, `_dbg/p4-f{00,20,22,24,28,30,36,40,47}.png` (proto4 frames for the same-method comparison) |
