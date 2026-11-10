# proto5 STAGE-LOG

Append-only log for `blender-moth-proto5`. Kept in the proto4 STAGE-LOG format:
one `## ` section per session, bullets + tables, numbers first, gaps stated.
proto5 starts as a verbatim copy of proto4 (frozen); proto4 and `moth-motion/`
are read-only reference.

## S1 — setup: copy + baseline parity (BRIEF-S1 Task 1)

Source: `robocopy blender-moth-proto4 blender-moth-proto5 /E` (incl. `tools/`,
`_dbg/`, `renders/`, `__pycache__`). Blender 4.2.3 LTS, Eevee NEXT 64 samples,
motion blur, 832x480, production camera. No ordering/clock dependence added.

**Copy verification**
- proto4: **264 files / 85,020,447 bytes**; proto5: **264 files / 85,020,447 bytes**;
  Compare-Object over (relative path + size): **identical, 0 differences**.
- proto4 was not written to: 0 files with a timestamp newer than the copy start;
  newest proto4 timestamp 2026-11-10 21:19:20 (pre-copy).
- Code identity: `build_moth.py`, `assemble.py`, `run_blender.py` md5-equal to proto4
  (`ec4258f9…`, `c71f840a…`, `22330faa…`) — the baselines below were rendered from the
  unmodified copied code.

**Baseline stills** (route: `python run_blender.py build_moth.py -- still <name>.png` + flags)

| still | flags | exit | size (B) | wall |
|---|---|---|---|---|
| `base-f0.png` | — | 0 | 689,454 | 10.3 s |
| `base-ns.png` | `--noshadow` | 0 | 693,613 | 9.8 s |
| `base-sil.png` | `--silhouette` | 0 | 156,875 | 5.8 s |
| `base-obl.png` | `--oblique 35` | 0 | 223,092 | 11.0 s |

(proto4 final stills for comparison: stageB-f0 689,453 / stageB-ns 693,601 /
stageB-sil 156,875 / stageB-obl 223,076 — deltas are Eevee TAA noise.)

**Parity table — proto5 baselines vs proto4's reported values**

| metric | proto4 reported | proto5 measured | tool | parity |
|---|---|---|---|---|
| A1 forewing L/W | 1.434 | **1.434** | `tools/measure.py` | exact |
| A2 fan IoU (best aligned) | 0.777 | **0.777** | `tools/measure.py` | exact |
| A2 render IoU (best aligned) | 0.780 | **0.780** | `tools/measure.py` | exact |
| A2 photo own-mirror (reference) | 0.901 | **0.901** | `tools/measure.py` | exact |
| A3 designed / peaks / clear | 7 / 7 / 5 | **7 / 7 / 5** | `tools/measure.py` | exact |
| A4 body / wingspan | 0.345 | **0.345** | `tools/measure.py` | exact |
| D1 grain, segmented mask | 4.96 | **4.963** | `tools/measure.py` | +0.06 % |
| D1 grain, geometric outline | 4.62 | **4.620** | `tools/final_metrics.py` | exact |
| E3 geometric 0-3 px | -12.87 % | **-12.87 %** | `tools/final_metrics.py` | exact |
| E3 geometric 3-6 px | -23.48 % | **-23.48 %** | `tools/final_metrics.py` | exact |
| E3 geometric 6-12 px | -24.23 % | **-24.23 %** | `tools/final_metrics.py` | exact |
| E3 geometric 12-24 px | -14.59 % | **-14.59 %** | `tools/final_metrics.py` | exact |
| E3 geometric 24-60 px | -1.19 % | **-1.19 %** | `tools/final_metrics.py` | exact |
| E3 geometric 0-3 px darkest 5 % | -59.26 % | **-59.26 %** | `tools/final_metrics.py` | exact |
| E3 segmented 0-3 / 3-6 / 6-12 / 12-24 / 24-60 px | -3.43 / -16.53 / -27.49 / -17.53 / -1.74 | **-3.43 / -16.51 / -27.49 / -17.53 / -1.74** | `tools/e3seg.py` | 3-6 px -0.02 pt, rest exact |
| E3 segmented 0-3 px darkest 5 % | -55.75 % | **-55.75 %** | `tools/e3seg.py` | exact |
| wing median colour (n≈32 k) | 111,156,148 | **111,156,148** | `tools/measure.py` | exact |
| body mean colour | ~94,121,87 | **93.9,121.0,86.6** | `tools/thorax.py` | exact to 1 level |
| thorax blob area / bbox / aspect / PCA | 2,276 px / 59x53 / 1.113 / 1.358 | **2,273 px / 59x52 / 1.135 / 1.358** | `tools/thorax.py` | area -0.1 %, PCA exact, bbox 1 row shorter |
| C1 thorax-footprint overlap | 99.7 % | **99.7 %** | `tools/c1parity_s1.py` | exact |
| C1 dark-centre rate | 0.002 | **0.000** (0.55 and 0.65 × wing median) | `tools/c1parity_s1.py`, `tools/c1probe_s1.py` | **not reproduced — see below** |

**C1 dark-centre rate — honest gap.** The brief's 0.002 is not reproduced by any
threshold in either centre box. Proto5 measures **0.0000** of a 36x36 box at the
silhouette centroid below 0.55/0.60/0.65 × wing median, and **0.0000** of the 40x40
c1audit box at the same thresholds; the sweep only leaves zero at 0.70 × median
(0.0170 / 0.0213, i.e. 22 / 34 px of 1296 / 1600). 0.002 would be 2-3 px, so the
closest honest statement is: the C1 criterion (no dark centre) **passes more
strictly than proto4's reported 0.002**; the exact 0.002 definition is unknown and
is not claimed. Full sweep in `_dbg/s1c1probe.log`.

**Tooling note (deviation, stated).** `tools/c1audit.py` and `tools/c1shape.py`
hardcode `T = ...\blender-moth-proto4` and the filenames `base-f0.png` /
`stageB-f0.png`, so they cannot be aimed at the proto5 stills without editing them.
They were therefore **ported verbatim** (same thresholds, boxes, morphology,
definitions; only the folder and the input path are un-hardcoded) to
`tools/c1parity_s1.py`, and `tools/c1probe_s1.py` sweeps the C1 thresholds. The
proto4 scripts themselves are unmodified. All other tools used
(`measure.py`, `final_metrics.py`, `e3seg.py`, `thorax.py`) take their paths from
their own location or the command line and ran unmodified.

Logs: `_dbg/s1render-{f0,ns,sil,obl}.log`, `_dbg/s1measure-base.log`,
`_dbg/s1final-metrics.log`, `_dbg/s1e3seg.log`, `_dbg/s1thorax.log`,
`_dbg/s1c1.log`, `_dbg/s1c1probe.log`, `_dbg/measure-base.json`.

## S1 — Part A motion footage: abundant-footage species (BRIEF-S1 Task 2)

Folder `moth-motion2/` (sibling of the frozen `moth-motion/`). Full report:
`moth-motion2/analysis-report.md`, copied to `_dbg/motion2-analysis-report.md`.
Machine-readable headline numbers: `moth-motion2/events-summary.json`.

**Resolution — obtained, not assumed.** Re-tested the YouTube format ceiling with an
11-client matrix (`client_matrix2.py`, `cm2.txt`): only `android` succeeds and it exposes
**640×360 only** (format 18); `android_vr` demands sign-in, `ios/mweb/web/web_safari` report
no format, `tv`/`tv_embedded` are rejected. **8 of 9 clips are 640×360; `tX7XwPGmv0g` is
480×360 (its own source size).** No cookies/PO tokens were used, so this is the same ceiling
proto4 documented — recorded per clip rather than papered over.

**Clips (9 downloaded of 10 attempted; `N_0ODUyywO8` cabbage white refused: "Requested
format is not available")**

| id | channel | species (confidence) | fps | res | role |
|---|---|---|---|---|---|
| `I5L3_UW1flY` | Ant Lab | geometrid *Xanthotype* sp. (MEDIUM, channel's tentative list) | 30.00 | 640×360 | geometrid perched + takeoff |
| `JQL25_hoQ1k` | Ant Lab | blinded sphinx *Paonias exaecata* (HIGH by label) | 29.97 | 640×360 | sphinx raise + takeoff |
| `nzU2WPOBaAI` | Ant Lab | blue morpho (HIGH by chapter, species unverified) | 29.97 | 640×360 | basking cycles |
| `B1ntmKcGqYI` | Jake Tobin | monarch (HIGH) | 30.00 | 640×360 | takeoff (visual only) |
| `Z-FEfNpZUL0` | JcmdiStockFootage | monarch (HIGH by title) | 29.97 | 640×360 | takeoff, "300 fps" title only |
| `y-cZ3rZF8bA` | Kaleidoscope Garden | queen *Danaus gilippus* (MEDIUM-HIGH) | 30.00 | 640×360 | basking |
| `2dGUKCPC8AU` | JcmdiStockFootage | painted lady (MEDIUM-HIGH by title) | 29.97 | 640×360 | pre-takeoff fan + liftoff |
| `MhsvrLinGus` | NanoNaturalist | white-lined sphinx *Hyles lineata* (MEDIUM) | 29.97 | 640×360 | shiver → takeoff (bounded) |
| `tX7XwPGmv0g` | naturalistguy | sphinx, family only (LOW) | 30.00 | **480×360** | shiver onset/duration |

**Measured headline numbers** (all playback-time; no real-time conversion — Ant Lab's
6000 fps label is speed-ramped, proto4's trap re-applied):

- Geometrid perched, left half: silhouette area **7,840 → 19,855 px** (×2.53), vertical
  extent **116 → 217 px** (+87 %), one open-close cycle **≈8.3 s playback**; right half:
  area **8,767 → 23,031 px** (×2.63), open/closed/open runs 1.33/5.67/2.83 s, **no period**.
- Morpho basking: 3 open runs, periods **18.42 s / 24.82 s** (mean 21.62, CV 0.148, n=2
  intervals), open mean **12.55 s**, closed mean **8.54 s**, duty **0.59 : 0.41**; horizontal
  wing extent **405 → 596 px**.
- Sphinx takeoff: wing raise **288.15–288.56 s**, liftoff **≈288.81 s**, first full
  downstroke **289.09–289.16 s** (≈0.66 s raise-to-liftoff, visual + motion energy).
- Geometrid takeoff: left first airborne **j637 ≈289.20 s**, right **j635 ≈289.07 s**;
  wings sweep down together (visual).
- Painted lady: wing opening **66.13–66.53 s**, airborne **≈67.33 s**, first wingbeat
  **67.73–68.13 s** → fan-to-liftoff **≈1.2 s playback** (≈0.12 s only if the title's
  300 fps is literal — flagged, not verified).
- Shiver (`tX7XwPGmv0g`, 30 fps, no slowdown stated): onset **≈2.3 s**, peak
  **3.43–3.57 s**, departure **3.57–3.63 s** → **≈1.1–1.3 s onset-to-gone**, wing-region
  motion energy **≈6× baseline**, monotone growth. Frequency **not measurable** (30 fps alias).

**Honest gaps (nothing filled from general knowledge):** antenna behaviour (individual
antennae not resolvable in 2–3× studio crops), abdomen pumping, real-time wingbeat frequency,
3-D dihedral / stroke plane, short-interval perched flick rhythm, two-moth interaction
(clip 1 is split-screen of separate shots), leg kinematics, and species identity for the
geometrid/morpho/warm-up sphinx. Full list in report §7.

**Untouched, verified:** `blender-moth-proto4` (264 files, newest 21:19:20) and `moth-motion/`
(newest 20:35:40) — 0 files newer than the session start. All work is in `moth-motion2/` and
`blender-moth-proto5/`. No wall-clock dependence in any measurement (frame indices only).

## S2 — detail build (BRIEF-S2, Session 2 of 3)

Working state = the verified proto5 copy from S1. Every change below was applied
ONE element at a time, rendered, measured with the existing tools, and kept or
reverted on its gate. Baseline for all comparisons = the S1 parity values
(A1 1.434; A2 fan 0.777 / render 0.780; D1 segmented 4.963; E3 geometric bands
0-3 px -12.87 % and darkest-5 % -59.26 %; wing median (111,156,148); body mean
(93.9,121.0,86.6); C1 99.7 % / 0.000). All stills Eevee 64 samples, motion blur,
832x480 (the two 1664x960 stills are Change 3's 2x test, as the brief requires).

New debug views added to `build_moth.py` (default off, additive only):
`--nohair` (Change 1 measurement), `--noscale` (Change 3 before/after),
`--res2x` (1664x960 stills). A `--nocontact` view was used only for the
Change-4 mechanism test and was removed again with the revert.

### S2 Change 1 — thorax particle hair (DETAIL-RESEARCH §2, option a) — KEPT

Both moths, thorax + head (they share the body-material path).

| parameter | value | note |
|---|---|---|
| system | legacy particle hair (`PARTICLE_SYSTEM`, type HAIR) | fixed seeds: thorax 11, head 23; no children (see traps) |
| count | 7 000 (thorax) / 3 500 (head) | strands |
| length | 0.030 / 0.015 object units | = 0.43 x the part's mean local radius (thorax ~0.07, head ~0.035); inside the brief's 0.25-0.45x |
| measured strand | chord 0.0207 / 0.0104 obj = 0.024 / 0.012 world = **6.6 / 3.3 px** at 832 | depsgraph readback, `_dbg/probe_hairlen4.py` |
| strand radius | 0.0032 obj (**~1.0 px**) | `radius_scale`; root 1.0 -> tip 0.0 |
| taper | `shape = -0.35` | rooted thick, needle-thin tip |
| roughness / clumping | WAVE kink amp 0.006, freq 2.5, extra 2 steps; `roughness_1` 0.004 | reads as radiating pile, not spikes |
| colour | `HairThoraxMat` / `HairHeadMat`: body family, tips one step lighter | Eevee Next returns 0 from Hair Info `Intercept` (probed, `_dbg/probe_hair.py`), so root->tip is approximated by object-space radius from the part centre (stated approximation, not a true along-strand gradient) |

**Two Blender 4.2 traps, both found by probe rather than guessed**
1. Writing `normal_factor` / `tangent_factor` in advanced-hair mode RECOMPUTES
   `hair_length` (0.024 -> 4.0); length is therefore written LAST
   (`_dbg/probe_hairseq.py`).
2. `factor_random > 0` scales the per-strand random velocity by the 4.0 default,
   so the *evaluated* strand ignored `hair_length` and came out **0.88 object
   units (11x the thorax)** — a fuzz cloud that covered the whole frame — even
   though the RNA read back 0.024. With `factor_random = 0` the strand equals
   `hair_length` exactly. Interpolated children showed the same bug, so children
   are off and the clumping/roughness is the WAVE kink
   (`_dbg/probe_hair_iso2.py`, `_dbg/probe_hair_iso3.py`).

**Gate (a) silhouette edge reads fuzzy/soft — PASS.** `tools/thfuzz.py`
compares the render against the identical scene with `--nohair` (oblique-35
view, thorax edge against the flat grey world):

| measure | without hair | with hair |
|---|---|---|
| changed px outside the skin outline (pile standing proud) | — | **202 px**, mean 1.27 px out, p90 3.0, max 4.0 |
| changed px total (incl. over the skin) | — | 3 723 |
| straight-view changed px (f0, both moths) | — | 7 438 |
| ray 20-80 % transition width | 5.75 px | 4.79 px |

The ray-width number is **not** a usable discriminator (its rays cross wing
edges as well as the thorax, so the median is wing-dominated and moves with
noise); the direct fringe numbers plus the visual check are the evidence.
Visual, 4x/6x crops (`_dbg/h1-obl-th*.png`, `_dbg/h1-f0-th*.png`): with hair the
outline is "fuzzy and irregular, a dense coat of visible hairs radiating
outward... individual hairs appear as distinct short strands rather than only a
soft blur"; the same crop without hair is "broadly capsule/oval, reasonably
clean and continuous". It reads as pile, not spikes.

**Gate (b) C1 parity — PASS on the audit box; the overlap criterion became
vacuous.** The S1 audit (`tools/c1parity_s1.py`) defines the thorax footprint by
colour — `(g-b) > 14`, olive body vs cyan wing — and the dark patch by
`lum < 0.80 x wing median`. The pale pile removes both (blob 2 273 px -> 183 px,
then 65 px after Change 2; dark patch 2 034 px -> **0 px**), so it can no longer
compute the overlap. `tools/c1parity_s2.py` keeps every S1 rule (same 40x40 box
at (563,245), same 0.55/0.65/0.80 thresholds, same morphology) and freezes only
the reference region to the Session-1 footprint (`_dbg/c1_footprint_s1.npy`,
2 273 px, reproduces the S1 baseline mean (93.9,121.0,86.6) exactly):

| C1 quantity | S1 baseline | S2 final |
|---|---|---|
| dark-centre frac < 0.65 x wing median (audit box) | 0.000 | **0.000** |
| dark-centre frac < 0.55 x wing median (audit box) | 0.000 | **0.000** |
| same, frozen thorax footprint | 0.014 / 0.007 | **0.000 / 0.000** |
| dark patch area (< 0.80 x wing median) | 2 020 px | **0 px** |
| overlap dark patch / thorax footprint | 100.0 % | **n/a (0 px dark area)** |

So "no dark centre" is met more strictly than in S1; "overlap >= 99 %" is
vacuous because there is no dark region left to overlap. Stated as a gap, not
claimed as a pass.

**Gate (c) A2 render IoU — PASS.** 0.779 (baseline 0.780, floor 0.77);
fan IoU 0.777 (exact). Hair adds white pixels to the silhouette render and moves
IoU by 0.001-0.002.

**Gate (d) full 48-frame render — PASS.** `--render`, 48 frames, Eevee 64
samples, motion blur, 832x480: **183.6 s** wall (3.83 s/frame) at the Change-1
state, i.e. **1.11x proto4's 165.6 s**, against the brief's ~3x failure
threshold. Re-rendered on the final accepted build at the end of the session:
**182.7 s** (3.81 s/frame, **1.10x**), 48 frames rewritten into `renders/`.

**Mid-flight follow check (brief requirement).** Frame 24 of the 48-frame
render (`renders/frame_0024.png`) shows the left moth airborne, and a
hair-on/hair-off diff at that frame (`_dbg/f24-diff.png`) puts the entire pile
at the moth's CURRENT body position with nothing left behind at the launch
site — no weight/vertex-group work was needed, as the brief predicted.

**Side effects measured at the Change-1 state:** D1 sigma_flat 4.963 -> 5.031
(+1.4 %); A1/A4 unchanged (1.434 / 0.345); wing median (111,156,148) ->
(111,157,149); body footprint mean (93.9,121.0,86.6) -> (149.4,158.7,146.3) —
deliberately handled by Change 2.

### S2 Change 2 — body colour + antennae colour (DETAIL-RESEARCH §0-1) — KEPT

The materials are linear-space; the measurement is a lit, Standard-view-
transform render, so the calibration was done empirically
(render -> measure -> adjust), as the brief requires. Five iterations:

| # | `BODY_HEX` | body / shell / hair sheen | rendered body mean (fixed footprint, 2 273 px) | verdict |
|---|---|---|---|---|
| 0 (S1) | `#6f8f66` | 0.85 / 0.60 / 0.55 | (93.9,121.0,86.6) | olive, B -51 short |
| 1 | `#528964` | unchanged | (145.8,157.6,146.0) | R +31 over |
| 2 | `#396e5f` | 0.40 / 0.25 / 0.12 | (106.7,121.4,112.7) | R in, G/B short |
| 3 | `#3b8674` | same | (107.5,128.9,118.5) | G/B still short |
| 4 | `#46bda9` | same | (108.8,148.7,135.4) | B +2.6 short |
| **5 (kept)** | **`#46c3b2`** / tip `#58cfba` | same | **(109.0,151.2,138.7)** | **IN box** |

Iteration 1 is the interesting negative result: changing the base colour alone
moved the rendered mean by only 3.6 levels, because the thorax reads through the
fur shell + pile, whose S1 sheen (0.85 body / 0.60 shell / 0.55 hair) supplied a
large additive white term. Fitting `linear_render = k*albedo + b` on iteration 2
gave b ~= 0.132 linear (a grey wash), which is why the accepted albedo is a
saturated `#46c3b2`: the render adds ~40 % of the wing brightness on top of it.
The sheen was cut in the same step (0.40 / 0.25 / 0.12) so the thorax reads matte
pile instead of a white fuzz (the shell base factor 1.10 -> 1.00 and roughness
0.90 -> 0.92, 0.88 -> 0.95).

**Antennae:** split off the dark leg material into their own `AntennaMat`
(`mat_body`, rough 0.60, sheen 0.25, fur 0.10). Rendered core of the arcs at f0
(`tools/antenna.py`, brightest px in a 3x3 window along the projected tube path,
25 samples per side):

| material | left core (RGB) | right core (RGB) |
|---|---|---|
| S1 `LegMat` `#4a4030` | (66,65,55) | (67,67,60) |
| `#ab9a7e` | (120,111,89) | (125,115,95) |
| `#c4b393` | (132,128,102) | (137,131,106) |
| `#d8c6a4` | (144,137,114) | (148,142,118) |
| **`#e5cab1` (kept)** | **(153,144,122)** | **(157,149,127)** |

Photo antennae median is (157,146,131) — the right arc lands on it. Against the
surrounding render (wing median (111,155,148), local wing pixels ~(73,78,68) at
the head), the arcs now read clearly lighter AND warmer (R > G > B instead of
B > R). `antenna.py` also prints a centre-line median; for the left arc that
value stays dark (70,67,57) because the projected centre line partly misses the
~1 px tube — the core column is the meaningful one. Segmentation of the nodes
(§5) is not modelled: recorded as below-resolution at 832.

**Gates — PASS.** body mean (109.0,151.2,138.7) inside R105-115 / G148-158 /
B138-150; wing median (111,155,148) vs baseline (111,156,148) = -0/+1/0, inside
±3 and the wing materials were not touched; C1 as in the table above.

### S2 Change 3 — wing scale layer (DETAIL-RESEARCH §3, option a) — KEPT

An anisotropic Voronoi cell layer added ON TOP of the existing grain, mottle
and bands (nothing was removed). Nodes: `Mapping(UV) -> Voronoi(F1, 2D,
randomness 1.0)` with the mapping scale set from the target pitch, then the
Voronoi `Color.X` drives a per-cell brightness multiplier, a roughness offset
and (chained after the existing bump) a cell bump; `Distance` is the bump
height. (`voronoi_dimensions = "2D"` came from the S2b merge: cells in the fan
plane rather than a slice through 3D cells, which is the correct reading of a
scale layer on a flat wing.)

**Pitch.** The anatomical rule is wingspan/550. The perched moth's wingspan on
screen is ~343 px at 832 (fan extent 0.535 x root scale 1.157 x 277.3 px/unit),
so the anatomical pitch would be **0.62 px** — below the sampling limit, and
far from the "about 1.5 px" the brief/DETAIL-RESEARCH assume (they assume an
~800 px wingspan). The brief's stated number was followed: **1.5 px along the
span, 0.75 px along the chord (2:1 elongation along the radial direction)**,
i.e. mapping scale (114.5, 154.0) = 171.7 px / 1.5 px and 115.5 px / 0.75 px.
Recorded as a deliberate deviation from the anatomy at this resolution.

**Amplitude steps** (each re-measured). Two workers ran this change in parallel;
the merge (see the "S2b completion note" at the end of this section) accepted
step 4. Rows 1-3 are this session's 3D-cell steps, row 4 the merged 2D setting.

| step | colour jitter | bump strength / distance | rough offset | 2x D1 (block 48) | 2x patch sigma | verdict |
|---|---|---|---|---|---|---|
| 1 | ±2.5 % | 0.16 / 0.00035 | ±0.03 | 4.267 -> 4.454 (+4.4 %) | 5.656 -> 5.851 | not resolving, raised |
| 2 | ±5.5 % | 0.40 / 0.00080 | ±0.05 | 4.267 -> 5.190 (+21.6 %) | 5.656 -> 6.461 | resolves |
| 3 | ±8.5 % | 0.55 / 0.00120 | ±0.07 | 4.267 -> 6.192 (+45.1 %) | 5.656 -> 7.346 | resolves at 2x, but the 5x 1x-crop reads as a faint cellular network -> rejected by the merge |
| **4 (accepted)** | **±2.75 %** | **0.20 / 0.00080** | **±0.025** | **4.267 -> 4.465 (+4.6 %)** | **5.656 -> 5.885 (+4.1 %)** | **kept (Voronoi 2D)** |

**Gate: D1 not below baseline — PASS.** 832x480, segmented mask
(`tools/measure.py`): **4.963 -> 5.782**; geometric outline
(`tools/final_metrics.py`): 4.620 -> 4.962. A small rise, as the brief expected.

**Gate: B2 bands not masked — PASS.** `tools/bandcheck2.py` (same statistic as
`tools/bandcheck.py`, on the S2 images): med 138.3 both; p90/med 1.058 -> 1.058;
p99/med 1.108 -> 1.106; std 7.05 -> 7.06. Band contrast and spread unchanged.

**Gate: 2x comparison shows the layer resolving — PASS, but marginal and
stated as such.** From the same build with the layer bypassed (`--noscale`) vs
enabled, both at 1664x960: `scale2x-before.png` / `scale2x-after.png` (kept in
the folder). Matched *physical* scale (`tools/d1x2.py`, 48 px blocks at 2x = 24 px
at 1x): 4.267 -> **4.465 (+4.6 %)**; patch (546,258)x2 = (1092,516), 48 px:
5.656 -> **5.885 (+4.1 %)**. The S2b note reports the cell rows being visible in
that pair; this session's own read of the amplified 2x diff is that a fine
elongated-cell lattice is present but faint. **The +4.6 % is not the +45 % that
the rejected stronger setting produced** — so the honest statement is: the layer
is in, it is measurably stronger than no layer, and at 1664x960 it is subtle
rather than obvious. Flagged for Session 3 (this is what the 8-second master
will show).

**No tile / reptile reading at 1x — PASS at the accepted setting.** 10x zoom of
a wing patch at 832 (`_dbg/c3d-sbs1x.png` was the stronger setting; the accepted
setting is the smoother of the two): grainier than no layer, no resolved cells,
amplified diff is fine noise. At the rejected stronger setting the 5x crop
(`_dbg/c3d-5x.png`) does show a meandering cellular network, which is why the
merge halved it — this session agrees with that call on inspection.

### S2 Change 4 — contact shadow root fix (DETAIL-RESEARCH §4) — TRIED, REVERTED

**What was tried.** (i) The contact terms (`0.22*contact_edge + 0.08*contact_fill`
in the code; the S1 docstring's "0.60/0.18" was stale) were removed from the
offset cast decal, and (ii) a new per-moth decal was added at EXACTLY the moth's
rest transform (zero offset, same rotation/scale, +0.8 mm above the cast plane so
the two cannot z-fight) carrying only the silhouette with a tight box blur
(radius 3 at 512x296 = **5.4 px support** at render scale), alpha **0.40**, and
the same opacity keyframe ramp as the existing decals.

**Measured** (E3, geometric mask, same method as parity):

| state | 0-3 px mean | 0-3 px darkest 5 % |
|---|---|---|
| S1 baseline (proto4 contact copy inside the offset decal) | -12.87 % | -59.26 % |
| contact terms removed, no new layer (`--nocontact` equivalent) | -11.45 % | -59.22 % |
| new zero-offset layer, alpha 0.40 | **-11.82 %** | -59.22 % |
| new zero-offset layer, alpha **1.00** (mechanism probe, not a candidate) | -12.50 % | -60.00 % |

So the layer contributes **0.37 points** at the brief's starting alpha and only
**1.05 points** even at alpha 1.0. The gate wants <= -15 %, i.e. ~3.2 points
more than the no-contact state; even an infinite-alpha version of this
construction cannot get there (extrapolating the two points gives ~alpha 3.5).
Halo guard was never in danger (darkest 5 % stays at -59.2 to -60.0, well inside
the -68 % cap), and the 1x view shows no drawn-outline look (`_dbg/c4-diff.png`:
only a thin 1 px halo hugging the wing edge plus a speckled showing-through over
the pile).

**Why it cannot work in this scene (measured, not assumed).** The zero-offset
decal sits UNDER the moth, so it can only show where the moth does not cover it.
`tools/edgeoffset.py` on this build: the rendered moth covers **100 %** of the
geometric fan outline (0 px of the outline falls off the moth) and **only 21.9 %**
of the 0-3 px band is bare bark; `tools/contactprobe.py` (decals differing only
by this layer) finds just **380 of the 3 042** band pixels changed. The layer's
own alpha field (`tools/contactfield.py`, decal-only render vs plate-only) peaks
at 0.22-0.28 well *inside* the outline and is only 0.09 within 0-2 px outside it
— exactly the part the wing's own antialiased footprint covers. The E3 0-3 px
band is therefore darkenable by an *offset* soft shadow (whose gradient lands in
the band) but not by a coincident contact layer.

**Verdict: reverted** to the proto4 stack, per the brief ("if it still fails,
revert ... and log"). Verification that the revert is exact: `make_shadow_image`
is byte-identical to the S1 baseline except (a) its stale alpha docstring was
corrected to the values the code actually uses and (b) an unreachable duplicate
block after `return img` was deleted (`_dbg/cmp_shadow.py`); the zero-offset
decal, its texture function, its opacity keys and the `--nocontact` flag are
gone. The experiment is reproducible from this log: `tools/contactprobe.py`,
`tools/contactfield.py`, `tools/edgeoffset.py` and `tools/e3band.py` are kept.

**Honest consequence:** the S2 final 0-3 px band is **-11.70 %**, not the S1
-12.87 %. That difference is NOT the shadow stack (which is byte-identical) — it
is Change 2. Two independent readings agree on the mechanism:
`tools/e3band.py` shows the band's own luminance rising on both sides
(n 3 042 both; lum_ns 103.81 -> 107.57, lum_f0 90.45 -> 94.93), and the S2b
worker's `--nohair` A/B keeps the *absolute* drop identical (-13.37 levels vs
-13.36 at S1) while the relative figure moves because the `--noshadow` reference
the drop is divided by got brighter. So the E3 mean is not a clean shadow metric
across the S2 colour change and should be re-baselined before Session 3 uses it.

### S2 end-of-session parity table (final accepted state)

| metric | S1 baseline | S2 final | delta | tool |
|---|---|---|---|---|
| A1 forewing L/W | 1.434 | **1.434** | 0 | `measure.py` |
| A2 fan IoU | 0.777 | **0.777** | 0 | `measure.py` |
| A2 render IoU | 0.780 | **0.779** | -0.001 | `measure.py` |
| A3 designed/peaks/clear | 7/7/5 | **7/7/5** | 0 | `measure.py` |
| A4 body/wingspan | 0.345 | **0.345** | 0 | `measure.py` |
| D1 grain, segmented | 4.963 | **5.782** | +16.5 % (Change 3, merged setting) | `measure.py` |
| D1 grain, geometric | 4.620 | **4.962** | +7.4 % (S2b's own re-render: 4.992) | `final_metrics.py` |
| E3 geom 0-3 px | -12.87 % | **-11.70 %** | +1.17 pt (body brightness, not shadow) | `final_metrics.py` |
| E3 geom 3-6 px | -23.48 % | **-21.85 %** | +1.63 pt | `final_metrics.py` |
| E3 geom 6-12 px | -24.23 % | **-23.04 %** | +1.19 pt | `final_metrics.py` |
| E3 geom 12-24 px | -14.59 % | **-14.20 %** | +0.39 pt | `final_metrics.py` |
| E3 geom 24-60 px | -1.19 % | **-1.18 %** | +0.01 pt | `final_metrics.py` |
| E3 geom 0-3 px darkest 5 % | -59.26 % | **-58.52 %** | +0.74 pt | `final_metrics.py` |
| wing median | (111,156,148) | **(111,156,148)** | exact | `measure.py`, `bodymean.py` |
| body mean (fixed footprint) | (93.9,121.0,86.6) | **(109.2,151.2,138.7)** | Change 2 target box; ±0.2 run-to-run (hair noise) | `bodymean.py` |
| body mean (`thorax.py`, colour-mask method) | (93.9,121.0,86.6) | **not measurable** | see gap 1 | `thorax.py` |
| C1 thorax-footprint overlap | 99.7 % | **vacuous (0 px dark area)** | see Change 1 (b) | `c1parity_s2.py` |
| C1 dark-centre rate (0.55/0.65 x wing median) | 0.000 | **0.000 / 0.000** | passes more strictly | `c1parity_s2.py` |
| antennae, rendered core (R/L) | (67,67,60)/(66,65,55) | **(157,149,127)/(153,144,122)** | Change 2 | `antenna.py` |
| 48-frame render wall time | 165.6 s (proto4) | **182.7 s** on the final accepted build (3.81 s/frame, **1.10x**); earlier states measured 183.6 s (Change 1), 199.3 s (3D cells) and 208.1 s (S2b's independent run). All well inside the ~500 s failure threshold | 1.10-1.26x | `--render` |

### S2 gaps, deviations and environment notes

1. **`thorax.py` can no longer measure the body.** It isolates the body by
   `(g-b) > 14` (olive vs cyan). After Change 2 the body is teal like the wing,
   so it returns a 65 px fringe blob and a meaningless mean
   ((73.8,150.2,134.0)). The body gate is therefore reported with
   `tools/bodymean.py`, which keeps the S1 statistic (mean sRGB over the body
   region) on the S1 footprint frozen from `base-f0.png` (2 273 px). That
   2 273 px footprint reproduces the S1 body mean exactly, so it is the same
   region; only the colour-based *selector* is replaced. Stated, not hidden.
2. **C1 overlap criterion is vacuous** after the pile + colour change (no dark
   region < 0.80 x wing median exists). The dark-centre criterion passes more
   strictly than S1. `tools/c1parity_s1.py` (verbatim port) still runs on the
   final still but reports a 4 px "dark patch" and a 0.0 % overlap, which is
   meaningless at that size.
3. **E3's segmented variant drifted for the same reason** (-3.43 % -> -0.77 %):
   `perched_mask`'s green threshold now includes the brighter thorax/pile, so
   the band composition changed. Use the geometric variant for comparisons, and
   re-baseline E3 before Session 3.
4. **Change 3's pitch deviates from the anatomy** (1.5 px vs 0.62 px) — see the
   Change 3 section; the brief's stated number was followed, and the true
   anatomical pitch at this wingspan is under the sampling limit.
5. **Change 3 amplitude — two workers, one merge.** The brief's only stated
   amplitude rule is to *reduce* once if the cells read as tiles at 1x. This
   session first hit the opposite failure (no resolution at 2x) and raised the
   amplitude twice; the parallel S2b worker, running in a private copy, judged a
   5x crop of the 1x still to show a reptile-like cell network at the strong
   setting and halved it once. On inspection this session agrees with that call,
   so the merged (halved, Voronoi-2D) setting is the accepted one and every
   number in this section refers to it. Cost of the compromise, stated: the 2x
   resolution gain is **+4.6 %** instead of **+45 %** — the layer is present and
   measurable at 1664x960, but subtle. All four measured settings are tabulated
   in the Change 3 stages.
6. **Two writers, one folder — what happened and what it means.** A second
   completion worker ran S2b in `blender-moth-proto5-s2b/` and merged back while
   this session was still measuring: `voronoi_dimensions = "2D"` (adopted — it is
   the right choice for a planar cell pattern and measures slightly better),
   the halved amplitudes (adopted, see gap 5), an earlier `RES_MULT` /
   `--res2x` variant (superseded by this session's `--res2x`, no trace left), and
   the "S2b completion note" at the end of this section. One 2x still was briefly
   written at 3328x1920 during that overlap and was re-rendered. **All numbers in
   this section were re-measured on the final file** (md5 `248545747cd8`,
   snapshots `_dbg/build_moth.h1.py`, `_dbg/build_moth.s2-c3.py`,
   `_dbg/build_moth.s2-revert4.py`, `_dbg/build_moth.s2-accepted.py`), and the
   file's md5 was checked before and after each render pair. If a third writer
   appears, re-verify the md5 before trusting a measurement. Note also that two
   renders running at once inflate wall time: the same final build measured
   182.7 s alone vs 208.1 s in the S2b run that overlapped another render.
7. **Not done (out of scope by the brief):** Session 3 choreography changes and
   the final mp4 assembly. `proto4`/`moth-motion*` were not written to; all work
   stayed in `blender-moth-proto5`.

8. **The hair is not bit-deterministic between runs.** Two renders of the
   identical accepted state differ in 3 662 px (>2 levels, mean |d| 0.105,
   max 56, all inside the pile; bbox x 210-627 y 198-278). The same scene with
   `--nohair` differs in only 59 px (mean |d| 0.0015) — i.e. proto4's known
   Eevee TAA non-determinism, unchanged. Nothing drifts in the *parameters*:
   seeds are fixed, there is no runtime randomness, and `_dbg/motion-log.json`
   is byte-identical across builds, so this is render noise on sub-pixel
   strands, not parameter instability (`tools/detcheck.py`,
   `_dbg/nh-a.png` vs `_dbg/nh-b.png`, `final-f0.png` vs `_dbg/final-check.png`).
   It is ~0.9 % of the frame and invisible at 1x, but it means adjacent frames
   will shimmer very slightly in the pile — flagging it before Session 3
   assembles the mp4.

Logs: `_dbg/s2*.log`, `_dbg/measure-final.json`, `_dbg/measure-s2base.json`,
`_dbg/c1_footprint_s1.npy`, `_dbg/s2-log-draft.md`, `_dbg/probe_hair*.py`.

### S2b completion note — Change 3 amplitude correction + independent verification

A second completion worker ran S2b in a private copy (`blender-moth-proto5-s2b/`) while this
folder was being written by the main session. The two builds converged on the same Change 1-2
state and on the same Change-4 verdict (split tried, gate (a) failed, proto4 stack restored),
and differed only in the Change-3 cell amplitude. This note records the merged outcome so the
numbers elsewhere in this section and the accepted `build_moth.py` agree.

**Change 3 amplitude (merged, accepted).** The layer is kept with the halved setting:
brightness jitter +-2.75 % (`_sc1` 0.055), cell bump Strength 0.20 / Distance 0.0008,
roughness jitter +-0.025, Voronoi `2D`. At the stronger setting (0.17 / 0.55 / 0.14) a 5x crop
of `final-f0.png` reads as a **visible reptile-skin cell pattern at 1x** — the brief treats
that as a gate ("if the cells read as visible tiles at 1x, halve the amplitude once"). At the
reduced setting the same crop reads as natural grain, while the 2x pair still resolves cell
rows. Measured on the merged build with the stills re-rendered (`final-f0/ns/sil.png`):

| measure | stronger setting | merged (halved) | S1 baseline |
|---|---|---|---|
| D1 sigma_flat, segmented mask | 7.118 | **5.782** | 4.963 |
| D1 sigma_flat, geometric outline | — | **4.992** | 4.620 |
| 2x flat-patch high-frequency sigma (`scale2x-before/after.png`) | — | **4.242 -> 4.433 (+4.5 %)** | — |
| 2x flat-patch low-frequency sigma | — | 3.787 -> 3.809 (+0.6 %) | — |
| band proxy p90/med, p99/med | — | 1.058 / 1.106 | 1.058 / 1.108 |

**Change 4, numbers measured independently** (same-build A/B, `--nocontact` hides only the new
decals): 0-3 px band mean -11.46 % (contact off) -> **-11.97 %** (contact on, blur 3 / alpha
0.55 after the one allowed adjustment; blur 2 / alpha 0.40 gave -11.61 %), darkest 5 %
-58.52 -> -59.26 % (guard -68 %), halo band 0-3 px outside the drawn silhouette -0.77 ->
-0.81 % (S1 baseline -1.12 %, i.e. not worse), 3-6/6-12/12-24 unchanged within 0.2 pt. Gate
(a) needs <= -15 %, so the split is 3.3 pp short and the revert stands.
**Why:** 22.1 % of the 3042 px in the 0-3 px ring are covered by the drawn moth (AA + fringe +
pile) and cannot be darkened by any decal underneath; the 77.9 % that is bare ground is already
at the cast decal's 0.90 alpha cap at its inner edge. A zero-offset layer's blur ramp is centred
on the wing edge, so only its outer tail is exposed: the ring samples texture alpha
0.46 / 0.30 / 0.23 at 0-1 / 1-2 / 2-3 px, yet the rendered darkening there is only
-0.5 / -0.55 / -0.12 levels. Reaching -15 % would need the darkening spread 2-4 px further out,
i.e. a visible halo that the same gate's darkest-5 % guard forbids.

**End-of-session parity on the merged build** (stills re-rendered from it, same tools):
A1 1.434; A2 fan 0.777 / render 0.778; A3 7 / 7 / 5; A4 0.345; D1 segmented 5.782, geometric
4.992; E3 geometric 0-3 px **-11.73 %**, 3-6 -21.88, 6-12 -23.02, 12-24 -14.19, 24-60 -1.18,
darkest 5 % **-58.52 %**; wing median **(111,156,148)**; body footprint **(109.0,151.3,138.7)**
(IN); C1 dark-centre 0.000 at 0.55/0.65/0.75 x wing median and 0 px dark patch (overlap
criterion vacuous); antenna cores **R 157,149,127 / L 153,144,122**; 48-frame render wall
**199.3 s** (main session) / **208.1 s** (independent run) — both with the other worker
rendering at the same time, so quote the range: 199-208 s, i.e. +20 to +26 % over proto4's
165.6 s, well inside the ~500 s failure threshold.

**E3 attribution — why the 0-3 px band reads -11.73 % while the shadow stack is byte-identical
to proto4.** With the pile hidden (`--nohair`, same build) the band reads **-12.65 %** and its
*absolute* drop is **-13.37 levels** against -13.36 levels at S1: the shadow itself is unchanged.
The relative figure moves because the 0-3 px ring contains body/pile pixels that Change 1-2
made much lighter (footprint mean 93.9,121.0,86.6 -> 109.0,151.3,138.7), which raises the
`--noshadow` reference the drop is divided by. The remaining 0.2 pp is the pile fringe plus TAA
noise. Not a shadow regression; recorded as a measurement-definition artefact.

**Provenance / artefacts.** Private copy: `blender-moth-proto5-s2b/` (accepted file identical to
this folder's `build_moth.py` apart from comments). Full independent write-up of the same
session (all four changes, measurements, gates) is kept for review at
`_dbg/s2b-independent-log-draft.md`. Logs: `_dbg/s2b-*.log`, `_dbg/measure-s2bmerged.json`,
`_dbg/patch-{before,after}.png`, `_dbg/wing2x-{before,after}.png`.

## S3 — choreography v3 + final render + proto5.mp4 (BRIEF-S3b, Session 3 of 3)

Full report: `stageC-results.md` in this folder. Single writer; the folder was idle
(no Blender/ffmpeg running at start). Session-start `build_moth.py` md5 **248545747cd8**
verified byte-identical to `_dbg/build_moth.s2-accepted.py` (disk size 90 709 B — the
brief's 90 698 B is a typo in the count, the md5 matches). After the edit md5
**8858284bd1d3df4ffb6f5289aacc3fe0** (93 959 B), snapshot `_dbg/build_moth.s3b-v3.py`;
re-checked immediately before the 48-frame render and again after: unchanged.

**Changes — choreography only** (geometry / materials / shadow stack / S2 detail
layers untouched, confirmed by frame-0 identity: S2 `final-f0.png` vs this build's
`renders/frame_0000.png` mean |d| 0.0068, 285 px >2, all in the pile).

1. `BASE_HZ` 9.5 → **8.0** (period 3.0 frames at 24 fps), with the dossier §1
   transfer argument (Morpho 5.9 ± 1.7 Hz / winter geometrids ~4 Hz / f ∝ √(m/A_f)) in
   the comment. Beat starts became **18.0, 21.0, 23.7, 27.06, 29.88**.
2. **f21 right-antenna-only event deleted**; **both-antennae sweep** (28, 0.0)
   (29, 3.0) (30, 0.0) added on BOTH antenna channels of the perched moth, placed
   immediately before the wing events (transferred reflex order, Manjunath et al. 2025;
   amplitude invention).
3. Perched partial wing raise (31, 0) → (37, +9.0) → **(39, 0)**, rest from f39
   (6 f up / 2 f down = the Xanthotype 75:25 open:close asymmetry, amplitude reduced).
4. Everything else byte-identical in value: f5 antenna event, f13 flick, breathing
   keys, f8-15 / f16-17 / f18 liftoff keys, exit path, factors [1.00, 0.90, 1.12, 0.94],
   1-frame left lead.
5. Description string → `proto5 v3 (8.0 Hz transferred; Xanthotype-shaped perched
   events; dossier-2)`. The choreography comment above the table now carries the
   dossier §4 labels (MEASURED / TRANSFERRED / INVENTION) per event.

**Motion verification — PASS** (`_dbg/s3b_motion_check.py`, log `_dbg/s3b-motion-check.log`;
proto4's `tools/motion_verify.py` hard-codes the v2 event list and reads the removed
`R_ant_right_deg`, so its block 1-5 logic was re-implemented on the same log for v3 and
the proto4 tool was left unmodified):
beat starts equal the brief's list to 0.00000; sweep present at f28-30 on both
channels; **no antenna event at f21** (max 0.0° over f19-f27); raise peak 9.0° at f37
with 0 at f31 and f39 and rest from f39; f5 event and f13 flick unchanged; left moth
rest 3.0 / f6 dip 1.5 / raise 28.0 @f15 / hold f16-17 / liftoff f18; root fixed f0-f17
and the exit path bit-identical to v2 at f24/f28/f36/f47; 1-frame lead preserved.
Per-object check on `proto2.blend` (`_dbg/s3b_ant_probe.py`): `Moth_R_antenna_L` and
`_R` are identical at every sampled frame and both read +2.50° at f5, +3.00° at f29,
0.00° at f19/21/23/28/30/39. Log diff vs v2 (`_dbg/motion-log.v2-baseline.json`):
**f0-f18 identical**; differences confined to the left wing/twist series from f19
(new beat grid) and the perched antennae/wing at f21-23, f29, f31-38.

**Stage C (proto4 → v3), frame-0 stills** (`s3b-f0/ns/sil.png`, same tools):
A1 1.434 exact; A2 fan **0.777** / render **0.779** (0.778-0.780, hair noise); A3 7/7/5;
A4 0.345; D1 segmented 5.782 (S2 layer) / geometric **4.948** (S2 measured 4.962-4.992);
E3 geometric 0-3 px **-11.76 %**, 3-6 -21.80, 6-12 -23.02, 12-24 -14.20, 24-60 -1.18,
darkest 5 % **-58.52 %**; wing median **(111,156,148)** exact; body footprint
**(108.9,151.2,138.7)** IN box; antenna cores **R (157,149,127) / L (153,144,122)**;
C1 dark-centre **0.000** at 0.55/0.65/0.75 × wing median and **0 px** dark patch
(overlap criterion vacuous). Everything reproduces the merged-build baseline within
TAA noise; nothing moved materially.

**Exit check** (same silhouette detector on both builds): proto4 f24 25 883 px
(11.4 % "outside"), v3 f24 **27 976 px = 4.8 % outside**; f47 = **0 px** in both, and
also 0 px from f36 on. The f24 difference was investigated: root position and
pitch/roll/yaw at f24/f28/f36/f47 are **bit-identical** to v2 and no bbox is off-frame,
so this ratio is pose-sensitive — at f24 the new beat grid puts the wings at +27.1°
(top of stroke) where v2 had them at -28.4° (bottom). Path unchanged; not a slower exit.

**Render + assembly.** 48 frames, Eevee 64 samples, motion blur, 832x480, wall
**182.9 s** (3.81 s/frame, 1.10× proto4's 165.6 s; S2's accepted build measured 182.7 s
alone). `assemble.py` was re-pointed at `proto5.mp4` (it previously wrote proto4.mp4;
proto4.mp4 untouched) and committed after verification:

```
ffprobe: {"codec_name":"h264","width":832,"height":480,"r_frame_rate":"24/1",
          "avg_frame_rate":"24/1","duration":"2.000000","nb_read_frames":"48",
          "pix_fmt":"yuv420p","bit_rate":"995664"}  format size 250366 B (244.5 KB)
```

832x480 / 24 fps / 48 frames / 2.000 s all verified. `f00.png`, `f12.png`, `f24.png`,
`f47.png` copied from `renders/`.

**Remains** (detail in `stageC-results.md` §6): A2 still 0.777/0.779 vs the 0.90 gate
(structural cap — the photo's own mirror is 0.901); the contact layer stays reverted
(zero-offset decals cannot darken the 0-3 px band; 0-3 px now -11.76 %, ~1.1 pt of which
is the brighter body raising the no-shadow reference); the hair pile is not
bit-deterministic (same build, two renders: mean |d| 0.107, 3 778 px >2, parameters
byte-identical); perched antenna behaviour is still unmeasured in footage (both events
are TRANSFERRED order / invented amplitude); 8.0 Hz is a band transfer, not a
measurement of this species; `_dbg/build_params.json` still says `"proto": 4`
(unrelated lineage field, left untouched by design).

Logs: `_dbg/s3b-*.log`, `_dbg/s3b_*.py`, `_dbg/measure-s3b.json`,
`_dbg/motion-log.json` (v3) + `_dbg/motion-log.v2-baseline.json`.
