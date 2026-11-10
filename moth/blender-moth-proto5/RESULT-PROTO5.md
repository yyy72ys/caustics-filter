# RESULT — Blender moth proto5

Base: proto4 (`delivery/blender-moth-proto4/`). Yardstick: `MOTH-DOSSIER.md`.
User directives for this round (2026-11-10 22:44): (1) MOTION — learn from other
closely related, data-rich species and transfer; (2) DETAIL — research how the
fine modelling can be pushed further, then apply it.

Build: `build_moth.py` (this folder), final md5 `8858284bd1d3df4ffb6f5289aacc3fe0`
(S2 detail state md5 `248545747cd8` + choreography v3 only). Render: Eevee
64 samples, motion blur, 832×480, 48 f @ 24 fps, wall **182.9 s** (1.10× proto4).
`proto5.mp4` ffprobe-verified by the session and independently by the
coordinator: h264, 832×480, 24 fps, **48 frames, 2.000 s**, 250,366 B.
Full records: `STAGE-LOG.md` (S1/S2/S3 + merged note), `stageC-results.md`,
`MOTION-DOSSIER-2.md`, `DETAIL-RESEARCH.md`, footage data on PC in `moth-motion2/`.

## Part A — Motion research v2 → choreography v3

Footage: 9 clips obtained (resolution ceiling re-confirmed by an 11-client
matrix: 640×360 max, one 480×360). Literature cross-check in MOTION-DOSSIER-2.

What changed in the animation (everything else byte-identical in value —
motion log f0–f18 identical to proto4's v2):

| item | proto4 | proto5 | basis |
|---|---|---|---|
| wingbeat BASE_HZ | 9.5 | **8.0** (period 3.0 f; beats 18.0/21.0/23.7/27.06/29.88) | TRANSFERRED |
| perched antenna event | f21, right antenna only | **deleted** → both-antennae sweep (28,0)(29,3°)(30,0), immediately before the wing event | order TRANSFERRED |
| perched partial wing raise | (30,0)(34,9)(36,9)(39,0) with hold | **(31,0)(37,+9°)(39,0)**, rest from f39 — 6 f up / 2 f down | shape MEASURED (class) |
| takeoff / exit / lead | — | unchanged (verified bit-identical root path) | MEASURED (proto4 + S1) |

Motion verification: ALL CHECKS PASS (beat starts exact to 0.00000; per-object
probe: both right-moth antennae +2.50° @f5, +3.00° @f29, 0.00° across f19–f27;
raise peak 9.0° @f37; right root never translates; f47 left moth = 0 px).

### Motion checklist v2 (per-item source labels)

- Takeoff structure (wings first → grounded raise → short hold + pitch →
  downstroke liftoff; legs do not push first): **MEASURED** (geometrid, proto4;
  re-confirmed in S1 across *Xanthotype*, *Paonias*, *Vanessa*).
- Grounded raise ≈62–79 % of the pre-liftoff window: **MEASURED** (3 species).
- Wingbeat 8.0 Hz: **TRANSFERRED** — Morpho 5.9 ± 1.7 Hz anchor, winter
  geometrids ≈4 Hz lower bound, scaling f ∝ √(m/A_f); NOT measured for this
  species. (Manduca's 25 Hz / 115–120° amplitude explicitly NOT transferred:
  different size, loading and muscle type.)
- Irregular beat factors [1.00, 0.90, 1.12, 0.94]: kept, seeded (proto4).
- 1-frame left lead at the first downstroke: **INVENTION** (both geometrid
  datasets show simultaneity within sampling resolution).
- Both-antennae sweep before the wing event: order **TRANSFERRED** (hawkmoth
  rapid-initiation reflex sequence, Manjunath et al. 2025), amplitude
  **INVENTION** (antennae are not resolvable in any obtained footage —
  3 datasets in a row).
- Partial wing raise 6 f / 2 f shape: **MEASURED class+shape** (*Xanthotype*
  perched open/close, open:close ≈ 75:25, S1), amplitude **TRANSFERRED**-reduced
  (full cycle ×2.53 silhouette area NOT copied; +9° partial kept).
- f5 antenna event: **INVENTION** (kept from proto4, labelled).
- f13 wing flick: **TRANSFERRED class** (peacock startle flick; *Hypagyrtis*
  small raises).
- Breathing keys: **INVENTION** (abdomen pumping never established).
- Pre-flight shiver: **researched and REJECTED** — geometrids are
  thermoconformers with very low wing loading (Heinrich & Mommsen 1985), so
  hawkmoth warm-up does not transfer; also conflicts with the user's
  don't-delay-takeoff rule.
- Butterfly basking posture: **NOT transferred** (dorsal basking ≠ geometrid
  flat rest); Morpho cycles used only as irregularity statistics (periods
  18.4/24.8 s, CV 0.148 — not a metronome).
- No reaction of the staying moth to the other's takeoff: no data either way;
  implemented as none.

## Part B — Detail research → build (DETAIL-RESEARCH.md)

Photo ground truth first (coordinator-measured): the photo's thorax is teal
pile ≈ wing lightness (median (109,152,145)); the "brownish" elements are the
pale tan antennae (median (157,146,131)). The mismatch decomposed into body
hue/lightness + antenna colour + hair length, and was fixed in that form.

| change | verdict | key numbers |
|---|---|---|
| 1. Thorax particle hair (7,000 + 3,500 strands, children off — a measured Blender 4.2 `factor_random` trap) | **KEPT** | fringe 202 px standing proud of skin (mean 1.27 px, p90 3 px); C1 dark patch 2,020 → **0 px**; A2 0.779; pile follows the moth in flight (f24 diff) |
| 2. Body colour + antennae (BODY_HEX `#46c3b2` saturated teal compensating a measured ~40 % additive wash; antennae own material `#e5cab1`) | **KEPT** | body footprint **(109.0,151.2,138.7)** inside target box R105–115/G148–158/B138–150 (proto4: (93.9,121.0,86.6)); antenna cores **R (157,149,127) / L (153,144,122)** vs photo (157,146,131); wing median untouched (111,156,148) |
| 3. Wing scale layer (anisotropic Voronoi cells, 1.5 px pitch, 2:1 along span) | **KEPT at halved amplitude** | strong setting (jitter ±8.5 %) hit the brief's 1× tile gate on independent check (visible reptile-skin in a 5× crop) → halved per the rule (jitter ±2.75 %, `_sc1` 0.055): D1 segmented 4.963 → **5.782**, geometric → **4.948**; 2× patch HF σ +4.5 % (layer present and resolving at 1664×960, subtle); B2 band statistics unchanged. Pitch deviation recorded: anatomical pitch at this on-screen span (343 px) is 0.62 px, sub-sampling; the brief's 1.5 px was followed |
| 4. Contact-shadow root fix (zero-offset contact decal; contact terms removed from the offset decal) | **TRIED, REVERTED with measured proof** | zero-offset layer contributes only 0.37 pt at α0.40 and 1.05 pt even at α1.00 vs the −15 % gate: only 21.9 % of the 0–3 px band is bare bark (the moth covers the rest, and the layer sits under it). Both workers reached the same verdict independently. Proto4 shadow stack restored byte-identical |

## A–F re-scores (MOTH-DOSSIER checklist)

- **A shape**: A1 L/W **1.434 PASS** (1.4–1.5); A2 fan 0.777 / render 0.779 —
  **gate 0.90 still not met, structural** (the photo's own mirror IoU is 0.901;
  a symmetric asset cannot pass without warping to one photo's pose, judged
  wrong in proto3 and held); A3 scallops 7/7/5 PASS; A4 0.345 (proto4's
  explicit definition; unchanged).
- **B bands**: maintained under the scale layer — p90/med 1.058→1.060,
  p99/med 1.108→1.111; band contrast unchanged, only micro-grain rose. PASS.
- **C body**: C1 dark-centre **0.000** at 0.55/0.65/0.75 × wing median
  (passes more strictly than proto4); dark patch 0 px; the overlap criterion
  is **vacuous** (no dark region exists to overlap) — stated, not claimed.
  C3 fuzz: NEW capability — see Change 1; thorax colour = photo pile family,
  antennae pale tan and segmented-look at the head (annulation itself is
  sub-resolution at 832 and not modelled).
- **D surface**: D1 segmented **5.782** / geometric **4.948** (proto4
  4.963/4.620); the scale layer resolves at 2× and is invisible as tiles at 1×
  at the accepted amplitude.
- **E rest/contact**: E3 geometric 0–3 px **−11.76 %**, darkest-5 %
  **−58.52 %** (halo guard fine). The shadow stack is byte-identical to
  proto4; the apparent −1.1 pt vs proto4's −12.87 % is proven to be the
  brighter body inside the band raising the no-shadow reference (absolute
  drop identical: −13.37 vs −13.36 levels with the pile hidden). E3 must be
  re-baselined on the proto5 values from here on.
- **F motion**: checklist v2 above; all machine checks PASS.

## Corrections recorded this round

- proto4's "C1 dark-centre rate 0.002" was **not reproducible** — S1 measures
  0.000 at thresholds 0.55–0.70 in both centre boxes (full sweep in STAGE-LOG).
- C1 overlap 99.7 % (S1 script) vs 100.0 % (S2 frozen-footprint variant) is a
  definitional difference between the two audit scripts; both ≥ the gate's
  intent, and the criterion is vacuous in proto5 anyway.
- The f24 "outside fraction" (4.8 % vs proto4's 11.4 %) is a pose-sensitive
  area ratio, not an exit slowdown: root position and pitch/roll/yaw at
  f24/f28/f36/f47 are bit-identical to proto4; at f24 the 8 Hz grid puts the
  wings at the top of the stroke (+27.1°) where 9.5 Hz had them at −28.4°.
  Full exit is 0 px from f36 on in both builds.
- Change 3's accepted amplitude is the merged (halved) one; the main S2
  section's strong-setting numbers are preserved in STAGE-LOG next to the
  merged note that supersedes them.

## Honest remains

1. A2 silhouette IoU 0.78 vs the 0.90 gate — structural cap, unchanged since
   proto3; closing it needs a pose-asymmetric asset or a different gate.
2. The contact/AO root fix does not exist in this decal construction
   (measured). A future route would need the darkening spread *outside* the
   outline, which the halo guard forbids, or real geometry/AO in the 8s
   piece's own scene.
3. The hair pile is not bit-deterministic between renders (Eevee TAA on
   sub-pixel strands: same build twice → mean |d| 0.107, ~0.9 % of pixels;
   parameters and motion log byte-identical). Expect a very slight shimmer
   in the pile across frames.
4. Perched antenna behaviour remains unmeasured in footage (3 datasets);
   both antenna events are order-transferred / amplitude-invented.
5. 8.0 Hz is a literature band transfer, not a measurement of *G.
   papilionaria*; real-time wingbeat frequency was never established.
6. Visual (coordinator's eyeball of f00): the wings still read flatter and
   more uniformly mint than the photo's mottled, bark-blended wings, and the
   thorax rim is still a fairly clean disc boundary — the pile fixes the
   interior and the silhouette edge more than the rim. The moths still sit
   slightly *on* the bark rather than *in* it.
7. `build_params.json` still carries `"proto": 4` (lineage field, left
   untouched by the choreography-only rule).

## Process notes (for the record)

Two infrastructure incidents, both recovered without deliverable damage:
the Session-2 monitor was lost in a runtime restart while the PC-side
session lived on, and a parallel completion worker (S2b) ran against the
same build — the two Change-3 amplitudes were reconciled by measurement
against the brief's tile gate, and the accepted file was pinned by md5
(`248545747cd8`) before Session 3 built on it. Session 3's first attempt
died on a transport error after verification only; a completion session
(S3b) finished from the verified state. All PC sessions ran
deepseek-v4.1-flash only, per the project's unification rule.
