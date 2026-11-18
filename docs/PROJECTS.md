# PROJECTS — project-level map with full idea contents
# (user's taxonomy, fixed 2026-11-11: v1–v8 are all versions of the
#  filter and live under PJ1; PJ2+ are the separate bodies of work)
#
# VERSION NOTATION (user, 2026-11-12, unified from now on): v4-<stage>.<rev>
#   v4-1.x = stage 1 (surface) — v4-1.1 is its first build (flowing surface)
#   v4-2.x = stage 2 (bottom + objects 3D); v4-3.x = stage 3 (volumetric
#   water); v4-4.x = stage 4 (camera viewpoint only; the eye is PJ9).

## PJ1 — The caustics filter (ACTIVE, this chat)
The artwork: "make Earth's environment in code and express it visually."
Environment nodes (sun / water / wind & waves / weather, physically
coupled; real mode = place + date/time -> astronomical solar position +
live Open-Meteo weather) -> a simulated water state -> light transport
(ray focusing through the simulated surface; caustics must EMERGE, never
painted as procedural pattern layers) -> image. Forms: a web app for
real-time recipe building + an offline physical pipeline for the final
render; one recipe JSON (caustics-recipe) is the shared body. Main
output: seamless 8 s loop video (stills exportable). Current production:
the 8 s moth piece (Blender base + this filter composited over it).
Version line (all under this PJ):
- v1: real-time WebGL web app (FIX-1..16; published on GitHub Pages)
- v2: offline deterministic transport (all-preimage Jacobian, sun-disc
  quadrature; photon mapping kept as cross-check only) + composite v3
  frozen (ambient x2.6 low-pass "Blow", exposure 2.2 on the direct term);
  production 8 s loop delivered 2026-11-10 (user verdict: positive)
- v3.1 water + bottom: wind-readout H_s with Stokes cap, realised slope
  variance (Cox–Munk-calibrated), Smith G1 visibility on the direct
  term, substrate albedo presets, spectral Beer–Lambert depth
- v3.2 rain: rain rate as the only new input; Marshall–Palmer drop-size
  distribution; isotropic ring-wave band whose centre moves with rate;
  rain damping of longer waves (documented approximation)
- v3.3 retina/observer: post-stage AFTER the composite (master
  untouched): display-luminance premise, Moon–Spencer pupil (gated
  model; Watson–Yellott as reference readout), Weber/de Vries–Rose
  adaptation, defocus blur from pupil x dioptre error, CSF analysis
- v3.4 spectrum determination: full dispersion (gravity + capillarity,
  tanh(kh)); short-fetch slope anchor replacing the open-ocean Cox–Munk
  anchor (realised sigma^2 ~2x at pond scale); breaking clamp
- v3.5 sky split: Perez all-weather sky; direct fraction eta splits
  direct vs dome-integrated diffuse transport; measured: diffuse
  caustics persist at eta=0 and contrast is non-monotone near eta~0.25
- v3.6 in-water absorption/scattering split a(lambda)/b(lambda) + HG
  phase function — designed (PIPELINE-V8), implementation pending
- v4: caustics itself into 3D, 3 stages (2D media+surface -> 3D bottom
  and objects with organic simulation -> volumetric 3D water surface);
  viewpoint becomes a variable: human EYE (v3.3 retina observer-side +
  perception/sampling nodes: foveation, saccades/gaze path, peripheral
  desaturation, sample budget) vs CAMERA (geometry + imaging) —
  research filed (ROADMAP-V4), not started
- v5: environmental factors on surface/water/bottom — research filed:
  ADOPT water temperature (strongest), humidity+air-temperature as an
  evaporation pair, surface-film coverage, CDOM spectral absorption;
  recorded as negligible/derived: water pressure, mean air pressure
- v6: meta-map of everything that influences Earth's water — PARKED by
  the user ("a bit later"), do not start
- v7: aquatic organisms reproduced mathematically — living assets built
  from code/math, not mesh-gen AI (Tripo rejected for this line
  2026-11-10: reproducibility in code is the point); the Blender moth
  (proto1–5, frozen as moth v1, 2026-11-11) is the first asset
- v8: cosmic-level numbers (Earth/Sun, solar activity) as the upstream
  layer of v5 — research verdict: NO physical knobs at pond scale;
  log TSI/Kp/flare class as metadata only

## Sub-PJ (under PJ1) — The moth generation line (added 2026-11-11 by the user)
The full production series that generates the moth — PJ1 v7's first
asset line and the 8 s piece's base footage. The user's rules for the
line: the moth is MODELLED in code (no photo cutouts, no mesh-gen AI —
Tripo rejected 2026-11-10 for this moth because reproducibility in code
is the point); camera truly locked; the loop closes exactly by
keyframes; what made the early Wan takes feel alive was the motion's
irregularity, so motion is researched from real footage/related species
and rebuilt as seeded, deterministic irregularity.
Sequence:
- Wan 2.2 tests (2026-11-10): 1 s single-moth (3 takes; user verdict:
  exiting the frame by ~0.5 s is realistic and GOOD — the loop comes
  from the moth returning) and 1 s two-moth tests; route then
  SUPERSEDED by Blender (Wan re-generates every pixel, so the camera
  can never be truly locked)
- Species dossier (2026-11-10): the photo's moth identified as most
  likely Geometra papilionaria (emerald geometer), confidence ~70–80%;
  structural findings (abdomen hidden from above at rest, forewing
  L/W 1.4–1.5, pigment green, scale size) in MOTH-DOSSIER.md
- proto1: mechanics proven (zero drift, exact choreography), read as
  vector-illustration CG
- proto2: one fan mesh per side, photo-sampled colour, silhouette
  contact shadow; abdomen-on-wings identified as the top defect
- proto3: species-researched rebuild; abdomen made structurally
  invisible (centreline seam patch); silhouette IoU 0.74–0.79 vs the
  0.90 gate (the photo's own mirror IoU caps it)
- proto4: takeoff MEASURED from slow-motion geometrid footage (wings
  first, grounded wing-lift ~79%, hold+pitch, liftoff on the first
  downstroke; legs don't push); colour corrected by measurement
  (proto3 was ~23% darker than the photo); MOTION-DOSSIER.md
- proto5 = MOTH v1 (frozen 2026-11-11, pushed to GitHub): wingbeat
  9.5->8.0 Hz transferred from related-species data; thorax fur
  (10,500 hairs), body colour inside the photo's measured box, wing
  scale layer at half amplitude; MOTION-DOSSIER-2.md. Moth v2+ =
  future work ("someday", user's word)
- 8 s piece bases: base v1 (pixel-exact seam; user review: shadow
  wrong, opening too still, not enough flutter, out-and-back is lame)
  -> base v2 (in-frame fluttery circuit, ~1 s opening, large 8 Hz
  wingbeat, height-driven shadow; delivered 2026-11-11)
- USER VERDICT on base v2 (2026-11-11): still far from good — motion
  reads odd, modelling is naive. Standing challenges C-M1 (motion)
  and C-M2 (modelling) + the learn-while-improving fork rule live in
  MOTH-CHALLENGES.md; the challenges are also filed as this chat's
  "current challenges" in the global knowledge cycle, so matching
  knowledge forks bounded improvement rounds (user delegation).

## PJ2 — Water held in shape by a field (NEXT project, user-decided)
The user's own words, 2026-11-08 evening: "in Blender, for a solid —
could you make, in code, an object that is itself water, held by a
core's attraction?" Refined the same night: "not just a sphere — adjust
a magnetic field to hold a solid in a complex shape; but it is water
itself, in a state of swaying with the external environment."
Decided 2026-11-09 in his words: "after the filter is complete, make
this as the next project."
Idea content (from SIDE-IDEA-blender-water-object.md, aligned 2026-11-09):
- Replace one-directional gravity with a designed FIELD: multiple
  attracting/repelling poles + a shape (SDF) field carve the water's
  equilibrium shape. A point core alone settles to a sphere — the field
  design is the artwork's core problem.
- Physical anchor: water is diamagnetic (real magnetic levitation);
  visual cousin: ferrofluid spike formations.
- Shape = BOTH of the user's picks: (1) a clear sculptural form (like a
  statue) AND (2) organic amorphousness — his reason: "controllable,
  and an organic feel." Synthesis = two-layer field: macro layer (SDF
  field) controls the silhouette; meso layer (pole field) grows organic
  protrusions/undulations that move; micro layer: external-environment
  disturbance keeps the surface perpetually swaying, the field's
  restoring force pulls it back to shape. "Silhouette = sculpture,
  contents = living thing."
- Implementation routes noted: self-written SPH particles (pole field +
  pressure + surface tension) -> meshing; mid: particles + metaballs;
  light: Geometry Nodes field deformation + water shader.

## PJ3 — The paper (arose 2026-11-11; prior-art done; go/no-go = user's)
Origin: the user's question "couldn't these calculations become a
paper?" Prior-art research completed (PAPER-PRIOR-ART.md): all three
claims partially covered, cores untrodden.
- Claim 1: diffuse-light caustics PERSIST (at eta=0 the field keeps
  structure: p99 1.23, max 2.61) and contrast vs the direct fraction
  eta is NON-MONOTONE with a bump near eta~0.25, mechanism = Perez
  horizon brightening + Snell-window angular compression. "Diffuse
  caustics exist" itself is published (Stramski 1986; Hieronymi) — the
  bump, the mechanism, and spatial field statistics via deterministic
  transport are unclaimed.
- Claim 2: the pond-scale (fetch ~100 m) short-fetch slope correction's
  quantitative caustic consequence (sigma^2 0.021 vs 0.0107; peaks
  1.45 -> 2.11). Frame as a regime difference, not "Cox–Munk is wrong."
- Claim 3 (cleanest novelty): photon-mapping field-level validation at
  2048^2 is shot-noise-limited at any feasible budget — the deliverable
  image agrees (MS-SSIM 0.996) while the raw field does not (0.285);
  a validation-methodology claim.
Proposed split: one optics paper (claims 1+2) + one graphics methods
paper (method + claim 3). Strongest extra evidence would be a real
tank/pool diffuse-light experiment for the eta~0.25 bump.

## PJ4 — Physical painting (added 2026-11-11 by the user)
His words: "I want to take the output stills and make them into
physical paintings. So: layer-splitting for that — e.g. a layer that
supports the underdrawing relative to the size, then the next layer is
painted like this — a methodology that divides the painting into fine
steps and proceduralises it; and an app that helps mix the colours
with the paints that exist."
Contents:
(a) Layer-decomposition methodology: given a still and a physical
    canvas size, decompose the image into painting layers — a
    support/underdrawing layer scaled to the canvas, then successive
    layers each with a defined painting action (what to paint, in what
    order, covering what) — the craft of painting turned into an
    explicit, repeatable procedure.
(b) Paint-mixing support app: input = a target colour (picked from the
    still) + the inventory of paints he actually owns; output = mixing
    ratios/recipes to reproduce that colour with those paints.
Status: recorded only; design not started.

## PJ5 — Camera app (added 2026-11-11 by the user)
His words: "a camera app — use the waves and caustics we made as a
filter (here it's not about computing; make the computed filter usable,
prepare several presets, let simple nodes be tweaked) and take the
shot with it. It can shoot video too."
Contents: the filter in precomputed/baked, phone-runnable form (the
heavy physics stays in PJ1's pipeline); a set of presets; a small
number of simple node-style controls; a working shutter for stills;
video recording with the filter applied.
Status: recorded only; design not started.

## PJ6 — Music-reactive wave/sound filter work (added 2026-11-11)
His words: "I also want to make a work on the concept of filtering
sound into waves that move organically in reaction to music."
Contents: music drives the wave field — its spectrum/envelope/agitation
becomes the "wind" on the water, so the waves move organically with the
music; and that motion in turn filters/shapes the sound (and image).
Where PJ1's water is driven by environment parameters, here the driver
is music — a sibling artwork on the same body, entered through sound.
Status: recorded only; design not started.

## PJ7 — Tilt-linked viewing (added 2026-11-11 by the user)
Source (his pointer: "in the last two weeks' bookmarks"): @monoradio102,
2026-11-01 — already in GLOBAL-LINKS as **#38** (it arrived with the
moth-refs batch). "Tilt to Peek": from ONE still, generate a 360-degree
orbit video (MiniMax H3 + ORB360 LoRA), then scrub its frames back and
forth in sync with the phone's tilt, so the picture seems to rotate in
3D. Not true 3D — a pre-rendered frame sequence indexed by tilt.
Sibling reference: @L_go_mrk's "oil-motion" (frame-switching driven by
scroll / mouse / drag / touch / device tilt), bookmark box 2026-08-10.
STANDALONE project (user correction, 2026-11-11: not a layer of
PJ1/PJ5 — a PJ in its own right). The work itself: tilt-linked
viewing as a medium — take a work (an image, a loop, a computed
sequence) and make the phone's tilt an index into its frames, so
the piece responds to the hand that holds it.
Application sketches kept as material, not as subordination: with
caustics frames the tilt can scrub COMPUTED states (loop phase,
sun/wind azimuth) so every tilt position stays a physically exact
frame; tilt can also act as gravity/viewpoint input. The project
stands on its own and may draw frames from any other PJ.
User's criterion (2026-11-11): integrating results across PJs is
welcome — what matters is whether the integrated result is
interesting. Standalone bookkeeping, free combination of outputs.
Status: recorded only; source secured; design not started.

## PJ8 — Wave-image reverse engineering (added 2026-11-12 by the user)
His words: 「波の画像からリバースエンジニアリングするプログラムがほしい。リファレンスなげるー、それと同じ静止画を計算的に正しい波でつくる」.
The work: input is a reference image of a water surface (waves); the program works backwards from the image to the environment that would produce it — wind, fetch, depth, water, sun — and then renders the SAME still with computationally correct waves, using the PJ1 forward stack as the renderer. An inverse problem bolted onto the forward pipeline: the forward direction (environment -> image) is PJ1; this project is the return path (image -> environment -> image).
Head start already in PJ1's instruments: the displacement-field yardstick, slope statistics (the Cox-Munk anchor relates glitter statistics to slope variance), caustic-field statistics (C quantiles / p99), and colour-from-medium (a/b) — the inverse direction reads those same statistics off the reference and solves for the scene, then verifies by re-rendering and comparing statistics, not by eyeballing.
STANDALONE project in its own right, drawing its renderer from PJ1 (integration across PJs is welcome per the user's criterion).
Status: recorded 2026-11-12. The user will throw a reference image. Sequencing (user, 2026-11-12): PJ8 comes AFTER v4 — directionally, once v4 is done, this is the direction to take; work starts then, from the first reference.

## PJ9 — The human eye (vision responses) (split out of v4 by the user, 2026-11-12)
His words: 「視覚的部分と視点ね。人間の目の反応は全部まとめて pj にしよう」「v4 ではカメラの角度や機能を搭載するとこまで」.
The work: EVERYTHING about how the human eye responds, bundled as one project — the observer side of the renderer. Contents: the v3.3 retina line (pupil / adaptation / blur / CSF; already implemented as a post-stage inside the v3 stack, and its viewpoint connection moves here), and the perception/sampling node group (foveation and resolution falloff, gaze / saccades, peripheral desaturation and temporal-resolution falloff, predictive fill-in, the total sample budget allocated by attention). The premise recorded 2026-11-11 stands: resolution is a property of the observer's sampling, not of the world.
Boundary (user, 2026-11-12): v4's stage 4 keeps only the CAMERA side — camera angle and camera functions (position / distance / angle plus imaging characteristics). The EYE mode of that stage is this PJ.
STANDALONE project in its own right; it receives the same radiance field the other PJs compute and models the eye that receives it.

## PJ10 — Rain and weather (split out by the user, 2026-11-12)
His words: 「雨も天候として pj 化」.
The work: rain AS WEATHER, as its own project — not a surface knob. Starting content: the re-derivation of rain (the v3.2 rain implementation was rejected as REJECTED_INVALID_DERIVATION in the 2026-11-11 audit rally and its re-derivation was parked as a v4-era item; it moves here). Rain means the falling water itself, its impacts and ring waves on the surface, and the weather frame around it (how rain couples with wind, light, and the surface state the other PJs compute).
STANDALONE project in its own right, drawing its surface from PJ1 like the others.
Status: recorded 2026-11-12; design not started.
