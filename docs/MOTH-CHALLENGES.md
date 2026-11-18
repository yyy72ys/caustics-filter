# Sub-PJ moth line — standing challenges & improvement loop

Opened 2026-11-11 by the user's verdict on the 8 s piece base v2
(video delivered 05:44): **「まだまだだなぁ、動きもおかしいしモデリングも甘い。
今後これらを学びながらも改善してほしい」** — the moth is still far
from good: the motion reads odd and the modelling is naive. Improve
both, learning as we go.

## Challenges
- **C-M1 Motion naturalness.** The flight still reads as animated, not
  alive. Takeoff/landing are part-measured, part-invented; the circuit
  is choreographic invention. Needs: more learning from real footage
  and related species (the proto4/5 method), better seeded
  irregularity, re-check of wingbeat amplitude/body coupling in the
  full piece context (base v2: 8.07 Hz, stroke p2p 1.06x v1 — measured,
  but the user's eye still judges it odd).
- **C-M2 Modelling fidelity.** Wings still flatter / more uniform mint
  than the photo; the scale pattern reads regular/cellular at render
  scale; thorax outline still clean; silhouette A2 = 0.779 (gate 0.90,
  capped by the photo's own pose asymmetry). Needs: mottling/scale
  research (the 9-author moth-refs + DESIGN.md), per-part measured
  gates, one element at a time (proto3 rule).

## Method (standing)
Learn-while-improving: research first (footage, literature, community
practice — research-on-stall applies), then a bounded proto round with
numerical gates, delivered to the user for his eye. Moth v1 stays
frozen on GitHub; improvements accumulate as v2+ candidates, never
silently replacing v1.

## Fork rule (user delegation, 2026-11-11)
This chat takes part in the global knowledge cycle (slug `caustics`).
When the morning pickup (09:00) or the instant-match hook
(`knowledge-drop-match`) brings knowledge that fits C-M1/C-M2, do NOT
stop at a proposal for these two challenges: **fork a bounded
improvement round** (a research task or one proto step, PC via dsh,
deepseek-v4.1-flash), run it, and report the result in this chat.
Anything needing the user's taste judgement still stops at delivery
for his verdict. This delegation covers the moth line only.

## Blind-test gate (added 2026-11-11, from the 9時ラン bringback)
- Method: photo vs render, judged blind by an outside eye (Grok
  vision). Pass = the judge cannot identify the render.
- Round 1 (2026-11-11, delivery/moth-blind-test/VERDICT.md): FAIL —
  Grok identified proto5's render in 26 s. Ranked tells, in order:
  (1) wing colour too uniform (flat mint, #1 tell), (2) pale spots
  too regular/symmetrical, (3) no scale texture at all,
  (4) wing margins too evenly scalloped, (5) thorax a smooth
  featureless oval (no setae/relief), (6) antennae are uniform lines
  (no segmentation/taper), (7) "airbrushed" shading, no crisp
  highlights/occlusion. Render rated render-like on all five
  checked features; photo photo-like on all five.
- This list is the v2 checklist, in this order. Re-run the blind
  test after each proto round.

## Candidate measurement tool for C-M1 — REA-Video (memo, 2026-11-12)

- **What**: REA-Video ("Reverse Engineer Anything – Video") — an OSS
  agent skill that decomposes a local video file shot by shot: cut
  positions/tempo, subtitle timing (SRT/ASS out), telop behaviour, SE
  positions, BGM BPM, audio type; emits report.html plus a
  reproduction plan. Offline tool-only use is free.
- **Why it fits C-M1**: today "irregularity" of the moth's motion is
  judged by eye. REA-Video could extract wingbeat rhythm / motion
  intervals from REAL moth-flight footage as numbers — the
  "reference side" measurement, on the same mechanical-check footing
  as the video-rally-review M1–M5 checks on our own pieces.
- **Status**: user approved recording this memo 2026-11-12 ("おけ").
  NOT yet dispatched. Bounded first step when run: ONE real moth
  video; compare extracted numbers against eyeballing; drop the tool
  if they disagree badly. Extraction accuracy is unverified
  (evidence grade: existence/spec = A, accuracy = C, no third-party
  measurement) — treat any numbers as provisional until that check.
- **Sources （出典）**:
  - X post: @tetumemo, 2026-11-11, status 2120343261948813775
    (quoting Riley Brown @rileybrown, 2026-11-10) — card S3 in
    `~/workspace/youtube-history-day/link-batch-2026-11-11-evening.md`.
  - Repository (independently confirmed 2026-11-11):
    github.com/JuneYaooo/rea-video, Apache-2.0, needs ffmpeg /
    Python 3.10+.
  - Registry: GLOBAL-LINKS resource R22
    (`~/workspace/knowledge/global-links/resource-rea-video.md`);
    LIBRARY.md entry 2026-11-11「リンク束・夜バッチ」.
  - Brought to this chat by the library watch, 2026-11-11 22:39.
