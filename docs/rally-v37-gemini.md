# Rally v3.7 — Gemini (AI Studio) digest (2026-11-11, prompt 1qoJkIHSoA67, D1–D5 round)

Verdicts: D1 WEAKENED / D2 DEAD (without continuity proof) /
D3 DEAD / D4 SURVIVES / D5 SURVIVES + one inter-module failure claim.

Load-bearing points:
- D1: depth cap H_s,max ≈ 0.55h is 5–18x the fetch value at the pond
  cells — a depth-limited H_s correction CANNOT bind here ("physics
  theatre"). Partly right for the pure depth CAP; the regime-gated
  Young–Verhagen growth correction (GPT) is the defensible middle.
  Also: the H_s override is an "escape hatch" → design answer: override
  stays (it is an environment statement) but is LOGGED as
  hs_source=override in every numbers file, never silent.
  Focal-plane point (true H_s will weaken z=1 m caustics) = physics,
  not a defect; noted for expectation-setting on the windy re-render.
- D2: same splice attack as Grok (step/kink → Gibbs rings in the
  caustic web; Phi(k) extrapolation; double counting at k_m).
  → Answered by the additive ECKV construction (no splice).
- D3 DEAD: sharpest single point of the round — ECKV's OWN B_h at
  U=6 (u*≈0.22) integrated to ~1200 rad/m gives MSS ABOVE 0.034, so
  "accept ECKV ≈ Cox–Munk 0.034" is internally inconsistent; and the
  taper was the only closed-loop sigma^2 check — removing it without a
  replacement acceptance test leaves the pipeline open-loop.
  → Design answer: no MSS target at all (MSS is an output, reported
  against CM / old-target / full-ECKV references) AND a replacement
  acceptance-test battery (per-band shares, H_s-vs-law, realisation
  vs analytic integral, positivity, no-silent-noop).
- D4: two-cascade is MANDATORY, not an escape clause (at U=10,
  lambda_p 2.8 m vs a 1.5 m scene patch); cascades must share ONE
  quantised clock or cross-modulation beats appear.
- Inter-module claim: capillary Hessian variance (∫k⁴S dk) explodes
  1e4–1e6x vs the truncated baseline → det J zero-count multiplies →
  Q=25 sun-disc quadrature shredded → eta=1 route collapses into
  salt-and-pepper. PARTIALLY physical (real curvature is large), but
  the transport sees the FOOTPRINT-filtered field; design answer:
  photon cross-check + C-percentile gates after v3.7, with the
  research's sub-footprint micro-normal route armed as the remedy
  if (and only if) those gates fail.
