# RALLY — v4-1.1 design claims C1–C5 (2026-11-12)
Design: DESIGN-V4-1-1.md. Prompt: adversarial reviewer, per-claim verdicts STANDS/WEAKENED/WRONG + counterargument + failure condition, full text pasted (no attachments).

## Threads (one per AI, continue these)
- GPT: https://chatgpt.com/c/6af588e8-e5bc-83ee-bfab-e1320ae60012
- Grok: https://grok.com/c/cbdd84bd-f592-460b-90bf-e9f864219028
- Gemini (AI Studio): https://aistudio.google.com/prompts/1eUgwu8pGf3mXU4YYsfak33YC4L6bP-8v

## Grok (2026-11-12 18:38, 1 min)
- C1 WEAKENED — free gravity-capillary waves, near-resonance nonlinearity, boils/scars missing; fails when relative submergence <~3–4 or Fr → ~0.7.
- C2 WRONG — claims T(k) is algebraically wrong/incomplete, does not reproduce the long-wave limit, kx=0 zero and resonance handling wrong; fails for finite bed amplitude or Fr ≳ 0.5.
- C3 WEAKENED — Dolcetti is single-size fixed spheres at 4–10 cm depth; insufficient as the anchor for multiscale mobile natural gravel beds; fails outside its envelope.
- C4 WRONG — standing waves need an independent generation/amplitude model (shear-flow–roughness interaction, radial propagation), not just the transfer resonance; fails where surface velocity exceeds the gravity-capillary minimum phase speed.
- C5 WEAKENED — frozen advection OK only for the shortest waves; large structures evolve on the loop timescale; statistics degrade beyond a few seconds / integral length ~0.2–0.3 m.

## Muse's own check on C2 (2026-11-12)
Long-wave limit, flow-aligned case (kx=κ): cosh→1, sinh(κH)/(κH)→1 ⇒ T → [1 − 1/Fr²]⁻¹ = Fr²/(Fr²−1). The limit IS reproduced. kx→0: the term κg/(kx²U²) → ∞ ⇒ T → 0, as claimed. So Grok's stated ground for WRONG ("does not reproduce the long-wave limit") does not hold as stated — possibly aggravated by the ASCII transcription of the formula in the prompt (kappa written as k). To be re-tested in-thread with the derivation shown, after GPT/Gemini land.

## GPT (2026-11-12 18:39, 1 min)
- C1 WEAKENED — decomposition is a fair engineering approximation, but the components are not independent/additive: the stationary band overlaps the bed response spectrally (double-counting risk), and turbulence is not advected at one uniform surface velocity.
- C2 WRONG — on a different ground than Grok: the transfer function as written contains NO surface tension, so its pole is the gravity-only stationary condition (kxU)^2 = g*kappa*tanh(kappa*H); the claim's resonance condition including the capillary term (sigma*k^3/rho) does not match T's own pole. (Muse check: this attack LANDS — the claim conflated T's pole with the capillary-inclusive free-wave dispersion. Fix: state T's pole as gravity-only; surface tension belongs to the free-wave band (B) selection, or T must be re-derived with the capillary dynamic boundary condition.)
- C3 WEAKENED — the flume bed (25.4 mm spheres) is relevant but not equivalent to natural polydisperse gravel; matching sigma_eta and dominant wavenumber alone does not validate the generator for natural beds in general.
- C4 WRONG — a resonance condition sets neither amplitude, bandwidth, nor energy distribution; the undamped linear T is singular at resonance, so a damping / finite-band amplitude model is required as an explicit component (the design's |T| clamp must be promoted to a named, justified rule).
- C5 WEAKENED — advection at surface velocity has experimental support, but Us is not uniform, turbulent structures have finite lifetimes, and an 8 s loop additionally needs the periodicity consistency T_loop = Lx / Us.

## Gemini (AI Studio, Gemini Flash Latest / 3.8 Flash) (2026-11-12 18:39, 2 min)
- C1 WEAKENED — free dispersive gravity-capillary waves (turbulent-pressure forced, propagating at their own group velocity, not just advected) are missing; fails at low submergence H/d84 <~ 2 and irregular banks.
- C2 WRONG — same valid ground as GPT (capillary term absent from T's denominator; T's pole contradicts the claimed capillary-inclusive resonance), plus the pole is unregularized and U is uniform (no vertical shear). Fails in the resonant band, k >~ 100 rad/m, and where boundary-layer shear matters.
- C3 WEAKENED — monodisperse packed spheres give a narrow-band, regular bed signal; natural polydisperse beds (clusters, imbrication, protruding clasts with separation vortices / local jumps / aeration) are not represented. Fails for wide grain sorting and H/d84 < 2.
- C4 WRONG — the undamped linear T is singular on the resonance curve; treating the band as-is gives non-physical infinite spectral energy. Damping or an independent amplitude model is required at every resonant wavenumber.
- C5 WRONG — turbulence surface signatures decorrelate on the eddy turnover time tau ~ H/u* (typically 0.5-2.0 s in shallow flow); 8 s of frozen advection is a "sliding sheet" artefact that ignores several renewals. Fails for any rendering beyond ~1 s.

## Synthesis (Muse, 2026-11-12)
- Tally: C1 WEAKENED 3/3. C2 WRONG 3/3 (GPT+Gemini ground VALID: T's pole is gravity-only, the claim's capillary-inclusive resonance condition contradicts T; Grok's stated ground — long-wave limit — does NOT hold, see the check above; the formula itself stands once the resonance is stated as gravity-only and the pole is regularized). C3 WEAKENED 3/3. C4 WRONG 3/3 (unanimous: an explicit resonance amplitude/damping rule is a required component, not an implementation detail). C5 WEAKENED 2 + WRONG 1 (finite lifetime tau ~ H/u* must be modelled; loop consistency T_loop = Lx/Us required).
- Design changes adopted into DESIGN-V4-1-1 (v2) and BRIEF-V4-1-1:
  1. Resonance stated as the gravity stationary condition (kx*U)^2 = g*kappa*tanh(kappa*H). Surface tension appears ONLY in the free-wave dispersion used by component C.
  2. Component B is a diagnostic partition of A's near-pole band, NOT added energy (kills the double-counting GPT flagged).
  3. Resonance rule (new named component): breaking-limited slope cap on the resonant band, applied deterministically, with the capped energy share printed per scene.
  4. Component C modes carry phase omega = kx*Us + omega0(k) (gravity-capillary intrinsic frequency, Doppler form) and an Ornstein-Uhlenbeck amplitude/phase renewal with correlation time tau = H/u*, tau printed.
  5. Scope labels: the flume gate validates AT flume conditions only; natural scenes carry an extrapolation label (bed type). Validity envelope printed: Fr <= 0.7, H/D84 >= 2.
  6. Loop consistency: Lx and Us chosen so T_loop = Lx/Us holds for the advected part; checked and printed.

## Round 2 (v4-1.2 claims D1–D5, same threads, 2026-11-12 evening)
Design: DESIGN-V4-1-2.md. New evidence briefed: per-condition Table II, k0 is predicted not measured, ring holds 81–94% of variance, friction-regularised complex T, sphere-bed lattice spectrum, turbulence correlations.

### Grok (2026-11-12 19:34)
- D1 WEAKENED — better than the singular real T, but still linear potential flow over a gently varying bed; the 0→π/2→π phase is an idealised 1-D result; in 3-D the resonance broadens and becomes direction-dependent; fails at relative submergence <~4 or with sharp roughness elements.
- D2 WEAKENED — ring+tail fits high-U0 spectra, but the 81–94% share is an OBSERVED property of the stationary mechanism, not an a-priori closure for arbitrary beds; fails near/below c_min or Fr <~ 0.4.
- D3 STANDS — per-condition std within a factor 2, k0 not gating: correct as a minimum consistency check; it is NOT a sufficient validation for caustics statistics or appearance.
- D4 WEAKENED — switching off the lattice-bed response is pragmatically justified, but periodic pressure harmonics can leak to the surface via nonlinear rectification/separation; fails for imperfect packing or depth comparable to sphere diameter.
- D5 WEAKENED — the two closures reproduce orders of magnitude inside their data range; confounded, no demonstrated extrapolation in relative submergence or Reynolds number; fails on natural gravel beds outside the calibration range.

### GPT (2026-11-12 19:34)
- D1 WEAKENED — friction regularisation is well motivated, but "the correct replacement" overclaims: phase agreement does not guarantee amplitude agreement; basal friction is not a universal damping coefficient; the Fr≈1 response remains unverified.
- D2 WRONG — the 81–94% is NOT a measured variance share; it is model-derived (depends on directional distribution and alpha assumptions). Ring and tail are neither independent components nor turbulence-only.
- D3 WEAKENED — per-condition std validation is legitimate as an amplitude gate, but the 0.5–2.0 band is loose, RMS hides spectral error, and 13 conditions from one apparatus/bed prove no generality.
- D4 WEAKENED — scale separation is persuasive for the ideal lattice, but real beds carry low-wavenumber topography and defects; a ks-only representation discards coherent bed information.
- D5 WEAKENED — usable as initial empirical hypotheses, but extrapolation-sensitive (exponential), V undefined, sigma_C confounded with depth, and L0 alone does not fix the spectrum.

### Gemini (2026-11-12 19:34)
- D1 WRONG — basal-friction regularisation presumes shallow water (kH<<1), but the stationary waves live at k0H≈1.2–13 where the bed is exponentially decoupled from the surface; the true regularisation there is turbulent eddy viscosity, radiation damping, nonlinear steepening. Fails for kH≳1 stationary/resonant components.
- D2 WEAKENED — the share is an empirical observation on a narrow flume with uniform spheres, not a physical law; fixing the share while predicting total variance from turbulence scaling contradicts ring amplitude being set by resonant bed-spectrum power. Fails at bed transitions and where bed height variance does not scale with u*^2.
- D3 WEAKENED — for caustics, sigma_eta is nearly irrelevant; curvature is set by the 4th spectral moment m4. Passing sigma within 0.5–2.0 (a 16x energy band) can coexist with order-of-magnitude errors at high k; tolerates slope error at k>100 rad/m.
- D4 WRONG — with the bed spectrum zeroed, the linear transfer gives identically zero surface response; stationary waves cannot arise spontaneously — subharmonic injection to much lower k is needed; detached from the bed, the ring becomes ungrounded phenomenology inside a transfer pipeline.
- D5 WEAKENED — the L0/ks exponential explodes ~650x over V/u* in [3,15], predicting L0 far beyond depth/width, contradicting L0 <~ 2–3H. Fails for V/u* >~ 8.

## Synthesis round 2 (Muse, 2026-11-12) — design v4-1.2 revised again
- Tally: D1 WEAKENED 2 + WRONG 1; D2 WRONG 1 + WEAKENED 2; D3 STANDS 1 + WEAKENED 2; D4 WEAKENED 2 + WRONG 1; D5 WEAKENED 3.
- What lands and changes the design:
  1. (Gemini D1) Basal friction cannot regularise the k0H≳1 band. Regularisation is split by band: friction for kH<<1; turbulent diffusion for the stationary band; slope cap as final guard. In practice the forced part is NOT routed through the singular T at all (see 2).
  2. (GPT D2, Gemini D2/D4) The ring is not closed by a fixed share and not glued on: the surface model for the forced part IS Dolcetti & Krynkin's spectral model (ring at k0 + k^-alpha tail), variance-closed by the Table-II/Nichols scaling of sigma with U0 and H (empirical closures, ranges printed). The direct bed response T*zeta remains as the separate physical term for beds with low-k content (natural gravel); for the lattice bed it is computed and reported (expected ~0), not switched off by fiat.
  3. (Gemini D3) The gate gains spectral content: sigma per condition (measurement) + ring share and tail slope vs the 2017 model (labelled model-vs-model) + cond-1 ring-off behavioural check (U0<c_min) + slope variance reported per condition.
  4. (Gemini D5) L0 is capped at 2.5H; V is fixed as the depth-averaged velocity; closures used inside calibration ranges only, ranges printed; sigma_C sensitivity at +-2x reported.
