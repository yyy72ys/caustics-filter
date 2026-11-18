# Rally v3.7 — GPT digest (2026-11-11, thread 6af2ea4a, D1–D5 round)

Overall: "v3.7 not yet validated"; v3.4 REJECTED as a physically valid
environmental wave generator (transport machinery stays valid) — the old
gate could report success while the state was wrong.

Verdicts (final matrix): D1-replace SURVIVES / D2 WEAKENED /
D3-retire-taper SURVIVES / D3-adopt-0.034 WEAKENED / D4-patch SURVIVES /
D4-gates WEAKENED / D4-cascades SURVIVES / D5-dispersion SURVIVES /
D5-quantisation WEAKENED.

Load-bearing points:
- D1: fetch law needs DURATION labelling (t_min ≈ 7 min @F100/U6,
  34 min @F1000/U6 — label states fetch-equilibrium when duration
  unknown); freeze the wind convention (U10 everywhere; u* from U10).
  Depth correction is REGIME-gated: kh_p at h=0.3 is 11.0 (U1.5/F100)
  down to 0.67 (U10/F1000) — dispersion vs depth-limited GROWTH are
  different corrections; use Young–Verhagen regime, not blanket.
- D2: conventions are the trap — JONSWAP is an elevation spectrum,
  ECKV B is a curvature spectrum; document units/Jacobian/2π factors.
  Never renormalise B_h to a target MSS after assembly (that re-creates
  the old calibration under a new name). Phi(k,θ) must integrate to 1
  and the slope COVARIANCE (not just total MSS) is a failure surface.
  Freshwater transfer of ECKV = model prediction, not truth.
- D3: SPLIT the claim — taper retirement is decisive; swapping 0.0548
  for 0.034 is NOT settled (Donelan wave-age evidence is about
  aerodynamic roughness z0, not optical MSS — the three quantities must
  not be conflated). Keep both as competing reference hypotheses;
  MSS_total = MSS_gravity + MSS_capillary reported per band.
- D4: lambda_p <= S/4 is marginal (4 cycles) — add a patch-doubling
  CONVERGENCE test (peak, MSS, caustics drift). Cascades need a
  partition of unity + shared time/dispersion + phase consistency.
  Nyquist is not derivative accuracy — roll off before Nyquist.
- D5: quantisation error up to ±0.39 rad/s hits LONG waves hardest
  (39–79% at omega 0.5–1) — measure per-band phase error against an
  explicit tolerance; keep physical omega as source of truth.
- Compatibility invariants: Beer–Lambert unchanged for identical
  optical inputs (wave-depth and optical path are separate
  dependencies — never double-count); frozen composite bit-stable for
  identical buffers; eta=1 route golden test I_v3.7(eta=1)=I_frozen.
- New: N1 fetch is DIRECTIONAL F(theta_wind, shoreline) — a directional
  effective-fetch input suffices for v3. N2 LOW-WIND POSITIVITY:
  ECKV alpha_m = 0.01[1+ln(u*/c_m)] goes NEGATIVE at U=1.5
  (u*≈0.057 < c_m≈0.23) — positivity test across all cells is P0;
  never silently clamp. N3 cutoff needs a sensitivity test
  (800/1000/1200/1500) with a physical damping envelope.
- Completion conditions: spectral correctness + numerical convergence
  + compatibility golden tests; an 8-cell pass alone is not enough.
