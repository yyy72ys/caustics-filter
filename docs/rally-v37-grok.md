# Rally v3.7 — Grok digest (2026-11-11, thread c03959aa, round: D1–D5 design under test)

Verdicts: D1 WEAKENED / D2 DEAD (as specified) / D3 WEAKENED / D4 SURVIVES / D5 SURVIVES.

- D1 WEAKENED: replacing waveHeightM is mandatory, but the depth-limited
  correction is second-order at the pond cells (kh ≳ 1–2 at the energy
  peak); an unvalidated depth term risks over-correcting fetch-dominated
  cells. Magnitude must be quantified per cell before it can be trusted.
- D2 DEAD as specified: the splice wavenumber and B(k) continuity are
  unspecified — a jump/kink in the curvature spectrum lands directly in
  the slope PDF and prints artificial rings/gaps into the caustic web.
  Using ECKV only for B_h + spreading with JONSWAP k_p/amplitude is a
  reasonable split, but without an explicit matching condition the
  proposal is incomplete. → Design answer (mine): do not splice; use
  ECKV's ADDITIVE form S(k) = (B_l + B_h)/k^3, re-parameterising only
  B_l's k_p and level from the JONSWAP fetch laws. A sum has no splice
  discontinuity by construction.
- D3 WEAKENED: retiring the taper is correct, but the replacement target
  (ocean Cox–Munk ~0.034 @U6) is under-motivated for ponds — short-fetch
  data (Woźniak/Donelan) do show steeper MSS; lowering the target does
  not by itself restore variance in the right band.
- D4 SURVIVES: gates are the minimal conditions; two-cascade is the
  standard fix; residual risk is implementational (cascade boundary).
- D5 SURVIVES: dispersion already correct; loop quantisation is a frozen
  constraint; the remaining duty is verifying the quantised frequencies
  still sample the new peak + capillary range.
- Frozen-behaviour attack surface: Beer–Lambert and the frozen composite
  are safe if the height field stays stationary; the eta=1 route needs
  the amplitude path to be zero-able (D1's override preserves it).
  Highest silent-breakage risk: the D2 splice and D1's unquantified
  depth term.
