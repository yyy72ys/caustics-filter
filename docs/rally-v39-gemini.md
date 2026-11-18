# Rally v3.9 (low-wind floor) — Gemini digest (AI Studio, Gemini 3.8 Flash)
Thread: https://aistudio.google.com/prompts/1qoJkIHSoA67BrMjDktzEsFUE59oKUjmy (continued, 2026-11-11)

Verdicts: R1 SURVIVES. R2 SURVIVES. R3 WEAKENED. R4 SURVIVES. R5 SURVIVES. R6 SURVIVES (direction only above a steepness threshold). Decision (a)+(d): REJECTED — "undisguised goalpost-moving".

Key arguments (Gemini):
- R1: with a caveat — the zero crossing sits suspiciously near the aerodynamic wave-growth threshold (u*_c ≈ 0.04–0.08 m/s); a glassy regime at near-zero wind is physical, but placing the cutoff at U=2.7 m/s is the log fit's extrapolation failure.
- R2: the proof is sound — a share gate cannot be met by driving the integral past the measured total; C1 wrongly assumed capillary dominance is wind-invariant.
- R3: category error — the Cox–Munk 0.003 is open-ocean residual (swell, internal waves, 1954 photo blur), not pond capillary variance. "Free parameter wearing a measurement's clothes."
- R4: survives; Hara's excess came from tidal-estuary shear turbulence a stagnant pond lacks, but α_m=0.00107 sits in the fit's steep-gradient zone, so underestimation at U=3 is probable.
- R5: survives on physics — at U=1.5 there are no sharp crests/micro-breakers; the parasitic source collapses naturally (viscous β ∝ k² wins).
- R6: survives as direction only ABOVE the crest-steepness threshold (a·k ≈ 0.1–0.2, Ebuchi 1993); at U=1.5 on 100 m fetch the dominant wave never reaches it, so the increment there is ~0; direction supported for U ≳ 4 m/s.
- Counter-proposal: (1) replace the floor with a viscous-transition α_m(u*) = α_m0 · ½(1 + tanh((u* − u*_c)/Δu*)), u*_c ≈ 0.04 m/s, honestly labelled as an ambient/thermal noise roll-off, NOT the ocean intercept; (2) replace relative-share C1 with an ABSOLUTE saturation gate: B(k) ≤ max(B_ECKV(k;u*), B̂_Phillips) over 50–200 rad/m, smooth across the 370 junction; calm water is allowed to be gravity-dominated.
- Decisive measurement: lab flume optical slope measurement, short fetch (≤20 m), clean freshwater, U=1.5, k=10–1500 rad/m. σ²(k>370) ≥ 0.002 redeems R3; < 10⁻⁴ kills the floor entirely (pond = glassy surface).
