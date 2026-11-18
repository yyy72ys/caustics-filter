# Rally v3.8 — Gemini (AI Studio) digest

Thread: https://aistudio.google.com/prompts/1qoJkIHSoA67BrMjDktzEsFUE59oKUjmy (continued, 2026-11-11)
Claims D1–D5 from the gravity-tail-termination research (research_notes/gravity-tail-termination-spectrum-20261111-0528).

## Verdicts
- **D1 (level problem, not cutoff) — SURVIVES.** Finite-depth objection fails on its own numbers: at k=50, h=0.3, kh=15 → tanh≡1, sech²≈1e-13; the band is deep-water to double precision. Fails only if h<0.05 m.
- **D2 (envelope reset to ECKV's own >370 fraction, ~0.10–0.25) — WEAKENED.** Two attacks: (1) circularity — grading Elfouhaily against Elfouhaily validates the implementation, not nature. (2) Young-sea physics: Laxague's ~5% pure-capillary share was measured at Ω_c≤2; at Ω_c=4–6.3 steep crests (ak_p→0.3–0.4) generate parasitic capillary trains shunting energy from k_p≈5–15 directly into k≈300–800 — a source ECKV lacks entirely (its B_h is wind-stress equilibrium only). The ECKV-referenced target may therefore under-represent a real short-fetch pond.
- **D3 (replace S_g with B_l wholesale) — SURVIVES.** Attacks recorded: γ(Ω_c=6.3)=6.49 vs JONSWAP's usual ≤3.3 (≈5 in extreme tanks) may over-peak k_p; H_s-as-prediction can drift via α_p(Ω_c). Defenses: no seam anywhere (structural smoothness), B_l rolls off before k_m instead of a scale-free power law, α_p(6.3)=0.0169 sits inside measured short-fetch equilibrium α range. NOTE: Gemini's own sketch of B_l's rolloff factor is loose; use the research report's verified formula in the design, not Gemini's.
- **D4 (defer Kudryavtsev) — SURVIVES.** Deferral is the only defensible choice: B_p kernels vary ×4 across the 2003/2005/2013 papers; Bringer shows ECKV vs Kudryavtsev diverge up to 300% in the 100–1000 band. Fails if ECKV B_h totals disagree with young-basin data by >50%.
- **D5 (acceptance battery) — DEAD, on one real point and one misreading.**
  - REAL: the |ΔH_s|>15% per-cell rejection gate is a poison pill. Fetch-law H_s scaling and ECKV's α_p k_p^-2 scaling have incompatible exponents (σ_η² ∝ U^0.55 k_p^-1.725 under ECKV closure); the two cross at one fetch only, so across the 8-cell sweep the ratio swings >25% systematically and the gate would reject valid cells. Gemini proposes an empirical ±35% envelope instead. → Design must demote ΔH_s to report-only or a wide envelope, never a 15% rejection.
  - MISREADING (ours to fix by precision): Gemini attacked "bit-stability" as requiring the new-spectrum render to equal the old. The v3.7 G8a gate is bit-identity of the frozen composite ON IDENTICAL INPUT BUFFERS (code-path regression), not across the spectrum change; renders with the new state are expected to change. The design must state this explicitly.
  - PARTLY FAIR: the ratio-to-ECKV criteria measure implementation fidelity to ECKV given our inputs, not nature-accuracy. Design should label them as fidelity checks and keep nature-facing checks separate (total vs Cox–Munk, B̂ ceiling in 50–200, open band-share reporting).

## Net for the design
D1+D3+D4 stand. D2's target stays but is labelled a fidelity reference, with the parasitic-capillary gap recorded as a known missing source (future Kudryavtsev chunk). D5 rewritten: ΔH_s report-only (wide envelope), G8 bit-identity scoped to identical buffers.
