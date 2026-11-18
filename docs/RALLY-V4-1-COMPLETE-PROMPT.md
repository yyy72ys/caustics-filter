# RALLY-V4-1-COMPLETE — prompt（2026-11-14）

三者ラリー（GPT／Grok／Gemini＝AI Studio・Gemini 4 Hydrogen）。主題：v4-1（第1段・水面）は本当に完了か。
スレッド：流れ水面ラリーの各AIスレッドを継続（新規スレッド禁止）。

---

ROLE: You are an adversarial reviewer of a physics-first water-surface simulation project. Your job is to attack the claim below and find anything missed. Do not be agreeable; a missed item found now is worth more than politeness.

CLAIM UNDER TEST: "Stage 1 of v4 (the water surface) is complete. There is nothing left to do in v4-1."

WHAT THE PROJECT IS: A physically based caustics/water-ripple filter. The governing principle (non-negotiable): reproduce nature, do not approximate its look. A numerical device is admissible only if it is neutral numerics or corresponds to real physics. The full pipeline: environment → water state (simulation) → light transport → image. Stage 1 covers ONLY the water surface state (wind waves + flowing-water surface). Light transport/deposit is a separate lane already settled. Stage 2 = bed and 3D objects. Stage 3 = volumetric water (in-water scattering/absorption as 3D fields, density interfaces). Stage 4 = camera/lens. The human eye is a separate project (PJ9); rain is a separate project (PJ10); splash/bubbles were re-homed to stage 2.

WHAT STAGE 1 HAS SETTLED (measured, do not re-derive):
- Wind-wave surface: Elfouhaily-form spectrum is the default; phase function corrected to Petzold FF (v16); low-wind capillary floor uses an operational prior sigma^2 = 0.003, status formally UNRESOLVED (4 low-wind cells never passed, by design decision).
- Flowing surface (stream): phase-locked bed response (linear transfer, anti-phase subcritical) + stationary dispersive waves (Dolcetti-Krynkin ring+tail; full dispersion, not shallow-water) + Horoshenkov turbulence component. Budget accounting fixed so components cannot silently zero each other.
- Deposit/transport: area-conserving mesh projection proven; sun-disc angular spread measured — physical spread does NOT explain fine-line excess.
- Warp-lane audit: all numerical suspects (pixel aspect, resolution, gradient estimators) cleared; warp shifts existing texture band content, creates none.
- Groupiness diagnostic: wave groupiness (GF) of all scenes is inside the same-spectrum random-phase ensemble — groupiness is natural, not a defect. nu (spectral bandwidth) is a height-bandwidth ruler, not a predictor of slope/curvature groupiness.
- Line-depth judgement: single-seed gate replaced by a distribution diagnostic (realization-unit, cluster bootstrap); thresholds 1.05 (mean) / 1.25 (q90) validated on 15 realizations and promoted to operational.

CANDIDATE LEFTOVERS — for EACH, give a verdict: [MUST-DO in v4-1] / [belongs elsewhere — name where] / [not needed — say why]:
L1. Langmuir circulation: wind-parallel convergence streaks (spacing ~ water depth) + a material-line transport lane (foam/film/debris lines — the surface "seams"). Research says material lines are not a height field.
L2. Gust structure ("cat's paws") reworked properly: wave-action balance with band-dependent relaxation (local equilibrium was shown to fail: relaxation ~45 s at lambda = 0.2 m), and surface film as an INDEPENDENT scalar field (Marangoni damping is resonant-band, not the flip side of gusts).
L3. Bound-wave correction: second-order phase coupling (Sharma-Dean QTF) for the gravity band; parasitic capillaries as a separate mechanism.
L4. Stream surface with spatially varying current U(x,y): full dispersion with two branches + wave-action conservation (current implementation assumes uniform U).
L5. Rock/obstacle surface signatures: boils, dimples, scars as intermittent advected events (vortex-core dimple ~ Gamma^2/(8 pi^2 g r^2)); standing patterns over bed features are already in via the bed-response term.
L6. Pond temperature lane: thermal convection cells at the surface + their effect on the surface state (thermal schlieren optics were assigned to stage 3).
L7. The 4 UNRESOLVED low-wind cells: close them, or is "operational prior retained" an acceptable permanent state for stage 1?
L8. Benjamin-Feir / modulational focusing: research concluded out of scope for our scenes (broad bandwidth + wide directional spread + damping dominance) — attack or confirm this exclusion.

THEN: list anything NOT above that stage 1 still owes — physical mechanisms of a real wind/stream/pond surface not yet modelled, validation never done, or consistency checks skipped. For each: why it matters physically, and where it belongs if not v4-1.

FORMAT: per-item verdict lines first (L1-L8), then your "missed items" list, then a final one-line verdict on the claim: STANDS (nothing left) / FALLS (list the must-dos). No vague agreement; if you judge an item "not needed", give the physical reason.
