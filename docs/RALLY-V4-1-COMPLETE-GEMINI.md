# RALLY-V4-1-COMPLETE — Gemini leg（verbatim）

- 日時: 2026-11-14 16:28–16:33 JST
- スレッド: https://aistudio.google.com/prompts/1eUgwu8pGf3mXU4YYsfak33YC4L6bP-8v （流れ水面ラリー継続）
- モデル: Gemini 4 Hydrogen（AI Studio、指定どおり選択を確認）
- プロンプト: RALLY-V4-1-COMPLETE-PROMPT.md
- 最終判定: **FALLS**

---

Per-Item Verdicts (L1–L8)

L1. [belongs elsewhere — Stage 3 (3D subsurface convection rolls) & Stage 2 (floating foam/debris particulates)] — Foam and debris lines are discrete surface materials, not a continuous height field (consistent with splash/bubbles in Stage 2). Langmuir circulation rolls are 3D subsurface shear/Stokes-drift instabilities that belong in Stage 3; Stage 1 can only consume their surface convergence velocity/strain field to advect the independent surfactant film scalar field (L2).

L2. [MUST-DO in v4-1] — Local equilibrium demonstrably fails for short-gravity/capillary gusts (τ_relax ∼45 s at λ=0.2 m), and Marangoni damping from surface films acts as a sharp resonant-band dissipation rather than inverse gustiness. Without wave-action balance and an independent film scalar field, moderate- and low-wind pond surfaces cannot reproduce real cat's paws or natural slicks without faking their look.

L3. [MUST-DO in v4-1] — Linear random-phase Elfouhaily waves produce Gaussian symmetric crests/troughs and miss the intense high-curvature ripple trains of parasitic capillaries ahead of steep short-gravity crests. Because caustic topology and fine-line density are governed by surface curvature (∇²η), omitting second-order Sharma-Dean bound coupling and parasitic capillaries leaves the surface linearly distorted and is a prime candidate for the unexplained fine-line excess.

L4. [MUST-DO in v4-1] — Real streams do not have uniform velocity U; bank shear, thalweg acceleration, and bed-induced velocity variations refract surface waves, split dispersion into two branches, and cause wave blocking/steepening where U→−c_g. Assuming uniform U freezes the stationary dispersion ring to one global k_0 across the entire domain.

L5. [MUST-DO in v4-1] — Linear bed-response theory assumes irrotational, non-separated potential flow and misses vortex shedding behind bluff bed features and rocks. Vortex-core dimples (Δη ∼ −Γ²/(8π² g r²)) and upwelling boil mounds are free-surface height/slope deformations belonging to Stage 1; however, to obey the governing principle, they must be physically driven by wake vorticity/separation parameters from Stage 2 obstacle geometry rather than cosmetic ad-hoc stamps.

L6. [not needed — say why] — Surface height/slope perturbations driven by pond thermal convection are physically negligible: thermal expansion (βΔTH ∼10⁻⁵ m) and thermocapillary Marangoni pressure deflections (Δγ/(ρgL) ∼10⁻⁶ m) yield slopes ≲10⁻⁴ rad, orders of magnitude below even minimal capillary ripples. The dominant optical signature of pond thermal cells is refractive index gradients (∂n/∂T) already assigned to Stage 3, while bulk temperature dependence of γ(T) and ν(T) enters Stage 1 solely as environmental scalar parameters.

L7. [MUST-DO in v4-1] — Retaining an uncalibrated operational prior (σ²=0.003) where 4 low-wind cells formally failed directly violates both the claim of completion and the governing principle ("reproduce nature, do not approximate its look"). Low-wind failure reflects real physics missing from standard Elfouhaily: friction-velocity growth thresholds (u_* < u_*,crit), viscous/Marangoni cutoffs, and intermittent gust/film dynamics (L2).

L8. [not needed — say why] — In ponds, streams, and short-fetch scenes, broad frequency bandwidth, wide directional spreading (Δθ ≳45°), and viscous/turbulent damping suppress 4-wave Benjamin-Feir modulational growth long before the modulational timescale τ_BF ∼(ϵ²ω_0)⁻¹ can act (consistent with the verified random-phase groupiness diagnostics). Short gravity-capillary nonlinearity is governed by 3-wave/bound interactions (L3), not BFI.

Missed Items Still Owed by Stage 1

Unresolved "fine-line excess" & lack of fourth-moment (m_4) curvature validation:

Why it matters physically: Sun-disc spread, numerical warp artifacts, and wave groupiness were all cleared, meaning the observed fine-line excess in caustics remains unexplained. Caustic line density and focusing intensity depend on the fourth spectral moment (m_4 = ∫k⁴S(k)dk, curvature variance), yet validation has relied on height variance (σ) and line-depth distributions without verifying high-k slope/curvature spectra against real water surfaces (e.g., Cox-Munk / laser slope gauge data).

Where it belongs: v4-1.

Resolution-independent viscous/capillary cutoff & DK-tail vs. Horoshenkov double-counting check:

Why it matters physically: Dolcetti-Krynkin's model already contains a k^(−α) (α=5–7) turbulent/non-stationary tail alongside the k_0 ring. Adding Horoshenkov's turbulence spectrum risks double-counting high-k energy even if total height variance σ² is budgeted, distorting slope (k²S) and curvature (k⁴S) bands. Furthermore, power-law tails without a physical viscous dissipation cutoff (2νk² ≥ forcing) make ∫k⁴S(k)dk dependent on simulation grid Nyquist frequency k_max, violating neutral numerics.

Where it belongs: v4-1.

Deep-water Elfouhaily applied to shallow water (missing finite-depth shoaling & bottom-friction damping):

Why it matters physically: Caustic scenes are inherently shallow (H ∼0.05–1.0 m), yet standard Elfouhaily assumes deep water (kH ≫1). For low-to-mid wavenumbers where kH ≲2, neglecting tanh(kH) dispersion, TMA-style shallow-water spectral saturation, and bottom-boundary dissipation produces non-physical wave speeds and overestimates low-frequency energy.

Where it belongs: v4-1.

Finite-fetch spatial growth & bank/shoreline wave reflection:

Why it matters physically: Real ponds and streams are bounded, not infinite periodic domains. Wind waves grow from zero at the windward bank over finite fetch x, and waves reflect off banks, shorelines, and emergent rocks (forming cross-wave interference patterns). Periodic Fourier tiling suppresses windward calm zones and creates artificial wraparound coherence.

Where it belongs: v4-1 (consuming bank/shoreline geometry boundaries from Stage 2).

Combined wind + stream current coupling (non-linear superposition):

Why it matters physically: Wind over a flowing stream acts via the relative wind vector U_10 − U_surface, and wind-generated ripples are Doppler-shifted, compressed, and steepened (or blocked) by opposing currents. Simply adding independent wind and stream height fields violates the dispersion relation and wave-action conservation.

Where it belongs: v4-1.

Temporal decorrelation & dispersive propagation of stream turbulence:

Why it matters physically: If the Horoshenkov turbulence field is advected via frozen Taylor translation at U, turbulent boils never deform, decay, or regenerate over eddy turnover times (τ ∼H/u_*), and radiated gravity-capillary ripples fail to propagate at their intrinsic group velocities c_g(k)+U.

Where it belongs: v4-1.

Micro-breaking / steepness stability dissipation cap:

Why it matters physically: Without a physical wave-breaking/micro-breaking energy dissipation limit (capping local downward acceleration a_z ≲−0.5g and limiting wave steepness before overturning), linear/second-order wave superposition and wave-current compression (L4) can produce unphysically sharp cusps and extreme caustic spikes. (Air entrainment/splash remains in Stage 2; the surface-state dissipation cap belongs in Stage 1.)

Where it belongs: v4-1.

Capillary contact-line menisci at emergent Stage 2 boundaries:

Why it matters physically: Surface tension curves the free surface over the capillary length (l_c = √(γ/ρg) ≈2.7 mm) wherever water meets emergent rocks, banks, or vegetation, generating intense localized caustic halos.

Where it belongs: v4-1 (surface boundary condition coupled to Stage 2 emergent solid contact lines).

Final Verdict

FALLS — Stage 1 is incomplete; must-dos before closing v4-1 are: L2 (wave-action gust relaxation + independent Marangoni film field), L3 (Sharma-Dean bound waves + parasitic capillaries), L4 (spatially varying U(x,y) two-branch dispersion + wave-action conservation), L5 (obstacle wake vortex dimples/boils coupled to Stage 2 geometry), L7 (physical closure of the 4 unresolved low-wind cells replacing the σ²=0.003 prior), plus the missed items: resolving the fine-line excess & m_4 curvature validation, DK-tail/Horoshenkov consistency & physical viscous cutoff, shallow-water finite-depth dispersion, finite-fetch/shoreline reflection, wind–stream current coupling, stream turbulence temporal decorrelation, micro-breaking dissipation limit, and contact-line menisci.
