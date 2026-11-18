# RALLY-V4-1-14 — Gemini 脚

- 脚: Gemini（Google AI Studio、モデル **Gemini 4 Hydrogen** 選択を確認して投稿）
- プロンプトURL: https://aistudio.google.com/prompts/1eUgwu8pGf3mXU4YYsfak33YC4L6bP-8v
- 取得: 2026-11-14 16:53–17:03 JST
- 最終判定: **FALLS**
- 判定内訳: STANDS 1（C13）／WRONG 3（C3, C6, C9）／WEAKENED 14（残り）

---

Claim Verdicts (C1–C18)

C1. Verdict: WEAKENED

Strongest counterargument: Defining U_mat = surface flow velocity + U_S omits the wind-driven viscous/turbulent skin shear current (≈ 0.02–0.035 U_10), which in fetch-limited ponds often equals or exceeds Stokes drift U_S, and omits instantaneous wave-orbital convergence (∇_H ⋅ u_orb) that compresses and stretches surfactant films between crests and troughs (modulating Marangoni damping at the wave phase). Conversely, if surface flow velocity is ever evaluated using instantaneous surface u_orb(z=η), adding U_S double-counts Stokes drift.

Failure condition: A Lagrangian drifter or film-patch simulation in a wind-driven zero-stream pond whose drift speed systematically underpredicts 0.03 U_10 + U_S, or a particle trace advected on u_orb(z=η) + U_S that exhibits 2× Stokes drift.

C2. Verdict: WEAKENED

Strongest counterargument: Standard S_in (Miles/Plant), S_nl (4-wave Hasselmann quadruplets), and S_ds (whitecapping) are oceanic gravity-wave closures that physically break down in the capillary-gravity transition (k ∼ 100–1000 rad/m) where 3-wave interactions dominate and wind-input coefficients carry ±50–100% uncertainty. An imbalance scored by U2 in the gust band conflates genuine disequilibrium of Elfouhaily with the known invalidity of 4-wave S_nl and Plant S_in in that band.

Failure condition: Changing the diagnostic S_nl closure from 4-wave quadruplets to 3-wave gravity-capillary triads flips the sign or magnitude of the diagnosed imbalance in the gust band by >50%.

C3. Verdict: WRONG

Strongest counterargument: Plant's (1982) growth rate β(k) ∝ (u_*/c)^2 ω governs wind-forced growth (u_* > 0), not symmetric relaxation during a gust lull. When wind drops (u_* → 0), τ_Plant = 1/β → ∞, wrongly predicting that ripples never decay! Physical gust decay is governed by viscous damping (2νk^2), Marangoni film damping, and nonlinear cross-band transfer S_nl — mechanisms asymmetric to wind growth and completely missing from a single symmetric Plant τ(k).

Failure condition: A step-down wind test (U_10: 6 → 0.5 m/s) where the relaxation rate predicted by Plant rates misses the actual viscous/film decay timescale by an order of magnitude (e.g., τ_Plant → ∞ while λ = 2 cm waves decay in ∼5 s).

C4. Verdict: WEAKENED

Strongest counterargument: While Marangoni resonance indeed peaks near 3–10 Hz and uniform broadband damping is false, restricting film dissipation exclusively to that resonance band omits two major physical effects of real surfactant films: (1) above resonance (f > 10 Hz, pure capillary regime), an inextensible film shifts the free-surface boundary condition from free-slip to quasi-no-slip, enhancing viscous boundary-layer damping O(νωk) across all short capillary wavenumbers; and (2) damping 3–10 Hz roughness reduces aerodynamic roughness length z_0, suppressing wind stress u_* and wind input S_in to longer bands.

Failure condition: A measured slicked-pond spectrum showing >2× attenuation at f > 15 Hz and reduced wind input at f < 2 Hz that a 3–10 Hz-only dissipation term fails to reproduce.

C5. Verdict: WEAKENED

Strongest counterargument: The claim gets causal direction backwards for L1. While Stage 1 can compute wave/current-induced surface strain, the horizontal convergence field driving Langmuir material seams (∇_H ⋅ u_LC) is the z=0 surface trace of Stage 3's 3D Craik–Leibovich vortex rolls. Stage 1 cannot physically derive Stage 3's convective downwelling convergence from 2D surface height η; Stage 1 must receive that convergence velocity from Stage 3 as an input port, not originate it.

Failure condition: Attempting to generate hierarchical Langmuir convergence streaks (∼3× mixed-layer depth) inside Stage 1 without a Stage 3 u_LC input, which forces either unphysical procedural streak synthesis (forbidden by L1) or zero Langmuir convergence.

C6. Verdict: WRONG

Strongest counterargument: Wilton's 1:2 triad resonance (λ ≈ 2.42 cm, k ≈ 260 rad/m, where c(k) = c(2k)) and parasitic capillary waves are physically distinct mechanisms. Parasitic capillaries are short ripples (λ ∼ 4–8 mm) forced by localized surface-tension pressure gradients at the sharp crests of much longer short-gravity waves (λ ∼ 5–30 cm) matching phase speed c_cap = c_grav (k_cap ≫ k_grav). Parasitic capillaries do not resolve or cover the 1:2 Wilton self-harmonic pole or the continuous resonant-triad manifold spanning the k ∈ [50, 500] rad/m band.

Failure condition: Inspecting the spectrum around k ≈ 260 rad/m in a 3 m/s pond and finding either a spectral hole where QTF was excised without Wilton-triad energy, or mis-phased 4–8 mm parasitic ripples substituting for 2.4 cm Wilton dimpled profiles.

C7. Verdict: WEAKENED

Strongest counterargument: First, the math as stated is wrong: for η_1 = a cos(kx) and η_2 = (1/2) k a^2 cos(2kx), the height ratio is (1/2)ka and the curvature ratio is 2ka (which is 4× the height ratio, not "2ka relative to its height ratio"). Second, for gravity waves bounded by breaking (ka ≲ 0.3), second-order bound curvature adds only (2ka)^2 ≲ 36% to band m_4 — likely too small to explain a large fine-line caustic excess compared to parasitic capillaries (curvature ratio (k_cap/k_grav)^2 ∼ 10^2–10^3) or an unregularized high-k tail. Third, scalar m_4 alone cannot distinguish diffuse curvature variance from phase-locked caustic-forming crest curvature (curvature kurtosis/skewness).

Failure condition: Measuring m_4 before/after η_2 shows <20% change while fine-line caustic excess persists, OR fine-line density changes dramatically under bound phase-coupling at nearly fixed m_4.

C8. Verdict: WEAKENED

Strongest counterargument: Conditioning a fundamental physical mechanism (second-order wave skewness/asymmetry) on whether a naive SVD of the QTF kernel happens to be low-rank violates the governing principle ("reproduce nature"). Moreover, the QTF only acts on the low-k gravity band (k ≪ 260 rad/m), which fits on a coarse sub-grid (N ∼ 64–128) where exact O(N^2) or algebraic FFT-factored Sharma–Dean evaluation is cheap regardless of matrix rank.

Failure condition: An SVD test showing high rank across broadband directional spreads causes U5 to reject bound waves entirely, leaving gravity crests artificially symmetric even though coarse-grid exact QTF evaluation takes milliseconds.

C9. Verdict: WRONG

Strongest counterargument: In natural gravel/rock-bed streams (H ∼ 0.04–0.5 m, λ_0 ∼ 0.05–0.2 m), velocity variations over boulders, bars, and bank shear layers occur on scales L_U ∼ H ∼ d_84 ∼ 0.05–0.3 m. Thus L_U/λ_0 ∼ 0.5–3, directly violating the WKB requirement λ_0/(2πL_U) ≪ 1, and wave blocking/caustics occur around nearly every emergent obstacle. Declaring "out-of-domain" without interpolation or uniform asymptotics (Airy/Pearcey functions or wave-resolving propagation) leaves massive spatial holes in real stream scenes.

Failure condition: Applying the WKB check λ_0/(2πL_U) < 0.1 to any realistic rocky stream reach flags >20–50% of the stream surface around rocks and riffles as out-of-domain.

C10. Verdict: WEAKENED

Strongest counterargument: Evaluating aerodynamic stress and S_in with U_10 − U_surface is required by Galilean invariance of the air-water boundary layer, NOT by wave-action conservation. Wave-action conservation (∂_t N + ∇⋅[(c_g + U)N] = S_tot/ω_r) holds mathematically for any specified forcing term S_tot; using absolute wind produces the wrong forcing rate S_in and wrong wave-age c_r/u_*, not a violation of conservation.

Failure condition: Solving the action-balance equation with absolute-wind S_in and verifying that ∫N dk dx still exactly satisfies the integral action budget ΔN = ∫(S_tot/ω_r) dt, proving the error lies in S_in physics rather than action conservation.

C11. Verdict: WEAKENED

Strongest counterargument: While a band-wise shear criterion is necessary, U6's proposed formula — U_eff(k) using deep-water e^{2kz} weighting — is physically wrong for intermediate/shallow stream bands (kH ≲ 2), which require the finite-depth hyperbolic kernel 2k cosh[2k(z+H)]/sinh(2kH). Furthermore, the Stewart–Joy linear integral ignores profile curvature U′′(z) (Rayleigh equation) and critical-layer resonance where U(z) = c inside the sheared boundary layer.

Failure condition: Comparing the e^{2kz} Doppler shift against the exact Rayleigh dispersion relation on a logarithmic open-channel profile U(z) = (u_*/κ) ln(z/z_0) at kH = 0.5–1.5 shows >15% phase-speed discrepancy.

C12. Verdict: WEAKENED

Strongest counterargument: Under short-fetch, light-to-moderate winds in typical ponds (U_10 ≤ 4 m/s, F ≤ 20 m), the JONSWAP spectral peak sits at f_p ≳ 3 Hz (k_p ≳ 35 rad/m), so even at H = 0.1 m, k_p H ≳ 3.5 > 2 — meaning almost no wind-wave energy exists in kH ≲ 2. Even when peak waves do feel the bottom (k_p H ≲ 2 at higher winds/longer fetch), caustic-governing slope (m_2) and curvature (m_4) moments reside at kH ≫ 2 and remain virtually invariant unless bottom friction alters breaking or nonlinear cascades.

Failure condition: A depth sensitivity sweep (H = 0.1 m vs. ∞) for a U_10 = 3 m/s, F = 15 m pond shows <1% change in σ, m_2, and m_4.

C13. Verdict: STANDS

Strongest counterargument: Elfouhaily (1997) embeds an empirical Gaussian roll-off around k_m ≈ 370 rad/m that makes continuum m_4 finite; however, it is not ν(T)-dependent, typical cm-scale grids clip inside k_m, and the stream lane's k^{−α} (α ∈ [5, 6]) power-law tail has a strictly divergent continuum integral m_4 = ∫ k^{5−α} dk → ∞ as k_Nyq → ∞.

Failure condition: Grid refinement (2× and 4× finer Δx at operational domain sizes) yields <1% variation in m_4 across both wind and stream lanes without an explicit ν(T)-dependent dissipation cutoff.

C14. Verdict: WEAKENED

Strongest counterargument: Dolcetti & Krynkin's k^{−α} tail and Horoshenkov's turbulence field both model the exact same physical process (bed-generated non-stationary surface turbulence), so overlap is virtually certain — yet auditing only variance (σ^2) attribution under DESIGN-V4-1-4 is an inadequate test. Because Horoshenkov decays exponentially/Gaussianly while DK decays as a power law, variance budgeting can make σ^2 match while severely distorting or double-counting the m_2 (slope) and m_4 (curvature) spectral bands.

Failure condition: An overlap audit passes total σ^2 budgeting within tolerance while showing >2× discrepancy in k^4 S(k) across the shared k ∈ [30, 300] rad/m band.

C15. Verdict: WEAKENED

Strongest counterargument: At low winds (U_10 < 2 m/s), physics deficits (missing wind-input thresholds, U4 gusts, U4 films, U11 finite fetch), measurement uncertainty (uncontrolled biogenic slicks causing order-of-magnitude scatter in field/tank capillary variance), and gate mismatch (comparing steady-equilibrium models to intermittent data) are physically confounded. Forcing every cell into "a single dominant cause class" guarantees arbitrary misattribution if U4 and U11 are not yet active or if slick contamination in the reference dataset is unmeasured.

Failure condition: Two independent diagnostic metrics (e.g., spectral σ^2 vs. line-depth distribution) assign the same low-wind cell to different dominant cause classes, or sensitivity to unmeasured film coverage spans the entire gap between model and gate.

C16. Verdict: WEAKENED

Strongest counterargument: Banner–Phillips wave breaking caps the short-gravity band (k ∼ 10–150 rad/m), NOT the high-k capillary tail (k ≳ 250 rad/m), where surface tension prevents overturning and viscous dissipation (2νk^2), Marangoni damping, and 3-wave/parasitic transfer are the true binding limits. However, in the short-gravity gust band itself, U4's action balance does blow up without U10 breaking dissipation — refuting U10's late placement at position 10.

Failure condition: Raising wind input at fixed fetch in the k > 300 rad/m capillary band saturates via 2νk^2 and Marangoni dissipation well below any geometric Banner–Phillips breaking threshold, while k ∈ [20, 120] rad/m blows up when U10 is omitted.

C17. Verdict: WEAKENED

Strongest counterargument: In a fixed-point (Eulerian) frame, the decorrelation rate is τ_Euler^{−1} ≈ U_s/L_0 + u_*/H. Because gravel-bed streams have L_0/H ∼ 1.8–2.9 and U_s/u_* ∼ 5–15, advective sweeping U_s/L_0 is 2–8× faster than eddy turnover u_*/H! Fixed-point autocorrelation is therefore insensitive to τ ∼ H/u_* and a poor verification metric; the error of frozen advection shows up starkly in Lagrangian (flow-following) autocorrelation and ω–k Doppler-line broadening (Δω ∼ u_*/H).

Failure condition: Computing fixed-point Eulerian autocorrelation at U_s/u_* = 10, L_0/H = 2.5 shows <10% difference in 1/e lag between frozen advection and τ = H/u_*, even while Lagrangian eddies visually slide rigidly without boiling.

C18. Verdict: WEAKENED

Strongest counterargument: While static menisci decay over ℓ_c ≈ 2.7 mm, real contact lines on oscillating pond banks and flowing stream boulders are moving contact lines that undergo stick-slip hysteresis (θ_A ≠ θ_R) and radiate capillary waves far beyond 1 cm into the field. Numerically, on typical Stage 1 grids (Δx ∼ 3–10 mm ≳ ℓ_c), a 2.7 mm exponential boundary layer cannot be resolved on η's grid and aliases into wrong surface normals/caustic halos.

Failure condition: Evaluating U13 on a Δx = 5 mm grid produces >30% error in contact-line slope and caustic ring intensity relative to the analytic Young–Laplace slope, or an oscillating wave against a solid wall fails to radiate capillary waves.

Attacks on the Proposed Build Order (U1 → U2 → U3 → U4 → U5 → U6 → U7 → U8 → U9 → U10 → U11 → U12 → U13)

U3 (Environmental-Parameter Consistency Audit) at position 3 silently invalidates U1 and U2:

What breaks: U1 computes Stokes drift U_S (sensitive to U_10 reference height, u_*, and H) and U2 computes the wind energy-balance and τ(k) table using U_10, u_*, ν(T), γ(T), H. Placing the input-unification audit after U1 and U2 guarantees rework whenever U3 shifts an environmental input.

Fix: Move U3 to position 1 (before everything).

U10 (Breaking / Micro-Breaking Dissipation) at position 10 breaks U2, U4, U5, U8, and U9:

What breaks: U2 explicitly includes "viscous + breaking + parasitic-capillary dissipation" inside S_ds, and U4 builds the gust action-balance layer using U2's τ(k). Without U10 active, U4 cannot saturate short-gravity gusts (C16), U5 triggers parasitic capillaries on unphysically over-steep un-broken crests, and U8/U9 validate an unstable surface that changes once U10 is added at step 10.

Fix: Merge U10 into U2/U4 (before U5, U8, and U9).

U7 (Finite-Depth Correction) at position 7 breaks U1, U2, U5, and U6:

What breaks: U1's U_S, U2's c_g and bottom friction, U5's Sharma–Dean QTF kernel and shallow-water Ursell check (H ≈ 0.6 m), and U6's shear weighting (e^{2kz} vs. cosh/sinh) and dispersion all depend explicitly on depth H via tanh(kH). Building U1, U2, U5, and U6 with deep-water kernels and patching finite depth at U7 forces rewriting every dispersion-dependent unit downstream.

Fix: Move U7 immediately after U3 (before U1, U2, U4, U5, U6).

U8's numerical convergence, viscous cutoff, and DK/Horoshenkov audit at position 8 invalidates U5 and U6:

What breaks: As C13 establishes, stream m_4 is grid-dependent without a physical viscous cutoff, and cm-scale wind grids clip inside k_m. Evaluating U5's m_4 contribution before establishing a grid-converged viscous cutoff measures grid noise, not bound-wave physics. Similarly, U6 builds ray advection of the DK tail before U8.5 checks whether the DK tail and Horoshenkov double-count the same turbulence band!

Fix: Split U8: run U8a (grid convergence, viscous cutoff, DK/Horoshenkov attribution audit) before U5 and U6; run U8b (confronting U5's m_4, U6's k_0, and full spectral/curvature suite) after all physics units.

U11 (Finite Fetch + Shore Reflection) at position 11 invalidates U8 and U9:

What breaks: U8 validates pond spectra and U9 disposes of the 4 low-wind cells using an infinite periodic domain, even though low-wind pond cells are fetch-limited and boundary-reflected. Adding fetch boundary conditions to U4's action balance after U8/U9 forces re-running both validation and disposition.

Fix: Move U11's fetch BC into U4 (before U8b and U9), and move U9 to the very end after all physics units (U1–U7, U10–U13) close.

Missed Entirely

1. L5 Dropped Completely (Obstacle Boils, Vortex-Core Dimples, and Wake Scars):

Why it matters physically: L5 was confirmed as a split settlement where Stage 1 owns the three separate surface-deformation representations (boil dome, vortex dimple ∼ −Γ^2/(8π^2 g r^2), wake scar) paired with Stage 2 rock wake drivers. Yet none of U1–U13 implements L5. Without it, rocky streams exhibit only linear potential-flow bed waves and static menisci, omitting the strongest rotational caustic generators in natural streams.

Where it belongs: New Unit U14 in Stage 1 (with an explicit Stage 2 → Stage 1 wake-parameter port).

2. Missing Stage 3 → Stage 1 Input Port for Langmuir Convergence (u_LC) and Thermal Coupling:

Why it matters physically: U4 outputs ∇⋅u_s to advect films (L1), but Craik–Leibovich vortex rolls live in Stage 3. Without an explicit Stage 3 → Stage 1 surface-velocity forcing port, L1 Langmuir windrows cannot enter Stage 1's film transport without faking streaks inside Stage 1.

Where it belongs: U1 / U4 external forcing port (receiving u_LC from Stage 3).

3. Wind-Driven Surface Skin Shear Drift (u_drift ≈ 0.02–0.035 U_10):

Why it matters physically: Surface film advection and relative wind (U_10 − U_surface) depend on the wind-driven viscous/turbulent shear current at z=0, which in ponds is absent from U_S alone.

Where it belongs: U1 (part of U_surface and U_mat).

4. Stream Turbulence–Wave Damping and Ring Broadening:

Why it matters physically: Bottom-generated turbulence in U6/U12 does not merely superpose linearly onto waves; turbulent eddy viscosity damps short waves (∝ ν_t k^2) and turbulent velocity fluctuations scatter stationary waves, physically broadening the Dolcetti–Krynkin ring in k-space.

Where it belongs: U6 + U2 S_ds.

5. Transcritical Flow (Fr → 1) Regularization in Streams:

Why it matters physically: U6 declares blocking and caustics out-of-domain, but gravel-bed streams routinely pass through local Fr ≈ 1 over cobbles and riffles, forming smooth Airy/Pearcey transition zones and non-breaking capillary-gravity shocklets rather than spatial holes.

Where it belongs: U6.

6. Energy-Budget Coupling Between U2/U4 Parasitic-Capillary Dissipation and U5 Ripple Generation:

Why it matters physically: U2 deducts parasitic-capillary dissipation from the gravity band while U5 injects parasitic capillaries onto gravity crests, with no conservation link ensuring the energy dissipated from parent crests matches the energy radiated into parasitic ripples.

Where it belongs: Coupling contract between U2/U4 and U5.

7. Caustic-Critical Validation Gaps in U8 (Cox–Munk Slope, Curvature PDF/Kurtosis, and ω–k Dispersion):

Why it matters physically: Scalar m_4 cannot distinguish caustic-forming crest sharpness from Gaussian noise; U8 omits slope skewness, curvature kurtosis/Gaussian-curvature PDFs, empirical Cox–Munk mean-square-slope vs. U_10 confrontation, and animated ω–k dispersion verification (proving bound waves and parasitic capillaries propagate locked to parent crests rather than linear dispersion, and stream turbulence broadens around ω = k⋅U).

Where it belongs: U8.

Final Verdict

FALLS — The design drops confirmed leftover L5 (obstacle boils/dimples/scars) completely, mis-splits Wilton resonance vs. parasitic capillaries (U5), uses an invalid symmetric Plant growth rate for gust decay (U4) and out-of-domain hole-punching for WKB failures (U6), misses Stage 3 → Stage 1 Langmuir forcing, wind skin drift, turbulence wave-damping, and Cox–Munk / ω–k / curvature-kurtosis validation (U8), and orders U3, U7, U8a, U10, and U11 after the units that mathematically depend on them (reordering required: U3 → U7 → U1 → (U2+U10) → U8a → (U4+U11) → U5 → U6 → U12 → U13 → U14(L5) → U8b → U9).
