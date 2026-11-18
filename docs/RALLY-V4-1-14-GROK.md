# RALLY-V4-1-14 — Grok 脚

- 脚: Grok（grok.com、Fast モデル）
- 会話URL: https://grok.com/c/47942b46-175f-4eca-a7bc-6e2ff49c9b37 （会話タイトル "Adversarial Review of Water Surface Simulation Design v4-1.14"）
- 取得: 2026-11-15 00:03–00:09 JST
- 最終判定: **FALLS**
- 判定内訳: STANDS 15（C1, C2, C3, C4, C5, C7, C8, C9, C10, C11, C12, C13, C16, C17, C18）／WEAKENED 3（C6, C14, C15）／WRONG 0
- 経路注記: **当日例外で musespa29 アカウントの新規会話**（本人の 2026-11-15 指示「今日は一旦29に寄せて」による）。恒久分担では波の部屋の Grok は musesupa6 の既存スレッド https://grok.com/c/cbdd84bd-f592-460b-90bf-e9f864219028 が正であり、そちらにはプロンプトを 2026-11-14 に1回投稿済みだが無料プラン上限で返信が生成されなかった（20:11・23:11 の再試行でも上限のまま）。本脚はその代替として 29 側の新規会話で同一プロンプト全文を投稿して取得したもの。作業後、アカウントは musespa29 のまま（6 へは戻していない）。有料アップグレードは行っていない。

---

C1. STANDS.
Strongest counter: A pseudo-velocity reconstructed solely from ∂tη or from linear orbital formulae derived from η alone necessarily mixes Eulerian and Lagrangian frames and cannot recover the Stokes drift contribution that is second-order in amplitude and independent of the instantaneous height field. Film and material-line advection driven by that pseudo-velocity therefore transport mass at the wrong velocity.
Failure condition: A monochromatic progressive wave in a periodic domain with measured zero Eulerian mean flow but non-zero Stokes drift; the film centroid must still move at the analytic Stokes drift. If it does not, C1 is confirmed; if a pure-η reconstruction somehow recovers it, C1 falls.

C2. STANDS.
Strongest counter: Elfouhaily is a fully-developed deep-water spectral shape, not a solution of any closed action-balance equation under the standard Miles/Plant + Hasselmann + Banner–Phillips/viscous closures. Band-by-band residual S_in + S_nl + S_ds is known to be non-zero, especially in the high-frequency gust band where dissipation time-scales are short.
Failure condition: A deep-water, fully-developed audit that returns residuals inside a published tolerance for every band that the project actually uses; if that occurs, the "does not balance" claim is falsified for those bands.

C3. STANDS.
Strongest counter: The measured relaxation time ≈ 45 s at λ = 0.2 m already falsifies instantaneous equilibrium. Plant growth rates supply the only first-principles τ(k) that is consistent with that measurement.
Failure condition: A step-change wind experiment in which the band energy adjusts on a time-scale statistically indistinguishable from zero (or from the group-travel time across the domain) rather than from the Plant τ(k).

C4. STANDS.
Strongest counter: Marangoni dissipation is a resonant, elasticity-dependent mechanism peaked in the 3–10 Hz band; a uniform broadband coefficient neither conserves action correctly nor reproduces the measured frequency dependence of film damping.
Failure condition: Laboratory monochromatic gravity-capillary damping rates (with and without controlled film coverage) that are reproduced by a single constant coefficient to within measurement error across the resonance band.

C5. STANDS.
Strongest counter: Stage 1 already owns the film field and therefore must own the surface velocity field that advects it. The visual appearance of foam or debris is stage 2; the kinematic driver is not.
Failure condition: A demonstration that film transport can be driven solely by a stage-2 post-process without violating mass conservation or action conservation on the surface state itself.

C6. WEAKENED.
Strongest counter: Even with surface tension restored in the denominator, the classical Sharma–Dean kernel still diverges or becomes singular in a finite neighbourhood of the Wilton resonance; simply declaring that neighbourhood "outside the domain" is an ad-hoc cut, not a derivation. Parasitic capillaries are a distinct physical process, but the precise width of the excluded band is not fixed by the QTF analysis alone.
Failure condition: An eigendecomposition or direct residual evaluation of the tension-corrected kernel that remains well-behaved through the Wilton band for the project's typical steepness; if it does, the forced split is unnecessary.

C7. STANDS.
Strongest counter: Curvature ∇²η sets line density. Bound waves of amplitude ratio ∼ 2ka relative to the free wave therefore amplify m4 by a factor that is first-order in steepness and measurable.
Failure condition: Before/after m4 comparison on the same linear field that shows no statistically significant change once the high-k viscous cutoff is held fixed.

C8. STANDS.
Strongest counter: Low-rank structure is the only practical justification for adopting a full QTF on a production surface mesh. If reconstruction error after truncation exceeds the project's action-conservation or curvature tolerance, the cost premise collapses.
Failure condition: Measured reconstruction error (relative L2 or m4) that remains above tolerance for every rank that is computationally affordable.

C9. STANDS.
Strongest counter: Typical stream length scales L_U are metres while λ0 is centimetres; WKB holds. Near blocking or caustics the ray approximation fails and any interpolation fabricates energy.
Failure condition: A controlled opposing-current experiment in which a wave packet is observed to tunnel or scatter through the nominal blocking line at amplitudes inconsistent with the out-of-domain flag.

C10. STANDS.
Strongest counter: Absolute wind injects energy with the wrong relative velocity; the resulting action flux is inconsistent with the Doppler-shifted dispersion relation under opposing current.
Failure condition: An opposing-current run in which the absolute-wind formulation conserves wave action to the same precision as the relative-wind formulation.

C11. STANDS.
Strongest counter: Vertical shear produces a k-dependent effective Doppler velocity. Without a quantitative band-wise threshold there is no basis for claiming that the surface velocity is an adequate surrogate.
Failure condition: A measured shear profile for which the deviation U_eff(k) − U_s remains below the stated threshold for every band that carries significant variance, yet the surface statistics still diverge from the uniform-current prediction.

C12. STANDS.
Strongest counter: For kH ≲ 2 the finite-depth dispersion and the Kitaigorodskii/TMA saturation factor change both the spectral shape and the integrated slope variance by amounts that exceed ordinary numerical error.
Failure condition: A suite of scenes with H in the project's stated range whose headline statistics (σ, slope variance, peak frequency) change by less than the project's own acceptance tolerance when the finite-depth factors are switched on.

C13. STANDS.
Strongest counter: Without a physical viscous cutoff the high-k tail is set by the grid; m4 therefore inherits that arbitrary scale and cannot be compared with nature.
Failure condition: A grid-convergence study in which m4 asymptotes to a unique value independent of both grid spacing and the numerical high-k filter once the filter is placed beyond the viscous scale.

C14. WEAKENED.
Strongest counter: The DESIGN-V4-1-4 accounting already treats the bed-response term as a separate quadrature; an overlap audit may simply confirm that the residual shared variance is already zero by construction. The claim is therefore still open.
Failure condition: An attribution audit that recovers a statistically significant shared variance band after the DESIGN-V4-1-4 separation has been applied.

C15. WEAKENED.
Strongest counter: The four cells may each be a linear combination of physics deficit, measurement error and gate defect; forcing a single dominant class is an assumption, not a result.
Failure condition: Re-measurement under U8 rulers that leaves residual variance that cannot be attributed to any one of the three classes above a pre-declared threshold.

C16. STANDS.
Strongest counter: Phillips saturation and Banner–Phillips breaking criteria are the only mechanisms known to bound the high-k tail under increasing wind input. An action balance that omits them grows the short-wave energy without limit.
Failure condition: A fixed-fetch, increasing-wind experiment in which the high-k energy saturates at the Phillips level even when the breaking term is switched off.

C17. STANDS.
Strongest counter: Frozen Taylor advection preserves perfect temporal correlation along streamlines; real bottom turbulence decorrelates on the outer time τ ∼ H/u*. Fixed-point autocorrelation therefore differs.
Failure condition: Fixed-point temporal autocorrelation measured on the generated stream surface that is statistically indistinguishable from the frozen case for the project's range of H/u*.

C18. STANDS.
Strongest counter: The meniscus decays exponentially over ℓ_c ≈ 2.7 mm. Outside a few capillary lengths its contribution to slope and curvature is negligible compared with wind and stream waves.
Failure condition: A contact-line experiment in which the meniscus-induced slope remains above the project's capillary-wave detection threshold at distances ≳ 1 cm from the wall.

Attacks on the proposed order

U8 (especially the physical viscous cutoff for m4) must precede or be concurrent with U5. U5's bound-wave prediction is stated as a change in m4. Without a grid-independent viscous cutoff already in place, that change is meaningless: any measured Δm4 could be an artefact of the numerical filter. Building the prediction first and only later discovering that m4 itself is unphysical breaks the verification contract of U5. Move the viscous-cutoff decision (and the associated convergence study) into a pre-U5 foundation unit or make it the first deliverable of U8.

U10 cannot wait until position 10. Both U2's diagnostic balance and U4's time-evolution already list breaking/micro-breaking as an S_ds term. Leaving the concrete Banner–Phillips / Phillips-saturation implementation until after U4 has already been written means the τ(k) tables and the action-balance layer are first built with an incomplete or placeholder dissipation, then later retro-fitted. That is a classic source of silent inconsistency. U10 must be absorbed into U2/U4 (or executed immediately after U2) so that the balance that U4 inherits already contains the ceiling.

U7 after U5 and U6 is late. U5's Ursell-number check and the declaration of second-order domain limits are depth-dependent. U6's two-branch dispersion and blocking condition are also depth-dependent once H is no longer infinite. Performing the finite-depth correction only after both units have already chosen their kernels and out-of-domain flags forces either a re-work of those kernels or an inconsistent hybrid (deep-water QTF + finite-depth dispersion). U7 should sit immediately after the foundations and before U5/U6.

U3 third is too late for the numbers U1 and U2 already consume. U1's kinematics ports and U2's Plant rates depend on the precise interpretation of U10, u*, ν(T), γ(T) and depth. An audit that discovers a systematic mismatch after those modules have been coded does not merely report a difference; it invalidates the numerical values already embedded in the τ(k) tables and the Stokes-drift formulae. The correspondence table must be frozen before any generation or balance code is written.

U11's fetch boundary condition is a boundary condition of U4. Treating fetch-limited growth as a later add-on means U4 is first written for periodic or infinite-fetch domains and only later given an upwind shore. The growth integrals and the shore-reflection port should be part of the initial design of the action-balance layer, not a subsequent patch.

Missed entirely

- Explicit viscous cutoff as a physical high-k bound, not merely a numerical filter. Without it m4 remains grid-dependent (already noted in C13) and the parasitic-capillary generation of U5 has no natural high-wavenumber termination. Belongs inside U2/U8 as a required physical term before any curvature statistic is trusted.
- Directional relaxation and directionality of the gust response. Plant rates and the measured 45 s time-scale are usually quoted for the omnidirectional spectrum; the project never states how the directional distribution itself relaxes under a changing wind vector. A purely scalar τ(k) is incomplete. Belongs inside U4.
- Consistent treatment of the capillary-gravity transition under simultaneous current and finite depth. U6 and U7 are written separately; the combined dispersion relation (two-branch + tanh(kH)) and the associated blocking loci are never derived. A new joint unit or an explicit cross-check after both is required.
- Surface-tension gradient (thermocapillary) coupling to the film field even when the height deflection is negligible. L6 correctly dismisses the static deflection, but a horizontal γ(T) gradient still drives Marangoni stresses that alter the film and the short-wave dissipation. The design treats γ only as a scalar parameter. Belongs as a source term in U4's film equation, fed by stage-3 temperature.
- Validation of the material transport velocity itself against laboratory Stokes-drift measurements. U1 only checks internal consistency. An external confrontation with measured Lagrangian drift is absent. Belongs in U8.
- Coupling of the stage-1 convergence field to stage-2 foam/debris generation. The design correctly emits the strain field, but never states the quantitative map from surface convergence to foam coverage or lifetime. That map is a stage-2 responsibility, yet the interface contract is missing.
- Energy flux through the air–water interface under opposing current. Relative wind is used for S_in, but the corresponding work done on the air side and the consistency of the total energy budget are never audited. Belongs inside the U2 balance.

Final verdict on the claim under test: FALLS. The decomposition captures the confirmed leftovers, but the build order places critical physical constraints (viscous cutoff, breaking ceiling, finite-depth dispersion, parameter freeze) after the modules that already depend on them; those dependencies must be reordered before the sequence can be accepted as the right way to close stage 1.
