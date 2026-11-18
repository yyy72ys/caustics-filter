# RALLY-V4-1-14 — GPT 脚

- 脚: GPT（ChatGPT）
- スレッド: https://chatgpt.com/c/6af588e8-e5bc-83ee-bfab-e1320ae60012
- 取得: 2026-11-14 16:53–17:01 JST

---

Adversarial review — v4-1.14
FALLS — decomposition and build order

Central finding: The thirteen-unit decomposition is not yet physically closed, and the proposed sequence contains circular dependencies. The most serious problems are not missing cosmetic mechanisms. They are inconsistencies between the energy-balance model, the nonlinear surface model, the validation measurements, and the physical meaning of the emitted velocity fields.

Three problems deserve particular attention:

1. U2 attempts to establish a band-wise energy balance before several of the processes contributing to that balance exist.
2. U5 attempts to validate curvature-sensitive physics before the physical high-wavenumber limit is established.
3. U6 and U7 separate two pieces of dispersion physics that should be developed against a common finite-depth, current-dependent dispersion framework.

There is also a more fundamental issue: several units describe diagnostics or output ports as though they were physical closures. Emitting a velocity field does not establish that the velocity is correct. Measuring a balance residual does not close the balance.

Part I — Claim-by-claim verdicts

C1. η alone is insufficient; independent kinematic outputs are required

Verdict: STANDS

Strongest counterargument: The proposed definition U_mat = U_surface + U_S is not automatically correct merely because its terms have been given different names.

Stokes drift is a difference between Lagrangian-mean and Eulerian-mean transport. It is not an additional instantaneous velocity that can indiscriminately be added to every surface velocity.

Three quantities must be distinguished:

Instantaneous Eulerian surface velocity.
Oscillatory orbital velocity.
Lagrangian-mean transport velocity, including Stokes drift.

There is also a mass-conservation issue. In a periodic domain, wave-induced Stokes transport can require a compensating Eulerian return flow, depending on the imposed mass-flux constraint.

The proposed zero-net-Eulerian-mean test is therefore insufficient by itself. A wave can have zero imposed mean Eulerian flow while carrying nonzero Stokes transport.

Failure condition: The claim fails if the implementation adds Stokes drift to an already Lagrangian velocity, double-counts orbital motion, or emits a transport velocity inconsistent with the chosen mass-flux constraint.

The foundation is necessary. Its proposed verification needs revision.

C2. Elfouhaily spectrum fails the band-wise source-term balance

Verdict: WEAKENED

Strongest counterargument: The proposed audit risks confusing a stationary spectrum with zero net transfer in every spectral band.

The governing spectral balance has the form ∂N/∂t + ∇_x · (c_g N) + ∇_k · (k̇ N) = S_in + S_nl + S_ds.

Even for a statistically stationary spectrum, nonlinear transfer does not have to vanish in every band.

Nonlinear interactions can transfer energy between bands. A stationary cascade can have nonzero flux through a band while its energy remains constant.

Therefore: S_in + S_nl + S_ds = 0 must be interpreted with the spectral-transfer flux and the relevant transport terms included.

A band-wise residual is not automatically evidence that the prescribed spectrum is physically inconsistent.

A second problem: the Elfouhaily spectrum is a spectral parameterisation, not a guarantee that the chosen combination of Miles/Plant input, nonlinear transfer, viscous damping, breaking, and parasitic-capillary dissipation reproduces its exact equilibrium.

The result depends on the source-term closures.

Failure condition: The claim fails if U2 reports a physically valid nonlinear energy cascade as an imbalance simply because individual bands have nonzero net source terms.

The audit must distinguish local source balance, inter-band transfer, and spatial transport.

C3. Gust response uses band-wise relaxation times from Plant growth rates

Verdict: WRONG

Strongest counterargument: Plant growth rates are not automatically gust-relaxation times.

A wind-input growth rate describes the growth of wave energy under a specified wind-wave configuration.

For a simplified exponential growth law, dE/dt = 2βωE, the corresponding growth timescale is τ_growth = 1/(2βω).

But this is not necessarily the relaxation time after a wind-speed change.

Actual relaxation depends on:

Wind-input adjustment.
Nonlinear spectral transfer.
Dissipation.
Wave propagation out of the forcing region.
The pre-existing wave spectrum.
The gust's duration and spatial extent.

The approximately 45-second observation at λ = 0.2 m is evidence against instantaneous equilibrium. It does not establish that Plant's growth timescale predicts the measured relaxation.

There is another logical issue: U2's proposed diagnostic table must somehow produce a relaxation time that U4 can consume. A residual from a stationary spectral balance is not a relaxation eigenvalue.

Failure condition: The claim fails if the Plant-derived timescale differs materially from measured post-gust energy relaxation, or if the same τ(k) is used for physically different wind histories that produce different responses.

Required correction: Separate wind-input growth rate, spectral relaxation time, and spatial propagation time. Derive or measure the latter two rather than renaming the first.

C4. Film damping must be Marangoni-type and resonance-localised

Verdict: WEAKENED

Strongest counterargument: The rejection of arbitrary uniform broadband damping is correct, but the proposed replacement is too narrowly specified.

Film-induced damping depends on the film's constitutive behaviour, not just its concentration.

Relevant properties include:

Surface elasticity.
Surface shear and dilatational viscosity.
Soluble versus insoluble surfactant behaviour.
Adsorption and desorption.
Surface diffusion.
Bulk-to-surface exchange.
Frequency and wavelength.

The resonance-band observation around 3–10 Hz is not a universal fixed interval. Its location and width depend on wave conditions and film properties.

Moreover, film effects need not enter only as an additional dissipation term. A film can alter the dispersion relation and the effective interfacial boundary conditions.

A scalar concentration field with one damping-rate lookup is not yet a complete physical film model.

Failure condition: The claim fails if measured damping varies with frequency or concentration in a way that the proposed resonance-band model cannot reproduce, or if a nominally identical concentration produces different damping because the film's rheology differs.

C5. Surface convergence and strain are Stage 1 state quantities

Verdict: STANDS

Strongest counterargument: The proposed port is necessary, but divergence and strain must be calculated from the correct velocity field and on the actual interface.

For a deforming free surface, a planar approximation to ∇·u_s may not equal the full surface divergence.

The surface transport equation must account for the geometry of the interface.

For film concentration, the relevant conservation law is of the form ∂c/∂t + ∇_s · (c U_mat) = S_c.

This includes both material transport and surface-area change.

There is also an architectural dependency: Stage 3 owns the three-dimensional Langmuir circulation, but Stage 1 must receive the appropriate surface velocity and convergence forcing from that structure.

Simply emitting strain calculated from the existing wave field does not implement Langmuir circulation.

Failure condition: The claim fails if the emitted strain field is derived from an incomplete velocity field, or if Stage 3's circulation cannot influence Stage 1's surface transport.

C6. Surface-tension-corrected Sharma–Dean QTF covers gravity waves; Wilton region goes to parasitic capillaries

Verdict: WRONG

Strongest counterargument: The proposed separation confuses a mathematical resonance problem with a separate physical wave-generation mechanism.

Adding surface tension to a denominator is not sufficient to establish a valid gravity-capillary second-order QTF.

The complete interaction kernel must be derived consistently from the governing equations and boundary conditions, including:

Finite-depth effects.
Gravity-capillary dispersion.
Sum- and difference-frequency interactions.
Resonant triad conditions.
The appropriate solvability conditions near resonance.

The Wilton resonance is a resonant interaction involving harmonics whose phase speeds coincide under gravity-capillary dispersion.

It cannot simply be declared covered by parasitic capillaries.

Parasitic capillaries are a physical wave population associated with steep parent gravity-wave crests. They are not a generic mathematical replacement for a singular second-order kernel.

In fact, using the parasitic-capillary mechanism to cover the excluded QTF neighbourhood risks double-counting or missing the resonant energy-transfer mechanism.

The problem is especially serious because the resonance lies inside the stated pond wavelength band.

Failure condition: The claim fails if a numerically regularised QTF produces discontinuities, incorrect energy transfer, or divergent amplitudes near the Wilton condition, or if the separately generated parasitic capillaries cannot reproduce the measured resonant triad structure.

Required correction: The Wilton neighbourhood needs its own coupled resonant-interaction treatment or a formally justified exclusion. It cannot be filled by simply switching on parasitic capillaries.

C7. Bound waves are a leading candidate for the fine-line excess; m4 decides

Verdict: WEAKENED

Strongest counterargument: The mechanism is plausible, but the proposed curvature argument is not sufficiently established to make bound waves the leading explanation.

For a surface spectrum, curvature variance is strongly weighted toward high wavenumbers: m4 ∼ ∫ k^4 E(k) dk for an appropriately defined isotropic spectral convention.

That makes m4 a relevant diagnostic.

However, the claimed amplification of approximately 2ka is not a universal result. It depends on the actual QTF coefficient and the amplitude convention.

For a second harmonic, curvature introduces a factor of four relative to its height amplitude because its wavenumber is doubled. The resulting relative curvature contribution depends on the second-harmonic amplitude, not merely on its existence.

There are competing explanations for fine-line excess:

Parasitic capillaries.
Inadequate high-wavenumber dissipation.
Surface-film effects.
Spectral overlap or double-counting.
Insufficient grid convergence.
Incorrect phase relationships.
Breaking-related short-wave saturation.

An increase in m4 after adding bound waves would demonstrate that the correction changes curvature. It would not establish that the correction explains the observed excess.

Failure condition: The candidate is rejected if the corrected m4 moves away from independently measured values, or if the fine-line discrepancy remains after accounting for the bound-wave contribution.

Required correction: Require agreement in m4, spectral distribution, and bispectral phase coupling, not just an increase in curvature variance.

C8. QTF adoption requires measured low-rankness

Verdict: WEAKENED

Strongest counterargument: Low-rankness is a computational feasibility question, not a physical admissibility condition.

The project principle allows neutral numerical methods and physically valid mechanisms.

A QTF does not become physically inadmissible because its matrix representation is not low-rank.

If the direct computation is affordable, the QTF remains valid.

Conversely, a low-rank approximation can be computationally attractive while introducing errors in the triads that dominate curvature.

The reconstruction error must also be measured in a physically relevant norm. A small global matrix error can conceal large errors in high-wavenumber contributions to m4.

Failure condition: The claim fails if a high-rank QTF is rejected despite a feasible direct or FFT-based implementation, or if a low-rank approximation passes a matrix norm but fails the m4 and bispectral tests.

Low-rankness should determine the implementation route, not whether the physical mechanism is adopted.

C9. WKB and action conservation are valid for the stream scenes

Verdict: WEAKENED

Strongest counterargument: The criterion L_U ≫ λ0 is necessary but insufficient to establish validity.

WKB validity depends on the spatial variation of the complete local wave solution, not only the current magnitude.

Relevant quantities include:

Current gradients.
Depth gradients.
Wavenumber gradients.
Branch separation.
Proximity to blocking.
Variation in wave amplitude.
Changes in the local eigenstructure.

Near a blocking point, the group velocity approaches zero and the ordinary ray description becomes singular.

At caustics, ray amplitudes can also diverge even when the physical wave field remains finite.

Declaring such regions out-of-domain is honest, but it may remove precisely the flow features that the renderer is intended to represent.

There is also a physical issue with branch merging: two branches cannot simply be treated as independent populations through a branch-coalescence point.

Failure condition: The claim fails if a scene contains physically important blocking or caustic regions that are silently omitted from the verified scope, or if L_U / λ0 is large while other local variation parameters still invalidate WKB.

C10. Wind input uses relative wind U10 − U_surface

Verdict: WEAKENED

Strongest counterargument: Relative wind is necessary, but the proposed subtraction is not a complete wind-input model.

The air-side wind profile must be interpreted relative to the moving water boundary.

However, the wind input depends on the air-side shear and the wave's intrinsic phase speed, not merely on the difference between a reference-height wind speed and a surface-current vector.

The reference-height velocity must first be translated into the appropriate boundary-layer profile and stress.

There is also a distinction between:

Water-current-relative wind.
Wave-phase-relative wind.
Wind-induced air stress.

These are related but not identical.

The claim that absolute wind necessarily violates action conservation is also imprecise. Action conservation concerns transport and source terms. An incorrect wind-input model produces incorrect forcing; it does not automatically violate the conservation equation if that forcing is included consistently.

Failure condition: The claim fails if relative-wind subtraction predicts the wrong input under opposing currents, especially when wave age and air-side boundary-layer structure materially affect the stress.

C11. Uniform-current approximation requires a band-wise shear criterion

Verdict: WEAKENED

Strongest counterargument: The need for a shear criterion stands, but the proposed effective Doppler velocity is an approximation whose validity needs to be established.

The e^{2kz}-weighted velocity is associated with the depth sensitivity of a deep-water wave's orbital structure.

In finite depth, the vertical weighting changes.

In strong shear, the wave eigenfunction itself changes, so applying a fixed weighting to the undisturbed profile can become inconsistent.

Moreover, a small difference between effective Doppler velocity and surface velocity does not necessarily imply negligible shear effects near a critical level.

The criterion also needs to account for whether the vertical profile contains a level at which the intrinsic wave phase speed matches the local current.

Failure condition: The claim fails if the effective-velocity approximation classifies a materially sheared flow as uniform while the full sheared-flow dispersion relation predicts a significant change.

C12. Finite-depth correction changes headline statistics for kH ≲ 2

Verdict: WEAKENED

Strongest counterargument: Finite-depth physics matters, but the threshold kH ≲ 2 does not establish that every headline statistic changes measurably.

At kH = 2, tanh(2) ≈ 0.964.

The corresponding gravity-wave phase-speed correction relative to the deep-water approximation is only about 1.8%.

At lower kH, the correction grows substantially.

Therefore, the claim needs to be evaluated band by band, weighted by the actual spectral energy and by the statistic being measured.

A small dispersion correction in a low-energy band may have negligible influence on total variance but a measurable influence on a phase-sensitive or slope-sensitive quantity.

There is another problem: TMA is a depth-limited spectral formulation, not a universal physical finite-depth correction that can simply be applied to every arbitrary scene.

Failure condition: The claim fails if the proposed finite-depth correction produces no statistically distinguishable change in the headline statistics over the actual scene distribution, or if TMA disagrees with the appropriate finite-depth measurements.

The deep-water approximation is not universally neutral, but neither is the proposed threshold a universal error boundary.

C13. Without a physical viscous cutoff, m4 is grid-dependent

Verdict: STANDS

Strongest counterargument: A physical high-wavenumber limitation is necessary, but a single hard cutoff is not automatically the correct physical representation.

Because curvature variance weights the spectrum by k^4, even relatively weak high-wavenumber energy can dominate m4.

If the spectrum has no physically justified high-wavenumber roll-off, the computed m4 can increase as resolution increases.

This makes the numerical grid part of the physical result, which violates the project principle.

However, real dissipation is generally distributed across wavenumbers and mechanisms. Viscous damping is not necessarily a sharp terminal cutoff.

The effective high-wavenumber spectrum may also depend on:

Surface films.
Breaking.
Parasitic capillaries.
Nonlinear transfer.
Finite-depth effects.

Failure condition: The claim would be weakened if m4 converges across resolution and domain-size changes under an independently validated physical spectrum. Otherwise, the current curvature statistics are not established as physical.

C14. DK tail and Horoshenkov component double-count variance

Verdict: WEAKENED

Strongest counterargument: Spectral overlap does not prove physical double-counting.

Two independent physical mechanisms can contribute energy in the same wavenumber range.

Conversely, two nominally separate components can be statistically correlated and produce a cross term even if their individual spectra appear reasonable.

The relevant variance is Var(A+C) = Var(A) + Var(C) + 2Cov(A,C).

The existing quadrature accounting must therefore be tested against the physical cross-spectrum, not merely the overlap of the two marginal spectra.

The DK tail and Horoshenkov turbulence model may represent distinct wave populations, partially overlapping populations, or alternative descriptions of some of the same energy.

That question cannot be answered from their spectral slopes alone.

Failure condition: The claim is rejected if the measured cross-spectrum and variance attribution show no duplicate contribution. It is confirmed if the same physical population is counted in both components.

The audit must permit either result.

C15. Each low-wind failure has one dominant cause class

Verdict: WRONG

Strongest counterargument: The procedure cannot guarantee that a single dominant cause exists.

Physics, measurement, and gate defects can interact.

For example:

A deficient capillary spectrum can make a measurement appear inconsistent.
A measurement-resolution limitation can conceal a genuine physical deficit.
A gate can be too strict for the measurement uncertainty while still exposing a real modelling error.

There may be multiple contributing causes, or the evidence may be insufficient to distinguish them.

Forcing every cell into one dominant class creates false certainty.

Failure condition: The claim fails as soon as a cell has two materially contributing causes or insufficient evidence to distinguish them.

Required correction: Permit mixed-cause and indeterminate dispositions. "Cause unknown" is a legitimate scientific result, not a project-management failure.

The operational/physical/scope distinction can remain unchanged.

C16. Breaking dissipation is the binding constraint on the gust band's high-k tail

Verdict: WRONG

Strongest counterargument: The claim assigns universal dominance to one dissipation mechanism without establishing the regime in which it dominates.

High-wavenumber energy can be limited by:

Viscous damping.
Film-induced damping.
Breaking and microbreaking.
Nonlinear transfer.
Parasitic-capillary generation and dissipation.
Other surface contamination effects.

The Phillips saturation concept is not a universal hard ceiling imposed on every high-k band under every wind and depth condition.

Banner-type breaking criteria identify instability conditions; they do not automatically specify the complete energy dissipation rate.

Also, non-air-entraining microbreaking is not necessarily captured by a criterion calibrated to conventional breaking onset.

There is a deeper inconsistency: U10 proposes changing only the dissipation term while leaving the surface-generation stage unchanged. Active breaking also modifies the surface geometry and the wave's subsequent nonlinear evolution.

Failure condition: The claim fails if measured high-k saturation occurs in conditions where viscous or film damping dominates, or if the proposed breaking term predicts saturation at the wrong level.

C17. Stream turbulence decorrelation follows τ ∼ H/u*

Verdict: WEAKENED

Strongest counterargument: The outer eddy-turnover estimate is a plausible timescale, not a complete temporal correlation model.

τ_outer ∼ H / u* provides a dimensional estimate for large-scale turbulence evolution.

But the surface response is filtered by the vertical structure of the turbulence and by the free-surface boundary condition.

Different spatial scales can have different lifetimes.

The fixed-point Eulerian autocorrelation also contains the effects of advection. It is not equivalent to a Lagrangian eddy lifetime.

An exponential decay with one timescale may be a convenient statistical model, but real turbulence correlations need not be exponential.

Failure condition: The claim fails if measured temporal correlations exhibit multiple decay scales, non-exponential behaviour, or strong dependence on spatial scale and advection velocity that cannot be represented by the proposed closure.

C18. Meniscus is confined to a few capillary lengths

Verdict: WEAKENED

Strongest counterargument: The capillary length is the correct leading scale for a gravity-capillary meniscus, but the claim assumes a simple equilibrium geometry.

The exponential decay is associated with a linearised Young–Laplace problem under appropriate boundary conditions.

Real contact lines can exhibit:

Contact-angle hysteresis.
Contact-line pinning.
Dynamic contact angles.
Sharp geometric corners.
Nonlinear meniscus slopes.
Spatially varying wetting conditions.

These effects can change the local profile.

The claim also says "negligible elsewhere" without specifying a numerical tolerance.

A meniscus with a small residual elevation can still have a meaningful slope or curvature near the boundary.

Failure condition: The claim fails if measured or analytically calculated contact-line profiles extend beyond the chosen few-ℓ_c domain at a level that materially changes the surface geometry.

The boundary-condition approach is appropriate, but the contact-line model cannot be reduced universally to a fixed exponential profile.

Part II — Attack on the proposed build order

O1. U1 → U2 → U3: the foundations are in the wrong order
Reorder required

Breakage: U1 and U2 consume physical parameters before U3 has established what those parameters mean.

For example:

Is wind speed U10 or a velocity measured at another height?
Is the stream velocity surface velocity or depth-averaged velocity?
Which density and viscosity are used?
Is surface tension temperature-dependent?
Which depth is supplied to the dispersion calculation?

U1's orbital and Stokes-drift calculations depend on the wave state and its environmental parameters. U2's source-term balance depends on the same parameters.

If U3 subsequently changes their interpretation, the earlier outputs are no longer verified.

The proposed statement that changes will merely be reported does not repair the dependency.

Required order:

U3 first, establishing the environmental input contract.
U1 second, establishing the kinematic output contract.
U2 third, establishing the baseline energy-budget framework.

U3 should not silently change existing scene values. It should establish their physical interpretation and record discrepancies.

O2. U2 precedes U4, U5 and U10 despite depending on them
Critical circular dependency

Breakage: U2 cannot close the stated energy balance while several of its dissipation and transfer terms are still unspecified.

Its proposed balance contains: S_in + S_nl + S_ds.

But the later units introduce or modify:

U4: film-induced damping.
U5: parasitic-capillary generation and associated dissipation.
U10: breaking and microbreaking dissipation.
U7: finite-depth effects on dispersion and dissipation.
U11: fetch-dependent spatial transport.

Consequently, U2's relaxation table is calculated against an incomplete physical model.

U4 then consumes that table as though it were a validated physical quantity.

This is backwards.

There is also a deeper problem: Plant's growth rate is not a complete relaxation model, as identified in C3.

Required correction: Split U2 into two deliverables:

U2a — Baseline energy-budget architecture: Define the spectral balance equation, source-term interfaces, transfer fluxes, and diagnostic conventions.
U2b — Final energy-budget closure: Complete only after the relevant dissipation, film, breaking, parasitic-capillary, and environmental terms have been specified.

U4 may begin against a provisional relaxation model, but its physical verification cannot be declared complete until U2b is complete.

O3. U8 comes too late for U5
Critical validation dependency

Breakage: U5 proposes a prediction of m4 before the physical high-wavenumber domain has been established.

The problem is circular:

U5 adds second-order bound-wave curvature.
U8 measures m4.
U8 requires a physical high-k cutoff.
That cutoff determines how much curvature the model can contain.

Therefore, the cutoff influences the very prediction U5 is supposed to validate.

Without the cutoff, the predicted m4 may be a property of grid resolution rather than the physical model.

This makes the proposed before/after comparison ambiguous.

Required correction: Split U8 into two units.

U8a — Baseline validation and physical spectral-domain definition: Establish measurement conventions, physical high-k dissipation, grid convergence, and baseline m4.
U8b — Post-correction validation: Compare the U5 result against independent measurements, including m4 and bispectral phase coupling.

U8a must precede U5.
U8b must follow U5.

A single U8 cannot correctly occupy only position 8.

O4. U10 cannot wait until position 10
Critical energy-budget dependency

Breakage: U2 already names breaking dissipation in S_ds, and U4 consumes the resulting relaxation table.

If breaking is not yet implemented, then the U2 balance and U4 action evolution are both incomplete.

This creates a hidden inconsistency:

U2 calculates a balance without the final breaking law.
U4 uses the resulting timescales.
U10 later introduces breaking.
The original timescales are now invalid.

There is another issue. Breaking can change both energy and geometry. A dissipation-only term cannot automatically represent the complete surface state during breaking.

Required correction: Move the breaking-model specification ahead of final U2 closure.

At minimum:

Establish the breaking regime and onset criteria early.
Define the energy-dissipation contribution before U2b.
Define the surface-state consequences before U4 is declared complete.

Leave splash and bubble generation in Stage 2.

U10 may remain a separate implementation package, but its physical contract must exist before U2 and U4 are verified.

O5. U7 comes too late relative to U5 and U6
Critical dispersion dependency

Breakage: U5 and U6 use finite-depth physics before U7 has established the common finite-depth framework.

U5 explicitly requires an Ursell-number check in shallow water. Its QTF also requires a consistent dispersion relation.

U6 requires full dispersion for two branches under a spatially varying current.

U7 later changes the wind lane's dispersion and depth-dependent spectrum.

Although U7 is nominally a wind-lane unit, its underlying dispersion conventions are shared physical infrastructure.

There are two separate risks:

U5 uses a finite-depth criterion based on one dispersion implementation while U6 uses another.
U6 and U7 calculate wave properties using inconsistent depth-dependent conventions.

There is an additional issue: U7's TMA correction is not itself a universal finite-depth physical solver. The underlying finite-depth dispersion relation and the chosen spectral parameterisation must remain distinct.

Required correction: Create a common finite-depth dispersion foundation before U5, U6 and U7.

Then:

U5 consumes it for the QTF and Ursell checks.
U6 consumes it for current-dependent branches.
U7 consumes it for wind-wave depth limitation.

The finite-depth dispersion contract should not be invented independently by three units.

O6. U6 is placed after U5 despite a shared dispersion dependency
Reorder recommended

Breakage: U5's second-order interactions and U6's two-branch dispersion can use inconsistent wave branches and phase relationships if independently implemented.

This is particularly dangerous around gravity-capillary resonance.

The QTF's sum- and difference-frequency interactions must use the same physical dispersion relation that defines the participating wave modes.

U6 is not simply a stream-only implementation detail if its dispersion framework supplies the general gravity-capillary branch structure.

Required correction: Establish the common dispersion and branch-selection infrastructure before the specialised U5 and U6 physics.

U5's nonlinear interaction kernel and U6's spatial ray transport can then be implemented separately.

O7. U11 comes too late for U4
Dependency mismatch

Breakage: U11 supplies the finite-fetch boundary condition to U4's action-balance system.

But U4 is developed and verified before the spatial boundary condition exists.

A spatially uniform wind step-change test does not verify a fetch-limited wave field.

The upwind boundary determines the initial wave action. The downstream spectrum then evolves through propagation, input, nonlinear transfer, and dissipation.

The finite-fetch condition changes the spatial history of the spectrum.

There is also a problem with treating shore reflection as a single coefficient. In general, reflection is frequency-, direction-, and potentially phase-dependent.

Required correction: Define the action-balance transport architecture and its boundary-condition interface before U4's final verification.

The actual fetch and reflection implementation can be developed in U11, but the interface must exist earlier.

O8. U9 depends on U8, but U8's baseline must precede most modifications
Reorder required

Breakage: U9 is supposed to diagnose four failed low-wind cells using U8's rulers.

But if U8 is performed after U4, U5, U6 and U7, the measured baseline no longer cleanly identifies the original failure.

A changed low-wind result could be caused by:

The original capillary-floor prior.
New film damping.
The high-k cutoff.
Bound-wave corrections.
Finite-depth corrections.
Numerical changes.

The diagnosis becomes confounded by accumulated model changes.

Required correction: U8a must establish the baseline before new physics is introduced.

U9 should then use that baseline and the later controlled comparisons to determine each cell's disposition.

A final U9 decision may legitimately occur later, but the baseline measurement cannot.

O9. U12 comes too late to verify the temporal architecture
Reorder recommended

Breakage: U4 and U6 introduce time-dependent transport and action evolution, but U12 adds turbulence decorrelation only after most of the stream dynamics have been developed.

This creates two incompatible temporal assumptions:

The stream's stationary-wave field evolves through action transport.
Its turbulence component remains frozen or follows a separate temporal model.

The combined surface can have incorrect temporal cross-correlations even if each component's marginal variance is correct.

Required correction: Define the stream temporal-state architecture early.

U12's specific turbulence-decay law can be implemented later, but the interface for component lifetimes and cross-component temporal correlations must exist before U6's time-dependent stream validation.

O10. U13 is correctly late in implementation, but its boundary dependency is incomplete

Verdict: WEAKENED

Breakage: U13 assumes Stage 2 supplies a static meniscus shape.

But a moving contact line changes the boundary condition.

If the stage-2 object moves, wets, dries, or changes contact angle, the meniscus cannot remain a fixed shape.

Required correction: Keep the meniscus geometry in Stage 2, but define a dynamic contact-line boundary interface in U13.

U13 may remain late because it is a relatively local boundary effect. Its interface must nevertheless be designed before Stage 2 integration is declared complete.

Revised build order

The proposed order should not simply be reversed. The dependencies need to be split into physical foundations, baseline measurements, physical mechanisms, and final validation.

Recommended dependency structure

U3 — Environmental parameter contract
Establish input definitions and physical units.
U1 — Surface kinematics
Establish Eulerian, orbital, Stokes and material-transport ports.
U8a — Baseline validation and physical spectral limits
Establish baseline m4, high-k dissipation, convergence, spectra and independent measurement conventions.
Common dispersion foundation
Finite depth, gravity-capillary dispersion, branch structure and shear applicability.
U2a — Energy-balance architecture
Establish spectral fluxes and source-term interfaces.
U10 — Breaking physics
Establish the breaking and microbreaking energy and surface-state contributions.
U11 — Fetch and reflection boundary architecture
Establish spatial action transport and boundary conditions.
U4 — Gust action evolution and film state
Use a provisional but physically specified energy model.
U5 — Bound waves, resonant interactions and parasitic capillaries
Use the common dispersion foundation and physical spectral limits.
U6 — Spatial current, shear and two-branch action transport
Include wind-current coupling and the relevant propagation boundaries.
U7 — Finite-depth wind-spectrum completion
Complete the wind-specific depth-dependent spectral model.
U12 — Stream turbulence temporal evolution
Integrate with the already-defined temporal-state architecture.
U13 — Contact-line boundary integration
Connect Stage 2 geometry and dynamic wetting conditions.
U2b — Final energy-budget closure
Recalculate relaxation and verify the complete source-term model.
U8b — Final validation group
Validate the completed physical surface against independent measurements.
U9 — Final low-wind disposition
Close, exclude or explicitly leave indeterminate each failed cell.

This is a dependency structure, not a claim that every unit must be completed serially. Some work can proceed in parallel after its physical interfaces are fixed.

The essential change is that baseline measurement, physical cutoff establishment, and common dispersion foundations must precede the mechanisms they are supposed to validate.

Part III — Missed entirely

These are omissions or insufficiently specified couplings that remain even after accepting all thirteen units.

M1. Wave-induced mean return flow and mass conservation
New Stage 1 requirement

Why it matters: U1 emits Stokes drift but does not explicitly close the mean mass-transport budget.

A periodic wave field can have nonzero Stokes transport even when its imposed Eulerian mean is zero.

Depending on the physical boundary conditions, a compensating Eulerian return current or an imposed mass flux is required.

Without this, the emitted material transport velocity can be physically inconsistent.

Where: U1. Add a mass-flux convention and a return-flow consistency test.

Failure condition: Nonzero net transport appears in a domain whose prescribed physical mass flux is zero.

M2. Film physics lacks an equation of state
New Stage 1 requirement

Why it matters: U4 defines concentration, advection and material properties, but does not explicitly define the constitutive relationship between film concentration and surface tension.

It needs something equivalent to: γ = γ(c, T).

It also needs a consistent relationship between concentration, surface pressure, surface elasticity and the resulting tangential stress.

A concentration field alone cannot determine Marangoni forcing.

Soluble and insoluble films can behave differently because material can exchange with the bulk.

Where: U4.

Failure condition: Two physically different films with the same nominal concentration produce the same damping because the model has no representation of their different constitutive properties.

M3. Film feedback on the wind-input mechanism
New Stage 1 requirement

Why it matters: U4 inserts film damping into the action balance, but does not explicitly account for film-induced changes to wind-wave energy input.

A contaminated surface can modify interfacial stress and the air-water momentum-transfer process.

Treating film effects as an additional dissipation term while leaving wind input unchanged can produce an inconsistent energy budget.

Where: U4 and U2b.

Failure condition: Measurements show film-dependent wind input that cannot be explained by the added dissipation term.

M4. Breaking changes surface geometry, not just energy
New Stage 1 requirement

Why it matters: U10 is specified primarily as a dissipation mechanism.

But breaking changes the actual surface state:

Crest deformation.
Crest instability.
Energy redistribution.
Local surface steepness.
Subsequent short-wave generation and dissipation.

Removing energy from a spectral band without representing the corresponding surface-state consequences does not necessarily reproduce the physical surface.

Where: U10, with splash and entrained air remaining in Stage 2.

Failure condition: A breaking event dissipates energy numerically while the generated surface continues to exhibit an unbroken crest.

M5. The stream's turbulence forcing lacks a generation closure
New Stage 1 requirement

Why it matters: U12 addresses turbulence decay and regeneration, but does not sufficiently specify where the turbulence-generated surface spectrum comes from.

The Horoshenkov component supplies a statistical description. It does not automatically determine the correct forcing under every bed roughness, shear velocity, depth and flow condition.

The model needs to distinguish:

Turbulence production.
Transfer of turbulent fluctuations to the free surface.
Surface response.
Advection.
Dissipation and decorrelation.

Otherwise, U12 risks becoming a stochastic lifetime adjustment on a pre-existing pattern.

Where: U12, with bed-generated turbulence parameters supplied by Stage 2.

Failure condition: Two flows with materially different turbulence production produce the same surface turbulence component merely because they share the same nominal decay timescale.

M6. The reflection coefficient in U11 is under-specified
New Stage 1 requirement

Why it matters: A shoreline reflection coefficient is generally not one universal scalar.

It can depend on:

Frequency.
Incident angle.
Shore geometry.
Bottom slope.
Surface conditions.
Breaking and dissipation near the shore.

Reflection also changes phase, not merely amplitude.

A partial standing wave requires a coherent relationship between incident and reflected waves.

Where: U11, with shoreline geometry supplied by Stage 2.

Failure condition: A single coefficient reproduces reflected amplitude but predicts incorrect interference patterns or frequency-dependent reflection.

M7. Stage 3 Langmuir forcing is not fully connected to Stage 1
Cross-stage interface missing

Why it matters: The design correctly places three-dimensional Langmuir circulation in Stage 3, but U4's convergence and strain ports do not fully specify the incoming circulation state.

Stage 1 needs to receive the resulting surface velocity and relevant forcing fields.

It also needs to distinguish those fields from the wave-induced orbital and Stokes velocities.

Otherwise, the same physical transport can be counted twice or omitted.

Where: Stage 3 → U1/U4 interface.

Failure condition: Activating Stage 3 Langmuir circulation changes the material-line field but does not consistently change the surface transport state.

M8. Spatially varying surface tension and temperature
New Stage 1 requirement

Why it matters: U3 includes environmental γ(T) and ν(T), but the architecture appears to treat them primarily as scene-wide parameters.

Spatially varying temperature or film concentration can produce spatially varying surface tension.

That affects capillary-wave propagation and interfacial stresses.

Thermal surface deformation can be negligible while thermocapillary forcing remains relevant to the velocity field.

Where: U3 supplies the fields; U4 consumes the film and surface-stress fields. Temperature-field generation remains Stage 3.

Failure condition: Spatial gradients produce measurable capillary or tangential-stress effects that the uniform-parameter model misses.

M9. Finite-fetch growth and periodic-domain coherence are not jointly validated
New validation requirement

Why it matters: U11 tests empirical growth curves, but that does not establish that the periodic computational domain reproduces the physical spatial coherence of a finite-fetch surface.

A periodic tile can repeatedly introduce the same wave groups and phases.

Even if the one-dimensional spectrum matches the expected fetch-limited spectrum, the spatial and temporal correlations may be wrong.

Where: U11 + U8b.

Failure condition: The surface's correlation functions or wave-group recurrence change materially when the domain length is increased.

M10. No explicit interface contract for the three surface components
New Stage 1 architecture requirement

Why it matters: The established stream architecture has three components:

Bed response.
Stationary dispersive waves.
Turbulence.

U4–U6 introduce additional transport and forcing mechanisms, but the design does not fully specify how their phase, energy and covariance relationships interact with those three components.

This is particularly important because DESIGN-V4-1-4's quadrature accounting must remain intact.

A budget can prevent silent zeroing without preventing physically incorrect cross-correlations.

For example, an independently generated turbulence field can have the correct variance and still have the wrong phase relationship to the stationary-wave field.

Where: U2b, U4, U6 and U8b. Establish a common component-interaction contract.

Failure condition: The total spectrum or variance is correct while the cross-spectrum and phase relationships disagree with independent measurements.

M11. No explicit test of the combined surface's physical boundary conditions
New validation requirement

Why it matters: The units test individual mechanisms, but the final surface must satisfy the combined free-surface conditions.

At minimum, the combined solution must be consistent with:

The kinematic free-surface boundary condition.
The normal-stress boundary condition.
The tangential-stress boundary condition.
Mass conservation.
Appropriate bed and shoreline boundary conditions.

These are not independent decorations. Film stresses, waves, currents, and breaking all interact through these conditions.

The project does not need to solve the complete Navier–Stokes equations everywhere, but its reduced-order models must not contradict the governing boundary conditions.

Where: New U14 — combined boundary-condition consistency audit.

Failure condition: Individual units pass their isolated tests, but their combined surface implies incompatible velocity, pressure or stress conditions.

Part IV — What I would change before implementation

The most dangerous design decisions are concentrated in five places.

Priority / Design defect / Required action

Critical — U2 calculates relaxation before physical source terms are complete — Split U2a/U2b
Critical — U5 validates curvature before physical high-k limits — Split U8a/U8b
Critical — U5 substitutes parasitic capillaries for Wilton resonance — Develop a separate resonant-interaction treatment
Critical — U10 adds breaking dissipation without complete surface-state consequences — Define breaking geometry and energy effects
High — U1 treats Stokes drift as a simple additional velocity — Close the mass-flux and velocity definitions

Two additional decisions should be made before coding:

Whether the project intends to represent active breaking or explicitly exclude it from the verified domain.
Whether the QTF is computationally feasible and physically valid near the gravity-capillary resonance. Low-rankness alone cannot decide this.

The other units are worth retaining, but they should not be treated as already physically closed simply because their proposed verification tests have been written.

Final verdict
FALLS

v4-1.14 is not the right final decomposition or build order. Retain most of its thirteen physical work packages, but reorder the foundations, split U2 and U8 into baseline/final deliverables, move breaking and fetch boundary specifications ahead of energy-balance verification, establish common finite-depth dispersion before U5/U6/U7, replace the Plant-growth-rate shortcut, resolve Wilton resonance separately, and add mass-flux closure, film constitutive physics, breaking geometry, component cross-correlation and combined free-surface boundary-condition validation.
