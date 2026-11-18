# RALLY-V4-1-COMPLETE — GPT leg（verbatim）

- 日時: 2026-11-14 16:28–16:30 JST
- スレッド: https://chatgpt.com/c/6af588e8-e5bc-83ee-bfab-e1320ae60012 （流れ水面ラリー継続）
- プロンプト: RALLY-V4-1-COMPLETE-PROMPT.md
- 最終判定: **FALLS**

---

Adversarial review — v4-1 completion claim

Initial finding: the claim is not established. The fact that the major surface-generation mechanisms have been assigned does not prove that the surface-state model is physically closed.

The critical distinction is between:

A mechanism that changes the actual water-surface geometry.

A mechanism that changes a material property of the interface without necessarily changing its geometry.

A mechanism that transports material along the interface without being representable as a height field.

A mechanism that belongs to a different stage but still supplies boundary conditions to Stage 1.

A physics-first renderer must not silently turn the second and third categories into height perturbations. Equally, it must not omit them merely because the existing height-field architecture cannot represent them.

L1 — Langmuir circulation and material-line transport
MUST-DO in v4-1

Strongest counterargument: The fact that material lines are not a height field is not a valid reason to exclude their underlying transport physics from the surface-state model.

Langmuir circulation is a real wave-current interaction. The Craik–Leibovich mechanism couples Stokes drift and vorticity through the vortex force. It produces counter-rotating circulation cells and surface convergence zones.

The candidate correctly identifies two distinct outputs:

Surface convergence and divergence, which modify the surface flow.

Material-line accumulation, which describes the transport of foam, films, and floating debris.

However, there is an important correction: the spacing is not universally equal to water depth. Observed Langmuir streak spacing exhibits a hierarchy of scales; values of order three times the mixed-layer depth are commonly reported.

What must be implemented: A surface velocity/convergence field and a material-advection capability, with Langmuir forcing activated only under appropriate wave, wind, and turbulence conditions.

Do not draw streaks directly into the height field.

Failure condition: Omitting this mechanism fails for wind-driven water scenes in which Langmuir circulation is physically active and surface tracers or their influence on surface-wave statistics matter.

If the deliverable is strictly geometric caustics with no foam, film, debris, or wave-modulation coupling, the material-line renderer itself can be deferred. That does not justify claiming the underlying surface-flow physics has been implemented.

L2 — Gust structure and independent surface-film scalar
MUST-DO in v4-1

Strongest counterargument: The two mechanisms must remain separate, and the existing equilibrium-spectrum approach cannot substitute for either one.

The 45-second relaxation time at λ = 0.2 m is particularly damaging to any assumption that the local wave spectrum instantaneously follows the local wind.

A wind gust is not simply a spatially varying amplitude multiplier.

It changes wave energy through wind input, propagation, dissipation, and spectral redistribution. A wave-action balance is the appropriate physical framework when propagation and spatially varying forcing matter.

The film mechanism is independently important.

Surface films alter tangential stress through interfacial elasticity and Marangoni effects. Their damping can be strongly frequency-dependent and resonant. They are not equivalent to gust-induced changes in wind input.

Two implementation traps:

A band-dependent relaxation time must be derived or calibrated for the relevant conditions. Do not extrapolate the 45-second observation indiscriminately across all wavelengths.

A film scalar needs physical transport and constitutive properties. A grayscale film mask with an arbitrary damping multiplier is not a physical film model.

What must be implemented: A non-equilibrium wind-wave energy/action evolution mechanism and an independent surface-film state with concentration, transport, and physically grounded damping.

Failure condition: The claim fails for gusting wind or contaminated/film-covered surfaces if the spectrum is simply regenerated from instantaneous wind speed or if film damping is folded into an arbitrary wind-dependent coefficient.

For permanently steady, clean-water scenes, these mechanisms can be disabled by the scene conditions, but their absence cannot be presented as universal physical closure.

L3 — Second-order bound waves and parasitic capillaries
MUST-DO in v4-1

Strongest counterargument: A first-order random-phase spectrum does not reproduce the complete geometry of finite-amplitude gravity waves.

The Sharma–Dean second-order theory supplies sum- and difference-frequency interactions. These produce bound-wave components that are phase-correlated with the primary waves and affect crest shape, asymmetry, and spectral organization.

This is not a cosmetic correction.

A second-order surface has different crest elevations, curvature distributions, and local gradients from a first-order Gaussian superposition.

There is also a separate error to avoid: treating parasitic capillaries as merely the high-frequency tail of the same second-order correction.

Parasitic capillaries are a distinct physical process. They occur near sufficiently steep gravity-capillary wave crests and can remain phase-locked to the generating wave through phase-speed or group-speed matching.

The implementation should distinguish:

Bound gravity-wave corrections.

Free gravity-capillary waves.

Crest-generated parasitic capillary waves.

Their respective dissipation and energy transfer.

Do not add a generic high-frequency noise layer and call it parasitic capillarity.

What must be implemented: A second-order correction for the applicable gravity-wave conditions, plus a separately gated parasitic-capillary mechanism.

Failure condition: The first-order-only model fails when wave steepness is sufficient for bound-wave contributions to materially affect surface geometry. The absence of parasitic capillaries fails when the carrier-crest curvature exceeds the relevant generation threshold.

This does not mean every scene needs a fully nonlinear free-surface solver. It means the model needs a physically justified regime boundary and a valid correction within its intended operating range.

L4 — Spatially varying stream current
MUST-DO in v4-1

Strongest counterargument: A uniform-current dispersion relation is not a valid general description of a spatially varying stream.

The Doppler-shifted dispersion relation depends on the local current projection:

ω = σ(k, H) + k · U(x, y).

Spatial gradients change the local wavenumber and direction. For steady current fields, absolute frequency is conserved along rays, while wave action, rather than wave energy, is the relevant conserved quantity.

The proposed two-branch treatment is useful, but there is a danger of treating two algebraic roots as two universally valid independent wave populations.

Branches can approach one another, coalesce, and encounter blocking conditions. Their existence and admissibility depend on the local dispersion relation and propagation direction.

There is another limitation: conventional wave-action ray transport assumes slowly varying background conditions. If the current varies significantly over a wavelength, the WKB approximation becomes questionable.

What must be implemented:

Spatially varying current.

Local full dispersion.

Branch selection and propagation direction.

Wave-action transport.

Refraction and blocking treatment.

A regime check for the slowly varying-current approximation.

Failure condition: The uniform-current implementation fails around jets, eddies, shear regions, channel contractions, and spatially varying flow fields where current gradients materially alter wave propagation.

A particularly important additional issue is vertical shear: depth-averaged U(x, y) does not necessarily represent the velocity profile experienced by waves. That is a separate outstanding item below.

L5 — Rock and obstacle surface signatures
belongs elsewhere — Stage 2, with Stage 1 surface coupling

Strongest counterargument: The proposed events are real surface phenomena, but their forcing originates in the interaction between the flow and submerged geometry.

The distinction matters:

A bed-response wave is a relatively coherent response to bed topography.

A boil is an upward turbulent event.

A dimple is associated with a vortex approaching or attached to the free surface.

A scar is associated with surface convergence and downwelling.

These cannot all be represented by the existing stationary bed-response term.

The literature explicitly distinguishes dimples, boils, and scars as different free-surface imprints of coherent subsurface motions.

The proposed vortex-core formula,

h ∼ Γ² / (8π² g r²),

is a scaling estimate, not a universal surface-deformation solution.

Its use requires assumptions about vortex geometry, circulation distribution, core structure, depth below the interface, and the free-surface boundary condition.

Correct ownership: Stage 2 owns the rock/obstacle geometry and associated hydrodynamic forcing. Stage 1 must accept the resulting surface forcing and represent the resulting surface deformation.

Failure condition: Treating obstacle-induced turbulence as already covered by stationary bed waves fails whenever intermittent wakes and coherent vortex structures are present.

Do not implement a generic animated dimple stamp in Stage 1. That would replace a physical forcing mechanism with a visual effect.

L6 — Thermal convection in ponds
belongs elsewhere — Stage 3, with Stage 1 boundary coupling

Strongest counterargument: Assigning thermal schlieren optics to Stage 3 does not settle the thermal effect on surface geometry.

Temperature gradients can generate buoyancy-driven circulation. Surface heating and evaporative cooling can also drive horizontal convection and surface currents. Experimental work documents closed circulation patterns under horizontal thermal forcing.

There are two distinct effects:

Thermal fields change optical properties inside the water: Stage 3.

Thermal fields generate flow and modify the free-surface boundary conditions: Stage 3 forcing coupled into Stage 1.

The second effect cannot be dismissed as an optical-only phenomenon.

There is an additional mechanism that should not be overlooked: temperature-dependent surface tension can produce Marangoni stresses. A thermal gradient can therefore affect the interface through both buoyancy and tangential surface stress.

Correct ownership: Stage 3 owns temperature and volumetric convection. Stage 1 receives the resulting surface velocity, stress, and deformation boundary conditions.

Failure condition: The omission matters for ponds with substantial thermal gradients, surface heating, evaporative cooling, or thermal circulation.

For an isothermal, well-mixed pond, the thermal-convection contribution can legitimately be disabled.

L7 — Four unresolved low-wind cells
MUST-DO in v4-1

Strongest counterargument: An operational prior is acceptable as an engineering decision, but it is not evidence of physical correctness.

The existing status is:

σ² = 0.003

with four low-wind cells never passing the validation gate.

Calling the prior permanent does not make those cells resolved.

The distinction is:

Operational closure: The implementation has a defined behavior and does not fail.

Physical closure: The behavior is supported by an adequate physical model or measurements.

Scope exclusion: The project explicitly refuses to claim validity in those conditions.

These are three different statuses.

The prior can remain indefinitely if it is documented as an unresolved operational assumption. However, the project must not describe those conditions as physically validated.

The low-wind regime deserves particular attention because capillary-wave generation, wind input, and viscous damping can be comparable in importance. A small error in the short-wave spectrum can have a disproportionate effect on fine surface slopes and caustic structure.

What must be done:

Preserve the existing prior.

Record the four failed cells individually.

Establish whether the failure is caused by insufficient physical modelling, inadequate measurements, or a validation defect.

Define an explicit uncertainty envelope and operating-domain restriction.

Either close the physical question or formally exclude the affected cells from the validated domain.

Failure condition: The claim fails if the prior is silently promoted from an operational fallback to a physically established low-wind law.

The prior itself does not need to be deleted. The unresolved status does need to remain visible.

L8 — Benjamin–Feir / modulational focusing
not needed — retain the exclusion conditionally

Strongest counterargument: The exclusion is physically defensible for the stated scene ensemble, but the explanation must not become a universal claim that damping and bandwidth eliminate modulational instability.

The classical Benjamin–Feir mechanism requires a sufficiently narrow-banded, predominantly unidirectional wave train. It transfers energy into sidebands and can produce wave-group modulation.

Segur et al. (2005) demonstrated that particular forms of dissipation can stabilize the instability, including in laboratory experiments.

The project has three relevant pieces of evidence:

Broad spectral bandwidth.

Wide directional spreading.

Damping-dominated dynamics.

The existing groupiness result also does not establish Benjamin–Feir instability. Groupiness within a same-spectrum random-phase ensemble is compatible with ordinary linear stochastic wave superposition.

The important distinction is that groupiness is an observed spatial statistic, while modulational instability is a dynamical mechanism.

Why not needed: Under the established operating conditions, implementing a dedicated Benjamin–Feir evolution mechanism would add complexity without demonstrated physical necessity.

Failure condition: The exclusion fails if a scene enters a narrowband, steep, predominantly unidirectional regime with sufficiently weak damping and sufficient evolution time.

A local wave packet can satisfy those conditions even if the global spectrum is broad. Therefore, the exclusion should be conditional on local bandwidth, directionality, steepness, and growth timescale, not just global ν.

Missed items — not in L1–L8

These are additional gaps that prevent a defensible completion declaration.

M1. Wave breaking and microbreaking
MUST-DO in v4-1

Why it matters: The wind-wave model cannot claim physical coverage of natural wind-driven surfaces while omitting the breaking regime.

Breaking is not merely a splash or bubble effect. It changes the surface geometry and dissipates wave energy.

Importantly, microbreaking can occur without conspicuous whitecaps or substantial air entrainment. It affects short-wave energy and surface roughness.

Correct ownership:

Stage 1: breaking onset, surface deformation, wave-energy dissipation.

Stage 2: splash, bubbles, entrained air and associated objects.

Failure condition: Wind conditions permit breaking, but the generator continues producing an unrestricted linear or weakly nonlinear wave field.

If breaking is explicitly excluded from the operating domain, that restriction must be recorded.

M2. Complete wind-wave energy budget
MUST-DO in v4-1

Why it matters: An Elfouhaily-form spectrum is not equivalent to a complete dynamical model of wind-wave energy evolution.

The physical spectrum depends on competing mechanisms:

Wind input.

Viscous dissipation.

Nonlinear wave-wave transfer.

Breaking dissipation.

Generation and attenuation of parasitic capillaries.

These mechanisms dominate different spectral regions.

L2 addresses non-equilibrium gust response, but it does not automatically establish that the underlying energy budget is closed.

Correct ownership: Stage 1.

Failure condition: The implementation changes spectral amplitudes without accounting for where energy enters, moves between bands, and leaves the wave field.

A stationary statistical generator can use a validated equilibrium spectrum without integrating every source term. But it must not claim to simulate physically consistent temporal evolution unless those mechanisms are accounted for.

M3. Wind waves interacting with stream currents
MUST-DO in v4-1

Why it matters: L4 addresses spatially varying stream currents, but the architecture still needs an explicit answer to what happens when wind waves and flowing water coexist.

Adding an Elfouhaily wind field to a stream surface is not generally equivalent to simulating their combined physical state.

Current changes the intrinsic and absolute frequencies, wave propagation, wavenumbers, and energy distribution. Opposing currents can increase wave steepness and eventually block propagation.

This is a cross-component interaction, not a third independent height layer.

Correct ownership: Stage 1.

Failure condition: Wind-wave and stream components are independently generated and superimposed despite significant current-induced changes to their dispersion and propagation.

M4. Vertical shear in the stream velocity profile
MUST-DO in v4-1

Why it matters: L4 uses a spatially varying horizontal current U(x, y), but this does not settle the vertical velocity profile U(z).

A depth-averaged current is not necessarily the velocity experienced by every wave component.

Vertical shear changes the wave dispersion relation and can introduce important differences from the uniform-current solution.

Correct ownership: Stage 1.

Failure condition: The model applies uniform-current dispersion to a materially sheared flow.

The implementation need not solve arbitrary three-dimensional turbulence everywhere. It does need a justified criterion for when a depth-uniform current approximation is acceptable.

M5. Surface kinematics must be consistent with surface geometry
MUST-DO in v4-1

Why it matters: A height field alone is insufficient to describe all the physical state variables required by the other mechanisms.

The system needs to distinguish:

Surface elevation.

Surface velocity.

Wave orbital velocity.

Mean current.

Wave-induced mass transport.

Stokes drift.

Material transport velocity.

These velocities are not interchangeable.

This omission affects Langmuir circulation, material-line transport, wave-current interactions, and the evolution of turbulence-generated structures.

For example, using the same velocity for material advection and wave phase propagation would be physically unjustified.

Correct ownership: Stage 1.

Failure condition: A downstream mechanism needs surface kinematics, but the generator only supplies height and spatial gradients.

This is an architectural consistency requirement, not necessarily a demand for a full 3-D fluid solver.

M6. Spectral and temporal validation beyond RMS and groupiness
MUST-DO in v4-1

Why it matters: The current validation results establish important properties, but they do not establish complete physical agreement.

The existing tests cover:

Surface variance.

Groupiness.

Line-depth distributions.

Numerical warp behavior.

Still missing are independent validations of:

Directional spectrum — Correct directional energy distribution

Frequency spectrum — Correct temporal wave content

Wavenumber spectrum — Correct spatial scales

Frequency-wavenumber spectrum — Correct dispersion and propagation

Cross-spectrum between components — Correct phase relationships and coupling

Surface autocorrelation — Correct temporal coherence

Slope distribution — Correct local ray-direction statistics

Curvature distribution — Correct fine-scale geometry

Crest/trough statistics — Non-Gaussian and nonlinear wave structure

The stream validation has an especially important limitation: the reported k0 is predicted by the stationary-wave equation, not independently measured across all 13 conditions.

The surface-standard-deviation gate therefore cannot establish the correctness of the stationary-wave spectrum.

Correct ownership: Stage 1 validation.

Failure condition: The model passes RMS and groupiness tests while producing incorrect spectral, temporal, or geometric statistics.

M7. Wind forcing and environmental parameter consistency
MUST-DO in v4-1

Why it matters: A physically named spectrum is not enough. Its inputs must correspond to the physical scene.

The project needs a consistent mapping between:

Reference-height wind speed.

Wind friction velocity.

Wind direction.

Wave age or fetch assumptions.

Water depth.

Density.

Dynamic viscosity.

Surface tension.

Temperature-dependent water properties.

This is particularly important when moving between ocean-like wind waves, shallow ponds, and streams.

A spectrum calibrated for one environmental regime should not be assumed valid in another simply because the same nominal wind speed is supplied.

Correct ownership: Stage 1, with environmental parameters supplied by the scene/environment system.

Failure condition: Two physically different water states produce the same generated spectrum because the generator ignores a parameter that materially changes dispersion, forcing, or dissipation.

M8. Numerical convergence of the water-surface generator
MUST-DO in v4-1

Why it matters: The warp-lane audit does not establish convergence of the underlying surface simulation.

The following remain distinct numerical questions:

Grid resolution versus the highest represented wavenumber.

Spectral aliasing.

Domain-size dependence.

Periodic-boundary artifacts.

Phase evolution under discretization.

Energy conservation through spectral truncation.

Gradient and curvature convergence.

Sensitivity to the numerical cutoff.

The highest-frequency capillary components are particularly vulnerable.

A surface may have the correct RMS height while its slope and curvature statistics change substantially with resolution.

Correct ownership: Stage 1 numerical validation.

Failure condition: Increasing grid resolution or changing domain size materially changes the surface statistics or downstream ray geometry without a physically justified explanation.

Final verdict
FALLS

The statement "Stage 1 of v4 is complete; there is nothing left to do in v4-1" is not defensible.

The highest-priority remaining work is:

Wind-wave energy evolution, gust response, and breaking.

Second-order gravity waves and parasitic capillaries.

Wind/current coupling and nonuniform-current propagation.

Surface kinematic consistency.

Low-wind unresolved-cell disposition.

Independent spectral, temporal, and numerical validation.

L5 and L6 should be implemented through their respective Stage 2 and Stage 3 ownership, not duplicated as independent surface generators. L8 is a defensible exclusion under the established operating conditions.

The most important architectural objection: v4-1 currently has a collection of physically motivated surface components, but component coverage is not equivalent to physical closure. The project still needs to demonstrate that the components interact consistently, that their energy budgets do not contradict one another, and that the resulting surface reproduces independently measured geometry and dynamics.

Until those conditions are met, Stage 1 can be described as substantially specified with an operational implementation, not physically complete.
