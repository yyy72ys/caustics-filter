# RALLY-V4-3-1-GROK — 第3段設計ラリー Grok 脚（verbatim）

- 脚: Grok（端末ブラウザ経由、新規会話）
- 対象: RALLY-V4-3-1-PROMPT.md（DESIGN-V4-3-1-DRAFT の主張 C1–C10 への敵対的査読）
- 日時: 2026-11-14 16:39–16:43 JST
- 会話URL: https://grok.com/c/0d9742a6-7027-49bb-974c-5d89a950832c?rid=f74cbe4f-bdfc-45ec-a845-dbb5eb0f6f5e
- 注記: 以下は返答本文の verbatim 保存。改変なし。

---

C1. STANDS
Strongest counterargument: Incoming solar rays (sun → bed) and outgoing rays (bed → eye) are not symmetric under the same n-field; a pure vertical stratification affects Snell angles differently on the two legs, and a plume’s ∇n deflects the deposit location before the view-side ray even starts. Omitting the sun-side therefore systematically misplaces the energy that stage 1 and the deposit quadrature have already conserved.
Failure condition: A controlled experiment (or high-resolution BOS + caustic imaging) in which a known bed-heated plume produces measurable lateral shift of the caustic pattern on the bed while the camera-side image of a fixed bed mark remains unshifted; if that occurs, one-sided transport is incomplete.

C2. WEAKENED
Strongest counterargument: The 1–4 mm figure is a derivation check, not a measurement. Apparent displacement scales with path length × (dn/dT)×ΔT×geometric leverage; at the main verification depth (~0.6 m) and realistic cool-skin / bed-ΔT values (a few K), many published water schlieren / BOS results show sub-millimetre or sub-pixel shifts once the ray is nearly normal and the integration path is short. Pixel threshold therefore depends on the exact camera sampling and on whether the plume is viewed through the full depth or only a fraction of it. Choosing bed heating solely because the number exceeds a pixel is circular until the number is measured in the actual geometry.
Failure condition: Direct BOS or schlieren imaging of a bed-heated plume at 0.6 m depth under the project’s illumination geometry yields <1 pixel rms displacement; the “first-field” ranking collapses.

C3. WEAKENED
Strongest counterargument: Secondary-flow surface patterns were admitted as labelled inputs only after they were shown to be kinematically constrained by the free-surface boundary conditions already solved in stage 1. A volumetric plume field is not so constrained; its positions, spacing and ΔT are free parameters whose joint statistics are fixed by the Boussinesq equations that the project has already ruled out of scope. Labelling the field “input” does not convert an unconstrained statistical guess into physics.
Failure condition: Two different parametrisations of the same surface energy balance and bed albedo produce statistically distinguishable caustic deposits or veil images under otherwise identical transport; the field is then an arbitrary texture, not a labelled physical input.

C4. STANDS
Strongest counterargument: The FF phase function is a local material property of the scatterers. Once the local optical depth and the local b/c are known along a bent path, evaluating the same CDF kernel at those local values is the correct first-order generalisation; the kernel shape itself does not acquire new angular structure merely because a and b vary slowly in space. Path-integral bookkeeping (optical-depth accumulation, local sampling of the kernel, path-radiance accumulation) is the only necessary change.
Failure condition: A delta-tracking MC reference on a strongly inhomogeneous field shows that the angular shape of the single-scattering halo deviates from the local FF evaluation by more than the homogeneous-case residual already accepted in v16.

C5. STANDS
Strongest counterargument: Path radiance is an integral of the source function along the ray; the top-down view already defines a unique set of rays from bed to camera plane. Changing the camera later only re-points the terminal segment; the volume contribution accumulated up to the surface remains valid and can be transformed by the new camera mapping. Building the term against the fixed view is therefore not premature.
Failure condition: When stage 4 is introduced, the re-pointed eye-side rays produce a veil image whose spatial structure cannot be obtained by any re-mapping of the stage-3 veil field; the volume integral is then camera-dependent in a non-trivial way.

C6. STANDS
Strongest counterargument: HydroLight (and any plane-parallel invariant-imbedding solver) assumes horizontally infinite layers and a uniform bottom. The project geometry is a finite patch with spatially structured bed reflectance and finite lateral extent of the caustic deposit. Those features generate horizontal gradients that plane-parallel codes cannot represent; a delta-tracking MC on the actual 3-D fields is the appropriate reference.
Failure condition: A plane-parallel calculation on the horizontally averaged fields reproduces the MC deposit energies and halo containment to within the project’s existing tolerance; the geometric distinction is then irrelevant at the verification sizes.

C7. STANDS
Strongest counterargument: Optical depth cL ≈ 0.2 is the regime in which the single-scattering series is known to converge rapidly for ocean-water phase functions. At cL → 1 the second- and higher-order contributions become measurable; the only rigorous way to locate the crossover for the project’s specific fields and phase function is the controlled order-capped versus uncapped MC comparison already proposed. Drawing the boundary by “taste” would violate the measurement-first rule already applied to the phase function itself.
Failure condition: The order-capped MC at cL = 0.2 already differs from the uncapped run by more than the homogeneous-stack residual; single-scattering is then inadequate even for today’s scenes.

C8. STANDS
Strongest counterargument: A halocline maintained by continuous freshwater inflow is a genuine material discontinuity on the scale of the optical path; treating it as a step (Snell refraction at a sharp interface) is the continuum limit of the GRIN equation. A thermal plume is a continuous temperature field whose gradients are set by diffusion and advection; staircasing it replaces the true ∇n with artificial jumps and alters the ray curvature. The split is therefore forced by the different physical origin of the two interfaces.
Failure condition: High-resolution refractive-index tomography of a laboratory salt wedge shows a continuous transition zone whose thickness is optically significant at the project wavelengths; the step approximation then fails.

C9. WEAKENED
Strongest counterargument: The Ellison–Turner / Turner E(Ri) curves were obtained in laboratory flumes at Re ≲ 10³. Field observations of atmospheric and oceanic gravity currents at Re ~ 10⁷ show systematically higher entrainment rates and a different Ri dependence (approximately Ri^{-3/4}). A 0.6 m creek is intermediate; assuming the laboratory curve transfers without modification is an untested extrapolation, not the same status as Dolcetti’s surface data which were already at comparable scale.
Failure condition: Direct measurement of entrainment velocity or interface deepening in a 0.6 m-class salt-wedge flume or field creek yields E values outside the laboratory Ellison–Turner envelope; the closure must be re-calibrated.

C10. WEAKENED
Strongest counterargument: The governing principle requires that the fields themselves be products of actual physics or neutral numerics. Reduced models that supply time-varying (n,a,b) can satisfy the principle only if every free parameter is pinned by measurement or by a higher-fidelity closure that has itself been verified. Once the reduced models contain tunable spatial statistics (plume spacing, front sharpness, renewal times), the division of labour becomes an appearance-driven approximation unless those statistics are independently constrained. Holding the transport exact while the fields are free parameters does not automatically reproduce nature.
Failure condition: Two different reduced-model realisations that are observationally indistinguishable at the surface produce statistically different caustic or veil images; the volume state is then under-constrained and the division of labour fails the principle.

Missed entirely
Polarisation / birefringence of the scattering. Even for unpolarised sunlight the FF phase function is polarisation-dependent; the deposit and the veil acquire polarisation structure that a scalar RTE discards. Matters for any later polarimetric verification or for the correct energy in the halo. Belongs in 3.2 (or forces a vector RTE upgrade).
Surface-wave focusing / defocusing of the volume paths. Stage 1 already supplies a time-varying free surface; rays that enter through different surface facets sample different vertical paths through the same (n,a,b) field. The design treats the volume as if the surface were flat after the initial refraction. Physically first-order for any non-flat surface. Belongs in the coupling between stage 1 and 3.1/3.2.
Temporal correlation of the fields with the surface. Plume rise and front motion are not independent of the surface strain field already computed in stage 1; the design evolves them separately and only feeds a roughness modulation back. The missing two-way coupling can shift both the caustic locations and the veil on the 8 s loop time scale. Belongs in 3.4.
Finite-size effects of the sun disc inside the volume. The deposit already uses disc quadrature at the bed; once rays bend, the effective solid angle changes along the path and the path-radiance integral must be performed over the same finite disc. The design inherits the homogeneous-medium treatment. Belongs in 3.1/3.2.
Validation gap: absolute radiometric scale of the veil. Energy conservation and relative halo containment are planned, but no absolute comparison of path radiance against a calibrated in-water radiance measurement is required. Without it the veil term can be arbitrarily scaled while still passing the relative tests. Belongs in 3.4.

Final verdict on the claim under test: FALLS — the four-block shape is directionally correct, but C2, C3, C9 and C10 are already weakened by missing scale checks and unconstrained reduced-model parameters; the missed surface–volume coupling and polarisation make the decomposition incomplete for caustic image formation.
