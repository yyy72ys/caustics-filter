# RALLY-V4-1-14 — dsh leg (DeepSeek Harness, pc-dsh)

Seat: dsh (DeepSeek Harness headless, user PC DESKTOP-UPSCHGD)
Prompt: RALLY-V4-1-14-PROMPT.md (PC copy: D:\mutaro-work\rally-in\RALLY-V4-1-14-PROMPT.md, 24,063 bytes, size-verified)
Dispatched: 2026-11-14 16:51 JST; completed 2026-11-14 ~16:54 JST (156 s)
PC save: D:\yohey2026\dsh\muse\out\review_v4_1_14.md (24,438 bytes, fetched verbatim below)

---

# RALLY-V4-1-14 adversarial review

## Verdicts C1–C18

### C1 — η alone insufficient; material transport must be U_mat = surface flow + Stokes drift as independent field (U1)
**Verdict: STANDS**
**Strongest counterargument:** For linear / weakly nonlinear waves orbital velocity and Stokes drift are deterministic functionals of the directional spectrum, so emitting them as fields adds I/O, not physics; a pseudo-velocity from η via linear theory (u ≈ ∇φ) coincides with the Lagrangian mean to O((ka)²), and in stream-dominated scenes the Eulerian current dominates U_S by an order of magnitude, making the error from an η-derived proxy quantitatively negligible there. U1's real value is interface discipline, not new physics.
**Failure condition:** Kills it if a film-patch / drifter advection test in a wave-dominated, current-free pond shows centroid displacement predicted from an η-derived pseudo-velocity matches the U_mat prediction to within measurement error, or if monochromatic closed-form checks show U_S from the emitted field disagrees with a³k²σ theory.

### C2 — Prescribed Elfouhaily does not balance band-by-band under standard S_in/S_nl/S_ds; audit will show gust-band imbalance (U2)
**Verdict: STANDS**
**Strongest counterargument:** S_ds (breaking + viscous + parasitic) contains order-of-magnitude tunable coefficients, and S_nl has multiple competing forms (DIA vs exact); with tuned breaking coefficients almost any prescribed spectrum can be made to "close," so the predicted imbalance is a statement about the auditor's choice of closures, not about nature. Without frozen S_in/S_nl/S_ds forms the claim is unfalsifiable.
**Failure condition:** Kills it if U2, using documented ST4/ST6 + DIA + viscous closures with no ad-hoc tuning, closes within its own stated tolerance in every band including the gust band in the deep-water fully-developed limit.

### C3 — Gust response is band-wise relaxation with τ(k) from Plant; instantaneous equilibrium already refuted by ~45 s at λ=0.2 m (U4)
**Verdict: STANDS**
**Strongest counterargument:** Plant rates are scatter-factor-~2 empirical fits, and the 45 s number is one band / one facility / one wind history; the fast capillary bands relax in <1 s via viscosity while nonlinear transfer can shortcut linear growth, so a single linear τ(k) table over-predicts lag for the short waves and under-predicts memory for swell. Band-wise relaxation is the right skeleton but Plant-τ as the quantitative table is provisional.
**Failure condition:** Kills it if a step-wind test shows band energies reaching the new equilibrium in ≪ τ(k) (e.g. capillaries in <2 s where Plant predicts >10 s), or instant-modulation predicts the relaxation curve as well as τ(k) does.

### C4 — Film damping must be Marangoni resonance-band (3–10 Hz) action dissipation set by coverage × elasticity; uniform broadband damping is wrong (U4)
**Verdict: STANDS**
**Strongest counterargument:** Real films are mixed soluble/insoluble/biological slicks with elasticity varying over decades and patchy coverage; the resonance frequency and peak width move with elasticity, and surfactant-altered surface boundary conditions plus entrainment change viscous dissipation outside the textbook 3–10 Hz peak. A misplaced resonance with wrong elasticity is no more physical than a uniform coefficient — both are wrong without in-scene film characterisation.
**Failure condition:** Kills it if controlled film experiments show >2× damping persisting an octave or more outside the modelled resonance band at fixed coverage, or uniform damping fits measured fetch-growth suppression as well as the resonance model.

### C5 — Surface convergence/strain field is a stage-1 state output; material look is stage 2 but advection field is stage 1 (U4, L1 port)
**Verdict: STANDS**
**Strongest counterargument:** For near-incompressible surface flow ∇·u_s is zero at leading order (rigid-lid); non-zero divergence comes from second-order Stokes divergence, upwelling, and numerical differentiation noise. Emitting ∇·u_s and strain components on η's grid risks shipping differentiated noise as "physics," and stage 2 could compute the same derivatives from the velocity fields it already receives.
**Failure condition:** Kills it if the emitted ∇·u_s in representative wind + stream scenes is shown to be numerically dominated (changes sign / magnitude under grid refinement while velocities converge), or if film-advection driven by stage-2-computed derivatives is indistinguishable from stage-1-emitted fields.

### C6 — Sharma–Dean QTF with tension in denominator covers gravity band; Wilton neighbourhood excluded to parasitic capillaries (U5)
**Verdict: WEAKENED**
**Strongest counterargument:** Sharma–Dean is a pure-gravity second-order result; inserting tension into the denominator is an ad-hoc hybrid, not a derived capillary-gravity second-order kernel. The Wilton resonance (λ≈1.7 cm) sits inside the pond's 1–20 cm main band, so the "excluded neighbourhood" swallows the band centre the design most needs, and parasitic capillaries (phase-locked, crest-triggered, dissipative) are a different physical mechanism from resonant second-order bound harmonics — one cannot cover the other's excluded domain by declaration. Near Wilton second order itself breaks down and higher order / fully nonlinear treatment is required.
**Failure condition:** Kills it if the tension-augmented QTF compared against a fully nonlinear reference (high-order spectral or tank Stokes-wave + capillary-gravity triads) disagrees by more than tolerance even one octave away from Wilton, or if the declared Wilton-excluded bandwidth exceeds ~50% of the gravity-capillary overlap band.

### C7 — Bound waves are a leading candidate for the fine-line excess; m4 before/after decides (U5+U8)
**Verdict: WEAKENED**
**Strongest counterargument:** Curvature amplification ~2ka relative to height ratio cuts both ways: at typical steepness ka≈0.05–0.15 the bound height fraction is O((ka)) so the m4 fraction is O((ka)²) ≈ 0.25–2% — an order of magnitude too small to explain a factor-level line-density excess unless steepness is extreme everywhere. Parasitic capillaries and the unresolved high-k tail / viscous-cutoff choice move m4 far more. As a "candidate to be tested" the logic is sound; as "leading candidate" it over-claims.
**Failure condition:** Kills it (as leading candidate) if U8 measures the U5-induced Δm4 at <10% while the fine-line excess requires ≥50% more curvature variance, or if surfaces with U5 switched off but parasitic capillaries on reproduce the excess.

### C8 — QTF adoption conditional on measured low-rankness first (U5)
**Verdict: STANDS**
**Strongest counterargument:** The gate may be the wrong gate: broadband QTFs are generically not low-rank (rank grows with tile size and triad count), so the measurement is likely to fail and kill a route that FFT-based convolution (O(N log N) per triad class) could have delivered cheaply. Conditioning on low-rankness mistakes one implementation premise for physical admissibility.
**Failure condition:** Kills it if eigendecomposition shows reconstruction error above tolerance at any affordable rank while direct FFT convolution meets the cost budget — the QTF route then stands despite failing the low-rank gate, proving the gate mis-specified.

### C9 — WKB + action conservation valid for stream scenes (L_U≫λ0); out-of-domain declaration at blocking/caustics correct (U6)
**Verdict: WEAKENED**
**Strongest counterargument:** The scenes where U6 matters most — flow over rocks, boil/dimple/scar neighbourhoods, confluences — have L_U comparable to λ0 by construction (obstacle scale ~ decimetres, ring λ0 ~ centimetres to decimetres). WKB is valid in the open reach where the answer is boring and invalid exactly where amplitude closure is needed; flagging blocking zones and caustics "out of domain" leaves holes in the rendered surface rather than a surface. Correct behaviour for a diagnostic, insufficient behaviour for a generator without a stated infill (matched asymptotics, full-wave patch, or calibrated local model).
**Failure condition:** Kills it if mapping L_U/λ0 over representative stream scenes shows L_U<5λ0 over more than a few percent of visible water, or if declared out-of-domain patches coincide with the boil/scar footprints the design must render.

### C10 — Wind input with relative wind U10−U_surface; absolute wind violates action conservation measurably on opposing current (U6)
**Verdict: STANDS**
**Strongest counterargument:** Stress (hence u* and S_in magnitude) itself depends on relative wind via the drag law; subtracting U_surface only inside S_in's argument while taking u* from absolute-wind tables double-counts or under-counts. The measurable violation claim also requires a strong opposing current (U/c_p not small), rare in pond scenes and modest in streams except at boil cores.
**Failure condition:** Kills it if a uniform opposing-current test (U≈0.3–0.5 m/s against λ≈0.2 m waves) shows absolute-wind vs relative-wind action-balance difference below the design's own tolerance, or if coupled stress-corrected absolute wind closes action equally well.

### C11 — Uniform-current approximation only where band-wise U_eff(k) deviation below threshold; shear criterion required (U6)
**Verdict: STANDS**
**Strongest counterargument:** The e^{2kz}-weighted U_eff(k) is a monochromatic linear construct; for a broadband directional spectrum with wind drift shear concentrated in the top millimetres, every short band "fails" while the long bands pass, and the threshold itself is free. Without a prescribed threshold tied to a dispersion-error budget the criterion is a label, not a test; directional shear and inflection-point instability escape it entirely.
**Failure condition:** Kills it if bands passing the criterion still show >tolerance dispersion error against a full Euler shear solver (e.g. Kirby–Chen / Dong–Kirby), or if all capillary bands fail in every scene, making "uniform is fine" vacuous.

### C12 — Finite-depth correction changes wind-lane headlines measurably where kH≲2; deep-water form not neutral there (U7)
**Verdict: STANDS**
**Strongest counterargument:** In H=0.05–1.0 m scenes the waves that feel depth (kH≲2 ⇒ λ≳πH) are the long gravity tail, which carries variance σ but little slope/curvature variance; capillaries and the gust band (λ~cm–dm) remain deep-water even in the shallow pond. Depth-limited saturation (TMA φ(kH)) is further offset by fetch limitation in small ponds, so the net headline move in slope/m4 — the quantities tied to the open fine-line problem — may be within noise even while σ moves.
**Failure condition:** Kills it if H-sensitivity runs (H=0.05, 0.2, 0.6, 1.0 m, fixed wind/fetch) change σ, slope variance and m4 all by less than the design's significance threshold.

### C13 — Without physical viscous cutoff m4 is grid-dependent; current curvature statistics not established as physical (U8)
**Verdict: STANDS**
**Strongest counterargument:** If the prescribed Elfouhaily tail plus viscous dissipation already imposes an effective high-k roll-off inside the resolved band, m4 converges with refinement and grid-dependence is an implementation bug (cutoff placed beyond dissipation range), not a fundamental un-established status. The claim over-generalises from a missing-cutoff implementation to all current surfaces.
**Failure condition:** Kills it if a grid/domain/cutoff convergence study shows m4 stable to <10% under 2× refinement in each dimension with the existing tail, establishing current m4 as physical without any new cutoff.

### C14 — DK tail and Horoshenkov component double-count overlap-band variance (U8 audit item)
**Verdict: WEAKENED**
**Strongest counterargument:** Overlap in wavenumber is not double-counting: the DK tail is coherent stationary (phase-locked to bed/flow, narrow in frequency–wavenumber), Horoshenkov is stochastic turbulence-forced (broad, decorrelating). Under the DESIGN-V4-1-4 quadrature accounting they add as coherent + incoherent terms; variance adds only where both claim the same (k,ω) coherence cell. The audit may well find zero overlap once frequency–wavenumber coherence is used instead of k-magnitude.
**Failure condition:** Kills it (in the strong double-count form) if f–k coherence attribution shows the coherent fraction in the overlap band accounts for the DK variance with residual incoherent variance matching Horoshenkov independently, i.e. summed variance matches measurement with no correction.

### C15 — 4 low-wind-cell failures separable into one dominant cause class (physics/measurement/gate); disposition will not end "unknown" (U9)
**Verdict: WEAKENED**
**Strongest counterargument:** Low-wind scenes couple the three classes irreducibly: capillary-floor shape uncertainty, unmeasured surfactant contamination, anemometer / reference-height error at <2 m/s, and gate thresholds tuned on higher-wind data. Promising a single dominant cause per cell plus "never unknown" is a project-management pledge, not a physical claim; mixed-cause cells (e.g. 40/35/25 splits) are the expected outcome and forcing a single label manufactures false certainty.
**Failure condition:** Kills it if U8 re-measurement attributes two cells to mixed causes with no class exceeding ~50% share, or any cell remains within joint uncertainty of all three classes after the rulers are applied.

### C16 — Breaking dissipation is the binding constraint on the gust-band high-k tail; balance without it over-grows short waves (U10)
**Verdict: WEAKENED**
**Strongest counterargument:** At gust-band steepness below the Banner–Phillips threshold breaking probability is near zero; viscous + parasitic-capillary + S_nl drain dominates the short-wave sink. Phillips saturation describes the equilibrium high-wind tail, not a transient gust front, and micro-breaking rates are uncertain to an order of magnitude — invoking them as "the" binding constraint pre-empts the U2 audit that is supposed to decide the sink partition quantitatively.
**Failure condition:** Kills it if U2/U4 balances with viscous + parasitic + S_nl sinks (breaking set to zero) already saturate the tail at the Phillips level under gust forcing, or saturation tests show tail growth linear in wind input even with the breaking term on.

### C17 — Turbulence decorrelation τ∼H/u* measurably differs from frozen advection; frozen time statistics are wrong (U12)
**Verdict: STANDS**
**Strongest counterargument:** Fixed-point autocorrelation is dominated by wave orbital advection (periods ≤1 s) and mean-flow sweeping, both faster than outer-scale eddy turnover (H/u* ~ seconds to tens of seconds); the frozen-vs-decorrelating difference may sit below wave-induced variance. Keeping total variance stationary while imposing exp(−t/τ) decorrelation also imposes a specific regeneration spectrum that is itself unvalidated.
**Failure condition:** Kills it if fixed-point autocorrelation with frozen Taylor advection vs τ=H/u* decorrelation differ by less than measurement / sampling noise over the scenes' record lengths, with variance held equal.

### C18 — Meniscus confined within a few ℓc (≲1 cm); boundary-condition treatment correct (U13)
**Verdict: STANDS**
**Strongest counterargument:** Cumulative shoreline / rock / floating-debris contact perimeter in pond/stream scenes is large, and surfactant-dependent contact-angle hysteresis plus meniscus scattering of short capillary waves extend influence beyond the static Young–Laplace decay; slope-variance integrals over many metres of contact line can contribute non-negligibly even if each meniscus is locally confined.
**Failure condition:** Kills it if integrating the analytic meniscus slope field over a representative scene's contact perimeter changes scene slope variance or near-shore line density above tolerance, or measured near-shore profiles deviate from exponential ℓc decay (hysteresis / dynamic contact line).

## Attacks on the order

1. **U8's rulers — especially the viscous-cutoff decision — must precede U5's prediction, not follow it.** U5's headline deliverable is Δm4 and a bispectral triad test, but C13 concedes m4 is currently grid-dependent without a physical cutoff. A Δm4 computed on a grid-dependent ruler is itself grid-dependent: change the high-k cutoff and the "bound-wave contribution" changes sign/size arbitrarily. U8 item (1) even admits its outcome may force a viscous cutoff as a "separate rally-reviewed change" — i.e. the design schedules the ruler after the measurement it rules. Move the cutoff/convergence audit (U8 items 1+4) before U5, or split U8 into U8a (rulers+cutoff, before U5) and U8b (full suite, after).
2. **U10 cannot wait until position 10.** U2's S_ds and U4's action balance both already name breaking/micro-breaking dissipation as a term, and U4's gust-band evolution is untestable without its sink partition: an action balance with sources but deferred sinks either over-grows (C16's own warning) or is tuned with a placeholder sink that U10 later invalidates, forcing U4's verification (step-wind relaxation, saturation) to be redone. Either pull the breaking-dissipation term forward as part of U2/U4 (leaving only calibration to U10), or admit U4's verification is provisional until U10. As drawn, U4→…→U10 builds the house before pricing the roof.
3. **U7 after U6 and U5 creates rework on both neighbours.** U6 solves stationary-phase k0(x,y) and action closure under deep-water dispersion, then U7 changes the dispersion relation (tanh kH) and saturation underneath it — every k0 field, blocking point (U=c_g(H)), and amplitude must be recomputed. U5's Ursell check (H≈0.6 m class) already needs H, and the QTF denominator itself is depth-dependent in finite depth. U7 is independent of U6 only in deep water, which the scenes are not (kH≲2 admitted). Move U7 before U6 (or fuse U6+U7 as one dispersion/action unit with H as input from the start); at minimum move U7 before U5's shallow-water gate.
4. **U3 placed third can retroactively invalidate U1/U2 inputs.** U3 unifies reference height, u*, wave age, ν(T), γ(T) — every number U2's balance audit and U1's Stokes/transport magnitudes consume. If the audit moves U10→u* by 20% (typical log-profile correction) or changes γ/ν with temperature, τ(k), U_S magnitude, viscous cutoff, and Marangoni resonance all shift, and U1/U2 verification must rerun. Parameter correspondence is zeroth-order work: move U3 to position 0/1 (before U1/U2), with U1/U2 consuming only audited inputs.
5. **U11 needs more than U4's frame, and U6 needs U11's shore.** U11 integrates fetch growth from zero at the upwind shore, but U6's stream inlet and U4's gust initial/boundary conditions already assume developed spectra at domain edges; bolting fetch on at position 11 invalidates the boundary conditions U4/U6 were verified against. Specify fetch/boundary treatment jointly with U4 (wind) and U6 (stream inlet), not as a late add-on. Similarly U12's τ∼H/u* needs U3's u*/H (another reason U3 goes first), and U9's disposition needs U8a rulers plus the U10 sink decision for low-wind capillaries — U9 after full U8 is fine only if U8 is split as above.
6. **Proposed re-order (minimal change):** U3 (audit) → U1 (kinematics) → U8a (convergence + viscous-cutoff ruler) → U2 (balance closure on audited inputs with breaking term included) + U10-term forward → U4 (action+gust+film) → U7 (depth in dispersion/saturation) → U6 (two-branch + shear + action, on finite-depth dispersion) → U5 (bound+parasitic, on finite-depth QTF gate + cutoff ruler) → U8b (full spectral suite, k0, double-count audit) → U9 (disposition) → U11/U12/U13 (fetch+reflection, decorrelation, meniscus ports). U6/U7 may alternatively fuse.

## Missed entirely

1. **The S_nl operator has no home.** U2 diagnoses balance with S_nl, U4 evolves action with wind input + Marangoni + breaking sinks — but no unit implements the nonlinear transfer operator (DIA/exact, bound-vs-resonant partition, depth-scaled) inside the time-evolution layer. Without S_nl the gust band cannot downshift peak frequency (JONSWAP fetch laws cited in U2/U11 verification require it) and energy piles at the forced scale. Belongs in U4 (with U2 specifying the form); currently the action balance is sources minus sinks with no redistribution.
2. **Blocked-action sink.** U6 flags blocking/caustics out-of-domain but never says where blocked action goes (breaking, reflection, transfer to turbulence, upstream radiation). Action conservation along rays that terminate at U=c_g leaks energy unless a sink/source is specified — a conservation violation by omission. Belongs in U6 (+U10 for the breaking channel); verification should assert global action accounting across blocking lines.
3. **Wind-stress coupling, not just relative-wind argument.** U10−U_surface in S_in is half the coupling; u* itself (hence all of S_in's magnitude, τ(k), H/u*, TMA saturation) depends on wave-age/sea-state via the drag law, and gust advection at "its own speed" ignores boundary-layer eddy structure and pressure–slope coupling. A gust field advected as a frozen scalar plus local S_in misses the documented gust–wave phase correlation. Belongs in U4 (stress/gust structure) with inputs from U3; at minimum state the drag closure and its iteration with the wave field.
4. **Wave-averaged vortex-force / Stokes–shear feedback (stage-1↔3 boundary).** U6 advects waves on currents and U4 emits convergence/strain, but waves exert Craik–Leibovich vortex forces and Stokes-drift shear stresses back on the surface current that Langmuir (L1) and shear (U6) assume given. With L1's circulation deferred to stage 3, the design needs an explicit one-way forcing port (wave field → CL force terms for stage 3) or it silently assumes no feedback. Name it in U1/U6 as a stage-3 handoff.
5. **Depth-dependent QTF + bottom-friction partition in U5/U7 interface.** U5 checks Ursell at H≈0.6 m but evaluates a deep-water QTF; in kH≲2 second-order kernels, saturation (TMA φ), and bottom friction redistribute variance jointly — U7's bottom-friction term and U5's bound tail overlap in the long-wave shoaling regime with no stated partition (cf. the DK/Horoshenkov audit in U8, which has no finite-depth counterpart). Add to U7/U5 interface or U8 audit.
6. **Meniscus–film–contact dynamics, not just static meniscus.** U13 holds a static Young–Laplace port with contact angle as scene input, but surfactant (U4's film) changes contact angle/hysteresis, and moving contact lines with waves exhibit dynamic-angle dissipation that damps short waves at shorelines — the very low-wind pond scenes of U9. Couple U13↔U4 (film-dependent angle) and U13↔U11 (reflection coefficient already depends on shoreline dissipation); currently U11's reflection port and U13's meniscus port are two ports to the same shoreline that never meet.
7. **Temporal-validation gap.** U8's suite is spectral + autocorrelation + slope/curvature PDFs, but gust memory (U4), decorrelation (U12), and groupiness interact in time–frequency (bispectral nonstationarity, envelope statistics, peak/trough asymmetry after U5). No unit validates joint (x,t) statistics against measurements (e.g. time-resolved slope image sequences), only marginal spectra. Add time–frequency validation to U8b.
8. **Integration architecture / cost envelope.** U4 time-stepping + U6 ray tracing + U5 full-tile QTF convolution + U12 regeneration interact through resolution, CFL, and stochastic-seed discipline; there is no unit owning the coupled cost/stability budget (U5 gates only QTF rank; U8 checks convergence post-hoc). A generator that is physical per-term but unstable or unaffordable coupled fails the project principle's "neutral numerics" clause. Needs a named integration owner (new small unit or part of U8a).

## Final verdict

**FALLS — decomposition is largely the right inventory, but order must change (U3 first, cutoff/convergence ruler before U5, breaking sink forward into U2/U4, U7 before U6/U5) and the S_nl operator, blocked-action sink, wind-stress coupling, and meniscus–film–reflection joins must be added.**
