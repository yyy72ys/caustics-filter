# RALLY-V4-4-1 — Gemini leg（verbatim）

- 日時: 2026-11-14 16:39–16:44 JST
- 席: Google AI Studio、モデル **Gemini 4 Hydrogen**（gemini-4-hydrogen、Thinking level: High）選択を確認して投稿
- プロンプトURL: https://aistudio.google.com/prompts/1XVHqvrjAgQhRmqpdsOCz5Se4ms9ELrQ6
- プロンプト: RALLY-V4-4-1-PROMPT.md
- 最終判定: **FALLS**
- 判定内訳: STANDS 1（C8）／WEAKENED 5（C1, C4, C5, C6, C7）／WRONG 3（C2, C3, C9）

---

Per-Claim Verdicts (C1–C9)

C1 — WEAKENED

Strongest counterargument: The negative claim (do not substitute defocus by pre-blurring the C field or displacement field) holds, though even its stated coincidence condition is wrong: pre-blurring the displacement field fails even on a flat fronto-parallel bed because the caustic Jacobian is nonlinear (∇⋅(blur(d)) ≠ blur(∇⋅d)), and pre-blurring C ignores bed albedo texture and shadows. More fatally, the positive claim ("computed only as an aperture integral over the post-transport radiance field") contradicts v4-4.2/v4-4.3: if the post-transport output is a 2D single-view image plus a scalar apparent-distance map, a per-pixel PSF convolution is not a physical aperture integral. Across a finite aperture D = f/N, different sub-aperture rays pass through different regions of the refracting water surface (experiencing different wave curvatures), look around foreground edges of Stage 2's 3D bed objects (aperture parallax/partial occlusion), and integrate Stage 3's volumetric scattering over different depths. A 2D PSF over a single post-transport image cannot reconstruct occluded background radiance or multi-depth compositing.

Failure condition: A Stage 2 3D object suspended 0.3 m above the bed (or Stage 3 volumetric turbidity + Stage 1 sun glitter) imaged at f/2.8 focused on the surface: a 2D PSF over the post-transport buffer produces hard silhouette haloing/bleeding around the 3D object and forces surface glitter, volume scatter, and bed radiance to share one blur radius, whereas real aperture rays peek behind the silhouette and focus each depth component independently.

C2 — WRONG

Strongest counterargument: Both supporting assertions fail mathematically. First, using geometric distance does not miss the focal plane by a factor n in air above water: with camera height h_a and water depth d_w, s_geom = h_a + d_w while s_app = h_a + d_w/n. In C7's own main case (h_a = 1.0 m, d_w = 0.6 m, n = 1.333), s_geom/s_app = 1.60/1.45 = 1.103 (a 10% shift, not 33%); at h_a = 3.0 m it is 1.043. The ratio approaches n only as h_a → 0 (lens pressed flush to the water). Second, invoking "the image of the bed is not even a surface" to mandate an apparent-distance field s(u,v) is self-contradictory: a scalar per-pixel distance field s(u,v) defines a single surface by construction and cannot represent an astigmatic two-sheeted focal structure (tangential vs. sagittal focal lines), nor can it represent real-wave refraction where curved crests drive apparent depth through ±∞ into real images floating above the water.

Failure condition: Overhead or oblique capture of a caustic-forming wave crest with principal curvature κ > n(n−1)d_w / n (e.g., wave radius R < 0.15 m at d_w = 0.6 m, or R < 1.0 m in sub at d_w = 4.0 m): paraxially 1/d_app = n/d_w − (n−1)κ ≤ 0, so the apparent image splits by principal curvature into two orthogonal distances—one going through ∞ into a real image in air (potentially behind the lens, s < 0) while the other remains virtual beneath the surface.

C3 — WRONG

Strongest counterargument: E = 2.2 and A = 2.6 cannot both be "recording-side constants." By the settled facts, E = 2.2 is applied to the direct term only, while A = 2.6 is an ambient scale. A physical camera exposure (shutter time × aperture area × sensor gain) is blind to photon provenance—it multiplies total incident radiance (L_direct + L_ambient) by a single scalar. It cannot apply 2.2× gain to solar caustic photons and 2.6× gain to diffuse sky photons. While one common scale factor can represent recording exposure, the ratio A/E = 2.6/2.2 ≈ 1.182 alters the sun-to-sky irradiance balance, which is strictly an environment / light-transport parameter (or an arbitrary look-dev adjustment violating the project's governing principle).

Failure condition: Any physically calibrated radiometric check comparing the direct-to-ambient contrast ratio of the image against clear-sky solar vs. sky downwelling irradiance, or any aperture/shutter sweep: classifying both E and A as camera-side constants implies the camera sensor itself changes the weather's sun-to-sky ratio.

C4 — WEAKENED

Strongest counterargument: Even on its own aberration-free thin-lens terms, v4-4.3's model has three physical/mathematical defects:

Swapped denominator in the blur formula: v4-4.3 writes b = f^2/N |s − s_f| / (s_f (s − f)). The exact thin-lens image-plane blur circle is b = D |v_f − v| / v = f^2/N |s − s_f| / (s (s_f − f)) (with s and s_f swapped in the denominator!). v4-4.3's formula gives the wrong s → ∞ limit (f^2/(N s_f) instead of f^2/(N (s_f − f))) and diverges to ∞ as s → f+, whereas physically collimated rays at s = f yield a finite blur circle equal to the aperture diameter D = f/N. (C7's 0.139 mm figure used this buggy formula; the exact value is 0.141 mm.)

Invalid geometric-plus-Airy composition: Convolving a top-hat disc b with an infinity-focused Airy disc (2.44 λN, which also omits the finite-conjugate working f-number N_w = N s_f/(s_f − f)) violates wave optics. In the transition regime of 1–5 waves of defocus (precisely where C7's f/8–f/16 examples sit), coherent Fresnel diffraction over the exit pupil creates dark/bright central interference spots (e.g., on-axis intensity drops to zero at W_020 = λ) and ringing that incoherent Disk∗Airy smears out.

Omission of non-aberration off-axis thin-lens geometry: Natural cos^4 θ illumination falloff and entrance-pupil foreshortening (cos θ elliptical bokeh) are general geometry of a planar sensor and circular stop, not lens-specific aberrations; omitting them drops corner exposure by ∼48% at 32° half-field (main at 35 mm). Furthermore, "within an error range measurable by PSF-shape verification" is unfalsifiable without a specified tolerance or reference lens class, especially when 10^4× dynamic-range sun glitter amplifies tiny PSF tail differences.

Failure condition: Close-up focus or wave-refracted real image near s → f (v4-4.3's b → ∞ singularity); f/8–f/16 small-defocus regime (W_020 ≈ λ–3λ) where Fresnel on-axis nulls diverge from Disk∗Airy; or 10^4× sun-glitter bokeh in main corners where real elliptical pupil foreshortening and cos^4 θ falloff are absent.

C5 — WEAKENED

Strongest counterargument: Separating viewing-condition CoC conventions (PJ9) from lens blur physics (Stage 4) is architecturally sound in principle, but the design violates its own separation and oversimplifies the output. First, C7 immediately re-imports a hardcoded viewing CoC (c = 0.029 mm) into Stage 4 to declare focus pulls "invisible" and gate a "mandatory" regime map—contradicting C5 and failing whenever the 2048^2 artwork is inspected at pixel scale (1 px = 17.6 μm on 36 mm full frame < 29 μm) or projected large in a gallery. Second, Stage 4's sensor pitch and OLPF/diffraction (v4-4.5) impose a real camera-side spatial-frequency floor before PJ9 ever sees the image. Third, Stage 4 cannot output "the blur-circle diameter itself" as a scalar when the PSF is anisotropically astigmatic (C9) or diffraction-dominated.

Failure condition: Gallery-scale projection or 100% pixel inspection of main (2048^2, 17.6 μm pitch) where blur circles of 20–28 μm (1.1–1.6 px) are classified "in focus" by 0.029 mm yet visibly soften sharp 2–3 px caustic cusps—or an astigmatic surface point whose PSF has two orthogonal diameters (b1, b2) rather than a single b.

C6 — WEAKENED

Strongest counterargument: While internal-focusing multi-element changes in f_eff(s_f) are indeed lens-specific, a unit-focusing ideal thin lens inherently produces focus breathing because the image plane shifts to v_f = f s_f/(s_f − f) as focus distance s_f changes. The angular FOV scales with 1/v_f, meaning pulling focus from the water surface (s_f = 0.5 m) to the bed (s_f = 1.0 m) narrows FOV by 5.6% at f = 50 mm and 10.2% at f = 85 mm. Worse, v4-4.1's window formula W = W_s (s/f − 1) holds only on the focused plane (s = s_f); off the focal plane, a thin lens projects W(s; s_f) = W_s s/v_f = W_s s (s_f − f)/(f s_f). Therefore, if v1 "excludes focus breathing" by locking FOV constant during v4-4.4's focus pull, it breaks thin-lens projection geometry by creating an unphysical constant-FOV pinhole + thin-lens-DoF hybrid.

Failure condition: Executing v4-4.4's surface-to-bed focus pull with a 50 mm or 85 mm lens (s_f: 0.5 m → 1.0 m): either frame boundaries breathe by 5–10% (violating v1's "excluded" claim and shifting main's window), or FOV is artificially pinned by fixing v_f = f in projection while varying v_f in blur (violating thin-lens optics).

C7 — WEAKENED

Strongest counterargument: Beyond importing the viewing CoC 0.029 mm (violating C5) and using v4-4.3's buggy b formula, C7 contains three concrete quantitative and experimental flaws:

Pixel-pitch mismatch with settled 2048^2 master resolution: C7 uses 5.9 μm (23 px blur, 3.6 px Airy disc), which corresponds to a ∼6100 px (24 MP) sensor—not the settled 2048^2 master resolution (36 mm/2048 = 17.6 μm pitch), where 0.139 mm is only 7.9 px and a 21.5 μm Airy disc is 1.2 px.

Unfair regime comparison that breaks the PRESETS framing: At h_a = 1.0 m, d_w = 0.6 m, a 35 mm full-frame lens frames main (W = 1.4755 m at s = 1.45 m). Raising the camera to h_a = 3.0 m while keeping f = 35 mm widens the framed bed by 2.4×, abandoning the main preset! Holding the main preset fixed at h_a = 3.0 m requires f ≈ 84 mm, where f/8 DoF at s_f = 3.0 m is 2.74–3.32 m—leaving the bed (3.45 m) still out of focus. Likewise, in sub (d_w = 4.0 m, apparent depth 3.0 m), 35 mm framing puts h_a ≈ 1.70 m, where f/8 surface focus (D_far ≈ 2.51 m) leaves the bed (4.70 m) strongly out of focus.

Omission of caustic cusp scale & wave curvature: While caustic cells are 100–400 px wide, caustic filaments set by the 0.00395 rad sun disc are only ∼2.4 mm (3.3 px on main), so 7.9 px defocus noticeably blurs them. Moreover, "does nothing under most conditions" mistakes wide-DoF physical optics for a software defect: when real optics keep surface and bed sharp, reproducing nature means rendering them sharp.

Failure condition: Holding the main preset window constant across camera height 1 m → 3 m (f = 35 → 84 mm at f/8), evaluating sub at f/8, or measuring pixel blur at the actual 2048^2 pitch (17.6 μm): C7's pixel counts are off by 3× and its "focus pull is invisible at h = 3 m, f/8" conclusion reverses.

C8 — STANDS

Strongest counterargument: The physical observation is completely right: a curved wave acts as a mirror with focal length f_m = R/2, placing the sun's reflected image at s_glitter = h_a + R/2—spanning from h_a ± cm on steep capillary waves all the way to ±∞ at wave inflection points (R → ∞), matching neither the water surface nor the bed. However, the hinted fix ("so a separate layer may be needed") is deeply inadequate: (1) glitter itself has no single layer distance (it spans h_a ± R/2 → ∞ with strong wave-curvature astigmatism R1 ≠ R2); (2) a glitter highlight from a 0.00395 rad sun off a curved wave only illuminates the sub-set of the camera aperture that satisfies specular reflection, producing non-uniform/partial-aperture bokeh rather than a simple thin-lens disc; and (3) Stage 3 volumetric scattering has the exact same problem (radiance integrated continuously over [h_a, h_a + d_w/n]), showing that patching a 2D post-transport PSF with ad-hoc "layers" cannot replace multi-depth or ray-based aperture integration.

Failure condition: A sun glint near a wave inflection point (R → ∞, s_glitter → ∞) overlapping a sharp capillary trough (R = −0.04 m, real image 2 cm above surface) when the camera is focused on the surface: both belong to "the glitter layer," yet one requires the infinity blur circle and the other is nearly sharp or front-focused, and both partially transmit/overlay refracted bed radiance underneath.

C9 — WRONG

Strongest counterargument: C9 assumes apparent-image astigmatism is purely an oblique-viewing effect and harmless overhead. Both assumptions fail:

Wave-induced astigmatism at strictly 0° overhead capture: Real water is rippled, not flat. By Coddington's equations at normal incidence, a surface with principal curvatures (κ1, κ2) maps bed depth d_w to two distinct apparent depths: 1/d_app,i = n/d_w − (n−1)κ_i. By Helmholtz optical reciprocity, the exact curvature condition for a wave crest to focus parallel sunlight into a caustic on the bed (κ ≥ n(n−1)d_w / n) is the exact condition that drives d_app,1 through ±∞ into a real image above the water, while κ2 (e.g., along a wave ridge where κ2 ≤ 0) keeps d_app,2 ≤ d_w/n. Thus, every caustic-forming ripple produces severe direction-dependent astigmatism directly overhead.

Oblique astigmatism already enters the "nearly overhead" wide-angle presets: In main (35 mm full-frame, W = 1.4755 m at s = 1.45 m), corner rays enter at θ2 ≈ 32° (θ1 ≈ 23.4°), yielding even on flat water a 15% tangential-vs-sagittal apparent-depth split (d_t = d_w cos^3 θ2 / (n cos^3 θ1) = 0.355 m vs. d_s = d_w cos θ2 / (n cos θ1) = 0.416 m).

Failure condition: Any overhead frame of main or sub containing sharp caustics (especially sub at d_w = 4.0 m, where κ_caustic = 0.33 m^−1, R ≤ 3 m): directly over bright caustic cusps, a representative scalar distance cannot match both principal focal axes, missing anisotropic defocus streaking and hitting singularities where 1/d_app,1 = 0.

Open Question Q1: Fate of Frozen Constants E = 2.2 and A = 2.6

Attack on Branch 1 (REINTERPRET):

What breaks: It is physically impossible for a camera's recording exposure (shutter × aperture × gain) to apply E = 2.2 exclusively to the direct term and A = 2.6 exclusively to the ambient term; a sensor integrates all incident photons blindly without tagging them as "sun" vs. "sky." Furthermore, locking output level to fixed E = 2.2, A = 2.6 breaks physical exposure coupling when Stage 4 sweeps aperture N (for DoF in v4-4.3) or shutter time (v4-4.5)—forcing an implicit, unmodeled inverse gain compensation. Worse, if the frozen composite is already display-clipped or clamped after E = 2.2, running Stage 4 PSF convolution afterward blurs clamped values instead of linear HDR radiance, destroying caustic/glitter highlight blooming.

What it contaminates downstream: It smuggles a scene/lighting contrast parameter (A/E = 2.6/2.2) into Stage 4 disguised as camera behavior, corrupting Stage 4's physical isolation. It also starves PJ9 (eye adaptation / CSF) of absolute physical luminance (cd/m^2), since unitless frozen multipliers have no link to real solar irradiance or camera exposure value (EV).

Attack on Branch 2 (REPLACE):

What breaks: Retiring E = 2.2 and A = 2.6 inside Stage 4 leaves the upstream transport output uncalibrated. Unless upstream transport outputs absolute spectral radiance in physical units (W⋅sr^−1⋅m^−2⋅nm^−1) driven by physical solar and sky irradiances, feeding dimensionless transport numbers into a physical camera equation (t_shutter ⋅ π/(4N^2) ⋅ gain) produces arbitrary exposure scale (black or clipped). Crucially, retiring A = 2.6 also changes the direct-to-ambient ratio (2.2/2.6 → 1.0), altering shadow depth and scene contrast unless A/E is explicitly re-homed upstream as a physical sky-to-sun irradiance ratio. Varying N across C7's f/2.8 → f/16 test also drops image brightness by 5 stops (32×) unless shutter/gain actively tracks exposure.

What it contaminates downstream: It breaks visual and numerical parity with the frozen composite and all pre-Stage-4 PRESET baselines, and shifts highlight saturation knees relative to PJ9's adaptation/tone mapping.

Attack on the Framing of Q1 Itself:
The binary choice ("keep both as camera exposure" vs. "replace both with camera exposure") is ill-posed because E = 2.2 and A = 2.6 are not the same kind of quantity. Writing I = E⋅(L_direct + (A/E) L_ambient) exposes the split:

Only a single global scale (e.g., E or a radiometric calibration constant K) can belong to Stage-4 exposure/recording.

The ratio A/E = 2.6/2.2 is strictly an environment / light-transport parameter (sky-to-sun balance or look-dev tweak) and can never be produced or replaced by a camera recording model. Neither branch is coherent unless A/E is separated from total exposure scaling.

Evidence that would distinguish the branches (and resolve the split):

Audit the units and raw mean values of L_direct and L_ambient in the frozen composite prior to E = 2.2, A = 2.6: compare raw L_direct/L_ambient to physical clear-sky sun/sky irradiance ratios on Earth (∼4:1–8:1 at surface). If transport already uses physical units, E and A are unphysical look-dev hacks; if transport outputs normalized dimensionless terms (⟨C⟩ = 1, ambient = 1), then A/E is substituting for missing environment sky/sun balance and E for missing radiometric gain.

Check whether the frozen composite preserves unclipped linear HDR peaks (>1.0) on caustic cusps and glitter: if clipped, Branch 1 is immediately disqualified for physical PSF imaging.

Run a physical aperture sweep (f/2.8 → f/8 at fixed shutter/ISO): physical exposure scales uniformly as 1/N^2 on direct and ambient together, never differentially.

Missed Items (Physics, Validation, and Consistency Checks)

1. Pipeline ordering inversion: v4-4.1 finite-camera perspective rays affect upstream light transport.

Why it matters physically: Moving from orthographic ("photo = window") to a finite camera at h_a = 1.0 m introduces viewing angles up to ∼32° in main. Viewing ray inclination θ(u,v) alters water-surface Snell refraction (nonlinear barrel/pincushion magnification of the bed across the flat air-water interface: x_bed = h_a tan θ2 + d_w sin θ2 / sqrt(n^2 − sin^2 θ2)), Fresnel surface reflectance R_F(θ), and volumetric water path length d_w/cos θ1. Generating v4-4.1's ray-vector field after transport cannot retroactively fix orthographic transport rays.

Where it belongs: Capture geometry (v4-4.1) must precede and drive light transport (Stages 1–3), not sit post-transport inside Stage 4.

2. Aspect-ratio mismatch between PRESETS (15:11 = 1.3636) and master resolution (2048^2 = 1:1).

Why it matters physically: sub (4.8 × 3.52 m) and contrast (2.4 × 1.76 m) are 15:11 rectangles (main 1.4755/1.0828 = 1.3627), whereas 2048^2 is square. Without explicit cropping or sensor aspect matching, pixels are anisotropic (Δx/Δy = 1.364), which distorts any pixel-space circular blur disc, Airy PSF, or OLPF into a physical ellipse and breaks MTF Nyquist symmetry.

Where it belongs: Consistency check between v4-4.1 (sensor dimensions W_s × H_s) and v4-4.5 (sampling pitch p_x, p_y).

3. Temporal shutter integration (motion blur across the 8 s seamless loop).

Why it matters physically: v4-4.5 invokes physical shutter time solely as a brightness multiplier, while v4-4.4 declares "no new time-direction physics beyond s_f(t)." In nature, a finite shutter time t_shutter integrates moving capillary waves (≥0.23 m/s, several px/frame on main), moving caustic cusps, streaking sun glitter, and moving focus s_f(t) over [t, t + t_shutter]. Instantaneous δ(t) frames with a "shutter time" scalar are physically inconsistent, and rolling vs. global shutter affects the t = 0 ↔ 8 s loop boundary.

Where it belongs: Stage 4 v4-4.4 / v4-4.5 (temporal exposure integration wrapped periodically at T = 8 s).

4. Linear HDR ordering vs. full-well saturation, shot noise, and veiling glare on extreme-DR glitter/caustics.

Why it matters physically: Caustic cusps and sun glitter span 10^3–10^5× dynamic range. Optical PSF convolution (defocus, diffraction, internal lens veiling glare) occurs in linear optical irradiance before sensor full-well clipping, photon shot noise, read noise, and ADC quantization. The design specifies neither where clipping happens nor sensor noise/dynamic range.

Where it belongs: Strict HDR linearity enforced in v4-4.3; full-well capacity, shot/read noise, and encoding in v4-4.5.

5. Wavelength-dependent diffraction PSF(λ) and Bayer CFA motivation for OLPF.

Why it matters physically: v4-4.6 notes that water dispersion creates rainbow fringes on caustics. Across 450–650 nm, Airy disc diameter (2.44 λN_w) varies by 44%; a single monochromatic λ mis-models diffracted dispersion fringes. Furthermore, v4-4.5 includes an OLPF (which physically exists to suppress color moiré on a Bayer CFA grid) without modeling CFA sub-sampling or demosaicing.

Where it belongs: Spectral PSF integration in v4-4.3 coupled to v4-4.5's spectral sensitivity and explicit CFA/sensor-grid specification.

6. Polarization coupling at air-water Brewster angles and camera optics.

Why it matters physically: Near Brewster's angle (∼53° for water), surface glitter/sky reflection is strongly s-polarized while refracted bed light is p-enriched; camera polarizers or oblique optical coatings alter the glitter-to-caustic ratio dramatically.

Where it belongs: Water Fresnel polarization in Stage 1/transport; optional polarizer/coating transmission in Stage 4 (v4-4.5/v4-4.6).

7. Skipped validation & consistency checks:

- Non-flat energy conservation: v4-4.3 checks flux conservation only "on flat regions." A gather-based PSF using spatially varying b(u,v) violates energy conservation wherever b varies (across waves and 3D edges); only scatter/splat or true ray-sampled aperture integration conserves flux globally.
- Curved-wave analytic calibration: v4-4.2 calibrates only against flat water, never testing a curved sinusoidal/cylindrical wave where 1/s_app = n/d_w − (n−1)κ goes negative/singular.
- Transport-bandlimit audit: No check on whether the frozen 3 px displacement low-pass ("Blow") already cuts off caustic spatial frequencies below Stage 4's Airy/OLPF/defocus cutoff, masking Stage 4 PSF differences in PRESET tests.

Final Verdict on the Decomposition

FALLS — To make Stage 4 physically valid, (1) move v4-4.1's finite-camera ray geometry upstream to drive transport (or evaluate a true multi-ray aperture light field through the water surface rather than a 2D post-transport PSF), (2) replace v4-4.2's scalar apparent-distance field (which breaks on wave-induced caustic astigmatism, sun glitter, 3D occlusion, and Stage 3 volume), (3) fix v4-4.3's swapped-denominator blur formula, working-f-number Airy scaling, and wave-optics defocus/energy conservation, (4) add temporal shutter integration across the 8 s loop in v4-4.4/v4-4.5, and (5) separate the upstream direct-to-ambient ratio (A/E) from Stage 4's uniform physical exposure.
