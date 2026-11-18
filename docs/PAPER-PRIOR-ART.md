# Prior-art / novelty assessment — deterministic water caustics under a Perez sky

Research conducted 2026-11-10 to 2026-11-11. Sources are flagged **opened** (full text fetched via browser.open) or **index** (search-result / abstract-level). No live interactive-site verification was needed; all sources are static publications.

## Summary

| Claim | Verdict | One-line basis |
|---|---|---|
| 1 — η-sweep, non-monotonic p99 at η≈0.25, diffuse caustics persist | **(b) partially covered** | Diffuse/overcast weakening-but-not-vanishing and diffuseness sweeps are published in ocean optics (Stramski 1986; Hieronymi 2011–2013). The non-monotonic bump, the Perez-sky bottom-field statistic, and the horizon-brightening + Snell-compression mechanism were **not found anywhere**. |
| 2 — Cox–Munk underestimates pond-scale slope variance ≈2×, caustic peaks 1.45×→2.11× | **(b) partially covered** | Fetch-dependence and a fetch-aware Cox–Munk modification are published (Woźniak 1996; JONSWAP/Elfouhaily). No pond-scale (100 m fetch) slope/glitter/caustics study, no rendering correction, and no caustic-brightness consequence were found. |
| 3 — Photon validation of a 2048² field is shot-noise-limited (MS-SSIM 0.285 @4M, 0.46 @16M despite image agreement 0.996) | **(b) partially covered — specific result apparently unclaimed** | "Millions of photons still noisy for caustics" is published (Hachisuka et al. 2008). All published caustic comparisons are image-space. No photons-per-cell field-validation analysis or log-field MS-SSIM precedent was found. |

Confidence: moderate-high for the "not found" parts — searches covered SIGGRAPH/TOG, Eurographics/CGF, Pacific Graphics, JCGT, Applied Optics, JOSA, Optics Express, Oceanologia, Ocean Science, JGR and arXiv, in English plus a Japanese batch. Absence-of-evidence is never proof; wording in the paper should be "to our knowledge, first…" rather than absolute.

## Findings

### Claim 1 — Caustics under diffuse sky, direct-fraction sweep, non-monotonicity

**Closest prior work**

1. **Stramski, D. (1986), "The effect of daylight diffuseness on the focusing of sunlight by sea surface waves," Oceanologia 24, 11–27.** Opened: https://ecp.iopan.pl/oceanologia/OC_24/OC_24_11-27.pdf
   - This is the single closest paper to sub-questions (i) and (ii). Theory of underwater light "flashes" as a function of daylight diffuseness d = diffuse/global (= 1−η in the project's notation), assuming only the direct component varies in time, plus simultaneous Baltic Sea measurements of surface diffuseness and underwater flashing light (525 nm).
   - Establishes: flash characteristics vary little for d < 0.5–0.6, then decrease **drastically** in frequency and maximum intensity for d > 0.7–0.8. The dependence it reports is monotone/threshold-like — no intermediate bump.
   - Does **not** establish: spatial bottom-plane caustic-field statistics (it is temporal flash statistics at a point), a Perez sky, spatial p99/std of a field, or any non-monotonicity.

2. **Hieronymi, M. (2013), "Monte Carlo code for the study of the dynamic light field at the wavy atmosphere–ocean interface," J. Eur. Opt. Soc. – Rapid Publ. 8, 13039, DOI 10.2971/jeos.2013.13039.** Abstract, index: https://www.mendeley.com/catalogue/8d11885e-adf5-3ac4-a6e1-35fc312d4791/
   - Monte Carlo radiative transfer at the wavy air–sea interface; characterises irradiance/reflectance/radiance below the surface under changing sun zenith angles **and percentages of diffuse sky radiation**. Focusing most effective under clear skies; variability **significantly smaller under overcast skies**. Also notes transmission appears overestimated for solar zenith angles near the horizon — a neighbour observation to the grazing-angle regime, but framed as an error, not as a focusing enhancement.
   - Companion work: Hieronymi & Macke (2012), Ocean Science 8, 103–120 (clear-sky setup uses 10% diffuse / 90% direct, isotropic diffuse skylight; relations page opened: https://os.copernicus.org/articles/8/103/2012/os-8-103-2012-relations.html ; discussion paper: https://os.copernicus.org/preprints/8/2101/2011/osd-8-2101-2011.pdf) and Hieronymi & Macke (2012), Ocean Science 8, 455–471 (peaks >10× mean at winds ≤5 m/s).
   - A Hieronymi manuscript/chapter, "Variability of Downwelling Irradiance in the Upper Ocean in Response to Surface Waves and Diffuse Sky Radiation" (ResearchGate record, index: https://www.researchgate.net/publication/259338146_Variability_of_Downwelling_Irradiance_in_the_Upper_Ocean_in_Response_to_Surface_Waves_and_Diffuse_Sky_Radiation ), states verbatim in its abstract: underwater fluctuations are strongly reduced as the direct portion decreases, "**however, even under overcast sky with 100% diffuse light input, small wave-induced irradiance variability can be observed**," and its comparison figure shows that under overcast input the foci are much weaker, more distorted, and not directly beneath the wave crest. **This pre-publishes the qualitative half of "diffuse caustics weaken but do not vanish."** Its exact venue/year could not be verified (see Could not verify) — cite the JEOS 2013 / Ocean Science papers for the peer-reviewed version of the claim.

3. **Supporting ocean-optics literature (all index unless noted):**
   - Dera, J. & Stramski, D. (1986), Oceanologia 23, 15–42 — https://old.iopan.pl/oceanologia/OC_23/OC_23_15-42.pdf — measured flashes >1.5× mean up to ~200/min at 1 m; strongest >5× mean ≈1/min; durations 10–30 ms; focusing strongest in clear water, clear skies, high sun, winds 2–5 m/s.
   - Gernez, P. & Antoine, D. (2009), J. Geophys. Res. 114, C06025, DOI 10.1029/2008JC005059 — field characterisation; fluctuations maximal at ~0.5 m waves / 1–5 m/s wind, damped under overcast.
   - Stramska, M. & Dickey, T. D. (1998), Deep-Sea Res. I 45, 1393–1410, DOI 10.1016/S0967-0637(98)00020-X — measured wave-focusing + cloud-cover variability.
   - Gege, P. & Pinnel, N. (2011), Applied Optics 50(15), 2192–2203, DOI 10.1364/AO.50.002192 — https://elib.dlr.de/69924/ — downwelling irradiance "highly variable due to focusing and defocusing of sunlight **and skylight**"; validated on three German lakes (the closest inland-water setting found anywhere in this research); typical variability 20–40% near surface. Related Gege work models direct and diffuse components as affected differently / only weakly correlated by wave focusing — a mechanism neighbour, not an η-sweep of caustic contrast.
   - Wijesekera, H. W., Pegau, W. S. & Boyd, T. J. (2005), Optics Express 13(23), 9257–9264 — wave effects on upper-ocean irradiance distribution.
   - Sabbah, S. & Shashar, N. (2006), Applied Optics — radiance/polarisation fluctuations associated with diffuse light, variability highest at 70° viewing zenith and rising with solar zenith angle.

4. **The two ingredients of the proposed mechanism are individually standard:**
   - Perez, R., Seals, R. & Michalsky, J. (1993), Solar Energy 50, 235–245 — the Perez all-weather model's coefficients explicitly encode **horizon darkening/brightening** and the near-horizon luminance gradient (via comparison paper, index: https://mdpi-res.com/d_attachment/energies/energies-13-05455/article_deploy/energies-13-05455.pdf?version=1603106351 ).
   - Snell's window compression (whole sky into a ~97.2° cone) is textbook; see e.g. the Optica spotlight on Snell's window in wavy water, Appl. Opt. 54(4), B8–B11 (2015): https://opg.optica.org/spotlight/spotlight-modal.cfm?id=301426 . No source found connects Snell compression of grazing sky light to a caustic-contrast enhancement.

**Comparison — Claim 1 vs closest prior work**

| Aspect | Project claim | Closest prior work | Covered? |
|---|---|---|---|
| Caustics/fluctuations under diffuse or overcast light exist but weaken | p99 = 1.23× mean, max 2.6× mean at η=0 | Hieronymi: residual variability at 100% diffuse; Stramski 1986; Gege & Pinnel 2011 | **Yes, qualitatively** |
| Direct-fraction / diffuseness sweep | η = 1→0 sweep of field statistics | Stramski 1986 (flash stats vs d, theory + Baltic data); Hieronymi 2013 (% diffuse sky) | **Yes, in ocean-optics form** (point/temporal stats, not spatial field stats) |
| Non-monotonic p99 rise at η≈0.25 | p99 rises slightly above trend while std falls −59% | Nothing found; Stramski reports little change then drastic (monotone) fall | **No — unclaimed** |
| Mechanism: Perez horizon brightening + Snell compression (dθw/dθi≈0.19 @80°) focusing grazing diffuse light | Stated mechanism | Horizon brightening is in Perez 1993; Snell compression is textbook; their caustic consequence not found | **No — unclaimed as a mechanism** |
| Setting: deterministic Jacobian transport, bottom plane, Perez sky | — | Ocean-optics work is Monte Carlo RT or in-situ measurement, sun+isotropic-diffuse | **No** |

**Verdict Claim 1: (b) partially covered.** Already published in substance: that wave focusing persists but weakens under diffuse/overcast illumination, and that fluctuation statistics fall with diffuseness. **Remaining unclaimed:** (1) the non-monotonic p99 behaviour near η≈0.25; (2) the horizon-brightening/Snell-compression explanation; (3) the specific form of the result — spatial bottom-field contrast statistics (std, p99, max) swept over η under a full Perez sky with a deterministic transport method. The paper must cite Stramski 1986 and Hieronymi prominently and must **not** claim "first to show diffuse caustics persist."

### Claim 2 — Short fetch, Cox–Munk, slope variance and caustic peaks

**Closest prior work**

1. **Cox, C. & Munk, W. (1954), J. Opt. Soc. Am. 44, 838–850** — https://userpages.umbc.edu/~martins/PHYS650/Cox%20and%20Munk%20Glint%20paper.pdf — the baseline: Gaussian/Gram–Charlier slope pdf from sun glitter; mean square slope grows linearly with wind speed; an open-ocean, large-scale statistical description (its unsuitability for high-resolution/shallow settings is noted e.g. in IOCCG material, index: http://ioccg.org/wp-content/uploads/2018/07/hedley_shallow_water.pdf ).

2. **Woźniak, S. B. (1996), "Sea surface slope distribution and foam coverage as functions of the mean height of wind waves," Oceanologia 38(3), 317–332** — http://www.iopan.gda.pl/oceanologia/OC_38/OC38_3-2_Wozniak.pdf — **the closest Cox–Munk-correction prior work.** Explicitly proposes a modified Cox–Munk expression for optical modelling, replacing wind-speed-only parameterisation with mean wave height H̄ precisely because H̄ connects to "the wind fetch D, the sea depth h, the shape of the coastal region," and argues wind-only dispersions can differ enormously. Follow-ups (Oceanologia 38(4), 1996; 39(1), 1997) build a spectral reflectance/transmittance model on it with fetch, depth and coastline as inputs. So "correct Cox–Munk for fetch" is **not** new in ocean optics — but it is parameterised by wave height, aimed at reflectance/transmittance, and never carried through to caustic peak brightness.

3. **Fetch-limited spectra are standard oceanography:** JONSWAP (Hasselmann et al. 1973) is by construction fetch-limited; Pierson–Moskowitz is the fully-developed counterpart. Elfouhaily, T. et al. (1997), JGR 102, 15781–15796, DOI 10.1029/97JC00467 (PDF, index: https://web-wp.archive.org/web/20251219021744/http://archimer.ifremer.fr/doc/00091/20226/17877.pdf ) — unified spectrum designed for diverse fetch conditions, matching Cox–Munk MSS at high wavenumber and JONSWAP at long waves. In rendering, Tessendorf-lineage ocean work routinely uses JONSWAP/TMA spectra with wind, fetch and depth parameters — fetch-limited wave *spectra* in graphics are not new either.

4. **Counterweight — do not over-claim against Cox–Munk for the ocean:** Bréon, F.-M. & Henriot, N. (2006), JGR 111, C06005, DOI 10.1029/2005JC003343 (abstract, index: https://www.mendeley.com/catalogue/a6c4bfec-e9b7-3639-8472-7f1ebf6a8d3c/ ) — Cox–Munk fits POLDER ocean glint "excellently." A later IASI revisit (arXiv:2210.05456: https://ar5iv.labs.arxiv.org/html/2210.05456 ) finds CM-compatible MSS, with CM if anything slightly **over**estimating upwind MSS at moderate wind. Claim 2's novelty is therefore a **regime** claim (small fetch, inland, caustics-relevant band), not "Cox–Munk is wrong."

5. **Review-risk note:** slope variance / MSS is cutoff-dependent (filtered MSS with an upper cutoff frequency is standard, e.g. Sensors 2022, index: https://pmc.ncbi.nlm.nih.gov/articles/PMC9571747/ ). The σ² ≈ 0.021 vs 0.0107 comparison will be challenged unless integration limits/bands are matched and stated. Young, steep short-fetch seas are documented (e.g. Beal 2000: https://www.jhuapl.edu/techdigest/Content/techdigest/pdf/V21-N01/21-01-Beal2.pdf ; Donelan et al. 2012, JGR 117, C00J23), which supports the direction of the claim.

**Verdict Claim 2: (b) partially covered.** Already published: fetch-dependence of wave/slope statistics, and at least one explicit fetch-aware Cox–Munk modification (Woźniak 1996). **Remaining unclaimed:** (1) any slope/glitter/caustics study at pond scale (fetch ≈100 m, dimensionless fetch ≈39); (2) a short-fetch Cox–Munk correction in a *rendering* context; (3) the quantitative caustic consequence — realised σ² ≈2× the Cox–Munk value with peak brightness rising ≈1.45×→≈2.11× the mean. Gege & Pinnel's German-lake measurements (Claim 1 list) are the nearest inland-water optics precedent and should be cited.

### Claim 3 — Photon-noise limits in validating a high-resolution caustic field

**Closest prior work**

1. **Hachisuka, T., Ogaki, S. & Jensen, H. W. (2008), "Progressive Photon Mapping," ACM Trans. Graph. 27(5).** Opened: https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Hachisuka08.pdf — the strongest published anchor: photon mapping "becomes very costly for scenes dominated by caustics… **to avoid noise it is necessary to use a large number of photons**," and even with **20 million photons** "this is not enough to capture details in the illumination and the rendered image exhibits low-frequency noise." This is about rendered images and memory limits, not field-vs-field validation.
2. **Frisvad, J. R. et al. (2014), "Photon Differential Splatting for Rendering Caustics," Computer Graphics Forum 33(6), 252–263, DOI 10.1111/cgf.12347** — https://backend.orbit.dtu.dk/ws/portalfiles/portal/115332999/diffsplat.pdf — the closest "caustics-only metric" precedent: RMSE and SSIM against a path-traced reference, including a caustic-illumination-only comparison. Still **image-space**, at image resolutions, not a 2048² field at ~1 photon/cell.
3. **Huang, J. et al. (2024), "Online Photon Guiding with 3D Gaussians for Caustics Rendering," arXiv:2403.03641** — https://arxiv.org/html/2403.03641v2 — reports MSE and (1−SSIM) on rendered images; tellingly, its *reference* images are SPPM renders using **64K iterations × 1M photons (64 billion photons)** because noiseless caustics are so hard to obtain — a vivid published illustration of the budget problem, again image-space.
4. The density-estimation bias/variance trade-off (few photons → noise, large kernel → blur) is textbook photon-mapping material (Jensen's lineage; Hachisuka & Jensen 2009 SPPM).

**Verdict Claim 3: (b) partially covered, specific result apparently unclaimed.** Already published: that caustics demand enormous photon budgets and remain noisy at tens of millions of photons. **Not found anywhere:** (1) a photons-per-cell formulation of the validation limit for a 2048² caustic *field*; (2) a field-level comparison metric (log-field MS-SSIM) as opposed to image-space metrics; (3) the documented divergence pattern — means agreeing to 1e-12 and final images at MS-SSIM 0.996 / 0 px peak shift while the field metric sits at 0.285 (4M) / 0.46 (16M). This is the cleanest methodological novelty of the three claims, but also the one most likely to be read by graphics reviewers as "known photon noise" unless framed precisely as a *validation-methodology* result (field metrics vs image metrics disagree; image agreement does not certify field agreement).

### Closest works A — deterministic (non-photon) caustics transport via Jacobian / specular-manifold methods

| Work | What it establishes | Distance from the project's method |
|---|---|---|
| Heckbert & Hanrahan (1984), SIGGRAPH '84; Watt (1990), SIGGRAPH '90 | Beam tracing of light through water — early deterministic transport lineage | Polygonal beams, no Jacobian field statistics |
| Mitchell & Hanrahan (1992), SIGGRAPH '92 | Analytic caustics from curved reflectors via Fermat/Newton | Implicit surfaces only, costly 3D search |
| Walter, B., Zhao, S., Holzschuch, N. & Bala, K. (2009), ACM TOG 28(3), Art. 92 | **Fully deterministic Newton–Raphson search for all contributory paths** through one refraction, Jacobian pseudoinverse, per-triangle enumeration | Closest in spirit to "all-preimage" completeness — but point light / volumetric single scattering, not a wave-surface bottom field under a sky |
| Jakob & Marschner (2012); Hanika, Droske & Fascione (2015), CGF 34(4) — https://jo.dreggn.org/home/2015_mnee.pdf | Manifold exploration / Manifold Next Event Estimation: deterministic manifold walk to a valid caustic connection per shading point, no photon cache | Monte Carlo overall; typically one solution per connection; point/directional lights |
| Zeltner, T., Georgiev, I. & Jakob, W. (2020), ACM TOG 39(4), Art. 149, DOI 10.1145/3386569.3392408 | Specular Manifold Sampling — generalises MNEE with stochastic initialisation + unbiased weights for high-frequency caustics and glints | Sampling method, 1–2 specular vertices; not an exhaustive field transport |
| Loubet, G. et al. (2020), ACM TOG 39(6), Art. 239, DOI 10.1145/3414685.3417811 | Analytic slope-space integrals for specular NEE; Fig. 1(b) is literally underwater caustics from sunlight | Exactly one specular interaction |
| Yan, L.-Q. et al. (2014), SIGGRAPH 2014 — http://sites.cs.ucsb.edu/~lingqi/publications/paper_glints.pdf | Deterministic glint rendering; Jacobian determinant sets highlight intensity | Glints (directional caustics), normal-mapped surfaces |
| **Fan, Z. et al. (2024), "Specular Polynomials," ACM TOG 43(4), Art. 126, 13 pp., DOI 10.1145/3658132** — https://arxiv.org/pdf/2405.13409 | Newton-free, "**completely deterministic**" polynomial root-finding giving the **complete set of admissible specular paths** for one bounce; applied to glints and caustics; "free from stochastic sampling, thus avoiding noise" | The strongest "all solutions, deterministically" precedent — but a per-connection path solver in a Monte Carlo framework, not a surface-to-receiver Jacobian field transport. (Note: the author PDF header shows Art. 1 / DOI 10.1145/3618360; Crossref and the author's site give Art. 126 / DOI 10.1145/3658132 — cite the latter, flag the discrepancy if quoting the PDF.) |
| Fan, Z. et al. (2025), "Bernstein Bounds for Caustics," ACM TOG 44(4), arXiv:2504.19163 — https://arxiv.org/pdf/2504.19163 | Deterministic admissible-set search is slow/high-variance at low sample rates, so triangles are sampled by Bernstein irradiance bounds | Sampling on top of deterministic bounds |
| Yuksel, C. & Keyser, J., "Fast Real-time Caustics from Height Fields" — https://www.cemyuksel.com/research/heightfield_caustics/fastcaustics.pdf | Receiver-backward two-pass deterministic-ish algorithm, height-field water + planar receiver — geometrically the closest scene setup | Real-time approximation, sun light, no physical statistics. Lineage also includes Guardado & Sánchez-Crespo (2004), GPU Gems ch. 2 (explicitly aesthetics-driven: https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-2-rendering-water-caustics ), Shah et al. (2007) TVCG caustics mapping, Wyman & Davis (2006), Stam (1996) "Random Caustics" (https://www.dgp.toronto.edu/public_user/stam/INRIA/caustics.html ) |

Positioning: "deterministic + all solutions" is **not** by itself novel post-2024 (Fan et al.; Walter et al. 2009). The defensible novelty is the combination: exhaustive preimage/Jacobian transport specialised to a wave height-field → planar receiver, producing a full field whose *statistics* are the object of study, under area (sky) illumination.

### Closest works B — physically-based caustics under sky (not just sun) illumination in rendering

1. **Dammertz, H. & Hanika, J. (2009), "Plane Sampling for Light Paths from the Environment Map," Journal of Graphics Tools.** Opened: https://jo.dreggn.org/home/2009_plane_sampling.pdf — **the key graphics precedent.** Starts light paths on an emitting environment map "in the context of Monte Carlo global illumination with image-based lighting or **physical sky models**… we can efficiently render **unbiased caustics from these kinds of lights**," and makes photon mapping/MIS/instant radiosity correct with environment maps. Directly relevant caveat: "a mostly diffuse environment such as a cloudy sky is best sampled by a path tracer" — i.e. the diffuse-sky caustic regime was identified as the awkward limit, but its field statistics were never studied.
2. **Tandianus, B., Johan, H. & Seah, H. S. (2012), Entertainment Computing 3(4), 129–141, DOI 10.1016/j.entcom.2011.11.004** — real-time approximate caustics under environment illumination, approximated as a set of important directional lights with precomputed per-direction patterns; "visually similar," no physical statistics.
3. **Wand, M. & Straßer, W. (2003), Eurographics 2003, CGF 22(3)** — http://www.staff.uni-mainz.de/wandm/papers/eg03.pdf — specular surfaces as pinhole cameras projecting the environmental light onto receivers; dynamic environments.
4. **Premož, S. & Ashikhmin, M. (2000/2001), Pacific Graphics 2000 / CGF 20(4), 189–199** — physically-based light transport for natural waters (sun+sky context) with particle/photon-style caustics.
5. The ocean-optics literature (Hieronymi, Stramski, Gege — Claim 1) is the real "sky caustics" science literature, but as underwater irradiance science, not rendering, and none of it uses a Perez sky.

**Japanese work (asked for explicitly):** no Japanese-language paper on diffuse/Perez-sky caustics statistics was found. Major relevant work by Japanese authors is in English venues: Iwasaki, K., Dobashi, Y. & Nishita, T., Pacific Graphics 2001 (http://nishitalab.org/user/nis/cdrom/pg/pg2001_iwa.pdf ) and Iwasaki, Dobashi & Nishita (2003), CGF 22, 601–609 (https://onlinelibrary.wiley.com/doi/epdf/10.1111/1467-8659.t01-2-00708 ) — fast/real-time water caustics, sun-based. One Japanese-venue item: Guerin, K. & Izumi, M. (2017), IPSJ SIG Technical Report 2017-CG-166 No. 6 (https://ipsj.ixsq.nii.ac.jp/record/178414/files/IPSJ-CG17166006.pdf ) — real-time approximation for large marine environments, light-direction based, no diffuse-sky statistics.

### Publication venues — realistic assessment (simulation-internal measurements only, no photographic validation yet)

**Optics — best fit for Claims 1 and 2 as physics findings (recommended first targets):**
- **Applied Optics** — Gege & Pinnel's lake work and much of the fluctuation literature sits here / in sibling Optica venues; simulation + theory papers are in scope without field data.
- **JOSA A** — wave-optics / statistical optics framing suits the η-sweep and slope-statistics results.
- **Optics Express** — Wijesekera et al. 2005 precedent; fast, simulation-friendly.
- **Journal of the European Optical Society – Rapid Publications** — Hieronymi 2013 precedent, directly on-topic readership.
- **Ocean Science** — Hieronymi & Macke precedent; modelling papers without new field data are normal there.
- JGR Oceans is **not** realistic yet — it will expect field data.

**Graphics — realistic if framed around the method (Claims A/3):**
- **JCGT (Journal of Computer Graphics Techniques)** — explicitly values practicality and reproducibility over novelty/theory, and accepts "new and extensive analysis of previous work" (scope, index: https://www.jcgt.org/write.html ). A good home for the deterministic-transport technique + validation-methodology (Claim 3) paper.
- **Eurographics Symposium on Rendering (EGSR, proc. in Computer Graphics Forum)** and **Pacific Graphics (proc. in CGF)** — plausible for the transport method if generality beyond water is argued; Pacific Graphics also matches the Japanese water-rendering lineage.
- **Computers & Graphics / The Visual Computer** — realistic journal fallbacks for a scene-specialised method.
- **SIGGRAPH / ACM TOG: unlikely in current form.** Reviewers there will press on (i) generality of the transport method vs Fan et al. 2024 / Walter et al. 2009, and (ii) absence of any real-world validation. A TOG attempt becomes credible only with either a broader method contribution or photographic/experimental validation (e.g. a tank/pool experiment under controllable diffuse light — which would simultaneously be the strongest possible evidence for Claim 1's η≈0.25 bump).

Suggested split: one optics paper (Claims 1+2, η-sweep + short-fetch, citing Stramski/Hieronymi/Woźniak as the lineage) and one graphics techniques paper (deterministic transport + Claim 3 validation methodology, for JCGT/EGSR).

## Could not verify

- **Exact venue/year of Hieronymi's "Variability of Downwelling Irradiance… Diffuse Sky Radiation"** — found only as a ResearchGate-hosted manuscript (HZG Institute of Coastal Research); its content is quoted above but it should not be cited with a venue that was not verified. The peer-reviewed anchors for the same claims are Hieronymi 2013 (JEOS-RP) and Hieronymi & Macke 2012 (Ocean Science).
- **Fan et al. 2024 article metadata discrepancy** — author PDF header: TOG 43(4) Art. 1, DOI 10.1145/3618360; Crossref/author site: Art. 126, DOI 10.1145/3658132. Both are reported above; verify against the ACM DL before citing.
- **Negative results are search-bounded:** no precedent was *found* for the η≈0.25 non-monotonicity, the pond-scale caustics consequence, or field-level caustic metrics — within English + one Japanese search batch across the venues listed. This supports "to our knowledge, first" phrasing, not an absolute priority claim.
- The mission's own numbers (σ² values, MS-SSIM values, dθw/dθi≈0.19) are the project's internal results; this research assessed prior art only and did not independently verify them.

## Sources

Opened (full text), read 2026-11-10:
- Stramski 1986 — https://ecp.iopan.pl/oceanologia/OC_24/OC_24_11-27.pdf
- Hachisuka et al. 2008 — https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Hachisuka08.pdf
- Dammertz & Hanika 2009 — https://jo.dreggn.org/home/2009_plane_sampling.pdf
- Hieronymi & Macke 2012 relations — https://os.copernicus.org/articles/8/103/2012/os-8-103-2012-relations.html

Index / abstract-level, read 2026-11-10 – 2026-11-11 (selection; full list in notes/key-sources.md):
- Hieronymi 2013 abstract — https://www.mendeley.com/catalogue/8d11885e-adf5-3ac4-a6e1-35fc312d4791/
- Hieronymi diffuse-sky manuscript — https://www.researchgate.net/publication/259338146_Variability_of_Downwelling_Irradiance_in_the_Upper_Ocean_in_Response_to_Surface_Waves_and_Diffuse_Sky_Radiation
- Hieronymi & Macke discussion paper — https://os.copernicus.org/preprints/8/2101/2011/osd-8-2101-2011.pdf
- Dera & Stramski 1986 — https://old.iopan.pl/oceanologia/OC_23/OC_23_15-42.pdf
- Gege & Pinnel 2011 — https://elib.dlr.de/69924/
- Woźniak 1996 — http://www.iopan.gda.pl/oceanologia/OC_38/OC38_3-2_Wozniak.pdf
- Cox & Munk 1954 — https://userpages.umbc.edu/~martins/PHYS650/Cox%20and%20Munk%20Glint%20paper.pdf
- Bréon & Henriot 2006 abstract — https://www.mendeley.com/catalogue/a6c4bfec-e9b7-3639-8472-7f1ebf6a8d3c/
- Elfouhaily et al. 1997 — https://web-wp.archive.org/web/20251219021744/http://archimer.ifremer.fr/doc/00091/20226/17877.pdf
- Frisvad et al. 2014 — https://backend.orbit.dtu.dk/ws/portalfiles/portal/115332999/diffsplat.pdf
- Huang et al. 2024 — https://arxiv.org/html/2403.03641v2
- Fan et al. 2024 — https://arxiv.org/pdf/2405.13409 ; PDF — https://s.bingheai.cn/pdf/2405.13409
- Fan et al. 2025 — https://arxiv.org/pdf/2504.19163
- Hanika et al. 2015 (MNEE) — https://jo.dreggn.org/home/2015_mnee.pdf
- Yan et al. 2014 — http://sites.cs.ucsb.edu/~lingqi/publications/paper_glints.pdf
- Yuksel & Keyser — https://www.cemyuksel.com/research/heightfield_caustics/fastcaustics.pdf
- Wand & Straßer 2003 — http://www.staff.uni-mainz.de/wandm/papers/eg03.pdf
- Perez horizon-brightening summary — https://mdpi-res.com/d_attachment/energies/energies-13-05455/article_deploy/energies-13-05455.pdf?version=1603106351
- Snell's window spotlight — https://opg.optica.org/spotlight/spotlight-modal.cfm?id=301426
- Iwasaki et al. 2001 — http://nishitalab.org/user/nis/cdrom/pg/pg2001_iwa.pdf ; 2003 — https://onlinelibrary.wiley.com/doi/epdf/10.1111/1467-8659.t01-2-00708
- Guerin & Izumi 2017 — https://ipsj.ixsq.nii.ac.jp/record/178414/files/IPSJ-CG17166006.pdf
- JCGT scope — https://www.jcgt.org/write.html
