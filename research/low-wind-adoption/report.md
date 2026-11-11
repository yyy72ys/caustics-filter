# Adopting the low-wind capillary treatment with its basis stated: what can honestly be claimed, and on what

**Question (user decision 2026-11-12).** Below the ECKV α_m fit domain (u\* < c_m/e ≈ 0.085 m/s; pipeline U ≲ 2.7 m/s), the pipeline floors α_m per cell so the capillary part carries σ²_h = 0.003, provenance `low_wind_prior`. The v3.9 rally (3/3) ruled those cells UNRESOLVED, never PASS. The user now wants the matter *closed before v4* — not by claiming validation, but by adopting the treatment as a standing rule with its basis stated (「根拠つき採用」), decided in a fresh three-way rally, with the basis visible in the project's GitHub docs. This report surveys the published low-wind treatments, identifies the candidate basis with the least invented content, drafts the repo status wording, and writes the hardest counterarguments for the rally prompt.

**Bottom line.** The strongest honest basis is the one already shipped: the σ²_h = 0.003 prior (candidate A), adopted not as a point value but as the **centre of a swept interval** {0.0015, 0.003, 0.006}. It is the only candidate whose level is a measured number; every alternative either invents more constants (viscous-transition forms), imports a level 27× lower that the bracket evidence contradicts (Apel), cannot be implemented from open sources (Kudryavtsev 2013), or replaces the whole spectrum to fix a boundary problem (Plant 2002). What adoption may claim is narrow and must be said exactly: *the zero band is the one option measurements exclude; the adopted magnitude is the only measured wind-independent slope residual, used as a magnitude anchor with its mechanism explicitly not claimed; the resulting totals sit inside every published low-wind MSS envelope.* What it may never claim: that 0.003 is a measured pond capillary level, that any low-wind gate passes, or that the cells are validated.

All web facts cited with URL and read date 2026-11-12 unless flagged as carried from the 2026-11-11 research (`low-wind-capillary-floor-20261111/report.md`), whose citations are not re-verified here. Quantities computed from project files are labelled *derived*.

---

## 1. Scope discipline: what "adoption" changes and what it must not touch

This section is the load-bearing distinction of the whole proposal, stated before any literature, because the rally's first attack will be that adoption is the v3.9-rejected (a)+(d) returning in a new costume.

**What (a)+(d) proposed (rejected 3/3, 2026-11-11):** change the spectrum (install the floor *as physics*) *and* re-scope gate C1 so the failing check disappears. Two objects moved; the visible outcome would have moved green-ward.

**What adoption proposes now:**

| Object | Under (a)+(d) | Under ADOPTED_WITH_BASIS |
|---|---|---|
| Spectrum numbers | Floor installed as a physical claim | **Unchanged** — the v15 default already computes max(raw, floor) below the crossing (REPORT-V15-MERGE §1); adoption moves no number |
| Gate C1 at U ≤ 3 | Re-scoped away | **Unchanged** — still computed, still printed, still *not PASS* in every report |
| Two-column reporting (v41 condition) | — | **Unchanged** |
| Provenance label `low_wind_prior` on every affected state | — | **Unchanged** — travels with outputs |
| Cell status string | UNRESOLVED → effectively accepted | UNRESOLVED → **ADOPTED_WITH_BASIS**, a third status distinct from both PASS and UNRESOLVED |
| Repo docs | — | Basis, limits, and supersession condition written where a reader of the code meets them |

The only thing that changes is the *governance status of a disclosed convention*: from "open question that blocks closure" to "recorded decision with a stated basis". This is precisely the move GPT's v3.9 digest left open — re-scoping is defensible "as an explicit change of claim scope" — and it follows the precedent the rally itself set in the B-envelope round (v41): a status/gate change was ruled *principled correction, not goalpost-moving* when (i) the prior position's provenance was shown to be weaker than it looked, (ii) the replacement was literature-anchored, and (iii) reporting stayed two-column and uncollapsed. §4 states the admission criteria under which adoption is legitimate at all, so the rally can test this case against them — and so the same move cannot later be used to launder genuinely failed items (rain 3.2) into "adopted".

**One-sentence test the rally should apply:** name one number, one gate outcome, or one printed status that becomes *greener* under adoption than it is today. There is none — ADOPTED_WITH_BASIS is not a colour on the PASS axis at all; it sits on the separate axis v41 created (conformance/decision vs observation).

## 2. Survey: published low-wind treatments and measurements

### 2.1 The problem, restated in one paragraph (established; from the 2026-11-11 research)

ECKV eq. (44) is a two-regime logarithmic fit of α_m to Jähne & Riemer (1990) tank and Hara et al. (1994) field curvature levels; its lower branch crosses zero at u\* = c_m/e ≈ 0.085 m/s (*derived* pipeline crossing U ≈ 2.72 m/s) and prescribes nothing below. Gao et al. (2025) document the negative-value failure explicitly and note the low-wind fit data were "close to the noise level of the instrument". ECKV's own Cox–Munk conformance claim covers only 1 < U₁₀ < 13 m/s. Any spectrum used at U = 1.5–3 therefore contains *somebody's* completion of ECKV below its domain; the only question is whose, and how honestly labelled.

### 2.2 The Cox–Munk intercept and its three independent re-readings

The candidate anchor: σ²_cross = **0.003** + 1.92×10⁻³U (±0.002); the whole U = 0 intercept sits in the crosswind component; total 0.003 + 5.12×10⁻³U (±0.004), data spanning 1–14 m/s per the original abstract (Ocean Optics Book Cox–Munk page, *verified live* 2026-11-11). Mobley's interpretation: "a residual amount of slope that is not attributable to the local wind", "often ignored when modeling".

Three re-readings of low-wind MSS, which matter because they are *not* the same fit:

- **Wu's logarithmic reanalysis of the Cox–Munk data** (in ECKV §5.1, eq. 27): σ² = 0.0135 at U = 1.5, 0.0218 at U = 3 (*derived* from the published form). **Same underlying data** — GPT's v3.9 attack stands: agreement with Wu is calibration consistency, not independent confirmation.
- **Hu et al.'s CALIPSO lidar law**, s² = 0.0146·√U for U < 7 (via Gao et al. 2025, eq. 4): 0.0179 at U = 1.5, 0.0253 at U = 3. A *different instrument on a different platform* (spaceborne lidar vs aerial glitter photography), attributed to lidar sensitivity to finer capillary structure; joins Cox–Munk near U ≈ 7. This is the closest thing to an independent bracket at our winds, and it sits *above* the linear law at both.
- **Bréon & Henriot (2006)**, POLDER satellite glitter vs NSCAT winds (*verified* 2026-11-12, abstract): the Cox–Munk slope model "permits an excellent fit of the observation… the modeled mean square slope is in near perfect agreement with that derived from the satellite data", with wind retrievable to < 1 m/s error. Satellite-scale confirmation that CM statistics — intercept included in the component form they fitted — describe real glitter at low-to-moderate winds; not a band-resolved capillary measurement.

**Filtered-MSS residuals (Hwang).** Hwang's fits s²(k_s) = b₁(k_s)·U + b₀(k_s) carry a measured nonzero intercept b₀ at finite wavenumber cutoffs, and his spectrum integrated only to k = 316 rad/m already reproduces clean-water Cox–Munk (DTIC ADA442354; from the 2026-11-11 research). The residual is therefore not an artefact of integrating to infinite k: a wind-independent slope component exists *inside* our band of interest, at ocean scale. (Its partition between our B_l and B_h is not published — a limit, stated.)

### 2.3 Field roughness at near-zero forcing: Walsh et al. (2005) — evidence both ways

The single most on-point field study found in this round (SOWEX, NASA Scanning Radar Altimeter + optical imaging; *verified live* 2026-11-12, full text):

- When the gravity-capillary waves disappeared in the absence of wind forcing, the surface **maintained a minimum total mss of ≈ 0.0015** — total MSS did not vanish even with no local wind.
- **But its composition is the opposite of what a naive reading wants**: ≈ ⅔ of the 0.0015 minimum was contributed by the 1–10 m wavelength region and ≈ ⅓ by the dominant waves. The 0.0015 is a *gravity-scale persistence* floor, not a capillary floor. **Using Walsh's 0.0015 as a σ²_h anchor would double-count** — those scales already live in our B_l. This trap is named here because it is the most likely "better anchor" a reviewer (or a future us) reaches for.
- Donelan & Pierson thresholds, as reported there: at SOWEX water temperature, ≈ **2 m/s (10-m height) to generate 50-cm waves, rising to almost 6 m/s for 1-cm waves**. Viscous decay times (Kinsman, cited there): a 1.7-cm wave decays in ≈ 4 s unforced — the capillary band is forcing-following, with no persistence of its own.
- Near the shore of an Alaskan lake (sheltered, surfactants suspected): **small-scale mss ≈ 10⁻⁶** — a measured real-water case where the small-scale slope variance sat ~3000× below 0.003. Any honest basis must carry this as a stated limit: film-covered or sheltered water can be far smoother than the floor. (Our scenes are defined as clean water under sustained wind — the scene definition, not the floor, is what excludes that case; the rally should test whether that exclusion is legitimate.)

### 2.4 Tank evidence: inception has no minimum; effective growth has a threshold

Two tank results that appear to conflict and in fact bracket different claims:

- **Wu (1978, Tellus 30, 93–96; *verified live* 2026-11-12)**: "There is no minimum wind velocity found for the inception" — infinitesimal waves were generated even below U = 1.6 m/s; the most effective growth starts at U ≈ 1.6 m/s, at a friction velocity close to Miles' analytic minimum for wave growth, **u\* ≈ 4.4 cm/s**. Fetch at the test section 5.5 m.
- **Donelan & Plant (2009, "A threshold for wind-wave growth"; Plant's project summary *verified live* 2026-11-12)**: in a recirculating wind-wave tank with slowly ramped wind, gravity-capillary waves **do not begin to grow until wind input overcomes viscous damping**; the threshold is sharp in *friction velocity* (invariant to water temperature, current, and fetch changes of ×2) while the wind-speed threshold varies. NSCAT ocean data show the same threshold character.

Reconciliation (*inference*, stated as such): Wu's inception is a linear-instability statement (infinitesimal ripples exist at any wind); Donelan & Plant's threshold is about growth to a *sustained, measurable equilibrium band*. The prior claims neither inception nor a measured equilibrium — it claims a magnitude for a band the renderer must populate under sustained wind. The threshold result is the strongest physics *against* a high floor at U = 1.5 (a sustained capillary equilibrium may genuinely be very small there), and Wu + Hara are the strongest *for* "not identically zero". The floor at 0.003 sits between those: defensible as a magnitude anchor, indefensible as a measured equilibrium. The wording in §5 says exactly this.

**Jähne-group tank slope imaging (Aeolotron).** The modern descendant of the Jähne–Riemer measurements behind α_m: imaging slope gauges (Kiefhaber et al. 2014 CISG) measuring both slope components at 0.22 mm²/pixel, 1500 fps; **MSS detection limit ≈ 0.0003, uncertainty < 10 % for MSS > 0.002** (Engel et al. 2025, EGUsphere, *verified live* 2026-11-12, citing the clean-freshwater reference run Kunz 2017). Two uses: (i) the instrument class that *could* make the decisive measurement exists and resolves an order of magnitude below the floor; (ii) the same facility shows surface films at U < 6 m/s suppress MSS by **1–2 (locally up to 3) orders of magnitude** versus clean freshwater — the quantitative form of the §2.3 limit. Kunz 2017's clean-water MSS(U) values themselves could not be extracted in this round (thesis not openly readable) — see Could-not-verify.

### 2.5 Hara et al. (1994) — the field-vs-fit direction

Field SLSG + radar saturation compared against the Jähne–Riemer tank fit behind α_m: "reasonable agreement… but our data indicate a **higher degree of saturation at very low wind speeds**" (JGR 99, 12593; abstract *verified* 2026-11-11). The only direct statement in the literature that the fit behind α_m *under-reads the field* at exactly our winds. Qualitative (an abstract sentence; Gemini's v3.9 caveat — tidal-estuary shear turbulence a stagnant pond lacks — is recorded and fair), but it fixes the *direction* in which the fit is known to err, which is what an adoption basis needs: the floor corrects in the measured direction of error, not against it.

### 2.6 Viscous-transition forms: what exists, what is invented

- **Gemini's v3.9 counter-proposal**: α_m(u\*) = α_m0 · ½(1 + tanh((u\* − u\*_c)/Δu\*)), u\*_c ≈ 0.04 m/s, labelled an ambient/thermal noise roll-off. Physically shaped (a real capillary collapse exists near u\* ~ 0.04–0.08, far below the fit crossing — the rally's own note), but: α_m0 and Δu\* are **two further unmeasured constants**, and nothing in the literature fixes the level this form delivers at U = 1.5. It replaces one borrowed number with a borrowed function plus two borrowed constants.
- **Donelan & Pierson / Donelan & Plant threshold physics** (§2.4): real, measured, and *against* a smooth floor — it predicts a near-absence of sustained short waves below threshold, with hysteresis. As a *form* for adoption it would argue for σ²_h near zero at U = 1.5 rising sharply later — i.e., it is a counter-model, not a completion of ECKV, and adopting it would contradict Hara's field excess and every measured total at U = 1.5 (the gravity part alone in our cells is 0.0103–0.0150; a near-zero capillary part is compatible with totals — it is the *partition* claim that would change).
- **Kudryavtsev balance model**: the physically richest completion (parasitic-capillary source keeps the band from vanishing *in structure*), but *derived* in the 2026-11-11 research from its verified open terms: at U = 1.5, β_ν < 0 above k ≈ 50 rad/m and β_ν(k_g) < 0 starves the parasitic source — its wind-fed reference solution also delivers no capillary band at U = 1.5; and the current (2013) α/n functions are paywalled, the open (2003) ones superseded by their own authors. Not adoptable today; correctly named as a supersession route.

### 2.7 Plant (2002) and the "Plant limit"

Plant's unified directional spectrum (JGR 2002) is a composite built to reproduce measured spectra — the same Jähne/Hara anchor class as ECKV — plus Plant's (1982) momentum-flux bound on total MSS, mss_max = 0.08 ± 0.04 (upper value 0.12; *verified* 2026-11-12 via Davis et al. 2023's use of it). No published low-wind floor or below-domain completion specific to Plant's spectrum was found in this round (*index-level* — see Could-not-verify): it inherits the same data sparsity below u\* ≈ 0.1 m/s, and swapping the project's whole spectrum to solve a boundary labelling problem would import far more change than the problem warrants. Plant's role in an honest basis is as an *envelope* (totals far below his limit) and as an independent reminder that the field's canonical spectra are all anchored to the same sparse low-wind data.

### 2.8 Apel-form α_m(U) and operational practice

- **Apel α_m(U) = 10^(−4.95 + 3.45(1 − exp(−U/4.7)))** — the repair Gao et al. (2025) actually adopted for ocean-colour glint work (*verified* 2026-11-11): positive at all U, a single published function, zero invented constants *of ours*. But its level: α_m = 9.8×10⁻⁵ at U = 1.5 → σ²_h ≈ 0.00011, **27× below the residual anchor**, and at U = 3 it falls *below* the ECKV log value it replaces (floored totals 0.58× Cox–Munk vs status-quo 0.62×; *derived*, 2026-11-11 table). Adopting Apel = adopting "positive but nearly zero" — the completion the bracket evidence (Hu above CM; Hara above the fit) least supports. It also mixes definitions: Apel's spectrum peaks its capillary part at 750 rad/m, ECKV's B_h at k_m = 367 rad/m.
- **Practice**: operational glint/reflectance codes apply Cox–Munk statistics at near-zero wind without zeroing anything — Bréon–Henriot (2006) statistics are exercised at U = 0.2 m/s in the reference `coxmunk` implementation (*verified* 2026-11-12); SeaDAS Rayleigh LUTs are tabulated from W = 0 (*verified* 2026-11-12, Ocean Optics Book). **Honesty caveat**: the isotropic SeaDAS form uses mss = 0.00512W with *no* intercept — practice is mixed on the intercept itself. What practice uniformly does *not* do is what our original code did: let a fitted component go identically zero inside the working wind range while continuing to ship the product. The deviation was the zero, not the residual.

## 3. Candidate bases compared

Invented content is counted as: numbers or forms that appear in no measurement or published model and exist only because we chose them.

| Candidate | Adopted claim (exactly) | Explicitly NOT claimed | Invented content | Fatal weakness |
|---|---|---|---|---|
| **A. σ²_h = 0.003 prior, swept as an interval** (status quo numbers) | Below the fit domain the capillary band is assigned the magnitude of the only measured wind-independent slope residual (CM crosswind intercept), with ECKV's B_h shape; totals stay inside all published low-wind MSS envelopes | That 0.003 is a measured *pond capillary* level; any mechanism; any gate pass; validation | **One borrowed scalar** (its placement in B_h) | Mechanism mismatch is real (ocean residual, plausibly swell-borne; ponds have no swell); anchor uncertainty ±0.002 on 0.003 |
| B. Viscous-transition α_m(u\*) (tanh roll-off) | The band collapses smoothly across u\*_c ≈ 0.04 with a stated roll-off shape | That u\*_c, Δu\*, α_m0 are measured | Form + **two constants** + still a level choice | More invented content than A for a better story; level at U = 1.5 unconstrained by any number |
| C. Apel α_m(U) | A published, positive-everywhere completion, adopted wholesale as α_m | That its level matches pond reality | **Zero of ours** — but imports Apel's peak definition (750 vs 367 rad/m) | Level 27× below the residual scale; contradicts the direction of Hara and the Hu bracket; a sign fix, not a level |
| D. Kudryavtsev (2013) balance | The physically derived completion | — | None in principle | Unimplementable from open sources (paywalled, 2003 form superseded); reference solution also vanishes at U = 1.5 (*derived*) |
| E. Plant (2002) composite | Replace the spectrum with Plant's | — | None in principle | No published below-domain completion found; whole-spectrum swap to fix a boundary label |
| F. Zero / flat below domain (v14 fallback, already user-rejected) | Nothing below domain | — | Zero | Positively excluded: Wu inception, Hara excess, Hwang b₀, every measured total's partition |

**Least invented content: A**, and it is not close between A and B (the only other candidate that engages the level at all): A borrows one measured scalar; B borrows a function and two constants to reach an unmeasured level. C borrows nothing but answers a different question (sign, not level) with a level the evidence brackets against.

**Refinement that further reduces A's claim — adopt the interval, not the point.** The v13 design already mandates the sensitivity sweep σ²_h ∈ {0, 0.0015, 0.003, 0.006}. Adoption should name 0.003 the *centre* of the adopted interval [0.0015, 0.006] (factor-2 either way, matching the anchor's own ±0.002 uncertainty at the same order), with the sweep shipped in every report that uses low-wind cells. The adopted object is then "a bracketed magnitude with a stated centre", which is a strictly smaller claim than a point prior — and it is already the project's practice, so adoption again changes no computation.

## 4. Admission criteria for adoption as a governance act (so this cannot launder failures)

The rally should be asked to rule on these criteria explicitly, and on whether this case meets them. They are what separates this closure from a general amnesty:

1. **The null option is positively excluded by measurement.** (Zero band: excluded by Wu 1978 inception, Hara 1994 field excess, Hwang's b₀ residuals, and the partition arithmetic — a zero capillary part is consistent with no measured total's *shape* evidence. Contrast rain 3.2: its derivation was positively *wrong*, not merely unmeasured — adoption is unavailable there.)
2. **The adopted level is bracketed by at least two constraints that are not the same data.** (CM linear law and Hu lidar law — different instruments, different decades; Wu counts only as consistency, per GPT's v3.9 point, and is not used as a bracket leg.)
3. **No gate outcome, number, or printed status improves.** (Adoption is a claim-scope change only; C1 at U ≤ 3 remains not-PASS forever under this status.)
4. **A concrete, doable supersession measurement exists and is named.** (A slope spectrum at U = 1–3 m/s resolving k > 200 rad/m — the Aeolotron CISG class resolves MSS to 3×10⁻⁴; or a faithful Kudryavtsev-2013 implementation if its appendix is obtained. On either arriving, ADOPTED_WITH_BASIS expires automatically and the cells return to measurement judgement.)

If the rally finds any criterion unmet, the honest outcome is to keep UNRESOLVED — this report's recommendation is conditional on all four.

## 5. Draft status wording for the repo docs (verbatim candidate)

> **Status: ADOPTED_WITH_BASIS (not validated; not a gate pass).** For winds below the Elfouhaily et al. (1997) α_m fit domain (u\* < c_m/e ≈ 0.085 m/s; U ≲ 2.7 m/s in this pipeline), where the published logarithmic fit crosses zero and prescribes nothing, the capillary spectrum B_h is not set to zero: α_m is floored per cell so that the capillary slope variance equals σ²_h = 0.003 (α_m ≈ 0.0022–0.0027), retaining the ECKV B_h shape, and results are reported over the adopted interval σ²_h ∈ [0.0015, 0.006]. Basis: 0.003 is the Cox & Munk (1954) crosswind zero-wind intercept — the only measured wind-independent slope residual in the literature — used as a magnitude anchor only; the mechanism behind that residual (ocean background roughness) is explicitly **not** claimed to operate on a pond, and no published measurement resolves the capillary band of a pond at these winds. The zero alternative is the one option measurements exclude (tank inception at arbitrarily low wind, Wu 1978; field saturation above the tank fit at very low wind, Hara et al. 1994; nonzero filtered-MSS intercepts, Hwang), and the floored totals sit inside the measured low-wind MSS envelope (Cox–Munk 0.0107–0.0184, Hu et al. lidar 0.0179–0.0253, at U = 1.5–3 m/s). Known limits: the anchor carries ±0.002 fit uncertainty; the floor is wind- and fetch-independent where it acts while real levels are not; growth-threshold physics (Donelan & Plant) implies the true sustained capillary level at U = 1.5 m/s may be lower; film-covered or sheltered water can be far smoother (small-scale mss ≈ 10⁻⁶ measured). Band-share gate C1 at these winds remains **not passed** and is printed, not waived. This convention is superseded automatically by any published slope spectrum at U = 1–3 m/s resolving k > 200 rad/m, or by a faithful implementation of Kudryavtsev et al. (2013). Every state built with the floor active carries the provenance label `low_wind_prior`.

## 6. Hardest counterarguments — for the rally prompt, verbatim

Put these to the reviewers exactly as written:

1. **"This is the same goalpost-moving rejected in the v3.9 round. The floor was inadmissible as physics then; renaming its status to ADOPTED_WITH_BASIS changes no evidence. What, precisely, is different now — other than the project's desire to close the item before v4?"**
2. **"Name one number, one gate outcome, or one printed status that becomes greener under adoption. If none, this is pure bookkeeping — and bookkeeping does not need a three-way rally. If the rally is being used to lend scientific weight to an administrative act, that is worse than the original (a)+(d)."**
3. **"The adopted centre, 0.003, carries a ±0.002 fit uncertainty and is an extrapolation to U = 0 of a 1954 ocean fit whose data may start at 3 m/s (Karaev et al. 2019) rather than 1 m/s. You are adopting, as your 'basis', a number statistically consistent with zero. Why is an interval whose lower leg is 0.0015 — itself inside the noise — a basis rather than a confession?"**
4. **"Walsh et al. (2005) measured small-scale mss ≈ 10⁻⁶ on a real lake, and the Aeolotron programme shows surface films suppress MSS by 1–3 orders of magnitude below U = 6 m/s. Ponds — this project's headline scene — are precisely the waters with films, sheltering, and short fetch. Your floor may overshoot the truth for your own scenes by orders of magnitude. What in the scene definition, stated in advance and not invented for this defence, excludes film and shelter?"**
5. **"Donelan & Plant's threshold is measured physics: below a friction-velocity threshold, gravity-capillary waves do not grow; the threshold for 1-cm waves is ≈ 6 m/s wind. Viscous decay kills a 1.7-cm wave in ~4 s. At U = 1.5 m/s there may be no sustained capillary equilibrium band at all — Wu's 'infinitesimal inception' ripples are not a band carrying σ² = 0.003. Is the floor modelling a band that physics says is not there under sustained forcing?"**
6. **"Availability is not a basis. The Cox–Munk intercept is 'the only measured residual' largely because it is the only number of its kind anyone published. In what other branch of physics would borrowing an intercept across mechanisms — ocean swell-borne background roughness standing in for wind-driven pond capillaries — survive review? If the mechanism is disclaimed, the claim reduces to 'we needed a number and this was the nearest one.'"**
7. **"Governance precedent: if ADOPTED_WITH_BASIS is legitimate here, state the admission criteria that this case meets and that rain (3.2, REJECTED_INVALID_DERIVATION) and the unvalidated scattering halo (3.6) do not. If you cannot state criteria that exclude them, adoption is a general amnesty mechanism and every future UNRESOLVED will queue up to use it."**
8. **"Motivation audit: the sensitivity sweep exists because the floor restores fine texture — the warp the artwork needs. The level was chosen from the MSS literature, you say, before renders. But the *decision to adopt rather than keep UNRESOLVED* arrives exactly when an open status became inconvenient (v4 waiting). Why should a reviewer not read the basis as post-hoc cover for an aesthetic requirement plus schedule pressure?"**
9. **"The GitHub README will be read by people who skip caveats. A status containing the word ADOPTED, in a repo whose theme is 'physically correct water', will be read as endorsement. Two-column reporting protects the PDF reports; it does not protect the one-line status in the docs. Is ADOPTED_WITH_BASIS honest at the exact place it will be read most carelessly?"**

## 7. Could not verify

- Kunz (2017) clean-freshwater Aeolotron MSS(U) values (thesis; cited by Engel et al. 2025 for the reference curve and the 0.0003 detection limit, but the curve's low-wind values were not extractable).
- Plant (2002)'s own treatment below u\* ≈ 0.1 m/s (no full text read; statements about it are index-level, from its abstract, its data anchors, and its use by others).
- Whether Cox & Munk's *slope* dataset specifically begins at 1 or 3 m/s (original abstract says experiments spanned 1–14 m/s; Karaev et al. 2019 say 3–14 m/s) — flagged in argument 3 above because it moves the intercept's extrapolation length.
- Jähne & Riemer (1990)'s u\* data range (carried from 2026-11-11: unverified either way).
- Donelan & Plant (2009)'s numerical threshold u\* (the invariant-threshold *character* is verified via Plant's project summary and the journal abstract; the exact threshold value was not extracted).
- Any pond-scale (fetch 100–1000 m) wavenumber spectrum at U = 1–3 m/s — none found in either round; that absence remains the fact of the matter.

## Sources

Full list with URLs and verification levels: `notes/sources.md` in this directory. Project context: `~/workspace/research_notes/low-wind-capillary-floor-20261111/report.md`; `~/workspace/caustics-filter/PIPELINE-V11.md`; `~/workspace/caustics-filter/PIPELINE-V12.md` (incl. outcome addendum); rally digests `rally-v39-{gpt,grok,gemini}.md`, `rally-v41-{gpt,grok,gemini}.md`; `delivery/offline-v15/REPORT-V15-MERGE.md` §1.
