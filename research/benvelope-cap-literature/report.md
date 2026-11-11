# Is a hard cap B(k) ≤ 0.012 a literature-supported universal bound?

## Summary

**No. The universal cap is a category error.** The number 0.012 is the *mean fitted value of
Phillips' equilibrium-range coefficient β* (Juszko, Marsden & Waddell 1995: mean β = 0.012,
observed range 0.006–0.024), a coefficient that applies **only inside the equilibrium range**
(≈1.3 f_p ≲ f ≲ 3 f_p, behind the spectral peak) and only in the frequency-spectrum form
E(f) = β·I(p)·g·u\*·f⁻⁴. No published source found proposes B(k) ≤ 0.012 as a bound on the
**total** saturation spectrum B(k) = k³S(k) including its peaks.

Worse for the gate: Elfouhaily et al. (1997) **itself** exceeds 0.012 at the gravity peak by
construction for young seas. Its measurement-fitted peak formula gives B_l(k_p) = 0.0120 at
Ω = 4.9, 0.0176 at Ω = 7.3 and 0.0241 at Ω = 10.3 — and the model's flank maxima are ≈3.7 %
higher, i.e. ≈0.0183 and ≈0.0250, which are *exactly* the mission's two "failing" gravity-peak
cells. The third failing cell, 0.01249 at the capillary peak, is likewise Elfouhaily's own
B_h(k_m) ≈ ½α_m ≈ 0.0118–0.0125 at U₁₀ = 10 m/s. The gate therefore rejects the Elfouhaily
model for producing Elfouhaily values. All three failing cells should be judged against the
region-specific criteria in Q5, against which they pass.

Confidence: high for Q1, Q2, Q4 (primary paper read in full; arithmetic reproduced below);
medium for Q3 (measured peak values are mostly published as figures/fitted formulas rather
than tables — see flags).

---

## Findings

### Q1. Elfouhaily et al. (1997): B_l peak at k = k_p vs Ω — yes, it exceeds 0.012 by construction

Source: Elfouhaily, T., Chapron, B., Katsaros, K. & Vandemark, D. (1997), "A unified
directional spectrum for long and short wind-driven waves", *J. Geophys. Res.* 102(C7),
15781–15796 (primary; full text read 2026-11-11,
https://web-wp.archive.org/web/20251106103856/http://archimer.ifremer.fr/doc/00091/20226/17877.pdf).
Equation transcription cross-checked against Mobley, *Ocean Optics Web Book*, "Wave Variance
Spectra: Examples" (https://oceanopticsbook.info/index.php/view/surfaces/level-2/wave-variance-spectra-examples),
which reproduces the ECKV equations with their original numbers.

The model (ECKV Eqs. 30–32):

- S(k) = k⁻³ [B_l + B_h]
- B_l(k) = ½ α_p (c_p/c) F_p
- F_p = L_PM · J_p · exp[ −(Ω_c/√10)(√(k/k_p) − 1) ], with
  L_PM = exp[ −1.25 (k_p/k)² ], J_p = γ^Γ, Γ = exp{ −[√(k/k_p) − 1]² / (2σ²) },
  σ = 0.08 (1 + 4 Ω_c⁻³)

At k = k_p every k-dependent factor collapses: c_p/c = 1, the cutoff exponential = 1,
Γ = exp(0) = 1 so J_p = γ, and L_PM = exp(−1.25) = 0.2865. Hence the **peak saturation level is**

> **B_l(k_p) = ½ · α_p · L_PM(k_p) · γ = 0.1433 · α_p · γ**

with the measurement-fitted (Donelan et al. 1985, Lake Ontario) parameters Elfouhaily adopts
(ECKV §3.1, Eq. 34):

- α_p = 0.006 · Ω^0.55
- γ = 1.7 for Ω ≤ 1; γ = 1.7 + 6·log₁₀(Ω) for Ω > 1

Numeric values (my arithmetic from the formulas above):

| Ω (inverse wave age) | α_p | γ | **B_l(k_p)** |
|---|---|---|---|
| 0.84 (fully developed) | 0.00545 | 1.70 | **0.00133** |
| 3.4 | 0.01176 | 4.889 | **0.00824** |
| 4.9 | 0.01438 | 5.841 | **0.01203** |
| 7.3 | 0.01791 | 6.880 | **0.01765** |
| 10.3 | 0.02164 | 7.777 | **0.02411** |

**Answer: yes.** Elfouhaily's own B_l crosses 0.012 at the peak itself at Ω ≈ 4.9 and rises
steeply thereafter, because both α_p (∝ Ω^0.55) and the overshoot factor γ (∝ log Ω) grow with
inverse wave age — that growth *is* the measured "overshoot" of young, fetch-limited seas.
The actual maximum of B_l sits on the peak flank slightly above k_p (L_PM and (k/k_p)^½ rise
faster than J_p and the cutoff fall), ≈3.7 % above B_l(k_p): 0.01765 → ≈0.0183 at Ω = 7.3 and
0.02411 → ≈0.0250 at Ω = 10.3. These are, to three significant figures, the mission's failing
cells (1) and (2). The model is not misbehaving; the gate is mis-specified.

Two caveats, neither of which rescues the gate:
- Elfouhaily/Donelan state the α_p, γ fits for 0.84 < Ω < 5 (Ocean Optics summary: young seas
  Ω = 2–5, "maximum allowed value is 5"). Ω = 7.3 and 10.3 are extrapolations of the log-γ
  branch — the extrapolation the mission model itself (and most published implementations)
  uses. A gate may legitimately flag *out-of-range Ω* as an applicability warning; it may not
  claim a literature cap of 0.012 that Elfouhaily's in-range formula already reaches at Ω ≈ 4.9.
- At Ω = 0.84 the peak is far below 0.012 (0.0013), and the mission's passing cell
  (Ω_c = 3.40, B max = 0.0047) is *below* the Elfouhaily B_l(k_p) value of 0.0082 for that Ω —
  i.e. the passing cell is if anything low against Elfouhaily, not validated by the cap.

### Q2. Where does 0.012 come from? The Phillips equilibrium-range β — and no source caps the total spectrum with it

The 0.012 in the literature is the coefficient β in Phillips' equilibrium-range spectrum:

- Phillips, O. M. (1985), "Spectral and statistical properties of the equilibrium range in
  wind-generated gravity waves", *J. Fluid Mech.* 156, 505–531: in the **equilibrium range
  only** — the band behind the peak where wind input, nonlinear transfer and breaking
  dissipation balance — E(f) = β·I(p)·g·u\*·f⁻⁴ (equivalently φ(k) = β·u\*·g^(−1/2)·k^(−5/2)),
  with β an empirical constant and I(p) ≈ 1.9–3.1 a directional-spreading factor.
- Juszko, B. A., Marsden, R. F. & Waddell (1995), as reported in Thomson, J. et al. (2013),
  "Waves and the equilibrium range at Ocean Weather Station P", *JGR Oceans* 118,
  doi:10.1002/2013JC008837 (primary; http://faculty.washington.edu/jmt3rd/Publications/Thomson_etal_2013.pdf):
  Juszko et al. "obtained a **mean value for the constant β = 0.012 from a range of
  0.006 < β < 0.024**". So 0.012 is a *mean*, not even the upper end of its own range —
  measured β reaches 0.024, twice the gate value.
- Babanin et al. (2020), "Estimating wind speed and direction using wave spectra", *JGR
  Oceans*, doi:10.1029/2019JC015717 (primary): "the observed values of β in this study
  correspond well to those referred to by Phillips (1985), with **β = [0.006 − 0.011]** and
  found by Juszko et al. (1995) with **β = 0.0122 ± 3.6×10⁻³**"; their own estimate is
  β ≈ 0.009 for U₁₀ > 6 m/s. The equilibrium range there is explicitly 1.3 f_p ≲ f ≲ 3 f_p.
- Later wind-retrieval papers simply adopt the mean: "β is set as β = 0.012 following Juszko,
  Marsden, and Waddell (1995)" (*Coastal Engineering Journal* wind-velocity paper, 2024).

So the mission's hypothesis is *almost* right: 0.012 ≈ the top of the commonly quoted
Phillips range [0.006–0.011] / the Juszko mean 0.0122. But β is **not** B(k): it is a
coefficient in a power-law valid in one spectral band, entangled with I(p) and u\*. Converting
it into a pointwise cap on k³S(k) at the peaks is a dimensional/definitional error — at the
gravity peak the spectrum is deliberately *not* in equilibrium (it carries the JONSWAP/Donelan
overshoot γ), and at the capillary peak it is set by α_m(u\*), a different parameter entirely.

**Does any published source propose B(k) ≤ 0.012 for the total spectrum including peaks?**
**None found.** Searches across the Elfouhaily, Phillips, Thomson, Hwang, Lenain–Melville and
spectrum-review literature returned 0.012 only as (i) the equilibrium β above, and (ii) an
unrelated coincidence — the **Charnock constant α = 0.012** assumed for the roughness length
z₀ = αu\*²/g in Thomson et al. (2013) and elsewhere. (A negative search result cannot prove
no source exists anywhere; it is reported as "none found", not "none exists".)

For completeness, two neighbouring constants that are also *not* caps on B(k): the
Pierson–Moskowitz Phillips constant α_PM = 0.0081 (fully developed level), and Forristall's
(1981) saturation-range coefficient α_P = 0.019 — the latter *exceeds* 0.012 in its own band
(Resio et al. parametric-spectrum review, via arXiv).

### Q3. Measured peak saturations B(k_p) for young / fetch-limited seas

Directly tabulated B(k_p) values are rare in this literature — peaks are published as
measurement-*fitted formulas* and figures. The fitted formulas are the measurements, in
parametric form:

1. **Donelan, Hamilton, Hui & Stewart (1985)**, "Directional spectra of wind-generated
   waves", *Phil. Trans. R. Soc. Lond.* A315, 509–562 — 14-staff array, fetch-limited Lake
   Ontario. Their fits α_D = 0.006·Ω^0.55 and γ(Ω) (the very formulas Elfouhaily adopts, Q1)
   give B(k_p) rising from 0.0013 (Ω = 0.84) through **0.0120 at Ω = 4.9** — i.e. measured
   developing seas reach the gate value at the top of the fitted range, on a steeply rising
   trend (∝ Ω^0.55 · log Ω). Exceeding 0.012 at the peak is the measured behaviour, not an
   anomaly.
2. **Hasselmann et al. (1973)**, JONSWAP, *Measurements of wind-wave growth and swell decay
   during the Joint North Sea Wave Project* — measured fetch-limited North Sea spectra fit by
   γ = 3.3 and α_J = 0.076·(gX/U₁₀²)^(−0.22). In saturation form B(k_p) = ½·α_J·e^(−1.25)·γ
   (same structure as Q1), this measured fit gives:
   - X = 10 km, U = 10 m/s: α_J = 0.0167 → B(k_p) = **0.0079**
   - X = 1 km, U = 10 m/s: α_J = 0.0277 → B(k_p) = **0.0131** (> 0.012)
   - X = 100 m, U = 10 m/s: α_J = 0.0460 → B(k_p) = **0.0217** (> 0.018)
   - X = 100 m, U = 6 m/s: α_J = 0.0367 → B(k_p) = **0.0174**
   (Flag: the last two extrapolate the JONSWAP α fetch law below its km-scale measurement
   fetches; the 1 km value does not. Individual JONSWAP spectra had γ varying around the 3.3
   mean, so individual measured peaks scatter above these mean-fit values — with γ ≈ 4 the
   100 m / 10 m/s case reaches ≈0.026, i.e. > 0.025.)
3. **Hwang, P. A. (2005)**, "Wavenumber spectrum of intermediate-scale ocean surface waves",
   NRL/PP/7330-05-5202 (primary; https://dtic.minsky.ai/index/ADA442354/pdf) — field
   compilation, 291 wind-sea cases (U = 3.6–14.2 m/s) plus 106 mixed-sea cases: the measured
   dimensionless spectrum follows B(k) = A₀(k)·(u\*/c)^a₀(k) over 1 < k < 316 rad/m, and for
   k above the peak Hwang adopts the equilibrium form S(k) = b·u\*·g^(−1/2)·k^(−5/2) with
   spectral constant **b = 5.2×10⁻²** — a measured saturation *coefficient* more than four
   times 0.012 (it multiplies u\*/c < 1, so it is not itself a B value; it shows how far the
   "0.012 as universal constant" idea is from measured spectral levels in this band).
4. **Jähne & Riemer (1990)** laboratory tank data (laser slope gauge, fetch ≈ 90 m in the
   tank studies cited by Elfouhaily; reproduced as Elfouhaily Fig. 5a): measured frequency
   curvature spectra B(f) show the gravity peak and the secondary gravity-capillary peak at
   30–40 Hz, with peak levels rising with wind speed into the 10⁻² range and above.
   **UNVERIFIED as exact numbers**: the values exist only as figure curves (axis to 0.1);
   they were not digitised for this report.

**Answer:** Yes — measurements (via their published fits) exceed 0.012 at the gravity peak
for fetches ≲ a few km / Ω ≳ 5; exceed 0.018 at short fetch (100 m class, U = 10 m/s mean
fit 0.0217; Donelan/Elfouhaily trend reaches 0.0176 at Ω = 7.3); and reach ≈0.025 for the
youngest cases (Elfouhaily-form flank maximum 0.0250 at Ω = 10.3; JONSWAP mean fit plus
observed γ scatter). No measurement compilation supports a ceiling anywhere near 0.012 at
the peak.

### Q4. The capillary peak: B_h(k_m) in Elfouhaily, and measured gravity-capillary peaks

Source: Elfouhaily et al. (1997), §5.2.2, Eqs. 40–44 (primary, as Q1).

- B_h(k) = ½ α_m (c_m/c) F_m, F_m = exp[ −¼ (k/k_m − 1)² ] (the printed Eq. 41's L_PM·J_p
  factors do not change the value at k_m: L_PM(k_m) = exp[−1.25(k_p/k_m)²] ≈ 1 since
  k_p ≪ k_m, and J_p(k_m) = γ^Γ with Γ ≈ 0 at k ≫ k_p, so J_p ≈ 1).
- At k = k_m ≈ 370 rad/m: c_m/c = 1 (c_m = 0.23 m/s is the minimum phase speed, at k_m) and
  the exponential = 1. Hence

> **B_h(k_m) ≈ ½ α_m**

- α_m is the two-regime logarithmic fit (Eq. 44) to the Jähne & Riemer (1990) and Hara, Bock
  & Lyzenga (1994) measured secondary-peak levels:
  α_m = 0.01·[1 + ln(u\*/c_m)] for u\* ≤ c_m; α_m = 0.01·[1 + 3·ln(u\*/c_m)] for u\* > c_m.
- At U₁₀ = 10 m/s: u\* ≈ 0.36–0.38 m/s (Elfouhaily's Garratt/Wu drag, C_d10N ≈ 0.0013–0.00144;
  equivalently u\* ≈ 0.0378·U₁₀ as used in later reproductions), so u\*/c_m ≈ 1.57–1.65 and
  **α_m ≈ 0.0235–0.0250, giving B_h(k_m) ≈ 0.0118–0.0125.**

The mission's failing cell (3), B max = 0.01249 at k = 368 rad/m ≈ k_m, is therefore
Elfouhaily's own capillary peak at U = 10 m/s, at the top of the drag-law range. A 0.012 cap
rejects it by 4 % — inside the drag-coefficient uncertainty of u\* itself.
**Flag / inconsistency in the mission statement:** the stated α_m ≈ 0.016 at U = 10 m/s would
give B_h(k_m) ≈ 0.008, which does *not* reproduce the reported 0.01249; 0.01249 is consistent
with Eq. 44 (α_m ≈ 0.025), not with α_m = 0.016. One of the two mission numbers uses a
different α_m definition (e.g. the linear form α_m = 1.4×10⁻²·(u\*/c_m) of Eq. 43 gives
0.022–0.023, also ≈ Eq. 44). This should be reconciled in the model, but either way the
literature value of the peak is ≈0.012, not safely below it.

Measured gravity-capillary peaks: Elfouhaily's Fig. 7a plots the Jähne & Riemer (1990) and
Hara et al. (1994) **measured** secondary-peak levels at 400 rad/m against u\*, and Eq. 44 is
fitted to them — i.e. the measured peak level at u\* ≈ 0.36 m/s *is* ≈ α_m ≈ 0.024 in
α_m units (≈0.012 as B_h(k_m) = ½α_m), rising to α_m ≈ 0.033 (B_h ≈ 0.017) at u\* = 0.5 m/s
and α_m ≈ 0.01·[1 + 3 ln(0.85/0.23)] ≈ 0.049 (B_h ≈ 0.025) at u\* = 0.85 m/s on the fitted
curve. Munk (2009), "An inconvenient sea truth", *Annu. Rev. Mar. Sci.* 1, 377–415, Fig. 11,
likewise summarises the measured saturation spectrum at the gravity-capillary transition as
an "enhanced wind-dependent intensity centered at the transition" — a peak, not a plateau at
any universal level. **Yes: published measurements of the gravity-capillary peak exceed
0.012** (marginally at U₁₀ = 10 m/s, clearly at higher winds). Exact digitised values from
the Jähne & Riemer / Hara figures are **UNVERIFIED** (figures not tabulated; values above are
read through Elfouhaily's published fit to them, which is the citable form).

### Q5. Synthesis: defensible acceptance criteria, by spectral region

**Yes — the universal 0.012 cap is a category error**: it takes an equilibrium-range
*coefficient* (Phillips β, a mean fit over one band, in a frequency-domain power law) and
applies it as a pointwise *ceiling* on a different quantity (B(k) = k³S(k)) in bands — the two
spectral peaks — whose defining physics is precisely that they are **not** in that equilibrium
(peak overshoot γ at k_p; α_m(u\*) forcing at k_m). Judged region by region:

| Region | Defensible criterion | Source |
|---|---|---|
| **(a) Gravity peak** (k ≈ k_p, say k_p/2–2k_p) | No absolute cap. Compare with the measurement-fitted peak level for the same Ω (or fetch): B_l(k_p) = ½·α_p·e^(−1.25)·γ, α_p = 0.006Ω^0.55, γ = 1.7 + 6 log₁₀Ω (flank max ≈ 1.04× this); equivalently JONSWAP B(k_p) = ½·α_J·e^(−1.25)·γ for the same (U, X). A model within ~a factor 1.5 of the fitted value, with Ω inside (or flagged outside) 0.84–5, passes. | Donelan et al. (1985); Hasselmann et al. (1973); Elfouhaily et al. (1997) Eqs. 31–34 |
| **(b) Equilibrium range** (≈1.3 f_p–3 f_p; k_p < k < k_n, up to ~10k_p in Elfouhaily's long-wave part) | Fit the model's band level to E(f) = β·I(p)·g·u\*·f⁻⁴ and require the *fitted coefficient* β ∈ [0.006, 0.024] (observed range; central values 0.009–0.0122). This is a range check on a fitted coefficient — the only legitimate use of "0.012" (as a central value, never a cap). In the saturation range beyond k_n (φ ∝ k⁻³), the measured level is B_s ≈ 5–9×10⁻³ (definition-dependent). | Phillips (1985); Juszko et al. (1995) via Thomson et al. (2013); Babanin et al. (2020); Banner et al. (1989) / Lenain & Melville (2017) for B_s |
| **(c) Capillary peak** (k ≈ k_m ≈ 370 rad/m) | Compare with B_h(k_m) ≈ ½·α_m(u\*), α_m from Elfouhaily Eq. 44 (log form) or Eq. 43 (linear 1.4×10⁻²·u\*/c_m): expect ≈0.012 at U₁₀ = 10 m/s, ≈0.017 at u\* = 0.5 m/s, ≈0.025 at u\* = 0.85 m/s. Judge the cell against ½α_m for *its own* u\*, with drag-law uncertainty (±~10 %) allowed. | Elfouhaily et al. (1997) Eqs. 40–44, Fig. 7a; fitted to Jähne & Riemer (1990), Hara et al. (1994) |

Legitimate *global* sanity checks (integrals, not pointwise B) do exist and are the kind of
gate Elfouhaily himself used: total mean-square slope vs Cox & Munk (1954),
mss_clean = 10⁻³ + 5.12×10⁻³·U₁₀ ± 0.004 (Elfouhaily Eq. 26), and Plant's limit
mss_clean < 0.12 (Eq. 29).

**Re-judging the three failing cells against the correct criteria:**

| Cell | B max | Correct benchmark | Verdict |
|---|---|---|---|
| U = 6, F = 100 m, Ω_c = 7.30, max at k = 14.9 (gravity flank) | 0.0183 | Elfouhaily B_l flank max ≈ 0.0183 (B_l(k_p) = 0.01765 × ≈1.04) | **Pass** — exact Elfouhaily value (with Ω-extrapolation warning, Q1) |
| U = 10, F = 100 m, Ω_c = 10.26, max at k = 10.5 (gravity peak) | 0.0250 | Elfouhaily B_l flank max ≈ 0.0250 (B_l(k_p) = 0.0241 at Ω = 10.3 × ≈1.04); JONSWAP mean fit at this (U, X) = 0.0217, above it with observed γ scatter | **Pass** — Elfouhaily value; again flag Ω > 5 extrapolation |
| U = 10, F = 1000 m, Ω_c = 4.87, max at k = 368 (capillary peak) | 0.01249 | B_h(k_m) ≈ ½α_m ≈ 0.0118–0.0125 at U₁₀ = 10 m/s | **Pass** — Elfouhaily value, within drag-law uncertainty |

---

## Could not verify

- **No source for a universal B(k) ≤ 0.012 cap**: none found in the searched literature
  (negative result; reported as "none found", not proven non-existent).
- **Exact digitised measured B(k_p) / B(k_m) values** from Jähne & Riemer (1990), Hara et al.
  (1994), Hwang's Fig. 1c–d and Donelan et al. (1985): these exist as figures, not tables;
  the report therefore uses their *published fitted formulas* (which are exact) and marks
  figure-read levels UNVERIFIED rather than inventing digits.
- **The mission's α_m ≈ 0.016 at U = 10 m/s** is inconsistent with its own reported capillary
  B max = 0.01249 under B_h(k_m) ≈ ½α_m, and with Elfouhaily Eq. 44 (α_m ≈ 0.024–0.025);
  which α_m definition the model uses could not be verified from the literature side.
- Phillips (1985) was verified through its equations and β values as quoted in Thomson et
  al. (2013), Babanin et al. (2020) and Romero & Melville (2022); the 1985 paper itself was
  not available in full text during this run.
- Elfouhaily's γ branch for Ω > 5: secondary transcriptions differ (log-branch extrapolation
  vs a 2.7·Ω^0.57 branch in some reproductions). The mission's numbers match the log-branch
  extrapolation; the in-range (Ω ≤ 5) conclusions are unaffected.

## Sources

Primary papers read in full text (2026-11-11, index/full-text fetch; Elfouhaily also via the
archived publisher PDF):

1. Elfouhaily, T., Chapron, B., Katsaros, K., Vandemark, D. (1997). A unified directional
   spectrum for long and short wind-driven waves. *J. Geophys. Res.* 102(C7), 15781–15796.
   https://web-wp.archive.org/web/20251106103856/http://archimer.ifremer.fr/doc/00091/20226/17877.pdf
2. Thomson, J., D'Asaro, E. A., Cronin, M. F., Rogers, W. E., Harcourt, R. R., Shcherbina, A.
   (2013). Waves and the equilibrium range at Ocean Weather Station P. *J. Geophys. Res.
   Oceans* 118, doi:10.1002/2013JC008837.
   http://faculty.washington.edu/jmt3rd/Publications/Thomson_etal_2013.pdf
3. Babanin, A. V. et al. (2020). Estimating wind speed and direction using wave spectra.
   *J. Geophys. Res. Oceans*, doi:10.1029/2019JC015717.
   https://minerva-access.unimelb.edu.au/rest/bitstreams/044ca7e5-9548-5c79-a070-61a301bc843f/retrieve
4. Hwang, P. A. (2005). Wavenumber spectrum of intermediate-scale ocean surface waves.
   NRL/PP/7330–05–5202. https://dtic.minsky.ai/index/ADA442354/pdf
5. Kleiss, J. M., Melville, W. K. (2010). Observations of wave breaking kinematics in
   fetch-limited seas. *J. Phys. Oceanogr.* 40, 2575–2604.
   https://airsea.ucsd.edu/wp-content/uploads/sites/10/2019/06/2010_Kleiss_Melville-Journal_of_Physical_Oceanography_vol_40.pdf
6. Romero, L., Melville, W. K. (2022). Mass transfer at the ocean–atmosphere interface…
   *Annu. Rev. Fluid Mech.* 54, 191–224. https://par.nsf.gov/servlets/purl/10388564

Cited through the above / through scholarly secondary transcription (flagged in text where
relied on): Phillips (1985) *J. Fluid Mech.* 156, 505–531; Juszko, Marsden & Waddell (1995);
Hasselmann et al. (1973) JONSWAP; Donelan et al. (1985) *Phil. Trans. R. Soc.* A315, 509–562;
Jähne & Riemer (1990); Hara, Bock & Lyzenga (1994); Cox & Munk (1954); Banner et al. (1989);
Lenain & Melville (2017) *J. Phys. Oceanogr.*, doi:10.1175/JPO-D-17-0017.1; Munk (2009)
*Annu. Rev. Mar. Sci.* 1, 377–415 (https://airsea.ucsd.edu/wp-content/uploads/sites/10/2020/04/Munk-2009.pdf);
Mobley, Ocean Optics Web Book, Wave Variance Spectra: Examples
(https://oceanopticsbook.info/index.php/view/surfaces/level-2/wave-variance-spectra-examples);
Williams, Mueller & Deike, Observations and empirical functions for the ocean surface wave
spectrum (https://arxiv.org/html/2606.28165).

Working notes: notes/calculations.md (all arithmetic), notes/sources.md (source log).
