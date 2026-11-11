# Low-wind status — UNRESOLVED, operational prior retained

**Status: UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a gate pass).**

For winds below the Elfouhaily et al. (1997) α_m fit domain (u* < c_m/e ≈ 0.085 m/s; U ≲ 2.7 m/s in this pipeline), the published logarithmic fit crosses zero and does not prescribe a capillary-spectrum level. The renderer retains its existing operational convention: α_m is floored per cell to give capillary slope variance σ²_h = 0.003, retaining the ECKV B_h shape. This is a modelling prior, not a measured low-wind capillary spectrum; results using affected cells must report the sensitivity interval σ²_h ∈ [0.0015, 0.006].

**Basis and limitation.** 0.003 is borrowed as a magnitude anchor from the Cox & Munk (1954) crosswind zero-wind MSS intercept. The mechanism behind that residual is not claimed to operate on a pond, and the intercept does not directly measure the capillary-band component used here. The measured low-wind total-MSS ranges from Cox–Munk and Hu et al. provide broad context for total surface roughness; they do not independently bracket σ²_h. No published measurement currently resolves the capillary-band slope spectrum of a pond at U = 1–3 m/s.

**Interpretation of the zero alternative.** Measurements of low-wind wave inception and nonzero surface roughness argue against treating the entire water surface as perfectly smooth. They do not establish a nonzero sustained capillary-band variance of 0.003, or rule out a very small capillary contribution under sheltered, film-covered, or below-threshold conditions.

**Known limits.** The floor is wind- and fetch-independent where active, while real capillary levels may depend on wind, fetch, shelter, surface films and growth thresholds. The small-scale MSS measured on a sheltered lake demonstrates that substantially smoother conditions occur. Band-share gate C1 at these winds remains **not passed** and is printed, not waived. Every state built with the floor active carries the provenance label `low_wind_prior`.

**Supersession.** This operational convention must be reconsidered when a published slope spectrum at U = 1–3 m/s resolves k > 200 rad/m and distinguishes the capillary contribution from total MSS. A faithful implementation of Kudryavtsev et al. (2013) may trigger a model-comparison review, but does not by itself count as independent observational validation.

---

## How this status was decided (2026-11-12)

An earlier proposal would have labelled this prior "adopted with basis". A structured three-way review was run on that proposal (research note: `research_notes/low-wind-adoption-20261112` in the project records; the round is v42 in the project's rally registry). One reviewer rejected the adoption claim outright — the cited measurements constrain *total* surface roughness, not the capillary-band component this prior assigns — and the other allowed adoption only under a rename whose headline leads with UNRESOLVED. (The third reviewer's answer could not be collected; the project owner accepted a two-way record for this round as a one-time exception.) The wording above is the stricter reviewer's required replacement, adopted by the project owner. The interval [0.0015, 0.006] is a **declared sensitivity interval** around the operational centre, not an evidence bracket: the supplied measurements do not establish those as physical lower and upper bounds.

In code, this status lives in `lowwind_status.py` (`STATUS_ID = "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED"`), and an automated test (`test_lowwind_provenance.py`) checks that every low-wind output carries the `low_wind_prior` provenance and never reports an unqualified pass.
