"""lowwind_status -- the SETTLED low-wind status of the caustics pipeline.

NEW MODULE (v17).  Nothing in v11..v16 is imported, edited or re-derived here:
this file only records, in code, the wording the user settled on 2026-11-12
after rally v42 (Gemini + GPT; Grok skipped by one-time exception).

What this module is
-------------------
The pipeline floors alpha_m per cell below the Elfouhaily (1997) alpha_m fit
crossing so that the capillary slope variance is sigma2_h = 0.003, labelled
`low_wind_prior` (v13 `ALPHA_M_DOMAIN_PRIOR`, carried by the v15 cell records
under `alpha.low_wind_prior` / `provenance.alpha_m_domain`).  The rally did NOT
conclude that this prior is validated: the status is UNRESOLVED, and the prior
is retained only as an operational convention that must travel with a
sensitivity interval and must never be read as a gate pass.  The stronger
"validated" reading is deliberately absent from this module, and the phrase
that was rejected in the rally does not appear anywhere in it.

Known field-name facts established from the v15 family (v17 Stage 0 probe)
--------------------------------------------------------------------------
* the per-cell record's `provenance.provenance` carries the ARM label
  (`low_wind_absolute_gate` for the merge arm D), NOT the string
  `low_wind_prior`;
* the string `low_wind_prior` is carried at `provenance.alpha_m_domain` and
  `alpha.alpha_m_domain`;
* "the prior was active" is carried by `alpha.prior_active` (bool) and by
  `alpha.low_wind_prior` (the label dict, or None);
* cells with u* >= c_m/e (i.e. U ~ 3.0 m/s on this matrix) are low-wind cells
  in the U <= 3 sense but are NOT in the prior zone: their alpha_m is the raw
  ECKV value and their provenance is `eckv_fit_domain_raw`.

No gates, no tuning, no numbers are produced by this module.
"""

STATUS_ID = "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED"

STATUS_HEADING = ("UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; "
                  "not a gate pass)")

STATUS_TEXT = (
    "Status: UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a "
    "gate pass). For winds below the Elfouhaily et al. (1997) α_m fit domain "
    "(u* < c_m/e ≈ 0.085 m/s; U ≲ 2.7 m/s in this pipeline), the published "
    "logarithmic fit crosses zero and does not prescribe a capillary-spectrum "
    "level. The renderer retains its existing operational convention: α_m is "
    "floored per cell to give capillary slope variance σ²_h = 0.003, retaining "
    "the ECKV B_h shape. This is a modelling prior, not a measured low-wind "
    "capillary spectrum; results using affected cells must report the "
    "sensitivity interval σ²_h ∈ [0.0015, 0.006].\n"
    "\n"
    "Basis and limitation: 0.003 is borrowed as a magnitude anchor from the "
    "Cox & Munk (1954) crosswind zero-wind MSS intercept. The mechanism behind "
    "that residual is not claimed to operate on a pond, and the intercept does "
    "not directly measure the capillary-band component used here. The measured "
    "low-wind total-MSS ranges from Cox–Munk and Hu et al. provide broad "
    "context for total surface roughness; they do not independently bracket "
    "σ²_h. No published measurement currently resolves the capillary-band "
    "slope spectrum of a pond at U = 1–3 m/s.\n"
    "\n"
    "Interpretation of the zero alternative: Measurements of low-wind wave "
    "inception and nonzero surface roughness argue against treating the entire "
    "water surface as perfectly smooth. They do not establish a nonzero "
    "sustained capillary-band variance of 0.003, or rule out a very small "
    "capillary contribution under sheltered, film-covered, or below-threshold "
    "conditions.\n"
    "\n"
    "Known limits: The floor is wind- and fetch-independent where active, "
    "while real capillary levels may depend on wind, fetch, shelter, surface "
    "films and growth thresholds. The small-scale MSS measured on a sheltered "
    "lake demonstrates that substantially smoother conditions occur. "
    "Band-share gate C1 at these winds remains not passed and is printed, not "
    "waived. Every state built with the floor active carries the provenance "
    "label low_wind_prior.\n"
    "\n"
    "Supersession: This operational convention must be reconsidered when a "
    "published slope spectrum at U = 1–3 m/s resolves k > 200 rad/m and "
    "distinguishes the capillary contribution from total MSS. A faithful "
    "implementation of Kudryavtsev et al. (2013) may trigger a model-comparison "
    "review, but does not by itself count as independent observational "
    "validation."
)

STATUS_DICT = {
    "status_id": "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED",
    "sigma2_h": 0.003,
    "sensitivity_interval": [0.0015, 0.006],
    "interval_kind": "declared sensitivity interval (not an evidence bracket)",
    "provenance_label": "low_wind_prior",
    "gate_C1_low_wind": "NOT PASSED (printed, not waived)",
    "supersession": ("published slope spectrum at U = 1-3 m/s resolving "
                     "k > 200 rad/m that distinguishes the capillary "
                     "contribution from total MSS"),
}


def status_block():
    """STATUS_DICT plus the headline and the full settled wording."""
    block = dict(STATUS_DICT)
    block["heading"] = STATUS_HEADING
    block["text"] = STATUS_TEXT
    return block
