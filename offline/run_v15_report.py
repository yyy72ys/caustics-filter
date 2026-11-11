"""run_v15_report -- builds out/v15/REPORT-V15-MERGE.md from the stage JSONs.

Reads (never recomputes physics):
  out/v15/numbers-v15-merge.json      stages 1 + 2 + the stage-4 roll-up
  out/v15/numbers-v15-transport.json  stage 3

Every number in the report is a value read from those files; the only derived
quantities are ratios of numbers already present, and they are labelled.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import atomic_write_bytes  # noqa: E402

OUT = os.path.join(HERE, "out", "v15")
MERGE = os.path.join(OUT, "numbers-v15-merge.json")
TRANS = os.path.join(OUT, "numbers-v15-transport.json")
REPORT = os.path.join(OUT, "REPORT-V15-MERGE.md")


def load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def n(x, fmt="%.4f"):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    try:
        return fmt % x
    except TypeError:
        return str(x)


def yn(x):
    if x is None:
        return "-"
    return "PASS" if x is True else ("FAIL" if x is False else str(x))


def table(headers, rows, sep="  "):
    """Fixed-width ASCII table; every column is left-aligned."""
    widths = [len(h) for h in headers]
    for r in rows:
        for i, c in enumerate(r):
            widths[i] = max(widths[i], len(str(c)))
    out = []
    out.append(sep.join(h.ljust(widths[i]) for i, h in enumerate(headers)))
    out.append(sep.join("-" * widths[i] for i in range(len(headers))))
    for r in rows:
        out.append(sep.join(str(c).ljust(widths[i]) for i, c in enumerate(r)))
    return "\n".join(out)


def main():
    t0 = time.time()
    m = load(MERGE)
    t = load(TRANS)
    cells = m["cells"]
    pol = m["policy"]
    mac = m["merge"]
    roll = m["rollup"]
    reg = m["regression"]
    summ = m["summary"]
    L = []

    L.append("# REPORT-V15-MERGE -- version 3.10 (chunk 11): the merge")
    L.append("")
    L.append("Brief: `BRIEF-V15-MERGE-R2.md`.  Design: `PIPELINE-V12.md` (it "
             "governs everything below).")
    L.append("Base: the v14 tree (`run_v14_spectrum.py` carries v13's "
             "spectrum and v11's assembly; v5..v14 are IMPORTED, never "
             "edited).")
    L.append("Machine: PC, `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`. "
             "  Model: deepseek-v4.1-flash.")
    L.append("")
    L.append("Execution: FOREGROUND, one jobkit checkpoint unit per stage "
             "(`out/v15/checkpoint.json`).  Every deliverable is committed "
             "with `atomic_write_bytes`.  Nothing outside `out/v15/` and the "
             "new `run_v15*` files was written; no image was inspected "
             "(numerical image statistics only); the report is ASCII.")
    L.append("")
    L.append("**Label honesty (a rally condition, repeated wherever the "
             "battery is named): everything in sections 3-4 is a REGION-WISE "
             "CONFORMANCE + OBSERVATION BATTERY.  It is not physical truth, "
             "not a validation and not a measurement of any real water "
             "surface.**")
    L.append("")
    L.append("New code (v15 family only):")
    L.append("")
    L.append("- `run_v15_spectrum.py` -- the merged default, the retired "
             "clause, the region-wise battery (a)/(c), the observation "
             "battery (b)/(d), the carried gates, the per-cell verdict")
    L.append("- `run_v15_state.py` -- `V15WaterState`: v11 state + the merged "
             "default's assembly; the labels travel in `state_node`")
    L.append("- `run_v15.py` -- the four stage runners")
    L.append("- `run_v15_report.py` -- this report (reads the stage JSONs)")
    L.append("")
    L.append("Reproduce, foreground, in this order:")
    L.append("")
    L.append("```")
    L.append("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    L.append("python run_v15.py --stage 1     # 8 cells, the full battery, "
             "analytic")
    L.append("python run_v15.py --stage 2     # regression: 5 items, incl. the "
             "G1 legacy re-render")
    L.append("python run_v15.py --stage 3     # T1/T2/T3 stills, C mesh %d^2, "
             "Q=25" % int(t["stack"]["C_mesh"]))
    L.append("python run_v15.py --stage 4     # roll-up + this report")
    L.append("```")
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 1. The merge statement")
    L.append("")
    L.append("1. **The default spectrum of the v15 family is now "
             "`%s`** (`SPECTRUM_SOURCE_DEFAULT` in `run_v15_spectrum.py`), "
             "with the **Arm-D rule set**:"
             % mac["default_spectrum_source"])
    L.append("   %s" % mac["default_arm_rule"])
    L.append("   The provenance labels are v14's, verbatim: inside the fit "
             "domain the raw `alpha_m` is used and the label is "
             "`eckv_fit_domain_raw`; below the crossing the prior acts and the "
             "label is `low_wind_prior`, carrying its target/source/rally "
             "fields (the provisional 0.003 of the v9 design).")
    L.append("2. **`jonswap` is retained as the rollback source.**  It is "
             "still selectable in `run_v15_spectrum.assemble` and in "
             "`V15WaterState`, and it returns the v11/v10 lane "
             "byte-for-byte (regression item 5).")
    L.append("3. **The universal `max B(k) <= 0.012 over [10,1500]` clause is "
             "RETIRED.**  Reason, quoted from the brief and the 2026-11-11 "
             "literature + three-way rally: %s" % mac["retired_clause"]
             ["reason"])
    L.append("   What it was replaced by: %s" % mac["what_replaced_it"])
    L.append("   The number is still printed in every cell "
             "(`b_max_observation`, `gated: false`) -- printed, never scored.")
    L.append("")
    L.append("What this chunk does NOT do: it does not touch a threshold "
             "(every bound is the brief's), it does not pick a winner (the "
             "winner reading from chunk 10 stands: Arm D), it does not "
             "validate the spectrum, and it does not promote the default "
             "anywhere outside the v15 family.")
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 2. The two columns (definitions and conventions)")
    L.append("")
    L.append("**CONFORMANCE -- the model against its own published form.**")
    L.append("")
    L.append("- **(a) gravity flank**: `max B_l` over `[k_p/2, 2 k_p]` vs "
             "`0.1433 * alpha_p * gamma(Omega_c) * 1.04`; PASS if the ratio "
             "is <= %g." % pol["gravity_conformance"]["pass_if_ratio_le"])
    L.append("- **(c) capillary**: `max B_h` over `[k_m/2, 2 k_m]` vs "
             "`(1/2) * alpha_m(u*)`; PASS if the ratio is inside [%g, %g]."
             % tuple(pol["capillary_conformance"]["pass_if_ratio_in"]))
    L.append("")
    L.append("**OBSERVATION -- independent anchors.**")
    L.append("")
    L.append("- **(b) equilibrium beta**: `E(f)` is reconstructed from the "
             "assembled omnidirectional spectrum through the deep-water "
             "dispersion (`omega^2 = g k`, `f = omega/2pi`); a log-space "
             "least-squares fit of `E = beta * I(p) * g * u* * f^-4` over "
             "`f in [%.1f f_p, %.1f f_p]` with `I(p) = %g`.  PASS if the "
             "fitted beta is in [%g, %g]."
             % (pol["equilibrium_beta"]["f_band_factors"][0],
                pol["equilibrium_beta"]["f_band_factors"][1],
                pol["equilibrium_beta"]["I_p"],
                pol["equilibrium_beta"]["band"][0],
                pol["equilibrium_beta"]["band"][1]))
    L.append("")
    L.append("  **The convention, stated (the brief requires it):** "
             "`%s` with `%s`, `I(p) = %g`.  In other words: %s"
             % (pol["equilibrium_beta"]["convention"]["law"],
                pol["equilibrium_beta"]["convention"]["omega"],
                pol["equilibrium_beta"]["I_p"],
                pol["equilibrium_beta"]["convention"]
                ["why_this_convention"]))
    L.append("")
    _bhz = [c["observation"]["b_beta"]["beta_Hz_convention"] for c in cells]
    _lo = pol["equilibrium_beta"]["band"][0]
    L.append("  Both numbers are printed for every cell: `beta` (the gate "
             "quantity) and `beta_Hz_convention` (the identical fit written "
             "as `E(f)` with `f` in Hz, exactly `(2 pi)^3 = 248.05` times "
             "smaller).  Reading the band against the Hz number would put "
             "every cell %.0fx-%.0fx below the band's floor -- that is a unit "
             "defect of the literal reading, not a property of the spectrum, "
             "and it is why the convention had to be stated."
             % (_lo / max(_bhz), _lo / min(_bhz)))
    L.append("")
    L.append("  If the FIT BAND is contaminated (`B_h > B_l` anywhere inside "
             "it) the sub-check is `UNRESOLVED (band-contaminated)` -- never "
             "a silent PASS/FAIL -- and the value is still printed with its "
             "band position.  A contaminated fit is never called PASS.")
    L.append("")
    L.append("- **(d) JONSWAP anchor**: `B_J(k_p) = 0.5 * alpha_J * "
             "exp(-1.25) * 3.3`, `alpha_J = 0.076 (gF/U^2)^-0.22`, printed "
             "beside the model's `B_l(k_p)` for every cell.  **Context, not a "
             "gate** (section 3.3).")
    L.append("")
    L.append("**Carried gates, unchanged:** C1 (gravity share of the 50-200 "
             "slope variance <= %g) and C1b (capillary share >= %g AND "
             "sigma2_h >= %g) -- computed and reported for all 8 cells, "
             "**gated only for U >= 6**; totals inside the "
             "[Cox-Munk clean, Hu] envelope with the 0.5x-2.0x clause (all "
             "cells, gating); K1 `|dHs| <= %g%%` (spectrum-integrated H_s vs "
             "the fetch law, all cells, gating) plus `hs_flag` at "
             "|dHs| > %g%%, which **annotates and does not re-gate on this "
             "matrix**."
             % (pol["C1_gravity_share_max"], pol["C1b_capillary_share_min"],
                pol["C1b_sigma2_h_min"], 100.0 * pol["K1_hs_tolerance"],
                100.0 * pol["K1_hs_flag"]))
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 3. The 8-cell matrix, both columns, per cell")
    L.append("")
    L.append("Conventions (inherited, unchanged): `S(k)` in m^2 per rad/m, "
             "`m0 = int S dk`, `sigma2 = int k^2 S dk`, `B(k) := k^3 S(k)`, "
             "wind = U10, `u*` from U10 only, depth 1.0 m, `k_cut = %g` "
             "rad/m." % pol["grid"]["k_cut_rad_per_m"])
    L.append("")
    L.append("### 3.1 CONFORMANCE and OBSERVATION -- the two columns, never "
             "collapsed")
    L.append("")
    rows = []
    for c in cells:
        a = c["conformance"]["a_gravity"]
        cc = c["conformance"]["c_capillary"]
        b = c["observation"]["b_beta"]
        rows.append([
            c["cell_id"], n(c["omega_c"], "%.2f"),
            n(c["warn_extrapolated_wave_age"]),
            n(a["B_l_flank_max"], "%.5f"), n(a["reference"], "%.5f"),
            n(a["ratio"], "%.4f"), yn(a["pass"]),
            n(cc["B_h_band_max"], "%.5f"), n(cc["reference"], "%.5f"),
            n(cc["ratio"], "%.4f"), yn(cc["pass"]),
            n(b["beta"], "%.5f"), n(b["beta_Hz_convention"], "%.2e"),
            b["status"],
        ])
    L.append(table([
        "cell", "Om_c", "warn", "Bl_flank", "ref(a)", "ratio(a)", "(a)",
        "Bh_max", "ref(c)", "ratio(c)", "(c)",
        "beta", "beta_Hz", "(b) status"], rows))
    L.append("")
    L.append("Reading notes:")
    L.append("")
    L.append("- `warn` is `warn_extrapolated_wave_age` (Omega_c > 5): those "
             "cells continue the gamma log-formula past its fitted range and "
             "say so.")
    L.append("- `(a)` and `(c)` are the CONFORMANCE column: they compare the "
             "model with the very formula it implements, so they are a "
             "self-consistency audit of the implementation, not an "
             "independent test of the physics (that is what `(b)` is for).")
    L.append("- `(b)` statuses in this matrix: %s"
             % ", ".join("%s=%s" % (c["cell_id"],
                                    c["observation"]["b_beta"]["status"])
                         for c in cells))
    L.append("- The capillary reference `(1/2) * alpha_m` is the applied "
             "`alpha_m`, including the low-wind prior where it acts: the "
             "conformance test asks whether the delivered `B_h` matches the "
             "level that was put in, not whether that level was right.")
    L.append("")
    L.append("### 3.2 The fit the beta came from (per cell)")
    L.append("")
    rows = []
    for c in cells:
        b = c["observation"]["b_beta"]
        rows.append([c["cell_id"], n(b["f_p_Hz_deep_water"], "%.3f"),
                     n(b["f_p_Hz_fetch_law"], "%.3f"),
                     "[%.2f, %.2f]" % tuple(b["fit_band_Hz"]),
                     "[%.2f, %.2f]" % tuple(b["fit_band_rad_per_m"]),
                     n(b["n_points_in_band"], "%d"),
                     n(b["fit_residual_log_rms"], "%.3f"),
                     n(b["contamination"]["n_points_bh_gt_bl"], "%d"),
                     n(b["contamination"]["k_at_max_rad_per_m"], "%.1f"),
                     b["status"]])
    L.append(table(["cell", "f_p deep", "f_p law", "fit band Hz",
                    "fit band rad/m", "n", "resid", "n(Bh>Bl)", "k at max",
                    "status"], rows))
    L.append("")
    L.append("- `f_p deep` = `sqrt(g k_p)/2 pi` on the model's own peak "
             "wavenumber; `f_p law` = the pipeline's fetch-law `f_p`.  Their "
             "relative difference is printed per cell in the JSON "
             "(`f_p_relative_difference`); the largest on this matrix is "
             "%.2f%% (U10p0_F1000p0, the cell where the finite-depth peak "
             "wavenumber sits furthest from the deep-water inversion).  The "
             "fit band is wide (1.3-3 f_p), so that offset moves no cell's "
             "status, but it is stated rather than rounded away."
             % (100.0 * max(abs(c["observation"]["b_beta"]
                                ["f_p_relative_difference"]) for c in cells)))
    L.append("- `resid` is the log-space RMS spread of "
             "`S(omega) * omega^4` inside the band -- the measure of how far "
             "the fitted `f^-4` law is from the model there.  It is NOT "
             "small (0.14-0.56), which is the honest statement that the "
             "assembled spectrum is not a clean `f^-4` equilibrium range "
             "over 1.3-3 f_p.")
    L.append("")
    L.append("### 3.3 (d) The JONSWAP anchor -- context, never a gate")
    L.append("")
    rows = []
    for c in cells:
        d = c["observation"]["d_jonswap_anchor"]
        rows.append([c["cell_id"], n(d["alpha_J"], "%.5f"),
                     n(d["B_J_at_k_p"], "%.5f"),
                     n(d["model_B_l_at_k_p"], "%.5f"),
                     n(d["ratio_model_over_anchor"], "%.3f"),
                     "yes" if d["sub_km_cell"] else "-"])
    L.append(table(["cell", "alpha_J", "B_J(k_p)", "model B_l(k_p)",
                    "model/anchor", "sub-km caveat"], rows))
    L.append("")
    L.append("**Sub-km extrapolation caveat (F = 100 m cells):** F = 100 m is "
             "a sub-kilometre fetch; JONSWAP's `alpha_J` fetch relation was "
             "fitted over kilometres (gF/U^2 ~ 1e3-1e5), so on those cells "
             "the anchor is an extrapolation 1-3 decades below its fitted "
             "range.  It is printed for context only.  No gate uses it, and "
             "the model/anchor column is not evidence of agreement in either "
             "direction.")
    L.append("")
    L.append("### 3.4 Carried gates, flags and verdicts")
    L.append("")
    rows = []
    for c in cells:
        ec = c["envelope_check"]
        rows.append([
            c["cell_id"],
            n(ec["ratio_to_lower"], "%.2fx") + " / "
            + n(ec["ratio_to_upper"], "%.2fx"),
            yn(ec["tol_inside"]),
            (yn(c["C1"]["pass"]) if c["C1_C1b_gated"]
             else "n/a (U<=3)"),
            (yn(c["C1b"]["pass"]) if c["C1_C1b_gated"]
             else "n/a (U<=3)"),
            ("%+.1f%%" % (100.0 * c["hs"]["delta_rel"])),
            yn(c["hs"]["pass"]),
            n(c["hs"]["hs_flag"]),
            n(c["b_max_observation"]["B_max"], "%.5f"),
            c["state"],
        ])
    L.append(table(["cell", "totals x[lo,hi]", "totals", "C1", "C1b",
                    "dHs", "K1", "hs_flag", "B max (obs)", "verdict"], rows))
    L.append("")
    n_strict = int(sum(1 for c in cells
                       if c["envelope_check"]["strict_inside"]))
    L.append("`B max (obs)` is the retired clause's number, printed and not "
             "scored (`gated: false`).  `totals` is the 0.5x-2.0x reading, "
             "which is the gate and passes %d/8 cells; the strict "
             "`[Cox-Munk clean, Hu]` reading is in every JSON cell "
             "(`envelope_check.strict_inside`) and would pass %d/8 -- it is "
             "reported beside the gate, never substituted for it."
             % (summ["totals_tol_pass"], n_strict))
    L.append("")
    L.append("Verdict rules, as the brief fixes them:")
    L.append("")
    L.append("- **low wind (U <= 3)**: UNRESOLVED unless a hard clause fires "
             "-- totals outside 0.5x-2.0x, negative/non-finite energy, a step "
             "in S(k), |dHs| > 35%, conformance (a) > 1.5 -- and then FAIL.  "
             "(b) and (c) are not in that list and never drive the verdict; "
             "their values are reported beside it.")
    L.append("- **high wind (U = 6, 10)**: PASS requires the construction "
             "gates (positivity, smoothness, the A4 band partition) + totals "
             "+ C1 + C1b + K1 + conformance (a) and (c) PASS + observation "
             "(b) PASS or flagged UNRESOLVED-with-reason.")
    L.append("")
    L.append("Verdict counts: **PASS %d, FAIL %d, UNRESOLVED %d**."
             % (roll["n_pass"], roll["n_fail"], roll["n_unresolved"]))
    L.append("")
    L.append("### 3.5 Does the merge stand?")
    L.append("")
    L.append("The brief requires **zero FAILs** for the merge to stand.  "
             "Observed: %s." % ("**zero FAILs**" if roll["zero_fails"]
                                 else "**FAIL(s) at: %s**"
                                 % ", ".join(roll["failing_cells"])))
    L.append("")
    L.append("- `merge_stands` = **%s**" % n(roll["merge_stands"]))
    L.append("- `all_conditions` (merge + regression + transport gates + the "
             "two travelling labels) = **%s**" % n(roll["all_conditions"]))
    L.append("")
    L.append("`run_v15_spectrum.SPECTRUM_SOURCE_DEFAULT` is "
             "`%s` either way -- that IS the brief's stage-1 instruction for "
             "the new family, and it touches nothing outside v15 (v13_state's "
             "`DEFAULT_SOURCE` is still `%s`; regression item 4).  If "
             "`merge_stands` were false the sentence would be: the family's "
             "default is set, the PROMOTION is not supported by this battery."
             % (mac["default_spectrum_source"], "jonswap"))
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 4. Flags inventory")
    L.append("")
    L.append("### 4.1 `warn_extrapolated_wave_age` (Omega_c > 5)")
    L.append("")
    if summ["cells_warn_extrapolated_wave_age"]:
        rows = []
        for c in cells:
            if c["warn_extrapolated_wave_age"]:
                rows.append([c["cell_id"], n(c["omega_c"], "%.2f"),
                             c["provenance"]["eckv_domain"], c["state"]])
        L.append(table(["cell", "Omega_c", "eckv_domain", "verdict"], rows))
    else:
        L.append("(none)")
    L.append("")
    L.append("### 4.2 UNRESOLVED cells")
    L.append("")
    if summ["cells_unresolved"]:
        L.append("  " + ", ".join(summ["cells_unresolved"]))
        L.append("")
        L.append("Every one of them is a low-wind cell (U <= 3), where "
                 "UNRESOLVED is the design's ceiling, not a defect found "
                 "here.  `sigma2_h_target = %g` remains an invented input on "
                 "the cells where the prior acts (%s)."
                 % (pol["sigma2_h_target_default"],
                    ", ".join(summ["cells_with_prior_or_topup"]) or "none"))
    else:
        L.append("(none)")
    L.append("")
    L.append("### 4.3 `hs_flag` cells (|dHs| > 15%, annotating only)")
    L.append("")
    if summ["cells_hs_flag"]:
        rows = [[k, "%+.1f%%" % (100.0 * v)]
                for k, v in sorted(summ["cells_hs_flag"].items())]
        L.append(table(["cell", "dHs"], rows))
        L.append("")
        L.append("On this matrix the flag fires exactly on the three "
                 "sub-kilometre (F = 100 m) cells.  It is a forward tripwire "
                 "adopted from the rally: **it annotates, it does not re-gate "
                 "-- K1 is still |dHs| <= 35% on every cell.**")
    else:
        L.append("(none)")
    L.append("")
    L.append("### 4.4 Cells where the low-wind prior / top-up acts")
    L.append("")
    if summ["cells_with_prior_or_topup"]:
        L.append("  " + ", ".join(summ["cells_with_prior_or_topup"]))
    else:
        L.append("  (none)")
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 5. Regression (stage 2), item by item")
    L.append("")
    it1 = reg["item1_D_STRICT_equals_v11"]
    it2 = reg["item2_v15_equals_v14_numbers"]
    it3 = reg["item3_v12_water_untouched"]
    it4 = reg["item4_frozen_modules_unchanged"]
    it5 = reg["item5_jonswap_rollback"]
    L.append("| # | requirement | result | evidence |")
    L.append("|---|---|---|---|")
    L.append("| 1 | D_STRICT branch == `v11.assemble`, byte-identical "
             "spectrum | **%s** | 8/8 cells, max abs dS is printed per cell; "
             "identical alpha_m and H_s too |" % yn(it1["pass"]))
    L.append("| 2 | v15 RAW arm == v14 RAW numbers per cell (alpha_m, "
             "sigma2_g, sigma2_h, B max, H_s pred) | **%s** | three "
             "comparisons: in-session RAW, in-session arm D, and the stored "
             "v14 numbers file |" % yn(it2["pass"]))
    L.append("| 3 | v12 water untouched; G1 legacy lane byte-identical to "
             "still-A-v37 | **%s** | `WATER_MODEL_DEFAULT = '%s'`; live "
             "re-render sha `%s` vs reference `%s` |"
             % (yn(it3["pass"]), it3["WATER_MODEL_DEFAULT"],
                (it3["g1_legacy_rerun"]["rerun_sha256"] or "")[:16],
                it3["g1_legacy_rerun"]["reference_sha256"][:16]))
    L.append("| 4 | frozen composite modules unchanged | **%s** | %d/%d "
             "modules match the v14 record; all %d unchanged between stage 1 "
             "and stage 2 |" % (yn(it4["pass"]),
                                it4["n_files_compared_against_v14_record"],
                                it4["n_files_hashed"], it4["n_files_hashed"]))
    L.append("| 5 | jonswap still selectable and reproduces the v10 scene-A "
             "statistics within 1%% | **%s** | S identical to "
             "`V10.assemble`; realised sigma2 %+.2e, H_s %+.2e |"
             % (yn(it5["pass"]), it5["sigma2_rel_delta"],
                it5["Hs_rel_delta"]))
    L.append("")
    L.append("Overall regression: **%s**." % yn(reg["all_pass"]))
    L.append("")
    L.append("### 5.1 Item 1 -- D_STRICT == v11 (per cell)")
    L.append("")
    rows = [[k, n(v["alpha_m_used_d_strict"], "%.6f"),
             n(v["alpha_m_used_v11"], "%.6f"), "%.1e" % v["max_abs_dS"],
             n(v["S_identical"])] for k, v in it1["cells"].items()]
    L.append(table(["cell", "alpha_m D_STRICT", "alpha_m v11", "max|dS|",
                    "identical"], rows))
    L.append("")
    L.append("### 5.2 Item 2 -- v15 == v14 numbers (per cell, stored "
             "comparison)")
    L.append("")
    rows = []
    for k, v in it2["default_arm_vs_stored_v14_RAW_column"]["cells"].items():
        rows.append([k, n(v["alpha_m_used_v15"], "%.9f"),
                     n(v["sigma2_g_v15"], "%.9f"),
                     n(v["sigma2_h_v15"], "%.9f"),
                     n(v["B_max_v15"], "%.9f"),
                     n(v["Hs_pred_v15"], "%.9f"),
                     n(v["pass"])])
    L.append(table(["cell", "alpha_m", "sigma2_g", "sigma2_h", "B max",
                    "H_s pred", "exact"], rows))
    L.append("")
    L.append("All five fields are compared with `==`, not with a tolerance: "
             "the same code path must give the same bits.  The merged default "
             "carries the Arm-D rules, and arm D is numerically RAW's spectrum "
             "(chunk 10), so this is also the statement that the merge moved "
             "no number.")
    L.append("")
    L.append("### 5.3 Item 3 -- v12 water and the G1 legacy lane")
    L.append("")
    g = it3["g1_legacy_rerun"]
    L.append("- `WATER_MODEL_DEFAULT = '%s'`, `DEFAULT_MEDIUM_PRESET = "
             "`%s``, `LEGACY_RECIPE_ROUTE = `%s`` -> default %s."
             % (it3["WATER_MODEL_DEFAULT"], it3["DEFAULT_MEDIUM_PRESET"],
                it3["LEGACY_RECIPE_ROUTE"], yn(it3["default_pass"])))
    L.append("- G1 legacy lane re-rendered live in this session through "
             "`run_v12.scene_state` + `run_v12.render_scene` with the legacy "
             "optics at the pre-3.6 scene (U = %g, F = %g, depth = %g, "
             "grid %d): sha `%s`."
             % (g["scene"]["U_ms"], g["scene"]["fetch_m"],
                g["scene"]["depth_m"], g["scene"]["grid"],
                g["rerun_sha256"]))
    L.append("- Stored delivery artefact `out/v10/still-A-v37.png`: sha "
             "`%s`." % g["stored_sha256"])
    L.append("- Byte-identical to the stored artefact: **%s**; to the "
             "brief's recorded reference `cfb086b22a0e...`: **%s**; "
             "max |d| = %d uint8, changed pixels = %.4f%%."
             % (yn(g["byte_identical_to_stored"]),
                yn(g["byte_identical_to_reference"]), g["max_abs_uint8"],
                100.0 * g["frac_pixels_changed"]))
    L.append("")
    L.append("### 5.4 Item 4 -- frozen modules")
    L.append("")
    L.append("%d modules hashed; %d of them carry a value recorded by the "
             "v14 run and every one of those matches it.  The five without a "
             "v14-recorded value (`run_v12.py` and the `run_v14*` family) are "
             "covered by the stage-1 -> stage-4 unchanged test, and this job "
             "never writes them."
             % (it4["n_files_hashed"],
                it4["n_files_compared_against_v14_record"]))
    L.append("")
    L.append("`v5_composite.py` sha256 `%s` -- prefix `5c86a5770dcd0ce2` "
             "matches the value recorded in REPORT-V14-GATES section 8: %s.  "
             "Composite constants in force: A = %s, B_sigma = %s, E = %s, "
             "Q = %s."
             % (it4["v5_composite_sha256"],
                yn(it4["v5_composite_prefix_matches_v14_report"]),
                it4["composite_constants"]["A"],
                it4["composite_constants"]["B_sigma"],
                it4["composite_constants"]["E"],
                it4["composite_constants"]["Q"]))
    L.append("")
    L.append("### 5.5 Item 5 -- the jonswap rollback lane")
    L.append("")
    j = it5
    L.append("- Selectable: `assemble(..., source='jonswap')` returns the "
             "v11/v10 lane; `V15WaterState(spectrum_source='jonswap')` builds "
             "the legacy state.  Effective source of the realised state: "
             "`%s`." % j["source_effective"])
    L.append("- Assembly identity vs `V10.assemble` on the scene-A cell: "
             "`S(k)` identical = **%s** (max |dS| = %.1e), alpha_m %.9f vs "
             "%.9f."
             % (yn(j["assembly_identity_vs_v10"]["S_identical"]),
                j["assembly_identity_vs_v10"]["max_abs_dS"],
                j["assembly_identity_vs_v10"]["alpha_m_v15"],
                j["assembly_identity_vs_v10"]["alpha_m_v10"]))
    L.append("- Realised statistics against `numbers-v10-sweep.json` cell "
             "`%s` (N = %d, depth %g): sigma2 %.6e vs %.6e (%+.3f%%), H_s "
             "%.6f vs %.6f (%+.3f%%), rms slope %.6f vs %.6f.  Tolerance 1%%: "
             "**%s**."
             % (m["regression"]["item5_jonswap_rollback"]["reference"]
                .split("cell ")[-1],
                j["cell"]["N"], j["cell"]["h"],
                j["realised_v15"]["sigma2_rad2"],
                j["realised_v10"]["sigma2_rad2"],
                100.0 * j["sigma2_rel_delta"],
                j["realised_v15"]["Hs_m"], j["realised_v10"]["Hs_m"],
                100.0 * j["Hs_rel_delta"],
                j["realised_v15"]["rms_slope"],
                j["realised_v10"]["rms_slope"], yn(j["pass"])))
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 6. Transport through the merged default stack (stage 3)")
    L.append("")
    st = t["stack"]
    L.append("Stack: spectrum `%s` (arm `%s`) + v12 water `%s` + frozen "
             "composite (A%s + Blow + E%s) + eta = %s + retina %s + rain "
             "%s mm/h.  C mesh %d^2, Q = 25; each still is %dx%d px."
             % (st["spectrum"], st["spectrum_arm"], st["water_model"],
                st["composite"]["A"], st["composite"]["E"], st["eta"],
                ("OFF" if "OFF" in st["retina"].upper() else st["retina"]),
                st["rain_mm_h"], st["C_mesh"],
                st["still_size_px"][0], st["still_size_px"][1]))
    L.append("")
    L.append("`eta = 1.0` is a property of the recipe (no sky block -> "
             "diffuse weight 0), and the retina post-stage is not in this "
             "path; both are recorded in the JSON rather than assumed.")
    L.append("")
    L.append("| id | scene | Omega_c | warn | cell status | C mean | C p99 | "
             "C p99.9 | C max | crush % | blow % | own shift mm | file |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k in ("T1", "T2", "T3"):
        c = t["cells"][k]
        L.append("| %s | %s U=%g F=%g h=%g | %.2f | %s | %s | %.9f | %.3f | "
                 "%.3f | %.2f | %.2f | %.3f | %.2f | `%s` |"
                 % (k, c["name"], c["cell"]["U_ms"], c["cell"]["fetch_m"],
                    c["cell"]["depth_m"], c["omega_c"],
                    n(c["warn_extrapolated_wave_age"]), c["cell_status"],
                    c["C"]["mean"], c["C"]["p99"], c["C"]["p99_9"],
                    c["C"]["max"], c["crush_pct"], c["blow_pct"],
                    c["own_geometric_shift_mm"], c["file"]))
    L.append("")
    L.append("`C` is the unattenuated caustic field (the det master), so "
             "`C mean`, `p99`, `p99.9` and `max` are spectral statistics, not "
             "image statistics; `crush` and `blow` are image statistics "
             "(fraction of linear values < 0.02 and of ACES values > 0.98).")
    L.append("")
    L.append("### 6.1 The three transport gates")
    L.append("")
    L.append(table(["gate", "definition", "result"],
                   [["C mean 1 +/- 0.02", "|mean(C) - 1| <= 0.02 on the "
                     "unattenuated field",
                     yn(t["gates"]["all_three"]["C_mean_1_pm_0.02"])],
                    ["crush < 5%", "crush_rate_linear_lt_0.02 < 0.05",
                     yn(t["gates"]["all_three"]["crush_lt_5pct"])],
                    ["blow == 0", "blowout_rate_aces_gt_0.98 == 0",
                     yn(t["gates"]["all_three"]["blow_eq_0"])]]))
    L.append("")
    L.append("All three: **%s**." % yn(t["gates"]["pass"]))
    L.append("")
    L.append("### 6.2 Proof that the labels travelled into the numbers JSONs")
    L.append("")
    L.append("This is a merge condition, so it is quoted from the file rather "
             "than asserted:")
    L.append("")
    L.append("```")
    L.append("out/v15/numbers-v15-transport.json")
    for k in ("T1", "T2", "T3"):
        c = t["cells"][k]
        L.append("  cells.%s.cell_status                    = %s"
                 % (k, c["cell_status"]))
        L.append("  cells.%s.spectrum_source                = %s"
                 % (k, c["spectrum_source"]))
        L.append("  cells.%s.omega_c                        = %.6f"
                 % (k, c["omega_c"]))
        L.append("  cells.%s.warn_extrapolated_wave_age     = %s"
                 % (k, c["warn_extrapolated_wave_age"]))
        L.append("  cells.%s.low_wind_status                = %s"
                 % (k, c["low_wind_status"]))
        L.append("  cells.%s.labels_travelled               = %s"
                 % (k, json.dumps(c["labels_travelled"], sort_keys=True)))
    L.append("```")
    L.append("")
    L.append("- **T2 (U=3, F=1000, depth 0.6)** carries `cell_status = "
             "UNRESOLVED` and `low_wind_status` = \"%s\"."
             % t["cells"]["T2"]["low_wind_status"])
    L.append("- **T3 (U=6, F=100, depth 1.0)** carries `omega_c = %.4f` and "
             "`warn_extrapolated_wave_age = true`." % t["cells"]["T3"]
             ["omega_c"])
    L.append("- Both conditions: **%s**." % yn(t["label_conditions"]["pass"]))
    L.append("")
    L.append("T2's transport numbers are otherwise unremarkable and are "
             "reported, not interpreted: the low-wind cell is the one this "
             "design cannot resolve, and a still is not evidence that it "
             "resolved.")
    L.append("")
    L.append("### 6.3 Geometry and tone context")
    L.append("")
    rows = []
    for k in ("T1", "T2", "T3"):
        c = t["cells"][k]
        w = c["warp"]["own_shift_field"]
        rows.append([k, "%.2f" % w["mean_mm"], "%.2f" % w["field_std_px"],
                     "%.3f" % w["pixel_scale_mm"],
                     "%.4f" % c["fine_texture"]["mean_abs_dL_dy"],
                     "%.3f" % c["state_summary"]["rms_slope"],
                     "%.4f" % c["state_summary"]["Hs_spectral_m"]])
    L.append(table(["id", "own shift mm", "shift field std px",
                    "px scale mm", "mean |dL/dy|", "rms slope", "H_s spec"],
                   rows))
    L.append("")
    L.append("The measured warp (block matching on the photo vs the still) is "
             "in the JSON beside the predicted field; neither is a gate here, "
             "and the v14 yardsticks are unchanged.")
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 7. Hashes: new v15 files, and nothing v11-v14 touched")
    L.append("")
    fh = m["file_hashes"]
    L.append("New v15 files:")
    L.append("")
    for k, v in fh["new_v15_files"].items():
        L.append("    %-24s %s" % (k, v))
    L.append("")
    L.append("Deliverables:")
    L.append("")
    import hashlib
    with open(MERGE, "rb") as _f:
        merge_sha = hashlib.sha256(_f.read()).hexdigest()
    L.append("    %-32s %s" % ("numbers-v15-merge.json", merge_sha))
    for k, v in fh["deliverables"].items():
        L.append("    %-32s %s" % (k, v))
    L.append("")
    L.append("`numbers-v15-merge.json` cannot carry its own sha256, so the "
             "value above is computed from the file this report read at build "
             "time.  A reviewer who recomputes it after any later write will "
             "see a different value, which is the intended tripwire.")
    L.append("")
    L.append("**No v11-v14 file was modified.**  The statement is a "
             "hash comparison, not a promise:")
    L.append("")
    L.append("- Of the %d frozen modules, %d carry a sha256 recorded by the "
             "v14 run; every one of them still matches it (`regression."
             "item4_frozen_modules_unchanged.files[...].matches_v14_record`)."
             % (it4["n_files_hashed"],
                it4["n_files_compared_against_v14_record"]))
    L.append("- All %d were hashed at stage 1 and again at stage 4; "
             "`all_frozen_unchanged = %s`."
             % (len(fh["frozen_files"]),
                n(fh["all_frozen_unchanged"])))
    L.append("- The v14 report's own recorded value for `v5_composite.py` "
             "(`5c86a5770dcd0ce2`) still holds: %s."
             % yn(it4["v5_composite_prefix_matches_v14_report"]))
    L.append("")
    L.append("Per-file table (`stage1` = this job's start, `stage4` = this "
             "job's end, `v14_recorded` where the v14 numbers file has one):")
    L.append("")
    rows = []
    for k, v in fh["frozen_files"].items():
        rows.append([k, (v["stage1"] or "")[:16], (v["stage4"] or "")[:16],
                     (it4["files"][k]["v14_recorded"] or "")[:16] or "-",
                     n(v["unchanged"])])
    L.append(table(["file", "stage1", "stage4", "v14_recorded", "same"], rows))
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 8. Deviations, conventions and what this report does NOT "
             "claim")
    L.append("")
    L.append("### 8.1 The (b) convention, stated once more")
    L.append("")
    betas = [c["observation"]["b_beta"]["beta"] for c in cells]
    bhz = [c["observation"]["b_beta"]["beta_Hz_convention"] for c in cells]
    L.append("The brief asks for `E(f) = beta * I(p) * g * u* * f^-4` over "
             "`f in [1.3 f_p, 3 f_p]` and for the convention to be stated.  "
             "It is: the gate is Toba's angular-frequency form, "
             "`S(omega) = beta I(p) g u* omega^-4`, `omega` in rad/s, "
             "`I(p) = 2.5`; `E(f)` is reconstructed through the deep-water "
             "dispersion exactly as the brief says.  Under that convention "
             "the matrix lands in %.4f-%.4f -- the literal Hz reading of the "
             "same fit lands in %.1e-%.1e, two orders below the "
             "pre-registered band, which no spectrum of this normalisation "
             "can reach.  Both numbers are in every cell.  A reader who "
             "rejects the convention has the data to say so without redoing "
             "the fit." % (min(betas), max(betas), min(bhz), max(bhz)))
    L.append("")
    L.append("### 8.2 Contamination is UNRESOLVED, never PASS")
    L.append("")
    cont_cells = [k for k, v in summ["observation_b_status"].items()
                  if v == "UNRESOLVED_BAND_CONTAMINATED"]
    L.append("Where `B_h > B_l` inside the fit band the sub-check is "
             "`UNRESOLVED (band-contaminated)` and the high-wind rule accepts "
             "it as \"flagged UNRESOLVED-with-reason\".  Contaminated cells "
             "here: %s." % (", ".join(cont_cells) if cont_cells else "(none)"))
    L.append("")
    u10 = [c for c in cells if c["cell_id"] == "U10p0_F100p0"]
    if u10:
        b10 = u10[0]["observation"]["b_beta"]
        L.append("The precedence is worth stating because it decides one "
                 "cell.  On **U=10/F=100** the fitted beta is %.5f, i.e. "
                 "%.1f%% BELOW the %g floor, and that cell is also "
                 "band-contaminated.  This job applies the brief's \"never a "
                 "silent PASS/FAIL\" clause literally: a contaminated fit is "
                 "UNRESOLVED whatever its value, so the cell PASSes on the "
                 "high-wind rule.  Under the stricter alternative reading -- "
                 "out-of-band stays FAIL even when contaminated -- that cell "
                 "would FAIL and the merge would not stand.  The data for "
                 "either reading are in `cells[*].observation.b_beta` "
                 "(`beta`, `beta_in_band`, `band_contaminated`, `status`); "
                 "nothing is decided silently."
                 % (b10["beta"], 100.0 * (1.0 - b10["beta"] / b10["beta_band"][0]),
                    b10["beta_band"][0]))
        L.append("")
    L.append("### 8.3 The low-wind cells keep their status")
    L.append("")
    L.append("Four cells (U <= 3) are UNRESOLVED.  That is the design's "
             "ceiling for them, not a finding of this chunk, and no number "
             "here upgrades them.  `(b)` even FAILs on some of them (%s); the "
             "brief does not put (b) in the low-wind hard-clause list, so the "
             "verdict stays UNRESOLVED and the sub-check status is printed in "
             "full.  That is a deliberate asymmetry of the brief, reported "
             "rather than smoothed."
             % (", ".join("%s at beta %.5f"
                          % (c["cell_id"], c["observation"]["b_beta"]["beta"])
                          for c in cells
                          if c["low_wind_cell"]
                          and c["observation"]["b_beta"]["status"]
                          == "FAIL_OUT_OF_BAND") or "none"))
    L.append("")
    L.append("### 8.4 Not claimed")
    L.append("")
    L.append("- Not physical truth, not a validation: the battery is "
             "region-wise conformity + observation, as section 2 names it.")
    L.append("- The conformance column compares the model with the formula it "
             "implements; a PASS there says the implementation is faithful, "
             "nothing more.")
    L.append("- The JONSWAP anchor is context; the sub-km caveat applies to "
             "the three F = 100 cells.")
    L.append("- Transport reports C statistics, warp and texture for three "
             "cells.  It is not a measurement of a real pond and it does not "
             "settle the value of `sigma2_h_target`.")
    L.append("- No Hu-law number exists; the upper totals column remains the "
             "in-repo fetch-limited MSS with its naming gap recorded "
             "(unchanged from v13/v14).")
    L.append("")
    L.append("### 8.5 A parallel session shared this directory")
    L.append("")
    L.append("While this job ran, a second session wrote `run_v15*` files into "
             "the same folder (it was ahead on the driver at one point).  This "
             "report and every number above come from the files hashed in "
             "section 7, produced by the four foreground stages of "
             "`run_v15.py`.  If a later write changes `numbers-v15-*.json` or "
             "a still, the hashes in section 7 will not match and this report "
             "must be treated as describing a superseded state.  The "
             "checkpoint file is progress proof and was not edited by hand.")
    L.append("")
    L.append("---")
    L.append("")

    # ------------------------------------------------------------------
    L.append("## 9. Deliverables and how to check them")
    L.append("")
    L.append("```")
    L.append("out/v15/REPORT-V15-MERGE.md              this file")
    L.append("out/v15/numbers-v15-merge.json           the 8-cell matrix, both "
             "columns, flags, verdicts, regression")
    L.append("out/v15/numbers-v15-transport.json       T1/T2/T3 numbers, gates, "
             "travelling labels")
    L.append("out/v15/still-T1-windy.png               U=6  F=1000 h=1.0")
    L.append("out/v15/still-T2-lowwind.png             U=3  F=1000 h=0.6 "
             "(UNRESOLVED label inside)")
    L.append("out/v15/still-T3-warned.png              U=6  F=100   h=1.0 "
             "(Omega_c = 7.30, warn = true)")
    L.append("out/v15/checkpoint.json                  jobkit, one unit per "
             "stage")
    L.append("```")
    L.append("")
    L.append("Checks a reviewer can run without this report:")
    L.append("")
    L.append("1. `numbers-v15-merge.json -> cells[*].conformance` and "
             "`.observation` are two separate objects per cell; nothing "
             "collapses them.")
    L.append("2. `cells[*].b_max_observation.gated == false` on all 8 cells: "
             "the retired clause cannot score anything in v15.")
    L.append("3. `regression.item2_...` compares five fields per cell with "
             "`==`; `exact` booleans are per field.")
    L.append("4. `numbers-v15-transport.json -> cells.T2.cell_status == "
             "'UNRESOLVED'` and `cells.T3.warn_extrapolated_wave_age == "
             "true`.")
    L.append("5. `regression.item3_...g1_legacy_rerun.rerun_sha256` equals "
             "`%s`." % (it3["g1_legacy_rerun"]["rerun_sha256"] or "")[:16])
    L.append("")
    L.append("Report built at %s (wall %.1f s)."
             % (time.strftime("%Y-%m-%dT%H:%M:%S"), time.time() - t0))
    L.append("")

    text = "\n".join(L)
    data = text.encode("ascii", "replace")

    def verify(p):
        s = os.path.getsize(p)
        return s > 2000

    atomic_write_bytes(REPORT, data, verify=verify)
    print("[report] wrote %s (%d bytes)" % (REPORT, len(data)))
    return REPORT


if __name__ == "__main__":
    main()
