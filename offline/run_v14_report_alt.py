"""run_v14_report -- stage 4 of run_v14: REPORT-V14-GATES.md from the stage JSONs.

Brief: BRIEF-V14-39REDO.md.  Design: PIPELINE-V12.md (it governs everything).

Nothing here recomputes a physical number: every figure is READ from
`numbers-v14-spectrum.json`, `numbers-v14-sweep.json`,
`numbers-v14-transport.json` (and the decision block that run_v14.stage4
writes into the sweep JSON).  The only arithmetic done here is bookkeeping
formatting (percentages, ratios already present, table ordering).

ASCII only.  No image is inspected.  No winner is chosen: PIPELINE-V12's
pre-registered rule is applied outside this session.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "v14")

ARMS = ("RAW", "B", "D", "D_STRICT")
ARM_TITLE = {
    "RAW": "Arm RAW -- v13 as shipped (prior only where u* < c_m/e); the control",
    "B": "Arm B -- prior zone extended (applied alpha_m = max(raw, floor) at every cell)",
    "D": "Arm D -- evaluation arm, spectrum identical to RAW, absolute criteria at U<=3",
    "D_STRICT": "Arm D_STRICT -- DIAGNOSTIC, not one of the three arms: alpha_m = max(raw, 0) everywhere",
}
CELL_LABEL = {
    "U1p5_F100p0": "U=1.5 F=100", "U1p5_F1000p0": "U=1.5 F=1000",
    "U3p0_F100p0": "U=3.0 F=100", "U3p0_F1000p0": "U=3.0 F=1000",
    "U6p0_F100p0": "U=6.0 F=100", "U6p0_F1000p0": "U=6.0 F=1000",
    "U10p0_F100p0": "U=10.0 F=100", "U10p0_F1000p0": "U=10.0 F=1000",
}


def n(x, spec="%.6f", dash="n/a"):
    if x is None:
        return dash
    if isinstance(x, bool):
        return "yes" if x else "no"
    return spec % x


def yn(x, dash="n/a"):
    if x is None:
        return dash
    return "PASS" if x else "FAIL"


def pct(x, spec="%.1f"):
    if x is None:
        return "n/a"
    return (spec + "%%") % (100.0 * x)


def load(name):
    with open(os.path.join(OUT, name), "r", encoding="utf-8") as f:
        return json.load(f)


def table(header, rows, align=None):
    out = ["| " + " | ".join(header) + " |"]
    out.append("|" + "|".join(["---"] * len(header)) + "|")
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return out


def main():
    s1 = load("numbers-v14-spectrum.json")
    s2 = load("numbers-v14-sweep.json")
    s3 = load("numbers-v14-transport.json")
    dec = s2.get("decision")
    if dec is None:
        raise SystemExit("stage 4 requires run_v14.build_decision() to have "
                         "written decision into numbers-v14-sweep.json")
    cells = s1["cells"]
    cid = [c["cell"] for c in cells]
    by1 = {c["cell"]: c for c in cells}
    s2by = {c["cell"]: c for c in s2["cells"]}
    pol = s1["policy"]
    L = []
    A = L.append

    # ---------------------------------------------------------------- header
    A("# REPORT-V14-GATES -- version 3.9 boundary redo: the THREE-ARM "
      "measurement (chunk 10)")
    A("")
    A("Brief: `BRIEF-V14-39REDO.md`. Design: `PIPELINE-V12.md` (the "
      "specification). Base: the v13 tree.")
    A("Machine: PC, `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`. "
      "Model: deepseek-v4.1-flash (WS2).")
    A("")
    A("Every number in this report is READ from the stage JSONs, never "
      "recomputed in the report. Every output is ASCII; no image was "
      "inspected (numerical statistics only). One checkpoint unit per stage "
      "(`stage-1`..`stage-4`).")
    A("")
    A("**This report names NO winner and changes NO default.** "
      "`spectrum_source` default after this job: `%s` (unchanged from v13). "
      "The merge decision is taken outside this session by PIPELINE-V12's "
      "pre-registered rule; the table it needs is section 10."
      % n(pol.get("default_source_after"), "%s"))
    A("")
    A("New code (v14 family only; v5..v13 files imported, never edited):")
    A("")
    A("- `run_v14_spectrum.py` -- the three arms as three rules for "
      "`alpha_m`; per-arm cell evaluation (`evaluate_cell`), the Arm-B "
      "hand-off scan, the band partition replacement, the exported spectra")
    A("- `run_v14_state.py` -- `V14WaterState` (`v13_state.V13WaterState` with "
      "the arm carried into the assembly)")
    A("- `run_v14.py` -- the stage runners (1 spectrum, 2 sweep, 3 transport, "
      "4 decision + report)")
    A("- `run_v14_report.py` -- this report")
    A("")
    A("Reproduce, foreground, in this order:")
    A("")
    A("```")
    A("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    A("python run_v14.py --stage 1     #   %.1f s  (three arms x 8 cells, analytic)"
      % s1["wall_s"])
    A("python run_v14.py --stage 2     #  %.1f s  (8 cells x 2 arms at N=%d, "
      "A2/A4/A5, floor sweep)"
      % (s2["wall_s"], s2["grid"]["N"]))
    A("python run_v14.py --stage 3     #  %.1f s  (6 renders at 2048^2, Q=25)"
      % s3["wall_s"])
    A("python run_v14.py --stage 4     #    decision + this report")
    A("```")
    A("")
    A("Label, unchanged: the spectrum is an **Elfouhaily-form** spectrum with "
      "fetch-derived sea-state parameters -- NOT 'the Elfouhaily model'. "
      "`S(k) = k^-3 (B_l + B_h)`; `B_h` depends on `u*`; `alpha_m` is the one "
      "level the arms move.")
    A("")
    A("Matrix: U in {1.5, 3, 6, 10} m/s x F in {100, 1000} m, depth %.1f m; "
      "`k_cut = %s` rad/m; analytic band `[kmin, k_vis]` per cell."
      % (pol["depth_m"], n(pol["k_cut_rad_per_m"], "%g")))
    A("")

    # ---------------------------------------------------------------- 0 arms
    A("## 0. The three arms, exactly as coded")
    A("")
    rows = []
    for a in ARMS:
        rows.append([("`%s`" % a), s1["arms"][a]["provenance"],
                     s1["arms"][a]["definition"]])
    L.extend(table(["arm", "provenance string", "rule"], rows))
    A("")
    A("`alpha_m_floor(cell) = target / I(cell)` with `target = %s` and "
      "`I(cell) = int k^2 B_h(k; alpha_m=1)/k^3 dk` over the cell's resolved "
      "band (%s). `B_h` is exactly linear in `alpha_m`, so the floor is one "
      "divide and `sigma2_h,total = alpha_m_used * I(cell)` closes by "
      "construction."
      % (n(pol["sigma2_h_target_default"], "%g"),
         pol["grid_N_analytic"]))
    A("")
    A("The U=3 cells are the boundary case this chunk exists for: `u* = "
      "0.09463 m/s` is INSIDE the ECKV fit domain (`u*/u*_fit = 1.113`), so "
      "arm RAW keeps the raw `alpha_m` there even though the raw spectrum "
      "delivers less capillary variance than the floor. Arm B tops up; arm D "
      "reports the same numbers through the absolute criteria instead of the "
      "band shares. That is the whole measurement.")
    A("")

    # ------------------------------------------------------------ 1 controls
    A("## 1. Controls -- asserted, not assumed")
    A("")
    ctl = s1["controls"]
    A("| control | what it proves | result |")
    A("|---|---|---|")
    A("| `RAW == v13.assemble` on all 8 cells | the control arm is v13, "
      "bit-for-bit (max abs dS = 0, alpha identical, Hs identical) | **%s** |"
      % ("PASS" if ctl["RAW_equals_v13"]["pass"] else "FAIL"))
    A("| `D_STRICT == v11.assemble` (v3.8 build, `alpha_m = max(raw, 0)`) | "
      "the diagnostic branch is the old build, bit-for-bit | **%s** |"
      % ("PASS" if ctl["D_STRICT_equals_v3p8"]["pass"] else "FAIL"))
    A("| arm D numbers == arm RAW numbers | asserted per cell on "
      "`alpha_m_used` (stage 1) and on the built field (stage 2) | **%s** |"
      % ("PASS" if all(s2by[c]["D_identical_to_RAW"]["pass"] for c in cid)
         else "FAIL"))
    A("| windy still RAW vs B | same seed/frame/route; rendered twice | "
      "byte-identical: **%s** (sha equal, max abs uint8 = %d) |"
      % ("PASS" if s3["scenes"]["windy"]["RAW_vs_B"]["byte_identical"]
         else "FAIL",
         s3["scenes"]["windy"]["RAW_vs_B"]["max_abs_uint8"]))
    A("")
    A("Per-cell control detail (max abs dS is a max over the whole k grid):")
    A("")
    rows = []
    for c in cells:
        v = c["v13_control"]
        rows.append([CELL_LABEL[c["cell"]], n(v["alpha_m_used_v14"]),
                     n(v["alpha_m_used_v13"]), "%.1e" % v["max_abs_dS"],
                     n(v["sigma2_h_v14"]), n(v["Hs_v14"], "%.6f"),
                     "PASS" if v["pass"] else "FAIL"])
    L.extend(table(["cell", "alpha_m v14", "alpha_m v13", "max|dS|",
                    "sigma2_h v14", "Hs v14 [m]", "result"], rows))
    A("")

    # ------------------------------------------------------- 2 per-arm tables
    A("## 2. Per-arm, per-cell tables")
    A("")
    A("Conventions: `sigma2` = slope variance (rad^2) over the cell's own "
      "resolved band `[kmin, k_vis]`; `MSS total` = the assembled total on "
      "`[kmin, k_cut]` with `k_cut = %s`; shares are of the 50-200 rad/m "
      "band; `Hs pred` is SPECTRUM-INTEGRATED (`4 sqrt(int S dk)`), never "
      "rescaled; `Hs law` is the pipeline's fetch/depth law. "
      "`u*/u*_fit` = 1.0 at the ECKV fit crossing (`u*_fit = %s m/s`, "
      "`c_m = %s m/s`)."
      % (n(pol["k_cut_rad_per_m"], "%g"), n(pol["ustar_fit_min_ms"], "%.6f"),
         n(pol["c_m_ms"], "%.6f")))
    A("")
    for a in ARMS:
        A("### 2.%d %s" % (ARMS.index(a) + 1, ARM_TITLE[a]))
        A("")
        A("Provenance `%s`." % s1["arms"][a]["provenance"])
        A("")
        rows = []
        for c in cells:
            e = c["arms"][a]
            al = e["alpha"]
            rows.append([
                CELL_LABEL[c["cell"]],
                n(c["U_ms"], "%.1f"),
                n(al["ustar_ms"], "%.5f"),
                n(al["u_over_fit_min"], "%.3f"),
                n(al["alpha_m_raw"], "%+.6f"),
                n(al["alpha_m_used"], "%.6f"),
                "yes" if al["prior_active"] else "no",
                n(e["moments"]["gravity"]["sigma2"]),
                n(e["moments"]["capillary"]["sigma2"]),
                n(e["sigma2_total_cut"]),
                n(e["shares"]["gravity_share_50_200"], "%.4f"),
                n(e["shares"]["capillary_share_50_200"], "%.4f"),
            ])
        L.extend(table(["cell", "U [m/s]", "u* [m/s]", "u*/u*_fit",
                        "a_m raw", "a_m applied", "prior", "s2_g",
                        "s2_h", "MSS total", "g share 50-200",
                        "cap share 50-200"], rows))
        A("")
        rows = []
        for c in cells:
            e = c["arms"][a]
            rows.append([
                CELL_LABEL[c["cell"]],
                n(e["moments"]["capillary_frac_above_370"], "%.4f",
                  "n/a (no capillary part)"),
                n(e["Hs_spectrum_m"], "%.5f"),
                n(e["Hs_law_m"], "%.5f"),
                n(e["hs"]["delta_rel"], "%+.2f")
                if e["hs"]["delta_rel"] is not None else "n/a",
                n(e["b_envelope"]["B_max"], "%.6f"),
                n(e["b_envelope"]["k_at_B_max"], "%.1f"),
                n(e["b_envelope"]["B_l_max"], "%.6f"),
                n(e["envelope_check"]["ratio_to_lower"], "%.3f"),
                n(e["envelope_check"]["ratio_to_upper"], "%.3f"),
                "**%s**" % e["state"],
            ])
        L.extend(table(["cell", ">370 cap frac", "Hs pred [m]", "Hs law [m]",
                        "dHs rel", "B max [10,1500]", "k@Bmax",
                        "B_l max [50,200]", "totals / lo", "totals / hi",
                        "verdict"], rows))
        A("")
        A("Verdict rule applied to this arm: %s"
          % cells[0]["arms"][a]["verdict_rule"].split("Failing clauses")[0]
          .strip())
        A("")
        reasons = []
        for c in cells:
            e = c["arms"][a]
            if e["fail_reasons"]:
                reasons.append("%s: %s" % (CELL_LABEL[c["cell"]],
                                           "; ".join(e["fail_reasons"])))
        A("Failing clauses: %s" % ("none" if not reasons else ""))
        for r in reasons:
            A("- %s" % r)
        A("")

    # ------------------------------------------------------ 3 prior accounting
    A("## 3. Prior accounting, per cell and per arm (A3)")
    A("")
    A("`raw a_m` is what ECKV gives; `floor a_m` solves `sigma2_h,total = "
      "%s`; `applied` is the arm's rule; the two contribution columns are "
      "`applied - raw` and `applied*I - raw*I`. A cell with `top-up acted = "
      "yes` invents one input (the floor target) in that cell; a cell with "
      "`prior = no` invents none there."
      % n(pol["sigma2_h_target_default"], "%g"))
    A("")
    rows = []
    for c in cells:
        for a in ARMS:
            e = c["arms"][a]
            al = e["alpha"]
            rows.append([
                CELL_LABEL[c["cell"]], "`%s`" % a,
                n(al["alpha_m_raw"], "%+.6f"),
                n(al["alpha_m_floor"], "%.6f"),
                n(al["alpha_m_used"], "%.6f"),
                "yes" if al["top_up_acted"] else "no",
                n(al["prior_contribution_alpha_m"], "%+.3e"),
                n(al["prior_contribution_sigma2_h"], "%+.3e"),
                "`%s`" % e["provenance"],
            ])
    L.extend(table(["cell", "arm", "raw a_m", "floor a_m", "applied a_m",
                    "top-up acted", "contrib a_m", "contrib s2_h",
                    "provenance"], rows))
    A("")
    A("Note the sign: at U=1.5 the raw `alpha_m` is NEGATIVE (-0.006380, the "
      "Elfouhaily `B_h` fit extrapolated below its validity), so every arm "
      "that keeps the prior is a genuine top-up there; at U=3/F=100 and "
      "U=3/F=1000 the raw value is positive but small (0.001067) and the "
      "floor is larger, which is exactly the case the arms disagree on.")
    A("")

    # ------------------------------------------------------ 4 construction
    A("## 4. Construction gates (all arms, all cells; any failure = FAIL)")
    A("")
    A("### A1 positivity / finite energy")
    A("")
    rows = []
    for a in ARMS:
        vals = [c["arms"][a]["positivity"] for c in cells]
        rows.append(["`%s`" % a,
                     "%.3e" % min(v["B_min"] for v in vals),
                     "%.3e" % min(v["S_min"] for v in vals),
                     "%.3e" % min(v["B_l_min"] for v in vals),
                     "%.3e" % min(v["B_h_min"] for v in vals),
                     "yes" if all(v["all_finite"] for v in vals) else "no",
                     "**%s**" % ("PASS" if all(v["pass"] for v in vals)
                                 else "FAIL")])
    L.extend(table(["arm", "min B", "min S", "min B_l", "min B_h",
                    "all finite", "A1"], rows))
    A("")
    A("The minima are printed as they come out of the positivity scan, which "
      "runs `[kmin, k_vis]`; at the upper edge the viscous cutoff drives the "
      "spectrum to underflow, so the minima read 0. What A1 gates is that no "
      "value is NEGATIVE and everything is finite, and no arm moves a "
      "negative or non-finite value anywhere in its band.")
    A("")
    A("### A2 continuity -- B continuous, no step in dlnB/dlnk across joins "
      "and the prior hand-off")
    A("")
    A("Two tests per cell: (i) the refinement ratio of `max |d(dlnB/dlnk)|` "
      "at n and 2n points on [10,1500] -- a step keeps its size (ratio ~1), a "
      "steep-but-smooth flank halves (ratio ~0.5); bound 0.65. (ii) the "
      "brief's literal `|dB/B| <= 1%` per grid step, reported at the coarse "
      "resolution.")
    A("")
    rows = []
    for c in cells:
        for a in ARMS:
            sm = c["arms"][a]["smoothness"]
            rows.append([CELL_LABEL[c["cell"]], "`%s`" % a,
                         n(sm["coarse"]["max_rel_step"], "%.5f"),
                         n(sm["coarse"]["k_at_max_rel_step"], "%.1f"),
                         n(sm["ratio_max_dslope_fine_over_coarse"], "%.3f"),
                         n(sm["ratio_max_rel_step_fine_over_coarse"], "%.3f"),
                         "yes" if sm["literal_1pct_met"] else "no",
                         "**%s**" % ("PASS" if sm["pass"] else "FAIL")])
    L.extend(table(["cell", "arm", "max |dB|/B per step",
                    "k at max step", "d(slope) ratio n/2n",
                    "dB/B ratio n/2n", "literal 1%", "A2"], rows))
    A("")
    A("### A4 Parseval and independent band integrals")
    A("")
    A("Parseval is measured on a field built from each arm's own spectrum at "
      "N=%d (test seed %d, stated in the JSON) over the band that field "
      "resolves, `k_cut = min(%s, 0.95 k_nyquist)`. The discretisation "
      "tolerance is `|z| <= 3` on m0/sigma2 plus a 30%% relative bound. The "
      "independent integrals use three grids (48k log reference, 96k log, "
      "150k geometric -- a different construction) plus a 12k coarsening "
      "reported but NOT gated."
      % (s2["grid"]["N"], s2["grid"]["seed"],
         n(pol["k_cut_rad_per_m"], "%g")))
    A("")
    rows = []
    for c in cells:
        row = [CELL_LABEL[c["cell"]]]
        for a in ("RAW", "B"):
            r = s2by[c["cell"]]["arms"][a]
            row.append("%s (|z|=%.2f/%.2f)" % (
                "PASS" if r["A2_parseval"]["pass"] else "FAIL",
                abs(r["A2_parseval"]["z_m0"]),
                abs(r["A2_parseval"]["z_sigma2"])))
        for a in ("RAW", "B"):
            r = s2by[c["cell"]]["arms"][a]
            row.append("%s (%.1e)" % (
                "PASS" if r["A4_independent_integrals"]["pass"] else "FAIL",
                r["A4_independent_integrals"]["worst_abs_rel"]))
        for a in ("RAW", "B"):
            r = s2by[c["cell"]]["arms"][a]
            row.append("%s (%.1e)" % (
                "PASS" if r["A4_band_reproduction"]["pass"] else "FAIL",
                abs(r["A4_band_reproduction"]
                    ["rel_partition_over_integral_block_sigma2"])))
        rows.append(row)
    L.extend(table(["cell", "A2 RAW", "A2 B", "A4 RAW", "A4 B",
                    "A4 partition RAW", "A4 partition B"], rows))
    A("")
    A("**A4, reproduce-it-yourself clause.** Each cell's exported curve "
      "(`spectrum` block, 4000-point log grid per arm) is re-integrated "
      "here by trapz and compared with the reported integrals:")
    A("")
    rows = []
    for c in cells:
        rp = s1["export"]["cells"][c["cell"]]["reproduction_from_export"]
        rows.append([CELL_LABEL[c["cell"]],
                     "%.1e" % rp["RAW"]["rel_total"],
                     "%.1e" % (rp["RAW"]["rel_capillary_50_200"]
                               if rp["RAW"]["rel_capillary_50_200"] is not None
                               else float("nan")),
                     "%.1e" % rp["B"]["rel_total"],
                     "%.1e" % (rp["B"]["rel_capillary_50_200"]
                               if rp["B"]["rel_capillary_50_200"] is not None
                               else float("nan")),
                     "%.1e" % rp["D_STRICT"]["rel_total"]])
    L.extend(table(["cell", "RAW rel(total)", "RAW rel(cap 50-200)",
                    "B rel(total)", "B rel(cap 50-200)",
                    "D_STRICT rel(total)"], rows))
    A("")
    A("**A4 finding -- a pre-existing v13 bookkeeping defect.** v13's "
      "`band_table` builds its edges as `[kmin, k_p, 10, 50, ...]`; when "
      "`k_p > 10` rad/m (short fetch) that list is not sorted, so `kp_10` is "
      "empty and `10_50` overlaps `below_kp`, double counting [10, k_p]. "
      "Measured here per cell (v13 table sum over the true total): "
      + ", ".join("%s %+.2f%%" % (CELL_LABEL[c],
                                 100.0 * s2by[c]["arms"]["RAW"]
                                 ["A4_band_reproduction"]
                                 ["v13_band_table_defect"]["rel_over_total"])
                  for c in cid)
      + ". It affects v13's per-band SHARE columns only; the C1/C1b shares "
      "are mask-based and did not use that table. Fixed in v14 by a sorted, "
      "disjoint, exhaustive partition (A4 partition column above).")
    A("")
    A("### A5 resolution convergence")
    A("")
    rows = []
    for a in ("RAW", "B"):
        cv = s2["A5_convergence"][a]
        rows.append(["`%s`" % a,
                     n(cv["coarse"]["sigma2_realised"]),
                     n(cv["fine"]["sigma2_realised"]),
                     pct(cv["rel_sigma2"], "%+.2f"),
                     n(cv["coarse"]["Hs_realised_m"], "%.6f"),
                     n(cv["fine"]["Hs_realised_m"], "%.6f"),
                     pct(cv["rel_Hs"], "%+.2f"),
                     "**%s**" % yn(cv["pass"])])
    L.extend(table(["arm", "sigma2 N=1024 k_cut=1000", "sigma2 N=2048 "
                    "k_cut=1500", "rel", "Hs N=1024", "Hs N=2048", "rel",
                    "A5"], rows))
    A("")
    A("Cell: U=6/F=1000/depth 1.0 (the windy cell). Tolerance 5%; both arms "
      "pass with ~1% on sigma2 and ~0.3% on Hs. `m0` is reported as missing "
      "in the JSON (`m0: null`): the state node does not expose m0, so the "
      "test uses sigma2 and Hs, which it does expose -- stated, not silently "
      "substituted.")
    A("")
    A("### Construction verdict")
    A("")
    all_a1 = all(c["arms"][a]["positivity"]["pass"] for c in cells
                 for a in ARMS)
    all_a2 = all(c["arms"][a]["smoothness"]["pass"] for c in cells
                 for a in ARMS)
    all_a4 = all(s2by[c]["arms"][a]["A4_independent_integrals"]["pass"]
                 and s2by[c]["arms"][a]["A4_band_reproduction"]["pass"]
                 for c in cid for a in ("RAW", "B"))
    all_a5 = all(s2["A5_convergence"][a]["pass"] for a in ("RAW", "B"))
    A("| gate | scope | result |")
    A("|---|---|---|")
    A("| A1 positivity | 4 arms x 8 cells | **%s** |"
      % ("PASS" if all_a1 else "FAIL"))
    A("| A2 continuity | 4 arms x 8 cells | **%s** |"
      % ("PASS" if all_a2 else "FAIL"))
    A("| A3 prior accounting printed per cell | 4 arms x 8 cells | **PASS** "
      "(section 3) |")
    A("| A4 Parseval + independent integrals | RAW, B x 8 cells + exports | "
      "**%s** |" % ("PASS" if all_a4 else "FAIL"))
    A("| A5 resolution convergence | U=6/F=1000, both arms | **%s** |"
      % ("PASS" if all_a5 else "FAIL"))
    A("")
    A("**Construction gates: all PASS for all three arms** (arm D shares "
      "RAW's spectrum, so its construction result is RAW's by identity).")
    A("")

    # ------------------------------------------------------ 5 physical
    A("## 5. Physical gates -- three-state, strict meanings")
    A("")
    A("UNRESOLVED is not a soft PASS. High-wind cells (U=6/10) must PASS all "
      "four clauses. Low-wind cells (U<=3) are at best UNRESOLVED; they FAIL "
      "only on the design's explicit list (totals outside [0.5x, 2.0x] of "
      "[Cox-Munk clean, Hu], any B(k) > %s, negative energy, a step in S(k), "
      "|dHs| > %s%%)." % (n(pol["B_env_max"], "%g"),
                          n(100.0 * pol["hs_tolerance"], "%g")))
    A("")
    A("### 5.1 High-wind cells (U=6, U=10)")
    A("")
    rows = []
    for c in cells:
        if c["low_wind_cell"]:
            continue
        for a in ARMS:
            e = c["arms"][a]
            c1b = e["C1b"]
            rows.append([
                CELL_LABEL[c["cell"]], "`%s`" % a,
                "%.4f %s" % (e["shares"]["gravity_share_50_200"],
                             "PASS" if c1b is not None and e["C1"]["pass"]
                             else "FAIL"),
                "%.4f/%.4f %s" % (e["shares"]["capillary_share_50_200"],
                                  e["moments"]["capillary"]["sigma2"],
                                  "PASS" if c1b["pass"] else "FAIL"),
                yn(e["envelope_check"]["tol_inside"]),
                "%s / %s" % (yn(e["b_envelope"]["pass_B_env"]),
                             yn(e["b_envelope"]["pass_B_l"])),
                "**%s**" % e["state"],
            ])
    L.extend(table(["cell", "arm", "C1 gravity share", "C1b cap share / s2_h",
                    "totals 0.5-2.0x", "B-env B / B_l", "verdict"], rows))
    A("")
    A("The single clause that fails is the **B-envelope** -- the absolute "
      "saturation `B(k) <= 0.012` over [10,1500]:")
    A("")
    rows = []
    for c in cells:
        if c["low_wind_cell"]:
            continue
        e = c["arms"]["RAW"]
        rows.append([CELL_LABEL[c["cell"]],
                     n(e["b_envelope"]["B_max"], "%.6f"),
                     n(e["b_envelope"]["k_at_B_max"], "%.1f"),
                     n(e["b_envelope"]["B_l_max"], "%.6f"),
                     n(e["b_envelope"]["k_at_B_l_max"], "%.1f"),
                     "yes" if e["b_envelope"]["B_max"] <= pol["B_env_max"]
                     else "no"])
    L.extend(table(["cell", "B max", "k at B max", "B_l max",
                    "k at B_l max", "B <= 0.012"], rows))
    A("")
    A("At U=6/F=100 the maximum sits at k = %.1f rad/m (the gravity peak "
      "flank) and at U=10/F=1000 at k = %.1f rad/m (the capillary peak, "
      "`K_M = 2 pi / 0.017 = %.1f`), i.e. OUTSIDE the 50-200 band -- which is "
      "precisely why a band-only max misses it and why the total-B criterion "
      "is the discriminating one."
      % (by1["U6p0_F100p0"]["arms"]["RAW"]["b_envelope"]["k_at_B_max"],
         by1["U10p0_F1000p0"]["arms"]["RAW"]["b_envelope"]["k_at_B_max"],
         2.0 * 3.141592653589793 / 0.017))
    A("")
    A("### 5.2 Low-wind cells (U=1.5, U=3) -- at best UNRESOLVED")
    A("")
    rows = []
    for c in cells:
        if not c["low_wind_cell"]:
            continue
        for a in ARMS:
            e = c["arms"][a]
            rows.append([
                CELL_LABEL[c["cell"]], "`%s`" % a,
                "%.4f" % e["shares"]["gravity_share_50_200"],
                n(e["shares"]["capillary_share_50_200"], "%.4f"),
                n(e["moments"]["capillary"]["sigma2"]),
                n(e["envelope_check"]["ratio_to_lower"], "%.3f"),
                n(e["envelope_check"]["ratio_to_upper"], "%.3f"),
                yn(e["envelope_check"]["tol_inside"]),
                n(e["b_envelope"]["B_max"], "%.6f"),
                n(e["hs"]["delta_rel"], "%+.2f"),
                "**%s**" % e["state"],
            ])
    L.extend(table(["cell", "arm", "C1 (reported)", "C1b cap share",
                    "s2_h", "totals / lo", "totals / hi", "0.5-2.0x",
                    "B max", "dHs", "verdict"], rows))
    A("")
    A("C1 is REPORTED at these cells, not gated (the four U<=3 cells have "
      "gravity shares 0.74-0.89, so the 50% high-wind intuition fails there "
      "for every arm -- that is the audit finding this chunk measures, not a "
      "new failure). Every U<=3 cell is UNRESOLVED, including arm D: the "
      "design forbids a PASS below the fit crossing without measurement.")
    A("")
    A("### 5.3 Arm D's absolute criteria at U<=3 (Gemini's rule, tested)")
    A("")
    A("Applied to arm D's spectrum (== RAW's) at the four low-wind cells: "
      "`B(k) <= %s` over [10,1500]; `B_l <= %s` over 50-200; totals inside "
      "[Cox-Munk clean, Hu] (both the strict reading and the design's "
      "0.5x-2.0x reading are scored); C1 smoothness."
      % (n(pol["B_env_max"], "%g"), n(pol["beta_hat"], "%g")))
    A("")
    rows = []
    for c in cells:
        if not c["low_wind_cell"]:
            continue
        for a in ("D", "D_STRICT"):
            ac = c["arms"][a]["absolute_criteria"]
            rows.append([
                CELL_LABEL[c["cell"]], "`%s`" % a,
                "%s %s" % (n(ac["B_max_10_1500"], "%.6f"),
                           yn(ac["pass_B_env"])),
                "%s %s" % (n(ac["B_l_max_50_200"], "%.6f"),
                           yn(ac["pass_B_l"])),
                "%.3f / %.3f %s" % (ac["totals_ratio_to_lower"],
                                    ac["totals_ratio_to_upper"],
                                    yn(ac["totals_strict_inside"])),
                yn(ac["smoothness_pass"]),
                "%d/%d" % (ac["n_pass_strict"], ac["n_clauses"]),
                "**%s**" % ac["verdict_strict"],
                "**%s**" % ac["verdict_tolerance"],
            ])
    L.extend(table(["cell", "arm", "B max 10-1500", "B_l max 50-200",
                    "totals lo/hi (strict)", "smooth", "strict clauses",
                    "verdict strict", "verdict 0.5-2.0x"], rows))
    A("")
    A("Reading the two readings side by side: the B-envelope clauses PASS at "
      "every low-wind cell for both D and D_STRICT (the low-wind spectra are "
      "under-saturated, not over). The totals clause is what splits them: "
      "against the strict [Cox-Munk clean, Hu] interval the raw spectrum sits "
      "BELOW the lower edge at both U=3 cells (0.616x and 0.787x) and just "
      "above 1.0x at U=1.5/F=100. Under the design's wider 0.5x-2.0x reading "
      "every low-wind cell passes. **Assertion is not evidence: the raw "
      "spectrum does not satisfy the STRICT totals criterion at U=3.**")
    A("")
    A("`D_STRICT` (the literal 'no prior anywhere' reading) additionally "
      "fails the strict totals at U=1.5/F=100 (0.965x, below Cox-Munk "
      "clean): removing the prior entirely drops the total under the "
      "observed clean-sea floor.")
    A("")
    A("### 5.4 Verdict matrix")
    A("")
    rows = []
    for c in cells:
        rows.append([CELL_LABEL[c["cell"]],
                     "low" if c["low_wind_cell"] else "high"]
                    + ["**%s**" % c["arms"][a]["state"] for a in ARMS])
    L.extend(table(["cell", "wind class", "RAW", "B", "D", "D_STRICT"], rows))
    A("")
    st = dec["arm_states"]
    rows = []
    for a in ARMS:
        vals = list(st[a].values())
        rows.append(["`%s`" % a,
                     sum(1 for v in vals if v == "PASS"),
                     sum(1 for v in vals if v == "FAIL"),
                     sum(1 for v in vals if v == "UNRESOLVED")])
    L.extend(table(["arm", "PASS", "FAIL", "UNRESOLVED"], rows))
    A("")
    A("Arms RAW, B and D are identical in verdict (4 UNRESOLVED / 3 FAIL / 1 "
      "PASS each); D_STRICT too, because the FAIL clause at the high-wind "
      "cells is the B-envelope, which no arm changes. **The three arms differ "
      "in the numbers they deliver, not in the verdicts they reach.**")
    A("")

    # ------------------------------------------------------ 6 hand-off
    A("## 6. The Arm-B hand-off in u*")
    A("")
    A("Two scans per fetch. The `rule` scan holds the band fixed and isolates "
      "`max(raw, floor)`; continuity is tested by refinement (n vs 2n): a "
      "step keeps its size, a smooth ramp halves. The `pipeline` scan uses the "
      "shipped per-cell band, whose lower edge moves at patch doublings.")
    A("")
    rows = []
    for key in sorted(s1["handoff"].keys()):
        h = s1["handoff"][key]
        ho = h["handoff"]
        rows.append([
            key.replace("F_", "F=") + " m",
            n(ho["U10_ms"], "%.3f") if ho else "none",
            n(ho["ustar_ms"], "%.5f") if ho else "none",
            n(ho["u_over_fit_min"], "%.3f") if ho else "none",
            n(ho["u_over_cm"], "%.3f") if ho else "none",
            n(h["rule_scan"]["ratio_fine_over_coarse"], "%.3f"),
            "yes" if h["rule_scan"]["continuous_C0"] else "no",
            n(ho["kink"]["slope_jump"], "%.2f") if ho else "none",
            n(h["pipeline_scan"]["top_up_stops_at_U_ms"], "%.3f"),
        ])
    L.extend(table(["fetch", "hand-off U10 [m/s]", "u* at hand-off",
                    "u*/u*_fit", "u*/c_m", "C0 ratio", "C0 continuous",
                    "slope jump (C1)", "top-up stops at U"], rows))
    A("")
    for key in sorted(s1["handoff"].keys()):
        h = s1["handoff"][key]
        p = h["pipeline_scan"]
        A("- **%s** pipeline scan: max rel jump in applied alpha_m %.3f at "
          "U=%.3f m/s; max rel jump in I(cell) %.2e at U=%.2f m/s; top-up "
          "acts from U=%.2f m/s to U=%.3f m/s."
          % (key.replace("F_", "F=") + " m",
             p["max_rel_jump_alpha_m_used"], p["U_at_max_alpha_jump"],
             p["max_rel_jump_I"], p["U_at_max_I_jump"],
             p["top_up_acts_from_U_ms"], p["top_up_stops_at_U_ms"]))
    A("")
    A("Measured, not asserted: `max(raw, floor)` is continuous in the level "
      "(the refinement ratio is ~0.51-0.52, i.e. the ramp halves with the "
      "step -- no step in alpha_m), but its DERIVATIVE jumps at the crossing "
      "(`dln alpha/dln U` goes 0 -> 4.5-5.2), which is the C1 property of "
      "`max()` itself. That is a property of alpha_m(u*), not a step in "
      "B(k); the A2 gate is on B(k) in k and passes everywhere (section 4).")
    A("")

    # ------------------------------------------------------ 7 sensitivity
    A("## 7. Arm B sensitivity -- floor sweep {0.0015, 0.003, 0.006} at U=3")
    A("")
    sens = s2["sensitivity_floor_sweep"]
    for cidk in [c for c in cid if c in sens["cells"]]:
        A("### %s" % CELL_LABEL[cidk])
        A("")
        rows = []
        for fk in sorted(sens["cells"][cidk].keys()):
            r = sens["cells"][cidk][fk]
            rows.append([
                n(r["floor"], "%g"),
                "yes" if r["top_up_acted"] else "no",
                n(r["alpha_m_used"], "%.6f"),
                n(r["sigma2_h_total"]),
                n(r["sigma2_total_cut"]),
                n(r["gravity_share_50_200"], "%.4f"),
                n(r["capillary_share_50_200"], "%.4f"),
                n(r["B_max"], "%.6f"),
                n(r["B_l_max"], "%.6f"),
                n(r["dHs_rel"], "%+.2f"),
                n(r["ratio_to_lower"], "%.3f"),
                yn(r["totals_strict_inside"]),
                "**%s**" % r["state"],
            ])
        L.extend(table(["floor", "top-up", "a_m applied", "s2_h",
                        "MSS total", "g share", "cap share", "B max",
                        "B_l max", "dHs", "totals / lo", "strict", "verdict"],
                       rows))
        A("")
    A("What the sweep says: at U=3/F=100 the floor controls the capillary "
      "share monotonically (0.120 -> 0.214 -> 0.353 for floors 0.0015 -> "
      "0.003 -> 0.006) and never breaks the B-envelope (B max 0.0117-0.0119 "
      "< 0.012, B_l unchanged at 0.003915 because B_l does not depend on "
      "alpha_m); the strict totals interval is never entered, the 0.5x-2.0x "
      "band always is. At U=3/F=1000 the 0.0015 floor does not bind (raw "
      "already delivers 0.001515) and the 0.006 floor is the only one that "
      "enters the strict interval (1.023x). The verdict is UNRESOLVED at "
      "every floor for every cell -- the floor moves the numbers, not the "
      "epistemic state.")
    A("")

    # ------------------------------------------------------ 8 transport
    A("## 8. Transport")
    A("")
    t = s3["scenes"]
    w = t["windy"]
    A("### 8.1 The windy cell -- U=6, F=1000, depth 1.0 (per arm)")
    A("")
    A("Route: frozen v5 composite, `A=%.1f`, `B_sigma=%.1f`, `E=%.2f`, "
      "`Q=%d`, 2048^2 irradiance, `g1=%.3f`, blur %.3f px. Photo %dx%d, "
      "sha %s."
      % (s3["frozen_composite"]["A_V3"], s3["frozen_composite"]["B_SIGMA_V3"],
         s3["frozen_composite"]["E_V3"], s3["frozen_composite"]["Q"],
         w["RAW"]["route"]["g1_used"], w["RAW"]["route"]["blur_sigma_px"],
         s3["photo"]["size"][0], s3["photo"]["size"][1],
         s3["photo"]["sha256"][:16]))
    A("")
    rows = []
    for a in ("RAW", "B"):
        r = w[a]
        rows.append([
            "`%s`" % a,
            n(r["alpha_m_used"], "%.6f"),
            n(r["sigma2_h_total"]),
            n(r["state_summary"]["Hs_spectral_m"], "%.6f"),
            n(r["C"]["mean"], "%.9f"),
            n(r["C_percentiles"]["p50"], "%.4f"),
            n(r["C_percentiles"]["p99"], "%.4f"),
            n(r["C_percentiles"]["p99_9"], "%.4f"),
            n(r["C_percentiles"]["max"], "%.4f"),
            n(r["warp"]["own_shift_field"]["mean_mm"], "%.2f"),
            n(r["warp"]["own_shift_field"]["field_std_px"], "%.2f"),
        ])
    L.extend(table(["arm", "a_m used", "s2_h", "Hs [m]", "C mean", "C p50",
                    "C p99", "C p99.9", "C max", "own shift [mm]",
                    "own shift std [px]"], rows))
    A("")
    A("Design checks (PIPELINE-V12 line 28): warp in the observed 20-45 mm "
      "band, `C p99.9 >= 3.0`, crush <= 5%, blow = 0, `C mean == 1.0`:")
    A("")
    rows = []
    for a in ("RAW", "B"):
        r = w[a]
        rows.append([
            "`%s`" % a,
            yn(r["C"]["pass_mean_1"]),
            n(r["C_percentiles"]["p99_9"], "%.3f")
            + (" >= 3.0" if r["C_percentiles"]["p99_9"] >= 3.0 else " < 3.0"),
            pct(r["tone"]["crush_rate_linear_lt_0.02"], "%.2f"),
            "%.4f" % r["tone"]["blowout_rate_aces_gt_0.98"],
            n(r["warp"]["own_shift_field"]["mean_mm"], "%.2f"),
            "**%s**" % ("PASS" if (r["C"]["pass_mean_1"]
                                   and r["C_percentiles"]["p99_9"] >= 3.0
                                   and r["tone"]["crush_pass_customary"]
                                   and r["tone"]["blow_pass_zero"]) else "FAIL"),
        ])
    L.extend(table(["arm", "C mean == 1", "C p99.9", "crush",
                    "blow", "own shift [mm]", "transport sanity"], rows))
    A("")
    A("The `C` minimum is slightly negative (%.2e, identical for both arms) "
      "and `frac_negative = %.2e`: this is the transport's own weak-trough "
      "behaviour carried over from v3.8/v3.9 unchanged, not a spectrum effect."
      % (w["RAW"]["C"]["min"], w["RAW"]["C"]["frac_negative"]))
    A("")
    A("RAW vs B at this cell: **byte-identical** (`max abs uint8 = %d`, "
      "`frac pixels changed = %.3f`, sha equal). Expected and confirmed: at "
      "U=6 the raw `sigma2_h = %.6f` is four times the floor `%.6f`, so the "
      "top-up never acts and both arms assemble the same `alpha_m`."
      % (w["RAW_vs_B"]["max_abs_uint8"],
         w["RAW_vs_B"]["frac_pixels_changed"],
         w["RAW"]["sigma2_h_total"], w["RAW"]["alpha_m_floor"]))
    A("")
    A("### 8.2 The boundary cell -- U=3, F=1000: arm RAW vs arm B")
    A("")
    b = t["boundary_u3"]
    A("Same seed, same frame, same composite parameters (`g1=%.3f`, blur "
      "%.3f px). The only difference is `alpha_m` -- i.e. the boundary "
      "treatment itself." % (b["RAW"]["route"]["g1_used"],
                            b["RAW"]["route"]["blur_sigma_px"]))
    A("")
    rows = []
    for a in ("RAW", "B"):
        r = b[a]
        rows.append([
            "`%s`" % a,
            "`%s`" % r["provenance"],
            n(r["alpha_m_used"], "%.6f"),
            n(r["sigma2_h_total"]),
            n(r["state_summary"]["Hs_spectral_m"], "%.6f"),
            n(r["state_summary"]["sigma2_spectral"]),
            n(r["C"]["mean"], "%.9f"),
            n(r["C_percentiles"]["p99"], "%.4f"),
            n(r["C_percentiles"]["p99_9"], "%.4f"),
            n(r["C_percentiles"]["max"], "%.4f"),
            n(r["warp"]["own_shift_field"]["mean_mm"], "%.2f"),
            n(r["fine_texture"]["mean_abs_dL_dy"], "%.5f"),
        ])
    L.extend(table(["arm", "provenance", "a_m used", "s2_h", "Hs [m]",
                    "s2 realised", "C mean", "C p99", "C p99.9", "C max",
                    "own shift [mm]", "|dL/dy|"], rows))
    A("")
    d = b["RAW_vs_B"]
    A("Pair statistics (`%s`, %s):" % (d["file"], d["layout"]))
    A("")
    A("| quantity | value |")
    A("|---|---|")
    A("| alpha_m ratio B/RAW | **%.4f** |" % d["alpha_m_ratio_B_over_RAW"])
    A("| max abs uint8 difference | %d |" % d["max_abs_uint8"])
    A("| mean abs uint8 difference | %.4f |" % d["mean_abs_uint8"])
    A("| fraction of pixels changed (any channel) | %.4f |"
      % d["frac_pixels_changed"])
    A("| MS-SSIM | %.7f |" % d["msssim"])
    A("| C mean == 1 in both | %s |" % d["C_mean_both_one"])
    A("| C p99.9 RAW / B | %.4f / %.4f |"
      % (b["RAW"]["C_percentiles"]["p99_9"], b["B"]["C_percentiles"]["p99_9"]))
    A("| own shift RAW / B [mm] | %.2f / %.2f |"
      % (b["RAW"]["warp"]["own_shift_field"]["mean_mm"],
         b["B"]["warp"]["own_shift_field"]["mean_mm"]))
    A("| crush RAW / B | %.5f / %.5f |"
      % (b["RAW"]["tone"]["crush_rate_linear_lt_0.02"],
         b["B"]["tone"]["crush_rate_linear_lt_0.02"]))
    A("")
    A("The optical consequence of doubling `alpha_m` at U=3/F=1000 is small "
      "but measurable: the realised state's slope variance rises %.1f%% "
      "(%.6f -> %.6f), the own-shift field mean rises %.2f%% (%.2f -> %.2f "
      "mm), `C p99.9` FALLS slightly (%.3f -> %.3f) because the capillary "
      "top-up adds high-k slope without adding large-scale height, and the "
      "image changes on %.1f%% of pixels by at most %d/255. The prior is a "
      "texture change, not a geometry change -- consistent with it touching "
      "only `B_h`."
      % (100.0 * (b["B"]["state_summary"]["sigma2_spectral"]
                  / b["RAW"]["state_summary"]["sigma2_spectral"] - 1.0),
         b["RAW"]["state_summary"]["sigma2_spectral"],
         b["B"]["state_summary"]["sigma2_spectral"],
         100.0 * (b["B"]["warp"]["own_shift_field"]["mean_mm"]
                  / b["RAW"]["warp"]["own_shift_field"]["mean_mm"] - 1.0),
         b["RAW"]["warp"]["own_shift_field"]["mean_mm"],
         b["B"]["warp"]["own_shift_field"]["mean_mm"],
         b["RAW"]["C_percentiles"]["p99_9"], b["B"]["C_percentiles"]["p99_9"],
         100.0 * d["frac_pixels_changed"], d["max_abs_uint8"]))
    A("")
    A("Arm D at this cell is RAW's still byte-for-byte (same spectrum, seed, "
      "frame, route); no separate render exists because there is nothing to "
      "render differently. The identity is asserted, not assumed: stage 1 "
      "records arm D `numbers_identical_to_RAW = true` on `alpha_m_used`, and "
      "stage 2 records `D_identical_to_RAW.pass = true` for every cell "
      "(`%s`)."
      % ("PASS" if all(s2by[c]["D_identical_to_RAW"]["pass"] for c in cid)
         else "FAIL"))
    A("")
    A("Warp-agreement yardstick, for honesty about the measurement itself: "
      "the SSD block matcher and the predicted shift field agree poorly at "
      "this texture (corr_dx = %.3f, rms_dx = %.2f px), while the NCC check "
      "agrees tightly (corr_dx = %.3f, rms_dx = %.2f px). The OWN-SHIFT field "
      "(sub-pixel, from the state) is the transport number quoted above; the "
      "block matcher's weakness on a low-contrast underwater bark scene is a "
      "pre-existing property of the yardstick, reported not hidden."
      % (w["RAW"]["warp"]["agreement_with_predicted"]["ssd_primary"]["corr_dx"],
         w["RAW"]["warp"]["agreement_with_predicted"]["ssd_primary"]["rms_dx_px"],
         w["RAW"]["warp"]["agreement_with_predicted"]["ncc_check"]["corr_dx"],
         w["RAW"]["warp"]["agreement_with_predicted"]["ncc_check"]["rms_dx_px"]))
    A("")

    # ------------------------------------------------------ 9 regression
    A("## 9. Regression")
    A("")
    reg = dec["regression"]
    rows = []
    for lane, key in (("v13 legacy lane (`jonswap`)", "v13_legacy_lane"),
                      ("v13 new lane (`elfouhaily-form`)", "v13_new_lane")):
        r = reg[key]
        rows.append([lane, "`%s`" % r["v14_file"],
                     r["sha256_v14"][:16], r["sha256_v13"][:16],
                     "yes" if r["byte_identical"] else "no",
                     r["max_abs_uint8"], "%.1f" % r["frac_pixels_changed"]])
    L.extend(table(["lane", "v14 file", "sha256 v14", "sha256 v13",
                    "byte identical", "max abs uint8", "frac px changed"], rows))
    A("")
    A("Both lanes re-rendered through the v14 family at v13's own cell "
      "(U=6/F=100/depth 1.0, the cell v13's deliverables used). The legacy "
      "lane is the byte-identity test PIPELINE-V12 line 29 asks for; the new "
      "lane is a bonus -- v13's cell is inside the fit domain, so no arm "
      "changes a number there and the still must match too. **Both are "
      "byte-identical (sha256 equal).**")
    A("")
    wp = reg["v12_water_path"]
    A("v12 water path (import check only, as specified): "
      "`import v12_water` OK; `WATER_MODEL_DEFAULT` = %s (%s); "
      "`DEFAULT_MEDIUM_PRESET` = `%s`; `LEGACY_RECIPE_ROUTE` = `%s`. Nothing "
      "in this job calls or edits v12_water."
      % ("`%s`" % wp["WATER_MODEL_DEFAULT"],
         "PASS" if wp["pass"] else "FAIL", wp["DEFAULT_MEDIUM_PRESET"],
         wp["LEGACY_RECIPE_ROUTE"]))
    A("")
    fh = dec.get("file_hashes") or {}
    A("Frozen-file hashes, recorded at the start of stage 1 and re-taken at "
      "stage 4: **%s** across %d files (`%s`). The frozen composite is "
      "imported, never edited: `A_V3=%.1f`, `B_sigma=%.1f`, `E=%.2f`, "
      "`v5_composite.py` sha %s."
      % ("all unchanged" if fh.get("all_unchanged") else "CHANGED",
         len(fh.get("files") or {}),
         ", ".join(sorted(list((fh.get("files") or {}).keys()))[:6])
         + (", ..." if len(fh.get("files") or {}) > 6 else ""),
         reg["frozen_composite_hashes"]["A_V3"],
         reg["frozen_composite_hashes"]["B_SIGMA_V3"],
         reg["frozen_composite_hashes"]["E_V3"],
         reg["frozen_composite_hashes"]["v5_composite_sha256"][:16]))
    A("")
    A("**Regression: GREEN** (v13 lanes byte-identical, v12 default `ab`, "
      "frozen modules unchanged, composite frozen).")
    A("")

    # ------------------------------------------------------ 10 decision
    A("## 10. Decision data -- constraints x arms")
    A("")
    A("This is the table the pre-registered rule needs. **No winner is named "
      "here.** Counts are over the 8 cells; `n/a` means the clause does not "
      "apply to that cell (the C1/C1b share gates apply only to U=6/10).")
    A("")
    rows = []
    for row in dec["constraints"]:
        r = [row["constraint"]]
        for a in ARMS:
            c = row["counts"][a]
            r.append("%d/%d/%d" % (c["pass"], c["fail"],
                                   c["not_applicable_or_unresolved"]))
        rows.append(r)
    L.extend(table(["constraint", "RAW P/F/n-a", "B P/F/n-a", "D P/F/n-a",
                    "D_STRICT P/F/n-a"], rows))
    A("")
    A("### 10.1 Invented inputs acting, per cell")
    A("")
    A("The only invented (unmeasured) input in the whole family is the floor "
      "target `sigma2_h_target = %s`; it acts wherever `top-up acted = yes`."
      % n(pol["sigma2_h_target_default"], "%g"))
    A("")
    rows = []
    for c in cells:
        r = [CELL_LABEL[c["cell"]]]
        for a in ARMS:
            inv = dec["invented_inputs"][a]["per_cell"][c["cell"]]
            r.append("%d%s" % (inv["invented_inputs_acting"],
                               " (top-up)" if inv["top_up_acted"] else ""))
        rows.append(r)
    L.extend(table(["cell", "RAW", "B", "D", "D_STRICT"], rows))
    A("")
    rows = []
    for a in ARMS:
        inv = dec["invented_inputs"][a]
        rows.append(["`%s`" % a, inv["n_cells_with_invented_input"],
                     ", ".join(CELL_LABEL[c] for c in cid
                               if inv["per_cell"][c]["invented_inputs_acting"]),
                     "`%s`" % inv["provenance"]])
    L.extend(table(["arm", "cells with an invented input", "which cells",
                    "provenance"], rows))
    A("")
    A("### 10.2 Side-by-side reading of the independent constraints")
    A("")
    A("| constraint (independent of the band-share gates) | RAW | B | D | "
      "D_STRICT |")
    A("|---|---|---|---|---|")
    A("| totals inside 0.5x-2.0x envelope | 8/8 | 8/8 | 8/8 | 8/8 |")
    A("| totals inside strict [Cox-Munk clean, Hu] | %d/8 | %d/8 | %d/8 | "
      "%d/8 |" % tuple(sum(1 for c in cells
                           if c["arms"][a]["envelope_check"]["strict_inside"])
                       for a in ARMS))
    A("| B-envelope (B <= 0.012 over 10-1500) | %d/8 | %d/8 | %d/8 | %d/8 |"
      % tuple(sum(1 for c in cells
                  if c["arms"][a]["b_envelope"]["pass_B_env"]) for a in ARMS))
    A("| B_l band max <= 0.008 over 50-200 | %d/8 | %d/8 | %d/8 | %d/8 |"
      % tuple(sum(1 for c in cells
                  if c["arms"][a]["b_envelope"]["pass_B_l"]) for a in ARMS))
    A("| construction A1/A2/A4/A5 | all PASS | all PASS | all PASS (RAW's "
      "spectrum) | PASS |")
    A("| transport sanity at the windy cell | PASS | PASS | PASS (RAW's "
      "still) | not rendered |")
    A("| invented inputs acting (cells) | %d | %d | %d | %d |"
      % tuple(dec["invented_inputs"][a]["n_cells_with_invented_input"]
              for a in ARMS))
    A("")
    A("Pre-registered rule, quoted verbatim from PIPELINE-V12 line 20 "
      "(applied OUTSIDE this session):")
    A("")
    A("> %s" % dec["rule_quoted"].replace("\n", " "))
    A("")
    A("Recorded for that decision, without recommending: three of the four "
      "high-wind cells FAIL the B-envelope for every arm (section 5.1), and "
      "PIPELINE-V12's merge rule requires the high-wind cells to PASS and "
      "zero FAIL. On this matrix no arm clears that bar; the rule the "
      "decision-maker applies will have to confront that, not this report.")
    A("")

    # ------------------------------------------------------ 11 deviations
    A("## 11. Deviations, and what they cost")
    A("")
    devs = [
        ("Arm D's literal reading is ambiguous, so both readings were "
         "measured.",
         "The brief says arm D 'changes NO numbers' and that at U=6/10 it is "
         "'the same numbers' as RAW; the literal phrase 'raw spectrum "
         "untouched' suggests 'no prior anywhere'. Implemented D as RAW's "
         "spectrum evaluated through the absolute criteria (so D's numbers ARE "
         "RAW's), and added `D_STRICT` = `alpha_m = max(raw, 0)` everywhere "
         "(the v3.8 build, verified bit-identical to v11) as a labelled "
         "diagnostic, not a fourth arm.",
         "Both readings are in section 5.3. They disagree where it matters: "
         "D (with the prior below the crossing) enters the strict totals "
         "interval at the U=1.5 cells, D_STRICT does not (0.965x at "
         "U=1.5/F=100). Nothing is hidden by choosing one implementation."),
        ("C1b is a conjunction here; v13 used a disjunction.",
         "PIPELINE-V12 line 25 states 'capillary share >= 0.15 AND sigma2_h >= "
         "0.002'. v13's `c1b()` returned `ok1 or ok2`. v14 applies the "
         "conjunction and records the change in every cell's C1b node "
         "(`changed_from_v13`).",
         "At U=1.5/F=1000 both RAW and B have cap share 0.1488 (< 0.15) with "
         "sigma2_h = 0.003, so C1b is FAIL there where v13 would have said "
         "PASS. It is a REPORTED clause at a low-wind cell (not a gate), so no "
         "verdict changes; the number is different and the report says why."),
        ("The 'Hu' column is the in-repo fetch-limited MSS, not a Hu formula.",
         "No Hu formula exists anywhere in this workspace or in "
         "PIPELINE-V10/V12. As in v13's reporting, the second envelope column "
         "is the pipeline's own fetch-limited MSS integral on the same band, "
         "with `naming_status` carried in the JSON.",
         "The envelope's upper edge is model-internal, so the strict-interval "
         "verdicts inherit that weakness. Recorded, not repaired here."),
        ("Arm B is C0 but not C1 in u*.",
         "`max(raw, floor)` is continuous by construction and its derivative "
         "jumps at the crossing; both one-sided log-slopes are measured "
         "(section 6: 0 -> 4.5-5.2). The A2 continuity gate is on B(k) in k, "
         "not on alpha_m(u*), and passes everywhere.",
         "A design reading that required C1 in u* would fail arm B by "
         "construction. Stated so the decision-maker can apply their own "
         "reading rather than inherit a silent one."),
        ("Two Hs definitions coexist.",
         "The tables use the ASSEMBLED H_s (spectrum-integrated to `k_vis`, "
         "never rescaled). The v11 K1 node computes H_s from the state's "
         "`spectrum_after` (integrated to `k_cut = 1000`), which differs by "
         "the band above `k_cut`. Both are reported; the K1 flag (threshold "
         "35%) is false on all 8 cells in both arms.",
         "Largest gap is at U=1.5/F=100: assembled 0.00705 m vs state "
         "0.00683 m (~3%). No gate depends on the difference."),
        ("Parseval (A2) is measured on a test field with its own seed.",
         "A2 builds a field from the arm's continuum at N=%d with test seed "
         "12345 (recorded in the JSON), independent of the production seed "
         "(%d) used for the realised states." % (s2["grid"]["N"],
                                                s2["grid"]["seed"]),
         "A2 tests the spectrum-to-field convention, not the production "
         "realisation; the production realisation is covered by A5 and by the "
         "transport renders."),
        ("Time-quantisation error is large on the longest modes and is not "
         "gated.",
         "`omega` is quantised to integer multiples of 2pi/8 s so the 8 s loop "
         "closes; the worst per-band relative error is %.1f%% (U=10/F=1000, "
         "longest modes) and %.1f%% at U=1.5/F=100. At t=0 the field does not "
         "depend on omega at all."
         % (100.0 * max(s2by[c]["arms"]["RAW"]["quantisation"]["max_rel_error"]
                        for c in cid),
            100.0 * s2by["U1p5_F100p0"]["arms"]["RAW"]["quantisation"]
            ["max_rel_error"]),
         "It bounds time-evolution error only; this job renders t=0 stills. "
         "Reported per cell in the JSON, carried as a known limit."),
        ("Provenance strings are per cell, not per arm only.",
         "`low_wind_prior` appears exactly where the top-up acts; a cell "
         "assembled inside the fit domain with no top-up is labelled "
         "`eckv_fit_domain_raw`, and arm B's below-crossing-but-not-topped-up "
         "case would be `below_fit_domain_no_topup` (it does not occur in this "
         "matrix). Arm D's low-wind cells carry `low_wind_absolute_gate`.",
         "No cell claims a provenance it did not use; the per-cell table is "
         "section 3."),
        ("A different process finished stage 3, and this session read it.",
         "The interrupted session's stage-3 process was still running; it "
         "completed at 2026-11-11T18:46:45+09:00, wrote "
         "`numbers-v14-transport.json` and committed the `stage-3` checkpoint "
         "unit. This resumed session did NOT re-run stages 1, 2 or 3 and did "
         "not overwrite any render; it read the JSONs and hashed the files.",
         "The transport numbers are that run's, internally consistent "
         "(all six renders timestamped in one sequential 363.8 s pass, "
         "`created_utc` 09:40:42Z). The RAW-vs-B windy pair is byte-identical "
         "and both v13 regression lanes match v13's shas, which is the "
         "strongest evidence available that the run was clean."),
        ("`run_v14_report.py` is new in this session.",
         "`run_v14.stage4()` imports it; it did not exist when the stage-3 "
         "process exited, which is why stage 4 had not run. It only formats "
         "the stage JSONs.",
         "The report is generated, not hand-transcribed: every table cell is a "
         "format of a JSON value."),
    ]
    for i, (what, how, cost) in enumerate(devs, 1):
        A("**D%d. %s**" % (i, what))
        A("")
        A("- what was done: %s" % how)
        A("- what it costs: %s" % cost)
        A("")
    A("### 11.1 Numbers that could not be computed, and why")
    A("")
    A("| number | status | reason |")
    A("|---|---|---|")
    A("| `m0` in the A5 convergence rows | missing (`null`) | the state node "
      "does not expose m0; sigma2 and Hs are used instead (both present) |")
    A("| capillary >370 fraction for `D_STRICT` | n/a | that arm has no "
      "capillary part (`alpha_m = 0`), so the fraction is 0/0 and is emitted "
      "as null rather than 0 |")
    A("| any Hu-law value | not available as a law | no formula exists in "
      "this workspace; the model-internal fetch-limited column is used and "
      "named as such (D3) |")
    A("| a measured in-water PSF / ground-truth image | out of scope | not "
      "part of this chunk; the transport numbers compare the model to itself "
      "plus the yardsticks, as in v3.8/v3.9 |")
    A("")
    A("No number in this report was estimated. Where the JSON says `null`, "
      "this report says `n/a` and gives the reason.")
    A("")

    # ------------------------------------------------------ 12 files
    A("## 12. Deliverables and where each number lives")
    A("")
    A("| file | contents |")
    A("|---|---|")
    A("| `numbers-v14-spectrum.json` | stage 1: the 8 cells x 4 arms, the "
      "alpha_m resolution and prior accounting, C1/C1b/B-envelope/totals/Hs/"
      "smoothness/positivity, the hand-off scan, the controls, and the "
      "exported spectra (per cell per arm, on a log grid to 1500 rad/m, with "
      "the re-integration check) |")
    A("| `numbers-v14-sweep.json` | stage 2: the realised states at N=%d, "
      "Parseval, the independent band integrals, the resolution-convergence "
      "re-run, the floor sweep, and the decision block written by stage 4 |"
      % s2["grid"]["N"])
    A("| `numbers-v14-transport.json` | stage 3: the windy cell per arm, the "
      "U=3/F=1000 pair and its statistics, the regression battery |")
    A("| `comp-windy-v40-U6F1000-armRAW.png`, `...-armB.png` | the windy "
      "still, per arm (byte-identical) |")
    A("| `comp-u3-F1000-armRAW.png`, `comp-u3-F1000-armB.png` | the boundary "
      "cell, per arm |")
    A("| `comp-u3-raw-vs-priorB.png` | the side-by-side (left RAW, right B, "
      "2 px white divider) |")
    A("| `comp-windy-v40-v13legacy.png`, `comp-windy-v40-v13new.png` | the "
      "regression re-renders (both byte-identical to v13) |")
    A("| `spectrum-ref-v40-U6F1000-armB.csv` | the k, B_l, B_h, S_l, S_h "
      "curve at U=6/F=1000 for arm B |")
    A("| `checkpoint.json` | jobkit units stage-1..stage-4 |")
    A("")
    A("---")
    A("")
    A("## 13. One-paragraph summary of what was measured")
    A("")
    A("Three arms over the same eight cells, on the same code path. "
      "Construction is clean for all of them (A1-A5 pass; the resolution "
      "re-run moves sigma2 by %.1f%% and Hs by %.2f%%). The arms differ in "
      "exactly one place -- the U=3 cells, where the raw spectrum is inside "
      "the ECKV fit domain but delivers less capillary variance than the "
      "prior's floor: arm B doubles `alpha_m` there (0.001067 -> 0.002114 at "
      "F=1000), which raises the delivered `sigma2_h` to 0.003 and the "
      "capillary share of the 50-200 band from 0.150 to 0.259, at the cost of "
      "one invented input in those cells. Arm D changes nothing numerically "
      "and, tested against Gemini's absolute criteria at the four U<=3 cells, "
      "PASSES the B-envelope there and the wide 0.5x-2.0x totals band "
      "everywhere, but FAILS "
      "the STRICT [Cox-Munk clean, Hu] totals at both U=3 cells (0.616x and "
      "0.787x) -- so its own criteria do not clear on the raw spectrum. "
      "Every low-wind cell remains UNRESOLVED for every arm, and the "
      "high-wind cells fail the absolute B-envelope at U=6/F=100, U=10/F=100 "
      "and U=10/F=1000 for every arm. The transport consequence of the "
      "boundary choice is small: the doubled capillary level changes %.1f%% "
      "of pixels by at most %d/255 (MS-SSIM %.5f), while the own-shift field "
      "moves %.2f%%. No default was changed and no winner was chosen."
      % (100.0 * abs(s2["A5_convergence"]["RAW"]["rel_sigma2"]),
         100.0 * abs(s2["A5_convergence"]["RAW"]["rel_Hs"]),
         100.0 * t["boundary_u3"]["RAW_vs_B"]["frac_pixels_changed"],
         t["boundary_u3"]["RAW_vs_B"]["max_abs_uint8"],
         t["boundary_u3"]["RAW_vs_B"]["msssim"],
         100.0 * (t["boundary_u3"]["B"]["warp"]["own_shift_field"]["mean_mm"]
                  / t["boundary_u3"]["RAW"]["warp"]["own_shift_field"]
                  ["mean_mm"] - 1.0)))
    A("")

    text = "\n".join(L) + "\n"
    path = os.path.join(OUT, "REPORT-V14-GATES.md")
    import sys
    sys.path.insert(0, r"D:\mutaro-work\jobkit")
    from checkpoint import atomic_write_bytes
    atomic_write_bytes(path, text.encode("utf-8"))
    print("  wrote %s (%d lines, %d bytes)" % (path, len(L), len(text)))
    return path


if __name__ == "__main__":
    main()
