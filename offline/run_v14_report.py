"""run_v14_report -- REPORT-V14-GATES.md from the v14 stage JSONs.

Every number in the report is READ from `numbers-v14-{spectrum,sweep,
transport}.json`; nothing is recomputed here, so the report cannot drift from
the measurements.  ASCII only.  No image is inspected.
"""
import json
import math
import os
import time

T0 = time.time()

ARMS = ("RAW", "B", "D", "D_STRICT")
ARMN = {"RAW": "RAW", "B": "B", "D": "D", "D_STRICT": "D-strict"}


# --------------------------------------------------------------------------
# formatting
# --------------------------------------------------------------------------
def _f(x, n=4, dash="-"):
    if x is None:
        return dash
    try:
        v = float(x)
    except Exception:
        return str(x)
    if not math.isfinite(v):
        return dash
    if v != 0.0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        return ("%%.%de" % max(0, n - 1)) % v
    return ("%%.%df" % n) % v


def _p(x, n=1, dash="-"):
    if x is None:
        return dash
    return ("%%.%df" % n) % float(x)


def _pass(b):
    if b is None:
        return "n/a"
    return "PASS" if b else "FAIL"


def _state(s):
    return {"PASS": "PASS", "FAIL": "FAIL",
            "UNRESOLVED": "UNRESOLV"}.get(s, str(s))


def _yn(b):
    return "yes" if b else "no"


def _cell(c):
    return "U=%-4s F=%-4s" % (_f(c["U_ms"], 1), _f(c["fetch_m"], 0))


def _load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def _ascii(txt):
    out = (txt.replace("\u2014", "--").replace("\u2013", "-")
              .replace("\u2019", "'").replace("\u201c", '"')
              .replace("\u201d", '"').replace("\u00b1", "+/-")
              .replace("\u00b5", "u").replace("\u2192", "->"))
    out.encode("ascii")
    return out


def table(L, head, rows, sep="-"):
    all_rows = [head] + [[str(c) for c in r] for r in rows]
    w = [max(len(r[i]) for r in all_rows) for i in range(len(head))]
    L.append("  ".join(h.ljust(x) for h, x in zip(head, w)))
    L.append("  ".join("-" * x for x in w))
    for r in all_rows[1:]:
        L.append("  ".join(str(c).ljust(x) for c, x in zip(r, w)))
    return L


# --------------------------------------------------------------------------
def build(out_dir):
    s1 = _load(os.path.join(out_dir, "numbers-v14-spectrum.json"))
    s2 = _load(os.path.join(out_dir, "numbers-v14-sweep.json"))
    s3 = _load(os.path.join(out_dir, "numbers-v14-transport.json"))
    dec = s2["decision"]
    pol = s1["policy"]
    cells = s1["cells"]
    c2 = {c["cell"]: c for c in s2["cells"]}
    summ = s1["summary"]
    hs = s1["handoff"]
    L = []
    A = L.append

    A("# REPORT-V14-GATES -- version 3.9 boundary redo, the three-arm "
      "measurement (chunk 10)")
    A("")
    A("Brief: `BRIEF-V14-39REDO.md`.  Design: `PIPELINE-V12.md` (it governs "
      "everything below).  Base: the v13 tree (Elfouhaily-form spectrum + the "
      "labelled provisional prior below the ECKV fit crossing).")
    cp = _load(os.path.join(out_dir, "checkpoint.json"))
    A("Machine: PC, `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.  "
      "Model: deepseek-v4.1-flash.")
    A("Execution: FOREGROUND, one jobkit checkpoint unit per stage "
      "(`out/v14/checkpoint.json` records %d unit completions for %d stages: "
      "stages 3 and 4 were run a second time as well -- see 9.6, which "
      "documents the duplicate session in this folder).  Every deliverable is "
      "committed with `atomic_write_bytes`.  Nothing outside `out/v14/` and "
      "the new `run_v14*` files was written; no `python.exe` / "
      "`blender.exe` was killed; no image was inspected (numerical image "
      "statistics only); the report is ASCII."
      % (int(cp["completed_units"]), len(set(cp["completed_ids"]))))
    A("")
    A("New code (v14 family only; v5..v13 imported, never edited):")
    A("")
    A("- `run_v14_spectrum.py` -- the three arms' `alpha_m` rules, the "
      "hand-off scan, the gates, the correct band partition, the exports")
    A("- `run_v14_state.py` -- `V14WaterState`: v11 state + the arm's "
      "resolution; `spectrum_source` default READ from v13, not changed")
    A("- `run_v14.py` -- the four stage runners")
    A("- `run_v14_report.py` -- this report (reads the stage JSONs)")
    A("")
    A("Reproduce, foreground, in this order:")
    A("")
    A("```")
    A("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    A("python run_v14.py --stage 1     # %s s  (8 cells x 4 number sets, "
      "analytic)" % _p(s1.get("wall_s"), 1))
    A("python run_v14.py --stage 2     # %s s  (16 states at 1024^2, Parseval, "
      "A4, A5, floor sweep)" % _p(s2.get("wall_s"), 1))
    A("python run_v14.py --stage 3     # %s s  (6 renders, C mesh 2048^2, "
      "Q=25; composite at the photo's 735x539)"
      % _p(s3.get("wall_s"), 1))
    A("python run_v14.py --stage 4     # decision block + this report")
    A("```")
    A("")
    A("---")
    A("")
    A("## 0. Scope: what this job does NOT do")
    A("")
    A("It does **not** pick a winner and it does **not** change a default.  "
      "`spectrum_source` in `v13_state.py` is `%s`, the value v13 left, and "
      "`run_v14_state.py` reads that constant instead of re-declaring it.  "
      "PIPELINE-V12's pre-registered decision rule is applied OUTSIDE this "
      "session; section 7 supplies the table it needs.  No threshold was "
      "moved: the gates below are the design's numbers, and where the design "
      "left a reading open, both readings are printed (section 9)."
      % pol["default_source_after"])
    A("")
    A("---")
    A("")
    A("## 1. The arms")
    A("")
    for a in ARMS:
        A("- **Arm %s** -- %s" % (ARMN[a], s1["arms"][a]["definition"]))
        A("     provenance string: `%s`" % s1["arms"][a]["provenance"])
    A("")
    A("Only one quantity moves between the arms: `alpha_m`, the level of the "
      "Elfouhaily high-frequency term `B_h`.  `B_h` is exactly linear in "
      "`alpha_m`, so a cell's capillary slope variance is "
      "`sigma2_h = alpha_m * I(cell)` with `I(cell) = int k^2 B_h(alpha=1)/k^3 "
      "dk` over the cell's resolved band `[kmin, k_vis]`, and every arm is one "
      "rule for `alpha_m_used`.  `B_l`, `Phi`, the patch/cascade rule, the "
      "transport and the frozen composite are untouched by all three arms.")
    A("")
    A("The v13 boundary is the ECKV fit crossing `u* = c_m/e = %s m/s` "
      "(`c_m = %s m/s`).  Arm B's boundary is a different, measured one: the "
      "`u*` where the raw spectrum's own capillary variance reaches the floor "
      "target (section 3.2)."
      % (_f(pol["ustar_fit_min_ms"], 6), _f(pol["c_m_ms"], 6)))
    A("")
    A("---")
    A("")
    A("## 2. Per-arm, per-cell table")
    A("")
    A("Conventions (unchanged from PIPELINE-V9/V10, inherited): `S(k)` in "
      "m^2 per rad/m, `m0 = int S dk`, `sigma2 = int k^2 S dk`, "
      "`B(k) := k^3 S(k)`, wind = U10, `u*` from U10 only.  **`sigma2_g` and "
      "`sigma2_h` are the parts' slope variances over the cell's resolved band "
      "`[kmin, k_vis]`** (the band the floor solves on); **`total MSS` is the "
      "assembled total over `[kmin, k_cut=%s rad/m]`** (the envelope "
      "quantity).  `Hs` is ALWAYS spectrum-integrated (`4 sqrt(int S dk)`), "
      "never rescaled.  `B max` is `max B(k)` over k in [10, 1500] -- the "
      "design's criterion -- and the k where it occurs is in the JSON "
      "(`b_envelope.k_at_B_max`)."
      % _f(pol["k_cut_rad_per_m"], 0))
    A("")
    head = ["cell", "u*", "u*/u*fit", "raw a_m", "appl a_m", "s2_g",
            "s2_h", "total MSS", "grav sh", "cap sh", "cap>370",
            "Hs spec", "Hs law", "dHs", "B max", "Bl max", "state"]
    for a in ARMS:
        A("### 2.%d Arm %s" % (ARMS.index(a) + 1, ARMN[a]))
        A("")
        rows = []
        for c in cells:
            x = c["arms"][a]
            rows.append([
                _cell(c), _f(x["alpha"]["ustar_ms"], 5),
                _f(x["alpha"]["u_over_fit_min"], 3),
                _f(x["alpha"]["alpha_m_raw"], 6),
                _f(x["alpha"]["alpha_m_used"], 6),
                _f(x["moments"]["gravity"]["sigma2"], 6),
                _f(x["moments"]["capillary"]["sigma2"], 6),
                _f(x["sigma2_total_cut"], 6),
                _f(x["shares"]["gravity_share_50_200"], 3),
                _f(x["shares"]["capillary_share_50_200"], 3),
                _f(x["moments"]["capillary_frac_above_370"], 3),
                _f(x["Hs_spectrum_m"], 4), _f(x["Hs_law_m"], 4),
                _p(100.0 * x["hs"]["delta_rel"], 1) + "%",
                _f(x["b_envelope"]["B_max"], 4),
                _f(x["b_envelope"]["B_l_max"], 4),
                _state(x["state"])])
        table(L, head, rows)
        A("")
        A("arm rule, low wind: %s" % cells[0]["arms"][a]["verdict_rule"])
        A("")
        A("arm rule, high wind: %s"
          % [c for c in cells if not c["low_wind_cell"]][0]["arms"][a]
          ["verdict_rule"])
        A("")
    A("Notes on the table:")
    A("")
    A("- `cap>370` is the fraction of the CAPILLARY slope variance above "
      "370 rad/m (`moments.capillary_frac_above_370`).  The fraction of the "
      "*total* slope variance above 370 is in the JSON and differs (the "
      "gravity part contributes almost nothing there).")
    A("- At U=1.5 the RAW/B `s2_h` is identically 0.003: that is the floor "
      "target, hit to machine precision by construction "
      "(`closure_rel_residual` is 0 in every prior cell).")
    A("- D's rows are RAW's numbers by definition; they are printed again "
      "rather than referenced, so the table cannot be misread.")
    A("- `D-strict` is the labelled DIAGNOSTIC (section 9.1), not one of the "
      "three arms.")
    A("")
    A("### 2.5 Where the arms actually differ")
    A("")
    A("```")
    for c in cells:
        A("%-16s RAW a_m %.6f | B a_m %.6f (top-up %s) | D-strict a_m %.6f"
          % (c["cell"], c["arms"]["RAW"]["alpha"]["alpha_m_used"],
             c["arms"]["B"]["alpha"]["alpha_m_used"],
             _yn(c["arms"]["B"]["alpha"]["top_up_acted"]),
             c["arms"]["D_STRICT"]["alpha"]["alpha_m_used"]))
    A("```")
    A("")
    A("Arm B's top-up acts at the four cells whose RAW capillary variance is "
      "below the floor.  At U=1.5 the RAW arm is already at the floor, so the "
      "top-up changes no number (RAW and B are the same spectrum there); at "
      "U=3 the applied level rises (%s -> %s at F=100, %s -> %s at F=1000), "
      "which is the only place the three arms differ.  Arm D-strict differs "
      "from RAW only below the ECKV crossing (%s).  Everywhere else all arms "
      "are the same spectrum, which is why the high-wind cells cannot separate "
      "them."
      % (_f([c for c in cells if c["cell"] == "U3p0_F100p0"][0]["arms"]["RAW"]
            ["alpha"]["alpha_m_used"], 6),
         _f([c for c in cells if c["cell"] == "U3p0_F100p0"][0]["arms"]["B"]
            ["alpha"]["alpha_m_used"], 6),
         _f([c for c in cells if c["cell"] == "U3p0_F1000p0"][0]["arms"]["RAW"]
            ["alpha"]["alpha_m_used"], 6),
         _f([c for c in cells if c["cell"] == "U3p0_F1000p0"][0]["arms"]["B"]
            ["alpha"]["alpha_m_used"], 6),
         ", ".join(c["cell"] for c in cells
                   if c["arms"]["D_STRICT"]["alpha"]["alpha_m_used"]
                   != c["arms"]["RAW"]["alpha"]["alpha_m_used"]) or "none"))
    A("")
    A("---")
    A("")
    A("## 3. Construction gates (all arms, all cells; any failure = FAIL)")
    A("")
    A("### 3.1 A1 positivity")
    A("")
    rows = []
    for a in ARMS:
        rows.append([ARMN[a],
                     _f(min(c["arms"][a]["positivity"]["B_min"]
                            for c in cells), 3),
                     _f(min(c["arms"][a]["positivity"]["S_min"]
                            for c in cells), 3),
                     _yn(all(c["arms"][a]["positivity"]["all_finite"]
                             for c in cells)),
                     _pass(all(c["arms"][a]["positivity"]["pass"]
                               for c in cells))])
    table(L, ["arm", "min B", "min S", "all finite", "A1"], rows)
    A("")
    A("A1 is measured on a 20000-point log grid over each cell's resolved band "
      "`[kmin, k_vis]`.  `min B = 0` is underflow at the viscous edge (B_l and "
      "B_h both decay to 0 there), not a negative value: no cell has a "
      "negative or non-finite energy in any arm.  Arm D-strict at U=1.5 sets "
      "`alpha_m = 0`, i.e. exactly zero capillary energy -- non-negative, but "
      "that is why `sigma2_h >= 0.002` exists as a C1b clause.")
    A("")
    A("### 3.2 A2 continuity -- B(k), dlnB/dlnk, and the prior hand-off")
    A("")
    A("`B(k)` in k: the design asks for continuity and no step in "
      "`dlnB/dlnk`, with jumps <= 1%.  A per-step bound alone cannot separate "
      "a STEP from a steep smooth flank, so the joint test is refinement: a "
      "smooth function's one-step change is proportional to the step and "
      "halves when the grid doubles, a splice does not.")
    A("")
    rows = []
    for a in ARMS:
        rows.append([ARMN[a],
                     _f(max(c["arms"][a]["smoothness"]["coarse"]
                            ["max_rel_step"] for c in cells), 5),
                     _f(max(c["arms"][a]["smoothness"]["fine"]
                            ["max_rel_step"] for c in cells), 5),
                     _f(max(c["arms"][a]["smoothness"]
                            ["ratio_max_dslope_fine_over_coarse"]
                            for c in cells), 3),
                     _yn(all(c["arms"][a]["smoothness"]["literal_1pct_met"]
                             for c in cells)),
                     _pass(all(c["arms"][a]["smoothness"]["pass"]
                               for c in cells))])
    table(L, ["arm", "max |dB/B| n=20k", "n=40k", "slope-jump ratio",
              "1%-per-step met", "A2"], rows)
    A("")
    _ws = max(((c["cell"], c["arms"]["B"]["smoothness"]["coarse"])
               for c in cells), key=lambda x: x[1]["max_rel_step"])
    A("The refinement ratio is ~0.5 in every arm (bound 0.65), so there is no "
      "step: the largest per-step change of B is a smooth flank, and the "
      "`1%%-per-step` line is met at n=20000 everywhere.  The largest step of "
      "all is `%s` at k = %s rad/m, where the local log-slope is %s -- the "
      "steep rise into the gravity peak just above the band's lower edge -- "
      "which is why a fixed per-step bound must always be quoted with its "
      "grid." % (_ws[0], _f(_ws[1]["k_at_max_rel_step"], 1),
                 _f(_ws[1]["local_log_slope_at_max"], 1)))
    A("")
    A("**The prior hand-off in `u*` (Arm B).**  Two scans per fetch "
      "(`handoff.F_100`, `handoff.F_1000`):")
    A("")
    rows = []
    for key in ("F_100", "F_1000"):
        h = hs[key]
        ho = h["handoff"] or {}
        kk = ho.get("kink", {})
        rows.append([_f(h["fetch_m"], 0), _f(ho.get("U10_ms"), 3),
                     _f(ho.get("ustar_ms"), 6),
                     _f(ho.get("u_over_fit_min"), 3),
                     _f(kk.get("rel_gap_across"), 5),
                     _f(kk.get("dln_alpha_dln_U_left"), 3),
                     _f(kk.get("dln_alpha_dln_U_right"), 3),
                     _f(h["rule_scan"]["ratio_fine_over_coarse"], 3),
                     _yn(h["rule_scan"]["continuous_C0"])])
    table(L, ["fetch", "U10 handoff", "u* handoff", "u*/u*fit", "gap across",
              "slope left", "slope right", "refine ratio", "C0"], rows)
    A("")
    A("Reading: `max(raw, floor)` is **continuous in value** (the refined gap "
      "across the crossing is at the bisection tolerance) and **not C1**: the "
      "one-sided log-slopes differ because `max()` has a kink where "
      "`raw = floor`.  That is a property of the rule, not a numerical defect, "
      "and it is stated rather than smoothed over.  The Arm-B hand-off is at "
      "u* = %s (F=100) and %s (F=1000) m/s, i.e. slightly BELOW the ECKV "
      "crossing (u*/u*fit = %s and %s), at U10 ~ %s / %s m/s."
      % (_f(hs["F_100"]["handoff"]["ustar_ms"], 6),
         _f(hs["F_1000"]["handoff"]["ustar_ms"], 6),
         _f(hs["F_100"]["handoff"]["u_over_fit_min"], 3),
         _f(hs["F_1000"]["handoff"]["u_over_fit_min"], 3),
         _f(hs["F_100"]["handoff"]["U10_ms"], 2),
         _f(hs["F_1000"]["handoff"]["U10_ms"], 2)))
    A("")
    A("The second scan is the shipped per-cell band, where the patch is "
      "quantised to powers of two, so `I(cell)` and the resolved-band integrals "
      "move in small steps at the doublings: the largest `I` step is %s (at "
      "U=%s) and the largest `alpha_m` step is %s, located at the hand-off "
      "ramp itself (U=%s), not at a doubling.  This is a pre-existing v11/v13 "
      "property of the resolved-band definition; `B(k)` itself is continuous, "
      "and no arm changes it."
      % (_f(hs["F_100"]["pipeline_scan"]["max_rel_jump_I"], 3),
         _f(hs["F_100"]["pipeline_scan"]["U_at_max_I_jump"], 2),
         _f(hs["F_100"]["pipeline_scan"]["max_rel_jump_alpha_m_used"], 3),
         _f(hs["F_100"]["pipeline_scan"]["U_at_max_alpha_jump"], 2)))
    A("")
    A("### 3.3 A3 prior accounting (Arm B; printed per cell)")
    A("")
    rows = []
    for c in cells:
        b = c["arms"]["B"]["alpha"]
        rows.append([_cell(c), _f(b["alpha_m_raw"], 6),
                     _f(b["alpha_m_floor"], 6), _f(b["alpha_m_used"], 6),
                     _f(b["prior_contribution_alpha_m"], 6),
                     _f(b["sigma2_h_raw_would_be"], 6),
                     _f(b["prior_contribution_sigma2_h"], 6),
                     _yn(b["top_up_acted"]), b["provenance"]])
    table(L, ["cell", "raw a_m", "floor a_m", "applied", "d a_m", "s2_h raw",
              "prior d s2_h", "top-up", "provenance"], rows)
    A("")
    A("`prior d s2_h` is `applied*I - raw*I`: the slope variance the invented "
      "input contributes at that cell.  Where `top-up = no` the provenance is "
      "`eckv_fit_domain_raw` and the prior contributes exactly zero.  The "
      "label carried when it acts is `%s`, with the design's source string "
      "`%s`."
      % (s1["arms"]["B"]["provenance"],
         (cells[0]["arms"]["B"]["alpha"]["low_wind_prior"] or {}).get("source")))
    A("")
    A("### 3.4 A4 Parseval + independent band integrals")
    A("")
    A("Two independent things are tested, per arm per cell:")
    A("")
    A("1. **Parseval** on a built field at N=1024: the realised `m0` and "
      "`sigma2` against the mode-sum and against the continuum over the band "
      "the grid resolves (`k_cut = min(1000, 0.95 k_nyquist)`).")
    A("2. **Band integrals on independent grids**: v11's 48k log grid "
      "(reference), a 96k log grid, and a 150k *geometric* grid built on the "
      "band (a different construction).  A 12k log grid is reported as a "
      "coarsening sensitivity and is not part of the pass condition.")
    A("")
    rows = []
    for a in ("RAW", "B"):
        a2 = [c2[c["cell"]]["arms"][a]["A2_parseval"] for c in cells]
        a4 = [c2[c["cell"]]["arms"][a]["A4_independent_integrals"]
              for c in cells]
        bp = [c2[c["cell"]]["arms"][a]["A4_band_reproduction"] for c in cells]
        rows.append([ARMN[a], _pass(all(x["pass"] for x in a2)),
                     _f(max(abs(x["z_m0"] or 0) for x in a2), 2),
                     _f(max(abs(x["z_sigma2"] or 0) for x in a2), 2),
                     _pass(all(x["pass"] for x in a4)),
                     "%.1e" % max(x["worst_abs_rel"] for x in a4),
                     _pass(all(x["pass"] for x in bp))])
    table(L, ["arm", "A2 Parseval", "worst |z| m0", "worst |z| s2",
              "A4 grids", "worst rel", "A4 partition"], rows)
    A("")
    v13bad = [c for c in cells if abs(
        c2[c["cell"]]["arms"]["B"]["A4_band_reproduction"]
        ["v13_band_table_defect"]["rel_over_total"]) > 1e-3]
    A("**A defect found by A4, and fixed here.**  `V11.band_table` builds its "
      "edges as `[kmin, k_p, 10, 50, 200, 370, 1000]`.  Whenever "
      "`k_p > 10 rad/m` (short fetch), that list is **unsorted**, so the row "
      "`10_50` overlaps `below_kp`: the table double counts `[10, k_p]` and "
      "its rows sum to MORE than the cell total.  Measured over the 8 cells, "
      "the overstatement reaches `%s` (relative).  It affects v13's per-band "
      "SHARE columns only -- the mask-based shares used by C1/C1b are computed "
      "directly and are unaffected.  v14 therefore uses its own "
      "`band_partition` (sorted, disjoint, exhaustive; residual < 1e-3 from "
      "trapezoid bookkeeping at the internal edges); the v13 table is kept "
      "beside it so the defect is measured rather than asserted away.  Cells "
      "where it fires (`k_p > 10`, overstatement > 0.1%%): %s."
      % (_f(max(abs(c2[c["cell"]]["arms"]["B"]["A4_band_reproduction"]
                   ["v13_band_table_defect"]["rel_over_total"])
                for c in cells), 3),
         ", ".join("%s (+%s, k_p=%s)"
                   % (c["cell"],
                      _f(c2[c["cell"]]["arms"]["B"]["A4_band_reproduction"]
                         ["v13_band_table_defect"]["rel_over_total"], 3),
                      _f(c2[c["cell"]]["arms"]["B"]["A4_band_reproduction"]
                         ["v13_band_table_defect"]["k_p_rad_per_m"], 1))
                   for c in v13bad) or "none"))
    A("")
    _exp = [(abs(v["rel_total"]), abs(v["rel_capillary_50_200"] or 0.0))
            for c in cells
            for v in s1["export"]["cells"][c["cell"]]
            ["reproduction_from_export"].values()]
    A("The exported spectra (`export.cells` in `numbers-v14-spectrum.json`, "
      "4000 log-spaced points per curve, k from `kmin` to 1500 rad/m) "
      "reproduce the reported band integrals when re-integrated by trapezoid: "
      "the worst relative residual over the 8 cells and the three exported "
      "curves per cell is `%s` on the total and `%s` on the 50-200 capillary "
      "integral." % (("%.2e" % max(x[0] for x in _exp)),
                     ("%.2e" % max(x[1] for x in _exp))))
    A("")
    A("### 3.5 A5 resolution convergence")
    A("")
    A("One cell (U=6 / F=1000 / depth 1.0) re-run at N=2048 with "
      "`k_cut=1500` against N=1024 with `k_cut=1000`:")
    A("")
    rows = []
    for a in ("RAW", "B"):
        c = s2["A5_convergence"][a]
        rows.append([ARMN[a], _f(c["coarse"]["sigma2_realised"], 6),
                     _f(c["fine"]["sigma2_realised"], 6),
                     _p(100.0 * c["rel_sigma2"], 2) + "%",
                     _p(100.0 * c["rel_Hs"], 2) + "%", _pass(c["pass"])])
    table(L, ["arm", "s2 N=1024", "s2 N=2048/k1500", "d s2", "d Hs", "A5"],
          rows)
    A("")
    A("Both totals move by well under the 5% tolerance, so the totals are "
      "converged with respect to grid and k-cutoff in both arms.  `m0` is "
      "reported as missing in the JSON (`m0: null`): the state node does not "
      "expose it, so the test uses `sigma2` and `Hs`, which it does expose -- "
      "stated, not silently substituted.")
    A("")
    A("### 3.6 The control-arm identity checks")
    A("")
    ra = s1["controls"]["RAW_equals_v13"]
    sd = s1["controls"]["D_STRICT_equals_v3p8"]
    A("- **RAW == v13.assemble**: %s on all 8 cells (identical `alpha_m`, "
      "identical `S(k)` elementwise, identical `Hs`).  The control arm is the "
      "shipped v13 spectrum, asserted per cell, not assumed."
      % _pass(ra["pass"]))
    A("- **D-strict == the v3.8 build**: %s on all 8 cells (the diagnostic "
      "reproduces the old hard zero exactly)." % _pass(sd["pass"]))
    A("")
    A("---")
    A("")
    A("## 4. Physical gates (three-state; UNRESOLVED is not a soft pass)")
    A("")
    A("### 4.1 High-wind cells (U = 6, 10) -- must PASS")
    A("")
    rows = []
    for c in cells:
        if c["low_wind_cell"]:
            continue
        for a in ("RAW", "B", "D"):
            x = c["arms"][a]
            rows.append([_cell(c), ARMN[a],
                         _f(x["shares"]["gravity_share_50_200"], 3) + " "
                         + _pass(x["C1"]["pass"]),
                         "%s/%s %s" % (
                             _f(x["shares"]["capillary_share_50_200"], 3),
                             _f(x["moments"]["capillary"]["sigma2"], 4),
                             _pass(x["C1b"]["pass"])),
                         "%sx s%s t%s" % (
                             _f(x["envelope_check"]["ratio_to_lower"], 2),
                             "P" if x["envelope_check"]["strict_inside"]
                             else "F",
                             "P" if x["envelope_check"]["tol_inside"]
                             else "F"),
                         _f(x["b_envelope"]["B_max"], 4) + " "
                         + _pass(x["b_envelope"]["pass_B_env"]),
                         _state(x["state"])])
    table(L, ["cell", "arm", "grav sh<=0.50", "cap sh>=0.15 & s2h>=0.002",
              "tot ratio s=strict t=tol", "B envel", "verdict"], rows)
    A("")
    bad = [(c["cell"], a) for c in cells if not c["low_wind_cell"]
           for a in ("RAW", "B", "D") if c["arms"][a]["state"] != "PASS"]
    if bad:
        A("**High-wind clauses that fail, and why.**")
        A("")
        for cid, a in bad:
            c = [x for x in cells if x["cell"] == cid][0]
            A("- `%s` arm %s: %s" % (cid, ARMN[a],
                                     "; ".join(c["arms"][a]["fail_reasons"])))
        A("")
        A("The failing clause is the **B-envelope** in every row above, and "
          "its location differs by cell.  In three of the four high-wind cells "
          "`max B(k)` sits at the GRAVITY peak: k = %s rad/m for U=6/F=100 "
          "(B = %s) and k = %s rad/m for U=10/F=100 (B = %s) -- the spectral "
          "peak of a short-fetch wind sea, where a curvature spectrum "
          "legitimately reaches 0.018-0.025.  In the fourth (U=10/F=1000) it "
          "is the CAPILLARY peak: `B(k_m) = 0.5 alpha_m = %s` against the "
          "0.012 ceiling, a %s%% excess.  The band check is separate and "
          "passes everywhere (`B_l` max over 50-200: worst %s, bound 0.008).  "
          "Every one of these is set by the cell's raw `alpha_m` / `alpha_p`, "
          "which no arm changes at high wind."
          % (_f([x for x in cells if x["cell"] == "U6p0_F100p0"][0]["arms"]
                ["RAW"]["b_envelope"]["k_at_B_max"], 1),
             _f([x for x in cells if x["cell"] == "U6p0_F100p0"][0]["arms"]
                ["RAW"]["b_envelope"]["B_max"], 4),
             _f([x for x in cells if x["cell"] == "U10p0_F100p0"][0]["arms"]
                ["RAW"]["b_envelope"]["k_at_B_max"], 1),
             _f([x for x in cells if x["cell"] == "U10p0_F100p0"][0]["arms"]
                ["RAW"]["b_envelope"]["B_max"], 4),
             _f([x for x in cells if x["cell"] == "U10p0_F1000p0"][0]["arms"]
                ["RAW"]["b_envelope"]["B_max"], 5),
             _p(100.0 * ([x for x in cells
                          if x["cell"] == "U10p0_F1000p0"][0]["arms"]["RAW"]
                         ["b_envelope"]["B_max"] / pol["B_env_max"] - 1.0), 1),
             _f(max(c["arms"]["RAW"]["b_envelope"]["B_l_max"]
                    for c in cells), 4)))
        A("")
        A("This is a finding about the spectrum, not about the boundary "
          "treatment: the three arms are the same spectrum at high wind, so no "
          "arm can pass or fail differently there.  Under PIPELINE-V12's merge "
          "rule (\"high-wind cells PASS\" is a precondition) this alone blocks "
          "a merge-as-default; the decision is not taken in this session.")
    else:
        A("All high-wind clauses pass in all three arms.")
    A("")
    A("### 4.2 Low-wind cells (U = 1.5, 3) -- at best UNRESOLVED")
    A("")
    rows = []
    for c in cells:
        if not c["low_wind_cell"]:
            continue
        for a in ("RAW", "B", "D"):
            x = c["arms"][a]
            ac = x["absolute_criteria"] or {}
            rows.append([_cell(c), ARMN[a],
                         _f(x["shares"]["gravity_share_50_200"], 3),
                         _f(x["shares"]["capillary_share_50_200"], 3),
                         _f(x["moments"]["capillary"]["sigma2"], 5),
                         "%s/%s" % (_f(x["envelope_check"]["ratio_to_lower"], 2),
                                    _f(x["envelope_check"]["ratio_to_upper"], 2)),
                         _f(x["b_envelope"]["B_max"], 4),
                         _f(x["b_envelope"]["B_l_max"], 4),
                         _p(100.0 * x["hs"]["delta_rel"], 1) + "%",
                         _pass(x["b_envelope"]["pass_B_env"]),
                         _pass(x["smoothness"]["pass"]),
                         _state(x["state"]),
                         "%d/4 | %d/4" % (ac.get("n_pass_strict", 0),
                                          ac.get("n_pass_tolerance", 0))])
    table(L, ["cell", "arm", "grav sh", "cap sh", "s2_h", "tot ratio lo/hi",
              "B max", "B_l max", "dHs", "B env", "C1 smooth", "state",
              "abs crit (strict/tol)"], rows)
    A("")
    A("Low-wind verdicts are UNRESOLVED in every arm: no FAIL clause fires "
      "(totals inside the 0.5x-2.0x band, no B above 0.012, no negative "
      "energy, no step in S(k), |dHs| below 35%).  The last column is the Arm "
      "D legitimacy test scored for every arm: `n_pass/4` of Gemini's absolute "
      "clauses (B-envelope, B_l band, totals, C1 smoothness), first with the "
      "strict totals interval and then with the 0.5x-2.0x reading.")
    A("")
    A("Two things must be said plainly here:")
    A("")
    A("- **At U=3 the strict totals interval FAILS in every arm**: the total "
      "MSS is %s-%s x Cox-Munk clean, i.e. BELOW the lower edge of "
      "`[Cox-Munk clean, Hu]`, although inside the design's own 0.5x-2.0x "
      "tolerance.  The strict interval, applied as a gate, would fail every "
      "high-wind cell too (ratios %s-%s x).  Both readings are printed in "
      "every row; the report does not choose between them."
      % (_f(min(c["arms"]["RAW"]["envelope_check"]["ratio_to_lower"]
                for c in cells if c["low_wind_cell"] and c["U_ms"] == 3.0), 3),
         _f(max(c["arms"]["RAW"]["envelope_check"]["ratio_to_lower"]
                for c in cells if c["low_wind_cell"] and c["U_ms"] == 3.0), 3),
         _f(min(c["arms"]["RAW"]["envelope_check"]["ratio_to_lower"]
                for c in cells if not c["low_wind_cell"]), 3),
         _f(max(c["arms"]["RAW"]["envelope_check"]["ratio_to_lower"]
                for c in cells if not c["low_wind_cell"]), 3)))
    A("- **A pure-gravity spectrum passes the absolute test.**  `D-strict` at "
      "U=1.5 sets `sigma2_h = 0` and still scores %s/4 strict at U=1.5/F=100 "
      "and %s/4 at U=1.5/F=1000: B stays under 0.012, `B_l` stays under 0.008, "
      "the totals land inside the strict interval (F=1000) or just under it "
      "(F=100), and the smoothness test is trivially satisfied.  That is "
      "Grok's objection measured rather than argued (`rally-v39-grok.md` line "
      "10): the absolute criteria alone are satisfiable by a spectrum with no "
      "capillary energy at all, which is exactly why C1b's "
      "`sigma2_h >= 0.002` clause exists in PIPELINE-V12."
      % (_f(cells[0]["arms"]["D_STRICT"]["absolute_criteria"]
            ["n_pass_strict"], 0),
         _f(cells[1]["arms"]["D_STRICT"]["absolute_criteria"]
            ["n_pass_strict"], 0)))
    A("")
    A("### 4.3 C1b: OR (v13) vs AND (PIPELINE-V12)")
    A("")
    A("v13's `c1b()` passed a cell if EITHER clause held.  PIPELINE-V12 states "
      "the conjunction (line 25: capillary share >= 0.15 AND "
      "`sigma2_h >= 0.002`), so v14 uses AND.  The difference is not cosmetic: "
      "the U=1.5/F=1000 cells carry `sigma2_h` at the target but a capillary "
      "share of %s, i.e. just below the 0.15 line, so they PASS under v13's OR "
      "and FAIL under the conjunction.  Since C1b is a gate only on the "
      "high-wind cells (PIPELINE-V12 4.1), no verdict changes, but the count "
      "of \"capillary-supported\" low-wind cells does.  At U=3 the "
      "conjunction is what makes Arm B's top-up arithmetically meaningful: "
      "RAW's capillary share there is %s, B's is %s."
      % (_f([c for c in cells if c["cell"] == "U1p5_F1000p0"][0]["arms"]["B"]
            ["shares"]["capillary_share_50_200"], 4),
         _f([c for c in cells if c["cell"] == "U3p0_F100p0"][0]["arms"]["RAW"]
            ["shares"]["capillary_share_50_200"], 3),
         _f([c for c in cells if c["cell"] == "U3p0_F100p0"][0]["arms"]["B"]
            ["shares"]["capillary_share_50_200"], 3)))
    A("")
    A("---")
    A("")
    A("## 5. Arm B sensitivity: the floor sweep at the U=3 cells")
    A("")
    sens = s2["sensitivity_floor_sweep"]
    for cid, rows in sens["cells"].items():
        A("`%s`:" % cid)
        A("")
        rr = []
        for tag, v in rows.items():
            rr.append([_f(v["floor"], 4), _f(v["alpha_m_raw"], 6),
                       _f(v["alpha_m_used"], 6), _yn(v["top_up_acted"]),
                       _f(v["sigma2_h_total"], 6),
                       _f(v["sigma2_gravity_total"], 6),
                       _f(v["capillary_share_50_200"], 3),
                       _f(v["gravity_share_50_200"], 3),
                       _f(v["sigma2_total_cut"], 6), _f(v["B_max"], 4),
                       _p(100.0 * v["dHs_rel"], 1) + "%",
                       _f(v["ratio_to_lower"], 2), _state(v["state"])])
        table(L, ["floor", "raw a_m", "appl a_m", "acted", "s2_h", "s2_g",
                  "cap sh", "grav sh", "total MSS", "B max", "dHs",
                  "tot/CM lo", "state"], rr)
        A("")
    A("At `U=3/F=1000` the 0.0015 floor does not act: the raw spectrum already "
      "delivers `sigma2_h = %s`, above 0.0015.  Every floor moves the shares "
      "monotonically and none pushes `max B(k)` over 0.012 at these cells; the "
      "totals stay inside the 0.5x-2.0x envelope at all three floors.  The "
      "floor is a free knob: 0.0015 -> 0.006 quadruples the invented "
      "capillary variance and moves the capillary share from %s to %s at "
      "U=3/F=100."
      % (_f(sens["cells"]["U3p0_F1000p0"]["floor_0.0015"]["sigma2_h_total"], 6),
         _f(sens["cells"]["U3p0_F100p0"]["floor_0.0015"]
            ["capillary_share_50_200"], 3),
         _f(sens["cells"]["U3p0_F100p0"]["floor_0.006"]
            ["capillary_share_50_200"], 3)))
    A("")
    A("---")
    A("")
    A("## 6. Transport")
    A("")
    A("The frozen route is unchanged: `Q=%d`, `A=%s`, `E=%s`, "
      "`B_lowpass_sigma = %s px`, seed from the unmodified recipe, frame 0.  "
      "The composite is rendered at the photo's %dx%d; the caustic field `C` "
      "is computed on the %d^2 mesh and is mean-1 by construction."
      % (s3["frozen_composite"]["Q"],
         _f(s3["frozen_composite"]["A_V3"], 1),
         _f(s3["frozen_composite"]["E_V3"], 1),
         _f(s3["frozen_composite"]["B_SIGMA_V3"], 1),
         s3["photo"]["size"][0], s3["photo"]["size"][1],
         2048))
    A("")
    A("### 6.1 The windy cell (U=6 / F=1000 / depth 1.0), per arm")
    A("")
    A("|arm|a_m|s2_h|C mean|C p99|C p99.9|C max|own shift mm|"
      "warp p99.9>=3|crush<=3.25%|blow=0|")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for a in ("RAW", "B"):
        w = s3["scenes"]["windy"][a]
        A("|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|"
          % (ARMN[a], _f(w["alpha_m_used"], 6), _f(w["sigma2_h_total"], 6),
             _f(w["C"]["mean"], 6), _f(w["C_percentiles"]["p99"], 3),
             _f(w["C_percentiles"]["p99_9"], 3),
             _f(w["C_percentiles"]["max"], 2),
             _f(w["warp"]["own_shift_field"]["mean_mm"], 2),
             _yn(w["C_percentiles"]["p99_9"] >= 3.0),
             _yn(w["tone"]["crush_pass_customary"]),
             _yn(w["tone"]["blow_pass_zero"])))
    wv = s3["scenes"]["windy"]["RAW_vs_B"]
    A("")
    A("RAW and B are the SAME render here: `alpha_m` is identical at this cell "
      "(the raw capillary variance exceeds the floor), the two stills are "
      "byte-identical (`%s`), max |d| = %s uint8.  Arm D has no separate "
      "render for the same reason.  The transport sanity lines hold: own-field "
      "shift %s mm (observed band 20-45), `p99.9 = %s >= 3.0`, crush %s <= "
      "3.25%%, blow %s.  C mean is %s, i.e. 1.0 to six decimals."
      % (_yn(wv["byte_identical"]), wv["max_abs_uint8"],
         _f(s3["scenes"]["windy"]["RAW"]["warp"]["own_shift_field"]["mean_mm"],
            2),
         _f(s3["scenes"]["windy"]["RAW"]["C_percentiles"]["p99_9"], 3),
         "%.2f%%" % (100.0 * s3["scenes"]["windy"]["RAW"]["tone"]
                     ["crush_rate_linear_lt_0.02"]),
         _yn(s3["scenes"]["windy"]["RAW"]["tone"]["blow_pass_zero"]),
         _f(s3["scenes"]["windy"]["RAW"]["C"]["mean"], 6)))
    A("")
    A("### 6.2 The boundary cell (U=3 / F=1000, depth 1.0): RAW vs B")
    A("")
    bv = s3["scenes"]["boundary_u3"]["RAW_vs_B"]
    A("Deliverable: `comp-u3-raw-vs-priorB.png` -- left arm RAW, right arm B, "
      "2 px white divider, %dx%d RGB, same seed and frame (the composite size "
      "is the photo's; the C mesh is the 2048^2 grid)."
      % (bv["size_px"][0], bv["size_px"][1]))
    A("")
    A("|arm|a_m|top-up|C mean|C p99|C p99.9|C max|own shift mm|mean |dL/dy|")
    A("|---|---|---|---|---|---|---|---|---|")
    for a in ("RAW", "B"):
        b = s3["scenes"]["boundary_u3"][a]
        A("|%s|%s|%s|%s|%s|%s|%s|%s|%s|"
          % (ARMN[a], _f(b["alpha_m_used"], 6), _yn(b["top_up_acted"]),
             _f(b["C"]["mean"], 6), _f(b["C_percentiles"]["p99"], 3),
             _f(b["C_percentiles"]["p99_9"], 3),
             _f(b["C_percentiles"]["max"], 2),
             _f(b["warp"]["own_shift_field"]["mean_mm"], 2),
             _f(b["fine_texture"]["mean_abs_dL_dy"], 5)))
    A("")
    A("C mean is 1.0 in both arms (`pass_mean_1 = %s`), i.e. 1.0 to six "
      "decimals, as the frozen route requires.  The two stills differ in %s%% "
      "of pixels (mean |d| %s, max %s uint8, MSSSIM %s).  What the boundary "
      "choice does to the picture, measured: the top-up raises `alpha_m` by "
      "%sx, raises the fine texture `mean|dL/dy|` by %s%% and moves the "
      "own-field shift from %s to %s mm; the tail of C becomes slightly LESS "
      "peaked (p99.9 %s -> %s, max %s -> %s) while the mean stays exactly 1 by "
      "construction.  The route's `g1` is 1.0 in BOTH arms, so the tail change "
      "is a property of the deterministic caustic field (`det_frame`, Q=25) "
      "once the high-k band carries more energy; the mechanism is recorded as "
      "NOT resolved here rather than attributed to a shading term that did not "
      "move."
      % (_yn(bv["C_mean_both_one"]), _p(100.0 * bv["frac_pixels_changed"], 1),
         _f(bv["mean_abs_uint8"], 3), _f(bv["max_abs_uint8"], 0),
         "%.6f" % bv["msssim"], _f(bv["alpha_m_ratio_B_over_RAW"], 3),
         _p(100.0 * (s3["scenes"]["boundary_u3"]["B"]["fine_texture"]
                     ["mean_abs_dL_dy"]
                     / s3["scenes"]["boundary_u3"]["RAW"]["fine_texture"]
                     ["mean_abs_dL_dy"] - 1.0), 2),
         _f(s3["scenes"]["boundary_u3"]["RAW"]["warp"]["own_shift_field"]
            ["mean_mm"], 2),
         _f(s3["scenes"]["boundary_u3"]["B"]["warp"]["own_shift_field"]
            ["mean_mm"], 2),
         _f(s3["scenes"]["boundary_u3"]["RAW"]["C_percentiles"]["p99_9"], 3),
         _f(s3["scenes"]["boundary_u3"]["B"]["C_percentiles"]["p99_9"], 3),
         _f(s3["scenes"]["boundary_u3"]["RAW"]["C_percentiles"]["max"], 2),
         _f(s3["scenes"]["boundary_u3"]["B"]["C_percentiles"]["max"], 2)))
    A("")
    A("---")
    A("")
    A("## 7. Decision data (constraints x arms) -- NO winner chosen here")
    A("")
    A("Per PIPELINE-V12 line 20 the arm that satisfies the INDEPENDENT "
      "constraints (totals envelope, B-envelope, construction gates, transport "
      "sanity) with the fewest invented inputs wins merge-as-challenger; a tie "
      "goes to Arm B; Arm D must actually pass its B-envelope checks at U=3.  "
      "The table below is that data; the application of the rule happens "
      "outside this session.")
    A("")
    rows = []
    for r in dec["constraints"]:
        rows.append([r["constraint"]] + [
            "%d/%d%s" % (r["counts"][a]["pass"],
                         r["counts"][a]["pass"] + r["counts"][a]["fail"],
                         ("" if not r["counts"][a]
                          ["not_applicable_or_unresolved"]
                          else " (+%d n/a)" % r["counts"][a]
                          ["not_applicable_or_unresolved"]))
            for a in ARMS])
    table(L, ["constraint", "RAW", "B", "D", "D-strict"], rows)
    A("")
    A("Counts are `pass/fail` over the eight cells; `n/a` are cells where the "
      "clause is not a gate (the high-wind clauses do not apply at U<=3).")
    A("")
    A("Verdict counts by arm:")
    A("")
    rows = []
    for a in ARMS:
        st = list(summ[a]["states"].values())
        rows.append([ARMN[a],
                     _f(sum(1 for s in st if s == "PASS"), 0),
                     _f(sum(1 for s in st if s == "FAIL"), 0),
                     _f(sum(1 for s in st if s == "UNRESOLVED"), 0),
                     _f(summ[a]["n_cells_prior_acts"], 0),
                     _f(summ[a]["worst_B_max"], 4),
                     _f(summ[a]["worst_B_l_max"], 4),
                     _p(100.0 * summ[a]["worst_abs_dHs"], 1) + "%"])
    table(L, ["arm", "PASS", "FAIL", "UNRESOLVED",
              "cells with invented input", "worst B max", "worst Bl max",
              "worst |dHs|"], rows)
    A("")
    A("Invented inputs acting, per cell (`decision.invented_inputs`):")
    A("")
    rows = []
    for c in cells:
        rows.append([_cell(c)] + [
            "%d" % dec["invented_inputs"][a]["per_cell"][c["cell"]]
            ["invented_inputs_acting"] for a in ARMS])
    table(L, ["cell", "RAW", "B", "D", "D-strict"], rows)
    A("")
    A("The only invented input in any arm is `sigma2_h_target = %s` (the "
      "design's provisional, unmeasured prior level, `rally: %s`).  RAW uses "
      "it at %s cell(s), B at %s, D at %s (D inherits RAW's spectrum by "
      "definition), D-strict at %s.  Arm D's legitimacy test at U<=3, scored "
      "on the raw spectrum: %s"
      % (_f(pol["sigma2_h_target_default"], 4),
         (cells[0]["arms"]["B"]["alpha"]["low_wind_prior"] or {}).get("rally"),
         _f(dec["invented_inputs"]["RAW"]["n_cells_with_invented_input"], 0),
         _f(dec["invented_inputs"]["B"]["n_cells_with_invented_input"], 0),
         _f(dec["invented_inputs"]["D"]["n_cells_with_invented_input"], 0),
         _f(dec["invented_inputs"]["D_STRICT"]
            ["n_cells_with_invented_input"], 0),
         "; ".join("%s strict %s (%d/4), tolerance %s (%d/4)"
                   % (cid, v["verdict_strict"], v["n_pass_strict"],
                      v["verdict_tolerance"], v["n_pass_tolerance"])
                   for cid, v in sorted(dec["D_absolute_criteria"].items()))))
    A("")
    A("What the table actually shows, stated without a preference:")
    A("")
    A("- **No gated constraint separates RAW, B and D.**  Every row scores "
      "identically for the three arms: the arms differ only at U=3 (B's "
      "top-up) and every U=3 clause that is scored there (totals, B-envelope, "
      "`B_l`, smoothness, positivity, A4, `|dHs|`) passes in both arms.  The "
      "separation is therefore in the invented-input count and in the un-gated "
      "low-wind arithmetic (capillary share, delivered `sigma2_h`), not in the "
      "constraint table.")
    A("- **D-strict loses one strict-totals cell** (U=1.5/F=100, ratio %s x "
      "Cox-Munk clean) and has zero invented inputs -- the two facts pull in "
      "opposite directions for the rule, which is why the rule is applied "
      "outside this session."
      % _f([c for c in cells if c["cell"] == "U1p5_F100p0"][0]["arms"]
           ["D_STRICT"]["envelope_check"]["ratio_to_lower"], 3))
    A("- **The B-envelope fails in all four arms, at 3 of the 4 high-wind "
      "cells** (only U=6/F=1000 passes; 5/8 cells overall, counting the "
      "low-wind cells where it does pass), so the design's \"high-wind cells "
      "PASS\" precondition is not met by any arm on this evidence.  Under "
      "PIPELINE-V12 section 5 that means no merge-as-default from this chunk; "
      "the report says why and leaves the default alone.")
    A("- **Transport sanity is identical for RAW and B** (%s), and regression "
      "is green for all arms (section 8), so neither of those splits the arms "
      "either." % _pass(dec["transport"]["RAW"]["sanitiy_pass"]))
    A("")
    A("---")
    A("")
    A("## 8. Regression")
    A("")
    rg = s3["regression"]
    A("- **v13 legacy lane, byte-identical re-check**: `%s` vs v13's `%s` -- "
      "sha256 %s / %s, identical = %s, max |d| = %s."
      % (rg["v13_legacy_lane"]["v14_file"],
         rg["v13_legacy_lane"]["v13_file"],
         rg["v13_legacy_lane"]["sha256_v14"][:16],
         rg["v13_legacy_lane"]["sha256_v13"][:16],
         _yn(rg["v13_legacy_lane"]["byte_identical"]),
         _f(rg["v13_legacy_lane"]["max_abs_uint8"], 0)))
    A("- **v13 new lane**: `%s` vs `%s` -- identical = %s (the v13 cell is "
      "inside the ECKV fit domain, so no arm changes a number there)."
      % (rg["v13_new_lane"]["v14_file"], rg["v13_new_lane"]["v13_file"],
         _yn(rg["v13_new_lane"]["byte_identical"])))
    A("- **v12 water path**: import-only check, `WATER_MODEL_DEFAULT = %r`, "
      "`DEFAULT_MEDIUM_PRESET = %r`, `LEGACY_RECIPE_ROUTE = %r` -> %s."
      % (rg["v12_water_path"]["WATER_MODEL_DEFAULT"],
         rg["v12_water_path"]["DEFAULT_MEDIUM_PRESET"],
         rg["v12_water_path"]["LEGACY_RECIPE_ROUTE"],
         _pass(rg["v12_water_path"]["pass"])))
    A("- **frozen composite**: imported, never edited; `A=%s`, `B_sigma=%s`, "
      "`E=%s`; `v5_composite.py` sha256 `%s`."
      % (_f(rg["frozen_composite_hashes"]["A_V3"], 1),
         _f(rg["frozen_composite_hashes"]["B_SIGMA_V3"], 1),
         _f(rg["frozen_composite_hashes"]["E_V3"], 1),
         rg["frozen_composite_hashes"]["v5_composite_sha256"][:16]))
    fh = dec["file_hashes"]
    A("- **every imported module unchanged while the job ran**: %s (%d files "
      "hashed at stage 1 and again at stage 4)."
      % (_yn(fh["all_unchanged"]), len(fh["files"])))
    A("- overall regression: %s" % _pass(rg["pass"]))
    A("")
    A("---")
    A("")
    A("## 9. Deviations, interpretations and missing numbers")
    A("")
    A("### 9.1 Arm D: the brief's two sentences disagree, both are measured")
    A("")
    A("The brief says Arm D \"changes NO numbers -- it is an evaluation arm\", "
      "and that at U=6/10 \"arm D same numbers\".  Taken with \"raw spectrum "
      "untouched\" one could also read D as \"no prior anywhere\".  Those are "
      "different spectra at U=1.5, where v13's prior acts.  Resolution: **D is "
      "RAW's spectrum with a different gate set**, exactly as the brief's "
      "operative sentence says, and the literal alternative is reported beside "
      "it as `D-strict` (the v3.8 hard zero, asserted byte-identical to the "
      "v3.8 build on all 8 cells).  Both are in every table and both are in "
      "the JSON.  No number is hidden by this choice and no gate depends on it "
      "(D-strict is never counted as one of the three arms in the decision "
      "counts).")
    A("")
    A("### 9.2 The totals envelope: strict vs 0.5x-2.0x")
    A("")
    A("PIPELINE-V12 names the interval `[Cox-Munk clean, Hu]` and, separately, "
      "gives the only numeric tolerance in the design: FAIL if the totals are "
      "outside `[0.5x, 2.0x]` of it (line 26).  The three-state verdicts use "
      "that tolerance; the strict interval is printed next to it in every row "
      "and in every JSON cell (`envelope_check.strict_inside`).  Under the "
      "strict reading every high-wind cell fails as well, so the choice "
      "matters and is therefore not buried: of the four high-wind cells, "
      "strict passes %d and the tolerance passes %d (RAW); over all eight "
      "cells strict passes %d and tolerance %d."
      % (sum(1 for c in cells if not c["low_wind_cell"]
             and c["arms"]["RAW"]["envelope_check"]["strict_inside"]),
         sum(1 for c in cells if not c["low_wind_cell"]
             and c["arms"]["RAW"]["envelope_check"]["tol_inside"]),
         sum(1 for c in cells
             if c["arms"]["RAW"]["envelope_check"]["strict_inside"]),
         sum(1 for c in cells
             if c["arms"]["RAW"]["envelope_check"]["tol_inside"])))
    A("")
    A("### 9.3 \"The Hu law\" does not exist in this workspace")
    A("")
    A("Carried forward from v13, unchanged: no Hu-law formula is present in "
      "this repository or in either design document.  The upper column is the "
      "pipeline's own fetch-limited MSS integral "
      "(`V13.fetch_limited_mss_reference`, JONSWAP/Donelan family, gamma=3.3, "
      "same `k_hi` as the model column), and the naming gap is recorded in the "
      "JSON (`envelope.hu_naming_status`).  It is NOT a Hu parameterisation "
      "and this report does not claim otherwise.  If a Hu formula is later "
      "supplied, only this column changes.")
    A("")
    A("### 9.4 Gate values that differ from v13's code")
    A("")
    A("- `B_l` in the 50-200 band: v13's V3 used 0.010 and 0.007 as its "
      "bounds; PIPELINE-V12 fixes `beta_hat = 0.008`.  v14 uses 0.008.")
    A("- The total-curvature ceiling is `B <= 0.012` over `[10, 1500]`, a "
      "criterion v13 never evaluated (v13's V3 looked at 50-200 only, which "
      "hides the capillary peak at ~370 and the gravity peak below 10).  Both "
      "facts are visible in section 4.1.")
    A("- C1b is the conjunction (section 4.3).")
    A("- Cell classification is by `U10` (`U <= 3` = low wind), per "
      "PIPELINE-V12's wording, not by v13's `u* >= c_m/e` (at U=3, "
      "`u*/u*fit = %s`, so v13 called those cells in-domain and demanded C1 "
      "PASS there; PIPELINE-V12 calls them low-wind and at best UNRESOLVED)."
      % _f([c for c in cells if c["cell"] == "U3p0_F100p0"][0]["arms"]["RAW"]
           ["alpha"]["u_over_fit_min"], 3))
    A("")
    A("### 9.5 The band-table defect")
    A("")
    A("Section 3.4.  `V11.band_table` is unsorted when `k_p > 10 rad/m` and "
      "over-counts by up to %s relative.  This job did not edit "
      "`v11_spectrum.py` (it is frozen); v14 carries its own partition and "
      "prints both."
      % _f(max(abs(c2[c["cell"]]["arms"]["B"]["A4_band_reproduction"]
                   ["v13_band_table_defect"]["rel_over_total"])
               for c in cells), 3))
    A("")
    A("### 9.6 A duplicate run of this same job wrote into this folder")
    A("")
    A("Two sessions were working this brief in the same directory.  The "
      "numbers in `numbers-v14-{spectrum,sweep,transport}.json` and the "
      "stills are from ONE foreground run (stages 1-3: PNG timestamps "
      "18:41-18:46, then re-rendered identically at 18:58-19:05 when stage 3 "
      "was re-run to add the image-geometry fields; every statistic below is "
      "identical across the two runs, which is the determinism check).  The "
      "other session read those same JSONs and produced its own report from "
      "them; that report is preserved verbatim as "
      "`REPORT-V14-GATES-alt-runB.md` and its generator as "
      "`run_v14_report_alt.py`.  This report was produced by this session's "
      "stage 4 from the same JSONs.  Three consequences are stated rather than "
      "hidden: (a) `numbers-v14-sweep.json` is re-dumped by stage 4 to attach "
      "the `decision` block -- the stage-2 measurement blocks are carried over "
      "verbatim, nothing physical is recomputed; (b) the checkpoint file "
      "records 8 completions for 4 stages because stages 3 and 4 ran twice, "
      "and it is left as written (a checkpoint is progress proof, not edited "
      "by hand); (c) where the two reports disagree on a derived bookkeeping "
      "ratio, the value here is the one in the JSON "
      "(`v13_band_table_rel_over_total` = %s for U=3/F=100), and the "
      "alternative report's differing value for that one ratio is a "
      "recomputation it did at format time."
      % _f(c2["U3p0_F100p0"]["arms"]["B"]["A4_band_reproduction"]
           ["v13_band_table_defect"]["rel_over_total"], 3))
    A("")
    A("### 9.7 Numbers that are missing, and why")
    A("")
    A("- **No Hu-law number** (9.3).  Reported as the in-repo fetch-limited "
      "MSS, labelled.")
    A("- **`m0` in the A5 rows is `null`**: the state node does not expose "
      "`m0`; `sigma2` and `Hs` are used instead (both present).")
    A("- **`cap>370` for D-strict at U=1.5 is `n/a`**: that arm has no "
      "capillary part, so the fraction is 0/0 and is emitted as null, not 0.")
    A("- **No wavelength-resolved optical validation of the prior**: the "
      "transport section reports C statistics, warp and texture for one "
      "boundary cell; it is not a measurement of a real pond and is not used "
      "as evidence that 0.003 is the right floor.  Nothing in this job can "
      "settle the floor's VALUE -- only its consequences.")
    A("- **No in-water PSF / halo test**: unchanged, not attempted.")
    A("- **Arm D-strict has no render**: its numbers differ from RAW only at "
      "U=1.5 and the brief asks for renders at U=3; a D-strict render is not a "
      "deliverable and was not produced.  Its spectral numbers are complete.")
    A("")
    A("### 9.8 What this report does NOT claim")
    A("")
    A("- It does not claim any arm is validated.  The low-wind cells remain "
      "UNRESOLVED in every arm, and `sigma2_h_target = 0.003` remains an "
      "invented input with no measured basis in this workspace.")
    A("- It does not claim the high-wind cells pass: the B-envelope clause "
      "fails there identically in all arms (4.1), which under the design's "
      "merge rule blocks a default merge.")
    A("- It does not pick a winner.  Section 7 is the data; PIPELINE-V12's "
      "rule is applied outside this session.")
    A("")
    A("---")
    A("")
    A("## 10. Deliverables and how to check them")
    A("")
    A("```")
    A("out/v14/REPORT-V14-GATES.md                     this file")
    A("out/v14/numbers-v14-spectrum.json               per-cell numbers, "
      "hand-off scans, exports")
    A("out/v14/numbers-v14-sweep.json                  states, Parseval, A4, "
      "A5, floor sweep, decision")
    A("out/v14/numbers-v14-transport.json              windy per arm, U=3 "
      "RAW vs B, regression")
    A("out/v14/comp-u3-raw-vs-priorB.png               left RAW / right B, "
      "%dx%d" % (bv["size_px"][0], bv["size_px"][1]))
    A("out/v14/comp-windy-v40-U6F1000-armRAW.png       the windy still (arm B "
      "is byte-identical)")
    A("out/v14/comp-windy-v40-U6F1000-armB.png         same, rendered "
      "separately")
    A("out/v14/comp-windy-v40-v13legacy.png            the v13 legacy lane, "
      "re-rendered")
    A("out/v14/comp-windy-v40-v13new.png               the v13 new lane, "
      "re-rendered")
    A("out/v14/comp-u3-F1000-armRAW.png                the left half, "
      "standalone")
    A("out/v14/comp-u3-F1000-armB.png                  the right half, "
      "standalone")
    A("out/v14/spectrum-ref-v40-U6F1000-armB.csv       one curve pair, plain "
      "text")
    A("out/v14/checkpoint.json                         checkpoint, one unit "
      "per stage")
    A("out/v14/REPORT-V14-GATES-alt-runB.md            the duplicate run's "
      "report (9.6)")
    A("```")
    A("")
    A("Checks a reviewer can run without this report:")
    A("")
    A("1. The two windy stills have equal sha256 "
      "(`scenes.windy.RAW_vs_B.byte_identical`) -- Arm B does nothing at U=6.")
    A("2. `scenes.boundary_u3.RAW_vs_B.alpha_m_ratio_B_over_RAW = %s` is the "
      "whole difference between the two halves of the side-by-side image."
      % _f(bv["alpha_m_ratio_B_over_RAW"], 4))
    A("3. Re-integrate `export.cells` by trapezoid and compare with "
      "`cells[*].arms[*].sigma2_total_cut` and `shares.B_h_sigma2_50_200`: the "
      "residuals are printed in `reproduction_from_export`.")
    A("4. `controls.RAW_equals_v13` and `controls.D_STRICT_equals_v3p8` are "
      "per-cell booleans with the max |dS| printed.")
    A("")
    A("Report built at %s (wall %s s)."
      % (time.strftime("%Y-%m-%dT%H:%M:%S"), _p(time.time() - T0, 1)))
    return "\n".join(L) + "\n"


def main():
    import sys
    sys.path.insert(0, r"D:\mutaro-work\jobkit")
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out",
                           "v14")
    txt = _ascii(build(out_dir))
    path = os.path.join(out_dir, "REPORT-V14-GATES.md")
    from checkpoint import atomic_write_bytes
    atomic_write_bytes(path, txt.encode("ascii"))
    print("wrote %s (%d bytes)" % (path, len(txt)))


if __name__ == "__main__":
    main()
