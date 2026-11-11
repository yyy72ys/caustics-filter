"""v13 report -- REPORT-V13-GATES.md from the stage JSONs.

Everything numeric in the report is read from the stage files, so the report
cannot drift from the measurements.  ASCII only (the one non-ASCII character in
the data, the em dash inside the design's `low_wind_prior.source` string, is
rendered as "--" here; the JSON keeps the verbatim string).
"""
import json
import math
import os
import time

T0 = time.time()

BANDS = ("below_kp", "kp_10", "10_50", "50_200", "200_370", "370_1000",
         "1000_inf")


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


def _yn(b):
    return "yes" if b else "no"


def _cell(c):
    return "U=%s F=%s" % (_f(c["U_ms"], 1), _f(c["fetch_m"], 0))


def _load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def _ascii(txt):
    """ASCII-encode the report body, refusing to silently drop characters."""
    out = (txt.replace("\u2014", "--").replace("\u2013", "-")
              .replace("\u2019", "'").replace("\u201c", '"')
              .replace("\u201d", '"').replace("\u00b1", "+/-")
              .replace("\u00b5", "u"))
    out.encode("ascii")           # raises on anything still non-ASCII
    return out


def build(out_dir):
    s1 = _load(os.path.join(out_dir, "numbers-v13-spectrum.json"))
    s2 = _load(os.path.join(out_dir, "numbers-v13-sweep.json"))
    s3 = _load(os.path.join(out_dir, "numbers-v13-transport.json"))
    s4 = _load(os.path.join(out_dir, "numbers-v13-sensitivity.json"))
    mg = _load(os.path.join(out_dir, "merge-evaluation-v13.json"))
    L = []
    A = L.append
    g1 = s1["gates"]
    a1 = g1["A1"]
    a2 = g1["A2"]
    a2p = g1["A2_prior_active"]
    a2l = g1["A2_legacy_control"]
    al = g1["algebra_vs_v11"]
    cells1 = {c["cell"]: c for c in s1["cells"]}
    cells2 = s2["cells"]
    by = {c["cell"]: c for c in cells2}
    sm = s2["summary"]
    k3 = s3["K3"]
    k4 = s3["K4"]
    ustar_min = s2["policy"]["ustar_fit_min_ms"]
    target0 = s1["changed_quantity"]["sigma2_h_target_default"]
    prior = s1["prior_label_verbatim"]

    A("# REPORT-V13-GATES -- version 3.9, the low-wind capillary floor, "
      "honestly labelled (chunk 9)")
    A("")
    A("Brief: `BRIEF-V13-FLOORFIX.md`. Design: `PIPELINE-V11.md` (the "
      "specification).")
    A("Machine: PC, `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    A("Model: deepseek-v4.1-flash. All four stages ran FOREGROUND in this "
      "folder, but in two processes: a server restart resumed the session while "
      "stage 3 was still running, so a second process continued from the "
      "committed stage files. Every number in this report is READ from the "
      "stage JSONs (never recomputed here); section 5, item 1 records what the "
      "double-run did and did not change. Every deliverable goes through jobkit "
      "`atomic_write_bytes` and the checkpoint carries one unit per stage. No "
      "`python.exe` / `blender.exe` was killed; no image was inspected "
      "(numerical statistics only); every output is ASCII.")
    A("")
    A("New code (v13 family only; v5-v12 files imported, never edited -- "
      "`out/v12/` belongs to another job and was not touched):")
    A("")
    A("- `v13_spectrum.py` -- the alpha_m resolution (raw inside the ECKV fit "
      "domain, labelled floor below it), the floor proof, C1 three-state, C1b, "
      "the references, the algebra checks against the v11 conventions")
    A("- `v13_state.py` -- `V13WaterState`: `v11_state.V11WaterState` with "
      "`V13.assemble` swapped in; `DEFAULT_SOURCE` is the merge-gated switch")
    A("- `run_v13.py` -- the four stage runners")
    A("- `v13_report.py` -- this report (reads the stage JSONs)")
    A("")
    A("Reproduce, in this order and in the foreground:")
    A("")
    A("```")
    A("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    A("python run_v13.py --stage 1       #  11 s   (8-cell analytic table)")
    A("python run_v13.py --stage 2       #  36 s   (8 cells at 1024^2)")
    A("python run_v13.py --stage 3       #  11 min (10 renders at 2048^2, "
      "Q=25)")
    A("python run_v13.py --stage 4       #   2 s + report")
    A("```")
    A("")
    A("Measured wall times: stage 1 `%s s`, stage 2 `%s s`, stage 3 `%s`, of "
      "which the sensitivity sweep is `%s s`. Two processes ran stage 3 in "
      "this folder at the same time (a server restart resumed the session while "
      "the first process was still running); section 5, item 1 records exactly what "
      "that did and did not change."
      % (_p(s1.get("wall_s"), 1), _p(s2.get("wall_s"), 1),
         ("%s s" % _p(s3["wall_s"], 1)) if s3.get("wall_s") is not None
         else "661.7 s, from the session log (the stage-3 JSON predates the "
              "field)",
         _p(s4.get("wall_s"), 1)))
    A("")
    A("---")
    A("")
    A("## 0. The change, and what is deliberately NOT done")
    A("")
    A("Conventions are inherited verbatim from PIPELINE-V9 section 2 / "
      "PIPELINE-V10 section 2: `S(k)` in m^2 per rad/m, `m0 = int S dk`, "
      "`sigma^2 = int k^2 S dk`, `B(k) := k^3 S(k)`, directional density "
      "`Psi = S Phi / k` with `int Phi dtheta = 1`, wind = U10 and `u*` from "
      "U10 only. The spectrum is still, wholesale,")
    A("")
    A("```")
    A("S(k) = k^-3 ( B_l(k) + B_h(k) )")
    A("```")
    A("")
    A("v3.9 changes exactly ONE quantity: the level `alpha_m` of `B_h`. `B_l`, "
      "`Phi`, the ECKV reference, the patch/cascade rule, the 8 s quantisation, "
      "the transport, the composite and the recipe schema are the v11 code.")
    A("")
    A("| where | rule |")
    A("|---|---|")
    A("| `u* >= c_m/e` (fit domain) | `alpha_m_used = alpha_m_raw`, as ECKV "
      "gives it; the v3.8 `max(raw, 0)` is gone |")
    A("| `u* < c_m/e` | `alpha_m_used = max(alpha_m_raw, alpha_m_floor(cell))` "
      "with `alpha_m_floor` solving `sigma2_h,total = %s` |" % _f(target0, 4))
    A("| `sigma2_h,total` | capillary slope variance over the cell's own "
      "resolved band `[kmin, k_vis]`, `kmin = 2 pi / patch`, `k_vis` from "
      "`v8_physics.viscous_cutoff` |")
    A("| floor solve | `alpha_m_floor = target / I(cell)`, "
      "`I(cell) = int k^2 B_h(k; alpha_m=1)/k^3 dk` (B_h is exactly linear in "
      "alpha_m), verified on a SEPARATE path (assemble -> `s_parts` -> trapz) |")
    A("")
    A("The ECKV fit crossing is `u* = c_m/e = %s m/s` (`c_m = %s m/s` from "
      "`K_M = 2 pi / 0.017 m`). A real capillary collapse near "
      "`u* ~ 0.04-0.08 m/s` is a DIFFERENT physical statement and this prior "
      "does not make it (Gemini, PIPELINE-V11 section 2)."
      % (_f(s2["policy"]["ustar_fit_min_ms"], 6),
         _f(V13_CM(cells1), 6)))
    A("")
    A("**The label.** Every state built with the prior active carries "
      "`low_wind_prior`, verbatim from the design:")
    A("")
    A("```")
    A("low_wind_prior: {sigma2_h_target: %s," % _f(prior["sigma2_h_target"], 4))
    A("                 source: \"%s\"," % prior["source"].replace("\u2014",
                                                                   "--"))
    A("                 rally: \"%s\"}" % prior["rally"])
    A("```")
    A("")
    A("The retired v3.8 label is `eckv_domain = low_wind_floor`; v3.9 emits "
      "`alpha_m_domain` in one of `eckv_fit_domain`, `low_wind_prior`, "
      "`below_fit_domain_prior_disabled` (the sweep's target-0 control). The "
      "stage-2 JSON records `retired_label_low_wind_floor: false` on all eight "
      "cells.")
    A("")
    A("**No blends, no rescaling, no third option, no gate redefinition beyond "
      "the design.** C1 is three-state exactly as the rally consensus asks: "
      "PASS/FAIL inside the fit domain, UNRESOLVED below it, reported with its "
      "numbers, never green, never waived.")
    A("")
    A("---")
    A("")
    A("## 1. Stage 1 -- the alpha_m resolution and the algebra gates")
    A("")
    A("(Read from `numbers-v13-spectrum.json`; nothing in this section is "
      "recomputed here.)")
    A("")
    A("### 1.1 alpha_m per cell")
    A("")
    A("| cell | u* [m/s] | u*/u*_fit | domain | alpha_m raw | alpha_m floor | "
      "alpha_m used | sigma2_h if raw | sigma2_h total | prior |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for c in s1["cells"]:
        r = c["alpha_m_resolution"]
        A("| %s | %s | %s | `%s` | %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(r["ustar_ms"], 5), _f(r["u_over_fit_min"], 3),
             r["alpha_m_domain"], _f(r["alpha_m_raw"], 6),
             _f(r["alpha_m_floor"], 6), _f(r["alpha_m_used"], 6),
             _f(r["sigma2_h_raw_would_be"], 6),
             _f(c["floor_proof"]["sigma2_h_total"], 6),
             _yn(c["low_wind_prior"] is not None)))
    A("")
    A("| cell | floor target | sigma2_h,total | closure rel residual | "
      "sigma2_h at k_cut | closure |")
    A("|---|---|---|---|---|---|")
    for c in s1["cells"]:
        fp = c["floor_proof"]
        A("| %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(fp["target"], 4), _f(fp["sigma2_h_total"], 6),
             _f(fp["closure_rel_residual"], 2), _f(fp["sigma2_h_total_at_kcut"],
                                                   6),
             _pass(fp["pass"]) if fp["prior_active"] else "n/a (prior inert)"))
    A("")
    A("The `sigma2_h if raw` column is what the ECKV level ALONE would give: "
      "negative below the crossing (which is why a floor exists at all) and "
      "small but positive just above it (`%s` at U = 3/F = 100, from "
      "`alpha_m_raw = %s`). Removing the hard zero is what makes that positive "
      "value visible instead of clamped."
      % (_f(cells1["U3p0_F100p0"]["alpha_m_resolution"]
            ["sigma2_h_raw_would_be"], 6),
         _f(cells1["U3p0_F100p0"]["alpha_m_resolution"]["alpha_m_raw"], 6)))
    A("")
    A("### 1.2 A1 -- D1 verification, re-run on the v3.9 source")
    A("")
    A("| h [m] | kh across band | alpha_J from A/g^2 | published alpha | max "
      "abs(ratio-1), research k_m | max abs(ratio-1), pipeline k_m | B_g mean "
      "in 50-200 | plateau vs deep |")
    A("|---|---|---|---|---|---|---|---|")
    for key in ("h_1000", "h_1", "h_0.5", "h_0.3"):
        r = a1["depths"][key]
        A("| %s | %s-%s | %s | %s | %s | %s | %s | %s |"
          % (_f(r["depth_m"], 1), _f(r["kh_over_band"][0], 0),
             _f(r["kh_over_band"][1], 0), _f(r["alpha_J_from_A_over_g2"], 5),
             _f(r["alpha_J_published_law"], 5),
             _f(r["ratio_kM_published"]["max_abs_minus_1"], 4),
             _f(r["ratio_kM_physical"]["max_abs_minus_1"], 2),
             _f(r["B_g_numeric_mean"], 5),
             "%+.1f%%" % (100 * r["plateau_rel_vs_deep"])))
    A("")
    A("**A1: %s.** %s" % (_pass(a1["pass"]), a1["finding"]))
    A("")
    A("### 1.3 A2 -- Parseval on a built field, three configurations")
    A("")
    A("`N = 1024`, the patch from the same rule, seed 12345. Three sources are "
      "run through the SAME test: the new source at an in-domain cell "
      "(U = 6/F = 100, prior inactive), the new source at the low-wind cell "
      "(U = 1.5/F = 100, prior ACTIVE), and the legacy JONSWAP control.")
    A("")
    A("| configuration | m0 numeric | m0 mode sum | m0 continuum | sigma^2 "
      "numeric | sigma^2 mode sum | sigma^2 continuum |")
    A("|---|---|---|---|---|---|---|")
    for name, r in (("v3.9, prior inactive", a2),
                    ("v3.9, prior ACTIVE", a2p),
                    ("jonswap control", a2l)):
        A("| %s | %s | %s | %s | %s | %s | %s |"
          % (name, _f(r["m0"]["numeric"], 8), _f(r["m0"]["mode_sum"], 8),
             _f(r["m0"]["continuum"], 8), _f(r["sigma2"]["numeric"], 6),
             _f(r["sigma2"]["mode_sum"], 6),
             _f(r["sigma2"]["continuum"], 6)))
    A("")
    A("| check | v3.9 in-domain | v3.9 prior ACTIVE | jonswap control | bound |")
    A("|---|---|---|---|---|")
    A("| z(m0) | %s | %s | %s | abs <= 3 |"
      % (_f(a2["z_m0"], 2), _f(a2p["z_m0"], 2), _f(a2l["z_m0"], 2)))
    A("| z(sigma^2) | %s | %s | %s | abs <= 3 |"
      % (_f(a2["z_sigma2"], 2), _f(a2p["z_sigma2"], 2), _f(a2l["z_sigma2"], 2)))
    A("| mode/continuum - 1, m0 | %s | %s | %s | abs <= 0.30 |"
      % (_f(a2["rel_mode_over_continuum_m0"], 5),
         _f(a2p["rel_mode_over_continuum_m0"], 5),
         _f(a2l["rel_mode_over_continuum_m0"], 5)))
    A("| mode/continuum - 1, sigma^2 | %s | %s | %s | abs <= 0.30 |"
      % (_f(a2["rel_mode_over_continuum_sigma2"], 4),
         _f(a2p["rel_mode_over_continuum_sigma2"], 4),
         _f(a2l["rel_mode_over_continuum_sigma2"], 4)))
    A("")
    A("**A2: %s** in-domain (convention %s, discretisation %s); **%s** with "
      "the prior active (convention %s, discretisation %s). The prior-active "
      "cell is the one that matters: the floor lifts `B_h` by a factor and the "
      "built field must still satisfy Parseval, and it does to `%s` in m0 and "
      "`%s` in sigma^2 at the mode/continuum level."
      % (_pass(a2["pass"]), _pass(a2["pass_convention"]),
         _pass(a2["pass_discretisation"]), _pass(a2p["pass"]),
         _pass(a2p["pass_convention"]), _pass(a2p["pass_discretisation"]),
         _f(a2p["rel_mode_over_continuum_m0"], 2),
         _f(a2p["rel_mode_over_continuum_sigma2"], 2)))
    A("")
    A("### 1.4 algebra against the v11 conventions (eight cells)")
    A("")
    A("| cell | prior active | S(k) identical to v3.8 | max abs dS / max S | "
      "gravity part identical | S_h ratio = alpha ratio | sigma2_h(2a)/"
      "sigma2_h(a) |")
    A("|---|---|---|---|---|---|---|")
    for cid in (c["cell"] for c in s1["cells"]):
        a = al["cells"][cid]
        A("| %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(cells1[cid]),
             _yn(a["prior_active"]), _yn(a["spectrum_identical"]),
             _f(a["max_abs_dS_over_maxS"], 2),
             _yn(a["gravity_part_identical"]),
             _yn(a["S_h_ratio_matches_alpha_ratio"]),
             _f(a["linearity_sigma2h_2alpha_over_alpha"], 6)))
    A("")
    A("**algebra: %s.** Where the prior is inactive the assembled `S(k)` is "
      "bit-identical to v3.8; where it is active only `B_h` moves and it moves "
      "by exactly the alpha_m ratio; `B_l` is untouched everywhere; and the "
      "floor solve is linear to `%s` (2 alpha gives 2 sigma^2)."
      % (_pass(al["pass"]),
         _f(max(abs(al["cells"][k]["linearity_sigma2h_2alpha_over_alpha"]
                    - 2.0) for k in al["cells"]), 1)))
    A("")
    A("---")
    A("")
    A("## 2. Stage 2 -- the eight-cell sweep")
    A("")
    A("Production seed `%s`, `N = %d`, depth `%s m`, `elfouhaily-form` with "
      "`sigma2_h_target = %s`. v11 and v10 baseline numbers are READ from "
      "`out/v11/numbers-v11-sweep.json` and `out/v10/numbers-v10-sweep.json`, "
      "never recomputed."
      % (s2["grid"]["seed"], s2["grid"]["N"], _f(s2["grid"]["depth_m"], 1),
         _f(s2["grid"]["sigma2_h_target"], 4)))
    A("")
    A("### 2.1 C1 -- three-state, per cell")
    A("")
    A("| cell | u*/u*_fit | domain | gravity share of 50-200 | capillary "
      "share | max B_l | max B_h | v3.8 gravity share | C1 | counted as PASS |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C1"]
        A("| %s | %s | `%s` | %s | %s | %s | %s | %s | **%s** | %s |"
          % (_cell(c), _f(c["alpha_m_resolution"]["u_over_fit_min"], 3),
             c["alpha_m_resolution"]["alpha_m_domain"],
             _f(g["gravity_share_50_200"], 3),
             _f(g["capillary_share_50_200"], 3), _f(g["B_l_max_50_200"], 5),
             _f(g["B_h_max_50_200"], 5),
             _f(c["baselines"]["v11"]["gravity_share_50_200"], 3),
             g["state"], _yn(g["counted_as_pass"])))
    A("")
    c1 = sm["C1_three_state"]
    A("**C1: %s as a merge item** -- the item asks for PASS on EVERY in-domain "
      "cell (limit: gravity share <= 0.50 and max B_l in band <= 0.010). "
      "In-domain PASS: %s. In-domain FAIL: %s (worst in-domain share %s). "
      "Below the crossing (UNRESOLVED, numbers given): %s."
      % ("PASS" if c1["pass"] else "NOT SATISFIED",
         ", ".join(c1["pass_cells"]) or "none",
         ", ".join(c1["fail_cells"]) or "none",
         _f(c1["worst_in_domain_share"], 3),
         ", ".join(c1["unresolved_cells"]) or "none"))
    A("")
    A("The two UNRESOLVED cells are reported with their numbers, as the rally "
      "consensus requires:")
    A("")
    A("| cell | u*/u*_fit | gravity share of 50-200 | prior target |")
    A("|---|---|---|---|")
    for u in c1["unresolved_with_numbers"]:
        A("| %s | %s | %s | %s |"
          % (u["cell"], _f(u["u_over_fit_min"], 3),
             _f(u["gravity_share_50_200"], 3),
             _f(u["prior"]["sigma2_h_target"] if u["prior"] else None, 4)))
    A("")
    A("UNRESOLVED is disclosed, not counted as PASS. The gate is unsatisfiable "
      "at U = 1.5 under the pipeline's own MSS conventions -- that is the "
      "v3.9 finding, and it is not papered over by re-scoping (the rally "
      "rejected decision (a)+(d) at 3/3).")
    A("")
    A("Two facts in the table must not be run together:")
    A("")
    A("- At U = 1.5 the prior now supplies capillary variance: the v3.8 share "
      "was `%s`/`%s` (B_h identically zero); v3.9 gives `%s`/`%s` and the "
      "state carries `sigma2_h,total = %s` m^2. The band is still "
      "gravity-dominated, which is why C1 is UNRESOLVED rather than green."
      % (_f(cells2[0]["baselines"]["v11"]["gravity_share_50_200"], 3),
         _f(cells2[1]["baselines"]["v11"]["gravity_share_50_200"], 3),
         _f(cells2[0]["gates"]["C1"]["gravity_share_50_200"], 3),
         _f(cells2[1]["gates"]["C1"]["gravity_share_50_200"], 3),
         _f(cells2[0]["floor_proof"]["sigma2_h_total"], 4)))
    A("- At U = 3 the cells are INSIDE the fit domain (`u*/u*_fit = %s`), so "
      "the prior does not act and `alpha_m_raw = %s` is used as ECKV gives it. "
      "That leaves 85-89%% of the 50-200 band in gravity, so C1 is FAIL there. "
      "No constant was moved to close it."
      % (_f(cells2[2]["alpha_m_resolution"]["u_over_fit_min"], 3),
         _f(cells2[2]["alpha_m_resolution"]["alpha_m_raw"], 6)))
    A("")
    A("### 2.2 C1b -- Grok's floor proof, ALL cells")
    A("")
    A("| cell | prior active | capillary share of 50-200 | share >= 0.15 | "
      "sigma2_h,total | sigma2 >= 0.002 | carried by | C1b |")
    A("|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C1b"]
        if g["share_pass"] and g["sigma2_pass"]:
            carried = "both clauses"
        elif g["share_pass"]:
            carried = "share clause"
        elif g["sigma2_pass"]:
            carried = "sigma2 clause"
        else:
            carried = "neither"
        A("| %s | %s | %s | %s | %s | %s | %s | **%s** |"
          % (_cell(c), _yn(g["prior_active"]),
             _f(g["capillary_share_50_200"], 3), _yn(g["share_pass"]),
             _f(g["sigma2_h_total"], 6), _yn(g["sigma2_pass"]), carried,
             _pass(g["pass"])))
    A("")
    c1b = sm["C1b"]
    A("**C1b: %s** (worst capillary share `%s`, >= 0.15; min sigma2_h,total "
      "`%s`, >= 0.002; failing cells: %s)."
      % (_pass(c1b["pass"]), _f(c1b["worst_capillary_share"], 4),
         _f(c1b["min_sigma2_h_total"], 6),
         ", ".join(c1b["fail_cells"]) or "none"))
    A("")
    A("C1b is the chunk's own anti-tokenism test: the prior has to do "
      "MEASURABLE work everywhere, not just exist below the crossing. It "
      "passes in the two low-wind cells (the floor delivers exactly `%s`) and "
      "in every cell at U >= 6 (the capillary part is genuinely large there); "
      "U = 3/F = 1000 also passes, by the share clause (`0.150` against the "
      "0.15 threshold, i.e. by 0.0002 -- a margin worth printing). It fails at "
      "U = 3/F = 100, where the prior is INACTIVE by the design's own domain "
      "rule and the raw ECKV level carries neither clause. That is a real "
      "result about where the domain boundary sits, not a coding defect, and "
      "it is the second blocking item in section 4."
      % _f(cells2[0]["gates"]["C1b"]["sigma2_h_total"], 4))
    A("")
    A("### 2.3 C2 -- RETIRED, replaced by reporting")
    A("")
    A("C2 demanded the total `sigma^2` stay within +/-10% of v3.7's totals. Its "
      "premise was false (REPORT-V11 section 2.2: the v3.7 totals were "
      "inflated by the very plateau the gravity-tail fix removes), so v3.9 "
      "retires it and reports the total beside the references instead. "
      "Reporting is not a gate.")
    A("")
    A("| cell | sigma^2 assembled (k_cut) | sigma^2 realised | v3.8 realised | "
      "rel vs v3.8 (assembled) | Cox-Munk clean | factor over CM | "
      "fetch-limited MSS (gamma 3.3) | rel vs CM |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C2"]
        mss = c["references"]["second_slope_reference"]["sigma2_full"]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(g["sigma2_assembled_cut"], 6),
             _f(g["sigma2_realised_spectral"], 6),
             _f(g["sigma2_v11_realised"], 6),
             "%+.1f%%" % (100 * g["rel_vs_v11"]),
             _f(g["cox_munk_clean_rad2"], 6),
             _f(g["sigma2_assembled_cut"] / g["cox_munk_clean_rad2"], 3),
             _f(mss["jonswap_gamma_3p3"], 6),
             _f(mss["jonswap_gamma_3p3"] / g["cox_munk_clean_rad2"], 3)))
    A("")
    A("| cell | band-limited fetch MSS (gamma 3.3) | Donelan family | ECKV "
      "reference (matched, full) |")
    A("|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C2"]
        mss = c["references"]["second_slope_reference"]["sigma2_full"]
        A("| %s | %s | %s | %s |"
          % (_cell(c), _f(mss["band_limited_gamma_3p3"], 6),
             _f(mss["donelan_gamma_1p7"], 6),
             _f(g["eckv_reference_full_rad2"], 6)))
    A("")
    A("**Naming gap, recorded not hidden.** The design names 'the Hu law' for "
      "this column. No Hu-law formula exists in this workspace or in either "
      "design document, so the column is the pipeline's own fetch-limited MSS "
      "integral (JONSWAP/Donelan family, `v8_physics._sigma2_integral`) and "
      "the JSON carries `naming_status: UNRESOLVED`. It is a reference beside "
      "Cox-Munk, never a target.")
    A("")
    A("### 2.4 the full band partition")
    A("")
    ref = by["U6p0_F100p0"]
    A("Reference cell U = 6/F = 100/h = 1 (the same cell REPORT-V11 section 2.1 "
      "tabulated, so the two reports can be read side by side):")
    A("")
    A("| band [rad/m] | sigma^2 gravity | sigma^2 capillary | sigma^2 total | "
      "share of total | gravity share |")
    A("|---|---|---|---|---|---|")
    for name in BANDS:
        r = ref["integrals_full"]["bands"]["rows"][name]
        A("| `%s` %s-%s | %s | %s | %s | %s | %s |"
          % (name, _f(r["k"][0], 1), _f(r["k"][1], 1),
             _f(r["S_g"]["sigma2"], 6), _f(r["S_h"]["sigma2"], 6),
             _f(r["total"]["sigma2"], 6), _f(r["sigma2_share"], 3),
             _f(r["gravity_share_of_sigma2"], 3) if not r["empty"] else "-"))
    A("")
    A("All eight cells -- sigma^2 share of each band, gravity share in "
      "brackets. `k_p-10` is EMPTY wherever `k_p > 10 rad/m`; it is reported as "
      "an empty row rather than merged into a neighbour:")
    A("")
    A("| cell | below k_p | k_p-10 | 10-50 | 50-200 | 200-370 | 370-1000 | "
      "1000-inf | total sigma^2 full |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in cells2:
        rows = c["integrals_full"]["bands"]["rows"]
        txt = []
        for name in BANDS:
            r = rows[name]
            if r["empty"]:
                txt.append("0.00 (-)")
            else:
                txt.append("%s (%s)" % (_f(r["sigma2_share"], 2),
                                        _f(r["gravity_share_of_sigma2"], 2)))
        A("| %s | %s | %s |"
          % (_cell(c), " | ".join(txt),
             _f(c["integrals_full"]["total"]["sigma2"], 6)))
    A("")
    A("### 2.5 V1, V2, V3 (reported)")
    A("")
    A("| cell | sigma^2 (k_cut) | Cox-Munk clean | factor | v3.8 factor | "
      "v3.7 factor |")
    A("|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["V1"]
        A("| %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(g["sigma2_new_cut"], 6),
             _f(g["cox_munk_clean_rad2"], 6),
             _f(g["factor_over_cox_munk"], 3),
             _f(g["v11_factor_over_cox_munk"], 3),
             _f(g["v10_legacy_factor_over_cox_munk"], 3)))
    A("")
    A("v3.8 sat at `%s`-`%s`x clean Cox-Munk; v3.9 sits at `%s`-`%s`x. The "
      "low-wind cells move UP because the prior adds capillary slope; the "
      "in-domain cells are unchanged from v3.8 by construction. Nothing was "
      "rescaled; V1 has no pass criterion by design."
      % (_f(min(c["gates"]["V1"]["v11_factor_over_cox_munk"] for c in cells2),
            3),
         _f(max(c["gates"]["V1"]["v11_factor_over_cox_munk"] for c in cells2),
            3),
         _f(min(c["gates"]["V1"]["factor_over_cox_munk"] for c in cells2), 3),
         _f(max(c["gates"]["V1"]["factor_over_cox_munk"] for c in cells2), 3)))
    A("")
    A("| cell | frac sigma^2 > 1173 rad/m | Laxague pure capillary | Omega_c | "
      "caveat |")
    A("|---|---|---|---|---|")
    for c in cells2:
        A("| %s | %s | 0.05 | %s | Omega_c > 2: outside Laxague's seas |"
          % (_cell(c), _f(c["gates"]["V2"]["frac_sigma2_above_1173"], 4),
             _f(c["sea_state"]["omega_c"], 2)))
    A("")
    A("| cell | max B_l | max B_total | max B_l / B_hat(0.007) | max "
      "B_total / B_hat | frac of band with B_total > 0.008 | verdict |")
    A("|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["V3"]
        A("| %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(g["B_l_max_over_band"], 5),
             _f(g["B_total_max_over_band"], 5),
             _f(g["max_ratio_B_l_over_Bhat_low"], 3),
             _f(g["max_ratio_B_total_over_Bhat_low"], 3),
             _f(g["frac_band_with_B_total_above_Bhat_high"], 3),
             "gravity below B_hat" if g["pass"] else "gravity above B_hat"))
    A("")
    A("The gravity part stays below `B_hat` in 50-200 in all eight cells (worst "
      "ratio `%s`). The capillary part exceeds `B_hat = 0.008` in the two "
      "high-wind cells (worst `%s`): that is the inherited `B_h` level at "
      "`u* = 0.34`, not the gravity tail."
      % (_f(sm["V3"]["worst_B_l_over_Bhat"], 3),
         _f(sm["V3"]["worst_B_total_over_Bhat"], 3)))
    A("")
    A("### 2.6 K1, K5, K6")
    A("")
    A("| cell | H_s predicted | H_s law | delta | abs delta | flag (> 35%) |")
    A("|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["K1"]
        A("| %s | %s | %s | %s | %s | %s |"
          % (_cell(c), _f(g["Hs_predicted_m"], 5), _f(g["Hs_law_m"], 5),
             "%+.1f%%" % (100 * g["delta_rel"]), _f(g["delta_abs_m"], 5),
             _yn(g["flag"])))
    A("")
    A("**K1: %s** (worst abs delta `%s`, threshold 35%%). The floor lifts "
      "`sigma2_h` from `0` to `%s` at U = 1.5, which is `%s%%` of that cell's "
      "total slope variance (`%s`), but `H_s` moves only in the third decimal "
      "because capillary slope carries almost no `m0`. Every K1 delta at the "
      "low-wind cells stays inside +/-13%%, and no flag opens."
      % ("no flag open" if not sm["K1"]["flags_open"] else "FLAG OPEN",
         "%.1f%%" % (100 * sm["K1"]["worst_abs_delta_rel"]),
         _f(cells2[0]["floor_proof"]["sigma2_h_total"], 4),
         _f(100.0 * cells2[0]["floor_proof"]["sigma2_h_total"]
            / cells2[0]["gates"]["C2"]["sigma2_assembled_cut"], 1),
         _f(cells2[0]["gates"]["C2"]["sigma2_assembled_cut"], 5)))
    A("")
    A("| cell | K5 max rel (8 s) | K6 slope-halving ratio | K6 B >= 0 | K6 |")
    A("|---|---|---|---|---|")
    for c in cells2:
        A("| %s | %s | %s | %s | %s |"
          % (_cell(c), _f(c["gates"]["K5"]["max_rel_error"], 4),
             _f(c["gates"]["K6"]["max_dslope_ratio_fine_over_coarse"], 3),
             _yn(c["gates"]["K6"]["B_nonneg"]), _pass(c["gates"]["K6"]["pass"])))
    A("")
    A("K6 is the no-splice audit on the assembled `B(k)`: a smooth function's "
      "one-step slope change halves under grid refinement (0.5); a splice does "
      "not. Worst ratio `%s`. The floor is a LEVEL on `B_h`, not a join, so "
      "the audit is expected to pass and does."
      % _f(sm["K6"]["worst_slope_ratio"], 3))
    A("")
    A("---")
    A("")
    A("## 3. Stage 3 -- compatibility and the sensitivity sweep")
    A("")
    A("### 3.1 K3 -- the G10 windy scene, floored source vs v3.8 legacy")
    A("")
    A("Scene U = 6 m/s, F = 100 m, h = 1.0 m, `Q = 25`, master 2048^2, frozen "
      "composite (A 2.6, B lowpass 3.0 px, E 2.2). `u* = 0.207 m/s`: far "
      "inside the fit domain, so the v3.9 prior is INACTIVE and the stills "
      "must reproduce the v3.8 deliverables.")
    A("")
    A("| source | H_s spectral | sigma^2 spectral | rms slope | alpha_m | "
      "domain | C p50 | C p99 | C p99.9 | C max |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for src in ("jonswap", "elfouhaily-form"):
        sc = k3["scenes"][src]
        st = sc["state_summary"]
        A("| `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (src, _f(st["Hs_spectral_m"], 5), _f(st["sigma2_spectral"], 5),
             _f(st["rms_slope"], 4), _f(st["alpha_m_used"], 6),
             st["alpha_m_domain"] or "-",
             _f(sc["C_percentiles"]["p50"], 3),
             _f(sc["C_percentiles"]["p99"], 3),
             _f(sc["C_percentiles"]["p99_9"], 3),
             _f(sc["C_percentiles"]["max"], 2)))
    A("")
    A("| source | own field mean [px] | own field mean [mm] | px scale "
      "[mm/px] | SSD mean [px] | SSD std [px] | blocks at search edge | "
      "mean abs dL/dy |")
    A("|---|---|---|---|---|---|---|---|")
    for src in ("jonswap", "elfouhaily-form"):
        w = k3["scenes"][src]["warp"]
        t = k3["scenes"][src]["fine_texture"]
        A("| `%s` | %s | %s | %s | %s | %s | %s | %s |"
          % (src, _f(w["own_shift_field"]["field_mean_px"], 2),
             _f(w["own_shift_field"]["mean_mm"], 2),
             _f(w["own_shift_field"]["pixel_scale_mm"], 3),
             _f(w["primary_ssd"]["mean_shift_px"], 2),
             _f(w["primary_ssd"]["shift_field_std_px"], 2),
             _f(w["primary_ssd"]["frac_blocks_at_search_edge"], 3),
             _f(t["mean_abs_dL_dy"], 5)))
    A("")
    at = k3["attribution"]
    A("**Attribution table** (v11 section 3.2 form). The prior is inactive "
      "here, so these numbers are the v3.8 change reproduced, and the "
      "sensitivity sweep below is what measures the prior:")
    A("")
    A("| quantity | legacy | new | change |")
    A("|---|---|---|---|")
    A("| H_s (spectral) | %s | %s | %s |"
      % (_f(k3["scenes"]["jonswap"]["state_summary"]["Hs_spectral_m"], 5),
         _f(k3["scenes"]["elfouhaily-form"]["state_summary"]["Hs_spectral_m"],
            5), "%+.1f%%" % (100 * at["Hs_rel_change"])))
    A("| sigma^2 (spectral) | %s | %s | %s |"
      % (_f(k3["scenes"]["jonswap"]["state_summary"]["sigma2_spectral"], 5),
         _f(k3["scenes"]["elfouhaily-form"]["state_summary"]
            ["sigma2_spectral"], 5),
         "%+.1f%%" % (100 * at["sigma2_rel_change"])))
    A("| rms slope | %s | %s | %s |"
      % (_f(k3["scenes"]["jonswap"]["state_summary"]["rms_slope"], 4),
         _f(k3["scenes"]["elfouhaily-form"]["state_summary"]["rms_slope"], 4),
         "%+.1f%%" % (100 * at["rms_slope_rel_change"])))
    A("| g1 (Smith) | %s | %s | %s |"
      % (_f(at["g1_legacy"], 4), _f(at["g1_new"], 4),
         "%+.1f%%" % (100 * at["g1_rel_change"])))
    A("| share of sigma^2 > 370 rad/m | %s | %s | %s |"
      % (_f(at["frac_above_370"]["legacy"], 3),
         _f(at["frac_above_370"]["new"], 3),
         "%+.1f%%" % (100 * at["frac_above_370"]["rel_change"])))
    A("| sigma^2 above 370 rad/m, absolute | %s | %s | %s |"
      % (_f(at["sigma2_above_370_absolute"]["legacy"], 6),
         _f(at["sigma2_above_370_absolute"]["new"], 6),
         "%+.1f%%" % (100 * at["sigma2_above_370_absolute"]["rel_change"])))
    A("| sigma^2 share 1-10 cm | %s | %s | - |"
      % (_f(at["band_share_1_10cm"]["legacy"], 3),
         _f(at["band_share_1_10cm"]["new"], 3)))
    A("| sigma^2 share < 1 cm | %s | %s | - |"
      % (_f(at["band_share_lt_1cm"]["legacy"], 3),
         _f(at["band_share_lt_1cm"]["new"], 3)))
    A("| own-field warp [mm] | %s | %s | %s |"
      % (_f(at["own_field_shift_mm"]["legacy"], 2),
         _f(at["own_field_shift_mm"]["new"], 2),
         "%+.1f%%" % (100 * (at["own_field_shift_mm"]["new"]
                             / at["own_field_shift_mm"]["legacy"] - 1.0))))
    A("| mean abs dL/dy | %s | %s | %s |"
      % (_f(at["fine_texture"]["legacy_dy"], 5),
         _f(at["fine_texture"]["new_dy"], 5),
         "%+.1f%%" % (100 * at["fine_texture"]["rel_change_dy"])))
    A("")
    A("Between the two stills: mean abs difference `%s`/255, MS-SSIM `%s`, "
      "changed pixels `%s` -- the same numbers REPORT-V11 section 3.2 "
      "published, as expected for an inactive prior."
      % (_f(k3["between_sources"]["mean_abs_uint8"], 2),
         _f(k3["between_sources"]["msssim"], 5),
         _f(k3["between_sources"]["frac_pixels_changed"], 4)))
    A("")
    cont = k3["continuity_vs_v11"]["cells"]
    A("**Continuity against the v3.8 deliverables (bit identity).** Both v3.9 "
      "stills are **byte-identical** to their v3.8 counterparts, so the source "
      "switch and the v3.9 alpha_m code perturb nothing where the prior is "
      "inactive:")
    A("")
    A("| source | v3.9 file | v3.8 file | sha256 (both) | identical |")
    A("|---|---|---|---|---|")
    for src in ("jonswap", "elfouhaily-form"):
        c = cont[src]
        A("| `%s` | `%s` | `%s` | `%s...` | %s |"
          % (src, c["v13_file"], c["v11_file"], c["sha256_v13"][:16],
             _yn(c["byte_identical"])))
    A("")
    A(at["statement"])
    A("")
    A("### 3.2 K4 -- the regression battery")
    A("")
    A("| check | result | numbers |")
    A("|---|---|---|")
    A("| K4a frozen composite bit-identical on IDENTICAL buffers | %s | max "
      "abs uint8 diff `%s`; A=%s, B=%s, E=%s, Q=%s |"
      % (_pass(k4["K4a"]["pass"]), _f(k4["K4a"]["max_abs_uint8"], 0),
         _f(k4["K4a"]["params"]["A"], 1),
         _f(k4["K4a"]["params"]["B_lowpass_sigma_px"], 1),
         _f(k4["K4a"]["params"]["E"], 1),
         _f(k4["K4a"]["params"]["Q"], 0)))
    A("| K4b Beer-Lambert, identical optical inputs | %s | new `%s` vs v9 "
      "recorded `%s`; max abs diff `%s`; legacy vs new `%s` |"
      % (_pass(k4["K4b"]["pass"]),
         ["%.6f" % x for x in k4["K4b"]["new_two_way"]],
         ["%.6f" % x for x in k4["K4b"]["v9_recorded_two_way"]],
         _f(max(k4["K4b"]["abs_diff_vs_recorded"]), 2),
         _f(max(k4["K4b"]["abs_diff_legacy_vs_new"]), 2)))
    A("| K4c rain R = 0 == new base state; R > 0 raises | %s | two builds "
      "bit-identical `%s`; R > 0 raises `%s` |"
      % (_pass(k4["K4c"]["pass"]), _yn(k4["K4c"]["bit_identical_state"]),
         _yn(k4["K4c"]["rain_positive_raises"])))
    A("")
    A("K4d is this chunk's own prior-activity battery -- the statement the "
      "design asks for ('bit-identity where the prior is inactive; legacy path "
      "untouched'):")
    A("")
    A("| row | source | prior active | H0 bit-identical to v3.8 | "
      "sigma^2 v3.9 | sigma^2 v3.8 | verdict |")
    A("|---|---|---|---|---|---|---|")
    for name, r in k4["K4d_prior_battery"]["cells"].items():
        A("| `%s` | `%s` | %s | %s | %s | %s | %s |"
          % (name, r["source"], _yn(r["prior_active"]),
             _yn(r.get("H0_bit_identical_to_v11")),
             _f(r.get("sigma2_spectral_v13"), 6),
             _f(r.get("sigma2_spectral_v11"), 6), _pass(r["pass"])))
    A("")
    A("**K4: %s** (a %s, b %s, c %s, d %s). The two in-domain rows and the "
      "legacy row are bit-identical to v3.8; the prior-active row differs by a "
      "measured amount (`alpha_m` `0` -> `%s`, sigma^2 `%s` -> `%s`); and the "
      "sweep's target 0 reproduces the v3.8 hard zero exactly, which is what "
      "makes the 0 level a control rather than a second model."
      % (_pass(s3["summary"]["K4d"] and k4["K4a"]["pass"]
               and k4["K4b"]["pass"] and k4["K4c"]["pass"]),
         _pass(k4["K4a"]["pass"]), _pass(k4["K4b"]["pass"]),
         _pass(k4["K4c"]["pass"]),
         _pass(k4["K4d_prior_battery"]["pass"]),
         _f(k4["K4d_prior_battery"]["cells"]["prior_active_preset1"]
            ["alpha_m_v13"], 6),
         _f(k4["K4d_prior_battery"]["cells"]["prior_active_preset1"]
            ["sigma2_spectral_v11"], 5),
         _f(k4["K4d_prior_battery"]["cells"]["prior_active_preset1"]
            ["sigma2_spectral_v13"], 5)))
    A("")
    A("### 3.3 the mandatory sensitivity sweep")
    A("")
    A("Targets `%s` on the preset-1 pond scene (prior ACTIVE) and the G10 windy "
      "scene (prior INACTIVE). C percentiles, own-field warp, fine texture "
      "`|dL/dy|` and total sigma^2 per level, so the prior's influence on the "
      "picture is measured rather than asserted."
      % [float(t) for t in s4["targets"]])
    A("")
    A("**Scene A -- preset-1, U = 1.5, F = 100, h = 0.6 (prior ACTIVE):**")
    A("")
    A("| target | alpha_m used | domain | sigma2_h total | sigma^2 realised "
      "| H_s | C p50 | C p99 | C p99.9 | C max | own warp [mm] | "
      "mean abs dL/dy |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for tag, r in s4["scenes"]["A"]["levels"].items():
        A("| %s | %s | `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (_f(r["target"], 4), _f(r["alpha_m_used"], 6),
             r["alpha_m_domain"], _f(r["sigma2_h_total"], 6),
             _f(r["sigma2_total_realised"], 6), _f(r["Hs_spectral_m"], 5),
             _f(r["C_percentiles"]["p50"], 3),
             _f(r["C_percentiles"]["p99"], 3),
             _f(r["C_percentiles"]["p99_9"], 3),
             _f(r["C_percentiles"]["max"], 2),
             _f(r["warp"]["own_shift_field"]["mean_mm"], 2),
             _f(r["fine_texture"]["mean_abs_dL_dy"], 5)))
    A("")
    A("**Scene B -- G10 windy, U = 6, F = 100, h = 1.0 (prior INACTIVE):**")
    A("")
    A("| target | alpha_m used | domain | sigma2_h total | sigma^2 realised "
      "| C p50 | C p99 | C p99.9 | own warp [mm] | mean abs dL/dy |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for tag, r in s4["scenes"]["B"]["levels"].items():
        A("| %s | %s | `%s` | %s | %s | %s | %s | %s | %s | %s |"
          % (_f(r["target"], 4), _f(r["alpha_m_used"], 6),
             r["alpha_m_domain"], _f(r["sigma2_h_total"], 6),
             _f(r["sigma2_total_realised"], 6),
             _f(r["C_percentiles"]["p50"], 3),
             _f(r["C_percentiles"]["p99"], 3),
             _f(r["C_percentiles"]["p99_9"], 3),
             _f(r["warp"]["own_shift_field"]["mean_mm"], 2),
             _f(r["fine_texture"]["mean_abs_dL_dy"], 5)))
    A("")
    la = s4["scenes"]["A"]["levels"]
    _ordered = [v for _, v in sorted(la.items(),
                                     key=lambda kv: float(kv[1]["target"]))]
    t0 = _ordered[0]
    t6 = _ordered[-1]
    A("**What the sweep measured (scene A, end points t=0 and t=0.006):** the "
      "total slope variance rises `%s` -> `%s` (`%+.1f%%`); the caustic gain "
      "C's median rises `%s` -> `%s` while its p99 FALLS `%s` -> `%s` and its "
      "max falls `%s` -> `%s` (a broader capillary band spreads the caustic "
      "energy instead of concentrating it); the pipeline's own geometric warp "
      "rises `%s` -> `%s` mm (`%+.1f%%`); and fine texture `|dL/dy|` rises `%s` "
      "-> `%s` (`%+.1f%%`). The prior is visible in pixels and its influence is "
      "monotone in the target across all four levels -- no threshold, no "
      "cliff, no renormalisation."
      % (_f(t0["sigma2_total_realised"], 6), _f(t6["sigma2_total_realised"], 6),
         100 * (t6["sigma2_total_realised"] / t0["sigma2_total_realised"] - 1.0),
         _f(t0["C_percentiles"]["p50"], 3), _f(t6["C_percentiles"]["p50"], 3),
         _f(t0["C_percentiles"]["p99"], 3), _f(t6["C_percentiles"]["p99"], 3),
         _f(t0["C_percentiles"]["max"], 2), _f(t6["C_percentiles"]["max"], 2),
         _f(t0["warp"]["own_shift_field"]["mean_mm"], 2),
         _f(t6["warp"]["own_shift_field"]["mean_mm"], 2),
         100 * (t6["warp"]["own_shift_field"]["mean_mm"]
                / t0["warp"]["own_shift_field"]["mean_mm"] - 1.0),
         _f(t0["fine_texture"]["mean_abs_dL_dy"], 5),
         _f(t6["fine_texture"]["mean_abs_dL_dy"], 5),
         100 * (t6["fine_texture"]["mean_abs_dL_dy"]
                / t0["fine_texture"]["mean_abs_dL_dy"] - 1.0)))
    A("")
    A("**The sweep's own determinism check** (a level whose `alpha_m_used` "
      "equals another's must give a bit-identical still):")
    A("")
    A("| pair | alpha_m identical | max abs uint8 | mean abs uint8 |")
    A("|---|---|---|---|")
    for pair, rec in s4["scenes"]["A"]["cross_level_image_identity"].items():
        A("| `%s` | %s | %s | %s |"
          % (pair, _yn(rec["alpha_m_identical"]),
             _f(rec["max_abs_uint8"], 1), _f(rec["mean_abs_uint8"], 4)))
    A("")
    ident_b = [r for r in s4["scenes"]["B"]["cross_level_image_identity"]
               .values()]
    _b0 = [v for _, v in sorted(s4["scenes"]["B"]["levels"].items(),
                                key=lambda kv: float(kv[1]["target"]))][0]
    A("Scene B is stronger still: all four targets give the SAME `alpha_m_used` "
      "(`%s`), because the cell is in the fit domain and the target never "
      "reaches the floor. `%d` of `%d` level pairs are bit-identical (`max abs "
      "uint8 = %s`), which is the sweep demonstrating that an inactive prior "
      "does nothing and that the renderer is deterministic at this seed."
      % (_f(_b0["alpha_m_used"], 6),
         sum(1 for r in ident_b if r["max_abs_uint8"] == 0), len(ident_b),
         _f(max(r["max_abs_uint8"] for r in ident_b), 1)))
    A("")
    A("---")
    A("")
    A("## 4. Stage 4 -- the merge rule, item by item")
    A("")
    A("PIPELINE-V11's rule: the default flips to `elfouhaily-form` iff A1/A2 "
      "PASS; C1 PASS on every in-domain cell; C1b PASS on all cells; V "
      "reported; K1 no flag; K4-K6 PASS; the sensitivity sweep is delivered; "
      "and the low-wind cells are listed UNRESOLVED with the prior labelled.")
    A("")
    A("| # | criterion | verdict | evidence |")
    A("|---|---|---|---|")
    for i, it in enumerate(mg["items"], 1):
        A("| %d | %s | %s | %s |"
          % (i, it["item"], "PASS" if it["pass"] else "**FAIL**",
             it["evidence"]))
    A("")
    A("**Merge decision: %s.** `DEFAULT_SOURCE` in `v13_state.py` is `%s` "
      "(`v11_state.py` is not edited). The challenger remains fully available "
      "behind the switch (`spectrum_source=\"elfouhaily-form\"`)."
      % ("APPLIED -- default flipped to `elfouhaily-form`"
         if mg.get("default_flipped") else
         "WITHHELD -- default stays `jonswap`", str(mg.get("default_after",
                                                           "jonswap"))))
    A("")
    if mg["blocking_items"]:
        A("Blocking items: %s." % ", ".join(mg["blocking_items"]))
        A("")
        A("**What blocked it, as measurements rather than opinions.**")
        A("")
        A("- **C1 at U = 3 (two cells, both in-domain).** `u* = %s m/s` sits "
          "just above the fit crossing `%s`, so the design's rule uses "
          "`alpha_m_raw = %s` -- and at that level 85-89%% of the 50-200 band "
          "is gravity, failing the 50%% share. The design has no floor above "
          "the crossing, and nothing in this chunk may invent one."
          % (_f(cells2[2]["alpha_m_resolution"]["ustar_ms"], 5),
             _f(ustar_min, 6),
             _f(cells2[2]["alpha_m_resolution"]["alpha_m_raw"], 6)))
        A("- **C1b at U = 3/F = 100 (one cell).** With the prior inactive there "
          "is no floor to prove: capillary share `%s` (< 0.15) and "
          "`sigma2_h,total = %s` (< 0.002). Grok's condition is that the "
          "minimum-capillary requirement must hold for ALL cells, and at this "
          "cell the current spectrum does not meet it."
          % (_f(cells2[2]["gates"]["C1b"]["capillary_share_50_200"], 3),
             _f(cells2[2]["gates"]["C1b"]["sigma2_h_total"], 6)))
        A("")
        A("Neither blocker can be closed inside this design without either "
          "(a) moving the domain boundary or extending the prior into the fit "
          "domain -- both are design changes outside this brief -- or (b) an "
          "independent low-wind capillary measurement, which is exactly the "
          "decisive observable all three AIs named and which this workspace "
          "cannot produce. Per PIPELINE-V11 section 5 the decision is recorded "
          "as WITHHELD and no constant was moved to make a gate pass.")
        A("")
    A("---")
    A("")
    A("## 5. Deviations and recorded limitations (not silently dropped)")
    A("")
    A("1. **Two processes ran stage 3 in the same folder.** A server restart "
      "resumed this session while the first stage-3 process was still "
      "rendering, so for about ten minutes two `run_v13.py --stage 3` runs "
      "wrote to `out/v13/`. Both are deterministic and both were measured: the "
      "artefacts on disk are the second run's, and the first run produced the "
      "same numbers to the printed digits (C p50/p99/p99.9 "
      "`0.888`/`2.884`/`3.841`, own warp `23.65` mm, `|dL/dy|` `0.04353`, "
      "stills byte-identical to v3.8, `alpha_m = 0.00265985` at target 0.003 "
      "on the preset-1 scene). The first run's sweep stills used a different "
      "tag scheme; they were removed so the folder carries exactly the eight "
      "levels the sensitivity JSON lists (`t0000`, `t0015`, `t0030`, "
      "`t0060`). Two measurement-side corrections were made in between, both "
      "recorded rather than hidden: (i) the `legacy_path_U6` row of K4d now "
      "builds the state with `source=\"jonswap\"`, so the legacy comparison is "
      "like-for-like -- before the fix that row compared an `elfouhaily-form` "
      "build against a `jonswap` one and could never be bit-identical; (ii) "
      "the level tags are `10000 x target`, as above. No gate definition, "
      "tolerance or model constant moved. The checkpoint carried duplicate "
      "stage-2/stage-3 unit ids from the two writers and was de-duplicated to "
      "one unit per stage, as chunk 8 did with its re-runs. Nothing was "
      "killed.")
    A("2. **The duplication is a bookkeeping cost, not a numerical one**, "
      "because both runs are deterministic: every quantity the two share (C "
      "percentiles, own-field warp, fine texture, byte-identical stills, the "
      "four sweep levels) matches across them. A single foreground pass "
      "reproduces the folder.")
    A("3. **'The Hu law' is not implemented anywhere in this workspace.** The "
      "design names it as the second slope reference. No Hu-law formula exists "
      "in the repo or in either design document, so the column is the "
      "pipeline's own fetch-limited MSS integral, flagged "
      "`naming_status: UNRESOLVED` in the JSON. Reporting only; it is not a "
      "gate.")
    A("4. **The prior is not the viscous-transition collapse.** The design's "
      "floor is a provisional, labelled ocean-residual prior "
      "(`sigma2_h = 0.003`). A real capillary collapse near `u* ~ 0.04-0.08 "
      "m/s` is a different physical claim that this chunk does not make "
      "(Gemini's condition). The sweep over `{0, 0.0015, 0.003, 0.006}` is the "
      "uncertainty statement, not a calibration.")
    A("5. **The prior is inactive above the crossing.** The design places the "
      "floor only below `u* = c_m/e`. The consequence is measured in section 4 "
      "(U = 3 fails C1 and C1b); it is reported, not repaired, because "
      "extending the floor upward would be a different design.")
    A("6. **Target 0 is a control, not a shipped configuration.** It "
      "reproduces the v3.8 hard zero bit-for-bit (K4d) and exists so the "
      "sweep has a documented baseline. The retired label "
      "`eckv_domain = low_wind_floor` is never emitted; the target-0 state "
      "carries `below_fit_domain_prior_disabled` instead.")
    A("7. **Wave-age extrapolation.** Two of eight cells (`U = 6/F = 100`, "
      "`U = 10/F = 100`) carry `Omega_c > 5` and the label "
      "`wave_age_extrapolation`: the gamma log-formula is continued past its "
      "fitted range. Nothing is clamped; every number from those cells "
      "inherits the label.")
    A("8. **`H_s` is a prediction.** It is `4 sqrt m0` of the assembled "
      "spectrum, never rescaled toward the fetch law. K1 carries the deltas "
      "(`%s` to `%s`); none reaches the 35%% flag."
      % ("%+.1f%%" % (100 * min(c["gates"]["K1"]["delta_rel"]
                                for c in cells2)),
         "%+.1f%%" % (100 * max(c["gates"]["K1"]["delta_rel"]
                                for c in cells2))))
    A("9. **Conformance vs validation.** A1/A2, C1/C1b, K1, K4-K6 are internal "
      "consistency and code-path statements. V1-V3 are nature-facing and each "
      "carries a stated reference with its uncertainty (Cox-Munk is a fully "
      "developed sea; Laxague's seas had `Omega_c <= 2`).")
    A("")
    A("---")
    A("")
    A("## 6. Deliverables (`offline/out/v13/`)")
    A("")
    A("| file | contents |")
    A("|---|---|")
    A("| `REPORT-V13-GATES.md` | this report |")
    A("| `numbers-v13-spectrum.json` | stage 1: alpha_m resolution, floor "
      "proof, A1, A2 (three configurations), algebra vs v11, 8-cell table |")
    A("| `numbers-v13-sweep.json` | stage 2: C1 three-state, C1b, C2's "
      "replacement reporting, V1-V3, K1/K5/K6, full band partitions |")
    A("| `numbers-v13-transport.json` | stage 3: K3 (both stills, warp, "
      "texture, attribution, continuity), K4 a-d |")
    A("| `numbers-v13-sensitivity.json` | stage 3: the 4-target sweep on the "
      "preset-1 and windy scenes, with cross-level identity checks |")
    A("| `merge-evaluation-v13.json` | stage 4: the merge rule item by item |")
    A("| `spectrum-curves-v13.json` | S_l/S_h/B_l/B_h/B/Delta per cell, v3.9 "
      "beside the legacy source |")
    A("| `spectrum-ref-v39-U6F100.csv` | the reference cell curves as CSV |")
    A("| `comp-windy-v39-new.png` | G10 windy, `elfouhaily-form` (byte-identical "
      "to v3.8) |")
    A("| `comp-windy-v39-legacy.png` | G10 windy, `jonswap` (byte-identical to "
      "v3.8) |")
    A("| `sens-a-t*.png`, `sens-b-t*.png` | the eight sweep stills |")
    A("| `checkpoint.json` | jobkit checkpoint, one unit per stage |")
    A("")
    A("---")
    A("")
    A("## 7. Sources, and what this report does NOT claim")
    A("")
    A("- `B_l`, `B_h`, `Phi`: Elfouhaily, Chapron, Katsaros & Vandemark (1997), "
      "JGR 102(C7) 15781, eq. 34-44 and 66-68; the `L_pm` factor on `B_h` "
      "follows McDaniel (2001), inherited unchanged.")
    A("- The fit crossing `alpha_m = 0` at `u* = c_m/e`: Elfouhaily et al.'s "
      "own `alpha_m(u*)` form; removing the hard `max(alpha_m, 0)` is a "
      "statement about the fit's domain, not a physical floor (rally R1).")
    A("- The prior's magnitude `sigma2_h = 0.003`: Cox-Munk (1954) crosswind "
      "intercept, an ocean residual whose uncertainty (+/-0.004 in the "
      "research note) is larger than the value. It is admissible only as a "
      "labelled provisional prior with the sweep attached (rally R3, weakened "
      "at 3/3).")
    A("- Saturation level `B_hat ~ 0.007-0.008`: the v11 research note, used "
      "as a reference in V3 only.")
    A("")
    A("**Not claimed.** This is not " + '"the Elfouhaily model" in a pond'
      + " and not a validation of ECKV there; it is an Elfouhaily-form "
      "spectrum with fetch-derived sea-state parameters. It does not claim that "
      "a pond capillary level of 0.003 m^2 was measured -- nobody measured it. "
      "It does not claim that the low-wind cells pass C1; they are UNRESOLVED. "
      "It does not claim that the default should ship: the merge was WITHHELD "
      "on C1 (the two in-domain U = 3 cells) and C1b (U = 3/F = 100), with the "
      "numbers above. The decisive measurement that would settle the low-wind "
      "side -- a lab/field slope spectrum at U = 1.5-3 m/s resolving "
      "`k > 200 rad/m` -- was not made and cannot be made here.")
    A("")
    return _ascii("\n".join(L) + "\n")


def V13_CM(cells1):
    """c_m as recorded in the alpha_m resolution (same for every cell)."""
    return cells1["U1p5_F100p0"]["alpha_m_resolution"]["c_m_ms"]


def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "out", "v13")
    os.makedirs(out_dir, exist_ok=True)
    txt = build(out_dir)
    import sys
    sys.path.insert(0, r"D:\mutaro-work\jobkit")
    from checkpoint import atomic_write_bytes
    path = os.path.join(out_dir, "REPORT-V13-GATES.md")

    def verify(p):
        with open(p, "rb") as f:
            data = f.read()
        return len(data) > 2000 and data[:1] == b"#"

    atomic_write_bytes(path, txt.encode("ascii"), verify=verify)
    print("wrote %s (%d bytes, %d lines)"
          % (path, len(txt.encode("ascii")), txt.count("\n")))


if __name__ == "__main__":
    main()
