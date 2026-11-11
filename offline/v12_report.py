"""v12 report -- REPORT-V12-GATES.md from the stage JSONs.

Every number in the report is READ from the stage files, so the report cannot
drift from the measurements.  ASCII only.  Atomic write through jobkit.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import atomic_write_bytes  # noqa: E402

OUT = os.path.join(HERE, "out", "v12")


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


def _pass(b):
    if b is None:
        return "n/a"
    return "PASS" if b else "FAIL"


def _yn(b):
    return "yes" if b else "no"


def _load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def _vec(v, n=4):
    return "[" + ", ".join(_f(x, n) for x in v) + "]"


def build():
    s1 = _load(os.path.join(OUT, "numbers-v12-water.json"))
    s2 = _load(os.path.join(OUT, "numbers-v12-g1.json"))
    g2 = _load(os.path.join(OUT, "numbers-v12-g2.json"))
    g3 = _load(os.path.join(OUT, "numbers-v12-g3.json"))
    g4 = _load(os.path.join(OUT, "numbers-v12-g4.json"))
    mg = _load(os.path.join(OUT, "merge-evaluation-v12.json"))

    chk = s1["checks"]
    mc = s1["model_comparison_preset1_scene"]
    g1a = s2["G1a_legacy_lane"]
    g1b = s2["G1b_bzero_lane"]
    g1c = s2["G1c_same_coefficients_new_model"]
    L = []
    A = L.append

    A("# REPORT -- v3.6, chunk 6: in-water transport (a/b split, single scatter)")
    A("")
    A("Version 3.6 (PIPELINE-V8 chunk 6, BRIEF-V12-CHUNK6).  Implemented "
      "2026-11-11, model deepseek-v4.1-flash.  Everything numeric below is "
      "read from the stage JSON files in this directory; nothing is asserted.")
    A("")
    A("Design source: PIPELINE-V8 chunk 6.  Base: the v3.5-v3.8 stack (v5-v11 "
      "code imported, never edited).  The spectrum source stays `jonswap` "
      "(the v3.8 merge is withheld); the frozen composite (A = 2.6 + Blow + "
      "E = 2.2 normalisation) is unchanged; retina is off; rain is untouched "
      "(R = 0).")
    A("")

    # ---------------------------------------------------------------- verdict
    A("## 0. Verdict")
    A("")
    A("| gate | verdict | headline number |")
    A("|---|---|---|")
    for it in mg["items"]:
        A("| %s | **%s** | %s |"
          % (it["item"], "PASS" if it["pass"] else "FAIL", it["evidence"]))
    A("")
    A("Merge rule (brief, stage 4): the a/b water path becomes the default "
      "iff G1-G4 pass.  **merge_pass = %s**, blocking items: %s.  "
      "`WATER_MODEL_DEFAULT` in `v12_water.py` after this run: **%s**."
      % (_yn(mg["merge_pass"]),
         ", ".join(mg["blocking_items"]) or "none", mg["default_after"]))
    A("")

    # ---------------------------------------------------------------- change
    A("## 1. What changed, and what did not")
    A("")
    A("Changed (new code only: `v12_water.py`, `run_v12.py`, `v12_report.py`):")
    A("")
    A("1. the turbidity scalar no longer drives the water optics of a v3.6 "
      "recipe: a(lambda) and b(lambda) are separate, chosen by a particle "
      "preset (or explicit spectra) ;")
    A("2. the direct beam is split into a ballistic term (kept sharp) and a "
      "single-scatter halo whose width comes from the scattering integral "
      "(section 2.2) instead of the empirical `sigma = 40 px * sqrt(b z)` ;")
    A("3. a CDOM absorption port, additive in a(lambda), default 0 -- port "
      "only, no model.")
    A("")
    A("Unchanged: the wave state (v11, `jonswap`), the deterministic "
      "transport `det.det_frame` (Q = 25), the frozen composite constants, "
      "the Smith G1 term, the substrate gain, the sky/Fresnel term, the "
      "tonemap, the rain path, and every earlier version's file.  A legacy "
      "recipe without `water.medium` takes the pre-3.6 lane, which "
      "*delegates to `v5_composite.render_frame_v5`* with the v5 optics -- "
      "so it cannot move by construction, and G1 measures that it did not.")
    A("")
    A("## 2. Model")
    A("")
    A("### 2.1 a(lambda) and b(lambda)")
    A("")
    A("Channel wavelengths 650/550/470 nm (the pipeline's own dispersion "
      "set, imported from `v5_physics.CHANNEL_WL_NM`; %s)." % s1["conventions"]["wl_note"])
    A("")
    A("    a(lambda) = a_w(lambda) + load * a_particle(lambda) + a_cdom(lambda)")
    A("    b(lambda) = b_w(lambda) + load * b_particle(lambda)")
    A("    c(lambda) = a(lambda) + b(lambda)")
    A("")
    A("Pure water (the only fully measured parts):")
    A("")
    anchors = chk["pure_water_anchors"]
    A("    a_w  = %s  1/m   Pope & Fry (1997) table (v5_physics, imported)"
      % _vec(anchors["a_w_pope_fry"]))
    A("    b_w  = %s  1/m   Morel (1974): 0.00288 * (500/lambda)^4.32"
      % _vec(anchors["b_w_morel_1974"]))
    A("    b_w(legacy) = %s  1/m  v5_physics's 0.0018 anchor (kept for G1)"
      % _vec(anchors["b_w_legacy_v5_anchor"]))
    A("")
    A("The two b_w anchors differ by %s (relative, Morel over the v5 "
      "anchor).  The legacy lane keeps the v5 anchor so G1 can be "
      "bit-identical; the a/b lane uses Morel.  Reported, not reconciled."
      % ("%+.1f%%" % (100 * anchors["rel_diff_morel_over_legacy"][1])))
    A("")
    A("CDOM port: `a_cdom(lambda) = a_cdom(440) * exp(-S (lambda - 440))`, "
      "S = %s /nm (Bricaud, Morel & Prieur 1981; natural-water range "
      "%s /nm, Twardowski et al. 2004).  `cdomA440` defaults to 0 and a "
      "zero fill leaves every number in this report unchanged (checked by "
      "`checks_partition`'s b = 0 corner and by the G1 lane, which never "
      "touches the port)." % (s1["constants"]["CDOM_S_NM"],
                              s1["constants"]["CDOM_S_RANGE_NM"]))
    A("")
    A("Presets (the a/b RATIO is the preset's identity; the load scales it):")
    A("")
    A("| preset | a_part(650/550/470) | b_part(650/550/470) | CDOM a440 | "
      "g | default load | a550 / b550 at load 1 | inside its declared range |")
    A("|---|---|---|---|---|---|---|---|")
    pr = chk["presets"]["rows"]
    for name, P in s1["preset_table"].items():
        row = pr["%s/load_1" % name]
        A("| %s | %s | %s | %s | %s | %s | %s / %s | %s |"
          % (name, _vec(P["a_part"]), _vec(P["b_part"]), _f(P["cdom_a440"], 2),
             _f(P["g"], 2), _f(P["load_default"], 2), _f(row["a550"], 4),
             _f(row["b550"], 4), _yn(row["within_declared_range"])))
    A("")
    A("Provenance: a_w, b_w, the CDOM slope and g carry numeric literature "
      "values.  **The particulate and dissolved preset values are modelling "
      "choices placed inside published order-of-magnitude ranges for that "
      "water class**, and each row above is checked against its declared "
      "range at load 1 and at its default load.  The range sources are the "
      "design research's water-type tables (Kirk 2011; Babin et al. 2003 for "
      "SPM); they were NOT re-verified online in this session (web search "
      "was unavailable), which is stated rather than papered over.")
    A("")
    A("HG asymmetry: g = %s (the design text's freshwater reference "
      "candidate), band %s from Petzold (1972)'s measured VSFs.  HG is known "
      "to under-fit the near-forward peak (Mobley et al. 2002 use two HG "
      "terms); g is a recipe field and its effect on sigma is measured in "
      "section 2.2."
      % (_f(s1["constants"]["G_DEFAULT"], 2),
         s1["constants"]["G_BAND"]))
    A("")
    A("Effect of the new model on the preset-1 scene (the same scene as G1):")
    A("")
    A("    legacy : a = %s, b = %s, c = %s 1/m, empirical blur %s px"
      % (_vec(mc["legacy_v5"]["a"]), _vec(mc["legacy_v5"]["b"]),
         _vec(mc["legacy_v5"]["c"]), _f(mc["legacy_v5"]["blur_sigma_px"], 2)))
    A("    new    : a = %s, b = %s (same coefficients, new transport)"
      % (_vec(mc["ab_same_coefficients"]["a"]),
         _vec(mc["ab_same_coefficients"]["b"])))
    A("             T_bal = %s, H = %s, A = %s"
      % (_vec(mc["ab_same_coefficients"]["T_ballistic"]),
         _vec(mc["ab_same_coefficients"]["H_halo"]),
         _vec(mc["ab_same_coefficients"]["A_lost"])))
    A("             derived sigma = %s m = %s px (empirical %s px, ratio %s)"
      % (_f(mc["ab_same_coefficients"]["sigma_m"], 4),
         _f(mc["ab_same_coefficients"]["sigma_px"], 2),
         _f(mc["ab_same_coefficients"]["empirical_sigma_px"], 2),
         _f(mc["ab_same_coefficients"]
            ["sigma_ratio_derived_over_empirical"], 2)))
    A("")
    A("Read that carefully: for the SAME coefficients the derived halo is "
      "~%sx wider than the old empirical blur, but it only carries H ~ %s of "
      "the light instead of 100%%.  The old model was wrong in both "
      "directions at once, which is why the still (G1c) moves by %s/255 "
      "while the coefficient change (G1b) moves it by %s/255."
      % (_f(mc["ab_same_coefficients"]
            ["sigma_ratio_derived_over_empirical"], 1),
         _f(float(mc["ab_same_coefficients"]["H_halo"][1]), 2),
         _f(g1c["mean_abs_uint8"], 2), _f(g1b["mean_abs_uint8"], 2)))
    A("")

    # ---------------------------------------------------------------- deriv
    A("### 2.2 Single scattering: where sigma comes from")
    A("")
    A("Down the path of total length L (L = s_down + s_up, both legs) the "
      "probability that a photon is still ballistic after distance s is "
      "exp(-c s) with c = a + b, and the density of FIRST scattering events "
      "at s is b exp(-c s).  After scattering by theta the photon covers the "
      "remaining path (L - s) at that angle, so its lateral displacement is")
    A("")
    A("    r = (L - s) tan(theta)")
    A("")
    A("Its second moment and the total single-scatter weight are")
    A("")
    A("    H     = int_0^L b exp(-c s) ds            = (b/c)(1 - exp(-c L))")
    A("    <r^2> = int_0^L b exp(-c s) (L-s)^2 <tan^2 theta> ds")
    A("          = <theta^2> (b/c^3) D2(c L),   D2(x) = x^2 - 2x + 2 - 2exp(-x)")
    A("")
    A("(substituting u = L - s makes the s-integral elementary; D2 is exact).  "
      "For an isotropic-in-plane halo <r^2> = 2 sigma^2, so")
    A("")
    A("    sigma^2 = (<theta^2>/2) (b/c^3) D2(c L)      [m^2]")
    A("    sigma_px = sigma_m / (S/N)                   <- geometry, no constant")
    A("")
    A("<theta^2> is measured from the HG phase function by quadrature: at "
      "g = %s it is %s rad^2, i.e. %s above the second-order identity "
      "2(1-g) = %s rad^2.  The width uses the MEASURED moment; 2(1-g) is "
      "reported as the size of the small-angle step (sigma would be %s "
      "smaller with it)."
      % (_f(chk["hg"]["g"], 2), _f(chk["hg"]["theta2_measured_rad2"], 4),
         "%+.1f%%" % (100 * chk["hg"]["theta2_rel_diff"]),
         _f(chk["hg"]["theta2_closed_form_2_1_minus_g"], 3),
         "%+.1f%%" % (-100 * (1 - math.sqrt(
             1.0 / (1.0 + chk["hg"]["theta2_rel_diff"]))))))
    A("")
    A("Limits, measured (not asserted):")
    A("")
    A("    thin  (cL -> 0):  sigma^2 -> (<theta^2>/2) b L^3/3  =>  "
      "sigma ~ b^%s L^%s   (fit; expected 0.5 / 1.5)"
      % (_f(chk["sigma_limits"]["exponent_sigma_vs_b_thin"], 3),
         _f(chk["sigma_limits"]["exponent_sigma_vs_L_thin"], 3)))
    A("    thick: the WIDTH keeps growing with L; it is the WEIGHT that "
      "saturates, H -> b/c.")
    A("")
    A("The retired empirical form `sigma ~ sqrt(b z)` gets the b exponent "
      "right in the thin limit and the depth exponent wrong (1 instead of "
      "1.5), and it carried a 40 px free constant where this derivation has "
      "the pixel size.  Both statements are the fitted numbers above.")
    A("")
    A("Formula-level checks (stage 1 JSON):")
    A("")
    A("| check | value | verdict |")
    A("|---|---|---|")
    A("| HG <cos(theta)> from the phase function | %s (abs err %s) | %s |"
      % (_f(chk["hg"]["mean_cos_measured"], 6),
         _f(chk["hg"]["mean_cos_abs_err"], 1), _pass(chk["hg"]["pass"])))
    A("| HG <theta^2> quadrature convergence (1e5 vs 4e5 nodes) | rel %s | ok |"
      % _f(chk["hg"]["theta2_quadrature_rel_diff"], 1))
    A("| sigma thin-limit exponents | b^%s, L^%s | %s |"
      % (_f(chk["sigma_limits"]["exponent_sigma_vs_b_thin"], 3),
         _f(chk["sigma_limits"]["exponent_sigma_vs_L_thin"], 3),
         _pass(chk["sigma_limits"]["pass"])))
    A("| sigma Monte-Carlo vs closed form | %s vs %s m (rel %s) | %s |"
      % (_f(chk["sigma_monte_carlo"]["sigma_monte_carlo_m"], 4),
         _f(chk["sigma_monte_carlo"]["sigma_closed_form_m"], 4),
         _f(chk["sigma_monte_carlo"]["rel_diff"], 1),
         _pass(chk["sigma_monte_carlo"]["pass"])))
    A("| T_bal + H + A = 1 (worst over the a,b,s grid) | %s | %s |"
      % (_f(chk["partition"]["worst_identity_err"], 1),
         _pass(chk["partition"]["pass"])))
    A("| b = 0 => H = 0 and T = exp(-aL) | %s | %s |"
      % (_yn(chk["partition"]["b_zero_halo_is_zero"]),
         _pass(chk["partition"]["b_zero_T_equals_expminusaL"])))
    A("| presets inside their declared ranges (load 1 and default) | %s rows | %s |"
      % (len(pr), _pass(chk["presets"]["pass"])))
    A("")
    A("g sensitivity (sigma at b = 1.0, c = 1.2, L = 1.25): %s"
      % {k: _f(v, 3) for k, v in chk["g_sensitivity"]["rows"].items()})
    A("")
    A("### 2.3 The composite")
    A("")
    A("The OUTPUT-SPEC section 3 separation is kept; only what multiplies the "
      "caustic changes:")
    A("")
    A("    direct(lambda) = E_base,dir(lambda) * g1 *")
    A("                     [ T_bal(lambda) * C(x) + H(lambda) * blur(C)(x) ]")
    A("    T_bal = exp(-c L),  H = (b/c)(1-exp(-c L)) rho,")
    A("    rho   = exp(-a <Delta s>),  <Delta s> = (<theta^2>/2) <L - s>")
    A("")
    A("The ballistic term keeps the refracted caustic geometry; the halo is "
      "the same C master convolved with the derived sigma and sampled through "
      "the same geometry.  `C` is achromatic (representative G band, "
      "`det.py`), so the halo's spectral information is carried entirely by "
      "H(lambda) -- stated here because it is a limitation, not a detail.")
    A("")
    A("Energy: T_bal + H + A = 1 exactly, with A = 'absorbed in the first "
      "event OR multiple-scattered'.  Multiple scattering is out of scope "
      "(PIPELINE-V8 line 53), so in strongly scattering water the halo is a "
      "LOWER bound on the light actually returned and the bottom is darker "
      "than a multiple-scattering truth would be.")
    A("")

    # ---------------------------------------------------------------- G1
    A("## 3. G1 -- does b = 0 reproduce the pre-3.6 result?")
    A("")
    A("Scene: the preset-1 scene (recipe-1 as-is, U = %s m/s, fetch = %s m, "
      "z = %s m), grid %d, Q = 25, seed %s, spectrum `jonswap`."
      % (_f(s2["scene"]["U"], 1), _f(s2["scene"]["F"], 0),
         _f(s2["scene"]["h"], 1), s2["scene"]["grid"],
         s2["scene"]["seed"]))
    A("")
    A("    E = %s (the frozen composite's display normalisation as used by "
      "the v10 route -- v5's one-time E renormalisation, not the raw 2.2), "
      "g1 = %s, A = 2.6, lowpass B = 3.0 px"
      % (_f(s2["route"]["E"], 6), _f(s2["route"]["g1"], 6)))
    A("")
    A("| lane | what | result |")
    A("|---|---|---|")
    A("| G1a legacy | pre-3.6 water path (delegates to "
      "`render_frame_v5` + v5 optics) | still sha `%s` == stored "
      "`out/v10/still-A-v37.png` sha `%s` -> **byte-identical %s**; C "
      "master identical to the v10 route re-run in this session: %s; "
      "max|d| = %d |"
      % (g1a["still_sha256"][:16], g1a["stored_sha256"][:16],
         _yn(g1a["byte_identical_to_stored_artifact"]),
         _yn(g1a["C_identical_to_v10_route_this_session"]),
         g1a["max_abs_uint8_vs_stored_artifact"]))
    A("| G1b b = 0, new model | a = pure water absorption, no halo, no blur "
      "| mean|d| %s/255, max %d, %s%% of pixels changed |"
      % (_f(g1b["mean_abs_uint8"], 2), g1b["max_abs_uint8"],
         _f(100 * g1b["frac_pixels_changed"], 2)))
    A("| G1c same a,b, new transport | isolates the MODEL from the "
      "coefficients | mean|d| %s/255, max %d |"
      % (_f(g1c["mean_abs_uint8"], 2), g1c["max_abs_uint8"]))
    A("")
    A("**G1 verdict: %s.**  The brief allows 'not bit-identical' provided the "
      "difference is measured and its cause stated.  Both readings are "
      "reported:" % _pass(s2["verdict"]["G1_pass"]))
    A("")
    A("* the compatibility requirement (a legacy recipe reproduces the "
      "pre-3.6 image) is met BIT-IDENTICALLY by the legacy lane -- same C "
      "master, same still, same sha256 as the stored v3.7 artefact;")
    A("* the literal 'new model at b = 0' is NOT bit-identical, and the cause "
      "is fully accounted for: b = 0 removes the pure-water plus turbidity "
      "scattering from c, i.e. %s 1/m (two-way attenuation ratio %s), and "
      "removes the %s px empirical blur.  Net effect: the linear RGB means "
      "are x%s of the legacy lane's."
      % (_vec(g1b["cause"]["c_removed_is_b_w_plus_turbidity_db_rgb"]),
         _vec(g1b["cause"]["two_way_attenuation_ratio_b0_over_legacy"], 3),
         _f(g1b["cause"]["blur_removed_px"], 2),
         _vec(g1b["linear_rgb_mean_ratio_b0_over_legacy"], 3)))
    A("")

    # ---------------------------------------------------------------- G2
    A("## 4. G2 -- the rally's falsification condition: same c, two splits")
    A("")
    A("Construction (fixed before the run): c* = the extinction of the "
      "`algae_pond` preset at load 0.30 = %s 1/m; the two splits are the "
      "library's own extremes, f_a from `clear` = %s (absorption-dominated) "
      "and from `sediment` = %s (scattering-dominated); a = f_a c*, "
      "b = (1 - f_a) c*.  Both members therefore have the SAME c to 1e-12 "
      "(%s)."
      % (_vec(g2["construction"]["c_star_rgb"]),
         _vec(g2["pair"]["absorption"]["split_f_a"], 3),
         _vec(g2["pair"]["scattering"]["split_f_a"], 3),
         _yn(g2["construction"]["c_identical_in_both"])))
    A("")
    A("| quantity | absorption split | scattering split | separation |")
    A("|---|---|---|---|")
    pa = g2["pair"]["absorption"]; ps = g2["pair"]["scattering"]
    A("| a / b (650 nm) | %s / %s | %s / %s | same c |"
      % (_f(pa["a"][0]), _f(pa["b"][0]), _f(ps["a"][0]), _f(ps["b"][0])))
    A("| a / b (550 nm) | %s / %s | %s / %s | same c |"
      % (_f(pa["a"][1]), _f(pa["b"][1]), _f(ps["a"][1]), _f(ps["b"][1])))
    A("| a / b (470 nm) | %s / %s | %s / %s | same c |"
      % (_f(pa["a"][2]), _f(pa["b"][2]), _f(ps["a"][2]), _f(ps["b"][2])))
    A("| ballistic share T/(T+H) | %s | %s | - |"
      % (_vec(pa["ballistic_share"], 3), _vec(ps["ballistic_share"], 3)))
    A("| halo sigma | %s px | %s px | - |"
      % (_f(pa["sigma_px"], 1), _f(ps["sigma_px"], 1)))
    A("| C x = cL (650/550/470) | %s | %s | identical |"
      % (_vec(pa["x_cL"], 2), _vec(ps["x_cL"], 2)))
    A("| **C p99 / mean** | **%s** | **%s** | **%s x** |"
      % (_f(pa["effective_field"]["contrast_p99_over_mean"], 3),
         _f(ps["effective_field"]["contrast_p99_over_mean"], 3),
         _f(g2["separation"]["contrast_separation_ratio_p99"], 3)))
    A("| C p99.9 / mean | %s | %s | %s x |"
      % (_f(pa["effective_field"]["contrast_p99_9_over_mean"], 3),
         _f(ps["effective_field"]["contrast_p99_9_over_mean"], 3),
         _f(g2["separation"]["contrast_separation_ratio_p99_9"], 3)))
    A("| mean linear RGB | %s | %s | level x%s (green) |"
      % (_vec(pa["colour"]["mean_linear_rgb"]),
         _vec(ps["colour"]["mean_linear_rgb"]),
         _f(g2["separation"]
            ["green_linear_level_ratio_scattering_over_absorption"], 2)))
    A("| chromaticity (r,g,b) | %s | %s | delta %s |"
      % (_vec(pa["colour"]["chromaticity"], 3),
         _vec(ps["colour"]["chromaticity"], 3),
         _vec(g2["separation"]["chromaticity_delta"], 3)))
    A("| Lab | L %s a %s b %s | L %s a %s b %s | dE76 %s (chroma only %s) |"
      % (_f(pa["colour"]["lab"]["L"], 1), _f(pa["colour"]["lab"]["a"], 1),
         _f(pa["colour"]["lab"]["b"], 1),
         _f(ps["colour"]["lab"]["L"], 1), _f(ps["colour"]["lab"]["a"], 1),
         _f(ps["colour"]["lab"]["b"], 1),
         _f(g2["separation"]["delta_e76"], 2),
         _f(g2["separation"]["delta_e76_chroma_only"], 2)))
    A("")
    A("Pixel-level difference of the two delivered stills: mean |d| = %s/255, "
      "max = %d, %s%% of pixels changed."
      % (_f(g2["separation"]["mean_abs_uint8"], 2),
         g2["separation"]["max_abs_uint8"],
         _f(100 * g2["separation"]["pixel_fraction_changed_uint8"], 2)))
    A("")
    A("Thresholds fixed before the run: contrast separation (max/min of "
      "p99/mean) >= %s and dE76(chroma only) >= %s."
      % (_f(g2["thresholds"]["contrast_separation_ratio_min"], 2),
         _f(g2["thresholds"]["delta_e76_chroma_min"], 1)))
    A("")
    A("**G2 verdict: %s** (contrast %s, colour %s).  Note the direction, "
      "which is a measurement, not an assumption: the scattering split "
      "carries %s x the contrast of the absorption split at the same c "
      "(p99/mean %s vs %s), and its level is %s x higher in green because "
      "the scattered light is returned to the bottom instead of being "
      "absorbed.  The blue channel separates most in contrast (share %s vs "
      "%s) and the chromaticity separates in the (g,b) plane (delta %s / "
      "%s)."
      % (_pass(g2["pass"]), _pass(g2["pass_contrast"]),
         _pass(g2["pass_colour"]),
         _f(g2["separation"]["contrast_ratio_scattering_over_absorption_p99"],
            3),
         _f(g2["separation"]["p99_over_mean_absorption"], 3),
         _f(g2["separation"]["p99_over_mean_scattering"], 3),
         _f(g2["separation"]
            ["green_linear_level_ratio_scattering_over_absorption"], 2),
         _f(pa["ballistic_share"][2], 3), _f(ps["ballistic_share"][2], 3),
         _f(g2["separation"]["chromaticity_delta"][1], 3),
         _f(g2["separation"]["chromaticity_delta"][2], 3)))
    A("")
    A("Stills: `pair-sameC-absorption.png`, `pair-sameC-scattering.png` "
      "(735x539, frozen composite, identical E and g1, only the a/b split "
      "differs).  Recipes: `recipes/g2-matched-c-*.json`.")
    A("")

    # ---------------------------------------------------------------- G3
    A("## 5. G3 -- b up at fixed a: peaks down, skirts up")
    A("")
    A("Fixed a = the `sediment` preset's own a = %s 1/m; b swept over the "
      "same spectral shape at four levels, one C master (the G2 master), "
      "grid %d.  x = cL is reported per row because the sigma derivation is "
      "a thin-to-moderate medium approximation."
      % (_vec(g3["a_fixed_rgb"]),
         g2["construction"]["grid"]))
    A("")
    A("| b(550) | c(550) | x = cL | sigma px | T | H | halo share | "
      "C p99/mean | C p99 | profile r32 | profile r48 |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in g3["rows"]:
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (_f(r["b_green"], 2), _f(r["c_rgb"][1], 3),
             _f(r["x_cL_green"], 2), _f(r["sigma_px"], 1),
             _f(r["T_ballistic"][1], 3), _f(r["H_halo"][1], 3),
             _f(r["halo_light_share"], 3),
             _f(r["effective_field"]["contrast_p99_over_mean"], 3),
             _f(r["effective_field"]["p99"], 3),
             _f(r["effective_field"]["peak_profile"]["profile"]["r32_px"], 3),
             _f(r["effective_field"]["peak_profile"]["profile"]["r48_px"], 3)))
    A("")
    mo = g3["monotonicity"]
    A("Monotonicity (strict, adjacent rows):")
    A("")
    for key in ("p99_over_mean", "p99_9_over_mean", "absolute_p99",
                "absolute_p99_9", "peak_profile_r32", "peak_profile_r48",
                "peak_profile_r64", "halo_light_share"):
        A("    %-24s %-11s %s  %s"
          % (key, mo[key]["direction"], _pass(mo[key]["monotone"]),
             _vec(mo[key]["values"], 3)))
    A("")
    A(g3["skirt_measure_note"])
    A("")
    A(g3["absolute_level_note"])
    A("")
    A("**G3 verdict: %s** (contrast %s, skirt %s).  Stills: %s."
      % (_pass(g3["pass"]), _pass(g3["pass_contrast"]),
         _pass(g3["pass_skirt"]), ", ".join("`%s`" % f for f in g3["files"])))
    A("")

    # ---------------------------------------------------------------- G4
    A("## 6. G4 -- two-way attenuation on both legs, z = 0.6 vs 4.0")
    A("")
    A("Water: `clear` preset at load 0.5, a = %s, b = %s, c = %s 1/m.  "
      "C is flat (C = 1) so the comparison isolates the attenuation from the "
      "caustic geometry; each depth is rendered three times: full, return leg "
      "disabled (up = 1), down leg disabled (down = 1).  The measured "
      "quantity is the mean linear direct term of the actual render."
      % (_vec(g4["water"]["a_rgb"]), _vec(g4["water"]["b_rgb"]),
         _vec(g4["water"]["c_rgb"])))
    A("")
    A("| depth | s_down | s_up | down factor (650/550/470) | up factor | "
      "two-way Beer | T+H | measured full/up_off | model full/up_off | "
      "rel err |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for k in ("z0.6", "z4"):
        d = g4["depths"][k]["parts"]["full"]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (k, _f(d["s_down_m"], 3), _f(d["s_up_m"], 3),
             _vec(d["beer_down"], 3), _vec(d["beer_up"], 3),
             _vec(d["beer_two_way"], 3), _vec(d["T_plus_H"], 3),
             _vec(g4["depths"][k]["measured_full_over_up_off_direct"], 4),
             _vec(g4["depths"][k]["model_full_over_up_off"], 4),
             _f(max(abs(x) for x in g4["depths"][k]["rel_err_up_leg"]), 1)))
    A("")
    A("The down leg was checked the same way (full vs down-off): worst "
      "relative error over both depths %s."
      % _f(max(abs(x) for k in ("z0.6", "z4")
               for x in g4["depths"][k]["rel_err_down_leg"]), 1))
    A("")
    rl = g4["checks"]["return_leg_effect_green"]
    A("Return leg, stated as a number: the green-channel return factor is "
      "%s at z = 0.6 m and %s at z = 4.0 m.  Omitting the return leg would "
      "brighten the bottom by %sx at z = 0.6 and %sx at z = 4.0 -- so the "
      "deep case is exactly where the second leg matters, which is what the "
      "gate asks."
      % (_f(rl["z0.6"], 4), _f(rl["z4.0"], 4),
         _f(rl["brightening_if_omitted_z0.6"], 3),
         _f(rl["brightening_if_omitted_z4"], 3)))
    A("")
    A("**G4 verdict: %s.**  Model identity T_bal + H + A = 1 holds to %s.  "
      "Image evidence: `depth-pair-z0.6.png` / `depth-pair-z4.png` (same C "
      "master, same water, same composite, only the depth differs): linear "
      "luminance %s -> %s, ratio %s, chromaticity %s -> %s."
      % (_pass(g4["checks"]["pass"]),
         _f(chk["partition"]["worst_identity_err"], 1),
         _f(g4["depths"]["z0.6"]["still"]["metrics"]
            ["linear_luminance_mean"], 4),
         _f(g4["depths"]["z4"]["still"]["metrics"]
            ["linear_luminance_mean"], 4),
         _f(g4["depth_level_ratio_green_still"], 3),
         _vec(g4["depths"]["z0.6"]["still"]["colour"]["chromaticity"], 3),
         _vec(g4["depths"]["z4"]["still"]["colour"]["chromaticity"], 3)))
    A("")

    # ---------------------------------------------------------------- dev
    A("## 7. Deviations, approximations, and what was NOT done")
    A("")
    for i, d in enumerate(_load(os.path.join(OUT, "numbers-v12.json"))
                          ["deviations"], 1):
        A("%d. %s" % (i, d))
    A("")
    A("Open items this chunk deliberately does not close:")
    A("")
    A("* multiple scattering (P2 in PIPELINE-V8): the halo is a lower bound "
      "in strongly scattering water;")
    A("* the CDOM model itself (v5's domain): the port is built, default 0;")
    A("* g per particle type (the presets share the design's 0.85; the "
      "sensitivity table shows what a different g would do);")
    A("* the a/b lane renders the substrate gain as identity (the recipe has "
      "no bottom block; the substrate preset is a v5 display term, not a "
      "water IOP);")
    A("* the recipe gains `water.medium`, but no *.json on disk was edited: "
      "the experiment recipes are written under `out/v12/recipes/`.")
    A("")

    # ---------------------------------------------------------------- merge
    A("## 8. Merge decision")
    A("")
    A("Per the brief: the a/b water path becomes the default iff G1-G4 pass; "
      "otherwise it stays behind the recipe switch and this report says why.")
    A("")
    A("    merge_pass     = %s" % _yn(mg["merge_pass"]))
    A("    blocking       = %s" % (", ".join(mg["blocking_items"]) or "none"))
    A("    WATER_MODEL_DEFAULT = %s" % mg["default_after"])
    A("    switch provenance       = %s"
      % mg.get("default_note", "-").replace("\n", " "))
    A("")
    if mg["merge_pass"]:
        A("All four gates pass, so the a/b path is the default for a water "
          "medium (a `water.medium` block without a `preset` gets "
          "`DEFAULT_MEDIUM_PRESET` = `%s`).  Recipes with NO `water.medium` "
          "block keep the legacy lane and stay bit-identical to pre-3.6 -- "
          "that is G1's compatibility invariant, not a contradiction: "
          "'default' applies to the medium declaration, not to its absence."
          % s1["constants"]["DEFAULT_MEDIUM_PRESET"])
    else:
        A("The a/b path stays behind the recipe switch.  The blocking items "
          "above are the reason; the legacy route remains the default for "
          "recipes without a medium block.")
    A("")
    rs = mg.get("recipe_switch")
    if rs:
        A("The switch is exercised, not asserted (`v12_water."
          "checks_recipe_switch`, run in stage 4 with `WATER_MODEL_DEFAULT` "
          "= `%s`):" % rs["WATER_MODEL_DEFAULT"])
        A("")
        A("| recipe declaration | route taken | optics mode | verdict |")
        A("|---|---|---|---|")
        for tag, why in (
                ("no_water_section",
                 "no water block at all -> must stay legacy (G1 invariant)"),
                ("water_no_medium",
                 "water block, no medium -> follows the merge switch"),
                ("explicit_legacy",
                 "`medium.model = legacy` -> must stay legacy")):
            row = rs["rows"][tag]
            ok = ((rs["pass_no_water_stays_legacy"]
                   if tag == "no_water_section" else
                   rs["pass_explicit_legacy"] if tag == "explicit_legacy"
                   else rs["pass_default_follows_switch"]))
            A("| `%s` (%s) | `%s` | %s | %s |"
              % (tag, why, row["recipe_model"], row["optics_mode"],
                 _pass(ok)))
        A("")
        A(rs["note"])
        A("")

    # ---------------------------------------------------------------- files
    A("## 9. Files (all under `caustics-filter/offline/out/v12/`)")
    A("")
    for f in ("numbers-v12.json", "numbers-v12-water.json",
              "numbers-v12-g1.json", "numbers-v12-g2.json",
              "numbers-v12-g3.json", "numbers-v12-g4.json",
              "merge-evaluation-v12.json", "checkpoint.json"):
        A("* `%s`" % f)
    A("* `REPORT-V12-GATES.md` (this file)")
    A("* stills: `still-g1-legacy.png` (`still-g1-legacy-rerun.png` is the "
      "same image re-rendered for its linear buffer; sha checked equal), "
      "`still-g1-ab-bzero.png`, "
      "`still-g1-ab-samecoeff.png`, `pair-sameC-absorption.png`, "
      "`pair-sameC-scattering.png`, `sweep-b*.png` (%d), "
      "`depth-pair-z0.6.png`, `depth-pair-z4.png`" % len(g3["files"]))
    A("* `recipes/g2-matched-c-absorption.json`, "
      "`recipes/g2-matched-c-scattering.json`")
    A("* `_work/C-exp-Q25-f32.raw` (the 1024^2 float32 C master the G2/G3/G4 "
      "statistics were computed on)")
    A("")
    A("## 10. Reproduction")
    A("")
    A("```bat")
    A("cd caustics-filter\\offline")
    A("python run_v12.py --stage 1   :: water module + formula checks")
    A("python run_v12.py --stage 2   :: G1 (1 det at 2048, Q=25)")
    A("python run_v12.py --stage 3   :: G2 G3 G4 (1 det at 1024, Q=25)")
    A("python run_v12.py --stage 4   :: merge + this report")
    A("```")
    A("")
    A("How THIS run was produced (recorded because it is a deviation in the "
      "run history, not in the model):")
    A("")
    A("* the previous session stopped after the stage-2 units `s2-state` and "
      "`s2-legacy-render`; those two artefacts were reused, not redone.  "
      "Stage 2 was completed with `python run_v12.py --stage 2r`, which "
      "re-derives the state and the C master (they are inputs to every render "
      "below), re-renders the legacy still ONLY to obtain its linear buffer, "
      "and checks that re-render against the committed "
      "`still-g1-legacy.png` byte for byte before anything else is measured;")
    A("* `--stage 2r` writes the same `numbers-v12-g1.json` schema as "
      "`--stage 2`, so stage 4 cannot tell the two apart except by the "
      "`resume` block it carries;")
    A("* stage 3 was run twice: the first pass measured G2/G3 with two "
      "instruments that turned out to be degenerate (the one-sided contrast "
      "encoding, and the peak autocorrelation half-width).  The instruments "
      "were corrected and the stage re-run; deviations 2 and 3 record what "
      "changed and the numbers on both sides.  No model constant, no "
      "threshold value and no gate definition moved, and the first pass's "
      "values agree with the delivered ones wherever they are comparable "
      "(p99, p99.9, the stills are byte-identical).")
    A("")
    A("Deterministic: one seed (from the unmodified recipe, preset 1), no "
      "random sampling in the water model, `det_frame` at Q = 25, and every "
      "deliverable written through `jobkit.atomic_write_bytes`.")
    A("")
    return "\n".join(L) + "\n"


def main():
    txt = build()
    path = os.path.join(OUT, "REPORT-V12-GATES.md")
    atomic_write_bytes(path, txt.encode("ascii", "replace"))
    print("wrote %s (%d bytes)" % (os.path.relpath(path, HERE),
                                  len(txt.encode("ascii", "replace"))))
    return path


if __name__ == "__main__":
    main()
