"""run_v14_spectrum -- version 3.9 boundary redo, THREE-ARM measurement.

Brief: BRIEF-V14-39REDO.md.  Design: PIPELINE-V12.md.  Base: the v13 tree
(Elfouhaily-form spectrum + labelled provisional prior below the ECKV fit
crossing, `spectrum_source` switch, frozen composite, WATER_MODEL_DEFAULT=ab
untouched).  v5..v13 files are IMPORTED, never edited.

The one quantity that moves between arms is the level `alpha_m` of the
Elfouhaily high-frequency part B_h (V10.elfouhaily_bh).  B_h is exactly linear
in alpha_m, so a cell's capillary slope variance is

    sigma2_h(cell) = alpha_m * I(cell),
    I(cell) = int_{kmin}^{k_vis} k^2 * B_h(k; alpha_m=1) / k^3 dk

and every arm is one rule for alpha_m_used:

  RAW        (v13 as shipped, the control)
             u* >= c_m/e : alpha_m_used = alpha_m_raw
             u* <  c_m/e : alpha_m_used = max(alpha_m_raw, alpha_m_floor)
             list of invented inputs: the floor, acting only below the crossing.

  B          (prior zone extended -- GPT/Grok rank 1)
             alpha_m_used = max(alpha_m_raw, alpha_m_floor) at EVERY cell
             with alpha_m_floor = sigma2_h_target / I(cell), i.e. the applied
             level is max(raw, floor): continuous in u* by construction (the
             two branches meet where raw = floor).  Provenance `low_wind_prior`
             wherever the top-up acts (floor > raw), including U=3 when it acts.

  D          (evaluation arm -- Gemini rank 1)
             changes NO numbers: it is RAW's spectrum, evaluated at U<=3 with
             Gemini's ABSOLUTE criteria instead of the band-share gates.
             Provenance `low_wind_absolute_gate` (a gate label, not a spectrum
             edit).  Rationale for D == RAW numerically: the brief says D
             "changes NO numbers" and the U=6/10 rows are "arm D same numbers".
             The literal alternative ("no prior anywhere") is measured too and
             reported as D_STRICT, a labelled diagnostic -- see the report's
             Deviations section.

The ECKV fit crossing (v13's boundary) is u* = c_m/e = 0.085048 m/s.  The
Arm-B hand-off is the u* where raw sigma2_h = target; it is solved per fetch
in `handoff_scan` and reported, and the scan proves is applies continuously.
"""
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10
import v11_spectrum as V11
import v13_spectrum as V13

G = P8.G
C_M = V11.C_M
K_M = V11.K_M
USTAR_FIT_MIN = V13.USTAR_FIT_MIN

ARMS = ("RAW", "B", "D")
D_STRICT = "D_STRICT"

#: the design's provisional prior target and the mandatory sweep points.
SIGMA2_H_TARGET_DEFAULT = V13.SIGMA2_H_TARGET_DEFAULT
SIGMA2_H_TARGETS = V13.SIGMA2_H_TARGETS

#: gate constants, PIPELINE-V12 (values verbatim; v13's V3 used 0.010/0.007).
B_ENV_MAX = 0.012          # max B(k) over k in [10, 1500]
BETA_HAT = 0.008           # max B_l(k) over the 50-200 band
GRAVITY_SHARE_MAX = 0.50   # C1
CAP_SHARE_MIN = 0.15       # C1b clause 1
SIGMA2_H_MIN = 0.002       # C1b clause 2
BAND = (50.0, 200.0)
B_K_RANGE = (10.0, 1500.0)
HS_TOL = 0.35              # |Hs_spectrum / Hs_law - 1|
TOTALS_TOL = (0.5, 2.0)    # PIPELINE-V12 line 26

PROVENANCE = {
    "RAW": "v13_shipped_prior_below_fit_crossing",
    "B": "low_wind_prior",
    "D": "low_wind_absolute_gate",
    D_STRICT: "raw_strict_no_prior",
}


def arm_definition(arm):
    return {
        "RAW": ("v13 as shipped: alpha_m raw inside the ECKV fit domain, "
                "max(raw, floor) below it (prior only where u* < c_m/e)"),
        "B": ("prior zone extended: alpha_m_used = max(raw, floor) at every "
              "cell, floor solving sigma2_h,total = target; continuous in u* "
              "by construction; provenance low_wind_prior where it acts"),
        "D": ("evaluation arm, NO spectrum change (identical to RAW): at U<=3 "
              "the band-share gates are replaced by the absolute criteria "
              "(B<=0.012 over 10-1500, B_l<=0.008 over 50-200, totals in "
              "[Cox-Munk clean, Hu], C1 smoothness)"),
        D_STRICT: ("DIAGNOSTIC, not one of the three arms: the literal 'raw "
                   "spectrum untouched' reading -- alpha_m = max(raw, 0) "
                   "everywhere (the v3.8 hard zero below the crossing), no "
                   "floor anywhere"),
    }[arm]


def provenance_label(arm, top_up_acted, in_domain):
    """The provenance string this cell carries under this arm."""
    if arm == "B":
        if top_up_acted:
            return PROVENANCE["B"]
        return ("eckv_fit_domain_raw" if in_domain
                else "below_fit_domain_no_topup")
    if arm == "D":
        return PROVENANCE["D"] if not in_domain else "eckv_fit_domain_raw"
    if arm == D_STRICT:
        return (PROVENANCE[D_STRICT] if not in_domain
                else "eckv_fit_domain_raw")
    if arm == "RAW":
        return (PROVENANCE["RAW"] if not in_domain else "eckv_fit_domain_raw")
    raise ValueError("unknown arm %r" % (arm,))


# --------------------------------------------------------------------------
# 1. assembly
# --------------------------------------------------------------------------
def assemble(U, F, h, arm="RAW", hs_override=None, k_cut=None, gamma=None,
             source="elfouhaily-form",
             sigma2_h_target=SIGMA2_H_TARGET_DEFAULT):
    """V11 (v3.8) spectrum + this arm's alpha_m resolution.

    Returns the v11/v13 ctx with the same key set v13 emits, plus
    `v14_arm` and the arm-specific `alpha_m_resolution` block.  source
    "jonswap" is returned BYTE-FOR-BYTE untouched (the legacy lane).
    """
    if arm not in ARMS and arm != D_STRICT:
        raise ValueError("unknown arm %r" % (arm,))
    ctx = V11.assemble(U, F, h, hs_override=hs_override, k_cut=k_cut,
                       gamma=gamma, source=source)
    if source != "elfouhaily-form":
        return ctx

    target = float(sigma2_h_target)
    raw = float(ctx["alpha_m"]["alpha_m_raw"])
    unit = V13.alpha_m_unit_integral(ctx)
    I = float(unit["integral"])
    floor = (target / I) if (target > 0.0 and I > 0.0) else 0.0
    ustar = float(ctx["ustar"]["ustar_ms"])
    in_dom = bool(ustar >= USTAR_FIT_MIN)

    if arm in ("RAW", "D"):
        if in_dom:
            used, domain, prior, top_up = raw, V13.ALPHA_M_DOMAIN_FIT, None, False
        else:
            used = max(raw, floor)
            domain, prior = V13.ALPHA_M_DOMAIN_PRIOR, V13.low_wind_prior_label(
                target)
            top_up = bool(floor > raw)
    elif arm == "B":
        used = max(raw, floor)
        top_up = bool(floor > raw)
        if top_up:
            domain, prior = V13.ALPHA_M_DOMAIN_PRIOR, V13.low_wind_prior_label(
                target)
        else:
            domain, prior = (V13.ALPHA_M_DOMAIN_FIT if in_dom
                             else V13.ALPHA_M_DOMAIN_PRIOR_OFF), None
    else:                                   # D_STRICT (diagnostic)
        used = max(raw, 0.0)
        top_up = False
        domain, prior = (V13.ALPHA_M_DOMAIN_FIT if in_dom
                         else V13.ALPHA_M_DOMAIN_PRIOR_OFF), None

    res = {
        "arm": arm,
        "u_ms": float(ctx["laws"]["U10_ms"]),
        "ustar_ms": ustar,
        "c_m_ms": float(C_M),
        "ustar_fit_min_ms": USTAR_FIT_MIN,
        "u_over_cm": float(ustar / C_M),
        "u_over_fit_min": float(ustar / USTAR_FIT_MIN),
        "in_eckv_fit_domain": in_dom,
        "alpha_m_raw": raw,
        "alpha_m_floor": float(floor),
        "alpha_m_used": float(used),
        "alpha_m_domain": domain,
        "prior_active": bool(prior is not None),
        "low_wind_prior": prior,
        "top_up_acted": bool(top_up),
        "prior_contribution_alpha_m": float(used - raw),
        "sigma2_h_target": target,
        "sigma2_h_unit_integral": I,
        "sigma2_h_total_closure": float(used * I),
        "sigma2_h_raw_would_be": float(raw * I),
        "prior_contribution_sigma2_h": float(used * I - raw * I),
        "provenance": provenance_label(arm, top_up, in_dom),
        "rule": arm_definition(arm),
    }
    ctx["version"] = "3.9+v14"
    ctx["v14_arm"] = arm
    ctx["alpha_m_value"] = float(used)
    ctx["alpha_m_used"] = float(used)
    ctx["alpha_m_raw"] = raw
    ctx["alpha_m_floor"] = float(floor)
    ctx["alpha_m_domain"] = domain
    ctx["alpha_m_resolution"] = res
    ctx["low_wind_prior"] = prior
    ctx["sigma2_h_target"] = target
    ctx["sigma2_h_total"] = float(used * I)
    ctx["provenance"] = res["provenance"]
    ctx["label"] = (ctx["label"] + " -- v3.9 boundary redo, arm %s" % arm)
    # H_s is a PREDICTION of the assembled spectrum (PIPELINE-V10 section 4).
    kcut = float(ctx["k_cut_rad_per_m"])
    k = V11._k_grid(kcut, k_hi=float(P8.viscous_cutoff(float(h))))
    s = V11.s_parts(k, ctx)
    kmin = V13.cell_kmin(ctx["laws"]["lambda_p_m"])[0]
    m = k >= kmin
    ctx["Hs_m"] = float(4.0 * math.sqrt(
        max(float(np.trapz(s["S"][m], k[m])), 0.0)))
    ctx["Hs_pred_rel_law"] = float(ctx["Hs_m"] / ctx["Hs_law_m"] - 1.0)
    return ctx


def s_parts(k, ctx):
    return V11.s_parts(k, ctx)


def cell_geometry(ctx):
    kmin, pm = V13.cell_kmin(ctx["laws"]["lambda_p_m"])
    return {"kmin_rad_per_m": float(kmin), "patch_m": float(pm),
            "k_viscous_rad_per_m": float(P8.viscous_cutoff(
                float(ctx["depth_m"])))}


# --------------------------------------------------------------------------
# 2. the numbers per cell
# --------------------------------------------------------------------------
def moments(ctx, k_cut=None):
    """sigma2 / m0 per part, on the cell's resolved band [kmin, k_vis]."""
    geo = cell_geometry(ctx)
    kmin = geo["kmin_rad_per_m"]
    khi = geo["k_viscous_rad_per_m"]
    k = V11._k_grid(khi, k_hi=khi)
    s = V11.s_parts(k, ctx)
    m = k >= kmin
    kk = k[m]
    out = {}
    for tag, arr in (("gravity", s["S_l"][m]), ("capillary", s["S_h"][m]),
                     ("total", s["S"][m])):
        out[tag] = {"sigma2": float(np.trapz(kk ** 2 * arr, kk)),
                    "m0": float(np.trapz(arr, kk))}
    # capillary share of the slope variance above 370 rad/m
    mh = m & (k >= 370.0)
    cap_all = out["capillary"]["sigma2"]
    cap_hi = float(np.trapz(k[mh] ** 2 * s["S_h"][mh], k[mh]))
    tot_hi = float(np.trapz(k[mh] ** 2 * s["S"][mh], k[mh]))
    out["frac_sigma2_above_370_total"] = (
        float(tot_hi / out["total"]["sigma2"])
        if out["total"]["sigma2"] else None)
    out["capillary_frac_above_370"] = (float(cap_hi / cap_all)
                                       if cap_all else None)
    out["band_rad_per_m"] = [kmin, khi]
    return out


def b_envelope(ctx, k_lo=B_K_RANGE[0], k_hi=B_K_RANGE[1], n=4000):
    """Max B(k) over [10,1500] and max B_l over the 50-200 band."""
    k = np.logspace(math.log10(float(k_lo)), math.log10(float(k_hi)), int(n))
    s = V11.s_parts(k, ctx)
    B = s["B_l"] + s["B_h"]
    i = int(np.argmax(B))
    kb = np.logspace(math.log10(BAND[0]), math.log10(BAND[1]), 3000)
    sb = V11.s_parts(kb, ctx)
    j = int(np.argmax(sb["B_l"]))
    Bt = sb["B_l"] + sb["B_h"]
    return {
        "k_range_rad_per_m": [float(k_lo), float(k_hi)],
        "B_max": float(B[i]), "k_at_B_max": float(k[i]),
        "B_max_where": "total B = B_l + B_h",
        "B_l_max": float(sb["B_l"][j]), "k_at_B_l_max": float(kb[j]),
        "B_total_max_50_200": float(np.max(Bt)),
        "k_at_B_total_max_50_200": float(kb[int(np.argmax(Bt))]),
        "B_env_max": B_ENV_MAX, "beta_hat": BETA_HAT,
        "pass_B_env": bool(B[i] <= B_ENV_MAX),
        "pass_B_l": bool(sb["B_l"][j] <= BETA_HAT),
        "pass": bool(B[i] <= B_ENV_MAX and sb["B_l"][j] <= BETA_HAT),
        "note": ("B_env_max applies to the TOTAL curvature everywhere in "
                 "[10,1500] (Gemini's absolute criterion, PIPELINE-V12 lines "
                 "15/26); beta_hat applies to B_l in the 50-200 band.  The "
                 "capillary peak sits at k ~ K_M = %.1f rad/m, i.e. outside "
                 "the 50-200 band, so a band-only max misses it."
                 % K_M),
    }


def envelope(U, F, h, kmin, k_cut):
    """The [Cox-Munk clean, Hu] bracket, matched definitions (v13's columns).

    "The Hu law" has no formula anywhere in this workspace or in either design
    document; as in v13's reporting the second column is the pipeline's own
    fetch-limited MSS integral and the naming gap is recorded, not hidden.
    """
    cm = V10.cox_munk_reference(float(U))
    hu = V13.fetch_limited_mss_reference(float(U), float(F), float(h),
                                         kmin=kmin, k_hi=k_cut)
    hu_v = float(hu["sigma2_full"]["jonswap_gamma_3p3"])
    return {"cox_munk_clean_rad2": float(cm["total_clean_rad2"]),
            "cox_munk": cm,
            "hu_named_column_rad2": hu_v,
            "hu_naming_status": hu["naming_status"],
            "hu_reference": hu,
            "k_hi_rad_per_m": float(k_cut),
            "kmin_rad_per_m": float(kmin),
            "definitions_matched": ("both columns are slope variances in "
                                    "rad^2; the model column is the assembled "
                                    "S(k) integrated over [kmin, k_cut]; the "
                                    "Hu-named column is the in-repo "
                                    "fetch-limited MSS integrated to the same "
                                    "k_cut")}


def envelope_check(total, env):
    """Three readings of the totals gate -- all reported, none hidden."""
    lo = float(env["cox_munk_clean_rad2"])
    hi = float(env["hu_named_column_rad2"])
    if hi < lo:
        lo, hi = hi, lo
    t = float(total)
    return {
        "total_sigma2_rad2": t,
        "lower_cox_munk_clean": lo, "upper_hu_named": hi,
        "ratio_to_lower": float(t / lo) if lo else None,
        "ratio_to_upper": float(t / hi) if hi else None,
        "strict_inside": bool(lo <= t <= hi),
        "strict_below_lower": bool(t < lo),
        "tolerance_band": [TOTALS_TOL[0] * lo, TOTALS_TOL[1] * hi],
        "tol_inside": bool(TOTALS_TOL[0] * lo <= t <= TOTALS_TOL[1] * hi),
        "note": ("strict = the interval as written in the design ([Cox-Munk "
                 "clean, Hu]); tol = PIPELINE-V12 line 26's explicit "
                 "[0.5x, 2.0x] band.  The gate used for the three-state "
                 "verdict is `tol_inside`; `strict_inside` is reported "
                 "beside it for every cell."),
    }


def hs_node(ctx):
    hs = float(ctx["Hs_m"])
    law = float(ctx["Hs_law_m"])
    d = float(hs / law - 1.0) if law else None
    return {"Hs_spectrum_m": hs, "Hs_fetch_law_m": law,
            "delta_rel": d, "abs_delta_rel": (None if d is None else abs(d)),
            "tolerance": HS_TOL,
            "pass": (None if d is None else bool(abs(d) <= HS_TOL)),
            "note": ("H_s is always SPECTRUM-INTEGRATED here (4 sqrt(int S "
                     "dk) over the cell's resolved band); it is never "
                     "rescaled to the law")}


def shares(ctx, integrals):
    o = integrals.get("overlap_50_200") or {}
    g = o.get("gravity_share_of_sigma2")
    return {"gravity_share_50_200": g,
            "capillary_share_50_200": (None if g is None else float(1.0 - g)),
            "B_l_sigma2_50_200": o.get("B_l_sigma2"),
            "B_h_sigma2_50_200": o.get("B_h_sigma2")}


def gate_c1(sh):
    g = sh["gravity_share_50_200"]
    return {"gate": "C1 -- gravity share of 50-200 slope variance <= 0.50",
            "value": g, "limit": GRAVITY_SHARE_MAX,
            "pass": (None if g is None else bool(g <= GRAVITY_SHARE_MAX))}


def gate_c1b(sh, mom):
    cap = sh["capillary_share_50_200"]
    s2h = float(mom["capillary"]["sigma2"])
    ok1 = None if cap is None else bool(cap >= CAP_SHARE_MIN)
    ok2 = bool(s2h >= SIGMA2_H_MIN)
    return {"gate": ("C1b -- capillary share 50-200 >= 0.15 AND "
                     "sigma2_h,total >= 0.002 (PIPELINE-V12 line 25; v13 "
                     "used OR -- the conjunction is this job's rule)"),
            "capillary_share_50_200": cap, "share_threshold": CAP_SHARE_MIN,
            "share_pass": ok1,
            "sigma2_h_total": s2h, "sigma2_threshold": SIGMA2_H_MIN,
            "sigma2_pass": ok2,
            "pass": (None if ok1 is None else bool(ok1 and ok2)),
            "changed_from_v13": ("v13's c1b() returned ok1 OR ok2; PIPELINE-V12 "
                                 "states the conjunction, so both clauses are "
                                 "required here")}


def smoothness(ctx, n=20000):
    """C1-in-k: B(k) continuous, dlnB/dlnk without steps.

    Two numbers per resolution, because a per-step bound alone cannot tell a
    STEP from a steep-but-smooth flank:
      * max |dB/B| in one grid step at n and 2n points on [10, 1500];
      * max |d(dlnB/dlnk)| at n and 2n, and the ratio fine/coarse.
    A smooth function's step change is proportional to the step, so it halves
    when the grid doubles (ratio ~0.5); a splice or a step does not (ratio ~1).
    The literal "<= 1% per step" line of the brief is reported as its own
    boolean, at the resolution it is met, with the location and the local
    log-slope of the largest step.
    """
    klo, khi = 10.0, 1500.0
    out = {}
    for tag, nn in (("coarse", int(n)), ("fine", int(n) * 2)):
        k = np.logspace(math.log10(klo), math.log10(khi), nn)
        s = V11.s_parts(k, ctx)
        B = s["B_l"] + s["B_h"]
        rel = np.abs(np.diff(B)) / np.maximum(B[:-1], 1e-300)
        i = int(np.argmax(rel))
        lk = np.log(k)
        slope = np.diff(np.log(np.maximum(B, 1e-300))) / np.diff(lk)
        dslope = np.abs(np.diff(slope))
        j = int(np.argmax(dslope))
        out[tag] = {
            "n": nn, "log_step": float(lk[1] - lk[0]),
            "max_rel_step": float(rel[i]), "k_at_max_rel_step": float(k[i]),
            "local_log_slope_at_max": float(slope[i]),
            "max_dslope": float(dslope[j]), "k_at_max_dslope": float(k[j + 1]),
            "frac_steps_above_1pct": float(np.mean(rel > 0.01)),
            "B_min": float(np.min(B)), "B_nonneg": bool(np.min(B) >= 0.0),
        }
    ratio_s = (float(out["fine"]["max_rel_step"]
                     / out["coarse"]["max_rel_step"])
               if out["coarse"]["max_rel_step"] > 0 else None)
    ratio_d = (float(out["fine"]["max_dslope"] / out["coarse"]["max_dslope"])
               if out["coarse"]["max_dslope"] > 0 else None)
    return {
        "gate": ("C1 smoothness -- B(k) continuous and dlnB/dlnk without "
                 "steps on [10,1500], tested by grid refinement"),
        "coarse": out["coarse"], "fine": out["fine"],
        "ratio_max_rel_step_fine_over_coarse": ratio_s,
        "ratio_max_dslope_fine_over_coarse": ratio_d,
        "refinement_bound": 0.65,
        "no_slope_discontinuity": bool(ratio_d is not None
                                       and ratio_d <= 0.65),
        "literal_1pct_met": bool(out["coarse"]["max_rel_step"] <= 0.01),
        "literal_1pct_violation_at": (
            None if out["coarse"]["max_rel_step"] <= 0.01 else {
                "k_rad_per_m": out["coarse"]["k_at_max_rel_step"],
                "max_rel_step": out["coarse"]["max_rel_step"],
                "local_log_slope": out["coarse"]["local_log_slope_at_max"],
                "explanation": ("the step is %0.2f%% of B in a %0.3f%% step of "
                                "ln k because |dlnB/dlnk| = %.1f there (the "
                                "steep flank below the gravity peak); the "
                                "step halves with the grid, so it is steepness "
                                "and not a discontinuity"
                                % (100.0 * out["coarse"]["max_rel_step"],
                                   100.0 * out["coarse"]["log_step"],
                                   -out["coarse"]["local_log_slope_at_max"]))}),
        "pass": bool(out["coarse"]["B_nonneg"] and out["fine"]["B_nonneg"]
                     and ratio_d is not None and ratio_d <= 0.65),
        "note": ("the absolute per-step change of B is resolution-dependent "
                 "for any steep smooth spectrum; the refinement ratio is the "
                 "discriminating test and both numbers are printed"),
    }


def positivity(ctx, n=20000, k_lo=None, k_hi=None):
    """A1: B >= 0, S >= 0, everything finite on a dense grid."""
    geo = cell_geometry(ctx)
    lo = float(geo["kmin_rad_per_m"] if k_lo is None else k_lo)
    hi = float(geo["k_viscous_rad_per_m"] if k_hi is None else k_hi)
    k = np.logspace(math.log10(lo), math.log10(hi), int(n))
    s = V11.s_parts(k, ctx)
    B = s["B_l"] + s["B_h"]
    bad = {
        "B_min": float(np.min(B)), "S_min": float(np.min(s["S"])),
        "B_l_min": float(np.min(s["B_l"])),
        "B_h_min": float(np.min(s["B_h"])),
        "all_finite": bool(np.all(np.isfinite(B))
                           and np.all(np.isfinite(s["S"]))),
    }
    bad["pass"] = bool(bad["B_min"] >= 0.0 and bad["S_min"] >= 0.0
                       and bad["all_finite"]
                       and bad["B_l_min"] >= 0.0 and bad["B_h_min"] >= 0.0)
    bad["band_rad_per_m"] = [lo, hi]
    bad["n"] = int(n)
    return bad


def independent_integrals(ctx, k_cut=None):
    """A4: independent grids for the same band integrals.

    The reference grid is v11's own 48k log grid; the checks are a finer log
    grid and a 150k *geometric* grid on the band (a different construction).
    A 12k log grid is reported as a coarsening sensitivity and is NOT part of
    the pass condition -- 12000 points over four decades is simply too coarse
    for the 50-200 band and says nothing about the spectrum.
    """
    geo = cell_geometry(ctx)
    kmin = geo["kmin_rad_per_m"]
    kc = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])

    def blocks(k):
        s = V11.s_parts(k, ctx)
        m = (k >= kmin) & (k <= kc)
        kk = k[m]
        tot = float(np.trapz(kk ** 2 * s["S"][m], kk))
        mo = m & (k >= BAND[0]) & (k <= BAND[1])
        ko = k[mo]
        g2 = float(np.trapz(ko ** 2 * s["B_l"][mo] / ko ** 3, ko))
        h2 = float(np.trapz(ko ** 2 * s["B_h"][mo] / ko ** 3, ko))
        return {"sigma2_total": tot, "sigma2_gravity_50_200": g2,
                "sigma2_capillary_50_200": h2, "N": int(k.size)}

    out = {"log_48k": blocks(V11._k_grid(kc)),
           "log_96k": blocks(V11._k_grid(kc, n=96000)),
           "geom_150k": blocks(np.geomspace(kmin, kc, 150001)),
           "log_12k": blocks(V11._k_grid(kc, n=12000))}
    ref = out["log_48k"]
    rel = {}
    for tag in ("log_96k", "geom_150k", "log_12k"):
        rel[tag] = {
            "total": float(out[tag]["sigma2_total"]
                           / ref["sigma2_total"] - 1.0),
            "g_50_200": float(out[tag]["sigma2_gravity_50_200"]
                              / max(ref["sigma2_gravity_50_200"], 1e-300) - 1.0),
            "h_50_200": float(out[tag]["sigma2_capillary_50_200"]
                              / max(ref["sigma2_capillary_50_200"], 1e-300)
                              - 1.0)}
    worst = max(abs(rel[t][k]) for t in ("log_96k", "geom_150k")
                for k in ("total", "g_50_200", "h_50_200"))
    return {"grids": out, "rel_vs_log_48k": rel, "worst_abs_rel": float(worst),
            "checks": ["log_96k", "geom_150k"],
            "reported_not_gated": ["log_12k"],
            "tolerance": 1.0e-3, "pass": bool(worst <= 1.0e-3),
            "band_rad_per_m": [kmin, kc],
            "note": ("three grids integrate the same assembled parts; the "
                     "acceptance applies to the finer log grid and the "
                     "independent geometric grid only")}


def band_partition(k, s, kmin, k_hi, k_p=None):
    """A correct band partition of [kmin, k_hi] (v14's replacement).

    v13 used `V11.band_table`, whose edge list is [kmin, k_p, 10, 50, ...]:
    when k_p > 10 rad/m (which happens at short fetch) that list is NOT sorted,
    so the rows overlap and the table sums to MORE than the cell total.  The
    v14 partition sorts the design's fixed edges and drops any that fall
    outside [kmin, k_hi], so the rows are disjoint and exhaustive.
    """
    edges = sorted({float(kmin), float(k_hi)}
                   | {float(e) for e in V11.BAND_EDGES_V11
                      if float(kmin) < float(e) < float(k_hi)})
    rows = {}
    tot = {"m0": 0.0, "sigma2": 0.0}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (k >= lo) & (k < hi)
        kk = k[m]
        r = {"k": [float(lo), float(hi)], "empty": bool(not m.any())}
        if m.any():
            r["gravity"] = {"m0": float(np.trapz(s["B_l"][m] / kk ** 3, kk)),
                            "sigma2": float(np.trapz(kk ** 2 * s["B_l"][m]
                                                     / kk ** 3, kk))}
            r["capillary"] = {"m0": float(np.trapz(s["B_h"][m] / kk ** 3, kk)),
                              "sigma2": float(np.trapz(kk ** 2 * s["B_h"][m]
                                                       / kk ** 3, kk))}
            r["total"] = {kk2: r["gravity"][kk2] + r["capillary"][kk2]
                          for kk2 in ("m0", "sigma2")}
            r["B_l_max"] = float(np.max(s["B_l"][m]))
            r["B_h_max"] = float(np.max(s["B_h"][m]))
        else:
            r["gravity"] = r["capillary"] = {"m0": 0.0, "sigma2": 0.0}
            r["total"] = {"m0": 0.0, "sigma2": 0.0}
            r["B_l_max"] = r["B_h_max"] = None
        rows["%g_%g" % (lo, hi)] = r
        tot["m0"] += r["total"]["m0"]
        tot["sigma2"] += r["total"]["sigma2"]
    # the total on the same grid, and v13's table for the defect record
    m_all = (k >= float(kmin)) & (k <= float(k_hi))
    kk = k[m_all]
    S = s["S"][m_all]
    tot_ref = {"m0": float(np.trapz(S, kk)),
               "sigma2": float(np.trapz(kk ** 2 * S, kk))}
    v13tab = V11.band_table(k, s, float(k_p if k_p is not None else k_hi),
                            float(kmin), float(k_hi))
    v13_sum = sum(r["total"]["sigma2"] for r in v13tab["rows"].values())
    for name, r in rows.items():
        r["sigma2_share"] = (float(r["total"]["sigma2"]
                                   / tot_ref["sigma2"])
                             if tot_ref["sigma2"] else 0.0)
        r["m0_share"] = (float(r["total"]["m0"] / tot_ref["m0"])
                         if tot_ref["m0"] else 0.0)
        g = r["total"]["sigma2"]
        r["gravity_share"] = (float(r["gravity"]["sigma2"] / g) if g else None)
    return {
        "rows": rows, "edges_rad_per_m": edges,
        "sum_rows": tot, "reference_total": tot_ref,
        "rel_rows_over_total_sigma2": float(tot["sigma2"]
                                            / tot_ref["sigma2"] - 1.0),
        "rel_rows_over_total_m0": float(tot["m0"] / tot_ref["m0"] - 1.0),
        "pass": bool(abs(tot["sigma2"] / tot_ref["sigma2"] - 1.0) <= 1e-3
                     and abs(tot["m0"] / tot_ref["m0"] - 1.0) <= 1e-3),
        "v13_band_table_sum_sigma2": float(v13_sum),
        "v13_band_table_rel_over_total": float(v13_sum
                                               / tot_ref["sigma2"] - 1.0),
        "note": ("disjoint, sorted, exhaustive; the residual is trapezoid "
                 "bookkeeping at the internal edges, not a gap"),
    }


def band_reproduction(ctx, k_cut=None):
    """A4 (second clause): the band partition vs the cell total."""
    geo = cell_geometry(ctx)
    kmin = geo["kmin_rad_per_m"]
    kc = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    k = V11._k_grid(kc)
    s = V11.s_parts(k, ctx)
    bp = band_partition(k, s, kmin, kc, float(ctx["k_p_rad_per_m"]))
    ib = V11.integral_block(k, s, kmin, kc, float(ctx["k_p_rad_per_m"]))
    return {
        "partition": bp,
        "total_sigma2_integral_block": float(ib["total"]["sigma2"]),
        "rel_partition_over_integral_block_sigma2": float(
            bp["sum_rows"]["sigma2"] / ib["total"]["sigma2"] - 1.0)
        if ib["total"]["sigma2"] else None,
        "pass": bool(bp["pass"] and abs(
            bp["sum_rows"]["sigma2"] / ib["total"]["sigma2"] - 1.0) <= 1e-3),
        "v13_band_table_defect": {
            "rel_over_total": bp["v13_band_table_rel_over_total"],
            "k_p_rad_per_m": float(ctx["k_p_rad_per_m"]),
            "fired": bool(abs(bp["v13_band_table_rel_over_total"]) > 1e-6),
            "cause": ("V11.band_table's edge list is [kmin, k_p, 10, 50, ...]; "
                      "when k_p > 10 rad/m it is unsorted, so kp_10 is empty "
                      "and 10_50 overlaps below_kp, double counting "
                      "[10, k_p]"),
            "note": ("found by this job's A4 bookkeeping test; it affects "
                     "v13's per-band SHARE columns, not the mask-based "
                     "gravity/capillary shares used by C1/C1b")},
        "note": ("the partition is the v14 replacement; the v13 table is kept "
                 "only to measure its own defect"),
    }


# --------------------------------------------------------------------------
# 3. the arm-B hand-off in u*
# --------------------------------------------------------------------------
def _rule_scan(F, h, target, Us):
    """max(raw, floor) with the BAND HELD FIXED -- the rule alone.

    The fixed band removes the patch-quantisation effect (the resolved band's
    lower edge moves in discrete doublings), so what is left is the rule.
    """
    ctx = assemble(3.0, float(F), float(h), arm="D")
    I = float(ctx["alpha_m_resolution"]["sigma2_h_unit_integral"])
    floor = float(target) / I
    rows = []
    for U in Us:
        ustar = float(V10.ustar_from_u10(float(U))["ustar_ms"])
        raw = float(V10.alpha_m_elfouhaily(ustar)["alpha_m_raw"])
        rows.append((float(U), ustar, raw, max(raw, floor),
                     max(raw, floor) * I))
    return {"I_fixed": I, "floor_fixed": floor,
            "band_from_cell": {"U_ref_ms": 3.0, "fetch_m": float(F),
                               "depth_m": float(h)},
            "rows": rows}


def handoff_scan(F, h, target=SIGMA2_H_TARGET_DEFAULT, n=120,
                 U_lo=0.8, U_hi=12.0):
    """Where Arm B's top-up starts/stops acting, and how smooth the join is.

    Two scans:
      * `rule`      -- max(raw, floor) with the band held fixed: this isolates
                       the hand-off RULE.  Continuity is tested by refinement
                       (n vs 2n): a step keeps its size, a steep smooth ramp
                       halves with the step.
      * `pipeline`  -- the shipped per-cell band (patch quantised to powers of
                       two).  Any jump there is the resolved band's lower edge
                       moving at a patch doubling, NOT a step in B(k); both are
                       reported so the two mechanisms cannot be confused.
    """
    def run(nn):
        Us = np.logspace(math.log10(U_lo), math.log10(U_hi), int(nn))
        r = _rule_scan(F, h, target, Us)
        a = np.array([(x[0], x[1], x[3], x[4]) for x in r["rows"]])
        d_am = np.abs(np.diff(a[:, 2])) / np.maximum(np.abs(a[:, 2][:-1]), 1e-300)
        d_s2 = np.abs(np.diff(a[:, 3])) / np.maximum(np.abs(a[:, 3][:-1]), 1e-300)
        return {"n": int(nn), "I_fixed": r["I_fixed"],
                "floor_fixed": r["floor_fixed"],
                "max_rel_jump_alpha_m_used": float(np.max(d_am)),
                "max_rel_jump_sigma2_h": float(np.max(d_s2)),
                "rows": r["rows"]}

    coarse = run(n)
    fine = run(int(n) * 2)
    ratio = (float(fine["max_rel_jump_alpha_m_used"]
                   / coarse["max_rel_jump_alpha_m_used"])
             if coarse["max_rel_jump_alpha_m_used"] > 0 else None)

    # the hand-off itself: raw alpha_m * I = target, bisected in U10
    def gap(U):
        r = _rule_scan(F, h, target, [U])
        U_, us, raw, used, s2 = r["rows"][0]
        return raw * r["I_fixed"] - float(target), us, raw, r["I_fixed"]

    lo, hi = float(U_lo), float(U_hi)
    g_lo = gap(lo)[0]
    g_hi = gap(hi)[0]
    handoff = None
    if g_lo < 0.0 < g_hi:
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if gap(mid)[0] < 0.0:
                lo = mid
            else:
                hi = mid
        U_h = 0.5 * (lo + hi)
        g, us, raw, I_f = gap(U_h)
        handoff = {"U10_ms": float(U_h), "ustar_ms": float(us),
                   "u_over_fit_min": float(us / USTAR_FIT_MIN),
                   "u_over_cm": float(us / C_M),
                   "raw_alpha_m_at_handoff": float(raw),
                   "alpha_m_used_at_handoff": float(max(raw, float(target) / I_f)),
                   "sigma2_h_raw": float(raw * I_f),
                   "target": float(target),
                   "band": "fixed at the U=3 cell's patch (rule-only scan)",
                   "method": ("bisection on U10 at fixed (F,h) for "
                              "sigma2_h,raw(u*) = target")}
        # one-sided log-slopes: max(a,b) is C0 but not C1 at the crossing
        eps = 2e-4
        pts = {}
        for tag, uu in (("a", U_h * (1.0 - 2 * eps)),
                        ("b", U_h * (1.0 - eps)),
                        ("c", U_h * (1.0 + eps)),
                        ("d", U_h * (1.0 + 2 * eps))):
            r = _rule_scan(F, h, target, [uu])["rows"][0]
            pts[tag] = (math.log(r[0]), math.log(r[3]), r[3])
        sl_l = (pts["b"][1] - pts["a"][1]) / (pts["b"][0] - pts["a"][0])
        sl_r = (pts["d"][1] - pts["c"][1]) / (pts["d"][0] - pts["c"][0])
        handoff["kink"] = {
            "alpha_at_handoff_left": float(pts["b"][2]),
            "alpha_at_handoff_right": float(pts["c"][2]),
            "rel_gap_across": float(pts["c"][2] / pts["b"][2] - 1.0),
            "dln_alpha_dln_U_left": float(sl_l),
            "dln_alpha_dln_U_right": float(sl_r),
            "slope_jump": float(sl_r - sl_l),
            "statement": ("max(raw, floor) is continuous (C0: the two sides "
                          "meet at the crossing) and its derivative jumps "
                          "there (C1 fails by construction of max()).  Both "
                          "one-sided log-slopes and the value gap are "
                          "measured, not asserted."),
        }

    # the shipped per-cell band (patch quantised)
    Us = np.logspace(math.log10(U_lo), math.log10(U_hi), int(n) * 2)
    pr = []
    for U in Us:
        ctx = assemble(float(U), float(F), float(h), arm="B",
                       sigma2_h_target=target)
        r = ctx["alpha_m_resolution"]
        pr.append((float(U), float(r["ustar_ms"]),
                   float(r["alpha_m_used"]),
                   float(r["sigma2_h_total_closure"]),
                   float(r["sigma2_h_unit_integral"]),
                   bool(r["top_up_acted"])))
    pa = np.array([(x[0], x[2], x[3], x[4]) for x in pr])
    pd_am = np.abs(np.diff(pa[:, 1])) / np.maximum(np.abs(pa[:, 1][:-1]), 1e-300)
    pd_s2 = np.abs(np.diff(pa[:, 2])) / np.maximum(np.abs(pa[:, 2][:-1]), 1e-300)
    pd_I = np.abs(np.diff(pa[:, 3])) / np.maximum(np.abs(pa[:, 3][:-1]), 1e-300)
    i_am = int(np.argmax(pd_am))
    i_I = int(np.argmax(pd_I))
    acted = [x[5] for x in pr]
    idx = [i for i, v in enumerate(acted) if v]

    return {
        "fetch_m": float(F), "depth_m": float(h), "target": float(target),
        "U_range_ms": [float(U_lo), float(U_hi)],
        "rule_scan": {
            "what": ("max(raw, floor) with the band held fixed: the hand-off "
                     "rule alone"),
            "coarse": {"n": coarse["n"],
                       "max_rel_jump_alpha_m_used":
                           coarse["max_rel_jump_alpha_m_used"],
                       "max_rel_jump_sigma2_h":
                           coarse["max_rel_jump_sigma2_h"]},
            "fine": {"n": fine["n"],
                     "max_rel_jump_alpha_m_used":
                         fine["max_rel_jump_alpha_m_used"],
                     "max_rel_jump_sigma2_h": fine["max_rel_jump_sigma2_h"]},
            "ratio_fine_over_coarse": ratio,
            "refinement_bound": 0.65,
            "continuous_C0": bool(ratio is not None and ratio <= 0.65),
            "I_fixed": fine["I_fixed"], "floor_fixed": fine["floor_fixed"],
            "note": ("the ratio ~0.5 says the ramp is smooth, not a step; the "
                     "residual per-step change is the ramp's own steepness")},
        "handoff": handoff,
        "pipeline_scan": {
            "what": ("the shipped per-cell band: patch quantised to powers of "
                     "two, so I(cell) and the resolved-band integrals move in "
                     "steps at the doublings"),
            "n": int(n) * 2,
            "max_rel_jump_alpha_m_used": float(pd_am.max()),
            "U_at_max_alpha_jump": float(pa[i_am][0]),
            "max_rel_jump_sigma2_h": float(pd_s2.max()),
            "U_at_max_sigma2_jump": float(pa[int(np.argmax(pd_s2))][0]),
            "max_rel_jump_I": float(pd_I.max()),
            "U_at_max_I_jump": float(pa[i_I][0]),
            "top_up_acts_from_U_ms": (float(pr[idx[0]][0]) if idx else None),
            "top_up_stops_at_U_ms": (float(pr[idx[-1]][0]) if idx else None),
            "note": ("these steps are the patch doubling (the resolved band's "
                     "lower edge), a pre-existing v11/v13 property; B(k) "
                     "itself is continuous and no arm changes this"),
        },
        "pass_C0": bool(ratio is not None and ratio <= 0.65),
        "scan": [{"U_ms": x[0], "ustar_ms": x[1], "alpha_m_used": x[2],
                  "sigma2_h": x[3], "I_cell": x[4], "top_up_acted": x[5]}
                 for x in pr],
    }


# --------------------------------------------------------------------------
# 4. per-arm cell evaluation
# --------------------------------------------------------------------------
def evaluate_cell(U, F, h, arm, k_cut=None, sigma2_h_target=None, ctx=None):
    """One arm at one cell: every number the report's rows need."""
    tgt = (SIGMA2_H_TARGET_DEFAULT if sigma2_h_target is None
           else float(sigma2_h_target))
    if ctx is None:
        ctx = assemble(float(U), float(F), float(h), arm=arm, k_cut=k_cut,
                       sigma2_h_target=tgt)
    geo = cell_geometry(ctx)
    kmin = geo["kmin_rad_per_m"]
    kvisc = geo["k_viscous_rad_per_m"]
    kc = float(ctx["k_cut_rad_per_m"])
    kf = V11._k_grid(kvisc, k_hi=kvisc)
    an_full = V11.integral_block(kf, V11.s_parts(kf, ctx), kmin, kvisc,
                                 float(ctx["k_p_rad_per_m"]))
    kcut = V11._k_grid(kc)
    an_cut = V11.integral_block(kcut, V11.s_parts(kcut, ctx), kmin, kc,
                                float(ctx["k_p_rad_per_m"]))
    mom = moments(ctx)
    sh = shares(ctx, an_cut)
    env = envelope(U, F, h, kmin, kc)
    ec = envelope_check(an_cut["total"]["sigma2"], env)
    be = b_envelope(ctx)
    sm = smoothness(ctx)
    pos = positivity(ctx)
    hs = hs_node(ctx)
    c1 = gate_c1(sh)
    c1b = gate_c1b(sh, mom)
    U_f = float(U)
    low_wind = bool(U_f <= 3.0 + 1e-12)

    # ---- three-state verdict, per PIPELINE-V12 -----------------------------
    fails = []
    if low_wind:
        if not ec["tol_inside"]:
            fails.append("totals outside [0.5x, 2.0x] of the envelope")
        if not be["pass_B_env"]:
            fails.append("B(k) > 0.012 somewhere in [10,1500]")
        if not pos["pass"]:
            fails.append("negative or non-finite energy")
        if not sm["pass"]:
            fails.append("step in S(k) / B(k)")
        if hs["pass"] is False:
            fails.append("|dHs| > 35%")
        state = "FAIL" if fails else "UNRESOLVED"
        rule = ("low-wind cell (U<=3): at best UNRESOLVED.  FAIL only on the "
                "design's explicit list (totals outside 0.5x-2.0x of the "
                "envelope, B > 0.012, negative energy, a step in S(k), "
                "|dHs| > 35%).  UNRESOLVED is not a soft pass.")
    else:
        checks = {"C1": c1["pass"], "C1b": c1b["pass"],
                  "totals_tol": ec["tol_inside"], "B_envelope": be["pass"]}
        bad = [k for k, v in checks.items() if v is not True]
        state = "PASS" if not bad else "FAIL"
        rule = ("high-wind cell (U=6/10): must PASS -- C1 (gravity share "
                "<= 0.50), C1b (capillary share >= 0.15 AND sigma2_h >= "
                "0.002), totals inside the 0.5x-2.0x envelope, B-envelope "
                "(B <= 0.012 over 10-1500, B_l <= 0.008 over 50-200).  "
                "Failing clauses: %s" % (", ".join(bad) or "none"))
        fails = bad
    absolute = None
    if low_wind:
        absolute = {
            "label": PROVENANCE["D"],
            "B_max_10_1500": be["B_max"], "B_env_max": B_ENV_MAX,
            "pass_B_env": be["pass_B_env"],
            "B_l_max_50_200": be["B_l_max"], "beta_hat": BETA_HAT,
            "pass_B_l": be["pass_B_l"],
            "totals_strict_inside": ec["strict_inside"],
            "totals_ratio_to_lower": ec["ratio_to_lower"],
            "totals_ratio_to_upper": ec["ratio_to_upper"],
            "totals_tol_inside": ec["tol_inside"],
            "smoothness_pass": sm["pass"],
            "n_clauses": 4,
            "n_pass_strict": int(sum([bool(be["pass_B_env"]),
                                      bool(be["pass_B_l"]),
                                      bool(ec["strict_inside"]),
                                      bool(sm["pass"])])),
            "n_pass_tolerance": int(sum([bool(be["pass_B_env"]),
                                         bool(be["pass_B_l"]),
                                         bool(ec["tol_inside"]),
                                         bool(sm["pass"])])),
            "verdict_strict": ("PASS" if all([be["pass_B_env"],
                                              be["pass_B_l"],
                                              ec["strict_inside"],
                                              sm["pass"]]) else "FAIL"),
            "verdict_tolerance": ("PASS" if all([be["pass_B_env"],
                                                 be["pass_B_l"],
                                                 ec["tol_inside"],
                                                 sm["pass"]])
                                  else "FAIL"),
            "note": ("Gemini's absolute criteria, tested on this arm's "
                     "spectrum.  Both the strict totals reading and the "
                     "0.5x-2.0x reading are scored; neither is substituted "
                     "for the other."),
        }
    return {
        "arm": arm,
        "cell": {"U_ms": U_f, "fetch_m": float(F), "depth_m": float(h)},
        "low_wind_cell": low_wind,
        "geometry": geo,
        "sea_state": {k: ctx["sea_state"][k] for k in (
            "k_p_rad_per_m", "c_p_ms", "omega_c", "alpha_p", "gamma_bl",
            "sigma_bl", "eckv_domain", "k_p_dispersion_rel_diff",
            "kp_times_h")},
        "alpha": ctx["alpha_m_resolution"],
        "provenance": ctx["provenance"],
        "moments": mom,
        "integrals_cut": an_cut,
        "integrals_full": an_full,
        "shares": sh,
        "b_envelope": be,
        "envelope": env,
        "envelope_check": ec,
        "hs": hs,
        "C1": c1, "C1b": c1b,
        "smoothness": sm,
        "positivity": pos,
        "absolute_criteria": absolute,
        "state": state,
        "fail_reasons": fails,
        "verdict_rule": rule,
        "Hs_spectrum_m": float(ctx["Hs_m"]),
        "Hs_law_m": float(ctx["Hs_law_m"]),
        "sigma2_total_cut": float(an_cut["total"]["sigma2"]),
        "sigma2_total_full": float(an_full["total"]["sigma2"]),
    }


def export_curves(ctx, n=4000):
    """S/B on a dense log grid -- enough to redo every band integral."""
    geo = cell_geometry(ctx)
    k = np.logspace(math.log10(geo["kmin_rad_per_m"]),
                    math.log10(B_K_RANGE[1]), int(n))
    s = V11.s_parts(k, ctx)
    return {"k_rad_per_m": [float(x) for x in k],
            "B_l": [float(x) for x in s["B_l"]],
            "B_h": [float(x) for x in s["B_h"]],
            "S_l_m2_per_rad_m": [float(x) for x in s["S_l"]],
            "S_h_m2_per_rad_m": [float(x) for x in s["S_h"]],
            "band_rad_per_m": [geo["kmin_rad_per_m"], float(B_K_RANGE[1])],
            "n": int(n)}


def reproduce_from_export(curve, k_cut):
    """Re-integrate an exported curve -- the A4 'reproduce it yourself' check."""
    k = np.asarray(curve["k_rad_per_m"])
    S = np.asarray(curve["S_l_m2_per_rad_m"]) + np.asarray(
        curve["S_h_m2_per_rad_m"])
    m = k <= float(k_cut)
    kk = k[m]
    tot = float(np.trapz(kk ** 2 * S[m], kk))
    mo = m & (k >= BAND[0]) & (k <= BAND[1])
    ko = k[mo]
    Bl = np.asarray(curve["B_l"])[mo]
    Bh = np.asarray(curve["B_h"])[mo]
    g2 = float(np.trapz(ko ** 2 * Bl / ko ** 3, ko))
    h2 = float(np.trapz(ko ** 2 * Bh / ko ** 3, ko))
    return {"sigma2_total": tot, "sigma2_gravity_50_200": g2,
            "sigma2_capillary_50_200": h2}
