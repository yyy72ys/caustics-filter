"""v13 spectrum -- version 3.9 (PIPELINE-V11 chunk 9): the low-wind capillary
floor, honestly labelled.

Base: `v11_spectrum` (v3.8) is IMPORTED and never edited.  v13 changes exactly
one quantity -- the level `alpha_m` of the Elfouhaily high-frequency part B_h:

    inside the ECKV fit domain (u* >= c_m/e):
        alpha_m_used = alpha_m_raw            (as ECKV gives it, no max())
    below it (u* < c_m/e), with a positive target:
        alpha_m_used = max(alpha_m_raw, alpha_m_floor(cell))
        alpha_m_floor(cell) solves sigma2_h,total = sigma2_h_target for this cell

`sigma2_h,total` is the capillary part's slope variance over the cell's own
resolved band [kmin, k_vis], with kmin = 2 pi / patch and k_vis from
v8_physics.viscous_cutoff.  Because B_h is exactly linear in alpha_m
(V10.elfouhaily_bh), the solve is closed form

    alpha_m_floor = sigma2_h_target / I(cell),
    I(cell) = int_{kmin}^{k_vis} k^2 * B_h(k; alpha_m = 1) / k^3 dk,

and the closure sigma2_h(alpha_m_floor) = target is verified on a SEPARATE code
path (assemble -> V11.s_parts -> trapz) and reported as `floor_proof`.

What is retired: the v3.8 label `eckv_domain = "low_wind_floor"` (the hard
alpha_m zero as a physical claim).  v13 emits `alpha_m_domain` in

    "eckv_fit_domain"                    u* >= c_m/e, raw used
    "low_wind_prior"                     u* <  c_m/e, prior active (target > 0)
    "below_fit_domain_prior_disabled"    u* <  c_m/e, target == 0 (sweep control)

and never emits "low_wind_floor".  A state built with the prior active carries
`low_wind_prior` with the design's target/source/rally fields verbatim.

Nothing else is touched: S(k) = (B_l + B_h)/k^3, the band table, the ECKV
reference, Phi, the patch/cascade rule and the transport are the v11 code.
"""
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10
import v11_spectrum as V11

G = P8.G
TF = P8.TF
C_M = V11.C_M
K_M = V11.K_M

#: the ECKV fit crossing: alpha_m_raw = 0 at u* = c_m/e (0.08505 m/s).
USTAR_FIT_MIN = float(C_M / math.e)

#: the design's provisional prior and its mandatory sweep points.
SIGMA2_H_TARGET_DEFAULT = 0.003
SIGMA2_H_TARGETS = (0.0, 0.0015, 0.003, 0.006)

#: the design's label, verbatim (PIPELINE-V11 section 2).
LOW_WIND_PRIOR_SOURCE = ("Cox-Munk crosswind intercept \u2014 ocean residual, "
                         "mechanism NOT claimed")
LOW_WIND_PRIOR_RALLY = "v3.9 R3 weakened at 3/3"

#: alpha_m_domain labels (the retired one is kept here only as a constant so a
#: test can prove it is never emitted).
ALPHA_M_DOMAIN_FIT = "eckv_fit_domain"
ALPHA_M_DOMAIN_PRIOR = "low_wind_prior"
ALPHA_M_DOMAIN_PRIOR_OFF = "below_fit_domain_prior_disabled"
RETIRED_ALPHA_M_DOMAIN = "low_wind_floor"

SOURCES = V11.SOURCES
BAND_EDGES_V13 = V11.BAND_EDGES_V11
BAND_NAMES_V13 = V11.BAND_NAMES_V11


def low_wind_prior_label(sigma2_h_target):
    """The design's `low_wind_prior` label, target set to the value in use.

    At the design's provisional value (0.003) this is the label verbatim; the
    sensitivity sweep re-uses the same structure with its own target, because a
    label that printed 0.003 while a different target was applied would be a
    false label.
    """
    return {"sigma2_h_target": float(sigma2_h_target),
            "source": LOW_WIND_PRIOR_SOURCE,
            "rally": LOW_WIND_PRIOR_RALLY}


def cell_kmin(lambda_p_m):
    """The cell's own resolved-band lower edge (v11 patch rule)."""
    pm = V10.patch_for(float(lambda_p_m))["patch_m"]
    return 2.0 * math.pi / float(pm), float(pm)


def _k_band(ctx, k_hi=None):
    khi = float(k_hi if k_hi is not None else
                P8.viscous_cutoff(float(ctx["depth_m"])))
    return V11._k_grid(khi, k_hi=khi), khi


def alpha_m_unit_integral(ctx, kmin=None, k_hi=None):
    """I(cell) = int k^2 B_h(k; alpha_m = 1)/k^3 dk over [kmin, k_vis]."""
    k, khi = _k_band(ctx, k_hi)
    pm = cell_kmin(ctx["laws"]["lambda_p_m"])[0] if kmin is None else kmin
    m = k >= float(pm)
    kk = k[m]
    bh1 = V10.elfouhaily_bh(kk, float(ctx["depth_m"]),
                            float(ctx["ustar"]["ustar_ms"]),
                            float(ctx["k_p_rad_per_m"]), 1.0)
    s_h1 = bh1 / kk ** 3
    return {"integral": float(np.trapz(kk ** 2 * s_h1, kk)),
            "band_rad_per_m": [float(pm), float(khi)],
            "k_vis_rad_per_m": float(khi),
            "note": ("B_h is exactly linear in alpha_m, so this single "
                     "integral holds the whole dependence: "
                     "sigma2_h = alpha_m * I(cell)")}


def alpha_m_resolve(ctx, sigma2_h_target=SIGMA2_H_TARGET_DEFAULT,
                    unit=None):
    """raw / floor / used alpha_m for one cell, with the domain label."""
    target = float(sigma2_h_target)
    depth = float(ctx["depth_m"])
    ustar = float(ctx["ustar"]["ustar_ms"])
    kp = float(ctx["k_p_rad_per_m"])
    v10 = V10.alpha_m_elfouhaily(ustar)
    raw = float(v10["alpha_m_raw"])
    u = unit if unit is not None else alpha_m_unit_integral(ctx)
    I = float(u["integral"])
    floor = (target / I) if (target > 0.0 and I > 0.0) else 0.0
    in_domain = bool(ustar >= USTAR_FIT_MIN)
    if in_domain:
        used = raw
        domain = ALPHA_M_DOMAIN_FIT
        prior = None
        rule = ("u* >= c_m/e: alpha_m is used RAW as ECKV gives it "
                "(the v3.8 hard max(raw, 0) is gone)")
    elif target > 0.0:
        used = max(raw, floor)
        domain = ALPHA_M_DOMAIN_PRIOR
        prior = low_wind_prior_label(target)
        rule = ("u* < c_m/e: alpha_m_used = max(alpha_m_raw, "
                "alpha_m_floor(cell)) with alpha_m_floor solving "
                "sigma2_h,total = %.6g" % target)
    else:
        used = max(raw, 0.0)
        domain = ALPHA_M_DOMAIN_PRIOR_OFF
        prior = None
        rule = ("u* < c_m/e, target = 0: the sweep's control point -- the "
                "v3.8 hard zero, kept ONLY here as the sensitivity "
                "baseline")
    if domain == RETIRED_ALPHA_M_DOMAIN:      # pragma: no cover - defensive
        raise RuntimeError("the retired v3.8 label was emitted")
    return {
        "u_ms": float(ctx["laws"]["U10_ms"]),
        "ustar_ms": ustar,
        "c_m_ms": float(C_M),
        "ustar_fit_min_ms": USTAR_FIT_MIN,
        "u_over_cm": float(ustar / C_M),
        "u_over_fit_min": float(ustar / USTAR_FIT_MIN),
        "in_eckv_fit_domain": in_domain,
        "alpha_m_raw": raw,
        "alpha_m_floor": float(floor),
        "alpha_m_used": float(used),
        "alpha_m_domain": domain,
        "prior_active": bool(prior is not None),
        "low_wind_prior": prior,
        "sigma2_h_target": target,
        "sigma2_h_unit_integral": I,
        "sigma2_h_total_closure": float(used * I),
        "sigma2_h_raw_would_be": float(raw * I),
        "alpha_m_v11_hard_zero": float(v10["alpha_m"]),
        "rule": rule,
        "domain_rule": ("the ECKV fit crossing is u* = c_m/e = %.6f m/s; "
                        "c_m = %.6f m/s from K_M = 2 pi/0.017 m.  A real "
                        "capillary collapse near u* ~ 0.04-0.08 m/s is a "
                        "DIFFERENT physical statement this prior does not "
                        "make (Gemini, PIPELINE-V11 section 2)."
                        % (USTAR_FIT_MIN, C_M)),
    }


def assemble(U, F, h, hs_override=None, k_cut=None, gamma=None,
             source="elfouhaily-form",
             sigma2_h_target=SIGMA2_H_TARGET_DEFAULT):
    """The v11 spectrum with the v3.9 alpha_m resolution applied.

    source="jonswap" returns the v11 legacy ctx BYTE-FOR-BYTE (the legacy path
    is untouched); source="elfouhaily-form" patches alpha_m_value and the
    labels, then re-predicts H_s from the assembled spectrum.
    """
    ctx = V11.assemble(U, F, h, hs_override=hs_override, k_cut=k_cut,
                       gamma=gamma, source=source)
    if source != "elfouhaily-form":
        return ctx
    res = alpha_m_resolve(ctx, sigma2_h_target=sigma2_h_target)
    ctx["version"] = "3.9"
    ctx["alpha_m_value"] = float(res["alpha_m_used"])
    ctx["alpha_m_used"] = float(res["alpha_m_used"])
    ctx["alpha_m_raw"] = float(res["alpha_m_raw"])
    ctx["alpha_m_floor"] = float(res["alpha_m_floor"])
    ctx["alpha_m_domain"] = res["alpha_m_domain"]
    ctx["alpha_m_resolution"] = res
    ctx["low_wind_prior"] = res["low_wind_prior"]
    ctx["sigma2_h_target"] = float(res["sigma2_h_target"])
    ctx["sigma2_h_total"] = float(res["sigma2_h_total_closure"])
    ctx["label"] = (ctx["label"] + " -- v3.9: low-wind capillary prior below "
                                    "the ECKV fit domain, labelled, swept")
    # H_s is a PREDICTION of the assembled spectrum (PIPELINE-V10 section 4):
    # re-predict it with the resolved alpha_m, same rule as v11.
    kcut = float(ctx["k_cut_rad_per_m"])
    k = V11._k_grid(kcut, k_hi=float(P8.viscous_cutoff(float(h))))
    s = V11.s_parts(k, ctx)
    kmin = cell_kmin(ctx["laws"]["lambda_p_m"])[0]
    m = k >= kmin
    ctx["Hs_m"] = float(4.0 * math.sqrt(
        max(float(np.trapz(s["S"][m], k[m])), 0.0)))
    ctx["Hs_pred_rel_law"] = float(ctx["Hs_m"] / ctx["Hs_law_m"] - 1.0)
    return ctx


def s_parts(k, ctx):
    return V11.s_parts(k, ctx)


def s_total(k, ctx):
    return V11.s_total(k, ctx)


def mode_density_2d(*a, **kw):
    return V11.mode_density_2d(*a, **kw)


def curves(ctx, patch_m, k_cut=None, n=420):
    return V11.curves(ctx, patch_m, k_cut=k_cut, n=n)


def gate_a2(ctx, tag="elfouhaily-form (v3.9)", **kw):
    return V11.gate_a2(ctx, tag=tag, **kw)


def gate_k6(ctx, **kw):
    return V11.gate_k6(ctx, **kw)


# --------------------------------------------------------------------------
# the floor proof (C1b's second clause) and the C1 three-state
# --------------------------------------------------------------------------
def floor_proof(ctx, k_cut=None):
    """sigma2_h,total of the ASSEMBLED spectrum, recomputed independently.

    Path: assemble -> V11.s_parts (B_h -> B_h/k^3) -> trapz.  This is a
    different code path from alpha_m_unit_integral, so the residual between
    sigma2_h_total and the target is a real closure test of the floor.
    """
    h = float(ctx["depth_m"])
    khi = float(P8.viscous_cutoff(h))
    k = V11._k_grid(khi, k_hi=khi)
    s = V11.s_parts(k, ctx)
    kmin = cell_kmin(ctx["laws"]["lambda_p_m"])[0]
    m = k >= kmin
    kk = k[m]
    cap = float(np.trapz(kk ** 2 * s["S_h"][m], kk))
    grav = float(np.trapz(kk ** 2 * s["S_l"][m], kk))
    target = float(ctx.get("sigma2_h_target") or 0.0)
    active = bool(ctx.get("low_wind_prior") is not None)
    rel = (None if (not active or target <= 0.0)
           else float(cap / target - 1.0))
    kc = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    kcut = V11._k_grid(kc)
    sc = V11.s_parts(kcut, ctx)
    mc = kcut >= kmin
    kkc = kcut[mc]
    cap_cut = float(np.trapz(kkc ** 2 * sc["S_h"][mc], kkc))
    return {
        "band_rad_per_m": [float(kmin), float(khi)],
        "sigma2_h_total": cap,
        "sigma2_gravity_total": grav,
        "sigma2_h_total_at_kcut": cap_cut,
        "k_cut_rad_per_m": kc,
        "target": target,
        "prior_active": active,
        "closure_rel_residual": rel,
        "closure_abs_residual": (None if rel is None
                                 else float(cap - target)),
        "pass": (bool(rel is not None and abs(rel) <= 1.0e-6)
                 if active else None),
        "clause2_reference": 0.002,
        "supplies_short_wave_variance": bool(cap >= 0.002),
        "note": ("sigma2_h,total = the capillary part's slope variance over "
                 "the cell's resolved band [kmin, k_vis]; the prior claims "
                 "this number, not a measured pond level"),
    }


def c1_three_state(ctx, k_cut=None):
    """C1 per PIPELINE-V11: PASS/FAIL in the fit domain, UNRESOLVED below."""
    kc = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    k = V11._k_grid(kc)
    s = V11.s_parts(k, ctx)
    kmin = cell_kmin(ctx["laws"]["lambda_p_m"])[0]
    ib = V11.integral_block(k, s, kmin, kc, float(ctx["k_p_rad_per_m"]))
    o = ib["overlap_50_200"]
    grav_share = float(o["gravity_share_of_sigma2"])
    cap_share = float(1.0 - grav_share)
    in_domain = bool(float(ctx["ustar"]["ustar_ms"]) >= USTAR_FIT_MIN)
    pass_share = bool(grav_share <= 0.50)
    pass_maxb = bool(o["B_l_max"] <= 0.010)
    if in_domain:
        state = "PASS" if (pass_share and pass_maxb) else "FAIL"
    else:
        state = "UNRESOLVED"
    return {
        "state": state,
        "in_eckv_fit_domain": in_domain,
        "gravity_share_50_200": grav_share,
        "capillary_share_50_200": cap_share,
        "B_l_max_50_200": float(o["B_l_max"]),
        "B_h_max_50_200": float(o["B_h_max"]),
        "B_l_sigma2_50_200": float(o["B_l_sigma2"]),
        "B_h_sigma2_50_200": float(o["B_h_sigma2"]),
        "limits": {"gravity_share_max": 0.50, "B_l_max_in_band": 0.010},
        "pass_share": pass_share,
        "pass_maxB": pass_maxb,
        "counted_as_pass": bool(state == "PASS"),
        "note": ("three-state per the v3.9 rally consensus: below u* = "
                 "c_m/e = %.6f m/s C1 is UNRESOLVED -- reported with its "
                 "numbers, never green, never waived" % USTAR_FIT_MIN),
    }


def c1b(ctx, fp=None):
    """Grok's floor proof: capillary share of 50-200 >= 15% OR
    sigma2_h,total >= 0.002, on ALL cells."""
    c1 = c1_three_state(ctx)
    f = fp if fp is not None else floor_proof(ctx)
    cap_share = float(c1["capillary_share_50_200"])
    ok1 = bool(cap_share >= 0.15)
    ok2 = bool(f["sigma2_h_total"] >= 0.002)
    return {
        "gate": ("C1b -- one wind-independent minimum-capillary test for ALL "
                 "cells (Grok's condition)"),
        "capillary_share_50_200": cap_share,
        "share_threshold": 0.15,
        "share_pass": ok1,
        "sigma2_h_total": float(f["sigma2_h_total"]),
        "sigma2_threshold": 0.002,
        "sigma2_pass": ok2,
        "prior_active": bool(f["prior_active"]),
        "alpha_m_domain": ctx["alpha_m_domain"],
        "pass": bool(ok1 or ok2),
        "note": ("PASS if EITHER clause holds.  The prior must do measurable "
                 "work everywhere or the chunk fails; the numbers say which "
                 "clause carried each cell."),
    }


# --------------------------------------------------------------------------
# references reported beside the total (C2's replacement; not gates)
# --------------------------------------------------------------------------
def fetch_limited_mss_reference(U, F, h, kmin=None, k_hi=None):
    """In-repo fetch-limited slope variance (JONSWAP-family, same band).

    The design names 'the Hu law' for this column.  No Hu-law formula exists
    in this workspace or in either design document, so the column is filled
    with the pipeline's own fetch-limited MSS integral and the naming gap is
    recorded as UNRESOLVED (reporting is not a gate).
    """
    rel = P8.jonswap_fetch_relations(float(U), float(F))
    khi = float(k_hi if k_hi is not None else V10.KCUT_DEFAULT)
    lo = float(kmin if kmin is not None else 1.0e-3)
    out = {}
    for name, gam in (("jonswap_gamma_3p3", 3.3), ("donelan_gamma_1p7", 1.7)):
        r = P8._sigma2_integral(float(U), float(F), float(h), gam, khi,
                                n=40000)
        out[name] = r["sigma2"]
    # band-limited variant on the cell's own resolved band
    k = np.logspace(math.log10(max(lo, 1e-3)), math.log10(khi), 24000)
    om = P8.dispersion_omega(k, float(h))
    f = om / (2.0 * math.pi)
    fp = rel["f_p_Hz"]
    a = rel["alpha"]
    sj = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-((f - fp) ** 2) / (2.0 * sj ** 2 * fp ** 2))
    sf = (a * G ** 2 / (2.0 * math.pi) ** 4 * f ** -5
          * np.exp(-1.25 * (fp / np.maximum(f, 1e-12)) ** 4) * 3.3 ** r)
    sk = sf * P8.dispersion_domega_dk(k, float(h)) / (2.0 * math.pi)
    out["band_limited_gamma_3p3"] = float(np.trapz(k * k * sk, k))
    return {
        "what": ("the pipeline's own fetch-limited MSS integral "
                 "(v8_physics._sigma2_integral, JONSWAP/Donelan family)"),
        "design_name": "the Hu law",
        "naming_status": ("UNRESOLVED -- the design names 'the Hu law'; no "
                          "Hu-law formula is present in this workspace or in "
                          "the two design documents, so this column is the "
                          "in-repo fetch-limited MSS, NOT a Hu "
                          "parameterisation"),
        "k_hi_rad_per_m": khi,
        "band_limited_kmin_rad_per_m": lo,
        "sigma2_full": out,
        "note": ("a reference beside Cox-Munk, never a target; the full-band "
                 "numbers use the same k_hi as the model column"),
    }


def cox_munk_reference(U):
    return V10.cox_munk_reference(U)


# --------------------------------------------------------------------------
# algebra checks against the v11 conventions
# --------------------------------------------------------------------------
def algebra_vs_v11(U, F, h, kmin=None, k_cut=None, sigma2_h_target=None):
    """v13 vs v11 on the same cell: identical where the prior is inactive.

    Three statements, all measured:
      1. raw alpha_m: v13 == v10.alpha_m_elfouhaily's raw value (same formula).
      2. S(k): elementwise IDENTICAL to v11 where the prior is inactive;
         where it is active, the ratio S_h_v13/S_h_v11 equals alpha_used/alpha_v11
         and S_l is untouched.
      3. linearity: sigma2_h(2 alpha) = 2 sigma2_h(alpha) to machine precision.
    """
    tgt = float(SIGMA2_H_TARGET_DEFAULT if sigma2_h_target is None
                else sigma2_h_target)
    ctx13 = assemble(U, F, h, k_cut=k_cut, sigma2_h_target=tgt)
    ctx11 = V11.assemble(U, F, h, k_cut=k_cut, source="elfouhaily-form")
    khi = float(P8.viscous_cutoff(float(h)))
    k = V11._k_grid(khi, k_hi=khi)
    s13 = V11.s_parts(k, ctx13)
    s11 = V11.s_parts(k, ctx11)
    active = ctx13["low_wind_prior"] is not None
    am11 = float(ctx11["alpha_m_value"])
    ratio = (float(ctx13["alpha_m_value"] / am11) if am11 > 0 else None)
    dS = np.abs(s13["S"] - s11["S"])
    dSl = np.abs(s13["S_l"] - s11["S_l"])
    dSh = np.abs(s13["S_h"] - s11["S_h"])
    scale = max(float(np.max(s13["S"])), 1e-300)
    # linearity, on the v13 ctx
    bh1 = V10.elfouhaily_bh(k, float(h), float(ctx13["ustar"]["ustar_ms"]),
                            float(ctx13["k_p_rad_per_m"]), 1.0)
    km = cell_kmin(ctx13["laws"]["lambda_p_m"])[0] if kmin is None else kmin
    m = k >= float(km)
    kk = k[m]
    cap1 = float(np.trapz(kk ** 2 * (bh1[m] / kk ** 3), kk))
    cap2 = float(np.trapz(kk ** 2 * (2.0 * bh1[m] / kk ** 3), kk))
    return {
        "cell": {"U_ms": float(U), "fetch_m": float(F), "depth_m": float(h)},
        "sigma2_h_target": tgt,
        "prior_active": bool(active),
        "alpha_m_raw_v10": float(V10.alpha_m_elfouhaily(
            float(ctx11["ustar"]["ustar_ms"]))["alpha_m_raw"]),
        "alpha_m_raw_v13": float(ctx13["alpha_m_raw"]),
        "alpha_m_used_v13": float(ctx13["alpha_m_value"]),
        "alpha_m_v11": float(ctx11["alpha_m_value"]),
        "alpha_ratio_v13_over_v11": ratio,
        "max_abs_dS_over_maxS": float(dS.max() / scale),
        "max_abs_dS_gravity": float(dSl.max()),
        "spectrum_identical": bool(not dS.any()) if not active else False,
        "gravity_part_identical": bool(not dSl.any()),
        "S_h_ratio_matches_alpha_ratio": bool(
            not active or ratio is None or np.allclose(
                s13["S_h"][s11["S_h"] > 0] / s11["S_h"][s11["S_h"] > 0],
                ratio, rtol=1e-12)),
        "linearity_sigma2h_2alpha_over_alpha": float(cap2 / cap1),
        "linearity_pass": bool(abs(cap2 / cap1 - 2.0) <= 1e-9),
        "max_abs_dS_h": float(dSh.max()),
        "note": ("where the prior is inactive the v3.9 source reproduces the "
                 "v3.8 source bit-for-bit in the assembled S(k); where it is "
                 "active only B_h moves, by the alpha_m ratio"),
    }
