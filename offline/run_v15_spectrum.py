"""run_v15_spectrum -- version 3.10 (chunk 11): THE MERGE.

Brief: BRIEF-V15-MERGE-R2.md.  Base: `run_v14_spectrum` (which carries v13's
spectrum and v11's assembly).  v5..v14 are IMPORTED, never edited.

What changes, and what does not
--------------------------------
1. The default spectrum source of THIS family is exposed as

       SPECTRUM_SOURCE_DEFAULT = "elfouhaily-form"

   with the Arm-D rule set (v13's raw alpha_m inside the ECKV fit domain, and
   `max(raw, floor)` below it; `floor` solves sigma2_h,total = 0.003).  The
   provenance strings `eckv_fit_domain_raw` / `low_wind_prior` are the v14
   ones, verbatim.  `jonswap` remains selectable (the rollback lane) and is
   returned BYTE-FOR-BYTE untouched.

2. The universal `max B(k) <= 0.012 over [10,1500]` clause is RETIRED.  It was
   a category error: 0.012 is the equilibrium-range coefficient (a MEAN, with
   a measured range 0.006-0.024), never a cap on the total saturation
   spectrum.  `b_max_observation` still PRINTS max B and max B_l -- as an
   observation, with `gated: false` in every record.

3. It is replaced by a REGION-WISE battery, reported in TWO columns per cell
   and never collapsed:

   CONFORMANCE (the model vs its own published form)
     (a) gravity:   max B_l over [k_p/2, 2 k_p]
                    vs 0.1433 * alpha_p * gamma(Omega_c) * 1.04; PASS <= 1.5
     (c) capillary: max B_h over [k_m/2, 2 k_m]
                    vs (1/2) * alpha_m(u*); PASS inside [0.8, 1.25]

   OBSERVATION (independent anchors)
     (b) equilibrium beta: E(f) reconstructed from the assembled
         omnidirectional spectrum through the deep-water dispersion
         (omega^2 = g k, f = omega / 2 pi); log-space least squares of
         E = beta * I(p) * g * u* * f^-4 over f in [1.3 f_p, 3 f_p], I(p) = 2.5.
         PASS if fitted beta in [0.006, 0.024].
         CONVENTION (stated, as the brief requires): the pre-registered band is
         the literature range of the equilibrium-range coefficient in TOBA's
         form, i.e. of S(omega) = beta * I(p) * g * u* * omega^-4 with omega in
         rad/s.  The gate uses that convention and the SAME fit under the
         literal `E(f), f in Hz` reading is printed in every cell beside it:
         it is exactly (2 pi)^3 = 248.05 times smaller and is NOT reachable by
         the pre-registered band for any spectrum of this normalisation.
         If the fit band is contaminated (B_h > B_l anywhere inside it) the
         sub-check is UNRESOLVED(band-contaminated) -- never a silent
         PASS/FAIL -- and the value is still printed with its band position.
     (d) JONSWAP anchor (context, never a gate): B_J(k_p) =
         0.5 * alpha_J * exp(-1.25) * 3.3, alpha_J = 0.076 (g F / U^2)^-0.22,
         printed beside the model's B_l(k_p) for every cell; the sub-km
         extrapolation caveat is printed for the F = 100 m cells.

4. C1 / C1b band shares are computed and reported for ALL cells, GATED only
   for U >= 6 (as in v13/v14); at U <= 3 they are reported-not-gated.
5. Totals envelope [Cox-Munk clean, Hu] with the 0.5x-2.0x clause: unchanged,
   all cells, gating.
6. K1 |dHs| (spectrum-integrated H_s vs the fetch law): FAIL if > 35% (the
   pre-registered threshold, unchanged) AND `hs_flag = true` when |dHs| > 15%
   (a forward tripwire: on this matrix it ANNOTATES, it does not re-gate).

Verdicts
--------
  low wind (U <= 3): UNRESOLVED unless a hard clause fires -- totals outside
     the 0.5x-2.0x envelope, negative / non-finite energy, a step in S(k),
     |dHs| > 35%, conformance (a) ratio > 1.5 -- and then FAIL.  (b) and (c)
     are NOT in the brief's low-wind hard list and therefore never drive that
     verdict; their values are reported and the status is printed.
  high wind (U = 6, 10): PASS requires the construction gates (positivity,
     smoothness, the A4 band partition) + totals + C1 + C1b + K1 +
     conformance (a) and (c) PASS + observation (b) PASS or flagged
     UNRESOLVED-with-reason.  Anything else is FAIL.
  Zero FAILs are required for the merge to stand (the brief's rule; the report
  applies it, this module only measures).
"""
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10
import v11_spectrum as V11
import v13_spectrum as V13
import run_v14_spectrum as V14

G = P8.G
C_M = V11.C_M
K_M = V11.K_M
USTAR_FIT_MIN = V13.USTAR_FIT_MIN

#: THE MERGE: the default source of this family.  jonswap stays selectable.
SPECTRUM_SOURCE_DEFAULT = "elfouhaily-form"
SPECTRUM_SOURCES = ("elfouhaily-form", "jonswap")

#: the merged default's alpha_m rule set (Arm D of chunk 10).
MERGE_ARM = "D"
MERGE_RULE = ("v13 as shipped, evaluated as the Arm-D rule set: alpha_m raw "
              "inside the ECKV fit domain (u* >= c_m/e = 0.085050 m/s), "
              "max(raw, floor) below it, with the floor solving "
              "sigma2_h,total = sigma2_h_target (default 0.003).  Provenance "
              "'eckv_fit_domain_raw' / 'low_wind_prior' as in v14.  No number "
              "moves relative to v14's Arm D / RAW spectrum.")
ARMS = V14.ARMS
D_STRICT = V14.D_STRICT
SIGMA2_H_TARGET_DEFAULT = V14.SIGMA2_H_TARGET_DEFAULT
SIGMA2_H_TARGETS = V14.SIGMA2_H_TARGETS

#: ---- the RETIRED clause (kept only so a reader can see it is not used) ---
RETIRED_B_ENV_MAX = V14.B_ENV_MAX
RETIRED_CLAUSE = {
    "clause": "max B(k) <= 0.012 over k in [10, 1500]",
    "retired": True,
    "reason": ("0.012 is the equilibrium-range coefficient beta (a MEAN, "
               "measured range 0.006-0.024), never a cap on the total "
               "saturation spectrum; Elfouhaily's own B_l(k_p) = "
               "0.1433 alpha_p gamma reproduces the three 'failing' v14 "
               "values exactly (0.0183 @ Omega 7.3, 0.0250 @ Omega 10.3, "
               "0.01249 = (1/2) alpha_m at U = 10)"),
    "where_it_still_appears": ("as an OBSERVATION only: b_max_observation() "
                               "prints max B and max B_l with gated=false"),
}

#: ---- CONFORMANCE constants (model vs its own published form) -------------
CONF_A_COEF = 0.1433
CONF_A_FLANK_FACTOR = 1.04
CONF_A_RATIO_MAX = 1.5
CONF_C_RANGE = (0.8, 1.25)

#: ---- OBSERVATION constants (independent anchors) -------------------------
BETA_I_P = 2.5
BETA_F_BAND = (1.3, 3.0)
BETA_BAND = (0.006, 0.024)
BETA_SLOPE = -4.0
BETA_CONVENTION = {
    "law": "S(omega) = beta * I(p) * g * u* * omega^-4",
    "omega": "angular frequency, rad/s (Toba's equilibrium-range form)",
    "I(p)": BETA_I_P,
    "fit": ("log-space least squares with the slope fixed at -4 over the "
            "intercept: log beta_hat = mean(log S(omega) + 4 log omega) - "
            "log(I(p) g u*)"),
    "band": "omega in 2 pi [1.3 f_p, 3 f_p], f_p = sqrt(g k_p) / 2 pi",
    "why_this_convention": ("the pre-registered band [0.006, 0.024] is the "
                            "literature range of THIS coefficient; the "
                            "identical fit written as E(f) with f in Hz "
                            "returns beta / (2 pi)^3 (248.05x smaller), which "
                            "the band cannot reach for any spectrum of this "
                            "normalisation -- both numbers are printed"),
    "literal_hz_reading": ("E(f) = beta I(p) g u* f^-4 with f in Hz is "
                           "printed as beta_Hz_convention in every cell"),
}

#: ---- carried, unchanged gates -------------------------------------------
HS_TOL = V14.HS_TOL                 # K1: |dHs| <= 35% (pre-registered)
HS_FLAG = 0.15                      # forward tripwire: annotates, no re-gate
TOTALS_TOL = V14.TOTALS_TOL
GRAVITY_SHARE_MAX = V14.GRAVITY_SHARE_MAX
CAP_SHARE_MIN = V14.CAP_SHARE_MIN
SIGMA2_H_MIN = V14.SIGMA2_H_MIN
BAND = V14.BAND
B_K_RANGE = V14.B_K_RANGE

PROVENANCE = V14.PROVENANCE

#: labels that must travel with every number and every output
WARN_OMEGA_C = 5.0


# --------------------------------------------------------------------------
# 1. assembly (identical code path to v14; only the labels are added)
# --------------------------------------------------------------------------
def assemble(U, F, h, arm=MERGE_ARM, hs_override=None, k_cut=None,
             gamma=None, source=SPECTRUM_SOURCE_DEFAULT,
             sigma2_h_target=SIGMA2_H_TARGET_DEFAULT):
    """`V14.assemble` with this family's default source and Arm-D rules.

    The seasonal numbers are produced by the SAME code path as v14 (v14's
    `assemble` is called, not copied), so `arm="RAW"` here is bit-identical to
    `arm="RAW"` in v14 by construction -- asserted per cell in the regression
    stage, not assumed.  `source="jonswap"` is returned by v14/v11 untouched
    (the rollback lane).
    """
    if source not in SPECTRUM_SOURCES:
        raise ValueError("unknown spectrum_source %r" % (source,))
    ctx = V14.assemble(U, F, h, arm=arm, hs_override=hs_override,
                       k_cut=k_cut, gamma=gamma, source=source,
                       sigma2_h_target=sigma2_h_target)
    if source != "elfouhaily-form":
        return ctx
    ctx["version"] = "3.10+v15"
    ctx["spectrum_source_default"] = SPECTRUM_SOURCE_DEFAULT
    ctx["merge_arm"] = MERGE_ARM
    ctx["label"] = (ctx["label"] + " -- v3.10 merge: elfouhaily-form is the "
                                    "default source, Arm-D rule set")
    prov = provenance(ctx)
    ctx["v15_provenance"] = prov
    ctx["omega_c"] = prov["omega_c"]
    ctx["warn_extrapolated_wave_age"] = prov["warn_extrapolated_wave_age"]
    return ctx


def provenance(ctx):
    """The provenance block every per-cell record and every output carries."""
    om = float(ctx["omega_c"])
    u = float(ctx["laws"]["U10_ms"])
    low = bool(u <= 3.0 + 1e-12)
    return {
        "omega_c": om,
        "warn_extrapolated_wave_age": bool(om > WARN_OMEGA_C),
        "warn_rule": ("Omega_c > 5: the gamma log-formula is continued past "
                      "its fitted range; the label says so and the warning "
                      "travels with every output"),
        "eckv_domain": ctx["eckv_domain"],
        "provenance": ctx["provenance"],
        "alpha_m_domain": ctx["alpha_m_domain"],
        "low_wind_cell": low,
        "low_wind_status": ("UNRESOLVED_AT_BEST (U <= 3: low-wind cells are "
                            "never a soft pass)" if low
                            else "gated (U >= 6 cells must PASS)"),
        "spectrum_source": ctx["spectrum_source"],
        "spectrum_source_default": SPECTRUM_SOURCE_DEFAULT,
        "merge_arm": MERGE_ARM,
    }


def s_parts(k, ctx):
    return V11.s_parts(k, ctx)


def cell_id(U, F):
    """The matrix cell label used by every v-family numbers file."""
    return "U%s_F%s" % (str(float(U)).replace(".", "p"),
                        str(float(F)).replace(".", "p"))


# --------------------------------------------------------------------------
# 2. the retired clause, as an OBSERVATION
# --------------------------------------------------------------------------
def b_max_observation(ctx, k_lo=B_K_RANGE[0], k_hi=B_K_RANGE[1], n=4000):
    """max B and max B_l -- printed, never gated (the clause is retired)."""
    k = np.logspace(math.log10(float(k_lo)), math.log10(float(k_hi)), int(n))
    s = V11.s_parts(k, ctx)
    B = s["B_l"] + s["B_h"]
    i = int(np.argmax(B))
    kb = np.logspace(math.log10(BAND[0]), math.log10(BAND[1]), 3000)
    sb = V11.s_parts(kb, ctx)
    j = int(np.argmax(sb["B_l"]))
    return {
        "k_range_rad_per_m": [float(k_lo), float(k_hi)],
        "B_max": float(B[i]), "k_at_B_max": float(k[i]),
        "B_l_max_50_200": float(sb["B_l"][j]),
        "k_at_B_l_max_50_200": float(kb[j]),
        "B_total_max_50_200": float(np.max(sb["B_l"] + sb["B_h"])),
        "retired_clause_value": RETIRED_B_ENV_MAX,
        "gated": False,
        "note": ("the 0.012 clause is RETIRED (see RETIRED_CLAUSE): this is an "
                 "observation, so nothing in this record turns B_max into a "
                 "verdict.  It is printed because a reader must be able to see "
                 "the number the retired clause used to score."),
    }


# --------------------------------------------------------------------------
# 3. CONFORMANCE -- the model against its own published form
# --------------------------------------------------------------------------
def conformance_gravity(ctx, n=6000):
    """(a) gravity flank: max B_l over [k_p/2, 2 k_p] vs 0.1433 a_p g 1.04."""
    kp = float(ctx["k_p_rad_per_m"])
    ap = float(ctx["alpha_p"])
    gam = float(ctx["gamma_bl"])
    lo, hi = float(kp) / 2.0, 2.0 * float(kp)
    k = np.logspace(math.log10(lo), math.log10(hi), int(n))
    s = V11.s_parts(k, ctx)
    bl = s["B_l"]
    i = int(np.argmax(bl))
    ref = CONF_A_COEF * ap * gam * CONF_A_FLANK_FACTOR
    ratio = float(bl[i] / ref) if ref else None
    # the peak value itself, for the report's cross-check with the literature
    kpk = np.array([kp])
    bl_at_kp = float(V11.s_parts(kpk, ctx)["B_l"][0])
    return {
        "what": ("max B_l over the gravity flank [k_p/2, 2 k_p] against "
                 "Elfouhaily's own B_l(k_p) = 0.1433 alpha_p gamma, times the "
                 "flank factor 1.04"),
        "band_rad_per_m": [lo, hi], "n": int(n),
        "k_p_rad_per_m": kp, "alpha_p": ap, "gamma_omega_c": gam,
        "B_l_flank_max": float(bl[i]), "k_at_flank_max": float(k[i]),
        "B_l_at_k_p": bl_at_kp,
        "reference": float(ref),
        "reference_form": "0.1433 * alpha_p * gamma(Omega_c) * 1.04",
        "ratio": ratio,
        "ratio_max": CONF_A_RATIO_MAX,
        "pass": (None if ratio is None else bool(ratio <= CONF_A_RATIO_MAX)),
        "column": "CONFORMANCE",
    }


def conformance_capillary(ctx, n=6000):
    """(c) capillary: max B_h over [k_m/2, 2 k_m] vs (1/2) alpha_m(u*)."""
    lo, hi = float(K_M) / 2.0, 2.0 * float(K_M)
    k = np.logspace(math.log10(lo), math.log10(hi), int(n))
    s = V11.s_parts(k, ctx)
    bh = s["B_h"]
    i = int(np.argmax(bh))
    am = float(ctx["alpha_m_value"])
    ref = 0.5 * am
    ratio = float(bh[i] / ref) if ref else None
    ok = (None if ratio is None else bool(
        CONF_C_RANGE[0] <= ratio <= CONF_C_RANGE[1]))
    return {
        "what": ("max B_h over the capillary band [k_m/2, 2 k_m] against "
                 "Elfouhaily's own B_h(k_m) = (1/2) alpha_m"),
        "band_rad_per_m": [lo, hi], "n": int(n),
        "k_m_rad_per_m": float(K_M),
        "B_h_band_max": float(bh[i]), "k_at_band_max": float(k[i]),
        "alpha_m_used": am, "reference": float(ref),
        "reference_form": "(1/2) * alpha_m(u*)",
        "ratio": ratio, "ratio_range": list(CONF_C_RANGE),
        "pass": ok,
        "column": "CONFORMANCE",
    }


# --------------------------------------------------------------------------
# 4. OBSERVATION -- the equilibrium-range coefficient
# --------------------------------------------------------------------------
def observation_beta(ctx, n=20000):
    """(b) the equilibrium-range beta, Toba convention (see BETA_CONVENTION).

    Reconstruction (the brief's): omega = sqrt(g k) (deep water), f = omega/2pi,
    E(f) = S(k) dk/df, dk/df = 2 (2 pi)^2 f / g.  The gate quantity is the
    coefficient of the omega-convention law; the Hz-convention number is the
    same fit divided by (2 pi)^3 and is printed beside it.
    """
    kp = float(ctx["k_p_rad_per_m"])
    us = float(ctx["ustar"]["ustar_ms"])
    fp = math.sqrt(G * kp) / (2.0 * math.pi)      # deep-water peak frequency
    khi = float(P8.viscous_cutoff(float(ctx["depth_m"])))
    k = np.logspace(math.log10(max(kp / 60.0, 1.0e-4)), math.log10(khi),
                    int(n))
    s = V11.s_parts(k, ctx)
    S = s["S"]
    om = np.sqrt(G * k)                            # deep water, per the brief
    f = om / (2.0 * math.pi)
    E = S * 2.0 * (2.0 * math.pi) ** 2 * f / G     # m^2 per Hz
    Sw = S * 2.0 * om / G                          # m^2 s per rad
    lo, hi = BETA_F_BAND[0] * fp, BETA_F_BAND[1] * fp
    m = (f >= lo) & (f <= hi) & np.isfinite(E) & (E > 0.0)
    if not m.any():
        return {"status": "FAIL_OUT_OF_BAND", "pass": False,
                "reason": "the fit band is empty on this cell",
                "column": "OBSERVATION"}
    fm, omm, Em, Swm = f[m], om[m], E[m], Sw[m]
    den = BETA_I_P * G * us
    c_om = math.exp(float(np.mean(np.log(Swm) + 4.0 * np.log(omm))))
    c_hz = math.exp(float(np.mean(np.log(Em) + 4.0 * np.log(fm))))
    beta_om = float(c_om / den) if den else None
    beta_hz = float(c_hz / den) if den else None
    resid = float(np.std(np.log(Swm) + 4.0 * np.log(omm)))
    in_band = (None if beta_om is None
               else bool(BETA_BAND[0] <= beta_om <= BETA_BAND[1]))
    # contamination of the FIT BAND: does the capillary part exceed gravity?
    Bl, Bh = s["B_l"][m], s["B_h"][m]
    excess = Bh - Bl
    cont = bool(np.any(excess > 0.0))
    j = int(np.argmax(excess))
    if cont:
        status, pass_flag = "UNRESOLVED_BAND_CONTAMINATED", None
        reason = ("B_h exceeds B_l inside the fit band at "
                  "k = %.2f rad/m (f = %.2f Hz, B_h - B_l = %+.3e): the fitted "
                  "beta is band-contaminated, so it is reported as UNRESOLVED, "
                  "never a silent PASS/FAIL"
                  % (float(k[m][j]), float(fm[j]), float(excess[j])))
    elif in_band:
        status, pass_flag = "PASS", True
        reason = "beta inside [0.006, 0.024] and the fit band is uncontaminated"
    else:
        status, pass_flag = "FAIL_OUT_OF_BAND", False
        reason = ("beta outside [0.006, 0.024] on an uncontaminated fit band "
                  "(%.5f)" % beta_om)
    return {
        "what": ("equilibrium-range coefficient from the assembled "
                 "omnidirectional spectrum: log-space LSQ of "
                 "S(omega) = beta I(p) g u* omega^-4 over f in [1.3 f_p, "
                 "3 f_p]"),
        "column": "OBSERVATION",
        "convention": dict(BETA_CONVENTION),
        "I_p": BETA_I_P, "u_star_ms": us, "g": G,
        "f_p_Hz_deep_water": float(fp),
        "f_p_Hz_fetch_law": float(ctx["laws"]["f_p_Hz"]),
        "f_p_relative_difference": float(fp / float(ctx["laws"]["f_p_Hz"]) - 1.0),
        "fit_band_Hz": [float(lo), float(hi)],
        "fit_band_rad_per_s": [float(2 * math.pi * lo),
                               float(2 * math.pi * hi)],
        "fit_band_rad_per_m": [float(k[m].min()), float(k[m].max())],
        "n_points_in_band": int(m.sum()),
        "beta": beta_om,
        "beta_Hz_convention": beta_hz,
        "beta_ratio_omega_over_hz": (None if not beta_hz else
                                     float(beta_om / beta_hz)),
        "beta_band": list(BETA_BAND),
        "beta_in_band": in_band,
        "fit_residual_log_rms": resid,
        "band_contaminated": cont,
        "contamination": {"n_points_bh_gt_bl": int(np.sum(excess > 0.0)),
                          "max_bh_minus_bl": float(float(excess[j])),
                          "k_at_max_rad_per_m": float(k[m][j]),
                          "f_at_max_Hz": float(fm[j])},
        "status": status, "pass": pass_flag, "reason": reason,
        "note": ("the gate quantity is beta in the stated Toba/omega "
                 "convention; beta_Hz_convention is the identical fit written "
                 "as E(f) with f in Hz and is exactly (2 pi)^3 smaller"),
    }


# --------------------------------------------------------------------------
# 5. OBSERVATION -- the JONSWAP anchor (context only)
# --------------------------------------------------------------------------
def jonswap_anchor(U, F, ctx):
    """(d) B_J(k_p) from JONSWAP's alpha_J, beside the model's B_l(k_p)."""
    U = float(U)
    F = float(F)
    gF_over_U2 = G * F / max(U, 1e-9) ** 2
    a_j = 0.076 * gF_over_U2 ** -0.22
    b_j = 0.5 * a_j * math.exp(-1.25) * 3.3
    kp = float(ctx["k_p_rad_per_m"])
    b_l_kp = float(V11.s_parts(np.array([kp]), ctx)["B_l"][0])
    sub_km = bool(F < 1000.0)
    return {
        "what": ("JONSWAP's own peak curvature anchor, printed beside the "
                 "model -- CONTEXT, never a gate"),
        "alpha_J": float(a_j),
        "alpha_J_form": "0.076 * (g F / U^2)^-0.22",
        "gF_over_U2": float(gF_over_U2),
        "B_J_at_k_p": float(b_j),
        "B_J_form": "0.5 * alpha_J * exp(-1.25) * 3.3",
        "model_B_l_at_k_p": b_l_kp,
        "ratio_model_over_anchor": (float(b_l_kp / b_j) if b_j else None),
        "sub_km_cell": sub_km,
        "sub_km_caveat": (
            "F = %.0f m is a SUB-KILOMETRE fetch: JONSWAP's alpha_J fetch "
            "relation was fitted over kilometres (gF/U^2 ~ 1e3-1e5), so on "
            "this cell the anchor is an extrapolation 1-3 decades below its "
            "fitted range.  It is printed for context only; no gate uses it."
            % F) if sub_km else None,
        "gated": False,
    }


# --------------------------------------------------------------------------
# 6. K1 -- |dHs| with the forward tripwire
# --------------------------------------------------------------------------
def hs_node(ctx):
    """K1: spectrum-integrated H_s vs the fetch law, plus the hs_flag."""
    node = V14.hs_node(ctx)
    d = node["delta_rel"]
    node["hs_flag"] = (None if d is None else bool(abs(d) > HS_FLAG))
    node["hs_flag_threshold"] = HS_FLAG
    node["hs_flag_note"] = ("forward tripwire adopted from the rally: it "
                            "ANNOTATES a cell (|dHs| > 15%); it does not "
                            "re-gate anything on this matrix -- the gate is "
                            "still |dHs| <= 35%")
    return node


# --------------------------------------------------------------------------
# 7. one cell, both columns, all flags, one verdict
# --------------------------------------------------------------------------
def evaluate_cell(U, F, h, k_cut=None, sigma2_h_target=None, ctx=None):
    """The v15 record for one cell: CONFORMANCE and OBSERVATION, never merged."""
    tgt = (SIGMA2_H_TARGET_DEFAULT if sigma2_h_target is None
           else float(sigma2_h_target))
    if ctx is None:
        ctx = assemble(float(U), float(F), float(h), k_cut=k_cut,
                       sigma2_h_target=tgt)
    geo = V14.cell_geometry(ctx)
    kmin = geo["kmin_rad_per_m"]
    kvisc = geo["k_viscous_rad_per_m"]
    kc = float(ctx["k_cut_rad_per_m"])
    kf = V11._k_grid(kvisc, k_hi=kvisc)
    an_full = V11.integral_block(kf, V11.s_parts(kf, ctx), kmin, kvisc,
                                 float(ctx["k_p_rad_per_m"]))
    kcut = V11._k_grid(kc)
    an_cut = V11.integral_block(kcut, V11.s_parts(kcut, ctx), kmin, kc,
                                float(ctx["k_p_rad_per_m"]))
    mom = V14.moments(ctx)
    sh = V14.shares(ctx, an_cut)
    env = V14.envelope(U, F, h, kmin, kc)
    ec = V14.envelope_check(an_cut["total"]["sigma2"], env)
    obs_b = observation_beta(ctx)
    conf_a = conformance_gravity(ctx)
    conf_c = conformance_capillary(ctx)
    anchor = jonswap_anchor(U, F, ctx)
    bobs = b_max_observation(ctx)
    sm = V14.smoothness(ctx)
    pos = V14.positivity(ctx)
    hs = hs_node(ctx)
    c1 = V14.gate_c1(sh)
    c1b = V14.gate_c1b(sh, mom)
    br = V14.band_reproduction(ctx, k_cut=kc)
    U_f = float(U)
    low_wind = bool(U_f <= 3.0 + 1e-12)
    construction = {
        "positivity": pos["pass"],
        "smoothness_no_step": sm["pass"],
        "a4_band_partition": br["pass"],
        "pass": bool(pos["pass"] and sm["pass"] and br["pass"]),
    }

    # ---- verdict --------------------------------------------------------
    hard = []
    if not ec["tol_inside"]:
        hard.append("totals outside [0.5x, 2.0x] of the envelope")
    if not pos["pass"]:
        hard.append("negative or non-finite energy")
    if not sm["pass"]:
        hard.append("a step in S(k) / B(k)")
    if hs["pass"] is False:
        hard.append("|dHs| > 35%")
    if conf_a["pass"] is False:
        hard.append("conformance (a) ratio > 1.5")
    fails = []
    if low_wind:
        fails = list(hard)
        state = "FAIL" if fails else "UNRESOLVED"
        rule = ("low-wind cell (U <= 3): UNRESOLVED unless a hard clause fires "
                "(totals outside 0.5x-2.0x, negative energy, a step in S(k), "
                "|dHs| > 35%, conformance (a) ratio > 1.5) -- then FAIL.  "
                "(b) and (c) are not in this list and never drive the verdict. "
                "UNRESOLVED is not a soft pass.")
    else:
        checks = {
            "construction": construction["pass"],
            "totals_tol": ec["tol_inside"],
            "C1": c1["pass"],
            "C1b": c1b["pass"],
            "K1_dHs": hs["pass"],
            "conformance_a": conf_a["pass"],
            "conformance_c": conf_c["pass"],
            "observation_b": (obs_b["pass"] is not False),
        }
        fails = [k for k, v in checks.items() if v is not True]
        state = "PASS" if not fails else "FAIL"
        rule = ("high-wind cell (U >= 6): PASS requires construction gates + "
                "totals inside 0.5x-2.0x + C1 + C1b + K1 + conformance (a),(c) "
                "PASS + observation (b) PASS or flagged UNRESOLVED-with-reason. "
                "Failing clauses: %s" % (", ".join(fails) or "none"))
    prov = provenance(ctx)
    return {
        "arm": MERGE_ARM,
        "cell": {"U_ms": U_f, "fetch_m": float(F), "depth_m": float(h)},
        "cell_id": cell_id(U_f, F),
        "spectrum_source": ctx["spectrum_source"],
        "low_wind_cell": low_wind,
        "provenance": prov,
        "omega_c": prov["omega_c"],
        "warn_extrapolated_wave_age": prov["warn_extrapolated_wave_age"],
        "geometry": geo,
        "sea_state": {k: ctx["sea_state"][k] for k in (
            "k_p_rad_per_m", "c_p_ms", "omega_c", "alpha_p", "gamma_bl",
            "sigma_bl", "eckv_domain", "k_p_dispersion_rel_diff",
            "kp_times_h")},
        "alpha": ctx["alpha_m_resolution"],
        "moments": mom,
        "integrals_cut": an_cut,
        "integrals_full": an_full,
        "shares": sh,
        "b_max_observation": bobs,
        "envelope": env,
        "envelope_check": ec,
        "hs": hs,
        "C1": c1, "C1b": c1b,
        "C1_C1b_gated": bool(not low_wind),
        "construction": construction,
        "smoothness": sm,
        "positivity": pos,
        "band_partition": br,
        "conformance": {"a_gravity": conf_a, "c_capillary": conf_c},
        "observation": {"b_beta": obs_b, "d_jonswap_anchor": anchor},
        "state": state,
        "fail_reasons": fails,
        "hard_clauses_fired": hard,
        "verdict_rule": rule,
        "Hs_spectrum_m": float(ctx["Hs_m"]),
        "Hs_law_m": float(ctx["Hs_law_m"]),
        "sigma2_total_cut": float(an_cut["total"]["sigma2"]),
        "sigma2_total_full": float(an_full["total"]["sigma2"]),
    }


def export_curves(ctx, n=4000):
    return V14.export_curves(ctx, n=n)


def reproduce_from_export(curve, k_cut):
    return V14.reproduce_from_export(curve, k_cut)


def cell_geometry(ctx):
    return V14.cell_geometry(ctx)
