"""v11 spectrum -- version 3.8 (PIPELINE-V10 chunk 8): gravity-tail termination.

One convention only (PIPELINE-V10 section 2, inherited verbatim from V9):

    S(k)  [m^2 per rad/m]   omnidirectional wavenumber spectrum
    m0    = int_0^inf S(k) dk          (height variance)
    sigma^2 = int_0^inf k^2 S(k) dk    (slope variance)
    Psi(k,theta) = S(k) * Phi(k,theta) / k,   int Phi dtheta = 1
    B(k) := k^3 S(k)                    (curvature spectrum, first-class here)

The gravity part is REPLACED WHOLESALE (PIPELINE-V10 section 1):

    S(k) = k^-3 * ( B_l(k) + B_h(k) )

    B_l  Elfouhaily et al. (1997) low-frequency part, PIPELINE-V10 section 3,
         driven by THIS pipeline's fetch laws:
             k_p from the fetch-law f_p through the existing finite-depth
             gravity-capillary dispersion (v8_physics),  c_p = omega(k_p)/k_p,
             Omega_c = U10 / c_p.  ECKV's own fetch->Omega_c fit is NOT used.
    B_h  v3.7's Elfouhaily high-frequency part, imported unchanged
         (v10_spectrum.elfouhaily_bh: L_pm-carried, McDaniel 2001 form).

No join, no blend, no splice, no renormalisation of either part.

HONEST LABEL (rally D3): the result is NOT "the Elfouhaily model".  It is an
**Elfouhaily-form spectrum with fetch-derived sea-state parameters**.  Code,
numbers files and reports use that label.

Domain labels (PIPELINE-V10 section 3), per cell:
    eckv_domain = "published"              for Omega_c <= 5
                = "wave_age_extrapolation" for Omega_c  > 5
The gamma log-formula is continued past Omega_c = 5, nothing is clamped and
nothing diverges; the label says plainly that those cells are uncalibrated
extrapolation.  The v3.7 low-wind B_h floor keeps its "low_wind_floor" label
in the separate field alpha_m_domain.

Nothing in this module is tuned.  There are no free parameters.
"""
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10

G = P8.G
TF = P8.TF
OMEGA0 = P8.OMEGA0

K_M = V10.K_M                      # 369.6 rad/m  (g/capillary crossover)
C_M = V10.C_M                      # 0.2310 m/s   minimum phase speed
SQRT10 = math.sqrt(10.0)

SOURCES = ("jonswap", "elfouhaily-form")

#: PIPELINE-V10 section 5-C band table, k in rad/m.  The first lower edge is
#: k_p (per cell); everything below k_p is reported separately as "below_kp"
#: so the rows of any cell sum to the cell total (no hidden bucket).
BAND_EDGES_V11 = (10.0, 50.0, 200.0, 370.0, 1000.0)
BAND_NAMES_V11 = ("kp_10", "10_50", "50_200", "200_370", "370_1000",
                  "1000_inf")

#: gravity / capillary split of the assembled spectrum, NEVER applied to the
#: curvature B (B is reported per part on the full band)
PART_GRAVITY = "gravity"
PART_CAPILLARY = "capillary"


# --------------------------------------------------------------------------
# 1. fetch laws -> sea-state parameters for B_l
# --------------------------------------------------------------------------
def peak_wavenumber_from_dispersion(omega_p, h):
    """Invert the pipeline's finite-depth gravity-capillary dispersion.

    omega^2 = (g k + (gamma/rho) k^3) tanh(k h).  The fetch law gives f_p
    (= omega_p / 2 pi) in the FREQUENCY domain; the peak wavenumber of the
    pipeline's own dispersion is the root of omega(k) = omega_p.  Monotonic
    in k, so a log-bisection is exact to machine precision.
    """
    w = float(omega_p)
    lo, hi = 1.0e-6, 1.0e6
    flo = float(P8.dispersion_omega(np.array([lo]), h)[0]) - w
    fhi = float(P8.dispersion_omega(np.array([hi]), h)[0]) - w
    if flo > 0.0 or fhi < 0.0:
        raise ValueError("omega_p outside the dispersion's bracket")
    for _ in range(220):
        mid = math.sqrt(lo * hi)
        fm = float(P8.dispersion_omega(np.array([mid]), h)[0]) - w
        if fm < 0.0:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < 1.0e-15:
            break
    k = math.sqrt(lo * hi)
    c = w / k
    return {"k_p_rad_per_m": float(k), "c_p_ms": float(c),
            "omega_residual_rel": float(
                float(P8.dispersion_omega(np.array([k]), h)[0]) / w - 1.0),
            "method": ("root of omega(k; h) = omega_p on the pipeline "
                       "dispersion (v8_physics.dispersion_omega)"),
            }


def sea_state(U, F, h):
    """fetch laws -> (k_p, c_p, Omega_c, alpha_p, gamma, sigma, domain)."""
    laws = V10.fetch_laws(U, F, h)
    pk = peak_wavenumber_from_dispersion(laws["omega_p_rad_s"], h)
    k_p = pk["k_p_rad_per_m"]
    c_p = pk["c_p_ms"]
    om_c = float(U) / c_p
    alpha_p = 6.0e-3 * om_c ** 0.55
    if om_c <= 1.0:
        gamma = 1.7
        gamma_form = "1.7"
    else:
        gamma = 1.7 + 6.0 * math.log10(om_c)
        gamma_form = "1.7 + 6 log10(Omega_c)"
    sigma = 0.08 * (1.0 + 4.0 * om_c ** -3.0)
    domain = ("published" if om_c <= 5.0 else "wave_age_extrapolation")
    return {
        "laws": laws,
        "k_p_rad_per_m": float(k_p), "c_p_ms": float(c_p),
        "omega_c": float(om_c),
        "alpha_p": float(alpha_p), "gamma_bl": float(gamma),
        "gamma_form": gamma_form, "sigma_bl": float(sigma),
        "eckv_domain": domain,
        "eckv_domain_rule": ("published for Omega_c <= 5; "
                             "wave_age_extrapolation for Omega_c > 5 "
                             "(gamma log-form continued, nothing clamped)"),
        "k_p_deep_water_rad_per_m": float(laws["k_p_rad_per_m"]),
        "k_p_dispersion_rel_diff": float(k_p / laws["k_p_rad_per_m"] - 1.0),
        "omega_residual_rel": pk["omega_residual_rel"],
        "kp_times_h": float(k_p * float(h)),
        "note": ("k_p is the ROOT of the pipeline dispersion at the fetch-law "
                 "omega_p (PIPELINE-V10 section 3), not the deep-water "
                 "lambda_p inversion v3.7 used; c_p = omega_p/k_p from the "
                 "same dispersion; Omega_c = U10/c_p.  ECKV's own ocean-fetch "
                 "relation is not used."),
    }


# --------------------------------------------------------------------------
# 2. B_l (PIPELINE-V10 section 3, verbatim formula)
# --------------------------------------------------------------------------
def b_l(k, h, k_p, c_p, omega_c, alpha_p, gamma_bl, sigma_bl):
    """Elfouhaily et al. (1997) low-frequency curvature spectrum B_l(k)."""
    k = np.asarray(k, dtype=np.float64)
    kk = np.maximum(k, 1.0e-12)
    om = P8.dispersion_omega(kk, h)
    c = om / kk
    rp = np.sqrt(kk / float(k_p))
    l_pm = np.exp(-1.25 * (float(k_p) / kk) ** 2)
    gam_arg = -((rp - 1.0) ** 2) / (2.0 * float(sigma_bl) ** 2)
    j_p = np.power(float(gamma_bl), np.exp(gam_arg))
    expo = np.exp(-(float(omega_c) / SQRT10) * (rp - 1.0))
    f_p = l_pm * j_p * expo
    bl = 0.5 * float(alpha_p) * (float(c_p) / c) * f_p
    bl = np.where(np.isfinite(bl), bl, 0.0)
    return np.maximum(bl, 0.0)


def b_h(k, h, ustar, k_p, alpha_m):
    """v3.7's high-frequency part, imported -- not re-derived here."""
    return V10.elfouhaily_bh(k, h, ustar, k_p, alpha_m)


# --------------------------------------------------------------------------
# 3. assembly
# --------------------------------------------------------------------------
def assemble(U, F, h, hs_override=None, k_cut=None, gamma=None,
             source="elfouhaily-form"):
    """The spectrum of the selected source, parameterised and ready.

    source="jonswap"          -> v10's ctx, untouched (the legacy source).
    source="elfouhaily-form"  -> the v11 ctx below.
    """
    kcut = float(k_cut if k_cut is not None else V10.KCUT_DEFAULT)
    if source not in SOURCES:
        raise ValueError("unknown spectrum_source %r" % (source,))
    if source == "jonswap":
        ctx = V10.assemble(U, F, h, hs_override=hs_override, k_cut=kcut,
                           gamma=(V10.JONSWAP_GAMMA if gamma is None
                                  else gamma))
        ctx["spectrum_source"] = "jonswap"
        ctx["label"] = "v3.7 two-regime (JONSWAP tail + Elouhaily B_h)"
        return ctx

    ss = sea_state(U, F, h)
    us = V10.ustar_from_u10(U)
    amr = V10.alpha_m_elfouhaily(us["ustar_ms"])
    ctx = {
        "version": "3.8",
        "spectrum_source": "elfouhaily-form",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters (NOT 'the Elfouhaily model')"),
        "laws": ss["laws"],
        "sea_state": ss,
        "ustar": us,
        "alpha_m": amr,
        "alpha_m_value": float(amr["alpha_m"]),
        "alpha_m_domain": amr["eckv_domain"],
        "eckv_domain": ss["eckv_domain"],
        "k_p_rad_per_m": float(ss["k_p_rad_per_m"]),
        "c_p_ms": float(ss["c_p_ms"]),
        "omega_c": float(ss["omega_c"]),
        "alpha_p": float(ss["alpha_p"]),
        "gamma_bl": float(ss["gamma_bl"]),
        "sigma_bl": float(ss["sigma_bl"]),
        "depth_m": float(h),
        "k_cut_rad_per_m": kcut,
        "hs_override_m": (None if hs_override is None else float(hs_override)),
        "hs_source": "predicted",
        "gamma": float(V10.JONSWAP_GAMMA),
    }
    # H_s is a PREDICTION under this source (PIPELINE-V10 section 4)
    k = _k_grid(kcut, k_hi=float(P8.viscous_cutoff(float(h))))
    s = s_parts(k, ctx)
    kmin = 2.0 * math.pi / V10.patch_for(ss["laws"]["lambda_p_m"])["patch_m"]
    m = k >= kmin
    ctx["Hs_m"] = float(4.0 * math.sqrt(
        max(float(np.trapz(s["S"][m], k[m])), 0.0)))
    ctx["Hs_law_m"] = float(ss["laws"]["Hs_depth_limited_m"])
    ctx["Hs_pred_rel_law"] = float(ctx["Hs_m"] / ctx["Hs_law_m"] - 1.0)
    return ctx


def s_parts(k, ctx):
    """S_g/S_l/S_h/B and the total at k, per the ctx's source."""
    k = np.maximum(np.asarray(k, dtype=np.float64), 0.0)
    if ctx.get("spectrum_source", "jonswap") == "jonswap":
        p = V10.s_parts(k, ctx)
        return {"k": k, "S_g": p["S_g"], "S_l": np.zeros_like(p["S_g"]),
                "S_h": p["S_h"], "B_l": np.zeros_like(p["S_g"]),
                "B_h": p["B_h"], "S": p["S"], "omega": p["omega"]}
    h = float(ctx["depth_m"])
    om = P8.dispersion_omega(np.maximum(k, 1.0e-12), h)
    bl = b_l(k, h, ctx["k_p_rad_per_m"], ctx["c_p_ms"], ctx["omega_c"],
             ctx["alpha_p"], ctx["gamma_bl"], ctx["sigma_bl"])
    bh = b_h(k, h, ctx["ustar"]["ustar_ms"], ctx["k_p_rad_per_m"],
             ctx["alpha_m_value"])
    kk = np.maximum(k, 1.0e-12)
    s_l = np.where(k > 0.0, bl / kk ** 3, 0.0)
    s_h = np.where(k > 0.0, bh / kk ** 3, 0.0)
    s_l = np.where(np.isfinite(s_l), s_l, 0.0)
    s_h = np.where(np.isfinite(s_h), s_h, 0.0)
    return {"k": k, "S_g": np.zeros_like(s_l), "S_l": s_l, "S_h": s_h,
            "B_l": bl, "B_h": bh, "S": s_l + s_h, "omega": om}


def s_total(k, ctx):
    return s_parts(k, ctx)["S"]


def mode_density_2d(kx, ky, ctx, weight=None, k_cut=None, theta_w=0.0,
                    delta=None):
    """2-D spectral density E(kx,ky) [m^2 per (rad/m)^2].

    E = S(k) Phi(k,theta) / k.  Phi is the UNCHANGED ECKV unified function on
    the total spectrum (PIPELINE-V10 section 1).
    """
    if ctx.get("spectrum_source", "jonswap") == "jonswap":
        return V10.mode_density_2d(kx, ky, ctx, weight=weight, k_cut=k_cut,
                                   theta_w=theta_w, delta=delta)
    kx = np.asarray(kx, dtype=np.float64)
    ky = np.asarray(ky, dtype=np.float64)
    k = np.sqrt(kx * kx + ky * ky)
    kk = np.maximum(k, 1.0e-12)
    S = s_total(kk.ravel(), ctx).reshape(k.shape)
    if delta is None:
        d = V10.eckv_delta(kk.ravel(), float(ctx["depth_m"]),
                           ctx["ustar"]["ustar_ms"],
                           ctx["k_p_rad_per_m"]).reshape(k.shape)
    else:
        d = np.asarray(delta, dtype=np.float64).reshape(k.shape)
    th = np.arctan2(ky, kx) - float(theta_w)
    phi = V10.eckv_phi(th, d)
    E = S * phi / kk
    if k_cut is not None:
        E = np.where(k <= float(k_cut), E, 0.0)
    if weight is not None:
        E = E * np.asarray(weight, dtype=np.float64)
    E = np.where(np.isfinite(E), E, 0.0)
    E[k == 0.0] = 0.0
    return E


# --------------------------------------------------------------------------
# 4. analytic integrals + the C1 band table
# --------------------------------------------------------------------------
K_INT_MIN = 1.0e-3
N_INT = 48000


def _k_grid(k_cut, k_hi=None, n=N_INT):
    hi = float(k_cut if k_hi is None else k_hi)
    return np.logspace(math.log10(K_INT_MIN), math.log10(hi), int(n))


def band_table(k, s, k_p, kmin, k_hi):
    """PIPELINE-V10 section 5-C table, split by part, per band.

    Rows: below_kp, then the six design bands.  Every row carries m0 and
    sigma^2 for the gravity part (B_l), the capillary part (B_h) and the
    total, plus the shares.  The rows sum to the cell total over
    [kmin, k_hi] (no hidden bucket).
    """
    edges = [float(kmin), float(k_p)] + list(BAND_EDGES_V11) + [float(k_hi)]
    names = ["below_kp"] + list(BAND_NAMES_V11)
    rows = {}
    tot = {"m0": 0.0, "sigma2": 0.0}
    for name, lo, hi in zip(names, edges[:-1], edges[1:]):
        if hi <= lo:
            rows[name] = {"k": [lo, hi], "empty": True,
                          "S_g": {"m0": 0.0, "sigma2": 0.0},
                          "S_h": {"m0": 0.0, "sigma2": 0.0},
                          "total": {"m0": 0.0, "sigma2": 0.0},
                          "sigma2_share": 0.0, "m0_share": 0.0,
                          "B_l_max": None, "B_h_max": None}
            continue
        m = (k >= lo) & (k < hi)
        kk = k[m]
        row = {"k": [float(lo), float(hi)], "empty": False}
        for tag, sl in (("S_g", s["B_l"][m] / kk ** 3),
                        ("S_h", s["B_h"][m] / kk ** 3)):
            row[tag] = {"m0": float(np.trapz(sl, kk)),
                        "sigma2": float(np.trapz(kk ** 2 * sl, kk))}
        row["total"] = {
            "m0": row["S_g"]["m0"] + row["S_h"]["m0"],
            "sigma2": row["S_g"]["sigma2"] + row["S_h"]["sigma2"]}
        row["B_l_max"] = (float(np.max(s["B_l"][m])) if m.any() else None)
        row["B_h_max"] = (float(np.max(s["B_h"][m])) if m.any() else None)
        row["B_l_mean"] = (float(np.trapz(s["B_l"][m], np.log(kk))
                                 / max(np.log(hi / lo), 1e-30))
                           if m.any() else None)
        rows[name] = row
        tot["m0"] += row["total"]["m0"]
        tot["sigma2"] += row["total"]["sigma2"]
    for name in rows:
        r = rows[name]
        r["sigma2_share"] = (float(r["total"]["sigma2"] / tot["sigma2"])
                             if tot["sigma2"] > 0 else 0.0)
        r["m0_share"] = (float(r["total"]["m0"] / tot["m0"])
                         if tot["m0"] > 0 else 0.0)
        if r["sigma2_share"] and r["total"]["sigma2"] > 0:
            r["gravity_share_of_sigma2"] = float(
                r["S_g"]["sigma2"] / r["total"]["sigma2"])
        else:
            r["gravity_share_of_sigma2"] = None
    return {"rows": rows, "total": tot, "range": [float(kmin), float(k_hi)],
            "note": ("bands per PIPELINE-V10 section 5-C; 'below_kp' carries "
                     "everything from kmin to k_p so the rows are exhaustive; "
                     "per-part integrals are of B_l/k^3 and B_h/k^3")}


def integral_block(k, s, kmin, k_hi, k_p):
    """m0 / sigma^2 / band table / part shares / B envelopes."""
    S = s["S"]
    m = (k >= kmin) & (k <= k_hi)
    kk = k[m]
    out = {}
    for tag, arr in (("total", S[m]), ("gravity", s["B_l"][m] / kk ** 3),
                     ("capillary", s["B_h"][m] / kk ** 3)):
        out[tag] = {"m0": float(np.trapz(arr, kk)),
                    "sigma2": float(np.trapz(kk * kk * arr, kk))}
    out["bands"] = band_table(k, s, k_p, kmin, k_hi)
    sig = out["total"]["sigma2"]
    for kx in (100.0, 370.0, 1173.0):
        mh = m & (k >= kx)
        s2 = float(np.trapz(k[mh] ** 2 * S[mh], k[mh]))
        out["frac_sigma2_above_%d" % int(kx)] = (float(s2 / sig) if sig
                                                 else None)
    mo = m & (k >= 50.0) & (k <= 200.0)
    ko = k[mo]
    if ko.size:
        g2 = float(np.trapz(ko ** 2 * s["B_l"][mo] / ko ** 3, ko))
        h2 = float(np.trapz(ko ** 2 * s["B_h"][mo] / ko ** 3, ko))
        out["overlap_50_200"] = {
            "B_l_sigma2": g2, "B_h_sigma2": h2,
            "gravity_share_of_sigma2": float(g2 / max(g2 + h2, 1e-30)),
            "B_l_max": float(np.max(s["B_l"][mo])),
            "B_h_max": float(np.max(s["B_h"][mo])),
            "B_total_max": float(np.max(s["B_l"][mo] + s["B_h"][mo])),
            "B_l_mean": float(np.trapz(s["B_l"][mo], np.log(ko))
                              / math.log(200.0 / 50.0)),
            "note": ("PIPELINE-V10 5-C gate quantities: gravity share of the "
                     "50-200 slope variance (<= 0.50) and max B_gravity in "
                     "50-200 (<= 0.010)"),
        }
    else:
        out["overlap_50_200"] = None
    return out


def anisotropy_prediction(ctx, k, kmin, k_cut):
    """Same estimator as v10, on the v11 spectrum (used by the state)."""
    s = s_parts(k, ctx)
    m = (k >= kmin) & (k <= k_cut)
    kk = k[m]
    w = kk * kk * s["S"][m]
    d = V10.eckv_delta(kk, float(ctx["depth_m"]), ctx["ustar"]["ustar_ms"],
                       ctx["k_p_rad_per_m"])
    tot = float(np.trapz(w, kk))
    de = float(np.trapz(w * d, kk) / tot) if tot > 0 else 0.0
    sa = tot * 0.5 * (1.0 + de / 2.0)
    sc = tot * 0.5 * (1.0 - de / 2.0)
    return {"delta_effective": de, "sigma2_along": float(sa),
            "sigma2_cross": float(sc),
            "ratio_along_over_cross": (float(sa / sc) if sc > 0 else None)}


# --------------------------------------------------------------------------
# 5. the ECKV reference at MATCHED inputs (C3 / C2 comparison)
# --------------------------------------------------------------------------
def eckv_reference(k, h, k_p, c_p, omega_c, alpha_p, gamma_bl, sigma_bl,
                   ustar, alpha_m, deep_water=True):
    """ECKV-as-published shape at the SAME (k_p, Omega_c, u*).

    The only difference to the implemented source is the phase speed: ECKV as
    published uses deep-water c(k) = sqrt(g/k + (gamma/rho) k).  That is the
    stated deviation of PIPELINE-V10 section 3, not a free parameter.
    """
    k = np.maximum(np.asarray(k, dtype=np.float64), 1.0e-12)
    if deep_water:
        om = np.sqrt(G * k + TF * k ** 3)
    else:
        om = P8.dispersion_omega(k, h)
    c = om / k
    rp = np.sqrt(k / float(k_p))
    l_pm = np.exp(-1.25 * (float(k_p) / k) ** 2)
    j_p = np.power(float(gamma_bl),
                   np.exp(-((rp - 1.0) ** 2) / (2.0 * float(sigma_bl) ** 2)))
    bl = 0.5 * float(alpha_p) * (float(c_p) / c) * l_pm * j_p \
        * np.exp(-(float(omega_c) / SQRT10) * (rp - 1.0))
    fm = np.exp(-0.25 * (k / K_M - 1.0) ** 2)
    bh = 0.5 * float(alpha_m) * (C_M / c) * l_pm * fm
    bl = np.where(np.isfinite(bl), bl, 0.0)
    bh = np.where(np.isfinite(bh), bh, 0.0)
    b = np.maximum(bl, 0.0) + np.maximum(bh, 0.0)
    return {"k": k, "B": b, "B_l": bl, "B_h": bh, "c": c}


def eckv_reference_totals(ctx, k_hi=None):
    """Total sigma^2 / m0 / frac>370 of ECKV at matched inputs (no cut)."""
    khi = float(k_hi if k_hi is not None else P8.viscous_cutoff(
        float(ctx["depth_m"])))
    k = _k_grid(khi, k_hi=khi)
    ref = eckv_reference(k, float(ctx["depth_m"]), ctx["k_p_rad_per_m"],
                         ctx["c_p_ms"], ctx["omega_c"], ctx["alpha_p"],
                         ctx["gamma_bl"], ctx["sigma_bl"],
                         ctx["ustar"]["ustar_ms"], ctx["alpha_m_value"],
                         deep_water=True)
    S = ref["B"] / k ** 3
    m0 = float(np.trapz(S, k))
    s2 = float(np.trapz(k ** 2 * S, k))
    mh = k >= 370.0
    s2h = float(np.trapz(k[mh] ** 2 * S[mh], k[mh]))
    return {"k_hi": khi, "m0_m2": m0, "sigma2_rad2": s2,
            "Hs_m": float(4.0 * math.sqrt(max(m0, 0.0))),
            "frac_sigma2_above_370": (float(s2h / s2) if s2 else None),
            "phase_speed": "deep water (ECKV as published)",
            "note": ("ECKV shape at the SAME (k_p, Omega_c, u*) supplied by "
                     "this pipeline's fetch laws; no ECKV fetch relation.  "
                     "Compared with the implemented spectrum this isolates "
                     "(a) the finite-depth phase speed and (b) the k_cut.")}


# --------------------------------------------------------------------------
# 6. gate A1 -- is the v3.7 B_g really the deep alpha_J/2 plateau?
# --------------------------------------------------------------------------
def bg_closed_form(k, h, k_p, omega_p, alpha_j, k_m=None):
    """B_g ~ (alpha_J/2) (1+3k^2/k_m^2)/(1+k^2/k_m^2)^3 exp[-1.25 (w_p/w)^4].

    Exact for the pipeline's finite-depth gravity-capillary dispersion in the
    DEEP limit (tanh -> 1); the exp uses the pipeline's own omega(k; h).

    k_m defaults to the research's k_m = 2 pi / 0.017 m = 369.6 rad/m.  The
    pipeline's dispersion has its own crossover sqrt(g/(gamma/rho)) = 367.0
    rad/m; the difference is a 0.7% convention difference and both are
    reported (ratio_kM_published / ratio_kM_physical) rather than hidden.
    """
    k = np.asarray(k, dtype=np.float64)
    kap2 = (k / float(K_M if k_m is None else k_m)) ** 2
    om = P8.dispersion_omega(k, h)
    b = (0.5 * float(alpha_j) * (1.0 + 3.0 * kap2) / (1.0 + kap2) ** 3
         * np.exp(-1.25 * (float(omega_p) / om) ** 4))
    return np.where(np.isfinite(b), b, 0.0)


def gate_a1(h_list=(1.0e3, 1.0, 0.5, 0.3), U=6.0, F=100.0,
            k_band=(50.0, 200.0), n=2000):
    """A1 (D1 verification): numeric v3.7 B_g vs the research closed form."""
    k = np.logspace(math.log10(k_band[0]), math.log10(k_band[1]), int(n))
    rows = {}
    deep_plateau = None
    for h in h_list:
        ctx = V10.assemble(U, F, h)
        s = V10.s_parts(k, ctx)
        bnum = k ** 3 * s["S_g"]
        alpha_j = float(ctx["jonswap_norm"]["A"]) / (G * G)
        bcf = bg_closed_form(k, h, ctx["k_p_rad_per_m"],
                             ctx["omega_p_rad_s"], alpha_j)
        k_m_phys = math.sqrt(G / TF)
        bcf_phys = bg_closed_form(k, h, ctx["k_p_rad_per_m"],
                                  ctx["omega_p_rad_s"], alpha_j,
                                  k_m=k_m_phys)
        ratio = np.where(bcf > 0, bnum / np.maximum(bcf, 1e-300), np.nan)
        ratio_phys = np.where(bcf_phys > 0,
                              bnum / np.maximum(bcf_phys, 1e-300), np.nan)
        plateau = float(np.exp(np.mean(np.log(np.maximum(bnum, 1e-300)))))
        row = {
            "depth_m": float(h),
            "kh_over_band": [float(k[0] * h), float(k[-1] * h)],
            "alpha_J_from_A_over_g2": float(alpha_j),
            "alpha_J_published_law": float(ctx["laws"]["alpha"]),
            "B_g_numeric_mean": plateau,
            "B_g_numeric_min": float(np.min(bnum)),
            "B_g_numeric_max": float(np.max(bnum)),
            "ratio_numeric_over_closed": {
                "min": float(np.nanmin(ratio)), "max": float(np.nanmax(ratio)),
                "mean": float(np.nanmean(ratio)),
                "max_abs_minus_1": float(np.nanmax(np.abs(ratio - 1.0))),
                "k_m_used_rad_per_m": float(K_M),
                "k_m_source": "research closed form: k_m = 2 pi / 0.017 m",
            },
            "ratio_kM_published": {
                "min": float(np.nanmin(ratio)),
                "max": float(np.nanmax(ratio)),
                "max_abs_minus_1": float(np.nanmax(np.abs(ratio - 1.0))),
            },
            "ratio_kM_physical": {
                "min": float(np.nanmin(ratio_phys)),
                "max": float(np.nanmax(ratio_phys)),
                "max_abs_minus_1": float(np.nanmax(np.abs(ratio_phys - 1.0))),
                "k_m_used_rad_per_m": float(k_m_phys),
                "k_m_source": "sqrt(g/(gamma/rho)) of the pipeline dispersion",
            },
            "k_m_convention_rel_diff": float(k_m_phys / K_M - 1.0),
            "ratio_using_published_alpha": {
                "min": float(np.nanmin(ratio * ctx["laws"]["alpha"]
                                       / alpha_j)),
                "max": float(np.nanmax(ratio * ctx["laws"]["alpha"]
                                       / alpha_j)),
            },
            "k_at_max_dev": float(k[int(np.nanargmax(np.abs(ratio - 1.0)))]),
            "lambda_peak_m": float(ctx["laws"]["lambda_p_m"]),
            "k_p_rad_per_m": float(ctx["k_p_rad_per_m"]),
        }

        if deep_plateau is None:
            deep_plateau = plateau
        row["plateau_rel_vs_deep"] = float(plateau / deep_plateau - 1.0)
        rows["h_%g" % h] = row
    worst = 0.0
    for h, r in rows.items():
        if h != "h_%g" % h_list[0]:
            worst = max(worst, abs(r["plateau_rel_vs_deep"]))
    return {
        "gate": "A1 -- D1 verification (PIPELINE-V10 section 5-A)",
        "what": ("B_g(k) computed NUMERICALLY from the current v3.7 code path "
                 "(V10.s_parts: S_omega -> exact domega/dk transform -> "
                 "B = k^3 S) against the research closed form with "
                 "kappa = k/k_m, k_m = 2 pi/0.017 m"),
        "cell": {"U_ms": float(U), "fetch_m": float(F), "k_band": list(k_band),
                 "n_points": int(n)},
        "deep_water_convention": {
            "h_m": float(h_list[0]),
            "why": ("h = 1000 m puts kh >= 5e4 across the band, i.e. "
                    "tanh(kh) = 1 to machine precision: the deep-water "
                    "closed form is exact there and any residual is numerical"),
        },
        "depths": rows,
        "worst_finite_depth_plateau_rel": float(worst),
        "failure_condition": 0.30,
        "finding": None,
        "pass": None,
    }


# --------------------------------------------------------------------------
# 7. gate A2 -- Parseval on a built test field (v10 G0 analogue)
# --------------------------------------------------------------------------
def build_test_field(ctx, N, patch_m, k_cut=None, seed=12345, theta_w=0.0):
    """Single-grid realisation of the assembled spectrum (G0's test field)."""
    N = int(N)
    S_patch = float(patch_m)
    dx = S_patch / N
    fx = np.fft.fftfreq(N, d=dx)
    kx1 = 2.0 * np.pi * fx
    kx, ky = np.meshgrid(kx1, kx1)
    kcut = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    E = mode_density_2d(kx, ky, ctx, k_cut=kcut, theta_w=theta_w)
    k = np.sqrt(kx * kx + ky * ky)
    dk = 2.0 * np.pi / S_patch
    D = E * dk * dk
    rng = np.random.default_rng(int(seed))
    a = rng.standard_normal((N, N))
    b = rng.standard_normal((N, N))
    H0 = (np.sqrt(2.0 * N ** 4 * D) * (a + 1j * b) / np.sqrt(2.0)
          ).astype(np.complex64)
    H0[0, 0] = 0.0
    h = np.real(np.fft.ifft2(H0.astype(np.complex128)))
    gx = np.real(np.fft.ifft2(1j * kx * H0.astype(np.complex128)))
    gy = np.real(np.fft.ifft2(1j * ky * H0.astype(np.complex128)))
    m0_num = float(np.mean(h * h))
    sig_num = float(np.mean(gx * gx + gy * gy))
    m0_mode = float(D.sum())
    sig_mode = float((k * k * D).sum())
    dof_m0 = float(D.sum() ** 2 / max((D * D).sum(), 1e-300))
    Dg = k * k * D
    dof_sig = float(Dg.sum() ** 2 / max((Dg * Dg).sum(), 1e-300))
    return {"m0_numeric": m0_num, "sigma2_numeric": sig_num,
            "m0_mode_sum": m0_mode, "sigma2_mode_sum": sig_mode,
            "dof_m0": dof_m0, "dof_sigma2": dof_sig,
            "patch_m": S_patch, "N": N, "k_cut": kcut}


def gate_a2(ctx, N=1024, patch_m=None, k_cut=None, seed=12345, theta_w=0.0,
            tag="elfouhaily-form"):
    """A2: the built field's m0 / sigma^2 match the assembled S(k)."""
    pm = float(patch_m if patch_m is not None else
               V10.patch_for(ctx["laws"]["lambda_p_m"])["patch_m"])
    t = build_test_field(ctx, N, pm, k_cut=k_cut, seed=seed, theta_w=theta_w)
    kcut = float(t["k_cut"])
    kmin = 2.0 * math.pi / pm
    k = _k_grid(kcut)
    s = s_parts(k, ctx)
    m = k >= kmin
    m0_cont = float(np.trapz(s["S"][m], k[m]))
    sig_cont = float(np.trapz(k[m] ** 2 * s["S"][m], k[m]))

    def rel(a, b):
        return None if b == 0 else float(a / b - 1.0)

    def z(num, exp, dof):
        sd = exp * math.sqrt(2.0 / dof) if dof > 0 else 0.0
        return None if sd <= 0 else float((num - exp) / sd)

    out = {
        "gate": "A2 -- Parseval on a built field from the %s spectrum" % tag,
        "spectrum_source": tag,
        "convention": ("m0 = int S dk, sigma^2 = int k^2 S dk, "
                       "Psi = S*Phi/k with int Phi dtheta = 1"),
        "test_field": {"N": int(N), "patch_m": pm, "k_cut": kcut,
                       "kmin": kmin, "seed": int(seed),
                       "k_nyquist": float(math.pi * int(N) / pm)},
        "m0": {"numeric": t["m0_numeric"], "mode_sum": t["m0_mode_sum"],
               "continuum": m0_cont},
        "sigma2": {"numeric": t["sigma2_numeric"],
                   "mode_sum": t["sigma2_mode_sum"],
                   "continuum": sig_cont},
        "dof": {"m0": t["dof_m0"], "sigma2": t["dof_sigma2"]},
        "rel_numeric_over_mode_m0": rel(t["m0_numeric"], t["m0_mode_sum"]),
        "rel_numeric_over_mode_sigma2": rel(t["sigma2_numeric"],
                                            t["sigma2_mode_sum"]),
        "rel_mode_over_continuum_m0": rel(t["m0_mode_sum"], m0_cont),
        "rel_mode_over_continuum_sigma2": rel(t["sigma2_mode_sum"], sig_cont),
        "z_m0": z(t["m0_numeric"], t["m0_mode_sum"], t["dof_m0"]),
        "z_sigma2": z(t["sigma2_numeric"], t["sigma2_mode_sum"],
                      t["dof_sigma2"]),
        "tolerance": {"z_abs_max": 3.0, "discretisation_rel_max": 0.30},
    }
    ok_a = (abs(out["z_m0"] or 0.0) <= 3.0
            and abs(out["z_sigma2"] or 0.0) <= 3.0)
    ok_b = (abs(out["rel_mode_over_continuum_m0"] or 0.0) <= 0.30
            and abs(out["rel_mode_over_continuum_sigma2"] or 0.0) <= 0.30)
    out["pass_convention"] = bool(ok_a)
    out["pass_discretisation"] = bool(ok_b)
    out["pass"] = bool(ok_a and ok_b)
    return out


# --------------------------------------------------------------------------
# 8. gate K6 -- continuity audit of the assembled B(k)
# --------------------------------------------------------------------------
def _k6_pass(ctx, k_lo, k_hi, n, floor_rel):
    k = np.logspace(math.log10(k_lo), math.log10(k_hi), int(n))
    s = s_parts(k, ctx)
    B = s["B_l"] + s["B_h"]
    step = np.diff(np.log(k))
    ok = B >= float(floor_rel) * float(np.max(B))
    i0 = int(np.argmax(ok)) if ok.any() else 0
    kb, Bb = k[i0:], B[i0:]
    rel_step = np.diff(Bb) / np.maximum(Bb[:-1], 1e-300)
    slope = np.diff(np.log(np.maximum(Bb, 1e-300))) / step[i0:]
    dslope = np.abs(np.diff(slope))
    return {
        "n": int(n), "log_step": float(step[0]),
        "B_min": float(np.min(B)), "B_max": float(np.max(B)),
        "B_all_finite": bool(np.all(np.isfinite(B))),
        "B_nonneg": bool(np.min(B) >= 0.0),
        "resolvable_from_rad_per_m": float(k[i0]),
        "max_rel_step": float(np.max(np.abs(rel_step))),
        "k_at_max_rel_step": float(kb[int(np.argmax(np.abs(rel_step)))]),
        "max_dslope": float(np.max(dslope)),
        "k_at_max_dslope": float(kb[int(np.argmax(dslope)) + 1]),
        "slope_range": [float(np.min(slope)), float(np.max(slope))],
    }


def gate_k6(ctx, k_lo=1.0, k_hi=3000.0, n=4000, floor_rel=1e-12):
    """K6: B(k) and d ln B / d ln k on a log grid -- no jumps, B >= 0.

    The discontinuity test is a RESOLUTION test, not a magic threshold: the
    largest one-step change of the log-slope is measured at n and 2n points
    on the same band.  A smooth function (however sharp its peak) has
    |d slope| proportional to the step, so it HALVES when the step halves; a
    splice or a renormalisation step has a fixed absolute slope jump and does
    not.  Both numbers and the ratio are reported.

    Below floor_rel * max(B) the spectrum has underflowed (B_l ~
    exp(-1.25 (k_p/k)^2) is 1e-117 at k = 1 rad/m): B is still reported there
    (finite, non-negative) but excluded from the step tests, because a
    relative-step test in denormal territory measures rounding, not physics.
    """
    a = _k6_pass(ctx, k_lo, k_hi, n, floor_rel)
    b = _k6_pass(ctx, k_lo, k_hi, int(n) * 2, floor_rel)
    ratio = (float(b["max_dslope"] / a["max_dslope"])
             if a["max_dslope"] > 0 else None)
    step_bound = 0.50          # |dB/B| in one step (smooth here: <= 0.12)
    ratio_bound = 0.65         # smooth => 0.5 ; a fixed jump => 1.0
    return {
        "gate": "K6 -- continuity audit of the assembled B(k)",
        "k_range_rad_per_m": [float(k_lo), float(k_hi)],
        "resolvable_floor_rel": float(floor_rel),
        "coarse": a, "fine": b,
        "max_dslope_ratio_fine_over_coarse": ratio,
        "step_threshold_rel": float(step_bound),
        "no_step_jump": bool(a["max_rel_step"] <= step_bound
                             and b["max_rel_step"] <= step_bound),
        "slope_halving_ratio_bound": float(ratio_bound),
        "no_slope_discontinuity": bool(ratio is not None
                                       and ratio <= ratio_bound),
        "B_nonneg": bool(a["B_nonneg"] and b["B_nonneg"]),
        "domain_label_present": bool(
            ctx.get("eckv_domain") in
            ("published", "wave_age_extrapolation")),
        "eckv_domain": ctx.get("eckv_domain"),
        "alpha_m_domain": ctx.get("alpha_m_domain"),
        "B_min_over_grid": a["B_min"],
        "excluded_below_rad_per_m": [float(k_lo),
                                     a["resolvable_from_rad_per_m"]],
        "excluded_note": ("B there falls to %.1e: underflow, not a jump; "
                          "still finite and non-negative."
                          % a["B_min"]),
        "pass": bool(a["B_nonneg"] and a["B_all_finite"]
                     and a["max_rel_step"] <= step_bound
                     and b["max_rel_step"] <= step_bound
                     and ratio is not None and ratio <= ratio_bound
                     and ctx.get("eckv_domain") in
                     ("published", "wave_age_extrapolation")),
        "note": ("a splice would be a step of order 1 in B (>= 50% in one "
                 "0.2% step of ln k) or a slope jump that survives grid "
                 "refinement; neither is present"),
    }


# --------------------------------------------------------------------------
# 9. curves for the deliverable
# --------------------------------------------------------------------------
def curves(ctx, patch_m, k_cut=None, n=420):
    kcut = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    kmin = 2.0 * math.pi / float(patch_m)
    k = np.logspace(math.log10(kmin), math.log10(kcut), int(n))
    s = s_parts(k, ctx)
    d = V10.eckv_delta(k, float(ctx["depth_m"]), ctx["ustar"]["ustar_ms"],
                       ctx["k_p_rad_per_m"])
    return {"patch_m": float(patch_m),
            "spectrum_source": ctx.get("spectrum_source", "jonswap"),
            "k_rad_per_m": [float(x) for x in k],
            "S_l_m2_per_rad_m": [float(x) for x in s["S_l"]],
            "S_g_m2_per_rad_m": [float(x) for x in s["S_g"]],
            "S_h_m2_per_rad_m": [float(x) for x in s["S_h"]],
            "S_total_m2_per_rad_m": [float(x) for x in s["S"]],
            "B_l": [float(x) for x in s["B_l"]],
            "B_h": [float(x) for x in s["B_h"]],
            "B_total": [float(x) for x in (s["B_l"] + s["B_h"])],
            "delta_eckv": [float(x) for x in d],
            "omega_rad_s": [float(x) for x in s["omega"]]}
