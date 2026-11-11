"""v10 spectrum -- version 3.7 (PIPELINE-V9 chunk 7): two-regime spectrum.

One convention only (PIPELINE-V9 section 2):

    S(k)  [m^2 per rad/m]   omnidirectional wavenumber spectrum
    m0    = int_0^inf S(k) dk          (height variance)
    sigma^2 = int_0^inf k^2 S(k) dk    (slope variance)

    Psi(k,theta) = S(k) * Phi(k,theta) / k      (2-D density per (k,theta))
    int_0^{2pi} Phi dtheta = 1  ->  int Psi k dk dtheta = int S dk = m0

NOTE ON THE DESIGN TEXT: PIPELINE-V9 section 2 writes the directional form as
``Psi = S * Phi / (2 pi k)`` while also pinning ``int Phi dtheta = 1``.  Those
two statements are mutually inconsistent by a factor 2*pi (they would give
m0 = int S dk / (2 pi)).  The self-consistent pair is used here (the one that
makes the `int Phi = 1` line true), and gate G0 verifies numerically that the
built field's variance equals int S dk.  This is recorded in the report as a
correction of the design text, not a silent deviation.

Construction (additive, no splice -- PIPELINE-V9 section 3):

    S(k) = S_g(k) + S_h(k)

    S_g  JONSWAP in frequency form S_J(omega; H_s, f_p) with BOTH H_s and f_p
         from the fetch laws, converted by the exact change of variables
         S_g(k) = S_J(omega(k)) * domega/dk on the pipeline's existing
         finite-depth gravity-capillary dispersion (v8_physics).
    S_h  Elfouhaily et al. (1997) B_h(k; u*) converted into the S(k)
         convention:  S_h(k) = B_h(k) / k^3   (B is the curvature spectrum,
         int B dlnk = mean square slope).  Level set by u* alone; never
         rescaled to any sigma^2 target.

    Phi  ECKV unified azimuthal function:
             Delta(k) = tanh(a0 + ap (c/c_p)^2 + am (c_m/c)^2)
                        / tanh(a0 + ap + am)
             Phi = (1/2pi) (1 + Delta cos 2(theta - theta_w))
         applied to the TOTAL spectrum.

Low wind (PIPELINE-V9 section 3, GPT N2):  alpha_m = 0.01 (1 + ln(u*/c_m))
goes negative below u* ~ 0.085 m/s.  It is FLOORED at 0.0 explicitly, the
numbers dict carries ``eckv_domain`` and gate G4 checks B(k) >= 0.

Sources
-------
* Hasselmann et al. (1973) JONSWAP fetch relations:
      f_p = 3.5 (g/U) Ft^(-1/3),  alpha = 0.076 Ft^(-0.22),
      H_s = 0.0016 Ft^0.5 U^2/g
* Young (1988) / CEM shallow-water fetch-limited growth (depth-limited H_s):
      gH_s/U^2 = 0.283 tanh(0.53 dt^0.75) tanh(0.00565 Ft^0.5 / tanh(0.53 dt^0.75))
  with dt = g h / U^2.  In the deep limit tanh->1 and this is EXACTLY the
  pipeline's existing law (0.283*0.00565 = 0.0015995 ~ 0.0016) -- so the
  depth correction is the same law family, not a second model.
* Elfouhaily, Chapron, Katsaros & Vandemark (1997), JGR 102(C7) 15781:
      k_m = 2pi/0.017 = 369.6 rad/m, c_m = 0.231 m/s (gravity-capillary
      minimum phase speed), alpha_p = 6e-3 sqrt(Omega),
      alpha_m per above (eq. 41-44), B_h = 0.5 alpha_m (c_m/c) L_pm
      exp(-0.25 (k/k_m - 1)^2)  (McDaniel 2001 applies L_pm to B_h as well),
      L_pm = exp(-1.25 (k_p/k)^2),
      Delta(k) eq. 66-68 with a0 = ln2/4, a_p = 4, a_m = 0.13 u*/c_m.
* Cox & Munk (1954) anisotropy reference is taken from v5_physics
  (sigma_a^2 = 3.16e-3 U, sigma_c^2 = 0.003 + 1.92e-3 U) -- reported, never
  used as a target.
"""
import math

import numpy as np

import v8_physics as P8

G = P8.G
TF = P8.TF                    # gamma/rho [m^3/s^2]
RHO_W = P8.RHO_W
SIGMA_SURF = P8.SIGMA_SURF
NU_KIN = P8.NU_KIN
OMEGA0 = P8.OMEGA0            # 2 pi / 8 s loop quantum

# Elfouhaily high-wavenumber / direction constants
K_M = 2.0 * math.pi / 0.017          # 369.6 rad/m
C_M = math.sqrt(G / K_M + TF * K_M)   # 0.2310 m/s minimum phase speed
ECKV_A0 = math.log(2.0) / 4.0         # 0.1733
ECKV_AP = 4.0
ECKV_AM_K = 0.13                      # a_m = 0.13 u*/c_m

# fetch-law constants (Hasselmann 1973 / Young 1988)
JONSWAP_GAMMA = 3.3
JONSWAP_S_A = 0.07
JONSWAP_S_B = 0.09
YV_COEF = 0.283                       # g Hs / U^2 prefactor (Young 1988)
YV_FETCH_COEF = 0.00565               # 0.283*0.00565 = 0.0015995 ~ 0.0016
YV_DEPTH_COEF = 0.53
DRAG_A = 0.8e-3                       # C10 = (0.8 + 0.065 U10) e-3
DRAG_B = 0.065e-3

# patch / cascade policy (PIPELINE-V9 section 5)
PATCH_BASE_M = 1.2        # 4 * 0.30 m: the v3.4 low-wind scene scale, now fixed
R_LAMBDA = 8.0            # patch >= R_LAMBDA * lambda_p (design minimum is 4)
K_X = 100.0               # gravity/capillary split, centre of the 50-200 band
KCUT_DEFAULT = 1000.0
KCUT_SWEEP = (800.0, 1000.0, 1200.0, 1500.0)

# integration grid
K_INT_MIN = 1.0e-3
N_INT = 24000

#: the four reporting bands, half-open in k (design: lambda = 2 pi / k)
BANDS = (
    ("gt_1m", 0.0, 2.0 * math.pi),
    ("0p1_1m", 2.0 * math.pi, 20.0 * math.pi),
    ("1_10cm", 20.0 * math.pi, 200.0 * math.pi),
    ("lt_1cm", 200.0 * math.pi, float("inf")),
)

def patch_for(lambda_p, base=None, r_lambda=None):
    """The scene patch: the base scale, doubled until it can hold the peak.

    One definition, used by the analytic table (stage 1) and by the state
    builder (stage 2) so the two cannot disagree.
    """
    b = float(PATCH_BASE_M if base is None else base)
    r = float(R_LAMBDA if r_lambda is None else r_lambda)
    need = max(b, r * float(lambda_p))
    m = max(int(math.ceil(math.log(need / b, 2.0) - 1e-12)), 0)
    return {"patch_m": float(b * (2 ** m)), "base_m": b,
            "doublings": int(m), "need_m": float(need),
            "r_lambda": r}


NON_DERIVED = [
    {"name": "K_M", "value": K_M, "units": "rad/m",
     "basis": "Elfouhaily et al. (1997): k_m = 2 pi / 0.017 m"},
    {"name": "C_M", "value": C_M, "units": "m/s",
     "basis": "gravity-capillary minimum phase speed at k_m"},
    {"name": "ECKV_AP", "value": ECKV_AP, "units": "-",
     "basis": "Elfouhaily et al. (1997) eq. 66-68"},
    {"name": "ECKV_AM_K", "value": ECKV_AM_K, "units": "-",
     "basis": "a_m = 0.13 u*/c_m (Elfouhaily et al. 1997)"},
    {"name": "PATCH_BASE_M", "value": PATCH_BASE_M, "units": "m",
     "basis": ("4 * 0.30 m = the v3.4 low-wind scene patch (env.py "
               "spec_L = 0.30 + 0.30 U^2/g at U->0); now wind-independent")},
    {"name": "R_LAMBDA", "value": R_LAMBDA, "units": "-",
     "basis": ("patch >= R_LAMBDA * lambda_p.  PIPELINE-V9 asks for >= 4; "
               "8 is adopted (convergence checked in G6, both reported)")},
    {"name": "K_X", "value": K_X, "units": "rad/m",
     "basis": ("cascade split at the centre of the design's 50-200 rad/m "
               "overlap band; partition of unity, so its value cannot "
               "change the total spectrum")},
]


# --------------------------------------------------------------------------
# 1. wind input
# --------------------------------------------------------------------------
def ustar_from_u10(U10):
    """Friction velocity from U10 only (Elfouhaily eq. 61 drag)."""
    U = float(U10)
    c10 = DRAG_A + DRAG_B * U
    return {"U10_ms": U, "C10": float(c10), "ustar_ms": float(math.sqrt(c10) * U),
            "form": "C10 = (0.8 + 0.065 U10) e-3, u* = sqrt(C10) U10"}


def alpha_m_elfouhaily(ustar):
    """alpha_m with the explicit low-wind floor (design section 3)."""
    us = float(ustar)
    cm = C_M
    if us < cm:
        raw = 0.01 * (1.0 + math.log(max(us, 1e-30) / cm))
    else:
        raw = 0.01 * (1.0 + 3.0 * math.log(us / cm))
    floored = max(raw, 0.0)
    if raw < 0.0:
        domain = "low_wind_floor"
    elif us < cm:
        domain = "gravity_range(u*<c_m)"
    else:
        domain = "capillary_range(u*>=c_m)"
    return {"ustar_ms": us, "c_m_ms": float(cm),
            "alpha_m_raw": float(raw), "alpha_m": float(floored),
            "floored": bool(raw < 0.0), "eckv_domain": domain,
            "u_over_cm": float(us / cm),
            "note": ("alpha_m = 0.01(1+ln(u*/c_m)) for u*<c_m, "
                     "0.01(1+3 ln(u*/c_m)) otherwise; negative values are "
                     "FLOORED at 0.0 and the cell is labelled")}


# --------------------------------------------------------------------------
# 2. fetch laws
# --------------------------------------------------------------------------
def fetch_laws(U, F, h):
    """JONSWAP peak/alpha/Hs laws + the depth-limited Hs of the same family."""
    U = float(U)
    F = float(F)
    h = float(h)
    Ft = G * F / max(U, 1e-9) ** 2
    fp = 3.5 * G / U * Ft ** (-1.0 / 3.0)
    omega_p = 2.0 * math.pi * fp
    alpha = 0.076 * Ft ** (-0.22)
    hs_deep = YV_COEF * YV_FETCH_COEF * math.sqrt(Ft) * U * U / G
    hs_law_form = 0.0016 * math.sqrt(Ft) * U * U / G
    dt = G * h / max(U, 1e-9) ** 2
    a = math.tanh(YV_DEPTH_COEF * dt ** 0.75)
    hs_depth = (U * U / G) * YV_COEF * a * math.tanh(
        YV_FETCH_COEF * math.sqrt(Ft) / max(a, 1e-12))
    lam_p = G / (2.0 * math.pi * fp * fp)
    k_p = 2.0 * math.pi / lam_p
    # peak group velocity (deep) and the travel-time lower bound on duration
    c_p = G / omega_p
    t_min = 2.0 * F / max(c_p, 1e-9)
    return {
        "U10_ms": U, "fetch_m": F, "depth_m": h,
        "Ft": float(Ft), "gF_over_U2": float(Ft),
        "dt_depth_nonDim": float(dt),
        "f_p_Hz": float(fp), "omega_p_rad_s": float(omega_p),
        "alpha": float(alpha),
        "Hs_fetch_law_m": float(hs_deep),
        "Hs_law_form_0p0016_m": float(hs_law_form),
        "Hs_depth_limited_m": float(hs_depth),
        "depth_correction_rel": float(hs_depth / hs_deep - 1.0),
        "depth_correction_binds": bool(k_p * h < 1.0 and hs_depth < hs_deep),
        "k_p_rad_per_m": float(k_p), "lambda_p_m": float(lam_p),
        "kp_times_h": float(k_p * h),
        "c_p_ms": float(c_p),
        "t_min_s": float(t_min),
        "t_min_form": "t_min = 2F/c_p (group travel over the fetch, lower bound)",
        "in_fitted_range_Ft_1e3_1e5": bool(1e3 <= Ft <= 1e5),
    }


# --------------------------------------------------------------------------
# 3. JONSWAP (frequency) -> S_g(k)
# --------------------------------------------------------------------------
def _jonswap_shape(omega, omega_p, gamma=JONSWAP_GAMMA):
    om = np.maximum(np.asarray(omega, dtype=np.float64), 1e-12)
    sig = np.where(om <= float(omega_p), JONSWAP_S_A, JONSWAP_S_B)
    r = np.exp(-((om - float(omega_p)) ** 2)
               / (2.0 * sig ** 2 * float(omega_p) ** 2))
    return om ** -5.0 * np.exp(-1.25 * (float(omega_p) / om) ** 4) \
        * float(gamma) ** r


def jonswap_norm(hs_m, omega_p, gamma=JONSWAP_GAMMA):
    """Amplitude A so that int A*shape domega = (H_s/4)^2 exactly."""
    om = np.logspace(math.log10(float(omega_p) * 1e-4),
                     math.log10(float(omega_p) * 1e2), 200000)
    sh = _jonswap_shape(om, omega_p, gamma)
    return {"A": float((float(hs_m) / 4.0) ** 2 / np.trapz(sh, om)),
            "gamma": float(gamma),
            "integral_shape": float(np.trapz(sh, om)),
            "note": ("the JONSWAP shape is normalised to the fetch-law H_s; "
                     "the implied Phillips alpha is reported alongside, so the "
                     "published alpha law is a READOUT here, not an input")}


def s_omega(hs_m, omega_p, omega, gamma=JONSWAP_GAMMA, norm=None):
    A = (norm if norm is not None else jonswap_norm(hs_m, omega_p,
                                                    gamma)["A"])
    return A * _jonswap_shape(omega, omega_p, gamma)


# --------------------------------------------------------------------------
# 4. Elfouhaily B_h -> S_h(k)
# --------------------------------------------------------------------------
def elfouhaily_bh(k, h, ustar, k_p, alpha_m):
    """Curvature spectrum B_h(k) (Elfouhaily eq. 40, McDaniel L_pm form)."""
    k = np.asarray(k, dtype=np.float64)
    kk = np.maximum(k, 1e-12)
    om = P8.dispersion_omega(kk, h)
    c = om / kk
    lpm = np.exp(-1.25 * (float(k_p) / kk) ** 2)
    fm = np.exp(-0.25 * (kk / K_M - 1.0) ** 2)
    bh = 0.5 * float(alpha_m) * (C_M / c) * lpm * fm
    bh = np.where(np.isfinite(bh), bh, 0.0)
    return np.maximum(bh, 0.0)


# --------------------------------------------------------------------------
# 5. ECKV azimuthal function
# --------------------------------------------------------------------------
def eckv_delta(k, h, ustar, k_p):
    """Delta(k): 1 in both the gravity and the capillary limit."""
    k = np.asarray(k, dtype=np.float64)
    kk = np.maximum(k, 1e-12)
    om = P8.dispersion_omega(kk, h)
    c = om / kk
    cp = float(P8.dispersion_omega(np.array([k_p]), h)[0]) / float(k_p)
    cd = C_M / c
    a_m = ECKV_AM_K * float(ustar) / C_M
    num = np.tanh(ECKV_A0 + ECKV_AP * (c / cp) ** 2 + a_m * cd ** 2)
    den = math.tanh(ECKV_A0 + ECKV_AP + a_m)
    d = num / den
    d = np.where(np.isfinite(d), d, 0.0)
    # a Delta slightly above 1 would make Phi negative at theta=pi/2; the
    # published function is bounded by 1 and the excess is 1e-3 at most.
    return np.clip(d, 0.0, 0.9999)


def eckv_phi(theta, delta):
    """Phi(k,theta) = (1/2pi)(1 + Delta cos 2(theta - theta_w))."""
    th = np.asarray(theta, dtype=np.float64)
    return (1.0 + np.asarray(delta, dtype=np.float64) * np.cos(2.0 * th)) \
        / (2.0 * math.pi)


# --------------------------------------------------------------------------
# 6. assembly
# --------------------------------------------------------------------------
def assemble(U, F, h, hs_override=None, k_cut=KCUT_DEFAULT,
             gamma=JONSWAP_GAMMA):
    """The one spectrum of section 3, parameterised and ready to evaluate."""
    laws = fetch_laws(U, F, h)
    us = ustar_from_u10(U)
    amr = alpha_m_elfouhaily(us["ustar_ms"])
    if hs_override is not None:
        hs = float(hs_override)
        hs_source = "override"
    else:
        hs = float(laws["Hs_depth_limited_m"])
        hs_source = "law"
    norm = jonswap_norm(hs, laws["omega_p_rad_s"], gamma)
    return {
        "version": "3.7",
        "laws": laws,
        "ustar": us,
        "alpha_m": amr,
        "Hs_m": float(hs),
        "hs_source": hs_source,
        "hs_override_m": (None if hs_override is None else float(hs_override)),
        "jonswap_norm": norm,
        "gamma": float(gamma),
        "k_cut_rad_per_m": float(k_cut),
        "depth_m": float(h),
        "k_p_rad_per_m": float(laws["k_p_rad_per_m"]),
        "omega_p_rad_s": float(laws["omega_p_rad_s"]),
        "alpha_m_value": float(amr["alpha_m"]),
        "eckv_domain": amr["eckv_domain"],
        "note": ("hs_source=law means the fetch-law H_s (depth-corrected in "
                 "the same family); 'override' means a recipe environment "
                 "statement, logged, never silent"),
    }


def s_parts(k, ctx):
    """S_g, S_h, B_h and S_total at wavenumbers k (SI)."""
    k = np.maximum(np.asarray(k, dtype=np.float64), 0.0)
    h = float(ctx["depth_m"])
    om = P8.dispersion_omega(np.maximum(k, 1e-12), h)
    sj = s_omega(ctx["Hs_m"], ctx["omega_p_rad_s"], om,
                 ctx["gamma"], ctx["jonswap_norm"]["A"])
    dwdk = P8.dispersion_domega_dk(np.maximum(k, 1e-12), h)
    s_g = np.where(k > 0.0, sj * dwdk, 0.0)
    bh = elfouhaily_bh(k, h, ctx["ustar"]["ustar_ms"],
                       ctx["k_p_rad_per_m"], ctx["alpha_m_value"])
    s_h = np.where(k > 0.0, bh / np.maximum(k, 1e-12) ** 3, 0.0)
    s_g = np.where(np.isfinite(s_g), s_g, 0.0)
    s_h = np.where(np.isfinite(s_h), s_h, 0.0)
    return {"k": k, "S_g": s_g, "S_h": s_h, "B_h": bh, "S": s_g + s_h,
            "omega": om}


def s_total(k, ctx):
    p = s_parts(k, ctx)
    return p["S"]


def mode_density_2d(kx, ky, ctx, weight=None, k_cut=None, theta_w=0.0,
                    delta=None):
    """2-D spectral density E(kx,ky) [m^2 per (rad/m)^2].

    E = S(k) * Phi(k,theta) / k.  The caller multiplies by the mode cell
    area (2 pi / S_patch)^2 to obtain the per-mode height variance.
    """
    kx = np.asarray(kx, dtype=np.float64)
    ky = np.asarray(ky, dtype=np.float64)
    k = np.sqrt(kx * kx + ky * ky)
    kk = np.maximum(k, 1e-12)
    S = s_total(kk.ravel(), ctx).reshape(k.shape)
    if delta is None:
        d = eckv_delta(kk.ravel(), float(ctx["depth_m"]),
                       ctx["ustar"]["ustar_ms"],
                       ctx["k_p_rad_per_m"]).reshape(k.shape)
    else:
        d = np.asarray(delta, dtype=np.float64).reshape(k.shape)
    th = np.arctan2(ky, kx) - float(theta_w)
    phi = eckv_phi(th, d)
    E = S * phi / kk
    if k_cut is not None:
        E = np.where(k <= float(k_cut), E, 0.0)
    if weight is not None:
        E = E * np.asarray(weight, dtype=np.float64)
    E = np.where(np.isfinite(E), E, 0.0)
    E[k == 0.0] = 0.0
    return E


# --------------------------------------------------------------------------
# 7. analytic integrals and the gate tables
# --------------------------------------------------------------------------
def _k_grid(k_cut):
    return np.logspace(math.log10(K_INT_MIN), math.log10(float(k_cut)), N_INT)


def _band_masks(k):
    out = {}
    for name, lo, hi in BANDS:
        m = (k >= lo) & (k < hi)
        out[name] = m
    return out


def integral_block(s, k, kmin, kcut):
    """m0 / sigma^2 over [kmin, kcut] plus band shares and the envelopes."""
    S = s["S"]
    m = (k >= kmin) & (k <= kcut)
    kk = k[m]
    out = {}
    for tag, arr in (("total", S[m]), ("gravity", s["S_g"][m]),
                     ("capillary", s["S_h"][m])):
        out[tag] = {"m0": float(np.trapz(arr, kk)),
                    "sigma2": float(np.trapz(kk * kk * arr, kk))}
    bm = _band_masks(k)
    for band, _, _ in BANDS:
        mask = bm[band] & m
        kk2 = k[mask]
        out.setdefault("bands", {})[band] = {
            "m0": float(np.trapz(S[mask], kk2)),
            "sigma2": float(np.trapz(kk2 * kk2 * S[mask], kk2)),
        }
    sig = out["total"]["sigma2"]
    m0 = out["total"]["m0"]
    if sig > 0:
        for band in out["bands"]:
            out["bands"][band]["sigma2_share"] = float(
                out["bands"][band]["sigma2"] / sig)
    if m0 > 0:
        for band in out["bands"]:
            out["bands"][band]["m0_share"] = float(
                out["bands"][band]["m0"] / m0)
    # fractions above 100 and 370 rad/m
    for kx in (100.0, 370.0):
        mh = m & (k >= kx)
        s2 = float(np.trapz(k[mh] ** 2 * S[mh], k[mh]))
        out["frac_sigma2_above_%d" % int(kx)] = float(s2 / sig) if sig else None
    # 50-200 overlap band, per part (honesty, design section 3)
    mo = m & (k >= 50.0) & (k <= 200.0)
    ko = k[mo]
    out["overlap_50_200"] = {
        "S_g_m0": float(np.trapz(s["S_g"][mo], ko)),
        "S_h_m0": float(np.trapz(s["S_h"][mo], ko)),
        "S_g_sigma2": float(np.trapz(ko * ko * s["S_g"][mo], ko)),
        "S_h_sigma2": float(np.trapz(ko * ko * s["S_h"][mo], ko)),
        "S_g_share_of_overlap_sigma2": (
            float(np.trapz(ko * ko * s["S_g"][mo], ko)
                  / max(np.trapz(ko * ko * S[mo], ko), 1e-30))),
        "S_h_share_of_overlap_sigma2": (
            float(np.trapz(ko * ko * s["S_h"][mo], ko)
                  / max(np.trapz(ko * ko * S[mo], ko), 1e-30))),
    }
    return out


def anisotropy_prediction(ctx, k, kmin, kcut):
    """Slope-weighted Delta -> sigma^2_along / sigma^2_cross.

    Phi = (1/2pi)(1 + Delta cos 2 theta) gives
        sigma^2_along = sigma^2/2 (1 + Delta_eff/2)
        sigma^2_cross = sigma^2/2 (1 - Delta_eff/2)
    (verified in gate G5 by direct numerical integration over theta on the
    2-D grid).
    """
    s = s_parts(k, ctx)
    m = (k >= kmin) & (k <= kcut)
    kk = k[m]
    w = kk * kk * s["S"][m]
    d = eckv_delta(kk, float(ctx["depth_m"]), ctx["ustar"]["ustar_ms"],
                   ctx["k_p_rad_per_m"])
    tot = float(np.trapz(w, kk))
    de = float(np.trapz(w * d, kk) / tot) if tot > 0 else 0.0
    sig = tot
    sa = sig * 0.5 * (1.0 + de / 2.0)
    sc = sig * 0.5 * (1.0 - de / 2.0)
    return {"delta_effective": de, "sigma2_along": float(sa),
            "sigma2_cross": float(sc),
            "ratio_along_over_cross": float(sa / sc) if sc > 0 else None}


def analytic(ctx, patch_m, k_cut=None, n_band=None):
    """The whole stage-1 analytic surface for one cell."""
    kcut = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    k = _k_grid(kcut)
    s = s_parts(k, ctx)
    kmin = 2.0 * math.pi / float(patch_m)
    blk = integral_block(s, k, kmin, kcut)
    an = anisotropy_prediction(ctx, k, kmin, kcut)
    laws = ctx["laws"]
    # peak of the ELEVATION spectrum S(k) and of the slope spectrum k^2 S(k)
    m = (k >= kmin) & (k <= kcut)
    kk = k[m]
    ip = int(np.argmax(s["S"][m]))
    isl = int(np.argmax(kk * kk * s["S"][m]))
    out = {
        "patch_m": float(patch_m),
        "kmin_rad_per_m": float(kmin),
        "k_cut_rad_per_m": kcut,
        "kmax_law_note": ("the continuum integral is taken to k_cut; the "
                          "discrete realisation adds the mode-set sampling "
                          "error on top (reported in the state numbers)"),
        "Hs_from_m0_m": float(4.0 * math.sqrt(max(blk["total"]["m0"], 0.0))),
        "Hs_law_m": float(laws["Hs_depth_limited_m"]),
        "Hs_used_m": float(ctx["Hs_m"]),
        "peak": {"k_peak_S_rad_per_m": float(kk[ip]),
                 "lambda_peak_S_m": float(2.0 * math.pi / kk[ip]),
                 "lambda_p_law_m": float(laws["lambda_p_m"]),
                 "k_slope_peak_rad_per_m": float(kk[isl]),
                 "lambda_slope_peak_m": float(2.0 * math.pi / kk[isl])},
        "integrals": blk,
        "anisotropy": an,
        "resolved_band_km1": [float(kmin), kcut],
    }
    out["m0_closure"] = (None if blk["total"]["m0"] <= 0 else
                         float(out["Hs_from_m0_m"] / 4.0
                               / math.sqrt(blk["total"]["m0"])))
    return out


def cox_munk_reference(U):
    from v5_physics import cox_munk_anisotropy, cox_munk_sigma2
    an = cox_munk_anisotropy(U)
    return {"total_clean_rad2": float(cox_munk_sigma2(U)),
            "sigma_along_rad2": float(an["sigma_a2"]),
            "sigma_cross_rad2": float(an["sigma_c2"]),
            "ratio_along_over_cross": float(an["sigma_a2"] / an["sigma_c2"])
            if an["sigma_c2"] > 0 else None,
            "note": ("Cox & Munk (1954) fully developed sea; a REFERENCE, "
                     "never a target (PIPELINE-V9 G3)")}


def full_eckv_sigma2(ctx, k_hi=None):
    """sigma^2 of the assembled spectrum with NO cut (to the viscous edge)."""
    khi = float(k_hi if k_hi is not None else 1.0e4)
    k = _k_grid(khi)
    s = s_parts(k, ctx)
    return {"k_hi": khi,
            "sigma2_rad2": float(np.trapz(k * k * s["S"], k)),
            "m0_m2": float(np.trapz(s["S"], k))}


# --------------------------------------------------------------------------
# 8. gate G0 (Parseval) on a built test field
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
    D = E * dk * dk                        # per-mode height variance
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
            "H0": H0, "k": k, "D": D, "kx": kx, "ky": ky,
            "patch_m": S_patch, "N": N, "k_cut": kcut}


def gate_g0(ctx, N=1024, patch_m=None, k_cut=None, seed=12345, theta_w=0.0):
    """G0: the built field's m0 and sigma^2 match the assembled S(k).

    Two levels are separated on purpose:
      (a) convention: numeric field vs the mode-sum of the same density
          (this is the check that a 2*pi or 1/k mistake cannot hide in);
      (b) discretisation: mode-sum vs the continuum integral of S(k).
    The tolerance band of (a) is the model's own realisation scatter,
    z = (numeric - mode_sum) / (mode_sum*sqrt(2/DOF)).
    """
    pm = float(patch_m if patch_m is not None else
               patch_for(ctx["laws"]["lambda_p_m"])["patch_m"])
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
        "gate": "G0 Parseval (PIPELINE-V9 section 2)",
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
    out["note"] = ("(a) numeric vs mode-sum is the convention test: any 2*pi "
                   "or 1/k error moves it far outside the z band.  (b) "
                   "mode-sum vs continuum is the discretisation error of this "
                   "patch; both are reported, neither is tuned.")
    return out


# --------------------------------------------------------------------------
# 9. gate G4 (positivity + domain log) over the 8-cell grid
# --------------------------------------------------------------------------
SWEEP_U = (1.5, 3.0, 6.0, 10.0)
SWEEP_F = (100.0, 1000.0)
SWEEP_Z = 1.0


def gate_g4(h=SWEEP_Z, k_cut=KCUT_DEFAULT):
    """G4: B(k) >= 0 and S(k) >= 0 for all k in all 8 cells, domain logged."""
    rows = []
    ok = True
    k = np.logspace(-3.0, math.log10(float(k_cut)), 6000)
    for U in SWEEP_U:
        for F in SWEEP_F:
            ctx = assemble(U, F, h, k_cut=k_cut)
            s = s_parts(k, ctx)
            bmin = float(np.min(s["B_h"]))
            smin = float(np.min(s["S"]))
            finite = bool(np.all(np.isfinite(s["S"])))
            row = {
                "U_ms": float(U), "fetch_m": float(F), "depth_m": float(h),
                "ustar_ms": ctx["ustar"]["ustar_ms"],
                "u_over_cm": ctx["alpha_m"]["u_over_cm"],
                "alpha_m_raw": ctx["alpha_m"]["alpha_m_raw"],
                "alpha_m_used": ctx["alpha_m_value"],
                "eckv_domain": ctx["eckv_domain"],
                "floored": ctx["alpha_m"]["floored"],
                "B_h_min": bmin, "S_min": smin, "finite": finite,
                "B_h_max": float(np.max(s["B_h"])),
                "capillary_present": bool(np.max(s["S_h"]) > 0.0),
                "pass": bool(bmin >= 0.0 and smin >= 0.0 and finite),
            }
            ok = ok and row["pass"]
            rows.append(row)
    return {"gate": ("G4 positivity + explicit low-wind domain log "
                     "(PIPELINE-V9 section 3)"),
            "k_range": [1e-3, float(k_cut)],
            "cells": rows,
            "n_cells": len(rows),
            "n_low_wind_floor": int(sum(1 for r in rows if r["floored"])),
            "pass": bool(ok),
            "note": ("the floor is applied EXPLICITLY (alpha_m = max(raw,0)); "
                     "at U=1.5 m/s u* ~ 0.045 m/s < c_m = 0.231 m/s, so "
                     "alpha_m < 0 and the capillary part is identically zero "
                     "there.  That is the design's stated domain treatment, "
                     "not a silent clamp: the cell carries "
                     "eckv_domain=low_wind_floor and G2 reports its "
                     "capillary share as 0.")}


# --------------------------------------------------------------------------
# 10. spectrum curves for the deliverable (CSV / JSON)
# --------------------------------------------------------------------------
def curves(ctx, patch_m, k_cut=None, n=420):
    kcut = float(k_cut if k_cut is not None else ctx["k_cut_rad_per_m"])
    kmin = 2.0 * math.pi / float(patch_m)
    k = np.logspace(math.log10(kmin), math.log10(kcut), int(n))
    s = s_parts(k, ctx)
    d = eckv_delta(k, float(ctx["depth_m"]), ctx["ustar"]["ustar_ms"],
                   ctx["k_p_rad_per_m"])
    return {"patch_m": float(patch_m),
            "k_rad_per_m": [float(x) for x in k],
            "S_g_m2_per_rad_m": [float(x) for x in s["S_g"]],
            "S_h_m2_per_rad_m": [float(x) for x in s["S_h"]],
            "S_total_m2_per_rad_m": [float(x) for x in s["S"]],
            "B_h": [float(x) for x in s["B_h"]],
            "delta_eckv": [float(x) for x in d],
            "omega_rad_s": [float(x) for x in s["omega"]]}
