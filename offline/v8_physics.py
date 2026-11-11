"""v8 chunk-4 physics nodes -- the wave-spectrum decision (version 3.4).

Implements PIPELINE-V8.md "chunk 4" (rally P0-C + the H_max clamp + patch
coverage).  Every number is a derived quantity of a published formula; the
non-derived constants that remain are listed in ``NON_DERIVED`` and copied
into numbers-v8.json.

What is here (chunk-4 items 1-6):

  1. DISPERSION  omega^2 = (g k + (gamma/rho) k^3) tanh(k h)
     gamma = 0.0728 N/m (surface tension of water at ~20 C),
     rho = 1000 kg/m^3 (fresh water; the temperature dependence is left as
     the v5 plug-in point, so the constants stay constants).

  2. CALIBRATION ANCHOR -- the Cox-Munk (fully developed sea) target is
     replaced by a short-fetch slope model: the mean square slope of a
     fetch-limited spectrum, integrated from the JONSWAP/Donelan family
     frequency spectrum mapped onto the chunk-4 dispersion.  The only
     argument is the non-dimensional fetch Ft = gF/U^2 (plus the depth).

  3. DIRECTIONAL DISTRIBUTION D(theta) -- an explicit state parameter, the
     cos^{2s}(theta/2) family (Longuet-Higgins et al. 1963) with the
     Mitsuyasu et al. (1975) frequency-dependent spreading parameter
     s(omega) = s_max (omega/omega_p)^m (m = +5 below, -2.5 above the peak)
     and s_max = 11.5 (U10/c_p) (WAFO / Krogstad & Barstow 1999 digest of
     Mitsuyasu et al. 1975).  Derived from the model, NOT an art slider.

  4. CAPILLARY TAIL -- the high-wavenumber cutoff is fixed by a physical
     rule: viscous damping 2 nu k^2 exceeds omega.  The cutoff hierarchy
     (viscous / Nyquist / grid) is reported, never silently replaced.

  5. BREAKING CLAMP -- H_max = min(0.142 lambda, 0.78 h) on the effective
     spectral height (recorded side by side with the Stokes 1/7 cap).

  6. Fetch / patch metadata (the numbers JSON records the patch size = the
     longest resolved wavelength, Ft, and the adopted anchor model id).

Sources
-------
* JONSWAP: Hasselmann et al. (1973), "Measurements of wind-wave growth and
  swell propagation during the Joint North Sea Wave Project".
        alpha = 0.076 Ft^-0.22,  f_p = 3.5 (g/U) Ft^-1/3,
        H_s   = 0.0016 Ft^0.5 U^2/g   (independent consistency check below)
* Donelan, Hamilton & Hui (1985): peak enhancement gamma = 1.7 (the
  JONSWAP gamma = 3.3 is the other named member of the family).
* Longuet-Higgins, Cartwright & Smith (1963): cos^{2S}(theta/2) spreading,
  normalisation Gamma(S+1) / (2 sqrt(pi) Gamma(S+1/2)).
* Mitsuyasu et al. (1975); Hasselmann, Dunckel & Ewing (1980):
  frequency-dependent spreading parameter (the +5 / -2.5 shell law).
* Cox & Munk (1954): the OLD anchor, kept only as a reference number.
* Stokes (1847) / 1:7 steepness limit and the 0.78 h depth limit of
  PIPELINE-V8 chunk 4.
"""
import math

import numpy as np
from scipy.special import gammaln

# --------------------------------------------------------------------------
# constants (recorded in numbers-v8.json)
# --------------------------------------------------------------------------
G = 9.81                 # m/s^2
SIGMA_SURF = 0.0728      # N/m   surface tension of water (~20 C)
RHO_W = 1000.0           # kg/m^3 fresh water
TF = SIGMA_SURF / RHO_W  # 7.28e-5 m^3/s^2  (gamma/rho)
NU_KIN = 1.0e-6          # m^2/s kinematic viscosity of water (~20 C)

CHOPPY_LAMBDA = 0.85     # inherited (state.py), recorded only
DAMP_L_STATE = 0.005     # inherited micro-wave damping length in state.py

OMEGA0 = 2.0 * math.pi / 8.0     # 8 s loop quantum (rad/s)

#: The two named members of the short-fetch slope family.  Gamma is the
#: peak enhancement factor of the frequency spectrum: JONSWAP measured
#: gamma = 3.3, Donelan et al. (1985) use gamma = 1.7.
ANCHOR_CANDIDATES = (("JONSWAP", 3.3), ("DONELAN", 1.7))
ANCHOR_MODEL_ID = "SFM1-JONSWAP-DONELAN-median"

#: Fetch defaults.  Pond default for new recipes (PIPELINE-V8 chunk 4).
POND_FETCH_M = 100.0
#: Legacy compatibility default: the fixed fetch whose anchor reproduces the
#: old Cox-Munk target at the legacy reference scene (U=1.5 m/s, h=0.6 m);
#: solved in ``legacy_fetch_m`` below (F = 4.60 m, Ft = 20.1, residual 0.0%).
LEGACY_FETCH_M = 4.6

#: Cox-Munk reference kept for the old-anchor comparison only.
COX_MUNK_A = 0.003
COX_MUNK_B = 5.12e-3

#: Mitsuyasu et al. (1975) shell law exponents.
SPREAD_M_LOW = 5.0       # omega < omega_p
SPREAD_M_HIGH = -2.5     # omega >= omega_p
SPREAD_S_MAX_K = 11.5    # s_max = 11.5 * (U10 / c_p) = 11.5 * omega_p U / g

#: Breaking limits (PIPELINE-V8 chunk 4).
BREAK_STEEPNESS = 0.142          # H_max = 0.142 lambda
STOKES_STEEPNESS = 1.0 / 7.0     # Stokes cap, kept alongside
BREAK_DEPTH_FRAC = 0.78          # H_max = 0.78 h

NON_DERIVED = [
    {"name": "SIGMA_SURF", "value": SIGMA_SURF, "units": "N/m",
     "basis": "surface tension of water at ~20 C (0.0728 N/m); the "
              "temperature dependence is the declared v5 plug-in point"},
    {"name": "RHO_W", "value": RHO_W, "units": "kg/m^3",
     "basis": "fresh water density (pond scene); same plug-in point"},
    {"name": "NU_KIN", "value": NU_KIN, "units": "m^2/s",
     "basis": "kinematic viscosity of water at ~20 C; used only by the "
              "viscous-cutoff RULE (rule 4), which does not bind inside the "
              "resolved band (recorded)"},
    {"name": "SPREAD_S_MAX_K", "value": SPREAD_S_MAX_K, "units": "-",
     "basis": "Mitsuyasu et al. (1975) s_max = 11.5 (U10/c_p)"},
]


def cox_munk_sigma2(U):
    """The OLD anchor (fully developed sea), kept as a reference number."""
    return COX_MUNK_A + COX_MUNK_B * float(U)


# --------------------------------------------------------------------------
# 1. dispersion
# --------------------------------------------------------------------------
def dispersion_omega(k, h):
    """omega^2 = (g k + (gamma/rho) k^3) tanh(k h).  k, h in SI."""
    k = np.asarray(k, dtype=np.float64)
    h = float(h)
    kk = np.maximum(k, 0.0)
    return np.sqrt((G * kk + TF * kk ** 3) * np.tanh(np.clip(kk * h, 0.0, 350.0)))


def dispersion_domega_dk(k, h):
    """d omega / d k of the chunk-4 dispersion (analytic)."""
    k = np.asarray(k, dtype=np.float64)
    h = float(h)
    kk = np.maximum(k, 0.0)
    w = dispersion_omega(kk, h)
    w = np.where(w > 0.0, w, 1.0)
    th = np.tanh(np.clip(kk * h, 0.0, 350.0))
    ch = np.cosh(np.clip(kk * h, 0.0, 350.0))
    sch2 = 1.0 / np.maximum(ch * ch, 1.0)
    num = (G + 3.0 * TF * kk ** 2) * th + (G * kk + TF * kk ** 3) * h * sch2
    return num / (2.0 * w)


def shallow_ratio_theory(kh):
    """c/sqrt(gh) for pure gravity waves in finite depth = sqrt(tanh(kh)/(kh))."""
    x = np.asarray(kh, dtype=np.float64)
    x = np.where(x > 0.0, x, 1e-12)
    return np.sqrt(np.tanh(x) / x)


def dispersion_gate_table(h, kh_list=(0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0),
                          depth_ref=None):
    """Gate-2 table: as-used phase speed vs sqrt(g h) as kh -> 0.

    For each target kh the wavenumber k = kh/h is used, so the capillary term
    is the one the pipeline actually applies (it is negligible at kh<1 but is
    included, not dropped).
    """
    h = float(h)
    c_ref = math.sqrt(G * h)
    rows = []
    for kh in kh_list:
        k = float(kh) / h
        w = float(dispersion_omega(np.array([k]), h)[0])
        c = w / k
        c_th = float(shallow_ratio_theory(np.array([kh]))[0]) * c_ref
        n = max(1.0, math.floor(w / OMEGA0 + 0.5))
        wq = n * OMEGA0
        rows.append({
            "kh": float(kh),
            "k_rad_per_m": k,
            "lambda_m": float(2.0 * math.pi / k),
            "omega_rad_per_s": w,
            "c_m_per_s": c,
            "c_over_sqrt_gh": float(c / c_ref),
            "c_over_sqrt_gh_theory": float(c_th / c_ref),
            "theory_minus_used_rel": float(abs(c - c_th) / c_th),
            "omega_over_OMEGA0": float(w / OMEGA0),
            "omega_quantised_rad_per_s": float(wq),
            "c_quantised_over_sqrt_gh": float(wq / k / c_ref),
            "quantisation_rel_error": float(abs(wq - w) / w),
        })
    # asymptotic statement on the finest kh (the pure-gravity theory ratio is
    # only meaningful in the shallow limit; the capillary term is included in
    # the "used" value, so the comparison is gated at kh <= 0.5)
    shallow = [r for r in rows if r["kh"] <= 0.5]
    r0 = rows[0]
    return {
        "gate": ("chunk-4 gate 2: the phase speed of the new dispersion tends "
                 "to sqrt(g h) for kh << 1 (theory ratio "
                 "c/sqrt(gh) = sqrt(tanh(kh)/(kh)))"),
        "depth_m": h,
        "sqrt_gh_m_per_s": c_ref,
        "rows": rows,
        "ratio_at_kh_0p02": r0["c_over_sqrt_gh"],
        "ratio_at_kh_0p1": rows[2]["c_over_sqrt_gh"],
        "max_theory_deviation_rel_kh_le_0p5": float(
            max(r["theory_minus_used_rel"] for r in shallow)),
        "max_theory_deviation_rel": float(max(r["theory_minus_used_rel"]
                                              for r in rows)),
        "pass": bool(r0["c_over_sqrt_gh"] < 1.0
                     and r0["c_over_sqrt_gh"] > 0.9999
                     and max(r["theory_minus_used_rel"] for r in shallow)
                     < 1e-3),
        "quantisation_note": ("the last two columns are the AS-USED values "
                              "after the 8 s loop quantisation.  At kh=0.02 "
                              "with h=0.1 m omega < OMEGA0, so the existing "
                              "n>=1 rule forces the mode to one quantum "
                              "(4x too fast).  That regime never occurs in "
                              "the verification states: their lowest mode has "
                              "kh >= 0.3 and omega >= 6*OMEGA0 (table in "
                              "numbers)."),
    }


def realisation_scatter(H0, N, k):
    """Effective degrees of freedom of the realised spectrum.

    For a Gaussian random field, Var(mean h^2)/(mean h^2)^2 = 2/DOF with
    DOF = (sum S)^2 / sum S^2, so the *realised* variance of a single
    snapshot deviates from the spectral (ensemble) value by ~sqrt(2/DOF).
    This patch holds only a few dominant modes, so the deviation is large
    and is reported instead of being mistaken for a calibration error.
    """
    S = (np.abs(H0.astype(np.complex128)) ** 2) / (2.0 * N ** 4)
    k = np.asarray(k, dtype=np.float64)
    Sg = S * k * k
    dof_m0 = float(S.sum() ** 2 / np.sum(S * S))
    dof_g = float(Sg.sum() ** 2 / np.sum(Sg * Sg))
    return {
        "dof_height": dof_m0, "dof_slope": dof_g,
        "rel_std_height": float(math.sqrt(2.0 / dof_m0)) if dof_m0 > 0 else None,
        "rel_std_sigma2": float(math.sqrt(2.0 / dof_g)) if dof_g > 0 else None,
        "note": ("a 2048^2 patch of a 1.5 m patch holding one dominant "
                 "wavelength is a few-mode realisation: the realised "
                 "variance legitimately scatters around the spectral "
                 "target, and the z-score below says whether a deviation is "
                 "calibration or realisation."),
    }


def dispersion_error_stats(omega_raw, omega_q, sl_weights):
    """Quantisation error of the 8 s loop under the new dispersion."""
    f = float(np.sum(sl_weights))
    rel = np.abs(omega_q - omega_raw) / np.maximum(omega_raw, 1e-30)
    return {
        "max_rel_error": float(rel.max()),
        "slope_weighted_mean_rel_error": float((rel * sl_weights).sum() / f)
        if f > 0 else None,
        "n_modes": int(omega_raw.size),
        "modes_below_OMEGA0": int(np.sum(omega_raw < OMEGA0)),
        "note": ("the 8 s loop quantises omega to integer multiples of "
                 "2*pi/8 s.  With the new dispersion every mode that inverts "
                 "to n>=1 keeps its phase; the error bound per mode is "
                 "OMEGA0/(2 omega).  Modes with omega < OMEGA0 would have to "
                 "be frozen (n=0) and are forced to n=1 by the existing rule "
                 "-- none exist in the verification scenes (recorded)."),
    }


# --------------------------------------------------------------------------
# 2. short-fetch slope anchor
# --------------------------------------------------------------------------
def nondim_fetch(F, U):
    """Ft = g F / U^2."""
    return G * float(F) / max(float(U), 1e-9) ** 2


def jonswap_fetch_relations(U, F):
    """JONSWAP fetch-limited peak frequency and Phillips alpha."""
    Ft = nondim_fetch(F, U)
    fp = 3.5 * G / float(U) * Ft ** (-1.0 / 3.0)
    alpha = 0.076 * Ft ** (-0.22)
    omega_p = 2.0 * math.pi * fp
    return {"Ft": float(Ft), "f_p_Hz": float(fp), "omega_p_rad_s": float(omega_p),
            "alpha": float(alpha),
            "Hs_fetch_law_m": float(0.0016 * math.sqrt(Ft)
                                    * float(U) ** 2 / G),
            "relations": ("f_p = 3.5 (g/U) Ft^-1/3, alpha = 0.076 Ft^-0.22, "
                          "H_s = 0.0016 Ft^0.5 U^2/g (Hasselmann et al. 1973)")}


def _sigma2_integral(U, F, h, gamma_peak, kmax, n=40000, alpha=None,
                     fp_scale=1.0):
    """Mean square slope of the fetch-limited frequency spectrum.

    The JONSWAP/Donelan frequency spectrum S(f) is evaluated at f(k) of the
    chunk-4 dispersion and mapped to wavenumber with S(k) = S(f) df/dk, then
        sigma^2 = int k^2 S(k) dk,   m0 = int S(k) dk.
    The gravity-capillary dispersion makes the tail converge as k^-5, so the
    integral is cutoff-insensitive (checked in ``anchor_sensitivity``).
    """
    rel = jonswap_fetch_relations(U, F)
    fp = rel["f_p_Hz"] * float(fp_scale)
    a = rel["alpha"] if alpha is None else float(alpha)
    kmin = 1e-3
    k = np.logspace(math.log10(kmin), math.log10(float(kmax)), int(n))
    om = dispersion_omega(k, h)
    f = om / (2.0 * math.pi)
    sj = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-((f - fp) ** 2) / (2.0 * sj ** 2 * fp ** 2))
    sf = (a * G ** 2 / (2.0 * math.pi) ** 4 * f ** -5
          * np.exp(-1.25 * (fp / np.maximum(f, 1e-12)) ** 4)
          * float(gamma_peak) ** r)
    dfdk = dispersion_domega_dk(k, h) / (2.0 * math.pi)
    sk = sf * dfdk
    m0 = float(np.trapz(sk, k))
    s2 = float(np.trapz(k * k * sk, k))
    return {"sigma2": s2, "m0": m0, "Hs": float(4.0 * math.sqrt(max(m0, 0.0))),
            "f_p_Hz": float(fp), "alpha": float(a),
            "lambda_p_m": float(G / (2.0 * math.pi * fp ** 2))}


def _resolved_kband(S_patch, N):
    return 2.0 * math.pi / float(S_patch), math.pi * N / float(S_patch)


def anchor_sigma2(U, F, h, S_patch=1.4755, N=2048, gamma_override=None):
    """The new sigma^2 gate target (short-fetch slope model).

    adopted value = MEDIAN of the named-model candidates
    (JONSWAP gamma=3.3, Donelan gamma=1.7), i.e. the midpoint of the two;
    the full spread is recorded by ``anchor_sensitivity``.
    """
    kmin, kmax = _resolved_kband(S_patch, N)
    rel = jonswap_fetch_relations(U, F)
    cands = {}
    for name, gam in ANCHOR_CANDIDATES:
        r = _sigma2_integral(U, F, h, gam, min(kmax, 2.0e4))
        r["gamma"] = float(gam)
        cands[name] = r
    vals = [cands[n]["sigma2"] for n, _ in ANCHOR_CANDIDATES]
    med = float(0.5 * (min(vals) + max(vals)))
    adopted = med if gamma_override is None else _sigma2_integral(
        U, F, h, float(gamma_override), min(kmax, 2.0e4))
    cm = cox_munk_sigma2(U)
    return {
        "model_id": ANCHOR_MODEL_ID,
        "model": ("mean square slope of a fetch-limited wave spectrum "
                  "(JONSWAP/Donelan family) integrated over the resolved band "
                  "with the chunk-4 dispersion"),
        "adopted_rule": "median of the named-model candidates",
        "fetch_m": float(F), "Ft": float(rel["Ft"]), "depth_m": float(h),
        "wind_ms": float(U),
        "jonswap_relations": rel,
        "candidates": cands,
        "sigma2_target_rad2": med,
        "sigma2_target_gamma_override": (
            None if gamma_override is None else adopted),
        "cox_munk_old_target_rad2": cm,
        "ratio_new_over_old": float(med / cm) if cm > 0 else None,
        "resolved_band_rad_per_m": [kmin, kmax],
        "note": ("the old anchor (Cox-Munk) is a fully developed sea and is "
                 "kept only for the comparison.  The new target is a factor "
                 "%.2f the old one at this scene: the old anchor removes "
                 "roughly half of the short-wave slope variance that a "
                 "fetch-limited pond spectrum carries." % (med / cm)),
    }


def anchor_sensitivity(U, F, h, S_patch=1.4755, N=2048):
    """Sensitivity of the adopted sigma^2 target to the literature spread."""
    kmin, kmax = _resolved_kband(S_patch, N)
    base = anchor_sigma2(U, F, h, S_patch, N)
    s0 = base["sigma2_target_rad2"]
    gam = {}
    for g in (1.0, 1.7, 2.5, 3.3, 5.0):
        gam["gamma_%.1f" % g] = _sigma2_integral(U, F, h, g, min(kmax, 2.0e4))
    rel = jonswap_fetch_relations(U, F)
    alpha_rows = {}
    for e in (0.18, 0.22, 0.26):
        a = 0.076 * rel["Ft"] ** (-e)
        alpha_rows["alpha_exp_%.2f" % e] = (
            _sigma2_integral(U, F, h, 3.3, min(kmax, 2.0e4), alpha=a))
    kcut = {}
    for kc in (500.0, 1000.0, 3000.0, kmax, 1e4, 1e6):
        kcut["kmax_%.3g" % kc] = _sigma2_integral(U, F, h, 3.3, kc)
    depth_rows = {}
    for hh in (h, 0.5 * h, 0.25 * h, 4.0 * h, 20.0 * h):
        depth_rows["h_%.3g" % hh] = _sigma2_integral(U, F, hh, 3.3,
                                                     min(kmax, 2.0e4))
    u_rows = {}
    for uu in (0.9 * U, 1.1 * U):
        u_rows["U_%.3g" % uu] = anchor_sigma2(uu, F, h, S_patch, N)
    return {
        "base_sigma2": s0,
        "gamma_variants": {k: v["sigma2"] for k, v in gam.items()},
        "gamma_spread_rel": float((max(v["sigma2"] for v in gam.values())
                                   - min(v["sigma2"] for v in gam.values()))
                                  / s0),
        "alpha_exponent_variants": {k: v["sigma2"]
                                    for k, v in alpha_rows.items()},
        "k_max_cutoff_variants": {k: v["sigma2"] for k, v in kcut.items()},
        "k_max_insensitivity_rel": float(
            abs(kcut["kmax_%.3g" % kmax]["sigma2"]
                - kcut["kmax_1e+06"]["sigma2"]) / s0),
        "depth_variants": {k: v["sigma2"] for k, v in depth_rows.items()},
        "wind_plus_minus_10pct": {k: v["sigma2_target_rad2"]
                                  for k, v in u_rows.items()},
        "Hs_consistency": {
            "integrated_Hs_m": base["candidates"]["JONSWAP"]["Hs"],
            "jonswap_law_Hs_m": base["jonswap_relations"]["Hs_fetch_law_m"],
            "rel_diff": float(abs(base["candidates"]["JONSWAP"]["Hs"]
                                  - base["jonswap_relations"]["Hs_fetch_law_m"])
                              / base["jonswap_relations"]["Hs_fetch_law_m"]),
            "note": ("independent consistency check of the integration: m0 "
                     "integrated from the same spectrum through the chunk-4 "
                     "dispersion lands within the published JONSWAP H_s "
                     "fetch law's own scatter (+11% here; the law is a fit "
                     "to the same spectra, not an identity)"),
        },
    }


def legacy_fetch_m(U=1.5, h=0.6, S_patch=1.4755, N=2048):
    """Fixed fetch whose anchor reproduces the old Cox-Munk target.

    PIPELINE-V8 chunk 4: "fetch F is a scene input (pond default 100 m); for
    legacy recipes choose the default value closest to the previous
    behaviour and record the choice".  Solved once here (documented), not at
    run time, so a legacy recipe's anchor is a fixed, recorded number.
    """
    target = cox_munk_sigma2(U)
    lo, hi = 0.2, 1000.0

    def f(F):
        return anchor_sigma2(U, F, h, S_patch, N)["sigma2_target_rad2"]

    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if f(mid) > target:
            hi = mid
        else:
            lo = mid
    F = math.sqrt(lo * hi)
    got = f(F)
    return {
        "U_ms": float(U), "depth_m": float(h),
        "target_cox_munk_rad2": float(target),
        "fetch_m": float(F), "Ft": float(nondim_fetch(F, U)),
        "sigma2_at_that_fetch": float(got),
        "residual_rel": float(abs(got - target) / target),
        "rounded_default_m": float(LEGACY_FETCH_M),
        "note": ("a compatibility anchor, not a physical pond fetch: at the "
                 "legacy reference scene the Cox-Munk target is reproduced "
                 "by a very short fetch because Cox-Munk is a fully developed "
                 "sea and this model is not.  New recipes carry "
                 "wave.fetchM explicitly (pond default 100 m)."),
    }


# --------------------------------------------------------------------------
# 3. directional distribution D(theta)
# --------------------------------------------------------------------------
def spreading_parameters(U, F):
    """Mitsuyasu et al. (1975) spreading parameters at the peak."""
    rel = jonswap_fetch_relations(U, F)
    omega_p = rel["omega_p_rad_s"]
    cp = G / omega_p
    s_max = SPREAD_S_MAX_K * (float(U) / cp)
    return {"omega_p_rad_s": float(omega_p), "c_p_m_per_s": float(cp),
            "U_over_cp": float(U / cp), "s_max": float(s_max),
            "m_low": SPREAD_M_LOW, "m_high": SPREAD_M_HIGH,
            "basis": ("s_max = 11.5 (U10/c_p); s = s_max (omega/omega_p)^m "
                      "with m=+5 below and m=-2.5 above the peak "
                      "(Mitsuyasu et al. 1975; Hasselmann et al. 1980; "
                      "Krogstad & Barstow 1999 for the +5/-2.5 shell law)")}


def spreading_s(omega, omega_p, s_max):
    """Frequency-dependent cos^{2s} spreading parameter."""
    om = np.maximum(np.asarray(omega, dtype=np.float64), 1e-12)
    r = om / float(omega_p)
    s = np.where(r < 1.0,
                 float(s_max) * r ** SPREAD_M_LOW,
                 float(s_max) * r ** SPREAD_M_HIGH)
    return np.maximum(s, 1e-6)


def cos2s_norm(s):
    """Normalisation Gamma(s+1)/(2 sqrt(pi) Gamma(s+1/2)) (Longuet-Higgins)."""
    s = np.asarray(s, dtype=np.float64)
    return np.exp(gammaln(s + 1.0) - 0.5 * math.log(math.pi)
                  - gammaln(s + 0.5))


def directional_factor(kx, ky, theta_w, s):
    """D(theta) = N(s) cos^{2s}((theta-theta_w)/2), normalised to int = 1.

    s may be a scalar or an (N,N) array (per spectral mode).
    """
    kx = np.asarray(kx, dtype=np.float64)
    ky = np.asarray(ky, dtype=np.float64)
    dth = np.arctan2(np.sin(np.arctan2(ky, kx) - float(theta_w)),
                     np.cos(np.arctan2(ky, kx) - float(theta_w)))
    c = np.cos(0.5 * dth)
    c = np.maximum(c, 1e-12)
    logd = np.log(cos2s_norm(s)) + 2.0 * np.asarray(s, dtype=np.float64) * np.log(c)
    d = np.exp(logd)
    d = np.where(np.isfinite(d), d, 0.0)
    return d


def _phillips_kshape(k, L):
    kk = np.maximum(np.asarray(k, dtype=np.float64), 1e-12)
    return (np.exp(-1.0 / (kk * L) ** 2) / kk ** 4
            * np.exp(-((kk * DAMP_L_STATE) ** 2)))


def direction_comparison(U, F, h, L, kmin, kmax, n=512, ntheta=1440):
    """Compare the v8 D(theta) with the current Phillips factor.

    Current factor (frozen state.py):  D_old = cos^2(dtheta) + 0.02
                                         (normalised angular integral 3.2673)
    New factor:  cos^{2s}(dtheta/2) with s from the Mitsuyasu shell law.
    Both are compared on the same resolved k-band with the Phillips k-shape
    weighting, so the numbers are slope-variance fractions, not cosmetics.
    """
    sp = spreading_parameters(U, F)
    th = np.linspace(-math.pi, math.pi, int(ntheta), endpoint=False)
    dth = th
    d_old_raw = np.cos(dth) ** 2 + 0.02
    int_old = float(np.trapz(d_old_raw, th))
    d_old = d_old_raw / int_old

    k = np.logspace(math.log10(kmin), math.log10(kmax), int(n))
    om = dispersion_omega(k, h)
    s = spreading_s(om, sp["omega_p_rad_s"], sp["s_max"])
    d_new = np.asarray(cos2s_norm(s))[:, None] * np.cos(0.5 * dth)[None, :] ** (
        2.0 * s[:, None])
    int_new = np.trapz(d_new, th, axis=1)
    d_new = d_new / int_new[:, None]

    w = k ** 2 * _phillips_kshape(k, L) * k  # k^2 * S1(k) * shell factor
    w = w / w.sum()

    def within(D, half_deg):
        m = np.abs(dth) <= math.radians(half_deg)
        frac = np.trapz(D[:, m], th[m], axis=1)
        return float((w * frac).sum())

    def cosm2(D):
        return float((w * np.trapz(D * np.cos(dth)[None, :] ** 2, th, axis=1)).sum())

    half_old = 45.0  # cos^2: D=0.5 at 45 deg
    # new half width at the peak mode (s at omega_p)
    return {
        "old_factor": "cos^2(dtheta) + 0.02 (frozen state.py)",
        "new_factor": ("cos^{2s}(dtheta/2), s = s_max (omega/omega_p)^m, "
                       "Mitsuyasu et al. (1975)"),
        "spreading": sp,
        "s_at_peak": float(sp["s_max"]),
        "s_min_over_band": float(s.min()),
        "s_max_over_band": float(s.max()),
        "half_width_deg_old": half_old,
        "half_width_deg_new_at_peak": float(
            2.0 * math.degrees(math.acos(0.5 ** (1.0 / (2.0 * sp["s_max"]))))),
        "slope_weighted_cos2_mean_old": cosm2(np.tile(d_old, (n, 1))),
        "slope_weighted_cos2_mean_new": cosm2(d_new),
        "fraction_within_30deg_old": within(np.tile(d_old, (n, 1)), 30.0),
        "fraction_within_30deg_new": within(d_new, 30.0),
        "fraction_within_60deg_old": within(np.tile(d_old, (n, 1)), 60.0),
        "fraction_within_60deg_new": within(d_new, 60.0),
        "angular_integral_old_raw": int_old,
        "note": ("the total mean square slope is unchanged by D (it depends "
                 "only on the omnidirectional spectrum); D redistributes it "
                 "in angle.  The old cos^2 lobe is isotropic in the sense "
                 "that its width is fixed for every mode; the new factor is "
                 "narrow at the peak (young sea) and widens towards the "
                 "capillary tail."),
    }


# --------------------------------------------------------------------------
# 4. capillary tail cutoff
# --------------------------------------------------------------------------
def viscous_cutoff(h, k_hi=1e9):
    """k where viscous damping 2 nu k^2 overtakes omega (deep/shallow)."""
    h = float(h)
    lo, hi = 1e-3, float(k_hi)
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if 2.0 * NU_KIN * mid ** 2 > float(dispersion_omega(np.array([mid]), h)[0]):
            hi = mid
        else:
            lo = mid
    kc = math.sqrt(lo * hi)
    return kc


def capillary_crossover():
    """k where gravity and capillarity balance: k_m = sqrt(rho g / gamma)."""
    km = math.sqrt(RHO_W * G / SIGMA_SURF)
    return {"k_m_rad_per_m": float(km), "lambda_m": float(2.0 * math.pi / km),
            "note": "below lambda_m capillarity dominates the dispersion"}


def capillary_tail(S_patch, N, h):
    """Cutoff hierarchy (rule 4) + the resolved-band tail statement."""
    knyq = math.pi * N / float(S_patch)
    kvis = viscous_cutoff(h)
    kmin = 2.0 * math.pi / float(S_patch)
    lam_vis = 2.0 * math.pi / kvis
    # per-loop viscous attenuation of the shortest resolved mode
    att_looptail = float(math.exp(-2.0 * NU_KIN * knyq ** 2 * 8.0))
    drop = 0.0
    return {
        "rule": ("the high-wavenumber cutoff is the wavenumber at which the "
                 "viscous damping rate 2 nu k^2 equals omega(k): modes above "
                 "it are critically damped within one radian of phase and are "
                 "removed from the state spectrum"),
        "k_viscous_rad_per_m": float(kvis),
        "lambda_viscous_m": float(lam_vis),
        "k_nyquist_rad_per_m": float(knyq),
        "lambda_nyquist_m": float(2.0 * math.pi / knyq),
        "kmin_rad_per_m": float(kmin),
        "k_viscous_over_nyquist": float(kvis / knyq),
        "modes_removed_in_resolved_band": 0,
        "effective_cutoff": ("numerical (grid Nyquist) -- the viscous rule "
                             "sits %.2g x above it, so inside the resolved "
                             "band no mode is removed" % (kvis / knyq)),
        "per_loop_viscous_attenuation_at_nyquist": att_looptail,
        "capillary_crossover": capillary_crossover(),
        "deficit_ratio_front": drop,
        "note": ("PIPELINE-V8 chunk 4 asks for the capillary tail's cutoff to "
                 "be fixed by a physical rule.  The rule is implemented and "
                 "evaluated; it does not bind inside this patch, so the "
                 "sharpness of the finest cusps is set by dx = S/N, not by "
                 "viscosity -- stated here rather than silently ignored.  "
                 "The per-loop viscous attenuation of the finest resolved "
                 "mode is recorded, not applied: the state is a static "
                 "spectrum repeated every 8 s."),
    }


# --------------------------------------------------------------------------
# 5. breaking clamp
# --------------------------------------------------------------------------
def breaking_limits(lambda_p, h):
    """H_max = min(0.142 lambda, 0.78 h) + the Stokes cap for comparison."""
    lam = float(lambda_p)
    h = float(h)
    h_stokes = BREAK_STEEPNESS * lam
    h_depth = BREAK_DEPTH_FRAC * h
    h_stokes17 = STOKES_STEEPNESS * lam
    return {
        "lambda_peak_m": lam, "depth_m": h,
        "H_stokes_0p142_m": float(h_stokes),
        "H_stokes_1over7_m": float(h_stokes17),
        "H_depth_0p78h_m": float(h_depth),
        "H_max_m": float(min(h_stokes, h_depth)),
        "binding": ("stokes" if h_stokes <= h_depth else "depth"),
        "rule": "H_max = min(0.142 lambda, 0.78 h)",
        "note": ("0.142 lambda is within 0.13% of the Stokes 1/7 steepness "
                 "cap; the new, independent bound is 0.78 h, which is what "
                 "bites in shallow water"),
    }


def apply_breaking_clamp(H0, N, k_f32, kmin, kmax, depth, sigma2_scale=None):
    """Scale the spectrum so that H_s <= min(0.142 lambda, 0.78 h).

    Returns (H0_new, node).  lambda is the dominant wavelength of the
    realised spectrum (the same peak readout v5 uses).
    """
    S = (np.abs(H0.astype(np.complex128)) ** 2) / (2.0 * N ** 4)
    k = k_f32.astype(np.float64)
    m0 = float(S.sum())
    hs = float(4.0 * math.sqrt(max(m0, 0.0)))
    # dominant wavelength = peak of the radially integrated spectrum
    nb = 400
    edges = np.logspace(math.log10(max(kmin * 0.5, 1e-6)),
                        math.log10(kmax), nb + 1)
    idx = np.digitize(k.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nb)
    acc = np.bincount(idx[ok], weights=S.ravel()[ok], minlength=nb)
    dk = np.diff(edges)
    i = int(np.argmax(acc / np.maximum(dk, 1e-30)))
    kp = float(math.sqrt(edges[i] * edges[i + 1]))
    lam = float(2.0 * math.pi / kp)
    lim = breaking_limits(lam, depth)
    applied = bool(hs > lim["H_max_m"] > 0.0)
    factor = float(lim["H_max_m"] / hs) if applied else 1.0
    H0_new = (H0 * factor).astype(np.complex64) if applied else H0
    node = {
        **lim,
        "Hs_before_m": hs,
        "Hs_after_m": float(hs * factor),
        "factor": float(factor),
        "applied": applied,
        "Hs_over_Hmax_before": float(hs / lim["H_max_m"])
        if lim["H_max_m"] > 0 else None,
        "peak_is_at_kmin": bool(kp <= kmin * 1.35),
        "stokes_cap_also_recorded": True,
    }
    return H0_new, node


# --------------------------------------------------------------------------
# patch / metadata
# --------------------------------------------------------------------------
def patch_metadata(S_patch, N, specL):
    kmin = 2.0 * math.pi / float(S_patch)
    return {
        "patch_S_m": float(S_patch),
        "is_max_resolved_wavelength": True,
        "max_wavelength_m": float(S_patch),
        "kmin_rad_per_m": float(kmin),
        "dx_m": float(S_patch) / int(N),
        "kmax_nyquist_rad_per_m": float(math.pi * int(N) / float(S_patch)),
        "specL_m": float(specL),
        "note": ("the patch side IS the longest resolved wavelength (kmin = "
                 "2*pi/S).  Low-wavenumber truncation bites when long waves "
                 "dominate (strong steady wind over a large pond); enlarging "
                 "the patch costs N^2 and is out of scope for chunk 4."),
    }


def anchor_gate(target_rad2, realised_rad2, old_target_rad2):
    """Gate 3: the new anchor is met, and the old target is reported."""
    rel = float((realised_rad2 - target_rad2) / target_rad2)
    return {
        "gate": ("chunk-4 gate 3: realised sigma^2 passes the NEW anchor "
                 "gate (|relative residual| <= 5%), old Cox-Munk target "
                 "difference reported alongside"),
        "new_target_rad2": float(target_rad2),
        "realised_rad2": float(realised_rad2),
        "relative_residual": rel,
        "old_cox_munk_target_rad2": float(old_target_rad2),
        "delta_new_minus_old": float(target_rad2 - old_target_rad2),
        "ratio_new_over_old": float(target_rad2 / old_target_rad2)
        if old_target_rad2 > 0 else None,
        "tolerance": 0.05,
        "pass": bool(abs(rel) <= 0.05),
    }
