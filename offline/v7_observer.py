"""v7 observer: human-retina post-stage (report.md section 4).

Everything in sections 1-3 of the report produces a radiance field.  This
module is the stage the report's section 4 asks for: the eye that looks at
the display.  It runs AFTER the frozen composite (v5_composite), on the
final linear image, and it is driven by one physical input only -- the
adapting luminance L_a [cd/m^2] -- plus the stated display/viewing geometry
without which pupil diameters and Weber fractions are numbers without a
referent (report section 4, prerequisite).

Nodes (report section 4.1-4.4 / section 5 variable map):

  pupil D(L_a)      Moon & Spencer 1944 (report 4.1, "simpler form").
                    Implemented exactly as published; clamped to the
                    physiological 2-8 mm band.  Watson & Yellott 2012 is
                    provided as a REFERENCE readout only (it needs the
                    adapting field area, which is a viewing-geometry input).
  adaptation gain   Steady-state divisive gain from JND integration:
                      threshold T(L) = a.(L + sqrt(L.L_t))  (Weber T~L high,
                        de Vries-Rose T~sqrt(L) low, report 4.2),
                      response  R(L) = ln(1 + sqrt(L/L_t)),
                      map       out = L_a . R(L)/R(L_a),
                    with L_t = 0.3 cd/m^2 (report 4.2: de Vries-Rose blend
                    below ~0.3).  Neutral point exact by construction:
                    out(L_a) = L_a.  Near L_a (photopic) the map is Weber
                    (log-like, local exponent << 1); for L -> 0 it is
                    proportional to sqrt(L) (Rose asymptote).
                    This is a steady-state gain, NOT adaptation dynamics:
                    dark adaptation takes ~30 min and does not apply to a
                    still (report 4.2, "what cannot be included honestly").
  defocus blur      Uniform defocus only (report 4.3 option (i)): the display
                    is flat, the eye accommodates to the display distance, so
                    there is no honest per-pixel depth blur.  Geometric blur
                    disk diameter d = p.DeltaD [rad] (p = pupil [m],
                    DeltaD = focus error [diopters]); realised as the
                    equal-variance Gaussian sigma = d/4.  DeltaD = 0 returns
                    the input untouched (bit-identity gate).  The labelled
                    artistic pseudo-depth blur (report 4.3 option (ii)) is
                    NOT implemented -- it is art, not physics.
  CSF weighting     Mannos & Sakrison 1974 (report 4.4):
                      H(f) = 2.6.(0.0192+0.114f).exp(-(0.114f)^1.1),
                    f in cycles/degree, used normalised to peak 1.  Used as
                    ANALYSIS (weight the still's spatial-frequency content,
                    report 4.4 honest use (i)) and as an OPTIONAL labelled
                    shaping filter (strength 0 = off, report 4.4 (ii)).
                    Without pixels-per-degree the CSF has no frequency axis
                    (report 4.4 prerequisite); ppd is an explicit input.

Nothing else in this module is a free parameter.  Display scale
(peak_cd_m2), adapting level (L_a, default = image mean), viewing geometry
(ppd), age (WY reference only) and defocus are the recipe inputs.
"""

import numpy as np

# --------------------------------------------------------------------------
# 4.1 pupil
# --------------------------------------------------------------------------

PUPIL_MIN_MM = 2.0
PUPIL_MAX_MM = 8.0

# Moon & Spencer 1944, as quoted in report 4.1 (D in mm, L in cd/m^2).
MS_BASE_MM = 4.9
MS_AMP_MM = 3.0
MS_SLOPE = 0.4
MS_OFFSET = 0.00114

# Report 4.1 quoted values (computed from the formula above): the gate checks
# the implementation against THESE rounded numbers, so a transcription error
# fails loudly instead of comparing the formula to itself.
MS_ANCHORS = (
    (0.01, 6.9),
    (1.0, 4.9),
    (100.0, 2.9),
    (1.0e4, 2.1),
)
MS_ANCHOR_TOL_MM = 0.5

# Watson & Yellott 2012 core (report 4.1).  L = adapting luminance [cd/m^2],
# a = adapting field area [deg^2].  Reference readout only.
WY_C = 7.75
WY_K = 5.75
WY_REF = 846.0
WY_EXP = 0.41
WY_REF_AGE = 28.58


def pupil_diameter_ms(L_cd_m2):
    """Moon & Spencer pupil diameter [mm], clamped to [2, 8].

    D = 4.9 - 3.tanh(0.4.log10(L) - 0.00114).  Scalar or array in, same out.
    L <= 0 (or non-finite) maps to the dark limit 7.9 mm, then clamped.
    Strictly decreasing in L on (0, +inf).
    """
    L = np.asarray(L_cd_m2, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = MS_SLOPE * np.log10(L) - MS_OFFSET
        D = MS_BASE_MM - MS_AMP_MM * np.tanh(x)
    dark_limit = MS_BASE_MM + MS_AMP_MM  # tanh -> -1 as L -> 0
    D = np.where(np.isfinite(D) & (L > 0.0), D, dark_limit)
    D = np.clip(D, PUPIL_MIN_MM, PUPIL_MAX_MM)
    if D.ndim == 0:
        return float(D)
    return D


def pupil_diameter_wy(L_cd_m2, field_area_deg2, age_yr=WY_REF_AGE):
    """Watson & Yellott 2012 unified pupil [mm].  REFERENCE readout, not gated.

    D_SD(L,a) = 7.75 - 5.75.[(L.a/846)^0.41 / ((L.a/846)^0.41 + 2)],
    D_U = D_SD + (y - 28.58).[0.02132 - 0.009562.D_SD].
    The eccentricity/duty factor M(e) in the full form is taken as 1
    (stated simplification); field_area_deg2 must be supplied from the
    viewing geometry -- without it the value has no referent (report 4.1).
    """
    L = float(L_cd_m2)
    a = float(field_area_deg2)
    y = float(age_yr)
    if not (L > 0.0 and a > 0.0):
        raise ValueError("WY needs L > 0 and field area > 0 (got %r, %r)"
                         % (L_cd_m2, field_area_deg2))
    q = (L * a / WY_REF) ** WY_EXP
    d_sd = WY_C - WY_K * (q / (q + 2.0))
    return float(d_sd + (y - WY_REF_AGE) * (0.02132 - 0.009562 * d_sd))


def pupil_gate(L_lo=1.0e-3, L_hi=1.0e5, n=200):
    """Gate: MS pupil strictly decreasing, in [2,8], anchors +-0.5 mm."""
    Ls = np.logspace(np.log10(L_lo), np.log10(L_hi), int(n))
    Ds = pupil_diameter_ms(Ls)
    diffs = np.diff(Ds)
    anchors = []
    for L_ref, D_ref in MS_ANCHORS:
        Dgot = pupil_diameter_ms(L_ref)
        anchors.append({"L_cd_m2": float(L_ref),
                        "D_mm": float(Dgot),
                        "quoted_mm": float(D_ref),
                        "abs_err_mm": float(abs(Dgot - D_ref)),
                        "pass": bool(abs(Dgot - D_ref) <= MS_ANCHOR_TOL_MM)})
    out = {"L_grid": [float(L_lo), float(L_hi), int(n)],
           "D_min_mm": float(Ds.min()), "D_max_mm": float(Ds.max()),
           "strictly_decreasing": bool(np.all(diffs < 0.0)),
           "within_2_8_mm": bool(np.all((Ds >= PUPIL_MIN_MM)
                                        & (Ds <= PUPIL_MAX_MM))),
           "anchors": anchors,
           "anchors_pass": bool(all(a["pass"] for a in anchors))}
    out["pass"] = bool(out["strictly_decreasing"] and out["within_2_8_mm"]
                       and out["anchors_pass"])
    return out


# --------------------------------------------------------------------------
# 4.2 adaptation
# --------------------------------------------------------------------------

# de Vries-Rose / Weber transition (report 4.2: blend below ~0.3 cd/m^2).
ADAPT_TRANSITION_CD_M2 = 0.3


def _adapt_response(L_cd_m2, Lt=ADAPT_TRANSITION_CD_M2):
    """R(L) = ln(1 + sqrt(L/Lt)); R(0) = 0.  JND-integral of the
    Weber/Rose-blended threshold T(L) = a.(L + sqrt(L.Lt))."""
    L = np.asarray(L_cd_m2, dtype=np.float64)
    return np.log1p(np.sqrt(np.clip(L, 0.0, None) / float(Lt)))


def adaptation_map(L_cd_m2, La_cd_m2, Lt=ADAPT_TRANSITION_CD_M2):
    """Steady-state adaptation map.  out = La.R(L)/R(La).

    Exact neutral point out(La) = La; out(0) = 0; monotone increasing.
    Weber-form (log-like) for L >> Lt, sqrt-law for L -> 0.
    """
    La = float(La_cd_m2)
    Lt = float(Lt)
    if not (La > 0.0):
        raise ValueError("adapting luminance La must be > 0 (got %r)"
                         % (La_cd_m2,))
    if not (Lt > 0.0):
        raise ValueError("transition Lt must be > 0 (got %r)" % (Lt,))
    R = _adapt_response(L_cd_m2, Lt)
    R_a = _adapt_response(La, Lt)
    out = La * (R / R_a)
    if np.ndim(out) == 0:
        return float(out)
    return out


def _loglog_exponent(Ls, outs):
    """Local power-law exponent of out vs L (polyfit in log10-log10)."""
    x = np.log10(np.asarray(Ls, dtype=np.float64))
    y = np.log10(np.maximum(np.asarray(outs, dtype=np.float64),
                            1.0e-300))
    return float(np.polyfit(x, y, 1)[0])


def adaptation_gate(La_cd_m2, Lt=ADAPT_TRANSITION_CD_M2):
    """Gate: neutral point, Weber form (photopic), Rose asymptote (scotopic).

    Weber-form band: [La, 4.La] must be compressive log-like
      (exponent in [0.03, 0.35]; pure log at La=200/Lt=0.3 gives ~0.15).
    Rose band: [La.1e-5, La.1e-4] (<< Lt) must be sqrt-like (0.5 +- 0.05).
    """
    La = float(La_cd_m2)
    R_a = _adapt_response(La, Lt)
    neutral = adaptation_map(La, La, Lt)
    hi = np.logspace(np.log10(La), np.log10(4.0 * La), 25)
    lo = np.logspace(np.log10(La * 1.0e-5), np.log10(La * 1.0e-4), 25)
    exp_weber = _loglog_exponent(hi, adaptation_map(hi, La, Lt))
    exp_rose = _loglog_exponent(lo, adaptation_map(lo, La, Lt))
    zero = adaptation_map(0.0, La, Lt)
    mono = np.all(np.diff(adaptation_map(
        np.logspace(-5.0, 5.0, 400), La, Lt)) > 0.0)
    out = {"La_cd_m2": La, "Lt_cd_m2": float(Lt),
           "neutral_out": float(neutral),
           "neutral_exact": bool(neutral == La),
           "out_zero": float(zero),
           "monotone_increasing": bool(mono),
           "weber_band_L": [float(hi[0]), float(hi[-1])],
           "weber_exponent": float(exp_weber),
           "weber_pass": bool(0.03 <= exp_weber <= 0.35),
           "rose_band_L": [float(lo[0]), float(lo[-1])],
           "rose_exponent": float(exp_rose),
           "rose_pass": bool(abs(exp_rose - 0.5) <= 0.05),
           "R_La": float(R_a)}
    out["pass"] = bool(out["neutral_exact"] and out["monotone_increasing"]
                       and out["weber_pass"] and out["rose_pass"]
                       and zero == 0.0)
    return out


# --------------------------------------------------------------------------
# 4.3 defocus
# --------------------------------------------------------------------------

def defocus_sigma_px(pupil_m, delta_D, pixels_per_degree):
    """Equal-variance Gaussian sigma [px] for a focus error delta_D [diopter].

    Geometric disk diameter d = p.DeltaD [rad]; sigma_rad = d/4 (a uniform
    disk of diameter d has 1-D marginal variance d^2/16).  DeltaD = 0 gives
    exactly 0.0 (the caller must bypass filtering at 0 for bit-identity).
    """
    p = float(pupil_m)
    dD = float(delta_D)
    ppd = float(pixels_per_degree)
    if not (p > 0.0):
        raise ValueError("pupil must be > 0 (got %r)" % (pupil_m,))
    if not (dD >= 0.0):
        raise ValueError("defocus must be >= 0 (got %r)" % (delta_D,))
    if not (ppd > 0.0):
        raise ValueError("ppd must be > 0 (got %r)" % (pixels_per_degree,))
    if dD == 0.0:
        return 0.0
    sigma_rad = (p * dD) / 4.0
    return float(sigma_rad * (180.0 / np.pi) * ppd)


def apply_defocus(linear_rgb, sigma_px):
    """Gaussian defocus on linear RGB.  sigma_px <= 0 returns input unchanged
    (same object; the DeltaD = 0 identity gate relies on this bypass)."""
    img = np.asarray(linear_rgb, dtype=np.float64)
    if float(sigma_px) <= 0.0:
        return linear_rgb
    from scipy import ndimage as _ndi
    out = np.empty_like(img)
    for b in range(img.shape[2]):
        out[:, :, b] = _ndi.gaussian_filter(img[:, :, b], float(sigma_px),
                                            mode="reflect")
    return out


def _impulse_width_px(size, sigma_px):
    """Second-moment sigma [px] of the impulse response of apply_defocus."""
    img = np.zeros((int(size), int(size)), dtype=np.float64)
    img[int(size) // 2, int(size) // 2] = 1.0
    from scipy import ndimage as _ndi
    bl = _ndi.gaussian_filter(img, float(sigma_px), mode="reflect")
    n = bl.shape[0]
    ax = np.arange(n) - n // 2
    tot = bl.sum()
    if not (tot > 0.0):
        return 0.0
    marg = bl.sum(axis=0) / tot
    var = float((marg * (ax ** 2)).sum())
    return float(np.sqrt(max(var, 0.0)))


def defocus_gate(pupil_mm, pixels_per_degree, delta_Ds=(0.5, 1.0, 2.0),
                 size=256):
    """Gate: DeltaD = 0 identity; measured width = 4.sigma_meas tracks
    p.DeltaD linearly within +-10% (slope 1 +- 0.1, intercept ~ 0)."""
    p_m = float(pupil_mm) / 1000.0
    ppd = float(pixels_per_degree)
    px_per_rad = ppd * (180.0 / np.pi)
    # identity at zero (small-array bit check)
    probe = np.arange(48, dtype=np.float64).reshape(6, 8, 1).repeat(3, axis=2)
    idem = apply_defocus(probe, defocus_sigma_px(p_m, 0.0, ppd))
    ident = {"max_abs_diff": float(np.max(np.abs(
        np.asarray(idem, dtype=np.float64) - probe))),
             "identical": bool(np.array_equal(np.asarray(idem), probe))}
    rows = []
    for dD in delta_Ds:
        sig = defocus_sigma_px(p_m, float(dD), ppd)
        meas = _impulse_width_px(size, sig)
        w_rad = 4.0 * meas / px_per_rad
        pred = p_m * float(dD)
        rows.append({"delta_D": float(dD), "sigma_px": float(sig),
                     "measured_width_rad": float(w_rad),
                     "predicted_p_dD_rad": float(pred),
                     "ratio": float(w_rad / pred) if pred > 0 else None})
    xs = np.array([r["predicted_p_dD_rad"] for r in rows])
    ys = np.array([r["measured_width_rad"] for r in rows])
    slope, intercept = (float(a) for a in np.polyfit(xs, ys, 1))
    ratios_ok = all(abs(r["ratio"] - 1.0) <= 0.10 for r in rows)
    out = {"pupil_mm": float(pupil_mm), "ppd": float(ppd),
           "identity_D0": ident, "rows": rows, "slope": slope,
           "intercept_rad": float(intercept),
           "slope_pass": bool(abs(slope - 1.0) <= 0.10),
           "intercept_pass": bool(abs(intercept) <= 0.05 * xs.min()),
           "rows_pass": bool(ratios_ok)}
    out["pass"] = bool(ident["identical"] and out["slope_pass"]
                       and out["intercept_pass"] and out["rows_pass"])
    return out


# --------------------------------------------------------------------------
# 4.4 CSF
# --------------------------------------------------------------------------

def csf_mannos_sakrison(f_cpd):
    """Unnormalised Mannos & Sakrison 1974 CSF at f [cycles/degree]."""
    f = np.asarray(f_cpd, dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        H = 2.6 * (0.0192 + 0.114 * f) * np.exp(-((0.114 * f) ** 1.1))
    H = np.where(f > 0.0, H, 0.0)
    if H.ndim == 0:
        return float(H)
    return H


def csf_peak(step=0.001, f_max=60.0):
    """Peak of the MS1974 CSF: (f_peak_cpd, H_peak).  Deterministic scan."""
    f = np.arange(step, f_max + step, step)
    H = csf_mannos_sakrison(f)
    i = int(np.argmax(H))
    return float(f[i]), float(H[i])


def csf_normalized(f_cpd):
    """H(f)/H_peak (peak = 1 by construction)."""
    _, Hp = csf_peak()
    return csf_mannos_sakrison(f_cpd) / Hp


def csf_gate(lo_cpd=6.0, hi_cpd=10.0):
    """Gate: CSF peak in [6, 10] cpd (report 4.4: ~8, analysis band 6-10);
    normalised peak == 1; response ~0 by 60 cpd."""
    fp, Hp = csf_peak()
    h60 = float(csf_normalized(60.0))
    out = {"f_peak_cpd": float(fp), "H_peak": float(Hp),
           "peak_in_band": bool(lo_cpd <= fp <= hi_cpd),
           "normalised_peak": float(csf_normalized(fp)),
           "H_60_over_peak": h60,
           "tail_pass": bool(h60 < 0.05)}
    out["pass"] = bool(out["peak_in_band"] and out["tail_pass"]
                       and out["normalised_peak"] == 1.0)
    return out


def csf_analysis(gray_linear, pixels_per_degree, patch_S_m=None,
                 image_W_px=None):
    """Honest use (i) of the CSF: weight the still's frequency content.

    Input: 2-D linear luminance in [0, 1].  Returns radially averaged
    amplitude, CSF weighting, weighted/unweighted energy, peak response
    frequency, and (when the patch scale is given) the share of the
    CSF-weighted energy inside the rain-ring band 310-630 rad/m
    (report 4.4: check whether the ring band survives perception).
    Frequencies: cycles/px -> cpd via ppd.  DC excluded from energies.
    """
    g = np.asarray(gray_linear, dtype=np.float64)
    Hh, Ww = g.shape
    ppd = float(pixels_per_degree)
    F = np.fft.rfft2(g - g.mean())
    amp = np.abs(F)
    fy = np.fft.fftfreq(Hh)[:, None]
    fx = np.fft.rfftfreq(Ww)[None, :]
    cyc_px = np.sqrt(fx ** 2 + fy ** 2)
    cpd = cyc_px * ppd
    H = csf_normalized(np.maximum(cpd, 1.0e-6))
    wenergy = float((np.abs(F * H) ** 2).sum())
    energy = float((np.abs(F) ** 2).sum())
    nbins = 240
    edges = np.linspace(0.0, float(cpd.max()), nbins + 1)
    idx = np.clip(np.digitize(cpd.ravel(), edges) - 1, 0, nbins - 1)
    num = np.bincount(idx, weights=(np.abs(F * H) ** 2).ravel(),
                      minlength=nbins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    peak_resp = float(centres[int(np.argmax(num[1:])) + 1]) if nbins > 1 else 0.0
    ring = None
    if patch_S_m is not None and image_W_px is not None:
        m_per_deg = (float(patch_S_m) / float(image_W_px)) * ppd
        k = cpd / max(m_per_deg, 1.0e-12) * (2.0 * np.pi)  # rad/m
        in_band = (k >= 310.0) & (k <= 630.0)
        ring = {"m_per_deg": float(m_per_deg),
                "weighted_share_in_310_630_rad_per_m":
                    float((np.abs(F * H)[in_band] ** 2).sum()
                          / max(wenergy, 1.0e-300))}
    return {"ppd": ppd, "image_hw": [int(Hh), int(Ww)],
            "energy_unweighted": energy,
            "energy_csf_weighted": wenergy,
            "weighted_over_unweighted":
                float(wenergy / max(energy, 1.0e-300)),
            "peak_response_cpd": peak_resp,
            "ring_band": ring}


def apply_csf(linear_gray, pixels_per_degree, strength=0.0):
    """OPTIONAL labelled shaping (report 4.4 (ii)): blend the amplitude
    spectrum toward the CSF-weighted one.  strength = 0 returns the input
    unchanged (same object; identity relied on by the OFF gate path)."""
    s = float(strength)
    if not (0.0 <= s <= 1.0):
        raise ValueError("csf strength must be in [0,1] (got %r)"
                         % (strength,))
    if s == 0.0:
        return linear_gray
    g = np.asarray(linear_gray, dtype=np.float64)
    Hh, Ww = g.shape
    F = np.fft.rfft2(g)
    fy = np.fft.fftfreq(Hh)[:, None]
    fx = np.fft.rfftfreq(Ww)[None, :]
    H = csf_normalized(np.maximum(np.sqrt(fx ** 2 + fy ** 2)
                                  * float(pixels_per_degree), 1.0e-6))
    G = np.fft.irfft2(F * ((1.0 - s) + s * H), s=g.shape)
    return G


# --------------------------------------------------------------------------
# full post-stage on a still
# --------------------------------------------------------------------------

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)


def observer_still(linear_rgb, peak_cd_m2=200.0, La_cd_m2=None,
                   pupil_mm=None, delta_D=0.0, pixels_per_degree=60.0,
                   csf_strength=0.0):
    """Apply the observer post-stage to linear RGB in [0, 1].

    Display scale: L_display = Y_rec709 . peak_cd_m2 (peak white = peak).
    L_a default = mean displayed luminance (stated steady-state choice).
    pupil default = MS D(L_a).  Returns dict(out=linear RGB [0,1 clipped],
    Y_in, Y_out, L_a used, pupil, sigma_px).  Pure function; OFF bypass
    (bit-identity) is the caller's job -- see run_v7.
    """
    img = np.asarray(linear_rgb, dtype=np.float64)
    peak = float(peak_cd_m2)
    ppd = float(pixels_per_degree)
    if not (peak > 0.0):
        raise ValueError("peak_cd_m2 must be > 0")
    Y = (img * REC709).sum(axis=-1)
    L = Y * peak
    La = float(np.mean(L)) if La_cd_m2 is None else float(La_cd_m2)
    if not (La > 0.0):
        raise ValueError("adapting luminance must be > 0")
    p_mm = pupil_diameter_ms(La) if pupil_mm is None else float(pupil_mm)
    sig = defocus_sigma_px(p_mm / 1000.0, float(delta_D), ppd)
    work = apply_defocus(img, sig)
    Yw = (np.asarray(work) * REC709).sum(axis=-1)
    Lw = Yw * peak
    Yp = adaptation_map(np.clip(Lw, 0.0, None), La) / peak
    with np.errstate(divide="ignore", invalid="ignore"):
        gain = np.where(Yw > 1.0e-12, Yp / np.maximum(Yw, 1.0e-12), 0.0)
    out = np.asarray(work) * gain[:, :, None]
    if float(csf_strength) > 0.0:
        oy = (np.asarray(out) * REC709).sum(axis=-1)
        shaped = apply_csf(np.clip(oy, 0.0, None), ppd,
                           strength=float(csf_strength))
        with np.errstate(divide="ignore", invalid="ignore"):
            g2 = np.where(oy > 1.0e-12,
                          shaped / np.maximum(oy, 1.0e-12), 0.0)
        out = np.asarray(out) * np.clip(g2, 0.0, None)[:, :, None]
    return {"out": np.clip(out, 0.0, 1.0), "Y_in": Y, "L_a": La,
            "pupil_mm": float(p_mm), "sigma_px": float(sig)}
