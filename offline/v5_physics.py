"""v5 physics nodes (PIPELINE-V5 chunk 1, report.md sections 1, 2, 5).

Every number here is a *derived* quantity from a published formula in
research-v3-variables/report.md.  No aesthetic slider is introduced.

Sections implemented
--------------------
1.1  Significant wave height H_s = 4*sqrt(m0) (readout) and the height-driven
     mode (scalar multiplier on the spectrum) with the Stokes cap H/lambda<=1/7.
1.2  "Calmness" = mean square slope sigma^2, Cox-Munk reference and the
     anisotropic glitter slope variances (sigma_a^2, sigma_c^2).
1.3  Smith G1 masking-shadowing, statistical (Beckmann rational
     approximation, Walter et al. 2007) -- explicitly an approximation.
2.1  Bottom albedo rho_b as per-channel substrate presets (broadband
     literature ranges, report section 2.1).
2.2  Beer-Lambert c(lambda)=a(lambda)+b(lambda) with the Pope & Fry pure-water
     absorption table at the pipeline's channel wavelengths ~650/550/470 nm.
2.3  Turbidity scattering blur sigma_blur ~ sqrt(b*z).

All calibrations that are NOT derivable from the report are named
``*_CALIBRATION`` and documented in the returned dicts / numbers-v5.json.
"""
import math

import numpy as np

# --------------------------------------------------------------------------
# Channel model
# --------------------------------------------------------------------------
# Effective channel wavelengths of the pipeline's IOR dispersion set
# (R 1.33132 / G 1.33434 / B 1.33925).  report section 2.2: "use the R/G/B
# effective wavelengths implied by its dispersion set, ~650/550/470 nm as an
# approximation -- and say so".
CHANNEL_WL_NM = (650.0, 550.0, 470.0)
CHANNEL_WL_NOTE = ("effective channel wavelengths ~650/550/470 nm "
                   "(approximation stated in report section 2.2)")

# Pope & Fry (1997) pure-water absorption, NASA Ocean Optics Protocols
# Table 1.1 / Kirk Table 3.1 (report section 2.2 table).
POPE_FRY_NM = (400.0, 420.0, 440.0, 500.0, 550.0, 600.0, 650.0, 700.0)
POPE_FRY_AW = (0.0066, 0.0045, 0.0064, 0.0204, 0.0565, 0.2224, 0.340, 0.624)

# Pure-water scattering: report section 2.2 gives "order 0.001-0.01 /m,
# Rayleigh-like, ~lambda^-4.32" without a table.  One documented anchor at
# 550 nm is used, everything else follows the published exponent.
BW_550 = 0.0018          # 1/m  (documented order-of-magnitude anchor)
BW_EXP = 4.32            # Rayleigh-like exponent (report section 2.2)

# Turbidity node (0..1, existing recipe field water.turbidity) -> added beam
# attenuation.  The report states turbidity is "the addition da+db on top of
# the pure-water floors" but gives no scale, so the scale is a NAMED
# calibration constant (not derived).
TURBID_DC_MAX = 0.30     # 1/m at turbidity = 1.0
TURBID_SCATTER_FRAC = 0.7  # db/(da+db); suspended particles scatter mostly

# Bottom albedo presets (broadband, report section 2.1):
#   sand 0.24-0.30 (bright carbonate, rises toward red)
#   rock 0.10-0.25 (spectrally flatter, less when wet/submerged)
#   vegetation ~0.12 (seagrass; below 0.10 across 500-600 nm, green peak,
#                     red absorption)
#   mud ~0.05 (dark organic sediment, no canonical table)
# Values are per-channel spectral albedos; the *broadband* albedo quoted in
# the report is the rec709-weighted mean of these (recorded in numbers).
SUBSTRATE_ALBEDO = {
    "sand":       (0.300, 0.270, 0.220),
    "rock":       (0.190, 0.185, 0.175),
    "vegetation": (0.070, 0.130, 0.045),
    "mud":        (0.055, 0.050, 0.042),
}
# The photo *is* the bottom albedo in this pipeline, so a substrate preset is
# a *relative* spectral gain.  Sand is the reference substrate (gain 1.0), so
# "sand" and a legacy recipe (no bottom block) leave the photo untouched.
SUBSTRATE_REFERENCE = "sand"

# Scattering blur: sigma_blur ~ sqrt(b*z) (report section 2.3).  The
# proportionality constant absorbs the forward-scatter angle and the caustic
# pattern scale; it is a NAMED calibration constant.
BLUR_K_PX = 40.0         # px per sqrt(1/m * m) at the master grid 2048^2

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)

STOKES_STEEPNESS = 1.0 / 7.0   # H/lambda limit (report section 1.1)


def _interp_table(x, xs, ys):
    return float(np.interp(float(x), np.asarray(xs, np.float64),
                           np.asarray(ys, np.float64)))


def pure_water_absorption(wl_nm):
    """Pope & Fry (1997) a_w(lambda), linearly interpolated. Returns (3,)."""
    return np.array([_interp_table(w, POPE_FRY_NM, POPE_FRY_AW)
                     for w in wl_nm], dtype=np.float64)


def pure_water_scattering(wl_nm):
    """b_w(lambda) = BW_550 * (550/lambda)^4.32. Returns (3,)."""
    return np.array([BW_550 * (550.0 / float(w)) ** BW_EXP
                     for w in wl_nm], dtype=np.float64)


def turbidity_terms(turbid):
    """Turbidity node -> (da, db) per channel (spectrally flat, report 2.2)."""
    dc = TURBID_DC_MAX * max(0.0, min(1.0, float(turbid)))
    db = TURBID_SCATTER_FRAC * dc
    da = (1.0 - TURBID_SCATTER_FRAC) * dc
    return np.full(3, da, np.float64), np.full(3, db, np.float64), dc


def water_attenuation(turbid):
    """c(lambda) = a_w + b_w + da + db at the three channel wavelengths."""
    a_w = pure_water_absorption(CHANNEL_WL_NM)
    b_w = pure_water_scattering(CHANNEL_WL_NM)
    da, db, dc = turbidity_terms(turbid)
    return {
        "a_w": a_w, "b_w": b_w, "da": da, "db": db,
        "a": a_w + da, "b": b_w + db, "c": a_w + da + b_w + db,
        "dc_total": float(dc),
        "wl_nm": list(CHANNEL_WL_NM),
        "note": CHANNEL_WL_NOTE,
    }


def refraction_angle_water(sun_alt_deg, n_water=1.333):
    """In-water sun zenith angle theta_w after Snell refraction (<=48.6 deg).

    ``sun_alt_deg`` is the altitude above the horizon, so the incidence angle
    from the surface normal is theta_i = 90 deg - altitude.
    """
    alt = math.radians(max(0.0, min(90.0, float(sun_alt_deg))))
    theta_i = 0.5 * math.pi - alt
    s = math.sin(theta_i) / float(n_water)
    s = max(-1.0, min(1.0, s))
    return math.asin(s)


def beer_lambert(depth_m, sun_alt_deg, turbid, view_zenith_deg=0.0,
                 n_water=1.333):
    """Two-way Beer-Lambert factors for the direct (caustic-carrying) beam.

    down-path s_down = z / cos(theta_w)  (theta_w = refracted sun angle)
    up-path   s_up   = z / cos(theta_v)  (theta_v = view zenith; the offline
                     composite is a rectified top view, theta_v = 0)
    Returns dict with per-channel factors (3,) and the path lengths.
    """
    z = max(0.0, float(depth_m))
    th_w = refraction_angle_water(sun_alt_deg, n_water)
    th_v = math.radians(max(0.0, min(48.6, float(view_zenith_deg))))
    s_down = z / max(math.cos(th_w), 1e-6)
    s_up = z / max(math.cos(th_v), 1e-6)
    att = water_attenuation(turbid)
    down = np.exp(-att["c"] * s_down)
    up = np.exp(-att["c"] * s_up)
    return {
        "down": down, "up": up, "two_way": down * up,
        "s_down_m": float(s_down), "s_up_m": float(s_up),
        "theta_w_deg": float(math.degrees(th_w)),
        "theta_v_deg": float(math.degrees(th_v)),
        "c": att["c"], "a": att["a"], "b": att["b"],
        "turbid": float(turbid),
    }


def beer_lambert_gate_example():
    """PIPELINE-V5 gate: clear water, z=2 m, vertical sun, two-way.

    550 nm ~ 0.80, 650 nm ~ 0.26 (+-0.02).  Clear water = turbidity 0,
    vertical sun = altitude 90 deg -> theta_w = 0 -> s_down = s_up = 2 m.
    """
    r = beer_lambert(2.0, 90.0, 0.0)
    tw = r["two_way"]
    return {
        "condition": {"depth_m": 2.0, "sun_alt_deg": 90.0, "turbidity": 0.0,
                      "view_zenith_deg": 0.0, "water": "clear"},
        "s_down_m": r["s_down_m"], "s_up_m": r["s_up_m"],
        "c_rgb": [float(x) for x in r["c"]],
        "two_way_rgb": [float(x) for x in tw],
        "two_way_550nm": float(tw[1]), "two_way_650nm": float(tw[0]),
        "expected_550nm": 0.80, "expected_650nm": 0.26, "tol": 0.02,
        "pass_550nm": bool(abs(float(tw[1]) - 0.80) <= 0.02),
        "pass_650nm": bool(abs(float(tw[0]) - 0.26) <= 0.02),
    }


# --------------------------------------------------------------------------
# 1.1 / 1.2  wave statistics
# --------------------------------------------------------------------------
def cox_munk_sigma2(wind_ms):
    """Cox-Munk mean square slope: 0.003 + 5.12e-3 * U  (rad^2)."""
    return 0.003 + 5.12e-3 * float(wind_ms)


def cox_munk_anisotropy(wind_ms):
    """Anisotropic Cox-Munk slope variances (along / cross wind), rad^2.

    sigma_a^2 = 3.16e-3 U, sigma_c^2 = 0.003 + 1.92e-3 U.
    NOTE: their sum is 0.003 + 5.08e-3 U, while the report's total is
    0.003 + 5.12e-3 U (the report's own rounding); the total formula is used
    as the calibration target and the discrepancy is recorded.
    """
    u = float(wind_ms)
    sa = 3.16e-3 * u
    sc = 0.003 + 1.92e-3 * u
    return {"sigma_a2": float(sa), "sigma_c2": float(sc),
            "sum": float(sa + sc),
            "total_formula": float(cox_munk_sigma2(u)),
            "sum_minus_total": float(sa + sc - cox_munk_sigma2(u))}


def stokes_cap_hs(lambda_p_m):
    """Stokes limiting significant height for a dominant wavelength."""
    return float(STOKES_STEEPNESS * float(lambda_p_m))


# --------------------------------------------------------------------------
# 1.3  Smith G1 (statistical masking-shadowing)
# --------------------------------------------------------------------------
def smith_g1_beckmann(alpha, theta_deg):
    """Smith G1 = 1/(1+Lambda) with the Walter et al. (2007) rational fit.

    alpha = rms surface slope (dimensionless), theta = incidence angle from
    the surface normal (degrees).  a = 1/(alpha*tan theta);
    Lambda ~ (1 - 1.259a + 0.396a^2)/(3.535a + 2.181a^2) for a<1.6 else 0.
    This is the STATISTICAL approximation, not a ray-exact visibility test
    (report section 1.3; do not combine the two).
    """
    th = math.radians(max(0.0, min(90.0, float(theta_deg))))
    tan_t = math.tan(th)
    if alpha <= 0.0 or tan_t <= 1e-9:
        return {"G1": 1.0, "Lambda": 0.0, "a": float("inf"),
                "theta_deg": float(theta_deg), "alpha": float(alpha),
                "regime": "Lambda=0 (a>=1.6 or degenerate)"}
    a = 1.0 / (float(alpha) * tan_t)
    if a >= 1.6:
        lam = 0.0
    else:
        lam = ((1.0 - 1.259 * a + 0.396 * a * a)
               / (3.535 * a + 2.181 * a * a))
    return {"G1": float(1.0 / (1.0 + lam)), "Lambda": float(lam),
            "a": float(a), "theta_deg": float(theta_deg),
            "alpha": float(alpha),
            "regime": ("Lambda=0 (a>=1.6)" if a >= 1.6 else "rational fit")}


# --------------------------------------------------------------------------
# 2.1  bottom albedo presets
# --------------------------------------------------------------------------
def substrate_albedo(name):
    """Per-channel spectral albedo of a substrate preset, or None (legacy)."""
    if name is None:
        return None
    key = str(name).strip().lower()
    if key not in SUBSTRATE_ALBEDO:
        raise ValueError("unknown bottom.substrate %r (expected one of %s)"
                         % (name, sorted(SUBSTRATE_ALBEDO)))
    return np.array(SUBSTRATE_ALBEDO[key], dtype=np.float64)


def substrate_gain(name, albedo_scale=1.0):
    """Channel gain applied to the photo's linear values (direct term only).

    gain = rho_b(lambda) / rho_sand(lambda) * albedoScale.  Legacy recipes
    (substrate None) return an identity gain, i.e. the old behaviour.
    """
    rho = substrate_albedo(name)
    if rho is None:
        g = np.ones(3, np.float64)
        ref = None
    else:
        ref = substrate_albedo(SUBSTRATE_REFERENCE)
        g = rho / ref
    g = g * float(albedo_scale)
    return {
        "substrate": name, "albedo_scale": float(albedo_scale),
        "rho_rgb": (None if rho is None else [float(x) for x in rho]),
        "rho_broadband_rec709": (None if rho is None else
                                 float((rho * REC709).sum())),
        "reference": SUBSTRATE_REFERENCE,
        "reference_rho_rgb": (None if ref is None else
                              [float(x) for x in ref]),
        "gain_rgb": [float(x) for x in g],
        "applies_to": ("photo linear values on the E_base,dir * C term "
                       "(direct-bottom term only, report section 2.1)"),
    }


# --------------------------------------------------------------------------
# 2.3  turbidity scattering blur
# --------------------------------------------------------------------------
def scattering_blur_sigma_px(att, depth_m, master_grid=2048):
    """sigma_blur ~ sqrt(b*z); b = scattering coefficient (green channel).

    Returns the Gaussian sigma in master-grid pixels at ``master_grid``.
    BLUR_K_PX is a NAMED calibration constant (report section 2.3 gives the
    proportionality only).
    """
    b_sca = float(att["b"][1])  # green band, achromatic master C
    bz = max(0.0, b_sca * float(depth_m))
    sig = BLUR_K_PX * math.sqrt(bz) * (float(master_grid) / 2048.0)
    return {"sigma_px": float(sig), "b_scatter_green": b_sca,
            "b_z": float(bz), "k_px_calibration": BLUR_K_PX,
            "form": "sigma_px = k_px * sqrt(b_sca*z)"}


# --------------------------------------------------------------------------
# 1.2  anisotropic glitter slope distribution
# --------------------------------------------------------------------------
def glitter_slope_sigmas(wind_ms, wind_dir_rad=0.0):
    """Anisotropic Cox-Munk slope variances for glitter sampling.

    The isotropic case is sigma_a = sigma_c (special case, as in report 1.2).
    """
    an = cox_munk_anisotropy(wind_ms)
    return {"sigma_a2": an["sigma_a2"], "sigma_c2": an["sigma_c2"],
            "sigma_a": math.sqrt(an["sigma_a2"]),
            "sigma_c": math.sqrt(an["sigma_c2"]),
            "wind_dir_rad": float(wind_dir_rad),
            "isotropic_special_case": {
                "sigma2": 0.5 * (an["sigma_a2"] + an["sigma_c2"]),
                "note": "sigma_a^2 = sigma_c^2 = sigma^2/2"}}
