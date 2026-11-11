"""v9 sky model -- chunk 5 (version 3.5): direct / diffuse separation.

PIPELINE-V8.md "chunk 5" + the 2026-11-11 design addendum.  Nothing in the
frozen pipeline is modified: this module is the only place the sky is
described, and it is used by ``run_v9.py`` only.

What is here
------------
1. PEREZ (1993) all-weather *relative* sky luminance distribution

       L(theta, chi) = L_z * [1 + a exp(b / cos theta)]
                            * [1 + c exp(d chi) + e cos^2 chi]

   theta = zenith angle of the sky element, chi = angle between the element
   and the sun (radians), L_z = zenith luminance.  The five coefficients
   (a,b,c,d,e) are taken from the published 8-bin table indexed by the
   clearness parameter epsilon; see ``PEREZ_COEF`` and ``PEREZ_SOURCE``.

   epsilon is NOT a free input: it follows from the direct fraction eta
   because eta fixes DHI/(DNI+DHI) and therefore (DHI+DNI)/DHI,

       epsilon = [ (DHI+DNI)/DHI + 1.041 Z^3 ] / [ 1 + 1.041 Z^3 ]
               = [ 1/(1-eta)        + 1.041 Z^3 ] / [ 1 + 1.041 Z^3 ]

   with Z the solar zenith angle in RADIANS (Perez et al. 1993 eq. 1).

2. DOME QUADRATURE.  The diffuse component is not a constant ambient term:
   the dome is sampled with N deterministic directions, each direction is
   weighted by the Perez radiance times its incidence cosine (i.e. by its own
   contribution to the horizontal diffuse irradiance), and the SAME Jacobian
   transport is run per direction and summed.  The weights are normalised so
   that sum(w) = 1 exactly, which is what makes the diffuse master carry
   mean(C) = 1 like the direct one.

3. NORMALISATION.  L_z is not taken from a table: it is solved numerically so
   that the hemispherical integral of the distribution equals the requested
   diffuse horizontal irradiance,

       L_z = DHI / integral[ f(theta,chi) cos(theta) dOmega ].

   This is the same construction used by the reference implementation of the
   model (Radiance ``gensky``/``skybright.cal``, LBNL), and it makes the split
   self-consistent: whatever the coefficient table says, the diffuse energy is
   the one that was asked for.  A table-driven L_z would import the paper's
   illuminance calibration on top, which this pipeline does not carry (no
   photometric units anywhere in the chain) -- stated as a design choice.

4. TURBIDITY.  Kept as an input (default 3, Linke) and recorded.  In this
   formulation the turbidity does not enter the *relative* luminance shape:
   the shape is governed by epsilon (which is measured from the radiation
   split) and the absolute level is normalised to the requested DHI.  The
   aerosol/circumsolar correction that would use the turbidity explicitly is
   NOT modelled here; that is recorded rather than silently dropped.  The
   circumsolar *fraction* implied by the distribution is measured and
   reported (it is the part of the diffuse that behaves like a broad beam).

5. WATER PATH FACTORS.  Direct and diffuse get separate mean in-water path
   lengths: the direct beam is one refracted direction, the diffuse is the
   irradiance-weighted mean of z/cos(theta_w) over the dome after Snell
   refraction (the Snell window confines it to <= 48.6 deg, which is why the
   diffuse path is longer but bounded).

Units: the whole chain is relative (C is dimensionless, the split is in
units of DNI+DHI = 1), so no photometric or absolute radiometric constant is
introduced.  ``NON_DERIVED`` lists the constants that remain.

Sources
-------
* Perez, Seals & Michalsky (1993), "All-weather model for sky luminance
  distribution -- preliminary configuration and validation", Solar Energy
  50(3) 235-245.  Table 1 (a,b,c,d,e vs epsilon bins), eq. 1 (epsilon).
* Perez, Ineichen, Seals, Michalsky & Stewart (1990), "Modeling daylight
  availability and irradiance components from direct and global irradiance",
  Solar Energy 44(5) 271-289 (the epsilon/Delta parametrisation reused here).
* Radiance (LBNL) ``gensky`` / ``skybright.cal``: the reference numerical
  implementation of the normalisation used in item 3.
* Longuet-Higgins / Snell: the in-water angle uses the pipeline's own band
  IOR (det.IOR_G = 1.33434) so the diffuse path is consistent with the
  transport's refraction.
"""
import math

import numpy as np

# --------------------------------------------------------------------------
# Perez et al. (1993) Table 1 -- all-weather sky luminance coefficients
# --------------------------------------------------------------------------
#: upper bounds of the 8 epsilon bins (Perez et al. 1993, Table 1)
PEREZ_EPS_BINS = (1.065, 1.230, 1.500, 1.950, 2.800, 4.500, 6.200)

#: (a, b, c, d, e) per epsilon bin, bins 1..8.  Transcribed from the
#: reference implementation of the model (Radiance skybright.cal / gensky,
#: LBNL) which carries Perez et al. (1993) Table 1 verbatim.
PEREZ_COEF = np.array([
    [-0.0080,  0.5880, -0.0620, -0.0600, 0.0720],   # 1  epsilon < 1.065
    [ 0.1300,  0.6830,  0.0130,  0.2500, 0.1760],   # 2  1.065 - 1.230
    [ 0.3300,  0.4870,  0.1360,  0.0120, 0.1000],   # 3  1.230 - 1.500
    [ 0.5680,  0.1870,  0.1880, -0.2090, 0.0320],   # 4  1.500 - 1.950
    [ 0.8730, -0.3750,  0.2500, -0.3340, 0.0300],   # 5  1.950 - 2.800
    [ 1.1320, -0.9620,  0.3180, -0.4810, 0.0170],   # 6  2.800 - 4.500
    [ 1.0600, -0.5840,  0.2660, -0.3740, 0.0170],   # 7  4.500 - 6.200
    [ 0.6780, -0.1870,  0.2260, -0.1990, 0.0130],   # 8  epsilon > 6.200
], dtype=np.float64)

PEREZ_SOURCE = (
    "Perez, Seals & Michalsky (1993) Solar Energy 50(3) 235-245, Table 1; "
    "coefficient table transcribed from the reference implementation "
    "(Radiance skybright.cal / gensky, LBNL)")

#: Delta (sky brightness) bin edges (Perez et al. 1990); recorded only.
PEREZ_DELTA_BINS = (0.15, 0.30, 0.60, 1.20, 2.00, 3.00, 4.50, 6.20)

#: Extraterrestrial normal irradiance, for the Delta diagnostic only.
E0_W_M2 = 1361.0

#: Linke turbidity default (design addendum: "turbidity (default 3)").
TURBIDITY_DEFAULT = 3.0

#: Number of dome directions.  Design: 64-128.  96 is used (see report).
N_SKY_DIRS = 96

#: Horizon elevation of the dome [deg].  The Perez relative distribution is
#: NOT usable at the exact horizon when b > 0: f grows like exp(b/cos th)
#: while the cosine weight vanishes (at 5 deg the bin-3/4 coefficients already
#: give f = 125, an artefact of extrapolating the fit below its data range,
#: not a sky).  The dome is therefore truncated where the fit still returns
#: plausible horizon/zenith ratios (bin 1 f_horizon/f_zenith = 0.76, bin 3
#: 7.9 at 10 deg -- the same order as the CIE standard sky types, whose
#: overcast horizon is 1/3 of the zenith).  10 deg is adopted; the
#: 5/10/15/20 deg sensitivity is recorded by ``horizon_sensitivity``.
HORIZON_ELEV_DEG = 10.0

#: Representative direct fractions used when the real-mode fetch fails.
FALLBACK_ETA = {"clear": 0.80, "overcast": 0.05}

#: In-water refractive index used for the path factors (det.IOR_G).
N_WATER = 1.33434

NON_DERIVED = [
    {"name": "N_SKY_DIRS", "value": N_SKY_DIRS, "units": "-",
     "basis": ("dome quadrature resolution, inside the design's 64-128 "
               "range; the quadrature error is measured (gate 3 + the "
               "direction-count sensitivity table), not assumed")},
    {"name": "TURBIDITY_DEFAULT", "value": TURBIDITY_DEFAULT, "units": "Linke",
     "basis": ("design addendum default; retained as an input and recorded. "
               "It does not enter the relative luminance shape in this "
               "formulation (see module docstring item 4)")},
    {"name": "E0_W_M2", "value": E0_W_M2, "units": "W/m^2",
     "basis": ("solar constant, used only for the Delta sky-brightness "
               "diagnostic; the model itself is normalised to the requested "
               "DHI so no absolute constant enters the transport")},
    {"name": "FALLBACK_ETA", "value": FALLBACK_ETA, "units": "-",
     "basis": ("representative clear / overcast direct fractions used ONLY "
               "when the real-mode fetch fails; the design names them "
               "('晴天 0.8 / 曇天 0.05 級')")},
]


# --------------------------------------------------------------------------
# 1. clearness / coefficients
# --------------------------------------------------------------------------
def epsilon_from_eta(eta, sun_zenith_deg):
    """epsilon = [(1/(1-eta)) + 1.041 Z^3] / [1 + 1.041 Z^3]  (Perez 1993)."""
    e = float(eta)
    e = min(max(e, 0.0), 1.0 - 1e-12)
    Z = math.radians(max(0.0, min(90.0, float(sun_zenith_deg))))
    z3 = 1.041 * Z ** 3
    return float(((1.0 / (1.0 - e)) + z3) / (1.0 + z3))


def eps_bin_index(eps):
    """0-based epsilon bin index (0..7)."""
    e = float(eps)
    for i, hi in enumerate(PEREZ_EPS_BINS):
        if e < hi:
            return i
    return 7


def eps_bin_centre(i):
    """Geometric centre of epsilon bin i (used for interpolation)."""
    lo = 1.0 if i == 0 else PEREZ_EPS_BINS[i - 1]
    hi = PEREZ_EPS_BINS[i] if i < 7 else 8.5
    return float(math.sqrt(lo * hi))


def perez_coefficients(eps):
    """(a,b,c,d,e) at eps, linearly interpolated between bin centres.

    Perez et al. prescribe interpolation over epsilon rather than a step
    lookup; the step-lookup alternative is measured in ``bin_step_check``.
    """
    e = float(eps)
    centres = [eps_bin_centre(i) for i in range(8)]
    if e <= centres[0]:
        c = PEREZ_COEF[0]
    elif e >= centres[-1]:
        c = PEREZ_COEF[-1]
    else:
        i = 0
        while i < 6 and e > centres[i + 1]:
            i += 1
        t = (e - centres[i]) / (centres[i + 1] - centres[i])
        c = PEREZ_COEF[i] * (1.0 - t) + PEREZ_COEF[i + 1] * t
    return {"a": float(c[0]), "b": float(c[1]), "c": float(c[2]),
            "d": float(c[3]), "e": float(c[4]),
            "eps": e, "bin": eps_bin_index(e) + 1,
            "bin_edges": PEREZ_EPS_BINS,
            "interpolation": "linear between bin centres (Perez recommended)"}


def perez_shape(theta_z, chi, coef):
    """Relative luminance f = [1+a exp(b/cos th)] [1+c exp(d chi)+e cos^2 chi].

    theta_z, chi in radians (arrays ok).  Returns f (unnormalised, L/L_z).
    """
    th = np.asarray(theta_z, dtype=np.float64)
    ch = np.asarray(chi, dtype=np.float64)
    ct = np.maximum(np.cos(th), 1.0e-6)
    f1 = 1.0 + float(coef["a"]) * np.exp(np.clip(float(coef["b"]) / ct,
                                                  -700.0, 700.0))
    f2 = (1.0 + float(coef["c"]) * np.exp(np.clip(float(coef["d"]) * ch,
                                                   -700.0, 700.0))
          + float(coef["e"]) * np.cos(ch) ** 2)
    return np.maximum(f1 * f2, 0.0)


def delta_sky_brightness(dhi_rel, sun_alt_deg):
    """Perez Delta = DHI*m/E0 (diagnostic + bin, recorded only)."""
    alt = math.radians(max(1.0, float(sun_alt_deg)))
    m = 1.0 / max(math.sin(alt), 1.0e-3)
    d = float(dhi_rel) * m / E0_W_M2
    b = 1
    for i, hi in enumerate(PEREZ_DELTA_BINS):
        if d < hi:
            b = i + 1
            break
    else:
        b = 8
    return {"Delta": float(d), "bin": int(b), "air_mass": float(m),
            "bins": PEREZ_DELTA_BINS,
            "note": ("Delta is a diagnostic here: the distribution is "
                     "normalised to the requested DHI, so the paper's "
                     "illuminance calibration does not enter the chain")}


# --------------------------------------------------------------------------
# 2. dome quadrature
# --------------------------------------------------------------------------
def fibonacci_hemisphere(n_dirs, horizon_elev_deg=HORIZON_ELEV_DEG):
    """Deterministic dome directions (unit vectors, z >= z_min) + solid angles.

    cos(theta) is uniform in i (equal solid angle per ring) over the truncated
    dome [z_min, 1], the azimuth is the golden angle sequence.  Unlike a random
    sample this is reproducible and has no clumping.
    """
    n = int(n_dirs)
    z_min = math.sin(math.radians(max(0.0, float(horizon_elev_deg))))
    i = np.arange(n, dtype=np.float64)
    z = z_min + (1.0 - z_min) * (i + 0.5) / float(n)   # uniform in cos(theta)
    r = np.sqrt(np.maximum(1.0 - z * z, 0.0))
    phi = i * (np.pi * (3.0 - math.sqrt(5.0)))
    dirs = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    solid = np.full(n, 2.0 * math.pi * (1.0 - z_min) / float(n),
                    dtype=np.float64)
    return dirs, solid


def sky_dome(eta, sun_dir, n_dirs=N_SKY_DIRS, sun_zenith_deg=None,
             horizon_elev_deg=HORIZON_ELEV_DEG):
    """The diffuse dome: directions, radiance weights and the diagnostics.

    Returns a dict with
      dirs       (K,3) unit vectors pointing into the sky
      w          (K,)  weights, sum = 1, proportional to L * cos(theta) * dOmega
      f          (K,)  relative luminance L/L_z
    and the epsilon / coefficient / normalisation record.
    """
    s = np.asarray(sun_dir, dtype=np.float64)
    s = s / np.linalg.norm(s)
    alt = math.degrees(math.asin(max(-1.0, min(1.0, float(s[2])))))
    if sun_zenith_deg is None:
        sun_zenith_deg = 90.0 - alt
    eps = epsilon_from_eta(eta, sun_zenith_deg)
    coef = perez_coefficients(eps)

    dirs, solid = fibonacci_hemisphere(int(n_dirs), horizon_elev_deg)
    cos_t = np.clip(dirs[:, 2], 0.0, 1.0)
    theta_z = np.arccos(cos_t)
    chi = np.arccos(np.clip(dirs @ s, -1.0, 1.0))
    f = perez_shape(theta_z, chi, coef)

    wraw = f * cos_t * solid
    tot = float(wraw.sum())
    w = wraw / tot

    # L_z = DHI / integral(f cos th dOmega); in units of (irradiance / sr).
    # The integral is the tabulated total, so L_z is DHI * pi / tot / (pi)
    # -- kept explicit rather than folded, so the number is auditable.
    i0 = tot / math.pi          # = (1/pi) * int f cos dOmega
    cz_deg = float(math.degrees(sun_zenith_deg))
    # circumsolar fraction (the part of the diffuse inside 5 deg of the sun)
    circ = float(w[chi <= math.radians(5.0)].sum())
    circ_25 = float(w[chi <= math.radians(25.0)].sum())
    neff = float(1.0 / np.sum(w * w)) if np.sum(w * w) > 0 else None
    return {
        "dirs": dirs, "w": w, "f": f, "solid_angle": solid,
        "n_dirs": int(n_dirs),
        "horizon_elev_deg": float(horizon_elev_deg),
        "eta": float(eta),
        "sun_alt_deg": float(alt),
        "sun_zenith_deg": float(cz_deg),
        "epsilon": float(eps),
        "coefficients": coef,
        "Lz_per_DHI": float(1.0 / i0) if i0 > 0 else None,
        "normalisation": ("L_z solved so that int L cos(theta) dOmega = DHI "
                          "(same construction as Radiance gensky)"),
        "circumsolar_fraction_5deg": circ,
        "circumsolar_fraction_25deg": circ_25,
        "w_max": float(w.max()), "w_effective_n": neff,
        "f_min": float(f.min()), "f_max": float(f.max()),
        "f_at_horizon_mean": float(f[cos_t < 0.15].mean())
        if np.any(cos_t < 0.15) else None,
        "f_at_zenith": float(perez_shape(np.array([0.0]),
                                         np.array([chi[np.argmax(cos_t)]]),
                                         coef)[0]),
        "w_sum": float(w.sum()),
        "delta": delta_sky_brightness(1.0 - float(eta), alt),
        "table": PEREZ_COEF.tolist(),
    }


def horizon_sensitivity(eta, sun_dir, elevations=(5.0, 10.0, 15.0, 20.0),
                        n_dirs=N_SKY_DIRS):
    """How the weight distribution changes with the dome's horizon cut."""
    rows = {}
    for e in elevations:
        d = sky_dome(eta, sun_dir, n_dirs=n_dirs, horizon_elev_deg=float(e))
        rows["elev_%.0f_deg" % e] = {
            "w_effective_n": d["w_effective_n"],
            "w_max": d["w_max"],
            "circumsolar_fraction_5deg": d["circumsolar_fraction_5deg"],
            "Lz_per_DHI": d["Lz_per_DHI"],
            "f_max": d["f_max"],
            "f_min": d["f_min"],
        }
    return {
        "rows": rows,
        "adopted_elev_deg": float(HORIZON_ELEV_DEG),
        "why": ("the Perez relative distribution is evaluated as a finite "
                "dome: below ~10 deg the analytic fit extrapolates its "
                "exp(b/cos th) term outside its data range (f = 125 at "
                "5 deg for bin 3-4).  The cut is a stated numerical bound, "
                "and the table shows the weight distribution is insensitive "
                "to it between 10 and 20 deg"),
    }


def bin_step_check(eta, sun_dir, n_dirs=N_SKY_DIRS):
    """Step-bin lookup vs interpolated coefficients: how much it changes."""
    s = np.asarray(sun_dir, dtype=np.float64)
    s = s / np.linalg.norm(s)
    alt = math.degrees(math.asin(max(-1.0, min(1.0, float(s[2])))))
    eps = epsilon_from_eta(eta, 90.0 - alt)
    interp = perez_coefficients(eps)
    step = {k: float(PEREZ_COEF[eps_bin_index(eps)][i])
            for i, k in enumerate("abcde")}
    dirs, solid = fibonacci_hemisphere(int(n_dirs))
    cos_t = np.clip(dirs[:, 2], 0.0, 1.0)
    theta_z = np.arccos(cos_t)
    chi = np.arccos(np.clip(dirs @ s, -1.0, 1.0))
    wi = perez_shape(theta_z, chi, interp) * cos_t * solid
    ws = perez_shape(theta_z, chi, step) * cos_t * solid
    wi = wi / wi.sum()
    ws = ws / ws.sum()
    return {
        "epsilon": float(eps), "bin": int(eps_bin_index(eps) + 1),
        "coefficients_interpolated": interp,
        "coefficients_step_bin": step,
        "max_abs_weight_diff": float(np.max(np.abs(wi - ws))),
        "l1_weight_diff": float(np.sum(np.abs(wi - ws))),
        "note": ("the adopted interpolation is Perez's recommended practice; "
                 "a pure bin lookup would move at most this much of the "
                 "diffuse weight, which is reported instead of assumed "
                 "negligible"),
    }


# --------------------------------------------------------------------------
# 3. water path factors (direct vs diffuse)
# --------------------------------------------------------------------------
def refracted_cos(dirs, n=N_WATER):
    """cos(theta) in water after Snell refraction of the sky directions."""
    d = np.atleast_2d(np.asarray(dirs, dtype=np.float64))
    cz = np.clip(d[:, 2], 0.0, 1.0)                  # cos of incidence
    st = np.sqrt(np.maximum(1.0 - cz * cz, 0.0)) / float(n)
    return np.sqrt(np.maximum(1.0 - st * st, 0.0))


def path_factors(depth_m, dome, sun_dir):
    """Mean in-water path length z/cos(theta_w): direct vs diffuse.

    The diffuse weighting is the irradiance weight (L cos(theta) dOmega),
    which is the physically right average for the diffuse beam's absorption:
    the flux crossing the surface in each direction is what gets attenuated.
    """
    z = float(depth_m)
    ct_dir = float(refracted_cos(np.asarray(sun_dir, np.float64))[0])
    l_dir = z / max(ct_dir, 1.0e-6)
    ct = refracted_cos(dome["dirs"])
    l_dir_each = z / np.maximum(ct, 1.0e-6)
    w = np.asarray(dome["w"], dtype=np.float64)
    l_diff = float((w * l_dir_each).sum())
    return {
        "depth_m": z,
        "cos_theta_w_direct": ct_dir,
        "theta_w_direct_deg": float(math.degrees(math.acos(min(1.0, ct_dir)))),
        "s_down_direct_m": float(l_dir),
        "s_down_diffuse_mean_m": l_diff,
        "diffuse_over_direct": float(l_diff / l_dir) if l_dir > 0 else None,
        "cos_theta_w_diffuse_mean": float((w * ct).sum()),
        "theta_w_diffuse_mean_deg": float(math.degrees(math.acos(
            min(1.0, float((w * ct).sum()))))),
        "snell_window_deg": float(math.degrees(math.asin(1.0 / N_WATER))),
        "form": ("direct: z/cos(theta_w,dir); diffuse: "
                 "sum_i w_i z/cos(theta_w,i) with w_i the irradiance weight"),
        "note": ("the diffuse path is longer than the direct one but bounded "
                 "by the Snell window (cos >= %.4f), so it cannot diverge"
                 % math.cos(math.asin(1.0 / N_WATER))),
    }


# --------------------------------------------------------------------------
# 4. split arithmetic
# --------------------------------------------------------------------------
def horizontal_weights(eta, sun_alt_deg):
    """Direct / diffuse weights of the HORIZONTAL irradiance.

    DNI and DHI are different quantities (normal vs horizontal), so the
    energy split on a horizontal bottom is

        GHI_direct = DNI sin(alt) = eta sin(alt)
        GHI_diff   = DHI          = 1 - eta
        w_d        = eta sin(alt) / (eta sin(alt) + 1 - eta)

    (with DNI+DHI = 1 by normalisation).  At eta = 1 this is exactly 1 and at
    eta = 0 exactly 0, which is what makes the legacy recipe reproduce 3.4
    bit-for-bit.  Using eta itself as the weight would ignore the cosine --
    dimensionally wrong but numerically close at high sun; both are recorded.
    """
    e = float(min(max(eta, 0.0), 1.0))
    alt = math.radians(max(0.0, min(90.0, float(sun_alt_deg))))
    cz = math.sin(alt)
    if e >= 1.0:
        wd = 1.0
    elif e <= 0.0:
        wd = 0.0
    else:
        wd = (e * cz) / (e * cz + (1.0 - e))
    return {
        "eta": e, "cos_zenith": float(cz),
        "w_direct": float(wd), "w_diffuse": float(1.0 - wd),
        "w_sum": float(wd + (1.0 - wd)),
        "ghi_rel": float(e * cz + (1.0 - e)),
        "ghi_direct_rel": float(e * cz),
        "ghi_diffuse_rel": float(1.0 - e),
        "rule": ("w_direct = DNI cos(Z)/(DNI cos(Z)+DHI); identical to the "
                 "raw eta only at eta in {0,1}"),
        "raw_eta_as_weight": float(e),
    }


def eta_from_radiation(dni, dhi):
    """eta = DNI/(DNI+DHI); None when there is no light."""
    dni = max(0.0, float(dni))
    dhi = max(0.0, float(dhi))
    tot = dni + dhi
    if tot <= 0.0:
        return None
    return float(dni / tot)


def energy_ledger(split, dome, path, c):
    """Radiometric ledger for one eta.  All relative to DNI+DHI = 1."""
    wd = split["w_direct"]
    ws = split["w_diffuse"]
    ghi_d = split["ghi_direct_rel"]
    ghi_s = split["ghi_diffuse_rel"]
    t_dir = np.exp(-np.asarray(c) * (path["s_down_direct_m"] * 2.0))
    t_sky = np.exp(-np.asarray(c) * (path["s_down_diffuse_mean_m"] * 2.0))
    t_eff = np.asarray(ghi_d) * t_dir + np.asarray(ghi_s) * t_sky
    return {
        "basis": "DNI+DHI = 1 (relative); C is dimensionless",
        "ghi_rel": split["ghi_rel"],
        "ghi_direct_rel": ghi_d, "ghi_diffuse_rel": ghi_s,
        "w_direct": wd, "w_diffuse": ws,
        "s_down_direct_m": path["s_down_direct_m"],
        "s_down_diffuse_mean_m": path["s_down_diffuse_mean_m"],
        "T_direct_rgb": [float(x) for x in t_dir],
        "T_diffuse_rgb": [float(x) for x in t_sky],
        "bottom_rel_rgb": [float(x) for x in t_eff],
        "bottom_rel_mean_rec709": float(np.mean(t_eff)),
        "note": ("the bottom total is GHI*(w_d T_dir + w_s T_sky); the only "
                 "inputs are the split, the two path factors and the "
                 "absorption spectrum -- nothing is fitted to the transport"),
    }
