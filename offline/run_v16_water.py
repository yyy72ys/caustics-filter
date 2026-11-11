"""run_v16_water -- chunk 12 (BRIEF-V16-SCATTER): the Fournier-Forand phase
function and a halo kernel built from its angular CDF.

What this file is
-----------------
v12_water.py (version 3.6) models single scattering with a Henyey-Greenstein
phase function at g = 0.85 and represents the halo as a Gaussian blur whose
sigma is the second moment <theta^2> of that phase function.  The 2026-11-12
measurement research refutes HG-0.85 for the algae/detritus/clear presets:
Petzold's average-particle phase function has <cos> = 0.924, B = 0.0183 and
<theta^2> = 0.177 rad^2, while HG-0.85 has B = 0.0361 and <theta^2> = 0.368.
No value of g fixes HG's SHAPE (a moment-matched Gaussian puts 4.2% of the
scattered light inside 5 deg where the measurement puts 55.9%).

This module adds, WITHOUT editing any v5-v15 file:

  1. the Fournier-Forand (Fournier & Forand 1994; Fournier & Jonasz 1999)
     phase function in the Ocean Optics Web Book form (two terms, Eq. 1 of
     https://oceanopticsbook.info/view/scattering/the-fournier-forand-phase-function),
     normalised by quadrature on a log-spaced theta grid [1e-6, pi];
  2. the halo kernel built from that phase function's angular CDF used
     directly in the v12 displacement integral (replacing the <theta^2>
     moment-matched Gaussian), selectable beside `kernel: gaussian`;
  3. the scorecard / moment machinery of the brief's Stage A.

Energy bookkeeping is DELIBERATELY unchanged from v12
------------------------------------------------------
The brief pre-registers that H (the halo energy fraction) is phase-function
independent in this model class and must not move.  This module therefore
takes the weights T_bal / H / A from ``v12_water.partition`` VERBATIM (same
code, same floats, same argument g = the preset's v12 g), and the phase
function acts only on the ANGULAR SHAPE of the halo (the kernel).  For the
same a, b, L and g this makes H bit-identical to the v12 run by construction,
which is exactly what Stage B tests: any movement would be a code leak.

The <theta^2> moment of the selected phase function is still computed and
reported (``sigma_diag``); it is a DIAGNOSTIC in this family, not the kernel
width, because the FF-CDF kernel is not a Gaussian.

Kernel construction (the displacement integral, stated once)
------------------------------------------------------------
v12:  sigma^2 = (<theta^2>/2) (b/c^3) D2(cL),  D2(x) = x^2-2x+2-2exp(-x).
Here the displacement is sampled instead of moment-matched.  A photon that
scatters once at path position s (density b exp(-c s) on [0, L]) travels the
remaining path at angle theta and lands at

    r = (L - s) tan(theta)          [the v12 convention, in the image plane]

with theta drawn from the selected phase function.  The radial mass profile
of r is the kernel; it is radially symmetric, so it is applied as a periodic
FFT convolution exactly like ``transport.gaussian_blur_fft`` (same convention,
so old vs new differs only in the kernel SHAPE).  The kernel is truncated at
r_max px and renormalised to unit mass over the captured region; the captured
fraction and the truncation loss are printed for every state.
"""
import math

import numpy as np

import v12_water as W
from transport import gaussian_blur_fft

# --------------------------------------------------------------------------
# phase-function parameters (the published fits -- NOT knobs)
# --------------------------------------------------------------------------
#: theta grid of the quadrature.  Log-spaced, 1e-6 rad to pi: the FF function
#: varies over ~6 decades near the forward peak.
THETA_MIN_RAD = 1e-6
N_THETA = 60001

#: per-preset Fournier-Forand parameters.  ``algae_pond`` is the default
#: medium of v12 (DEFAULT_MEDIUM_PRESET); ``clear`` is the Jonasz clear-water
#: fit; ``sediment`` has NO published FF fit for a measured sediment-dominated
#: freshwater VSF, so its parameters are an ASSUMPTION and every sediment
#: state carries phase_assumption = True plus the envelope.
FF_PHASE = {
    "algae_pond": {
        "n": 1.10, "mu": 3.5835,
        "label": "Petzold average-particle fit",
        "provenance": ("n = 1.10, mu = 3.5835 is the Fournier-Forand pair "
                       "that best fits Petzold's average-particle phase "
                       "function (Ocean Optics Web Book, Fig. 5, B = 0.0183 "
                       "contour); a MEASUREMENT fit, not a knob"),
        "phase_assumption": False,
    },
    "clear": {
        "n": 1.098, "mu": 3.37,
        "label": "Jonasz clear-water fit",
        "provenance": ("n = 1.098, mu = 3.37 is the FF fit quoted for clear "
                       "ocean water (Jonasz family); a published fit, not a "
                       "knob"),
        "phase_assumption": False,
    },
    "sediment": {
        "n": 1.15, "mu": 4.0,
        "label": "nominal mineral sediment (ASSUMPTION)",
        "provenance": ("NO published FF fit for a measured sediment-"
                       "dominated freshwater VSF was found; n = 1.15, "
                       "mu = 4.0 is a NOMINAL choice placed between the "
                       "biogenic and mineral FF families.  This is an "
                       "ASSUMPTION, not a measurement, and every state that "
                       "uses it prints phase_assumption = True and the "
                       "envelope below"),
        "phase_assumption": True,
    },
}

#: the sediment envelope the brief requires: [Petzold fit, nominal sediment].
SEDIMENT_ENVELOPE = (
    {"preset": "algae_pond", "n": 1.10, "mu": 3.5835,
     "label": "Petzold average-particle fit (measurement)"},
    {"preset": "sediment", "n": 1.15, "mu": 4.0,
     "label": "nominal mineral sediment (ASSUMPTION)"},
)

#: the two lanes this family exposes.
PHASES = ("ff", "hg")
KERNELS = ("ff-cdf", "gaussian")
PHASE_DEFAULT = "ff"          # the candidate default under test
KERNEL_DEFAULT = "ff-cdf"     # the candidate default under test
HG_G_DEFAULT = W.G_DEFAULT    # 0.85, the v12 legacy lane

#: Stage A's declared beta(90 deg) reference (the brief's pre-registration):
#: Petzold's flat-back level, NOT the FF value.
BETA90_REFERENCE_SR = 0.00291
BETA90_TOL_FACTOR = 1.5


# --------------------------------------------------------------------------
# the Fournier-Forand phase function
# --------------------------------------------------------------------------
def ff_phase(theta, n, mu):
    """Fournier & Jonasz (1999) FF phase function in 1/sr (OOBN Eq. 1).

    beta(theta) = 1/(4 pi (1-d)^2 d^nu) [ nu(1-d) - (1-d^nu)
                   + (d(1-d^nu) - nu(1-d)) / sin^2(theta/2) ]
                 + (1 - d180^nu) / (16 pi (d180 - 1) d180^nu) (3 cos^2 - 1)

    d = 4/(3 (n-1)^2) sin^2(theta/2),  nu = (3 - mu)/2,  d180 = d(pi).
    The SECOND term is part of the Fournier & Jonasz form; omitting it moves
    beta(90 deg) from 0.004197 to 0.004670 sr^-1.  Values are NOT normalised
    here -- see ``phase_table``, which divides by the quadrature norm.
    """
    th = np.asarray(theta, dtype=np.float64)
    s2 = np.sin(0.5 * th) ** 2
    k = 4.0 / (3.0 * (n - 1.0) ** 2)
    delta = k * s2
    nu = (3.0 - float(mu)) / 2.0
    dnu = np.exp(nu * np.log(np.maximum(delta, 1e-300)))
    one_minus = 1.0 - delta
    pref = 1.0 / (4.0 * math.pi * one_minus ** 2 * dnu)
    bracket = (nu * one_minus - (1.0 - dnu)
               + (delta * (1.0 - dnu) - nu * one_minus)
               / np.maximum(s2, 1e-300))
    d180 = k
    d180nu = d180 ** nu
    second = ((1.0 - d180nu) / (16.0 * math.pi * (d180 - 1.0) * d180nu)
              * (3.0 * np.cos(th) ** 2 - 1.0))
    return pref * bracket + second


def hg_phase_theta(theta, g):
    """HG on the theta grid, through v12's own function (same code path)."""
    return W.hg_phase(np.cos(np.asarray(theta, dtype=np.float64)), float(g))


def phase_table(phase="ff", preset=None, g=None, n_theta=N_THETA):
    """theta grid + analytic phase + quadrature normalisation.

    The normalisation constant is the integral of the analytic form over the
    sphere; it is 1.0 for HG by construction and 0.988-1.000 for the FF fits
    (the analytic FF is not exactly normalised).  Every printed beta(theta)
    divides by it.
    """
    theta = np.geomspace(THETA_MIN_RAD, math.pi, int(n_theta))
    if phase == "hg":
        gg = HG_G_DEFAULT if g is None else float(g)
        ph = hg_phase_theta(theta, gg)
        out = {"phase": "hg", "g": gg, "n": None, "mu": None,
               "label": "Henyey-Greenstein g=%.4g (v12 legacy lane)" % gg,
               "provenance": "v12 legacy lane (v12_water.hg_phase, unchanged)"}
    elif phase == "ff":
        P = FF_PHASE[preset]
        ph = ff_phase(theta, P["n"], P["mu"])
        out = {"phase": "ff", "g": None, "n": P["n"], "mu": P["mu"],
               "label": P["label"], "provenance": P["provenance"],
               "preset": preset}
    else:
        raise ValueError("unknown phase %r (expected one of %s)"
                         % (phase, PHASES))
    domega = 2.0 * math.pi * np.sin(theta)
    norm = float(np.trapz(ph * domega, theta))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("phase function has a bad norm (%r)" % (norm,))
    out.update({"theta_rad": theta, "phase_1_per_sr_raw": ph,
                "normalisation": norm, "domega": domega,
                "phase_1_per_sr": ph / norm})
    return out


def phase_moments(table):
    """<cos>, B, <theta^2>, F(<1,5,10 deg), beta(1,5,10,90 deg), F(>30 deg).

    All by quadrature on the table's log grid; the integrals are over the
    sphere with dOmega = 2 pi sin(theta) dtheta.  ``F(<x)`` is the fraction of
    the scattered light inside x degrees; B is the fraction beyond 90 deg.
    """
    th = table["theta_rad"]
    p = table["phase_1_per_sr"]
    w = table["domega"]
    norm = float(np.trapz(p * w, th))

    def _frac(mask):
        return float(np.trapz(np.where(mask, p * w, 0.0), th) / norm)

    def _beta(deg):
        v = np.interp(math.radians(deg), th, p)
        return float(v)

    out = {
        "label": table["label"], "phase": table["phase"],
        "n": table["n"], "mu": table["mu"], "g": table["g"],
        "normalisation_quadrature": table["normalisation"],
        "mean_cos": float(np.trapz(p * w * np.cos(th), th) / norm),
        "theta2_rad2": float(np.trapz(p * w * th * th, th) / norm),
        "backscatter_fraction_B": _frac(th >= 0.5 * math.pi),
        "F_lt_1deg": _frac(th <= math.radians(1.0)),
        "F_lt_5deg": _frac(th <= math.radians(5.0)),
        "F_lt_10deg": _frac(th <= math.radians(10.0)),
        "F_gt_30deg": _frac(th > math.radians(30.0)),
        "beta_1deg_1_per_sr": _beta(1.0),
        "beta_5deg_1_per_sr": _beta(5.0),
        "beta_10deg_1_per_sr": _beta(10.0),
        "beta_90deg_1_per_sr": _beta(90.0),
    }
    out["beta90_over_reference"] = (out["beta_90deg_1_per_sr"]
                                    / BETA90_REFERENCE_SR)
    out["provenance"] = table.get("provenance")
    return out


def scorecard_clauses(row, reference=BETA90_REFERENCE_SR,
                      tol_beta=BETA90_TOL_FACTOR):
    """The brief's five Stage-A clauses for the algae/generic preset."""
    lo_b = reference / tol_beta
    hi_b = reference * tol_beta
    t2 = 0.177
    clauses = {
        "mean_cos_0.924_pm_0.02": {
            "target": [0.904, 0.944], "value": row["mean_cos"],
            "pass": bool(abs(row["mean_cos"] - 0.924) <= 0.02)},
        "B_0.0183_pm_0.005": {
            "target": [0.0133, 0.0233],
            "value": row["backscatter_fraction_B"],
            "pass": bool(abs(row["backscatter_fraction_B"] - 0.0183)
                         <= 0.005)},
        "theta2_0.177_pm_15pct": {
            "target": [t2 * 0.85, t2 * 1.15], "value": row["theta2_rad2"],
            "pass": bool(abs(row["theta2_rad2"] - t2) <= 0.15 * t2)},
        "F_lt_5deg_0.559_pm_0.10": {
            "target": [0.459, 0.659], "value": row["F_lt_5deg"],
            "pass": bool(abs(row["F_lt_5deg"] - 0.559) <= 0.10)},
        "beta90_within_x1.5_of_%.5f" % reference: {
            "target": [lo_b, hi_b], "value": row["beta_90deg_1_per_sr"],
            "pass": bool(lo_b <= row["beta_90deg_1_per_sr"] <= hi_b),
            "ratio_to_reference": (row["beta_90deg_1_per_sr"] / reference)},
    }
    n_pass = sum(1 for c in clauses.values() if c["pass"])
    return {"clauses": clauses, "n_clauses": len(clauses),
            "n_pass": int(n_pass), "pass_all": bool(n_pass == len(clauses))}


def phase_provenance(preset=None, phase="ff"):
    """The provenance block that must travel with every v16 water state."""
    if phase == "hg":
        return {"phase": "hg", "g": HG_G_DEFAULT,
                "provenance": ("v12 legacy lane: HG g=%.2f, the v3.6 default; "
                               "REFUTED for algae/detritus/clear by the "
                               "2026-11-12 research, retained for regression"
                               % HG_G_DEFAULT),
                "phase_assumption": False, "envelope": None}
    P = FF_PHASE[preset]
    out = {"phase": "ff", "preset": preset, "n": P["n"], "mu": P["mu"],
           "label": P["label"], "provenance": P["provenance"],
           "phase_assumption": bool(P["phase_assumption"])}
    if P["phase_assumption"]:
        out["envelope"] = [dict(e) for e in SEDIMENT_ENVELOPE]
        out["envelope_note"] = ("sediment's phase function is an ASSUMPTION: "
                                "the state is bracketed by the two FF rows "
                                "above, and no measurement is claimed")
    else:
        out["envelope"] = None
    return out


# --------------------------------------------------------------------------
# the halo kernel: FF angular CDF in the displacement integral
# --------------------------------------------------------------------------
def _phase_bin_masses(table, n_bins=2400):
    """Mass per log-spaced theta bin (sums to 1) from the fine quadrature."""
    th = table["theta_rad"]
    p = table["phase_1_per_sr"]
    w = table["domega"]
    norm = float(np.trapz(p * w, th))
    cum_fine = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] * w[1:] + p[:-1] * w[:-1])
                                               * np.diff(th))]) / norm
    edges = np.geomspace(th[0], th[-1], int(n_bins) + 1)
    cum_edges = np.interp(edges, th, cum_fine)
    masses = np.diff(cum_edges)
    masses = np.clip(masses, 0.0, None)
    return edges, masses / masses.sum()


def _scatter_positions(c, L, n_s=192):
    """Gauss-Legendre nodes on [0, L] with density proportional to exp(-c s)."""
    x, wq = np.polynomial.legendre.leggauss(int(n_s))
    s = 0.5 * float(L) * (x + 1.0)
    w = 0.5 * float(L) * wq * np.exp(-float(c) * s)
    return s, w / w.sum()


def radial_profile(table, b, c, L, dx_m, r_max_px, n_bins=1500, n_s=192,
                   n_theta_bins=2400):
    """Radial mass profile of the single-scatter displacement, in pixels.

    r = (L - s) |tan(theta)| / dx, theta ~ phase, s ~ truncated exponential.
    Returns the radial bin edges, the mass per bin (summing to ``captured``),
    the cumulative fraction of the FULL halo mass, and diagnostics.
    """
    dx = float(dx_m)
    if dx <= 0.0:
        raise ValueError("dx_m must be > 0")
    edges_th, mass_th = _phase_bin_masses(table, n_theta_bins)
    th_rep = np.sqrt(edges_th[:-1] * edges_th[1:])          # log mid-point
    s, ws = _scatter_positions(c, L, n_s=n_s)
    tan_th = np.abs(np.tan(th_rep))
    r_px = np.outer(L - s, tan_th) / dx                     # (n_s, n_theta)
    wgt = np.outer(ws, mass_th)
    edges = np.geomspace(1e-3, float(r_max_px), int(n_bins) + 1)
    mass, _ = np.histogram(r_px.ravel(), bins=edges,
                           weights=wgt.ravel())
    captured = float(mass.sum())
    cum = np.cumsum(mass)
    centres = np.sqrt(edges[:-1] * edges[1:])
    area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        density = np.where(area > 0, mass / area, 0.0)
    p = {
        "edges_px": edges, "centres_px": centres, "mass": mass,
        "density_per_px2": density, "cumulative": cum,
        "captured_fraction": captured,
        "r_max_px": float(r_max_px),
        "rms_radius_px": float(math.sqrt(
            max(float((mass * centres ** 2).sum() / max(captured, 1e-300)),
                0.0))),
        "n_s_nodes": int(n_s), "n_theta_bins": int(n_theta_bins),
        "b": float(b), "c": float(c), "L_m": float(L), "dx_m": dx,
        "truncation_loss_fraction": float(1.0 - captured),
    }
    return p


def encircled_energy(profile, radii_px):
    """Fraction of the FULL halo mass inside r (linear interp of the CDF)."""
    edges = profile["edges_px"]
    cdf = np.concatenate([[0.0], profile["cumulative"]])
    return [float(np.interp(float(r), edges, cdf)) for r in radii_px]


def kernel_2d(profile, n):
    """The periodic convolution kernel on an n x n grid (offsets wrapped).

    Value = radial density per px^2 at the pixel's radius, normalised so the
    discrete kernel sums to 1 (the same unit-mass convention
    ``transport.gaussian_blur_fft`` obeys).  The normalisation factor is
    returned so the report can show how far the discrete grid is from the
    analytic annulus mass.
    """
    n = int(n)
    dy = np.arange(n)
    dy = np.where(dy < n // 2, dy, dy - n)
    dx = np.arange(n)
    dx = np.where(dx < n // 2, dx, dx - n)
    R = np.sqrt(dy[:, None].astype(np.float64) ** 2
                + dx[None, :].astype(np.float64) ** 2)
    K = np.interp(R.ravel(), profile["centres_px"],
                  profile["density_per_px2"], left=0.0, right=0.0)
    K = K.reshape((n, n))
    raw = float(K.sum())
    if raw <= 0.0:
        raise ValueError("empty kernel")
    K = K / raw
    return K, {"raw_sum_before_normalisation": raw,
               "n": n, "sum_after": float(K.sum())}


def ff_blur_fft(C, Kp):
    """Periodic FFT convolution with a unit-mass kernel (gaussian_blur_fft
    convention: mean preserved, wrap-around)."""
    C = np.asarray(C)
    if C.ndim == 3:
        out = np.empty_like(C, dtype=np.float64)
        for ch in range(C.shape[2]):
            out[:, :, ch] = ff_blur_fft(C[:, :, ch], Kp)
        return out
    if Kp.shape != C.shape:
        raise ValueError("kernel shape %s != field shape %s"
                         % (Kp.shape, C.shape))
    F = np.fft.rfft2(np.asarray(C, dtype=np.float64))
    Kf = np.fft.rfft2(Kp)
    return np.fft.irfft2(F * Kf, s=C.shape)


def gaussian_kernel_profile(sigma_px, radii_px):
    """Encircled energy of a 2-D Gaussian blur kernel: 1-exp(-r^2/2 sigma^2)."""
    return [float(1.0 - math.exp(-0.5 * (float(r) / float(sigma_px)) ** 2))
            for r in radii_px]


# --------------------------------------------------------------------------
# optics: v12 weights (verbatim) + the selected phase/kernel
# --------------------------------------------------------------------------
def ff_params_for_preset(preset):
    key = preset if preset in FF_PHASE else W.DEFAULT_MEDIUM_PRESET
    return key, FF_PHASE[key]


def sigma_diagnostic(b_green, c_green, L, theta2):
    """sigma^2 = (<theta^2>/2) (b/c^3) D2(cL) -- the v12 formula, any moment."""
    b = max(float(b_green), 0.0)
    c = float(c_green)
    if b <= 0.0 or c <= 0.0:
        return {"sigma_m": 0.0, "sigma2_m2": 0.0, "theta2_rad2": float(theta2),
                "D2": 0.0, "x": 0.0}
    x = c * float(L)
    D2 = x * x - 2.0 * x + 2.0 - 2.0 * math.exp(-x)
    sig2 = 0.5 * float(theta2) * (b / (c ** 3)) * D2
    return {"sigma_m": float(math.sqrt(max(sig2, 0.0))),
            "sigma2_m2": float(sig2), "theta2_rad2": float(theta2),
            "D2": float(D2), "x": float(x), "b": b, "c": c,
            "L_m": float(L)}


def kernel_r_max_px(sigma_px, grid, cap_frac=0.45, n_sigma=8.0,
                    floor=48.0):
    grid = int(grid)
    r = max(float(floor), float(n_sigma) * float(sigma_px))
    return float(min(r, cap_frac * grid))


def optics_v16(spec, env, depth_m, grid, patch_m, phase=PHASE_DEFAULT,
               kernel=KERNEL_DEFAULT, preset=None, r_max_px=None,
               view_zenith_deg=0.0):
    """Water optics with the v12 energy budget and the selected halo shape.

    The weights (T_bal, H, A) come from ``v12_water.ab_optics`` unchanged, so
    H is bit-identical to the v12 run for the same a, b, L and g.  The phase
    function changes (i) the moment diagnostic and (ii) the kernel.
    """
    base = W.ab_optics(spec, env, depth_m, grid, patch_m,
                       view_zenith_deg=view_zenith_deg)
    key = spec.get("preset") or preset
    pk, P = ff_params_for_preset(key)
    if phase == "ff":
        table = phase_table("ff", preset=pk)
        mom = phase_moments(table)
        prov = phase_provenance(pk, "ff")
    elif phase == "hg":
        table = phase_table("hg", g=float(spec["g"]))
        mom = phase_moments(table)
        prov = phase_provenance(pk, "hg")
    else:
        raise ValueError("unknown phase %r" % (phase,))

    b_g = float(np.asarray(spec["b"], np.float64)[1])
    c_g = float(base["partition"]["c"][1])
    L = float(base["partition"]["L_m"])
    dx = float(base["dx_m"])
    if phase == "hg":
        # the legacy lane takes v12's own sigma object verbatim (same code
        # path as v12_water.sigma_from_scattering), so it is bit-identical
        sig = dict(base["sigma"])
        sig_px = float(base["sigma_px"])
    else:
        sig = sigma_diagnostic(b_g, c_g, L, mom["theta2_rad2"])
        sig_px = sig["sigma_m"] / dx if dx > 0 else 0.0

    out = dict(base)
    out["kernel"] = kernel
    out["phase_kind"] = phase
    out["phase"] = {"label": mom["label"], "n": mom["n"], "mu": mom["mu"],
                    "g": mom["g"],
                    # the tables are big; keep only scalar moments on the
                    # optics object and let the driver recompute tables.
                    "mean_cos": mom["mean_cos"],
                    "theta2_rad2": mom["theta2_rad2"],
                    "backscatter_fraction_B": mom["backscatter_fraction_B"],
                    "beta_90deg_1_per_sr": mom["beta_90deg_1_per_sr"]}
    out["phase_provenance"] = prov
    out["sigma_v12_hg"] = base["sigma"]
    out["sigma_px_v12_hg"] = float(base["sigma_px"])
    out["sigma"] = sig
    out["sigma_px"] = float(sig_px)
    out["sigma_diag_m"] = sig["sigma_m"]
    out["sigma_diag_px"] = float(sig_px)
    out["sigma_diag_source"] = ("v12_water.sigma_from_scattering (verbatim)"
                               if phase == "hg" else
                               "sigma_diagnostic with the FF quadrature <theta^2>")

    if kernel == "ff-cdf":
        r_max = (float(r_max_px) if r_max_px is not None
                 else kernel_r_max_px(sig_px, grid))
        prof = radial_profile(table, b_g, c_g, L, dx, r_max)
        Kp, kinfo = kernel_2d(prof, int(grid))
        out["kernel_profile"] = prof
        out["kernel_2d"] = Kp
        out["kernel_info"] = kinfo
        out["kernel_encircled_px"] = encircled_energy(prof, (8.5, 32.0, 79.0))
    else:
        out["kernel_profile"] = None
        out["kernel_2d"] = None
    return out


# --------------------------------------------------------------------------
# the composite: v12 body with the kernel substitution
# --------------------------------------------------------------------------
def render_frame_v16(fields, C, photo_lin, env, Wd, Hd, optics,
                     ambient_scale=None, c_layer_exposure=None,
                     substrate_gain_rgb=(1.0, 1.0, 1.0), g1=1.0,
                     lowpass_sigma_px=None, glitter_gain=0.0,
                     return_debug=False):
    """One v16 still.

    (a) phase == hg AND kernel == gaussian  ->  DELEGATES to
        ``v12_water.render_frame_v12`` with the v12 optics, so the legacy lane
        cannot move by construction (Stage D item i).
    (b) otherwise this is the v12 a/b composite body with ONE substitution:
        the halo field is the FF-CDF convolution instead of the Gaussian blur.
        Stage D asserts that forcing the kernel to gaussian here reproduces
        v12 byte for byte, so the copied body is proven faithful.
    """
    from v5_composite import A_V3, B_SIGMA_V3, E_V3
    A = A_V3 if ambient_scale is None else float(ambient_scale)
    E = E_V3 if c_layer_exposure is None else float(c_layer_exposure)
    lp = B_SIGMA_V3 if lowpass_sigma_px is None else float(lowpass_sigma_px)

    if optics.get("phase_kind") == "hg" and optics.get("kernel") == "gaussian":
        return W.render_frame_v12(
            fields, C, photo_lin, env, Wd, Hd, optics,
            ambient_scale=A, c_layer_exposure=E,
            substrate_gain_rgb=substrate_gain_rgb, g1=g1,
            lowpass_sigma_px=lp, glitter_gain=glitter_gain,
            return_debug=return_debug)

    from image import (IOR, SKY_REF, aces_narkowicz, e_base_dir_from_env,
                       linear_to_srgb)
    from state import sample_bicubic
    from v5_composite import (_sample_c, _state_gradients, _wobbled_alb,
                              v3_beer)

    S = float(fields["S"])
    N_buf = int(np.asarray(C).shape[0])
    depth = float(env["depthM"])
    diffuse = float(env["diffuse"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    Ebase = e_base_dir_from_env(env)
    gain = np.asarray(substrate_gain_rgb, dtype=np.float64)

    w_bal = np.asarray(optics["w_ballistic"], dtype=np.float64)
    w_halo = np.asarray(optics["w_halo"], dtype=np.float64)
    Cc = np.asarray(C, dtype=np.float64)
    if float(np.max(w_halo)) > 0.0:
        if optics["kernel"] == "ff-cdf":
            Ch = ff_blur_fft(Cc, optics["kernel_2d"]).astype(np.float64)
        else:
            sig_px = float(optics["sigma_px"])
            Ch = (gaussian_blur_fft(Cc, sig_px).astype(np.float64)
                  if sig_px > 0.0 else Cc)
    else:
        Ch = Cc

    gx_s, gy_s = _state_gradients(fields, lp)
    iy, ix = np.mgrid[0:Hd, 0:Wd]
    u = (ix.astype(np.float64) + 0.5) / Wd
    v = (iy.astype(np.float64) + 0.5) / Hd
    px = u * S
    py = v * S
    gx = sample_bicubic(gx_s.astype(np.float32), px, py, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nz = 1.0 / nl

    beer_v3 = v3_beer(env)
    scat_amt = float(1.0 - beer_v3.mean())
    amb_gain = (0.35 + 0.65 * diffuse) * (0.3 + 0.7 * sunI)
    E_amb = scatter * (scat_amt * amb_gain) * A

    col = np.zeros((Hd, Wd, 3), dtype=np.float64)
    alb_all = np.zeros((Hd, Wd, 3), dtype=np.float64)
    direct_all = np.zeros((Hd, Wd, 3), dtype=np.float64)
    irr_all = np.zeros((Hd, Wd, 3), dtype=np.float64)
    for bb in range(3):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        alb = _wobbled_alb(photo_lin[:, :, bb:bb + 1], u, v, shx, shy,
                           S)[:, :, 0]
        alb_all[:, :, bb] = alb
        c_bal = _sample_c(Cc, u, v, shx, shy, S, N_buf)
        c_hal = _sample_c(Ch, u, v, shx, shy, S, N_buf)
        irr = float(w_bal[bb]) * c_bal + float(w_halo[bb]) * c_hal
        irr_all[:, :, bb] = irr
        direct = float(Ebase[bb]) * irr * float(g1)
        direct_all[:, :, bb] = alb * float(gain[bb]) * direct \
            * strength * exposure
        bottom = alb * float(gain[bb]) * (E * direct) + alb * float(E_amb[bb])
        col[:, :, bb] = bottom * strength

    cos_t = nz
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - diffuse) \
        + SKY_REF[None, None, :] * diffuse
    sky_gain = (0.10 + 0.30 * (1.0 - diffuse)) \
        * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * (sky_gain * A)

    if glitter_gain and float(glitter_gain) > 0.0:
        from v5_physics import cox_munk_anisotropy
        an = cox_munk_anisotropy(float(env["windSpeedMS"]))
        s = np.asarray(env["sunDir"], dtype=np.float64)
        s = s / np.linalg.norm(s)
        Hv = s + np.array([0.0, 0.0, 1.0])
        Hv = Hv / np.linalg.norm(Hv)
        ndh = np.clip(-nx * Hv[0] - ny * Hv[1] + nz * Hv[2], 0.0, 1.0)
        ma = np.sqrt(np.maximum(1.0 - ndh * ndh, 0.0))
        mc = np.zeros_like(ma)
        p = (np.exp(-0.5 * (ma * ma / max(an["sigma_a2"], 1e-9)
                            + mc * mc / max(an["sigma_c2"], 1e-9)))
             / (2.0 * np.pi * math.sqrt(an["sigma_a2"] * an["sigma_c2"])))
        spark = (np.power(ndh, 1400.0) * 2.0 + np.power(ndh, 180.0) * 0.35)
        col += sunColor[None, None, :] * (spark * p
                                          * float(glitter_gain))[:, :, None]

    lin = col * exposure
    aces = aces_narkowicz(np.maximum(lin, 0.0))
    out = (np.clip(linear_to_srgb(aces), 0.0, 1.0) * 255.0
           + 0.5).astype(np.uint8)
    if return_debug:
        return out, {"lin": lin, "aces": aces, "alb": alb_all,
                     "direct_lin": direct_all,
                     "direct_irradiance": irr_all,
                     "w_ballistic": [float(x) for x in w_bal],
                     "w_halo": [float(x) for x in w_halo],
                     "beer_down": [float(x) for x in optics["beer_down"]],
                     "beer_up": [float(x) for x in optics["beer_up"]],
                     "g1": float(g1), "E": E, "A": A,
                     "kernel": optics.get("kernel"),
                     "blur_sigma_px": float(optics.get("sigma_px", 0.0)),
                     "substrate_gain_rgb": [float(x) for x in gain]}
    return out


def effective_c_field_v16(C, optics):
    """The C multiplier the transport applies (green-channel weights)."""
    C = np.asarray(C, dtype=np.float64)
    wb = float(optics["w_ballistic"][1])
    wh = float(optics["w_halo"][1])
    if wh <= 0.0:
        return wb * C
    if optics.get("kernel") == "ff-cdf":
        ch = ff_blur_fft(C, optics["kernel_2d"])
    else:
        sig = float(optics["sigma_px"])
        ch = gaussian_blur_fft(C, sig) if sig > 0.0 else C
    return wb * C + wh * ch
