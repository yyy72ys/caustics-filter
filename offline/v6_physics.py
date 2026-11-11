"""v6 physics: rain statistics -> ring-wave spectral band (PIPELINE-V6, chunk 2).

Every number here is a DERIVED quantity from a published formula in
report.md section 3 (rain on windless water) / section 5 (variable map).
The single exception, named explicitly, is the impact->ring-wave energy
conversion efficiency ``RING_ETA`` (PIPELINE-V6 section 2: "導出できない定数は
η ... の1つだけ").

Chain (PIPELINE-V6 section 1-2, report section 3.1-3.5):

  R [mm/h]
    -> Marshall-Palmer  N(D) = N0 exp(-Lambda D), Lambda = 4.1 R^-0.21
    -> Atlas-Ulbrich    v(D) = 9.65 - 10.3 exp(-0.6 D)          [m/s]
    -> F(R) = integral N(D) v(D) dD                             [drops m^-2 s^-1]
    -> flux-weighted mean impact kinetic energy Ebar            [J/drop]
    -> D0(R) = 3.67 / Lambda                                    [mm]
    -> k_ring(R) from D0 (bigger drops -> longer rings)         [rad/m]
    -> ring band variance  sigma_ring^2 = eta F Ebar tau / (rho g + sigma_s k^2)
    -> second effect: low-k (gravity band) damping factor

Nothing else in this module is a free parameter: N0, Lambda(R), v(D), F(R),
Ebar(R), D0(R), the 310-630 rad/m band, the gravity-capillary dispersion and
the ~150 ms impact-process lifetime all come from the report.

Ring placement (k_ring vs R)
---------------------------
report section 3.4 states the observed ring band (lambda 1-2 cm), and that the
ring peak frequency decreases with drop diameter while being nearly
independent of rain rate.  No closed-form frequency law is published in the
report, so the placement inside the observed band is an explicit MONOTONE
PLACEMENT RULE (not a calibration constant, not a new fit):

    lambda_ring(D0) = 10 mm + 10 mm * (D0 - D0(R=1)) / (D0(R=50) - D0(R=1))

i.e. the two band limits (1 cm, 2 cm) are anchored at the two ends of the
rain-rate range of the verification table, and the interpolation variable is
the mass-median drop diameter D0(R) (which is what report section 3.4 says
sets the ring scale).  Consequence at the reference rate: lambda(R=10) =
14.9 mm -> k = 422 rad/m, matching PIPELINE-V6's "k_ring ~ 420 rad/m
(lambda ~ 1.5 cm)" to 0.5%.  k_ring is monotone decreasing in R by
construction and stays inside 310-630 rad/m for R in [1, 50].
"""
import math

import numpy as np

# --------------------------------------------------------------------------
# physical constants (report sections 3.2, 3.4, 6.1)
# --------------------------------------------------------------------------
G = 9.81                 # m/s^2
RHO_W = 1000.0           # kg/m^3 (water, report section 3.4)
SIGMA_SURF = 0.072       # N/m   surface tension (report section 3.4)

# Marshall-Palmer (report section 3.1)
MP_N0 = 8000.0           # m^-3 mm^-1
MP_LAMBDA_A = 4.1        # mm^-1
MP_LAMBDA_P = -0.21

# Atlas-Ulbrich terminal velocity (report section 3.2), D in mm, v in m/s
AU_V0 = 9.65
AU_V1 = 10.3
AU_K = 0.6
AU_VALID_MM = (0.5, 5.0)
DSD_DMAX_MM = 6.0        # drops >~6 mm break up in fall (report section 3.1)

# Ring waves (report section 3.4)
RING_LAMBDA_BAND_MM = (10.0, 20.0)   # observed 1-2 cm band
RING_K_BAND = (310.0, 630.0)         # PIPELINE-V6 section 2 gate band [rad/m]
RING_K_DESIGN_REF = 420.0            # PIPELINE-V6 section 2 reference centre
RING_TAU_S = 0.15                    # impact process ~150 ms (report 3.4)
RING_RATE_ANCHORS = (1.0, 50.0)      # rain-rate ends of the verification table
RING_SIGMA_K_FRAC = 1.0 / 6.0        # bump width, as a fraction of the band

# THE one calibration constant of this chunk (PIPELINE-V6 section 2/6):
RING_ETA = 0.10
RING_ETA_BASIS = (
    "impact -> ring-wave energy conversion efficiency.  The report derives "
    "everything else in the rain chain (Lambda(R), v(D), F(R), D0(R), the "
    "k band, the lifetime) but not this.  Value 0.10: about a tenth of the "
    "impact kinetic energy survives into the ring band after the cavity, "
    "crown and Worthington-jet/splash phase -- which PIPELINE-V6 section 8 "
    "explicitly excludes from the model, so that energy is lost by "
    "construction.  Order-of-magnitude calibration, NOT derived; the "
    "sensitivity is recorded in numbers-v6.json (eta 0.03/0.10/0.30).")

# Verification table of report section 3.3 [drops m^-2 s^-1]
F_TABLE = {1.0: 1.3e3, 10.0: 4.2e3, 25.0: 6.4e3, 50.0: 8.6e3}
F_TABLE_TOL_FACTOR = 1.5   # PIPELINE-V6 section 6 gate 2


def dsd_lambda(rate_mm_h):
    """Marshall-Palmer slope Lambda = 4.1 R^-0.21  [1/mm]."""
    R = float(rate_mm_h)
    if R <= 0.0:
        raise ValueError("rain rate must be > 0 for DSD evaluation")
    return float(MP_LAMBDA_A * R ** MP_LAMBDA_P)


def dsd_number_density(rate_mm_h, D_mm):
    """N(D) = N0 exp(-Lambda D)  [m^-3 mm^-1]."""
    return MP_N0 * np.exp(-dsd_lambda(rate_mm_h) * np.asarray(D_mm, np.float64))


def fall_speed(D_mm):
    """Atlas-Ulbrich terminal velocity [m/s]; D in mm (valid 0.5-5 mm)."""
    return AU_V0 - AU_V1 * np.exp(-AU_K * np.asarray(D_mm, np.float64))


def mass_median_diameter(rate_mm_h):
    """D0 = 3.67 / Lambda  [mm] (report section 3.1)."""
    return float(3.67 / dsd_lambda(rate_mm_h))


def mean_diameter(rate_mm_h):
    """DSD mean diameter 1/Lambda [mm] (report section 3.1); readout only."""
    return float(1.0 / dsd_lambda(rate_mm_h))


def _dsd_grid(n=240001):
    return np.linspace(0.0, DSD_DMAX_MM, int(n))


def _trapz(y, x):
    """Trapezoidal integral, working on numpy 1.x (np.trapz) and 2.x."""
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(f(y, x))


def impact_flux(rate_mm_h, n=240001):
    """F = integral N(D) v(D) dD  [drops m^-2 s^-1] (report section 3.3).

    Three variants are returned because the Atlas-Ulbrich fit is only valid
    for D ~ 0.5-5 mm while Marshall-Palmer extends below 0.5 mm:

      F_clamped      : v = max(v, 0) over 0..6 mm  -- PHYSICAL variant used
                       for the ring energy budget (a negative terminal
                       velocity would mean upward rain flux).
      F_unclamped    : the raw fit over 0..6 mm; this is what reproduces the
                       report's own table (report: "integral taken over the
                       validity band and slightly beyond").
      F_validity_band: integration restricted to 0.5..5 mm -- the strict
                       validity band; much smaller because ~87% (R=1) of the
                       Marshall-Palmer drops sit below 0.5 mm.

    The gate (PIPELINE-V6 section 6 item 2) uses F_clamped against the table.
    """
    R = float(rate_mm_h)
    D = _dsd_grid(n)
    N = dsd_number_density(R, D)
    v_raw = fall_speed(D)
    v_cl = np.maximum(v_raw, 0.0)
    integ_cl = float(_trapz(N * v_cl, D))
    integ_raw = float(_trapz(N * v_raw, D))
    band = (D >= AU_VALID_MM[0]) & (D <= AU_VALID_MM[1])
    integ_band = float(_trapz((N * v_cl)[band], D[band]))
    vol_flux = float(_trapz(N * (np.pi / 6.0) * D ** 3 * v_cl, D))
    # volume flux [mm^3 m^-2 s^-1] -> mm/h: 1 mm^3 / m^2 = 1e-6 m^3/1e4 cm^2...
    # direct check: mm/h = 3.6e6 * (m^3 m^-2 s^-1) with D in mm -> 1e-9 m^3
    rate_check = float(vol_flux * 1e-9 * 3.6e6)
    return {
        "rate_mm_h": R,
        "lambda_per_mm": dsd_lambda(R),
        "D0_mm": mass_median_diameter(R),
        "F_clamped": integ_cl,
        "F_unclamped": integ_raw,
        "F_validity_band_0p5_5mm": integ_band,
        "F_table": F_TABLE.get(R),
        "F_ratio_clamped_to_table": (integ_cl / F_TABLE[R]) if R in F_TABLE else None,
        "rate_back_check_mm_h": rate_check,
        "rate_back_check_rel": (abs(rate_check - R) / R),
        "note": ("F_clamped is the physical variant and the gate input; "
                 "F_unclamped reproduces the report table; the 0.5-5 mm band "
                 "restriction gives the report's own 'factor ~2' caveat "
                 "concrete"),
        "quadrature": {"n": int(n), "D_range_mm": [0.0, DSD_DMAX_MM]},
    }


def mean_impact_energy(rate_mm_h, n=240001):
    """Flux-weighted mean impact kinetic energy  E[1/2 m v^2]  [J/drop].

    m = (pi/6) rho_w D^3 (spherical drop, report section 3.2).  Flux
    weighting (weight N v dD) is the correct average: an impact count is
    weighted by how many drops of that size actually arrive per second.
    """
    R = float(rate_mm_h)
    D = _dsd_grid(n)
    N = dsd_number_density(R, D)
    v = np.maximum(fall_speed(D), 0.0)
    m = (np.pi / 6.0) * RHO_W * (D * 1e-3) ** 3
    e = 0.5 * m * v * v
    w = N * v
    flux_mean = float(_trapz(w * e, D) / _trapz(w, D))
    num_mean = float(_trapz(N * e, D) / _trapz(N, D))
    return {
        "rate_mm_h": R,
        "E_flux_weighted_J": flux_mean,
        "E_number_weighted_J": num_mean,
        "E_flux_weighted_uJ": flux_mean * 1e6,
        "mean_fall_speed_m_s_flux_weighted": float(
            _trapz(w * v, D) / _trapz(w, D)),
        "note": ("E_flux_weighted is the mean energy per IMPACT and is what "
                 "the variance budget uses"),
    }


def ring_wavelength_mm(rate_mm_h):
    """Monotone ring-scale placement anchored to the observed 1-2 cm band.

    See the module docstring: lambda = 10 mm + 10 mm * t, with
    t = (D0(R) - D0(R=1)) / (D0(R=50) - D0(R=1)).  Monotone increasing in R
    (bigger drops -> longer rings, report section 3.4).
    """
    R = float(rate_mm_h)
    d0 = mass_median_diameter(R)
    d_lo = mass_median_diameter(RING_RATE_ANCHORS[0])
    d_hi = mass_median_diameter(RING_RATE_ANCHORS[1])
    t = (d0 - d_lo) / (d_hi - d_lo)
    lo, hi = RING_LAMBDA_BAND_MM
    lam = lo + (hi - lo) * t
    clipped = bool(lam < lo or lam > hi)
    lam_c = min(max(lam, lo), hi)
    return {
        "rate_mm_h": R, "D0_mm": d0, "t_placement": float(t),
        "lambda_ring_mm": float(lam),
        "lambda_ring_mm_clipped": float(lam_c),
        "clipped_to_band": clipped,
        "rule": ("lambda = 10mm + 10mm * (D0-D0(R=1))/(D0(R=50)-D0(R=1)); "
                 "the two observed band limits are the anchors, D0(R) is the "
                 "interpolation variable (report section 3.4: bigger drops "
                 "-> longer rings).  Modelling placement, not a calibration."),
    }


def ring_wavenumber(rate_mm_h):
    """k_ring = 2 pi / lambda_ring  [rad/m]; monotone decreasing in R."""
    w = ring_wavelength_mm(rate_mm_h)
    lam_m = w["lambda_ring_mm_clipped"] * 1e-3
    k = 2.0 * np.pi / lam_m
    return {
        "rate_mm_h": float(rate_mm_h),
        "k_ring_rad_per_m": float(k),
        "lambda_ring_m": float(lam_m),
        "in_band_310_630": bool(RING_K_BAND[0] <= k <= RING_K_BAND[1]),
        "design_reference_k_420": RING_K_DESIGN_REF,
        "placement": w,
    }


def gravity_capillary_omega(k_rad_per_m):
    """omega^2 = g k + (sigma/rho) k^3  [rad/s] (report sections 3.4 / 6.1)."""
    k = np.asarray(k_rad_per_m, np.float64)
    return np.sqrt(np.maximum(G * k + (SIGMA_SURF / RHO_W) * k ** 3, 0.0))


def capillary_crossover():
    """l_c = sqrt(sigma/(rho g)), k_c, lambda_c (report sections 3.4 / 6.1)."""
    l_c = math.sqrt(SIGMA_SURF / (RHO_W * G))
    k_c = 1.0 / l_c
    return {"l_c_m": l_c, "k_c_rad_per_m": k_c,
            "lambda_c_m": 2.0 * math.pi * l_c,
            "c_min_m_s": (4.0 * G * SIGMA_SURF / RHO_W) ** 0.25}


def ring_wave_energy_density_per_variance(k_rad_per_m):
    """(rho g + sigma k^2): energy per unit area per unit height variance.

    For a linear wave train the energy density is E = (rho g + sigma k^2) a^2
    / 2 with <h^2> = a^2/2, so <h^2> = E / (rho g + sigma k^2).
    """
    k = float(k_rad_per_m)
    return RHO_W * G + SIGMA_SURF * k * k


def ring_band_variance(rate_mm_h, eta=RING_ETA):
    """Total ring-band height variance sigma_ring^2 [m^2] (PIPELINE-V6 sec. 2).

    variance ~ eta * F(R) * Ebar(R) * tau  is an energy per unit area
    [J m^-2]; dividing by (rho g + sigma k_ring^2) turns it into a height
    variance.  tau is the ~150 ms impact-process lifetime of report 3.4.
    """
    k = ring_wavenumber(rate_mm_h)["k_ring_rad_per_m"]
    fl = impact_flux(rate_mm_h)
    en = mean_impact_energy(rate_mm_h)
    e_ring = float(eta) * fl["F_clamped"] * en["E_flux_weighted_J"] * RING_TAU_S
    denom = ring_wave_energy_density_per_variance(k)
    var = e_ring / denom
    return {
        "rate_mm_h": float(rate_mm_h),
        "eta": float(eta),
        "F_clamped": fl["F_clamped"],
        "E_flux_weighted_J": en["E_flux_weighted_J"],
        "tau_s": RING_TAU_S,
        "power_in_W_m2": float(fl["F_clamped"] * en["E_flux_weighted_J"]),
        "energy_per_area_J_m2": e_ring,
        "rho_g_plus_sigma_k2": float(denom),
        "ring_variance_m2": float(var),
        "ring_rms_height_mm": float(math.sqrt(max(var, 0.0)) * 1e3),
        "ring_rms_slope_rad": float(math.sqrt(max(var, 0.0)) * k),
        "k_ring_rad_per_m": float(k),
    }


def ring_sigma_k(k_ring, frac=RING_SIGMA_K_FRAC):
    """Gaussian half-width of the bump in k, as a fraction of the band width."""
    lo, hi = RING_K_BAND
    return float(max((hi - lo) * float(frac), 1e-6))


def gravity_band_damping(rate_mm_h, m0_lowk_m2, tau_s=RING_TAU_S):
    """Second rain effect (PIPELINE-V6 section 4): long-wave damping factor.

    The design asks only for the MECHANISM, with a coefficient close to 1 and
    the shape documented as an approximation.  No new calibration constant is
    introduced: the rate is the first-order estimate

        gamma = (rain impact power in [W/m^2]) / (gravity-band energy [J/m^2])
              = F * Ebar / (rho g m0_lowk)

    i.e. "how long it would take the rain's kinetic-energy flux to replace the
    energy of the sub-band (gravity) waves", and the factor is evaluated over
    ONE impact-process lifetime tau (the same 150 ms of report 3.4), which is
    the process time the ring band itself is built on:

        gain(k) = 1 / (1 + gamma * tau * w(k)),  w(k) = 1/(1+(k/310)^4)

    w(k) -> 0 above the ring band, so the capillary band is untouched
    (ring band up, gravity band down -- report section 3.4).  The per-8 s-loop
    value of the same rate is also returned but is NOT applied: the state is a
    static spectrum, and the design's "close to 1" coefficient refers to the
    in-scene effect.
    """
    fl = impact_flux(rate_mm_h)
    en = mean_impact_energy(rate_mm_h)
    p_in = float(fl["F_clamped"] * en["E_flux_weighted_J"])
    e_grav = float(RHO_W * G * max(float(m0_lowk_m2), 1e-15))
    gamma = p_in / e_grav
    return {
        "rate_mm_h": float(rate_mm_h),
        "power_in_W_m2": p_in,
        "gravity_band_energy_J_m2": e_grav,
        "gamma_per_s": float(gamma),
        "tau_s": float(tau_s),
        "gamma_tau": float(gamma * tau_s),
        "gain_applied": float(1.0 / (1.0 + gamma * tau_s)),
        "gain_per_8s_loop_not_applied": float(math.exp(-gamma * 8.0)),
        "weight": "w(k)=1/(1+(k/310)^4) on the sub-band side only",
        "note": ("approximate FORM permitted by PIPELINE-V6 section 4; the "
                 "magnitude is the first-order rate above, with no tuned "
                 "constant.  In the windless verification scene there is "
                 "nothing but the Cox-Munk floor to damp, so this is close "
                 "to a no-op (recorded)."),
    }


def rain_state_summary(rate_mm_h, m0_lowk_m2=1.0e-4, eta=RING_ETA):
    """Everything the state builder needs for one rain rate."""
    R = float(rate_mm_h)
    if R <= 0.0:
        return {"rate_mm_h": 0.0, "active": False,
                "note": "R=0 -> no ring band, no damping: v5 behaviour exactly"}
    k = ring_wavenumber(R)
    var = ring_band_variance(R, eta=eta)
    fl = impact_flux(R)
    en = mean_impact_energy(R)
    damp = gravity_band_damping(R, m0_lowk_m2)
    cw = capillary_crossover()
    return {
        "rate_mm_h": R,
        "active": True,
        "lambda_per_mm": fl["lambda_per_mm"],
        "lambda_per_m": fl["lambda_per_mm"] * 1e3,
        "D_mean_mm": mean_diameter(R),
        "D0_mm": fl["D0_mm"],
        "k_ring_rad_per_m": k["k_ring_rad_per_m"],
        "lambda_ring_m": k["lambda_ring_m"],
        "k_ring_in_band": k["in_band_310_630"],
        "placement": k["placement"],
        "flux": fl,
        "energy": en,
        "band": var,
        "damping": damp,
        "capillary_crossover": cw,
        "physics_chain": ("R -> Lambda(R) -> N(D) -> v(D) -> F -> Ebar -> "
                          "D0 -> k_ring -> sigma_ring^2 (eta is the only "
                          "calibrated factor)"),
    }


def eta_sensitivity(rates=(2.0, 10.0, 50.0), etas=(0.03, 0.10, 0.30)):
    """Ring-band slope variance for a few eta values (honest sensitivity)."""
    out = {}
    for R in rates:
        out["R%g" % R] = {
            "eta_%g" % e: {
                "ring_variance_m2": ring_band_variance(R, eta=e)["ring_variance_m2"],
                "ring_rms_slope_rad": ring_band_variance(R, eta=e)["ring_rms_slope_rad"],
                "ring_rms_height_mm": ring_band_variance(R, eta=e)["ring_rms_height_mm"],
            } for e in etas}
    return out


def flux_gate(rates=(1.0, 10.0, 25.0, 50.0)):
    """PIPELINE-V6 section 6 gate 2: F(R) vs the report section 3.3 table."""
    rows = {}
    ok = True
    for R in rates:
        fl = impact_flux(R)
        ratio = fl["F_clamped"] / F_TABLE[R]
        passed = bool(1.0 / F_TABLE_TOL_FACTOR <= ratio <= F_TABLE_TOL_FACTOR)
        ok = ok and passed
        rows["R%g" % R] = {
            "F_clamped": fl["F_clamped"],
            "F_unclamped": fl["F_unclamped"],
            "F_table": F_TABLE[R],
            "ratio": float(ratio),
            "dev_rel": float(ratio - 1.0),
            "pass_factor_1p5": passed,
            "F_validity_band_0p5_5mm": fl["F_validity_band_0p5_5mm"],
            "rate_back_check_mm_h": fl["rate_back_check_mm_h"],
        }
    return {"gate": "F(R) integral vs report section 3.3 table, factor 1.5",
            "tolerance_factor": F_TABLE_TOL_FACTOR,
            "rows": rows, "pass": bool(ok)}


def wavenumber_gate(rates=None):
    """PIPELINE-V6 section 6 gate 3: k_ring monotone decreasing, in band."""
    if rates is None:
        rates = [1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0,
                 30.0, 40.0, 50.0]
    rows = []
    prev = None
    mono = True
    inside = True
    for R in rates:
        k = ring_wavenumber(R)
        kk = k["k_ring_rad_per_m"]
        if prev is not None and not (kk < prev):
            mono = False
        prev = kk
        if not k["in_band_310_630"]:
            inside = False
        rows.append({"rate_mm_h": R, "D0_mm": k["placement"]["D0_mm"],
                     "lambda_ring_mm": k["lambda_ring_m"] * 1e3,
                     "k_ring_rad_per_m": kk,
                     "in_band": k["in_band_310_630"]})
    return {"gate": "k_ring(R) strictly decreasing and inside 310-630 rad/m",
            "band_rad_per_m": list(RING_K_BAND), "rows": rows,
            "strictly_decreasing": bool(mono), "all_inside_band": bool(inside),
            "pass": bool(mono and inside)}


if __name__ == "__main__":
    cw = capillary_crossover()
    print("crossover l_c=%.4f m  k_c=%.1f  lambda_c=%.4f m  c_min=%.4f m/s"
          % (cw["l_c_m"], cw["k_c_rad_per_m"], cw["lambda_c_m"],
             cw["c_min_m_s"]))
    for R in (1.0, 2.0, 10.0, 25.0, 50.0):
        s = rain_state_summary(R)
        print("R=%-5g Lambda=%.3f D0=%.3fmm F=%.0f (tab %s) Ebar=%.2f uJ "
              "k=%.1f lam=%.2fcm var=%.3e m2 rms_h=%.3f mm rms_slope=%.4f"
              % (R, s["lambda_per_mm"], s["D0_mm"], s["flux"]["F_clamped"],
                 s["flux"]["F_table"], s["energy"]["E_flux_weighted_uJ"],
                 s["k_ring_rad_per_m"], s["lambda_ring_m"] * 100,
                 s["band"]["ring_variance_m2"],
                 s["band"]["ring_rms_height_mm"],
                 s["band"]["ring_rms_slope_rad"]))
    g2 = flux_gate()
    print("F gate:", g2["pass"],
          {k: round(v["ratio"], 3) for k, v in g2["rows"].items()})
    g3 = wavenumber_gate()
    print("k gate:", g3["pass"], "mono=%s inside=%s"
          % (g3["strictly_decreasing"], g3["all_inside_band"]))
    print("eta sensitivity:", eta_sensitivity())
