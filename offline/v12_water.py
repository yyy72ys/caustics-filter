"""v12_water -- version 3.6 chunk 6: in-water transport, a(lambda)/b(lambda)
split, Henyey-Greenstein single scattering.

What this file is
-----------------
The water path of the pipeline before this chunk attenuated the direct beam
with ONE per-channel Beer-Lambert factor driven by a turbidity scalar, and
softened the caustic master with an empirical Gaussian of
``sigma = BLUR_K_PX * sqrt(b * z)`` (v5_physics, imported unchanged).  This
module replaces that model with

  1. a(lambda) and b(lambda) SEPARATED, driven by a particle-type preset
     (clear / algae-rich pond / sediment) or by explicit spectra.  Every
     value carries its provenance; values that are modelling choices are
     labelled ``modelling_choice`` and are checked against the published
     range they were placed inside.
  2. SINGLE SCATTERING in the transport: the ballistic beam keeps the
     refracted caustic and the attenuation exp(-c L); the single-scattered
     light is a halo whose spatial spread is the second moment of the
     scattering integral (the derivation of ``sigma`` is below and is
     repeated in REPORT-V12-GATES.md).  The empirical blur becomes a
     numerical representation whose sigma comes out of the physics.
  3. a CDOM PORT: an additive term in a(lambda), default zero.  The port is
     built; the model that would fill it is v5's domain and is not attempted.

Nothing in v5-v11 is edited; this module imports them.

Geometry / conventions
----------------------
* channel wavelengths 650/550/470 nm (the pipeline's own dispersion set,
  v5_physics.CHANNEL_WL_NM -- imported, not re-declared)
* depths in m, coefficients in 1/m, phase-function angles in rad
* L = s_down + s_up is the total in-water path of the direct beam
  (s_down = z/cos(theta_w) after Snell, s_up = z/cos(theta_v); the offline
  composite is a rectified top view so theta_v = 0)
* the master C is ACHROMATIC (representative G band, det.py); it is blurred
  into the halo, and the halo's per-channel weight carries the spectral
  information.  That is stated where it matters.

Derivation of the single-scatter halo width (formula-level checks live in
``checks_*`` below)
-----------------------------------------------------------------------
Down the path of length L the probability that a photon is still ballistic
after a distance s is exp(-c s), c = a + b.  The density of FIRST scattering
events at s is b exp(-c s).  After scattering by theta the photon covers the
remaining path (L - s) at that angle, so its lateral displacement is

    r = (L - s) * tan(theta),        r^2 = (L - s)^2 tan^2(theta).

The total single-scatter weight and the mean square displacement are

    H  = int_0^L b exp(-c s) ds = (b/c) (1 - exp(-c L))
    <r^2> = int_0^L b exp(-c s) (L - s)^2 <tan^2 theta> ds
          = <theta^2> * (b/c^3) * D2(c L),
    D2(x) = x^2 - 2x + 2 - 2 exp(-x)

(the substitution u = L - s makes the integral elementary).  For an
isotropic-in-plane halo <r^2> = 2 sigma^2, and tan^2 theta -> theta^2 in the
small-angle limit, so

    sigma^2 = (<theta^2>/2) * (b / c^3) * D2(c L)          [m^2]

<theta^2> is MEASURED from the Henyey-Greenstein phase function by
quadrature (``theta2_hg``); at g = 0.85 it is 0.368 rad^2, i.e. 23% above the
second-order identity 2(1-g) = 0.300 rad^2, and that measured number is what
the width uses.  Limits:

    sigma_thin  (c L -> 0):  sigma^2 -> (<theta^2>/2) b L^3 / 3
    sigma_thick (c L large): the WIDTH keeps growing like L (the halo of a
                             thick medium is the remaining path times the
                             angular spread); what saturates is the WEIGHT,
                             H -> b/c.

The empirical form it replaces, sigma ~ sqrt(b z), has the right b-scaling
in the thin limit (sigma^2 ~ b) but the wrong z-scaling (z^3, not z^2) and
a free constant (BLUR_K_PX = 40 px) where this derivation has geometry
(sigma_px = sigma_m / dx).  Both statements are measured in stage 1.

Energy bookkeeping (single scattering only)
-------------------------------------------
    T_bal = exp(-c L)                          ballistic
    H     = (b/c)(1 - T_bal) * rho             scattered, returned
    rho   = exp(-a * <Delta s>), <Delta s> = (<theta^2>/2) * <L - s>
    A     = 1 - T_bal - H                      absorbed in the first event
                                               OR multiple-scattered (not
                                               modelled here)
The identity T_bal + H + A = 1 holds exactly by construction and is checked.
``A`` mixes absorption with the multiple-scattering loss: PIPELINE-V8 line 53
keeps multiple scattering out of this chunk, so the halo is deliberately a
LOWER bound on the returned light in strongly scattering water.
"""
import math

import numpy as np

import v5_physics as V5
import v5_composite as V5C
from transport import gaussian_blur_fft

# --------------------------------------------------------------------------
# channel model and literature anchors
# --------------------------------------------------------------------------
WL = V5.CHANNEL_WL_NM                      # (650, 550, 470) nm -- not re-declared
WL_NOTE = V5.CHANNEL_WL_NOTE

#: Morel (1974) pure-water scattering: b_w(500) = 0.00288 /m, exponent 4.32.
#: The legacy path (v5_physics.BW_550 = 0.0018 at 550 nm) is kept untouched for
#: the b = 0 regression; the two anchors differ by 6% at 550 nm and that
#: difference is reported, never silently swapped.
BW_MOREL_500 = 0.00288
BW_MOREL_EXP = 4.32

#: CDOM absorption slope (Bricaud, Morel & Prieur 1981: S ~ 0.014 /nm;
#: Twardowski et al. 2004 put the natural-water range at 0.010-0.020 /nm).
CDOM_S_NM = 0.014
CDOM_S_RANGE_NM = (0.010, 0.020)

#: Henyey-Greenstein asymmetry.  The design text (PIPELINE-V8 chunk 6) names
#: g ~ 0.85 for freshwater forward scattering as the reference candidate;
#: Petzold (1972) measured VSFs of natural waters imply g ~ 0.85-0.95 and HG
#: is known to under-fit the near-forward peak (Mobley et al. 2002 use two
#: HG terms).  g is a recipe field; the sensitivity of sigma to it is measured
#: in stage 1 and reported rather than tuned.
G_DEFAULT = 0.85
G_BAND = (0.80, 0.95)

# --------------------------------------------------------------------------
# particle-type presets
# --------------------------------------------------------------------------
# Each preset stores the PARTICULATE and DISSOLVED parts at a reference load
# of 1.0, so the recipe can scale the load without changing the a/b RATIO
# (the ratio is the preset's identity -- PIPELINE-V8 chunk 6).  Pure water is
# added by ``preset_spectra`` from the literature tables, never stored here.
#
# ``a_part`` / ``b_part`` are modelling choices PLACED INSIDE the published
# order-of-magnitude ranges for that water class (recorded per preset in
# ``ranges``); they are not measured spectra.  The only measured numbers in
# the whole table are the pure-water parts, the CDOM exponential slope (S) and
# the phase-function asymmetry (g).  ``checks_presets`` verifies that every
# preset lands inside its declared range at load 1 and at its default load,
# and that a >= 0, b >= 0 over the usable load interval.
PRESETS = {
    "clear": {
        "label": "clear freshwater (oligotrophic lake)",
        "a_part": [0.015, 0.020, 0.035],
        "b_part": [0.045, 0.055, 0.075],
        "cdom_a440": 0.10,
        "g": G_DEFAULT,
        "load_default": 1.0,
        "ranges": {
            "a550_range": (0.05, 0.50),
            "b550_range": (0.02, 0.50),
            "a550_range_source": ("order-of-magnitude range for clear "
                                  "lakes, design research (Kirk 2011 ch.3 "
                                  "water-type tables), NOT re-verified in "
                                  "this session"),
            "b550_range_source": "same",
        },
    },
    "algae_pond": {
        "label": "algae-rich pond (eutrophic, green water)",
        "a_part": [0.35, 0.55, 1.80],
        "b_part": [1.20, 1.40, 1.80],
        "cdom_a440": 1.20,
        "g": G_DEFAULT,
        "load_default": 0.6,
        "ranges": {
            "a550_range": (0.50, 10.0),
            "b550_range": (0.50, 10.0),
            "a550_range_source": ("eutrophic / productive inland water, "
                                  "design research (Kirk 2011 ch.3), NOT "
                                  "re-verified in this session"),
            "b550_range_source": "same",
        },
    },
    "sediment": {
        "label": "sediment-laden (resuspended mineral particles)",
        "a_part": [0.10, 0.12, 0.18],
        "b_part": [2.50, 3.00, 3.50],
        "cdom_a440": 0.05,
        "g": G_DEFAULT,
        "load_default": 0.25,
        "ranges": {
            "a550_range": (0.03, 3.0),
            "b550_range": (0.50, 100.0),
            "a550_range_source": ("mineral-dominated turbid inland water, "
                                  "design research (Babin et al. 2003 SPM "
                                  "IOP relations; Kirk 2011), NOT "
                                  "re-verified in this session"),
            "b550_range_source": "same",
        },
    },
}

#: Which preset a medium declaration gets when it names none.  Stage 4 of
#: run_v12 flips this ONLY if G1-G4 pass (PIPELINE-V8 chunk 6 merge rule).
DEFAULT_MEDIUM_PRESET = "algae_pond"

#: The route taken by a recipe WITHOUT a medium block.  This is a
#: compatibility invariant (G1), not a merge candidate: it is never flipped.
LEGACY_RECIPE_ROUTE = "legacy"

#: THE MERGE SWITCH (PIPELINE-V8 chunk 6 merge rule, stage 4 of run_v12).
#: "legacy" = a recipe that carries a water block but names no medium keeps
#: the v5 path; "ab" = such a recipe gets DEFAULT_MEDIUM_PRESET through the a/b
#: transport.  Recipes with NO water section at all are NOT affected: they stay
#: on the v5 path bit for bit (that is the G1 compatibility invariant, and it
#: is why the merge can be flipped without touching any existing recipe).
#: Stage 4 rewrites this one line, atomically, iff G1-G4 all pass.
WATER_MODEL_DEFAULT = "ab"


# --------------------------------------------------------------------------
# pure-water and dissolved parts
# --------------------------------------------------------------------------
def a_water(wl=WL):
    """Pope & Fry (1997) a_w(lambda) -- v5_physics's table, imported."""
    return np.asarray(V5.pure_water_absorption(wl), dtype=np.float64)


def b_water_morel(wl=WL):
    """Morel (1974): b_w(lambda) = 0.00288 * (500/lambda)^4.32 (1/m)."""
    lam = np.asarray(wl, dtype=np.float64)
    return BW_MOREL_500 * (500.0 / lam) ** BW_MOREL_EXP


def b_water_legacy(wl=WL):
    """The pre-3.6 anchor, for the b = 0 regression only (v5_physics)."""
    return np.asarray(V5.pure_water_scattering(wl), dtype=np.float64)


def a_cdom(a440, wl=WL, slope_nm=CDOM_S_NM):
    """CDOM absorption port: a_g(lambda) = a_g(440) exp(-S (lambda-440)).

    The port is additive and defaults to ZERO.  Filling it is v5's business;
    this chunk only guarantees that the slot exists, is spectral, and that a
    zero fill leaves everything else untouched bit for bit.
    """
    lam = np.asarray(wl, dtype=np.float64)
    return float(a440) * np.exp(-float(slope_nm) * (lam - 440.0))


# --------------------------------------------------------------------------
# preset -> a/b spectra
# --------------------------------------------------------------------------
def resolve_preset(name):
    key = str(name or DEFAULT_MEDIUM_PRESET)
    if key not in PRESETS:
        raise ValueError("unknown water preset %r (expected one of %s)"
                         % (name, sorted(PRESETS)))
    return key


def preset_spectra(name, load=None, g=None, cdom_a440=None, wl=WL):
    """a(lambda), b(lambda) for a particle preset at a load.  Returns a dict
    with the parts kept separate so the report can attribute every 1/m."""
    key = resolve_preset(name)
    P = PRESETS[key]
    ld = float(P["load_default"] if load is None else load)
    if ld < 0.0:
        raise ValueError("water load must be >= 0 (got %r)" % (load,))
    a_p = np.asarray(P["a_part"], dtype=np.float64) * ld
    b_p = np.asarray(P["b_part"], dtype=np.float64) * ld
    c440 = float(P["cdom_a440"] if cdom_a440 is None else cdom_a440) * ld
    a_w = a_water(wl)
    b_w = b_water_morel(wl)
    a_g = a_cdom(c440, wl)
    a = a_w + a_p + a_g
    b = b_w + b_p
    gg = float(P["g"] if g is None else g)
    return {
        "preset": key, "label": P["label"], "load": ld,
        "g": gg, "cdom_a440": c440, "cdom_s_nm": CDOM_S_NM,
        "wl_nm": list(wl),
        "a_w": a_w, "a_particle": a_p, "a_cdom": a_g, "a": a,
        "b_w": b_w, "b_particle": b_p, "b": b,
        "c": a + b,
        "a_over_c": (a / np.maximum(a + b, 1e-30)).tolist(),
        "g_band": list(G_BAND),
        "provenance": {
            "a_w": "Pope & Fry (1997) table, v5_physics (literature)",
            "b_w": ("Morel (1974) b_w(500)=0.00288 /m, exponent 4.32 "
                    "(literature)"),
            "cdom_slope": ("Bricaud et al. (1981) S=0.014 /nm, natural-water "
                           "range 0.010-0.020 /nm (literature)"),
            "a_particle": "modelling_choice inside the preset's declared range",
            "b_particle": "modelling_choice inside the preset's declared range",
            "cdom_a440": ("modelling_choice (dissolved load); the MODEL that "
                          "would predict it is v5's domain -- port only"),
            "g": ("design text reference value 0.85, inside the Petzold "
                  "(1972) band 0.85-0.95; HG under-fits the forward peak"),
        },
    }


def explicit_spectra(a_rgb, b_rgb, g=G_DEFAULT, cdom_a440=0.0, wl=WL):
    """A recipe that states the spectra directly (absolute values, 1/m).

    ``a_rgb``/``b_rgb`` REPLACE a and b wholesale: pure water is the caller's
    responsibility.  The values are clipped at 0 and reported as clipped.
    """
    a = np.clip(np.asarray(a_rgb, dtype=np.float64), 0.0, None)
    b = np.clip(np.asarray(b_rgb, dtype=np.float64), 0.0, None)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("explicit spectra need 3 channels (R,G,B)")
    return {
        "preset": None, "label": "explicit spectra", "load": 1.0,
        "g": float(g), "cdom_a440": float(cdom_a440), "cdom_s_nm": CDOM_S_NM,
        "wl_nm": list(wl),
        "a_w": a_water(wl), "a_particle": a, "a_cdom": a_cdom(cdom_a440, wl),
        "a": a, "b_w": b_water_morel(wl), "b_particle": b, "b": b,
        "c": a + b,
        "a_over_c": (a / np.maximum(a + b, 1e-30)).tolist(),
        "g_band": list(G_BAND),
        "provenance": {
            "a": "explicit recipe spectra (absolute)",
            "b": "explicit recipe spectra (absolute)",
            "note": ("explicit spectra include (or omit) pure water by the "
                     "caller's choice; nothing is added behind the scenes"),
        },
    }


def parse_medium(raw):
    """Read ``water.medium`` from a recipe.  Absent -> the legacy route.

    Schema (all optional except ``model``)::

        water.medium = {
          "model": "ab" | "legacy",       # default: absent -> legacy
          "preset": "clear"|"algae_pond"|"sediment",
          "load": 0.25,                   # particle-load scalar
          "g": 0.85,                      # HG asymmetry
          "cdomA440": 0.0,                # CDOM port, default 0
          "spectra": {"a": [..3..], "b": [..3..]}   # absolute, replaces preset
        }

    A legacy recipe has no ``water.medium`` at all and takes the unchanged
    v5 route (G1).  A recipe with a water block but no medium block follows
    ``WATER_MODEL_DEFAULT`` -- the merge switch.
    """
    water = raw.get("water")
    node = (water or {}).get("medium") if isinstance(water, dict) else None
    if not isinstance(node, dict):
        if not isinstance(water, dict):
            # No water section at all -> the compatibility lane, never flipped.
            return {"model": LEGACY_RECIPE_ROUTE, "declared": False,
                    "route": "no water section (compatibility invariant, G1)"}
        if WATER_MODEL_DEFAULT != "ab":
            return {"model": LEGACY_RECIPE_ROUTE, "declared": False,
                    "route": "water section, no medium block, default=legacy"}
        spec = preset_spectra(DEFAULT_MEDIUM_PRESET)
        spec["declared"] = False
        spec["source_kind"] = "default preset (merge default)"
        return {"model": "ab", "declared": False, "spec": spec,
                "route": "water section, no medium block, default=ab"}
    model = str(node.get("model", "ab")).strip().lower()
    if model not in ("ab", "legacy"):
        raise ValueError("water.medium.model must be 'ab' or 'legacy'")
    if model == "legacy":
        return {"model": "legacy", "declared": True,
                "note": "explicitly declared legacy route"}
    sp = node.get("spectra")
    if isinstance(sp, dict) and "a" in sp and "b" in sp:
        spec = explicit_spectra(sp["a"], sp["b"],
                                g=node.get("g", G_DEFAULT),
                                cdom_a440=node.get("cdomA440", 0.0))
        spec["declared"] = True
        spec["source_kind"] = "explicit spectra"
        return {"model": "ab", "declared": True, "spec": spec}

    spec = preset_spectra(node.get("preset", DEFAULT_MEDIUM_PRESET),
                          load=node.get("load"),
                          g=node.get("g"),
                          cdom_a440=node.get("cdomA440"))
    spec["declared"] = True
    spec["source_kind"] = "particle preset"
    return {"model": "ab", "declared": True, "spec": spec}


def optics_from_recipe(raw, env, depth_m, grid, patch_m, g=None):
    """The ONE place a recipe is turned into water optics (merge switch).

    ``parse_medium`` decides; this function dispatches.  Callers that want a
    specific lane for an experiment bypass it and call ``legacy_optics`` /
    ``ab_optics`` directly (the gates do), but any production recipe goes
    through here, so flipping ``WATER_MODEL_DEFAULT`` moves every recipe that
    declares a water block and leaves untouched every recipe that does not.
    """
    med = parse_medium(raw)
    out = {"route": med.get("route"), "declared": med.get("declared"),
           "recipe_model": med["model"]}
    if med["model"] == "ab":
        o = ab_optics(med["spec"], env, depth_m, grid, patch_m, g=g)
        o["recipe_model"] = "ab"
        return o
    o = legacy_optics(env, depth_m, grid)
    o["recipe_model"] = "legacy"
    return o


def checks_recipe_switch(grid=1024, patch_m=2.4, default=None):
    """Stage-4 side: the merge switch does what the merge rule says.

    Three recipes are pushed through ``optics_from_recipe`` with a synthetic
    env: (a) no water block, (b) water block without a medium, (c) water block
    with an explicit ``model: legacy``.  The gate is that (a) and (c) are the
    v5 lane whatever ``WATER_MODEL_DEFAULT`` is, and that (b) follows the
    switch.  Nothing here renders.

    ``default`` lets stage 4 test the value it has just written to the file:
    the module-level constant in the running process was bound at import time
    and would otherwise describe the state BEFORE the flip.
    """
    global WATER_MODEL_DEFAULT
    saved = WATER_MODEL_DEFAULT
    if default is not None:
        WATER_MODEL_DEFAULT = str(default)
    env = {"sun": {"altitudeDeg": 59.0265833497501}, "depthM": 0.6,
           "turbid": 0.35, "sunI": 0.8, "diffuse": 0.3, "strength": 1.0,
           "exposure": 1.0, "scatter": [0.02, 0.03, 0.05],
           "sunColor": [1.0, 0.98, 0.94], "absorb": [0.4, 0.2, 0.1],
           "windSpeedMS": 1.5,
           "sunDir": [0.0, 0.0, 1.0]}
    rows = {}
    try:
        for tag, raw in (
                ("no_water_section", {"name": "a"}),
                ("water_no_medium", {"name": "b", "water": {"depthM": 0.6}}),
                ("explicit_legacy", {"name": "c",
                                     "water": {"medium": {"model":
                                                          "legacy"}}})):
            o = optics_from_recipe(raw, env, 0.6, grid, patch_m)
            rows[tag] = {"recipe_model": o["recipe_model"],
                         "optics_mode": o["mode"],
                         "spec_preset": (o["spec"]["preset"]
                                         if o["mode"] == "ab" else None)}
        expected = "legacy" if WATER_MODEL_DEFAULT != "ab" else "ab"
        out = {
            "WATER_MODEL_DEFAULT": WATER_MODEL_DEFAULT,
            "tested_against": ("file value written by stage 4"
                               if default is not None else "module value"),
            "rows": rows,
            "expected_water_no_medium": expected,
            "pass_no_water_stays_legacy": bool(
                rows["no_water_section"]["optics_mode"] == "legacy"),
            "pass_explicit_legacy": bool(
                rows["explicit_legacy"]["optics_mode"] == "legacy"),
            "pass_default_follows_switch": bool(
                rows["water_no_medium"]["optics_mode"] == expected),
            "pass": bool(rows["no_water_section"]["optics_mode"] == "legacy"
                         and rows["explicit_legacy"]["optics_mode"] == "legacy"
                         and rows["water_no_medium"]["optics_mode"] == expected),
            "note": ("the merge rule can only move recipes that declare a "
                     "water block; a recipe without one keeps the v5 path bit "
                     "for bit (G1), which is what makes the flip safe"),
        }
    finally:
        WATER_MODEL_DEFAULT = saved
    return out


# --------------------------------------------------------------------------
# Henyey-Greenstein phase function and its moments
# --------------------------------------------------------------------------
def hg_phase(mu, g):
    """HG phase function p(mu) in 1/sr, int over the sphere = 1."""
    g = float(g)
    mu = np.clip(np.asarray(mu, dtype=np.float64), -1.0, 1.0)
    return (1.0 - g * g) / (4.0 * math.pi
                            * np.power(1.0 + g * g - 2.0 * g * mu, 1.5))


def hg_moments(g, n_mu=400000):
    """Measure <cos>, <theta^2> and <tan^2>_capped from the phase function.

    The integration is over mu = cos(theta) with dOmega = 2 pi dmu.
    ``tan^2`` has a logarithmic divergence at theta = 90 deg for HG (finite
    phase-function density at mu = 0), so it is reported with a cap; that
    divergence is exactly why the halo width uses the small-angle second
    moment <theta^2> = 2(1-g) instead of the exact tangent.
    """
    g = float(g)
    mu = np.linspace(-1.0, 1.0, int(n_mu))
    p = hg_phase(mu, g)
    w = 2.0 * math.pi * p
    norm = float(np.trapz(w, mu))
    mean_mu = float(np.trapz(w * mu, mu) / norm)
    theta = np.arccos(mu)
    th2 = float(np.trapz(w * theta * theta, mu) / norm)
    caps = {}
    for cap_deg in (20.0, 30.0, 45.0):
        cap = math.radians(cap_deg)
        mask = theta <= cap
        with np.errstate(divide="ignore", invalid="ignore"):
            t2 = np.where(mask, np.tan(theta) ** 2, 0.0)
        caps["tan2_capped_%gdeg" % cap_deg] = float(
            np.trapz(w * t2, mu) / norm)
    return {
        "g": g, "norm_over_sphere": norm,
        "mean_cos_measured": mean_mu,
        "mean_cos_expected": g,
        "mean_cos_abs_err": abs(mean_mu - g),
        "theta2_measured_rad2": th2,
        "theta2_closed_form_2_1_minus_g": 2.0 * (1.0 - g),
        "theta2_rel_diff": float(th2 / (2.0 * (1.0 - g)) - 1.0),
        "tan2_capped": caps,
        "note": ("<theta^2> = 2(1-g) is the small-angle identity used by the "
                 "halo width; <tan^2> is divergent at 90 deg for HG, so the "
                 "capped values are reported as the size of the approximation "
                 "error, not as a usable moment"),
    }


def hg_sample_cos(rng, g, n):
    """Inverse-CDF sample of cos(theta) from HG (standard closed form)."""
    g = float(g)
    u = rng.random(int(n))
    if abs(g) < 1e-9:
        return 2.0 * u - 1.0
    return (1.0 + g * g - ((1.0 - g * g) / (1.0 - g + 2.0 * g * u)) ** 2) \
        / (2.0 * g)


# --------------------------------------------------------------------------
# the single-scatter halo: weight and width
# --------------------------------------------------------------------------
def scatter_position_moments(c, L):
    """Mean scattering position along the path, truncated exponential.

    p(s) = c exp(-c s) / (1 - exp(-c L)),  s in [0, L].
    """
    x = float(c) * float(L)
    if abs(x) < 1e-12:
        return {"s_mean_m": 0.5 * float(L), "L_minus_s_mean_m": 0.5 * float(L),
                "x": x}
    e = math.exp(-x)
    s_mean = (1.0 - e * (1.0 + x)) / (float(c) * (1.0 - e))
    return {"s_mean_m": float(s_mean),
            "L_minus_s_mean_m": float(L) - float(s_mean),
            "x": float(x)}


def theta2_hg(g, n_mu=200000):
    """<theta^2> of the HG phase function, by quadrature (no constant).

    The second-order identity <theta^2> ~ 2(1-g) is only the small-angle
    expansion: at g = 0.85 the measured moment is 23% larger (theta^4/12
    term).  The halo width uses the MEASURED moment; 2(1-g) is reported
    beside it as the size of the small-angle step, not as the value.
    """
    g = float(g)
    mu = np.linspace(-1.0, 1.0, int(n_mu))
    w = 2.0 * math.pi * hg_phase(mu, g)
    norm = float(np.trapz(w, mu))
    theta = np.arccos(mu)
    return float(np.trapz(w * theta * theta, mu) / norm)


def sigma_from_scattering(b, c, L, g):
    """sigma of the single-scatter halo, from the scattering integral.

    sigma^2 = (<theta^2>/2) * (b/c^3) * D2(c L),
    D2(x) = x^2 - 2x + 2 - 2exp(-x),
    <theta^2> measured from the HG phase function (``theta2_hg``).
    Returns m and the dimensionless pieces so the report can print them.
    """
    b = max(float(b), 0.0)
    c = float(c)
    if b <= 0.0 or c <= 0.0:
        return {"sigma_m": 0.0, "b": b, "c": c, "L_m": float(L), "g": float(g),
                "D2": 0.0, "x": 0.0,
                "closed_form": "(<theta^2>/2)(b/c^3) D2(cL) -> 0 for b = 0"}
    x = c * float(L)
    D2 = x * x - 2.0 * x + 2.0 - 2.0 * math.exp(-x)
    th2 = theta2_hg(g)
    sig2 = 0.5 * th2 * (b / (c ** 3)) * D2
    return {"sigma_m": float(math.sqrt(max(sig2, 0.0))), "sigma2_m2": float(sig2),
            "b": b, "c": c, "L_m": float(L), "g": float(g),
            "D2": float(D2), "x": float(x),
            "theta2_rad2": th2,
            "theta2_small_angle_2_1_minus_g": 2.0 * (1.0 - g),
            "theta2_small_angle_rel_diff": float(
                th2 / (2.0 * (1.0 - g)) - 1.0),
            "closed_form": "sigma^2 = (<theta^2>/2)(b/c^3)(x^2-2x+2-2exp(-x))"}


def sigma_thin_limit(b, L, g):
    """sigma^2 -> (<theta^2>/2) b L^3 / 3 as c L -> 0 (scaling check)."""
    return math.sqrt(0.5 * theta2_hg(g) * float(b) * float(L) ** 3 / 3.0)


def sigma_thin_limit(b, L, g):
    """sigma^2 -> (1-g) b L^3 / 3 as c L -> 0 (used by the scaling check)."""
    return math.sqrt((1.0 - float(g)) * float(b) * float(L) ** 3 / 3.0)


def partition(a, b, s_down, s_up, g):
    """Ballistic / single-scatter-halo / lost weights of the direct beam.

    Every channel is independent; a and b are (3,) arrays in 1/m.
    T_bal + H + A = 1 exactly (checked in ``checks_partition``).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    c = np.maximum(a + b, 1e-30)
    L = float(s_down) + float(s_up)
    x = c * L
    T = np.exp(-x)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_b = np.where(c > 0, b / c, 0.0)
    H_raw = frac_b * (1.0 - T)
    # per-channel mean scattering position (truncated exponential at rate c)
    half_th2 = 0.5 * theta2_hg(float(g))
    ds = np.zeros(3, np.float64)
    s_mean = np.zeros(3, np.float64)
    for i in range(3):
        pm = scatter_position_moments(c[i], L)
        s_mean[i] = pm["s_mean_m"]
        ds[i] = half_th2 * pm["L_minus_s_mean_m"]
    rho = np.exp(-a * ds)
    H = H_raw * rho
    A = 1.0 - T - H
    return {
        "a": a, "b": b, "c": c, "g": float(g),
        "s_down_m": float(s_down), "s_up_m": float(s_up), "L_m": L, "x": x,
        "T_ballistic": T, "H_halo": H, "A_lost": A,
        "frac_b_over_c": frac_b,
        "rho_halo_extra_absorption": rho,
        "ds_extra_path_m": ds, "s_mean_m": s_mean,
        "theta2_half": half_th2,
        "down_factor": np.exp(-c * float(s_down)),
        "up_factor": np.exp(-c * float(s_up)),
        "two_way_beer": np.exp(-c * L),
        "identity_max_abs_err": float(np.max(np.abs(T + H + A - 1.0))),
    }


# --------------------------------------------------------------------------
# optics assembly for the composite
# --------------------------------------------------------------------------
def legacy_optics(env, depth_m, grid):
    """The pre-3.6 water optics, byte-for-byte the v5 recipe (G1 lane)."""
    att = V5.water_attenuation(env["turbid"])
    beer = V5.beer_lambert(depth_m, env["sun"]["altitudeDeg"], env["turbid"])
    blur = V5.scattering_blur_sigma_px(att, depth_m, grid)
    return {"mode": "legacy", "att": att, "beer": beer, "blur": blur,
            "depth_m": float(depth_m), "grid": int(grid),
            "code_path": ("v5_physics.water_attenuation + beer_lambert + "
                          "scattering_blur_sigma_px (imported, unchanged)")}


def ab_optics(spec, env, depth_m, grid, patch_m, view_zenith_deg=0.0,
              g=None, force_up_one=False, force_down_one=False):
    """The v3.6 water optics: a/b split + single-scatter halo + derived sigma."""
    g = float(spec["g"] if g is None else g)
    a = np.asarray(spec["a"], dtype=np.float64)
    b = np.asarray(spec["b"], dtype=np.float64)
    z = max(0.0, float(depth_m))
    th_w = V5.refraction_angle_water(env["sun"]["altitudeDeg"])
    th_v = math.radians(max(0.0, min(48.6, float(view_zenith_deg))))
    s_down = z / max(math.cos(th_w), 1e-6)
    s_up = z / max(math.cos(th_v), 1e-6)
    part = partition(a, b, s_down, s_up, g)
    if force_down_one:
        part_down_off = partition(a, b, 0.0, s_up, g)
        w_bal, w_halo = part_down_off["T_ballistic"], part_down_off["H_halo"]
    elif force_up_one:
        part_up_off = partition(a, b, s_down, 0.0, g)
        w_bal, w_halo = part_up_off["T_ballistic"], part_up_off["H_halo"]
    else:
        w_bal, w_halo = part["T_ballistic"], part["H_halo"]
    sig = sigma_from_scattering(float(b[1]), float(part["c"][1]),
                                float(part["L_m"]), g)
    dx = float(patch_m) / float(grid)
    sig_px = sig["sigma_m"] / dx if dx > 0 else 0.0
    return {
        "mode": "ab",
        "spec": spec, "a": a, "b": b, "c": part["c"], "g": g,
        "depth_m": z, "grid": int(grid), "patch_m": float(patch_m),
        "dx_m": dx, "s_down_m": s_down, "s_up_m": s_up,
        "theta_w_deg": float(math.degrees(th_w)),
        "theta_v_deg": float(math.degrees(th_v)),
        "partition": part,
        "w_ballistic": w_bal, "w_halo": w_halo,
        "beer_down": part["down_factor"], "beer_up": part["up_factor"],
        "beer_two_way": part["two_way_beer"],
        "sigma": sig, "sigma_px": float(sig_px),
        "sigma_px_per_channel": [sig["sigma_m"] / dx] * 3,
        "channel_master": "G (the master C is achromatic, det.py)",
        "legs_forced": ("down=1" if force_down_one else
                        ("up=1" if force_up_one else None)),
    }


# --------------------------------------------------------------------------
# the composite: ballistic caustic + single-scatter halo
# --------------------------------------------------------------------------
def render_frame_v12(fields, C, photo_lin, env, W, H, optics,
                     ambient_scale=V5C.A_V3, c_layer_exposure=V5C.E_V3,
                     substrate_gain_rgb=(1.0, 1.0, 1.0), g1=1.0,
                     lowpass_sigma_px=V5C.B_SIGMA_V3, glitter_gain=0.0,
                     return_debug=False):
    """One v3.6 still.

    The legacy lane DELEGATES to v5_composite.render_frame_v5 with the v5
    optics, so a legacy recipe cannot move by construction (G1).

    The a/b lane keeps the OUTPUT-SPEC section 3 separation and only changes
    what multiplies the caustic:

        direct(lambda) = E_base,dir(lambda) * g1 *
                         [ T_bal(lambda) * C(x) + H(lambda) * blur(C)(x) ]

    with T_bal = exp(-c L), H = (b/c)(1-exp(-c L)) rho -- the weights of
    ``partition``.  The blur is the numerical representation of the
    single-scatter halo; its sigma comes from ``sigma_from_scattering``.
    """
    if optics.get("mode") == "legacy":
        return V5C.render_frame_v5(
            fields, C, photo_lin, env, W, H,
            ambient_scale=ambient_scale, c_layer_exposure=c_layer_exposure,
            substrate_gain_rgb=substrate_gain_rgb,
            beer_down=optics["beer"]["down"], beer_up=optics["beer"]["up"],
            g1=g1, blur_sigma_px=optics["blur"]["sigma_px"],
            lowpass_sigma_px=lowpass_sigma_px, glitter_gain=glitter_gain,
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
    A = float(ambient_scale)
    E = float(c_layer_exposure)
    gain = np.asarray(substrate_gain_rgb, dtype=np.float64)

    # ---- single-scatter halo: the ONLY use of the blur ----
    w_bal = np.asarray(optics["w_ballistic"], dtype=np.float64)
    w_halo = np.asarray(optics["w_halo"], dtype=np.float64)
    sig_px = float(optics["sigma_px"])
    Cc = np.asarray(C, dtype=np.float64)
    if sig_px > 0.0 and float(np.max(w_halo)) > 0.0:
        Ch = gaussian_blur_fft(Cc, sig_px).astype(np.float64)
    else:
        Ch = Cc

    gx_s, gy_s = _state_gradients(fields, lowpass_sigma_px)
    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
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

    PH, PW, _ = photo_lin.shape
    col = np.zeros((H, W, 3), dtype=np.float64)
    alb_all = np.zeros((H, W, 3), dtype=np.float64)
    direct_all = np.zeros((H, W, 3), dtype=np.float64)
    irr_all = np.zeros((H, W, 3), dtype=np.float64)
    for bb in range(3):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        alb = _wobbled_alb(photo_lin[:, :, bb:bb + 1], u, v, shx, shy, S)[:, :, 0]
        alb_all[:, :, bb] = alb
        c_bal = _sample_c(Cc, u, v, shx, shy, S, N_buf)
        c_hal = _sample_c(Ch, u, v, shx, shy, S, N_buf)
        irr = float(w_bal[bb]) * c_bal + float(w_halo[bb]) * c_hal
        irr_all[:, :, bb] = irr
        direct = float(Ebase[bb]) * irr * float(g1)
        direct_all[:, :, bb] = alb * float(gain[bb]) * direct * strength * exposure
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
        col += sunColor[None, None, :] * (spark * p * float(glitter_gain))[:, :, None]

    lin = col * exposure
    aces = aces_narkowicz(np.maximum(lin, 0.0))
    out = (np.clip(linear_to_srgb(aces), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    if return_debug:
        return out, {"lin": lin, "aces": aces, "alb": alb_all,
                     "direct_lin": direct_all, "direct_irradiance": irr_all,
                     "w_ballistic": [float(x) for x in w_bal],
                     "w_halo": [float(x) for x in w_halo],
                     "beer_down": [float(x) for x in optics["beer_down"]],
                     "beer_up": [float(x) for x in optics["beer_up"]],
                     "g1": float(g1), "E": E, "A": A,
                     "blur_sigma_px": sig_px,
                     "substrate_gain_rgb": [float(x) for x in gain]}
    return out


# --------------------------------------------------------------------------
# field statistics used by the gates
# --------------------------------------------------------------------------
def effective_c_field(C, optics):
    """The C multiplier the transport actually applies (G channel weights)."""
    C = np.asarray(C, dtype=np.float64)
    if optics.get("mode") != "ab":
        sig = float(optics["blur"]["sigma_px"])
        return gaussian_blur_fft(C, sig) if sig > 0 else C
    wb = float(optics["w_ballistic"][1])
    wh = float(optics["w_halo"][1])
    ch = gaussian_blur_fft(C, float(optics["sigma_px"])) \
        if float(optics["sigma_px"]) > 0 and wh > 0 else C
    return wb * C + wh * ch


def peak_profile_lags(F, lags=(16, 32, 48, 64), top_frac=1e-3):
    """Mean normalized field at a fixed lag from the brightest peaks.

    The classical half-width of a caustic peak is set by the unresolved
    ballistic core at every b level (measured: 4 px in all four G3 cases), so
    it cannot discriminate.  What does discriminate is how much light sits at
    a FIXED distance from the peaks: the halo floods the troughs, so this
    profile rises monotonically toward the mean as b grows.
    """
    f = np.asarray(F, dtype=np.float64)
    m = float(f.mean())
    g = f / m
    n = g.shape[0]
    thr = np.percentile(g, 100.0 * (1.0 - top_frac))
    mask = (g >= thr).astype(np.float64)
    cnt = max(float(mask.sum()), 1.0)
    cor = np.fft.irfft2(np.fft.rfft2(g) * np.conj(np.fft.rfft2(mask)),
                        s=g.shape) / cnt
    prof = {}
    for r in lags:
        prof["r%d_px" % r] = float(0.5 * (cor[0, r] + cor[r, 0]))
    return {"profile": prof, "lag_px": [int(r) for r in lags],
            "n_peaks": int(cnt), "top_frac": float(top_frac),
            "note": ("g = F/mean; 1.0 is the field mean.  A rise means the "
                     "troughs around the peaks are being flooded by the halo")}


def field_stats(F):
    """Peak / contrast / skirt numbers of a caustic field."""
    f = np.asarray(F, dtype=np.float64)
    m = float(f.mean())
    s = float(f.std())
    p = {q: float(np.percentile(f, q)) for q in (50, 90, 99, 99.9, 99.99)}
    kurt = float(np.mean((f - m) ** 4) / max(s ** 4, 1e-30))
    # skirt width: radial autocorrelation half width (periodic, axis-averaged)
    fc = f - m
    n = fc.shape[0]
    Fp = np.abs(np.fft.fft2(fc)) ** 2
    ac = np.real(np.fft.ifft2(Fp))
    ac = ac / ac[0, 0]
    a = 0.5 * (np.abs(ac[0, :n // 2]) + np.abs(ac[:n // 2, 0]))
    lags = np.arange(a.size)
    def _half(level):
        idx = np.where(a < level)[0]
        if idx.size:
            return float(idx[0])
        return None
    # area share above a low shoulder threshold (skirt growth)
    thr = 1.0 + 0.10 * (p[99] / m - 1.0)
    return {
        "mean": m, "std": s, "p50": p[50], "p90": p[90], "p99": p[99],
        "p99_9": p[99.9], "p99_99": p[99.99],
        "max": float(f.max()), "min": float(f.min()),
        "kurtosis": kurt,
        "contrast_p99_over_mean": (float(p[99] / m) if m != 0 else None),
        "contrast_p99_9_over_mean": (float(p[99.9] / m) if m != 0 else None),
        "excess_over_mean_p99": float(p[99] - m),
        "skirt_r50_px": _half(0.5), "skirt_r1e_px": _half(math.exp(-1.0)),
        "skirt_r50_frac_of_halfgrid": (None if _half(0.5) is None
                                       else float(_half(0.5) / (n / 2.0))),
        "shoulder_area_frac": float(np.mean(f > thr)),
        "shoulder_threshold": float(thr),
        "grid": int(n),
        "peak_profile": peak_profile_lags(f),
    }


def chromaticity(rgb):
    x = np.asarray(rgb, dtype=np.float64)
    s = float(x.sum())
    return (x / s).tolist() if s > 0 else [None, None, None]


def lab_from_linear_rgb(rgb):
    """Mean-colour CIELAB (D65) of a linear-RGB triple or an image array."""
    x = np.asarray(rgb, dtype=np.float64)
    if x.ndim >= 2:
        x = x.reshape(-1, 3).mean(axis=0)
    M = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = M @ x
    wn = np.array([0.95047, 1.0, 1.08883])
    t = xyz / wn

    def f(u):
        return np.where(u > (6.0 / 29.0) ** 3,
                        np.cbrt(u), u / (3.0 * (6.0 / 29.0) ** 2) + 4.0 / 29.0)
    fx, fy, fz = f(t)
    return {"L": float(116.0 * fy - 16.0), "a": float(500.0 * (fx - fy)),
            "b": float(200.0 * (fy - fz))}


def delta_e76(lab1, lab2, chroma_only=False):
    d = [lab2[k] - lab1[k] for k in ("L", "a", "b")]
    if chroma_only:
        d[0] = 0.0
    return float(math.sqrt(sum(x * x for x in d)))


# --------------------------------------------------------------------------
# formula-level checks (stage 1 -- every number is measured here)
# --------------------------------------------------------------------------
def checks_hg(g):
    m = hg_moments(g)
    m["theta2_quadrature_convergence"] = {
        "n_1e5": theta2_hg(g, n_mu=100000),
        "n_4e5": theta2_hg(g, n_mu=400000),
    }
    a = m["theta2_quadrature_convergence"]["n_1e5"]
    bq = m["theta2_quadrature_convergence"]["n_4e5"]
    m["theta2_quadrature_rel_diff"] = float(a / bq - 1.0)
    m["pass_mean_cos"] = bool(m["mean_cos_abs_err"] < 1e-4)
    m["pass_quadrature"] = bool(abs(m["theta2_quadrature_rel_diff"]) < 1e-6)
    m["small_angle_step"] = (
        "<theta^2> = %.6f rad^2 vs the second-order identity 2(1-g) = %.6f "
        "rad^2: %+.1f%%.  The width uses the measured moment; 2(1-g) is "
        "reported as the size of the small-angle step."
        % (m["theta2_measured_rad2"],
           m["theta2_closed_form_2_1_minus_g"],
           100.0 * m["theta2_rel_diff"]))
    m["pass"] = bool(m["pass_mean_cos"] and m["pass_quadrature"])
    return m


def checks_sigma_limits(g=G_DEFAULT):
    """Thin-limit scaling exponents and the thick-limit weight."""
    rows = []
    for L in (0.5, 1.25, 4.0, 8.336):
        for b in (1e-3, 1e-2, 0.05, 0.2, 1.0, 4.0):
            c = 0.02 + b
            s = sigma_from_scattering(b, c, L, g)
            thin = sigma_thin_limit(b, L, g)
            rows.append({"L_m": L, "b": b, "c": c, "x": s["x"],
                         "sigma_m": s["sigma_m"],
                         "sigma_thin_limit_m": thin,
                         "ratio_to_thin": (s["sigma_m"] / thin
                                           if thin > 0 else None)})
    # fitted exponents: sigma ~ b^p at fixed L, sigma ~ L^q at fixed b
    def fit(xs, ys):
        xs = np.log(np.asarray(xs)); ys = np.log(np.asarray(ys))
        A = np.vstack([xs, np.ones_like(xs)]).T
        p, _ = np.linalg.lstsq(A, ys, rcond=None)[0]
        return float(p)
    b_lo = [r for r in rows if r["L_m"] == 1.25 and r["b"] <= 0.05]
    p_b = fit([r["b"] for r in b_lo], [r["sigma_m"] for r in b_lo])
    thin_rows = [{"L": L,
                  "sigma_thin": sigma_thin_limit(1e-3, L, g),
                  "sigma_full": sigma_from_scattering(
                      1e-3, 0.02 + 1e-3, L, g)["sigma_m"]}
                 for L in (0.25, 0.5, 1.0, 2.0, 4.0)]
    q_L = fit([r["L"] for r in thin_rows], [r["sigma_thin"] for r in thin_rows])
    return {
        "g": g, "rows": rows,
        "exponent_sigma_vs_b_thin": p_b,
        "exponent_sigma_vs_b_expected": 0.5,
        "exponent_sigma_vs_L_thin": q_L,
        "exponent_sigma_vs_L_expected": 1.5,
        "thin_limit_table": thin_rows,
        "pass_b": bool(abs(p_b - 0.5) < 0.02),
        "pass_L": bool(abs(q_L - 1.5) < 0.02),
        "pass": bool(abs(p_b - 0.5) < 0.02 and abs(q_L - 1.5) < 0.02),
        "note": ("the retired empirical form had sigma ~ sqrt(b z): the b "
                 "exponent matches (0.5, thin limit) but the depth exponent "
                 "is 1.5, not 1 -- that difference is the physics this "
                 "chunk adds"),
    }


def checks_sigma_mc(b=1.0, c=1.2, L=1.25, g=G_DEFAULT, n=400000, seed=12345):
    """Monte-Carlo check of the closed form (scattering position integral).

    theta^2 is drawn from its HG value 2(1-g) (the small-angle model the
    closed form uses); the sample counts the single-scatter weight too:
    <r^2> = H * <(L-s)^2 tan^2> / ... implemented as
    <r^2> = <(L-s)^2> * 2(1-g) * H  with s from the truncated exponential.
    """
    rng = np.random.default_rng(int(seed))
    x = float(c) * float(L)
    e = math.exp(-x)
    u = rng.random(int(n))
    s = -np.log(1.0 - u * (1.0 - e)) / float(c)
    H = (float(b) / float(c)) * (1.0 - e)
    th2 = theta2_hg(float(g))
    r2 = ((float(L) - s) ** 2) * th2
    mean_r2 = float(np.mean(r2)) * H
    sig_mc = math.sqrt(mean_r2 / 2.0)
    sig_cf = sigma_from_scattering(b, c, L, float(g))["sigma_m"]
    return {
        "b": b, "c": c, "L_m": L, "g": float(g), "n": int(n),
        "H_halo": H, "theta2_rad2": th2,
        "sigma_monte_carlo_m": sig_mc, "sigma_closed_form_m": sig_cf,
        "rel_diff": float(sig_mc / sig_cf - 1.0),
        "pass": bool(abs(sig_mc / sig_cf - 1.0) < 0.01),
        "note": ("MC samples the scattering position from the truncated "
                 "exponential and uses the measured <theta^2>; it therefore "
                 "validates the position integral and the weight algebra, "
                 "not the small-angle assumption itself"),
    }


def checks_partition():
    """T + H + A = 1 exactly, and the b = 0 and a = 0 corners."""
    rows = []
    for a in (0.0, 0.05, 0.3, 1.0, 3.0):
        for b in (0.0, 0.01, 0.3, 1.0, 5.0):
            for (sd, su) in ((0.65, 0.6), (4.336, 4.0)):
                p = partition(np.full(3, a), np.full(3, b), sd, su, G_DEFAULT)
                rows.append({"a": a, "b": b, "s_down": sd, "s_up": su,
                             "T": float(p["T_ballistic"][1]),
                             "H": float(p["H_halo"][1]),
                             "A": float(p["A_lost"][1]),
                             "sum": float((p["T_ballistic"]
                                           + p["H_halo"] + p["A_lost"])[1]),
                             "identity_err": p["identity_max_abs_err"]})
    worst = max(abs(r["identity_err"]) for r in rows)
    b0 = [r for r in rows if r["b"] == 0.0]
    a0 = [r for r in rows if r["a"] == 0.0]
    return {
        "rows": rows, "worst_identity_err": float(worst),
        "b_zero_halo_is_zero": bool(all(r["H"] == 0.0 for r in b0)),
        "b_zero_T_equals_expminusaL": bool(all(
            abs(r["T"] - math.exp(-r["a"] * (r["s_down"] + r["s_up"]))) < 1e-12
            for r in b0)),
        "a_zero_T_plus_H_is_one": bool(all(
            abs(r["T"] + r["H"] - 1.0) < 2e-3 for r in a0)),
        "pass": bool(worst < 1e-12),
        "note": ("A is 'absorbed in the first event OR multiple-scattered'; "
                 "with a = 0 the halo returns what the ballistic beam loses, "
                 "up to the rho correction (2e-3 band)"),
    }


def checks_legacy_vs_morel(wl=WL):
    """The two pure-water scattering anchors, side by side (never swapped)."""
    bl = b_water_legacy(wl)
    bm = b_water_morel(wl)
    return {
        "wl_nm": list(wl),
        "b_w_legacy_v5_anchor": bl.tolist(),
        "b_w_morel_1974": bm.tolist(),
        "rel_diff_morel_over_legacy": (bm / bl - 1.0).tolist(),
        "a_w_pope_fry": a_water(wl).tolist(),
        "note": ("the legacy route keeps v5's 0.0018 anchor (G1 must be bit "
                 "identical); the a/b route uses Morel (1974). The 6% "
                 "difference at 550 nm is reported, not reconciled"),
    }


def checks_presets():
    """Every preset inside its declared published range; a,b >= 0."""
    rows = {}
    ok = True
    for name, P in PRESETS.items():
        for tag, load in (("load_1", 1.0), ("load_default", P["load_default"])):
            s = preset_spectra(name, load=load)
            a550 = float(s["a"][1]); b550 = float(s["b"][1])
            lo_a, hi_a = P["ranges"]["a550_range"]
            lo_b, hi_b = P["ranges"]["b550_range"]
            within = bool(lo_a <= a550 <= hi_a and lo_b <= b550 <= hi_b)
            ok = ok and within
            rows["%s/%s" % (name, tag)] = {
                "load": load,
                "a550": a550, "b550": b550,
                "c_rgb": s["c"].tolist(),
                "a_over_c": s["a_over_c"],
                "a550_range": [lo_a, hi_a], "b550_range": [lo_b, hi_b],
                "within_declared_range": within,
                "nonneg": bool(np.all(s["a"] >= 0) and np.all(s["b"] >= 0)),
            }
    return {"rows": rows, "pass": bool(ok),
            "note": ("the ranges are the order-of-magnitude published ranges "
                     "recorded in the design research and NOT re-verified "
                     "online in this session; the preset values themselves "
                     "are modelling choices inside them")}


def checks_beer_two_way():
    """Down/up legs and the a/b model on the z = 0.6 vs z = 4.0 pair (G4
    model side; the render side is measured in run_v12 stage 3)."""
    spec = preset_spectra("clear", load=0.5)
    out = {"preset": "clear", "load": 0.5, "depths": {}}
    fake_env = {"sun": {"altitudeDeg": 67.29}}  # preset-1 sun altitude
    for z in (0.6, 4.0):
        o = ab_optics(spec, fake_env, z, 1024, 4.8)
        p = o["partition"]
        out["depths"]["z%g" % z] = {
            "s_down_m": o["s_down_m"], "s_up_m": o["s_up_m"],
            "down": [float(x) for x in p["down_factor"]],
            "up": [float(x) for x in p["up_factor"]],
            "two_way_beer": [float(x) for x in p["two_way_beer"]],
            "T_ballistic": [float(x) for x in p["T_ballistic"]],
            "H_halo": [float(x) for x in p["H_halo"]],
            "T_plus_H": [float(x) for x in
                         (p["T_ballistic"] + p["H_halo"])],
            "ratio_up_over_one": [float(x) for x in p["up_factor"]],
            "sigma_m": o["sigma"]["sigma_m"],
            "sigma_px_at_1024": o["sigma_px"],
        }
    d06 = out["depths"]["z0.6"]; d40 = out["depths"]["z4"]
    out["up_leg_effect"] = {
        "green_up_factor_z0.6": d06["up"][1],
        "green_up_factor_z4": d40["up"][1],
        "green_full_over_down_only_z0.6": float(
            d06["two_way_beer"][1] / d06["down"][1]),
        "green_full_over_down_only_z4": float(
            d40["two_way_beer"][1] / d40["down"][1]),
        "statement": ("the return leg multiplies the ballistic beam by "
                      "exp(-c s_up); at z = 4 m that factor is %.3g in the "
                      "green channel, so omitting it would brighten the "
                      "bottom by %.1fx"
                      % (d40["up"][1], 1.0 / max(d40["up"][1], 1e-30))),
    }
    return out


def checks_all(g=G_DEFAULT):
    return {
        "hg": checks_hg(g),
        "sigma_limits": checks_sigma_limits(g),
        "sigma_monte_carlo": checks_sigma_mc(g=g),
        "partition": checks_partition(),
        "presets": checks_presets(),
        "pure_water_anchors": checks_legacy_vs_morel(),
        "beer_two_way_model": checks_beer_two_way(),
        "g_sensitivity": {
            "note": ("sigma scales as sqrt(1-g); the band is the Petzold "
                     "(1972) 0.85-0.95 range around the design's 0.85"),
            "rows": {("%.2f" % gg): sigma_from_scattering(
                1.0, 1.2, 1.25, gg)["sigma_m"] for gg in
                (0.80, 0.85, 0.90, 0.92, 0.95)},
        },
    }
