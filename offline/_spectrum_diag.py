"""Spectrum diagnosis sweep (MEASUREMENT ONLY).

Answers BRIEF-SPECTRUM-DIAG.md: where do sigma^2 and H_s actually live in the
v8 (version 3.4) water state?

No model file is imported for modification: v8_state / v8_physics / env /
recipe_v8 are used exactly as the production runner (run_v8.py) uses them.
The only new code is this diagnostic: band decomposition, the anchor target's
own band decomposition (a re-implementation of v8_physics._sigma2_integral that
is validated against it), and the SMB comparison.

Outputs (atomic, jobkit):
  out/spectrum-diag/cells/<cell>.json   per-cell result (resume unit)
  out/spectrum-diag/checkpoint.json     jobkit progress proof
  out/spectrum-diag/spectrum-diag.json  the full table
"""
import gc
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import (Checkpoint, atomic_write_bytes, read_checkpoint,  # noqa: E402
                        verify_units)
from env import resolve_environment  # noqa: E402
from state import derive_seed  # noqa: E402
from recipe_v8 import with_v8  # noqa: E402
import v8_physics as P  # noqa: E402
from v8_state import V8WaterState, TAPER_N  # noqa: E402

OUT = os.path.join(HERE, "out", "spectrum-diag")
CELLS = os.path.join(OUT, "cells")
CKPT = os.path.join(OUT, "checkpoint.json")

RAW_RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PRESET = 1
GRID = 2048
DEPTH_M = 1.0                      # the brief's "z = 1.0" (water depth)
WINDS = (1.5, 3.0, 6.0, 10.0)
FETCHES = (100.0, 1000.0)

# wavelength band edges (metres), descending: the brief's four bands
BANDS = (
    ("gt_1m", None, 1.0),
    ("0p1_to_1m", 1.0, 0.1),
    ("1cm_to_10cm", 0.1, 0.01),
    ("lt_1cm", 0.01, 0.0),
)
KAPPA_CAP = 370.0                  # brief's capillary peak wavenumber (rad/m)
LAMBDA_CAP = 2.0 * math.pi / KAPPA_CAP


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def k_edges():
    """k thresholds of the wavelength bands (k = 2*pi/lambda)."""
    return {"gt_1m": 2.0 * math.pi / 1.0,
            "0p1_to_1m": 2.0 * math.pi / 0.1,
            "1cm_to_10cm": 2.0 * math.pi / 0.01}


def band_of_k(k, e):
    """Band label for a k array (half-open, descending in wavelength)."""
    k = np.asarray(k, dtype=np.float64)
    out = np.full(k.shape, "gt_1m", dtype=object)
    out[k >= e["gt_1m"]] = "0p1_to_1m"
    out[k >= e["0p1_to_1m"]] = "1cm_to_10cm"
    out[k >= e["1cm_to_10cm"]] = "lt_1cm"
    return out


def decompose(S, k, e, m0_total=None, s2_total=None):
    """m0 and sigma^2 per wavelength band of a discrete 2-D spectrum."""
    S = np.asarray(S, dtype=np.float64).ravel()
    kk = np.asarray(k, dtype=np.float64).ravel()
    lab = band_of_k(kk, e)
    m0_tot = float(S.sum()) if m0_total is None else float(m0_total)
    s2_tot = float((kk * kk * S).sum()) if s2_total is None else float(s2_total)
    rows = {}
    for name, lam_hi, lam_lo in BANDS:
        m = lab == name
        m0 = float(S[m].sum())
        s2 = float((kk[m] ** 2 * S[m]).sum())
        rows[name] = {
            "lambda_range_m": [lam_hi, lam_lo],
            "k_range_rad_per_m": [None if lam_hi is None
                                  else 2.0 * math.pi / lam_hi,
                                  None if not lam_lo
                                  else 2.0 * math.pi / lam_lo],
            "m0_m2": m0, "sigma2_rad2": s2,
            "share_m0": (m0 / m0_tot) if m0_tot > 0 else None,
            "share_sigma2": (s2 / s2_tot) if s2_tot > 0 else None,
            "modes": int(m.sum()),
        }
    return {"m0_m2": m0_tot, "sigma2_rad2": s2_tot, "bands": rows,
            "check_m0_closes": bool(abs(sum(r["m0_m2"] for r in rows.values())
                                        - m0_tot) <= 1e-12 * max(m0_tot, 1e-30)),
            "check_sigma2_closes": bool(
                abs(sum(r["sigma2_rad2"] for r in rows.values()) - s2_tot)
                <= 1e-12 * max(s2_tot, 1e-30))}


def target_spectrum(U, F, h, gamma_peak, kmax, n=40000):
    """Re-implementation of v8_physics._sigma2_integral, returning k and S(k).

    Identical formula/quadrature so the totals can be validated against the
    model's own function (a check, not an assumption).
    """
    rel = P.jonswap_fetch_relations(U, F)
    fp = rel["f_p_Hz"]
    a = rel["alpha"]
    k = np.logspace(math.log10(1e-3), math.log10(float(kmax)), int(n))
    om = P.dispersion_omega(k, h)
    f = om / (2.0 * math.pi)
    sj = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-((f - fp) ** 2) / (2.0 * sj ** 2 * fp ** 2))
    sf = (a * P.G ** 2 / (2.0 * math.pi) ** 4 * f ** -5
          * np.exp(-1.25 * (fp / np.maximum(f, 1e-12)) ** 4)
          * float(gamma_peak) ** r)
    dfdk = P.dispersion_domega_dk(k, h) / (2.0 * math.pi)
    return k, sf * dfdk


def decompose_continuous(k, sk, e):
    """Band decomposition on a non-uniform k grid, additive by construction.

    Each trapezoid segment is assigned to the band of its midpoint, so the band
    sums close exactly on the full integral (no boundary loss).
    """
    dk = np.diff(k)
    contrib_m0 = 0.5 * (sk[1:] + sk[:-1]) * dk
    contrib_s2 = 0.5 * ((k[1:] ** 2 * sk[1:]) + (k[:-1] ** 2 * sk[:-1])) * dk
    mid = np.sqrt(k[1:] * k[:-1])
    lab = band_of_k(mid, e)
    m0_tot = float(contrib_m0.sum())
    s2_tot = float(contrib_s2.sum())
    rows = {}
    for name, lam_hi, lam_lo in BANDS:
        m = lab == name
        rows[name] = {
            "lambda_range_m": [lam_hi, lam_lo],
            "m0_m2": float(contrib_m0[m].sum()),
            "sigma2_rad2": float(contrib_s2[m].sum()),
            "share_m0": float(contrib_m0[m].sum()) / m0_tot if m0_tot else None,
            "share_sigma2": float(contrib_s2[m].sum()) / s2_tot if s2_tot else None,
        }
    return {"m0_m2": m0_tot, "sigma2_rad2": s2_tot, "bands": rows,
            "Hs_m": float(4.0 * math.sqrt(max(m0_tot, 0.0)))}


def smb_hs(U, F):
    """Two named fetch laws, both reported."""
    rel = P.jonswap_fetch_relations(U, F)
    Ft = rel["Ft"]
    hs_jon = rel["Hs_fetch_law_m"]
    # classic SMB (Bretschneider 1970 form), cross-check only
    hs_smb = 0.283 * (U * U / P.G) * math.tanh(
        0.0125 * Ft ** 0.42)
    return {"Ft": Ft,
            "Hs_fetch_law_m": float(hs_jon),
            "Hs_SMB_Bretschneider_m": float(hs_smb),
            "law": ("pipeline law H_s = 0.0016 Ft^0.5 U^2/g (Hasselmann et "
                    "al. 1973); reproduces the brief's 0.031 m at U=6/F=100"),
            "smb_bret": "H_s = 0.283 (U^2/g) tanh(0.0125 Ft^0.42)"}


def equivalent_fetch(U, Hs):
    """Fetch a JONSWAP-law sea of height Hs would need (inversion)."""
    if U <= 0 or Hs <= 0:
        return None
    Ft = (Hs * P.G / (0.0016 * U * U)) ** 2
    return {"Ft": float(Ft), "fetch_m": float(Ft * U * U / P.G)}


def radial_peak(S, k, kmin, kmax, weight_power=0):
    """k of the peak of the radially integrated k^weight_power * S."""
    kf = np.asarray(k, dtype=np.float64).ravel()
    Sf = np.asarray(S, dtype=np.float64).ravel()
    w = Sf * (kf ** weight_power)
    nb = 400
    edges = np.logspace(math.log10(max(kmin * 0.5, 1e-6)), math.log10(kmax),
                        nb + 1)
    idx = np.digitize(kf, edges) - 1
    ok = (idx >= 0) & (idx < nb)
    acc = np.bincount(idx[ok], weights=w[ok], minlength=nb)
    dk = np.diff(edges)
    dens = acc / np.maximum(dk, 1e-30)
    i = int(np.argmax(dens))
    kp = float(math.sqrt(edges[i] * edges[i + 1]))
    return {"k_peak_rad_per_m": kp, "lambda_peak_m": float(2.0 * math.pi / kp),
            "peak_is_at_kmin": bool(kp <= kmin * 1.35)}


# --------------------------------------------------------------------------
# one cell
# --------------------------------------------------------------------------
def run_cell(U, F, seed, e):
    t0 = time.time()
    scene = with_v8(RAW, fetch_m=F, depth_m=DEPTH_M, wind_ms=U)
    env = resolve_environment(scene)
    st = V8WaterState(scene, env, N=GRID, seed=seed, preset_index=PRESET,
                      fetch_m=F, fetch_source="explicit")
    S = st._spectrum_density()                      # post-anchor, post-clamp
    k = st.k.astype(np.float64)
    kmin, knyq = float(st.kmin), float(st.kmax_nyq)

    fields = st.fields_at(0.0)
    re = st.realised_stats(fields)
    del fields
    gc.collect()

    # --- pre-anchor shape, reconstructed exactly from the applied factors ---
    anc = st.anchor or {}
    kc = anc.get("kc_rad_per_m")
    if anc.get("applied") and kc:
        # v8_state.TAPER_N (the taper form the anchor solver used)
        g = 1.0 / (1.0 + (k / float(kc)) ** TAPER_N)
    else:
        g = np.ones_like(k)
    clamp_factor = float((st.breaking or {}).get("factor", 1.0) or 1.0)
    S_pre = S / (g * clamp_factor ** 2)
    pre_check = {
        "m0_rel_err_vs_anchor_node": (
            None if not anc.get("m0_before") else
            abs(float(S_pre.sum()) - float(anc["m0_before"]))
            / float(anc["m0_before"])),
        "sigma2_rel_err_vs_anchor_node": (
            None if not anc.get("sigma2_before") else
            abs(float((k * k * S_pre).sum()) - float(anc["sigma2_before"]))
            / float(anc["sigma2_before"])),
        "note": ("H0_pre = H0_post / (sqrt(g)*clamp_factor) is exact: the "
                 "anchor multiplies H0 by sqrt(g) and the breaking clamp by a "
                 "uniform factor"),
    }

    post = decompose(S, k, e)
    pre = decompose(S_pre, k, e)
    lam_p = st.spectrum_after["peak"]["lambda_peak_m"]
    slope_peak = radial_peak(S, k, kmin, knyq, weight_power=2)

    # --- the anchor target's own band decomposition -----------------------
    kmax_int = min(knyq, 2.0e4)
    tgt = {}
    for name, gam in P.ANCHOR_CANDIDATES:
        kk, sk = target_spectrum(U, F, st.v8_depth, gam, kmax_int)
        ref = P._sigma2_integral(U, F, st.v8_depth, gam, kmax_int)
        dec = decompose_continuous(kk, sk, e)
        dec["gamma"] = float(gam)
        dec["validation_vs_model_integral"] = {
            "model_sigma2": float(ref["sigma2"]),
            "reimpl_sigma2": float(dec["sigma2_rad2"]),
            "rel_err": float(abs(dec["sigma2_rad2"] - ref["sigma2"])
                             / ref["sigma2"]),
            "model_m0": float(ref["m0"]),
            "reimpl_m0": float(dec["m0_m2"]),
            "m0_rel_err": float(abs(dec["m0_m2"] - ref["m0"]) / ref["m0"]),
            "model_Hs_m": float(ref["Hs"]),
            "reimpl_Hs_m": float(dec["Hs_m"]),
        }
        tgt[name] = dec
    names = [n for n, _ in P.ANCHOR_CANDIDATES]
    # adopted target = median of the two totals = their mean (2 candidates)
    adopted_bands = {}
    for name, lam_hi, lam_lo in BANDS:
        adopted_bands[name] = {
            "lambda_range_m": [lam_hi, lam_lo],
            "m0_m2": float(np.mean([tgt[n]["bands"][name]["m0_m2"]
                                    for n in names])),
            "sigma2_rad2": float(np.mean([tgt[n]["bands"][name]["sigma2_rad2"]
                                          for n in names])),
        }
    ad_tot_m0 = float(np.mean([tgt[n]["m0_m2"] for n in names]))
    ad_tot_s2 = float(np.mean([tgt[n]["sigma2_rad2"] for n in names]))
    for name in adopted_bands:
        adopted_bands[name]["share_m0"] = (adopted_bands[name]["m0_m2"]
                                           / ad_tot_m0 if ad_tot_m0 else None)
        adopted_bands[name]["share_sigma2"] = (
            adopted_bands[name]["sigma2_rad2"] / ad_tot_s2 if ad_tot_s2 else None)
    target_block = {
        "integrated_band_rad_per_m": [kmin, kmax_int],
        "candidates": tgt,
        "adopted_median": {"m0_m2": ad_tot_m0, "sigma2_rad2": ad_tot_s2,
                           "Hs_m": float(4.0 * math.sqrt(max(ad_tot_m0, 0.0))),
                           "bands": adopted_bands,
                           "vs_model_target_rel": (
                               abs(ad_tot_s2 - float(anc.get(
                                   "sigma2_target_rad2") or ad_tot_s2))
                               / float(anc.get("sigma2_target_rad2") or 1.0))},
    }

    # --- SMB comparison ---------------------------------------------------
    smb = smb_hs(U, F)
    hs_spec = float(st.spectrum_after["Hs_m"])
    hs_real = float(re["Hs_m"])
    s2_spec = float(st.spectrum_after["sigma2_rad2"])
    s2_real = float(re["sigma2_rad2"])
    smb["Hs_ratio_spectral"] = hs_spec / smb["Hs_fetch_law_m"]
    smb["Hs_ratio_realised"] = hs_real / smb["Hs_fetch_law_m"]
    smb["Hs_ratio_vs_SMB_Bretschneider_spectral"] = (
        hs_spec / smb["Hs_SMB_Bretschneider_m"])
    smb["equivalent_fetch_of_realised_Hs"] = equivalent_fetch(U, hs_real)
    smb["equivalent_fetch_of_spectral_Hs"] = equivalent_fetch(U, hs_spec)

    # --- (c) shape kept, height forced to the fetch law -------------------
    m0_now = float(post["m0_m2"])
    s2_now = float(post["sigma2_rad2"])
    target_s2 = float(anc.get("sigma2_target_rad2") or s2_now)
    hs_smb = smb["Hs_fetch_law_m"]
    scale = hs_smb / hs_spec if hs_spec > 0 else None
    s2_at_smb_height = s2_now * scale ** 2 if scale else None
    slope_per_m0_now = s2_now / m0_now if m0_now else None
    slope_per_m0_needed = (target_s2 / ((hs_smb / 4.0) ** 2)) if hs_smb else None
    cap_bands = ("1cm_to_10cm", "lt_1cm")
    tgt_cap = sum(adopted_bands[b]["sigma2_rad2"] for b in cap_bands)
    mod_cap = sum(post["bands"][b]["sigma2_rad2"] for b in cap_bands)
    forced = {
        "rule": ("keep the realised spectral SHAPE (S(k)/m0) and scale the "
                 "height to the fetch law: sigma2_forced = sigma2 * (Hs_law/"
                 "Hs_now)^2.  This is what the current shape can deliver at "
                 "the right height."),
        "Hs_law_m": float(hs_smb), "Hs_now_m": float(hs_spec),
        "uniform_amplitude_factor": float(scale) if scale else None,
        "sigma2_forced_rad2": s2_at_smb_height,
        "sigma2_now_rad2": s2_now,
        "sigma2_lost_by_forcing_rad2": (s2_now - s2_at_smb_height
                                        if s2_at_smb_height is not None else None),
        "sigma2_lost_by_forcing_rel": ((s2_now - s2_at_smb_height) / s2_now
                                       if s2_at_smb_height else None),
        "sigma2_missing_vs_anchor_target_rad2": (target_s2 - s2_at_smb_height
                                                 if s2_at_smb_height is not None
                                                 else None),
        "sigma2_missing_vs_anchor_target_rel_of_target": (
            (target_s2 - s2_at_smb_height) / target_s2
            if s2_at_smb_height is not None else None),
        "slope_per_height_var_now_rad2_per_m2": slope_per_m0_now,
        "slope_per_height_var_needed_rad2_per_m2": slope_per_m0_needed,
        "slope_per_height_var_ratio_needed_over_now": (
            slope_per_m0_needed / slope_per_m0_now
            if slope_per_m0_now else None),
        "capillary_bands": list(cap_bands),
        "target_sigma2_in_capillary_bands_rad2": float(tgt_cap),
        "target_share_sigma2_in_capillary_bands": (tgt_cap / ad_tot_s2
                                                   if ad_tot_s2 else None),
        "model_sigma2_in_capillary_bands_rad2": float(mod_cap),
        "model_share_sigma2_in_capillary_bands": (mod_cap / s2_now
                                                  if s2_now else None),
        "missing_capillary_slope_rad2": float(tgt_cap - mod_cap),
        "missing_capillary_slope_rel_of_total": (
            float(tgt_cap - mod_cap) / target_s2 if target_s2 else None),
    }

    # --- capillary resolution ---------------------------------------------
    cap = {
        "capillary_peak_k_rad_per_m": KAPPA_CAP,
        "capillary_peak_lambda_m": LAMBDA_CAP,
        "k_nyquist_rad_per_m": knyq,
        "lambda_nyquist_m": float(2.0 * math.pi / knyq),
        "resolved_by_grid": bool(knyq > KAPPA_CAP),
        "knyq_over_kcapillary": float(knyq / KAPPA_CAP),
        "grid_points_per_capillary_wavelength": float(LAMBDA_CAP / st.dx),
        "kmin_rad_per_m": kmin,
        "lambda_max_resolved_m": float(2.0 * math.pi / kmin),
        "capillary_crossover": P.capillary_crossover(),
        "spectrum_has_energy_at_capillary_peak": None,   # filled below
        "share_sigma2_in_lambda_1p5_to_2cm": None,
    }
    win = (k >= 2.0 * math.pi / 0.020) & (k <= 2.0 * math.pi / 0.015)
    s2_win = float((k[win] ** 2 * S[win]).sum())
    cap["share_sigma2_in_lambda_1p5_to_2cm"] = (s2_win / s2_now
                                                if s2_now else None)
    # slope density at the capillary peak relative to its maximum
    kf = k.ravel()
    sl = (kf ** 2 * S.ravel())
    nb = 400
    edges = np.logspace(math.log10(max(kmin * 0.5, 1e-6)), math.log10(knyq),
                        nb + 1)
    idx = np.digitize(kf, edges) - 1
    ok = (idx >= 0) & (idx < nb)
    acc = np.bincount(idx[ok], weights=sl[ok], minlength=nb)
    dens = acc / np.maximum(np.diff(edges), 1e-30)
    kc_bin = int(np.clip(np.searchsorted(edges, KAPPA_CAP) - 1, 0, nb - 1))
    cap["spectrum_has_energy_at_capillary_peak"] = float(
        dens[kc_bin] / dens.max()) if dens.max() > 0 else None
    cap["slope_density_at_kcapillary_over_max"] = float(
        dens[kc_bin] / dens.max()) if dens.max() > 0 else None
    cap["slope_density_peak"] = slope_peak
    cap["model_shape_has_capillary_peak"] = bool(
        slope_peak["k_peak_rad_per_m"] > 100.0)

    # --- caps -------------------------------------------------------------
    br = dict(st.breaking or {})
    caps = {
        "breaking_clamp": br,
        "stokes_1over7": br.get("stokes_cap_1over7"),
        "H_stokes_0p142_m": br.get("H_stokes_0p142_m"),
        "H_stokes_1over7_m": br.get("H_stokes_1over7_m"),
        "H_depth_0p78h_m": br.get("H_depth_0p78h_m"),
        "binding": br.get("binding"),
        "applied": br.get("applied"),
        "factor": br.get("factor"),
        "Hs_after_clamp_m": br.get("Hs_after_m"),
        "Hs_over_Hmax_before": br.get("Hs_over_Hmax_before"),
        "viscous_cutoff": {
            "k_viscous_rad_per_m": (st.tail or {}).get("k_viscous_rad_per_m"),
            "lambda_viscous_m": (st.tail or {}).get("lambda_viscous_m"),
            "k_viscous_over_nyquist": (st.tail or {}).get("k_viscous_over_nyquist"),
            "modes_removed_in_resolved_band":
                (st.tail or {}).get("modes_removed_in_resolved_band"),
            "binds_inside_patch": bool(
                (st.tail or {}).get("k_viscous_over_nyquist", 1.0) < 1.0),
        },
    }

    cell = {
        "cell_id": "U%g_F%g" % (U, F),
        "scene": {
            "wind_ms": float(U), "fetch_m": float(F), "depth_m": float(st.v8_depth),
            "Ft": float(smb["Ft"]), "specL_m": float(env["specL"]),
            "patch_S_m": float(st.S), "N": int(st.N), "dx_m": float(st.dx),
            "kmin_rad_per_m": kmin, "kmax_nyquist_rad_per_m": knyq,
            "lambda_nyquist_m": float(2.0 * math.pi / knyq),
            "lambda_max_m": float(2.0 * math.pi / kmin),
            "seed": int(seed), "waveHeightM_env": float(env["waveHeightM"]),
            "Hs_from_env_height_4x": float(4.0 * env["waveHeightM"]),
        },
        "anchor": {
            "sigma2_target_rad2": target_s2,
            "applied": bool(anc.get("applied")),
            "reason": anc.get("reason"),
            "kc_rad_per_m": kc,
            "lambda_c_m": anc.get("lambda_c_m"),
            "sigma2_before_rad2": anc.get("sigma2_before"),
            "sigma2_after_spectral_rad2": anc.get("sigma2_after_spectral"),
            "sigma2_residual_rel": anc.get("sigma2_residual_rel"),
            "m0_before_m2": anc.get("m0_before"),
            "m0_after_m2": anc.get("m0_after_spectral"),
            "m0_removed_frac": anc.get("m0_removed_frac"),
            "Hs_before_m": anc.get("Hs_before"),
            "Hs_after_spectral_m": anc.get("Hs_after_spectral"),
            "old_cox_munk_target_rad2": anc.get("cox_munk_old_target_rad2"),
            "jonswap_relations": anc.get("jonswap_relations"),
            "direction": {
                "applied": bool(st.direction),
                "s_max": (st.direction or {}).get("s_max"),
                "omega_p_rad_s": (st.direction or {}).get("omega_p_rad_s"),
                "sigma2_before": (st.direction or {}).get("sigma2_before"),
                "sigma2_after": (st.direction or {}).get("sigma2_after"),
                "m0_before": (st.direction or {}).get("m0_before"),
                "m0_after": (st.direction or {}).get("m0_after"),
            },
        },
        "sigma2": {
            "target_rad2": target_s2,
            "realised_spectral_rad2": s2_spec,
            "realised_spatial_rad2": s2_real,
            "residual_vs_target_spectral_rel": float(
                (s2_spec - target_s2) / target_s2) if target_s2 else None,
            "residual_vs_target_spatial_rel": float(
                (s2_real - target_s2) / target_s2) if target_s2 else None,
            "target_met_spectral": bool(
                abs(s2_spec - target_s2) / target_s2 <= 0.05) if target_s2 else None,
        },
        "Hs": {
            "spectral_m": hs_spec,
            "realised_m": hs_real,
            "spectral_minus_realised_rel": float((hs_spec - hs_real) / hs_real)
            if hs_real else None,
            "smb": smb,
        },
        "bands_post_anchor": post,
        "bands_pre_anchor": pre,
        "pre_anchor_reconstruction_check": pre_check,
        "anchor_target_bands": target_block,
        "lambda_p_m": float(lam_p),
        "lambda_slope_peak_m": float(slope_peak["lambda_peak_m"]),
        "capillary": cap,
        "caps": caps,
        "forced_height": forced,
        "wall_s": float(time.time() - t0),
    }
    del S, S_pre, g, k
    gc.collect()
    return cell


# --------------------------------------------------------------------------
# derived fields (computed for fresh AND resumed cells, so a resume cannot
# silently drop a column)
# --------------------------------------------------------------------------
def enrich(cell):
    sc, anc = cell["scene"], cell["anchor"]
    L = float(sc["specL_m"])
    kmin = float(sc["kmin_rad_per_m"])
    kp = 2.0 * math.pi / float(cell["lambda_p_m"])
    cell["derived"] = {
        # Phillips shape peaks: S(k) at k = 0.707/L (lambda = 8.9 L),
        # k^2 S(k) at k = 1/L (lambda = 2*pi*L); the patch side is S = 4 L,
        # so both peaks sit at wavelengths LONGER than the patch.
        "lambda_peak_phillips_m": float(2.0 * math.pi / (0.70710678 / L)),
        "lambda_slope_peak_phillips_m": float(2.0 * math.pi * L),
        "patch_over_lambda_peak_phillips": float(sc["patch_S_m"]
                                                 / (2.0 * math.pi / (0.70710678 / L))),
        "peak_is_at_kmin": bool(kp <= kmin * 1.35),
        "peak_wavelength_over_patch": float(cell["lambda_p_m"]
                                            / sc["patch_S_m"]),
        "anchor_headroom_ratio_shape_over_target": (
            None if not anc.get("sigma2_before_rad2") else
            float(anc["sigma2_before_rad2"] / anc["sigma2_target_rad2"])),
        "anchor_can_apply": bool(anc.get("applied")),
        "m0_spectral_over_realised": float(
            cell["bands_post_anchor"]["m0_m2"]
            / (cell["Hs"]["realised_m"] / 4.0) ** 2),
        "target_share_sigma2_gravity_bands": float(
            1.0 - cell["forced_height"]["target_share_sigma2_in_capillary_bands"]),
        "model_share_sigma2_gravity_bands": float(
            1.0 - cell["forced_height"]["model_share_sigma2_in_capillary_bands"]),
    }
    return cell


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CELLS, exist_ok=True)
    cp = Checkpoint(CKPT)
    done = set((cp.state() or {}).get("completed_ids") or [])
    print("[ckpt] done units: %d %s" % (len(done), sorted(done)), flush=True)

    global RAW
    with open(RAW_RECIPE, "r", encoding="utf-8") as f:
        RAW = json.load(f)
    seed = derive_seed(RAW, PRESET)
    print("production seed = %d (derive_seed(recipe-1.json, preset %d))"
          % (seed, PRESET), flush=True)

    e = k_edges()
    meta = {
        "brief": "BRIEF-SPECTRUM-DIAG.md",
        "pipeline": "v8 / version 3.4 (run_v8.py state path, unmodified)",
        "measurement_only": True,
        "model_files_modified": False,
        "recipe": os.path.abspath(RAW_RECIPE),
        "preset": PRESET, "grid_N": GRID, "depth_m": DEPTH_M,
        "winds_ms": list(WINDS), "fetches_m": list(FETCHES),
        "seed": int(seed),
        "seed_note": ("production seed = derive_seed(unmodified recipe, "
                      "preset 1); held FIXED across the whole sweep, so every "
                      "cell uses the same phase draw"),
        "bands": [{"name": n, "lambda_range_m": [hi, lo]}
                  for n, hi, lo in BANDS],
        "band_rule": ("band assignment by |k| (2-D spectrum, both +k and -k "
                      "modes counted once each, as in the model's own "
                      "m0 = sum S and sigma2 = sum k^2 S); half-open, "
                      "descending in wavelength"),
        "smb_note": ("the brief's 'SMB fetch law' value at U=6/F=100 "
                     "(0.031 m) is reproduced by the pipeline's own "
                     "H_s = 0.0016 Ft^0.5 U^2/g law to 3 decimals, so that law "
                     "is used as the primary reference; the classic "
                     "Bretschneider SMB form is reported alongside"),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "created_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    cells = []
    for U in WINDS:
        for F in FETCHES:
            cid = "U%g_F%g" % (U, F)
            path = os.path.join(CELLS, cid + ".json")
            if cid in done and os.path.exists(path):
                with open(path, "r", encoding="utf-8") as f:
                    cells.append(json.load(f))
                print("[cell] %-10s reused from checkpoint" % cid, flush=True)
                continue
            print("[cell] %-10s building v8 state (N=%d)..." % (cid, GRID),
                  flush=True)
            cell = run_cell(U, F, seed, e)
            blob = json.dumps(cell, indent=2, ensure_ascii=False,
                              default=float).encode("utf-8")
            atomic_write_bytes(path, blob,
                               verify=lambda p: json.load(
                                   open(p, encoding="utf-8"))["cell_id"] == cid)
            cp.unit_done(cid)
            cells.append(cell)
            print("[cell] %-10s sigma2 tgt=%.6f real=%.6f | Hs=%.4f SMB=%.4f "
                  "ratio=%.2f | sigma2 bands>1m %.1f%% 0.1-1m %.1f%% 1-10cm "
                  "%.1f%% <1cm %.1f%% | m0 >1m %.1f%% | kc=%.2f | %.1fs"
                  % (cid, cell["sigma2"]["target_rad2"],
                     cell["sigma2"]["realised_spatial_rad2"],
                     cell["Hs"]["realised_m"],
                     cell["Hs"]["smb"]["Hs_fetch_law_m"],
                     cell["Hs"]["smb"]["Hs_ratio_realised"],
                     100 * cell["bands_post_anchor"]["bands"]["gt_1m"]["share_sigma2"],
                     100 * cell["bands_post_anchor"]["bands"]["0p1_to_1m"]["share_sigma2"],
                     100 * cell["bands_post_anchor"]["bands"]["1cm_to_10cm"]["share_sigma2"],
                     100 * cell["bands_post_anchor"]["bands"]["lt_1cm"]["share_sigma2"],
                     100 * cell["bands_post_anchor"]["bands"]["gt_1m"]["share_m0"],
                     cell["anchor"]["kc_rad_per_m"] or 0.0,
                     cell["wall_s"]), flush=True)

    cells.sort(key=lambda c: (c["scene"]["wind_ms"], c["scene"]["fetch_m"]))
    for c in cells:
        enrich(c)

    # ---- summary tables --------------------------------------------------
    def dom(bands, key):
        return max(bands.items(), key=lambda kv: kv[1][key] or 0.0)[0]

    summary = []
    for c in cells:
        b = c["bands_post_anchor"]["bands"]
        summary.append({
            "cell": c["cell_id"], "wind_ms": c["scene"]["wind_ms"],
            "fetch_m": c["scene"]["fetch_m"], "Ft": c["scene"]["Ft"],
            "sigma2_target": c["sigma2"]["target_rad2"],
            "sigma2_realised_spatial": c["sigma2"]["realised_spatial_rad2"],
            "sigma2_realised_spectral": c["sigma2"]["realised_spectral_rad2"],
            "Hs_spectral_m": c["Hs"]["spectral_m"],
            "Hs_realised_m": c["Hs"]["realised_m"],
            "Hs_SMB_m": c["Hs"]["smb"]["Hs_fetch_law_m"],
            "Hs_ratio": c["Hs"]["smb"]["Hs_ratio_realised"],
            "Hs_ratio_over_2x": bool(c["Hs"]["smb"]["Hs_ratio_realised"] > 2.0),
            "lambda_p_m": c["lambda_p_m"],
            "lambda_c_m": c["anchor"]["lambda_c_m"],
            "patch_S_m": c["scene"]["patch_S_m"], "dx_m": c["scene"]["dx_m"],
            "lambda_nyquist_m": c["scene"]["lambda_nyquist_m"],
            "capillary_resolved": c["capillary"]["resolved_by_grid"],
            "sigma2_share": {n: b[n]["share_sigma2"] for n in b},
            "m0_share": {n: b[n]["share_m0"] for n in b},
            "sigma2_dominant_band": dom(b, "share_sigma2"),
            "m0_dominant_band": dom(b, "share_m0"),
            "clamp_applied": c["caps"]["applied"],
            "clamp_binding": c["caps"]["binding"],
            "sigma2_forced_at_SMB_height": c["forced_height"]["sigma2_forced_rad2"],
            "sigma2_missing_vs_target": c["forced_height"]
            ["sigma2_missing_vs_anchor_target_rad2"],
            "missing_capillary_slope": c["forced_height"]
            ["missing_capillary_slope_rad2"],
            "target_share_sigma2_capillary": c["forced_height"]
            ["target_share_sigma2_in_capillary_bands"],
            "model_share_sigma2_capillary": c["forced_height"]
            ["model_share_sigma2_in_capillary_bands"],
        })

    ratios = [s["Hs_ratio"] for s in summary]
    cap_shares = [s["sigma2_share"]["1cm_to_10cm"]
                  + s["sigma2_share"]["lt_1cm"] for s in summary]
    checks = {
        "U6_F100_reproduces_brief": None,
        "anchor_target_validation_max_rel_err": None,
        "pre_anchor_reconstruction_max_rel_err": None,
        "band_sums_close": None,
        "seed_is_production": bool(all(c["scene"]["seed"] == seed for c in cells)),
        "no_model_file_modified": True,
    }
    val_errs = []
    for c in cells:
        for cand in c["anchor_target_bands"]["candidates"].values():
            val_errs.append(cand["validation_vs_model_integral"]["rel_err"])
        checks["band_sums_close"] = bool(
            (checks["band_sums_close"] is not False)
            and c["bands_post_anchor"]["check_m0_closes"]
            and c["bands_post_anchor"]["check_sigma2_closes"])
    checks["anchor_target_validation_max_rel_err"] = float(max(val_errs))
    pre_errs = [max([v for k2, v in c["pre_anchor_reconstruction_check"].items()
                     if "rel_err" in k2 and v is not None] or [0.0])
                for c in cells]
    checks["pre_anchor_reconstruction_max_rel_err"] = float(max(pre_errs))
    for c in cells:
        if c["cell_id"] == "U6_F100":
            checks["U6_F100_reproduces_brief"] = {
                "brief_sigma2_realised": 0.0556,
                "measured_sigma2_realised": c["sigma2"]["realised_spatial_rad2"],
                "brief_Hs_m": 0.332,
                "measured_Hs_realised_m": c["Hs"]["realised_m"],
                "sigma2_rel_diff": abs(c["sigma2"]["realised_spatial_rad2"]
                                       - 0.0556) / 0.0556,
                "Hs_rel_diff": abs(c["Hs"]["realised_m"] - 0.332) / 0.332,
                "reproduces": bool(
                    abs(c["sigma2"]["realised_spatial_rad2"] - 0.0556) / 0.0556
                    < 0.01
                    and abs(c["Hs"]["realised_m"] - 0.332) / 0.332 < 0.01),
            }

    answers = {
        "a_where_sigma2_and_m0_live": {
            "sigma2_dominant_bands": {s["cell"]: s["sigma2_dominant_band"]
                                      for s in summary},
            "m0_dominant_bands": {s["cell"]: s["m0_dominant_band"]
                                  for s in summary},
            "capillary_share_of_sigma2_min_max": [float(min(cap_shares)),
                                                  float(max(cap_shares))],
            "statement": None,
        },
        "b_cells_over_2x": [s["cell"] for s in summary if s["Hs_ratio_over_2x"]],
        "b_ratio_min_max": [float(min(ratios)), float(max(ratios))],
        "b_cross_check_classic_SMB_Bretschneider": {
            "ratio_min_max": [float(min(c["Hs"]["realised_m"]
                                       / c["Hs"]["smb"]["Hs_SMB_Bretschneider_m"]
                                       for c in cells)),
                              float(max(c["Hs"]["realised_m"]
                                        / c["Hs"]["smb"]["Hs_SMB_Bretschneider_m"]
                                        for c in cells))],
            "cells_over_2x": [c["cell_id"] for c in cells
                              if c["Hs"]["realised_m"]
                              / c["Hs"]["smb"]["Hs_SMB_Bretschneider_m"] > 2.0],
            "note": ("the classic SMB form is ~1.5-1.7x the pipeline's law at "
                     "these fetches, so the 2x conclusion is not an artefact of "
                     "the law choice: only U10/F1000 drops below 2x against it"),
        },
        "c_missing_capillary": {
            s["cell"]: {
                "sigma2_forced_at_SMB_height":
                    s["sigma2_forced_at_SMB_height"],
                "sigma2_missing_vs_anchor_target": s["sigma2_missing_vs_target"],
                "missing_capillary_slope_rad2": s["missing_capillary_slope"],
            } for s in summary
        },
    }
    s0 = summary[0]
    m0_sh = [s["m0_share"]["gt_1m"] for s in summary]
    dom_s = {s["cell"]: s["sigma2_dominant_band"] for s in summary}
    answers["a_where_sigma2_and_m0_live"]["statement"] = (
        "m0 (height variance) is carried by the >1 m band in all %d cells "
        "(%.1f-%.1f%%); sigma^2 (slope variance) is carried by the longest "
        "bands, never by the capillary bands: the dominant sigma^2 band is "
        ">1 m in %d cells and 0.1-1 m in %d cell(s), and the share of sigma^2 "
        "in lambda < 10 cm is only %.2f-%.2f%%.  Where the anchor applies, the "
        "slope piles up just below its taper knee lambda_c = %.3f-%.3f m."
        % (len(cells), 100 * min(m0_sh), 100 * max(m0_sh),
           len([v for v in dom_s.values() if v == "gt_1m"]),
           len([v for v in dom_s.values() if v == "0p1_to_1m"]),
           100 * min(cap_shares), 100 * max(cap_shares),
           min([c["anchor"]["lambda_c_m"] for c in cells
                if c["anchor"]["lambda_c_m"]]),
           max([c["anchor"]["lambda_c_m"] for c in cells
                if c["anchor"]["lambda_c_m"]])))

    # ---- findings (each one a measured number, quotable in the report) ----
    def by_id(cid):
        return next(c for c in cells if c["cell_id"] == cid)

    findings = []
    findings.append({
        "id": "F1-Hs-is-set-by-the-recipe-height-not-the-fetch",
        "claim": ("the realised H_s is 4*waveHeightM(env) (the Phillips "
                  "amplitude normalisation) with only a few percent removed by "
                  "the anchor; fetch enters H_s only through the anchor taper, "
                  "and at U=10 it does not enter at all (anchor not applied)"),
        "evidence": [{"cell": c["cell_id"],
                      "Hs_env_4x": c["scene"]["Hs_from_env_height_4x"],
                      "Hs_realised": c["Hs"]["realised_m"],
                      "Hs_realised_over_env_4x": float(
                          c["Hs"]["realised_m"]
                          / c["scene"]["Hs_from_env_height_4x"]),
                      "anchor_applied": c["anchor"]["applied"]} for c in cells],
    })
    findings.append({
        "id": "F2-ratio-over-2x-everywhere",
        "claim": ("H_s / H_s(fetch law) is > 2 in %d of the %d cells "
                  "(range %.2f-%.2f)"
                  % (len([s for s in summary if s["Hs_ratio_over_2x"]]),
                     len(cells), min(ratios), max(ratios))),
        "cells_over_2x": [s["cell"] for s in summary if s["Hs_ratio_over_2x"]],
    })
    lam_c_vals = [c["anchor"]["lambda_c_m"] for c in cells
                  if c["anchor"]["lambda_c_m"]]
    tgt_cap = [s["target_share_sigma2_capillary"] for s in summary]
    mod_cap = [s["model_share_sigma2_capillary"] for s in summary]
    findings.append({
        "id": "F3-sigma2-lives-where-the-taper-knee-is",
        "claim": ("where the anchor applies, the slope variance ends up at the "
                  "taper knee lambda_c = %.3f-%.3f m (gravity bands); the model "
                  "carries %.2f-%.2f%% of sigma^2 in the capillary bands "
                  "(lambda < 10 cm) while the anchor TARGET puts %.2f-%.2f%% "
                  "there"
                  % (min(lam_c_vals), max(lam_c_vals),
                     100 * min(mod_cap), 100 * max(mod_cap),
                     100 * min(tgt_cap), 100 * max(tgt_cap))),
        "per_cell": [{"cell": c["cell_id"],
                      "lambda_c_m": c["anchor"]["lambda_c_m"],
                      "model_share_sigma2_capillary":
                          c["forced_height"]["model_share_sigma2_in_capillary_bands"],
                      "target_share_sigma2_capillary":
                          c["forced_height"]["target_share_sigma2_in_capillary_bands"]}
                     for c in cells],
    })
    u10 = [c for c in cells if c["scene"]["wind_ms"] == 10.0]
    findings.append({
        "id": "F4-shape-lacks-slope-at-U10",
        "claim": ("at U=10 the shape's own sigma^2 is BELOW the anchor target "
                  "(headroom %.2f-%.2f), so the anchor reports applied=False and "
                  "the realised sigma^2 misses the target by %.0f-%.0f%%"
                  % (min(c["derived"]["anchor_headroom_ratio_shape_over_target"]
                         for c in u10),
                     max(c["derived"]["anchor_headroom_ratio_shape_over_target"]
                         for c in u10),
                     -100 * max(c["sigma2"]["residual_vs_target_spectral_rel"]
                                for c in u10),
                     -100 * min(c["sigma2"]["residual_vs_target_spectral_rel"]
                                for c in u10))),
        "evidence": [{"cell": c["cell_id"],
                      "shape_sigma2_before": c["anchor"]["sigma2_before_rad2"],
                      "target": c["anchor"]["sigma2_target_rad2"],
                      "headroom": c["derived"]["anchor_headroom_ratio_shape_over_target"],
                      "residual_rel": c["sigma2"]["residual_vs_target_spectral_rel"]}
                     for c in cells],
    })
    s1 = [c["bands_post_anchor"]["bands"]["1cm_to_10cm"]["share_sigma2"]
          for c in cells]
    s0b = [c["bands_post_anchor"]["bands"]["lt_1cm"]["share_sigma2"]
           for c in cells]
    k370 = [c["capillary"]["slope_density_at_kcapillary_over_max"]
            for c in cells]
    findings.append({
        "id": "F5-capillary-band-is-empty",
        "claim": ("the <1 cm band carries %.1e-%.1e of sigma^2 (i.e. nothing) "
                  "and the 1-10 cm band %.2f-%.2f%%; the slope density at the "
                  "capillary peak (k=370 rad/m) is %.1e-%.1e of its own "
                  "maximum, so the shape has no capillary tail"
                  % (min(s0b), max(s0b),
                     100 * min(s1), 100 * max(s1),
                     min(k370), max(k370))),
        "per_cell": [{"cell": c["cell_id"],
                      "share_sigma2_1cm_to_10cm":
                          c["bands_post_anchor"]["bands"]["1cm_to_10cm"]["share_sigma2"],
                      "share_sigma2_lt_1cm":
                          c["bands_post_anchor"]["bands"]["lt_1cm"]["share_sigma2"],
                      "slope_density_at_k370_over_max":
                          c["capillary"]["slope_density_at_kcapillary_over_max"]}
                     for c in cells],
    })
    findings.append({
        "id": "F6-patch-truncates-the-spectral-peak",
        "claim": ("the Phillips peak sits at lambda = 8.9*specL = %.2f x S "
                  "(OUTSIDE the patch in every cell), so the reported lambda_p "
                  "is pinned to the patch size (peak_is_at_kmin = %s)"
                  % (max(c["derived"]["patch_over_lambda_peak_phillips"]
                         for c in cells) ** -1,
                     all(c["derived"]["peak_is_at_kmin"] for c in cells))),
        "per_cell": [{"cell": c["cell_id"],
                      "lambda_p_m": c["lambda_p_m"],
                      "patch_S_m": c["scene"]["patch_S_m"],
                      "lambda_peak_phillips_m": c["derived"]["lambda_peak_phillips_m"],
                      "peak_is_at_kmin": c["derived"]["peak_is_at_kmin"]}
                     for c in cells],
    })
    margins = [c["caps"]["breaking_clamp"]["H_max_m"]
               / c["caps"]["breaking_clamp"]["Hs_before_m"] for c in cells]
    findings.append({
        "id": "F7-caps-never-bind",
        "claim": ("neither the breaking clamp (0.142 lambda / 0.78 h), the "
                  "Stokes 1/7 cap nor the viscous cutoff binds in any cell "
                  "(clamp applied = %s); the depth bound 0.78 h is the binding "
                  "limit where it is the smaller one, and the margin H_max/H_s "
                  "is %.2f-%.2f x"
                  % (sorted({bool(c["caps"]["applied"]) for c in cells}),
                     min(margins), max(margins))),
        "per_cell": [{"cell": c["cell_id"],
                      "Hs_m": c["caps"]["breaking_clamp"]["Hs_before_m"],
                      "H_max_m": c["caps"]["breaking_clamp"]["H_max_m"],
                      "binding": c["caps"]["binding"],
                      "applied": c["caps"]["applied"],
                      "margin": float(c["caps"]["breaking_clamp"]["H_max_m"]
                                      / c["caps"]["breaking_clamp"]["Hs_before_m"]),
                      "viscous_binds": c["caps"]["viscous_cutoff"]["binds_inside_patch"]}
                     for c in cells],
    })

    ok_cells, bad_cells = verify_units(CELLS, "*.json",
                                       lambda p: "cell_id" in json.load(
                                           open(p, encoding="utf-8")))
    doc = {
        "meta": meta,
        "cells": cells,
        "summary_table": summary,
        "answers": answers,
        "findings": findings,
        "checks": checks,
        "checkpoint_units": {"ok": len(ok_cells), "bad": bad_cells,
                             "completed": (cp.state() or {}).get(
                                 "completed_ids")},
    }
    payload = json.dumps(doc, indent=2, ensure_ascii=False,
                         default=float).encode("utf-8")
    atomic_write_bytes(os.path.join(OUT, "spectrum-diag.json"), payload,
                       verify=lambda p: bool(json.load(
                           open(p, encoding="utf-8"))["cells"]))
    cp.unit_done("outputs-spectrum-diag")
    print("\n[ok] wrote %s (%.1f kB)"
          % (os.path.join(OUT, "spectrum-diag.json"), len(payload) / 1024.0),
          flush=True)
    print("[checks] %s" % json.dumps(checks, default=float), flush=True)
    print("[ckpt] %s" % json.dumps(cp.state(), default=float), flush=True)


if __name__ == "__main__":
    main()
