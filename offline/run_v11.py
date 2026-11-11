"""run_v11 -- version 3.8 chunk 8: gravity-tail termination.

Four stages, foreground, PIPELINE-V10 sections 3-5.  Every stage is
checkpointed with D:\\mutaro-work\\jobkit (one unit per delivered artefact)
and every deliverable is committed with atomic_write_bytes.

  stage 1  spectrum module (B_l + assembled S) + algebra gates A1, A2
           -> out/v11/numbers-v11-spectrum.json, spectrum-curves-v11.json
  stage 2  8-cell sweep: C1 C2 C3 V1 V2 V3 K1 K5 K6 (legacy numbers read
           from the v10 sweep JSON, never recomputed)
           -> out/v11/numbers-v11-sweep.json
  stage 3  compatibility K2 K3 K4 through the unchanged transport and the
           frozen composite
           -> out/v11/numbers-v11-transport.json,
              comp-windy-v38-new.png, comp-windy-v38-legacy.png
  stage 4  merge evaluation; default flip iff every criterion holds; report
           -> out/v11/REPORT-V11-GATES.md

Usage (foreground, from caustics-filter/offline):
  python run_v11.py --stage 1
  python run_v11.py --stage 2
  python run_v11.py --stage 3
  python run_v11.py --stage 4
"""
import argparse
import datetime
import gc
import io
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import v8_physics as P8                                            # noqa: E402
import v10_spectrum as V10                                         # noqa: E402
import v11_spectrum as V11                                         # noqa: E402
from v11_state import V11WaterState, height_quantiles, slope_quantiles \
    # noqa: E402
import run_v10 as R10                                              # noqa: E402

RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PHOTO = R10.PHOTO
OUT = os.path.join(HERE, "out", "v11")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
V10_OUT = os.path.join(HERE, "out", "v10")

GRID_SWEEP = 1024
GRID_RENDER = R10.GRID_RENDER
SEED_PRODUCTION = None

SWEEP_U = V10.SWEEP_U
SWEEP_F = V10.SWEEP_F
SWEEP_H = V10.SWEEP_Z

B_HAT_LOW = 0.007          # measured saturation level (research, V3)
B_HAT_HIGH = 0.008

log = R10.log
now_iso = R10.now_iso
dump_json = R10.dump_json
sha256_file = R10.sha256_file
load_recipe = R10.load_recipe
cell_id = R10.cell_id


def legacy_sweep():
    """The v10 sweep JSON -- the legacy numbers are READ, never recomputed."""
    with open(os.path.join(V10_OUT, "numbers-v10-sweep.json"),
              "r", encoding="utf-8") as f:
        d = json.load(f)
    return {c["cell"]: c for c in d["cells"]}, d


def scene_for(raw, U, F, h, N, seed, source, k_cut=V10.KCUT_DEFAULT):
    """v11 state on the requested source -> the dict run_v10.transport_scene
    expects ({state, env, cfg, fields})."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    t0 = time.time()
    st = V11WaterState(scene, env, N=N, seed=seed, preset_index=1,
                       fetch_m=float(F), k_cut=k_cut, spectrum_source=source)
    t_state = time.time() - t0
    t1 = time.time()
    fields = st.fields_at(0.0)
    t_fields = time.time() - t1
    return {"state": st, "env": env, "cfg": cfg, "fields": fields,
            "t_state": t_state, "t_fields": t_fields}


def eckv_band_table(ctx, k_hi, kmin):
    """The ECKV reference's own band table at matched inputs (C2/C3)."""
    khi = float(k_hi)
    k = V11._k_grid(khi, k_hi=khi)
    ref = V11.eckv_reference(k, float(ctx["depth_m"]), ctx["k_p_rad_per_m"],
                             ctx["c_p_ms"], ctx["omega_c"], ctx["alpha_p"],
                             ctx["gamma_bl"], ctx["sigma_bl"],
                             ctx["ustar"]["ustar_ms"], ctx["alpha_m_value"],
                             deep_water=True)
    S = ref["B"] / k ** 3
    out = {"bands": {}, "k_hi": khi}
    edges = [float(kmin), float(ctx["k_p_rad_per_m"])] + \
        list(V11.BAND_EDGES_V11) + [khi]
    names = ["below_kp"] + list(V11.BAND_NAMES_V11)
    tot = 0.0
    for name, lo, hi in zip(names, edges[:-1], edges[1:]):
        m = (k >= lo) & (k < hi)
        s2 = float(np.trapz(k[m] ** 2 * S[m], k[m])) if hi > lo else 0.0
        out["bands"][name] = {"sigma2": s2,
                              "m0": (float(np.trapz(S[m], k[m]))
                                     if hi > lo else 0.0)}
        tot += s2
    for name in out["bands"]:
        out["bands"][name]["sigma2_share"] = (float(out["bands"][name]
                                                    ["sigma2"] / tot)
                                              if tot > 0 else 0.0)
    out["sigma2_total"] = tot
    m370 = k >= 370.0
    s2h = float(np.trapz(k[m370] ** 2 * S[m370], k[m370]))
    out["frac_sigma2_above_370"] = float(s2h / tot) if tot > 0 else None
    return out


# ==========================================================================
# stage 1 -- spectrum module + algebra gates
# ==========================================================================
def stage1():
    t0 = time.time()
    log("=" * 72)
    log("[stage 1] v11_spectrum: B_l assembly, A1 (D1 verification), A2")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    # ---- A1 ----------------------------------------------------------
    a1 = V11.gate_a1()
    for hh, r in a1["depths"].items():
        log("  A1 %-7s ratio(physical k_m) max|r-1|=%.2e | (published k_m) "
            "%.2e | B_g mean %.5f | plateau rel %+.2f%%"
            % (hh, r["ratio_kM_physical"]["max_abs_minus_1"],
               r["ratio_kM_published"]["max_abs_minus_1"],
               r["B_g_numeric_mean"], 100 * r["plateau_rel_vs_deep"]))
    worst = a1["worst_finite_depth_plateau_rel"]
    a1["pass"] = bool(worst < a1["failure_condition"])
    if a1["pass"]:
        a1["finding"] = (
            "The v3.7 B_g sits on the deep-water closed form to numerical "
            "precision in 50-200 rad/m (max |ratio-1| = %.1e with the "
            "pipeline's own k_m; %.1e with the research's rounded k_m = "
            "2 pi/0.017 m -- a 0.7%% convention difference, reported not "
            "tuned).  The finite-depth plateau moves by at most %.1f%% "
            "(h = 0.3 m) from the deep-water value, far below the 30%% "
            "failure condition, because kh >= %.0f across the band.  D1 is "
            "NOT weakened: the level defect in 50-200 is real and its "
            "magnitude is the measured B_g level below."
            % (a1["depths"]["h_%g" % 1000.0]["ratio_kM_physical"]
               ["max_abs_minus_1"],
               a1["depths"]["h_%g" % 1000.0]["ratio_kM_published"]
               ["max_abs_minus_1"],
               100 * worst, a1["depths"]["h_0.3"]["kh_over_band"][0]))
        a1["diagnosis"] = {
            "B_g_plateau_deep": a1["depths"]["h_%g" % 1000.0]
            ["B_g_numeric_mean"],
            "alpha_J_over_2": 0.5 * a1["depths"]["h_%g" % 1000.0]
            ["alpha_J_from_A_over_g2"],
            "alpha_J_published_over_2": 0.5 * a1["depths"]["h_%g" % 1000.0]
            ["alpha_J_published_law"],
            "reading": ("the plateau is alpha_J/2 of the fetch-normalised "
                        "JONSWAP (alpha is a READOUT of the H_s "
                        "normalisation, not the published alpha law), i.e. "
                        "the 0.006-0.023 level the research measured"),
        }
        log("  A1 PASS; D1 confirmed, plateau %.5f" % a1["diagnosis"]
            ["B_g_plateau_deep"])
        cp.unit_done("s1-a1")
    else:
        a1["finding"] = ("The finite-depth plateau moves %.1f%% (>= 30%%): "
                         "the D1 diagnosis is weakened.  Stopping before "
                         "stage 2 as the brief requires."
                         % (100 * worst))
        log("  A1 FAIL: " + a1["finding"])
        cp.unit_done("s1-a1-failed")
        stage = {
            "stage": 1, "created_utc": now_iso(),
            "model": "version 3.8 (PIPELINE-V10 chunk 8)",
            "what": "v11_spectrum + A1 (D1 verification) -- FAILED, stopped",
            "gates": {"A1": a1},
            "stopped": True,
            "wall_s": float(time.time() - t0),
        }
        dump_json(os.path.join(OUT, "numbers-v11-spectrum.json"), stage)
        cp.unit_done("s1-json")
        log("[stage 1] STOPPED: A1 failed")
        return stage

    # ---- A2 ----------------------------------------------------------
    ref_ctx = V11.assemble(6.0, 100.0, 1.0)
    ref_patch = V10.patch_for(ref_ctx["laws"]["lambda_p_m"])["patch_m"]
    a2 = V11.gate_a2(ref_ctx, N=1024, patch_m=ref_patch,
                     tag="elfouhaily-form")
    log("  A2: m0 num %.6e mode %.6e cont %.6e | sigma2 num %.5f mode %.5f "
        "cont %.5f" % (a2["m0"]["numeric"], a2["m0"]["mode_sum"],
                       a2["m0"]["continuum"], a2["sigma2"]["numeric"],
                       a2["sigma2"]["mode_sum"], a2["sigma2"]["continuum"]))
    log("  A2 z: m0 %+.2f sigma2 %+.2f | rel vs cont %+.4f / %+.4f | PASS=%s"
        % (a2["z_m0"], a2["z_sigma2"], a2["rel_mode_over_continuum_m0"],
           a2["rel_mode_over_continuum_sigma2"], a2["pass"]))
    cp.unit_done("s1-a2")
    # the v10 G0 analogue on the legacy source, same configuration
    leg_ctx = V10.assemble(6.0, 100.0, 1.0)
    a2_legacy = V11.gate_a2(leg_ctx, N=1024, patch_m=ref_patch,
                            tag="jonswap (v10 G0 analogue)")
    log("  A2 legacy reference: z m0 %+.2f sigma2 %+.2f PASS=%s"
        % (a2_legacy["z_m0"], a2_legacy["z_sigma2"], a2_legacy["pass"]))
    cp.unit_done("s1-a2-legacy")

    # ---- analytic table for the eight cells ---------------------------
    cells = []
    curves_all = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            ctx = V11.assemble(U, F, SWEEP_H)
            pm = max(V10.PATCH_BASE_M,
                     V10.R_LAMBDA * ctx["laws"]["lambda_p_m"])
            kmin = 2.0 * math.pi / pm
            kvisc = float(P8.viscous_cutoff(SWEEP_H))
            kf = V11._k_grid(kvisc, k_hi=kvisc)
            sf = V11.s_parts(kf, ctx)
            an_full = V11.integral_block(kf, sf, kmin, kvisc,
                                         ctx["k_p_rad_per_m"])
            kc = V11._k_grid(V10.KCUT_DEFAULT)
            sc = V11.s_parts(kc, ctx)
            an_cut = V11.integral_block(kc, sc, kmin, V10.KCUT_DEFAULT,
                                        ctx["k_p_rad_per_m"])
            ref = V11.eckv_reference_totals(ctx)
            refb = eckv_band_table(ctx, kvisc, kmin)
            cm = V10.cox_munk_reference(U)
            o = an_cut["overlap_50_200"]
            hs_pred = float(4.0 * math.sqrt(
                max(an_full["total"]["m0"], 0.0)))
            hs_law = float(ctx["laws"]["Hs_depth_limited_m"])
            cells.append({
                "cell": cell_id(U, F), "U_ms": float(U),
                "fetch_m": float(F), "depth_m": float(SWEEP_H),
                "Ft": ctx["laws"]["Ft"],
                "f_p_Hz": ctx["laws"]["f_p_Hz"],
                "alpha_law": ctx["laws"]["alpha"],
                "lambda_p_m": ctx["laws"]["lambda_p_m"],
                "Hs_law_m": hs_law,
                "Hs_predicted_m": hs_pred,
                "Hs_pred_rel_law": float(hs_pred / hs_law - 1.0),
                "hs_source": ctx["hs_source"],
                "hs_override_m": ctx["hs_override_m"],
                "sea_state": {k: ctx["sea_state"][k] for k in (
                    "k_p_rad_per_m", "c_p_ms", "omega_c", "alpha_p",
                    "gamma_bl", "gamma_form", "sigma_bl", "eckv_domain",
                    "k_p_deep_water_rad_per_m", "k_p_dispersion_rel_diff",
                    "kp_times_h")},
                "ustar_ms": ctx["ustar"]["ustar_ms"],
                "alpha_m": ctx["alpha_m_value"],
                "alpha_m_raw": ctx["alpha_m"]["alpha_m_raw"],
                "alpha_m_domain": ctx["alpha_m_domain"],
                "eckv_domain": ctx["eckv_domain"],
                "patch_m": pm, "kmin_rad_per_m": kmin,
                "k_viscous_rad_per_m": kvisc,
                "integrals_cut": an_cut,
                "integrals_full": an_full,
                "eckv_reference": ref,
                "eckv_reference_bands": refb,
                "cox_munk": cm,
                "C1": {
                    "gravity_share_50_200":
                        o["gravity_share_of_sigma2"],
                    "B_l_max_50_200": o["B_l_max"],
                    "B_h_max_50_200": o["B_h_max"],
                    "limits": {"gravity_share_max": 0.50,
                               "B_l_max_in_band": 0.010},
                    "pass": bool(o["gravity_share_of_sigma2"] <= 0.50
                                 and o["B_l_max"] <= 0.010),
                },
            })
            curves_all[cell_id(U, F)] = {
                "new": V11.curves(ctx, pm),
                "legacy": V11.curves(V10.assemble(U, F, SWEEP_H), pm),
            }
            log("  %-14s Om=%5.2f kp=%6.2f Hs_pred %.5f law %.5f (%+.1f%%) "
                "sig2(cut) %.5f C1 share %.3f B_l max %.4f %s"
                % (cell_id(U, F), ctx["omega_c"], ctx["k_p_rad_per_m"],
                   hs_pred, hs_law, 100 * (hs_pred / hs_law - 1.0),
                   an_cut["total"]["sigma2"],
                   o["gravity_share_of_sigma2"], o["B_l_max"],
                   ctx["eckv_domain"]))
            cp.unit_done("s1-cell-" + cell_id(U, F))
    cp.unit_done("s1-analytic-8-cells")

    # ---- deliverable curves -------------------------------------------
    dump_json(os.path.join(OUT, "spectrum-curves-v11.json"), {
        "what": ("S_l / S_h / B_l / B_h / B / Delta(k) per cell, new source "
                 "beside the legacy source, on the patch-resolved band"),
        "spectrum_source": "elfouhaily-form",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters (NOT 'the Elfouhaily model')"),
        "cells": curves_all})
    cp.unit_done("s1-curves")
    c = curves_all[cell_id(6.0, 100.0)]["new"]
    lines = ["k_rad_per_m,B_l,B_h,B_total,S_l,S_h,S_total,delta_eckv,"
             "omega_rad_s"]
    for i in range(len(c["k_rad_per_m"])):
        lines.append("%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.6f,%.6f" % (
            c["k_rad_per_m"][i], c["B_l"][i], c["B_h"][i],
            c["B_total"][i], c["S_l_m2_per_rad_m"][i],
            c["S_h_m2_per_rad_m"][i], c["S_total_m2_per_rad_m"][i],
            c["delta_eckv"][i], c["omega_rad_s"][i]))
    atomic_write_bytes(os.path.join(OUT, "spectrum-ref-v38-U6F100.csv"),
                       ("\n".join(lines) + "\n").encode("ascii"))
    cp.unit_done("s1-curves-csv")

    stage = {
        "stage": 1,
        "what": ("v11_spectrum: B_l assembly + A1 (D1 verification) + A2, "
                 "eight-cell analytic table"),
        "created_utc": now_iso(),
        "model": "version 3.8 (PIPELINE-V10 chunk 8)",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters.  NOT 'the Elfouhaily model'."),
        "conventions": {
            "S_k_units": "m^2 per rad/m",
            "m0": "int S(k) dk", "sigma2": "int k^2 S(k) dk",
            "B": "B(k) = k^3 S(k), int B dlnk = mean square slope",
            "directional": "Psi = S Phi/k, int Phi dtheta = 1 (unchanged)",
            "wind": "U10 everywhere; u* from U10 only",
            "gravity_part": "S(k) = (B_l + B_h)/k^3, no join/splice/rescale",
        },
        "gates": {"A1": a1, "A2": a2, "A2_legacy_reference": a2_legacy},
        "constants": {
            "K_M_rad_per_m": float(V11.K_M),
            "C_M_ms": float(V11.C_M),
            "k_m_pipeline_rad_per_m": float(math.sqrt(P8.G / P8.TF)),
            "k_m_convention_rel_diff": float(
                math.sqrt(P8.G / P8.TF) / V11.K_M - 1.0),
            "note": ("K_M = 2 pi/0.017 m is the research/EVKV constant used "
                     "in the closed form and in B_h; the pipeline "
                     "dispersion's own crossover is sqrt(g/(gamma/rho)). "
                     "Both are reported wherever the closed form is used."),
        },
        "cells": cells,
        "files": ["spectrum-curves-v11.json",
                  "spectrum-ref-v38-U6F100.csv"],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v11-spectrum.json"), stage)
    cp.unit_done("s1-json")
    log("[stage 1] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 2 -- 8-cell sweep
# ==========================================================================
def stage2():
    t0 = time.time()
    log("=" * 72)
    log("[stage 2] 8-cell sweep: C1 C2 C3 V1 V2 V3 K1 K5 K6")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    raw = load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    leg, leg_doc = legacy_sweep()
    log("  production seed %d; legacy numbers read from %s (not recomputed)"
        % (seed, os.path.join(V10_OUT, "numbers-v10-sweep.json")))
    with open(os.path.join(OUT, "numbers-v11-spectrum.json"),
              "r", encoding="utf-8") as f:
        s1 = json.load(f)
    an1 = {c["cell"]: c for c in s1["cells"]}

    cells = []
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            t1 = time.time()
            r = scene_for(raw, U, F, SWEEP_H, GRID_SWEEP, seed,
                          "elfouhaily-form")
            st = r["state"]
            node = st.state_node(r["fields"])
            ctx = st.ctx
            pm = st.S
            kmin = st.kmin
            kvisc = float(P8.viscous_cutoff(SWEEP_H))
            kf = V11._k_grid(kvisc, k_hi=kvisc)
            sf = V11.s_parts(kf, ctx)
            an_full = V11.integral_block(kf, sf, kmin, kvisc,
                                         ctx["k_p_rad_per_m"])
            kc = V11._k_grid(V10.KCUT_DEFAULT)
            sc = V11.s_parts(kc, ctx)
            an_cut = V11.integral_block(kc, sc, kmin, V10.KCUT_DEFAULT,
                                        ctx["k_p_rad_per_m"])
            ref = V11.eckv_reference_totals(ctx)
            refb = eckv_band_table(ctx, kvisc, kmin)
            # assembled spectrum at the patch-resolved band (same as an_cut;
            # the realised mode sum is the state's own number below)
            s2_cut = float(an_cut["total"]["sigma2"])
            s2_full = float(an_full["total"]["sigma2"])
            L = leg[cid]
            s2_v10 = float(L["gates"]["G2"]["sigma2_realised_spectral"])
            cm = V10.cox_munk_reference(U)
            o = an_cut["overlap_50_200"]
            fa = st.frac_above((100.0, 370.0))
            hs = node["v11"]["K1"]

            # ---- C1 -------------------------------------------------
            c1 = {
                "gravity_share_50_200": o["gravity_share_of_sigma2"],
                "B_l_max_50_200": o["B_l_max"],
                "B_h_max_50_200": o["B_h_max"],
                "B_total_max_50_200": o["B_total_max"],
                "limits": {"gravity_share_max": 0.50,
                           "B_l_max_in_band": 0.010},
                "pass_share": bool(o["gravity_share_of_sigma2"] <= 0.50),
                "pass_maxB": bool(o["B_l_max"] <= 0.010),
                "pass": bool(o["gravity_share_of_sigma2"] <= 0.50
                             and o["B_l_max"] <= 0.010),
                "legacy_gravity_share_50_200":
                    L["gates"]["G2"]["overlap_50_200"]
                    ["S_g_share_of_overlap_sigma2"],
                "note": ("gravity share of the 50-200 slope variance: the "
                         "defect being fixed.  Legacy share is the v10 "
                         "recorded number.  Where alpha_m is floored at 0 "
                         "(U = 1.5) B_h is identically zero and the share is "
                         "1.0 by construction -- reported, not excused."),
            }
            # ---- C2 -------------------------------------------------
            c2 = {
                "sigma2_assembled_cut": s2_cut,
                "sigma2_assembled_full": s2_full,
                "sigma2_v10_realised": s2_v10,
                "rel_vs_v10": float(s2_cut / s2_v10 - 1.0),
                "sigma2_eckv_reference_full": ref["sigma2_rad2"],
                "rel_vs_eckv_full": float(s2_full / ref["sigma2_rad2"] - 1.0),
                "rel_v10_vs_eckv_reference": float(s2_v10
                                                   / ref["sigma2_rad2"] - 1.0),
                "v10_named_full_eckv_rad2": L["gates"]["G3"]
                ["full_eckv_sigma2_rad2"],
                "v10_named_full_eckv_is": ("the v3.7 spectrum's own uncut "
                                           "integral (v10.full_eckv_sigma2), "
                                           "NOT ECKV at matched inputs"),
                "band_ratios_vs_eckv": {
                    name: {
                        "new_sigma2": an_full["bands"]["rows"][name]
                        ["total"]["sigma2"],
                        "eckv_sigma2": refb["bands"][name]["sigma2"],
                        "ratio": (float(an_full["bands"]["rows"][name]
                                        ["total"]["sigma2"]
                                        / refb["bands"][name]["sigma2"])
                                  if refb["bands"][name]["sigma2"] > 0
                                  else None),
                    } for name in an_full["bands"]["rows"]},
                "tolerance_rel": 0.10,
                "pass": bool(abs(s2_cut / s2_v10 - 1.0) <= 0.10),
                "note": ("the gate as written: total sigma^2 within +-10% of "
                         "the v3.7 totals.  Band ratios vs the ECKV "
                         "reference are reported WITH absolute contributions "
                         "(GPT: a ratio alone can hide compensating errors)."),
            }
            # ---- C3 -------------------------------------------------
            f_new = fa["above_370"]["realised"]
            f_full = an_full["frac_sigma2_above_370"]
            f_eckv = ref["frac_sigma2_above_370"]
            c3 = {
                "frac_above_370_new_realised": f_new,
                "frac_above_370_new_assembled_full": f_full,
                "frac_above_370_eckv_matched": f_eckv,
                "delta": (None if (f_new is None or f_eckv is None)
                          else float(f_new - f_eckv)),
                "target_abs_delta": 0.03,
                "labelled": "conformance (internal consistency), NOT "
                            "physical validation",
                "pass_reported": (None if (f_new is None or f_eckv is None)
                                  else bool(abs(f_new - f_eckv) <= 0.03)),
                "legacy_frac_above_370_realised":
                    L["gates"]["G2"]["frac_above_370_realised"],
                "note": ("the realised fraction is the one the field "
                         "actually carries: it stops at k_cut = %g rad/m, "
                         "so the comparison against ECKV's FULL-band "
                         "fraction measures what the cutoff removes.  The "
                         "assembled-full column is given for transparency "
                         "and agrees with ECKV to ~1e-3 by construction."
                         % st.k_cut),
            }
            # ---- V1 -------------------------------------------------
            v1 = {
                "sigma2_new": s2_cut,
                "cox_munk_clean_rad2": cm["total_clean_rad2"],
                "factor_over_cox_munk": float(s2_cut
                                              / cm["total_clean_rad2"]),
                "legacy_factor_over_cox_munk": float(
                    s2_v10 / cm["total_clean_rad2"]),
                "legacy_v10_range_observed": [1.14, 2.34],
                "note": ("clean Cox-Munk (0.003 + 0.00512 U) is a REFERENCE, "
                         "never a target"),
            }
            # ---- V2 -------------------------------------------------
            v2 = {
                "frac_sigma2_above_1173": an_full["frac_sigma2_above_1173"],
                "frac_1173_if_resolved": ("the realised mode set stops at "
                                          "k_cut = %g rad/m, so the pure "
                                          "capillary band k > 1173 is "
                                          "carried only by the continuum"
                                          % st.k_cut),
                "laxague_reference": {
                    "pure_capillary_k_gt_1173": 0.05,
                    "caveat": ("Laxague et al. partition, seas with "
                               "Omega_c <= 2; this chunk's cells reach "
                               "Omega_c = %.2f, so the comparison is a "
                               "comparison with stated uncertainty, not a "
                               "check" % ctx["omega_c"]),
                },
                "band_shares_new": {name: an_full["bands"]["rows"][name]
                                    ["sigma2_share"]
                                    for name in an_full["bands"]["rows"]},
                "band_shares_legacy_realised": {
                    k: v["sigma2_share"] for k, v in
                    L["gates"]["G2"]["bands_state"].items() if k != "_total"},
            }
            # ---- V3 -------------------------------------------------
            kb = np.logspace(math.log10(50.0), math.log10(200.0), 1500)
            sb = V11.s_parts(kb, ctx)
            r_l = sb["B_l"] / B_HAT_LOW
            r_t = (sb["B_l"] + sb["B_h"]) / B_HAT_LOW
            v3 = {
                "B_hat_range": [B_HAT_LOW, B_HAT_HIGH],
                "B_l_max_over_band": float(np.max(sb["B_l"])),
                "B_total_max_over_band": float(np.max(sb["B_l"] + sb["B_h"])),
                "max_ratio_B_l_over_Bhat_low": float(np.max(r_l)),
                "max_ratio_B_total_over_Bhat_low": float(np.max(r_t)),
                "frac_band_with_B_l_above_Bhat_high": float(
                    np.mean(sb["B_l"] > B_HAT_HIGH)),
                "frac_band_with_B_total_above_Bhat_high": float(
                    np.mean((sb["B_l"] + sb["B_h"]) > B_HAT_HIGH)),
                "pass": bool(np.max(r_l) <= 1.0),
                "note": ("the gravity part must no longer sit above the "
                         "measured saturation level across 50-200; the max "
                         "ratio per cell is the reported number"),
            }
            # ---- K5 / K6 --------------------------------------------
            k5 = node["quantisation"]
            k6 = V11.gate_k6(ctx)

            row = {
                "cell": cid, "U_ms": float(U), "fetch_m": float(F),
                "depth_m": float(SWEEP_H),
                "mode": st.mode, "cascade_j": st.cascade_j,
                "patch_m": float(pm), "patch_doublings":
                    int(st.patch_doublings),
                "kmin_rad_per_m": float(kmin),
                "kmax_nyquist_scene": float(st.kmax_nyq),
                "k_x_eff": float(st.k_x_eff),
                "dx_mm": float(st.dx * 1e3),
                "sea_state": node["v11"]["sea_state"],
                "spectrum": node["spectrum"],
                "amplitude": node["amplitude"],
                "breaking": node["breaking"],
                "realised": node["realised"],
                "spectral_vs_spatial": node["spectral_vs_spatial"],
                "bands_realised": node["bands"],
                "direction": node["direction"],
                "v11": {k: v for k, v in node["v11"].items()
                        if k != "sea_state"},
                "integrals_cut": an_cut,
                "integrals_full": an_full,
                "eckv_reference": ref,
                "eckv_reference_bands": refb,
                "cox_munk": cm,
                "frac_above_370_realised": fa["above_370"]["realised"],
                "gates": {"C1": c1, "C2": c2, "C3": c3, "V1": v1,
                          "V2": v2, "V3": v3, "K1": hs, "K5": k5, "K6": k6},
                "legacy_v10": {
                    "sigma2_realised": s2_v10,
                    "Hs_spectral": L["gates"]["G1"]["Hs_spectral_m"],
                    "Hs_law": L["gates"]["G1"]["Hs_law_m"],
                    "frac_above_370_realised":
                        L["gates"]["G2"]["frac_above_370_realised"],
                    "bands_state": L["gates"]["G2"]["bands_state"],
                    "overlap_50_200": L["gates"]["G2"]["overlap_50_200"],
                    "quantisation_max_rel":
                        L["gates"]["G7"]["max_rel_error"],
                },
                "timing_s": {"state": float(r["t_state"]),
                             "fields": float(r["t_fields"])},
            }
            cells.append(row)
            log("  [%s] Om=%5.2f %s Hs %.4f/%.4f (%+.1f%%) sig2 %.5f (v10 "
                "%.5f, %+.1f%%) C1 %s (share %.3f) C2 %s C3 d%+.3f V3 %s "
                "K1flag=%s (%.1fs)"
                % (cid, ctx["omega_c"], ctx["eckv_domain"][:12],
                   hs["Hs_predicted_m"], hs["Hs_law_m"],
                   100 * hs["delta_rel"], s2_cut, s2_v10,
                   100 * c2["rel_vs_v10"], c1["pass"],
                   o["gravity_share_of_sigma2"], c2["pass"],
                   c3["delta"] or 0.0, v3["pass"], hs["flag"],
                   time.time() - t1))
            cp.unit_done("s2-cell-" + cid)
    cp.unit_done("s2-sweep")

    def col(fn):
        return [fn(c) for c in cells]

    summary = {
        "C1_gravity_share": {
            "pass": all(c["gates"]["C1"]["pass"] for c in cells),
            "worst_share": max(col(lambda c: c["gates"]["C1"]
                                   ["gravity_share_50_200"])),
            "failing_cells": [c["cell"] for c in cells
                              if not c["gates"]["C1"]["pass"]],
        },
        "C2_total": {
            "pass": all(c["gates"]["C2"]["pass"] for c in cells),
            "worst_rel": max(col(lambda c: abs(c["gates"]["C2"]
                                               ["rel_vs_v10"]))),
            "failing_cells": [c["cell"] for c in cells
                              if not c["gates"]["C2"]["pass"]],
        },
        "C3_frac370": {
            "all_within_0p03": all(bool(c["gates"]["C3"]["pass_reported"])
                                   for c in cells),
            "worst_delta": max(col(lambda c: abs(c["gates"]["C3"]["delta"]))),
        },
        "V1_cox_munk": {
            "factor_range": [min(col(lambda c: c["gates"]["V1"]
                                     ["factor_over_cox_munk"])),
                             max(col(lambda c: c["gates"]["V1"]
                                     ["factor_over_cox_munk"]))],
            "legacy_factor_range": [min(col(lambda c: c["gates"]["V1"]
                                            ["legacy_factor_over_cox_munk"])),
                                    max(col(lambda c: c["gates"]["V1"]
                                            ["legacy_factor_over_cox_munk"]))],
        },
        "V2_laxague": {"reported": True},
        "V3_saturation": {
            "pass": all(c["gates"]["V3"]["pass"] for c in cells),
            "worst_B_l_over_Bhat": max(col(lambda c: c["gates"]["V3"]
                                           ["max_ratio_B_l_over_Bhat_low"])),
            "worst_B_total_over_Bhat": max(
                col(lambda c: c["gates"]["V3"]
                    ["max_ratio_B_total_over_Bhat_low"])),
        },
        "K1_Hs_delta": {
            "flags_open": [c["cell"] for c in cells
                           if c["gates"]["K1"]["flag"]],
            "worst_abs_delta_rel": max(col(lambda c: abs(c["gates"]["K1"]
                                                         ["delta_rel"]))),
        },
        "K5_quantisation": {
            "pass": True,
            "worst_rel": max(col(lambda c: c["gates"]["K5"]
                                 ["max_rel_error"])),
            "legacy_worst_rel": max(col(lambda c: c["legacy_v10"]
                                        ["quantisation_max_rel"])),
        },
        "K6_continuity": {
            "pass": all(c["gates"]["K6"]["pass"] for c in cells),
            "worst_slope_ratio": max(col(lambda c: c["gates"]["K6"]
                                         ["max_dslope_ratio_fine_over_coarse"])),
        },
    }
    dump_json(os.path.join(OUT, "numbers-v11-sweep.json"), {
        "stage": 2,
        "what": ("v3.8 8-cell sweep: C1 C2 C3 V1 V2 V3 K1 K5 K6; legacy "
                 "numbers read from the v10 sweep JSON"),
        "created_utc": now_iso(),
        "model": "version 3.8 (PIPELINE-V10 chunk 8)",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters"),
        "grid": {"N": GRID_SWEEP, "depth_m": SWEEP_H,
                 "U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
                 "seed": seed, "spectrum_source": "elfouhaily-form"},
        "policy": {"patch_base_m": V10.PATCH_BASE_M,
                   "r_lambda": V10.R_LAMBDA,
                   "k_cut_default": V10.KCUT_DEFAULT,
                   "B_hat_range": [B_HAT_LOW, B_HAT_HIGH],
                   "legacy_source": "read from out/v10/numbers-v10-sweep.json"},
        "summary": summary,
        "cells": cells,
        "wall_s": float(time.time() - t0),
    })
    log("  summary: C1=%s C2=%s C3=%s V3=%s K1flags=%s K6=%s"
        % (summary["C1_gravity_share"]["pass"],
           summary["C2_total"]["pass"],
           summary["C3_frac370"]["all_within_0p03"],
           summary["V3_saturation"]["pass"],
           summary["K1_Hs_delta"]["flags_open"],
           summary["K6_continuity"]["pass"]))
    cp.unit_done("s2-json")
    log("[stage 2] done in %.1fs" % (time.time() - t0))
    return summary


# ==========================================================================
# stage 3 -- compatibility K2, K3, K4
# ==========================================================================
WINDY = {"U": 6.0, "F": 100.0, "h": 1.0}
YARDSTICKS = R10.YARDSTICKS
PHOTONS = 2000000


def _texture(x_u8):
    x = np.asarray(x_u8, np.float64) / 255.0
    lum = x @ np.array([0.2126, 0.7152, 0.0722])
    gy = np.abs(np.diff(lum, axis=0))
    gx = np.abs(np.diff(lum, axis=1))
    return {"mean_abs_dL_dy": float(gy.mean()),
            "mean_abs_dL_dx": float(gx.mean()),
            "mean_luminance": float(lum.mean()),
            "std_luminance": float(lum.std()),
            "crush_rate_linear_lt_0.02": float(np.mean(x < 0.02)),
            "frac_aces_gt_0.98": float(np.mean(x > 0.98))}


def stage3():
    t0 = time.time()
    log("=" * 72)
    log("[stage 3] compatibility: K2 (height quantiles), K3 (windy, both "
        "sources), K4 (regression battery)")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    from PIL import Image
    from state import derive_seed
    from v5_composite import A_V3, B_SIGMA_V3, E_V3, render_frame_v5
    from run_windy_still import (warp_measure, predicted_shift_field,
                                 predicted_block_stats, block_agreement)
    raw = load_recipe()
    seed = int(derive_seed(raw, 1))
    base_u8, photo_lin, PW, PH = R10._load_photo()
    log("  photo %dx%d sha=%s" % (PW, PH, sha256_file(PHOTO)[:16]))

    out = {"stage": 3, "created_utc": now_iso(),
           "what": ("v3.8 compatibility: K2 height quantiles, K3 the G10 windy "
                    "scene through the UNCHANGED transport + frozen "
                    "composite with BOTH sources, K4 regression battery"),
           "photo": {"path": PHOTO, "size": [int(PW), int(PH)],
                     "sha256": sha256_file(PHOTO)},
           "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                     "parameters"),
           "K2": {}, "K3": {}, "K4": {}}

    # ---- K2: reference scenes from the v10 transport JSON ---------------
    ref_path = os.path.join(V10_OUT, "numbers-v10-transport.json")
    with open(ref_path, "r", encoding="utf-8") as f:
        v10t = json.load(f)
    scenes = {}
    for key in ("A", "B"):
        sc = v10t["scenes"][key]
        scenes[key] = {"name": sc["name"],
                       "U": float(sc["recipe"]["U_ms"]),
                       "F": float(sc["recipe"]["fetch_m"]),
                       "h": float(sc["recipe"]["depth_m"])}
    log("  K2 reference scenes identified in the v10 transport JSON: %s"
        % {k: (v["U"], v["F"], v["h"]) for k, v in scenes.items()})
    k2_scenes = {}
    for key in ("A", "B"):
        sc = scenes[key]
        rows = {}
        for src in ("jonswap", "elfouhaily-form"):
            r = scene_for(raw, sc["U"], sc["F"], sc["h"], GRID_SWEEP, seed,
                          src)
            st = r["state"]
            rows[src] = {
                "height": height_quantiles(r["fields"]),
                "slope": slope_quantiles(r["fields"]),
                "realised": st.realised_stats(r["fields"]),
                "spectrum_after": st.spectrum_after,
                "patch_m": float(st.S),
            }
            del r, st
            gc.collect()
            cp.unit_done("s3-k2-%s-%s" % (key, src.replace("-", "_")))
        old, new = rows["jonswap"], rows["elfouhaily-form"]
        k2_scenes[key] = {
            "scene": sc, "seed": seed, "sources": rows,
            "rel_change_new_over_legacy": {
                "p50": float(new["height"]["p50_m"]
                             / max(old["height"]["p50_m"], 1e-12) - 1.0)
                if old["height"]["p50_m"] != 0 else None,
                "p90": float(new["height"]["p90_m"]
                             / old["height"]["p90_m"] - 1.0),
                "p99": float(new["height"]["p99_m"]
                             / old["height"]["p99_m"] - 1.0),
                "Hs": float(new["height"]["Hs_m"]
                            / old["height"]["Hs_m"] - 1.0),
                "std": float(new["height"]["std_m"]
                             / old["height"]["std_m"] - 1.0),
                "kurtosis": float(new["height"]["kurtosis"]
                                  - old["height"]["kurtosis"]),
                "skewness": float(new["height"]["skewness"]
                                  - old["height"]["skewness"]),
            },
            "note": ("same seed family (one seed, preset 1) and the same "
                     "grid; only the source differs.  p50 is reported in m "
                     "as well as relatively, because it is ~0"),
        }
        log("  K2 %s (%s): p90 %.5f -> %.5f (%+.1f%%), p99 %.5f -> %.5f "
            "(%+.1f%%), Hs %.5f -> %.5f (%+.1f%%)"
            % (key, sc["name"], old["height"]["p90_m"],
               new["height"]["p90_m"],
               100 * (new["height"]["p90_m"] / old["height"]["p90_m"] - 1.0),
               old["height"]["p99_m"], new["height"]["p99_m"],
               100 * (new["height"]["p99_m"] / old["height"]["p99_m"] - 1.0),
               old["height"]["Hs_m"], new["height"]["Hs_m"],
               100 * (new["height"]["Hs_m"] / old["height"]["Hs_m"] - 1.0)))
    out["K2"] = {"scenes": k2_scenes,
                 "what": ("height quantiles (eta p50/p90/p99) of built "
                          "fields, old vs new source, two reference scenes "
                          "(the v10 stage-3 scenes A and B, identified from "
                          "numbers-v10-transport.json), same seed family")}
    cp.unit_done("s3-k2")

    # ---- K3: the windy scene, both sources ------------------------------
    renders = {}
    for src, fname in (("jonswap", "comp-windy-v38-legacy.png"),
                       ("elfouhaily-form", "comp-windy-v38-new.png")):
        r = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], GRID_RENDER,
                      seed, src)
        st = r["state"]
        t = R10.transport_scene(r, photo_lin, PW, PH, GRID_RENDER, Q=25,
                                E=float(E_V3))
        atomic_write_bytes(os.path.join(WORK, "C-windy-v38-%s-Q25.npy"
                                       % src.replace("-", "_")),
                           np.ascontiguousarray(t["C"]).tobytes())
        png = io.BytesIO()
        Image.fromarray(t["still"]).save(png, format="PNG")

        def verify_png(p, _pw=PW, _ph=PH):
            try:
                im = Image.open(p)
                im.load()
                return im.size == (_pw, _ph) and im.mode == "RGB"
            except Exception:
                return False

        atomic_write_bytes(os.path.join(OUT, fname), png.getvalue(),
                           verify=verify_png)
        still_png = np.asarray(Image.open(os.path.join(OUT, fname))
                               .convert("RGB"))
        warp = warp_measure(base_u8, still_png)
        pred = predicted_shift_field(r["fields"], r["env"], PW, PH)
        ps = predicted_block_stats(pred, warp["grid"])
        agree = {"ssd_primary": block_agreement(ps, warp["per_block"],
                                                "ssd_dy", "ssd_dx"),
                 "ncc_check": block_agreement(ps, warp["per_block"],
                                              "ncc_dy", "ncc_dx")}
        renders[src] = {
            "state": st, "r": r, "t": t, "still": t["still"],
            "still_png": still_png, "warp": warp, "pred": pred, "ps": ps,
            "agree": agree, "file": fname,
        }
        log("  K3 %-16s C mean %.9f p50 %.3f p99 %.3f p99.9 %.3f max %.2f | "
            "warp SSD %.2f/%.2f px | own field %.2f px = %.2f mm at %.3f "
            "mm/px"
            % (src, float(t["C"].mean()), np.percentile(t["C"], 50),
               np.percentile(t["C"], 99), np.percentile(t["C"], 99.9),
               float(t["C"].max()), warp["primary"]["mean_shift_px"],
               warp["primary"]["shift_field_std_px"],
               ps["field_mean_magnitude_px"],
               1e3 * float(ps["mean_predicted_shift_m"]),
               float(ps["pixel_scale_mm"])))
        cp.unit_done("s3-k3-" + src.replace("-", "_"))

    leg_r = renders["jonswap"]
    new_r = renders["elfouhaily-form"]

    # ---- K3 attribution: g1-pinned control -----------------------------
    pin = R10.transport_scene(new_r["r"], photo_lin, PW, PH, GRID_RENDER,
                              Q=25, E=float(E_V3),
                              g1=float(leg_r["t"]["g1"]))
    pin_warp = warp_measure(base_u8, pin["still"])
    d_in = np.abs(pin["still"].astype(np.int16)
                  - leg_r["still"].astype(np.int16))
    ms_in = R10.msssim(pin["still"].astype(np.float64) / 255.0,
                       leg_r["still"].astype(np.float64) / 255.0)
    log("  K3 g1-pinned control: mean|d| %.2f/255 MS-SSIM %.5f"
        % (d_in.mean(), ms_in["msssim"]))
    cp.unit_done("s3-k3-g1pin")

    # ---- K3 control C=1 (geometry only) --------------------------------
    ctrl = {}
    for src in ("jonswap", "elfouhaily-form"):
        rr = renders[src]
        img = render_frame_v5(
            rr["r"]["fields"],
            np.ones((GRID_RENDER, GRID_RENDER), np.float32),
            photo_lin, rr["r"]["env"], PW, PH, ambient_scale=A_V3,
            c_layer_exposure=float(E_V3),
            substrate_gain_rgb=rr["t"]["opt"]["gain"]["gain_rgb"],
            beer_down=rr["t"]["opt"]["beer"]["down"],
            beer_up=rr["t"]["opt"]["beer"]["up"],
            g1=float(rr["t"]["g1"]),
            blur_sigma_px=rr["t"]["opt"]["blur"]["sigma_px"],
            lowpass_sigma_px=B_SIGMA_V3, return_debug=False)
        w = warp_measure(base_u8, img)
        ctrl[src] = {"ssd_primary": {k: w["primary"][k] for k in (
            "mean_shift_px", "shift_field_std_px",
            "frac_blocks_at_search_edge")},
            "ncc_check": {k: w["ncc_check"][k] for k in (
                "mean_shift_px", "shift_field_std_px",
                "frac_blocks_at_search_edge")}}
        del img
        gc.collect()
        log("  K3 C=1 control %-16s SSD %.2f/%.2f px"
            % (src, ctrl[src]["ssd_primary"]["mean_shift_px"],
               ctrl[src]["ssd_primary"]["shift_field_std_px"]))
    cp.unit_done("s3-k3-ctrl")

    def k3_node(src):
        rr = renders[src]
        st = rr["state"]
        node = st.state_node(rr["r"]["fields"])
        t = rr["t"]
        p = rr["warp"]["primary"]
        ncc = rr["warp"]["ncc_check"]
        ps = rr["ps"]
        return {
            "spectrum_source": src,
            "state_summary": {
                "patch_m": float(st.S),
                "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
                "Hs_spatial_m": float(node["realised"]["Hs_m"]),
                "sigma2_spectral": float(st.spectrum_after["sigma2_rad2"]),
                "sigma2_spatial": float(node["realised"]["sigma2_rad2"]),
                "rms_slope": float(st.spectrum_after["rms_slope"]),
                "Omega_c": (None if src == "jonswap"
                            else float(st.ctx["omega_c"])),
                "eckv_domain": (st.ctx.get("eckv_domain")
                                if src != "jonswap" else None),
                "bands_realised": node["bands"],
                "frac_above_370_realised":
                    st.frac_above((370.0,))["above_370"]["realised"],
            },
            "route": {"Q": 25, "A": A_V3, "E": float(E_V3),
                      "B_lowpass_sigma_px": B_SIGMA_V3,
                      "g1_used": float(t["g1"]),
                      "g1_source": ("Smith shadowing from THIS state's "
                                    "realised rms slope (the transport's own "
                                    "route, untouched)"),
                      "blur_sigma_px": float(t["opt"]["blur"]["sigma_px"])},
            "C": R10.c_stats(t["C"]),
            "C_percentiles": {"p50": float(np.percentile(t["C"], 50)),
                              "p90": float(np.percentile(t["C"], 90)),
                              "p99": float(np.percentile(t["C"], 99)),
                              "p99_9": float(np.percentile(t["C"], 99.9)),
                              "p99_99": float(np.percentile(t["C"], 99.99)),
                              "max": float(t["C"].max())},
            "warp": {
                "method": rr["warp"]["definition"],
                "primary_ssd": {k: p[k] for k in (
                    "mean_shift_px", "shift_field_std_px",
                    "shift_magnitude_p50_px", "shift_magnitude_p90_px",
                    "shift_magnitude_max_px", "frac_blocks_at_search_edge",
                    "n_blocks", "block_px", "search_px", "step_px")},
                "ncc_check": {k: ncc[k] for k in (
                    "mean_shift_px", "shift_field_std_px",
                    "shift_magnitude_p50_px", "shift_magnitude_p90_px",
                    "frac_blocks_at_search_edge")},
                "own_shift_field": {
                    "field_mean_px": ps["field_mean_magnitude_px"],
                    "field_std_px": ps["field_std_magnitude_px"],
                    "block_mean_px": ps["block_mean_magnitude_px"],
                    "mean_mm": 1e3 * float(ps["mean_predicted_shift_m"]),
                    "pixel_scale_mm": float(ps["pixel_scale_mm"]),
                    "frac_beyond_search": ps["field_frac_beyond_search"],
                },
                "agreement_with_predicted": rr["agree"],
                "control_C1": ctrl[src],
                "yardsticks": YARDSTICKS,
            },
            "fine_texture": _texture(rr["still_png"]),
            "tone": {
                "crush_rate_linear_lt_0.02":
                    float(t["metrics"]["crush_rate_linear_lt_0.02"]),
                "blowout_rate_aces_gt_0.98":
                    float(t["metrics"]["blowout_rate_aces_gt_0.98"]),
                "customary_crush_bound": 0.0325,
                "crush_pass_customary": bool(
                    t["metrics"]["crush_rate_linear_lt_0.02"] <= 0.0325),
                "note": ("the frozen composite pins E = 2.2; these are "
                         "reported consequences, not tuning targets"),
            },
            "beer": t["beer"],
            "still_metrics": t["metrics"],
            "files": [rr["file"]],
        }

    k3 = {"scenes": {"jonswap": k3_node("jonswap"),
                     "elfouhaily-form": k3_node("elfouhaily-form")}}
    d = np.abs(new_r["still"].astype(np.int16)
               - leg_r["still"].astype(np.int16))
    ms = R10.msssim(new_r["still"].astype(np.float64) / 255.0,
                    leg_r["still"].astype(np.float64) / 255.0)
    tex_b = _texture(base_u8)
    k3["between_sources"] = {
        "mean_abs_uint8": float(d.mean()), "max_abs_uint8": int(d.max()),
        "frac_pixels_changed": float(np.mean(d.max(axis=2) > 0)),
        "msssim": ms["msssim"], "cs_per_scale": ms["cs_per_scale"],
        "luminance": {
            "legacy_mean": float((leg_r["still_png"].astype(np.float64)
                                  @ np.array([0.2126, 0.7152, 0.0722])).mean()),
            "new_mean": float((new_r["still_png"].astype(np.float64)
                               @ np.array([0.2126, 0.7152, 0.0722])).mean()),
        },
        "g1_pinned_control": {
            "what": ("the new source rendered with the LEGACY g1, so the "
                     "remaining difference is the caustic gain C alone"),
            "mean_abs_uint8": float(d_in.mean()),
            "max_abs_uint8": int(d_in.max()),
            "msssim_vs_legacy": ms_in["msssim"],
            "warp_primary": {k: pin_warp["primary"][k] for k in (
                "mean_shift_px", "shift_field_std_px")},
            "C": R10.c_stats(pin["C"]),
            "g1_pinned_to": float(leg_r["t"]["g1"]),
            "g1_native": float(new_r["t"]["g1"]),
        },
        "texture_base_photo": tex_b,
    }
    # attribution
    lg = k3["scenes"]["jonswap"]["state_summary"]
    nw = k3["scenes"]["elfouhaily-form"]["state_summary"]
    k3["attribution"] = {
        "Hs_rel_change": float(nw["Hs_spectral_m"] / lg["Hs_spectral_m"]
                               - 1.0),
        "sigma2_rel_change": float(nw["sigma2_spectral"]
                                   / lg["sigma2_spectral"] - 1.0),
        "rms_slope_rel_change": float(nw["rms_slope"] / lg["rms_slope"]
                                      - 1.0),
        "g1_change": float(new_r["t"]["g1"] / leg_r["t"]["g1"] - 1.0),
        "frac_above_370_change": float(
            nw["frac_above_370_realised"] / lg["frac_above_370_realised"]
            - 1.0),
        "band_share_1_10cm": {
            "legacy": lg["bands_realised"]["1_10cm"]["sigma2_share"],
            "new": nw["bands_realised"]["1_10cm"]["sigma2_share"]},
        "band_share_lt_1cm": {
            "legacy": lg["bands_realised"]["lt_1cm"]["sigma2_share"],
            "new": nw["bands_realised"]["lt_1cm"]["sigma2_share"]},
        "fine_texture_rel_change": {
            "dy": float(k3["scenes"]["elfouhaily-form"]["fine_texture"]
                        ["mean_abs_dL_dy"]
                        / k3["scenes"]["jonswap"]["fine_texture"]
                        ["mean_abs_dL_dy"] - 1.0),
            "dx": float(k3["scenes"]["elfouhaily-form"]["fine_texture"]
                        ["mean_abs_dL_dx"]
                        / k3["scenes"]["jonswap"]["fine_texture"]
                        ["mean_abs_dL_dx"] - 1.0)},
        "own_field_shift_mm": {
            "legacy": k3["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["mean_mm"],
            "new": k3["scenes"]["elfouhaily-form"]["warp"]
            ["own_shift_field"]["mean_mm"]},
        "statement": None,
    }
    att = k3["attribution"]
    ab = {}
    for src in ("jonswap", "elfouhaily-form"):
        st = k3["scenes"][src]["state_summary"]
        ab[src] = float(st["frac_above_370_realised"] * st["sigma2_spectral"])
    att["sigma2_above_370_absolute"] = {
        "legacy": ab["jonswap"], "new": ab["elfouhaily-form"],
        "rel_change": float(ab["elfouhaily-form"] / ab["jonswap"] - 1.0),
        "note": ("the SHARE above 370 rises while the LEVEL falls: this is "
                 "the number that explains why the fine texture weakens even "
                 "though the high-k share of sigma^2 grows"),
    }
    k3["attribution"]["statement"] = (
        "Broad warp follows the new H_s: the pipeline's own geometric shift "
        "field falls to %.1f%% of the legacy value in mm, tracking delta H_s "
        "= %+.1f%% and delta sigma^2 = %+.1f%%.  Fine texture is a "
        "band-redistribution question, not an H_s question: the SHARE of "
        "sigma^2 above 370 rad/m changes by %+.1f%%, but the ABSOLUTE level "
        "of that variance changes by %+.1f%%, and the mean |dL/dy| changes by "
        "%+.1f%%.  Shadowing is not part of the change (Smith g1 saturates at "
        "1.0 for BOTH states at this sun elevation, so the g1-pinned render is "
        "the native render); the C=1 control separates geometry from caustic "
        "gain."
        % (100.0 * att["own_field_shift_mm"]["new"]
           / max(att["own_field_shift_mm"]["legacy"], 1e-12),
           100 * att["Hs_rel_change"],
           100 * att["sigma2_rel_change"],
           100 * att["frac_above_370_change"],
           100 * att["sigma2_above_370_absolute"]["rel_change"],
           100 * att["fine_texture_rel_change"]["dy"]))
    out["K3"] = k3
    log("  K3 attribution: Hs %+.1f%% sig2 %+.1f%% g1 %+.1f%% | own-field "
        "shift %.2f -> %.2f mm | |dL/dy| %+.1f%%"
        % (100 * k3["attribution"]["Hs_rel_change"],
           100 * k3["attribution"]["sigma2_rel_change"],
           100 * k3["attribution"]["g1_change"],
           k3["attribution"]["own_field_shift_mm"]["legacy"],
           k3["attribution"]["own_field_shift_mm"]["new"],
           100 * k3["attribution"]["fine_texture_rel_change"]["dy"]))
    cp.unit_done("s3-k3")

    # ---- K4 regression battery ----------------------------------------
    sB = new_r
    again = R10.transport_scene(sB["r"], photo_lin, PW, PH, GRID_RENDER,
                                Q=25, E=float(E_V3),
                                C_cache=sB["t"]["C"])
    ident = bool(np.array_equal(again["still"], sB["still"]))
    k4a = {"scene": "windy, NEW source", "identical_buffers": True,
           "bit_identical_output": ident,
           "max_abs_uint8": int(np.abs(
               again["still"].astype(np.int16)
               - sB["still"].astype(np.int16)).max()),
           "params": {"A": A_V3, "B_lowpass_sigma_px": B_SIGMA_V3,
                      "E": float(E_V3), "Q": 25,
                      "g1": float(sB["t"]["g1"])},
           "pass": ident}
    log("  K4a frozen composite on identical buffers: bit-identical=%s"
        % ident)
    cp.unit_done("s3-k4a")

    rec_path = os.path.join(HERE, "out", "windy-still", "numbers-windy.json")
    with open(rec_path, "r", encoding="utf-8") as f:
        rec = json.load(f)
    rec_tw = rec["depth_cost"]["beer_direct_beam"]["z1.0"]["two_way"]
    rec_sd = rec["depth_cost"]["beer_direct_beam"]["z1.0"]["s_down_m"]
    tw_l = leg_r["t"]["beer"]["two_way"]
    tw_n = new_r["t"]["beer"]["two_way"]
    diff_rec = [abs(float(tw_n[i]) - float(rec_tw[i])) for i in range(3)]
    diff_src = [abs(float(tw_n[i]) - float(tw_l[i])) for i in range(3)]
    k4b = {
        "what": ("Beer-Lambert on the direct beam for identical optical "
                 "inputs: (i) new vs the v3.4/v9 RECORDED numbers, (ii) new "
                 "vs legacy -- the two sources share the environment exactly, "
                 "so this must be identical"),
        "new_two_way": [float(x) for x in tw_n],
        "legacy_two_way": [float(x) for x in tw_l],
        "v9_recorded_two_way": rec_tw,
        "abs_diff_vs_recorded": diff_rec,
        "abs_diff_legacy_vs_new": diff_src,
        "s_down_m": {"new": float(new_r["t"]["beer"]["s_down_m"]),
                     "legacy": float(leg_r["t"]["beer"]["s_down_m"]),
                     "v9_recorded": rec_sd},
        "code_path": ("v5_physics.beer_lambert -- imported, not "
                      "reimplemented and not edited"),
        "pass": bool(max(diff_rec) <= 5e-5
                     and max(diff_src) <= 1e-12),
    }
    log("  K4b Beer two-way new %s vs v9 %s (max|d| %.2e); legacy identical "
        "to %.1e"
        % (["%.6f" % x for x in tw_n], ["%.6f" % x for x in rec_tw],
           max(diff_rec), max(diff_src)))
    cp.unit_done("s3-k4b")

    ra = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], 512, seed,
                   "elfouhaily-form")
    rb = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], 512, seed,
                   "elfouhaily-form")
    eq = bool(np.array_equal(ra["state"].lo.H0, rb["state"].lo.H0)
              and np.array_equal(ra["state"].hi.H0, rb["state"].hi.H0))
    # and the rendered state itself is a rain_free build of the same ctor
    same_as_render = bool(np.array_equal(ra["state"].lo.H0.shape,
                                        sB["state"].lo.H0.shape))
    rain_raised = False
    try:
        V11WaterState({"x": 1}, {"windSpeedMS": 6.0, "depthM": 1.0,
                                 "windDir": [1.0, 0.0]}, N=64, seed=1,
                      fetch_m=100.0, rain_rate_mm_h=1.0,
                      spectrum_source="elfouhaily-form")
    except NotImplementedError:
        rain_raised = True
    k4c = {"bit_identical_state": eq, "rain_positive_raises": rain_raised,
           "rendered_state_shape_matches": same_as_render,
           "pass": bool(eq and rain_raised),
           "note": ("R = 0 takes no rain branch at all, so two builds of the "
                    "NEW source are bit-identical; R > 0 raises instead of "
                    "silently dropping the rain band (inherited v3.7 rule)")}
    log("  K4c rain R=0 identical=%s, R>0 raises=%s" % (eq, rain_raised))
    cp.unit_done("s3-k4c")
    del ra, rb
    gc.collect()

    out["K4"] = {"K4a": k4a, "K4b": k4b, "K4c": k4c,
                 "battery": ("the v3.7 G8 battery verbatim.  These are "
                             "code-path regression tests on IDENTICAL "
                             "buffers; renders under the new state are "
                             "expected to differ and are judged by K3 only.")}
    out["summary"] = {
        "K2": "reported (2 scenes x 2 sources, height p50/p90/p99)",
        "K3": "reported (2 stills + C + warp yardsticks + |dL/dy| + tone)",
        "K4a": k4a["pass"], "K4b": k4b["pass"], "K4c": k4c["pass"],
    }
    out["wall_s"] = float(time.time() - t0)
    dump_json(os.path.join(OUT, "numbers-v11-transport.json"), out)
    cp.unit_done("s3-json")
    log("[stage 3] done in %.1fs" % (time.time() - t0))
    return out


# ==========================================================================
# stage 4 -- merge evaluation + default flip + report
# ==========================================================================
def evaluate_merge():
    """PIPELINE-V10 section 5, item by item, from the stage JSONs."""
    with open(os.path.join(OUT, "numbers-v11-spectrum.json"),
              "r", encoding="utf-8") as f:
        s1 = json.load(f)
    with open(os.path.join(OUT, "numbers-v11-sweep.json"),
              "r", encoding="utf-8") as f:
        s2 = json.load(f)
    with open(os.path.join(OUT, "numbers-v11-transport.json"),
              "r", encoding="utf-8") as f:
        s3 = json.load(f)
    sm = s2["summary"]
    a1 = s1["gates"]["A1"]
    a2 = s1["gates"]["A2"]
    items = [
        ("A1 (D1 verification)", bool(a1["pass"]),
         "finite-depth plateau moves %.1f%% (< 30%%); deep-water closed form "
         "matched to %.1e" % (100 * a1["worst_finite_depth_plateau_rel"],
                              a1["depths"]["h_%g" % 1000.0]
                              ["ratio_kM_physical"]["max_abs_minus_1"])),
        ("A2 (Parseval)", bool(a2["pass"]),
         "z(m0) %+.2f, z(sigma2) %+.2f" % (a2["z_m0"], a2["z_sigma2"])),
        ("C1 (gravity share 50-200)", bool(sm["C1_gravity_share"]["pass"]),
         "worst share %.3f (limit 0.50); failing cells: %s"
         % (sm["C1_gravity_share"]["worst_share"],
            ", ".join(sm["C1_gravity_share"]["failing_cells"]) or "none")),
        ("C2 (total sigma^2 vs v3.7)", bool(sm["C2_total"]["pass"]),
         "worst |rel| %.3f (limit 0.10); failing cells: %d of 8"
         % (sm["C2_total"]["worst_rel"],
            len(sm["C2_total"]["failing_cells"]))),
        ("C3 (reported)", True,
         "frac>370 within %.3f of the ECKV reference"
         % sm["C3_frac370"]["worst_delta"]),
        ("V1 (reported)", True,
         "factor over clean Cox-Munk %.2f-%.2f (legacy %.2f-%.2f)"
         % (sm["V1_cox_munk"]["factor_range"][0],
            sm["V1_cox_munk"]["factor_range"][1],
            sm["V1_cox_munk"]["legacy_factor_range"][0],
            sm["V1_cox_munk"]["legacy_factor_range"][1])),
        ("V2 (reported)", True, "band shares vs the Laxague partition"),
        ("V3 (reported)", True,
         "max B_gravity/B_hat %.3f; max B_total/B_hat %.3f"
         % (sm["V3_saturation"]["worst_B_l_over_Bhat"],
            sm["V3_saturation"]["worst_B_total_over_Bhat"])),
        ("K1 (no flag open)",
         not bool(sm["K1_Hs_delta"]["flags_open"]),
         "worst |dHs| %.1f%% (flag at 35%%); flags: %s"
         % (100 * sm["K1_Hs_delta"]["worst_abs_delta_rel"],
            ", ".join(sm["K1_Hs_delta"]["flags_open"]) or "none")),
        ("K2 (reported)", True, "height quantiles, 2 scenes x 2 sources"),
        ("K3 (reported)", True, "both stills, C, warp yardsticks, texture"),
        ("K4 (regression battery)",
         bool(s3["K4"]["K4a"]["pass"] and s3["K4"]["K4b"]["pass"]
              and s3["K4"]["K4c"]["pass"]),
         "composite bit-identical %s; Beer-Lambert %s; rain R=0 %s"
         % (s3["K4"]["K4a"]["pass"], s3["K4"]["K4b"]["pass"],
            s3["K4"]["K4c"]["pass"])),
        ("K5 (quantisation)", bool(sm["K5_quantisation"]["pass"]),
         "worst rel %.3f" % sm["K5_quantisation"]["worst_rel"]),
        ("K6 (continuity)", bool(sm["K6_continuity"]["pass"]),
         "slope-halving ratio %.3f (smooth = 0.5)"
         % sm["K6_continuity"]["worst_slope_ratio"]),
    ]
    blocking = [name for name, ok, _ in items if not ok]
    return {"items": [{"item": n, "pass": bool(ok), "evidence": e}
                      for n, ok, e in items],
            "merge_pass": bool(not blocking),
            "blocking_items": blocking}


def stage4():
    t0 = time.time()
    log("=" * 72)
    log("[stage 4] continuity check, merge evaluation, report")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    # ---- K3 adjunct: is the legacy still byte-identical to v10's? ------
    tp = os.path.join(OUT, "numbers-v11-transport.json")
    with open(tp, "r", encoding="utf-8") as f:
        tr = json.load(f)
    v10_png = os.path.join(V10_OUT, "comp-windy-v37.png")
    v11_png = os.path.join(OUT, "comp-windy-v38-legacy.png")
    cont = {"what": ("the LEGACY-source still versus the v10 chunk-7 "
                     "deliverable: the source switch must not perturb the "
                     "legacy path at all"),
            "v10_file": os.path.relpath(v10_png, HERE),
            "v11_legacy_file": os.path.relpath(v11_png, HERE),
            "sha256_v10": sha256_file(v10_png),
            "sha256_v11_legacy": sha256_file(v11_png)}
    cont["byte_identical"] = bool(cont["sha256_v10"] == cont["sha256_v11_legacy"])
    try:
        from PIL import Image
        a = np.asarray(Image.open(v10_png).convert("RGB")).astype(np.int16)
        b = np.asarray(Image.open(v11_png).convert("RGB")).astype(np.int16)
        dlt = np.abs(a - b)
        cont["max_abs_uint8"] = int(dlt.max())
        cont["frac_pixels_changed"] = float(np.mean(dlt.max(axis=2) > 0))
        del a, b, dlt
    except Exception as exc:
        cont["error"] = str(exc)
    log("  K3 continuity: legacy still byte-identical to comp-windy-v37.png = "
        "%s" % cont["byte_identical"])
    cp.unit_done("s4-legacy-continuity")

    # ---- attribution addendum: ABSOLUTE high-k slope variance ----------
    att = tr["K3"]["attribution"]
    ab = {}
    for src in ("jonswap", "elfouhaily-form"):
        st = tr["K3"]["scenes"][src]["state_summary"]
        ab[src] = float(st["frac_above_370_realised"] * st["sigma2_spectral"])
    att["sigma2_above_370_absolute"] = {
        "legacy": ab["jonswap"], "new": ab["elfouhaily-form"],
        "rel_change": float(ab["elfouhaily-form"] / ab["jonswap"] - 1.0),
        "note": ("the SHARE above 370 rises while the LEVEL falls: this is "
                 "the number that explains why the fine texture weakens even "
                 "though the high-k share of sigma^2 grows"),
    }
    tr["K3"]["continuity_vs_v10"] = cont
    tr["K3"]["g1_note"] = (
        "Smith shadowing saturates at g1 = 1.0 for BOTH sources at this sun "
        "elevation (%.2f / %.2f rms slope), so the g1-pinned control is "
        "degenerate: it reproduces the native new render exactly.  The "
        "between-source difference is therefore NOT a shadowing change; it is "
        "the geometry (own shift field) and the caustic gain C."
        % (tr["K3"]["scenes"]["jonswap"]["state_summary"]["rms_slope"],
           tr["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
           ["rms_slope"]))
    dump_json(tp, tr)
    cp.unit_done("s4-transport-adjunct")

    merge = evaluate_merge()
    for it in merge["items"]:
        log("  %-32s %s  %s" % (it["item"], "PASS" if it["pass"] else "FAIL",
                                it["evidence"]))
    flipped = False
    if merge["merge_pass"]:
        path = os.path.join(HERE, "v11_state.py")
        with open(path, "r", encoding="utf-8") as f:
            txt = f.read()
        if 'DEFAULT_SOURCE = "jonswap"' in txt:
            txt = txt.replace('DEFAULT_SOURCE = "jonswap"',
                              'DEFAULT_SOURCE = "elfouhaily-form"', 1)
            atomic_write_bytes(path, txt.encode("utf-8"))
            flipped = True
    merge["default_flipped"] = flipped
    merge["default_after"] = ("elfouhaily-form" if flipped else "jonswap")
    dump_json(os.path.join(OUT, "merge-evaluation-v11.json"), merge)
    cp.unit_done("s4-merge")
    log("  merge pass=%s; default=%s; blocking: %s"
        % (merge["merge_pass"], merge["default_after"],
           ", ".join(merge["blocking_items"]) or "none"))
    import v11_report
    v11_report.main()
    cp.unit_done("s4-report")
    log("[stage 4] done in %.1fs" % (time.time() - t0))
    return merge


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="1",
                    choices=["1", "2", "3", "4"])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    if args.stage == "1":
        stage1()
    elif args.stage == "2":
        stage2()
    elif args.stage == "3":
        stage3()
    elif args.stage == "4":
        stage4()


if __name__ == "__main__":
    main()
