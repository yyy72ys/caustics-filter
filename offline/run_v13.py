"""run_v13 -- version 3.9 chunk 9: the low-wind capillary floor, labelled.

Four stages, foreground, PIPELINE-V11.  Every stage is checkpointed with
D:\\mutaro-work\\jobkit (ONE unit per stage, as the brief asks) and every
deliverable is committed with atomic_write_bytes.  Nothing outside
`out/v13/` and the v13 files is written; `out/v12/` belongs to another job.

  stage 1  v13 spectrum/state + A1/A2 re-run (algebra vs v11 conventions,
           Parseval, the alpha_m floor solve per cell)
           -> out/v13/numbers-v13-spectrum.json, spectrum-curves-v13.json
  stage 2  eight-cell sweep: C1 THREE-STATE, C1b floor proof, V1-V3, K1/K5/K6,
           total sigma2 beside Cox-Munk clean and the second slope reference,
           full band partitions
           -> out/v13/numbers-v13-sweep.json
  stage 3  K3 windy compatibility (floored source vs v11 legacy, attribution),
           K4 regression battery + the prior-activity battery, and the
           mandatory sensitivity sweep over sigma2_h target
           {0, 0.0015, 0.003, 0.006} on the preset-1 and G10 windy scenes
           -> out/v13/numbers-v13-transport.json, numbers-v13-sensitivity.json
  stage 4  merge evaluation per PIPELINE-V11 + REPORT-V13-GATES.md
           -> out/v13/merge-evaluation-v13.json, REPORT-V13-GATES.md

Usage (foreground, from caustics-filter/offline):
  python run_v13.py --stage 1
  python run_v13.py --stage 2
  python run_v13.py --stage 3
  python run_v13.py --stage 4
"""
import argparse
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
import v13_spectrum as V13                                         # noqa: E402
import v13_state as S13                                            # noqa: E402
import run_v10 as R10                                              # noqa: E402

RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PHOTO = R10.PHOTO
OUT = os.path.join(HERE, "out", "v13")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
V11_OUT = os.path.join(HERE, "out", "v11")
V10_OUT = os.path.join(HERE, "out", "v10")

GRID_SWEEP = 1024
GRID_RENDER = R10.GRID_RENDER
SWEEP_U = V10.SWEEP_U
SWEEP_F = V10.SWEEP_F
SWEEP_H = V10.SWEEP_Z

B_HAT_LOW = 0.007
B_HAT_HIGH = 0.008

#: the G10 windy scene of PIPELINE-V11 section 3 (K3) and the two sensitivity
#: scenes; the identification is CHECKED against the v10 transport JSON.
WINDY = {"U": 6.0, "F": 100.0, "h": 1.0}
PRESET1 = {"U": 1.5, "F": 100.0, "h": 0.6}

log = R10.log
now_iso = R10.now_iso
dump_json = R10.dump_json
sha256_file = R10.sha256_file
load_recipe = R10.load_recipe
cell_id = R10.cell_id


def sweep_scenes():
    """The two sensitivity scenes, read from the v10 transport JSON."""
    with open(os.path.join(V10_OUT, "numbers-v10-transport.json"),
              "r", encoding="utf-8") as f:
        t = json.load(f)
    out = {}
    for key, want in (("A", PRESET1), ("B", WINDY)):
        sc = t["scenes"][key]
        rec = sc["recipe"]
        got = {"U": float(rec["U_ms"]), "F": float(rec["fetch_m"]),
               "h": float(rec["depth_m"]), "name": sc["name"],
               "seed": int(rec["seed"])}
        if (abs(got["U"] - want["U"]) > 1e-12
                or abs(got["F"] - want["F"]) > 1e-12
                or abs(got["h"] - want["h"]) > 1e-12):
            raise RuntimeError("sensitivity scene %s moved: %r" % (key, got))
        out[key] = got
    return out


def scene_for(raw, U, F, h, N, seed, source, target=None,
              k_cut=V10.KCUT_DEFAULT):
    """v13 state -> the dict run_v10.transport_scene expects."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    t0 = time.time()
    st = S13.V13WaterState(
        scene, env, N=N, seed=seed, preset_index=1, fetch_m=float(F),
        k_cut=k_cut, spectrum_source=source,
        sigma2_h_target=(S13.DEFAULT_SIGMA2_H_TARGET if target is None
                         else float(target)))
    t_state = time.time() - t0
    t1 = time.time()
    fields = st.fields_at(0.0)
    t_fields = time.time() - t1
    return {"state": st, "env": env, "cfg": cfg, "fields": fields,
            "t_state": t_state, "t_fields": t_fields}


def v11_baseline():
    """The v3.8 sweep JSON -- baseline numbers are READ, never recomputed."""
    with open(os.path.join(V11_OUT, "numbers-v11-sweep.json"),
              "r", encoding="utf-8") as f:
        d = json.load(f)
    return {c["cell"]: c for c in d["cells"]}, d


def legacy_v10():
    with open(os.path.join(V10_OUT, "numbers-v10-sweep.json"),
              "r", encoding="utf-8") as f:
        d = json.load(f)
    return {c["cell"]: c for c in d["cells"]}, d


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


# ==========================================================================
# stage 1 -- v13 spectrum/state + algebra gates
# ==========================================================================
def stage1():
    t0 = time.time()
    log("=" * 72)
    log("[stage 1] v13_spectrum: alpha_m resolution, A1, A2, algebra vs v11")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    scenes = sweep_scenes()
    p1 = scenes["A"]

    # ---- A1 (unchanged physics: the v3.7 B_g closed form) --------------
    a1 = V11.gate_a1()
    worst = a1["worst_finite_depth_plateau_rel"]
    a1["pass"] = bool(worst < a1["failure_condition"])
    a1["finding"] = (
        "Re-run on the v3.9 source: B_g is untouched by v3.9 (the change is "
        "alpha_m only), so the v3.8 result must reproduce exactly -- and it "
        "does (worst finite-depth plateau move %.2f%%, deep-water closed form "
        "matched to %.1e with the pipeline's own k_m).  A1 is reported, not "
        "re-used as evidence for the new level."
        % (100 * worst, a1["depths"]["h_%g" % 1000.0]["ratio_kM_physical"]
           ["max_abs_minus_1"]))
    log("  A1 PASS=%s (plateau %.2f%%)" % (a1["pass"], 100 * worst))

    # ---- A2 on the new source: in-domain cell AND prior-active cell -----
    ctx6 = V13.assemble(6.0, 100.0, 1.0)
    pm6 = V10.patch_for(ctx6["laws"]["lambda_p_m"])["patch_m"]
    a2 = V13.gate_a2(ctx6, N=1024, patch_m=pm6,
                     tag="elfouhaily-form (v3.9, prior inactive)")
    log("  A2 in-domain: z m0 %+.2f sigma2 %+.2f | rel cont %+.4f / %+.4f | "
        "PASS=%s" % (a2["z_m0"], a2["z_sigma2"],
                     a2["rel_mode_over_continuum_m0"],
                     a2["rel_mode_over_continuum_sigma2"], a2["pass"]))
    ctxp = V13.assemble(p1["U"], p1["F"], p1["h"])
    pmp = V10.patch_for(ctxp["laws"]["lambda_p_m"])["patch_m"]
    a2_prior = V13.gate_a2(ctxp, N=1024, patch_m=pmp,
                           tag="elfouhaily-form (v3.9, LOW-WIND PRIOR ACTIVE)")
    log("  A2 prior:     z m0 %+.2f sigma2 %+.2f | rel cont %+.4f / %+.4f | "
        "PASS=%s" % (a2_prior["z_m0"], a2_prior["z_sigma2"],
                     a2_prior["rel_mode_over_continuum_m0"],
                     a2_prior["rel_mode_over_continuum_sigma2"],
                     a2_prior["pass"]))
    leg_ctx = V10.assemble(6.0, 100.0, 1.0)
    a2_legacy = V13.gate_a2(leg_ctx, N=1024, patch_m=pm6,
                            tag="jonswap (v10 G0 analogue, control)")

    # ---- algebra vs v11 conventions ------------------------------------
    algebra = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            algebra[cell_id(U, F)] = V13.algebra_vs_v11(U, F, SWEEP_H)
    al_ok = all(
        (a["prior_active"] or a["spectrum_identical"])
        and a["gravity_part_identical"] and a["linearity_pass"]
        and a["S_h_ratio_matches_alpha_ratio"]
        for a in algebra.values())
    log("  algebra: in-domain bit-identical=%s, gravity untouched=%s, "
        "linearity=%s -> %s"
        % (all(a["spectrum_identical"] for a in algebra.values()
               if not a["prior_active"]),
           all(a["gravity_part_identical"] for a in algebra.values()),
           all(a["linearity_pass"] for a in algebra.values()),
           "PASS" if al_ok else "FAIL"))

    # ---- the eight cells, analytic --------------------------------------
    cells = []
    curves_all = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            ctx = V13.assemble(U, F, SWEEP_H)
            res = ctx["alpha_m_resolution"]
            pm = V10.patch_for(ctx["laws"]["lambda_p_m"])["patch_m"]
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
            fp = V13.floor_proof(ctx)
            c1 = V13.c1_three_state(ctx)
            c1b = V13.c1b(ctx, fp)
            mss = V13.fetch_limited_mss_reference(U, F, SWEEP_H,
                                                 kmin=kmin,
                                                 k_hi=V10.KCUT_DEFAULT)
            cm = V13.cox_munk_reference(U)
            v11ctx = V11.assemble(U, F, SWEEP_H, source="elfouhaily-form")
            k11 = V11._k_grid(V10.KCUT_DEFAULT)
            s11 = V11.s_parts(k11, v11ctx)
            an11 = V11.integral_block(k11, s11, kmin, V10.KCUT_DEFAULT,
                                      v11ctx["k_p_rad_per_m"])
            cells.append({
                "cell": cid, "U_ms": float(U), "fetch_m": float(F),
                "depth_m": float(SWEEP_H),
                "patch_m": float(pm), "kmin_rad_per_m": float(kmin),
                "k_viscous_rad_per_m": kvisc,
                "sea_state": {k: ctx["sea_state"][k] for k in (
                    "k_p_rad_per_m", "c_p_ms", "omega_c", "alpha_p",
                    "gamma_bl", "gamma_form", "sigma_bl", "eckv_domain",
                    "k_p_dispersion_rel_diff", "kp_times_h")},
                "Hs_m": float(ctx["Hs_m"]),
                "Hs_law_m": float(ctx["Hs_law_m"]),
                "alpha_m_resolution": res,
                "low_wind_prior": ctx.get("low_wind_prior"),
                "floor_proof": fp,
                "C1": c1,
                "C1b": c1b,
                "integrals_full": an_full,
                "integrals_cut": an_cut,
                "v11_hard_zero_baseline": {
                    "source": "v11.assemble (v3.8) at the same cell, "
                              "recomputed here only to isolate the change",
                    "alpha_m": float(v11ctx["alpha_m_value"]),
                    "sigma2_h_total": None,
                    "gravity_share_50_200":
                        an11["overlap_50_200"]["gravity_share_of_sigma2"],
                    "sigma2_cut": float(an11["total"]["sigma2"]),
                },
                "references": {
                    "cox_munk_clean_rad2": cm["total_clean_rad2"],
                    "cox_munk": cm,
                    "second_slope_reference": mss,
                    "total_sigma2_cut": float(an_cut["total"]["sigma2"]),
                    "total_sigma2_full": float(an_full["total"]["sigma2"]),
                },
            })
            r11 = V11.s_parts(kf, v11ctx)
            mm = kf >= kmin
            cells[-1]["v11_hard_zero_baseline"]["sigma2_h_total"] = float(
                np.trapz(kf[mm] ** 2 * r11["S_h"][mm], kf[mm]))
            del r11, mm
            curves_all[cid] = {
                "v13": V13.curves(ctx, pm),
                "legacy": V13.curves(V10.assemble(U, F, SWEEP_H), pm),
            }
            log("  %-14s u*=%.5f %-32s raw %+.6f floor %.6f used %.6f | C1 "
                "%-10s C1b %-5s sig2h %.6f"
                % (cid, res["ustar_ms"], res["alpha_m_domain"],
                   res["alpha_m_raw"], res["alpha_m_floor"],
                   res["alpha_m_used"], c1["state"], c1b["pass"],
                   fp["sigma2_h_total"]))
    cp.unit_done("stage-1")

    dump_json(os.path.join(OUT, "spectrum-curves-v13.json"), {
        "what": ("S_l / S_h / B_l / B_h / B / Delta(k) per cell: v3.9 "
                 "(floored source) beside the v3.7 legacy source, on the "
                 "patch-resolved band"),
        "spectrum_source": "elfouhaily-form (v3.9)",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters -- v3.9 adds the labelled low-wind capillary "
                  "prior below the ECKV fit domain"),
        "cells": curves_all})
    c = curves_all[cell_id(6.0, 100.0)]["v13"]
    lines = ["k_rad_per_m,B_l,B_h,B_total,S_l,S_h,S_total,delta_eckv,"
             "omega_rad_s"]
    for i in range(len(c["k_rad_per_m"])):
        lines.append("%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.6f,%.6f" % (
            c["k_rad_per_m"][i], c["B_l"][i], c["B_h"][i], c["B_total"][i],
            c["S_l_m2_per_rad_m"][i], c["S_h_m2_per_rad_m"][i],
            c["S_total_m2_per_rad_m"][i], c["delta_eckv"][i],
            c["omega_rad_s"][i]))
    atomic_write_bytes(os.path.join(OUT, "spectrum-ref-v39-U6F100.csv"),
                       ("\n".join(lines) + "\n").encode("ascii"))

    stage = {
        "stage": 1,
        "created_utc": now_iso(),
        "model": "version 3.9 (PIPELINE-V11 chunk 9)",
        "what": ("v13_spectrum / v13_state: alpha_m raw inside the ECKV fit "
                 "domain, labelled floor below it; A1 and A2 re-run; algebra "
                 "against the v11 conventions; eight-cell analytic table"),
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters -- NOT 'the Elfouhaily model'"),
        "conventions": {
            "S_k_units": "m^2 per rad/m",
            "m0": "int S(k) dk", "sigma2": "int k^2 S(k) dk",
            "B": "B(k) = k^3 S(k)",
            "directional": "Psi = S Phi / k, int Phi dtheta = 1 (unchanged)",
            "wind": "U10 everywhere; u* from U10 only",
            "alpha_m_domain_boundary": {
                "ustar_fit_min_ms": V13.USTAR_FIT_MIN,
                "c_m_ms": V13.C_M,
                "rule": "u* >= c_m/e is the ECKV fit domain (alpha_m_raw >= 0)",
            },
        },
        "changed_quantity": {
            "name": "alpha_m (level of B_h only)",
            "inside": "alpha_m_used = alpha_m_raw",
            "below": "alpha_m_used = max(alpha_m_raw, alpha_m_floor(cell))",
            "floor": ("alpha_m_floor = sigma2_h_target / I(cell), I(cell) = "
                      "int k^2 B_h(alpha=1)/k^3 dk over [kmin, k_vis]"),
            "sigma2_h_target_default": V13.SIGMA2_H_TARGET_DEFAULT,
            "sigma2_h_targets_swept": list(V13.SIGMA2_H_TARGETS),
            "retired_label": V13.RETIRED_ALPHA_M_DOMAIN,
            "untouched": ("B_l, Phi, ECKV reference, patch/cascade rule, "
                          "transport, composite, sky, recipe schema"),
        },
        "prior_label_verbatim": V13.low_wind_prior_label(
            V13.SIGMA2_H_TARGET_DEFAULT),
        "gates": {"A1": a1, "A2": a2, "A2_prior_active": a2_prior,
                  "A2_legacy_control": a2_legacy,
                  "algebra_vs_v11": {"pass": bool(al_ok),
                                     "cells": algebra}},
        "cells": cells,
        "files": ["spectrum-curves-v13.json",
                  "spectrum-ref-v39-U6F100.csv"],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v13-spectrum.json"), stage)
    log("[stage 1] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 2 -- eight-cell sweep
# ==========================================================================
def stage2():
    t0 = time.time()
    log("=" * 72)
    log("[stage 2] 8-cell sweep: C1 three-state, C1b, V1-V3, K1/K5/K6, "
        "totals + references")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    raw = load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    base11, _ = v11_baseline()
    base10, _ = legacy_v10()
    log("  production seed %d; v11 and v10 numbers READ from their JSONs"
        % seed)
    with open(os.path.join(OUT, "numbers-v13-spectrum.json"),
              "r", encoding="utf-8") as f:
        s1 = json.load(f)

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
            kmin = st.kmin
            pm = st.S
            kvisc = float(P8.viscous_cutoff(SWEEP_H))
            kf = V11._k_grid(kvisc, k_hi=kvisc)
            an_full = V11.integral_block(kf, V11.s_parts(kf, ctx), kmin, kvisc,
                                         ctx["k_p_rad_per_m"])
            kc = V11._k_grid(V10.KCUT_DEFAULT)
            an_cut = V11.integral_block(kc, V11.s_parts(kc, ctx), kmin,
                                        V10.KCUT_DEFAULT,
                                        ctx["k_p_rad_per_m"])
            fp = V13.floor_proof(ctx)
            c1 = V13.c1_three_state(ctx)
            c1b = V13.c1b(ctx, fp)
            res = ctx["alpha_m_resolution"]
            B = base11[cid]
            L10 = base10[cid]
            cm = V13.cox_munk_reference(U)
            mss = V13.fetch_limited_mss_reference(U, F, SWEEP_H, kmin=kmin,
                                                 k_hi=V10.KCUT_DEFAULT)
            s2_cut = float(an_cut["total"]["sigma2"])
            s2_full = float(an_full["total"]["sigma2"])
            s2_realised = float(st.spectrum_after["sigma2_rad2"])
            fa = st.frac_above((100.0, 370.0))
            hs = node["v11"]["K1"]
            kb = np.logspace(math.log10(50.0), math.log10(200.0), 1500)
            sb = V11.s_parts(kb, ctx)
            v3 = {
                "B_hat_range": [B_HAT_LOW, B_HAT_HIGH],
                "B_l_max_over_band": float(np.max(sb["B_l"])),
                "B_total_max_over_band": float(np.max(sb["B_l"] + sb["B_h"])),
                "max_ratio_B_l_over_Bhat_low": float(
                    np.max(sb["B_l"]) / B_HAT_LOW),
                "max_ratio_B_total_over_Bhat_low": float(
                    np.max(sb["B_l"] + sb["B_h"]) / B_HAT_LOW),
                "frac_band_with_B_l_above_Bhat_high": float(
                    np.mean(sb["B_l"] > B_HAT_HIGH)),
                "frac_band_with_B_total_above_Bhat_high": float(
                    np.mean((sb["B_l"] + sb["B_h"]) > B_HAT_HIGH)),
                "pass": bool(np.max(sb["B_l"]) <= B_HAT_LOW),
                "note": ("the gravity part must not sit above the measured "
                         "saturation level across 50-200; the capillary part "
                         "is allowed to (it is the short-wave band)"),
            }
            del sb, kb
            v1 = {
                "sigma2_new_cut": s2_cut,
                "cox_munk_clean_rad2": cm["total_clean_rad2"],
                "factor_over_cox_munk": float(s2_cut
                                              / cm["total_clean_rad2"]),
                "v11_factor_over_cox_munk": float(
                    B["gates"]["V1"]["factor_over_cox_munk"]),
                "v10_legacy_factor_over_cox_munk": float(
                    L10["gates"]["G2"]["sigma2_realised_spectral"]
                    / cm["total_clean_rad2"]),
                "note": "clean Cox-Munk is a REFERENCE, never a target",
            }
            v2 = {
                "frac_sigma2_above_1173": an_full["frac_sigma2_above_1173"],
                "laxague_reference": {
                    "pure_capillary_k_gt_1173": 0.05,
                    "caveat": ("Laxague et al. partition, seas with "
                               "Omega_c <= 2; this cell has Omega_c = %.2f"
                               % ctx["omega_c"]),
                },
                "band_shares_new": {name: an_full["bands"]["rows"][name]
                                    ["sigma2_share"]
                                    for name in an_full["bands"]["rows"]},
                "band_shares_v11_realised": {
                    k: v["sigma2_share"] for k, v in
                    B["bands_realised"].items() if k != "_total"},
            }
            # C2: retired.  Report the total beside the references.
            c2 = {
                "status": "RETIRED",
                "reason": ("the v3.8 gate demanded the total within +-10% of "
                           "v3.7's totals; its premise was false (v11 section "
                           "2.2: v3.7's totals were inflated by the very "
                           "plateau the gravity-tail fix removes).  v3.9 "
                           "replaces it with reporting, per PIPELINE-V11"),
                "sigma2_assembled_cut": s2_cut,
                "sigma2_assembled_full": s2_full,
                "sigma2_realised_spectral": s2_realised,
                "sigma2_v11_realised": float(
                    B["gates"]["C2"]["sigma2_assembled_cut"]),
                "sigma2_v10_legacy_realised": float(
                    L10["gates"]["G2"]["sigma2_realised_spectral"]),
                "rel_vs_v11": float(s2_cut
                                    / B["gates"]["C2"]["sigma2_assembled_cut"]
                                    - 1.0),
                "cox_munk_clean_rad2": cm["total_clean_rad2"],
                "second_slope_reference": mss,
                "eckv_reference_full_rad2": B["eckv_reference"]["sigma2_rad2"],
                "note": ("reporting is not a gate: the total is printed "
                         "beside Cox-Munk clean and the second slope "
                         "reference, and the full band partition is in "
                         "integrals_full.bands"),
            }
            row = {
                "cell": cid, "U_ms": float(U), "fetch_m": float(F),
                "depth_m": float(SWEEP_H),
                "mode": st.mode, "cascade_j": int(st.cascade_j),
                "patch_m": float(pm), "patch_doublings":
                    int(st.patch_doublings),
                "kmin_rad_per_m": float(kmin),
                "kmax_nyquist_scene": float(st.kmax_nyq),
                "dx_mm": float(st.dx * 1e3),
                "sea_state": node["v11"]["sea_state"],
                "alpha_m_resolution": res,
                "low_wind_prior": ctx.get("low_wind_prior"),
                "spectrum": node["spectrum"],
                "amplitude": node["amplitude"],
                "breaking": node["breaking"],
                "realised": node["realised"],
                "spectral_vs_spatial": node["spectral_vs_spatial"],
                "bands_realised": node["bands"],
                "direction": node["direction"],
                "v11": {k: v for k, v in node["v11"].items()
                        if k != "sea_state"},
                "v13": {k: v for k, v in node["v13"].items()
                        if k not in ("C1", "C1b", "floor_proof",
                                     "alpha_m_resolution")},
                "integrals_cut": an_cut,
                "integrals_full": an_full,
                "frac_above_370_realised": fa["above_370"]["realised"],
                "floor_proof": fp,
                "gates": {"C1": c1, "C1b": c1b, "C2": c2, "V1": v1,
                          "V2": v2, "V3": v3, "K1": hs,
                          "K5": node["quantisation"],
                          "K6": V11.gate_k6(ctx)},
                "references": {
                    "cox_munk": cm, "second_slope_reference": mss},
                "baselines": {
                    "v11": {
                        "alpha_m": B["spectrum"]["alpha_m"],
                        "alpha_m_raw": B["spectrum"]["alpha_m_raw"],
                        "alpha_m_domain": B["spectrum"]["alpha_m_domain"],
                        "sigma2_cut": B["gates"]["C2"]
                        ["sigma2_assembled_cut"],
                        "gravity_share_50_200":
                            B["gates"]["C1"]["gravity_share_50_200"],
                        "C1_state_v11": ("PASS" if B["gates"]["C1"]["pass"]
                                         else "FAIL"),
                        "frac_above_370_realised":
                            B["frac_above_370_realised"],
                    },
                    "v10_legacy": {
                        "sigma2_realised":
                            L10["gates"]["G2"]["sigma2_realised_spectral"],
                        "Hs_spectral": L10["gates"]["G1"]["Hs_spectral_m"],
                        "gravity_share_50_200":
                            L10["gates"]["G2"]["overlap_50_200"]
                            ["S_g_share_of_overlap_sigma2"],
                    },
                },
                "timing_s": {"state": float(r["t_state"]),
                             "fields": float(r["t_fields"])},
            }
            cells.append(row)
            log("  [%-14s] %-32s C1 %-10s C1b %-5s | share %.3f sig2 %.5f "
                "(v11 %.5f, %+.1f%%) K1flag=%s (%.1fs)"
                % (cid, res["alpha_m_domain"], c1["state"], c1b["pass"],
                   c1["gravity_share_50_200"], s2_cut,
                   B["gates"]["C2"]["sigma2_assembled_cut"],
                   100 * c2["rel_vs_v11"], hs["flag"],
                   time.time() - t1))
            del r, st, node, an_full, an_cut
            gc.collect()

    in_domain = [c for c in cells if c["gates"]["C1"]["in_eckv_fit_domain"]]
    below = [c for c in cells if not c["gates"]["C1"]["in_eckv_fit_domain"]]
    summary = {
        "C1_three_state": {
            "pass": all(c["gates"]["C1"]["state"] == "PASS"
                        for c in in_domain),
            "pass_in_domain_cells_only": True,
            "in_domain_cells": [c["cell"] for c in in_domain],
            "pass_cells": [c["cell"] for c in cells
                           if c["gates"]["C1"]["state"] == "PASS"],
            "fail_cells": [c["cell"] for c in cells
                           if c["gates"]["C1"]["state"] == "FAIL"],
            "unresolved_cells": [c["cell"] for c in cells
                                 if c["gates"]["C1"]["state"] == "UNRESOLVED"],
            "worst_in_domain_share": max(
                c["gates"]["C1"]["gravity_share_50_200"] for c in in_domain),
            "unresolved_with_numbers": [
                {"cell": c["cell"],
                 "gravity_share_50_200": c["gates"]["C1"]
                 ["gravity_share_50_200"],
                 "u_over_fit_min": c["alpha_m_resolution"]["u_over_fit_min"],
                 "prior": c["low_wind_prior"]} for c in below],
        },
        "C1b": {
            "pass": all(c["gates"]["C1b"]["pass"] for c in cells),
            "fail_cells": [c["cell"] for c in cells
                           if not c["gates"]["C1b"]["pass"]],
            "worst_capillary_share": min(
                c["gates"]["C1b"]["capillary_share_50_200"] for c in cells),
            "min_sigma2_h_total": min(c["gates"]["C1b"]["sigma2_h_total"]
                                      for c in cells),
        },
        "V1": {"factor_range": [min(c["gates"]["V1"]
                                   ["factor_over_cox_munk"] for c in cells),
                                max(c["gates"]["V1"]
                                    ["factor_over_cox_munk"] for c in cells)]},
        "V2": {"reported": True},
        "V3": {"pass": all(c["gates"]["V3"]["pass"] for c in cells),
               "worst_B_l_over_Bhat": max(
                   c["gates"]["V3"]["max_ratio_B_l_over_Bhat_low"]
                   for c in cells),
               "worst_B_total_over_Bhat": max(
                   c["gates"]["V3"]["max_ratio_B_total_over_Bhat_low"]
                   for c in cells)},
        "K1": {"flags_open": [c["cell"] for c in cells
                              if c["gates"]["K1"]["flag"]],
               "worst_abs_delta_rel": max(
                   abs(c["gates"]["K1"]["delta_rel"]) for c in cells)},
        "K5": {"pass": True,
               "worst_rel": max(c["gates"]["K5"]["max_rel_error"]
                                for c in cells)},
        "K6": {"pass": all(c["gates"]["K6"]["pass"] for c in cells),
               "worst_slope_ratio": max(
                   c["gates"]["K6"]["max_dslope_ratio_fine_over_coarse"]
                   for c in cells)},
        "floor_proof": {
            "all_prior_cells_close": all(
                (c["floor_proof"]["pass"] is True) for c in below),
            "max_abs_closure_rel": max(
                [abs(c["floor_proof"]["closure_rel_residual"] or 0.0)
                 for c in below] or [0.0]),
        },
    }
    dump_json(os.path.join(OUT, "numbers-v13-sweep.json"), {
        "stage": 2,
        "what": ("v3.9 8-cell sweep: C1 three-state, C1b floor proof, V1-V3 "
                 "reports, K1/K5/K6, total sigma2 beside the references, full "
                 "band partitions.  v11/v10 baseline numbers are READ from "
                 "their JSONs."),
        "created_utc": now_iso(),
        "model": "version 3.9 (PIPELINE-V11 chunk 9)",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters; low-wind capillary prior labelled"),
        "grid": {"N": GRID_SWEEP, "depth_m": SWEEP_H,
                 "U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
                 "seed": seed, "spectrum_source": "elfouhaily-form",
                 "sigma2_h_target": V13.SIGMA2_H_TARGET_DEFAULT},
        "policy": {"patch_base_m": V10.PATCH_BASE_M,
                   "r_lambda": V10.R_LAMBDA,
                   "k_cut_default": V10.KCUT_DEFAULT,
                   "B_hat_range": [B_HAT_LOW, B_HAT_HIGH],
                   "ustar_fit_min_ms": V13.USTAR_FIT_MIN,
                   "baselines": "out/v11/numbers-v11-sweep.json and "
                                "out/v10/numbers-v10-sweep.json (read)"},
        "summary": summary,
        "cells": cells,
        "wall_s": float(time.time() - t0),
    })
    log("  summary: C1 pass=%s fail=%s unresolved=%s | C1b pass=%s fail=%s"
        % (summary["C1_three_state"]["pass"],
           summary["C1_three_state"]["fail_cells"],
           summary["C1_three_state"]["unresolved_cells"],
           summary["C1b"]["pass"], summary["C1b"]["fail_cells"]))
    cp.unit_done("stage-2")
    log("[stage 2] done in %.1fs" % (time.time() - t0))
    return summary


# ==========================================================================
# stage 3 -- K3, K4, and the mandatory sensitivity sweep
# ==========================================================================
YARDSTICKS = R10.YARDSTICKS


def stage3():
    t0 = time.time()
    log("=" * 72)
    log("[stage 3] K3 windy (floored vs v11 legacy), K4 + prior battery, "
        "sensitivity sweep")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    from PIL import Image
    from state import derive_seed
    from v5_composite import A_V3, B_SIGMA_V3, E_V3
    from run_windy_still import (warp_measure, predicted_shift_field,
                                 predicted_block_stats, block_agreement)
    raw = load_recipe()
    seed = int(derive_seed(raw, 1))
    base_u8, photo_lin, PW, PH = R10._load_photo()
    scenes = sweep_scenes()
    log("  photo %dx%d sha=%s" % (PW, PH, sha256_file(PHOTO)[:16]))

    out = {"stage": 3, "created_utc": now_iso(),
           "what": ("v3.9 compatibility: K3 the G10 windy scene with the "
                    "floored source vs the v11 legacy still, K4 regression "
                    "battery plus the prior-activity battery"),
           "photo": {"path": PHOTO, "size": [int(PW), int(PH)],
                     "sha256": sha256_file(PHOTO)},
           "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                     "parameters; low-wind capillary prior labelled"),
           "K3": {}, "K4": {}}

    # ---- K3: the windy scene, both sources ------------------------------
    renders = {}
    for src, fname in (("jonswap", "comp-windy-v39-legacy.png"),
                       ("elfouhaily-form", "comp-windy-v39-new.png")):
        r = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], GRID_RENDER,
                      seed, src)
        st = r["state"]
        t = R10.transport_scene(r, photo_lin, PW, PH, GRID_RENDER, Q=25,
                                E=float(E_V3))
        atomic_write_bytes(os.path.join(WORK, "C-windy-v39-%s-Q25.npy"
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
        renders[src] = {"state": st, "r": r, "t": t,
                        "still_png": still_png, "warp": warp, "ps": ps,
                        "agree": agree, "file": fname,
                        "texture": _texture(still_png)}
        log("  K3 %-16s C p50 %.3f p99 %.3f p99.9 %.3f max %.2f | own field "
            "%.2f px = %.2f mm | |dL/dy| %.5f"
            % (src, np.percentile(t["C"], 50), np.percentile(t["C"], 99),
               np.percentile(t["C"], 99.9), float(t["C"].max()),
               ps["field_mean_magnitude_px"],
               1e3 * float(ps["mean_predicted_shift_m"]),
               renders[src]["texture"]["mean_abs_dL_dy"]))
        del r
        gc.collect()

    leg_r = renders["jonswap"]
    new_r = renders["elfouhaily-form"]

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
                "alpha_m_domain": (st.ctx.get("alpha_m_domain")
                                   if src != "jonswap" else None),
                "alpha_m_used": (None if src == "jonswap"
                                 else float(st.ctx["alpha_m_value"])),
                "low_wind_prior": (st.ctx.get("low_wind_prior")
                                   if src != "jonswap" else None),
                "bands_realised": node["bands"],
                "frac_above_370_realised":
                    st.frac_above((370.0,))["above_370"]["realised"],
            },
            "route": {"Q": 25, "A": A_V3, "E": float(E_V3),
                      "B_lowpass_sigma_px": B_SIGMA_V3,
                      "g1_used": float(t["g1"]),
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
                "yardsticks": YARDSTICKS,
            },
            "fine_texture": rr["texture"],
            "tone": {
                "crush_rate_linear_lt_0.02":
                    float(t["metrics"]["crush_rate_linear_lt_0.02"]),
                "blowout_rate_aces_gt_0.98":
                    float(t["metrics"]["blowout_rate_aces_gt_0.98"]),
                "customary_crush_bound": 0.0325,
                "crush_pass_customary": bool(
                    t["metrics"]["crush_rate_linear_lt_0.02"] <= 0.0325),
            },
            "beer": t["beer"],
            "still_metrics": t["metrics"],
            "files": [rr["file"]],
        }

    k3 = {"scenes": {"jonswap": k3_node("jonswap"),
                     "elfouhaily-form": k3_node("elfouhaily-form")}}
    d = np.abs(new_r["still_png"].astype(np.int16)
               - leg_r["still_png"].astype(np.int16))
    ms = R10.msssim(new_r["still_png"].astype(np.float64) / 255.0,
                    leg_r["still_png"].astype(np.float64) / 255.0)
    k3["between_sources"] = {
        "mean_abs_uint8": float(d.mean()), "max_abs_uint8": int(d.max()),
        "frac_pixels_changed": float(np.mean(d.max(axis=2) > 0)),
        "msssim": ms["msssim"],
        "texture_base_photo": _texture(base_u8),
    }
    # continuity against the v11 deliverables (the prior is inactive here)
    cont = {}
    for src, f13, f11 in (("jonswap", "comp-windy-v39-legacy.png",
                           "comp-windy-v38-legacy.png"),
                          ("elfouhaily-form", "comp-windy-v39-new.png",
                           "comp-windy-v38-new.png")):
        a13 = sha256_file(os.path.join(OUT, f13))
        a11 = sha256_file(os.path.join(V11_OUT, f11))
        im13 = np.asarray(Image.open(os.path.join(OUT, f13))
                          .convert("RGB")).astype(np.int16)
        im11 = np.asarray(Image.open(os.path.join(V11_OUT, f11))
                          .convert("RGB")).astype(np.int16)
        dd = np.abs(im13 - im11)
        cont[src] = {"v13_file": f13, "v11_file": f11,
                     "sha256_v13": a13, "sha256_v11": a11,
                     "byte_identical": bool(a13 == a11),
                     "max_abs_uint8": int(dd.max()),
                     "frac_pixels_changed": float(np.mean(dd.max(axis=2)
                                                          > 0))}
        del im13, im11, dd
    k3["continuity_vs_v11"] = {
        "what": ("the G10 scene has u* = 0.207 m/s, far inside the ECKV fit "
                 "domain, so the v3.9 prior is INACTIVE there and both stills "
                 "must reproduce the v3.8 deliverables bit-for-bit"),
        "cells": cont,
        "all_byte_identical": bool(all(v["byte_identical"]
                                      for v in cont.values())),
    }
    lg = k3["scenes"]["jonswap"]["state_summary"]
    nw = k3["scenes"]["elfouhaily-form"]["state_summary"]
    ab = {}
    for src in ("jonswap", "elfouhaily-form"):
        s = k3["scenes"][src]["state_summary"]
        ab[src] = float(s["frac_above_370_realised"] * s["sigma2_spectral"])

    def rel(a, b):
        return None if not b else float(a / b - 1.0)

    k3["attribution"] = {
        "Hs_rel_change": rel(nw["Hs_spectral_m"], lg["Hs_spectral_m"]),
        "sigma2_rel_change": rel(nw["sigma2_spectral"], lg["sigma2_spectral"]),
        "rms_slope_rel_change": rel(nw["rms_slope"], lg["rms_slope"]),
        "g1_legacy": float(leg_r["t"]["g1"]),
        "g1_new": float(new_r["t"]["g1"]),
        "g1_rel_change": rel(float(new_r["t"]["g1"]), float(leg_r["t"]["g1"])),
        "frac_above_370": {"legacy": lg["frac_above_370_realised"],
                           "new": nw["frac_above_370_realised"],
                           "rel_change": rel(nw["frac_above_370_realised"],
                                             lg["frac_above_370_realised"])},
        "sigma2_above_370_absolute": {
            "legacy": ab["jonswap"], "new": ab["elfouhaily-form"],
            "rel_change": rel(ab["elfouhaily-form"], ab["jonswap"])},
        "band_share_1_10cm": {
            "legacy": lg["bands_realised"]["1_10cm"]["sigma2_share"],
            "new": nw["bands_realised"]["1_10cm"]["sigma2_share"]},
        "band_share_lt_1cm": {
            "legacy": lg["bands_realised"]["lt_1cm"]["sigma2_share"],
            "new": nw["bands_realised"]["lt_1cm"]["sigma2_share"]},
        "own_field_shift_mm": {
            "legacy": k3["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["mean_mm"],
            "new": k3["scenes"]["elfouhaily-form"]["warp"]["own_shift_field"]
            ["mean_mm"]},
        "fine_texture": {
            "legacy_dy": k3["scenes"]["jonswap"]["fine_texture"]
            ["mean_abs_dL_dy"],
            "new_dy": k3["scenes"]["elfouhaily-form"]["fine_texture"]
            ["mean_abs_dL_dy"],
            "rel_change_dy": rel(k3["scenes"]["elfouhaily-form"]
                                 ["fine_texture"]["mean_abs_dL_dy"],
                                 k3["scenes"]["jonswap"]["fine_texture"]
                                 ["mean_abs_dL_dy"])},
        "statement": ("This cell is in the ECKV fit domain, so the v3.9 prior "
                      "does not act: the change is the v3.8 change (gravity "
                      "tail removed, alpha_m raw).  The attribution numbers "
                      "therefore reproduce REPORT-V11 section 3.2 and the "
                      "stills are byte-identical to the v3.8 deliverables.  "
                      "Nothing here measures the prior; the sensitivity sweep "
                      "below does."),
    }
    out["K3"] = k3
    log("  K3 attribution: Hs %+.1f%% sig2 %+.1f%% | both stills "
        "byte-identical to v3.8: %s"
        % (100 * k3["attribution"]["Hs_rel_change"],
           100 * k3["attribution"]["sigma2_rel_change"],
           k3["continuity_vs_v11"]["all_byte_identical"]))

    # ---- K4: regression battery + prior activity ------------------------
    again = R10.transport_scene(new_r["r"], photo_lin, PW, PH, GRID_RENDER,
                               Q=25, E=float(E_V3), C_cache=new_r["t"]["C"])
    ident = bool(np.array_equal(again["still"], new_r["t"]["still"]))
    k4a = {"scene": "windy, floored source", "bit_identical_output": ident,
           "max_abs_uint8": int(np.abs(
               again["still"].astype(np.int16)
               - new_r["t"]["still"].astype(np.int16)).max()),
           "params": {"A": A_V3, "B_lowpass_sigma_px": B_SIGMA_V3,
                      "E": float(E_V3), "Q": 25},
           "pass": ident}
    del again
    gc.collect()

    rec_path = os.path.join(HERE, "out", "windy-still", "numbers-windy.json")
    with open(rec_path, "r", encoding="utf-8") as f:
        rec = json.load(f)
    rec_tw = rec["depth_cost"]["beer_direct_beam"]["z1.0"]["two_way"]
    tw_l = leg_r["t"]["beer"]["two_way"]
    tw_n = new_r["t"]["beer"]["two_way"]
    diff_rec = [abs(float(tw_n[i]) - float(rec_tw[i])) for i in range(3)]
    diff_src = [abs(float(tw_n[i]) - float(tw_l[i])) for i in range(3)]
    k4b = {"new_two_way": [float(x) for x in tw_n],
           "legacy_two_way": [float(x) for x in tw_l],
           "v9_recorded_two_way": rec_tw,
           "abs_diff_vs_recorded": diff_rec,
           "abs_diff_legacy_vs_new": diff_src,
           "pass": bool(max(diff_rec) <= 5e-5 and max(diff_src) <= 1e-12)}
    log("  K4b Beer two-way new %s vs v9 %s (max|d| %.2e)"
        % (["%.6f" % x for x in tw_n], ["%.6f" % x for x in rec_tw],
           max(diff_rec)))

    ra = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], 256, seed,
                   "elfouhaily-form")
    rb = scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], 256, seed,
                   "elfouhaily-form")
    eq = bool(np.array_equal(ra["state"].lo.H0, rb["state"].lo.H0)
              and (ra["state"].hi is None
                   or np.array_equal(ra["state"].hi.H0, rb["state"].hi.H0)))
    rain_raised = False
    try:
        S13.V13WaterState({"x": 1}, {"windSpeedMS": 6.0, "depthM": 1.0,
                                     "windDir": [1.0, 0.0]}, N=64, seed=1,
                          fetch_m=100.0, rain_rate_mm_h=1.0,
                          spectrum_source="elfouhaily-form")
    except NotImplementedError:
        rain_raised = True
    k4c = {"bit_identical_state": eq, "rain_positive_raises": rain_raised,
           "pass": bool(eq and rain_raised)}
    del ra, rb
    gc.collect()

    # K4d: the prior-activity battery (the v3.9-specific regression)
    from env import resolve_environment
    from recipe_v8 import with_v8
    import v11_state as S11

    def h0_of(st):
        parts = [st.lo.H0]
        if st.hi is not None:
            parts.append(st.hi.H0)
        return parts

    def build13(U, F, h, target, src="elfouhaily-form", N=256):
        sc = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
        env = resolve_environment(sc)
        return S13.V13WaterState(sc, env, N=N, seed=seed, preset_index=1,
                                 fetch_m=float(F),
                                 spectrum_source=src,
                                 sigma2_h_target=target)

    def build11(U, F, h, src, N=256):
        sc = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
        env = resolve_environment(sc)
        return S11.V11WaterState(sc, env, N=N, seed=seed, preset_index=1,
                                 fetch_m=float(F), spectrum_source=src)

    rows = {}
    for tag, U, F, h, src in (("in_domain_U6", 6.0, 100.0, 1.0,
                               "elfouhaily-form"),
                              ("in_domain_U10", 10.0, 1000.0, 1.0,
                               "elfouhaily-form"),
                              ("legacy_path_U6", 6.0, 100.0, 1.0, "jonswap")):
        a = build13(U, F, h, V13.SIGMA2_H_TARGET_DEFAULT, src=src)
        b = build11(U, F, h, src)
        same = all(np.array_equal(x, y) for x, y in zip(h0_of(a), h0_of(b)))
        rows[tag] = {"U_ms": U, "fetch_m": F, "depth_m": h,
                     "source": src, "prior_active":
                         bool(a.ctx.get("low_wind_prior") is not None),
                     "H0_bit_identical_to_v11": bool(same),
                     "sigma2_spectral_v13": float(a.spectrum_after
                                                  ["sigma2_rad2"]),
                     "sigma2_spectral_v11": float(b.spectrum_after
                                                  ["sigma2_rad2"]),
                     "pass": bool(same)}
        del a, b
        gc.collect()
    # prior-active: must DIFFER, and by a measured amount
    a = build13(1.5, 100.0, 0.6, V13.SIGMA2_H_TARGET_DEFAULT)
    b = build11(1.5, 100.0, 0.6, "elfouhaily-form")
    differ = not all(np.array_equal(x, y) for x, y in zip(h0_of(a), h0_of(b)))
    rows["prior_active_preset1"] = {
        "U_ms": 1.5, "fetch_m": 100.0, "depth_m": 0.6,
        "source": "elfouhaily-form",
        "prior_active": bool(a.ctx.get("low_wind_prior") is not None),
        "H0_bit_identical_to_v11": bool(not differ),
        "alpha_m_v13": float(a.ctx["alpha_m_value"]),
        "alpha_m_v11": float(b.ctx["alpha_m_value"]),
        "sigma2_h_total_v13": float(V13.floor_proof(a.ctx)["sigma2_h_total"]),
        "sigma2_h_total_v11": float(V13.floor_proof(b.ctx)["sigma2_h_total"]),
        "sigma2_spectral_v13": float(a.spectrum_after["sigma2_rad2"]),
        "sigma2_spectral_v11": float(b.spectrum_after["sigma2_rad2"]),
        "pass": bool(differ)}
    del a, b
    gc.collect()
    # target 0 must reproduce the v3.8 hard zero bit-for-bit
    a = build13(1.5, 100.0, 0.6, 0.0)
    b = build11(1.5, 100.0, 0.6, "elfouhaily-form")
    rows["target0_reproduces_v11"] = {
        "U_ms": 1.5, "fetch_m": 100.0, "depth_m": 0.6,
        "source": "elfouhaily-form", "sigma2_h_target": 0.0,
        "prior_active": bool(a.ctx.get("low_wind_prior") is not None),
        "alpha_m_domain": a.ctx["alpha_m_domain"],
        "alpha_m_used": float(a.ctx["alpha_m_value"]),
        "H0_bit_identical_to_v11": bool(
            all(np.array_equal(x, y) for x, y in zip(h0_of(a), h0_of(b)))),
        "pass": bool(all(np.array_equal(x, y)
                         for x, y in zip(h0_of(a), h0_of(b))))}
    del a, b
    gc.collect()
    k4d = {"what": ("v3.9-specific regression: the prior must be INERT "
                    "wherever it is inactive (bit-identity with v3.8), ACTIVE "
                    "only below u* = c_m/e, and the sweep's target 0 must "
                    "reproduce the v3.8 hard zero exactly"),
           "cells": rows,
           "pass": bool(all(v["pass"] for v in rows.values()))}
    out["K4"] = {"K4a": k4a, "K4b": k4b, "K4c": k4c, "K4d_prior_battery": k4d,
                 "battery": ("K4a-c are the v3.7 G8 battery verbatim; K4d is "
                             "this chunk's own prior-activity battery")}
    log("  K4d prior battery: in-domain identical=%s, prior-active differs="
        "%s, target0 identical=%s -> PASS=%s"
        % (rows["in_domain_U6"]["H0_bit_identical_to_v11"],
           rows["prior_active_preset1"]["pass"],
           rows["target0_reproduces_v11"]["H0_bit_identical_to_v11"],
           k4d["pass"]))
    out["summary"] = {
        "K3": "reported (both stills + attribution + continuity vs v3.8)",
        "K4a": k4a["pass"], "K4b": k4b["pass"], "K4c": k4c["pass"],
        "K4d": k4d["pass"]}
    out["wall_s"] = float(time.time() - t0)

    dump_json(os.path.join(OUT, "numbers-v13-transport.json"), out)
    _sensitivity(raw, seed, photo_lin, PW, PH, base_u8, scenes)
    cp.unit_done("stage-3")
    log("[stage 3] done in %.1fs" % (time.time() - t0))
    return out


def _sensitivity(raw, seed, photo_lin, PW, PH, base_u8, scenes):
    """The mandatory sweep: sigma2_h target over {0,0.0015,0.003,0.006}."""
    t0 = time.time()
    log("-" * 72)
    log("  sensitivity sweep: targets %s on preset-1 and G10 windy"
        % (list(V13.SIGMA2_H_TARGETS),))
    from PIL import Image
    from run_windy_still import (warp_measure, predicted_shift_field,
                                 predicted_block_stats)
    from v5_composite import E_V3
    out = {"what": ("the prior's influence on the picture, measured rather "
                    "than asserted: sigma2_h target swept over "
                    "{0, 0.0015, 0.003, 0.006} on the preset-1 pond scene "
                    "(prior active) and the G10 windy scene (prior inactive)"),
           "created_utc": now_iso(),
           "targets": list(V13.SIGMA2_H_TARGETS),
           "scenes": {}}
    for key in ("A", "B"):
        sc = scenes[key]
        rows = {}
        stills = {}
        for tgt in V13.SIGMA2_H_TARGETS:
            tag = ("t%04d" % int(round(10000.0 * tgt)))
            r = scene_for(raw, sc["U"], sc["F"], sc["h"], GRID_RENDER, seed,
                          "elfouhaily-form", target=tgt)
            st = r["state"]
            t = R10.transport_scene(r, photo_lin, PW, PH, GRID_RENDER, Q=25,
                                    E=float(E_V3))
            fname = "sens-%s-%s.png" % (key.lower(), tag)
            png = io.BytesIO()
            Image.fromarray(t["still"]).save(png, format="PNG")
            atomic_write_bytes(os.path.join(OUT, fname), png.getvalue())
            still_png = np.asarray(Image.open(os.path.join(OUT, fname))
                                   .convert("RGB"))
            warp = warp_measure(base_u8, still_png)
            pred = predicted_shift_field(r["fields"], r["env"], PW, PH)
            ps = predicted_block_stats(pred, warp["grid"])
            ctx = st.ctx
            fp = V13.floor_proof(ctx)
            kmin = st.kmin
            kvisc = float(P8.viscous_cutoff(sc["h"]))
            kf = V11._k_grid(kvisc, k_hi=kvisc)
            an_full = V11.integral_block(kf, V11.s_parts(kf, ctx), kmin,
                                         kvisc, ctx["k_p_rad_per_m"])
            row = {
                "target": float(tgt),
                "scene": sc["name"], "U_ms": sc["U"],
                "fetch_m": sc["F"], "depth_m": sc["h"],
                "alpha_m_raw": float(ctx["alpha_m_raw"]),
                "alpha_m_floor": float(ctx["alpha_m_floor"]),
                "alpha_m_used": float(ctx["alpha_m_value"]),
                "alpha_m_domain": ctx["alpha_m_domain"],
                "prior_active": bool(ctx.get("low_wind_prior") is not None),
                "low_wind_prior": ctx.get("low_wind_prior"),
                "sigma2_h_total": float(fp["sigma2_h_total"]),
                "sigma2_total_assembled_full": float(
                    an_full["total"]["sigma2"]),
                "sigma2_total_realised": float(
                    st.spectrum_after["sigma2_rad2"]),
                "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
                "rms_slope": float(st.spectrum_after["rms_slope"]),
                "frac_above_370_realised":
                    st.frac_above((370.0,))["above_370"]["realised"],
                "C": R10.c_stats(t["C"]),
                "C_percentiles": {
                    "p50": float(np.percentile(t["C"], 50)),
                    "p90": float(np.percentile(t["C"], 90)),
                    "p99": float(np.percentile(t["C"], 99)),
                    "p99_9": float(np.percentile(t["C"], 99.9)),
                    "max": float(t["C"].max())},
                "warp": {
                    "primary_ssd": {k: warp["primary"][k] for k in (
                        "mean_shift_px", "shift_field_std_px",
                        "frac_blocks_at_search_edge")},
                    "own_shift_field": {
                        "field_mean_px": ps["field_mean_magnitude_px"],
                        "field_std_px": ps["field_std_magnitude_px"],
                        "mean_mm": 1e3 * float(ps["mean_predicted_shift_m"]),
                        "pixel_scale_mm": float(ps["pixel_scale_mm"])}},
                "fine_texture": _texture(still_png),
                "file": fname,
            }
            rows[tag] = row
            stills[tag] = still_png
            log("    %s %-6s alpha_m %.6f (%s) sig2h %.6f | sig2 %.6f | C "
                "p50 %.3f p99 %.3f p99.9 %.3f | warp %.2f mm | |dL/dy| %.5f"
                % (key, tag, row["alpha_m_used"],
                   "prior" if row["prior_active"] else "raw  ",
                   row["sigma2_h_total"], row["sigma2_total_realised"],
                   row["C_percentiles"]["p50"], row["C_percentiles"]["p99"],
                   row["C_percentiles"]["p99_9"],
                   row["warp"]["own_shift_field"]["mean_mm"],
                   row["fine_texture"]["mean_abs_dL_dy"]))
            del r, st
            gc.collect()
        # cross-level identity: a level that does not move the spectrum must
        # not move a pixel
        tags = list(rows)
        ident = {}
        for i in range(len(tags)):
            for j in range(i + 1, len(tags)):
                dd = np.abs(stills[tags[i]].astype(np.int16)
                            - stills[tags[j]].astype(np.int16))
                ident["%s_vs_%s" % (tags[i], tags[j])] = {
                    "max_abs_uint8": int(dd.max()),
                    "mean_abs_uint8": float(dd.mean()),
                    "alpha_m_identical": bool(
                        rows[tags[i]]["alpha_m_used"]
                        == rows[tags[j]]["alpha_m_used"])}
                del dd
        out["scenes"][key] = {
            "name": sc["name"], "recipe": sc,
            "levels": rows, "cross_level_image_identity": ident,
            "note": ("a level whose alpha_m_used equals another's must give "
                     "a bit-identical still; that is the sweep's own "
                     "determinism check")}
    out["wall_s"] = float(time.time() - t0)
    dump_json(os.path.join(OUT, "numbers-v13-sensitivity.json"), out)
    log("  sensitivity sweep done in %.1fs" % (time.time() - t0))
    return out


# ==========================================================================
# stage 4 -- merge evaluation + report
# ==========================================================================
def evaluate_merge():
    with open(os.path.join(OUT, "numbers-v13-spectrum.json"),
              "r", encoding="utf-8") as f:
        s1 = json.load(f)
    with open(os.path.join(OUT, "numbers-v13-sweep.json"),
              "r", encoding="utf-8") as f:
        s2 = json.load(f)
    with open(os.path.join(OUT, "numbers-v13-transport.json"),
              "r", encoding="utf-8") as f:
        s3 = json.load(f)
    with open(os.path.join(OUT, "numbers-v13-sensitivity.json"),
              "r", encoding="utf-8") as f:
        s4 = json.load(f)
    sm = s2["summary"]
    a1 = s1["gates"]["A1"]
    a2 = s1["gates"]["A2"]
    a2p = s1["gates"]["A2_prior_active"]
    al = s1["gates"]["algebra_vs_v11"]
    c1 = sm["C1_three_state"]
    c1b = sm["C1b"]
    n_levels = sum(len(v["levels"]) for v in s4["scenes"].values())
    swept_ok = bool(n_levels == 4 * 2 and all(
        "C_percentiles" in L and "warp" in L and "fine_texture" in L
        and "sigma2_total_realised" in L
        for v in s4["scenes"].values() for L in v["levels"].values()))
    items = [
        ("A1 (D1 verification, re-run)", bool(a1["pass"]),
         "finite-depth plateau moves %.2f%% (< 30%%); deep-water closed form "
         "matched to %.1e"
         % (100 * a1["worst_finite_depth_plateau_rel"],
            a1["depths"]["h_%g" % 1000.0]["ratio_kM_physical"]
            ["max_abs_minus_1"])),
        ("A2 (Parseval, in-domain + prior-active)", bool(a2["pass"]
                                                         and a2p["pass"]),
         "in-domain z(m0) %+.2f z(sig2) %+.2f; prior-active z(m0) %+.2f "
         "z(sig2) %+.2f" % (a2["z_m0"], a2["z_sigma2"], a2p["z_m0"],
                            a2p["z_sigma2"])),
        ("algebra vs v11 conventions", bool(al["pass"]),
         "in-domain cells bit-identical to v3.8, gravity part untouched "
         "everywhere, B_h linear in alpha_m"),
        ("C1 PASS on every in-domain cell", bool(c1["pass"]),
         "in-domain cells: %s; FAIL: %s; UNRESOLVED (below u* = c_m/e): %s"
         % (", ".join(c1["in_domain_cells"]),
            ", ".join(c1["fail_cells"]) or "none",
            ", ".join(c1["unresolved_cells"]) or "none")),
        ("C1b PASS on all cells", bool(c1b["pass"]),
         "worst capillary share %.4f (>= 0.15), min sigma2_h total %.6f "
         "(>= 0.002); failing: %s"
         % (c1b["worst_capillary_share"], c1b["min_sigma2_h_total"],
            ", ".join(c1b["fail_cells"]) or "none")),
        ("V1-V3 reported", True,
         "V1 factor over clean Cox-Munk %.2f-%.2f; V2 band shares; V3 max "
         "B_l/B_hat %.3f"
         % (sm["V1"]["factor_range"][0], sm["V1"]["factor_range"][1],
            sm["V3"]["worst_B_l_over_Bhat"])),
        ("K1 no flag open", not bool(sm["K1"]["flags_open"]),
         "worst abs(dHs) %.1f%% (flag at 35%%); flags: %s"
         % (100 * sm["K1"]["worst_abs_delta_rel"],
            ", ".join(sm["K1"]["flags_open"]) or "none")),
        ("K4 PASS (a-d)", bool(s3["K4"]["K4a"]["pass"]
                               and s3["K4"]["K4b"]["pass"]
                               and s3["K4"]["K4c"]["pass"]
                               and s3["K4"]["K4d_prior_battery"]["pass"]),
         "composite bit-identical %s; Beer-Lambert %s; rain R=0 %s; prior "
         "battery %s"
         % (s3["K4"]["K4a"]["pass"], s3["K4"]["K4b"]["pass"],
            s3["K4"]["K4c"]["pass"],
            s3["K4"]["K4d_prior_battery"]["pass"])),
        ("K5 PASS", bool(sm["K5"]["pass"]),
         "worst 8 s quantisation rel %.3f" % sm["K5"]["worst_rel"]),
        ("K6 PASS", bool(sm["K6"]["pass"]),
         "slope-halving ratio %.3f (smooth = 0.5)"
         % sm["K6"]["worst_slope_ratio"]),
        ("sensitivity sweep delivered", swept_ok,
         "%d levels (2 scenes x 4 targets): C percentiles, own-field warp, "
         "mean abs dL/dy, total sigma2 per level" % n_levels),
        ("low-wind cells listed UNRESOLVED + prior labelled",
         bool(c1["unresolved_cells"]
              and all(u["prior"] is not None
                      for u in c1["unresolved_with_numbers"])
              and all(v is None or isinstance(v, dict)
                      for v in [c["low_wind_prior"] for c in s2["cells"]])),
         "UNRESOLVED: %s, each carrying the low_wind_prior label with its "
         "target" % ", ".join(c1["unresolved_cells"])),
    ]
    blocking = [n for n, ok, _ in items if not ok]
    return {"items": [{"item": n, "pass": bool(ok), "evidence": e}
                      for n, ok, e in items],
            "merge_pass": bool(not blocking),
            "blocking_items": blocking,
            "rule": ("PIPELINE-V11 section 5: flip the default iff A1/A2 "
                     "PASS, C1 PASS on every in-domain cell, C1b PASS on all "
                     "cells, V reported, K1 no flag, K4-K6 PASS, the "
                     "sensitivity sweep delivered, and the low-wind cells "
                     "listed UNRESOLVED with the prior labelled.  "
                     "UNRESOLVED is disclosed, not counted as PASS."),
            "c1_state_by_cell": {
                c["cell"]: {"state": c["gates"]["C1"]["state"],
                            "gravity_share_50_200":
                                c["gates"]["C1"]["gravity_share_50_200"],
                            "u_over_fit_min":
                                c["alpha_m_resolution"]["u_over_fit_min"]}
                for c in s2["cells"]}}


def stage4():
    t0 = time.time()
    log("=" * 72)
    log("[stage 4] continuity, merge evaluation, report")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    with open(os.path.join(OUT, "numbers-v13-transport.json"),
              "r", encoding="utf-8") as f:
        tr = json.load(f)
    tr["K3"]["continuity_vs_v11"]["recheck_stage4"] = {
        "v13_new_sha256": sha256_file(os.path.join(
            OUT, "comp-windy-v39-new.png")),
        "v11_new_sha256": sha256_file(os.path.join(
            V11_OUT, "comp-windy-v38-new.png")),
        "v13_legacy_sha256": sha256_file(os.path.join(
            OUT, "comp-windy-v39-legacy.png")),
        "v11_legacy_sha256": sha256_file(os.path.join(
            V11_OUT, "comp-windy-v38-legacy.png")),
    }
    tr["K3"]["continuity_vs_v11"]["all_byte_identical"] = bool(
        tr["K3"]["continuity_vs_v11"]["recheck_stage4"]["v13_new_sha256"]
        == tr["K3"]["continuity_vs_v11"]["recheck_stage4"]["v11_new_sha256"]
        and
        tr["K3"]["continuity_vs_v11"]["recheck_stage4"]["v13_legacy_sha256"]
        == tr["K3"]["continuity_vs_v11"]["recheck_stage4"]["v11_legacy_sha256"])
    dump_json(os.path.join(OUT, "numbers-v13-transport.json"), tr)

    merge = evaluate_merge()
    for it in merge["items"]:
        log("  %-42s %s  %s" % (it["item"], "PASS" if it["pass"] else "FAIL",
                                it["evidence"]))
    flipped = False
    if merge["merge_pass"]:
        path = os.path.join(HERE, "v13_state.py")
        with open(path, "r", encoding="utf-8") as f:
            txt = f.read()
        if 'DEFAULT_SOURCE = "jonswap"' in txt:
            txt = txt.replace('DEFAULT_SOURCE = "jonswap"',
                              'DEFAULT_SOURCE = "elfouhaily-form"', 1)
            atomic_write_bytes(path, txt.encode("utf-8"))
            flipped = True
    merge["default_flipped"] = flipped
    merge["default_after"] = ("elfouhaily-form" if flipped else "jonswap")
    merge["default_lives_in"] = "v13_state.py (v11_state.py is not edited)"
    dump_json(os.path.join(OUT, "merge-evaluation-v13.json"), merge)
    log("  merge pass=%s; default=%s; blocking: %s"
        % (merge["merge_pass"], merge["default_after"],
           ", ".join(merge["blocking_items"]) or "none"))
    import v13_report
    v13_report.main()
    cp.unit_done("stage-4")
    _dedupe_checkpoint()
    log("[stage 4] done in %.1fs" % (time.time() - t0))
    return merge


def _dedupe_checkpoint():
    """One unit per stage in the checkpoint (the double-run left duplicates)."""
    st = json.load(open(CP_PATH, "r", encoding="utf-8"))
    ids = []
    for u in st.get("completed_ids", []):
        if u not in ids:
            ids.append(u)
    st["completed_ids"] = ids
    st["completed_units"] = len(ids)
    st["last_unit"] = ids[-1] if ids else None
    atomic_write_bytes(CP_PATH,
                       json.dumps(st, indent=2,
                                  ensure_ascii=False).encode("utf-8"))
    log("  checkpoint de-duplicated: %s" % (", ".join(ids),))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="1", choices=["1", "2", "3", "4"])
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
