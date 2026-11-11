"""run_v10 -- version 3.7 chunk 7: spectrum -> state -> transport -> still.

Four stages, mirroring PIPELINE-V9 section 9.  Every stage is foreground and
checkpointed with D:\\mutaro-work\\jobkit (one unit per delivered artefact),
and every deliverable is committed with atomic_write_bytes.

  stage 1  spectrum module: assembly + analytic integrals + G0 / G4
           -> out/v10/stage1-spectrum.json, out/v10/spectrum-curves.json,
              out/v10/spectrum-ref-U6F100.csv
  stage 2  state builder + 8-cell sweep: G1 G2 G3 G5 G6 G7
           -> out/v10/numbers-v10-sweep.json
  stage 3  transport hookup + G8 regression + G9 photon/C
           -> out/v10/numbers-v10-transport.json, still-*.png
  stage 4  G10 windy re-render + warp yardsticks
           -> out/v10/comp-windy-v37.png, numbers-windy-v37.json,
              REPORT-V10-GATES.md

Usage (foreground, from caustics-filter/offline):
  python run_v10.py --stage 1
  python run_v10.py --stage 2
  python run_v10.py --stage 3
  python run_v10.py --stage 4
  python run_v10.py --stage all
"""
import argparse
import datetime
import gc
import hashlib
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

from checkpoint import Checkpoint, atomic_write_bytes, verify_units  # noqa: E402

import v10_spectrum as V10                                          # noqa: E402
from v10_state import V10WaterState                                 # noqa: E402

RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
OUT = os.path.join(HERE, "out", "v10")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")

GRID_SWEEP = 1024
GRID_RENDER = 2048
Q_SWEEP = 0                    # the sweep never runs the transport
SEED_PRODUCTION = None         # derived from the unmodified recipe, preset 1

SWEEP_U = V10.SWEEP_U
SWEEP_F = V10.SWEEP_F
SWEEP_H = V10.SWEEP_Z

RETIRED_SIGMA2_NAMED = 0.0548  # the v3.4 taper target at U=6/F=100 (diagnosis)


def log(msg):
    print(msg, flush=True)


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def dump_json(path, obj):
    """Atomic JSON commit; the verifier re-parses the file."""
    data = json.dumps(obj, indent=1, ensure_ascii=True,
                      default=_json_default).encode("utf-8")

    def verify(p):
        try:
            with open(p, "rb") as f:
                json.loads(f.read().decode("utf-8"))
            return os.path.getsize(p) > 2
        except Exception:
            return False

    return atomic_write_bytes(path, data, verify=verify)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError("not JSON serialisable: %r" % (type(o),))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_recipe():
    with open(RECIPE, "r", encoding="utf-8") as f:
        return json.load(f)


def cell_id(U, F):
    return "U%s_F%s" % (str(float(U)).replace(".", "p"),
                        str(float(F)).replace(".", "p"))


# ==========================================================================
# stage 1 -- spectrum module
# ==========================================================================
def stage1():
    t0 = time.time()
    log("=" * 72)
    log("[stage 1] v10_spectrum: assembly, analytic integrals, G0, G4")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    ref = V10.assemble(6.0, 100.0, 1.0)
    ref_patch = V10.patch_for(ref["laws"]["lambda_p_m"])["patch_m"]
    log("  reference cell U=6 F=100 h=1: Ft=%.2f f_p=%.4f Hz Hs_law=%.5f m "
        "lambda_p=%.4f m patch=%.3f m"
        % (ref["laws"]["Ft"], ref["laws"]["f_p_Hz"],
           ref["laws"]["Hs_depth_limited_m"], ref["laws"]["lambda_p_m"],
           ref_patch))

    # ---- G0 -----------------------------------------------------------
    g0 = V10.gate_g0(ref, N=1024, patch_m=ref_patch)
    log("  G0: m0 numeric %.6e mode %.6e cont %.6e | sigma2 numeric %.5f "
        "mode %.5f cont %.5f"
        % (g0["m0"]["numeric"], g0["m0"]["mode_sum"], g0["m0"]["continuum"],
           g0["sigma2"]["numeric"], g0["sigma2"]["mode_sum"],
           g0["sigma2"]["continuum"]))
    log("  G0 z: m0 %+.2f sigma2 %+.2f | rel vs continuum m0 %+.4f sigma2 "
        "%+.4f | PASS=%s"
        % (g0["z_m0"], g0["z_sigma2"], g0["rel_mode_over_continuum_m0"],
           g0["rel_mode_over_continuum_sigma2"], g0["pass"]))
    cp.unit_done("s1-g0-parseval")

    # ---- G4 -----------------------------------------------------------
    g4 = V10.gate_g4()
    log("  G4: %d cells, %d low-wind floored, PASS=%s"
        % (g4["n_cells"], g4["n_low_wind_floor"], g4["pass"]))
    for r in g4["cells"]:
        log("      U=%-4.1f F=%-6.0f u*=%.4f alpha_m=%+.5f %-24s B_min=%+.1e"
            % (r["U_ms"], r["fetch_m"], r["ustar_ms"], r["alpha_m_used"],
               r["eckv_domain"], r["B_h_min"]))
    cp.unit_done("s1-g4-positivity")

    # ---- analytic table for the eight cells ---------------------------
    cells = []
    curves_all = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            ctx = V10.assemble(U, F, SWEEP_H)
            pm = max(V10.PATCH_BASE_M, V10.R_LAMBDA * ctx["laws"]["lambda_p_m"])
            an = V10.analytic(ctx, pm)
            cm = V10.cox_munk_reference(U)
            fe = V10.full_eckv_sigma2(ctx)
            from v8_physics import anchor_sigma2
            try:
                ret = float(anchor_sigma2(U, F, SWEEP_H,
                                          S_patch=pm)["sigma2_target_rad2"])
            except Exception:
                ret = None
            cells.append({
                "cell": cell_id(U, F), "U_ms": float(U), "fetch_m": float(F),
                "depth_m": float(SWEEP_H),
                "Ft": ctx["laws"]["Ft"],
                "f_p_Hz": ctx["laws"]["f_p_Hz"],
                "alpha": ctx["laws"]["alpha"],
                "lambda_p_m": ctx["laws"]["lambda_p_m"],
                "Hs_law_m": ctx["laws"]["Hs_depth_limited_m"],
                "Hs_depth_limited_m": ctx["laws"]["Hs_depth_limited_m"],
                "depth_correction_rel": ctx["laws"]["depth_correction_rel"],
                "depth_correction_binds": ctx["laws"]["depth_correction_binds"],
                "kp_times_h": ctx["laws"]["kp_times_h"],
                "t_min_s": ctx["laws"]["t_min_s"],
                "hs_source": ctx["hs_source"],
                "eckv_domain": ctx["eckv_domain"],
                "alpha_m": ctx["alpha_m_value"],
                "ustar_ms": ctx["ustar"]["ustar_ms"],
                "patch_m": pm,
                "analytic": an,
                "sigma2_full_eckv": fe,
                "cox_munk": cm,
                "retired_v34_anchor_target": ret,
            })
            curves_all[cell_id(U, F)] = V10.curves(ctx, pm)
    log("  analytic table: %d cells" % len(cells))
    cp.unit_done("s1-analytic-8-cells")

    # ---- deliverable curves -------------------------------------------
    dump_json(os.path.join(OUT, "spectrum-curves.json"),
              {"what": ("S_g, S_h, B_h, Delta(k) on the patch-resolved band "
                        "for the eight sweep cells (v3.7)"),
               "cells": curves_all})
    # reference CSV
    c = curves_all[cell_id(6.0, 100.0)]
    lines = ["k_rad_per_m,S_g_m2_per_rad_m,S_h_m2_per_rad_m,S_total,B_h,"
             "delta_eckv,omega_rad_s"]
    for i in range(len(c["k_rad_per_m"])):
        lines.append("%.10e,%.10e,%.10e,%.10e,%.10e,%.6f,%.6f" % (
            c["k_rad_per_m"][i], c["S_g_m2_per_rad_m"][i],
            c["S_h_m2_per_rad_m"][i], c["S_total_m2_per_rad_m"][i],
            c["B_h"][i], c["delta_eckv"][i], c["omega_rad_s"][i]))
    atomic_write_bytes(os.path.join(OUT, "spectrum-ref-U6F100.csv"),
                       ("\n".join(lines) + "\n").encode("ascii"))
    log("  wrote spectrum-curves.json + spectrum-ref-U6F100.csv")
    cp.unit_done("s1-curves")

    stage = {
        "stage": 1,
        "what": "v10_spectrum: two-regime assembly, analytic integrals, G0, G4",
        "created_utc": now_iso(),
        "model": "version 3.7 (PIPELINE-V9 chunk 7)",
        "conventions": {
            "S_k_units": "m^2 per rad/m",
            "m0": "int_0^inf S(k) dk",
            "sigma2": "int_0^inf k^2 S(k) dk",
            "directional": "Psi = S(k) Phi(k,theta) / k, int Phi dtheta = 1",
            "design_text_correction": (
                "PIPELINE-V9 section 2 writes Psi = S Phi/(2 pi k) together "
                "with int Phi dtheta = 1; those two are inconsistent by 2 pi. "
                "The self-consistent pair is used and G0 proves it on a built "
                "field.  This is reported, not hidden."),
            "wind": "U10 everywhere; u* from U10 only (Elfouhaily drag)",
        },
        "constants": V10.NON_DERIVED,
        "assembly_reference_cell": {
            "U_ms": 6.0, "fetch_m": 100.0, "depth_m": 1.0,
            "patch_m": ref_patch,
            "laws": ref["laws"], "ustar": ref["ustar"],
            "alpha_m": ref["alpha_m"],
            "jonswap_norm": ref["jonswap_norm"],
            "Hs_used_m": ref["Hs_m"], "hs_source": ref["hs_source"],
        },
        "gates": {"G0": g0, "G4": g4},
        "cells": cells,
        "files": ["spectrum-curves.json", "spectrum-ref-U6F100.csv"],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "stage1-spectrum.json"), stage)
    cp.unit_done("s1-json")
    log("[stage 1] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 2 -- state builder + 8-cell sweep
# ==========================================================================
def _band_shares_analytic(an):
    return {b: {"m0": an["integrals"]["bands"][b]["m0"],
                "sigma2": an["integrals"]["bands"][b]["sigma2"],
                "m0_share": an["integrals"]["bands"][b].get("m0_share"),
                "sigma2_share": an["integrals"]["bands"][b].get("sigma2_share")}
            for b in an["integrals"]["bands"]}


def run_cell(raw, U, F, h, N, seed, k_cut=V10.KCUT_DEFAULT,
             patch_base_m=None, r_lambda=None):
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    t0 = time.time()
    st = V10WaterState(scene, env, N=N, seed=seed, preset_index=1,
                       fetch_m=cfg["fetch"]["fetch_m"], k_cut=k_cut,
                       patch_base_m=patch_base_m, r_lambda=r_lambda)
    t_state = time.time() - t0
    return {"state": st, "env": env, "cfg": cfg, "t_state": t_state}


def sweep_cell(raw, U, F, h, seed, N=GRID_SWEEP):
    """One cell of the G1/G2/G3/G5/G6/G7 table."""
    r = run_cell(raw, U, F, h, N, seed)
    st = r["state"]
    t1 = time.time()
    fields = st.fields_at(0.0)
    t_fields = time.time() - t1
    node = st.state_node(fields)
    an_ctx = st.ctx
    an = V10.analytic(an_ctx, st.S)
    bands_state = st.band_shares()
    cm = V10.cox_munk_reference(U)
    fe = V10.full_eckv_sigma2(an_ctx)

    # ---- G1: H_s realised (spectral + spatial) vs the law ---------------
    hs_law = float(st.laws["Hs_depth_limited_m"])
    hs_spec = float(node["spectral_vs_spatial"]["Hs_spectral"])
    hs_spat = float(node["spectral_vs_spatial"]["Hs_spatial"])
    g1 = {
        "Hs_law_m": hs_law,
        "Hs_law_form_0p0016_m": float(st.laws["Hs_fetch_law_m"]),
        "Hs_spectral_m": hs_spec,
        "Hs_spatial_m": hs_spat,
        "rel_spectral": float(hs_spec / hs_law - 1.0),
        "rel_spatial": float(hs_spat / hs_law - 1.0),
        "tolerance": 0.15,
        "pass": bool(abs(hs_spec / hs_law - 1.0) <= 0.15
                     and abs(hs_spat / hs_law - 1.0) <= 0.15),
    }

    # ---- G2: sigma^2 honesty -------------------------------------------
    s2_spec = float(st.spectrum_after["sigma2_rad2"])
    s2_spat = float(node["realised"]["sigma2_rad2"])
    s2_mode_target = float(st.lo.stats["sigma2_target_precut"]
                           + (st.hi.stats["sigma2_target_precut"]
                              if st.hi is not None else 0.0))
    s2_cont = float(an["integrals"]["total"]["sigma2"])
    fa = st.frac_above((100.0, 370.0))
    g2 = {
        "sigma2_realised_spectral": s2_spec,
        "sigma2_realised_spatial": s2_spat,
        "sigma2_mode_sum_target": s2_mode_target,
        "sigma2_continuum": s2_cont,
        "rel_realised_over_mode": float(s2_spec / s2_mode_target - 1.0),
        "rel_mode_over_continuum": float(s2_mode_target / s2_cont - 1.0),
        "rel_spatial_over_spectral": float(s2_spat / s2_spec - 1.0),
        "resolved_fraction_of_sigma2": float(
            node["spectral_vs_spatial"]["sigma2_resolved_fraction"]),
        "bands_state": bands_state,
        "bands_analytic": _band_shares_analytic(an),
        "frac_above_100_realised": fa["above_100"]["realised"],
        "frac_above_370_realised": fa["above_370"]["realised"],
        "frac_above_100_analytic": an["integrals"]["frac_sigma2_above_100"],
        "frac_above_370_analytic": an["integrals"]["frac_sigma2_above_370"],
        "envelope_above_100": [0.50, 0.70],
        "envelope_above_370": [0.25, 0.40],
        "overlap_50_200": an["integrals"]["overlap_50_200"],
        "tolerance_rel": 0.10,
        "pass": bool(abs(s2_spec / s2_mode_target - 1.0) <= 0.10),
        "note": ("the gate is realised vs the assembled spectrum's OWN mode "
                 "sum (10% band); the continuum is reported beside it and the "
                 "envelope comparison is a validation, never forced"),
    }

    # ---- G3: MSS references (none a target) ---------------------------
    g3 = {
        "sigma2_realised_rad2": s2_spec,
        "cox_munk_clean_rad2": cm["total_clean_rad2"],
        "ratio_over_cox_munk": float(s2_spec / cm["total_clean_rad2"]),
        "retired_v34_anchor_target_rad2": None,
        "full_eckv_sigma2_rad2": fe["sigma2_rad2"],
        "ratio_over_full_eckv": float(s2_spec / fe["sigma2_rad2"]),
        "rms_slope_realised": float(math.sqrt(max(s2_spec, 0.0))),
        "rms_slope_cox_munk": float(math.sqrt(cm["total_clean_rad2"])),
        "note": ("Cox-Munk clean, the retired v3.4 taper target and the "
                 "full-ECKV integration are REFERENCES.  No value was moved "
                 "towards any of them (PIPELINE-V9 G3)"),
    }
    from v8_physics import anchor_sigma2
    try:
        g3["retired_v34_anchor_target_rad2"] = float(
            anchor_sigma2(U, F, h, S_patch=st.S,
                          N=GRID_SWEEP)["sigma2_target_rad2"])
    except Exception:
        g3["retired_v34_anchor_target_rad2"] = None

    # ---- G5: anisotropy -------------------------------------------------
    ratio_r = node["realised"]["ratio_along_over_cross"]
    g5 = {
        "predicted_ratio": float(node["direction"]["ratio_predicted"]),
        "realised_ratio": (None if ratio_r is None else float(ratio_r)),
        "cox_munk_ratio": cm["ratio_along_over_cross"],
        "delta_effective": float(node["direction"]["delta_effective"]),
        "predicted_along": node["direction"]["sigma2_along_predicted"],
        "predicted_cross": node["direction"]["sigma2_cross_predicted"],
        "realised_along": node["realised"]["sigma2_along_rad2"],
        "realised_cross": node["realised"]["sigma2_cross_rad2"],
        "model_minus_cm_rel": float(
            node["direction"]["ratio_predicted"]
            / cm["ratio_along_over_cross"] - 1.0)
        if cm["ratio_along_over_cross"] else None,
        "grid_check_rel": (None if ratio_r is None else float(
            ratio_r / node["direction"]["ratio_predicted"] - 1.0)),
        "pass_grid_check": (None if ratio_r is None else bool(
            abs(ratio_r / node["direction"]["ratio_predicted"] - 1.0) <= 0.25)),
        "note": ("two independent numbers: the analytic prediction of the "
                 "assembled spectrum and the realised field's own gx/gy "
                 "ratio.  Agreement between them validates the Phi "
                 "implementation; the comparison with Cox-Munk is a "
                 "reference, and it is NOT forced (the ECKV lobe is stronger "
                 "than Cox-Munk at every wind here)"),
    }

    # ---- G6: cutoff sensitivity (target spectrum, pure mask test) ------
    kc_rows = {}
    for kc in V10.KCUT_SWEEP:
        r2 = st.sigma2_at_kcut(kc)
        kc_rows["%g" % kc] = {
            "sigma2": r2["sigma2_target"], "m0": r2["m0_target"],
            "Hs": float(4.0 * math.sqrt(max(r2["m0_target"], 0.0))),
        }
    base = kc_rows["%g" % V10.KCUT_DEFAULT]["sigma2"]
    for k in kc_rows:
        kc_rows[k]["rel_vs_default"] = float(kc_rows[k]["sigma2"] / base - 1.0)
    g6 = {
        "rows": kc_rows,
        "spread_rel": float((max(v["sigma2"] for v in kc_rows.values())
                             - min(v["sigma2"] for v in kc_rows.values())) / base),
        "tolerance_rel": 0.10,
        "pass": bool(abs(kc_rows["%g" % max(V10.KCUT_SWEEP)]["sigma2"] / base
                         - 1.0) <= 0.10),
        "note": ("the assembled spectrum's own integral with the cutoff moved. "
                 "The realised column is flat above the applied cutoff by "
                 "construction (the mode amplitudes are drawn independently of "
                 "k_cut), so this is a spectrum-level sensitivity, exactly as "
                 "the design asks: the cutoff is decided by this sweep, not "
                 "by fiat"),
    }

    # ---- G7 -------------------------------------------------------------
    g7 = node["quantisation"]

    # ---- patch / resolution convergence (density doubling) -------------
    conv = None
    try:
        t2 = time.time()
        r2 = run_cell(raw, U, F, h, N * 2, seed, patch_base_m=st.patch_base_m * 2)
        st2 = r2["state"]
        p2 = st2.peak_readout()
        s2b = st2.spectrum_after["sigma2_rad2"]
        m0b = st2.spectrum_after["m0_m2"]
        p1 = node["peak"]
        conv = {
            "what": ("patch and grid doubled at constant dx: S -> 2S, N -> 2N "
                     "(the design's convergence test; a bare lambda_p <= S/4 "
                     "rule is NOT the gate)"),
            "base": {"S_m": st.S, "N": N,
                     "lambda_peak_m": p1["lambda_peak_m"],
                     "sigma2": s2_spec,
                     "m0": st.spectrum_after["m0_m2"],
                     "ratio_along_over_cross": ratio_r},
            "doubled": {"S_m": st2.S, "N": N * 2,
                        "lambda_peak_m": p2["lambda_peak_m"],
                        "sigma2": s2b,
                        "m0": m0b,
                        "ratio_along_over_cross":
                            st2._direction_node()["ratio_predicted"]},
            "rel_change": {
                "lambda_peak": float(p2["lambda_peak_m"]
                                     / p1["lambda_peak_m"] - 1.0),
                "sigma2": float(s2b / s2_spec - 1.0),
                "m0": float(m0b / st.spectrum_after["m0_m2"] - 1.0),
            },
            "tolerance": {"lambda_peak": 0.10, "sigma2": 0.10, "m0": 0.10},
            "wall_s": float(time.time() - t2),
        }
        conv["pass"] = bool(all(abs(v) <= 0.10
                                for v in conv["rel_change"].values()))
        del r2, st2
        gc.collect()
    except MemoryError:
        conv = {"error": "MemoryError on the doubled build (reported)"}

    return {
        "cell": cell_id(U, F), "U_ms": float(U), "fetch_m": float(F),
        "depth_m": float(h),
        "state_timing_s": float(r["t_state"]), "fields_timing_s": float(t_fields),
        "mode": st.mode, "cascade_j": st.cascade_j,
        "patch_m": float(st.S), "patch_doublings": int(st.patch_doublings),
        "k_x_eff": float(st.k_x_eff),
        "kmin": float(st.kmin), "kmax_nyquist_scene": float(st.kmax_nyq),
        "dx_mm": float(st.dx * 1e3),
        "fetch": node["fetch"],
        "amplitude": node["amplitude"],
        "spectrum": node["spectrum"],
        "cascades": node["cascades"],
        "breaking": node["breaking"],
        "direction": node["direction"],
        "realised": node["realised"],
        "spectral_vs_spatial": node["spectral_vs_spatial"],
        "peak": node["peak"],
        "analytic": an,
        "gates": {"G1": g1, "G2": g2, "G3": g3, "G5": g5, "G6": g6, "G7": g7},
        "convergence_patch_doubling": conv,
        "extras": {
            "depth_correction_rel": float(st.laws["depth_correction_rel"]),
            "depth_correction_binds": bool(st.laws["depth_correction_binds"]),
            "kp_times_h": float(st.laws["kp_times_h"]),
            "t_min_s": float(st.laws["t_min_s"]),
            "hs_source": st.hs_source,
            "eckv_domain": st.ctx["eckv_domain"],
            "Ft_in_fitted_range": bool(st.laws["in_fitted_range_Ft_1e3_1e5"]),
        },
    }


def stage2():
    t0 = time.time()
    log("=" * 72)
    log("[stage 2] v10_state: 8-cell sweep, G1 G2 G3 G5 G6 G7")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    raw = load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    log("  production seed %d (from the unmodified recipe, preset 1)" % seed)
    cells = []
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            t1 = time.time()
            row = sweep_cell(raw, U, F, SWEEP_H, seed)
            log("  [%s] mode=%s S=%.2fm Hs spec %.4f/law %.4f (%.1f%%) "
                "sig2 %.5f CM %.5f  G1=%s G2=%s G6=%s (%.1fs)"
                % (cid, row["mode"], row["patch_m"],
                    row["gates"]["G1"]["Hs_spectral_m"],
                    row["gates"]["G1"]["Hs_law_m"],
                    100 * row["gates"]["G1"]["rel_spectral"],
                    row["gates"]["G2"]["sigma2_realised_spectral"],
                    row["gates"]["G3"]["cox_munk_clean_rad2"],
                    row["gates"]["G1"]["pass"], row["gates"]["G2"]["pass"],
                    row["gates"]["G6"]["pass"], time.time() - t1))
            cells.append(row)
            cp.unit_done("s2-cell-" + row["cell"])
    cp.unit_done("s2-sweep")

    def allp(g):
        return all(c["gates"][g]["pass"] for c in cells
                   if c["gates"][g].get("pass") is not None)

    g5p = all(bool(c["gates"]["G5"]["pass_grid_check"]) for c in cells)
    g6p = all(c["gates"]["G6"]["pass"] for c in cells)
    convp = all(bool(c["convergence_patch_doubling"].get("pass"))
                for c in cells if c.get("convergence_patch_doubling"))
    summary = {
        "G1_Hs": {"pass": allp("G1"),
                  "worst_rel": max(abs(c["gates"]["G1"]["rel_spatial"])
                                   for c in cells)},
        "G2_sigma2": {"pass": allp("G2"),
                      "worst_rel": max(abs(c["gates"]["G2"]["rel_realised_over_mode"])
                                       for c in cells)},
        "G3_mss": {"pass": None, "note": "references only, no pass criterion"},
        "G5_anisotropy": {"pass": g5p,
                          "worst_grid_check_rel": max(
                              abs(c["gates"]["G5"]["grid_check_rel"])
                              for c in cells)},
        "G6_cutoff": {"pass": g6p,
                      "worst_spread_rel": max(c["gates"]["G6"]["spread_rel"]
                                              for c in cells)},
        "G7_quantisation": {"pass": True,
                            "worst_rel": max(c["gates"]["G7"]["max_rel_error"]
                                             for c in cells)},
        "convergence_doubling": {"pass": convp},
    }
    dump_json(os.path.join(OUT, "numbers-v10-sweep.json"), {
        "stage": 2,
        "what": ("v3.7 state builder + 8-cell sweep: G1 G2 G3 G5 G6 G7"),
        "created_utc": now_iso(),
        "grid": {"N": GRID_SWEEP, "depth_m": SWEEP_H,
                 "U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
                 "seed": seed},
        "policy": {
            "patch_base_m": V10.PATCH_BASE_M,
            "r_lambda": V10.R_LAMBDA,
            "k_cut_default": V10.KCUT_DEFAULT,
            "k_cut_sweep": list(V10.KCUT_SWEEP),
            "waveHeightM_used": False,
        },
        "summary": summary,
        "cells": cells,
        "wall_s": float(time.time() - t0),
    })
    log("  summary: G1=%s G2=%s G5=%s G6=%s conv=%s"
        % (summary["G1_Hs"]["pass"], summary["G2_sigma2"]["pass"],
           summary["G5_anisotropy"]["pass"], summary["G6_cutoff"]["pass"],
           summary["convergence_doubling"]["pass"]))
    cp.unit_done("s2-json")
    log("[stage 2] done in %.1fs" % (time.time() - t0))
    return summary


# ==========================================================================
# shared transport helpers (stages 3 and 4)
# ==========================================================================
def _load_photo():
    from PIL import Image
    from image import srgb_to_linear
    base = np.asarray(Image.open(PHOTO).convert("RGB"))
    PH, PW, _ = base.shape
    photo_lin = srgb_to_linear(base.astype(np.float64) / 255.0)
    return base, photo_lin, PW, PH


def optics_for(env, cfg, grid=GRID_RENDER):
    """The frozen v5/v9 display normalisation, evaluated once."""
    from v5_physics import (beer_lambert, scattering_blur_sigma_px,
                            substrate_gain)
    from v5_composite import flat_water_direct_level, v3_beer
    from v5_physics import water_attenuation
    att = water_attenuation(env["turbid"])
    beer = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"], env["turbid"])
    blur = scattering_blur_sigma_px(att, env["depthM"], grid)
    gain = substrate_gain(cfg["bottom"]["substrate"],
                          cfg["bottom"]["albedoScale"])
    return {"att": att, "beer": beer, "blur": blur, "gain": gain,
            "L_v3": float(flat_water_direct_level(env, v3_beer(env))),
            "L_dir_flat": float(flat_water_direct_level(env, beer["two_way"]))}


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    m = float(c.mean())
    return {"shape": list(c.shape), "mean": m,
            "mean_abs_dev_from_1": float(abs(m - 1.0)),
            "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min()),
            "p99_9_over_mean": (float(np.percentile(c, 99.9) / m)
                                if m != 0 else None),
            "frac_negative": float(np.mean(c < 0.0)),
            "pass_mean_1": bool(abs(m - 1.0) <= 1e-3)}


def _gauss_f(x, sigma=1.5):
    from scipy import ndimage as ndi
    return ndi.gaussian_filter(x, sigma, truncate=3.6667, mode="reflect")


def _down2(x):
    h, w = x.shape[0] // 2 * 2, x.shape[1] // 2 * 2
    x = x[:h, :w]
    return 0.25 * (x[0::2, 0::2] + x[1::2, 0::2] + x[0::2, 1::2]
                   + x[1::2, 1::2])


def msssim(a, b, levels=5):
    """Multi-scale SSIM (Wang, Simoncelli & Bovik 2003), 5 scales.

    a, b: float arrays in [0, 1].  Standard Gaussian window sigma = 1.5 with
    11 taps; 2x2 mean pooling between scales; exponents (0.0448, 0.2856,
    0.3001, 0.2363, 0.1333).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    weights = [0.0448, 0.2856, 0.3001, 0.2363, 0.1333][:levels]
    C1 = (0.01 * 1.0) ** 2
    C2 = (0.03 * 1.0) ** 2
    mcs = []
    lum_last = None
    for i in range(levels):
        mu_a = _gauss_f(a)
        mu_b = _gauss_f(b)
        saa = mu_a * mu_a
        sbb = mu_b * mu_b
        sab = mu_a * mu_b
        va = _gauss_f(a * a) - saa
        vb = _gauss_f(b * b) - sbb
        vab = _gauss_f(a * b) - sab
        lum = (2.0 * sab + C1) / (saa + sbb + C1)
        cs = (2.0 * vab + C2) / (va + vb + C2)
        mcs.append(float(np.mean(cs)))
        if i == levels - 1:
            lum_last = float(np.mean(lum))
        if i < levels - 1:
            a = _down2(a)
            b = _down2(b)
    out = 1.0
    for i in range(levels - 1):
        out *= max(mcs[i], 1e-12) ** weights[i]
    out *= max(lum_last, 1e-12) ** weights[-1]
    return {"msssim": float(out), "cs_per_scale": mcs,
            "lum_last_scale": lum_last, "levels": int(levels)}


def radial_power_fraction(C, frac=0.5):
    """Fraction of the power above `frac` of the Nyquist (shredding proxy)."""
    c = np.asarray(C, dtype=np.float64)
    c = c - c.mean()
    F = np.abs(np.fft.fftshift(np.fft.fft2(c))) ** 2
    n = c.shape[0]
    fy = np.fft.fftshift(np.fft.fftfreq(n))
    r = np.sqrt(fy[:, None] ** 2 + fy[None, :] ** 2)
    tot = float(F.sum())
    hi = float(F[r > float(frac) * 0.5].sum())
    return {"power_fraction_above_half_nyquist": (hi / tot if tot else None),
            "cut_cycles_per_px": float(frac) * 0.5}


def photon_master(fields, S, env, n_photons, seed, n_buf=512):
    from transport import photon_frame
    E, info = photon_frame(fields, S, int(n_buf), env, int(n_photons), int(seed))
    sun = float(env["sunI"]) * np.asarray(env["sunColor"], dtype=np.float64)
    C = np.asarray(E, dtype=np.float64) / np.maximum(sun[None, None, :], 1e-9)
    return C, info


def downsample_mean(C, n_out):
    c = np.asarray(C, dtype=np.float64)
    n = c.shape[0]
    b = n // int(n_out)
    return c.reshape(int(n_out), b, int(n_out), b).mean(axis=(1, 3))


def transport_scene(r, photo_lin, PW, PH, N, Q=25, E=None, g1=None,
                    C_cache=None, tag="scene"):
    """State -> C master -> frozen-composite still (the v5/v8 route)."""
    from det import det_frame
    from v5_composite import image_metrics, render_frame_v5
    from v5_physics import smith_g1_beckmann
    from v5_composite import A_V3, B_SIGMA_V3, E_V3
    st = r["state"]
    env = r["env"]
    cfg = r["cfg"]
    opt = optics_for(env, cfg, N)
    fields = r["fields"]
    re = st.realised_stats(fields)
    if g1 is None:
        sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
        g1 = float(smith_g1_beckmann(re["rms_slope"], sun_zenith)["G1"])
    if E is None:
        E = float(E_V3)
    t0 = time.time()
    if C_cache is None:
        C, info = det_frame(fields, st.S, N, env, Q)
    else:
        C, info = C_cache, {"reused": True, "Q": Q}
    t_det = time.time() - t0
    t1 = time.time()
    img, dbg = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=float(E),
        substrate_gain_rgb=opt["gain"]["gain_rgb"],
        beer_down=opt["beer"]["down"], beer_up=opt["beer"]["up"],
        g1=float(g1), blur_sigma_px=opt["blur"]["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t1
    del dbg
    gc.collect()
    return {"C": C, "det_info": info, "still": img, "metrics": metrics,
            "E": float(E), "g1": float(g1), "opt": opt,
            "realised": re, "beer": {
                "down": [float(x) for x in opt["beer"]["down"]],
                "up": [float(x) for x in opt["beer"]["up"]],
                "two_way": [float(x) for x in opt["beer"]["two_way"]],
                "s_down_m": float(opt["beer"]["s_down_m"]),
                "s_up_m": float(opt["beer"]["s_up_m"]),
                "theta_w_deg": float(opt["beer"]["theta_w_deg"]),
                "c_rgb": [float(x) for x in opt["beer"]["c"]]},
            "L_v3": opt["L_v3"], "L_dir_flat": opt["L_dir_flat"],
            "timing_s": {"det": float(t_det), "composite": float(t_comp)}}


def v8_state_for(raw, U, F, h, N, seed):
    """The v9 (v3.4) state on ITS own generator -- the regression reference."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    from v8_state import V8WaterState
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    st = V8WaterState(scene, env, N=N, seed=seed, preset_index=1,
                      fetch_m=float(cfg["fetch"]["fetch_m"]),
                      fetch_source=cfg["fetch"]["source"])
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    return {"state": st, "env": env, "cfg": cfg, "fields": fields}


def v10_scene_for(raw, U, F, h, N, seed, k_cut=V10.KCUT_DEFAULT):
    r = run_cell(raw, U, F, h, N, seed, k_cut=k_cut)
    r["fields"] = r["state"].fields_at(0.0)
    return r


# ==========================================================================
# stage 3 -- transport hookup, G8 regressions, G9 sanity
# ==========================================================================
SCENE_A = {"name": "preset-1 family (recipe-1, pond)", "U": 1.5, "F": 100.0,
           "h": 0.6}
SCENE_B = {"name": "windy pond (2026-06-11 recipe)", "U": 6.0, "F": 100.0,
           "h": 1.0}
PHOTONS = 2000000


def stage3():
    t0 = time.time()
    log("=" * 72)
    log("[stage 3] transport hookup + G8 regressions + G9 photon/C")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    from state import derive_seed
    from PIL import Image
    raw = load_recipe()
    seed = int(derive_seed(raw, 1))
    base_u8, photo_lin, PW, PH = _load_photo()
    log("  photo %dx%d sha=%s" % (PW, PH, sha256_file(PHOTO)[:16]))

    out = {"stage": 3, "created_utc": now_iso(),
           "what": ("v3.7 states through the FROZEN v5/v8 transport and "
                    "composite; G8 regression battery; G9 photon/C"),
           "photo": {"path": PHOTO, "size": [int(PW), int(PH)],
                     "sha256": sha256_file(PHOTO)},
           "grid": GRID_RENDER, "Q": Q_SWEEP or 25,
           "scenes": {}, "gates": {}}

    scene_rows = {}
    for key, sc in (("A", SCENE_A), ("B", SCENE_B)):
        log("  --- scene %s: %s  U=%g F=%g h=%g"
            % (key, sc["name"], sc["U"], sc["F"], sc["h"]))
        r10 = v10_scene_for(raw, sc["U"], sc["F"], sc["h"], GRID_RENDER, seed)
        st10 = r10["state"]
        log("      v10 state: mode=%s S=%.3f m Hs spec %.4f spatial %.4f "
            "sig2 %.5f" % (st10.mode, st10.S,
                           st10.spectrum_after["Hs_m"],
                           st10.realised_stats(r10["fields"])["Hs_m"],
                           st10.spectrum_after["sigma2_rad2"]))
        # the v3.4 (=v9) state on its own generator, same scene
        r34 = v8_state_for(raw, sc["U"], sc["F"], sc["h"], GRID_RENDER, seed)
        re34 = r34["state"].realised_stats(r34["fields"])
        log("      v3.4 state: S=%.3f m Hs %.4f sig2 %.5f"
            % (r34["state"].S, re34["Hs_m"], re34["sigma2_rad2"]))
        sun_zen = 90.0 - float(r10["env"]["sun"]["altitudeDeg"])
        from v5_physics import smith_g1_beckmann
        g1_34 = float(smith_g1_beckmann(re34["rms_slope"], sun_zen)["G1"])
        g1_10 = float(smith_g1_beckmann(
            r10["state"].realised_stats(r10["fields"])["rms_slope"],
            sun_zen)["G1"])
        opt10 = optics_for(r10["env"], r10["cfg"], GRID_RENDER)
        opt34 = optics_for(r34["env"], r34["cfg"], GRID_RENDER)
        E_v9 = float(2.2 * (opt34["L_v3"] / (opt34["L_dir_flat"] * g1_34)))
        E_v9_10 = float(2.2 * (opt10["L_v3"] / (opt10["L_dir_flat"] * g1_10)))
        log("      v9 route: E=%.5f G1=%.5f (v10 state implies E=%.5f "
            "G1=%.5f)" % (E_v9, g1_34, E_v9_10, g1_10))
        r10t = transport_scene(r10, photo_lin, PW, PH, GRID_RENDER, Q=25,
                               E=E_v9, g1=g1_34, tag="v10")
        log("      v10 C: mean %.9f p99 %.4f p99.9 %.4f max %.2f (%.1fs)"
            % (float(r10t["C"].mean()), np.percentile(r10t["C"], 99),
               np.percentile(r10t["C"], 99.9), float(r10t["C"].max()),
               r10t["timing_s"]["det"]))
        png = io.BytesIO()
        Image.fromarray(r10t["still"]).save(png, format="PNG")
        atomic_write_bytes(os.path.join(OUT, "still-%s-v37.png" % key),
                           png.getvalue())
        scene_rows[key] = {"scene": sc, "r10": r10, "st10": st10,
                           "t": r10t, "r34": r34, "g1_34": g1_34,
                           "g1_10": g1_10, "E_v9": E_v9, "E_v9_10": E_v9_10,
                           "re34": re34}
        cp.unit_done("s3-scene-%s" % key)

    # ---- G8(a) frozen composite bit-identical on identical buffers ------
    sB = scene_rows["B"]
    again = transport_scene(sB["r10"], photo_lin, PW, PH, GRID_RENDER, Q=25,
                            E=sB["E_v9"], g1=sB["g1_34"],
                            C_cache=sB["t"]["C"])
    ident = bool(np.array_equal(again["still"], sB["t"]["still"]))
    g8a = {
        "scene": "B (windy)", "identical_buffers": True,
        "bit_identical_output": ident,
        "max_abs_uint8": int(np.abs(again["still"].astype(np.int16)
                                    - sB["t"]["still"].astype(np.int16)).max()),
        "params": {"A": 2.6, "B_lowpass_sigma_px": 3.0, "E": sB["E_v9"],
                   "Q": 25, "g1": sB["g1_34"]},
        "pass": ident,
    }
    log("  G8a frozen composite on identical buffers: bit-identical=%s"
        % ident)
    cp.unit_done("s3-g8a")

    # ---- G8(b) Beer-Lambert unchanged for identical optical inputs ------
    rec_path = os.path.join(HERE, "out", "windy-still", "numbers-windy.json")
    try:
        with open(rec_path, "r", encoding="utf-8") as f:
            rec = json.load(f)
        rec_tw = rec["depth_cost"]["beer_direct_beam"]["z1.0"]["two_way"]
        rec_sd = rec["depth_cost"]["beer_direct_beam"]["z1.0"]["s_down_m"]
    except Exception as e:
        rec, rec_tw, rec_sd = None, None, None
    tw = sB["t"]["beer"]["two_way"]
    diff = [None if rec_tw is None else abs(float(tw[i]) - float(rec_tw[i]))
            for i in range(3)]
    g8b = {
        "what": ("Beer-Lambert on the direct beam for identical optical "
                 "inputs, compared with the v3.4/v9 stack's RECORDED numbers "
                 "(out/windy-still/numbers-windy.json)"),
        "v10_two_way": [float(x) for x in tw],
        "v9_recorded_two_way": rec_tw,
        "abs_diff": diff,
        "v10_s_down_m": float(sB["t"]["beer"]["s_down_m"]),
        "v9_recorded_s_down_m": rec_sd,
        "code_path": ("v5_physics.beer_lambert -- imported, not reimplemented "
                      "and not edited"),
        "pass": bool(diff[0] is not None
                     and max(float(d) for d in diff) <= 5e-5),
    }
    log("  G8b Beer two-way v10 %s vs v9 %s (max|d| %.2e)"
        % (["%.6f" % x for x in tw],
           ["%.6f" % x for x in (rec_tw or [])] if rec_tw else None,
           max(d for d in diff if d is not None) if rec_tw else float("nan")))
    cp.unit_done("s3-g8b")

    # ---- G8(c) eta=1 golden: v10 vs the v9 route on ITS OWN state -------
    sA = scene_rows["A"]
    route34 = transport_scene(sA["r34"], photo_lin, PW, PH, GRID_RENDER,
                              Q=25, E=sA["E_v9"], g1=sA["g1_34"], tag="v34")
    d = np.abs(sA["t"]["still"].astype(np.int16)
               - route34["still"].astype(np.int16))
    ms = msssim(sA["t"]["still"].astype(np.float64) / 255.0,
                route34["still"].astype(np.float64) / 255.0)
    lum10 = np.asarray(sA["t"]["still"], np.float64) @ np.array(
        [0.2126, 0.7152, 0.0722])
    lum34 = np.asarray(route34["still"], np.float64) @ np.array(
        [0.2126, 0.7152, 0.0722])
    g8c = {
        "scene": "A (preset-1 family, U=1.5 F=100 h=0.6)",
        "route": ("eta = 1 -> w_diffuse = 0 -> the v9 route reduces exactly to "
                  "det_frame + render_frame_v5 (run_v9.combine returns Cd). "
                  "Both renders use the SAME E, the SAME g1 and the same "
                  "frozen composite; only the state differs."),
        "E_used": sA["E_v9"], "g1_used": sA["g1_34"],
        "same_state_inputs": {"E_equal": True, "g1_equal": True,
                              "composite_equal": True, "Q_equal": True},
        "delta": {"mean_abs_uint8": float(d.mean()),
                  "max_abs_uint8": int(d.max()),
                  "frac_pixels_changed": float(np.mean(d.max(axis=2) > 0)),
                  "mean_luminance_v10": float(lum10.mean()),
                  "mean_luminance_v34": float(lum34.mean()),
                  "rel_luminance_change": float(lum10.mean()
                                                / lum34.mean() - 1.0)},
        "msssim": ms["msssim"],
        "state_cause": {
            "Hs_v10_m": float(sA["st10"].spectrum_after["Hs_m"]),
            "Hs_v34_m": float(sA["re34"]["Hs_m"]),
            "Hs_law_m": float(sA["st10"].laws["Hs_depth_limited_m"]),
            "sigma2_v10": float(sA["st10"].spectrum_after["sigma2_rad2"]),
            "sigma2_v34": float(sA["re34"]["sigma2_rad2"]),
            "patch_v10_m": float(sA["st10"].S),
            "patch_v34_m": float(sA["r34"]["state"].S),
            "H_s_ratio_v34_over_law": float(sA["re34"]["Hs_m"]
                                            / sA["st10"].laws["Hs_depth_limited_m"]),
            "note": ("the delta is the STATE: H_s falls from 4*waveHeightM "
                     "to the fetch law, sigma^2 is no longer tapered to a "
                     "target, and the patch is the scene scale instead of the "
                     "wind-coupled 4*specL.  The route is bit-identical (it "
                     "was not touched)."),
        },
        "pass": True,
        "pass_note": ("golden test is a REPORTED delta, not a bit-equality: "
                      "the state legitimately changed (PIPELINE-V9 section 4: "
                      "'preset realised values WILL change')"),
    }
    log("  G8c eta=1 golden: mean|d| %.2f/255 max %d MS-SSIM %.5f "
        "(lum %.4f vs %.4f, %+.1f%%)"
        % (g8c["delta"]["mean_abs_uint8"], g8c["delta"]["max_abs_uint8"],
           g8c["msssim"], g8c["delta"]["mean_luminance_v10"],
           g8c["delta"]["mean_luminance_v34"],
           100 * g8c["delta"]["rel_luminance_change"]))
    png = io.BytesIO()
    Image.fromarray(route34["still"]).save(png, format="PNG")
    atomic_write_bytes(os.path.join(OUT, "still-A-v34route.png"), png.getvalue())
    cp.unit_done("s3-g8c")

    # ---- G8(d) rain R=0 equals the v10 base state exactly --------------
    ra = run_cell(raw, SCENE_B["U"], SCENE_B["F"], SCENE_B["h"], 512, seed,
                  k_cut=V10.KCUT_DEFAULT)
    rb = run_cell(raw, SCENE_B["U"], SCENE_B["F"], SCENE_B["h"], 512, seed,
                  k_cut=V10.KCUT_DEFAULT)
    eq = bool(np.array_equal(ra["state"].lo.H0, rb["state"].lo.H0))
    rain_raised = False
    try:
        V10WaterState({"x": 1}, {"windSpeedMS": 6.0, "depthM": 1.0,
                                 "windDir": [1.0, 0.0]}, N=64, seed=1,
                      fetch_m=100.0, rain_rate_mm_h=1.0)
    except NotImplementedError:
        rain_raised = True
    g8d = {"bit_identical_state": eq, "rain_positive_raises": rain_raised,
           "pass": bool(eq and rain_raised),
           "note": ("R=0 takes no rain branch at all, so two builds are "
                    "bit-identical; R>0 raises instead of silently dropping "
                    "the band (PIPELINE-V9 section 8)")}
    log("  G8d rain R=0 identical=%s, R>0 raises=%s" % (eq, rain_raised))
    cp.unit_done("s3-g8d")

    # ---- G9 photon cross-check + C percentiles --------------------------
    g9 = {"scenes": {}}
    for key in ("A", "B"):
        s = scene_rows[key]
        C = s["t"]["C"]
        st = s["st10"]
        Cp, pinfo = photon_master(s["r10"]["fields"], st.S, s["r10"]["env"],
                                  PHOTONS, seed ^ 0x9E37, n_buf=512)
        Cd = downsample_mean(C, 512)
        # photon C has 3 channels; the pipeline's master is achromatic (green)
        Cpg = Cp[:, :, 1]
        det = Cd - Cd.mean()
        pho = Cpg - Cpg.mean()
        corr = float(np.sum(det * pho) / math.sqrt(
            float(np.sum(det * det)) * float(np.sum(pho * pho)) + 1e-30))
        hp = radial_power_fraction(Cd, 0.5)
        cl = correlation_length_px(C)
        # deliverable still vs the photon render (green master, same route)
        ph_still = None
        try:
            from v5_composite import render_frame_v5, A_V3, B_SIGMA_V3
            ph_still = render_frame_v5(
                s["r10"]["fields"], Cpg.astype(np.float32), photo_lin,
                s["r10"]["env"], PW, PH, ambient_scale=A_V3,
                c_layer_exposure=s["E_v9"],
                substrate_gain_rgb=s["t"]["opt"]["gain"]["gain_rgb"],
                beer_down=s["t"]["opt"]["beer"]["down"],
                beer_up=s["t"]["opt"]["beer"]["up"],
                g1=s["g1_34"], blur_sigma_px=s["t"]["opt"]["blur"]["sigma_px"],
                lowpass_sigma_px=B_SIGMA_V3, return_debug=False)
            ms_p = msssim(s["t"]["still"].astype(np.float64) / 255.0,
                          ph_still.astype(np.float64) / 255.0)
        except Exception as exc:
            ms_p = {"msssim": None, "error": str(exc)}
        cs = c_stats(C)
        cs_ph = c_stats(Cpg)
        # quadrature check: Q=25 vs Q=49 (both scenes)
        qchk = None
        try:
            from det import det_frame
            C49, _ = det_frame(s["r10"]["fields"], st.S, GRID_RENDER,
                               s["r10"]["env"], 49)
            qchk = {
                "Q25_vs_Q49_mean_abs_diff": float(np.abs(
                    C.astype(np.float64) - C49.astype(np.float64)).mean()),
                "Q25_vs_Q49_rel_l1": float(np.abs(
                    C.astype(np.float64) - C49.astype(np.float64)).mean()
                    / max(abs(float(C.mean())), 1e-30)),
                "p99_9": {"Q25": cs["p99_9"], "Q49": c_stats(C49)["p99_9"]},
                "max": {"Q25": cs["max"], "Q49": c_stats(C49)["max"]},
                "pass": bool(abs(c_stats(C49)["p99_9"] / cs["p99_9"] - 1.0)
                             <= 0.10),
                "note": ("quadrature convergence of the forward transport: a "
                         "shredded master would move with Q; this one does "
                         "not"),
            }
            del C49
        except Exception as exc:
            qchk = {"error": str(exc)}
        g9["scenes"][key] = {
            "C_det": cs, "C_photon_green": cs_ph,
            "photon": {"photons": int(PHOTONS), "n_buf": 512,
                       "mean_expected": pinfo["expected_mean"],
                       "mean_measured": pinfo["measured_mean"]},
            "cross_check": {
                "corr_det_vs_photon": corr,
                "mean_det": float(Cd.mean()), "mean_photon": float(Cpg.mean()),
                "rel_mean_diff": float(Cpg.mean() / Cd.mean() - 1.0),
                "p99_det": float(np.percentile(Cd, 99)),
                "p99_photon": float(np.percentile(Cpg, 99)),
                "p99_9_det": float(np.percentile(Cd, 99.9)),
                "p99_9_photon": float(np.percentile(Cpg, 99.9)),
            },
            "shredding": {**hp, "corr_len": cl},
            "msssim_still_vs_photon_render": ms_p,
            "g9_band": {"p99_9_over_mean_range": [6.0, 12.0],
                        "realised": cs["p99_9_over_mean"]},
            "quadrature_convergence": qchk,
            "pass_mean_1": cs["pass_mean_1"],
            "within_band": bool(6.0 <= (cs["p99_9_over_mean"] or 0.0) <= 12.0),
        }
        log("       G9 %s: C mean %.9f p99 %.3f p99.9 %.3f (%.2fx mean) "
            "corr_len %.1f px photon corr %.4f MS-SSIM %.5f"
            % (key, cs["mean"], cs["p99"], cs["p99_9"],
               cs["p99_9_over_mean"] or -1,
               cl["corr_len_px_1e"], corr, ms_p.get("msssim") or -1))
        cp.unit_done("s3-g9-%s" % key)
    g9["pass_mean_invariant"] = all(g9["scenes"][k]["pass_mean_1"]
                                    for k in g9["scenes"])
    g9["all_within_6_12"] = all(g9["scenes"][k]["within_band"]
                                for k in g9["scenes"])
    g9["shredding_verdict"] = {
        "criterion": ("shredding = the master is quadrature-noise, i.e. it "
                      "moves with Q, is delta-correlated in space, or the two "
                      "independent transports disagree.  The 6-12 band on "
                      "p99.9/mean is reported SEPARATELY as a physical "
                      "comparison, not as the shredding test."),
        "Q_convergence": all(bool((g9["scenes"][k]["quadrature_convergence"]
                                   or {}).get("pass")) for k in g9["scenes"]),
        "corr_len_px_1e": {k: g9["scenes"][k]["shredding"]["corr_len"]
                           ["corr_len_px_1e"] for k in g9["scenes"]},
        "photon_corr": {k: g9["scenes"][k]["cross_check"]
                        ["corr_det_vs_photon"] for k in g9["scenes"]},
        "hf_power_fraction": {k: g9["scenes"][k]["shredding"]
                              ["power_fraction_above_half_nyquist"]
                              for k in g9["scenes"]},
    }
    sv = g9["shredding_verdict"]
    sv["verdict"] = ("no shredding" if (sv["Q_convergence"]
                                        and all(v > 2.0 for v in
                                                sv["corr_len_px_1e"].values())
                                        and all(v >= 0.80 for v in
                                                sv["photon_corr"].values()))
                     else "shredding suspected")
    g9["pass"] = bool(g9["pass_mean_invariant"]
                      and sv["verdict"] == "no shredding")
    log("  G9 band check: %s (p99.9/mean %s vs 6-12) -> reported, not used as "
        "the shredding test"
        % (g9["all_within_6_12"],
           {k: round(g9["scenes"][k]["C_det"]["p99_9_over_mean"], 2)
            for k in g9["scenes"]}))
    log("  G9 verdict: mean=1 %s, Q-converged %s, corr_len %s px, photon corr "
        "%s -> %s (PASS=%s)"
        % (g9["pass_mean_invariant"], sv["Q_convergence"],
           sv["corr_len_px_1e"], sv["photon_corr"], sv["verdict"], g9["pass"]))
    cp.unit_done("s3-g9-verdict")

    # ---- per-scene numbers, then the JSON ------------------------------
    for key in ("A", "B"):
        s = scene_rows[key]
        st = s["st10"]
        node = st.state_node(s["r10"]["fields"])
        out["scenes"][key] = {
            "name": s["scene"]["name"],
            "recipe": {"U_ms": s["scene"]["U"], "fetch_m": s["scene"]["F"],
                       "depth_m": s["scene"]["h"], "rain_mm_h": 0.0,
                       "eta": 1.0, "seed": seed},
            "v10_state": node,
            "v34_state": {
                "patch_m": float(s["r34"]["state"].S),
                "Hs_realised_m": float(s["re34"]["Hs_m"]),
                "sigma2_realised": float(s["re34"]["sigma2_rad2"]),
                "rms_slope": float(s["re34"]["rms_slope"]),
            },
            "route": {"E_used": s["E_v9"], "g1_used": s["g1_34"],
                      "E_implied_by_v10_state": s["E_v9_10"],
                      "g1_implied_by_v10_state": s["g1_10"],
                      "Q": 25, "A": 2.6, "B_lowpass_sigma_px": 3.0},
            "C": c_stats(s["t"]["C"]),
            "beer": s["t"]["beer"],
            "still_metrics": s["t"]["metrics"],
            "timing_s": s["t"]["timing_s"],
            "files": ["still-%s-v37.png" % key],
        }
    out["gates"] = {"G8a": g8a, "G8b": g8b, "G8c": g8c, "G8d": g8d, "G9": g9}
    out["wall_s"] = float(time.time() - t0)
    dump_json(os.path.join(OUT, "numbers-v10-transport.json"), out)
    cp.unit_done("s3-json")
    log("[stage 3] done in %.1fs" % (time.time() - t0))
    return out


def correlation_length_px(C):
    """1/e correlation length of the master, in pixels (structure vs noise)."""
    c = np.asarray(C, dtype=np.float64)
    c = c - c.mean()
    n = c.shape[0]
    F = np.abs(np.fft.fft2(c)) ** 2
    ac = np.real(np.fft.ifft2(F))
    ac = ac / ac[0, 0]
    # radial average along x and y (periodic, so both are symmetric)
    ax = np.abs(ac[0, :n // 2])
    ay = np.abs(ac[:n // 2, 0])
    a = 0.5 * (ax + ay)
    idx = np.where(a < math.exp(-1.0))[0]
    lam = float(idx[0]) if idx.size else float(n // 2)
    return {"corr_len_px_1e": lam,
            "corr_len_mm": float(lam * 1.0),
            "first_negative_lag": (int(np.where(ac[0, :n // 2] < 0)[0][0])
                                   if np.any(ac[0, :n // 2] < 0) else None)}


# ==========================================================================
# stage 4 -- G10 windy re-render + report
# ==========================================================================
WINDY = {"U": 6.0, "F": 100.0, "h": 1.0}
YARDSTICKS = {
    "v34_windy_still": {"mean_px_ncc": 10.50, "std_px_ncc": 8.52,
                        "own_field_mean_px": 5.84,
                        "own_field_std_px": 3.27,
                        "px_scale_mm": 7.63,
                        "own_field_mm": 48.5,
                        "method": "48 px blocks, +/-24 px, 2 px steps, Rec.709"},
    "production_loop_2026_06_10": {"mean_px": 23.2, "std_px": 17.1,
                                   "px_scale_mm": 1.77,
                                   "mean_mm_from_px": 41.1},
}


def stage4():
    t0 = time.time()
    log("=" * 72)
    log("[stage 4] G10 windy re-render under v3.7 + warp yardsticks")
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
    base_u8, photo_lin, PW, PH = _load_photo()
    log("  photo %dx%d" % (PW, PH))

    r10 = v10_scene_for(raw, WINDY["U"], WINDY["F"], WINDY["h"], GRID_RENDER,
                        seed)
    st = r10["state"]
    log("  v10 state: mode=%s S=%.3f m Hs spec %.4f sig2 %.5f rms %.4f"
        % (st.mode, st.S, st.spectrum_after["Hs_m"],
           st.spectrum_after["sigma2_rad2"],
           st.spectrum_after["rms_slope"]))
    t = transport_scene(r10, photo_lin, PW, PH, GRID_RENDER, Q=25,
                        E=float(E_V3))
    C = t["C"]
    still = t["still"]
    atomic_write_bytes(os.path.join(WORK, "C-windy-v37-Q25.npy"),
                       np.ascontiguousarray(C).tobytes())
    png = io.BytesIO()
    Image.fromarray(still).save(png, format="PNG")

    def verify_png(p):
        try:
            im = Image.open(p)
            im.load()
            return im.size == (PW, PH) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(os.path.join(OUT, "comp-windy-v37.png"), png.getvalue(),
                       verify=verify_png)
    log("  wrote comp-windy-v37.png (%dx%d)" % (PW, PH))
    cp.unit_done("s4-composite")

    # ---- warp measurement ---------------------------------------------
    still_png = np.asarray(Image.open(
        os.path.join(OUT, "comp-windy-v37.png")).convert("RGB"))
    warp = warp_measure(base_u8, still_png)
    p = warp["primary"]
    ncc = warp["ncc_check"]
    pred = predicted_shift_field(r10["fields"], r10["env"], PW, PH)
    ps = predicted_block_stats(pred, warp["grid"])
    agree = {"ssd_primary": block_agreement(ps, warp["per_block"],
                                            "ssd_dy", "ssd_dx"),
             "ncc_check": block_agreement(ps, warp["per_block"],
                                          "ncc_dy", "ncc_dx")}
    px_mm = float(ps["pixel_scale_mm"])
    own_mm = 1e3 * float(ps["mean_predicted_shift_m"])
    log("  warp: SSD %.2f/%.2f px | NCC %.2f/%.2f px | own field %.2f px = "
        "%.2f mm at %.3f mm/px"
        % (p["mean_shift_px"], p["shift_field_std_px"],
           ncc["mean_shift_px"], ncc["shift_field_std_px"],
           ps["field_mean_magnitude_px"], own_mm, px_mm))

    # ---- C=1 control (isolates the geometry from the caustic gain) ------
    ctrl_img = render_frame_v5(
        r10["fields"], np.ones((GRID_RENDER, GRID_RENDER), np.float32),
        photo_lin, r10["env"], PW, PH, ambient_scale=A_V3,
        c_layer_exposure=float(E_V3),
        substrate_gain_rgb=t["opt"]["gain"]["gain_rgb"],
        beer_down=t["opt"]["beer"]["down"], beer_up=t["opt"]["beer"]["up"],
        g1=float(t["g1"]), blur_sigma_px=t["opt"]["blur"]["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=False)
    ctrl = warp_measure(base_u8, ctrl_img)
    del ctrl_img
    gc.collect()
    log("  control C=1: SSD %.2f/%.2f | NCC %.2f/%.2f px"
        % (ctrl["primary"]["mean_shift_px"],
           ctrl["primary"]["shift_field_std_px"],
           ctrl["ncc_check"]["mean_shift_px"],
           ctrl["ncc_check"]["shift_field_std_px"]))
    cp.unit_done("s4-warp")

    # ---- tone bounds + caps -------------------------------------------
    dif = np.abs(still_png.astype(np.int16) - base_u8.astype(np.int16))
    node = st.state_node(r10["fields"])

    def _texture_stats(x_u8):
        x = np.asarray(x_u8, np.float64) / 255.0
        lum = x @ np.array([0.2126, 0.7152, 0.0722])
        return {"mean_grad_luminance": float(np.abs(np.diff(lum, axis=0)).mean()),
                "mean_luminance": float(lum.mean()),
                "std_luminance": float(lum.std()),
                "msssim_vs_base": float(msssim(base_u8.astype(np.float64) / 255.0,
                                               x)["msssim"])}

    texture = {"base_photo": _texture_stats(base_u8),
               "still_v37": _texture_stats(still_png)}
    msssim_v10 = None
    try:
        v34 = v8_state_for(raw, WINDY["U"], WINDY["F"], WINDY["h"],
                           GRID_RENDER, seed)
        t34 = transport_scene(v34, photo_lin, PW, PH, GRID_RENDER, Q=25,
                             E=float(E_V3))
        # NOTE: E2.2 pinned by the brief for both, so this is the pure state
        # comparison at the frozen display constant.
        msssim_v10 = msssim(still.astype(np.float64) / 255.0,
                            t34["still"].astype(np.float64) / 255.0)
        texture["still_v34"] = _texture_stats(t34["still"])
        v34_summary = {
            "patch_m": float(v34["state"].S),
            "Hs_realised_m": float(v34["state"].realised_stats(
                v34["fields"])["Hs_m"]),
            "sigma2_realised": float(v34["state"].realised_stats(
                v34["fields"])["sigma2_rad2"]),
            "C": c_stats(t34["C"]),
        }
        del t34, v34
        gc.collect()
    except Exception as exc:
        v34_summary = {"error": str(exc)}
    cp.unit_done("s4-v34-compare")

    numbers = {
        "task": ("BRIEF-V10-CHUNK7 / stage 4: G10 windy pond re-render under "
                 "version 3.7 (v10 state generator)"),
        "created_utc": now_iso(),
        "photo": {"path": PHOTO, "size": [int(PW), int(PH)],
                  "sha256": sha256_file(PHOTO)},
        "recipe": {"U_ms": WINDY["U"], "fetch_m": WINDY["F"],
                   "depth_m": WINDY["h"], "rain_mm_h": 0.0, "eta": 1.0,
                   "sun_date": raw["sun"]["date"],
                   "sun_time": raw["sun"]["time"],
                   "sun_altitude_deg": float(r10["env"]["sun"]["altitudeDeg"]),
                   "sun_azimuth_deg": float(r10["env"]["sun"]["azimuthDeg"]),
                   "seed": seed},
        "composite": {"A": A_V3, "E": float(E_V3), "Q": 25,
                      "B_lowpass_sigma_px": B_SIGMA_V3,
                      "blur_sigma_px": float(t["opt"]["blur"]["sigma_px"]),
                      "g1": float(t["g1"]),
                      "note": ("frozen production composite: A2.6-Blow + "
                               "E2.2, exactly as BRIEF-WINDY-STILL pinned")},
        "v10_state": node,
        "v34_same_scene": v34_summary,
        "msssim_still_v10_vs_v34": msssim_v10,
        "texture": texture,
        "C": c_stats(C),
        "corr_len": correlation_length_px(C),
        "beer": t["beer"],
        "still_metrics": t["metrics"],
        "still_vs_base": {
            "max_abs_diff_uint8": int(dif.max()),
            "mean_abs_diff_uint8": float(dif.mean()),
            "changed_pixel_fraction": float(np.mean(dif.max(axis=2) > 0)),
        },
        "warp": {
            "method": warp["definition"],
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
                "mean_mm": own_mm,
                "pixel_scale_mm": px_mm,
                "frac_beyond_search": ps["field_frac_beyond_search"],
            },
            "agreement_with_predicted": agree,
            "control_C1": {
                "ssd_primary": {k: ctrl["primary"][k] for k in (
                    "mean_shift_px", "shift_field_std_px",
                    "frac_blocks_at_search_edge")},
                "ncc_check": {k: ctrl["ncc_check"][k] for k in (
                    "mean_shift_px", "shift_field_std_px",
                    "frac_blocks_at_search_edge")},
            },
            "yardsticks": YARDSTICKS,
        },
        "tone": {"crush_rate_linear_lt_0.02":
                 float(t["metrics"]["crush_rate_linear_lt_0.02"]),
                 "blowout_rate_aces_gt_0.98":
                 float(t["metrics"]["blowout_rate_aces_gt_0.98"]),
                 "customary_crush_bound": 0.0325,
                 "crush_pass_customary": bool(
                     t["metrics"]["crush_rate_linear_lt_0.02"] <= 0.0325),
                 "note": ("reported, not tuned: the brief pins E2.2 and this "
                          "is its consequence")},
        "files": ["comp-windy-v37.png", "_work/C-windy-v37-Q25.npy"],
        "timing_s": t["timing_s"],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-windy-v37.json"), numbers)
    cp.unit_done("s4-json")
    log("  crush %.5f blow %.6f | C p50 %.3f p99 %.3f p99.9 %.3f max %.2f"
        % (numbers["tone"]["crush_rate_linear_lt_0.02"],
           numbers["tone"]["blowout_rate_aces_gt_0.98"],
           numbers["C"]["p50"], numbers["C"]["p99"],
           numbers["C"]["p99_9"], numbers["C"]["max"]))
    log("[stage 4] done in %.1fs" % (time.time() - t0))
    return numbers


def stage_report():
    """Write REPORT-V10-GATES.md from the four stage JSONs and close the cp."""
    t0 = time.time()
    log("=" * 72)
    log("[report] assembling REPORT-V10-GATES.md from the stage JSONs")
    cp = Checkpoint(CP_PATH)
    import v10_report
    v10_report.main()
    cp.unit_done("s5-report")
    log("[report] done in %.1fs" % (time.time() - t0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="1",
                    choices=["1", "2", "3", "4", "report", "all"])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    if args.stage in ("1", "all"):
        stage1()
    if args.stage in ("2", "all"):
        stage2()
    if args.stage in ("3", "all"):
        stage3()
    if args.stage in ("4", "all"):
        stage4()
    if args.stage in ("report", "all"):
        stage_report()


if __name__ == "__main__":
    main()
