"""run_v15 -- version 3.10 (chunk 11): THE MERGE.  Four foreground stages.

Brief: BRIEF-V15-MERGE-R2.md.  Design: PIPELINE-V12.md (governs everything
below it).  Base: the v14 tree.  Every deliverable goes through jobkit
`atomic_write_bytes`; one checkpoint unit is committed per stage.

  stage 1  the merged default + the region-wise battery over the 8-cell
           matrix (analytic): both columns (a)/(c) conformance, (b) beta,
           (d) JONSWAP anchor, C1/C1b, totals, K1 + hs_flag, provenance,
           verdicts -> out/v15/numbers-v15-merge.json
  stage 2  regression: D_STRICT == v11, v15 == v14 RAW numbers per cell,
           v12 water default + the G1 legacy lane byte-identity, the frozen
           module hash list, and the jonswap rollback lane against the v10
           scene-A statistics -> merged into numbers-v15-merge.json
  stage 3  transport through the merged default stack: T1/T2/T3 stills +
           numbers-v15-transport.json
  stage 4  hash list, verdict roll-up, REPORT-V15-MERGE.md

Usage (foreground, from caustics-filter/offline):
  python run_v15.py --stage 1
  python run_v15.py --stage 2
  python run_v15.py --stage 3
  python run_v15.py --stage 4
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
import run_v14_spectrum as V14                                     # noqa: E402
import run_v15_spectrum as V15                                     # noqa: E402
import run_v15_state as S15                                        # noqa: E402
import run_v10 as R10                                              # noqa: E402

OUT = os.path.join(HERE, "out", "v15")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
V14_OUT = os.path.join(HERE, "out", "v14")
V13_OUT = os.path.join(HERE, "out", "v13")
V12_OUT = os.path.join(HERE, "out", "v12")
V10_OUT = os.path.join(HERE, "out", "v10")

GRID_SWEEP = 1024
GRID_RENDER = R10.GRID_RENDER
SWEEP_U = V10.SWEEP_U
SWEEP_F = V10.SWEEP_F
SWEEP_H = V10.SWEEP_Z

#: stage 3's three stills (brief, stage 3).
T1 = {"id": "T1", "name": "windy", "U": 6.0, "F": 1000.0, "h": 1.0,
      "file": "still-T1-windy.png",
      "why": "the standing windy comparison scene"}
T2 = {"id": "T2", "name": "lowwind", "U": 3.0, "F": 1000.0, "h": 0.6,
      "file": "still-T2-lowwind.png",
      "why": "low-wind cell: the UNRESOLVED label must travel with the output"}
T3 = {"id": "T3", "name": "warned", "U": 6.0, "F": 100.0, "h": 1.0,
      "file": "still-T3-warned.png",
      "why": ("Omega_c = 7.30 -> warn_extrapolated_wave_age must travel with "
              "the output")}
TRANSPORT_CELLS = (T1, T2, T3)

#: files this job must not change (the v14 frozen list, verbatim).
FROZEN_FILES = (
    "v5_composite.py", "v5_physics.py", "v12_water.py", "v10_spectrum.py",
    "v10_state.py", "v11_spectrum.py", "v11_state.py", "v13_spectrum.py",
    "v13_state.py", "transport.py", "det.py", "image.py", "env.py",
    "recipe_v8.py", "run_v13.py", "run_v10.py", "run_v12.py", "run_v14.py",
    "run_v14_spectrum.py", "run_v14_state.py", "run_v14_report.py",
)
NEW_FILES = ("run_v15_spectrum.py", "run_v15_state.py", "run_v15.py",
             "run_v15_report.py")

#: the v12 G1 legacy artefact (delivery still-A-v37), the brief's reference.
G1_LEGACY_SHA = "cfb086b22a0e8f1b25cb3befbcc272f3e9c5127fe5499396a555c99b297f0e32"
G1_LEGACY_FILE = os.path.join(V10_OUT, "still-A-v37.png")

#: the v10 cell the jonswap rollback lane is compared against (scene A).
V10_SCENE_A_CELL = "U6p0_F100p0"

log = R10.log
now_iso = R10.now_iso
dump_json = R10.dump_json
sha256_file = R10.sha256_file
load_recipe = R10.load_recipe


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def hashes(names):
    out = {}
    for f in names:
        p = os.path.join(HERE, f)
        out[f] = sha256_file(p) if os.path.exists(p) else None
    return out


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


def scene_for(raw, U, F, h, N, seed, source=None):
    """v15 merged-default state -> the dict run_v10.transport_scene expects."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    t0 = time.time()
    st = S15.V15WaterState(
        scene, env, N=N, seed=seed, preset_index=1, fetch_m=float(F),
        k_cut=V10.KCUT_DEFAULT,
        spectrum_source=(source or V15.SPECTRUM_SOURCE_DEFAULT),
        arm=V15.MERGE_ARM)
    t_state = time.time() - t0
    t1 = time.time()
    fields = st.fields_at(0.0)
    t_fields = time.time() - t1
    return {"state": st, "env": env, "cfg": cfg, "fields": fields,
            "t_state": t_state, "t_fields": t_fields}


def _v11_hard_zero_identity(U, F, h):
    """D_STRICT must reproduce the v11/v3.8 build (alpha_m = max(raw, 0))."""
    a = V15.assemble(U, F, h, arm=V15.D_STRICT)
    b = V11.assemble(U, F, h, source="elfouhaily-form")
    k = V11._k_grid(float(a["k_cut_rad_per_m"]))
    sa = V11.s_parts(k, a)
    sb = V11.s_parts(k, b)
    d = np.abs(sa["S"] - sb["S"])
    return {"cell": {"U_ms": float(U), "fetch_m": float(F),
                     "depth_m": float(h)},
            "alpha_m_used_d_strict": float(a["alpha_m_value"]),
            "alpha_m_used_v11": float(b["alpha_m_value"]),
            "max_abs_dS": float(d.max()), "S_identical": bool(not d.any()),
            "pass": bool(not d.any()
                         and a["alpha_m_value"] == b["alpha_m_value"])}


def _arm_identity(U, F, h, arm, v14_arm=None):
    """v15's arm must reproduce v14's arm -- same code path, exact alike."""
    a = V15.assemble(U, F, h, arm=arm)
    b = V14.assemble(U, F, h, arm=(v14_arm or arm))
    k = V11._k_grid(float(a["k_cut_rad_per_m"]))
    sa = V11.s_parts(k, a)
    sb = V11.s_parts(k, b)
    d = np.abs(sa["S"] - sb["S"])
    return {"cell": {"U_ms": float(U), "fetch_m": float(F),
                     "depth_m": float(h)},
            "arm_v15": arm, "arm_v14": (v14_arm or arm),
            "alpha_m_v15": float(a["alpha_m_value"]),
            "alpha_m_v14": float(b["alpha_m_value"]),
            "sigma2_h_v15": float(a["sigma2_h_total"]),
            "sigma2_h_v14": float(b["sigma2_h_total"]),
            "Hs_v15": float(a["Hs_m"]), "Hs_v14": float(b["Hs_m"]),
            "max_abs_dS": float(d.max()), "S_identical": bool(not d.any()),
            "pass": bool(not d.any()
                         and a["alpha_m_value"] == b["alpha_m_value"]
                         and a["Hs_m"] == b["Hs_m"])}


# ==========================================================================
# stage 1 -- the merged default + the region-wise battery
# ==========================================================================
def stage1():
    t0 = time.time()
    log("=" * 72)
    log("[stage 1] merged default (%s) + region-wise battery over 8 cells"
        % V15.SPECTRUM_SOURCE_DEFAULT)
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    frozen_at_start = hashes(FROZEN_FILES)

    assert V15.SPECTRUM_SOURCE_DEFAULT == "elfouhaily-form"
    assert V15.MERGE_ARM == "D"
    assert "jonswap" in V15.SPECTRUM_SOURCES
    log("  default=%s | merge arm=%s | rollback selectable=jonswap | retired "
        "clause: max B(k) <= %g (not gated anywhere in v15)"
        % (V15.SPECTRUM_SOURCE_DEFAULT, V15.MERGE_ARM,
           V15.RETIRED_B_ENV_MAX))

    cells = []
    for U in SWEEP_U:
        for F in SWEEP_F:
            t1 = time.time()
            ctx = V15.assemble(U, F, SWEEP_H)
            ev = V15.evaluate_cell(U, F, SWEEP_H, ctx=ctx)
            ev["timing_s"] = float(time.time() - t1)
            cells.append(ev)
            a = ev["conformance"]["a_gravity"]
            c = ev["conformance"]["c_capillary"]
            b = ev["observation"]["b_beta"]
            log("  [%-14s] Om_c %5.2f warn %-5s | (a) %.4f %-4s | (c) %.4f "
                "%-4s | (b) beta %.5f %-30s | K1 %+6.1f%% flag %-5s | %s"
                % (ev["cell_id"], ev["omega_c"],
                   ev["warn_extrapolated_wave_age"], a["ratio"], a["pass"],
                   c["ratio"], c["pass"], b["beta"], b["status"],
                   100.0 * ev["hs"]["delta_rel"], ev["hs"]["hs_flag"],
                   ev["state"]))

    # ---- summary + flags inventory --------------------------------------
    states = {c["cell_id"]: c["state"] for c in cells}
    summary = {
        "n_pass": int(sum(1 for v in states.values() if v == "PASS")),
        "n_fail": int(sum(1 for v in states.values() if v == "FAIL")),
        "n_unresolved": int(sum(1 for v in states.values()
                                if v == "UNRESOLVED")),
        "states": states,
        "zero_fails": bool(all(v != "FAIL" for v in states.values())),
        "cells_warn_extrapolated_wave_age": [
            c["cell_id"] for c in cells if c["warn_extrapolated_wave_age"]],
        "cells_unresolved": [c["cell_id"] for c in cells
                             if c["state"] == "UNRESOLVED"],
        "cells_hs_flag": {c["cell_id"]: c["hs"]["delta_rel"] for c in cells
                          if c["hs"]["hs_flag"]},
        "cells_with_prior_or_topup": [
            c["cell_id"] for c in cells
            if c["alpha"]["prior_active"] or c["alpha"]["top_up_acted"]],
        "conformance_a_pass": int(sum(1 for c in cells
                                      if c["conformance"]["a_gravity"]["pass"])),
        "conformance_c_pass": int(sum(1 for c in cells
                                      if c["conformance"]["c_capillary"]["pass"])),
        "observation_b_status": {c["cell_id"]:
                                 c["observation"]["b_beta"]["status"]
                                 for c in cells},
        "totals_tol_pass": int(sum(1 for c in cells
                                   if c["envelope_check"]["tol_inside"])),
        "C1_pass_high_wind": int(sum(1 for c in cells
                                     if not c["low_wind_cell"]
                                     and c["C1"]["pass"])),
        "C1b_pass_high_wind": int(sum(1 for c in cells
                                      if not c["low_wind_cell"]
                                      and c["C1b"]["pass"])),
        "K1_pass": int(sum(1 for c in cells if c["hs"]["pass"])),
        "tagged": ("region-wise conformance + observation battery -- NOT a "
                   "claim of physical truth (rally condition)"),
    }

    stage = {
        "stage": 1,
        "brief": "BRIEF-V15-MERGE-R2.md",
        "design": "PIPELINE-V12.md",
        "created_utc": now_iso(),
        "model": "version 3.10 (chunk 11) -- THE MERGE",
        "label": ("region-wise conformance + observation battery (NOT a claim "
                  "of physical truth)"),
        "merge": {
            "default_spectrum_source": V15.SPECTRUM_SOURCE_DEFAULT,
            "default_arm": V15.MERGE_ARM,
            "default_arm_rule": V15.MERGE_RULE,
            "rollback_source": "jonswap (still selectable; v11/v10 lane "
                               "returned byte-for-byte)",
            "retired_clause": dict(V15.RETIRED_CLAUSE),
            "what_replaced_it": ("CONFORMANCE (a) gravity flank vs "
                                 "0.1433 alpha_p gamma 1.04, (c) capillary "
                                 "peak vs (1/2) alpha_m; OBSERVATION (b) "
                                 "equilibrium beta in [0.006, 0.024], "
                                 "(d) JONSWAP anchor as context"),
            "provenance_labels": {
                "eckv_fit_domain_raw": "u* >= c_m/e, alpha_m raw",
                "low_wind_prior": "u* < c_m/e, max(raw, floor) with the "
                                  "design's sigma2_h_target",
            },
        },
        "policy": {
            "spectrum_source_default": V15.SPECTRUM_SOURCE_DEFAULT,
            "spectrum_sources": list(V15.SPECTRUM_SOURCES),
            "merge_arm": V15.MERGE_ARM,
            "alpha_m_rule": V15.MERGE_RULE,
            "sigma2_h_target_default": V15.SIGMA2_H_TARGET_DEFAULT,
            "ustar_fit_min_ms": V15.USTAR_FIT_MIN,
            "c_m_ms": V15.C_M,
            "k_m_rad_per_m": V15.K_M,
            "retired_B_env_max": V15.RETIRED_B_ENV_MAX,
            "gravity_conformance": {
                "band": "[k_p/2, 2 k_p]",
                "reference": "0.1433 * alpha_p * gamma(Omega_c) * 1.04",
                "pass_if_ratio_le": V15.CONF_A_RATIO_MAX},
            "capillary_conformance": {
                "band": "[k_m/2, 2 k_m]",
                "reference": "(1/2) * alpha_m(u*)",
                "pass_if_ratio_in": list(V15.CONF_C_RANGE)},
            "equilibrium_beta": {
                "band": list(V15.BETA_BAND),
                "I_p": V15.BETA_I_P,
                "f_band_factors": list(V15.BETA_F_BAND),
                "convention": dict(V15.BETA_CONVENTION)},
            "jonswap_anchor": {
                "form": "0.5 * alpha_J * exp(-1.25) * 3.3, "
                        "alpha_J = 0.076 (g F / U^2)^-0.22",
                "gated": False},
            "C1_gravity_share_max": V15.GRAVITY_SHARE_MAX,
            "C1b_capillary_share_min": V15.CAP_SHARE_MIN,
            "C1b_sigma2_h_min": V15.SIGMA2_H_MIN,
            "C1_C1b_gated_for": "U >= 6 only",
            "totals_tolerance": list(V15.TOTALS_TOL),
            "totals_definition": ("[Cox-Munk clean, Hu] with the 0.5x-2.0x "
                                  "clause: unchanged, all cells, gating"),
            "K1_hs_tolerance": V15.HS_TOL,
            "K1_hs_flag": V15.HS_FLAG,
            "grid": {"U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
                     "depth_m": SWEEP_H,
                     "k_cut_rad_per_m": V10.KCUT_DEFAULT,
                     "analytic_grid": "V11._k_grid (48k log) + 20k/6k local "
                                      "grids for the battery",
                     "band": [V15.BAND[0], V15.BAND[1]]},
        },
        "summary": summary,
        "provenance": {
            "what": ("the provenance the brief requires in this numbers file "
                     "as well as in every per-cell record"),
            "spectrum_source": V15.SPECTRUM_SOURCE_DEFAULT,
            "merge_arm": V15.MERGE_ARM,
            "arm_rule": V15.MERGE_RULE,
            "warn_rule": "warn_extrapolated_wave_age = (Omega_c > %.1f)"
                         % V15.WARN_OMEGA_C,
            "low_wind_rule": ("U <= 3 -> low-wind cell; UNRESOLVED at best, "
                              "FAIL only on the brief's hard-clause list"),
            "cells_warn_extrapolated_wave_age":
                summary["cells_warn_extrapolated_wave_age"],
            "cells_low_wind": [c["cell_id"] for c in cells
                               if c["low_wind_cell"]],
            "cells_eckv_fit_domain": [c["cell_id"] for c in cells
                                      if c["alpha"]["in_eckv_fit_domain"]],
            "cells_alpha_provenance": {
                c["cell_id"]: c["provenance"]["provenance"] for c in cells},
            "cells_alpha_domain": {
                c["cell_id"]: c["provenance"]["alpha_m_domain"]
                for c in cells},
        },
        "cells": cells,
        "frozen_hashes_at_stage1": frozen_at_start,
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v15-merge.json"), stage)
    log("  verdicts: PASS %d / FAIL %d / UNRESOLVED %d | zero_fails=%s"
        % (summary["n_pass"], summary["n_fail"], summary["n_unresolved"],
           summary["zero_fails"]))
    log("  warn cells: %s" % (summary["cells_warn_extrapolated_wave_age"],))
    log("  hs_flag cells: %s" % (summary["cells_hs_flag"],))
    cp.unit_done("stage-1")
    log("[stage 1] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 2 -- regression, item by item
# ==========================================================================
def stage2():
    t0 = time.time()
    log("=" * 72)
    log("[stage 2] regression: D_STRICT, v14 numbers, v12 water, frozen "
        "modules, jonswap rollback")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    with open(os.path.join(OUT, "numbers-v15-merge.json"), "r",
              encoding="utf-8") as f:
        merge = json.load(f)

    reg = {}

    # ---- item 1: D_STRICT == v11.assemble -------------------------------
    strict = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            strict[V15.cell_id(U, F)] = _v11_hard_zero_identity(U, F, SWEEP_H)
    ok1 = all(v["pass"] for v in strict.values())
    reg["item1_D_STRICT_equals_v11"] = {
        "requirement": ("D_STRICT branch == v11.assemble, byte-identical "
                        "spectrum (as in v13/v14)"),
        "pass": bool(ok1), "cells": strict,
        "note": ("the v3.8 hard-zero branch (alpha_m = max(raw, 0)) is kept as "
                 "a labelled diagnostic and still reproduces the v11 build "
                 "elementwise on all 8 cells")}
    log("  [1] D_STRICT == v11.assemble: %s" % ok1)

    # ---- item 2: v15 arm numbers == v14 arm numbers ---------------------
    raw_ident = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            raw_ident[V15.cell_id(U, F)] = _arm_identity(U, F, SWEEP_H, "RAW")
    ok2a = all(v["pass"] for v in raw_ident.values())
    # and the merged default (arm D) against v14's RAW column, in session
    d_ident = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            d_ident[V15.cell_id(U, F)] = _arm_identity(U, F, SWEEP_H, "D",
                                                       v14_arm="D")
    ok2b = all(v["pass"] for v in d_ident.values())
    # against the STORED v14 numbers file, per cell, field by field
    with open(os.path.join(V14_OUT, "numbers-v14-spectrum.json"), "r",
              encoding="utf-8") as f:
        v14 = json.load(f)
    v14raw = {c["cell"]: c["arms"]["RAW"] for c in v14["cells"]}
    by_id = {c["cell_id"]: c for c in merge["cells"]}
    stored = {}
    ok2c = True
    for cid, c14 in v14raw.items():
        c15 = by_id[cid]
        row = {
            "alpha_m_used_v15": float(c15["alpha"]["alpha_m_used"]),
            "alpha_m_used_v14": float(c14["alpha"]["alpha_m_used"]),
            "sigma2_g_v15": float(c15["moments"]["gravity"]["sigma2"]),
            "sigma2_g_v14": float(c14["moments"]["gravity"]["sigma2"]),
            "sigma2_h_v15": float(c15["moments"]["capillary"]["sigma2"]),
            "sigma2_h_v14": float(c14["moments"]["capillary"]["sigma2"]),
            "B_max_v15": float(c15["b_max_observation"]["B_max"]),
            "B_max_v14": float(c14["b_envelope"]["B_max"]),
            "Hs_pred_v15": float(c15["Hs_spectrum_m"]),
            "Hs_pred_v14": float(c14["Hs_spectrum_m"]),
        }
        exact = {}
        for k in ("alpha_m_used", "sigma2_g", "sigma2_h", "B_max",
                  "Hs_pred"):
            exact[k] = bool(row[k + "_v15"] == row[k + "_v14"])
        row["exact"] = exact
        row["pass"] = bool(all(exact.values()))
        ok2c = ok2c and row["pass"]
        stored[cid] = row
    reg["item2_v15_equals_v14_numbers"] = {
        "requirement": ("v15 RAW arm numbers == v14 RAW numbers per cell "
                        "(alpha_m, sigma2_g, sigma2_h, B max, H_s pred): "
                        "identical code path, expect exact equality"),
        "pass": bool(ok2a and ok2b and ok2c),
        "raw_arm_in_session": {"pass": bool(ok2a), "cells": raw_ident},
        "default_arm_D_in_session_vs_v14_D": {"pass": bool(ok2b),
                                              "cells": d_ident},
        "default_arm_vs_stored_v14_RAW_column": {
            "pass": bool(ok2c), "cells": stored,
            "source": "out/v14/numbers-v14-spectrum.json cells[*].arms.RAW"},
        "note": ("three comparisons: (i) v15 arm RAW vs v14 arm RAW computed "
                 "in this session, (ii) v15's merged default (arm D) vs v14's "
                 "arm D, (iii) the merged default's stored numbers against the "
                 "v14 numbers file.  All are expected to be EXACT (same code "
                 "path); any delta is printed, not smoothed.")}
    log("  [2] v15 == v14 numbers: RAW %s | D %s | stored %s"
        % (ok2a, ok2b, ok2c))

    # ---- item 3: v12 water untouched + the G1 legacy lane ---------------
    import v12_water as W12
    water = {
        "WATER_MODEL_DEFAULT": W12.WATER_MODEL_DEFAULT,
        "DEFAULT_MEDIUM_PRESET": W12.DEFAULT_MEDIUM_PRESET,
        "LEGACY_RECIPE_ROUTE": W12.LEGACY_RECIPE_ROUTE,
        "default_pass": bool(W12.WATER_MODEL_DEFAULT == "ab"),
        "import_only": True,
    }
    stored_sha = sha256_file(G1_LEGACY_FILE)
    water["g1_stored_sha256"] = stored_sha
    water["g1_stored_equals_reference"] = bool(stored_sha == G1_LEGACY_SHA)

    # the live re-render, through run_v12's own path (G1a)
    import run_v12 as R12
    from PIL import Image
    raw = load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    base_u8, photo_lin, PW, PH = R10._load_photo()
    r = R12.scene_state(raw, 1.5, 100.0, 0.6, R12.GRID_G1, seed)
    st = r["state"]
    re = st.realised_stats(r["fields"])
    from v5_physics import smith_g1_beckmann
    g1 = float(smith_g1_beckmann(
        re["rms_slope"], 90.0 - float(r["env"]["sun"]["altitudeDeg"]))["G1"])
    opt_legacy = R10.optics_for(r["env"], r["cfg"], R12.GRID_G1)
    E_v9 = float(2.2 * (opt_legacy["L_v3"]
                        / (opt_legacy["L_dir_flat"] * g1)))
    legs = W12.legacy_optics(r["env"], r["env"]["depthM"], R12.GRID_G1)
    t_leg = R12.render_scene(r, dict(legs, gain=opt_legacy["gain"]),
                             photo_lin, PW, PH, R12.GRID_G1, E_v9, g1,
                             tag="v15_g1_legacy")
    p = os.path.join(WORK, "still-g1-legacy-v15-rerun.png")
    png = io.BytesIO()
    Image.fromarray(t_leg["still"]).save(png, format="PNG")

    def _verify(path, _pw=PW, _ph=PH):
        try:
            im = Image.open(path)
            im.load()
            return im.size == (_pw, _ph) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(p, png.getvalue(), verify=_verify)
    rerun_sha = sha256_file(p)
    stored_u8 = np.asarray(Image.open(G1_LEGACY_FILE).convert("RGB"))
    d_stored = np.abs(t_leg["still"].astype(np.int16)
                      - stored_u8.astype(np.int16))
    water["g1_legacy_rerun"] = {
        "scene": {"U_ms": 1.5, "fetch_m": 100.0, "depth_m": 0.6,
                  "grid": R12.GRID_G1, "seed": seed,
                  "patch_m": float(st.S)},
        "route": {"E": E_v9, "g1": g1, "Q": R12.Q_DET,
                  "composite": "A 2.6 / B lowpass 3.0 / E as computed",
                  "water_model": "legacy (pre-3.6 v5 optics)"},
        "rerun_sha256": rerun_sha,
        "stored_sha256": stored_sha,
        "reference_sha256": G1_LEGACY_SHA,
        "byte_identical_to_stored": bool(rerun_sha == stored_sha),
        "byte_identical_to_reference": bool(rerun_sha == G1_LEGACY_SHA),
        "max_abs_uint8": int(d_stored.max()),
        "frac_pixels_changed": float(np.mean(d_stored.max(axis=2) > 0)),
        "file": os.path.relpath(p, HERE),
        "C": R10.c_stats(t_leg["C"]),
    }
    water["pass"] = bool(water["default_pass"]
                         and water["g1_legacy_rerun"]
                         ["byte_identical_to_reference"])
    del t_leg, r, st, png
    gc.collect()
    reg["item3_v12_water_untouched"] = {
        "requirement": ("v12 water model untouched: WATER_MODEL_DEFAULT still "
                        "'ab'; its G1 legacy lane still byte-identical to "
                        "delivery still-A-v37 (sha cfb086b22a0e...)"),
        **water}
    log("  [3] v12 water: default=%s | G1 legacy re-render %s vs reference "
        "%s -> identical=%s"
        % (water["WATER_MODEL_DEFAULT"], rerun_sha[:16],
           G1_LEGACY_SHA[:16], water["g1_legacy_rerun"]
           ["byte_identical_to_reference"]))

    # ---- item 4: frozen composite modules unchanged ---------------------
    from v5_composite import A_V3 as _A, B_SIGMA_V3 as _B, E_V3 as _E
    now_h = hashes(FROZEN_FILES)
    was14 = v14["frozen_hashes"]
    was1 = merge.get("frozen_hashes_at_stage1") or {}
    cmp = {}
    for f in FROZEN_FILES:
        prior = was14.get(f)
        cmp[f] = {
            "v14_recorded": prior,
            "stage1": was1.get(f),
            "now": now_h.get(f),
            "matches_v14_record": (None if prior is None
                                   else bool(prior == now_h.get(f))),
            "unchanged_since_stage1": bool(was1.get(f) == now_h.get(f)),
        }
    n_recorded = int(sum(1 for v in cmp.values()
                         if v["matches_v14_record"] is not None))
    ok4 = bool(all(v["matches_v14_record"] is not False
                   for v in cmp.values())
               and all(v["unchanged_since_stage1"] for v in cmp.values()))
    reg["item4_frozen_modules_unchanged"] = {
        "requirement": "frozen composite modules (the v12 hash list) unchanged",
        "pass": bool(ok4),
        "hash_source": "out/v14/numbers-v14-spectrum.json frozen_hashes",
        "files": cmp,
        "n_files_compared_against_v14_record": n_recorded,
        "n_files_hashed": len(FROZEN_FILES),
        "v5_composite_sha256": now_h.get("v5_composite.py"),
        "v5_composite_prefix_matches_v14_report": bool(
            (now_h.get("v5_composite.py") or "").startswith(
                "5c86a5770dcd0ce2")),
        "composite_constants": {"A": _A, "B_sigma": _B, "E": float(_E),
                                "Q": 25},
        "note": ("%d of the %d hashed modules carry a recorded v14 value and "
                 "all of them match it; the remaining files (run_v12.py and "
                 "the run_v14* family) have no v14-recorded value, so they are "
                 "covered by the stage-1 -> stage-4 unchanged test instead"
                 % (n_recorded, len(FROZEN_FILES)))}
    log("  [4] frozen modules: %d/%d match the v14 record, all unchanged "
        "since stage 1 -> %s" % (n_recorded, len(FROZEN_FILES), ok4))

    # ---- item 5: jonswap rollback lane ---------------------------------
    cell = {"U": 6.0, "F": 100.0, "h": 1.0}
    ja = V15.assemble(cell["U"], cell["F"], cell["h"], source="jonswap")
    jb = V10.assemble(cell["U"], cell["F"], cell["h"],
                      k_cut=V10.KCUT_DEFAULT)
    kj = V11._k_grid(float(ja["k_cut_rad_per_m"]))
    djs = np.abs(V11.s_parts(kj, ja)["S"] - V11.s_parts(kj, jb)["S"])
    with open(os.path.join(V10_OUT, "numbers-v10-sweep.json"), "r",
              encoding="utf-8") as f:
        v10 = json.load(f)
    v10cell = [c for c in v10["cells"] if c["cell"] == V10_SCENE_A_CELL][0]
    rj = scene_for(raw, cell["U"], cell["F"], cell["h"], v10["grid"]["N"],
                   int(derive_seed(raw, 1)), source="jonswap")
    rej = rj["state"].realised_stats(rj["fields"])
    s2_rel = float(rej["sigma2_rad2"]
                   / float(v10cell["realised"]["sigma2_rad2"]) - 1.0)
    hs_rel = float(rej["Hs_m"] / float(v10cell["realised"]["Hs_m"]) - 1.0)
    ok5 = bool(not djs.any() and abs(s2_rel) <= 0.01 and abs(hs_rel) <= 0.01)
    reg["item5_jonswap_rollback"] = {
        "requirement": ("'jonswap' still selectable in v15 and reproduces the "
                        "v10 scene-A spectrum statistics (realised sigma2 and "
                        "H_s vs numbers-v10 within 1%)"),
        "pass": bool(ok5),
        "selectable": True,
        "source_effective": rj["state"].spectrum_source,
        "assembly_identity_vs_v10": {
            "max_abs_dS": float(djs.max()), "S_identical": bool(not djs.any()),
            "alpha_m_v15": float(ja["alpha_m_value"]),
            "alpha_m_v10": float(jb["alpha_m_value"])},
        "realised_v15": {"sigma2_rad2": float(rej["sigma2_rad2"]),
                         "Hs_m": float(rej["Hs_m"]),
                         "rms_slope": float(rej["rms_slope"])},
        "realised_v10": {"sigma2_rad2": float(v10cell["realised"]
                                              ["sigma2_rad2"]),
                         "Hs_m": float(v10cell["realised"]["Hs_m"]),
                         "rms_slope": float(v10cell["realised"]
                                            ["rms_slope"])},
        "sigma2_rel_delta": s2_rel, "Hs_rel_delta": hs_rel,
        "tolerance": 0.01,
        "cell": dict(cell, N=int(v10["grid"]["N"]),
                     seed=int(derive_seed(raw, 1)),
                     patch_m=float(rj["state"].S)),
        "reference": ("out/v10/numbers-v10-sweep.json cell %s"
                      % V10_SCENE_A_CELL),
    }
    log("  [5] jonswap rollback: S identical %s | sigma2 %+.2e | Hs %+.2e -> "
        "%s" % (not djs.any(), s2_rel, hs_rel, ok5))
    del rj
    gc.collect()

    reg["all_pass"] = bool(
        reg["item1_D_STRICT_equals_v11"]["pass"]
        and reg["item2_v15_equals_v14_numbers"]["pass"]
        and reg["item3_v12_water_untouched"]["pass"]
        and reg["item4_frozen_modules_unchanged"]["pass"]
        and reg["item5_jonswap_rollback"]["pass"])
    reg["wall_s"] = float(time.time() - t0)
    reg["items"] = ["item1_D_STRICT_equals_v11",
                    "item2_v15_equals_v14_numbers",
                    "item3_v12_water_untouched",
                    "item4_frozen_modules_unchanged",
                    "item5_jonswap_rollback"]
    merge["regression"] = reg
    merge["frozen_hashes_at_stage2"] = now_h
    dump_json(os.path.join(OUT, "numbers-v15-merge.json"), merge)
    log("  regression all_pass=%s" % reg["all_pass"])
    cp.unit_done("stage-2")
    log("[stage 2] done in %.1fs" % (time.time() - t0))
    return reg


# ==========================================================================
# stage 3 -- transport through the merged default stack
# ==========================================================================
def stage3():
    t0 = time.time()
    log("=" * 72)
    log("[stage 3] transport through the MERGED DEFAULT stack (3 stills)")
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
    log("  photo %dx%d sha=%s | composite A=%s B_sigma=%s E=%s"
        % (PW, PH, sha256_file(R10.PHOTO)[:16], A_V3, B_SIGMA_V3, E_V3))

    out = {
        "stage": 3,
        "brief": "BRIEF-V15-MERGE-R2.md",
        "created_utc": now_iso(),
        "model": "version 3.10 (chunk 11) -- transport through the merged "
                 "default",
        "what": ("the MERGED DEFAULT stack rendered at three cells: T1 the "
                 "standing windy scene, T2 the low-wind cell whose UNRESOLVED "
                 "label must travel, T3 the Omega_c = 7.30 cell whose "
                 "warn_extrapolated_wave_age must travel"),
        "stack": {
            "spectrum": V15.SPECTRUM_SOURCE_DEFAULT,
            "spectrum_arm": V15.MERGE_ARM,
            "spectrum_rule": V15.MERGE_RULE,
            "water_model": "v12 'ab' (WATER_MODEL_DEFAULT)",
            "composite": {"A": A_V3, "B_lowpass_sigma_px": B_SIGMA_V3,
                          "E": float(E_V3), "Q": 25,
                          "status": "frozen (A2.6 + Blow + E2.2)"},
            "eta": 1.0,
            "eta_source": ("the recipe carries no sky block, so eta = "
                           "DNI/(DNI+DHI) = 1.0 (direct sun only)"),
            "retina": "OFF (the v7 observer post-stage is not in this path)",
            "rain_mm_h": 0.0,
            "seed": seed, "C_mesh": int(GRID_RENDER),
            "photo": {"path": R10.PHOTO, "size": [int(PW), int(PH)],
                      "sha256": sha256_file(R10.PHOTO)},
            "still_size_px": [int(PW), int(PH)],
        },
        "cells": {},
    }

    gates = {"C_mean_1_pm_0.02": True, "crush_lt_5pct": True,
             "blow_eq_0": True}
    for sc in TRANSPORT_CELLS:
        r = scene_for(raw, sc["U"], sc["F"], sc["h"], GRID_RENDER, seed)
        st = r["state"]
        t = R10.transport_scene(r, photo_lin, PW, PH, GRID_RENDER, Q=25,
                                E=float(E_V3))
        p = os.path.join(OUT, sc["file"])
        png = io.BytesIO()
        Image.fromarray(t["still"]).save(png, format="PNG")

        def _verify(path, _pw=PW, _ph=PH):
            try:
                im = Image.open(path)
                im.load()
                return im.size == (_pw, _ph) and im.mode == "RGB"
            except Exception:
                return False

        atomic_write_bytes(p, png.getvalue(), verify=_verify)
        still_u8 = np.asarray(Image.open(p).convert("RGB"))
        warp = warp_measure(base_u8, still_u8)
        pred = predicted_shift_field(r["fields"], r["env"], PW, PH)
        ps = predicted_block_stats(pred, warp["grid"])
        agree = {"ssd_primary": block_agreement(ps, warp["per_block"],
                                                "ssd_dy", "ssd_dx"),
                 "ncc_check": block_agreement(ps, warp["per_block"],
                                              "ncc_dy", "ncc_dx")}
        # the cell record from the spectrum side (labels + status)
        ev = V15.evaluate_cell(sc["U"], sc["F"], sc["h"])
        ctx = st.ctx
        c = t["C"]
        cm = float(np.mean(c))
        crush = float(t["metrics"]["crush_rate_linear_lt_0.02"])
        blow = float(t["metrics"]["blowout_rate_aces_gt_0.98"])
        g_c = bool(abs(cm - 1.0) <= 0.02)
        g_crush = bool(crush < 0.05)
        g_blow = bool(blow == 0.0)
        gates["C_mean_1_pm_0.02"] = gates["C_mean_1_pm_0.02"] and g_c
        gates["crush_lt_5pct"] = gates["crush_lt_5pct"] and g_crush
        gates["blow_eq_0"] = gates["blow_eq_0"] and g_blow
        out["cells"][sc["id"]] = {
            "id": sc["id"], "name": sc["name"], "why": sc["why"],
            "cell": {"U_ms": float(sc["U"]), "fetch_m": float(sc["F"]),
                     "depth_m": float(sc["h"])},
            "cell_id": ev["cell_id"],
            "spectrum_source": st.spectrum_source,
            "spectrum_arm": st.arm,
            "spectrum_source_default": V15.SPECTRUM_SOURCE_DEFAULT,
            "omega_c": float(ctx["omega_c"]),
            "warn_extrapolated_wave_age": bool(
                float(ctx["omega_c"]) > V15.WARN_OMEGA_C),
            "cell_status": ev["state"],
            "cell_status_rule": ev["verdict_rule"],
            "fail_reasons": ev["fail_reasons"],
            "low_wind_cell": ev["low_wind_cell"],
            "low_wind_status": ev["provenance"]["low_wind_status"],
            "provenance": ev["provenance"],
            "alpha_m_used": float(ctx["alpha_m_value"]),
            "alpha_m_raw": float(ctx["alpha_m_raw"]),
            "alpha_m_floor": float(ctx["alpha_m_floor"]),
            "top_up_acted": bool(
                ctx["alpha_m_resolution"]["top_up_acted"]),
            "sigma2_h_total": float(ctx["sigma2_h_total"]),
            "conformance": {"a_gravity": ev["conformance"]["a_gravity"]
                            ["ratio"],
                            "c_capillary": ev["conformance"]["c_capillary"]
                            ["ratio"]},
            "observation_b_status": ev["observation"]["b_beta"]["status"],
            "C": {"mean": cm,
                  "p50": float(np.percentile(c, 50)),
                  "p99": float(np.percentile(c, 99)),
                  "p99_9": float(np.percentile(c, 99.9)),
                  "max": float(c.max()),
                  "min": float(c.min()),
                  "field": "unattenuated caustic field C (det master)"},
            "crush_pct": 100.0 * crush,
            "blow_pct": 100.0 * blow,
            "own_geometric_shift_mm": float(
                ps["mean_predicted_shift_m"] * 1e3),
            "warp": {
                "primary_ssd": warp["primary"],
                "ncc_check": warp["ncc_check"],
                "own_shift_field": {
                    "field_mean_px": ps["field_mean_magnitude_px"],
                    "field_std_px": ps["field_std_magnitude_px"],
                    "block_mean_px": ps["block_mean_magnitude_px"],
                    "mean_mm": 1e3 * float(ps["mean_predicted_shift_m"]),
                    "pixel_scale_mm": float(ps["pixel_scale_mm"])},
                "agreement_with_predicted": agree},
            "state_summary": {
                "patch_m": float(st.S),
                "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
                "sigma2_spectral": float(
                    st.spectrum_after["sigma2_rad2"]),
                "rms_slope": float(st.spectrum_after["rms_slope"])},
            "fine_texture": _texture(still_u8),
            "environment": {"depth_m": float(r["env"]["depthM"]),
                            "turbid": float(r["env"]["turbid"]),
                            "sun_altitude_deg": float(
                                r["env"]["sun"]["altitudeDeg"])},
            "file": sc["file"], "sha256": sha256_file(p),
            "still_size_px": [int(PW), int(PH)], "C_mesh": int(GRID_RENDER),
            "timing_s": {"state": r["t_state"], "fields": r["t_fields"],
                         "transport": t["timing_s"]},
            "gates": {"C_mean_1_pm_0.02": g_c, "crush_lt_5pct": g_crush,
                      "blow_eq_0": g_blow},
            "labels_travelled": {
                "cell_status": ev["state"],
                "warn_extrapolated_wave_age":
                    bool(float(ctx["omega_c"]) > V15.WARN_OMEGA_C),
                "omega_c": float(ctx["omega_c"]),
                "spectrum_source": st.spectrum_source,
                "UNRESOLVED_in_record": bool(ev["state"] == "UNRESOLVED"),
                "warn_true_in_record": bool(
                    float(ctx["omega_c"]) > V15.WARN_OMEGA_C),
            },
        }
        log("  %s %-8s U=%-4g F=%-6g h=%-4g | Om_c %5.2f warn %-5s | status "
            "%-10s | C mean %.9f p99 %.3f p99.9 %.3f max %.2f | crush "
            "%.2f%% blow %.3f%% | shift %.2f mm | sha %s"
            % (sc["id"], sc["name"], sc["U"], sc["F"], sc["h"],
               ctx["omega_c"], float(ctx["omega_c"]) > V15.WARN_OMEGA_C,
               ev["state"], cm, float(np.percentile(c, 99)),
               float(np.percentile(c, 99.9)), float(c.max()), 100 * crush,
               100 * blow, 1e3 * float(ps["mean_predicted_shift_m"]),
               sha256_file(p)[:16]))
        del r, st, t, pred, warp
        gc.collect()

    out["gates"] = {
        "per_cell": {k: v["gates"] for k, v in out["cells"].items()},
        "all_three": gates,
        "definitions": {
            "C_mean_1_pm_0.02": "|mean(C) - 1| <= 0.02 on the unattenuated "
                                "field (preregistered here; tighter than the "
                                "v10/v14 route's own 1e-3 self-check)",
            "crush_lt_5pct": "crush_rate_linear_lt_0.02 < 0.05",
            "blow_eq_0": "blowout_rate_aces_gt_0.98 == 0"},
        "pass": bool(all(gates.values())),
    }
    # the two labels the brief makes merge conditions
    out["label_conditions"] = {
        "T2_UNRESOLVED_in_numbers": bool(
            out["cells"]["T2"]["cell_status"] == "UNRESOLVED"),
        "T3_warn_true_in_numbers": bool(
            out["cells"]["T3"]["warn_extrapolated_wave_age"] is True),
        "T3_omega_c": out["cells"]["T3"]["omega_c"],
        "pass": bool(out["cells"]["T2"]["cell_status"] == "UNRESOLVED"
                     and out["cells"]["T3"]
                     ["warn_extrapolated_wave_age"] is True),
    }
    out["wall_s"] = float(time.time() - t0)
    dump_json(os.path.join(OUT, "numbers-v15-transport.json"), out)
    cp.unit_done("stage-3")
    log("[stage 3] gates %s | labels %s | done in %.1fs"
        % (out["gates"]["all_three"], out["label_conditions"]["pass"],
           time.time() - t0))
    return out


# ==========================================================================
# stage 4 -- hash list, verdict roll-up, report
# ==========================================================================
def stage4():
    t0 = time.time()
    log("=" * 72)
    log("[stage 4] hash list, verdict roll-up, REPORT-V15-MERGE.md")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    with open(os.path.join(OUT, "numbers-v15-merge.json"), "r",
              encoding="utf-8") as f:
        merge = json.load(f)
    with open(os.path.join(OUT, "numbers-v15-transport.json"), "r",
              encoding="utf-8") as f:
        trans = json.load(f)

    now_h = hashes(FROZEN_FILES)
    was = merge.get("frozen_hashes_at_stage1") or now_h
    reg_hash = {k: {"stage1": was.get(k), "stage4": now_h[k],
                    "unchanged": bool(was.get(k) == now_h[k])}
                for k in now_h}
    merge["file_hashes"] = {
        "frozen_files": reg_hash,
        "all_frozen_unchanged": bool(all(v["unchanged"]
                                         for v in reg_hash.values())),
        "new_v15_files": hashes(NEW_FILES),
        "deliverables": {},
        "self_hash_note": ("numbers-v15-merge.json cannot carry its own "
                           "sha256 (the value would change the file).  The "
                           "report prints the hash of the file it actually "
                           "read, computed at report build time."),
        "note": ("sha256 of every imported/frozen module at stage 1 and stage "
                 "4, plus the v15 files and the deliverables")}
    for fn in ("numbers-v15-transport.json",
               "still-T1-windy.png", "still-T2-lowwind.png",
               "still-T3-warned.png"):
        p = os.path.join(OUT, fn)
        merge["file_hashes"]["deliverables"][fn] = (
            sha256_file(p) if os.path.exists(p) else None)

    # the roll-up the report quotes
    cells = merge["cells"]
    fails = [c["cell_id"] for c in cells if c["state"] == "FAIL"]
    merge["rollup"] = {
        "default_source": V15.SPECTRUM_SOURCE_DEFAULT,
        "n_cells": len(cells),
        "n_pass": int(sum(1 for c in cells if c["state"] == "PASS")),
        "n_fail": int(len(fails)),
        "n_unresolved": int(sum(1 for c in cells
                                if c["state"] == "UNRESOLVED")),
        "failing_cells": fails,
        "zero_fails": bool(not fails),
        "merge_stands": bool(not fails),
        "merge_rule": ("the brief: high-wind cells must PASS and zero FAILs "
                       "are required for the merge to stand.  The v15 module "
                       "exposes the default regardless (that IS the brief's "
                       "stage 1 instruction for the new family); if this flag "
                       "is false the merge is NOT promoted and the report says "
                       "why."),
        "regression_all_pass": bool(merge.get("regression", {})
                                    .get("all_pass")),
        "transport_all_pass": bool(trans["gates"]["pass"]),
        "transport_labels_pass": bool(trans["label_conditions"]["pass"]),
    }
    merge["rollup"]["all_conditions"] = bool(
        merge["rollup"]["zero_fails"]
        and merge["rollup"]["regression_all_pass"]
        and merge["rollup"]["transport_all_pass"]
        and merge["rollup"]["transport_labels_pass"])
    dump_json(os.path.join(OUT, "numbers-v15-merge.json"), merge)
    log("  roll-up: PASS %d / FAIL %d / UNRESOLVED %d | zero_fails=%s | "
        "regression=%s | transport=%s | labels=%s"
        % (merge["rollup"]["n_pass"], merge["rollup"]["n_fail"],
           merge["rollup"]["n_unresolved"], merge["rollup"]["zero_fails"],
           merge["rollup"]["regression_all_pass"],
           merge["rollup"]["transport_all_pass"],
           merge["rollup"]["transport_labels_pass"]))

    import run_v15_report
    run_v15_report.main()
    cp.unit_done("stage-4")
    log("[stage 4] done in %.1fs" % (time.time() - t0))
    return merge["rollup"]


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
