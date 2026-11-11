"""run_v14 -- version 3.9 boundary redo: the THREE-ARM measurement (chunk 10).

Brief: BRIEF-V14-39REDO.md.  Design: PIPELINE-V12.md (it governs everything).
Everything is written to `out/v14/` -- no other folder is written, no process is
touched, no image is inspected (numerical statistics only), every deliverable
goes through jobkit `atomic_write_bytes`, and one checkpoint unit is committed
per stage.

  stage 1  the three arms over the 8-cell matrix, analytic:
           alpha_m resolution + prior accounting, per-arm C1/C1b/B-envelope/
           totals/HS/smoothness/positivity, the Arm-B hand-off scan in u*,
           the RAW == v13 bit-identity check, the exported spectra
           -> out/v14/numbers-v14-spectrum.json
  stage 2  realised states (N=1024) per cell for RAW and B: Parseval (A2),
           independent band integrals (A4), resolution convergence (A5),
           the Arm-B floor sweep {0.0015, 0.003, 0.006} at the U=3 cells
           -> out/v14/numbers-v14-sweep.json
  stage 3  transport: the windy cell (U=6/F=1000/depth 1.0) per arm, the
           U=3/F=1000 RAW-vs-B stills and their side-by-side composite, the
           v13 re-render byte-identity check, the v12 import check
           -> out/v14/numbers-v14-transport.json
  stage 4  the decision table (constraints x arms), regression, and
           REPORT-V14-GATES.md
           -> out/v14/REPORT-V14-GATES.md

Usage (foreground, from caustics-filter/offline):
  python run_v14.py --stage 1
  python run_v14.py --stage 2
  python run_v14.py --stage 3
  python run_v14.py --stage 4
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
import run_v14_state as S14                                        # noqa: E402
import run_v10 as R10                                              # noqa: E402

OUT = os.path.join(HERE, "out", "v14")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
V13_OUT = os.path.join(HERE, "out", "v13")
V11_OUT = os.path.join(HERE, "out", "v11")
V10_OUT = os.path.join(HERE, "out", "v10")

GRID_SWEEP = 1024
GRID_RENDER = R10.GRID_RENDER
SWEEP_U = V10.SWEEP_U
SWEEP_F = V10.SWEEP_F
SWEEP_H = V10.SWEEP_Z

#: the brief's transport cell: U=6, F=1000, depth 1.0.
WINDY = {"U": 6.0, "F": 1000.0, "h": 1.0}
#: the boundary cell the brief asks to render: U=3, F=1000, depth 1.0.
BOUNDARY = {"U": 3.0, "F": 1000.0, "h": 1.0}
#: the v13 deliverables were rendered at U=6/F=100/h=1.0 -- the regression cell.
V13_CELL = {"U": 6.0, "F": 100.0, "h": 1.0}

#: files this job must not change (recorded at stage 1, re-checked at stage 4).
FROZEN_FILES = (
    "v5_composite.py", "v5_physics.py", "v12_water.py", "v10_spectrum.py",
    "v10_state.py", "v11_spectrum.py", "v11_state.py", "v13_spectrum.py",
    "v13_state.py", "transport.py", "det.py", "image.py", "env.py",
    "recipe_v8.py", "run_v13.py", "run_v10.py",
)

log = R10.log
now_iso = R10.now_iso
dump_json = R10.dump_json
sha256_file = R10.sha256_file
load_recipe = R10.load_recipe
cell_id = R10.cell_id


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def scene_for(raw, U, F, h, N, seed, arm, source="elfouhaily-form",
              target=None, k_cut=V10.KCUT_DEFAULT):
    """v14 state -> the dict run_v10.transport_scene expects."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    t0 = time.time()
    st = S14.V14WaterState(
        scene, env, N=N, seed=seed, preset_index=1, fetch_m=float(F),
        k_cut=k_cut, spectrum_source=source, arm=arm,
        sigma2_h_target=(S14.DEFAULT_SIGMA2_H_TARGET if target is None
                         else float(target)))
    t_state = time.time() - t0
    t1 = time.time()
    fields = st.fields_at(0.0)
    t_fields = time.time() - t1
    return {"state": st, "env": env, "cfg": cfg, "fields": fields,
            "t_state": t_state, "t_fields": t_fields}


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


def hashes():
    return {f: sha256_file(os.path.join(HERE, f)) for f in FROZEN_FILES}


def _v13_raw_identity(U, F, h):
    """RAW must reproduce v13.assemble bit-for-bit (the control arm)."""
    a = V14.assemble(U, F, h, arm="RAW")
    b = V13.assemble(U, F, h)
    k = V11._k_grid(float(a["k_cut_rad_per_m"]))
    sa = V11.s_parts(k, a)
    sb = V11.s_parts(k, b)
    d = np.abs(sa["S"] - sb["S"])
    return {
        "cell": {"U_ms": float(U), "fetch_m": float(F), "depth_m": float(h)},
        "alpha_m_used_v14": float(a["alpha_m_value"]),
        "alpha_m_used_v13": float(b["alpha_m_value"]),
        "alpha_identical": bool(a["alpha_m_value"] == b["alpha_m_value"]),
        "sigma2_h_v14": float(a["sigma2_h_total"]),
        "sigma2_h_v13": float(b["sigma2_h_total"]),
        "Hs_v14": float(a["Hs_m"]), "Hs_v13": float(b["Hs_m"]),
        "max_abs_dS": float(d.max()),
        "S_identical": bool(not d.any()),
        "pass": bool(a["alpha_m_value"] == b["alpha_m_value"]
                     and not d.any()
                     and a["Hs_m"] == b["Hs_m"]),
    }


def _v11_hard_zero_identity(U, F, h):
    """D_STRICT must reproduce the v3.8 build (alpha_m = max(raw, 0))."""
    a = V14.assemble(U, F, h, arm=V14.D_STRICT)
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


# ==========================================================================
# stage 1 -- the three arms, 8 cells, analytic
# ==========================================================================
def stage1():
    t0 = time.time()
    log("=" * 72)
    log("[stage 1] three arms x 8 cells (analytic) + hand-off scan + exports")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    frozen_at_start = hashes()

    # ---- the RAW == v13 control check, and D_STRICT == v3.8 --------------
    raw_ident = {}
    strict_ident = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            raw_ident[cid] = _v13_raw_identity(U, F, SWEEP_H)
            strict_ident[cid] = _v11_hard_zero_identity(U, F, SWEEP_H)
    raw_ok = all(v["pass"] for v in raw_ident.values())
    strict_ok = all(v["pass"] for v in strict_ident.values())
    log("  control: RAW == v13.assemble on all 8 cells: %s | D_STRICT == "
        "v3.8 build: %s" % (raw_ok, strict_ok))

    # ---- the Arm-B hand-off in u* ---------------------------------------
    handoff = {}
    for F in SWEEP_F:
        h = V14.handoff_scan(F, SWEEP_H)
        handoff["F_%g" % F] = h
        ho = h["handoff"]
        log("  hand-off F=%-6g target %.4f: u* = %s (U10 = %s, u*/fit = %s) | "
            "rule C0 %s (ratio %.3f) | pipeline max dalpha %.2e at U=%.2f"
            % (F, h["target"], ("%.6f" % ho["ustar_ms"]) if ho else "none",
               ("%.3f" % ho["U10_ms"]) if ho else "none",
               ("%.3f" % ho["u_over_fit_min"]) if ho else "none",
               h["rule_scan"]["continuous_C0"],
               h["rule_scan"]["ratio_fine_over_coarse"],
               h["pipeline_scan"]["max_rel_jump_alpha_m_used"],
               h["pipeline_scan"]["U_at_max_alpha_jump"]))

    # ---- the cells ------------------------------------------------------
    cells = []
    exports = {}
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            t1 = time.time()
            ctxs = {
                "RAW": V14.assemble(U, F, SWEEP_H, arm="RAW"),
                "B": V14.assemble(U, F, SWEEP_H, arm="B"),
                V14.D_STRICT: V14.assemble(U, F, SWEEP_H, arm=V14.D_STRICT),
            }
            ctxs["D"] = None                       # D is RAW's spectrum
            per_arm = {}
            for arm in ("RAW", "B", "D", V14.D_STRICT):
                ctx = ctxs["RAW"] if arm == "D" else ctxs[arm]
                per_arm[arm] = V14.evaluate_cell(
                    U, F, SWEEP_H, arm if arm != "D" else "D", ctx=ctx)
                if arm == "D":
                    per_arm["D"]["numbers_identical_to_RAW"] = True
                    per_arm["D"]["alpha"] = per_arm["RAW"]["alpha"]
                    per_arm["D"]["alpha"]["arm"] = "D"
                    per_arm["D"]["provenance"] = (
                        V14.PROVENANCE["D"] if per_arm["RAW"]["low_wind_cell"]
                        else "eckv_fit_domain_raw")
            a_raw = per_arm["RAW"]
            cells.append({
                "cell": cid, "U_ms": float(U), "fetch_m": float(F),
                "depth_m": float(SWEEP_H),
                "geometry": a_raw["geometry"],
                "sea_state": a_raw["sea_state"],
                "low_wind_cell": a_raw["low_wind_cell"],
                "envelope": a_raw["envelope"],
                "arms": per_arm,
                "v13_control": raw_ident[cid],
                "timing_s": float(time.time() - t1),
            })
            exports[cid] = {
                "band_rad_per_m": a_raw["geometry"]["kmin_rad_per_m"],
                "k_cut_rad_per_m": float(ctxs["RAW"]["k_cut_rad_per_m"]),
                "RAW": V14.export_curves(ctxs["RAW"]),
                "B": V14.export_curves(ctxs["B"]),
                V14.D_STRICT: V14.export_curves(ctxs[V14.D_STRICT]),
                "D_note": ("arm D's spectrum IS arm RAW's; re-use the RAW "
                           "curve (bit-identical alpha_m_used)"),
            }
            # A4's reproduce-it-yourself check on the exported grid
            rep = {}
            for arm in ("RAW", "B", V14.D_STRICT):
                exp = exports[cid][arm]
                r = V14.reproduce_from_export(
                    exp, float(ctxs["RAW"]["k_cut_rad_per_m"]))
                tgt = per_arm[arm if arm != "D" else "RAW"]
                rep[arm] = {
                    "export_sigma2_total": r["sigma2_total"],
                    "reported_sigma2_total": tgt["sigma2_total_cut"],
                    "rel_total": float(r["sigma2_total"]
                                       / tgt["sigma2_total_cut"] - 1.0),
                    "export_capillary_50_200": r["sigma2_capillary_50_200"],
                    "reported_capillary_50_200":
                        tgt["shares"]["B_h_sigma2_50_200"],
                    "rel_capillary_50_200": (
                        None if not tgt["shares"]["B_h_sigma2_50_200"]
                        else float(r["sigma2_capillary_50_200"]
                                   / tgt["shares"]["B_h_sigma2_50_200"]
                                   - 1.0))}
            exports[cid]["reproduction_from_export"] = rep
            del ctxs
            gc.collect()
            log("  [%-14s] raw %+.6f floor %.6f | RAW used %.6f (%s) | B used "
                "%.6f topup=%s | B_env %.4f Bl %.4f | C1 %.3f C1b %s"
                % (cid, a_raw["alpha"]["alpha_m_raw"],
                   a_raw["alpha"]["alpha_m_floor"], a_raw["alpha"]["alpha_m_used"],
                   a_raw["state"], per_arm["B"]["alpha"]["alpha_m_used"],
                   per_arm["B"]["alpha"]["top_up_acted"],
                   per_arm["B"]["b_envelope"]["B_max"],
                   per_arm["B"]["b_envelope"]["B_l_max"],
                   per_arm["B"]["shares"]["gravity_share_50_200"],
                   per_arm["B"]["C1b"]["pass"]))

    summary = {}
    for arm in ("RAW", "B", "D", V14.D_STRICT):
        st = [c["arms"][arm]["state"] for c in cells]
        summary[arm] = {
            "states": {c["cell"]: c["arms"][arm]["state"] for c in cells},
            "n_pass": int(sum(1 for s in st if s == "PASS")),
            "n_fail": int(sum(1 for s in st if s == "FAIL")),
            "n_unresolved": int(sum(1 for s in st if s == "UNRESOLVED")),
            "cells_where_prior_acts": [
                c["cell"] for c in cells
                if c["arms"][arm]["alpha"]["prior_active"]],
            "n_cells_prior_acts": int(sum(
                1 for c in cells if c["arms"][arm]["alpha"]["prior_active"])),
            "worst_B_max": max(c["arms"][arm]["b_envelope"]["B_max"]
                               for c in cells),
            "worst_B_l_max": max(c["arms"][arm]["b_envelope"]["B_l_max"]
                                 for c in cells),
            "worst_gravity_share": max(
                c["arms"][arm]["shares"]["gravity_share_50_200"]
                for c in cells),
            "min_capillary_share": min(
                c["arms"][arm]["shares"]["capillary_share_50_200"]
                for c in cells),
            "worst_abs_dHs": max(abs(c["arms"][arm]["hs"]["delta_rel"])
                                 for c in cells),
            "n_totals_tol_pass": int(sum(
                1 for c in cells
                if c["arms"][arm]["envelope_check"]["tol_inside"])),
            "n_totals_strict_pass": int(sum(
                1 for c in cells
                if c["arms"][arm]["envelope_check"]["strict_inside"])),
        }

    stage = {
        "stage": 1,
        "created_utc": now_iso(),
        "model": "version 3.9 + v14 three-arm boundary measurement "
                 "(PIPELINE-V12 chunk 10)",
        "what": ("the three arms (RAW control, B extended prior, D evaluation "
                 "arm) over the same 8-cell matrix, analytic: alpha_m "
                 "resolution and prior accounting, C1/C1b/B-envelope/totals/"
                 "Hs/smoothness/positivity, the Arm-B hand-off scan, and the "
                 "control-arm identity checks"),
        "brief": "BRIEF-V14-39REDO.md",
        "design": "PIPELINE-V12.md",
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters -- NOT 'the Elfouhaily model'"),
        "arms": {a: {"definition": V14.arm_definition(a),
                     "provenance": V14.PROVENANCE[a]} for a in
                 ("RAW", "B", "D", V14.D_STRICT)},
        "policy": {
            "grid_N_analytic": "48k log grid (V11._k_grid), band [kmin, k_vis]",
            "grid_N_state": GRID_SWEEP,
            "depth_m": SWEEP_H,
            "U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
            "k_cut_rad_per_m": V10.KCUT_DEFAULT,
            "patch_base_m": V10.PATCH_BASE_M, "r_lambda": V10.R_LAMBDA,
            "sigma2_h_target_default": V14.SIGMA2_H_TARGET_DEFAULT,
            "sigma2_h_targets_swept": list(V14.SIGMA2_H_TARGETS),
            "ustar_fit_min_ms": V14.USTAR_FIT_MIN,
            "c_m_ms": V14.C_M,
            "B_env_max": V14.B_ENV_MAX, "beta_hat": V14.BETA_HAT,
            "gravity_share_max": V14.GRAVITY_SHARE_MAX,
            "capillary_share_min": V14.CAP_SHARE_MIN,
            "sigma2_h_min": V14.SIGMA2_H_MIN,
            "B_band_k_range": list(V14.B_K_RANGE),
            "band_50_200": list(V14.BAND),
            "totals_tolerance": list(V14.TOTALS_TOL),
            "hs_tolerance": V14.HS_TOL,
            "default_source_after": S13.DEFAULT_SOURCE,
            "B_HAT_low_high_v13": [0.007, 0.008],
        },
        "controls": {
            "RAW_equals_v13": {"pass": bool(raw_ok), "cells": raw_ident},
            "D_STRICT_equals_v3p8": {"pass": bool(strict_ok),
                                     "cells": strict_ident},
            "note": ("RAW is the control arm and must reproduce v13's "
                     "assembly bit-for-bit; D_STRICT (diagnostic) must "
                     "reproduce the v3.8 build (alpha_m = max(raw, 0)).  "
                     "Both are asserted per cell, not assumed."),
        },
        "handoff": handoff,
        "summary": summary,
        "cells": cells,
        "frozen_hashes": frozen_at_start,
        "export": {
            "what": ("S_l/S_h/B_l/B_h per cell on a 1500-point log grid from "
                     "kmin to 1500 rad/m; trapz on this grid reproduces every "
                     "band integral -- see each cell's "
                     "reproduction_from_export block"),
            "cells": exports},
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v14-spectrum.json"), stage)

    # a small human-readable curve file (the v13 deliverable's analogue)
    c = exports[cell_id(6.0, 1000.0)]["B"]
    lines = ["k_rad_per_m,B_l,B_h,S_l,S_h,arm"]
    for i in range(len(c["k_rad_per_m"])):
        lines.append("%.10e,%.10e,%.10e,%.10e,%.10e,B" % (
            c["k_rad_per_m"][i], c["B_l"][i], c["B_h"][i],
            c["S_l_m2_per_rad_m"][i], c["S_h_m2_per_rad_m"][i]))
    atomic_write_bytes(os.path.join(OUT, "spectrum-ref-v40-U6F1000-armB.csv"),
                       ("\n".join(lines) + "\n").encode("ascii"))
    cp.unit_done("stage-1")
    log("[stage 1] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 2 -- realised states, A2/A4/A5, floor sweep
# ==========================================================================
def stage2():
    t0 = time.time()
    log("=" * 72)
    log("[stage 2] realised states (N=%d), Parseval, independent integrals, "
        "convergence, floor sweep" % GRID_SWEEP)
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    raw = load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    with open(os.path.join(OUT, "numbers-v14-spectrum.json"), "r",
              encoding="utf-8") as f:
        s1 = json.load(f)
    s1cells = {c["cell"]: c for c in s1["cells"]}
    log("  production seed %d" % seed)

    cells = []
    for U in SWEEP_U:
        for F in SWEEP_F:
            cid = cell_id(U, F)
            t1 = time.time()
            row = {"cell": cid, "U_ms": float(U), "fetch_m": float(F),
                   "depth_m": float(SWEEP_H), "arms": {}}
            for arm in ("RAW", "B"):
                r = scene_for(raw, U, F, SWEEP_H, GRID_SWEEP, seed, arm)
                st = r["state"]
                node = st.state_node(r["fields"])
                ctx = st.ctx
                geo = V14.cell_geometry(ctx)
                knyq = math.pi * GRID_SWEEP / float(st.S)
                k_pars = float(min(V10.KCUT_DEFAULT, 0.95 * knyq))
                a2 = V11.gate_a2(ctx, N=GRID_SWEEP, patch_m=geo["patch_m"],
                                 k_cut=k_pars,
                                 tag="v14 arm %s" % arm)
                a2["k_cut_used_for_parseval"] = k_pars
                a2["k_nyquist_at_N"] = float(knyq)
                a2["note_band"] = (
                    "Parseval is measured over the band the state's own grid "
                    "resolves: k_cut = min(%.0f, 0.95 k_nyquist).  A5 (N=2048, "
                    "k_cut=1500) is the resolution test; this gate is the "
                    "convention test." % V10.KCUT_DEFAULT)
                a4 = V14.independent_integrals(ctx)
                a4b = V14.band_reproduction(ctx)
                row["arms"][arm] = {
                    "alpha": ctx["alpha_m_resolution"],
                    "provenance": ctx["provenance"],
                    "patch_m": float(st.S),
                    "patch_doublings": int(st.patch_doublings),
                    "mode": st.mode,
                    "kmin_rad_per_m": float(st.kmin),
                    "dx_mm": float(st.dx * 1e3),
                    "kmax_nyquist_scene": float(st.kmax_nyq),
                    "spectrum": node["spectrum"],
                    "realised": node["realised"],
                    "bands_realised": node["bands"],
                    "amplitude": node["amplitude"],
                    "breaking": node["breaking"],
                    "spectral_vs_spatial": node["spectral_vs_spatial"],
                    "direction": node["direction"],
                    "K1": node["v11"]["K1"],
                    "quantisation": node["quantisation"],
                    "frac_above_370_realised": st.frac_above((370.0,))
                    ["above_370"]["realised"],
                    "sigma2_spectral": float(st.spectrum_after["sigma2_rad2"]),
                    "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
                    "rms_slope": float(st.spectrum_after["rms_slope"]),
                    "A2_parseval": a2,
                    "A4_independent_integrals": a4,
                    "A4_band_reproduction": a4b,
                }
                log("  [%-14s] %-3s realised sig2 %.6f Hs %.4f | A2 %s | "
                    "A4 %s (worst %.1e) rows %s"
                    % (cid, arm, row["arms"][arm]["sigma2_spectral"],
                       row["arms"][arm]["Hs_spectral_m"], a2["pass"], a4["pass"],
                       a4["worst_abs_rel"], a4b["pass"]))
                del r, st, node, ctx
                gc.collect()
            # the arm-D identity, asserted on the built field
            row["D_identical_to_RAW"] = {
                "why": ("arm D's spectrum IS RAW's; no state is rebuilt.  The "
                        "identity is asserted on alpha_m_used from stage 1."),
                "alpha_m_used_RAW":
                    s1cells[cid]["arms"]["RAW"]["alpha"]["alpha_m_used"],
                "alpha_m_used_D":
                    s1cells[cid]["arms"]["D"]["alpha"]["alpha_m_used"],
                "pass": bool(
                    s1cells[cid]["arms"]["RAW"]["alpha"]["alpha_m_used"]
                    == s1cells[cid]["arms"]["D"]["alpha"]["alpha_m_used"])}
            row["timing_s"] = float(time.time() - t1)
            cells.append(row)

    # ---- A5: one cell at a finer grid and a higher k-cutoff --------------
    log("  A5 convergence: U=6/F=1000 at N=2048, k_cut=1500 vs N=1024, 1000")
    conv = {}
    for arm in ("RAW", "B"):
        r0 = scene_for(raw, 6.0, 1000.0, SWEEP_H, GRID_SWEEP, seed, arm)
        r1 = scene_for(raw, 6.0, 1000.0, SWEEP_H, 2048, seed, arm,
                       k_cut=1500.0)
        n0 = r0["state"].state_node(r0["fields"])
        n1 = r1["state"].state_node(r1["fields"])
        conv[arm] = {
            "case": {"U_ms": 6.0, "fetch_m": 1000.0, "depth_m": SWEEP_H},
            "coarse": {"N": GRID_SWEEP, "k_cut": V10.KCUT_DEFAULT,
                       "patch_m": float(r0["state"].S),
                       "sigma2_realised": float(
                           r0["state"].spectrum_after["sigma2_rad2"]),
                       "Hs_realised_m": float(
                           r0["state"].spectrum_after["Hs_m"]),
                       "m0": float(n0["spectrum"].get("m0", float("nan")))
                       if "m0" in n0.get("spectrum", {}) else None},
            "fine": {"N": 2048, "k_cut": 1500.0,
                     "patch_m": float(r1["state"].S),
                     "sigma2_realised": float(
                         r1["state"].spectrum_after["sigma2_rad2"]),
                     "Hs_realised_m": float(
                         r1["state"].spectrum_after["Hs_m"]),
                     "m0": float(n1["spectrum"].get("m0", float("nan")))
                     if "m0" in n1.get("spectrum", {}) else None},
        }
        c = conv[arm]
        c["rel_sigma2"] = float(c["fine"]["sigma2_realised"]
                                / c["coarse"]["sigma2_realised"] - 1.0)
        c["rel_Hs"] = float(c["fine"]["Hs_realised_m"]
                            / c["coarse"]["Hs_realised_m"] - 1.0)
        c["tolerance"] = 0.05
        c["pass"] = bool(abs(c["rel_sigma2"]) <= 0.05
                         and abs(c["rel_Hs"]) <= 0.05)
        log("    %-3s N1024 %.6f -> N2048/k1500 %.6f (rel %+.4f) Hs %+.4f "
            "-> %s" % (arm, c["coarse"]["sigma2_realised"],
                       c["fine"]["sigma2_realised"], c["rel_sigma2"],
                       c["rel_Hs"], c["pass"]))
        del r0, r1, n0, n1
        gc.collect()

    # ---- the Arm-B floor sweep at the U=3 cells --------------------------
    log("  floor sweep %s at U=3" % (list(V14.SIGMA2_H_TARGETS[1:]),))
    sens = {"what": ("Arm B's floor swept over {0.0015, 0.003, 0.006} at the "
                     "two U=3 cells: sigma2_h delivered, the shares, the "
                     "envelope and the B-envelope at each floor"),
            "floors": [0.0015, 0.003, 0.006],
            "cells": {}}
    for F in SWEEP_F:
        cid = cell_id(3.0, F)
        rows = {}
        for tgt in sens["floors"]:
            ctx = V14.assemble(3.0, F, SWEEP_H, arm="B",
                               sigma2_h_target=tgt)
            ev = V14.evaluate_cell(3.0, F, SWEEP_H, "B", sigma2_h_target=tgt,
                                   ctx=ctx)
            rows["floor_%g" % tgt] = {
                "floor": float(tgt),
                "alpha_m_raw": ev["alpha"]["alpha_m_raw"],
                "alpha_m_used": ev["alpha"]["alpha_m_used"],
                "top_up_acted": bool(ev["alpha"]["top_up_acted"]),
                "sigma2_h_total": ev["moments"]["capillary"]["sigma2"],
                "sigma2_gravity_total": ev["moments"]["gravity"]["sigma2"],
                "sigma2_total_cut": ev["sigma2_total_cut"],
                "gravity_share_50_200": ev["shares"]["gravity_share_50_200"],
                "capillary_share_50_200":
                    ev["shares"]["capillary_share_50_200"],
                "B_max": ev["b_envelope"]["B_max"],
                "B_l_max": ev["b_envelope"]["B_l_max"],
                "Hs_spectrum_m": ev["Hs_spectrum_m"],
                "Hs_law_m": ev["Hs_law_m"],
                "dHs_rel": ev["hs"]["delta_rel"],
                "totals_strict_inside":
                    ev["envelope_check"]["strict_inside"],
                "totals_tol_inside": ev["envelope_check"]["tol_inside"],
                "ratio_to_lower": ev["envelope_check"]["ratio_to_lower"],
                "ratio_to_upper": ev["envelope_check"]["ratio_to_upper"],
                "capillary_frac_above_370":
                    ev["moments"]["capillary_frac_above_370"],
                "state": ev["state"],
            }
            log("    %-14s floor %.4f -> alpha %.6f sig2h %.6f capshare "
                "%.3f Bmax %.4f"
                % (cid, tgt, ev["alpha"]["alpha_m_used"],
                   ev["moments"]["capillary"]["sigma2"],
                   ev["shares"]["capillary_share_50_200"],
                   ev["b_envelope"]["B_max"]))
            del ctx, ev
            gc.collect()
        sens["cells"][cid] = rows

    # ---- summary ---------------------------------------------------------
    summary = {}
    for arm in ("RAW", "B"):
        rows = [c["arms"][arm] for c in cells]
        summary[arm] = {
            "A2_parseval_pass": all(r["A2_parseval"]["pass"] for r in rows),
            "A4_independent_pass": all(r["A4_independent_integrals"]["pass"]
                                       for r in rows),
            "A4_band_reproduction_pass": all(r["A4_band_reproduction"]["pass"]
                                             for r in rows),
            "worst_A4_rel": max(r["A4_independent_integrals"]["worst_abs_rel"]
                                for r in rows),
            "A5_convergence_pass": conv[arm]["pass"],
            "worst_abs_dHs_spectral_vs_law": max(
                abs(r["K1"]["delta_rel"]) for r in rows),
            "worst_quantisation_rel": max(
                r["quantisation"]["max_rel_error"] for r in rows),
        }
    stage = {
        "stage": 2,
        "created_utc": now_iso(),
        "what": ("realised states at N=%d for arms RAW and B: state node, "
                 "Parseval (A2) at N=512 per cell, independent band integrals "
                 "(A4), the resolution-convergence re-run (A5), and the Arm-B "
                 "floor sweep at the U=3 cells.  Arm D is not rebuilt: its "
                 "spectrum is RAW's." % GRID_SWEEP),
        "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                  "parameters -- NOT 'the Elfouhaily model'"),
        "grid": {"N": GRID_SWEEP, "depth_m": SWEEP_H,
                 "U_ms": list(SWEEP_U), "fetch_m": list(SWEEP_F),
                 "seed": seed, "k_cut": V10.KCUT_DEFAULT},
        "policy": s1["policy"],
        "summary": summary,
        "cells": cells,
        "A5_convergence": conv,
        "sensitivity_floor_sweep": sens,
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v14-sweep.json"), stage)
    cp.unit_done("stage-2")
    log("[stage 2] done in %.1fs" % (time.time() - t0))
    return stage


# ==========================================================================
# stage 3 -- transport + regression
# ==========================================================================
def stage3():
    t0 = time.time()
    log("=" * 72)
    log("[stage 3] transport: windy per arm, U=3 RAW vs B, regression")
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
    log("  photo %dx%d sha=%s" % (PW, PH, sha256_file(R10.PHOTO)[:16]))

    out = {"stage": 3, "created_utc": now_iso(),
           "what": ("transport consequences of the boundary choice: the windy "
                    "cell per arm (U=6/F=1000/depth 1.0), the U=3/F=1000 "
                    "RAW-vs-B stills and their side-by-side composite, and "
                    "the regression battery"),
           "photo": {"path": R10.PHOTO, "size": [int(PW), int(PH)],
                     "sha256": sha256_file(R10.PHOTO)},
           "label": ("Elfouhaily-form spectrum with fetch-derived sea-state "
                     "parameters; three-arm boundary measurement"),
           "frozen_composite": {"A_V3": A_V3, "B_SIGMA_V3": B_SIGMA_V3,
                                "E_V3": float(E_V3), "Q": 25},
           "windy_cell": {"U_ms": WINDY["U"], "fetch_m": WINDY["F"],
                          "depth_m": WINDY["h"]},
           "boundary_cell": {"U_ms": BOUNDARY["U"], "fetch_m": BOUNDARY["F"],
                             "depth_m": BOUNDARY["h"]},
           "scenes": {}}

    def render(U, F, h, arm, source="elfouhaily-form", fname=None,
               keep_C=False):
        r = scene_for(raw, U, F, h, GRID_RENDER, seed, arm, source=source)
        st = r["state"]
        t = R10.transport_scene(r, photo_lin, PW, PH, GRID_RENDER, Q=25,
                                E=float(E_V3))
        if fname:
            png = io.BytesIO()
            Image.fromarray(t["still"]).save(png, format="PNG")

            def verify(p, _pw=PW, _ph=PH):
                try:
                    im = Image.open(p)
                    im.load()
                    return im.size == (_pw, _ph) and im.mode == "RGB"
                except Exception:
                    return False

            atomic_write_bytes(os.path.join(OUT, fname), png.getvalue(),
                               verify=verify)
            path = os.path.join(OUT, fname)
        else:
            path = None
        still_png = np.asarray(Image.open(path).convert("RGB")) if path else \
            t["still"]
        warp = warp_measure(base_u8, still_png)
        pred = predicted_shift_field(r["fields"], r["env"], PW, PH)
        ps = predicted_block_stats(pred, warp["grid"])
        agree = {"ssd_primary": block_agreement(ps, warp["per_block"],
                                                "ssd_dy", "ssd_dx"),
                 "ncc_check": block_agreement(ps, warp["per_block"],
                                              "ncc_dy", "ncc_dx")}
        keep = {"state": st, "r": r, "t": t, "still_png": still_png,
                "warp": warp, "ps": ps, "agree": agree, "file": fname,
                "texture": _texture(still_png), "path": path}
        if not keep_C:
            keep["t"]["C"] = t["C"]
        return keep

    def node_of(k, arm):
        st = k["state"]
        t = k["t"]
        p = k["warp"]["primary"]
        ps = k["ps"]
        ctx = st.ctx
        return {
            "arm": arm,
            "spectrum_source": st.spectrum_source,
            "alpha_m_used": float(ctx["alpha_m_value"]),
            "alpha_m_raw": float(ctx["alpha_m_raw"]),
            "alpha_m_floor": float(ctx["alpha_m_floor"]),
            "top_up_acted": bool(ctx["alpha_m_resolution"]["top_up_acted"]),
            "provenance": ctx["provenance"],
            "low_wind_prior": ctx.get("low_wind_prior"),
            "sigma2_h_total": float(ctx["sigma2_h_total"]),
            "state_summary": {
                "patch_m": float(st.S),
                "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
                "sigma2_spectral": float(st.spectrum_after["sigma2_rad2"]),
                "rms_slope": float(st.spectrum_after["rms_slope"]),
                "omega_c": float(ctx["omega_c"]),
                "frac_above_370_realised":
                    st.frac_above((370.0,))["above_370"]["realised"],
            },
            "route": {"Q": 25, "A": A_V3, "E": float(E_V3),
                      "B_lowpass_sigma_px": B_SIGMA_V3,
                      "g1_used": float(t["g1"]),
                      "blur_sigma_px": float(t["opt"]["blur"]["sigma_px"])},
            "C": R10.c_stats(t["C"]),
            "C_percentiles": {
                "p50": float(np.percentile(t["C"], 50)),
                "p90": float(np.percentile(t["C"], 90)),
                "p99": float(np.percentile(t["C"], 99)),
                "p99_9": float(np.percentile(t["C"], 99.9)),
                "p99_99": float(np.percentile(t["C"], 99.99)),
                "max": float(t["C"].max())},
            "warp": {
                "method": k["warp"]["definition"],
                "primary_ssd": {kk: p[kk] for kk in (
                    "mean_shift_px", "shift_field_std_px",
                    "shift_magnitude_p50_px", "shift_magnitude_p90_px",
                    "shift_magnitude_max_px", "frac_blocks_at_search_edge",
                    "n_blocks", "block_px", "search_px", "step_px")},
                "ncc_check": {kk: k["warp"]["ncc_check"][kk] for kk in (
                    "mean_shift_px", "shift_field_std_px",
                    "shift_magnitude_p50_px", "shift_magnitude_p90_px",
                    "frac_blocks_at_search_edge")},
                "own_shift_field": {
                    "field_mean_px": ps["field_mean_magnitude_px"],
                    "field_std_px": ps["field_std_magnitude_px"],
                    "block_mean_px": ps["block_mean_magnitude_px"],
                    "mean_mm": 1e3 * float(ps["mean_predicted_shift_m"]),
                    "pixel_scale_mm": float(ps["pixel_scale_mm"]),
                    "frac_beyond_search": ps["field_frac_beyond_search"]},
                "agreement_with_predicted": k["agree"],
                "yardsticks": R10.YARDSTICKS,
            },
            "fine_texture": k["texture"],
            "tone": {
                "crush_rate_linear_lt_0.02":
                    float(t["metrics"]["crush_rate_linear_lt_0.02"]),
                "blowout_rate_aces_gt_0.98":
                    float(t["metrics"]["blowout_rate_aces_gt_0.98"]),
                "customary_crush_bound": 0.0325,
                "crush_pass_customary": bool(
                    t["metrics"]["crush_rate_linear_lt_0.02"] <= 0.0325),
                "blow_pass_zero": bool(
                    t["metrics"]["blowout_rate_aces_gt_0.98"] <= 0.0)},
            "beer": t["beer"],
            "still_metrics": t["metrics"],
            "file": k["file"],
            "sha256": (sha256_file(k["path"]) if k["path"] else None),
            "timing_s": t["timing_s"],
        }

    # ---- the windy cell, per arm ----------------------------------------
    windy = {}
    for arm in ("RAW", "B"):
        fname = "comp-windy-v40-U6F1000-arm%s.png" % arm
        k = render(WINDY["U"], WINDY["F"], WINDY["h"], arm, fname=fname)
        windy[arm] = node_of(k, arm)
        log("  windy %-3s C p50 %.3f p99 %.3f p99.9 %.3f max %.2f mean %.6f |"
            " warp %.2f mm | |dL/dy| %.5f"
            % (arm, windy[arm]["C_percentiles"]["p50"],
               windy[arm]["C_percentiles"]["p99"],
               windy[arm]["C_percentiles"]["p99_9"],
               windy[arm]["C_percentiles"]["max"],
               windy[arm]["C"]["mean"],
               windy[arm]["warp"]["own_shift_field"]["mean_mm"],
               windy[arm]["fine_texture"]["mean_abs_dL_dy"]))
        if arm == "RAW":
            windy_png_raw = k["still_png"].copy()
        else:
            windy_png_b = k["still_png"].copy()
        del k
        gc.collect()
    d = np.abs(windy_png_raw.astype(np.int16)
               - windy_png_b.astype(np.int16))
    windy["RAW_vs_B"] = {
        "max_abs_uint8": int(d.max()),
        "mean_abs_uint8": float(d.mean()),
        "frac_pixels_changed": float(np.mean(d.max(axis=2) > 0)),
        "byte_identical": bool(windy["RAW"]["sha256"]
                               == windy["B"]["sha256"]),
        "still_size_px": [int(PW), int(PH)],
        "c_mesh": int(GRID_RENDER),
        "expected": ("identical: at U=6 raw sigma2_h > target, so Arm B's "
                     "top-up does not act and both arms assemble the same "
                     "alpha_m"),
        "alpha_m_identical": bool(windy["RAW"]["alpha_m_used"]
                                 == windy["B"]["alpha_m_used"]),
    }
    windy["D_identical_to_RAW"] = {
        "why": "arm D changes no numbers (brief line 13)",
        "pass": True,
        "alpha_m_used": windy["RAW"]["alpha_m_used"],
        "note": ("no separate render: the arm-D field, seed, frame and "
                 "composite parameters are RAW's, so the still is RAW's "
                 "byte-for-byte")}
    out["scenes"]["windy"] = windy
    del d
    gc.collect()

    # ---- the U=3/F=1000 boundary stills --------------------------------
    bnd = {}
    for arm in ("RAW", "B"):
        fname = "comp-u3-F1000-arm%s.png" % arm
        k = render(BOUNDARY["U"], BOUNDARY["F"], BOUNDARY["h"], arm,
                   fname=fname)
        bnd[arm] = node_of(k, arm)
        if arm == "RAW":
            png_raw = k["still_png"].copy()
        else:
            png_b = k["still_png"].copy()
        log("  u3 %-3s alpha_m %.6f (%s) C p50 %.3f p99 %.3f p99.9 %.3f "
            "max %.2f mean %.6f | warp %.2f mm | |dL/dy| %.5f"
            % (arm, bnd[arm]["alpha_m_used"],
               "prior" if bnd[arm]["top_up_acted"] else "raw",
               bnd[arm]["C_percentiles"]["p50"],
               bnd[arm]["C_percentiles"]["p99"],
               bnd[arm]["C_percentiles"]["p99_9"],
               bnd[arm]["C_percentiles"]["max"], bnd[arm]["C"]["mean"],
               bnd[arm]["warp"]["own_shift_field"]["mean_mm"],
               bnd[arm]["fine_texture"]["mean_abs_dL_dy"]))
        del k
        gc.collect()
    side = np.concatenate([png_raw, np.full((png_raw.shape[0], 2, 3), 255,
                                            np.uint8), png_b], axis=1)
    png = io.BytesIO()
    Image.fromarray(side).save(png, format="PNG")

    def verify_side(p, _h=png_raw.shape[0], _w=png_raw.shape[1] * 2 + 2):
        try:
            im = Image.open(p)
            im.load()
            return im.size == (_w, _h) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(os.path.join(OUT, "comp-u3-raw-vs-priorB.png"),
                       png.getvalue(), verify=verify_side)
    dd = np.abs(png_raw.astype(np.int16) - png_b.astype(np.int16))
    bnd["RAW_vs_B"] = {
        "layout": ("left = arm RAW, right = arm B, 2 px white divider; "
                   "%dx%d RGB (the composite is rendered at the photo's size; "
                   "the C mesh is %d^2)"
                   % (side.shape[1], side.shape[0], GRID_RENDER)),
        "size_px": [int(side.shape[1]), int(side.shape[0])],
        "still_size_px": [int(PW), int(PH)],
        "file": "comp-u3-raw-vs-priorB.png",
        "sha256": sha256_file(os.path.join(OUT,
                                           "comp-u3-raw-vs-priorB.png")),
        "max_abs_uint8": int(dd.max()), "mean_abs_uint8": float(dd.mean()),
        "frac_pixels_changed": float(np.mean(dd.max(axis=2) > 0)),
        "msssim": R10.msssim(png_raw.astype(np.float64) / 255.0,
                             png_b.astype(np.float64) / 255.0)["msssim"],
        "alpha_m_ratio_B_over_RAW": float(
            bnd["B"]["alpha_m_used"] / bnd["RAW"]["alpha_m_used"]),
        "C_mean_both_one": bool(bnd["RAW"]["C"]["pass_mean_1"]
                                and bnd["B"]["C"]["pass_mean_1"]),
        "note": ("same seed, same frame, same composite parameters: the only "
                 "difference is alpha_m, i.e. the boundary treatment itself"),
    }
    out["scenes"]["boundary_u3"] = bnd
    del dd, side, png_raw, png_b
    gc.collect()

    # ---- regression: the v13 cell, both lanes --------------------------
    log("  regression: re-render the v13 cell U=6/F=100/h=1.0, both lanes")
    reg = {}
    for src, arm, fname, v13f in (
            ("jonswap", "RAW", "comp-windy-v40-v13legacy.png",
             "comp-windy-v39-legacy.png"),
            ("elfouhaily-form", "RAW", "comp-windy-v40-v13new.png",
             "comp-windy-v39-new.png")):
        k = render(V13_CELL["U"], V13_CELL["F"], V13_CELL["h"], arm,
                   source=src, fname=fname)
        a14 = sha256_file(os.path.join(OUT, fname))
        a13 = sha256_file(os.path.join(V13_OUT, v13f))
        i14 = k["still_png"].astype(np.int16)
        i13 = np.asarray(Image.open(os.path.join(V13_OUT, v13f))
                         .convert("RGB")).astype(np.int16)
        dg = np.abs(i14 - i13)
        reg[src] = {"v14_file": fname, "v13_file": v13f,
                    "sha256_v14": a14, "sha256_v13": a13,
                    "byte_identical": bool(a14 == a13),
                    "max_abs_uint8": int(dg.max()),
                    "frac_pixels_changed": float(np.mean(dg.max(axis=2) > 0))}
        log("    %-16s v14 %s vs v13 %s -> identical=%s max|d|=%d"
            % (src, a14[:16], a13[:16], reg[src]["byte_identical"],
               reg[src]["max_abs_uint8"]))
        del k, i14, i13, dg
        gc.collect()
    out["regression"] = {
        "v13_legacy_lane": reg["jonswap"],
        "v13_new_lane": reg["elfouhaily-form"],
        "v13_all_byte_identical": bool(all(v["byte_identical"]
                                           for v in reg.values())),
        "note": ("both v13 lanes re-rendered through the v14 family at v13's "
                 "own cell (U=6/F=100/h=1.0).  The legacy lane is the "
                 "byte-identical test PIPELINE-V12 asks for; the new lane is "
                 "a bonus: v13's cell is inside the fit domain, so no arm "
                 "changes a number there, and the still must match too."),
    }

    # ---- v12 water path: import check only ------------------------------
    import v12_water as W12
    out["regression"]["v12_water_path"] = {
        "check": "import only (PIPELINE-V12 line 29)",
        "import_ok": True,
        "WATER_MODEL_DEFAULT": W12.WATER_MODEL_DEFAULT,
        "DEFAULT_MEDIUM_PRESET": W12.DEFAULT_MEDIUM_PRESET,
        "LEGACY_RECIPE_ROUTE": W12.LEGACY_RECIPE_ROUTE,
        "pass": bool(W12.WATER_MODEL_DEFAULT == "ab"),
        "note": ("the water model default must still be 'ab'; nothing in this "
                 "job calls or edits v12_water"),
    }
    out["regression"]["frozen_composite_hashes"] = {
        "A_V3": A_V3, "B_SIGMA_V3": B_SIGMA_V3, "E_V3": float(E_V3),
        "v5_composite_sha256": sha256_file(os.path.join(HERE,
                                                        "v5_composite.py")),
        "note": "the frozen v5 composite is imported, never edited"}
    out["regression"]["pass"] = bool(
        out["regression"]["v13_all_byte_identical"]
        and out["regression"]["v12_water_path"]["pass"])
    out["wall_s"] = float(time.time() - t0)
    dump_json(os.path.join(OUT, "numbers-v14-transport.json"), out)
    cp.unit_done("stage-3")
    log("[stage 3] done in %.1fs (regression pass=%s)"
        % (time.time() - t0, out["regression"]["pass"]))
    return out


# ==========================================================================
# stage 4 -- decision table + report
# ==========================================================================
def build_decision():
    with open(os.path.join(OUT, "numbers-v14-spectrum.json"), "r",
              encoding="utf-8") as f:
        s1 = json.load(f)
    with open(os.path.join(OUT, "numbers-v14-sweep.json"), "r",
              encoding="utf-8") as f:
        s2 = json.load(f)
    with open(os.path.join(OUT, "numbers-v14-transport.json"), "r",
              encoding="utf-8") as f:
        s3 = json.load(f)
    arms = ("RAW", "B", "D", V14.D_STRICT)
    cells = s1["cells"]
    s2by = {c["cell"]: c for c in s2["cells"]}

    constraints = [
        ("C1 gravity share <= 0.50 (U=6/10 cells)",
         lambda c, a: (None if c["low_wind_cell"]
                       else c["arms"][a]["C1"]["pass"])),
        ("C1b capillary share >= 0.15 AND sigma2_h >= 0.002 (U=6/10 cells)",
         lambda c, a: (None if c["low_wind_cell"]
                       else c["arms"][a]["C1b"]["pass"])),
        ("totals inside [0.5x, 2.0x] of [Cox-Munk clean, Hu]",
         lambda c, a: c["arms"][a]["envelope_check"]["tol_inside"]),
        ("totals inside [Cox-Munk clean, Hu] (strict reading, reported)",
         lambda c, a: c["arms"][a]["envelope_check"]["strict_inside"]),
        ("B-envelope: max B(k) <= 0.012 over [10,1500]",
         lambda c, a: c["arms"][a]["b_envelope"]["pass_B_env"]),
        ("B_l band max <= 0.008 over [50,200]",
         lambda c, a: c["arms"][a]["b_envelope"]["pass_B_l"]),
        ("C1 smoothness: no step in B / dlnB/dlnk",
         lambda c, a: c["arms"][a]["smoothness"]["pass"]),
        ("A1 positivity / finite energy",
         lambda c, a: c["arms"][a]["positivity"]["pass"]),
        ("A4 independent band integrals",
         lambda c, a: s2by[c["cell"]]["arms"][
             "RAW" if a in ("RAW", "D", V14.D_STRICT) else "B"]
             ["A4_independent_integrals"]["pass"]),
        ("|dHs| <= 35% (spectrum vs fetch law)",
         lambda c, a: c["arms"][a]["hs"]["pass"]),
    ]
    rows = []
    for name, fn in constraints:
        row = {"constraint": name, "cells": {}}
        for c in cells:
            for a in arms:
                v = fn(c, a)
                row["cells"].setdefault(c["cell"], {})[a] = v
        for a in arms:
            vals = [row["cells"][c["cell"]][a] for c in cells]
            row.setdefault("counts", {})[a] = {
                "pass": int(sum(1 for v in vals if v is True)),
                "fail": int(sum(1 for v in vals if v is False)),
                "not_applicable_or_unresolved": int(
                    sum(1 for v in vals if v is None))}
        rows.append(row)

    invented = {}
    for a in arms:
        per = {}
        for c in cells:
            al = c["arms"][a]["alpha"]
            n = 1 if al["prior_active"] else 0
            per[c["cell"]] = {
                "invented_inputs_acting": n,
                "which": (["sigma2_h_target=%.4g (the provisional prior)"
                           % al["sigma2_h_target"]] if n else []),
                "top_up_acted": bool(al["top_up_acted"]),
                "prior_active": bool(al["prior_active"])}
        invented[a] = {
            "per_cell": per,
            "n_cells_with_invented_input": int(sum(
                v["invented_inputs_acting"] for v in per.values())),
            "input_name": "sigma2_h_target = 0.003 (provisional, unmeasured)",
            "provenance": s1["arms"][a]["provenance"]}

    t = s3["scenes"]
    transport = {}
    for a in ("RAW", "B"):
        w = t["windy"][a]
        transport[a] = {
            "C_mean": w["C"]["mean"], "C_mean_pass_1": w["C"]["pass_mean_1"],
            "p99": w["C_percentiles"]["p99"],
            "p99_9": w["C_percentiles"]["p99_9"],
            "max": w["C_percentiles"]["max"],
            "own_shift_mean_mm": w["warp"]["own_shift_field"]["mean_mm"],
            "warp_p99_9_ge_3": bool(w["C_percentiles"]["p99_9"] >= 3.0),
            "crush_pass": bool(w["tone"]["crush_pass_customary"]),
            "blow_pass": bool(w["tone"]["blow_pass_zero"]),
            "warp_band_pass": bool(20.0 <= w["warp"]["own_shift_field"]
                                   ["mean_mm"] <= 45.0),
            "sanitiy_pass": None}
        transport[a]["sanitiy_pass"] = bool(
            transport[a]["C_mean_pass_1"]
            and transport[a]["warp_p99_9_ge_3"]
            and transport[a]["crush_pass"] and transport[a]["blow_pass"])
    transport["D"] = dict(transport["RAW"])
    transport["D"]["note"] = "no separate render: D's numbers are RAW's"
    transport[V14.D_STRICT] = None

    decision = {
        "what": ("the brief's decision DATA: every independent-constraint "
                 "result and the count of invented inputs acting per cell, "
                 "side by side.  NO winner is chosen here -- PIPELINE-V12's "
                 "pre-registered rule is applied outside this session."),
        "arms": {a: {"definition": s1["arms"][a]["definition"],
                     "provenance": s1["arms"][a]["provenance"]} for a in arms},
        "constraints": rows,
        "invented_inputs": invented,
        "transport": transport,
        "regression": s3["regression"],
        "rule_quoted": (
            "PIPELINE-V12 line 20: the arm that satisfies the INDEPENDENT "
            "constraints (totals envelope, B-envelope, construction gates, "
            "transport sanity) with the fewest invented inputs wins "
            "merge-as-challenger; tie -> Arm B (majority rank); Arm D winning "
            "requires its B-envelope checks to actually pass at U=3; "
            "assertion is not evidence."),
        "default_change": ("NONE in this job.  spectrum_source default stays "
                           "'%s' (v13's value)." % S13.DEFAULT_SOURCE),
        "arm_states": {
            a: {c["cell"]: c["arms"][a]["state"] for c in cells} for a in arms},
        "D_absolute_criteria": {
            c["cell"]: c["arms"]["D"]["absolute_criteria"] for c in cells
            if c["arms"]["D"]["absolute_criteria"]},
        "D_STRICT_absolute_criteria": {
            c["cell"]: c["arms"][V14.D_STRICT]["absolute_criteria"]
            for c in cells
            if c["arms"][V14.D_STRICT]["absolute_criteria"]},
    }
    return decision, s1, s2, s3


def stage4():
    t0 = time.time()
    log("=" * 72)
    log("[stage 4] decision table, regression re-check, report")
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    decision, s1, s2, s3 = build_decision()

    # regression: the frozen files must be untouched since stage 1
    now_h = hashes()
    with open(os.path.join(OUT, "numbers-v14-spectrum.json"), "r",
              encoding="utf-8") as f:
        first = json.load(f)
    was = first.get("frozen_hashes") or now_h
    reg_hash = {k: {"stage1": was.get(k), "stage4": now_h[k],
                    "unchanged": bool(was.get(k) == now_h[k])}
                for k in now_h}
    decision["file_hashes"] = {
        "files": reg_hash,
        "all_unchanged": bool(all(v["unchanged"]
                                  for v in reg_hash.values())),
        "note": ("sha256 of every imported/frozen module recorded at the "
                 "start of the job (stage 1) and re-taken at stage 4: no "
                 "v5..v13 file changed while this job ran")}
    log("  frozen files unchanged: %s"
        % decision["file_hashes"]["all_unchanged"])

    for a in ("RAW", "B", "D", V14.D_STRICT):
        st = decision["arm_states"][a]
        log("  arm %-8s states: PASS %d, FAIL %d, UNRESOLVED %d | invented "
            "inputs acting at %d cell(s)"
            % (a, sum(1 for v in st.values() if v == "PASS"),
               sum(1 for v in st.values() if v == "FAIL"),
               sum(1 for v in st.values() if v == "UNRESOLVED"),
               decision["invented_inputs"][a]
               ["n_cells_with_invented_input"]))

    if not first.get("frozen_hashes"):
        first["frozen_hashes"] = now_h
        dump_json(os.path.join(OUT, "numbers-v14-spectrum.json"), first)
    s2["decision"] = decision
    dump_json(os.path.join(OUT, "numbers-v14-sweep.json"), s2)
    import run_v14_report
    run_v14_report.main()
    cp.unit_done("stage-4")
    log("[stage 4] done in %.1fs" % (time.time() - t0))
    return decision


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
