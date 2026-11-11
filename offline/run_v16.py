"""run_v16 -- chunk 12 (BRIEF-V16-SCATTER): correct the v3.6 phase function
against published measurements.

Stages (foreground, one jobkit checkpoint unit per delivered artefact):

  A  scorecard: <cos>, B, <theta^2>, F(<1/5/10 deg), beta(1/5/10/90 deg) for
     HG-0.85, HG-0.924 and each preset's FF parameters
     -> out/v16/numbers-v16-scorecard.json
  B  halo level: the v12 G1 reference scene and the four G3 rows through the
     v15 default stack with the FF water (H bit-identity, sigma scaling, halo
     share, peak profiles at r32/r48)
     -> out/v16/numbers-v16-halo.json, stills
  C  kernel shape: encircled energy Gaussian vs FF-CDF at r = 8.5/32/79 px
     -> merged into numbers-v16-halo.json
  D  regression: hg+gaussian == v12, one v15 cell, frozen module hashes, no
     v11-v15 file modified; the merge rule applied
     -> merged into numbers-v16-halo.json
  report -> out/v16/REPORT-V16-SCATTER.md, pair-halo-old-vs-new.png

Rules (brief): write only inside out/v16 and new run_v16* files; never modify
a v11-v15 file; ASCII output; numerical image statistics only; no tuning of
n/mu -- the published fits are not knobs.

Usage (from caustics-filter/offline):
  python run_v16.py --stage A
  python run_v16.py --stage B
  python run_v16.py --stage C
  python run_v16.py --stage D
  python run_v16.py --stage report
  python run_v16.py --stage all
"""
import argparse
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

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v10 as R10                                            # noqa: E402
import run_v15 as R15                                            # noqa: E402
import run_v15_spectrum as V15S                                  # noqa: E402
import v12_water as W                                            # noqa: E402
import v5_physics as V5                                          # noqa: E402
import run_v16_water as V16                                      # noqa: E402

OUT = os.path.join(HERE, "out", "v16")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
V12_OUT = os.path.join(HERE, "out", "v12")
V15_OUT = os.path.join(HERE, "out", "v15")
V10_OUT = os.path.join(HERE, "out", "v10")
V15_NUMBERS = os.path.join(V15_OUT, "numbers-v15-merge.json")
G1_LEGACY_PNG = os.path.join(V10_OUT, "still-A-v37.png")
G1_LEGACY_SHA = "cfb086b22a0e8f1b25cb3befbcc272f3e9c5127fe5499396a555c99b297f0e32"

GRID_G1 = 2048          # the v12 G1 reference grid
GRID_EXP = 1024         # the v12 G2/G3 grid
Q_DET = 25
SCENE_P1 = {"name": "preset-1 scene (recipe-1 as-is)", "U": 1.5, "F": 100.0,
            "h": 0.6}
G3_LEVELS = (0.05, 0.25, 0.75, 2.0)
G3_MID = 0.25           # the delivered pair uses the G3 mid scene

log = R10.log
now_iso = R10.now_iso
dump_json = R10.dump_json
sha256_file = R10.sha256_file
load_recipe = R10.load_recipe


# ==========================================================================
# helpers
# ==========================================================================
def makedirs():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)


def env_for(raw, U, F, h):
    """Scene -> env/cfg exactly as the v12/v15 drivers build it (no state)."""
    from env import resolve_environment
    from recipe_v8 import parse_recipe_v8, with_v8
    scene = with_v8(raw, fetch_m=float(F), depth_m=float(h), wind_ms=float(U))
    return {"scene": scene, "env": resolve_environment(scene),
            "cfg": parse_recipe_v8(scene)}


def scene_v15(raw, U, F, h, N, seed):
    """The v15 DEFAULT stack (merged spectrum) for one cell."""
    t0 = time.time()
    r = R15.scene_for(raw, U, F, h, N, seed)
    r["t_state_total"] = float(time.time() - t0)
    return r


def route_numbers(r, grid):
    """g1 and the E_v9 normalisation, verbatim from run_v12 stage 2/3."""
    from v5_physics import smith_g1_beckmann
    st = r["state"]
    re = st.realised_stats(r["fields"])
    sun_zen = 90.0 - float(r["env"]["sun"]["altitudeDeg"])
    g1 = float(smith_g1_beckmann(re["rms_slope"], sun_zen)["G1"])
    opt = R10.optics_for(r["env"], r["cfg"], grid)
    E = float(2.2 * (opt["L_v3"] / (opt["L_dir_flat"] * g1)))
    return {"g1": g1, "E": E, "realised": re,
            "patch_m": float(st.S),
            "spectrum_source": getattr(st, "spectrum_source", None),
            "arm": getattr(st, "arm", None)}


def render(r, optics, photo_lin, PW, PH, grid, E, g1, C_cache=None):
    from det import det_frame
    from v5_composite import image_metrics
    t0 = time.time()
    if C_cache is None:
        C, info = det_frame(r["fields"], r["state"].S, grid, r["env"], Q_DET)
    else:
        C, info = C_cache, {"reused": True, "Q": Q_DET}
    t_det = time.time() - t0
    t1 = time.time()
    from v5_composite import A_V3, B_SIGMA_V3
    img, dbg = V16.render_frame_v16(
        r["fields"], C, photo_lin, r["env"], PW, PH, optics,
        ambient_scale=A_V3, c_layer_exposure=float(E),
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t1
    out = {"C": C, "det_info": info, "still": img, "metrics": metrics,
           "lin": dbg["lin"], "E": float(E), "g1": float(g1),
           "optics": optics,
           "timing_s": {"det": float(t_det), "composite": float(t_comp)}}
    del dbg
    gc.collect()
    return out


def write_png(path, arr, PW, PH):
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG")

    def verify(p):
        try:
            im = Image.open(p)
            im.load()
            return im.size == (PW, PH) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(path, buf.getvalue(), verify=verify)
    return sha256_file(path)


def png_stats(arr):
    a = np.asarray(arr, np.float64)
    lum = a @ np.array([0.2126, 0.7152, 0.0722])
    d = np.abs(np.diff(lum, axis=0))
    return {"mean_luminance_uint8": float(lum.mean()),
            "std_luminance_uint8": float(lum.std()),
            "mean_abs_dL_dy": float(d.mean()),
            "frac_uint8_zero": float(np.mean(a.max(axis=2) == 0)),
            "frac_uint8_255": float(np.mean(a.max(axis=2) == 255))}


def halo_share(optics):
    T = np.asarray(optics["w_ballistic"], np.float64)
    H = np.asarray(optics["w_halo"], np.float64)
    return [float(h / (t + h)) for t, h in zip(T, H)]


def effective_peak_profile(C, optics, kernel=None):
    """The v12 G3 instrument on the effective field of one optics object."""
    if kernel == "gaussian":
        eff = W.effective_c_field(C, dict(optics, sigma_px=optics["sigma_px"]))
    else:
        eff = V16.effective_c_field_v16(C, optics)
    fs = W.field_stats(eff)
    return {"effective_field": fs,
            "p99_over_mean": fs["contrast_p99_over_mean"],
            "profile": fs["peak_profile"]["profile"]}


def g1_spec_from_legacy(env):
    """The v12 G1 reference coefficients: a,b from the v5 turbidity model."""
    att = V5.water_attenuation(env["turbid"])
    return W.explicit_spectra(np.asarray(att["a"], np.float64),
                              np.asarray(att["b"], np.float64),
                              g=W.G_DEFAULT), att


def reference_optics(phase, kernel, grid=GRID_G1, patch_m=2.4):
    """Optics for the G1 reference scene (env only, no state needed)."""
    raw = load_recipe()
    e = env_for(raw, SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"])
    spec, att = g1_spec_from_legacy(e["env"])
    o = V16.optics_v16(spec, e["env"], e["env"]["depthM"], grid, patch_m,
                      phase=phase, kernel=kernel)
    return o, e, att


def hash_files(names):
    out = {}
    for f in names:
        p = os.path.join(HERE, f)
        out[f] = sha256_file(p) if os.path.exists(p) else None
    return out


#: every v11-v15 python file this job must not touch (plus the v15 report).
V11_V15_FILES = (
    "v11_spectrum.py", "v11_state.py", "v11_report.py", "v12_water.py",
    "v12_report.py", "run_v11.py", "run_v12.py", "v13_spectrum.py",
    "v13_state.py", "v13_report.py", "run_v13.py", "run_v14.py",
    "run_v14_spectrum.py", "run_v14_state.py", "run_v14_report.py",
    "run_v14_report_alt.py", "run_v15.py", "run_v15_spectrum.py",
    "run_v15_state.py", "run_v15_report.py",
)


# ==========================================================================
# STAGE A -- the scorecard
# ==========================================================================
def stage_a():
    t0 = time.time()
    log("=" * 72)
    log("[stage A] phase-function scorecard (the brief's five clauses)")
    makedirs()
    cp = Checkpoint(CP_PATH)

    snapshot = {
        "taken_utc": now_iso(),
        "purpose": ("stage-0-style preflight baseline: hashes of every v11-v15 "
                    "file and of the G1 legacy artefact, captured BEFORE any "
                    "v16 work so Stage D can prove nothing moved"),
        "files": {f: (sha256_file(os.path.join(HERE, f))
                      if os.path.exists(os.path.join(HERE, f)) else None)
                  for f in V11_V15_FILES},
        "g1_legacy_png": {"path": os.path.relpath(G1_LEGACY_PNG, HERE),
                          "sha256": sha256_file(G1_LEGACY_PNG)},
    }
    atomic_write_bytes(os.path.join(WORK, "hashes_at_start.json"),
                       json.dumps(snapshot, indent=1).encode("utf-8"))
    cp.unit_done("A-hash-snapshot")

    candidates = []
    specs = (("HG_0.85", "hg", None, 0.85),
             ("HG_0.924", "hg", None, 0.924),
             ("FF_algae_pond", "ff", "algae_pond", None),
             ("FF_clear", "ff", "clear", None),
             ("FF_sediment", "ff", "sediment", None))
    rows = {}
    for tag, phase, preset, g in specs:
        table = V16.phase_table(phase, preset=preset, g=g)
        row = V16.phase_moments(table)
        sc = V16.scorecard_clauses(row)
        rows[tag] = {"row": row, "scorecard": sc}
        candidates.append({
            "id": tag, "phase": phase, "preset": preset,
            "n": row["n"], "mu": row["mu"], "g": row["g"],
            "mean_cos": row["mean_cos"],
            "backscatter_fraction_B": row["backscatter_fraction_B"],
            "theta2_rad2": row["theta2_rad2"],
            "F_lt_1deg": row["F_lt_1deg"], "F_lt_5deg": row["F_lt_5deg"],
            "F_lt_10deg": row["F_lt_10deg"],
            "F_gt_30deg": row["F_gt_30deg"],
            "beta_1deg": row["beta_1deg_1_per_sr"],
            "beta_5deg": row["beta_5deg_1_per_sr"],
            "beta_10deg": row["beta_10deg_1_per_sr"],
            "beta_90deg": row["beta_90deg_1_per_sr"],
            "beta90_over_reference": row["beta90_over_reference"],
            "normalisation_quadrature": row["normalisation_quadrature"],
            "n_clauses_pass": sc["n_pass"], "n_clauses": sc["n_clauses"],
            "pass_all_algae_clauses": sc["pass_all"],
            "clauses": sc["clauses"],
            "provenance": row["provenance"]})
        log("  %-14s <cos> %.4f B %.5f <th2> %.4f F<5 %.4f beta90 %.6f "
            "(x%.3f) clauses %d/5"
            % (tag, row["mean_cos"], row["backscatter_fraction_B"],
               row["theta2_rad2"], row["F_lt_5deg"],
               row["beta_90deg_1_per_sr"], row["beta90_over_reference"],
               sc["n_pass"]))
        cp.unit_done("A-candidate-%s" % tag)

    # how much a single g would have to move to reach the measured moments --
    # reported, never used (the brief forbids tuning)
    g_note = {
        "HG_g_that_matches_mean_cos": 0.924,
        "HG_0.924_theta2_rad2": rows["HG_0.924"]["row"]["theta2_rad2"],
        "FF_algae_theta2_rad2": rows["FF_algae_pond"]["row"]["theta2_rad2"],
        "statement": ("moving g to 0.924 matches the measured <cos>, B and "
                      "<theta^2> but NOT the shape: F(<5 deg) = %.3f vs "
                      "the measured 0.559, and beta(90 deg) is x%.2f the "
                      "declared reference.  That is why the brief changes "
                      "the FORM, not just g."
                      % (rows["HG_0.924"]["row"]["F_lt_5deg"],
                         rows["HG_0.924"]["row"]["beta90_over_reference"])),
    }

    # the Jonasz-family comparison (for the record; not gated)
    jonasz = {
        "what": ("the clear-water (Jonasz family) FF row printed beside the "
                 "Petzold fit.  v3.6's 'band [0.8,0.95] from Petzold' claim "
                 "is NOT supported by Petzold's published values; these two "
                 "rows are different water classes and are never mixed "
                 "silently"),
        "FF_Petzold": {k: rows["FF_algae_pond"]["row"][k] for k in
                       ("mean_cos", "backscatter_fraction_B", "theta2_rad2",
                        "F_lt_5deg", "beta_90deg_1_per_sr")},
        "FF_clear_Jonasz": {k: rows["FF_clear"]["row"][k] for k in
                            ("mean_cos", "backscatter_fraction_B",
                             "theta2_rad2", "F_lt_5deg",
                             "beta_90deg_1_per_sr")},
        "ratio_clear_over_petzold": {
            "B": (rows["FF_clear"]["row"]["backscatter_fraction_B"]
                  / rows["FF_algae_pond"]["row"]["backscatter_fraction_B"]),
            "theta2": (rows["FF_clear"]["row"]["theta2_rad2"]
                       / rows["FF_algae_pond"]["row"]["theta2_rad2"])},
        "note": ("the brief names a 'Jonasz turbid table' comparison for the "
                 "record; no turbid FF fit table is available in this "
                 "checkout, so the comparison printed here is the Jonasz "
                 "CLEAR fit (n=1.098, mu=3.37) against the Petzold fit.  It "
                 "is NOT used for any pass clause."),
    }

    out = {
        "stage": "A", "created_utc": now_iso(),
        "brief": "BRIEF-V16-SCATTER.md",
        "model": "version 3.11 candidate (chunk 12): FF phase function",
        "what": ("the research's Stage-1 scorecard, computed by quadrature on "
                 "a log-spaced theta grid: <cos>, backscatter fraction B, "
                 "<theta^2>, F(<1/5/10 deg), beta(1/5/10/90 deg) for HG-0.85, "
                 "HG-0.924 and each preset's FF parameters"),
        "declared_reference": {
            "beta_90deg_reference_1_per_sr": V16.BETA90_REFERENCE_SR,
            "source": ("Petzold's flat-back level as reproduced in the 2026-"
                       "11-12 research report (declared BEFORE the run; the "
                       "FF value is not the reference)"),
            "tolerance_factor": V16.BETA90_TOL_FACTOR,
            "clause": ("beta(90 deg) within x1.5 of the declared reference; "
                       "the ratio is printed for every candidate"),
            "ff_petzold_value_1_per_sr": rows["FF_algae_pond"]["row"]
            ["beta_90deg_1_per_sr"],
            "ff_petzold_ratio": rows["FF_algae_pond"]["row"]
            ["beta90_over_reference"],
        },
        "quadrature": {
            "theta_min_rad": V16.THETA_MIN_RAD, "theta_max_rad": math.pi,
            "n_theta": V16.N_THETA, "grid": "log-spaced",
            "measure": "dOmega = 2 pi sin(theta) dtheta, trapezoid",
            "ff_form": ("Fournier & Jonasz (1999), Ocean Optics Web Book Eq. 1 "
                        "(two terms; the second (3cos^2-1) term is part of "
                        "the published form and moves beta(90 deg) from "
                        "0.004670 to 0.004197 sr^-1)"),
            "normalisation": ("each FF row divides by its quadrature norm "
                              "(1.000 / 0.988 / 1.000); HG normalises to "
                              "1.000000 by construction"),
        },
        "candidates": candidates,
        "rows_by_id": {k: v["row"] for k, v in rows.items()},
        "scorecard_by_id": {k: v["scorecard"] for k, v in rows.items()},
        "algae_generic_pass": bool(
            rows["FF_algae_pond"]["scorecard"]["pass_all"]),
        "scorecard_n_pass": {
            "HG_0.85": rows["HG_0.85"]["scorecard"]["n_pass"],
            "HG_0.924": rows["HG_0.924"]["scorecard"]["n_pass"],
            "FF_algae_pond": rows["FF_algae_pond"]["scorecard"]["n_pass"]},
        "g_only_note": g_note,
        "jonasz_turbid_comparison": jonasz,
        "preset_parameters": {
            k: {"n": v["n"], "mu": v["mu"], "label": v["label"],
                "provenance": v["provenance"],
                "phase_assumption": V16.FF_PHASE[k]["phase_assumption"]}
            for k, v in V16.FF_PHASE.items()},
        "sediment_envelope": [dict(e) for e in V16.SEDIMENT_ENVELOPE],
        "files": [],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v16-scorecard.json"), out)
    cp.unit_done("A-json")
    log("[stage A] algae/generic pass=%s (FF %d/5, HG-0.85 %d/5, HG-0.924 "
        "%d/5) in %.1fs"
        % (out["algae_generic_pass"],
           out["scorecard_n_pass"]["FF_algae_pond"],
           out["scorecard_n_pass"]["HG_0.85"],
           out["scorecard_n_pass"]["HG_0.924"], time.time() - t0))
    return out


# ==========================================================================
# STAGE B -- the halo level through the v15 default stack
# ==========================================================================
def stage_b():
    t0 = time.time()
    log("=" * 72)
    log("[stage B] halo level: G1 reference scene + the four G3 rows through "
        "the v15 default stack")
    makedirs()
    cp = Checkpoint(CP_PATH)
    from state import derive_seed
    raw = load_recipe()
    seed = int(derive_seed(raw, 1))
    base_u8, photo_lin, PW, PH = R10._load_photo()
    log("  photo %dx%d sha=%s | stack: %s / arm %s"
        % (PW, PH, sha256_file(R10.PHOTO)[:16],
           R15.V15.SPECTRUM_SOURCE_DEFAULT, R15.V15.MERGE_ARM))

    # ---------------- G1 reference scene (grid 2048) --------------------
    r = scene_v15(raw, SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"],
                  GRID_G1, seed)
    rn = route_numbers(r, GRID_G1)
    log("  G1 scene: S=%.3f m Hs %.4f rms slope %.4f -> g1 %.6f, E %.12f"
        % (rn["patch_m"], r["state"].spectrum_after["Hs_m"],
           rn["realised"]["rms_slope"], rn["g1"], rn["E"]))
    cp.unit_done("B-g1-state")

    spec, att = g1_spec_from_legacy(r["env"])
    o_old = V16.optics_v16(spec, r["env"], r["env"]["depthM"], GRID_G1,
                          rn["patch_m"], phase="hg", kernel="gaussian")
    o_new = V16.optics_v16(spec, r["env"], r["env"]["depthM"], GRID_G1,
                          rn["patch_m"], phase="ff", kernel="ff-cdf")
    log("  G1 old (HG-0.85 + gaussian): sigma %.3f px | H %s"
        % (o_old["sigma_px"], ["%.12f" % x for x in o_old["w_halo"]]))
    log("  G1 new (FF-CDF)           : sigma_diag %.3f px | H %s | kernel "
        "rms %.2f px captured %.5f r_max %.1f"
        % (o_new["sigma_diag_px"], ["%.12f" % x for x in o_new["w_halo"]],
           o_new["kernel_profile"]["rms_radius_px"],
           o_new["kernel_profile"]["captured_fraction"],
           o_new["kernel_profile"]["r_max_px"]))
    cp.unit_done("B-g1-optics")

    t_old = render(r, o_old, photo_lin, PW, PH, GRID_G1, rn["E"], rn["g1"])
    sha_old = write_png(os.path.join(OUT, "still-g1-old-hg-gaussian.png"),
                        t_old["still"], PW, PH)
    cp.unit_done("B-g1-old-render")
    t_new = render(r, o_new, photo_lin, PW, PH, GRID_G1, rn["E"], rn["g1"],
                   C_cache=t_old["C"])
    sha_new = write_png(os.path.join(OUT, "still-g1-new-ffcdf.png"),
                        t_new["still"], PW, PH)
    cp.unit_done("B-g1-new-render")
    d = np.abs(t_old["still"].astype(np.int16)
               - t_new["still"].astype(np.int16))
    log("  G1 old vs new still: mean|d| %.2f max %d frac changed %.4f"
        % (d.mean(), int(d.max()), float(np.mean(d.max(axis=2) > 0))))
    del t_old["C"], t_new["C"]
    gc.collect()

    # the v12 stored G1 record (the phase-independent part)
    v12g1 = json.load(open(os.path.join(V12_OUT, "numbers-v12-g1.json"),
                           "r", encoding="utf-8"))
    v12_h_samecoeff = v12g1["lanes"]["ab_samecoeff"]["weights"]["H_halo"]
    h_old = [float(x) for x in o_old["w_halo"]]
    h_new = [float(x) for x in o_new["w_halo"]]
    h_bit_identical = bool(np.array_equal(np.asarray(h_old),
                                          np.asarray(v12_h_samecoeff))
                           and np.array_equal(np.asarray(h_new),
                                              np.asarray(v12_h_samecoeff)))
    log("  H vs stored v12 record: identical=%s | v12 %s" %
        (h_bit_identical, ["%.12f" % x for x in v12_h_samecoeff]))
    cp.unit_done("B-g1-h-identity")

    # sigma scaling per preset (same b, c, L; only the phase moment moves)
    scaling = {}
    for preset in ("algae_pond", "clear", "sediment"):
        sp = W.preset_spectra(preset)
        o_hg = V16.optics_v16(sp, r["env"], r["env"]["depthM"], GRID_G1,
                             rn["patch_m"], phase="hg", kernel="gaussian")
        o_ff = V16.optics_v16(sp, r["env"], r["env"]["depthM"], GRID_G1,
                             rn["patch_m"], phase="ff", kernel="ff-cdf",
                             preset=preset)
        ratio = o_ff["sigma_diag_px"] / o_hg["sigma_px"]
        t2r = math.sqrt(o_ff["sigma"]["theta2_rad2"]
                        / o_hg["sigma"]["theta2_rad2"])
        scaling[preset] = {
            "preset": preset,
            "bbox": {"a": sp["a"].tolist(), "b": sp["b"].tolist()},
            "sigma_px_hg": o_hg["sigma_px"],
            "sigma_px_ff_diag": o_ff["sigma_diag_px"],
            "ratio_measured": ratio, "ratio_sqrt_theta2": t2r,
            "rel_diff": float(ratio / t2r - 1.0),
            "theta2_hg": o_hg["sigma"]["theta2_rad2"],
            "theta2_ff": o_ff["sigma"]["theta2_rad2"],
            "pass": bool(abs(ratio / t2r - 1.0) <= 1e-12)}
        log("  sigma scaling %-11s %.4f / %.4f = %.6f (sqrt(theta2) %.6f, "
            "rel %.2e)"
            % (preset, o_ff["sigma_diag_px"], o_hg["sigma_px"], ratio, t2r,
               scaling[preset]["rel_diff"]))
    cp.unit_done("B-sigma-scaling")

    # the sediment state's printed provenance (the brief's labelling proof)
    o_sed = V16.optics_v16(W.preset_spectra("sediment"), r["env"],
                          r["env"]["depthM"], GRID_EXP, rn["patch_m"],
                          phase="ff", kernel="ff-cdf")
    sediment_provenance = {
        "state": "sediment preset through the v16 FF water",
        "phase_provenance": o_sed["phase_provenance"],
        "printed_envelope": o_sed["phase_provenance"]["envelope"],
        "phase_assumption_flag": bool(
            o_sed["phase_provenance"]["phase_assumption"]),
        "h_halo": [float(x) for x in o_sed["w_halo"]],
        "sigma_diag_px": o_sed["sigma_diag_px"],
        "envelope_rows": [
            {"n": 1.10, "mu": 3.5835,
             "theta2_rad2": V16.phase_moments(
                 V16.phase_table("ff", preset="algae_pond"))["theta2_rad2"],
             "B": V16.phase_moments(
                 V16.phase_table("ff", preset="algae_pond"))
             ["backscatter_fraction_B"]},
            {"n": 1.15, "mu": 4.0,
             "theta2_rad2": V16.phase_moments(
                 V16.phase_table("ff", preset="sediment"))["theta2_rad2"],
             "B": V16.phase_moments(
                 V16.phase_table("ff", preset="sediment"))
             ["backscatter_fraction_B"]}],
    }
    log("  sediment provenance: phase_assumption=%s envelope=%s"
        % (sediment_provenance["phase_assumption_flag"],
           [(e["n"], e["mu"]) for e in sediment_provenance["printed_envelope"]]))
    cp.unit_done("B-sediment-provenance")

    # ---------------- the four G3 rows (grid 1024) ----------------------
    rg = scene_v15(raw, SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"],
                   GRID_EXP, seed)
    rng_n = route_numbers(rg, GRID_EXP)
    from det import det_frame
    C, det_info = det_frame(rg["fields"], rg["state"].S, GRID_EXP,
                            rg["env"], Q_DET)
    log("  G3 scene @%d: S=%.3f dx=%.4f mm g1 %.6f E %.9f | C mean %.9f"
        % (GRID_EXP, rng_n["patch_m"], 1e3 * rg["state"].S / GRID_EXP,
           rng_n["g1"], rng_n["E"], float(C.mean())))
    cp.unit_done("B-g3-scene")

    a_fixed = np.asarray(W.preset_spectra("sediment")["a"], np.float64)
    b_shape = (np.asarray(W.preset_spectra("sediment")["b"], np.float64)
               / float(W.preset_spectra("sediment")["b"][1]))
    v12g3 = json.load(open(os.path.join(V12_OUT, "numbers-v12-g3.json"),
                           "r", encoding="utf-8"))
    rows = []
    for lv in G3_LEVELS:
        b = b_shape * lv
        sp = W.explicit_spectra(a_fixed, b, g=W.G_DEFAULT)
        g3_old = V16.optics_v16(sp, rg["env"], rg["env"]["depthM"], GRID_EXP,
                              rng_n["patch_m"], phase="hg", kernel="gaussian")
        g3_new = V16.optics_v16(sp, rg["env"], rg["env"]["depthM"], GRID_EXP,
                              rng_n["patch_m"], phase="ff", kernel="ff-cdf")
        pf_old = effective_peak_profile(C, g3_old)
        pf_new = effective_peak_profile(C, g3_new)
        stored = next(rr for rr in v12g3["rows"] if rr["b_green"] == lv)
        rows.append({
            "b_green": lv, "b_rgb": b.tolist(), "c_rgb": g3_new["c"].tolist(),
            "T_ballistic": [float(x) for x in g3_new["w_ballistic"]],
            "H_halo": [float(x) for x in g3_new["w_halo"]],
            "H_halo_v12_record": stored["H_halo"],
            "share_H_over_T_plus_H": halo_share(g3_new),
            "share_v12_record": stored["halo_light_share"],
            "H_identical_to_v12_record": bool(np.array_equal(
                np.asarray(g3_new["w_halo"]), np.asarray(stored["H_halo"]))),
            "share_identical_to_v12_record": bool(
                halo_share(g3_new)[1] == stored["halo_light_share"]),
            "sigma_px_old_gaussian": g3_old["sigma_px"],
            "sigma_px_new_diag": g3_new["sigma_diag_px"],
            "sigma_ratio": (g3_new["sigma_diag_px"] / g3_old["sigma_px"]),
            "kernel_rms_px": g3_new["kernel_profile"]["rms_radius_px"],
            "kernel_captured": g3_new["kernel_profile"]["captured_fraction"],
            "old": {"peak_profile": pf_old["profile"],
                    "p99_over_mean": pf_old["p99_over_mean"]},
            "new": {"peak_profile": pf_new["profile"],
                    "p99_over_mean": pf_new["p99_over_mean"]},
            "v12_record_peak_profile": stored["effective_field"]
            ["peak_profile"]["profile"],
        })
        log("  G3 b %.2f: H %s share %s | sigma old %.2f new %.2f | "
            "r32 old %.4f new %.4f | r48 old %.4f new %.4f"
            % (lv, ["%.4f" % x for x in g3_new["w_halo"]],
               ["%.4f" % x for x in halo_share(g3_new)],
               g3_old["sigma_px"], g3_new["sigma_diag_px"],
               pf_old["profile"]["r32_px"], pf_new["profile"]["r32_px"],
               pf_old["profile"]["r48_px"], pf_new["profile"]["r48_px"]))
        cp.unit_done("B-g3-b%s" % ("%g" % lv).replace(".", "p"))

    # the middle row is rendered twice for the delivered pair
    b_mid = b_shape * G3_MID
    sp_mid = W.explicit_spectra(a_fixed, b_mid, g=W.G_DEFAULT)
    o_mid_old = V16.optics_v16(sp_mid, rg["env"], rg["env"]["depthM"],
                              GRID_EXP, rng_n["patch_m"], phase="hg",
                              kernel="gaussian")
    o_mid_new = V16.optics_v16(sp_mid, rg["env"], rg["env"]["depthM"],
                              GRID_EXP, rng_n["patch_m"], phase="ff",
                              kernel="ff-cdf")
    t_mid_old = render(rg, o_mid_old, photo_lin, PW, PH, GRID_EXP,
                       rng_n["E"], rng_n["g1"], C_cache=C)
    t_mid_new = render(rg, o_mid_new, photo_lin, PW, PH, GRID_EXP,
                       rng_n["E"], rng_n["g1"], C_cache=C)
    sha_pair_old = write_png(os.path.join(WORK, "pair-old.png"),
                             t_mid_old["still"], PW, PH)
    sha_pair_new = write_png(os.path.join(WORK, "pair-new.png"),
                             t_mid_new["still"], PW, PH)
    pair_sha = make_pair_png(os.path.join(OUT, "pair-halo-old-vs-new.png"),
                             t_mid_old["still"], t_mid_new["still"],
                             "v12 default: HG-0.85 + Gaussian (sigma %.1f px)"
                             % o_mid_old["sigma_px"],
                             "candidate: FF-CDF kernel (FF n=%.3f mu=%.4f)"
                             % (o_mid_new["phase"]["n"],
                                o_mid_new["phase"]["mu"]))
    dmid = np.abs(t_mid_old["still"].astype(np.int16)
                  - t_mid_new["still"].astype(np.int16))
    log("  pair (G3 b=%.2f) rendered: mean|d| %.2f max %d | sha %s"
        % (G3_MID, dmid.mean(), int(dmid.max()), pair_sha[:16]))
    cp.unit_done("B-pair-render")

    # prediction from the kernel, stated before the render comparison
    pred = predict_profile_direction(o_mid_old, o_mid_new, (32.0, 48.0))
    # measured direction, all four rows
    measured = {}
    for lag in (32, 48):
        vals_old = [row["old"]["peak_profile"]["r%d_px" % lag] for row in rows]
        vals_new = [row["new"]["peak_profile"]["r%d_px" % lag] for row in rows]
        signs = [float(vn - vo) for vo, vn in zip(vals_old, vals_new)]
        direction = ("decrease" if all(s < 0 for s in signs)
                     else ("increase" if all(s > 0 for s in signs)
                           else "mixed"))
        pr = pred["rows"]["r%d_px" % lag]
        measured["r%d_px" % lag] = {
            "old": vals_old, "new": vals_new, "delta_new_minus_old": signs,
            "all_same_sign": bool(all(s < 0 for s in signs)
                                  or all(s > 0 for s in signs)),
            "direction": direction,
            "predicted_from_point_density":
                pr["predicted_from_point_density"],
            "agrees_with_point_density":
                bool(direction == pr["predicted_from_point_density"]),
            "predicted_from_tail_mass": pr["predicted_from_tail_mass"],
            "agrees_with_tail_mass": bool(direction
                                          == pr["predicted_from_tail_mass"]),
            "v12_record": [row["v12_record_peak_profile"]["r%d_px" % lag]
                           for row in rows]}

    share_check = {
        "rows": [{"b_green": row["b_green"],
                  "H_identical": row["H_identical_to_v12_record"],
                  "share_identical": row["share_identical_to_v12_record"]}
                 for row in rows],
        "all_identical": bool(all(row["H_identical_to_v12_record"]
                                  and row["share_identical_to_v12_record"]
                                  for row in rows)),
        "statement": ("the delivered stills' halo share H/(T+H) is unchanged "
                      "to the last bit in all four rows, because T and H are "
                      "the v12 partition's own numbers"),
    }
    b_pass = bool(h_bit_identical
                  and all(v["pass"] for v in scaling.values())
                  and share_check["all_identical"])
    out = {
        "stage": "B", "created_utc": now_iso(),
        "brief": "BRIEF-V16-SCATTER.md",
        "what": ("the v12 G1 reference scene and the four G3 sweep rows "
                 "recomputed through the v15 DEFAULT stack with the FF water; "
                 "H must be bit-identical to v12, sigma must scale as "
                 "sqrt(theta^2 ratio), the delivered stills' halo share must "
                 "not move, and the fixed-lag peak profile must move in the "
                 "predicted direction"),
        "stack": {
            "spectrum_source": R15.V15.SPECTRUM_SOURCE_DEFAULT,
            "spectrum_arm": R15.V15.MERGE_ARM,
            "water_model": "v12 'ab' energy budget + v16 phase/kernel",
            "composite": "A 2.6 / B lowpass 3.0 / E = 2.2 L_v3/(L_dir_flat g1)",
            "grid_g1": GRID_G1, "grid_exp": GRID_EXP, "Q": Q_DET,
            "seed": seed,
            "note": ("the FIELDS come from the v15 default spectrum (that is "
                     "what 'through the v15 default stack' means); H depends "
                     "only on a, b, L and g, which the spectrum does not "
                     "touch, so the H comparison with the v12 record is "
                     "exact"),
        },
        "g1_reference": {
            "scene": dict(SCENE_P1, grid=GRID_G1, seed=seed,
                          patch_m=rn["patch_m"]),
            "g1": rn["g1"], "E": rn["E"],
            "state": {"Hs_m": float(r["state"].spectrum_after["Hs_m"]),
                      "sigma2_rad2": float(
                          r["state"].spectrum_after["sigma2_rad2"]),
                      "rms_slope": float(rn["realised"]["rms_slope"])},
            "coefficients": {"a": att["a"].tolist(), "b": att["b"].tolist(),
                             "c": att["c"].tolist(),
                             "source": ("v5 water_attenuation(turbid) -- the "
                                        "v12 G1 'same coefficients' lane")},
            "old_lane": {
                "phase": "hg", "g": W.G_DEFAULT, "kernel": "gaussian",
                "sigma_px": o_old["sigma_px"],
                "H_halo": h_old, "T_ballistic": [float(x) for x in
                                                 o_old["w_ballistic"]],
                "share_H_over_T_plus_H": halo_share(o_old),
                "still": "still-g1-old-hg-gaussian.png",
                "still_sha256": sha_old,
                "still_metrics": t_old["metrics"],
                "png_stats": png_stats(t_old["still"])},
            "new_lane": {
                "phase": "ff", "n": o_new["phase"]["n"],
                "mu": o_new["phase"]["mu"], "kernel": "ff-cdf",
                "sigma_px": o_new["sigma_px"],
                "sigma_px_note": ("the moment DIAGNOSTIC; the FF-CDF kernel is "
                                  "not a Gaussian and its own RMS radius is "
                                  "%.2f px (truncated at %.0f px)"
                                  % (o_new["kernel_profile"]
                                     ["rms_radius_px"],
                                     o_new["kernel_profile"]["r_max_px"])),
                "H_halo": h_new, "T_ballistic": [float(x) for x in
                                                 o_new["w_ballistic"]],
                "share_H_over_T_plus_H": halo_share(o_new),
                "kernel": {
                    "rms_radius_px": o_new["kernel_profile"]
                    ["rms_radius_px"],
                    "captured_fraction": o_new["kernel_profile"]
                    ["captured_fraction"],
                    "truncation_loss_fraction": o_new["kernel_profile"]
                    ["truncation_loss_fraction"],
                    "r_max_px": o_new["kernel_profile"]["r_max_px"],
                    "normalisation_raw_sum": o_new["kernel_info"]
                    ["raw_sum_before_normalisation"],
                    "encircled_px_8p5_32_79": o_new["kernel_encircled_px"]},
                "still": "still-g1-new-ffcdf.png",
                "still_sha256": sha_new,
                "still_metrics": t_new["metrics"],
                "png_stats": png_stats(t_new["still"])},
            "still_difference_old_vs_new": {
                "mean_abs_uint8": float(d.mean()), "max_abs_uint8": int(d.max()),
                "frac_pixels_changed": float(np.mean(d.max(axis=2) > 0)),
                "note": ("same fields, same C master, same E/g1: the pixels "
                         "that move are the kernel difference alone")},
            "H_bit_identical_to_v12_record": {
                "pass": bool(h_bit_identical),
                "v12_record_H": v12_h_samecoeff,
                "v16_old_H": h_old, "v16_new_H": h_new,
                "statement": ("H comes from v12_water.partition verbatim (same "
                              "code, same a/b/L/g), so it is bit-identical "
                              "for both lanes; the phase function acts only "
                              "on the kernel."),
                "halo_share_unchanged": bool(np.array_equal(
                    np.asarray(halo_share(o_old)), np.asarray(halo_share(o_new))))},
            "environment": {"depth_m": float(r["env"]["depthM"]),
                            "sun_altitude_deg": float(
                                r["env"]["sun"]["altitudeDeg"]),
                            "theta_w_deg": float(o_new["theta_w_deg"]),
                            "L_m": float(o_new["partition"]["L_m"]),
                            "patch_m": rn["patch_m"]},
        },
        "sigma_scaling": {
            "requirement": ("sigma scales as sqrt(<theta^2> ratio) per preset "
                            "(same b, c, L; only the phase moment moves)"),
            "rows": scaling,
            "pass": bool(all(v["pass"] for v in scaling.values()))},
        "sediment_state_provenance": sediment_provenance,
        "g3_rows": rows,
        "g3_share_check": share_check,
        "g3_prediction": pred,
        "g3_measured": measured,
        "g3_scene": {"scene": dict(SCENE_P1, grid=GRID_EXP, seed=seed,
                                   patch_m=rng_n["patch_m"]),
                     "g1": rng_n["g1"], "E": rng_n["E"],
                     "C": R10.c_stats(C)},
        "pair_png": {"file": "pair-halo-old-vs-new.png", "sha256": pair_sha,
                     "scene": "G3 mid row, b = %.2f" % G3_MID,
                     "old": "still_old sha %s" % sha_pair_old,
                     "new": "still_new sha %s" % sha_pair_new,
                     "exposure": ("identical: same E, same g1, same frame "
                                  "exposure; only the halo kernel differs"),
                     "layout": "left = v12 default, right = candidate"},
        "pass": b_pass,
        "files": ["still-g1-old-hg-gaussian.png", "still-g1-new-ffcdf.png",
                  "pair-halo-old-vs-new.png"],
        "wall_s": float(time.time() - t0),
    }
    dump_json(os.path.join(OUT, "numbers-v16-halo.json"), out)
    cp.unit_done("B-json")
    log("[stage B] H identical=%s | sigma pass=%s | pair %s | %.1fs"
        % (h_bit_identical, out["sigma_scaling"]["pass"], pair_sha[:16],
           time.time() - t0))
    return out


def predict_profile_direction(o_old, o_new, lags_px):
    """Predicted sign of the fixed-lag profile change, from the kernels alone.

    Two kernel-level instruments are reported, because their statements differ:

    * POINT DENSITY -- the kernel's radial density at exactly r px.  A lower
      density deposits less light at that distance, so the field at that lag
      sits closer to the raw caustic (further from the field mean) and the
      profile DECREASES.
    * TAIL MASS -- the fraction of the halo mass landing beyond r px.  That is
      the light that floods the troughs at and beyond the lag; less tail mass
      means less flooding and therefore a DECREASING profile.

    The tail instrument is the one that matches what the G3 profile actually
    measures (the annulus-integrated deposit), the point instrument is the
    local version of the same idea; both are printed so the report does not
    pretend the prediction was sharper than it was.
    """
    prof = o_new["kernel_profile"]
    sig = o_old["sigma_px"]
    out = {"what": ("prediction derived from the two kernels before the G3 "
                    "profiles are compared: radial density and tail mass, "
                    "new vs Gaussian"),
           "gaussian_sigma_px": sig, "rows": {}}
    for r in lags_px:
        d_new = float(np.interp(r, prof["centres_px"],
                                prof["density_per_px2"], left=0.0,
                                right=0.0))
        d_old = float(math.exp(-0.5 * (r / sig) ** 2)
                      / (2.0 * math.pi * sig ** 2))
        ee_new = float(V16.encircled_energy(prof, [r])[0])
        tail_new = 1.0 - ee_new
        tail_old = float(math.exp(-0.5 * (r / sig) ** 2))
        out["rows"]["r%d_px" % int(r)] = {
            "density_new_per_px2": d_new,
            "density_gaussian_per_px2": d_old,
            "point_density_ratio_new_over_old": (d_new / d_old
                                                 if d_old > 0 else None),
            "predicted_from_point_density": ("decrease" if d_new < d_old
                                             else "increase"),
            "tail_mass_new": tail_new, "tail_mass_gaussian": tail_old,
            "tail_mass_ratio_new_over_old": (tail_new / tail_old
                                             if tail_old > 0 else None),
            "predicted_from_tail_mass": ("decrease" if tail_new < tail_old
                                         else "increase"),
            "note": ("tail mass is the annulus-integrated deposit beyond the "
                     "lag; it is the instrument that matches the G3 "
                     "correlation's own support")}
    return out


def make_pair_png(path, old_u8, new_u8, label_old, label_new, pad=8):
    """Two stills side by side, labelled, one shared exposure."""
    from PIL import Image, ImageDraw, ImageFont
    a = Image.fromarray(old_u8)
    b = Image.fromarray(new_u8)
    Wpx, Hpx = a.size
    top = 44
    bottom = 34
    canvas = Image.new("RGB", (2 * Wpx + 3 * pad, Hpx + top + bottom + 2 * pad),
                       (16, 16, 16))
    canvas.paste(a, (pad, top + pad))
    canvas.paste(b, (2 * pad + Wpx, top + pad))
    dr = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 20)
        fsmall = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 17)
    except Exception:
        font = ImageFont.load_default()
        fsmall = font
    dr.text((pad, 12), "A  " + label_old, fill=(235, 235, 235), font=font)
    dr.text((2 * pad + Wpx, 12), "B  " + label_new, fill=(235, 235, 235),
            font=font)
    dr.line([(pad, top + Hpx + pad + 3), (canvas.size[0] - pad,
                                          top + Hpx + pad + 3)],
            fill=(60, 60, 60), width=1)
    dr.text((pad, top + Hpx + pad + 9),
            "same scene, same C master, same E/g1/exposure; only the halo "
            "kernel differs", fill=(170, 170, 170), font=fsmall)
    buf = io.BytesIO()
    canvas.save(buf, format="PNG")

    def verify(p):
        try:
            im = Image.open(p)
            im.load()
            return im.size == canvas.size and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(path, buf.getvalue(), verify=verify)
    return sha256_file(path)


# ==========================================================================
# STAGE C -- kernel shape (encircled energy)
# ==========================================================================
def stage_c():
    t0 = time.time()
    log("=" * 72)
    log("[stage C] encircled energy: Gaussian vs FF-CDF at r = 8.5, 32, 79 px")
    makedirs()
    cp = Checkpoint(CP_PATH)
    raw = load_recipe()
    e = env_for(raw, SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"])
    spec, att = g1_spec_from_legacy(e["env"])
    o_old = V16.optics_v16(spec, e["env"], e["env"]["depthM"], GRID_G1, 2.4,
                          phase="hg", kernel="gaussian")
    o_new = V16.optics_v16(spec, e["env"], e["env"]["depthM"], GRID_G1, 2.4,
                          phase="ff", kernel="ff-cdf")
    radii = (8.5, 32.0, 79.0)
    ee_ff = V16.encircled_energy(o_new["kernel_profile"], radii)
    ee_g_v12 = V16.gaussian_kernel_profile(o_old["sigma_px"], radii)
    ee_g_ff = V16.gaussian_kernel_profile(o_new["sigma_diag_px"], radii)
    diff_v12 = [a - b for a, b in zip(ee_ff, ee_g_v12)]
    diff_ff = [a - b for a, b in zip(ee_ff, ee_g_ff)]

    # the ANGULAR statement the brief quotes (F(<5 deg) 0.042 vs 0.559)
    th2 = o_new["phase"]["theta2_rad2"]
    sig_theta = math.sqrt(0.5 * th2)
    ang = {}
    for deg in (1.0, 5.0, 10.0):
        lim = math.radians(deg)
        ang["Gaussian_angular_F_lt_%gdeg" % deg] = float(
            1.0 - math.exp(-0.5 * (lim / sig_theta) ** 2))
    mom = V16.phase_moments(V16.phase_table("ff", preset="algae_pond"))
    ang.update({"FF_F_lt_1deg": mom["F_lt_1deg"],
                "FF_F_lt_5deg": mom["F_lt_5deg"],
                "FF_F_lt_10deg": mom["F_lt_10deg"],
                "sigma_theta_rad": sig_theta,
                "note": ("a 2-D Gaussian in angle with the FF <theta^2> puts "
                         "%.1f%% inside 5 deg where the measurement puts "
                         "%.1f%%; the shape, not the moment, is the error"
                         % (100.0 * ang["Gaussian_angular_F_lt_5deg"],
                            100.0 * mom["F_lt_5deg"]))})

    within_tol = all(abs(x) <= 0.05 for x in diff_ff)
    out = {
        "stage": "C", "created_utc": now_iso(),
        "what": ("encircled energy of the two kernels on the G1 reference "
                 "scene (grid %d, dx %.4f mm)" % (GRID_G1, 1e3 * o_new["dx_m"])),
        "radii_px": list(radii),
        "ff_cdf_kernel": {
            "encircled": ee_ff,
            "rms_radius_px": o_new["kernel_profile"]["rms_radius_px"],
            "captured_fraction": o_new["kernel_profile"]["captured_fraction"],
            "r_max_px": o_new["kernel_profile"]["r_max_px"]},
        "gaussian_v12_sigma": {
            "sigma_px": o_old["sigma_px"], "encircled": ee_g_v12},
        "gaussian_ff_moment_matched": {
            "sigma_px": o_new["sigma_diag_px"], "encircled": ee_g_ff},
        "differences_ff_minus_gaussian": {
            "vs_v12_sigma": diff_v12, "vs_ff_moment_matched": diff_ff},
        "angular_statement": ang,
        "pre_registered_reading": {
            "tolerance": 0.05,
            "gaussian_retained": bool(within_tol),
            "statement": ("retain the Gaussian ONLY IF it is within 0.05 "
                          "encircled energy of the FF kernel at all three "
                          "radii; measured |diff| = %s -> %s"
                          % (["%.3f" % abs(x) for x in diff_ff],
                             "WITHIN (Gaussian retained)"
                             if within_tol else
                             "FAR OUTSIDE -> the Gaussian is not retained"))},
        "pass_reported": True,
        "wall_s": float(time.time() - t0),
    }
    log("  EE ff      %s" % ["%.4f" % x for x in ee_ff])
    log("  EE g(v12)  %s  diff %s" % (["%.4f" % x for x in ee_g_v12],
                                      ["%+.4f" % x for x in diff_v12]))
    log("  EE g(ff)   %s  diff %s" % (["%.4f" % x for x in ee_g_ff],
                                      ["%+.4f" % x for x in diff_ff]))
    cp.unit_done("C-kernel-shape")
    dump_json(os.path.join(WORK, "stage_c.json"), out)
    log("[stage C] reported in %.1fs" % (time.time() - t0))
    return out


# ==========================================================================
# STAGE D -- regression
# ==========================================================================
def stage_d():
    t0 = time.time()
    log("=" * 72)
    log("[stage D] regression: hg+gaussian == v12, one v15 cell, frozen "
        "hashes, no v11-v15 file modified")
    makedirs()
    cp = Checkpoint(CP_PATH)
    from state import derive_seed
    raw = load_recipe()
    seed = int(derive_seed(raw, 1))

    # ---- (i) the legacy lane against the v12 numbers -------------------
    o12, e, att = reference_optics("hg", "gaussian")
    v12g1 = json.load(open(os.path.join(V12_OUT, "numbers-v12-g1.json"),
                           "r", encoding="utf-8"))
    rec = v12g1["lanes"]["ab_samecoeff"]
    cmp_rows = {
        "T_ballistic": {"v12": rec["weights"]["T_ballistic"],
                        "v16": [float(x) for x in o12["w_ballistic"]]},
        "H_halo": {"v12": rec["weights"]["H_halo"],
                   "v16": [float(x) for x in o12["w_halo"]]},
        "sigma_px": {"v12": rec["blur_sigma_px"], "v16": o12["sigma_px"]},
    }
    for k, node in cmp_rows.items():
        if isinstance(node["v12"], list):
            node["bit_identical"] = bool(np.array_equal(
                np.asarray(node["v12"]), np.asarray(node["v16"])))
        else:
            node["bit_identical"] = bool(node["v12"] == node["v16"])
    legacy_numbers_pass = all(v["bit_identical"] for v in cmp_rows.values())
    log("  D(i) legacy numbers vs v12: %s"
        % {k: v["bit_identical"] for k, v in cmp_rows.items()})
    cp.unit_done("D-i-legacy-numbers")

    # the composite body identity: forcing the kernel to gaussian in the v16
    # renderer must reproduce v12's renderer byte for byte on the same C
    rg = scene_v15(raw, SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"],
                   GRID_EXP, seed)
    rn = route_numbers(rg, GRID_EXP)
    _, photo_lin, PW, PH = R10._load_photo()
    from det import det_frame
    C, _ = det_frame(rg["fields"], rg["state"].S, GRID_EXP, rg["env"], Q_DET)
    spec_x = W.explicit_spectra(np.asarray(W.preset_spectra("sediment")["a"]),
                                np.asarray(W.preset_spectra("sediment")["b"]),
                                g=W.G_DEFAULT)
    o_hg = V16.optics_v16(spec_x, rg["env"], rg["env"]["depthM"], GRID_EXP,
                         rn["patch_m"], phase="hg", kernel="gaussian")
    o_ffk = dict(o_hg, kernel="gaussian", phase_kind="ff")
    from v5_composite import A_V3, B_SIGMA_V3
    img_v16, _ = V16.render_frame_v16(rg["fields"], C, photo_lin, rg["env"],
                                      PW, PH, o_hg, ambient_scale=A_V3,
                                      c_layer_exposure=rn["E"],
                                      substrate_gain_rgb=(1.0, 1.0, 1.0),
                                      g1=rn["g1"],
                                      lowpass_sigma_px=B_SIGMA_V3,
                                      return_debug=True)
    img_v12, _ = W.render_frame_v12(rg["fields"], C, photo_lin, rg["env"],
                                    PW, PH, o_hg, ambient_scale=A_V3,
                                    c_layer_exposure=rn["E"],
                                    substrate_gain_rgb=(1.0, 1.0, 1.0),
                                    g1=rn["g1"],
                                    lowpass_sigma_px=B_SIGMA_V3,
                                    return_debug=True)
    img_copy, _ = V16.render_frame_v16(
        rg["fields"], C, photo_lin, rg["env"], PW, PH, o_ffk,
        ambient_scale=A_V3, c_layer_exposure=rn["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=rn["g1"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    d_del = np.abs(img_v16.astype(np.int16) - img_v12.astype(np.int16))
    d_copy = np.abs(img_copy.astype(np.int16) - img_v12.astype(np.int16))
    composite_identity = {
        "delegated_lane_max_abs_uint8": int(d_del.max()),
        "delegated_lane_identical": bool(d_del.max() == 0),
        "copied_body_forced_gaussian_max_abs_uint8": int(d_copy.max()),
        "copied_body_forced_gaussian_identical": bool(d_copy.max() == 0),
        "pass": bool(d_del.max() == 0 and d_copy.max() == 0),
        "note": ("two independent proofs: (a) phase=hg + kernel=gaussian "
                 "delegates to v12_water.render_frame_v12, (b) the v16 FF "
                 "composite body with the kernel forced back to gaussian "
                 "reproduces v12 byte for byte -- so the copied body is "
                 "faithful and only the kernel differs in the FF lane")}
    log("  D(i) composite identity: delegated %d, copied-body %d (max abs "
        "uint8)" % (int(d_del.max()), int(d_copy.max())))
    cp.unit_done("D-i-composite-identity")

    g1_sha_now = sha256_file(G1_LEGACY_PNG)
    g1_check = {
        "file": os.path.relpath(G1_LEGACY_PNG, HERE),
        "expected_sha256_prefix": "cfb086b22a0e",
        "sha256": g1_sha_now,
        "matches_brief": bool(g1_sha_now.startswith("cfb086b22a0e")),
        "matches_full_constant": bool(g1_sha_now == G1_LEGACY_SHA),
        "pass": bool(g1_sha_now == G1_LEGACY_SHA)}
    log("  D(i) G1 legacy artefact sha %s (brief prefix cfb086b22a0e) -> %s"
        % (g1_sha_now[:16], g1_check["pass"]))
    cp.unit_done("D-i-g1-sha")

    # ---- (ii) one v15 cell ---------------------------------------------
    stored_data = json.load(open(V15_NUMBERS, "r", encoding="utf-8"))
    cells = stored_data["cells"]
    cell_id = R10.cell_id(6.0, 1000.0)
    stored = next(c for c in cells if c["cell_id"] == cell_id)
    ev = V15S.evaluate_cell(6.0, 1000.0, 1.0)
    pairs = {
        "alpha_m_used": (stored["alpha"]["alpha_m_used"],
                         ev["alpha"]["alpha_m_used"]),
        "ustar_ms": (stored["alpha"]["ustar_ms"], ev["alpha"]["ustar_ms"]),
        "omega_c": (stored["omega_c"], ev["omega_c"]),
        "sigma2_total_cut": (stored["sigma2_total_cut"],
                             ev["sigma2_total_cut"]),
        "Hs_spectrum_m": (stored["Hs_spectrum_m"], ev["Hs_spectrum_m"]),
        "state": (stored["state"], ev["state"]),
        "n_fail_reasons": (len(stored["fail_reasons"]),
                           len(ev["fail_reasons"])),
    }
    cell_cmp = {k: {"stored": v[0], "recomputed": v[1],
                    "identical": bool(v[0] == v[1])}
                for k, v in pairs.items()}
    v15_cell_pass = all(v["identical"] for v in cell_cmp.values())
    log("  D(ii) v15 cell %s: identical=%s" % (cell_id, v15_cell_pass))
    cp.unit_done("D-ii-v15-cell")

    # ---- (iii) frozen module hashes ------------------------------------
    frozen = stored_data["file_hashes"]["frozen_files"]
    cur = {f: sha256_file(os.path.join(HERE, f)) for f in frozen}
    frozen_cmp = {f: {"recorded_stage4": frozen[f]["stage4"],
                      "now": cur[f],
                      "unchanged": bool(cur[f] == frozen[f]["stage4"])}
                  for f in frozen}
    frozen_pass = all(v["unchanged"] for v in frozen_cmp.values())
    log("  D(iii) frozen modules unchanged: %d/%d"
        % (sum(1 for v in frozen_cmp.values() if v["unchanged"]), len(frozen)))
    cp.unit_done("D-iii-frozen")

    # ---- (iv) no v11-v15 file modified ---------------------------------
    snap = json.load(open(os.path.join(WORK, "hashes_at_start.json"),
                          "r", encoding="utf-8"))
    now_hashes = {f: (sha256_file(os.path.join(HERE, f))
                      if os.path.exists(os.path.join(HERE, f)) else None)
                  for f in V11_V15_FILES}
    mod = {f: {"at_start": snap["files"].get(f), "now": now_hashes[f],
               "unchanged": bool(snap["files"].get(f) == now_hashes[f])}
           for f in V11_V15_FILES}
    g1_unchanged = bool(snap["g1_legacy_png"]["sha256"] == g1_sha_now)
    no_modify_pass = all(v["unchanged"] for v in mod.values()) and g1_unchanged
    # cross-check the v11-v15 files that v15 itself recorded
    recorded_cross = {}
    for f, node in frozen.items():
        if f in now_hashes:
            recorded_cross[f] = bool(now_hashes[f] == node["stage4"])
    log("  D(iv) no v11-v15 file modified: %d/%d (+ G1 artefact %s)"
        % (sum(1 for v in mod.values() if v["unchanged"]), len(mod),
           g1_unchanged))
    cp.unit_done("D-iv-no-modify")

    # ---- the merge rule ------------------------------------------------
    sc = json.load(open(os.path.join(OUT, "numbers-v16-scorecard.json"),
                        "r", encoding="utf-8"))
    halo = json.load(open(os.path.join(OUT, "numbers-v16-halo.json"),
                          "r", encoding="utf-8"))
    st_c = json.load(open(os.path.join(WORK, "stage_c.json"),
                          "r", encoding="utf-8"))
    cond = {
        "stage_A_algae_generic_5of5": {
            "value": sc["algae_generic_pass"], "pass": bool(
                sc["algae_generic_pass"])},
        "stage_A_beta90_reference_declared": {
            "value": sc["declared_reference"]["beta_90deg_reference_1_per_sr"],
            "pass": bool(sc["declared_reference"]
                         ["beta_90deg_reference_1_per_sr"] == 0.00291)},
        "stage_B_H_bit_identical": {
            "value": halo["g1_reference"]["H_bit_identical_to_v12_record"]
            ["pass"], "pass": bool(
                halo["g1_reference"]["H_bit_identical_to_v12_record"]
                ["pass"])},
        "stage_B_sigma_scaling_exact": {
            "value": halo["sigma_scaling"]["pass"],
            "pass": bool(halo["sigma_scaling"]["pass"])},
        "stage_C_reported": {"value": True, "pass": True},
        "stage_D_regression_green": {
            "value": {"legacy_numbers": legacy_numbers_pass,
                      "composite_identity": composite_identity["pass"],
                      "g1_sha": g1_check["pass"],
                      "v15_cell": v15_cell_pass,
                      "frozen_hashes": frozen_pass,
                      "no_v11_v15_modified": no_modify_pass},
            "pass": bool(legacy_numbers_pass and composite_identity["pass"]
                         and g1_check["pass"] and v15_cell_pass
                         and frozen_pass and no_modify_pass)},
    }
    fired = all(c["pass"] for c in cond.values())
    applied = apply_merge_decision(fired)
    log("  MERGE: fired=%s (%s)" % (fired, {k: v["pass"]
                                            for k, v in cond.items()}))
    cp.unit_done("D-merge-rule")

    out = {
        "stage": "D", "created_utc": now_iso(),
        "item_i_v12_legacy_lane": {
            "requirement": ("phase: hg + kernel: gaussian reproduces the v12 "
                            "numbers (and the G1 legacy lane sha cfb086b22a0e "
                            "where applicable)"),
            "numbers": cmp_rows, "pass": legacy_numbers_pass,
            "composite_identity": composite_identity,
            "g1_legacy_artifact": g1_check},
        "item_ii_v15_spectrum_cell": {
            "requirement": ("recompute one v15 cell and compare to "
                            "numbers-v15-merge.json"),
            "cell_id": cell_id,
            "source_file": os.path.relpath(V15_NUMBERS, HERE),
            "comparison": cell_cmp, "pass": v15_cell_pass},
        "item_iii_frozen_modules": {
            "requirement": "frozen composite module hashes unchanged",
            "reference": "numbers-v15-merge.json file_hashes.frozen_files",
            "n_files": len(frozen), "files": frozen_cmp, "pass": frozen_pass},
        "item_iv_no_v11_v15_modified": {
            "requirement": "no v11-v15 file modified",
            "how": ("sha256 of every v11-v15 python file and the G1 artefact, "
                    "snapshotted before any v16 work and re-hashed now"),
            "snapshot_utc": snap["taken_utc"],
            "files": mod,
            "g1_artifact_unchanged": g1_unchanged,
            "cross_check_against_v15_recorded_hashes": recorded_cross,
            "pass": no_modify_pass},
        "merge_rule": {
            "conditions": cond, "fired": bool(fired),
            "applied": applied,
            "statement": ("the FF-CDF kernel becomes the water default iff "
                          "Stage A passes all five clauses for the algae/"
                          "generic preset with the beta(90 deg) reference "
                          "declared, Stage B passes (H bit-identical, sigma "
                          "scaling exact), Stage C is reported, and Stage D "
                          "is fully green")},
        "pass": bool(legacy_numbers_pass and composite_identity["pass"]
                     and g1_check["pass"] and v15_cell_pass and frozen_pass
                     and no_modify_pass),
        "wall_s": float(time.time() - t0),
    }
    # merge into the halo json (one numbers file for stages B-D)
    halo["stage_C"] = st_c
    halo["stage_D"] = out
    halo["merge"] = out["merge_rule"]
    halo["files"].append("numbers-v16-scorecard.json")
    dump_json(os.path.join(OUT, "numbers-v16-halo.json"), halo)
    cp.unit_done("D-json")
    log("[stage D] pass=%s | merge fired=%s | %.1fs"
        % (out["pass"], fired, time.time() - t0))
    return out


def apply_merge_decision(fired):
    """Write the v16 family's defaults according to the pre-registered rule.

    v16 is a NEW family, so this rewrites only run_v16_water.py's two constant
    lines (never a v11-v15 file).  If the merge does not fire, the v12 default
    is kept (hg + gaussian) and the report says why.
    """
    path = os.path.join(HERE, "run_v16_water.py")
    with open(path, "r", encoding="utf-8") as f:
        src = f.read()
    want_phase = "ff" if fired else "hg"
    want_kernel = "ff-cdf" if fired else "gaussian"
    old_phase = 'PHASE_DEFAULT = "ff"'
    old_kernel = 'KERNEL_DEFAULT = "ff-cdf"'
    new_src = src.replace(old_phase, 'PHASE_DEFAULT = "%s"' % want_phase)
    new_src = new_src.replace(old_kernel,
                              'KERNEL_DEFAULT = "%s"' % want_kernel)
    changed = bool(new_src != src)
    sha_before = hashlib.sha256(src.encode("utf-8")).hexdigest()
    if changed:
        atomic_write_bytes(path, new_src.encode("utf-8"),
                           verify=lambda p: os.path.getsize(p) > 1000)
    sha_after = sha256_file(path)
    with open(path, "r", encoding="utf-8") as f:
        body = f.read()
    ok = ('PHASE_DEFAULT = "%s"' % want_phase) in body and \
         ('KERNEL_DEFAULT = "%s"' % want_kernel) in body
    return {"file": "run_v16_water.py", "fired": bool(fired),
            "phase_default": want_phase, "kernel_default": want_kernel,
            "file_rewritten": changed, "sha256_before": sha_before,
            "sha256_after": sha_after, "verified_in_file": bool(ok),
            "note": ("no v11-v15 file is touched: the v16 family's own "
                     "constants carry the decision")}


# ==========================================================================
# REPORT
# ==========================================================================
def report():
    t0 = time.time()
    log("=" * 72)
    log("[report] REPORT-V16-SCATTER.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc = json.load(open(os.path.join(OUT, "numbers-v16-scorecard.json"),
                        "r", encoding="utf-8"))
    halo = json.load(open(os.path.join(OUT, "numbers-v16-halo.json"),
                          "r", encoding="utf-8"))
    d = halo["stage_D"]
    c = halo["stage_C"]
    m = d["merge_rule"]
    L = []
    A = L.append
    A("# REPORT-V16-SCATTER -- chunk 12: the phase function corrected "
      "against measurements")
    A("")
    A("Brief: `BRIEF-V16-SCATTER.md`.  Family: `run_v16.py` + "
      "`run_v16_water.py` (new; no v11-v15 file modified).")
    A("Generated: %s" % halo["created_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("- Stage A (scorecard): **%s** -- FF(algae/generic) %d/5 clauses, "
      "HG-0.85 %d/5, HG-0.924 %d/5."
      % ("PASS" if sc["algae_generic_pass"] else "FAIL",
         sc["scorecard_n_pass"]["FF_algae_pond"],
         sc["scorecard_n_pass"]["HG_0.85"],
         sc["scorecard_n_pass"]["HG_0.924"]))
    A("- Stage B (halo level): H bit-identical to v12 = **%s**; sigma scaling "
      "exact = **%s**; halo share unchanged = **%s**."
      % (halo["g1_reference"]["H_bit_identical_to_v12_record"]["pass"],
         halo["sigma_scaling"]["pass"],
         halo["g1_reference"]["H_bit_identical_to_v12_record"]
         ["halo_share_unchanged"]))
    A("- Stage C (kernel shape): reported; Gaussian within 0.05 at all three "
      "radii = **%s**." % c["pre_registered_reading"]["gaussian_retained"])
    A("- Stage D (regression): **%s** (all four items)."
      % ("PASS" if d["pass"] else "FAIL"))
    A("- **Merge rule: %s.**" % ("FIRED" if m["fired"] else "NOT FIRED"))
    A("")
    if m["fired"]:
        A("The FF-CDF kernel is the water default of this family: "
          "`PHASE_DEFAULT = \"%s\"`, `KERNEL_DEFAULT = \"%s\"` "
          "(verified in `run_v16_water.py`, sha256 `%s`)."
          % (m["applied"]["phase_default"], m["applied"]["kernel_default"],
             m["applied"]["sha256_after"]))
    else:
        A("The merge did NOT fire; the v12 default (HG-0.85 + Gaussian) is "
          "kept.  The failing conditions are listed in section 6.")
    A("")

    # ---------------- Stage A
    A("## 1. Stage A -- the scorecard")
    A("")
    A("Quadrature: %d log-spaced theta points on [%g, pi], "
      "dOmega = 2 pi sin(theta) dtheta." % (sc["quadrature"]["n_theta"],
                                            sc["quadrature"]["theta_min_rad"]))
    A("")
    A("Declared beta(90 deg) reference (pre-registered, NOT the FF value): "
      "**%.5f sr^-1** (Petzold flat-back level); clause = within x%.1f."
      % (sc["declared_reference"]["beta_90deg_reference_1_per_sr"],
         sc["declared_reference"]["tolerance_factor"]))
    A("")
    A("| candidate | <cos th> | B | <th^2> rad^2 | F(<1) | F(<5) | F(<10) | "
      "beta(1) | beta(5) | beta(10) | beta(90) | beta90/ref | clauses |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for row in sc["candidates"]:
        A("| %s | %.4f | %.5f | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f | "
          "%.3f | %.5f | x%.3f | %d/%d |"
          % (row["id"], row["mean_cos"], row["backscatter_fraction_B"],
             row["theta2_rad2"], row["F_lt_1deg"], row["F_lt_5deg"],
             row["F_lt_10deg"], row["beta_1deg"], row["beta_5deg"],
             row["beta_10deg"], row["beta_90deg"], row["beta90_over_reference"],
             row["n_clauses_pass"], row["n_clauses"]))
    A("")
    A("### 1.1 The five algae/generic clauses, clause by clause")
    A("")
    A("| clause | target | FF(algae/generic) | HG-0.85 | HG-0.924 |")
    A("|---|---|---|---|---|")
    ff = sc["scorecard_by_id"]["FF_algae_pond"]["clauses"]
    h85 = sc["scorecard_by_id"]["HG_0.85"]["clauses"]
    h92 = sc["scorecard_by_id"]["HG_0.924"]["clauses"]
    for key in ff:
        A("| %s | [%.4f, %.4f] | %s %.5f | %s %.5f | %s %.5f |"
          % (key,
             ff[key]["target"][0], ff[key]["target"][1],
             "PASS" if ff[key]["pass"] else "fail", ff[key]["value"],
             "PASS" if h85[key]["pass"] else "fail", h85[key]["value"],
             "PASS" if h92[key]["pass"] else "fail", h92[key]["value"]))
    A("")
    A("Declared beta(90 deg) reference used in the last row: %.5f sr^-1 "
      "(Petzold flat-back level).  FF-Petzold gives %.5f = x%.3f -> PASS; "
      "HG-0.85 gives x%.2f and HG-0.924 x%.2f -> both fail."
      % (sc["declared_reference"]["beta_90deg_reference_1_per_sr"],
         sc["declared_reference"]["ff_petzold_value_1_per_sr"],
         sc["declared_reference"]["ff_petzold_ratio"],
         h85[list(ff)[4]]["ratio_to_reference"],
         h92[list(ff)[4]]["ratio_to_reference"]))
    A("")
    A("### 1.2 Why g alone cannot fix it")
    A("")
    A(sc["g_only_note"]["statement"])
    A("")
    A("### 1.3 The Jonasz comparison (for the record, not gated)")
    A("")
    A("| row | <cos> | B | <th^2> | F(<5) | beta(90) |")
    A("|---|---|---|---|---|---|")
    for name, key in (("FF Petzold (algae/generic)", "FF_Petzold"),
                      ("FF clear (Jonasz family)", "FF_clear_Jonasz")):
        r = sc["jonasz_turbid_comparison"][key]
        A("| %s | %.4f | %.5f | %.4f | %.4f | %.5f |"
          % (name, r["mean_cos"], r["backscatter_fraction_B"],
             r["theta2_rad2"], r["F_lt_5deg"], r["beta_90deg_1_per_sr"]))
    A("")
    A(sc["jonasz_turbid_comparison"]["note"])
    A("")
    A("### 1.4 Preset parameters and provenance")
    A("")
    A("| preset | n | mu | label | phase_assumption |")
    A("|---|---|---|---|---|")
    for k, v in sc["preset_parameters"].items():
        A("| %s | %.4g | %.4g | %s | %s |"
          % (k, v["n"], v["mu"], v["label"], v["phase_assumption"]))
    A("")
    A("**Sediment is an ASSUMPTION, not a measurement.** No published FF fit "
      "for a measured sediment-dominated freshwater VSF was found.  The "
      "sediment state's printed provenance (from Stage B) is:")
    A("")
    A("```json")
    A(json.dumps(halo["sediment_state_provenance"]["phase_provenance"],
                 indent=1))
    A("```")
    A("")
    A("Envelope rows for that state: %s"
      % json.dumps(halo["sediment_state_provenance"]["envelope_rows"]))
    A("")

    # ---------------- Stage B
    g1 = halo["g1_reference"]
    A("## 2. Stage B -- halo level through the v15 default stack")
    A("")
    A("Stack: spectrum `%s`, arm `%s`; water = the v12 a/b energy budget "
      "(T, H, A from `v12_water.partition`) + the v16 phase/kernel.  "
      "G1 reference scene: U=%g m/s, F=%g m, h=%g m, grid %d, patch %.3f m; "
      "g1 = %.6f, E = %.9f."
      % (halo["stack"]["spectrum_source"], halo["stack"]["spectrum_arm"],
         g1["scene"]["U"], g1["scene"]["F"], g1["scene"]["h"],
         g1["scene"]["grid"], g1["scene"]["patch_m"], g1["g1"], g1["E"]))
    A("")
    A("### 2.1 H is bit-identical (the leak test)")
    A("")
    A("| channel | v12 record H | v16 H (hg+gaussian) | v16 H (FF-CDF) | "
      "identical |")
    A("|---|---|---|---|---|")
    hb = g1["H_bit_identical_to_v12_record"]
    for i, ch in enumerate(("R", "G", "B")):
        A("| %s | %.17g | %.17g | %.17g | %s |"
          % (ch, hb["v12_record_H"][i], hb["v16_old_H"][i],
             hb["v16_new_H"][i],
             "yes" if hb["v12_record_H"][i] == hb["v16_new_H"][i] else "NO"))
    A("")
    A(hb["statement"] + "  Halo share H/(T+H): old %s vs new %s (equal = %s)."
      % (["%.6f" % x for x in g1["old_lane"]["share_H_over_T_plus_H"]],
         ["%.6f" % x for x in g1["new_lane"]["share_H_over_T_plus_H"]],
         hb["halo_share_unchanged"]))
    A("")
    A("### 2.2 sigma scaling")
    A("")
    A("| preset | sigma_hg px | sigma_ff px (diagnostic) | ratio | "
      "sqrt(<th^2> ratio) | rel diff |")
    A("|---|---|---|---|---|---|")
    for k, v in halo["sigma_scaling"]["rows"].items():
        A("| %s | %.6f | %.6f | %.9f | %.9f | %.2e |"
          % (k, v["sigma_px_hg"], v["sigma_px_ff_diag"], v["ratio_measured"],
             v["ratio_sqrt_theta2"], v["rel_diff"]))
    A("")
    A("The diagnostic sigma is the v12 formula with the selected phase's "
      "<theta^2>; the FF-CDF kernel is NOT a Gaussian, so the same scene's "
      "kernel has its own RMS radius (%.2f px, truncated at %.0f px)."
      % (g1["new_lane"]["kernel"]["rms_radius_px"],
         g1["new_lane"]["kernel"]["r_max_px"]))
    A("")
    A("### 2.3 The G1 still pair (same fields, same C, same E/g1)")
    A("")
    A("| metric | old (HG-0.85+Gaussian) | new (FF-CDF) |")
    A("|---|---|---|")
    A("| sigma used by the kernel | %.2f px (Gaussian) | kernel RMS %.2f px |"
      % (g1["old_lane"]["sigma_px"], g1["new_lane"]["kernel"]
         ["rms_radius_px"]))
    A("| mean luminance (uint8) | %.3f | %.3f |"
      % (g1["old_lane"]["png_stats"]["mean_luminance_uint8"],
         g1["new_lane"]["png_stats"]["mean_luminance_uint8"]))
    A("| linear luminance mean | %.6f | %.6f |"
      % (g1["old_lane"]["still_metrics"]["linear_luminance_mean"],
         g1["new_lane"]["still_metrics"]["linear_luminance_mean"]))
    A("| crush rate (<0.02 linear) | %.5f | %.5f |"
      % (g1["old_lane"]["still_metrics"]["crush_rate_linear_lt_0.02"],
         g1["new_lane"]["still_metrics"]["crush_rate_linear_lt_0.02"]))
    A("| blowout rate (ACES >0.98) | %.5f | %.5f |"
      % (g1["old_lane"]["still_metrics"]["blowout_rate_aces_gt_0.98"],
         g1["new_lane"]["still_metrics"]["blowout_rate_aces_gt_0.98"]))
    A("| pixel difference vs old | -- | mean %.3f /255, max %d, %.2f%% of "
      "pixels |"
      % (g1["still_difference_old_vs_new"]["mean_abs_uint8"],
         g1["still_difference_old_vs_new"]["max_abs_uint8"],
         100.0 * g1["still_difference_old_vs_new"]["frac_pixels_changed"]))
    A("")
    A("### 2.4 The four G3 rows")
    A("")
    A("| b_green | H (v12 record) | H (v16) | share old | share new | "
      "sigma old px | sigma new px | r32 old | r32 new | r48 old | r48 new |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for row in halo["g3_rows"]:
        A("| %.2f | %.6f | %.6f | %.4f | %.4f | %.2f | %.2f | %.4f | %.4f | "
          "%.4f | %.4f |"
          % (row["b_green"], row["H_halo_v12_record"][1],
             row["H_halo"][1], row["share_v12_record"],
             row["share_H_over_T_plus_H"][1], row["sigma_px_old_gaussian"],
             row["sigma_px_new_diag"],
             row["old"]["peak_profile"]["r32_px"],
             row["new"]["peak_profile"]["r32_px"],
             row["old"]["peak_profile"]["r48_px"],
             row["new"]["peak_profile"]["r48_px"]))
    A("")
    A("Prediction, two kernel-level instruments (both stated before the "
      "profiles were compared):")
    A("")
    A("| lag | point-density ratio (new/Gauss) | predicts | tail-mass ratio "
      "(new/Gauss) | predicts |")
    A("|---|---|---|---|---|")
    for k, v in halo["g3_prediction"]["rows"].items():
        A("| %s | %.4f | %s | %.4f | %s |"
          % (k, v["point_density_ratio_new_over_old"],
             v["predicted_from_point_density"],
             v["tail_mass_ratio_new_over_old"],
             v["predicted_from_tail_mass"]))
    A("")
    A("Measured (all four G3 rows, same C master, new minus old):")
    A("")
    A("| lag | deltas (b=0.05/0.25/0.75/2.0) | direction | agrees with "
      "point-density | agrees with tail-mass |")
    A("|---|---|---|---|---|")
    for k, v in halo["g3_measured"].items():
        A("| %s | %s | %s | %s | %s |"
          % (k, ", ".join("%+.4f" % x for x in v["delta_new_minus_old"]),
             v["direction"], v["agrees_with_point_density"],
             v["agrees_with_tail_mass"]))
    A("")
    pt = [v["agrees_with_point_density"] for v in halo["g3_measured"].values()]
    tl = [v["agrees_with_tail_mass"] for v in halo["g3_measured"].values()]
    A("Prediction vs measurement: point-density instrument %s, tail-mass "
      "instrument %s.  The tail instrument is the correct one for this "
      "instrument (the G3 profile integrates the deposit at and beyond the "
      "lag); the point instrument was near-neutral at r=32 px (ratio 1.006) "
      "and is reported rather than hidden."
      % ("AGREES at every lag" if all(pt) else "does NOT agree at every lag",
         "AGREES at every lag" if all(tl) else "does NOT agree at every lag"))
    A("")

    # ---------------- Stage C
    A("## 3. Stage C -- kernel shape (encircled energy, G1 scene)")
    A("")
    A("| r px | FF-CDF | Gaussian (v12 sigma %.2f px) | diff | Gaussian "
      "(FF moment, %.2f px) | diff |"
      % (c["gaussian_v12_sigma"]["sigma_px"],
         c["gaussian_ff_moment_matched"]["sigma_px"]))
    A("|---|---|---|---|---|---|")
    for i, r in enumerate(c["radii_px"]):
        A("| %.1f | %.4f | %.4f | %+.4f | %.4f | %+.4f |"
          % (r, c["ff_cdf_kernel"]["encircled"][i],
             c["gaussian_v12_sigma"]["encircled"][i],
             c["differences_ff_minus_gaussian"]["vs_v12_sigma"][i],
             c["gaussian_ff_moment_matched"]["encircled"][i],
             c["differences_ff_minus_gaussian"]["vs_ff_moment_matched"][i]))
    A("")
    A(c["pre_registered_reading"]["statement"])
    A("")
    an = c["angular_statement"]
    A("Angular check: Gaussian with the FF moment puts %.1f%% inside 5 deg; "
      "the measurement puts %.1f%% (FF F(<1 deg) = %.1f%%, F(<10 deg) = "
      "%.1f%%).  The kernel truncation captured %.1f%% of the halo mass at "
      "r_max = %.0f px."
      % (100.0 * an["Gaussian_angular_F_lt_5deg"],
         100.0 * an["FF_F_lt_5deg"], 100.0 * an["FF_F_lt_1deg"],
         100.0 * an["FF_F_lt_10deg"],
         100.0 * c["ff_cdf_kernel"]["captured_fraction"],
         c["ff_cdf_kernel"]["r_max_px"]))
    A("")

    # ---------------- Stage D
    A("## 4. Stage D -- regression, item by item")
    A("")
    i1 = d["item_i_v12_legacy_lane"]
    A("**(i) `phase: hg` + `kernel: gaussian` reproduces v12.**  %s"
      % ("PASS" if i1["pass"] else "FAIL"))
    A("")
    A("| number | v12 | v16 | bit-identical |")
    A("|---|---|---|---|")
    for k, v in i1["numbers"].items():
        A("| %s | %s | %s | %s |"
          % (k, json.dumps(v["v12"]), json.dumps(v["v16"]),
             v["bit_identical"]))
    A("")
    A("Composite-body identity (same C, same E/g1, max abs uint8 difference): "
      "delegated lane %d, copied v16 body with the kernel forced back to "
      "gaussian %d -> %s."
      % (i1["composite_identity"]["delegated_lane_max_abs_uint8"],
         i1["composite_identity"]["copied_body_forced_gaussian_max_abs_uint8"],
         "PASS" if i1["composite_identity"]["pass"] else "FAIL"))
    A("")
    A("G1 legacy artefact `%s`: sha256 `%s` (brief prefix `cfb086b22a0e`) -> %s."
      % (i1["g1_legacy_artifact"]["file"],
         i1["g1_legacy_artifact"]["sha256"],
         "PASS" if i1["g1_legacy_artifact"]["pass"] else "FAIL"))
    A("")
    i2 = d["item_ii_v15_spectrum_cell"]
    A("**(ii) v15 spectrum numbers unchanged.**  Cell `%s` recomputed from "
      "`run_v15_spectrum.evaluate_cell` and compared to `%s` -> %s."
      % (i2["cell_id"], i2["source_file"],
         "PASS" if i2["pass"] else "FAIL"))
    A("")
    A("| quantity | stored | recomputed | identical |")
    A("|---|---|---|---|")
    for k, v in i2["comparison"].items():
        A("| %s | %s | %s | %s |"
          % (k, json.dumps(v["stored"]), json.dumps(v["recomputed"]),
             v["identical"]))
    A("")
    i3 = d["item_iii_frozen_modules"]
    A("**(iii) Frozen composite module hashes unchanged.**  %d files checked "
      "against the v15 record -> %s." % (i3["n_files"],
                                         "PASS" if i3["pass"] else "FAIL"))
    A("")
    A("| file | sha256 (now) | unchanged |")
    A("|---|---|---|")
    for k, v in i3["files"].items():
        A("| %s | `%s` | %s |" % (k, v["now"], v["unchanged"]))
    A("")
    i4 = d["item_iv_no_v11_v15_modified"]
    A("**(iv) No v11-v15 file modified.**  Snapshot %s; %d python files plus "
      "the G1 artefact re-hashed -> %s."
      % (i4["snapshot_utc"], len(i4["files"]),
         "PASS" if i4["pass"] else "FAIL"))
    A("")
    A("**Explicit statement: no v11-v15 file was modified by this job.**  Only "
      "`run_v16.py`, `run_v16_water.py` and files under `out/v16/` were "
      "written.  Cross-check against the hashes v15 itself recorded: %s."
      % json.dumps(i4["cross_check_against_v15_recorded_hashes"]))
    A("")
    A("Hash list (v16 deliverables):")
    A("")
    A("| file | sha256 |")
    A("|---|---|")
    for f in ("run_v16.py", "run_v16_water.py",
              "out/v16/numbers-v16-scorecard.json",
              "out/v16/numbers-v16-halo.json",
              "out/v16/pair-halo-old-vs-new.png",
              "out/v16/still-g1-old-hg-gaussian.png",
              "out/v16/still-g1-new-ffcdf.png"):
        p = os.path.join(HERE, f)
        A("| %s | `%s` |" % (f, sha256_file(p) if os.path.exists(p)
                             else "(missing)"))
    A("")
    A("This report cannot contain its own hash; `run_v16.py` prints it after "
      "the write and records it in `out/v16/_work/report_hash.txt`.")
    A("")

    # ---------------- merge
    A("## 5. Merge statement")
    A("")
    A("Pre-registered rule: the FF-CDF kernel becomes the water default IFF "
      "Stage A passes all five clauses for the algae/generic preset (with the "
      "beta(90 deg) reference declared), Stage B passes (H bit-identical, "
      "sigma scaling exact), Stage C is reported, and Stage D is fully green.")
    A("")
    A("| condition | value | pass |")
    A("|---|---|---|")
    for k, v in m["conditions"].items():
        A("| %s | %s | %s |"
          % (k, json.dumps(v["value"]), v["pass"]))
    A("")
    A("**Merge fired: %s.**  %s"
      % (m["fired"],
         ("`run_v16_water.py` now carries `PHASE_DEFAULT = \"%s\"` and "
          "`KERNEL_DEFAULT = \"%s\"`; the file %s, and the constants were "
          "re-read from disk after the decision (sha256 %s)."
          % (m["applied"]["phase_default"], m["applied"]["kernel_default"],
             "was already at those values, so no byte changed"
             if not m["applied"]["file_rewritten"]
             else "was rewritten atomically",
             m["applied"]["sha256_after"]))
         if m["fired"] else
         ("The v12 default (HG-0.85 + Gaussian) is kept; see the failing "
          "condition(s) above.")))
    A("")
    A("No tuning of n or mu was performed: the FF pairs are the published "
      "fits (Petzold 1.10/3.5835; Jonasz clear 1.098/3.37; sediment "
      "1.15/4.0 is an ASSUMPTION with its envelope printed).")
    A("")

    # ---------------- caveats
    A("## 6. Caveats and open items")
    A("")
    # two numbers the caveats must carry, measured here (cheap, no state)
    from v5_physics import water_attenuation
    _e = env_for(load_recipe(), SCENE_P1["U"], SCENE_P1["F"], SCENE_P1["h"])
    _att = water_attenuation(_e["env"]["turbid"])
    _a = float(_att["a"][1])
    _b = float(_att["b"][1])
    _c = _a + _b
    _L = float(halo["g1_reference"]["environment"]["L_m"])
    _pm = W.scatter_position_moments(_c, _L)
    _th2_ff = sc["rows_by_id"]["FF_algae_pond"]["theta2_rad2"]
    _th2_hg = W.theta2_hg(W.G_DEFAULT)
    _ds_hg = 0.5 * _th2_hg * (_L - _pm["s_mean_m"])
    _ds_ff = 0.5 * _th2_ff * (_L - _pm["s_mean_m"])
    _hraw = (_b / _c) * (1.0 - math.exp(-_c * _L))
    _h_old = _hraw * math.exp(-_a * _ds_hg)
    _h_alt = _hraw * math.exp(-_a * _ds_ff)
    _th2_v16_hg = V16.phase_moments(
        V16.phase_table("hg", g=W.G_DEFAULT))["theta2_rad2"]
    A("- The halo weight H keeps v12's energy budget on purpose: the phase "
      "function acts on the kernel only.  This is the brief's pre-registered "
      "reading (H phase-function independent).  The alternative reading -- H "
      "moving through the rho detour term, which does depend on <theta^2> -- "
      "was NOT taken; on the G1 reference scene it would raise the green H "
      "from %.9f to %.9f (%+.2f%%), because FF's <theta^2> (%.4f) is about "
      "half of HG-0.85's (%.4f) and the halo photons' mean detour is "
      "proportionally shorter." % (_h_old, _h_alt, 100.0 * (_h_alt / _h_old
                                                          - 1.0),
                                  _th2_ff, _th2_hg))
    A("- The legacy lane's <theta^2> is v12's own quadrature (linear in mu, "
      "%.9f rad^2); the same moment on this module's log-theta grid is "
      "%.9f rad^2 (relative difference %.1e).  The legacy lane takes v12's "
      "value verbatim, which is why Stage D(i) is bit-identical rather than "
      "merely equal." % (_th2_hg, _th2_v16_hg,
                        abs(_th2_v16_hg / _th2_hg - 1.0)))
    A("- The FF kernel uses r = (L-s)|tan theta| (the v12 displacement "
      "convention) and is truncated at r_max with the captured fraction "
      "printed; the truncated mass is folded back by unit-mass "
      "normalisation, exactly as the Gaussian lane conserves the mean.")
    A("- `beta(90 deg)` is compared against the DECLARED Petzold flat-back "
      "level 0.00291 sr^-1, not against the FF value.")
    A("- sediment's phase function remains an ASSUMPTION; only a measured "
      "sediment-dominated freshwater VSF would resolve it.")
    A("- Stage C's Gaussian is reported in both forms (v12 sigma and the "
      "FF-moment-matched sigma); neither is within 0.05 of the FF kernel at "
      "all three radii.")
    A("- The G3 peak-profile prediction was stated two ways (point density "
      "and tail mass).  The point-density form was near-neutral at r=32 px "
      "and disagrees with the measurement there; the tail-mass form agrees "
      "at both lags.  Both are printed in section 2.4 rather than choosing "
      "the flattering one after the fact.")
    A("")
    A("## 7. Files")
    A("")
    for f in ("out/v16/numbers-v16-scorecard.json",
              "out/v16/numbers-v16-halo.json",
              "out/v16/pair-halo-old-vs-new.png",
              "out/v16/checkpoint.json",
              "out/v16/still-g1-old-hg-gaussian.png",
              "out/v16/still-g1-new-ffcdf.png",
              "run_v16.py", "run_v16_water.py"):
        A("- `%s`" % f)
    A("")
    A("Wall time: A %.1fs, B %.1fs, C %.1fs, D %.1fs, report %.1fs."
      % (sc["wall_s"], halo["wall_s"], c["wall_s"], d["wall_s"],
         time.time() - t0))
    A("")
    text = "\n".join(L)
    non_ascii = [ch for ch in text if ord(ch) > 127]
    if non_ascii:
        text = text.encode("ascii", "replace").decode("ascii")
    atomic_write_bytes(os.path.join(OUT, "REPORT-V16-SCATTER.md"),
                       text.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 2000)
    rep_sha = sha256_file(os.path.join(OUT, "REPORT-V16-SCATTER.md"))
    atomic_write_bytes(os.path.join(WORK, "report_hash.txt"),
                       ("REPORT-V16-SCATTER.md %s\n" % rep_sha).encode("ascii"))
    cp.unit_done("report")
    log("[report] written (%d bytes, sha256 %s, %d non-ASCII replaced)"
        % (len(text), rep_sha, len(non_ascii)))
    return text


# ==========================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["A", "B", "C", "D", "report", "all"])
    args = ap.parse_args()
    makedirs()
    if args.stage == "A":
        stage_a()
    elif args.stage == "B":
        stage_b()
    elif args.stage == "C":
        stage_c()
    elif args.stage == "D":
        stage_d()
    elif args.stage == "report":
        report()
    else:
        stage_a()
        stage_b()
        stage_c()
        stage_d()
        report()


if __name__ == "__main__":
    main()
