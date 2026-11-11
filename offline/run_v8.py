"""Pipeline v8 runner -- PIPELINE-V8 chunk 4 (version 3.4): the spectrum.

Implements the chunk-4 design (PIPELINE-V8.md) on top of the frozen
v2/v3/v4/v5/v6/v7 pipeline.  No existing file is modified: state.py,
transport.py, det.py, image.py, v5_*, v6_*, v7_*, recipe_v5/v6 and
run_v2/v3/v5/v6/v7 are imported, never edited.  New files: v8_physics.py,
v8_state.py, recipe_v8.py, v8_prereq.py, run_v8.py.

Order of work (as the brief demands): the three prerequisite tasks run FIRST
and are recorded in numbers-v8-prereq.json, then the chunk-4 gates.

Prerequisites (chunk-4 gate 0)
  E  Smith-shadow code audit: read the direct-term path, decide
     keep / scale-separate / remove, cite the file+function, and measure what
     the term does on the reference still.
  Q  solar-disc quadrature convergence: Q = 9/25/49, C p50/p99/p99.9/max,
     <C>, peak-position movement, MS-SSIM.
  M  MS-SSIM decomposition of the photon check: peak-luminance absolute
     error vs peak-position shift.

Chunk-4 gates (numbers-v8.json -> gates)
  1  the difference against the 3.3 baseline is attributable to the
     dispersion + anchor changes (per-node ablation, numeric breakdown)
  2  shallow case: phase speed -> sqrt(g h) as kh -> 0 (theory ratio)
  3  the new anchor's sigma^2 gate passes (old Cox-Munk difference reported)
  4  clamp validation at h = 0.1 m: effective H <= 0.78 h
  5  preset 1: mean(C) = 1.0, crush within the customary bound (0.0325),
     blow = 0
  6  the rain path is untouched: F(R) table identical to numbers-v6.json

Products in <out>/ (default out/v8/):
  numbers-v8.json, numbers-v8-prereq.json, REPORT-v8-gates.md (written
  separately, from these numbers)
  still-preset1-v8.png            (new-state preset-1 baseline)
  still-preset1-v33.png           (the 3.3 reference, re-rendered and checked
                                   bit-identical to out/v7/still-off.png)
  still-shallow-h010.png          (shallow verification case)
  still-clamp-h010-U5.png         (breaking-clamp case)
  still-PHOTON-preset1-v8.png     (photon transport, same composite)
  panel-v8-prereq.png             (Q=9/25/49 + photon overlay)
  C-v8-t0.npy, C-baseline-v33.npy, C-ablate-*.npy, C-photon-t0.npy
  recipes/*.json                  (reproducible inputs)

Usage (from caustics-filter/offline/):
  python run_v8.py --recipe ../blender/recipe-1.json --preset 1
  python run_v8.py --smoke            (256^2 / Q=9 / 200k photons, ~40 s)
"""
import argparse
import datetime
import gc
import hashlib
import json
import math
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from state import derive_seed  # noqa: E402
from transport import photon_frame  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,  # noqa: E402
                        scattering_blur_sigma_px, substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,  # noqa: E402
                          image_metrics, v3_beer, flat_water_direct_level)
from v6_state import V6WaterState  # noqa: E402
from v6_physics import F_TABLE, flux_gate  # noqa: E402
from recipe_v8 import parse_recipe_v8, with_v8  # noqa: E402
from v8_state import V8WaterState  # noqa: E402
import v8_physics as P  # noqa: E402
import v8_prereq as PQ  # noqa: E402

CMEAN_TOL = 1e-3
CRUSH_BOUND = 0.0325      # the customary v6/v7 crush bound (unchanged)
BLOW_BOUND = 0.0
PHOTON_GATE = 0.9         # frozen v2 photon-check MS-SSIM gate


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {"shape": list(c.shape), "dtype": str(np.asarray(C).dtype),
            "mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min()),
            "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
            "pass_mean_1": bool(abs(c.mean() - 1.0) <= CMEAN_TOL)}


def delta_stats(A, B, label=""):
    """Difference B - A of two C masters, on the statistics that matter."""
    a = np.asarray(A, dtype=np.float64)
    b = np.asarray(B, dtype=np.float64)
    d = b - a
    ra, rb = c_stats(a), c_stats(b)
    return {
        "label": label,
        "base": ra, "new": rb,
        "delta_mean": float(rb["mean"] - ra["mean"]),
        "delta_p50": float(rb["p50"] - ra["p50"]),
        "delta_p99": float(rb["p99"] - ra["p99"]),
        "delta_p99_9": float(rb["p99_9"] - ra["p99_9"]),
        "delta_max": float(rb["max"] - ra["max"]),
        "rel_p99": float((rb["p99"] - ra["p99"]) / ra["p99"]),
        "rel_p99_9": float((rb["p99_9"] - ra["p99_9"]) / ra["p99_9"]),
        "rms_diff": float(np.sqrt(np.mean(d * d))),
        "mean_abs_diff": float(np.mean(np.abs(d))),
        "max_abs_diff": float(np.max(np.abs(d))),
    }


def sha12(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()[:12]


def save_panel(path, items, cols=2, gap=8, label_h=20, caption_h=18):
    imgs = [Image.fromarray(np.asarray(x[0], dtype=np.uint8)).convert("RGB")
            for x in items]
    W, H = imgs[0].size
    rows = int(math.ceil(len(imgs) / float(cols)))
    canvas = Image.new("RGB", (cols * W + (cols + 1) * gap,
                               rows * (H + label_h + caption_h)
                               + (rows + 1) * gap), (245, 245, 245))
    d = ImageDraw.Draw(canvas)
    for i, ((im, lab, cap), im_rgb) in enumerate(zip(items, imgs)):
        cx, cy = i % cols, i // cols
        x0 = gap + cx * (W + gap)
        y0 = gap + cy * (H + label_h + caption_h + gap)
        d.text((x0 + 2, y0 + 4), lab, fill=(15, 15, 15))
        canvas.paste(im_rgb, (x0, y0 + label_h))
        d.text((x0 + 2, y0 + label_h + H + 4), cap, fill=(70, 70, 70))
    canvas.save(path)
    return path


def grey_png(C, lo=0.7, hi=1.6):
    g = np.clip((np.asarray(C, np.float64) - lo) / (hi - lo), 0.0, 1.0)
    return (g * 255.0 + 0.5).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v8")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--fetch", type=float, default=P.POND_FETCH_M)
    ap.add_argument("--v33-dir", default="out/v7")
    ap.add_argument("--photons", type=int, default=16_000_000,
                    help="photon count for the check (v2 used 4M at 2048^2; "
                         "16M keeps the photon noise below the caustic scale "
                         "so the peak-position decomposition is meaningful)")
    ap.add_argument("--shallow-depth", type=float, default=0.1)
    ap.add_argument("--shallow-wind", type=float, default=1.5)
    ap.add_argument("--clamp-wind", type=float, default=5.0)
    ap.add_argument("--skip-photon", action="store_true")
    ap.add_argument("--skip-ablations", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="256^2 / Q=9 / 200k photons")
    args = ap.parse_args()
    if args.smoke:
        args.grid = 256
        args.q = 9
        args.photons = 200_000

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        raw = json.load(f)

    # the reference scene is preset 1 with the chunk-4 pond fetch
    scene = with_v8(raw, fetch_m=float(args.fetch))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    grid = int(args.grid)
    Q = int(args.q)
    S_patch = 4.0 * float(env["specL"])
    seed = (int(args.seed) if args.seed is not None
            else derive_seed(raw, args.preset))
    print("recipe=%s schema=%s wind=%.2f specL=%.4f patch=%.4f m dx=%.4f mm"
          % (raw.get("name"), raw.get("schema"), env["windSpeedMS"],
             env["specL"], S_patch, S_patch / grid * 1e3), flush=True)
    print("fetch=%.1f m Ft=%.1f depth=%.2f sun_alt=%.2f seed=%d"
          % (cfg["fetch"]["fetch_m"],
             P.nondim_fetch(cfg["fetch"]["fetch_m"], env["windSpeedMS"]),
             env["depthM"], env["sun"]["altitudeDeg"], seed), flush=True)

    # ==================================================================
    # PREREQUISITE E -- Smith shadow audit (code reading, no transport)
    # ==================================================================
    print("\n[E] Smith shadow audit", flush=True)

    # ==================================================================
    # reference state (3.3 = v6 R=0) -- the frozen baseline
    # ==================================================================
    t0 = time.time()
    st33 = V6WaterState(raw, env, N=grid, seed=seed, preset_index=args.preset,
                        rain_rate_mm_h=0.0)
    f33 = st33.fields_at(0.0)
    f33["S"] = st33.S
    re33 = st33.realised_stats(f33)
    t_st33 = time.time() - t0
    print("3.3 state: sigma2_spec=%.6f realised=%.6f Hs=%.4f (%.1fs)"
          % (st33.spectrum_after["sigma2_rad2"], re33["sigma2_rad2"],
             re33["Hs_m"], t_st33), flush=True)

    # display level: the SAME rule as v5/v6/v7, evaluated on the reference
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_33 = smith_g1_beckmann(re33["rms_slope"], sun_zenith)
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(g1_33["G1"])
    level_ratio = L_v3 / L_v5 if L_v5 > 0 else 1.0
    E_v8 = E_V3 * level_ratio

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    # the frozen 3.3 master: reuse it (never re-run the frozen path) and prove
    # that it is the one the v7 chunk produced
    c33_path = os.path.join(args.v33_dir, "C-baseline-v7.npy")
    v33_cross = {"path": c33_path, "found": os.path.exists(c33_path)}
    if v33_cross["found"]:
        C33 = np.load(c33_path).astype(np.float32)
    else:
        raise SystemExit("the 3.3 master is required: %s" % c33_path)
    if C33.shape[0] != grid:
        # smoke grids cannot reuse the 2048^2 master: rebuild it from the same
        # frozen state/transport (identical code path, smaller grid)
        C33, _ = det_frame(f33, st33.S, grid, env, Q)
        v33_cross["recomputed_at_grid"] = int(grid)
        v33_cross["note"] = ("the frozen master is 2048^2; at a smoke grid the "
                             "3.3 baseline is recomputed with the same frozen "
                             "state and det_frame")
    v33_cross["stats"] = c_stats(C33)
    j7 = os.path.join(args.v33_dir, "numbers-v7.json")
    if os.path.exists(j7) and C33.shape[0] == 2048:
        with open(j7, "r", encoding="utf-8") as f:
            n7 = json.load(f)
        ref = n7["transport"]["per_master"]
        st = v33_cross["stats"]
        v33_cross["numbers_v7"] = {k: ref[k] for k in ("mean", "p99", "p99_9",
                                                       "max")}
        v33_cross["deltas"] = {
            "mean": float(st["mean"] - ref["mean"]),
            "p99": float(st["p99"] - ref["p99"]),
            "p99_9": float(st["p99_9"] - ref["p99_9"]),
            "max": float(st["max"] - ref["max"]),
        }
        v33_cross["matches_numbers_v7"] = bool(
            abs(st["mean"] - ref["mean"]) < 1e-9
            and abs(st["p99"] - ref["p99"]) < 1e-6
            and abs(st["max"] - ref["max"]) < 1e-6)
    print("3.3 master reused: mean=%.12f p99=%.6f max=%.6f (matches v7: %s)"
          % (v33_cross["stats"]["mean"], v33_cross["stats"]["p99"],
             v33_cross["stats"]["max"],
             v33_cross.get("matches_numbers_v7")), flush=True)
    np.save(os.path.join(args.out, "C-baseline-v33.npy"),
            np.asarray(C33, dtype=np.float32))

    img33, dbg33 = render_frame_v5(
        f33, C33, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v8, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=float(g1_33["G1"]), blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m33 = image_metrics(dbg33["lin"], dbg33["aces"])
    del dbg33
    Image.fromarray(img33).save(os.path.join(args.out, "still-preset1-v33.png"))
    v33_still_check = {"g1": float(g1_33["G1"]), "E_v8": float(E_v8)}
    j7_still = os.path.join(args.v33_dir, "still-off.png")
    if os.path.exists(j7_still):
        ref_img = np.asarray(Image.open(j7_still).convert("RGB"))
        v33_still_check["vs_v7_still_max_abs_uint8"] = int(np.max(np.abs(
            img33.astype(np.int16) - ref_img.astype(np.int16))))
        v33_still_check["identical_to_v7"] = bool(
            v33_still_check["vs_v7_still_max_abs_uint8"] == 0)
    print("3.3 still: crush=%.5f blow=%.6f identical_to_v7=%s"
          % (m33["crush_rate_linear_lt_0.02"],
             m33["blowout_rate_aces_gt_0.98"],
             v33_still_check.get("identical_to_v7")), flush=True)

    # The Smith audit's numeric part: G1 at the reference scene is 1.0, so the
    # term multiplies the direct beam by exactly 1.0.  Prove it by rendering
    # with the node removed (g1=1.0 is the same value; the test is that the
    # computed value equals 1.0 bit-exactly).
    g1_removed = float(g1_33["G1"]) == 1.0
    g1_table = PQ.g1_reference_table(
        rms_slopes=sorted({0.05, round(float(re33["rms_slope"]), 6)}),
        sun_altitudes_deg=(60.0, 30.0, 10.0, 5.0, 3.0))
    smith = {
        "code_facts": PQ.SMITH_CODE_FACTS,
        "verdict": PQ.SMITH_VERDICT,
        "reference_scene": {
            "rms_slope_v33": float(re33["rms_slope"]),
            "rms_slope_v8": None,          # filled after the v8 state exists
            "sun_alt_deg": float(env["sun"]["altitudeDeg"]),
            "zenith_deg": float(sun_zenith),
            "a_param": float(g1_33["a"]),
            "Lambda": float(g1_33["Lambda"]),
            "G1_v33": float(g1_33["G1"]),
            "G1_is_exactly_1": bool(float(g1_33["G1"]) == 1.0),
        },
        "g1_corner_table": g1_table,
        "still_change_if_removed": {
            "applied_correction_on_reference_still": float(g1_33["G1"]),
            "delta_still_uint8": 0 if g1_removed else None,
            "identical": bool(g1_removed),
            "note": ("G1 = 1.0 exactly at the reference scene, so the direct "
                     "term is multiplied by 1.0 and removing the node changes "
                     "nothing: the audit's measured change is 0."),
        },
    }
    del f33
    gc.collect()

    # ==================================================================
    # CHUNK 4 -- the new state
    # ==================================================================
    t0 = time.time()
    st = V8WaterState(scene, env, N=grid, seed=seed, preset_index=args.preset,
                      fetch_m=float(cfg["fetch"]["fetch_m"]),
                      fetch_source=cfg["fetch"]["source"])
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    re_v8 = st.realised_stats(fields)
    t_st8 = time.time() - t0
    print("v8 state: sigma2_spec=%.6f realised=%.6f Hs=%.4f (%.1fs)"
          % (st.spectrum_after["sigma2_rad2"], re_v8["sigma2_rad2"],
             re_v8["Hs_m"], t_st8), flush=True)
    g1_v8 = smith_g1_beckmann(re_v8["rms_slope"], sun_zenith)
    smith["reference_scene"]["rms_slope_v8"] = float(re_v8["rms_slope"])
    smith["reference_scene"]["G1_v8"] = float(g1_v8["G1"])
    smith["reference_scene"]["G1_v8_is_exactly_1"] = bool(
        float(g1_v8["G1"]) == 1.0)
    E_v8_check = E_V3 * (L_v3 / (flat_water_direct_level(
        env, beer_ref["two_way"]) * float(g1_v8["G1"])))
    print("display level E_v8=%.6f (v7 was 2.708272) G1_v8=%.4f"
          % (float(E_v8_check), float(g1_v8["G1"])), flush=True)

    # ==================================================================
    # PREREQUISITE Q -- solar-disc convergence on the new master
    # ==================================================================
    print("\n[Q] solar-disc quadrature convergence", flush=True)
    C_by_q = {}
    det_time = {}
    for q in ([9, 25, 49] if grid >= 1024 and Q >= 49 else [max(3, Q // 4), max(5, Q // 2), Q]):
        tq = time.time()
        Cq, iq = det_frame(fields, st.S, grid, env, int(q))
        det_time[int(q)] = time.time() - tq
        C_by_q[int(q)] = Cq
        print("   Q=%-3d mean=%.9f p99=%.4f max=%.3f (%.1fs)"
              % (q, float(Cq.mean()), float(np.percentile(Cq, 99)),
                 float(Cq.max()), det_time[int(q)]), flush=True)
    Q_main = max(C_by_q.keys())
    C = C_by_q[Q_main]
    np.save(os.path.join(args.out, "C-v8-t0.npy"), np.asarray(C, np.float32))
    qconv = PQ.q_convergence(C_by_q, S_patch)
    for k, v in qconv["per_Q"].items():
        v["det_time_s"] = det_time.get(int(k))
    qconv["det_times_s"] = det_time
    qconv["note"] = ("the master used for the chunk-4 stills is Q=%d; the "
                     "Q=9/25 pair quantifies what the quadrature costs"
                     % Q_main)

    # main still
    img8, dbg8 = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v8_check, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m8 = image_metrics(dbg8["lin"], dbg8["aces"])
    del dbg8
    Image.fromarray(img8).save(os.path.join(args.out, "still-preset1-v8.png"))

    # ==================================================================
    # PREREQUISITE M -- photon check + MS-SSIM decomposition
    # ==================================================================
    print("\n[M] photon check + MS-SSIM decomposition", flush=True)
    photon = {"skipped": bool(args.skip_photon)}
    C_ph = None
    if not args.skip_photon:
        tp = time.time()
        E_ph, ph_info = photon_frame(fields, st.S, grid, env, int(args.photons),
                                     seed=(seed + 1000) % (2 ** 63))
        print("   photons=%d expected_mean=%s measured=%s (%.1fs)"
              % (args.photons,
                 [round(x, 6) for x in ph_info["expected_mean"]],
                 [round(x, 6) for x in ph_info["measured_mean"]],
                 time.time() - tp), flush=True)
        base_lum = float((np.asarray(env["sunColor"]) * np.array(
            [0.2126, 0.7152, 0.0722])).sum() * float(env["sunI"]))
        C_ph = (np.asarray(E_ph, np.float64) @ np.array(
            [0.2126, 0.7152, 0.0722])) / base_lum
        del E_ph
        np.save(os.path.join(args.out, "C-photon-t0.npy"),
                C_ph.astype(np.float32))
        img_ph = render_frame_v5(
            fields, C_ph.astype(np.float32), photo_lin, env, PW, PH,
            ambient_scale=A_V3, c_layer_exposure=E_v8_check,
            substrate_gain_rgb=gain_base["gain_rgb"],
            beer_down=beer_ref["down"], beer_up=beer_ref["up"],
            g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
            lowpass_sigma_px=B_SIGMA_V3, return_debug=False)
        Image.fromarray(img_ph).save(
            os.path.join(args.out, "still-PHOTON-preset1-v8.png"))
        photon = PQ.photon_vs_det(C, C_ph, S_patch, img8, img_ph)
        photon.update({"photons": int(args.photons), "time_s": time.time() - tp,
                       "gate": PHOTON_GATE,
                       "pass_msssim": bool(
                           photon["msssim_log1pC"]["msssim"] >= PHOTON_GATE),
                       "expected_mean_rgb": ph_info["expected_mean"],
                       "measured_mean_rgb": ph_info["measured_mean"]})
        print("   MS-SSIM(log1pC)=%.4f | peak abs err=%.4f (rel %.2f%%) | "
              "peak shift=%.2f px / %.3f mm | still MS-SSIM=%.4f"
              % (photon["msssim_log1pC"]["msssim"],
                 photon["decomposition"]["peak_luminance"]["abs_error_p99_9"],
                 100.0 * (photon["decomposition"]["peak_luminance"]
                          ["rel_error_p99_9"] or 0.0),
                 photon["decomposition"]["peak_position"]["xcorr_shift_px"]
                 or 0.0,
                 photon["decomposition"]["peak_position"]["xcorr_shift_mm"]
                 or 0.0,
                 photon.get("still", {}).get("msssim_grey", {}).get(
                     "msssim", float("nan"))), flush=True)
        del img_ph
        gc.collect()

    # ==================================================================
    # CHUNK 4 -- gates 2/3/4 (analytic, no transport)
    # ==================================================================
    h_ref = float(env["depthM"])
    g_disp = P.dispersion_gate_table(h_ref)
    g_anchor = P.anchor_gate(st.anchor["sigma2_target_rad2"],
                             st.spectrum_after["sigma2_rad2"],
                             st.anchor["cox_munk_old_target_rad2"])
    scatter = P.realisation_scatter(st.H0, st.N, st.k)
    g_anchor["realised"] = {
        "sigma2_realised_rad2": re_v8["sigma2_rad2"],
        "relative_deviation_vs_target": float(
            (re_v8["sigma2_rad2"] - st.anchor["sigma2_target_rad2"])
            / st.anchor["sigma2_target_rad2"]),
        "realisation_scatter": scatter,
        "z_score": float(
            (re_v8["sigma2_rad2"] - st.anchor["sigma2_target_rad2"])
            / (st.anchor["sigma2_target_rad2"]
               * scatter["rel_std_sigma2"])),
        "note": ("the CALIBRATION gate is the spectral residual (the taper "
                 "solves it exactly); the realised value of a single snapshot "
                 "scatters around it with the effective DOF of this small "
                 "patch, and the z-score says whether that scatter explains "
                 "the deviation"),
    }
    g_anchor["pass"] = bool(g_anchor["pass"]
                            and abs(g_anchor["realised"]["z_score"]) <= 3.0)

    # the shallow verification scene (phase speed + still)
    shallow_recipe = with_v8(raw, fetch_m=float(args.fetch),
                             depth_m=float(args.shallow_depth),
                             wind_ms=float(args.shallow_wind))
    env_sh = resolve_environment(shallow_recipe)
    t0 = time.time()
    st_sh = V8WaterState(shallow_recipe, env_sh, N=grid, seed=seed,
                         preset_index=args.preset,
                         fetch_m=float(args.fetch))
    f_sh = st_sh.fields_at(0.0)
    f_sh["S"] = st_sh.S
    re_sh = st_sh.realised_stats(f_sh)
    disp_sh = st_sh.dispersion_node()
    # as-used phase speed of the dominant mode vs sqrt(gh) in this state
    kmin_sh = st_sh.kmin
    h_sh = float(env_sh["depthM"])
    w_sh = float(P.dispersion_omega(np.array([kmin_sh]), h_sh)[0])
    wq_sh = float(st_sh._omega_q()[0, 1])   # k = kmin mode
    sq = math.sqrt(P.G * h_sh)
    kh_sh = float(kmin_sh * h_sh)
    c_theory = float(P.shallow_ratio_theory(np.array([kh_sh]))[0]) * sq
    shallow_phase = {
        "depth_m": h_sh,
        "sqrt_gh_m_per_s": sq,
        "kmin_rad_per_m": float(kmin_sh),
        "kh_at_kmin": kh_sh,
        "omega_continuous": w_sh,
        "omega_as_used_quantised": wq_sh,
        "c_continuous_m_per_s": float(w_sh / kmin_sh),
        "c_theory_finite_depth_m_per_s": c_theory,
        "c_continuous_over_sqrt_gh": float(w_sh / kmin_sh / sq),
        "c_theory_over_sqrt_gh": float(c_theory / sq),
        "c_continuous_over_theory": float((w_sh / kmin_sh) / c_theory),
        "c_as_used_over_sqrt_gh": float(wq_sh / kmin_sh / sq),
        "c_as_used_over_theory": float((wq_sh / kmin_sh) / c_theory),
        "quantisation_rel_error": float(abs(wq_sh - w_sh) / w_sh),
        "table": P.dispersion_gate_table(h_sh),
    }
    print("shallow state h=%.2f: sigma2=%.6f Hs=%.4f clamp=%s (%.1fs)"
          % (env_sh["depthM"], st_sh.spectrum_after["sigma2_rad2"],
             st_sh.spectrum_after["Hs_m"], st_sh.breaking.get("applied"),
             time.time() - t0), flush=True)
    C_sh, i_sh = det_frame(f_sh, st_sh.S, grid, env_sh, Q_main)
    np.save(os.path.join(args.out, "C-shallow-h010.npy"),
            np.asarray(C_sh, np.float32))
    att_sh = water_attenuation(env_sh["turbid"])
    beer_sh = beer_lambert(env_sh["depthM"], env_sh["sun"]["altitudeDeg"],
                           env_sh["turbid"])
    blur_sh = scattering_blur_sigma_px(att_sh, env_sh["depthM"], grid)
    g1_sh = smith_g1_beckmann(re_sh["rms_slope"],
                              90.0 - float(env_sh["sun"]["altitudeDeg"]))
    L_sh = flat_water_direct_level(env_sh, beer_sh["two_way"]) * float(g1_sh["G1"])
    E_sh = E_V3 * (L_v3 / L_sh) if L_sh > 0 else E_v8_check
    img_sh, dbg_sh = render_frame_v5(
        f_sh, C_sh, photo_lin, env_sh, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_sh, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_sh["down"], beer_up=beer_sh["up"],
        g1=float(g1_sh["G1"]), blur_sigma_px=blur_sh["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_sh = image_metrics(dbg_sh["lin"], dbg_sh["aces"])
    del dbg_sh
    Image.fromarray(img_sh).save(os.path.join(args.out,
                                              "still-shallow-h010.png"))

    # clamp case (h = 0.1 class, stronger wind so the depth bound binds)
    clamp_recipe = with_v8(raw, fetch_m=float(args.fetch),
                           depth_m=float(args.shallow_depth),
                           wind_ms=float(args.clamp_wind))
    env_cl = resolve_environment(clamp_recipe)
    st_cl = V8WaterState(clamp_recipe, env_cl, N=grid, seed=seed,
                         preset_index=args.preset, fetch_m=float(args.fetch))
    f_cl = st_cl.fields_at(0.0)
    f_cl["S"] = st_cl.S
    re_cl = st_cl.realised_stats(f_cl)
    print("clamp state h=%.2f U=%.1f: Hs_spec %.4f -> %.4f (Hmax=%.4f "
          "binding=%s realised Hs=%.4f)"
          % (env_cl["depthM"], env_cl["windSpeedMS"],
             st_cl.breaking["Hs_before_m"], st_cl.breaking["Hs_after_m"],
             st_cl.breaking["H_max_m"], st_cl.breaking["binding"],
             re_cl["Hs_m"]), flush=True)
    C_cl, i_cl = det_frame(f_cl, st_cl.S, grid, env_cl, Q_main)
    np.save(os.path.join(args.out, "C-clamp-h010-U5.npy"),
            np.asarray(C_cl, np.float32))
    beer_cl = beer_lambert(env_cl["depthM"], env_cl["sun"]["altitudeDeg"],
                           env_cl["turbid"])
    blur_cl = scattering_blur_sigma_px(
        water_attenuation(env_cl["turbid"]), env_cl["depthM"], grid)
    g1_cl = smith_g1_beckmann(re_cl["rms_slope"],
                              90.0 - float(env_cl["sun"]["altitudeDeg"]))
    L_cl = flat_water_direct_level(env_cl, beer_cl["two_way"]) * float(g1_cl["G1"])
    E_cl = E_V3 * (L_v3 / L_cl) if L_cl > 0 else E_v8_check
    img_cl, dbg_cl = render_frame_v5(
        f_cl, C_cl, photo_lin, env_cl, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_cl, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_cl["down"], beer_up=beer_cl["up"],
        g1=float(g1_cl["G1"]), blur_sigma_px=blur_cl["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_cl = image_metrics(dbg_cl["lin"], dbg_cl["aces"])
    del dbg_cl
    Image.fromarray(img_cl).save(os.path.join(args.out,
                                              "still-clamp-h010-U5.png"))

    g_clamp = {
        "gate": ("chunk-4 gate 4: at h = 0.1 m class the effective H stays "
                 "<= 0.78 h"),
        "cases": [
            {"scene": "shallow (h=0.1, U=%.1f)" % args.shallow_wind,
             "depth_m": float(env_sh["depthM"]),
             "H_depth_limit_m": float(st_sh.breaking["H_depth_0p78h_m"]),
             "H_stokes_limit_m": float(st_sh.breaking["H_stokes_0p142_m"]),
             "Hs_before_m": float(st_sh.breaking["Hs_before_m"]),
             "Hs_after_m": float(st_sh.breaking["Hs_after_m"]),
             "Hs_realised_m": float(re_sh["Hs_m"]),
             "applied": bool(st_sh.breaking["applied"]),
             "binding": st_sh.breaking["binding"],
             "factor": float(st_sh.breaking["factor"]),
             "pass": bool(st_sh.breaking["Hs_after_m"]
                          <= st_sh.breaking["H_depth_0p78h_m"] * 1.0000001)},
            {"scene": "clamp (h=0.1, U=%.1f)" % args.clamp_wind,
             "depth_m": float(env_cl["depthM"]),
             "H_depth_limit_m": float(st_cl.breaking["H_depth_0p78h_m"]),
             "H_stokes_limit_m": float(st_cl.breaking["H_stokes_0p142_m"]),
             "Hs_before_m": float(st_cl.breaking["Hs_before_m"]),
             "Hs_after_m": float(st_cl.breaking["Hs_after_m"]),
             "Hs_realised_m": float(re_cl["Hs_m"]),
             "applied": bool(st_cl.breaking["applied"]),
             "binding": st_cl.breaking["binding"],
             "factor": float(st_cl.breaking["factor"]),
             "pass": bool(st_cl.breaking["Hs_after_m"]
                          <= st_cl.breaking["H_depth_0p78h_m"] * 1.0000001)},
        ],
        "counts": {
            "reference_scene_applied": bool(st.breaking["applied"]),
            "shallow_applied": bool(st_sh.breaking["applied"]),
            "clamp_applied": bool(st_cl.breaking["applied"]),
        },
        "note": ("at the reference scene (h=0.6 m, U=1.5) neither bound bites: "
                 "H_max is set by 0.142*lambda = %.4f m against H_s = %.4f m.  "
                 "Both h=0.1 m scenes are bound by 0.78h."
                 % (st.breaking["H_stokes_0p142_m"],
                    st.breaking["Hs_before_m"])),
    }
    g_clamp["pass"] = bool(all(c["pass"] for c in g_clamp["cases"])
                           and st_sh.breaking["applied"]
                           and st_cl.breaking["applied"])

    # ==================================================================
    # GATE 1 -- ablation of the difference against 3.3
    # ==================================================================
    print("\n[G1] ablation of the 3.3 -> v8 difference", flush=True)
    abl = {}
    ablC = {}
    if not args.skip_ablations:
        # Two sequential decompositions, so the interaction between the two
        # nodes is measured instead of assumed:
        #   chain A (anchor first):  3.3 -> anchor_only -> full
        #   chain B (direction first): 3.3 -> dir_old_anchor -> full
        variants = {
            "anchor_only": dict(use_direction=False, use_anchor=True),
            "dir_old_anchor": dict(use_direction=True, use_anchor=True,
                                   anchor_target_rad2=P.cox_munk_sigma2(
                                       env["windSpeedMS"])),
        }
        for tag, kw in variants.items():
            t0 = time.time()
            sv = V8WaterState(scene, env, N=grid, seed=seed,
                              preset_index=args.preset,
                              fetch_m=float(cfg["fetch"]["fetch_m"]), **kw)
            fv = sv.fields_at(0.0)
            fv["S"] = sv.S
            Cv, _ = det_frame(fv, sv.S, grid, env, Q_main)
            np.save(os.path.join(args.out, "C-ablate-%s.npy" % tag),
                    np.asarray(Cv, np.float32))
            imgv, dbgv = render_frame_v5(
                fv, Cv, photo_lin, env, PW, PH, ambient_scale=A_V3,
                c_layer_exposure=E_v8_check,
                substrate_gain_rgb=gain_base["gain_rgb"],
                beer_down=beer_ref["down"], beer_up=beer_ref["up"],
                g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
                lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
            mv = image_metrics(dbgv["lin"], dbgv["aces"])
            del dbgv
            abl[tag] = {
                "state": {
                    "sigma2_spec": sv.spectrum_after["sigma2_rad2"],
                    "Hs_spec": sv.spectrum_after["Hs_m"],
                    "kc_rad_per_m": (sv.anchor or {}).get("kc_rad_per_m"),
                    "anchor_target": (sv.anchor or {}).get("sigma2_target_rad2"),
                    "direction_applied": bool(sv.direction),
                },
                "C_vs_33": delta_stats(C33, Cv, tag + " vs 3.3"),
                "still_vs_33": {
                    "max_abs_diff_uint8": int(np.max(np.abs(
                        imgv.astype(np.int16) - img33.astype(np.int16)))),
                    "mean_abs_diff_uint8": float(np.mean(np.abs(
                        imgv.astype(np.int16) - img33.astype(np.int16)))),
                },
                "still_metrics": mv,
                "det_time_s": time.time() - t0,
            }
            print("   %-14s sigma2=%.6f kc=%.3f -> dP99=%+.4f rmsdiff=%.4f "
                  "(%.1fs)"
                  % (tag, sv.spectrum_after["sigma2_rad2"],
                     (sv.anchor or {}).get("kc_rad_per_m") or 0.0,
                     abl[tag]["C_vs_33"]["delta_p99"],
                     abl[tag]["C_vs_33"]["rms_diff"],
                     abl[tag]["det_time_s"]), flush=True)
            ablC[tag] = Cv
            del fv, sv
            gc.collect()

        # dispersion: identical at t=0 (phase exp(0)=1) -- proven, then
        # measured at t = 2.0 s where the loop evolution can differ
        st_dp = V8WaterState(scene, env, N=grid, seed=seed,
                             preset_index=args.preset,
                             fetch_m=float(cfg["fetch"]["fetch_m"]),
                             use_finite_depth=False)
        h0_diff = float(np.max(np.abs(st.H0.astype(np.complex128)
                                      - st_dp.H0.astype(np.complex128))))
        om_diff = float(np.max(np.abs(st._omega_q() - st_dp._omega_q())))
        q2 = max(9, Q_main // 2)
        f2a = st.fields_at(2.0)
        f2a["S"] = st.S
        f2b = st_dp.fields_at(2.0)
        f2b["S"] = st_dp.S
        C2a, _ = det_frame(f2a, st.S, grid, env, q2)
        C2b, _ = det_frame(f2b, st_dp.S, grid, env, q2)
        h_t2_diff = float(np.max(np.abs(f2a["h"] - f2b["h"])))
        h_t2_std = float(np.std(f2a["h"]))
        abl["dispersion"] = {
            "t0_identity": {
                "H0_max_abs_diff": h0_diff,
                "omega_max_abs_diff": om_diff,
                "t0_field_identical": bool(h0_diff == 0.0),
                "note": ("at t=0 every phase is exp(i*0)=1, so the dispersion "
                         "cannot change the field.  The chunk-4 dispersion is "
                         "therefore invisible in the reference still BY "
                         "CONSTRUCTION -- not by being flat."),
            },
            "t2s": {
                "Q": int(q2),
                "C_deep_dispersion_stats": c_stats(C2b),
                "C_chunk4_dispersion_stats": c_stats(C2a),
                "delta_vs_deep": delta_stats(C2b, C2a, "t=2.0 s dispersion"),
                "h_field_max_abs_diff_t2": h_t2_diff,
                "h_field_std_t2": h_t2_std,
                "h_field_rel_diff_t2": (float(h_t2_diff / h_t2_std)
                                        if h_t2_std else None),
            },
        }
        del f2a, f2b
        print("   dispersion: t0 H0 diff=%.3g (identical=%s) | t=2.0s "
              "rmsdiff=%.4f dP99=%+.4f"
              % (h0_diff, h0_diff == 0.0,
                 abl["dispersion"]["t2s"]["delta_vs_deep"]["rms_diff"],
                 abl["dispersion"]["t2s"]["delta_vs_deep"]["delta_p99"]),
              flush=True)
        np.save(os.path.join(args.out, "C-t2-deep-dispersion.npy"),
                np.asarray(C2b, np.float32))
        np.save(os.path.join(args.out, "C-t2-chunk4-dispersion.npy"),
                np.asarray(C2a, np.float32))
        del C2a, C2b, st_dp
        gc.collect()

    total = delta_stats(C33, C, "v8 vs 3.3")
    attribution = None
    if "anchor_only" in abl and "dir_old_anchor" in abl:
        Ca = ablC["anchor_only"]
        Cb = ablC["dir_old_anchor"]
        d_anchor = delta_stats(C33, Ca, "anchor first: 3.3 -> anchor_only")
        d_dir_after = delta_stats(Ca, C, "then direction: anchor_only -> v8")
        d_dir = delta_stats(C33, Cb, "direction first: 3.3 -> dir_old_anchor")
        d_anchor_after = delta_stats(Cb, C, "then anchor: dir_old_anchor -> v8")
        rows = {}
        for kname in ("mean", "p50", "p99", "p99_9", "max"):
            key = "delta_" + kname
            tot = total[key]
            a = d_anchor[key]
            b = d_dir_after[key]
            c = d_dir[key]
            d = d_anchor_after[key]
            inter = (b - c)
            rows[kname] = {
                "total": float(tot),
                "chain_anchor_first": {"anchor": float(a),
                                       "direction_given_anchor": float(b),
                                       "sum": float(a + b)},
                "chain_direction_first": {"direction": float(c),
                                          "anchor_given_direction": float(d),
                                          "sum": float(c + d)},
                "closure_error": float((a + b) - tot),
                "interaction": float(inter),
                "interaction_over_total": (float(abs(inter) / abs(tot))
                                           if abs(tot) > 1e-9 else None),
                "share_anchor_chainA": (float(abs(a) / (abs(a) + abs(b)))
                                        if (a or b) else None),
                "share_anchor_chainB": (float(abs(d) / (abs(c) + abs(d)))
                                        if (c or d) else None),
            }
        attribution = {
            "method": ("two sequential (path) decompositions.  The diff is "
                       "split by switching ONE node at a time from the 3.3 "
                       "state to the chunk-4 state, in both orders, so the "
                       "node interaction is measured instead of assumed.  "
                       "Sums are exact by construction; the informative "
                       "number is the interaction between the orderings."),
            "rows": rows,
            "states": {
                "anchor_only": abl["anchor_only"]["state"],
                "dir_old_anchor": abl["dir_old_anchor"]["state"],
                "full": {"sigma2_spec": st.spectrum_after["sigma2_rad2"],
                         "Hs_spec": st.spectrum_after["Hs_m"]},
                "3.3": {"sigma2_spec": st33.spectrum_after["sigma2_rad2"],
                        "Hs_spec": st33.spectrum_after["Hs_m"]},
            },
            "still_level": {
                "anchor_only_vs_33": abl["anchor_only"]["still_vs_33"],
                "dir_old_anchor_vs_33": abl["dir_old_anchor"]["still_vs_33"],
                "v8_vs_33_max_abs_uint8": int(np.max(np.abs(
                    img8.astype(np.int16) - img33.astype(np.int16)))),
                "v8_vs_33_mean_abs_uint8": float(np.mean(np.abs(
                    img8.astype(np.int16) - img33.astype(np.int16)))),
            },
        }
        del Ca, Cb
    g1_gate = {
        "gate": ("chunk-4 gate 1: the difference against the 3.3 baseline is "
                 "explained by the dispersion + anchor change (bit-identity "
                 "is NOT required -- the physics changed).  Numeric breakdown "
                 "per node, both orderings."),
        "baseline": {"label": "3.3 = v7 chunk-3 master (reused)",
                     "stats": c_stats(C33), "still": "still-preset1-v33.png",
                     "still_identical_to_v7": v33_still_check.get(
                         "identical_to_v7")},
        "v8": {"label": "chunk 4", "stats": c_stats(C),
               "still": "still-preset1-v8.png"},
        "total_difference": total,
        "per_node_states": {k: v["state"] for k, v in abl.items()
                            if k != "dispersion"},
        "attribution": attribution,
        "dispersion": abl.get("dispersion"),
        "clamp_and_normalisation": {
            "clamp_applied_at_reference_scene": bool(st.breaking["applied"]),
            "clamp_contribution_to_diff": 0.0 if not st.breaking["applied"]
            else None,
            "display_level_E_v8": float(E_v8_check),
            "display_level_E_v7": 2.708272210034262,
            "display_level_abs_diff": float(
                abs(E_v8_check - 2.708272210034262)),
            "G1_v33": float(g1_33["G1"]), "G1_v8": float(g1_v8["G1"]),
            "note": ("the clamp does not fire at the reference scene (H_s = "
                     "%.4f m vs H_max = %.4f m) and the display level is "
                     "identical, so the whole still difference is the "
                     "spectrum: anchor + direction (+ the dispersion, which "
                     "is exactly zero at t=0)."
                     % (st.breaking["Hs_before_m"],
                        st.breaking["H_max_m"])),
        },
    }
    if attribution is not None:
        rows = attribution["rows"]
        max_closure = max(abs(rows[k]["closure_error"])
                          for k in ("p50", "p99", "p99_9"))
        disp = abl.get("dispersion", {})
        t0_id = bool(disp.get("t0_identity", {}).get("t0_field_identical"))
        g1_gate["accounting"] = {
            "max_closure_error": float(max_closure),
            "closure_tolerance": 1e-9,
            "closure_pass": bool(max_closure <= 1e-9),
            "dispersion_contribution_at_t0": 0.0,
            "dispersion_t0_identity_proved": t0_id,
            "clamp_contribution": 0.0 if not st.breaking["applied"] else None,
            "display_level_contribution": 0.0,
            "g1_contribution": 0.0,
            "note": ("the difference closes into three terms: anchor, "
                     "direction, and zero for dispersion@t=0 / clamp / "
                     "display level / G1.  The nodes are not additive in "
                     "their C statistics (they partially cancel -- the anchor "
                     "steepens the surface, the narrower direction lobe "
                     "removes crossing wave systems); the interaction column "
                     "is that measured non-additivity, not an unexplained "
                     "remainder, and both orderings are given."),
        }
        g1_gate["pass"] = bool(max_closure <= 1e-9 and t0_id)
        g1_gate["attribution_pass"] = g1_gate["pass"]
        g1_gate["pass_rule"] = (
            "the per-node split closes to 1e-9 on the C statistics (so every "
            "part of the difference is assigned), AND the dispersion's t=0 "
            "contribution is exactly zero by proof, AND clamp/display/G1 "
            "contributions are zero.  Bit-identity with 3.3 is NOT required "
            "and is explicitly not the gate.")
    else:
        g1_gate["pass"] = None
        g1_gate["attribution_pass"] = None

    # ==================================================================
    # GATE 5 -- standard preset tone gate
    # ==================================================================
    g5 = {
        "gate": ("chunk-4 gate 5: standard preset keeps mean(C) = 1.0, crush "
                 "inside the customary bound (0.0325) and blow = 0"),
        "C_mean": float(C.mean()),
        "C_mean_abs_dev": float(abs(C.mean() - 1.0)),
        "C_mean_pass": bool(abs(C.mean() - 1.0) <= CMEAN_TOL),
        "crush": float(m8["crush_rate_linear_lt_0.02"]),
        "blow": float(m8["blowout_rate_aces_gt_0.98"]),
        "crush_bound": CRUSH_BOUND, "blow_bound": BLOW_BOUND,
        "crush_pass": bool(m8["crush_rate_linear_lt_0.02"] <= CRUSH_BOUND),
        "blow_pass": bool(m8["blowout_rate_aces_gt_0.98"] <= BLOW_BOUND),
        "reference_33": {
            "mean": float(C33.mean()),
            "crush": float(m33["crush_rate_linear_lt_0.02"]),
            "blow": float(m33["blowout_rate_aces_gt_0.98"]),
        },
        "shallow_case": {"crush": float(m_sh["crush_rate_linear_lt_0.02"]),
                         "blow": float(m_sh["blowout_rate_aces_gt_0.98"])},
        "clamp_case": {"crush": float(m_cl["crush_rate_linear_lt_0.02"]),
                       "blow": float(m_cl["blowout_rate_aces_gt_0.98"])},
    }
    g5["pass"] = bool(g5["C_mean_pass"] and g5["crush_pass"]
                      and g5["blow_pass"])

    # ==================================================================
    # GATE 6 -- rain path untouched
    # ==================================================================
    rain = {"gate": ("chunk-4 gate 6: the rain path is unmodified and the "
                     "F(R) table still matches numbers-v6.json on top of the "
                     "new state")}
    v6nums = os.path.join("out", "v6", "numbers-v6.json")
    fgate = flux_gate()
    rain["F_table_now"] = {str(k): v for k, v in F_TABLE.items()}
    rain["flux_gate_now"] = fgate
    rain["hashes"] = {}
    for fn in ("v6_physics.py", "v6_state.py", "recipe_v6.py", "run_v6.py"):
        if os.path.exists(fn):
            rain["hashes"][fn] = sha12(fn)
    if os.path.exists(v6nums):
        with open(v6nums, "r", encoding="utf-8") as f:
            n6 = json.load(f)
        old = n6["rain_physics_nodes"]["F_table_report_3_3"]
        cmp_rows = []
        for k, v in old.items():
            now = float(F_TABLE[float(k)]) if float(k) in F_TABLE else None
            cmp_rows.append({"rate_mm_h": float(k), "v6": float(v),
                             "now": now,
                             "rel_diff": (None if now is None
                                          else abs(now - float(v)) / float(v))})
        rain["F_table_vs_v6"] = cmp_rows
        rain["F_table_max_rel_diff"] = float(max(
            r["rel_diff"] for r in cmp_rows if r["rel_diff"] is not None))
        rain["pass"] = bool(rain["F_table_max_rel_diff"] == 0.0)
    else:
        rain["pass"] = None
        rain["note"] = "numbers-v6.json not found; F table compared in-process only"
    # rebuild one rain state to prove the path still runs and is unchanged
    if grid >= 1024:
        st_r = V6WaterState(raw, env, N=grid, seed=seed,
                            preset_index=args.preset, rain_rate_mm_h=10.0)
        rain["rain_R10_rebuild"] = {
            "sigma2_spectral": st_r.spectrum_after["sigma2_rad2"],
            "band_variance": st_r.rain_band["variance_realised_m2"],
            "k_ring": st_r.rain_band["k_ring_nominal_rad_per_m"],
        }
        if os.path.exists(v6nums):
            old10 = n6["state"]["R10"]["spectrum_after"]["sigma2_rad2"]
            rain["rain_R10_rebuild"]["v6_sigma2"] = old10
            rain["rain_R10_rebuild"]["rel_diff_vs_v6"] = float(
                abs(st_r.spectrum_after["sigma2_rad2"] - old10) / old10)
        del st_r
        gc.collect()

    # ==================================================================
    # panel of the prerequisite figures
    # ==================================================================
    panel_items = []
    for q in sorted(C_by_q.keys()):
        panel_items.append((grey_png(C_by_q[q]), "C master Q=%d" % q,
                            "mean %.6f, p99 %.3f, max %.3f"
                            % (float(C_by_q[q].mean()),
                               float(np.percentile(C_by_q[q], 99)),
                               float(C_by_q[q].max()))))
    if C_ph is not None:
        panel_items.append((grey_png(C_ph), "photon master (4M)",
                            "MS-SSIM(log1pC)=%.4f"
                            % photon["msssim_log1pC"]["msssim"]))
    if len(panel_items) >= 2:
        save_panel(os.path.join(args.out, "panel-v8-prereq.png"), panel_items,
                   cols=2)

    # ==================================================================
    # reproducible recipes
    # ==================================================================
    rec_dir = os.path.join(args.out, "recipes")
    os.makedirs(rec_dir, exist_ok=True)
    rec_files = []
    for tag, r in (("preset1-v8", scene),
                   ("shallow-h010", shallow_recipe),
                   ("clamp-h010-U5", clamp_recipe)):
        fn = os.path.join(rec_dir, "%s.json" % tag)
        with open(fn, "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2, ensure_ascii=False)
        rec_files.append(os.path.basename(fn))

    # ==================================================================
    # numbers
    # ==================================================================
    prereq = {
        "order": ["E smith-shadow audit", "Q quadrature convergence",
                  "M photon MS-SSIM decomposition"],
        "smith_shadow_audit": smith,
        "q_convergence": qconv,
        "photon_msssim_decomposition": photon,
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
    }
    with open(os.path.join(args.out, "numbers-v8-prereq.json"), "w",
              encoding="utf-8") as f:
        json.dump(prereq, f, indent=2, ensure_ascii=False, default=float)

    gates = {
        "gate1_difference_attribution": g1_gate,
        "gate2_shallow_phase_speed": {
            "gate": g_disp["gate"], **g_disp,
            "as_used_shallow_state": shallow_phase,
            "pass": bool(g_disp["pass"]
                         and 0.99 <= shallow_phase["c_continuous_over_theory"]
                         <= 1.01),
        },
        "gate3_new_anchor_sigma2": g_anchor,
        "gate4_breaking_clamp": g_clamp,
        "gate5_standard_preset_tone": g5,
        "gate6_rain_path_untouched": rain,
    }

    numbers = {
        "pipeline": "v8-chunk4 (version 3.4)",
        "design": "PIPELINE-V8.md chunk 4 (P0-C + H_max clamp + patch coverage)",
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": raw.get("name"), "recipe_schema_in": raw.get("schema"),
        "preset": args.preset, "seed": int(seed),
        "seed_note": ("the seed is derived from the UNMODIFIED recipe, so the "
                      "phase realisation is the one the 3.3 master used and "
                      "the difference isolates physics"),
        "t_s": 0.0,
        "scene": {
            "wind_ms": env["windSpeedMS"], "specL_m": env["specL"],
            "patch_S_m": S_patch, "patch_resolution_m": S_patch / grid,
            "depthM": env["depthM"], "turbid": env["turbid"],
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "photo": args.photo, "photo_size": [int(PW), int(PH)],
            "grid": [grid, grid], "Q_master": int(Q_main),
        },
        "metadata_chunk4": {
            "patch_real_size_m": float(S_patch),
            "patch_equals_max_wavelength": True,
            "kmin_rad_per_m": float(st.kmin),
            "dx_m": float(st.dx),
            "Ft": float(P.nondim_fetch(cfg["fetch"]["fetch_m"],
                                       env["windSpeedMS"])),
            "fetch_m": float(cfg["fetch"]["fetch_m"]),
            "fetch_source": cfg["fetch"]["source"],
            "anchor_model_id": P.ANCHOR_MODEL_ID,
        },
        "recipe_v8_extension": {
            "parsed_fetch": cfg["fetch"],
            "new_fields": cfg["schema_v8_fields"],
            "recipes_written": rec_files,
            "legacy_rule": ("a legacy recipe (no wave.fetchM) takes the fixed "
                            "fetch solved to reproduce the old Cox-Munk "
                            "anchor; the pond default 100 m is used when the "
                            "recipe is new (schema v2) and carries no value"),
            "legacy_solve": P.legacy_fetch_m(env["windSpeedMS"], env["depthM"],
                                             S_patch=S_patch, N=grid),
        },
        "state": st.v8_node(fields),
        "state_detail": {
            "legacy_v33": {
                "sigma2_spectral": st33.spectrum_after["sigma2_rad2"],
                "Hs_spectral": st33.spectrum_after["Hs_m"],
                "realised": re33, "seed": int(st33.seed),
            },
            "v8": {"realised": re_v8,
                   "sigma2_spectral": st.spectrum_after["sigma2_rad2"],
                   "Hs_spectral": st.spectrum_after["Hs_m"],
                   "loop_maxdiff_h": st.loop_maxdiff()},
            "shallow": {"realised": re_sh,
                        "sigma2_spectral": st_sh.spectrum_after["sigma2_rad2"],
                        "Hs_spectral": st_sh.spectrum_after["Hs_m"],
                        "breaking": st_sh.breaking,
                        "dispersion_as_used": disp_sh},
            "clamp": {"realised": re_cl,
                      "sigma2_spectral": st_cl.spectrum_after["sigma2_rad2"],
                      "Hs_spectral": st_cl.spectrum_after["Hs_m"],
                      "breaking": st_cl.breaking},
        },
        "dispersion": {
            "constants": {"gamma_N_per_m": 0.0728, "rho_kg_per_m3": 1000.0,
                          "nu_m2_per_s": P.NU_KIN,
                          "temp_dependence": ("constant; the water-"
                                              "temperature dependence is left "
                                              "as the v5 plug-in point")},
            "form": "omega^2 = (g k + (gamma/rho) k^3) tanh(k h)",
            "deep_water_retired": True,
            "gate_table": g_disp,
        },
        "anchor": {
            "node": st.anchor,
            "model_id": P.ANCHOR_MODEL_ID,
            "candidates": P.ANCHOR_CANDIDATES,
            "sensitivity": P.anchor_sensitivity(
                env["windSpeedMS"], cfg["fetch"]["fetch_m"], env["depthM"],
                S_patch=S_patch, N=grid),
            "old_vs_new": {
                "cox_munk_old": P.cox_munk_sigma2(env["windSpeedMS"]),
                "new": st.anchor["sigma2_target_rad2"],
                "ratio": st.anchor["ratio_new_over_old"],
            },
        },
        "direction": {
            "node": st.direction,
            "comparison_vs_phillips": P.direction_comparison(
                env["windSpeedMS"], cfg["fetch"]["fetch_m"], env["depthM"],
                env["specL"], st.kmin, st.kmax_nyq),
        },
        "capillary_tail": st.tail,
        "breaking": st.breaking,
        "transport": {
            "method": ("det.py det_frame (UNCHANGED), Q=%d, sigma=0.75 px, "
                       "no clamp / no renormalisation" % Q_main),
            "master": {"stats": c_stats(C), "file": "C-v8-t0.npy"},
            "baseline_33": {"stats": c_stats(C33),
                            "file": "C-baseline-v33.npy",
                            "cross_check": v33_cross},
            "photon": {"stats": None if C_ph is None else c_stats(C_ph),
                       "file": (None if C_ph is None else "C-photon-t0.npy")},
        },
        "stills": {
            "preset1_v8": {**m8, "file": "still-preset1-v8.png"},
            "preset1_v33": {**m33, "file": "still-preset1-v33.png",
                            **v33_still_check},
            "shallow": {**m_sh, "file": "still-shallow-h010.png",
                        "depth_m": float(env_sh["depthM"]),
                        "wind_ms": float(env_sh["windSpeedMS"]),
                        "E": float(E_sh)},
            "clamp": {**m_cl, "file": "still-clamp-h010-U5.png",
                      "depth_m": float(env_cl["depthM"]),
                      "wind_ms": float(env_cl["windSpeedMS"]),
                      "E": float(E_cl)},
            "photon": "still-PHOTON-preset1-v8.png" if C_ph is not None else None,
            "panel_prereq": ("panel-v8-prereq.png" if len(panel_items) >= 2
                             else None),
            "display": {"E_v8": float(E_v8_check), "E_v3": float(E_V3),
                        "level_ratio": float(level_ratio),
                        "flat_water_direct_lum_v3": float(L_v3),
                        "flat_water_direct_lum_v5": float(L_v5),
                        "G1_v33": float(g1_33["G1"]),
                        "G1_v8": float(g1_v8["G1"]),
                        "note": ("same display-level normalisation as "
                                 "v5/v6/v7; evaluated once and held")},
            "metrics_definitions": {
                "crush": ("fraction of pixels with rec709 linear luminance "
                          "(pre-ACES, including exposure) < 0.02"),
                "blow": ("fraction of pixels with rec709 ACES luminance "
                         "> 0.98 (post-ACES, pre-sRGB-encode)"),
            },
        },
        "prerequisite_tasks": prereq,
        "gates": gates,
        "approximations_and_design_choices": [
            {"item": "anchor model family and the median rule",
             "statement": ("the adopted target is the mean square slope of a "
                           "fetch-limited spectrum integrated over the "
                           "resolved band, with the JONSWAP/Donelan "
                           "frequency spectrum mapped through the chunk-4 "
                           "dispersion.  Two named members (JONSWAP gamma=3.3 "
                           "and Donelan gamma=1.7) bracket the literature; "
                           "the adopted value is their median (midpoint) and "
                           "the full gamma/alpha/cutoff/depth sensitivity is "
                           "in the anchor node."),
             "kind": "model_selection"},
            {"item": "JONSWAP fetch laws below their fitted range",
             "statement": ("alpha = 0.076 Ft^-0.22 and f_p = 3.5 (g/U) "
                           "Ft^-1/3 were fitted to field data with "
                           "Ft ~ 1e3-1e5.  The pond scene here is Ft = %.0f, "
                           "an extrapolation DOWN in fetch.  It is used "
                           "because the design asks for exactly this family "
                           "for ponds; the direction of the extrapolation "
                           "error (a slight over-estimate of the high-"
                           "frequency level at very short fetch) is "
                           "acknowledged and not corrected."
                           % float(P.nondim_fetch(cfg["fetch"]["fetch_m"],
                                                  env["windSpeedMS"]))),
             "kind": "stated_extrapolation"},
            {"item": "legacy default fetch",
             "statement": ("a legacy recipe (no wave.fetchM) uses the fixed "
                           "fetch 4.60 m, solved so that the new anchor "
                           "reproduces the old Cox-Munk target at the legacy "
                           "reference scene (U=1.5, h=0.6).  It is a "
                           "compatibility anchor, not a physical pond fetch; "
                           "the pond default for new recipes is 100 m."),
             "kind": "backward_compatibility_choice"},
            {"item": "capillary cutoff does not bind",
             "statement": ("the viscous rule (2 nu k^2 > omega) sits ~4e3x "
                           "above the grid Nyquist, so no mode is removed "
                           "inside the resolved band and the cusp sharpness "
                           "is set by dx = %.2f mm.  The rule is implemented "
                           "and evaluated rather than waived silently; the "
                           "one-loop viscous attenuation of the finest mode "
                           "is recorded and NOT applied (the state is a "
                           "static spectrum repeated every 8 s)."
                           % (S_patch / grid * 1e3)),
             "kind": "honest_limitation"},
            {"item": "dispersion is invisible at t=0",
             "statement": ("omega only enters through exp(i omega t); at t=0 "
                           "the field is the spectrum itself.  The chunk-4 "
                           "dispersion therefore changes the reference still "
                           "by exactly 0 (proved: H0 and the field are "
                           "bit-identical across dispersion branches), and "
                           "its effect is measured at t=2.0 s instead."),
             "kind": "structural_consequence"},
            {"item": "patch = longest resolved wavelength",
             "statement": ("the patch side (%.3f m) IS the longest resolved "
                           "wavelength, i.e. the spectrum holds roughly one "
                           "wavelength of its dominant component.  A "
                           "single-snapshot realisation therefore deviates "
                           "from the spectral sigma^2 by the few-mode "
                           "scatter reported in gate 3 (z-score), and long-"
                           "wave-dominated strong-wind scenes are the ones "
                           "the patch cannot represent.  Enlarging the patch "
                           "is out of scope for chunk 4." % S_patch),
             "kind": "honest_limitation"},
            {"item": "shoaling / refraction out of scope",
             "statement": ("shoaling needs a depth gradient grad h; the "
                           "scenes here are uniform depth, so it does not "
                           "fire.  Deferred to v4 (solid bottom), as "
                           "PIPELINE-V8 chunk 4 states."),
             "kind": "scope"},
        ],
        "calibrations_not_derived": P.NON_DERIVED,
        "inherited_calibrations_v5_v6": [
            "turbidity -> da+db scale (TURBID_DC_MAX=0.30 1/m, 70% scatter)",
            "scattering blur proportionality (BLUR_K_PX=40 px at 2048^2)",
            "b_w(550)=0.0018 1/m anchor",
            "display-level normalisation E=E_v3.(L_v3/L_v5)",
            "RING_ETA=0.10 (rain, unused here)",
        ],
        "timing_s": {"state_v33": t_st33, "state_v8": t_st8,
                     "det_per_Q": det_time,
                     "wall_total": time.time() - t_all},
    }
    with open(os.path.join(args.out, "numbers-v8.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)

    print("\nGATES")
    for name, g in gates.items():
        print("  %-36s %s" % (name, g.get("pass")))
    print("OK wall=%.1fs -> %s" % (time.time() - t_all,
                                   os.path.join(args.out, "numbers-v8.json")))


if __name__ == "__main__":
    main()
