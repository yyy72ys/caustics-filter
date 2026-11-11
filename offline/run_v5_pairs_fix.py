"""v5 verification pairs, rebuild (chunk 1: water & bottom) -- pair fix pass.

Why this file exists
--------------------
The first v5 pair set (run_v5.py) put the *physically exact low range* on the
left of every pair:

    pair2  z = 0.5 m   (sand, wind mode)
    pair3  sand        (the reference substrate, gain 1.0 == untouched photo)
    pair4  wind mode   (H_s readout 0.084 m)

At those settings the master ratio C barely modulates (measured p1/p99 =
0.84/1.20 at z=0.5 m), so the left panel is indistinguishable from the
original photo and the pairs are not comparisons.  This pass rebuilds the
three pairs with BOTH sides inside the range where the term is visible.
Nothing about the design is changed: same recipe, preset 1, same derived
seed, same v5 physics/composite code, same display normalisation E_v5.

Pairs produced (new file names; nothing existing is overwritten):
  pair2b-depth-z1.5-vs-z4.png      z = 1.5 m vs z = 4.0 m   (sand, wind mode)
  pair3b-substrate-rock-vs-vegetation.png   rock vs vegetation
  pair4b-wave-fixed0.12-vs-fixed0.30.png    fixed H_s = 0.12 m vs 0.30 m
                                            (0.30 m is Stokes-capped: the
                                             effective height is in the label)

Singles, numbers-v5-pairs-fix.json and C evidence arrays go to the same out dir.

Usage (from caustics-filter/offline/):
  python run_v5_pairs_fix.py
"""
import argparse
import datetime
import gc
import json
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,
                        substrate_gain, scattering_blur_sigma_px,
                        glitter_slope_sigmas, cox_munk_sigma2,
                        CHANNEL_WL_NM)  # noqa: E402
from v5_state import V5WaterState  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,
                          image_metrics, v3_beer,
                          flat_water_direct_level)  # noqa: E402

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)

# Pair settings (task: pairs-fix.txt).  Everything else stays at the v5
# baseline (recipe-1, preset 1, derived seed, E_v5 display normalisation).
DEPTH_PAIR = (1.5, 4.0)
SUBSTRATE_PAIR = ("rock", "vegetation")
WAVE_PAIR_M = (0.12, 0.30)


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
def band_mod(a, k=8, win=9):
    """Relative modulation at the caustic scale (8x downsample, local mean).

    Dimensionless: std of (coarse local value / local mean).  Used as the
    quantitative stand-in for "can you see a light net at all".  Reference:
    the untouched photo measures ~0.26 on its own texture, so this metric is
    only used for the *master C* (which has no photo texture in it).
    """
    a = np.asarray(a, dtype=np.float64)
    if a.ndim == 3:
        a = (a * REC709).sum(-1)
    H, W = a.shape[:2]
    Hc, Wc = (H // k) * k, (W // k) * k
    m = a[:Hc, :Wc].reshape(Hc // k, k, Wc // k, k).mean(axis=(1, 3))
    from scipy import ndimage as ndi
    loc = ndi.uniform_filter(m, win)
    return float(np.std(m / np.maximum(loc, 1e-9)))


def c_vis(C):
    c = np.asarray(C, dtype=np.float64)
    return {"mean": float(c.mean()),
            "p1": float(np.percentile(c, 1)),
            "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "max": float(c.max()),
            "mod8": band_mod(c)}


def photo_diff(img_u8, photo_srgb):
    """How far the render is from the unmodified photo (sRGB, 0-1)."""
    d = np.abs(img_u8.astype(np.float64) / 255.0 - photo_srgb)
    return {"mean_abs": float(d.mean()),
            "p95_abs": float(np.percentile(d, 95))}


def save_pair(path, left, right, left_lines, right_lines, gap=10, label_h=30):
    """Side-by-side with two-line labels (new helper; run_v5.py untouched)."""
    a = Image.fromarray(np.asarray(left))
    b = Image.fromarray(np.asarray(right))
    W, H = a.size
    canvas = Image.new("RGB", (2 * W + gap, H + label_h), (250, 250, 250))
    canvas.paste(a, (0, label_h))
    canvas.paste(b, (W + gap, label_h))
    d = ImageDraw.Draw(canvas)
    for i, t in enumerate(left_lines):
        d.text((3, 3 + 13 * i), t, fill=(20, 20, 20))
    for i, t in enumerate(right_lines):
        d.text((W + gap + 3, 3 + 13 * i), t, fill=(20, 20, 20))
    canvas.save(path)


def fmt3(v):
    return "/".join("%.2f" % x for x in v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v5")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    env = resolve_environment(recipe)
    grid, Q = int(args.grid), int(args.q)

    pair_files = ["pair2b-depth-z1.5-vs-z4.png",
                  "pair3b-substrate-rock-vs-vegetation.png",
                  "pair4b-wave-fixed0.12-vs-fixed0.30.png"]
    still_files = ["still-2b-depth-z1.5.png", "still-2b-depth-z4.png",
                   "still-3b-substrate-rock.png",
                   "still-3b-substrate-vegetation.png",
                   "still-4b-wave-fixed-hs0.12.png",
                   "still-4b-wave-fixed-hs0.30.png"]
    npy_files = ["C-depth-z1.5.npy", "C-depth-z4.npy",
                 "C-wave-fixed0.12.npy", "C-wave-fixed0.30.npy"]
    num_file = "numbers-v5-pairs-fix.json"
    targets = [os.path.join(args.out, f) for f in
               pair_files + still_files + npy_files + [num_file]]
    existing = [p for p in targets if os.path.exists(p)]
    if existing:
        raise SystemExit("refusing to overwrite existing files: %s" % existing)

    # ------------------------------------------------------------- state
    t0 = time.time()
    state = V5WaterState(recipe, env, N=grid, seed=args.seed,
                         preset_index=args.preset)
    fields = state.fields_at(0.0)
    fields["S"] = state.S
    realised = state.realised_stats(fields)
    t_state = time.time() - t0
    print("state seed=%d Hs_spec=%.4f Hs_real=%.4f rms_slope=%.4f (%.1fs)"
          % (int(state.seed), state.spectrum_after["Hs_m"],
             realised["Hs_m"], realised["rms_slope"], t_state), flush=True)

    # ----------------------------------------------------------- physics
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_node = smith_g1_beckmann(realised["rms_slope"], sun_zenith)
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    # identical display normalisation to run_v5.py (E_v5 depends only on the
    # reference condition, so every panel keeps the baseline level reference)
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(g1_node["G1"])
    E_v5 = E_V3 * (L_v3 / L_v5)
    print("E_v5=%.4f G1=%.4f ref blur=%.2fpx" %
          (E_v5, g1_node["G1"], blur_ref["sigma_px"]), flush=True)

    # reference master C (same surface, reference depth): reuse the frozen
    # evidence array when present, otherwise recompute from this state.
    ref_c_path = os.path.join(args.out, "C-baseline-v5.npy")
    if os.path.exists(ref_c_path):
        C_ref = np.load(ref_c_path)
        if C_ref.shape != (grid, grid):
            raise SystemExit("C-baseline-v5.npy shape %s != grid %d"
                             % (C_ref.shape, grid))
        c_source = {"file": ref_c_path, "reused": True}
    else:
        C_ref, _ = det_frame(fields, state.S, grid, env, Q)
        np.save(os.path.join(args.out, "C-baseline-v5-pairsfix.npy"),
                np.asarray(C_ref, dtype=np.float32))
        c_source = {"file": "C-baseline-v5-pairsfix.npy", "reused": False}
    c_ref_stats = c_vis(C_ref)

    photo_u8 = np.asarray(Image.open(args.photo).convert("RGB"),
                          dtype=np.uint8)
    photo_srgb = photo_u8.astype(np.float64) / 255.0
    photo_lin = srgb_to_linear(photo_srgb)
    PH, PW, _ = photo_lin.shape
    gain_base = substrate_gain("sand", 1.0)

    metrics = {}

    def render(tag, f_use, C, env_use, gain_rgb, beer, blur_px, label):
        t = time.time()
        img, dbg = render_frame_v5(
            f_use, C, photo_lin, env_use, PW, PH, ambient_scale=A_V3,
            c_layer_exposure=E_v5, substrate_gain_rgb=gain_rgb,
            beer_down=beer["down"], beer_up=beer["up"], g1=g1_node["G1"],
            blur_sigma_px=blur_px, lowpass_sigma_px=B_SIGMA_V3,
            return_debug=True)
        fn = "still-%s.png" % tag
        Image.fromarray(img).save(os.path.join(args.out, fn))
        m = image_metrics(dbg["lin"], dbg["aces"])
        m.update({"tag": tag, "file": fn, "label": label,
                  "time_s": time.time() - t,
                  "beer_two_way_rgb": [float(x) for x in
                                       np.asarray(beer["down"])
                                       * np.asarray(beer["up"])],
                  "blur_sigma_px": float(blur_px),
                  "substrate_gain_rgb": [float(x) for x in gain_rgb],
                  "c_layer_exposure": float(E_v5),
                  "photo_diff_srgb": photo_diff(img, photo_srgb)})
        m["img_mean_srgb"] = [float(x) for x in img.reshape(-1, 3).mean(0)]
        metrics[tag] = m
        print("%-30s crush=%.4f blow=%.6f lum=%.4f (%.1fs) -> %s"
              % (tag, m["crush_rate_linear_lt_0.02"],
                 m["blowout_rate_aces_gt_0.98"], m["linear_luminance_mean"],
                 m["time_s"], fn), flush=True)
        return img

    # ------------------------------------------------------ 2b depth pair
    depth_info = []
    for z in DEPTH_PAIR:
        env_z = json.loads(json.dumps(env))
        env_z["depthM"] = float(z)
        t = time.time()
        C_z, info_z = det_frame(fields, state.S, grid, env_z, Q)
        cst = c_vis(C_z)
        cst["det_time_s"] = time.time() - t
        cst["det_info"] = info_z
        npy = "C-depth-z%g.npy" % z
        np.save(os.path.join(args.out, npy), np.asarray(C_z, dtype=np.float32))
        beer_z = beer_lambert(float(z), env["sun"]["altitudeDeg"], env["turbid"])
        blur_z = scattering_blur_sigma_px(att, float(z), grid)
        lab = "depth z=%.2f m, sand, wind mode" % z
        tag = "2b-depth-z%g" % z
        img = render(tag, fields, C_z, env_z, gain_base["gain_rgb"], beer_z,
                     blur_z["sigma_px"], lab)
        depth_info.append({
            "depth_m": float(z), "npy": npy, "C": cst, "still": "still-%s.png" % tag,
            "beer": {k: ([float(x) for x in v] if isinstance(v, np.ndarray) else v)
                     for k, v in beer_z.items()},
            "blur": blur_z, "label": lab, "metrics": metrics[tag]})
        del C_z
        gc.collect()

    # -------------------------------------------------- 3b substrate pair
    sub_info = []
    for sub in SUBSTRATE_PAIR:
        g = substrate_gain(sub, 1.0)
        tag = "3b-substrate-%s" % sub
        lab = "substrate %s (rho %s)" % (sub, fmt3(g["rho_rgb"]))
        img = render(tag, fields, C_ref, env, g["gain_rgb"], beer_ref,
                     blur_ref["sigma_px"], lab)
        sub_info.append({"substrate": sub, "gain": g, "label": lab,
                         "still": "still-%s.png" % tag,
                         "metrics": metrics[tag]})

    # ------------------------------------------------------- 4b wave pair
    wave_info = []
    del fields
    gc.collect()
    for hs in WAVE_PAIR_M:
        t = time.time()
        st = V5WaterState(recipe, env, N=grid, seed=args.seed,
                          preset_index=args.preset, height_mode="fixed",
                          height_m=float(hs))
        f_h = st.fields_at(0.0)
        f_h["S"] = st.S
        rl = st.realised_stats(f_h)
        t_s = time.time() - t
        t = time.time()
        C_h, info_h = det_frame(f_h, st.S, grid, env, Q)
        cst = c_vis(C_h)
        cst["det_time_s"] = time.time() - t
        cst["det_info"] = info_h
        npy = "C-wave-fixed%.2f.npy" % hs
        np.save(os.path.join(args.out, npy), np.asarray(C_h, dtype=np.float32))
        cap = st.height["stokes_cap"]
        eff = st.height["Hs_after"]
        if cap["applied"]:
            lab = ("fixed H_s=%.2f m -> effective %.3f m (Stokes cap "
                   "H/lam=1/7 at lam_p=%.2f m)" % (hs, eff, cap["lambda_peak_m"]))
        else:
            lab = "fixed H_s=%.2f m -> effective %.3f m (no cap)" % (hs, eff)
        tag = "4b-wave-fixed-hs%.2f" % hs
        img = render(tag, f_h, C_h, env, gain_base["gain_rgb"], beer_ref,
                     blur_ref["sigma_px"], lab)
        wave_info.append({
            "requested_m": float(hs), "effective_m": float(eff),
            "capped": bool(cap["applied"]),
            "Hs_cap_m": float(cap["Hs_cap_m"]),
            "lambda_peak_m": float(cap["lambda_peak_m"]),
            "factor": float(st.height["factor"]),
            "realised": rl, "state_build_s": t_s, "npy": npy,
            "C": cst, "still": "still-%s.png" % tag, "label": lab,
            "metrics": metrics[tag]})
        del f_h, C_h, st
        gc.collect()

    # -------------------------------------------------------------- pairs
    d0, d1 = depth_info
    save_pair(os.path.join(args.out, pair_files[0]),
              Image.open(os.path.join(args.out, d0["still"])),
              Image.open(os.path.join(args.out, d1["still"])),
              [d0["label"], "two-way beer %s" % fmt3(d0["beer"]["two_way"]) +
               "  blur %.1f px" % d0["blur"]["sigma_px"]],
              [d1["label"], "two-way beer %s" % fmt3(d1["beer"]["two_way"]) +
               "  blur %.1f px" % d1["blur"]["sigma_px"]])
    s0, s1 = sub_info
    save_pair(os.path.join(args.out, pair_files[1]),
              Image.open(os.path.join(args.out, s0["still"])),
              Image.open(os.path.join(args.out, s1["still"])),
              [s0["label"], "gain vs sand %s" % fmt3(s0["gain"]["gain_rgb"])],
              [s1["label"], "gain vs sand %s" % fmt3(s1["gain"]["gain_rgb"])])
    w0, w1 = wave_info
    save_pair(os.path.join(args.out, pair_files[2]),
              Image.open(os.path.join(args.out, w0["still"])),
              Image.open(os.path.join(args.out, w1["still"])),
              [w0["label"], "realised H_s %.3f m  slope rms %.3f" %
               (w0["realised"]["Hs_m"], w0["realised"]["rms_slope"])],
              [w1["label"], "realised H_s %.3f m  slope rms %.3f" %
               (w1["realised"]["Hs_m"], w1["realised"]["rms_slope"])])

    # ------------------------------------- reference (the complained panels)
    old_ref = {}
    for f, key in (("C-depth-z0.5.npy", "old_left_pair2_z0.5"),
                   ("C-baseline-v5.npy", "old_left_pair3_sand_and_pair4_wind")):
        p = os.path.join(args.out, f)
        if os.path.exists(p):
            c = np.load(p)
            old_ref[key] = {"file": f, "p1_p99": [float(np.percentile(c, 1)),
                                                  float(np.percentile(c, 99))],
                            "mod8": band_mod(c)}

    # ------------------------------------------------------------ numbers
    numbers = {
        "pipeline": "v5-chunk1",
        "purpose": ("rebuild of the three verification pairs so that BOTH "
                    "panels sit where the term is visible; the first set's "
                    "left panels (z=0.5 m / sand / wind mode) were inside the "
                    "range where C barely modulates and read as the original "
                    "photo"),
        "design": "PIPELINE-V5.md (unchanged)",
        "based_on": {"runner": "run_v5.py (unchanged, imported constants)",
                     "numbers": "numbers-v5.json"},
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": recipe.get("name"), "preset": args.preset,
        "seed": int(state.seed), "grid": [grid, grid], "Q": Q,
        "t_s": 0.0,
        "baseline_conditions": {
            "recipe": args.recipe, "photo": args.photo,
            "photo_size": [int(PW), int(PH)],
            "depth_reference_m": env["depthM"], "turbid": env["turbid"],
            "wind_ms": env["windSpeedMS"], "diffuse": env["diffuse"],
            "sunI": env["sunI"], "sun_alt_deg": env["sun"]["altitudeDeg"],
            "A_ambient": A_V3, "B_lowpass_sigma_px": B_SIGMA_V3,
            "c_layer_exposure_E_v5": float(E_v5),
            "E_v5_note": ("display normalisation recomputed with the same "
                          "formula as run_v5.py; identical value 2.708272"),
            "G1_smith": g1_node,
            "reference_blur": blur_ref,
            "reference_C": {**c_source, **c_ref_stats},
        },
        "wind_state": {
            "Hs_spec_m": state.spectrum_after["Hs_m"],
            "Hs_realised_m": realised["Hs_m"],
            "rms_slope": realised["rms_slope"],
            "sigma2_rad2": realised["sigma2_rad2"],
            "cox_munk_sigma2": cox_munk_sigma2(env["windSpeedMS"]),
            "calmness": state.calm,
            "spectrum_before": state.spectrum_before,
            "spectrum_after": state.spectrum_after,
        },
        "pairs": {
            "depth": {"panels": depth_info,
                      "note": ("same sand substrate, same wind-driven state, "
                               "same seed; only depth differs.  E_v5 is held "
                               "at the reference-condition value for both "
                               "panels, so the spectral Beer-Lambert darkening "
                               "with depth is shown, not compensated")},
            "substrate": {"panels": sub_info,
                          "note": ("rock and vegetation are both non-sand, so "
                                   "neither panel is the untouched photo; the "
                                   "reference master C (reference depth) is "
                                   "shared")},
            "wave": {"panels": wave_info,
                     "note": ("both panels are height-driven (fixed mode); "
                              "0.30 m is Stokes-capped and the label carries "
                              "the effective height.  The states are built "
                              "from the same recipe/seed; calmness "
                              "calibration is not applied in fixed mode "
                              "(as designed)")},
        },
        "old_left_panel_reference": old_ref,
        "metric_definitions": {
            "mod8": ("std of (8x-downsampled local value / local mean) of the "
                     "master C -- relative caustic modulation; the old left "
                     "panels measured 0.021-0.025"),
            "crush": "fraction of pixels with linear luminance < 0.02",
            "blow": "fraction with ACES luminance > 0.98",
            "photo_diff_srgb": ("mean / p95 absolute difference between the "
                                "rendered sRGB still and the unmodified photo"),
        },
        "timing_s": {"state_build_wind": t_state,
                     "wall_total": time.time() - t_all},
        "generated_files": pair_files + still_files + npy_files + [num_file],
    }
    with open(os.path.join(args.out, num_file), "w", encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)

    print("\nPAIRS")
    for f in pair_files:
        print("  %s (%.0f KB)" % (f, os.path.getsize(os.path.join(args.out, f)) / 1024))
    print("STILLS")
    for f in still_files:
        print("  %s (%.0f KB)" % (f, os.path.getsize(os.path.join(args.out, f)) / 1024))
    print("C EVIDENCE")
    for f in npy_files:
        print("  %s (%.0f KB)" % (f, os.path.getsize(os.path.join(args.out, f)) / 1024))
    print("NUMBERS %s" % num_file)
    print("\nMOD8 (master C): ref=%.4f | old z0.5=%.4f | z1.5=%.4f z4=%.4f | "
          "hs0.12=%.4f hs0.30=%.4f"
          % (c_ref_stats["mod8"],
             old_ref.get("old_left_pair2_z0.5", {}).get("mod8", float("nan")),
             depth_info[0]["C"]["mod8"], depth_info[1]["C"]["mod8"],
             wave_info[0]["C"]["mod8"], wave_info[1]["C"]["mod8"]))
    print("wall %.1fs" % (time.time() - t_all))


if __name__ == "__main__":
    main()
