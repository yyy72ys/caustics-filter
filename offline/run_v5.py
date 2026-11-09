"""Pipeline v5 runner -- PIPELINE-V5 chunk 1 (water & bottom variables).

Implements the chunk-1 variable set from PIPELINE-V5.md on top of the frozen
v2/v3/v4 pipeline.  Nothing existing is modified: state.py, transport.py,
det.py, image.py, run_v2/v3/prod are imported, never edited.

New nodes / terms (all derived from report.md sections 1, 2, 5):
  1  H_s = 4*sqrt(m0): read-only node in wind mode; height-driven mode scales
     the spectrum so m0=(H_s/4)^2 and caps at the Stokes limit H/lambda<=1/7.
  2  calmness sigma^2 = integral k^2 S dk: realised readout + Cox-Munk
     reference; the high-k side is calibrated to the Cox-Munk value in wind
     mode.  Anisotropic glitter slope sigmas (sigma_a^2, sigma_c^2) are
     exposed (the offline v3 composite has no glitter term, so the sampler is
     implemented but not applied to the verification stills -- recorded).
  3  Smith G1(sun) on the direct term (statistical approximation).
  4  bottom albedo rho_b substrate presets as per-channel gains on the photo's
     linear values (direct-bottom term only).
  5  Beer-Lambert exp(-c s_down) * exp(-c s_up) with Pope & Fry a_w at
     ~650/550/470 nm + turbidity da/db, replacing the ad-hoc depth darkening.
  6  depth/contrast only through (i) the spectral attenuation and (ii) the
     turbidity scattering blur sigma ~ sqrt(b*z) on C.

Products in <out>/ (default out/v5/):
  numbers-v5.json
  still-v3-final-repro.png, still-v5-baseline.png
  still-depth-z0.5.png, still-depth-z3.0.png
  still-substrate-sand.png, still-substrate-vegetation.png
  still-wave-wind.png, still-wave-fixed-hs0.30.png
  pair1-baseline-v3-vs-v5.png ... pair4-wave-wind-vs-fixed.png
  C-depth-z0.5.npy, C-depth-z3.0.npy, C-wave-fixed.npy (transport evidence)

Usage (from caustics-filter/offline/):
  python run_v5.py --recipe ../blender/recipe-1.json --preset 1
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
from state import WaterState  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert,
                        beer_lambert_gate_example, smith_g1_beckmann,
                        substrate_gain, scattering_blur_sigma_px,
                        cox_munk_sigma2, cox_munk_anisotropy,
                        glitter_slope_sigmas, CHANNEL_WL_NM)  # noqa: E402
from v5_state import V5WaterState  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,
                          render_v3_final_repro, image_metrics, v3_beer,
                          flat_water_direct_level)  # noqa: E402
from recipe_v5 import parse_recipe_v5  # noqa: E402

# Reference numbers from PIPELINE-V5.md ("v3 final 7.9% / 0%").
V3_FINAL_CRUSH_REF = 0.079
V3_FINAL_BLOW_REF = 0.0
# Height-driven gate case from PIPELINE-V5.md.
WAVE_PAIR_TARGET_M = 0.30
# Additional (non-gated) demonstration case that stays under the Stokes cap.
WAVE_UNCAPPED_DEMO_M = 0.20
HS_GATE_TOL = 0.05


def to_lum(rgb):
    rec709 = np.array([0.2126, 0.7152, 0.0722])
    return (np.asarray(rgb, dtype=np.float64) * rec709).sum(axis=-1)


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {"shape": list(c.shape), "dtype": str(np.asarray(C).dtype),
            "mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min()),
            "mean_abs_dev_from_1": float(abs(c.mean() - 1.0))}


def save_side_by_side(path, left, right, left_label, right_label, gap=10,
                      label_h=18):
    a = Image.fromarray(left)
    b = Image.fromarray(right)
    W, H = a.size
    canvas = Image.new("RGB", (2 * W + gap, H + label_h), (250, 250, 250))
    canvas.paste(a, (0, label_h))
    canvas.paste(b, (W + gap, label_h))
    d = ImageDraw.Draw(canvas)
    d.text((3, 4), left_label, fill=(20, 20, 20))
    d.text((W + gap + 3, 4), right_label, fill=(20, 20, 20))
    canvas.save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v5")
    ap.add_argument("--v2dir", default="out/v2")
    ap.add_argument("--v3dir", default="out/v3")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--substrate-pair", default="sand,vegetation")
    ap.add_argument("--depth-pair", default="0.5,3.0")
    ap.add_argument("--wave-height", type=float, default=WAVE_PAIR_TARGET_M)
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    cfg = parse_recipe_v5(recipe)
    env = resolve_environment(recipe)
    grid = int(args.grid)
    Q = int(args.q)
    depth_pair = [float(x) for x in str(args.depth_pair).split(",") if x.strip()]
    sub_pair = [s.strip().lower() for s in str(args.substrate_pair).split(",")
                if s.strip()]

    print("recipe=%s schema=%s legacy=%s wave=%s bottom=%s"
          % (recipe.get("name"), cfg["schema"], cfg["legacy_recipe"],
             cfg["wave"], cfg["bottom"]), flush=True)

    # ---------------------------------------------------------------- state
    t0 = time.time()
    state = V5WaterState(recipe, env, N=grid, seed=args.seed,
                         preset_index=args.preset,
                         height_mode=cfg["wave"]["heightMode"],
                         height_m=cfg["wave"]["height_m"])
    t_state = time.time() - t0
    fields = state.fields_at(0.0)
    fields["S"] = state.S
    realised = state.realised_stats(fields)
    print("state seed=%d S=%.6f Hs_spec=%.4f Hs_real=%.4f sigma2_spec=%.5f "
          "sigma2_real=%.5f rms_slope=%.4f (%.1fs)"
          % (int(state.seed), state.S, state.spectrum_after["Hs_m"],
             realised["Hs_m"], state.spectrum_after["sigma2_rad2"],
             realised["sigma2_rad2"], realised["rms_slope"], t_state),
          flush=True)

    # ------------------------------------------------------------ physics
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_node = smith_g1_beckmann(realised["rms_slope"], sun_zenith)
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)

    # Display-level normalisation: keep the flat-water direct-term luminance of
    # the v3 final composition so the crush/blow comparison is about
    # distribution, not about an accidental exposure shift.  Display transform,
    # not physics.
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(g1_node["G1"])
    level_ratio = L_v3 / L_v5 if L_v5 > 0 else 1.0
    E_v5 = E_V3 * level_ratio
    print("beer two_way=%s g1=%.4f E_v3=%.2f level_ratio=%.4f E_v5=%.4f "
          "blur=%.2fpx" % ([round(float(x), 4) for x in beer_ref["two_way"]],
                           g1_node["G1"], E_V3, level_ratio, E_v5,
                           blur_ref["sigma_px"]), flush=True)

    # ------------------------------------------------------------- master C
    c_path = os.path.join(args.v2dir, "C-t0.npy")
    C_base = np.load(c_path)
    with open(os.path.join(args.v2dir, "numbers-v2.json"), "r",
              encoding="utf-8") as f:
        v2n = json.load(f)
    seed_match_v2 = bool(int(state.seed) == int(v2n["seed"]))
    if not seed_match_v2:
        print("WARNING: v5 state seed %d != v2 seed %s"
              % (int(state.seed), v2n.get("seed")), flush=True)
    if C_base.shape != (grid, grid):
        raise SystemExit("v2 C master %s does not match --grid %d"
                         % (C_base.shape, grid))
    C_stats = {"v2_master_v3_state": c_stats(C_base)}

    # The v5 calibrated state needs its own transport master: C is a property
    # of the surface, so reusing the v3-state master with v5 fields would be
    # physically inconsistent.  det.py is UNCHANGED; only the input state is.
    t1 = time.time()
    C_v5, info_v5 = det_frame(fields, state.S, grid, env, Q)
    C_stats["v5_baseline"] = {**c_stats(C_v5),
                              "det_time_s": time.time() - t1,
                              "det_info": info_v5}
    np.save(os.path.join(args.out, "C-baseline-v5.npy"),
            np.asarray(C_v5, dtype=np.float32))
    print("v5 C master mean=%.9f p99=%.3f max=%.3f (%.1fs)"
          % (float(C_v5.mean()), float(np.percentile(C_v5, 99)),
             float(C_v5.max()), time.time() - t1), flush=True)

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    stills = {}
    metrics = {}

    def render(tag, fields_use, C, env_use, gain_rgb, beer, g1, blur_px,
               E_use, label):
        tv = time.time()
        img, dbg = render_frame_v5(fields_use, C, photo_lin, env_use, PW, PH,
                                   ambient_scale=A_V3, c_layer_exposure=E_use,
                                   substrate_gain_rgb=gain_rgb,
                                   beer_down=beer["down"],
                                   beer_up=beer["up"], g1=g1,
                                   blur_sigma_px=blur_px,
                                   lowpass_sigma_px=B_SIGMA_V3,
                                   return_debug=True)
        dt = time.time() - tv
        fn = "still-%s.png" % tag
        Image.fromarray(img).save(os.path.join(args.out, fn))
        m = image_metrics(dbg["lin"], dbg["aces"])
        m["time_s"] = dt
        m["label"] = label
        m["file"] = fn
        m["beer_down_rgb"] = [float(x) for x in beer["down"]]
        m["beer_up_rgb"] = [float(x) for x in beer["up"]]
        m["G1"] = float(g1)
        m["blur_sigma_px"] = float(blur_px)
        m["c_layer_exposure"] = float(E_use)
        m["substrate_gain_rgb"] = [float(x) for x in gain_rgb]
        stills[tag] = img
        metrics[tag] = m
        print("%-22s crush=%.4f blow=%.6f lum=%.4f (%.1fs) -> %s"
              % (tag, m["crush_rate_linear_lt_0.02"],
                 m["blowout_rate_aces_gt_0.98"], m["linear_luminance_mean"],
                 dt, fn), flush=True)
        return img

    # ---- still 1a: v3 final reproduction (same conditions) ----
    # NOTE: the v3 final used the *v3 (uncalibrated) state*.  Build it
    # separately; the v5 baseline below uses the v5 calibrated state.
    t0 = time.time()
    state_v3 = WaterState(recipe, env, N=grid, seed=args.seed,
                          preset_index=args.preset)
    fields_v3 = state_v3.fields_at(0.0)
    fields_v3["S"] = state_v3.S
    v3_img, v3_dbg = render_v3_final_repro(fields_v3, C_base, photo_lin, env,
                                           PW, PH)
    del fields_v3, state_v3
    gc.collect()
    Image.fromarray(v3_img).save(os.path.join(args.out,
                                              "still-v3-final-repro.png"))
    v3_metrics = image_metrics(v3_dbg["lin"], v3_dbg["aces"])
    v3_metrics["file"] = "still-v3-final-repro.png"
    stills["v3_final_repro"] = v3_img
    metrics["v3_final_repro"] = v3_metrics
    v3_ref_path = os.path.join(args.v3dir, "still-A2.6-Blow-E2.2.png")
    v3_ref_diff = None
    if os.path.exists(v3_ref_path):
        ref = np.asarray(Image.open(v3_ref_path).convert("RGB"),
                         dtype=np.int16)
        v3_ref_diff = {"file": v3_ref_path,
                       "max_abs_diff_uint8": int(np.max(np.abs(ref.astype(np.int16)
                                                             - v3_img.astype(np.int16)))),
                       "mean_abs_diff_uint8": float(np.mean(np.abs(ref.astype(np.int16)
                                                                  - v3_img.astype(np.int16))))}
    print("v3 final repro crush=%.4f blow=%.6f (ref png diff %s)"
          % (v3_metrics["crush_rate_linear_lt_0.02"],
             v3_metrics["blowout_rate_aces_gt_0.98"], v3_ref_diff), flush=True)

    # ---- still 1b: v5 baseline ----
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    render("v5-baseline", fields, C_v5, env, gain_base["gain_rgb"], beer_ref,
           g1_node["G1"], blur_ref["sigma_px"], E_v5,
           "v5 baseline")

    # diagnostic: same v5 physics at the frozen v3 exposure E=2.2 (no display
    # normalisation) -- recorded in numbers, not a verification still
    _diag_img, _diag_dbg = render_frame_v5(
        fields, C_v5, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_V3, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=g1_node["G1"], blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    baseline_at_E_v3 = image_metrics(_diag_dbg["lin"], _diag_dbg["aces"])
    img_level = {
        "direct_lin_mean_v3": float(np.mean(v3_dbg["direct_lin"])),
        "direct_lin_mean_v5": float(np.mean(_diag_dbg["direct_lin"])),
    }
    img_level["ratio"] = (img_level["direct_lin_mean_v3"]
                          / img_level["direct_lin_mean_v5"])
    del _diag_img, _diag_dbg
    print("v5 baseline at frozen E=2.2 crush=%.4f blow=%.6f"
          % (baseline_at_E_v3["crush_rate_linear_lt_0.02"],
             baseline_at_E_v3["blowout_rate_aces_gt_0.98"]), flush=True)

    # ---- still 2: depth pair ----
    for z in depth_pair:
        env_z = json.loads(json.dumps(env))
        env_z["depthM"] = float(z)
        tz = time.time()
        C_z, info_z = det_frame(fields, state.S, grid, env_z, Q)
        C_stats["depth_z%g" % z] = {**c_stats(C_z),
                                    "det_time_s": time.time() - tz,
                                    "det_info": info_z}
        np.save(os.path.join(args.out, "C-depth-z%g.npy" % z),
                np.asarray(C_z, dtype=np.float32))
        beer_z = beer_lambert(float(z), env["sun"]["altitudeDeg"],
                              env["turbid"])
        blur_z = scattering_blur_sigma_px(att, float(z), grid)
        render("depth-z%g" % z, fields, C_z, env_z, gain_base["gain_rgb"],
               beer_z, g1_node["G1"], blur_z["sigma_px"], E_v5,
               "depth z=%.2f m" % z)
        del C_z
        gc.collect()

    # ---- still 3: substrate pair ----
    for sub in sub_pair:
        g = substrate_gain(sub, cfg["bottom"]["albedoScale"])
        render("substrate-%s" % sub, fields, C_v5, env, g["gain_rgb"],
               beer_ref, g1_node["G1"], blur_ref["sigma_px"], E_v5,
               "substrate %s" % sub)

    # ---- still 4a: wind-driven wave panel (readout H_s) ----
    render("wave-wind", fields, C_v5, env, gain_base["gain_rgb"], beer_ref,
           g1_node["G1"], blur_ref["sigma_px"], E_v5,
           "wind mode (H_s readout)")

    wave_wind_node = state.wave_node()
    calm_node = state.calmness_node(fields)
    wind_state_info = {
        "seed": int(state.seed), "grid": [grid, grid], "S_m": float(state.S),
        "dx_m": float(state.dx), "realised": realised,
        "Hs_spec": state.spectrum_after["Hs_m"],
    }

    # free the wind state before building the height-driven state
    del fields, state
    gc.collect()

    # ---- still 4b: height-driven state (same seed) ----
    t1 = time.time()
    state_h = V5WaterState(recipe, env, N=grid, seed=args.seed,
                           preset_index=args.preset, height_mode="fixed",
                           height_m=float(args.wave_height))
    t_state_h = time.time() - t1
    fields_h = state_h.fields_at(0.0)
    fields_h["S"] = state_h.S
    realised_h = state_h.realised_stats(fields_h)
    print("fixed H_s: requested=%.3f after=%.4f cap=%.4f realised=%.4f"
          % (float(args.wave_height), state_h.height["Hs_after"],
             state_h.height["stokes_cap"]["Hs_cap_m"], realised_h["Hs_m"]),
          flush=True)
    t1 = time.time()
    C_h, info_h = det_frame(fields_h, state_h.S, grid, env, Q)
    C_stats["wave_fixed"] = {**c_stats(C_h), "det_time_s": time.time() - t1,
                             "det_info": info_h}
    np.save(os.path.join(args.out, "C-wave-fixed.npy"),
            np.asarray(C_h, dtype=np.float32))
    render("wave-fixed-hs%.2f" % float(args.wave_height), fields_h, C_h, env,
           gain_base["gain_rgb"], beer_ref, g1_node["G1"],
           blur_ref["sigma_px"], E_v5,
           "fixed mode H_s=%.2f m" % float(args.wave_height))
    fixed_state_info = {
        "seed": int(state_h.seed), "grid": [grid, grid],
        "S_m": float(state_h.S), "dx_m": float(state_h.dx),
        "realised": realised_h,
        "height": state_h.height,
        "spectrum_before": state_h.spectrum_before,
        "spectrum_after": state_h.spectrum_after,
        "calmness_node": state_h.calmness_node(fields_h),
    }
    del fields_h, C_h, state_h
    gc.collect()

    # ---- H_s gate: an uncapped height-driven case (state only, no transport)
    t1 = time.time()
    state_u = V5WaterState(recipe, env, N=grid, seed=args.seed,
                           preset_index=args.preset, height_mode="fixed",
                           height_m=float(WAVE_UNCAPPED_DEMO_M))
    fields_u = state_u.fields_at(0.0)
    realised_u = state_u.realised_stats(fields_u)
    uncapped_info = {
        "requested_m": float(WAVE_UNCAPPED_DEMO_M),
        "after_m": state_u.height["Hs_after"],
        "realised_m": realised_u["Hs_m"],
        "capped": state_u.height["stokes_cap"]["applied"],
        "time_s": time.time() - t1,
    }
    del fields_u, state_u
    gc.collect()

    # ------------------------------------------------------------- gates
    def dev(a, b):
        return float((a - b) / b) if b else None

    hs_fixed = {
        "requested_m": float(args.wave_height),
        "realised_m": realised_h["Hs_m"],
        "uncapped_request_m": fixed_state_info["height"]["Hs_uncapped_requested"],
        "cap_m": fixed_state_info["height"]["stokes_cap"]["Hs_cap_m"],
        "cap_applied": fixed_state_info["height"]["stokes_cap"]["applied"],
        "lambda_peak_m": fixed_state_info["height"]["stokes_cap"]["lambda_peak_m"],
    }
    # realised vs the requested 0.30 m
    hs_fixed["dev_vs_requested"] = dev(hs_fixed["realised_m"],
                                       hs_fixed["requested_m"])
    hs_fixed["pass_vs_requested"] = bool(
        abs(hs_fixed["dev_vs_requested"]) <= HS_GATE_TOL)
    # realised vs the post-cap effective target (the physically reachable one)
    eff = min(hs_fixed["requested_m"], hs_fixed["cap_m"])
    hs_fixed["effective_target_m"] = eff
    hs_fixed["dev_vs_effective"] = dev(hs_fixed["realised_m"], eff)
    hs_fixed["pass_vs_effective"] = bool(
        abs(hs_fixed["dev_vs_effective"]) <= HS_GATE_TOL)
    hs_uncapped = {
        **uncapped_info,
        "dev_vs_requested": dev(uncapped_info["realised_m"],
                                uncapped_info["requested_m"]),
        "pass_vs_requested": bool(
            abs(dev(uncapped_info["realised_m"],
                    uncapped_info["requested_m"])) <= HS_GATE_TOL),
    }

    cmean_gate = {}
    for k, v in C_stats.items():
        cmean_gate[k] = {"mean": v["mean"],
                         "abs_dev": v["mean_abs_dev_from_1"],
                         "pass": bool(abs(v["mean"] - 1.0) <= 1e-3)}

    base_m = metrics["v5-baseline"]
    crush_gate = {
        "reference": {"crush": V3_FINAL_CRUSH_REF, "blow": V3_FINAL_BLOW_REF,
                      "source": "PIPELINE-V5.md line 48 (v3 final 7.9% / 0%)"},
        "v3_final_repro": {
            "crush": v3_metrics["crush_rate_linear_lt_0.02"],
            "blow": v3_metrics["blowout_rate_aces_gt_0.98"]},
        "v5_baseline": {
            "crush": base_m["crush_rate_linear_lt_0.02"],
            "blow": base_m["blowout_rate_aces_gt_0.98"]},
        "v5_baseline_at_frozen_E2.2": baseline_at_E_v3,
        "pass_vs_design_ref": bool(
            base_m["crush_rate_linear_lt_0.02"] <= V3_FINAL_CRUSH_REF
            and base_m["blowout_rate_aces_gt_0.98"] <= V3_FINAL_BLOW_REF),
        "pass_vs_v3_repro": bool(
            base_m["crush_rate_linear_lt_0.02"]
            <= v3_metrics["crush_rate_linear_lt_0.02"]
            and base_m["blowout_rate_aces_gt_0.98"]
            <= max(v3_metrics["blowout_rate_aces_gt_0.98"], 0.0)),
        "pass": bool(base_m["crush_rate_linear_lt_0.02"] <= V3_FINAL_CRUSH_REF
                     and base_m["blowout_rate_aces_gt_0.98"] <= V3_FINAL_BLOW_REF),
        "note": ("the v3 final reproduction measures 8.71% crush (the design "
                 "document quotes 7.9%); both references are reported.  "
                 "E_v5 is a display-level normalisation, so the un-normalised "
                 "E=2.2 diagnostic is also recorded."),
    }

    beer_gate = beer_lambert_gate_example()
    beer_gate["pass"] = bool(beer_gate["pass_550nm"] and beer_gate["pass_650nm"])

    gates = {
        "H_s_height_mode": {
            "gate": "realised H_s within +-5% of the target",
            "case_fixed_0.30m": hs_fixed,
            "case_uncapped_0.20m": hs_uncapped,
            "pass": bool(hs_fixed["pass_vs_requested"]),
            "pass_uncapped_demo": bool(hs_uncapped["pass_vs_requested"]),
            "pass_vs_effective_target": bool(hs_fixed["pass_vs_effective"]),
            "pass_note": ("the specified 0.30 m case is Stokes-capped "
                          "(H/lambda<=1/7) at the realised dominant "
                          "wavelength, so the literal +-5% gate on 0.30 m "
                          "FAILS by construction (realised 0.2102 m, "
                          "-29.9%).  The cap is the designed physical "
                          "behaviour and is recorded, not hidden.  The same "
                          "mode meets +-5% of the post-cap effective target "
                          "and an uncapped 0.20 m case meets +-5% of its "
                          "target.  All three numbers are reported."),
        },
        "beer_lambert_example": beer_gate,
        "C_mean_1": {"gate": "mean(C) = 1.0 (+-1e-3, energy conservation)",
                     "per_master": cmean_gate,
                     "pass": bool(all(v["pass"] for v in cmean_gate.values()))},
        "crush_blow_vs_v3_final": crush_gate,
    }

    # ------------------------------------------------------------- pairs
    pair_files = []
    p = "pair1-baseline-v3-vs-v5.png"
    save_side_by_side(os.path.join(args.out, p), stills["v3_final_repro"],
                      stills["v5-baseline"], "v3 final (repro)", "v5 baseline")
    pair_files.append(p)
    p = "pair2-depth-z%g-vs-z%g.png" % (depth_pair[0], depth_pair[1])
    save_side_by_side(os.path.join(args.out, p),
                      stills["depth-z%g" % depth_pair[0]],
                      stills["depth-z%g" % depth_pair[1]],
                      "depth z=%.2f m" % depth_pair[0],
                      "depth z=%.2f m" % depth_pair[1])
    pair_files.append(p)
    p = "pair3-substrate-%s-vs-%s.png" % (sub_pair[0], sub_pair[1])
    save_side_by_side(os.path.join(args.out, p),
                      stills["substrate-%s" % sub_pair[0]],
                      stills["substrate-%s" % sub_pair[1]],
                      "substrate %s" % sub_pair[0],
                      "substrate %s" % sub_pair[1])
    pair_files.append(p)
    p = "pair4-wave-wind-vs-fixed.png"
    save_side_by_side(os.path.join(args.out, p), stills["wave-wind"],
                      stills["wave-fixed-hs%.2f" % float(args.wave_height)],
                      "wind mode H_s=%.3f m" % realised["Hs_m"],
                      "fixed H_s=%.2f m (cap %.3f)"
                      % (float(args.wave_height),
                         fixed_state_info["height"]["stokes_cap"]["Hs_cap_m"]))
    pair_files.append(p)

    # ---- reproducible caustics-recipe/v2 copies (new fields only) ----
    rec_dir = os.path.join(args.out, "recipes")
    os.makedirs(rec_dir, exist_ok=True)

    def write_recipe(name, wave, bottom, depth=None):
        r2 = json.loads(json.dumps(recipe))
        r2["schema"] = "caustics-recipe/v2"
        r2["wave"] = wave
        r2["bottom"] = bottom
        if depth is not None:
            r2["water"]["depthM"] = float(depth)
        with open(os.path.join(rec_dir, name), "w", encoding="utf-8") as f:
            json.dump(r2, f, indent=2, ensure_ascii=False)

    write_recipe("baseline-v2.json", {"heightMode": "wind"},
                 {"substrate": None, "albedoScale": 1.0})
    for sub in sub_pair:
        write_recipe("substrate-%s-v2.json" % sub, {"heightMode": "wind"},
                     {"substrate": sub,
                      "albedoScale": cfg["bottom"]["albedoScale"]})
    write_recipe("wave-wind-v2.json", {"heightMode": "wind"},
                 {"substrate": None, "albedoScale": 1.0})
    write_recipe("wave-fixed-hs%.2f-v2.json" % float(args.wave_height),
                 {"heightMode": "fixed",
                  "height_m": float(args.wave_height)},
                 {"substrate": None, "albedoScale": 1.0})
    for z in depth_pair:
        write_recipe("depth-z%g-v2.json" % z, {"heightMode": "wind"},
                     {"substrate": None, "albedoScale": 1.0}, depth=z)

    # ------------------------------------------------------------ numbers
    numbers = {
        "pipeline": "v5-chunk1",
        "design": "PIPELINE-V5.md (caustics v3 chunk 1: water & bottom)",
        "evidence": "report.md sections 1, 2, 5",
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": recipe.get("name"), "recipe_schema": recipe.get("schema"),
        "preset": args.preset, "seed": int(wind_state_info["seed"]),
        "t_s": 0.0,
        "recipe_v2_extension": {
            "parsed": cfg,
            "legacy_recipe": cfg["legacy_recipe"],
            "note": ("v1 recipes without wave/bottom blocks read exactly as "
                     "before: identity substrate gain, wind height mode"),
        },
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "depthM": env["depthM"], "turbid": env["turbid"],
            "clarity": env.get("clarity"), "wind_ms": env["windSpeedMS"],
            "photo": args.photo, "photo_size": [int(PW), int(PH)],
        },
        "state": {
            "grid": [grid, grid], "patch_S_m": float(wind_state_info["S_m"]),
            "dx_m": float(wind_state_info["dx_m"]),
            "seed_match_v2": seed_match_v2,
            "wind_mode": wind_state_info,
            "calmness_node": calm_node,
            "H_s_node": wave_wind_node,
            "fixed_mode": fixed_state_info,
            "uncapped_demo": uncapped_info,
        },
        "physics": {
            "channel_wavelengths_nm": list(CHANNEL_WL_NM),
            "water_attenuation": {
                "a_w": [float(x) for x in att["a_w"]],
                "b_w": [float(x) for x in att["b_w"]],
                "da": [float(x) for x in att["da"]],
                "db": [float(x) for x in att["db"]],
                "c": [float(x) for x in att["c"]],
                "turbidity_calibration": {
                    "dc_max_per_m": 0.30,
                    "scatter_fraction": 0.7,
                    "note": "named calibration: the report gives the form "
                            "(turbidity adds da+db) but no scale"},
                "source": "Pope & Fry 1997 table (report section 2.2); "
                          "b_w(550)=0.0018 1/m, Rayleigh lambda^-4.32",
            },
            "beer_lambert_reference_condition": {
                **{k: ([float(x) for x in v] if isinstance(v, np.ndarray)
                       else v) for k, v in beer_ref.items()},
                "display_level_normalisation": {
                    "flat_water_direct_lum_v3": float(L_v3),
                    "flat_water_direct_lum_v5": float(L_v5),
                    "ratio": float(level_ratio),
                    "image_level_check": img_level,
                    "E_v3": float(E_V3), "E_v5": float(E_v5),
                    "note": "display transform, not physics: keeps the "
                            "flat-water direct-term luminance of the v3 "
                            "final composition so the crush/blow gate "
                            "compares distributions, not an exposure shift",
                },
            },
            "smith_G1": {**g1_node,
                         "note": "statistical Smith approximation (report "
                                 "1.3); ray-exact visibility NOT implemented; "
                                 "the two are never combined"},
            "scattering_blur": blur_ref,
            "substrates": {
                name: substrate_gain(name, 1.0)
                for name in ("sand", "rock", "vegetation", "mud")},
            "glitter": {
                "slope_sigmas": glitter_slope_sigmas(env["windSpeedMS"]),
                "cox_munk_total": cox_munk_sigma2(env["windSpeedMS"]),
                "applied_to_stills": False,
                "note": ("the offline v3 composite has NO glitter term (only "
                         "the Fresnel sky reflection); the anisotropic "
                         "(sigma_a^2, sigma_c^2) Gaussian sampler is "
                         "implemented in v5_physics/v5_composite but is off "
                         "in the verification stills so that still 1 stays "
                         "the v3 final condition.  Wiring it into a glitter "
                         "term is the composite-stage follow-up."),
            },
        },
        "transport": {
            "method": "det.py det_frame (UNCHANGED), Q=%d, sigma=0.75px, "
                      "no clamp / no renormalisation" % Q,
            "v3_master_for_repro": {
                "file": c_path, "source": "v2 frozen master (v3 state)",
                "used_for": "still 1a v3-final reproduction only"},
            "v5_master": {
                "file": "C-baseline-v5.npy",
                "source": "recomputed from the v5 calibrated state",
                "note": "C is a property of the surface; reusing the v3-state "
                        "master with v5 fields would be inconsistent"},
            "per_master": C_stats,
        },
        "stills": {
            "individual": {k: v["file"] for k, v in metrics.items()},
            "pairs": pair_files,
            "size": [int(PW), int(PH)],
            "metrics": metrics,
            "metrics_definitions": {
                "crush": "fraction of pixels with rec709 linear luminance "
                         "(pre-ACES, incl. exposure) < 0.02",
                "blow": "fraction of pixels with rec709 ACES luminance "
                        "> 0.98 (post-ACES, pre-sRGB-encode)",
            },
            "v3_final_repro_vs_reference_png": v3_ref_diff,
            "notes": {
                "depth_pair_display_gain": ("E_v5 is held at the "
                                            "reference-condition value for "
                                            "both depth panels, so the "
                                            "spectral Beer-Lambert darkening "
                                            "with depth is visible, not "
                                            "compensated"),
                "substrate_reference": ("sand is the reference substrate "
                                        "(gain 1.0), so the sand panel is "
                                        "identical to the v5 baseline"),
                "wave_pair": ("wind mode shows the H_s readout; fixed mode "
                              "requests 0.30 m and is Stokes-capped to the "
                              "reported effective height"),
            },
        },
        "gates": gates,
        "timing_s": {
            "state_build_wind": t_state, "state_build_fixed": t_state_h,
            "wall_total": time.time() - t_all,
        },
        "calibrations_not_derived": [
            "turbidity -> da+db scale (TURBID_DC_MAX=0.30 1/m, 70% scattering)",
            "scattering blur proportionality (BLUR_K_PX=40 px at 2048^2)",
            "b_w(550)=0.0018 1/m anchor (report gives order 0.001-0.01 only)",
            "display-level normalisation E_v5=E_v3*(L_v3/L_v5) (display "
            "transform, not physics)",
        ],
    }
    with open(os.path.join(args.out, "numbers-v5.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)
    print("OK wall=%.1fs gates: Hs_0.30m(literal)=%s Hs_0.20m=%s beer=%s "
          "Cmean=%s crushblow=%s"
          % (time.time() - t_all,
             gates["H_s_height_mode"]["pass"],
             gates["H_s_height_mode"]["pass_uncapped_demo"],
             gates["beer_lambert_example"]["pass"],
             gates["C_mean_1"]["pass"],
             gates["crush_blow_vs_v3_final"]["pass"]), flush=True)


if __name__ == "__main__":
    main()
