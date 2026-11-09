"""Pipeline v3 runner (PIPELINE-V3, composition-adjustment round).

Frozen: state generation and transport are NOT rerun. The v2 master
C-t0.npy (and its recipe/seed/env) is reused; the water state is rebuilt
deterministically (same recipe+preset+seed) only to obtain the identical
t=0 fields for image formation. det.py / transport.py are never called.

Adjusted (image.py composition stage only): 6 still variants, preset 1 t=0:
  A (ambient scale): 1.0 / 1.7 / 2.6 applied to E_amb AND the camera-side
     additive (Fresnel sky term). The direct term E_base,dir*C is unchanged.
     NOTE: this path has no separate L_in-scatter buffer (underwater
     scattering is folded into E_amb); the sky term is the only other
     C-independent additive, so A scales both additives. Recorded in numbers.
  B (displacement band limit): full (as v2) / low (state gx,gy smoothed with
     a periodic Gaussian sigma=3.0 px at 2048^2 BEFORE shift computation;
     albedo wobble and C sampling share the smoothed shifts).

Per variant metrics -> numbers-v3.json:
  linear mean (pre-tonemap), crush rate (linear lum < 0.02),
  blowout rate (ACES lum > 0.98), fracture metric (Laplacian-energy ratio
  of wobbled albedo vs original photo; 1.0 = no shattering), run time.

Scope: t=0 stills only. No full-frame production. web/ untouched.

Usage (from caustics-filter/offline/):
  python run_v3.py --recipe ../blender/recipe-1.json --preset 1
Smoke:
  python run_v3.py --recipe ../blender/recipe-1.json --preset 1 --grid 256 --out out/v3-smoke --v2dir out/v2
  (smoke reuses the FULL-RES v2 master for stills; --grid only rebuilds the
  state fields -- for a like-for-like smoke use --grid 2048.)
"""
import argparse
import datetime
import json
import os
import sys
import time

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from state import WaterState  # noqa: E402
from image import render_frame_v3, srgb_to_linear  # noqa: E402

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)
A_LEVELS = (1.0, 1.7, 2.6)
B_LEVELS = (("full", 0.0), ("low", 3.0))


def to_lum(rgb):
    rgb = np.asarray(rgb, dtype=np.float64)
    return (rgb * REC709).sum(axis=-1)


def laplacian_energy(Y):
    """Mean squared 4-neighbour Laplacian (edge-replicated borders)."""
    Y = np.asarray(Y, dtype=np.float64)
    P = np.pad(Y, 1, mode="edge")
    lap = 4.0 * P[1:-1, 1:-1] - P[:-2, 1:-1] - P[2:, 1:-1] - P[1:-1, :-2] - P[1:-1, 2:]
    return float(np.mean(lap * lap))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--out", default="out/v3")
    ap.add_argument("--v2dir", default="out/v2")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    if recipe.get("schema") != "caustics-recipe/v1":
        raise SystemExit("recipe schema mismatch")
    env = resolve_environment(recipe)

    # ---- frozen master reuse (no transport rerun) ----
    c_path = os.path.join(args.v2dir, "C-t0.npy")
    C = np.load(c_path)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise SystemExit("v2 C-t0.npy must be a square 2D master")
    v2_numbers_path = os.path.join(args.v2dir, "numbers-v2.json")
    with open(v2_numbers_path, "r", encoding="utf-8") as f:
        v2n = json.load(f)

    t0 = time.time()
    state = WaterState(recipe, env, N=args.grid, seed=args.seed,
                       preset_index=args.preset)
    t_state = time.time() - t0
    fields = state.fields_at(0.0)
    fields["S"] = state.S

    if int(args.seed) if args.seed is not None else int(state.seed) != int(v2n["seed"]):
        print("WARNING: v3 state seed %d != v2 seed %s (recipe/preset drift?)"
              % (int(state.seed), v2n.get("seed")), flush=True)
    if C.shape != (args.grid, args.grid):
        print("NOTE: C master is %s but state grid is %d; image formation "
              "resamples both (v2 convention). For the canonical v3 run use "
              "--grid 2048." % (C.shape, args.grid), flush=True)

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape
    E_orig = laplacian_energy(to_lum(photo_lin))
    print("photo %dx%d E_orig(lap)=%.6g seed=%d S=%.4f"
          % (PW, PH, E_orig, int(state.seed), float(state.S)), flush=True)

    variants = []
    for A in A_LEVELS:
        for bname, sigma in B_LEVELS:
            key = "A%.1f-B%s" % (A, bname)
            fname = "still-A%.1f-B%s.png" % (A, bname)
            tv = time.time()
            still, dbg = render_frame_v3(fields, C, photo_lin, env, PW, PH,
                                         ambient_scale=A,
                                         lowpass_sigma_px=sigma,
                                         return_debug=True)
            dt = time.time() - tv
            lin, aces, alb = dbg["lin"], dbg["aces"], dbg["alb"]
            lin_lum = to_lum(lin)
            aces_lum = to_lum(np.clip(aces, 0.0, 1.0))
            linear_mean = float(np.mean(lin))
            linear_lum_mean = float(np.mean(lin_lum))
            crush = float(np.mean(lin_lum < 0.02))
            blow = float(np.mean(aces_lum > 0.98))
            E_wob = laplacian_energy(to_lum(alb))
            fracture = float(E_wob / E_orig) if E_orig > 0 else 0.0
            Image.fromarray(still).save(os.path.join(args.out, fname))
            print("%s: linmean=%.4f crush=%.4f blow=%.4f frac=%.4f (%.1fs) -> %s"
                  % (key, linear_mean, crush, blow, fracture, dt, fname),
                  flush=True)
            variants.append({
                "key": key, "A_ambient_scale": A,
                "B_band": bname, "lowpass_sigma_px_state": sigma,
                "file": fname, "size": [int(PW), int(PH)],
                "metrics": {
                    "linear_mean_pre_tonemap": linear_mean,
                    "linear_luminance_mean": linear_lum_mean,
                    "crush_rate_linear_lt_0.02": crush,
                    "blowout_rate_aces_gt_0.98": blow,
                    "fracture_laplacian_energy_ratio": fracture,
                    "fracture_E_wobbled": E_wob,
                    "fracture_E_original": E_orig,
                },
                "definitions": {
                    "linear_mean": "mean of lin (pre-ACES linear incl. exposure, all channels)",
                    "crush_rate": "fraction of pixels with rec709 lin luminance < 0.02",
                    "blowout_rate": "fraction of pixels with rec709 ACES luminance > 0.98 (post-ACES, pre-sRGB-encode)",
                    "fracture": "mean(Lap(Ywob)^2)/mean(Lap(Yorig)^2), 4-neighbour Laplacian, edge-replicated; Ywob = rec709 luminance of the wobbled albedo sampling (pre-lighting), Yorig = photo_lin luminance; 1.0 = no shattering",
                },
                "time_s": dt,
            })

    numbers = {
        "pipeline": "v3",
        "recipe": recipe.get("name"), "recipe_schema": recipe.get("schema"),
        "preset": args.preset, "seed": int(state.seed), "t_s": 0.0,
        "frozen": {
            "note": "state generation and transport frozen; det.py/transport.py never called; v2 master reused; state rebuilt deterministically for identical t=0 fields",
            "v2_C_file": os.path.join(args.v2dir, "C-t0.npy"),
            "v2_C_shape": list(C.shape), "v2_C_dtype": str(C.dtype),
            "v2_C_mean": float(C.astype(np.float64).mean()),
            "v2_Q_adopted": v2n["transport_main"]["Q_adopted"],
            "v2_E_base_dir_rgb": v2n["transport_main"]["E_base_dir_rgb"],
            "v2_E_base_dir_luminance": v2n["transport_main"]["E_base_dir_luminance"],
            "seed_match_v2": bool(int(state.seed) == int(v2n["seed"])),
            "state": {
                "grid": [args.grid, args.grid], "patch_S_m": state.S,
                "dx_m": state.dx,
            },
        },
        "variant_axes": {
            "A_ambient_scale": {
                "levels": list(A_LEVELS),
                "applies_to": "E_amb AND the Fresnel sky additive (camera-side additive); direct term E_base,dir*C is unchanged",
                "in_scatter_note": "this offline path has no separate L_in-scatter buffer (underwater scattering is folded into E_amb); the sky term is the only other C-independent additive, so A scales both additives",
            },
            "B_band": {
                "levels": [{"name": n, "sigma_px_state": s} for n, s in B_LEVELS],
                "method": "periodic (wrap) Gaussian on state gx,gy BEFORE normal/shift computation (sigma quoted at 2048^2, rescaled with grid); albedo wobble and C sampling share the smoothed shifts",
            },
        },
        "metrics_definitions": variants[0]["definitions"] if variants else {},
        "variants": variants,
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "depthM": env["depthM"], "clarity": env.get("clarity"),
        },
        "timing_s": {
            "state_build": t_state,
            **{"image_%s" % v["key"]: v["time_s"] for v in variants},
            "wall_total": time.time() - t_all,
        },
        "outputs": {
            "numbers": "numbers-v3.json",
            "stills": [v["file"] for v in variants],
            "still_size": [int(PW), int(PH)],
            "note": "photo-native size 735x539 (PIPELINE-V3 text says 736x540; outputs follow the input photo like v2)",
        },
    }
    with open(os.path.join(args.out, "numbers-v3.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)
    print("OK wall=%.1fs stills=%d size=%dx%d"
          % (time.time() - t_all, len(variants), PW, PH), flush=True)


if __name__ == "__main__":
    main()
