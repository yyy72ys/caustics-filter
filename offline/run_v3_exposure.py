"""Pipeline v3 exposure probe (preset 1 / t=0, A2.6-Blow base).

Task: caustics-filter/offline v3 code, variant A2.6-Blow
(ambient x2.6, displacement lowpass sigma=3.0px) linear synthesis result,
apply exposure E=1.5/2.2/3.0 to the C-layer (display layer) ONLY.

Definition (this probe):
  base lin (A=2.6, Blow) via image.render_frame_v3(return_debug=True):
    col = alb*(direct+E_amb_A)*strength + sky_A
    lin_base = col * exposure_env
  where direct = Ebase[b]*C_sampled*beer[b] (the ONLY C-modulated term),
  E_amb_A = raw E_amb * 2.6, sky_A = sky*F*sky_gain * 2.6 (both fixed),
  sampling/shifts identical to v3 (Blow lowpass shared by albedo+C).
  probe:
    lin_E = (alb*(E*direct + E_amb_A)*strength + sky_A) * exposure_env
          = lin_base + (E-1) * alb*direct*strength*exposure_env
  aces_E = aces_narkowicz(max(lin_E,0)); sRGB PNG output.
  Ambient/sky stay at A=2.6; direct (C layer) alone carries E.
  Master (C-t0.npy, flow-t0.npy, state) is read-only, never written.

Outputs (out/v3/):
  still-A2.6-Blow-E1.5.png, -E2.2.png, -E3.0.png
  numbers-v3-exposure.json (append/merge by E; base record included)

Scope: preset 1, t=0 stills only. web/ untouched.

Usage (from caustics-filter/offline/):
  python run_v3_exposure.py --recipe ../blender/recipe-1.json --preset 1
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
from state import WaterState, sample_bicubic  # noqa: E402
from image import (render_frame_v3, srgb_to_linear, linear_to_srgb,  # noqa: E402
                   aces_narkowicz, e_base_dir_from_env, IOR, SKY_REF)

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)
A_BASE = 2.6
B_SIGMA = 3.0
E_LEVELS = (1.5, 2.2, 3.0)


def to_lum(rgb):
    return (np.asarray(rgb, dtype=np.float64) * REC709).sum(axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--out", default="out/v3")
    ap.add_argument("--v2dir", default="out/v2")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--exposures", default="1.5,2.2,3.0")
    args = ap.parse_args()
    exposures = tuple(float(x) for x in str(args.exposures).split(",") if x.strip())

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    if recipe.get("schema") != "caustics-recipe/v1":
        raise SystemExit("recipe schema mismatch")
    env = resolve_environment(recipe)

    # ---- frozen master reuse (read-only; never write to v2dir) ----
    c_path = os.path.join(args.v2dir, "C-t0.npy")
    c_stat_before = os.stat(c_path)
    C = np.load(c_path)  # read-only load; no save anywhere in this script
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise SystemExit("v2 C-t0.npy must be a square 2D master")
    with open(os.path.join(args.v2dir, "numbers-v2.json"), "r",
              encoding="utf-8") as f:
        v2n = json.load(f)

    t0 = time.time()
    state = WaterState(recipe, env, N=args.grid, seed=args.seed,
                       preset_index=args.preset)
    t_state = time.time() - t0
    fields = state.fields_at(0.0)
    fields["S"] = state.S
    if int(args.seed) if args.seed is not None else int(state.seed) != int(v2n["seed"]):
        print("WARNING: v3-exposure state seed %d != v2 seed %s"
              % (int(state.seed), v2n.get("seed")), flush=True)

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    # Base A2.6-Blow via canonical v3 path (bit-identical composition).
    tv = time.time()
    still_base, dbg = render_frame_v3(fields, C, photo_lin, env, PW, PH,
                                      ambient_scale=A_BASE,
                                      lowpass_sigma_px=B_SIGMA,
                                      return_debug=True)
    dt_base = time.time() - tv
    lin_base = dbg["lin"]
    print("base A2.6-Blow via render_frame_v3: %.1fs linmean=%.6f"
          % (dt_base, float(np.mean(lin_base))), flush=True)

    # ---- C-layer-only direct component map (identical sampling) ----
    # Recompute alb*direct*strength*exposure per-pixel (float64 HxWx3),
    # i.e. the linear contribution of the C-modulated direct term only.
    S = float(fields["S"])
    N_buf = C.shape[0]
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    Ebase = e_base_dir_from_env(env)
    absorb = np.asarray(env["absorb"], dtype=np.float64)

    gx_in = np.asarray(fields["gx"], dtype=np.float64)
    gy_in = np.asarray(fields["gy"], dtype=np.float64)
    n_state = gx_in.shape[0]
    sig_eff = B_SIGMA * (float(n_state) / 2048.0)
    from scipy import ndimage as _ndi
    gx_s = _ndi.gaussian_filter(gx_in, sig_eff, mode="wrap")
    gy_s = _ndi.gaussian_filter(gy_in, sig_eff, mode="wrap")

    iy, ix = np.mgrid[0:PH, 0:PW]
    # NOTE render_frame_v3 uses (W,H)=(PW,PH) and HxW lattices; keep same.
    H, W = PH, PW
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S
    gx = sample_bicubic(gx_s.astype(np.float32), px, py, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    beer = np.exp(-absorb * path_len * 0.22)
    Cc = np.asarray(C, dtype=np.float64)

    direct_lin = np.zeros((H, W, 3), dtype=np.float64)
    alb_all = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        su = np.clip(u + shx / S, 0.0, 1.0) * (PW - 1.0)
        sv = np.clip(v + shy / S, 0.0, 1.0) * (PH - 1.0)
        x0 = np.floor(su).astype(np.int64)
        y0 = np.floor(sv).astype(np.int64)
        x1 = np.clip(x0 + 1, 0, PW - 1)
        y1 = np.clip(y0 + 1, 0, PH - 1)
        x0 = np.clip(x0, 0, PW - 1)
        y0 = np.clip(y0, 0, PH - 1)
        fx = su - np.floor(su)
        fy = sv - np.floor(sv)
        alb = (photo_lin[y0, x0, b] * (1 - fx) * (1 - fy)
               + photo_lin[y0, x1, b] * fx * (1 - fy)
               + photo_lin[y1, x0, b] * (1 - fx) * fy
               + photo_lin[y1, x1, b] * fx * fy)
        alb_all[:, :, b] = alb
        bu = np.mod(u + shx / S, 1.0) * N_buf - 0.5
        bv = np.mod(v + shy / S, 1.0) * N_buf - 0.5
        bx0 = np.floor(bu).astype(np.int64) % N_buf
        by0 = np.floor(bv).astype(np.int64) % N_buf
        bx1 = (bx0 + 1) % N_buf
        by1 = (by0 + 1) % N_buf
        fxb = bu - np.floor(bu)
        fyb = bv - np.floor(bv)
        c = (Cc[by0, bx0] * (1 - fxb) * (1 - fyb)
             + Cc[by0, bx1] * fxb * (1 - fyb)
             + Cc[by1, bx0] * (1 - fxb) * fyb
             + Cc[by1, bx1] * fxb * fyb)
        direct = float(Ebase[b]) * c * float(beer[b])
        direct_lin[:, :, b] = alb * direct * strength * exposure

    # Consistency: lin_base must equal alb*(direct+E_amb_A)*strength*exposure+sky_A*exposure.
    # Recompute check (max abs diff) for the log, not enforced strictly.
    maxdiff_alb = float(np.max(np.abs(alb_all - dbg["alb"])))
    print("alb-map maxdiff vs render_frame_v3 debug: %.3g" % maxdiff_alb, flush=True)

    probes = []
    for E in exposures:
        tv = time.time()
        lin_E = lin_base + (float(E) - 1.0) * direct_lin
        aces_E = aces_narkowicz(np.maximum(lin_E, 0.0))
        out = (np.clip(linear_to_srgb(aces_E), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
        lin_lum = to_lum(lin_E)
        aces_lum = to_lum(np.clip(aces_E, 0.0, 1.0))
        linear_mean = float(np.mean(lin_E))
        linear_lum_mean = float(np.mean(lin_lum))
        crush = float(np.mean(lin_lum < 0.02))
        blow = float(np.mean(aces_lum > 0.98))
        dt = time.time() - tv
        # Filename: one decimal as specified (-E1.5 / -E2.2 / -E3.0).
        e_tag = ("%.1f" % float(E))
        fname = "still-A2.6-Blow-E%s.png" % e_tag
        Image.fromarray(out).save(os.path.join(args.out, fname))
        print("E=%s: linmean=%.6f crush=%.6f blow=%.6f (%.1fs) -> %s"
              % (e_tag, linear_mean, crush, blow, dt, fname), flush=True)
        probes.append({
            "E_C_layer": float(E),
            "base": "A2.6-Blow",
            "file": fname,
            "size": [int(PW), int(PH)],
            "metrics": {
                "linear_mean_post_exposure": linear_mean,
                "linear_luminance_mean": linear_lum_mean,
                "crush_rate_linear_lt_0.02": crush,
                "blowout_rate_aces_gt_0.98": blow,
            },
            "definitions": {
                "linear_mean": "mean of lin_E (pre-ACES linear incl. env exposure; C-modulated direct term xE, ambient/sky fixed at A=2.6; all channels)",
                "crush_rate": "fraction of pixels with rec709 lin_E luminance < 0.02",
                "blowout_rate": "fraction of pixels with rec709 ACES(lin_E) luminance > 0.98 (post-ACES, pre-sRGB-encode)",
            },
            "time_s": dt,
        })

    # ---- append/merge numbers file ----
    numbers_path = os.path.join(args.out, "numbers-v3-exposure.json")
    prev = None
    if os.path.exists(numbers_path):
        with open(numbers_path, "r", encoding="utf-8") as f:
            try:
                prev = json.load(f)
            except Exception:
                prev = None
    prev_probes = []
    if isinstance(prev, dict) and isinstance(prev.get("probes"), list):
        prev_probes = [p for p in prev["probes"]
                       if float(p.get("E_C_layer", -1)) not in
                       set(float(E) for E in exposures)]
    merged = prev_probes + probes
    merged.sort(key=lambda p: float(p["E_C_layer"]))
    numbers = {
        "pipeline": "v3-exposure-probe",
        "recipe": recipe.get("name"), "recipe_schema": recipe.get("schema"),
        "preset": args.preset, "seed": int(state.seed), "t_s": 0.0,
        "base": {
            "key": "A2.6-Blow", "A_ambient_scale": A_BASE,
            "B_band": "low", "lowpass_sigma_px_state": B_SIGMA,
            "file": "still-A2.6-Blow.png",
            "note": "canonical v3 still reused as base; this probe recomputes its lin via render_frame_v3 (read-only) and scales ONLY the C-modulated direct term",
        },
        "C_layer_exposure": {
            "levels": [float(E) for E in exposures],
            "applies_to": "direct term E_base,dir*C only (C-layer/display-layer gain); E_amb and Fresnel sky additive stay at A=2.6; env strength/exposure unchanged",
            "formula": "lin_E = (alb*(E*direct + E_amb_A2.6)*strength + sky_A2.6) * exposure_env",
            "master_untouched": {
                "v2_C_file": os.path.join(args.v2dir, "C-t0.npy"),
                "v2_C_shape": list(C.shape), "v2_C_dtype": str(C.dtype),
                "note": "master C/flow/state files are only read, never written by this probe",
            },
        },
        "metrics_definitions": probes[0]["definitions"] if probes else {},
        "probes": merged,
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "depthM": env["depthM"], "clarity": env.get("clarity"),
        },
        "timing_s": {
            "state_build": t_state,
            "base_render_v3_A2.6_Blow": dt_base,
            **{"probe_E%g" % float(p["E_C_layer"]): p["time_s"] for p in probes},
            "wall_total": time.time() - t_all,
        },
        "outputs": {
            "numbers": "numbers-v3-exposure.json",
            "stills": [p["file"] for p in merged],
            "still_size": [int(PW), int(PH)],
            "note": "filenames use one decimal (E1.5/E2.2/E3.0); timing keys use %g-style tags",
        },
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
    }
    with open(numbers_path, "w", encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)

    # Master-untouched guard: C file stat must be unchanged.
    c_stat_after = os.stat(c_path)
    assert (c_stat_after.st_mtime_ns == c_stat_before.st_mtime_ns
            and c_stat_after.st_size == c_stat_before.st_size), \
        "v2 master C-t0.npy was modified!"
    print("OK wall=%.1fs probes=%s"
          % (time.time() - t_all, [p["file"] for p in probes]), flush=True)


if __name__ == "__main__":
    main()
