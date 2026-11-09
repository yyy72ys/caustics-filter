"""Offline pipeline v1 entry point.

Usage (from caustics-filter/offline/):
  python run.py --recipe ../blender/recipe-1.json --preset 1
  python run.py --recipe ../blender/recipe-1.json --preset 1 --grid 256 --photons 20000 --frames 2 --out out-smoke

Re-executable: same recipe + preset (+ optional seed) reproduces the outputs.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from state import WaterState  # noqa: E402
from transport import photon_frame  # noqa: E402
from image import render_frame, srgb_to_linear  # noqa: E402

FRAME_DT = 0.05  # 20 fps-equivalent stepping
VIDEO_FPS = 20.0


def lum_stats(a):
    a = np.asarray(a, dtype=np.float64).ravel()
    mean = float(np.mean(a))
    std = float(np.std(a))
    return {
        "method": "rec709 luminance of linear values, flattened over crop/field",
        "mean": mean,
        "cv": float(std / mean) if mean > 0 else 0.0,
        "p50": float(np.percentile(a, 50)),
        "p99": float(np.percentile(a, 99)),
        "max": float(np.max(a)),
    }


def to_lum(rgb):
    rgb = np.asarray(rgb, dtype=np.float64)
    return 0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--photons", type=int, default=4000000)
    ap.add_argument("--frames", type=int, default=16)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--out", default="out")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--ref1", default="../docs/refs/ref1-green-ripples.jpg")
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    if recipe.get("schema") != "caustics-recipe/v1":
        raise SystemExit("recipe schema mismatch")
    env = resolve_environment(recipe)

    t0 = time.time()
    state = WaterState(recipe, env, N=args.grid, seed=args.seed,
                       preset_index=args.preset)
    t_state = time.time() - t0

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    loop_diff = state.loop_maxdiff()

    E0 = None
    transport_info0 = None
    t_transport = 0.0
    t_image = 0.0
    frames_u8 = []
    VW, VH = (PW + 1) // 2 * 2, (PH + 1) // 2 * 2  # even dims for mp4 (736x540)
    for i in range(args.frames):
        t = i * FRAME_DT
        fa = time.time()
        fields = state.fields_at(t)
        fields["S"] = state.S
        fb = time.time()
        E, info = photon_frame(fields, state.S, args.grid, env,
                               args.photons,
                               seed=(state.seed + 1000 + i) % (2 ** 63))
        t_transport += time.time() - fb
        if i == 0:
            E0 = E.copy()
            transport_info0 = info
        fc = time.time()
        W, H = (PW, PH) if i == 0 else (VW, VH)
        img = render_frame(fields, E, photo_lin, env, W, H)
        t_image += time.time() - fc
        frames_u8.append(img)
        del fields, E
        print("frame %d/%d t=%.2fs (fields %.1fs)" % (i + 1, args.frames, t, fb - fa),
              flush=True)

    # still (native photo size, t=0 frame)
    still = frames_u8[0]
    Image.fromarray(still).save(os.path.join(args.out, "preset1-still.png"))

    # video (even dims, 20 fps)
    import cv2
    vpath = os.path.join(args.out, "preset1-test16.mp4")
    vw = cv2.VideoWriter(vpath, cv2.VideoWriter_fourcc(*"mp4v"),
                         VIDEO_FPS, (VW, VH))
    if not vw.isOpened():
        raise SystemExit("VideoWriter failed to open")
    for im in frames_u8:
        if im.shape[1] != VW or im.shape[0] != VH:
            im2 = np.asarray(Image.fromarray(im).resize((VW, VH),
                                                        Image.LANCZOS))
        else:
            im2 = im
        vw.write(cv2.cvtColor(im2, cv2.COLOR_RGB2BGR))
    vw.release()

    # field stats (linear, pre-tonemap) on still frame buffer
    field_stats = lum_stats(to_lum(E0))

    # ref1 water-region crop stats (linear luminance)
    ref = np.asarray(Image.open(args.ref1).convert("RGB"),
                     dtype=np.float64) / 255.0
    RH, RW, _ = ref.shape
    CW, CH = 512, 512
    cx0, cy0 = (RW - CW) // 2, (RH - CH) // 2
    crop = ref[cy0:cy0 + CH, cx0:cx0 + CW, :]
    ref1_stats = lum_stats(to_lum(srgb_to_linear(crop)))
    ref1_stats["file"] = os.path.basename(args.ref1)
    ref1_stats["crop_xywh"] = [cx0, cy0, CW, CH]

    numbers = {
        "recipe": recipe.get("name"),
        "recipe_schema": recipe.get("schema"),
        "preset": args.preset,
        "seed": int(state.seed),
        "state": {
            "grid": [args.grid, args.grid],
            "patch_S_m": state.S,
            "dx_m": state.dx,
            "kmin_rad_per_m": float(state.kmin),
            "kmax_nyquist_rad_per_m": float(state.kmax_nyq),
            "band_decades": float(state.decades),
            "specL_m": float(state.specL),
            "choppy_lambda": 0.85,
            "omega0_rad_per_s": float(2.0 * np.pi / 8.0),
            "loop_T_s": 8.0,
            "target_std_m": float(state.target_std),
            "amplitude_scale": float(state.amplitude_scale),
            "damping_l_m": 0.005,
        },
        "transport": {
            "photons_per_frame": args.photons,
            "frames": args.frames,
            "frame_dt_s": FRAME_DT,
            "sun_disc_rad": float(__import__("transport").sun_disc_radius(
                env.get("clarity", 0.8))),
            "buffer": [args.grid, args.grid],
            "sigma_px": 0.75,
            "ior_rgb": [1.33132, 1.33434, 1.33925],
            "expected_mean_rgb": transport_info0["expected_mean"],
            "measured_mean_rgb_frame0": transport_info0["measured_mean"],
        },
        "field_stats_frame0_linear": field_stats,
        "ref1_stats_linear": ref1_stats,
        "loop_check": {"t0_s": 0.0, "t1_s": 8.0,
                       "max_abs_diff_h_m": float(loop_diff)},
        "timing_s": {
            "state_build": t_state,
            "transport_total": t_transport,
            "image_total": t_image,
            "wall_total": time.time() - t_all,
        },
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"],
            "depthM": env["depthM"],
        },
        "outputs": {
            "still": "preset1-still.png",
            "video": "preset1-test16.mp4",
            "video_fps": VIDEO_FPS,
            "still_size": [int(PW), int(PH)],
            "video_size": [int(VW), int(VH)],
        },
    }
    with open(os.path.join(args.out, "preset1-numbers.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)
    print("OK wall=%.1fs" % (time.time() - t_all), flush=True)


if __name__ == "__main__":
    main()
