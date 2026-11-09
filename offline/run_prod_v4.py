"""Pipeline v4 production loop render (PIPELINE-V4, all parameters frozen).

Frozen: preset 1 (recipe-1, deterministic seed), state 2048^2 Phillips FFT
(omega quantized to 2pi/8s multiples), transport det.py Q=25 sigma 0.75px
(float32 master, unclamped), composition = render_frame_v3 A=2.6 +
displacement lowpass sigma 3.0px@2048^2, plus C-layer exposure E=2.2
(direct term only, same formula as run_v3_exposure.py).

Frames: 24fps x 192 (8 s loop), t = frame/24 s.
Idempotent per frame: if out/C-XXXX.npy + flow-XXXX.npy + frame-XXXX.png
all exist, the frame is skipped (logged as skip).
Progress appended per frame to out/progress.log.

Usage (from caustics-filter/offline/):
  python run_prod_v4.py --recipe ../blender/recipe-1.json --preset 1
All logging goes to files (progress.log, prod-run.log); stdout is minimal
so the script survives stdio-restricted sandboxes.
"""
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from state import WaterState, sample_bicubic  # noqa: E402
from det import det_frame, flow_backward_master  # noqa: E402
from image import (render_frame_v3, srgb_to_linear, linear_to_srgb,  # noqa: E402
                   aces_narkowicz, e_base_dir_from_env, IOR)

FPS = 24
NFRAMES = 192
GRID = 2048
Q_FROZEN = 25
A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2
SEED_V2 = 3255045276382562994

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)


def frame_files(out, i):
    tag = "%04d" % i
    return (os.path.join(out, "C-%s.npy" % tag),
            os.path.join(out, "flow-%s.npy" % tag),
            os.path.join(out, "frame-%s.png" % tag))


def to_lum(rgb):
    return (np.asarray(rgb, dtype=np.float64) * REC709).sum(axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--out", default="out/prod")
    ap.add_argument("--v2dir", default="out/v2")
    ap.add_argument("--photo", default="../samples/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--frames", default="0-191")
    args = ap.parse_args()

    if "-" in str(args.frames):
        lo, hi = str(args.frames).split("-")
        frange = range(int(lo), int(hi) + 1)
    else:
        frange = [int(x) for x in str(args.frames).split(",") if x.strip()]

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    runlog = os.path.join(args.out, "prod-run.log")
    prolog = os.path.join(args.out, "progress.log")

    def rlog(msg):
        with open(runlog, "a", encoding="utf-8") as f:
            f.write(msg + "\n")

    def prog(msg):
        with open(prolog, "a", encoding="utf-8") as f:
            f.write(msg + "\n")

    rlog("v4 run start %s frames=%s"
         % (datetime.datetime.now(datetime.timezone.utc).isoformat(),
            args.frames))

    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    if recipe.get("schema") != "caustics-recipe/v1":
        rlog("ERROR recipe schema mismatch")
        raise SystemExit("recipe schema mismatch")
    env = resolve_environment(recipe)

    t0 = time.time()
    state = WaterState(recipe, env, N=GRID, seed=args.seed,
                       preset_index=args.preset)
    t_state = time.time() - t0
    if int(state.seed) != SEED_V2:
        rlog("ERROR seed %d != v2 seed %d" % (int(state.seed), SEED_V2))
        raise SystemExit("seed drift vs v2")
    rlog("state built N=%d S=%.10f seed=%d (%.2fs)"
         % (GRID, float(state.S), int(state.seed), t_state))

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape
    rlog("photo %dx%d" % (PW, PH))

    # ---- frozen composition constants (same as run_v3_exposure.py) ----
    S = float(state.S)
    N_buf = GRID
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    Ebase = e_base_dir_from_env(env)
    absorb = np.asarray(env["absorb"], dtype=np.float64)
    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    beer = np.exp(-absorb * path_len * 0.22)
    Cc_dtype_note = "float32"
    rlog("Ebase=%s beer=%s strength=%.4g exposure=%.4g" %
         (list(float(x) for x in Ebase), list(float(x) for x in beer),
          strength, exposure))

    H, W = PH, PW
    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S

    results = []
    cum = 0.0
    for i in frange:
        c_path, f_path, p_path = frame_files(args.out, i)
        if (os.path.exists(c_path) and os.path.exists(f_path)
                and os.path.exists(p_path)):
            prog("frame=%04d t=%.4fs status=skip (artifacts exist)"
                 % (i, i / float(FPS)))
            results.append({"frame": i, "t_s": i / float(FPS),
                            "status": "skip"})
            continue
        tv = time.time()
        t = i / float(FPS)
        fields = state.fields_at(t)
        fields["S"] = state.S
        C, info = det_frame(fields, S, N_buf, env, Q_FROZEN)
        flow = flow_backward_master(fields, S, N_buf, env)
        np.save(c_path, C.astype(np.float32))
        np.save(f_path, flow.astype(np.float32))
        still, dbg = render_frame_v3(fields, C, photo_lin, env, PW, PH,
                                     ambient_scale=A_FROZEN,
                                     lowpass_sigma_px=B_SIGMA_FROZEN,
                                     return_debug=True)
        lin_base = dbg["lin"]
        # ---- C-layer-only direct component (identical to exposure probe)
        n_state = np.asarray(fields["gx"]).shape[0]
        sig_eff = B_SIGMA_FROZEN * (float(n_state) / 2048.0)
        from scipy import ndimage as _ndi
        gx_s = _ndi.gaussian_filter(
            np.asarray(fields["gx"], dtype=np.float64), sig_eff, mode="wrap")
        gy_s = _ndi.gaussian_filter(
            np.asarray(fields["gy"], dtype=np.float64), sig_eff, mode="wrap")
        gx = sample_bicubic(gx_s.astype(np.float32), px, py,
                            S).astype(np.float64)
        gy = sample_bicubic(gy_s.astype(np.float32), px, py,
                            S).astype(np.float64)
        nl = np.sqrt(gx * gx + gy * gy + 1.0)
        nx = -gx / nl
        ny = -gy / nl
        Cc = np.asarray(C, dtype=np.float64)
        direct_lin = np.zeros((H, W, 3), dtype=np.float64)
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
        alb_re = np.zeros_like(dbg["alb"])
        for b in range(3):
            Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
            shx = -nx * depth * Cb
            shy = -ny * depth * Cb
            su = np.clip(u + shx / S, 0.0, 1.0) * (PW - 1.0)
            sv = np.clip(v + shy / S, 0.0, 1.0) * (PH - 1.0)
            x0 = np.clip(np.floor(su).astype(np.int64), 0, PW - 1)
            y0 = np.clip(np.floor(sv).astype(np.int64), 0, PH - 1)
            x1 = np.clip(x0 + 1, 0, PW - 1)
            y1 = np.clip(y0 + 1, 0, PH - 1)
            fx = su - np.floor(su)
            fy = sv - np.floor(sv)
            alb_re[:, :, b] = (photo_lin[y0, x0, b] * (1 - fx) * (1 - fy)
                               + photo_lin[y0, x1, b] * fx * (1 - fy)
                               + photo_lin[y1, x0, b] * (1 - fx) * fy
                               + photo_lin[y1, x1, b] * fx * fy)
        # direct_lin alb consistency vs render_frame_v3 debug alb:
        # recover alb from direct_lin via direct map is exact by
        # construction (same loop); check dbg alb vs alb_re instead.
        maxdiff_alb = float(np.max(np.abs(alb_re - dbg["alb"])))
        lin_E = lin_base + (E_FROZEN - 1.0) * direct_lin
        aces_E = aces_narkowicz(np.maximum(lin_E, 0.0))
        out = (np.clip(linear_to_srgb(aces_E), 0.0, 1.0) * 255.0
               + 0.5).astype(np.uint8)
        Image.fromarray(out).save(p_path)
        lin_lum = to_lum(lin_E)
        aces_lum = to_lum(np.clip(aces_E, 0.0, 1.0))
        dt = time.time() - tv
        cum += dt
        C64 = np.asarray(C, dtype=np.float64)
        rec = {"frame": i, "t_s": t, "status": "rendered",
               "C_mean": float(C64.mean()), "C_max": float(C64.max()),
               "C_min": float(C64.min()),
               "lin_mean_post_E": float(np.mean(lin_E)),
               "crush": float(np.mean(lin_lum < 0.02)),
               "blow": float(np.mean(aces_lum > 0.98)),
               "alb_maxdiff_vs_v3dbg": maxdiff_alb,
               "time_s": dt}
        results.append(rec)
        prog("frame=%04d t=%.4fs status=rendered dt=%.1fs cum=%.0fs "
             "Cmean=%.6f Cmax=%.3f crush=%.4f blow=%.6f albDiff=%.2g" %
             (i, t, dt, cum, rec["C_mean"], rec["C_max"], rec["crush"],
              rec["blow"], maxdiff_alb))
        rlog("frame %04d done dt=%.1f cum=%.0f" % (i, dt, cum))
        del fields, C, flow, still, dbg, lin_base, direct_lin, lin_E

    # ---- seam verification: state t=8s vs t=0 (C max diff) ----
    rlog("seam check start")
    f0 = state.fields_at(0.0)
    f0["S"] = state.S
    C0_chk, _ = det_frame(f0, S, N_buf, env, Q_FROZEN)
    fl0 = flow_backward_master(f0, S, N_buf, env)
    f8 = state.fields_at(8.0)
    f8["S"] = state.S
    C8, _ = det_frame(f8, S, N_buf, env, Q_FROZEN)
    fl8 = flow_backward_master(f8, S, N_buf, env)
    C0_master = np.load(frame_files(args.out, 0)[0]).astype(np.float64)
    seam_C = float(np.max(np.abs(C8.astype(np.float64)
                                 - C0_master)))
    seam_C_self = float(np.max(np.abs(C8.astype(np.float64)
                                      - C0_chk.astype(np.float64))))
    seam_flow = float(np.max(np.abs(fl8.astype(np.float64)
                                    - fl0.astype(np.float64))))
    rlog("seam C(t=8 vs master C-0000) maxdiff=%.6g" % seam_C)
    rlog("seam C(t=8 vs fresh t=0) maxdiff=%.6g" % seam_C_self)
    rlog("seam flow(t=8 vs t=0) maxdiff=%.6g px" % seam_flow)

    # ---- representative still: frame 96 ----
    still_src = frame_files(args.out, 96)[2]
    still_dst = os.path.join(args.out, "still-f96.png")
    if os.path.exists(still_src):
        shutil.copyfile(still_src, still_dst)
        rlog("still copied frame-0096 -> still-f96.png")
    else:
        rlog("ERROR still source missing: frame 96")
        raise SystemExit("frame 96 missing")

    numbers = {
        "pipeline": "v4",
        "recipe": recipe.get("name"),
        "recipe_schema": recipe.get("schema"),
        "preset": args.preset,
        "seed": int(state.seed),
        "frames": {"fps": FPS, "count": NFRAMES,
                   "t_s": [i / float(FPS) for i in frange],
                   "loop_T_s": 8.0},
        "frozen": {
            "state": {"grid": [GRID, GRID], "patch_S_m": float(state.S),
                      "dx_m": float(state.dx)},
            "transport": {"method": "det.py det_frame",
                          "Q": Q_FROZEN, "pixel_filter_sigma_px": 0.75,
                          "master_dtype": Cc_dtype_note,
                          "processing": "none (no clamp/renorm)"},
            "composition": {
                "base": "render_frame_v3",
                "A_ambient_scale": A_FROZEN,
                "B_band": "low",
                "lowpass_sigma_px_state": B_SIGMA_FROZEN,
                "C_layer_exposure": E_FROZEN,
                "formula": ("lin_E = (alb*(E*direct + E_amb_A2.6)*strength"
                            " + sky_A2.6) * exposure_env, "
                            "E=2.2 on the C-modulated direct term only"),
            },
        },
        "seam": {
            "C_t8_vs_C_frame0000_maxdiff": seam_C,
            "C_t8_vs_fresh_t0_maxdiff": seam_C_self,
            "flow_t8_vs_t0_maxdiff_px": seam_flow,
        },
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "depthM": env["depthM"], "clarity": env.get("clarity"),
        },
        "photo_size": [int(PW), int(PH)],
        "padspec_video": {
            "note": ("photo-native 735x539 has odd dims; right+bottom "
                     "edge-replicated 1px pad to 736x540 for H.264 encode"),
        },
        "timing_s": {"state_build": t_state,
                     "wall_render_total": time.time() - t_all},
        "outputs": {
            "C": ["C-%04d.npy" % i for i in frange],
            "flow": ["flow-%04d.npy" % i for i in frange],
            "frames": ["frame-%04d.png" % i for i in frange],
            "still": "still-f96.png",
        },
        "per_frame": results,
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
    }
    with open(os.path.join(args.out, "numbers-prod.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)
    rlog("numbers-prod.json written frames=%d wall=%.0fs"
         % (len(results), time.time() - t_all))
    prog("run-phase1(render+seam+still+numbers) done wall=%.0fs"
         % (time.time() - t_all))


if __name__ == "__main__":
    main()
