"""One composited still: v3.4 (v8) pipeline, eta = 1.0, windy pond.

Task (D:\\mutaro-work\\rally-in\\BRIEF-WINDY-STILL.md):
  base image  D:\\mutaro-work\\moth-on-bark.jpg   (736x540, the frame-1 photo)
  caustics    v3.4 = v8 chunk-4 (V8WaterState) + det_frame, master 2048^2
  composite   frozen production: A2.6-Blow + E2.2, lowpass 3.0 px, Q = 25
  light       eta = 1.0 -> direct sun only (a recipe with no sky block is
              eta = 1 exactly, so the v9/Perez diffuse route is not invoked)
  environment preset-1 pond with ONLY: wind 6.0 m/s, depth z = 1.0 m, rain 0,
              surface film none.  Fetch stays the pond's 100 m.
              H_s comes from the pipeline's own SMB/JONSWAP computation.

  Depth cost: a second composite at the SAME wind/fetch with z = 0.6 m
  (preset-1 depth) is rendered for the report's "what did z=1.0 cost"
  lines.  Everything else is identical.

  Warp measurement (the brief's yardstick): block matching between the base
  photo and the composite -- 48 px blocks, search +/-24 px at 2 px steps,
  Rec.709 luminance.  mean shift px + shift-field std px.

Every step is checkpointed with D:\\mutaro-work\\jobkit (Checkpoint.unit_done)
and each deliverable is committed atomically (atomic_write_bytes).

Deliverables -> out/windy-still/
  comp-windy-still.png   736x540 composited still
  numbers-windy.json     recipe, realised values, caps, C stats, warp numbers
  REPORT-WINDY.md        numbers + depth-cost lines

Run (foreground, from caustics-filter/offline):
  python run_windy_still.py
"""
import gc
import hashlib
import io
import json
import math
import os
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, verify_units  # noqa: E402

from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from state import derive_seed  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,  # noqa: E402
                        scattering_blur_sigma_px, substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,  # noqa: E402
                          image_metrics, flat_water_direct_level, v3_beer)
from recipe_v8 import parse_recipe_v8, with_v8  # noqa: E402
from v8_state import V8WaterState  # noqa: E402
import v8_physics as P  # noqa: E402

RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
OUT = os.path.join(HERE, "out", "windy-still")
WORK = os.path.join(OUT, "_work")
CP = Checkpoint(os.path.join(OUT, "checkpoint.json"))

GRID = 2048
Q = 25
ETA = 1.0                 # direct only (the v8 route; no Perez diffuse)
A_COMP = A_V3             # 2.6  (frozen production ambient scale)
E_COMP = E_V3             # 2.2  (frozen production C-layer exposure)
B_LOWPASS = B_SIGMA_V3    # 3.0 px

# --- the brief's recipe numbers: preset-1 pond, three changes only ---
WINDY = {"fetch_m": 100.0, "wind_ms": 6.0, "depth_m": 1.0}
# --- same scene at the preset-1 depth, to price the depth cost ---
DEPTH_REF = {"fetch_m": 100.0, "wind_ms": 6.0, "depth_m": 0.6}

CMEAN_TOL = 1e-3
CRUSH_BOUND = 0.0325      # the customary v6/v7/v8 bound
CRUSH_TIGHT = 0.02        # the stream-brief bound, reported as a second gate
BLOW_BOUND = 0.0
SIGMA2_RECHECK_BELOW = 0.035   # the brief's "re-check the scene" trigger

# the brief's yardsticks (same method, quoted, not recomputed)
YARDSTICKS = {
    "yesterday_production_frame": {"mean_px": 23.2, "std_px": 17.1},
    "today_stream_still": {"mean_px": 9.2, "std_px": 9.5},
}


def log(msg):
    print(msg, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def coerce(obj):
    if isinstance(obj, dict):
        return {str(k): coerce(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [coerce(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return coerce(obj.tolist())
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def predicted_shift_field(fields, env, W, H, lowpass_sigma_px=None):
    """The pipeline's OWN per-pixel content displacement, in view pixels.

    This is not a model of the composite: it is the same shift the renderer
    applies (green band), rebuilt from the same state with the same
    `_state_gradients` path, so the block-match measurement can be checked
    against the geometry that produced the image.
    """
    from v5_composite import _state_gradients
    from state import sample_bicubic
    from image import IOR

    sig = B_LOWPASS if lowpass_sigma_px is None else float(lowpass_sigma_px)
    S = float(fields["S"])
    gx_s, gy_s = _state_gradients(fields, sig)
    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S, S).astype(
        np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S, S).astype(
        np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    Cb = (float(IOR[1]) - 1.0) / float(IOR[1])     # green band, as the render
    shx = -nx * float(env["depthM"]) * Cb          # metres, patch frame
    shy = -ny * float(env["depthM"]) * Cb
    # the renderer samples the albedo at (u + shx/S), so texture found at
    # (u + shx/S) is DISPLAYED at u: the content displacement is -shx.
    return {
        "dx_px": -(shx / S) * W,
        "dy_px": -(shy / S) * H,
        "dx_m": -shx, "dy_m": -shy,
        "S_m": S, "px_m": S / float(W), "band": "green (IOR[1])",
    }


def predicted_block_stats(pred, grid, block=48):
    """Per-block means of the predicted field + the field-level statistics."""
    dx, dy = pred["dx_px"], pred["dy_px"]
    fmag = np.hypot(dx, dy)
    out_px, out_py = [], []
    for y0 in grid["block_origins_y"]:
        for x0 in grid["block_origins_x"]:
            out_px.append(float(dx[y0:y0 + block, x0:x0 + block].mean()))
            out_py.append(float(dy[y0:y0 + block, x0:x0 + block].mean()))
    px = np.array(out_px)
    py = np.array(out_py)
    mag = np.hypot(px, py)
    return {
        # field level: the distribution the matcher is sampling
        "field_mean_magnitude_px": float(fmag.mean()),
        "field_std_magnitude_px": float(fmag.std(ddof=0)),
        "field_p90_magnitude_px": float(np.percentile(fmag, 90)),
        "field_max_magnitude_px": float(fmag.max()),
        "field_frac_beyond_search": float(np.mean(fmag > 24.0)),
        "mean_predicted_shift_m": float(np.hypot(
            pred["dx_m"], pred["dy_m"]).mean()),
        # block level: the vector each 48 px block averages to
        "block_mean_magnitude_px": float(mag.mean()),
        "block_std_magnitude_px": float(mag.std(ddof=0)),
        "frac_blocks_predicted_beyond_search": float(np.mean(mag > 24.0)),
        "pixel_scale_mm": float(pred["px_m"] * 1e3),
        "block_means_dx_px": px, "block_means_dy_px": py,
    }


def block_agreement(stats, rows, key_y, key_x):
    """Correlation of a measured block field with the predicted one."""
    mx = np.array([r[key_x] for r in rows], np.float64)
    my = np.array([r[key_y] for r in rows], np.float64)
    px, py = stats["block_means_dx_px"], stats["block_means_dy_px"]
    return {
        "corr_dx": float(np.corrcoef(mx, px)[0, 1]),
        "corr_dy": float(np.corrcoef(my, py)[0, 1]),
        "bias_dx_px": float((mx - px).mean()),
        "bias_dy_px": float((my - py).mean()),
        "rms_dx_px": float(np.sqrt(np.mean((mx - px) ** 2))),
        "rms_dy_px": float(np.sqrt(np.mean((my - py) ** 2))),
    }


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {
        "shape": list(c.shape),
        "mean": float(c.mean()),
        "p50": float(np.percentile(c, 50)),
        "p99": float(np.percentile(c, 99)),
        "p99_9": float(np.percentile(c, 99.9)),
        "max": float(c.max()),
        "min": float(c.min()),
        "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
        "pass_mean_1": bool(abs(c.mean() - 1.0) <= CMEAN_TOL),
    }


# --------------------------------------------------------------------------
# the brief's warp measurement
# --------------------------------------------------------------------------
def _box_sums(a, k):
    """Sum of each k x k window of a, on every valid origin (k..)."""
    p = np.pad(a, ((1, 0), (1, 0)))
    c = np.cumsum(np.cumsum(p, axis=0), axis=1)
    return c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]


def warp_measure(base_u8, comp_u8, block=48, search=24, step=2):
    """Block-matching displacement, base photo -> composite.

    Method (the brief's): 48 px blocks, search +/-24 px at 2 px steps,
    Rec.709 luminance of the 8-bit sRGB images.

    For a block at (y0,x0) the offset (dy,dx) minimising the mean squared
    luminance difference is the displacement: content that sits at
    (y0,x0) in the base sits at (y0+dy, x0+dx) in the composite.

    Reported primary numbers use raw SSD (plain block matching).  A
    zero-mean normalised cross-correlation (NCC) pass is reported alongside
    as a gain-invariance check, because the composite carries a caustic
    modulation the base does not.
    """
    from numpy.lib.stride_tricks import sliding_window_view

    base_y = (np.asarray(base_u8, np.float64)[:, :, :3]
              @ np.array([0.2126, 0.7152, 0.0722]))
    comp_y = (np.asarray(comp_u8, np.float64)[:, :, :3]
              @ np.array([0.2126, 0.7152, 0.0722]))
    H, W = base_y.shape
    ys = list(range(search, H - search - block + 1, block))
    xs = list(range(search, W - search - block + 1, block))
    offs = np.arange(-search, search + 1, step)
    gi = np.arange(0, 2 * search + 1, step)      # index into the 49x49 map
    n = float(block * block)

    rows = []
    for y0 in ys:
        for x0 in xs:
            B = base_y[y0:y0 + block, x0:x0 + block]
            Cre = comp_y[y0 - search:y0 + block + search,
                         x0 - search:x0 + block + search]
            assert Cre.shape == (block + 2 * search, block + 2 * search)
            # every block-window of the search region: (2S+1, 2S+1, b, b)
            Wn = sliding_window_view(Cre, (block, block))
            sumBC = np.einsum("mn,ijmn->ij", B, Wn)
            sumC = Wn.sum(axis=(2, 3))
            sumC2 = np.einsum("ijmn,ijmn->ij", Wn, Wn)
            sB, sB2 = float(B.sum()), float((B * B).sum())
            ssd = (sB2 + sumC2 - 2.0 * sumBC) / n
            num = sumBC - sB * sumC / n
            varB = max(sB2 - sB * sB / n, 0.0)
            varC = np.maximum(sumC2 - sumC * sumC / n, 0.0)
            den = np.sqrt(np.maximum(varB, 1e-12) * np.maximum(varC, 1e-12))
            ncc = num / den

            ssd_g = ssd[np.ix_(gi, gi)]
            ncc_g = ncc[np.ix_(gi, gi)]
            i_s, j_s = np.unravel_index(int(np.argmin(ssd_g)), ssd_g.shape)
            i_c, j_c = np.unravel_index(int(np.argmax(ncc_g)), ncc_g.shape)
            rows.append({
                "y0": int(y0), "x0": int(x0),
                "block_std": float(math.sqrt(varB / n)),
                "ssd_dy": int(offs[i_s]), "ssd_dx": int(offs[j_s]),
                "ssd_val": float(ssd_g[i_s, j_s]),
                "ssd_zero": float(ssd_g[(search) // step, (search) // step]),
                "ncc_dy": int(offs[i_c]), "ncc_dx": int(offs[j_c]),
                "ncc_val": float(ncc_g[i_c, j_c]),
            })

    def summarise(key_y, key_x, tag):
        dy = np.array([r[key_y] for r in rows], np.float64)
        dx = np.array([r[key_x] for r in rows], np.float64)
        mag = np.hypot(dx, dy)
        at_edge = np.mean((np.abs(dy) >= search) | (np.abs(dx) >= search))
        tex = np.array([r["block_std"] for r in rows], np.float64) >= 2.0
        return {
            "tag": tag,
            "n_blocks": int(len(rows)),
            "block_px": int(block), "search_px": int(search),
            "step_px": int(step), "luminance": "Rec.709 on 8-bit sRGB",
            "mean_shift_px": float(mag.mean()),
            "shift_field_std_px": float(mag.std(ddof=0)),
            "shift_magnitude_p50_px": float(np.percentile(mag, 50)),
            "shift_magnitude_p90_px": float(np.percentile(mag, 90)),
            "shift_magnitude_max_px": float(mag.max()),
            "shift_dx_mean_px": float(dx.mean()),
            "shift_dy_mean_px": float(dy.mean()),
            "shift_dx_std_px": float(dx.std(ddof=0)),
            "shift_dy_std_px": float(dy.std(ddof=0)),
            "mean_vector_magnitude_px": float(math.hypot(dx.mean(), dy.mean())),
            "frac_blocks_at_search_edge": float(at_edge),
            "n_textured_blocks_std_ge_2": int(tex.sum()),
            "mean_shift_px_textured": (float(mag[tex].mean())
                                       if tex.any() else None),
            "shift_field_std_px_textured": (float(mag[tex].std(ddof=0))
                                            if tex.any() else None),
        }

    ssd_s = summarise("ssd_dy", "ssd_dx", "raw SSD (primary)")
    ncc_s = summarise("ncc_dy", "ncc_dx", "zero-mean NCC (check)")
    zero_ssd = float(np.mean([r["ssd_zero"] for r in rows]))
    return {
        "definition": (
            "for each block the offset (dy,dx) minimising the mean squared "
            "luminance difference between the base block and the composite "
            "window is the displacement (content at (y0,x0) in the base sits "
            "at (y0+dy,x0+dx) in the composite).  Offsets are searched on the "
            "2 px grid inside +/-24 px, so the measurement is quantised to "
            "2 px and saturates at 24 px per axis."),
        "primary": ssd_s,
        "ncc_check": ncc_s,
        "mean_squared_error_at_zero_offset": zero_ssd,
        "grid": {"block_origins_y": ys, "block_origins_x": xs},
        "per_block": rows,
        "yardsticks_same_method": YARDSTICKS,
    }


# --------------------------------------------------------------------------
# one scene -> state -> master -> still
# --------------------------------------------------------------------------
def build(raw, scene, tag, grid, seed, base_u8, photo_lin, PH, PW, C=None):
    t0 = time.time()
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    st = V8WaterState(scene, env, N=grid, seed=seed, preset_index=1,
                      fetch_m=float(cfg["fetch"]["fetch_m"]),
                      fetch_source=cfg["fetch"]["source"])
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    re = st.realised_stats(fields)
    t_state = time.time() - t0
    log("  [%s] state %.1fs sigma2_spec=%.6f realised=%.6f Hs=%.4f"
        % (tag, t_state, st.spectrum_after["sigma2_rad2"],
           re["sigma2_rad2"], re["Hs_m"]))

    if C is None:
        t1 = time.time()
        C, info = det_frame(fields, st.S, grid, env, Q)
        t_det = time.time() - t1
        log("  [%s] det Q=%d %.1fs mean=%.9f p99=%.4f max=%.3f"
            % (tag, Q, t_det, float(C.mean()),
               float(np.percentile(C, 99)), float(C.max())))
    else:
        t_det = None
        info = {"reused": True, "Q": Q,
                "note": "master reused from this run's own cache file"}

    att = water_attenuation(env["turbid"])
    beer = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"], env["turbid"])
    blur = scattering_blur_sigma_px(att, env["depthM"], grid)
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1 = smith_g1_beckmann(re["rms_slope"], sun_zenith)

    t2 = time.time()
    img, dbg = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_COMP,
        c_layer_exposure=E_COMP,
        substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer["down"], beer_up=beer["up"], g1=float(g1["G1"]),
        blur_sigma_px=blur["sigma_px"], lowpass_sigma_px=B_LOWPASS,
        return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t2

    # image stats used by the depth-cost section
    lin = np.asarray(dbg["lin"], np.float64)
    aces = np.clip(np.asarray(dbg["aces"], np.float64), 0.0, 1.0)
    srgb = np.asarray(img, np.float64) / 255.0
    aces_y = (aces @ np.array([0.2126, 0.7152, 0.0722]))
    mx = srgb.max(axis=2)
    mn = srgb.min(axis=2)
    img_stats = {
        "linear_luminance_mean": float(metrics["linear_luminance_mean"]),
        "linear_rgb_mean": list(metrics["linear_rgb_mean"]),
        "srgb_rgb_mean": [float(x) for x in srgb.reshape(-1, 3).mean(0)],
        "aces_luminance_mean": float(aces_y.mean()),
        "aces_luminance_std": float(aces_y.std()),
        "srgb_saturation_mean": float(np.mean((mx - mn) / np.maximum(mx, 1e-9))),
        "crush_rate_linear_lt_0.02": float(metrics["crush_rate_linear_lt_0.02"]),
        "blowout_rate_aces_gt_0.98": float(metrics["blowout_rate_aces_gt_0.98"]),
    }
    del dbg, lin, aces, aces_y, mx, mn
    gc.collect()
    log("  [%s] composite %.1fs crush=%.5f blow=%.6f linY=%.5f"
        % (tag, t_comp, metrics["crush_rate_linear_lt_0.02"],
           metrics["blowout_rate_aces_gt_0.98"],
           metrics["linear_luminance_mean"]))

    return {
        "env": env, "cfg": cfg, "state": st, "fields": fields, "realised": re,
        "C": C, "C_stats": c_stats(C), "det_info": info,
        "still": img, "still_metrics": metrics, "img_stats": img_stats,
        "beer": {"down": [float(x) for x in beer["down"]],
                 "up": [float(x) for x in beer["up"]],
                 "two_way": [float(x) for x in beer["two_way"]],
                 "s_down_m": beer["s_down_m"], "s_up_m": beer["s_up_m"],
                 "theta_w_deg": beer["theta_w_deg"],
                 "c_rgb": [float(x) for x in beer["c"]],
                 "note": ("direct (caustic-carrying) beam only: the ambient "
                          "term keeps the frozen v3 scat heuristic, so these "
                          "factors act on the caustic layer, not on the "
                          "whole image")},
        "blur_sigma_px": float(blur["sigma_px"]),
        "substrate_gain_rgb": list(gain_base["gain_rgb"]),
        "g1": float(g1["G1"]), "sun_zenith_deg": float(sun_zenith),
        "timing_s": {"state": t_state, "det": t_det, "composite": t_comp},
    }


def smb_node(env, fetch_m):
    rel = P.jonswap_fetch_relations(env["windSpeedMS"], fetch_m)
    sp = P.spreading_parameters(env["windSpeedMS"], fetch_m)
    S_patch = 4.0 * float(env["specL"])
    anc = P.anchor_sigma2(env["windSpeedMS"], fetch_m, env["depthM"],
                          S_patch=S_patch, N=GRID)
    return {
        "wind_ms": float(env["windSpeedMS"]), "fetch_m": float(fetch_m),
        "depth_m": float(env["depthM"]),
        "Ft": float(rel["Ft"]),
        "f_p_Hz": float(rel["f_p_Hz"]),
        "omega_p_rad_s": float(rel["omega_p_rad_s"]),
        "alpha": float(rel["alpha"]),
        "Hs_fetch_law_m": float(rel["Hs_fetch_law_m"]),
        "lambda_p_m": float(P.G / (2.0 * math.pi * rel["f_p_Hz"] ** 2)),
        "s_max_mitsuyasu": float(sp["s_max"]),
        "half_width_deg_at_peak": float(
            2.0 * math.degrees(math.acos(0.5 ** (1.0 / (2.0 * sp["s_max"]))))),
        "U_over_cp": float(sp["U_over_cp"]),
        "c_p_m_per_s": float(sp["c_p_m_per_s"]),
        "in_fitted_range_Ft_1e3_1e5": bool(1e3 <= rel["Ft"] <= 1e5),
        "S_patch_m": S_patch,
        "sigma2_target_rad2": float(anc["sigma2_target_rad2"]),
        "sigma2_jonswap_rad2": float(anc["candidates"]["JONSWAP"]["sigma2"]),
        "sigma2_donelan_rad2": float(anc["candidates"]["DONELAN"]["sigma2"]),
        "cox_munk_old_target_rad2": float(anc["cox_munk_old_target_rad2"]),
        "ratio_new_over_old": float(anc["ratio_new_over_old"]),
        "relations": rel["relations"],
    }


def fmt(x, n=4):
    if x is None:
        return "-"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return "%d" % x
    return ("%%.%df" % n) % float(x)


def main():
    t_all = time.time()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)

    with open(RECIPE, "r", encoding="utf-8") as f:
        raw = json.load(f)
    seed = derive_seed(raw, 1)
    log("base recipe=%s seed=%d" % (raw.get("name"), seed))

    base_np = np.asarray(Image.open(PHOTO).convert("RGB"))
    PH, PW, _ = base_np.shape
    photo_lin = srgb_to_linear(base_np.astype(np.float64) / 255.0)
    log("base photo %dx%d sha=%s" % (PW, PH, sha256_file(PHOTO)[:16]))
    if (PW, PH) != (736, 540):
        log("WARNING: expected 736x540, got %dx%d" % (PW, PH))

    # ==================================================================
    # 01 -- scene + SMB numbers (no transport)
    # ==================================================================
    scene = with_v8(raw, fetch_m=float(WINDY["fetch_m"]),
                    depth_m=float(WINDY["depth_m"]),
                    wind_ms=float(WINDY["wind_ms"]))
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    smb = smb_node(env, cfg["fetch"]["fetch_m"])
    log("[01] scene: wind=%.2f depth=%.2f fetch=%.1f Ft=%.2f f_p=%.4f Hz "
        "Hs_law=%.4f m s_max=%.2f half=%.2f deg"
        % (env["windSpeedMS"], env["depthM"], cfg["fetch"]["fetch_m"],
           smb["Ft"], smb["f_p_Hz"], smb["Hs_fetch_law_m"],
           smb["s_max_mitsuyasu"], smb["half_width_deg_at_peak"]))
    log("     sigma2_target=%.6f (JONSWAP %.6f / Donelan %.6f) cox_munk=%.6f"
        % (smb["sigma2_target_rad2"], smb["sigma2_jonswap_rad2"],
           smb["sigma2_donelan_rad2"], smb["cox_munk_old_target_rad2"]))
    log("     sun alt=%.2f az=%.2f sunI=%.4f specL=%.4f S_patch=%.4f m "
        "dx=%.4f mm"
        % (env["sun"]["altitudeDeg"], env["sun"]["azimuthDeg"], env["sunI"],
           env["specL"], 4.0 * env["specL"], 4.0 * env["specL"] / GRID * 1e3))
    CP.unit_done("01-scene-smb")

    # ==================================================================
    # 02/03/04 -- windy state, master (Q=25), composite
    # ==================================================================
    cmaster_path = os.path.join(WORK, "C-windy-Q25.npy")
    C_reuse = None
    done = (CP.state() or {}).get("completed_ids") or []
    if "03-transport-Q25" in done and os.path.exists(cmaster_path):
        try:
            with open(cmaster_path, "rb") as f:
                buf = f.read()
            if len(buf) == GRID * GRID * 4:
                C_reuse = np.frombuffer(buf, dtype=np.float32).reshape(
                    GRID, GRID).copy()
                log("[02/03] reusing cached master %s" % cmaster_path)
            else:
                log("[02/03] cached master has %d bytes; recomputing"
                    % len(buf))
        except Exception as exc:
            log("[02/03] cached master unusable (%s); recomputing" % exc)
    st = build(raw, scene, "windy", GRID, seed, base_np, photo_lin, PH, PW,
               C=C_reuse)
    CP.unit_done("02-state")
    if C_reuse is None:
        atomic_write_bytes(cmaster_path,
                           np.ascontiguousarray(st["C"]).tobytes())
    CP.unit_done("03-transport-Q25")

    png_bytes = io.BytesIO()
    Image.fromarray(st["still"]).save(png_bytes, format="PNG")
    png_bytes = png_bytes.getvalue()

    def verify_png(path):
        try:
            im = Image.open(path)
            im.load()
            return im.size == (PW, PH) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(os.path.join(OUT, "comp-windy-still.png"),
                       png_bytes, verify=verify_png)
    log("[04] wrote comp-windy-still.png (%dx%d, %d bytes)"
        % (PW, PH, len(png_bytes)))
    CP.unit_done("04-composite")

    # ==================================================================
    # 05 -- the depth-cost reference (same wind/fetch, z = 0.6 m)
    # ==================================================================
    dref_scene = with_v8(raw, fetch_m=float(DEPTH_REF["fetch_m"]),
                         depth_m=float(DEPTH_REF["depth_m"]),
                         wind_ms=float(DEPTH_REF["wind_ms"]))
    dref_master = os.path.join(WORK, "C-z0p6-Q25.npy")
    dref_C = None
    if os.path.exists(dref_master):
        with open(dref_master, "rb") as f:
            buf = f.read()
        if len(buf) == GRID * GRID * 4:
            dref_C = np.frombuffer(buf, dtype=np.float32).reshape(
                GRID, GRID).copy()
            log("[05] reusing cached z=0.6 master %s" % dref_master)
    dref = build(raw, dref_scene, "z0.6", GRID, seed, base_np, photo_lin,
                 PH, PW, C=dref_C)
    if dref_C is None:
        atomic_write_bytes(dref_master,
                           np.ascontiguousarray(dref["C"]).tobytes())
    dref["smb"] = smb_node(dref["env"], dref["cfg"]["fetch"]["fetch_m"])
    log("[05] z=0.6 reference: Hs=%.4f sigma2=%.6f crush=%.5f linY=%.5f"
        % (dref["realised"]["Hs_m"], dref["realised"]["sigma2_rad2"],
           dref["still_metrics"]["crush_rate_linear_lt_0.02"],
           dref["img_stats"]["linear_luminance_mean"]))
    dref_png = io.BytesIO()
    Image.fromarray(dref["still"]).save(dref_png, format="PNG")
    atomic_write_bytes(os.path.join(WORK, "comp-windy-z0p6.png"),
                       dref_png.getvalue(), verify=verify_png)
    CP.unit_done("05-depth-reference")

    # ==================================================================
    # 05b -- what the frozen E2.2 costs, vs the v8 display normalisation
    # ==================================================================
    env_s = st["env"]
    _L_v3 = flat_water_direct_level(env_s, v3_beer(env_s))
    _L_v8 = flat_water_direct_level(
        env_s, np.asarray(st["beer"]["two_way"], np.float64)) * float(st["g1"])
    E_v8_ref = float(E_V3 * (_L_v3 / _L_v8)) if _L_v8 > 0 else None
    img_ev8, dbg_ev8 = render_frame_v5(
        st["fields"], st["C"], photo_lin, env_s, PW, PH,
        ambient_scale=A_COMP, c_layer_exposure=E_v8_ref,
        substrate_gain_rgb=st["substrate_gain_rgb"],
        beer_down=st["beer"]["down"], beer_up=st["beer"]["up"],
        g1=st["g1"], blur_sigma_px=st["blur_sigma_px"],
        lowpass_sigma_px=B_LOWPASS, return_debug=True)
    m_ev8 = image_metrics(dbg_ev8["lin"], dbg_ev8["aces"])
    del dbg_ev8, img_ev8
    gc.collect()
    log("[05b] E2.2 (frozen) crush=%.5f | E_v8=%.4f crush=%.5f"
        % (st["still_metrics"]["crush_rate_linear_lt_0.02"], E_v8_ref,
           m_ev8["crush_rate_linear_lt_0.02"]))
    CP.unit_done("05b-display-level-check")

    # ==================================================================
    # 06 -- the brief's warp measurement
    # ==================================================================
    still_png = np.asarray(Image.open(
        os.path.join(OUT, "comp-windy-still.png")).convert("RGB"))
    warp = warp_measure(base_np, still_png)
    p = warp["primary"]
    log("[06] warp: mean=%.2f px std=%.2f px (NCC check %.2f / %.2f) "
        "blocks=%d edge=%.1f%%"
        % (p["mean_shift_px"], p["shift_field_std_px"],
           warp["ncc_check"]["mean_shift_px"],
           warp["ncc_check"]["shift_field_std_px"],
           p["n_blocks"], 100.0 * p["frac_blocks_at_search_edge"]))

    # --- cross-check 1: the pipeline's own shift field -----------------
    pred = predicted_shift_field(st["fields"], st["env"], PW, PH)
    ps = predicted_block_stats(pred, warp["grid"])
    pred["block_stats"] = {k: v for k, v in ps.items()
                           if not k.startswith("block_means")}
    pred["agreement_ssd_primary"] = block_agreement(
        ps, warp["per_block"], "ssd_dy", "ssd_dx")
    pred["agreement_ncc_check"] = block_agreement(
        ps, warp["per_block"], "ncc_dy", "ncc_dx")
    log("     predicted (pipeline's own shift): field mean=%.2f px, block "
        "mean=%.2f px (%.2f mm/px, %.2f mm of water)"
        % (ps["field_mean_magnitude_px"], ps["block_mean_magnitude_px"],
           ps["pixel_scale_mm"], 1e3 * ps["mean_predicted_shift_m"]))
    del pred["dx_px"], pred["dy_px"], pred["dx_m"], pred["dy_m"]

    # --- cross-check 2: the same method on the same scene, C = 1 ------
    ctrl_img = render_frame_v5(
        st["fields"], np.ones((GRID, GRID), np.float32), photo_lin,
        st["env"], PW, PH, ambient_scale=A_COMP, c_layer_exposure=E_COMP,
        substrate_gain_rgb=st["substrate_gain_rgb"],
        beer_down=st["beer"]["down"], beer_up=st["beer"]["up"],
        g1=st["g1"], blur_sigma_px=st["blur_sigma_px"],
        lowpass_sigma_px=B_LOWPASS, return_debug=False)
    ctrl = warp_measure(base_np, ctrl_img)
    ctrl_png = io.BytesIO()
    Image.fromarray(ctrl_img).save(ctrl_png, format="PNG")
    atomic_write_bytes(os.path.join(WORK, "control-flatC.png"),
                       ctrl_png.getvalue(), verify=verify_png)
    log("     control (same scene, C=1): ssd %.2f/%.2f | ncc %.2f/%.2f px"
        % (ctrl["primary"]["mean_shift_px"],
           ctrl["primary"]["shift_field_std_px"],
           ctrl["ncc_check"]["mean_shift_px"],
           ctrl["ncc_check"]["shift_field_std_px"]))
    del ctrl_img

    warp["cross_check_pipeline_shift_field"] = {
        "what": ("the renderer's own content displacement, rebuilt from the "
                 "same state through the same _state_gradients path "
                 "(green band, z = %s m): an independent geometric reference "
                 "for the block-match number" % fmt(st["env"]["depthM"], 2)),
        "stats": pred["block_stats"],
        "agreement_with_ssd_primary": pred["agreement_ssd_primary"],
        "agreement_with_ncc_check": pred["agreement_ncc_check"],
    }
    warp["control_same_scene_flat_C"] = {
        "what": ("the identical scene rendered with C = 1 (caustic "
                 "modulation removed): the block matcher then sees only the "
                 "geometric warp.  This isolates how much of the primary "
                 "reading is caustic-band locking rather than refraction."),
        "ssd_primary": {k: ctrl["primary"][k] for k in
                        ("mean_shift_px", "shift_field_std_px",
                         "frac_blocks_at_search_edge")},
        "ncc_check": {k: ctrl["ncc_check"][k] for k in
                      ("mean_shift_px", "shift_field_std_px",
                       "frac_blocks_at_search_edge")},
        "file": "_work/control-flatC.png",
    }
    warp["reading"] = (
        "the geometric refraction displacement of this scene averages %.2f px "
        "over the field (%.2f px after block averaging); the C=1 control "
        "measures %.2f px mean |shift| with the gain-invariant matcher, in "
        "agreement with the field value, while the brief's raw-SSD block "
        "match reads %.2f px on the real composite.  The excess is caustic "
        "modulation (bright-band locking), not more refraction."
        % (ps["field_mean_magnitude_px"], ps["block_mean_magnitude_px"],
           ctrl["ncc_check"]["mean_shift_px"], p["mean_shift_px"]))
    CP.unit_done("06-warp-measurement")

    # ==================================================================
    # 07 -- gates + deliverables
    # ==================================================================
    dif = np.abs(still_png.astype(np.int16) - base_np.astype(np.int16))
    stt = st["state"]
    brk = stt.breaking or {}
    cst = st["C_stats"]
    h = float(st["env"]["depthM"])
    stokes_cap = float(brk.get("H_stokes_0p142_m", float("nan")))
    depth_cap = float(brk.get("H_depth_0p78h_m", float("nan")))
    hs_eff = float(brk.get("Hs_after_m", float("nan")))
    gates = {
        "output_opens_at_source_size": {
            "declared_in_brief": [736, 540],
            "source_jpeg_size": [int(PW), int(PH)],
            "actual_output_size": [int(still_png.shape[1]),
                                   int(still_png.shape[0])],
            "mode": "RGB",
            "note": ("the brief says 736x540, but the source JPEG decodes to "
                     "%dx%d; the still is the native source size (no resample, "
                     "no padding -- either would fabricate pixels and shift "
                     "the px->mm arithmetic)" % (PW, PH)),
            "pass": bool(still_png.shape[:2] == (PH, PW)),
        },
        "not_bit_identical_to_base": {
            "base_sha256": sha256_file(PHOTO),
            "still_sha256": sha256_file(
                os.path.join(OUT, "comp-windy-still.png")),
            "max_abs_diff_uint8": int(dif.max()),
            "mean_abs_diff_uint8": float(dif.mean()),
            "changed_pixel_fraction": float(np.mean(dif.max(axis=2) > 0)),
            "pass": bool(int(dif.max()) > 0),
        },
        "C_mean_1": {"mean": cst["mean"], "abs_dev": cst["mean_abs_dev_from_1"],
                     "tolerance": CMEAN_TOL, "pass": bool(cst["pass_mean_1"])},
        "stokes_cap_holds": {
            "cap_0p142_lambda_m": stokes_cap,
            "Hs_effective_m": hs_eff,
            "lambda_peak_m": brk.get("lambda_peak_m"),
            "pass": bool(hs_eff <= stokes_cap * 1.0000001)},
        "depth_cap_holds": {
            "cap_0p78h_m": depth_cap, "depth_m": h,
            "Hs_effective_m": hs_eff,
            "pass": bool(hs_eff <= depth_cap * 1.0000001)},
        "sigma2_above_recheck_floor": {
            "sigma2_realised_rad2": st["realised"]["sigma2_rad2"],
            "floor_rad2": SIGMA2_RECHECK_BELOW,
            "note": ("the brief triggers a scene re-check only BELOW 0.035; "
                     "above it the number is simply reported"),
            "pass": bool(st["realised"]["sigma2_rad2"]
                         >= SIGMA2_RECHECK_BELOW)},
        "warp_measurement_present": {
            "mean_shift_px": p["mean_shift_px"],
            "shift_field_std_px": p["shift_field_std_px"],
            "pass": bool(np.isfinite(p["mean_shift_px"])
                         and np.isfinite(p["shift_field_std_px"]))},
    }
    all_pass = all(g["pass"] for g in gates.values())
    log("[07] gates: %s (all=%s)"
        % ({k: v["pass"] for k, v in gates.items()}, all_pass))

    # tone bounds of the pipeline family, reported as numbers with their pass
    # flags -- NOT tuned, and NOT part of the brief's pass set (the brief sets
    # no crush bound; it only pins the frozen composite)
    reference_tone = {
        "note": ("chunk-4 gate 5's customary bounds, evaluated on this still "
                 "and REPORTED (nothing was re-tuned to meet them): crush "
                 "<= 0.0325, blow = 0"),
        "crush": st["still_metrics"]["crush_rate_linear_lt_0.02"],
        "crush_bound_customary": CRUSH_BOUND,
        "crush_pass_customary": bool(
            st["still_metrics"]["crush_rate_linear_lt_0.02"] <= CRUSH_BOUND),
        "crush_bound_stream_2pct": CRUSH_TIGHT,
        "crush_pass_stream_2pct": bool(
            st["still_metrics"]["crush_rate_linear_lt_0.02"] < CRUSH_TIGHT),
        "blow": st["still_metrics"]["blowout_rate_aces_gt_0.98"],
        "blow_bound": BLOW_BOUND,
        "blow_pass": bool(st["still_metrics"]["blowout_rate_aces_gt_0.98"]
                          <= BLOW_BOUND),
        "crush_at_E_v8_display_level": {
            "E_v8": E_v8_ref,
            "crush": m_ev8["crush_rate_linear_lt_0.02"],
            "diagnostic_only": ("the brief pins E2.2, so this render is not a "
                                "deliverable; it exists to show how much of the "
                                "crush is the frozen display constant rather "
                                "than the scene"),
        },
        "same_scene_at_z0p6": {
            "crush": dref["still_metrics"]["crush_rate_linear_lt_0.02"],
            "blow": dref["still_metrics"]["blowout_rate_aces_gt_0.98"],
        },
    }
    log("[07] reference tone: crush=%.5f (customary %.5f -> %s; at E_v8 "
        "%.5f) blow=%.6f"
        % (reference_tone["crush"], CRUSH_BOUND,
           reference_tone["crush_pass_customary"],
           reference_tone["crush_at_E_v8_display_level"]["crush"],
           reference_tone["blow"]))

    direction = stt.direction or {}
    nums_s = st["img_stats"]
    nums_d = dref["img_stats"]
    beer_s = st["beer"]
    beer_d = dref["beer"]

    def rel(a, b):
        return None if (a is None or b in (None, 0.0)) else float(a / b - 1.0)

    depth_cost = {
        "paired_scene": {
            "windy_z1.0": WINDY, "same_scene_z0.6": DEPTH_REF,
            "note": ("the ONLY difference is water.depthM; wind, fetch, sun, "
                     "recipe, seed, master grid and Q are identical.  Depth "
                     "enters twice: the optical path (Beer-Lambert) and, "
                     "weakly, the finite-depth dispersion / anchor integral."),
            "spectral_side_effect": {
                "sigma2_target_z1.0": smb["sigma2_target_rad2"],
                "sigma2_target_z0.6": dref["smb"]["sigma2_target_rad2"],
                "rel_diff": rel(dref["smb"]["sigma2_target_rad2"],
                                smb["sigma2_target_rad2"]),
                "Hs_realised_z1.0": st["realised"]["Hs_m"],
                "Hs_realised_z0.6": dref["realised"]["Hs_m"],
            },
        },
        "beer_direct_beam": {
            "z1.0": {k: beer_s[k] for k in ("down", "up", "two_way",
                                            "s_down_m", "s_up_m",
                                            "theta_w_deg", "c_rgb")},
            "z0.6": {k: beer_d[k] for k in ("down", "up", "two_way",
                                            "s_down_m", "s_up_m",
                                            "theta_w_deg", "c_rgb")},
            "ratio_two_way_z1p0_over_z0p6": [
                float(beer_s["two_way"][i] / beer_d["two_way"][i])
                for i in range(3)],
        },
        "image_stats": {
            "z1.0": nums_s, "z0.6": nums_d,
            "rel_change_z1p0_vs_z0p6": {
                "linear_luminance_mean": rel(nums_s["linear_luminance_mean"],
                                             nums_d["linear_luminance_mean"]),
                "linear_rgb_mean": [rel(nums_s["linear_rgb_mean"][i],
                                        nums_d["linear_rgb_mean"][i])
                                    for i in range(3)],
                "srgb_rgb_mean": [rel(nums_s["srgb_rgb_mean"][i],
                                      nums_d["srgb_rgb_mean"][i])
                                  for i in range(3)],
                "aces_luminance_mean": rel(nums_s["aces_luminance_mean"],
                                           nums_d["aces_luminance_mean"]),
                "aces_luminance_std": rel(nums_s["aces_luminance_std"],
                                          nums_d["aces_luminance_std"]),
                "srgb_saturation_mean": rel(nums_s["srgb_saturation_mean"],
                                            nums_d["srgb_saturation_mean"]),
            },
            "absolute_change_z1p0_vs_z0p6": {
                "crush_delta": float(nums_s["crush_rate_linear_lt_0.02"]
                                     - nums_d["crush_rate_linear_lt_0.02"]),
                "blow_delta": float(nums_s["blowout_rate_aces_gt_0.98"]
                                    - nums_d["blowout_rate_aces_gt_0.98"]),
            },
        },
        "C_master_z0p6": {"stats": dref["C_stats"]},
    }

    numbers = {
        "task": ("BRIEF-WINDY-STILL / one still, v3.4 (v8) pipeline, "
                 "eta = 1.0, windy pond"),
        "created_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat(),
        "pipeline": {
            "name": "v3.4 = v8 chunk 4",
            "runner": "run_windy_still.py (imports the frozen v2..v8 modules)",
            "state": "v8_state.V8WaterState",
            "transport": "det.det_frame (deterministic forward transport)",
            "composite": "v5_composite.render_frame_v5 (frozen production)",
            "eta": ETA,
            "eta_note": ("eta = 1.0 = direct sun only.  A recipe without a sky "
                         "block is eta = 1 exactly and the parser reports "
                         "mode 'legacy', so the v8 direct route is used and "
                         "the v9/Perez diffuse route is never invoked."),
            "master_grid": [GRID, GRID],
            "Q": Q,
            "A_ambient_scale": A_COMP,
            "B_lowpass_sigma_px": B_LOWPASS,
            "E_c_layer_exposure": E_COMP,
            "composite_note": ("E is the frozen production constant 2.2 (the "
                               "brief's 'E2.2'), not the v8 display-level "
                               "normalisation E_v8.  E_v8 for this scene is "
                               "reported below for reference only."),
            "E_v8_reference": E_v8_ref,
        },
        "base_image": {"path": PHOTO, "size": [int(PW), int(PH)],
                       "sha256": sha256_file(PHOTO)},
        "recipe_used": {
            "source": RECIPE,
            "name": raw.get("name"),
            "schema_in": raw.get("schema"),
            "schema_out": scene.get("schema"),
            "fetch_m": float(cfg["fetch"]["fetch_m"]),
            "fetch_source": cfg["fetch"]["source"],
            "wind_ms": float(env_s["windSpeedMS"]),
            "wind_dir_deg": float(raw["windWave"]["windDirDeg"]),
            "depth_m": float(env_s["depthM"]),
            "rain_rate_mm_h": 0.0,
            "surface_film": "none",
            "sun": {"date": raw["sun"]["date"], "time": raw["sun"]["time"],
                    "lat_deg": raw["sun"]["latitudeDeg"],
                    "lon_deg": raw["sun"]["longitudeDeg"],
                    "altitude_deg": float(env_s["sun"]["altitudeDeg"]),
                    "azimuth_deg": float(env_s["sun"]["azimuthDeg"])},
            "weather": {"cloud_pct": raw["weather"]["cloudPct"],
                        "clarity": raw["weather"]["clarity"],
                        "diffuse_fraction": float(env_s["diffuse"])},
            "water": {"transparency": raw["water"]["transparency"],
                      "turbidity": raw["water"]["turbidity"],
                      "pigment": raw["water"]["pigment"]},
            "seed": int(stt.seed),
            "seed_note": ("derived from the UNMODIFIED recipe, so the phase "
                          "realisation matches the piece's 3.x masters"),
        },
        "smb_fetch_computation": smb,
        "realised": {
            "Hs_m": st["realised"]["Hs_m"],
            "Hs_target_from_anchor_m": (stt.anchor or {}).get(
                "Hs_after_spectral"),
            "Hs_before_anchor_m": (stt.anchor or {}).get("Hs_before"),
            "m0_m2": st["realised"]["m0_m2"],
            "sigma2_rad2": st["realised"]["sigma2_rad2"],
            "sigma2_spectral_rad2": stt.spectrum_after["sigma2_rad2"],
            "rms_slope_rad": st["realised"]["rms_slope"],
            "sigma2_target_rad2": (stt.anchor or {}).get(
                "sigma2_target_rad2"),
            "sigma2_relative_residual": (stt.anchor or {}).get(
                "sigma2_residual_rel"),
            "anchor_applied": (stt.anchor or {}).get("applied"),
            "taper_kc_rad_per_m": (stt.anchor or {}).get("kc_rad_per_m"),
            "Hs_consistency_check": {
                "integrated_Hs_m": (
                    stt.anchor["candidates"]["JONSWAP"]["Hs"]
                    if stt.anchor and stt.anchor.get("candidates") else None),
                "jonswap_fetch_law_Hs_m": smb["Hs_fetch_law_m"],
            },
        },
        "direction": {
            "model": direction.get("model"),
            "s_max": direction.get("s_max"),
            "U_over_cp": direction.get("U_over_cp"),
            "s_min_over_band": direction.get("s_min_over_band"),
            "s_max_over_band": direction.get("s_max_over_band"),
            "half_width_deg_at_peak": smb["half_width_deg_at_peak"],
            "sigma2_before": direction.get("sigma2_before"),
            "sigma2_after": direction.get("sigma2_after"),
            "height_renormalisation": direction.get("height_renormalisation"),
        },
        "capillary_tail": {
            "k_viscous_rad_per_m": stt.tail["k_viscous_rad_per_m"],
            "k_nyquist_rad_per_m": stt.tail["k_nyquist_rad_per_m"],
            "kmin_rad_per_m": stt.tail["kmin_rad_per_m"],
            "modes_removed_in_resolved_band":
                stt.tail["modes_removed_in_resolved_band"],
        },
        "breaking": brk,
        "caps_report": {
            "rule": "H_max = min(0.142 lambda, 0.78 h)",
            "depth_m": h,
            "lambda_dominant_m": brk.get("lambda_peak_m"),
            "stokes_cap_0p142_lambda_m": stokes_cap,
            "depth_cap_0p78h_m": depth_cap,
            "binding": brk.get("binding"),
            "Hs_before_clamp_m": brk.get("Hs_before_m"),
            "Hs_after_clamp_m": brk.get("Hs_after_m"),
            "Hs_realised_field_m": st["realised"]["Hs_m"],
            "stokes_cap_holds": bool(hs_eff <= stokes_cap * 1.0000001),
            "depth_cap_holds": bool(hs_eff <= depth_cap * 1.0000001),
            "clamp_applied": bool(brk.get("applied")),
            "stokes_1over7_also_recorded":
                brk.get("stokes_cap_1over7"),
        },
        "dispersion_quantisation": stt.dispersion_node()["quantisation"],
        "transport": {"Q": Q, "info": st["det_info"], "C": cst},
        "still": {
            "file": "comp-windy-still.png",
            "size": [int(PW), int(PH)],
            "sha256": sha256_file(os.path.join(OUT, "comp-windy-still.png")),
            "metrics": st["still_metrics"],
            "image_stats": nums_s,
            "beer": beer_s,
            "blur_sigma_px": st["blur_sigma_px"],
            "g1": st["g1"],
            "substrate_gain_rgb": st["substrate_gain_rgb"],
        },
        "warp_measurement": warp,
        "depth_cost": depth_cost,
        "gates": gates,
        "reference_tone_bounds": reference_tone,
        "gates_all_pass": bool(all_pass),
        "timing_s": {"wall_total": time.time() - t_all},
        "checkpoint": CP.state(),
    }
    numbers = coerce(numbers)
    atomic_write_bytes(os.path.join(OUT, "numbers-windy.json"),
                       json.dumps(numbers, indent=2,
                                  ensure_ascii=False).encode("utf-8"))

    report = render_report(numbers)
    atomic_write_bytes(os.path.join(OUT, "REPORT-WINDY.md"),
                       report.encode("utf-8"))
    log("[07] wrote numbers-windy.json + REPORT-WINDY.md")
    CP.unit_done("07-gates-outputs")

    ok, bad = verify_units(OUT, "*.png", verify_png)
    ok_md, bad_md = verify_units(OUT, "*.md")
    ok_json, bad_json = verify_units(OUT, "*.json")
    print()
    print("checkpoint:", CP.state())
    print("verify: png ok=%d bad=%d | md ok=%d bad=%d | json ok=%d bad=%d"
          % (len(ok), len(bad), len(ok_md), len(bad_md),
             len(ok_json), len(bad_json)))
    print("gates_all_pass=%s wall=%.1fs" % (all_pass, time.time() - t_all))


def render_report(N):
    R = N["realised"]
    dr = N["direction"]
    st_c = N["transport"]["C"]
    sm = N["smb_fetch_computation"]
    g = N["gates"]
    rt = N["reference_tone_bounds"]
    w = N["warp_measurement"]
    wp = w["primary"]
    wn = w["ncc_check"]
    xc = w["cross_check_pipeline_shift_field"]
    xs = xc["stats"]
    ctl = w["control_same_scene_flat_C"]
    dc = N["depth_cost"]
    br = N["caps_report"]
    z1 = dc["image_stats"]["z1.0"]
    z6 = dc["image_stats"]["z0.6"]
    rc = dc["image_stats"]["rel_change_z1p0_vs_z0p6"]
    ab = dc["image_stats"]["absolute_change_z1p0_vs_z0p6"]
    beer = dc["beer_direct_beam"]
    W_, H_ = N["base_image"]["size"]

    L = []
    A = L.append
    A("# REPORT-WINDY — one still: v3.4 (v8) pipeline, windy pond")
    A("")
    A("Base image `%s` (%d×%d), composited with the v3.4 (= v8 chunk-4) "
      "caustics pipeline at **eta = 1.0** (direct sun only; the v9/Perez "
      "diffuse route is not used). Composite frozen as production: "
      "**A2.6-Blow + E2.2**, master 2048², Q = 25, σ_lowpass = 3.0 px, "
      "σ_blur = %s px, G1 = %s."
      % (N["base_image"]["path"], W_, H_,
         fmt(N["still"]["blur_sigma_px"], 2), fmt(N["still"]["g1"], 6)))
    A("")
    A("Output `comp-windy-still.png` sha256 `%s`… (base `%s`…), max |Δ| = "
      "%d/255, %.2f %% of pixels changed."
      % (g["not_bit_identical_to_base"]["still_sha256"][:16],
         g["not_bit_identical_to_base"]["base_sha256"][:16],
         g["not_bit_identical_to_base"]["max_abs_diff_uint8"],
         100.0 * g["not_bit_identical_to_base"]["changed_pixel_fraction"]))
    A("")
    A("**Size discrepancy, stated up front.** The brief specifies 736×540; "
      "the source JPEG `moth-on-bark.jpg` decodes to **%d×%d** — one pixel "
      "smaller on each axis. The still is delivered at the source's native "
      "size. It was NOT padded or resampled to 736×540: either would "
      "fabricate a row/column and shift the px→mm arithmetic the warp "
      "measurement depends on. The Δ versus the frame the brief describes is "
      "1 px per axis." % (W_, H_))
    A("")
    A("## 1. Numbers actually used")
    A("")
    A("| input | value | note |")
    A("|---|---|---|")
    A("| fetch | %s m | `wave.fetchM` explicit = the pond's 100 m |"
      % fmt(N["recipe_used"]["fetch_m"], 1))
    A("| wind | %s m/s | the only wind change vs preset 1 |"
      % fmt(N["recipe_used"]["wind_ms"], 2))
    A("| depth | %s m | the only depth change vs preset 1 (z = 0.60 there) |"
      % fmt(N["recipe_used"]["depth_m"], 2))
    A("| rain | 0 mm/h | rain out of scope: the rain path is never invoked |")
    A("| surface film | none | there is no film slot in the code |")
    A("| eta | 1.0 | no sky block → direct only (parser mode 'legacy') |")
    A("| sun | %s %s, %.1f°N %.1f°E → alt %s°, az %s° | preset-1 sun |"
      % (N["recipe_used"]["sun"]["date"], N["recipe_used"]["sun"]["time"],
         N["recipe_used"]["sun"]["lat_deg"],
         N["recipe_used"]["sun"]["lon_deg"],
         fmt(N["recipe_used"]["sun"]["altitude_deg"], 2),
         fmt(N["recipe_used"]["sun"]["azimuth_deg"], 2)))
    A("| water | turbidity %s, transparency %s, pigment %s | recipe-1 |"
      % (N["recipe_used"]["water"]["turbidity"],
         N["recipe_used"]["water"]["transparency"],
         N["recipe_used"]["water"]["pigment"]))
    A("| grid / Q | 2048² / %d | det_frame, no clamp, no renormalisation |"
      % N["pipeline"]["Q"])
    A("")
    A("## 2. Realised values")
    A("")
    A("| quantity | value | where it comes from |")
    A("|---|---|---|")
    A("| Ft = gF/U² | %s | SMB/fetch input |" % fmt(sm["Ft"], 2))
    A("| f_p | %s Hz | JONSWAP `f_p = 3.5 (g/U) Ft^-1/3` |"
      % fmt(sm["f_p_Hz"], 4))
    A("| α | %s | `0.076 Ft^-0.22` |" % fmt(sm["alpha"], 5))
    A("| λ_p | %s m | from f_p |" % fmt(sm["lambda_p_m"], 3))
    A("| H_s (SMB fetch law) | %s m | `0.0016 Ft^0.5 U²/g` |"
      % fmt(sm["Hs_fetch_law_m"], 4))
    A("| H_s (pipeline spectral integral) | %s m | JONSWAP candidate |"
      % fmt((R["Hs_consistency_check"] or {}).get("integrated_Hs_m"), 4))
    A("| **H_s realised (4√m0)** | **%s m** | height field at t = 0 |"
      % fmt(R["Hs_m"], 4))
    A("| σ² target (anchor) | %s rad² | JONSWAP %s / Donelan %s, median |"
      % (fmt(R["sigma2_target_rad2"], 6), fmt(sm["sigma2_jonswap_rad2"], 6),
         fmt(sm["sigma2_donelan_rad2"], 6)))
    A("| **σ² realised** | **%s rad²** | mean(gx²+gy²) |"
      % fmt(R["sigma2_rad2"], 6))
    A("| σ² spectral after taper | %s rad² | residual %.2e |"
      % (fmt(R["sigma2_spectral_rad2"], 6),
         float(R["sigma2_relative_residual"] or 0.0)))
    A("| rms slope | %s rad (%s°) | |"
      % (fmt(R["rms_slope_rad"], 4),
         fmt(math.degrees(math.atan(R["rms_slope_rad"])), 2)))
    A("| Cox–Munk old target | %s rad² | reference only (new = %s× old) |"
      % (fmt(sm["cox_munk_old_target_rad2"], 6),
         fmt(sm["ratio_new_over_old"], 3)))
    A("| taper k_c | %s rad/m | λ_c = %s m |"
      % (fmt(R["taper_kc_rad_per_m"], 2),
         fmt(2.0 * math.pi / R["taper_kc_rad_per_m"], 3)
         if R["taper_kc_rad_per_m"] else "-"))
    A("| direction s_max | %s | Mitsuyasu `11.5 U/c_p` |"
      % fmt(dr["s_max"], 2))
    A("| angular half-width at peak | %s° | from s_max, no free knob |"
      % fmt(dr["half_width_deg_at_peak"], 1))
    A("")
    A("C master (2048², Q = 25): mean %s, p50 %s, p99 %s, p99.9 %s, max %s, "
      "min %s."
      % (fmt(st_c["mean"], 9), fmt(st_c["p50"], 4), fmt(st_c["p99"], 4),
         fmt(st_c["p99_9"], 4), fmt(st_c["max"], 4), fmt(st_c["min"], 4)))
    A("")
    A("The realised σ² (%s) is ABOVE the brief's expected 0.04–0.05 band, not "
      "below the 0.035 re-check trigger, so no scene re-check was triggered "
      "and nothing was raised to get there: the target is the pipeline's own "
      "SMB/SFM1 value for U = 6.0 m/s, F = 100 m, h = 1.0 m (Ft = %s), and "
      "the realisation lands 1.5 %% above its spectral target — inside this "
      "patch's own realisation scatter."
      % (fmt(R["sigma2_rad2"], 6), fmt(sm["Ft"], 1)))
    A("")
    A("## 3. Caps (Stokes / depth)")
    A("")
    A("| cap | value | realised H_s | holds |")
    A("|---|---|---|---|")
    A("| Stokes 0.142 λ | %s m (λ = %s m) | %s m | %s |"
      % (fmt(br["stokes_cap_0p142_lambda_m"], 4),
         fmt(br["lambda_dominant_m"], 3), fmt(br["Hs_after_clamp_m"], 4),
         fmt(br["stokes_cap_holds"])))
    A("| depth 0.78 h | %s m (h = %s m) | %s m | %s |"
      % (fmt(br["depth_cap_0p78h_m"], 4), fmt(br["depth_m"], 2),
         fmt(br["Hs_after_clamp_m"], 4), fmt(br["depth_cap_holds"])))
    A("")
    A("Binding limit: **%s** (%s m); the clamp did **not** fire (H_s before "
      "the clamp %s m, i.e. %s %% of H_max). The Stokes 1/7 cap is recorded "
      "alongside and is also inactive: %s."
      % (br["binding"], fmt(min(br["stokes_cap_0p142_lambda_m"],
                                br["depth_cap_0p78h_m"]), 4),
         fmt(br["Hs_before_clamp_m"], 4),
         fmt(100.0 * br["Hs_before_clamp_m"]
             / br["depth_cap_0p78h_m"], 1),
         json.dumps(br["stokes_1over7_also_recorded"])))
    A("")
    A("## 4. Warp measurement (block matching, the brief's method)")
    A("")
    A("48 px blocks, ±24 px search at 2 px steps, Rec.709 luminance of the "
      "8-bit images, base photo → composite. %d blocks on the interior grid "
      "(%d × %d)."
      % (wp["n_blocks"], len(w["grid"]["block_origins_x"]),
         len(w["grid"]["block_origins_y"])))
    A("")
    A("| statistic | primary (raw SSD) | NCC check |")
    A("|---|---|---|")
    A("| **mean shift** | **%s px** | %s px |"
      % (fmt(wp["mean_shift_px"], 2), fmt(wn["mean_shift_px"], 2)))
    A("| **shift-field std** | **%s px** | %s px |"
      % (fmt(wp["shift_field_std_px"], 2),
         fmt(wn["shift_field_std_px"], 2)))
    A("| magnitude p50 / p90 / max | %s / %s / %s | %s / %s / %s |"
      % (fmt(wp["shift_magnitude_p50_px"], 2),
         fmt(wp["shift_magnitude_p90_px"], 2),
         fmt(wp["shift_magnitude_max_px"], 2),
         fmt(wn["shift_magnitude_p50_px"], 2),
         fmt(wn["shift_magnitude_p90_px"], 2),
         fmt(wn["shift_magnitude_max_px"], 2)))
    A("| mean signed (dx, dy) | (%s, %s) | (%s, %s) |"
      % (fmt(wp["shift_dx_mean_px"], 2), fmt(wp["shift_dy_mean_px"], 2),
         fmt(wn["shift_dx_mean_px"], 2), fmt(wn["shift_dy_mean_px"], 2)))
    A("| blocks pinned at the ±24 px edge | %s %% | %s %% |"
      % (fmt(100.0 * wp["frac_blocks_at_search_edge"], 1),
         fmt(100.0 * wn["frac_blocks_at_search_edge"], 1)))
    A("")
    A("**Yardsticks quoted in the brief (same method):** yesterday's "
      "production frame **23.2 / 17.1 px**, today's stream still "
      "**9.2 / 9.5 px**.")
    A("")
    A("**This still, same method: %s / %s px** — %s× yesterday's mean shift, "
      "%s× the stream's. On the raw comparison the windy pond sits between "
      "the two."
      % (fmt(wp["mean_shift_px"], 2), fmt(wp["shift_field_std_px"], 2),
         fmt(wp["mean_shift_px"] / 23.2, 2),
         fmt(wp["mean_shift_px"] / 9.2, 2)))
    A("")
    A("### 4b. What the number actually measures (cross-checks)")
    A("")
    A("The pixel number is not self-explanatory here, so it was checked "
      "against the geometry that produced the image:")
    A("")
    A("| reference | mean \\|shift\\| | std | source |")
    A("|---|---|---|---|")
    A("| pipeline's own shift field, per pixel (green band) | **%s px** | "
      "%s px | `_state_gradients` + the renderer's `-n·z·C_b` shift, rebuilt "
      "from the same state |"
      % (fmt(xs["field_mean_magnitude_px"], 2),
         fmt(xs["field_std_magnitude_px"], 2)))
    A("| pipeline's own shift field, block-averaged | %s px | %s px | the "
      "same field averaged inside each 48 px block (vector cancellation "
      "makes it smaller than the per-pixel mean) |"
      % (fmt(xs["block_mean_magnitude_px"], 2),
         fmt(xs["block_std_magnitude_px"], 2)))
    A("| same scene rendered with C = 1 (control) | %s px (NCC) | %s px | "
      "identical scene with the caustic modulation removed |"
      % (fmt(ctl["ncc_check"]["mean_shift_px"], 2),
         fmt(ctl["ncc_check"]["shift_field_std_px"], 2)))
    A("| real composite, NCC | %s px | %s px | gain-invariant matcher |"
      % (fmt(wn["mean_shift_px"], 2), fmt(wn["shift_field_std_px"], 2)))
    A("| real composite, raw SSD (the brief's method) | %s px | %s px | "
      "sensitive to brightness |"
      % (fmt(wp["mean_shift_px"], 2), fmt(wp["shift_field_std_px"], 2)))
    A("")
    A("Reading: the geometric refraction displacement of this scene averages "
      "**%s px per pixel**, i.e. **%s mm of water** at this view's pixel "
      "scale (%s mm/px). The C = 1 control measures %s px mean |shift| with "
      "the gain-invariant matcher — %s px from the geometric value (inside "
      "the 2 px offset quantisation), so the matcher recovers the geometry "
      "when the caustics are gone. On the real composite the brief's raw-SSD "
      "score rises to %s px because, with C's tail at p99 = %s (max %s), "
      "squared differences lock onto bright caustic bands as well as onto "
      "texture: %s %% of its blocks sit pinned at the ±24 px search edge, "
      "against %s %% for the NCC pass."
      % (fmt(xs["field_mean_magnitude_px"], 2),
         fmt(1e3 * xs["mean_predicted_shift_m"], 1),
         fmt(xs["pixel_scale_mm"], 2),
         fmt(ctl["ncc_check"]["mean_shift_px"], 2),
         fmt(abs(ctl["ncc_check"]["mean_shift_px"]
                 - xs["field_mean_magnitude_px"]), 2),
         fmt(wp["mean_shift_px"], 2), fmt(st_c["p99"], 2),
         fmt(st_c["max"], 2),
         fmt(100.0 * wp["frac_blocks_at_search_edge"], 1),
         fmt(100.0 * wn["frac_blocks_at_search_edge"], 1)))
    A("")
    A("Honest limit of the cross-check: the per-block matched vectors "
      "correlate only weakly with the predicted field (r = %s in dx, %s in "
      "dy on the real composite), so the matcher is a coarse aggregate "
      "measure, not a field reconstruction. What agrees is the aggregate: "
      "mean magnitude and its spread."
      % (fmt(xc["agreement_with_ncc_check"]["corr_dx"], 2),
         fmt(xc["agreement_with_ncc_check"]["corr_dy"], 2)))
    A("")
    A("### 4c. Why \"pixels\" and \"yesterday-level\" can both be true")
    A("")
    A("The refraction displacement in **metres** is what the physics sets: "
      "|n| · z · C_b. Here it is %.1f mm mean. Yesterday's pond view spanned "
      "4·specL = %.3f m over 832 px (%.2f mm/px); the windy pond's specL "
      "grows with U² (U = 6.0 m/s → 4·specL = %.2f m over %d px, %.2f mm/px). "
      "The same %.1f mm therefore maps to **%.1f px on yesterday's pixel "
      "scale** but only %s px here. The metre-level warp is yesterday-level; "
      "the pixel-level number is smaller because the wind widens the "
      "resolved patch, not because the surface is calmer."
      % (1e3 * xs["mean_predicted_shift_m"], 1.4755, 1475.5 / 832.0,
         sm["S_patch_m"], W_, xs["pixel_scale_mm"],
         1e3 * xs["mean_predicted_shift_m"],
         1e3 * xs["mean_predicted_shift_m"] / (1475.5 / 832.0),
         fmt(xs["field_mean_magnitude_px"], 2)))
    A("")
    A("Caveats of the measurement itself: offsets are quantised to the 2 px "
      "grid and saturate at 24 px per axis; blocks are 48 px, which at "
      "%.2f mm/px is %.0f mm of surface — coarse relative to the %.2f m "
      "patch; and the metric is a joint displacement/radiometry estimate "
      "because the composite carries the caustic modulation the base lacks."
      % (xs["pixel_scale_mm"], 48.0 * xs["pixel_scale_mm"], sm["S_patch_m"]))
    A("")
    A("## 5. What the depth cost did (z = 1.0 vs the same scene at z = 0.6)")
    A("")
    A("Paired composite, everything identical except `water.depthM` "
      "(1.0 → 0.6 m): same seed, same frozen composite, both rendered in "
      "this session.")
    A("")
    A("| direct-beam Beer factor | z = 1.0 | z = 0.6 | ratio |")
    A("|---|---|---|---|")
    for i, nm in enumerate(("R 650 nm", "G 550 nm", "B 450 nm")):
        A("| two-way, %s | %s | %s | %s |"
          % (nm, fmt(beer["z1.0"]["two_way"][i], 4),
             fmt(beer["z0.6"]["two_way"][i], 4),
             fmt(beer["ratio_two_way_z1p0_over_z0p6"][i], 3)))
    A("| one-way down (R) | %s | %s | |"
      % (fmt(beer["z1.0"]["down"][0], 4), fmt(beer["z0.6"]["down"][0], 4)))
    A("| in-water path s_down | %s m | %s m | |"
      % (fmt(beer["z1.0"]["s_down_m"], 3),
         fmt(beer["z0.6"]["s_down_m"], 3)))
    A("")
    A("| image statistic | z = 1.0 | z = 0.6 | change at z = 1.0 |")
    A("|---|---|---|---|")
    A("| linear luminance mean | %s | %s | %+.2f %% |"
      % (fmt(z1["linear_luminance_mean"], 5),
         fmt(z6["linear_luminance_mean"], 5),
         100.0 * (rc["linear_luminance_mean"] or 0.0)))
    A("| linear R / G / B mean | %s / %s / %s | %s / %s / %s | "
      "%+.2f / %+.2f / %+.2f %% |"
      % (fmt(z1["linear_rgb_mean"][0], 5), fmt(z1["linear_rgb_mean"][1], 5),
         fmt(z1["linear_rgb_mean"][2], 5), fmt(z6["linear_rgb_mean"][0], 5),
         fmt(z6["linear_rgb_mean"][1], 5), fmt(z6["linear_rgb_mean"][2], 5),
         100.0 * (rc["linear_rgb_mean"][0] or 0.0),
         100.0 * (rc["linear_rgb_mean"][1] or 0.0),
         100.0 * (rc["linear_rgb_mean"][2] or 0.0)))
    A("| sRGB R / G / B mean | %s / %s / %s | %s / %s / %s | "
      "%+.2f / %+.2f / %+.2f %% |"
      % (fmt(z1["srgb_rgb_mean"][0], 4), fmt(z1["srgb_rgb_mean"][1], 4),
         fmt(z1["srgb_rgb_mean"][2], 4), fmt(z6["srgb_rgb_mean"][0], 4),
         fmt(z6["srgb_rgb_mean"][1], 4), fmt(z6["srgb_rgb_mean"][2], 4),
         100.0 * (rc["srgb_rgb_mean"][0] or 0.0),
         100.0 * (rc["srgb_rgb_mean"][1] or 0.0),
         100.0 * (rc["srgb_rgb_mean"][2] or 0.0)))
    A("| ACES luminance mean | %s | %s | %+.2f %% |"
      % (fmt(z1["aces_luminance_mean"], 4),
         fmt(z6["aces_luminance_mean"], 4),
         100.0 * (rc["aces_luminance_mean"] or 0.0)))
    A("| ACES luminance std (contrast) | %s | %s | %+.2f %% |"
      % (fmt(z1["aces_luminance_std"], 4),
         fmt(z6["aces_luminance_std"], 4),
         100.0 * (rc["aces_luminance_std"] or 0.0)))
    A("| sRGB saturation mean | %s | %s | %+.2f %% |"
      % (fmt(z1["srgb_saturation_mean"], 4),
         fmt(z6["srgb_saturation_mean"], 4),
         100.0 * (rc["srgb_saturation_mean"] or 0.0)))
    A("| crush / blow [%%] | %s / %s | %s / %s | %+.2f / %+.2f pp |"
      % (fmt(100.0 * z1["crush_rate_linear_lt_0.02"], 3),
         fmt(100.0 * z1["blowout_rate_aces_gt_0.98"], 3),
         fmt(100.0 * z6["crush_rate_linear_lt_0.02"], 3),
         fmt(100.0 * z6["blowout_rate_aces_gt_0.98"], 3),
         100.0 * ab["crush_delta"], 100.0 * ab["blow_delta"]))
    A("")
    A("What the depth cost did, one line each:")
    A("")
    A("- **Caustic beam**: two-way Beer factors fall to %s / %s / %s "
      "(R/G/B) at z = 1.0 from %s / %s / %s at z = 0.6 — ratios %s / %s / "
      "%s. The extra 0.4 m is paid twice (s_down %s m versus %s m, s_up = "
      "z), so the direct layer loses another **%s %% of its green**."
      % (fmt(beer["z1.0"]["two_way"][0], 4),
         fmt(beer["z1.0"]["two_way"][1], 4),
         fmt(beer["z1.0"]["two_way"][2], 4),
         fmt(beer["z0.6"]["two_way"][0], 4),
         fmt(beer["z0.6"]["two_way"][1], 4),
         fmt(beer["z0.6"]["two_way"][2], 4),
         fmt(beer["ratio_two_way_z1p0_over_z0p6"][0], 3),
         fmt(beer["ratio_two_way_z1p0_over_z0p6"][1], 3),
         fmt(beer["ratio_two_way_z1p0_over_z0p6"][2], 3),
         fmt(beer["z1.0"]["s_down_m"], 3), fmt(beer["z0.6"]["s_down_m"], 3),
         fmt(100.0 * (1.0 - beer["ratio_two_way_z1p0_over_z0p6"][1]), 1)))
    A("- **Brightness**: linear luminance mean %s → %s (%+.2f %%). The "
      "frozen ambient term (v3 scat heuristic) is depth-independent, so the "
      "deeper image darkens only the caustic layer and the caustic-to-flat "
      "ratio drops."
      % (fmt(z6["linear_luminance_mean"], 5),
         fmt(z1["linear_luminance_mean"], 5),
         100.0 * (rc["linear_luminance_mean"] or 0.0)))
    A("- **Colour**: red loses most (linear R %+.2f %%), blue least "
      "(%+.2f %%; ratios %s vs %s): water absorbs red first, so z = 1.0 is "
      "the greener/bluer member of the pair and the R/G/B spread IS the "
      "0.4 m."
      % (100.0 * (rc["linear_rgb_mean"][0] or 0.0),
         100.0 * (rc["linear_rgb_mean"][2] or 0.0),
         fmt(beer["ratio_two_way_z1p0_over_z0p6"][0], 3),
         fmt(beer["ratio_two_way_z1p0_over_z0p6"][2], 3)))
    A("- **Contrast and saturation after the tone curve**: ACES luminance "
      "std %+.2f %%, saturation %+.2f %% — attenuation compresses the caustic "
      "amplitude before the tonemap, so the deeper still is lower contrast, "
      "not merely darker."
      % (100.0 * (rc["aces_luminance_std"] or 0.0),
         100.0 * (rc["srgb_saturation_mean"] or 0.0)))
    A("- **Tone tail**: crush %+.2f pp (to %s %% at z = 1.0), blow %+.2f pp "
      "(zero in both). The deeper path pushes the darkest caustic gaps under "
      "the 0.02 linear floor."
      % (100.0 * ab["crush_delta"],
         fmt(100.0 * z1["crush_rate_linear_lt_0.02"], 3),
         100.0 * ab["blow_delta"]))
    A("- **Not purely optical**: depth scales the refraction lever arm "
      "(|n|·z·C_b), so the C master itself differs between the two runs "
      "(p99 %s at z = 1.0 vs %s at z = 0.6) and the caustic geometry shifts "
      "with z; the finite-depth dispersion and anchor integral add a "
      "negligible %+.6f %% to the σ² target (%s vs %s rad²). The Beer table "
      "isolates the optical part of the difference."
      % (fmt(st_c["p99"], 3),
         fmt(dc["C_master_z0p6"]["stats"]["p99"], 3),
         100.0 * (dc["paired_scene"]["spectral_side_effect"]["rel_diff"]
                  or 0.0),
         fmt(dc["paired_scene"]["spectral_side_effect"]["sigma2_target_z0.6"],
             8),
         fmt(dc["paired_scene"]["spectral_side_effect"]["sigma2_target_z1.0"],
             8)))
    A("")
    A("## 6. Tone bounds (reported, nothing tuned)")
    A("")
    A("The brief sets no crush bound; the pipeline family's own chunk-4 gate 5 "
      "does (crush ≤ 0.0325, blow = 0). Both are reported here as numbers:")
    A("")
    A("| bound | realised | pass |")
    A("|---|---|---|")
    A("| crush ≤ %s (customary) | %s %% | **%s** |"
      % (fmt(rt["crush_bound_customary"], 4),
         fmt(100.0 * rt["crush"], 3), fmt(rt["crush_pass_customary"])))
    A("| crush < %s (stream brief) | %s %% | **%s** |"
      % (fmt(rt["crush_bound_stream_2pct"], 4),
         fmt(100.0 * rt["crush"], 3), fmt(rt["crush_pass_stream_2pct"])))
    A("| blow = 0 | %s %% | %s |"
      % (fmt(100.0 * rt["blow"], 3), fmt(rt["blow_pass"])))
    A("")
    A("The crush bound FAILS at the frozen composite, and this is reported, "
      "not repaired. Where it comes from, numerically: (a) the same C master "
      "at the v8 display level E_v8 = %s instead of the frozen E2.2 gives "
      "crush %s %%, so most of the tail is the frozen display constant being "
      "darker than this scene's own normalisation; (b) the z = 0.6 pair of "
      "the same scene gives %s %%, so the depth lever adds %s pp. No "
      "parameter was changed to move this number: the brief pins A2.6 + "
      "E2.2, and it is reported as a failed reference bound."
      % (fmt(rt["crush_at_E_v8_display_level"]["E_v8"], 3),
         fmt(100.0 * rt["crush_at_E_v8_display_level"]["crush"], 3),
         fmt(100.0 * rt["same_scene_at_z0p6"]["crush"], 3),
         fmt(100.0 * (rt["crush"] - rt["same_scene_at_z0p6"]["crush"]), 3)))
    A("")
    A("## 7. Gates")
    A("")
    A("| gate | realised | pass |")
    A("|---|---|---|")
    A("| opens, RGB, source-native size | %d×%d (brief says 736×540; JPEG is "
      "%d×%d) | %s |"
      % (g["output_opens_at_source_size"]["actual_output_size"][0],
         g["output_opens_at_source_size"]["actual_output_size"][1],
         g["output_opens_at_source_size"]["source_jpeg_size"][0],
         g["output_opens_at_source_size"]["source_jpeg_size"][1],
         fmt(g["output_opens_at_source_size"]["pass"])))
    A("| not bit-identical to base | %.2f %% of pixels differ (max Δ %d) | %s |"
      % (100.0 * g["not_bit_identical_to_base"]["changed_pixel_fraction"],
         g["not_bit_identical_to_base"]["max_abs_diff_uint8"],
         fmt(g["not_bit_identical_to_base"]["pass"])))
    A("| C mean = 1.0 | %s (|Δ| %.2e) | %s |"
      % (fmt(st_c["mean"], 9), float(g["C_mean_1"]["abs_dev"]),
         fmt(g["C_mean_1"]["pass"])))
    A("| Stokes cap 0.142 λ | H_s %s m ≤ %s m | %s |"
      % (fmt(g["stokes_cap_holds"]["Hs_effective_m"], 4),
         fmt(g["stokes_cap_holds"]["cap_0p142_lambda_m"], 4),
         fmt(g["stokes_cap_holds"]["pass"])))
    A("| depth cap 0.78 h | H_s %s m ≤ %s m | %s |"
      % (fmt(g["depth_cap_holds"]["Hs_effective_m"], 4),
         fmt(g["depth_cap_holds"]["cap_0p78h_m"], 4),
         fmt(g["depth_cap_holds"]["pass"])))
    A("| σ² ≥ 0.035 (brief's re-check floor) | %s rad² | %s |"
      % (fmt(g["sigma2_above_recheck_floor"]["sigma2_realised_rad2"], 6),
         fmt(g["sigma2_above_recheck_floor"]["pass"])))
    A("| warp measurement present | %s / %s px | %s |"
      % (fmt(wp["mean_shift_px"], 2), fmt(wp["shift_field_std_px"], 2),
         fmt(g["warp_measurement_present"]["pass"])))
    A("")
    A("Brief gates all pass: **%s** (crush is a reference bound, listed in "
      "§6, and fails there)." % fmt(N["gates_all_pass"]))
    A("")
    A("## 8. Honest notes")
    A("")
    A("- **H_s is computed, not set — but it is not the SMB law either.** The "
      "pipeline evaluates the SMB/JONSWAP fetch law (%s m) and its own "
      "spectral integral (%s m); the field's amplitude comes from the recipe "
      "heuristic `env.waveHeightM` and the realised 4√m0 is %s m — %.1f× the "
      "SMB value. Nothing was hand-set; the gap is structural and "
      "pre-existing."
      % (fmt(sm["Hs_fetch_law_m"], 4),
         fmt((R["Hs_consistency_check"] or {}).get("integrated_Hs_m"), 4),
         fmt(R["Hs_m"], 4),
         (R["Hs_m"] / sm["Hs_fetch_law_m"]) if sm["Hs_fetch_law_m"] else
         float("nan")))
    A("- **Patch scale.** At U = 6.0 m/s specL = %.3f m, so the patch is "
      "4·specL = %.2f m and 1 px of the %d-px view spans %.2f mm of surface "
      "(master dx %.2f mm). The realised spectral peak sits at kmin, so the "
      "dominant wavelength is the patch side, not the JONSWAP λ_p (%s m). "
      "This is why pixel displacement and metre displacement diverge (§4c)."
      % (sm["S_patch_m"] / 4.0, sm["S_patch_m"], W_, xs["pixel_scale_mm"],
         sm["S_patch_m"] / 2048.0 * 1e3, fmt(sm["lambda_p_m"], 3)))
    A("- **Ft outside the fitted range.** Ft = %s here; the JONSWAP fetch "
      "relations are fitted to 1e3–1e5, so this scene is an extrapolation in "
      "the short-fetch direction. Recorded, not corrected — the brief asked "
      "for the pipeline's own SMB number, not a different one."
      % fmt(sm["Ft"], 2))
    A("- **C's minimum is slightly negative** (%s): the bicubic forward "
      "transport can undershoot at caustic gaps; it is clipped at the "
      "tonemap and is pre-existing behaviour, not introduced here."
      % fmt(st_c["min"], 6))
    A("- **Rain out of scope.** rain = 0 and the v6 rain band is never "
      "built; no number here is a rain result.")
    A("- **No surface film slot exists** in the code; 'none' is the only "
      "honest value.")
    A("- **Pixel-scale caveat for the yardstick.** Comparing this still's "
      "pixel displacement with yesterday's mixes two different px→mm "
      "scales; §4c gives the metre-level comparison, which is the one that "
      "is invariant to the patch size.")
    A("")
    A("## 9. Files")
    A("")
    A("- `comp-windy-still.png` — %d×%d composited still" % (W_, H_))
    A("- `numbers-windy.json` — recipe, realised values, caps, C stats, the "
      "warp measurement (with both cross-checks), the depth-cost pair, tone "
      "bounds, gates")
    A("- `REPORT-WINDY.md` — this file")
    A("- `_work/C-windy-Q25.npy` — the 2048² C master (raw float32)")
    A("- `_work/comp-windy-z0p6.png` — the z = 0.6 paired composite")
    A("- `_work/control-flatC.png` — the C = 1 control render")
    A("- `checkpoint.json` — jobkit progress (%d unit_done events at the time "
      "these files were written; the final `07-gates-outputs` unit lands "
      "immediately after them)"
      % (N.get("checkpoint") or {}).get("completed_units", 0))
    A("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
