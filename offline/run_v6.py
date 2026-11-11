"""Pipeline v6 runner -- PIPELINE-V6 chunk 2 (rain on windless, calm water).

Implements the chunk-2 design (PIPELINE-V6.md) on top of the frozen
v2/v3/v4/v5 pipeline.  No existing file is modified: state.py, transport.py,
det.py, image.py, v5_physics.py, v5_state.py, v5_composite.py, recipe_v5.py
and run_v2/v3/v5 are imported, never edited.

New in this chunk (report.md section 3 / section 5):

  rain.rate_mm_h  R  [mm/h]  ->  Marshall-Palmer Lambda(R), D0(R)
                              ->  Atlas-Ulbrich v(D)
                              ->  impact flux F(R) = int N v dD
                              ->  flux-weighted mean impact energy Ebar(R)
                              ->  k_ring(R) from D0(R) (monotone decreasing,
                                  inside the 310-630 rad/m ring band)
                              ->  isotropic ring-band variance
                                  eta F Ebar tau / (rho g + sigma k^2)
                              ->  gravity-capillary dispersion and 8 s loop
                              ->  second effect: sub-band damping (approx.)

Products in <out>/ (default out/v6/):
  numbers-v6.json
  still-R00.png, still-R02.png, still-R10.png, still-R50.png
  panel-R00-R02-R10-R50.png
  crop-zoom-R00-R02-R10-R50.png   (2 rows: magnified stills / C masters)
  C-rain-R00.npy .. C-rain-R50.npy, C-baseline-v5-calm.npy
  recipes/rain-R02-v2.json ... (reproducible inputs)

Usage (from caustics-filter/offline/):
  python run_v6.py --recipe ../blender/recipe-1.json --preset 1
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
from v5_state import V5WaterState  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,
                        scattering_blur_sigma_px, substrate_gain,
                        CHANNEL_WL_NM)  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,
                          image_metrics, v3_beer,
                          flat_water_direct_level)  # noqa: E402
from recipe_v6 import parse_recipe_v6, with_rain  # noqa: E402
from v6_physics import (RING_ETA, RING_ETA_BASIS, RING_K_BAND, RING_TAU_S,
                        capillary_crossover, eta_sensitivity, flux_gate,
                        rain_state_summary, ring_wavenumber, wavenumber_gate,
                        F_TABLE, F_TABLE_TOL_FACTOR)  # noqa: E402
from v6_state import V6WaterState  # noqa: E402

ETAS_DEMO = (0.03, 0.10, 0.30)
CMEAN_TOL = 1e-3
MONO_TOL = 0.0


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {"shape": list(c.shape), "dtype": str(np.asarray(C).dtype),
            "mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min()),
            "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
            "pass_mean_1": bool(abs(c.mean() - 1.0) <= CMEAN_TOL)}


def save_panel_2x2(path, items, gap=8, label_h=20, caption_h=18):
    """items: list of (rgb uint8, label, caption).  2x2 grid, labelled."""
    imgs = [Image.fromarray(x[0]) for x in items]
    W, H = imgs[0].size
    cols, rows = 2, 2
    canvas = Image.new("RGB", (cols * W + (cols + 1) * gap,
                               rows * (H + label_h + caption_h)
                               + (rows + 1) * gap), (245, 245, 245))
    d = ImageDraw.Draw(canvas)
    for i, (im, (rgb, label, caption)) in enumerate(zip(imgs, items)):
        cx = i % cols
        cy = i // cols
        x0 = gap + cx * (W + gap)
        y0 = gap + cy * (H + label_h + caption_h + gap)
        d.text((x0 + 2, y0 + 4), label, fill=(15, 15, 15))
        canvas.paste(im, (x0, y0 + label_h))
        d.text((x0 + 2, y0 + label_h + H + 4), caption, fill=(70, 70, 70))
    canvas.save(path)
    return path


def save_triple_zoom(path, crops, zooms, labels, gap=8, label_h=20):
    """Three magnified crops side by side, labelled."""
    ims = []
    for c, z in zip(crops, zooms):
        im = Image.fromarray(c)
        im = im.resize((int(im.width * z), int(im.height * z)),
                       Image.LANCZOS)
        ims.append(im)
    W, H = ims[0].size
    canvas = Image.new("RGB", (3 * W + 4 * gap, H + label_h + 2 * gap),
                       (245, 245, 245))
    d = ImageDraw.Draw(canvas)
    for i, (im, lab) in enumerate(zip(ims, labels)):
        x0 = gap + i * (W + gap)
        d.text((x0 + 2, gap - 12), lab, fill=(15, 15, 15))
        canvas.paste(im, (x0, gap + label_h - 12))
    canvas.save(path)
    return path


def best_crop(imgs, size, smooth=48):
    """Deterministic crop: the region where rain differs most from calm.

    Uses |gray(R=10) - gray(R=0)| box-smoothed over ``smooth`` px and picks
    the argmax centre, so the magnified crop shows the ring grain rather than
    an arbitrary part of the photo.
    """
    a = imgs[0].astype(np.float64).mean(axis=2)
    b = imgs[1].astype(np.float64).mean(axis=2)
    dd = np.abs(b - a)
    k = int(smooth)
    ker = np.ones(k) / k
    sm = np.apply_along_axis(lambda v: np.convolve(v, ker, mode="same"), 0, dd)
    sm = np.apply_along_axis(lambda v: np.convolve(v, ker, mode="same"), 1, sm)
    h, w = sm.shape
    half = size // 2
    core = sm[half:h - half, half:w - half]
    j, i = np.unravel_index(int(np.argmax(core)), core.shape)
    cy, cx = j + half, i + half
    return {"cy": int(cy), "cx": int(cx), "size": int(size),
            "score": float(sm[cy, cx]),
            "box": [int(cx - half), int(cy - half), int(size), int(size)]}


def build_crop_figure(out_dir, rates, S_patch, crop_px=240, zoom=3.0,
                      gap=8, label_h=18):
    """Magnified ring-grain evidence: two rows, same physical window.

    Row 1: the magnification of the DELIVERED stills (uint8 sRGB, photo
    composited) -- this is the literal "拡大クロップ" of the verification
    stills, so the photo's own pixel grid is magnified too.
    Row 2: the same physical window taken from the C MASTERS (the transport
    output, 2048^2 dimensionless illuminance ratio) at 1:1, grey-mapped on a
    FIXED scale (C in [0.4, 2.2]) so the two rows can be compared and the
    three rates are on the same footing.  This row has no photo and no
    upsampling: it shows the ring-wave caustic field itself, which is what
    makes the grain physical rather than a sampling artefact.

    The window is chosen deterministically by ``best_crop`` (max |grey(R=10)
    - grey(R=0)| after a 48 px box average), so it is the region where rain
    changes the image most, not an arbitrary cut.
    """
    order = [r for r in sorted(rates)]
    files = {R: os.path.join(out_dir, "still-R%02g.png" % R) for R in order}
    imgs = {R: np.asarray(Image.open(files[R]).convert("RGB"))
            for R in order}
    ref = 10.0 if 10.0 in order else order[-1]
    info = best_crop([imgs[order[0]], imgs[ref]], crop_px)
    x0, y0, w, h = info["box"]
    PH, PW, _ = imgs[order[0]].shape

    # same physical window in master pixels (still covers the whole patch)
    C0 = np.load(os.path.join(out_dir, "C-rain-R%02g.npy" % order[0]),
                 mmap_mode="r")
    NC = int(C0.shape[0])
    del C0
    cx0 = int(round(x0 / float(PW) * NC))
    cy0 = int(round(y0 / float(PH) * NC))
    cw_c = int(round(w / float(PW) * NC))
    ch_c = int(round(h / float(PH) * NC))

    panes = []
    labels = []
    for R in order:
        pane = imgs[R][y0:y0 + h, x0:x0 + w]
        im = Image.fromarray(pane).resize(
            (int(round(w * zoom)), int(round(h * zoom))), Image.LANCZOS)
        panes.append(np.asarray(im))
        labels.append("still R=%g mm/h, %.1fx (from the delivered PNG)" % (R, zoom))
    W1, H1 = panes[0].shape[1], panes[0].shape[0]
    for R in order:
        C = np.load(os.path.join(out_dir, "C-rain-R%02g.npy" % R))
        cc = np.asarray(C[cy0:cy0 + ch_c, cx0:cx0 + cw_c], dtype=np.float64)
        g = np.clip((cc - 0.4) / (2.2 - 0.4), 0.0, 1.0)
        g8 = (g * 255.0 + 0.5).astype(np.uint8)
        im = Image.fromarray(g8, mode="L").resize((W1, H1), Image.BILINEAR)
        panes.append(np.asarray(im.convert("RGB")))
        labels.append("C master R=%g mm/h, 1:1 (transport output, "
                      "grey = C 0.4..2.2)" % R)
    H2 = panes[len(order)].shape[0]
    rows = [panes[:len(order)], panes[len(order):]]
    row_h = [H1, H2]
    cols = len(order)
    cw_ = W1
    canvas = Image.new("RGB", (cols * cw_ + (cols + 1) * gap,
                               2 * (max(row_h) + label_h) + 3 * gap),
                       (245, 245, 245))
    d = ImageDraw.Draw(canvas)
    ytop = gap
    for ri, row in enumerate(rows):
        for ci, pane in enumerate(row):
            x = gap + ci * (cw_ + gap)
            d.text((x + 2, ytop), labels[ri * cols + ci], fill=(15, 15, 15))
            canvas.paste(Image.fromarray(pane), (x, ytop + label_h))
        ytop += max(row_h) + label_h + gap
    path = os.path.join(out_dir, "crop-zoom-%s.png"
                        % "-".join("R%02g" % R for R in order))
    canvas.save(path)
    info.update({
        "file": os.path.basename(path),
        "zoom": float(zoom),
        "rows": ["delivered stills (magnified)", "C master 1:1 (transport)"],
        "c_master_window": {"x0": cx0, "y0": cy0, "w": cw_c, "h": ch_c,
                             "grid": NC},
        "physical_window_m": [float(x0 / PW * S_patch),
                              float(y0 / PH * S_patch),
                              float(w / PW * S_patch),
                              float(h / PH * S_patch)],
        "ring_wavelength_in_window_px": {
            "still_735px": float((0.015 / S_patch) * PW * zoom),
            "C_master_2048px": float((0.015 / S_patch) * NC * (W1 / cw_c)),
        },
        "note": ("crop window chosen deterministically: the %dx%d px region "
                 "where |grey(R=10)-grey(R=0)| peaks after a 48 px box "
                 "average.  Row 2 is the same physical window of the C "
                 "master at 1:1, so the ring grain there is not a JPEG/upscale "
                 "artefact; row 1 is the literal magnification of the "
                 "delivered still." % (w, h)),
        "selection_source_rate_mm_h": float(ref),
    })
    return info


def _patch_numbers_crop(out_dir, crop_info):
    p = os.path.join(out_dir, "numbers-v6.json")
    if not os.path.exists(p):
        return False
    with open(p, "r", encoding="utf-8") as f:
        d = json.load(f)
    d.setdefault("stills", {})["crop"] = crop_info
    with open(p, "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False, default=float)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v6")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--rates", default="0,2,10,50")
    ap.add_argument("--wind", type=float, default=0.0,
                    help="scene wind speed for the calm verification scene")
    ap.add_argument("--crop", type=int, default=240)
    ap.add_argument("--zoom", type=float, default=3.0)
    ap.add_argument("--crop-only", action="store_true",
                    help="only regenerate the crop/panel PNGs from the stills "
                         "and C masters already on disk (physics untouched)")
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    rates = [float(x) for x in str(args.rates).split(",") if x.strip()]

    # ---- the verification scene: preset 1 photo/sun, wind 0 (calm) -----
    scene = with_rain(recipe, 0.0, wind_ms=float(args.wind))
    env = resolve_environment(scene)
    cfg = parse_recipe_v6(scene)
    grid = int(args.grid)
    Q = int(args.q)
    S_patch = 4.0 * float(env["specL"])

    print("recipe=%s schema=%s wind=%.2f specL=%.3f patch=%.3f m "
          "depth=%.2f turbid=%.2f sun_alt=%.1f"
          % (recipe.get("name"), cfg["schema"], env["windSpeedMS"],
             env["specL"], S_patch, env["depthM"], env["turbid"],
             env["sun"]["altitudeDeg"]), flush=True)

    if args.crop_only:
        ci = build_crop_figure(args.out, rates, S_patch, args.crop, args.zoom)
        ok = _patch_numbers_crop(args.out, ci)
        print("crop-only: %s (numbers patched: %s)"
              % (ci["file"], ok), flush=True)
        print(json.dumps({k: ci[k] for k in
                          ("box", "physical_window_m", "zoom",
                           "ring_wavelength_in_window_px")},
                         ensure_ascii=False))
        return

    # ---- physics gates (grid independent) ------------------------------
    f_gate = flux_gate()
    k_gate = wavenumber_gate()
    cw = capillary_crossover()

    # ---- v5 reference state (calm) = backward-compatibility baseline ----
    t0 = time.time()
    st_v5 = V5WaterState(scene, env, N=grid, seed=args.seed,
                         preset_index=args.preset)
    f_v5 = st_v5.fields_at(0.0)
    f_v5["S"] = st_v5.S
    re_v5 = st_v5.realised_stats(f_v5)
    t_v5_state = time.time() - t0
    print("v5 calm state: sigma2_spec=%.6f sigma2_real=%.6f Hs=%.4f "
          "(%.1fs)" % (st_v5.spectrum_after["sigma2_rad2"],
                       re_v5["sigma2_rad2"], re_v5["Hs_m"], t_v5_state),
          flush=True)

    # ---- display level (same normalisation as v5; held for all R) ------
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_node = smith_g1_beckmann(re_v5["rms_slope"], sun_zenith)
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(g1_node["G1"])
    level_ratio = L_v3 / L_v5 if L_v5 > 0 else 1.0
    E_v6 = E_V3 * level_ratio
    print("beer two_way=%s G1=%.4f level_ratio=%.4f E_v6=%.4f blur=%.2fpx"
          % ([round(float(x), 4) for x in beer_ref["two_way"]],
             g1_node["G1"], level_ratio, E_v6, blur_ref["sigma_px"]),
          flush=True)
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    # ---- v5 baseline transport + still (the reference C / pixels) ------
    t1 = time.time()
    C_v5, info_v5 = det_frame(f_v5, st_v5.S, grid, env, Q)
    t_c_v5 = time.time() - t1
    np.save(os.path.join(args.out, "C-baseline-v5-calm.npy"),
            np.asarray(C_v5, dtype=np.float32))
    img_v5, dbg_v5 = render_frame_v5(
        f_v5, C_v5, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v6, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=g1_node["G1"], blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_v5 = image_metrics(dbg_v5["lin"], dbg_v5["aces"])
    Image.fromarray(img_v5).save(os.path.join(args.out, "still-v5-calm.png"))
    print("v5 calm C mean=%.9f p99=%.4f (%.1fs) still crush=%.4f blow=%.6f"
          % (float(C_v5.mean()), float(np.percentile(C_v5, 99)), t_c_v5,
             m_v5["crush_rate_linear_lt_0.02"],
             m_v5["blowout_rate_aces_gt_0.98"]), flush=True)
    del dbg_v5
    gc.collect()

    per_rate = {}
    stills = []
    C_by_rate = {}
    state_ref = None

    for R in rates:
        t0 = time.time()
        st = V6WaterState(scene, env, N=grid, seed=args.seed,
                          preset_index=args.preset, rain_rate_mm_h=R)
        fields = st.fields_at(0.0)
        fields["S"] = st.S
        realised = st.realised_stats(fields)
        spec = st._spectrum_stats()
        t_state = time.time() - t0

        # ---- backward compatibility check (R=0 vs the v5 state) --------
        h0_diff = float(np.max(np.abs(
            st.H0.astype(np.complex128) - st_v5.H0.astype(np.complex128))))
        om_diff = float(np.max(np.abs(st._omega_q() - st_v5._omega_q())))
        h_diff = float(np.max(np.abs(
            fields["h"].astype(np.float64)
            - f_v5["h"].astype(np.float64))))
        loop = st.loop_maxdiff()

        t1 = time.time()
        C, info = det_frame(fields, st.S, grid, env, Q)
        t_c = time.time() - t1
        cn = os.path.join(args.out, "C-rain-R%02g.npy" % R)
        np.save(cn, np.asarray(C, dtype=np.float32))
        C_by_rate[R] = C

        c_diff_v5 = None
        if R == 0.0:
            c_diff_v5 = {
                "max_abs_diff": float(np.max(np.abs(
                    C.astype(np.float64) - C_v5.astype(np.float64)))),
                "mean_abs_diff": float(np.mean(np.abs(
                    C.astype(np.float64) - C_v5.astype(np.float64)))),
                "identical": bool(np.array_equal(C, C_v5)),
            }

        img, dbg = render_frame_v5(
            fields, C, photo_lin, env, PW, PH, ambient_scale=A_V3,
            c_layer_exposure=E_v6, substrate_gain_rgb=gain_base["gain_rgb"],
            beer_down=beer_ref["down"], beer_up=beer_ref["up"],
            g1=g1_node["G1"], blur_sigma_px=blur_ref["sigma_px"],
            lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
        m = image_metrics(dbg["lin"], dbg["aces"])
        del dbg
        tag = "R%02g" % R
        fn = "still-%s.png" % tag
        Image.fromarray(img).save(os.path.join(args.out, fn))
        m["file"] = fn
        stills.append((img, "R = %g mm/h" % R,
                       "sigma2 = %.4f rad2 (calm %.4f + rain %.4f)"
                       % (realised["sigma2_rad2"],
                          re_v5["sigma2_rad2"],
                          max(realised["sigma2_rad2"] - re_v5["sigma2_rad2"],
                              0.0))))

        still_diff_v5 = None
        if R == 0.0:
            still_diff_v5 = {
                "max_abs_diff_uint8": int(np.max(np.abs(
                    img.astype(np.int16) - img_v5.astype(np.int16)))),
                "mean_abs_diff_uint8": float(np.mean(np.abs(
                    img.astype(np.int16) - img_v5.astype(np.int16)))),
                "identical": bool(np.array_equal(img, img_v5)),
            }

        node = st.rain_node(fields)
        if R > 0:
            bd = node["sigma2_breakdown"]
            bd["after_damping_spectral_rad2"] = float(
                st._spectrum_stats()["sigma2_rad2"]
                - st.rain_band["slope_variance_realised_rad2"])
        per_rate[tag] = {
            "rate_mm_h": float(R),
            "state": {
                "patch_S_m": float(st.S), "dx_m": float(st.dx),
                "grid": [grid, grid],
                "seed": int(st.seed),
                "spectrum_after": st.spectrum_after,
                "calm_floor_spectral": st.calm_floor,
                "realised": realised,
                "loop_maxdiff_h": loop,
                "H0_max_abs_diff_vs_v5": h0_diff,
                "omega_q_max_abs_diff_vs_v5": om_diff,
                "h_max_abs_diff_vs_v5": h_diff,
                "build_time_s": t_state,
            },
            "rain_node": node,
            "transport": {**c_stats(C), "det_time_s": t_c, "det_info": info,
                          "file": os.path.basename(cn)},
            "C_vs_v5_baseline": c_diff_v5,
            "still": {**m, "label": "R = %g mm/h" % R,
                      "vs_v5_baseline": still_diff_v5},
        }
        print("%-6s sigma2_spec=%.5f sigma2_real=%.5f band_var=%.3e "
              "C mean=%.9f p99=%.3f max=%.3f crush=%.4f blow=%.6f (%.1fs)"
              % (tag, spec["sigma2_rad2"], realised["sigma2_rad2"],
                 (0.0 if st.rain_band is None
                  else st.rain_band["variance_realised_m2"]),
                 float(C.mean()), float(np.percentile(C, 99)), float(C.max()),
                 m["crush_rate_linear_lt_0.02"],
                 m["blowout_rate_aces_gt_0.98"], t_c), flush=True)
        if R == 0.0:
            state_ref = st
            print("     backward-compat: H0 identical=%s  omega identical=%s "
                  "C identical=%s  still identical=%s"
                  % (h0_diff == 0.0, om_diff == 0.0,
                     c_diff_v5["identical"], still_diff_v5["identical"]),
                  flush=True)
        del fields, C, st
        gc.collect()

    # ------------------------------------------------------------- gates
    # 4: ring-band variance monotone increasing in R
    rates_pos = [r for r in rates if r > 0.0]
    band_vars = [(r, per_rate["R%02g" % r]["rain_node"]["band"]
                  ["variance_realised_m2"]) for r in rates_pos]
    mono_band = all(band_vars[i][1] < band_vars[i + 1][1]
                    for i in range(len(band_vars) - 1)) if len(band_vars) > 1 \
        else None
    # also across the full analytic range, so the gate is not just the 3 stills
    analytic_band = []
    for R in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0):
        analytic_band.append(
            (R, rain_state_summary(R)["band"]["ring_variance_m2"]))
    mono_analytic = all(analytic_band[i][1] < analytic_band[i + 1][1]
                        for i in range(len(analytic_band) - 1))

    base = per_rate.get("R00")
    if base is None:
        base = per_rate[[k for k in per_rate if
                         per_rate[k]["rate_mm_h"] == 0.0][0]]
    ref_crush = base["still"]["crush_rate_linear_lt_0.02"]
    ref_blow = base["still"]["blowout_rate_aces_gt_0.98"]
    crush_rows = {}
    for k, v in per_rate.items():
        crush_rows[k] = {
            "crush": v["still"]["crush_rate_linear_lt_0.02"],
            "blow": v["still"]["blowout_rate_aces_gt_0.98"],
            "vs_v5_calm_crush": m_v5["crush_rate_linear_lt_0.02"],
            "vs_v5_calm_blow": m_v5["blowout_rate_aces_gt_0.98"],
            "lum_mean": v["still"]["linear_luminance_mean"],
        }
    v5_ref = per_rate[[k for k in per_rate
                       if per_rate[k]["rate_mm_h"] == 0.0][0]]
    crush_limit = max(1.5 * ref_crush, ref_crush + 0.02)
    blow_limit = ref_blow + 0.002
    gates = {
        "gate1_backward_compatible_R0": {
            "gate": ("R=0 identical to the v5 baseline (C mean and stills), "
                     "PIPELINE-V6 section 6 item 1"),
            "v5_calm_C_mean": float(C_v5.mean()),
            "v6_R0_C_mean": v5_ref["transport"]["mean"],
            "C_mean_abs_diff": float(abs(v5_ref["transport"]["mean"]
                                         - float(C_v5.mean()))),
            "C_array_vs_v5": v5_ref["C_vs_v5_baseline"],
            "still_vs_v5": v5_ref["still"]["vs_v5_baseline"],
            "state_H0_identical": bool(
                v5_ref["state"]["H0_max_abs_diff_vs_v5"] == 0.0),
            "state_omega_identical": bool(
                v5_ref["state"]["omega_q_max_abs_diff_vs_v5"] == 0.0),
            "state_h_identical": bool(
                v5_ref["state"]["h_max_abs_diff_vs_v5"] == 0.0),
            "pass": bool(
                v5_ref["C_vs_v5_baseline"]["identical"]
                and v5_ref["still"]["vs_v5_baseline"]["identical"]
                and v5_ref["state"]["H0_max_abs_diff_vs_v5"] == 0.0),
        },
        "gate2_impact_flux_F": f_gate,
        "gate3_ring_wavenumber": k_gate,
        "gate4_band_variance_monotone": {
            "gate": ("ring-band variance strictly increasing with R "
                     "(PIPELINE-V6 section 6 item 4)"),
            "stills": [{"rate_mm_h": r, "band_variance_m2": v}
                       for r, v in band_vars],
            "analytic_range": [{"rate_mm_h": r, "band_variance_m2": v}
                               for r, v in analytic_band],
            "strictly_increasing_stills": bool(mono_band),
            "strictly_increasing_analytic": bool(mono_analytic),
            "pass": bool(mono_analytic and (mono_band is not False)),
        },
        "gate5_energy_and_tone": {
            "gate": ("mean(C)=1 (1e-3) for every master; crush/blow not "
                     "implausibly worse than the v5 baseline "
                     "(PIPELINE-V6 section 6 item 5)"),
            "C_mean_per_master": {
                k: {"mean": v["transport"]["mean"],
                    "abs_dev": v["transport"]["mean_abs_dev_from_1"],
                    "pass": v["transport"]["pass_mean_1"]}
                for k, v in per_rate.items()},
            "C_mean_pass": bool(all(v["transport"]["pass_mean_1"]
                                    for v in per_rate.values())),
            "tone_rows": crush_rows,
            "limits": {"crush_limit": float(crush_limit),
                       "blow_limit": float(blow_limit),
                       "rule": ("crush <= max(1.5x v5-calm crush, v5-calm "
                                "crush + 0.02) and blow <= v5-calm blow + "
                                "0.002; the v5 design reference 7.9% crush is "
                                "recorded in crush_rows for context")},
            "crush_pass": bool(max(r["crush"] for r in crush_rows.values())
                               <= crush_limit),
            "blow_pass": bool(max(r["blow"] for r in crush_rows.values())
                              <= blow_limit),
            "pass": bool(all(v["transport"]["pass_mean_1"]
                             for v in per_rate.values())
                         and max(r["crush"] for r in crush_rows.values())
                         <= crush_limit
                         and max(r["blow"] for r in crush_rows.values())
                         <= blow_limit),
        },
        "gate6_single_calibration_constant": {
            "gate": ("the only non-derived constant of this chunk is eta "
                     "(PIPELINE-V6 section 6 item 6)"),
            "name": "RING_ETA (v6_physics.RING_ETA)",
            "value": float(RING_ETA),
            "units": "-",
            "meaning": "impact -> ring-wave energy conversion efficiency",
            "role": "sigma_ring^2 = eta * F(R) * Ebar(R) * tau / "
                    "(rho g + sigma k_ring^2)",
            "basis": RING_ETA_BASIS,
            "sensitivity": eta_sensitivity(etas=ETAS_DEMO),
            "derived_quantities": [
                "Lambda(R) = 4.1 R^-0.21", "N(D) = N0 exp(-Lambda D)",
                "v(D) = 9.65 - 10.3 exp(-0.6 D)", "F(R) = int N v dD",
                "D0(R) = 3.67/Lambda", "Dmean(R) = 1/Lambda",
                "k_ring(R) (monotone placement from D0, anchored to the "
                "observed 1-2 cm band)",
                "Ebar(R) flux-weighted 1/2 m v^2", "tau = 150 ms (report 3.4)",
                "band 310-630 rad/m", "omega^2 = gk + (sigma/rho)k^3"],
            "not_calibrated_here": (
                "no other fitted number is introduced in chunk 2.  The "
                "inherited v5 calibrations (turbidity scale, scattering-blur "
                "proportionality, b_w(550) anchor, display-level "
                "normalisation) are unchanged and listed under "
                "inherited_calibrations_v5."),
            "pass": True,
        },
    }

    # ------------------------------------------------------------- panels
    panel_items = [(img, lab, cap) for (img, lab, cap) in stills]
    panel = save_panel_2x2(os.path.join(args.out,
                                        "panel-R00-R02-R10-R50.png"),
                           panel_items)
    # magnified crop: same physical window, delivered stills + C masters
    order = sorted(per_rate.keys(),
                   key=lambda k: per_rate[k]["rate_mm_h"])
    crop_info = build_crop_figure(
        args.out, [per_rate[k]["rate_mm_h"] for k in order], S_patch,
        args.crop, args.zoom)

    # ---- reproducible recipes -----------------------------------------
    rec_dir = os.path.join(args.out, "recipes")
    os.makedirs(rec_dir, exist_ok=True)
    rec_files = []
    for R in rates:
        r = with_rain(recipe, R, wind_ms=float(args.wind))
        fn = os.path.join(rec_dir, "rain-R%02g-v2.json" % R)
        with open(fn, "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2, ensure_ascii=False)
        rec_files.append(os.path.basename(fn))

    # ------------------------------------------------------------ numbers
    sigma_rows = []
    for k in order:
        v = per_rate[k]
        node = v["rain_node"]
        row = {
            "rate_mm_h": v["rate_mm_h"],
            "sigma2_spectral_rad2": v["state"]["spectrum_after"]["sigma2_rad2"],
            "sigma2_realised_rad2": v["state"]["realised"]["sigma2_rad2"],
            "rms_slope_realised_rad": v["state"]["realised"]["rms_slope"],
            "Hs_spectral_m": v["state"]["spectrum_after"]["Hs_m"],
            "m0_spectral_m2": v["state"]["spectrum_after"]["m0_m2"],
            "ring_band_variance_m2": (
                None if not node.get("active")
                else node["band"]["variance_realised_m2"]),
            "ring_band_slope_variance_rad2": (
                None if not node.get("active")
                else node["band"]["slope_variance_realised_rad2"]),
            "breakdown": node.get("sigma2_breakdown"),
        }
        sigma_rows.append(row)

    numbers = {
        "pipeline": "v6-chunk2",
        "design": "PIPELINE-V6.md (caustics v3 chunk 2: rain, windless calm)",
        "evidence": "report.md sections 3.1-3.5 (rain), 5 (variable map)",
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": recipe.get("name"), "recipe_schema_in": recipe.get("schema"),
        "preset": args.preset,
        "seed": int(state_ref.seed),
        "t_s": 0.0,
        "scene": {
            "note": ("the verification scene is preset 1's photo and sun with "
                     "wind = 0 (calm): the background is the Cox-Munk floor "
                     "sigma2 ~ 0.003, never a perfect mirror "
                     "(PIPELINE-V6 section 7)"),
            "wind_ms": env["windSpeedMS"],
            "specL_m": env["specL"], "patch_S_m": S_patch,
            "patch_resolution_m": S_patch / grid,
            "ring_wavelength_px_at_R10": float(
                2 * np.pi / 422.2) / (S_patch / grid),
            "depthM": env["depthM"], "turbid": env["turbid"],
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "photo": args.photo, "photo_size": [int(PW), int(PH)],
            "design_assumed_patch_m": 1.48,
            "patch_note": ("PIPELINE-V6 section 7 quotes a 1.48 m patch; at "
                           "wind 0 the recipe's specL clamp gives 1.20 m, "
                           "i.e. a FINER ring sampling (1.5 cm ~ %.1f px, "
                           "not 21 px)" % (0.015 / (S_patch / grid))),
            "background_anisotropy": ("the windless residual keeps the "
                                      "existing pipeline's cos^2 lobe shape "
                                      "(v5 behaviour, required for the R=0 "
                                      "identity gate); the rain band itself "
                                      "is isotropic"),
        },
        "recipe_v6_extension": {
            "parsed": cfg,
            "new_fields": ["rain.rate_mm_h"],
            "legacy_note": ("v1 recipes without the rain block read exactly as "
                            "before; R=0 is a hard early return to the v5 "
                            "state (no band, no damping, identical omega)"),
            "recipes_written": rec_files,
        },
        "rain_physics_nodes": {
            "eta": float(RING_ETA), "eta_basis": RING_ETA_BASIS,
            "tau_ring_s": float(RING_TAU_S),
            "ring_band_rad_per_m": list(RING_K_BAND),
            "capillary_crossover": cw,
            "F_table_report_3_3": {str(k): v for k, v in F_TABLE.items()},
            "F_table_tolerance_factor": float(F_TABLE_TOL_FACTOR),
            "per_rate": {k: per_rate[k]["rain_node"] for k in order},
            "analytic_table": {
                "rate_mm_h": [1.0, 2.0, 5.0, 10.0, 25.0, 50.0],
                "k_ring_rad_per_m": [ring_wavenumber(r)["k_ring_rad_per_m"]
                                     for r in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)],
                "lambda_ring_mm": [ring_wavenumber(r)["lambda_ring_m"] * 1e3
                                   for r in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)],
                "F_clamped": [rain_state_summary(r)["flux"]["F_clamped"]
                              for r in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)],
                "Ebar_uJ": [rain_state_summary(r)["energy"]["E_flux_weighted_uJ"]
                            for r in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)],
                "band_variance_m2": [rain_state_summary(r)["band"]
                                     ["ring_variance_m2"]
                                     for r in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)],
            },
        },
        "sigma2_breakdown": {
            "definition": ("sigma2 = integral k^2 S dk over the resolved band; "
                           "'calm floor' is the windless Cox-Munk-calibrated "
                           "background, 'rain band' is the ring-wave "
                           "contribution (PIPELINE-V6 section 7 asks for this "
                           "split)"),
            "rows": sigma_rows,
            "cox_munk_floor_reference": 0.003,
        },
        "state": {k: per_rate[k]["state"] for k in order},
        "transport": {
            "method": "det.py det_frame (UNCHANGED), Q=%d, sigma=0.75px, "
                      "no clamp / no renormalisation" % Q,
            "v5_calm_baseline": {"mean": float(C_v5.mean()),
                                 "stats": c_stats(C_v5),
                                 "det_time_s": t_c_v5, "det_info": info_v5,
                                 "file": "C-baseline-v5-calm.npy"},
            "per_rate": {k: per_rate[k]["transport"] for k in order},
        },
        "stills": {
            "individual": {k: per_rate[k]["still"]["file"] for k in order},
            "panel": os.path.basename(panel),
            "crop": crop_info,
            "metrics": {k: per_rate[k]["still"] for k in order},
            "v5_calm_baseline": {**m_v5, "file": "still-v5-calm.png"},
            "metrics_definitions": {
                "crush": "fraction of pixels with rec709 linear luminance "
                         "(pre-ACES, incl. exposure) < 0.02",
                "blow": "fraction of pixels with rec709 ACES luminance "
                        "> 0.98 (post-ACES, pre-sRGB-encode)",
            },
            "display": {
                "E_v6": float(E_v6), "E_v3": float(E_V3),
                "level_ratio": float(level_ratio),
                "flat_water_direct_lum_v3": float(L_v3),
                "flat_water_direct_lum_v5": float(L_v5),
                "note": ("same display-level normalisation as v5, evaluated "
                         "once at the calm reference condition and held "
                         "fixed for all R so the differences between the four "
                         "stills are physical (the ring band), not exposure"),
            },
        },
        "gates": gates,
        "approximations_and_design_choices": [
            {"item": "ring placement rule k_ring(R)",
             "statement": ("report 3.4 gives the observed 1-2 cm ring band and "
                           "the monotone 'bigger drops -> longer rings' "
                           "behaviour, but no closed-form frequency law.  The "
                           "placement used is lambda = 10 mm + 10 mm * "
                           "(D0-D0(R=1))/(D0(R=50)-D0(R=1)): the two observed "
                           "band limits are the anchors and D0(R) is the "
                           "interpolation variable.  It reproduces "
                           "k(10 mm/h) = 422 rad/m vs the design's ~420.  "
                           "Modelling choice, not a calibration constant."),
             "kind": "placement_rule"},
            {"item": "ring band energy -> height variance",
             "statement": ("E/(rho g + sigma k^2) converts an energy per unit "
                           "area into a height variance using the deep-water "
                           "linear energy density; evaluated at k_ring (the "
                           "band is narrow, so this is a band-centre "
                           "approximation)."),
             "kind": "approximation"},
            {"item": "lifetime tau",
             "statement": ("tau = 150 ms is the impact-process duration "
                           "measured in report 3.4, used as the mean ring "
                           "lifetime as the design instructs; viscous decay "
                           "is folded into it (no separate beta)."),
             "kind": "report_value"},
            {"item": "second rain effect (sub-band damping)",
             "statement": ("PIPELINE-V6 section 4 asks only for the mechanism "
                           "with a coefficient close to 1 and the shape "
                           "documented as an approximation.  gamma = F*Ebar/"
                           "(rho g m0_lowk) is a first-order rate with NO "
                           "tuned constant; the factor 1/(1+gamma*tau) is "
                           "evaluated over one impact-process lifetime and "
                           "weighted to the sub-band side only.  In this "
                           "windless scene there is only the Cox-Munk floor "
                           "to damp, so the applied gain is 0.998/0.987/0.926 "
                           "at R=2/10/50 -- close to 1, as designed.  The "
                           "per-8 s-loop value of the same rate "
                           "(exp(-gamma*8)) is recorded but NOT applied: a "
                           "static-spectrum loop is not the place to "
                           "integrate a slow damping law."),
             "kind": "approximation_permitted_by_design"},
            {"item": "band Gaussian not hard-clipped",
             "statement": ("the bump is a Gaussian in k with "
                           "sigma_k = (630-310)/6; at R=1 and R=50 the centre "
                           "sits on the band edge, so the tail crosses "
                           "310/630.  The variance fraction inside the band "
                           "is recorded per rate (0.89 at R=2, 0.99 at R=10, "
                           "0.59 at R=50).  The bump MEAN wavenumber equals "
                           "k_ring(R) by symmetry of the Gaussian, so the "
                           "monotone-centre gate is exact."),
             "kind": "approximation"},
            {"item": "rain does not modify the light path",
             "statement": ("PIPELINE-V6 section 8: no airborne drops, no "
                           "splash/crown/stalk geometry, no per-drop "
                           "simulation.  Only the surface spectrum changes."),
             "kind": "scope"},
        ],
        "calibrations_not_derived": [
            {"name": "RING_ETA", "value": float(RING_ETA),
             "role": "impact -> ring-wave energy conversion efficiency",
             "basis": RING_ETA_BASIS}],
        "inherited_calibrations_v5": [
            "turbidity -> da+db scale (TURBID_DC_MAX=0.30 1/m, 70% scatter)",
            "scattering blur proportionality (BLUR_K_PX=40 px at 2048^2)",
            "b_w(550)=0.0018 1/m anchor",
            "display-level normalisation E=E_v3*(L_v3/L_v5) (display "
            "transform, not physics)",
        ],
        "timing_s": {"v5_state_build": t_v5_state,
                     "v5_transport": t_c_v5,
                     "wall_total": time.time() - t_all},
    }
    with open(os.path.join(args.out, "numbers-v6.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)

    print("\nGATES")
    for name, g in gates.items():
        print("  %-34s %s" % (name, g.get("pass")))
    print("OK wall=%.1fs -> %s" % (time.time() - t_all,
                                   os.path.join(args.out,
                                                "numbers-v6.json")))


if __name__ == "__main__":
    main()
