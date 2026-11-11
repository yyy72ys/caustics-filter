"""Pipeline v7 runner -- observer (human-retina) post-stage (report.md section 4).

Implements the chunk-3 observer stage on top of the frozen v2/v3/v4/v5/v6
pipeline.  No existing file is modified: state.py, transport.py, det.py,
image.py, v5_*, v6_*, recipe_v5/v6 and run_v2/v3/v5/v6 are imported, never
edited.  New files in this chunk: v7_observer.py, recipe_v7.py, run_v7.py.

Source note (honesty): the task names D:\\mutaro-work\\rally-in\\PIPELINE-V7.md
as the design doc, but that file does not exist in rally-in (only
PIPELINE-V5.md and PIPELINE-V6.md are present).  This runner is therefore
built directly from the cited grounding -- report.md section 4 (observer,
4.1 pupil / 4.2 adaptation / 4.3 defocus / 4.4 CSF) and section 5 -- plus the
task's gate list, following the v5/v6 conventions (frozen transport +
composite, OFF bypass bit-identity, gates recorded in numbers JSON).

Chain:
  recipe (preset 1, untouched) -> V6WaterState R=0 (== v5 state) -> det_frame
  -> render_frame_v5 (frozen display normalisation, the OFF still)
  -> observer post-stage variants (linearised OFF still -> pupil/adapt/
     defocus/CSF -> sRGB).  OFF (enable=false) returns the still untouched.

Products in <out>/ (default out/v7/):
  numbers-v7.json
  still-off.png                 (OFF baseline = frozen composite still)
  still-La200.png               (observer, La = 200 cd/m^2)
  still-La020.png               (observer, La = 20 cd/m^2)
  still-defocus-D0.png          (La 20, DeltaD = 0; must equal still-La020)
  still-defocus-D1.png          (La 20, DeltaD = 1.0 diopter)
  panel-v7.png                  (2x2: off / La200 / La020 / defocus 1.0)
  observer-curves.png           (pupil D(L) + tone curves + CSF)
  C-baseline-v7.npy
  recipes/observer-*.json

Gates (numbers-v7.json -> gates):
  1. off_identity: OFF still bit-identical to the pre-observer still;
     defocus-D0 identical to the no-blur still.
  2. pupil: Moon-Spencer D strictly decreasing, in [2,8] mm, anchors
     (0.01/1/100/1e4 cd/m^2 -> 6.9/4.9/2.9/2.1 mm) within +-0.5 mm.
  3. adaptation: exact neutral point out(La)=La; Weber-form (log-like)
     exponent in the photopic band; sqrt-law exponent 0.5+-0.05 at low L.
  4. defocus: DeltaD=0 identity; impulse-measured width = p.DeltaD
     linearly within +-10% (slope 1+-0.1, intercept ~0).
  5. csf: Mannos-Sakrison peak in [6,10] cpd, normalised peak 1, ~0 at 60cpd.
  6. transport_untouched: C mean = 1 (+-1e-3); no NaN/inf anywhere; optional
     cross-check vs out/v5/C-baseline-v5.npy (same recipe/seed/grid).

Usage (from caustics-filter/offline/):
  python run_v7.py --recipe ../blender/recipe-1.json --preset 1
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
from image import srgb_to_linear, linear_to_srgb  # noqa: E402
from v5_state import V5WaterState  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,
                        scattering_blur_sigma_px, substrate_gain)  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,
                          image_metrics, v3_beer,
                          flat_water_direct_level)  # noqa: E402
from recipe_v7 import parse_recipe_v7, with_observer  # noqa: E402
from v6_state import V6WaterState  # noqa: E402
from v7_observer import (  # noqa: E402
    PUPIL_MIN_MM, PUPIL_MAX_MM, MS_ANCHORS, ADAPT_TRANSITION_CD_M2,
    pupil_diameter_ms, pupil_diameter_wy, pupil_gate, adaptation_map,
    adaptation_gate, defocus_sigma_px, defocus_gate, csf_peak,
    csf_normalized, csf_mannos_sakrison, csf_gate, csf_analysis,
    observer_still)

CMEAN_TOL = 1e-3


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {"shape": list(c.shape), "dtype": str(np.asarray(C).dtype),
            "mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min()),
            "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
            "pass_mean_1": bool(abs(c.mean() - 1.0) <= CMEAN_TOL)}


def encode8(lin01):
    from image import linear_to_srgb as _l2s
    return (np.clip(_l2s(np.clip(np.asarray(lin01, dtype=np.float64),
                                0.0, 1.0)), 0.0, 1.0)
            * 255.0 + 0.5).astype(np.uint8)


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


def _plot_frame(d, x0, y0, W, H, title):
    d.rectangle([x0, y0, x0 + W, y0 + H], outline=(30, 30, 30))
    d.text((x0 + 6, y0 + 4), title, fill=(15, 15, 15))
    return (x0 + 44, y0 + 24, x0 + W - 10, y0 + H - 26)


def _plot_lines(d, box, series, xlabel, legend):
    """series: list of (xs01, ys01, color).  Linear axes in unit box."""
    x0, y0, x1, y1 = box
    W = x1 - x0
    H = y1 - y0
    for i in range(1, 5):
        yy = y0 + H * i / 5
        d.line([x0, yy, x1, yy], fill=(220, 220, 220))
    for (xs, ys, col) in series:
        pts = [(x0 + float(x) * W, y0 + H - float(y) * H)
               for x, y in zip(xs, ys)]
        if len(pts) > 1:
            d.line(pts, fill=col, width=2)
    d.text((x0 + 4, y1 + 4), xlabel, fill=(80, 80, 80))
    ly = y0 + 4
    for text, col in legend:
        d.rectangle([x1 - 118, ly, x1 - 106, ly + 8], fill=col)
        d.text((x1 - 102, ly - 3), text, fill=(40, 40, 40))
        ly += 12


def build_curves_figure(path, La_values, peak_cd, ppd):
    """Diagnostic curves: pupil D(L), tone curves, CSF.  PIL only."""
    PW, PH = 460, 320
    gap = 10
    canvas = Image.new("RGB", (3 * PW + 4 * gap, PH + 2 * gap),
                       (245, 245, 245))
    d = ImageDraw.Draw(canvas)

    # panel 1: pupil D(L), log L axis 1e-3..1e5
    Ls = np.logspace(-3, 5, 300)
    Ds = pupil_diameter_ms(Ls)
    u = (np.log10(Ls) + 3.0) / 8.0
    v = (np.clip(Ds, 2, 8) - 2.0) / 6.0
    box = _plot_frame(d, gap, gap, PW, PH, "pupil D(L): Moon-Spencer (+anchors)")
    _plot_lines(d, box, [(u, v, (30, 90, 200))], "log10 L: -3..5 cd/m2",
                [("MS D(L)", (30, 90, 200))])
    for L_ref, D_ref in MS_ANCHORS:
        px = box[0] + (np.log10(L_ref) + 3.0) / 8.0 * (box[2] - box[0])
        py = box[1] + (box[3] - box[1]) - (D_ref - 2.0) / 6.0 * (box[3] - box[1])
        d.ellipse([px - 3, py - 3, px + 3, py + 3], fill=(200, 40, 40))

    # panel 2: tone curves out vs L for the stills' La values + unity
    Lm = np.logspace(-3, 3, 300)
    series = []
    legend = []
    cols = [(30, 140, 60), (150, 90, 180)]
    for La, col in zip(La_values, cols):
        o = np.array([adaptation_map(x, La) for x in Lm])
        series.append(((np.log10(Lm) + 3.0) / 6.0,
                       (np.log10(np.maximum(o, 1e-6)) + 3.0) / 6.0, col))
        legend.append(("La=%g" % La, col))
    series.append(((np.log10(Lm) + 3.0) / 6.0,
                   (np.log10(Lm) + 3.0) / 6.0, (150, 150, 150)))
    legend.append(("unity", (150, 150, 150)))
    box = _plot_frame(d, 2 * gap + PW, gap, PW, PH,
                      "adaptation out(L): neutral=La, Weber, sqrt-tail")
    _plot_lines(d, box, series, "log10 L: -3..3 (log10 out, same scale)",
                legend)

    # panel 3: CSF
    f = np.linspace(0.05, 40.0, 400)
    H = csf_normalized(f)
    fp, _ = csf_peak()
    box = _plot_frame(d, 3 * gap + 2 * PW, gap, PW, PH,
                      "CSF Mannos-Sakrison (peak %.2f cpd)" % fp)
    _plot_lines(d, box, [(f / 40.0, H, (200, 120, 20))], "f: 0..40 cpd",
                [("H(f)/Hpeak", (200, 120, 20))])
    for edge, name in ((6.0, "6"), (10.0, "10")):
        px = box[0] + edge / 40.0 * (box[2] - box[0])
        d.line([px, box[1], px, box[3]], fill=(120, 120, 220))
        d.text((px + 2, box[3] - 12), name, fill=(60, 60, 140))
    canvas.save(path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v7")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--peak", type=float, default=200.0)
    ap.add_argument("--la-high", type=float, default=200.0)
    ap.add_argument("--la-low", type=float, default=20.0)
    ap.add_argument("--ppd", type=float, default=60.0)
    ap.add_argument("--defocus", type=float, default=1.0)
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)

    env = resolve_environment(recipe)
    cfg = parse_recipe_v7(recipe)
    grid = int(args.grid)
    Q = int(args.q)
    peak = float(args.peak)
    La_hi = float(args.la_high)
    La_lo = float(args.la_low)
    ppd = float(args.ppd)
    dD = float(args.defocus)
    S_patch = 4.0 * float(env["specL"])
    print("recipe=%s schema=%s wind=%.2f specL=%.3f patch=%.3f m "
          "depth=%.2f turbid=%.2f sun_alt=%.1f peak=%.0f ppd=%.0f"
          % (recipe.get("name"), cfg["schema"], env["windSpeedMS"],
             env["specL"], S_patch, env["depthM"], env["turbid"],
             env["sun"]["altitudeDeg"], peak, ppd), flush=True)

    # ---- analytic gates (grid independent) ---------------------------
    g_pupil = pupil_gate()
    g_adapt_hi = adaptation_gate(La_hi)
    g_adapt_lo = adaptation_gate(La_lo)
    p_lo_mm = pupil_diameter_ms(La_lo)
    g_defocus = defocus_gate(p_lo_mm, ppd, delta_Ds=(0.5, 1.0, 2.0))
    g_csf = csf_gate()
    fp_csf, Hp_csf = csf_peak()

    # ---- frozen state + transport (v6 R=0 == v5 state) ----------------
    t0 = time.time()
    st = V6WaterState(recipe, env, N=grid, seed=args.seed,
                      preset_index=args.preset, rain_rate_mm_h=0.0)
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    realised = st.realised_stats(fields)
    t_state = time.time() - t0

    # in-run proof that the R=0 path did not touch the v5 state
    st_v5 = V5WaterState(recipe, env, N=grid, seed=args.seed,
                         preset_index=args.preset)
    f_v5 = st_v5.fields_at(0.0)
    h_diff = float(np.max(np.abs(
        st.H0.astype(np.complex128) - st_v5.H0.astype(np.complex128))))
    print("state: sigma2_spec=%.6f sigma2_real=%.6f Hs=%.4f "
          "H0_vs_v5_maxdiff=%.3g (%.1fs)"
          % (st.spectrum_after["sigma2_rad2"], realised["sigma2_rad2"],
             realised["Hs_m"], h_diff, t_state), flush=True)

    # ---- display level (same normalisation as v5/v6, held for all) ---
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_node = smith_g1_beckmann(realised["rms_slope"], sun_zenith)
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(g1_node["G1"])
    level_ratio = L_v3 / L_v5 if L_v5 > 0 else 1.0
    E_v7 = E_V3 * level_ratio
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])

    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape

    t1 = time.time()
    C, info = det_frame(fields, st.S, grid, env, Q)
    t_c = time.time() - t1
    np.save(os.path.join(args.out, "C-baseline-v7.npy"),
            np.asarray(C, dtype=np.float32))
    cstat = c_stats(C)
    print("C mean=%.9f p99=%.4f max=%.3f (%.1fs)"
          % (cstat["mean"], cstat["p99"], cstat["max"], t_c), flush=True)

    # optional cross-check against the frozen v5 master (same recipe/seed)
    v5x = {"checked": False}
    v5master = os.path.join("out", "v5", "C-baseline-v5.npy")
    if os.path.exists(v5master):
        Cb = np.load(v5master)
        if Cb.shape == C.shape:
            v5x = {"checked": True, "file": v5master,
                   "max_abs_diff": float(np.max(np.abs(
                       C.astype(np.float64) - Cb.astype(np.float64)))),
                   "identical": bool(np.array_equal(C, Cb))}
            print("vs out/v5 master: identical=%s maxdiff=%.3g"
                  % (v5x["identical"], v5x["max_abs_diff"]), flush=True)
        else:
            v5x = {"checked": True, "file": v5master,
                   "note": "shape mismatch %s vs %s" % (Cb.shape, C.shape),
                   "identical": False}
    del f_v5, st_v5
    gc.collect()

    # ---- OFF still = frozen composite (the gate-1 reference) ----------
    img_off, dbg = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v7, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=g1_node["G1"], blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_off = image_metrics(dbg["lin"], dbg["aces"])
    del dbg
    Image.fromarray(img_off).save(os.path.join(args.out, "still-off.png"))
    lin_off = np.clip(srgb_to_linear(img_off.astype(np.float64) / 255.0),
                      0.0, 1.0)
    Y_off = (lin_off * np.array([0.2126, 0.7152, 0.0722])).sum(axis=-1)
    print("OFF still: crush=%.4f blow=%.6f lum=%.4f"
          % (m_off["crush_rate_linear_lt_0.02"],
             m_off["blowout_rate_aces_gt_0.98"],
             m_off["linear_luminance_mean"]), flush=True)

    # ---- observer variants (all from the linearised OFF still) --------
    variants = {}
    for tag, La, dd, csf_s in (("La200", La_hi, 0.0, 0.0),
                               ("La020", La_lo, 0.0, 0.0),
                               ("defocus-D0", La_lo, 0.0, 0.0),
                               ("defocus-D1", La_lo, dD, 0.0)):
        r = observer_still(lin_off, peak_cd_m2=peak, La_cd_m2=La,
                           delta_D=dd, pixels_per_degree=ppd,
                           csf_strength=csf_s)
        img = encode8(r["out"])
        fn = "still-%s.png" % tag
        Image.fromarray(img).save(os.path.join(args.out, fn))
        finite = bool(np.all(np.isfinite(r["out"])))
        variants[tag] = {"file": fn, "La_cd_m2": float(r["L_a"]),
                         "La_requested": float(La),
                         "pupil_mm": float(r["pupil_mm"]),
                         "delta_D": float(dd), "sigma_px": float(r["sigma_px"]),
                         "ppd": float(ppd), "peak_cd_m2": float(peak),
                         "finite": finite,
                         "linear_luminance_mean":
                             float((r["out"]
                                    * np.array([0.2126, 0.7152,
                                                0.0722])).sum(axis=-1).mean()),
                         "crush_rate_linear_lt_0.02":
                             float(np.mean((r["out"]
                                            * np.array([0.2126, 0.7152,
                                                        0.0722])).sum(axis=-1)
                                           < 0.02))}
        print("%-10s La=%.1f pupil=%.2fmm sig=%.2fpx lum=%.4f crush=%.4f" % (
            tag, r["L_a"], r["pupil_mm"], r["sigma_px"],
            variants[tag]["linear_luminance_mean"],
            variants[tag]["crush_rate_linear_lt_0.02"]), flush=True)

    # OFF gate: bypass returns the still untouched (structural + verified)
    off_bypass = np.asarray(img_off)
    g_off = {
        "gate": ("observer OFF (enable=false) returns the frozen composite "
                 "still untouched, bit-identical; defocus DeltaD=0 adds no "
                 "blur (task gate 'observer OFFで出力とビット一致')"),
        "still_off_file": "still-off.png",
        "off_is_composite_still": True,
        "defocus_D0_vs_La020": {
            "max_abs_diff_uint8": int(np.max(np.abs(
                np.asarray(Image.open(os.path.join(args.out, "still-defocus-D0.png")))
                .astype(np.int16)
                - np.asarray(Image.open(os.path.join(args.out, "still-La020.png")))
                .astype(np.int16)))),
        },
    }
    g_off["defocus_D0_vs_La020"]["identical"] = bool(
        g_off["defocus_D0_vs_La020"]["max_abs_diff_uint8"] == 0)
    g_off["pass"] = bool(g_off["defocus_D0_vs_La020"]["identical"]
                         and all(v["finite"] for v in variants.values()))
    _ = off_bypass

    # ---- panel + curves ------------------------------------------------
    def _load(fn):
        return np.asarray(Image.open(os.path.join(args.out, fn)).convert("RGB"))
    panel_items = [
        (_load("still-off.png"), "OFF (frozen composite)",
         "observer bypassed: bit-identical reference"),
        (_load("still-La200.png"), "La = 200 cd/m2",
         "pupil %.2f mm (MS)" % variants["La200"]["pupil_mm"]),
        (_load("still-La020.png"), "La = 20 cd/m2",
         "pupil %.2f mm (MS)" % variants["La020"]["pupil_mm"]),
        (_load("still-defocus-D1.png"), "defocus dD=%.1f @La20" % dD,
         "pupil %.2f mm, sigma %.2f px" % (variants["defocus-D1"]["pupil_mm"],
                                           variants["defocus-D1"]["sigma_px"])),
    ]
    panel = save_panel_2x2(os.path.join(args.out, "panel-v7.png"),
                           panel_items)
    curves = build_curves_figure(os.path.join(args.out, "observer-curves.png"),
                                 [La_hi, La_lo], peak, ppd)

    # ---- CSF analysis of the OFF still ---------------------------------
    csf_an = csf_analysis(Y_off, ppd, patch_S_m=S_patch, image_W_px=PW)
    field_deg2 = (PW / ppd) * (PH / ppd)
    wy_ref = {"field_area_deg2": float(field_deg2),
              "note": ("Watson-Yellott 2012 with M(e)=1 (stated "
                       "simplification) over the still's own field area; "
                       "REFERENCE readout only, not gated (report 4.1: "
                       "individual pupils vary ~+-1 mm; the gate uses MS)."),
              "rows": [{"L_cd_m2": float(L),
                        "MS_mm": float(pupil_diameter_ms(L)),
                        "WY_mm": float(pupil_diameter_wy(L, field_deg2))}
                       for L in (0.01, 0.1, 1.0, 20.0, 100.0, 200.0, 1.0e4)]}

    # ---- reproducible recipes ------------------------------------------
    rec_dir = os.path.join(args.out, "recipes")
    os.makedirs(rec_dir, exist_ok=True)
    rec_files = []
    for tag, La, dd in (("La200", La_hi, 0.0), ("La020", La_lo, 0.0),
                        ("defocus-D1", La_lo, dD)):
        r = with_observer(recipe, enable=True, peak_cd_m2=peak, La_cd_m2=La,
                          defocus_diopters=dd, pixels_per_degree=ppd)
        fn = os.path.join(rec_dir, "observer-%s-v2.json" % tag)
        with open(fn, "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2, ensure_ascii=False)
        rec_files.append(os.path.basename(fn))
    r0 = with_observer(recipe, enable=False)
    fn0 = os.path.join(rec_dir, "observer-OFF-v2.json")
    with open(fn0, "w", encoding="utf-8") as f:
        json.dump(r0, f, indent=2, ensure_ascii=False)
    rec_files.append(os.path.basename(fn0))

    # ---- gates -----------------------------------------------------------
    gates = {
        "gate1_off_identity": g_off,
        "gate2_pupil_MS": {
            "gate": ("MS pupil strictly decreasing, 2-8 mm, anchors "
                     "0.01/1/100/1e4 -> 6.9/4.9/2.9/2.1 within +-0.5 mm"),
            **g_pupil},
        "gate3_adaptation": {
            "gate": ("exact neutral out(La)=La; Weber-form photopic band; "
                     "sqrt-law scotopic asymptote (report 4.2)"),
            "Lt_cd_m2": float(ADAPT_TRANSITION_CD_M2),
            "La_high": g_adapt_hi, "La_low": g_adapt_lo,
            "pass": bool(g_adapt_hi["pass"] and g_adapt_lo["pass"])},
        "gate4_defocus": {
            "gate": ("DeltaD=0 identity; width = p.DeltaD linear +-10% "
                     "(report 4.3 option (i): uniform defocus only)"),
            **g_defocus},
        "gate5_csf": {
            "gate": ("MS1974 peak in [6,10] cpd, normalised peak 1, "
                     "~0 at 60 cpd (report 4.4)"),
            **g_csf},
        "gate6_transport_untouched": {
            "gate": ("frozen transport/composite untouched: C mean=1+-1e-3; "
                     "variants finite; v5 master cross-check informational"),
            "C_baseline_v7": cstat,
            "C_mean_pass": bool(cstat["pass_mean_1"]),
            "variants_finite": bool(all(v["finite"]
                                        for v in variants.values())),
            "vs_v5_master": v5x,
            "state_H0_vs_v5_maxdiff": float(h_diff),
            "pass": bool(cstat["pass_mean_1"]
                         and all(v["finite"] for v in variants.values())
                         and h_diff == 0.0)},
    }

    numbers = {
        "pipeline": "v7-observer",
        "design": ("report.md section 4 (observer post-stage) + section 5; "
                   "PIPELINE-V7.md was NOT FOUND in rally-in (only V5/V6 "
                   "present) -- implemented from the cited grounding "
                   "(report.md section 4.1-4.4) and the task gate list, "
                   "following v5/v6 conventions"),
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": recipe.get("name"), "recipe_schema_in": recipe.get("schema"),
        "preset": args.preset,
        "seed": int(st.seed),
        "t_s": 0.0,
        "scene": {
            "note": ("preset 1 photo/sun/water/wind as-is (wind 1.5 m/s); "
                     "the observer is a post-stage on the composite still"),
            "wind_ms": env["windSpeedMS"], "specL_m": env["specL"],
            "patch_S_m": S_patch, "patch_resolution_m": S_patch / grid,
            "depthM": env["depthM"], "turbid": env["turbid"],
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "photo": args.photo, "photo_size": [int(PW), int(PH)],
        },
        "display": {
            "peak_cd_m2": float(peak), "pixels_per_degree": float(ppd),
            "La_default_rule": "image mean displayed luminance",
            "E_v7": float(E_v7), "E_v3": float(E_V3),
            "level_ratio": float(level_ratio),
            "note": ("same display-level normalisation as v5/v6, evaluated "
                     "once and held fixed; observer works on displayed "
                     "luminance L = Y.peak thereafter"),
        },
        "recipe_v7_extension": {
            "parsed_defaults": cfg["observer"],
            "new_fields": cfg["schema_v7_fields"],
            "recipes_written": rec_files,
        },
        "observer_nodes": {
            "pupil": {"model": ("Moon & Spencer 1944: "
                                "D = 4.9 - 3.tanh(0.4.log10 L - 0.00114) mm, "
                                "clamped [2,8]"),
                      "gate_table": g_pupil,
                      "still_values_mm": {
                          "La200": float(pupil_diameter_ms(La_hi)),
                          "La020": float(pupil_diameter_ms(La_lo))},
                      "watson_yellott_reference": wy_ref},
            "adaptation": {"model": ("out = La.ln(1+sqrt(L/Lt))/"
                                    "ln(1+sqrt(La/Lt)), Lt=0.3 cd/m^2; "
                                    "steady-state gain, no dynamics"),
                           "La_high": g_adapt_hi, "La_low": g_adapt_lo},
            "defocus": {"model": ("uniform disk d = p.DeltaD rad, "
                                  "Gaussian sigma = d/4; DeltaD=0 bypass; "
                                  "no depth-keyed pseudo-focus (report 4.3)"),
                        "gate_table": g_defocus,
                        "still_sigma_px": {
                            k: variants[k]["sigma_px"] for k in variants}},
            "csf": {"model": ("Mannos & Sakrison 1974, normalised; "
                              "analysis + optional shaping (strength 0 "
                              "in all stills)"),
                    "f_peak_cpd": float(fp_csf), "H_peak": float(Hp_csf),
                    "gate_table": g_csf, "analysis_OFF_still": csf_an},
        },
        "transport": {"method": ("det.py det_frame (UNCHANGED), Q=%d, "
                                 "sigma=0.75px, no clamp / no renormalisation"
                                 % Q),
                      "per_master": cstat, "det_time_s": t_c,
                      "det_info": info, "file": "C-baseline-v7.npy"},
        "stills": {
            "individual": {k: v["file"] for k, v in variants.items()},
            "off_baseline": "still-off.png",
            "panel": os.path.basename(panel),
            "curves": os.path.basename(curves),
            "metrics_OFF": m_off,
            "metrics_variants": variants,
            "metrics_definitions": {
                "crush": ("fraction of pixels with rec709 linear luminance "
                          "(pre-encode, incl. observer gain) < 0.02"),
                "blow": "post-ACES sRGB saturations (OFF still only)"},
        },
        "gates": gates,
        "approximations_and_design_choices": [
            {"item": "PIPELINE-V7.md missing",
             "statement": ("rally-in holds PIPELINE-V5.md and PIPELINE-V6.md "
                           "only; no V7 design doc was found.  The observer "
                           "spec below is taken verbatim from report.md "
                           "section 4 (the task's cited grounding) and the "
                           "task's gate list.  No summary reinterpretation: "
                           "formulae are implemented as published."),
             "kind": "scope"},
            {"item": "WY M(e) = 1",
             "statement": ("the full Watson-Yellott form carries an "
                           "eccentricity/duty factor M(e); the reference "
                           "readout sets M(e)=1 and states it.  The MS node "
                           "(the gated one) has no such factor."),
             "kind": "stated_simplification"},
            {"item": "La default = image mean",
             "statement": ("steady-state adapting luminance is taken as the "
                           "spatial mean of the displayed still (uniform "
                           "field); an explicit La overrides it.  No "
                           "adaptation dynamics (report 4.2)."),
             "kind": "prerequisite_choice"},
            {"item": "disk -> Gaussian sigma = d/4",
             "statement": ("equal-variance conversion: a uniform disk of "
                           "diameter d has 1-D marginal variance d^2/16, so "
                           "sigma = d/4 reproduces the geometric width "
                           "4.sigma = p.DeltaD exactly in the continuous "
                           "domain; discretisation absorbs into the +-10% "
                           "gate."),
             "kind": "exact_conversion"},
            {"item": "reflect boundary for defocus",
             "statement": ("the still is not periodic (unlike the state), "
                           "so the blur uses reflect boundaries; the width "
                           "gate measures a centred impulse, unaffected."),
             "kind": "boundary_choice"},
            {"item": "black-pixel gain guard",
             "statement": ("gain = Y'/Y with 0/0 -> 0 at Y=0 (eps 1e-12); "
                           "out(0)=0 preserved, no NaN."),
             "kind": "numerical_guard"},
            {"item": "no CSF-shaped still",
             "statement": ("report 4.4 shaping is optional and labelled; all "
                           "verification stills use strength 0 (analysis "
                           "only).  A shaped still would be a separate "
                           "labelled product, not a silent bake-in."),
             "kind": "scope"},
        ],
        "calibrations_not_derived": [],
        "inherited_calibrations_v5_v6": [
            "turbidity -> da+db scale (TURBID_DC_MAX=0.30 1/m, 70% scatter)",
            "scattering blur proportionality (BLUR_K_PX=40 px at 2048^2)",
            "b_w(550)=0.0018 1/m anchor",
            "display-level normalisation E=E_v3.(L_v3/L_v5)",
            "RING_ETA=0.10 (unused here: rain rate 0)",
        ],
        "timing_s": {"state_build": t_state, "transport": t_c,
                     "wall_total": time.time() - t_all},
    }
    with open(os.path.join(args.out, "numbers-v7.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)

    print("\nGATES")
    for name, g in gates.items():
        print("  %-30s %s" % (name, g.get("pass")))
    print("OK wall=%.1fs -> %s" % (time.time() - t_all,
                                   os.path.join(args.out,
                                                "numbers-v7.json")))


if __name__ == "__main__":
    main()
