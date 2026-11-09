"""Diagnostic: does the v5 physics still produce caustics at higher wind?

PIPELINE-V5's calmness calibration forces the realised mean square slope to
the Cox-Munk value.  At preset 1 (U = 1.5 m/s) that makes the surface very
calm, so the caustic contrast drops a lot.  This is the designed physical
correction (the v3 state was rougher than a 14 m/s ocean at 1.5 m/s wind).

This script is NOT part of the verification stills.  It re-runs the same v5
pipeline at higher wind speeds to show that the caustic contrast is restored
by the wind node (a scene input), not by weakening the physics.

Usage (from caustics-filter/offline/):
  python diag_v5_wind.py --winds 5,8
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
from det import det_frame  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert,
                        smith_g1_beckmann, scattering_blur_sigma_px,
                        substrate_gain)  # noqa: E402
from v5_state import V5WaterState  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,
                          image_metrics, v3_beer,
                          flat_water_direct_level)  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=25)
    ap.add_argument("--winds", default="5,8")
    ap.add_argument("--out", default="out/v5")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    args = ap.parse_args()
    winds = [float(x) for x in str(args.winds).split(",") if x.strip()]
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape
    recs = []
    for U in winds:
        r = json.loads(json.dumps(recipe))
        r["windWave"]["windSpeedMS"] = float(U)
        env = resolve_environment(r)
        t0 = time.time()
        st = V5WaterState(r, env, N=args.grid, preset_index=args.preset)
        fields = st.fields_at(0.0)
        fields["S"] = st.S
        real = st.realised_stats(fields)
        att = water_attenuation(env["turbid"])
        beer = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
        g1 = smith_g1_beckmann(real["rms_slope"],
                               90.0 - env["sun"]["altitudeDeg"])["G1"]
        blur = scattering_blur_sigma_px(att, env["depthM"], args.grid)["sigma_px"]
        L_v3 = flat_water_direct_level(env, v3_beer(env))
        L_v5 = flat_water_direct_level(env, beer["two_way"]) * g1
        E_v5 = E_V3 * (L_v3 / L_v5 if L_v5 > 0 else 1.0)
        C, info = det_frame(fields, st.S, args.grid, env, int(args.q))
        gain = substrate_gain(None)["gain_rgb"]
        img, dbg = render_frame_v5(fields, C, photo_lin, env, PW, PH,
                                   ambient_scale=A_V3, c_layer_exposure=E_v5,
                                   substrate_gain_rgb=gain,
                                   beer_down=beer["down"], beer_up=beer["up"],
                                   g1=g1, blur_sigma_px=blur,
                                   lowpass_sigma_px=B_SIGMA_V3,
                                   return_debug=True)
        fn = "diag-wind%g.png" % U
        Image.fromarray(img).save(os.path.join(args.out, fn))
        m = image_metrics(dbg["lin"], dbg["aces"])
        rec = {
            "wind_ms": float(U), "file": fn,
            "Hs_spec_m": st.spectrum_after["Hs_m"],
            "Hs_realised_m": real["Hs_m"],
            "sigma2_spectral": st.spectrum_after["sigma2_rad2"],
            "sigma2_realised": real["sigma2_rad2"],
            "rms_slope_deg": float(np.degrees(np.arctan(real["rms_slope"]))),
            "cox_munk_target": st.calmness_node(fields)["cox_munk_reference"],
            "calm_applied": st.calm["applied"] if st.calm else None,
            "calm_kc": (st.calm or {}).get("kc_rad_per_m"),
            "C_mean": float(np.mean(C)), "C_p99": float(np.percentile(C, 99)),
            "C_max": float(np.max(C)),
            "E_v5": float(E_v5), "blur_sigma_px": float(blur),
            "G1": float(g1), "metrics": m, "time_s": time.time() - t0,
        }
        recs.append(rec)
        print("U=%g Hs=%.3f sigma2=%.5f(real %.5f) Cmean=%.6f C p99=%.2f "
              "max=%.2f crush=%.4f blow=%.6f (%.0fs)"
              % (U, rec["Hs_realised_m"], rec["sigma2_spectral"],
                 rec["sigma2_realised"], rec["C_mean"], rec["C_p99"],
                 rec["C_max"], m["crush_rate_linear_lt_0.02"],
                 m["blowout_rate_aces_gt_0.98"], rec["time_s"]), flush=True)
        del fields, C, st
    with open(os.path.join(args.out, "diag-wind.json"), "w",
              encoding="utf-8") as f:
        json.dump({"purpose": "diagnostic only: caustic contrast vs wind "
                              "under the v5 Cox-Munk calibration (not one of "
                              "the four verification pairs)",
                   "recipe": recipe.get("name"), "preset": args.preset,
                   "Q": int(args.q), "grid": int(args.grid),
                   "runs": recs}, f, indent=2, ensure_ascii=False, default=float)


if __name__ == "__main__":
    main()
