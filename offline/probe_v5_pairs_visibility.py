"""Visibility A/B probe for the v5 verification panels.

For each panel it renders the still TWICE -- once with the real master ratio C
and once with a flat C = 1 -- and reports the distribution of the per-pixel
ratio of the two sRGB luminances.  That ratio is the perceptually relevant
quantity: a p5..p95 band of 0.94..1.06 is invisible against photo texture,
0.78..1.25 is an obvious light net.

Reason it exists: "does the panel show the effect" is otherwise argued from
opinion.  The first v5 pair set's left panels (z = 0.5 m / sand / wind mode)
measured 0.94..1.06 and read as the untouched photo, which is why the pair
set was rebuilt (run_v5_pairs_fix.py).

Prints only; writes no files.  Run from caustics-filter/offline/:
  python probe_v5_pairs_visibility.py
"""
import gc
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,
                        substrate_gain, scattering_blur_sigma_px)  # noqa: E402
from v5_state import V5WaterState  # noqa: E402
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5, v3_beer,
                          flat_water_direct_level)  # noqa: E402

REC709 = np.array([0.2126, 0.7152, 0.0722])
OUT = "out/v5"
recipe = json.load(open("../blender/recipe-1.json", encoding="utf-8"))
env = resolve_environment(recipe)
att = water_attenuation(env["turbid"])
photo = np.asarray(Image.open("../test-inputs/moth-on-bark.jpg").convert("RGB"),
                   dtype=np.float64) / 255.0
photo_lin = srgb_to_linear(photo)
PH, PW, _ = photo_lin.shape

st = V5WaterState(recipe, env, N=2048, seed=None, preset_index=1)
fields = st.fields_at(0.0)
fields["S"] = st.S
r = st.realised_stats(fields)
g1 = smith_g1_beckmann(r["rms_slope"], 90 - env["sun"]["altitudeDeg"])["G1"]
beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"], env["turbid"])
blur_ref = scattering_blur_sigma_px(att, env["depthM"], 2048)
E_v5 = E_V3 * flat_water_direct_level(env, v3_beer(env)) \
    / (flat_water_direct_level(env, beer_ref["two_way"]) * g1)
gain_sand = substrate_gain("sand", 1.0)["gain_rgb"]


def srgb_lum(img):
    return (img.astype(np.float64) / 255.0 * REC709).sum(-1)


def ab(tag, C, f, e, gain, beer, blur):
    kw = dict(ambient_scale=A_V3, c_layer_exposure=E_v5, substrate_gain_rgb=gain,
              beer_down=beer["down"], beer_up=beer["up"], g1=g1,
              blur_sigma_px=blur, lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    a, _ = render_frame_v5(f, C, photo_lin, e, PW, PH, **kw)
    b, _ = render_frame_v5(f, np.ones_like(np.asarray(C, np.float32)),
                           photo_lin, e, PW, PH, **kw)
    la, lb = srgb_lum(a), srgb_lum(b)
    ratio = la / np.maximum(lb, 1e-6)
    print("%-24s sRGB ratio p1=%.3f p5=%.3f p50=%.3f p95=%.3f p99=%.3f | "
          "frac|dL|>0.05: %.1f%%  max|dL|=%.3f"
          % (tag, *np.percentile(ratio, [1, 5, 50, 95, 99]),
             100 * np.mean(np.abs(la - lb) > 0.05), np.abs(la - lb).max()))


for z in (0.5, 1.5, 3.0, 4.0):
    C = np.load(os.path.join(OUT, "C-depth-z%g.npy" % z))
    ez = json.loads(json.dumps(env))
    ez["depthM"] = float(z)
    ab("depth z=%.1f" % z, C, fields, ez, gain_sand,
       beer_lambert(z, env["sun"]["altitudeDeg"], env["turbid"]),
       scattering_blur_sigma_px(att, z, 2048)["sigma_px"])
C = np.load(os.path.join(OUT, "C-baseline-v5.npy"))
ab("ref 0.6 sand", C, fields, env, gain_sand, beer_ref, blur_ref["sigma_px"])

del fields
gc.collect()
for hs, npy in ((0.12, "C-wave-fixed0.12.npy"), (0.30, "C-wave-fixed0.30.npy")):
    s = V5WaterState(recipe, env, N=2048, seed=None, preset_index=1,
                     height_mode="fixed", height_m=hs)
    f = s.fields_at(0.0)
    f["S"] = s.S
    ab("wave fixed %.2f" % hs, np.load(os.path.join(OUT, npy)), f, env,
       gain_sand, beer_ref, blur_ref["sigma_px"])
    del s, f
    gc.collect()
