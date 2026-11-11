import json
import time
import numpy as np

from env import resolve_environment
from recipe_v8 import parse_recipe_v8, with_v8
from v6_state import V6WaterState
from v8_state import V8WaterState
import v8_physics as P

with open("../blender/recipe-1.json", "r", encoding="utf-8") as f:
    raw = json.load(f)

env = resolve_environment(raw)
cfg = parse_recipe_v8(raw)
print("fetch parse (legacy):", cfg["fetch"])
print("env specL=%.4f S=%.4f depth=%.2f U=%.2f" % (env["specL"], 4 * env["specL"],
                                                   env["depthM"], env["windSpeedMS"]))

print("\nanchor:", json.dumps(P.anchor_sigma2(1.5, 100.0, 0.6,
                                             S_patch=4 * env["specL"]), indent=1))
print("\nlegacy fetch:", json.dumps(P.legacy_fetch_m(1.5, 0.6,
                                                    S_patch=4 * env["specL"]), indent=1))
print("\nsens:", json.dumps(P.anchor_sensitivity(1.5, 100.0, 0.6,
                                                S_patch=4 * env["specL"]), indent=1))
print("\ndispersion gate h=0.1:", json.dumps(P.dispersion_gate_table(0.1), indent=1))
print("\ntail:", json.dumps(P.capillary_tail(1.4755, 2048, 0.6), indent=1))
print("\ndirection:", json.dumps(P.direction_comparison(
    1.5, 100.0, 0.6, env["specL"], 2 * np.pi / (4 * env["specL"]),
    np.pi * 512 / (4 * env["specL"])), indent=1))
print("\nbreaking preset:", json.dumps(P.breaking_limits(1.1969, 0.6), indent=1))
print("breaking h=0.1:", json.dumps(P.breaking_limits(1.1969, 0.1), indent=1))

t0 = time.time()
st0 = V6WaterState(raw, env, N=512, preset_index=1)
f0 = st0.fields_at(0.0)
print("\nbaseline(512) build %.1fs sigma2=%.6f Hs=%.4f" % (
    time.time() - t0, st0.spectrum_after["sigma2_rad2"], st0.spectrum_after["Hs_m"]))

t0 = time.time()
st = V8WaterState(raw, env, N=512, preset_index=1, fetch_m=100.0)
f = st.fields_at(0.0)
print("v8(512) build %.1fs" % (time.time() - t0))
print("  anchor applied=%s target=%.6f before=%.6f after=%.6f kc=%.2f" % (
    st.anchor["applied"], st.anchor["sigma2_target_rad2"],
    st.anchor["sigma2_before"], st.anchor["sigma2_after_spectral"],
    st.anchor["kc_rad_per_m"]))
print("  direction s_max=%.2f sigma2 %.6f -> %.6f" % (
    st.direction["s_max"], st.direction["sigma2_before"],
    st.direction["sigma2_after"]))
print("  breaking=%s" % json.dumps(st.breaking))
print("  realised:", st.realised_stats(f))
print("  quant mean rel err:", st.dispersion_node()["quantisation"]["slope_weighted_mean_rel_error"])
print("  spectrum_after:", st.spectrum_after["sigma2_rad2"], st.spectrum_after["Hs_m"])

# t=0 field equality between dispersion branches
st2 = V8WaterState(raw, env, N=512, preset_index=1, fetch_m=100.0,
                   use_finite_depth=False)
f2 = st2.fields_at(0.0)
print("  t0 h identical across dispersion branches:",
      float(np.max(np.abs(f["h"] - f2["h"]))),
      "H0:", float(np.max(np.abs(st.H0 - st2.H0))))
print("  omega max diff:", float(np.max(np.abs(st._omega_q() - st2._omega_q()))))
f_2s = st.fields_at(2.0)
f2_2s = st2.fields_at(2.0)
print("  t2 h max diff:", float(np.max(np.abs(f_2s["h"] - f2_2s["h"]))),
      "std h:", float(np.std(f_2s["h"])))
print("  loop maxdiff v8:", st.loop_maxdiff())
