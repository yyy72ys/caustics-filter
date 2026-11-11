import json

d = json.load(open("out/v8/numbers-v8.json", encoding="utf-8"))
pre = d["prerequisite_tasks"]
g = d["gates"]


def f(x, n=4):
    if x is None:
        return "n/a"
    return ("%%.%df" % n) % x


print("=== scene / metadata")
print(json.dumps(d["scene"], indent=1))
print(json.dumps(d["metadata_chunk4"], indent=1))
print("seed", d["seed"], "wall", d["timing_s"]["wall_total"])
print("states:", json.dumps(d["state_detail"], indent=1))
print("\n=== dispersion quantisation")
print(json.dumps(d["state"]["dispersion"]["quantisation"], indent=1))
print("\n=== anchor")
a = d["anchor"]
print("target", a["node"]["sigma2_target_rad2"], "old", a["old_vs_new"])
print(json.dumps(a["node"], indent=1)[:2500])
print("sensitivity:", json.dumps(a["sensitivity"], indent=1))
print("\n=== direction")
print(json.dumps(d["direction"]["node"], indent=1))
print(json.dumps(d["direction"]["comparison_vs_phillips"], indent=1)[:2000])
print("\n=== tail / breaking / patch")
print(json.dumps(d["capillary_tail"], indent=1))
print(json.dumps(d["breaking"], indent=1))
print("\n=== stills")
for k in ("preset1_v8", "preset1_v33", "shallow", "clamp"):
    print(k, json.dumps(d["stills"][k], indent=1))
print(json.dumps(d["stills"]["display"], indent=1))
print("\n=== transport")
print(json.dumps(d["transport"]["master"], indent=1))
print(json.dumps(d["transport"]["baseline_33"], indent=1)[:900])
print("\n=== PREREQ E")
print(json.dumps(pre["smith_shadow_audit"]["reference_scene"], indent=1))
print(json.dumps(pre["smith_shadow_audit"]["still_change_if_removed"], indent=1))
print(json.dumps(pre["smith_shadow_audit"]["g1_corner_table"], indent=0))
print("\n=== PREREQ Q")
q = pre["q_convergence"]
for k, v in q["per_Q"].items():
    print("Q", k, json.dumps(v, indent=0))
for p in q["pairs"]:
    print("pair", p["Qa"], p["Qb"], json.dumps(
        {kk: p[kk] for kk in ("quantile_rel", "mean_rel", "rms_diff",
                              "msssim_log1pC", "peak_luminance_abs_err")}, indent=0))
    print("   peak:", json.dumps(p["peak"], indent=0))
print("\n=== PREREQ M")
m = pre["photon_msssim_decomposition"]
print(json.dumps({k: m[k] for k in ("photons", "gamma", "quantiles_det",
                                    "quantiles_photon", "quantile_rel",
                                    "fwhm_det", "fwhm_photon",
                                    "mean_abs_diff_log1pC", "pass_msssim",
                                    "expected_mean_rgb", "measured_mean_rgb")},
                 indent=1))
print("msssim C:", json.dumps(m["msssim_log1pC"], indent=1))
print("decomp:", json.dumps(m["decomposition"], indent=1))
print("still:", json.dumps(m.get("still"), indent=1))
print("\n=== GATES")
for k, v in g.items():
    print(k, "pass=", v.get("pass"))
print(json.dumps(g["gate1_difference_attribution"], indent=1))
print(json.dumps(g["gate2_shallow_phase_speed"]["as_used_shallow_state"], indent=1))
print(json.dumps({kk: vv for kk, vv in g["gate2_shallow_phase_speed"].items()
                  if kk not in ("rows", "table", "as_used_shallow_state")}, indent=1))
print(json.dumps(g["gate3_new_anchor_sigma2"], indent=1))
print(json.dumps(g["gate4_breaking_clamp"], indent=1))
print(json.dumps(g["gate5_standard_preset_tone"], indent=1))
print(json.dumps({kk: vv for kk, vv in g["gate6_rain_path_untouched"].items()
                  if kk not in ("flux_gate_now", "F_table_now")}, indent=1))
print("=== report")
print(json.dumps(d["prerequisite_tasks"]["q_convergence"]["per_Q"], indent=1)[:400])
