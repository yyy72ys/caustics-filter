import json
p = "out/v8-smoke/numbers-v8.json"
d = json.load(open(p, encoding="utf-8"))
g1 = d["gates"]["gate1_difference_attribution"]
print("total:", json.dumps({k: g1["total_difference"][k] for k in
      ("delta_mean","delta_p50","delta_p99","delta_p99_9","delta_max",
       "rms_diff","mean_abs_diff")}, indent=1))
print("p99/99.9/max/mean rows:")
for k in ("mean","p50","p99","p99_9","max"):
    print(" ", k, json.dumps(g1["attribution"]["rows"][k], indent=None))
print("states:", json.dumps(g1["attribution"]["states"], indent=1))
print("clamp:", json.dumps(g1["clamp_and_normalisation"], indent=1))
print("\ngate3:", json.dumps({k: d["gates"]["gate3_new_anchor_sigma2"][k] for k in
      ("new_target_rad2","realised_rad2","relative_residual","pass")}, indent=1))
print("realised detail:", json.dumps(d["gates"]["gate3_new_anchor_sigma2"]["realised"], indent=1))
print("\ngate2 shallow as used:", json.dumps(d["gates"]["gate2_shallow_phase_speed"]["as_used_shallow_state"], indent=1)[:1200])
print("\nq conv pairs:", json.dumps(d["prerequisite_tasks"]["q_convergence"]["pairs"], indent=1)[:2500])
print("\nphoton:", json.dumps(d["prerequisite_tasks"]["photon_msssim_decomposition"]["decomposition"], indent=1))
print("\nstill photon:", json.dumps(d["prerequisite_tasks"]["photon_msssim_decomposition"].get("still"), indent=1))
print("\nrain:", json.dumps(d["gates"]["gate6_rain_path_untouched"], indent=1)[:1500])
print("\ndirection cmp:", json.dumps(d["direction"]["comparison_vs_phillips"], indent=1)[:1500])
