import json
d = json.load(open("out/v8-smoke/numbers-v8.json", encoding="utf-8"))
p = d["prerequisite_tasks"]["photon_msssim_decomposition"]
print(json.dumps(p["decomposition"]["peak_position"], indent=1))
print(json.dumps(p.get("still"), indent=1))
g1 = d["gates"]["gate1_difference_attribution"]
print(json.dumps(g1["accounting"], indent=1))
print(json.dumps(g1["attribution"]["rows"]["p99"], indent=1))
print("gates pass:", {k: v.get("pass") for k, v in d["gates"].items()})
