import json, os
OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
keys = {}
tot = 0.0
for i in range(62, 86):
    j = json.load(open(os.path.join(OUT, "parts", "frame-%04d.json" % i)))
    w = j["wall_s"]
    tot += w["total"]
    for k, v in w.items():
        if k == "total":
            continue
        keys[k] = keys.get(k, 0.0) + v
print("sum of components:")
for k, v in keys.items():
    print("  %-14s %8.2f s  %6.2f min" % (k, v, v/60))
s = sum(keys.values())
print("components sum %.2f, total sum %.2f, residual %.2f" % (s, tot, tot - s))
