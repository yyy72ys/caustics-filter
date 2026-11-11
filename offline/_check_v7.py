import json
import numpy as np

C = np.load("out/v7/C-baseline-v7.npy")
c = np.asarray(C, np.float64)
print(C.dtype, C.shape)
print("mean", repr(float(c.mean())), "p99", repr(float(np.percentile(c, 99))),
      "p99.9", repr(float(np.percentile(c, 99.9))), "max", repr(float(c.max())))
d = json.load(open("out/v7/numbers-v7.json", encoding="utf-8"))
p = d["transport"]["per_master"]
print("v7 ", repr(p["mean"]), repr(p["p99"]), repr(p["p99_9"]), repr(p["max"]))
print("diffs", float(c.mean()) - p["mean"],
      float(np.percentile(c, 99)) - p["p99"],
      float(c.max()) - p["max"])
