import numpy as np
from PIL import Image
OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
pairs = [("eta100", "still-eta100.png"), ("eta080", "still-eta080.png")]
for tag, still in pairs:
    d = np.asarray(Image.open(OUT + "\\_dec-%s-f13.png" % tag).convert("RGB")).astype(np.float64)
    s = np.asarray(Image.open(OUT + "\\" + still).convert("RGB")).astype(np.float64)
    print("%s decoded frame13 mean %.3f | still mean %.3f | mean|d| %.3f | max %.0f"
          % (tag, d.mean(), s.mean(), np.abs(d - s).mean(), np.abs(d - s).max()))
