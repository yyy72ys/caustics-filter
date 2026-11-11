import os, json, glob, hashlib
import numpy as np
from PIL import Image

OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
res = {}

def sha12(p):
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return h[:12]

for tag in ("eta100", "eta080"):
    d = os.path.join(OUT, "frames", tag)
    files = sorted(glob.glob(os.path.join(d, "frame-*.png")))
    dims = []
    bad = []
    for f in files:
        im = Image.open(f)
        dims.append((im.size, im.mode))
        im.load()
    uniq = sorted(set(dims))
    res[tag] = {"count": len(files), "first": os.path.basename(files[0]) if files else None,
                "last": os.path.basename(files[-1]) if files else None,
                "sizes": [list(u[0]) for u in uniq], "modes": [u[1] for u in uniq],
                "expected_name_seq_ok": [os.path.basename(f) for f in files] ==
                                        ["frame-%04d.png" % i for i in range(1, 25)]}

parts = []
for i in range(62, 86):
    p = os.path.join(OUT, "parts", "frame-%04d.json" % i)
    j = json.load(open(p))
    assert j["frame"] == i, (p, j["frame"])
    parts.append(j)
res["parts_frames"] = [p["frame"] for p in parts]
res["n_dirs_values"] = sorted(set(p["n_dirs"] for p in parts))
res["Q_direct_values"] = sorted(set(p["Q_direct"] for p in parts))

w = [p["wall_s"] for p in parts]
a_only = sum(x["fields"] + x["direct"] + x["composite_A"] for x in w)
b_extra = sum(x["sky"] + x["composite_B"] for x in w)
total = sum(x["total"] for x in w)
res["wall"] = {
    "sum_total_s": total, "sum_total_min": total / 60.0,
    "A_only_s": a_only, "A_only_min": a_only / 60.0,
    "B_extra_s": b_extra, "B_extra_min": b_extra / 60.0,
    "A_cumulative_s": a_only, "B_cumulative_s": a_only + b_extra,
    "mean_s_per_frame": total / len(parts),
    "mean_A_s_per_frame": a_only / len(parts),
    "mean_B_s_per_frame": (a_only + b_extra) / len(parts),
    "max_frame_total_s": max(x["total"] for x in w),
    "min_frame_total_s": min(x["total"] for x in w),
}
res["crush_blow"] = {
    "A_crush_mean": float(np.mean([p["still_A"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "B_crush_mean": float(np.mean([p["still_B"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "A_blow_mean": float(np.mean([p["still_A"]["blowout_rate_aces_gt_0.98"] for p in parts])),
    "B_blow_mean": float(np.mean([p["still_B"]["blowout_rate_aces_gt_0.98"] for p in parts])),
    "A_crush_min": float(np.min([p["still_A"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "A_crush_max": float(np.max([p["still_A"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "B_crush_min": float(np.min([p["still_B"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "B_crush_max": float(np.max([p["still_B"]["crush_rate_linear_lt_0.02"] for p in parts])),
    "A_blow_max": float(np.max([p["still_A"]["blowout_rate_aces_gt_0.98"] for p in parts])),
    "B_blow_max": float(np.max([p["still_B"]["blowout_rate_aces_gt_0.98"] for p in parts])),
}
res["png_AB"] = {
    "mean_abs_uint8_all_frames": float(np.mean([p["A_vs_B_png"]["mean_abs_uint8"] for p in parts])),
    "max_abs_uint8_max": int(max(p["A_vs_B_png"]["max_abs_uint8"] for p in parts)),
}

# stills f0074 == frame-0013 (74-62+1 = 13)
idx = 74 - 62 + 1
sa = Image.open(os.path.join(OUT, "still-eta100.png"))
sb = Image.open(os.path.join(OUT, "still-eta080.png"))
fa = Image.open(os.path.join(OUT, "frames", "eta100", "frame-%04d.png" % idx))
fb = Image.open(os.path.join(OUT, "frames", "eta080", "frame-%04d.png" % idx))
aa = np.asarray(sa.convert("RGB")).astype(np.int16)
ab = np.asarray(sb.convert("RGB")).astype(np.int16)
d = np.abs(aa - ab)
res["stills"] = {
    "frame_index": idx,
    "still_A_size": list(sa.size), "still_B_size": list(sb.size),
    "sha12_A": sha12(os.path.join(OUT, "still-eta100.png")),
    "sha12_B": sha12(os.path.join(OUT, "still-eta080.png")),
    "identical_bytes": open(os.path.join(OUT, "still-eta100.png"), "rb").read() ==
                       open(os.path.join(OUT, "still-eta080.png"), "rb").read(),
    "mean_abs_diff_uint8": float(d.mean()),
    "max_abs_diff_uint8": int(d.max()),
    "pct_pixels_gt2": float((d.max(axis=2) > 2).mean() * 100),
    "stillA_matches_frame0013": bool(np.array_equal(aa, np.asarray(fa.convert("RGB")).astype(np.int16))),
    "stillB_matches_frame0013": bool(np.array_equal(ab, np.asarray(fb.convert("RGB")).astype(np.int16))),
}
res["mp4_bytes"] = {f: os.path.getsize(os.path.join(OUT, f)) for f in
                    ("comp-eta100-1s.mp4", "comp-eta080-1s.mp4")}
res["mp4_sha12"] = {f: sha12(os.path.join(OUT, f)) for f in
                    ("comp-eta100-1s.mp4", "comp-eta080-1s.mp4")}

print(json.dumps(res, indent=2))
json.dump(res, open(os.path.join(OUT, "_verify_now.json"), "w"), indent=2)
