"""Independent verification of the 1-second samples (no pipeline imports)."""
import glob
import hashlib
import json
import os
import subprocess

import numpy as np
from PIL import Image

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "v9samples")


def probe(p):
    cmd = ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
           "-show_entries",
           "stream=codec_name,width,height,nb_read_frames,r_frame_rate,"
           "avg_frame_rate,duration,pix_fmt", "-of", "json", p]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(r.stderr)
    return json.loads(r.stdout)["streams"][0]


def h(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


res = {"dir": OUT, "mp4": {}, "frames": {}}
for tag in ("eta100", "eta080"):
    p = os.path.join(OUT, "comp-%s-1s.mp4" % tag)
    info = probe(p)
    res["mp4"][tag] = {
        "bytes": os.path.getsize(p), "sha12": h(p), "ffprobe": info,
        "pass": bool(int(info["width"]) == 832 and int(info["height"]) == 480
                     and int(info["nb_read_frames"]) == 24
                     and abs(float(info["duration"]) - 1.0) < 1e-9)}
for tag in ("eta100", "eta080"):
    fs = sorted(glob.glob(os.path.join(OUT, "frames", tag, "frame-*.png")))
    hs = [h(f) for f in fs]
    arrs = [np.asarray(Image.open(f).convert("RGB")) for f in fs]
    shapes = sorted({a.shape for a in arrs})
    diffs = [float(np.abs(arrs[i].astype(np.int16)
                         - arrs[i + 1].astype(np.int16)).mean())
             for i in range(len(arrs) - 1)]
    res["frames"][tag] = {
        "count": len(fs), "unique_sha": len(set(hs)), "shapes": shapes,
        "identical_within_series": len(set(hs)) != len(fs),
        "consecutive_mean_abs_diff": {"min": min(diffs), "max": max(diffs),
                                      "mean": float(np.mean(diffs))},
        "frame_mean_uint8": {"min": float(min(a.mean() for a in arrs)),
                             "max": float(max(a.mean() for a in arrs))},
    }
a = np.asarray(Image.open(os.path.join(OUT, "still-eta100.png")).convert("RGB"))
b = np.asarray(Image.open(os.path.join(OUT, "still-eta080.png")).convert("RGB"))
res["stills"] = {"shape_A": list(a.shape), "shape_B": list(b.shape),
                 "identical": bool(h(os.path.join(OUT, "still-eta100.png"))
                                   == h(os.path.join(OUT, "still-eta080.png"))),
                 "mean_abs_diff": float(np.abs(a.astype(np.int16)
                                               - b.astype(np.int16)).mean()),
                 "max_abs_diff": int(np.abs(a.astype(np.int16)
                                            - b.astype(np.int16)).max())}
res["tmp_files"] = glob.glob(os.path.join(OUT, "**", "*.tmp"), recursive=True)
res["report_exists"] = os.path.exists(os.path.join(OUT, "SAMPLE-REPORT.md"))
res["report_bytes"] = (os.path.getsize(os.path.join(OUT, "SAMPLE-REPORT.md"))
                       if res["report_exists"] else None)
parts = [json.load(open(p, encoding="utf-8")) for p in
         sorted(glob.glob(os.path.join(OUT, "parts", "frame-*.json")))]
res["parts"] = {
    "count": len(parts),
    "frames": [p["frame"] for p in parts],
    "configs": sorted({(p["Q_direct"], p["n_dirs"], round(p["w_direct"], 6),
                        round(p["w_diffuse"], 6)) for p in parts}),
    "per_frame_wall_s": {"mean": float(np.mean([p["wall_s"]["total"]
                                                for p in parts])),
                         "sum_min": float(np.sum([p["wall_s"]["total"]
                                                  for p in parts]) / 60.0)},
    "A_vs_B_png_mean_abs_uint8": {
        "min": min(p["A_vs_B_png"]["mean_abs_uint8"] for p in parts),
        "max": max(p["A_vs_B_png"]["mean_abs_uint8"] for p in parts)},
}
res["pass"] = bool(all(v["pass"] for v in res["mp4"].values())
                   and all(v["count"] == 24 and not v["identical_within_series"]
                           and v["shapes"] == [(480, 832, 3)]
                           for v in res["frames"].values())
                   and not res["stills"]["identical"]
                   and not res["tmp_files"] and res["report_exists"])
with open(os.path.join(OUT, "verification.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, ensure_ascii=False)
print(json.dumps(res, indent=2, ensure_ascii=False))
