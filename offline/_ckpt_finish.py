import sys, json
sys.path.insert(0, r"D:\mutaro-work\jobkit")
from checkpoint import Checkpoint

OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
cp = Checkpoint(OUT + r"\checkpoint.json")
print("before:", json.dumps(cp.state(), ensure_ascii=False)[:200])
for uid in ("verify (24x2 PNGs + ffprobe + stills diff)",
            "encode-nvenc (h264_nvenc mp4s + ffprobe 24fr/1.000s)",
            "report (SAMPLE-REPORT.md, dome=48 recorded)"):
    cp.unit_done(uid)
st = cp.state()
print("completed_units:", st["completed_units"])
print("last_unit:", st["last_unit"])
print("updated_at:", st["updated_at"])
