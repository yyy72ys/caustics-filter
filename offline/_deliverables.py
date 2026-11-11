import os, json, hashlib, subprocess
OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
files = ["comp-eta100-1s.mp4", "comp-eta080-1s.mp4", "still-eta100.png",
         "still-eta080.png", "SAMPLE-REPORT.md"]
for f in files:
    p = os.path.join(OUT, f)
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    print("%-22s %9d  %s  %s" % (f, os.path.getsize(p), h[:16], "OK" if os.path.exists(p) else "MISSING"))
