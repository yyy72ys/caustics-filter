import hashlib, os
DST = r"D:\yohey2026\dsh\muse\caustics-filter\delivery\offline-v9\samples-1s"
for f in ("comp-eta100-1s.mp4","comp-eta080-1s.mp4","still-eta100.png","still-eta080.png"):
    p = os.path.join(DST,f)
    print(f, hashlib.sha256(open(p,"rb").read()).hexdigest())
