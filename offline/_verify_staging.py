import os, hashlib, shutil
SRC = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
DST = r"D:\yohey2026\dsh\muse\caustics-filter\delivery\offline-v9\samples-1s"
files = ["comp-eta100-1s.mp4", "comp-eta080-1s.mp4", "still-eta100.png",
         "still-eta080.png", "SAMPLE-REPORT.md"]
os.makedirs(DST, exist_ok=True)
shutil.copyfile(os.path.join(SRC, "SAMPLE-REPORT.md"), os.path.join(DST, "SAMPLE-REPORT.md"))

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

ok = True
for f in files:
    a, b = os.path.join(SRC, f), os.path.join(DST, f)
    ha, hb = sha(a), sha(b)
    same = ha == hb
    ok &= same
    print("%-22s %9d %9d  %s  %s" % (f, os.path.getsize(a), os.path.getsize(b),
                                     ha[:16], "MATCH" if same else "MISMATCH"))
print("ALL_MATCH", ok)
print("staged files:", sorted(os.listdir(DST)))
