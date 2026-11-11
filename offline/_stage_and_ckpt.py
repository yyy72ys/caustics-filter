import sys, os, hashlib, shutil, json
sys.path.insert(0, r"D:\mutaro-work\jobkit")
from checkpoint import Checkpoint

SRC = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"
DST = r"D:\yohey2026\dsh\muse\caustics-filter\delivery\offline-v9\samples-1s"
files = ["comp-eta100-1s.mp4", "comp-eta080-1s.mp4", "still-eta100.png",
         "still-eta080.png", "SAMPLE-REPORT.md"]
os.makedirs(DST, exist_ok=True)
for f in files:
    shutil.copyfile(os.path.join(SRC, f), os.path.join(DST, f))

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

checks = {}
for f in files:
    checks[f] = sha(os.path.join(SRC, f)) == sha(os.path.join(DST, f))
assert all(checks.values()), checks

cp = Checkpoint(os.path.join(SRC, "checkpoint.json"))
cp.unit_done("delivery-staged (delivery/offline-v9/samples-1s, 5 files byte-identical; no scp route exists)")
print(json.dumps({"staged": checks, "units": cp.state()["completed_units"],
                  "last": cp.state()["last_unit"]}, indent=2))
