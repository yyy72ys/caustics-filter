import os, shutil, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
V14 = os.path.join(HERE, "out", "v14")
W = os.path.join(V14, "_work")
src = os.path.join(W, "SNAPSHOT-run_v14_report.py")
dst = os.path.join(W, "run_v14_report.py.concurrent-writer")
if os.path.exists(src):
    shutil.move(src, dst)
    print("renamed builder snapshot ->", os.path.basename(dst),
          "(foreign implementation, kept as evidence)")
for f in ("_v14_done.py",):
    p = os.path.join(HERE, f)
    if os.path.exists(p):
        shutil.move(p, os.path.join(W, f))
shutil.copy2(os.path.join(V14, "REPORT-V14-GATES.md"),
             os.path.join(W, "SNAPSHOT-REPORT-V14-GATES.md"))
rep = os.path.join(V14, "REPORT-V14-GATES.md")
print("report", os.path.getsize(rep),
      hashlib.sha256(open(rep, "rb").read()).hexdigest()[:16])
print("\n_work contents:")
for f in sorted(os.listdir(W)):
    print("  " + f)
