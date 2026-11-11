import os, sys, json
paths = [
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\run_v14_spectrum.py',
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\out\v14\REPORT-V14-GATES.md',
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\out\v14\numbers-v14-spectrum.json',
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\out\v10',
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\run_v12.py',
 r'D:\yohey2026\dsh\muse\caustics-filter\offline\v12_water.py',
]
out = []
allok = True
for p in paths:
    ok = os.path.exists(p)
    allok = allok and ok
    line = ('OK   ' if ok else 'MISS ') + p
    if ok and os.path.isdir(p):
        line += '  (dir, %d entries)' % len(os.listdir(p))
    out.append(line)
    print(line)
print('CWD', os.getcwd())
print('PREFLIGHT', 'PASS' if allok else 'FAIL')
sys.exit(0 if allok else 1)
