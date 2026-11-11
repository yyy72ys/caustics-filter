import sys, json, os
sys.path.insert(0, r"D:\mutaro-work\jobkit")
from checkpoint import atomic_write_bytes

p = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples\checkpoint.json"
st = json.load(open(p, encoding="utf-8"))
seen, ids = set(), []
for i in st["completed_ids"]:
    if i not in seen:
        seen.add(i); ids.append(i)
st["completed_ids"] = ids
st["completed_units"] = len(ids)
st["last_unit"] = ids[-1]
out = json.dumps(st, ensure_ascii=False, indent=2) + "\n"
atomic_write_bytes(p, out.encode("utf-8"))
chk = json.load(open(p, encoding="utf-8"))
print("completed_units:", chk["completed_units"], "| unique ids:", len(set(chk["completed_ids"])))
print("last_unit:", chk["last_unit"])
