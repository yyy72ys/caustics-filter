import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()
for p in ['out/v10/still-A-v37.png', 'out/v10/still-B-v37.png', 'out/v10/comp-windy-v37.png']:
    fp = os.path.join(HERE, p)
    print(sha(fp)[:24], os.path.getsize(fp), p)
print('--- v10 numbers keys')
d = json.load(open(os.path.join(HERE, 'out/v10/numbers-v10-sweep.json'), encoding='utf-8'))
print(list(d.keys()))
print(json.dumps(d, indent=1)[:1500])
