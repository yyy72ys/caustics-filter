import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
for name in ['numbers-v12.json', 'numbers-v12-g1.json', 'numbers-v12-g4.json']:
    p = os.path.join(HERE, 'out/v12', name)
    d = json.load(open(p, encoding='utf-8'))
    print('==', name, list(d.keys()))

def walk(o, path=''):
    if isinstance(o, dict):
        for k, v in o.items():
            if isinstance(v, str) and len(v) == 64 and all(c in '0123456789abcdef' for c in v):
                print('  HASH', path + '/' + k, v[:16])
            else:
                walk(v, path + '/' + k)
    elif isinstance(o, list):
        for i, v in enumerate(o[:6]):
            walk(v, path + '[%d]' % i)

for name in ['numbers-v12.json', 'numbers-v12-g1.json', 'numbers-v12-g4.json']:
    d = json.load(open(os.path.join(HERE, 'out/v12', name), encoding='utf-8'))
    walk(d, name)
