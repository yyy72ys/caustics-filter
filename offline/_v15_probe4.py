import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'out/v10/numbers-windy-v37.json'), encoding='utf-8'))
print(list(d.keys()))
sc = d.get('scenes', {})
print(list(sc.keys()))
for k, v in sc.items():
    print('==', k, list(v.keys()))
    for kk in ('state_v10', 'v10', 'state', 'spectrum', 'realised'):
        if kk in v:
            print('  ', kk, json.dumps(v[kk], indent=1)[:800])
