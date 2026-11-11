import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'out/v10/numbers-v10-transport.json'), encoding='utf-8'))
print(list(d.keys()))
sc = d.get('scenes', {})
print('scenes:', list(sc.keys()))
for k, v in sc.items():
    print('==', k, list(v.keys())[:30])
    print(json.dumps({kk: v.get(kk) for kk in ('name','U','F','h','v10_state','spectrum','realised') if kk in v}, indent=1)[:1800])
