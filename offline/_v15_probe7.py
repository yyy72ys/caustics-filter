import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'out/v10/numbers-v10-transport.json'), encoding='utf-8'))
a = d['scenes']['A']['v10_state']
print('KEYS:', list(a.keys()))
for k in a:
    if k in ('realised', 'spectrum_after', 'spectral_vs_spatial', 'quantisation', 'amplitude'):
        print('==', k)
        print(json.dumps(a[k], indent=1)[:1200])
