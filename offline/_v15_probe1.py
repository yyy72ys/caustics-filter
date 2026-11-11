import json, sys
p = sys.argv[1] if len(sys.argv) > 1 else 'out/v14/numbers-v14-spectrum.json'
d = json.load(open(p, encoding='utf-8'))
print('KEYS:', list(d.keys()))
print(json.dumps(d.get('policy'), indent=1)[:2500])
c = d['cells'][0]
print('CELL KEYS:', list(c.keys()))
print('ARM KEYS:', list(c['arms']['RAW'].keys()))
print(json.dumps(c['arms']['RAW']['b_envelope'], indent=1))
