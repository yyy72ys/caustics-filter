import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'out/v10/numbers-windy-v37.json'), encoding='utf-8'))
print(json.dumps(d.get('task'), indent=1))
print(json.dumps(d.get('v10_state'), indent=1)[:2500])
print('--- files', d.get('files'))
