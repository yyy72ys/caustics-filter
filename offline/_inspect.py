import json

def load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)

v7 = load("out/v7/numbers-v7.json")
print("V7 keys:", list(v7.keys()))
print("V7 transport:", json.dumps(v7["transport"]["per_master"], indent=1)[:1200])
print("V7 scene:", json.dumps(v7["scene"], indent=1))
print("V7 display:", json.dumps(v7["display"], indent=1))
print("V7 stills metrics OFF:", json.dumps(v7["stills"]["metrics_OFF"], indent=1))
print("V7 gate6:", json.dumps(v7["gates"]["gate6_transport_untouched"], indent=1)[:800])
print("V7 seed:", v7["seed"], "Q", v7["transport"]["det_info"]["Q"])

v6 = load("out/v6/numbers-v6.json")
print("\nV6 transport v5 baseline:", json.dumps(v6["transport"]["v5_calm_baseline"], indent=1)[:900])
print("V6 gate5 limits:", json.dumps(v6["gates"]["gate5_energy_and_tone"]["limits"], indent=1))
print("V6 gate5 crush rows:", json.dumps(v6["gates"]["gate5_energy_and_tone"]["tone_rows"], indent=1))
print("V6 state R00:", json.dumps(v6["state"]["R00"], indent=1)[:1500])
print("V6 still metrics R00:", json.dumps(v6["stills"]["metrics"]["R00"], indent=1)[:900])
print("V6 v5 calm baseline metrics:", json.dumps(v6["stills"]["v5_calm_baseline"], indent=1))
print("V6 display:", json.dumps(v6["stills"]["display"], indent=1))
