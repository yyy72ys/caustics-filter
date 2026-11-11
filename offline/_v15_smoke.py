import json
import run_v15_spectrum as V15

ev = V15.evaluate_cell(6.0, 100.0, 1.0)
print("state", ev["state"], ev["fail_reasons"])
print("prov", json.dumps(ev["provenance"], indent=1))
print("a", json.dumps({k: ev["conformance"]["a_gravity"][k] for k in
                       ("model_B_l_max", "reference_B_l_kp", "ratio", "pass")}, indent=1))
print("c", json.dumps({k: ev["conformance"]["c_capillary"][k] for k in
                       ("model_B_h_max", "reference_half_alpha_m_used", "ratio", "pass")}, indent=1))
print("b", json.dumps({k: ev["observation"]["b_equilibrium_beta"][k] for k in
                       ("beta", "in_band", "contaminated", "status",
                        "diagnostic_standard_phillips_beta")}, indent=1))
print("d", json.dumps({k: ev["observation"]["d_jonswap_anchor"][k] for k in
                       ("alpha_J", "B_J_kp", "model_B_l_kp", "caveat")}, indent=1))
print("hs", ev["hs"]["delta_rel"], "flag", ev["hs_flag"])
print("C1", ev["C1"]["pass"], "C1b", ev["C1b"]["pass"], "tot", ev["envelope_check"]["tol_inside"])
