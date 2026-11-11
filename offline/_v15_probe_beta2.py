import math
import numpy as np
import v11_spectrum as V11
import run_v14_spectrum as V14

G = 9.81
print("cell        fp     u*      brief_beta  phillips_beta  mean_k3S   B_l(kp)  B_h_max")
for U in (1.5, 3.0, 6.0, 10.0):
    for F in (100.0, 1000.0):
        ctx = V14.assemble(U, F, 1.0, arm="RAW")
        fp = float(ctx["laws"]["f_p_Hz"])
        f = np.linspace(1.3 * fp, 3.0 * fp, 400)
        k = (2 * math.pi * f) ** 2 / G
        s = V11.s_parts(k, ctx)
        dkdf = 8 * math.pi ** 2 * f / G
        Sf = s["S"] * dkdf
        ust = float(ctx["ustar"]["ustar_ms"])
        beta_brief = math.exp(float(np.mean(np.log(Sf) - np.log(2.5 * G * ust * f ** -4))))
        beta_phil = math.exp(float(np.mean(np.log(Sf * (2 * math.pi) ** 4 * f ** 5 / G ** 2))))
        # mean curvature in the fit band
        k3S = s["B_l"] + s["B_h"]
        print("U=%-4g F=%-6g %.3f %.4f  %.3e   %.4f      %.4f   %.5f  %.5f"
              % (U, F, fp, ust, beta_brief, beta_phil, float(np.mean(k3S)),
                 float(np.interp(float(ctx["k_p_rad_per_m"]),
                                 V11.s_parts(np.linspace(0.9*ctx["k_p_rad_per_m"], 1.1*ctx["k_p_rad_per_m"], 50), ctx)["k"],
                                 V11.s_parts(np.linspace(0.9*ctx["k_p_rad_per_m"], 1.1*ctx["k_p_rad_per_m"], 50), ctx)["B_l"])),
                 float(np.max(s["B_h"]))))
