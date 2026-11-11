import math
import numpy as np
import v11_spectrum as V11
import run_v14_spectrum as V14

G = 9.81
print("cell           fp     brief(m2/Hz)  Eomega  phillips  B_l_flank/ref  Bh/0.5am")
for U in (1.5, 3.0, 6.0, 10.0):
    for F in (100.0, 1000.0):
        ctx = V14.assemble(U, F, 1.0, arm="RAW")
        fp = float(ctx["laws"]["f_p_Hz"])
        f = np.linspace(1.3 * fp, 3.0 * fp, 400)
        om = 2 * math.pi * f
        k = om ** 2 / G
        s = V11.s_parts(k, ctx)
        Sf = s["S"] * (8 * math.pi ** 2 * f / G)
        ust = float(ctx["ustar"]["ustar_ms"])
        b1 = math.exp(float(np.mean(np.log(Sf) - np.log(2.5 * G * ust * f ** -4))))
        Somega = s["S"] * (2 * om / G)
        b2 = math.exp(float(np.mean(np.log(Somega) - np.log(2.5 * G * ust * om ** -4))))
        b3 = math.exp(float(np.mean(np.log(Sf * (2 * math.pi) ** 4 * f ** 5 / G ** 2))))
        # (a) gravity flank
        kp = float(ctx["k_p_rad_per_m"])
        kk = np.logspace(math.log10(kp / 2), math.log10(2 * kp), 4000)
        ss = V11.s_parts(kk, ctx)
        bl = float(np.max(ss["B_l"]))
        refa = 0.1433 * float(ctx["alpha_p"]) * float(ctx["gamma_bl"]) * 1.04
        # (c) capillary
        km = 369.6
        kc = np.logspace(math.log10(km / 2), math.log10(2 * km), 4000)
        sc = V11.s_parts(kc, ctx)
        bh = float(np.max(sc["B_h"]))
        refc = 0.5 * float(ctx["alpha_m_value"])
        print("U=%-4g F=%-6g %.3f  %.3e     %.5f  %.4f     %.3f          %.3f"
              % (U, F, fp, b1, b2, b3, bl / refa, bh / refc))
