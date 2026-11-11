import math, json
import numpy as np
import v8_physics as P8
import v10_spectrum as V10
import v11_spectrum as V11
import run_v14_spectrum as V14

G = 9.81
for U in (1.5, 3.0, 6.0, 10.0):
    for F in (100.0, 1000.0):
        ctx = V14.assemble(U, F, 1.0, arm="RAW")
        fp = float(ctx["laws"]["f_p_Hz"])
        f = np.linspace(1.3 * fp, 3.0 * fp, 400)
        k = (2 * math.pi * f) ** 2 / G
        s = V11.s_parts(k, ctx)
        dkdf = 8 * math.pi ** 2 * f / G
        E = s["S"] * dkdf
        ust = float(ctx["ustar"]["ustar_ms"])
        ref = 2.5 * G * ust * f ** -4
        lb = float(np.mean(np.log(E) - np.log(ref)))
        beta = math.exp(lb)
        # contamination: B_h > B_l in the fit band?
        cont = bool(np.any(s["B_h"] > s["B_l"]))
        print("U=%-5g F=%-6g fp=%.3f Om=%.2f u*=%.4f beta=%.6f cont=%s Bmax=%.5f"
              % (U, F, fp, float(ctx["omega_c"]), ust, beta, cont,
                 float(np.max(s["B_l"] + s["B_h"]))))
