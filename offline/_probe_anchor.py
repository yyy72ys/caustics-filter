import numpy as np

G = 9.81
GAMMA = 0.0728
RHO = 1000.0
TF = GAMMA / RHO  # 7.28e-5 m^3/s^2

def omega(k, h):
    return np.sqrt((G * k + TF * k ** 3) * np.tanh(k * h))

def domegadk(k, h):
    w = omega(k, h)
    th = np.tanh(k * h)
    sch = 1.0 / np.cosh(np.clip(k * h, 0, 30)) ** 2
    num = (G + 3 * TF * k ** 2) * th + (G * k + TF * k ** 3) * h * sch
    return num / (2.0 * w)

def jonswap_sigma2(U, F, h, gammas=(3.3,), kmax=6000.0, n=400000):
    Ft = G * F / U ** 2
    fp = 3.5 * G / U * Ft ** (-1.0 / 3.0)
    alpha = 0.076 * Ft ** (-0.22)
    out = {}
    for gam in gammas:
        k = np.logspace(-3, np.log10(kmax), n)
        f = omega(k, h) / (2 * np.pi)
        sf = 0.07 if True else None
        sj = np.where(f <= fp, 0.07, 0.09)
        r = np.exp(-((f - fp) ** 2) / (2 * sj ** 2 * fp ** 2))
        Sf = (alpha * G ** 2 / (2 * np.pi) ** 4 * f ** -5
              * np.exp(-1.25 * (fp / f) ** 4) * gam ** r)
        dfdk = domegadk(k, h) / (2 * np.pi)
        s2 = np.trapz(k ** 2 * Sf * dfdk, k)
        m0 = np.trapz(Sf * dfdk, k)
        out[gam] = dict(Ft=Ft, fp=fp, alpha=alpha, sigma2=s2, m0=m0,
                        Hs=4*np.sqrt(m0) if m0>0 else 0,
                        lam_p=G/(2*np.pi*fp**2))
    return out

for U, F, h in ((1.5, 100.0, 0.6), (1.5, 8.9, 0.6), (5.0, 100.0, 1.0), (8.0, 100.0, 0.6)):
    print("U=%.1f F=%.1f h=%.2f  Ft=%.1f  CM=%.5f" % (U, F, h, G*F/U**2, 0.003+5.12e-3*U))
    r = jonswap_sigma2(U, F, h, gammas=(1.0, 1.7, 3.3, 5.0))
    for gam, v in r.items():
        print("   gam=%.1f sigma2=%.5f m0=%.3e Hs=%.4f fp=%.3f lam_p=%.3f alpha=%.4f"
              % (gam, v["sigma2"], v["m0"], v["Hs"], v["fp"], v["lam_p"], v["alpha"]))
    print("   ratio to CM:", {gam: round(v["sigma2"]/(0.003+5.12e-3*U), 3) for gam, v in r.items()})
