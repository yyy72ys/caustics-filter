import numpy as np

G = 9.81
GAMMA = 0.0728
RHO = 1000.0
TF = GAMMA / RHO


def omega(k, h):
    return np.sqrt((G * k + TF * k ** 3) * np.tanh(k * h))


def domegadk(k, h):
    w = omega(k, h)
    th = np.tanh(k * h)
    sch = 1.0 / np.cosh(np.clip(k * h, 0, 30)) ** 2
    return ((G + 3 * TF * k ** 2) * th + (G * k + TF * k ** 3) * h * sch) / (2.0 * w)


def sigma2(U, F, h, gam=2.5, kmax=5361.65, kmin=4.26, n=200000, full=False):
    Ft = G * F / U ** 2
    fp = 3.5 * G / U * Ft ** (-1.0 / 3.0)
    alpha = 0.076 * Ft ** (-0.22)
    k = np.logspace(np.log10(0.05), np.log10(kmax), n)
    f = omega(k, h) / (2 * np.pi)
    sj = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-((f - fp) ** 2) / (2 * sj ** 2 * fp ** 2))
    Sf = (alpha * G ** 2 / (2 * np.pi) ** 4 * f ** -5
          * np.exp(-1.25 * (fp / f) ** 4) * gam ** r)
    dfdk = domegadk(k, h) / (2 * np.pi)
    m = k >= kmin
    s2 = np.trapz((k ** 2 * Sf * dfdk)[m], k[m])
    m0 = np.trapz((Sf * dfdk)[m], k[m])
    return dict(Ft=Ft, fp=fp, alpha=alpha, sigma2=s2, m0=m0,
                Hs=4 * np.sqrt(m0), lam_p=G / (2 * np.pi * fp ** 2))


U, h = 1.5, 0.6
for F in (0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100, 200, 500, 1000):
    r = sigma2(U, F, h)
    print("F=%7.1f Ft=%8.2f sigma2=%.5f (CM ratio %.2f) Hs=%.4f lam_p=%.3f"
          % (F, r["Ft"], r["sigma2"], r["sigma2"]/(0.003+5.12e-3*U), r["Hs"], r["lam_p"]))

print()
print("kmax sensitivity at U=1.5 F=100 h=0.6:")
for kmax in (300, 1000, 3000, 5361.65, 10000, 1e5, 1e6, 1.8e7):
    r = sigma2(U, 100.0, h, kmax=kmax)
    print("  kmax=%.3g sigma2=%.6f  m0=%.4e Hs=%.4f" % (kmax, r["sigma2"], r["m0"], r["Hs"]))

print()
print("kmin sensitivity:")
for kmin in (0.5, 2.0, 4.26, 8.0):
    r = sigma2(U, 100.0, h, kmin=kmin)
    print("  kmin=%.2f sigma2=%.6f" % (kmin, r["sigma2"]))

print()
print("MoC: equivalent fetch for CM target at U=1.5 h=0.6")
from math import log
target = 0.003 + 5.12e-3*U
lo, hi = 0.05, 1000.0
for _ in range(80):
    mid = np.sqrt(lo*hi)
    if sigma2(U, mid, h)["sigma2"] > target:
        hi = mid
    else:
        lo = mid
print("  F=%.3f m  Ft=%.2f" % (np.sqrt(lo*hi), G*np.sqrt(lo*hi)/U**2))
