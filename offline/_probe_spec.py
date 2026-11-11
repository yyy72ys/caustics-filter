import json
import numpy as np

from env import resolve_environment
from v6_state import V6WaterState
from v8_state import V8WaterState

with open("../blender/recipe-1.json", "r", encoding="utf-8") as f:
    raw = json.load(f)
env = resolve_environment(raw)

def check(tag, st):
    S = st._spectrum_density()
    k = st.k.astype(np.float64)
    s2 = float((k * k * S).sum())
    m0 = float(S.sum())
    h = np.real(np.fft.ifft2(st.H0.astype(np.complex128)))
    gx = np.real(np.fft.ifft2(1j * st.kx.astype(np.float64) * st.H0.astype(np.complex128)))
    gy = np.real(np.fft.ifft2(1j * st.ky.astype(np.float64) * st.H0.astype(np.complex128)))
    print("%-22s spec sigma2=%.6f m0=%.6e | direct field sigma2=%.6f m0=%.6e  ratio=%.3f"
          % (tag, s2, m0, float(np.mean(gx*gx+gy*gy)), float(np.mean(h*h)),
             float(np.mean(gx*gx+gy*gy))/s2))

st0 = V6WaterState(raw, env, N=256, preset_index=1)
check("v6 baseline", st0)
print("  N=%d kmin=%.4f kmax=%.2f" % (st0.N, st0.kmin, st0.kmax_nyq))
print("  k[0,:5]=", st0.k[0, :5], " kx[0,:5]=", st0.kx[0, :5])

st1 = V8WaterState(raw, env, N=256, preset_index=1, fetch_m=100.0)
check("v8 final", st1)
