"""Water state: Phillips-spectrum height field on a periodic patch via FFT.

Grid N x N (default 2048), patch side S = 4 * specL [m].
Spectrum: Phillips-type P(k) = A * exp(-1/(kL)^2) / k^4 * D(theta) * exp(-(k*l)^2)
  with directional distribution D = cos^2 type around wind direction,
  micro-wave damping length l = 0.005 m.
Choppy horizontal displacement with lambda = 0.85 (spectral:
  Dx(k) = -i kx/k H(k) lambda, same for y).
Time dependence: omega_raw = sqrt(g*k), g = 9.81; omega is quantized to integer
  multiples of omega0 = 2*pi/8 s so the state is exactly 8 s periodic.
Phase is fixed deterministically by the recipe-derived seed.
Height/gradient/displacement fields are sampled at arbitrary points with
periodic bicubic (Catmull-Rom) interpolation; both the photon side and the
image-formation side reference the same state object.
"""
import hashlib
import json
import numpy as np

G = 9.81
CHOPPY_LAMBDA = 0.85
LOOP_T = 8.0
OMEGA0 = 2.0 * np.pi / LOOP_T
DAMP_L = 0.005


def derive_seed(recipe, preset_index=1):
    blob = json.dumps(recipe, sort_keys=True, ensure_ascii=False)
    blob += "|preset%d" % int(preset_index)
    h = hashlib.sha256(blob.encode("utf-8")).digest()
    return int.from_bytes(h[:8], "little") % (2 ** 63)


class WaterState:
    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1):
        self.N = int(N)
        self.specL = float(recipe_env["specL"])
        self.S = 4.0 * self.specL
        self.dx = self.S / self.N
        self.seed = int(seed) if seed is not None else derive_seed(recipe, preset_index)
        rng = np.random.default_rng(self.seed)

        fx = np.fft.fftfreq(self.N, d=self.dx)  # cycles/m
        kx1 = 2.0 * np.pi * fx
        kx, ky = np.meshgrid(kx1.astype(np.float64), kx1.astype(np.float64))
        self.kx = kx.astype(np.float32)
        self.ky = ky.astype(np.float32)
        k = np.sqrt(kx * kx + ky * ky)
        self.k = k.astype(np.float32)

        L = self.specL
        wdx, wdy = recipe_env["windDir"]
        theta_w = np.arctan2(wdy, wdx)
        theta = np.arctan2(ky, kx)
        c = np.cos(theta - theta_w)
        D = np.where(c > 0.0, c * c, 0.0) + 0.02  # cos^2 lobe + isotropic floor

        kk = k.copy()
        kk[kk == 0.0] = 1.0
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            P = (np.exp(-1.0 / (kk * L) ** 2) / (kk ** 4)
                 * D * np.exp(-((kk * DAMP_L) ** 2)))
        P[k == 0.0] = 0.0
        P[~np.isfinite(P)] = 0.0
        self.P = P.astype(np.float64)

        a = rng.standard_normal((self.N, self.N))
        b = rng.standard_normal((self.N, self.N))
        F0_raw = ((a + 1j * b) / np.sqrt(2.0)) * np.sqrt(self.P)
        h_raw = np.real(np.fft.ifft2(F0_raw))
        std_raw = float(np.std(h_raw))
        target_std = float(recipe_env["waveHeightM"])
        if not np.isfinite(std_raw) or std_raw <= 0.0:
            raise RuntimeError("degenerate raw height field")
        self.amplitude_scale = target_std / std_raw
        self.H0 = (F0_raw * self.amplitude_scale).astype(np.complex64)
        self.target_std = target_std

        # band edges / diagnostics
        self.kmin = 2.0 * np.pi / self.S
        self.kmax_nyq = np.pi / self.dx
        self.decades = float(np.log10(self.kmax_nyq / self.kmin))

        # keep float64 k for omega computation, release mesh temps
        del kx, ky, k, kk, P, a, b, F0_raw, h_raw

    def _omega_q(self):
        raw = np.sqrt(G * np.maximum(self.k.astype(np.float64), 0.0))
        n = np.maximum(1.0, np.floor(raw / OMEGA0 + 0.5))
        return (n * OMEGA0).astype(np.float64)

    def fields_at(self, t):
        """Return dict of spatial fields at time t (seconds).

        Keys: h, gx, gy, Dx, Dy — each (N,N) float32 in meters (gradients
        dimensionless). Periodic in t with period 8 s by construction.
        """
        wq = self._omega_q()
        phase = np.exp((1j * wq * float(t))).astype(np.complex64)
        F = self.H0 * phase
        del phase
        h = np.real(np.fft.ifft2(F)).astype(np.float32)
        kx = self.kx.astype(np.complex64)
        ky = self.ky.astype(np.complex64)
        gx = np.real(np.fft.ifft2(1j * kx * F)).astype(np.float32)
        gy = np.real(np.fft.ifft2(1j * ky * F)).astype(np.float32)
        kinv = np.zeros_like(self.k)
        nz = self.k > 0
        kinv[nz] = 1.0 / self.k[nz]
        kinv = kinv.astype(np.complex64)
        lam = np.complex64(CHOPPY_LAMBDA)
        Dx = np.real(np.fft.ifft2(-1j * kx * kinv * lam * F)).astype(np.float32)
        Dy = np.real(np.fft.ifft2(-1j * ky * kinv * lam * F)).astype(np.float32)
        return {"h": h, "gx": gx, "gy": gy, "Dx": Dx, "Dy": Dy}

    def loop_maxdiff(self):
        f0 = self.fields_at(0.0)["h"].astype(np.float64)
        f8 = self.fields_at(LOOP_T)["h"].astype(np.float64)
        return float(np.max(np.abs(f8 - f0)))


def sample_bicubic(grid, xs_m, ys_m, S):
    """Periodic Catmull-Rom bicubic sampling.

    grid: (N,N) float32, rows = y. xs_m/ys_m: arrays in meters (any range,
    wrapped). Returns float32 array of same shape as xs_m.
    """
    N = grid.shape[0]
    dx = S / N
    gx = (np.asarray(xs_m, dtype=np.float64) / dx)
    gy = (np.asarray(ys_m, dtype=np.float64) / dx)
    ix0 = np.floor(gx).astype(np.int64)
    iy0 = np.floor(gy).astype(np.int64)
    tx = (gx - np.floor(gx)).astype(np.float64)
    ty = (gy - np.floor(gy)).astype(np.float64)

    def w_cubic(t):
        t = np.abs(t)
        a = -0.5
        w1 = (a + 2.0) * t ** 3 - (a + 3.0) * t ** 2 + 1.0
        w2 = a * t ** 3 - 5.0 * a * t ** 2 + 8.0 * a * t - 4.0 * a
        return np.where(t <= 1.0, w1, np.where(t < 2.0, w2, 0.0))

    wx = [w_cubic(tx + 1), w_cubic(tx), w_cubic(1.0 - tx), w_cubic(2.0 - tx)]
    wy = [w_cubic(ty + 1), w_cubic(ty), w_cubic(1.0 - ty), w_cubic(2.0 - ty)]
    g = grid.astype(np.float64)
    out = np.zeros_like(tx)
    for j in range(4):
        iym = (iy0 + (j - 1)) % N
        row = np.zeros_like(tx)
        for i in range(4):
            ixm = (ix0 + (i - 1)) % N
            row += wx[i] * g[iym, ixm]
        out += wy[j] * row
    return out.astype(np.float32)


def sample_bilinear(grid, xs_m, ys_m, S, wrap=True):
    N = grid.shape[0]
    dx = S / N
    gx = np.asarray(xs_m, dtype=np.float64) / dx - 0.5
    gy = np.asarray(ys_m, dtype=np.float64) / dx - 0.5
    ix0 = np.floor(gx).astype(np.int64)
    iy0 = np.floor(gy).astype(np.int64)
    tx = (gx - ix0).astype(np.float64)
    ty = (gy - iy0).astype(np.float64)
    if wrap:
        ix0 %= N
        iy0 %= N
        ix1 = (ix0 + 1) % N
        iy1 = (iy0 + 1) % N
    else:
        ix0 = np.clip(ix0, 0, N - 1)
        iy0 = np.clip(iy0, 0, N - 1)
        ix1 = np.clip(ix0 + 1, 0, N - 1)
        iy1 = np.clip(iy0 + 1, 0, N - 1)
    g = grid.astype(np.float64)
    v = (g[iy0, ix0] * (1 - tx) * (1 - ty) + g[iy0, ix1] * tx * (1 - ty)
         + g[iy1, ix0] * (1 - tx) * ty + g[iy1, ix1] * tx * ty)
    return v.astype(np.float32)
