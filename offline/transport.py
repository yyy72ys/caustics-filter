"""Light transport: forward photon mapping onto the bottom plane.

Per frame N_PHOTONS photons are sampled uniformly over the patch. Incident
direction = sun direction + solar-disc jitter (disc radius 0.00465 rad at
clarity 1, widened by (1-clarity): r = 0.00465*(0.7+0.6*(1-clarity))).
At the water surface each photon is Snell-refracted with band-dependent IOR
(R 1.33132 / G 1.33434 / B 1.33925); one photon carries 3 sub-rays.
Surface points use the choppy-displaced position and the bicubic-sampled
normal from the shared water state. Bottom intersection (z = -depth) is
splatted bilinearly into a 2048x2048 HDR buffer (periodic wrap conserves
energy), then smoothed with a sigma=0.75 px Gaussian (FFT).
Normalization is physical: photon flux = sunI*sunColor*patch_area/N, divided
by cell area, so the buffer mean equals the analytic expectation
(sunI*sunColor) up to float error. Entry surface reflectance is neglected
(energy-conserving into water; documented approximation).
"""
import numpy as np

IOR = np.array([1.33132, 1.33434, 1.33925], dtype=np.float64)

from state import sample_bicubic  # noqa: E402


def sun_disc_radius(clarity):
    return 0.00465 * (0.7 + 0.6 * (1.0 - float(clarity)))


def orthonormal_basis(s):
    s = np.asarray(s, dtype=np.float64)
    s = s / np.linalg.norm(s)
    up = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(s, up))) > 0.99:
        up = np.array([1.0, 0.0, 0.0])
    u = up - s * float(np.dot(up, s))
    u = u / np.linalg.norm(u)
    v = np.cross(s, u)
    return s, u, v


def snell(I, N, eta):
    """Refract incident I (into surface) with normals N, eta=n1/n2.

    I, N: (...,3). Returns T (...,3); TIR mask (all False here since eta<1).
    """
    cos1 = -(I * N).sum(axis=-1, keepdims=True)
    cos1 = np.clip(cos1, -1.0, 1.0)
    k = 1.0 - eta * eta * (1.0 - cos1 * cos1)
    k = np.maximum(k, 0.0)
    T = eta * I + (eta * cos1 - np.sqrt(k)) * N
    return T


def gaussian_blur_fft(img, sigma_px):
    """Separable-equivalent Gaussian blur via FFT (periodic). img: (H,W[,C])."""
    single = (img.ndim == 2)
    if single:
        img = img[:, :, None]
    H, W, C = img.shape
    fy = np.fft.fftfreq(H)[:, None]
    fx = np.fft.fftfreq(W)[None, :]
    Hf = np.exp(-2.0 * (np.pi ** 2) * (sigma_px ** 2) * (fx ** 2 + fy ** 2))
    out = np.empty_like(img)
    for c in range(C):
        F = np.fft.rfft2(img[:, :, c])
        Hfr = Hf[:, :F.shape[1]]
        out[:, :, c] = np.fft.irfft2(F * Hfr, s=(H, W)).astype(img.dtype)
    if single:
        return out[:, :, 0]
    return out


def photon_frame(fields, S, N_buf, env, n_photons, seed, chunk=250_000):
    """Trace one frame. Returns (E, info) with E (N_buf,N_buf,3) float32."""
    rng = np.random.default_rng(seed)
    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    depth = float(env["depthM"])
    sunI = float(env["sunI"])
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    patch_area = S * S
    cell = S / N_buf
    cell_area = cell * cell

    acc = np.zeros((N_buf * N_buf, 3), dtype=np.float64)
    n_done = 0
    while n_done < n_photons:
        m = min(chunk, n_photons - n_done)
        x0 = rng.random(m) * S
        y0 = rng.random(m) * S
        # solar disc jitter (uniform on disc)
        rr = disc * np.sqrt(rng.random(m))
        ph = rng.random(m) * 2.0 * np.pi
        jx = rr * np.cos(ph)
        jy = rr * np.sin(ph)
        d = (s[None, :] + u[None, :] * jx[:, None] + v[None, :] * jy[:, None])
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        I = -d  # incident, pointing down

        h = sample_bicubic(fields["h"], x0, y0, S).astype(np.float64)
        Dx = sample_bicubic(fields["Dx"], x0, y0, S).astype(np.float64)
        Dy = sample_bicubic(fields["Dy"], x0, y0, S).astype(np.float64)
        gx = sample_bicubic(fields["gx"], x0, y0, S).astype(np.float64)
        gy = sample_bicubic(fields["gy"], x0, y0, S).astype(np.float64)
        nl = np.sqrt(gx * gx + gy * gy + 1.0)
        N = np.stack([-gx / nl, -gy / nl, np.ones(m) / nl], axis=1)
        Psx = x0 + Dx
        Psy = y0 + Dy
        Psz = h
        for b in range(3):
            eta = 1.0 / float(IOR[b])
            T = snell(I, N, eta)
            Tz = T[:, 2]
            Tz = np.minimum(Tz, -1e-6)
            t = (Psz + depth) / (-Tz)
            xb = np.mod(Psx + T[:, 0] * t, S)
            yb = np.mod(Psy + T[:, 1] * t, S)
            fx = xb / cell
            fy = yb / cell
            ix0 = np.floor(fx).astype(np.int64) % N_buf
            iy0 = np.floor(fy).astype(np.int64) % N_buf
            tx = fx - np.floor(fx)
            ty = fy - np.floor(fy)
            ix1 = (ix0 + 1) % N_buf
            iy1 = (iy0 + 1) % N_buf
            w = sunI * float(sunColor[b]) * patch_area / n_photons
            w00 = w * (1 - tx) * (1 - ty)
            w10 = w * tx * (1 - ty)
            w01 = w * (1 - tx) * ty
            w11 = w * tx * ty
            i00 = iy0 * N_buf + ix0
            i10 = iy0 * N_buf + ix1
            i01 = iy1 * N_buf + ix0
            i11 = iy1 * N_buf + ix1
            nb = N_buf * N_buf
            acc[:, b] += np.bincount(i00, weights=w00, minlength=nb)
            acc[:, b] += np.bincount(i10, weights=w10, minlength=nb)
            acc[:, b] += np.bincount(i01, weights=w01, minlength=nb)
            acc[:, b] += np.bincount(i11, weights=w11, minlength=nb)
        n_done += m

    E = (acc / cell_area).reshape(N_buf, N_buf, 3)
    E = gaussian_blur_fft(E.astype(np.float64), 0.75).astype(np.float32)
    expected = sunI * sunColor
    # NOTE: the mean must be accumulated in float64. np.mean on the float32
    # buffer is systematically biased (observed -0.2..-0.4% on folded
    # caustics data at 2048^2) while float64 accumulation matches theory
    # to ~1e-10. Values themselves are unaffected (astype rounds unbiased).
    measured = E.astype(np.float64).reshape(-1, 3).mean(axis=0)
    info = {
        "expected_mean": [float(x) for x in expected],
        "measured_mean": [float(x) for x in measured],
    }
    return E, info
