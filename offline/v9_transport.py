"""v9 weighted-direction transport -- chunk 5 (version 3.5).

``det_frame`` (det.py, FROZEN) traces Q equal-weight rays per surface texel,
all near the solar direction, and returns one C master.  The diffuse component
needs the same Jacobian transport run along many sky directions with unequal
weights.  This module is that generalisation, and it is a *faithful copy* of
det_frame's geometry rather than a re-derivation:

  * same source grid, same choppy displacement (Dx, Dy), same height h,
  * same normal from (gx, gy), same Snell refraction with the same IOR,
  * same bottom intersection z = -depth, same periodic bilinear splat,
  * same analytical divisor (rays per bottom cell for a uniform field) so the
    mean is a conservation consequence, never a rescaling,
  * same sigma = 0.75 px FFT pixel filter, scaled by N_buf/2048 so the filter
    has the SAME PHYSICAL WIDTH at every master resolution (a 0.75 px filter
    at 512^2 would be a 3 px filter at 2048^2 and would silently smooth the
    diffuse field it is supposed to measure).

The fidelity is checked at run time: ``disc_dirs`` reproduces the solar-disc
quadrature of det.py exactly, and ``det_frame_dirs`` with equal weights on
those directions must agree with ``det_frame`` to float32 rounding.  The
measured max difference is written into numbers-v9.json.

No frozen file is imported-and-modified: det.py is only imported to READ
IOR_G/SIGMA_PX constants.
"""
import math

import numpy as np

from det import IOR_G, SIGMA_PX
from transport import gaussian_blur_fft, orthonormal_basis, sun_disc_radius


def disc_dirs(env, Q):
    """The solar-disc quadrature of det.det_frame, reproduced bit-for-bit.

    Returns (dirs (Q,3), weights (Q,)) with weights = 1/Q, i.e. exactly the
    condition det_frame uses.
    """
    from det import fibonacci_disc
    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    quads = fibonacci_disc(int(Q), disc)
    d = s[None, :] + u[None, :] * quads[:, 0:1] + v[None, :] * quads[:, 1:2]
    d = d / np.linalg.norm(d, axis=1, keepdims=True)
    return d, np.full(int(Q), 1.0 / float(Q), dtype=np.float64)


def det_frame_dirs(fields, S, N_buf, env, dirs, weights=None, ior=IOR_G,
                   row_block=256, sigma_px=None, progress=None, tag=""):
    """Weighted sum of the deterministic forward transport over directions.

    dirs    (K,3) unit vectors pointing TOWARD the sky (same convention as
            env["sunDir"]).
    weights (K,)  normalised (sum = 1) or None (uniform).
    Returns (C (N_buf,N_buf) float32, info dict).
    """
    h = fields["h"].astype(np.float64)
    gx = fields["gx"].astype(np.float64)
    gy = fields["gy"].astype(np.float64)
    Dx = fields["Dx"].astype(np.float64)
    Dy = fields["Dy"].astype(np.float64)
    N_s = h.shape[0]
    S = float(S)
    N_buf = int(N_buf)
    depth = float(env["depthM"])

    d_all = np.asarray(dirs, dtype=np.float64)
    K = int(d_all.shape[0])
    if weights is None:
        w_all = np.full(K, 1.0 / float(K), dtype=np.float64)
    else:
        w_all = np.asarray(weights, dtype=np.float64)
        w_all = w_all / w_all.sum()
    eta = 1.0 / float(ior)

    cell = S / N_buf
    nb = N_buf * N_buf
    acc = np.zeros(nb, dtype=np.float64)
    acc_dir = np.zeros(nb, dtype=np.float64)
    xs_1d = (np.arange(N_s, dtype=np.float64) + 0.5) * (S / N_s)
    per_dir_mean = []

    for q in range(K):
        wq = float(w_all[q])
        if wq == 0.0:
            per_dir_mean.append(None)
            continue
        I = -d_all[q]
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])
        acc_dir.fill(0.0)
        for r0 in range(0, N_s, row_block):
            r1 = min(r0 + row_block, N_s)
            yy = xs_1d[r0:r1]
            xx = xs_1d
            X0, Y0 = np.meshgrid(xx, yy)
            H = h[r0:r1, :]
            GX = gx[r0:r1, :]
            GY = gy[r0:r1, :]
            DX = Dx[r0:r1, :]
            DY = Dy[r0:r1, :]
            nl = np.sqrt(GX * GX + GY * GY + 1.0)
            Nx = -GX / nl
            Ny = -GY / nl
            Nz = 1.0 / nl
            cos1 = -(Ix * Nx + Iy * Ny + Iz * Nz)
            cos1 = np.clip(cos1, -1.0, 1.0)
            kk = 1.0 - eta * eta * (1.0 - cos1 * cos1)
            kk = np.maximum(kk, 0.0)
            sq = np.sqrt(kk)
            Tx = eta * Ix + (eta * cos1 - sq) * Nx
            Ty = eta * Iy + (eta * cos1 - sq) * Ny
            Tz = eta * Iz + (eta * cos1 - sq) * Nz
            Tz = np.minimum(Tz, -1e-6)
            Psx = X0 + DX
            Psy = Y0 + DY
            Psz = H
            t = (Psz + depth) / (-Tz)
            xb = np.mod(Psx + Tx * t, S)
            yb = np.mod(Psy + Ty * t, S)
            fx = (xb / cell).ravel()
            fy = (yb / cell).ravel()
            ix0 = np.floor(fx).astype(np.int64) % N_buf
            iy0 = np.floor(fy).astype(np.int64) % N_buf
            tx = fx - np.floor(fx)
            ty = fy - np.floor(fy)
            ix1 = (ix0 + 1) % N_buf
            iy1 = (iy0 + 1) % N_buf
            w00 = (1 - tx) * (1 - ty)
            w10 = tx * (1 - ty)
            w01 = (1 - tx) * ty
            w11 = tx * ty
            i00 = iy0 * N_buf + ix0
            i10 = iy0 * N_buf + ix1
            i01 = iy1 * N_buf + ix0
            i11 = iy1 * N_buf + ix1
            acc_dir += np.bincount(i00, weights=w00, minlength=nb)
            acc_dir += np.bincount(i10, weights=w10, minlength=nb)
            acc_dir += np.bincount(i01, weights=w01, minlength=nb)
            acc_dir += np.bincount(i11, weights=w11, minlength=nb)
            del X0, Y0, H, GX, GY, DX, DY, nl, Nx, Ny, Nz
        # one weighted add per direction (keeps the 2048^2 temporaries down)
        acc += wq * acc_dir
        per_dir_mean.append(None if progress is None else 0)
        if progress is not None and (q % 16 == 0 or q == K - 1):
            progress(tag, q + 1, K)

    # analytical divisor: rays per bottom cell for a uniform field
    divisor = float(N_s) * float(N_s) / float(nb)
    C = (acc / divisor).reshape(N_buf, N_buf)
    sig = (SIGMA_PX * float(N_buf) / 2048.0 if sigma_px is None
           else float(sigma_px))
    if sig > 0.0:
        C = gaussian_blur_fft(C, sig)
    info = {
        "n_dirs": K,
        "N_src": int(N_s), "N_buf": int(N_buf),
        "divisor_analytical": float(divisor),
        "sigma_px": float(sig),
        "sigma_px_equivalent_at_2048": float(sig * 2048.0 / float(N_buf)),
        "weights_sum": float(w_all.sum()),
        "clamp": "none",
        "renormalization": "none (mean~1 by energy conservation)",
    }
    return C.astype(np.float32), info


def mean_of(arr):
    return float(np.asarray(arr, dtype=np.float64).mean())


def fft_upsample(C, N):
    """Band-limited (Fourier zero-pad) upsample of a periodic master.

    Used to extend the 512^2 diffuse master to the 2048^2 transport grid.
    Unlike cubic interpolation this PRESERVES THE MEAN EXACTLY (the k=0 mode
    is untouched) and adds no energy above the source band -- which is the
    honest way to say "the field has no structure finer than this": nothing
    is invented, and the DC (energy) term cannot drift.
    """
    a = np.asarray(C, dtype=np.float64)
    n = a.shape[0]
    N = int(N)
    if n == N:
        return a.astype(np.float32)
    if N % n != 0 and N < n:
        raise ValueError("upsample target must be larger (got %d -> %d)"
                         % (n, N))
    F = np.fft.fftshift(np.fft.fft2(a))
    big = np.zeros((N, N), dtype=np.complex128)
    h = n // 2
    c = N // 2
    big[c - h:c - h + n, c - h:c - h + n] = F
    out = np.real(np.fft.ifft2(np.fft.ifftshift(big))) * (float(N) / n) ** 2
    return out.astype(np.float32)
