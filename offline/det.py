"""Deterministic Jacobian-equivalent forward transport (PIPELINE-V2 main path).

Method (OUTPUT-SPEC section 4, PIPELINE-V2 item 1):
- The water surface state (2048^2 grid, Choppy-displaced real vertices) is
  forward-mapped onto the bottom grid. Every surface texel emits exactly one
  stratified ray per solar-disc quadrature direction, Snell-refracted with a
  representative band IOR (G = 1.33434; the 1ch master C is achromatic, the
  residual chromatic spread is documented and left to the photon check).
- Finite solar disc: Q-point deterministic Fibonacci-disc quadrature over the
  disc (radius shared with transport.sun_disc_radius, i.e. identical
  condition to the photon path). Equal weights.
- Each ray is splatted bilinearly into the bottom buffer (periodic wrap, so
  energy is conserved) and the buffer is divided by the bottom-cell area.
  Bilinear splat + sigma=0.75 px FFT Gaussian = the same pixel filter as the
  photon path (transport.py), so the two fields are directly comparable.
- Because forward mapping collects every landing ray, all inverse images are
  summed structurally (no explicit branch search). Point-source singularities
  are NOT clamped.
- Output is the dimensionless illuminance ratio C = E_direct / E_base_dir,
  where E_base_dir is the flat-water bottom irradiance (definition and value
  recorded in numbers/metadata, not applied as a rescaling).

Energy bookkeeping: every ray carries the same analytical weight
w = dx_src^2 (source cell area). C_cell = acc_cell / (N_rays_total/N_buf^2).
With N_src == N_buf the divisor is exactly Q, a known constant -- the buffer
is NEVER rescaled to unit mean, clamped, or tone-mapped. The observed
mean(C) ~= 1 is a conservation consequence, reported as a diagnostic.

Companion product: flow_backward on the master bottom grid, the reverse
displacement used by image formation (G band), in master-grid pixels.
Convention is documented in metadata.json (origin / pixel-center / axes).
"""
import numpy as np

from transport import orthonormal_basis, sun_disc_radius, gaussian_blur_fft

IOR_G = 1.33434
SIGMA_PX = 0.75


def fibonacci_disc(Q, radius):
    """Deterministic equal-weight quadrature points on a disc.

    r_i = R*sqrt((i+0.5)/Q), theta_i = i * golden_angle. Returns (Q,2).
    """
    i = np.arange(Q, dtype=np.float64)
    r = float(radius) * np.sqrt((i + 0.5) / float(Q))
    th = i * (np.pi * (3.0 - np.sqrt(5.0)))
    return np.stack([r * np.cos(th), r * np.sin(th)], axis=1)


def snell_vec(I, N, eta):
    cos1 = -(I * N).sum(axis=-1, keepdims=True)
    cos1 = np.clip(cos1, -1.0, 1.0)
    k = 1.0 - eta * eta * (1.0 - cos1 * cos1)
    k = np.maximum(k, 0.0)
    return eta * I + (eta * cos1 - np.sqrt(k)) * N


def det_frame(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256):
    """Deterministic forward transport for one time slice.

    fields: dict with h/gx/gy/Dx/Dy (N_s,N_s) float32 + S not required.
    Returns (C, info): C (N_buf,N_buf) float32, dimensionless, unclamped.
    """
    h = fields["h"].astype(np.float64)
    gx = fields["gx"].astype(np.float64)
    gy = fields["gy"].astype(np.float64)
    Dx = fields["Dx"].astype(np.float64)
    Dy = fields["Dy"].astype(np.float64)
    N_s = h.shape[0]
    assert h.shape == (N_s, N_s)
    S = float(S)
    N_buf = int(N_buf)
    depth = float(env["depthM"])

    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    quads = fibonacci_disc(int(Q), disc)
    eta = 1.0 / float(ior)

    cell = S / N_buf
    acc = np.zeros(N_buf * N_buf, dtype=np.float64)

    xs_1d = (np.arange(N_s, dtype=np.float64) + 0.5) * (S / N_s)
    for q in range(int(Q)):
        jx, jy = float(quads[q, 0]), float(quads[q, 1])
        d = s + u * jx + v * jy
        d = d / np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])
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
            nb = N_buf * N_buf
            # each ray carries unit weight; physical area normalization below
            acc += np.bincount(i00, weights=w00, minlength=nb)
            acc += np.bincount(i10, weights=w10, minlength=nb)
            acc += np.bincount(i01, weights=w01, minlength=nb)
            acc += np.bincount(i11, weights=w11, minlength=nb)
            del X0, Y0, H, GX, GY, DX, DY, nl, Nx, Ny, Nz

    total_rays = float(int(Q)) * float(N_s) * float(N_s)
    # analytical divisor: rays per bottom cell for a uniform field.
    # NOT a measured rescaling (mean is conserved, not enforced).
    divisor = total_rays / float(N_buf * N_buf)
    C = (acc / divisor).reshape(N_buf, N_buf)
    C = gaussian_blur_fft(C, SIGMA_PX).astype(np.float32)
    info = {
        "Q": int(Q),
        "ior": float(ior),
        "disc_rad": float(disc),
        "sigma_px": float(SIGMA_PX),
        "total_rays": float(total_rays),
        "divisor_analytical": float(divisor),
        "clamp": "none",
        "renormalization": "none (mean~1 by energy conservation)",
    }
    return C, info


def flow_backward_master(fields, S, N_buf, env, ior=IOR_G):
    """Reverse displacement on the master bottom grid (G band).

    flow[j,i] = (dx, dy) in master-grid pixels such that the wobbled sample
    of a master-sized photo is at (i+0.5+dx, j+0.5+dy) (pixel-center
    convention). Origin: top-left corner of the patch (u=0,v=0); +x = +u
    (right), +y = +v (down, row index increasing). float32 (N_buf,N_buf,2).
    Same shift formula as image formation (image.render_frame_v2).
    """
    gx = fields["gx"].astype(np.float64)
    gy = fields["gy"].astype(np.float64)
    N_s = gx.shape[0]
    S = float(S)
    N_buf = int(N_buf)
    depth = float(env["depthM"])
    Cb = (float(ior) - 1.0) / float(ior)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    shx = -nx * depth * Cb  # meters, +u
    shy = -ny * depth * Cb  # meters, +v
    if N_s != N_buf:
        raise ValueError("flow grid requires N_s == N_buf (got %d vs %d)"
                         % (N_s, N_buf))
    px_per_m = N_buf / S
    flow = np.stack([shx * px_per_m, shy * px_per_m], axis=-1)
    return flow.astype(np.float32)
