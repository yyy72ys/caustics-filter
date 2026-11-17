"""det_fast -- faster det_frame with the same deposit (NEW file, BRIEF-DET-SPEEDUP).

Equivalence contract
--------------------
`det_frame_fast` reproduces `det.det_frame` **bit for bit** (float32 C,
verified on production frames 0 and 24 against `out/prod/C-*.npy` and
against `det.det_frame` itself).  Nothing in the physics or the pixel
filter is changed:

  * solar-disc quadrature: `det.fibonacci_disc(Q, sun_disc_radius(clarity))`
    (imported);
  * Snell refraction with eta = 1/ior, per-texel normals from gx/gy, the
    finite-depth hit t = (h + depth)/(-Tz), `np.mod`-equivalent periodic
    wrap (fmod + the negative branch shifted by +S, exactly numpy's
    remainder), the analytic divisor rays/(N_buf^2), the periodic bilinear
    splat, and the sigma = 0.75 px FFT Gaussian
    (`transport.gaussian_blur_fft`, imported);
  * the deposit keeps production's **addition order**: for every row
    block of 256 rows and every direction q (q outermost), the four
    per-bin sums -- i00, i10, i01, i11 -- are summed in ascending ray
    order and then added into the accumulator once per bin, exactly as
    `acc += np.bincount(idx, weights=w, minlength=nb)` four times would.
    That is what makes the output bit-identical, so it is preserved.

Where the speed comes from
--------------------------
1. Fusion without temporaries.  `det.det_frame` materialises ~30 four-MiB
   temporaries per (q, row block) and calls `np.bincount` four times per
   block, each zero-filling a 4 Mi-element accumulator that is then added
   element-wise into `acc`: ~6.8 G element touches per frame.
2. The per-ray geometry (normals, Snell, bottom hit, bilinear indices) is
   embarrassingly parallel *per ray* -- every ray's value depends on its
   own texel only.  It is computed by a `prange` kernel over row chunks,
   so the CPU's cores are used without touching the summation order.
   Bit-exactness is unaffected: no value depends on scheduling.
3. The deposit itself stays single-threaded and in production order.  The
   four splats of a block use one scratch buffer with a generation stamp,
   so only the bins that were actually hit are committed and cleared.

The measured split is in out/det-speedup/REPORT-DET-SPEEDUP.md.

Public API
----------
det_frame_fast(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
               chunk_rows=64, work=None, n_threads=None)
    -> (C float32 (N_buf,N_buf), info dict)      # same as det.det_frame
warmup(n_threads=None)   -- compile once; timing must exclude this
"""
import os

import numpy as np

from det import IOR_G, SIGMA_PX, fibonacci_disc
from transport import orthonormal_basis, sun_disc_radius, gaussian_blur_fft

try:                                    # numba is present on this machine
    import numba
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:                       # pragma: no cover - fallback path
    HAVE_NUMBA = False
    numba = None

    def njit(*a, **kw):
        def deco(f):
            return f
        if a and callable(a[0]):
            return a[0]
        return deco

    def prange(*a):
        return range(*a)


# --------------------------------------------------------------------------
# per-ray geometry (parallel over row chunks; every ray independent)
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _geom_chunk(h, gx, gy, Dx, Dy, r0, r1, N_s, off, S, s_over_ns, depth,
                cell, N_buf, Ix, Iy, Iz, eta,
                ix0_o, iy0_o, w00_o, w10_o, w01_o, w11_o):
    """One row chunk of det.det_frame lines 102-136.

    Operation order is that of the numpy original, expression by
    expression, so the values are bit-identical; only the storage of the
    results changes (int32 indices instead of int64 temporaries).
    """
    for r in range(r0, r1):
        y0 = (r + 0.5) * s_over_ns
        base = off + (r - r0) * N_s
        for c in range(N_s):
            gxv = gx[r, c]
            gyv = gy[r, c]
            nl = np.sqrt(gxv * gxv + gyv * gyv + 1.0)
            Nx = -gxv / nl
            Ny = -gyv / nl
            Nz = 1.0 / nl
            cos1 = -(Ix * Nx + Iy * Ny + Iz * Nz)
            if cos1 < -1.0:
                cos1 = -1.0
            elif cos1 > 1.0:
                cos1 = 1.0
            kk = 1.0 - eta * eta * (1.0 - cos1 * cos1)
            if kk < 0.0:
                kk = 0.0
            sq = np.sqrt(kk)
            Txa = eta * cos1 - sq
            Tx = eta * Ix + Txa * Nx
            Ty = eta * Iy + Txa * Ny
            Tz = eta * Iz + Txa * Nz
            if Tz > -1e-6:
                Tz = -1e-6
            Psx = (c + 0.5) * s_over_ns + Dx[r, c]
            Psy = y0 + Dy[r, c]
            Psz = h[r, c]
            t = (Psz + depth) / (-Tz)
            # numpy's np.mod(P, S): fmod, and the negative branch + S.
            # ~14 % of the x hits are negative, so this matters.
            xb = np.fmod(Psx + Tx * t, S)
            if xb < 0.0:
                xb += S
            yb = np.fmod(Psy + Ty * t, S)
            if yb < 0.0:
                yb += S
            fx = xb / cell
            fy = yb / cell
            fxi = np.floor(fx)
            fyi = np.floor(fy)
            i0 = np.int32(fxi)          # fx in [0, N_buf) by construction
            j0 = np.int32(fyi)
            tx = fx - fxi
            ty = fy - fyi
            k = base + c
            ix0_o[k] = i0
            iy0_o[k] = j0
            w00_o[k] = (1.0 - tx) * (1.0 - ty)
            w10_o[k] = tx * (1.0 - ty)
            w01_o[k] = (1.0 - tx) * ty
            w11_o[k] = tx * ty


@njit(cache=False, fastmath=False, parallel=True, nogil=True)
def _geom_all(h, gx, gy, Dx, Dy, N_s, S, s_over_ns, depth, cell, N_buf,
              Ix, Iy, Iz, eta, chunk_rows,
              ix0_o, iy0_o, w00_o, w10_o, w01_o, w11_o):
    nchunks = (N_s + chunk_rows - 1) // chunk_rows
    for ci in prange(nchunks):
        r0 = ci * chunk_rows
        r1 = r0 + chunk_rows
        if r1 > N_s:
            r1 = N_s
        _geom_chunk(h, gx, gy, Dx, Dy, r0, r1, N_s, ci * chunk_rows * N_s,
                    S, s_over_ns, depth, cell, N_buf, Ix, Iy, Iz, eta,
                    ix0_o, iy0_o, w00_o, w10_o, w01_o, w11_o)


# --------------------------------------------------------------------------
# deposit: production order, one scratch buffer, generation stamps
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _splat4(acc, buf, gen, g, ix0_b, iy0_b, w00_b, w10_b, w01_b, w11_b,
            n, N_buf):
    """Four bilinear splats, each summed per bin in ascending ray order and
    added into `acc` once per bin -- the same arithmetic as

        acc += np.bincount(i00, w00, minlength=nb)   # then i10, i01, i11

    for the same block, in the same order.  `buf` is zero on entry and is
    restored to zero on the bins it touched.  `gen` is a per-bin
    generation stamp, so the commit touches only the bins that were hit.
    """
    # ---- pass 1: i00 ----
    for k in range(n):
        i = iy0_b[k] * N_buf + ix0_b[k]
        buf[i] += w00_b[k]
    for k in range(n):
        i = iy0_b[k] * N_buf + ix0_b[k]
        if gen[i] != g:
            gen[i] = g
            acc[i] += buf[i]
            buf[i] = 0.0
    g += 1
    # ---- pass 2: i10 ----
    for k in range(n):
        ix = ix0_b[k] + 1
        if ix == N_buf:
            ix = 0
        i = iy0_b[k] * N_buf + ix
        buf[i] += w10_b[k]
    for k in range(n):
        ix = ix0_b[k] + 1
        if ix == N_buf:
            ix = 0
        i = iy0_b[k] * N_buf + ix
        if gen[i] != g:
            gen[i] = g
            acc[i] += buf[i]
            buf[i] = 0.0
    g += 1
    # ---- pass 3: i01 ----
    for k in range(n):
        iy = iy0_b[k] + 1
        if iy == N_buf:
            iy = 0
        i = iy * N_buf + ix0_b[k]
        buf[i] += w01_b[k]
    for k in range(n):
        iy = iy0_b[k] + 1
        if iy == N_buf:
            iy = 0
        i = iy * N_buf + ix0_b[k]
        if gen[i] != g:
            gen[i] = g
            acc[i] += buf[i]
            buf[i] = 0.0
    g += 1
    # ---- pass 4: i11 ----
    for k in range(n):
        iy = iy0_b[k] + 1
        if iy == N_buf:
            iy = 0
        ix = ix0_b[k] + 1
        if ix == N_buf:
            ix = 0
        i = iy * N_buf + ix
        buf[i] += w11_b[k]
    for k in range(n):
        iy = iy0_b[k] + 1
        if iy == N_buf:
            iy = 0
        ix = ix0_b[k] + 1
        if ix == N_buf:
            ix = 0
        i = iy * N_buf + ix
        if gen[i] != g:
            gen[i] = g
            acc[i] += buf[i]
            buf[i] = 0.0
    return g + 1


def _set_threads(n_threads):
    if not HAVE_NUMBA:
        return
    if n_threads is None:
        return
    num_threads = int(n_threads)
    if num_threads > 0:
        numba.set_num_threads(min(num_threads, numba.config.NUMBA_NUM_THREADS))


def _njit_check():
    if not HAVE_NUMBA:
        raise RuntimeError("det_fast needs numba for the fused kernels")


def det_frame_fast(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
                   chunk_rows=64, work=None, n_threads=None):
    """Same result as det.det_frame, computed with fused/parallel kernels.

    row_block must equal the production value (256): it defines the
    summation grouping of the deposit.
    chunk_rows only splits the (order-independent) geometry work.
    work: optional dict reused across frames (buffers stay allocated).
    """
    _njit_check()
    _set_threads(n_threads)

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
    s_over_ns = S / N_s
    nb = N_buf * N_buf
    n_rays = N_s * N_s

    if work is None:
        work = {}
    if work.get("acc") is None or work["acc"].size != nb:
        work["acc"] = np.zeros(nb, dtype=np.float64)
        work["buf"] = np.zeros(nb, dtype=np.float64)
        work["gen"] = np.zeros(nb, dtype=np.int32)
        work["g"] = 1                       # gen starts at 0, so do stamps
    if work.get("ix0") is None or work["ix0"].size != n_rays:
        work["ix0"] = np.empty(n_rays, dtype=np.int32)
        work["iy0"] = np.empty(n_rays, dtype=np.int32)
        work["w00"] = np.empty(n_rays, dtype=np.float64)
        work["w10"] = np.empty(n_rays, dtype=np.float64)
        work["w01"] = np.empty(n_rays, dtype=np.float64)
        work["w11"] = np.empty(n_rays, dtype=np.float64)
    acc = work["acc"]
    buf = work["buf"]
    gen = work["gen"]
    ix0 = work["ix0"]
    iy0 = work["iy0"]
    w00 = work["w00"]
    w10 = work["w10"]
    w01 = work["w01"]
    w11 = work["w11"]

    acc.fill(0.0)                           # per-frame reset, as np.zeros
    g = int(work["g"])
    nb_i = np.int32(N_buf)

    for q in range(int(Q)):
        jx, jy = float(quads[q, 0]), float(quads[q, 1])
        d = s + u * jx + v * jy
        d = d / np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])

        _geom_all(h, gx, gy, Dx, Dy, N_s, S, s_over_ns, depth, cell, N_buf,
                  Ix, Iy, Iz, eta, int(chunk_rows),
                  ix0, iy0, w00, w10, w01, w11)

        for r0 in range(0, N_s, row_block):
            r1 = min(r0 + row_block, N_s)
            a = r0 * N_s
            b = r1 * N_s
            g = _splat4(acc, buf, gen, g, ix0[a:b], iy0[a:b], w00[a:b],
                        w10[a:b], w01[a:b], w11[a:b], b - a, nb_i)

    work["g"] = g
    total_rays = float(int(Q)) * float(N_s) * float(N_s)
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


def warmup(n_threads=None):
    """Compile both kernels once (tiny data) so timing excludes JIT cost."""
    _njit_check()
    _set_threads(n_threads)
    n = 32
    z = np.zeros((n, n), dtype=np.float64)
    nn = n * n
    acc = np.zeros(nn, dtype=np.float64)
    buf = np.zeros(nn, dtype=np.float64)
    gen = np.zeros(nn, dtype=np.int32)
    ix0 = np.zeros(nn, dtype=np.int32)
    iy0 = np.zeros(nn, dtype=np.int32)
    w = np.full(nn, 0.25, dtype=np.float64)
    _geom_all(z, z, z, z, z, n, 1.0, 1.0 / n, 0.6, 1.0 / n, n,
              0.0, 0.0, -1.0, 1.0, 8,
              ix0, iy0, w, w, w, w)
    _splat4(acc, buf, gen, 1, ix0, iy0, w, w, w, w, nn, np.int32(n))


if __name__ == "__main__":
    import argparse
    import json
    import sys
    import time

    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from env import resolve_environment
    from state import WaterState
    from det import det_frame

    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", type=int, default=0)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--recipe", default=os.path.join(here, "..", "blender",
                                                    "recipe-1.json"))
    args = ap.parse_args()

    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    env = resolve_environment(recipe)
    state = WaterState(recipe, env, N=2048, seed=None, preset_index=1)
    fields = state.fields_at(args.frame / 24.0)
    fields["S"] = state.S

    t0 = time.perf_counter()
    warmup(args.threads)
    t_jit = time.perf_counter() - t0
    t0 = time.perf_counter()
    C, info = det_frame_fast(fields, state.S, 2048, env, 25,
                             n_threads=args.threads)
    t_fast = time.perf_counter() - t0
    t0 = time.perf_counter()
    C_ref, _ = det_frame(fields, state.S, 2048, env, 25)
    t_slow = time.perf_counter() - t0
    d = np.abs(C.astype(np.float64) - C_ref.astype(np.float64))
    print("jit_warmup_s=%.2f fast_s=%.3f base_s=%.3f speedup=%.2fx"
          % (t_jit, t_fast, t_slow, t_slow / t_fast))
    print("bit_identical=%s max_abs_diff=%.3g"
          % (np.array_equal(C, C_ref), float(d.max())))
