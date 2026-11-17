"""det_fast2 -- det_frame with a stable counting-sort binned deposit (NEW file).

BRIEF-DET-BINSORT, STAGE 1.  Candidate module only: nothing here is wired
into production (run_prod_v4.py / DET_IMPL are untouched); the gate and the
bench scripts call `det_frame_fast2` explicitly.

What is new
-----------
`det_fast.det_frame_fast` (the current production det) deposits with 4 passes
x (scatter + commit) straight into a 32 MiB bin buffer in **ray order**, so
every one of the 200 (q, row-block) calls scatters through 32 MiB.

`det_fast2` inserts one step before the deposit: for each (q, row block) the
rays are **stably sorted by their destination class (iy0, ix0)** with a
counting sort, and the deposit then walks the rays in bin order.  Consecutive
rays write neighbouring bins instead of jumping around the buffer, so a bin
that is hit by several rays is finished in one visit and the 32 MiB scratch
buffer + generation stamps are not needed at all.

Why the result is still bit-identical
-------------------------------------
Production's deposit tree, per (q, row block) and per bin i, is

    acc[i] = ((((acc[i] + s1) + s2) + s3) + s4)

where s_p is the sum of the corner-p weights of every ray that lands in bin i
in pass p, accumulated in **ascending ray index** (`np.bincount` order), and
acc[i] is one single add per (pass, block).

The sort has to preserve exactly two things:

1. **Per-bin summation order.**  A stable sort keeps rays with equal keys in
   ascending index order, so every bin's sum still runs in ascending ray
   order -> the same floating-point tree inside s_p.

2. **One acc add per (bin, pass, block), in pass order 1,2,3,4.**  The sort
   key is the class (iy0, ix0) the four bilinear corners are derived from.
   For a fixed corner (dx,dy) the map

       (iy0, ix0) -> (iy0 + dy, ix0 + dx) mod N_buf        (dx,dy) in {0,1}^2

   is injective, so bin i receives the rays of exactly one class, and each
   class is contiguous after the sort.  A class can therefore be closed at
   its run boundary: one acc[i] add per (class, corner) == one acc[i] add per
   (bin, pass) == production.

   Deposit "pass4" walks the four corners in four separate passes (1,2,3,4),
   which is production's order with no further argument.
   Deposit "fused" walks the sorted block **once, in descending class
   order**, and closes all four corners at the same run boundary.  Descending
   class order means bin (Y,X) is reached as corner 0 of class (Y,X), then
   corner 1 of class (Y,X-1), then corner 2 of class (Y-1,X), then corner 3
   of class (Y-1,X-1) -- exactly production's 1,2,3,4 order, for a single
   pass over the weights.  The two cyclic seams (destination row 0 and
   destination column 0, where "Y-1"/"X-1" wraps) cannot be ordered that way,
   so those 2*N-1 bins are staged and flushed in production order at the end
   of the block.  See `_deposit_runs4`.

Nothing in the physics, the geometry kernel, the weights, the fmod wrap, the
divisor or the FFT is changed: the geometry stage is `det_fast._geom_all`,
imported unchanged.

Sort variants (`sort_mode`) -- all produce the same bin-ordered block
--------------------------------------------------------------------
  "radix" (default) two-level counting sort (11+11 bits: classes by row, then
          each row group by column).  Same output as "bin", but the counting
          arrays stay cache-resident instead of being a 16 MiB memset+prefix
          per block.
  "bin"   single-level counting sort on the full 22-bit class key -- the
          literal mechanism of the brief, kept for the measurement.
  "tile"  coarse stable sort on a square tile key + production's
          scatter/commit arithmetic confined to a tile window (the sorting
          cost floor; it cannot close a bin in one visit).

Public API (same signature as det_fast.det_frame_fast)
-----------------------------------------------------
det_frame_fast2(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
                chunk_rows=64, work=None, n_threads=None,
                sort_mode="radix", deposit="fused", tile=64, timing=None)
    -> (C float32 (N_buf,N_buf), info dict)      # same as det.det_frame
warmup(n_threads=None)   -- compile once; timing must exclude this
"""
import os
import time

import numpy as np

from det import IOR_G, SIGMA_PX, fibonacci_disc
from transport import orthonormal_basis, sun_disc_radius, gaussian_blur_fft

import det_fast as _df

try:                                    # numba is present on this machine
    import numba
    from numba import njit
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


# --------------------------------------------------------------------------
# sort 1: single-level counting sort on the class key (brief mechanism)
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _offsets(ix0, iy0, cnt, n, nbins, shift, box_shift, ntx, mode):
    """Count the block's rays per bucket, then exclusive prefix offsets.

    mode 0: bucket = class key (iy0 << shift) | ix0           (nbins = N_buf^2)
    mode 1: bucket = (iy0 >> box) * ntx + (ix0 >> box)        (nbins = ntiles)
    """
    for i in range(nbins):
        cnt[i] = 0
    if mode == 0:
        for k in range(n):
            cnt[(iy0[k] << shift) | ix0[k]] += 1
    else:
        for k in range(n):
            cnt[((iy0[k] >> box_shift) * ntx) + (ix0[k] >> box_shift)] += 1
    s = 0
    for i in range(nbins):
        c = cnt[i]
        cnt[i] = s
        s += c


@njit(cache=False, fastmath=False, nogil=True)
def _scatter(ix0, iy0, w00, w10, w01, w11, cnt, key_dst, wdst, n,
             shift, box_shift, ntx, mode):
    """Place ray k at its bucket offset, in ascending k -> stable."""
    for k in range(n):
        key = (iy0[k] << shift) | ix0[k]
        if mode == 0:
            b = key
        else:
            b = ((iy0[k] >> box_shift) * ntx) + (ix0[k] >> box_shift)
        p = cnt[b]
        cnt[b] = p + 1
        key_dst[p] = key
        wdst[p, 0] = w00[k]
        wdst[p, 1] = w10[k]
        wdst[p, 2] = w01[k]
        wdst[p, 3] = w11[k]


# --------------------------------------------------------------------------
# sort 2: two-level counting sort, same output, cache-resident counters
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _radix_rows(ix0, iy0, w00, w10, w01, w11, cnt1, key_t, w_t, n,
                N_buf, shift):
    """MSD pass: bucket the block's rays by destination row iy0 (N_buf
    buckets, L1-resident counters), stable in ascending ray index.
    `cnt1` is left holding each row group's end offset."""
    for i in range(N_buf):
        cnt1[i] = 0
    for k in range(n):
        cnt1[iy0[k]] += 1
    s = 0
    for i in range(N_buf):
        c = cnt1[i]
        cnt1[i] = s
        s += c
    for k in range(n):
        y = iy0[k]
        p = cnt1[y]
        cnt1[y] = p + 1
        key_t[p] = (y << shift) | ix0[k]
        w_t[p, 0] = w00[k]
        w_t[p, 1] = w10[k]
        w_t[p, 2] = w01[k]
        w_t[p, 3] = w11[k]


@njit(cache=False, fastmath=False, nogil=True)
def _radix_cols(cnt1, cnt2, key_t, w_t, key_d, w_d, n, N_buf, mask):
    """LSD pass inside every row group: bucket by column ix0.  Reads and
    writes stay inside one group (a few KiB), so this pass runs out of L1.
    Output: the block sorted by (iy0, ix0), ascending ray index inside a
    class -- exactly what `_offsets` + `_scatter` produce."""
    a = 0
    for y in range(N_buf):
        b = cnt1[y]
        ng = b - a
        if ng > 0:
            # only the group's own column range needs counters: the group is
            # a few hundred rays, so this is ~8x less zeroing/prefix work
            # than sweeping all N_buf columns.
            cmin = N_buf
            cmax = -1
            for t in range(a, b):
                c = key_t[t] & mask
                if c < cmin:
                    cmin = c
                if c > cmax:
                    cmax = c
            for i in range(cmin, cmax + 1):
                cnt2[i] = 0
            for t in range(a, b):
                cnt2[key_t[t] & mask] += 1
            s = a
            for i in range(cmin, cmax + 1):
                c = cnt2[i]
                cnt2[i] = s
                s += c
            for t in range(a, b):
                key = key_t[t]
                p = cnt2[key & mask]
                cnt2[key & mask] = p + 1
                key_d[p] = key
                w_d[p, 0] = w_t[t, 0]
                w_d[p, 1] = w_t[t, 1]
                w_d[p, 2] = w_t[t, 2]
                w_d[p, 3] = w_t[t, 3]
        a = b


# --------------------------------------------------------------------------
# deposit: four bilinear corners of each class, one acc add per (bin, pass)
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _deposit_runs(acc, key_dst, wdst, n, N_buf, mask, shift, dx, dy, col):
    """Corner (dx,dy) of every class, closed once per class.

    Equivalent to `acc += np.bincount(i_corner, w_corner, minlength=nb)` for
    this block: one add per bin, weights summed in ascending ray order.
    """
    t = 0
    while t < n:
        key = key_dst[t]
        s = wdst[t, col]
        t += 1
        while t < n and key_dst[t] == key:
            s += wdst[t, col]
            t += 1
        ix = (key & mask) + dx
        if ix == N_buf:
            ix = 0
        iy = (key >> shift) + dy
        if iy == N_buf:
            iy = 0
        acc[iy * N_buf + ix] += s


@njit(cache=False, fastmath=False, nogil=True)
def _deposit_runs4(acc, key_dst, wdst, n, N_buf, mask, shift, st):
    """All four corners of every class, in one descending walk.

    Descending class order reaches bin (Y,X) as corner 0 of class (Y,X),
    corner 1 of class (Y,X-1), corner 2 of class (Y-1,X) and corner 3 of
    class (Y-1,X-1) -- production's pass order 1,2,3,4 -- so the four adds
    into one acc cell happen in production's order, with production's
    rounding after each add.

    The two cyclic seams (destination row 0 and column 0) cannot be ordered
    that way; those 2*N-1 bins are staged in `st` and flushed at the end, in
    the same 1,2,3,4 order.  Bin (Y,X) is staged in slot X if Y == 0, else in
    slot N_buf + Y (valid because the affected bins are exactly row 0 and
    column 0).
    """
    t = n - 1
    while t >= 0:
        key = key_dst[t]
        end = t + 1
        t -= 1
        while t >= 0 and key_dst[t] == key:
            t -= 1
        start = t + 1
        s0 = wdst[start, 0]
        s1 = wdst[start, 1]
        s2 = wdst[start, 2]
        s3 = wdst[start, 3]
        for u in range(start + 1, end):
            s0 += wdst[u, 0]
            s1 += wdst[u, 1]
            s2 += wdst[u, 2]
            s3 += wdst[u, 3]
        ix0 = key & mask
        iy0 = key >> shift
        ix1 = ix0 + 1
        if ix1 == N_buf:
            ix1 = 0
        iy1 = iy0 + 1
        if iy1 == N_buf:
            iy1 = 0
        # corner 0 -> (iy0, ix0)
        if iy0 == 0:
            st[ix0, 0] = s0
        elif ix0 == 0:
            st[N_buf + iy0, 0] = s0
        else:
            acc[iy0 * N_buf + ix0] += s0
        # corner 1 -> (iy0, ix1)
        if iy0 == 0:
            st[ix1, 1] = s1
        elif ix1 == 0:
            st[N_buf + iy0, 1] = s1
        else:
            acc[iy0 * N_buf + ix1] += s1
        # corner 2 -> (iy1, ix0)
        if iy1 == 0:
            st[ix0, 2] = s2
        elif ix0 == 0:
            st[N_buf + iy1, 2] = s2
        else:
            acc[iy1 * N_buf + ix0] += s2
        # corner 3 -> (iy1, ix1)
        if iy1 == 0:
            st[ix1, 3] = s3
        elif ix1 == 0:
            st[N_buf + iy1, 3] = s3
        else:
            acc[iy1 * N_buf + ix1] += s3
    # ---- flush the staged row 0 / column 0 bins in pass order ----
    # `st` is zero on entry (the caller keeps it zeroed), so an untouched
    # slot adds 0.0 -- which is exactly what production's `acc += zeros`
    # does for a bin that no ray hit in that pass.
    for j in range(N_buf):
        i = j
        acc[i] += st[i, 0]
        acc[i] += st[i, 1]
        acc[i] += st[i, 2]
        acc[i] += st[i, 3]
        st[i, 0] = 0.0
        st[i, 1] = 0.0
        st[i, 2] = 0.0
        st[i, 3] = 0.0
    for y in range(1, N_buf):
        i = y * N_buf
        acc[i] += st[N_buf + y, 0]
        acc[i] += st[N_buf + y, 1]
        acc[i] += st[N_buf + y, 2]
        acc[i] += st[N_buf + y, 3]
        st[N_buf + y, 0] = 0.0
        st[N_buf + y, 1] = 0.0
        st[N_buf + y, 2] = 0.0
        st[N_buf + y, 3] = 0.0


@njit(cache=False, fastmath=False, nogil=True)
def _deposit_tiles(acc, buf, gen, g, key_dst, wdst, n, N_buf, mask, shift,
                   dx, dy, col, box_shift, ntx, hit):
    """Corner (dx,dy) in tile order: production's scatter/commit arithmetic,
    but only a tile-sized window of `buf`/`gen` is touched at a time.

    A bin belongs to exactly one tile, so it is committed (one acc add) once
    per (pass, block) -- the same tree as production.
    """
    t = 0
    while t < n:
        tile = (((key_dst[t] >> shift) >> box_shift) * ntx
                + ((key_dst[t] & mask) >> box_shift))
        nh = 0
        while t < n:
            key = key_dst[t]
            tl = (((key >> shift) >> box_shift) * ntx
                  + ((key & mask) >> box_shift))
            if tl != tile:
                break
            ix = (key & mask) + dx
            if ix == N_buf:
                ix = 0
            iy = (key >> shift) + dy
            if iy == N_buf:
                iy = 0
            i = iy * N_buf + ix
            if gen[i] != g:
                gen[i] = g
                buf[i] = wdst[t, col]
                hit[nh] = i
                nh += 1
            else:
                buf[i] += wdst[t, col]
            t += 1
        for j in range(nh):
            i = hit[j]
            acc[i] += buf[i]
            buf[i] = 0.0
        g += 1
    return g


def _set_threads(n_threads):
    _df._set_threads(n_threads)


def _njit_check():
    if not HAVE_NUMBA or not _df.HAVE_NUMBA:
        raise RuntimeError("det_fast2 needs numba for the fused kernels")


def _pow2_shift(n):
    s = 0
    v = 1
    while v < n:
        v <<= 1
        s += 1
    return s if v == n else None


def det_frame_fast2(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
                    chunk_rows=64, work=None, n_threads=None,
                    sort_mode="radix", deposit="fused", tile=64, timing=None):
    """Same result as det.det_frame, with the deposit walked in bin order.

    sort_mode / deposit only choose how the order is built and walked; every
    combination here is bit-exact.  `timing` (dict) gets geom/sort/deposit
    seconds added in, so the net gain is visible stage by stage.
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

    shift = _pow2_shift(N_buf)
    box = int(tile)
    if sort_mode == "tile":
        deposit = "tiles"
    if shift is None or (deposit == "tiles"
                         and (box < 2 or (N_buf % box) != 0)):
        # generic fallback: production deposit, unsorted (never wrong, just
        # not localised).  Production is N_buf = 2048.
        C, info = _df.det_frame_fast(fields, S, N_buf, env, Q, ior=ior,
                                     row_block=row_block,
                                     chunk_rows=chunk_rows, work=work,
                                     n_threads=n_threads)
        if timing is not None:
            timing["fallback"] = True
        return C, info
    mask = (1 << shift) - 1
    box_shift = _pow2_shift(box) if deposit == "tiles" else 0
    ntx = (N_buf // box) if deposit == "tiles" else 0
    ntiles = ntx * ntx
    mode = 1 if sort_mode == "tile" else 0

    if work is None:
        work = {}
    if work.get("acc") is None or work["acc"].size != nb:
        work["acc"] = np.zeros(nb, dtype=np.float64)
    if work.get("ix0") is None or work["ix0"].size != n_rays:
        work["ix0"] = np.empty(n_rays, dtype=np.int32)
        work["iy0"] = np.empty(n_rays, dtype=np.int32)
        work["w00"] = np.empty(n_rays, dtype=np.float64)
        work["w10"] = np.empty(n_rays, dtype=np.float64)
        work["w01"] = np.empty(n_rays, dtype=np.float64)
        work["w11"] = np.empty(n_rays, dtype=np.float64)
    nblk = int(row_block) * N_s
    if work.get("key_dst") is None or work["key_dst"].size != nblk:
        work["key_dst"] = np.empty(nblk, dtype=np.int32)
        work["wdst"] = np.empty((nblk, 4), dtype=np.float64)
    if sort_mode == "radix" and (work.get("key_t") is None
                                 or work["key_t"].size != nblk):
        work["key_t"] = np.empty(nblk, dtype=np.int32)
        work["w_t"] = np.empty((nblk, 4), dtype=np.float64)
    if deposit == "fused" and (work.get("st") is None
                               or work["st"].shape[0] != 2 * N_buf):
        work["st"] = np.zeros((2 * N_buf, 4), dtype=np.float64)
    ncnt = nb if mode == 0 and sort_mode == "bin" else N_buf
    if mode == 1:
        ncnt = ntiles
    if work.get("cnt") is None or work["cnt"].size < ncnt:
        work["cnt"] = np.zeros(ncnt, dtype=np.int32)
    if sort_mode == "radix" and (work.get("cnt2") is None
                                 or work["cnt2"].size < N_buf):
        work["cnt2"] = np.zeros(N_buf, dtype=np.int32)
    if deposit == "tiles":
        if work.get("buf") is None or work["buf"].size != nb:
            work["buf"] = np.zeros(nb, dtype=np.float64)
            work["gen"] = np.zeros(nb, dtype=np.int32)
            work["hit"] = np.empty(box * box, dtype=np.int32)
            work["g"] = 1
        if int(work["g"]) > (1 << 30):
            work["gen"].fill(0)             # gen is 0 on entry, stamps >= 1
            work["g"] = 1

    acc = work["acc"]
    ix0 = work["ix0"]
    iy0 = work["iy0"]
    w00 = work["w00"]
    w10 = work["w10"]
    w01 = work["w01"]
    w11 = work["w11"]
    key_dst = work["key_dst"]
    wdst = work["wdst"]
    cnt = work["cnt"]

    acc.fill(0.0)                           # per-frame reset, as np.zeros

    t_geom = 0.0
    t_sort = 0.0
    t_sort2 = 0.0
    t_dep = 0.0
    for q in range(int(Q)):
        jx, jy = float(quads[q, 0]), float(quads[q, 1])
        d = s + u * jx + v * jy
        d = d / np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])

        t0 = time.perf_counter()
        _df._geom_all(h, gx, gy, Dx, Dy, N_s, S, s_over_ns, depth, cell,
                      N_buf, Ix, Iy, Iz, eta, int(chunk_rows),
                      ix0, iy0, w00, w10, w01, w11)
        t_geom += time.perf_counter() - t0

        for r0 in range(0, N_s, row_block):
            r1 = min(r0 + row_block, N_s)
            a = r0 * N_s
            b = r1 * N_s
            n = b - a

            t0 = time.perf_counter()
            if sort_mode == "radix":
                _radix_rows(ix0[a:b], iy0[a:b], w00[a:b], w10[a:b],
                            w01[a:b], w11[a:b], cnt, work["key_t"][:n],
                            work["w_t"][:n], n, N_buf, shift)
                t1 = time.perf_counter()
                _radix_cols(cnt, work["cnt2"], work["key_t"][:n],
                            work["w_t"][:n], key_dst[:n], wdst[:n], n,
                            N_buf, mask)
                t_sort2 += time.perf_counter() - t1
            else:
                _offsets(ix0[a:b], iy0[a:b], cnt, n, ncnt, shift, box_shift,
                         ntx, mode)
                _scatter(ix0[a:b], iy0[a:b], w00[a:b], w10[a:b], w01[a:b],
                         w11[a:b], cnt, key_dst[:n], wdst[:n], n, shift,
                         box_shift, ntx, mode)
            t_sort += time.perf_counter() - t0

            t0 = time.perf_counter()
            if deposit == "fused":
                _deposit_runs4(acc, key_dst[:n], wdst[:n], n, N_buf, mask,
                               shift, work["st"])
            elif deposit == "pass4":
                _deposit_runs(acc, key_dst[:n], wdst[:n], n, N_buf, mask,
                              shift, 0, 0, 0)
                _deposit_runs(acc, key_dst[:n], wdst[:n], n, N_buf, mask,
                              shift, 1, 0, 1)
                _deposit_runs(acc, key_dst[:n], wdst[:n], n, N_buf, mask,
                              shift, 0, 1, 2)
                _deposit_runs(acc, key_dst[:n], wdst[:n], n, N_buf, mask,
                              shift, 1, 1, 3)
            else:
                g = int(work["g"])
                for col, (dx, dy) in enumerate(((0, 0), (1, 0), (0, 1),
                                                (1, 1))):
                    g = _deposit_tiles(acc, work["buf"], work["gen"], g,
                                       key_dst[:n], wdst[:n], n, N_buf, mask,
                                       shift, dx, dy, col, box_shift, ntx,
                                       work["hit"])
                work["g"] = g
            t_dep += time.perf_counter() - t0

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
    if timing is not None:
        timing["geom_s"] = timing.get("geom_s", 0.0) + t_geom
        timing["sort_s"] = timing.get("sort_s", 0.0) + t_sort
        timing["sort_msd_s"] = timing.get("sort_msd_s", 0.0) + (t_sort - t_sort2)
        timing["sort_lsd_s"] = timing.get("sort_lsd_s", 0.0) + t_sort2
        timing["deposit_s"] = timing.get("deposit_s", 0.0) + t_dep
        timing["sort_mode"] = sort_mode
        timing["deposit"] = deposit
        timing["tile"] = box
        timing["blocks"] = int(Q) * ((N_s + row_block - 1) // row_block)
    return C, info


def warmup(n_threads=None):
    """Compile every kernel once (tiny data) so timing excludes JIT cost."""
    _njit_check()
    _set_threads(n_threads)
    n = 32
    nbuf = 64
    z = np.zeros((n, n), dtype=np.float64)
    nn = n * n
    acc = np.zeros(nbuf * nbuf, dtype=np.float64)
    ix0 = np.zeros(nn, dtype=np.int32)
    iy0 = np.zeros(nn, dtype=np.int32)
    w = np.full(nn, 0.25, dtype=np.float64)
    _df._geom_all(z, z, z, z, z, n, 1.0, 1.0 / n, 0.6, 1.0 / n, nbuf,
                  0.0, 0.0, -1.0, 1.0, 8,
                  ix0, iy0, w, w, w, w)
    shift = _pow2_shift(nbuf)
    mask = (1 << shift) - 1
    key_dst = np.zeros(nn, dtype=np.int32)
    wdst = np.zeros((nn, 4), dtype=np.float64)
    # single-level
    cnt = np.zeros(nbuf * nbuf, dtype=np.int32)
    _offsets(ix0, iy0, cnt, nn, nbuf * nbuf, shift, 0, 0, 0)
    _scatter(ix0, iy0, w, w, w, w, cnt, key_dst, wdst, nn, shift, 0, 0, 0)
    _deposit_runs(acc, key_dst, wdst, nn, nbuf, mask, shift, 0, 0, 0)
    # radix
    cnt1 = np.zeros(nbuf, dtype=np.int32)
    cnt2 = np.zeros(nbuf, dtype=np.int32)
    key_t = np.zeros(nn, dtype=np.int32)
    w_t = np.zeros((nn, 4), dtype=np.float64)
    _radix_rows(ix0, iy0, w, w, w, w, cnt1, key_t, w_t, nn, nbuf, shift)
    _radix_cols(cnt1, cnt2, key_t, w_t, key_dst, wdst, nn, nbuf, mask)
    st = np.zeros((2 * nbuf, 4), dtype=np.float64)
    _deposit_runs4(acc, key_dst, wdst, nn, nbuf, mask, shift, st)
    # tiles
    box_shift = 2                     # tile = 4 on a 64 grid
    ntx = nbuf // 4
    cnt3 = np.zeros(ntx * ntx, dtype=np.int32)
    _offsets(ix0, iy0, cnt3, nn, ntx * ntx, shift, box_shift, ntx, 1)
    _scatter(ix0, iy0, w, w, w, w, cnt3, key_dst, wdst, nn, shift,
             box_shift, ntx, 1)
    buf = np.zeros(nbuf * nbuf, dtype=np.float64)
    gen = np.zeros(nbuf * nbuf, dtype=np.int32)
    hit = np.empty(4 * 4, dtype=np.int32)
    _deposit_tiles(acc, buf, gen, 1, key_dst, wdst, nn, nbuf, mask, shift,
                   0, 0, 0, box_shift, ntx, hit)


if __name__ == "__main__":
    import argparse
    import json
    import sys

    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from env import resolve_environment
    from state import WaterState
    from det_fast import det_frame_fast

    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", type=int, default=0)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--sort-mode", default="radix",
                    choices=["radix", "bin", "tile"])
    ap.add_argument("--deposit", default="fused",
                    choices=["fused", "pass4", "tiles"])
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
    _df.warmup(args.threads)
    t_jit = time.perf_counter() - t0
    tm = {}
    t0 = time.perf_counter()
    C2, info2 = det_frame_fast2(fields, state.S, 2048, env, 25,
                                n_threads=args.threads,
                                sort_mode=args.sort_mode,
                                deposit=args.deposit, timing=tm)
    t2 = time.perf_counter() - t0
    t0 = time.perf_counter()
    C1, info1 = det_frame_fast(fields, state.S, 2048, env, 25,
                               n_threads=args.threads)
    t1 = time.perf_counter() - t0
    d = np.abs(C2.astype(np.float64) - C1.astype(np.float64))
    print("jit=%.2fs fast2=%.3fs fast=%.3fs ratio=%.3f"
          % (t_jit, t2, t1, t1 / t2))
    print("timing:", json.dumps(tm, indent=1))
    print("bit_identical=%s max_abs_diff=%.3g info_equal=%s"
          % (np.array_equal(C2, C1), float(d.max()), info1 == info2))
