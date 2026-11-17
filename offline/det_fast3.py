"""det_fast3 -- det_frame with the sort cost cut (NEW file, BRIEF-DET-SORT2).

det_fast2 (radix + fused, wired into production by cf-0063) spends 70 % of
the deposit stage on sorting: 2.311 s of 3.306 s, split MSD 1.079 s +
LSD 1.232 s (REPORT-DET-BINSORT ﾂｧ4.1).  This module attacks that sort.

The two changes the brief asks for
----------------------------------
(a) FINER MSD BUCKET -- bucket the block's rays by `(iy0 << kh) | (ix0 >> kl)`
    instead of by destination row alone.  A det_fast2 row group holds ~256
    rays whose destination columns are spread over all 2048 columns; with
    kh = 6 a bucket holds ~5.9 rays that also share their high column bits
    (measured: p99 = 20, max = 43, numbers-probe.json).  With kh = 0 the
    bucket is exactly det_fast2's row bucket, so the two are directly
    comparable.

(b) THE LSD IS ELIMINATED, not just shortened.  det_fast2's LSD does three
    scans per row group (min/max, count, scatter) and sweeps a 2048-wide
    counter window for every row group -- ~8.4 M counter operations per
    block whatever the ray count.  det_fast3 replaces it with a **stamped
    per-bucket accumulator**: every ray of a bucket is accumulated straight
    into `sacc[col_lo]` (one slot per low column, L1-resident), and the
    bucket is then emitted class by class.  There is no second scatter, no
    second block-sized buffer and no weight movement beyond the MSD:

        MSD scatter -> for each bucket: acc[col_lo] += weights (in ray order)
                    -> emit the bucket's classes, descending column order

    Why this is still bit-identical: a class is a single (iy0, ix0), so the
    rays of a class all live in one bucket and are visited in ascending ray
    index (the MSD scatter is stable).  `sacc` therefore holds exactly the
    per-class sums s0..s3 in ascending ray order, i.e. production's
    `np.bincount` sums.  Emitting the buckets in descending (row, col_hi)
    order and the classes inside a bucket in descending column order reaches
    bin (Y,X) as corner 0 of class (Y,X), then corner 1 of class (Y,X-1),
    then corner 2 of class (Y-1,X), then corner 3 of class (Y-1,X-1) --
    production's 1,2,3,4 order -- so `_emit_class` can add the four corners
    straight into `acc` exactly as det_fast2's `_deposit_runs4` does.  The
    cyclic seams (destination row 0 / column 0) keep det_fast2's staged
    `st` buffer and its flush order.

What is NOT changed
-------------------
  * geometry: `det_fast._geom_all` imported and called unchanged (fmod,
    clip order, divisor, FFT all as production);
  * `row_block = 256` (it defines the summation grouping);
  * the deposit statements (`_emit_class`, `_flush_st`) are det_fast2's
    `_deposit_runs4` statements, unchanged;
  * the MSD scatter still writes (key, w00, w10, w01, w11) = 36 B per ray,
    so (a) is measured on its own merits.

Variants kept for the measurement (`mode=`)
-------------------------------------------
  "acc"    (a)+(b): pristine accumulator, touched slots found by sweeping the
           bucket's whole low-column window (gen-stamped)
  "mask"   (a)+(b): same arithmetic, touched slots tracked in a bitset so the
           cost is proportional to the classes present
  "count"  (a) only: fine bucket + det_fast2's min/max/count/scatter LSD
  "micro"  (b) alternative: fine bucket + stable insertion sort per bucket
  "fused"  (b) alternative: fine bucket + det_fast2's LSD into an L1 scratch
           that the deposit consumes immediately

Measured outcome (REPORT-DET-SORT2): none of them beats det_fast2's
radix+fused on this machine.  The two reasons are quantified in the report;
short version: (a) makes the MSD scatter *slower* (its random-write bytes do
not change, but the counter working set grows from 8 KiB to 128 KiB-1 MiB),
and (b) cannot beat det_fast2's LSD because any per-class ordering work costs
`nbuckets * nslot` = 4.19 M counter slots per block, or - for the swept-free
accumulator - a per-ray accumulate that costs about the same as the LSD it
replaces.  Kept in the tree as the measured evidence, not as a candidate.

Public API
----------
det_frame_fast3(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
                chunk_rows=64, work=None, n_threads=None,
                kh=6, mode="acc", micro_max=64, timing=None)
    -> (C float32 (N_buf,N_buf), info dict)      # same as det.det_frame
warmup(n_threads=None)   -- compile once; timing must exclude this
"""
import os
import time

import numpy as np

from det import IOR_G, SIGMA_PX, fibonacci_disc
from transport import orthonormal_basis, sun_disc_radius, gaussian_blur_fft

import det_fast as _df
import det_fast2 as _df2

try:
    import numba
    from numba import njit
    HAVE_NUMBA = True
except Exception:                       # pragma: no cover
    HAVE_NUMBA = False
    numba = None

    def njit(*a, **kw):
        def deco(f):
            return f
        if a and callable(a[0]):
            return a[0]
        return deco


# --------------------------------------------------------------------------
# MSD: count + prefix on the composite bucket (row, high column bits)
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _msd_count(ix0, iy0, cnt, n, nbuckets, kh, kl):
    """Count the block's rays per bucket, then exclusive prefix offsets."""
    for i in range(nbuckets):
        cnt[i] = 0
    for k in range(n):
        cnt[(iy0[k] << kh) | (ix0[k] >> kl)] += 1
    s = 0
    for i in range(nbuckets):
        c = cnt[i]
        cnt[i] = s
        s += c


@njit(cache=False, fastmath=False, nogil=True)
def _msd_scatter(ix0, iy0, w00, w10, w01, w11, cnt, key_d, w_d, n,
                 shift, kh, kl):
    """Place ray k at its bucket offset, in ascending k -> stable."""
    for k in range(n):
        iy = iy0[k]
        ix = ix0[k]
        b = (iy << kh) | (ix >> kl)
        p = cnt[b]
        cnt[b] = p + 1
        key_d[p] = (iy << shift) | ix
        w_d[p, 0] = w00[k]
        w_d[p, 1] = w10[k]
        w_d[p, 2] = w01[k]
        w_d[p, 3] = w11[k]


# --------------------------------------------------------------------------
# deposit: statements of det_fast2._deposit_runs4, factored per class
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _emit_class(acc, key, s0, s1, s2, s3, N_buf, mask, shift, st):
    """The four corners of one class, exactly as det_fast2 writes them."""
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


@njit(cache=False, fastmath=False, nogil=True)
def _flush_st(acc, N_buf, st):
    """The staged cyclic seams (destination row 0 / column 0), in pass
    order; every slot is reset so an untouched slot adds exactly 0.0."""
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
def _deposit_range(acc, key_d, w_d, a, e, N_buf, mask, shift, st):
    """Descending walk of a sorted range, one add per class (det_fast2's
    run detection, unchanged)."""
    t = e - 1
    while t >= a:
        key = key_d[t]
        end = t + 1
        t -= 1
        while t >= a and key_d[t] == key:
            t -= 1
        start = t + 1
        s0 = w_d[start, 0]
        s1 = w_d[start, 1]
        s2 = w_d[start, 2]
        s3 = w_d[start, 3]
        for u in range(start + 1, end):
            s0 += w_d[u, 0]
            s1 += w_d[u, 1]
            s2 += w_d[u, 2]
            s3 += w_d[u, 3]
        _emit_class(acc, key, s0, s1, s2, s3, N_buf, mask, shift, st)


# --------------------------------------------------------------------------
# (a)+(b): stamped per-bucket accumulator, then emit class by class
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _acc_emit2(acc, key_d, w_d, ends, nbuckets, N_buf, mask, shift, kh, kl,
               sacc, gen, g0, st):
    """Walk the fine buckets descending; accumulate each bucket's rays into
    `sacc[col_lo]` and emit the classes in descending column order.

    The rays of a class all live in one bucket and are visited in ascending
    ray index (the MSD scatter is stable), so `sacc` holds production's
    per-class sums s0..s3.  `gen` is a generation stamp per low-column slot,
    so no clearing is needed between buckets; `g0` is returned.
    """
    nslot = 1 << kl
    slmask = nslot - 1
    g = g0
    for b in range(nbuckets - 1, -1, -1):
        a = 0 if b == 0 else ends[b - 1]
        e = ends[b]
        if e == a:
            continue
        g += 1
        for t in range(a, e):
            c = key_d[t] & slmask
            if gen[c] != g:
                gen[c] = g
                sacc[c, 0] = w_d[t, 0]
                sacc[c, 1] = w_d[t, 1]
                sacc[c, 2] = w_d[t, 2]
                sacc[c, 3] = w_d[t, 3]
            else:
                sacc[c, 0] += w_d[t, 0]
                sacc[c, 1] += w_d[t, 1]
                sacc[c, 2] += w_d[t, 2]
                sacc[c, 3] += w_d[t, 3]
        row = b >> kh
        chi = b & ((1 << kh) - 1)
        base = chi << kl
        krow = row << shift
        for c in range(nslot - 1, -1, -1):
            if gen[c] == g:
                _emit_class(acc, krow | (base | c), sacc[c, 0], sacc[c, 1],
                            sacc[c, 2], sacc[c, 3], N_buf, mask, shift, st)
    _flush_st(acc, N_buf, st)
    return g


# --------------------------------------------------------------------------
# (a) only: fine bucket + det_fast2's min/max/count/scatter LSD
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _lsd_count(ends, key_t, w_t, key_d, w_d, nbuckets, mask_lo, cnt2):
    a = 0
    for b in range(nbuckets):
        e = ends[b]
        if e - a == 1:
            key_d[a] = key_t[a]
            w_d[a, 0] = w_t[a, 0]
            w_d[a, 1] = w_t[a, 1]
            w_d[a, 2] = w_t[a, 2]
            w_d[a, 3] = w_t[a, 3]
        elif e - a > 1:
            cmin = mask_lo + 1
            cmax = -1
            for t in range(a, e):
                c = key_t[t] & mask_lo
                if c < cmin:
                    cmin = c
                if c > cmax:
                    cmax = c
            for i in range(cmin, cmax + 1):
                cnt2[i] = 0
            for t in range(a, e):
                cnt2[key_t[t] & mask_lo] += 1
            s = a
            for i in range(cmin, cmax + 1):
                c = cnt2[i]
                cnt2[i] = s
                s += c
            for t in range(a, e):
                key = key_t[t]
                c = key & mask_lo
                p = cnt2[c]
                cnt2[c] = p + 1
                key_d[p] = key
                w_d[p, 0] = w_t[t, 0]
                w_d[p, 1] = w_t[t, 1]
                w_d[p, 2] = w_t[t, 2]
                w_d[p, 3] = w_t[t, 3]
        a = e


# --------------------------------------------------------------------------
# (b) alternative kept on record: stable insertion sort per fine bucket
# --------------------------------------------------------------------------
@njit(cache=False, fastmath=False, nogil=True)
def _isort_range(key_d, w_d, a, e, mask_lo):
    """Stable insertion sort of [a, e) by (key & mask_lo)."""
    for i in range(a + 1, e):
        key = key_d[i]
        lo = key & mask_lo
        w0 = w_d[i, 0]
        w1 = w_d[i, 1]
        w2 = w_d[i, 2]
        w3 = w_d[i, 3]
        j = i - 1
        while j >= a and (key_d[j] & mask_lo) > lo:
            key_d[j + 1] = key_d[j]
            w_d[j + 1, 0] = w_d[j, 0]
            w_d[j + 1, 1] = w_d[j, 1]
            w_d[j + 1, 2] = w_d[j, 2]
            w_d[j + 1, 3] = w_d[j, 3]
            j -= 1
        key_d[j + 1] = key
        w_d[j + 1, 0] = w0
        w_d[j + 1, 1] = w1
        w_d[j + 1, 2] = w2
        w_d[j + 1, 3] = w3


@njit(cache=False, fastmath=False, nogil=True)
def _csort_range(key_d, w_d, a, e, mask_lo, cnt2, tk, tw):
    """Stable counting sort of [a, e) by (key & mask_lo) into tk/tw[0:r]."""
    cmin = mask_lo + 1
    cmax = -1
    for t in range(a, e):
        c = key_d[t] & mask_lo
        if c < cmin:
            cmin = c
        if c > cmax:
            cmax = c
    for i in range(cmin, cmax + 1):
        cnt2[i] = 0
    for t in range(a, e):
        cnt2[key_d[t] & mask_lo] += 1
    s = 0
    for i in range(cmin, cmax + 1):
        c = cnt2[i]
        cnt2[i] = s
        s += c
    for t in range(a, e):
        key = key_d[t]
        c = key & mask_lo
        p = cnt2[c]
        cnt2[c] = p + 1
        tk[p] = key
        tw[p, 0] = w_d[t, 0]
        tw[p, 1] = w_d[t, 1]
        tw[p, 2] = w_d[t, 2]
        tw[p, 3] = w_d[t, 3]


@njit(cache=False, fastmath=False, nogil=True)
def _acc_mask(acc, key_d, w_d, ends, nbuckets, N_buf, mask, shift, kh, kl,
              sacc, words, hi_tab, st):
    """(a)+(b) final form: stamped accumulator with a bitmask emit.

    Identical arithmetic to `_acc_emit2`, but the touched low-column slots
    are tracked in a bitset instead of being found by sweeping the whole
    window, so the cost is proportional to the classes actually present
    (`_acc_emit2` sweeps nbuckets * nslot = 4.19 M slots per block, which is
    what makes it a wash against det_fast2's LSD).
    """
    nslot = 1 << kl
    nwords = (nslot + 63) >> 6
    slmask = nslot - 1
    for b in range(nbuckets - 1, -1, -1):
        a = 0 if b == 0 else ends[b - 1]
        e = ends[b]
        if e == a:
            continue
        for i in range(nwords):
            words[i] = 0
        for t in range(a, e):
            c = key_d[t] & slmask
            wi = c >> 6
            bit = np.int64(1) << (c & 63)
            if words[wi] & bit:
                sacc[c, 0] += w_d[t, 0]
                sacc[c, 1] += w_d[t, 1]
                sacc[c, 2] += w_d[t, 2]
                sacc[c, 3] += w_d[t, 3]
            else:
                words[wi] |= bit
                sacc[c, 0] = w_d[t, 0]
                sacc[c, 1] = w_d[t, 1]
                sacc[c, 2] = w_d[t, 2]
                sacc[c, 3] = w_d[t, 3]
        row = b >> kh
        base = (b & ((1 << kh) - 1)) << kl
        krow = row << shift
        wi = nwords - 1
        while wi >= 0:
            w = words[wi]
            while w != 0:
                sh = 56
                nb = 0
                while sh >= 0:
                    by = (w >> sh) & 255
                    if by != 0:
                        nb = sh + hi_tab[by]
                        break
                    sh -= 8
                c = (wi << 6) | nb
                _emit_class(acc, krow | (base | c), sacc[c, 0], sacc[c, 1],
                            sacc[c, 2], sacc[c, 3], N_buf, mask, shift, st)
                w &= ~(np.int64(1) << nb)
            wi -= 1
    _flush_st(acc, N_buf, st)


@njit(cache=False, fastmath=False, nogil=True)
def _max_bucket(ends, nbuckets):
    m = 0
    a = 0
    for b in range(nbuckets):
        e = ends[b]
        if e - a > m:
            m = e - a
        a = e
    return m


@njit(cache=False, fastmath=False, nogil=True)
def _fused_count_deposit(acc, key_t, w_t, ends, nbuckets, N_buf, mask, shift,
                         mask_lo, cnt2, tk, tw, st):
    """(b) fused: sort each bucket into an L1 scratch and deposit it at once.

    The three LSD scans (min/max, count, scatter) are unchanged, but the
    scatter target is a small scratch that the deposit reads back while it
    is still in L1/L2 -- det_fast2 instead writes a second 18.9 MiB block
    buffer and reads it back from memory.  Buckets are visited descending,
    so the emitted corner order is still production's 1,2,3,4.
    """
    for b in range(nbuckets - 1, -1, -1):
        a = 0 if b == 0 else ends[b - 1]
        e = ends[b]
        r = e - a
        if r == 0:
            continue
        if r == 1:
            _emit_class(acc, key_t[a], w_t[a, 0], w_t[a, 1], w_t[a, 2],
                        w_t[a, 3], N_buf, mask, shift, st)
            continue
        cmin = mask_lo + 1
        cmax = -1
        for t in range(a, e):
            c = key_t[t] & mask_lo
            if c < cmin:
                cmin = c
            if c > cmax:
                cmax = c
        for i in range(cmin, cmax + 1):
            cnt2[i] = 0
        for t in range(a, e):
            cnt2[key_t[t] & mask_lo] += 1
        s = 0
        for i in range(cmin, cmax + 1):
            c = cnt2[i]
            cnt2[i] = s
            s += c
        for t in range(a, e):
            key = key_t[t]
            c = key & mask_lo
            p = cnt2[c]
            cnt2[c] = p + 1
            tk[p] = key
            tw[p, 0] = w_t[t, 0]
            tw[p, 1] = w_t[t, 1]
            tw[p, 2] = w_t[t, 2]
            tw[p, 3] = w_t[t, 3]
        _deposit_range(acc, tk, tw, 0, r, N_buf, mask, shift, st)
    _flush_st(acc, N_buf, st)


@njit(cache=False, fastmath=False, nogil=True)
def _sortdep_fused(acc, key_d, w_d, ends, nbuckets, n, N_buf, mask, shift,
                   mask_lo, micro_max, cnt2, tk, tw, st):
    """Sort + deposit in one descending walk (insertion sort variant)."""
    for b in range(nbuckets - 1, -1, -1):
        a = 0 if b == 0 else ends[b - 1]
        e = ends[b]
        if e > a:
            if e - a <= micro_max:
                _isort_range(key_d, w_d, a, e, mask_lo)
                _deposit_range(acc, key_d, w_d, a, e, N_buf, mask, shift, st)
            else:
                _csort_range(key_d, w_d, a, e, mask_lo, cnt2, tk, tw)
                _deposit_range(acc, tk, tw, 0, e - a, N_buf, mask, shift, st)
    _flush_st(acc, N_buf, st)


@njit(cache=False, fastmath=False, nogil=True)
def _sort_only(ends, key_d, w_d, nbuckets, mask_lo, micro_max, cnt2, tk, tw):
    """Sort every bucket in place (oversized ones through tk/tw and back)."""
    a = 0
    for b in range(nbuckets):
        e = ends[b]
        if e - a > 1:
            if e - a <= micro_max:
                _isort_range(key_d, w_d, a, e, mask_lo)
            else:
                r = e - a
                _csort_range(key_d, w_d, a, e, mask_lo, cnt2, tk, tw)
                for t in range(r):
                    key_d[a + t] = tk[t]
                    w_d[a + t, 0] = tw[t, 0]
                    w_d[a + t, 1] = tw[t, 1]
                    w_d[a + t, 2] = tw[t, 2]
                    w_d[a + t, 3] = tw[t, 3]
        a = e


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def _set_threads(n_threads):
    _df._set_threads(n_threads)


def _njit_check():
    if not HAVE_NUMBA or not _df.HAVE_NUMBA:
        raise RuntimeError("det_fast3 needs numba for the fused kernels")


def _pow2_shift(n):
    s = 0
    v = 1
    while v < n:
        v <<= 1
        s += 1
    return s if v == n else None


def det_frame_fast3(fields, S, N_buf, env, Q, ior=IOR_G, row_block=256,
                    chunk_rows=64, work=None, n_threads=None,
                    kh=6, mode="acc", micro_max=64, timing=None):
    """Same result as det.det_frame, with the sort cost cut.

    kh    high column bits folded into the MSD bucket (0 = det_fast2's row
          bucket, which is directly comparable).
    mode  "acc" (candidate: stamped accumulator), "count" ((a) only with
          det_fast2's LSD), "micro" (insertion sort per bucket, kept on
          record as the rejected alternative).
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
    kh = int(kh)
    mode = str(mode).lower()
    if shift is None or kh < 0 or kh >= shift or mode not in (
            "acc", "mask", "count", "micro", "fused"):
        C, info = _df2.det_frame_fast2(fields, S, N_buf, env, Q, ior=ior,
                                       row_block=row_block,
                                       chunk_rows=chunk_rows, work=work,
                                       n_threads=n_threads, timing=timing)
        if timing is not None:
            timing["fallback"] = True
        return C, info
    kl = shift - kh
    mask = (1 << shift) - 1
    mask_lo = (1 << kl) - 1
    nbuckets = N_buf << kh
    nslot = 1 << kl

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
    if work.get("key_t") is None or work["key_t"].size != nblk:
        work["key_t"] = np.empty(nblk, dtype=np.int32)
        work["w_t"] = np.empty((nblk, 4), dtype=np.float64)
    if work.get("cnt") is None or work["cnt"].size < nbuckets:
        work["cnt"] = np.zeros(nbuckets, dtype=np.int32)
    if work.get("cnt2") is None or work["cnt2"].size < nslot:
        work["cnt2"] = np.zeros(max(nslot, 8), dtype=np.int32)
    if work.get("st") is None or work["st"].shape[0] != 2 * N_buf:
        work["st"] = np.zeros((2 * N_buf, 4), dtype=np.float64)
    if work.get("sacc") is None or work["sacc"].shape[0] != nslot:
        work["sacc"] = np.zeros((nslot, 4), dtype=np.float64)
        work["gen"] = np.zeros(nslot, dtype=np.int64)
        work["g"] = 1
    if int(work["g"]) > (1 << 60):
        work["gen"].fill(0)
        work["g"] = 1
    if mode in ("count", "micro") and (work.get("key_d") is None
                                       or work["key_d"].size != nblk):
        work["key_d"] = np.empty(nblk, dtype=np.int32)
        work["w_d"] = np.empty((nblk, 4), dtype=np.float64)
    if work.get("tk") is None or work["tk"].size < 64:
        work["tk"] = np.empty(64, dtype=np.int32)
        work["tw"] = np.empty((64, 4), dtype=np.float64)
    if mode == "mask":
        if work.get("words") is None or work["words"].size < (nslot + 63) // 64:
            work["words"] = np.zeros((nslot + 63) // 64, dtype=np.int64)
        if work.get("hi_tab") is None:
            tab = np.zeros(256, dtype=np.int64)
            for i in range(256):
                j = 0
                while (1 << j) <= i:
                    if i & (1 << j):
                        tab[i] = j
                    j += 1
            work["hi_tab"] = tab

    acc = work["acc"]
    ix0 = work["ix0"]
    iy0 = work["iy0"]
    w00 = work["w00"]
    w10 = work["w10"]
    w01 = work["w01"]
    w11 = work["w11"]
    key_t = work["key_t"]
    w_t = work["w_t"]
    cnt = work["cnt"]
    st = work["st"]

    acc.fill(0.0)

    t_geom = 0.0
    t_msd = 0.0
    t_lsd = 0.0
    t_dep = 0.0
    t_acc = 0.0
    n_grow = 0
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
            _msd_count(ix0[a:b], iy0[a:b], cnt, n, nbuckets, kh, kl)
            _msd_scatter(ix0[a:b], iy0[a:b], w00[a:b], w10[a:b], w01[a:b],
                         w11[a:b], cnt, key_t[:n], w_t[:n], n, shift, kh, kl)
            t_msd += time.perf_counter() - t0

            if mode == "acc":
                t0 = time.perf_counter()
                work["g"] = _acc_emit2(acc, key_t[:n], w_t[:n], cnt, nbuckets,
                                       N_buf, mask, shift, kh, kl,
                                       work["sacc"], work["gen"],
                                       int(work["g"]), st)
                t_acc += time.perf_counter() - t0
            elif mode == "mask":
                t0 = time.perf_counter()
                _acc_mask(acc, key_t[:n], w_t[:n], cnt, nbuckets, N_buf, mask,
                          shift, kh, kl, work["sacc"], work["words"],
                          work["hi_tab"], st)
                t_acc += time.perf_counter() - t0
            elif mode == "count":
                t0 = time.perf_counter()
                _lsd_count(cnt, key_t[:n], w_t[:n], work["key_d"][:n],
                           work["w_d"][:n], nbuckets, mask_lo, work["cnt2"])
                t_lsd += time.perf_counter() - t0
                t0 = time.perf_counter()
                _df2._deposit_runs4(acc, work["key_d"][:n], work["w_d"][:n],
                                    n, N_buf, mask, shift, st)
                t_dep += time.perf_counter() - t0
            elif mode == "fused":
                need = int(_max_bucket(cnt, nbuckets))
                if need > work["tk"].size:
                    cap = 64
                    while cap < need:
                        cap *= 2
                    work["tk"] = np.empty(cap, dtype=np.int32)
                    work["tw"] = np.empty((cap, 4), dtype=np.float64)
                    n_grow += 1
                t0 = time.perf_counter()
                _fused_count_deposit(acc, key_t[:n], w_t[:n], cnt, nbuckets,
                                     N_buf, mask, shift, mask_lo,
                                     work["cnt2"], work["tk"], work["tw"], st)
                t_lsd += time.perf_counter() - t0
            else:                                   # "micro"
                need = int(_max_bucket(cnt, nbuckets))
                if need > work["tk"].size:
                    cap = 64
                    while cap < need:
                        cap *= 2
                    work["tk"] = np.empty(cap, dtype=np.int32)
                    work["tw"] = np.empty((cap, 4), dtype=np.float64)
                    n_grow += 1
                t0 = time.perf_counter()
                _sortdep_fused(acc, key_t[:n], w_t[:n], cnt, nbuckets, n,
                               N_buf, mask, shift, mask_lo, int(micro_max),
                               work["cnt2"], work["tk"], work["tw"], st)
                t_lsd += time.perf_counter() - t0

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
        timing["sort_msd_s"] = timing.get("sort_msd_s", 0.0) + t_msd
        timing["sort_lsd_s"] = timing.get("sort_lsd_s", 0.0) + t_lsd
        timing["acc_emit_s"] = timing.get("acc_emit_s", 0.0) + t_acc
        timing["deposit_s"] = timing.get("deposit_s", 0.0) + t_dep
        timing["sort_or_acc_s"] = t_msd + t_lsd + t_acc
        timing["kh"] = kh
        timing["kl"] = kl
        timing["mode"] = mode
        timing["micro_max"] = int(micro_max)
        timing["blocks"] = int(Q) * ((N_s + row_block - 1) // row_block)
        timing["micro_grows"] = int(n_grow)
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
    kh, kl = 2, shift - 2
    mask_lo = (1 << kl) - 1
    nslot = 1 << kl
    nbuckets = nbuf << kh
    cnt = np.zeros(nbuckets, dtype=np.int32)
    key_t = np.zeros(nn, dtype=np.int32)
    w_t = np.zeros((nn, 4), dtype=np.float64)
    cnt2 = np.zeros(nslot, dtype=np.int32)
    tk = np.zeros(nn, dtype=np.int32)
    tw = np.zeros((nn, 4), dtype=np.float64)
    st = np.zeros((2 * nbuf, 4), dtype=np.float64)
    sacc = np.zeros((nslot, 4), dtype=np.float64)
    gen = np.zeros(nslot, dtype=np.int64)
    _msd_count(ix0, iy0, cnt, nn, nbuckets, kh, kl)
    _msd_scatter(ix0, iy0, w, w, w, w, cnt, key_t, w_t, nn, shift, kh, kl)
    _acc_emit2(acc, key_t, w_t, cnt, nbuckets, nbuf, mask, shift, kh, kl,
               sacc, gen, 1, st)
    _flush_st(acc, nbuf, st)
    _lsd_count(cnt, key_t, w_t, key_t, w_t, nbuckets, mask_lo, cnt2)
    _isort_range(key_t, w_t, 0, nn, mask_lo)
    _csort_range(key_t, w_t, 0, nn, mask_lo, cnt2, tk, tw)
    _sort_only(cnt, key_t, w_t, nbuckets, mask_lo, 32, cnt2, tk, tw)
    _deposit_range(acc, key_t, w_t, 0, nn, nbuf, mask, shift, st)
    _df2._deposit_runs4(acc, key_t, w_t, nn, nbuf, mask, shift, st)
    _sortdep_fused(acc, key_t, w_t, cnt, nbuckets, nn, nbuf, mask, shift,
                   mask_lo, 32, cnt2, tk, tw, st)
    _max_bucket(cnt, nbuckets)
    _fused_count_deposit(acc, key_t, w_t, cnt, nbuckets, nbuf, mask, shift,
                         mask_lo, cnt2, tk, tw, st)
    _flush_st(acc, nbuf, st)


if __name__ == "__main__":
    import argparse
    import json
    import sys

    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from env import resolve_environment
    from state import WaterState
    from det_fast2 import det_frame_fast2

    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", type=int, default=0)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--kh", type=int, default=6)
    ap.add_argument("--mode", default="acc",
                    choices=["acc", "mask", "count", "micro", "fused"])
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
    _df2.warmup(args.threads)
    t_jit = time.perf_counter() - t0
    tm = {}
    t0 = time.perf_counter()
    C3, info3 = det_frame_fast3(fields, state.S, 2048, env, 25,
                                n_threads=args.threads, kh=args.kh,
                                mode=args.mode, timing=tm)
    t3 = time.perf_counter() - t0
    t0 = time.perf_counter()
    C2, info2 = det_frame_fast2(fields, state.S, 2048, env, 25,
                                n_threads=args.threads)
    t2 = time.perf_counter() - t0
    d = np.abs(C3.astype(np.float64) - C2.astype(np.float64))
    print("jit=%.2fs fast3=%.3fs fast2=%.3fs ratio=%.3f"
          % (t_jit, t3, t2, t2 / t3))
    print("timing:", json.dumps(tm, indent=1))
    print("bit_identical=%s max_abs_diff=%.3g info_equal=%s"
          % (np.array_equal(C3, C2), float(d.max()), info2 == info3))
