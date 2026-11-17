"""transport_mesh_fast -- BRIEF-V4-1-16-FASTDEPO: the SAME area-conserving
deposition as `transport_mesh.py`, only faster.

WHAT IS AND IS NOT ALLOWED HERE
-------------------------------
The brief allows only *neutral numerics*: the same area-preserving mesh
projection, computed faster.  Concretely, this module changes NOTHING about

  * the projection (`transport_mesh.landing_pixels` and
    `transport_mesh.cell_triangles` are imported and called verbatim -- the
    Snell refraction, the corner grid, the quad -> 2 triangle split, the
    per-cell energy w_cell = W*H/N_s^2 are untouched),
  * the deposit *mathematics*: every triangle's energy is still distributed
    in exact proportion to the analytic overlap area of the projected
    triangle with each pixel square, obtained by the same four half-plane
    Sutherland-Hodgman clips (x >= x0, x <= x0+1, y >= y0, y <= y0+1) and the
    same shoelace area.  No splat, no approximation, no clipping shortcut, no
    coverage approximation, no solar-disc change, no renormalisation.
  * the periodic bookkeeping (i mod W, j mod H), the bulk/tail split rule,
    the uniform backstop for undeposited tail area, and every reported
    diagnostic of `deposit`.

What changes is only HOW the same arithmetic is executed: the Python/numpy
element-wise clip stack of `_clip`/`tri_box_area` is replaced, box by box, by
a Numba (LLVM) loop that performs the identical operations.  Because the
individual sums are evaluated in a different order (grouping of triangles is
a vectorisation artefact and is not reproduced; per-pixel accumulation is
sequential in both), floating-point rounding can differ in the last bits.
That is what the pre-registered equivalence gates of the brief bound.

The public API mirrors the verbatim module:

    deposit_fast(T, A_tri, A_cell, w_cell, W, H)      -> acc (H,W), info
    mesh_C_fast(fields, S, W, H, env, ior=IOR_G)      -> C (H,W) float32, info
    tri_box_area_fast(T, x0, y0)                      -> exact areas (px^2)

and `BOX_BULK_MAX` is exposed as a module constant with the same meaning as
in `transport_mesh` (which triangles take the vectorised bulk path).  Both
paths compute exact areas; the constant is a performance switch only.

SELF TEST
---------
    python transport_mesh_fast.py --selftest

compares the fast kernels against `transport_mesh` on: the verbatim module's
own exactness cases (flat water, lattice shear, integer dilation, lattice
rotation, non-lattice rotation), random triangle/box overlaps, the tail path
forced on, and a random wave field through the whole mesh_C chain.
"""
import numpy as np
from numba import njit

import transport_mesh as TM

IOR_G = TM.IOR_G
K_CAP = TM.K_CAP                       # clipper vertex capacity (8)
UNIFORM_GUARD = TM.UNIFORM_GUARD

#: same documented rasterisation constant as the verbatim module: triangles
#: whose integer bbox holds <= BOX_BULK_MAX pixels go to the bulk path.  Both
#: paths are exact; this is a speed switch, nothing else.
BOX_BULK_MAX = TM.BOX_BULK_MAX         # 16

#: rows of candidate boxes per block in the verbatim bulk path.  The fast
#: bulk path is a per-triangle loop and does not need the block, but the name
#: is kept so callers can see the parallel.
CHUNK_ROWS = TM.CHUNK_ROWS


# ---------------------------------------------------------------------------
# the exact clip, one polygon at a time (identical arithmetic to TM._clip)
# ---------------------------------------------------------------------------
@njit(cache=True, fastmath=False)
def _clip_xy(sx, sy, n, a, b, c, dx, dy):
    """Clip the convex polygon (sx[0:n], sy[0:n]) by a*x+b*y-c >= 0.

    Vertex emission order is the verbatim module's: for every input slot, the
    intersection point with the previous edge comes first, then the vertex
    itself, packed in slot order; capacity K_CAP.  Written straight from
    `transport_mesh._clip` (same f, same edge from the previous slot, same
    1e-300 denominator guard, same interpolation f/(f-fp)).
    """
    m = 0
    for i in range(n):
        f = a * sx[i] + b * sy[i] - c
        ins = f >= 0.0
        j = i - 1
        if j < 0:
            j = n - 1
        fp = a * sx[j] + b * sy[j] - c
        if ins != (fp >= 0.0):
            den = f - fp
            if abs(den) < 1e-300:
                den = 1e-300
            t = f / den
            if m < K_CAP:
                dx[m] = sx[i] + (sx[j] - sx[i]) * t
                dy[m] = sy[i] + (sy[j] - sy[i]) * t
            m += 1
        if ins:
            if m < K_CAP:
                dx[m] = sx[i]
                dy[m] = sy[i]
            m += 1
    if m > K_CAP:
        m = K_CAP
    return m


@njit(cache=True, fastmath=False)
def _shoelace_xy(x, y, n):
    s = 0.0
    for i in range(n):
        j = i + 1
        if j >= n:
            j = 0
        s += x[i] * y[j] - x[j] * y[i]
    return 0.5 * abs(s)


@njit(cache=True, fastmath=False)
def _tri_box_area_one(T, x0, y0, ax, ay, bx, by):
    """Exact area of (triangle T ∩ [x0,x0+1]x[y0,y0+1]), same 4 clips."""
    ax[0] = T[0, 0]
    ay[0] = T[0, 1]
    ax[1] = T[1, 0]
    ay[1] = T[1, 1]
    ax[2] = T[2, 0]
    ay[2] = T[2, 1]
    n = _clip_xy(ax, ay, 3, 1.0, 0.0, x0, bx, by)
    n = _clip_xy(bx, by, n, -1.0, 0.0, -(x0 + 1.0), ax, ay)
    n = _clip_xy(ax, ay, n, 0.0, 1.0, y0, bx, by)
    n = _clip_xy(bx, by, n, 0.0, -1.0, -(y0 + 1.0), ax, ay)
    return _shoelace_xy(ax, ay, n)


@njit(cache=True, fastmath=False)
def _tri_box_area_many(T, x0, y0):
    """Vectorised entry point: exact box areas for many triangles at once."""
    M = T.shape[0]
    out = np.zeros(M, dtype=np.float64)
    ax = np.zeros(K_CAP, dtype=np.float64)
    ay = np.zeros(K_CAP, dtype=np.float64)
    bx = np.zeros(K_CAP, dtype=np.float64)
    by = np.zeros(K_CAP, dtype=np.float64)
    for i in range(M):
        out[i] = _tri_box_area_one(T[i], x0[i], y0[i], ax, ay, bx, by)
    return out


def tri_box_area_fast(T, x0, y0):
    """Exact area of (triangle ∩ unit pixel) for every triangle (same API)."""
    T = np.ascontiguousarray(np.asarray(T, dtype=np.float64))
    x0 = np.ascontiguousarray(np.asarray(x0, dtype=np.float64))
    y0 = np.ascontiguousarray(np.asarray(y0, dtype=np.float64))
    return _tri_box_area_many(T, x0, y0)


# ---------------------------------------------------------------------------
# deposition kernels
# ---------------------------------------------------------------------------
@njit(cache=True, fastmath=False)
def _deposit_bulk_njit(T, i0, j0, sw, sh, dens, acc, W, H):
    """Exact overlap on the integer bbox stencil of each triangle.

    Returns (area_in, area_out, energy_in, energy_out, bulk_boxes).  The
    arithmetic per box is `_tri_box_area_one`; the accumulation is
    `acc[(y mod H)*W + (x mod W)] += dens*a`, the same map `bincount` builds
    in the verbatim path.  Triangles are visited in index order and boxes in
    (row, column) order, so the per-pixel addition order is the verbatim one.
    """
    M = T.shape[0]
    area_in = 0.0
    area_out = 0.0
    energy_in = 0.0
    energy_out = 0.0
    boxes = 0
    ax = np.zeros(K_CAP, dtype=np.float64)
    ay = np.zeros(K_CAP, dtype=np.float64)
    bx = np.zeros(K_CAP, dtype=np.float64)
    by = np.zeros(K_CAP, dtype=np.float64)
    Wf = float(W)
    Hf = float(H)
    for t in range(M):
        swt = int(sw[t])
        sht = int(sh[t])
        dt = dens[t]
        for gy in range(sht):
            yi = j0[t] + gy
            iy = int(yi) % H
            inside_y = (yi >= 0.0) and (yi < Hf)
            for gx in range(swt):
                xi = i0[t] + gx
                a = _tri_box_area_one(T[t], xi, yi, ax, ay, bx, by)
                w = dt * a
                ix = int(xi) % W
                acc[iy * W + ix] += w
                if inside_y and (xi >= 0.0) and (xi < Wf):
                    area_in += a
                    energy_in += w
                else:
                    area_out += a
                    energy_out += w
                boxes += 1
    return area_in, area_out, energy_in, energy_out, boxes


@njit(cache=True, fastmath=False)
def _deposit_large_njit(P, dens_one, acc, W, H):
    """Tail path: huge / sliver triangle, rasterised row strip by row strip.

    Same as `transport_mesh._deposit_large`: clip to the row band, take the
    integer column span of the surviving polygon, clamp it to [-2W, 3W] and to
    4W per row, then clip every candidate column exactly.

    Returns (deposited_area, area_in, area_out, energy_in, energy_out,
             tail_boxes, did_work).
    """
    Wf = float(W)
    Hf = float(H)
    ymin = P[0, 1]
    ymax = P[0, 1]
    for k in range(1, 3):
        if P[k, 1] < ymin:
            ymin = P[k, 1]
        if P[k, 1] > ymax:
            ymax = P[k, 1]
    lo = ymin
    if -2.0 * Hf > lo:
        lo = -2.0 * Hf
    hi = ymax
    if 3.0 * Hf < hi:
        hi = 3.0 * Hf
    j0 = int(np.floor(lo))
    j1 = int(np.floor(hi))
    if j1 < j0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0
    area_in = 0.0
    area_out = 0.0
    energy_in = 0.0
    energy_out = 0.0
    boxes = 0
    total = 0.0
    ax = np.zeros(K_CAP, dtype=np.float64)
    ay = np.zeros(K_CAP, dtype=np.float64)
    bx = np.zeros(K_CAP, dtype=np.float64)
    by = np.zeros(K_CAP, dtype=np.float64)
    rx = np.zeros(K_CAP, dtype=np.float64)
    ry = np.zeros(K_CAP, dtype=np.float64)
    for row in range(j0, j1 + 1):
        rf = float(row)
        # clip the triangle to the row band [row, row+1]; the surviving row
        # polygon is kept in (rx, ry) because every column clip below must
        # start from it
        ax[0] = P[0, 0]
        ay[0] = P[0, 1]
        ax[1] = P[1, 0]
        ay[1] = P[1, 1]
        ax[2] = P[2, 0]
        ay[2] = P[2, 1]
        n = _clip_xy(ax, ay, 3, 0.0, 1.0, rf, bx, by)
        n = _clip_xy(bx, by, n, 0.0, -1.0, -(rf + 1.0), rx, ry)
        if n <= 0:
            continue
        xlo = rx[0]
        xhi = rx[0]
        for k in range(1, n):
            if rx[k] < xlo:
                xlo = rx[k]
            if rx[k] > xhi:
                xhi = rx[k]
        c0 = np.floor(xlo)
        c1 = np.floor(xhi)
        if c0 < -2.0 * Wf:
            c0 = -2.0 * Wf
        if c1 > 3.0 * Wf:
            c1 = 3.0 * Wf
        if (c1 - c0) > 4.0 * Wf:
            c1 = c0 + 4.0 * Wf
        cnt = int(c1 - c0 + 1.0)
        if cnt < 0:
            cnt = 0
        for k in range(cnt):
            col = c0 + float(k)
            m = _clip_xy(rx, ry, n, 1.0, 0.0, col, bx, by)
            m = _clip_xy(bx, by, m, -1.0, 0.0, -(col + 1.0), ax, ay)
            a = _shoelace_xy(ax, ay, m)
            w = a * dens_one
            iy = row % H
            ix = int(col) % W
            acc[iy * W + ix] += w
            if (col >= 0.0) and (col < Wf) and (rf >= 0.0) and (rf < Hf):
                area_in += a
                energy_in += w
            else:
                area_out += a
                energy_out += w
            total += a
            boxes += 1
    return total, area_in, area_out, energy_in, energy_out, boxes, 1


# ---------------------------------------------------------------------------
# deposit (same bookkeeping as the verbatim module)
# ---------------------------------------------------------------------------
def deposit_fast(T, A_tri, A_cell, w_cell, W, H):
    """Distribute every cell's energy by exact overlap area (same as deposit).

    Returns (acc (H,W) float64, info dict) with the verbatim keys.
    """
    W, H = int(W), int(H)
    acc = np.zeros(W * H, dtype=np.float64)
    T = np.ascontiguousarray(np.asarray(T, dtype=np.float64))
    A_tri = np.asarray(A_tri, dtype=np.float64)
    A_cell = np.asarray(A_cell, dtype=np.float64)
    w_cell = np.asarray(w_cell, dtype=np.float64)
    info = {"area_in": 0.0, "area_out": 0.0, "energy_in": 0.0, "energy_out": 0.0,
            "bulk_boxes": 0, "tail_boxes": 0, "n_bulk": 0, "n_tail": 0,
            "energy_guard": 0.0, "n_degenerate_cell": 0}
    ncell = A_cell.size
    dens_cell = np.zeros(ncell)
    ok = A_cell > UNIFORM_GUARD
    dens_cell[ok] = w_cell[ok] / A_cell[ok]
    info["n_degenerate_cell"] = int((~ok).sum())
    dens = np.concatenate([dens_cell, dens_cell])          # per triangle
    wtri = np.concatenate([w_cell * (A_tri[:ncell] / np.maximum(A_cell, UNIFORM_GUARD)),
                           w_cell * (A_tri[ncell:] / np.maximum(A_cell, UNIFORM_GUARD))])

    xmin = T[:, :, 0].min(1)
    xmax = T[:, :, 0].max(1)
    ymin = T[:, :, 1].min(1)
    ymax = T[:, :, 1].max(1)
    i0 = np.floor(xmin)
    j0 = np.floor(ymin)
    sw = np.floor(xmax) - i0 + 1
    sh = np.floor(ymax) - j0 + 1
    per = sw * sh
    bulk = np.nonzero(per <= BOX_BULK_MAX)[0]
    tail = np.nonzero(per > BOX_BULK_MAX)[0]
    info["n_bulk"] = int(bulk.size)
    info["n_tail"] = int(tail.size)

    if bulk.size:
        bi0 = np.ascontiguousarray(i0[bulk])
        bj0 = np.ascontiguousarray(j0[bulk])
        bsw = np.ascontiguousarray(sw[bulk])
        bsh = np.ascontiguousarray(sh[bulk])
        bdens = np.ascontiguousarray(dens[bulk])
        a_in, a_out, e_in, e_out, nbox = _deposit_bulk_njit(
            T[bulk], bi0, bj0, bsw, bsh, bdens, acc, W, H)
        info["area_in"] += float(a_in)
        info["area_out"] += float(a_out)
        info["energy_in"] += float(e_in)
        info["energy_out"] += float(e_out)
        info["bulk_boxes"] += int(nbox)

    # ---- tail: one triangle at a time (rare by construction) --------------
    tail_area = np.zeros(tail.size)
    if tail.size:
        for k, t in enumerate(tail):
            (tot, a_in, a_out, e_in, e_out, nbox,
             did) = _deposit_large_njit(np.ascontiguousarray(T[t]),
                                        float(dens[t]), acc, W, H)
            tail_area[k] = tot
            info["area_in"] += float(a_in)
            info["area_out"] += float(a_out)
            info["energy_in"] += float(e_in)
            info["energy_out"] += float(e_out)
            info["tail_boxes"] += int(nbox)

    # ---- backstop: undeposited tail area is spread uniformly and REPORTED -
    deficit_px = 0.0
    for k, t in enumerate(tail):
        if tail_area[k] < A_tri[t] * (1.0 - 1e-9):
            deficit_px += float(A_tri[t] - tail_area[k]) * float(dens[t])
    info["energy_guard"] = float(deficit_px)
    total_w = float(np.sum(wtri))
    info["energy_cell_total"] = total_w
    info["energy_deposited"] = float(acc.sum())
    if deficit_px > 0.0:
        acc += deficit_px / float(W * H)
        info["energy_deposited_after_guard"] = float(acc.sum())
    info["rel_err"] = float(abs(info.get("energy_deposited_after_guard",
                                          info["energy_deposited"]) - total_w)
                             / total_w)
    info["offframe_fraction_area"] = (info["area_out"]
                                      / max(info["area_in"] + info["area_out"], 1e-30))
    info["offframe_fraction_energy"] = (info["energy_out"]
                                        / max(info["energy_in"] + info["energy_out"], 1e-30))
    return acc.reshape(H, W), info


def mesh_C_fast(fields, S, W, H, env, ior=IOR_G, w_total=None):
    """Area-conserving C for one time slice, via the fast kernels.

    Identical chain to `transport_mesh.mesh_C`: the projection and the
    triangle decomposition are the verbatim module's own functions; only the
    deposition arithmetic is executed by the Numba kernels.
    """
    N = np.asarray(fields["h"]).shape[0]
    PXe, PYe = TM.landing_pixels(fields, S, env, ior=ior, W=W, H=H)
    T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
    wt = float(W) * float(H) / (float(N) * float(N)) if w_total is None else float(w_total)
    w_cell = np.full(A_cell.size, wt)
    acc, info = deposit_fast(T, A_tri, A_cell, w_cell, W, H)
    info.update({"N_s": int(N), "W": int(W), "H": int(H), "S_m": float(S),
                 "w_cell": wt, "ior": float(ior),
                 "method": "area push-forward (quad -> 2 triangles, exact "
                           "overlap with the frame pixels); numba kernel, "
                           "same arithmetic as transport_mesh.deposit",
                 "kernel": "none", "jitter": "none", "band_limit": "none",
                 "per_cell_energy_rule": "w_cell = W*H/N_s^2 (det.det_frame's "
                                         "divisor convention, rectangular "
                                         "frame analogue)"})
    if not np.all(np.isfinite(acc)):
        raise FloatingPointError("non-finite C")
    return acc.astype(np.float32), info


# ---------------------------------------------------------------------------
# self test
# ---------------------------------------------------------------------------
def _selftest():
    import time
    ok = True

    def chk(name, cond, detail=""):
        nonlocal ok
        ok = ok and bool(cond)
        print("[%s] %s %s" % ("ok" if cond else "NG", name, detail))

    rng = np.random.default_rng(20261116)

    # ---- 1. per-box areas, random triangles vs verbatim tri_box_area ------
    M = 4000
    c = rng.uniform(2, 600, size=(M, 2))
    th = rng.uniform(0, 2 * np.pi, size=(M, 3, 1))
    rr = rng.uniform(0.02, 60.0, size=(M, 3, 1))
    Tr = c[:, None, :] + rr * np.stack([np.cos(th[..., 0]), np.sin(th[..., 0])], axis=-1)
    x0 = np.floor(Tr[:, :, 0].min(1))
    y0 = np.floor(Tr[:, :, 1].min(1))
    dmax = 0.0
    drel = 0.0
    for dx in range(6):
        for dy in range(6):
            A_ref = TM.tri_box_area(Tr, x0 + dx, y0 + dy)
            A_fast = tri_box_area_fast(Tr, x0 + dx, y0 + dy)
            d = np.abs(A_ref - A_fast)
            dmax = max(dmax, float(d.max()))
            drel = max(drel, float((d / np.maximum(np.abs(A_ref), 1e-30)).max()))
    chk("fast per-box area == verbatim tri_box_area (4000 triangles x 36 boxes)",
        dmax < 1e-9, "max|diff| = %.3e px^2, max rel = %.2e" % (dmax, drel))

    # ---- 2. the verbatim module's own exactness cases ---------------------
    N, W, H = 64, 736, 540
    S = 1.4755102040816326
    z = np.zeros((N, N), dtype=np.float32)
    f = {"h": z, "gx": z, "gy": z, "Dx": z, "Dy": z, "S": S}
    env = {"depthM": 0.6, "sunDir": np.array([0.4965974, -0.1350766, 0.8574062])}
    Cv, iv = TM.mesh_C(f, S, W, H, env)
    Cf, iff = mesh_C_fast(f, S, W, H, env)
    chk("flat water: C == 1", np.max(np.abs(Cf - 1.0)) < 1e-9,
        "max|C-1| = %.3e (verbatim %.3e)"
        % (np.max(np.abs(Cf - 1.0)), np.max(np.abs(Cv - 1.0))))
    chk("flat water: fast C == verbatim C", np.max(np.abs(Cf - Cv)) < 1e-12,
        "max|diff| = %.3e" % np.max(np.abs(Cf - Cv)))

    axc = np.arange(N + 1) / N
    XC, YC = np.meshgrid(axc, axc, indexing="xy")
    accv = {}          # name -> (PXe, PYe, W, H, w_cell value)
    wt = W * H / float(N * N)
    accv["lattice shear"] = ((XC + 0.5 * YC) * W, YC * H, W, H, wt)
    accv["integer dilation J=2"] = ((0.5 + (XC - 0.5) * 2.0) * W,
                                    (0.5 + (YC - 0.5) * 2.0) * H, W, H, wt)
    Nsq, Wsq = 64, 368
    aqs = np.arange(Nsq + 1) / Nsq
    XS, YS = np.meshgrid(aqs, aqs, indexing="xy")
    wtsq = Wsq * Wsq / float(Nsq * Nsq)
    accv["90 deg rotation + shift"] = ((1.0 - YS + 0.37) * Wsq,
                                       (XS + 0.23) * Wsq, Wsq, Wsq, wtsq)
    th2 = 0.41
    XC0, YC0 = XS - 0.5, YS - 0.5
    accv["non-lattice rotation"] = (
        (0.5 + np.cos(th2) * XC0 - np.sin(th2) * YC0) * Wsq,
        (0.5 + np.sin(th2) * XC0 + np.cos(th2) * YC0) * Wsq, Wsq, Wsq, wtsq)
    for name, (PXe, PYe, w, hh, wval) in accv.items():
        T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
        wc = np.full(A_cell.size, wval)
        acc_v, info_v = TM.deposit(T, A_tri, A_cell, wc, w, hh)
        acc_f, info_f = deposit_fast(T, A_tri, A_cell, wc, w, hh)
        d = float(np.max(np.abs(acc_v - acc_f)))
        chk("%s: fast == verbatim" % name, d < 1e-9 and info_f["rel_err"] < 1e-11,
            "max|diff| = %.3e  rel_err %.2e (verbatim %.2e)  fold %.2f"
            % (d, info_f["rel_err"], info_v["rel_err"],
               info_f["offframe_fraction_energy"]))

    # ---- 3. the tail path, forced on by BOX_BULK_MAX = 0 ------------------
    PXe, PYe, w, hh, wval = accv["non-lattice rotation"]
    T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
    wc = np.full(A_cell.size, wval)
    global BOX_BULK_MAX
    old = BOX_BULK_MAX
    BOX_BULK_MAX = 0
    acc_v, info_v = TM.deposit(T, A_tri, A_cell, wc, w, hh)
    acc_f, info_f = deposit_fast(T, A_tri, A_cell, wc, w, hh)
    BOX_BULK_MAX = old
    d = float(np.max(np.abs(acc_v - acc_f)))
    chk("tail path only: fast == verbatim",
        d < 1e-9 and abs(info_f["rel_err"] - info_v["rel_err"]) < 1e-12,
        "max|diff| = %.3e  rel_err fast %.2e / verbatim %.2e  tail boxes %d/%d"
        % (d, info_f["rel_err"], info_v["rel_err"],
           info_f["tail_boxes"], info_v["tail_boxes"]))

    # ---- 4. a real wave field through mesh_C (128^2 -> 736x540) ----------
    Nf = 128
    af = np.arange(Nf) / float(Nf)
    Xf, Yf = np.meshgrid(af, af, indexing="xy")
    hh = (0.01 * np.sin(2 * np.pi * 5 * Xf) * np.cos(2 * np.pi * 3 * Yf)
          + 0.006 * np.sin(2 * np.pi * 11 * Xf + 0.7)).astype(np.float32)
    gxx = np.gradient(np.asarray(hh, dtype=np.float64), axis=1) * Nf / S
    gyy = np.gradient(np.asarray(hh, dtype=np.float64), axis=0) * Nf / S
    fr = {"h": hh, "gx": gxx.astype(np.float32), "gy": gyy.astype(np.float32),
          "Dx": np.zeros_like(hh), "Dy": np.zeros_like(hh), "S": S}
    envd = {"depthM": 4.0, "sunDir": np.array([0.4965974, -0.1350766, 0.8574062])}
    Cv, iv = TM.mesh_C(fr, S, W, H, envd)
    Cf, iff = mesh_C_fast(fr, S, W, H, envd)
    d = float(np.max(np.abs(Cv.astype(np.float64) - Cf.astype(np.float64))))
    chk("depth 4.0 m wave field: fast == verbatim",
        d < 1e-6,
        "max|diff| = %.3e  std %.6f/%.6f  p99 %.6f/%.6f  rel_err %.2e/%.2e"
        % (d, float(np.std(Cf)), float(np.std(Cv)),
           float(np.percentile(Cf, 99)), float(np.percentile(Cv, 99)),
           iff["rel_err"], iv["rel_err"]))

    # ---- 5. timing, one 512^2 field --------------------------------------
    t0 = time.perf_counter(); TM.mesh_C(fr, S, W, H, envd); tv = time.perf_counter() - t0
    t0 = time.perf_counter(); mesh_C_fast(fr, S, W, H, envd); tf = time.perf_counter() - t0
    print("[--] 128^2 field -> 736x540: verbatim %.2f s, fast %.2f s (%.1fx)"
          % (tv, tf, tv / max(tf, 1e-9)))
    print("-" * 60)
    print("SELFTEST_OK" if ok else "SELFTEST_FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    print(__doc__)
