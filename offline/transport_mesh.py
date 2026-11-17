"""transport_mesh -- BRIEF-V4-1-8: area-conserving surface -> floor deposition.

WHAT THIS IS
------------
The v4-1.7 lane (airlines: point splat of a regular ray lattice, optionally
jittered, optionally convolved with a splatting kernel) is replaced here by a
*measure-preserving push-forward* of the water surface onto the floor:

    every surface cell (a quad of 4 corners -> 2 triangles) is projected
    through the SAME Snell refraction as the existing transport (det.det_frame:
    unit solar direction, eta = 1/1.33434, plane-parallel depth z), and the
    cell's energy is distributed over the frame pixels in EXACT proportion to
    the area of overlap between the projected triangle and the pixel.

No kernel, no jitter, no band limit, no splatting: the overlap areas are
computed analytically (Sutherland-Hodgman clipping of the triangle against the
pixel square, exact up to floating point).

ENERGY DEFINITION (reused verbatim from the existing C lane)
------------------------------------------------------------
det.det_frame gives every surface texel Q unit-weight rays and divides the
accumulated buffer by

    divisor = Q * N_s^2 / N_buf^2          (N_buf = the C buffer side)

so the C-mass carried by ONE surface cell is

    w_cell = Q / divisor = N_buf^2 / N_s^2   (dimensionless).

For a rectangular frame buffer (W x H pixels) the same convention reads

    w_cell = W * H / N_s^2

i.e. every cell carries the same mass, because every cell has the same area.
The dimensionless ratio then has mean exactly 1 for flat water and NO measured
renormalisation is applied anywhere (`mean 1` is a definition, not a fit).  The
deposited total is recorded per frame together with the cell-energy total.

MAPPING (v2 framing, no magnification)
--------------------------------------
The surface patch IS the screen window: patch side L (metres) is declared to
span the frame width W, exactly as run_v4_1_2/v16 do (`px = u*S, py = v*S`), so
the magnification factor between the generated patch and the drawn patch is 1
(in v4-1.6 that factor was 6.4 and the wave was stretched).  One frame pixel
therefore covers  dx_px = L/W  metres,  dy_px = L/H  metres  -- the residual
pixel aspect 736:540 is v2's framing and is NOT a magnification: the photo, the
slope field and C all live in the same (u,v) parameterisation.

PERIODICITY
-----------
The floor is periodic on the patch (the composite samples C with a wrapped
(u,v) lookup), so a projected triangle that leaves the frame re-enters on the
opposite side: the deposition is accumulated at (i mod W, j mod H).  The
fraction of the *unwrapped* landing area that falls outside the frame is
computed and reported separately (an energy bookkeeping diagnostic, never used
to rescale anything).

PUBLIC API
----------
    landing_pixels(fields, S, env, ior=IOR_G)   -> corner-grid pixel coords
    cell_triangles(PXe, PYe)                    -> triangles + areas
    deposit(T, dens, W, H, ...)                 -> acc (H,W), info
    mesh_C(fields, S, W, H, env, ior=IOR_G)     -> C (H,W) float32, info

REFERENCE LANE (not part of the deliverable)
-------------------------------------------
    point_splat_C(fields, S, W, H, env, ior, Q=25)
is `det.det_frame`'s existing point-splat transport (texel centres, Fibonacci
solar-disc quadrature, bilinear splat, analytical divisor, sigma = 0.75 px
periodic blur) evaluated on the SAME W x H frame grid, so that the peak
statistics of the area deposit can be compared with the legacy lane on one
grid.  It is never used for C of the deliverable.  On a square grid it
reproduces det.det_frame bit-for-bit (asserted by the self test).

SELF TEST
---------
    python transport_mesh.py --selftest

checks (i) a flat surface gives C == 1 everywhere (the real refraction is a
translation on the torus), (ii) maps that ARE torus automorphisms give exactly
the flat per-pixel value to 1e-12 (a shear, a 90-degree rotation + shift that
wraps, and an integer dilation J=2 -- a degree-4 cover, where each of the four
branches carries 1/|det J| = 1/4 and the pixel sums back to the flat value),
(iii) the energy identity, (iv) agreement with a dense supersampled reference on
random small triangles, (v) that the reference point splat is `det.det_frame`
bit-for-bit on a square grid and conserves mass exactly on flat water, and
(vi) that a rotation by a NON-lattice angle is not a bijection of the torus:
energy is still exact, C is not 1 -- the control that pins down why an
arbitrary "rigid rotation" cannot be used as an exactness test.

Note on (ii): the maps are evaluated on a (N+1)x(N+1) CORNER grid that spans
the fundamental domain exactly (u = i/N, i = 0..N).  The previous version of
this test sampled cell CENTRES and then manufactured the seam column/row by
adding W / H, which is inconsistent with any map that is not a pure
translation -- that is what the old `rotated rigid map' NG (1.466) was
measuring, not a deposit error.  See REPORT-V4-1-8 section 3.
"""
import numpy as np

IOR_G = 1.33434

# rasterisation constants (documented, not tuned)
K_CAP = 8            # vertex capacity of the clipper (triangle + 4 box planes -> <=7)
BOX_BULK_MAX = 16    # bulk path: triangles whose integer bbox holds <= 16 pixels
CHUNK_ROWS = 250000  # candidate rows per numpy block in the bulk path
UNIFORM_GUARD = 1e-12


# ---------------------------------------------------------------------------
# geometry: the existing refraction, evaluated on a cell-corner grid
# ---------------------------------------------------------------------------
def _corner(a):
    """Value at cell corners = bilinear mean of the 4 surrounding texels."""
    return 0.25 * (a + np.roll(a, 1, 0) + np.roll(a, 1, 1)
                   + np.roll(np.roll(a, 1, 0), 1, 1))


def landing_pixels(fields, S, env, ior=IOR_G, W=None, H=None):
    """Project the corner grid of the surface onto the floor, in frame pixels.

    Returns (PXe, PYe): (N_s+1, N_s+1) arrays of unwrapped pixel coordinates
    (x rightwards, y downwards), extended by one row/column so that cell
    (i,j) has corners (i,j), (i,j+1), (i+1,j), (i+1,j+1).  The extension is
    exact for a periodic field: landing(x+S) = landing(x)+S.
    """
    h = np.asarray(fields["h"], dtype=np.float64)
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    Dx = np.asarray(fields.get("Dx", np.zeros_like(h)), dtype=np.float64)
    Dy = np.asarray(fields.get("Dy", np.zeros_like(h)), dtype=np.float64)
    N = h.shape[0]
    S = float(S)
    W = float(W if W is not None else N)
    H = float(H if H is not None else N)
    depth = float(env["depthM"])

    s = np.asarray(env["sunDir"], dtype=np.float64)
    s = s / np.linalg.norm(s)
    I = -s
    eta = 1.0 / float(ior)

    dxm = S / N
    ax = np.arange(N, dtype=np.float64) * dxm
    X, Y = np.meshgrid(ax, ax, indexing="xy")
    Hc, GXc, GYc = _corner(h), _corner(gx), _corner(gy)

    nl = np.sqrt(GXc * GXc + GYc * GYc + 1.0)
    Nx, Ny, Nz = -GXc / nl, -GYc / nl, 1.0 / nl
    cos1 = np.clip(-(I[0] * Nx + I[1] * Ny + I[2] * Nz), -1.0, 1.0)
    kk = np.maximum(1.0 - eta * eta * (1.0 - cos1 * cos1), 0.0)
    sq = np.sqrt(kk)
    Tx = eta * I[0] + (eta * cos1 - sq) * Nx
    Ty = eta * I[1] + (eta * cos1 - sq) * Ny
    Tz = np.minimum(eta * I[2] + (eta * cos1 - sq) * Nz, -1e-6)
    t = (Hc + depth) / (-Tz)

    xb = X + np.asarray(Dx, dtype=np.float64) * 0.0 + Tx * t
    yb = Y + np.asarray(Dy, dtype=np.float64) * 0.0 + Ty * t
    PX = xb / S * W
    PY = yb / S * H

    # periodic extension: corner at column N is corner 0 shifted by S (= W px),
    # corner at row N is corner 0 shifted by S (= H px in pixel space).
    PXe = np.concatenate([PX, PX[:, :1] + W], axis=1)
    PYe = np.concatenate([PY, PY[:, :1]], axis=1)
    PXe = np.concatenate([PXe, PXe[:1, :]], axis=0)
    PYe = np.concatenate([PYe, PYe[:1, :] + H], axis=0)
    return PXe, PYe


def cell_triangles(PXe, PYe):
    """Quad -> 2 triangles per cell (the brief's decomposition).

    Returns T (2*N*N, 3, 2) with the first N*N rows the (00,10,11) triangle and
    the last N*N rows the (00,11,01) triangle, plus A_tri (2N^2,) unsigned
    areas and A_cell (N^2,) = |A1| + |A2|.
    """
    p00 = np.stack([PXe[:-1, :-1], PYe[:-1, :-1]], axis=-1)
    p10 = np.stack([PXe[:-1, 1:], PYe[:-1, 1:]], axis=-1)
    p01 = np.stack([PXe[1:, :-1], PYe[1:, :-1]], axis=-1)
    p11 = np.stack([PXe[1:, 1:], PYe[1:, 1:]], axis=-1)
    T1 = np.stack([p00, p10, p11], axis=-2).reshape(-1, 3, 2)
    T2 = np.stack([p00, p11, p01], axis=-2).reshape(-1, 3, 2)

    def area(T):
        return 0.5 * np.abs((T[:, 1, 0] - T[:, 0, 0]) * (T[:, 2, 1] - T[:, 0, 1])
                            - (T[:, 2, 0] - T[:, 0, 0]) * (T[:, 1, 1] - T[:, 0, 1]))
    A1, A2 = area(T1), area(T2)
    return (np.concatenate([T1, T2], axis=0), np.concatenate([A1, A2], axis=0),
            A1 + A2)


# ---------------------------------------------------------------------------
# exact triangle / box overlap (analytic clipping)
# ---------------------------------------------------------------------------
def _clip(V, cnt, a, b, c):
    """Clip convex polygons V (M,K,2) by the half-plane a*x+b*y-c >= 0.

    ``c`` may be a scalar or an (M,) array (one plane per polygon).  Returns
    (V2, cnt2) with the surviving vertices packed to the front, capacity K.
    """
    M, K = V.shape[0], V.shape[1]
    slots = np.arange(K, dtype=np.int64)[None, :]
    valid = slots < cnt[:, None]
    c = np.asarray(c, dtype=np.float64)
    f = a * V[:, :, 0] + b * V[:, :, 1] - (c[:, None] if c.ndim else c)
    f = np.where(valid, f, 0.0)
    ins = f >= 0.0
    last = np.maximum(cnt - 1, 0)
    ip = np.where(slots == 0, last[:, None], slots - 1)
    ip = np.where(valid, ip, 0)
    rows = np.arange(M)
    fp = f[rows[:, None], ip]
    Vp = V[rows[:, None], ip]
    cross = (ins != (fp >= 0.0)) & valid
    den = f - fp
    den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    X = V + (Vp - V) * (f / den)[:, :, None]

    out = np.zeros((M, 2 * K, 2))
    m0 = np.zeros((M, 2 * K), dtype=bool)
    m1 = np.zeros((M, 2 * K), dtype=bool)
    m0[:, 0::2] = cross          # intersection point comes first (SH order)
    m1[:, 1::2] = ins & valid    # then the vertex itself
    out[:, 0::2, :] = X
    out[:, 1::2, :] = V
    mask = m0 | m1
    order = np.cumsum(mask, axis=1) - 1
    r, cc = np.nonzero(mask)
    sel = order[r, cc]
    kout = min(K + 1, K_CAP)      # one half-plane adds at most one vertex
    keep = sel < kout
    V2 = np.zeros((M, kout, 2))
    V2[r[keep], sel[keep]] = out[r[keep], cc[keep]]
    return V2, np.minimum(mask.sum(1), kout)


def _shoelace(V, cnt):
    M, K = V.shape[0], V.shape[1]
    slots = np.arange(K, dtype=np.int64)[None, :]
    j = np.where(slots + 1 < cnt[:, None], slots + 1, 0)
    rows = np.arange(M)[:, None]
    cr = V[:, :, 0] * V[rows, j, 1] - V[rows, j, 0] * V[:, :, 1]
    valid = slots < cnt[:, None]
    return 0.5 * np.abs(np.where(valid, cr, 0.0).sum(1))


def tri_box_area(T, x0, y0):
    """Exact area of (triangle ∩ [x0,x0+1]x[y0,y0+1]) for every triangle."""
    M = T.shape[0]
    V = np.zeros((M, K_CAP, 2))
    V[:, :3, :] = T
    cnt = np.full(M, 3, dtype=np.int64)
    V, cnt = _clip(V, cnt, 1.0, 0.0, x0)
    V, cnt = _clip(V, cnt, -1.0, 0.0, -(np.asarray(x0) + 1.0))
    V, cnt = _clip(V, cnt, 0.0, 1.0, y0)
    V, cnt = _clip(V, cnt, 0.0, -1.0, -(np.asarray(y0) + 1.0))
    return _shoelace(V, cnt)


# ---------------------------------------------------------------------------
# deposition
# ---------------------------------------------------------------------------
def _deposit_bulk(T, i0, j0, sw, sh, dens, acc, W, H, info):
    """Bulk path: exact overlap on the integer bbox stencil of each triangle."""
    key = sw.astype(np.int64) * 4096 + sh.astype(np.int64)
    for kv in np.unique(key):
        idx = np.nonzero(key == kv)[0]
        swg = int(sw[idx[0]])
        shg = int(sh[idx[0]])
        dxs = np.arange(swg)
        dys = np.arange(shg)
        gx_, gy_ = np.meshgrid(dxs, dys, indexing="xy")
        gx_ = gx_.ravel()
        gy_ = gy_.ravel()
        nb = gx_.size
        step = max(1, CHUNK_ROWS // nb)
        for s in range(0, idx.size, step):
            sub = idx[s:s + step]
            m = sub.size
            Ti = np.repeat(T[sub], nb, axis=0)
            xi = (i0[sub][:, None] + gx_[None, :]).ravel()
            yi = (j0[sub][:, None] + gy_[None, :]).ravel()
            a = tri_box_area(Ti, xi, yi)
            w = np.repeat(dens[sub], nb) * a
            flat = (yi.astype(np.int64) % H) * W + (xi.astype(np.int64) % W)
            acc += np.bincount(flat, weights=w, minlength=W * H).ravel()
            inside = ((xi >= 0) & (xi < W) & (yi >= 0) & (yi < H))
            info["area_in"] += float(a[inside].sum())
            info["area_out"] += float(a[~inside].sum())
            info["energy_in"] += float(w[inside].sum())
            info["energy_out"] += float(w[~inside].sum())
            info["bulk_boxes"] += int(m * nb)


def _deposit_large(T_one, dens_one, acc, W, H, info):
    """Tail path: huge / sliver triangles, rasterised row strip by row strip."""
    P = T_one
    ymin = float(P[:, 1].min())
    ymax = float(P[:, 1].max())
    j0 = int(np.floor(max(ymin, -2.0 * H)))
    j1 = int(np.floor(min(ymax, 3.0 * H)))
    if j1 < j0:
        return 0.0
    rows = np.arange(j0, j1 + 1, dtype=np.float64)
    R = rows.size
    V = np.repeat(P[None, :, :], R, axis=0)
    cnt = np.full(R, 3, dtype=np.int64)
    V, cnt = _clip(V, cnt, 0.0, 1.0, rows)
    V, cnt = _clip(V, cnt, 0.0, -1.0, -(rows + 1.0))
    slots = np.arange(V.shape[1], dtype=np.int64)[None, :]
    valid = slots < cnt[:, None]
    xs = np.where(valid, V[:, :, 0], np.inf)
    xe = np.where(valid, V[:, :, 0], -np.inf)
    c0 = np.floor(xs.min(1))
    c1 = np.floor(xe.max(1))
    c0 = np.maximum(c0, -2.0 * W)
    c1 = np.minimum(c1, 3.0 * W)
    # guard: a single row cannot need more than 4 frame widths
    big = (c1 - c0) > 4 * W
    c1 = np.where(big, c0 + 4 * W, c1)
    counts = np.maximum((c1 - c0 + 1).astype(np.int64), 0)
    if counts.sum() == 0:
        return 0.0
    ridx = np.repeat(np.arange(R), counts)
    start = np.repeat(np.cumsum(counts) - counts, counts)
    col = (c0[ridx] + (np.arange(counts.sum()) - start)).astype(np.float64)
    Vr = V[ridx]
    cntr = cnt[ridx]
    Vr, cntr = _clip(Vr, cntr, 1.0, 0.0, col)
    Vr, cntr = _clip(Vr, cntr, -1.0, 0.0, -(col + 1.0))
    a = _shoelace(Vr, cntr)
    w = a * dens_one
    flat = (rows[ridx].astype(np.int64) % H) * W + (col.astype(np.int64) % W)
    acc += np.bincount(flat, weights=w, minlength=W * H).ravel()
    inside = ((col >= 0) & (col < W) & (rows[ridx] >= 0) & (rows[ridx] < H))
    info["area_in"] += float(a[inside].sum())
    info["area_out"] += float(a[~inside].sum())
    info["energy_in"] += float(w[inside].sum())
    info["energy_out"] += float(w[~inside].sum())
    info["tail_boxes"] += int(a.size)
    return float(a.sum())


def deposit(T, A_tri, A_cell, w_cell, W, H, verbose=False):
    """Distribute every cell's energy by exact overlap area.

    T      (2N^2,3,2) projected triangles in unwrapped pixel coordinates
    A_tri  (2N^2,)    their |area| (px^2)
    A_cell (N^2,)     |A1|+|A2| per cell
    w_cell (N^2,)     energy mass per cell (dimensionless)

    Returns (acc (H,W) float64, info dict).  A cell whose two triangles have
    zero total area can carry no density and is recorded, not silently dropped:
    its energy is added to the uniform backstop.
    """
    W, H = int(W), int(H)
    acc = np.zeros(W * H, dtype=np.float64)
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

    _deposit_bulk(T[bulk], i0[bulk], j0[bulk], sw[bulk], sh[bulk],
                  dens[bulk], acc, W, H, info)

    # ---- tail: python loop, one triangle at a time (rare by construction) --
    tail_area = np.zeros(tail.size)
    for k, t in enumerate(tail):
        tail_area[k] = _deposit_large(T[t], float(dens[t]), acc, W, H, info)

    # ---- backstop: any triangle area not deposited (clipping guards) is
    # spread uniformly and REPORTED, never silently normalised away ----------
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


def mesh_C(fields, S, W, H, env, ior=IOR_G, w_total=None):
    """Area-conserving C for one time slice on the W x H frame grid."""
    N = np.asarray(fields["h"]).shape[0]
    PXe, PYe = landing_pixels(fields, S, env, ior=ior, W=W, H=H)
    T, A_tri, A_cell = cell_triangles(PXe, PYe)
    wt = float(W) * float(H) / (float(N) * float(N)) if w_total is None else float(w_total)
    w_cell = np.full(A_cell.size, wt)
    acc, info = deposit(T, A_tri, A_cell, w_cell, W, H)
    info.update({"N_s": int(N), "W": int(W), "H": int(H), "S_m": float(S),
                 "w_cell": wt, "ior": float(ior),
                 "method": "area push-forward (quad -> 2 triangles, exact "
                           "overlap with the frame pixels)",
                 "kernel": "none", "jitter": "none", "band_limit": "none",
                 "per_cell_energy_rule": "w_cell = W*H/N_s^2 (det.det_frame's "
                                         "divisor convention, rectangular "
                                         "frame analogue)"})
    if not np.all(np.isfinite(acc)):
        raise FloatingPointError("non-finite C")
    return acc.astype(np.float32), info


# ---------------------------------------------------------------------------
# reference lane: the pre-existing point splat, on the same frame grid
# ---------------------------------------------------------------------------
def point_splat_C(fields, S, W, H, env, ior=IOR_G, Q=25, sigma_px=0.75,
                  blur=True, row_block=256):
    """`det.det_frame`'s point-splat transport, evaluated on the W x H frame.

    Everything is det.det_frame's, verbatim in method: one ray per surface
    texel centre, Q equal-weight Fibonacci solar-disc directions, Snell
    refraction, bilinear splat with periodic wrap, analytical divisor
    (total_rays / pixels), then a sigma = 0.75 px periodic Gaussian.  The only
    difference is that the two frame axes are allowed to have different
    metre-per-pixel scales (dx = S/W, dy = S/H), which is the same
    parameterisation the area deposit uses; on a square grid this is
    numerically the legacy function (self test v).

    Returns (C (H,W) float32, info).  NOT used for the deliverable: it exists
    so the frame-24 peak statistics can be compared on one grid.
    """
    from det import fibonacci_disc
    from transport import orthonormal_basis, gaussian_blur_fft, sun_disc_radius

    h = np.asarray(fields["h"], dtype=np.float64)
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    Dx = np.asarray(fields.get("Dx", np.zeros_like(h)), dtype=np.float64)
    Dy = np.asarray(fields.get("Dy", np.zeros_like(h)), dtype=np.float64)
    N_s = h.shape[0]
    S = float(S)
    W, H = int(W), int(H)
    depth = float(env["depthM"])

    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    quads = fibonacci_disc(int(Q), disc)
    eta = 1.0 / float(ior)

    xs_1d = (np.arange(N_s, dtype=np.float64) + 0.5) * (S / N_s)
    acc = np.zeros(W * H, dtype=np.float64)
    for q in range(int(Q)):
        d = s + u * float(quads[q, 0]) + v * float(quads[q, 1])
        d = d / np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])
        for r0 in range(0, N_s, row_block):
            r1 = min(r0 + row_block, N_s)
            X0, Y0 = np.meshgrid(xs_1d, xs_1d[r0:r1])
            GX, GY = gx[r0:r1, :], gy[r0:r1, :]
            DX, DY = Dx[r0:r1, :], Dy[r0:r1, :]
            HH = h[r0:r1, :]
            nl = np.sqrt(GX * GX + GY * GY + 1.0)
            Nx, Ny, Nz = -GX / nl, -GY / nl, 1.0 / nl
            cos1 = np.clip(-(Ix * Nx + Iy * Ny + Iz * Nz), -1.0, 1.0)
            kk = np.maximum(1.0 - eta * eta * (1.0 - cos1 * cos1), 0.0)
            sq = np.sqrt(kk)
            Tx = eta * Ix + (eta * cos1 - sq) * Nx
            Ty = eta * Iy + (eta * cos1 - sq) * Ny
            Tz = np.minimum(eta * Iz + (eta * cos1 - sq) * Nz, -1e-6)
            t = (HH + depth) / (-Tz)
            xb = np.mod(X0 + DX + Tx * t, S)
            yb = np.mod(Y0 + DY + Ty * t, S)
            fx = (xb / S * W).ravel()
            fy = (yb / S * H).ravel()
            ix0 = np.floor(fx).astype(np.int64) % W
            iy0 = np.floor(fy).astype(np.int64) % H
            tx = fx - np.floor(fx)
            ty = fy - np.floor(fy)
            ix1 = (ix0 + 1) % W
            iy1 = (iy0 + 1) % H
            i00 = iy0 * W + ix0
            i10 = iy0 * W + ix1
            i01 = iy1 * W + ix0
            i11 = iy1 * W + ix1
            nb = W * H
            acc += np.bincount(i00, weights=(1.0 - tx) * (1.0 - ty), minlength=nb)
            acc += np.bincount(i10, weights=tx * (1.0 - ty), minlength=nb)
            acc += np.bincount(i01, weights=(1.0 - tx) * ty, minlength=nb)
            acc += np.bincount(i11, weights=tx * ty, minlength=nb)

    total_rays = float(int(Q)) * float(N_s) * float(N_s)
    divisor = total_rays / float(W * H)
    C = (acc / divisor).reshape(H, W)
    if blur:
        C = gaussian_blur_fft(C.astype(np.float64), float(sigma_px))
    info = {"method": "point splat (det.det_frame, same frame grid)",
            "Q": int(Q), "ior": float(ior), "disc_rad": float(disc),
            "sigma_px": float(sigma_px) if blur else 0.0,
            "total_rays": total_rays, "divisor_analytical": float(divisor),
            "N_s": int(N_s), "W": int(W), "H": int(H), "S_m": float(S),
            "kernel": "none (bilinear splat only)", "jitter": "none",
            "band_limit": "none", "normalisation": "none",
            "mean": float(C.mean())}
    return C.astype(np.float32), info


# ---------------------------------------------------------------------------
# self test
# ---------------------------------------------------------------------------
def _selftest():
    rng = np.random.default_rng(20261113)
    ok = True

    def chk(name, cond, detail=""):
        nonlocal ok
        ok = ok and bool(cond)
        print("[%s] %s %s" % ("ok" if cond else "NG", name, detail))

    # (i) flat surface: C must be exactly 1 everywhere
    N, W, H = 64, 736, 540
    S = 1.4755102040816326
    z = np.zeros((N, N), dtype=np.float32)
    f = {"h": z, "gx": z, "gy": z, "Dx": z, "Dy": z, "S": S}
    env = {"depthM": 0.6, "sunDir": np.array([0.4965974, -0.1350766, 0.8574062])}
    C, info = mesh_C(f, S, W, H, env)
    chk("flat surface -> C == 1", np.max(np.abs(C - 1.0)) < 1e-9,
        "max|C-1| = %.3e  rel_err = %.2e" % (np.max(np.abs(C - 1.0)), info["rel_err"]))
    chk("flat surface energy identity", info["rel_err"] < 1e-12,
        "rel_err = %.2e" % info["rel_err"])

    # (ii) a lattice-preserving map (shear): a torus bijection -> C == 1
    # exactly.  Every pixel is covered once and every cell has the same area,
    # so any deviation is a rasterisation error.
    axc = np.arange(N + 1) / N
    XC, YC = np.meshgrid(axc, axc, indexing="xy")
    PXe = (XC + 0.5 * YC) * W
    PYe = YC * H
    T, A_tri, A_cell = cell_triangles(PXe, PYe)
    wt = W * H / float(N * N)
    acc, info2 = deposit(T, A_tri, A_cell, np.full(A_cell.size, wt), W, H)
    chk("sheared lattice map -> C == 1", np.max(np.abs(acc - 1.0)) < 1e-9,
        "max|C-1| = %.3e" % np.max(np.abs(acc - 1.0)))
    chk("sheared lattice map energy identity", info2["rel_err"] < 1e-12,
        "rel_err = %.2e" % info2["rel_err"])

    # (ii-b) integer dilation J=2: the torus is covered 4x by cells of 4x the
    # area, so C is exactly 1 again (mean 1 is a conservation consequence).
    PXe = (0.5 + (XC - 0.5) * 2.0) * W
    PYe = (0.5 + (YC - 0.5) * 2.0) * H
    T, A_tri, A_cell = cell_triangles(PXe, PYe)
    acc2, info3 = deposit(T, A_tri, A_cell, np.full(A_cell.size, wt), W, H)
    chk("integer dilation J=2 -> C == 1", np.max(np.abs(acc2 - 1.0)) < 1e-9,
        "max|C-1| = %.3e" % np.max(np.abs(acc2 - 1.0)))

    # (ii-c) why the dilation test is stated as "== flat" and not "== 1/4":
    # C is a per-PIXEL mass, and the J=2 map is a degree-4 cover, so a pixel
    # sums four branches.  Each branch carries exactly 1/|det J| = 1/4 of the
    # flat density, the four sum back to the flat value.  A depositor that
    # forgot to divide by |det J| would give 4x here; one that reported one
    # branch would give 1/4.  Both are excluded by the two checks below.
    PXf, PYf = XC * W, YC * H
    Tf, Af_tri, Af_cell = cell_triangles(PXf, PYf)
    accf, _ = deposit(Tf, Af_tri, Af_cell, np.ones(Af_cell.size), W, H)
    chk("flat map, unit cell mass -> C == N^2/(W*H) (the flat unit)",
        np.max(np.abs(accf - N * N / float(W * H))) < 1e-9,
        "max|C-unit| = %.3e unit = %.9f"
        % (np.max(np.abs(accf - N * N / float(W * H))), N * N / float(W * H)))
    acc2b, _ = deposit(T, A_tri, A_cell, np.ones(A_cell.size), W, H)
    ratio = np.max(np.abs(acc2b / accf))
    chk("dilation J=2 (degree-4 cover): per-pixel mass == flat (not 4x, "
        "not 1/4x)", abs(ratio - 1.0) < 1e-9,
        "max C_J2/C_flat = %.12f ; per-branch density = 1/|det J| = %.3f"
        % (ratio, 1.0 / 4.0))
    chk("dilation J=2: A_cell == |det J| * A_flat",
        np.max(np.abs(A_cell / (W * H / float(N * N)) - 4.0)) < 1e-9,
        "max|A_cell/A_flat - 4| = %.3e"
        % np.max(np.abs(A_cell / (W * H / float(N * N)) - 4.0)))

    # (ii-d) RIGID rotation that IS a torus automorphism: 90 degrees about the
    # patch centre plus a shift (x' = 1-v+0.37, y' = u+0.23).  Square grid, so
    # the lattice is mapped to itself; the image tiles the torus once (and the
    # shift makes it wrap).  This replaces the old `rigid rotation by 0.41 rad'
    # check, which could never be exact: a rotation by a non-lattice angle does
    # not map Z^2 to Z^2, so it is NOT a bijection of the torus (see ii-e).
    Nsq, Wsq = 64, 368
    aqs = np.arange(Nsq + 1) / Nsq
    XS, YS = np.meshgrid(aqs, aqs, indexing="xy")
    PXs = (1.0 - YS + 0.37) * Wsq
    PYs = (XS + 0.23) * Wsq
    Ts, As_tri, As_cell = cell_triangles(PXs, PYs)
    accs, info4 = deposit(Ts, As_tri, As_cell,
                          np.full(As_cell.size, Wsq * Wsq / float(Nsq * Nsq)),
                          Wsq, Wsq)
    chk("90 deg rigid rotation + shift (torus automorphism) -> C == 1",
        np.max(np.abs(accs - 1.0)) < 1e-9,
        "max|C-1| = %.3e  rel_err = %.2e  wraps = %s"
        % (np.max(np.abs(accs - 1.0)), info4["rel_err"],
           info4["offframe_fraction_energy"] > 0))

    # (ii-e) control: a rotation by a NON-lattice angle.  Energy must still be
    # exact (mass is conserved by construction) but C must NOT be 1 -- the
    # rotated square is not a fundamental domain of the axis-aligned lattice,
    # so the pushed-forward measure overlaps itself.  Any `exactness' claim for
    # such a map is a statement about the map, not about the depositor.
    th = 0.41
    XC0, YC0 = XS - 0.5, YS - 0.5
    PXr = (0.5 + np.cos(th) * XC0 - np.sin(th) * YC0) * Wsq
    PYr = (0.5 + np.sin(th) * XC0 + np.cos(th) * YC0) * Wsq
    Tr2, Ar_tri, Ar_cell = cell_triangles(PXr, PYr)
    accr, infor = deposit(Tr2, Ar_tri, Ar_cell,
                          np.full(Ar_cell.size, Wsq * Wsq / float(Nsq * Nsq)),
                          Wsq, Wsq)
    dev = float(np.max(np.abs(accr - 1.0)))
    chk("non-lattice rotation: energy exact", infor["rel_err"] < 1e-12,
        "rel_err = %.2e  sum = %.6f (W*H = %d)"
        % (infor["rel_err"], accr.sum(), Wsq * Wsq))
    chk("non-lattice rotation is NOT a bijection (C != 1)", dev > 0.2,
        "max|C-1| = %.3f, empty pixels = %.4f, area excess = %.3e px^2"
        % (dev, float(np.mean(accr < 1e-12)),
           float(As_cell.sum() - Ar_cell.sum())))

    # (iii) random small triangles vs dense supersampling
    M = 400
    c = rng.uniform(2, 600, size=(M, 2))
    th2 = rng.uniform(0, 2 * np.pi, size=(M, 3, 1))
    rr = rng.uniform(0.05, 1.2, size=(M, 3, 1))
    Tr = c[:, None, :] + rr * np.stack([np.cos(th2[..., 0]), np.sin(th2[..., 0])], axis=-1)
    x0 = np.floor(Tr[:, :, 0].min(1)) - 1.0
    y0 = np.floor(Tr[:, :, 1].min(1)) - 1.0
    Ssub = 256
    sub = (np.arange(Ssub) + 0.5) / Ssub
    ref = np.zeros(M)
    for k in range(M):
        px = Tr[k]
        gx = x0[k] + 4.0 * sub[None, :]
        gy = y0[k] + 4.0 * sub[:, None]
        d1 = (px[1, 0] - px[0, 0]) * (gy - px[0, 1]) - (px[1, 1] - px[0, 1]) * (gx - px[0, 0])
        d2 = (px[2, 0] - px[1, 0]) * (gy - px[1, 1]) - (px[2, 1] - px[1, 1]) * (gx - px[1, 0])
        d3 = (px[0, 0] - px[2, 0]) * (gy - px[2, 1]) - (px[0, 1] - px[2, 1]) * (gx - px[2, 0])
        sg = np.sign(d1) + np.sign(d2) + np.sign(d3)
        ref[k] = (np.abs(sg) >= 2).sum() / (Ssub * Ssub) * 16.0
    act = np.zeros(M)
    for dx in range(4):
        for dy in range(4):
            act += tri_box_area(Tr, x0 + dx, y0 + dy)
    d = np.abs(act - ref)
    chk("exact clip vs 256^2 supersample", d.max() < 0.01,
        "max|diff| = %.2e px^2 (mean %.2e)" % (d.max(), d.mean()))

    # (iv) additivity: the two triangles of a cell sum to the cell area
    chk("triangle areas add to the cell area", True,
        "max|A1+A2-A_cell| = %.2e (by construction)"
        % np.max(np.abs((A_tri[:N * N] + A_tri[N * N:]) - A_cell)))

    # (v) the reference lane (point splat) is the legacy transport
    Ns = 64
    ax2 = np.arange(Ns) / float(Ns)
    X2, Y2 = np.meshgrid(ax2, ax2, indexing="xy")
    hh = (0.004 * np.sin(2 * np.pi * 3 * X2) * np.cos(2 * np.pi * 2 * Y2)).astype(np.float32)
    gxx = np.gradient(np.asarray(hh, dtype=np.float64), axis=1) * Ns / S
    gyy = np.gradient(np.asarray(hh, dtype=np.float64), axis=0) * Ns / S
    fr = {"h": hh, "gx": gxx.astype(np.float32), "gy": gyy.astype(np.float32),
          "Dx": np.zeros_like(hh), "Dy": np.zeros_like(hh), "S": S}
    Cps, _ = point_splat_C(fr, S, 256, 256, env, Q=25)
    from det import det_frame as _det_frame
    Cdet, _ = _det_frame(fr, S, 256, env, 25)
    chk("point splat == det.det_frame on a square grid",
        np.max(np.abs(Cps - Cdet)) < 1e-5,
        "max|diff| = %.3e (C mean %.6f)" % (np.max(np.abs(Cps - Cdet)),
                                            float(np.mean(Cps))))
    f512 = {"h": np.zeros((512, 512), dtype=np.float32),
            "gx": np.zeros((512, 512), dtype=np.float32),
            "gy": np.zeros((512, 512), dtype=np.float32),
            "Dx": np.zeros((512, 512), dtype=np.float32),
            "Dy": np.zeros((512, 512), dtype=np.float32), "S": S}
    Cpf, infof = point_splat_C(f512, S, W, H, env, Q=25)
    chk("point splat on flat water: mass exact (mean == 1)",
        abs(infof["mean"] - 1.0) < 1e-12,
        "mean = %.12f" % infof["mean"])
    print("[--] point splat on flat water: max|C-1| = %.3f  "
          "(lattice ripple: %d cells onto %d pixels; the area deposit gives "
          "exactly 0 on the same field, check (i))"
          % (np.max(np.abs(Cpf - 1.0)), 512 * 512, W * H))
    print("-" * 60)
    print("SELFTEST_OK" if ok else "SELFTEST_FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    print(__doc__)
