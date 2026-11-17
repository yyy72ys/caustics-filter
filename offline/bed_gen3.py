"""bed_gen3 -- the grain-layer mechanisms of v4-2.15 (mechanism A / mechanism B).

This module is the FIRST implementation of the design's grain layer.  It is a
library: every mechanism is a standalone function that takes a `Bed` and
returns a new `Bed` (never mutates the input) plus a diagnostics dict, so the
factorial driver (run_v4_2_15.py) can switch each step on/off independently.

A-0  initial_bed()      bed_gen2's placement, verbatim (patch box), or four
                        verbatim 0.6 m populations laid side by side on one
                        rectangular form (strip box; v4-2-3's precedent).
A-1  winnow()           tau(x,y) = rho g H S0 times a wake-deficit field behind
                        the coarse clasts; per-grain exposure p_i; critical
                        exposure p_c from a moment balance that involves
                        tau*, D_i, the pivot angle and the overburden; removal
                        with probability W*(phi); removed grains leave the
                        system (no replacement).
A-2  granular_sorting() size-selective vertical rearrangement (percolation of
                        d < Dc*/3, bridging of Dc*/3 <= d < Dc*, squeeze
                        expulsion of the coarse), mass kept in the box.
A-3  interlock()        imbrication (upstream edge buried), 3-point support
                        settling of the surface grains, cluster cover to ~10%.
B-1  facetise()         Voronoi-style fragmentation: each ellipsoid becomes a
                        convex polyhedron with 6-12 flat faces, same calipers
                        (b/a, Corey unchanged by construction), sharp corners.
B-2  planation()        local mean + k*sigma_z surface; 'clip' cuts the heads
                        (grains remain) or 'remove' deletes the protruding
                        grains and lets the bed below show.
B-3  pocket_fill()      fines (< Dc*/3) into the inter-grain depressions.

The 2-D height field is the top-most surface (vertical ray cast, periodic),
same z-buffer rule as bed_gen2.  Faceted grains are rasterised with their
exact half-space intersection (top face = min over the upper faces), which is
why the facet change is a change of *shape*, not a change of placement.

Declared constants (physics, not fitted) are at the top.  The one calibrated
number is TAU_STAR_DENSE (0.05), the Shields value at which a fully exposed
D50 grain starts to move; it anchors the moment balance and is stated in the
report.
"""
from __future__ import annotations

import copy
import math

import numpy as np

import bed_gen as BG
import bed_gen2 as B2

# --- water / sediment (declared) -------------------------------------------
RHO_W = 1000.0
RHO_S = 2650.0
DRHO = RHO_S - RHO_W
G = 9.81
KAPPA = 0.4
C_D = 1.2
C_L = 0.4
Z0_OVER_D50 = 1.0 / 30.0
H_FLOW = 0.25                     # v4-2.2 scene
S0 = 0.003                        # v4-2.2 scene
TAU0_SCENE = RHO_W * G * H_FLOW * S0          # 7.3575 Pa
TAU_STAR_SWEEP = (0.03, 0.045, 0.06, 0.08, 0.10)
TAU_STAR_DENSE = 0.05             # calibration anchor (fully exposed D50)
PHI_DEG = (25.0, 40.0)            # pivot angle range
GOLDEN = math.pi * (3.0 - math.sqrt(5.0))

# --- vector contract -------------------------------------------------------
PATCH = B2.PATCH
CELL = B2.PATCH / B2.HIRES        # 0.29296875 mm
INPAINT_PX = 60.0


# ===========================================================================
# geometry helpers
# ===========================================================================
def fibonacci_sphere(n):
    i = np.arange(n, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = np.mod(GOLDEN * np.arange(n), 2.0 * math.pi)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def _ellipsoid_radius(dirs, abc):
    A, B, C = abc
    q = (dirs[:, 0] / A) ** 2 + (dirs[:, 1] / B) ** 2 + (dirs[:, 2] / C) ** 2
    return 1.0 / np.sqrt(q)


def _plane_vertices(N, off):
    """Vertices of {N p <= off} for a bounded convex set (via the dual hull).

    q_k = n_k / off_k; the convex hull of q has one facet per vertex of P.
    """
    from scipy.spatial import ConvexHull
    q = N / off[:, None]
    hull = ConvexHull(q)
    verts = []
    for tri in hull.simplices:
        M = N[tri]
        if abs(np.linalg.det(M)) < 1e-14:
            continue
        verts.append(np.linalg.solve(M, off[tri]))
    if not verts:
        return np.zeros((0, 3))
    V = np.array(verts)
    # keep only points inside every half space
    keep = np.all(V @ N.T <= off[None, :] + 1e-9, axis=1)
    V = V[keep]
    # unique
    if len(V):
        _, idx = np.unique(np.round(V, 10), axis=0, return_index=True)
        V = V[np.sort(idx)]
    return V


def _face_polygons(V, N, off, tol=1e-7):
    """For each face, the polygon (ordered vertices) and its area/edge lengths."""
    out = []
    for k in range(len(N)):
        d = V @ N[k] - off[k]
        on = np.abs(d) < 1e-6
        P = V[on]
        if len(P) < 3:
            out.append({"area": 0.0, "edges": [], "n": int(len(P))})
            continue
        c = P.mean(axis=0)
        n = N[k] / np.linalg.norm(N[k])
        e1 = P[0] - c
        e1 = e1 - (e1 @ n) * n
        if np.linalg.norm(e1) < 1e-12:
            e1 = np.cross(n, np.array([1.0, 0.0, 0.0]))
        e1 = e1 / max(np.linalg.norm(e1), 1e-12)
        e2 = np.cross(n, e1)
        ang = np.arctan2((P - c) @ e2, (P - c) @ e1)
        o = np.argsort(ang)
        P = P[o]
        area = 0.0
        edges = []
        for i in range(len(P)):
            a = P[i] - c
            b = P[(i + 1) % len(P)] - c
            area += 0.5 * np.linalg.norm(np.cross(a, b))
            edges.append(float(np.linalg.norm(P[(i + 1) % len(P)] - P[i])))
        out.append({"area": float(area), "edges": edges, "n": int(len(P))})
    return out


def poly_metrics(V, N, off):
    """Volume, calipers, facet statistics and the two roundness proxies."""
    polys = _face_polygons(V, N, off)
    vol = sum(p["area"] * float(off[k]) for k, p in enumerate(polys)) / 3.0
    cal = np.max(np.abs(V), axis=0) if len(V) else np.zeros(3)
    area = sum(p["area"] for p in polys)
    # facet linear dimension: sqrt(face area) and mean edge
    lin = [math.sqrt(p["area"]) for p in polys if p["area"] > 0]
    edges = [e for p in polys for e in p["edges"]]
    Rc = float(np.max(np.linalg.norm(V, axis=1))) if len(V) else 0.0
    deq = (6.0 * vol / math.pi) ** (1.0 / 3.0) if vol > 0 else 0.0
    wadell_vol = deq / (2.0 * Rc) if Rc > 0 else 0.0
    # top-view silhouette corner angles -> Powers/Krumbein-style roundness
    ang = None
    if len(V) >= 4:
        from scipy.spatial import ConvexHull
        try:
            h = ConvexHull(V[:, :2])
            P = V[h.vertices][:, :2]
            th = []
            for i in range(len(P)):
                a = P[(i - 1) % len(P)] - P[i]
                b = P[(i + 1) % len(P)] - P[i]
                ca = a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-15)
                th.append(math.acos(max(-1.0, min(1.0, ca))))
            ang = float(np.mean(th))
        except Exception:
            ang = None
    return {"volume": float(vol), "calipers": [float(c) for c in cal],
            "area": float(area), "facet_areas": [p["area"] for p in polys],
            "facet_lin": lin, "edges": edges,
            "n_faces": int(len([p for p in polys if p["area"] > 0])),
            "circumradius": Rc, "deq": float(deq),
            "wadell_volume_proxy": float(wadell_vol),
            "corner_angle_mean_deg": (None if ang is None
                                      else float(math.degrees(ang))),
            "powers_roundness": (None if ang is None
                                 else float(0.6 * np.clip((math.degrees(ang) - 90.0) / 90.0, 0.0, 1.0)))}


def make_facet_faces(rng, abc, n_faces=None):
    """Voronoi-style fracturing: the ellipsoid cut by n flat fracture planes.

    Each plane is a shallow cut near a random surface direction
    (offset = u * r_ellipsoid(n) with u in [0.87, 0.96]); the grain is the
    ellipsoid intersected with the n half spaces.  The result is a pebble
    whose original rounded surface is mostly preserved and which carries
    6-12 flat, sharp-edged breakage faces -- the design's B-1 object.

    Returned planes are in the grain's LOCAL frame (metres).  Volume and facet
    statistics are computed here: the facet is the ellipse where the plane
    cuts the ellipsoid (60-point polygon, area and interior angles).
    """
    A, B, C = abc
    if n_faces is None:
        n_faces = int(rng.integers(6, 13))
    dirs = fibonacci_sphere(n_faces)
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    Rq = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                   [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                   [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    dirs = dirs @ Rq.T
    # break the Fibonacci symmetry of the offsets
    dirs = dirs + rng.normal(0.0, 0.06, size=dirs.shape)
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    r = _ellipsoid_radius(dirs, abc)
    u = rng.uniform(0.87, 0.96, size=n_faces)
    off = r * u
    N = dirs
    # ---- facet polygons: plane n.x = off  ∩  ellipsoid ---------------
    polys = []
    for k in range(n_faces):
        n = N[k]
        # affine map to the unit sphere: x = (A y0, B y1, C y2); the plane
        # n.x = off becomes  m.y = off  with  m = (A n0, B n1, C n2)
        m = np.array([n[0] * A, n[1] * B, n[2] * C])
        # plane in x-space: n.x = off  ->  m.(D x) = off, |Dx|=1
        Mm = np.linalg.norm(m)
        if Mm <= 1e-12:
            continue
        h = float(off[k]) / Mm
        if abs(h) >= 1.0:
            fran = 0.0
            continue
        e1 = np.array([1.0, 0, 0]) if abs(m[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = e1 - (e1 @ m) / Mm ** 2 * m
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(m / Mm, e1)
        rho = math.sqrt(max(1.0 - h * h, 0.0))
        th = np.linspace(0.0, 2 * math.pi, 60, endpoint=False)
        yy = (h / Mm) * m + rho * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
        xx = np.stack([A * yy[:, 0], B * yy[:, 1], C * yy[:, 2]], axis=1)
        c = xx.mean(axis=0)
        dd = xx - c
        cr = np.cross(dd, np.roll(dd, -1, axis=0))
        ar = 0.5 * float(np.linalg.norm(cr.sum(axis=0)))
        # dihedral angle between the flat facet and the smooth surface at
        # the cut boundary (0 = flush/tangent, pi/2 = knife edge)
        ne = np.stack([xx[:, 0] / A ** 2, xx[:, 1] / B ** 2,
                       xx[:, 2] / C ** 2], axis=1)
        ne /= np.linalg.norm(ne, axis=1)[:, None]
        cosd = np.clip(np.abs(ne @ n), 0.0, 1.0)
        dih = float(np.degrees(np.arccos(cosd).mean()))
        polys.append({"area": float(ar), "lin": math.sqrt(max(ar, 0.0)),
                      "dihedral_deg": dih,
                      "cap_height": float(off[k] / Mm)})
    # ---- volume: unit-sphere lattice, cut by the same planes ----------
    g1 = (np.arange(12) + 0.5) / 12.0 * 2 - 1
    GX, GY, GZ = np.meshgrid(g1, g1, g1)
    sel = GX ** 2 + GY ** 2 + GZ ** 2 <= 1.0
    Y = np.stack([GX[sel], GY[sel], GZ[sel]], axis=1)
    X = np.stack([A * Y[:, 0], B * Y[:, 1], C * Y[:, 2]], axis=1)
    inside = np.all(X @ N.T <= off[None, :] + 1e-12, axis=1)
    Vell = 4.0 / 3.0 * math.pi * A * B * C
    vol = float(Vell * inside.mean())
    lin = [p["lin"] for p in polys]
    dihs = [p["dihedral_deg"] for p in polys]
    Rc = float(max(A, B, C))
    deq = (6.0 * vol / math.pi) ** (1.0 / 3.0)
    # Knud Thomsen ellipsoid area (p = 1.6)
    Ap, Bp, Cp = A ** 1.6, B ** 1.6, C ** 1.6
    area_ell = 4.0 * math.pi * ((Ap * Bp + Ap * Cp + Bp * Cp) / 3.0) ** (1.0 / 1.6)
    area_frac = float(sum(p["area"] for p in polys) / max(area_ell, 1e-18))
    dm = float(np.mean(dihs)) if dihs else None
    m = {"volume": vol, "calipers": [2 * A, 2 * B, 2 * C],
         "facet_areas": [p["area"] for p in polys],
         "facet_lin": lin, "n_faces": int(len(polys)),
         "facet_area_fraction": area_frac,
         "facet_dihedral_deg": dm,
         "wadell_volume_proxy": float(deq / (2 * Rc)),
         "facet_edge_roundness": (None if dm is None
                                  else float(0.6 * np.clip(1.0 - dm / 90.0, 0.0, 1.0))),
         "facet_edge_roundness_definition":
             "0.6*clip(1 - mean(facet/surface dihedral angle)/90 deg, 0, 1); "
             "0 = knife-edge corner (angular), 0.6 = flush facet (rounded), "
             "Powers/Krumbein visual scale"}
    return {"N": N, "off": off, "metrics": m}


def _rot_from_dip_az(dip_deg, az_deg):
    """Local frame with the a-axis at az and the a-b plane dipping dip_deg.

    nrm = (-sin d, 0, cos d) after the azimuth rotation -> for flow +x the
    UPSTREAM edge of the a-b plane is the low one (imbrication).
    """
    d = math.radians(dip_deg)
    az = math.radians(az_deg)
    Rz = np.array([[math.cos(az), -math.sin(az), 0.0],
                   [math.sin(az), math.cos(az), 0.0],
                   [0.0, 0.0, 1.0]])
    nrm = np.array([-math.sin(d), 0.0, math.cos(d)])
    ea = np.array([math.cos(d), 0.0, math.sin(d)])
    eb = np.cross(nrm, ea)
    return np.stack([ea, eb, nrm], axis=1) @ Rz.T


def _world_faces(g):
    """(k,4) world-frame planes [nx, ny, nz, off] of a fractured grain."""
    F = g["frac"]
    Nw = F["N"] @ g["R"].T
    return np.column_stack([Nw, F["off"]])


# ===========================================================================
# bed container
# ===========================================================================
def clone_grains(grains):
    out = []
    for g in grains:
        h = dict(g)
        # immutable/np arrays are shared on purpose (never mutated)
        out.append(h)
    return out


class Bed:
    def __init__(self, Lx, Ly, cell, grains, form, meta=None):
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.cell = float(cell)
        self.nx = int(round(self.Lx / self.cell))
        self.ny = int(round(self.Ly / self.cell))
        self.grains = grains
        self.form = form
        self.meta = dict(meta or {})
        self._z = None
        self._gid = None
        self.clip = None            # optional planation surface (ny,nx)
        self.extra = []             # extra fill grains (kept in .grains too)
        self.tau = None             # last tau field (coarse grid) for the record

    def copy(self):
        b = Bed(self.Lx, self.Ly, self.cell, clone_grains(self.grains),
                self.form, self.meta)
        b.clip = None if self.clip is None else self.clip.copy()
        b.extra = list(self.extra)
        return b

    @property
    def xs(self):
        return (np.arange(self.nx) + 0.5) * self.cell

    @property
    def ys(self):
        return (np.arange(self.ny) + 0.5) * self.cell

    # ---- raster ----------------------------------------------------------
    def render(self, force=False):
        if self._z is not None and not force:
            return self._z, self._gid
        z, gid = raycast_zgid(self.grains, self.xs, self.ys, self.nx, self.ny)
        hole = gid < 0
        if hole.any():
            hf = float(hole.mean())
            if hf > 1e-4:
                sig = INPAINT_PX if hf > 0.02 else 15.0
                num = _gauss(np.where(hole, 0.0, z), sig)
                den = _gauss((~hole).astype(np.float64), sig)
                z = np.where(hole, num / np.maximum(den, 1e-9), z)
        if self.clip is not None:
            z = np.minimum(z, self.clip)
        self._z, self._gid = z, gid
        return z, gid


def _gauss(a, sigma_px):
    from scipy import ndimage
    return ndimage.gaussian_filter(a, sigma_px, mode="wrap")


def raycast_zgid(grains, xs, ys, nx, ny):
    """bed_gen2.raycast2's z-buffer, without the normals (same arithmetic)."""
    z = np.full((ny, nx), -np.inf)
    gid = np.full((ny, nx), -1, dtype=np.int32)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    for gi, g in enumerate(grains):
        if "frac" in g:
            _stamp_frac(z, gid, gi, g, xs, ys, nx, ny)
            continue
        A, B, C = g["abc"]
        R = g["R"]
        ax = np.array([A, B, C])
        hx = BG._half_extent(g, 0)
        hy = BG._half_extent(g, 1)
        i0 = int(math.floor((g["x"] - hx - xs[0]) / dx))
        i1 = int(math.ceil((g["x"] + hx - xs[0]) / dx)) + 1
        j0 = int(math.floor((g["y"] - hy - ys[0]) / dy))
        j1 = int(math.ceil((g["y"] + hy - ys[0]) / dy)) + 1
        ixw = np.arange(i0, i1) % nx
        iyw = np.arange(j0, j1) % ny
        DX = (xs[ixw] - g["x"])[None, :]
        DY = (ys[iyw] - g["y"])[:, None]
        u = np.stack([R[0, i] * DX + R[1, i] * DY for i in range(3)])
        w = R[2, :][:, None, None]
        aa = ((w / ax[:, None, None]) ** 2).sum(0)
        bb = ((u * w) / ax[:, None, None] ** 2).sum(0)
        cc = ((u / ax[:, None, None]) ** 2).sum(0) - 1.0
        disc = bb * bb - aa * cc
        ok = disc > 0.0
        if not ok.any():
            continue
        sq = np.sqrt(np.where(ok, disc, 0.0))
        s2 = (-bb + sq) / aa
        ztop = g["z"] + s2
        idx = np.ix_(iyw, ixw)
        sub = z[idx]
        hit = ok & (ztop > sub)
        if not hit.any():
            continue
        z[idx] = np.where(hit, ztop, sub)
        gsub = gid[idx]
        gid[idx] = np.where(hit, gi, gsub)
    if (gid < 0).any():
        zz = np.where(np.isfinite(z), z, np.nan)
        z = np.where(np.isfinite(z), z, np.nanmean(zz))
    return z, gid


def _stamp_frac(z, gid, gi, g, xs, ys, nx, ny):
    """Ellipsoid top, clipped by the fracture planes (exact vertical ray)."""
    A, B, C = g["abc"]
    R = g["R"]
    ax = np.array([A, B, C])
    hx = BG._half_extent(g, 0)
    hy = BG._half_extent(g, 1)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    i0 = int(math.floor((g["x"] - hx - xs[0]) / dx))
    i1 = int(math.ceil((g["x"] + hx - xs[0]) / dx)) + 1
    j0 = int(math.floor((g["y"] - hy - ys[0]) / dy))
    j1 = int(math.ceil((g["y"] + hy - ys[0]) / dy)) + 1
    ixw = np.arange(i0, i1) % nx
    iyw = np.arange(j0, j1) % ny
    DX = (xs[ixw] - g["x"])[None, :]
    DY = (ys[iyw] - g["y"])[:, None]
    u = np.stack([R[0, i] * DX + R[1, i] * DY for i in range(3)])
    w = R[2, :][:, None, None]
    aa = ((w / ax[:, None, None]) ** 2).sum(0)
    bb = ((u * w) / ax[:, None, None] ** 2).sum(0)
    cc = ((u / ax[:, None, None]) ** 2).sum(0) - 1.0
    disc = bb * bb - aa * cc
    ok = disc > 0.0
    if not ok.any():
        return
    sq = np.sqrt(np.where(ok, disc, 0.0))
    ztop = g["z"] + (-bb + sq) / aa
    zlow = g["z"] + (-bb - sq) / aa
    F = _world_faces(g)
    for k in range(len(F)):
        nz = F[k, 2]
        rhs = F[k, 3] - F[k, 0] * DX - F[k, 1] * DY
        if nz > 1e-9:
            ztop = np.minimum(ztop, g["z"] + rhs / nz)
        elif nz < -1e-9:
            zlow = np.maximum(zlow, g["z"] + rhs / nz)
    hit = ok & (ztop > zlow) & np.isfinite(ztop)
    if not hit.any():
        return
    idx = np.ix_(iyw, ixw)
    sub = z[idx]
    hit = hit & (ztop > sub)
    if not hit.any():
        return
    z[idx] = np.where(hit, ztop, sub)
    gsub = gid[idx]
    gid[idx] = np.where(hit, gi, gsub)


# ===========================================================================
# A-0 -- initial bed
# ===========================================================================
def sigma_form_ref(stats_path):
    import json
    st = json.load(open(stats_path, encoding="utf-8"))
    return float(st["sigma_z_fixpoint_final_parameters"]["sigma_form_m"])


def form_field_rect(nx, ny, Lx, Ly, sigma, seed, w_low=B2.FORM_W_LOW,
                    skew_c=B2.FORM_SKEW_C, H=B2.SUBSTRATE_H):
    """Rectangular version of bed_gen2.form_field (same two layers, same
    quadratic skew).  Used only for the strip box, which bed_gen2 cannot make.
    """
    rng = np.random.default_rng(int(seed) + 7)
    w = rng.standard_normal((ny, nx))
    kx = np.fft.fftfreq(nx, d=Lx / nx) * 2.0 * math.pi
    ky = np.fft.fftfreq(ny, d=Ly / ny) * 2.0 * math.pi
    KX, KY = np.meshgrid(kx, ky)
    K = np.hypot(KX, KY)
    amp = np.zeros_like(K)
    m = K > 0
    amp[m] = K[m] ** (-(H + 1.0))
    f = np.real(np.fft.ifft2(np.fft.fft2(w) * amp))
    f = f - f.mean()
    f = f / max(f.std(), 1e-12)
    # low modes, rectangular
    rng2 = np.random.default_rng(int(seed) + 11)
    ii = (np.arange(nx) + 0.5) * (Lx / nx)
    jj = (np.arange(ny) + 0.5) * (Ly / ny)
    X, Y = np.meshgrid(ii, jj)
    u = np.zeros((ny, nx))
    for i in range(-2, 3):
        for j in range(-2, 3):
            k2 = i * i + j * j
            if k2 == 0 or k2 > 4:
                continue
            a = k2 ** (-(H + 1.0) / 2.0)
            ph = rng2.uniform(0.0, 2.0 * math.pi)
            u += a * np.cos(2.0 * math.pi * (i * X / Lx + j * Y / Ly) + ph)
    u = u - u.mean()
    u = u / max(u.std(), 1e-12)
    g = w_low * u + math.sqrt(max(1.0 - w_low * w_low, 0.0)) * f
    g = g / max(g.std(), 1e-12)
    v = g + skew_c * (g * g - 1.0)
    v -= v.mean()
    return (v * (sigma / max(v.std(), 1e-12))).astype(np.float64)


def _place_on_form(grains, form, Lx, Ly):
    ny, nx = form.shape
    cellx = Lx / nx
    celly = Ly / ny
    Xg = np.array([g["x"] for g in grains])
    Yg = np.array([g["y"] for g in grains])
    fx = np.asarray(BG._bilinear(form, (np.mod(Xg, Lx) / cellx)[None, :],
                                 (np.mod(Yg, Ly) / celly)[None, :])[0])
    hz = np.array([BG._half_extent(g, 2) for g in grains])
    p = np.array([g["p"] for g in grains])
    z = fx + p - hz
    for i, g in enumerate(grains):
        g["z"] = float(z[i])


def initial_bed(seed, box="patch", sigma_form=0.0106289907,
                tile_seeds=None):
    """A-0.  bed_gen2's grain placement, unchanged.

    patch : B2.form_field(512, sigma_form, seed) + B2.get_population(seed,
            B2.P_SCALE0) + place(); the two mono boulders are dropped (they are
            shadow probes, not bed) and their footprint is inpainted, exactly
            the v4-2-2 rule.
    strip : 4 verbatim 0.6 m populations laid side by side (v4-2-3's
            precedent) on one rectangular form; the packing is not seam-aware
            (disclosed, same as v4-2-3).
    """
    if box == "patch":
        Lx = Ly = PATCH
        form = B2.form_field(512, sigma_form, seed)
        grains, _cm = B2.get_population(seed, B2.P_SCALE0)
        grains = [g for g in grains if not str(g["kind"]).startswith("mono")]
        B2.place(grains, form)
        meta = {"box": "patch", "form_seed": int(seed),
                "form_grid": [int(v) for v in form.shape],
                "populations": [int(seed)]}
    elif box == "strip":
        Lx, Ly = 2.4, 0.6
        nfx, nfy = 2048, 512
        seeds = list(tile_seeds) if tile_seeds else [
            int(seed) + 10 * k for k in range(4)]
        form = form_field_rect(nfx, nfy, Lx, Ly, sigma_form, seed)
        grains = []
        for k, sd in enumerate(seeds):
            gs, _cm = B2.get_population(int(sd), B2.P_SCALE0)
            gs = [g for g in gs if not str(g["kind"]).startswith("mono")]
            for g in gs:
                g["x"] = float((g["x"] + k * PATCH) % Lx)
                g["y"] = float(g["y"] % Ly)
            grains.extend(gs)
        _place_on_form(grains, form, Lx, Ly)
        meta = {"box": "strip", "form_seed": int(seed),
                "form_grid": [int(v) for v in form.shape],
                "populations": [int(s) for s in seeds],
                "seam_note": "packing is not aware of the tile seams; "
                             "grains may overlap across the three seam lines"}
    else:
        raise ValueError(box)
    b = Bed(Lx, Ly, CELL, grains, form, meta)
    b.meta["sigma_form_m"] = float(sigma_form)
    b.meta["seed"] = int(seed)
    b.meta["construction"] = "bed_gen2 population, verbatim"
    return b


# ===========================================================================
# shared measurements on a bed
# ===========================================================================
def grain_table(bed, z, gid):
    """Per-grain exposed flag, exposed px, top and depth below the surface."""
    n = len(bed.grains)
    valid = gid >= 0
    counts = np.bincount(gid[valid].ravel(), minlength=n).astype(np.int64)
    D = np.array([g["D"] for g in bed.grains], dtype=np.float64)
    kind = np.array([str(g["kind"]) for g in bed.grains])
    hz = np.array([BG._half_extent(g, 2) for g in bed.grains])
    top = np.array([g["z"] for g in bed.grains]) + hz
    cluster = np.array([k.startswith("cluster") for k in kind])
    exposed = counts > 0
    # local surface height at the grain centre
    xs = bed.cell
    zi = np.clip((np.array([g["x"] for g in bed.grains]) / xs).astype(int),
                 0, bed.nx - 1)
    yi = np.clip((np.array([g["y"] for g in bed.grains]) / xs).astype(int),
                 0, bed.ny - 1)
    zs = z[yi, zi]
    depth = zs - top
    return {"counts": counts, "D": D, "kind": kind, "cluster": cluster,
            "exposed": exposed, "top": top, "z_surface": zs,
            "depth_below_surface": depth}


def _wq(values, weights, q):
    o = np.argsort(values)
    v = values[o]
    w = weights[o]
    c = np.cumsum(w) / max(w.sum(), 1e-30)
    return float(np.interp(q, c, v))


def gsd_metrics(bed, z, gid):
    t = grain_table(bed, z, gid)
    cell2 = bed.cell * bed.cell
    ok = ~t["cluster"]
    D = t["D"]
    exposed = t["exposed"] & ok
    buried = (~t["exposed"]) & ok
    area = t["counts"] * cell2
    foot = np.array([BG.footprint_area(g) for g in bed.grains])
    out = {"n_exposed": int(exposed.sum()), "n_buried": int(buried.sum()),
           "exposed_area_m2": float(area[exposed].sum())}
    if exposed.sum() >= 10:
        out["surface_D50_m"] = _wq(D[exposed], area[exposed], 0.50)
        out["surface_D84_m"] = _wq(D[exposed], area[exposed], 0.84)
        out["surface_D16_m"] = _wq(D[exposed], area[exposed], 0.16)
        out["surface_D90_m"] = _wq(D[exposed], area[exposed], 0.90)
        out["surface_sigma_g"] = (out["surface_D84_m"] / out["surface_D50_m"]
                                  if out["surface_D50_m"] > 0 else None)
    if buried.sum() >= 10:
        out["buried_D50_m"] = _wq(D[buried], foot[buried], 0.50)
        out["buried_D84_m"] = _wq(D[buried], foot[buried], 0.84)
    if out.get("surface_D50_m") and out.get("buried_D50_m"):
        out["R_a"] = out["surface_D50_m"] / out["buried_D50_m"]
    # legacy, bed_gen2's own kind-based armour ratio
    issurf = t["kind"] == "surface"
    isbulk = t["kind"] == "bulk"
    if (issurf & t["exposed"]).sum() >= 10 and isbulk.sum() >= 10:
        s50 = _wq(D[issurf & t["exposed"]], area[issurf & t["exposed"]], 0.50)
        b50 = _wq(D[isbulk], foot[isbulk], 0.50)
        out["R_a_kind"] = s50 / b50
        out["surface_D50_kind_m"] = s50
        out["buried_D50_kind_m"] = b50
    out["method"] = ("exposed-area weighted D50 of every exposed non-cluster "
                     "grain / footprint weighted D50 of every never-exposed "
                     "non-cluster grain; R_a_kind is bed_gen2's kind-based "
                     "legacy ratio")
    return out


def fines_profile(bed, z, gid, Dc_star):
    """Volume fraction of fines (d < Dc*) in depth bins below the surface.

    The depth unit is the D90 of the FINE fraction itself (the grains that
    percolate / bridge), not the D90 of the framework: the bed is only ~1.6
    D50 deep, so a framework-D90 depth unit would put every grain in the top
    bin.
    """
    t = grain_table(bed, z, gid)
    ok = ~t["cluster"]
    D = t["D"]
    vol = np.array([4.0 / 3.0 * math.pi * g["abc"][0] * g["abc"][1]
                    * g["abc"][2] for g in bed.grains])
    fine = ok & (D < Dc_star)
    if fine.sum() >= 5:
        d90 = float(np.percentile(D[fine], 90))
    else:
        d90 = float(Dc_star / 2.0)
    dep = t["depth_below_surface"] / max(d90, 1e-9)
    lab = ["<0 (above surface)", "0-0.5", "0.5-1", "1-2", "2-3",
           "3-5", "5-8", ">8"]
    out = []
    for i in range(len(lab)):
        if i == 0:
            m = ok & (dep < 0)
        elif i == len(lab) - 1:
            m = ok & (dep > 8)
        else:
            lo, hi = [(0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0),
                      (3.0, 5.0), (5.0, 8.0)][i - 1]
            m = ok & (dep >= lo) & (dep < hi)
        vv = vol[m].sum()
        fv = vol[m & (D < Dc_star)].sum()
        vf = vol[m & (D < Dc_star / 3.0)].sum()
        out.append({"bin": lab[i], "n": int(m.sum()),
                    "solid_volume_m3": float(vv),
                    "fines_volume_fraction": (float(fv / vv) if vv > 0 else None),
                    "very_fine_volume_fraction": (float(vf / vv) if vv > 0 else None)})
    return {"bins": out, "Dc_star_m": float(Dc_star),
            "D90_fine_m": float(d90),
            "depth_unit": "D90 of the fine fraction (d < Dc*)"}


def cluster_coverage(bed, z, gid):
    t = grain_table(bed, z, gid)
    cell2 = bed.cell * bed.cell
    area = t["counts"] * cell2
    m = t["cluster"] & ~(t["kind"] == "mono")
    return float(area[m].sum() / (bed.nx * bed.ny * cell2))


def porosity_bed(bed, z, grid=None):
    """Voxel solid fraction in the top slab (same object as bed_gen.porosity,
    but it also knows the faceted half spaces)."""
    if grid is None:
        if bed.Lx <= 0.61:
            nx, ny, nz = 384, 384, 128
        else:
            nx, ny, nz = 768, 192, 96
    else:
        nx, ny, nz = grid
    slab = 0.05
    dx = bed.Lx / nx
    dy = bed.Ly / ny
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    ztop = BG._bilinear(z, xs[None, :] / bed.cell, ys[:, None] / bed.cell)
    z1 = float(np.median(ztop)) + 1.5 * B2.D50_IN
    z0 = z1 - slab
    dz = slab / (nz - 1)
    zs = z0 + (np.arange(nz) + 0.5) * dz
    occ = np.zeros((nz, ny, nx), dtype=bool)
    for g in bed.grains:
        hx = BG._half_extent(g, 0)
        hy = BG._half_extent(g, 1)
        hz = BG._half_extent(g, 2)
        i0 = max(0, int((g["x"] - hx) / dx))
        i1 = min(nx, int((g["x"] + hx) / dx) + 1)
        j0 = max(0, int((g["y"] - hy) / dy))
        j1 = min(ny, int((g["y"] + hy) / dy) + 1)
        k0 = max(0, int((g["z"] - hz - z0) / dz))
        k1 = min(nz, int((g["z"] + hz - z0) / dz) + 1)
        if i1 <= i0 or j1 <= j0 or k1 <= k0:
            continue
        DX = xs[i0:i1][None, None, :] - g["x"]
        DY = ys[j0:j1][None, :, None] - g["y"]
        DZ = zs[k0:k1][:, None, None] - g["z"]
        if "frac" in g:
            F = _world_faces(g)
            inside = np.ones(np.broadcast(DX, DY, DZ).shape, dtype=bool)
            A, B, C = g["abc"]
            R = g["R"]
            ax = np.array([A, B, C])
            q = np.zeros(inside.shape)
            for m in range(3):
                u = (R[0, m] * DX + R[1, m] * DY + R[2, m] * DZ) / ax[m]
                q += u * u
            inside &= (q <= 1.0)
            for k in range(len(F)):
                inside &= (F[k, 0] * DX + F[k, 1] * DY + F[k, 2] * DZ <= F[k, 3])
            occ[k0:k1, j0:j1, i0:i1] |= inside
        else:
            A, B, C = g["abc"]
            R = g["R"]
            ax = np.array([A, B, C])
            q = np.zeros(np.broadcast(DX, DY, DZ).shape)
            for m in range(3):
                u = (R[0, m] * DX + R[1, m] * DY + R[2, m] * DZ) / ax[m]
                q += u * u
            occ[k0:k1, j0:j1, i0:i1] |= (q <= 1.0)
    ktop = np.clip(((ztop - z0) / dz).astype(np.int64), 1, nz)
    idx = np.arange(nz)[:, None, None] < ktop[None, :, :]
    solid = int((occ & idx).sum())
    total = int(idx.sum())
    return {"porosity": float(1.0 - solid / max(total, 1)),
            "solid_fraction": float(solid / max(total, 1)),
            "slab_m": slab, "z0_m": z0, "z1_m": z1, "grid": [nx, ny, nz]}


def local_mean_map(z, half_px):
    """Periodic box mean (uniform filter with wrap)."""
    from scipy import ndimage
    k = int(max(3, 2 * half_px + 1))
    return ndimage.uniform_filter(z, size=k, mode="wrap")


# ===========================================================================
# A-1 winnow
# ===========================================================================
def tau_grid(bed, tau_star, rng, grid_px=4):
    """tau(x,y) = tau0 + downstream wake deficit of the coarse clasts.

    tau0 = tau* (rho_s-rho) g D50_ref; the wake is a Gaussian strip of length
    Lw in [2,4] D_clast and width ~0.4 D_clast behind every clast with
    D >= D84.  The field exists as a sampled grid (returned) and is sampled
    bilinearly at the grain centres.
    """
    d50 = bed.meta["D50_ref"]
    tau0 = tau_star * DRHO * G * d50
    gx = int(math.ceil(bed.nx / grid_px))
    gy = int(math.ceil(bed.ny / grid_px))
    cellg = bed.cell * grid_px
    M = np.ones((gy, gx))
    xs = (np.arange(gx) + 0.5) * cellg
    ys = (np.arange(gy) + 0.5) * cellg
    D84 = bed.meta.get("D84_ref", 2.0 * d50)
    wakes = []
    large = [g for g in bed.grains if g["D"] >= D84 and g["kind"] != "mono"]
    for g in large:
        Lw = rng.uniform(2.0, 4.0) * g["D"]
        amp = 0.75
        wid = 0.45 * g["D"]
        xs_i = xs - g["x"]
        xs_i = xs_i - bed.Lx * np.round(xs_i / bed.Lx)
        ys_i = ys - g["y"]
        ys_i = ys_i - bed.Ly * np.round(ys_i / bed.Ly)
        m = (xs_i[None, :] > 0) & (xs_i[None, :] < 4.0 * g["D"])
        if not m.any():
            continue
        defc = amp * np.exp(-xs_i[None, :] / Lw) * np.exp(-(ys_i[:, None] / wid) ** 2)
        M = np.where(m, M - defc, M)
        wakes.append({"x": float(g["x"]), "y": float(g["y"]),
                      "D": float(g["D"]), "Lw": float(Lw)})
    M = np.clip(M, 0.05, 1.0)
    tau = tau0 * M
    bed.tau = {"tau_star": float(tau_star), "tau0_Pa": float(tau0),
               "grid": [gy, gx], "cell_m": float(cellg),
               "min_Pa": float(tau.min()), "mean_Pa": float(tau.mean()),
               "max_Pa": float(tau.max()),
               "frac_below_0p8_tau0": float((tau < 0.8 * tau0).mean()),
               "n_wakes": len(wakes), "wakes": wakes[:200]}
    # bilinear sample at grain centres (periodic)
    px = np.mod(np.array([g["x"] for g in bed.grains]), bed.Lx) / cellg - 0.5
    py = np.mod(np.array([g["y"] for g in bed.grains]), bed.Ly) / cellg - 0.5
    i0 = np.floor(px).astype(int)
    j0 = np.floor(py).astype(int)
    fx = px - i0
    fy = py - j0
    i0 %= gx
    j0 %= gy
    i1 = (i0 + 1) % gx
    j1 = (j0 + 1) % gy
    tau_i = ((tau[j0, i0] * (1 - fx) * (1 - fy) + tau[j0, i1] * fx * (1 - fy)
              + tau[j1, i0] * (1 - fx) * fy + tau[j1, i1] * fx * fy))
    return tau, tau_i


def p_crit(tau_Pa, D, phi_deg, N_over, abc):
    """Critical relative exposure of a grain from a moment balance.

    Drag and lift act at the top of the exposed part and at mid-height; the
    submerged weight resists at half the short axis.  With
        u* = sqrt(tau/rho),  u(z) = (u*/kappa) ln(z/z0), z0 = D50/30,
        Fd = 0.5 rho C_D A_f u(D/2)^2,  Fl = 0.5 rho C_L A_f u(D/2)^2
    the balance Fd*(0.5C + p*D) + Fl*0.5C = (W' + N)*0.5B*sin(phi) gives
        p_c = [ (W'+N)*0.5B*sin(phi) - (Fd+Fl)*0.5C ] / (Fd*D)
    which depends on tau*, D, the pivot angle and the overburden (never on
    D_i/D50 alone).
    """
    A, B, C = abc
    tau = max(float(tau_Pa), 1e-9)
    ustar = math.sqrt(tau / RHO_W)
    z0 = Z0_OVER_D50 * max(D, 1e-4)
    zz = max(0.5 * D, z0 * 1.0001)
    u = (ustar / KAPPA) * math.log(zz / z0)
    Af = math.pi / 4.0 * B * C
    Fd = 0.5 * RHO_W * C_D * Af * u * u
    Fl = 0.5 * RHO_W * C_L * Af * u * u
    V = 4.0 / 3.0 * math.pi * A * B * C
    W = DRHO * G * V
    phi = math.radians(phi_deg)
    mom_res = (W + N_over) * 0.5 * B * math.sin(phi)
    mom_drv0 = (Fd + Fl) * 0.5 * C
    den = Fd * D
    if den <= 1e-30:
        return 1.5
    pc = (mom_res - mom_drv0) / den
    return float(np.clip(pc, 0.0, 1.5))


def _exposure(bed, z, gid, t):
    """p_i = (z_top,i - z_local_mean) / D_i on a 3 D_i square.

    z_local_mean is read from a periodic box filter at the nearest of 12
    log-spaced window sizes (within +-12 % of 3 D_i); a per-grain box at 2048^2
    or 8192^2 would cost O(N_grains * window^2).
    """
    from scipy import ndimage
    p = np.zeros(len(bed.grains))
    D = np.array([g["D"] for g in bed.grains])
    wx = np.clip(3.0 * D / bed.cell, 3, min(bed.nx, bed.ny) // 3)
    sizes = np.unique(np.round(np.geomspace(max(3, wx.min()), max(3, wx.max()),
                                            12)).astype(int))
    means = {int(s): ndimage.uniform_filter(z, size=int(s), mode="wrap")
             for s in sizes}
    for i, g in enumerate(bed.grains):
        if t["cluster"][i]:
            p[i] = np.nan
            continue
        s = int(sizes[np.argmin(np.abs(np.log(sizes / wx[i])))])
        cx = int(g["x"] / bed.cell) % bed.nx
        cy = int(g["y"] / bed.cell) % bed.ny
        p[i] = (t["top"][i] - means[s][cy, cx]) / g["D"]
    return p


def _remove(bed, idx):
    return [g for i, g in enumerate(bed.grains) if i not in idx]


def winnow_at(bed, tau_star, rng, max_passes=4, seed_note=""):
    """One tau* realisation of A-1.

    One pass = one transport step.  Eligible grains (exposed, non-cluster,
    p_i > p_c) are ordered by mobility W*_i / D_i (the most mobile fines
    first), and removed until the exposed-area D50 reaches the armour target
    2 x (never-exposed D50).  This bounds the over-shoot: a transport step
    does not remove every unstable grain at once, and winnow is size
    selective.  The bed is re-rasterised once per pass (a later pass picks up
    the grains that the removal exposed).
    """
    b = bed.copy()
    for g in b.grains:
        if "phi_deg" not in g:
            g["phi_deg"] = float(rng.uniform(*PHI_DEG))
    hist = []
    removed = 0
    stop = "max_passes"
    last_p = None
    n_elig = 0
    for it in range(max_passes):
        z, gid = b.render(force=True)
        t = grain_table(b, z, gid)
        m = gsd_metrics(b, z, gid)
        ra = m.get("R_a")
        hist.append({"pass": it, "R_a": ra, "R_a_kind": m.get("R_a_kind"),
                     "n_grains": len(b.grains),
                     "surface_D50_m": m.get("surface_D50_m")})
        if ra is not None and ra >= 2.0:
            stop = "R_a=2"
            break
        buried = m.get("buried_D50_m")
        if buried is None:
            stop = "no_buried_size"
            break
        target = 2.0 * buried
        tau, tau_i = tau_grid(
            b, tau_star, np.random.default_rng(int(rng.integers(1 << 62))),
            grid_px=4)   # 4 bed pixels per tau-grid cell (0.293 mm cell)
        p = _exposure(b, z, gid, t)
        last_p = p
        cell2 = b.cell * b.cell
        area = t["counts"] * cell2
        cand = []
        n_elig = 0
        for i, g in enumerate(b.grains):
            if not t["exposed"][i] or t["cluster"][i]:
                continue
            if not np.isfinite(p[i]):
                continue
            pc = p_crit(tau_i[i], g["D"], g["phi_deg"], 0.0, g["abc"])
            if p[i] <= pc:
                continue
            n_elig += 1
            w = min(1.0, max(0.0, (p[i] - pc) / 0.15)) ** 1.5
            cand.append((w / max(g["D"], 1e-6), w, i))
        if not cand:
            stop = "no_mover"
            break
        cand.sort(reverse=True)
        # exposed set (non-cluster) and its D50 as a function of removals
        exp = [i for i in range(len(b.grains))
               if t["exposed"][i] and not t["cluster"][i]]
        Dv = np.array([b.grains[i]["D"] for i in exp])
        Av = np.array([area[i] for i in exp])
        pos = {gi: k for k, gi in enumerate(exp)}
        drop = np.zeros(len(exp), dtype=bool)
        take = []
        for score, w, i in cand:
            if w < 0.02:
                break
            k = pos.get(i)
            if k is None:
                continue
            drop[k] = True
            take.append(i)
            keep = ~drop
            if keep.sum() < 10:
                break
            d50 = _wq(Dv[keep], Av[keep], 0.50)
            if d50 >= target:
                break
        if not take:
            stop = "no_mover"
            break
        b.grains = _remove(b, set(take))
        removed += len(take)
        b._z = None
        b._gid = None
        if len(take) == 0:
            stop = "no_mover"
            break
    z, gid = b.render(force=True)
    m = gsd_metrics(b, z, gid)
    diag = {"tau_star": float(tau_star), "passes": len(hist),
            "removed_grains": int(removed), "stop_reason": stop,
            "R_a": m.get("R_a"), "R_a_kind": m.get("R_a_kind"),
            "surface_D50_m": m.get("surface_D50_m"),
            "target_D50_m": (2.0 * m["buried_D50_m"]
                             if m.get("buried_D50_m") else None),
            "history": hist,
            "p_stats": ({"median": float(np.nanmedian(last_p)),
                         "p84": float(np.nanpercentile(last_p, 84)),
                         "n_eligible_last_pass": int(n_elig)}
                        if last_p is not None else None),
            "wake_field": b.tau,
            "note": seed_note}
    return b, diag


def winnow(bed, tau_list=TAU_STAR_SWEEP, seed=0):
    """The A-1 sweep.  The selected tau* is the first (ascending) value that
    reaches R_a >= 2; if none does, the highest tau* is kept and the record
    says 'no armour at any tau*' (the brief forbids raising tau* further)."""
    rng = np.random.default_rng(int(seed) + 915)
    rows = []
    selected = None
    for ts in tau_list:
        b, d = winnow_at(bed, ts, rng)
        rows.append(d)
        if selected is None and d["R_a"] is not None and d["R_a"] >= 2.0:
            selected = (ts, b)
    if selected is None:
        selected = (tau_list[-1], b)
    return selected[1], {"sweep": rows, "selected_tau_star": float(selected[0]),
                         "selected_from": "first tau* reaching R_a>=2, else "
                                          "the highest tau* (no armour)"}


# ===========================================================================
# A-2 granular sorting
# ===========================================================================
def granular_sorting(bed, rng, dc_factor=0.20):
    """Size-selective vertical rearrangement below the armour surface.

    Dc* = dc_factor * D50_surface (the brief's 0.15-0.25 D_framework band,
    D_framework = the surface framework D50; the midpoint 0.20 is fixed
    before the run).  In the sheared layer (down to 5 D90 of the fines):
    d < Dc*/3 percolates down (1-4 D90 per pass); Dc*/3 <= d < Dc* bridges
    and stops at 2-5 D90; the coarse is pushed up by squeeze expulsion.
    Solid volume is conserved (the residual is recorded).
    """
    b = bed.copy()
    z, gid = b.render(force=True)
    m = gsd_metrics(b, z, gid)
    D50s = m.get("surface_D50_m") or B2.D50_SURF
    Dc = dc_factor * D50s
    t = grain_table(b, z, gid)
    D = t["D"]
    fine = (~t["cluster"]) & (D < Dc)
    D90 = float(np.percentile(D[fine], 90)) if fine.sum() >= 5 else Dc / 2.0
    vol = np.array([4.0 / 3.0 * math.pi * g["abc"][0] * g["abc"][1]
                    * g["abc"][2] for g in b.grains])
    before = fines_profile(b, z, gid, Dc)
    moved = 0
    dz_sum = 0.0
    vdz = 0.0
    for i, g in enumerate(b.grains):
        if t["cluster"][i]:
            continue
        dep = t["depth_below_surface"][i]
        if dep < -2.0 * D90 or dep > 5.0 * D90:
            continue
        d = g["D"]
        if d < Dc / 3.0:
            step = rng.uniform(1.0, 4.0) * D90
            g["z"] -= step
            moved += 1
            vdz += vol[i] * step
        elif d < Dc:
            target_top = -(rng.uniform(2.0, 5.0) * D90)
            cur_top = -(g["z"] + BG._half_extent(g, 2) - t["z_surface"][i])
            dz = target_top - cur_top
            if dz < 0:
                g["z"] += dz
                moved += 1
                vdz += vol[i] * dz
    # squeeze expulsion: the coarse rises, the moved fines take the same
    # solid volume down (volume conservation in the vertical)
    coarse = [i for i, g in enumerate(b.grains)
              if not t["cluster"][i] and g["D"] >= Dc
              and t["exposed"][i]]
    up = 0.0
    if coarse and vdz < 0:
        target = -vdz
        weights = np.array([vol[i] for i in coarse])
        weights /= weights.sum()
        for k, i in enumerate(coarse):
            dzk = target * weights[k] / vol[i]
            dzk = min(dzk, 0.25 * b.grains[i]["D"])
            b.grains[i]["z"] += dzk
            up += vol[i] * dzk
    b._z = None
    b._gid = None
    z2, gid2 = b.render(force=True)
    after = fines_profile(b, z2, gid2, Dc)
    diag = {"Dc_star_m": float(Dc), "Dc_over_D50_surface": float(dc_factor),
            "D50_surface_m": float(D50s), "D90_fine_m": float(D90),
            "n_moved": int(moved), "n_coarse_lifted": int(len(coarse)),
            "volume_down_m3": float(vdz), "volume_up_m3": float(up),
            "volume_residual_m3": float(vdz + up),
            "profile_before": before, "profile_after": after}
    return b, diag


# ===========================================================================
# A-3 interlock
# ===========================================================================
def interlock(bed, rng, cluster_target=0.10):
    """Imbrication + 3-point support settling + cluster cover to 10 %."""
    b = bed.copy()
    z, gid = b.render(force=True)
    m = gsd_metrics(b, z, gid)
    D84 = m.get("surface_D84_m") or B2.D84_SURF
    D90 = m.get("surface_D90_m") or D84
    t = grain_table(b, z, gid)
    n0 = len(b.grains)
    poro0 = porosity_bed(b, z)
    cov0 = cluster_coverage(b, z, gid)
    # ---- clusters up to the target (added in batches, one re-raster each)
    added = 0
    cov = cov0
    for batch in range(3):
        if cov >= cluster_target:
            break
        area_box = b.Lx * b.Ly
        # area of one group ~ the anchor footprint + 2 mediums + 3 fines
        Dg = float(rng.uniform(1.0, 1.22) * D90)
        a_grp = math.pi * (0.5 * Dg) ** 2 * 2.4
        n_grp = int(math.ceil(max(1.0, (cluster_target - cov) * area_box / a_grp)))
        n_grp = min(n_grp, 24)
        for cg in range(n_grp):
            D = float(rng.uniform(1.0, 1.22) * D90)
            anchor = BG.make_grain(rng, D, kind="cluster_anchor", cluster=99)
            anchor["x"] = float(rng.uniform(0.05, 0.95) * b.Lx)
            anchor["y"] = float(rng.uniform(0.05, 0.95) * b.Ly)
            anchor["p"] = 0.34 * D
            fx = BG._bilinear(b.form,
                              np.array([[anchor["x"] / (b.Lx / b.form.shape[1])]]),
                              np.array([[anchor["y"] / (b.Ly / b.form.shape[0])]]))[0, 0]
            anchor["z"] = float(fx + anchor["p"] - BG._half_extent(anchor, 2))
            grp = [anchor]
            for k in range(2):
                Dm = D * rng.uniform(0.38, 0.55)
                gm = BG.make_grain(rng, Dm, kind="cluster_medium", cluster=99)
                gm["x"] = float((anchor["x"] + (0.55 + 0.32 * k)
                                 * BG._half_extent(anchor, 0)) % b.Lx)
                gm["y"] = float((anchor["y"] + rng.uniform(-0.35, 0.35)
                                 * BG._half_extent(anchor, 1)) % b.Ly)
                gm["z"] = anchor["z"] - 0.35 * BG._half_extent(anchor, 2)
                grp.append(gm)
            for k in range(3):
                Df = D * rng.uniform(0.16, 0.26)
                gf = BG.make_grain(rng, Df, kind="cluster_fine", cluster=99)
                gf["x"] = float((anchor["x"] - (0.55 + 0.42 * k)
                                 * BG._half_extent(anchor, 0)) % b.Lx)
                gf["y"] = float((anchor["y"] + rng.uniform(-0.3, 0.3)
                                 * BG._half_extent(anchor, 1)) % b.Ly)
                gf["z"] = anchor["z"] - 0.6 * BG._half_extent(anchor, 2)
                grp.append(gf)
            b.grains.extend(grp)
            added += len(grp)
        b._z = None
        b._gid = None
        z, gid = b.render(force=True)
        cov = cluster_coverage(b, z, gid)
    # ---- imbrication: re-orient every surface grain ----------------------
    t = grain_table(b, z, gid)
    dips = []
    for i, g in enumerate(b.grains):
        if t["cluster"][i] or not t["exposed"][i]:
            continue
        dip = float(np.clip(rng.normal(28.0, 6.5), 15.0, 40.0))
        az = float(g.get("az_deg", 90.0))
        g["R"] = _rot_from_dip_az(dip, az)
        g["dip_deg"] = dip
        dips.append(dip)
    # ---- 3-point support settling ---------------------------------------
    z, gid = b.render(force=True)
    t = grain_table(b, z, gid)
    P = np.array([[g["x"], g["y"]] for g in b.grains])
    rr = np.array([math.sqrt(BG.footprint_area(g) / math.pi)
                   for g in b.grains])
    hz = np.array([BG._half_extent(g, 2) for g in b.grains])
    from scipy.spatial import cKDTree
    tree = cKDTree(np.mod(P, [b.Lx, b.Ly]), boxsize=[b.Lx, b.Ly])
    supports = np.zeros(len(b.grains), dtype=np.int32)
    order = np.argsort(-t["top"])
    settled = 0
    for i in order:
        if t["cluster"][i] or not t["exposed"][i]:
            continue
        pos = np.mod(P[i], [b.Lx, b.Ly])
        nb = tree.query_ball_point(pos, 2.6 * rr[i])
        base = t["top"][i] - 2.0 * hz[i]
        # matrix support: the local surface is above the grain's base
        loc = t["z_surface"][i]
        mat = 1 if loc > base + 0.1 * b.grains[i]["D"] else 0
        sup = mat
        for j in nb:
            if j == i or t["cluster"][j]:
                continue
            dxy = np.abs(np.mod(np.array([b.grains[j]["x"] - b.grains[i]["x"],
                                          b.grains[j]["y"] - b.grains[i]["y"]]),
                                [b.Lx, b.Ly]))
            dxy = np.minimum(dxy, [b.Lx, b.Ly] - dxy)
            if np.hypot(*dxy) > 1.15 * (rr[i] + rr[j]):
                continue
            # a prop: the neighbour summit lies inside this grain's lower half
            if base - 0.3 * b.grains[i]["D"] < t["top"][j] < base + 0.6 * b.grains[i]["D"]:
                sup += 1
        for _ in range(6):
            if sup >= 3:
                break
            b.grains[i]["z"] -= 0.02 * b.grains[i]["D"]
            settled += 1
            t["top"][i] -= 0.02 * b.grains[i]["D"]
            base = t["top"][i] - 2.0 * hz[i]
            loc = t["z_surface"][i]
            mat = 1 if loc > base + 0.1 * b.grains[i]["D"] else 0
            sup = mat
            for j in nb:
                if j == i or t["cluster"][j]:
                    continue
                dxy = np.abs(np.mod(np.array([b.grains[j]["x"] - b.grains[i]["x"],
                                              b.grains[j]["y"] - b.grains[i]["y"]]),
                                    [b.Lx, b.Ly]))
                dxy = np.minimum(dxy, [b.Lx, b.Ly] - dxy)
                if np.hypot(*dxy) > 1.15 * (rr[i] + rr[j]):
                    continue
                if base - 0.3 * b.grains[i]["D"] < t["top"][j] < base + 0.6 * b.grains[i]["D"]:
                    sup += 1
        supports[i] = sup
    b._z = None
    b._gid = None
    z, gid = b.render(force=True)
    t = grain_table(b, z, gid)
    poro1 = porosity_bed(b, z)
    cov1 = cluster_coverage(b, z, gid)
    surf = np.array([bool(t["exposed"][i] and not t["cluster"][i])
                     for i in range(len(b.grains))])
    sup_s = supports[surf]
    diag = {"n_added_cluster_grains": int(added),
            "cluster_coverage_before": float(cov0),
            "cluster_coverage_after": float(cov1),
            "cluster_target": float(cluster_target),
            "n_surface": int(surf.sum()),
            "support_hist": {str(k): int((sup_s == k).sum())
                             for k in range(0, 7)},
            "support_median": float(np.median(sup_s)) if sup_s.size else None,
            "frac_support_ge3": float((sup_s >= 3).mean()) if sup_s.size else None,
            "n_sink_steps": int(settled),
            "dip_median_deg": float(np.median(dips)) if dips else None,
            "dip_p16_deg": float(np.percentile(dips, 16)) if dips else None,
            "dip_p84_deg": float(np.percentile(dips, 84)) if dips else None,
            "n_dip": int(len(dips)),
            "porosity_before": poro0, "porosity_after": poro1,
            "porosity_drop": float(poro0["porosity"] - poro1["porosity"])}
    return b, diag


# ===========================================================================
# B-1 faceted grains
# ===========================================================================
def facetise(bed, rng):
    """Replace every grain's smooth ellipsoid by a fractured ellipsoid."""
    b = bed.copy()
    stats = []
    rejected = 0
    for g in b.grains:
        if "frac" in g:
            continue
        chosen = None
        f = None
        for attempt in range(4):
            f = make_facet_faces(rng, g["abc"])
            if f is None:
                rejected += 1
                continue
            mm = f["metrics"]
            if mm["n_faces"] < 6:
                rejected += 1
                continue
            chosen = f
            break
        if chosen is None:
            chosen = f
        g["frac"] = chosen
        stats.append(chosen["metrics"])
    b._z = None
    b._gid = None
    all_lin = [v for s in stats for v in s["facet_lin"]]
    nfaces = np.array([s["n_faces"] for s in stats])
    ba = np.array([s["calipers"][1] / s["calipers"][0] for s in stats])
    corey = np.array([s["calipers"][2] / math.sqrt(s["calipers"][0]
                                                   * s["calipers"][1])
                      for s in stats])
    diag = {"n_grains": len(stats), "rejected_draws": int(rejected),
            "volume_over_ellipsoid": {
                "median": float(np.median([
                    s["volume"] / (4.0 / 3.0 * math.pi
                                   * g["abc"][0] * g["abc"][1] * g["abc"][2])
                    for s, g in zip(stats, [x for x in b.grains if "frac" in x])])),
                "definition": "solid volume of the fractured grain / smooth "
                              "ellipsoid volume at the same calipers"},
            "n_faces": {"min": int(nfaces.min()), "median": float(np.median(nfaces)),
                        "max": int(nfaces.max()),
                        "hist": {str(k): int((nfaces == k).sum())
                                 for k in range(6, 13)}},
            "facet_linear_m": {"n": len(all_lin),
                               "median": float(np.median(all_lin)),
                               "p16": float(np.percentile(all_lin, 16)),
                               "p84": float(np.percentile(all_lin, 84)),
                               "min": float(np.min(all_lin)),
                               "max": float(np.max(all_lin)),
                               "definition": "sqrt(face area)"},
            "facet_linear_mm": {"median": float(np.median(all_lin) * 1000),
                                "p16": float(np.percentile(all_lin, 16) * 1000),
                                "p84": float(np.percentile(all_lin, 84) * 1000)},
            "ba": {"median": float(np.median(ba)), "p16": float(np.percentile(ba, 16)),
                   "p84": float(np.percentile(ba, 84))},
            "corey": {"median": float(np.median(corey)),
                      "p16": float(np.percentile(corey, 16)),
                      "p84": float(np.percentile(corey, 84))},
            "wadell_volume_proxy": {"median": float(np.median([s["wadell_volume_proxy"] for s in stats])),
                                    "p16": float(np.percentile([s["wadell_volume_proxy"] for s in stats], 16)),
                                    "p84": float(np.percentile([s["wadell_volume_proxy"] for s in stats], 84)),
                                    "definition": "(6V/pi)^(1/3) / (2 R_circum)"},
            "facet_edge_roundness": {"median": float(np.median([s["facet_edge_roundness"] for s in stats])),
                                     "p16": float(np.percentile([s["facet_edge_roundness"] for s in stats], 16)),
                                     "p84": float(np.percentile([s["facet_edge_roundness"] for s in stats], 84)),
                                     "definition": "0.6*clip(1 - mean(dihedral)/90deg,0,1); "
                                                   "0 = angular knife edge, 0.6 = flush (Powers scale)"},
            "facet_dihedral_deg": {"median": float(np.median([s["facet_dihedral_deg"] for s in stats])),
                                   "p16": float(np.percentile([s["facet_dihedral_deg"] for s in stats], 16)),
                                   "p84": float(np.percentile([s["facet_dihedral_deg"] for s in stats], 84))},
            "facet_area_fraction": {"median": float(np.median([s["facet_area_fraction"] for s in stats])),
                                    "p16": float(np.percentile([s["facet_area_fraction"] for s in stats], 16)),
                                    "p84": float(np.percentile([s["facet_area_fraction"] for s in stats], 84)),
                                    "definition": "total flat-facet area / smooth ellipsoid area"},
            "volume_preserved": {"median": float(np.median([s["volume"] for s in stats])),
                                 "ellipsoid_reference": None,
                                 "definition": "Voxel/lattice volume of the fractured grain over "
                                               "(4/3) pi A B C; the reference is filled by the driver"}}
    return b, diag


def facet_stats(bed, z=None, gid=None, exposed_only=False):
    """Facet size / shape statistics over all grains or over the exposed set."""
    polys = []
    areas = []
    vols = []
    vells = []
    t = grain_table(bed, z, gid) if (exposed_only and z is not None) else None
    for i, g in enumerate(bed.grains):
        if "frac" not in g:
            continue
        if exposed_only and t is not None and not t["exposed"][i]:
            continue
        m = g["frac"]["metrics"]
        polys.extend(m["facet_lin"])
        areas.extend(m["facet_areas"])
        vols.append(m["volume"])
        vells.append(4.0 / 3.0 * math.pi * g["abc"][0] * g["abc"][1] * g["abc"][2])
    if not polys:
        return None
    lin = np.array(polys)
    return {"n_facets": int(len(lin)),
            "median_m": float(np.median(lin)),
            "p16_m": float(np.percentile(lin, 16)),
            "p84_m": float(np.percentile(lin, 84)),
            "median_mm": float(np.median(lin) * 1000),
            "area_weighted_median_mm": float(_wq(lin, np.array(areas),
                                                 0.5) * 1000),
            "volume_over_ellipsoid": float(np.mean(np.array(vols) / np.array(vells)))}


# ===========================================================================
# B-2 planation
# ===========================================================================
def planation_surface(bed, z, k=0.5):
    """Local mean + k * local sigma over a 3 D50 square."""
    from scipy import ndimage
    half = int(max(2, round(1.5 * bed.meta["D50_ref"] / bed.cell)))
    size = 2 * half + 1
    mu = ndimage.uniform_filter(z, size=size, mode="wrap")
    mu2 = ndimage.uniform_filter(z * z, size=size, mode="wrap")
    sd = np.sqrt(np.maximum(mu2 - mu * mu, 0.0))
    return mu + k * sd, {"k": float(k), "window_m": size * bed.cell,
                         "mean_m": float(mu.mean()), "sd_m": float(sd.mean())}


def ndimage_uniform(a, size):
    from scipy import ndimage
    return ndimage.uniform_filter(a, size=int(size), mode="wrap")


def planation(bed, k=0.5, mode="clip", rng=None, top_f=None):
    """B-2a ('clip') cuts the heads; B-2b ('remove') deletes the protruding
    grains and lets the bed below show.

    top_f (optional) is the equivalent quantile reading: the surface is the
    local (1-f) quantile of the summit heights, reported as a cross-check.
    """
    b = bed.copy()
    z, gid = b.render(force=True)
    P, meta = planation_surface(b, z, k=k)
    if top_f is not None:
        from scipy.stats import norm
        half = int(max(2, round(1.5 * b.meta["D50_ref"] / b.cell)))
        size = 2 * half + 1
        mu = ndimage_uniform(z, size)
        sd = np.sqrt(np.maximum(ndimage_uniform(z * z, size) - mu * mu, 0.0))
        Pq = mu + float(norm.ppf(1.0 - top_f)) * sd
    else:
        Pq = None
    above = z > P
    f_meas = float(above.mean())
    V_cut = float(np.maximum(z - P, 0.0).sum() * b.cell * b.cell)
    diag = {"mode": mode, "k": float(k), "window_m": meta["window_m"],
            "clip_area_fraction": f_meas,
            "equiv_top_f": (float(top_f) if top_f is not None else None),
            "equiv_top_f_measured": (float((z > Pq).mean()) if Pq is not None else None),
            "cut_volume_m3": V_cut,
            "surface_mean_m": float(P.mean()),
            "surface_min_m": float(P.min()), "surface_max_m": float(P.max())}
    if mode == "clip":
        b.clip = P
        b._z = None
        b._gid = None
        diag["note"] = "grain heads cut at the planation surface; grains stay; " \
                       "the cut volume leaves the system"
        # per-grain clip amount for the record
        t = grain_table(b, z, gid)
        clipped = []
        for i, g in enumerate(b.grains):
            if t["cluster"][i]:
                continue
            zi = int(np.clip(g["x"] / b.cell, 0, b.nx - 1))
            yi = int(np.clip(g["y"] / b.cell, 0, b.ny - 1))
            d = t["top"][i] - P[yi, zi]
            if d > 0:
                clipped.append(d)
        diag["n_clipped_grains"] = int(len(clipped))
        diag["clipped_median_m"] = float(np.median(clipped)) if clipped else 0.0
        diag["clipped_p84_m"] = float(np.percentile(clipped, 84)) if clipped else 0.0
    else:
        t = grain_table(b, z, gid)
        yidx = np.clip((np.array([g["y"] for g in b.grains]) / b.cell).astype(int), 0, b.ny - 1)
        xidx = np.clip((np.array([g["x"] for g in b.grains]) / b.cell).astype(int), 0, b.nx - 1)
        excess = t["top"] - P[yidx, xidx]
        cand = [i for i in range(len(b.grains))
                if not t["cluster"][i] and t["exposed"][i] and excess[i] > 0]
        cand.sort(key=lambda i: -excess[i])
        # B-2a and B-2b are the same cut: the removed grains cover the same
        # surface fraction f that B-2a clips.  Removing them by descending
        # excess until the exposed area removed reaches f * A keeps the pair
        # comparable and stops the removal from emptying the armour layer.
        cell2 = b.cell * b.cell
        a_target = f_meas * (b.nx * b.ny * cell2)
        rm = []
        a_cum = 0.0
        for i in cand:
            if a_cum >= a_target:
                break
            rm.append(i)
            a_cum += t["counts"][i] * cell2
        diag["n_removed_grains"] = int(len(rm))
        diag["removed_area_fraction"] = float(a_cum / (b.nx * b.ny * cell2))
        diag["removed_excess_median_m"] = (float(np.median(excess[rm])) if rm else 0.0)
        seats = [(b.grains[i]["x"], b.grains[i]["y"], b.grains[i]["D"],
                  P[yidx[i], xidx[i]]) for i in rm]
        b.grains = _remove(b, set(rm))
        # the seat is filled by the grains directly below (the brief: "その場所
        # を下の粒／ポケット細粒で埋める").  Up to three neighbours below each
        # seat are raised toward the planation surface; the raised volume is
        # recorded, and the removed grains' volume leaves the system.
        from scipy.spatial import cKDTree
        if b.grains:
            Pg = np.array([[g["x"], g["y"]] for g in b.grains])
            tree = cKDTree(np.mod(Pg, [b.Lx, b.Ly]), boxsize=[b.Lx, b.Ly])
            raised = 0
            vol_up = 0.0
            t2 = grain_table(b, b.render(force=True)[0], b.render()[1])
            for (sx, sy, sD, sP) in seats:
                nb = tree.query_ball_point(np.mod([sx, sy], [b.Lx, b.Ly]),
                                           1.5 * sD)
                cand = [j for j in nb
                        if not t2["cluster"][j]
                        and sP - 2.5 * sD < t2["top"][j] < sP + 0.05 * sD]
                cand.sort(key=lambda j: -t2["top"][j])
                for j in cand[:3]:
                    want = sP - 0.25 * sD
                    dz = min(want - t2["top"][j], 2.0 * b.grains[j]["D"])
                    if dz > 0:
                        b.grains[j]["z"] += dz
                        t2["top"][j] += dz
                        raised += 1
                        vol_up += (4.0 / 3.0 * math.pi
                                   * b.grains[j]["abc"][0]
                                   * b.grains[j]["abc"][1]
                                   * b.grains[j]["abc"][2]) * dz
            diag["seats"] = len(seats)
            diag["n_raised"] = int(raised)
            diag["raised_volume_m3"] = float(vol_up)
        b._z = None
        b._gid = None
        diag["note"] = ("protruding grains removed (volume leaves the system); "
                        "the seat is filled by raising up to three grains "
                        "directly below toward the planation surface; the "
                        "remaining pore space is B-3's")
    return b, diag


# ===========================================================================
# B-3 pocket fill
# ===========================================================================
def pocket_fill(bed, rng, Dc_star, fill_fraction=0.5, max_pores=6000):
    """Fines (d < Dc*/3) scattered into the inter-grain depressions.

    Candidates are the pixels where the surface lies below the planation
    surface by more than max(2 mm, 0.03 D50); they are taken deepest-first
    with a spacing so that at most `max_pores` grains are added (a declared
    computational cap, reported).  Each fine is seated at 60 % of the local
    depression depth.  The porosity is measured before and after.
    """
    b = bed.copy()
    z, gid = b.render(force=True)
    poro0 = porosity_bed(b, z)
    if b.clip is not None:
        P = b.clip
    else:
        P, _ = planation_surface(b, z, k=1.0)
    dep = np.maximum(P - z, 0.0)
    m = dep > max(0.002, 0.03 * b.meta["D50_ref"])
    if not m.any():
        return b, {"n_filled": 0, "filled_volume_m3": 0.0,
                   "depression_area_fraction": 0.0, "fill_fraction": fill_fraction,
                   "porosity_before": poro0, "porosity_after": poro0}
    target_vol = fill_fraction * float(dep[m].sum() * b.cell * b.cell)
    idx = np.argwhere(m)
    order = np.argsort(-dep[m])
    idx = idx[order]
    span = max(1, int(math.ceil(len(idx) / float(max_pores))))
    vol_added = 0.0
    n = 0
    for k in range(0, len(idx), span):
        yy, xx = int(idx[k][0]), int(idx[k][1])
        depth = float(dep[yy, xx])
        D = float(rng.uniform(0.2, 0.9) * Dc_star)
        g = BG.make_grain(rng, D, kind="pocket_fine", cluster=None)
        g["x"] = float((xx + 0.5 + rng.uniform(-0.5, 0.5)) * b.cell) % b.Lx
        g["y"] = float((yy + 0.5 + rng.uniform(-0.5, 0.5)) * b.cell) % b.Ly
        g["p"] = 0.0
        g["z"] = float(z[yy, xx] + 0.60 * depth - BG._half_extent(g, 2))
        b.grains.append(g)
        vol_added += 4.0 / 3.0 * math.pi * g["abc"][0] * g["abc"][1] * g["abc"][2]
        n += 1
        if vol_added >= target_vol:
            break
    b._z = None
    b._gid = None
    z2, gid2 = b.render(force=True)
    poro1 = porosity_bed(b, z2)
    diag = {"n_filled": int(n), "filled_volume_m3": float(vol_added),
            "target_volume_m3": float(target_vol),
            "fill_of_target": float(vol_added / max(target_vol, 1e-30)),
            "pore_cap": int(max_pores), "stride": int(span),
            "depression_area_fraction": float(m.mean()),
            "fill_fraction": float(fill_fraction),
            "Dc_star_m": float(Dc_star),
            "D_range_m": [float(Dc_star * 0.2), float(Dc_star * 0.9)],
            "porosity_before": poro0, "porosity_after": poro1}
    return b, diag
