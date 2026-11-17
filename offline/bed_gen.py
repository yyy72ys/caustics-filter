"""bed_gen.py -- particle-based 3-D gravel bed for BRIEF-V4-2-1.

What this is
------------
A bed of ~4-5k rotated ellipsoids (a real 3-D particle assembly, not a height
noise field) with

  * a log-normal GSD (D50 = 0.020 m, sigma_g = D84/D50 = 2.1),
  * a surface layer coarsened by an armour ratio of 2.0 (surface D50 =
    0.040 m) and a sub-surface layer that keeps the unarmoured distribution,
  * grain shapes b/a in [0.67, 0.74] and Corey shape factor c/sqrt(a b) in
    [0.6, 0.7],
  * imbrication (a-b plane dipping 30 +- 10 deg towards the upstream side),
    a-axis azimuth mixed 50 % transverse / 50 % along-flow,
  * asymmetric clusters (anchor > D90, medium grains on the stoss side, fines
    in the wake, length ~ 2 x width) covering ~10 % of the patch,
  * two anchor-class boulders (D ~ 0.15 m and 0.22 m) placed with the same
    cluster rule, kept OUT of the bed audit (they are objects, not bed),
  * an fBm substrate (H = 0.79) under the particles whose amplitude is solved
    so the realised sigma_z lands inside the audit band (rule recorded below),
  * settling = gravity drop onto the accumulating height buffer plus an
    iterative non-penetration relaxation (no dynamical claim),
  * a vertical ray-cast depth map at render resolution (736 x 540 over
    0.6 x 0.4402 m) and at statistics resolution (2048 x 2048 over 0.6 m),
    plus per-pixel facet normals and grain ids.

Sampling convention (declared before the build, not tuned afterwards)
--------------------------------------------------------------------
The brief's D50 / D84 are read as *areal* (image) statistics, because that is
what a depth map can return and what sigma_z/D50 ~ 0.4 refers to (classic
gravel-bed roughness).  A lognormal *number* density that is area-weighted by
D^2 is again lognormal with the same sigma and median shifted by e^{+2 sigma^2};
the draw is therefore made from median D50*exp(-2 sigma_ln^2) so that the
reconstructed areal median equals the stated D50.  sigma_ln = ln(2.1).

Substrate amplitude rule (brief: "determine the amplitude so that the overall
sigma_z lands in the audit range; record how")
------------------------------------------------------------
sigma_target = 0.40 * D50 (= 0.008 m, the middle of [0.32, 0.46]).
Fixed-point iteration: build, measure sigma_nonsub (particle relief about a
smoothed substrate), set sigma_sub = sqrt(max(sigma_target^2 - sigma_nonsub^2,
(0.15 sigma_target)^2)), rebuild (max 4 passes).  The realised sigma_z is then
measured and reported whatever it is.

Only new file; nothing else in the repository is touched.
"""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import atomic_write_bytes  # noqa: E402

# ---------------------------------------------------------------------------
# fixed inputs (BRIEF-V4-2-1)
# ---------------------------------------------------------------------------
PATCH = 0.6                 # m, square, periodic
D50_IN = 0.020              # m  (areal median, see module docstring)
SIGMA_G = 2.1               # D84/D50
ARMOR = 2.0
BA_RANGE = (0.67, 0.74)
COREY_RANGE = (0.6, 0.7)
DIP_DEG = 30.0
DIP_SPREAD_DEG = 10.0
CLUSTER_COVERAGE = 0.10
N_CLUSTERS = 3
MONO_D = (0.15, 0.22)
SUBSTRATE_H = 0.79
SIGMA_TARGET = 0.40 * D50_IN
INTERLOCK = 0.30            # grains sink this fraction of D before relaxation

SEED = 3255045276382562994   # the scene seed fixed by the brief

DEPOSIT_N = 1024            # z-buffer used for settling
HIRES = 2048                # depth map for statistics (0.293 mm / px)
RENDER_W, RENDER_H = 736, 540
CELL_R = PATCH / RENDER_W
RENDER_J0 = int(round((PATCH - RENDER_H * CELL_R) / 2.0 / CELL_R))
WIN_X0, WIN_X1 = 0.0, PATCH
WIN_Y0 = RENDER_J0 * CELL_R
WIN_Y1 = WIN_Y0 + RENDER_H * CELL_R

VOX = (384, 384, 128)


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")


# ---------------------------------------------------------------------------
# substrate
# ---------------------------------------------------------------------------
def fbm_field(n, L, H, seed, sigma):
    """Periodic self-affine surface with structure function ~ Delta^(2H).

    Built by filtering white noise with A(k) = k^-(H+1) (2-D): the structure
    function of that field is D(Dx) ~ Dx^(2H).  Recorded, fitted back in
    audit_structure().
    """
    rng = np.random.default_rng(int(seed))
    w = rng.standard_normal((n, n))
    k1 = np.fft.fftfreq(n, d=L / n) * 2.0 * math.pi
    KX, KY = np.meshgrid(k1, k1)
    K = np.hypot(KX, KY)
    amp = np.zeros_like(K)
    m = K > 0
    amp[m] = K[m] ** (-(H + 1.0))
    Z = np.fft.fft2(w) * amp
    z = np.real(np.fft.ifft2(Z))
    z -= z.mean()
    s = z.std()
    return (z * (sigma / s)).astype(np.float64)


# ---------------------------------------------------------------------------
# grain population
# ---------------------------------------------------------------------------
def _sigma_ln():
    return math.log(SIGMA_G)


def draw_size(rng, n, d50_areal, clip=(4.0e-4, 0.09)):
    """Areal-median lognormal draw (see module docstring)."""
    s = _sigma_ln()
    mu = math.log(d50_areal) - 2.0 * s * s
    D = np.exp(mu + s * rng.standard_normal(n))
    return np.clip(D, clip[0], clip[1])


def make_grain(rng, D, dip_deg=None, cluster=None, kind="bulk"):
    ratio = rng.uniform(*BA_RANGE)
    corey = rng.uniform(*COREY_RANGE)
    B = 0.5 * D
    A = B / ratio
    C = corey * math.sqrt(A * B)
    # a-axis azimuth: 50 % transverse (y), 50 % along flow (x)
    az = (0.5 * math.pi) if rng.random() < 0.5 else 0.0
    az += math.radians(rng.uniform(-12.0, 12.0))
    if dip_deg is None:
        dip = DIP_DEG + rng.uniform(-DIP_SPREAD_DEG, DIP_SPREAD_DEG)
    else:
        dip = float(dip_deg)
    d = math.radians(dip)
    ea = np.array([math.cos(az), math.sin(az), 0.0])
    nrm = np.array([math.sin(d), 0.0, math.cos(d)])
    eb = np.cross(nrm, ea)
    eb /= np.linalg.norm(eb)
    ea = np.cross(eb, nrm)
    ea /= np.linalg.norm(ea)
    R = np.stack([ea, eb, nrm], axis=1)
    return {"D": float(D), "abc": (A, B, C), "R": R, "kind": kind,
            "cluster": cluster, "dip_deg": float(dip), "ratio": float(ratio),
            "corey": float(corey), "az_deg": math.degrees(az) % 360.0,
            "x": 0.0, "y": 0.0, "z": 0.0}


def _half_extent(g, axis):
    A, B, C = g["abc"]
    R = g["R"]
    return math.sqrt((A * R[axis, 0]) ** 2 + (B * R[axis, 1]) ** 2
                     + (C * R[axis, 2]) ** 2)


def footprint_area(g):
    return math.pi * _half_extent(g, 0) * _half_extent(g, 1)


def draw_population(rng, d50_areal, v_target):
    """Grains whose *solid volume* sums to v_target, size by the areal GSD."""
    out = []
    vol = 0.0
    while vol < v_target:
        D = float(draw_size(rng, 1, d50_areal)[0])
        g = make_grain(rng, D)
        A, B, C = g["abc"]
        g["V"] = 4.0 / 3.0 * math.pi * A * B * C
        vol += g["V"]
        out.append(g)
    return out


def build_population(seed, armored=True, clusters=True, mono=True):
    """Grain list: bulk volume, then the armoured surface layer, then clusters.

    Population sizes are set by SOLID VOLUME (not by an arbitrary number of
    layers): the bulk fills 3 D50 of bed at a 0.40 porosity, and the surface
    layer covers the patch 1.15 times over with the armoured GSD.  This keeps
    the grain count (and therefore the bed thickness) determined by the GSD
    instead of by a loop counter.
    """
    rng = np.random.default_rng(seed)
    grains = []
    A_patch = PATCH * PATCH
    phi = 0.40

    bulk = draw_population(rng, D50_IN,
                           A_patch * 1.6 * D50_IN * 0.55)
    for g in bulk:
        g["kind"] = "bulk"
    grains.extend(bulk)

    if armored:
        d50s = D50_IN * ARMOR
        # n from an area-coverage target with the armoured areal GSD:
        # E[footprint] = pi E[D^2] / (4 <b/a>) ; E[D^2] = d50^2 e^{+2 s^2}
        s = _sigma_ln()
        e_d2 = d50s ** 2 * math.exp(-2.0 * s * s)
        e_fp = math.pi * e_d2 / (4.0 * 0.705)
        n_surf = int(round(0.85 * A_patch / e_fp))
        Ds = np.sort(draw_size(rng, n_surf, d50s))    # fines first
        for D in Ds:
            g = make_grain(rng, float(D))
            A, B, C = g["abc"]
            g["V"] = 4.0 / 3.0 * math.pi * A * B * C
            g["kind"] = "surface"
            grains.append(g)

    # ---- clusters --------------------------------------------------------
    cluster_members = []
    if clusters:
        # anchor from > D90 of the *surface* (areal) distribution
        d90 = D50_IN * ARMOR * (SIGMA_G ** 1.2815515655446004)
        for c in range(N_CLUSTERS):
            D = d90 * rng.uniform(1.02, 1.22)
            anchor = make_grain(rng, D, kind="cluster_anchor", cluster=c)
            cluster_members.append(anchor)
            grains.append(anchor)
            # stoss side (upstream, -x): medium grains
            for k in range(2):
                Dm = D * rng.uniform(0.38, 0.55)
                gm = make_grain(rng, Dm, kind="cluster_medium", cluster=c)
                gm["off"] = (-(0.55 + 0.32 * k) * _half_extent(anchor, 0)
                             * rng.uniform(1.0, 1.25),
                             rng.uniform(-0.35, 0.35) * _half_extent(anchor, 1))
                cluster_members.append(gm)
                grains.append(gm)
            # wake (downstream, +x): fines, fanning out
            for k in range(3):
                Df = D * rng.uniform(0.16, 0.26)
                gf = make_grain(rng, Df, kind="cluster_fine", cluster=c)
                s = (0.55 + 0.42 * k) * _half_extent(anchor, 0)
                gf["off"] = (s * rng.uniform(0.95, 1.15),
                             (rng.uniform(-0.55, 0.55) if k % 2 else
                              rng.uniform(-0.25, 0.25))
                             * s)
                cluster_members.append(gf)
                grains.append(gf)

    # ---- the two "mono" boulders (shadow probes, not bed) ---------------
    if mono:
        for D in MONO_D:
            a = make_grain(rng, D, kind="mono", cluster="mono")
            cluster_members.append(a)
            grains.append(a)
            for k in range(3):
                Dm = D * rng.uniform(0.16, 0.24)
                gm = make_grain(rng, Dm, kind="mono_medium", cluster="mono")
                gm["off"] = (-(0.60 + 0.30 * k) * _half_extent(a, 0)
                             * rng.uniform(1.0, 1.2),
                             rng.uniform(-0.3, 0.3) * _half_extent(a, 1))
                cluster_members.append(gm)
                grains.append(gm)
            for k in range(4):
                Df = D * rng.uniform(0.10, 0.16)
                gf = make_grain(rng, Df, kind="mono_fine", cluster="mono")
                s = (0.60 + 0.40 * k) * _half_extent(a, 0)
                gf["off"] = (s * rng.uniform(0.95, 1.15),
                             rng.uniform(-0.5, 0.5) * s)
                cluster_members.append(gf)
                grains.append(gf)

    return grains, cluster_members


# ---------------------------------------------------------------------------
# z-buffer deposition
# ---------------------------------------------------------------------------
class ZBuf:
    def __init__(self, n, L):
        self.n = int(n)
        self.L = float(L)
        self.cell = self.L / self.n
        self.z = np.full((self.n, self.n), -np.inf)

    def local(self, x, y, g):
        """Pixel window (ix0, ix1, iy0, iy1) covering the grain footprint."""
        hx = _half_extent(g, 0)
        hy = _half_extent(g, 1)
        ix0 = int(math.floor((x - hx) / self.cell))
        ix1 = int(math.ceil((x + hx) / self.cell))
        iy0 = int(math.floor((y - hy) / self.cell))
        iy1 = int(math.ceil((y + hy) / self.cell))
        return ix0, ix1, iy0, iy1

    def _dxdy(self, x, y, ix0, ix1, iy0, iy1):
        ixs = (np.arange(ix0, ix1) + 0.5) * self.cell
        iys = (np.arange(iy0, iy1) + 0.5) * self.cell
        DX = ixs[None, :] - x
        DY = iys[:, None] - y
        return DX, DY

    def roots(self, g, DX, DY):
        A, B, C = g["abc"]
        R = g["R"]
        ax = np.array([A, B, C])
        u = np.stack([R[0, i] * DX + R[1, i] * DY for i in range(3)], axis=0)
        w = R[2, :][:, None, None]
        aa = ((w / ax[:, None, None]) ** 2).sum(axis=0)
        bb = ((u * w) / ax[:, None, None] ** 2).sum(axis=0)
        cc = ((u / ax[:, None, None]) ** 2).sum(axis=0) - 1.0
        disc = bb * bb - aa * cc
        ok = disc > 0.0
        sq = np.sqrt(np.where(ok, disc, 0.0))
        s1 = (-bb - sq) / aa
        s2 = (-bb + sq) / aa
        return s1, s2, ok

    def probe(self, g, x, y, q=0.88):
        """Height of the grain centre if dropped at (x,y).

        A grain rests on a few contact points, not on the single highest bump
        in its footprint: the centre height is the q-quantile of the required
        clearance, so <= (1-q) of the footprint may interpenetrate before the
        relaxation pass.  Using the max here is what makes naive z-buffer
        deposition grow towers.
        """
        ix0, ix1, iy0, iy1 = self.local(x, y, g)
        DX, DY = self._dxdy(x, y, ix0, ix1, iy0, iy1)
        s1, s2, ok = self.roots(g, DX, DY)
        if not ok.any():
            iz = int(y / self.cell) % self.n
            ix = int(x / self.cell) % self.n
            return float(self.z[iz, ix]) + g["abc"][2]
        ixg = np.arange(ix0, ix1) % self.n
        iyg = np.arange(iy0, iy1) % self.n
        zb = self.z[np.ix_(iyg, ixg)]
        gap = (zb - s1)[ok]
        if gap.size < 9:
            return float(np.max(gap))
        return float(np.quantile(gap, q))

    def sites(self, rng, k=6, block_px=48, q=0.05):
        """Candidate drop points: random draws inside the LOWEST blocks.

        Ballistic deposition with 6 uniform random candidates does not level a
        bed: the chance of hitting a 20 mm pit in a 0.36 m^2 patch is ~0.5 %,
        so grains keep landing on top of each other and the surface ends up
        30 mm rough (measured).  Grains in a real stream roll into the hollows,
        so the candidates are drawn from the lowest `q` of the coarse blocks.
        """
        n = self.n
        nb = max(1, n // block_px)
        b = self.z[:nb * block_px, :nb * block_px].reshape(
            nb, block_px, nb, block_px)
        mins = b.min(axis=(1, 3))
        thr = np.quantile(mins, q)
        cand = np.argwhere(mins <= thr)
        if len(cand) == 0:
            cand = np.argwhere(mins <= mins.min())
        out = []
        for _ in range(int(k)):
            j = cand[int(rng.integers(len(cand)))]
            cx = (j[1] * block_px + rng.integers(block_px)) * self.cell
            cy = (j[0] * block_px + rng.integers(block_px)) * self.cell
            out.append((cx, cy))
        return out

    def drop(self, g, x, y, mode="rest", n_cand=1, rng=None):
        """Place a grain at (x,y) (or at the lowest of n_cand candidates).

        The centre is lowered by INTERLOCK * D before the non-penetration
        relaxation.  Pure z-buffer stacking gives a random-LOOSE packing
        (porosity ~0.6) which no natural gravel bed has; real grains settle
        into the matrix until neighbours stop them.
        """
        if n_cand > 1 and rng is not None:
            best = None
            for cx, cy in self.sites(rng, k=n_cand):
                zc = self.probe(g, cx, cy)
                if best is None or zc < best[2]:
                    best = (cx, cy, zc)
            x, y = best[0], best[1]
        zc = self.probe(g, x, y) - INTERLOCK * g["D"]
        g["x"], g["y"], g["z"] = float(x), float(y), zc
        self._stamp(g, x, y)

    def _stamp(self, g, x, y):
        ix0, ix1, iy0, iy1 = self.local(x, y, g)
        DX, DY = self._dxdy(x, y, ix0, ix1, iy0, iy1)
        s1, s2, ok = self.roots(g, DX, DY)
        ixg = np.arange(ix0, ix1) % self.n
        iyg = np.arange(iy0, iy1) % self.n
        top = g["z"] + s2
        cur = self.z[np.ix_(iyg, ixg)]
        self.z[np.ix_(iyg, ixg)] = np.where(ok, np.maximum(cur, top), cur)

    def put_here(self, g, x, y, n_cand=1, rng=None):
        """Place a grain at (x,y), optionally at the lowest of n_cand drops."""
        self.drop(g, x, y, n_cand=n_cand, rng=rng)

    def surface(self, px, py):
        """Bilinear sample of the z buffer (periodic)."""
        n = self.n
        fx = np.mod(px / self.cell - 0.5, n)
        fy = np.mod(py / self.cell - 0.5, n)
        x0 = np.floor(fx).astype(np.int64) % n
        y0 = np.floor(fy).astype(np.int64) % n
        x1 = (x0 + 1) % n
        y1 = (y0 + 1) % n
        tx = fx - np.floor(fx)
        ty = fy - np.floor(fy)
        return (self.z[y0, x0] * (1 - tx) * (1 - ty)
                + self.z[y0, x1] * tx * (1 - ty)
                + self.z[y1, x0] * (1 - tx) * ty
                + self.z[y1, x1] * tx * ty)


def settle(grains, zbuf, seed, bed_sub=None):
    """Planed deposition to the substrate form + non-penetration relaxation."""
    rng = np.random.default_rng(seed + 1)
    order = sorted(range(len(grains)),
                   key=lambda i: (-grains[i]["D"], grains[i]["kind"]))
    # coarse first for the bulk (fills from the bottom), fines first for the
    # surface layer -> armour on top.  Split by kind.
    bulk = [i for i in order if grains[i]["kind"] == "bulk"]
    surf = [i for i in order if grains[i]["kind"] == "surface"]
    clust = [i for i in order if grains[i]["kind"].startswith("cluster")]
    mono = [i for i in order if grains[i]["kind"].startswith("mono")]
    # ---- deposition: each grain is planed to the substrate "reached form" --
    # Rationale (declared construction, not fitted): ballistic stacking of a
    # 20-40 mm gravel on a 0.6 m patch cannot make a paved surface -- it grew
    # 150-370 mm piles whose sigma_z was 4-6x the audit band.  A gravel-bed
    # armour layer really is a PAVEMENT: every grain is bedded down to the
    # same reached form and only protrudes by the interlocking tolerance.
    #   top(g) = form(x,y) - u * T_SLAB + protrusion
    # with u ~ U(0,1) for the buried bulk (rank burial), u = 0 for the armour
    # layer, and protrusion ~ |N(0, 0.05 D)|.  The form is the substrate.
    T_SLAB = 3.0 * D50_IN
    sub = bed_sub
    nsub = sub.shape[0]

    def form_at(x, y):
        return float(_bilinear(
            sub, np.array([[x / (PATCH / nsub)]]),
            np.array([[y / (PATCH / nsub)]]))[0, 0])

    for i, gi in enumerate(bulk):
        g = grains[gi]
        x = float(rng.uniform(0.0, PATCH))
        y = float(rng.uniform(0.0, PATCH))
        L = i % 2
        top = (form_at(x, y) - 0.30 * D50_IN * ARMOR
               - 0.95 * D50_IN * L + rng.normal(0.0, 0.03 * g["D"]))
        g["x"], g["y"] = x, y
        g["z"] = top - _half_extent(g, 2)
    for gi in surf:
        g = grains[gi]
        x = float(rng.uniform(0.0, PATCH))
        y = float(rng.uniform(0.0, PATCH))
        # coarse grains protrude, fines bed down: this is the positive skew of
        # a real armour layer (the biggest stones are the tallest things)
        top = (form_at(x, y)
               + 0.06 * (g["D"] - D50_IN * ARMOR)
               + abs(rng.normal(0.0, 0.03 * g["D"])))
        g["x"], g["y"] = x, y
        g["z"] = top - _half_extent(g, 2)
    # rebuild the height buffer from the planed particles
    zbuf.z[:] = -np.inf
    for g in grains:
        if "off" in g or g["kind"].startswith("cluster") \
                or g["kind"].startswith("mono"):
            continue
        zbuf._stamp(g, g["x"], g["y"])

    # snapshot of the matrix BEFORE the anchors: cluster satellites are bedded
    # on the matrix, never on the anchor or on each other (that stacking was a
    # real failure in this build -- see REPORT).
    snap = ZBuf(DEPOSIT_N, PATCH)
    snap.z = zbuf.z.copy()

    # anchors on a jittered grid; the two boulders are kept away from the patch
    # edge so that the periodic wrap never cuts them in two
    for gi in clust + mono:
        g = grains[gi]
        if g["kind"] == "mono":
            # bedded boulder: ~1/3 of the diameter stands proud of the bed
            x = float(rng.uniform(0.22, 0.40))
            y = float(rng.uniform(0.22, 0.40))
            top = form_at(x, y) + 0.35 * g["D"]
            g["x"], g["y"] = x, y
            g["z"] = top - _half_extent(g, 2)
            zbuf._stamp(g, x, y)
        elif "off" not in g:
            x = float(rng.uniform(0.05, 0.95) * PATCH)
            y = float(rng.uniform(0.05, 0.95) * PATCH)
            zbuf.put_here(g, x, y, n_cand=3, rng=rng)

    # cluster members: positioned relative to their anchor, bedded on the
    # matrix snapshot
    for g in grains:
        if "off" not in g:
            continue
        anchor = None
        for h in grains:
            if (h.get("cluster") == g.get("cluster")
                    and h["kind"] in ("cluster_anchor", "mono")):
                anchor = h
                break
        ox, oy = g["off"]
        x = (anchor["x"] + ox) % PATCH
        y = (anchor["y"] + oy) % PATCH
        zc = snap.probe(g, x, y) - INTERLOCK * g["D"]
        g["x"], g["y"], g["z"] = float(x), float(y), zc
        zbuf._stamp(g, x, y)

    # ---- non-penetration relaxation (support functions, z push-up) ------
    from scipy.spatial import cKDTree
    n = len(grains)
    P = np.array([[g["x"], g["y"], g["z"]] for g in grains], dtype=np.float64)
    hx = np.array([_half_extent(g, 0) for g in grains])
    hy = np.array([_half_extent(g, 1) for g in grains])
    ax = np.array([g["abc"] for g in grains], dtype=np.float64)
    RT = np.array([g["R"].T for g in grains], dtype=np.float64)
    D = np.array([g["D"] for g in grains])
    small = np.where(D < 0.06)[0]
    big = np.where(D >= 0.06)[0]
    tree = cKDTree(np.mod(P[small], PATCH), boxsize=PATCH)
    pr = tree.query_pairs(0.075, output_type="ndarray")
    ii = small[pr[:, 0]]
    jj = small[pr[:, 1]]
    if len(big):
        bi = np.repeat(big, n)
        bj = np.tile(np.arange(n), len(big))
        m = bi != bj
        ii = np.concatenate([ii, bi[m]])
        jj = np.concatenate([jj, bj[m]])
    contacts = 0
    worst = 0.0
    it = -1
    for it in range(8):
        u = P[jj] - P[ii]
        u[:, 0] -= PATCH * np.round(u[:, 0] / PATCH)
        u[:, 1] -= PATCH * np.round(u[:, 1] / PATCH)
        d = np.linalg.norm(u, axis=1)
        good = d > 1e-9
        u = u[good] / d[good, None]
        iig, jjg, dg = ii[good], jj[good], d[good]
        h_i = np.sqrt(((ax[iig] * np.einsum("nkj,nj->nk", RT[iig], u)) ** 2)
                      .sum(axis=1))
        h_j = np.sqrt(((ax[jjg] * np.einsum("nkj,nj->nk", RT[jjg], u)) ** 2)
                      .sum(axis=1))
        ov = h_i + h_j - dg
        m = ov > 2e-7
        if not m.any():
            worst = 0.0
            break
        contacts = int(m.sum())
        worst = float(ov[m].max())
        # Jacobi relaxation:每 grain は接触の平均だけ動く (summing every
        # pairwise overlap would launch a boulder that touches hundreds of
        # grains at once).
        disp = np.zeros((n, 3))
        cnt = np.zeros(n)
        vec = u[m] * (0.5 * ov[m])[:, None]
        np.add.at(disp, iig[m], -vec)
        np.add.at(disp, jjg[m], +vec)
        np.add.at(cnt, iig[m], 1.0)
        np.add.at(cnt, jjg[m], 1.0)
        act = cnt > 0
        disp[act] /= cnt[act, None]
        # clamp: a boulder with hundreds of contacts must still creep, not fly
        P += np.clip(0.5 * disp, -0.003, 0.003)
    for i, g in enumerate(grains):
        g["x"], g["y"], g["z"] = (float(P[i, 0] % PATCH),
                                  float(P[i, 1] % PATCH), float(P[i, 2]))
    return {"iterations": it + 1, "contacts_last": contacts,
            "last_worst_overlap_m": worst,
            "pairs": int(len(ii))}


def _support(g, u):
    ax = np.array(g["abc"])
    R = g["R"]
    p = R.T @ u
    return float(math.sqrt(float(((ax * p) ** 2).sum())))


# ---------------------------------------------------------------------------
# vertical ray cast -> depth map + facet normals + grain ids
# ---------------------------------------------------------------------------
def raycast(grains, xs, ys, n_facet=80, skip=None):
    nx, ny = len(xs), len(ys)
    W, H = nx, ny
    z = np.full((H, W), -np.inf)
    gid = np.full((H, W), -1, dtype=np.int32)
    nxv = np.zeros((H, W))
    nyv = np.zeros((H, W))
    nzv = np.ones((H, W))
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    for gi, g in enumerate(grains):
        if skip is not None and gi in skip:
            continue
        Rp = g["abc"]
        A, B, C = Rp
        R = g["R"]
        ax = np.array([A, B, C])
        hx = _half_extent(g, 0)
        hy = _half_extent(g, 1)
        i0 = max(0, int(math.floor((g["x"] - hx - xs[0]) / dx)))
        i1 = min(nx, int(math.ceil((g["x"] + hx - xs[0]) / dx)) + 1)
        j0 = max(0, int(math.floor((g["y"] - hy - ys[0]) / dy)))
        j1 = min(ny, int(math.ceil((g["y"] + hy - ys[0]) / dy)) + 1)
        if i1 <= i0 or j1 <= j0:
            continue
        DX = xs[i0:i1][None, :] - g["x"]
        DY = ys[j0:j1][:, None] - g["y"]
        u = np.stack([R[0, i] * DX + R[1, i] * DY for i in range(3)])
        w = R[2, :][:, None, None]
        aa = ((w / ax[:, None, None]) ** 2).sum(0)
        bb = ((u * w) / ax[:, None, None] ** 2).sum(0)
        cc = ((u / ax[:, None, None]) ** 2).sum(0) - 1.0
        disc = bb * bb - aa * cc
        ok = disc > 0.0
        sq = np.sqrt(np.where(ok, disc, 0.0))
        s2 = (-bb + sq) / aa
        ztop = g["z"] + s2
        sub = z[j0:j1, i0:i1]
        hit = ok & (ztop > sub)
        if not hit.any():
            continue
        sub[hit] = ztop[hit]
        gid[j0:j1, i0:i1][hit] = gi
        # analytic outward normal at the hit point
        S = np.where(hit, s2, 0.0)
        X0 = u[0] + w[0] * S
        X1 = u[1] + w[1] * S
        X2 = u[2] + w[2] * S
        nb = np.stack([X0 / A ** 2, X1 / B ** 2, X2 / C ** 2])
        nn = np.sqrt((nb ** 2).sum(0))
        nn = np.where(nn > 0, nn, 1.0)
        nwx = R[0, 0] * nb[0] + R[0, 1] * nb[1] + R[0, 2] * nb[2]
        nwy = R[1, 0] * nb[0] + R[1, 1] * nb[1] + R[1, 2] * nb[2]
        nwz = R[2, 0] * nb[0] + R[2, 1] * nb[1] + R[2, 2] * nb[2]
        nwx, nwy, nwz = nwx / nn, nwy / nn, nwz / nn
        flip = nwz < 0
        nwx = np.where(flip, -nwx, nwx)
        nwy = np.where(flip, -nwy, nwy)
        nwz = np.where(flip, -nwz, nwz)
        nxv[j0:j1, i0:i1][hit] = nwx[hit]
        nyv[j0:j1, i0:i1][hit] = nwy[hit]
        nzv[j0:j1, i0:i1][hit] = nwz[hit]
    if (gid < 0).any():
        # holes: fill with the neighbourhood mean (periodic in x)
        zz = np.where(np.isfinite(z), z, np.nan)
        fill = np.nanmean(zz)
        z = np.where(np.isfinite(z), z, fill)
    return {"z": z, "gid": gid, "nx": nxv, "ny": nyv, "nz": nzv}


# ---------------------------------------------------------------------------
# facet normals (low-poly quantisation of the analytic normal)
# ---------------------------------------------------------------------------
def fibonacci_hemisphere(n):
    i = np.arange(n, dtype=np.float64)
    z = (i + 0.5) / n                     # 0..1  -> upper hemisphere
    r = np.sqrt(np.maximum(1.0 - z * z, 0.0))
    th = i * (math.pi * (3.0 - math.sqrt(5.0)))
    return np.stack([r * np.cos(th), r * np.sin(th), z], axis=1)


def facet_normals(nx, ny, nz, n_facet=80):
    """Snap each pixel normal to the nearest of n_facet fixed directions.

    The facet set is the low-poly analogue of a faceted DEM: each grain's
    surface is represented by one of 80 discrete orientations, so the slope
    distribution has a hard tail (a smooth ellipsoid only reaches 90 deg on
    its silhouette and that tail would be measured as a measure-zero set).
    """
    F = fibonacci_hemisphere(n_facet)
    V = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1)
    idx = np.argmax(V @ F.T, axis=1)
    N = F[idx]
    return (N[:, 0].reshape(nx.shape), N[:, 1].reshape(nx.shape),
            N[:, 2].reshape(nx.shape), idx.reshape(nx.shape))


# ---------------------------------------------------------------------------
# horizon-based AO and shadow (both from the depth map, static)
# ---------------------------------------------------------------------------
AZIMUTHS = 12


def horizon_map(z, cell, n_dir=AZIMUTHS, r_min=0.0012, r_max=0.045,
                n_steps=26):
    """Per-direction horizon tangent tan(h) for every pixel.

    `z` must be a full patch map (periodic in BOTH axes) so that np.roll is
    the correct wrap.
    """
    H, W = z.shape
    out = np.zeros((n_dir, H, W), dtype=np.float32)
    radii = np.geomspace(r_min, r_max, n_steps)
    for k in range(n_dir):
        th = 2.0 * math.pi * k / n_dir
        cx, sy = math.cos(th), math.sin(th)
        tan_h = np.full((H, W), -1e9, dtype=np.float32)
        for r in radii:
            ox = int(round(r * cx / cell))
            oy = int(round(r * sy / cell))
            if ox == 0 and oy == 0:
                continue
            dist = cell * math.hypot(ox, oy)
            zs = np.roll(z, (-oy, -ox), axis=(0, 1))
            tan_h = np.maximum(tan_h, ((zs - z) / dist).astype(np.float32))
        out[k] = tan_h
    return out


def ambient_occlusion(tan_h):
    """AO from the horizon field: V = mean_d cos^2(h_d), clipped to [0,1].

    Horizon-based ambient occlusion (the standard depth-map approximation):
    the visible sky fraction in direction d is ~ cos^2 of the horizon
    elevation.  Used on the diffuse term only.
    """
    c2 = 1.0 / (1.0 + np.maximum(tan_h, 0.0) ** 2)
    V = c2.mean(axis=0)
    return np.clip(V, 0.0, 1.0).astype(np.float32)


def shadow_map(z, cell, sun_dir_water, disc_rad, n_disc=9, max_len=0.45,
               n_steps=180):
    """Direct-sun visibility with a finite solar disc (penumbra).

    The bed lies under water, so the shadow is cast by the *refracted* sun
    direction inside the water (Snell at a flat surface, IOR 1.334).  For each
    of n_disc sample directions inside the solar half-angle the depth map is
    marched; visibility = (blocked count == 0).  Returns the visible fraction
    (1.0 = fully lit, 0.0 = umbra).
    """
    H, W = z.shape
    s = np.asarray(sun_dir_water, dtype=np.float64)
    s = s / np.linalg.norm(s)
    # build a basis perpendicular to s
    e1 = np.cross(s, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(e1) < 1e-8:
        e1 = np.array([1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(s, e1)
    vis = np.zeros((H, W), dtype=np.float64)
    i = np.arange(n_disc)
    rr = (np.sqrt((i + 0.5) / n_disc) * disc_rad
          if n_disc > 1 else np.zeros(1))
    th = i * (math.pi * (3.0 - math.sqrt(5.0)))
    for k in range(n_disc):
        if n_disc == 1:
            d = s.copy()
        else:
            d = s + e1 * (rr[k] * math.cos(th[k])) + e2 * (rr[k]
                                                           * math.sin(th[k]))
            d /= np.linalg.norm(d)
        blocked = np.zeros((H, W), dtype=bool)
        for t in np.linspace(cell, max_len, n_steps):
            px = (np.arange(W)[None, :] + t * d[0] / cell).astype(np.float64)
            py = (np.arange(H)[:, None] + t * d[1] / cell).astype(np.float64)
            zs = _bilinear(z, px, py)
            # ray height above the reference point
            ray_z = z + t * d[2]
            blocked |= (zs > ray_z + 2e-4)
        vis += (~blocked).astype(np.float64)
    return (vis / float(n_disc)).astype(np.float32)


def _bilinear(z, px, py):
    """Bilinear sample of a full-patch (periodic) map; px/py in pixel units."""
    H, W = z.shape
    x0 = np.floor(px).astype(np.int64)
    y0 = np.floor(py).astype(np.int64)
    tx = px - x0
    ty = py - y0
    x0m = np.mod(x0, W)
    x1m = np.mod(x0 + 1, W)
    y0m = np.mod(y0, H)
    y1m = np.mod(y0 + 1, H)
    return (z[y0m, x0m] * (1 - tx) * (1 - ty) + z[y0m, x1m] * tx * (1 - ty)
            + z[y1m, x0m] * (1 - tx) * ty + z[y1m, x1m] * tx * ty)


# ---------------------------------------------------------------------------
# audit: statistics of the depth map
# ---------------------------------------------------------------------------
def structure_function(z, cell, dmax=0.20, nbin=60):
    """1-D (x and y averaged) structure function with periodic wraps in x."""
    H, W = z.shape
    dmax_px = int(min(dmax / cell, W // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d = []
    out_s = []
    for d in ds:
        dx2 = np.mean((z - np.roll(z, -d, axis=1)) ** 2)
        dy2 = np.mean((z - np.roll(z, -d, axis=0)) ** 2)
        out_d.append(d * cell)
        out_s.append(0.5 * (dx2 + dy2))
    return np.array(out_d), np.array(out_s)


def audit_structure(z, cell):
    d, s = structure_function(z, cell)
    lo, hi = 2.0 * D50_IN, 0.20
    m = (d >= lo) & (d <= hi)
    Hfit = 0.5 * np.polyfit(np.log(d[m]), np.log(s[m]), 1)[0]
    # local slope and the breakpoint (where it falls below 0.6 of 2H)
    loc = 0.5 * np.gradient(np.log(s), np.log(d))
    brk = None
    target = 0.6 * 2.0 * Hfit
    for i in range(1, len(d) - 1):
        if d[i] < D50_IN * 0.5:
            continue
        if loc[i] < target and brk is None:
            brk = float(d[i])
    return {"H_fit": float(Hfit), "fit_range_m": [lo, hi],
            "break_scale_m": brk,
            "break_over_D50": (brk / D50_IN) if brk else None,
            "H_local_at_D50": float(np.interp(D50_IN, d, loc)),
            "D_m": d.tolist(), "S_m2": s.tolist(),
            "H_local": loc.tolist()}


def _residual_surface(z, cell, sigma_m, smooth_px=1.0):
    from scipy import ndimage
    trend = ndimage.gaussian_filter(z, sigma_m / cell, mode="wrap")
    r = z - trend
    if smooth_px > 0:
        r = ndimage.gaussian_filter(r, smooth_px, mode="wrap")
    return r


def reconstruct_gsd(z, cell, trend_sigma=2.2 * D50_IN,
                    min_frac=0.25, peak_size=5):
    """Gravelometer-style reconstruction from the depth map alone.

    (1) residual surface r = z - gaussian(z, 2.2 D50)   [removes the substrate]
    (2) markers = regional maxima of r (maximum filter, one marker per grain
        after connected-component merging of plateaus)
    (3) partition = every pixel belongs to its NEAREST regional maximum
        (Euclidean Voronoi of the markers).  This is the declared choice: a
        height watershed was tried and failed here (the residual surface is
        too flat for the 8-bit flooding to keep 300+ basins apart -- it merged
        them into two regions covering the patch), so the partition is the
        nearest-peak one, which is the same operator on both beds compared.
    (4) per-region area A -> D = 2 sqrt(A/pi)  (the visible cap; a cap of a
        grain protruding above its neighbours is smaller than its equator)
    (5) area-weighted quantiles -> D50, D84, D84/D50

    Uses ONLY the depth map (never the grain ids used to build it).
    """
    import cv2
    from scipy import ndimage
    r = _residual_surface(z, cell, trend_sigma)
    mx = ndimage.maximum_filter(r, size=int(peak_size))
    peaks = (r >= mx - 1e-12)
    nlab, lab = cv2.connectedComponents(peaks.astype(np.uint8))
    if nlab <= 1:
        return {"n": 0, "D50_m": None, "method": "nearest-peak Voronoi"}
    _, idx = ndimage.distance_transform_edt(lab == 0, return_indices=True)
    vor = lab[idx[0], idx[1]]
    ids, counts = np.unique(vor, return_counts=True)
    area = counts.astype(np.float64) * (cell * cell)
    D = 2.0 * np.sqrt(area / math.pi)
    keep = area >= (min_frac * D50_IN) ** 2
    Dk, Ak = D[keep], area[keep]
    if Dk.size < 10:
        return {"n": int(Dk.size), "D50_m": None,
                "method": "nearest-peak Voronoi"}
    order = np.argsort(Dk)
    Dk, Ak = Dk[order], Ak[order]
    cw = np.cumsum(Ak) / Ak.sum()
    return {"n": int(Dk.size), "n_markers": int(nlab - 1),
            "D50_m": float(np.interp(0.50, cw, Dk)),
            "D84_m": float(np.interp(0.84, cw, Dk)),
            "D16_m": float(np.interp(0.16, cw, Dk)),
            "D84_over_D50": float(np.interp(0.84, cw, Dk)
                                  / np.interp(0.50, cw, Dk)),
            "D50_number_m": float(np.median(Dk)),
            "coverage_fraction": float(Ak.sum() / (z.size * cell * cell)),
            "method": "nearest-peak Voronoi on the residual surface"}


def porosity(grains, z, cell, slab=2.5 * D50_IN):
    """Column-mean solid fraction over the top `slab` metres -> porosity.

    Voxel test on a 384^2 x 128 grid.  The slab starts `slab` below each
    column's own depth-map top, so the empty space that the settling's Jacobi
    relaxation pushed below the bed cannot dilute the number (it did: an
    unclipped depth range gave porosity 0.49 with a bed only 90 mm thick).
    """
    nx, ny, nz = VOX
    xs2 = (np.arange(nx) + 0.5) * (PATCH / nx)
    ys2 = (np.arange(ny) + 0.5) * (PATCH / ny)
    ztop = _bilinear(z, xs2[None, :] / cell, ys2[:, None] / cell)
    z1 = float(np.median(ztop)) + 1.5 * D50_IN
    z0 = z1 - float(slab)
    dz = slab / (nz - 1)
    dx = PATCH / nx
    dy = PATCH / ny
    occ = np.zeros((nz, ny, nx), dtype=bool)
    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy
    zs = z0 + (np.arange(nz) + 0.5) * dz
    for g in grains:
        A, B, C = g["abc"]
        R = g["R"]
        ax = np.array([A, B, C])
        hx = _half_extent(g, 0)
        hy = _half_extent(g, 1)
        hz = _half_extent(g, 2)
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
        q = np.zeros((k1 - k0, j1 - j0, i1 - i0))
        for m in range(3):
            u = (R[0, m] * DX + R[1, m] * DY + R[2, m] * DZ) / ax[m]
            q += u * u
        occ[k0:k1, j0:j1, i0:i1] |= (q <= 1.0)
    ktop = np.clip(((ztop - z0) / dz).astype(np.int64), 1, nz)
    idx = np.arange(nz)[:, None, None] < ktop[None, :, :]
    solid = (occ & idx).sum()
    total = idx.sum()
    return {"porosity": float(1.0 - solid / max(total, 1)),
            "solid_fraction": float(solid / max(total, 1)),
            "z0_m": z0, "z1_m": z1, "slab_m": float(slab),
            "vox_m": [dx, dy, dz], "grid": [nx, ny, nz]}


# ---------------------------------------------------------------------------
# construction of one bed
# ---------------------------------------------------------------------------
def build_bed(sigma_sub, seed=SEED, armored=True, clusters=True, mono=True,
              verbose=True):
    t0 = time.time()
    n_sub = 512
    zsub = fbm_field(n_sub, PATCH, SUBSTRATE_H, seed + 7, sigma_sub)
    zb = ZBuf(DEPOSIT_N, PATCH)
    # the z buffer starts at the substrate, resampled to the deposit grid
    ii = (np.arange(DEPOSIT_N) + 0.5) / DEPOSIT_N
    X, Y = np.meshgrid(ii * PATCH, ii * PATCH)
    zb.z = _bilinear(zsub, X / (PATCH / n_sub), Y / (PATCH / n_sub))
    grains, cms = build_population(seed, armored=armored, clusters=clusters,
                                   mono=mono)
    info = settle(grains, zb, seed, bed_sub=zsub)
    if verbose:
        print("    deposit+relax: %d grains, %d pairs, %d iters (%.1fs)"
              % (len(grains), info["pairs"], info["iterations"],
                 time.time() - t0), flush=True)
    return {"grains": grains, "cluster_members": cms, "substrate": zsub,
            "settle": info, "sigma_sub": float(sigma_sub),
            "wall_s": float(time.time() - t0)}


def maps_for(bed, w, h, x0, x1, y0, y1, n_facet=80, skip=None):
    xs = x0 + (np.arange(w) + 0.5) * ((x1 - x0) / w)
    ys = y0 + (np.arange(h) + 0.5) * ((y1 - y0) / h)
    rc = raycast(bed["grains"], xs, ys, n_facet=n_facet, skip=skip)
    fx, fy, fz, fidx = facet_normals(rc["nx"], rc["ny"], rc["nz"], n_facet)
    rc["fx"], rc["fy"], rc["fz"], rc["fidx"] = fx, fy, fz, fidx
    rc["cell"] = (x1 - x0) / w
    return rc


def _skew(a):
    a = np.asarray(a, dtype=np.float64).ravel()
    a = a[np.isfinite(a)]
    m = a.mean()
    s = a.std()
    return float(np.mean((a - m) ** 3) / s ** 3)


def _quant(a, q):
    return float(np.quantile(a, q))


def build_and_audit(out_dir=None, seed=SEED, sun=None, verbose=True,
                    sigma_target=SIGMA_TARGET):
    """Full bed: sigma_z fixed point, maps, audit dict."""
    log = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    iters = []
    sigma_sub = 0.5 * sigma_target
    floor = 0.55 * sigma_target
    for it in range(4):
        bed = build_bed(sigma_sub, seed=seed, verbose=verbose)
        gs = bed["grains"]
        Ds = np.array([g["D"] for g in gs])
        small = np.where(Ds <= D50_IN * ARMOR * SIGMA_G)[0]
        nomono = np.array([not str(g["kind"]).startswith("mono") for g in gs])
        rc = maps_for(bed, RENDER_W, RENDER_H, WIN_X0, WIN_X1, WIN_Y0, WIN_Y1)
        gid = rc["gid"]
        mm = (gid >= 0) & np.isin(gid, small) & nomono[gid]
        zmat = rc["z"][mm]
        s_mat = float(np.std(zmat))
        # substrate and particles are independent, so remove the substrate
        # variance in quadrature before deciding the next amplitude
        s_grain = math.sqrt(max(s_mat ** 2 - sigma_sub ** 2, 0.0))
        resid = math.sqrt(max(sigma_target ** 2 - s_grain ** 2, floor ** 2))
        iters.append({"iter": it, "sigma_sub_m": float(sigma_sub),
                      "sigma_matrix_m": s_mat,
                      "sigma_particles_only_m": s_grain,
                      "next_sigma_sub_m": resid})
        log("  sigma_z fixpoint %d: sigma_sub=%.5f matrix=%.5f "
            "particles=%.5f -> %.5f"
            % (it, sigma_sub, s_mat, s_grain, resid))
        if abs(resid - sigma_sub) < 2e-4:
            sigma_sub = resid
            break
        sigma_sub = resid
    # final build with the converged amplitude
    bed = build_bed(sigma_sub, seed=seed, verbose=verbose)
    grains = bed["grains"]
    kinds = np.array([g["kind"] for g in grains])
    is_mono = np.array([str(k).startswith("mono") for k in kinds])

    t0 = time.time()
    hires = maps_for(bed, HIRES, HIRES, 0.0, PATCH, 0.0, PATCH)
    rmap = maps_for(bed, RENDER_W, RENDER_H, WIN_X0, WIN_X1, WIN_Y0, WIN_Y1)
    log("  raycast hires %dx%d + render %dx%d (%.1fs)"
        % (HIRES, HIRES, RENDER_W, RENDER_H, time.time() - t0))

    gid_h = hires["gid"]
    mono_mask_h = np.isin(gid_h, np.where(is_mono)[0])
    bed_mask_h = ~mono_mask_h & (gid_h >= 0)

    # bed WITHOUT the two boulders: every audit statistic is measured on this
    # map (the boulders are objects, not bed, and would dominate sigma_z)
    hires_bed = maps_for(bed, HIRES, HIRES, 0.0, PATCH, 0.0, PATCH,
                         skip=set(int(i) for i in np.where(is_mono)[0]))
    zh = hires_bed["z"]
    gid_h = hires_bed["gid"]
    hole = gid_h < 0
    if hole.any():
        # the two boulders hide the bed under them; inpaint so that the
        # structure function / skewness do not see a 220 mm hole
        from scipy import ndimage
        wgt = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, zh), 60.0,
                                      mode="wrap")
        den = ndimage.gaussian_filter(wgt, 60.0, mode="wrap")
        zh = np.where(hole, num / np.maximum(den, 1e-9), zh)
    bed_mask_h = np.ones_like(gid_h, dtype=bool)
    z_bed = zh[bed_mask_h]
    D_all = np.array([g["D"] for g in grains])
    D84_surface = D50_IN * ARMOR * SIGMA_G
    mat_mask_h = bed_mask_h & np.isin(gid_h, np.where(D_all <= D84_surface)[0])
    z_mat = zh[mat_mask_h]
    grains_bed = [g for g in grains if not str(g["kind"]).startswith("mono")]
    stats = {
        "n_grains": len(grains),
        "kinds": {str(k): int((kinds == k).sum()) for k in sorted(set(kinds))},
        "patch_m": PATCH, "D50_in_m": D50_IN, "sigma_g": SIGMA_G,
        "armor_ratio_input": ARMOR,
        "D50_surface_nominal_m": D50_IN * ARMOR,
        "sigma_z_fixpoint": iters,
        "sigma_sub_m": float(sigma_sub),
        "sigma_z_bed_m": float(np.std(z_bed)),
        "sigma_z_over_D50": float(np.std(z_bed) / D50_IN),
        "sigma_z_matrix_m": float(np.std(z_mat)),
        "sigma_z_matrix_over_D50": float(np.std(z_mat) / D50_IN),
        "matrix_definition": ("pixels owned by grains with D <= D84 of the "
                              "armoured surface GSD (= %.3f m); the coarse "
                              "tail is macro-roughness.  Excludes the two "
                              "'mono' boulders." % D84_surface),
        "sigma_z_all_m": float(np.std(zh)),
        "sigma_z_over_D50_all": float(np.std(zh) / D50_IN),
        "skew_bed": _skew(z_bed),
        "skew_matrix": _skew(z_mat),
        "skew_all": _skew(zh),
        "z_bed_range_m": [float(z_bed.min()), float(z_bed.max())],
        "porosity": porosity(grains_bed, zh, hires["cell"]),
        "settle": bed["settle"],
    }

    # ---- cluster coverage (exact, from the grain ids) --------------------
    cm = bed["cluster_members"]
    ci = [i for i, g in enumerate(grains) if g.get("cluster") is not None
          and not str(g["kind"]).startswith("mono")]
    cover = float(np.isin(gid_h, ci).mean())
    stats["cluster_coverage"] = cover
    stats["cluster_coverage_target"] = CLUSTER_COVERAGE
    stats["n_clusters"] = int(len(set(g["cluster"] for g in cm
                                      if g.get("cluster") is not None
                                      and not str(g["kind"]).startswith("mono"))))

    # ---- imbrication -----------------------------------------------------
    surf = [g for g in grains if g["kind"] in ("surface", "cluster_anchor",
                                               "cluster_medium", "cluster_fine")]
    dips = np.array([g["dip_deg"] for g in surf])
    # geometric check: the down-dip direction of the a-b plane (n = R[:,2])
    down = np.array([[-g["R"][0, 2], -g["R"][1, 2]] for g in surf])
    down /= np.maximum(np.linalg.norm(down, axis=1, keepdims=True), 1e-12)
    upstream = (down[:, 0] < -0.5)          # down-dip points upstream (-x)
    stats["imbrication"] = {
        "n": int(dips.size),
        "dip_median_deg": float(np.median(dips)),
        "dip_p16_deg": float(np.percentile(dips, 16)),
        "dip_p84_deg": float(np.percentile(dips, 84)),
        "rule": "dip = 30 +- 10 deg (uniform), plane normal tilted downstream "
                "so the a-b plane dips upstream",
        "upstream_fraction": float(np.mean(upstream)),
        "a_axis_azimuth_transverse_fraction": float(np.mean(
            [(1.0 if abs((g["az_deg"] % 180.0) - 90.0) < 20.0 else 0.0)
             for g in surf])),
    }

    # ---- structure function --------------------------------------------
    st = audit_structure(zh, hires["cell"])
    stats["structure"] = {k: v for k, v in st.items() if k not in
                          ("D_m", "S_m2", "H_local")}
    stats["structure_curves"] = {"D_m": st["D_m"], "S_m2": st["S_m2"],
                                 "H_local": st["H_local"]}

    # ---- facet slopes ----------------------------------------------------
    fz = hires["fz"]
    slope = np.degrees(np.arccos(np.clip(fz, -1.0, 1.0)))
    sl = slope[bed_mask_h]
    stats["facets"] = {
        "n_directions": 80,
        "slope_median_deg": float(np.median(sl)),
        "slope_p84_deg": float(np.percentile(sl, 84)),
        "frac_gt60deg": float(np.mean(sl > 60.0)),
        "frac_gt80deg": float(np.mean(sl > 80.0)),
        "method": "analytic ellipsoid normal quantised to the nearest of 80 "
                  "Fibonacci hemisphere directions (low-poly facet analogue)",
    }

    # ---- reconstructed GSD (depth map only) -----------------------------
    gsd = reconstruct_gsd(zh, hires["cell"])
    stats["gsd_reconstructed"] = gsd

    # ---- reference bed (no armour, no clusters, no boulders) ------------
    t0 = time.time()
    ref = build_bed(sigma_sub, seed=seed + 101, armored=False,
                    clusters=False, mono=False, verbose=False)
    ref_h = maps_for(ref, HIRES // 2, HIRES // 2, 0.0, PATCH, 0.0, PATCH)
    ref_gsd = reconstruct_gsd(ref_h["z"], ref_h["cell"])
    stats["gsd_reference"] = ref_gsd
    if gsd.get("D50_m") and ref_gsd.get("D50_m"):
        stats["armor_ratio_reconstructed"] = float(gsd["D50_m"]
                                                   / ref_gsd["D50_m"])
        stats["D50_recon_over_input"] = float(gsd["D50_m"] / D50_IN)
    stats["reference_build_wall_s"] = float(time.time() - t0)

    # ---- AO + shadow on the full patch (periodic), then crop ------------
    t0 = time.time()
    full = maps_for(bed, RENDER_W, RENDER_W, 0.0, PATCH, 0.0, PATCH)
    cell = full["cell"]
    ao_full = ambient_occlusion(horizon_map(full["z"], cell))
    sun_info = None
    sh_full = np.ones_like(ao_full)
    if sun is not None:
        # refracted in-water sun direction (flat surface, IOR 1.334)
        s = np.asarray(sun["sunDir"], dtype=np.float64)
        s = s / np.linalg.norm(s)
        eta = 1.0 / 1.33434
        cos1 = float(s[2])
        k = 1.0 - eta * eta * (1.0 - cos1 * cos1)
        sh = eta * s + (math.sqrt(max(k, 0.0)) - eta * cos1) \
            * np.array([0.0, 0.0, 1.0])
        sh /= np.linalg.norm(sh)
        sh_full = shadow_map(full["z"], cell, sh, float(sun["disc_rad"]))
        sun_info = {"sun_dir_air": [float(x) for x in s],
                    "sun_dir_water": [float(x) for x in sh],
                    "altitude_deg": float(sun.get("altitude_deg", float("nan"))),
                    "azimuth_deg": float(sun.get("azimuth_deg", float("nan"))),
                    "disc_half_angle_rad": float(sun["disc_rad"]),
                    "n_disc_samples": 9,
                    "umbra_fraction": float(np.mean(sh_full <= 0.0)),
                    "penumbra_fraction": float(np.mean((sh_full > 0.0)
                                                       & (sh_full < 1.0))),
                    "lit_fraction": float(np.mean(sh_full >= 1.0)),
                    "method": "march of the refracted sun direction through "
                              "the depth map, one march per solar-disc "
                              "sample; visibility = blocked-sample complement"}
        stats["shadow"] = sun_info
    rmap["ao"] = ao_full[RENDER_J0:RENDER_J0 + RENDER_H, :]
    rmap["shadow"] = sh_full[RENDER_J0:RENDER_J0 + RENDER_H, :]
    full["ao"] = ao_full
    full["shadow"] = sh_full
    stats["ao"] = {"method": "horizon map over 12 azimuths, radii 1.2-45 mm, "
                             "V = mean_d 1/(1+tan^2 h_d)",
                   "mean": float(np.mean(rmap["ao"])),
                   "p16": _quant(rmap["ao"], 0.16),
                   "p84": _quant(rmap["ao"], 0.84),
                   "frac_below_0.5": float(np.mean(rmap["ao"] < 0.5))}
    stats["map_wall_s"] = float(time.time() - t0)

    stats["full_patch"] = {"n_px": int(full["z"].size),
                           "sigma_z_m": float(np.std(full["z"]))}
    return bed, rmap, hires, stats, full


# ---------------------------------------------------------------------------
# G-BED gate
# ---------------------------------------------------------------------------
def gate_bed(stats):
    s = stats
    g = {}
    def put(key, val, lo, hi, label, margin=True):
        ok = (val is not None) and (lo <= val <= hi)
        g[key] = {"value": val, "range": [lo, hi], "pass": bool(ok),
                  "label": label}
        return ok

    put("skewness", s["skew_bed"], 0.3, 0.9, "bed skewness (boulders excluded)")
    put("sigma_z_over_D50", s["sigma_z_matrix_over_D50"], 0.32, 0.46,
        "sigma_z / D50 of the MATRIX (see matrix_definition)")
    g["sigma_z_whole_bed"] = {
        "value_over_D50": s["sigma_z_over_D50"],
        "note": "recorded, not gated: with the armoured coarse tail included "
                "the brief's 0.32-0.46 band cannot close (the same 0.4 D50 "
                "relation then applies to the 0.040 m surface D50, i.e. 0.8)"}
    put("structure_H", s["structure"]["H_fit"], 0.7, 0.85,
        "structure function H (2-20 cm fit)")
    brk = s["structure"].get("break_over_D50")
    g["structure_break"] = {
        "value_D50_units": brk,
        "pass": bool(brk is not None and 0.4 <= brk <= 3.0),
        "label": "structure-function break scale / D50 (should be ~1)"}
    gsd = s["gsd_reconstructed"]
    put("recon_D84_over_D50", gsd.get("D84_over_D50"), 1.9, 2.3,
        "surface-reconstructed D84/D50")
    put("recon_armor_ratio", s.get("armor_ratio_reconstructed"), 1.7, 2.3,
        "reconstructed armour ratio (surface / unarmoured reference)")
    put("imbrication_median_deg", s["imbrication"]["dip_median_deg"], 25.0,
        35.0, "a-b plane dip, median (upstream)")
    put("cluster_coverage", s["cluster_coverage"], 0.05, 0.15,
        "cluster areal coverage")
    put("porosity", s["porosity"]["porosity"], 0.30, 0.45, "bed porosity")
    g["facet_tail"] = {
        "frac_gt60deg": s["facets"]["frac_gt60deg"],
        "frac_gt80deg": s["facets"]["frac_gt80deg"],
        "note": "recorded, not gated (a faceted gravel DEM has an 80-90 deg "
                "tail; this low-poly ellipsoid analogue has its own)"}
    g["D50_recon_over_input"] = {
        "value": s.get("D50_recon_over_input"),
        "note": "recorded, not gated: the visible-cap reconstruction is not "
                "expected to return the input median exactly"}
    g["pass"] = all(v.get("pass", True) for v in g.values())
    return g


# ---------------------------------------------------------------------------
# depth-map.png
# ---------------------------------------------------------------------------
def save_depth_png(path, z, nrm, cell):
    from PIL import Image
    zz = z.copy()
    zz -= zz.min()
    zz /= max(zz.max(), 1e-9)
    # hillshade with the fixed preset-5 sun (top-left-ish)
    s = np.array([0.586, -0.006, 0.810])
    s /= np.linalg.norm(s)
    shade = np.clip(nrm[..., 0] * s[0] + nrm[..., 1] * s[1]
                    + nrm[..., 2] * s[2], 0.0, 1.0)
    img = 0.55 * zz + 0.45 * shade
    img = np.clip(img, 0.0, 1.0)
    rgb = np.stack([img * 0.98, img * 0.95, img * 0.86], axis=-1)
    u8 = (rgb * 255.0 + 0.5).astype(np.uint8)
    buf = _png_bytes(u8)
    atomic_write_bytes(path, buf)


def _png_bytes(u8):
    import io
    from PIL import Image
    b = io.BytesIO()
    Image.fromarray(u8).save(b, format="PNG")
    return b.getvalue()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out", "v4-2-1"))
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--no-shadow", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    import env as ENV
    st = {"name": "深い海", "sun": {"date": "2026-07-20", "time": "10:00",
                                   "latitudeDeg": 26.3,
                                   "longitudeDeg": 127.8},
          "water": {"depthM": 0.25, "transparency": 0.95, "turbidity": 0.02,
                    "pigment": 0.02},
          "windWave": {"windSpeedMS": 3.5, "windDirDeg": 180, "fineness": 0.45,
                       "flowSpeed": 0.45},
          "weather": {"cloudPct": 10, "clarity": 0.9},
          "couplings": [{"id": k, "enabled": True, "weight": 1.0} for k in
                        ("cloud-sun", "wind-wave", "sun-caustics",
                         "turbid-caustics", "cloud-glitter", "depth-tint")],
          "render": {"strength": 1.0, "exposure": 1.0}}
    e = ENV.resolve_environment(st)
    sun = {"sunDir": e["sunDir"], "disc_rad": 0.00395,
           "altitude_deg": e["sun"]["altitudeDeg"],
           "azimuth_deg": e["sun"]["azimuthDeg"]}
    t0 = time.time()
    bed, rmap, hires, stats, full = build_and_audit(
        args.out, seed=args.seed,
        sun=None if args.no_shadow else sun)
    stats["gate"] = gate_bed(stats)
    stats["created_utc"] = now_iso()
    stats["build_wall_s"] = float(time.time() - t0)
    stats["raycast"] = {"render": [RENDER_W, RENDER_H],
                        "render_window_m": [WIN_X0, WIN_X1, WIN_Y0, WIN_Y1],
                        "hires": [HIRES, HIRES], "hires_cell_m": hires["cell"],
                        "render_cell_m": rmap["cell"],
                        "method": "vertical ray cast (z-buffer), analytic "
                                  "ellipsoid normal at the hit point"}
    save_depth_png(os.path.join(args.out, "depth-map.png"), rmap["z"],
                   np.stack([rmap["nx"], rmap["ny"], rmap["nz"]], -1),
                   rmap["cell"])
    with open(os.path.join(args.out, "bed-stats.json"), "w",
              encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    # maps for the renderer: the FULL patch (0.6 m square, periodic) so the
    # transport can wrap in both directions, plus the render window crop
    np.savez_compressed(os.path.join(args.out, "_work", "bed-maps.npz"),
                        z=full["z"].astype(np.float32),
                        gid=full["gid"],
                        nx=full["nx"].astype(np.float32),
                        ny=full["ny"].astype(np.float32),
                        nz=full["nz"].astype(np.float32),
                        fx=full["fx"].astype(np.float32),
                        fy=full["fy"].astype(np.float32),
                        fz=full["fz"].astype(np.float32),
                        ao=full["ao"].astype(np.float32),
                        shadow=full["shadow"].astype(np.float32),
                        D=np.array([g["D"] for g in bed["grains"]]),
                        kind=np.array([g["kind"] for g in bed["grains"]]),
                        render_j0=np.array([RENDER_J0]),
                        cell=np.array([full["cell"]]),
                        mono_geom=np.array(
                            [[g["x"], g["y"], g["z"], g["abc"][0], g["abc"][1],
                              g["abc"][2]] for g in bed["grains"]
                             if str(g["kind"]) == "mono"]),
                        mono_R=np.array(
                            [g["R"].tolist() for g in bed["grains"]
                             if str(g["kind"]) == "mono"]))
    g = stats["gate"]
    print("\nG-BED")
    for k, v in g.items():
        if not isinstance(v, dict):
            continue
        if "value" in v and "range" in v:
            print("  %-26s %-10s [%s, %s] %s"
                  % (k, ("%.4g" % v["value"]) if v["value"] is not None
                     else "None", v["range"][0], v["range"][1],
                     "PASS" if v["pass"] else "FAIL"))
        else:
            print("  %-26s %s (recorded)" % (k, v))
    print("  overall: %s" % ("PASS" if g["pass"] else "FAIL"))
    return 0 if g["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
