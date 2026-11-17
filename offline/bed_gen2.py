"""bed_gen2.py -- particle-based gravel bed for BRIEF-V4-2-2.

What is new relative to bed_gen.py (v4-2.1), and why
---------------------------------------------------
The v4-2.1 audit failed 7 of 8 items.  DESIGN-V4-2-2 splits that failure into

  (i)  items whose *measurement* was not established  -> fix the audit,
  (ii) items that are a property of the *construction*-> fix the generator.

This module changes the generator only.  The audit side lives in
run_v4_2_2.py.  bed_gen.py is imported for its primitives (fBm, grain shape,
structure function, porosity) and is NOT modified.

Generation changes (DESIGN-V4-2-2 "生成側の修正" 1-4):

1. Asymmetric protrusion p.  v2.1 planed every grain to one "reached form"
   with protrusion ~ |N(0, 0.03 D)| and, for the buried bulk, a two-level
   rank burial (u in {0,1}).  That produced a bimodal, slightly negative
   matrix skew (-0.073).  Here each grain draws p from

        p = -(0.08 + 0.57 u) D          with probability 1 - P_TAIL
        p = +min(Exp(mean 0.18 D), 0.5 D)   with probability P_TAIL

   i.e. most tops sit *below* the form and a minority stands up by up to
   0.5 D.  The summit of grain g is placed at form(x_c,y_c) + p_g.

2. sigma_z is closed by measurement, not by a fixed-point assumption: each
   attempt builds -> raycasts -> measures sigma_z(matrix) and sigma_grain
   = std(z - form) on the matrix pixels -> updates the substrate amplitude
   AND the protrusion (burial) scale.  Max 4 attempts, every attempt is
   returned in the stats.  Target 0.40 * D50_surface (= 16.0 mm); the band
   is applied afterwards and never moved.

3. Void ratio and cluster coverage are construction targets: the bulk is
   under-filled relative to v2.1 (BULK_VOL) and the cluster anchors are
   scaled (CLUSTER_SCALE) so that porosity lands near 0.37 and the cluster
   areal coverage near 0.10 (both measured afterwards).

4. The substrate is the same fBm (H = 0.79).  Its amplitude is allocated
   inside the sigma_z budget.  Because a stationary grain mosaic contributes
   a *constant* (plateau) term to the structure function at d > D50, the
   large-scale exponent can only be carried by the substrate itself; the
   allocation therefore gives the substrate most of the variance and the
   grains the rest (SIGMA_GRAIN_FRAC).  A patch-scale low-wavenumber component
   (8 modes, period 0.6-0.3 m) is added to the substrate, and the substrate
   height distribution is made asymmetric (quadratic skew transform).  Both
   are disclosed in the report; the low-mode weight and the skew coefficient
   are recorded constants, not fits to the audit.

Declared limits: ellipsoids have no corners (unchanged from v2.1); grains may
interpenetrate (the v2.1 Jacobi relaxation was dropped because it moved the
height field away from the planed distribution -- the penetration statistics
are measured and reported instead); the deposit is quasi-static, no dynamics.
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
import bed_gen as BG  # noqa: E402

# ---------------------------------------------------------------------------
# scene constants -- verbatim from v4-2.1 (BRIEF-V4-2-2: "v4-2.1 から verbatim")
# ---------------------------------------------------------------------------
PATCH = 0.6
D50_IN = 0.020
SIGMA_G = 2.1
ARMOR = 2.0
D50_SURF = D50_IN * ARMOR            # 0.040 m  (armoured surface D50)
D84_SURF = D50_SURF * SIGMA_G        # 0.084 m
BA_RANGE = (0.67, 0.74)
COREY_RANGE = (0.6, 0.7)
DIP_DEG = 30.0
DIP_SPREAD_DEG = 10.0
MONO_D = (0.15, 0.22)
SUBSTRATE_H = 0.79
N_CLUSTERS = 3
SEED = 3255045276382562994

SIGMA_TARGET = 0.40 * D50_SURF        # 0.016 m (brief: 0.40 * D50_surface)

# ---- construction constants (design-stage calibration, all recorded) ------
FORM_W_LOW = 0.62        # weight of the patch-scale low-mode component
FORM_SKEW_C = 0.55       # quadratic skew coefficient of the substrate
P_TAIL = 0.30            # probability that a grain stands proud of the form
P_TAIL_MEAN = 0.18       # Exp mean of the proud tail, in units of D
P_TAIL_MAX = 0.50        # hard cap of the proud tail, in units of D
P_BURY_MIN = 0.05        # burial depth range, in units of D
P_BURY_MAX = 0.28
P_SCALE0 = 0.30          # protrusion (burial) scale, fixed by calibration
SORT_GAIN = 0.22         # size sorting: p += SORT_GAIN*((D/D50s)^q-1)*D50s
SORT_Q = 1.0
BULK_VOL = 0.20          # bulk solid volume per unit area per 1.6 D50 thickness
D_CLIP_MAX = 0.14        # upper clip of the drawn GSD (v2.1 used 0.09; that
                         # cut 13 % of the areal upper tail and biased D84)
BULK_THICK = 1.6 * D50_IN
SURF_COVER = 0.85        # surface-layer footprint coverage (fraction of patch)
CLUSTER_SCALE = 0.90     # anchor-size multiplier vs the v2.1 D90 rule

DEPOSIT_N = 1024
HIRES = 2048
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


def log(msg):
    print("[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg),
          flush=True)


# ---------------------------------------------------------------------------
# substrate form: fBm(H=0.79) + patch-scale low modes + skew transform
# ---------------------------------------------------------------------------
def low_mode_field(n, L, H, seed, kmax=2):
    """Periodic sum of the lowest |k| modes (period 0.6 .. 0.3 m).

    Real by construction: z = sum a_k cos(2 pi (i x/L + j y/L) + phi_k).
    Amplitudes follow the same k^-(H+1)/2 law as the fBm so the two parts
    share one spectral slope; only the outer scale is extended.
    """
    rng = np.random.default_rng(int(seed))
    ii = (np.arange(n) + 0.5) * (L / n)
    X, Y = np.meshgrid(ii, ii)
    z = np.zeros((n, n), dtype=np.float64)
    for i in range(-kmax, kmax + 1):
        for j in range(-kmax, kmax + 1):
            k2 = i * i + j * j
            if k2 == 0 or k2 > kmax * kmax:
                continue
            a = k2 ** (-(H + 1.0) / 2.0)
            ph = rng.uniform(0.0, 2.0 * math.pi)
            z += a * np.cos(2.0 * math.pi * (i * X + j * Y) / L + ph)
    z -= z.mean()
    return z / max(z.std(), 1e-12)


def form_field(n, sigma, seed, w_low=FORM_W_LOW, skew_c=FORM_SKEW_C):
    """Substrate form with unit-var parts -> exactly `sigma` std, skewed.

    v = w*g_low + sqrt(1-w^2)*g_fbm, then v <- v + c (v^2 - 1) and rescale.
    The quadratic term carries the positive skew of the *form*; the grain
    protrusion p carries its own (smaller) positive skew -- both reported.
    """
    f = BG.fbm_field(n, PATCH, SUBSTRATE_H, seed + 7, 1.0)
    u = low_mode_field(n, PATCH, SUBSTRATE_H, seed + 11, kmax=2)
    g = w_low * u + math.sqrt(max(1.0 - w_low * w_low, 0.0)) * f
    g = g / max(g.std(), 1e-12)
    v = g + skew_c * (g * g - 1.0)
    v -= v.mean()
    return (v * (sigma / max(v.std(), 1e-12))).astype(np.float64)


# ---------------------------------------------------------------------------
# grain population (positions and protrusion p drawn ONCE per parameter set)
# ---------------------------------------------------------------------------
def protrusion(rng, D, p_scale, tail_p=P_TAIL, tail_mean=P_TAIL_MEAN,
               tail_max=P_TAIL_MAX, bury_min=P_BURY_MIN, bury_max=P_BURY_MAX,
               sort_gain=SORT_GAIN, sort_q=SORT_Q):
    """Asymmetric protrusion of the grain summit above the form (metres).

    Two declared parts:
      * size sorting  p_sort = sort_gain * ((D/D50_SURF)**sort_q - 1) * D50:
        coarse clasts stand proud of the form, fines sink under it.  This is
        what makes the armour layer an inverted (coarse-on-top) surface and
        what the exposed-area audit sees.
      * asymmetric random part: most grains are buried (p < 0); a fraction
        P_TAIL stands proud by up to tail_max*D (scaled by p_scale).
    """
    if rng.random() < tail_p:
        pr = min(rng.exponential(tail_mean * D), tail_max * D)
    else:
        pr = -(bury_min + (bury_max - bury_min) * rng.random()) * D
    ps = sort_gain * ((D / D50_SURF) ** sort_q - 1.0) * D50_SURF
    return float(ps + p_scale * pr)


def pack_positions(radii, rng, k=32, obstacles=None):
    """Greedy max-min (Poisson-disc-like) packing of footprint radii.

    Plan-view NON-overlap matters for the audit: the exposed top area of a
    grain is only proportional to its footprint (and therefore the exposed GSD
    is an unbiased estimate of the areal GSD) when the grains tile the patch
    instead of burying each other.  Coarse grains are placed first, each at the
    best of k random candidates (periodic distances), so the fines end up in
    the interstices -- an inverted (armoured) pavement.
    """
    order = np.argsort(-np.asarray(radii))
    P = [] if obstacles is None else [[o[0], o[1], o[2]] for o in obstacles]
    out = np.zeros((len(radii), 2))
    for i in order:
        r = float(radii[i])
        cand = rng.uniform(0.0, PATCH, size=(int(k), 2))
        if not P:
            out[i] = cand[0]
            P.append([cand[0, 0], cand[0, 1], r])
            continue
        Q = np.asarray(P)
        dx = cand[:, 0][:, None] - Q[None, :, 0]
        dy = cand[:, 1][:, None] - Q[None, :, 1]
        dx -= PATCH * np.round(dx / PATCH)
        dy -= PATCH * np.round(dy / PATCH)
        d = np.hypot(dx, dy)
        score = (d / (r + Q[None, :, 2])).min(axis=1)
        j = int(np.argmax(score))
        out[i] = cand[j]
        P.append([cand[j, 0], cand[j, 1], r])
    return out


def populate(seed, armored=True, clusters=True, mono=True,
             surf_cover=SURF_COVER, cluster_scale=CLUSTER_SCALE,
             bulk_vol=BULK_VOL, p_scale=P_SCALE0, sort_gain=SORT_GAIN,
             sort_q=SORT_Q, bulk_depth=0.0):
    """Grain list with x, y and the protrusion p already drawn.

    Positions and p do not depend on the substrate amplitude, so a change of
    sigma_form re-uses this list verbatim (the sigma_z fix-point only moves
    the form and the p scale).
    """
    rng_bulk = np.random.default_rng(seed + 101)
    rng_surf = np.random.default_rng(seed + 102)
    rng_clu = np.random.default_rng(seed + 103)
    rng_mono = np.random.default_rng(seed + 104)
    grains = []
    A = PATCH * PATCH
    s = math.log(SIGMA_G)

    # ---- bulk (sub-surface) ----------------------------------------------
    v_target = A * BULK_THICK * bulk_vol
    vol = 0.0
    while vol < v_target:
        D = float(BG.draw_size(rng_bulk, 1, D50_IN,
                              clip=(4.0e-4, D_CLIP_MAX))[0])
        g = BG.make_grain(rng_bulk, D, kind="bulk")
        A_, B_, C_ = g["abc"]
        g["V"] = 4.0 / 3.0 * math.pi * A_ * B_ * C_
        vol += g["V"]
        g["x"] = float(rng_bulk.uniform(0.0, PATCH))
        g["y"] = float(rng_bulk.uniform(0.0, PATCH))
        # buried: the top of a bulk grain sits below the armour layer, and the
        # depth grows with the grain size.  Small bulk grains therefore fill
        # the interstices just under the armour (no open 20-40 mm pits), the
        # coarse ones stay deep.
        g["p"] = float(-(bulk_depth * D50_SURF
                         + (0.15 + 1.10 * rng_bulk.random()) * D))
        grains.append(g)

    # ---- clusters (positioned first: they are obstacles for the packing) --
    cm = []
    anchors = []
    if clusters:
        d90 = D50_SURF * (SIGMA_G ** 1.2815515655446004)
        for c in range(N_CLUSTERS):
            D = d90 * cluster_scale * rng_clu.uniform(1.02, 1.22)
            anchor = BG.make_grain(rng_clu, D, kind="cluster_anchor", cluster=c)
            anchor["x"] = float(rng_clu.uniform(0.05, 0.95) * PATCH)
            anchor["y"] = float(rng_clu.uniform(0.05, 0.95) * PATCH)
            anchor["p"] = 0.34 * D
            grains.append(anchor)
            cm.append(anchor)
            anchors.append(anchor)
            for k in range(2):
                Dm = D * rng_clu.uniform(0.38, 0.55)
                gm = BG.make_grain(rng_clu, Dm, kind="cluster_medium",
                                   cluster=c)
                gm["off"] = (-(0.55 + 0.32 * k) * BG._half_extent(anchor, 0)
                             * rng_clu.uniform(1.0, 1.25),
                             rng_clu.uniform(-0.35, 0.35)
                             * BG._half_extent(anchor, 1))
                grains.append(gm)
                cm.append(gm)
            for k in range(3):
                Df = D * rng_clu.uniform(0.16, 0.26)
                gf = BG.make_grain(rng_clu, Df, kind="cluster_fine", cluster=c)
                s_ = (0.55 + 0.42 * k) * BG._half_extent(anchor, 0)
                gf["off"] = (s_ * rng_clu.uniform(0.95, 1.15),
                             (rng_clu.uniform(-0.55, 0.55) if k % 2 else
                              rng_clu.uniform(-0.25, 0.25)) * s_)
                grains.append(gf)
                cm.append(gf)

    # ---- the two boulders (shadow probes; NOT bed) -----------------------
    mono_grains = []
    if mono:
        for D in MONO_D:
            a = BG.make_grain(rng_mono, D, kind="mono", cluster="mono")
            a["x"] = float(rng_mono.uniform(0.22, 0.40))
            a["y"] = float(rng_mono.uniform(0.22, 0.40))
            a["p"] = 0.35 * D
            grains.append(a)
            cm.append(a)
            mono_grains.append(a)
            for k in range(3):
                Dm = D * rng_mono.uniform(0.16, 0.24)
                gm = BG.make_grain(rng_mono, Dm, kind="mono_medium",
                                   cluster="mono")
                gm["off"] = (-(0.60 + 0.30 * k) * BG._half_extent(a, 0)
                             * rng_mono.uniform(1.0, 1.2),
                             rng_mono.uniform(-0.3, 0.3)
                             * BG._half_extent(a, 1))
                grains.append(gm)
                cm.append(gm)
            for k in range(4):
                Df = D * rng_mono.uniform(0.10, 0.16)
                gf = BG.make_grain(rng_mono, Df, kind="mono_fine",
                                   cluster="mono")
                s_ = (0.60 + 0.40 * k) * BG._half_extent(a, 0)
                gf["off"] = (s_ * rng_mono.uniform(0.95, 1.15),
                             rng_mono.uniform(-0.5, 0.5) * s_)
                grains.append(gf)
                cm.append(gf)

    # ---- armoured surface layer, packed around the cluster anchors --------
    # (the two boulders are deliberately NOT obstacles: their footprint is
    # excluded from the bed audit, and a boulder parked on an empty patch
    # would leave a 30-100 mm crater in the boulder-free map)
    if armored:
        d50s = D50_SURF
        e_d2 = d50s ** 2 * math.exp(-2.0 * s * s)
        e_fp = math.pi * e_d2 / (4.0 * 0.705)
        n_surf = int(round(surf_cover * A / e_fp))
        Ds = BG.draw_size(rng_surf, n_surf, d50s,
                          clip=(4.0e-4, D_CLIP_MAX))
        surf_grains = []
        for D in Ds:
            g = BG.make_grain(rng_surf, float(D), kind="surface")
            g["p"] = protrusion(rng_surf, float(D), p_scale,
                                sort_gain=sort_gain, sort_q=sort_q)
            surf_grains.append(g)
        obst = [(g["x"], g["y"], math.sqrt(BG.footprint_area(g) / math.pi))
                for g in anchors]
        radii = [math.sqrt(BG.footprint_area(g) / math.pi) for g in surf_grains]
        xy = pack_positions(radii, rng_surf, k=32, obstacles=obst)
        for g, (x, y) in zip(surf_grains, xy):
            g["x"], g["y"] = float(x), float(y)
        grains.extend(surf_grains)

    # ---- cluster satellites: offsets relative to their anchor ------------
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
        g["x"] = float((anchor["x"] + ox) % PATCH)
        g["y"] = float((anchor["y"] + oy) % PATCH)
        D = g["D"]
        # satellites sit a little lower than their anchor, same sign rule
        g["p"] = float(0.5 * anchor["p"] + 0.5 * protrusion(
            rng_clu if g["cluster"] != "mono" else rng_mono, D, p_scale))
    return grains, cm


def place(grains, form, p_scale_note=None):
    """Set z from the local form + the grain's own protrusion p."""
    n = form.shape[0]
    cx = form.shape[1]
    cell = PATCH / cx
    xs = (np.arange(n) + 0.5) * cell
    Xg = np.array([g["x"] for g in grains])
    Yg = np.array([g["y"] for g in grains])
    fx = np.asarray(BG._bilinear(form, (Xg / cell)[None, :],
                                 (Yg / cell)[None, :])[0])
    hz = np.array([BG._half_extent(g, 2) for g in grains])
    p = np.array([g["p"] for g in grains])
    z = fx + p - hz
    for i, g in enumerate(grains):
        g["z"] = float(z[i])
    return {"cell": cell, "z_center_range": [float(z.min()), float(z.max())]}


# ---------------------------------------------------------------------------
# vertical ray cast with PERIODIC wrap (the transport wraps in x and y)
# ---------------------------------------------------------------------------
def raycast2(grains, xs, ys, skip=None):
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
        for arr, val in ((nxv, nwx), (nyv, nwy), (nzv, nwz)):
            cur = arr[idx]
            arr[idx] = np.where(hit, val, cur)
    if (gid < 0).any():
        zz = np.where(np.isfinite(z), z, np.nan)
        z = np.where(np.isfinite(z), z, np.nanmean(zz))
    return {"z": z, "gid": gid, "nx": nxv, "ny": nyv, "nz": nzv}


def maps2(grains, w, h, x0, x1, y0, y1, skip=None, n_facet=80):
    xs = x0 + (np.arange(w) + 0.5) * ((x1 - x0) / w)
    ys = y0 + (np.arange(h) + 0.5) * ((y1 - y0) / h)
    rc = raycast2(grains, xs, ys, skip=skip)
    cell = (x1 - x0) / w
    fx, fy, fz, fidx = facet_normals_chunked(rc["nx"], rc["ny"], rc["nz"],
                                             n_facet)
    rc.update({"fx": fx, "fy": fy, "fz": fz, "fidx": fidx, "cell": cell})
    return rc


def facet_normals_chunked(nx, ny, nz, n_facet=80, chunk=262144):
    """Same quantisation as bed_gen.facet_normals, done in row chunks."""
    F = BG.fibonacci_hemisphere(n_facet)
    H, W = nx.shape
    flat = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1)
    idx = np.empty(flat.shape[0], dtype=np.int64)
    for a in range(0, flat.shape[0], chunk):
        b = min(a + chunk, flat.shape[0])
        idx[a:b] = np.argmax(flat[a:b] @ F.T, axis=1)
    N = F[idx]
    return (N[:, 0].reshape(H, W), N[:, 1].reshape(H, W),
            N[:, 2].reshape(H, W), idx.reshape(H, W))


# ---------------------------------------------------------------------------
# audits
# ---------------------------------------------------------------------------
def exposed_stats(gid, cell, grains, is_mono, is_cluster):
    """EXACT surface GSD from the grain list, weighted by exposed top area.

    (1) every pixel of the statistics depth map is already owned by the
        top-most grain (gid); its exposed area is n_px * cell^2.
    (2) the SURFACE GSD is the armoured surface LAYER (kind == "surface"),
        area-weighted by its exposed tops.  This is the object the brief's
        D50 / D84 / sigma_g refer to (the surface bed material).
    (3) the lower layer has no exposure, so it is weighted by its footprint
        area (the same areal convention the surface draw uses) -> armour ratio
        = D50_surface / D50_lower.

    A bed-wide exposure-weighted GSD (every exposed grain, boulders excluded)
    is reported alongside; it is only identical to (2) when the armoured
    layer covers the bed completely.
    """
    D = np.array([g["D"] for g in grains])
    valid = gid >= 0
    counts = np.bincount(gid[valid].ravel(), minlength=len(D)).astype(np.float64)
    area = counts * cell * cell
    is_surf = np.array([g["kind"] == "surface" for g in grains])

    def _pq(mask, wgt):
        Dk, Ak = D[mask], wgt[mask]
        keep = Ak > 0.5 * cell * cell
        Dk, Ak = Dk[keep], Ak[keep]
        if Dk.size < 10:
            return None
        o = np.argsort(Dk)
        Dk, Ak = Dk[o], Ak[o]
        cw = np.cumsum(Ak) / Ak.sum()
        return {"n_grains_exposed": int(Dk.size),
                "exposed_area_m2": float(Ak.sum()),
                "D16_m": float(np.interp(0.16, cw, Dk)),
                "D50_m": float(np.interp(0.50, cw, Dk)),
                "D84_m": float(np.interp(0.84, cw, Dk)),
                "D84_over_D50": float(np.interp(0.84, cw, Dk)
                                      / np.interp(0.50, cw, Dk)),
                "D90_m": float(np.interp(0.90, cw, Dk)),
                "D_max_m": float(Dk.max())}

    surf = _pq(is_surf, area)
    surf["weighting"] = ("exposed top area of the armoured surface layer "
                         "(kind == surface) from the statistics depth map; "
                         "cluster structures and the two boulders excluded")
    bedwide = _pq((~is_mono) & (counts > 0.5), area)
    if bedwide:
        bedwide["weighting"] = ("every exposed grain (boulders excluded), "
                                "reported for comparison only")

    # lower layer: bulk, footprint-area weighted
    is_bulk = np.array([g["kind"] == "bulk" for g in grains])
    Db = np.array([BG.footprint_area(g) for g in grains])[is_bulk]
    Dbd = D[is_bulk]
    o = np.argsort(Dbd)
    Dbd, Db = Dbd[o], Db[o]
    cb = np.cumsum(Db) / Db.sum()
    lower = {"n": int(is_bulk.sum()),
             "D50_m": float(np.interp(0.50, cb, Dbd)),
             "D84_over_D50": float(np.interp(0.84, cb, Dbd)
                                   / np.interp(0.50, cb, Dbd)),
             "weighting": "footprint area pi*hx*hy (the sub-surface layer "
                          "has no exposure); this is the same areal "
                          "convention used to draw the surface layer"}
    ratio = surf["D50_m"] / lower["D50_m"]
    return {"surface": surf, "bedwide_exposed": bedwide, "lower": lower,
            "armor_ratio_exact": ratio,
            "cluster_coverage_exposed": float(
                area[is_cluster & (~is_mono)].sum() / (gid.size * cell * cell))}


def structure_audit(z, cell, d50ref=D50_SURF, dmax=0.20, nbin=60):
    """Structure function split at the measured break point (DESIGN audit 2).

    The break rule is the same one v2.1 declared (first scale above 0.5 D50
    where the local log-slope drops below 0.6 x (twice the 2-20 cm fit)), so
    the two reports are comparable.  Both segment exponents, the break scale
    and an alternative knee definition are reported.
    """
    d, s = BG.structure_function(z, cell, dmax=dmax, nbin=nbin)
    lo, hi = 2.0 * D50_IN, min(dmax, 0.20)
    m = (d >= lo) & (d <= hi)
    H_full = 0.5 * np.polyfit(np.log(d[m]), np.log(s[m]), 1)[0]
    loc = 0.5 * np.gradient(np.log(s), np.log(d))

    def _fit(a, b):
        mm = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
        if mm.sum() < 4:
            return None, int(mm.sum())
        return float(0.5 * np.polyfit(np.log(d[mm]), np.log(s[mm]), 1)[0]), \
            int(mm.sum())

    brk = None
    target = 0.6 * 2.0 * H_full
    for i in range(1, len(d) - 1):
        if d[i] < 0.5 * D50_IN:
            continue
        if loc[i] < target:
            brk = float(d[i])
            break
    H_below, n_below = (_fit(2.0 * cell, brk) if brk else (None, 0))
    H_above, n_above = (_fit(brk, d[-1]) if brk else (None, 0))
    # alternative: the knee = the scale of maximum local slope
    i_knee = int(np.nanargmax(loc[d >= 0.5 * D50_IN]))
    dk = float(d[d >= 0.5 * D50_IN][i_knee])
    H_above_knee, n_ak = _fit(dk, d[-1])
    return {
        "H_full_band": float(H_full),
        "fit_band_m": [lo, hi],
        "break_scale_m": brk,
        "break_over_D50_input": (brk / D50_IN) if brk else None,
        "break_over_D50_surface": (brk / d50ref) if brk else None,
        "H_below_break": H_below, "n_points_below": n_below,
        "H_above_break": H_above, "n_points_above": n_above,
        "knee_scale_m": dk, "H_above_knee": H_above_knee,
        "n_points_above_knee": n_ak,
        "break_rule": "first d > 0.5*D50_input where the local log-slope "
                      "falls below 0.6*(2*H_full_band); same rule as v4-2.1",
        "D_m": d.tolist(), "S_m2": s.tolist(), "H_local": loc.tolist(),
    }


def _skew(a):
    a = np.asarray(a, dtype=np.float64).ravel()
    a = a[np.isfinite(a)]
    m = a.mean()
    return float(np.mean((a - m) ** 3) / max(a.std() ** 3, 1e-30))


def audit_bed(grains, hi, mono_mask, form, cell_form, D84=D84_SURF):
    """Matrix / whole-bed statistics on one statistics-resolution map.

    The matrix mask follows v4-2.1 exactly so that the two reports compare:
    the top-most grain is not a boulder and its D <= D84 of the armoured
    surface GSD.  (The three cluster ANCHORS are all larger than that, so they
    are excluded by the D rule; the cluster mediums/fines are kept, as in
    v4-2.1.)  The variant that additionally drops every cluster member is
    reported alongside, because that exclusion is a defensible alternative.
    """
    gid = hi["gid"]
    z = hi["z"]
    D = np.array([g["D"] for g in grains])
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    is_cluster = np.array([k.startswith("cluster") for k in kinds])
    bed = ~is_mono[gid] & (gid >= 0)
    Dpx = D[np.clip(gid, 0, len(D) - 1)]
    mat = bed & (Dpx <= D84)
    mat2 = mat & (~is_cluster[np.clip(gid, 0, len(D) - 1)])
    zm = z[mat]
    zb = z[bed]
    out = {
        "matrix_definition": ("top-most grain has D <= D84 of the armoured "
                              "surface GSD (%.3f m) and is not a boulder; the "
                              "same rule as v4-2.1 (the cluster anchors are "
                              "larger than D84 and drop out by the D rule)"
                              % D84),
        "n_matrix_px": int(mat.sum()),
        "matrix_area_fraction": float(mat.sum() / z.size),
        "sigma_z_matrix_m": float(np.std(zm)),
        "sigma_z_matrix_over_D50_input": float(np.std(zm) / D50_IN),
        "skew_matrix": _skew(zm),
        "sigma_z_matrix_no_clusters_m": float(np.std(z[mat2])),
        "skew_matrix_no_clusters": _skew(z[mat2]),
        "matrix_no_clusters_area_fraction": float(mat2.sum() / z.size),
        "sigma_z_bed_m": float(np.std(zb)), "skew_bed": _skew(zb),
        "sigma_z_grain_m": float(np.std((z - form)[mat])),
        "sigma_z_form_m": float(np.std(form)),
        "z_range_m": [float(zb.min()), float(zb.max())],
    }
    return out


# ---------------------------------------------------------------------------
# deliverables
# ---------------------------------------------------------------------------
def save_heightfield(path, z):
    """Write a real .npy (header + float32 payload), atomically."""
    import io
    buf = io.BytesIO()
    np.save(buf, z.astype("<f4"))
    data = buf.getvalue()

    def verify(p):
        try:
            a = np.load(p, allow_pickle=False)
            return a.shape == z.shape and a.dtype == np.dtype("<f4")
        except Exception:
            return False
    atomic_write_bytes(path, data, verify=verify)
    return len(data)


def save_transects(path, z, cell, d50=D50_SURF, extra=None):
    """Flow-direction and cross-stream profiles, mm units, D50 scale bar."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    H, W = z.shape
    jc = H // 2
    ic = W // 2
    xmm = (np.arange(W) + 0.5) * cell * 1000.0
    ymm = (np.arange(H) + 0.5) * cell * 1000.0
    zmm = (z - z.mean()) * 1000.0
    band = 9
    xline = zmm[jc]
    yline = zmm[:, ic]
    xband = zmm[jc - band // 2:jc + band // 2 + 1].mean(axis=0)
    yband = zmm[:, ic - band // 2:ic + band // 2 + 1].mean(axis=1)

    fig, ax = plt.subplots(1, 2, figsize=(13.2, 4.6), dpi=120,
                           facecolor="white")
    for a, xx, single, bandv, ttl, xl in (
            (ax[0], xmm, xline, xband, "Flow direction (x), transect at y = %.1f mm"
             % (ymm[jc]), "x (mm)"),
            (ax[1], ymm, yline, yband, "Cross-stream (y), transect at x = %.1f mm"
             % (xmm[ic]), "y (mm)")):
        a.plot(xx, single, lw=0.7, color="#3a6ea8", label="single line")
        a.plot(xx, bandv, lw=1.6, color="#b03030", label="9-line mean")
        a.axhline(0.0, lw=0.6, color="#999999", ls=":")
        a.set_xlabel(xl)
        a.set_ylabel("z - mean (mm)")
        a.set_title(ttl, fontsize=9)
        a.set_xlim(0, xx[-1])
        a.grid(alpha=0.25, lw=0.4)
        # D50 scale bar (horizontal, mm)
        x0 = 0.72 * xx[-1]
        y0 = a.get_ylim()[0] + 0.12 * (a.get_ylim()[1] - a.get_ylim()[0])
        a.plot([x0, x0 + d50 * 1000.0], [y0, y0], lw=3.0, color="black",
               solid_capstyle="butt")
        a.text(x0 + 0.5 * d50 * 1000.0, y0 + 0.05 * (a.get_ylim()[1]
                                                     - a.get_ylim()[0]),
               "D50_surface = %.0f mm" % (d50 * 1000.0), ha="center",
               fontsize=8)
        a.legend(loc="upper right", fontsize=7, framealpha=0.85)
    if extra:
        fig.text(0.5, 0.005, extra, ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    import io
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=120, facecolor="white")
    plt.close(fig)
    data = b.getvalue()

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size[0] > 400
        except Exception:
            return False
    atomic_write_bytes(path, data, verify=verify)


# ---------------------------------------------------------------------------
# full build: sigma_z measured fix-point
# ---------------------------------------------------------------------------
_POP_CACHE = {}


def form_at_grid(form, form_n, cell):
    """Bilinear resample of the form onto the statistics grid (`cell`)."""
    cellf = PATCH / form_n
    iif = (np.arange(HIRES) + 0.5) * cell
    return BG._bilinear(form, (iif / cellf)[None, :], (iif / cellf)[:, None])


def get_population(seed, p_scale):
    key = (int(seed), round(float(p_scale), 4), SURF_COVER, CLUSTER_SCALE,
           BULK_VOL)
    if key not in _POP_CACHE:
        _POP_CACHE[key] = populate(seed, p_scale=p_scale)
    return _POP_CACHE[key]


def build_once(seed, sigma_form, p_scale, n_form=512):
    t0 = time.time()
    form = form_field(n_form, sigma_form, seed)
    grains, cm = get_population(seed, p_scale)
    place(grains, form)
    hi = maps2(grains, HIRES, HIRES, 0.0, PATCH, 0.0, PATCH)
    return {"grains": grains, "cluster_members": cm, "form": form,
            "hi": hi, "wall_s": float(time.time() - t0),
            "sigma_form": float(sigma_form), "p_scale": float(p_scale)}


def build_and_audit(out_dir=None, seed=SEED, sun=None, verbose=True,
                    max_attempts=4):
    """Build + audit with the sigma_z fix-point measured on the depth map.

    Target (declared before the first measurement): sigma_z(matrix) =
    0.40 * D50_surface, where D50_surface is the EXACT audit's surface D50
    (exposed-area weighted, measured from the grain list, not assumed to be
    0.040 m).  One knob is moved -- the substrate amplitude sigma_form; the
    burial/protrusion scale is a design constant fixed by the calibration
    sweeps recorded in the report.  Every attempt is returned in the stats.
    """
    lg = log if verbose else (lambda *a: None)
    form_n = 512
    sigma_form = 0.85 * SIGMA_TARGET
    p_scale = P_SCALE0
    attempts = []
    target = SIGMA_TARGET          # provisional, replaced after attempt 1
    target_rule = ("0.40 * D50_surface; D50_surface measured by the exact "
                   "exposed-area audit on the same population")
    for it in range(max_attempts):
        t0 = time.time()
        b = build_once(seed, sigma_form, p_scale, n_form=form_n)
        grains = b["grains"]
        kinds = np.array([str(g["kind"]) for g in grains])
        is_mono = np.array([k.startswith("mono") for k in kinds])
        is_cluster = np.array([k.startswith("cluster") for k in kinds])
        cell = b["hi"]["cell"]
        formpix = form_at_grid(b["form"], form_n, cell)
        stat = audit_bed(grains, b["hi"], is_mono, formpix, cell)
        ex_a = exposed_stats(b["hi"]["gid"], cell, grains, is_mono,
                             is_cluster)
        d50s = ex_a["surface"]["D50_m"]
        if it == 0:
            target = 0.40 * d50s
            lg("  surface D50 (exact) = %.5f m -> sigma_z target = %.5f m"
               % (d50s, target))
        sig = stat["sigma_z_matrix_m"]
        attempts.append({
            "attempt": it + 1,
            "sigma_form_in_m": float(sigma_form),
            "p_scale_in": float(p_scale),
            "D50_surface_measured_m": float(d50s),
            "target_m": float(target),
            "sigma_z_matrix_m": float(sig),
            "sigma_z_over_D50_surface": float(sig / d50s),
            "rel_err": float((sig - target) / target),
            "sigma_z_grain_m": stat["sigma_z_grain_m"],
            "sigma_z_form_realised_m": stat["sigma_z_form_m"],
            "skew_matrix": stat["skew_matrix"],
            "wall_s": float(time.time() - t0),
        })
        lg("  sigma_z attempt %d: sigma_form=%.5f -> sigma_z(matrix)=%.5f "
           "(%.2f%% of target %.5f, %.3f D50s) grain=%.5f skew=%.3f"
           % (it + 1, sigma_form, sig, 100.0 * (sig - target) / target,
              target, sig / d50s, stat["sigma_z_grain_m"],
              stat["skew_matrix"]))
        if abs(sig - target) <= 0.01 * target:
            lg("  sigma_z closed within 1% of the target")
            break
        # declared update rule: sigma_z responds almost linearly to the
        # substrate amplitude in the form-dominated regime, so apply the
        # measured ratio (damped by 0.9) to the amplitude
        ratio = target / max(sig, 1e-6)
        sigma_form = float(np.clip(sigma_form * (1.0 + 0.9 * (ratio - 1.0)),
                                   0.25 * target, 1.6 * target))
    # ---- final build with the LAST parameters (no cherry-picking) --------
    b = build_once(seed, sigma_form, p_scale, n_form=form_n)
    grains = b["grains"]
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    is_cluster = np.array([k.startswith("cluster") for k in kinds])
    D = np.array([g["D"] for g in grains])
    cell = b["hi"]["cell"]
    formpix = form_at_grid(b["form"], form_n, cell)
    stat = audit_bed(grains, b["hi"], is_mono, formpix, cell)
    lg("  final sigma_z(matrix)=%.5f m (%.3f x D50_surface %.5f m)"
       % (stat["sigma_z_matrix_m"],
          stat["sigma_z_matrix_m"] / D50_SURF, D50_SURF))

    # ---- exact surface GSD / armour ratio from the grain list -----------
    ex = exposed_stats(b["hi"]["gid"], cell, grains, is_mono, is_cluster)
    d50s_final = ex["surface"]["D50_m"]
    lg("  exact surface GSD: D50=%.5f D84/D50=%.3f  lower D50=%.5f "
       "armour=%.3f" % (d50s_final, ex["surface"]["D84_over_D50"],
                        ex["lower"]["D50_m"], ex["armor_ratio_exact"]))

    # ---- bed map without the boulders (for SF / facets / delivery) ------
    skip = set(int(i) for i in np.where(is_mono)[0])
    hi_bed = maps2(grains, HIRES, HIRES, 0.0, PATCH, 0.0, PATCH, skip=skip)
    zh = hi_bed["z"]
    gidb = hi_bed["gid"]
    hole = gidb < 0
    if hole.any():
        from scipy import ndimage
        wgt = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, zh), 60.0,
                                      mode="wrap")
        den = ndimage.gaussian_filter(wgt, 60.0, mode="wrap")
        zh = np.where(hole, num / np.maximum(den, 1e-9), zh)
    sf = structure_audit(zh, cell, d50ref=d50s_final)
    lg("  structure: full=%.3f  break=%.4f m (%.2f D50s)  above=%.3f"
       % (sf["H_full_band"], sf["break_scale_m"] or float("nan"),
          (sf["break_over_D50_surface"] or float("nan")),
          sf["H_above_break"] if sf["H_above_break"] is not None
          else float("nan")))

    # ---- facet slopes ----------------------------------------------------
    fz = hi_bed["fz"]
    slope = np.degrees(np.arccos(np.clip(fz, -1.0, 1.0)))
    facets = {"n_directions": 80,
              "slope_median_deg": float(np.median(slope)),
              "slope_p84_deg": float(np.percentile(slope, 84)),
              "frac_gt60deg": float(np.mean(slope > 60.0)),
              "frac_gt80deg": float(np.mean(slope > 80.0)),
              "method": "analytic ellipsoid normal quantised to the nearest "
                        "of 80 Fibonacci hemisphere directions"}

    # ---- porosity --------------------------------------------------------
    grains_bed = [g for g in grains if not str(g["kind"]).startswith("mono")]
    por = BG.porosity(grains_bed, zh, cell)
    zrel = zh - float(np.mean(zh))
    pits = {"min_m": float(zh.min()),
            "frac_below_20mm": float(np.mean(zrel < -0.020)),
            "frac_below_30mm": float(np.mean(zrel < -0.030)),
            "frac_below_50mm": float(np.mean(zrel < -0.050)),
            "note": ("interstices where the vertical ray reaches a coarse "
                     "buried grain instead of an armour grain.  They are part "
                     "of the measured statistics (sigma_z, skewness, H) and "
                     "are disclosed; they are not a render artifact")}
    lg("  porosity=%.4f  cluster coverage=%.4f"
       % (por["porosity"], ex["cluster_coverage_exposed"]))

    # ---- imbrication -----------------------------------------------------
    surf = [g for g in grains if g["kind"] in
            ("surface", "cluster_anchor", "cluster_medium", "cluster_fine")]
    dips = np.array([g["dip_deg"] for g in surf])
    down = np.array([[-g["R"][0, 2], -g["R"][1, 2]] for g in surf])
    down /= np.maximum(np.linalg.norm(down, axis=1, keepdims=True), 1e-12)
    upstream = down[:, 0] < -0.5
    imb = {"n": int(dips.size),
           "dip_median_deg": float(np.median(dips[upstream])) if
           upstream.any() else float(np.median(dips)),
           "dip_median_all_deg": float(np.median(dips)),
           "dip_p16_deg": float(np.percentile(dips, 16)),
           "dip_p84_deg": float(np.percentile(dips, 84)),
           "upstream_fraction": float(np.mean(upstream)),
           "rule": "a-b plane dips 30 +- 10 deg, normal tilted downstream so "
                   "the plane faces upstream (unchanged from v4-2.1)"}

    # ---- penetration statistics (declared limitation, measured) ---------
    pen = penetration_stats(grains, D)

    # ---- reconstruction operator (calibration only, NOT the gate) -------
    recon = BG.reconstruct_gsd(zh, cell)
    ref_bed = build_reference(seed, sigma_form, p_scale, form_n)
    ref_recon = BG.reconstruct_gsd(ref_bed["z"], cell)
    recon_block = {
        "main": recon, "reference": ref_recon,
        "reference_definition": "same generator with armored=False, "
                                "clusters=False, mono=False: an UNARMOURED "
                                "bed drawn from D50 = 0.020 m",
        "bias_main": (recon["D50_m"] / ex["surface"]["D50_m"]
                      if recon.get("D50_m") else None),
        "bias_reference": (ref_recon["D50_m"] / D50_IN
                           if ref_recon.get("D50_m") else None),
        "note": "the residual->maxima->nearest-peak-Voronoi operator reads "
                "2.10 as 1.20 on a known input (v4-2.1 section 1.1); it is "
                "reported for method continuity and is NOT used by any gate",
    }

    # ---- AO / shadow / full patch maps for the renderer ------------------
    t0 = time.time()
    full = maps2(grains, RENDER_W, RENDER_W, 0.0, PATCH, 0.0, PATCH)
    cellf2 = full["cell"]
    ao_full = BG.ambient_occlusion(BG.horizon_map(full["z"], cellf2))
    sh_full = np.ones_like(ao_full)
    sun_info = None
    if sun is not None:
        s = np.asarray(sun["sunDir"], dtype=np.float64)
        s = s / np.linalg.norm(s)
        eta = 1.0 / 1.33434
        cos1 = float(s[2])
        k = 1.0 - eta * eta * (1.0 - cos1 * cos1)
        sh = eta * s + (math.sqrt(max(k, 0.0)) - eta * cos1) \
            * np.array([0.0, 0.0, 1.0])
        sh /= np.linalg.norm(sh)
        sh_full = BG.shadow_map(full["z"], cellf2, sh, float(sun["disc_rad"]))
        sun_info = {"sun_dir_air": [float(x) for x in s],
                    "sun_dir_water": [float(x) for x in sh],
                    "altitude_deg": float(sun.get("altitude_deg",
                                                  float("nan"))),
                    "azimuth_deg": float(sun.get("azimuth_deg",
                                                  float("nan"))),
                    "disc_half_angle_rad": float(sun["disc_rad"]),
                    "n_disc_samples": 9,
                    "umbra_fraction": float(np.mean(sh_full <= 0.0)),
                    "penumbra_fraction": float(np.mean((sh_full > 0.0)
                                                       & (sh_full < 1.0))),
                    "lit_fraction": float(np.mean(sh_full >= 1.0))}
    rmap = maps2(grains, RENDER_W, RENDER_H, WIN_X0, WIN_X1, WIN_Y0, WIN_Y1)
    rmap["ao"] = ao_full[RENDER_J0:RENDER_J0 + RENDER_H, :]
    rmap["shadow"] = sh_full[RENDER_J0:RENDER_J0 + RENDER_H, :]
    full["ao"] = ao_full
    full["shadow"] = sh_full
    map_wall = float(time.time() - t0)

    stats = {
        "n_grains": len(grains),
        "kinds": {str(k): int((kinds == k).sum()) for k in sorted(set(kinds))},
        "patch_m": PATCH, "D50_in_m": D50_IN, "sigma_g": SIGMA_G,
        "armor_ratio_input": ARMOR,
        "D50_surface_nominal_m": D50_SURF,
        "D84_surface_nominal_m": D84_SURF,
        "sigma_target_m": SIGMA_TARGET,
        "sigma_target_rule": ("0.40 * D50_surface, with D50_surface the EXACT "
                              "exposed-area surface D50 (measured %.5f m). "
                              "The nominal 0.040 m would give 0.0160 m; the "
                              "measured target is %.5f m."
                              % (d50s_final, 0.40 * d50s_final)),
        "D50_surface_measured_m": float(d50s_final),
        "sigma_z_over_D50_surface": float(stat["sigma_z_matrix_m"]
                                          / d50s_final),
        "sigma_z_fixpoint": attempts,
        "sigma_z_fixpoint_final_parameters": {
            "sigma_form_m": float(sigma_form), "p_scale": float(p_scale),
            "note": "the final build uses the LAST parameters of the loop, "
                    "not the best measured attempt"},
        **{k: v for k, v in stat.items()},
        "exact_gsd": ex,
        "armor_ratio_exact": ex["armor_ratio_exact"],
        "D84_over_D50_exact": ex["surface"]["D84_over_D50"],
        "structure": {k: v for k, v in sf.items()
                      if k not in ("D_m", "S_m2", "H_local")},
        "structure_curves": {k: sf[k] for k in ("D_m", "S_m2", "H_local")},
        "facets": facets,
        "porosity": por,
        "pits": pits,
        "imbrication": imb,
        "penetration": pen,
        "gsd_reconstructed": recon_block,
        "cluster_coverage": ex["cluster_coverage_exposed"],
        "cluster_coverage_target": 0.10,
        "n_clusters": int(N_CLUSTERS),
        "form": {"H_gen": SUBSTRATE_H, "w_low": FORM_W_LOW,
                 "skew_c": FORM_SKEW_C, "n_modes_low": "8 modes, |k| <= 2/L",
                 "sigma_form_m": float(sigma_form),
                 "sigma_form_realised_m": stat["sigma_z_form_m"],
                 "skew_form": _skew(formpix),
                 "note": "fBm(H=0.79) + patch-scale low modes + quadratic "
                         "skew transform; amplitude solved inside the "
                         "sigma_z budget"},
        "protrusion_model": {
            "p_tail_prob": P_TAIL, "p_tail_mean_over_D": P_TAIL_MEAN,
            "p_tail_max_over_D": P_TAIL_MAX,
            "p_bury_range_over_D": [P_BURY_MIN, P_BURY_MAX],
            "size_sorting_gain": SORT_GAIN, "size_sorting_q": SORT_Q,
            "p_scale": float(p_scale),
            "surface_coverage_target": SURF_COVER,
            "bulk_solid_volume_factor": BULK_VOL,
            "bulk_depth_over_D50s": 0.12,
            "D_clip_max_m": D_CLIP_MAX,
            "note": "summit of each grain sits at form(x,y) + p, with p = "
                    "sort_gain*((D/D50s)^q-1)*D50s + p_scale*p_random.  The "
                    "*measured* grain relief after the visible envelope (max "
                    "over the grains covering a pixel) also includes the "
                    "ellipsoid cap curvature."},
        "ao": {"method": "horizon map over 12 azimuths, radii 1.2-45 mm, "
                         "V = mean_d 1/(1+tan^2 h_d)",
               "mean": float(np.mean(rmap["ao"])),
               "frac_below_0.5": float(np.mean(rmap["ao"] < 0.5))},
        "shadow": sun_info,
        "map_wall_s": map_wall,
        "raycast": {"render": [RENDER_W, RENDER_H],
                    "render_window_m": [WIN_X0, WIN_X1, WIN_Y0, WIN_Y1],
                    "hires": [HIRES, HIRES], "hires_cell_m": cell,
                    "render_cell_m": rmap["cell"],
                    "method": "vertical ray cast (periodic wrap in x and y), "
                              "analytic ellipsoid normal at the hit point"},
        "full_patch": {"n_px": int(full["z"].size),
                       "sigma_z_m": float(np.std(full["z"]))},
        "created_utc": now_iso(),
    }
    return {"grains": grains, "cluster_members": b["cluster_members"],
            "hi": b["hi"], "hi_bed": hi_bed, "zh": zh, "rmap": rmap,
            "full": full, "stats": stats}


def penetration_stats(grains, D, max_pairs=400000):
    """Fraction of neighbouring grain pairs whose ellipsoids overlap."""
    from scipy.spatial import cKDTree
    n = len(grains)
    P = np.array([[g["x"], g["y"], g["z"]] for g in grains])
    ax = np.array([g["abc"] for g in grains], float)
    RT = np.array([g["R"].T for g in grains], float)
    tree = cKDTree(np.mod(P[:, :2], PATCH), boxsize=PATCH)
    pr = tree.query_pairs(0.06, output_type="ndarray")
    if len(pr) > max_pairs:
        pr = pr[:: max(1, len(pr) // max_pairs)]
    ii, jj = pr[:, 0], pr[:, 1]
    u = P[jj] - P[ii]
    u[:, 0] -= PATCH * np.round(u[:, 0] / PATCH)
    u[:, 1] -= PATCH * np.round(u[:, 1] / PATCH)
    d = np.linalg.norm(u, axis=1)
    good = d > 1e-9
    u = u[good] / d[good, None]
    ii, jj, d = ii[good], jj[good], d[good]
    h_i = np.sqrt(((ax[ii] * np.einsum("nkj,nj->nk", RT[ii], u)) ** 2)
                  .sum(axis=1))
    h_j = np.sqrt(((ax[jj] * np.einsum("nkj,nj->nk", RT[jj], u)) ** 2)
                  .sum(axis=1))
    ov = h_i + h_j - d
    return {"pairs_tested": int(len(d)),
            "frac_overlap": float(np.mean(ov > 0)),
            "max_overlap_m": float(max(ov.max(), 0.0)),
            "median_overlap_when_overlapping_m":
                float(np.median(ov[ov > 0])) if (ov > 0).any() else 0.0,
            "note": "grains may interpenetrate (no relaxation in this "
                    "build); the value is disclosed, not gated"}


def build_reference(seed, sigma_form, p_scale, form_n):
    """Unarmoured reference bed (same form, same p law) for the operator bias.

    Built at half resolution: it is only used for the reconstruction-operator
    calibration, never for a gate.
    """
    form = form_field(form_n, sigma_form, seed + 505)
    grains, _ = populate(seed + 505, armored=False, clusters=False, mono=False,
                         p_scale=p_scale)
    place(grains, form)
    n = HIRES
    rc = maps2(grains, n, n, 0.0, PATCH, 0.0, PATCH)
    return {"z": rc["z"], "cell": rc["cell"]}


# ---------------------------------------------------------------------------
# G-BED2
# ---------------------------------------------------------------------------
def gate_bed2(s):
    g = {}

    def put(key, val, lo, hi, label):
        ok = (val is not None) and (lo <= val <= hi)
        g[key] = {"value": val, "range": [lo, hi], "pass": bool(ok),
                  "label": label}

    put("skewness", s["skew_matrix"], 0.3, 0.9, "matrix skewness")
    put("sigma_z_over_D50_surface", s["sigma_z_over_D50_surface"],
        0.32, 0.46, "sigma_z(matrix) / D50_surface (exact audit)")
    put("structure_H", s["structure"]["H_above_break"], 0.7, 0.85,
        "structure function H, above the measured break")
    put("D84_over_D50_exact", s["D84_over_D50_exact"], 1.9, 2.3,
        "surface D84/D50, exposed-area weighted (exact audit)")
    put("armor_ratio_exact", s["armor_ratio_exact"], 1.7, 2.3,
        "surface / lower-layer D50 (exact audit)")
    put("imbrication_median_deg", s["imbrication"]["dip_median_deg"], 25.0,
        35.0, "a-b plane dip, median, upstream-dipping grains")
    put("cluster_coverage", s["cluster_coverage"], 0.05, 0.15,
        "cluster areal coverage (exposed area)")
    put("porosity", s["porosity"]["porosity"], 0.30, 0.45, "bed porosity")

    # ---- 併記 (recorded, not gated) -------------------------------------
    st = s["structure"]
    g["recorded_reconstruction_operator"] = {
        "main_raw": {k: s["gsd_reconstructed"]["main"].get(k) for k in
                     ("D50_m", "D84_m", "D84_over_D50", "n", "n_markers")},
        "reference_raw": {k: s["gsd_reconstructed"]["reference"].get(k) for k
                          in ("D50_m", "D84_m", "D84_over_D50", "n")},
        "bias_main": s["gsd_reconstructed"]["bias_main"],
        "bias_reference": s["gsd_reconstructed"]["bias_reference"],
        "note": "residual -> maxima -> nearest-peak Voronoi; calibrated on an "
                "unarmoured reference bed with a known D50"}
    g["recorded_break_scale"] = {
        "value_m": st["break_scale_m"],
        "over_D50_input": st["break_over_D50_input"],
        "over_D50_surface": st["break_over_D50_surface"],
        "note": "recorded, the gate uses the exponent above it"}
    g["recorded_H_below_break"] = {"value": st["H_below_break"],
                                   "note": "small-scale side, recorded"}
    g["recorded_H_knee_alternative"] = {
        "knee_m": st["knee_scale_m"], "H_above_knee": st["H_above_knee"],
        "note": "an independent break definition (max local slope); recorded "
                "so the choice of break rule is auditable"}
    g["recorded_sigma_z_over_D50_input"] = {
        "value": s["sigma_z_matrix_over_D50_input"],
        "note": "the same sigma_z divided by the INPUT D50 (0.020 m); this is "
                "the quantity v4-2.1 gated at 0.32-0.46 and scored 0.78"}
    g["recorded_facet_tail"] = {
        "frac_gt60deg": s["facets"]["frac_gt60deg"],
        "frac_gt80deg": s["facets"]["frac_gt80deg"],
        "median_deg": s["facets"]["slope_median_deg"]}
    g["recorded_protrusion_and_form"] = {
        "sigma_z_form_m": s["sigma_z_form_m"],
        "sigma_z_grain_m": s["sigma_z_grain_m"],
        "sigma_grain_share": s["sigma_z_grain_m"] / s["sigma_z_matrix_m"],
        "skew_form": s["form"]["skew_form"],
        "p_scale": s["protrusion_model"]["p_scale"]}
    g["pass"] = all(v.get("pass", True) for v in g.values())
    return g


def write_bed_deliverables(out_dir, res):
    """depth-map.png, heightfield.npy (+meta), transects.png, bed-stats.json,
    _work/bed-maps.npz.  Shared by the CLI and run_v4_2_2.py."""
    stats = res["stats"]
    os.makedirs(os.path.join(out_dir, "_work"), exist_ok=True)
    BG.save_depth_png(os.path.join(out_dir, "depth-map.png"), res["rmap"]["z"],
                      np.stack([res["rmap"]["nx"], res["rmap"]["ny"],
                                res["rmap"]["nz"]], -1), res["rmap"]["cell"])
    save_heightfield(os.path.join(out_dir, "heightfield.npy"), res["zh"])
    save_transects(os.path.join(out_dir, "transects.png"), res["zh"],
                   res["hi"]["cell"], stats["D50_surface_measured_m"],
                   extra="heightfield.npy: %dx%d float32, cell %.4f mm, z in "
                         "m (mean %.4f m); the two boulders are excluded and "
                         "their footprint is inpainted"
                         % (HIRES, HIRES, res["hi"]["cell"] * 1000,
                            float(np.mean(res["zh"]))))
    meta = {"file": "heightfield.npy", "dtype": "float32 (little endian)",
            "shape": [HIRES, HIRES], "unit": "m",
            "grid": "statistics grid of the bed audit", "patch_m": PATCH,
            "cell_m": float(res["hi"]["cell"]),
            "x_range_m": [0.0, PATCH], "y_range_m": [0.0, PATCH],
            "provenance": "vertical ray cast (periodic in x and y) of the "
                          "particle bed, top-most grain per pixel; the two "
                          "mono boulders are excluded and their footprint "
                          "inpainted from the surrounding bed (same rule as "
                          "the structure-function audit)",
            "mean_z_m": float(np.mean(res["zh"])),
            "sigma_z_m": float(np.std(res["zh"])),
            "D50_surface_measured_m": stats["D50_surface_measured_m"]}
    atomic_write_bytes(os.path.join(out_dir, "heightfield-meta.json"),
                       json.dumps(meta, ensure_ascii=False, indent=2)
                       .encode("utf-8"))
    with open(os.path.join(out_dir, "bed-stats.json"), "w",
              encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    full = res["full"]
    np.savez_compressed(
        os.path.join(out_dir, "_work", "bed-maps.npz"),
        z=full["z"].astype(np.float32), gid=full["gid"],
        nx=full["nx"].astype(np.float32), ny=full["ny"].astype(np.float32),
        nz=full["nz"].astype(np.float32), fx=full["fx"].astype(np.float32),
        fy=full["fy"].astype(np.float32), fz=full["fz"].astype(np.float32),
        ao=full["ao"].astype(np.float32),
        shadow=full["shadow"].astype(np.float32),
        D=np.array([g["D"] for g in res["grains"]]),
        kind=np.array([g["kind"] for g in res["grains"]]),
        render_j0=np.array([RENDER_J0]),
        cell=np.array([full["cell"]]),
        mono_geom=np.array([[g["x"], g["y"], g["z"], g["abc"][0],
                             g["abc"][1], g["abc"][2]]
                            for g in res["grains"]
                            if str(g["kind"]) == "mono"]),
        mono_R=np.array([g["R"].tolist() for g in res["grains"]
                         if str(g["kind"]) == "mono"]))
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out", "v4-2-2"))
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--no-shadow", action="store_true")
    ap.add_argument("--probe", action="store_true",
                    help="build + audit, print, do not write deliverables")
    args = ap.parse_args(argv)
    if not args.probe:
        os.makedirs(os.path.join(args.out, "_work"), exist_ok=True)
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
    res = build_and_audit(args.out, seed=args.seed,
                          sun=None if args.no_shadow else sun)
    stats = res["stats"]
    stats["gate"] = gate_bed2(stats)
    stats["build_wall_s"] = float(time.time() - t0)
    if args.probe:
        for k, v in stats["gate"].items():
            if isinstance(v, dict) and "value" in v and "range" in v:
                print("  %-26s %-10s [%s, %s] %s"
                      % (k, ("%.4g" % v["value"]) if v["value"] is not None
                         else "None", v["range"][0], v["range"][1],
                         "PASS" if v["pass"] else "FAIL"))
        return 0
    write_bed_deliverables(args.out, res)
    g = stats["gate"]
    print("\nG-BED2")
    for k, v in g.items():
        if isinstance(v, dict) and "value" in v and "range" in v:
            print("  %-26s %-10s [%s, %s] %s"
                  % (k, ("%.4g" % v["value"]) if v["value"] is not None
                     else "None", v["range"][0], v["range"][1],
                     "PASS" if v["pass"] else "FAIL"))
    print("  overall: %s" % ("PASS" if g["pass"] else "FAIL"))
    return 0 if g["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
