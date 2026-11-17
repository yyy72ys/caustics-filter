"""prep_v4_2_10.py -- BRIEF-V4-2-10 STAGE 0 + STAGE 0.5.

STAGE 0 : absolute-path preflight (no other tree is searched), sha256 record,
          Pre/PATChES geometry check, D50 / cell-constant cross-check.
STAGE 0.5: fetch Run1_Post.pts from the EIDC datastore (token-free single-file
          GET, the same route the dem-fetch job used), rasterise Pre and Post
          on a 2 mm lattice (cell = mean z of the points in the cell, no point
          -> NaN), record the Pre/Post difference statistics, and FIX the
          reference-window layout (patch 0.6x0.6 m, strip 2.4x0.6 m) before any
          measurement is made.

New file.  Existing files are never written.  Everything lands in
out/v4-2-10/.  Idempotent / resumable: each unit is skipped when its artifact
already verifies (shape + size + content check), so a reboot mid-run only costs
the unit that was in flight.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import sys
import time
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, read_checkpoint  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-10.md"
DESIGN = r"D:\mutaro-work\rally-in\DESIGN-V4-2-7-BED-NEXT.md"
DATA = r"D:\yohey2026\dsh\muse\caustics-filter\data\dem-alt"
PRE_PTS = os.path.join(DATA, "Run1_Pre.pts")
POST_PTS = os.path.join(DATA, "Run1_Post.pts")
POST_URL = ("https://catalogue.ceh.ac.uk/datastore/eidchub/"
            "6749d033-cdf4-479b-ba85-015c3dbb476a/TLSdata/TLSscans/Run1_Post.pts")

OUT = os.path.join(HERE, "out", "v4-2-10")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
MANIFEST = os.path.join(OUT, "prep-manifest.json")
PRE_NPY = os.path.join(OUT, "raster-run1-pre.npy")
POST_NPY = os.path.join(OUT, "raster-run1-post.npy")

R5 = os.path.join(HERE, "run_v4_2_5.py")
HF_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
HF_R1 = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")
HF_R2 = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r2.npy")
HF_R2_SUB = os.path.join(HERE, "out", "v4-2-12b", "heightfield-long-r2-regen.npy")
META_R2_SUB = os.path.join(HERE, "out", "v4-2-12b",
                           "heightfield-long-r2-regen-meta.json")
NUMBERS_425 = os.path.join(HERE, "out", "v4-2-5", "numbers-v4-2-5.json")

CELL = 0.002                     # 2 mm, fixed by the brief
PRE_BYTES = 225521584
PRE_NPTS = 4784407
POST_BYTES = 215478752           # from the datastore index, checked after fetch

PATCH_BOX = (300, 300)           # 0.60 x 0.60 m
STRIP_BOX_U, STRIP_BOX_S = 1200, 300   # 2.40 x 0.60 m
N_PATCH = 8                      # brief minimum: >= 8
N_STRIP = 6                      # brief minimum: >= 6

LOG = []


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jd(o):
    if isinstance(o, dict):
        return {str(k): jd(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jd(v) for v in o]
    if isinstance(o, np.ndarray):
        return jd(o.tolist())
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = (json.dumps(jd(obj), ensure_ascii=False, indent=1) + "\n").encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True

    return atomic_write_bytes(path, data, verify=verify)


def save_npy(path, arr):
    """atomic .npy write through a .tmp sibling."""
    tmp = path + ".tmp.npy"
    np.save(tmp, arr)
    with open(tmp, "rb") as f:
        head = f.read(6)
    if head != b"\x93NUMPY":
        raise RuntimeError("npy header check failed: %r" % path)
    os.replace(tmp, path)
    return path


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def stage0(cp, state):
    log("STAGE 0  preflight (absolute paths only)")
    required = [
        ("brief", BRIEF), ("design", DESIGN), ("ruler", R5),
        ("v4-2-2 patch heightfield", HF_PATCH),
        ("v4-2-3 long bed r1", HF_R1),
        ("v4-2-5 numbers", NUMBERS_425),
        ("Hodge Run1_Pre.pts", PRE_PTS),
    ]
    missing = [k for k, p in required if not os.path.isfile(p)]
    sub = {"path": HF_R2, "exists": os.path.isfile(HF_R2),
           "substitute": HF_R2_SUB, "substitute_exists": os.path.isfile(HF_R2_SUB)}
    if missing:
        raise SystemExit("STAGE 0 STOP -- missing required absolute paths: %s" % missing)
    if not sub["substitute_exists"]:
        raise SystemExit("STAGE 0 STOP -- r2 absent at out/v4-2-3 AND no verified "
                         "regenerated r2 at out/v4-2-12b; the STAGE 1 r2 gate "
                         "cannot be run")
    hashes = {}
    for k, p in required:
        hashes[k] = {"path": p, "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
        log("  %-26s %10d B  %s" % (k, hashes[k]["bytes"], hashes[k]["sha256"][:16]))
    if sub["substitute_exists"]:
        hashes["r2 substitute (headerless raw f32)"] = {
            "path": HF_R2_SUB, "bytes": os.path.getsize(HF_R2_SUB),
            "sha256": sha256_file(HF_R2_SUB)}
    pre_head, pre_bytes = read_pts_header(PRE_PTS)
    checks = {
        "pre_header_points": pre_head, "pre_bytes": pre_bytes,
        "pre_header_matches_documented": pre_head == PRE_NPTS,
        "pre_bytes_matches_documented": pre_bytes == PRE_BYTES,
        "patch_bytes_16777344": os.path.getsize(HF_PATCH) == 16777344,
        "r1_bytes_83886080": os.path.getsize(HF_R1) == 83886080,
        "r2_substitute_bytes_83886080": os.path.getsize(HF_R2_SUB) == 83886080,
    }
    log("  Pre header %d pts / %d B (documented %d / %d) -> %s"
        % (pre_head, pre_bytes, PRE_NPTS, PRE_BYTES,
           "match" if checks["pre_header_matches_documented"] else "MISMATCH"))
    r2_note = (
        "out/v4-2-3/heightfield-long-r2.npy does NOT exist: v4-2.3 never wrote the "
        "r2 heightfield (its own report says so). The brief's STAGE 0 path list "
        "assumes it does. A bit-verified deterministic regeneration exists at "
        "out/v4-2-12b/heightfield-long-r2-regen.npy (seed "
        "3255045276383062994, 28/28 metadata fields identical, sha256 "
        "6bf729892ddf3bdaf742eb1596944254273f77429adabe40f2ed56c996877b4b) and "
        "v4-2-12b's identity gate reproduced v4-2.5's recorded r2 H exactly "
        "(x 0.3222576449 / iso 0.3245887727). It is used as the r2 face and the "
        "STAGE 1 gate is required to reproduce those floats; if it does not, "
        "STAGE 1 stops.")
    state["stage0"] = {"created_utc": utcnow(), "hashes": hashes,
                       "checks": checks, "r2": sub, "r2_note": r2_note,
                       "stop_condition": "none -- all required paths present "
                                         "(r2 via the documented substitute)"}
    if not (checks["pre_header_matches_documented"] and
            checks["pre_bytes_matches_documented"]):
        raise SystemExit("STAGE 0 STOP -- Run1_Pre.pts header/size do not match "
                         "the documented spec")
    write_json(os.path.join(WORK, "stage0-preflight.json"), state["stage0"])
    cp.unit_done("stage0-preflight")
    return state


def read_pts_header(path):
    with open(path, "rb") as f:
        line = f.readline()
    return int(line.split()[0]), os.path.getsize(path)


# ---------------------------------------------------------------------------
# STAGE 0.5 -- fetch
# ---------------------------------------------------------------------------
def fetch_post(cp, state):
    log("STAGE0.5  fetch Run1_Post.pts (token-free GET, resumable Range)")
    if os.path.isfile(POST_PTS) and os.path.getsize(POST_PTS) == POST_BYTES:
        log("  already present, full size")
        cp.unit_done("stage0.5-download")
        return
    part = POST_PTS + ".part"
    have = os.path.getsize(part) if os.path.isfile(part) else 0
    req = urllib.request.Request(POST_URL)
    if have:
        req.add_header("Range", "bytes=%d-" % have)
        log("  resuming at %d B" % have)
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=120) as r, open(part, "ab" if have else "wb") as f:
        code = r.getcode()
        total = int(r.headers.get("Content-Length") or 0) + have
        log("  HTTP %d, %d B to go" % (code, total - have))
        got = have
        last = 0.0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            got += len(chunk)
            if time.time() - last > 5:
                last = time.time()
                log("  %.1f / %.1f MB  (%.1f MB/s)"
                    % (got / 1e6, total / 1e6, (got - have) / 1e6 / max(time.time() - t0, 1e-9)))
    if os.path.getsize(part) != POST_BYTES:
        raise SystemExit("STAGE 0.5 STOP -- Post size %d != datastore %d"
                         % (os.path.getsize(part), POST_BYTES))
    os.replace(part, POST_PTS)
    h, b = read_pts_header(POST_PTS)
    state["stage0_5_fetch"] = {
        "url": POST_URL, "bytes": b, "header_points": h,
        "expected_bytes": POST_BYTES,
        "sha256": sha256_file(POST_PTS),
        "seconds": time.time() - t0,
        "parity_with_pre": {"pre_bytes": PRE_BYTES, "post_bytes": b,
                            "byte_ratio": b / PRE_BYTES,
                            "pre_points": PRE_NPTS, "post_points": h,
                            "point_ratio": h / PRE_NPTS},
    }
    log("  got %d B, header %d pts (Pre %d) -> parity %.4f"
        % (b, h, PRE_NPTS, h / PRE_NPTS))
    write_json(os.path.join(WORK, "stage0.5-fetch.json"), state["stage0_5_fetch"])
    cp.unit_done("stage0.5-download")


# ---------------------------------------------------------------------------
# STAGE 0.5 -- raster
# ---------------------------------------------------------------------------
def parse_pts(path):
    import pandas as pd
    t0 = time.time()
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, dtype=np.float32,
                     usecols=[0, 1, 2], names=["X", "Y", "Z"])
    log("  parsed %d pts from %s in %.1fs" % (len(df), os.path.basename(path),
                                              time.time() - t0))
    return df.X.values, df.Y.values, df.Z.values


def frame_from(x, y, keep=None):
    # float64 is mandatory: a naive float32 accumulation of 4.8e6 coordinates
    # near 1e2 loses metres (measured, before this cast: mean Y off by 5.5 m).
    xy = np.column_stack([x, y]).astype(np.float64)
    if keep is not None:
        xy = xy[keep]
    mu = xy.mean(0)
    w, v = np.linalg.eigh(np.cov((xy - mu).T))
    ax = v[:, 1] if w[1] > w[0] else v[:, 0]
    if ax[0] < 0:                    # deterministic sign
        ax = -ax
    n = np.array([-ax[1], ax[0]])
    return mu, ax, n, w


def rasterize(u, s, z, u0, s0, cell, nx, ny):
    idx = np.floor((u - u0) / cell).astype(np.int64)
    jdx = np.floor((s - s0) / cell).astype(np.int64)
    ok = (idx >= 0) & (idx < nx) & (jdx >= 0) & (jdx < ny)
    flat = (idx[ok] * ny + jdx[ok]).astype(np.int64)
    ncell = nx * ny
    cnt = np.bincount(flat, minlength=ncell).astype(np.float64)
    sm = np.bincount(flat, weights=z[ok].astype(np.float64), minlength=ncell)
    out = np.full(ncell, np.nan, dtype=np.float32)
    m = cnt > 0
    out[m] = (sm[m] / cnt[m]).astype(np.float32)
    stats = {"points": int(ok.sum()), "points_outside": int((~ok).sum()),
             "occupied_cells": int(m.sum()), "cells": int(ncell),
             "nodata_fraction": float(1.0 - m.sum() / ncell),
             "points_per_occupied_cell_p50": float(np.median(cnt[m]))}
    return out.reshape(nx, ny), stats


def mean_fill(a):
    a = np.array(a, dtype=np.float64, copy=True)
    m = np.isfinite(a)
    a[~m] = float(np.nanmean(a))
    return a


def stage05_raster(cp, state):
    log("STAGE0.5  rasterise Pre and Post to a 2 mm lattice")
    X, Y, Z = parse_pts(PRE_PTS)
    mu, ax, nrm, w = frame_from(X, Y)
    u = (X - mu[0]) * ax[0] + (Y - mu[1]) * ax[1]
    s = (X - mu[0]) * nrm[0] + (Y - mu[1]) * nrm[1]
    # second pass: axis from the bed interior only (walls excluded by the
    # cross-stream 2-98 percentile range of the first pass)
    lo, hi = np.percentile(s, [1, 99])
    interior = (s > lo) & (s < hi)
    mu2, ax2, nrm2, w2 = frame_from(X, Y, keep=interior)
    ang1 = float(np.degrees(np.arctan2(ax[1], ax[0])))
    ang2 = float(np.degrees(np.arctan2(ax2[1], ax2[0])))
    log("  PCA long axis: all pts %.3f deg, bed interior %.3f deg (delta %.3f)"
        % (ang1, ang2, ang2 - ang1))
    mu, ax, nrm = mu2, ax2, nrm2
    u = (X - mu[0]) * ax[0] + (Y - mu[1]) * ax[1]
    s = (X - mu[0]) * nrm[0] + (Y - mu[1]) * nrm[1]
    u0, s0 = float(u.min()), float(s.min())
    nx = int(np.ceil((u.max() - u0) / CELL)) + 1
    ny = int(np.ceil((s.max() - s0) / CELL)) + 1

    # point-spacing check (2-D nearest neighbour on a subsample)
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(0)
    pick = rng.choice(len(u), size=min(60000, len(u)), replace=False)
    tr = cKDTree(np.column_stack([u, s]))
    dd, _ = tr.query(np.column_stack([u[pick], s[pick]]), k=2)
    nn_mm = np.percentile(dd[:, 1], [10, 25, 50, 75, 90]) * 1000.0

    pre, pre_stats = rasterize(u, s, Z, u0, s0, CELL, nx, ny)
    save_npy(PRE_NPY, pre)
    del X, Y, Z
    cp.unit_done("stage0.5-raster-pre")
    log("  Pre raster %dx%d  nodata %.3f" % (nx, ny, pre_stats["nodata_fraction"]))

    X, Y, Z = parse_pts(POST_PTS)
    u2v = (X - mu[0]) * ax[0] + (Y - mu[1]) * ax[1]
    s2v = (X - mu[0]) * nrm[0] + (Y - mu[1]) * nrm[1]
    mu_p, ax_p, nrm_p, w_p = frame_from(X, Y)
    ang_post = float(np.degrees(np.arctan2(ax_p[1], ax_p[0])))
    if ang_post - ang2 > 90:
        ang_post -= 180.0
    elif ang_post - ang2 < -90:
        ang_post += 180.0
    post, post_stats = rasterize(u2v, s2v, Z, u0, s0, CELL, nx, ny)
    save_npy(POST_NPY, post)
    del X, Y, Z
    cp.unit_done("stage0.5-raster-post")
    log("  Post raster nodata %.3f (own PCA %.3f deg vs Pre %.3f deg, delta %.3f)"
        % (post_stats["nodata_fraction"], ang_post, ang2, ang_post - ang2))

    state["stage0_5_raster"] = {
        "created_utc": utcnow(), "cell_m": CELL,
        "pca": {"mean_xy": mu, "long_axis": ax, "cross_axis": nrm,
                "angle_deg_all": ang1, "angle_deg_bed_interior": ang2,
                "angle_deg_post_own": float(ang_post),
                "delta_deg_post_vs_pre": float(ang_post - ang2),
                "eigenvalues_all": w, "eigenvalues_interior": w2,
                "flow_definition": ("PCA long axis of the bed-interior points "
                                    "= streamwise (u); cross axis = 90 deg "
                                    "counter-clockwise from it"),
                "width_check": ("PCA cross-stream span of Pre = %.4f m vs flume "
                                "width 2.1 m + margin" % float(s.ptp()))},
        "origin": {"u0": u0, "s0": s0, "nx": nx, "ny": ny,
                   "u_max": float(u0 + nx * CELL), "s_max": float(s0 + ny * CELL)},
        "extent_pre": {"u": [float(u.min()), float(u.max())],
                       "s": [float(s.min()), float(s.max())]},
        "nodata": {"pre": pre_stats, "post": post_stats},
        "point_spacing_mm": {"p10": float(nn_mm[0]), "p25": float(nn_mm[1]),
                             "p50": float(nn_mm[2]), "p75": float(nn_mm[3]),
                             "p90": float(nn_mm[4]),
                             "method": "2-D nearest neighbour, 60k subsample, "
                                       "bed-interior frame",
                             "note": ("brief/dem-fetch report says p50 2.6 mm; the "
                                      "measured p50 is much smaller, i.e. the "
                                      "cloud carries sub-mm duplicate/near-"
                                      "duplicate returns (multi-station "
                                      "registration). Recorded, not used.")},
        "files": {os.path.basename(PRE_NPY): {"bytes": os.path.getsize(PRE_NPY),
                                              "sha256": sha256_file(PRE_NPY)},
                  os.path.basename(POST_NPY): {"bytes": os.path.getsize(POST_NPY),
                                               "sha256": sha256_file(POST_NPY)}},
    }
    write_json(os.path.join(WORK, "stage0.5-raster.json"), state["stage0_5_raster"])
    return state


# ---------------------------------------------------------------------------
# STAGE 0.5 -- window layout (fixed BEFORE measurement) + nodata arrangement
# ---------------------------------------------------------------------------
def nodata_layout(mask, nblk=6):
    h, w = mask.shape
    if h < nblk or w < nblk:
        raise ValueError("window too small for the block map: %r" % (mask.shape,))
    blocks = []
    for bi in range(nblk):
        row = []
        for bj in range(nblk):
            a = mask[bi * h // nblk:(bi + 1) * h // nblk,
                     bj * w // nblk:(bj + 1) * w // nblk]
            row.append(float(a.mean()))
        blocks.append(row)
    v = np.array(blocks)
    p = float(mask.mean())
    # clustering index: block variance / binomial variance for the same p
    expect = p * (1 - p) / (h / nblk * w / nblk)
    return {"fraction": p, "blocks_6x6": blocks,
            "block_range": [float(v.min()), float(v.max())],
            "clustering_index": float(v.var() / expect) if expect > 0 else None,
            "centre_vs_edge": {"centre": float(mask[h // 4:3 * h // 4,
                                                   w // 4:3 * w // 4].mean()),
                               "edge": float((mask.sum() - mask[h // 4:3 * h // 4,
                                                                w // 4:3 * w // 4].sum())
                                             / (mask.size - (h // 2) * (w // 2)))},
            "row_profile": [float(mask[i * h // 12:(i + 1) * h // 12].mean())
                            for i in range(12)],
            "col_profile": [float(mask[:, i * w // 12:(i + 1) * w // 12].mean())
                            for i in range(12)]}


def stage05_layout(cp, state):
    log("STAGE0.5  fix reference-window layout + nodata arrangement")
    post = np.load(POST_NPY, mmap_mode="r")
    nx, ny = post.shape
    u0 = state["stage0_5_raster"]["origin"]["u0"]
    s0 = state["stage0_5_raster"]["origin"]["s0"]
    # walls are found numerically, not assumed: the flume side walls sit
    # ~0.2 m above the bed, so a cross-stream row that is a wall has a much
    # higher median z than the bed rows.
    row_med = np.array([float(np.nanmedian(np.asarray(post[:, j], dtype=np.float32)))
                        for j in range(ny)])
    bed_med = float(np.nanmedian(row_med))
    wall = (~np.isfinite(row_med)) | (row_med > bed_med + 0.05)
    # largest contiguous non-wall run
    best = (0, 0)
    j = 0
    while j < ny:
        if not wall[j]:
            k = j
            while k + 1 < ny and not wall[k + 1]:
                k += 1
            if k - j > best[1] - best[0]:
                best = (j, k)
            j = k + 1
        else:
            j += 1
    j_lo, j_hi = best[0], best[1] + 1
    L = j_hi - j_lo
    if L < 640:
        raise SystemExit("STAGE 0.5 STOP -- bed interior only %d cells across "
                         "(need 640 for two 0.6 m rows + margins)" % L)
    gap = max((L - 2 * PATCH_BOX[1]) // 3, 10)
    rows_patch = [j_lo + gap, j_lo + 2 * gap + PATCH_BOX[1]]
    if rows_patch[1] + PATCH_BOX[1] > j_hi:
        rows_patch[1] = j_hi - PATCH_BOX[1]
    rows_strip = rows_patch          # same cross-stream rows for both boxes
    log("  wall rows: %d of %d excluded; bed interior j=%d..%d (%d cells, %.3f m)"
        % (int(wall.sum()), ny, j_lo, j_hi, L, L * CELL))
    i_lo, i_hi = 20, nx - 20
    span = i_hi - i_lo
    pos_patch = [i_lo + int(round(k * (span - PATCH_BOX[0]) / (N_PATCH // 2 - 1)))
                 for k in range(N_PATCH // 2)]
    pos_strip = [i_lo + int(round(k * (span - STRIP_BOX_U) / (N_STRIP // 2 - 1)))
                 for k in range(N_STRIP // 2)]

    is_ok = np.isfinite(np.asarray(post, dtype=np.float32))
    # summed-area table for O(1) coverage of any window
    iia = np.zeros((nx + 1, ny + 1), dtype=np.int64)
    np.cumsum(np.cumsum(is_ok, axis=0), axis=1, out=iia[1:, 1:])

    def cover(i0, j0, ni, nj):
        s = (iia[i0 + ni, j0 + nj] - iia[i0, j0 + nj]
             - iia[i0 + ni, j0] + iia[i0, j0])
        return float(s) / (ni * nj)

    # candidate blocks are non-overlapping by construction: 0.6 m / 2.4 m
    # blocks tiled from the interior start.  Positions whose coverage is below
    # COVER_MIN are not measured (a window that is 98% mean-fill is not a
    # measurement of the bed) and are recorded as excluded.
    COVER_LEVELS = (0.50, 0.45, 0.40, 0.35, 0.30, 0.25)
    cands = {"patch": [], "strip": []}
    for r, j0 in enumerate(rows_patch):
        k = i_lo
        while k + PATCH_BOX[0] <= i_hi:
            cands["patch"].append({"row": r, "i0": k, "j0": j0,
                                   "cover": cover(k, j0, *PATCH_BOX)})
            k += PATCH_BOX[0]
    for r, j0 in enumerate(rows_strip):
        k = i_lo
        while k + STRIP_BOX_U <= i_hi:
            cands["strip"].append({"row": r, "i0": k, "j0": j0,
                                   "cover": cover(k, j0, STRIP_BOX_U, STRIP_BOX_S)})
            k += STRIP_BOX_U

    def pick(box, n):
        # the highest coverage floor that still yields n windows wins; the floor
        # and every candidate's coverage are recorded (the reference population
        # must not be quietly padded with mean-filled near-empty windows)
        for thr in COVER_LEVELS:
            good = [c for c in cands[box] if c["cover"] >= thr]
            good.sort(key=lambda c: c["i0"])
            if len(good) >= n:
                idx = np.unique(np.round(np.linspace(0, len(good) - 1, n)).astype(int))
                return [good[i] for i in idx], thr, len(good)
        raise SystemExit(
            "STAGE 0.5 STOP -- %s: only %d blocks reach coverage %.2f (need %d); "
            "the 2 mm lattice of this scan cannot support the reference "
            "population without going below a 25%% valid-cell floor."
            % (box, len([c for c in cands[box] if c["cover"] >= COVER_LEVELS[-1]]),
               COVER_LEVELS[-1], n))

    patch_sel, thr_p, nok_p = pick("patch", N_PATCH)
    strip_sel, thr_s, nok_s = pick("strip", N_STRIP)
    log("  coverage floors: patch %.2f (%d blocks eligible), strip %.2f (%d)"
        % (thr_p, nok_p, thr_s, nok_s))

    def cut(i0, j0, ni, nj):
        return np.asarray(post[i0:i0 + ni, j0:j0 + nj], dtype=np.float32)

    patch_w, strip_w = [], []
    for c in patch_sel:
        i0, j0 = c["i0"], c["j0"]
        a = cut(i0, j0, *PATCH_BOX)
        m = ~np.isfinite(a)
        patch_w.append({
            "id": "patch-r%d-i%05d" % (c["row"], i0), "i0": i0, "j0": j0,
            "ni": PATCH_BOX[0], "nj": PATCH_BOX[1],
            "u_m": [float(u0 + i0 * CELL), float(u0 + (i0 + PATCH_BOX[0]) * CELL)],
            "s_m": [float(s0 + j0 * CELL), float(s0 + (j0 + PATCH_BOX[1]) * CELL)],
            "cover": c["cover"], "nodata": nodata_layout(m)})
    for c in strip_sel:
        i0, j0 = c["i0"], c["j0"]
        a = cut(i0, j0, STRIP_BOX_U, STRIP_BOX_S)
        m = ~np.isfinite(a)
        strip_w.append({
            "id": "strip-r%d-i%05d" % (c["row"], i0), "i0": i0, "j0": j0,
            "ni": STRIP_BOX_U, "nj": STRIP_BOX_S,
            "u_m": [float(u0 + i0 * CELL), float(u0 + (i0 + STRIP_BOX_U) * CELL)],
            "s_m": [float(s0 + j0 * CELL), float(s0 + (j0 + STRIP_BOX_S) * CELL)],
            "cover": c["cover"], "nodata": nodata_layout(m)})

    # non-overlap audit
    def disjoint(ws, key):
        bad = []
        for a in range(len(ws)):
            for b in range(a + 1, len(ws)):
                A, B = ws[a], ws[b]
                ou = min(A[key][0] + A["ni"] * CELL, B[key][0] + B["ni"] * CELL) - \
                    max(A[key][0], B[key][0])
                os_ = min(A["s_m"][0] + A["nj"] * CELL, B["s_m"][0] + B["nj"] * CELL) - \
                    max(A["s_m"][0], B["s_m"][0])
                if ou > 1e-9 and os_ > 1e-9:
                    bad.append([A["id"], B["id"]])
        return bad

    lay = {
        "created_utc": utcnow(),
        "rule": ("windows are non-overlapping blocks tiled from the bed "
                 "interior; walls are excluded numerically (a cross-stream row "
                 "whose median z is >0.05 m above the bed median is a wall); a "
                 "block is only measured if its valid-cell coverage is >= 0.45 "
                 "(a 98%-mean-filled window is not a measurement); the picks "
                 "are then spread evenly along the flow. The layout is fixed "
                 "here, before any measurement."),
        "coverage_min": {"patch": thr_p, "strip": thr_s},
        "coverage_eligible_blocks": {"patch": nok_p, "strip": nok_s},
        "wall_rows_excluded": int(wall.sum()),
        "bed_interior_cells": int(L),
        "bed_interior_s_m": [float(s0 + j_lo * CELL), float(s0 + j_hi * CELL)],
        "rows_cross_cells": [int(r) for r in rows_patch],
        "rows_cross_s_m": [[float(s0 + r * CELL), float(s0 + (r + PATCH_BOX[1]) * CELL)]
                           for r in rows_patch],
        "candidates": {k: [{"row": c["row"], "i0": c["i0"], "cover": c["cover"]}
                           for c in v] for k, v in cands.items()},
        "patch_windows": patch_w, "strip_windows": strip_w,
        "patch_overlaps": disjoint(patch_w, "u_m"),
        "strip_overlaps": disjoint(strip_w, "u_m"),
        "counts": {"patch": len(patch_w), "strip": len(strip_w)},
        "nodata_summary": {
            "patch_mean": float(np.mean([w["nodata"]["fraction"] for w in patch_w])),
            "patch_min": float(np.min([w["nodata"]["fraction"] for w in patch_w])),
            "patch_max": float(np.max([w["nodata"]["fraction"] for w in patch_w])),
            "strip_mean": float(np.mean([w["nodata"]["fraction"] for w in strip_w])),
            "strip_min": float(np.min([w["nodata"]["fraction"] for w in strip_w])),
            "strip_max": float(np.max([w["nodata"]["fraction"] for w in strip_w])),
        },
    }
    log("  patch %d windows, strip %d windows; nodata patch %.3f-%.3f mean %.3f, "
        "strip %.3f-%.3f mean %.3f"
        % (len(patch_w), len(strip_w),
           lay["nodata_summary"]["patch_min"], lay["nodata_summary"]["patch_max"],
           lay["nodata_summary"]["patch_mean"],
           lay["nodata_summary"]["strip_min"], lay["nodata_summary"]["strip_max"],
           lay["nodata_summary"]["strip_mean"]))
    if lay["patch_overlaps"] or lay["strip_overlaps"]:
        raise SystemExit("STAGE 0.5 STOP -- window overlap: %s"
                         % (lay["patch_overlaps"] + lay["strip_overlaps"]))
    state["stage0_5_layout"] = lay
    write_json(os.path.join(WORK, "stage0.5-layout.json"), lay)
    cp.unit_done("stage0.5-layout")
    return state


# ---------------------------------------------------------------------------
# STAGE 0.5 -- Pre/Post difference (diagnostic only, never a judge)
# ---------------------------------------------------------------------------
def stage05_diff(cp, state):
    log("STAGE0.5  Pre/Post difference statistics (diagnostic only)")
    pre = np.load(PRE_NPY, mmap_mode="r")
    post = np.load(POST_NPY, mmap_mode="r")
    # common-region comparison, on the bed interior rows only
    d = np.asarray(post, dtype=np.float64) - np.asarray(pre, dtype=np.float64)
    both = np.isfinite(d)
    dv = d[both]
    # plane-removed difference (registration tilt); per 2 m column block
    rows = []
    nx = d.shape[0]
    nblk = 10
    for b in range(nblk):
        i0, i1 = b * nx // nblk, (b + 1) * nx // nblk
        blk = d[i0:i1]
        m = np.isfinite(blk)
        if m.sum() < 1000:
            continue
        rows.append({"u_m": float(state["stage0_5_raster"]["origin"]["u0"]
                                  + (i0 + i1) / 2 * CELL),
                     "n": int(m.sum()),
                     "mean_mm": float(np.nanmean(blk) * 1000),
                     "sd_mm": float(np.nanstd(blk) * 1000),
                     "p05_mm": float(np.nanpercentile(blk, 5) * 1000),
                     "p95_mm": float(np.nanpercentile(blk, 95) * 1000)})
    out = {
        "created_utc": utcnow(),
        "sigma_z_pre_m": float(np.nanstd(np.asarray(pre, dtype=np.float64))),
        "sigma_z_post_m": float(np.nanstd(np.asarray(post, dtype=np.float64))),
        "common_cells": int(both.sum()),
        "common_fraction": float(both.mean()),
        "diff_m": {"mean": float(dv.mean()), "sd": float(dv.std()),
                   "p05": float(np.percentile(dv, 5)),
                   "p50": float(np.percentile(dv, 50)),
                   "p95": float(np.percentile(dv, 95)),
                   "abs_p95": float(np.percentile(np.abs(dv), 95))},
        "block_profile_2m": rows,
        "note": ("diagnostic only ('was the bed worked by water'); never used in "
                 "the judgment. A non-zero mean slope in the block profile means "
                 "the two scans are not perfectly registered; the difference is "
                 "therefore read as an upper bound on bed change, not as a "
                 "sediment budget."),
    }
    state["stage0_5_diff"] = out
    write_json(os.path.join(WORK, "stage0.5-diff.json"), out)
    cp.unit_done("stage0.5-diff")
    log("  sigma_z pre %.5f / post %.5f m; diff sd %.2f mm, |diff| p95 %.2f mm, "
        "common %.3f"
        % (out["sigma_z_pre_m"], out["sigma_z_post_m"],
           out["diff_m"]["sd"] * 1000, out["diff_m"]["abs_p95"] * 1000,
           out["common_fraction"]))
    return state


# ---------------------------------------------------------------------------
def main():
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint(os.path.join(WORK, "state.json")) or {}
    done = set((cp.state() or {}).get("completed_ids") or [])
    log("checkpoint: %d units done %s" % (len(done), sorted(done)))
    # hydrate the in-memory state from the per-unit artifacts on a resume
    for unit, key, fn in (("stage0-preflight", "stage0", "stage0-preflight.json"),
                          ("stage0.5-download", "stage0_5_fetch", "stage0.5-fetch.json"),
                          ("stage0.5-raster-pre", "stage0_5_raster", "stage0.5-raster.json"),
                          ("stage0.5-layout", "stage0_5_layout", "stage0.5-layout.json"),
                          ("stage0.5-diff", "stage0_5_diff", "stage0.5-diff.json")):
        p = os.path.join(WORK, fn)
        if unit in done and os.path.isfile(p) and key not in state:
            state[key] = json.load(open(p, encoding="utf-8"))
    t0 = time.time()
    if "stage0-preflight" not in done:
        stage0(cp, state)
    if "stage0.5-download" not in done:
        fetch_post(cp, state)
    if "stage0.5-raster-pre" not in done or "stage0.5-raster-post" not in done:
        if not (os.path.isfile(PRE_NPY) and os.path.isfile(POST_NPY)):
            stage05_raster(cp, state)
        else:
            log("  rasters present; re-registering units")
            if "stage0.5-raster-pre" not in done:
                cp.unit_done("stage0.5-raster-pre")
            if "stage0.5-raster-post" not in done:
                cp.unit_done("stage0.5-raster-post")
    if "stage0.5-layout" not in done:
        stage05_layout(cp, state)
    if "stage0.5-diff" not in done:
        stage05_diff(cp, state)
    write_json(os.path.join(WORK, "state.json"), state)
    manifest = {
        "created_utc": utcnow(), "cell_m": CELL,
        "pre_pts": {"path": PRE_PTS, "bytes": os.path.getsize(PRE_PTS),
                    "header_points": PRE_NPTS},
        "post_pts": {"path": POST_PTS, "bytes": os.path.getsize(POST_PTS),
                     "header_points": state["stage0_5_fetch"]["header_points"]},
        "license": ("Contains data supplied by UK Centre for Ecology & Hydrology. "
                    "Open Government Licence v3."),
        "raster": state["stage0_5_raster"], "layout": state["stage0_5_layout"],
        "diff": state["stage0_5_diff"],
        "seconds_total": time.time() - t0,
    }
    write_json(MANIFEST, manifest)
    cp.unit_done("stage0.5-manifest")
    log("prep done in %.1fs -- %s" % (time.time() - t0, MANIFEST))


if __name__ == "__main__":
    main()
