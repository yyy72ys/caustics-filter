"""run_v4_2_12 -- BRIEF-V4-2-12: re-measure the existing generated bed with
Nikora's procedure, plus an inexpensive synthetic-fill (filleting) probe.

MEASUREMENT ONLY.  No generation, no new bed, no rendering, no parameter
change.  The only files this job writes are this .py, everything under
out/v4-2-12/, and its own jobkit checkpoint inside that directory.

What the brief asks
-------------------
STAGE 0  absolute-path preflight (stop if anything is missing -- do NOT go
         looking in other trees), geometry / D50 recorded from the generator
         side, and a reproduction gate: the v4-2.5 published numbers must come
         back within +-0.01 using run_v4_2_5's own functions (verbatim import).
STAGE 1  Nikora's procedure on the same three faces, two directions each:
         raw field (no detrending), axis cuts D(lx,0) / D(0,ly) on a 60-bin
         log ladder, break chosen by brute force over ladder points inside
         [0.25 D50, 2 D50] (two-segment log-log least squares, minimum total
         residual), then H_Nikora = single fit [min lag, break], the two-band
         slopes [min lag, break] and [break, 4 D50], the break sensitivity
         (+-0.8 / 1.25), and the whole H_loc(l) ladder.
STAGE 2  patch only: priority-flood depression fill F, then
         z'(d) = min(F, z + d) for d in {0.25, 0.5, 0.75} sigma_z, re-measured
         with both instruments (the fixed 0.4-3.6 cm band H_fit, and the
         Nikora procedure with the break re-identified on the filled surface).
         Prediction, recorded before measuring: at d = 0.5 sigma_z,
         dH_fit <= +0.1.

Conventions fixed BEFORE the first measurement (not tuned afterwards)
---------------------------------------------------------------------
* cell         : 0.6 m / 2048 px = 0.29296875 mm, exactly v4-2.2 / v4-2.3.
* ladder       : 60 log-spaced lags, rounded to whole pixels, np.unique'd
                 (identical arithmetic to run_v4_2_5.curve_dir, imported).
* dmax (axis cut):
    patch  x and y : 0.30 m (half of the 0.6 m box, the v4-2.5 choice)
    long   flow  x : 1.00 m (1/3 of the 3.0 m bed, the v4-2.3 choice)
    long   cross y : 0.30 m (half of the 0.6 m width)
* flow axis     : arrays are (ny, nx); flow / long axis = axis 1 = x,
                 cross-stream = axis 0 = y (same as run_v4_2_3/2_5).
* periodicity   : our fields are periodic by construction, so np.roll is the
                 correct wrap for the lag differences; the report says so.
* NO detrending : the brief forbids it for STAGE 1 (Nikora measures the bed
                 surface, including its form).
* break search  : candidates = ladder lags inside [0.25 D50, 2 D50]; each
                 candidate splits the curve into [first lag .. brk] and
                 [brk .. last lag]; both segments need >= 4 points; the split
                 with the smallest total log-log residual wins.
* fill          : Barnes priority-flood (min-heap from the open border,
                 level = max(z, level)); open (non-periodic) boundary is
                 declared, since priority-flood is an open-domain algorithm.

Stages are checkpointed one jobkit unit at a time (stage0 / stage1 / stage2 /
report).  The run is foreground: it either finishes all deliverables or stops
at a gate and writes the stop record.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import heapq
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-12.md"

OUT = os.path.join(HERE, "out", "v4-2-12")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-12.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-12.md")
PREFLIGHT = os.path.join(OUT, "stage0-preflight.json")

# --- STAGE 0: the brief's absolute paths, verbatim -------------------------
# One missing file is a STOP.  Nothing here is searched for elsewhere.
REQUIRED = [
    ("ruler (source of the verbatim functions)", os.path.join(HERE, "run_v4_2_5.py")),
    ("patch bed", os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")),
    ("long bed r1", os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")),
    ("long bed r2", os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r2.npy")),
]
GEN_SCRIPTS = [
    ("run_v4_2_2.py", os.path.join(HERE, "run_v4_2_2.py")),
    ("run_v4_2_3.py", os.path.join(HERE, "run_v4_2_3.py")),
    ("bed_gen2.py", os.path.join(HERE, "bed_gen2.py")),
]

HF_PATCH = REQUIRED[1][1]
HF_R1 = REQUIRED[2][1]
HF_R2 = REQUIRED[3][1]
META_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
META_R1 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta.json")
META_R2 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta-2.json")

THIS_FILE = "run_v4_2_12.py"

# --- the measurement contract (fixed before the first measurement) ---------
CELL = 0.6 / 2048.0                   # 0.00029296875 m, exactly v4-2.2/2.3
LONG_NY, LONG_NX = 2048, 10240        # 0.6 m x 3.0 m
NBIN = 60
DM_PATCH = 0.30
DM_LONG_X = 1.00
DM_LONG_Y = 0.30
FIXED_BAND = (0.004, 0.036)           # the ruler: 0.4-3.6 cm, unchanged
BREAK_BAND = (0.25, 2.0)              # break search window, x D50
TWO_BAND_TOP = 4.0                    # second band: break .. 4 x D50
MIN_SEG_POINTS = 4

# reproduction gate: v4-2.5's published scaling-band values (brief, STAGE 0)
GATE_EXPECT = {
    "patch_x": 0.323, "patch_y": 0.373, "patch_iso": 0.349,
    "long_r1_x": 0.294, "long_r1_y": 0.297,
    "long_r2_x": 0.322, "long_r2_y": 0.325,
}
GATE_TOL = 0.01

# STAGE 2 prediction, recorded before measuring
PREDICTION = {"delta_over_sigma": 0.5, "dH_fit_max": 0.10}
DELTAS = (0.25, 0.50, 0.75)

LOG = []


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    LOG.append(line)
    print(line, flush=True)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = (json.dumps(obj, ensure_ascii=False, indent=1) + "\n").encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            return json.loads(f.read().decode("utf-8")) is not None

    atomic_write_bytes(path, data, verify=verify)
    return path


class Blocker(Exception):
    """A gate failed.  Carries the partial context so a stop report can be
    written honestly from what was actually measured."""

    def __init__(self, reason, ctx):
        super().__init__(reason)
        self.reason = reason
        self.ctx = ctx


# ---------------------------------------------------------------------------
# estimators -- imported verbatim from run_v4_2_5 (the ruler's正本)
# ---------------------------------------------------------------------------
def ruler():
    """Import run_v4_2_5 and hand back its own estimators, unmodified."""
    import run_v4_2_5 as R5
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band"):
        assert hasattr(R5, name), "run_v4_2_5 lost %s" % name
    return R5


def axis_cut(z, cell, axis, dmax, R5, nbin=NBIN):
    """D(lx, 0) or D(0, ly): one axis, 60-bin log ladder, periodic roll.

    Same arithmetic as run_v4_2_5.curve_dir (imported and called, not copied).
    """
    return R5.curve_dir(z, cell, dmax, nbin=nbin, axis=axis)


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight, geometry, D50, reproduction gate
# ---------------------------------------------------------------------------
def probe_inputs():
    """Format / shape / lattice / spread of whatever is on disk, measured.

    Runs before the missing-file gate on purpose: the brief asks STAGE 0 to
    record the format, shape and lattice spacing of the beds, and that part is
    answerable for the files that do exist even when another one is missing.
    """
    out = {}
    expected_long = LONG_NY * LONG_NX * 4
    if os.path.isfile(HF_PATCH):
        a = np.load(HF_PATCH, mmap_mode="r")
        meta = {}
        if os.path.isfile(META_PATCH):
            meta = json.load(open(META_PATCH, encoding="utf-8"))
        z = np.asarray(a, dtype=np.float64)
        out["patch"] = {
            "present": True, "container": "npy (with header)",
            "dtype": str(a.dtype), "shape": list(a.shape),
            "cell_m_from_generator": CELL,
            "cell_m_meta": meta.get("cell_m"),
            "size_m": [a.shape[1] * CELL, a.shape[0] * CELL],
            "sigma_z_m": float(np.std(z)),
            "skew": float(((z - z.mean()) ** 3).mean() / max(np.std(z), 1e-30) ** 3),
            "D50_surface_m_meta": meta.get("D50_surface_measured_m"),
        }
        del z
    for tag, path in (("long_r1", HF_R1), ("long_r2", HF_R2)):
        if os.path.isfile(path):
            n = os.path.getsize(path)
            head = np.fromfile(path, dtype="<f4", count=LONG_NY * 64).astype(np.float64)
            out[tag] = {
                "present": True, "container": "headerless raw float32 (<f4)",
                "bytes": n, "bytes_expected": expected_long,
                "bytes_ok": bool(n == expected_long),
                "shape_expected": [LONG_NY, LONG_NX],
                "sample_stats": {"n_cells": int(head.size),
                                 "min_m": float(head.min()),
                                 "max_m": float(head.max()),
                                 "std_m": float(head.std())},
            }
        else:
            out[tag] = {"present": False, "container": "headerless raw float32 (<f4)",
                        "bytes": None, "bytes_expected": expected_long,
                        "shape_expected": [LONG_NY, LONG_NX]}
    for tag, meta_path in (("long_r1_meta", META_R1), ("long_r2_meta", META_R2)):
        if os.path.isfile(meta_path):
            m = json.load(open(meta_path, encoding="utf-8"))
            out[tag] = {"present": True, "path": meta_path,
                        "shape": m.get("shape"), "cell_m": m.get("cell_m"),
                        "length_m": m.get("length_m"), "width_m": m.get("width_m"),
                        "surface_D50_m": m.get("surface_D50_m"),
                        "sigma_z_m": m.get("sigma_z_m"),
                        "role": m.get("role", "primary (r1)")}
    return out


def stage0(cp):
    log("STAGE0  absolute-path preflight (brief list only; no other trees)")
    rows, missing = [], []
    for label, path in REQUIRED:
        ok = os.path.isfile(path)
        row = {"label": label, "path": path, "exists": bool(ok),
               "bytes": os.path.getsize(path) if ok else None,
               "sha256": sha256_file(path) if ok else None}
        rows.append(row)
        log("  %-46s %s  %s" % (label, "OK" if ok else "MISSING",
                                row["bytes"] if ok else path))
        if not ok:
            missing.append(row)

    gen_rows, gen_present = [], []
    for label, path in GEN_SCRIPTS:
        ok = os.path.isfile(path)
        gen_rows.append({"label": label, "path": path, "exists": bool(ok)})
        if ok:
            gen_present.append(label)
    log("  generator scripts present: %s" % (gen_present or "NONE"))

    ctx = {"checked_utc": utcnow(), "required": rows, "generator_scripts": gen_rows,
           "missing": missing, "probe": probe_inputs(),
           "geometry": None, "d50": None, "repro": None}

    if missing:
        names = ", ".join(os.path.basename(r["path"]) for r in missing)
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: missing input(s) %s -- STOP per brief"
                     % names)
        raise Blocker("STAGE 0 preflight failed: %s missing (%s)"
                      % (names, "; ".join(r["path"] for r in missing)), ctx)
    if not gen_present:
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: no generator script to date the geometry")
        raise Blocker("STAGE 0: none of run_v4_2_2.py / run_v4_2_3.py / "
                      "bed_gen2.py is present, so the physical size cannot be "
                      "confirmed", ctx)

    # --- geometry, from the generator side + the saved metadata -------------
    log("STAGE0  geometry and D50 from the generator side")
    import bed_gen2 as B2
    patch_m = float(B2.PATCH)
    hires = int(B2.HIRES)
    cell_gen = patch_m / hires
    zp = np.load(HF_PATCH, mmap_mode="r")
    meta_p = json.load(open(META_PATCH, encoding="utf-8"))
    meta_r1 = json.load(open(META_R1, encoding="utf-8"))
    meta_r2 = json.load(open(META_R2, encoding="utf-8")) if os.path.isfile(META_R2) else {}

    geo = {
        "cell_m_from_generator": cell_gen,
        "cell_m_patch_meta": meta_p["cell_m"],
        "cell_check": bool(abs(cell_gen - CELL) < 1e-15
                           and abs(meta_p["cell_m"] - CELL) < 1e-15),
        "patch": {
            "shape": list(zp.shape), "cell_m": CELL,
            "size_m": [float(zp.shape[1] * CELL), float(zp.shape[0] * CELL)],
            "patch_m_from_bed_gen2": patch_m,
            "hires_from_bed_gen2": hires,
            "size_source": "bed_gen2.PATCH / BED_GEN2.HIRES + heightfield-meta.json",
        },
        "long_common": {
            "shape_expected": [LONG_NY, LONG_NX],
            "cell_m": CELL,
            "size_m": [LONG_NX * CELL, LONG_NY * CELL],
            "length_m_meta": meta_r1.get("length_m"),
            "width_m_meta": meta_r1.get("width_m"),
            "size_source": "run_v4_2_3.py (LONG_LX=3.0, LONG_LY=0.6 at the "
                           "v4-2.2 cell) + longbed-meta*.json",
        },
        "bytes_check": {
            os.path.basename(HF_R1): os.path.getsize(HF_R1) == LONG_NY * LONG_NX * 4,
            os.path.basename(HF_R2): os.path.getsize(HF_R2) == LONG_NY * LONG_NX * 4,
        },
    }
    d50 = {
        "patch_m": float(meta_p["D50_surface_measured_m"]),
        "long_r1_m": float(meta_r1["surface_D50_m"]),
        "long_r2_m": float(meta_r2["surface_D50_m"]),
    }
    unique = (geo["cell_check"]
              and geo["patch"]["size_m"] == [patch_m, patch_m]
              and all(geo["bytes_check"].values())
              and abs(geo["long_common"]["length_m_meta"] - 3.0) < 1e-12
              and abs(geo["long_common"]["width_m_meta"] - 0.6) < 1e-12)
    geo["unique"] = bool(unique)
    log("  cell %.10f m (unique=%s); D50 patch %.6f m, r1 %.6f m, r2 %.6f m"
        % (CELL, unique, d50["patch_m"], d50["long_r1_m"], d50["long_r2_m"]))
    if not unique:
        ctx.update({"geometry": geo, "d50": d50})
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: physical size / grid not uniquely determined")
        raise Blocker("STAGE 0: the physical size of the beds is not uniquely "
                      "determined (geometry check failed)", ctx)

    # --- reproduction gate: v4-2.5's published numbers, +-0.01 --------------
    log("STAGE0  reproduction gate vs v4-2.5 (band %.3f-%.3f m, tol %.2f)"
        % (FIXED_BAND[0], FIXED_BAND[1], GATE_TOL))
    R5 = ruler()
    got = {}

    def band_H(z, cell, axis, dmax, lo, hi):
        d, s = axis_cut(z, cell, axis, dmax, R5)
        return d, s, R5.fit_band(d, s, lo, hi)

    zp64 = np.asarray(zp, dtype=np.float64)
    dx, sx, fx = band_H(zp64, CELL, 1, DM_PATCH, *FIXED_BAND)
    dy, sy, fy = band_H(zp64, CELL, 0, DM_PATCH, *FIXED_BAND)
    di, si = R5.curve_iso(zp64, CELL, 0.20, nbin=NBIN)
    fi = R5.fit_band(di, si, *FIXED_BAND)
    got.update({"patch_x": fx["H"], "patch_y": fy["H"], "patch_iso": fi["H"]})
    log("  patch    x %.4f  y %.4f  iso %.4f" % (fx["H"], fy["H"], fi["H"]))

    for tag, path in (("long_r1", HF_R1), ("long_r2", HF_R2)):
        z = np.fromfile(path, dtype="<f4").reshape(LONG_NY, LONG_NX).astype(np.float64)
        dz, sz, fz = band_H(z, CELL, 1, DM_LONG_X, *FIXED_BAND)
        dw, sw, fw = band_H(z, CELL, 0, DM_LONG_Y, *FIXED_BAND)
        got.update({tag + "_x": fz["H"], tag + "_y": fw["H"]})
        log("  %-7s x %.4f  y %.4f" % (tag, fz["H"], fw["H"]))
        del z

    dev = {k: {"expected": GATE_EXPECT[k], "got": got[k],
               "abs_diff": abs(got[k] - GATE_EXPECT[k])} for k in GATE_EXPECT}
    gate_ok = all(v["abs_diff"] <= GATE_TOL for v in dev.values())
    ctx.update({"geometry": geo, "d50": d50,
                "repro": {"expect": GATE_EXPECT, "got": got, "deviation": dev,
                          "tol": GATE_TOL, "pass": bool(gate_ok)}})
    write_json(PREFLIGHT, ctx)
    if not gate_ok:
        bad = ", ".join("%s %.4f vs %.3f" % (k, v["got"], v["expected"])
                        for k, v in dev.items() if v["abs_diff"] > GATE_TOL)
        cp.unit_done("stage0 BLOCKED: reproduction gate failed (%s)" % bad)
        raise Blocker("STAGE 0 reproduction gate failed: %s" % bad, ctx)
    cp.unit_done("stage0 (preflight OK, geometry unique, repro gate PASS)")
    return ctx


# ---------------------------------------------------------------------------
# break identification -- the Nikora step, fixed before measuring
# ---------------------------------------------------------------------------
def break_search(d, s, d50, lo_mult=BREAK_BAND[0], hi_mult=BREAK_BAND[1],
                 min_pts=MIN_SEG_POINTS):
    """Brute-force two-segment log-log split over ladder lags in
    [lo_mult*D50, hi_mult*D50]; the split with the smallest total residual
    wins.  Returns the break index and the per-candidate residuals."""
    ok = np.isfinite(s) & (s > 0) & (d > 0)
    dd, ss = d[ok], s[ok]
    x, y = np.log(dd), np.log(ss)
    cand = np.where((dd >= lo_mult * d50) & (dd <= hi_mult * d50))[0]
    trials = []
    best = None
    for k in cand:
        n1, n2 = k + 1, len(dd) - k
        if n1 < min_pts or n2 < min_pts:
            continue
        p1 = np.polyfit(x[:k + 1], y[:k + 1], 1)
        p2 = np.polyfit(x[k:], y[k:], 1)
        r1 = y[:k + 1] - np.polyval(p1, x[:k + 1])
        r2 = y[k:] - np.polyval(p2, x[k:])
        rss = float(np.sum(r1 ** 2) + np.sum(r2 ** 2))
        trials.append({"index": int(k), "lag_m": float(dd[k]), "rss": rss})
        if best is None or rss < best[1]:
            best = (int(k), rss)
    if best is None:
        return None
    return {"index": best[0], "lag_m": float(dd[best[0]]),
            "lag_over_d50": float(dd[best[0]] / d50),
            "rss": best[1], "n_candidates": len(trials),
            "candidates": trials}


def single_fit(d, s, a, b):
    m = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
    if m.sum() < 3:
        return {"H": None, "n_points": int(m.sum()), "note": "too few lags"}
    p = np.polyfit(np.log(d[m]), np.log(s[m]), 1)
    pred = np.polyval(p, np.log(d[m]))
    resid = np.log(s[m]) - pred
    sst = float(np.sum((np.log(s[m]) - np.log(s[m]).mean()) ** 2))
    return {"H": float(0.5 * p[0]), "n_points": int(m.sum()),
            "lag_lo_used_m": float(d[m][0]), "lag_hi_used_m": float(d[m][-1]),
            "rms_log": float(np.sqrt(np.mean(resid ** 2))),
            "r2": (None if sst <= 0 else float(1.0 - np.sum(resid ** 2) / sst))}


def nikora_measure(d, s, d50, label):
    """The full Nikora read for one (face, direction)."""
    brk = break_search(d, s, d50)
    if brk is None:
        return {"face_direction": label, "break": None,
                "H_Nikora": None, "note": "no valid break candidate"}
    d_min = float(d[d > 0].min())
    b = brk["lag_m"]
    h_main = single_fit(d, s, d_min, b)
    out = {"face_direction": label, "break": brk,
           "H_Nikora": h_main,
           "two_band": {
               "below": single_fit(d, s, d_min, b),
               "above": single_fit(d, s, b, TWO_BAND_TOP * d50),
           },
           "sensitivity": {}}
    # the break nudged, but read on the same ladder (nearest lag), so that the
    # instrument is unchanged and only the split scale moves
    for tag, mult in (("x0.8", 0.8), ("x1.25", 1.25)):
        target = b * mult
        idx = int(np.argmin(np.abs(d - target)))
        out["sensitivity"][tag] = {
            "target_m": float(target), "used_m": float(d[idx]),
            "H_Nikora": single_fit(d, s, d_min, float(d[idx]))}
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- Nikora procedure on the three faces
# ---------------------------------------------------------------------------
def stage1(cp, ctx):
    log("STAGE1  Nikora procedure: axis cuts, break, H_Nikora, sensitivity")
    R5 = ruler()
    d50p = ctx["d50"]["patch_m"]
    faces = []

    zp = np.load(HF_PATCH).astype(np.float64)
    faces.append(("patch", zp, d50p, [(1, DM_PATCH, "x (flow)"),
                                      (0, DM_PATCH, "y (cross)")]))
    for tag, path, d50, dx_m, dy_m in (
            ("long_r1", HF_R1, ctx["d50"]["long_r1_m"], DM_LONG_X, DM_LONG_Y),
            ("long_r2", HF_R2, ctx["d50"]["long_r2_m"], DM_LONG_X, DM_LONG_Y)):
        z = np.fromfile(path, dtype="<f4").reshape(LONG_NY, LONG_NX).astype(np.float64)
        faces.append((tag, z, d50, [(1, dx_m, "x (flow)"),
                                    (0, dy_m, "y (cross)")]))

    results = {}
    for name, z, d50, dirs in faces:
        results[name] = {"d50_m": d50, "shape": list(z.shape),
                         "sigma_z_m": float(np.std(z)),
                         "skew": float(((z - z.mean()) ** 3).mean()
                                       / max(np.std(z), 1e-30) ** 3),
                         "directions": {}}
        for axis, dmax, dlabel in dirs:
            d, s = axis_cut(z, CELL, axis, dmax, R5)
            m = nikora_measure(d, s, d50, "%s %s" % (name, dlabel))
            m.update({
                "axis": int(axis), "dmax_m": dmax, "n_lags": int(d.size),
                "lag_min_m": float(d.min()), "lag_max_m": float(d.max()),
                "lag_min_over_d50": float(d.min() / d50),
                "record_size_over_d50": float(z.shape[axis] * CELL / d50),
                "cell_mm": CELL * 1000.0,
                "H_loc": [float(v) if np.isfinite(v) else None
                          for v in R5.local_H(d, s)],
                "curve": {"D_m": [float(v) for v in d],
                          "S_m2": [float(v) for v in s]},
            })
            results[name]["directions"][dlabel] = m
            h = (m.get("H_Nikora") or {}).get("H")
            br = (m.get("break") or {}).get("lag_over_d50")
            log("  %-22s H_Nikora=%s  break=%s"
                % ("%s %s" % (name, dlabel),
                   "%.3f" % h if h is not None else "n/a",
                   "%.2f D50" % br if br else "n/a"))
        del z
    cp.unit_done("stage1 (3 faces x 2 directions, Nikora break + sensitivity)")
    return results


# ---------------------------------------------------------------------------
# STAGE 2 -- synthetic fill, patch only
# ---------------------------------------------------------------------------
def priority_flood(z):
    """Barnes priority-flood depression fill, open borders.

    z is (ny, nx).  Returns F with F >= z everywhere and F flat over every
    depression.  Pure heapq; 4.2e6 cells for the 0.6 m patch.
    """
    ny, nx = z.shape
    zz = np.ascontiguousarray(z, dtype=np.float64).ravel()
    n = zz.size
    F = np.empty(n, dtype=np.float64)
    seen = bytearray(n)
    heap = []

    def seed(i):
        if seen[i]:
            return
        seen[i] = 1
        F[i] = zz[i]
        heap.append((zz[i], i))

    for j in range(nx):
        seed(j)
        seed((ny - 1) * nx + j)
    for i in range(ny):
        seed(i * nx)
        seed(i * nx + nx - 1)
    heapq.heapify(heap)

    while heap:
        level, i = heapq.heappop(heap)
        r, c = divmod(i, nx)
        if c > 0:
            nb = i - 1
            if not seen[nb]:
                seen[nb] = 1
                v = zz[nb] if zz[nb] > level else level
                F[nb] = v
                heapq.heappush(heap, (v, nb))
        if c < nx - 1:
            nb = i + 1
            if not seen[nb]:
                seen[nb] = 1
                v = zz[nb] if zz[nb] > level else level
                F[nb] = v
                heapq.heappush(heap, (v, nb))
        if r > 0:
            nb = i - nx
            if not seen[nb]:
                seen[nb] = 1
                v = zz[nb] if zz[nb] > level else level
                F[nb] = v
                heapq.heappush(heap, (v, nb))
        if r < ny - 1:
            nb = i + nx
            if not seen[nb]:
                seen[nb] = 1
                v = zz[nb] if zz[nb] > level else level
                F[nb] = v
                heapq.heappush(heap, (v, nb))
    return F.reshape(ny, nx)


def stage2(cp, ctx):
    log("STAGE2  synthetic fill (priority-flood), patch only; prediction %s"
        % PREDICTION)
    R5 = ruler()
    zp = np.load(HF_PATCH).astype(np.float64)
    d50 = ctx["d50"]["patch_m"]
    sigma0 = float(np.std(zp))

    t0 = time.time()
    F = priority_flood(zp)
    log("  priority-flood: %.1f s, filled fraction %.4f"
        % (time.time() - t0, float(np.mean(F > zp + 1e-12))))

    d0, s0 = axis_cut(zp, CELL, 1, DM_PATCH, R5)
    base_fit = {"x": R5.fit_band(d0, s0, *FIXED_BAND)}
    d0y, s0y = axis_cut(zp, CELL, 0, DM_PATCH, R5)
    base_fit["y"] = R5.fit_band(d0y, s0y, *FIXED_BAND)
    d0i, s0i = R5.curve_iso(zp, CELL, 0.20, nbin=NBIN)
    base_fit["iso"] = R5.fit_band(d0i, s0i, *FIXED_BAND)
    base_nk = {"x": nikora_measure(d0, s0, d50, "patch x"),
               "y": nikora_measure(d0y, s0y, d50, "patch y"),
               "iso": nikora_measure(d0i, s0i, d50, "patch iso")}

    out = {"prediction": PREDICTION, "sigma_z_before_m": sigma0,
           "filled_fraction": float(np.mean(F > zp + 1e-12)),
           "baseline": {"H_fit": base_fit, "H_Nikora": base_nk},
           "deltas": {}}
    curves = {"x": (1, DM_PATCH), "y": (0, DM_PATCH)}
    for delta in DELTAS:
        tag = "%.2f_sigma" % delta
        zp2 = np.minimum(F, zp + delta * sigma0)
        rec = {"delta_over_sigma": delta, "delta_m": delta * sigma0,
               "sigma_z_after_m": float(np.std(zp2)),
               "d_sigma_z_m": float(np.std(zp2) - sigma0), "H_fit": {}, "H_Nikora": {}}
        for k, (axis, dmax) in curves.items():
            d, s = axis_cut(zp2, CELL, axis, dmax, R5)
            rec["H_fit"][k] = R5.fit_band(d, s, *FIXED_BAND)
            rec["H_Nikora"][k] = nikora_measure(d, s, d50, "patch %s %s" % (k, tag))
        d, s = R5.curve_iso(zp2, CELL, 0.20, nbin=NBIN)
        rec["H_fit"]["iso"] = R5.fit_band(d, s, *FIXED_BAND)
        rec["H_Nikora"]["iso"] = nikora_measure(d, s, d50, "patch iso " + tag)
        rec["dH_fit"] = {k: (rec["H_fit"][k]["H"] - base_fit[k]["H"]) for k in base_fit}
        rec["dH_Nikora"] = {
            k: ((rec["H_Nikora"][k].get("H_Nikora") or {}).get("H") or np.nan)
            - ((base_nk[k].get("H_Nikora") or {}).get("H") or np.nan)
            for k in base_nk}
        out["deltas"][tag] = rec
        log("  delta=%s: dH_fit iso %+0.4f  dH_Nikora iso %+0.4f  dsigma %+0.5f m"
            % (tag, float(rec["dH_fit"]["iso"]),
               float(rec["dH_Nikora"]["iso"]), rec["d_sigma_z_m"]))
    cp.unit_done("stage2 (priority-flood + 3 fills, both instruments)")
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _fmt(v, nd=4):
    if v is None:
        return "n/a"
    if isinstance(v, float) and not np.isfinite(v):
        return "n/a"
    return ("%." + str(nd) + "f") % v


def build_report(ctx, stop_reason=None):
    L = []
    A = L.append
    A("# REPORT-V4-2-12 — 既存の生成床を Nikora 手順で測り直す＋合成充填テスト")
    A("")
    A("- brief: `D:\\mutaro-work\\rally-in\\BRIEF-V4-2-12.md`")
    A("- 実行者: `run_v4_2_12.py`（測定のみ。新規生成なし）")
    A("- 作成: %s UTC" % ctx["created_utc"])
    A("")
    if stop_reason:
        A("> **この実行は STAGE 0 のゲートで停止した。STAGE 1 / STAGE 2 は測定していない。**")
        A("> 理由: %s" % stop_reason)
        A("> 以下の表は「測った値」ではなく「測る前に固定した手順」と「事前確認の結果」である。")
        A("")

    # --- stage 0 ---
    A("## STAGE 0 — 絶対パス事前確認")
    A("")
    A("brief に列挙された絶対パスのみを確認した（別ツリーは探していない）。")
    A("")
    A("| 対象 | 絶対パス | 存在 | bytes | sha256（先頭16） |")
    A("|---|---|---|---|---|")
    for r in ctx["stage0"]["required"]:
        A("| %s | `%s` | %s | %s | %s |"
          % (r["label"], r["path"], "OK" if r["exists"] else "**MISSING**",
             r["bytes"] if r["bytes"] else "-",
             (r["sha256"][:16] + "…") if r["sha256"] else "-"))
    A("")
    A("生成側スクリプト: " + ", ".join(
        "%s=%s" % (g["label"], "OK" if g["exists"] else "MISSING")
        for g in ctx["stage0"]["generator_scripts"]))
    A("")
    if ctx["stage0"]["missing"]:
        A("### 欠けている入力（STOP 条件）")
        A("")
        for r in ctx["stage0"]["missing"]:
            A("- **%s**: `%s` — 存在しない。" % (r["label"], r["path"]))
        A("")
        A("brief の STAGE 0 は「1つでも欠けたら STOP して報告。別ツリーを探さない」と")
        A("定めているため、ここで停止した。STAGE 1 / STAGE 2 は実行していない。")
        A("")
        A("#### 欠落の由来（既存記録からの事実。新たな生成はしていない）")
        A("")
        A("- `out/v4-2-3/run.log` は帯状床を2本作っている（r1 seed 3255045276382562994、")
        A("  r2 seed 3255045276383062994）が、`saved ...heightfield-long-r1.npy` は")
        A("  1行だけで、r2 の保存行は無い（run.log 15行目が r1 の保存、16行目から r2 の生成）。")
        A("- `run_v4_2_5.py` 495行目は r2 を「its heightfield was not written to disk,")
        A("  so no fresh recompute is possible」と明記している。")
        A("- `out/v4-2-3/longbed-meta-2.json` は r2 のメタデータ（`role: replication`）で、")
        A("  床は作られたが場が保存されなかった、という上記2点と整合する。")
        A("- したがって r2 の場は既存ツリーのどこにも無い（brief の指示どおり別ツリーは")
        A("  探していない）。復元には再生成が必要で、それは brief の「新規生成なし」に反する。")
        A("")

    probe = ctx["stage0"].get("probe") or {}
    if probe:
        A("### 入力の形式・形状・格子（実測）")
        A("")
        A("| 入力 | 形式 | shape | 格子間隔 | 物理サイズ | 備考 |")
        A("|---|---|---|---|---|---|")
        p = probe.get("patch")
        if p:
            A("| パッチ床 | %s %s | %s | %.8f mm | %.3f × %.3f m | σz %.6f m, 歪度 %.3f, D50 %.6f m |"
              % (p["container"], p["dtype"], p["shape"], p["cell_m_from_generator"] * 1000,
                 p["size_m"][0], p["size_m"][1], p["sigma_z_m"], p["skew"],
                 p["D50_surface_m_meta"] or float("nan")))
        for tag, label in (("long_r1", "帯状 r1"), ("long_r2", "帯状 r2")):
            q = probe.get(tag) or {}
            if q.get("present"):
                A("| %s | %s | %s | %.8f mm | %.3f × %.3f m | bytes %s/%s %s |"
                  % (label, q["container"], q["shape_expected"], CELL * 1000,
                     LONG_NX * CELL, LONG_NY * CELL, q["bytes"], q["bytes_expected"],
                     "OK" if q["bytes_ok"] else "**不一致**"))
            else:
                A("| %s | %s | %s | %.8f mm | %.3f × %.3f m | **ファイル無し**（期待 %d bytes） |"
                  % (label, q.get("container", "headerless raw float32"), q.get("shape_expected"),
                     CELL * 1000, LONG_NX * CELL, LONG_NY * CELL, q.get("bytes_expected")))
        A("")
        r1m, r2m = probe.get("long_r1_meta"), probe.get("long_r2_meta")
        if r1m and r2m:
            A("- 帯状メタ: r1 %.3f × %.3f m, D50 %.6f m, σz %.6f m / r2 %.3f × %.3f m, "
              "D50 %.6f m, σz %.6f m (role=%s)"
              % (r1m["length_m"], r1m["width_m"], r1m["surface_D50_m"], r1m["sigma_z_m"],
                 r2m["length_m"], r2m["width_m"], r2m["surface_D50_m"], r2m["sigma_z_m"],
                 r2m["role"]))
            A("")

    geo = ctx["stage0"].get("geometry")
    d50 = ctx["stage0"].get("d50")
    if geo:
        A("### 形状・格子・物理サイズ（生成側から確認）")
        A("")
        A("| 床 | shape | 格子間隔 | 物理サイズ | 根拠 |")
        A("|---|---|---|---|---|")
        A("| パッチ | %s | %.8f mm | %.3f × %.3f m | %s |"
          % (geo["patch"]["shape"], geo["patch"]["cell_m"] * 1000,
             geo["patch"]["size_m"][0], geo["patch"]["size_m"][1],
             geo["patch"]["size_source"]))
        A("| 帯状 r1/r2 | %s | %.8f mm | %.3f × %.3f m | %s |"
          % (geo["long_common"]["shape_expected"], geo["long_common"]["cell_m"] * 1000,
             geo["long_common"]["size_m"][0], geo["long_common"]["size_m"][1],
             geo["long_common"]["size_source"]))
        A("")
        A("- 一意性: **%s**" % ("OK" if geo["unique"] else "NG"))
    if d50:
        A("- 表面 D50: パッチ %.6f m / 帯状 r1 %.6f m / 帯状 r2 %.6f m"
          % (d50["patch_m"], d50["long_r1_m"], d50["long_r2_m"]))
    A("")
    rp = ctx["stage0"].get("repro")
    A("### 再現ゲート（v4-2.5 の既報値、±%.2f）" % GATE_TOL)
    A("")
    if not rp:
        A("**未実行**（入力欠落のため）。")
    else:
        A("| 面・方向 | 既報 | 実測 | |差| | 判定 |")
        A("|---|---|---|---|---|")
        for k, v in rp["deviation"].items():
            A("| %s | %.3f | %s | %s | %s |"
              % (k, v["expected"], _fmt(v["got"]), _fmt(v["abs_diff"]),
                 "PASS" if v["abs_diff"] <= GATE_TOL else "FAIL"))
        A("")
        A("再現ゲート: **%s**" % ("PASS" if rp["pass"] else "FAIL"))
    A("")

    # --- stage 1 ---
    A("## STAGE 1 — Nikora 手順での測り直し")
    A("")
    A("手順（測定前に固定。途中で変えない）:")
    A("")
    A("1. 原フィールドそのまま（detrending なし）。軸カット D(lx,0) / D(0,ly) を")
    A("   60-bin 対数ラダーで計算。うちの生成場は周期場なので `np.roll` の")
    A("   ラップはその前提で正当。")
    A("2. 折れの同定: log–log の2区間最小二乗。候補はラダー点のうち")
    A("   [0.25 × D50, 2 × D50]。各区間 4 点以上。総残差最小の折れを採る。")
    A("3. 指数: H_Nikora = [最小ラグ, 折れ] の単一フィット。二帯域 =")
    A("   [最小ラグ, 折れ] と [折れ, 4 × D50]。感度 = 折れ ×0.8 / ×1.25。")
    A("   H_loc(ℓ) は全ラダー値。")
    A("")
    if not ctx.get("stage1"):
        A("**未測定**（STAGE 0 で停止したため）。この節には実測値を書かない。")
        A("")
    else:
        A("| 面 | 方向 | 格子 mm | 記録サイズ ×D50 | 最小ラグ /D50 | 折れ mm | 折れ ×D50 | H_Nikora | 帯1 | 帯2 | ×0.8 | ×1.25 |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for face, fd in ctx["stage1"].items():
            for dl, m in fd["directions"].items():
                nk = m.get("H_Nikora") or {}
                tb = m.get("two_band") or {}
                se = m.get("sensitivity") or {}
                br = m.get("break") or {}
                A("| %s | %s | %.4f | %.1f | %.3f | %s | %s | %s | %s | %s | %s | %s |"
                  % (face, dl, m["cell_mm"], m["record_size_over_d50"],
                     m["lag_min_over_d50"],
                     _fmt((br.get("lag_m") or 0) * 1000, 3),
                     _fmt(br.get("lag_over_d50"), 2),
                     _fmt(nk.get("H"), 3),
                     _fmt(((tb.get("below") or {}).get("H")), 3),
                     _fmt(((tb.get("above") or {}).get("H")), 3),
                     _fmt(((se.get("x0.8") or {}).get("H_Nikora") or {}).get("H"), 3),
                     _fmt(((se.get("x1.25") or {}).get("H_Nikora") or {}).get("H"), 3)))
        A("")

    # --- stage 2 ---
    A("## STAGE 2 — 合成充填テスト（filleting の先行テスト）")
    A("")
    A("事前予測（測定前に記録）: **δ = 0.5 σz の充填で ΔH_fit ≤ +%.2f**。"
      % PREDICTION["dH_fit_max"])
    A("")
    if not ctx.get("stage2"):
        A("**未測定**（STAGE 0 で停止したため）。")
        A("")
    else:
        s2 = ctx["stage2"]
        A("priority-flood 充填率 %.4f、σz(充填前) = %.6f m。" % (s2["filled_fraction"], s2["sigma_z_before_m"]))
        A("")
        A("| δ/σz | δ (mm) | Δσz (mm) | ΔH_fit x | ΔH_fit y | ΔH_fit iso | ΔH_Nikora iso | 予測照合 |")
        A("|---|---|---|---|---|---|---|---|")
        for tag, r in s2["deltas"].items():
            A("| %s | %.3f | %+.4f | %+0.4f | %+0.4f | %+0.4f | %+0.4f | %s |"
              % (tag, r["delta_m"] * 1000, r["d_sigma_z_m"] * 1000,
                 r["dH_fit"]["x"], r["dH_fit"]["y"], r["dH_fit"]["iso"],
                 r["dH_Nikora"]["iso"],
                 ("該当" if tag.startswith("0.50") else "—")))
        A("")

    # --- disclosures ---
    A("## 開示事項")
    A("")
    A("- 測定のみ。既存ファイルは変更していない（`run_v4_2_12.py` と `out/v4-2-12/` のみ）。")
    A("- `np.roll` は周期場前提。軸カットの最大ラグはドメインの半分以下に取った。")
    A("- detrending なし（brief の指定）。したがって form（低波数）も H に入る。")
    A("- 折れ探索はラダー点の総当たりで、候補は [0.25, 2] × D50。折れの位置は")
    A("  ラダーの粒度（60 bin）に量子化される。感度 ×0.8 / ×1.25 は最近傍の")
    A("  ラダー点で読み、使用した実ラグも numbers に残した。")
    A("- STAGE 2 の充填は開境界の priority-flood（Barnes 型）。場は周期場だが、")
    A("  充填アルゴリズムは開領域で定義されるため境界は開で計算した。")
    A("- `run_v4_2_5.py` の関数は import して verbatim で使った（複製していない）。")
    A("")
    A("## チェックポイント")
    A("")
    st = ctx.get("checkpoint") or {}
    A("- units=%s, last=%s, updated=%s"
      % (st.get("completed_units"), st.get("last_unit"), st.get("updated_at")))
    A("")
    A("## ログ")
    A("")
    A("```")
    L.extend(LOG[-40:])
    A("```")
    A("")
    return "\n".join(L) + "\n"


def write_report(ctx, stop_reason=None):
    text = build_report(ctx, stop_reason)
    atomic_write_bytes(REPORT, text.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 0)
    log("  wrote %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="BRIEF-V4-2-12 measurement job")
    ap.add_argument("--stage0", action="store_true", help="preflight only")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    created = utcnow()
    t_run = time.time()

    ctx = {"created_utc": created, "brief": BRIEF, "this_file": THIS_FILE,
           "deliverables": {"report": REPORT, "numbers": NUMBERS,
                            "preflight": PREFLIGHT, "checkpoint": CP_PATH}}

    try:
        ctx["stage0"] = stage0(cp)
    except Blocker as e:
        ctx["stage0"] = e.ctx
        ctx["stage1"] = None
        ctx["stage2"] = None
        ctx["status"] = "STOPPED_AT_STAGE0"
        ctx["blocker"] = e.reason
        ctx["timings_s"] = {"total": time.time() - t_run}
        # one jobkit unit for the stop record; the checkpoint is captured after
        # the unit is committed so the report and the json agree with it
        cp.unit_done("report (STOP record: stopped at stage0 gate)")
        ctx["checkpoint"] = cp.state()
        numbers = {
            "brief": BRIEF, "created_utc": created, "status": "STOPPED_AT_STAGE0",
            "deliverables_complete": False, "blocker": e.reason,
            "stage0": ctx["stage0"], "stage1": None, "stage2": None,
            "contract": {"cell_m": CELL, "nbin": NBIN, "fixed_band_m": FIXED_BAND,
                         "break_band_over_d50": BREAK_BAND, "two_band_top_over_d50": TWO_BAND_TOP,
                         "prediction": PREDICTION, "gate_tol": GATE_TOL},
            "checkpoint": ctx["checkpoint"], "timings_s": ctx["timings_s"], "log": LOG,
        }
        write_json(NUMBERS, numbers)
        write_report(ctx, stop_reason=e.reason)
        log("STOP  %s" % e.reason)
        print("\n".join(LOG[-6:]))
        return 2

    ctx["stage1"] = stage1(cp, ctx["stage0"])
    ctx["stage2"] = stage2(cp, ctx["stage0"])
    ctx["status"] = "COMPLETE"
    ctx["timings_s"] = {"total": time.time() - t_run}
    cp.unit_done("report (REPORT-V4-2-12.md + numbers-v4-2-12.json)")
    ctx["checkpoint"] = cp.state()
    numbers = {
        "brief": BRIEF, "created_utc": created, "status": "COMPLETE",
        "deliverables_complete": True,
        "stage0": ctx["stage0"], "stage1": ctx["stage1"], "stage2": ctx["stage2"],
        "contract": {"cell_m": CELL, "nbin": NBIN, "fixed_band_m": FIXED_BAND,
                     "break_band_over_d50": BREAK_BAND,
                     "two_band_top_over_d50": TWO_BAND_TOP,
                     "prediction": PREDICTION, "gate_tol": GATE_TOL},
        "checkpoint": ctx["checkpoint"], "timings_s": ctx["timings_s"], "log": LOG,
    }
    write_json(NUMBERS, numbers)
    write_report(ctx)
    log("DONE  %s" % REPORT)
    log("      %s" % NUMBERS)
    print("\n".join(LOG[-6:]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
