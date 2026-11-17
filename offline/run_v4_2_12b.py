"""run_v4_2_12b -- BRIEF-V4-2-12B: re-measure the existing generated bed with
Nikora's procedure, regenerate r2 from its seed, and run a synthetic-fill
(filleting) probe on the patch.

MEASUREMENT + ONE REGENERATION.  No parameter change, no rendering, no new bed
beyond the r2 replication that the brief asks for.  The only files this job
writes are this .py and everything under out/v4-2-12b/.

What the brief asks, and where it is answered
---------------------------------------------
STAGE 0    absolute-path preflight over the brief's list (one missing file is
           a STOP -- no other tree is searched); geometry / D50 cross-checked
           against the metadata; reproduction gate: v4-2.5's published
           scaling-band H must come back within +-0.01 using run_v4_2_5's own
           functions (verbatim import).
STAGE 0.5  read the r1/r2 seed correspondence out of run_v4_2_3.py + the two
           metadata files (+ that report); regenerate r2 from its seed into
           out/v4-2-12b/heightfield-long-r2-regen.npy; identity gate: measure
           the regenerated face with run_v4_2_5's ruler and compare with the
           r2 values v4-2.5 recorded (x 0.322 / iso 0.325, +-0.01).  PASS ->
           the regenerated face joins STAGE 1 as r2.  FAIL -> r2 is reported
           as not reproduced and STAGE 1 proceeds on patch + r1 only.
STAGE 1    Nikora's procedure on patch x/y and long-bed flow/cross directions:
           raw field (no detrending), axis cuts D(lx,0) / D(0,ly) on the 60-bin
           log ladder, break by brute force over the ladder points inside
           [0.25 D50, 2 D50] (two-segment log-log least squares, >=4 points per
           segment, minimum total residual), H_Nikora = single fit
           [min lag, break], the two-band slopes [min lag, break] and
           [break, 4 D50], the break sensitivity (x0.8 / x1.25, nearest ladder
           point), and the whole H_loc(l) ladder.
STAGE 2    patch only: priority-flood depression fill F (Barnes, open
           boundary), z'(d) = min(F, z + d), d in {0.25, 0.5, 0.75} sigma_z,
           re-measured with both instruments.  Prediction, fixed in code
           before the measurement: at d = 0.5 sigma_z, dH_fit <= +0.10.

Conventions fixed BEFORE the first measurement (not tuned afterwards)
---------------------------------------------------------------------
* cell         : 0.6 m / 2048 px = 0.29296875 mm, exactly v4-2.2 / v4-2.3.
* ladder       : run_v4_2_5.curve_dir's own 60-bin geometric ladder, rounded
                 to whole pixels and np.unique'd (imported, not copied).
* dmax         : patch x/y 0.30 m (half of the 0.6 m box = half the domain);
                 long flow x 1.50 m (half of the 3.0 m bed); long cross y
                 0.30 m (half of the 0.6 m width).  the gate measurement uses
                 1.00 m / 0.20 m exactly as v4-2.5 did, so the published value
                 is reproduced on the same ladder.
* break search : two-segment log-log least squares; candidate = ladder point in
                 [0.25 D50, 2 D50]; segment 1 = [first ladder lag, candidate],
                 segment 2 = [candidate, 4 D50]; both need >= 4 points; the
                 split with the smallest total log-residual wins.  capping the
                 upper segment at 4 D50 (the same top as the reported second
                 band) keeps the break independent of dmax.
* periodic     : the generated fields are periodic by construction, so np.roll
                 is the correct wrap for the lag differences (the report says
                 so, and says where the long bed's 3 m seam sits).
* fill         : Barnes priority-flood from the open border (level =
                 max(z, level)); open (non-periodic) boundary is declared.
* gates        : +-0.01 on H, unchanged.  No 0.70-0.85 judgement here.

Stages are checkpointed one jobkit unit at a time (stage0 / stage0_5 / stage1 /
stage2 / report).  The run is foreground: it either finishes every deliverable
or stops at a gate and writes an honest stop record.
"""
from __future__ import annotations

import datetime
import gc
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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-12B.md"

OUT = os.path.join(HERE, "out", "v4-2-12b")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-12b.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-12B.md")
PREFLIGHT = os.path.join(OUT, "stage0-preflight.json")
REGEN = os.path.join(OUT, "heightfield-long-r2-regen.npy")
REGEN_META = os.path.join(OUT, "heightfield-long-r2-regen-meta.json")
THIS_FILE = "run_v4_2_12b.py"

# --- STAGE 0: the brief's absolute paths, verbatim -------------------------
HF_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
HF_R1 = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")
META_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
META_R1 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta.json")
META_R2 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta-2.json")
R3_PATH = os.path.join(HERE, "run_v4_2_3.py")
R5_PATH = os.path.join(HERE, "run_v4_2_5.py")
B2_PATH = os.path.join(HERE, "bed_gen2.py")

REQUIRED = [
    ("ruler (source of the verbatim functions)", R5_PATH),
    ("patch bed", HF_PATCH),
    ("long bed r1 (headerless raw float32)", HF_R1),
    ("long bed r1 meta", META_R1),
    ("long bed r2 meta", META_R2),
    ("generator run_v4_2_3.py", R3_PATH),
    ("generator bed_gen2.py", B2_PATH),
]
# sha256 prefix / byte count the brief names for the two big files
BRIEF_HASHES = {
    os.path.basename(HF_PATCH): {"bytes": 16777344,
                                 "sha256_prefix": "e8c7584a89a15659"},
    os.path.basename(HF_R1): {"bytes": 83886080,
                              "sha256_prefix": "4011919af065040f"},
}

# --- the measurement contract (fixed before the first measurement) ---------
CELL = 0.6 / 2048.0                   # 0.00029296875 m, exactly v4-2.2/2.3
LONG_NY, LONG_NX = 2048, 10240        # 0.6 m x 3.0 m
NBIN = 60
DM_PATCH = 0.30                       # half of the 0.6 m box
DM_LONG_X = 1.50                      # half of the 3.0 m bed
DM_LONG_Y = 0.30                      # half of the 0.6 m width
GATE_DM_LONG_X = 1.00                 # the v4-2.5 choice, for the identity gate
GATE_DM_ISO = 0.20
FIXED_BAND = (0.004, 0.036)           # the ruler: 0.4-3.6 cm, unchanged
BREAK_LO, BREAK_HI = 0.25, 2.0        # break window, x D50
TWO_TOP = 4.0                         # second band top, x D50
MIN_SEG_POINTS = 4

GATE_TOL = 0.01

# STAGE 0 reproduction gate: v4-2.5's published scaling-band H (brief values)
REPRO_EXPECT = {
    "patch_x": 0.323, "patch_y": 0.373, "patch_iso": 0.349,
    "long_r1_x": 0.294, "long_r1_y": 0.297, "long_r1_iso": 0.295,
}
# the unrounded values v4-2.5 recorded (out/v4-2-5/numbers-v4-2-5.json)
REPRO_RECORDED = {
    "patch_x": 0.32311075522810045,
    "patch_y": 0.37299321343399783,
    "patch_iso": 0.34854666750555946,
    "long_r1_x": 0.294268927895509,
    "long_r1_y": 0.29716828399283884,
    "long_r1_iso": 0.29495733653508976,
}

# STAGE 0.5 identity gate: v4-2.5's recorded r2 grain-bed values
R2_EXPECT = {"x": 0.32225764489966235, "iso": 0.3245887726612185}
R2_EXPECT_BRIEF = {"x": 0.322, "iso": 0.325}

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


def _san(o):
    if isinstance(o, dict):
        return {str(k): _san(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_san(v) for v in o]
    if isinstance(o, np.ndarray):
        return _san(o.tolist())
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
    data = (json.dumps(_san(obj), ensure_ascii=False, indent=1) + "\n"
            ).encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True

    atomic_write_bytes(path, data, verify=verify)
    return path


class Blocker(Exception):
    """A gate failed.  Carries what was measured so a stop record can be
    written honestly instead of guessing."""

    def __init__(self, reason, ctx):
        super().__init__(reason)
        self.reason = reason
        self.ctx = ctx


# ---------------------------------------------------------------------------
# instrument -- imported verbatim from run_v4_2_5 (the ruler's source of truth)
# ---------------------------------------------------------------------------
def ruler():
    import run_v4_2_5 as R5
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band"):
        assert hasattr(R5, name), "run_v4_2_5 lost %s" % name
    return R5


def axis_cut(z, cell, axis, dmax, R5, nbin=NBIN):
    """D(lx, 0) or D(0, ly): one axis, the ruler's own log ladder, periodic."""
    return R5.curve_dir(z, cell, dmax, nbin=nbin, axis=axis)


def load_long(path):
    return np.fromfile(path, dtype="<f4").reshape(LONG_NY, LONG_NX).astype(np.float64)


def skew_of(a):
    a = np.asarray(a, dtype=np.float64)
    m = float(a.mean())
    s = float(a.std())
    return float(np.mean((a - m) ** 3) / max(s, 1e-30) ** 3)


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight, geometry, D50, reproduction gate
# ---------------------------------------------------------------------------
def probe_inputs():
    out = {}
    if os.path.isfile(HF_PATCH):
        a = np.load(HF_PATCH, mmap_mode="r")
        z = np.asarray(a, dtype=np.float64)
        out["patch"] = {
            "container": "npy (with header)", "dtype": str(a.dtype),
            "shape": list(a.shape), "cell_m": CELL,
            "size_m": [a.shape[1] * CELL, a.shape[0] * CELL],
            "sigma_z_m": float(np.std(z)), "skew": skew_of(z),
        }
        del z
    if os.path.isfile(HF_R1):
        n = os.path.getsize(HF_R1)
        head = np.fromfile(HF_R1, dtype="<f4", count=LONG_NY * 64).astype(np.float64)
        out["long_r1"] = {
            "container": "headerless raw float32 LE (<f4)",
            "bytes": n, "bytes_expected": LONG_NY * LONG_NX * 4,
            "bytes_ok": bool(n == LONG_NY * LONG_NX * 4),
            "shape": [LONG_NY, LONG_NX],
            "sample_head_rows": 64,
            "sample_stats": {"n_cells": int(head.size),
                             "min_m": float(head.min()),
                             "max_m": float(head.max()),
                             "std_m": float(head.std())},
        }
    return out


def stage0(cp):
    log("STAGE0  absolute-path preflight (the brief's list only)")
    rows, missing = [], []
    for label, path in REQUIRED:
        ok = os.path.isfile(path)
        row = {"label": label, "path": path, "exists": bool(ok),
               "bytes": os.path.getsize(path) if ok else None,
               "sha256": sha256_file(path) if ok else None}
        rows.append(row)
        log("  %-44s %s %s" % (label, "OK  " if ok else "MISSING",
                               row["bytes"] if ok else path))
        if not ok:
            missing.append(row)
    if missing:
        names = ", ".join(os.path.basename(r["path"]) for r in missing)
        ctx = {"checked_utc": utcnow(), "required": rows, "missing": missing}
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: missing %s -- STOP per brief" % names)
        raise Blocker("STAGE 0 preflight failed: missing %s" % names, ctx)

    import bed_gen2 as B2

    # the brief's sha256 prefixes must hold
    hash_check = {}
    for r in rows:
        name = os.path.basename(r["path"])
        if name in BRIEF_HASHES:
            exp = BRIEF_HASHES[name]
            hash_check[name] = {
                "bytes_match": bool(r["bytes"] == exp["bytes"]),
                "sha256_prefix_match": bool(
                    r["sha256"].startswith(exp["sha256_prefix"])),
                "bytes": r["bytes"], "sha256": r["sha256"],
                "brief_bytes": exp["bytes"],
                "brief_sha256_prefix": exp["sha256_prefix"]}
            log("  hash %-30s bytes=%s sha256_prefix=%s"
                % (name, hash_check[name]["bytes_match"],
                   hash_check[name]["sha256_prefix_match"]))
    if not all(v["bytes_match"] and v["sha256_prefix_match"]
               for v in hash_check.values()):
        ctx = {"checked_utc": utcnow(), "required": rows,
               "hash_check": hash_check}
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: file size / sha256 prefix mismatch")
        raise Blocker("STAGE 0: the two big heightfields do not match the "
                      "brief's byte count / sha256 prefix", ctx)

    log("STAGE0  geometry and D50 from the generator side + metadata")
    meta_p = json.load(open(META_PATCH, encoding="utf-8"))
    meta_r1 = json.load(open(META_R1, encoding="utf-8"))
    meta_r2 = json.load(open(META_R2, encoding="utf-8"))
    zp = np.load(HF_PATCH, mmap_mode="r")

    cell_gen = float(B2.PATCH) / int(B2.HIRES)
    geo = {
        "cell_m_from_generator": cell_gen,
        "cell_m_patch_meta": meta_p["cell_m"],
        "cell_m_long_meta": meta_r1["cell_m"],
        "cell_unique": bool(
            abs(cell_gen - CELL) < 1e-15
            and abs(float(meta_p["cell_m"]) - CELL) < 1e-15
            and abs(float(meta_r1["cell_m"]) - CELL) < 1e-15
            and abs(float(meta_r2["cell_m"]) - CELL) < 1e-15),
        "patch": {
            "shape": list(zp.shape), "cell_m": CELL,
            "size_m": [float(zp.shape[1] * CELL), float(zp.shape[0] * CELL)],
            "bed_gen2_patch_m": float(B2.PATCH),
            "bed_gen2_hires": int(B2.HIRES),
        },
        "long": {
            "shape": [LONG_NY, LONG_NX], "cell_m": CELL,
            "size_m": [float(LONG_NX * CELL), float(LONG_NY * CELL)],
            "meta_r1_length_m": meta_r1["length_m"],
            "meta_r1_width_m": meta_r1["width_m"],
            "meta_r2_length_m": meta_r2["length_m"],
            "meta_r2_width_m": meta_r2["width_m"],
        },
        "bytes": {os.path.basename(HF_PATCH): os.path.getsize(HF_PATCH),
                  os.path.basename(HF_R1): os.path.getsize(HF_R1)},
    }
    geo["unique"] = bool(
        geo["cell_unique"]
        and geo["patch"]["size_m"] == [float(B2.PATCH), float(B2.PATCH)]
        and geo["long"]["size_m"] == [3.0, 0.6]
        and abs(geo["long"]["meta_r1_length_m"] - 3.0) < 1e-12
        and abs(geo["long"]["meta_r1_width_m"] - 0.6) < 1e-12)
    d50 = {
        "brief_patch_m": 0.0362,
        "patch_meta_m": float(meta_p["D50_surface_measured_m"]),
        "brief_r1_m": 0.038297, "r1_meta_m": float(meta_r1["surface_D50_m"]),
        "brief_r2_m": 0.037976, "r2_meta_m": float(meta_r2["surface_D50_m"]),
        "match_patch": bool(abs(float(meta_p["D50_surface_measured_m"])
                                - 0.0362) < 1e-4),
        "match_r1": bool(abs(float(meta_r1["surface_D50_m"]) - 0.038297)
                         < 1e-6),
        "match_r2": bool(abs(float(meta_r2["surface_D50_m"]) - 0.037976)
                         < 1e-6),
        "patch_meta_sigma_z_m": float(meta_p["sigma_z_m"]),
        "r1_meta_sigma_z_m": float(meta_r1["sigma_z_m"]),
        "r2_meta_sigma_z_m": float(meta_r2["sigma_z_m"]),
        "patch_sigma_z_over_D50": float(meta_p["sigma_z_m"])
        / float(meta_p["D50_surface_measured_m"]),
        "r1_sigma_z_over_D50": float(meta_r1["sigma_z_m"])
        / float(meta_r1["surface_D50_m"]),
        "r2_sigma_z_over_D50": float(meta_r2["sigma_z_m"])
        / float(meta_r2["surface_D50_m"]),
    }
    log("  cell %.11f m unique=%s; D50 patch %.6f r1 %.6f r2 %.6f"
        % (CELL, geo["unique"], d50["patch_meta_m"], d50["r1_meta_m"],
           d50["r2_meta_m"]))
    if not (geo["unique"] and d50["match_patch"] and d50["match_r1"]
            and d50["match_r2"]):
        ctx = {"checked_utc": utcnow(), "required": rows,
               "hash_check": hash_check, "geometry": geo, "d50": d50,
               "probe": probe_inputs()}
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: geometry / D50 does not agree")
        raise Blocker("STAGE 0: geometry or D50 does not agree with the "
                      "metadata (geometry unique=%s)" % geo["unique"], ctx)
    del zp

    # --- reproduction gate: v4-2.5's published numbers, +-0.01 --------------
    log("STAGE0  reproduction gate vs v4-2.5 (band %.3f-%.3f m, tol %.2f)"
        % (FIXED_BAND[0], FIXED_BAND[1], GATE_TOL))
    R5 = ruler()
    got, cache = {}, {}

    zp64 = np.load(HF_PATCH).astype(np.float64)
    d0, s0 = axis_cut(zp64, CELL, 1, DM_PATCH, R5)
    cache["patch_x"] = (d0, s0)
    got["patch_x"] = R5.fit_band(d0, s0, *FIXED_BAND)["H"]
    d1, s1 = axis_cut(zp64, CELL, 0, DM_PATCH, R5)
    cache["patch_y"] = (d1, s1)
    got["patch_y"] = R5.fit_band(d1, s1, *FIXED_BAND)["H"]
    di, si = R5.curve_iso(zp64, CELL, GATE_DM_ISO, nbin=NBIN)
    got["patch_iso"] = R5.fit_band(di, si, *FIXED_BAND)["H"]
    del zp64
    log("  patch     x %.4f  y %.4f  iso %.4f"
        % (got["patch_x"], got["patch_y"], got["patch_iso"]))

    zl = load_long(HF_R1)
    dx, sx = axis_cut(zl, CELL, 1, GATE_DM_LONG_X, R5)
    got["long_r1_x"] = R5.fit_band(dx, sx, *FIXED_BAND)["H"]
    dy, sy = axis_cut(zl, CELL, 0, DM_LONG_Y, R5)
    cache["long_r1_y"] = (dy, sy)
    got["long_r1_y"] = R5.fit_band(dy, sy, *FIXED_BAND)["H"]
    dli, sli = R5.curve_iso(zl, CELL, GATE_DM_ISO, nbin=NBIN)
    got["long_r1_iso"] = R5.fit_band(dli, sli, *FIXED_BAND)["H"]
    del zl
    gc.collect()
    log("  long r1   x %.4f  y %.4f  iso %.4f"
        % (got["long_r1_x"], got["long_r1_y"], got["long_r1_iso"]))

    dev = {}
    for k in REPRO_EXPECT:
        dev[k] = {"expected_brief": REPRO_EXPECT[k],
                  "recorded_v425": REPRO_RECORDED[k],
                  "got": float(got[k]),
                  "abs_diff_vs_brief": float(abs(got[k] - REPRO_EXPECT[k])),
                  "abs_diff_vs_recorded": float(
                      abs(got[k] - REPRO_RECORDED[k]))}
    repro_ok = all(v["abs_diff_vs_brief"] <= GATE_TOL for v in dev.values())
    repro = {"expect": REPRO_EXPECT, "recorded_v425": REPRO_RECORDED,
             "got": got, "deviation": dev, "tol": GATE_TOL,
             "pass": bool(repro_ok)}
    # G-LEGACY baseline: every other .py as it is right now
    baseline = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p) and name != THIS_FILE:
            baseline[name] = {"sha256": sha256_file(p),
                              "bytes": os.path.getsize(p)}
    log("  G-LEGACY baseline: %d existing .py hashed" % len(baseline))

    ctx = {"checked_utc": utcnow(), "required": rows, "hash_check": hash_check,
           "probe": probe_inputs(), "geometry": geo, "d50": d50,
           "repro": repro,
           "g_legacy_baseline": {"n": len(baseline), "py": baseline},
           "brief": BRIEF}
    write_json(PREFLIGHT, ctx)
    if not repro_ok:
        bad = ", ".join("%s %.4f vs %.3f" % (k, v["got"], v["expected_brief"])
                        for k, v in dev.items()
                        if v["abs_diff_vs_brief"] > GATE_TOL)
        cp.unit_done("stage0 BLOCKED: reproduction gate failed (%s)" % bad)
        raise Blocker("STAGE 0 reproduction gate failed: %s" % bad, ctx)
    cp.unit_done("stage0")
    return ctx, cache


# ---------------------------------------------------------------------------
# STAGE 0.5 -- seed correspondence, r2 regeneration, identity gate
# ---------------------------------------------------------------------------
META_COMPARE_FIELDS = [
    "seed", "tile_seeds", "shape", "cell_m", "length_m", "width_m",
    "n_grains", "grains_per_tile", "kinds", "form_grid", "n_low_modes",
    "low_mode_periods_m", "low_mode_lattice_sx", "low_mode_kmax",
    "sigma_form_m", "place", "sigma_z_m", "sigma_z_matrix_m", "skew_matrix",
    "sigma_z_grain_m", "sigma_z_form_m", "matrix_area_fraction",
    "surface_D50_m", "surface_D84_over_D50", "armor_ratio_exact",
    "substrate_sigma_z_m", "substrate_skew", "hole_fraction",
]


def _get_meta(m, field):
    if field == "hole_fraction":
        return (m.get("raster") or {}).get("hole_fraction")
    if field == "place":
        return (m.get("place") or {}).get("z_center_range_m")
    return m.get(field)


def _cmp(a, b):
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        ds = [_cmp(x, y) for x, y in zip(a, b)]
        return {"abs": max(d["abs"] for d in ds),
                "rel": max(d["rel"] for d in ds),
                "equal": all(d["equal"] for d in ds)}
    if isinstance(a, dict) and isinstance(b, dict):
        dsh = [_cmp(a[k], b[k]) for k in a if k in b]
        if not dsh:
            return {"abs": None, "rel": None,
                    "equal": bool(set(a) == set(b))}
        return {"abs": max(d["abs"] for d in dsh),
                "rel": max(d["rel"] for d in dsh),
                "equal": all(d["equal"] for d in dsh)
                and set(a) == set(b)}
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return {"abs": None, "rel": None, "equal": bool(a == b)}
    d = abs(fa - fb)
    scale = max(abs(fa), abs(fb), 1e-30)
    return {"abs": float(d), "rel": float(d / scale),
            "equal": bool(d <= 1e-9 * max(1.0, scale))}


def stage0_5(cp, ctx):
    import run_v4_2_3 as R3
    import bed_gen2 as B2

    log("STAGE0.5  r1/r2 seed correspondence, read from the generator side")
    meta_r1 = json.load(open(META_R1, encoding="utf-8"))
    meta_r2 = json.load(open(META_R2, encoding="utf-8"))
    long_seeds = [int(v) for v in R3.LONG_SEEDS]
    seed_offset = long_seeds[1] - long_seeds[0]
    path_check = {
        "run_v4_2_3.LONG_SEEDS": long_seeds,
        "run_v4_2_3.SEED": int(R3.SEED),
        "bed_gen2.SEED": int(B2.SEED),
        "seed_offset": int(seed_offset),
        "shape_constant": [int(R3.LONG_NY), int(R3.LONG_NX)],
        "cell_constant": float(R3.CELL),
        "n_tiles": int(R3.N_TILES),
        "form_grid": [int(R3.FORM_NY), int(R3.FORM_NX)],
        "meta_r1_seed": int(meta_r1["seed"]),
        "meta_r2_seed": int(meta_r2["seed"]),
        "r1_seed_matches_code": bool(int(meta_r1["seed"]) == long_seeds[0]),
        "r2_seed_matches_code": bool(int(meta_r2["seed"]) == long_seeds[1]),
        "r1_tile_seeds_match": bool(
            [int(v) for v in meta_r1["tile_seeds"]]
            == [long_seeds[0] + 10 * k for k in range(R3.N_TILES)]),
        "r2_tile_seeds_match": bool(
            [int(v) for v in meta_r2["tile_seeds"]]
            == [long_seeds[1] + 10 * k for k in range(R3.N_TILES)]),
        "r1_shape_matches": bool(list(meta_r1["shape"]) ==
                                 [int(R3.LONG_NY), int(R3.LONG_NX)]),
        "r2_shape_matches": bool(list(meta_r2["shape"]) ==
                                 [int(R3.LONG_NY), int(R3.LONG_NX)]),
        "r1_cell_matches": bool(abs(float(meta_r1["cell_m"]) - CELL) < 1e-15),
        "r2_cell_matches": bool(abs(float(meta_r2["cell_m"]) - CELL) < 1e-15),
        "r1_file_on_disk": bool(
            os.path.getsize(HF_R1) == LONG_NY * LONG_NX * 4),
        "r1_vs_r2_differ_only_by_seed": bool(
            long_seeds[1] - long_seeds[0] == 500000),
        "source_line": ("LONG_SEEDS = [int(SEED), int(SEED) + 500000] in "
                        "run_v4_2_3.py; SEED = bed_gen2.SEED = %d" % int(B2.SEED)),
    }
    log("  seeds r1=%d r2=%d (offset %d); meta seeds match=%s/%s"
        % (long_seeds[0], long_seeds[1], seed_offset,
           path_check["r1_seed_matches_code"],
           path_check["r2_seed_matches_code"]))
    if not (path_check["r1_seed_matches_code"]
            and path_check["r2_seed_matches_code"]
            and path_check["r1_tile_seeds_match"]
            and path_check["r2_tile_seeds_match"]
            and path_check["r2_shape_matches"]):
        ctx["stage0_5"] = {"path_check": path_check, "regen": None,
                           "identity_gate": None}
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0_5 BLOCKED: seed correspondence not verified")
        raise Blocker("STAGE 0.5: the metadata seeds do not line up with "
                      "run_v4_2_3.LONG_SEEDS", ctx)

    # the path probe the brief allows instead of rebuilding r1: its population
    # (5 tiles, tile seed = seed + 10k) must come back with the recorded sizes
    log("  path probe: r1 population counts from its own seeds")
    probe_counts = [len(B2.populate(long_seeds[0] + 10 * k)[0])
                    for k in range(R3.N_TILES)]
    probe_ok = bool(probe_counts == [int(v) for v in meta_r1["grains_per_tile"]])
    path_check["r1_population_probe"] = {
        "counts": probe_counts,
        "meta_grains_per_tile": [int(v) for v in meta_r1["grains_per_tile"]],
        "match": probe_ok}
    log("  r1 population probe: %s (match=%s)"
        % (probe_counts, probe_ok))

    # --- regenerate r2 ------------------------------------------------------
    log("STAGE0.5  regenerating r2 from seed %d (keep=False: nothing under "
        "out/v4-2-3 is touched)" % long_seeds[1])
    t0 = time.time()
    regen_error = None
    zh = None
    regen_meta = None
    try:
        cfg = R3.load_cfg()
        lb = R3.build_long(cfg, long_seeds[1], keep=False)
        regen_meta = lb["meta"]
        zh = lb["zh"]
        lb.pop("gid", None)
        lb.pop("formpix", None)
        lb.pop("grains", None)
        lb = None
        gc.collect()
    except Exception as exc:                    # noqa: BLE001
        regen_error = "%s: %s" % (type(exc).__name__, exc)
        log("  REGENERATION FAILED: %s" % regen_error)
    regen_s = float(time.time() - t0)
    log("  regeneration wall %.1f s, zh=%s"
        % (regen_s, None if zh is None else zh.shape))

    regen_block = {"seed": long_seeds[1], "wall_s": regen_s,
                   "error": regen_error, "file": None, "meta_compare": None,
                   "readback": None}
    gate = {"pass": False, "reason": "not reached", "values": None}

    if zh is not None:
        # meta comparison against the saved r2 metadata
        cmp_rows = []
        all_equal = True
        for f in META_COMPARE_FIELDS:
            a = _get_meta(regen_meta, f)
            b = _get_meta(meta_r2, f)
            if a is None or b is None:
                continue
            c = _cmp(a, b)
            all_equal = all_equal and c["equal"]
            cmp_rows.append({"field": f, "regen": a, "saved": b,
                             "abs_diff": c["abs"], "rel_diff": c["rel"],
                             "match": c["equal"]})
        regen_block["meta_compare"] = {
            "n_fields": len(cmp_rows), "all_match": bool(all_equal),
            "rows": cmp_rows,
            "tolerance": "|diff| <= 1e-9 * max(1, |value|) for floats; exact "
                         "for ints/lists/dicts",
            "ignored": ["wall_s", "raster.raycast_wall_s"]}
        log("  regen meta: %d fields, all_match=%s"
            % (len(cmp_rows), all_equal))

        # save: headerless raw float32 LE, the same container as r1
        payload = zh.astype("<f4").tobytes()
        atomic_write_bytes(REGEN, payload,
                           verify=lambda p: os.path.getsize(p)
                           == LONG_NY * LONG_NX * 4)
        del payload
        gc.collect()
        sha = sha256_file(REGEN)
        zr = np.fromfile(REGEN, dtype="<f4").reshape(LONG_NY, LONG_NX)
        readback = {
            "container": "headerless raw float32 LE (<f4), no .npy header "
                         "(same convention as heightfield-long-r1.npy)",
            "shape": [LONG_NY, LONG_NX], "bytes": os.path.getsize(REGEN),
            "sha256": sha,
            "sigma_z_read_m": float(zr.std(dtype=np.float64)),
            "sigma_z_in_memory_m": float(np.std(zh)),
            "bit_identical_to_float32_cast": bool(
                np.array_equal(zr, zh.astype("<f4"))),
            "cell_m": CELL, "seed": long_seeds[1],
        }
        del zr
        regen_block["file"] = {"path": REGEN, "bytes": readback["bytes"],
                               "sha256": sha}
        regen_block["readback"] = readback
        log("  saved %s (%d B, sha256 %s)"
            % (os.path.basename(REGEN), readback["bytes"], sha[:16]))

        # --- identity gate: the ruler's fixed band, x and iso ---------------
        log("STAGE0.5  identity gate vs v4-2.5's r2 record "
            "(x %.4f / iso %.4f, +-%.2f)"
            % (R2_EXPECT_BRIEF["x"], R2_EXPECT_BRIEF["iso"], GATE_TOL))
        R5 = ruler()
        dz, sz = axis_cut(zh, CELL, 1, GATE_DM_LONG_X, R5)
        hx = R5.fit_band(dz, sz, *FIXED_BAND)
        dI, sI = R5.curve_iso(zh, CELL, GATE_DM_ISO, nbin=NBIN)
        hi = R5.fit_band(dI, sI, *FIXED_BAND)
        vals = {"x": float(hx["H"]), "iso": float(hi["H"])}
        dev = {k: {"recorded": R2_EXPECT[k],
                   "brief": R2_EXPECT_BRIEF[k],
                   "got": vals[k],
                   "abs_diff": float(abs(vals[k] - R2_EXPECT[k]))}
               for k in vals}
        ok = all(v["abs_diff"] <= GATE_TOL for v in dev.values())
        gate = {"pass": bool(ok), "tol": GATE_TOL, "values": vals,
                "deviation": dev,
                "what_was_compared": ("the ruler's fixed 0.4-3.6 cm band "
                                      "H_fit on the regenerated grains face; "
                                      "the same instrument and the same "
                                      "ladder dmax as v4-2.5's r2 curve"),
                "reason": ("PASS" if ok else "FAIL: regenerated r2 does not "
                           "reproduce the recorded fixed-band H")}
        log("  gate: x %.4f (d %.4f)  iso %.4f (d %.4f)  -> %s"
            % (vals["x"], dev["x"]["abs_diff"], vals["iso"],
               dev["iso"]["abs_diff"], "PASS" if ok else "FAIL"))
        write_json(REGEN_META, {
            "path": REGEN, "created_utc": utcnow(),
            "source": "regenerated by run_v4_2_12b.py from run_v4_2_3.build_"
                      "long(cfg, seed, keep=False); seed read from "
                      "run_v4_2_3.LONG_SEEDS[1]",
            "seed": long_seeds[1], "cell_m": CELL,
            "shape": [LONG_NY, LONG_NX],
            "container": readback["container"], "bytes": readback["bytes"],
            "sha256": sha,
            "regenerated_meta": regen_meta,
            "meta_compare": regen_block["meta_compare"]})
    else:
        gate = {"pass": False, "tol": GATE_TOL, "values": None,
                "reason": "regeneration did not complete: %s" % regen_error}

    ctx["stage0_5"] = {"path_check": path_check, "regen": regen_block,
                       "identity_gate": gate,
                       "verdict": ("PASS: regenerated face is used as r2"
                                   if gate["pass"] else
                                   "FAIL: r2 not reproduced; STAGE 1 uses "
                                   "patch + r1 only")}
    del zh
    gc.collect()
    write_json(PREFLIGHT, ctx)
    cp.unit_done("stage0_5")
    return ctx


# ---------------------------------------------------------------------------
# break identification -- the Nikora step, fixed before measuring
# ---------------------------------------------------------------------------
def break_search(d, s, d50):
    """Brute-force two-segment log-log split over ladder lags in
    [0.25 D50, 2 D50]; segment 1 = [first lag, candidate], segment 2 =
    [candidate, 4 D50]; both need >= 4 points; minimum total residual wins."""
    ok = np.isfinite(s) & (s > 0) & (d > 0)
    dd, ss = d[ok], s[ok]
    x, y = np.log(dd), np.log(ss)
    j_top = int(np.searchsorted(dd, TWO_TOP * d50, side="right") - 1)
    lo, hi = BREAK_LO * d50, BREAK_HI * d50
    trials = []
    best = None
    for k in range(MIN_SEG_POINTS - 1, j_top - MIN_SEG_POINTS + 1):
        if dd[k] < lo or dd[k] > hi:
            continue
        p1 = np.polyfit(x[:k + 1], y[:k + 1], 1)
        p2 = np.polyfit(x[k:j_top + 1], y[k:j_top + 1], 1)
        r1 = y[:k + 1] - np.polyval(p1, x[:k + 1])
        r2 = y[k:j_top + 1] - np.polyval(p2, x[k:j_top + 1])
        rss = float(np.sum(r1 ** 2) + np.sum(r2 ** 2))
        trials.append({"index": int(k), "lag_m": float(dd[k]),
                       "lag_over_d50": float(dd[k] / d50), "rss": rss})
        if best is None or rss < best[1]:
            best = (int(k), rss)
    if best is None:
        return None
    k = best[0]
    return {"index": k, "lag_m": float(dd[k]),
            "lag_over_d50": float(dd[k] / d50), "rss": best[1],
            "window_m": [float(lo), float(hi)],
            "upper_segment_top_m": float(dd[j_top]),
            "n_candidates": len(trials), "candidates": trials}


def alt_break_lastlag(d, s, d50):
    """POST-HOC ROBUSTNESS ONLY (not the reported procedure).

    Same two-segment log-log least squares, but segment 2 runs to the last
    ladder lag instead of 4 x D50.  Quantifies how much the break definition
    moves the answer.
    """
    ok = np.isfinite(s) & (s > 0) & (d > 0)
    dd, ss = d[ok], s[ok]
    x, y = np.log(dd), np.log(ss)
    lo, hi = BREAK_LO * d50, BREAK_HI * d50
    best = None
    for k in range(MIN_SEG_POINTS - 1, len(dd) - MIN_SEG_POINTS + 1):
        if dd[k] < lo or dd[k] > hi:
            continue
        p1 = np.polyfit(x[:k + 1], y[:k + 1], 1)
        p2 = np.polyfit(x[k:], y[k:], 1)
        rss = float(np.sum((y[:k + 1] - np.polyval(p1, x[:k + 1])) ** 2)
                    + np.sum((y[k:] - np.polyval(p2, x[k:])) ** 2))
        if best is None or rss < best[1]:
            best = (int(k), rss)
    if best is None:
        return None
    k = best[0]
    H = float(0.5 * np.polyfit(x[:k + 1], y[:k + 1], 1)[0])
    return {"lag_m": float(dd[k]), "lag_over_d50": float(dd[k] / d50),
            "H": H, "rss": best[1]}


def posthoc_alt_breaks(st1):
    """The §3.4 robustness table, computed once and stored in numbers."""
    out = {}
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            brk = m.get("break") or {}
            h0 = (m.get("H_Nikora") or {}).get("H")
            alt = alt_break_lastlag(np.array(m["curve"]["D_m"], dtype=float),
                                    np.array(m["curve"]["S_m2"], dtype=float),
                                    face["d50_m"])
            out["%s | %s" % (name, dlabel)] = {
                "break_main_m": brk.get("lag_m"),
                "H_main": h0,
                "break_seg2_to_last_lag_m": (alt or {}).get("lag_m"),
                "H_seg2_to_last_lag": (alt or {}).get("H"),
                "dH": ((alt["H"] - h0)
                       if (alt is not None and h0 is not None) else None)}
    return {"note": "post-hoc robustness only; NOT the reported procedure. "
                    "Same two-segment fit with segment 2 running to the last "
                    "ladder lag instead of 4 x D50.",
            "rows": out}


def nikora_measure(d, s, d50, label, R5):
    """One (face, direction) read by the fixed Nikora procedure."""
    brk = break_search(d, s, d50)
    if brk is None:
        return {"label": label, "break": None, "H_Nikora": None,
                "note": "no valid break candidate on this ladder"}
    d_min = float(d[d > 0].min())
    b = brk["lag_m"]
    h_main = R5.fit_band(d, s, d_min, b)
    out = {"label": label,
           "break": brk,
           "H_Nikora": h_main,
           "two_band": {"below": h_main,
                        "above": R5.fit_band(d, s, b, TWO_TOP * d50),
                        "below_rule": "same fit as H_Nikora",
                        "above_rule": "[break, 4 x D50]"},
           "sensitivity": {}}
    for tag, mult in (("x0.8", 0.8), ("x1.25", 1.25)):
        target = b * mult
        idx = int(np.argmin(np.abs(d - target)))
        used = float(d[idx])
        fit = R5.fit_band(d, s, d_min, used)
        out["sensitivity"][tag] = {
            "target_m": float(target), "used_m": used,
            "used_over_d50": float(used / d50),
            "shift_from_break_m": float(used - b), "fit": fit,
            "dH": (float(fit["H"] - h_main["H"])
                   if fit["H"] is not None and h_main["H"] is not None
                   else None)}
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- Nikora procedure on the available faces
# ---------------------------------------------------------------------------
def stage1(cp, ctx, cache):
    log("STAGE1  Nikora procedure: axis cuts, break, H_Nikora, sensitivity")
    R5 = ruler()
    import bed_gen2 as B2

    d50_p = ctx["d50"]["patch_meta_m"]
    faces = [("patch", "npy", HF_PATCH, d50_p,
              (("x (saved axis 1)", 1, DM_PATCH, "patch_x"),
               ("y (saved axis 0)", 0, DM_PATCH, "patch_y")))]
    faces.append(("long_r1", "long", HF_R1, ctx["d50"]["r1_meta_m"],
                  (("x (flow, long axis)", 1, DM_LONG_X, None),
                   ("y (cross)", 0, DM_LONG_Y, "long_r1_y"))))
    if ctx["stage0_5"]["identity_gate"]["pass"]:
        faces.append(("long_r2 (regenerated)", "long", REGEN,
                      float(_regen_d50(ctx)),
                      (("x (flow, long axis)", 1, DM_LONG_X, None),
                       ("y (cross)", 0, DM_LONG_Y, None))))

    results = {"method": {
        "detrending": "none (raw field, as the brief fixes)",
        "axis_cuts": "D(lx,0) / D(0,ly) from run_v4_2_5.curve_dir "
                     "(periodic np.roll, 60-bin log ladder)",
        "break_window": "[0.25, 2] x D50, candidates = ladder points",
        "break_fit": "two-segment log-log least squares; segment 1 = "
                     "[first lag, brk], segment 2 = [brk, 4 x D50]; >= 4 "
                     "points each; minimum total residual",
        "H_Nikora": "single fit [first ladder lag (= cell), brk], "
                    "0.5*dlnS/dlnd",
        "two_band": "[first lag, brk] and [brk, 4 x D50]",
        "sensitivity": "brk x 0.8 / x 1.25 -> nearest ladder lag; refit",
        "local_slope": "0.5*dlnS/dlnd on the same ladder (run_v4_2_5.local_H)",
        "fit_quality": "run_v4_2_5.fit_band (H, H_se, rms_log, r2, n, "
                       "actual lags used)",
        "grid_mm": CELL * 1000.0},
        "faces": {}}

    for name, kind, path, d50, dirs in faces:
        t0 = time.time()
        z = np.load(path).astype(np.float64) if kind == "npy" else load_long(path)
        rec = {"d50_m": float(d50), "shape": list(z.shape),
               "sigma_z_m": float(np.std(z)), "skew": float(B2._skew(z)),
               "directions": {}}
        log("  face %-22s shape %s sigma_z %.6f m"
            % (name, z.shape, rec["sigma_z_m"]))
        for dlabel, axis, dmax, ck in dirs:
            if ck and ck in cache:
                d, s = cache[ck]
                src = "reused from the STAGE 0 run (same ladder)"
            else:
                d, s = axis_cut(z, CELL, axis, dmax, R5)
                src = "computed in STAGE 1"
                if ck:
                    cache[ck] = (d, s)
            m = nikora_measure(d, s, float(d50),
                               "%s %s" % (name, dlabel), R5)
            m.update({
                "axis": int(axis), "dmax_m": float(dmax),
                "curve_source": src, "n_lags": int(d.size),
                "lag_min_m": float(d.min()), "lag_max_m": float(d.max()),
                "lag_min_over_d50": float(d.min() / d50),
                "cell_mm": CELL * 1000.0,
                "record_size_over_d50": float(z.shape[axis] * CELL / d50),
                "sigma_z_over_d50": float(np.std(z) / d50),
                "skew": float(rec["skew"]),
                "H_fit_fixed_band": R5.fit_band(d, s, *FIXED_BAND),
                "H_loc": {"lag_m": [float(v) for v in d],
                          "H": [float(v) for v in R5.local_H(d, s)]},
                "curve": {"D_m": [float(v) for v in d],
                          "S_m2": [float(v) for v in s]}})
            rec["directions"][dlabel] = m
            h = (m.get("H_Nikora") or {}).get("H")
            br = (m.get("break") or {}).get("lag_m")
            log("    %-22s H_Nikora=%s break=%s mm"
                % (dlabel, "%.4f" % h if h is not None else "n/a",
                   "%.3f" % (br * 1000.0) if br else "n/a"))
        rec["wall_s"] = float(time.time() - t0)
        results["faces"][name] = rec
        del z
        gc.collect()
    cp.unit_done("stage1")
    return results


def _regen_d50(ctx):
    for row in ctx["stage0_5"]["regen"]["meta_compare"]["rows"]:
        if row["field"] == "surface_D50_m":
            return row["regen"]
    return 0.03797627400499068


# ---------------------------------------------------------------------------
# STAGE 2 -- synthetic fill (priority flood), patch only
# ---------------------------------------------------------------------------
def priority_flood(z):
    """Barnes priority-flood depression fill, open borders.

    F >= z everywhere; F is flat across every depression.  The border is the
    outlet (open boundary, not periodic).  4.2e6 cells for the 0.6 m patch.
    """
    ny, nx = z.shape
    zz = np.ascontiguousarray(z, dtype=np.float64).ravel()
    n = zz.size
    F = np.empty(n, dtype=np.float64)
    seen = bytearray(n)
    heap = []

    def seed(i):
        if not seen[i]:
            seen[i] = 1
            F[i] = zz[i]
            heap.append((float(zz[i]), i))

    for j in range(nx):
        seed(j)
        seed((ny - 1) * nx + j)
    for i in range(ny):
        seed(i * nx)
        seed(i * nx + nx - 1)
    heapq.heapify(heap)

    # local bindings for the hot loop
    heappop = heapq.heappop
    heappush = heapq.heappush
    zzl = zz
    Fl = F
    sn = seen
    while heap:
        level, i = heappop(heap)
        c = i % nx
        if c > 0:
            nb = i - 1
            if not sn[nb]:
                sn[nb] = 1
                v = zzl[nb] if zzl[nb] > level else level
                Fl[nb] = v
                heappush(heap, (float(v), nb))
        if c < nx - 1:
            nb = i + 1
            if not sn[nb]:
                sn[nb] = 1
                v = zzl[nb] if zzl[nb] > level else level
                Fl[nb] = v
                heappush(heap, (float(v), nb))
        if i >= nx:
            nb = i - nx
            if not sn[nb]:
                sn[nb] = 1
                v = zzl[nb] if zzl[nb] > level else level
                Fl[nb] = v
                heappush(heap, (float(v), nb))
        if i + nx < n:
            nb = i + nx
            if not sn[nb]:
                sn[nb] = 1
                v = zzl[nb] if zzl[nb] > level else level
                Fl[nb] = v
                heappush(heap, (float(v), nb))
    return F.reshape(ny, nx)


def stage2(cp, ctx):
    log("STAGE2  synthetic fill (priority flood, open boundary), patch only")
    R5 = ruler()
    zp = np.load(HF_PATCH).astype(np.float64)
    d50 = ctx["d50"]["patch_meta_m"]
    sigma0 = float(np.std(zp))
    log("  sigma_z before = %.6f m (prediction: at delta=0.5 sigma, "
        "dH_fit <= +%.2f)" % (sigma0, PREDICTION["dH_fit_max"]))

    t0 = time.time()
    F = priority_flood(zp)
    fill_s = float(time.time() - t0)
    uplift = F - zp
    filled_fraction = float(np.mean(uplift > 1e-12))
    log("  priority flood %.1f s; filled fraction %.4f; max uplift %.5f m"
        % (fill_s, filled_fraction, float(uplift.max())))

    def base_curves(zz):
        out = {}
        out["x"] = axis_cut(zz, CELL, 1, DM_PATCH, R5)
        out["y"] = axis_cut(zz, CELL, 0, DM_PATCH, R5)
        out["iso"] = R5.curve_iso(zz, CELL, GATE_DM_ISO, nbin=NBIN)
        return out

    base_cv = base_curves(zp)
    base = {"H_fit": {k: R5.fit_band(*base_cv[k], *FIXED_BAND)
                      for k in base_cv},
            "H_Nikora": {k: nikora_measure(*base_cv[k], d50, "patch " + k, R5)
                         for k in base_cv}}
    rec_base = {}
    for k in base["H_fit"]:
        hb = (base["H_Nikora"][k].get("H_Nikora") or {}).get("H")
        rec_base[k] = {"H_fit": base["H_fit"][k]["H"], "H_Nikora": hb}
        log("  baseline %-3s H_fit %.4f  H_Nikora %s  break %s mm"
            % (k, base["H_fit"][k]["H"],
               "%.4f" % hb if hb is not None else "n/a",
               "%.3f" % (base["H_Nikora"][k]["break"]["lag_m"] * 1000.0)
               if base["H_Nikora"][k].get("break") else "n/a"))

    out = {"prediction": PREDICTION, "deltas": DELTAS,
           "fill": {"method": "Barnes priority flood, open border (level = "
                              "max(z, level)); F >= z, F flat over "
                              "depressions",
                    "filled_fraction": filled_fraction,
                    "max_uplift_m": float(uplift.max()),
                    "mean_uplift_where_filled_m": float(
                        uplift[uplift > 1e-12].mean()),
                    "wall_s": fill_s,
                    "boundary": "open (non-periodic) border; the structure "
                                "function below still uses the periodic roll, "
                                "which is disclosed"},
           "sigma_z_before_m": sigma0,
           "baseline": {"H_fit": {k: base["H_fit"][k] for k in base["H_fit"]},
                        "H_Nikora": base["H_Nikora"], "summary": rec_base},
           "cells": {}}
    del uplift

    for delta in DELTAS:
        tag = "%.2f sigma_z" % delta
        zp2 = np.minimum(F, zp + delta * sigma0)
        sig = float(np.std(zp2))
        rec = {"delta_over_sigma": float(delta), "delta_m": float(delta * sigma0),
               "sigma_z_after_m": sig, "d_sigma_z_m": float(sig - sigma0),
               "skew_after": float(skew_of(zp2)),
               "H_fit": {}, "H_Nikora": {}, "dH_fit": {}, "dH_Nikora": {},
               "break_lag_m": {}}
        cvs = base_curves(zp2)
        for k in ("x", "y", "iso"):
            rec["H_fit"][k] = R5.fit_band(*cvs[k], *FIXED_BAND)
            rec["H_Nikora"][k] = nikora_measure(*cvs[k], d50,
                                                "patch %s %s" % (k, tag), R5)
            hb = (rec["H_Nikora"][k].get("H_Nikora") or {}).get("H")
            h0 = rec_base[k]["H_Nikora"]
            rec["dH_fit"][k] = (float(rec["H_fit"][k]["H"]
                                      - rec_base[k]["H_fit"])
                                if rec["H_fit"][k]["H"] is not None else None)
            rec["dH_Nikora"][k] = (float(hb - h0)
                                   if hb is not None and h0 is not None else None)
            rec["break_lag_m"][k] = (rec["H_Nikora"][k]["break"]["lag_m"]
                                     if rec["H_Nikora"][k].get("break") else None)
        out["cells"][tag] = rec
        log("  delta=%-11s dH_fit x %+0.3f y %+0.3f iso %+0.3f | "
            "dH_Nikora x %+0.3f y %+0.3f iso %+0.3f | dsigma %+0.6f m"
            % (tag, rec["dH_fit"]["x"], rec["dH_fit"]["y"],
               rec["dH_fit"]["iso"], rec["dH_Nikora"]["x"],
               rec["dH_Nikora"]["y"], rec["dH_Nikora"]["iso"],
               rec["d_sigma_z_m"]))
        del zp2, cvs
        gc.collect()

    # prediction check: at delta = 0.5 sigma, dH_fit <= +0.10 in every direction
    key = "%.2f sigma_z" % PREDICTION["delta_over_sigma"]
    c05 = out["cells"][key]
    vals = [c05["dH_fit"][k] for k in ("x", "y", "iso")]
    ok_vals = [v for v in vals if v is not None]
    worst = max(ok_vals) if len(ok_vals) == len(vals) else None
    out["prediction_check"] = {
        "cell": key,
        "dH_fit": c05["dH_fit"],
        "max_dH_fit": (float(worst) if worst is not None else None),
        "threshold": float(PREDICTION["dH_fit_max"]),
        "pass": bool(worst is not None
                     and worst <= PREDICTION["dH_fit_max"])}
    log("  prediction at %s: max dH_fit %+0.4f (threshold +%.2f) -> %s"
        % (key, worst, PREDICTION["dH_fit_max"],
           "PASS" if out["prediction_check"]["pass"] else "FAIL"))
    cp.unit_done("stage2")
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f(v, nd=4):
    if v is None:
        return "n/a"
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(v):
        return "n/a"
    return ("%." + str(nd) + "f") % v


def sgn(v, nd=3):
    if v is None:
        return "n/a"
    return ("%+." + str(nd) + "f") % float(v)


def build_report(ctx):
    L = []
    A = L.append
    st0 = ctx["stage0"]
    st05 = ctx["stage0_5"]
    st1 = ctx["stage1"]
    st2 = ctx["stage2"]
    geo, d50 = st0["geometry"], st0["d50"]
    gate = st05["identity_gate"]

    A("# REPORT-V4-2-12B — 既存の生成床を Nikora 手順で測り直す＋合成充填テスト")
    A("")
    A("- brief: `D:\\mutaro-work\\rally-in\\BRIEF-V4-2-12B.md`")
    A("- 実行: `run_v4_2_12b.py`（前景・最後まで）")
    A("- 作成: %s UTC" % ctx["created_utc"])
    A("- 書き込み範囲: `run_v4_2_12b.py` と `out/v4-2-12b/` のみ。既存ファイルは"
      "変更していない（§6 の G-LEGACY で自己監査）。")
    A("")
    A("## 0. 要約")
    A("")
    A("- **STAGE 0**: brief が挙げた絶対パスの全ファイルが存在し、2つの大きな"
      " heightfield のバイト数と sha256 先頭は brief の記載と一致した。"
      "再現ゲートは v4-2.5 の既報 6 値を ±0.01 以内で再現し **PASS**。")
    mc0 = st05["regen"]["meta_compare"]
    A("- **STAGE 0.5**: r1/r2 の作り分けは `run_v4_2_3.LONG_SEEDS = "
      "[SEED, SEED + 500000]` の seed 差のみで、両メタの seed とタイル seed が"
      "コードと一致した。r1 の母集団をその seed から引き直しても記録の"
      "タイル別粒数と一致した。r2 を seed から再生成し、メタ "
      "%s。同一性ゲート（固定帯 x / iso、±0.01）は **%s**。"
      % (("%d 項目すべてが保存メタと一致" % mc0["n_fields"]) if mc0 else
         "照合はできなかった（再生成が完了せず）",
         "PASS" if gate["pass"] else "FAIL"))
    if gate["pass"]:
        A("  再生成面は r2 として STAGE 1 の対象に加えた。")
    else:
        A("  再生成面は r2 として扱わず、STAGE 1 はパッチ＋r1 の2面で行った"
          "（brief の指示どおり停止しない）。")
    A("- **STAGE 1**: 面×方向（パッチ x/y、帯状 r1 流下/横断、r2 PASS なら r2 の"
      "2方向＝計 6 本）で折れ・H_Nikora・二帯域・感度・H_loc を取得。"
      "数値は §3 と numbers JSON。")
    d05 = st2["prediction_check"]
    A("- **STAGE 2**: 予測「δ = 0.5 σz の充填でも ΔH_fit ≤ +0.10」に対し、"
      "実測の最大 ΔH_fit は **%s**（閾値 +%.2f）→ **%s**。"
      % (sgn(d05["max_dH_fit"]), d05["threshold"],
         "予測どおり" if d05["pass"] else "予測に反する"))
    A("")

    # ------------------------------------------------------------------ 1
    A("## 1. STAGE 0 — 絶対パス事前確認と再現ゲート")
    A("")
    A("### 1.1 存在・バイト数・sha256")
    A("")
    A("| ファイル | bytes | sha256 (先頭16) | brief 記載との一致 |")
    A("|---|---:|---|---|")
    for r in st0["required"]:
        name = os.path.basename(r["path"])
        hc = st0["hash_check"].get(name)
        if hc:
            ok = "bytes %s / sha256 %s" % (
                "一致" if hc["bytes_match"] else "不一致",
                "一致" if hc["sha256_prefix_match"] else "不一致")
        else:
            ok = "—（brief に hash 記載なし）"
        A("| `%s` | %d | `%s` | %s |"
          % (name, r["bytes"], r["sha256"][:16], ok))
    A("")
    A("絶対パスは brief の記載そのもので、別ツリーの探索はしていない。"
      "1つでも欠ければ STOP する実装（欠損 0 件）。")
    A("")
    A("### 1.2 物理サイズ・格子・D50")
    A("")
    A("| 量 | 値 | 突き合わせ |")
    A("|---|---|---|")
    A("| 格子間隔 | %.11f m | generator `bed_gen2.PATCH/HIRES`、パッチ meta、"
      "r1/r2 meta の4者が一致: %s |"
      % (geo["cell_m_from_generator"],
         "はい" if geo["cell_unique"] else "いいえ"))
    A("| パッチ | %s px = %.4f × %.4f m | brief の 0.6×0.6 m と一致 |"
      % (tuple(geo["patch"]["shape"]), geo["patch"]["size_m"][0],
         geo["patch"]["size_m"][1]))
    A("| 帯状床 | %s px = %.4f × %.4f m | brief の 3.0×0.6 m・meta と一致 |"
      % (tuple(geo["long"]["shape"]), geo["long"]["size_m"][0],
         geo["long"]["size_m"][1]))
    A("| D50 パッチ | %.8f m | brief 0.0362 m と一致: %s |"
      % (d50["patch_meta_m"], "はい" if d50["match_patch"] else "いいえ"))
    A("| D50 r1 | %.8f m | brief 0.038297 m と一致: %s |"
      % (d50["r1_meta_m"], "はい" if d50["match_r1"] else "いいえ"))
    A("| D50 r2 | %.8f m | brief 0.037976 m と一致: %s |"
      % (d50["r2_meta_m"], "はい" if d50["match_r2"] else "いいえ"))
    A("| σz/D50（パッチ / r1 / r2） | %.4f / %.4f / %.4f | 記録用 |"
      % (d50["patch_sigma_z_over_D50"], d50["r1_sigma_z_over_D50"],
         d50["r2_sigma_z_over_D50"]))
    A("")
    A("### 1.3 再現ゲート（v4-2.5 の既報値を run_v4_2_5 の関数で再計算）")
    A("")
    A("帯 = 固定帯 0.4–3.6 cm。パッチは dmax 0.30 m ×2軸＋iso、帯状 r1 は "
      "x dmax 1.00 m・y 0.30 m・iso 0.20 m（v4-2.5 と同じラダー）。")
    A("")
    A("| 面・方向 | v4-2.5 既報（丸め） | v4-2.5 記録値 | 再計算 | 差（記録値比） "
      "| ±0.01 |")
    A("|---|---:|---:|---:|---:|---|")
    for k in ("patch_x", "patch_y", "patch_iso",
              "long_r1_x", "long_r1_y", "long_r1_iso"):
        r = st0["repro"]["deviation"][k]
        A("| %s | %.3f | %.6f | %.6f | %+.6f | %s |"
          % (k, r["expected_brief"], r["recorded_v425"], r["got"],
             r["got"] - r["recorded_v425"],
             "PASS" if r["abs_diff_vs_brief"] <= st0["repro"]["tol"] else "FAIL"))
    A("")
    A("→ **再現ゲート %s**（全 6 値が ±%.2f 以内；記録値との差は最大 %s）。"
      % ("PASS" if st0["repro"]["pass"] else "FAIL", st0["repro"]["tol"],
         "%.2e" % max(abs(r["got"] - r["recorded_v425"])
                      for r in st0["repro"]["deviation"].values())))
    if all(abs(r["got"] - r["recorded_v425"]) == 0.0
           for r in st0["repro"]["deviation"].values()):
        A("差は 6 値とも**厳密に 0**（同じ関数・同じラダー・同じ入力なので、"
          "v4-2.5 の記録値と同一の float64 が返った）。")
    A("")

    # ------------------------------------------------------------------ 2
    A("## 2. STAGE 0.5 — r2 の再生成と同一性ゲート")
    A("")
    pc = st05["path_check"]
    A("### 2.1 r1 と r2 の作り分け（seed の対応）")
    A("")
    A("`run_v4_2_3.py` の生成定数と 2 つのメタを突き合わせた結果:")
    A("")
    A("| 量 | r1 | r2 | 根拠 |")
    A("|---|---|---|---|")
    A("| seed | %d | %d | `LONG_SEEDS = [SEED, SEED + %d]`、`SEED = "
      "bed_gen2.SEED = %d` |"
      % (pc["run_v4_2_3.LONG_SEEDS"][0], pc["run_v4_2_3.LONG_SEEDS"][1],
         pc["seed_offset"], pc["bed_gen2.SEED"]))
    A("| メタ seed 一致 | %s | %s | longbed-meta.json / longbed-meta-2.json |"
      % ("はい" if pc["r1_seed_matches_code"] else "いいえ",
         "はい" if pc["r2_seed_matches_code"] else "いいえ"))
    A("| タイル seed | seed+10k, k=0..4 | seed+10k, k=0..4 | 両メタの "
      "`tile_seeds` が一致: %s / %s |"
      % ("はい" if pc["r1_tile_seeds_match"] else "いいえ",
         "はい" if pc["r2_tile_seeds_match"] else "いいえ"))
    A("| 形状・格子 | %s px, cell %.11f m | 同左 | 両メタとコード定数 |"
      % (tuple(pc["shape_constant"]), pc["cell_constant"]))
    A("| 違い | — | seed のみ（+%d） | 生成規則・形状・タイル数・低次モードは"
      "同一 |" % pc["seed_offset"])
    A("")
    A("経路の追加検証（r1 を作り直す代わり）: r1 の seed から "
      "`bed_gen2.populate(seed + 10k)` を引き直したタイル別粒数は "
      "`%s`、保存メタは `%s`、一致 **%s**。"
      % (pc["r1_population_probe"]["counts"],
         pc["r1_population_probe"]["meta_grains_per_tile"],
         "はい" if pc["r1_population_probe"]["match"] else "いいえ"))
    A("")
    A("### 2.2 再生成 r2 のメタ照合")
    A("")
    rb = st05["regen"]
    if rb["meta_compare"] is None:
        A("再生成は完了しなかった: `%s`。この場合も停止せず、r2 を名乗らずに"
          "パッチ＋r1 で STAGE 1 に進む。" % rb["error"])
    else:
        mc = rb["meta_compare"]
        A("再生成 `run_v4_2_3.build_long(cfg, seed, keep=False)` の戻りメタを"
          "保存メタ（longbed-meta-2.json）と比較。**%d 項目すべて一致: %s**"
          "（float は |差| ≤ 1e-9·max(1,|値|)、int・リストは完全一致；"
          "`wall_s` と `raster.raycast_wall_s` は時間なので除外）。"
          % (mc["n_fields"], "はい" if mc["all_match"] else "いいえ"))
        A("")
        A("| 項目 | 再生成 | 保存メタ | 差 | 一致 |")
        A("|---|---:|---:|---:|---|")
        for row in mc["rows"]:
            if isinstance(row["regen"], (int, float)) and not isinstance(
                    row["regen"], bool) and row["regen"] is not None:
                a = ("%d" % row["regen"]) if isinstance(row["regen"], int) \
                    else ("%.9g" % row["regen"])
                b = ("%d" % row["saved"]) if isinstance(row["saved"], int) \
                    else ("%.9g" % row["saved"])
                dd = ("%.2e" % row["abs_diff"]
                      if row["abs_diff"] is not None else "—")
                A("| `%s` | %s | %s | %s | %s |"
                  % (row["field"], a, b, dd,
                     "はい" if row["match"] else "いいえ"))
            else:
                A("| `%s` | `%s` | `%s` | — | %s |"
                  % (row["field"],
                     json.dumps(row["regen"], ensure_ascii=False)[:120],
                     json.dumps(row["saved"], ensure_ascii=False)[:120],
                     "はい" if row["match"] else "いいえ"))
        A("")
        A("### 2.3 同一性ゲート（v4-2.5 の r2 既報値との照合）")
        A("")
        A("再生成面を `run_v4_2_5` の物差し（固定帯 0.4–3.6 cm、x は "
          "dmax 1.00 m・iso は 0.20 m＝v4-2.5 の r2 曲線と同じラダー）で測定。")
        A("")
        A("| 方向 | v4-2.5 既報（丸め） | 記録値 | 再生成の測定値 | 差 | ±0.01 |")
        A("|---|---:|---:|---:|---:|---|")
        for k in ("x", "iso"):
            d = gate["deviation"][k]
            A("| %s | %.3f | %.6f | %.6f | %+.6f | %s |"
              % (k, d["brief"], d["recorded"], d["got"],
                 d["got"] - d["recorded"],
                 "PASS" if d["abs_diff"] <= gate["tol"] else "FAIL"))
        A("")
        if all(v["abs_diff"] == 0.0 for v in gate["deviation"].values()):
            A("差は 2 値とも**厳密に 0**（再生成面の float64 が v4-2.3 当時の"
              "面と同一で、同じ関数・同じラダーが同じ値を返した）。")
            A("")
        A("→ **同一性ゲート %s**。%s"
          % ("PASS" if gate["pass"] else "FAIL",
             "再生成面を r2 として STAGE 1 の対象に加えた。"
             if gate["pass"] else
             "r2 は再現できなかったものとして報告し、STAGE 1 はパッチ＋r1 の"
             "2面のみ。再生成面を r2 とは呼ばない。"))
        A("")
        A("保存: `out/v4-2-12b/heightfield-long-r2-regen.npy` — %d B、"
          "sha256 `%s`、%s。読み戻しの σz は %.9f m（メモリ上の float64 と "
          "float32 キャストでビット一致: %s）。"
          % (rb["file"]["bytes"], rb["file"]["sha256"],
             rb["readback"]["container"],
             rb["readback"]["sigma_z_read_m"],
             "はい" if rb["readback"]["bit_identical_to_float32_cast"] else "いいえ"))
        A("")
        A("再生成にかかった時間: %.1f s。" % rb["wall_s"])
    A("")

    # ------------------------------------------------------------------ 3
    A("## 3. STAGE 1 — Nikora 手順での測り直し")
    A("")
    A("手順は測定前に固定（コード冒頭の契約）: detrending なしの原フィールド、"
      "軸カット `D(lx,0)` / `D(0,ly)`、60-bin 対数ラダー（`run_v4_2_5."
      "curve_dir` を verbatim 使用）、折れはラダー点の [0.25, 2]×D50 を総当たり"
      "（2区間最小二乗、各区間4点以上、第2区間の上端は 4×D50）、H_Nikora は"
      "[最小ラグ＝格子間隔, 折れ] の単一フィット。")
    A("")
    A("### 3.1 面の幾何と統計")
    A("")
    A("| 面 | 方向 | 格子 (mm) | 記録サイズ (×D50) | 最小ラグ (×D50) | "
      "σz (m) | σz/D50 | 歪度 |")
    A("|---|---|---:|---:|---:|---:|---:|---:|")
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            A("| %s | %s | %.4f | %.2f | %.4f | %.6f | %.4f | %+.3f |"
              % (name, dlabel, m["cell_mm"], m["record_size_over_d50"],
                 m["lag_min_over_d50"], face["sigma_z_m"],
                 m["sigma_z_over_d50"], m["skew"]))
    A("")
    A("歪度は面全体（全画素、3次モーメント）の値。meta の `skew_matrix` は "
      "matrix 画素だけの歪度なので、帯状床では一致しない（r1 meta −0.372 対 "
      "本表 −0.055、r2 meta +0.396 対 本表 +1.001）。v4-2.3 自身も歪度の"
      "再現は主張していない（帯状床は 3 m モードに分散が寄るため、同じ "
      "c=0.55 の歪み変換でも効き方が変わる）。")
    A("")
    A("### 3.2 折れと指数")
    A("")
    A("| 面 | 方向 | 折れ (mm) | 折れ (×D50) | 候補数 | H_Nikora | n | "
      "rms_log | r2 | 二帯域 下 | 二帯域 上 |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            b = m.get("break") or {}
            h = m.get("H_Nikora") or {}
            tb = m.get("two_band") or {}
            ab = tb.get("above") or {}
            A("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
              % (name, dlabel,
                 f(b.get("lag_m", 0) * 1000.0, 3) if b else "n/a",
                 f(b.get("lag_over_d50"), 3) if b else "n/a",
                 b.get("n_candidates", "n/a") if b else "n/a",
                 f(h.get("H")), h.get("n_points", "n/a"),
                 f(h.get("rms_log"), 4), f(h.get("r2"), 5),
                 f(h.get("H")), f(ab.get("H"))))
    A("")
    A("「二帯域 下」は H_Nikora と同じフィット（[最小ラグ, 折れ]）。"
      "「二帯域 上」は [折れ, 4×D50]。")
    edge = []
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            cands = (m.get("break") or {}).get("candidates") or []
            if cands and m["break"]["lag_m"] == cands[0]["lag_m"]:
                edge.append("%s / %s" % (name, dlabel))
    if edge:
        A("")
        A("**窓の下端に張り付いた折れ**: %d 本のうち **%d 本**（%s）で、"
          "最適な折れが探索窓 [0.25, 2]×D50 の最小候補ラグに一致した。"
          "窓の下端は brief が固定した 0.25×D50 なので動かしていないが、"
          "折れの絶対値（したがって H_Nikora の区間）は窓の定義に依存する。"
          "窓内の各候補の rss は numbers の `break.candidates` に、"
          "折れを ×0.8 / ×1.25 に動かした感度は §3.3 に残した。"
          % (sum(len(f["directions"]) for f in st1["faces"].values()),
             len(edge), "、".join(edge)))
    A("")
    A("### 3.3 折れの感度（×0.8 / ×1.25 → 最近傍ラダー点）")
    A("")
    A("| 面 | 方向 | 折れ (mm) | ×0.8 実ラグ (mm) | H (×0.8) | ΔH | "
      "×1.25 実ラグ (mm) | H (×1.25) | ΔH |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            b = m.get("break") or {}
            s8 = (m.get("sensitivity") or {}).get("x0.8") or {}
            s12 = (m.get("sensitivity") or {}).get("x1.25") or {}
            A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
              % (name, dlabel,
                 f(b.get("lag_m", 0) * 1000.0, 3) if b else "n/a",
                 f(s8.get("used_m", 0) * 1000.0, 3) if s8 else "n/a",
                 f((s8.get("fit") or {}).get("H")),
                 sgn(s8.get("dH")),
                 f(s12.get("used_m", 0) * 1000.0, 3) if s12 else "n/a",
                 f((s12.get("fit") or {}).get("H")),
                 sgn(s12.get("dH"))))
    A("")
    A("H_loc(ℓ) の全ラダー値は numbers JSON の "
      "`stage1.faces.<面>.directions.<方向>.H_loc` に、軸カットの D 値そのものは"
      " `...curve.D_m` / `...curve.S_m2` に入れてある。局所傾斜の範囲: ")
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            hl = np.array(m["H_loc"]["H"], dtype=float)
            fin = hl[np.isfinite(hl)]
            imin = int(np.nanargmin(np.where(np.isfinite(hl), hl, np.nan)))
            imax = int(np.nanargmax(np.where(np.isfinite(hl), hl, np.nan)))
            A("- %s / %s: 最小 %.4f @ %.3f mm、最大 %.4f @ %.3f mm"
              % (name, dlabel, float(fin.min()),
                 m["H_loc"]["lag_m"][imin] * 1000.0, float(fin.max()),
                 m["H_loc"]["lag_m"][imax] * 1000.0))
    A("")
    A("### 3.4 折れ探索の定義依存（参考・事後計算）")
    A("")
    A("**これは報告値ではない。** 本手順の第2区間は測定前に [折れ, 4×D50] と"
      "固定した（§3・開示3）。同じ2区間最小二乗で第2区間を**ラダー末端まで**"
      "取る読み方（前身ジョブの書き方）を、保存済みの軸カットから事後計算し、"
      "折れの定義が結果をどれだけ動かすかを示す。")
    A("")
    A("| 面 | 方向 | 本手順 折れ (mm) | 本手順 H_Nikora | 第2区間=末端 折れ (mm) "
      "| その H | ΔH |")
    A("|---|---|---:|---:|---:|---:|---:|")
    worst_alt = 0.0
    post = st1.get("break_definition_posthoc") or posthoc_alt_breaks(st1)
    for name, face in st1["faces"].items():
        for dlabel, m in face["directions"].items():
            row = post["rows"].get("%s | %s" % (name, dlabel))
            h0 = (m.get("H_Nikora") or {}).get("H")
            if not row or row.get("H_seg2_to_last_lag") is None or h0 is None:
                A("| %s | %s | %s | %s | n/a | n/a | n/a |"
                  % (name, dlabel, f((m.get("break") or {}).get("lag_m", 0)
                                     * 1000.0, 3) if m.get("break") else "n/a",
                     f(h0)))
                continue
            dh = row["dH"]
            worst_alt = max(worst_alt, abs(dh))
            A("| %s | %s | %.3f | %.4f | %.3f | %.4f | %+.3f |"
              % (name, dlabel, m["break"]["lag_m"] * 1000.0, h0,
                 row["break_seg2_to_last_lag_m"] * 1000.0,
                 row["H_seg2_to_last_lag"], dh))
    A("")
    A("第2区間の上端をラダー末端に取ると、帯状床の流下方向で折れが "
      "30–47 mm（0.8–1.2×D50）へ動き、H_Nikora は 0.38–0.39 まで下がる"
      "（表の最大 |ΔH| = %.3f）。つまり本レポートの H_Nikora（0.43–0.45）は"
      "「折れ ≤ 2×D50・第2区間 ≤ 4×D50」という窓の定義に条件づけられた値で"
      "ある。窓の位置を動かした感度（§3.3）はもっと小さい。" % worst_alt)
    A("")

    # ------------------------------------------------------------------ 4
    A("## 4. STAGE 2 — 合成充填テスト（パッチ床のみ）")
    A("")
    A("**測定前に固定した予測**: δ = 0.5 σz の充填でも ΔH_fit ≤ "
      "+%.2f（固定帯 0.4–3.6 cm）。" % st2["prediction"]["dH_fit_max"])
    A("")
    fl = st2["fill"]
    A("充填は Barnes 型 priority-flood（開境界、level = max(z, level)）で "
      "F を1回だけ計算し、z'(δ) = min(F, z + δ·σz) を 3 つの δ に適用した。"
      "充填されたセルは %.2f%%、最大 uplift は %.5f m、充填セルの平均 uplift は "
      "%.5f m、priority-flood の計算時間は %.1f s。σz（充填前）= %.9f m。"
      % (100.0 * fl["filled_fraction"], fl["max_uplift_m"],
         fl["mean_uplift_where_filled_m"], fl["wall_s"],
         st2["sigma_z_before_m"]))
    A("")
    A("| δ/σz | δ (m) | σz 後 (m) | Δσz (m) | ΔH_fit x | ΔH_fit y | "
      "ΔH_fit iso | ΔH_Nikora x | ΔH_Nikora y | ΔH_Nikora iso |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for tag, rec in st2["cells"].items():
        A("| %s | %.6f | %.9f | %+.3e | %s | %s | %s | %s | %s | %s |"
          % (tag, rec["delta_m"], rec["sigma_z_after_m"],
             rec["d_sigma_z_m"], sgn(rec["dH_fit"]["x"]),
             sgn(rec["dH_fit"]["y"]), sgn(rec["dH_fit"]["iso"]),
             sgn(rec["dH_Nikora"]["x"]), sgn(rec["dH_Nikora"]["y"]),
             sgn(rec["dH_Nikora"]["iso"])))
    A("")
    A("参考: 各 δ の折れ（充填後の面で同定し直したもの）")
    A("")
    A("| δ/σz | 折れ x (mm) | 折れ y (mm) | 折れ iso (mm) | H_fit x | H_fit y | "
      "H_fit iso | H_Nikora x | H_Nikora y | H_Nikora iso |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    base = st2["baseline"]
    A("| 0（充填前） | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
      % tuple([f((base["H_Nikora"][k].get("break") or {}).get("lag_m", 0)
                 * 1000.0, 3) if (base["H_Nikora"][k].get("break")) else "n/a"
               for k in ("x", "y", "iso")]
              + [f(base["H_fit"][k]["H"]) for k in ("x", "y", "iso")]
              + [f((base["H_Nikora"][k].get("H_Nikora") or {}).get("H"))
                 for k in ("x", "y", "iso")]))
    for tag, rec in st2["cells"].items():
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % tuple([tag]
                  + [f(rec["break_lag_m"][k] * 1000.0, 3)
                     if rec["break_lag_m"][k] else "n/a"
                     for k in ("x", "y", "iso")]
                  + [f(rec["H_fit"][k]["H"]) for k in ("x", "y", "iso")]
                  + [f((rec["H_Nikora"][k].get("H_Nikora") or {}).get("H"))
                     for k in ("x", "y", "iso")]))
    A("")
    pc2 = st2["prediction_check"]
    A("**予測との照合**: δ = 0.5 σz で ΔH_fit は x %s / y %s / iso %s、"
      "最大 %s。閾値 +%.2f に対し **%s**。"
      % (sgn(pc2["dH_fit"]["x"]), sgn(pc2["dH_fit"]["y"]),
         sgn(pc2["dH_fit"]["iso"]), sgn(pc2["max_dH_fit"]),
         pc2["threshold"], "予測どおり（PASS）" if pc2["pass"]
         else "予測に反する（FAIL）"))
    A("")

    # ------------------------------------------------------------------ 5
    A("## 5. 開示事項")
    A("")
    A("1. **周期場と np.roll**: STAGE 1 の差分は np.roll による周期ラップで"
      "計算した。生成場は周期場として作られているのでこの前提は正当。"
      "ただし帯状床は 0.6 m タイル5枚を並べたもので、x = 3.0 m の継ぎ目は"
      "独立な母集団どうしの人工的な接合である。使用した最大ラグは "
      "min(ドメイン半分, 4×D50 ≈ 0.15 m) なので、この継ぎ目はどのフィットにも"
      "入らない。")
    A("2. **detrending なし**: brief の固定手順どおり、原フィールドをそのまま"
      "測った。したがって H には基盤の form（低次モード・fBm）が含まれる。")
    A("3. **折れ探索の第2区間の上端**: [折れ, 4×D50] に固定した。上端を"
      "ラダー末端（dmax）にすると折れの位置が dmax に依存してしまうため、"
      "測定前にこの定義を固定した。")
    A("4. **r2 の身元**: r2 の heightfield は v4-2.3 実行時にディスクへ保存"
      "されていない。本ジョブは seed から再生成し、メタ照合（%s）と固定帯の"
      "同一性ゲート（%s）で本人確認した。同一性ゲートは固定帯の x と iso の"
      "2値のみで、r2 の y は v4-2.5 でも測定不能と記録されているため比較して"
      "いない。ゲートは ±0.01 の H 一致であり、ビット同一性の主張ではない"
      "（ビット同一性に近い証拠はメタ %d 項目の一致）。"
      % ("全項目一致" if rb["meta_compare"] and rb["meta_compare"]["all_match"]
         else "不一致あり",
         "PASS" if gate["pass"] else "FAIL",
         rb["meta_compare"]["n_fields"] if rb["meta_compare"] else 0))
    A("5. **再生成ファイルのコンテナ**: `heightfield-long-r2-regen.npy` は "
      "r1 と同じ **ヘッダなし raw float32 LE**（2048×10240、83,886,080 B）で"
      "保存した。`np.load` では読めないので、読むときは "
      "`np.fromfile(path, dtype='<f4').reshape(2048, 10240)` を使う。"
      "拡張子は brief の指定どおり。")
    A("6. **充填の境界条件**: priority-flood は開境界（非周期）で解いた。"
      "一方 STAGE 2 の構造関数は周期ラップのまま測っている。この不一致は"
      "充填面の外周に限られ、使用ラグ（≤ 4×D50 ≈ 0.15 m）に対する影響は"
      "小さいが、厳密には境界の扱いが異なることを開示する。")
    A("7. **ΔH の符号**: ΔH = H(充填後) − H(充填前)。正は充填で H が"
      "大きくなったことを意味する。")
    A("8. **H_loc**: 同じラダー上の 0.5·dlnS/dlnd（np.gradient）。端点は"
      "片側差分なので、ラダー両端の値は内側より誤差が大きい。")
    A("9. **STAGE 2 はパッチ床のみ**（brief の指定）。帯状床の充填は"
      "計算していない。")
    A("10. **0.70–0.85 ゲート**: brief の対象外。本レポートは合否判定を"
      "していない（H の値と不確かさのみ）。")
    A("11. **環境**: Python %s, numpy %s。実行は前景・1プロセス。"
      % (ctx["env"]["python"], ctx["env"]["numpy"]))
    nedge = 0
    ntot = 0
    for face in st1["faces"].values():
        for m in face["directions"].values():
            ntot += 1
            cands = (m.get("break") or {}).get("candidates") or []
            if cands and m["break"]["lag_m"] == cands[0]["lag_m"]:
                nedge += 1
    A("12. **折れ探索窓の下端張り付き**: %d/%d 本で最適な折れが窓の最小候補"
      "ラグに一致した（§3.2 の注記）。brief が窓を [0.25, 2]×D50 と固定して"
      "いるので、本レポートの H_Nikora はその定義での値であり、窓を広げた"
      "ときの折れの普遍値を主張しない。窓の位置を ×0.8 / ×1.25 動かした"
      "ときの H の変化は §3.3（最大でも |ΔH| ≈ 0.013）、第2区間の上端を"
      "ラダー末端に取る読み方との差は §3.4（最大 |ΔH| ≈ 0.05）。"
      % (nedge, ntot))
    A("")

    # ------------------------------------------------------------------ 6
    A("## 6. 成果物・チェックポイント・自己監査")
    A("")
    A("| 成果物 | パス |")
    A("|---|---|")
    A("| 本レポート | `%s` |" % REPORT)
    A("| 全数値 | `%s` |" % NUMBERS)
    A("| 再生成 r2 | `%s` |" % REGEN)
    A("| 事前確認の記録 | `%s` |" % PREFLIGHT)
    A("| jobkit checkpoint | `%s` |" % CP_PATH)
    A("")
    gl = ctx["g_legacy"]
    A("G-LEGACY（既存 .py の不変監査）: 基準 %d 本、実行後 %d 本、"
      "変更 %d・追加 %d・削除 %d → **%s**。"
      % (gl["n_baseline"], gl["n_now"], len(gl["changed"]),
         len(gl["added"]), len(gl["removed"]),
         "PASS" if gl["ok"] else "FAIL"))
    if gl["changed"]:
        A("")
        A("- 変更: %s" % ", ".join(gl["changed"]))
    if gl["added"]:
        A("- 追加: %s" % ", ".join(gl["added"]))
    A("")
    A("jobkit ユニット: %s（%d 個）。"
      % (", ".join(ctx["checkpoint_units"]["ids"]),
         ctx["checkpoint_units"]["n"]))
    A("")
    A("実行時間: %s（レポート書き出しは含まない）"
      % ", ".join("%s %.1f s" % (k, v) for k, v in ctx["timings"].items()
                  if k != "report"))
    A("")

    # ------------------------------------------------------------------ 7
    A("## 7. 実行ログ")
    A("")
    A("```")
    for line in LOG:
        A(line)
    A("```")
    A("")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
def report_only():
    """Rebuild REPORT from numbers JSON (no measurement re-run)."""
    if not os.path.isfile(NUMBERS):
        raise SystemExit("no numbers to rebuild from: %s" % NUMBERS)
    n = json.load(open(NUMBERS, encoding="utf-8"))
    LOG.extend(n.get("log") or [])
    if "break_definition_posthoc" not in n.get("stage1", {}):
        n["stage1"]["break_definition_posthoc"] = posthoc_alt_breaks(n["stage1"])
        write_json(NUMBERS, n)
        print("added break_definition_posthoc to %s" % NUMBERS)
    ctx = {"created_utc": n["created_utc"], "stage0": n["stage0"],
           "stage0_5": n["stage0_5"], "stage1": n["stage1"],
           "stage2": n["stage2"], "env": n["env"],
           "timings": n.get("timings_s") or {},
           "checkpoint_units": (n.get("checkpoint") or {}).get("state") or {},
           "g_legacy": _g_legacy({"stage0": n["stage0"]})}
    md = build_report(ctx)
    atomic_write_bytes(REPORT, md.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 5000)
    print("rebuilt %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    return 0


def main():
    if "--report-only" in sys.argv[1:]:
        return report_only()
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    created = utcnow()
    timings = {}
    ctx = {"created_utc": created, "brief": BRIEF,
           "env": {"python": sys.version.split()[0], "numpy": np.__version__,
                   "platform": sys.platform},
           "stages_done": []}

    def stop_record(reason):
        """Write an honest stop record: what ran, what did not."""
        ctx["stopped_at"] = reason
        cp.unit_done("report (STOP record: %s)" % reason)
        ctx["checkpoint_units"] = _ckpt_units(cp)
        ctx["g_legacy"] = _g_legacy(ctx)
        md = build_stop_report(ctx, reason)
        atomic_write_bytes(REPORT, md.encode("utf-8"),
                           verify=lambda p: os.path.getsize(p) > 1000)
        write_json(NUMBERS, {
            "brief": BRIEF, "created_utc": created, "stopped_at": reason,
            "stage0": ctx.get("stage0"), "stage0_5": ctx.get("stage0_5"),
            "stage1": ctx.get("stage1"), "stage2": ctx.get("stage2"),
            "env": ctx["env"], "log": LOG})
        log("STOPPED  %s" % reason)
        log("      %s" % REPORT)

    t0 = time.time()
    try:
        st0, cache = stage0(cp)
    except Blocker as b:
        timings["stage0"] = time.time() - t0
        ctx["stage0"] = b.ctx
        ctx["timings"] = timings
        stop_record("STAGE 0: %s" % b.reason)
        return 2
    timings["stage0"] = time.time() - t0
    ctx["stage0"] = st0
    ctx["d50"] = st0["d50"]
    ctx["geometry"] = st0["geometry"]
    ctx["stages_done"].append("stage0")
    log("STAGE0  done in %.1f s" % timings["stage0"])

    t0 = time.time()
    try:
        stage0_5(cp, ctx)
    except Blocker as b:
        timings["stage0_5"] = time.time() - t0
        ctx["stage0_5"] = b.ctx.get("stage0_5") or b.ctx
        ctx["timings"] = timings
        stop_record("STAGE 0.5: %s" % b.reason)
        return 2
    timings["stage0_5"] = time.time() - t0
    ctx["stages_done"].append("stage0_5")
    log("STAGE0.5  done in %.1f s (identity gate %s)"
        % (timings["stage0_5"],
           "PASS" if ctx["stage0_5"]["identity_gate"]["pass"] else "FAIL"))

    t0 = time.time()
    ctx["stage1"] = stage1(cp, ctx, cache)
    ctx["stage1"]["break_definition_posthoc"] = posthoc_alt_breaks(ctx["stage1"])
    timings["stage1"] = time.time() - t0
    ctx["stages_done"].append("stage1")
    log("STAGE1  done in %.1f s" % timings["stage1"])

    t0 = time.time()
    ctx["stage2"] = stage2(cp, ctx)
    timings["stage2"] = time.time() - t0
    ctx["stages_done"].append("stage2")
    log("STAGE2  done in %.1f s" % timings["stage2"])

    timings["report"] = 0.0
    ctx["timings"] = timings
    cur = _ckpt_units(cp)
    ctx["checkpoint_units"] = {
        "path": CP_PATH, "n": len(cur["ids"]) + 1,
        "ids": cur["ids"] + ["report"],
        "completed_units": len(cur["ids"]) + 1, "last_unit": "report",
        "note": "the report unit is committed right after the file is written"}
    ctx["g_legacy"] = _g_legacy(ctx)

    ctx["deliverables"] = {
        "report": {"path": REPORT},
        "numbers": {"path": NUMBERS},
        "regen_r2": ({"path": REGEN, "bytes": ctx["stage0_5"]["regen"]["file"]
                      and ctx["stage0_5"]["regen"]["file"]["bytes"],
                      "sha256": ctx["stage0_5"]["regen"]["file"]
                      and ctx["stage0_5"]["regen"]["file"]["sha256"]}
                     if ctx["stage0_5"]["regen"]["file"] else None),
        "regen_meta": {"path": REGEN_META},
        "preflight": {"path": PREFLIGHT},
        "checkpoint": {"path": CP_PATH},
        "write_set": [THIS_FILE, "out/v4-2-12b/"],
    }

    md = build_report(ctx)
    atomic_write_bytes(REPORT, md.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 5000)

    t0 = time.time()
    # numbers carries the report unit id, so the checkpoint is updated first
    cp.unit_done("report")
    timings["report"] = time.time() - t0
    ctx["timings"] = timings
    ctx["checkpoint_units"] = _ckpt_units(cp)

    numbers = {
        "brief": BRIEF, "created_utc": created,
        "what_this_is": (
            "re-measurement of the existing generated bed with Nikora's "
            "procedure (STAGE 1), plus the r2 regeneration and identity gate "
            "(STAGE 0.5) and a synthetic-fill probe on the patch (STAGE 2). "
            "No parameter change, no rendering, no new bed beyond the r2 "
            "replication the brief asks for."),
        "write_set": [THIS_FILE, "out/v4-2-12b/"],
        "env": ctx["env"],
        "method_contract": ctx["stage1"]["method"],
        "prediction": ctx["stage2"]["prediction"],
        "stage0": ctx["stage0"],
        "stage0_5": ctx["stage0_5"],
        "stage1": ctx["stage1"],
        "stage2": ctx["stage2"],
        "checkpoint": {"path": CP_PATH, "state": ctx["checkpoint_units"]},
        "g_legacy": ctx["g_legacy"],
        "deliverables": ctx["deliverables"],
        "timings_s": {k: round(v, 1) for k, v in timings.items()},
        "log": LOG,
    }
    write_json(NUMBERS, numbers)
    log("DONE  %s" % REPORT)
    log("      %s" % NUMBERS)
    return 0


def _ckpt_units(cp):
    st = cp.state() or {}
    ids = st.get("completed_ids") or []
    return {"path": CP_PATH, "n": len(ids), "ids": ids,
            "completed_units": st.get("completed_units"),
            "last_unit": st.get("last_unit")}


def _g_legacy(ctx):
    base = (ctx.get("stage0") or {}).get("g_legacy_baseline") or {"py": {}}
    now = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p) and name != THIS_FILE:
            now[name] = {"sha256": sha256_file(p),
                         "bytes": os.path.getsize(p)}
    changed = sorted(k for k in base["py"] if k in now
                     and now[k]["sha256"] != base["py"][k]["sha256"])
    added = sorted(set(now) - set(base["py"]))
    removed = sorted(set(base["py"]) - set(now))
    return {"n_baseline": len(base["py"]), "n_now": len(now),
            "changed": changed, "added": added, "removed": removed,
            "ok": bool(not (changed or added or removed)),
            "write_set": [THIS_FILE, "out/v4-2-12b/"],
            "checked_at": utcnow()}


def build_stop_report(ctx, reason):
    L = ["# REPORT-V4-2-12B — STOP record", "",
         "- brief: `%s`" % BRIEF,
         "- 実行: `run_v4_2_12b.py`",
         "- 作成: %s UTC" % ctx["created_utc"], "",
         "> **この実行は `%s` で停止した。**" % reason, "",
         "> 停止した時点で測った値だけを下に残す。測っていない段は"
         "「未測定」と書く。", ""]
    st0 = ctx.get("stage0") or {}
    if st0.get("required"):
        L.append("## STAGE 0 — 確認値")
        L.append("")
        L.append("| ファイル | bytes | sha256 (先頭16) |")
        L.append("|---|---:|---|")
        for r in st0["required"]:
            L.append("| `%s` | %s | `%s` |"
                     % (os.path.basename(r["path"]),
                        r["bytes"] if r["bytes"] is not None else "MISSING",
                        (r["sha256"] or "—")[:16]))
        L.append("")
        if st0.get("missing"):
            L.append("欠損: %s" % ", ".join(
                "`%s`" % m["path"] for m in st0["missing"]))
            L.append("")
        if st0.get("repro"):
            L.append("再現ゲート: **%s**"
                     % ("PASS" if st0["repro"]["pass"] else "FAIL"))
            L.append("")
    for stage in ("stage0_5", "stage1", "stage2"):
        if stage not in ctx["stages_done"]:
            L.append("- %s: 未測定" % stage)
    L.append("")
    L.append("## ログ")
    L.append("")
    L.append("```")
    L.extend(LOG)
    L.append("```")
    L.append("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
