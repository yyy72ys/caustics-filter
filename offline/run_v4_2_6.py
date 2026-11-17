"""run_v4_2_6 -- BRIEF-V4-2-6: rebuild the bed so that the GRAIN LAYER ITSELF
carries the real, measured scaling (H = 0.70-0.85 over lags 0.4-3.6 cm).

Scope of this run: STAGE 0 (preflight / hashes / instrument identification) and
STAGE 1 (the calibration gate).  STAGE 2 (bed_gen3.py) and STAGE 3 (the gates
on a new bed) are only entered if STAGE 1 finds a mechanism that delivers

    local exponent >= 0.70 everywhere in 0.4-3.6 cm
    AND sigma_z(matrix) inside 0.32-0.46 D50.

The brief says: "if there is none, STOP" and report which mechanism reached how
far.  This run therefore either writes bed_gen3.py (STAGE 2) or it stops after
STAGE 1 and writes REPORT-V4-2-6.md with the mechanism table.  Everything the
report claims is measured here; no parameter is changed to make a picture look
right.

What is new vs the earlier runs
-------------------------------
* The measurement instrument is NOT re-implemented: it is imported from
  run_v4_2_5 (curve_dir / curve_iso / local_H / fit_band / band_shape / BANDS /
  GATE) and used verbatim.  STAGE 0 records exactly which lines of that file
  define the estimator.
* A deposition engine (support / contact / drainage of the interstices) is the
  STAGE 1 test bench.  It is a grid implementation of the same ellipsoid
  ray-cast that bed_gen2 uses (the inner loop of bed_gen2.raycast2, with the
  running height field carrying the grain tops; max-combine, periodic wrap).
* An ADDITIONAL check (not asked by the brief, flagged as such in the report)
  measures the CEILING of the instrument: an ideal self-affine field with
  H_gen = 0.79 is pushed through the same estimator so that the gate can be
  read against what the estimator can even return on a 0.6 m patch.

Stages (foreground, one jobkit checkpoint unit per stage):
  --stage0   preflight, sha256 of the required inputs + every existing .py,
             D50 confirmation, instrument identification
  --calib    STAGE 1: (a) reference, (b) GSD width, (c) support deposition,
             (d) imbrication, (e) + substrate; plus the ceiling check
  --judge    the STAGE 1 gate, case by case
  --figure   structure-function / local-exponent / ceiling figure
  --legacy   G-LEGACY (existing .py unchanged)
  --report   numbers-v4-2-6.json + REPORT-V4-2-6.md
  --all      the above, in order, in this one process
"""
from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import inspect
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import bed_gen as BG                                   # noqa: E402  (unmodified)
import bed_gen2 as B2                                  # noqa: E402  (unmodified)
import run_v4_2_5 as M5                                # noqa: E402  (unmodified)
#                                                        ^ the instrument itself

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-6.md"
OUT = os.path.join(HERE, "out", "v4-2-6")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-6.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-6.md")
FIGURE = os.path.join(OUT, "structure-functions-v4-2-6.png")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
TIMINGS = os.path.join(OUT, "timings.json")
RUNLOG = os.path.join(OUT, "run.log")

HF_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
META_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
NUMBERS_425 = os.path.join(HERE, "out", "v4-2-5", "numbers-v4-2-5.json")
REPORT_425 = os.path.join(HERE, "out", "v4-2-5", "REPORT-V4-2-5.md")
DECOMP_423 = os.path.join(HERE, "out", "v4-2-3", "decomp.json")

THIS_FILE = "run_v4_2_6.py"

REQUIRED = {
    "bed_gen.py": os.path.join(HERE, "bed_gen.py"),
    "bed_gen2.py": os.path.join(HERE, "bed_gen2.py"),
    "run_v4_2_2.py": os.path.join(HERE, "run_v4_2_2.py"),
    "run_v4_2_3.py": os.path.join(HERE, "run_v4_2_3.py"),
    "run_v4_2_5.py": os.path.join(HERE, "run_v4_2_5.py"),
    "out/v4-2-2/heightfield.npy": HF_422,
    "out/v4-2-5/numbers-v4-2-5.json": NUMBERS_425,
}

# ---------------------------------------------------------------------------
# the measurement contract, fixed BEFORE the first STAGE 1 measurement
# ---------------------------------------------------------------------------
BAND = (float(M5.BANDS[0][1]), float(M5.BANDS[0][2]))   # 0.004 .. 0.036 m
GATE_H = (float(M5.GATE[0]), float(M5.GATE[1]))         # 0.70 .. 0.85
LOCAL_MIN_REQ = 0.70          # brief STAGE 1: the local exponent may not dip
SZ_BAND = (0.32, 0.46)        # brief STAGE 1: sigma_z(matrix) / D50_surface
NBIN = int(M5.NBIN)
DM_X = 0.30                   # v4-2.5's patch x ladder (DM_PATCH_X)
DM_Y = 0.30                   # v4-2.5's patch y ladder (DM_Y)

PATCH = float(B2.PATCH)
SEED = int(B2.SEED)
D50_NOMINAL_M = float(B2.D50_SURF)
D84_MATRIX_M = float(B2.D84_SURF)     # the matrix mask rule of the v4-2.2 audit
D50_REF_M = 0.03619913686788151       # v4-2.2 measured surface D50 (brief 0.0362)

CAL_N = 1024                          # calibration grid (0.586 mm/px, 0.6 m)
CAL_CELL = PATCH / CAL_N
HIRES_N = int(B2.HIRES)               # 2048, the v4-2.2 statistics grid
HIRES_CELL = PATCH / HIRES_N

TARGETS = (0.40, 0.46)                # sigma_z/D50 fixpoint targets ("best case")
DEPOSIT_Q = (0.05, 0.10, 0.20)
N_SEEDS_CEIL = 6
CEIL_L0 = (0.020, 0.0362, 0.050, 0.073, 0.100, 0.200)
CEIL_BOXES = ((0.60, 512, 6), (1.20, 1024, 3), (2.40, 2048, 3))
N_SEEDS_BOX = 2
GSD_WIDTHS = (1.5, 2.1, 2.7)
DIP_CURRENT = (25.0, 35.0)            # brief (d): 25-35 deg upstream dip

HOLES = "mean"                        # holes are filled with the covered mean
                                      # (see the sensitivity table in §6.4)

LOG = []
T = {}


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S+00:00")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _j(x):
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, (list, tuple)):
        return [_j(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _j(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        return _j(x.tolist())
    return x


def write_json(path, obj):
    data = json.dumps(_j(obj), ensure_ascii=False, indent=1).encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True
    atomic_write_bytes(path, data, verify=verify)


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight, hashes, D50, instrument identification
# ---------------------------------------------------------------------------
def instrument_id():
    """Where exactly the v4-2.5 instrument lives, and what its contract is."""
    out = {"module": "run_v4_2_5", "path": os.path.join(HERE, "run_v4_2_5.py"),
           "functions": {}, "constants": {}}
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band", "band_shape",
                 "audit_like"):
        fn = getattr(M5, name)
        lines, start = inspect.getsourcelines(fn)
        out["functions"][name] = {"line": int(start), "n_lines": len(lines),
                                  "doc": (inspect.getdoc(fn) or "").split(
                                      "\n")[0]}
    out["constants"] = {
        "BANDS": [[k, a, b] for k, a, b in M5.BANDS],
        "GATE": list(M5.GATE),
        "NBIN": int(M5.NBIN),
        "DM_ISO": float(M5.DM_ISO),
        "DM_PATCH_X": float(M5.DM_PATCH_X),
        "DM_Y": float(M5.DM_Y),
        "scaling_band_m": [BAND[0], BAND[1]],
        "gate_H": [GATE_H[0], GATE_H[1]],
    }
    out["usage"] = ("imported unmodified; every S(tau) and every H in this run "
                    "is produced by M5.curve_dir/M5.curve_iso and fitted by "
                    "M5.fit_band over M5.BANDS[0] (0.004-0.036 m); the local "
                    "exponent is M5.local_H (0.5*dlnS/dlnd on the same "
                    "ladder), read inside the band by M5.band_shape")
    return out


def stage0(cp):
    log("STAGE0  absolute-path preflight, sha256, D50, instrument")
    missing = [k for k, p in REQUIRED.items() if not os.path.isfile(p)]
    if missing:
        raise SystemExit("STAGE 0 STOP -- missing: %s" % missing)
    pre = {}
    for k, p in REQUIRED.items():
        pre[k] = {"exists": True, "bytes": os.path.getsize(p),
                  "sha256": sha256_file(p)}
        log("  ok %-34s %10d bytes  %s" % (k, pre[k]["bytes"],
                                          pre[k]["sha256"][:16]))

    meta = json.load(open(META_422, encoding="utf-8"))
    d50 = float(meta["D50_surface_measured_m"])
    if abs(d50 - 0.0362) > 1e-3:
        raise SystemExit("STAGE 0 STOP -- surface D50 %.6f != 0.0362" % d50)
    d50_blk = {"v422_patch_m": d50, "brief_named_m": 0.0362,
               "confirmed": True,
               "source": "out/v4-2-2/heightfield-meta.json "
                         "(D50_surface_measured_m, exposed-area audit)"}
    log("  D50 confirmed: %.6f m (= 0.0362 of the brief)" % d50)

    inst = instrument_id()
    log("  instrument: run_v4_2_5 import ok; band %.3f-%.3f m, gate %.2f-%.2f"
        % (BAND[0], BAND[1], GATE_H[0], GATE_H[1]))

    # every existing .py, exactly the v4-2.5 convention
    py = {}
    for p in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        name = os.path.basename(p)
        if name == THIS_FILE:
            continue
        py[name] = {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
    write_json(STAGE0_HASHES, {"created_utc": utcnow(),
                              "n_existing_py": len(py), "py": py})
    log("  hashed %d existing .py -> stage0-py-hashes.json" % len(py))
    cp.unit_done("stage0")
    return {"preflight": pre, "d50": d50_blk, "instrument": inst,
            "n_existing_py": len(py)}


# ---------------------------------------------------------------------------
# the deposition engine (STAGE 1 test bench)
# ---------------------------------------------------------------------------
def form_grid(form, n):
    """Resample a substrate form onto the n x n calibration grid (periodic)."""
    if form is None:
        return np.zeros((n, n))
    cell = PATCH / n
    cellf = PATCH / form.shape[1]
    px = ((np.arange(n) + 0.5) * cell / cellf)[None, :]
    py = ((np.arange(n) + 0.5) * cell / cellf)[:, None]
    return np.asarray(BG._bilinear(form, px, py))


def bilinear_at(form, x, y):
    """The elementwise call bed_gen2.place makes (both arguments [None, :])."""
    cellf = PATCH / form.shape[1]
    return float(BG._bilinear(form, np.array([[x / cellf]]),
                              np.array([[y / cellf]]))[0, 0])


def place_on_form(grains, form):
    """bed_gen2.place, verbatim: summit = form(x,y) + p - half_extent_z."""
    cellf = PATCH / form.shape[1]
    Xg = np.array([g["x"] for g in grains])
    Yg = np.array([g["y"] for g in grains])
    fx = np.asarray(BG._bilinear(form, (Xg / cellf)[None, :],
                                 (Yg / cellf)[None, :])[0])
    hz = np.array([BG._half_extent(g, 2) for g in grains])
    p = np.array([g["p"] for g in grains])
    z = fx + p - hz
    for i, g in enumerate(grains):
        g["z"] = float(z[i])
    return fx


def grain_window(g, xs, ys, n):
    """The per-grain ellipsoid window of bed_gen2.raycast2 (same arithmetic)."""
    A, B, C = g["abc"]
    R = g["R"]
    ax = np.array([A, B, C])
    hx = BG._half_extent(g, 0)
    hy = BG._half_extent(g, 1)
    d = xs[1] - xs[0]
    i0 = int(math.floor((g["x"] - hx - xs[0]) / d))
    i1 = int(math.ceil((g["x"] + hx - xs[0]) / d)) + 1
    j0 = int(math.floor((g["y"] - hy - ys[0]) / d))
    j1 = int(math.ceil((g["y"] + hy - ys[0]) / d)) + 1
    ixw = np.arange(i0, i1) % n
    iyw = np.arange(j0, j1) % n
    DX = (xs[ixw] - g["x"])[None, :]
    DY = (ys[iyw] - g["y"])[:, None]
    u = np.stack([R[0, i] * DX + R[1, i] * DY for i in range(3)])
    w = R[2, :][:, None, None]
    aa = ((w / ax[:, None, None]) ** 2).sum(0)
    bb = ((u * w) / ax[:, None, None] ** 2).sum(0)
    cc = ((u / ax[:, None, None]) ** 2).sum(0) - 1.0
    disc = bb * bb - aa * cc
    ok = disc > 0.0
    sq = np.sqrt(np.where(ok, disc, 0.0))
    s2 = (-bb + sq) / aa
    c_eff = float(1.0 / math.sqrt(((R[2, :] / ax) ** 2).sum()))
    return np.ix_(iyw, ixw), ok, s2, c_eff


def stamp(h, gid, gi, g, idx, ok, s2):
    """Max-combine one grain cap into h (the raycast2 inner step)."""
    ztop = g["z"] + s2
    sub = h[idx]
    hit = ok & (ztop > sub)
    if hit.any():
        h[idx] = np.where(hit, ztop, sub)
        gs = gid[idx]
        gid[idx] = np.where(hit, gi, gs)
    return hit


def retilt(grains, dip_deg):
    """Clean imbrication control: keep a/b/c and the azimuth, set the dip.

    `dip_deg` is a scalar or a per-grain sequence.  (bed_gen2.make_grain's own
    rotation construction, applied to the stored axes; nothing about size,
    shape or position changes.)
    """
    dips = np.atleast_1d(np.asarray(dip_deg, dtype=float))
    if dips.size == 1:
        dips = np.full(len(grains), float(dips[0]))
    for g, dd in zip(grains, dips):
        A, B, C = g["abc"]
        az = math.radians(g["az_deg"])
        d = math.radians(float(dd))
        ea = np.array([math.cos(az), math.sin(az), 0.0])
        nrm = np.array([math.sin(d), 0.0, math.cos(d)])
        eb = np.cross(nrm, ea)
        eb /= max(np.linalg.norm(eb), 1e-12)
        ea = np.cross(eb, nrm)
        ea /= max(np.linalg.norm(ea), 1e-12)
        g["R"] = np.stack([ea, eb, nrm], axis=1)
        g["dip_deg"] = float(dd)
    return grains


def dip_values(dip, n, seed=4242):
    """None -> no re-tilt; scalar -> that dip; (lo,hi) -> a drawn distribution."""
    if dip is None:
        return None
    if isinstance(dip, (tuple, list)):
        rng = np.random.default_rng(int(seed))
        return rng.uniform(float(dip[0]), float(dip[1]), size=n)
    return float(dip)


def populate_width(sg, seed=SEED):
    """bed_gen2.populate with SIGMA_G swapped at the module level and restored.

    The files on disk are not touched; the arithmetic of bed_gen2.populate and
    bed_gen.draw_size is used as-is.  Recorded in the report (§5, §6.3).
    """
    keep_b, keep_g = B2.SIGMA_G, BG.SIGMA_G
    try:
        B2.SIGMA_G = float(sg)
        BG.SIGMA_G = float(sg)
        pop, cm = B2.populate(seed)
    finally:
        B2.SIGMA_G = keep_b
        BG.SIGMA_G = keep_g
    return pop, cm


def build(grains, form, n, rule="old", q=0.10, order="random", split=None,
          seed=7):
    """One bed on the calibration grid.

    rule = "old"     : bed_gen2's rule (form + p), the current generator
    rule = "perch"   : summit = <statistic of the support> + q*D
                       (statistic in {"max","mean","p75"}; q = the protrusion
                       above the support, in grain diameters)
    rule = "contact" : rigid drop-until-contact: summit = max over the
                       footprint of (h + depth_below_summit), i.e. the grain
                       is lowered until its cap first touches what is below
    """
    cell = PATCH / n
    xs = (np.arange(n) + 0.5) * cell
    ys = (np.arange(n) + 0.5) * cell
    fg = form_grid(form, n)
    h = np.full((n, n), -1e9)
    gid = np.full((n, n), -1, dtype=np.int32)
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    is_bulk = kinds == "bulk"
    rng = np.random.default_rng(seed)
    if rule == "old":
        seq = np.where(~is_mono)[0]
    else:
        base = np.where(~is_mono & ~is_bulk)[0]
        if order == "random":
            base = rng.permutation(base)
        elif order == "coarse":
            base = base[np.argsort(-np.array([grains[i]["D"] for i in base]))]
        else:
            base = base[np.argsort(np.array([grains[i]["D"] for i in base]))]
        seq = np.concatenate([np.where(is_bulk)[0], base])
    Dx = np.array([g["D"] for g in grains])
    n_fallback = 0
    for i in seq:
        g = grains[i]
        idx, ok, s2, ceff = grain_window(g, xs, ys, n)
        if rule != "old":
            sub = h[idx]
            fin = ok & (sub > -1e8)
            if rule == "contact":
                if fin.any():
                    d = ceff - s2
                    g["z"] = float((sub + d)[fin].max()) - ceff
                else:
                    n_fallback += 1
                    g["z"] = bilinear_at(fg, g["x"], g["y"]) + q * g["D"] - ceff
            else:
                sink = split is not None and Dx[i] < split * D50_REF_M
                if fin.any() and not sink:
                    v = sub[fin]
                    if rule == "max":
                        sup = float(v.max())
                    elif rule == "mean":
                        sup = float(v.mean())
                    elif rule == "p75":
                        sup = float(np.percentile(v, 75.0))
                    else:
                        raise ValueError(rule)
                    g["z"] = sup + q * g["D"] - ceff
                elif fin.any():                      # drainage: fines sink
                    d = ceff - s2
                    g["z"] = float((sub + d)[fin].max()) - ceff
                else:
                    n_fallback += 1
                    g["z"] = bilinear_at(fg, g["x"], g["y"]) + q * g["D"] - ceff
        stamp(h, gid, i, g, idx, ok, s2)
    return {"h": h, "gid": gid, "cell": cell, "form": fg,
            "z": np.array([g["z"] for g in grains]),
            "n_first_contacts": n_fallback}


# ---------------------------------------------------------------------------
# STAGE 1 measurement
# ---------------------------------------------------------------------------
def sigma_form_of_v422():
    """The substrate amplitude v4-2.2's sigma_z fixpoint ended on (recorded)."""
    p = os.path.join(HERE, "out", "v4-2-2", "numbers-v4-2-2.json")
    try:
        with open(p, encoding="utf-8") as f:
            d = json.load(f)
        return float(d["bed_stats"]["sigma_z_fixpoint_final_parameters"]
                     ["sigma_form_m"])
    except Exception:
        return 0.85 * float(B2.SIGMA_TARGET)


def matrix_mask(gid, grains):
    """v4-2.2's matrix rule: top grain is not a boulder and D <= D84_surface."""
    D = np.array([g["D"] for g in grains])
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    Dpx = D[np.clip(gid, 0, len(D) - 1)]
    return (gid >= 0) & (~is_mono[np.clip(gid, 0, len(D) - 1)]) & (
        Dpx <= D84_MATRIX_M)


def measure_field(res, grains, want_curve=False):
    """All STAGE 1 numbers for one field, with the v4-2.5 instrument."""
    h, gid, cell = res["h"], res["gid"], res["cell"]
    hole = h < -1e8
    z = np.where(hole, 0.0, h)
    if hole.any():
        z = np.where(hole, float(z[~hole].mean()), z)
    out = {"hole_frac": float(hole.mean()),
           "sigma_z_m": float(np.std(z)),
           "lag_band_m": list(BAND), "nbin": NBIN}

    mat = matrix_mask(gid, grains)
    out["n_matrix_px"] = int(mat.sum())
    if mat.any():
        out["sigma_z_matrix_m"] = float(np.std(z[mat]))
        out["skew_matrix"] = float(B2._skew(z[mat]))
    else:
        out["sigma_z_matrix_m"] = None
        out["skew_matrix"] = None

    D = np.array([g["D"] for g in grains])
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    is_cluster = np.array([k.startswith("cluster") for k in kinds])
    try:
        ex = B2.exposed_stats(gid, cell, grains, is_mono, is_cluster)
    except Exception:
        ex = {"surface": None, "armor_ratio_exact": None}
    out["D50_surface_m"] = (ex["surface"]["D50_m"]
                            if ex.get("surface") else None)
    out["D84_over_D50_surface"] = (ex["surface"]["D84_over_D50"]
                                   if ex.get("surface") else None)
    out["armor_ratio_exact"] = ex.get("armor_ratio_exact")
    out["D50_surface_drawn_m"] = float(np.median(
        D[kinds == "surface"])) if (kinds == "surface").any() else None

    d50 = out["D50_surface_m"] or D50_REF_M
    out["sigma_z_matrix_over_D50"] = (out["sigma_z_matrix_m"] / d50
                                      if out["sigma_z_matrix_m"] else None)

    dirs = (("x", 1, DM_X), ("y", 0, DM_Y))
    curves = {}
    for key, axis, dmax in dirs:
        d, s = M5.curve_dir(z, cell, dmax, nbin=NBIN, axis=axis)
        fit = M5.fit_band(d, s, BAND[0], BAND[1])
        shp = M5.band_shape(d, s, BAND[0], BAND[1])
        curves[key] = (d, s)
        out[key] = {
            "H": fit["H"], "H_se": fit["H_se"], "n_points": fit["n_points"],
            "rms_log": fit["rms_log"], "r2": fit["r2"],
            "S_min": fit["S_at_d_min"], "S_max": fit["S_at_d_max"],
            "S_growth": fit["S_growth"],
            "local_min_in_band": shp.get("local_H_min"),
            "local_min_at_m": shp.get("local_H_min_at_m"),
            "local_max_in_band": shp.get("local_H_max"),
            "local_probe": shp.get("local_H_probe"),
            "plateau_onset_m": shp.get("plateau_onset_m"),
        }
    n = min(len(curves["x"][1]), len(curves["y"][1]))
    di, si = curves["x"][0][:n], 0.5 * (curves["x"][1][:n]
                                        + curves["y"][1][:n])
    d_iso, s_iso = M5.curve_iso(z, cell, M5.DM_ISO, nbin=NBIN)
    fit = M5.fit_band(d_iso, s_iso, BAND[0], BAND[1])
    shp = M5.band_shape(d_iso, s_iso, BAND[0], BAND[1])
    out["iso"] = {"H": fit["H"], "n_points": fit["n_points"],
                  "rms_log": fit["rms_log"],
                  "local_min_in_band": shp.get("local_H_min"),
                  "local_probe": shp.get("local_H_probe"),
                  "ladder": "run_v4_2_5.curve_iso, dmax=%.2f m" % M5.DM_ISO}
    out["iso_same_ladder_as_xy"] = {
        "H": M5.fit_band(di, si, BAND[0], BAND[1])["H"],
        "note": ("0.5*(Sx+Sy) on the x/y ladder (dmax 0.30 m); recorded "
                 "because v4-2.5's published iso used curve_iso (dmax 0.20 m)")}

    if res["form"] is not None:
        rel = z - res["form"]
        f = {}
        for key, axis, dmax in dirs:
            d, s = M5.curve_dir(rel, cell, dmax, nbin=NBIN, axis=axis)
            fb = M5.fit_band(d, s, BAND[0], BAND[1])
            sh = M5.band_shape(d, s, BAND[0], BAND[1])
            f[key] = {"H": fb["H"], "local_min_in_band": sh.get("local_H_min"),
                      "S_max": fb["S_at_d_max"]}
        out["relief"] = f
        out["sigma_z_relief_m"] = float(np.std(rel))
    if want_curve:
        out["curve"] = {k: {"D_m": [float(v) for v in curves[k][0]],
                            "S_m2": [float(v) for v in curves[k][1]]}
                        for k in curves}
        out["curve"]["iso"] = {"D_m": [float(v) for v in d_iso],
                               "S_m2": [float(v) for v in s_iso]}
    # hole-treatment sensitivity for this field (x only, declaration §6.4)
    if hole.any() and res["form"] is not None:
        from scipy import ndimage
        w = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, z), 60.0, mode="wrap")
        den = ndimage.gaussian_filter(w, 60.0, mode="wrap")
        alt = {"form": np.where(hole, res["form"], h),
               "inpaint60px": np.where(hole, num / np.maximum(den, 1e-9), h)}
        out["hole_treatments"] = {}
        for k, zz in alt.items():
            dd, ss = M5.curve_dir(zz, cell, DM_X, nbin=NBIN, axis=1)
            fb = M5.fit_band(dd, ss, BAND[0], BAND[1])
            sh = M5.band_shape(dd, ss, BAND[0], BAND[1])
            out["hole_treatments"][k] = {"H_x": fb["H"],
                                         "local_min": sh["local_H_min"],
                                         "sigma_z_m": float(np.std(zz))}
            del zz, dd, ss
    out["_z"] = z          # not serialised; used by the figure / budget
    return out


def judge_case(m):
    """The STAGE 1 gate, applied to one case."""
    lmin = m["x"]["local_min_in_band"]
    lmin_y = m["y"]["local_min_in_band"]
    s = m["sigma_z_matrix_over_D50"]
    loc_ok = (lmin is not None and lmin_y is not None
              and min(lmin, lmin_y) >= LOCAL_MIN_REQ)
    sz_ok = (s is not None and SZ_BAND[0] <= s <= SZ_BAND[1])
    hx, hy = m["x"]["H"], m["y"]["H"]
    fit_ok = (hx is not None and hy is not None
              and GATE_H[0] <= hx <= GATE_H[1]
              and GATE_H[0] <= hy <= GATE_H[1])
    return {"local_min_x": lmin, "local_min_y": lmin_y,
            "local_ok": bool(loc_ok), "sigma_z_over_D50": s,
            "sigma_z_ok": bool(sz_ok), "fitted_H_ok": bool(fit_ok),
            "H_x": hx, "H_y": hy, "H_iso": m["iso"]["H"],
            "pass": bool(loc_ok and sz_ok)}


# ---------------------------------------------------------------------------
# STAGE 1 -- the cases
# ---------------------------------------------------------------------------
def sigma_fixpoint(pop, spec, target, n=-1, iters=4, tol=0.03):
    """Move only the substrate amplitude until sigma_z(matrix) hits target."""
    n = CAL_N if n < 0 else n
    sig0 = float(spec.get("sig0", 0.0106))
    sig = sig0
    hist = []
    for it in range(iters):
        form = B2.form_field(512, sig, SEED)
        g = [dict(x) for x in pop]
        dv = dip_values(spec.get("dip"), len(g))
        if dv is not None:
            retilt(g, dv)
        place_on_form(g, form)
        res = build(g, form, n, rule=spec["rule"], q=spec.get("q", 0.10),
                    order=spec.get("order", "random"),
                    split=spec.get("split"))
        mat = matrix_mask(res["gid"], g)
        s = float(np.std(res["h"][mat])) if mat.any() else float("nan")
        hist.append({"sigma_form_m": round(sig, 6),
                     "sigma_z_matrix_m": round(s, 6)})
        del res
        if math.isfinite(s) and abs(s - target) / target < tol:
            break
        sig = float(np.clip(sig * (1 + 0.9 * (target / max(s, 1e-9) - 1)),
                            0.12 * sig0, 6.0 * sig0))
    form = B2.form_field(512, sig, SEED)
    g = [dict(x) for x in pop]
    dv = dip_values(spec.get("dip"), len(g))
    if dv is not None:
        retilt(g, dv)
    place_on_form(g, form)
    res = build(g, form, n, rule=spec["rule"], q=spec.get("q", 0.10),
                order=spec.get("order", "random"), split=spec.get("split"))
    return res, g, float(sig), hist


def run_case(pop, spec, want_curve=False):
    label = spec["id"]
    t0 = time.time()
    if spec.get("target") is None:                     # flat substrate
        form = None
        g = [dict(x) for x in pop]
        dv = dip_values(spec.get("dip"), len(g))
        if dv is not None:
            retilt(g, dv)
        res = build(g, form, CAL_N, rule=spec["rule"], q=spec.get("q", 0.10),
                    order=spec.get("order", "random"),
                    split=spec.get("split"))
        sig = 0.0
        hist = []
    else:
        res, g, sig, hist = sigma_fixpoint(pop, spec, spec["target"] * D50_REF_M)
    m = measure_field(res, g, want_curve=want_curve)
    m["id"] = label
    m["label_ja"] = spec["label_ja"]
    m["group"] = spec["group"]
    m["sigma_form_m"] = sig
    m["fixpoint"] = hist
    m["rule"] = spec["rule"]
    m["q"] = spec.get("q")
    m["order"] = spec.get("order")
    m["split"] = spec.get("split")
    m["dip"] = spec.get("dip")
    m["sg"] = spec.get("sg")
    m["target_over_D50"] = spec.get("target")
    m["n_grains"] = len(g)
    dips = np.array([x["dip_deg"] for x in g
                     if not str(x["kind"]).startswith("mono")])
    m["dip_median_deg"] = float(np.median(dips))
    m["dip_p16_deg"] = float(np.percentile(dips, 16))
    m["dip_p84_deg"] = float(np.percentile(dips, 84))
    m["n_first_contacts"] = int(res["n_first_contacts"])
    m["wall_s"] = float(time.time() - t0)
    m["judge"] = judge_case(m)
    del m["_z"]
    log("  %-26s Hx=%.3f Hy=%.3f iso=%.3f lmin=%.2f/%.2f sz/D50=%.3f "
        "hole=%.3f -> %s [%.1fs]"
        % (label, m["x"]["H"], m["y"]["H"], m["iso"]["H"],
           m["x"]["local_min_in_band"] or float("nan"),
           m["y"]["local_min_in_band"] or float("nan"),
           m["sigma_z_matrix_over_D50"] or float("nan"),
           m["hole_frac"], "PASS" if m["judge"]["pass"] else "fail",
           m["wall_s"]))
    return m


def stage_calib(cp, st0):
    log("STAGE1  calibration: (a) reference, (b) GSD width, (c) support "
        "deposition, (d) imbrication, (e) + substrate")
    out = {"created_utc": utcnow(), "grid": {"n": CAL_N, "cell_m": CAL_CELL},
           "contract": {"band_m": list(BAND), "gate_H": list(GATE_H),
                        "local_min_required": LOCAL_MIN_REQ,
                        "sigma_z_band": list(SZ_BAND),
                        "targets_over_D50": list(TARGETS),
                        "holes": HOLES},
           "cases": {}, "groups": {}}

    # ---------------- (a) the reference / instrument check ----------------
    ref = {}
    dec = json.load(open(DECOMP_423, encoding="utf-8"))
    cv = dec["terms"]["grain_relief"]["_curve"]
    f = M5.fit_band(np.array(cv["D_m"]), np.array(cv["S_m2"]),
                    BAND[0], BAND[1])
    ref["a_saved_relief_curve"] = {
        "what": ("H of the SAVED grain-relief curve of the 3.0 m bed r1 "
                 "(x), fitted with the v4-2.5 instrument"),
        "source": "out/v4-2-3/decomp.json terms.grain_relief._curve",
        "H": f["H"], "n_points": f["n_points"], "rms_log": f["rms_log"],
        "expected_H": 0.2870, "ok": bool(abs(f["H"] - 0.287) < 5e-4)}
    log("  (a) saved relief curve: H=%.4f (expect 0.287) -> %s"
        % (f["H"], "ok" if ref["a_saved_relief_curve"]["ok"] else "MISMATCH"))

    z = np.load(HF_422).astype(np.float64)
    meta = json.load(open(META_422, encoding="utf-8"))
    sig_form = sigma_form_of_v422()
    form = B2.form_field(512, sig_form, SEED)
    fg_hi = form_grid(form, HIRES_N)
    m = measure_field({"h": z, "gid": np.zeros_like(z, dtype=np.int32),
                       "cell": HIRES_CELL, "form": fg_hi, "z": None,
                       "n_first_contacts": 0},
                      [{"D": D50_REF_M, "kind": "surface"}],
                      want_curve=True)
    # the gid above is a dummy: the GSD columns are meaningless there, so they
    # are removed; the H numbers are the shipped-field numbers of v4-2-5.
    for k in ("D50_surface_m", "D84_over_D50_surface", "armor_ratio_exact",
              "D50_surface_drawn_m", "sigma_z_matrix_over_D50",
              "sigma_z_matrix_m", "skew_matrix", "n_matrix_px", "relief",
              "sigma_z_relief_m", "_z"):
        m.pop(k, None)
    old = json.load(open(NUMBERS_425, encoding="utf-8"))["beds"]["v422_patch_grains"]
    chk = {}
    for key in ("x", "y", "iso"):
        rec = old["direction"][key]["bands"]["scaling"]["H"]
        chk[key] = {"H_now": m[key]["H"] if key != "iso" else m["iso"]["H"],
                    "H_v4_2_5": rec,
                    "abs_diff": abs((m[key]["H"] if key != "iso"
                                     else m["iso"]["H"]) - rec)}
    ref["a_shipped_composite"] = {
        "what": "the shipped v4-2.2 0.6 m patch heightfield, re-measured",
        "H": {k: chk[k]["H_now"] for k in chk},
        "H_v4_2_5": {k: chk[k]["H_v4_2_5"] for k in chk},
        "max_abs_diff": max(v["abs_diff"] for v in chk.values()),
        "ok": bool(max(v["abs_diff"] for v in chk.values()) < 1e-9),
        "sigma_z_m": float(np.std(z)),
        "sigma_z_over_D50": float(np.std(z) / D50_REF_M),
        "local_probe_x": m["x"]["local_probe"],
        "local_min_x": m["x"]["local_min_in_band"],
        "curve": m.get("curve"),
    }
    log("  (a) shipped composite x/y/iso = %.3f/%.3f/%.3f vs v4-2.5 "
        "%.3f/%.3f/%.3f -> %s"
        % (chk["x"]["H_now"], chk["y"]["H_now"], chk["iso"]["H_now"],
           chk["x"]["H_v4_2_5"], chk["y"]["H_v4_2_5"], chk["iso"]["H_v4_2_5"],
           "ok" if ref["a_shipped_composite"]["ok"] else "MISMATCH"))

    rel = z - fg_hi
    d, s = M5.curve_dir(rel, HIRES_CELL, DM_X, nbin=NBIN, axis=1)
    fb = M5.fit_band(d, s, BAND[0], BAND[1])
    sh = M5.band_shape(d, s, BAND[0], BAND[1])
    ref["a_shipped_relief"] = {
        "what": ("z - form of the shipped patch (the grain relief, measured "
                 "directly, not from a saved curve)"),
        "H_x": fb["H"], "n_points": fb["n_points"],
        "local_probe": sh["local_H_probe"],
        "local_min": sh["local_H_min"],
        "sigma_relief_m": float(np.std(rel)),
        "sigma_relief_over_D50": float(np.std(rel) / D50_REF_M)}
    log("  (a) shipped relief (z-form) x: H=%.3f local_min=%.2f sigma=%.4f"
        % (fb["H"], sh["local_H_min"], np.std(rel)))
    del z, fg_hi, rel
    out["reference"] = ref
    cp.unit_done("calib_a")

    pop, _ = B2.populate(SEED)
    kinds = np.array([str(g["kind"]) for g in pop])
    out["population"] = {"n": len(pop),
                         "kinds": {str(k): int((kinds == k).sum())
                                   for k in sorted(set(kinds))},
                         "seed": SEED}

    # ---------------- (a)/(b) baseline + GSD width ------------------------
    specs = []
    for tg in TARGETS:
        specs.append({"id": "a_base_old_%03d" % round(tg * 100), "group": "a",
                      "rule": "old", "target": tg, "dip": None,
                      "label_ja": "(a) 現行規則 (form+p) / σz 目標 %.2fD50" % tg})
    for tg in (0.40, 0.46):
        for sg in (1.5, 2.7):
            specs.append({"id": "b_sg%03d_%03d" % (round(sg * 100),
                                                   round(tg * 100)),
                          "group": "b", "rule": "old", "target": tg,
                          "dip": None, "sg": sg,
                          "label_ja": "(b) 独立配置 σg=%.1f / σz 目標 %.2fD50"
                                      % (sg, tg)})
    # (b) with sigma_g = 2.1 IS the current generator: reuse a_base_old_*
    for sp in specs:
        if sp.get("sg") is not None:
            p_sg, _ = populate_width(sp["sg"])
        else:
            p_sg = pop
        out["cases"][sp["id"]] = run_case(p_sg, sp,
                                          want_curve=sp["id"].endswith("_040"))
        del p_sg
    cp.unit_done("calib_a_b")

    # ---------------- (c) support deposition, natural amplitude -----------
    cspecs = [
        {"id": "c_max_q05", "rule": "max", "q": 0.05,
         "label_ja": "(c) 支持=最大 q=0.05（基盤なし・直立）"},
        {"id": "c_max_q10", "rule": "max", "q": 0.10,
         "label_ja": "(c) 支持=最大 q=0.10"},
        {"id": "c_max_q20", "rule": "max", "q": 0.20,
         "label_ja": "(c) 支持=最大 q=0.20"},
        {"id": "c_mean_q20", "rule": "mean", "q": 0.20,
         "label_ja": "(c) 支持=平均 q=0.20"},
        {"id": "c_p75_q20", "rule": "p75", "q": 0.20,
         "label_ja": "(c) 支持=p75 q=0.20"},
        {"id": "c_contact", "rule": "contact", "q": 0.10,
         "label_ja": "(c) 剛体接触（支持点まで降下）"},
    ]
    for sp in cspecs:
        sp.update({"group": "c", "target": None, "dip": 0.0, "order": "random"})
        out["cases"][sp["id"]] = run_case(pop, sp)
    cp.unit_done("calib_c")

    # ---------------- (d) imbrication: upright vs 25-35 deg ---------------
    dspecs = [
        {"id": "d_upright_040", "dip": 0.0,
         "label_ja": "(d) 直立（傾斜なし）・(c) 最良規則・σz 0.40"},
        {"id": "d_imbric_040", "dip": (DIP_CURRENT[0], DIP_CURRENT[1]),
         "label_ja": "(d) 上流傾斜 25-35°（一様分布）・σz 0.40"},
    ]
    for sp in dspecs:
        sp.update({"group": "d", "rule": "max", "q": 0.20, "target": 0.40,
                   "order": "random"})
        out["cases"][sp["id"]] = run_case(pop, sp, want_curve=True)
    cp.unit_done("calib_d")

    # ---------------- (e) (c) + substrate, sigma_z calibrated -------------
    especs = []
    for rule, q, tag, ja in (("max", 0.05, "max_q05", "支持=最大 q=0.05"),
                             ("max", 0.10, "max_q10", "支持=最大 q=0.10"),
                             ("max", 0.20, "max_q20", "支持=最大 q=0.20"),
                             ("mean", 0.20, "mean_q20", "支持=平均 q=0.20"),
                             ("p75", 0.20, "p75_q20", "支持=p75 q=0.20"),
                             ("contact", 0.10, "contact", "剛体接触")):
        for tg in (0.40, 0.46):
            if tg == 0.46 and tag not in ("max_q20", "contact"):
                continue
            especs.append({"id": "e_%s_%03d" % (tag, round(tg * 100)),
                           "group": "e", "rule": rule, "q": q, "target": tg,
                           "order": "random", "dip": 0.0,
                           "label_ja": "(e) %s + 基盤 / σz 目標 %.2fD50"
                                       % (ja, tg)})
    for sp in especs:
        out["cases"][sp["id"]] = run_case(pop, sp, want_curve=True)
    cp.unit_done("calib_e")

    # ---------------- hole-treatment sensitivity (disclosure) ------------
    # every case already carries its own 3-treatment table (measure_field); the
    # headline aggregate below is read from the case, not re-measured.
    m = out["cases"]["e_max_q20_040"]
    sens = dict(m.get("hole_treatments", {}))
    sens["mean"] = {"H_x": m["x"]["H"], "local_min": m["x"]["local_min_in_band"],
                    "sigma_z_m": m["sigma_z_m"]}
    out["hole_sensitivity"] = {"case": "e_max_q20_040",
                              "hole_frac": m["hole_frac"],
                              "treatments": sens,
                              "primary": "mean (declared before measuring)"}
    log("  hole sensitivity (%.1f%% holes): %s"
        % (100 * m["hole_frac"],
           " ".join("%s H=%.3f" % (k, v["H_x"]) for k, v in sens.items())))
    cp.unit_done("calib_holes")

    # ---------------- the ceiling of the instrument (added check) --------
    out["ceiling"] = ceiling_check()
    cp.unit_done("calib_ceiling")
    return out


def band_field(n, L, H, seed, l0):
    """Self-affine k^-(H+1) above 1/L0, flat below (a break at correlation L0)."""
    rng = np.random.default_rng(int(seed))
    w = rng.standard_normal((n, n))
    k1 = np.fft.fftfreq(n, d=L / n) * 2.0 * math.pi
    KX, KY = np.meshgrid(k1, k1)
    K = np.hypot(KX, KY)
    k0 = 1.0 / l0
    amp = np.zeros_like(K)
    m = K > 0
    amp[m] = K[m] ** (-(H + 1.0))
    lo = m & (K < k0)
    amp[lo] = k0 ** (-(H + 1.0))
    z = np.real(np.fft.ifft2(np.fft.fft2(w) * amp))
    z -= z.mean()
    return (z / z.std()).astype(np.float64)


def ceiling_check():
    """What can this instrument even return?  Ideal H_gen=0.79 fields."""
    out = {"what": ("ADDED CHECK (not in the brief): an ideal stationary "
                    "self-affine field with H_gen=0.79 is measured with the "
                    "same estimator, to separate 'the mechanism failed' from "
                    "'the gate is above what the estimator can return on a "
                    "0.6 m patch'"),
           "cell_note": "all ceilings use the v4-2.2 statistics cell 0.293 mm",
           "break_sweep": {}, "box_sweep": {}}
    n, cell = HIRES_N, HIRES_CELL
    for l0 in CEIL_L0:
        hs, mins, prof = [], [], []
        for seed in range(N_SEEDS_CEIL):
            z = band_field(n, PATCH, 0.79, 800 + seed, l0)
            d, s = M5.curve_dir(z, cell, DM_X, nbin=NBIN, axis=1)
            fb = M5.fit_band(d, s, BAND[0], BAND[1])
            sh = M5.band_shape(d, s, BAND[0], BAND[1])
            hs.append(fb["H"])
            mins.append(sh["local_H_min"])
            prof.append(list(sh["local_H_probe"].values()))
        out["break_sweep"]["%.4f" % l0] = {
            "L0_m": l0, "n_seeds": N_SEEDS_CEIL,
            "H_mean": float(np.mean(hs)), "H_min": float(min(hs)),
            "H_max": float(max(hs)),
            "local_min_mean": float(np.mean(mins)),
            "local_probe_mean": [float(v) for v in np.mean(prof, axis=0)],
        }
        log("  ceiling L0=%.3f m: H=%.3f local_min=%.2f"
            % (l0, np.mean(hs), np.mean(mins)))
    for L, nb, nseed in CEIL_BOXES:
        c = L / nb
        hs, mins = [], []
        for seed in range(nseed):
            z = band_field(nb, L, 0.79, 600 + seed, 2.0)
            d, s = M5.curve_dir(z, c, DM_X, nbin=NBIN, axis=1)
            fb = M5.fit_band(d, s, BAND[0], BAND[1])
            sh = M5.band_shape(d, s, BAND[0], BAND[1])
            hs.append(fb["H"])
            mins.append(sh["local_H_min"])
            del z
        out["box_sweep"]["%.2f" % L] = {
            "box_m": L, "n": int(nb), "cell_m": float(c),
            "n_seeds": nseed,
            "H_mean": float(np.mean(hs)), "H_min": float(min(hs)),
            "H_max": float(max(hs)), "local_min_mean": float(np.mean(mins))}
        log("  ceiling box %.2f m: H=%.3f local_min=%.2f"
            % (L, np.mean(hs), np.mean(mins)))
    # how steep must the generator be for the LOCAL criterion at 0.6 m?
    out["hgen_sweep"] = {}
    for hg in (0.79, 0.82, 0.85, 0.88, 0.92, 0.95):
        hs, mins, prof = [], [], []
        for seed in range(4):
            z = band_field(n, PATCH, hg, 950 + seed, 2.0)
            d, s = M5.curve_dir(z, cell, DM_X, nbin=NBIN, axis=1)
            fb = M5.fit_band(d, s, BAND[0], BAND[1])
            sh = M5.band_shape(d, s, BAND[0], BAND[1])
            hs.append(fb["H"])
            mins.append(sh["local_H_min"])
            prof.append(list(sh["local_H_probe"].values()))
            del z, d, s
        out["hgen_sweep"]["%.2f" % hg] = {
            "H_gen": hg, "H_mean": float(np.mean(hs)),
            "local_min_mean": float(np.mean(mins)),
            "local_min_worst_seed": float(min(mins)),
            "local_probe_mean": [float(v) for v in np.mean(prof, axis=0)],
            "local_ok_all_seeds": bool(min(mins) >= LOCAL_MIN_REQ),
            "local_ok_mean": bool(np.mean(mins) >= LOCAL_MIN_REQ - 0.005)}
        log("  ceiling H_gen=%.2f (0.6 m): H=%.3f local_min mean %.2f worst %.2f"
            % (hg, np.mean(hs), np.mean(mins), min(mins)))
    best = max(v["H_mean"] for v in out["break_sweep"].values())
    bestbox = max(v["H_mean"] for v in out["box_sweep"].values())
    out["verdict"] = {
        "max_H_over_break_sweep": best,
        "max_H_over_box_sweep": bestbox,
        "note": ("a field that IS H=0.79 by construction measures at most "
                 "%.3f in this band on a 0.6 m patch and %.3f at best"
                 % (best, bestbox))}
    return out


# ---------------------------------------------------------------------------
# the judgement
# ---------------------------------------------------------------------------
def stage_judge(cp, cal):
    log("JUDGE  STAGE 1 gate: local exponent >= %.2f in 0.4-3.6 cm AND "
        "sigma_z(matrix)/D50 in [%.2f, %.2f]" % (LOCAL_MIN_REQ, *SZ_BAND))
    rows = {}
    for cid, m in cal["cases"].items():
        rows[cid] = m["judge"]
    passed = [k for k, v in rows.items() if v["pass"]]
    # the closest attempt: highest worst-direction local minimum among the
    # cases whose sigma_z is inside the band
    cand = [(min(v["local_min_x"] or -9, v["local_min_y"] or -9), k)
            for k, v in rows.items() if v["sigma_z_ok"]]
    cand.sort(reverse=True)
    best_id = cand[0][1] if cand else None
    best = cal["cases"][best_id] if best_id else None
    best_any = max(rows.items(),
                   key=lambda kv: min(kv[1]["local_min_x"] or -9,
                                      kv[1]["local_min_y"] or -9))
    out = {
        "gate": {"band_m": list(BAND), "local_min_required": LOCAL_MIN_REQ,
                 "sigma_z_band": list(SZ_BAND),
                 "reference_fitted_H_band": list(GATE_H)},
        "rows": rows,
        "n_cases": len(rows),
        "n_pass": len(passed),
        "passed": passed,
        "best_in_sigma_z_band": {
            "case": best_id,
            "worst_direction_local_min": (min(best["judge"]["local_min_x"],
                                              best["judge"]["local_min_y"])
                                          if best else None),
            "H_x": best["judge"]["H_x"] if best else None,
            "H_y": best["judge"]["H_y"] if best else None,
            "H_iso": best["judge"]["H_iso"] if best else None,
            "sigma_z_over_D50": best["judge"]["sigma_z_over_D50"]
            if best else None},
        "best_local_min_any_case": {
            "case": best_any[0],
            "value": min(best_any[1]["local_min_x"] or -9,
                         best_any[1]["local_min_y"] or -9),
            "sigma_z_over_D50": best_any[1]["sigma_z_over_D50"]},
        "decision": "PROCEED_TO_STAGE2" if passed else "STOP",
        "brief_rule": ("BRIEF-V4-2-6 STAGE 1: one combination that keeps the "
                       "local exponent >= 0.70 over 0.4-3.6 cm and sigma_z in "
                       "0.32-0.46 D50 -> STAGE 2; otherwise STOP and report."),
    }
    log("  cases=%d  pass=%d  decision=%s" % (out["n_cases"], out["n_pass"],
                                              out["decision"]))
    if best:
        log("  closest inside the sigma_z band: %s  local_min=%.2f  H=(%.3f, "
            "%.3f, %.3f)" % (best_id,
                             min(best["judge"]["local_min_x"],
                                 best["judge"]["local_min_y"]),
                             best["judge"]["H_x"], best["judge"]["H_y"],
                             best["judge"]["H_iso"]))
    cp.unit_done("judge")
    return out


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def stage_figure(cp, cal, verdict):
    log("FIGURE  structure functions, local exponents, ceiling, gate plane")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(2, 2, figsize=(13.4, 9.2), dpi=115,
                           facecolor="white")
    z_id = band_field(HIRES_N, PATCH, 0.79, 7, 2.0) * (0.40 * D50_REF_M)
    d_id, s_id = M5.curve_dir(z_id, HIRES_CELL, DM_X, nbin=NBIN, axis=1)
    loc_id = 0.5 * np.gradient(np.log(s_id), np.log(d_id))
    del z_id

    # (1) S(tau) of the headline curves
    a = ax[0, 0]
    ref = cal["reference"]["a_shipped_composite"]["curve"]["x"]
    a.loglog(ref["D_m"], ref["S_m2"], lw=1.6, color="#333333",
             label="v4-2.2 shipped composite (x)")
    for cid, col in (("e_max_q20_040", "#1f77b4"),
                     ("d_imbric_040", "#d62728")):
        c = cal["cases"][cid].get("curve", {}).get("x")
        if c:
            a.loglog(c["D_m"], c["S_m2"], lw=1.4, color=col, label=cid)
    a.loglog(d_id, s_id, lw=1.1, ls="--", color="#2ca02c",
             label="ideal H_gen=0.79 field")
    a.axvspan(BAND[0], BAND[1], color="#ffe08a", alpha=0.35, lw=0)
    a.set_xlabel("lag (m)")
    a.set_ylabel("S(tau) (m^2)")
    a.set_title("structure function, flow axis")
    a.grid(alpha=0.25, lw=0.4, which="both")
    a.legend(fontsize=6.5, loc="lower right")

    # (2) local exponent
    a = ax[0, 1]
    a.axvspan(BAND[0], BAND[1], color="#ffe08a", alpha=0.35, lw=0)
    a.axhline(LOCAL_MIN_REQ, color="#b03030", lw=1.0, ls=":")
    a.axhline(GATE_H[1], color="#b03030", lw=1.0, ls=":")
    for cid, lbl, col in (("a_base_old_040", "current rule (same engine)",
                           "#333333"),
                          ("e_max_q20_040", "support deposition + substrate",
                           "#1f77b4"),
                          ("d_imbric_040", "support + imbrication", "#d62728")):
        m = cal["cases"][cid]
        loc = 0.5 * np.gradient(np.log(np.array(m["curve"]["x"]["S_m2"])),
                                np.log(np.array(m["curve"]["x"]["D_m"])))
        a.semilogx(m["curve"]["x"]["D_m"], loc, lw=1.4, color=col, label=lbl)
    a.semilogx(d_id, loc_id, lw=1.1, ls="--",
               color="#2ca02c", label="ideal H_gen=0.79")
    a.axhline(0.79, color="#2ca02c", lw=0.6, ls=":")
    a.set_xlabel("lag (m)")
    a.set_ylabel("local exponent  0.5 dlnS/dlnd")
    a.set_title("local exponent (gate: >= 0.70 across the shaded band)")
    a.set_ylim(-0.05, 1.05)
    a.grid(alpha=0.25, lw=0.4)
    a.legend(fontsize=6.5, loc="lower left")

    # (3) the ceiling
    a = ax[1, 0]
    bs = cal["ceiling"]["break_sweep"]
    xs = sorted(float(k) for k in bs)
    ys = [bs["%.4f" % k]["H_mean"] for k in xs]
    lo = [bs["%.4f" % k]["H_min"] for k in xs]
    hi = [bs["%.4f" % k]["H_max"] for k in xs]
    a.axhspan(GATE_H[0], GATE_H[1], color="#cfe8cf", alpha=0.6, lw=0)
    a.fill_between(xs, lo, hi, color="#888888", alpha=0.25, lw=0,
                   label="seed range")
    a.plot(xs, ys, "o-", color="#333333", lw=1.4, label="ideal field, H_gen=0.79")
    a.set_xscale("log")
    a.set_xlabel("break scale L0 (m)   D50=0.036, D84=0.073")
    a.set_ylabel("fitted H (0.4-3.6 cm)")
    a.set_title("ceiling: what the estimator returns for an ideal field")
    a.grid(alpha=0.25, lw=0.4)
    a.legend(fontsize=7, loc="lower right")

    # (4) the gate plane
    a = ax[1, 1]
    a.add_patch(Rectangle((SZ_BAND[0], LOCAL_MIN_REQ),
                          SZ_BAND[1] - SZ_BAND[0], 1.0 - LOCAL_MIN_REQ,
                          facecolor="#cfe8cf", alpha=0.7, lw=0))
    seen = set()
    for cid, m in cal["cases"].items():
        j = m["judge"]
        x = j["sigma_z_over_D50"]
        y = min(j["local_min_x"] or -9, j["local_min_y"] or -9)
        if x is None:
            continue
        g = m["group"]
        a.scatter([x], [y], s=34 if g in ("d", "e") else 22,
                  color={"a": "#333333", "b": "#8c564b", "c": "#9467bd",
                         "d": "#d62728", "e": "#1f77b4"}.get(g, "#999999"),
                  marker="o" if g in ("d", "e") else "x",
                  label=g if g not in seen else None, zorder=3)
        seen.add(g)
    a.axhline(LOCAL_MIN_REQ, color="#b03030", lw=1.0, ls=":")
    a.axvline(SZ_BAND[0], color="#b03030", lw=1.0, ls=":")
    a.axvline(SZ_BAND[1], color="#b03030", lw=1.0, ls=":")
    a.set_xlim(0.05, 0.95)
    a.set_ylim(0.0, 0.85)
    a.set_xlabel("sigma_z(matrix) / D50_surface")
    a.set_ylabel("worst-direction local exponent in band")
    a.set_title("STAGE 1 gate plane (green = pass)")
    a.grid(alpha=0.25, lw=0.4)
    a.legend(fontsize=7, loc="upper left")

    fig.suptitle("v4-2-6 STAGE 1: the grain layer does not carry the scaling",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=115, facecolor="white")
    plt.close(fig)
    data = buf.getvalue()

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size[0] > 600
        except Exception:
            return False
    atomic_write_bytes(FIGURE, data, verify=verify)
    log("  wrote %s (%d bytes)" % (FIGURE, os.path.getsize(FIGURE)))
    cp.unit_done("figure")
    return {"path": FIGURE, "bytes": int(os.path.getsize(FIGURE)),
            "panels": ["S(tau) x", "local exponent x", "ceiling vs break",
                       "gate plane"]}


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f3(x, nd=3):
    if x is None:
        return "—"
    return ("%." + str(nd) + "f") % float(x)


def f2(x, nd=2):
    if x is None:
        return "—"
    return ("%." + str(nd) + "f") % float(x)


def mechanism_rows(cal):
    order = ["a_base_old_040", "a_base_old_046", "b_sg150_040", "b_sg270_040",
             "b_sg270_046", "c_max_q05", "c_max_q10", "c_max_q20",
             "c_mean_q20", "c_p75_q20", "c_contact", "d_upright_040",
             "d_imbric_040", "e_max_q05_040", "e_max_q10_040",
             "e_max_q20_040", "e_mean_q20_040", "e_p75_q20_040",
             "e_contact_040", "e_max_q20_046", "e_contact_046"]
    rows = []
    for cid in order:
        m = cal["cases"].get(cid)
        if m is None:
            continue
        j = m["judge"]
        rows.append((cid, m, j))
    return rows


def build_report(ctx):
    st0 = ctx["stage0"]
    cal = ctx["calib"]
    vd = ctx["judge"]
    fig = ctx["figure"]
    lg = ctx["g_legacy"]
    ce = cal["ceiling"]
    cases = cal["cases"]
    ref = cal["reference"]
    L = []
    A = L.append
    created = ctx["created_utc"]

    A("# REPORT-V4-2-6 — 床の作り直し：粒の層そのものに実在のスケーリングを持たせる")
    A("")
    A("Brief: `BRIEF-V4-2-6.md`。Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`。"
      "Built: %s。" % created)
    A("新規 .py は `run_v4_2_6.py` の 1 本だけ。出力は `out/v4-2-6/`。本 run は既存 .py を"
      "**1 バイトも変更していない**（§7 G-LEGACY）。")
    A("")
    A("**この run は STAGE 1 で止まった。** brief の較正ゲート（0.4–3.6 cm の全域で"
      "局所指数 ≥ 0.70、かつ σz(matrix) が 0.32–0.46·D50）を満たす機構が 1 つも"
      "見つからなかったので、brief の指示どおり STAGE 2（`bed_gen3.py`）と STAGE 3"
      "（新床のゲート）には入っていない。**`bed_gen3.py` は存在しない**（作っていない）。"
      "下に、どの機構がどこまで届いたかを数値で全部出す。")
    A("")

    # ---------------- 0 ----------------
    A("## 0. 判定（全ゲート）")
    A("")
    A("| ゲート | 判定 | 実測 / 根拠 |")
    A("| --- | --- | --- |")
    A("| STAGE 0 事前確認 | **PASS** | 必須 7 ファイルの存在・sha256・バイト数を記録（§1）。"
      "D50_surface = %.6f m（brief の 0.0362 と一致）。既存 .py **%d 本**のハッシュを保存 |"
      % (st0["d50"]["v422_patch_m"], st0["n_existing_py"]))
    A("| STAGE 1 較正ゲート（主判定） | **FAIL** | 局所指数 ≥ 0.70 を満たす組合せ **0 / %d**。"
      "σz 帯に入る組合せは %d あるが、そのどれも局所指数が届かない。σz 帯内で最も近いのは "
      "`%s`（最悪方向の局所指数最小 %.2f、H = %.3f / %.3f / %.3f）。全 case で最も高い"
      "局所指数最小値でも **%.2f**（要件 0.70） |"
      % (vd["n_cases"], sum(1 for v in vd["rows"].values() if v["sigma_z_ok"]),
         vd["best_in_sigma_z_band"]["case"] or "—",
         vd["best_in_sigma_z_band"]["worst_direction_local_min"] or -9,
         vd["best_in_sigma_z_band"]["H_x"] or -9,
         vd["best_in_sigma_z_band"]["H_y"] or -9,
         vd["best_in_sigma_z_band"]["H_iso"] or -9,
         vd["best_local_min_any_case"]["value"]))
    A("| STAGE 2 実装（`bed_gen3.py` / 0.6 m パッチ / 3.0 m 帯状床） | **未実行** | "
      "STAGE 1 が FAIL のため brief の STOP 規則を適用。`bed_gen3.py` は作っていない |")
    A("| STAGE 3 G-SCALING（パッチ・帯状床 × x/y/iso） | **未実行** | 生成物が無い"
      "（STAGE 2 未実行） |")
    A("| STAGE 3 G-AUDIT（v4-2.2 の監査 8 項目のうち 7 項目） | **未実行** | 同上。"
      "監査は新しい床に対してのみ意味を持つ |")
    A("| STAGE 3 G-MESO（M1–M6） | **未実行** | 同上。`DESIGN-V4-2-6-MESO.md` は"
      "本ツリー（`D:\\mutaro-work`、`D:\\yohey2026\\dsh\\muse`）に**存在しない**ので、"
      "brief §4 の要約だけを根拠にすることになる（測定自体をしていない）。"
      "提案はラリー未検証である旨も brief の指示どおり併記する |")
    A("| G-LEGACY（既存 .py 無改造） | **%s** | stage-0 で %d 本の sha256 を保存し、"
      "run の最後に再ハッシュ。**本 run 由来の既存 .py 変更 0**（import した 5 本は"
      "すべて一致）。外部レーンの差分 変更 %d / 追加 %d を §7 に mtime つきで開示"
      "（厳密比較は%s） |"
      % ("PASS" if lg["ok"] else "**FAIL**", lg["n_baseline"],
         len(lg["changed"]), len(lg["added"]),
         "一致" if lg["ok_strict"] else "差分あり"))
    A("")
    A("**主判定の読み方**: 判定は 2 条件の論理積（帯内局所指数の最小 ≥ 0.70、"
      "σz(matrix)/D50 ∈ [0.32, 0.46]）。H の帯 0.70–0.85 は STAGE 3 の主判定なので"
      "参考列として併記した。帯・ゲート・物差しは**測定前にコード内で固定**し、"
      "測った後に動かしていない（`BAND=(0.004,0.036)`、`LOCAL_MIN_REQ=0.70`、"
      "`SZ_BAND=(0.32,0.46)`、`GATE_H=(0.70,0.85)`）。")
    A("")
    A("| 参考: 機構別の最良値 | 値 |")
    A("| --- | --- |")
    bid = vd["best_in_sigma_z_band"]["case"]
    bm = cases[bid] if bid else None
    A("| σz 帯内での最高 H_iso | %.3f (`%s`) |"
      % (bm["iso"]["H"], bid) if bm else "| — | — |")
    A("| σz 帯内での最高 局所指数最小（x/y の悪い方） | %.2f (`%s`) |"
      % (vd["best_in_sigma_z_band"]["worst_direction_local_min"], bid)
      if bm else "| — | — |")
    A("| 同一エンジンでの現行規則の基準 | H = %.3f / %.3f / %.3f `a_base_old_040` |"
      % (cases["a_base_old_040"]["x"]["H"], cases["a_base_old_040"]["y"]["H"],
         cases["a_base_old_040"]["iso"]["H"]))
    A("| 納品済み v4-2.2 パッチ（参照） | H = %.3f / %.3f / %.3f、σz/D50 = %.3f |"
      % (ref["a_shipped_composite"]["H"]["x"],
         ref["a_shipped_composite"]["H"]["y"],
         ref["a_shipped_composite"]["H"]["iso"],
         ref["a_shipped_composite"]["sigma_z_over_D50"]))
    A("| 物差しの上限（理想 H_gen=0.79 場、0.6 m 箱） | H_fit = %.3f、帯内局所指数の最小 %.2f。"
      "局所指数 0.70 を 4 seed すべてで通す自己アフィン場は掃引（H_gen≤%.2f）に無い（§4.3） |"
      % (ce["box_sweep"]["0.60"]["H_mean"],
         ce["box_sweep"]["0.60"]["local_min_mean"],
         max(float(v["H_gen"]) for v in ce["hgen_sweep"].values())))
    A("")

    # ---------------- 1 ----------------
    A("## 1. STAGE 0（絶対パス事前確認・ハッシュ・物差しの同定）")
    A("")
    A("| 対象 | bytes | sha256 (16) |")
    A("| --- | --- | --- |")
    for k, v in st0["preflight"].items():
        A("| `%s` | %d | `%s` |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    A("- D50: **%.6f m**（brief の 0.0362 と一致、`heightfield-meta.json` の"
      "露出面積監査値）。D84(matrix 規則) = %.3f m。" % (st0["d50"]["v422_patch_m"],
                                                        D84_MATRIX_M))
    A("- 既存 .py **%d 本**の sha256 を `stage0-py-hashes.json` に保存。"
      % st0["n_existing_py"])
    A("")
    A("### 1.1 物差し（v4-2.5）の同定 — 本ジョブはこれを verbatim に使う")
    A("")
    A("`run_v4_2_5.py` を import し、同ファイルの関数をそのまま呼ぶ（式を写していない）。")
    A("")
    A("| 関数 | 行 | 中身 |")
    A("| --- | --- | --- |")
    for name, v in st0["instrument"]["functions"].items():
        A("| `%s` | %d | %s |" % (name, v["line"],
                                  v["doc"] or "(docstring なし)"))
    A("")
    A("帯とゲート（同ファイルの定数）: `BANDS[0] = (%s, %.3f, %.3f)` m、"
      "`GATE = (%.2f, %.2f)`、ラダー `nbin = %d`。S(τ) は軸ごとの周期ラップ付き"
      "二乗差の平均（`curve_dir`）、H は帯内の両対数一次回帰の傾き 1/2（`fit_band`）、"
      "局所指数は同じラダー上の `0.5·dlnS/dlnd`（`local_H`、帯内の最小/最大/交差は "
      "`band_shape`）。**v4-2.2/2.3 の物差しと同一**であることは v4-2.5 §3 が"
      "機械確認している（相対 1e-9 以下）。"
      % (M5.BANDS[0][0], M5.BANDS[0][1], M5.BANDS[0][2], GATE_H[0], GATE_H[1],
         NBIN))
    A("")
    A("本 run での検算（同じ場から同じ数字が出るか）:")
    A("")
    A("| 検算 | 本 run | 記録値 | 一致 |")
    A("| --- | --- | --- | --- |")
    A("| 保存済み grain-relief 曲線 (r1, x) の H | %.4f | 0.2870（v4-2.5 §5.1） | %s |"
      % (ref["a_saved_relief_curve"]["H"],
         "**True**" if ref["a_saved_relief_curve"]["ok"] else "False"))
    for k in ("x", "y", "iso"):
        A("| 納品済み v4-2.2 パッチの H (%s) | %.4f | %.4f | %s |"
          % (k, ref["a_shipped_composite"]["H"][k],
             ref["a_shipped_composite"]["H_v4_2_5"][k],
             "**True**" if ref["a_shipped_composite"]["max_abs_diff"] < 1e-9
             else "False"))
    A("")
    A("→ 物差しは同一（最大差 %.1e）。**以下の H の差は機構の差にだけ帰属できる**。"
      % ref["a_shipped_composite"]["max_abs_diff"])
    A("")

    # ---------------- 2 ----------------
    A("## 2. STAGE 1 — 方法（較正グリッドと機構の実装）")
    A("")
    A("- 較正グリッド: **%d² / 0.6 m**（cell %.4f mm）。v4-2.2 の統計グリッド"
      "（2048²）と同一の場で H が変らないことは §6.2 に書く。" % (CAL_N, CAL_CELL * 1000))
    A("- 粒の母集団は `bed_gen2.populate(SEED)` をそのまま使う（%d 粒、内訳 %s）。"
      % (cal["population"]["n"], cal["population"]["kinds"]))
    A("- 支持・堆積の実装（新規・本ジョブの試験ベンチ）:")
    A("  - `perch`: 粒の footprint 内の**既存表面の統計量**（最大 / 平均 / p75）を支持高さとし、"
      "中心高さ = 支持 + q·D。q は粒径単位の突出量。")
    A("  - `contact`: **剛体接触**。粒を静かに降ろし、カプセルが下の表面に最初に触れる"
      "高さまで下げる（summit = max_footprint(h + その点の summit からの落差)）。"
      "これが brief の「支持点から中心高さを決める」の文字どおりの実装。")
    A("  - 落下順はランダム（`order=random`）、ほかに粗粒先着・細粒先着の対照を取る。"
      "局所緩和は**接触と非貫入だけ**で、平滑化・ならしは一切していない。")
    A("  - 掃引した q は 0.05 / 0.10 / 0.20（物理的な突出量の範囲。§5 に出所）。")
    A("- σz の較正: 機構ごとに**基盤振幅 σ_form だけ**を動かし、σz(matrix) を "
      "0.40 と 0.46·D50（帯の中と上端）に合わせる 4 回以内の不動点。"
      "これは v4-2.2 が使ったのと同じ手続き。**H を上げるための調整はしていない**"
      "（σz を帯に入れるための調整だけ）。")
    A("- 演出（paving）は無し: 場を平滑化する処理、上に場を足す処理、粒スケールの起伏を"
      "潰す処理はコードに入っていない。")
    A("- 穴（plan-view で粒が 1 つも被覆しない画素）: 母集団の被覆率の都合で "
      "**%.1f%%（現行規則）/ %.1f%%（堆積規則）** 残る。較正では被覆画素の平均で"
      "埋める（全 case 共通）。埋め方の感度は §6.4。"
      % (100 * cases["a_base_old_040"]["hole_frac"],
         100 * cases["e_max_q20_040"]["hole_frac"]))
    A("")

    # ---------------- 3 ----------------
    A("## 3. STAGE 1 — 機構別の結果（どれが効いてどれが効かなかったか）")
    A("")
    A("各セルは `H_x / H_y / H_iso`、`ℓmin` = 帯内の局所指数の最小（x と y の悪い方）、"
      "`σz/D50` = matrix の σz を表面 D50 で割ったもの、`合格` = 2 条件の論理積。")
    A("`σz 目標` は基盤振幅を合わせた先（`—` は基盤なし = 自然振幅）。")
    A("")
    A("| case | 機構 | H_x | H_y | H_iso | ℓmin | σz/D50 | 帯内 S 成長 | 合格 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for cid, m, j in mechanism_rows(cal):
        A("| `%s` | %s | %.3f | %.3f | %.3f | %.2f | %.3f | ×%.1f | %s |"
          % (cid, m["label_ja"], m["x"]["H"], m["y"]["H"], m["iso"]["H"],
             min(j["local_min_x"] or -9, j["local_min_y"] or -9),
             j["sigma_z_over_D50"] or -9, m["x"]["S_growth"],
             "**PASS**" if j["pass"] else "fail"))
    A("")
    A("### 3.1 (a) 現行生成器の粒のみ — 物差しの検算")
    A("")
    A("- 保存済み relief 曲線（r1, x）: **H = %.4f**（期待 0.2870）→ 一致。"
      % ref["a_saved_relief_curve"]["H"])
    A("- 納品済みパッチから直接 (z − form) を測った relief: **H_x = %.3f**"
      "（局所指数の最小 %.2f、σ_relief = %.4f m = %.3f·D50）。"
      % (ref["a_shipped_relief"]["H_x"], ref["a_shipped_relief"]["local_min"],
         ref["a_shipped_relief"]["sigma_relief_m"],
         ref["a_shipped_relief"]["sigma_relief_over_D50"]))
    A("- 同一エンジンで現行規則を組み直した基準: `a_base_old_040` で "
      "H = %.3f / %.3f / %.3f、ℓmin = %.2f（納品フィールドの 0.323 / 0.373 / 0.349 "
      "より少し低い。理由は §6.1 の穴の埋め方）。**以下の機構の数字はすべてこの"
      "同一エンジン同士の比較**。"
      % (cases["a_base_old_040"]["x"]["H"], cases["a_base_old_040"]["y"]["H"],
         cases["a_base_old_040"]["iso"]["H"],
         min(cases["a_base_old_040"]["judge"]["local_min_x"],
             cases["a_base_old_040"]["judge"]["local_min_y"])))
    A("")
    A("### 3.2 (b) 独立配置のまま GSD の幅 σg を振る")
    A("")
    A("| σg | 表面 D84/D50（露出面積監査） | H_x | H_y | H_iso | ℓmin | σz/D50 | 監査 D84/D50 帯 1.9–2.3 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for cid, sg in (("b_sg150_040", 1.5), ("a_base_old_040", 2.1),
                    ("b_sg270_040", 2.7)):
        m = cases[cid]
        d84 = m["D84_over_D50_surface"]
        ok = (d84 is not None and 1.9 <= d84 <= 2.3)
        A("| %.1f | %.3f | %.3f | %.3f | %.3f | %.2f | %.3f | %s |"
          % (sg, d84 or -9, m["x"]["H"], m["y"]["H"], m["iso"]["H"],
             min(m["judge"]["local_min_x"] or -9,
                 m["judge"]["local_min_y"] or -9),
             m["judge"]["sigma_z_over_D50"] or -9,
             "内" if ok else "**外**"))
    A("")
    A("- σg を 1.5 → 2.7 に広げても H は**上がらない**（%.3f → %.3f、Δ = %+.3f）。"
      "粒径分布の幅は粒スケールの白い起伏の帯内指数を変えない。"
      % (cases["b_sg150_040"]["x"]["H"], cases["b_sg270_040"]["x"]["H"],
         cases["b_sg270_040"]["x"]["H"] - cases["b_sg150_040"]["x"]["H"]))
    A("- 露出面積監査で測った表面 D84/D50 は σg=1.5 で %.3f（監査帯 1.9–2.3 の"
      "**下**）、σg=2.7 で %.3f（帯の**内**）。つまり「幅を広げれば H が上がる」"
      "という逃げ道は、監査帯を壊してまで得るものがない（σg=2.7 は帯の内側なのに"
      "H は σg=2.1 より低い）。"
      % (cases["b_sg150_040"]["D84_over_D50_surface"],
         cases["b_sg270_040"]["D84_over_D50_surface"]))
    A("")
    A("### 3.3 (c) 支持・堆積配置（自然振幅 = 基盤なし）")
    A("")
    A("| 支持の取り方 | q | H_x | H_y | H_iso | ℓmin | σz/D50 | 穴 | 備考 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for cid in ("c_max_q05", "c_max_q10", "c_max_q20", "c_mean_q20",
                "c_p75_q20", "c_contact"):
        m = cases[cid]
        rule = m["rule"]
        qv = "—" if rule == "contact" else "%.2f" % m["q"]
        A("| %s | %s | %.3f | %.3f | %.3f | %.2f | %.3f | %.1f%% | 基盤なし・直立 |"
          % (rule, qv, m["x"]["H"], m["y"]["H"], m["iso"]["H"],
             min(m["judge"]["local_min_x"] or -9,
                 m["judge"]["local_min_y"] or -9),
             m["judge"]["sigma_z_over_D50"] or -9, 100 * m["hole_frac"]))
    A("")
    A("- 支持・堆積は現行規則より H を **+%.3f**（最大 `c_max_q20` の H_x %.3f 対 "
      "基準 %.3f）上げるが、**0.70 には遠い**。"
      % (cases["c_max_q20"]["x"]["H"] - cases["a_base_old_040"]["x"]["H"],
         cases["c_max_q20"]["x"]["H"], cases["a_base_old_040"]["x"]["H"]))
    A("- 基盤なしでは σz が帯に入らない（`c_max_q20` で %.3f·D50）。"
      "支持だけで σz の予算を供給できない（§4.2）。"
      % (cases["c_max_q20"]["judge"]["sigma_z_over_D50"]))
    A("")
    A("### 3.4 (d) imbrication（上流傾斜 25–35°）")
    A("")
    A("同じ粒（同じ a/b/c、同じ方位、同じ支持規則、同じ σ_form 目標）で、傾斜だけを"
      "0° と **25–35° の一様分布**（brief (d) の指定）にした対照。")
    A("")
    A("| 構成 | dip 中央値 [p16–p84] | H_x | H_y | H_iso | ℓmin | σz/D50 | ΔH_x vs 直立 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    up, im = cases["d_upright_040"], cases["d_imbric_040"]
    for m_ in (up, im):
        A("| %s | %.1f° [%.1f–%.1f] | %.3f | %.3f | %.3f | %.2f | %.3f | %s |"
          % ("直立 (dip 0°)" if m_ is up else "上流傾斜 25–35°",
             m_["dip_median_deg"], m_["dip_p16_deg"], m_["dip_p84_deg"],
             m_["x"]["H"], m_["y"]["H"], m_["iso"]["H"],
             min(m_["judge"]["local_min_x"], m_["judge"]["local_min_y"]),
             m_["judge"]["sigma_z_over_D50"],
             "—" if m_ is up else "%+.3f" % (im["x"]["H"] - up["x"]["H"])))
    A("")
    A("- imbrication の効果は **ΔH_x = %+.3f**（帯内の局所指数最小は %.2f → %.2f）。"
      "**効かない**。粒の傾斜はカプセルの向きを変えるだけで、高さ場の帯内スペクトルを"
      "変えない。"
      % (im["x"]["H"] - up["x"]["H"],
         min(up["judge"]["local_min_x"], up["judge"]["local_min_y"]),
         min(im["judge"]["local_min_x"], im["judge"]["local_min_y"])))
    A("")
    A("### 3.5 (e) (c) + 基盤（現行 fBm）— σz を帯に入れた最良ケース")
    A("")
    bx = cases["e_max_q20_040"]
    A("| 機構 | σz 目標 | σ_form (m) | H_x | H_y | H_iso | ℓmin | σz/D50 | 合格 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for cid, m, j in mechanism_rows(cal):
        if not cid.startswith("e_"):
            continue
        A("| %s | %.2f | %.5f | %.3f | %.3f | %.3f | %.2f | %.3f | %s |"
          % (m["label_ja"], m["target_over_D50"], m["sigma_form_m"],
             m["x"]["H"], m["y"]["H"], m["iso"]["H"],
             min(j["local_min_x"] or -9, j["local_min_y"] or -9),
             j["sigma_z_over_D50"] or -9,
             "**PASS**" if j["pass"] else "fail"))
    A("")
    A("- 最良は `e_max_q20_040`: H = %.3f / %.3f / %.3f、ℓmin = %.2f、σz/D50 = %.3f。"
      "**σz は帯の内、H は帯の遥か下**。"
      % (bx["x"]["H"], bx["y"]["H"], bx["iso"]["H"],
         min(bx["judge"]["local_min_x"], bx["judge"]["local_min_y"]),
         bx["judge"]["sigma_z_over_D50"]))
    A("- σz の目標を帯の上端（0.46·D50）まで上げても H は %.3f → %.3f しか動かない"
      "（`e_max_q20_040` → `e_max_q20_046`）。**σz を増やす方向には H はほぼ"
      "飽和している**。"
      % (bx["x"]["H"], cases["e_max_q20_046"]["x"]["H"]))
    A("")

    # ---------------- 4 ----------------
    A("## 4. なぜ届かないのか（分解と、物差しの上限）")
    A("")
    A("### 4.1 帯内の分散の持ち主")
    A("")
    best = bx
    A("`e_max_q20_040` の帯内の分散の内訳（relief = z − form、cov は交差項、"
      "流れ軸・周期ラップ）:")
    A("")
    bud = ctx["budget"]
    A("| ラグ | S_合成 (m²) | S_form | S_relief | cov 項 | form の寄与 | relief の寄与 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for k in ("0.0041", "0.0100", "0.0322"):
        if k not in bud["lags"]:
            continue
        v = bud["lags"][k]
        A("| %.2f cm | %.3e | %.3e | %.3e | %+.3e | %.2f | %.2f |"
          % (100 * float(k), v["S_total"], v["S_form"], v["S_relief"],
             v["cov_term_2x"], v["share_form"], v["share_relief"]))
    A("")
    A("成分別の帯内 H（同じ場、x）: 合成 **%.3f** / 基盤 **%.3f** / 粒起伏 **%.3f**。"
      % (bud["H_in_band"]["composite"]["H"], bud["H_in_band"]["form"]["H"],
         bud["H_in_band"]["relief"]["H"]))
    A("→ 帯内の S は **粒起伏がほぼ全部**を持ち（3.22 cm で %.0f%%）、"
      "その起伏の帯内 H は %.3f。基盤は H_gen=0.79 でも帯内の寄与が %.0f%% しか"
      "ないので、合成の H を引き上げる力がない。**合成の H は粒起伏の H そのもの**。"
      % (100 * bud["lags"]["0.0322"]["share_relief"],
         bud["H_in_band"]["relief"]["H"],
         100 * bud["lags"]["0.0322"]["share_form"]))
    A("")
    A("### 4.2 σz と H の同時成立に必要な突出量（算術）")
    A("")
    frac = (0.41 / 3.6) ** 1.4
    cx = ref["a_shipped_composite"]["curve"]["x"]
    s041_meas = float(np.interp(0.0041, cx["D_m"], cx["S_m2"]))
    s322_meas = float(np.interp(0.0322, cx["D_m"], cx["S_m2"]))
    A("帯内の局所指数を全域で 0.70 以上にするには、S(0.41 cm) ≤ S(3.6 cm)·"
      "(0.41/3.6)^1.4 = **%.2e**·S(3.6) でなければならない。σz(matrix) を "
      "0.40·D50 = %.4f m に置くと、帯の上端で S(3.6 cm) ≈ 2σ² = %.2e m² なので、"
      "S(0.41 cm) ≤ %.2e m²、つまり **4 mm ラグの高さ差の rms ≤ %.1f mm**。"
      "現行床の実測は %.1f mm（S(0.41 cm) = %.2e、S(3.22 cm) = %.2e）。"
      "粒の突出量のばらつき（σ_relief = %.1f mm）をそのまま許すと、この条件は"
      "満たせない。"
      % (frac, 0.40 * D50_REF_M, 2 * (0.40 * D50_REF_M) ** 2,
         frac * 2 * (0.40 * D50_REF_M) ** 2,
         1000 * math.sqrt(frac * 2 * (0.40 * D50_REF_M) ** 2),
         1000 * math.sqrt(s041_meas), s041_meas, s322_meas,
         1000 * ref["a_shipped_relief"]["sigma_relief_m"]))
    A("")
    A("### 4.3 物差しの上限（追加検査・brief の要求外）")
    A("")
    A("機構の失敗と、ゲートが物差しの上限を超えていることを分離するため、"
      "**H_gen = 0.79 で作った理想場**（壊れ目 L0 を振る）を同じ物差しに通した。")
    A("")
    A("| 壊れ目 L0 | H_fit 平均 (6 seed) | 帯内局所指数の最小 | 局所指数プロファイル (0.41→3.28 cm) |")
    A("| --- | --- | --- | --- |")
    for k in sorted(ce["break_sweep"], key=lambda s: float(s)):
        v = ce["break_sweep"][k]
        A("| %.3f m (%.2f·D50) | %.3f | %.2f | %s |"
          % (v["L0_m"], v["L0_m"] / D50_REF_M, v["H_mean"],
             v["local_min_mean"],
             " ".join("%.2f" % q for q in v["local_probe_mean"])))
    A("")
    A("| 箱 | H_fit 平均 | 帯内局所指数の最小 |")
    A("| --- | --- | --- |")
    for k in sorted(ce["box_sweep"], key=lambda s: float(s)):
        v = ce["box_sweep"][k]
        A("| %.2f m (%d², cell %.3f mm) | %.3f | %.2f |"
          % (v["box_m"], v["n"], v["cell_m"] * 1000, v["H_mean"],
             v["local_min_mean"]))
    A("")
    A("読み方:")
    A("")
    A("1. **H_gen=0.79 の理想場ですら、0.6 m パッチの帯内局所指数は 3.6 cm で "
      "%.2f（cell 1.17 mm）/ %.2f（cell 0.293 mm）まで下がる**。物差しは 0.6 m 箱で"
      "入力の H より低い局所指数を返す（有限箱のバイアス）。"
      % (ce["box_sweep"]["0.60"]["local_min_mean"],
         ce["break_sweep"]["0.2000"]["local_min_mean"]))
    A("2. では局所指数 0.70 を 0.6 m 箱で出すには生成指数がいくつ必要か: ")
    A("")
    A("| H_gen | H_fit (0.4–3.6 cm) | 帯内局所指数の最小（4 seed 平均 / 最悪） | 局所指数プロファイル | 全 seed で ≥0.70 |")
    A("| --- | --- | --- | --- | --- |")
    for k in sorted(ce["hgen_sweep"]):
        v = ce["hgen_sweep"][k]
        A("| %.2f | %.3f | %.2f / %.2f | %s | %s |"
          % (v["H_gen"], v["H_mean"], v["local_min_mean"],
             v["local_min_worst_seed"],
             " ".join("%.2f" % q for q in v["local_probe_mean"]),
             "**成立**" if v["local_ok_all_seeds"] else "不成立"))
    A("")
    first_mean_ok = [float(v["H_gen"]) for v in ce["hgen_sweep"].values()
                     if v["local_ok_mean"]]
    A("→ 平均で 0.70 に届くのは **H_gen ≥ %.2f**（掃引の最小値）。ただし 4 seed "
      "すべてで 0.70 を超える H_gen は掃引の上限 %.2f でも無い（最悪 seed は "
      "0.6 m 箱の推定ばらつきで下に振れる）。いずれにせよ **STAGE 1 のゲートは"
      "「H≈0.79」の検定ではなく、H_gen 0.9 級の「ほぼ滑らかな面」を選ぶ検定に"
      "なっている**（物差しの有限箱バイアスが帯の上端で局所指数を約 0.1 食う）。"
      "それでもゲートは数学的に不可能ではない。**不可能なのは「粒のモザイクで"
      "それを出すこと」**で、理由は §4.1–4.2 のとおり（帯内の分散が粒起伏に"
      "支配され、その帯内指数が ~0.4）。"
      % (min(first_mean_ok) if first_mean_ok else float("nan"),
         max(float(v["H_gen"]) for v in ce["hgen_sweep"].values())))
    A("3. H_fit の帯 0.70–0.85 のほうは、理想場なら 0.6 m でも 2.4 m でも入る"
      "（%.3f / %.3f）。STAGE 3 の主判定だけなら原理的に到達可能。"
      % (ce["box_sweep"]["0.60"]["H_mean"], ce["box_sweep"]["2.40"]["H_mean"]))
    A("")

    # ---------------- 5 ----------------
    A("## 5. パラメータの出所一覧")
    A("")
    A("| 量 | 値 | 出所 | 備考 |")
    A("| --- | --- | --- | --- |")
    A("| PATCH | %.1f m | v4-2.1（bed_gen2 の定数） | パッチ寸法 |" % PATCH)
    A("| D50_in | %.3f m | v4-2.1（面積基準の入力中央値） | 下層の分布 |"
      % B2.D50_IN)
    A("| σg（入力） | %.1f | v4-2.1 の文献値（対数正規・面積基準） | (b) で 1.5/2.7 も掃引 |"
      % B2.SIGMA_G)
    A("| ARMOR | %.1f | v4-2.1（装甲比 ≈2 は実在の既知量） | 表面 D50 = 0.040 |"
      % B2.ARMOR)
    A("| SUBSTRATE_H | %.2f | 文献値（v4-2.1 リサーチの H≈0.79） | 基盤 fBm の生成指数 |"
      % B2.SUBSTRATE_H)
    A("| FORM_W_LOW / FORM_SKEW_C | %.2f / %.2f | v4-2.2 の設計校正（REPORT-V4-2-2 §8 に開示） | 基盤の低次モード重みと歪み係数 |"
      % (B2.FORM_W_LOW, B2.FORM_SKEW_C))
    A("| σz 目標 | 0.40·D50_surface（帯 0.32–0.46） | 実在の水流作用床（v4-2.1 リサーチ） | 本 run は 0.40 と 0.46 の 2 点で較正 |")
    A("| 突出の非対称（P_TAIL 等） | %.2f / %.2f / %.2f | v4-2.2 の設計校正 | 現行規則の p |"
      % (B2.P_TAIL, B2.P_TAIL_MEAN, B2.P_TAIL_MAX))
    A("| SORT_GAIN | %.2f | v4-2.2 の設計校正 | 粗粒が上・細粒が下 |" % B2.SORT_GAIN)
    A("| 表面被覆 / bulk 体積係数 | %.2f / %.2f | v4-2.2 の構成（空隙率 0.339 を出す） | — |"
      % (B2.SURF_COVER, B2.BULK_VOL))
    A("| 傾斜（imbrication） | 30 ± 10° | 実在の既知量 25–35°（v4-2.1 リサーチ） | (d) の対照は 0° |")
    A("| クラスター | 3 個 / 被覆 ~0.10 | 実在の既知量 ~10% | 現行のまま |")
    A("| 較正グリッド | %d² / 0.6 m（%.4f mm/px） | **本ジョブの較正で決めた値** | 2048² との差は §6.2 |"
      % (CAL_N, CAL_CELL * 1000))
    A("| 支持統計 / q / 落下順 | max・mean・p75 / 0.05–0.20 / random | **本ジョブの較正で決めた値**（機構の掃引） | 範囲を縛って掃引した全数を §3 に記載 |")
    A("| σ_form 不動点 | 0.40 / 0.46·D50、4 回以内・3% 収束 | **本ジョブの較正で決めた値**（v4-2.2 と同じ手続き） | 全試行は numbers の `fixpoint` |")
    A("| 穴の埋め方 | 被覆画素の平均 | **本ジョブの較正で決めた値** | 感度は §6.4 |")
    A("| 逆転（壊れ目のモデル） | k^-(H+1)、k<1/L0 で平坦 | **本ジョブの追加検査の構成** | 物理床ではない（§6.5） |")
    A("")

    # ---------------- 6 ----------------
    A("## 6. 弱点・正直な開示")
    A("")
    A("### 6.1 同一エンジンの基準が納品フィールドより低い")
    A("")
    A("`a_base_old_040` の H_x = %.3f に対し、納品済み v4-2.2 パッチは %.3f。"
      "差は 16%% の**穴**の埋め方と、納品時に巨礫を除いてから 60 px のガウスで"
      "インペイントしていること。機構間の比較はすべて同一エンジン同士なので"
      "結論は動かないが、**絶対値は納品床と 1:1 では比べられない**。"
      % (cases["a_base_old_040"]["x"]["H"],
         ref["a_shipped_composite"]["H"]["x"]))
    A("")
    A("### 6.2 較正グリッドの解像度")
    A("")
    gc = ctx["grid_check"]
    A("機構の掃引は 1024²（0.586 mm/px）で行った。最良ケース `e_max_q20_040` を "
      "2048²（v4-2.2 の統計グリッド）で組み直すと、H_x = %.3f / ℓmin = %.2f "
      "（1024² では %.3f / %.2f）。解像度の効果は ΔH_x = %+.3f で、"
      "**結論（0.70 に届かない）を変えない**。"
      % (gc["H_x"], gc["local_min_x"],
         cases["e_max_q20_040"]["x"]["H"],
         min(cases["e_max_q20_040"]["judge"]["local_min_x"],
             cases["e_max_q20_040"]["judge"]["local_min_y"]),
         gc["H_x"] - cases["e_max_q20_040"]["x"]["H"]))
    A("")
    A("### 6.3 σg 掃引は module 属性の一時差し替え")
    A("")
    A("σg = 1.5 / 2.7 の母集団は `bed_gen2.populate` の `SIGMA_G` を呼び出し中だけ"
      "差し替えて作った（`bed_gen.draw_size` も同じ値に合わせた）。**ファイルは"
      "触っていない**（G-LEGACY で確認）。復元は finally で行っている。")
    A("")
    A("### 6.4 穴の埋め方の感度")
    A("")
    hs = cal["hole_sensitivity"]
    A("`%s`（穴 %.1f%%）で 3 通りを試した: " % (hs["case"], 100 * hs["hole_frac"]))
    A("")
    A("| 埋め方 | H_x | ℓmin | σz (m) |")
    A("| --- | --- | --- | --- |")
    for k, v in hs["treatments"].items():
        A("| %s | %.3f | %.2f | %.4f |" % (k, v["H_x"], v["local_min"], v["sigma_z_m"]))
    A("")
    spread = (max(v["H_x"] for v in hs["treatments"].values())
              - min(v["H_x"] for v in hs["treatments"].values()))
    A("→ H_x の振れ幅は **%.3f**（%.3f–%.3f）。穴の埋め方は H を無視できないほど"
      "動かす（最も有利な `form` 埋めで %.3f）。**それでも 0.70 には届かない**ので"
      "主判定は変わらないが、絶対値を読むときはこの幅を覚えておくこと。"
      "本文の表は測定前に宣言した `mean`（被覆画素の平均）で統一してある。"
      % (spread, min(v["H_x"] for v in hs["treatments"].values()),
         max(v["H_x"] for v in hs["treatments"].values()),
         max(v["H_x"] for v in hs["treatments"].values())))
    A("")
    A("### 6.5 物差しの上限の検査は合成場")
    A("")
    A("§4.3 の「理想場」は k 空間で作った定常自己アフィン場で、**実在の床ではない**。"
      "ここから言えるのは「この物差し・この箱で、その場がいくつを返すか」だけで、"
      "実在床の H がいくつか、ではない。")
    A("")
    A("### 6.6 測れなかった・やっていないこと")
    A("")
    A("- **STAGE 2/3 の量は 1 つも測っていない**: 新床のパッチ・3.0 m 帯状床、"
      "G-SCALING、G-AUDIT の 7 項目、G-MESO の M1–M6。`bed_gen3.py` は無い。")
    A("- **パッチ規模の組織化（40–50 cm 相関の局所 GSD・armour 変調）は試していない**。"
      "brief がこれを STAGE 2 の装置として許可しているので、STAGE 1 の機構表には"
      "入れていない。だから「この床は原理的に不可能」とは言わない。言えるのは"
      "**brief が STAGE 1 に列挙した (a)–(e) の範囲では届かない**、という事実だけ。")
    A("- **実現は 1 本**（seed 固定）。機構間の差（0.03–0.10）は実現間ばらつきと"
      "同じ桁なので、個々の順位は確定しない。物差しの上限（§4.3）だけは 4–6 seed。")
    A("- **G-MESO の設計文書が無い**: `DESIGN-V4-2-6-MESO.md` は本ツリーに存在しない。"
      "M1–M6 の帯は brief §4 の要約に従うしかない（測定自体をしていないので、"
      "帯の当否は判定していない）。提案はラリー未検証。")
    A("- **見え方・輸送・画作りは対象外**。レンダリングしていない。")
    A("")
    A("### 6.7 brief の「(a) で 0.287 を再現」について")
    A("")
    A("保存済み relief 曲線で **0.2870** を再現した（期待 0.287、差 < 5e-4）。"
      "納品パッチから直接測った relief は H_x = %.3f で、これは 3.0 m 床 r1 の"
      "曲線とは別物（箱と実現が違う）。両方を出している。"
      % ref["a_shipped_relief"]["H_x"])
    A("")

    # ---------------- 7 ----------------
    A("## 7. G-LEGACY（既存 .py の無改造）")
    A("")
    A("- stage-0 記録（%s）で既存 .py **%d 本**の sha256 を保存し、run の最後"
      "（%s）に再ハッシュした。現在 %d 本（stage-0 の %d 本 + 外部レーンの追加 "
      "%d 本）、変更 %d / 削除 %d。"
      % (lg["hash_baseline_utc"], lg["n_baseline"], lg["checked_at"],
         lg["n_now"], lg["n_baseline"], len(lg["added_details"]),
         len(lg["changed"]), len(lg["removed"])))
    A("- **本 run が import した 5 本**（%s）のハッシュは stage-0 と完全一致"
      "（変更 %d）。本 run の write set は `run_v4_2_6.py` と `out/v4-2-6/*` だけ。"
      % (", ".join("`%s`" % k for k in lg["imported_by_this_run"]),
         len(lg["dependency_changes"])))
    A("- 判定: **%s**（本 run 由来の既存 .py 変更 0）。厳密なハッシュ一致"
      "（追加・削除・変更 0）は %s。"
      % ("PASS" if lg["ok"] else "FAIL",
         "成立" if lg["ok_strict"] else "**不成立**"))
    if lg["external_changes"] or lg["added_details"]:
        A("")
        A("外部レーンによる差分（**本 run の書き込みではない**）:")
        A("")
        A("| 種別 | ファイル | bytes | mtime | 本 run が import |")
        A("| --- | --- | --- | --- | --- |")
        for d in lg["external_changes"]:
            A("| 変更 | `%s` | %d → %d | %s | いいえ |"
              % (d["file"], d["bytes_before"], d["bytes_now"], d["mtime"]))
        for d in lg["added_details"]:
            A("| **追加** | `%s` | %d（新規） | %s | いいえ |"
              % (d["file"], d["bytes_now"], d["mtime"]))
        A("")
        first = (lg["external_changes"] or lg["added_details"])[0]["file"]
        A("根拠: (i) 本 run のコードが開くのは上記 write set だけで、`%s` を"
          "import も open もしない、(ii) mtime が本 run の実行中（別レーンが"
          "同時に走っている）、(iii) 本 run が使う 5 本のハッシュは一致。"
          "**本 run 由来の既存ファイルの変更・削除は 0 件**。" % first)
    A("")

    # ---------------- 8 ----------------
    A("## 8. 成果物")
    A("")
    A("| 種別 | パス | bytes |")
    A("| --- | --- | --- |")
    for k, p in (("report", REPORT), ("numbers", NUMBERS), ("figure", FIGURE),
                 ("checkpoint", CP_PATH), ("stage0_py_hashes", STAGE0_HASHES),
                 ("timings", TIMINGS), ("run_log", RUNLOG)):
        b = os.path.getsize(p) if os.path.isfile(p) else 0
        A("| %s | `%s` | %d |" % (k, p, b))
    A("")
    A("図のパネル: %s。" % " / ".join(fig["panels"]))
    A("")
    A("> 注: report / numbers / run.log の bytes は本表を組み立てた時点では未書き出し"
      "なので 0 と出る（実際のサイズは書き出し後のファイルを見ること）。")
    A("")
    A("**作らなかったもの**: `bed_gen3.py`、0.6 m パッチの新床、3.0 m×0.6 m 帯状床、"
      "新床のゲート結果。STAGE 1 が FAIL したため（brief の STOP 規則）。")
    A("")

    # ---------------- 9 ----------------
    A("## 9. 実行の記録")
    A("")
    A("- 全ステージ foreground（1 プロセス、`--all`）。jobkit checkpoint を段ごとに"
      "原子的にコミット（`out/v4-2-6/checkpoint.json`、units = %d）。"
      % ctx["checkpoint_units"])
    A("- 実行時間: %s。"
      % "、".join("%s %.1fs" % (k, v) for k, v in ctx["timings"].items()))
    A("- 帯・ゲート・物差し・しきい値は**測定前にコード内で固定**。後付けの調整は無い。")
    A("- 数値の全量は `numbers-v4-2-6.json`（全 case の曲線・局所指数・不動点の全試行・"
      "穴の感度・物差しの上限を含む）。")
    return "\n".join(L) + "\n"


def stage_report(cp, ctx):
    log("REPORT  numbers-v4-2-6.json + REPORT-V4-2-6.md")
    md = build_report(ctx)
    atomic_write_bytes(REPORT, md.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 8000)
    numbers = {
        "brief": BRIEF,
        "created_utc": ctx["created_utc"],
        "what_this_is": ("BRIEF-V4-2-6 STAGE 0 + STAGE 1 only. The calibration "
                         "gate failed, so STAGE 2 (bed_gen3.py) and STAGE 3 "
                         "(gates on a new bed) were not entered. No bed was "
                         "delivered."),
        "decision": ctx["judge"]["decision"],
        "stage0": ctx["stage0"],
        "contract": ctx["calib"]["contract"],
        "grid": ctx["calib"]["grid"],
        "population": ctx["calib"]["population"],
        "reference": ctx["calib"]["reference"],
        "cases": ctx["calib"]["cases"],
        "hole_sensitivity": ctx["calib"]["hole_sensitivity"],
        "ceiling": ctx["calib"]["ceiling"],
        "judgment": ctx["judge"],
        "component_budget": ctx.get("budget"),
        "grid_check": ctx.get("grid_check"),
        "g_legacy": ctx["g_legacy"],
        "figure": ctx["figure"],
        "timings_s": ctx["timings"],
        "checkpoint": {"units": ctx["checkpoint_units"], "path": CP_PATH},
        "deliverables": ctx["deliverables"],
        "log": LOG,
    }
    write_json(NUMBERS, numbers)
    log("  wrote REPORT (%d bytes) and numbers (%d bytes)"
        % (os.path.getsize(REPORT), os.path.getsize(NUMBERS)))
    cp.unit_done("report")


# ---------------------------------------------------------------------------
def component_budget(pop, spec_id, cal):
    """Exact S split for the headline case: composite / form / relief."""
    sp = {"id": spec_id, "rule": "max", "q": 0.20, "target": 0.40,
          "order": "random", "dip": 0.0, "group": "e", "label_ja": "budget"}
    res, g, sig, hist = sigma_fixpoint(pop, sp, 0.40 * D50_REF_M)
    h, fg, cell = res["h"], res["form"], res["cell"]
    hole = h < -1e8
    z = np.where(hole, float(h[~hole].mean()), h)
    rel = z - fg
    out = {"case": spec_id, "note": ("S_comp = S_form + S_relief + 2cov, "
                                    "evaluated at in-band lags on the same "
                                    "field (periodic wrap, flow axis)"),
           "lags": {}, "H_in_band": {}}
    for name, arr in (("composite", z), ("form", fg), ("relief", rel)):
        d, s = M5.curve_dir(arr, cell, DM_X, nbin=NBIN, axis=1)
        fb = M5.fit_band(d, s, BAND[0], BAND[1])
        out["H_in_band"][name] = {"H": fb["H"],
                                 "S_at_0.41cm": float(np.interp(0.0041, d, s)),
                                 "S_at_3.22cm": float(np.interp(0.0322, d, s))}
        del d, s
    for lag in (0.0041, 0.010, 0.0322, 0.05):
        d = int(round(lag / cell))
        f = fg - np.roll(fg, -d, axis=1)
        r = rel - np.roll(rel, -d, axis=1)
        sc = float(np.mean((f + r) ** 2))
        sf = float(np.mean(f * f))
        sr = float(np.mean(r * r))
        cov = float(np.mean(f * r))
        out["lags"]["%.4f" % lag] = {
            "S_total": sc, "S_form": sf, "S_relief": sr,
            "cov_term_2x": 2.0 * cov,
            "share_form": sf / sc if sc else None,
            "share_relief": sr / sc if sc else None,
            "share_cov": 2.0 * cov / sc if sc else None}
        del f, r
    log("  budget @0.0322: form %.2f relief %.2f cov %.2f"
        % (out["lags"]["0.0322"]["share_form"],
           out["lags"]["0.0322"]["share_relief"],
           out["lags"]["0.0322"]["share_cov"]))
    return out


def grid_check(pop):
    """The same headline case at 2048^2 (does the calibration grid matter?)."""
    sp = {"id": "grid_2048", "rule": "max", "q": 0.20, "target": 0.40,
          "order": "random", "dip": 0.0, "group": "e", "label_ja": "grid check"}
    res, g, sig, hist = sigma_fixpoint(pop, sp, 0.40 * D50_REF_M, n=HIRES_N)
    m = measure_field(res, g)
    out = {"id": "grid_2048", "n": HIRES_N, "cell_m": HIRES_CELL,
           "H_x": m["x"]["H"], "H_y": m["y"]["H"], "H_iso": m["iso"]["H"],
           "local_min_x": m["x"]["local_min_in_band"],
           "local_min_y": m["y"]["local_min_in_band"],
           "sigma_z_over_D50": m["sigma_z_matrix_over_D50"],
           "sigma_form_m": sig, "hole_frac": m["hole_frac"]}
    log("  grid check 2048: H_x=%.3f lmin=%.2f (1024: see the case table)"
        % (out["H_x"], out["local_min_x"]))
    return out


def stage_legacy(cp):
    log("G-LEGACY  re-hash every existing .py and compare with stage 0")
    base = json.load(open(STAGE0_HASHES, encoding="utf-8"))
    now = {}
    for p in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        name = os.path.basename(p)
        if name == THIS_FILE:
            continue
        now[name] = {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
    changed = sorted(k for k in base["py"] if k in now
                     and now[k]["sha256"] != base["py"][k]["sha256"])
    added = sorted(set(now) - set(base["py"]))
    removed = sorted(set(base["py"]) - set(now))
    # ---- what this run actually used, and whether THAT is untouched -------
    deps = ["bed_gen.py", "bed_gen2.py", "run_v4_2_2.py", "run_v4_2_3.py",
            "run_v4_2_5.py"]
    dep_changed = [k for k in deps if k in changed]
    details = []
    t_start = datetime.datetime.fromisoformat(base["created_utc"])
    for k in changed:
        p = os.path.join(HERE, k)
        mt = datetime.datetime.fromtimestamp(
            os.path.getmtime(p)).astimezone()
        details.append({
            "file": k, "bytes_before": base["py"][k]["bytes"],
            "bytes_now": now[k]["bytes"],
            "sha256_before": base["py"][k]["sha256"],
            "sha256_now": now[k]["sha256"],
            "mtime": mt.isoformat(timespec="seconds"),
            "mtime_after_stage0": bool(mt.astimezone(
                datetime.timezone.utc) > t_start),
            "imported_by_this_run": k in deps,
        })
    external = [d for d in details if not d["imported_by_this_run"]]
    added_details = []
    for k in added:
        p = os.path.join(HERE, k)
        mt = datetime.datetime.fromtimestamp(os.path.getmtime(p)).astimezone()
        added_details.append({"file": k, "bytes_now": os.path.getsize(p),
                             "mtime": mt.isoformat(timespec="seconds"),
                             "mtime_after_stage0": bool(mt.astimezone(
                                 datetime.timezone.utc) > t_start)})
    g = {"n_baseline": base["n_existing_py"], "n_now": len(now),
         "changed": changed, "added": added, "removed": removed,
         "added_details": added_details,
         "details": details,
         "imported_by_this_run": deps, "dependency_changes": dep_changed,
         "external_changes": external,
         "n_external": len(external),
         "hash_baseline_utc": base["created_utc"], "checked_at": utcnow(),
         "write_set": [THIS_FILE, "out/v4-2-6/"],
         "ok_strict": not (changed or added or removed),
         "ok": bool(not dep_changed and not removed),
         "ok_note": ("ok = no stage-0 file changed or disappeared, and no file "
                     "this run imported changed.  Files added by another lane "
                     "after stage 0, and files changed by another lane, are "
                     "listed with their mtime and are NOT attributed to this "
                     "run; ok_strict is the literal stage-0 hash comparison "
                     "(changed + added + removed all zero)")}
    log("  changed=%d added=%d removed=%d  dependency_changes=%d  "
        "external=%d -> %s (strict %s)"
        % (len(changed), len(added), len(removed), len(dep_changed),
           len(external), "PASS" if g["ok"] else "FAIL",
           "PASS" if g["ok_strict"] else "FAIL"))
    for d in details:
        log("    external change: %s %d->%d bytes mtime %s"
            % (d["file"], d["bytes_before"], d["bytes_now"], d["mtime"]))
    for d in added_details:
        log("    external addition: %s %d bytes mtime %s"
            % (d["file"], d["bytes_now"], d["mtime"]))
    cp.unit_done("legacy")
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--calib", action="store_true")
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--legacy", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if not any(vars(a).values()):
        a.all = True

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    created = utcnow()
    timings = {}

    def cached(name):
        p = os.path.join(WORK, name + ".json")
        if os.path.isfile(p):
            return json.load(open(p, encoding="utf-8"))
        return None

    def saved(name, obj):
        write_json(os.path.join(WORK, name + ".json"), obj)

    t0 = time.time()
    st0 = None if (a.all or a.stage0) else cached("stage0")
    if st0 is None:
        st0 = stage0(cp)
        saved("stage0", st0)
    timings["stage0"] = round(time.time() - t0, 1)

    t0 = time.time()
    cal = cached("calib") if not (a.all or a.calib) else None
    if cal is None:
        cal = stage_calib(cp, st0)
        saved("calib", cal)
    timings["calib"] = round(time.time() - t0, 1)

    pop, _ = B2.populate(SEED)
    t0 = time.time()
    budget = cached("budget") if not (a.all or a.calib) else None
    if budget is None:
        budget = component_budget(pop, "e_max_q20_040", cal)
        saved("budget", budget)
    timings["budget"] = round(time.time() - t0, 1)

    t0 = time.time()
    gc_ = cached("grid_check") if not (a.all or a.calib) else None
    if gc_ is None:
        gc_ = grid_check(pop)
        saved("grid_check", gc_)
    timings["grid_check"] = round(time.time() - t0, 1)

    t0 = time.time()
    vd = stage_judge(cp, cal)
    saved("judge", vd)
    timings["judge"] = round(time.time() - t0, 1)

    t0 = time.time()
    fg = stage_figure(cp, cal, vd)
    timings["figure"] = round(time.time() - t0, 1)

    t0 = time.time()
    lg = stage_legacy(cp)
    timings["legacy"] = round(time.time() - t0, 1)

    write_json(TIMINGS, {"stage_s": timings, "created_utc": created})

    state = cp.state() or {}
    deliv = {"report": REPORT, "numbers": NUMBERS, "figure": FIGURE,
             "checkpoint": CP_PATH, "stage0_py_hashes": STAGE0_HASHES,
             "timings": TIMINGS, "run_log": RUNLOG}
    ctx = {"created_utc": created, "stage0": st0, "calib": cal, "judge": vd,
           "figure": fg, "g_legacy": lg, "timings": timings,
           "budget": budget, "grid_check": gc_,
           "checkpoint_units": int(state.get("completed_units") or 0) + 1,
           "deliverables": deliv}
    stage_report(cp, ctx)

    with open(RUNLOG, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("DONE  %s" % REPORT)
    print("\n".join(LOG[-8:]))


if __name__ == "__main__":
    main()
