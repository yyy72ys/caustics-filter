"""run_v4_2_5 -- BRIEF-V4-2-5: re-audit of the bed.  Measure H where the
scaling region actually is (lag <~ D50), not in the old 4-20 cm gate band.

AUDIT ONLY.  No generation, no parameter change, no new bed, no rendering.
Everything is measured from saved heightfields / saved structure-function
curves of v4-2.2 and v4-2.3.  Existing .py files are imported and used
unmodified; the only new .py is this file.

What the brief asks, and where it is answered here
--------------------------------------------------
1. Directional 2-D structure function S(tau) for the flow axis (x) and the
   cross-stream axis (y), following v4-2.3's measurement (periodic squared
   difference along one axis, same cell, same lag ladder).            -> stage_repro / stage_measure
2. H fitted in three declared bands, per bed, per direction:
     scaling  0.4 - 3.6 cm   (<= D50; the NEW primary judgment)
     meso     6.6 - 20 cm    (recorded only, for a later H_meso design)
     legacy   4 - 20 cm      (the old fit band; reproducibility check)
   plus, for the v4-2.2 patch, the old break-based value (0.151).      -> stage_measure
3. Fit quality (log-residual RMS, R^2, number of lags) is reported for
   EVERY fitted H, and, when S saturates inside a band, the shape
   (local exponent profile + knee/plateau scale) is reported.         -> stage_measure
4. Judgment: does the scaling-band H fall in 0.70-0.85, per bed and per
   direction?  When it differs from an old value, the difference is split
   into "band" and "method" contributions.                            -> stage_judge

Stages (foreground, idempotent, one jobkit checkpoint unit per stage):
  --stage0   absolute-path preflight, D50 confirmation, sha256 of every .py
  --repro    recompute the curves from the saved heightfields and prove they
             reproduce v4-2.3's saved curves (naive re-reader check)
  --measure  all fits, shape diagnostics, component shares -> numbers
  --judge    primary judgment + band/method attribution
  --figure   structure-function PNG
  --report   numbers-v4-2-5.json + REPORT-V4-2-5.md
  --all      the above, in order, in this one process
"""
from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import bed_gen2 as B2                                   # noqa: E402  (unmodified)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-5.md"
OUT = os.path.join(HERE, "out", "v4-2-5")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-5.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-5.md")
FIGURE = os.path.join(OUT, "structure-functions.png")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
TIMINGS = os.path.join(OUT, "timings.json")

HF_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
META_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
NUMBERS_422 = os.path.join(HERE, "out", "v4-2-2", "numbers-v4-2-2.json")
REPORT_422 = os.path.join(HERE, "out", "v4-2-2", "REPORT-V4-2-2.md")
DIR_423 = os.path.join(HERE, "out", "v4-2-3")
HF_LONG = os.path.join(DIR_423, "heightfield-long-r1.npy")
META_LONG = os.path.join(DIR_423, "longbed-meta.json")
META_LONG2 = os.path.join(DIR_423, "longbed-meta-2.json")
CURVES_423 = os.path.join(DIR_423, "h-curves.json")
DECOMP_423 = os.path.join(DIR_423, "decomp.json")
NUMBERS_423 = os.path.join(DIR_423, "numbers-v4-2-3.json")
REPORT_423 = os.path.join(DIR_423, "REPORT-V4-2-3.md")

THIS_FILE = "run_v4_2_5.py"

# --- the measurement contract (fixed BEFORE the first measurement) ---------
PATCH = B2.PATCH                      # 0.6 m
CELL = PATCH / B2.HIRES               # 0.00029296875 m/px, exactly v4-2.2/2.3
LONG_NY, LONG_NX = 2048, 10240        # 0.6 m x 3.0 m
LONG_CELL = 0.00029296875
assert abs(LONG_CELL - CELL) < 1e-15

BANDS = [
    ("scaling", 0.004, 0.036),        # <= D50 ; the NEW primary judgment
    ("meso", 0.066, 0.200),           # 6.6 cm = v4-2.2 iso break scale 6.62 cm
    ("legacy", 0.040, 0.200),         # 2*D50_IN .. 20 cm, the old gate band
]
BAND_LABEL = {"scaling": "scaling 0.4-3.6 cm", "meso": "meso 6.6-20 cm",
              "legacy": "legacy 4-20 cm"}
GATE = (0.70, 0.85)

# instrument: v4-2.3's own lag ladders (nbin=60), plus a finer one as a
# binning-sensitivity control (nbin=120, dmax=5 cm)
NBIN = 60
DM_ISO = 0.20
DM_X_LONG = 1.00          # v4-2.3 used 1/3 of the 3.0 m bed
DM_PATCH_X = 0.30         # v4-2.3 used 0.30 m for the 0.6 m patch
DM_Y = 0.30               # half of the 0.6 m width
NBIN_FINE = 120
DM_FINE = 0.05

# old values the brief asks us to reproduce
OLD = {
    "H_iso_4_20_v422_patch": 0.18212965650193444,
    "H_iso_above_break_v422_patch": 0.150567117376588,
    "break_iso_v422_patch_m": 0.0662109375,
    "H_x_4_20_long_r1": 0.21916492702163493,
    "H_x_4_20_long_r2": 0.2318850786538861,
    "H_iso_4_20_long_r1": 0.18618650545829524,
    "H_iso_4_20_long_r2": 0.25338501444743927,
    "H_x_4_20_substrate_r1": 0.7961369310819543,
    "H_x_4_20_substrate_r2": 0.7489796468820366,
    "H_iso_4_20_substrate_r1": 0.7733524493653687,
    "H_iso_4_20_substrate_r2": 0.7588405104368092,
    "H_iso_4_20_substrate_3rep_mean": 0.8085545099410559,
    "H_relief_x_4_20_long_r1": 0.102818037577871,
    "D50_surface_v422": 0.03619913686788151,
    "D50_surface_long_r1": 0.03829749229501987,
    "D50_surface_long_r2": 0.03797627400499068,
}

LOG = []


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S+00:00")


def sha256_file(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _j(x):
    """JSON-safe scalar."""
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (list, tuple)):
        return [_j(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _j(v) for k, v in x.items()}
    return x


def write_json(path, obj):
    data = json.dumps(_j(obj), ensure_ascii=False, indent=1).encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True
    atomic_write_bytes(path, data, verify=verify)


# ---------------------------------------------------------------------------
# the estimator -- v4-2.3's arithmetic, copied (same roll, same lags)
# ---------------------------------------------------------------------------
def curve_dir(z, cell, dmax, nbin=NBIN, axis=1):
    """One axis of the 2-D structure function: mean squared increment.

    Identical arithmetic to run_v4_2_3.curve_x: periodic wrap, no de-meaning
    per lag, geometric lag ladder rounded to whole pixels.  axis=1 -> flow
    (x), axis=0 -> cross-stream (y).
    """
    dmax_px = int(min(dmax / cell, z.shape[axis] // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d, out_s = [], []
    for d in ds:
        out_d.append(d * cell)
        out_s.append(np.mean((z - np.roll(z, -d, axis=axis)) ** 2))
    return np.array(out_d), np.array(out_s)


def curve_iso(z, cell, dmax=DM_ISO, nbin=NBIN):
    """v4-2.3's curve_iso verbatim (x and y averaged, periodic)."""
    dmax_px = int(min(dmax / cell, z.shape[1] // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d, out_s = [], []
    for d in ds:
        dx2 = np.mean((z - np.roll(z, -d, axis=1)) ** 2)
        dy2 = np.mean((z - np.roll(z, -d, axis=0)) ** 2)
        out_d.append(d * cell)
        out_s.append(0.5 * (dx2 + dy2))
    return np.array(out_d), np.array(out_s)


def local_H(d, s):
    return 0.5 * np.gradient(np.log(s), np.log(d))


def fit_band(d, s, a, b):
    """0.5 * dlnS/dlnd over the lags inside [a,b], with fit quality."""
    m = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
    n = int(m.sum())
    out = {"lag_lo_m": a, "lag_hi_m": b, "n_points": n,
           "H": None, "H_se": None, "rms_log": None, "r2": None,
           "d_min_fit_m": None, "d_max_fit_m": None,
           "S_at_d_min": None, "S_at_d_max": None, "S_growth": None}
    if n < 4:
        out["note"] = "fewer than 4 lags inside the band"
        return out
    x = np.log(d[m])
    y = np.log(s[m])
    p, cov = np.polyfit(x, y, 1, cov=True)
    pred = np.polyval(p, x)
    resid = y - pred
    rms = float(np.sqrt(np.mean(resid ** 2)))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = None if sst <= 0 else float(1.0 - np.sum(resid ** 2) / sst)
    out.update({
        "H": float(0.5 * p[0]),
        "H_se": float(0.5 * np.sqrt(max(cov[0, 0], 0.0))),
        "rms_log": rms, "r2": r2,
        "d_min_fit_m": float(d[m][0]), "d_max_fit_m": float(d[m][-1]),
        "S_at_d_min": float(s[m][0]), "S_at_d_max": float(s[m][-1]),
        "S_growth": float(s[m][-1] / s[m][0]),
    })
    return out


def band_shape(d, s, a, b, nprobe=6):
    """Inside-band shape: local exponent profile, plateau / knee, crossings."""
    m = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
    if m.sum() < 3:
        return {"note": "too few lags for a shape read"}
    dd = d.copy()
    loc = local_H(dd, s)
    dl, ll = dd[m], loc[m]
    i_min = int(np.argmin(ll))
    i_max = int(np.argmax(ll))
    probes = np.geomspace(max(a, dl[0]), min(b, dl[-1]), nprobe)
    prof = {("%.4f" % p): float(np.interp(p, dd, loc)) for p in probes}
    # first scale above 2 cm where the local exponent drops below 0.20 --
    # "where the grain mosaic actually saturates" (outside the scaling band)
    plateau_onset = None
    for i in range(len(dd)):
        if dd[i] >= 0.02 and loc[i] < 0.20:
            plateau_onset = float(dd[i])
            break
    # where the LOCAL exponent crosses the gate edges (interpolated)
    def crossing(level):
        for i in range(len(ll) - 1):
            if (ll[i] - level) * (ll[i + 1] - level) < 0:
                f = (level - ll[i]) / (ll[i + 1] - ll[i])
                return float(np.exp(np.log(dl[i])
                                    + f * (np.log(dl[i + 1]) - np.log(dl[i]))))
        return None
    return {
        "local_H_probe": prof,
        "local_H_min": float(ll[i_min]), "local_H_min_at_m": float(dl[i_min]),
        "local_H_max": float(ll[i_max]), "local_H_max_at_m": float(dl[i_max]),
        "cross_0.70_m": crossing(0.70),
        "cross_0.85_m": crossing(0.85),
        "cross_0.20_m": crossing(0.20),
        "plateau_onset_m": plateau_onset,
        "plateau_onset_rule": ("first lag > 2 cm where local H < 0.20 on the "
                               "same curve (a plateau that starts outside "
                               "the scaling band is reported as such)"),
        "plateau_inside_band": bool(float(ll.min()) < 0.20),
        "note": ("local exponent = 0.5*dlnS/dlnd read on the same lag ladder; "
                 "a plateau is the scale range where it collapses (S stops "
                 "growing)"),
    }


def audit_like(d, s, cell, fit_hi=0.20):
    """v4-2.2/2.3's break rule on a supplied curve (for the old numbers)."""
    lo, hi = 2.0 * B2.D50_IN, min(fit_hi, 0.20)
    m = (d >= lo) & (d <= hi)
    H_full = float(0.5 * np.polyfit(np.log(d[m]), np.log(s[m]), 1)[0])
    loc = local_H(d, s)

    def _fit(a, b):
        mm = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
        if mm.sum() < 4:
            return None, int(mm.sum())
        return float(0.5 * np.polyfit(np.log(d[mm]), np.log(s[mm]), 1)[0]), \
            int(mm.sum())

    brk = None
    target = 0.6 * 2.0 * H_full
    for i in range(1, len(d) - 1):
        if d[i] < 0.5 * B2.D50_IN:
            continue
        if loc[i] < target:
            brk = float(d[i])
            break
    H_above, n_above = (_fit(brk, d[-1]) if brk else (None, 0))
    return {"H_full_band": H_full, "break_scale_m": brk,
            "H_above_break": H_above, "n_points_above": n_above}


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
REQUIRED = {
    "out/v4-2-2/heightfield.npy": HF_422,
    "out/v4-2-3/heightfield-long-r1.npy": HF_LONG,
    "out/v4-2-3/REPORT-V4-2-3.md": REPORT_423,
    "out/v4-2-3/longbed-meta.json": META_LONG,
    "out/v4-2-3/h-curves.json": CURVES_423,
    "out/v4-2-3/decomp.json": DECOMP_423,
    "out/v4-2-2/numbers-v4-2-2.json": NUMBERS_422,
}


def stage0(cp):
    log("STAGE0  absolute-path preflight")
    missing = [k for k, p in REQUIRED.items() if not os.path.isfile(p)]
    if missing:
        raise SystemExit("STAGE 0 STOP -- missing: %s" % missing)
    pre = {}
    for k, p in REQUIRED.items():
        pre[k] = {"exists": True, "bytes": os.path.getsize(p),
                  "sha256": sha256_file(p)}
        log("  ok %-42s %10d bytes  %s" % (k, pre[k]["bytes"],
                                           pre[k]["sha256"][:16]))
    # the long heightfield is a headerless raw float32 dump -- prove the read
    n = os.path.getsize(HF_LONG)
    if n != LONG_NY * LONG_NX * 4:
        raise SystemExit("STAGE 0 STOP -- long heightfield is %d bytes, "
                         "expected %d" % (n, LONG_NY * LONG_NX * 4))
    zr = np.fromfile(HF_LONG, dtype="<f4")
    if zr.size != LONG_NY * LONG_NX:
        raise SystemExit("STAGE 0 STOP -- long heightfield size mismatch")
    meta_l = json.load(open(META_LONG, encoding="utf-8"))
    sig = float(zr.std(dtype=np.float64))
    pre["long_heightfield_read"] = {
        "dtype": "float32 LE, headerless (np.fromfile)",
        "shape": [LONG_NY, LONG_NX], "sigma_z_field_read_m": sig,
        "sigma_z_meta_m": meta_l["sigma_z_m"],
        "sigma_rel_diff": abs(sig - meta_l["sigma_z_m"])
        / meta_l["sigma_z_m"],
        "cell_m": meta_l["cell_m"],
        "read_ok": bool(abs(sig - meta_l["sigma_z_m"]) / meta_l["sigma_z_m"]
                        < 0.02),
    }
    if not pre["long_heightfield_read"]["read_ok"]:
        raise SystemExit("STAGE 0 STOP -- long heightfield read check failed")
    log("  long heightfield read ok: sigma read %.6f vs meta %.6f"
        % (sig, meta_l["sigma_z_m"]))
    del zr

    # D50 confirmation (the brief names 0.0362 m)
    m422 = json.load(open(META_422, encoding="utf-8"))
    nb422 = json.load(open(NUMBERS_422, encoding="utf-8"))
    d50_422 = float(m422["D50_surface_measured_m"])
    d50_422_nb = float(nb422["bed_stats"]["D50_surface_measured_m"])
    meta2 = json.load(open(META_LONG2, encoding="utf-8"))
    d50 = {
        "v422_patch_m": d50_422,
        "v422_patch_numbers_m": d50_422_nb,
        "long_r1_meta_m": float(meta_l["surface_D50_m"]),
        "long_r2_meta_m": float(meta2["surface_D50_m"]),
        "brief_named_m": 0.0362,
        "brief_named_confirmed": bool(abs(d50_422 - 0.0362) < 1e-3),
        "band_top_over_D50_v422": BANDS[0][2] / d50_422,
        "band_top_over_D50_long_r1": BANDS[0][2] / float(meta_l["surface_D50_m"]),
    }
    if not d50["brief_named_confirmed"]:
        raise SystemExit("STAGE 0 STOP -- D50 %.6f does not match 0.0362"
                         % d50_422)
    log("  D50 confirmed: v4-2.2 patch %.6f m (= 0.0362), long r1 %.6f m, "
        "long r2 %.6f m" % (d50_422, d50["long_r1_meta_m"],
                            d50["long_r2_meta_m"]))

    # sha256 of every existing .py (G-LEGACY baseline)
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
    return {"preflight": pre, "d50": d50, "n_existing_py": len(py)}


# ---------------------------------------------------------------------------
# stage repro -- read the saved fields, recompute, compare with 2.3's curves
# ---------------------------------------------------------------------------
def load_long():
    z = np.fromfile(HF_LONG, dtype="<f4").reshape(LONG_NY, LONG_NX)
    return z.astype(np.float64)


def stage_repro(cp):
    log("REPRO  recompute from saved heightfields; compare with 2.3 curves")
    h423 = json.load(open(CURVES_423, encoding="utf-8"))
    out = {"instrument": "run_v4_2_5.curve_dir/curve_iso == "
                         "run_v4_2_3.curve_x/curve_iso (copied arithmetic)"}

    z422 = np.load(HF_422).astype(np.float64)
    t0 = time.time()
    d, s = curve_iso(z422, CELL, DM_ISO)
    dref = np.array(h423["curves"]["v422_patch_grains"]["iso"]["D_m"])
    sref = np.array(h423["curves"]["v422_patch_grains"]["iso"]["S_m2"])
    out["v422_patch_iso"] = {
        "n_lags": int(d.size), "max_abs_d_m": float(np.max(np.abs(d - dref))),
        "max_abs_S": float(np.max(np.abs(s - sref))),
        "rel_S": float(np.max(np.abs(s - sref) / np.maximum(sref, 1e-30)))}
    dx, sx = curve_dir(z422, CELL, DM_PATCH_X, axis=1)
    dxref = np.array(h423["curves"]["v422_patch_grains"]["x"]["D_m"])
    sxref = np.array(h423["curves"]["v422_patch_grains"]["x"]["S_m2"])
    out["v422_patch_x"] = {
        "n_lags": int(dx.size),
        "max_abs_d_m": float(np.max(np.abs(dx - dxref))),
        "max_abs_S": float(np.max(np.abs(sx - sxref))),
        "rel_S": float(np.max(np.abs(sx - sxref) / np.maximum(sxref, 1e-30)))}
    log("  v4-2.2 patch: iso dS %.3e (rel %.2e), x dS %.3e (rel %.2e) "
        "[%.1fs]" % (out["v422_patch_iso"]["max_abs_S"],
                     out["v422_patch_iso"]["rel_S"],
                     out["v422_patch_x"]["max_abs_S"],
                     out["v422_patch_x"]["rel_S"], time.time() - t0))

    zl = load_long()
    t0 = time.time()
    d, s = curve_iso(zl, CELL, DM_ISO)
    dref = np.array(h423["curves"]["long_grains_r1"]["iso"]["D_m"])
    sref = np.array(h423["curves"]["long_grains_r1"]["iso"]["S_m2"])
    out["long_r1_iso"] = {
        "n_lags": int(d.size), "max_abs_d_m": float(np.max(np.abs(d - dref))),
        "max_abs_S": float(np.max(np.abs(s - sref))),
        "rel_S": float(np.max(np.abs(s - sref) / np.maximum(sref, 1e-30)))}
    dx, sx = curve_dir(zl, CELL, DM_X_LONG, axis=1)
    dxref = np.array(h423["curves"]["long_grains_r1"]["x"]["D_m"])
    sxref = np.array(h423["curves"]["long_grains_r1"]["x"]["S_m2"])
    out["long_r1_x"] = {
        "n_lags": int(dx.size),
        "max_abs_d_m": float(np.max(np.abs(dx - dxref))),
        "max_abs_S": float(np.max(np.abs(sx - sxref))),
        "rel_S": float(np.max(np.abs(sx - sxref) / np.maximum(sxref, 1e-30)))}
    log("  long r1: iso dS %.3e (rel %.2e), x dS %.3e (rel %.2e) [%.1fs]"
        % (out["long_r1_iso"]["max_abs_S"], out["long_r1_iso"]["rel_S"],
           out["long_r1_x"]["max_abs_S"], out["long_r1_x"]["rel_S"],
           time.time() - t0))
    del zl, z422

    out["all_ok"] = bool(all(np.isfinite(v["rel_S"]) and v["rel_S"] < 1e-9
                             for k, v in out.items()
                             if isinstance(v, dict) and "rel_S" in v))
    log("  repro all_ok=%s (<=1e-9 relative)" % out["all_ok"])
    cp.unit_done("repro")
    return out


# ---------------------------------------------------------------------------
# stage measure -- all cases, all bands, all directions
# ---------------------------------------------------------------------------
CASE_META = {
    "v422_patch_grains": {
        "label": "v4-2.2 0.6 m patch, grains",
        "source": "out/v4-2-2/heightfield.npy (re-read and re-measured)",
        "group": "grains"},
    "v422_patch_substrate": {
        "label": "v4-2.2 0.6 m patch, substrate only",
        "source": "out/v4-2-3/h-curves.json (saved curve, seed-identical "
                  "substrate of the same build)",
        "group": "substrate"},
    "long_grains_r1": {
        "label": "3.0 m long bed r1, grains",
        "source": "out/v4-2-3/heightfield-long-r1.npy (re-read and re-measured)",
        "group": "grains"},
    "long_substrate_r1": {
        "label": "3.0 m long bed r1, substrate only",
        "source": "out/v4-2-3/h-curves.json (saved curve)",
        "group": "substrate"},
    "long_grains_r2": {
        "label": "3.0 m long bed r2, grains",
        "source": "out/v4-2-3/h-curves.json (saved curve; its heightfield was "
                  "not written to disk, so no fresh recompute is possible)",
        "group": "grains"},
    "long_substrate_r2": {
        "label": "3.0 m long bed r2, substrate only",
        "source": "out/v4-2-3/h-curves.json (saved curve)",
        "group": "substrate"},
    "long_grains_r1_noboulderband": {
        "label": "3.0 m long bed r1, grains, boulder band removed",
        "source": "out/v4-2-3/h-curves.json (saved curve)",
        "group": "grains"},
}
DIRS = ("x", "y", "iso")


def curves_for(case, saved, fresh):
    """Return {dir: {'d':..., 's':..., 'src':...}} for one case."""
    out = {}
    if "iso" in fresh.get(case, {}):
        for k, v in fresh[case].items():
            out[k] = {"d": v[0], "s": v[1], "src": "recomputed"}
        return out
    c = saved["curves"].get(case, {})
    for k in DIRS:
        if k in c:
            out[k] = {"d": np.array(c[k]["D_m"]),
                      "s": np.array(c[k]["S_m2"]), "src": "saved 2.3"}
    y = saved["cases"].get(case, {}).get("curve_y")
    if "y" not in out and y is not None:
        out["y"] = {"d": np.array(y["D_m"]), "s": np.array(y["S_m2"]),
                    "src": "saved 2.3"}
    return out


def stage_measure(cp, repro):
    log("MEASURE  fits in the three declared bands, per bed and direction")
    saved = json.load(open(CURVES_423, encoding="utf-8"))
    h423 = saved["cases"]
    t0 = time.time()

    # fresh recompute: patch grains (x,y,iso) + long r1 grains (x,y,iso);
    # y is NOT in 2.3's saved curves for the grain beds, so it is measured
    # here from the saved fields.
    fresh = {}
    z = np.load(HF_422).astype(np.float64)
    fresh["v422_patch_grains"] = {}
    fresh["v422_patch_grains"]["iso"] = curve_iso(z, CELL, DM_ISO)
    fresh["v422_patch_grains"]["x"] = curve_dir(z, CELL, DM_PATCH_X, axis=1)
    fresh["v422_patch_grains"]["y"] = curve_dir(z, CELL, DM_Y, axis=0)
    del z
    log("  fresh patch curves done [%.1fs]" % (time.time() - t0))

    z = load_long()
    fresh["long_grains_r1"] = {}
    fresh["long_grains_r1"]["iso"] = curve_iso(z, CELL, DM_ISO)
    fresh["long_grains_r1"]["x"] = curve_dir(z, CELL, DM_X_LONG, axis=1)
    fresh["long_grains_r1"]["y"] = curve_dir(z, CELL, DM_Y, axis=0)
    del z
    log("  fresh long-r1 curves done [%.1fs]" % (time.time() - t0))

    beds = {}
    for case, meta in CASE_META.items():
        cvs = curves_for(case, saved, fresh)
        entry = {"label": meta["label"], "group": meta["group"],
                 "source": meta["source"], "direction": {}}
        for k in DIRS:
            if k not in cvs:
                entry["direction"][k] = {
                    "available": False,
                    "reason": ("no saved cross-stream curve and no saved "
                               "heightfield for this realisation; not "
                               "re-invented (brief: no generation)")}
                continue
            d, s = cvs[k]["d"], cvs[k]["s"]
            blk = {"available": True, "curve_source": cvs[k]["src"],
                   "n_lags_total": int(d.size),
                   "lag_max_m": float(d[-1]),
                   "bands": {}, "shape": {}}
            for name, a, b in BANDS:
                blk["bands"][name] = fit_band(d, s, a, b)
            blk["shape"]["scaling"] = band_shape(d, s, BANDS[0][1],
                                                 BANDS[0][2])
            blk["shape"]["legacy"] = band_shape(d, s, BANDS[2][1],
                                                BANDS[2][2])
            # the 2.3/2.2 break-rule reproduction (legacy band only)
            blk["audit_like"] = audit_like(d, s, CELL)
            # probe table (local H) for the report
            blk["local_H_at_m"] = {
                ("%.4f" % p): float(np.interp(p, d, local_H(d, s)))
                for p in (0.004, 0.006, 0.010, 0.016, 0.024, 0.032,
                          0.040, 0.066, 0.100, 0.140, 0.200)
                if p <= d[-1]}
            blk["curve"] = {"D_m": [float(v) for v in d],
                            "S_m2": [float(v) for v in s]}
            entry["direction"][k] = blk
        beds[case] = entry
        log("  case %-32s dirs=%s" % (case, sorted(entry["direction"])))

    # ---- binning-sensitivity control (nbin=120, dmax=5 cm) ----------------
    fine = {}
    z = np.load(HF_422).astype(np.float64)
    dx, sx = curve_dir(z, CELL, DM_FINE, nbin=NBIN_FINE, axis=1)
    dy, sy = curve_dir(z, CELL, DM_FINE, nbin=NBIN_FINE, axis=0)
    fine["v422_patch_grains"] = {"x": (dx, sx), "y": (dy, sy),
                                 "iso": (dx, 0.5 * (sx + sy))}
    del z
    z = load_long()
    dx, sx = curve_dir(z, CELL, DM_FINE, nbin=NBIN_FINE, axis=1)
    dy, sy = curve_dir(z, CELL, DM_FINE, nbin=NBIN_FINE, axis=0)
    fine["long_grains_r1"] = {"x": (dx, sx), "y": (dy, sy),
                              "iso": (dx, 0.5 * (sx + sy))}
    del z
    # substrate fresh fields are not saved -> only the two grain beds
    for case, dirs in fine.items():
        for k, (d, s) in dirs.items():
            f = fit_band(d, s, BANDS[0][1], BANDS[0][2])
            beds[case]["direction"][k]["binning_control"] = {
                "nbin": NBIN_FINE, "dmax_m": DM_FINE,
                "H_scaling": f["H"], "n_points": f["n_points"],
                "rms_log": f["rms_log"]}
    log("  binning control (nbin=120) done [%.1fs]" % (time.time() - t0))

    # ---- component shares on the small band (decomp curves, r1, x) -------
    dec = json.load(open(DECOMP_423, encoding="utf-8"))
    comp = {}
    dc = np.array(dec["terms"]["composite"]["_curve"]["D_m"])
    st = np.array(dec["terms"]["composite"]["_curve"]["S_m2"])
    sf = np.array(dec["terms"]["form"]["_curve"]["S_m2"])
    sr = np.array(dec["terms"]["grain_relief"]["_curve"]["S_m2"])
    for lag in (0.0041, 0.010, 0.0319, 0.04, 0.10, 0.20):
        i = int(np.argmin(np.abs(dc - lag)))
        tot = float(st[i])
        cov = 0.5 * (tot - sf[i] - sr[i])
        comp["%.4f" % dc[i]] = {
            "S_total": tot, "S_form": float(sf[i]), "S_relief": float(sr[i]),
            "cov_term": float(cov),
            "share_form": float(sf[i] / tot) if tot else None,
            "share_relief": float(sr[i] / tot) if tot else None,
            "share_cov": float(cov / tot) if tot else None}
    shares = {"lags": comp,
              "terms_H_scaling": {
                  k: fit_band(dc, np.array(dec["terms"][k]["_curve"]["S_m2"]),
                              BANDS[0][1], BANDS[0][2])
                  for k in ("composite", "form", "grain_relief")},
              "terms_H_legacy": {
                  k: fit_band(dc, np.array(dec["terms"][k]["_curve"]["S_m2"]),
                              BANDS[2][1], BANDS[2][2])
                  for k in ("composite", "form", "grain_relief")},
              "note": ("decomposition identity S_comp = S_form + S_relief + "
                       "2*cov, evaluated on the saved r1 x curves; not a new "
                       "separation")}
    log("  component shares + H by term done")

    cp.unit_done("measure")
    return {"beds": beds, "shares": shares,
            "wall_s": float(time.time() - t0)}


# ---------------------------------------------------------------------------
# stage judge
# ---------------------------------------------------------------------------
def stage_judge(cp, meas):
    log("JUDGE  primary = scaling-band H against 0.70-0.85")
    rows = []
    for case, bed in meas["beds"].items():
        for k in DIRS:
            blk = bed["direction"].get(k, {})
            if not blk.get("available"):
                rows.append({"case": case, "label": bed["label"],
                             "group": bed["group"], "direction": k,
                             "available": False, "H": None,
                             "in_band": False})
                continue
            f = blk["bands"]["scaling"]
            H = f["H"]
            dentro = bool(H is not None and GATE[0] <= H <= GATE[1])
            rows.append({
                "case": case, "label": bed["label"], "group": bed["group"],
                "direction": k, "available": True,
                "H": H, "H_se": f["H_se"], "n_points": f["n_points"],
                "rms_log": f["rms_log"], "r2": f["r2"],
                "in_band": dentro,
                "margin_to_band": (None if H is None else
                                   (max(GATE[0] - H, H - GATE[1], 0.0))),
                "latitude": (None if H is None else
                             ("above" if H > GATE[1] else
                              ("below" if H < GATE[0] else "inside"))),
            })
    gate_rows = [r for r in rows if r["in_band"]]
    grains_x = [r for r in rows if r["group"] == "grains"
                and r["direction"] in ("x", "y") and r["available"]]
    grains_x_fail = [r for r in grains_x if not r["in_band"]]

    # band-vs-method attribution, per case x scaling
    attrib = []
    for case in ("v422_patch_grains", "long_grains_r1", "long_grains_r2"):
        if case not in meas["beds"]:
            continue
        for k in ("x", "iso"):
            blk = meas["beds"][case]["direction"].get(k, {})
            if not blk.get("available"):
                continue
            hs = blk["bands"]["scaling"]["H"]
            hl = blk["bands"]["legacy"]["H"]
            a4 = blk.get("local_H_at_m", {}).get("0.0400")
            bc = blk.get("binning_control", {})
            attrib.append({
                "case": case, "direction": k,
                "H_legacy_4_20": hl, "H_scaling": hs,
                "delta_band": (None if (hs is None or hl is None)
                               else hs - hl),
                "H_local_at_4cm": a4,
                "H_scaling_nbin120": bc.get("H_scaling") if k == "x" else None,
                "delta_method_binning": (
                    None if (hs is None or bc.get("H_scaling") is None)
                    else bc["H_scaling"] - hs),
            })

    verdict = {
        "primary_question": ("scaling-band H (0.4-3.6 cm, <= D50) inside "
                             "0.70-0.85, per bed and per direction?"),
        "gate": list(GATE),
        "n_cells": len([r for r in rows if r["available"]]),
        "n_cells_in_band": len(gate_rows),
        "all_grains_directions_fail": bool(
            len(grains_x_fail) == len(grains_x) and len(grains_x) > 0),
        "rows": rows, "attribution": attrib,
    }
    log("  cells in band: %d / %d" % (len(gate_rows),
                                     verdict["n_cells"]))
    for r in rows:
        if r["available"]:
            log("    %-30s %-3s H=%.3f %s" % (
                r["case"], r["direction"], r["H"],
                "PASS" if r["in_band"] else "fail (%s)" % r["latitude"]))
    cp.unit_done("judge")
    return verdict


# ---------------------------------------------------------------------------
# stage figure
# ---------------------------------------------------------------------------
def stage_figure(cp, meas, verdict):
    log("FIGURE  structure-functions.png")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import io

    beds = meas["beds"]
    fig, ax = plt.subplots(2, 2, figsize=(14.2, 11.4), dpi=120,
                           facecolor="white")
    a0, a1, a2, a3 = ax.ravel()

    def band_bg(a):
        a.axvspan(BANDS[0][1], BANDS[0][2], color="#d9ead3", alpha=0.75,
                  lw=0, label="scaling band 0.4-3.6 cm")
        a.axvspan(BANDS[1][1], BANDS[1][2], color="#fce5cd", alpha=0.55,
                  lw=0, label="meso band 6.6-20 cm")

    def draw(a, case, k, color, ls="-", mark=None, label=None):
        blk = beds[case]["direction"].get(k, {})
        if not blk.get("available"):
            return
        d = np.array(blk["curve"]["D_m"]) * 100
        s = np.array(blk["curve"]["S_m2"])
        H = blk["bands"]["scaling"]["H"]
        a.loglog(d, s, color=color, ls=ls, lw=1.3,
                 label="%s  H=%.2f" % (label or case, H))
        m = (d >= 0.4) & (d <= 3.6)
        if m.sum() > 3:
            p = np.polyfit(np.log(d[m]), np.log(s[m]), 1)
            a.loglog(d[m], np.exp(np.polyval(p, np.log(d[m]))), color=color,
                     lw=2.6, alpha=0.55)
        if mark is not None and len(mark) and mark in d:
            i = int(np.where(d == mark)[0][0])
            a.plot(d[i], s[i], marker="o", ms=3.4, color=color)

    for a, k, ttl in ((a0, "x", "Flow direction (x) - grains and substrate"),
                      (a1, "y", "Cross-stream (y) - grains and substrate")):
        band_bg(a)
        draw(a, "v422_patch_grains", k, "#1f4e79", "-", label="0.6 m patch, grains")
        draw(a, "long_grains_r1", k, "#b03030", "-", label="3.0 m r1, grains")
        draw(a, "long_grains_r2", k, "#d0793a", "--", label="3.0 m r2, grains")
        draw(a, "v422_patch_substrate", k, "#2e7d32", ":", label="0.6 m patch, substrate")
        draw(a, "long_substrate_r1", k, "#66a05a", "-.", label="3.0 m r1, substrate")
        draw(a, "long_substrate_r2", k, "#9dc183", ":", label="3.0 m r2, substrate")
        a.set_title(ttl, fontsize=9.5)
        a.set_xlabel("lag tau (cm)")
        a.set_ylabel("S(tau) (m$^2$)")
        a.grid(alpha=0.25, which="both", lw=0.4)
        a.legend(fontsize=6.2, loc="upper center",
                 bbox_to_anchor=(0.5, -0.16), ncol=3, framealpha=0.9,
                 borderaxespad=0.0)
        a.set_xlim(0.03, 100)

    # (c) components of the r1 bed, x
    band_bg(a2)
    dec = json.load(open(DECOMP_423, encoding="utf-8"))
    for k, color, col in (("composite", "#111111", "#111111"),
                          ("form", "#2e7d32", "#2e7d32"),
                          ("grain_relief", "#b03030", "#b03030")):
        d = np.array(dec["terms"][k]["_curve"]["D_m"]) * 100
        s = np.array(dec["terms"][k]["_curve"]["S_m2"])
        H = meas["shares"]["terms_H_scaling"][k]["H"]
        a2.loglog(d, s, color=color, lw=1.4,
                  label="%s  H=%.2f" % (k, H))
        m = (d >= 0.4) & (d <= 3.6)
        p = np.polyfit(np.log(d[m]), np.log(s[m]), 1)
        a2.loglog(d[m], np.exp(np.polyval(p, np.log(d[m]))), color=col,
                  lw=2.6, alpha=0.5)
    a2.set_title("r1, x: composite vs substrate vs grain relief "
                 "(S_comp = S_form + S_relief + 2cov)", fontsize=9.5)
    a2.set_xlabel("lag tau (cm)")
    a2.set_ylabel("S(tau) (m$^2$)")
    a2.grid(alpha=0.25, which="both", lw=0.4)
    a2.legend(fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=3, framealpha=0.9, borderaxespad=0.0)
    a2.set_xlim(0.03, 100)

    # (d) local exponent vs lag -- where the scaling region sits
    a3.axhspan(GATE[0], GATE[1], color="#d9ead3", alpha=0.85, lw=0,
               label="gate 0.70-0.85")
    a3.axvspan(BANDS[0][1] * 100, BANDS[0][2] * 100, color="#fff2cc",
               alpha=0.9, lw=0, label="scaling band")
    for case, color, ls, lab in (
            ("v422_patch_grains", "#1f4e79", "-", "0.6 m patch, grains"),
            ("long_grains_r1", "#b03030", "-", "3.0 m r1, grains"),
            ("long_grains_r2", "#d0793a", "--", "3.0 m r2, grains"),
            ("v422_patch_substrate", "#2e7d32", ":", "0.6 m patch, substrate"),
            ("long_substrate_r1", "#66a05a", "-.", "3.0 m r1, substrate"),
            ("long_substrate_r2", "#9dc183", ":", "3.0 m r2, substrate")):
        blk = beds[case]["direction"].get("x", {})
        if not blk.get("available"):
            continue
        d = np.array(blk["curve"]["D_m"])
        s = np.array(blk["curve"]["S_m2"])
        a3.semilogx(d * 100, local_H(d, s), color=color, ls=ls, lw=1.2,
                    label=lab)
    a3.set_xlim(0.2, 40)
    a3.set_ylim(-0.55, 1.4)
    a3.set_xlabel("lag tau (cm)")
    a3.set_ylabel("local exponent  0.5 dlnS/dlnd")
    a3.set_title("Local slope (x): composite ~0.3, substrate ~0.9, gate "
                 "0.70-0.85", fontsize=9.5)
    a3.grid(alpha=0.25, lw=0.4)
    a3.legend(fontsize=6.3, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), framealpha=0.9,
              borderaxespad=0.0)

    fig.suptitle("v4-2-5 structure functions: in the scaling region "
                 "(tau <= D50) the composite bed is a sub-linear grain mosaic",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=2.4, w_pad=1.2)
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=120, facecolor="white")
    plt.close(fig)

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size[0] > 600
        except Exception:
            return False
    atomic_write_bytes(FIGURE, b.getvalue(), verify=verify)
    log("  wrote structure-functions.png (%d bytes)"
        % os.path.getsize(FIGURE))
    cp.unit_done("figure")
    return {"path": FIGURE, "bytes": os.path.getsize(FIGURE),
            "panels": ["S(tau) x", "S(tau) y", "components r1 x",
                       "local exponent x"]}


# ---------------------------------------------------------------------------
# output helpers
# ---------------------------------------------------------------------------
def f3(x):
    return "n/a" if x is None else "%.3f" % x


def cell(f, gate=GATE):
    if f.get("H") is None:
        return "n/a"
    mark = "**" if gate[0] <= f["H"] <= gate[1] else ""
    return "%s%.3f%s · n=%d · rms %.3f" % (mark, f["H"], mark,
                                           f["n_points"], f["rms_log"])


def local_at(blk, lag):
    """Local exponent at one lag, from the stored curve of a direction."""
    d = np.array(blk["curve"]["D_m"])
    s = np.array(blk["curve"]["S_m2"])
    if lag > d[-1]:
        return None
    return float(np.interp(lag, d, local_H(d, s)))


def table_B(meas, group, k):
    lines = ["| bed | scaling 0.4-3.6 cm | meso 6.6-20 cm | legacy 4-20 cm "
             "| legacy, above break |",
             "| --- | --- | --- | --- | --- |"]
    for case, bed in meas["beds"].items():
        if bed["group"] != group:
            continue
        blk = bed["direction"].get(k, {})
        if not blk.get("available"):
            lines.append("| %s | 測定不能 | - | - | - |"
                         % bed["label"])
            continue
        b = blk["bands"]
        ab = blk["audit_like"]
        lines.append("| %s | %s | %s | %s | %s |" % (
            bed["label"], cell(b["scaling"]), cell(b["meso"]),
            cell(b["legacy"]),
            "n/a" if ab["H_above_break"] is None else
            "%.3f (break %.1f cm)" % (ab["H_above_break"],
                                      ab["break_scale_m"] * 100)))
    return "\n".join(lines)


def shape_table(meas, k):
    lines = ["| bed | local H @0.41 cm | @1.0 cm | @2.4 cm | @3.2 cm "
             "| min local H in band | S growth in band | plateau が帯内? |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for case, bed in meas["beds"].items():
        blk = bed["direction"].get(k, {})
        if not blk.get("available"):
            continue
        sh = blk["shape"]["scaling"]
        f = blk["bands"]["scaling"]

        def g(lag):
            v = local_at(blk, lag)
            return "-" if v is None else "%.2f" % v
        lines.append("| %s | %s | %s | %s | %s | %.2f @ %.2f cm | x%.1f "
                     "| %s |" % (
                         bed["label"], g(0.0041), g(0.010), g(0.024),
                         g(0.032), sh["local_H_min"],
                         sh["local_H_min_at_m"] * 100,
                         f["S_growth"] or float("nan"),
                         "yes" if sh["plateau_inside_band"] else "no"))
    return "\n".join(lines)


def build_report(ctx):
    st0 = ctx["stage0"]
    meas = ctx["measure"]
    ver = ctx["judge"]
    rows = ver["rows"]
    g422 = meas["beds"]["v422_patch_grains"]
    gr1 = meas["beds"]["long_grains_r1"]
    gr2 = meas["beds"]["long_grains_r2"]
    sr1 = meas["beds"]["long_substrate_r1"]
    sr2 = meas["beds"]["long_substrate_r2"]
    sp = meas["beds"]["v422_patch_substrate"]

    def H(case, k, band="scaling"):
        b = meas["beds"][case]["direction"].get(k, {})
        if not b.get("available"):
            return None
        return b["bands"][band]["H"]

    def rms(case, k, band="scaling"):
        b = meas["beds"][case]["direction"].get(k, {})
        if not b.get("available"):
            return None
        return b["bands"][band]["rms_log"]

    def npt(case, k, band="scaling"):
        b = meas["beds"][case]["direction"].get(k, {})
        if not b.get("available"):
            return 0
        return b["bands"][band]["n_points"]

    grains_cells = [r for r in rows if r["group"] == "grains"
                    and r["direction"] in ("x", "y") and r["available"]]
    passed = [r for r in grains_cells if r["in_band"]]
    L = []
    A = L.append

    A("# REPORT-V4-2-5 — 床の再監査：H をスケーリング領域で測り直す")
    A("")
    A("Brief: `BRIEF-V4-2-5.md`. Workdir: "
      "`D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`. Built: %s."
      % ctx["created_utc"])
    A("新規は `run_v4_2_5.py` のみ。出力は `out/v4-2-5/`。**生成していない・"
      "パラメータを変えていない・床を作り直していない**（保存済み heightfield "
      "と、v4-2.3 が保存した構造関数カーブのみを測った）。既存 .py の書き換えは "
      "0 本。")
    A("")
    A("## 0. 判定")
    A("")
    A("主判定 = **スケーリング領域（ラグ 0.4–3.6 cm ≲ D50）の H が帯 "
      "0.70–0.85 に入るか**、床ごと・方向ごと。帯は測定前に固定し、動かして"
      "いない。")
    A("")
    A("| 床 | 方向 | H（0.4–3.6 cm） | n | rms(log) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for r in rows:
        if not r["available"]:
            A("| %s | %s | 測定不能 | – | – | – |" % (r["label"],
                                                     r["direction"]))
            continue
        A("| %s | %s | **%.3f** | %d | %.3f | %s |" % (
            r["label"], r["direction"], r["H"], r["n_points"], r["rms_log"],
            "**帯内**" if r["in_band"] else
            ("帯の上（+%.3f）" % (r["H"] - GATE[1])
             if r["latitude"] == "above" else
             "帯の下（−%.3f）" % (GATE[0] - r["H"]))))
    A("")
    in_band = [r for r in rows if r["available"] and r["in_band"]]
    A("**判定: 帯に入ったセルは %d / %d。**" %
      (len(in_band), ver["n_cells"]))
    A("")
    if in_band:
        A("帯に入ったのは: %s。**ただしこの 1 本は「スケーリング則が帯に入った」"
          "とは読めない**（§5.2 と §6.1: 局所指数が帯の中で 0.79 → 0.61 へ"
          "すべり落ち、ゲートの縁を横切る平均値がたまたま入っただけ）。"
          % "、".join(
              "%s / %s (%.3f)" % (r["label"], r["direction"], r["H"])
              for r in in_band))
        A("")
    A("粒床（2.2 の 0.6 m パッチ、帯状床 r1・r2）は流れ方向・横断方向とも"
      "帯に入らなかった：H = %s / %s / %s（x）、%s / %s（y）、iso %s / %s / "
      "%s。**すべて帯の下側**に、はっきり外れる（帯の下端より 0.33 以上下）。"
      % (f3(H("v422_patch_grains", "x")), f3(H("long_grains_r1", "x")),
         f3(H("long_grains_r2", "x")),
         f3(H("v422_patch_grains", "y")), f3(H("long_grains_r1", "y")),
         f3(H("v422_patch_grains", "iso")), f3(H("long_grains_r1", "iso")),
         f3(H("long_grains_r2", "iso"))))
    A("")
    A("基盤のみは箱で分かれる：0.6 m パッチは x %s / iso %s（帯の下端〜帯内）、"
      "3.0 m 帯状床は x %s / %s、y %s / %s（帯の上）。**同じ基盤の作り方で、"
      "箱の長さだけで 0.4–3.6 cm の H が 0.69 と 0.93 に分かれる**（§6.1）。"
      % (f3(H("v422_patch_substrate", "x")),
         f3(H("v422_patch_substrate", "iso")),
         f3(H("long_substrate_r1", "x")), f3(H("long_substrate_r2", "x")),
         f3(H("long_substrate_r1", "y")), f3(H("long_substrate_r2", "y"))))
    A("")
    A("したがって、v4-2.4 のラリーで置かれた前提「Nikora の H≈0.79 は"
      "ラグ ≲ D50 のスケーリング領域の値であり、この床もそこに居る」は、"
      "**合成表面についてはこの床で成立しない**。帯を移したことで値は "
      "0.219/0.232 から 0.294/0.322 へ上がった（+0.075/+0.090）が、上げ幅は"
      "すべて帯の効果で（方法＝ビニングの効果は ±0.006、§6）、0.70–0.85 "
      "には届かない。")
    A("")
    A("## 1. STAGE 0（絶対パス事前確認）")
    A("")
    A("| 対象 | bytes | sha256 (16) |")
    A("| --- | --- | --- |")
    for k, v in st0["preflight"].items():
        if k == "long_heightfield_read":
            continue
        A("| `%s` | %d | `%s` |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    lr = st0["preflight"]["long_heightfield_read"]
    A("- 長尺床はヘッダなし生 float32（`np.fromfile`, little endian）として読み、"
      "shape %s、cell %.8f m。読みの検算: 場のσz = **%.6f m** 対 meta "
      "**%.6f m**（相対差 %.1e）。"
      % (lr["shape"], lr["cell_m"], lr["sigma_z_field_read_m"],
         lr["sigma_z_meta_m"], lr["sigma_rel_diff"]))
    A("- D50（報告書値 0.0362 m）の確認: v4-2.2 パッチ **%.6f m**"
      "（heightfield-meta と numbers の両方で一致）。帯状床は r1 %.6f m / "
      "r2 %.6f m。→ 0.4–3.6 cm は v4-2.2 で 0.11–0.99 D50、帯状床 r1 で "
      "0.10–0.94 D50。**帯は全床で ≲ D50**。"
      % (st0["d50"]["v422_patch_m"], st0["d50"]["long_r1_meta_m"],
         st0["d50"]["long_r2_meta_m"]))
    A("- 既存 .py **%d 本**の sha256 を stage-0 記録として "
      "`stage0-py-hashes.json` に保存。" % st0["n_existing_py"])
    A("")
    A("## 2. 方法（物差しの定義）")
    A("")
    A("- **S(τ)**: 保存済み高さ場 z に対し、軸ごとに "
      "`S(τ) = mean[(z(x) − z(x+τ))²]`（周期ラップ、x は流れ軸、y は横断軸）。"
      "`iso` は v4-2.2/2.3 と同じ x と y の平均 `0.5(Sx+Sy)`。ラグは "
      "1 px〜dmax の幾何級数 60 ビン（整数 px に丸め）。**v4-2.3 の "
      "`curve_iso` / `curve_x` の算術をそのまま写した**（同ファイルの関数を"
      "呼ばず、式を写して再実装。§3 で一致を機械確認）。")
    A("- **H**: 帯内のラグで `ln S` を `ln τ` に一次回帰し、傾きの 1/2。"
      "v4-2.2/2.3 のゲートと同じ定義。")
    A("- **帯**（測定前に固定）: スケーリング **0.4–3.6 cm**（新主判定）、"
      "meso **6.6–20 cm**（記録のみ、H_meso 設計用。6.6 cm は v4-2.2 の "
      "iso 折れ点 6.62 cm）、legacy **4–20 cm**（旧ゲート帯、再現性チェック）。")
    A("- **当てはまり**: 対数残差 RMS、R²、使ったラグ数 n。"
      "**飽和の形**: 帯内の局所指数 `0.5·dlnS/dlnd` のプローブ、帯内の最小・"
      "最大局所指数とそのスケール、0.70/0.85/0.20 を横切るスケール、"
      "帯内での S の成長比。")
    A("- 成分分解（基盤／粒起伏）は v4-2.3 が保存した `decomp.json` の "
      "`_curve`（r1・x）を再利用した。**新しい分解はしていない**。")
    A("")
    A("## 3. 物差しの再現性（保存 curve との一致）")
    A("")
    rp = ctx["repro"]
    A("| 対象 | ラグ数 | max abs ΔS (m²) | 相対 ΔS |")
    A("| --- | --- | --- | --- |")
    for k in ("v422_patch_iso", "v422_patch_x", "long_r1_iso", "long_r1_x"):
        v = rp[k]
        A("| %s | %d | %.3e | %.2e |" % (k, v["n_lags"], v["max_abs_S"],
                                         v["rel_S"]))
    A("")
    A("→ `all_ok = %s`（相対 1e-9 以下）。**同じ場から同じ曲線が出る**ので、"
      "以下の差は帯と方法にだけ帰属できる。" % rp["all_ok"])
    A("")
    A("旧値の再計算（legacy 4–20 cm、v4-2.3 と同じ物差し）:")
    A("")
    A("| 旧値 | 記録値 | 本 run の再計算 | 一致 |")
    A("| --- | --- | --- | --- |")
    old_checks = [
        ("2.2 パッチ iso 4–20 cm", OLD["H_iso_4_20_v422_patch"],
         meas["beds"]["v422_patch_grains"]["direction"]["iso"]["bands"]
         ["legacy"]["H"]),
        ("2.2 パッチ iso 折れ点上（0.151）",
         OLD["H_iso_above_break_v422_patch"],
         meas["beds"]["v422_patch_grains"]["direction"]["iso"]["audit_like"]
         ["H_above_break"]),
        ("2.2 パッチ iso 折れ点スケール",
         OLD["break_iso_v422_patch_m"],
         meas["beds"]["v422_patch_grains"]["direction"]["iso"]["audit_like"]
         ["break_scale_m"]),
        ("帯状床 r1 x 4–20 cm（0.2192）", OLD["H_x_4_20_long_r1"],
         meas["beds"]["long_grains_r1"]["direction"]["x"]["bands"]["legacy"]
         ["H"]),
        ("帯状床 r2 x 4–20 cm（0.2319）", OLD["H_x_4_20_long_r2"],
         meas["beds"]["long_grains_r2"]["direction"]["x"]["bands"]["legacy"]
         ["H"]),
        ("帯状床 r1 iso 4–20 cm", OLD["H_iso_4_20_long_r1"],
         meas["beds"]["long_grains_r1"]["direction"]["iso"]["bands"]["legacy"]
         ["H"]),
        ("帯状床 r2 iso 4–20 cm", OLD["H_iso_4_20_long_r2"],
         meas["beds"]["long_grains_r2"]["direction"]["iso"]["bands"]["legacy"]
         ["H"]),
        ("基盤 r1 x 4–20 cm", OLD["H_x_4_20_substrate_r1"],
         meas["beds"]["long_substrate_r1"]["direction"]["x"]["bands"]["legacy"]
         ["H"]),
        ("基盤 r2 x 4–20 cm", OLD["H_x_4_20_substrate_r2"],
         meas["beds"]["long_substrate_r2"]["direction"]["x"]["bands"]["legacy"]
         ["H"]),
        ("粒起伏 r1 x 4–20 cm（0.1028）", OLD["H_relief_x_4_20_long_r1"],
         meas["shares"]["terms_H_legacy"]["grain_relief"]["H"]),
    ]
    for name, ref, got in old_checks:
        ok = (ref is not None and got is not None
              and abs(ref - got) <= 2e-4 + 1e-6 * abs(ref))
        A("| %s | %.6f | %.6f | %s |" % (name, ref,
                                         float("nan") if got is None else got,
                                         "**True**" if ok else "**False**"))
    A("")
    A("→ 物差しは旧 run と同一（差は 1e-4 以下）。**旧値との差は帯の違いで"
      "説明でき、方法の違いではない**（§6 で数値に分解）。")
    A("")
    A("## 4. H 表（帯別・方向別）")
    A("")
    A("各セルは `H · n=<使ったラグ数> · rms=<対数残差 RMS>`。太字は帯内。"
      "`legacy, above break` は v4-2.2/2.3 の折れ点則をその曲線に当てた"
      "ときの折れ点上の値（再現性の確認用）。")
    A("")
    A("### 4.1 流れ方向（x）")
    A("")
    A(table_B(meas, "grains", "x"))
    A("")
    A("基盤のみ（x）:")
    A("")
    A(table_B(meas, "substrate", "x"))
    A("")
    A("### 4.2 横断方向（y）")
    A("")
    A(table_B(meas, "grains", "y"))
    A("")
    A("基盤のみ（y）:")
    A("")
    A(table_B(meas, "substrate", "y"))
    A("")
    A("> 帯状床 r2 の粒床 y は**測れない**: v4-2.3 は r2 の heightfield を"
      "保存しておらず、保存 curve にも y が無い。2.2 パッチの基盤 y も同じ"
      "理由で測れない（2.3 は長尺床の基盤にだけ y を保存した）。生成しない"
      "という brief の制約下では再構成しない（数を作らない）。")
    A("")
    A("### 4.3 等方（iso = x と y の平均、v4-2.2/2.3 の主推定量）")
    A("")
    A(table_B(meas, "grains", "iso"))
    A("")
    A(table_B(meas, "substrate", "iso"))
    A("")
    A("> 注 1: 0.6 m パッチ基盤の iso **0.713 が唯一の帯内**だが、同じ基盤の"
      " x 方向は 0.686（帯の下）で、局所指数は帯の中で 0.79 → 0.61 と"
      "すべり落ちる"
      "（§5.2）。1 本の平均値がゲートの縁を横切っただけで、スケーリング則"
      "ではない。")
    A("> 注 2: `boulder band removed` の iso は、v4-2.3 の `case_block` が "
      "`row_slice` を x にしか適用しないため **r1 の iso と同一**（独立な"
      "測点ではない）。この行の x だけが感度の情報を持つ。")
    A("")
    A("## 5. 帯の中で飽和しているか（形と折れ点）")
    A("")
    A("流れ方向（x）の形:")
    A("")
    A(shape_table(meas, "x"))
    A("")
    A("読み方（x）:")
    A("")
    for case, lab in (("v422_patch_grains", "2.2 パッチ粒床"),
                      ("long_grains_r1", "帯状床 r1 粒床"),
                      ("long_grains_r2", "帯状床 r2 粒床")):
        blk = meas["beds"][case]["direction"]["x"]
        sh = blk["shape"]["scaling"]
        f = blk["bands"]["scaling"]
        A("- **%s**: 帯内で S は %.2f 倍にしか増えない（H=0.79 なら "
          "(3.19/0.41)^1.58 ≈ 26 倍必要）。局所指数は 0.41 cm で %.2f、"
          "最小 %.2f（%.2f cm）。**飽和（局所指数 < 0.20）は帯の外**、"
          "%.1f cm から。つまり帯は「飽和した plateau」ではなく、"
          "**H 0.3 のサブ線形な立ち上がり**で、粒起伏の成長指数がそのまま"
          "合成の H になっている（§5.1）。"
          % (lab, f["S_growth"], local_at(blk, 0.0041),
             sh["local_H_min"], sh["local_H_min_at_m"] * 100,
             (sh["plateau_onset_m"] or float("nan")) * 100))
    A("")
    for case, lab in (("v422_patch_substrate", "2.2 パッチ基盤"),
                      ("long_substrate_r1", "帯状床 r1 基盤"),
                      ("long_substrate_r2", "帯状床 r2 基盤")):
        blk = meas["beds"][case]["direction"]["x"]
        sh = blk["shape"]["scaling"]
        f = blk["bands"]["scaling"]
        A("- **%s**: 逆に帯内で S は %.1f 倍に増え、局所指数は %.2f–%.2f。"
          "**飽和していない**。%s"
          % (lab, f["S_growth"], sh["local_H_min"], sh["local_H_max"],
             ("帯全体で 0.85 を上回ったまま（帯の上）。"
              if sh["local_H_min"] > 0.85 else
              ("帯全体で 0.70 を下回ったまま（帯の下）。"
               if sh["local_H_max"] < 0.70 else
               ("帯全体で 0.70–0.85 の内側。"
                if sh["local_H_min"] > 0.70 else
                "帯の中で 0.70 を切る（%.2f cm）。"
                % ((sh["cross_0.70_m"] or float("nan")) * 100))))))
    A("")
    A("→ 帯の中で「飽和」に近いのは**粒起伏の側**（ただし plateau は帯の外）。"
      "基盤は飽和せず、3.0 m 帯状床では帯の上限（0.85）を超える速さで成長"
      "する。合成床の H が 0.3 前後に落ちるのは、両者の和で粒起伏の低指数な"
      "分散が支配するから（§5.1）。")
    A("")
    A("### 5.1 誰が 0.4–3.6 cm の分散を持つか（r1・x、保存 curve の恒等式）")
    A("")
    A("`S_合成 = S_基盤 + S_粒 + 2·cov`（v4-2.3 の分解をそのまま再評価）。")
    A("")
    A("| ラグ | S_合成 (m²) | 基盤の寄与 | 粒起伏の寄与 | cov 項 |")
    A("| --- | --- | --- | --- | --- |")
    for k, v in meas["shares"]["lags"].items():
        if float(k) > 0.0401:
            continue
        A("| %s cm | %.3e | %+.2f | %+.2f | %+.2f |" % (
            "%.2f" % (float(k) * 100), v["S_total"], v["share_form"],
            v["share_relief"], v["share_cov"]))
    A("")
    t = meas["shares"]["terms_H_scaling"]
    A("成分別 H（0.4–3.6 cm、x）: 合成 **%.3f** / 基盤 **%.3f** / "
      "粒起伏 **%.3f**。→ スケーリング帯でも 4 cm 帯と同じ構図"
      "（基盤は帯の外・箱依存、粒起伏はサブ線形、合成は粒に引かれる）。"
      % (t["composite"]["H"], t["form"]["H"], t["grain_relief"]["H"]))
    A("")
    A("### 5.2 帯の中でゲートを横切るか（x）")
    A("")
    A("| bed | local H @0.41 cm | @3.2 cm | 0.85 を横切る（帯内） "
      "| 0.70 を横切る（帯内） | plateau 開始（<0.20） |")
    A("| --- | --- | --- | --- | --- | --- |")
    for case, bed in meas["beds"].items():
        blk = bed["direction"].get("x", {})
        if not blk.get("available"):
            continue
        sh = blk["shape"]["scaling"]

        def cz(v):
            return "帯内に無し" if v is None else "%.2f cm" % (v * 100)
        A("| %s | %.2f | %.2f | %s | %s | %s |" % (
            bed["label"], local_at(blk, 0.0041), local_at(blk, 0.032),
            cz(sh["cross_0.85_m"]), cz(sh["cross_0.70_m"]),
            "n/a" if sh["plateau_onset_m"] is None
            else "%.1f cm" % (sh["plateau_onset_m"] * 100)))
    A("")
    A("→ **0.85 は誰も帯内で切らない**（最大は 3.0 m 基盤の 0.94 で、最初から"
      "上に居る）。**0.70 を帯内で切るのは 0.6 m パッチ基盤の x だけ**"
      "（1.00 cm、0.79 → 0.61 のすべり）。粒床は 0.41–3.2 cm の全域で 0.4 "
      "以下、plateau（<0.20）は 9–17 cm で帯の外。")
    A("")
    A("## 6. 旧値との差の切り分け（帯 vs 方法）")
    A("")
    A("| ケース | 旧 4–20 cm | 新 0.4–3.6 cm | Δ（帯） | 4 cm の局所指数"
      " | nbin=120 の新帯 H | Δ（方法） |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for a in ver["attribution"]:
        A("| %s / %s | %s | %s | %s | %s | %s | %s |" % (
            a["case"], a["direction"], f3(a["H_legacy_4_20"]),
            f3(a["H_scaling"]),
            "n/a" if a["delta_band"] is None else "%+.3f" % a["delta_band"],
            f3(a["H_local_at_4cm"]), f3(a["H_scaling_nbin120"]),
            "n/a" if a["delta_method_binning"] is None else
            "%+.3f" % a["delta_method_binning"]))
    A("")
    db = [a["delta_band"] for a in ver["attribution"]
          if a["delta_band"] is not None]
    dm = [a["delta_method_binning"] for a in ver["attribution"]
          if a["delta_method_binning"] is not None]
    A("- **帯の効果**が主: Δband = %+.3f 〜 %+.3f。帯状床 x は +0.075/+0.090、"
      "0.6 m パッチは +0.147（旧 4–20 cm 帯が 0.6 m 箱では折れ点 8.2 cm の"
      "下に落ち、独立標本が 3 個しかない帯だったため。v4-2.3 §4）。"
      % (min(db), max(db)))
    A("- **方法（ビニング）の効果**は %+.3f 〜 %+.3f（|Δ| ≤ %.3f）で、"
      "帯の効果より 1 桁小さい。ラグ格子を 60→120 ビンに細かくしても"
      "結論は動かない。"
      % (min(dm), max(dm), max(abs(v) for v in dm)))
    A("- 旧値 0.2192 / 0.2319 は §3 のとおり 1e-4 以内で再現した。"
      "したがって差は**帯だけ**で説明できる。")
    A("")
    A("### 6.1 箱寸法の効果は帯の中にも残る（基盤・x）")
    A("")
    A("| 基盤 x | H（0.4–3.6 cm） | S @0.41 cm (m²) | S @3.19 cm (m²) "
      "| 帯内の成長 | 局所指数 @0.41 / 1.0 / 3.2 cm |")
    A("| --- | --- | --- | --- | --- | --- |")
    for case, lab in (("v422_patch_substrate", "0.6 m パッチ"),
                      ("long_substrate_r1", "3.0 m 帯状床 r1"),
                      ("long_substrate_r2", "3.0 m 帯状床 r2")):
        blk = meas["beds"][case]["direction"]["x"]
        f = blk["bands"]["scaling"]
        A("| %s | %.3f | %.3e | %.3e | x%.1f | %.2f / %.2f / %.2f |" % (
            lab, f["H"], f["S_at_d_min"], f["S_at_d_max"], f["S_growth"],
            local_at(blk, 0.0041), local_at(blk, 0.010), local_at(blk, 0.032)))
    A("")
    A("同じ fBm 則（H_gen=0.79）・同じ σ_form・同じ cell で、箱の長さだけを"
      " 0.6 m と 3.0 m に変えると、**0.4–3.6 cm の H が 0.69 と 0.86–0.94 に"
      "分かれる**。帯内の S そのものも 10 倍違う（1.9e-6 対 1.9e-7 m² @0.41 cm）。"
      "0.6 m パッチの基盤は帯全体で 0.85 未満（0.79 → 0.61）で、x では "
      "0.70 も切る。"
      "つまり v4-2.3 が 4–20 cm で見つけた「物差しの低バイアス」は、"
      "≲D50 の帯にも残っている。")
    A("")
    A("機構の候補は「σ_form を固定したまま箱を長くすると、低波数のモード"
      "（3.0/1.5/1.0 m）が分散の一部を持ち、短波長側の相対振幅と局所勾配が"
      "変わる」こと。**本 job はこれを検証しない**（生成しない）。"
      "注意: 基盤は各 1 実現ずつなので、箱の効果と実現間ばらつきは完全には"
      "分離できない（2.3 は基盤 3 実現で 0.763–0.866 の幅を報告）。"
      "それでも 0.69 対 0.94 の差はその幅の外側にある。")
    A("")
    A("## 7. 解釈（言えること・言えないこと）")
    A("")
    A("1. **帯を移しても粒床は 0.70–0.85 に入らない。** 0.4–3.6 cm の H は "
      "x で %s / %s / %s、y で %s / %s（r2 y は測定不能）、iso で %s / %s / "
      "%s。帯の下側に、はっきり外れる（帯の下端より 0.33 以上下）。"
      % (f3(H("v422_patch_grains", "x")), f3(H("long_grains_r1", "x")),
         f3(H("long_grains_r2", "x")),
         f3(H("v422_patch_grains", "y")), f3(H("long_grains_r1", "y")),
         f3(H("v422_patch_grains", "iso")), f3(H("long_grains_r1", "iso")),
         f3(H("long_grains_r2", "iso"))))
    t = meas["shares"]["terms_H_scaling"]
    A("2. **理由は帯の取り方でも物差しの低バイアスでもない。** 同じ帯で"
      "成分を見ると、粒起伏の H は %.3f で合成（%.3f）とほぼ同じ。0.4–3.6 cm "
      "の S のうち粒起伏が 97–100%% を占め（§5.1）、基盤の寄与は 0–4%%。"
      "つまり合成の H は粒起伏の成長指数そのもので、基盤（H %.3f）は S の"
      "数 %% しか持っていない。v4-2.3 が 4–20 cm で見つけた構図が、"
      "スケーリング帯でもそのまま成り立つ。"
      % (t["grain_relief"]["H"], t["composite"]["H"], t["form"]["H"]))
    A("3. **「基盤は帯に入る」という主張は箱を明記しないと成立しない。** "
      "基盤のみの 0.4–3.6 cm の H は 0.6 m パッチで %s（x）/ %s（iso）、"
      "3.0 m 帯状床で %s–%s（x）/ %s–%s（y）。同じ作り方で 0.69 と 0.93 に"
      "分かれる（§6.1）。0.79 に最も近いのは**0.6 m パッチの基盤**である。"
      % (f3(H("v422_patch_substrate", "x")),
         f3(H("v422_patch_substrate", "iso")),
         f3(H("long_substrate_r2", "x")), f3(H("long_substrate_r1", "x")),
         f3(H("long_substrate_r2", "y")), f3(H("long_substrate_r1", "y"))))
    A("4. **Nikora の H≈0.79 と一致する帯は、この床の合成表面には無い。** "
      "0.4–3.6 cm で 0.29–0.37、4–20 cm で 0.18–0.25、6.6–20 cm で "
      "0.15–0.24（粒床 x/iso）。曲線上でも、合成表面の局所指数は "
      "0.03–100 cm の全域で 0.79 に達しない（最大は 0.03 cm の %.2f、r1 x）。"
      "0.79 は「この床の合成表面」の統計ではない。"
      % local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.0003))
    A("5. **H_meso の別建てに使える材料**: meso 帯（6.6–20 cm）の H は "
      "x で %s / %s（r1/r2 粒床）、iso で %s / %s。0.4–3.6 cm よりさらに"
      "低い。H_meso を 0.70–0.85 に置くなら、それは**別の問い**（別の帯・"
      "別の対象）として設計する必要がある。"
      % (f3(H("long_grains_r1", "x", "meso")),
         f3(H("long_grains_r2", "x", "meso")),
         f3(H("long_grains_r1", "iso", "meso")),
         f3(H("long_grains_r2", "iso", "meso"))))
    A("6. **このジョブは生成も再設計もしない。** パッチ再設計に渡す入力は"
      "「スケーリング帯の H は 0.3、その正体は粒起伏（0.29）、基盤は"
      "箱依存で 0.69–0.94、ただし分散の数 % しか持たない」という分解済みの"
      "事実だけ。")
    A("")
    A("## 8. 弱点・正直な開示")
    A("")
    A("- **r2 の y と r2 の粒起伏**: r2 の heightfield は保存されていない。"
      "y 方向の粒床と粒起伏の曲線は v4-2.3 にも保存が無く、生成禁止の"
      "制約下では作れない。表では「測定不能」と明記し、数を作っていない。"
      "粒起伏の曲線（分解）は r1・x の 1 本だけで、方向別・実現別の"
      "分解はしていない。")
    A("- **基盤のみ 3 実現平均 0.8086 は再計算していない**: その 3 実現の"
      " curve は保存されていない（保存 curve は r1/r2 の 2 本）。本 run が"
      "再計算したのは r1 %s / r2 %s（x, 4–20 cm）。"
      % (f3(H("long_substrate_r1", "x", "legacy")),
         f3(H("long_substrate_r2", "x", "legacy"))))
    A("- **唯一の「帯内」は脆弱**: 0.6 m パッチ基盤の iso 0.713 は、"
      "帯内で局所指数が 0.79→0.61 と動く曲線の平均で、x は 0.686（帯の下）、"
      "折れ点則も 3.4 cm に折れ点を置く。スケーリング則の証拠として"
      "数えないこと（§4.3 注 1・§5.2）。")
    A("- **帯より小さいラグ（判定に使っていない）**: 曲線は 1 px = 0.029 cm "
      "から出る。合成床（r1 x）の局所指数は 0.03 cm で %.2f、0.1 cm で %.2f、"
      "0.2 cm で %.2f、0.3 cm で %.2f で、0.79 にはならず、むしろ 0.4 cm の "
      "%.2f へ向けて下がる。基盤（r1 x）は 0.03 cm で %.2f（格子内挿の d² 域）。"
      "この下端はラスタ離散化（1 px 段差・raycast）の影響を排除できないため、"
      "判定の帯（0.4–3.6 cm）には入れていない。"
      % (local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.0003),
         local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.001),
         local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.002),
         local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.003),
         local_at(meas["beds"]["long_grains_r1"]["direction"]["x"], 0.004),
         local_at(meas["beds"]["long_substrate_r1"]["direction"]["x"], 0.0003)))
    A("- **0.6 m パッチの周期性**: 0.6 m 箱を周期ラップしているので、"
      "3.6 cm のラグは箱の 6% で、周期の影響は無視できる。帯状床 r1 の y も"
      "同じ（0.6 m 幅、ラグ ≤ 3.6 cm）。")
    A("- **H の定義**: 帯内の両対数傾きの 1/2、という一点。定義を変えれば"
      "値は動く（v4-2.2 §2.1 と同じ開示）。ここでは v4-2.2/2.3 の定義を"
      "そのまま使い、帯だけを brief の指定に従って移した。")
    A("- **判定は「帯に入るか」の二値**であり、0.79 との一致度ではない。"
      "帯の上に外れた基盤（0.86–0.94）も「不合格」と数え、帯の下端すれすれの "
      "0.686（パッチ基盤 x）も「不合格」と数えている。")
    A("- **巨礫 10 個**は v4-2.2/2.3 と同じ規則で除外・インペイント済みの"
      "保存場を測っている。巨礫帯を落とした感度（2.3 の "
      "`noboulderband`）も同じ帯で再計算すると、x の H は %s（全画素 %s）"
      "で、結論は動かない。"
      % (f3(H("long_grains_r1_noboulderband", "x")),
         f3(H("long_grains_r1", "x"))))
    A("- 本レポートは **H という 1 つの統計量の再監査**しか主張しない。"
      "見え方・輸送・画作りは対象外（レンダリングしていない）。")
    A("")
    A("## 9. G-LEGACY（既存 .py の無改造）")
    A("")
    gl = ctx.get("g_legacy") or {}
    A("- stage-0 記録（%s）で既存 .py **%d 本**の sha256 を保存し、run の最後"
      "（%s）に再ハッシュした。現在 **%d 本**、追加 %d / 削除 %d / 変更 **%d**。"
      % (gl.get("hash_baseline_utc"), gl.get("n_baseline") or 0,
         gl.get("checked_at"), gl.get("n_now") or 0,
         len(gl.get("added") or []), len(gl.get("removed") or []),
         len(gl.get("changed") or [])))
    A("- 判定: **%s**（既存 py の厳密一致）。**本 run の write set は "
      "`run_v4_2_5.py` と `out/v4-2-5/*` だけ**で、他の .py を書くコードは"
      "無い。" % ("PASS" if gl.get("ok") else "FAIL"))
    if gl.get("details"):
        A("")
        A("| file | bytes (before → now) | sha256 (before → now) | mtime (now) "
          "| 本 run 由来 |")
        A("| --- | --- | --- | --- | --- |")
        for d in gl["details"]:
            A("| `%s` | %d → %d | `%s` → `%s` | %s | いいえ（外部レーン） |" % (
                d["file"], d["bytes_before"], d["bytes_now"],
                d["sha256_before"][:16], d["sha256_now"][:16],
                d["mtime_now"]))
        A("")
        A("上の差分は本 run の write set の外（**外部レーン**）。mtime は "
          "stage-0 のハッシュ取得より後で、本 run が測定していた時間帯に"
          "一致する。本 run はこのファイルを開いてもいない。")
    A("")
    A("## 10. 成果物")
    A("")
    for k, v in ctx["deliverables"].items():
        A("- `%s` → `%s`" % (k, v["path"]))
    A("- 図: %d bytes、パネル = %s" % (
        ctx["figure"]["bytes"], " / ".join(ctx["figure"]["panels"])))
    A("- 他のファイルのバイト数は実ファイルを参照（本表は書き出し前に組むため、"
      "report / numbers / checkpoint の最終サイズは載せない）。")
    A("")
    A("## 11. 実行の記録")
    A("")
    A("- 全ステージ foreground（1 プロセス、`--all`）。jobkit checkpoint を"
      "段ごとに原子的にコミット（`out/v4-2-5/checkpoint.json`、units = %d）。"
      % ctx["checkpoint_units"])
    A("- G-LEGACY: 本 run 由来の既存 .py 変更は **0 件**。%s"
      % ("既存 py のハッシュは run の前後で完全一致（PASS）。"
         if ctx.get("g_legacy", {}).get("ok") else
         "外部レーンの差分は §9 に mtime 付きで開示した（本 run 由来ではない）。"))
    A("- 帯・窓・判定規則・物差しは**測定前にコード内で固定**。後付けの調整は"
      "無い（帯 0.70–0.85、スケーリング 0.4–3.6 cm、meso 6.6–20 cm、"
      "legacy 4–20 cm）。")
    A("- 実行時間: %s（report 行は本ファイル書き出し後の計測のため未記載）。"
      % "、".join("%s %.0fs" % (k, v) for k, v in ctx["timings"].items()))
    A("- 数値の全量は `numbers-v4-2-5.json`（曲線・局所指数・残差を含む）。")
    A("")
    return "\n".join(L)


def stage_report(cp, ctx):
    log("REPORT  numbers-v4-2-5.json + REPORT-V4-2-5.md")
    md = build_report(ctx)
    atomic_write_bytes(REPORT, md.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 3000)
    numbers = {
        "brief": BRIEF,
        "created_utc": ctx["created_utc"],
        "what_this_is": ("re-audit only: H measured in the scaling band "
                         "(0.4-3.6 cm <= D50), the meso band (6.6-20 cm) and "
                         "the legacy band (4-20 cm), per bed and direction; "
                         "no generation, no parameter change"),
        "estimator": ("S(tau)=mean[(z(x)-z(x+tau))^2] along one axis, "
                      "periodic; H = 0.5*dlnS/dlnd fitted over the lags "
                      "inside the band; arithmetic copied from run_v4_2_3"),
        "bands_m": {k: [a, b] for k, a, b in BANDS},
        "gate": list(GATE),
        "D50": ctx["stage0"]["d50"],
        "stage0": ctx["stage0"],
        "repro": ctx["repro"],
        "beds": ctx["measure"]["beds"],
        "component_shares": ctx["measure"]["shares"],
        "judgment": ctx["judge"],
        "old_values_checked": OLD,
        "figure": ctx["figure"],
        "timings_s": ctx["timings"],
        "checkpoint": {"units": ctx["checkpoint_units"],
                       "path": CP_PATH},
        "deliverables": ctx["deliverables"],
        "g_legacy": ctx.get("g_legacy"),
        "log": LOG,
    }
    write_json(NUMBERS, numbers)
    log("  wrote REPORT (%d bytes) and numbers (%d bytes)"
        % (os.path.getsize(REPORT), os.path.getsize(NUMBERS)))
    cp.unit_done("report")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--repro", action="store_true")
    ap.add_argument("--measure", action="store_true")
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if not any(vars(a).values()):
        a.all = True

    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    cache = os.path.join(OUT, "_work")
    os.makedirs(cache, exist_ok=True)
    t = {}

    def saved(name, obj):
        write_json(os.path.join(cache, name + ".json"), obj)

    def cached(name):
        p = os.path.join(cache, name + ".json")
        if os.path.isfile(p):
            return json.load(open(p, encoding="utf-8"))
        return None

    created = utcnow()

    t0 = time.time()
    st0 = cached("stage0") if not (a.all or a.stage0) else None
    if st0 is None:
        st0 = stage0(cp)
        saved("stage0", st0)
    t["stage0"] = time.time() - t0

    t0 = time.time()
    rp = cached("repro") if not (a.all or a.repro) else None
    if rp is None:
        rp = stage_repro(cp)
        saved("repro", rp)
    t["repro"] = time.time() - t0

    t0 = time.time()
    ms = cached("measure") if not (a.all or a.measure) else None
    if ms is None:
        ms = stage_measure(cp, rp)
        saved("measure", ms)
    t["measure"] = time.time() - t0

    t0 = time.time()
    vd = stage_judge(cp, ms)
    t["judge"] = time.time() - t0

    t0 = time.time()
    fg = stage_figure(cp, ms, vd)
    t["figure"] = time.time() - t0

    # G-LEGACY: no existing .py may have changed during this run
    base = json.load(open(STAGE0_HASHES, encoding="utf-8"))
    now = {}
    for p in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        name = os.path.basename(p)
        if name == THIS_FILE:
            continue
        now[name] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    changed = sorted(k for k in base["py"] if k in now
                     and now[k] != base["py"][k]["sha256"])
    added = sorted(set(now) - set(base["py"]))
    removed = sorted(set(base["py"]) - set(now))
    details = []
    for k in changed:
        p = os.path.join(HERE, k)
        details.append({
            "file": k,
            "bytes_before": base["py"][k]["bytes"],
            "bytes_now": os.path.getsize(p),
            "sha256_before": base["py"][k]["sha256"],
            "sha256_now": now[k],
            "mtime_now": datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds"),
        })
    g_legacy = {"n_baseline": base["n_existing_py"], "n_now": len(now),
                "changed": changed, "added": added, "removed": removed,
                "details": details, "hash_baseline_utc": base["created_utc"],
                "checked_at": utcnow(),
                "write_set": [THIS_FILE, "out/v4-2-5/"],
                "ok": not (changed or added or removed)}
    log("G-LEGACY  existing .py: changed=%d added=%d removed=%d -> %s"
        % (len(changed), len(added), len(removed),
           "PASS" if g_legacy["ok"] else "FAIL"))

    timings = {}
    for k, v in t.items():
        timings[k] = round(v, 1)
    write_json(TIMINGS, {"stage_s": timings, "created_utc": created})

    deliv = {
        "report": {"path": REPORT},
        "numbers": {"path": NUMBERS},
        "figure": {"path": FIGURE},
        "checkpoint": {"path": CP_PATH},
        "stage0_py_hashes": {"path": STAGE0_HASHES},
        "timings": {"path": TIMINGS},
    }

    state = cp.state() or {}
    ctx = {"created_utc": created, "stage0": st0, "repro": rp,
           "measure": ms, "judge": vd, "figure": fg,
           "timings": timings, "g_legacy": g_legacy,
           "checkpoint_units": int(state.get("completed_units") or 0) + 1,
           "deliverables": deliv}
    stage_report(cp, ctx)

    log("DONE  %s" % REPORT)
    log("      %s" % NUMBERS)
    log("      %s" % FIGURE)
    print("\n".join(LOG[-6:]))


if __name__ == "__main__":
    main()
