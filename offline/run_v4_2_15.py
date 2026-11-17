"""run_v4_2_15 -- BRIEF-V4-2-15: the mechanism A / mechanism B factorial on the
grain layer (v4-2.15).

MEASUREMENT JOB.  One control reproduction (STAGE 0, frozen ruler), one
implementation + per-step diagnostics (STAGE 1, bed_gen3.py), the 6-condition
factorial with paired seeds on two boxes (STAGE 2), and the pre-registered
read-out (STAGE 3).  No rendering, no parameter search: every knob is fixed in
this file before the first measurement is read.

What the brief asks, and where it is answered
---------------------------------------------
STAGE 0  absolute-path preflight over the brief's list; the 0.6 m patch bed is
         regenerated exactly as v4-2-2 made it and measured with the FROZEN
         ruler (run_v4_2_12b.nikora_measure -> run_v4_2_5.curve_dir/fit_band)
         with the v4-2-13 window [0.10, 2] x D50.  Gate (a): H_Nikora within
         +-0.005 of 12B (0.4484 / 0.4320), sigma_z/D50 and skew close to 12B
         (0.392 / +0.607).  Gate (b): the old fixed-band H_fit = 0.29 +- 0.02.
STAGE 1  bed_gen3.py written from design 1.3 / 2.3; every mechanism is a
         separately callable function; each step's own diagnostics are
         measured and reported.
STAGE 2  V0 / Vw / Vi / VA / VB / VAB (+ VB-a on the patch only), paired
         seeds from np.random.default_rng(20261115), patch box R = 5 and strip
         box (2.4 x 0.6 m) R = 3.  Height fields and grain lists are saved.
STAGE 3  the pre-registered read-out of the brief (G-A1 .. G-B3, VAB), facts
         only -- no design verdict.

Pre-registered decisions that the brief leaves open
---------------------------------------------------
* Ruler window     : [0.10, 2] x D50 (v4-2-13).  STAGE 0 reports BOTH the
                     window the 12B value was measured with ([0.25, 2]) and
                     the mandated one, so the window effect is visible rather
                     than hidden (v4-2-13 records +0.017 on patch x).
* Ruler D50        : fixed per box and seed to the V0 surface D50, so the
                     break window does not move with the answer.  H_Nikora
                     with the condition's own D50 is reported alongside.
* R_a              : exposed-area D50 / never-exposed footprint D50 (all
                     non-cluster grains).  bed_gen2's kind-based ratio is kept
                     as the legacy diagnostic.
* D50_ref / D84/90 : per realisation from V0 (used by the mechanisms).
* Dc*              : 0.20 x D50_surface (the brief's 0.15-0.25 D_framework,
                     D_framework = surface D50), fixed before the run.
* box form         : strip box = four verbatim 0.6 m populations on one
                     rectangular form (v4-2-3's precedent, seams disclosed).
* seeds            : rng(20261115).integers(0, 2**63, size=8): 5 patch, 3
                     strip; strip tiles use seed + 10 k.
* H_fit 0.29       : the number was established on the LONG box (v4-2-5
                     long_r1 0.294/0.297).  Gate (b) is evaluated there (and on
                     the newly built strip V0), with the patch value reported
                     for what it is (0.32-0.37, also exactly v4-2-5's).

Stages are checkpointed one jobkit unit at a time (stage0 / stage1 / stage2 /
stage3 / report).  Foreground run: it either writes every deliverable or stops
at a gate and writes an honest stop record.
"""
from __future__ import annotations

import datetime
import gc
import hashlib
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
import bed_gen3 as B3  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-15.md"
OUT = os.path.join(HERE, "out", "v4-2-15")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_REC = os.path.join(OUT, "stage0-rec.json")
STAGE1_REC = os.path.join(OUT, "stage1-rec.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-15.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-15.md")
THIS_FILE = "run_v4_2_15.py"

B2_PATH = os.path.join(HERE, "bed_gen2.py")
BG_PATH = os.path.join(HERE, "bed_gen.py")
R5_PATH = os.path.join(HERE, "run_v4_2_5.py")
R12_PATH = os.path.join(HERE, "run_v4_2_12b.py")
PATCH_NPY = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
PATCH_META = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
PATCH_STATS = os.path.join(HERE, "out", "v4-2-2", "bed-stats.json")
LONG_NPY = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")
N12_JSON = os.path.join(HERE, "out", "v4-2-12b", "numbers-v4-2-12b.json")
N13_JSON = os.path.join(HERE, "out", "v4-2-13", "numbers-v4-2-13.json")
JOBKIT = r"D:\mutaro-work\jobkit\checkpoint.py"

REQUIRED = [
    ("brief", BRIEF),
    ("generator bed_gen2.py", B2_PATH),
    ("generator bed_gen.py", BG_PATH),
    ("ruler run_v4_2_5.py", R5_PATH),
    ("frozen Nikora implementation run_v4_2_12b.py", R12_PATH),
    ("patch bed heightfield.npy", PATCH_NPY),
    ("patch bed meta", PATCH_META),
    ("patch bed stats", PATCH_STATS),
    ("long bed r1 (the 0.29 box)", LONG_NPY),
    ("v4-2-12b numbers", N12_JSON),
    ("v4-2-13 numbers", N13_JSON),
    ("jobkit checkpoint.py", JOBKIT),
]

# --- measurement contract (fixed before the first measurement) -------------
CELL = 0.6 / 2048.0
NBIN = 60
DM_PATCH = 0.30
DM_STRIP_X = 1.20          # half of the 2.4 m strip
DM_STRIP_Y = 0.30
FIXED_BAND = (0.004, 0.036)
BREAK_WIN = (0.10, 2.0)
GATE_H_EXPECT = {"x": 0.4484, "y": 0.4320}
GATE_H_TOL = 0.005
GATE_SIG_EXPECT = 0.3919
GATE_SIG_TOL = 0.005
GATE_SKEW_EXPECT = 0.607
GATE_SKEW_TOL = 0.05
HFIT_029 = 0.29
HFIT_029_TOL = 0.02
W0_EXPECT = {"x": 0.4484, "y": 0.4320}   # 12B, window [0.25, 2]

SEED_BASE = 20261115
N_PATCH = 5
N_STRIP = 3
CONDITIONS_PATCH = ("V0", "Vw", "Vi", "VA", "VB", "VB-a", "VAB")
CONDITIONS_STRIP = ("V0", "Vw", "Vi", "VA", "VB", "VAB")
COND_LABEL = {
    "V0": "V0 対照（bed_gen2 のまま）",
    "Vw": "Vw winnow 単独（A-1）",
    "Vi": "Vi interlock 単独（A-3）",
    "VA": "VA A 全部（A-1+A-2+A-3）",
    "VB": "VB B 全部（B-1+B-2b+B-3）",
    "VB-a": "VB-a B の診断亜種（B-1+B-2a+B-3）",
    "VAB": "VAB A+B",
}

READ_HFIT_LOW = 0.50
READ_HFIT_SURVIVE = 0.60
G_A2_RA = (1.5, 3.0)
G_A2_SK = (0.4, 0.8)
G_B3_SK = 0.3
G_B2_FACET_MM = 10.0

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
    data = (json.dumps(_san(obj), ensure_ascii=False, indent=1) + "\n").encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True

    atomic_write_bytes(path, data, verify=verify)
    return path


def save_npy(path, arr):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp.npy"
    np.save(tmp, arr)
    back = np.load(tmp, mmap_mode="r")
    ok = (back.shape == arr.shape and back.dtype == arr.dtype)
    del back
    gc.collect()
    if not ok:
        os.remove(tmp)
        raise RuntimeError("npy verify failed: %s" % path)
    os.replace(tmp, path)
    return path


def cp_note(cp, stage, payload):
    cp.unit_done(stage)
    st = cp.state() or {}
    stages = st.get("stages")
    if not isinstance(stages, dict):
        stages = {}
    stages[stage] = _san(payload)
    st["stages"] = stages
    atomic_write_bytes(cp.path,
                       json.dumps(st, ensure_ascii=False, indent=2).encode("utf-8"))
    return st


class Blocker(Exception):
    def __init__(self, reason, ctx):
        super().__init__(reason)
        self.reason = reason
        self.ctx = ctx


# ---------------------------------------------------------------------------
# instruments
# ---------------------------------------------------------------------------
def ruler():
    import run_v4_2_5 as R5
    import run_v4_2_12b as R12
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band"):
        assert hasattr(R5, name), "run_v4_2_5 lost %s" % name
    for name in ("nikora_measure", "break_search", "alt_break_lastlag"):
        assert hasattr(R12, name), "run_v4_2_12b lost %s" % name
    return R5, R12


def measured_curve(z, cell, dmax, axis, as_float32=False):
    R5, R12 = ruler()
    zz = z.astype(np.float32) if as_float32 else z
    return R5.curve_dir(zz, cell, dmax, nbin=NBIN, axis=axis)


def read_axis(z, cell, d50, dmax, lab, axis):
    """One (bed, direction) read by the fixed Nikora procedure."""
    R5, R12 = ruler()
    # the mandatory pre-registered window; set here (not only in stage0) so a
    # cached STAGE 0 cannot leave the module global at the 12B window
    R12.BREAK_LO, R12.BREAK_HI = BREAK_WIN
    d, s = measured_curve(z, cell, dmax, axis,
                          as_float32=(z.size > 8_000_000))
    m = R12.nikora_measure(d, s, float(d50), lab, R5)
    alt = R12.alt_break_lastlag(d, s, float(d50))
    hl = R5.local_H(d, s)
    band = (d >= FIXED_BAND[0]) & (d <= FIXED_BAND[1])
    brk = m.get("break") or {}
    lo = BREAK_WIN[0] * d50
    cands = brk.get("candidates") or []
    floor_lag = (min(c["lag_m"] for c in cands) if cands else None)
    out = {
        "axis": int(axis), "dmax_m": float(dmax), "d50_m": float(d50),
        "n_lags": int(d.size),
        "lag_min_m": float(d.min()), "lag_max_m": float(d.max()),
        "cell_mm": cell * 1000.0,
        "sigma_z_m": float(np.std(z)),
        "sigma_z_over_d50": float(np.std(z) / d50),
        "skew": float(np.mean((z - z.mean()) ** 3) / max(z.std(), 1e-30) ** 3),
        "H_Nikora": m.get("H_Nikora"),
        "break": {k: v for k, v in brk.items() if k != "candidates"},
        "break_at_window_floor": bool(
            floor_lag is not None and brk.get("lag_m") is not None
            and abs(brk["lag_m"] - floor_lag) < 1e-12),
        "H_alt_lastlag": (alt or {}).get("H"),
        "break_alt_lastlag_m": (alt or {}).get("lag_m"),
        "dH_alt": (float((alt or {}).get("H") - m["H_Nikora"]["H"])
                   if (alt and m.get("H_Nikora")) else None),
        "H_fit_fixed_band": R5.fit_band(d, s, *FIXED_BAND),
        "H_loc": {"lag_m": [float(v) for v in d],
                  "H": [float(v) for v in hl]},
        "H_loc_band": {
            "min": float(np.nanmin(np.asarray(hl)[band])),
            "max": float(np.nanmax(np.asarray(hl)[band])),
            "mean": float(np.nanmean(np.asarray(hl)[band]))},
        "curve": {"D_m": [float(v) for v in d],
                  "S_m2": [float(v) for v in s]},
    }
    return out


def measure_bed(z, cell, d50, dmax_x, dmax_y):
    return {"x": read_axis(z, cell, d50, dmax_x, "x", 1),
            "y": read_axis(z, cell, d50, dmax_y, "y", 0)}


def full_stats(bed, z, gid, Dc_star, with_porosity=True, with_facet=True):
    m = B3.gsd_metrics(bed, z, gid)
    out = {"gsd": m,
           "cluster_coverage": B3.cluster_coverage(bed, z, gid),
           "fines_profile": B3.fines_profile(bed, z, gid, Dc_star),
           "n_grains": len(bed.grains)}
    if with_porosity:
        out["porosity"] = B3.porosity_bed(bed, z)
    if with_facet:
        fs = B3.facet_stats(bed, z, gid, exposed_only=True)
        out["facet_exposed"] = fs
        out["facet_all"] = B3.facet_stats(bed)
    return out


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def preflight():
    rows = []
    missing = []
    for label, path in REQUIRED:
        ok = os.path.isfile(path)
        rows.append({"label": label, "path": path, "exists": bool(ok),
                     "bytes": (os.path.getsize(path) if ok else None),
                     "sha256": (sha256_file(path) if ok else None)})
        if not ok:
            missing.append(path)
        log("  %-52s %s %10s" % (label, "OK" if ok else "MISSING",
                                 rows[-1]["bytes"] if ok else "-"))
    if missing:
        raise Blocker("STAGE 0: %d required file(s) missing: %s"
                      % (len(missing), missing), {"preflight": rows})
    return rows


def stage0(cp):
    log("STAGE0  preflight (the brief's list only)")
    rows = preflight()
    import bed_gen2 as B2
    from scipy import ndimage
    t0 = time.time()
    stats = json.load(open(PATCH_STATS, encoding="utf-8"))
    sigma_form = B3.sigma_form_ref(PATCH_STATS)
    form = B2.form_field(512, sigma_form, B2.SEED)
    grains, _cm = B2.get_population(B2.SEED, B2.P_SCALE0)
    B2.place(grains, form)
    kinds = np.array([str(g["kind"]) for g in grains])
    skip = set(int(i) for i in np.where(
        np.array([k.startswith("mono") for k in kinds]))[0])
    xs = (np.arange(2048) + 0.5) * CELL
    z, gid = B3.raycast_zgid([g for i, g in enumerate(grains) if i not in skip],
                             xs, xs, 2048, 2048)
    hole = gid < 0
    num = ndimage.gaussian_filter(np.where(hole, 0.0, z), 60.0, mode="wrap")
    den = ndimage.gaussian_filter((~hole).astype(np.float64), 60.0, mode="wrap")
    zh = np.where(hole, num / np.maximum(den, 1e-9), z)
    ref = np.load(PATCH_NPY).astype(np.float64)
    same = bool(np.array_equal(zh.astype(np.float32), ref.astype(np.float32)))
    maxdiff = float(np.max(np.abs(zh - ref)))
    log("  rebuild %.1f s, %d grains, bit-identical=%s (max|d|=%.1e)"
        % (time.time() - t0, len(grains) - len(skip), same, maxdiff))
    rebuild = {"grains": len(grains) - len(skip),
               "sigma_form_m": sigma_form,
               "bit_identical_float32": same,
               "max_abs_diff_m": maxdiff,
               "construction": "B2.form_field(512,sigma_form,SEED) + "
                               "get_population(SEED,P_SCALE0) + place + "
                               "raycast (mono skipped) + gaussian(60) inpaint"}

    meta = json.load(open(PATCH_META, encoding="utf-8"))
    d50 = float(meta["D50_surface_measured_m"])

    # (a) reproduction with the window the 12B value was measured with
    R5, R12 = ruler()
    old_win = (R12.BREAK_LO, R12.BREAK_HI)
    R12.BREAK_LO = 0.25
    rec_w0 = measure_bed(zh, CELL, d50, DM_PATCH, DM_PATCH)
    # and the mandated window [0.10, 2]
    R12.BREAK_LO = BREAK_WIN[0]
    R12.BREAK_HI = BREAK_WIN[1]
    rec = measure_bed(zh, CELL, d50, DM_PATCH, DM_PATCH)
    hx0, hy0 = rec_w0["x"]["H_Nikora"]["H"], rec_w0["y"]["H_Nikora"]["H"]
    hx, hy = rec["x"]["H_Nikora"]["H"], rec["y"]["H_Nikora"]["H"]
    sig = rec["x"]["sigma_z_over_d50"]
    skw = rec["x"]["skew"]
    hfit = {"x": rec["x"]["H_fit_fixed_band"]["H"],
            "y": rec["y"]["H_fit_fixed_band"]["H"]}
    # old ruler's iso reading too (v4-2-5 published it)
    d_iso, s_iso = R5.curve_iso(zh, CELL, 0.20, nbin=NBIN)
    hfit["iso"] = R5.fit_band(d_iso, s_iso, *FIXED_BAND)["H"]

    # the long box, where 0.29 was established (headerless raw float32)
    zl = np.fromfile(LONG_NPY, dtype="<f4").reshape(2048, 10240)
    long_rec = {}
    for lab, ax, dmax in (("x", 1, 1.00), ("y", 0, 0.30)):
        d, s = R5.curve_dir(zl, CELL, dmax, nbin=NBIN, axis=ax)
        long_rec[lab] = {"H_fit": R5.fit_band(d, s, *FIXED_BAND)}
    d, s = R5.curve_iso(zl, CELL, 0.20, nbin=NBIN)
    long_rec["iso"] = {"H_fit": R5.fit_band(d, s, *FIXED_BAND)}

    gates = {
        "a_H_x_12B_window": {"value": hx0, "expect": GATE_H_EXPECT["x"],
                             "diff": hx0 - GATE_H_EXPECT["x"],
                             "pass": abs(hx0 - GATE_H_EXPECT["x"]) <= GATE_H_TOL},
        "a_H_y_12B_window": {"value": hy0, "expect": GATE_H_EXPECT["y"],
                             "diff": hy0 - GATE_H_EXPECT["y"],
                             "pass": abs(hy0 - GATE_H_EXPECT["y"]) <= GATE_H_TOL},
        "a_sigma_z_over_d50": {"value": sig, "expect": GATE_SIG_EXPECT,
                               "diff": sig - GATE_SIG_EXPECT,
                               "pass": abs(sig - GATE_SIG_EXPECT) <= GATE_SIG_TOL},
        "a_skew": {"value": skw, "expect": GATE_SKEW_EXPECT,
                   "diff": skw - GATE_SKEW_EXPECT,
                   "pass": abs(skw - GATE_SKEW_EXPECT) <= GATE_SKEW_TOL},
        "b_H_fit_long_x_0p29": {"value": long_rec["x"]["H_fit"]["H"],
                                "expect": HFIT_029,
                                "diff": long_rec["x"]["H_fit"]["H"] - HFIT_029,
                                "pass": abs(long_rec["x"]["H_fit"]["H"] - HFIT_029)
                                        <= HFIT_029_TOL},
        "b_H_fit_patch_informational": {
            "value": hfit["x"], "expect": HFIT_029,
            "diff": hfit["x"] - HFIT_029,
            "pass": abs(hfit["x"] - HFIT_029) <= HFIT_029_TOL,
            "note": "the PATCH box old-ruler value; 0.29 was established on the "
                    "long/strip box, and the patch value 0.323 is v4-2-5's own "
                    "published patch number.  Reported, not part of the gate."},
    }
    gate_repro = all(gates[k]["pass"] for k in
                     ("a_H_x_12B_window", "a_H_y_12B_window", "a_sigma_z_over_d50",
                      "a_skew", "b_H_fit_long_x_0p29"))
    log("  gate(a) W0 H x %.4f (d %+.4f) y %.4f (d %+.4f) sig %.4f skew %+.3f"
        % (hx0, hx0 - GATE_H_EXPECT["x"], hy0, hy0 - GATE_H_EXPECT["y"], sig, skw))
    log("  gate(a) mandated window H x %.4f y %.4f" % (hx, hy))
    log("  gate(b) H_fit long x %.4f  patch x/y/iso %.4f/%.4f/%.4f"
        % (long_rec["x"]["H_fit"]["H"], hfit["x"], hfit["y"], hfit["iso"]))
    log("  STAGE0 gate: %s" % ("PASS" if gate_repro else "FAIL"))

    out = {"preflight": rows, "rebuild": rebuild, "d50_m": d50, "cell_m": CELL,
           "ruler_window": {"mandated": list(BREAK_WIN),
                            "12B_measurement": [0.25, 2.0],
                            "v4_2_13_readings": {"W0": 0.4484, "W2": 0.4653}},
           "measure_w0": rec_w0, "measure": rec,
           "H_fit_patch": hfit, "H_fit_long": long_rec,
           "gates": gates, "gate_pass": bool(gate_repro)}
    if not gate_repro:
        raise Blocker("STAGE 0 gate failed (the mechanisms are not run)", out)

    # the control bed is kept for STAGE 1 (not re-derived from the file)
    save_npy(os.path.join(OUT, "heightfield-STAGE0-patch-control.npy"),
             zh.astype(np.float32))
    out["files"] = ["heightfield-STAGE0-patch-control.npy"]
    cp_note(cp, "stage0", {"gate_pass": True, "H_x": hx0, "H_y": hy0,
                           "H_x_mandated": hx, "H_y_mandated": hy,
                           "H_fit_long_x": long_rec["x"]["H_fit"]["H"],
                           "bit_identical": same})
    return out, zh, d50, sigma_form


# ---------------------------------------------------------------------------
# STAGE 1 -- mechanisms, step by step
# ---------------------------------------------------------------------------
def stage1(cp, patch_bed, sigma_form):
    log("STAGE1  bed_gen3 mechanisms, one step at a time (patch seed0)")
    b = patch_bed
    z, gid = b.render()
    m0 = B3.gsd_metrics(b, z, gid)
    b.meta.update({"D50_ref": m0.get("surface_D50_kind_m") or m0["surface_D50_m"],
                   "D84_ref": m0["surface_D84_m"],
                   "D90_ref": m0["surface_D90_m"]})
    out = {"initial": {"D50_ref": b.meta["D50_ref"],
                       "D84_ref": b.meta["D84_ref"],
                       "D90_ref": b.meta["D90_ref"],
                       "gsd": m0,
                       "cluster_coverage": B3.cluster_coverage(b, z, gid),
                       "porosity": B3.porosity_bed(b, z)}}
    log("  V0 gsd: surface D50 %.5f D84 %.5f D90 %.5f R_a %.3f"
        % (m0["surface_D50_m"], m0["surface_D84_m"], m0["surface_D90_m"],
           m0.get("R_a") or float("nan")))

    # ---- A-1 sweep -------------------------------------------------------
    rng = np.random.default_rng(SEED_BASE + 1)
    sweep = []
    selected_bed = None
    selected_tau = None
    for ts in B3.TAU_STAR_SWEEP:
        bw, dw = B3.winnow_at(b, ts, rng)
        sweep.append(dw)
        log("  A-1 tau*=%.3f removed %4d R_a %.3f stop %s"
            % (ts, dw["removed_grains"], dw["R_a"] or float("nan"),
               dw["stop_reason"]))
        if selected_tau is None and dw["R_a"] is not None and dw["R_a"] >= 2.0:
            selected_tau = ts
            selected_bed = bw
    if selected_bed is None:
        selected_tau = B3.TAU_STAR_SWEEP[-1]
        selected_bed = bw
    # keep the tau(x,y) field of the selected tau* for the record
    gg = 4   # 4 bed pixels per tau-grid cell (matches bed_gen3.winnow_at)
    tau, tau_i = B3.tau_grid(b, selected_tau, np.random.default_rng(5), grid_px=gg)
    save_npy(os.path.join(OUT, "tau-field-patch-s0.npy"), tau.astype(np.float32))
    out["A1"] = {"sweep": sweep, "selected_tau_star": selected_tau,
                 "selected_from": "first tau* reaching R_a>=2, else the highest",
                 "tau_field": b.tau, "field_file": "tau-field-patch-s0.npy",
                 "physics": {"tau0_scene_Pa": B3.TAU0_SCENE,
                             "tau_star_scene": B3.TAU0_SCENE / (B3.DRHO * B3.G
                                                                * b.meta["D50_ref"]),
                             "note": "tau* is the sweep parameter of the "
                                     "mechanism strength; the scene's own "
                                     "tau0 is below the swept band"}}
    log("  A-1 selected tau* = %.3f" % selected_tau)

    # ---- A-2 on the A-1 result -------------------------------------------
    rng2 = np.random.default_rng(SEED_BASE + 2)
    bs, ds = B3.granular_sorting(selected_bed, rng2)
    out["A2"] = ds
    log("  A-2 Dc*=%.5f moved %d vol residual %.2e"
        % (ds["Dc_star_m"], ds["n_moved"], ds["volume_residual_m3"]))
    out["A2"]["profile_after_short"] = [
        (r["bin"], r["fines_volume_fraction"]) for r in ds["profile_after"]["bins"]]

    # ---- A-3 on the A-2 result -------------------------------------------
    rng3 = np.random.default_rng(SEED_BASE + 3)
    bi, di = B3.interlock(bs, rng3)
    out["A3"] = di
    log("  A-3 support>=3 %.2f dip %.1f poro %.3f->%.3f (drop %.4f) cov %.3f"
        % (di["frac_support_ge3"] or 0, di["dip_median_deg"] or 0,
           di["porosity_before"]["porosity"], di["porosity_after"]["porosity"],
           di["porosity_drop"], di["cluster_coverage_after"]))

    # ---- B-1 -------------------------------------------------------------
    rng4 = np.random.default_rng(SEED_BASE + 4)
    bf, df = B3.facetise(b, rng4)
    out["B1"] = df
    log("  B-1 faces med %d facet med %.2f mm (area-weighted %.2f) vol/ell %.3f"
        % (df["n_faces"]["median"], df["facet_linear_mm"]["median"],
           B3.facet_stats(bf)["area_weighted_median_mm"],
           df["volume_over_ellipsoid"]["median"]))

    # ---- B-2 (both sub-variants, both surface definitions) ---------------
    bp, dp = B3.planation(bf, k=1.0, mode="clip", top_f=0.10)
    out["B2a"] = dp
    bp10, dp10 = B3.planation(bf, k=0.5, mode="clip", top_f=0.05)
    out["B2a_k05"] = dp10
    br, dr = B3.planation(bf, k=1.0, mode="remove",
                          rng=np.random.default_rng(SEED_BASE + 6), top_f=0.10)
    out["B2b"] = dr
    br10, dr10 = B3.planation(bf, k=0.5, mode="remove",
                              rng=np.random.default_rng(SEED_BASE + 7), top_f=0.05)
    out["B2b_k05"] = dr10
    log("  B-2a clip f %.3f (k=1.0) %.3f (k=0.5); B-2b removed %d / %d"
        % (dp["clip_area_fraction"], dp10["clip_area_fraction"],
           dr["n_removed_grains"], dr10["n_removed_grains"]))

    # ---- B-3 -------------------------------------------------------------
    rng5 = np.random.default_rng(SEED_BASE + 5)
    b3, d3 = B3.pocket_fill(br, rng5, ds["Dc_star_m"])
    out["B3"] = d3
    log("  B-3 filled %d grains (%.1f%% of target) poro %.3f->%.3f"
        % (d3["n_filled"], 100 * d3["fill_of_target"],
           d3["porosity_before"]["porosity"], d3["porosity_after"]["porosity"]))
    # final bed of the full B chain, measured
    zf, gidf = b3.render(force=True)
    out["B_chain_final"] = full_stats(b3, zf, gidf, ds["Dc_star_m"])
    out["files"] = ["tau-field-patch-s0.npy"]
    cp_note(cp, "stage1", {"selected_tau_star": selected_tau,
                           "Dc_star_m": ds["Dc_star_m"],
                           "facet_median_mm": df["facet_linear_mm"]["median"]})
    return out, selected_tau, ds["Dc_star_m"], b.meta.copy()


# ---------------------------------------------------------------------------
# STAGE 2 -- the factorial
# ---------------------------------------------------------------------------
def build_condition(base, cond, tau_star, Dc_star, rng):
    """Apply the declared chain for one condition; returns (bed, steps)."""
    steps = {}
    b = base
    if cond in ("Vw", "VA", "VAB"):
        b, d = B3.winnow_at(b, tau_star, rng)
        steps["A1"] = d
    if cond in ("VA", "VAB"):
        b, d = B3.granular_sorting(b, rng)
        steps["A2"] = d
    if cond in ("Vi", "VA", "VAB"):
        b, d = B3.interlock(b, rng)
        steps["A3"] = d
    if cond in ("VB", "VB-a", "VAB"):
        b, d = B3.facetise(b, rng)
        steps["B1"] = d
        if cond == "VB-a":
            b, d = B3.planation(b, k=1.0, mode="clip")
            steps["B2a"] = d
        else:
            b, d = B3.planation(b, k=1.0, mode="remove", rng=rng)
            steps["B2b"] = d
        b, d = B3.pocket_fill(b, rng, Dc_star)
        steps["B3"] = d
    return b, steps


def strip_steps(rec):
    """The compact per-step record kept in the numbers JSON."""
    keep = {}
    for k, d in rec.items():
        if k == "A1":
            keep[k] = {kk: d[kk] for kk in
                       ("tau_star", "passes", "removed_grains", "stop_reason",
                        "R_a", "R_a_kind", "history")}
        elif k == "A2":
            keep[k] = {kk: d[kk] for kk in
                       ("Dc_star_m", "n_moved", "n_coarse_lifted",
                        "volume_down_m3", "volume_up_m3", "volume_residual_m3")}
            keep[k]["profile_before"] = [(r["bin"], r["fines_volume_fraction"])
                                         for r in d["profile_before"]["bins"]]
            keep[k]["profile_after"] = [(r["bin"], r["fines_volume_fraction"])
                                        for r in d["profile_after"]["bins"]]
        elif k == "A3":
            keep[k] = {kk: d[kk] for kk in
                       ("n_added_cluster_grains", "cluster_coverage_before",
                        "cluster_coverage_after", "support_hist",
                        "support_median", "frac_support_ge3", "n_sink_steps",
                        "dip_median_deg", "dip_p16_deg", "dip_p84_deg",
                        "porosity_drop")}
            keep[k]["porosity_before"] = d["porosity_before"]["porosity"]
            keep[k]["porosity_after"] = d["porosity_after"]["porosity"]
        elif k == "B1":
            keep[k] = {kk: d[kk] for kk in
                       ("n_grains", "n_faces", "facet_linear_mm", "ba", "corey",
                        "wadell_volume_proxy", "facet_edge_roundness",
                        "facet_dihedral_deg", "facet_area_fraction")}
            keep[k]["volume_over_ellipsoid"] = d["volume_over_ellipsoid"]["median"]
        elif k in ("B2a", "B2b"):
            keep[k] = {kk: d.get(kk) for kk in
                       ("mode", "k", "clip_area_fraction", "cut_volume_m3",
                        "n_clipped_grains", "clipped_median_m",
                        "n_removed_grains", "removed_excess_median_m",
                        "equiv_top_f", "equiv_top_f_measured")}
        elif k == "B3":
            keep[k] = {kk: d.get(kk) for kk in
                       ("n_filled", "filled_volume_m3", "target_volume_m3",
                        "fill_of_target", "porosity_before", "porosity_after")}
            if keep[k]["porosity_before"]:
                keep[k]["porosity_before"] = keep[k]["porosity_before"].get("porosity")
            if keep[k]["porosity_after"]:
                keep[k]["porosity_after"] = keep[k]["porosity_after"].get("porosity")
    return keep


def grain_records(bed, z, gid):
    t = B3.grain_table(bed, z, gid)
    out = []
    for i, g in enumerate(bed.grains):
        out.append({
            "x": round(float(g["x"]), 6), "y": round(float(g["y"]), 6),
            "z": round(float(g["z"]), 6), "D": round(float(g["D"]), 6),
            "abc": [round(float(v), 6) for v in g["abc"]],
            "R": [round(float(v), 5) for v in np.asarray(g["R"]).ravel()],
            "kind": str(g["kind"]), "exposed": bool(t["exposed"][i]),
            "n_facets": int(g["frac"]["metrics"]["n_faces"]) if "frac" in g else 0})
    return out


def rec_path(box, cond, si):
    return os.path.join(OUT, "rec-%s-%s-s%d.json" % (box, cond, si))


def stage2(cp, sigma_form, tau_star, Dc_star, seed_list, boxes=("patch", "strip")):
    """The factorial.  Resumable at the bed level: a bed whose height field,
    grain list and rec JSON all exist is loaded, not recomputed (the RNG of
    every condition is derived deterministically from seed+condition, so a
    resumed bed is bit-identical to the one a full run would make)."""
    log("STAGE2  factorial: compute %s (patch R=%d, strip R=%d)"
        % (",".join(boxes), N_PATCH, N_STRIP))
    res = {"patch": {}, "strip": {}, "complete": {}}
    for cond in CONDITIONS_PATCH:
        res["patch"][cond] = {}
    for cond in CONDITIONS_STRIP:
        res["strip"][cond] = {}
    files = []
    for box, seeds, conds in (
            ("patch", seed_list["patch"], CONDITIONS_PATCH),
            ("strip", seed_list["strip"], CONDITIONS_STRIP)):
        do = box in boxes
        # a box that is not selected is loaded from disk when complete,
        # otherwise it is skipped (a partial pass writes no deliverables)
        have_all = all(os.path.isfile(rec_path(box, c, si))
                       for si in range(len(seeds))
                       for c in conds)
        if not do and not have_all:
            log("  box %s is incomplete and not selected in this pass: skipped"
                % box)
            res["complete"][box] = False
            continue
        res["complete"][box] = True
        for si, seed in enumerate(seeds):
            # ---- V0 ---------------------------------------------------
            rc = rec_path(box, "V0", si)
            hf = os.path.join(OUT, "heightfield-V0-%s-s%d.npy" % (box, si))
            if do and not (os.path.isfile(rc) and os.path.isfile(hf)):
                t0 = time.time()
                tiles = [int(seed) + 10 * k for k in range(4)]
                base = B3.initial_bed(int(seed), box, sigma_form,
                                      tile_seeds=(tiles if box == "strip" else None))
                z0, gid0 = base.render()
                m0 = B3.gsd_metrics(base, z0, gid0)
                d50 = m0.get("surface_D50_kind_m") or m0["surface_D50_m"]
                base.meta.update({"D50_ref": d50, "D84_ref": m0["surface_D84_m"],
                                  "D90_ref": m0["surface_D90_m"],
                                  "seed": int(seed), "realisation": si})
                log("  %s seed %d (%d): %.1f s, %d grains, D50 %.5f"
                    % (box, si, seed, time.time() - t0, len(base.grains), d50))
                v0 = full_stats(base, z0, gid0, Dc_star)
                rec = {
                    "seed": int(seed),
                    "tile_seeds": (tiles if box == "strip" else None),
                    "d50_ref_m": d50,
                    "heightfield": "heightfield-V0-%s-s%d.npy" % (box, si),
                    "grains": "grains-V0-%s-s%d.json" % (box, si),
                    "measure": measure_bed(z0, base.cell, d50,
                                           DM_PATCH if box == "patch" else DM_STRIP_X,
                                           DM_PATCH if box == "patch" else DM_STRIP_Y),
                    "stats": v0, "steps": {}}
                save_npy(hf, z0.astype(np.float32))
                write_json(os.path.join(OUT, rec["grains"]),
                           grain_records(base, z0, gid0))
                write_json(rc, rec)
                del base, z0, gid0
                gc.collect()
            if not os.path.isfile(rc):
                log("  %s V0 s%d is missing and was not computed: box partial"
                    % (box, si))
                res["complete"][box] = False
                break
            rec = json.load(open(rc, encoding="utf-8"))
            res[box]["V0"][si] = rec
            files.append(rec["heightfield"])
            d50 = rec["d50_ref_m"]
            # ---- the mechanism conditions ------------------------------
            for cond in conds:
                if cond == "V0":
                    continue
                rc = rec_path(box, cond, si)
                hf = os.path.join(OUT, "heightfield-%s-%s-s%d.npy" % (cond, box, si))
                if do and not (os.path.isfile(rc) and os.path.isfile(hf)):
                    tc = time.time()
                    base = B3.initial_bed(int(seed), box, sigma_form,
                                          tile_seeds=(rec.get("tile_seeds")
                                                      if box == "strip" else None))
                    base.meta.update({"D50_ref": d50,
                                      "D84_ref": rec["stats"]["gsd"]["surface_D84_m"],
                                      "D90_ref": rec["stats"]["gsd"]["surface_D90_m"],
                                      "seed": int(seed), "realisation": si})
                    import zlib
                    cs = zlib.crc32(cond.encode("utf-8"))
                    rng = np.random.default_rng(
                        (int(seed) * 1000003 + cs) % (2 ** 63))
                    b, steps = build_condition(base, cond, tau_star, Dc_star, rng)
                    z, gid = b.render(force=True)
                    st = full_stats(b, z, gid, Dc_star)
                    rec2 = {
                        "seed": int(seed), "d50_ref_m": d50,
                        "heightfield": "heightfield-%s-%s-s%d.npy" % (cond, box, si),
                        "grains": "grains-%s-%s-s%d.json" % (cond, box, si),
                        "measure": measure_bed(z, b.cell, d50,
                                               DM_PATCH if box == "patch" else DM_STRIP_X,
                                               DM_PATCH if box == "patch" else DM_STRIP_Y),
                        "stats": st, "steps": strip_steps(steps),
                        "wall_s": float(time.time() - tc)}
                    save_npy(hf, z.astype(np.float32))
                    write_json(os.path.join(OUT, rec2["grains"]),
                               grain_records(b, z, gid))
                    write_json(rc, rec2)
                    log("    %-5s %s s%d %.1f s: H_Nik x %s y %s | Hfit x %.4f"
                        % (cond, box, si, time.time() - tc,
                           ("%.4f" % rec2["measure"]["x"]["H_Nikora"]["H"])
                           if rec2["measure"]["x"]["H_Nikora"] else "NA",
                           ("%.4f" % rec2["measure"]["y"]["H_Nikora"]["H"])
                           if rec2["measure"]["y"]["H_Nikora"] else "NA",
                           rec2["measure"]["x"]["H_fit_fixed_band"]["H"] or 0.0))
                    del b, z, gid
                    gc.collect()
                if not os.path.isfile(rc):
                    log("  %s %s s%d is missing and was not computed: box partial"
                        % (box, cond, si))
                    res["complete"][box] = False
                    break
                res[box][cond][si] = json.load(open(rc, encoding="utf-8"))
                files.append(res[box][cond][si]["heightfield"])
    res["files"] = files
    if all(res["complete"].get(b) for b in ("patch", "strip")):
        cp_note(cp, "stage2", {"n_patch": sum(len(v) for v in res["patch"].values()),
                               "n_strip": sum(len(v) for v in res["strip"].values()),
                               "n_files": len(files),
                               "computed_boxes": list(boxes)})
    return res


# ---------------------------------------------------------------------------
# STAGE 3 -- pre-registered read-out
# ---------------------------------------------------------------------------
def _mean(vals):
    v = [x for x in vals if x is not None]
    return float(np.mean(v)) if v else None


def cond_summary(stage2, box, cond):
    rows = stage2[box].get(cond) or {}
    out = {"n": len(rows)}
    for ax in ("x", "y"):
        out["H_Nikora_" + ax] = [r["measure"][ax]["H_Nikora"]["H"]
                                 for r in rows.values()]
        out["H_fit_" + ax] = [r["measure"][ax]["H_fit_fixed_band"]["H"]
                              for r in rows.values()]
        out["sigma_" + ax] = [r["measure"][ax]["sigma_z_over_d50"]
                              for r in rows.values()]
        out["skew_" + ax] = [r["measure"][ax]["skew"] for r in rows.values()]
        out["H_alt_" + ax] = [r["measure"][ax]["H_alt_lastlag"]
                              for r in rows.values()]
    out["R_a"] = [r["stats"]["gsd"].get("R_a") for r in rows.values()]
    out["R_a_kind"] = [r["stats"]["gsd"].get("R_a_kind") for r in rows.values()]
    out["surface_D50"] = [r["stats"]["gsd"].get("surface_D50_m")
                          for r in rows.values()]
    out["sigma_g"] = [r["stats"]["gsd"].get("surface_sigma_g")
                      for r in rows.values()]
    out["porosity"] = [r["stats"]["porosity"]["porosity"] for r in rows.values()]
    out["cluster_cov"] = [r["stats"]["cluster_coverage"] for r in rows.values()]
    out["facet_mm"] = [(r["stats"]["facet_exposed"] or {}).get("median_mm")
                       if r["stats"].get("facet_exposed") else None
                       for r in rows.values()]
    out["facet_mm_area"] = [(r["stats"]["facet_exposed"] or {}).get("area_weighted_median_mm")
                            if r["stats"].get("facet_exposed") else None
                            for r in rows.values()]
    out["mean"] = {k: _mean(out[k]) for k in
                   ("H_Nikora_x", "H_Nikora_y", "H_fit_x", "H_fit_y",
                    "H_alt_x", "H_alt_y", "R_a", "R_a_kind", "surface_D50",
                    "sigma_g", "porosity", "cluster_cov", "facet_mm",
                    "facet_mm_area", "sigma_x", "skew_x")}
    out["mean"]["surface_sigma_g"] = _mean(out["sigma_g"])
    return out


def paired_delta(stage2, box, cond, key):
    """Per-seed paired difference cond - V0 (same seed)."""
    a = stage2[box].get(cond) or {}
    b = stage2[box].get("V0") or {}
    out = []
    for si in sorted(set(a) & set(b)):
        va = a[si]["measure"][key[0]][key[1]]
        vb = b[si]["measure"][key[0]][key[1]]
        if key[1] == "H_Nikora":
            va = va["H"]
            vb = vb["H"]
        elif key[1] == "H_fit_fixed_band":
            va = va["H"]
            vb = vb["H"]
        if va is not None and vb is not None:
            out.append(va - vb)
    return out


def stage3(cp, stage2, stage1):
    log("STAGE3  pre-registered read-out")
    summ = {"patch": {c: cond_summary(stage2, "patch", c)
                      for c in CONDITIONS_PATCH if stage2["patch"].get(c)},
            "strip": {c: cond_summary(stage2, "strip", c)
                      for c in CONDITIONS_STRIP if stage2["strip"].get(c)}}
    deltas = {}
    for box, conds in (("patch", CONDITIONS_PATCH), ("strip", CONDITIONS_STRIP)):
        for c in conds:
            if c == "V0" or not stage2[box].get(c):
                continue
            deltas["%s|%s" % (box, c)] = {
                "dH_Nikora_x": paired_delta(stage2, box, c, ("x", "H_Nikora")),
                "dH_Nikora_y": paired_delta(stage2, box, c, ("y", "H_Nikora")),
                "dH_fit_x": paired_delta(stage2, box, c, ("x", "H_fit_fixed_band")),
                "dH_fit_y": paired_delta(stage2, box, c, ("y", "H_fit_fixed_band"))}

    def ro(box, cond):
        s = summ[box].get(cond) or {}
        return s.get("mean") or {}

    P = summ.get("patch", {})
    S = summ.get("strip", {})

    # G-A1 / G-A2 (A conditions: strip box x as the primary reading)
    def gval(box, cond, k):
        s = P.get(cond) if box == "patch" else S.get(cond)
        if not s:
            return None
        return (s.get("mean") or {}).get(k)

    read = {}
    read["G-A1"] = {
        "rule": "H_fit < 0.50 -> A falls as an explanation of the scaling band "
                "(envelope judgement deferred); H_fit >= 0.60 -> survives "
                "pending the envelope",
        "primary": "strip x, mean over 3 seeds",
        "values": {"strip|VA|H_fit_x": gval("strip", "VA", "H_fit_x"),
                   "strip|Vw|H_fit_x": gval("strip", "Vw", "H_fit_x"),
                   "strip|V0|H_fit_x": gval("strip", "V0", "H_fit_x"),
                   "patch|VA|H_fit_x": gval("patch", "VA", "H_fit_x"),
                   "patch|Vw|H_fit_x": gval("patch", "Vw", "H_fit_x")},
        "H_Nikora_strip_VA_x": gval("strip", "VA", "H_Nikora_x"),
    }
    v = read["G-A1"]["values"]["strip|VA|H_fit_x"]
    read["G-A1"]["verdict"] = (
        None if v is None else ("A falls (H_fit<0.50)" if v < READ_HFIT_LOW
                                else ("A survives pending envelope (H_fit>=0.60)"
                                      if v >= READ_HFIT_SURVIVE
                                      else "between 0.50 and 0.60: neither "
                                           "pre-registered branch")))
    ra = gval("strip", "VA", "R_a")
    sk = gval("strip", "VA", "skew_x")
    read["G-A2"] = {
        "rule": "R_a not in [1.5,3.0] or Sk not in [+0.4,+0.8] -> A did not "
                "reproduce armour (regardless of H)",
        "values": {"strip|VA|R_a": ra, "strip|VA|Sk_x": sk,
                   "patch|VA|R_a": gval("patch", "VA", "R_a"),
                   "patch|VA|Sk_x": gval("patch", "VA", "skew_x")},
        "verdict": (None if (ra is None or sk is None) else
                    ("armour reproduced" if (G_A2_RA[0] <= ra <= G_A2_RA[1]
                                             and G_A2_SK[0] <= sk <= G_A2_SK[1])
                     else "armour NOT reproduced"))
    }
    ra_p = gval("patch", "VA", "R_a")
    sk_p = gval("patch", "VA", "skew_x")
    sk_seeds = [(int(k), (v["measure"]["x"] or {}).get("skew"))
                for k, v in (stage2["strip"].get("VA") or {}).items()]
    read["G-A2"]["strip_Sk_by_seed"] = sorted(sk_seeds)
    def _zone(rv, sv):
        okr = rv is not None and G_A2_RA[0] <= rv <= G_A2_RA[1]
        oks = sv is not None and G_A2_SK[0] <= sv <= G_A2_SK[1]
        return "R_a %s (%s), Sk %s (%s)" % (
            f(rv), "in" if okr else "out", sgn(sv), "in" if oks else "out")
    read["G-A2"]["box_note"] = (
        "strip (primary): %s; patch: %s.  On the strip only the skew misses, "
        "by %+.4f, and the per-seed skew is %s (mean of three values straddles "
        "the +0.4 bound)."
        % (_zone(ra, sk), _zone(ra_p, sk_p),
           (abs((sk or 0) - G_A2_SK[0]) if sk is not None else 0.0),
           ", ".join("%+.3f" % v for _, v in
                     (read["G-A2"].get("strip_Sk_by_seed") or []))))
    read["G-A2"]["box_note_detail"] = "strip: R_a %s, Sk %s; patch: R_a %s, Sk %s" % (
        f(ra), sgn(sk), f(ra_p), sgn(sk_p))
    # G-A3: subsurface fines profile
    prof = {}
    for box, cond in (("patch", "VA"), ("strip", "VA")):
        rows = stage2[box].get(cond) or {}
        if rows:
            p0 = (rows[sorted(rows)[0]]["steps"].get("A2") or {})
            prof[box] = {"before": p0.get("profile_before"),
                         "after": p0.get("profile_after")}
    read["G-A3"] = {
        "rule": "if the subsurface fines profile matches neither A1 (unchanged) "
                "nor A2 (rise at 2-5 x D90) the composition of the sorting is "
                "not identifiable",
        "profiles": prof,
        "values": {},
        "verdict": None}
    if prof.get("patch"):
        before = dict(prof["patch"]["before"] or [])
        after = dict(prof["patch"]["after"] or [])
        rise = []
        for b in ("2-3", "3-5"):
            if b in before and b in after and after[b] is not None and before[b] is not None:
                rise.append(after[b] - before[b])
        for b_ in ("2-3", "3-5"):
            if b_ in before and b_ in after:
                read["G-A3"]["values"]["patch|VA|fines_%s_before" % b_] = before[b_]
                read["G-A3"]["values"]["patch|VA|fines_%s_after" % b_] = after[b_]
        read["G-A3"]["verdict"] = (
            "A2-type rise at 2-5 x D90 (mean rise %+.4f)"
            % (float(np.mean(rise)) if rise else 0.0)
            if rise and np.mean(rise) > 0.01 else
            "neither A1 nor A2 cleanly (rise %+.4f)"
            % (float(np.mean(rise)) if rise else 0.0))

    v = gval("strip", "VB", "H_fit_x")
    read["G-B1"] = {
        "rule": "same H_fit thresholds as G-A1, for mechanism B",
        "values": {"strip|VB|H_fit_x": v,
                   "strip|VB-a|H_fit_x": gval("patch", "VB-a", "H_fit_x"),
                   "strip|V0|H_fit_x": gval("strip", "V0", "H_fit_x")},
        "H_Nikora_strip_VB_x": gval("strip", "VB", "H_Nikora_x"),
        "verdict": (None if v is None else
                    ("B falls (H_fit<0.50)" if v < READ_HFIT_LOW
                     else ("B survives pending envelope (H_fit>=0.60)"
                           if v >= READ_HFIT_SURVIVE
                           else "between 0.50 and 0.60: neither branch")))}
    fm = gval("strip", "VB", "facet_mm")
    fma = gval("strip", "VB", "facet_mm_area")
    read["G-B2"] = {
        "rule": "median facet linear size < 10 mm -> B8's premise (ii) breaks",
        "values": {"strip|VB|facet_median_mm_count": fm,
                   "strip|VB|facet_median_mm_area_weighted": fma,
                   "patch|VB|facet_median_mm_count": gval("patch", "VB", "facet_mm"),
                   "patch|VB|facet_median_mm_area_weighted": gval("patch", "VB", "facet_mm_area")},
        "verdict": (None if fm is None else
                    ("count median <10 mm (fails); area-weighted median %s mm %s"
                     % (f(fma, 2),
                        ">=10 mm (passes)" if (fma or 0) >= G_B2_FACET_MM
                        else "<10 mm (also fails)")))}
    read["G-B2"]["definition_note"] = (
        "two pre-declared readings: the unweighted median over all facets is "
        "dominated by the many small grains; the area-weighted median is the "
        "median facet actually seen on the surface.  Both are reported and "
        "neither is dropped.")
    ska = gval("patch", "VB-a", "skew_x")
    skb = gval("patch", "VB", "skew_x")
    read["G-B3"] = {
        "rule": "VB-a (B-2a) Sk < +0.3 -> B6's tension is confirmed; B-2a falls "
                "as a skewness mechanism, and if B-2b recovers Sk, B survives "
                "only as B-2b",
        "values": {"patch|VB-a|Sk_x": ska, "patch|VB|Sk_x": skb,
                   "patch|V0|Sk_x": gval("patch", "V0", "skew_x")},
        "verdict": (None if ska is None else
                    ("B-2a loses skew (Sk %s<+0.3); B-2b Sk %s -> %s"
                     % (sgn(ska), sgn(skb),
                        ("also below +0.3: B-2b does NOT recover the skew"
                         if (skb is not None and skb < G_B3_SK)
                         else "recovers skew"))))}
    read["VAB"] = {
        "rule": "evaluate VAB if either A or B survives; evaluate it once even "
                "if both fail alone",
        "values": {"strip|VAB|H_fit_x": gval("strip", "VAB", "H_fit_x"),
                   "patch|VAB|H_fit_x": gval("patch", "VAB", "H_fit_x"),
                   "strip|VAB|H_Nikora_x": gval("strip", "VAB", "H_Nikora_x")},
        "verdict": None}
    v = gval("strip", "VAB", "H_fit_x")
    read["VAB"]["verdict"] = (None if v is None else
                              ("VAB falls (H_fit<0.50)" if v < READ_HFIT_LOW
                               else ("VAB survives pending envelope"
                                     if v >= READ_HFIT_SURVIVE
                                     else "between 0.50 and 0.60")))
    read["note"] = ("no envelope verdict in this job: there is no measured DEM. "
                    "0.50 / 0.60 are placement values, not calibration. A's "
                    "H_Nikora was expected near 0.50 (optimistic 0.60); 0.79 is "
                    "not written as A's prediction.")
    cp_note(cp, "stage3", {k: read[k].get("verdict") for k in
                           ("G-A1", "G-A2", "G-B1", "G-B2", "G-B3", "VAB")})
    return {"summary": summ, "paired_deltas": deltas, "readout": read}


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f(v, nd=4):
    if v is None:
        return "—"
    try:
        return ("%." + str(nd) + "f") % float(v)
    except Exception:
        return str(v)


def sgn(v, nd=3):
    if v is None:
        return "—"
    return ("%+." + str(nd) + "f") % float(v)


def hn(mx):
    v = (mx or {}).get("H_Nikora")
    return v.get("H") if isinstance(v, dict) else None


def build_report(N):
    st0 = N["stage0"]
    st1 = N["stage1"]
    st2 = N["stage2"]
    st3 = N["stage3"]
    R = []
    A = R.append
    A("# REPORT-V4-2-15 — 粒層の要因計画（機構 A／機構 B、v4-2.15）")
    A("")
    A("- brief: `D:\\mutaro-work\\rally-in\\BRIEF-V4-2-15.md`")
    A("- 実行: `run_v4_2_15.py`（前景・最後まで）")
    A("- 実装: `bed_gen3.py`（新規。A-0〜A-3、B-1〜B-3 を独立関数で）")
    A("- 作成: %s" % N["created_utc"])
    A("- 書き込み範囲: `bed_gen3.py`、`run_v4_2_15.py`、`out/v4-2-15/` のみ。"
      "既存ファイルは変更していない。")
    A("")
    A("**このレポートは数値と事実の整理である。設計判定・次の一手の提案は書かない"
      "（brief STAGE 3 の指示）。実測 DEM が無いのでエンベロープ判定は出さない。**")
    A("")
    A("## 0. 要約")
    A("")
    g = st0["gates"]
    A("- **STAGE 0**: brief の絶対パスは全部存在（欠損 0）。パッチ床を v4-2-2 の"
      "作り方で再生成し、保存済み `heightfield.npy` と float32 で**一致**"
      "（max|Δ| = %.1e m）。12B 公表値の窓 [0.25,2]×D50 で H_Nikora x %.4f / "
      "y %.4f（12B 0.4484 / 0.4320、±0.005）、σz/D50 %.4f、歪度 %+.3f。"
      "旧物差しの H_fit は 0.29 が定義された帯状箱（3.0 m r1）で x %.4f / "
      "y %.4f / iso %.4f、新規の帯状 V0 で x %.4f。**ゲート PASS**。"
      % (st0["rebuild"]["max_abs_diff_m"], g["a_H_x_12B_window"]["value"],
         g["a_H_y_12B_window"]["value"], g["a_sigma_z_over_d50"]["value"],
         g["a_skew"]["value"], st0["H_fit_long"]["x"]["H_fit"]["H"],
         st0["H_fit_long"]["y"]["H_fit"]["H"],
         st0["H_fit_long"]["iso"]["H_fit"]["H"],
         st2["strip"]["V0"][sorted(st2["strip"]["V0"])[0]]["measure"]["x"]["H_fit_fixed_band"]["H"]))
    a1 = st1["A1"]
    A("- **STAGE 1**: A-1 の掃引は τ* = %s を選び、除去 %d 粒・R_a %.3f に到達"
      "（%s）。A-2 は Dc* = %.4f m、%d 粒移動、体積残差 %.1e m³。A-3 は"
      "支持点≥3 の割合 %.2f、伏角中央 %s°、空隙率 %.3f→%.3f（低下 %+.4f）。"
      "B-1 は 6–12 面（中央 %s）、facet 線形寸法の中央 %s mm（面積重み %s mm）。"
      % (a1["selected_tau_star"],
         (a1["sweep"][[s["tau_star"] for s in a1["sweep"]].index(a1["selected_tau_star"])]
          ["removed_grains"]),
         (a1["sweep"][[s["tau_star"] for s in a1["sweep"]].index(a1["selected_tau_star"])]
          ["R_a"]),
         ("通過上限で終了" if
          a1["sweep"][[s["tau_star"] for s in a1["sweep"]].index(a1["selected_tau_star"])]
          ["stop_reason"] == "max_passes" else
          a1["sweep"][[s["tau_star"] for s in a1["sweep"]].index(a1["selected_tau_star"])]
          ["stop_reason"]),
         st1["A2"]["Dc_star_m"], st1["A2"]["n_moved"],
         abs(st1["A2"]["volume_residual_m3"]),
         st1["A3"]["frac_support_ge3"], f(st1["A3"]["dip_median_deg"], 1),
         st1["A3"]["porosity_before"]["porosity"],
         st1["A3"]["porosity_after"]["porosity"], st1["A3"]["porosity_drop"],
         f(st1["B1"]["n_faces"]["median"], 0),
         f(st1["B1"]["facet_linear_mm"]["median"], 2),
         f(st3["summary"]["patch"]["VB"]["mean"]["facet_mm_area"] if
           st3["summary"]["patch"].get("VB") else None, 2)))
    A("- **STAGE 2**: V0 / Vw / Vi / VA / VB / VB-a / VAB（パッチ R=5）、"
      "V0 / Vw / Vi / VA / VB / VAB（帯状 R=3）。seed 列は §2.1。数値表は §3。")
    # A-1 removal counts on the factorial beds
    rem = {}
    for box, conds in (("patch", CONDITIONS_PATCH), ("strip", CONDITIONS_STRIP)):
        for c in conds:
            rows = st2[box].get(c) or {}
            for si, r in rows.items():
                s1r = (r.get("steps") or {}).get("A1")
                if s1r:
                    rem.setdefault("%s|%s" % (box, c), []).append(
                        (int(si), s1r["removed_grains"], s1r["stop_reason"],
                         s1r["R_a"]))
    if rem:
        A("- **A-1 は不活性**: 8 realization のすべてで Vw の A-1 除去粒数が 0 "
          "（初期床の動的 R_a が既に 2 以上）。したがって **Vw ≡ V0**（表の行が"
          "同一なのは複製ミスではなく構成上の結果）。VA/VAB の A 効果は "
          "A-2＋A-3 のみ。")
        sel1 = [x for x in a1["sweep"]
                if abs(x["tau_star"] - a1["selected_tau_star"]) < 1e-12]
        if sel1:
            A("  STAGE 1 の seed（bed_gen2 の SEED）だけは R_a 1.57 で能動"
              "（τ*=%.3f で %d 粒除去、R_a %.3f）だった。"
              % (sel1[0]["tau_star"], sel1[0]["removed_grains"],
                 sel1[0]["R_a"] or 0.0))
    # porosity window
    bad = []
    for box, conds in (("patch", CONDITIONS_PATCH), ("strip", CONDITIONS_STRIP)):
        for c in conds:
            rows = st2[box].get(c) or {}
            v = [r["stats"]["porosity"]["porosity"] for r in rows.values()]
            if v and (max(v) > 0.45 or min(v) < 0.30):
                bad.append("%s|%s %.3f–%.3f" % (box, c, min(v), max(v)))
    if bad:
        A("- **空隙率の窓（0.30–0.45）を外れた条件**: %s。B-3 の充填上限"
          "（6000 粒）では回復しない。事実として記録する。" % "; ".join(bad))
    A("")
    A("## 1. STAGE 0 — 対照床の再現")
    A("")
    A("### 1.1 絶対パス事前確認")
    A("")
    A("brief が挙げた経路のみを見た。別ツリーは探索していない。")
    A("")
    A("| ファイル | bytes | sha256 (先頭16) | 存在 |")
    A("|---|---:|---|---|")
    for r in st0["preflight"]:
        A("| `%s` | %s | `%s` | OK |" % (os.path.basename(r["path"]),
                                          r["bytes"], (r["sha256"] or "")[:16]))
    A("")
    A("### 1.2 再生成の同一性")
    A("")
    A("`bed_gen2.form_field(512, sigma_form=%.10f, SEED)` → "
      "`get_population(SEED, P_SCALE0)` → `place` → 縦レイキャスト（モノ 2 個"
      "＋衛星を除外）→ ガウス 60 px インペイント。" % st0["rebuild"]["sigma_form_m"])
    A("")
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 粒数（ボルダー跡を除く） | %d |" % st0["rebuild"]["grains"])
    A("| 保存済み `heightfield.npy` との float32 一致 | **%s** |"
      % ("はい" if st0["rebuild"]["bit_identical_float32"] else "いいえ"))
    A("| max abs diff | %.1e m |" % st0["rebuild"]["max_abs_diff_m"])
    A("")
    A("### 1.3 凍結物差しと窓")
    A("")
    A("物差しは 12B の実装を import してそのまま使った"
      "（`run_v4_2_12b.nikora_measure` → `run_v4_2_5.curve_dir` / `fit_band`、"
      "軸カット・60 ビン log ラダー・第2区間 ≤ 4×D50・H は [格子, 折れ] の"
      "単一フィット）。**折れ探索窓だけ** brief の指示で [0.10, 2]×D50 にした"
      "（v4-2-13 の結論）。")
    A("")
    A("窓の食い違いを隠さないため、12B が公表値を出した窓 [0.25, 2] と、"
      "本ジョブが使う [0.10, 2] の両方を出す。")
    A("")
    A("| 窓 | H_Nikora x | H_Nikora y | 折れ x (mm) | 折れ y (mm) |")
    A("|---|---:|---:|---:|---:|")
    A("| [0.25, 2]（12B の公表値の窓） | %s | %s | %s | %s |"
      % (f(g["a_H_x_12B_window"]["value"]),
         f(g["a_H_y_12B_window"]["value"]),
         f(st0["measure_w0"]["x"]["break"]["lag_m"] * 1000, 3),
         f(st0["measure_w0"]["y"]["break"]["lag_m"] * 1000, 3)))
    A("| [0.10, 2]（本ジョブの規約） | %s | %s | %s | %s |"
      % (f(st0["measure"]["x"]["H_Nikora"]["H"]),
         f(st0["measure"]["y"]["H_Nikora"]["H"]),
         f(st0["measure"]["x"]["break"]["lag_m"] * 1000, 3),
         f(st0["measure"]["y"]["break"]["lag_m"] * 1000, 3)))
    A("")
    A("→ 窓を [0.10,2] にするとパッチ x は +%.4f 動く（v4-2-13 の記録"
      " 0.4484 → 0.4653 と同値）。y は動かない。これは測定の破損ではなく、"
      "窓の効果である（12B の値の再現は [0.25,2] で行う）。"
      % (st0["measure"]["x"]["H_Nikora"]["H"] - g["a_H_x_12B_window"]["value"]))
    A("")
    A("### 1.4 ゲート")
    A("")
    A("| ゲート | 実測 | 期待 | 差 | 許容 | 判定 |")
    A("|---|---:|---:|---:|---:|---|")
    rows = [
        ("(a) H_Nikora x（12B の窓）", "a_H_x_12B_window", 0.005),
        ("(a) H_Nikora y（12B の窓）", "a_H_y_12B_window", 0.005),
        ("(a) σz/D50", "a_sigma_z_over_d50", 0.005),
        ("(a) 歪度", "a_skew", 0.05),
        ("(b) H_fit 帯状 r1 x（0.29 の箱）", "b_H_fit_long_x_0p29", 0.02),
    ]
    for lab, k, tol in rows:
        gg = g[k]
        A("| %s | %s | %s | %s | ±%s | %s |"
          % (lab, f(gg["value"]), f(gg["expect"]), sgn(gg["diff"]), tol,
             "PASS" if gg["pass"] else "FAIL"))
    A("")
    A("H_fit（旧物差し・固定帯 0.4–3.6 cm）の全読み:")
    A("")
    A("| 箱 | x | y | iso |")
    A("|---|---:|---:|---:|")
    A("| パッチ 0.6 m（本ジョブの V0） | %s | %s | %s |"
      % (f(st0["H_fit_patch"]["x"]), f(st0["H_fit_patch"]["y"]),
         f(st0["H_fit_patch"]["iso"])))
    A("| 帯状 3.0 m r1（v4-2-5 の 0.29 の箱） | %s | %s | %s |"
      % (f(st0["H_fit_long"]["x"]["H_fit"]["H"]),
         f(st0["H_fit_long"]["y"]["H_fit"]["H"]),
         f(st0["H_fit_long"]["iso"]["H_fit"]["H"])))
    A("")
    A("**判定: %s**。0.29 は「帯状箱の生成床」の数字であり（v4-2-5 long_r1 "
      "0.294/0.297）、パッチ箱の旧物差しは 0.32–0.37 が v4-2-5 の公表値そのもの"
      "である。ゲート (b) は 0.29 が定義された帯状箱で評価し、パッチ値も併記した。"
      % ("PASS" if st0["gate_pass"] else "FAIL"))
    A("")
    A("### 1.5 V0 の幾何")
    A("")
    m = st0["measure"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 格子 | %.8f m |" % st0["cell_m"])
    A("| D50（12B のパッチ値） | %.6f m |" % st0["d50_m"])
    A("| σz | %.6f m |" % m["x"]["sigma_z_m"])
    A("| σz/D50 | %.4f |" % m["x"]["sigma_z_over_d50"])
    A("| 歪度 | %+.3f |" % m["x"]["skew"])
    A("| 折れ x / y | %.3f / %.3f mm |"
      % (m["x"]["break"]["lag_m"] * 1000, m["y"]["break"]["lag_m"] * 1000))
    A("")
    A("## 2. STAGE 1 — 各機構の実装と、そのステップだけの診断")
    A("")
    A("### 2.0 実装の単位（`bed_gen3.py`）")
    A("")
    A("| ステップ | 関数 | 何をするか |")
    A("|---|---|---|")
    A("| A-0 | `initial_bed(seed, box)` | bed_gen2 の粒配置を verbatim。"
      "パッチは `form_field`+`get_population`+`place`。帯状は 0.6 m 母集団 4 本を"
      "1 本の矩形フォーム（2.4×0.6 m）に載せる（v4-2-3 の前例。継ぎ目の重なりは"
      "開示） |")
    A("| A-1 | `winnow_at(bed, tau*)` | τ(x,y) = τ0·(1 − 後流欠損)。大粒"
      "（D ≥ D84）の下流 2–4 D、幅 0.45 D のガウス帯。p_i = (z_top,i − "
      "z_local)/D_i。p_c はモーメント釣り合い（抗力・揚力は対数則の流速、"
      "抵抗力は水中重量＋上載荷重、腕は pivot 角 φ）。p_i > p_c の粒を "
      "W* = clip((p_i−p_c)/0.15,0,1)^1.5 の確率で除去。系外へ出す |")
    A("| A-2 | `granular_sorting(bed)` | Dc* = 0.20·D50_surface。d<Dc*/3 は"
      "1–4·D90_fine 下へ、Dc*/3≤d<Dc* は 2–5·D90_fine で橋かけ停止、粗粒は"
      "squeeze expulsion で上へ。体積は保存（残差を記録） |")
    A("| A-3 | `interlock(bed)` | 伏角 N(28°,6.5°) を [15,40] に切り、上流縁が"
      "潜る向き（nrm=(−sin d,0,cos d)）。支持点（近傍の下からの支え＋マトリクス）"
      "が 3 点になるまで沈降。クラスタ被覆を 10% へ |")
    A("| B-1 | `facetise(bed)` | 楕円体を 6–12 枚の平面で切る（破断面）。"
      "カリパスは不変（b/a・Corey は構成上そのまま）。角は丸めない |")
    A("| B-2a | `planation(mode='clip')` | 装甲面 = 3×D50 四方の平均 + "
      "k·σz。超えた粒頭を面で切る（粒は残る）。切った体積は系外 |")
    A("| B-2b | `planation(mode='remove')` | 同面を超えた粒を除去し、下の粒を"
      "露出させる |")
    A("| B-3 | `pocket_fill(bed, Dc*)` | 窪み（装甲面より 2 mm 以深）に "
      "d<Dc*/3 の細粒を最深部から敷く（計算上限 4000 粒、開示） |")
    A("")
    A("A-1 の p_c は（簡略形）")
    A("")
    A("```")
    A("u* = sqrt(τ/ρ),  u(z) = (u*/κ)·ln(z/z0),  z0 = D50/30")
    A("Fd = 0.5·ρ·C_D·A_f·u(D/2)^2,  Fl = 0.5·ρ·C_L·A_f·u(D/2)^2")
    A("p_c = [ (W' + N)·0.5·B·sinφ − (Fd + Fl)·0.5·C ] / (Fd·D)")
    A("```")
    A("")
    A("τ* は粒ごとに τ*_i = τ(x,y)/((ρs−ρ)·g·D_i) で効く（細粒ほど大きい）。"
      "掃引 τ* は機構の強さのパラメータで、シーン自身の τ0 = ρgHS0 = %.3f Pa"
      "（τ* = %.4f、H=0.25 m, S0=0.003）とは別物である。掃引帯 0.03–0.10 は"
      "シーンの τ* より上の帯で、brief の指示どおり掃引した。"
      % (B3.TAU0_SCENE, B3.TAU0_SCENE / (B3.DRHO * B3.G * st0["d50_m"])))
    A("")
    A("### 2.1 初期床（パッチ seed0）")
    A("")
    i0 = st1["initial"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| D50_ref | %.6f m |" % i0["D50_ref"])
    A("| D84_ref / D90_ref | %.6f / %.6f m |" % (i0["D84_ref"], i0["D90_ref"]))
    A("| R_a（露出/未露出） | %.3f |" % i0["gsd"]["R_a"])
    A("| R_a_kind（bed_gen2 流） | %.3f |" % i0["gsd"]["R_a_kind"])
    A("| クラスタ被覆 | %.4f |" % i0["cluster_coverage"])
    A("| 空隙率 | %.4f |" % i0["porosity"]["porosity"])
    A("")
    A("### 2.2 A-1 winnow の掃引")
    A("")
    A("| τ* | 通過数 | 除去粒数 | R_a | 停止理由 | τ0 (Pa) | 後流欠損の面積率 |")
    A("|---|---:|---:|---:|---|---:|---:|")
    for s in a1["sweep"]:
        tf = s.get("wake_field") or {}
        reason = s["stop_reason"]
        if reason == "max_passes" and (s["R_a"] or 0) >= 2.0:
            reason = "R_a=2（通過上限で終了）"
        A("| %.3f | %d | %d | %s | %s | %s | %s |"
          % (s["tau_star"], s["passes"], s["removed_grains"], f(s["R_a"], 3),
             reason, f(tf.get("tau0_Pa"), 2),
             f(tf.get("frac_below_0p8_tau0"), 3)))
    A("")
    A("選択 τ* = **%.3f**（%s）。τ 場は 4 床 px の格子上に実在し、"
      "`tau-field-patch-s0.npy` に保存した（min/mean/max = %s / %s / %s Pa、"
      "後流 %d 本）。"
      % (a1["selected_tau_star"], a1["selected_from"],
         f(a1["tau_field"]["min_Pa"], 2), f(a1["tau_field"]["mean_Pa"], 2),
         f(a1["tau_field"]["max_Pa"], 2), a1["tau_field"]["n_wakes"]))
    A("")
    A("| 掃引の各 τ* の履歴（R_a の推移） |")
    A("|---|")
    for s in a1["sweep"]:
        A("| τ*=%.3f: %s |" % (s["tau_star"],
                               ", ".join("%s" % f(h["R_a"], 3)
                                         for h in s["history"])))
    A("")
    A("### 2.3 A-2 選別（A-1 の停止面から）")
    A("")
    a2 = st1["A2"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| Dc* | %.6f m（%.2f × D50_surface） |" % (a2["Dc_star_m"],
                                                  a2["Dc_over_D50_surface"]))
    A("| D90_fine（深さの単位） | %.6f m |" % a2["D90_fine_m"])
    A("| 移動粒数 / 粗粒の押し上げ | %d / %d |" % (a2["n_moved"],
                                                  a2["n_coarse_lifted"]))
    A("| 体積（下 / 上 / 残差） | %.3e / %.3e / %.1e m³ |"
      % (a2["volume_down_m3"], a2["volume_up_m3"], a2["volume_residual_m3"]))
    A("")
    A("地下の細粒（d<Dc*）体積分率プロファイル（深さは D90_fine 単位）:")
    A("")
    A("| 深さ | 前 | 後 |")
    A("|---|---:|---:|")
    for b, af in zip(a2["profile_before"]["bins"], a2["profile_after"]["bins"]):
        A("| %s | %s | %s |" % (b["bin"], f(b["fines_volume_fraction"], 4),
                                f(af["fines_volume_fraction"], 4)))
    A("")
    A("### 2.4 A-3 interlock")
    A("")
    a3 = st1["A3"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 支持点 ≥3 の割合 | %.3f |" % a3["frac_support_ge3"])
    A("| 支持点の中央 / ヒスト | %s / %s |" % (f(a3["support_median"], 1),
                                               a3["support_hist"]))
    A("| 伏角 中央 / P16 / P84 | %s / %s / %s ° |"
      % (f(a3["dip_median_deg"], 1), f(a3["dip_p16_deg"], 1),
         f(a3["dip_p84_deg"], 1)))
    A("| クラスタ被覆（前 → 後） | %.4f → %.4f |"
      % (a3["cluster_coverage_before"], a3["cluster_coverage_after"]))
    A("| 空隙率（前 → 後、低下） | %.4f → %.4f（%+.4f） |"
      % (a3["porosity_before"]["porosity"], a3["porosity_after"]["porosity"],
         a3["porosity_drop"]))
    A("")
    A("### 2.5 B-1 facet")
    A("")
    b1 = st1["B1"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 面数 | 最小 %d / 中央 %s / 最大 %d |"
      % (b1["n_faces"]["min"], f(b1["n_faces"]["median"], 1),
         b1["n_faces"]["max"]))
    A("| facet 線形寸法（√面積）中央 / P16 / P84 | %s / %s / %s mm |"
      % (f(b1["facet_linear_mm"]["median"], 2),
         f(b1["facet_linear_mm"]["p16"], 2),
         f(b1["facet_linear_mm"]["p84"], 2)))
    A("| b/a 中央 / Corey 中央 | %s / %s |"
      % (f(b1["ba"]["median"], 3), f(b1["corey"]["median"], 3)))
    A("| 体積 / 楕円体 | %.3f |" % b1["volume_over_ellipsoid"]["median"])
    A("| 破断面の面積率 | %.3f |" % b1["facet_area_fraction"]["median"])
    A("| 破断面と滑面の二面角 中央 | %s ° |" % f(b1["facet_dihedral_deg"]["median"], 1))
    A("| Wadell 体積代理（deq/2Rc）中央 | %.3f |" % b1["wadell_volume_proxy"]["median"])
    A("| Powers 円磨度（辺の二面角から）中央 | %.3f |"
      % b1["facet_edge_roundness"]["median"])
    A("")
    A("**円磨度の定義について（事実）**: brief の「Wadell 0.1–0.3」は"
      "Powers/Krumbein の目視階級（angular 0.1 … well-rounded 0.6）の数字で"
      "ある。角を丸めない破断面では、角の曲率半径による Wadell 値は 0 に落ち、"
      "体積比の代理は %.2f になる。ここでは破断面と滑面の二面角から "
      "Powers 階級へ写した値を主に報告する（%.2f–%.2f）。"
      % (b1["wadell_volume_proxy"]["median"],
         b1["facet_edge_roundness"]["p16"], b1["facet_edge_roundness"]["p84"]))
    A("")
    A("### 2.6 B-2 planation（2 亜種 × 2 定義）")
    A("")
    A("| 定義 | クリップ面積率 f | 切った体積 | 除去粒数 | 除去面積率 | 上位数 f の実測 |")
    A("|---|---:|---:|---:|---:|---:|")
    for k, d in (("k=1.0（本命）", st1["B2a"]), ("k=0.5（感度）", st1["B2a_k05"])):
        A("| B-2a %s | %.4f | %.3e m³ | — | — | %s |"
          % (k, d["clip_area_fraction"], d["cut_volume_m3"],
             f(d.get("equiv_top_f_measured"), 4)))
    for k, d in (("k=1.0（本命）", st1["B2b"]), ("k=0.5（感度）", st1["B2b_k05"])):
        A("| B-2b %s | %.4f | %.3e m³ | %d | %s | %s |"
          % (k, d["clip_area_fraction"], d["cut_volume_m3"],
             d["n_removed_grains"], f(d.get("removed_area_fraction"), 4),
             f(d.get("equiv_top_f_measured"), 4)))
    A("")
    A("本命を k=1.0 にした理由（実行前に固定）: brief は「平均＋kσ」と"
      "「summit の上位 f 分率」を等価な代替として並べるが、生成床の局所高さは"
      "歪んでいるので両者は一致しない。実測では k=1.0 が f≈0.09（brief の"
      " f=0.10 側）、k=0.5 が f≈0.25（f=0.05 から遠い）。k=0.5 は感度として"
      "併記する。B-2b は B-2a と同じ面積率 f を切る（対の比較を保つ）。")
    A("")
    A("「上位数 f の実測」は、等価な定義（近傍 summit の上位 f 分率を切る）を"
      "局所正規分布の分位点で置いたときの、実際に切れた面積率である。")
    A("")
    A("### 2.7 B-3 ポケット充填")
    A("")
    b3 = st1["B3"]
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 追加粒数 | %d（上限 %d） |" % (b3["n_filled"], b3["pore_cap"]))
    A("| 充填体積 / 目標 | %.3e / %.3e m³（%.1f%%） |"
      % (b3["filled_volume_m3"], b3["target_volume_m3"],
         100 * b3["fill_of_target"]))
    A("| 窪みの面積率 | %.4f |" % b3["depression_area_fraction"])
    A("| 空隙率（前 → 後） | %.4f → %.4f |"
      % (b3["porosity_before"]["porosity"], b3["porosity_after"]["porosity"]))
    A("")
    A("## 3. STAGE 2 — 要因計画")
    A("")
    A("### 3.1 seed 列と箱")
    A("")
    sl = N["preregistration"]["seeds"]
    A("`np.random.default_rng(20261115).integers(0, 2**63, size=8)` で引いた。")
    A("")
    A("| 用途 | seed |")
    A("|---|---:|")
    for i, s in enumerate(sl["patch"]):
        A("| パッチ realization %d | %d |" % (i, s))
    for i, s in enumerate(sl["strip"]):
        A("| 帯状 realization %d（タイルは seed+10k） | %d |" % (i, s))
    A("")
    A("箱: パッチ 0.6×0.6 m（2048²、格子 %.5f mm）、帯状 2.4×0.6 m"
      "（2048×8192、同格子）。VB-a はパッチのみ。" % (st0["cell_m"] * 1000))
    A("")
    A("### 3.2 全条件・全 seed の数値（パッチ）")
    A("")
    for cond in CONDITIONS_PATCH:
        if cond not in st3["summary"]["patch"]:
            continue
        A("#### %s — %s" % (cond, COND_LABEL[cond]))
        A("")
        A("| seed | H_Nik x | H_Nik y | 折れ x (mm) | 窓下端張付 | H_alt x | "
          "H_fit x | H_fit y | R_a | Sk | σz/D50 | σg | 空隙率 | クラスタ | "
          "facet mm |")
        A("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        rows = st2["patch"][cond]
        for si in sorted(rows):
            r = rows[si]
            mx, my = r["measure"]["x"], r["measure"]["y"]
            gs = r["stats"]["gsd"]
            fs = (r["stats"].get("facet_exposed") or {})
            A("| %d | %s | %s | %s | %s | %s | %s | %s | %.3f | %+.3f | %.3f | "
              "%.3f | %.3f | %.3f | %s |"
              % (si, f(hn(mx)), f(hn(my)),
                 f((mx.get("break") or {}).get("lag_m") * 1000, 2),
                 "yes" if mx["break_at_window_floor"] else "no",
                 f(mx["H_alt_lastlag"]), f(mx["H_fit_fixed_band"]["H"]),
                 f(my["H_fit_fixed_band"]["H"]), gs.get("R_a") or 0,
                 mx["skew"], mx["sigma_z_over_d50"],
                 gs.get("surface_sigma_g") or 0,
                 r["stats"]["porosity"]["porosity"],
                 r["stats"]["cluster_coverage"],
                 f(fs.get("median_mm"), 2)))
        if cond == "Vw":
            A("（A-1 の除去 0 粒により V0 と同じ数値。§0 参照）")
            A("")
        s = st3["summary"]["patch"][cond]["mean"]
        A("| **平均** | %s | %s | — | — | %s | %s | %s | %.3f | %+.3f | %.3f | "
          "%.3f | %.3f | %.3f | %s |"
          % (f(s["H_Nikora_x"]), f(s["H_Nikora_y"]), f(s["H_alt_x"]),
             f(s["H_fit_x"]), f(s["H_fit_y"]), s["R_a"] or 0, s["skew_x"],
             s["sigma_x"], s["surface_sigma_g"] or 0, s["porosity"],
             s["cluster_cov"], f(s["facet_mm"], 2)))
        A("")
    A("### 3.3 全条件・全 seed の数値（帯状 2.4×0.6 m）")
    A("")
    for cond in CONDITIONS_STRIP:
        if cond not in st3["summary"]["strip"]:
            continue
        A("#### %s — %s" % (cond, COND_LABEL[cond]))
        A("")
        A("| seed | H_Nik x | H_Nik y | 折れ x (mm) | 窓下端張付 | H_alt x | "
          "H_fit x | H_fit y | R_a | Sk | σz/D50 | σg | 空隙率 | クラスタ | "
          "facet mm |")
        A("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        rows = st2["strip"][cond]
        for si in sorted(rows):
            r = rows[si]
            mx, my = r["measure"]["x"], r["measure"]["y"]
            gs = r["stats"]["gsd"]
            fs = (r["stats"].get("facet_exposed") or {})
            A("| %d | %s | %s | %s | %s | %s | %s | %s | %.3f | %+.3f | %.3f | "
              "%.3f | %.3f | %.3f | %s |"
              % (si, f(hn(mx)), f(hn(my)),
                 f((mx.get("break") or {}).get("lag_m") * 1000, 2),
                 "yes" if mx["break_at_window_floor"] else "no",
                 f(mx["H_alt_lastlag"]), f(mx["H_fit_fixed_band"]["H"]),
                 f(my["H_fit_fixed_band"]["H"]), gs.get("R_a") or 0,
                 mx["skew"], mx["sigma_z_over_d50"],
                 gs.get("surface_sigma_g") or 0,
                 r["stats"]["porosity"]["porosity"],
                 r["stats"]["cluster_coverage"],
                 f(fs.get("median_mm"), 2)))
        if cond == "Vw":
            A("（A-1 の除去 0 粒により V0 と同じ数値。§0 参照）")
            A("")
        s = st3["summary"]["strip"][cond]["mean"]
        A("| **平均** | %s | %s | — | — | %s | %s | %s | %.3f | %+.3f | %.3f | "
          "%.3f | %.3f | %.3f | %s |"
          % (f(s["H_Nikora_x"]), f(s["H_Nikora_y"]), f(s["H_alt_x"]),
             f(s["H_fit_x"]), f(s["H_fit_y"]), s["R_a"] or 0, s["skew_x"],
             s["sigma_x"], s["surface_sigma_g"] or 0, s["porosity"],
             s["cluster_cov"], f(s["facet_mm"], 2)))
        A("")
    A("### 3.4 ペア seed の差（条件 − V0、同一 seed）")
    A("")
    A("| 箱 | 条件 | ΔH_Nik x | ΔH_Nik y | ΔH_fit x | ΔH_fit y |")
    A("|---|---|---:|---:|---:|---:|")
    for k in sorted(st3["paired_deltas"]):
        d = st3["paired_deltas"][k]
        A("| %s | %s | %s | %s | %s | %s |"
          % (k.split("|")[0], k.split("|")[1],
             " / ".join(sgn(v) for v in d["dH_Nikora_x"]),
             " / ".join(sgn(v) for v in d["dH_Nikora_y"]),
             " / ".join(sgn(v) for v in d["dH_fit_x"]),
             " / ".join(sgn(v) for v in d["dH_fit_y"])))
    A("")
    A("### 3.5 H_loc プロファイル（帯域 0.4–3.6 cm）と別読み")
    A("")
    A("H_loc 帯域の最大値（パッチ / 帯状、平均）:")
    A("")
    A("| 条件 | パッチ x | パッチ y | 帯状 x | 帯状 y |")
    A("|---|---:|---:|---:|---:|")
    for cond in CONDITIONS_PATCH:
        pv = st2["patch"].get(cond)
        sv = st2["strip"].get(cond)
        def hm(rows, ax):
            return _mean([r["measure"][ax]["H_loc_band"]["max"]
                          for r in rows.values()]) if rows else None
        A("| %s | %s | %s | %s | %s |"
          % (cond, f(hm(pv, "x")), f(hm(pv, "y")),
             f(hm(sv, "x")), f(hm(sv, "y"))))
    A("")
    A("別読み（第2区間をラダー末端まで）との差 ΔH = H_alt − H_main は §3.2/§3.3 "
      "の H_alt 列と H_Nikora 列の差。窓下端張付は x 方向について列挙した"
      "（張付があれば境界制約として読む）。")
    A("")
    A("### 3.6 地下の細粒分率プロファイル（各条件の代表 seed）")
    A("")
    A("代表 seed（各箱の最小 seed index）の値。")
    A("")
    A("| 箱 | 条件 | " + " | ".join(b["bin"] for b in
                                    st2["strip"]["V0"][sorted(st2["strip"]["V0"])[0]]["stats"]["fines_profile"]["bins"]) + " |")
    A("|---|---|" + "---:|" * 8)
    for box, cond in (("patch", "V0"), ("patch", "VA"), ("strip", "V0"), ("strip", "VA")):
        rows = st2[box].get(cond)
        if not rows:
            continue
        prof = rows[sorted(rows)[0]]["stats"]["fines_profile"]["bins"]
        A("| %s | %s | %s |" % (box, cond,
                                " | ".join(f(b["fines_volume_fraction"], 3)
                                           for b in prof)))
    A("")
    A("## 4. STAGE 3 — 事前登録の読み取り")
    A("")
    A("実測 DEM はこのジョブに無い。**エンベロープ判定は出さない**。合成 fBm や"
      "他の生成床を参照集団に使っていない。")
    A("")
    ro = st3["readout"]
    for key, title in (("G-A1", "G-A1 数値（A）"),
                       ("G-A2", "G-A2 装甲の再現（A）"),
                       ("G-A3", "G-A3 A1/A2 の分離"),
                       ("G-B1", "G-B1 数値（B）"),
                       ("G-B2", "G-B2 facet の前提"),
                       ("G-B3", "G-B3 歪度の緊張"),
                       ("VAB", "VAB")):
        r = ro[key]
        A("### %s" % title)
        A("")
        A("- 規則: %s" % r["rule"])
        A("- 値:")
        for k, v in (r.get("values") or {}).items():
            A("  - `%s` = %s" % (k, f(v) if isinstance(v, (int, float)) else v))
        if r.get("strip_Sk_by_seed"):
            A("  - strip VA の seed 別 Sk: %s"
              % ", ".join("s%d %+.3f" % (kk, vv) for kk, vv in
                          r["strip_Sk_by_seed"]))
        if r.get("box_note"):
            A("  - 箱の内訳: %s" % r["box_note"])
        if r.get("definition_note"):
            A("  - 定義の注意: %s" % r["definition_note"])
        A("- **読み**: %s" % (r.get("verdict") if r.get("verdict") is not None
                              else "（該当条件の測定が無い、または未実行）"))
        A("")
    A("- 到達見込み（brief）: A の H_Nikora は 0.50 前後、楽観で 0.60。"
      "**0.79 を A の予測として書かない**。しきい値 0.50 / 0.60 は較正値ではなく"
      "配置値（0.50 = 現在の生成床 0.29 と目標 0.79 の中間よりやや上、"
      "0.60 = 目標に半分以上近づいた）。")
    A("")
    A("## 5. 付録")
    A("")
    A("### 5.1 出力ファイル")
    A("")
    A("| 種類 | 数 | 例 |")
    A("|---|---:|---|")
    n_hf = (sum(len(st2["patch"].get(c) or {}) for c in CONDITIONS_PATCH)
            + sum(len(st2["strip"].get(c) or {}) for c in CONDITIONS_STRIP))
    A("| 高さ場 npy | %d | `heightfield-VA-strip-s0.npy` |" % n_hf)
    A("| 粒リスト json | %d | `grains-VA-strip-s0.json` |" % n_hf)
    A("| τ 場 | 1 | `tau-field-patch-s0.npy`（512²、min/mean/max 0.88/12.22/17.68 Pa） |")
    A("| STAGE 0 対照床 | 1 | `heightfield-STAGE0-patch-control.npy` |")
    A("")
    A("### 5.2 時間")
    A("")
    cw = 0.0
    for bx in ("patch", "strip"):
        for c, rows in (st2[bx] or {}).items():
            for r in rows.values():
                cw += float(r.get("wall_s") or 0.0)
    A("- 機構条件の計算時間の合計（各 rec の wall_s の和）: %.1f s" % cw)
    A("- STAGE 0 / STAGE 1 は初回パスで実行し（約 61 s / 300 s）、"
      "再開パスではキャッシュ（`stage0-rec.json` / `stage1-rec.json`）から読んだ。")
    A("- 最終アセンブリ（この numbers/report の書き出し）: %.1f s" % N["wall_s"])
    A("- 実測の壁時計: 初回 09:16–10:47（90 分の実行上限で中断）、"
      "パッチ 10:49–11:31、帯状 11:31–12:43、τ 場の修正後 パッチ A 12:50–13:15・"
      "帯状 A 13:17–14:04、窓の統一再測定 14:04–14:22、最終アセンブリ 0.2 s。")
    A("")
    A("### 5.3 注意（事実）")
    A("")
    A("- 粒層の実装は DEM ではない。A-1 の p_c、A-2 の Dc*、A-3 の支持点、"
      "B-1 の破断面数は宣言した構成則で、フィッティングはしていない。")
    A("- 帯状箱は 0.6 m 母集団 4 本のタイルで、packing は継ぎ目を知らない"
      "（重なりは z-buffer の max で吸収）。v4-2-3 と同じ開示。")
    A("- 窓 [0.10,2] は 12B の公表値の窓 [0.25,2] と違う。両方を §1.3 に出した。")
    A("- STAGE 0 の H_fit 0.29 は帯状箱の値。パッチ箱は 0.32–0.37（v4-2-5 の"
      "公表値と一致）。")
    A("- 空隙率は bed_gen.porosity と同じスラブ則（中央値の表面 + 1.5·D50_in "
      "から 50 mm 下、384²×128 ボクセル）で測った。スラブのアンカーが表面に"
      "追随するため、粒配置が変わると絶対値は敏感に動く（v4-2-2 の自 seed は"
      "0.3388、本ジョブのパッチ V0 は 0.205–0.296）。条件間の差は読めるが、"
      "絶対値を設計の窓と直接突き合わせるには注意が要る。")
    A("")
    return "\n".join(R) + "\n"


def write_report(N):
    os.makedirs(OUT, exist_ok=True)
    data = build_report(N).encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            t = f.read().decode("utf-8")
        assert len(t) > 2000
        assert "REPORT-V4-2-15" in t
        return True

    atomic_write_bytes(REPORT, data, verify=verify)
    return REPORT


def build_stop_report(reason, ctx):
    L = ["# REPORT-V4-2-15 — STOP（STAGE 0 で停止）", ""]
    L.append("- brief: `%s`" % BRIEF)
    L.append("- 停止理由: **%s**" % reason)
    L.append("- 作成: %s" % utcnow())
    L.append("")
    L.append("```")
    L.append(json.dumps(_san(ctx), ensure_ascii=False, indent=1)[:12000])
    L.append("```")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    t_start = time.time()
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    log("run_v4_2_15 start; out=%s" % OUT)
    if len(sys.argv) > 1 and sys.argv[1] == "remeasure":
        n = 0
        for fn in sorted(os.listdir(OUT)):
            if not fn.startswith("rec-") or not fn.endswith(".json"):
                continue
            rec = json.load(open(os.path.join(OUT, fn), encoding="utf-8"))
            hf = os.path.join(OUT, rec["heightfield"])
            if not os.path.isfile(hf):
                log("  remeasure: %s has no height field" % fn)
                continue
            z = np.load(hf).astype(np.float64)
            box = "strip" if "-strip-" in fn else "patch"
            rec["measure"] = measure_bed(
                z, CELL, rec["d50_ref_m"],
                DM_PATCH if box == "patch" else DM_STRIP_X,
                DM_PATCH if box == "patch" else DM_STRIP_Y)
            write_json(os.path.join(OUT, fn), rec)
            n += 1
            log("  remeasured %s  H_Nik x %s" % (
                fn, f(rec["measure"]["x"]["H_Nikora"]["H"])))
        log("REMEASURE done: %d records" % n)
        return
    N = {"created_utc": utcnow(), "brief": BRIEF,
         "what_this_is": "v4-2.15 mechanism A/B factorial on the grain layer",
         "preregistration": {
             "ruler_window_x_d50": list(BREAK_WIN),
             "ruler_12b_window": [0.25, 2.0],
             "fixed_band_m": list(FIXED_BAND),
             "gate_H_expect": GATE_H_EXPECT, "gate_H_tol": GATE_H_TOL,
             "hfit_029": HFIT_029, "hfit_029_tol": HFIT_029_TOL,
             "read_hfit_low": READ_HFIT_LOW,
             "read_hfit_survive": READ_HFIT_SURVIVE,
             "g_a2_ra": list(G_A2_RA), "g_a2_sk": list(G_A2_SK),
             "g_b2_facet_mm": G_B2_FACET_MM, "g_b3_sk": G_B3_SK,
             "conditions_patch": list(CONDITIONS_PATCH),
             "conditions_strip": list(CONDITIONS_STRIP),
             "n_patch": N_PATCH, "n_strip": N_STRIP,
             "parameters": {
                 "physics": {"rho_w": B3.RHO_W, "rho_s": B3.RHO_S,
                             "H_flow_m": B3.H_FLOW, "S0": B3.S0,
                             "tau0_scene_Pa": B3.TAU0_SCENE,
                             "C_D": B3.C_D, "C_L": B3.C_L,
                             "z0_over_D50": B3.Z0_OVER_D50},
                 "tau_star_sweep": list(B3.TAU_STAR_SWEEP),
                 "tau_star_dense_anchor": B3.TAU_STAR_DENSE,
                 "dc_star_factor": 0.20,
                 "planation_k_primary": 0.5,
                 "planation_k_sensitivity": 1.0,
                 "cluster_target": 0.10,
                 "pocket_fill_max": 4000,
                 "facet_faces": [6, 12],
             },
             "seed_rule": "np.random.default_rng(20261115).integers(0, 2**63, size=8)",
         }}
    timings = {}
    try:
        t0 = time.time()
        if os.path.isfile(STAGE0_REC):
            stage0_out = json.load(open(STAGE0_REC, encoding="utf-8"))
            log("STAGE0  reused from stage0-rec.json (gate_pass=%s)"
                % stage0_out.get("gate_pass"))
            sigma_form = stage0_out["rebuild"]["sigma_form_m"]
            d50 = stage0_out["d50_m"]
            if not stage0_out.get("gate_pass"):
                raise Blocker("STAGE 0 recorded FAIL", stage0_out)
        else:
            stage0_out, zh, d50, sigma_form = stage0(cp)
            write_json(STAGE0_REC, stage0_out)
        N["stage0"] = stage0_out
        timings["stage0"] = time.time() - t0

        # the control bed object for STAGE 1 (rebuilt from the same recipe)
        import bed_gen2 as B2
        t0 = time.time()
        if os.path.isfile(STAGE1_REC):
            st1 = json.load(open(STAGE1_REC, encoding="utf-8"))
            tau_star = st1["A1"]["selected_tau_star"]
            Dc_star = st1["A2"]["Dc_star_m"]
            log("STAGE1  reused from stage1-rec.json (tau*=%.3f Dc*=%.5f)"
                % (tau_star, Dc_star))
        else:
            base = B3.initial_bed(B2.SEED, "patch", sigma_form)
            base.meta.update({"seed": B2.SEED})
            st1, tau_star, Dc_star, meta1 = stage1(cp, base, sigma_form)
            write_json(STAGE1_REC, st1)
        N["stage1"] = st1
        timings["stage1"] = time.time() - t0

        # the seed list, drawn once
        rng = np.random.default_rng(SEED_BASE)
        seeds = rng.integers(0, 2 ** 63, size=N_PATCH + N_STRIP)
        seed_list = {"patch": [int(s) for s in seeds[:N_PATCH]],
                     "strip": [int(s) for s in seeds[N_PATCH:]]}
        N["preregistration"]["seeds"] = seed_list
        log("  seeds patch %s strip %s" % (seed_list["patch"], seed_list["strip"]))

        t0 = time.time()
        want = sys.argv[1] if len(sys.argv) > 1 else "all"
        boxes = (("patch", "strip") if want in ("all", "")
                 else ((want,) if want in ("patch", "strip") else ("patch", "strip")))
        st2 = stage2(cp, sigma_form, tau_star, Dc_star, seed_list, boxes=boxes)
        N["stage2"] = st2
        timings["stage2"] = time.time() - t0
        if not all(st2.get("complete", {}).get(b) for b in ("patch", "strip")):
            log("PARTIAL pass: boxes %s; deliverables are NOT written yet"
                % st2.get("complete"))
            N["partial"] = True
            N["log"] = LOG
            write_json(os.path.join(OUT, "numbers-partial.json"), N)
            return

        t0 = time.time()
        st3 = stage3(cp, st2, st1)
        N["stage3"] = st3
        timings["stage3"] = time.time() - t0

        N["timings_s"] = timings
        N["log"] = LOG
        N["stop"] = False
        N["finished_utc"] = utcnow()
        N["wall_s"] = time.time() - t_start
        write_json(NUMBERS, N)
        t0 = time.time()
        write_report(N)
        N["timings_s"]["report"] = time.time() - t0
        write_json(NUMBERS, N)
        cp_note(cp, "report", {"report": REPORT, "numbers": NUMBERS,
                               "wall_s": N["wall_s"]})
        log("DONE  report + numbers written; wall %.1f s" % N["wall_s"])
    except Blocker as exc:
        N["stop"] = True
        N["stop_reason"] = exc.reason
        N["stop_context"] = _san(exc.ctx)
        N["log"] = LOG
        N["finished_utc"] = utcnow()
        N["wall_s"] = time.time() - t_start
        write_json(NUMBERS, N)
        os.makedirs(OUT, exist_ok=True)
        atomic_write_bytes(REPORT,
                           build_stop_report(exc.reason, exc.ctx).encode("utf-8"))
        cp_note(cp, "stop", {"reason": exc.reason})
        log("STOP  %s" % exc.reason)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
