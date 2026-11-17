"""run_v4_2_13 -- BRIEF-V4-2-13: break-window sensitivity of the Nikora ruler,
and a bias calibration of that ruler on synthetic fields with a known H.

NUMBERS ONLY.  No verdict, no design opinion, no parameter change, no render.

What is measured
----------------
STAGE 0  absolute-path preflight (the brief's list only; one missing file is a
         STOP -- no other tree is searched), then the 12B reproduction gate:
         with the frozen window [0.25, 2] x D50 the 12B H_Nikora must come back
         (patch x 0.4484 / y 0.4320, r1 x 0.4350 / y 0.4379, r2 x 0.4358 /
         y 0.4401, +-0.002).  The r2 face is the one 12B regenerated from seed.
STAGE 1  break-window sensitivity on the three generated faces, x and y:
         W0 = [0.25, 2] x D50 (the 12B reference)
         W1 = [0.15, 2] x D50
         W2 = [0.10, 2] x D50
         W3 = [0.05, 2] x D50
         only the lower edge moves; everything else is 12B's implementation
         verbatim.  Per (face, direction, window): H_Nikora, break lag in mm and
         x D50, whether the break stuck to the window's lower edge, the second
         segment's slope.  Plus 12B's alternate reading (segment 2 running to
         the last ladder lag) re-reported at W0.
STAGE 2  known-H synthetic calibration: fBm fields in the patch bed's own box
         (0.6 m, cell 0.293 mm, D50 = 0.0362 m for the window conversion),
         H_gen in {0.50, 0.65, 0.79}, two realizations each, sigma_z matched to
         the generated beds.  Each field is read at W0-W3; bias = estimate -
         H_gen is tabulated per window, and lower-edge sticking is recorded.

Conventions fixed BEFORE the first measurement
----------------------------------------------
* ruler      : run_v4_2_5's curve_dir / curve_iso / fit_band / local_H, imported
               verbatim (the same source 12B imported).
* cell       : 0.6 m / 2048 px = 0.29296875 mm.
* ladder     : curve_dir's own 60-bin geometric ladder, rounded to whole pixels
               and np.unique'd -- identical to 12B.
* break      : two-segment log-log least squares; candidate = ladder point in
               [lo, 2 D50]; segment 1 = [first ladder lag, candidate]; segment 2
               = [candidate, 4 D50]; >= 4 points per segment; minimum total
               residual wins.  Only `lo` is varied.
* H_Nikora   : single fit [first ladder lag (= cell), break], 0.5*dlnS/dlnd.
* sticking   : the chosen break equals the smallest candidate ladder lag that
               the window admits (i.e. the optimiser picked the lowest value it
               was allowed to pick).  Reported as a boolean plus the distance of
               the break above the window edge in x D50.
* sigma_z    : the synthetic fields use the mean of the three generated faces'
               sigma_z; the slope is scale-free so this only fixes the absolute
               S level.
* synthetic  : isotropic periodic spectral synthesis, amplitude |k|^-(H+1)
               (2-D PSD ~ |k|^-(2H+2) -> D(r) ~ r^(2H)), float32 container to
               match the stored beds, own seed per realization.

Stages are checkpointed one jobkit unit at a time (stage0 / stage1 / stage2 /
report).  Foreground; every deliverable is written before exit.
"""
from __future__ import annotations

import datetime
import gc
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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-13.md"

OUT = os.path.join(HERE, "out", "v4-2-13")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-13.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-13.md")
PREFLIGHT = os.path.join(OUT, "stage0-preflight.json")
THIS_FILE = "run_v4_2_13.py"

# --- STAGE 0: the brief's faces, resolved to the files 12B actually read -----
# The brief's parenthetical puts the patch bed under out/v4-2-3/.  That folder
# holds the long bed (heightfield-long-r1.npy); the patch bed 12B measured is
# out/v4-2-2/heightfield.npy (run_v4_2_12b.py: HF_PATCH).  "12B が使ったもの"
# is the authority, so the patch face is that file.  Recorded, not guessed.
HF_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
HF_R1 = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")
HF_R2 = os.path.join(HERE, "out", "v4-2-12b", "heightfield-long-r2-regen.npy")
META_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
META_R1 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta.json")
META_R2 = os.path.join(HERE, "out", "v4-2-3", "longbed-meta-2.json")
R5_PATH = os.path.join(HERE, "run_v4_2_5.py")
R12B_PATH = os.path.join(HERE, "run_v4_2_12b.py")

REQUIRED = [
    ("12B measurement implementation", R12B_PATH),
    ("ruler (verbatim functions)", R5_PATH),
    ("patch bed (the one 12B measured)", HF_PATCH),
    ("long bed r1", HF_R1),
    ("long bed r2 (12B seed regeneration)", HF_R2),
    ("patch meta", META_PATCH),
    ("long bed r1 meta", META_R1),
    ("long bed r2 meta", META_R2),
]

# --- the measurement contract ----------------------------------------------
CELL = 0.6 / 2048.0                   # 0.00029296875 m
LONG_NY, LONG_NX = 2048, 10240
NBIN = 60
DM_PATCH = 0.30                       # half of the 0.6 m box
DM_LONG_X = 1.50                      # half of the 3.0 m bed
DM_LONG_Y = 0.30                      # half of the 0.6 m width
HI_MULT = 2.0                         # fixed upper edge of the break window
TWO_TOP = 4.0                         # fixed top of the second segment
MIN_SEG_POINTS = 4

WINDOWS = [("W0", 0.25), ("W1", 0.15), ("W2", 0.10), ("W3", 0.05)]

# STAGE 0 gate: the 12B published H_Nikora, +-0.002 (brief)
GATE12B = {
    "patch | x": 0.4484, "patch | y": 0.4320,
    "long_r1 | x": 0.4350, "long_r1 | y": 0.4379,
    "long_r2 | x": 0.4358, "long_r2 | y": 0.4401,
}
GATE_TOL = 0.002

# STAGE 2 contract, fixed before generating anything
SYNTH_H = (0.50, 0.65, 0.79)
SYNTH_REAL = 2
SYNTH_SEEDS = {  # explicit, recorded
    (0.50, 0): 1300500, (0.50, 1): 1300501,
    (0.65, 0): 1300650, (0.65, 1): 1300651,
    (0.79, 0): 1300790, (0.79, 1): 1300791,
}
SYNTH_D50_M = 0.0362                  # the brief's conversion value
BAND_REF = (4 * CELL, 32 * CELL)      # generator check band, 1.2-9.4 mm (4-32 px)

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
    def __init__(self, reason, ctx):
        super().__init__(reason)
        self.reason = reason
        self.ctx = ctx


def ruler():
    import run_v4_2_5 as R5
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band"):
        assert hasattr(R5, name), "run_v4_2_5 lost %s" % name
    return R5


def axis_cut(z, cell, axis, dmax, R5, nbin=NBIN):
    return R5.curve_dir(z, cell, dmax, nbin=nbin, axis=axis)


def load_long(path):
    return np.fromfile(path, dtype="<f4").reshape(LONG_NY, LONG_NX).astype(np.float64)


# ---------------------------------------------------------------------------
# the Nikora step -- 12B's break_search with one parameter opened up
# ---------------------------------------------------------------------------
def break_search(d, s, d50, lo_mult):
    """Brute-force two-segment log-log split over ladder lags in
    [lo_mult*D50, HI_MULT*D50]; segment 1 = [first lag, candidate], segment 2 =
    [candidate, TWO_TOP*D50]; both need >= 4 points; minimum residual wins.

    Identical to run_v4_2_12b.break_search except that the lower edge is an
    argument instead of the literal BREAK_LO = 0.25.
    """
    ok = np.isfinite(s) & (s > 0) & (d > 0)
    dd, ss = d[ok], s[ok]
    x, y = np.log(dd), np.log(ss)
    j_top = int(np.searchsorted(dd, TWO_TOP * d50, side="right") - 1)
    lo, hi = lo_mult * d50, HI_MULT * d50
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
            "lower_edge_over_d50": float(lo_mult),
            "upper_segment_top_m": float(dd[j_top]),
            "n_candidates": len(trials), "candidates": trials}


def alt_break_lastlag(d, s, d50, lo_mult):
    """12B's post-hoc alternate: same two-segment fit, but segment 2 runs to the
    last ladder lag instead of 4 x D50.  Robustness reading only."""
    ok = np.isfinite(s) & (s > 0) & (d > 0)
    dd, ss = d[ok], s[ok]
    x, y = np.log(dd), np.log(ss)
    lo, hi = lo_mult * d50, HI_MULT * d50
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
    p1 = np.polyfit(x[:k + 1], y[:k + 1], 1)
    p2 = np.polyfit(x[k:], y[k:], 1)
    return {"lag_m": float(dd[k]), "lag_over_d50": float(dd[k] / d50),
            "H": float(0.5 * p1[0]), "slope2": float(p2[0]),
            "H2": float(0.5 * p2[0]), "rss": best[1],
            "n2": int(len(dd) - k)}


def nikora(d, s, d50, lo_mult, R5, label):
    """One (face, direction, window) read."""
    brk = break_search(d, s, d50, lo_mult)
    if brk is None:
        return {"label": label, "window_lo_over_d50": float(lo_mult),
                "break": None, "H_Nikora": None,
                "note": "no valid break candidate on this ladder"}
    d_min = float(d[d > 0].min())
    b = brk["lag_m"]
    h_main = R5.fit_band(d, s, d_min, b)
    up = R5.fit_band(d, s, b, TWO_TOP * d50)
    i_min_cand = min(c["index"] for c in brk["candidates"])
    i_max_cand = max(c["index"] for c in brk["candidates"])
    return {
        "label": label,
        "window_lo_over_d50": float(lo_mult),
        "window_m": brk["window_m"],
        "break": {k: brk[k] for k in ("index", "lag_m", "lag_over_d50",
                                      "rss", "n_candidates",
                                      "upper_segment_top_m")},
        "n_candidates": brk["n_candidates"],
        "candidate_index_span": [int(i_min_cand), int(i_max_cand)],
        "stuck_lower_edge": bool(brk["index"] == i_min_cand),
        "break_above_edge_over_d50": float(brk["lag_over_d50"] - lo_mult),
        "H_Nikora": h_main,
        "two_band": {"below": h_main, "above": up,
                     "below_rule": "same fit as H_Nikora",
                     "above_rule": "[break, 4 x D50]",
                     "second_segment_H2": up.get("H"),
                     "second_segment_loglog_slope": (None if up.get("H") is None
                                                     else 2.0 * up["H"])},
    }


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight + 12B reproduction
# ---------------------------------------------------------------------------
def stage0(cp):
    log("STAGE0  absolute-path preflight (the brief's faces only)")
    rows, missing = [], []
    for label, path in REQUIRED:
        ok = os.path.isfile(path)
        row = {"label": label, "path": path, "exists": bool(ok),
               "bytes": os.path.getsize(path) if ok else None,
               "sha256": sha256_file(path) if ok else None}
        rows.append(row)
        log("  %-40s %s" % (label, ("OK  %d" % row["bytes"]) if ok else "MISSING"))
        if not ok:
            missing.append(row)
    if missing:
        ctx = {"checked_utc": utcnow(), "required": rows, "missing": missing}
        write_json(PREFLIGHT, ctx)
        cp.unit_done("stage0 BLOCKED: missing %d file(s) -- STOP per brief"
                     % len(missing))
        raise Blocker("STAGE 0 preflight failed: missing %s"
                      % ", ".join(os.path.basename(m["path"]) for m in missing),
                      ctx)

    # what out/v4-2-3/ actually holds (the brief's parenthetical, checked)
    d3 = os.path.join(HERE, "out", "v4-2-3")
    v423 = sorted(n for n in os.listdir(d3)
                  if n.endswith(".npy")) if os.path.isdir(d3) else []
    path_resolution = {
        "brief_says_patch_bed_parent": "out\\v4-2-3\\",
        "npy_files_in_out_v4-2-3": v423,
        "patch_bed_used": HF_PATCH,
        "why": ("out/v4-2-3/ contains only the long bed "
                "(heightfield-long-r1.npy); the patch bed 12B measured is "
                "out/v4-2-2/heightfield.npy (run_v4_2_12b.py HF_PATCH), so "
                "'12B が使ったもの' resolves to that file."),
        "r2_face": HF_R2,
        "r2_why": "12B regenerated r2 from its seed into out/v4-2-12b/",
    }
    log("  patch face -> %s" % HF_PATCH)

    meta_p = json.load(open(META_PATCH, encoding="utf-8"))
    meta_r1 = json.load(open(META_R1, encoding="utf-8"))
    meta_r2 = json.load(open(META_R2, encoding="utf-8"))
    d50 = {"patch_m": float(meta_p["D50_surface_measured_m"]),
           "r1_m": float(meta_r1["surface_D50_m"]),
           "r2_m": float(meta_r2["surface_D50_m"])}
    geo = {"cell_m": CELL,
           "patch_shape": list(np.load(HF_PATCH, mmap_mode="r").shape),
           "r1_shape": [LONG_NY, LONG_NX],
           "r2_shape": [LONG_NY, LONG_NX],
           "d50_m": d50}

    log("STAGE0  reproduction gate: frozen window [0.25, 2] x D50, tol %.3f"
        % GATE_TOL)
    R5 = ruler()
    faces = [
        ("patch", "npy", HF_PATCH, d50["patch_m"],
         (("x", 1, DM_PATCH), ("y", 0, DM_PATCH))),
        ("long_r1", "long", HF_R1, d50["r1_m"],
         (("x", 1, DM_LONG_X), ("y", 0, DM_LONG_Y))),
        ("long_r2", "long", HF_R2, d50["r2_m"],
         (("x", 1, DM_LONG_X), ("y", 0, DM_LONG_Y))),
    ]
    curves, faces_meta, got = {}, {}, {}
    for name, kind, path, d50f, dirs in faces:
        t0 = time.time()
        z = np.load(path).astype(np.float64) if kind == "npy" else load_long(path)
        faces_meta[name] = {"path": path, "shape": list(z.shape),
                            "d50_m": float(d50f), "sigma_z_m": float(z.std()),
                            "skew": float(np.mean((z - z.mean()) ** 3)
                                          / max(float(z.std()), 1e-30) ** 3)}
        log("  face %-8s shape %s sigma_z %.6f m"
            % (name, z.shape, faces_meta[name]["sigma_z_m"]))
        for dlabel, axis, dmax in dirs:
            d, s = axis_cut(z, CELL, axis, dmax, R5)
            curves[(name, dlabel)] = (d, s)
            m = nikora(d, s, d50f, 0.25, R5, "%s | %s" % (name, dlabel))
            h = (m.get("H_Nikora") or {}).get("H")
            got["%s | %s" % (name, dlabel)] = h
            log("    %-8s %-2s H_Nikora %.4f  break %.4f mm (%.3f xD50)%s"
                % (name, dlabel, h,
                   m["break"]["lag_m"] * 1000.0, m["break"]["lag_over_d50"],
                   "  [lower-edge] " if m["stuck_lower_edge"] else ""))
        faces_meta[name]["wall_s"] = float(time.time() - t0)
        del z
        gc.collect()

    dev = {}
    for k, exp in GATE12B.items():
        g = got.get(k)
        dev[k] = {"expected_12b_brief": exp,
                  "measured_now": g,
                  "abs_diff": (None if g is None else float(abs(g - exp))),
                  "pass": bool(g is not None and abs(g - exp) <= GATE_TOL)}
    repro = {"expect_12b": GATE12B, "tol": GATE_TOL, "measured": got,
             "deviation": dev,
             "pass": bool(all(v["pass"] for v in dev.values())),
             "source_of_expectation": "BRIEF-V4-2-13 STAGE 0 (12B published)"}

    ctx = {"checked_utc": utcnow(), "brief": BRIEF, "required": rows,
           "path_resolution": path_resolution, "geometry": geo,
           "d50": d50, "faces": faces_meta, "repro": repro,
           "windows": {k: [v, HI_MULT] for k, v in WINDOWS}}
    write_json(PREFLIGHT, ctx)
    if not repro["pass"]:
        bad = ", ".join("%s %.4f vs %.4f" % (k, v["measured_now"], v["expected_12b_brief"])
                        for k, v in dev.items() if not v["pass"])
        cp.unit_done("stage0 BLOCKED: 12B reproduction failed -- STOP per brief")
        raise Blocker("STAGE 0 reproduction gate failed: %s" % bad, ctx)
    log("STAGE0  reproduction PASS (max |diff| %.5f)"
        % max(v["abs_diff"] for v in dev.values()))
    cp.unit_done("stage0")
    return ctx, curves


# ---------------------------------------------------------------------------
# STAGE 1 -- window sensitivity on the three generated faces
# ---------------------------------------------------------------------------
def stage1(cp, ctx, curves):
    log("STAGE1  lower-edge sensitivity W0-W3 on 3 faces x 2 directions")
    R5 = ruler()
    d50 = ctx["d50"]
    d50_of = {"patch": d50["patch_m"], "long_r1": d50["r1_m"],
              "long_r2": d50["r2_m"]}
    out = {
        "method": {
            "detrending": "none (raw field)",
            "axis_cuts": "run_v4_2_5.curve_dir (periodic np.roll, 60-bin "
                         "log ladder rounded to whole pixels)",
            "break_window": "lower edge = {0.25, 0.15, 0.10, 0.05} x D50; "
                            "upper edge fixed 2 x D50; candidates = ladder points",
            "break_fit": "two-segment log-log least squares; seg1 = [first ladder "
                         "lag, brk]; seg2 = [brk, 4 x D50]; >=4 points each; "
                         "minimum total residual",
            "H_Nikora": "single fit [first ladder lag (= cell), brk]",
            "second_segment": "fit_band(brk, 4 x D50); H2 = 0.5*dlnS/dlnd, "
                              "log-log slope = 2*H2",
            "stuck_rule": "break index == smallest candidate index in the window",
            "grid_mm": CELL * 1000.0,
            "unchanged_from_12b": [
                "axis cuts / ladder", "break objective and segment rules",
                "H_Nikora fit bounds (first ladder lag to break)",
                "segment-2 top = 4 x D50", "D50 per face from metadata"],
        },
        "windows": {k: [v, HI_MULT] for k, v in WINDOWS},
        "faces": {},
        "alt_reading_at_W0": {},
    }
    for name in ("patch", "long_r1", "long_r2"):
        d50f = float(d50_of[name])
        fm = {"d50_m": d50f, "sigma_z_m": ctx["faces"][name]["sigma_z_m"],
              "shape": ctx["faces"][name]["shape"], "directions": {}}
        for dlabel in ("x", "y"):
            d, s = curves[(name, dlabel)]
            rows = {}
            for wname, wlo in WINDOWS:
                m = nikora(d, s, d50f, wlo, R5, "%s | %s | %s" % (name, dlabel, wname))
                rows[wname] = m
                h = (m.get("H_Nikora") or {}).get("H")
                log("  %-8s %-2s %s H %.4f brk %.3f mm (%.3f D50)%s H2 %s"
                    % (name, dlabel, wname, h,
                       m["break"]["lag_m"] * 1000.0, m["break"]["lag_over_d50"],
                       " STUCK" if m["stuck_lower_edge"] else "",
                       "%.4f" % m["two_band"]["second_segment_H2"]
                       if m["two_band"]["second_segment_H2"] is not None else "n/a"))
            alt = alt_break_lastlag(d, s, d50f, 0.25)
            fm["directions"][dlabel] = {
                "windows": rows,
                "alt_seg2_to_last_lag_at_W0": alt,
                "curve": {"D_m": [float(v) for v in d],
                          "S_m2": [float(v) for v in s]},
            }
            h0 = rows["W0"]["H_Nikora"]["H"]
            out["alt_reading_at_W0"]["%s | %s" % (name, dlabel)] = {
                "main_break_m": rows["W0"]["break"]["lag_m"],
                "H_main": h0,
                "alt_break_m": (alt or {}).get("lag_m"),
                "H_alt": (alt or {}).get("H"),
                "dH_alt_minus_main": ((alt["H"] - h0)
                                      if alt is not None and h0 is not None
                                      else None),
            }
        out["faces"][name] = fm

    # window-wise fact table (6 face-direction readings per window)
    summ = {}
    for wname, wlo in WINDOWS:
        vals, stuck, ab = [], 0, []
        for name in out["faces"]:
            for dlabel in ("x", "y"):
                m = out["faces"][name]["directions"][dlabel]["windows"][wname]
                h = m["H_Nikora"]["H"]
                vals.append({"face": name, "dir": dlabel, "H": h,
                             "break_m": m["break"]["lag_m"],
                             "break_over_d50": m["break"]["lag_over_d50"],
                             "stuck": m["stuck_lower_edge"]})
                stuck += int(m["stuck_lower_edge"])
                ab.append(m["break_above_edge_over_d50"])
        hs = [v["H"] for v in vals]
        summ[wname] = {
            "lo_over_d50": wlo, "n_readings": len(vals),
            "n_stuck": stuck,
            "H_min": float(min(hs)), "H_max": float(max(hs)),
            "H_mean": float(np.mean(hs)),
            "diff_vs_0.79_min": float(min(0.79 - h for h in hs)),
            "diff_vs_0.79_max": float(max(0.79 - h for h in hs)),
            "break_mm_min": float(min(v["break_m"] for v in vals) * 1000.0),
            "break_mm_max": float(max(v["break_m"] for v in vals) * 1000.0),
            "break_above_edge_min_over_d50": float(min(ab)),
            "readings": vals,
        }
    out["window_summary"] = summ
    cp.unit_done("stage1")
    return out


# ---------------------------------------------------------------------------
# STAGE 2 -- known-H synthetic fields
# ---------------------------------------------------------------------------
def fbm_field(n, h, sigma, seed):
    """Isotropic periodic fBm by spectral synthesis.

    amplitude |k|^-(H+1)  <=>  2-D PSD ~ |k|^-(2H+2)  <=>  D(r) ~ r^(2H).
    Returned as float32 (the container the stored beds use).
    """
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((n, n))
    ky = np.fft.fftfreq(n)[:, None]
    kx = np.fft.rfftfreq(n)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    amp = np.zeros_like(k)
    nz = k > 0
    amp[nz] = k[nz] ** (-(h + 1.0))
    z = np.fft.irfft2(np.fft.rfft2(white) * amp, s=(n, n))
    z = z - z.mean()
    z = z / z.std() * sigma
    return z.astype(np.float32)


def stage2(cp, ctx, R5):
    log("STAGE2  known-H synthetic calibration (fBm, 2 realizations each)")
    sig = float(np.mean([ctx["faces"][n]["sigma_z_m"]
                         for n in ("patch", "long_r1", "long_r2")]))
    n = 2048
    contract = {
        "box_m": 0.6, "grid": [n, n], "cell_m": CELL,
        "cell_mm": CELL * 1000.0, "d50_m_for_windows": SYNTH_D50_M,
        "sigma_z_m": sig,
        "sigma_z_source": "mean of the three generated faces' sigma_z",
        "generator": "isotropic periodic spectral synthesis, amplitude "
                     "|k|^-(H+1) (2-D PSD ~ |k|^-(2H+2) => D(r) ~ r^(2H)); "
                     "float32 container",
        "seeds": {"H=%.2f r%d" % (h, i): s for (h, i), s in SYNTH_SEEDS.items()},
        "gen_check_band_m": list(BAND_REF),
        "windows_all": {k: [v, HI_MULT] for k, v in WINDOWS},
        "bias_definition": "estimate - H_gen",
    }
    out = {"contract": contract, "realizations": {}, "bias": {},
           "stuck_summary": {}}

    for hgen in SYNTH_H:
        for i in range(SYNTH_REAL):
            key = "H=%.2f r%d" % (hgen, i)
            seed = SYNTH_SEEDS[(hgen, i)]
            z = fbm_field(n, hgen, sig, seed).astype(np.float64)
            rec = {"H_gen": float(hgen), "realization": i, "seed": int(seed),
                   "sigma_z_m": float(z.std()),
                   "skew": float(np.mean((z - z.mean()) ** 3) / z.std() ** 3),
                   "directions": {}}
            for axis, dlabel in ((1, "x"), (0, "y")):
                d, s = axis_cut(z, CELL, axis, DM_PATCH, R5)
                chk = R5.fit_band(d, s, *BAND_REF)
                rows = {}
                for wname, wlo in WINDOWS:
                    m = nikora(d, s, SYNTH_D50_M, wlo, R5,
                               "%s | %s | %s" % (key, dlabel, wname))
                    rows[wname] = m
                    log("  %-12s %-2s %s H %.4f brk %.3f mm (%.3f D50)%s"
                        % (key, dlabel, wname, m["H_Nikora"]["H"],
                           m["break"]["lag_m"] * 1000.0,
                           m["break"]["lag_over_d50"],
                           " STUCK" if m["stuck_lower_edge"] else ""))
                rec["directions"][dlabel] = {
                    "gen_check_H": chk["H"], "gen_check_n": chk["n_points"],
                    "gen_check_band_m": list(BAND_REF),
                    "windows": rows,
                }
            out["realizations"][key] = rec
            del z
            gc.collect()

    # bias table: per (H_gen, direction, window) mean over the 2 realizations
    bias_rows = []
    for hgen in SYNTH_H:
        for dlabel in ("x", "y"):
            for wname in [w[0] for w in WINDOWS]:
                vals, stuck, brks = [], 0, []
                for i in range(SYNTH_REAL):
                    m = out["realizations"]["H=%.2f r%d" % (hgen, i)][
                        "directions"][dlabel]["windows"][wname]
                    vals.append(m["H_Nikora"]["H"])
                    stuck += int(m["stuck_lower_edge"])
                    brks.append(m["break"]["lag_over_d50"])
                bias_rows.append({
                    "H_gen": float(hgen), "dir": dlabel, "window": wname,
                    "H_mean": float(np.mean(vals)),
                    "bias_mean": float(np.mean(vals) - hgen),
                    "H_each": [float(v) for v in vals],
                    "bias_each": [float(v - hgen) for v in vals],
                    "n_stuck": int(stuck), "n": len(vals),
                    "break_over_d50_mean": float(np.mean(brks)),
                })
    out["bias_rows"] = bias_rows
    out["bias"] = {}
    for wname in [w[0] for w in WINDOWS]:
        node = {}
        for hgen in SYNTH_H:
            rs = [r for r in bias_rows
                  if r["window"] == wname and r["H_gen"] == hgen]
            bs = [b for r in rs for b in r["bias_each"]]
            node["H=%.2f" % hgen] = {
                "bias_mean_over_all_4_readings": float(np.mean(bs)),
                "bias_min": float(min(bs)), "bias_max": float(max(bs)),
                "bias_abs_max": float(max(abs(b) for b in bs)),
                "H_est_mean": float(np.mean([b + hgen for b in bs])),
                "n_stuck": int(sum(r["n_stuck"] for r in rs)),
                "n_readings": len(bs),
            }
        out["bias"][wname] = node
    for wname in [w[0] for w in WINDOWS]:
        out["stuck_summary"][wname] = int(sum(
            v["n_stuck"] for v in out["bias"][wname].values()))
    cp.unit_done("stage2")
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f(v, nd=4):
    if v is None:
        return "n/a"
    return ("%." + str(nd) + "f") % float(v)


def build_report(ctx):
    L = []
    A = L.append
    A("# REPORT-V4-2-13 — 折れ窓感度と既知 H 合成場での偏り")
    A("")
    A("- brief: `%s`" % BRIEF)
    A("- 実行: `run_v4_2_13.py`（前景、成果物を書くまで終えない）")
    A("- 作成: %s UTC" % ctx["created_utc"])
    A("- 役割: **数値を出す仕事。判定・設計への意見は書かない。**")
    A("")
    A("論点は2つ：(1) 折れ探索窓の下端を 0.25→0.15→0.10→0.05 ×D50 と下げると"
      "生成床3面の推定値と折れ位置がどう動くか、下端張り付きは解消するか。"
      "(2) 既知 H の合成場でこの物差し（各窓設定）にどれだけ偏りがあるか。")
    A("")

    # ---------------- STAGE 0
    st0 = ctx["stage0"]
    A("## STAGE 0 — 前提確認と 12B 再現")
    A("")
    A("絶対パスのみを確認した（別ツリーは探していない）。")
    A("")
    A("| 種別 | ファイル | bytes | sha256 (先頭16) |")
    A("|---|---|---:|---|")
    for r in st0["required"]:
        A("| %s | `%s` | %s | `%s` |"
          % (r["label"], os.path.basename(r["path"]),
             r["bytes"] if r["bytes"] is not None else "MISSING",
             (r["sha256"] or "—")[:16]))
    A("")
    pr = st0["path_resolution"]
    A("**パッチ床のパス解決（事実）**: brief の括弧は `out\\v4-2-3\\` と書くが、"
      "そのフォルダの .npy は `%s` のみで、12B が実際に読んだパッチ床は "
      "`out\\v4-2-2\\heightfield.npy`（`run_v4_2_12b.py` の HF_PATCH）。"
      "「12B が使ったもの」を優先し、パッチ面はそのファイルとした。"
      % (", ".join(pr["npy_files_in_out_v4-2-3"]) or "なし"))
    A("")
    A("| 面 | shape | D50 (m) | sigma_z (m) |")
    A("|---|---|---:|---:|")
    for name in ("patch", "long_r1", "long_r2"):
        fm = st0["faces"][name]
        A("| %s | %s | %.6f | %.6f |"
          % (name, "×".join(str(v) for v in fm["shape"]), fm["d50_m"],
             fm["sigma_z_m"]))
    A("")
    A("### 凍結窓 [0.25, 2]×D50 での 12B H_Nikora 再現（許容 ±%.3f）" % GATE_TOL)
    A("")
    A("| 面 | 方向 | 12B 公表値 | 本測定 | 差 | 判定 |")
    A("|---|---|---:|---:|---:|---|")
    for k in ("patch | x", "patch | y", "long_r1 | x", "long_r1 | y",
              "long_r2 | x", "long_r2 | y"):
        v = st0["repro"]["deviation"][k]
        A("| %s | %s | %s | %s | %s | %s |"
          % (k.split(" | ")[0], k.split(" | ")[1],
             f(v["expected_12b_brief"]), f(v["measured_now"]),
             f(v["abs_diff"]), "PASS" if v["pass"] else "FAIL"))
    A("")
    A("再現ゲート: **%s**（最大 |差| %s、許容 %.3f）。"
      % ("PASS" if st0["repro"]["pass"] else "FAIL",
         f(max(v["abs_diff"] for v in st0["repro"]["deviation"].values()), 5),
         GATE_TOL))
    A("")

    # ---------------- STAGE 1
    st1 = ctx["stage1"]
    A("## STAGE 1 — 窓感度（生成床 3 面 × x/y）")
    A("")
    A("下端だけを動かした。他は 12B の実装と同一：軸カット・60ビンlogラダー・"
      "2区間log-log最小残差・各区間4点以上・第2区間の上端 4×D50・"
      "H_Nikora は [ラダー最小ラグ(=格子), 折れ] の単一フィット。")
    A("")
    A("| 窓 | 下端 (×D50) |")
    A("|---|---:|")
    for wname, wlo in WINDOWS:
        A("| %s | %.2f |" % (wname, wlo))
    A("")
    A("「下端張付」= 選ばれた折れが、その窓が許す最小のラダー点であること。"
      "「第2区間 H2」= fit_band(折れ, 4×D50) の H（傾きは 2×H2）。")
    A("")

    for name in ("patch", "long_r1", "long_r2"):
        fm = st1["faces"][name]
        A("### %s（D50 = %.6f m, sigma_z = %.6f m）"
          % (name, fm["d50_m"], fm["sigma_z_m"]))
        A("")
        A("| 方向 | 窓 | H_Nikora | H_se | 折れ (mm) | 折れ (×D50) | 下端張付 | "
          "第2区間 H2 | 傾き 2H2 | n1 | n2 | 候補数 | rss |")
        A("|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|")
        for dlabel in ("x", "y"):
            for wname, wlo in WINDOWS:
                m = fm["directions"][dlabel]["windows"][wname]
                h = m["H_Nikora"]
                tb = m["two_band"]
                A("| %s | %s | %s | %s | %.4f | %.3f | %s | %s | %s | %d | %d | "
                  "%d | %s |"
                  % (dlabel, wname, f(h.get("H")), f(h.get("H_se"), 4),
                     m["break"]["lag_m"] * 1000.0,
                     m["break"]["lag_over_d50"],
                     "**yes**" if m["stuck_lower_edge"] else "no",
                     f(tb.get("second_segment_H2")),
                     f(tb.get("second_segment_loglog_slope")),
                     h.get("n_points", 0),
                     (tb.get("above") or {}).get("n_points", 0),
                     m["n_candidates"], f(m["break"]["rss"], 3)))
        A("")

    A("### 窓ごとの要約（6 読み = 3 面 × x/y）")
    A("")
    A("| 窓 | 下端 ×D50 | 張付数 | H 最小 | H 最大 | H 平均 | "
      "0.79−H 最小 | 0.79−H 最大 | 折れ mm 最小 | 折れ mm 最大 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for wname, wlo in WINDOWS:
        s = st1["window_summary"][wname]
        A("| %s | %.2f | %d/6 | %s | %s | %s | %s | %s | %.3f | %.3f |"
          % (wname, wlo, s["n_stuck"], f(s["H_min"]), f(s["H_max"]),
             f(s["H_mean"]), f(s["diff_vs_0.79_min"]), f(s["diff_vs_0.79_max"]),
             s["break_mm_min"], s["break_mm_max"]))
    A("")

    A("### 12B の別読み（第2区間をラダー末端まで）— W0 で再掲")
    A("")
    A("| 面 | 方向 | 本手順 折れ (mm) | 本手順 H | 別読み 折れ (mm) | 別読み H | 差 |")
    A("|---|---|---:|---:|---:|---:|---:|")
    for k in ("patch | x", "patch | y", "long_r1 | x", "long_r1 | y",
              "long_r2 | x", "long_r2 | y"):
        v = st1["alt_reading_at_W0"][k]
        A("| %s | %s | %.4f | %s | %.4f | %s | %s |"
          % (k.split(" | ")[0], k.split(" | ")[1],
             v["main_break_m"] * 1000.0, f(v["H_main"]),
             (v["alt_break_m"] or 0) * 1000.0, f(v["H_alt"]),
             f(v["dH_alt_minus_main"], 4)))
    A("")

    # ---------------- STAGE 2
    st2 = ctx["stage2"]
    c = st2["contract"]
    A("## STAGE 2 — 既知 H 合成場での偏り較正")
    A("")
    A("| 項目 | 値 |")
    A("|---|---|")
    A("| 箱 | %.1f m × %.1f m、%d×%d、格子 %.5f mm |"
      % (c["box_m"], c["box_m"], c["grid"][0], c["grid"][1], c["cell_mm"]))
    A("| 窓換算の D50 | %.4f m（brief の換算値、生成床と同じ） |" % c["d50_m_for_windows"])
    A("| sigma_z | %.6f m（生成床3面の平均） |" % c["sigma_z_m"])
    A("| 生成法 | 等方・周期スペクトル合成、振幅 \\|k\\|^-(H+1)"
      "（2次元PSD ~ \\|k\\|^-(2H+2)、D(r) ~ r^(2H)）、float32 |")
    A("| 実現数 | H ごとに 2（seed 明示） |")
    A("")
    A("生成した場が狙った H を持つことの確認。場は構成上 H_gen を持つ"
      "（振幅 \\|k\\|^-(H+1)）。下は**物差し自身**で 4–32 px（%.1f–%.1f mm）を"
      "測った値で、ここにも格子離散化・有限領域の読みが乗る。"
      % (BAND_REF[0] * 1000.0, BAND_REF[1] * 1000.0))
    A("")
    A("| H_gen | 実現 | seed | sigma_z (m) | 4–32 px H (x) | 4–32 px H (y) |")
    A("|---|---:|---:|---:|---:|---:|")
    for hgen in SYNTH_H:
        for i in range(SYNTH_REAL):
            k = "H=%.2f r%d" % (hgen, i)
            r = st2["realizations"][k]
            A("| %.2f | %d | %d | %.6f | %s | %s |"
              % (hgen, i, r["seed"], r["sigma_z_m"],
                 f(r["directions"]["x"]["gen_check_H"]),
                 f(r["directions"]["y"]["gen_check_H"])))
    A("")

    A("### 実現ごとの測定（H_gen / 方向 / 窓）")
    A("")
    A("| H_gen | 実現 | 方向 | 窓 | H_Nikora | 偏り | 折れ (mm) | 折れ (×D50) | 下端張付 | 第2区間 H2 |")
    A("|---|---:|---|---|---:|---:|---:|---:|---|---:|")
    for hgen in SYNTH_H:
        for i in range(SYNTH_REAL):
            k = "H=%.2f r%d" % (hgen, i)
            for dlabel in ("x", "y"):
                for wname, wlo in WINDOWS:
                    m = st2["realizations"][k]["directions"][dlabel]["windows"][wname]
                    A("| %.2f | %d | %s | %s | %s | %+0.4f | %.4f | %.3f | %s | %s |"
                      % (hgen, i, dlabel, wname, f(m["H_Nikora"]["H"]),
                         m["H_Nikora"]["H"] - hgen,
                         m["break"]["lag_m"] * 1000.0,
                         m["break"]["lag_over_d50"],
                         "**yes**" if m["stuck_lower_edge"] else "no",
                         f(m["two_band"].get("second_segment_H2"))))
    A("")

    A("### 偏り表（推定 − 真値、2 実現 × x/y = 4 読みの平均）")
    A("")
    A("| 窓 | 下端 ×D50 | " + " | ".join("H=%.2f 偏り" % h for h in SYNTH_H)
      + " | " + " | ".join("H=%.2f 推定平均" % h for h in SYNTH_H) + " | 張付数 |")
    A("|---|---:|" + "---:|" * len(SYNTH_H) + "---:|" * len(SYNTH_H) + "---:|")
    for wname, wlo in WINDOWS:
        node = st2["bias"][wname]
        cells = []
        for hgen in SYNTH_H:
            v = node["H=%.2f" % hgen]
            cells.append("%+0.4f（%+0.4f〜%+0.4f）"
                         % (v["bias_mean_over_all_4_readings"],
                            v["bias_min"], v["bias_max"]))
        means = [f(node["H=%.2f" % h]["H_est_mean"]) for h in SYNTH_H]
        A("| %s | %.2f | %s | %s | %d/12 |"
          % (wname, wlo, " | ".join(cells), " | ".join(means),
             st2["stuck_summary"][wname]))
    A("")
    A("（張付数は 3 H_gen × 2 実現 × 2 方向 = 12 読み中の数。）")
    A("")

    # ---------------- facts
    A("## 事実のまとめ")
    A("")
    w0 = st1["window_summary"]["W0"]
    w3 = st1["window_summary"]["W3"]
    A("1. **窓を下げたときの H_Nikora の動き。** 6 読みの平均は W0 %s → W1 %s → "
      "W2 %s → W3 %s。面・方向ごとの W0→W3 の差は次のとおり。"
      % (f(w0["H_mean"]), f(st1["window_summary"]["W1"]["H_mean"]),
         f(st1["window_summary"]["W2"]["H_mean"]), f(w3["H_mean"])))
    A("")
    A("   | 面 | 方向 | W0 | W1 | W2 | W3 | W3−W0 |")
    A("   |---|---|---:|---:|---:|---:|---:|")
    for name in ("patch", "long_r1", "long_r2"):
        for dlabel in ("x", "y"):
            ws = st1["faces"][name]["directions"][dlabel]["windows"]
            h0 = ws["W0"]["H_Nikora"]["H"]
            h3 = ws["W3"]["H_Nikora"]["H"]
            A("   | %s | %s | %s | %s | %s | %s | %+0.4f |"
              % (name, dlabel, f(h0), f(ws["W1"]["H_Nikora"]["H"]),
                 f(ws["W2"]["H_Nikora"]["H"]), f(h3), h3 - h0))
    A("")
    A("2. **下端張り付き。** W0 で %d/6、W1 で %d/6、W2 で %d/6、W3 で %d/6。"
      % (w0["n_stuck"], st1["window_summary"]["W1"]["n_stuck"],
         st1["window_summary"]["W2"]["n_stuck"], w3["n_stuck"]))
    A("   折れ位置の範囲は W0 %.3f–%.3f mm → W3 %.3f–%.3f mm。"
      % (w0["break_mm_min"], w0["break_mm_max"],
         w3["break_mm_min"], w3["break_mm_max"]))
    A("")
    A("3. **0.79 との差。** 生成床 6 読みの (0.79 − H_Nikora) は、"
      "W0 で %.3f–%.3f、W1 で %.3f–%.3f、W2 で %.3f–%.3f、W3 で %.3f–%.3f。"
      % (w0["diff_vs_0.79_min"], w0["diff_vs_0.79_max"],
         st1["window_summary"]["W1"]["diff_vs_0.79_min"],
         st1["window_summary"]["W1"]["diff_vs_0.79_max"],
         st1["window_summary"]["W2"]["diff_vs_0.79_min"],
         st1["window_summary"]["W2"]["diff_vs_0.79_max"],
         w3["diff_vs_0.79_min"], w3["diff_vs_0.79_max"]))
    A("   全 24 読み（6 面方向 × 4 窓）で差が 0.25 以上だったか: **%s**"
      "（最小 %s）。"
      % ("はい" if min(s["diff_vs_0.79_min"] for s in st1["window_summary"].values())
         >= 0.25 else "いいえ",
         f(min(s["diff_vs_0.79_min"] for s in st1["window_summary"].values()))))
    A("")
    bs0 = st2["bias"]["W0"]
    bs3 = st2["bias"]["W3"]
    A("4. **既知 H の偏り。** 4 読み平均の偏り（推定−真値）は、"
      "W0 で H=0.50 %+0.4f / 0.65 %+0.4f / 0.79 %+0.4f、"
      "W3 で H=0.50 %+0.4f / 0.65 %+0.4f / 0.79 %+0.4f。"
      % (bs0["H=0.50"]["bias_mean_over_all_4_readings"],
         bs0["H=0.65"]["bias_mean_over_all_4_readings"],
         bs0["H=0.79"]["bias_mean_over_all_4_readings"],
         bs3["H=0.50"]["bias_mean_over_all_4_readings"],
         bs3["H=0.65"]["bias_mean_over_all_4_readings"],
         bs3["H=0.79"]["bias_mean_over_all_4_readings"]))
    A("   偏りの絶対値の最大は W0 %s、W1 %s、W2 %s、W3 %s（窓ごと・全 H の全読み）。"
      % (f(max(v["bias_abs_max"] for v in bs0.values())),
         f(max(v["bias_abs_max"] for v in st2["bias"]["W1"].values())),
         f(max(v["bias_abs_max"] for v in st2["bias"]["W2"].values())),
         f(max(v["bias_abs_max"] for v in bs3.values()))))
    A("")
    A("5. **合成場の下端張り付き。** W0 %d/12、W1 %d/12、W2 %d/12、W3 %d/12。"
      % (st2["stuck_summary"]["W0"], st2["stuck_summary"]["W1"],
         st2["stuck_summary"]["W2"], st2["stuck_summary"]["W3"]))
    A("")
    A("6. **第2区間の傾き（生成床）。** 6 読みの第2区間 H2 は、"
      "W0 で %.3f–%.3f、W3 で %.3f–%.3f。"
      % (min(st1["faces"][n]["directions"][d]["windows"]["W0"]["two_band"]["second_segment_H2"]
             for n in st1["faces"] for d in ("x", "y")),
         max(st1["faces"][n]["directions"][d]["windows"]["W0"]["two_band"]["second_segment_H2"]
             for n in st1["faces"] for d in ("x", "y")),
         min(st1["faces"][n]["directions"][d]["windows"]["W3"]["two_band"]["second_segment_H2"]
             for n in st1["faces"] for d in ("x", "y")),
         max(st1["faces"][n]["directions"][d]["windows"]["W3"]["two_band"]["second_segment_H2"]
             for n in st1["faces"] for d in ("x", "y"))))
    A("")

    # ---------------- appendix
    A("## 付録 — 実行と成果物")
    A("")
    A("- 書き込み先: `%s` と `out/v4-2-13/`" % THIS_FILE)
    A("- 環境: Python %s / numpy %s / %s"
      % (ctx["env"]["python"], ctx["env"]["numpy"], ctx["env"]["platform"]))
    A("- 所要: " + ", ".join("%s %.1fs" % (k, v)
                            for k, v in (ctx.get("timings") or {}).items()))
    ck = ctx.get("checkpoint_units") or {}
    A("- checkpoint: `%s` — units %s（last: %s）"
      % (ck.get("path"), ck.get("ids"), ck.get("last_unit")))
    A("- 数値: `out/v4-2-13/numbers-v4-2-13.json`（全表・全曲線を含む）")
    A("")
    A("### 数値の読み方")
    A("")
    A("- H_Nikora は 12B と同じ定義（[ラダー最小ラグ, 折れ] の単一 log-log フィット）"
      "であり、窓を変えてもフィット下端は動かない。動くのは折れ候補の下限だけ。")
    A("- 折れは 60 ビンの log ラダー上の離散点なので、窓の下端を下げても"
      "候補点はラダー点の間隔ぶん飛ぶ。")
    A("- STAGE 2 の fBm は粒径スケールのカットオフを持たない（純 fBm）。"
      "生成床が持つ粒スケールの平坦化は含まれない。")
    A("")
    return "\n".join(L) + "\n"


def build_stop_report(ctx, reason):
    L = ["# REPORT-V4-2-13 — STOP record", "",
         "- brief: `%s`" % BRIEF,
         "- 実行: `run_v4_2_13.py`",
         "- 作成: %s UTC" % ctx.get("created_utc", utcnow()), "",
         "> **この実行は `%s` で停止した。**" % reason, "",
         "> 停止時点で確認・測定できた値だけを残す。", ""]
    st0 = ctx.get("stage0") or {}
    if st0.get("required"):
        L.append("## STAGE 0 — 確認値")
        L.append("")
        L.append("| 種別 | ファイル | bytes | sha256 (先頭16) |")
        L.append("|---|---|---:|---|")
        for r in st0["required"]:
            L.append("| %s | `%s` | %s | `%s` |"
                     % (r["label"], os.path.basename(r["path"]),
                        r["bytes"] if r["bytes"] is not None else "MISSING",
                        (r["sha256"] or "—")[:16]))
        L.append("")
        if st0.get("missing"):
            L.append("欠損: " + ", ".join("`%s`" % m["path"]
                                          for m in st0["missing"]))
            L.append("")
        if st0.get("repro"):
            L.append("再現ゲート: **%s**"
                     % ("PASS" if st0["repro"]["pass"] else "FAIL"))
            L.append("")
            L.append("| 面 | 方向 | 12B | 本測定 | 差 |")
            L.append("|---|---|---:|---:|---:|")
            for k, v in st0["repro"]["deviation"].items():
                L.append("| %s | %s | %s | %s | %s |"
                         % (k.split(" | ")[0], k.split(" | ")[1],
                            f(v["expected_12b_brief"]), f(v["measured_now"]),
                            f(v["abs_diff"])))
            L.append("")
    L.append("STAGE 1 / STAGE 2 は未測定。")
    L.append("")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
def _ckpt_units(cp):
    st = cp.state() or {}
    ids = st.get("completed_ids") or []
    return {"path": CP_PATH, "n": len(ids), "ids": ids,
            "completed_units": st.get("completed_units"),
            "last_unit": st.get("last_unit")}


def main():
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    created = utcnow()
    timings = {}
    ctx = {"created_utc": created, "brief": BRIEF,
           "env": {"python": sys.version.split()[0], "numpy": np.__version__,
                   "platform": sys.platform},
           "stages_done": []}

    def stop_record(reason):
        ctx["stopped_at"] = reason
        cp.unit_done("report (STOP: %s)" % reason)
        ctx["checkpoint_units"] = _ckpt_units(cp)
        md = build_stop_report(ctx, reason)
        atomic_write_bytes(REPORT, md.encode("utf-8"),
                           verify=lambda p: os.path.getsize(p) > 500)
        write_json(NUMBERS, {
            "brief": BRIEF, "created_utc": created, "stopped_at": reason,
            "stage0": ctx.get("stage0"), "stage1": ctx.get("stage1"),
            "stage2": ctx.get("stage2"), "env": ctx["env"], "log": LOG})
        log("STOPPED  %s" % reason)
        log("      %s" % REPORT)

    t0 = time.time()
    try:
        st0, curves = stage0(cp)
    except Blocker as b:
        timings["stage0"] = time.time() - t0
        ctx["stage0"] = b.ctx
        ctx["timings"] = timings
        stop_record("STAGE 0: %s" % b.reason)
        return 2
    timings["stage0"] = time.time() - t0
    ctx["stage0"] = st0
    ctx["d50"] = st0["d50"]
    ctx["faces"] = st0["faces"]
    ctx["stages_done"].append("stage0")
    log("STAGE0  done in %.1f s" % timings["stage0"])

    R5 = ruler()
    t0 = time.time()
    ctx["stage1"] = stage1(cp, ctx, curves)
    timings["stage1"] = time.time() - t0
    ctx["stages_done"].append("stage1")
    log("STAGE1  done in %.1f s" % timings["stage1"])

    t0 = time.time()
    ctx["stage2"] = stage2(cp, ctx, R5)
    timings["stage2"] = time.time() - t0
    ctx["stages_done"].append("stage2")
    log("STAGE2  done in %.1f s" % timings["stage2"])

    ctx["timings"] = timings
    # the report names the checkpoint state as it will be once the report unit
    # is committed; the file is written before that commit (fail-closed order)
    ids = (cp.state() or {}).get("completed_ids") or []
    ctx["checkpoint_units"] = {"path": CP_PATH, "ids": ids + ["report"],
                               "last_unit": "report",
                               "completed_units": len(ids) + 1,
                               "note": "report unit committed right after the "
                                       "file is written"}
    md = build_report(ctx)
    atomic_write_bytes(REPORT, md.encode("utf-8"),
                       verify=lambda p: os.path.getsize(p) > 6000)
    cp.unit_done("report")
    ctx["checkpoint_units"] = _ckpt_units(cp)

    numbers = {
        "brief": BRIEF, "created_utc": created,
        "what_this_is": (
            "break-window lower-edge sensitivity of the Nikora ruler on the "
            "three generated faces (STAGE 1), plus a known-H fBm bias "
            "calibration of the same ruler at each window (STAGE 2). "
            "Numbers only: no verdict, no design opinion, no parameter change."),
        "write_set": [THIS_FILE, "out/v4-2-13/"],
        "env": ctx["env"],
        "windows": {k: [v, HI_MULT] for k, v in WINDOWS},
        "stage0": ctx["stage0"],
        "stage1": ctx["stage1"],
        "stage2": ctx["stage2"],
        "checkpoint": {"path": CP_PATH, "state": ctx["checkpoint_units"]},
        "deliverables": {"report": REPORT, "numbers": NUMBERS,
                         "preflight": PREFLIGHT, "checkpoint": CP_PATH},
        "timings_s": {k: round(v, 1) for k, v in timings.items()},
        "log": LOG,
    }
    write_json(NUMBERS, numbers)
    log("DONE  %s" % REPORT)
    log("      %s" % NUMBERS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
