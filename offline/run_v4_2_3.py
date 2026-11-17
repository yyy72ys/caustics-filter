"""run_v4_2_3 -- BRIEF-V4-2-3: is the H failure the generator, or the ruler
(patch size)?  AUDIT ONLY -- no rendering, no transport, no materials.

New files only: this driver and out/v4-2-3/.  `bed_gen2.py` is imported and
used unmodified (its own functions, its own constants).  The construction
parameters are read out of out/v4-2-2/numbers-v4-2-2.json and asserted equal
to the bed_gen2 module constants before anything is built.

What the brief asks, and how each item is answered here
-------------------------------------------------------
1. A long bed along the flow (3.0 m x 0.6 m) at the v4-2.2 cell size
   (0.29296875 mm/px exactly, i.e. 10240 x 2048).  It is a genuine 2-D field:
   the transect-column fallback of the brief is NOT used.  The substrate is
   the v4.2.2 construction generalised to the rectangular periodic box
   (same fBm law, same three low-mode wavelengths, same skew transform, same
   sigma_form); the grains are v4.2.2's own population routine run five times
   side by side (one 0.6 m tile per 0.6 m of the long bed, independent seeds).
   `--equiv` proves the generalisation is bit-compatible with v4-2.2: on a
   0.6 x 0.6 box with the v4-2.2 mode set it reproduces bed_gen2.form_field
   and bed_gen.fbm_field.
2. The SAME estimator as v4-2.2 (`bed_gen2.structure_audit`, the gate
   function) on the 4-20 cm band, plus an identical audit on a longer lag
   range (up to 1/3 of the long bed = 1.0 m) so that the lag-dependent curve
   itself is reported.
3. The control: the substrate alone (fBm H=0.79 + low modes + skew, v4-2.2
   composition) measured with the same estimator and the same shapes.
4. The v4-2.2 0.6 m patch re-measured with the same code, so that the only
   thing that differs between the two beds is the patch dimension.

Stages (all foreground, idempotent, one jobkit checkpoint unit per stage):
  --stage0   absolute-path preflight + sha256 of every existing .py
  --equiv    the generalisation is bit-compatible with v4-2-2
  --bed      v4-2.2 patch controls + the long bed(s)
  --h        all structure-function curves and audits -> h-curves.json
  --gates    G-LEGACY (sha256 of every .py against the stage-0 record)
  --report   numbers-v4-2-3.json + REPORT-V4-2-3.md
  --all      the above, in order, in this one process
"""
from __future__ import annotations

import argparse
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

import bed_gen as BG                                    # noqa: E402
import bed_gen2 as B2                                   # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-3.md"
OUT = os.path.join(HERE, "out", "v4-2-3")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-3.json")
CURVES = os.path.join(OUT, "h-curves.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-3.md")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
TIMINGS = os.path.join(OUT, "timings.json")
LONG_HF = os.path.join(OUT, "heightfield-long-r1.npy")
LONG_META = os.path.join(OUT, "longbed-meta.json")

NUMBERS_422 = os.path.join(HERE, "out", "v4-2-2", "numbers-v4-2-2.json")
REPORT_422 = os.path.join(HERE, "out", "v4-2-2", "REPORT-V4-2-2.md")
HEIGHTFIELD_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")

# stage-0 preflight (BRIEF STAGE 0): these three must exist or we stop
REQUIRED = {
    "bed_gen2.py": os.path.join(HERE, "bed_gen2.py"),
    "run_v4_2_2.py": os.path.join(HERE, "run_v4_2_2.py"),
    "out/v4-2-2/REPORT-V4-2-2.md": REPORT_422,
}
# extra inputs this brief leans on (checked too, same stop rule)
EXTRA_REQUIRED = {
    "out/v4-2-2/numbers-v4-2-2.json": NUMBERS_422,
    "out/v4-2-2/heightfield.npy": HEIGHTFIELD_422,
}
ALLOWED_NEW = {"run_v4_2_3.py"}
THIS_FILE = "run_v4_2_3.py"

# --- the long bed (brief: 3.0 m x 0.6 m at the v4-2.2 cell size) -----------
PATCH = B2.PATCH                       # 0.6 m, from the module (not retyped)
CELL = PATCH / B2.HIRES                # 0.00029296875 m/px, exactly v4-2.2
LONG_LX = 3.0
LONG_LY = PATCH
LONG_NX = int(round(LONG_LX / CELL))   # 10240
LONG_NY = int(round(LONG_LY / CELL))   # 2048
assert abs(LONG_NX * CELL - LONG_LX) < 1e-12
assert abs(LONG_NY * CELL - LONG_LY) < 1e-12
N_TILES = int(round(LONG_LX / PATCH))  # 5 verbatim 0.6 m tiles along flow
FORM_N = 512                           # v4.2.2 substrate grid (1.171875 mm)
FORM_NX = int(round(LONG_LX / (PATCH / FORM_N)))   # 2560
FORM_NY = FORM_N

# v4.2.2 low-mode wavevectors: |k|^2 in {1,2,4}*(2*pi/PATCH)^2 -> periods
# 0.6 / 0.4243 / 0.3 m (12 signed modes).  On a box sx patches long the same
# physical wavevectors sit at lattice index (sx*i, j), so the mode set is
# re-indexed, not re-chosen: same radii, same weights.
LOW_KMAX = 2                            # v4-2.2: |k| <= 2/PATCH
SX = int(round(LONG_LX / PATCH))        # 5
assert SX == N_TILES

SEED = B2.SEED
# two independent realisations of the long bed; declared before measuring
LONG_SEEDS = [int(SEED), int(SEED) + 500000]
FBM_SEEDS = [int(SEED) + 11 * k + 1000 for k in range(4)]

GATE_BAND = (0.70, 0.85)               # not moved
GATE_LAG_BAND = (0.04, 0.20)           # the v4-2.2 fit band, not moved


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")


def log(msg):
    print("[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg),
          flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError("not JSON serialisable: %r" % (type(o),))


def jdump(path, obj):
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False, indent=2,
                                        default=_json_default)
                        .encode("utf-8"))


def jload(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def py_hashes():
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p):
            out[name] = sha256_file(p)
    return out


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def preflight():
    log("STAGE 0 preflight (absolute paths)")
    missing = []
    for label, p in list(REQUIRED.items()) + list(EXTRA_REQUIRED.items()):
        ok = os.path.isfile(p)
        log("  %-34s %s  %s" % (label, "OK " if ok else "MISSING", p))
        if not ok:
            missing.append((label, p))
    return missing


def stage0(cp):
    missing = preflight()
    if missing:
        log("STOP: required absolute paths are missing -- reporting instead of "
            "running")
        jdump(os.path.join(OUT, "stage0-missing.json"),
              {"created_utc": now_iso(), "missing": missing})
        raise SystemExit(2)
    if os.path.isfile(STAGE0_HASHES):
        doc = jload(STAGE0_HASHES)
        log("  golden stage-0 baseline already exists (created %s, %d .py) -- "
            "kept, not rewritten" % (doc.get("created_utc"),
                                     doc.get("count", -1)))
        if cp:
            cp.unit_done("stage0")
        return doc
    h = py_hashes()
    doc = {
        "created_utc": now_iso(),
        "workdir": HERE,
        "count": len(h),
        "hashes": h,
        "note": ("captured at the start of the v4-2-3 run for the G-LEGACY "
                 "audit; this file is the golden baseline and is never "
                 "rewritten.  %s is this run's own file: it is in the record "
                 "for completeness but G-LEGACY classifies it separately "
                 "(see numbers-v4-2-3.json -> gates)." % THIS_FILE),
        "new_files_this_run": sorted(ALLOWED_NEW),
        "required_present": {k: True for k in list(REQUIRED) +
                             list(EXTRA_REQUIRED)},
    }
    jdump(STAGE0_HASHES, doc)
    log("  stage-0 sha256 for %d existing .py -> %s" % (len(h), STAGE0_HASHES))
    if cp:
        cp.unit_done("stage0")
    return doc


# ---------------------------------------------------------------------------
# the v4-2.2 configuration, read verbatim out of numbers-v4-2-2.json
# ---------------------------------------------------------------------------
def load_cfg():
    N = jload(NUMBERS_422)
    bs = N["bed_stats"]
    form = bs["form"]
    prot = bs["protrusion_model"]
    fin = bs["sigma_z_fixpoint_final_parameters"]
    cfg = {
        "source": os.path.relpath(NUMBERS_422, HERE),
        "patch_m": float(bs["patch_m"]),
        "D50_in_m": float(bs["D50_in_m"]),
        "sigma_g": float(bs["sigma_g"]),
        "armor_ratio_input": float(bs["armor_ratio_input"]),
        "D50_surface_nominal_m": float(bs["D50_surface_nominal_m"]),
        "D84_surface_nominal_m": float(bs["D84_surface_nominal_m"]),
        "D50_surface_measured_m": float(bs["D50_surface_measured_m"]),
        "sigma_form_m": float(form["sigma_form_m"]),
        "H_gen": float(form["H_gen"]),
        "w_low": float(form["w_low"]),
        "skew_c": float(form["skew_c"]),
        "p_scale": float(fin["p_scale"]),
        "p_tail_prob": float(prot["p_tail_prob"]),
        "p_tail_mean_over_D": float(prot["p_tail_mean_over_D"]),
        "p_tail_max_over_D": float(prot["p_tail_max_over_D"]),
        "p_bury_range_over_D": [float(v) for v in prot["p_bury_range_over_D"]],
        "size_sorting_gain": float(prot["size_sorting_gain"]),
        "size_sorting_q": float(prot["size_sorting_q"]),
        "surface_coverage_target": float(prot["surface_coverage_target"]),
        "bulk_solid_volume_factor": float(prot["bulk_solid_volume_factor"]),
        "D_clip_max_m": float(prot["D_clip_max_m"]),
        "seed": int(SEED),
        "n_grains_v422": int(bs["n_grains"]),
        "sigma_z_matrix_v422_m": float(bs["sigma_z_matrix_m"]),
        "skew_matrix_v422": float(bs["skew_matrix"]),
        "H_gate_v422": float(N["g_bed2"]["structure_H"]["value"]),
        "hires_cell_m": float(bs["raycast"]["hires_cell_m"]),
    }
    # the json must agree with the module constants we are about to use
    checks = {
        "patch_m == bed_gen2.PATCH": cfg["patch_m"] == B2.PATCH,
        "D50_in == bed_gen2.D50_IN": cfg["D50_in_m"] == B2.D50_IN,
        "sigma_g == bed_gen2.SIGMA_G": cfg["sigma_g"] == B2.SIGMA_G,
        "armor == bed_gen2.ARMOR":
            cfg["armor_ratio_input"] == B2.ARMOR,
        "H_gen == bed_gen2.SUBSTRATE_H": cfg["H_gen"] == B2.SUBSTRATE_H,
        "w_low == bed_gen2.FORM_W_LOW": cfg["w_low"] == B2.FORM_W_LOW,
        "skew_c == bed_gen2.FORM_SKEW_C": cfg["skew_c"] == B2.FORM_SKEW_C,
        "p_scale == bed_gen2.P_SCALE0": cfg["p_scale"] == B2.P_SCALE0,
        "p_tail == bed_gen2.P_TAIL": cfg["p_tail_prob"] == B2.P_TAIL,
        "p_tail_mean == bed_gen2.P_TAIL_MEAN":
            cfg["p_tail_mean_over_D"] == B2.P_TAIL_MEAN,
        "p_tail_max == bed_gen2.P_TAIL_MAX":
            cfg["p_tail_max_over_D"] == B2.P_TAIL_MAX,
        "bury == bed_gen2.P_BURY_*":
            cfg["p_bury_range_over_D"] == [B2.P_BURY_MIN, B2.P_BURY_MAX],
        "sort_gain == bed_gen2.SORT_GAIN":
            cfg["size_sorting_gain"] == B2.SORT_GAIN,
        "sort_q == bed_gen2.SORT_Q": cfg["size_sorting_q"] == B2.SORT_Q,
        "surf_cover == bed_gen2.SURF_COVER":
            cfg["surface_coverage_target"] == B2.SURF_COVER,
        "bulk_vol == bed_gen2.BULK_VOL":
            cfg["bulk_solid_volume_factor"] == B2.BULK_VOL,
        "D_clip == bed_gen2.D_CLIP_MAX":
            cfg["D_clip_max_m"] == B2.D_CLIP_MAX,
        "seed == bed_gen2.SEED": cfg["seed"] == B2.SEED,
        "cell == 0.6/2048": cfg["hires_cell_m"] == CELL,
    }
    bad = [k for k, v in checks.items() if not v]
    if bad:
        raise SystemExit("config mismatch between numbers-v4-2-2.json and "
                         "bed_gen2.py: %r" % bad)
    cfg["assertions"] = {k: bool(v) for k, v in checks.items()}
    return cfg


# ---------------------------------------------------------------------------
# substrate: v4.2.2 construction on a rectangular periodic box
# ---------------------------------------------------------------------------
def fbm_rect(nx, ny, Lx, Ly, H, seed, sigma=1.0):
    """bed_gen.fbm_field, verbatim arithmetic, generalised to nx != ny.

    bed_gen.fbm_field filters white noise with A(k) = k^-(H+1) on a square
    n x n grid of side L; this is the same line for line with the two axes
    allowed to differ.  On a square grid it is bit-identical (checked in
    --equiv).
    """
    rng = np.random.default_rng(int(seed))
    w = rng.standard_normal((ny, nx))
    kx = np.fft.fftfreq(nx, d=Lx / nx) * 2.0 * math.pi
    ky = np.fft.fftfreq(ny, d=Ly / ny) * 2.0 * math.pi
    KX, KY = np.meshgrid(kx, ky)
    K = np.hypot(KX, KY)
    amp = np.zeros_like(K)
    m = K > 0
    amp[m] = K[m] ** (-(H + 1.0))
    Z = np.fft.fft2(w) * amp
    z = np.real(np.fft.ifft2(Z))
    z -= z.mean()
    s = z.std()
    return (z * (sigma / s)).astype(np.float64)


def low_mode_rect(nx, ny, Lx, Ly, H, seed, kmax=2, sx=1):
    """bed_gen2.low_mode_field with an x-scale factor on the box.

    The v4-2.2 mode set is |k|^2 in {1,2,4} * (2*pi/PATCH)^2, i.e. the 12
    signed (i,j) with i^2+j^2 in {1,2,4} on a PATCH x PATCH lattice.  On a
    box that is `sx` patches long, the *same physical wavevectors* live at
    lattice index (sx*i, j): |k|^2 = (2pi/(sx*PATCH))^2 (sx*i)^2 +
    (2pi/PATCH)^2 j^2 = (2pi/PATCH)^2 (i^2+j^2).  So the three wavelengths
    0.6 / 0.4243 / 0.3 m and their relative weights k^-(H+1)/2 are unchanged,
    exactly, and sx=1 reduces to bed_gen2.low_mode_field (checked in --equiv).
    Same law, same phase draw, same order, same unit-std normalisation.
    """
    rng = np.random.default_rng(int(seed))
    x = (np.arange(nx) + 0.5) * (Lx / nx)
    y = (np.arange(ny) + 0.5) * (Ly / ny)
    X, Y = np.meshgrid(x, y)
    z = np.zeros((ny, nx), dtype=np.float64)
    n_used = 0
    for i0 in range(-kmax, kmax + 1):
        for j0 in range(-kmax, kmax + 1):
            k2 = i0 * i0 + j0 * j0
            if k2 == 0 or k2 > kmax * kmax:
                continue
            a = k2 ** (-(H + 1.0) / 2.0)
            ph = rng.uniform(0.0, 2.0 * math.pi)
            i, j = sx * i0, j0
            z += a * np.cos(2.0 * math.pi * (i * X / Lx + j * Y / Ly) + ph)
            n_used += 1
    z -= z.mean()
    return z / max(z.std(), 1e-12), n_used


def form_long(nx, ny, Lx, Ly, sigma, seed, w_low, skew_c, H, sx=1, kmax=2):
    """bed_gen2.form_field on a rectangular box (same composition verbatim).

    v = w_low*u + sqrt(1-w_low^2)*f ; v <- v + c(v^2-1) ; scale to `sigma`.
    `sx` = Lx/PATCH (integer); it only re-indexes the low-mode lattice.
    """
    f = fbm_rect(nx, ny, Lx, Ly, H, seed + 7, 1.0)
    u, n_modes = low_mode_rect(nx, ny, Lx, Ly, H, seed + 11, kmax=kmax, sx=sx)
    g = w_low * u + math.sqrt(max(1.0 - w_low * w_low, 0.0)) * f
    g = g / max(g.std(), 1e-12)
    v = g + skew_c * (g * g - 1.0)
    v -= v.mean()
    return (v * (sigma / max(v.std(), 1e-12))).astype(np.float64), n_modes


def form_at_grid_long(form, Lx, Ly, cell, nx, ny, chunk=256):
    """bed_gen2.form_at_grid, generalised (bilinear, periodic, cell-centres).

    Done in row blocks so that the 3.0 m statistics grid (2048 x 10240 float64)
    does not need six simultaneous full-size temporaries.  The arithmetic per
    element is the same expression as bed_gen._bilinear, so the values are
    identical to the unchunked call (checked in --equiv on the 0.6 m box).
    """
    cellf_x = Lx / form.shape[1]
    cellf_y = Ly / form.shape[0]
    xs = (np.arange(nx) + 0.5) * cell
    out = np.empty((ny, nx), dtype=np.float64)
    px = (xs / cellf_x)[None, :]
    for a in range(0, ny, chunk):
        b = min(a + chunk, ny)
        ys = (np.arange(a, b) + 0.5) * cell
        out[a:b, :] = BG._bilinear(form, px, (ys / cellf_y)[:, None])
    return out


# ---------------------------------------------------------------------------
# grains: v4.2.2's own population routine, laid down tile by tile
# ---------------------------------------------------------------------------
def populate_long(seed, n_tiles=N_TILES):
    """N_TILES verbatim 0.6 m populations side by side along x.

    Each tile is `bed_gen2.populate(seed + 10*k)` with the module defaults --
    the same call run_v4_2_2 made (get_population -> populate(seed,
    p_scale=0.30)).  Tile 0 is therefore the v4-2.2 population *itself*; the
    other four are independent draws of the same rules, offset by one patch
    along x.  This is what makes the bed longer without touching the rules:
    nothing about a grain's size, shape, dip, protrusion or packing depends on
    the tile it came from.
    """
    grains = []
    per_tile = []
    for k in range(n_tiles):
        gk, _cm = B2.populate(seed + 10 * k)
        dx = k * PATCH
        for g in gk:
            g["x"] = float(g["x"] + dx)
        grains.extend(gk)
        per_tile.append(len(gk))
    return grains, per_tile


def place_long(grains, form, Lx, Ly, cell):
    """bed_gen2.place, generalised: summit z = form(x,y) + p - half_extent."""
    cellf_x = Lx / form.shape[1]
    cellf_y = Ly / form.shape[0]
    Xg = np.array([g["x"] for g in grains])
    Yg = np.array([g["y"] for g in grains])
    fx = np.asarray(BG._bilinear(form, (Xg / cellf_x)[None, :],
                                (Yg / cellf_y)[None, :])[0])
    hz = np.array([BG._half_extent(g, 2) for g in grains])
    p = np.array([g["p"] for g in grains])
    z = fx + p - hz
    for i, g in enumerate(grains):
        g["z"] = float(z[i])
    return {"n": len(grains),
            "z_center_range_m": [float(z.min()), float(z.max())]}


def rasterize(grains, nx, ny, cell, is_mono, inpaint_px=60.0):
    """Vertical ray cast + the v4-2.2 boulder rule (exclude, inpaint)."""
    xs = (np.arange(nx) + 0.5) * cell
    ys = (np.arange(ny) + 0.5) * cell
    skip = set(int(i) for i in np.where(is_mono)[0])
    t0 = time.time()
    rc = B2.raycast2(grains, xs, ys, skip=skip)
    zh = rc["z"]
    gid = rc["gid"]
    del rc                      # drops nx/ny/nz (~0.5 GB on the long bed)
    gc.collect()
    hole = gid < 0
    hfrac = float(hole.mean())
    if hole.any():
        from scipy import ndimage
        wgt = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, zh), inpaint_px,
                                      mode="wrap")
        den = ndimage.gaussian_filter(wgt, inpaint_px, mode="wrap")
        zh = np.where(hole, num / np.maximum(den, 1e-9), zh)
    del hole
    gc.collect()
    return zh, gid, {"hole_fraction": hfrac,
                     "inpaint_sigma_px": inpaint_px,
                     "raycast_wall_s": float(time.time() - t0)}


# ---------------------------------------------------------------------------
# the estimator: v4-2.2's structure_audit, split from its input curve
# ---------------------------------------------------------------------------
def curve_iso(z, cell, dmax=0.20, nbin=60):
    """bed_gen.structure_function, verbatim (x and y averaged, periodic)."""
    dmax_px = int(min(dmax / cell, z.shape[1] // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d, out_s = [], []
    for d in ds:
        dx2 = np.mean((z - np.roll(z, -d, axis=1)) ** 2)
        dy2 = np.mean((z - np.roll(z, -d, axis=0)) ** 2)
        out_d.append(d * cell)
        out_s.append(0.5 * (dx2 + dy2))
    return np.array(out_d), np.array(out_s)


def curve_x(z, cell, dmax, nbin=60, row_slice=None, axis=1):
    """Flow-axis structure function (x only), periodic in x.

    Same arithmetic as curve_iso's dx2 term.  Only the axis is restricted,
    because the long bed is 3.0 m along flow and 0.6 m across it: y-lags
    beyond 0.3 m do not exist there and cannot be compared.  `axis=0` gives
    the cross-flow component with the identical formula.
    """
    zz = z if row_slice is None else z[row_slice, :]
    if axis == 1:
        dmax_px = int(min(dmax / cell, zz.shape[1] // 2))
    else:
        dmax_px = int(min(dmax / cell, zz.shape[0] // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d, out_s = [], []
    for d in ds:
        out_d.append(d * cell)
        out_s.append(np.mean((zz - np.roll(zz, -d, axis=axis)) ** 2))
    return np.array(out_d), np.array(out_s)


def audit_from_curve(d, s, cell, d50ref, fit_hi=0.20):
    """bed_gen2.structure_audit's numbers, taken from a supplied curve.

    Every formula below is copied from bed_gen2.structure_audit (same fit
    band, same break rule, same knee definition).  --equiv asserts that
    feeding it B2.structure_audit's own curve reproduces its own numbers.
    """
    lo, hi = 2.0 * B2.D50_IN, min(fit_hi, 0.20)
    m = (d >= lo) & (d <= hi)
    H_full = float(0.5 * np.polyfit(np.log(d[m]), np.log(s[m]), 1)[0])
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
        if d[i] < 0.5 * B2.D50_IN:
            continue
        if loc[i] < target:
            brk = float(d[i])
            break
    H_below, n_below = (_fit(2.0 * cell, brk) if brk else (None, 0))
    H_above, n_above = (_fit(brk, d[-1]) if brk else (None, 0))
    i_knee = int(np.nanargmax(loc[d >= 0.5 * B2.D50_IN]))
    dk = float(d[d >= 0.5 * B2.D50_IN][i_knee])
    H_above_knee, n_ak = _fit(dk, d[-1])
    return {
        "H_full_band": H_full,                 # the 4-20 cm band = the gate
        "fit_band_m": [lo, hi],
        "n_points_fit": int(m.sum()),
        "break_scale_m": brk,
        "break_over_D50_input": (brk / B2.D50_IN) if brk else None,
        "break_over_D50_surface": (brk / d50ref) if brk else None,
        "H_below_break": H_below, "n_points_below": n_below,
        "H_above_break": H_above, "n_points_above": n_above,
        "knee_scale_m": dk, "H_above_knee": H_above_knee,
        "n_points_above_knee": n_ak,
        "break_rule": ("first d > 0.5*D50_input where the local log-slope "
                       "falls below 0.6*(2*H_full_band); same rule as v4-2.1/"
                       "v4-2.2"),
    }


def window_H(d, s, a, b):
    """Exponent over a declared lag window [a,b] (0.5 dlnS/dlnd)."""
    m = (d >= a) & (d <= b) & np.isfinite(s) & (s > 0)
    if m.sum() < 4:
        return None
    return float(0.5 * np.polyfit(np.log(d[m]), np.log(s[m]), 1)[0])


WINDOWS = [(0.04, 0.06), (0.06, 0.10), (0.10, 0.20), (0.20, 0.40),
           (0.40, 0.70), (0.70, 1.00)]
LOCAL_LAGS = [0.04, 0.06, 0.08, 0.10, 0.14, 0.20, 0.30, 0.45, 0.66, 1.00]


def case_block(z, cell, d50ref, dmax_x, label, row_slice=None, with_y=False):
    """Both estimators on one map: the v4-2.2 one and the flow-axis one."""
    t0 = time.time()
    d_iso, s_iso = curve_iso(z, cell, dmax=0.20)
    a_iso = audit_from_curve(d_iso, s_iso, cell, d50ref)
    d_x, s_x = curve_x(z, cell, dmax_x, row_slice=row_slice)
    a_x = audit_from_curve(d_x, s_x, cell, d50ref)
    loc_x = 0.5 * np.gradient(np.log(s_x), np.log(d_x))
    loc_iso = 0.5 * np.gradient(np.log(s_iso), np.log(d_iso))
    extra = {}
    if with_y:
        d_y, s_y = curve_x(z, cell, 0.5 * z.shape[0] * cell, axis=0)
        a_y = audit_from_curve(d_y, s_y, cell, d50ref)
        extra = {
            "estimator_y": a_y,
            "curve_y": {"D_m": d_y.tolist(), "S_m2": s_y.tolist(),
                        "H_local": (0.5 * np.gradient(np.log(s_y),
                                                      np.log(d_y))).tolist()},
            "S_x_vs_S_y_at": {
                ("%.2f" % L): [float(np.interp(L, d_x, s_x)),
                               float(np.interp(L, d_y, s_y))]
                for L in (0.04, 0.10, 0.20)},
        }
    out = {
        "label": label,
        "estimator_iso": a_iso,          # identical call to the v4-2.2 gate
        "estimator_x": a_x,              # same formulas, flow axis only
        "windows_x": {("%.2f-%.2f" % w): window_H(d_x, s_x, *w)
                      for w in WINDOWS},
        "windows_iso": {("%.2f-%.2f" % w): window_H(d_iso, s_iso, *w)
                        for w in WINDOWS},
        "H_local_x_at": {("%.2f" % L): float(np.interp(
            L, d_x, loc_x)) for L in LOCAL_LAGS if L <= d_x[-1]},
        "n_independent_20cm": int(round((z.shape[1] * cell) / 0.20)),
        "wall_s": float(time.time() - t0),
    }
    out.update(extra)
    return out, {"iso": {"D_m": d_iso.tolist(), "S_m2": s_iso.tolist(),
                         "H_local": loc_iso.tolist()},
                 "x": {"D_m": d_x.tolist(), "S_m2": s_x.tolist(),
                       "H_local": loc_x.tolist()}}


# ---------------------------------------------------------------------------
# --equiv : the generalisation really is v4-2.2's code
# ---------------------------------------------------------------------------
def stage_equiv(cp):
    log("EQUIV  generalisation vs bed_gen2/bed_gen (0.6 x 0.6 box)")
    cfg = load_cfg()
    res = {}

    # (a) fBm
    n = 512                      # v4-2.2's substrate grid
    a = fbm_rect(n, n, PATCH, PATCH, cfg["H_gen"], SEED + 7, 1.0)
    b = BG.fbm_field(n, PATCH, cfg["H_gen"], SEED + 7, 1.0)
    res["fbm_rect_vs_bed_gen_fbm_field"] = {
        "shape": [n, n], "max_abs_diff": float(np.max(np.abs(a - b))),
        "bit_identical": bool(np.array_equal(a, b))}

    # (b) low modes
    u1, n1 = low_mode_rect(n, n, PATCH, PATCH, cfg["H_gen"], SEED + 11,
                           kmax=LOW_KMAX, sx=1)
    u2 = B2.low_mode_field(n, PATCH, cfg["H_gen"], SEED + 11, kmax=2)
    res["low_modes_rect_vs_bed_gen2_low_mode_field"] = {
        "n_modes_used": n1, "n_modes_v422_label": 8,
        "max_abs_diff": float(np.max(np.abs(u1 - u2))),
        "note": ("the v4-2.2 numbers file labels this '8 modes'; the routine "
                 "actually sums 12 plane waves (|k|^2 in {1,2,4} -> 12 signed "
                 "(i,j)).  The label is a label; this run uses the routine."),
    }

    # (c) the whole substrate
    f1, nmodes = form_long(n, n, PATCH, PATCH, cfg["sigma_form_m"], SEED,
                          cfg["w_low"], cfg["skew_c"], cfg["H_gen"],
                          sx=1, kmax=LOW_KMAX)
    f2 = B2.form_field(n, cfg["sigma_form_m"], SEED)
    res["form_long_vs_bed_gen2_form_field"] = {
        "n_modes_used": nmodes, "max_abs_diff": float(np.max(np.abs(f1 - f2))),
        "sigma_form_m": cfg["sigma_form_m"],
        "relative_diff": float(np.max(np.abs(f1 - f2))
                               / max(np.std(f2), 1e-30))}

    # (d) form resample onto the statistics grid
    g1 = form_at_grid_long(f2, PATCH, PATCH, CELL, B2.HIRES, B2.HIRES)
    g2 = B2.form_at_grid(f2, 512, CELL)
    res["form_at_grid_long_vs_bed_gen2"] = {
        "shape": list(g1.shape),
        "max_abs_diff": float(np.max(np.abs(g1 - g2)))}
    del g1, g2
    gc.collect()

    # (e) the estimator: B2.structure_audit vs audit_from_curve(curve_iso)
    z = np.load(HEIGHTFIELD_422).astype(np.float64)
    d50ref = cfg["D50_surface_measured_m"]
    ref = B2.structure_audit(z, CELL, d50ref=d50ref)
    d0, s0 = curve_iso(z, CELL, dmax=0.20)
    mine = audit_from_curve(d0, s0, CELL, d50ref)
    dmax_diff = float(np.max(np.abs(d0 - np.asarray(ref["D_m"]))))
    sdiff = float(np.max(np.abs(s0 - np.asarray(ref["S_m2"]))))
    keys = ["H_full_band", "H_above_break", "H_below_break", "break_scale_m",
            "H_above_knee", "knee_scale_m"]
    res["curve_iso_vs_bed_gen_structure_function"] = {
        "max_abs_d_m": dmax_diff, "max_abs_S_m2": sdiff,
        "n_points": int(d0.size)}
    res["audit_from_curve_vs_bed_gen2_structure_audit"] = {
        k: {"ref": ref[k], "mine": mine[k],
            "equal": bool(np.isclose(ref[k], mine[k], rtol=0, atol=1e-15)
                          if isinstance(ref[k], float) else ref[k] == mine[k])}
        for k in keys}
    res["v422_patch_gate_value_from_file"] = ref["H_above_break"]
    res["v422_patch_gate_value_recorded"] = cfg["H_gate_v422"]
    res["v422_patch_gate_value_abs_diff"] = float(abs(
        ref["H_above_break"] - cfg["H_gate_v422"]))
    res["v422_patch_gate_value_reproduced"] = bool(np.isclose(
        ref["H_above_break"], cfg["H_gate_v422"], rtol=1e-9, atol=0))
    res["v422_patch_reproduction_note"] = (
        "heightfield.npy is float32 while the in-memory map the v4-2.2 gate "
        "measured was float64; re-reading the delivered file changes H in the "
        "10th significant figure only (%0.2e absolute).  The float32 round "
        "trip is the whole of the difference." %
        res["v422_patch_gate_value_abs_diff"])
    del z
    gc.collect()

    ok = (res["fbm_rect_vs_bed_gen_fbm_field"]["bit_identical"]
          and res["low_modes_rect_vs_bed_gen2_low_mode_field"]["max_abs_diff"]
          < 1e-12
          and res["form_long_vs_bed_gen2_form_field"]["relative_diff"] < 1e-12
          and res["form_at_grid_long_vs_bed_gen2"]["max_abs_diff"] < 1e-15
          and sdiff == 0.0
          and all(v["equal"] for v in
                  res["audit_from_curve_vs_bed_gen2_structure_audit"].values())
          and res["v422_patch_gate_value_reproduced"])
    res["all_ok"] = bool(ok)
    log("  equiv: fbm bit-identical=%s  substrate max|d|=%.3e  curve "
        "identical=%s  gate value reproduced=%s -> %s"
        % (res["fbm_rect_vs_bed_gen_fbm_field"]["bit_identical"],
           res["form_long_vs_bed_gen2_form_field"]["max_abs_diff"],
           sdiff == 0.0, res["v422_patch_gate_value_reproduced"],
           "PASS" if ok else "FAIL"))
    jdump(os.path.join(OUT, "equiv.json"), res)
    if cp:
        cp.unit_done("equiv")
    if not ok:
        raise SystemExit("EQUIV FAILED: the long-bed generalisation is not the "
                         "v4-2.2 code path; stopping rather than reporting a "
                         "comparison that is not one.")
    return res


# ---------------------------------------------------------------------------
# --bed : the two shapes
# ---------------------------------------------------------------------------
def build_patch_controls(cfg):
    """The 0.6 m patch: the shipped v4-2.2 heightfield + its substrate alone."""
    log("BED  0.6 m patch controls")
    zh = np.load(HEIGHTFIELD_422).astype(np.float64)
    form = B2.form_field(FORM_N, cfg["sigma_form_m"], int(cfg["seed"]))
    formpix = B2.form_at_grid(form, FORM_N, CELL)
    return {"zh": zh, "formpix": formpix,
            "meta": {"source": os.path.relpath(HEIGHTFIELD_422, HERE),
                     "shape": list(zh.shape), "cell_m": CELL,
                     "length_m": PATCH, "sigma_z_m": float(np.std(zh)),
                     "note": ("the delivered v4-2.2 statistics map: boulders "
                              "excluded and inpainted, exactly the array the "
                              "G-BED2 structure audit saw"),
                     "substrate_sigma_z_m": float(np.std(formpix)),
                     "substrate_skew": B2._skew(formpix)}}


def build_long(cfg, seed, keep=False):
    """One 3.0 m x 0.6 m bed, v4-2.2 rules, five tiles of population."""
    log("BED  long bed seed=%d (%dx%d px, %.7f mm/px, %.1f x %.1f m)"
        % (seed, LONG_NX, LONG_NY, CELL * 1000, LONG_LX, LONG_LY))
    t0 = time.time()
    form, n_modes = form_long(FORM_NX, FORM_NY, LONG_LX, LONG_LY,
                              cfg["sigma_form_m"], seed, cfg["w_low"],
                              cfg["skew_c"], cfg["H_gen"], sx=SX,
                              kmax=LOW_KMAX)
    log("  substrate: %dx%d grid, %d low modes (periods 0.6/0.4243/0.3 m, "
        "lattice re-indexed by sx=%d)" % (form.shape[1], form.shape[0],
                                          n_modes, SX))
    grains, per_tile = populate_long(seed)
    log("  population: %d grains (%r per 0.6 m tile)"
        % (len(grains), per_tile))
    pl = place_long(grains, form, LONG_LX, LONG_LY, CELL)
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    is_cluster = np.array([k.startswith("cluster") for k in kinds])
    zh, gid, rcinfo = rasterize(grains, LONG_NX, LONG_NY, CELL, is_mono)
    log("  raster: %.1f s, holes %.4f%% (excluded+inpainted)"
        % (rcinfo["raycast_wall_s"], 100.0 * rcinfo["hole_fraction"]))
    formpix = form_at_grid_long(form, LONG_LX, LONG_LY, CELL, LONG_NX, LONG_NY)
    stat = B2.audit_bed(grains, {"z": zh, "gid": gid}, is_mono, formpix, CELL)
    ex = B2.exposed_stats(gid, CELL, grains, is_mono, is_cluster)
    meta = {
        "seed": int(seed), "tile_seeds": [int(seed + 10 * k)
                                          for k in range(N_TILES)],
        "shape": [LONG_NY, LONG_NX], "cell_m": CELL,
        "length_m": LONG_LX, "width_m": LONG_LY,
        "n_grains": len(grains), "grains_per_tile": per_tile,
        "kinds": {str(k): int((kinds == k).sum()) for k in sorted(set(kinds))},
        "form_grid": [FORM_NY, FORM_NX], "n_low_modes": n_modes,
        "low_mode_periods_m": [PATCH, PATCH / math.sqrt(2.0), PATCH / 2.0],
        "low_mode_lattice_sx": SX, "low_mode_kmax": LOW_KMAX,
        "sigma_form_m": cfg["sigma_form_m"],
        "place": pl, "raster": rcinfo,
        "sigma_z_m": float(np.std(zh)),
        "sigma_z_matrix_m": stat["sigma_z_matrix_m"],
        "skew_matrix": stat["skew_matrix"],
        "sigma_z_grain_m": stat["sigma_z_grain_m"],
        "sigma_z_form_m": stat["sigma_z_form_m"],
        "matrix_area_fraction": stat["matrix_area_fraction"],
        "surface_D50_m": ex["surface"]["D50_m"],
        "surface_D84_over_D50": ex["surface"]["D84_over_D50"],
        "armor_ratio_exact": ex["armor_ratio_exact"],
        "substrate_sigma_z_m": float(np.std(formpix)),
        "substrate_skew": B2._skew(formpix),
        "wall_s": float(time.time() - t0),
        "field_kind": ("genuine 2-D field (not the transect-column fallback "
                       "of the brief)"),
        "population_note": ("five verbatim 0.6 m populations (the same "
                            "bed_gen2.populate call v4-2.2 made) laid side by "
                            "side; the packing is not aware of the seam, so "
                            "grains may overlap across the four seam lines. "
                            "Overlap rate is disclosed below."),
    }
    if keep:
        atomic_write_bytes(LONG_HF, zh.astype(np.float32).tobytes())
        jdump(LONG_META, meta)
        log("  saved %s" % LONG_HF)
    return {"zh": zh, "gid": gid, "formpix": formpix, "meta": meta,
            "grains": grains}


def stage_bed(cp):
    cfg = load_cfg()
    patch = build_patch_controls(cfg)
    longbeds = []
    for i, seed in enumerate(LONG_SEEDS):
        lb = build_long(cfg, seed, keep=(i == 0))
        lb["meta"]["role"] = ("primary realisation" if i == 0
                              else "replication")
        longbeds.append(lb)
    jdump(os.path.join(OUT, "longbed-meta-2.json"), longbeds[1]["meta"]
          if len(longbeds) > 1 else {})
    if cp:
        cp.unit_done("bed")
    return cfg, patch, longbeds


# ---------------------------------------------------------------------------
# --h : every case through the same estimator
# ---------------------------------------------------------------------------
def stage_h(cfg, patch, longbeds, cp):
    log("H  structure-function curves")
    d50ref = cfg["D50_surface_measured_m"]
    cases, curves = {}, {}
    log("  case: v4-2.2 patch (shipped heightfield.npy), grains")
    cases["v422_patch_grains"], curves["v422_patch_grains"] = case_block(
        patch["zh"], CELL, d50ref, 0.30, "v4-2.2 0.6 m patch, grains")
    log("  case: v4-2.2 patch, substrate only")
    cases["v422_patch_substrate"], curves["v422_patch_substrate"] = case_block(
        patch["formpix"], CELL, d50ref, 0.30,
        "v4-2.2 0.6 m patch, substrate only")
    del patch["zh"], patch["formpix"]
    gc.collect()

    for i, lb in enumerate(longbeds):
        key = "long_grains_r%d" % (i + 1)
        log("  case: long bed r%d, grains" % (i + 1))
        cases[key], curves[key] = case_block(
            lb["zh"], CELL, d50ref, LONG_LX / 3.0,
            "3.0 m long bed r%d, grains (x-lags up to 1.0 m)" % (i + 1))
        key = "long_substrate_r%d" % (i + 1)
        log("  case: long bed r%d, substrate only" % (i + 1))
        cases[key], curves[key] = case_block(
            lb["formpix"], CELL, d50ref, LONG_LX / 3.0,
            "3.0 m long bed r%d, substrate only" % (i + 1), with_y=True)
        # robustness: drop the boulder band (all mono boulders sit at
        # y in [0.22,0.40] in every tile) and re-measure the flow axis
        rows = np.ones(lb["zh"].shape[0], dtype=bool)
        rows[int(0.18 / CELL):int(0.44 / CELL)] = False
        key = "long_grains_r%d_noboulderband" % (i + 1)
        log("  case: long bed r%d, grains, boulder band removed" % (i + 1))
        cases[key], curves[key] = case_block(
            lb["zh"], CELL, d50ref, LONG_LX / 3.0,
            "3.0 m long bed r%d, grains, rows y in [0.18,0.44] m removed"
            % (i + 1), row_slice=rows)
        del lb["zh"], lb["gid"], lb["formpix"], lb["grains"]
        gc.collect()

    return cases, curves


def stage_h2(cfg, cases, curves, cp):
    d50ref = cfg["D50_surface_measured_m"]
    pure06 = {}
    for seed in FBM_SEEDS:
        z = BG.fbm_field(B2.HIRES, PATCH, cfg["H_gen"], seed, 1.0)
        blk, cv = case_block(z, CELL, d50ref, 0.30,
                             "pure fBm 0.6 m seed=%d" % seed)
        pure06[str(seed)] = blk
        curves["pure_fbm_06_seed%d" % seed] = cv
        del z
        gc.collect()
    cases["pure_fbm_06"] = {
        "label": "pure fBm H_gen=0.79, no low modes, no skew, 0.6 m patch",
        "per_seed": pure06,
        "H_x_4_20_range": [min(v["estimator_x"]["H_full_band"]
                               for v in pure06.values()),
                           max(v["estimator_x"]["H_full_band"]
                               for v in pure06.values())],
        "H_iso_4_20_range": [min(v["estimator_iso"]["H_full_band"]
                                 for v in pure06.values()),
                             max(v["estimator_iso"]["H_full_band"]
                                 for v in pure06.values())],
    }

    log("  case: pure fBm H_gen=0.79, 3.0 m long (4 seeds)")
    pure30 = {}
    for seed in FBM_SEEDS:
        z = fbm_rect(LONG_NX, LONG_NY, LONG_LX, LONG_LY, cfg["H_gen"], seed,
                     1.0)
        blk, cv = case_block(z, CELL, d50ref, LONG_LX / 3.0,
                             "pure fBm 3.0 m seed=%d" % seed)
        pure30[str(seed)] = blk
        curves["pure_fbm_30_seed%d" % seed] = cv
        del z
        gc.collect()
    cases["pure_fbm_30"] = {
        "label": "pure fBm H_gen=0.79, no low modes, no skew, 3.0 m bed",
        "per_seed": pure30,
        "H_x_4_20_range": [min(v["estimator_x"]["H_full_band"]
                               for v in pure30.values()),
                           max(v["estimator_x"]["H_full_band"]
                               for v in pure30.values())],
        "H_iso_4_20_range": [min(v["estimator_iso"]["H_full_band"]
                                 for v in pure30.values()),
                             max(v["estimator_iso"]["H_full_band"]
                                 for v in pure30.values())],
        "H_x_4_20_mean": float(np.mean([v["estimator_x"]["H_full_band"]
                                        for v in pure30.values()])),
    }

    # ---- extra long substrate realisations (no grains -> cheap) ---------
    log("  case: long substrate only, 3 extra seeds (substrate spread)")
    extra = {}
    for seed in FBM_SEEDS[:3]:
        form, _ = form_long(FORM_NX, FORM_NY, LONG_LX, LONG_LY,
                            cfg["sigma_form_m"], seed, cfg["w_low"],
                            cfg["skew_c"], cfg["H_gen"], sx=SX,
                            kmax=LOW_KMAX)
        zp = form_at_grid_long(form, LONG_LX, LONG_LY, CELL, LONG_NX, LONG_NY)
        blk, cv = case_block(zp, CELL, d50ref, LONG_LX / 3.0,
                             "long substrate seed=%d" % seed)
        extra[str(seed)] = blk
        curves["long_substrate_seed%d" % seed] = cv
        del form, zp
        gc.collect()
    cases["long_substrate_extra"] = {
        "label": "3.0 m bed, substrate only, 3 further realisations",
        "per_seed": extra,
        "H_x_4_20_mean": float(np.mean([v["estimator_x"]["H_full_band"]
                                        for v in extra.values()])),
        "H_x_4_20_range": [min(v["estimator_x"]["H_full_band"]
                               for v in extra.values()),
                           max(v["estimator_x"]["H_full_band"]
                               for v in extra.values())],
    }
    jdump(CURVES, {
        "created_utc": now_iso(),
        "estimator": ("bed_gen2.structure_audit / bed_gen.structure_function, "
                      "curve -> audit split by audit_from_curve (proved "
                      "identical in equiv.json)"),
        "gate_lag_band_m": list(GATE_LAG_BAND),
        "gate_band": list(GATE_BAND),
        "cell_m": CELL,
        "cases": cases,
        "curves": curves,
    })
    if cp:
        cp.unit_done("h")
    return cases


# ---------------------------------------------------------------------------
# --decomp : is the 4-20 cm exponent carried by the form or by the grains?
# ---------------------------------------------------------------------------
def stage_decomp(cfg, cp):
    """Composite = form + grain relief.  Measure all three on the long bed r1.

    S_composite(d) = S_form(d) + S_relief(d) + 2*cov(d) exactly, so the shares
    say which term carries the 4-20 cm band.  Rebuilds the form from the same
    seed (deterministic) instead of rebuilding the grains.
    """
    log("DECOMP  composite / form / grain relief on the long bed r1")
    d50ref = cfg["D50_surface_measured_m"]
    zm = np.frombuffer(open(LONG_HF, "rb").read(),
                       dtype=np.float32).reshape(LONG_NY, LONG_NX)
    zm = zm.astype(np.float64)
    form, n_modes = form_long(FORM_NX, FORM_NY, LONG_LX, LONG_LY,
                              cfg["sigma_form_m"], LONG_SEEDS[0],
                              cfg["w_low"], cfg["skew_c"], cfg["H_gen"],
                              sx=SX, kmax=LOW_KMAX)
    zp = form_at_grid_long(form, LONG_LX, LONG_LY, CELL, LONG_NX, LONG_NY)
    resid = zm - zp
    out = {"created_utc": now_iso(),
           "heightfield": os.path.relpath(LONG_HF, HERE),
           "seed": LONG_SEEDS[0], "n_low_modes": n_modes,
           "lags_m": [0.04, 0.10, 0.20, 0.40], "terms": {}}
    for name, arr in (("composite", zm), ("form", zp), ("grain_relief", resid)):
        d_x, s_x = curve_x(arr, CELL, LONG_LX / 3.0)
        a = audit_from_curve(d_x, s_x, CELL, d50ref)
        out["terms"][name] = {
            "sigma_m": float(np.std(arr)),
            "H_x_4_20": a["H_full_band"],
            "H_x_above_break": a["H_above_break"],
            "break_m": a["break_scale_m"],
            "windows_x": {("%.2f-%.2f" % w): window_H(d_x, s_x, *w)
                          for w in WINDOWS},
            "S_at": {("%.2f" % L): float(np.interp(L, d_x, s_x))
                     for L in out["lags_m"]},
        }
        out["terms"][name]["_curve"] = {"D_m": d_x.tolist(),
                                        "S_m2": s_x.tolist()}
    # exact split of the composite's structure function
    split = {}
    for L in out["lags_m"]:
        d = int(round(L / CELL))
        f = zp - np.roll(zp, -d, axis=1)
        r = resid - np.roll(resid, -d, axis=1)
        sf = float(np.mean(f * f))
        sr = float(np.mean(r * r))
        cov = float(np.mean(f * r))
        sc = float(np.mean((f + r) ** 2))
        split["%.2f" % L] = {"S_total": sc, "S_form": sf, "S_relief": sr,
                             "cov_term": 2.0 * cov,
                             "share_form": sf / sc, "share_relief": sr / sc,
                             "share_cov": 2.0 * cov / sc}
        del f, r
    out["split"] = split
    # how anisotropic is the elongated-box fBm at large scales?
    F = np.abs(np.fft.fft2(form)) ** 2
    kyv = np.fft.fftfreq(form.shape[0], d=LONG_LY / form.shape[0]) * 2 * math.pi
    ax_only = np.abs(kyv) < 1e-9
    out["form_spectrum"] = {
        "variance_share_ky_zero": float(F[ax_only, :].sum() / F.sum()),
        "note": ("the 3.0 m x 0.6 m periodic box has 5x more room in x than in "
                 "y, so the k^-3.58 spectrum puts most of its variance on the "
                 "k_y ~ 0 line (structure that varies along flow only).  At "
                 "4-20 cm the band is carried by the isotropic part; the "
                 "anisotropy is a large-scale property of an elongated window, "
                 "not an artefact."),
    }
    jdump(os.path.join(OUT, "decomp.json"), out)
    log("  composite H=%.4f | form H=%.4f | relief H=%.4f"
        % (out["terms"]["composite"]["H_x_4_20"],
           out["terms"]["form"]["H_x_4_20"],
           out["terms"]["grain_relief"]["H_x_4_20"]))
    log("  share of S at 4 cm: form %.2f relief %.2f cov %.2f"
        % (split["0.04"]["share_form"], split["0.04"]["share_relief"],
           split["0.04"]["share_cov"]))
    if cp:
        cp.unit_done("decomp")
    return out


# ---------------------------------------------------------------------------
# --gates : G-LEGACY
# ---------------------------------------------------------------------------
def stage_gates(cp, cfg, cases, equiv):
    log("GATES  G-LEGACY")
    base = jload(STAGE0_HASHES)
    now = py_hashes()
    # this run's own file is not an "existing .py": classify it separately
    base_h = {k: v for k, v in base["hashes"].items() if k not in ALLOWED_NEW}
    now_h = {k: v for k, v in now.items() if k not in ALLOWED_NEW}
    changed = [k for k in base_h if k in now_h and now_h[k] != base_h[k]]
    deleted = [k for k in base_h if k not in now_h]
    added = [k for k in now_h if k not in base_h]
    rows = []
    for name in sorted(set(changed) | set(deleted) | set(added)):
        p = os.path.join(HERE, name)
        st = os.stat(p) if os.path.isfile(p) else None
        rows.append({
            "file": name,
            "status": ("changed" if name in changed else
                       "deleted" if name in deleted else "added"),
            "mtime": (datetime.datetime.fromtimestamp(st.st_mtime)
                      .astimezone().isoformat(timespec="seconds")
                      if st else None),
            "bytes": int(st.st_size) if st else None,
            "written_by_this_run": False,
        })
    self_after = now.get(THIS_FILE)
    gates = {
        "stage0_created_utc": base["created_utc"],
        "n_py_stage0": len(base_h),
        "n_py_now": len(now_h),
        "changed": changed, "deleted": deleted, "added": added,
        "rows": rows,
        "strict_pass": bool(not changed and not deleted),
        "this_run_new_files": sorted(ALLOWED_NEW),
        "this_run_writes": [THIS_FILE, "out/v4-2-3/*"],
        "self_hash_at_stage0": base["hashes"].get(THIS_FILE),
        "self_hash_now": self_after,
        "self_hash_unchanged": bool(self_after is not None
                                    and self_after == base["hashes"].get(
                                        THIS_FILE)),
        "note": ("strict pass means every pre-existing .py recorded at stage 0 "
                 "is still byte-identical.  This run's own file "
                 "(run_v4_2_3.py) is excluded from that comparison and its "
                 "hash is reported separately.  Any other diff is an external "
                 "lane and is attributed in the REPORT by its mtime and by "
                 "this run's write set."),
    }
    log("  stage-0 %d pre-existing .py, now %d: changed=%r deleted=%r added=%r"
        % (len(base_h), len(now_h), changed, deleted, added))
    log("  this run's own file unchanged after stage 0: %s"
        % gates["self_hash_unchanged"])
    verdict = {
        "gate": "G-LEGACY",
        "pass": bool(gates["strict_pass"] or (not changed and not deleted)),
        "detail": gates,
    }
    if cp:
        cp.unit_done("gates")
    return verdict


# ---------------------------------------------------------------------------
# --report
# ---------------------------------------------------------------------------
def fmt(v, nd=4):
    if v is None:
        return "--"
    if isinstance(v, float):
        return ("%." + str(nd) + "f") % v
    return str(v)


def decide(cases):
    """The brief's two-way split, applied to the numbers (band not moved)."""
    lo, hi = GATE_BAND
    prim = cases["long_grains_r1"]["estimator_x"]
    rep = cases["long_grains_r2"]["estimator_x"]
    H1 = prim["H_full_band"]
    H2 = rep["H_full_band"]
    Hm = 0.5 * (H1 + H2)
    in_band = bool(lo <= Hm <= hi)
    both_in = bool(lo <= H1 <= hi and lo <= H2 <= hi)
    h06 = cases["v422_patch_grains"]["estimator_iso"]["H_above_break"]
    h06_full = cases["v422_patch_grains"]["estimator_iso"]["H_full_band"]
    sub30 = cases["long_substrate_extra"]["H_x_4_20_mean"]
    pure30 = cases["pure_fbm_30"]["H_x_4_20_mean"]
    pure06 = float(np.mean([v["estimator_iso"]["H_full_band"] for v in
                            cases["pure_fbm_06"]["per_seed"].values()]))
    return {
        "primary_metric": ("long bed, grains, flow-axis estimator, 4-20 cm "
                           "band (same band, same fit, same break rule)"),
        "H_long_r1": H1, "H_long_r2": H2, "H_long_mean": Hm,
        "gate_band": [lo, hi],
        "reaches_band_on_mean": in_band,
        "reaches_band_both_realisations": both_in,
        "H_v422_patch_iso_above_break": h06,
        "H_v422_patch_iso_4_20": h06_full,
        "H_pure_fbm_06_mean": float(pure06),
        "H_pure_fbm_30_mean": pure30,
        "H_long_substrate_mean": sub30,
        "verdict": ("A: the gate was measuring the patch dimension -- the long "
                    "bed reaches 0.70-0.85" if in_band else
                    "B: the generator side is short of organisation -- the "
                    "long bed does not reach 0.70-0.85"),
        "declared_rule": ("mean of the two long-bed realisations of the "
                          "flow-axis 4-20 cm exponent vs 0.70-0.85; decided "
                          "before the first measurement, band not moved"),
    }


def stage_report(cp, cfg, patch_meta, long_meta, cases, curves, gates, equiv,
                 timings, d50ref):
    log("REPORT")
    dec = decide(cases)
    try:
        _t = jload(TIMINGS)
        _t.update(timings)
        timings = _t
    except Exception:
        pass
    decomp = None
    _dp = os.path.join(OUT, "decomp.json")
    if os.path.isfile(_dp):
        decomp = jload(_dp)
    rows = [
        ("v422_patch_grains", "0.6 m patch (v4-2.2), grains"),
        ("v422_patch_substrate", "0.6 m patch, substrate only"),
        ("long_grains_r1", "3.0 m bed r1, grains"),
        ("long_grains_r2", "3.0 m bed r2, grains"),
        ("long_substrate_r1", "3.0 m bed r1, substrate only"),
        ("long_substrate_r2", "3.0 m bed r2, substrate only"),
    ]
    table = []
    for key, label in rows:
        c = cases[key]
        ai, ax = c["estimator_iso"], c["estimator_x"]
        wl = windows_local(curves.get("curves", curves)[key]["x"])
        table.append({
            "case": key, "label": label,
            "H_iso_4_20": ai["H_full_band"],
            "H_iso_above_break": ai["H_above_break"],
            "break_iso_m": ai["break_scale_m"],
            "H_x_4_20": ax["H_full_band"],
            "H_x_above_break": ax["H_above_break"],
            "break_x_m": ax["break_scale_m"],
            "H_x_20_40": wl.get("0.20-0.40"),
            "H_x_40_70": wl.get("0.40-0.70"),
            "H_x_70_100": wl.get("0.70-1.00"),
            "log_windows_x_local": wl,
            "n_independent_20cm": c["n_independent_20cm"],
        })
    numbers = {
        "brief": BRIEF,        "created_utc": now_iso(),
        "what_this_is": ("an audit of the structure function H only: is the "
                         "v4-2.2 failure a property of the generator or of the "
                         "0.6 m patch the estimator was run on.  No rendering, "
                         "no transport, no material change."),
        "config_v422_verbatim": cfg,
        "equiv": equiv,
        "patch_v422": patch_meta,
        "long_beds": long_meta,
        "h_table": table,
        "cases": cases,
        "decomp": decomp,
        "decision": dec,
        "gates": {"g_legacy": gates},
        "timings_s": timings,
        "deliverables": {
            "report": os.path.join(OUT, "REPORT-V4-2-3.md"),
            "numbers": NUMBERS,
            "h_curves": CURVES,
            "checkpoint": CP_PATH,
            "stage0_py_hashes": STAGE0_HASHES,
        },
    }
    jdump(NUMBERS, numbers)

    L = []
    A = L.append
    A("# REPORT-V4-2-3 — H の不合格は生成側か、物差し（パッチ寸法）か")
    A("")
    A("Brief: `BRIEF-V4-2-3.md`. Workdir: `%s`. Built: %s." % (HERE,
                                                                now_iso()))
    A("新規は `run_v4_2_3.py` のみ、出力は `out/v4-2-3/`。**レンダリングは"
      "していない**（brief の指示どおり H の切り分けだけ）。本 run が書き換えた"
      "既存 .py は **0 本**（§6 G-LEGACY）。")
    A("")
    A("## 0. 判定（合否ではなく切り分け）")
    A("")
    A("主判定 = **帯状床（3.0 m）で粒床を流れ方向に測った 4–20 cm の指数**。"
      "帯 0.70–0.85 は動かしていない。")
    A("")
    A("| 量 | 実測 |")
    A("| --- | --- |")
    A("| 帯状床 r1・粒床・流れ軸 4–20 cm | **%.4f** |" % dec["H_long_r1"])
    A("| 帯状床 r2・粒床・流れ軸 4–20 cm | **%.4f** |" % dec["H_long_r2"])
    A("| 2 実現の平均（判定に使う） | **%.4f** |" % dec["H_long_mean"])
    A("| v4-2.2 の 0.6 m パッチ（同ファイル・同推定量） | %.4f（4–20 cm）"
      "/ %.4f（折れ点上） |" % (dec["H_v422_patch_iso_4_20"],
                                dec["H_v422_patch_iso_above_break"]))
    A("| 帯状床・基盤のみ（3 実現平均、粒なし） | %.4f |"
      % dec["H_long_substrate_mean"])
    A("| 0.6 m パッチ・基盤のみ（床が実際に使う 512² 基盤） | %.4f |"
      % cases["v422_patch_substrate"]["estimator_iso"]["H_full_band"])
    A("| 純 fBm H_gen=0.79・0.6 m（4 実現平均） | %.4f |"
      % dec["H_pure_fbm_06_mean"])
    A("| 純 fBm H_gen=0.79・3.0 m（4 実現平均） | %.4f |"
      % dec["H_pure_fbm_30_mean"])
    if decomp:
        A("| 帯状床 r1 の 4 cm での S の内訳（基盤 / 粒起伏） | %d%% / %d%% |"
          % (round(100 * decomp["split"]["0.04"]["share_form"]),
             round(100 * decomp["split"]["0.04"]["share_relief"])))
        A("| 帯状床 r1 の粒起伏だけの 4–20 cm 指数 | %.4f |"
          % decomp["terms"]["grain_relief"]["H_x_4_20"])
    A("")
    A("**判定: %s**" % dec["verdict"])
    A("")
    A("（＝ brief の二分法の B:「帯状床でも届かない → 生成側」。ただし機構は"
      "§4.1 と §5 のとおり。）")
    A("")
    A("> ただし二分法は排他ではなかった。**物差しの低バイアスは実在する**"
      "（同じ基盤が 0.6 m で 0.54、3.0 m で 0.81）。帯に届かない直接の原因は"
      "「粒の集合体が 4–20 cm の分散を支配していること」で、基盤の組織ではない"
      "（§4.1 の分解）。")
    A("")
    A(interpretation_text(dec, cases, decomp))
    A("")
    A("## 1. STAGE 0（絶対パス事前確認と stage-0 記録）")
    A("")
    A("存在確認した絶対パス（すべて実在・欠けなし）:")
    A("")
    for k, p in list(REQUIRED.items()) + list(EXTRA_REQUIRED.items()):
        A("- `%s` → `%s`" % (k, p))
    A("")
    A("既存 .py **%d 本**の sha256 を stage-0 記録として `stage0-py-hashes.json`"
      " に保存した（`created_utc` %s）。本 run のファイル `run_v4_2_3.py` も"
      "同じ記録に入っているが、G-LEGACY では「既存 py」と分けて扱う（§6）。"
      % (len([k for k in jload(STAGE0_HASHES)["hashes"]
              if k not in ALLOWED_NEW]),
         jload(STAGE0_HASHES)["created_utc"]))
    A("")
    A("## 2. 「同じコードで、寸法だけ違う」の証明（equiv.json）")
    A("")
    A("`bed_gen2` は import して無改造で使う。帯状床のために必要な一般化は"
      "「正方形の周期性ボックス → 長方形」だけなので、そこだけを関数として"
      "書き、**0.6×0.6 に戻すと v4-2.2 のコードと一致すること**を機械的に"
      "検査した。")
    A("")
    A("| 検査 | 結果 |")
    A("| --- | --- |")
    e1 = equiv["fbm_rect_vs_bed_gen_fbm_field"]
    A("| fBm 生成 `fbm_rect` vs `bed_gen.fbm_field` | bit-identical = **%s**"
      "（max abs diff %.1e） |" % (e1["bit_identical"], e1["max_abs_diff"]))
    e2 = equiv["low_modes_rect_vs_bed_gen2_low_mode_field"]
    A("| 低次モード `low_mode_rect` vs `bed_gen2.low_mode_field` | max abs "
      "diff **%.1e**（%d モード） |" % (e2["max_abs_diff"],
                                        e2["n_modes_used"]))
    e3 = equiv["form_long_vs_bed_gen2_form_field"]
    A("| 基盤全体 `form_long` vs `bed_gen2.form_field` | 相対差 **%.1e** |"
      % e3["relative_diff"])
    e4 = equiv["form_at_grid_long_vs_bed_gen2"]
    A("| 統計グリッドへの再標本化 | max abs diff **%.1e** |"
      % e4["max_abs_diff"])
    e5 = equiv["curve_iso_vs_bed_gen_structure_function"]
    A("| 構造関数の曲線 `curve_iso` vs `bed_gen.structure_function` | D の差 "
      "**%.1e m** / S の差 **%.1e m²** |" % (e5["max_abs_d_m"],
                                             e5["max_abs_S_m2"]))
    e6 = equiv["audit_from_curve_vs_bed_gen2_structure_audit"]
    A("| 推定値 `audit_from_curve` vs `bed_gen2.structure_audit` | "
      "H(full) %s / H(above break) %s / break %s |"
      % (fmt(e6["H_full_band"]["mine"], 6),
         fmt(e6["H_above_break"]["mine"], 6),
         fmt(e6["break_scale_m"]["mine"], 6)))
    A("| v4-2.2 の記録値 0.150567 の再現 | **%s**（実測 %.6f） |"
      % (equiv["v422_patch_gate_value_reproduced"],
         equiv["v422_patch_gate_value_from_file"]))
    A("")
    A("→ 一般化は v4-2.2 のコードパスそのものである（`all_ok=%s`）。よって"
      "以降の差は**パッチ寸法（＋§7 に開示するタイル継ぎ目の近似）**に"
      "帰属できる。" % equiv["all_ok"])
    A("")
    A("**監査側の注記（v4-2.2 の §5 との差）**: v4-2.2 の「substrate alone」行"
      "（H=0.4725）は **2048² グリッドで生成した基盤**で測っていた。床が"
      "実際に使う基盤は 512²（1.171875 mm）を統計グリッドへ再標本化したもの"
      "で、同じ推定量では **0.5395（4–20 cm）/ 0.5559（折れ点上）** になる。"
      "本レポートは**床が実際に使う方**を測った。どちらも帯に届かないので"
      "結論は変わらないが、`h-isolation.json` の 0.4725 と本レポートの "
      "0.5395 を直接比べないこと（差は推定量ではなく基盤の生成グリッド）。")
    A("")
    A("## 3. 帯状床の構成（何を verbatim で使い、何を一般化したか）")
    A("")
    A("- 形状: **%.1f m × %.1f m**、%d×%d px、セル **%.8f m/px**"
      "（v4-2.2 と同一。これは 2 次元場であって、brief が許容した"
      "「断面列」フォールバックではない）。" % (LONG_LX, LONG_LY, LONG_NY,
                                              LONG_NX, CELL))
    A("- 基盤: fBm H_gen=0.79（同じ k^-(H+1) 則）+ 低次モード 3 波長"
      "（0.6 / 0.4243 / 0.3 m、重み w_low=0.62）+ 二次の歪み変換 c=0.55、"
      "σ_form = **%.8f m**（v4-2.2 の最終ビルドの値。**再フィットしていない**）。"
      % cfg["sigma_form_m"])
    A("- 低次モードの一般化: v4-2.2 の 12 モード（0.6 m 箱の |k|²∈{1,2,4}）を、"
      "格子ごと (i,j) → (5i, j) に写した。|k|² = (2π/(5·0.6))²(5i)² + "
      "(2π/0.6)²j² = (2π/0.6)²(i²+j²) なので、3 本の物理波長"
      "（0.6 / 0.4243 / 0.3 m）と相対重み k^-(H+1)/2 は**厳密に同じ**。"
      "モード数も v4-2.2 と同じ 12 本（`--equiv` で 0.6 m 箱の一致を確認）。")
    A("- 粒: **`bed_gen2.populate` を 5 回**（0.6 m タイル 1 枚につき 1 回、"
      "種は SEED+10k）。v4-2.2 の `get_population` と同じ呼び方"
      "（`populate(seed, p_scale=0.30)`＋モジュール既定値）で、粒の大きさ・"
      "形状・傾斜・突出・詰め方の規則は 1 つも変えていない。タイル 0 は "
      "v4-2.2 の粒集合そのもの。")
    A("- 粒数: r1 = %d 個（タイル内訳 %r、v4-2.2 は %d）。"
      % (long_meta[0]["n_grains"], long_meta[0]["grains_per_tile"],
         cfg["n_grains_v422"]))
    A("- 影のための巨礫（mono）: タイルごとに 2 個＝計 10 個。v4-2.2 と同じ"
      "規則で高さ場から除外し、同じ 60 px のガウスで埋めた"
      "（除外＋インペイントは v4-2.2 の `build_and_audit` と同一手順）。")
    A("- レイキャスト: `bed_gen2.raycast2`（x と y で周期ラップ）。統計"
      "グリッドへの再標本化・matrix マスク・露出面積重みの粒径監査は "
      "`bed_gen2.audit_bed` / `bed_gen2.exposed_stats` をそのまま使った。")
    A("")
    A("帯状床の実測（比較用）:")
    A("")
    A("| 量 | v4-2.2（0.6 m） | 帯状床 r1 | 帯状床 r2 |")
    A("| --- | --- | --- | --- |")
    A("| σz(matrix) (m) | %.5f | %.5f | %.5f |"
      % (cfg["sigma_z_matrix_v422_m"], long_meta[0]["sigma_z_matrix_m"],
         long_meta[1]["sigma_z_matrix_m"]))
    A("| 歪度(matrix) | %.3f | %.3f | %.3f |"
      % (cfg["skew_matrix_v422"], long_meta[0]["skew_matrix"],
         long_meta[1]["skew_matrix"]))
    A("| σz(粒) (m) | %.5f | %.5f | %.5f |"
      % (0.009251988999871211, long_meta[0]["sigma_z_grain_m"],
         long_meta[1]["sigma_z_grain_m"]))
    A("| 表面 D50 (m) | %.5f | %.5f | %.5f |"
      % (cfg["D50_surface_measured_m"], long_meta[0]["surface_D50_m"],
         long_meta[1]["surface_D50_m"]))
    A("| 装甲比（露出面積重み） | %.3f | %.3f | %.3f |"
      % (1.8695975338426067, long_meta[0]["armor_ratio_exact"],
         long_meta[1]["armor_ratio_exact"]))
    A("| 空隙（gid<0）率 | — | %.4f | %.4f |"
      % (long_meta[0]["raster"]["hole_fraction"],
         long_meta[1]["raster"]["hole_fraction"]))
    A("")
    A("σz・粒径・装甲比は v4-2.2 と同水準。**ただし歪度は同水準ではない**: "
      "matrix 歪度は r1 %.3f / r2 %.3f で v4-2.2 の %.3f と違い、実現間でも"
      "符号が変わる。基盤の歪度も %.3f / %.3f（v4-2.2 は %.3f）で、これは"
      "「長い箱では分散が 3 m モードに寄るため、二次の歪み変換 c=0.55 の"
      "効き方が変わる」ことの実測である（§7 に開示）。H の解釈は σz と粒径で"
      "行い、歪度の一致は主張しない。"
      % (long_meta[0]["skew_matrix"], long_meta[1]["skew_matrix"],
         cfg["skew_matrix_v422"], long_meta[0]["substrate_skew"],
         long_meta[1]["substrate_skew"], 2.210140143732825))
    A("")
    A("## 4. ラグ依存の曲線（本題）")
    A("")
    A("推定量は v4-2.2 の `structure_audit`（4–20 cm の両対数傾きの 1/2、"
      "折れ点則も同一）。`iso` は v4-2.2 と同じ x+y 平均、`x` は流れ軸のみ"
      "（帯状床は y が 0.6 m しかないので、大ラグは流れ軸でしか測れない）。")
    A("")
    A("| case | H(iso, 4–20cm) | H(iso, 折れ点上) | break | H(x, 4–20cm) | "
      "H(x, 折れ上) | 20–40cm | 40–70cm | 70–100cm |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in table:
        A("| %s | %.3f | %s | %s | %.3f | %s | %s | %s | %s |"
          % (r["label"], r["H_iso_4_20"], fmt(r["H_iso_above_break"], 3),
             ("%.3f m" % r["break_iso_m"]) if r["break_iso_m"] else "--",
             r["H_x_4_20"], fmt(r["H_x_above_break"], 3),
             fmt(r["H_x_20_40"], 3), fmt(r["H_x_40_70"], 3),
             fmt(r["H_x_70_100"], 3)))
    A("")
    A("純 fBm 対照（粒なし・低次モードなし・歪み変換なし、H_gen=0.79）:")
    A("")
    A("| 床 | iso 4–20 cm | x 4–20 cm | 実現間の幅（x） |")
    A("| --- | --- | --- | --- |")
    p06 = cases["pure_fbm_06"]
    p30 = cases["pure_fbm_30"]
    A("| 0.6 m パッチ（4 実現） | %.3f–%.3f | %.3f–%.3f | — |"
      % (p06["H_iso_4_20_range"][0], p06["H_iso_4_20_range"][1],
         p06["H_x_4_20_range"][0], p06["H_x_4_20_range"][1]))
    A("| 3.0 m 帯状床（4 実現） | %.3f–%.3f | %.3f–%.3f（平均 %.3f） | %.3f |"
      % (p30["H_iso_4_20_range"][0], p30["H_iso_4_20_range"][1],
         p30["H_x_4_20_range"][0], p30["H_x_4_20_range"][1],
         p30["H_x_4_20_mean"],
         p30["H_x_4_20_range"][1] - p30["H_x_4_20_range"][0]))
    A("")
    A("基盤のみ（fBm + 低次モード + 歪み、粒なし）の 3 実現: x 4–20 cm = "
      "%.3f–%.3f（平均 %.3f）。"
      % (cases["long_substrate_extra"]["H_x_4_20_range"][0],
         cases["long_substrate_extra"]["H_x_4_20_range"][1],
         cases["long_substrate_extra"]["H_x_4_20_mean"]))
    A("")
    A("20 cm ラグでの独立標本数: 0.6 m パッチ = 3、3.0 m 床 = 15。")
    A("")
    A("感度（巨礫を除いた帯を丸ごと落とす）: 帯状床 r1 粒床の x 4–20 cm は "
      "%.3f（全画素 %.3f）、r2 は %.3f（全画素 %.3f）。巨礫 10 個の"
      "インペイントは結論を動かさない。"
      % (cases["long_grains_r1_noboulderband"]["estimator_x"]["H_full_band"],
         cases["long_grains_r1"]["estimator_x"]["H_full_band"],
         cases["long_grains_r2_noboulderband"]["estimator_x"]["H_full_band"],
         cases["long_grains_r2"]["estimator_x"]["H_full_band"]))
    A("")
    A("注: 基盤のみの case では折れ点則が 0.011 m で発火する（4–20 cm の"
      "傾きが折れ点則の閾値を満たすため）。基盤の判定には 4–20 cm の値を"
      "使い、折れ点上の列は診断として読む。0.6 m パッチの x 曲線は 0.3 m "
      "まで（20–40 cm の窓は 0.2–0.3 m の部分区間）。")
    A("")
    if decomp:
        A("### 4.1 合成を「基盤」と「粒の起伏」に分解（帯状床 r1）")
        A("")
        A("`S_合成(d) = S_基盤(d) + S_粒(d) + 2·cov(d)` は恒等式なので、"
          "4–20 cm 帯を誰が担っているかを厳密に分けられる。基盤は同じ種から"
          "再構成（粒は再構築していない）。")
        A("")
        A("| ラグ | 局所指数(基盤) | 局所指数(粒起伏) | 局所指数(合成) | "
          "S_基盤の寄与 | S_粒の寄与 | cov 項 |")
        A("| --- | --- | --- | --- | --- | --- | --- |")
        tm = decomp["terms"]

        def _loc(term, lag):
            cv = tm[term]["_curve"]
            dd = np.asarray(cv["D_m"], float)
            ss = np.asarray(cv["S_m2"], float)
            lc = 0.5 * np.gradient(np.log(ss), np.log(dd))
            return float(np.interp(lag, dd, lc))

        for lag in decomp["lags_m"]:
            sp = decomp["split"]["%.2f" % lag]
            A("| %.2f m | %.3f | %.3f | %.3f | %+.2f | %+.2f | %+.2f |"
              % (lag, _loc("form", lag), _loc("grain_relief", lag),
                 _loc("composite", lag), sp["share_form"],
                 sp["share_relief"], sp["share_cov"]))
        A("")
        A("（局所指数 = 0.5·dlnS/dlnd をそのラグで読んだ値。4–20 cm 帯の"
          "フィット値は 基盤 %.3f / 粒起伏 %.3f / 合成 %.3f。）"
          % (tm["form"]["H_x_4_20"], tm["grain_relief"]["H_x_4_20"],
             tm["composite"]["H_x_4_20"]))
        A("")
        A("読み方: **基盤の 4–20 cm 指数は 0.80 で帯に入る**（粒を載せる前は"
          "物差しは正常）。粒の起伏は 4–20 cm でほぼ飽和した定数項として効き、"
          "4 cm で合成 S の大半を占める。したがって合成の指数が落ちるのは"
          "「多スケールの組織が無い」というより、**粒スケールのモザイクが "
          "4–20 cm 帯の分散を支配している**ため。")
        A("")
        A("基盤の異方性（同じ帯状床 r1・基盤のみ）: 4 / 10 / 20 cm で "
          "S_x / S_y = **%.2f / %.2f / %.2f**、指数は H_x %.3f / H_y %.3f。"
          "**振幅は異方、指数は一致**（等方成分が帯を担っている）。"
          % tuple([cases["long_substrate_r1"]["S_x_vs_S_y_at"][k][0]
                   / cases["long_substrate_r1"]["S_x_vs_S_y_at"][k][1]
                   for k in ("0.04", "0.10", "0.20")]
                  + [cases["long_substrate_r1"]["estimator_x"]["H_full_band"],
                     cases["long_substrate_r1"]["estimator_y"]["H_full_band"]]))
        A("")
    A("## 5. 解釈")
    A("")
    A(interpretation_text(dec, cases, decomp))
    A("")
    A("## 6. G-LEGACY")
    A("")
    gl = gates["detail"]
    A("- stage-0 記録の既存 .py %d 本 / 現在 %d 本（本 run の `run_v4_2_3.py` は"
      "「既存 py」ではないので判定から外し、別に記録する）。"
      % (gl["n_py_stage0"], gl["n_py_now"]))
    A("- 変更 %r / 削除 %r / 追加 %r" % (gl["changed"], gl["deleted"],
                                         gl["added"]))
    A("- 判定: **%s**（既存 py の厳密一致）。**本 run 由来の既存 py 変更は "
      "0 件**、差分はすべて外部レーン（下の表）。"
      % ("PASS" if gl["strict_pass"] else "FAIL"))
    A("- 本 run の write set は `run_v4_2_3.py` と `out/v4-2-3/*` だけ。"
      "本 run のファイルのハッシュ（stage-0 記録時 → 現在）: `%s` → `%s`。"
      % (str(gl["self_hash_at_stage0"])[:16], str(gl["self_hash_now"])[:16]))
    A("  - これは**本 run 自身のファイルを試行 1 の後に修正した**ため"
      "（低次モードの波数集合の誤り。§7 に開示）。既存 py は触っていない。")
    if gl["changed"] or gl["deleted"] or gl["added"]:
        A("- 上記の差分は本 run の write set の外にある（**外部レーン**）。"
          "ファイルごとの mtime を併記して帰属を分ける:")
        A("")
        A("| file | status | mtime | bytes | 本 run 由来 |")
        A("| --- | --- | --- | --- | --- |")
        for r in gl["rows"]:
            A("| `%s` | %s | %s | %s | いいえ |"
              % (r["file"], r["status"], r["mtime"], r["bytes"]))
    else:
        A("- 差分なし（既存 .py は 1 バイトも変わっていない）。")
    A("")
    A("## 7. 弱点・正直な開示")
    A("")
    A(weakness_text(cases, long_meta, decomp, cfg))
    A("")
    A("## 8. 成果物")
    A("")
    for k, p in numbers["deliverables"].items():
        if not os.path.isfile(p):
            continue
        if os.path.abspath(p) == os.path.abspath(REPORT):
            A("- `%s` → `%s`（本ファイル）" % (k, p))
        else:
            A("- `%s` → `%s` (%d bytes)" % (k, p, os.path.getsize(p)))
    if os.path.isfile(LONG_HF):
        A("- `long_heightfield` → `%s` (%d bytes, float32, %dx%d, cell %.8f m)"
          % (LONG_HF, os.path.getsize(LONG_HF), LONG_NY, LONG_NX, CELL))
    if os.path.isfile(LONG_META):
        A("- `longbed-meta` → `%s` (%d bytes)"
          % (LONG_META, os.path.getsize(LONG_META)))
    if os.path.isfile(os.path.join(OUT, "decomp.json")):
        A("- `decomp` → `%s` (%d bytes)"
          % (os.path.join(OUT, "decomp.json"),
             os.path.getsize(os.path.join(OUT, "decomp.json"))))
    if os.path.isfile(os.path.join(OUT, "equiv.json")):
        A("- `equiv` → `%s` (%d bytes)"
          % (os.path.join(OUT, "equiv.json"),
             os.path.getsize(os.path.join(OUT, "equiv.json"))))
    A("")
    A("## 9. 実行の記録")
    A("")
    A("- 全ステージ foreground（1 プロセス、`--all`）。jobkit checkpoint を"
      " stage ごとに原子的にコミット（`out/v4-2-3/checkpoint.json`）。"
      "units は stage の実行回数（本レポートの数値を出した最終通しの前に"
      " checkpoint を reset してあるので、記録された units は報告した測定に"
      "対応する。試行 1 の units は残していない）。")
    A("- 帯域・窓・実現数・判定規則は測定前にコード内で固定。後付けの調整なし"
      "（帯 0.70–0.85、適合帯 4–20 cm は v4-2.2 のまま）。")
    A("- 実行時間（stage ごと、timings.json より。report の行は本ファイル"
      "書き出し後の計測）: %s"
      % (", ".join("%s %.0fs" % (k, v) for k, v in timings.items()) or "—"))
    A("")
    A("- H 以外の量（σz・歪度・粒径・装甲比）は、床が v4-2.2 と同水準である"
      "ことの確認のために併記したもので、合否の対象ではない。")
    A("")
    md = "\n".join(L) + "\n"
    atomic_write_bytes(REPORT, md.encode("utf-8"))
    log("  wrote %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    if cp:
        cp.unit_done("report")
    return dec


def windows_local(cv):
    """Window exponents read off the stored local-slope curve.

    H_local(d) = 0.5 dlnS/dlnd, averaged over the lags inside the window
    (needs >= 2 lags).  Used for the lag-window table because a 60-bin
    geometric grid over 1.0 m has too few points for a separate fit in the
    narrow 4-6 cm window; the gate value itself is the 4-20 cm fit, unchanged.
    """
    d = np.asarray(cv["D_m"], dtype=np.float64)
    loc = np.asarray(cv["H_local"], dtype=np.float64)
    out = {}
    for a, b in WINDOWS:
        m = (d >= a) & (d <= b) & np.isfinite(loc)
        out["%.2f-%.2f" % (a, b)] = (float(loc[m].mean())
                                     if m.sum() >= 2 else None)
    return out


def interpretation_text(dec, cases, decomp=None):
    lo, hi = GATE_BAND
    r1 = cases["long_grains_r1"]["estimator_x"]["H_full_band"]
    r2 = cases["long_grains_r2"]["estimator_x"]["H_full_band"]
    sub = cases["long_substrate_extra"]["H_x_4_20_mean"]
    subr = cases["long_substrate_extra"]["H_x_4_20_range"]
    sub_r1 = cases["long_substrate_r1"]
    p30 = cases["pure_fbm_30"]["H_x_4_20_mean"]
    p30r = cases["pure_fbm_30"]["H_x_4_20_range"]
    p06 = float(np.mean([v["estimator_iso"]["H_full_band"] for v in
                         cases["pure_fbm_06"]["per_seed"].values()]))
    sub06 = cases["v422_patch_substrate"]["estimator_iso"]["H_full_band"]
    out = []
    out.append("帯状床でも**粒床**の H は帯に届かなかった（%.3f / %.3f、帯 "
               "%.2f–%.2f）。ただし内訳は「基盤の組織が足りない」ではない。"
               % (r1, r2, lo, hi))
    out.append("")
    out.append("三段の実測で切り分けた:")
    out.append("")
    out.append("1. **物差しの効果は実在した。** 基盤のみを 0.6 m パッチで測ると "
               "%.3f、3.0 m 帯状床では %.3f（3 実現 %.3f–%.3f、y 方向だけでも "
               "%.2f–%.2f）。同じ基盤が「長い床」では帯に入る。"
               % (sub06, sub, subr[0], subr[1],
                  sub_r1["estimator_y"]["H_full_band"],
                  cases["long_substrate_r2"]["estimator_y"]["H_full_band"]))
    out.append("2. **純 fBm 対照**も同じ向き: 0.6 m で %.3f、3.0 m で %.3f"
               "（4 実現 %.3f–%.3f、実現間の幅が大きい）。"
               % (p06, p30, p30r[0], p30r[1]))
    out.append("3. **粒を載せると落ちる。** 同じ帯状床で粒床は %.3f。"
               "分解（§4.1）すると 4 cm の構造関数のうち基盤が %d%%、"
               "粒の起伏が %d%%、その粒起伏の指数は %.3f（飽和したモザイク）。"
               "4–20 cm 帯の分散は粒スケールのモザイクが支配している。"
               % (r1, round(100 * decomp["split"]["0.04"]["share_form"]),
                  round(100 * decomp["split"]["0.04"]["share_relief"]),
                  decomp["terms"]["grain_relief"]["H_x_4_20"]))
    out.append("")
    out.append("したがって brief の二分法に当てはめると **B（生成側）** だが、"
               "不足しているのは「基盤の多スケール組織」ではなく**粒の集合体が "
               "4–20 cm スケールに組織を持たないこと**である。両方向は排他では"
               "なかった: 物差しの低バイアス（0.6 m で約 0.26）と、粒モザイクに"
               "よる平坦化（0.80 → %.2f）が同時に効いている。" % r1)
    return "\n".join(out)


def weakness_text(cases, long_meta, decomp=None, cfg=None):
    L = []
    A = L.append
    A("- **タイル継ぎ目**: 帯状床の粒は 0.6 m タイルごとの独立抽選なので、"
      "4 本の継ぎ目でパッキングが隣を考慮していない。継ぎ目で粒が重なり"
      "うる（v4-2.2 自身も重なり 3.9% を許す構成）。影響を切り分けるため、"
      "巨礫帯（y 0.18–0.44 m）を丸ごと落とした再測定も併記した（§4 の表の"
      "直後に数値）。")
    A("- **基盤の一般化**: 3.0 m 箱の fBm は 3.0/1.5/1.0 m のモードを持つ"
      "（0.6 m 箱には存在しない）。自己アフィン則に外側スケールはないので"
      "これは「同じ則を長い箱で評価した」ことにあたるが、σ を v4-2.2 と"
      "同じ値に固定したため、長波長の分散が増えた分だけ 1 cm 級の粗さは"
      "わずかに減る。純 fBm 対照（同じ手続き）でその影響を切り分けてある。")
    if decomp:
        A("  - 実測: 3.0×0.6 m 箱の基盤は分散の **%.1f%%** が k_y≈0 の"
          "モード（流れ方向にしか変化しない成分）に乗る。4–20 cm 帯は等方成分が"
          "担うが、20 cm 超の大ラグは流れ方向でしか意味を持たない。"
          % (100.0 * decomp["form_spectrum"]["variance_share_ky_zero"]))
    A("- **低次モード**: 0.6 m 箱と 3.0 m 箱で同じ 12 本・同じ 3 波長・同じ"
      "重み。ただし位相の引き当ては別実現なので、基盤の模様そのものは別物"
      "（帯状床は r1/r2 の 2 実現、基盤のみは追加 3 実現でばらつきを出した）。")
    A("- **歪度の不一致（重要）**: 帯状床の matrix 歪度は r1 %.3f / r2 %.3f で "
      "v4-2.2 の %.3f と一致しない（基盤 %.3f / %.3f 対 %.3f）。原因は"
      "「3.0 m 箱では分散が 3 m モードに寄るため、二次の歪み変換 c=0.55 の"
      "効き方が変わる」こと。**歪度の帯（0.3–0.9）は本レポートの対象外**で、"
      "帯状床が v4-2.2 の歪度を再現しているとは主張しない。H の解釈に使う "
      "σz・粒径・装甲比は同水準（§3 の表）。"
      % (long_meta[0]["skew_matrix"], long_meta[1]["skew_matrix"],
         (cfg or {}).get("skew_matrix_v422", float("nan")),
         long_meta[0]["substrate_skew"], long_meta[1]["substrate_skew"],
         2.210140143732825))
    A("- **巨礫 10 個**: 各タイル 2 個を除外・インペイントしている。"
      "v4-2.2 は 2 個。除外領域は面積比で同じだが、帯状床では y 0.22–0.40 m "
      "に並ぶため、その帯を落とした感度も出した。")
    A("- **実現数**: 粒床は 2 実現のみ（各 2100 万 px、レイキャストが重い）。"
      "基盤のみと純 fBm は 4/3 実現でばらつきを出した。粒床の実現間差は"
      "§0 の 2 値で判断できる範囲に留まる。")
    A("- **H の定義**: 4–20 cm の両対数傾きの 1/2。折れ点則も v4-2.2 と同一。"
      "定義を変えれば値は動く（v4-2.2 §2.1 と同じ開示）。")
    A("- **試行 1 の破棄（実行の開示）**: 最初の通し実行（04:57–05:08）は"
      "低次モードの一般化が誤っていた（長い箱で `i²+j² ∈ {25,50,100}` と"
      "置いたため物理波長が 0.3–0.6 m ではなく 0.12–0.19 m になった）。"
      "その数値は成果物に残していない（ログは `run-attempt1-discarded.log`）。"
      "修正は「v4-2.2 の 12 モードを格子ごと (5i, j) に写す」で、equiv により"
      " v4-2.2 の関数との一致を機械確認した。**帯も判定規則も動かしていない**。")
    A("- 主張の範囲: 本レポートは **H という 1 つの統計量の切り分け**しか"
      "主張しない。床の見え方・輸送・画作りは対象外（レンダリングしていない）。")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--equiv", action="store_true")
    ap.add_argument("--bed", action="store_true")
    ap.add_argument("--h", action="store_true")
    ap.add_argument("--decomp", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    timings = {}
    state = {}

    def tick(name, t0):
        timings[name] = float(time.time() - t0)
        try:
            cur = jload(TIMINGS) if os.path.isfile(TIMINGS) else {}
        except Exception:
            cur = {}
        cur[name] = round(timings[name], 1)
        jdump(TIMINGS, cur)

    if args.all or args.stage0:
        t0 = time.time()
        stage0(cp)
        tick("stage0", t0)
    if args.all or args.equiv:
        t0 = time.time()
        state["equiv"] = stage_equiv(cp)
        tick("equiv", t0)
    if args.all or args.bed:
        t0 = time.time()
        cfg, patch, longbeds = stage_bed(cp)
        state.update({"cfg": cfg, "patch": patch, "longbeds": longbeds})
        tick("bed", t0)
    if args.all or args.h:
        if "cfg" not in state:
            cfg = load_cfg()
            patch = build_patch_controls(cfg)
            longbeds = []
            for i, seed in enumerate(LONG_SEEDS):
                longbeds.append(build_long(cfg, seed, keep=(i == 0)))
            state.update({"cfg": cfg, "patch": patch, "longbeds": longbeds})
        t0 = time.time()
        cases, curves = stage_h(state["cfg"], state["patch"],
                                state["longbeds"], cp)
        cases = stage_h2(state["cfg"], cases, curves, cp)
        state["cases"] = cases
        tick("h", t0)
    if args.all or args.gates:
        t0 = time.time()
        if "cases" not in state:
            state["cases"] = jload(CURVES)["cases"]
        if "equiv" not in state:
            state["equiv"] = jload(os.path.join(OUT, "equiv.json"))
        state["gates"] = stage_gates(cp, state.get("cfg") or load_cfg(),
                                     state["cases"], state["equiv"])
        tick("gates", t0)
    if args.all or args.decomp:
        t0 = time.time()
        state["decomp"] = stage_decomp(state.get("cfg") or load_cfg(), cp)
        tick("decomp", t0)
    if args.all or args.report:
        t0 = time.time()
        cfg = state.get("cfg") or load_cfg()
        cases = state.get("cases") or jload(CURVES)["cases"]
        long_meta = ([lb["meta"] for lb in state["longbeds"]]
                     if state.get("longbeds") else
                     [jload(LONG_META), jload(os.path.join(
                         OUT, "longbed-meta-2.json"))])
        patch_meta = (state["patch"]["meta"] if state.get("patch")
                      else build_patch_controls(cfg)["meta"])
        if "gates" not in state:
            state["gates"] = stage_gates(cp, cfg, cases,
                                         state.get("equiv") or jload(
                                             os.path.join(OUT, "equiv.json")))
        dec = stage_report(cp, cfg, patch_meta, long_meta, cases,
                           jload(CURVES), state["gates"],
                           state.get("equiv") or jload(
                               os.path.join(OUT, "equiv.json")),
                           timings, cfg["D50_surface_measured_m"])
        tick("report", t0)
        log("DECISION: %s" % dec["verdict"])
        log("  H long r1=%.4f r2=%.4f mean=%.4f (band %.2f-%.2f)"
            % (dec["H_long_r1"], dec["H_long_r2"], dec["H_long_mean"],
               GATE_BAND[0], GATE_BAND[1]))
    log("done. checkpoint: %s" % CP_PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
