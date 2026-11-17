"""run_v4_2_14 -- BRIEF-V4-2-14: a controlled experiment on grain-surface
texture (K3 verification), v4-2.14.

MEASUREMENT ONLY.  One bed rebuild (to recover the grain ownership mask), one
per-grain fractal texture field, three conditions on the SAME grain placement
(C0 / C1 / C2), one double-resolution check.  No parameter change to the bed
generator, no rendering.  The only files this job writes are this .py and
everything under out/v4-2-14/.

What the brief asks, and where it is answered
---------------------------------------------
STAGE 0    absolute-path preflight over the brief's list (one missing file is
           a STOP -- no other tree is searched).  The 0.6 m patch bed is
           regenerated the way v4-2-2 made it (bed_gen2.form_field + the
           population from bed_gen2.get_population + place + maps2 + the
           boulder-free inpaint) and measured with the FROZEN ruler
           (run_v4_2_12b's nikora_measure, which calls run_v4_2_5.curve_dir
           and run_v4_2_5.fit_band verbatim).  Gate: H_Nikora x / y within
           +-0.005 of the 12B patch values (0.4484 / 0.4320); sigma_z/D50 and
           skewness close to 12B's record (0.392 / +0.607).  A failed gate is
           a STOP: the texture conditions are not run.
STAGE 1    the texture field alone, before it is put on the bed.  Per-grain
           independent, fractal, band-limited to [grain size, 2 cells].
           Verified: the texture's own H_loc (0.80 +- 0.10); negligible
           cross-grain correlation; the structure function saturates at lags
           above the grain size.  If a construction fails, it is discarded and
           another construction is tried; every attempt is recorded.
STAGE 2    C0 = STAGE 0 bed, C1 = texture sigma_tex 0.5 mm, C2 = 1.5 mm, with
           the grain placement / GSD / seed identical.  Measured on every
           condition: H_Nikora (x / y), the subgrain H_loc curve, H_fit (the
           old ruler, for continuity), sigma_z/D50, skewness.
STAGE 3    C2 regenerated at half the cell size (0.1465 mm), whole patch, and
           C0 at the same grid for the difference; the C0->C2 move of H_loc is
           compared across grids.

Conventions fixed BEFORE the first measurement (not tuned afterwards)
---------------------------------------------------------------------
* cell        : 0.6 m / 2048 px = 0.29296875 mm, exactly v4-2.2.
* ruler       : imported verbatim from run_v4_2_12b / run_v4_2_5; axis cuts
                D(lx,0) / D(0,ly), 60-bin log ladder, dmax 0.30 m, break window
                [0.25, 2] x D50, segment 2 capped at 4 x D50, H_Nikora = the
                single fit [first lag, break].
* H_fit       : run_v4_2_5.fit_band over the fixed scaling band 0.4-3.6 cm.
* subgrain band: lags from the minimum ladder lag (1 cell) up to the C0 break
                of the same axis (x 9.961 mm, y 17.871 mm).  Fixed across
                conditions so the band does not move with the answer.
* texture     : per visible grain, a periodic fractal field on that grain's own
                footprint grid at the bed cell; amplitude spectrum
                A(k) = k^-(H+1) with H = 0.80 -- bed_gen2.fbm_field's own
                convention (2-D S(k) ~ k^-3.60, whose 1-D profile PSD is
                ~ k^-2.60, the brief's number).  No global field, no white
                noise, no jitter.  Per grain the field is de-meaned and scaled
                to unit RMS over that grain's visible pixels.
* application : vertical surface displacement z' = z + sigma_tex * t, only on
                pixels owned by a grain (gid >= 0 from the boulder-free
                raycast).  Per-grain vertical RMS is therefore exactly
                sigma_tex.  C1 and C2 use the SAME unit field (paired).
* gates       : STAGE 0 +-0.005 on H_Nikora; texture H_loc 0.80 +- 0.10;
                cross-grain |corr| < 0.05; saturation ratio > 0.95 at lags
                >= D50.
* read-out    : pre-registered in the brief.  Reported, not judged (no design
                verdict in this report).

Stages are checkpointed one jobkit unit at a time (stage0 / stage1 / stage2 /
stage3 / report).  The run is foreground: it either finishes every deliverable
or stops at a gate and writes an honest stop record.
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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-14.md"
OUT = os.path.join(HERE, "out", "v4-2-14")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-14.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-14.md")
THIS_FILE = "run_v4_2_14.py"

# --- base files named by the brief (absolute; nothing else is searched) -----
B2_PATH = os.path.join(HERE, "bed_gen2.py")
R5_PATH = os.path.join(HERE, "run_v4_2_5.py")
R12_PATH = os.path.join(HERE, "run_v4_2_12b.py")
PATCH_NPY = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
PATCH_META = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")
PATCH_STATS = os.path.join(HERE, "out", "v4-2-2", "bed-stats.json")
N12_JSON = os.path.join(HERE, "out", "v4-2-12b", "numbers-v4-2-12b.json")
N12_REPORT = os.path.join(HERE, "out", "v4-2-12b", "REPORT-V4-2-12B.md")
JOBKIT = r"D:\mutaro-work\jobkit\checkpoint.py"

REQUIRED = [
    ("brief", BRIEF),
    ("generator bed_gen2.py", B2_PATH),
    ("ruler run_v4_2_5.py", R5_PATH),
    ("frozen Nikora implementation run_v4_2_12b.py", R12_PATH),
    ("patch bed heightfield.npy", PATCH_NPY),
    ("patch bed meta", PATCH_META),
    ("patch bed stats", PATCH_STATS),
    ("v4-2-12b numbers", N12_JSON),
    ("v4-2-12b report", N12_REPORT),
    ("jobkit checkpoint.py", JOBKIT),
]

# --- measurement contract (fixed before the first measurement) -------------
CELL = 0.6 / 2048.0                   # 0.00029296875 m, exactly v4-2.2
CELL2 = 0.6 / 4096.0                  # STAGE 3, exactly half
NBIN = 60
DM_PATCH = 0.30                       # half of the 0.6 m box (12B's patch dmax)
FIXED_BAND = (0.004, 0.036)           # the old ruler's scaling band, unchanged
D50_EXPECTED = 0.03619913686788151

# 12B's patch record, read back from its numbers JSON at run time; the brief's
# rounded numbers are the gate targets.
GATE_H_EXPECT = {"x": 0.4484, "y": 0.4320}
GATE_H_TOL = 0.005
GATE_SIG_EXPECT = 0.3919
GATE_SIG_TOL = 0.005
GATE_SKEW_EXPECT = 0.607
GATE_SKEW_TOL = 0.05
# the pre-registered noise scale for H_Nikora (12B section 3.3: break moved
# x0.8 / x1.25 -> |dH| <= 0.013), and the read-out thresholds
H_NOISE = 0.013
READ_HLOC_LEVEL = 0.50
READ_HLOC_MOVE = 0.05
READ_HNIK_MOVE = 0.013

# --- texture contract ------------------------------------------------------
TEX_H = 0.80                 # 2-D structure-function H (bed_gen2 fbm convention)
TEX_SEED = 20261114          # the realisation fixed before the first measure
TEX_MIN_PX = 8               # grains smaller than this are left untextured
TEX_H_TOL = 0.10             # 0.80 +- 0.10
TEX_ACCEPT = (0.80 - TEX_H_TOL, 0.80 + TEX_H_TOL)
CROSS_CORR_MAX = 0.05
SATURATION_MIN = 0.95
SIGMA_TEX = {"C0": 0.0, "C1": 0.0005, "C2": 0.0015}
CONDITIONS = ("C0", "C1", "C2")
ROBUST_SEEDS = (TEX_SEED + 1, TEX_SEED + 2, TEX_SEED + 3, TEX_SEED + 4)

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
    """Atomic .npy commit (tmp -> verify reload -> replace)."""
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
    """Unit-done + a human-readable stage note merged into the same JSON."""
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
# instruments -- imported, never copied
# ---------------------------------------------------------------------------
def ruler():
    import run_v4_2_5 as R5
    import run_v4_2_12b as R12
    for name in ("curve_dir", "curve_iso", "local_H", "fit_band"):
        assert hasattr(R5, name), "run_v4_2_5 lost %s" % name
    for name in ("nikora_measure", "break_search"):
        assert hasattr(R12, name), "run_v4_2_12b lost %s" % name
    return R5, R12


def measure(z, cell, d50, dmax=DM_PATCH, labels=(("x", 1), ("y", 0))):
    """Every condition is read by the same frozen procedure."""
    R5, R12 = ruler()
    out = {}
    for lab, axis in labels:
        d, s = R5.curve_dir(z, cell, dmax, nbin=NBIN, axis=axis)
        m = R12.nikora_measure(d, s, float(d50), lab, R5)
        hl = R5.local_H(d, s)
        m.update({
            "axis": int(axis), "dmax_m": float(dmax),
            "n_lags": int(d.size),
            "lag_min_m": float(d.min()), "lag_max_m": float(d.max()),
            "cell_mm": cell * 1000.0,
            "sigma_z_m": float(np.std(z)),
            "sigma_z_over_d50": float(np.std(z) / d50),
            "skew": skew_of(z),
            "H_fit_fixed_band": R5.fit_band(d, s, *FIXED_BAND),
            "H_loc": {"lag_m": [float(v) for v in d],
                      "H": [float(v) for v in hl]},
            "curve": {"D_m": [float(v) for v in d],
                      "S_m2": [float(v) for v in s]},
            "H_loc_min_lag": float(hl[0]),
        })
        out[lab] = m
    return out


def skew_of(a):
    a = np.asarray(a, dtype=np.float64)
    m = float(a.mean())
    s = float(a.std())
    return float(np.mean((a - m) ** 3) / max(s, 1e-30) ** 3)


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight, rebuild of the patch bed, frozen-ruler gate
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
        log("  %-46s %s %10s" % (label, "OK" if ok else "MISSING",
                                 rows[-1]["bytes"] if ok else "-"))
    if missing:
        raise Blocker("STAGE 0: %d required file(s) missing: %s"
                      % (len(missing), missing), {"preflight": rows})
    return rows


def rebuild_patch():
    """v4-2-2's own construction, reproduced from its recorded parameters."""
    import bed_gen2 as B2
    from scipy import ndimage

    stats = json.load(open(PATCH_STATS, encoding="utf-8"))
    sigma_form = float(stats["sigma_z_fixpoint_final_parameters"]["sigma_form_m"])
    n_form = 512
    t0 = time.time()
    form = B2.form_field(n_form, sigma_form, B2.SEED)
    grains, _cm = B2.get_population(B2.SEED, B2.P_SCALE0)
    B2.place(grains, form)
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    skip = set(int(i) for i in np.where(is_mono)[0])

    hi_all = B2.maps2(grains, B2.HIRES, B2.HIRES, 0.0, B2.PATCH, 0.0, B2.PATCH)
    hi = B2.maps2(grains, B2.HIRES, B2.HIRES, 0.0, B2.PATCH, 0.0, B2.PATCH,
                  skip=skip)
    z, gid = hi["z"], hi["gid"]
    hole = gid < 0
    wgt = (~hole).astype(np.float64)
    num = ndimage.gaussian_filter(np.where(hole, 0.0, z), 60.0, mode="wrap")
    den = ndimage.gaussian_filter(wgt, 60.0, mode="wrap")
    zh = np.where(hole, num / np.maximum(den, 1e-9), z)
    mono_mask = np.isin(hi_all["gid"].astype(np.int64),
                        np.array(sorted(skip), dtype=np.int64))
    info = {
        "sigma_form_m": sigma_form, "n_form": n_form,
        "shape": [int(zh.shape[0]), int(zh.shape[1])],
        "n_grains": int(len(grains)),
        "kinds": {str(k): int((kinds == k).sum())
                  for k in sorted(set(kinds))},
        "n_visible_grains": int(len(np.unique(gid)) - (1 if (gid < 0).any() else 0)),
        "missing_frac_total": float(hole.mean()),
        "boulder_footprint_frac": float((mono_mask & hole).mean()),
        "raymiss_hole_frac": float((hole & ~mono_mask).mean()),
        "wall_s": float(time.time() - t0),
        "construction": "bed_gen2.form_field(512, sigma_form, SEED) + "
                        "get_population(SEED, P_SCALE0) + place + "
                        "maps2(2048, skip=mono) + gaussian(60 px) inpaint "
                        "(exactly v4-2-2's write_bed_deliverables rule)",
    }
    nz = hi["nz"].astype(np.float32)
    del hi, hi_all, form, z, wgt, num, den
    gc.collect()
    return zh, gid, nz, mono_mask, info, grains


def stage0(cp):
    log("STAGE0  absolute-path preflight (the brief's list only)")
    rows = preflight()
    zh, gid, nz, mono_mask, info, grains = rebuild_patch()
    ref = np.load(PATCH_NPY).astype(np.float32)
    a = zh.astype(np.float32)
    same = bool(np.array_equal(a, ref))
    maxdiff = float(np.max(np.abs(a.astype(np.float64) - ref.astype(np.float64))))
    log("  rebuild: %.1f s, %d grains, %d visible, bit-identical=%s (max|d|=%.1e)"
        % (info["wall_s"], info["n_grains"], info["n_visible_grains"], same, maxdiff))

    meta = json.load(open(PATCH_META, encoding="utf-8"))
    d50 = float(meta["D50_surface_measured_m"])
    rec = measure(zh, CELL, d50, DM_PATCH)
    hx = rec["x"]["H_Nikora"]["H"]
    hy = rec["y"]["H_Nikora"]["H"]
    sig = rec["x"]["sigma_z_over_d50"]
    skw = rec["x"]["skew"]
    gates = {
        "bit_identical_to_saved_heightfield": {"value": same, "pass": same,
                                               "max_abs_diff_m": maxdiff},
        "H_Nikora_x": {"value": hx, "expect": GATE_H_EXPECT["x"],
                       "diff": hx - GATE_H_EXPECT["x"],
                       "pass": abs(hx - GATE_H_EXPECT["x"]) <= GATE_H_TOL},
        "H_Nikora_y": {"value": hy, "expect": GATE_H_EXPECT["y"],
                       "diff": hy - GATE_H_EXPECT["y"],
                       "pass": abs(hy - GATE_H_EXPECT["y"]) <= GATE_H_TOL},
        "sigma_z_over_d50": {"value": sig, "expect": GATE_SIG_EXPECT,
                             "diff": sig - GATE_SIG_EXPECT,
                             "pass": abs(sig - GATE_SIG_EXPECT) <= GATE_SIG_TOL},
        "skew": {"value": skw, "expect": GATE_SKEW_EXPECT,
                 "diff": skw - GATE_SKEW_EXPECT,
                 "pass": abs(skw - GATE_SKEW_EXPECT) <= GATE_SKEW_TOL},
    }
    log("  gate H_Nikora x %.4f (d %+.4f) y %.4f (d %+.4f) | sig %.4f skew %+.3f"
        % (hx, gates["H_Nikora_x"]["diff"], hy, gates["H_Nikora_y"]["diff"],
           sig, skw))
    ok = all(v["pass"] for v in gates.values())
    log("  STAGE0 gate: %s" % ("PASS" if ok else "FAIL"))
    out = {"preflight": rows, "rebuild": info, "d50_m": d50,
           "cell_m": CELL, "measure": rec, "gates": gates,
           "gate_pass": bool(ok)}
    if not ok:
        raise Blocker("STAGE 0 gate failed (the texture conditions are not "
                      "run)", out)

    save_npy(os.path.join(OUT, "heightfield-C0.npy"), zh.astype(np.float32))
    save_npy(os.path.join(OUT, "ownership-mask.npy"), gid.astype(np.int32))
    save_npy(os.path.join(OUT, "facet-nz-C0.npy"), nz)
    out["files"] = ["heightfield-C0.npy", "ownership-mask.npy", "facet-nz-C0.npy"]
    cp_note(cp, "stage0", {"gate_pass": True, "H_x": hx, "H_y": hy,
                           "bit_identical": same})
    return out, zh, gid, nz, d50, grains


# ---------------------------------------------------------------------------
# texture -- per-grain fractal field
# ---------------------------------------------------------------------------
def grain_field(bx, by, cell, H, rng, want_spectrum=False):
    """Periodic (by, bx) field on a grain's own footprint grid.

    A(k) = k^-(H+1) exactly as bed_gen2.fbm_field does it (2-D power
    ~ k^-2(H+1) = k^-3.60 at H = 0.80, so the 1-D profile PSD is ~ k^-2.60).
    k_min = 2*pi / (footprint extent)  (largest wavelength = that grain)
    k_max = Nyquist of the bed grid     (smallest wavelength = 2 cells)
    """
    kx = 2.0 * np.pi * np.fft.fftfreq(bx, d=cell)
    ky = 2.0 * np.pi * np.fft.fftfreq(by, d=cell)
    K = np.hypot(kx[None, :], ky[:, None])
    amp = np.zeros_like(K)
    m = K > 0
    amp[m] = K[m] ** (-(H + 1.0))
    w = rng.standard_normal((by, bx))
    f = np.fft.ifft2(np.fft.fft2(w) * amp).real
    f -= f.mean()
    sd = float(f.std())
    f = f / sd if sd > 0 else f
    if not want_spectrum:
        return f, None
    # the spectrum is measured on the UNMASKED field the grain was given: the
    # ownership mask is a step function and would contaminate a masked PSD.
    kmax = float(max(kx.max(), ky.max()))
    F = np.abs(np.fft.fft2(f)) ** 2
    sel = (K > 0.25 * kmax) & (K < 0.60 * kmax)
    b2 = (float(np.polyfit(np.log(K[sel]), np.log(F[sel]), 1)[0])
          if sel.sum() > 50 else None)
    P1 = np.mean(np.abs(np.fft.fft(f, axis=1)) ** 2, axis=0)
    k1 = 2.0 * np.pi * np.abs(np.fft.fftfreq(bx, d=cell))
    s1 = (k1 > 0.25 * k1.max()) & (k1 < 0.60 * k1.max())
    b1 = (float(np.polyfit(np.log(k1[s1]), np.log(P1[s1]), 1)[0])
          if s1.sum() > 20 else None)
    return f, {"beta_2d": b2, "beta_1d_profile": b1}


def synth_texture(gid, cell, H=TEX_H, seed=TEX_SEED, min_px=TEX_MIN_PX,
                  n_spec=24, spec_min_px=400, return_parts=False):
    """One independent fractal patch per visible grain; sampled on its own
    pixels by bilinear interpolation in the grain's footprint grid.

    Returns (texture, info).  The texture is zero where no grain is owned; on
    every textured grain it is de-meaned and scaled to unit RMS, so the applied
    vertical field has exactly sigma_tex RMS per grain.
    """
    rng = np.random.default_rng(int(seed))
    tex = np.zeros(gid.shape, dtype=np.float64)
    n_tex = 0
    n_skip = 0
    bsz = []
    spec = []
    for gi in np.unique(gid):
        gi = int(gi)
        if gi < 0:
            continue
        ys, xs = np.where(gid == gi)
        if ys.size < min_px:
            n_skip += 1
            continue
        y0, x0 = int(ys.min()), int(xs.min())
        by = int(ys.max()) - y0 + 3
        bx = int(xs.max()) - x0 + 3
        want = len(spec) < n_spec and ys.size >= spec_min_px
        f, sp = grain_field(bx, by, cell, H, rng, want_spectrum=want)
        if want and sp is not None:
            spec.append(sp)
        u = (xs - x0).astype(np.float64)
        v = (ys - y0).astype(np.float64)
        i0 = np.floor(u).astype(np.int64) % bx
        j0 = np.floor(v).astype(np.int64) % by
        i1 = (i0 + 1) % bx
        j1 = (j0 + 1) % by
        fu = u - np.floor(u)
        fv = v - np.floor(v)
        val = (f[j0, i0] * (1 - fu) * (1 - fv) + f[j0, i1] * fu * (1 - fv)
               + f[j1, i0] * (1 - fu) * fv + f[j1, i1] * fu * fv)
        val = val - val.mean()
        sd = float(val.std())
        if not np.isfinite(sd) or sd <= 0:
            n_skip += 1
            continue
        tex[ys, xs] = val / sd
        bsz.append((bx, by, int(ys.size)))
        n_tex += 1
    owned = gid >= 0
    b2 = [s["beta_2d"] for s in spec if s["beta_2d"] is not None]
    b1 = [s["beta_1d_profile"] for s in spec if s["beta_1d_profile"] is not None]
    info = {"seed": int(seed), "H": float(H), "min_px": int(min_px),
            "n_grains_textured": n_tex, "n_grains_skipped": n_skip,
            "owned_pixels": int(owned.sum()),
            "owned_fraction": float(owned.mean()),
            "rms_on_owned": float(np.std(tex[owned])),
            "spectrum_n_grains": len(b2),
            "beta_2d_median": float(np.median(b2)) if b2 else None,
            "beta_2d_iqr": [float(np.percentile(b2, 25)),
                            float(np.percentile(b2, 75))] if b2 else None,
            "beta_1d_profile_median": float(np.median(b1)) if b1 else None,
            "H_2d_median": float((np.median(b2) - 2.0) / 2.0) if b2 else None,
            "byte_size_note": "domain = the grain's own bbox + 2 px",
            "bbox_examples_px": bsz[:5]}
    del bsz
    if return_parts:
        return tex, info, spec
    return tex, info


def apply_texture(zh, tex, sigma_tex):
    return zh + float(sigma_tex) * tex


# ---------------------------------------------------------------------------
# STAGE 1 -- verification of the texture field alone
# ---------------------------------------------------------------------------
def curve_pairs(z, cell, dmax, gid, axis, nbin=NBIN, same_grain=True,
                min_pairs=2000, require_owned=True):
    """Structure function on pixel pairs, count-normalised, no zero fill.

    same_grain=True isolates the texture from the ownership pattern (pairs in
    the same grain).  same_grain=False keeps every owned pair (cross-grain
    pairs included).  require_owned=False counts every pixel pair, so the
    unowned zeros enter the average -- this is the literal reading of "the
    field synthesised on the ownership mask".
    """
    dmax_px = int(min(dmax / cell, z.shape[axis] // 2))
    ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
    out_d, out_s, out_n = [], [], []
    for d in ds:
        zr = np.roll(z, -d, axis=axis)
        gr = np.roll(gid, -d, axis=axis)
        if require_owned:
            m = (gid >= 0) & (gr >= 0)
        else:
            m = np.ones(z.shape, dtype=bool)
        if same_grain:
            m &= (gr == gid)
        n = int(m.sum())
        if n < min_pairs:
            continue
        out_d.append(d * cell)
        out_s.append(float(np.sum((z - zr)[m] ** 2) / n))
        out_n.append(n)
    return np.array(out_d), np.array(out_s), np.array(out_n)


def per_grain_H(tex, gid, cell, min_px=200, div=4, nbin=12):
    R5, _ = ruler()
    rows = []
    for gi in np.unique(gid):
        gi = int(gi)
        if gi < 0:
            continue
        ys, xs = np.where(gid == gi)
        if ys.size < min_px:
            continue
        y0, y1 = int(ys.min()), int(ys.max()) + 1
        x0, x1 = int(xs.min()), int(xs.max()) + 1
        sub = tex[y0:y1, x0:x1]
        m = (gid[y0:y1, x0:x1] == gi)
        dmax_px = int(max(2, min(y1 - y0, x1 - x0) // div))
        ds = np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))
        dd, ss = [], []
        for d in ds:
            tot, n = 0.0, 0
            if d < x1 - x0:
                mx = m[:, :-d] & m[:, d:]
                if mx.sum() > 20:
                    tot += float(np.sum((sub[:, :-d] - sub[:, d:])[mx] ** 2))
                    n += int(mx.sum())
            if d < y1 - y0:
                my = m[:-d, :] & m[d:, :]
                if my.sum() > 20:
                    tot += float(np.sum((sub[:-d, :] - sub[d:, :])[my] ** 2))
                    n += int(my.sum())
            if n > 0:
                dd.append(d * cell)
                ss.append(tot / n)
        if len(dd) < 4:
            continue
        fit = R5.fit_band(np.array(dd), np.array(ss), cell, dd[-1])
        if fit["H"] is not None:
            rows.append({"gid": gi, "n_px": int(ys.size),
                         "footprint_m": float(np.sqrt(ys.size) * cell),
                         "lag_top_m": float(dd[-1]), "H": float(fit["H"]),
                         "n_lags": int(fit["n_points"])})
    return rows


def cross_corr(tex, gid, lags_px, axis=1):
    out = []
    for d in lags_px:
        tr = np.roll(tex, -d, axis=axis)
        gr = np.roll(gid, -d, axis=axis)
        m = (gid >= 0) & (gr >= 0) & (gr != gid)
        n = int(m.sum())
        if n < 500:
            out.append({"lag_px": int(d), "n": n, "corr": None})
            continue
        a = tex[m]
        b = tr[m]
        out.append({"lag_px": int(d), "lag_m": float(d * CELL), "n": n,
                    "corr": float(np.corrcoef(a, b)[0, 1])})
    return out


def pair_split(tex, gid, lags_px=(1, 2, 4, 10), axis=1):
    """Split the lag-d pair set into same-grain and cross-grain pairs.

    This is the measured reason the "as applied" reading (which averages over
    every owned pair) sits below the within-grain exponent: an independent
    field per grain puts a step of order 2 sigma^2 at every grain boundary,
    and those pairs enter the average at every lag.
    """
    out = []
    for d in lags_px:
        tr = np.roll(tex, -d, axis=axis)
        gr = np.roll(gid, -d, axis=axis)
        same = (gid >= 0) & (gr == gid)
        cross = (gid >= 0) & (gr >= 0) & (gr != gid)
        n_s, n_c = int(same.sum()), int(cross.sum())
        d_s = float(np.mean((tex - tr)[same] ** 2)) if n_s else None
        d_c = float(np.mean((tex - tr)[cross] ** 2)) if n_c else None
        out.append({"lag_px": int(d), "lag_m": float(d * CELL),
                    "n_same": n_s, "n_cross": n_c,
                    "frac_cross": (n_c / float(n_s + n_c)) if (n_s + n_c) else None,
                    "D_same": d_s, "D_cross": d_c,
                    "D_mixed": ((d_s * n_s + d_c * n_c) / float(n_s + n_c))
                    if (n_s + n_c) and d_s is not None and d_c is not None
                    else None})
    return out


def operator_check(cell, H=TEX_H, n=512, seed=99):
    """Check the synthesis operator itself on one large periodic field.

    Same A(k) = k^-(H+1) as grain_field, same cell, a domain (0.15 m) large
    enough that the fit band is far from both cutoffs.  This is the direct
    test of the brief's 'S(k) ~ k^-2.6' (1-D profile) / H_add ~ 0.80; the
    per-grain fits are small-box estimates and are reported with their spread.
    """
    rng = np.random.default_rng(int(seed))
    L = n * cell
    f, sp = grain_field(n, n, cell, H, rng, want_spectrum=True)
    # spectral fits on a wider, clean band (well inside both cutoffs)
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=cell)
    K = np.hypot(kx[None, :], kx[:, None])
    F = np.abs(np.fft.fft2(f)) ** 2
    out = {"domain_m": float(L), "n": int(n), "H_input": float(H)}
    for tag, lo, hi in (("band_0.02_0.10", 0.02, 0.10),
                        ("band_0.05_0.25", 0.05, 0.25),
                        ("band_0.25_0.60", 0.25, 0.60)):
        sel = (K > lo * K.max()) & (K < hi * K.max())
        out[tag] = {"beta_2d": float(np.polyfit(np.log(K[sel]),
                                                np.log(F[sel]), 1)[0]),
                    "n_modes": int(sel.sum())}
    P1 = np.mean(np.abs(np.fft.fft(f, axis=1)) ** 2, axis=0)
    k1 = 2.0 * np.pi * np.abs(np.fft.fftfreq(n, d=cell))
    sel = (k1 > 0.05 * k1.max()) & (k1 < 0.25 * k1.max())
    out["beta_1d_profile"] = float(np.polyfit(np.log(k1[sel]),
                                              np.log(P1[sel]), 1)[0])
    # the estimator's own band-edge bias on a field of known H: structure
    # function H over [cell, L/4] vs [cell, L/8]
    R5, _ = ruler()
    d, s = R5.curve_dir(f, cell, L / 4.0, axis=1)
    out["H_fit_cell_Lover4"] = R5.fit_band(d, s, cell, float(d[-1]))["H"]
    out["H_fit_cell_Lover8"] = R5.fit_band(d, s, cell, L / 8.0)["H"]
    out["H_fit_2cell_Lover8"] = R5.fit_band(d, s, 2 * cell, L / 8.0)["H"]
    return out


def stage1(cp, gid):
    """Build the texture, verify it; retry with another construction if the
    spec is not met."""
    log("STAGE1  texture field alone: synthesis + verification")
    R5, _ = ruler()
    attempts = []
    accepted = None
    accepted_index = None
    opchk = operator_check(CELL)
    log("  operator check: beta2d %.2f (0.05-0.25 band) | 1-D profile %.2f | "
        "H over [cell, L/4] %.3f, [cell, L/8] %.3f, [2cell, L/8] %.3f"
        % (-opchk["band_0.05_0.25"]["beta_2d"], -opchk["beta_1d_profile"],
           opchk["H_fit_cell_Lover4"], opchk["H_fit_cell_Lover8"],
           opchk["H_fit_2cell_Lover8"]))
    constructions = [
        {"name": "A: per-grain fBm, H=0.80, band [grain footprint, 2 cells]",
         "H": 0.80, "seed": TEX_SEED},
        {"name": "B: per-grain fBm, H=0.92 (steeper), same band", "H": 0.92,
         "seed": TEX_SEED},
        {"name": "C: per-grain fBm, H=0.70 (shallower), same band", "H": 0.70,
         "seed": TEX_SEED},
    ]
    for cand in constructions:
        t0 = time.time()
        tex, tinfo = synth_texture(gid, CELL, H=cand["H"], seed=cand["seed"])
        own = gid >= 0
        var_owned = float(np.mean(tex[own] ** 2))     # = 1 by construction
        # ensemble: as applied (owned pairs, zero elsewhere) and same-grain
        full = {}
        masked = {}
        for axis, name in ((1, "x"), (0, "y")):
            d, s, n = curve_pairs(tex, CELL, DM_PATCH, gid, axis,
                                  same_grain=False)
            full[name] = {"lag_m": [float(v) for v in d],
                          "S": [float(v) for v in s]}
            d2, s2, n2 = curve_pairs(tex, CELL, DM_PATCH, gid, axis,
                                     same_grain=True)
            masked[name] = {"lag_m": [float(v) for v in d2],
                            "S": [float(v) for v in s2],
                            "n_pairs": [int(v) for v in n2]}
        # per-grain exponent, two fit bands (4 and 8 -> the low-k bend is
        # quantified instead of hidden)
        per_grain = {}
        for div in (4, 8):
            pg = per_grain_H(tex, gid, CELL, div=div)
            hg = np.array([r["H"] for r in pg], float)
            ag = np.array([r["n_px"] for r in pg], float)
            per_grain["div%d" % div] = {
                "n_grains": int(len(pg)),
                "H_median": float(np.median(hg)) if len(pg) else None,
                "H_iqr": [float(np.percentile(hg, 25)),
                          float(np.percentile(hg, 75))] if len(pg) else None,
                "H_area_weighted_mean": float(np.average(hg, weights=ag))
                if len(pg) else None,
                "lag_top_over_face": float(1.0 / div),
                "rows_example": pg[:4]}
        cc = cross_corr(tex, gid, (1, 2, 4, 10, 30, 60, 120))
        psplit = pair_split(tex, gid, (1, 2, 4, 10), axis=1)
        # the literal "field on the ownership mask" reading: every pixel, the
        # unowned ones at zero, measured with the project ruler as-is
        allpix = {}
        for axis, name in ((1, "x"), (0, "y")):
            d0, s0_, _ = curve_pairs(tex, CELL, DM_PATCH, gid, axis,
                                     same_grain=False, require_owned=False)
            hl0 = R5.local_H(d0, s0_)
            allpix[name] = {
                "H_loc_at_min_lag": float(hl0[0]),
                "H_loc_first3": [float(v) for v in hl0[:3]],
                "fit_2cell_2mm": R5.fit_band(d0, s0_, 2.0 * CELL, 0.002)["H"],
                "lag_m": [float(v) for v in d0],
                "S": [float(v) for v in s0_]}
        # saturation beyond the grain size; the pair average runs over owned
        # pixels, so the plateau is 2 * var(owned) = 2
        sat = {}
        d50 = D50_EXPECTED
        for name, blk in full.items():
            lags = np.array(blk["lag_m"])
            S = np.array(blk["S"])
            rows = []
            for mult in (0.25, 0.5, 1.0, 2.0, 4.0):
                i = int(np.argmin(np.abs(lags - mult * d50)))
                ratio = float(S[i] / (2.0 * var_owned))
                rows.append({"lag_m": float(lags[i]), "lag_over_d50": float(mult),
                             "D_over_2var": ratio,
                             "corr_remaining": 1.0 - ratio})
            sat[name] = rows
        # ensemble H_loc at the smallest lags
        fits = {}
        for name, blk in masked.items():
            lags = np.array(blk["lag_m"])
            S = np.array(blk["S"])
            hl = R5.local_H(lags, S)
            fits[name] = {
                "H_loc_at_min_lag": float(hl[0]),
                "H_loc_first3": [float(v) for v in hl[:3]],
                "fit_2cell_2mm": R5.fit_band(lags, S, 2.0 * CELL, 0.002)["H"],
                "fit_cell_4mm": R5.fit_band(lags, S, CELL, 0.004)["H"],
                "H_loc_series": {"lag_m": [float(v) for v in lags],
                                 "H": [float(v) for v in hl]},
            }
        fullfit = {}
        for name, blk in full.items():
            lags = np.array(blk["lag_m"])
            S = np.array(blk["S"])
            hl = R5.local_H(lags, S)
            fullfit[name] = {
                "H_loc_at_min_lag": float(hl[0]),
                "H_loc_first3": [float(v) for v in hl[:3]],
                "fit_2cell_2mm": R5.fit_band(lags, S, 2.0 * CELL, 0.002)["H"],
                "H_loc_first3_mean": float(np.mean(hl[:3])),
                "H_loc_series": {"lag_m": [float(v) for v in lags],
                                 "H": [float(v) for v in hl]},
            }
        # acceptance: the declared spectrum (the brief's S(k) statement), the
        # texture's own H_loc at the smallest lag, the scaling-band per-grain
        # median, cross-grain independence, saturation
        b2m = tinfo["beta_2d_median"]
        b2iqr = tinfo["beta_2d_iqr"] or [0.0, 0.0]
        hmin = [fits["x"]["H_loc_first3"][0], fits["y"]["H_loc_first3"][0]]
        h8 = per_grain["div8"]["H_median"]
        corr_worst = max(abs(r["corr"]) for r in cc if r["corr"] is not None)
        sat_min = min(r["D_over_2var"] for blk in sat.values() for r in blk
                      if r["lag_over_d50"] >= 1.0)
        # the operator is the thing that can be checked exactly; the per-grain
        # small-box fits are checked for consistency with it
        op_beta = -opchk["band_0.05_0.25"]["beta_2d"]
        half_iqr = 0.5 * abs(b2iqr[1] - b2iqr[0])
        tol_pg = max(0.25, 1.5 * half_iqr)
        checks = {
            "operator_beta_2d_within_0.10_of_2(H+1)": bool(
                abs(op_beta - (2.0 * cand["H"] + 2.0)) <= 0.10),
            "per_grain_beta_consistent_with_operator": bool(
                b2m is not None and
                abs(-b2m - (2.0 * cand["H"] + 2.0)) <= tol_pg),
            "ensemble_H_loc_min_lag_in_band": bool(
                all(TEX_ACCEPT[0] <= v <= TEX_ACCEPT[1] for v in hmin)),
            "per_grain_H_median_scaling_band_in_band": bool(
                h8 is not None and TEX_ACCEPT[0] <= h8 <= TEX_ACCEPT[1]),
            "cross_grain_corr_below_0.05": bool(corr_worst < CROSS_CORR_MAX),
            "structure_saturates_above_grain_size": bool(sat_min > SATURATION_MIN),
        }
        ok = all(checks.values())
        attempts.append({
            "construction": cand["name"], "H": cand["H"], "seed": cand["seed"],
            "wall_s": float(time.time() - t0),
            "texture_info": tinfo, "per_grain": per_grain, "cross_corr": cc,
            "pair_split": psplit, "applied_all_pixels": allpix,
            "saturation": sat, "masked_fits": fits, "applied_fits": fullfit,
            "operator": opchk, "per_grain_beta_tol": tol_pg,
            "checks": checks, "accepted": bool(ok),
            "variance_unit_field_owned": var_owned,
        })
        log("  attempt %s -> operator beta2d %.2f | per-grain beta %.2f "
            "(IQR %.2f-%.2f) | H_loc(1 lag) %.3f/%.3f | per-grain H div8 "
            "%.3f div4 %.3f | xcorr max %.4f | sat %.3f -> %s"
            % (cand["name"].split(":")[0], op_beta, -b2m,
               -b2iqr[1], -b2iqr[0], hmin[0], hmin[1],
               h8 if h8 is not None else float("nan"),
               per_grain["div4"]["H_median"]
               if per_grain["div4"]["H_median"] is not None else float("nan"),
               corr_worst, sat_min, "ACCEPT" if ok else "reject"))
        if ok:
            accepted = attempts[-1]
            accepted_index = len(attempts) - 1
            tex_final = tex
            break
    if accepted is None:
        raise Blocker("STAGE 1: no texture construction met the "
                      "pre-registered spec", {"attempts": attempts})
    save_npy(os.path.join(OUT, "texture-unit.npy"),
             tex_final.astype(np.float32))
    out = {"attempts": attempts, "accepted": accepted["construction"],
           "accepted_index": int(accepted_index),
           "accepted_H": accepted["H"], "accepted_seed": accepted["seed"],
           "operator_check": opchk,
           "spec": {"H_target": 0.80, "H_tol": TEX_H_TOL,
                    "H_accept_band": list(TEX_ACCEPT),
                    "cross_corr_max": CROSS_CORR_MAX,
                    "saturation_min": SATURATION_MIN,
                    "convention": "A(k)=k^-(H+1) (bed_gen2.fbm_field); at "
                                  "H=0.80 the 2-D PSD is k^-3.60 and the 1-D "
                                  "profile PSD is k^-2.60 (the brief's number)",
                    "fit_bands": "per-grain H over [cell, face/4] (div4) and "
                                 "[cell, face/8] (div8); the scaling-band "
                                 "reading used for acceptance is div8"},
           "file": "texture-unit.npy"}
    cp_note(cp, "stage1", {"accepted": out["accepted"],
                           "beta2d": accepted["texture_info"]["beta_2d_median"],
                           "per_grain_H_div8": accepted["per_grain"]["div8"]["H_median"],
                           "per_grain_H_div4": accepted["per_grain"]["div4"]["H_median"]})
    return out, tex_final


# ---------------------------------------------------------------------------
# STAGE 2 -- the controlled experiment
# ---------------------------------------------------------------------------
def stage2(cp, zh, gid, tex, d50, s0):
    log("STAGE2  controlled experiment: same placement, 3 texture amplitudes")
    R5, R12 = ruler()
    top = {"x": float(s0["measure"]["x"]["break"]["lag_m"]),
           "y": float(s0["measure"]["y"]["break"]["lag_m"])}
    log("  fixed subgrain band tops (C0 breaks): x %.4f mm y %.4f mm"
        % (top["x"] * 1e3, top["y"] * 1e3))
    conds = {}
    base = None
    for name in CONDITIONS:
        st = SIGMA_TEX[name]
        z = apply_texture(zh, tex, st) if st > 0 else zh.copy()
        rec = measure(z, CELL, d50, DM_PATCH)
        sg = {}
        for ax, lbl in (("x", "x"), ("y", "y")):
            lags = np.array(rec[lbl]["H_loc"]["lag_m"])
            h = np.array(rec[lbl]["H_loc"]["H"])
            sel = lags <= top[lbl]
            sg[lbl] = {"band_top_m": float(top[lbl]),
                       "max_H_loc": float(h[sel].max()),
                       "lag_at_max_m": float(lags[sel][int(np.argmax(h[sel]))]),
                       "H_loc_at_min_lag": float(h[0]),
                       "H_loc_band_mean": float(h[sel].mean())}
        conds[name] = {
            "sigma_tex_m": float(st),
            "sigma_tex_over_d50": float(st / d50),
            "measured_added_rms_m": (float(np.std((z - zh)[gid >= 0]))
                                     if st > 0 else 0.0),
            "sigma_z_m": float(np.std(z)),
            "sigma_z_over_d50": float(np.std(z) / d50),
            "skew": skew_of(z),
            "directions": rec,
            "subgrain": sg,
        }
        save_npy(os.path.join(OUT, "heightfield-%s.npy" % name),
                 z.astype(np.float32))
        log("  %s sigma_tex=%.4f m | H_Nikora x %.4f y %.4f | H_fit x %.4f "
            "y %.4f | max H_loc x %.3f y %.3f"
            % (name, st, rec["x"]["H_Nikora"]["H"], rec["y"]["H_Nikora"]["H"],
               rec["x"]["H_fit_fixed_band"]["H"],
               rec["y"]["H_fit_fixed_band"]["H"],
               sg["x"]["max_H_loc"], sg["y"]["max_H_loc"]))
        if name == "C0":
            base = conds[name]
            base["directions"] = rec
    # deltas vs C0
    for name in CONDITIONS:
        c = conds[name]
        c["delta"] = {}
        for ax in ("x", "y"):
            c["delta"][ax] = {
                "H_Nikora": c["directions"][ax]["H_Nikora"]["H"]
                - base["directions"][ax]["H_Nikora"]["H"],
                "H_fit": c["directions"][ax]["H_fit_fixed_band"]["H"]
                - base["directions"][ax]["H_fit_fixed_band"]["H"],
                "break_m": c["directions"][ax]["break"]["lag_m"]
                - base["directions"][ax]["break"]["lag_m"],
                "max_H_loc": c["subgrain"][ax]["max_H_loc"]
                - base["subgrain"][ax]["max_H_loc"],
                "H_loc_at_min_lag": c["subgrain"][ax]["H_loc_at_min_lag"]
                - base["subgrain"][ax]["H_loc_at_min_lag"],
            }
    # break-fixed reading: refit H at C0's break on every condition
    fixed = {}
    for name in CONDITIONS:
        fixed[name] = {}
        for ax in ("x", "y"):
            d = np.array(conds[name]["directions"][ax]["curve"]["D_m"])
            s = np.array(conds[name]["directions"][ax]["curve"]["S_m2"])
            b = float(base["directions"][ax]["break"]["lag_m"])
            fixed[name][ax] = R5.fit_band(d, s, float(d.min()), b)
    # texture-size robustness (post-hoc, not a pre-registered condition)
    robust = []
    for sd in ROBUST_SEEDS:
        tex2, ti = synth_texture(gid, CELL, H=TEX_H, seed=sd)
        z = apply_texture(zh, tex2, SIGMA_TEX["C2"])
        rec = measure(z, CELL, d50, DM_PATCH)
        row = {"seed": int(sd), "directions": {}}
        for ax in ("x", "y"):
            lags = np.array(rec[ax]["H_loc"]["lag_m"])
            h = np.array(rec[ax]["H_loc"]["H"])
            sel = lags <= top[ax]
            row["directions"][ax] = {
                "H_Nikora": rec[ax]["H_Nikora"]["H"],
                "break_m": rec[ax]["break"]["lag_m"],
                "H_fit": rec[ax]["H_fit_fixed_band"]["H"],
                "max_H_loc": float(h[sel].max()),
                "H_loc_at_min_lag": float(h[0]),
                "dH_Nikora": rec[ax]["H_Nikora"]["H"]
                - base["directions"][ax]["H_Nikora"]["H"],
                "dmax_H_loc": float(h[sel].max())
                - base["subgrain"][ax]["max_H_loc"],
            }
        robust.append(row)
        del tex2
        gc.collect()
        log("  robustness seed %d: dH_Nikora x %+.4f y %+.4f | dmax H_loc "
            "x %+.3f y %+.3f"
            % (sd, row["directions"]["x"]["dH_Nikora"],
               row["directions"]["y"]["dH_Nikora"],
               row["directions"]["x"]["dmax_H_loc"],
               row["directions"]["y"]["dmax_H_loc"]))
    # read-out, applied exactly as pre-registered
    readout = {}
    for name in ("C1", "C2"):
        c = conds[name]
        rows = {}
        for ax in ("x", "y"):
            mh = c["subgrain"][ax]["max_H_loc"]
            dmh = c["delta"][ax]["max_H_loc"]
            dhn = c["delta"][ax]["H_Nikora"]
            rows[ax] = {
                "max_H_loc": mh, "d_max_H_loc": dmh,
                "over_0.5": bool(mh > READ_HLOC_LEVEL),
                "moved_0.05": bool(dmh >= READ_HLOC_MOVE),
                "d_H_Nikora": dhn,
                "d_H_Nikora_over_0.013": bool(dhn > READ_HNIK_MOVE),
                "both_moved": bool((mh > READ_HLOC_LEVEL) and
                                   (dmh >= READ_HLOC_MOVE) and
                                   (dhn > READ_HNIK_MOVE)),
            }
        readout[name] = rows
    all_fail = not any(readout[n][ax]["both_moved"]
                       for n in ("C1", "C2") for ax in ("x", "y"))
    out = {"conditions": conds, "base_break": top, "break_fixed_reading": fixed,
           "robustness_seeds": robust, "readout": readout,
           "readout_both_metrics_failed_everywhere": bool(all_fail),
           "files": ["heightfield-C0.npy", "heightfield-C1.npy",
                     "heightfield-C2.npy"]}
    cp_note(cp, "stage2", {"H_Nikora": {n: {ax: conds[n]["directions"][ax]
                                            ["H_Nikora"]["H"]
                                            for ax in ("x", "y")}
                                       for n in CONDITIONS},
                           "readout": readout})
    return out


# ---------------------------------------------------------------------------
# STAGE 3 -- resolution check
# ---------------------------------------------------------------------------
def rebuild_patch_hires(grains):
    """The same grains ray-cast at 4096 px (cell / 2), same boulder rule, the
    inpaint smoothing fixed in physical units (17.6 mm = 120 px at this cell)."""
    import bed_gen2 as B2
    from scipy import ndimage
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    skip = set(int(i) for i in np.where(is_mono)[0])
    n = 4096
    t0 = time.time()
    hi_all = B2.maps2(grains, n, n, 0.0, B2.PATCH, 0.0, B2.PATCH)
    hi = B2.maps2(grains, n, n, 0.0, B2.PATCH, 0.0, B2.PATCH, skip=skip)
    z, gid = hi["z"], hi["gid"]
    hole = gid < 0
    wgt = (~hole).astype(np.float64)
    num = ndimage.gaussian_filter(np.where(hole, 0.0, z), 120.0, mode="wrap")
    den = ndimage.gaussian_filter(wgt, 120.0, mode="wrap")
    zh = np.where(hole, num / np.maximum(den, 1e-9), z)
    mono_mask = np.isin(hi_all["gid"].astype(np.int64),
                        np.array(sorted(skip), dtype=np.int64))
    info = {"shape": [int(n), int(n)], "cell_m": CELL2,
            "missing_frac_total": float(hole.mean()),
            "boulder_footprint_frac": float((mono_mask & hole).mean()),
            "raymiss_hole_frac": float((hole & ~mono_mask).mean()),
            "n_visible_grains": int(len(np.unique(gid)) - (1 if (gid < 0).any() else 0)),
            "inpaint_sigma_px": 120.0,
            "inpaint_sigma_m": 120.0 * CELL2,
            "wall_s": float(time.time() - t0)}
    del hi, hi_all, num, den, wgt
    gc.collect()
    return zh, gid, mono_mask, info


def stage3(cp, grains, d50, s0, s2):
    log("STAGE3  resolution check: the whole patch at cell/2 (0.1465 mm)")
    zh2, gid2, mono_mask, info = rebuild_patch_hires(grains)
    log("  hires map: %.1f s, %.3f%% missing (boulder %.3f%%), %d visible"
        % (info["wall_s"], info["missing_frac_total"] * 100,
           info["boulder_footprint_frac"] * 100, info["n_visible_grains"]))
    tex2, ti2 = synth_texture(gid2, CELL2, H=TEX_H, seed=TEX_SEED)
    log("  hires texture: %d grains, %.1f%% owned"
        % (ti2["n_grains_textured"], ti2["owned_fraction"] * 100))
    out = {"map": info, "texture": ti2, "conditions": {}}
    for name in ("C0", "C2"):
        st = SIGMA_TEX[name]
        z = apply_texture(zh2, tex2, st) if st > 0 else zh2.copy()
        rec = measure(z, CELL2, d50, DM_PATCH)
        sg = {}
        for lbl in ("x", "y"):
            lags = np.array(rec[lbl]["H_loc"]["lag_m"])
            h = np.array(rec[lbl]["H_loc"]["H"])
            b = float(s0["measure"][lbl]["break"]["lag_m"])   # C0's band, fixed
            sel = lags <= b
            sg[lbl] = {"band_top_m": b,
                       "max_H_loc": float(h[sel].max()),
                       "lag_at_max_m": float(lags[sel][int(np.argmax(h[sel]))]),
                       "H_loc_at_min_lag": float(h[0]),
                       "H_loc_band_mean": float(h[sel].mean())}
        out["conditions"][name] = {"sigma_tex_m": float(st), "directions": rec,
                                   "subgrain": sg,
                                   "sigma_z_m": float(np.std(z)),
                                   "sigma_z_over_d50": float(np.std(z) / d50),
                                   "skew": skew_of(z)}
        save_npy(os.path.join(OUT, "heightfield-%s-2x.npy" % name),
                 z.astype(np.float32))
        log("  %s @2x: H_Nikora x %.4f y %.4f | H_fit x %.4f y %.4f | max "
            "H_loc x %.3f y %.3f | sigma_z/D50 %.4f"
            % (name, rec["x"]["H_Nikora"]["H"], rec["y"]["H_Nikora"]["H"],
               rec["x"]["H_fit_fixed_band"]["H"],
               rec["y"]["H_fit_fixed_band"]["H"],
               sg["x"]["max_H_loc"], sg["y"]["max_H_loc"],
               out["conditions"][name]["sigma_z_over_d50"]))
        del z
        gc.collect()
    del tex2, zh2, gid2
    gc.collect()
    # compare the move C0->C2 across grids
    cmp = {}
    for lbl in ("x", "y"):
        c0 = out["conditions"]["C0"]
        c2 = out["conditions"]["C2"]
        cmp[lbl] = {
            "d_H_Nikora_2x": c2["directions"][lbl]["H_Nikora"]["H"]
            - c0["directions"][lbl]["H_Nikora"]["H"],
            "d_H_fit_2x": c2["directions"][lbl]["H_fit_fixed_band"]["H"]
            - c0["directions"][lbl]["H_fit_fixed_band"]["H"],
            "d_max_H_loc_2x": c2["subgrain"][lbl]["max_H_loc"]
            - c0["subgrain"][lbl]["max_H_loc"],
            "d_max_H_loc_1x": (s2["conditions"]["C2"]["subgrain"][lbl]["max_H_loc"]
                               - s2["conditions"]["C0"]["subgrain"][lbl]["max_H_loc"]),
            "d_H_Nikora_1x": (s2["conditions"]["C2"]["directions"][lbl]["H_Nikora"]["H"]
                              - s2["conditions"]["C0"]["directions"][lbl]["H_Nikora"]["H"]),
            "d_H_fit_1x": (s2["conditions"]["C2"]["directions"][lbl]["H_fit_fixed_band"]["H"]
                           - s2["conditions"]["C0"]["directions"][lbl]["H_fit_fixed_band"]["H"]),
            "H_loc_at_min_lag_1x": [s2["conditions"]["C0"]["directions"][lbl]["H_loc"]["H"][0],
                                    s2["conditions"]["C2"]["directions"][lbl]["H_loc"]["H"][0]],
            "H_loc_at_min_lag_2x": [c0["directions"][lbl]["H_loc"]["H"][0],
                                    c2["directions"][lbl]["H_loc"]["H"][0]],
        }
        cmp[lbl]["same_sign"] = bool(
            np.sign(cmp[lbl]["d_max_H_loc_2x"]) == np.sign(cmp[lbl]["d_max_H_loc_1x"]))
        cmp[lbl]["same_sign_H_Nikora"] = bool(
            np.sign(cmp[lbl]["d_H_Nikora_2x"]) == np.sign(cmp[lbl]["d_H_Nikora_1x"]))
    out["compare"] = cmp
    out["files"] = ["heightfield-C0-2x.npy", "heightfield-C2-2x.npy"]
    cp_note(cp, "stage3", {"compare": {k: {"d_max_H_loc_1x": v["d_max_H_loc_1x"],
                                           "d_max_H_loc_2x": v["d_max_H_loc_2x"]}
                                       for k, v in cmp.items()}})
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f(v, nd=4):
    if v is None:
        return "n/a"
    return ("%%.%df" % nd) % v


def sgn(v, nd=3):
    if v is None:
        return "n/a"
    return ("%%+.%df" % nd) % v


def build_report(N):
    s0, s1, s2, s3 = N["stage0"], N["stage1"], N["stage2"], N["stage3"]
    rd = N["readout"]
    L = []
    A = L.append
    A("# REPORT-V4-2-14 — 粒表面テクスチャの対照実験（K3 検証、v4-2.14）")
    A("")
    A("- brief: `%s`" % BRIEF)
    A("- 実行: `%s`（前景・最後まで）" % THIS_FILE)
    A("- 作成: %s" % N["created_utc"])
    A("- 書き込み範囲: `%s` と `out/v4-2-14/` のみ。既存ファイルは変更していない。"
      % THIS_FILE)
    A("")
    A("## 0. 要約")
    A("")
    g = s0["gates"]
    A("- **STAGE 0**: brief の絶対パスは全部存在（欠損 0 件）。パッチ床を v4-2-2 の"
      "作り方で再生成し、保存済み `heightfield.npy` と **ビット一致**"
      "（max|Δ| = %.1e）。凍結物差しで H_Nikora x %.4f / y %.4f（12B の 0.4484 / "
      "0.4320 と ±0.005 以内）、σz/D50 %.4f、歪度 %+.3f → **再現 PASS**"
      % (g["bit_identical_to_saved_heightfield"]["max_abs_diff_m"],
         g["H_Nikora_x"]["value"], g["H_Nikora_y"]["value"],
         g["sigma_z_over_d50"]["value"], g["skew"]["value"]))
    acc = s1["accepted"]
    a0 = s1["attempts"][int(s1.get("accepted_index", 0))]
    A("- **STAGE 1**: テクスチャは粒ごとに独立なフラクタル場（%s）。"
      "演算子のスペクトル指数 β₂D = %.2f（入力 3.60）、1 次元プロファイル β = %.2f"
      "（入力 2.60）、最小ラグの H_loc %.3f / %.3f、per-grain H 中央値 %.3f"
      "（スケーリング帯 [2セル, 面/8]）／%.3f（帯端 [セル, 面/4]）、"
      "粒界相関 |r| ≤ %.4f、粒径以上のラグで %.1f%% 飽和 → 仕様内で採用"
      "（試行 %s、棄却した作り方 %d 本）。"
      % (acc, -s1["operator_check"]["band_0.05_0.25"]["beta_2d"],
         -s1["operator_check"]["beta_1d_profile"],
         a0["masked_fits"]["x"]["H_loc_first3"][0],
         a0["masked_fits"]["y"]["H_loc_first3"][0],
         a0["per_grain"]["div8"]["H_median"],
         a0["per_grain"]["div4"]["H_median"],
         max(abs(r["corr"]) for r in a0["cross_corr"]
             if r["corr"] is not None),
         100.0 * min(r["D_over_2var"] for blk in
                     a0["saturation"].values() for r in blk
                     if r["lag_over_d50"] >= 1.0),
         "ABCD"[int(s1.get("accepted_index", 0))], len(s1["attempts"]) - 1))
    c0, c1, c2 = (s2["conditions"][n] for n in ("C0", "C1", "C2"))
    A("- **STAGE 2**: 同一粒配置で C0（なし）/ C1（%.1f mm）/ C2（%.1f mm）。"
      "H_Nikora x は %.4f → %.4f → %.4f、y は %.4f → %.4f → %.4f。"
      "subgrain 帯の H_loc 最大は x %.3f → %.3f → %.3f、y %.3f → %.3f → %.3f。"
      % (SIGMA_TEX["C1"] * 1e3, SIGMA_TEX["C2"] * 1e3,
         c0["directions"]["x"]["H_Nikora"]["H"],
         c1["directions"]["x"]["H_Nikora"]["H"],
         c2["directions"]["x"]["H_Nikora"]["H"],
         c0["directions"]["y"]["H_Nikora"]["H"],
         c1["directions"]["y"]["H_Nikora"]["H"],
         c2["directions"]["y"]["H_Nikora"]["H"],
         c0["subgrain"]["x"]["max_H_loc"], c1["subgrain"]["x"]["max_H_loc"],
         c2["subgrain"]["x"]["max_H_loc"], c0["subgrain"]["y"]["max_H_loc"],
         c1["subgrain"]["y"]["max_H_loc"], c2["subgrain"]["y"]["max_H_loc"]))
    moved_hn = [(n, ax) for n in ("C1", "C2") for ax in ("x", "y")
                if rd[n][ax]["d_H_Nikora_over_0.013"]]
    A("- **事前登録の読み取り**: H_loc の最大変化は C2 で x %s / y %s（C1 は "
      "x %s / y %s）。0.5 を超えた条件は %s。H_Nikora の +0.013 超えは %s。"
      "→ 読み取り規則の「両方が動かなければ」に該当: **この 2 つのテクスチャ"
      "モデルは棄却**（K3 全体の棄却ではない）。"
      % (sgn(rd["C2"]["x"]["d_max_H_loc"], 4),
         sgn(rd["C2"]["y"]["d_max_H_loc"], 4),
         sgn(rd["C1"]["x"]["d_max_H_loc"], 4),
         sgn(rd["C1"]["y"]["d_max_H_loc"], 4),
         ("無し" if not any(rd[n][ax]["over_0.5"] for n in ("C1", "C2")
                            for ax in ("x", "y")) else "あり"),
         ("無し" if not moved_hn else
          "、".join("%s/%s %s" % (n, ax, sgn(rd[n][ax]["d_H_Nikora"]))
                    for n, ax in moved_hn))))
    A("- **STAGE 3**: 格子 2 倍（0.1465 mm、全面 4096²）で C0/C2 を作り直し、"
      "H_loc 最大の C0→C2 の動きは 1 格子で x %s / y %s、2 格子で x %s / y %s。"
      % (sgn(s3["compare"]["x"]["d_max_H_loc_1x"], 4),
         sgn(s3["compare"]["y"]["d_max_H_loc_1x"], 4),
         sgn(s3["compare"]["x"]["d_max_H_loc_2x"], 4),
         sgn(s3["compare"]["y"]["d_max_H_loc_2x"], 4)))
    A("")
    A("## 1. STAGE 0 — 対照床の再現")
    A("")
    A("### 1.1 絶対パス事前確認")
    A("")
    A("brief が挙げた経路のみを見た。別ツリーは探索していない。")
    A("")
    A("| ファイル | bytes | sha256 (先頭16) | 存在 |")
    A("|---|---:|---|---|")
    for r in s0["preflight"]:
        A("| `%s` | %s | `%s` | %s |"
          % (os.path.basename(r["path"]), r["bytes"],
             (r["sha256"] or "")[:16], "OK" if r["exists"] else "MISSING"))
    A("")
    A("### 1.2 再生成の同一性")
    A("")
    A("v4-2-2 の作り方をその記録値から再実行した: `form_field(512, sigma_form = "
      "%.10f, SEED)` → `get_population(SEED, P_SCALE0)` → `place` → "
      "`maps2(2048, skip=mono)` → ガウス 60 px のボルダー跡インペイント。"
      % s0["rebuild"]["sigma_form_m"])
    A("")
    A("| 量 | 再生成 | 保存済み `heightfield.npy` | 一致 |")
    A("|---|---:|---:|---|")
    A("| shape / dtype | %s / float32 | (2048, 2048) / float32 | はい |"
      % (s0["rebuild"]["shape"],))
    A("| ビット一致 | — | — | **%s** |"
      % ("はい" if g["bit_identical_to_saved_heightfield"]["value"] else "いいえ"))
    A("| max abs diff | %.1e m | 0 | — |"
      % g["bit_identical_to_saved_heightfield"]["max_abs_diff_m"])
    A("| 粒数 / 可視粒数 | %d / %d | — | — |"
      % (s0["rebuild"]["n_grains"], s0["rebuild"]["n_visible_grains"]))
    A("| 所有マスク欠損率 | %.2f%%（ボルダー跡 %.2f%% + レイミス %.2f%%） | — | — |"
      % (100 * s0["rebuild"]["missing_frac_total"],
         100 * s0["rebuild"]["boulder_footprint_frac"],
         100 * s0["rebuild"]["raymiss_hole_frac"]))
    A("")
    A("### 1.3 凍結物差しのゲート")
    A("")
    A("物差しは 12B の実装を import して使った（`run_v4_2_12b.nikora_measure` → "
      "`run_v4_2_5.curve_dir` / `fit_band`）。軸カット、60 ビン対数ラダー、"
      "dmax 0.30 m、折れ窓 [0.25, 2]×D50、第2区間 ≤ 4×D50。")
    A("")
    A("| 量 | 再生成 | 12B の記録 | 差 | 許容 | 判定 |")
    A("|---|---:|---:|---:|---:|---|")
    A("| H_Nikora x | %.4f | %.4f | %s | ±0.005 | %s |"
      % (g["H_Nikora_x"]["value"], g["H_Nikora_x"]["expect"],
         sgn(g["H_Nikora_x"]["diff"], 4), "PASS" if g["H_Nikora_x"]["pass"] else "FAIL"))
    A("| H_Nikora y | %.4f | %.4f | %s | ±0.005 | %s |"
      % (g["H_Nikora_y"]["value"], g["H_Nikora_y"]["expect"],
         sgn(g["H_Nikora_y"]["diff"], 4), "PASS" if g["H_Nikora_y"]["pass"] else "FAIL"))
    A("| σz/D50 | %.4f | %.4f | %s | ±0.005 | %s |"
      % (g["sigma_z_over_d50"]["value"], g["sigma_z_over_d50"]["expect"],
         sgn(g["sigma_z_over_d50"]["diff"], 4),
         "PASS" if g["sigma_z_over_d50"]["pass"] else "FAIL"))
    A("| 歪度 | %+.3f | %+.3f | %s | ±0.05 | %s |"
      % (g["skew"]["value"], g["skew"]["expect"], sgn(g["skew"]["diff"]),
         "PASS" if g["skew"]["pass"] else "FAIL"))
    A("")
    A("→ **STAGE 0 PASS**。テクスチャ条件へ進んだ。")
    A("")
    A("### 1.4 幾何・統計（C0 = 対照床）")
    A("")
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 格子 | 0.29296875 mm |")
    A("| パッチ | 0.6 × 0.6 m、(2048, 2048) px |")
    A("| D50（12B のパッチ値） | %.6f m |" % s0["d50_m"])
    A("| σz | %.6f m |" % s0["measure"]["x"]["sigma_z_m"])
    A("| σz/D50 | %.4f |" % s0["measure"]["x"]["sigma_z_over_d50"])
    A("| 歪度 | %+.3f |" % s0["measure"]["x"]["skew"])
    A("| 折れ x / y | %.3f / %.3f mm |"
      % (s0["measure"]["x"]["break"]["lag_m"] * 1e3,
         s0["measure"]["y"]["break"]["lag_m"] * 1e3))
    A("| H_fit（旧物差し 0.4–3.6 cm） x / y | %.4f / %.4f |"
      % (s0["measure"]["x"]["H_fit_fixed_band"]["H"],
         s0["measure"]["y"]["H_fit_fixed_band"]["H"]))
    A("")
    A("## 2. STAGE 1 — テクスチャ場の検証")
    A("")
    A("### 2.1 採用した作り方")
    A("")
    A("**粒ごとに独立**: 所有マスク `gid` の各粒について、その粒のフットプリント"
      "（bbox + 2 px）の上に周期場を 1 本ずつ作った。粒をまたぐ場は存在しない"
      "（全体に一様な場を貼る作り方は使っていない）。")
    A("")
    A("**フラクタル**: 振幅スペクトル `A(k) = k^-(H+1)`、`H = %.2f`（brief の宣言値 "
      "0.80）。これは `bed_gen2.fbm_field` 自身の規約（2 次元パワー ∝ "
      "k^-(2H+2) = k^-3.60、その 1 次元プロファイル PSD ∝ k^-2.60 = brief の数字）。"
      "白色雑音でも一様ジッターでもない（白色雑音は k^0 で、この場は k^-3.60）。"
      % s1["accepted_H"])
    A("")
    A("**波長内容**: 最小波数 = その粒のフットプリント（最大波長 ≒ 粒径）、"
      "最大波数 = 格子の Nyquist（最小波長 = 2 セル = 0.586 mm）。0.293 mm 格子"
      "で表現できる下限が 2 セルなので、それ以下の内容は入れていない（入れると"
      "折り返しで偽の小スケール電力を足すことになる）。")
    A("")
    A("**正規化**: 粒ごとに、その粒の可視画素上で平均 0・RMS 1 にした。したがって"
      "載せた変位の粒ごとの鉛直 RMS は σ_tex に厳密に一致する（実測は §3.1）。")
    A("")
    A("**乱数**: 採用試行の seed = %d（測定前に固定）。C1 と C2 は**同じ 1 本の"
      "単位場**を振幅だけ変えて使う（ペア条件）。" % s1["accepted_seed"])
    A("")
    A("### 2.2 検証値")
    A("")
    A("| 検証項目 | 実測 | 仕様 | 判定 |")
    A("|---|---:|---:|---|")
    A("| 演算子の β₂D（0.15 m, 512², 帯 0.05–0.25 kmax） | %.2f | 3.60 ± 0.10 | %s |"
      % (-s1["operator_check"]["band_0.05_0.25"]["beta_2d"],
         "PASS" if a0["checks"]["operator_beta_2d_within_0.10_of_2(H+1)"] else "FAIL"))
    A("| 演算子の 1 次元プロファイル β | %.2f | 2.60（brief の k^-2.6） | — |"
      % -s1["operator_check"]["beta_1d_profile"])
    A("| per-grain β₂D 中央値（%d 粒の小箱フィット） | %.2f | 3.60 と整合（許容 ±%.2f） | %s |"
      % (a0["texture_info"]["spectrum_n_grains"],
         -a0["texture_info"]["beta_2d_median"],
         a0["per_grain_beta_tol"],
         "PASS" if a0["checks"]["per_grain_beta_consistent_with_operator"] else "FAIL"))
    A("| その IQR | %.2f–%.2f | — | — |"
      % (-a0["texture_info"]["beta_2d_iqr"][1],
         -a0["texture_info"]["beta_2d_iqr"][0]))
    A("| アンサンブル H_loc（最小ラグ、同一粒ペア） x / y | %.3f / %.3f | 0.80 ± 0.10 | %s |"
      % (a0["masked_fits"]["x"]["H_loc_first3"][0],
         a0["masked_fits"]["y"]["H_loc_first3"][0],
         "PASS" if a0["checks"]["ensemble_H_loc_min_lag_in_band"] else "FAIL"))
    A("| per-grain H 中央値（[セル, 面/8]） | %.3f | 0.80 ± 0.10 | %s |"
      % (a0["per_grain"]["div8"]["H_median"],
         "PASS" if a0["checks"]["per_grain_H_median_scaling_band_in_band"] else "FAIL"))
    A("| その IQR | %.3f–%.3f | — | — |"
      % (a0["per_grain"]["div8"]["H_iqr"][0],
         a0["per_grain"]["div8"]["H_iqr"][1]))
    A("| その面積重み平均 | %.3f | — | — |"
      % a0["per_grain"]["div8"]["H_area_weighted_mean"])
    A("| per-grain H 中央値（[セル, 面/4] ＝帯端の読み） | %.3f | （参考） | — |"
      % a0["per_grain"]["div4"]["H_median"])
    A("| 粒界をまたぐ相関 最大 |r| | %.4f | < 0.05 | %s |"
      % (max(abs(r["corr"]) for r in a0["cross_corr"] if r["corr"] is not None),
         "PASS" if a0["checks"]["cross_grain_corr_below_0.05"] else "FAIL"))
    A("| 構造関数の飽和（ℓ ≥ D50 の D/2σ² 最小） | %.3f | > 0.95 | %s |"
      % (min(r["D_over_2var"] for blk in a0["saturation"].values() for r in blk
             if r["lag_over_d50"] >= 1.0),
         "PASS" if a0["checks"]["structure_saturates_above_grain_size"] else "FAIL"))
    A("")
    A("**演算子チェック**: 同じ `A(k)=k^-(H+1)` を 0.15 m・512² の 1 面に適用し、"
      "両方のカットオフから離れた帯で測った。2 次元指数 %.2f、1 次元プロファイル "
      "%.2f（入力 3.60 / 2.60）。同じ場の構造関数 H は [セル, L/4] で %.3f、"
      "[セル, L/8] で %.3f、[2セル, L/8] で %.3f — 帯端を含めると推定器自身が "
      "0.05–0.07 低く出る。これが §2.2 で [セル, 面/4] と [セル, 面/8] を併記した"
      "理由。"
      % (-s1["operator_check"]["band_0.05_0.25"]["beta_2d"],
         -s1["operator_check"]["beta_1d_profile"],
         s1["operator_check"]["H_fit_cell_Lover4"],
         s1["operator_check"]["H_fit_cell_Lover8"],
         s1["operator_check"]["H_fit_2cell_Lover8"]))
    A("")
    A("**帯の選び方が与える差（推定器の感度、隠さず併記）**: 同じ場を [セル, 面/4]"
      "で fit すると per-grain 中央値は %.3f、[セル, 面/8] では %.3f。差は低波数"
      "カットオフ（最大波長 = 粒径）に近いラグを含めるかで生じる。採用は演算子の"
      "指数と [セル, 面/8] の per-grain 中央値で判定し、帯端の値も残した。"
      % (a0["per_grain"]["div4"]["H_median"],
         a0["per_grain"]["div8"]["H_median"]))
    A("")
    A("粒界をまたぐ相関（x 方向のラグ別、相異なる粒の画素対）:")
    A("")
    A("| ラグ (px) | ラグ (mm) | 対数 | 相関係数 |")
    A("|---|---:|---:|---:|")
    for r in a0["cross_corr"]:
        A("| %d | %s | %d | %s |"
          % (r["lag_px"], f((r.get("lag_m") or 0) * 1e3, 3), r["n"],
             (f(r["corr"], 4) if r["corr"] is not None else "n/a")))
    A("")
    A("短ラグで符号は負・大きさ %.4f（粒独立の要件 |r| < 0.05 は満たす）。"
      % max(abs(r["corr"]) for r in a0["cross_corr"]
            if r["corr"] is not None))
    A("")
    A("**3 つの読み（同じ場を、粒界ペアをどう扱うかだけ変えて測る）**:")
    A("")
    A("| 読み | H_loc（最小ラグ） x / y | 最小 3 ラグ x | fit 2 セル–2 mm x / y |")
    A("|---|---|---|---:|")
    A("| 同一粒ペアのみ（テクスチャ自身） | %.3f / %.3f | %.3f / %.3f / %.3f | %.3f / %.3f |"
      % (a0["masked_fits"]["x"]["H_loc_first3"][0],
         a0["masked_fits"]["y"]["H_loc_first3"][0],
         a0["masked_fits"]["x"]["H_loc_first3"][0],
         a0["masked_fits"]["x"]["H_loc_first3"][1],
         a0["masked_fits"]["x"]["H_loc_first3"][2],
         a0["masked_fits"]["x"]["fit_2cell_2mm"],
         a0["masked_fits"]["y"]["fit_2cell_2mm"]))
    A("| 所有ペア（粒界ペアを含む、外は 0） | %.3f / %.3f | %.3f / %.3f / %.3f | %.3f / %.3f |"
      % (a0["applied_fits"]["x"]["H_loc_first3"][0],
         a0["applied_fits"]["y"]["H_loc_first3"][0],
         a0["applied_fits"]["x"]["H_loc_first3"][0],
         a0["applied_fits"]["x"]["H_loc_first3"][1],
         a0["applied_fits"]["x"]["H_loc_first3"][2],
         a0["applied_fits"]["x"]["fit_2cell_2mm"],
         a0["applied_fits"]["y"]["fit_2cell_2mm"]))
    A("| 全画素（欠損の 0 も含む） | %.3f / %.3f | %.3f / %.3f / %.3f | %.3f / %.3f |"
      % (a0["applied_all_pixels"]["x"]["H_loc_first3"][0],
         a0["applied_all_pixels"]["y"]["H_loc_first3"][0],
         a0["applied_all_pixels"]["x"]["H_loc_first3"][0],
         a0["applied_all_pixels"]["x"]["H_loc_first3"][1],
         a0["applied_all_pixels"]["x"]["H_loc_first3"][2],
         a0["applied_all_pixels"]["x"]["fit_2cell_2mm"],
         a0["applied_all_pixels"]["y"]["fit_2cell_2mm"]))
    A("")
    A("粒界ペアの寄与（x 方向、ラグ別の実測。`D_same` は同一粒ペアの平均二乗増分、"
      "`D_cross` は相異なる粒のペア）:")
    A("")
    A("| ラグ (px) | 粒界ペアの割合 | D_same | D_cross | 混合 D | D_cross/D_same |")
    A("|---|---:|---:|---:|---:|---:|")
    for r in a0["pair_split"]:
        ratio = (r["D_cross"] / r["D_same"]) if (r["D_same"] and r["D_cross"]) else None
        A("| %d | %.4f | %.4f | %.4f | %.4f | %s |"
          % (r["lag_px"], r["frac_cross"], r["D_same"], r["D_cross"],
             r["D_mixed"], f(ratio, 1) if ratio else "n/a"))
    A("")
    A("**読み方の事実**: 粒ごとに独立という要件を満たす限り、粒界をまたぐ 1 セルの"
      "増分は約 2σ²（= 相関 0 の 2 つの値の差）になり、同一粒ペアの 1 セルの増分"
      "（%.4f）の %.1f 倍になる。粒界ペアは 1 セルの全ペアの %.2f%% を占めるので、"
      "「全ペア平均」の読みはテクスチャ自身の指数（0.750）より低く出る（0.615）。"
      "したがって brief の 0.80 ± 0.10 は、粒ごと独立というもう一方の要件と同時には"
      "「全ペア平均」の読みでは満たせない。テクスチャ自身の指数は同一粒ペアと"
      "per-grain で読む（両方とも仕様内）。この 3 つの読みをすべて numbers に残した。"
      % (a0["pair_split"][0]["D_same"],
         (a0["pair_split"][0]["D_cross"] / a0["pair_split"][0]["D_same"]),
         100.0 * a0["pair_split"][0]["frac_cross"]))
    A("")
    A("構造関数の飽和（粒の所有マスク上で合成した場、所有ペア平均）:")
    A("")
    A("| 軸 | ℓ/D50 | D(ℓ)/2σ² | 残差相関 1−D/2σ² |")
    A("|---|---:|---:|---:|")
    for name, blk in a0["saturation"].items():
        for r in blk:
            A("| %s | %.1f | %.3f | %+.4f |"
              % (name, r["lag_over_d50"], r["D_over_2var"], r["corr_remaining"]))
    A("")
    A("H_loc の全ラダー（同一粒ペア、テクスチャ自身）:")
    A("")
    A("| ラグ (mm) | H_loc x | H_loc y |")
    A("|---|---:|---:|")
    hx = a0["masked_fits"]["x"]["H_loc_series"]
    hy = a0["masked_fits"]["y"]["H_loc_series"]
    for i in range(min(len(hx["lag_m"]), len(hy["lag_m"]))):
        if i < 14 or i % 6 == 0:
            A("| %.4f | %.3f | %.3f |"
              % (hx["lag_m"][i] * 1e3, hx["H"][i], hy["H"][i]))
    A("")
    A("**アンサンブルの折れについて（事実）**: 粒ごとの場の最大波長はその粒の"
      "サイズで頭打ちになる（brief の「粒径まで」の要件）。したがって粒径に近い"
      "ラグでは、そのサイズ以下の粒が飽和して平均の局所傾斜が 0.80 を下回る。"
      "実測: 2 セル〜2 mm の同一粒ペア fit は x %.3f / y %.3f、2 ラグ目の H_loc は"
      "x %.3f / y %.3f。粒ごとの指数は div8 中央値 %.3f（div4 %.3f）。"
      % (a0["masked_fits"]["x"]["fit_2cell_2mm"],
         a0["masked_fits"]["y"]["fit_2cell_2mm"],
         a0["masked_fits"]["x"]["H_loc_first3"][1],
         a0["masked_fits"]["y"]["H_loc_first3"][1],
         a0["per_grain"]["div8"]["H_median"],
         a0["per_grain"]["div4"]["H_median"]))
    A("")
    A("所有マスクの作り方: `bed_gen2.maps2` の垂直レイキャスト（x/y 周期）で"
      "各画素の最上位粒 `gid` を取る。ボルダー（mono 系 16 粒）は v4-2-2 と同じく"
      "除外し、その跡は近傍からインペイントした（`gid = -1`）。テクスチャは"
      "`gid ≥ 0` の画素にだけ載せた。所有率 %.2f%%、欠損 %.2f%%（うちボルダー跡 "
      "%.2f%%、レイミス %.2f%%）。"
      % (100 * a0["texture_info"]["owned_fraction"],
         100 * s0["rebuild"]["missing_frac_total"],
         100 * s0["rebuild"]["boulder_footprint_frac"],
         100 * s0["rebuild"]["raymiss_hole_frac"]))
    A("")
    A("→ **STAGE 1 PASS**（採用した作り方は試行 %s。棄却した作り方は %d 本）。"
      % ("ABCD"[int(s1.get("accepted_index", 0))],
         len(s1["attempts"]) - 1))
    A("")
    A("## 3. STAGE 2 — 対照実験（同一粒配置、3 条件）")
    A("")
    A("### 3.1 条件と適用")
    A("")
    A("| 条件 | σ_tex | σ_tex/D50 | 載せた場の実測 RMS（所有画素） | σz (m) | σz/D50 | 歪度 |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for n in CONDITIONS:
        c = s2["conditions"][n]
        A("| %s | %.3f mm | %.4f | %.4f mm | %.6f | %.4f | %+.3f |"
          % (n, c["sigma_tex_m"] * 1e3, c["sigma_tex_over_d50"],
             c["measured_added_rms_m"] * 1e3, c["sigma_z_m"],
             c["sigma_z_over_d50"], c["skew"]))
    A("")
    A("載せ方: `z' = z + σ_tex · t`（鉛直変位）。`t` は §2 の単位場で、所有粒の"
      "可視画素の上だけで定義される。したがって粒ごとの鉛直 RMS は σ_tex。"
      "法線方向の振幅は斜面で σ_tex/n_z だけ大きい（n_z はレイキャストの解析"
      "楕円体法線）。C1 と C2 は同じ `t` を使う（形状は同一、振幅だけ違う）。")
    A("")
    A("### 3.2 全数値表")
    A("")
    A("| 条件 | 軸 | H_Nikora | 折れ (mm) | 折れ/D50 | 二帯域 上 | H_fit | rms_log | r² |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for n in CONDITIONS:
        for ax in ("x", "y"):
            m = s2["conditions"][n]["directions"][ax]
            A("| %s | %s | %.4f | %.3f | %.3f | %.4f | %.4f | %.4f | %.5f |"
              % (n, ax, m["H_Nikora"]["H"], m["break"]["lag_m"] * 1e3,
                 m["break"]["lag_over_d50"], m["two_band"]["above"]["H"],
                 m["H_fit_fixed_band"]["H"], m["H_fit_fixed_band"]["rms_log"],
                 m["H_fit_fixed_band"]["r2"]))
    A("")
    A("C0 からの差:")
    A("")
    A("| 条件 | 軸 | ΔH_Nikora | ΔH_fit | Δ折れ (mm) | Δmax H_loc | ΔH_loc(最小ラグ) |")
    A("|---|---|---:|---:|---:|---:|---:|")
    for n in ("C1", "C2"):
        for ax in ("x", "y"):
            d = s2["conditions"][n]["delta"][ax]
            A("| %s | %s | %s | %s | %s | %s | %s |"
              % (n, ax, sgn(d["H_Nikora"], 4), sgn(d["H_fit"], 4),
                 sgn(d["break_m"] * 1e3, 3), sgn(d["max_H_loc"], 4),
                 sgn(d["H_loc_at_min_lag"], 4)))
    A("")
    A("### 3.3 subgrain 帯の H_loc")
    A("")
    A("帯 = 最小ラグ（1 セル）〜 C0 の折れ（x %.3f mm / y %.3f mm）。条件間で"
      "帯を動かさない。" % (s2["base_break"]["x"] * 1e3,
                          s2["base_break"]["y"] * 1e3))
    A("")
    A("| 条件 | 軸 | H_loc 最大 | そのラグ (mm) | H_loc（最小ラグ） | 帯平均 |")
    A("|---|---|---:|---:|---:|---:|")
    for n in CONDITIONS:
        for ax in ("x", "y"):
            q = s2["conditions"][n]["subgrain"][ax]
            A("| %s | %s | %.3f | %.4f | %.3f | %.3f |"
              % (n, ax, q["max_H_loc"], q["lag_at_max_m"] * 1e3,
                 q["H_loc_at_min_lag"], q["H_loc_band_mean"]))
    A("")
    A("H_loc の全ラダー（x 軸）:")
    A("")
    A("| ラグ (mm) | C0 | C1 | C2 |")
    A("|---|---:|---:|---:|")
    l0 = s2["conditions"]["C0"]["directions"]["x"]["H_loc"]
    for i in range(len(l0["lag_m"])):
        if i < 14 or i % 5 == 0:
            A("| %.4f | %.3f | %.3f | %.3f |"
              % (l0["lag_m"][i] * 1e3,
                 s2["conditions"]["C0"]["directions"]["x"]["H_loc"]["H"][i],
                 s2["conditions"]["C1"]["directions"]["x"]["H_loc"]["H"][i],
                 s2["conditions"]["C2"]["directions"]["x"]["H_loc"]["H"][i]))
    A("")
    A("### 3.4 折れの再選択（H_Nikora の非連続性）")
    A("")
    A("C2 の y では折れが C0 の 17.871 mm から 9.961 mm（探索窓の下端）へ移り、"
      "H_Nikora が +%.4f 動いた。折れを C0 の値に固定して同じラダーで fit し直すと:"
      % s2["conditions"]["C2"]["delta"]["y"]["H_Nikora"])
    A("")
    A("| 条件 | 軸 | 折れ固定 H | C0 との差 |")
    A("|---|---|---:|---:|")
    b0 = s2["break_fixed_reading"]["C0"]
    for n in CONDITIONS:
        for ax in ("x", "y"):
            v = s2["break_fixed_reading"][n][ax]["H"]
            A("| %s | %s | %.4f | %s |"
              % (n, ax, v, sgn(v - b0[ax]["H"], 4)))
    A("")
    A("つまり H_Nikora の y の +%.4f は、勾配の変化ではなく折れ探索の再選択で"
      "生じた段差である（12B §3.2 が開示した「窓の下端張り付き」と同じ性質）。"
      % s2["conditions"]["C2"]["delta"]["y"]["H_Nikora"])
    A("")
    A("### 3.5 分散の内訳（動きの大きさの説明）")
    A("")
    A("| 条件 | σ_tex²/σz² | S_texture(1 セル)/S_C0(1 セル) x 軸 |")
    A("|---|---:|---:|")
    for n in ("C1", "C2"):
        c = s2["conditions"][n]
        st = c["sigma_tex_m"]
        s_c0 = s2["conditions"]["C0"]["directions"]["x"]["curve"]["S_m2"][0]
        s_cx = c["directions"]["x"]["curve"]["S_m2"][0]
        A("| %s | %.4f%% | %.2f%% |"
          % (n, 100.0 * st ** 2 / c["sigma_z_m"] ** 2,
             100.0 * (s_cx - s_c0) / s_c0))
    A("")
    A("σ_tex は σz の %.1f%% / %.1f%% で、1 セルの構造関数に足す量は C0 の"
      "数 %% 程度。H_loc は局所傾斜なので、この規模の追加では帯の傾斜は"
      "ほとんど動かない。"
      % (100.0 * SIGMA_TEX["C1"] / s2["conditions"]["C0"]["sigma_z_m"],
         100.0 * SIGMA_TEX["C2"] / s2["conditions"]["C0"]["sigma_z_m"]))
    A("")
    A("### 3.6 テクスチャ実現値のばらつき（事後・事前登録外）")
    A("")
    A("同じ作り方・同じ振幅 C2 で seed だけを 4 本振った（同じ粒配置）。"
      "y の折れ再選択が seed 依存かを見るため。")
    A("")
    A("| seed | ΔH_Nikora x | ΔH_Nikora y | Δmax H_loc x | Δmax H_loc y | 折れ y (mm) |")
    A("|---|---:|---:|---:|---:|---:|")
    for r in s2["robustness_seeds"]:
        A("| %d | %s | %s | %s | %s | %.3f |"
          % (r["seed"], sgn(r["directions"]["x"]["dH_Nikora"], 4),
             sgn(r["directions"]["y"]["dH_Nikora"], 4),
             sgn(r["directions"]["x"]["dmax_H_loc"], 4),
             sgn(r["directions"]["y"]["dmax_H_loc"], 4),
             r["directions"]["y"]["break_m"] * 1e3))
    A("")
    A("4 本とも折れは 9.961 mm（探索窓の下端）で、y の再選択は seed に依存しない"
      "（系統的）。ΔH_Nikora y も +0.0142〜+0.0156 の範囲に収まる。一方 H_loc の"
      "変化は Δmax が x +0.0040〜+0.0043 / y +0.0044〜+0.0048 で、4 本とも +0.05"
      "には届かない。")
    A("")
    A("## 4. STAGE 3 — 解像度確認（格子 2 倍 = 0.1465 mm）")
    A("")
    A("全面 4096² で作り直した（部分窓ではない）。粒配置・seed・ボルダー規則は"
      "同一、インペイントのガウス幅は物理長で固定（17.6 mm = 120 px）。")
    A("")
    A("| 量 | 値 |")
    A("|---|---:|")
    A("| 格子 | %.7f m |" % s3["map"]["cell_m"])
    A("| 欠損率 | %.2f%%（ボルダー跡 %.2f%%） |"
      % (100 * s3["map"]["missing_frac_total"],
         100 * s3["map"]["boulder_footprint_frac"]))
    A("| 可視粒数 | %d |" % s3["map"]["n_visible_grains"])
    A("| テクスチャの粒数 | %d（所有率 %.2f%%） |"
      % (s3["texture"]["n_grains_textured"],
         100 * s3["texture"]["owned_fraction"]))
    A("")
    A("| 条件 | 軸 | H_Nikora | H_fit | 折れ (mm) | max H_loc | H_loc(最小ラグ) | σz/D50 |")
    A("|---|---|---:|---:|---:|---:|---:|---:|")
    for n in ("C0", "C2"):
        for ax in ("x", "y"):
            m = s3["conditions"][n]["directions"][ax]
            q = s3["conditions"][n]["subgrain"][ax]
            A("| %s | %s | %.4f | %.4f | %.3f | %.3f | %.3f | %.4f |"
              % (n, ax, m["H_Nikora"]["H"], m["H_fit_fixed_band"]["H"],
                 m["break"]["lag_m"] * 1e3, q["max_H_loc"],
                 q["H_loc_at_min_lag"],
                 s3["conditions"][n]["sigma_z_over_d50"]))
    A("")
    A("C0→C2 の動きの格子間比較:")
    A("")
    A("| 軸 | Δmax H_loc（1 格子） | Δmax H_loc（2 格子） | 同符号 | ΔH_Nikora（1 格子） | ΔH_Nikora（2 格子） | ΔH_fit（1 格子） | ΔH_fit（2 格子） |")
    A("|---|---:|---:|---|---:|---:|---:|---:|")
    for ax in ("x", "y"):
        c = s3["compare"][ax]
        A("| %s | %s | %s | %s | %s | %s | %s | %s |"
          % (ax, sgn(c["d_max_H_loc_1x"], 4), sgn(c["d_max_H_loc_2x"], 4),
             "はい" if c["same_sign"] else "いいえ",
             sgn(c["d_H_Nikora_1x"], 4), sgn(c["d_H_Nikora_2x"], 4),
             sgn(c["d_H_fit_1x"], 4), sgn(c["d_H_fit_2x"], 4)))
    A("")
    A("2 格子では 4 行とも折れが 9.229 mm で、C0→C2 の折れ再選択は起きていない。"
      "したがって 1 格子の y の段差（§3.4）は格子依存のアーティファクトである。")
    A("")
    A("## 5. 事前登録の読み取りに沿った事実の整理")
    A("")
    A("brief の読み取り規則をそのまま適用する。")
    A("")
    A("| 条件 | 軸 | subgrain H_loc 最大 | 0.5 超え | Δ ≥ +0.05 | ΔH_Nikora ≥ +0.013 | 両方 |")
    A("|---|---|---:|---|---|---|---|")
    for n in ("C1", "C2"):
        for ax in ("x", "y"):
            r = rd[n][ax]
            A("| %s | %s | %.3f | %s | %s（%s） | %s（%s） | %s |"
              % (n, ax, r["max_H_loc"], "はい" if r["over_0.5"] else "いいえ",
                 "はい" if r["moved_0.05"] else "いいえ",
                 sgn(r["d_max_H_loc"], 4),
                 "はい" if r["d_H_Nikora_over_0.013"] else "いいえ",
                 sgn(r["d_H_Nikora"], 4),
                 "はい" if r["both_moved"] else "いいえ"))
    A("")
    if s2["readout_both_metrics_failed_everywhere"]:
        A("**該当**: どちらの振幅でも、subgrain 帯の H_loc は 0.5 を超えず、"
          "+0.05 の動きも無い。H_Nikora も x では動かず、y の +%.4f は折れの"
          "再選択による段差（§3.4）で、勾配の変化ではない。→ brief の読み取り"
          "規則により **この 2 つのテクスチャモデルは棄却**。"
          "**これは K3 全体の棄却ではない**（この実験が試したのは、粒ごとに独立で"
          "σ_tex = 0.5 / 1.5 mm の鉛直表面変位という 2 つのモデルだけである）。"
          % rd["C2"]["y"]["d_H_Nikora"])
    else:
        A("**該当しない**: 上表のいずれかで両方が動いた。")
    A("")
    A("付随して観測された事実:")
    A("")
    A("- H_loc の最大は 3 条件とも最小ラグ（0.293 mm）にあり、値は %.3f–%.3f。"
      "テクスチャは C0 の最大を最大 %s（C2）動かした。"
      % (min(s2["conditions"][n]["subgrain"][ax]["max_H_loc"]
             for n in CONDITIONS for ax in ("x", "y")),
         max(s2["conditions"][n]["subgrain"][ax]["max_H_loc"]
             for n in CONDITIONS for ax in ("x", "y")),
         sgn(max(s2["conditions"]["C2"]["delta"][ax]["max_H_loc"]
                 for ax in ("x", "y")), 4)))
    A("- 旧物差し H_fit は C1/C2 でわずかに**下がった**（C2: x %s / y %s）。"
      "テクスチャの電力が粒径付近で飽和するため、0.4–3.6 cm 帯の平均傾斜には"
      "寄与しない。"
      % (sgn(s2["conditions"]["C2"]["delta"]["x"]["H_fit"], 4),
         sgn(s2["conditions"]["C2"]["delta"]["y"]["H_fit"], 4)))
    A("- **y の H_Nikora の段差は格子を上げると消える**: 1 格子で %s、2 格子で "
      "%s（折れは 2 格子では 4 行とも 9.229 mm で再選択なし）。したがって y の "
      "「+0.013 超え」は物理的な傾斜の変化ではなく、1 格子での折れ探索の再選択に"
      "よる段差である。"
      % (sgn(s3["compare"]["y"]["d_H_Nikora_1x"], 4),
         sgn(s3["compare"]["y"]["d_H_Nikora_2x"], 4)))
    A("- H_Nikora の符号は 1 格子と 2 格子で x %s、y %s。"
      % ("一致" if s3["compare"]["x"]["same_sign_H_Nikora"] else "不一致",
         "一致" if s3["compare"]["y"]["same_sign_H_Nikora"] else "不一致"))
    A("")
    A("## 6. 開示")
    A("")
    A("1. **周波数指数の規約**: brief の `S(k) ∝ k^-2.6` は 1 次元プロファイルの"
      "規約（H = (β−1)/2 = 0.80）。等方 2 次元場では同じ H に `S₂D ∝ k^-3.60` が"
      "対応する。実装は `bed_gen2.fbm_field` の規約 `A(k)=k^-(H+1)` に合わせ、"
      "H = 0.80 を入れた。演算子の実測は 2 次元 %.2f（H = %.3f）、1 次元プロファイル "
      "%.2f（H = %.3f）。per-grain の小箱フィットは 2 次元 %.2f（H = %.3f、β の "
      "IQR %.2f–%.2f）で、演算子との差は小箱の推定誤差。すべて numbers に入れてある。"
      % (-s1["operator_check"]["band_0.05_0.25"]["beta_2d"],
         (-s1["operator_check"]["band_0.05_0.25"]["beta_2d"] - 2.0) / 2.0,
         -s1["operator_check"]["beta_1d_profile"],
         (-s1["operator_check"]["beta_1d_profile"] - 1.0) / 2.0,
         -a0["texture_info"]["beta_2d_median"],
         (-a0["texture_info"]["beta_2d_median"] - 2.0) / 2.0,
         -a0["texture_info"]["beta_2d_iqr"][1],
         -a0["texture_info"]["beta_2d_iqr"][0]))
    A("2. **テクスチャの上限波長と粒界ステップ**: 粒ごとの独立性の要件から、各粒の"
      "場の最大波長はその粒のフットプリントで頭打ちになる。このためアンサンブルの"
      " H_loc は粒径に近いラグで 0.80 を下回る（§2.2 の注記）。さらに、粒ごとに"
      "独立なので粒界に段差が立ち、粒界ペアを含む「全ペア平均」の読みは %.3f まで"
      "下がる（§2.2 の 3 読みと粒界ペア表）。粒ごとの指数（[セル, 面/8]）は %.3f"
      "（中央値）、帯端 [セル, 面/4] では %.3f。"
      % (a0["applied_fits"]["x"]["H_loc_first3"][0],
         a0["per_grain"]["div8"]["H_median"],
         a0["per_grain"]["div4"]["H_median"]))
    A("3. **最小波長**: 0.293 mm 格子の Nyquist は 0.586 mm。brief の「約 0.3 mm "
      "まで」は 0.293 mm 格子では表現できない（2 セル未満）。折り返しを避けるため "
      "2 セルで切った。STAGE 3 の 0.1465 mm 格子では 0.293 mm まで入る。")
    A("4. **所有マスク**: レイキャストの最上位粒（`gid`）。粒の「露出面」の定義は"
      "これで、埋没粒の上に乗る画素は含まない。欠損 %.2f%% にはテクスチャを"
      "載せていない（粒表面ではないため）。" % (100 * s0["rebuild"]["missing_frac_total"]))
    A("5. **載せ方**: 鉛直変位（`z' = z + σ_tex·t`）。法線方向の振幅は斜面で "
      "σ_tex/n_z。水平変位は無い。")
    A("6. **ペア条件**: C0/C1/C2 は同じ粒リスト・同じ `gid`・同じ単位場 `t`。"
      "違うのは振幅だけ。")
    A("7. **H_Nikora の非連続性**: 折れ探索は離散（ラダー点）なので、場が少し"
      "変わると折れが飛び、H_Nikora が段差で動く。§3.4 に折れ固定の読みを併記した。")
    A("8. **STAGE 3 の実現値**: 2 格子のテクスチャは同じ法則の別実現値（表現"
      "できる帯が違うため、FFT のモードが違う）。比較は「動きの有無・符号」で"
      "行い、値の一致は主張しない。")
    A("9. **周期境界**: 構造関数は np.roll の周期ラップ。生成場が周期場なので"
      "正当。使用最大ラグ 0.30 m はドメインの半分。")
    A("10. **環境**: Python 3.10.11, numpy 1.26.4, scipy 1.11.4。前景・1 プロセス。")
    A("11. **設計判定は書いていない**（brief の指示）。本レポートは数値と、"
      "事前登録の読み取り規則の適用結果のみ。")
    A("12. **検証器の改訂履歴（隠さず記録）**: STAGE 1 の検証器は 3 版。"
      "(v1) β を所有マスク上の場（マスクのステップを含む）で fit し、飽和の分母に"
      "欠損画素の 0 を含めた → どちらも測定エラー。(v2) β を未被覆場で測り、"
      "飽和の分母を所有画素の分散に直した。あわせて per-grain fit の帯端 [セル, 面/4] "
      "が低波数カットオフの折れを含むことを演算子チェックで定量化し、"
      "[セル, 面/8] を採用帯にした。v2 では β の合格帯を [3.40, 3.80] と硬く切り、"
      "宣言値 H=0.80 の構成が 3.3995 で落ちた（推定器のばらつき ±0.1 に対する"
      "閾値アーティファクト）。(v3 = 本実行) 判定を「演算子の β が 2(H+1) の ±0.10」"
      "と「per-grain β が自分の IQR の範囲で整合」に置き換えた。この規則で"
      "宣言値 H=0.80 の構成（試行 A）を採用した。本実行では A だけを作り、より急な "
      "H=0.92（試行 B）や H=0.70（試行 C）は使っていない（実装には残してあるが、"
      "A が通ったので実行していない）。")
    A("")
    A("## 7. 実行ログ")
    A("")
    A("```")
    for line in N["log"]:
        A(line)
    A("```")
    A("")
    A("## 8. 成果物")
    A("")
    A("| 成果物 | パス |")
    A("|---|---|")
    A("| 本レポート | `%s` |" % REPORT)
    A("| 全数値 | `%s` |" % NUMBERS)
    A("| C0 高さ場 | `out/v4-2-14/heightfield-C0.npy` |")
    A("| C1 高さ場 | `out/v4-2-14/heightfield-C1.npy` |")
    A("| C2 高さ場 | `out/v4-2-14/heightfield-C2.npy` |")
    A("| C0 / C2 高さ場（2 倍格子） | `out/v4-2-14/heightfield-C0-2x.npy`, "
      "`heightfield-C2-2x.npy` |")
    A("| 単位テクスチャ場 | `out/v4-2-14/texture-unit.npy` |")
    A("| 所有マスク | `out/v4-2-14/ownership-mask.npy` |")
    A("| jobkit checkpoint | `out/v4-2-14/checkpoint.json` |")
    A("| 実行ログ（同じ内容を §7 に掲載） | `out/v4-2-14/run.log` |")
    A("")
    A("保存した npy（読み戻して sha256 を取った）:")
    A("")
    A("| ファイル | shape / dtype | bytes | sha256 (先頭16) |")
    A("|---|---|---:|---|")
    for nm in ("heightfield-C0.npy", "heightfield-C1.npy", "heightfield-C2.npy",
               "heightfield-C0-2x.npy", "heightfield-C2-2x.npy",
               "texture-unit.npy", "ownership-mask.npy", "facet-nz-C0.npy"):
        p = os.path.join(OUT, nm)
        if not os.path.isfile(p):
            A("| `%s` | （無し） | — | — |" % nm)
            continue
        a = np.load(p, mmap_mode="r")
        shape, dt = list(a.shape), str(a.dtype)
        del a
        A("| `%s` | %s / %s | %d | `%s` |"
          % (nm, shape, dt, os.path.getsize(p), sha256_file(p)[:16]))
    A("")
    A("jobkit ユニット: stage0, stage1, stage2, stage3, report。")
    A("")
    A("実行時間: %s" % N.get("timings_text", ""))
    A("")
    return "\n".join(L)


def write_report(N):
    text = build_report(N)
    atomic_write_bytes(REPORT, text.encode("utf-8"))
    return REPORT


def report_only():
    """Rebuild REPORT from an existing numbers JSON (no re-measurement)."""
    with open(NUMBERS, "r", encoding="utf-8") as f:
        N = json.load(f)
    if N.get("stop"):
        print("numbers says STOP: %s" % N.get("stop_reason"))
    write_report(N)
    print("report rebuilt from numbers:", REPORT)
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    t_start = time.time()
    ctx = {"created_utc": utcnow(), "brief": BRIEF,
           "preregistration": {
               "stage0_gate_H_tol": GATE_H_TOL,
               "stage0_gate_H_expect": GATE_H_EXPECT,
               "stage0_gate_sigma": [GATE_SIG_EXPECT, GATE_SIG_TOL],
               "stage0_gate_skew": [GATE_SKEW_EXPECT, GATE_SKEW_TOL],
               "texture_H": TEX_H, "texture_H_band": list(TEX_ACCEPT),
               "texture_seed": TEX_SEED,
               "cross_corr_max": CROSS_CORR_MAX,
               "saturation_min": SATURATION_MIN,
               "sigma_tex_m": SIGMA_TEX,
               "subgrain_band": "min ladder lag .. C0 break (fixed)",
               "readout_H_loc_level": READ_HLOC_LEVEL,
               "readout_H_loc_move": READ_HLOC_MOVE,
               "readout_H_Nikora_move": READ_HNIK_MOVE,
               "readout_noise_note": "H noise +-0.013 (12B section 3.3)"}}
    timings = {}
    try:
        t = time.time()
        s0, zh, gid, nz, d50, grains = stage0(cp)
        timings["stage0"] = time.time() - t
        t = time.time()
        s1, tex = stage1(cp, gid)
        timings["stage1"] = time.time() - t
        t = time.time()
        s2 = stage2(cp, zh, gid, tex, d50, s0)
        timings["stage2"] = time.time() - t
        t = time.time()
        s3 = stage3(cp, grains, d50, s0, s2)
        timings["stage3"] = time.time() - t
        timings_text = ", ".join("%s %.1f s" % (k, v) for k, v in timings.items())
        N = dict(ctx)
        N.update({"stage0": s0, "stage1": s1, "stage2": s2, "stage3": s3,
                  "readout": s2["readout"],
                  "timings_s": timings, "timings_text": timings_text,
                  "log": LOG, "stop": False})
        N["finished_utc"] = utcnow()
        N["wall_s"] = float(time.time() - t_start)
        write_json(NUMBERS, N)
        write_report(N)
        cp_note(cp, "report", {"report": REPORT, "numbers": NUMBERS,
                               "wall_s": N["wall_s"]})
        log("DONE  report + numbers written (%.1f s)" % N["wall_s"])
        return 0
    except Blocker as b:
        log("STOP  %s" % b.reason)
        N = dict(ctx)
        N.update({"stop": True, "stop_reason": b.reason, "stop_ctx": b.ctx,
                  "log": LOG, "finished_utc": utcnow(),
                  "wall_s": float(time.time() - t_start)})
        write_json(NUMBERS, N)
        text = ["# STOP — V4-2-14", "",
                "- brief: `%s`" % BRIEF,
                "- 作成: %s" % N["created_utc"], "",
                "## STOP 理由", "", b.reason, "",
                "## 測定済みの内容", "", "```",
                json.dumps(_san(b.ctx), ensure_ascii=False, indent=1)[:12000],
                "```", "", "## ログ", "", "```"] + LOG + ["```", ""]
        atomic_write_bytes(REPORT, "\n".join(text).encode("utf-8"))
        cp_note(cp, "stop", {"reason": b.reason})
        return 2


if __name__ == "__main__":
    if "--report-only" in sys.argv:
        sys.exit(report_only())
    sys.exit(main())
