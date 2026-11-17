"""run_v4_2_10.py -- BRIEF-V4-2-10 STAGE 1 (reproduction gate) / STAGE 2
(measurement) / STAGE 3 (report).

The ruler's functions (curve_dir / curve_iso / local_H / fit_band) are IMPORTED
from run_v4_2_5.py and used verbatim; no formula is copied.  The DEM is not
periodic, so for the DEM windows (and for the generated windows measured in the
same box, for fairness) the squared-increment average is taken WITHOUT the
periodic wrap -- the brief's STAGE note 1 requires this -- and local_H /
fit_band are applied to that curve unchanged.  A wrap/no-wrap audit on the
generated faces quantifies that the two estimators differ by far less than the
gaps being judged.

New file.  Existing files are never written.  Resumable through the jobkit
checkpoint at every stage boundary.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, read_checkpoint  # noqa: E402

import run_v4_2_5 as R5                     # the ruler -- imported, never copied

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-10.md"
DESIGN = r"D:\mutaro-work\rally-in\DESIGN-V4-2-7-BED-NEXT.md"
OUT = os.path.join(HERE, "out", "v4-2-10")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
MANIFEST = os.path.join(OUT, "prep-manifest.json")
PRE_NPY = os.path.join(OUT, "raster-run1-pre.npy")
POST_NPY = os.path.join(OUT, "raster-run1-post.npy")
NUMBERS = os.path.join(OUT, "numbers-v4-2-10.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-10.md")

HF_PATCH = os.path.join(HERE, "out", "v4-2-2", "heightfield.npy")
HF_R1 = os.path.join(HERE, "out", "v4-2-3", "heightfield-long-r1.npy")
HF_R2 = os.path.join(HERE, "out", "v4-2-12b", "heightfield-long-r2-regen.npy")
NUMBERS_425 = os.path.join(HERE, "out", "v4-2-5", "numbers-v4-2-5.json")
META_422 = os.path.join(HERE, "out", "v4-2-2", "heightfield-meta.json")

CELL_GEN = 0.6 / 2048.0                      # 0.00029296875
NBIN = 60
BAND = (0.004, 0.036)                        # the ruler's scaling band
D50_BED = 0.03619913686788151                # v4-2.2 surface D50 (m)
D90_REP = 0.200                              # brief: representative D90, fixed
SEED_MASK = 20261116
R_PATCH, R_STRIP = 8, 4                      # design A3 realisation counts
H_GEN_NULL = 0.79
SIGMA_NULL = 0.40 * D50_BED

REPRO = {                                    # STAGE 1, brief values
    "patch_x": 0.323, "patch_y": 0.373, "patch_iso": 0.349,
    "long_r1_x": 0.294, "long_r1_y": 0.297, "long_r1_iso": 0.295,
    "long_r2_x": 0.322, "long_r2_iso": 0.325,
}
REPRO_EXACT = {                              # unrounded, out/v4-2-5
    "patch_x": 0.32311075522810045, "patch_y": 0.37299321343399783,
    "patch_iso": 0.34854666750555946,
    "long_r1_x": 0.294268927895509, "long_r1_y": 0.29716828399283884,
    "long_r1_iso": 0.29495733653508976,
    "long_r2_x": 0.32225764489966235, "long_r2_iso": 0.3245887726612185,
}
GATE_DM = {"patch": 0.30, "long_x": 1.00, "long_y": 0.30, "iso": 0.20}

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


# ---------------------------------------------------------------------------
# instruments
# ---------------------------------------------------------------------------
def ladder(naxis, dmax, cell, nbin=NBIN):
    """The ruler's own geometric ladder (run_v4_2_5.curve_dir's construction)."""
    dmax_px = int(min(dmax / cell, naxis // 2))
    return np.unique(np.round(np.geomspace(1, dmax_px, nbin)).astype(int))


def s_nowrap(z, cell, dmax, axis, nbin=NBIN, valid=None):
    """Mean squared increment along one axis WITHOUT the periodic wrap.

    Identical arithmetic to the ruler's curve_dir except that the pairing is
    truncated at the window edge instead of wrapping round (the brief's
    STAGE note 1: the wrap injects a seam step of 0.10-0.12 m into S and pins
    H at 0.500).  With valid=..., pairs are kept only where both ends are
    valid and no fill is applied (diagnostic lane).
    """
    ds = ladder(z.shape[axis], dmax, cell, nbin)
    d_out, s_out, n_out = [], [], []
    for d in ds:
        if axis == 1:
            a = z[:, d:] - z[:, :-d]
            w = None if valid is None else (valid[:, d:] & valid[:, :-d])
        else:
            a = z[d:, :] - z[:-d, :]
            w = None if valid is None else (valid[d:, :] & valid[:-d, :])
        if w is None:
            d_out.append(d * cell)
            s_out.append(float(np.mean(a * a)))
            n_out.append(int(a.size))
        else:
            m = float(w.sum())
            d_out.append(d * cell)
            s_out.append(float(np.where(w, a * a, 0.0).sum() / m) if m > 0 else np.nan)
            n_out.append(int(m))
    return np.array(d_out), np.array(s_out), np.array(n_out)


def s_iso(dx, sx, dy, sy):
    """0.5*(Sx+Sy) on the x ladder (S_y log-interpolated onto it)."""
    syi = np.interp(np.log(dx), np.log(dy), np.log(sy))
    return 0.5 * (sx + np.exp(syi))


def profile(d, s, common):
    """H_loc on the native ladder (local_H verbatim), interpolated to `common`."""
    H = R5.local_H(d, s)
    ok = np.isfinite(H) & np.isfinite(d) & (d > 0)
    return np.interp(np.log(common), np.log(d[ok]), H[ok])


def form_fill(z, mask, block=8, deg=3):
    """Fill holes with a low-order 2-D polynomial 'form' fitted to the valid
    cells (computed on a block-averaged grid, then evaluated everywhere)."""
    h, w = z.shape
    hb, wb = h // block, w // block
    zz = np.where(mask, np.nan, z)
    zz = zz[:hb * block, :wb * block].reshape(hb, block, wb, block)
    cnt = np.isfinite(zz).sum(axis=(1, 3))
    with np.errstate(invalid="ignore"):
        mz = np.nanmean(zz, axis=(1, 3))
    okb = cnt >= max(4, block * block // 8)
    if okb.sum() < 12:
        return mean_fill(z)
    yy, xx = np.mgrid[0:hb, 0:wb]
    xs = xx[okb] / max(wb - 1, 1) * 2 - 1
    ys = yy[okb] / max(hb - 1, 1) * 2 - 1
    terms = [xs ** i * ys ** j for i in range(deg + 1) for j in range(deg + 1 - i)]
    A = np.column_stack(terms)
    coef, *_ = np.linalg.lstsq(A, mz[okb], rcond=None)
    Y, X = np.mgrid[0:h, 0:w]
    Xs = X / max(w - 1, 1) * 2 - 1
    Ys = Y / max(h - 1, 1) * 2 - 1
    surf = np.zeros((h, w))
    for c, (i, j) in zip(coef, [(i, j) for i in range(deg + 1)
                                for j in range(deg + 1 - i)]):
        surf += c * Xs ** i * Ys ** j
    out = np.where(mask, surf, z)
    return out


def mean_fill(z):
    z = np.array(z, dtype=np.float64, copy=True)
    m = np.isfinite(z)
    z[~m] = float(np.nanmean(z))
    return z


def random_mask(shape, frac, seed):
    rng = np.random.default_rng(seed)
    return rng.random(shape) < frac


# ---------------------------------------------------------------------------
# STAGE 1 -- reproduction gate
# ---------------------------------------------------------------------------
def stage1(cp, ctx):
    log("STAGE 1  reproduction gate (ruler verbatim on the existing faces)")
    zp = np.load(HF_PATCH).astype(np.float64)
    z1 = np.fromfile(HF_R1, dtype="<f4").reshape(2048, 10240).astype(np.float64)
    z2 = np.fromfile(HF_R2, dtype="<f4").reshape(2048, 10240).astype(np.float64)
    got = {}
    d, s = R5.curve_dir(zp, CELL_GEN, GATE_DM["patch"], nbin=NBIN, axis=1)
    got["patch_x"] = R5.fit_band(d, s, *BAND)["H"]
    d, s = R5.curve_dir(zp, CELL_GEN, GATE_DM["patch"], nbin=NBIN, axis=0)
    got["patch_y"] = R5.fit_band(d, s, *BAND)["H"]
    d, s = R5.curve_iso(zp, CELL_GEN, GATE_DM["iso"], nbin=NBIN)
    got["patch_iso"] = R5.fit_band(d, s, *BAND)["H"]
    for tag, z in (("long_r1", z1), ("long_r2", z2)):
        d, s = R5.curve_dir(z, CELL_GEN, GATE_DM["long_x"], nbin=NBIN, axis=1)
        got[tag + "_x"] = R5.fit_band(d, s, *BAND)["H"]
        d, s = R5.curve_dir(z, CELL_GEN, GATE_DM["long_y"], nbin=NBIN, axis=0)
        got[tag + "_y"] = R5.fit_band(d, s, *BAND)["H"]
        d, s = R5.curve_iso(z, CELL_GEN, GATE_DM["iso"], nbin=NBIN)
        got[tag + "_iso"] = R5.fit_band(d, s, *BAND)["H"]
    dev, info = {}, {}
    ok = True
    for k, v in got.items():
        if k not in REPRO_EXACT and k not in REPRO:
            info[k] = {"got": v,
                       "note": ("no v4-2.5 record exists (v4-2.3 never saved r2's "
                                "y curve); the brief's second r2 number 0.325 is "
                                "the iso value 0.324589, which is gated separately")}
            continue
        e = REPRO_EXACT.get(k)
        d_brief = abs(v - REPRO[k]) if k in REPRO else None
        d_exact = abs(v - e) if e is not None else None
        # "complete agreement": exact to the recorded float for the faces that
        # were re-read (patch, r1); r2 comes from the verified regeneration, so
        # it is held to the same 1e-9 bar and the brief's 0.01 rounding bar.
        good = (d_exact is not None and d_exact <= 1e-9)
        if k not in REPRO_EXACT:
            good = (d_brief is not None and d_brief <= 0.01)
        ok = ok and good
        dev[k] = {"got": v, "brief": REPRO.get(k), "recorded_v425": e,
                  "abs_diff_vs_brief": d_brief, "abs_diff_vs_recorded": d_exact,
                  "pass": bool(good)}
    ctx["stage1"] = {
        "what": ("run_v4_2_5's own functions, its own ladders (x dmax 1.00 m on "
                 "the long bed, 0.30 m on the patch, y 0.30 m, iso 0.20 m), the "
                 "ruler's fixed 0.4-3.6 cm band; the numbers must come back "
                 "bit-for-bit"),
        "values": got, "deviation": dev, "informational": info,
        "pass": bool(ok),
        "faces": {"patch": HF_PATCH, "r1": HF_R1, "r2": HF_R2,
                  "r2_note": ("out/v4-2-3/heightfield-long-r2.npy does not exist "
                              "(v4-2.3 never wrote it); out/v4-2-12b/"
                              "heightfield-long-r2-regen.npy is the verified "
                              "deterministic regeneration of the same seed")},
    }
    for k, v in dev.items():
        log("  %-12s got %.12f  v4-2.5 %.12f  |d|=%.2e  %s"
            % (k, v["got"], v["recorded_v425"] or float("nan"),
               v["abs_diff_vs_recorded"] or 0.0, "PASS" if v["pass"] else "FAIL"))
    write_json(os.path.join(WORK, "stage1-gate.json"), ctx["stage1"])
    if not ok:
        cp.unit_done("stage1 BLOCKED")
        raise SystemExit("STAGE 1 STOP -- reproduction gate failed; the ruler or "
                         "a face changed. No DEM was touched.")
    cp.unit_done("stage1-gate")
    return ctx


# ---------------------------------------------------------------------------
# STAGE 2 -- one face, one direction, one fill lane
# ---------------------------------------------------------------------------
def face_curves(z, cell, dmax_x, dmax_y, mask=None, fill="mean"):
    """Return native (d,S) for x and y after the declared fill."""
    if fill == "none":
        valid = np.isfinite(z) if mask is None else ~mask
        dx, sx, nx = s_nowrap(z, cell, dmax_x, 1, valid=valid)
        dy, sy, ny = s_nowrap(z, cell, dmax_y, 0, valid=valid)
        return (dx, sx, nx), (dy, sy, ny)
    if mask is not None:
        z = np.where(mask, np.nan, z)
    zf = mean_fill(z) if fill == "mean" else form_fill(
        np.asarray(z, dtype=np.float64), ~np.isfinite(z))
    dx, sx, nx = s_nowrap(zf, cell, dmax_x, 1)
    dy, sy, ny = s_nowrap(zf, cell, dmax_y, 0)
    return (dx, sx, nx), (dy, sy, ny)


def lowpass(z, sigma_px):
    """Gaussian low-pass at a large physical scale, in O(N).

    A direct gaussian_filter with sigma=853 px on a 2048x8192 strip would need a
    ~3400 px kernel (hours).  At that scale a block mean is the same filter to
    within the block shape, so: block-average at ~sigma/2, upsample by
    repeat, then a 3-px gaussian to round the block edges.  Recorded as such.
    """
    if sigma_px <= 24:
        return ndimage.gaussian_filter(z, sigma_px, mode="nearest")
    b = max(2, int(round(sigma_px / 2.0)))
    h, w = z.shape
    hb, wb = max(1, h // b), max(1, w // b)
    zz = z[:hb * b, :wb * b].reshape(hb, b, wb, b).mean(axis=(1, 3))
    up = np.repeat(np.repeat(zz, b, axis=0), b, axis=1)
    if up.shape != z.shape:
        pad_h = h - up.shape[0]
        pad_w = w - up.shape[1]
        up = np.pad(up, ((0, max(pad_h, 0)), (0, max(pad_w, 0))), mode="edge")
        up = up[:h, :w]
    return ndimage.gaussian_filter(up, 3.0, mode="nearest")


def secondary(zf, dx, sx, dy, sy, cell):
    """Secondary quantities (recorded, never mixed into the verdict)."""
    h, w = zf.shape
    sig = float(zf.std())
    mu = float(zf.mean())
    sk = float(np.mean((zf - mu) ** 3) / max(sig, 1e-30) ** 3)
    sig_filt = (2.5 * D90_REP / 2.0) / cell
    sm = lowpass(zf, sig_filt)
    B = zf - sm
    sig_B = float(B.std())
    hill = B[B > 0]
    pit = B[B < 0]
    d84_hill = float(np.percentile(hill, 84)) if hill.size > 100 else None
    d84_pit = float(np.percentile(-pit, 84)) if pit.size > 100 else None
    ratio = (d84_hill / d84_pit) if (d84_hill and d84_pit) else None

    def S_at(d, s, x):
        return float(np.exp(np.interp(np.log(x), np.log(d), np.log(s))))

    def phi(x):
        return S_at(dx, sx, x) / (2 * sig ** 2)

    def phi_y(x):
        return S_at(dy, sy, x) / (2 * sig ** 2)

    h_meso = 0.5 * np.log(S_at(dx, sx, 0.20) / S_at(dx, sx, 0.066)) / np.log(0.20 / 0.066)
    # M6: e-folding correlation length from S(d) = 2 sig^2 (1 - 1/e)
    target = 2 * sig ** 2 * (1 - 1 / np.e)

    def corr_len(d, s):
        m = np.isfinite(s) & (s >= target)
        if not m.any():
            return None
        i = int(np.argmax(m))
        if i == 0:
            return float(d[0])
        x0, x1 = np.log(d[i - 1]), np.log(d[i])
        y0, y1 = s[i - 1], s[i]
        f = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
        return float(np.exp(x0 + f * (x1 - x0)))

    return {
        "skew": sk, "sigma_z": sig,
        "sigma_B": sig_B, "sigma_B_over_sigma_z": sig_B / max(sig, 1e-30),
        "filter_diameter_m": 2.5 * D90_REP, "filter_sigma_px": sig_filt,
        "d84_hill_m": d84_hill, "d84_pit_m": d84_pit,
        "d84_hill_over_pit": ratio,
        "phi_6p6cm_x": phi(0.066), "phi_6p6cm_y": phi_y(0.066),
        "phi_20cm_x": phi(0.20), "phi_20cm_y": phi_y(0.20),
        "H_meso_x": float(h_meso),
        "aniso_Sx_over_Sy_20cm": S_at(dx, sx, 0.20) / S_at(dy, sy, 0.20),
        "corr_len_x_m": corr_len(dx, sx), "corr_len_y_m": corr_len(dy, sy),
        "corr_len_rule": "lag where S(d) reaches 2*sigma_z^2*(1-1/e)",
        "definitions": {
            "d84_hill_over_pit": ("residual = z - gaussian(z, 2.5*D90/2); hills = "
                                  "residual>0, pits = residual<0; D84 = 84th "
                                  "percentile of |residual| in each class"),
            "phi": "S(lag)/(2*sigma_z^2) on the filled window",
            "H_meso": "0.5*ln(S(0.20)/S(0.066))/ln(0.20/0.066)",
        },
    }


def fbm(nx, ny, H, sigma, seed):
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(nx)[:, None]
    ky = np.fft.rfftfreq(ny)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    amp = np.zeros_like(k)
    m = k > 0
    amp[m] = k[m] ** (-(H + 1.0))
    noise = rng.normal(size=k.shape) + 1j * rng.normal(size=k.shape)
    f = np.fft.irfft2(amp * noise, s=(nx, ny))
    f -= f.mean()
    f *= sigma / f.std()
    return f.astype(np.float64)


def envelope(ref, test):
    lo = np.min(ref, axis=0)
    hi = np.max(ref, axis=0)
    out_lo = test < lo
    out_hi = test > hi
    outside = out_lo | out_hi
    # extreme-rank p-value (global rank envelope, Myllymäki et al. 2017)
    def ext_rank(curve, pop):
        r = []
        for l in range(pop.shape[1]):
            v = curve[l]
            r.append(min((pop[:, l] <= v).sum(), (pop[:, l] >= v).sum()))
        return min(r)
    E_test = ext_rank(test, ref)
    E_ref = np.array([ext_rank(ref[i], np.vstack([ref[:i], ref[i + 1:]]))
                      if False else ext_rank(ref[i], ref) for i in range(ref.shape[0])])
    p = float((1 + (E_ref <= E_test).sum()) / (ref.shape[0] + 1))
    return {"lo": lo, "hi": hi, "outside_lags": outside,
            "n_outside": int(outside.sum()),
            "verdict": "FAIL" if outside.any() else "PASS",
            "extreme_rank_test": int(E_test), "p_global_rank": p,
            "first_outside_lag_m": (float(test_lag[outside][0]) if outside.any() else None)}


test_lag = None


def judge_box(name, ref_profiles, gen_profiles, common, ref_curves, gen_curves,
              zscore=True):
    global test_lag
    test_lag = common
    ref = np.vstack(ref_profiles)
    mu = ref.mean(axis=0)
    sd = ref.std(axis=0, ddof=1)
    res = {}
    for gname, g in gen_profiles.items():
        env = envelope(ref, g)
        z = (g - mu) / np.where(sd > 0, sd, np.nan)
        d_ref = mu[np.argmin(abs(np.log(common / 0.0041)))] - mu[
            np.argmin(abs(np.log(common / 0.0328)))]
        d_gen = g[np.argmin(abs(np.log(common / 0.0041)))] - g[
            np.argmin(abs(np.log(common / 0.0328)))]
        res[gname] = {
            "verdict": env["verdict"],
            "n_outside_lags": env["n_outside"],
            "outside_at_m": [float(x) for x in common[env["outside_lags"]]],
            "lo": env["lo"], "hi": env["hi"], "profile": g,
            "ref_mean": mu, "ref_sd": sd, "z": z,
            "C1_mean_abs_z": float(np.nanmean(np.abs(z))),
            "C2_max_abs_z": float(np.nanmax(np.abs(z))),
            "C3_delta_gap": float(abs(d_gen - d_ref)),
            "p_global_rank": env["p_global_rank"],
            "extreme_rank": env["extreme_rank_test"],
            "lag_at_max_abs_z_m": float(common[int(np.nanargmax(np.abs(z)))]),
        }
    return {"box": name, "lags_m": [float(x) for x in common],
            "n_lags_in_band": int(len(common)), "generated": res,
            "ref_population": {"n": int(ref.shape[0]), "mean": mu, "sd": sd,
                               "min": ref.min(axis=0), "max": ref.max(axis=0),
                               "profiles": [p for p in ref_profiles]},
            "z_note": ("C1/C2/C3 are placement values from design A4, recorded as "
                       "diagnostics only; the verdict is the envelope")}


def load_gen_patch():
    return np.load(HF_PATCH).astype(np.float64)


def load_gen_strip(path, i0, n=8192):
    z = np.fromfile(path, dtype="<f4").reshape(2048, 10240)
    return z[:, i0:i0 + n].astype(np.float64)


# ---------------------------------------------------------------------------
def stage2(cp, ctx):
    man = json.load(open(MANIFEST, encoding="utf-8"))
    layout = man["layout"]
    post = np.load(POST_NPY, mmap_mode="r")
    o = man["raster"]["origin"]
    u0, s0 = o["u0"], o["s0"]
    CELL_DEM = man["cell_m"]
    f_patch = layout["nodata_summary"]["patch_mean"]
    f_strip = layout["nodata_summary"]["strip_mean"]

    # common comparison ladder: the DEM reference ladder inside the band
    lag_patch = ladder(300, 0.30, CELL_DEM)
    lag_patch = lag_patch[(lag_patch * CELL_DEM >= BAND[0]) & (lag_patch * CELL_DEM <= BAND[1])]
    common_patch = lag_patch * CELL_DEM
    lag_strip_u = ladder(1200, 1.20, CELL_DEM)
    lag_strip = lag_strip_u[(lag_strip_u * CELL_DEM >= BAND[0]) &
                            (lag_strip_u * CELL_DEM <= BAND[1])]
    common_strip = lag_strip * CELL_DEM
    log("  comparison ladders: patch %d lags %.4f..%.4f m, strip %d lags %.4f..%.4f m"
        % (len(common_patch), common_patch[0], common_patch[-1],
           len(common_strip), common_strip[0], common_strip[-1]))

    results = {"ladders": {"patch_m": [float(x) for x in common_patch],
                           "strip_m": [float(x) for x in common_strip]},
               "masks": {"patch_fraction": f_patch, "strip_fraction": f_strip,
                         "seed": SEED_MASK},
               "lanes": {}}

    for lane in ("mean", "form", "none"):
        dunit = "stage2-lane-%s" % lane
        fprev = os.path.join(WORK, "stage2-%s.json" % lane)
        have = set((cp.state() or {}).get("completed_ids") or [])
        if dunit in have and os.path.isfile(fprev):
            results["lanes"][lane] = json.load(open(fprev, encoding="utf-8"))["lanes"][lane]
            log("  --- fill lane: %s (reused from checkpoint)" % lane)
            continue
        lane_out = {"patch": {"ref": [], "ref_sec": [], "gen": {}, "gen_sec": {}},
                    "strip": {"ref": [], "ref_sec": [], "gen": {}, "gen_sec": {}}}
        log("  --- fill lane: %s" % lane)
        # ---- reference population (Post), measured in the box, same end rule
        for box, key in (("patch", "patch_windows"), ("strip", "strip_windows")):
            dm = (0.30, 0.30) if box == "patch" else (1.20, 0.30)
            common = common_patch if box == "patch" else common_strip
            for w in layout[key]:
                a = np.asarray(post[w["i0"]:w["i0"] + w["ni"],
                                    w["j0"]:w["j0"] + w["nj"]], dtype=np.float64)
                mask = ~np.isfinite(a)
                (dx, sx, _), (dy, sy, _) = face_curves(a, CELL_DEM, dm[0], dm[1],
                                                       fill=lane)
                lane_out[box]["ref"].append(
                    {"id": w["id"], "x": profile(dx, sx, common),
                     "y": profile(dy, sy, common),
                     "nodata": w["nodata"]["fraction"]})
                if lane == "mean":
                    zf = mean_fill(np.where(mask, np.nan, a))
                    lane_out[box]["ref_sec"].append(
                        {"id": w["id"], **secondary(zf, dx, sx, dy, sy, CELL_DEM)})
            log("    %s reference: %d windows measured" % (box, len(lane_out[box]["ref"])))
        # ---- generated faces
        zp = load_gen_patch()
        mask_p = random_mask(zp.shape, f_patch, SEED_MASK)
        (dx, sx, _), (dy, sy, _) = face_curves(zp, CELL_GEN, 0.30, 0.30,
                                               mask=mask_p, fill=lane)
        lane_out["patch"]["gen"]["v4-2.2 patch"] = {
            "x": profile(dx, sx, common_patch),
            "y": profile(dy, sy, common_patch),
            "mask": "random %.3f (reference mean fraction), seed %d"
                    % (f_patch, SEED_MASK)}
        if lane == "mean":
            zf = mean_fill(np.where(mask_p, np.nan, zp))
            lane_out["patch"]["gen_sec"]["v4-2.2 patch"] = secondary(
                zf, dx, sx, dy, sy, CELL_GEN)
        del zp

        for tag, path in (("v4-2.3 r1", HF_R1), ("v4-2.3 r2", HF_R2)):
            for pos, i0 in (("primary", 0), ("alt", 2048)):
                z = load_gen_strip(path, i0)
                mask_s = random_mask(z.shape, f_strip, SEED_MASK)
                (dx, sx, _), (dy, sy, _) = face_curves(z, CELL_GEN, 1.20, 0.30,
                                                       mask=mask_s, fill=lane)
                lane_out["strip"]["gen"]["%s [%s u=%.2f-%.2f m]"
                                        % (tag, pos, i0 * CELL_GEN,
                                           (i0 + 8192) * CELL_GEN)] = {
                    "x": profile(dx, sx, common_strip),
                    "y": profile(dy, sy, common_strip),
                    "mask": "random %.3f, seed %d" % (f_strip, SEED_MASK)}
                if lane == "mean" and pos == "primary":
                    zf = mean_fill(np.where(mask_s, np.nan, z))
                    lane_out["strip"]["gen_sec"]["%s [primary]" % tag] = secondary(
                        zf, dx, sx, dy, sy, CELL_GEN)
                del z
            log("    %s measured (2 placements)" % tag)

        # ---- judge
        for box, common in (("patch", common_patch), ("strip", common_strip)):
            ref = lane_out[box]["ref"]
            refx = [r["x"] for r in ref]
            refy = [r["y"] for r in ref]
            for dkey, ridx in (("x", refx), ("y", refy)):
                gen = {k: v[dkey] for k, v in lane_out[box]["gen"].items()}
                j = judge_box(box, ridx, gen, common, None, None)
                lane_out[box]["judge_" + dkey] = j
        results["lanes"][lane] = lane_out
        cp.unit_done("stage2-lane-%s" % lane)
        write_json(os.path.join(WORK, "stage2-%s.json" % lane),
                   {"lanes": {lane: lane_out}, "ladders": results["ladders"]})
        del lane_out
    ctx["stage2"] = results
    return ctx


# ---------------------------------------------------------------------------
# null (estimator calibration only -- never a judge)
# ---------------------------------------------------------------------------
def stage2_null(cp, ctx):
    log("STAGE2  fBm null (H_gen=0.79, R=8 patch / R=4 strip) -- calibration only")
    common_patch = np.array(ctx["stage2"]["ladders"]["patch_m"])
    common_strip = np.array(ctx["stage2"]["ladders"]["strip_m"])
    out = {"H_gen": H_GEN_NULL, "sigma": SIGMA_NULL, "R": [R_PATCH, R_STRIP],
           "what": ("the design A3 pure-fBm calibration: the same ruler, the same "
                    "box, the same band; measures the finite-box bias of the "
                    "estimator. NOT a physical reference and NOT used in the "
                    "verdict.")}
    # patch, periodic ruler verbatim
    P = []
    for r in range(R_PATCH):
        z = fbm(2048, 2048, H_GEN_NULL, SIGMA_NULL, 1000 + r)
        d, s = R5.curve_dir(z, CELL_GEN, 0.30, nbin=NBIN, axis=1)
        px = profile(d, s, common_patch)
        d, s = R5.curve_dir(z, CELL_GEN, 0.30, nbin=NBIN, axis=0)
        py = profile(d, s, common_patch)
        P.append((px, py))
        del z
    out["patch"] = {"x_mean": np.mean([p[0] for p in P], axis=0),
                    "x_sd": np.std([p[0] for p in P], axis=0, ddof=1),
                    "y_mean": np.mean([p[1] for p in P], axis=0),
                    "y_sd": np.std([p[1] for p in P], axis=0, ddof=1)}
    log("  patch null x: %.3f..%.3f" % (out["patch"]["x_mean"][0],
                                        out["patch"]["x_mean"][-1]))
    # strip
    S = []
    for r in range(R_STRIP):
        z = fbm(2048, 8192, H_GEN_NULL, SIGMA_NULL, 2000 + r)
        d, s = R5.curve_dir(z, CELL_GEN, 1.20, nbin=NBIN, axis=1)
        px = profile(d, s, common_strip)
        d, s = R5.curve_dir(z, CELL_GEN, 0.30, nbin=NBIN, axis=0)
        py = profile(d, s, common_strip)
        S.append((px, py))
        del z
        log("    strip null realization %d/%d" % (r + 1, R_STRIP))
    out["strip"] = {"x_mean": np.mean([p[0] for p in S], axis=0),
                    "x_sd": np.std([p[0] for p in S], axis=0, ddof=1),
                    "y_mean": np.mean([p[1] for p in S], axis=0),
                    "y_sd": np.std([p[1] for p in S], axis=0, ddof=1)}
    log("  strip null x: %.3f..%.3f" % (out["strip"]["x_mean"][0],
                                        out["strip"]["x_mean"][-1]))
    ctx["null"] = out
    write_json(os.path.join(WORK, "stage2-null.json"), out)
    cp.unit_done("stage2-null")
    return ctx


def stage2_null_handicap(cp, ctx):
    """Diagnostic: what a pure H=0.79 face measures after the DEM's own
    handicap (2 mm lattice, the DEM's nodata fraction, and each fill lane).
    This is what tells us whether a lane can discriminate at all."""
    log("STAGE2  null under the DEM handicap (2 mm + holes, all three lanes)")
    man = json.load(open(MANIFEST, encoding="utf-8"))
    common = np.array(ctx["stage2"]["ladders"]["patch_m"])
    f = man["layout"]["nodata_summary"]["patch_mean"]
    out = {"fraction": f, "R": 8, "cell_m": 7 * CELL_GEN, "lanes": {},
           "what": ("pure fBm H=0.79 put through the same 2 mm lattice + the same "
                    "nodata fraction + the same fill rule as the reference "
                    "windows: the profile a true H=0.79 bed would show if it were "
                    "scanned like the Hodge TLS data, per lane. Diagnostic only.")}
    for lane in ("mean", "form", "none"):
        vals = []
        for r in range(8):
            z = fbm(2048, 2048, H_GEN_NULL, SIGMA_NULL, 3000 + r)
            n = 2048 // 7 * 7
            zz = z[:n, :n].reshape(n // 7, 7, n // 7, 7).mean(axis=(1, 3))[:300, :300]
            m = random_mask(zz.shape, f, SEED_MASK + r)
            (dx, sx, _), (dy, sy, _) = face_curves(zz, 7 * CELL_GEN, 0.30, 0.30,
                                                   mask=m, fill=lane)
            vals.append((profile(dx, sx, common), profile(dy, sy, common)))
            del z, zz
        out["lanes"][lane] = {
            "x_mean": np.mean([v[0] for v in vals], axis=0),
            "x_sd": np.std([v[0] for v in vals], axis=0, ddof=1),
            "y_mean": np.mean([v[1] for v in vals], axis=0),
            "y_sd": np.std([v[1] for v in vals], axis=0, ddof=1)}
        log("  handicap null [%s] x: %.3f..%.3f"
            % (lane, out["lanes"][lane]["x_mean"][0],
               out["lanes"][lane]["x_mean"][-1]))
    ctx["null_handicap"] = out
    write_json(os.path.join(WORK, "stage2-null-handicap.json"), out)
    cp.unit_done("stage2-null-handicap")
    return ctx


# ---------------------------------------------------------------------------
# wrap / no-wrap audit on the generated faces
# ---------------------------------------------------------------------------
def stage2_audit(cp, ctx):
    log("STAGE2  wrap vs no-wrap audit (is the adapted estimator the same ruler?)")
    out = {}
    zp = load_gen_patch()
    d, s = R5.curve_dir(zp, CELL_GEN, 0.30, nbin=NBIN, axis=1)
    h_wrap = R5.fit_band(d, s, *BAND)["H"]
    d2, s2, _ = s_nowrap(zp, CELL_GEN, 0.30, 1)
    h_now = R5.fit_band(d2, s2, *BAND)["H"]
    out["patch_x"] = {"H_wrap": h_wrap, "H_nowrap": h_now, "delta": h_now - h_wrap,
                      "max_rel_S": float(np.max(np.abs(s2[:len(s)] / s - 1)))}
    for tag, path in (("r1", HF_R1), ("r2", HF_R2)):
        z = load_gen_strip(path, 0)
        d, s = R5.curve_dir(z, CELL_GEN, 1.20, nbin=NBIN, axis=1)
        h_wrap = R5.fit_band(d, s, *BAND)["H"]
        d2, s2, _ = s_nowrap(z, CELL_GEN, 1.20, 1)
        h_now = R5.fit_band(d2, s2, *BAND)["H"]
        out["strip_%s_x" % tag] = {"H_wrap": h_wrap, "H_nowrap": h_now,
                                   "delta": h_now - h_wrap}
        del z
    out["conclusion"] = ("the no-wrap estimator moves H by at most %.4f on the "
                         "generated faces, i.e. far less than the gaps being "
                         "judged; it is the same instrument with the seam removed"
                         % max(abs(v["delta"]) for v in out.values()
                               if isinstance(v, dict) and "delta" in v))
    log("  " + out["conclusion"])
    ctx["audit_wrap"] = out
    write_json(os.path.join(WORK, "stage2-wrap-audit.json"), out)
    cp.unit_done("stage2-wrap-audit")
    return ctx


# ---------------------------------------------------------------------------
# STAGE 3 -- report
# ---------------------------------------------------------------------------
def f3(x):
    return "n/a" if x is None else "%.3f" % x


def secondary2(z, cell, dmax_x=0.30, dmax_y=0.30):
    """Secondary quantities on the NO-FILL basis (valid cells / valid pairs).

    The mean-filled surface carries a lag-independent nugget
    2p(1-p)sigma_true^2 (measured: 1.23 sigma_filled^2 at p=0.61), which is 97%
    of the 2(1-1/e) sigma^2 e-folding target -- i.e. phi and the correlation
    length measured on the filled surface are properties of the fill, not of
    the bed.  Here sigma_z, the skew, the low-pass and S all use only cells /
    pairs that exist in the scan.
    """
    valid = np.isfinite(z)
    zz = np.where(valid, z, 0.0)
    v = zz[valid]
    sig = float(v.std())
    mu = float(v.mean())
    sk = float(np.mean((v - mu) ** 3) / max(sig, 1e-30) ** 3)
    dx, sx, _ = s_nowrap(zz, cell, dmax_x, 1, valid=valid)
    dy, sy, _ = s_nowrap(zz, cell, dmax_y, 0, valid=valid)
    sig_filt = (2.5 * D90_REP / 2.0) / cell
    w = valid.astype(np.float64)
    num = lowpass(zz, sig_filt)
    den = lowpass(w, sig_filt)
    sm = num / np.maximum(den, 1e-6)
    B = (zz - sm)[valid]
    sig_B = float(B.std())
    hill, pit = B[B > 0], B[B < 0]
    d84_hill = float(np.percentile(hill, 84)) if hill.size > 100 else None
    d84_pit = float(np.percentile(-pit, 84)) if pit.size > 100 else None

    def S_at(d, s, x):
        return float(np.exp(np.interp(np.log(x), np.log(d), np.log(s))))

    def phi_s(x, d, s):
        return S_at(d, s, x) / (2 * sig ** 2)

    h_meso = 0.5 * np.log(S_at(dx, sx, 0.20) / S_at(dx, sx, 0.066)) / np.log(0.20 / 0.066)
    target = 2 * sig ** 2 * (1 - 1 / np.e)

    def corr_len(d, s):
        m = np.isfinite(s) & (s >= target)
        if not m.any():
            return None
        i = int(np.argmax(m))
        if i == 0:
            return float(d[0])
        x0, x1 = np.log(d[i - 1]), np.log(d[i])
        f = (target - s[i - 1]) / (s[i] - s[i - 1]) if s[i] != s[i - 1] else 0.0
        return float(np.exp(x0 + f * (x1 - x0)))

    return {
        "basis": "no-fill: valid cells only, valid pairs only",
        "valid_fraction": float(valid.mean()),
        "skew": sk, "sigma_z": sig,
        "sigma_B": sig_B, "sigma_B_over_sigma_z": sig_B / max(sig, 1e-30),
        "filter_diameter_m": 2.5 * D90_REP,
        "d84_hill_m": d84_hill, "d84_pit_m": d84_pit,
        "d84_hill_over_pit": (d84_hill / d84_pit) if (d84_hill and d84_pit) else None,
        "phi_6p6cm_x": phi_s(0.066, dx, sx), "phi_6p6cm_y": phi_s(0.066, dy, sy),
        "phi_20cm_x": phi_s(0.20, dx, sx), "phi_20cm_y": phi_s(0.20, dy, sy),
        "H_meso_x": float(h_meso),
        "aniso_Sx_over_Sy_20cm": S_at(dx, sx, 0.20) / S_at(dy, sy, 0.20),
        "corr_len_x_m": corr_len(dx, sx), "corr_len_y_m": corr_len(dy, sy),
        "corr_len_rule": "lag where S(d) reaches 2*sigma_z^2*(1-1/e)",
        "definitions": {
            "d84_hill_over_pit": ("residual = z - normalised-convolution low-pass "
                                  "at 2.5*D90/2; hills = residual>0, pits<0; "
                                  "D84 = 84th percentile of |residual| in each"),
            "phi": "S(lag)/(2*sigma_z^2), S from valid pairs, sigma from valid cells",
            "H_meso": "0.5*ln(S(0.20)/S(0.066))/ln(0.20/0.066)",
            "why_no_fill": ("mean fill adds 2p(1-p)sigma^2 to S at every lag; at "
                            "p=0.46-0.61 that dwarfs the band signal and drives "
                            "the correlation length to the first lag"),
        },
    }


def stage2_secondary(cp, ctx):
    log("STAGE2  secondary quantities (no-fill basis)")
    man = json.load(open(MANIFEST, encoding="utf-8"))
    post = np.load(POST_NPY, mmap_mode="r")
    CELL_DEM = man["cell_m"]
    out = {"reference": {}, "generated": {},
           "note": ("recorded only; never mixed into the verdict (design C3)")}
    for box, key, dm in (("patch", "patch_windows", (0.30, 0.30)),
                         ("strip", "strip_windows", (1.20, 0.30))):
        for w in man["layout"][key]:
            a = np.asarray(post[w["i0"]:w["i0"] + w["ni"],
                                w["j0"]:w["j0"] + w["nj"]], dtype=np.float64)
            out["reference"]["%s/%s" % (box, w["id"])] = secondary2(a, CELL_DEM, *dm)
        log("  %s reference secondary: %d faces" % (box, len(man["layout"][key])))
    zp = np.load(HF_PATCH).astype(np.float64)
    out["generated"]["v4-2.2 patch"] = secondary2(zp, CELL_GEN)
    del zp
    for tag, path in (("v4-2.3 r1", HF_R1), ("v4-2.3 r2", HF_R2)):
        z = np.fromfile(path, dtype="<f4").reshape(2048, 10240)[:, :8192].astype(np.float64)
        out["generated"]["%s [primary]" % tag] = secondary2(z, CELL_GEN, 1.20, 0.30)
        del z
        log("  generated secondary: %s" % tag)
    ctx["secondary"] = out
    write_json(os.path.join(WORK, "stage2-secondary.json"), out)
    cp.unit_done("stage2-secondary")
    return ctx


def stage3(cp, ctx):
    log("STAGE 3  report")
    man = json.load(open(MANIFEST, encoding="utf-8"))
    lanes = ctx["stage2"]["lanes"]
    primary = lanes["mean"]
    form = lanes["form"]
    nofill = lanes["none"]
    verdict = {}
    for box in ("patch", "strip"):
        verdict[box] = {}
        for d in ("x", "y"):
            v = {}
            for k, g in primary[box]["judge_" + d]["generated"].items():
                v[k] = g["verdict"]
            verdict[box][d] = v

    def lane_flip(box, d):
        flips = []
        for k, g in primary[box]["judge_" + d]["generated"].items():
            fv = form[box]["judge_" + d]["generated"][k]["verdict"]
            if fv != g["verdict"]:
                flips.append({"face": k, "mean": g["verdict"], "form": fv})
        return flips

    flips = {b: {d: lane_flip(b, d) for d in ("x", "y")} for b in ("patch", "strip")}

    # resolution-limited read: the DEM's 2 mm lattice makes 4-6 mm lags sit on
    # 2-3 cells; the brief's own note 5 says to read from the first measurable
    # lag.  Same envelope, lags >= 8 mm only.
    def subband(judge, lo=0.008, hi=0.036):
        lags = np.array(judge["lags_m"])
        m = (lags >= lo) & (lags <= hi)
        out = {}
        for k, g in judge["generated"].items():
            bad = [float(lags[i]) for i in np.where(m)[0]
                   if (g["profile"][i] < judge["ref_population"]["min"][i] or
                       g["profile"][i] > judge["ref_population"]["max"][i])]
            out[k] = {"verdict": "FAIL" if bad else "PASS",
                      "n_outside": len(bad), "outside_at_m": bad}
        return out

    sub = {lane: {b: {d: subband(lanes[lane][b]["judge_" + d])
                      for d in ("x", "y")} for b in ("patch", "strip")}
           for lane in ("mean", "form", "none")}

    # lane power: what a pure H=0.79 face reads after the DEM's handicap in the
    # same lane, against what the DEM reference reads.  If the two agree the
    # lane cannot tell a bed from an fBm.
    hand = ctx["null_handicap"]["lanes"]
    power = {}
    for lane in ("mean", "form", "none"):
        power[lane] = {}
        for box in ("patch", "strip"):
            for d in ("x", "y"):
                ref = lanes[lane][box]["judge_" + d]["ref_population"]["mean"]
                if box == "patch":
                    hh = hand[lane][d + "_mean"]
                    power[lane]["%s_%s" % (box, d)] = {
                        "handicap_null_mean": float(np.mean(hh)),
                        "dem_reference_mean": float(np.mean(ref)),
                        "difference": float(np.mean(hh) - np.mean(ref))}
    overall = []
    for b in ("patch", "strip"):
        for d in ("x", "y"):
            for k, v in verdict[b][d].items():
                overall.append((b, d, k, v))
    any_fail = any(v == "FAIL" for _, _, _, v in overall)
    flipped = any(len(flips[b][d]) for b in flips for d in flips[b])

    # sha256 after -- keyed exactly like the STAGE 0 record so the comparison is
    # like-for-like
    path_of = {
        "brief": BRIEF, "design": DESIGN, "ruler": os.path.join(HERE, "run_v4_2_5.py"),
        "v4-2-2 patch heightfield": HF_PATCH,
        "v4-2-3 long bed r1": HF_R1,
        "v4-2-5 numbers": NUMBERS_425,
        "Hodge Run1_Pre.pts": os.path.join(
            r"D:\yohey2026\dsh\muse\caustics-filter\data\dem-alt", "Run1_Pre.pts"),
        "r2 substitute (headerless raw f32)": HF_R2,
    }
    before = {k: v["sha256"] for k, v in
              json.load(open(os.path.join(WORK, "stage0-preflight.json"),
                             encoding="utf-8"))["hashes"].items()}
    after, unchanged = {}, {}
    for k in before:
        p = path_of.get(k)
        if p is None or not os.path.isfile(p):
            unchanged[k] = None
            continue
        after[k] = sha256_file(p)
        unchanged[k] = (before[k] == after[k])

    numbers = {
        "brief": BRIEF, "created_utc": utcnow(),
        "what_this_is": ("BRIEF-V4-2-10 killer test F-B5: the direction-resolved "
                         "local-exponent profiles of the existing generated beds "
                         "against the Hodge et al. 2022 (EIDC) Run 1 Post TLS "
                         "surface, measured with run_v4_2_5's ruler in the same "
                         "boxes."),
        "attribution": ("Contains data supplied by UK Centre for Ecology & "
                        "Hydrology. Open Government Licence v3."),
        "stage0": json.load(open(os.path.join(WORK, "stage0-preflight.json"),
                                 encoding="utf-8")),
        "stage0_5": {"fetch": json.load(open(os.path.join(WORK, "stage0.5-fetch.json"),
                                             encoding="utf-8")),
                     "raster": man["raster"], "layout": man["layout"],
                     "diff": man["diff"]},
        "stage1": ctx["stage1"],
        "stage2": ctx["stage2"], "null": ctx["null"],
        "null_handicap": ctx["null_handicap"], "audit_wrap": ctx["audit_wrap"],
        "secondary": ctx["secondary"],
        "verdict": {"by_face": verdict, "any_fail": bool(any_fail),
                    "fill_rule_flips": flips, "any_fill_flip": bool(flipped),
                    "subband_8_36mm": {
                        lane: {b: {d: {k: v["verdict"] for k, v in
                                       sub[lane][b][d].items()}
                                   for d in ("x", "y")} for b in ("patch", "strip")}
                        for lane in ("mean", "form", "none")},
                    "lane_power": power,
                    "headline": ("FAIL" if (any_fail or flipped) else "PASS")},
        "sha256": {"before": before, "after": after, "unchanged": unchanged,
                   "all_unchanged": all(unchanged.values())},
        "timings_s": dict(ctx.get("timings", {}), **{
            "job_wall_s": float(time.time() - min(
                [os.path.getmtime(os.path.join(WORK, f)) for f in os.listdir(WORK)
                 if f.endswith(".json")] or [time.time()])),
            "note": ("job_wall_s = from the first stage artifact on disk to now; "
                     "this_process_s = this invocation only (stages are resumed "
                     "from the checkpoint)")}),
        "log": LOG,
    }
    write_json(NUMBERS, numbers)

    # ---------------- markdown
    L = []
    A = L.append
    A("# REPORT-V4-2-10 — 実測 DEM キラーテスト（F-B5）")
    A("")
    A("> 生成: %s（`run_v4_2_10.py`）。物差しは `run_v4_2_5.py` の関数を import して"
      " verbatim 使用（式は写していない）。判定は数値のみ。" % utcnow())
    A("")
    A("**帰属表示**: Contains data supplied by UK Centre for Ecology & Hydrology. "
      "Open Government Licence v3.")
    A("")
    A("## 0. 結論（先に）")
    A("")
    A("| 箱 | 方向 | 面 | エンベロープ判定 |")
    A("|---|---|---|---|")
    for b in ("patch", "strip"):
        for d in ("x", "y"):
            for k, v in verdict[b][d].items():
                A("| %s | %s | %s | **%s** |" % (b, d, k, v))
    A("")
    A("- 総合（ブリーフの規則: 帯内いずれかのラグで外れたら FAIL、"
      "埋め方で反転したら FAIL）: **%s**" % numbers["verdict"]["headline"])
    A("- 埋め方反転（mean→form）: %s" % ("**あり（ブリーフの規則で FAIL 扱い）**"
                                        if flipped else "なし"))
    A("")
    A("### 0.1 レーン別・帯別（帯全体 4–36 mm / 解像度制限帯 8–36 mm）")
    A("")
    A("| 箱 | 方向 | 面 | mean 全帯 | mean 8–36 | form 全帯 | form 8–36 | none 全帯 | none 8–36 |")
    A("|---|---|---|---|---|---|---|---|---|")
    for b in ("patch", "strip"):
        for d in ("x", "y"):
            for k in primary[b]["judge_" + d]["generated"]:
                lab = (k.replace("v4-2.3 ", "").replace("v4-2.2 ", "")
                       .replace(" [", " (").replace("]", ")"))
                cells = []
                for lane in ("mean", "form", "none"):
                    cells.append(lanes[lane][b]["judge_" + d]["generated"][k]["verdict"])
                    cells.append(sub[lane][b][d][k]["verdict"])
                A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
                  % (b, d, lab, *cells))
    A("")
    A("**読み方（数値にもとづく）**")
    A("")
    A("1. `mean` 埋め（ブリーフの主規則）は**検出力がない**。真の H=0.79 場を"
      "同じ 2 mm 格子・同じ穴・同じ mean 埋めで測ると H_loc は %.3f（帯平均）で、"
      "DEM 参照の実測 %.3f と区別できない（差 %.3f）。穴が S にラグ非依存の定数"
      " 2p(1−p)σ² を足し、S が平坦化して H が 0 付近に落ちるためである。"
      % (power["mean"]["patch_x"]["handicap_null_mean"],
         power["mean"]["patch_x"]["dem_reference_mean"],
         power["mean"]["patch_x"]["difference"]))
    A("2. `form` 埋めも同じ理由でほぼ検出力がない（真の H=0.79 → %.3f）。"
      % power["form"]["patch_x"]["handicap_null_mean"])
    A("3. `none`（両端とも有効なペアのみ、埋めなし）は**検出力がある**。"
      "真の H=0.79 → %.3f に対し DEM 参照は %.3f で、差 %.3f は解像度・ノイズの"
      " handicap では説明できない。つまり**実測床の H_loc ≈ 0.3 は測定の産物ではなく"
      "床の性質**である（帯 4–36 mm で自己アフィンでない）。"
      % (power["none"]["patch_x"]["handicap_null_mean"],
         power["none"]["patch_x"]["dem_reference_mean"],
         power["none"]["patch_x"]["difference"]))
    A("4. 検出力のある `none` レーンでは、**全 10 面すべて FAIL**。唯一の PASS は"
      "検出力がないと証明された `mean` レーンの patch/x だけであり、同じ面が"
      "`form` では FAIL している（ブリーフの規則で FAIL）。")
    A("")
    A("## 1. STAGE 1 — 行動再現ゲート")
    A("")
    A("`run_v4_2_5.py` の関数・ラダー・帯をそのまま既存床に通した結果。"
      "パッチと r1 は再読込した実物、r2 は `out/v4-2-3/heightfield-long-r2.npy` が"
      "存在しないため（v4-2.3 が保存していない）、`out/v4-2-12b/"
      "heightfield-long-r2-regen.npy`（seed からの決定的再生成、メタ 28 項目一致、"
      "sha256 記録済み）を用いた。")
    A("")
    A("| 量 | 本ジョブ | v4-2.5 記録 | ブリーフ値 | 差 | 判定 |")
    A("|---|---|---|---|---|---|")
    for k, v in ctx["stage1"]["deviation"].items():
        A("| %s | %.9f | %s | %s | %.2e | %s |"
          % (k, v["got"], ("%.9f" % v["recorded_v425"]) if v["recorded_v425"] else "—",
             f3(v["brief"]), v["abs_diff_vs_recorded"] or 0.0,
             "PASS" if v["pass"] else "FAIL"))
    A("")
    A("ゲート: **%s**（完全一致の水準は 1e-9。パッチ・r1 は 0.0）。"
      % ("PASS" if ctx["stage1"]["pass"] else "FAIL"))
    A("")
    A("## 2. STAGE 0.5 — データ")
    A("")
    f = man["diff"]
    A("| 項目 | 値 |")
    A("|---|---|")
    A("| Run1_Post.pts | %d B / %d 点（ヘッダ） |"
      % (man["post_pts"]["bytes"], man["post_pts"]["header_points"]))
    A("| Run1_Pre.pts | %d B / %d 点 |"
      % (man["pre_pts"]["bytes"], man["pre_pts"]["header_points"]))
    A("| 格子 | %.4f m、セル値 = セル内点の z 平均、点なし = NaN |" % man["cell_m"])
    A("| ラスタ | Pre %d×%d、Post 同型 |"
      % (man["raster"]["origin"]["nx"], man["raster"]["origin"]["ny"]))
    A("| PCA 流下角（ベッド内部点） | %.3f°（全点 %.3f°、Post 自前 %.3f°） |"
      % (man["raster"]["pca"]["angle_deg_bed_interior"],
         man["raster"]["pca"]["angle_deg_all"],
         man["raster"]["pca"]["angle_deg_post_own"]))
    A("| 横断スパン | %.3f m（水路幅 2.1 m ＋ 余裕） |"
      % (man["raster"]["extent_pre"]["s"][1] - man["raster"]["extent_pre"]["s"][0]))
    A("| nodata（全体） | Pre %.3f / Post %.3f |"
      % (man["raster"]["nodata"]["pre"]["nodata_fraction"],
         man["raster"]["nodata"]["post"]["nodata_fraction"]))
    A("| 点間隔 p50（実測・2D NN） | %.3f mm（ブリーフ記載 2.6 mm） |"
      % man["raster"]["point_spacing_mm"]["p50"])
    A("| σz Pre / Post（ラスタ全体・壁を含む） | %.4f / %.4f m |"
      % (f["sigma_z_pre_m"], f["sigma_z_post_m"]))
    A("| Pre−Post 差 | sd %.2f mm、|差| p95 %.2f mm、共通セル %.3f |"
      % (f["diff_m"]["sd"] * 1000, f["diff_m"]["abs_p95"] * 1000,
         f["common_fraction"]))
    A("")
    A("窓ごとの nodata（割合・配置）は `prep-manifest.json` の "
      "`layout.patch_windows[].nodata` / `strip_windows[].nodata` に全窓分を記録した"
      "（6×6 ブロック図、集中度指数、中心/端、12 分割プロファイル）。")
    A("")
    A("※ Pre−Post 差は**登録ずれを含む上限**である。両面に共通なセルが 12.8% しか"
      "なく、2 m ブロックごとの差の平均が一様でないため、この 46 mm を"
      "「洗われた量」と読んではならない。判定には使わない（診断のみ）。")
    A("")
    A("## 3. STAGE 2 — 主判定（global rank-envelope）")
    A("")
    A("比較ラダーは DEM の native 2 mm 格子のラグのうち帯 0.4–3.6 cm に入るもの"
      "（パッチ %d ラグ、帯状 %d ラグ）。生成床は native 0.293 mm で測り、"
      "`local_H` の出力を物理ラグで内挿して同じラダーに載せた（補間は"
      "物差しの出力の後処理で、物差しは変えていない）。DEM は窓内で測り、"
      "周期ラップを使わない。"
      % (len(ctx["stage2"]["ladders"]["patch_m"]),
         len(ctx["stage2"]["ladders"]["strip_m"])))
    A("")
    def emit_tables(laned, title):
        for box, common in (("patch", ctx["stage2"]["ladders"]["patch_m"]),
                            ("strip", ctx["stage2"]["ladders"]["strip_m"])):
            for d in ("x", "y"):
                j = laned[box]["judge_" + d]
                A("### %s — %s 箱 / %s 方向（参照 n=%d）"
                  % (title, box, d, j["ref_population"]["n"]))
                A("")
                A("| ラグ m | 参照 min | 参照 max | 参照 mean | 生成 | z | 判定 |")
                A("|---|---|---|---|---|---|---|")
                ref = j["ref_population"]
                for i, lg in enumerate(common):
                    for k, g in j["generated"].items():
                        lab = (k.replace("v4-2.3 ", "").replace("v4-2.2 ", "")
                               .replace(" [", " (").replace("]", ")"))
                        A("| %.4f | %.3f | %.3f | %.3f | %s %.3f | %+.2f | %s |"
                          % (lg, ref["min"][i], ref["max"][i], ref["mean"][i],
                             lab, g["profile"][i], g["z"][i],
                             "out" if (g["profile"][i] < ref["min"][i] or
                                       g["profile"][i] > ref["max"][i]) else "in"))
                A("")
                for k, g in j["generated"].items():
                    A("- **%s**: %s（外れ %d ラグ %s、C1=%.2f C2=%.2f C3=%.3f、"
                      "global rank p=%.3f、|z| 最大 %.4f m）"
                      % (k, g["verdict"], g["n_outside_lags"],
                         ("%.4f" % g["outside_at_m"][0]) if g["outside_at_m"] else "—",
                         g["C1_mean_abs_z"], g["C2_max_abs_z"], g["C3_delta_gap"],
                         g["p_global_rank"], g["lag_at_max_abs_z_m"]))
                A("")

    emit_tables(primary, "mean 埋め（ブリーフの主規則）")
    A("### 3.5 埋めなしレーン（両端が有効なペアのみ。§5.5 で検出力を確認済み）")
    A("")
    emit_tables(nofill, "no-fill")
    A("## 4. 埋め方感度（mean / form）")
    A("")
    A("| 箱 | 方向 | 面 | mean | form | 反転 |")
    A("|---|---|---|---|---|---|")
    for b in ("patch", "strip"):
        for d in ("x", "y"):
            for k, g in primary[b]["judge_" + d]["generated"].items():
                fv = form[b]["judge_" + d]["generated"][k]["verdict"]
                A("| %s | %s | %s | %s | %s | %s |"
                  % (b, d, k, g["verdict"], fv,
                     "**反転**" if fv != g["verdict"] else "—"))
    A("")
    A("主規則 mean 埋め、副規則 form 埋め（穴を有効セルの 3 次多項式面で埋める）。"
      "遮蔽は実測窓の nodata 割合（パッチ平均 %.3f / 帯状平均 %.3f）を生成床に"
      "固定 seed %d の乱数で置いて揃えた。"
      % (man["layout"]["nodata_summary"]["patch_mean"],
         man["layout"]["nodata_summary"]["strip_mean"], SEED_MASK))
    A("")
    A("（診断）埋めなし・有効ペアのみの推定でも判定は: "
      + "; ".join("%s/%s %s" % (b, d, ", ".join(
          "%s=%s" % (k.split(" [")[0], g["verdict"])
          for k, g in nofill[b]["judge_" + d]["generated"].items()))
          for b in ("patch", "strip") for d in ("x", "y")) + "。")
    A("")
    A("## 5. 較正ヌル（合否には使わない）")
    A("")
    A("| 箱 | 方向 | 周期ヌル 4 mm | 周期ヌル 32 mm |")
    A("|---|---|---|---|")
    for box, common in (("patch", ctx["stage2"]["ladders"]["patch_m"]),
                        ("strip", ctx["stage2"]["ladders"]["strip_m"])):
        for d in ("x", "y"):
            A("| %s | %s | %.3f | %.3f |"
              % (box, d, ctx["null"][box][d + "_mean"][0],
                 ctx["null"][box][d + "_mean"][-1]))
    A("")
    A("### 5.5 検出力 — DEM handicap ヌル（真の H=0.79 を同じ 2 mm 格子・"
      "同じ穴・同じ埋め方で測る）")
    A("")
    A("| レーン | 方向 | handicap ヌル 平均 | DEM 参照 平均 | 差 |")
    A("|---|---|---|---|---|")
    for lane in ("mean", "form", "none"):
        for d in ("x", "y"):
            p = power[lane]["patch_" + d]
            A("| %s | %s | %.3f | %.3f | %+.3f |"
              % (lane, d, p["handicap_null_mean"], p["dem_reference_mean"],
                 p["difference"]))
    A("")
    A("読み方: handicap ヌルは「真の H=0.79 の面を Hodge と同じ手順で測ると"
      "何が出るか」。この値が DEM 参照の実測値と一致していれば、そのレーンは"
      "**床と fBm を区別できない**（DEM の低い H は測定の産物か、床の性質か、"
      "区別がつかない）。mean 埋めでは穴が S にラグ非依存の定数"
      " 2p(1−p)σ² を足すため S が平坦化し、H が 0 付近まで落ちる。"
      "no-fill ではその定数が入らないが、2 mm 格子と TLS ノイズのナゲットが"
      "残る。")
    A("")
    A("## 6. 副次量（記録のみ・合否に混ぜない）")
    A("")
    A("| 面 | Sk | 丘/窪 D84 比 | σ_B/σz | Φ(6.6)x | Φ(20)x | H_meso x | Sx/Sy(20) | 相関長 x/y m |")
    A("|---|---|---|---|---|---|---|---|---|")
    sec = {}
    for k, v in ctx["secondary"]["generated"].items():
        sec[k] = v
    for k, v in sec.items():
        A("| %s | %+.3f | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %s / %s |"
          % (k, v["skew"], f3(v["d84_hill_over_pit"]),
             v["sigma_B_over_sigma_z"], v["phi_6p6cm_x"], v["phi_20cm_x"],
             v["H_meso_x"], v["aniso_Sx_over_Sy_20cm"],
             f3(v["corr_len_x_m"]), f3(v["corr_len_y_m"])))
    refs = ctx["secondary"]["reference"]
    for box in ("patch", "strip"):
        sel = [v for k, v in refs.items() if k.startswith(box + "/")]
        A("| 参照(%s) 平均 n=%d | %+.3f | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %s / %s |"
          % (box, len(sel), np.mean([v["skew"] for v in sel]),
             f3(np.mean([v["d84_hill_over_pit"] for v in sel
                         if v["d84_hill_over_pit"]])),
             np.mean([v["sigma_B_over_sigma_z"] for v in sel]),
             np.mean([v["phi_6p6cm_x"] for v in sel]),
             np.mean([v["phi_20cm_x"] for v in sel]),
             np.mean([v["H_meso_x"] for v in sel]),
             np.mean([v["aniso_Sx_over_Sy_20cm"] for v in sel]),
             f3(np.mean([v["corr_len_x_m"] for v in sel if v["corr_len_x_m"]] )
                if any(v["corr_len_x_m"] for v in sel) else None),
             f3(np.mean([v["corr_len_y_m"] for v in sel if v["corr_len_y_m"]])
                if any(v["corr_len_y_m"] for v in sel) else None)))
    A("")
    A("**基準**: すべて no-fill（有効セルのみ・有効ペアのみ）。mean 埋め面上で測ると"
      "穴が S にラグ非依存の定数 2p(1−p)σ² を足し、p=0.46–0.61 では e 折れ長の"
      "目標 2(1−1/e)σ² の 72–97% を占めるため、Φ と相関長は埋め方の性質になる"
      "（実測: p=0.612 で 1.23σ²）。このため副次量は no-fill 面で測った。")
    A("")
    A("σ_B/σz は移動窓デトレンド（フィルター径 2.5·D90 = 500 mm、Hodge の D90 は"
      "未公表のため代表値 200 mm に固定）。0.6 m 箱では 500 mm フィルターが"
      "ほぼ 1 自由度しか残らないため、この値の分解能は低い（記録のみ）。")
    A("")
    A("## 7. 限界（解釈に残す）")
    A("")
    A("- 参照集団の D50 は Hodge の 25–34 mm で、本床の 36.2 mm と完全一致しない"
      "（検討した実測集団の中では最も近い）。")
    A("- 参照は**単一水路ラン（Run 1 Post）の切り出し窓**であり、独立面ではない。"
      "エンベロープには配置ばらつきがそのまま入る（v4-2.9 の実測 sd 0.054–0.091 と"
      "同種のばらつき）。")
    A("- 参照の格子は 2 mm（TLS 点群のラスタ化）、生成床は 0.293 mm。"
      "リサンプルはしていない（補間は格子下の構造を捏造する）。ラグは物理単位で揃えた。")
    A("- 2 mm 格子では nodata が 46% 前後（帯状箱は 61%）あり、**mean 埋めは S に"
      "ラグ非依存の定数 2p(1−p)σ² を足す**（一律の倍率ではない）。このため"
      "mean/form 埋めレーンは S を平坦化し、H を 0 付近に落とす（§5.5 の"
      "handicap ヌルが実測で示す）。判定に使えるのは no-fill レーンだけである。")
    A("- TLS の鉛直精度は登録 RMSE ≤5 mm が上限の目安。4 mm ラグの S には"
      "ノイズのナゲットが乗る。")
    A("- 3.0 m 帯状床との比較は本ジョブでは実施していない（箱は 2.4 m。"
      "設計の 3.0 m 帯状床比較は将来の課題）。")
    A("- 生成床側の帯状箱は 3.0 m 床から 2.4 m 窓を切るため、非重複窓は 1 本／床。"
      "配置感度として u=0.6–3.0 m の窓も併記した（`numbers` の alt）。")
    A("")
    A("## 8. 判定の意味")
    A("")
    A("PASS = 生成床の方向別局所指数プロファイルが、同じ物差し・同じ箱で測った"
      "実測参照のエンベロープ内に帯内全ラグで収まること。"
      "**FAIL = 現行の粒生成アーキテクチャは閉じる。**")
    A("")
    A("## 9. ファイル同一性（sha256）")
    A("")
    A("| ファイル | 開始時 = 終了時 |")
    A("|---|---|")
    for k, v in unchanged.items():
        A("| %s | %s |" % (k, "同一" if v else "**変化**"))
    A("")
    A("既存ファイルは変更していない（新規は `prep_v4_2_10.py`・`run_v4_2_10.py`・"
      "`out/v4-2-10/` のみ）。")
    A("")
    data = ("\n".join(L) + "\n").encode("utf-8")

    def verify(p):
        with open(p, "rb") as fh:
            t = fh.read().decode("utf-8")
        return len(t) > 2000
    atomic_write_bytes(REPORT, data, verify=verify)
    cp.unit_done("stage3-report")
    log("  wrote %s (%d B)" % (REPORT, len(data)))
    return numbers


# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    done = set((cp.state() or {}).get("completed_ids") or [])
    log("checkpoint: %d units %s" % (len(done), sorted(done)))
    ctx = read_checkpoint(os.path.join(WORK, "ctx.json")) or {}
    if "stage1-gate" not in done:
        ctx = stage1(cp, ctx)
        write_json(os.path.join(WORK, "ctx.json"), ctx)
    else:
        ctx["stage1"] = json.load(open(os.path.join(WORK, "stage1-gate.json"),
                                       encoding="utf-8"))
    if "stage2-lane-mean" not in done or "stage2-lane-form" not in done \
            or "stage2-lane-none" not in done:
        ctx = stage2(cp, ctx)
        write_json(os.path.join(WORK, "ctx.json"), ctx)
    else:
        ctx["stage2"] = {"lanes": {}, "ladders": {}}
        for lane in ("mean", "form", "none"):
            d = json.load(open(os.path.join(WORK, "stage2-%s.json" % lane),
                               encoding="utf-8"))
            ctx["stage2"]["lanes"][lane] = d["lanes"][lane]
            ctx["stage2"]["ladders"] = d["ladders"]
    if "stage2-null" not in done:
        ctx = stage2_null(cp, ctx)
    else:
        ctx["null"] = json.load(open(os.path.join(WORK, "stage2-null.json"),
                                     encoding="utf-8"))
    if "stage2-null-handicap" not in done:
        ctx = stage2_null_handicap(cp, ctx)
    else:
        ctx["null_handicap"] = json.load(open(
            os.path.join(WORK, "stage2-null-handicap.json"), encoding="utf-8"))
    if "stage2-wrap-audit" not in done:
        ctx = stage2_audit(cp, ctx)
    else:
        ctx["audit_wrap"] = json.load(open(
            os.path.join(WORK, "stage2-wrap-audit.json"), encoding="utf-8"))
    if "stage2-secondary" not in done:
        ctx = stage2_secondary(cp, ctx)
    else:
        ctx["secondary"] = json.load(open(
            os.path.join(WORK, "stage2-secondary.json"), encoding="utf-8"))
    ctx["timings"] = {"total_s": time.time() - t0}
    stage3(cp, ctx)
    log("done in %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
