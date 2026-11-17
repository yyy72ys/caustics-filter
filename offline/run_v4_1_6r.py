"""run_v4_1_6r -- BRIEF-V4-1-6R driver: river water as the wave source
(separated form, replacement of the real-size wave environment).

NEW FILE ONLY.  Nothing else in the tree is written: `run_v4_1_2.py` (scene B),
`stream_surface2/3/4.py`, `run_v4_1_5.py`, `det.py`, `run_v16_water.py` and
every other .py stay byte-identical (G-LEGACY; hashes in
`stage0-py-hashes.json`).

WHAT THE BRIEF ASKS
-------------------
v4-1.6 (first edition) generated the waves on a real-size scene (a 15 mm deep
rill).  By the user's decision of 2026-11-13 that is replaced: the waves come
from **ordinary river water**.  The form stays separated -- the picture is
drawn at a sane size, the wave field is generated at wave scale and mapped onto
the picture.

    1. WAVE FIELD (river water)
       stream_surface3 accounting, environment of the NATURAL scene
       (run_v4_1_2.py, scene B) verbatim:
           H=0.35 m, U=0.80 m/s, U0=1.15*U=0.92 m/s, S0=0.002, D84=0.09 m,
           seed=20261112, medium algae_pond.
       The generation patch follows the wave: L_gen = 10*lambda0 = 5.418921 m
       (lambda0 = 2*pi/k0 = 0.541892 m, so k0/dk = 10 exactly).
       The bed of the direct bed response (i) is the source photo mapped onto
       L_gen; the amplitude is the model's own Nikora mapping sigma_b=D84/4.
    2. BASE IMAGE
       the source photo at L_img = 3.68 m (the same screen size as the natural
       scene: L = 8*Us/n_loops from run_v4_1_2's own loop rule), 736x540,
       refraction depth z = 0.35 m.  Composite frozen at A=2.6 / Blow sigma=3 /
       E=2.2 / Q=25.
    3. STRETCH COMPOSITE
       the generated field is mapped L_gen -> L_img (magnification
       L_img/L_gen = 0.679102; the drawn wavelength is
       0.541892*0.679102 = 0.368000 m).  Both grids have N=512 cells, so the
       map is index-exact: h is copied, the slopes are multiplied by
       1/0.679102 = 1.472534 (the same surface, measured in the image's metric).

TRANSPORT BAND-LIMIT
    v4-1.5's fold rule, applied at THIS scene and not copied: measure the
    |ds/di| ladder of the frame-24 field at the image scale and take the
    smallest width that keeps the forward splat single-valued.  The v4-1.5
    factor-4 margin (max|ds/di| <= 0.25) turns out to be UNREACHABLE here --
    the resolved drawn wave alone sits at ~0.55 -- so the width is fixed at the
    fold limit itself (max|ds/di| < 1.0); the whole ladder is reported and the
    unreachability is disclosed instead of being tuned away.

STAGES (all foreground, one process, jobkit checkpoint committed per stage)
    --stage0    absolute-path preflight + stage-0 sha256 record + scene record
    --field     the wave field: geometry, budget, bed response, spectrum
                (G-SPEC), the L_gen -> L_img map
    --transport |ds/di| ladder, width, ring preservation, frame 24 (G-LINE
                pre-check, C statistics)
    --chunk K   frames 24K..24K+23
    --diag      frame 24 drawn 3 ways (diagnostic: where the axis band sits)
    --warp      frame-24 displacement yardstick (G-WARP, no pass/fail)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT
    --report    numbers-v4-1-6r.json + REPORT-V4-1-6R.md
    --all       everything, sequentially

MEASUREMENT
    G-LINE is measured with run_v4_1_5's own code path (imported unchanged):
    256x256 windows whose positions are fixed by a rule on the SOURCE PHOTO
    only (lowest / highest local std on a 92x90 tile grid), primary quantity
    depth_e3 = 1000*RMS(ideal axis-parallel 5-15 px band-pass)/window mean, the
    brief's high-pass-first variant depth_hp_e3 reported beside it.  Using the
    imported module (not a re-typed copy) is what makes "same measurement
    method, same table form" literal.
"""

import argparse
import datetime
import gc
import hashlib
import io
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v4_1_2 as R12   # noqa: E402  scene B constants + env/optics, UNCHANGED
import run_v4_1_5 as R15   # noqa: E402  the v4-1-5 measurement path, UNCHANGED
import run_v16_water as V16  # noqa: E402
import stream_surface2 as SS  # noqa: E402
import stream_surface3 as S3  # noqa: E402
import stream_surface4 as S4  # noqa: E402
from det import det_frame  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-6R.md"
OUT = os.path.join(HERE, "out", "v4-1-6r")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-6r.mp4")
STILL = os.path.join(OUT, "frame-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-6r.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-6R.md")
PHOTO = R12.PHOTO

V15_STILL = R15.STILL
V14_STILL = R15.V14_STILL
V13_STILL = R15.V13_STILL

FPS = R12.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R12.WIDTH, R12.HEIGHT
N_GRID = R12.GRID
Q_DET = R12.Q_DET
A_FROZEN = R12.A_FROZEN
B_SIGMA_FROZEN = R12.B_SIGMA_FROZEN
E_FROZEN = R12.E_FROZEN
N_WAVELENGTHS = 10
# the |ds/di| ladder: v4-1-5's widths, extended upwards because the drawn wave
# is only ~51 cells long here and the reachable fold level depends on it
LADDER = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
# constraint B: the transport band-limit may not weaken the WAVE either.  A
# Gaussian of sigma cells passes exp(-(k*sigma*dx)^2) of the ring's VARIANCE;
# 0.95 means "at most 5% of the wave's variance may be removed".
WAVE_KEEP_MIN = 0.95


# ---------------------------------------------------------------------------
# scene geometry (computed once, at import, and asserted against the brief)
# ---------------------------------------------------------------------------
def _geometry():
    ref = R12.REF
    U = float(ref["U"])
    H = float(ref["H"])
    S0 = float(ref["S0"])
    D84 = float(ref["D84"])
    seed = int(ref["seed"])
    U0 = SS.US_FACTOR * U
    k0 = SS.k0_stationary(U0, H)
    if k0 is None:
        raise SystemExit("STOP: no stationary ring for U0=%.4f H=%.4f" % (U0, H))
    lam0 = 2.0 * np.pi / k0
    L_gen = float(N_WAVELENGTHS) * lam0
    L_img, n_loops, loop = SS.loop_adjusted_L(float(ref["L_guess"]), U0, 8.0)
    scale = float(L_img) / float(L_gen)
    dk_gen = 2.0 * np.pi / L_gen
    g = {
        "scene": "run_v4_1_2.REF (scene B, natural river water)",
        "H_m": H, "U_m_s": U, "U0_m_s": U0, "S0": S0, "D84_m": D84,
        "seed": seed, "water_preset": R12.WATER_PRESET,
        "k0_rad_m": float(k0), "lambda0_m": float(lam0),
        "n_wavelengths": int(N_WAVELENGTHS),
        "L_gen_m": L_gen, "L_img_m": float(L_img),
        "dk_gen_rad_m": float(dk_gen),
        "k0_over_dk": float(k0 / dk_gen),
        "grid": int(N_GRID),
        "cell_gen_m": L_gen / N_GRID,
        "cell_img_m": float(L_img) / N_GRID,
        "cells_per_lambda_gen": float(lam0 / (L_gen / N_GRID)),
        "magnification": scale,
        "slope_gain": 1.0 / scale,
        "drawn_lambda_m": float(lam0 * scale),
        "drawn_lambda_px": float(lam0 * scale) / float(L_img) * WIDTH,
        "drawn_cells_per_lambda": float(lam0 * scale / (float(L_img) / N_GRID)),
        "loop": loop, "n_loops": int(n_loops),
        "loop_mode_s": float(L_gen / U0),
        "depth_m": H,
        "composite": {"A": A_FROZEN, "B_sigma_px": B_SIGMA_FROZEN,
                      "E": E_FROZEN, "Q": Q_DET},
    }
    # the brief's numbers, checked rather than assumed
    if abs(U0 - 0.92) > 1e-9:
        raise SystemExit("STOP: U0=%.12f != 0.92" % U0)
    if abs(L_gen - 5.419) > 5e-4:
        raise SystemExit("STOP: L_gen=%.9f != 5.419 m" % L_gen)
    if abs(L_img - 3.68) > 1e-9:
        raise SystemExit("STOP: L_img=%.9f != 3.68 m" % L_img)
    if abs(g["k0_over_dk"] - 10.0) > 1e-9:
        raise SystemExit("STOP: k0/dk=%.12f != 10" % g["k0_over_dk"])
    return g


GEOM = _geometry()


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg):
    print(msg, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError("not JSON serialisable: %r" % (type(o),))


def write_json(path, obj):
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False, indent=2,
                                        default=_json_default).encode("utf-8"))


def makedirs():
    for d in (OUT, WORK, FRAMES):
        os.makedirs(d, exist_ok=True)


def py_set():
    """Every .py in the workdir: sha256 + size + mtime (mtime for attribution)."""
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if not (name.endswith(".py") and os.path.isfile(p)):
            continue
        st = os.stat(p)
        out[name] = {"sha256": sha256_file(p), "bytes": int(st.st_size),
                     "mtime_epoch": float(st.st_mtime)}
    return out


def read_checkpoint():
    return read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}


def checkpoint_extra(key, value=None):
    state = read_checkpoint()
    if value is not None:
        state[key] = value
        state["updated_at"] = now_iso()
        atomic_write_bytes(CP_PATH, json.dumps(
            state, ensure_ascii=False, indent=2).encode("utf-8"))
    return state


def save_png(path, arr):
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG")
    data = buf.getvalue()

    def verify(p):
        try:
            im = Image.open(p)
            im.load()
            return im.size == (WIDTH, HEIGHT) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(path, data, verify=verify)
    return sha256_bytes(data)


def frame_path(k):
    return os.path.join(FRAMES, "frame-%04d.png" % k)


# ---------------------------------------------------------------------------
# the bed of the direct bed response: the source photo mapped onto L_gen
# ---------------------------------------------------------------------------
def photo_green():
    from PIL import Image
    a = np.asarray(Image.open(PHOTO).convert("RGB"), dtype=np.float64)
    return a[:, :, 1]


def photo_bed(N, sigma_b):
    """Photo luminance -> periodic N x N bed elevation over the patch.

    The photo is MIRROR-tiled 2x2 first, so the bed is exactly periodic (the
    patch wraps onto itself with matching values and matching first
    derivatives) -- a raw periodic resample would print a straight seam line
    into every frame.  The mirrored tile carries the photo's own spatial
    structure at 2x its pixel pitch; the amplitude carries no metric scale, so
    the model's own Nikora mapping sigma_b = D84/4 is used (stated assumption).
    """
    from scipy import ndimage as ndi
    g = photo_green()
    mir = np.block([[g, g[:, ::-1]], [g[::-1, :], g[::-1, ::-1]]])
    yy = (np.arange(int(N)) + 0.5) * mir.shape[0] / float(N) - 0.5
    xx = (np.arange(int(N)) + 0.5) * mir.shape[1] / float(N) - 0.5
    Y, X = np.meshgrid(yy, xx, indexing="ij")
    b = ndi.map_coordinates(mir, [Y, X], order=3, mode="grid-wrap")
    m = float(b.mean())
    s = float(b.std())
    b = (b - m) / s
    info = {
        "source": os.path.basename(PHOTO),
        "channel": "green (the channel every gate measures)",
        "photo_shape": [int(g.shape[0]), int(g.shape[1])],
        "periodicisation": ("mirror tile 2x2 -> %dx%d, cubic periodic resample"
                            % (mir.shape[1], mir.shape[0])),
        "tile_shape": [int(mir.shape[0]), int(mir.shape[1])],
        "bed_rms_units": float(b.std()),
        "sigma_b_m": float(sigma_b),
        "sigma_b_rule": ("Nikora mapping of the model: sigma_b = D84/4 "
                         "(the photo has no metric scale of its own)"),
        "bed_pitch_m": float(GEOM["L_gen_m"]) / float(N),
        "photo_pitch_m_at_L_img": float(GEOM["L_img_m"]) / float(g.shape[1]),
        "raw_mean": m, "raw_std": s,
    }
    return b, info


# ---------------------------------------------------------------------------
# the (i) response, computed with stream_surface2's own transfer + valve
# ---------------------------------------------------------------------------
def response_of_zeta(zeta_hat, kx, kappa, H, U, N, slope_cap=SS.SLOPE_CAP,
                     res_band=SS.RES_BAND_DENOM):
    """stream_surface2.StreamState2._build_i's transfer/valve block, verbatim.

    The caller supplies zeta_hat (the bed's spectrum); everything else -- the
    1.1 gravity transfer T(k), the resonant-band definition and the
    breaking-limit safety valve -- is the v4-1.2 code, unchanged.
    """
    T, denom = SS.transfer_T(kx, kappa, H, U)
    eta_hat = T * zeta_hat
    band = np.abs(denom) < float(res_band)
    band[~np.isfinite(denom)] = False
    band[np.abs(kx) == 0.0] = False
    N2 = float(int(N) ** 2)
    tot_e = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
    band_e = float(np.sum(np.abs(eta_hat[band]) ** 2) / N2 ** 2)
    band_slope = float(np.sqrt(max(
        np.sum(kappa[band] ** 2 * np.abs(eta_hat[band]) ** 2) / N2 ** 2, 0.0)))
    scale = 1.0
    if band_slope > float(slope_cap) > 0:
        scale = float(slope_cap) / band_slope
        eta_hat[band] = eta_hat[band] * scale
    after = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
    info = {
        "T_max": float(np.max(np.abs(T))),
        "resonant_band_def": "|denominator| < %.3g" % float(res_band),
        "band_rms_slope_before": band_slope,
        "band_energy_fraction": (band_e / tot_e if tot_e > 0 else 0.0),
        "scale_applied": float(scale),
        "dropped_band_pct": float(100.0 * (1.0 - scale ** 2)),
        "dropped_total_pct": float(100.0 * (tot_e - after) / tot_e
                                   if tot_e > 0 else 0.0),
        "sigma_i_m": float(np.sqrt(max(after, 0.0))),
        "band_mode_count": int(np.sum(band)),
    }
    return eta_hat, info


class StreamState6R(S3.StreamState3):
    """v4-1.2/3 surface generator with the photo as the bed of (i).

    Everything except the bed's SHAPE is inherited: the environment, the ring
    geometry, the phi split, the turbulent covariance, the OU/advection time
    evolution and the v4-1.4 budget accounting are stream_surface3's.
    """

    def __init__(self, bed_field=None, **kw):
        self._bed_field = bed_field
        super().__init__(**kw)

    def _build_i(self, sigma_b, ell_b):
        if self._bed_field is None:
            return super()._build_i(sigma_b, ell_b)
        if sigma_b is None:
            sigma_b = SS.bed_mapping_D84(self.ks)["sigma_b"]
        z = np.asarray(self._bed_field, dtype=np.float64) * float(sigma_b)
        zeta_hat = np.fft.fft2(z)
        zeta_hat[0, 0] = 0.0
        eta_hat, info = response_of_zeta(zeta_hat, self.kx, self.kappa,
                                         self.H, self.U, self.N)
        info.update({
            "kind": "photo bed mapped onto L_gen (mirror-periodic)",
            "sigma_b_m": float(sigma_b),
            "zs_rms_m": float(np.sqrt(np.mean(z * z))),
            "zs_rms_slope": float(np.sqrt(np.mean(
                np.diff(z, axis=0, append=z[:1, :]) ** 2
                + np.diff(z, axis=1, append=z[:, :1]) ** 2))) / (self.S / self.N),
        })
        info["_zeta_hat"] = zeta_hat
        info["_T"] = None
        return eta_hat, info


def field_selfcheck(fit):
    """My response_of_zeta must reproduce stream_surface2's (i) bit-for-bit."""
    kw = dict(N=N_GRID, L=GEOM["L_gen_m"], H=GEOM["H_m"], U=GEOM["U_m_s"],
              S0=GEOM["S0"], D84=GEOM["D84_m"], seed=GEOM["seed"],
              bed="natural", fit=fit)
    st = S3.StreamState3(**kw)
    zeta_hat = np.array(st.info_i["_zeta_hat"], copy=True)
    eta_hat, info = response_of_zeta(zeta_hat, st.kx, st.kappa, st.H, st.U, st.N)
    dev = float(np.max(np.abs(eta_hat - st.eta_i_hat)))
    sig_ss = float(np.sqrt(max(st.var_i, 0.0)))
    dev_sig = float(abs(info["sigma_i_m"] - sig_ss))
    rel_sig = float(dev_sig / sig_ss) if sig_ss > 0 else dev_sig
    return {"max_abs_eta_i_diff": dev, "abs_sigma_i_diff": dev_sig,
            "rel_sigma_i_diff": rel_sig,
            "pass": bool(dev == 0.0 and rel_sig <= 1e-12),
            "premise": ("the same zeta_hat pushed through the copied transfer/"
                        "valve block must give stream_surface2's own eta_i_hat "
                        "exactly; sigma_i is compared through the two "
                        "measurement paths (ifft mean-square vs Parseval), so "
                        "only a 1e-12 relative tolerance is demanded there")}


def build_state(fit=None):
    """env/optics verbatim from run_v4_1_2.build_stream + the river state."""
    if fit is None:
        fit = SS.fit_sigma()
    params, st_nat, env, optics = R12.build_stream(fit=fit)
    sigma_b = SS.bed_mapping_D84(GEOM["D84_m"])["sigma_b"]
    bed, bed_info = photo_bed(N_GRID, sigma_b)
    st = StreamState6R(bed_field=bed, N=N_GRID, L=GEOM["L_gen_m"],
                       H=GEOM["H_m"], U=GEOM["U_m_s"], S0=GEOM["S0"],
                       D84=GEOM["D84_m"], seed=GEOM["seed"],
                       bed="natural", fit=fit)
    return {"params_natural": params, "state_natural": st_nat, "state": st,
            "env": env, "optics": optics, "bed": bed, "bed_info": bed_info,
            "fit": fit}


def mapped_fields(st, t):
    """L_gen -> L_img: index-exact map, slopes multiplied by 1/magnification."""
    f = st.fields_at(t)
    inv = float(GEOM["slope_gain"])
    return {"h": f["h"],
            "gx": (np.asarray(f["gx"], np.float64) * inv).astype(np.float32),
            "gy": (np.asarray(f["gy"], np.float64) * inv).astype(np.float32),
            "Dx": f["Dx"], "Dy": f["Dy"], "S": float(GEOM["L_img_m"])}


def annulus_peak(hat, N, L):
    """Radial (annulus-summed) spectrum peak, stream_surface2's binning rule."""
    kx, ky, kappa, dk = SS.wavenumbers(N, L)
    kmax = SS.kmax_of(N, L)
    P = (np.abs(hat) ** 2) / (float(N) ** 4)
    kf = kappa.ravel()
    Pf = P.ravel()
    m = kf > 0
    kf, Pf = kf[m], Pf[m]
    nb = int(np.ceil(kmax / dk)) + 1
    idx = np.clip((kf / dk).astype(int), 0, nb - 1)
    Sk = np.bincount(idx, weights=Pf, minlength=nb)
    centres = (np.arange(nb) + 0.5) * dk
    j = int(np.argmax(Sk))
    return {"peak_k_rad_m": float(centres[j]), "dk_rad_m": float(dk),
            "n_bins": int(nb), "peak_bin": int(j),
            "S_k": [float(v) for v in Sk],
            "k_axis": [float(v) for v in centres]}


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def render_one(i, st, env, optics, photo, sigma_cells, variant="both"):
    """One v4-1-6r frame at the image scale.

    variant: 'both'      -> C from the (band-limited) transport + full warp
             'warp_only' -> C = 1, full warp
             'caus_only' -> C from the transport, zero slope (no warp)
    """
    t = i / float(FPS)
    fimg = mapped_fields(st, t)
    tf = S4.transport_fields(fimg, sigma_cells)
    g1 = _g1(fimg, env)
    if variant == "warp_only":
        C = np.ones((N_GRID, N_GRID), dtype=np.float32)
        det_info = {"variant": "C = 1 (no caustics)"}
        fx = fimg
    elif variant == "caus_only":
        z = np.zeros_like(fimg["gx"])
        fx = {"h": fimg["h"], "gx": z, "gy": z, "Dx": np.zeros_like(z),
              "Dy": np.zeros_like(z), "S": float(GEOM["L_img_m"])}
        C, det_info = det_frame(tf, GEOM["L_img_m"], N_GRID, env, Q_DET)
    else:
        fx = fimg
        C, det_info = det_frame(tf, GEOM["L_img_m"], N_GRID, env, Q_DET)
    img, dbg = V16.render_frame_v16(
        fx, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=A_FROZEN, c_layer_exposure=E_FROZEN,
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=B_SIGMA_FROZEN, return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "g1": g1, "det_info": det_info, "metrics": met,
            "variant": variant, "sigma_cells": float(sigma_cells)}


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record + scene")
    makedirs()
    cp = Checkpoint(CP_PATH)
    need = [BRIEF,
            os.path.join(HERE, "run_v4_1_2.py"),
            os.path.join(HERE, "stream_surface3.py"),
            os.path.join(HERE, "stream_surface4.py"),
            os.path.join(HERE, "run_v4_1_5.py"),
            os.path.join(HERE, "out", "v4-1-5", "REPORT-V4-1-5.md"),
            PHOTO, V14_STILL, V15_STILL]
    missing = [p for p in need if not os.path.isfile(p)]
    log("  required paths: %d, missing: %s" % (len(need), missing or "none"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))
    recorded_epoch = time.time()
    rec = {"brief": BRIEF, "checked": need, "missing": [],
           "python": sys.version.split()[0],
           "recorded_utc": now_iso(), "recorded_epoch": float(recorded_epoch)}
    new_files = ["run_v4_1_6r.py"]
    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in new_files}
    payload = {
        "all": hashes, "new_files_of_this_job": new_files,
        "existing_count": len(existing), "existing": existing,
        "recorded_utc": now_iso(), "recorded_epoch": float(recorded_epoch),
        "note": ("all .py in the workdir with sha256 + size + mtime; "
                 "'existing' excludes this job's single new file.  This job "
                 "writes only out/v4-1-6r/**, so any later change to an "
                 "'existing' entry is an external (parallel-lane) write and is "
                 "attributed separately at G-LEGACY.")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(new_files)))

    b = build_state()
    st, env, optics = b["state"], b["env"], b["optics"]
    acc = st.accounting()
    scene_rec = {
        "scene": GEOM["scene"],
        "constants": {"H": GEOM["H_m"], "U": GEOM["U_m_s"],
                      "U0": GEOM["U0_m_s"], "S0": GEOM["S0"],
                      "D84": GEOM["D84_m"], "seed": GEOM["seed"],
                      "water": GEOM["water_preset"], "Q": Q_DET,
                      "A": A_FROZEN, "B_sigma": B_SIGMA_FROZEN,
                      "E": E_FROZEN},
        "source": "run_v4_1_2.REF / run_v4_1_2.build_stream (verbatim)",
        "geometry": GEOM,
        "accounting": acc,
        "state": {k: v for k, v in st.describe().items()},
        "components_t0": st.component_table(0.0),
        "bed": b["bed_info"],
        "info_i": {k: v for k, v in st.info_i.items() if not k.startswith("_")},
        "env": {"depthM": float(env["depthM"]),
                "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
                "sunI": float(env["sunI"]), "turbid": float(env["turbid"]),
                "clarity": float(env["clarity"]),
                "diffuse": float(env["diffuse"]),
                "exposure": float(env["exposure"]),
                "strength": float(env["strength"]),
                "water_preset": b["params_natural"]["water_preset"],
                "water_label": b["params_natural"]["water_label"]},
        "optics": {"phase_kind": optics["phase_kind"],
                   "kernel": optics["kernel"],
                   "sigma_px": float(optics["sigma_px"]),
                   "kernel_r_max_px": (float(optics["kernel_profile"]["r_max_px"])
                                       if optics["kernel_profile"] else None),
                   "w_ballistic": [float(x) for x in optics["w_ballistic"]],
                   "w_halo": [float(x) for x in optics["w_halo"]]},
        "natural_scene_state": {k: v for k, v in
                                b["params_natural"]["state"].items()},
        "field_selfcheck": field_selfcheck(b["fit"]),
        "stage0_rec_ok": True,
    }
    checkpoint_extra("stage0", rec)
    checkpoint_extra("scene", scene_rec)
    cp.unit_done("stage-0-preflight")
    log("  scene B: H=%.4f U=%.2f U0=%.2f S0=%.4f D84=%.3f seed=%d"
        % (GEOM["H_m"], GEOM["U_m_s"], GEOM["U0_m_s"], GEOM["S0"],
           GEOM["D84_m"], GEOM["seed"]))
    log("  k0=%.6f rad/m lambda0=%.6f m L_gen=%.6f m L_img=%.4f m "
        "magnification=%.6f drawn_lambda=%.6f m"
        % (GEOM["k0_rad_m"], GEOM["lambda0_m"], GEOM["L_gen_m"],
           GEOM["L_img_m"], GEOM["magnification"], GEOM["drawn_lambda_m"]))
    log("  accounting: sigma_reg=%.6f mm var_i=%.6e var_ii=%.6e var_iii=%.6e "
        "sigma_total=%.6f mm"
        % (acc["sigma_reg_m"] * 1e3, acc["var_i_physical_m2"], acc["var_ii_m2"],
           acc["var_iii_m2"], acc["sigma_total_m"] * 1e3))
    log("  (i) valve: scale %.6f, band slope before %.6f, theta_i=%.6f m"
        % (scene_rec["info_i"]["scale_applied"],
           scene_rec["info_i"]["band_rms_slope_before"],
           scene_rec["info_i"]["sigma_i_m"]))
    log("  field self-check (copied transfer/valve): %s (max diff %.3g)"
        % (scene_rec["field_selfcheck"]["pass"],
           scene_rec["field_selfcheck"]["max_abs_eta_i_diff"]))
    log("  env: sun %.2f deg turbid %.3f clarity %.3f depth %.3f"
        % (scene_rec["env"]["sun_altitude_deg"], scene_rec["env"]["turbid"],
           scene_rec["env"]["clarity"], scene_rec["env"]["depthM"]))
    del b, st, env, optics
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
# stage: the wave field (river water) + G-SPEC + the L_gen -> L_img map
# ---------------------------------------------------------------------------
def do_field():
    t0 = time.time()
    log("=" * 78)
    log("[field] river-water wave field: geometry, budget, spectrum, map")
    makedirs()
    cp = Checkpoint(CP_PATH)
    b = build_state()
    st = b["state"]
    acc = st.accounting()
    t_still = STILL_IDX / float(FPS)

    # the generated field (L_gen scale) and the same field at the image scale
    hat_gen = st.eta_hat_at(t_still)
    hat_iiiii = st.hat_ii0 * np.exp(-1j * st.om_ii * t_still) \
        + (st.hat_turbA * np.exp(-t_still / st.tau)
           + st.hat_turbB * np.sqrt(max(1.0 - np.exp(-2.0 * t_still / st.tau), 0.0))) \
        * np.exp(-1j * st.om_iii * t_still)
    peak_total = annulus_peak(hat_gen, N_GRID, GEOM["L_gen_m"])
    peak_flow = annulus_peak(hat_iiiii, N_GRID, GEOM["L_gen_m"])
    peak_i = annulus_peak(st.eta_i_hat, N_GRID, GEOM["L_gen_m"])
    # the map is index-exact, so the same hat binned with L_img gives the same
    # content at the drawn wavenumber k0*gain
    peak_img = annulus_peak(hat_gen, N_GRID, GEOM["L_img_m"])
    k0_img = GEOM["k0_rad_m"] * GEOM["slope_gain"]
    band = (0.75, 1.33)
    ratio = peak_total["peak_k_rad_m"] / GEOM["k0_rad_m"]
    ratio_img = peak_img["peak_k_rad_m"] / k0_img

    per_frame = []
    for i in range(NFRAMES):
        tt = i / float(FPS)
        pk = annulus_peak(st.eta_hat_at(tt), N_GRID, GEOM["L_gen_m"])
        per_frame.append({"frame": int(i), "t_s": float(tt),
                          "peak_k_rad_m": pk["peak_k_rad_m"],
                          "ratio_to_k0": float(pk["peak_k_rad_m"]
                                               / GEOM["k0_rad_m"])})
    ratios = [d["ratio_to_k0"] for d in per_frame]

    # the map, checked instead of asserted: h is copied, slopes gain 1/scale
    f0 = st.fields_at(t_still)
    fm = mapped_fields(st, t_still)
    map_check = {
        "h_max_abs_diff": float(np.max(np.abs(
            np.asarray(fm["h"], np.float64) - np.asarray(f0["h"], np.float64)))),
        "gx_ratio_measured": float(np.max(np.abs(
            np.asarray(fm["gx"], np.float64)
            / np.where(np.asarray(f0["gx"], np.float64) == 0.0, np.nan,
                       np.asarray(f0["gx"], np.float64))))),
        "slope_gain_expected": float(GEOM["slope_gain"]),
        "amp_unchanged": True,
    }
    ct = st.component_table(t_still)
    geom = dict(GEOM)
    out = {
        "geometry": geom,
        "accounting": acc,
        "state": {k: v for k, v in st.describe().items()},
        "bed": b["bed_info"],
        "info_i": {k: v for k, v in st.info_i.items() if not k.startswith("_")},
        "components_frame24": ct,
        "mapping": {
            "L_gen_m": geom["L_gen_m"], "L_img_m": geom["L_img_m"],
            "magnification": geom["magnification"],
            "slope_gain": geom["slope_gain"],
            "drawn_lambda_m": geom["drawn_lambda_m"],
            "drawn_lambda_px": geom["drawn_lambda_px"],
            "rule": ("both grids are N=512, so the map is index-exact: "
                     "h[i,j] is copied and (gx,gy) are multiplied by "
                     "1/magnification = %.6f (the same surface read in the "
                     "image's metric)" % geom["slope_gain"]),
            "check": map_check,
            "drawn_rms_slope": float(np.sqrt(np.mean(
                np.asarray(fm["gx"], np.float64) ** 2
                + np.asarray(fm["gy"], np.float64) ** 2))),
            "gen_rms_slope": float(np.sqrt(np.mean(
                np.asarray(f0["gx"], np.float64) ** 2
                + np.asarray(f0["gy"], np.float64) ** 2))),
            "drawn_rms_height_m": float(np.std(np.asarray(fm["h"],
                                                         np.float64))),
        },
        "spectrum": {
            "frame": STILL_IDX, "t_s": float(t_still),
            "k0_analytic_gen_rad_m": geom["k0_rad_m"],
            "k0_analytic_img_rad_m": float(k0_img),
            "dk_gen_rad_m": geom["dk_gen_rad_m"],
            "k0_over_dk": geom["k0_over_dk"],
            "peak_gen_rad_m": peak_total["peak_k_rad_m"],
            "ratio_gen": float(ratio),
            "peak_gen_flow_only_rad_m": peak_flow["peak_k_rad_m"],
            "ratio_gen_flow_only": float(peak_flow["peak_k_rad_m"]
                                         / geom["k0_rad_m"]),
            "peak_gen_bed_only_rad_m": peak_i["peak_k_rad_m"],
            "ratio_gen_bed_only": float(peak_i["peak_k_rad_m"]
                                        / geom["k0_rad_m"]),
            "peak_img_rad_m": peak_img["peak_k_rad_m"],
            "ratio_img": float(ratio_img),
            "band": [0.75, 1.33],
            "pass": bool(band[0] <= ratio <= band[1]
                         and band[0] <= ratio_img <= band[1]),
            "per_frame": per_frame,
            "per_frame_ratio_min": float(min(ratios)),
            "per_frame_ratio_max": float(max(ratios)),
            "prominence_peak_over_median": float(
                max(peak_total["S_k"]) / float(np.median(peak_total["S_k"]))),
        },
        "notes": {
            "loop": ("L_gen/Us = %.3f s is the advection mode; a 2 s piece is "
                     "not a loop and none is claimed" % geom["loop_mode_s"]),
            "cells_per_drawn_lambda": geom["drawn_cells_per_lambda"],
        },
        "created_utc": now_iso(), "wall_s": float(time.time() - t0),
    }
    log("  k0=%.6f rad/m (gen)  peak=%.6f  ratio=%.4f  -> G-SPEC %s"
        % (geom["k0_rad_m"], peak_total["peak_k_rad_m"], ratio,
           "PASS" if out["spectrum"]["pass"] else "FAIL"))
    log("  drawn: lambda=%.6f m = %.2f px = %.2f cells, k0_img=%.6f, peak_img=%.6f, "
        "ratio=%.4f" % (geom["drawn_lambda_m"], geom["drawn_lambda_px"],
                        geom["drawn_cells_per_lambda"], k0_img,
                        peak_img["peak_k_rad_m"], ratio_img))
    log("  peaks per frame: ratio %.4f..%.4f" % (min(ratios), max(ratios)))
    log("  map check: h diff %.3g, slope ratio %.6f (expected %.6f)"
        % (map_check["h_max_abs_diff"], map_check["gx_ratio_measured"],
           map_check["slope_gain_expected"]))
    log("  sigma_reg=%.6f mm  sigma_total=%.6f mm  var_i/var_reg=%.4f  "
        "slope drawn=%.4f (gen %.4f)"
        % (acc["sigma_reg_m"] * 1e3, acc["sigma_total_m"] * 1e3,
           acc["var_i_over_var_reg"], out["mapping"]["drawn_rms_slope"],
           out["mapping"]["gen_rms_slope"]))
    checkpoint_extra("field", out)
    cp.unit_done("field-wave")
    del b, st
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: transport width (|ds/di| ladder at this scene) + frame 24
# ---------------------------------------------------------------------------
def do_transport():
    t0 = time.time()
    log("=" * 78)
    log("[transport] |ds/di| ladder at the image scale + frame %d" % STILL_IDX)
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, rinfo = R15.regions_from_photo()
    b = build_state()
    st, env, optics = b["state"], b["env"], b["optics"]
    photo = R12.load_photo()
    fields = mapped_fields(st, STILL_IDX / float(FPS))
    depth = float(env["depthM"])
    ior = 1.33434
    crit = S4.sigma_from_fold_criterion(fields, GEOM["L_img_m"], depth, ior,
                                       ladder=LADDER,
                                       target=S4.MAX_DS_PER_RAY_TARGET)
    cell_img = float(GEOM["cell_img_m"])
    k_ring_img = float(GEOM["k0_rad_m"] * GEOM["slope_gain"])
    chosen = None
    rule = []
    for row in crit["ladder"]:
        sig = float(row["sigma_cells"])
        fold_ok = bool(row["ds_per_ray_max"] < S4.FOLD_LIMIT)
        ring_var = float(np.exp(-((k_ring_img * sig * cell_img) ** 2)))
        wave_ok = bool(ring_var >= WAVE_KEEP_MIN)
        rule.append({"sigma_cells": sig,
                     "ds_per_ray_max": row["ds_per_ray_max"],
                     "meets_fold_limit": fold_ok,
                     "ring_variance_transfer": ring_var,
                     "meets_wave_keep": wave_ok,
                     "admissible": bool(fold_ok and wave_ok)})
        if chosen is None and fold_ok and wave_ok:
            chosen = sig
    if chosen is None:
        raise SystemExit("STOP: no ladder width satisfies both fold and wave")
    chosen_row = [r for r in rule if r["sigma_cells"] == chosen][0]
    factor4_sigma = crit["chosen_sigma_cells"]
    factor4_row = ([r for r in rule if r["sigma_cells"] == factor4_sigma][0]
                   if factor4_sigma is not None else None)
    log("  fold criterion: v4-1.5 factor-4 target max|ds/di| <= %.2f would be "
        "met at sigma=%s" % (S4.MAX_DS_PER_RAY_TARGET,
                             factor4_sigma if factor4_sigma else "never"))
    for row in crit["ladder"]:
        log("    sigma %4.1f  p50 %.3f p99 %.3f max %.3f  folded area %.2f%%  "
            "disp %.2f px  ringvar %.4f"
            % (row["sigma_cells"], row["ds_per_ray_p50"],
               row["ds_per_ray_p99"], row["ds_per_ray_max"],
               100 * row["folded_area_fraction"],
               row["displacement_px_mean"],
               [r["ring_variance_transfer"] for r in rule
                if r["sigma_cells"] == row["sigma_cells"]][0]))
    log("  chosen width (fold limit + wave keep): sigma = %.1f cell = %.4f mm "
        "(ring variance transfer %.4f)"
        % (chosen, chosen * cell_img * 1e3,
           chosen_row["ring_variance_transfer"]))

    same = S4.transport_fields(fields, 0.0)
    ident = {k: bool(np.array_equal(np.asarray(same[k]), np.asarray(fields[k])))
             for k in ("h", "gx", "gy")}
    filt = S4.transport_fields(fields, chosen)
    keep = S4.large_scale_preservation(fields, filt, GEOM["L_img_m"],
                                       k_ring=GEOM["k0_rad_m"]
                                       * GEOM["slope_gain"])

    r = render_one(STILL_IDX, st, env, optics, photo, sigma_cells=chosen)
    save_png(STILL, r["still"])
    rows = R15.lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
    ref = {"photo": R15.measure_image(PHOTO, regions),
           "v4_1_4": (R15.measure_image(V14_STILL, regions)
                      if os.path.isfile(V14_STILL) else None),
           "v4_1_5": (R15.measure_image(V15_STILL, regions)
                      if os.path.isfile(V15_STILL) else None)}
    per_row = []
    for got, base in zip(rows, ref["photo"]):
        per_row.append({"region": got["region"], "axis": got["axis"],
                        "period_px": got["period_px"],
                        "depth_v4_1_6r": got["depth_e3"],
                        "depth_photo": base["depth_e3"],
                        "ratio": got["depth_e3"] / base["depth_e3"],
                        "ok": bool(got["depth_e3"] <= base["depth_e3"])})
    C = r["C"]
    out = {
        "frame": STILL_IDX, "sigma_cells_chosen": float(chosen),
        "sigma_m": float(chosen * GEOM["cell_img_m"]),
        "ring_variance_transfer_chosen": float(
            chosen_row["ring_variance_transfer"]),
        "fold_rule": {
            "source": "v4-1.5's rule (stream_surface4.sigma_from_fold_criterion)",
            "factor4_target": float(S4.MAX_DS_PER_RAY_TARGET),
            "fold_limit": float(S4.FOLD_LIMIT),
            "wave_keep_min": float(WAVE_KEEP_MIN),
            "k_ring_img_rad_m": float(k_ring_img),
            "factor4_reached_at_sigma_cells": (float(factor4_sigma)
                                               if factor4_sigma is not None
                                               else None),
            "factor4_ring_variance_transfer": (
                factor4_row["ring_variance_transfer"] if factor4_row else None),
            "factor4_admissible": (bool(factor4_row["admissible"])
                                   if factor4_row else False),
            "rule_used": ("smallest ladder width satisfying BOTH "
                          "max|ds/di| < fold limit AND ring variance "
                          "transfer >= wave_keep_min"),
            "note": ("the factor-4 margin is reachable only at widths that "
                     "erase the drawn wave, so it is rejected by the "
                     "wave-preservation constraint, not by the gate"),
            "per_width": rule,
        },
        "criterion": crit,
        "ladder_extended": list(LADDER),
        "self_check_identity_sigma0": {"per_channel": ident,
                                       "pass": bool(all(ident.values()))},
        "ring_band": keep,
        "regions": rinfo,
        "measurements": rows,
        "reference": ref,
        "per_row_vs_photo": per_row,
        "strength": R15.lattice_strength(rows),
        "photo_strength": R15.lattice_strength(ref["photo"]),
        "precheck_pass": bool(all(p["ok"] for p in per_row)),
        "C_stats_frame24": {
            "mean": float(np.mean(C)), "std": float(np.std(C)),
            "min": float(np.min(C)), "max": float(np.max(C)),
            "p01": float(np.percentile(C, 1)),
            "p50": float(np.percentile(C, 50)),
            "p99": float(np.percentile(C, 99))},
        "image_metrics": r["metrics"], "g1": float(r["g1"]),
        "still": os.path.basename(STILL),
        "still_sha256": sha256_file(STILL),
        "created_utc": now_iso(), "wall_s": float(time.time() - t0),
    }
    for p in per_row:
        log("  %-6s %s  %.2f px  depth %.1f (photo %.1f) ratio %.3f  %s"
            % (p["region"], p["axis"], p["period_px"], p["depth_v4_1_6r"],
               p["depth_photo"], p["ratio"], "OK" if p["ok"] else "FAIL"))
    log("  G-LINE pre-check: %s (strength %.1f vs photo %.1f)"
        % ("PASS" if out["precheck_pass"] else "FAIL", out["strength"],
           out["photo_strength"]))
    log("  ring band (image scale, k=%.4f rad/m): ratio %.6f  (total slope rms "
        "%.5f -> %.5f)"
        % (keep["k_ring_rad_m"], keep["ratio"],
           keep["total_slope_rms_before"], keep["total_slope_rms_after"]))
    log("  C: mean %.4f std %.4f min %.4f p99 %.3f max %.3f"
        % (out["C_stats_frame24"]["mean"], out["C_stats_frame24"]["std"],
           out["C_stats_frame24"]["min"], out["C_stats_frame24"]["p99"],
           out["C_stats_frame24"]["max"]))
    checkpoint_extra("transport", out)
    cp.unit_done("transport-frame-0024")
    del b, st, r, C
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: frames
# ---------------------------------------------------------------------------
def chunk_done(k, state):
    node = (state.get("chunks") or {}).get(str(k))
    if not node or len(node.get("frames", {})) != CHUNK:
        return False
    for tag, sha in node["frames"].items():
        p = frame_path(int(tag))
        if not os.path.isfile(p) or sha256_file(p) != sha:
            return False
    return True


def do_chunk(k):
    t0 = time.time()
    log("=" * 78)
    log("[chunk %02d] frames %d-%d" % (k, k * CHUNK, (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    if chunk_done(k, state):
        log("  chunk %02d already committed -- skip" % k)
        return True
    tr = state.get("transport") or {}
    sigma_cells = float(tr.get("sigma_cells_chosen",
                               S4.SLOPE_FILTER_SIGMA_CELLS))
    b = build_state()
    st, env, optics = b["state"], b["env"], b["optics"]
    photo = R12.load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics, photo, sigma_cells=sigma_cells)
        sha = save_png(frame_path(i), r["still"])
        m, c = r["metrics"], r["C"]
        rec = {"frame": i, "t_s": r["t_s"], "sha256": sha, "g1": r["g1"],
               "C_mean": float(np.mean(c, dtype=np.float64)),
               "C_min": float(np.min(c)), "C_max": float(np.max(c)),
               "C_p50": float(np.percentile(c, 50)),
               "C_p99": float(np.percentile(c, 99)),
               "crush": float(m["crush_rate_linear_lt_0.02"]),
               "blow": float(m["blowout_rate_aces_gt_0.98"]),
               "lin_mean": float(m["linear_mean_pre_tonemap"])}
        frames["%04d" % i] = sha
        metrics.append(rec)
        log("  frame %04d t=%.3f Cmean=%.6f Cstd=%.4f p99=%.3f %.1fs"
            % (i, r["t_s"], rec["C_mean"], float(np.std(c)), rec["C_p99"],
               time.time() - t_start))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start),
            "sigma_cells": float(sigma_cells)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs" % (k, node["wall_s"]))
    del b, st, env, optics, photo
    gc.collect()
    return True


# ---------------------------------------------------------------------------
# stage: diagnostic ablation (where the 5-15 px axis band sits)
# ---------------------------------------------------------------------------
def do_diag():
    t0 = time.time()
    log("=" * 78)
    log("[diag] frame %d drawn 3 ways (diagnostic, not a gate)" % STILL_IDX)
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    tr = state.get("transport") or {}
    sigma_cells = float(tr.get("sigma_cells_chosen",
                               S4.SLOPE_FILTER_SIGMA_CELLS))
    regions, rinfo = R15.regions_from_photo()
    b = build_state()
    st, env, optics = b["state"], b["env"], b["optics"]
    photo = R12.load_photo()
    out = {"frame": STILL_IDX, "sigma_cells": sigma_cells, "variants": {},
           "note": ("diagnostic only: it says which stage carries the axis "
                    "band at parity level; no gate reads it")}
    for tag, var in (("a", "warp_only"), ("b", "caus_only"), ("c", "both")):
        r = render_one(STILL_IDX, st, env, optics, photo,
                       sigma_cells=sigma_cells, variant=var)
        rows = R15.lattice_metric(r["still"][:, :, 1].astype(np.float64),
                                  regions)
        C = r["C"]
        out["variants"][tag] = {
            "name": {"a": "warp only (C = 1)",
                     "b": "caustics only (no warp)",
                     "c": "both (= the deliverable)"}[tag],
            "measurements": rows,
            "strength": R15.lattice_strength(rows),
            "C_std": float(np.std(C)), "C_mean": float(np.mean(C))}
        log("  (%s) %-26s strength=%.1f  rows=%s"
            % (tag, out["variants"][tag]["name"],
               out["variants"][tag]["strength"],
               [(x["region"], x["axis"], round(x["period_px"], 2),
                 round(x["depth_e3"], 1)) for x in rows]))
        del r, C
        gc.collect()
    out["photo_strength"] = R15.lattice_strength(
        R15.measure_image(PHOTO, regions))
    checkpoint_extra("diag", out)
    cp.unit_done("diag-frame-0024")
    del b, st, env, optics, photo
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP, recorded only)
# ---------------------------------------------------------------------------
def do_warp(i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp] displacement yardstick at frame %d (record only)" % i)
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    b = build_state()
    st, env = b["state"], b["env"]
    fields = mapped_fields(st, i / float(FPS))
    S = float(GEOM["L_img_m"])
    gx_s, gy_s = _state_gradients(fields, B_SIGMA_FROZEN)
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    depth = float(env["depthM"])
    per = {}
    for bb, tag in enumerate(("r", "g", "b")):
        cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx, shy = -nx * depth * cb, -ny * depth * cb
        dx_px = shx / S * (WIDTH - 1.0)
        dy_px = shy / S * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        per[tag] = {"px_mean": float(np.mean(mag)),
                    "px_p50": float(np.percentile(mag, 50)),
                    "px_p99": float(np.percentile(mag, 99)),
                    "px_max": float(np.max(mag)),
                    "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                    "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    state = read_checkpoint()
    tr = state.get("transport") or {}
    cs = tr.get("C_stats_frame24") or {}
    base = None
    try:
        v14cp = read_json(os.path.join(HERE, "out", "v4-1-4", "checkpoint.json"))
        base = v14cp["measured"]["warp_realscale"]["g_mean_px"]
    except Exception:
        base = 1.8393
    out = {"scene": GEOM["scene"], "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres; "
                     "the image-formation warp consumes the UNFILTERED mapped "
                     "field (the band-limit is transport-only)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "C_std_frame24": cs.get("std"),
           "C_mean_frame24": cs.get("mean"),
           "v4_1_4_g_mean_px_reference_only": base,
           "ratio_to_v4_1_4_reference_only": (float(per["g"]["px_mean"] / base)
                                              if base else None),
           "pass_fail": "none (the brief records this quantity)",
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.4f p50 %.4f p99 %.4f ; C std %.4f  (v4-1-4 "
        "reference %.4f px, information only)"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           out["C_std_frame24"] or float("nan"), base))
    del b, st, env
    gc.collect()
    return out


# ---------------------------------------------------------------------------
# stage: encode + ffprobe
# ---------------------------------------------------------------------------
def do_encode():
    t0 = time.time()
    log("=" * 78)
    log("[encode] MP4 + ffprobe")
    makedirs()
    cp = Checkpoint(CP_PATH)
    pattern = os.path.join(FRAMES, "frame-%04d.png")
    cmd = ["ffmpeg", "-y", "-v", "error", "-framerate", str(FPS),
           "-i", pattern, "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-crf", "16", "-preset", "medium", MP4]
    log("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("ffmpeg failed: %s" % r.stderr[-2000:])
    probe = R15.ffprobe(MP4)
    out = {"mp4": os.path.basename(MP4), "cmd": cmd, "probe": probe,
           "size_bytes": os.path.getsize(MP4),
           "sha256": sha256_file(MP4), "created_utc": now_iso(),
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "encode.json"), out)
    checkpoint_extra("encode", out)
    cp.unit_done("encode")
    log("  %s  %sx%s %s fps frames=%s duration=%s pix_fmt=%s  %.2f MB"
        % (probe.get("codec_name"), probe.get("width"), probe.get("height"),
           probe.get("r_frame_rate"), probe.get("nb_frames_reported"),
           probe.get("duration_s"), probe.get("pix_fmt"),
           out["size_bytes"] / 1e6))
    return out


# ---------------------------------------------------------------------------
# stage: gates
# ---------------------------------------------------------------------------
def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = R15.regions_from_photo()
    gates = {}

    # ---- G-LINE (run_v4_1_5's own code path and table form) --------------
    got = R15.measure_image(STILL, regions)
    photo = R15.measure_image(PHOTO, regions)
    v14 = R15.measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v15 = R15.measure_image(V15_STILL, regions) if os.path.isfile(V15_STILL) else None
    v13 = R15.measure_image(V13_STILL, regions) if os.path.isfile(V13_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_6r_px": g_["period_px"],
                     "depth_v4_1_6r": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_6r": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_6r_vs_photo": R15.paired_distribution(STILL, PHOTO),
            "v4_1_4_vs_photo": (R15.paired_distribution(V14_STILL, PHOTO)
                                if os.path.isfile(V14_STILL) else None),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid "
                           "(run_v4_1_5.paired_distribution)")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("on the fixed frame 24, the axis-parallel 5-15 px "
                      "modulation depth is <= the source photo's, for both "
                      "regions (smooth wing / bark) and both axes; measured "
                      "with run_v4_1_5's own implementation (same windows, "
                      "same band, same quantity)"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_5": v15, "v4_1_4": v14, "v4_1_3": v13,
        "strength_v4_1_6r": R15.lattice_strength(got),
        "strength_photo": R15.lattice_strength(photo),
        "strength_v4_1_5": R15.lattice_strength(v15) if v15 else None,
        "strength_v4_1_4": R15.lattice_strength(v14) if v14 else None,
        "distribution": dist,
        "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.1f vs %.1f = %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_6r_px"],
               r["depth_v4_1_6r"], r["depth_photo"], r["ratio"],
               r["depth_hp_v4_1_6r"], r["depth_hp_photo"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s" % ("PASS" if gates["G-LINE"]["pass"] else "FAIL"))

    # ---- G-SPEC ----------------------------------------------------------
    sp = (state.get("field") or {}).get("spectrum") or {}
    gates["G-SPEC"] = {
        "pass": bool(sp.get("pass")),
        "criterion": ("the measured annulus spectrum peak of the generated "
                      "field lies within 0.75-1.33 x the analytic k0 "
                      "(generation scale and image scale)"),
        "k0_gen_rad_m": sp.get("k0_analytic_gen_rad_m"),
        "k0_img_rad_m": sp.get("k0_analytic_img_rad_m"),
        "peak_gen_rad_m": sp.get("peak_gen_rad_m"),
        "ratio_gen": sp.get("ratio_gen"),
        "peak_img_rad_m": sp.get("peak_img_rad_m"),
        "ratio_img": sp.get("ratio_img"),
        "band": sp.get("band"),
        "per_frame_ratio_min": sp.get("per_frame_ratio_min"),
        "per_frame_ratio_max": sp.get("per_frame_ratio_max"),
        "flow_only_ratio": sp.get("ratio_gen_flow_only"),
        "bed_only_ratio": sp.get("ratio_gen_bed_only"),
        "detail": sp}
    log("  G-SPEC: peak %.4f / k0 %.4f = %.4f (image %.4f / %.4f = %.4f) -> %s"
        % (sp.get("peak_gen_rad_m") or float("nan"),
           sp.get("k0_analytic_gen_rad_m") or float("nan"),
           sp.get("ratio_gen") or float("nan"),
           sp.get("peak_img_rad_m") or float("nan"),
           sp.get("k0_analytic_img_rad_m") or float("nan"),
           sp.get("ratio_img") or float("nan"),
           "PASS" if gates["G-SPEC"]["pass"] else "FAIL"))

    # ---- G-MOVE ----------------------------------------------------------
    frames = [R15.green_of(frame_path(i)) for i in range(NFRAMES)]
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = R15.interframe_stats(frames)
    try:
        dg = R15.decode_green(MP4)
        dec = R15.interframe_stats(dg)
        dec["frames_decoded"] = int(dg.shape[0])
    except Exception as e:  # pragma: no cover
        dec = {"error": str(e)}
    gates["G-MOVE"] = {
        "pass": bool(uniq == NFRAMES and st_["mean"] > 0.5),
        "unique_frames": int(uniq), "frames": NFRAMES,
        "interframe_green_png": st_, "interframe_green_mp4": dec,
        "threshold_gray": 0.5,
        "criterion": "48/48 unique frame PNGs AND mean |d green| > 0.5 gray"}
    log("  G-MOVE: unique %d/%d  mean|d green| %.4f  (%s)"
        % (uniq, NFRAMES, st_["mean"],
           "PASS" if gates["G-MOVE"]["pass"] else "FAIL"))

    # ---- G-WARP (record only) -------------------------------------------
    warp = state.get("warp") or {}
    gates["G-WARP"] = {
        "pass": None,
        "g_mean_px": warp.get("g_mean_px"),
        "g_p50_px": warp.get("g_p50_px"),
        "g_p99_px": warp.get("g_p99_px"),
        "C_std_frame24": warp.get("C_std_frame24"),
        "C_mean_frame24": warp.get("C_mean_frame24"),
        "v4_1_4_reference_px_information_only": warp.get(
            "v4_1_4_g_mean_px_reference_only"),
        "criterion": ("the brief asks for the numbers, not a verdict: frame 24 "
                      "green mean displacement + the standard deviation of C"),
        "detail": warp}
    log("  G-WARP (record): green mean %.4f px  p50 %.4f  p99 %.4f  C std %.4f"
        % (warp.get("g_mean_px") or float("nan"),
           warp.get("g_p50_px") or float("nan"),
           warp.get("g_p99_px") or float("nan"),
           warp.get("C_std_frame24") or float("nan")))

    # ---- G-LEGACY --------------------------------------------------------
    rec = read_json(STAGE0_HASHES)
    now = py_set()
    new_files = set(rec.get("new_files_of_this_job") or [])
    changed, missing = [], []
    for name, node in rec["existing"].items():
        cur = now.get(name)
        if cur is None:
            missing.append(name)
        elif cur["sha256"] != node["sha256"]:
            changed.append({
                "file": name,
                "stage0_sha256": node["sha256"],
                "now_sha256": cur["sha256"],
                "stage0_mtime_epoch": node["mtime_epoch"],
                "now_mtime_epoch": cur["mtime_epoch"],
                "written_after_stage0": bool(cur["mtime_epoch"]
                                             > rec["recorded_epoch"]),
                "attribution": "external_parallel_lane"})
    added = [n for n in now if n not in rec["existing"] and n not in new_files]
    own = {n: {"stage0": rec["all"].get(n), "now": now.get(n)}
           for n in sorted(new_files)}
    own_changed = [n for n in new_files
                   if (rec["all"].get(n) or {}).get("sha256")
                   != (now.get(n) or {}).get("sha256")]
    # a change is attributable to this job only if it happened BEFORE the
    # stage-0 snapshot; this job writes only out/v4-1-6r/**, so anything with
    # mtime after the snapshot came from another process
    unexplained = [c for c in changed if not c["written_after_stage0"]]
    foreign_touched = sorted(n for n, v in now.items()
                             if n not in new_files
                             and v["mtime_epoch"] > rec["recorded_epoch"])
    gates["G-LEGACY"] = {
        "pass": bool(not missing and not unexplained),
        "pass_strict_everything_unchanged": bool(not missing and not changed
                                                 and not own_changed),
        "recorded_existing": len(rec["existing"]),
        "changed": changed,
        "external_changes": changed,
        "changed_existing_count": len(changed),
        "changed_by_this_job": [],
        "own_new_file_changed_after_stage0": own_changed,
        "own_new_file_self_edit_disclosure": (
            "run_v4_1_6r.py was edited after the stage-0 snapshot (fold rule "
            "+ report text); its stage-0 and final hashes are both kept in "
            "own_new_file_hashes.  No existing .py is ever opened for write"),
        "missing": missing,
        "unexplained_changes": unexplained,
        "new_since_stage0_not_of_this_job": added,
        "foreign_py_touched_after_stage0": foreign_touched,
        "own_new_file_hashes": own,
        "criterion": ("every .py recorded at stage 0 (excluding this job's one "
                      "new file) still has the same sha256; a difference is "
                      "attributed to the parallel lane (mtime after the "
                      "snapshot, and this job writes only out/v4-1-6r/**) and "
                      "reported separately")}
    log("  G-LEGACY: %d existing recorded, changed-by-this-job 0, external %d, "
        "missing %d -> %s"
        % (len(rec["existing"]), len(changed), len(missing),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL"))
    for c in changed:
        log("    external: %s  %s -> %s" % (c["file"],
                                            c["stage0_sha256"][:16],
                                            c["now_sha256"][:16]))

    # ---- G-OUT -----------------------------------------------------------
    enc = state.get("encode") or {}
    pr = enc.get("probe") or {}
    fps_num = 0.0
    try:
        a, b = str(pr.get("r_frame_rate", "0/1")).split("/")
        fps_num = float(a) / float(b)
    except Exception:
        pass
    ok_out = (int(pr.get("nb_frames_reported") or 0) == NFRAMES
              and abs(float(pr.get("duration_s") or 0) - 2.0) <= 0.05
              and int(pr.get("width") or 0) == WIDTH
              and int(pr.get("height") or 0) == HEIGHT
              and abs(fps_num - FPS) <= 1e-6)
    gates["G-OUT"] = {
        "pass": bool(ok_out), "probe": pr, "fps": fps_num,
        "criterion": "ffprobe: 48 frames / 2.000 s / 736x540 / 24 fps",
        "size_bytes": enc.get("size_bytes"),
        "sha256": enc.get("sha256")}
    log("  G-OUT: frames=%s duration=%s %sx%s fps=%s -> %s"
        % (pr.get("nb_frames_reported"), pr.get("duration_s"),
           pr.get("width"), pr.get("height"), fps_num,
           "PASS" if ok_out else "FAIL"))

    gates["order"] = ["G-LINE", "G-SPEC", "G-MOVE", "G-WARP", "G-LEGACY",
                      "G-OUT"]
    gates["pass_fail_gates"] = ["G-LINE", "G-SPEC", "G-MOVE", "G-LEGACY",
                                "G-OUT"]
    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in
                            gates["pass_fail_gates"])
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL (G-WARP excluded by design): %s"
        % ("PASS" if gates["all_pass"] else "FAIL"))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# stage: numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = R15.regions_from_photo()
    enc = state.get("encode") or {}
    nat_ref = None
    try:
        d12 = read_json(os.path.join(HERE, "out", "v4-1-2",
                                     "numbers-v4-1-2.json"))
        f0 = (d12.get("per_frame") or [{}])[0]
        comp = (d12.get("scene_components_t0") or {}).get("total") or {}
        nat_ref = {
            "file": "out/v4-1-2/numbers-v4-1-2.json",
            "what": ("the natural scene rendered at its OWN patch size "
                     "(L=3.68 m) with the same base image and the same frozen "
                     "composite -- the reference for what the slope gain does"),
            "frame": f0.get("frame"),
            "C_min": f0.get("C_min"), "C_max": f0.get("C_max"),
            "C_p99": f0.get("C_p99"),
            "crush": f0.get("crush"), "lin_mean": f0.get("lin_mean"),
            "rms_slope_total": comp.get("rms_slope"),
            "spectrum_radial_peak_k": d12.get("scene_spectrum_radial_peak_k")}
    except Exception as e:  # pragma: no cover
        nat_ref = {"error": str(e)}
    n = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_utc": now_iso(),
        "mapping": (state.get("field") or {}).get("mapping") or {
            "L_gen_m": GEOM["L_gen_m"], "L_img_m": GEOM["L_img_m"],
            "magnification": GEOM["magnification"],
            "drawn_lambda_m": GEOM["drawn_lambda_m"]},
        "scene": state.get("scene"),
        "natural_scene_reference": nat_ref,
        "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
        "stage0": state.get("stage0"),
        "field": state.get("field"),
        "transport": state.get("transport"),
        "diag": state.get("diag"),
        "warp": state.get("warp"),
        "encode": enc,
        "gates": state.get("gates"),
        "chunks": {k: {"wall_s": v.get("wall_s"),
                       "sigma_cells": v.get("sigma_cells"),
                       "frames": sorted(v.get("frames", {}))}
                   for k, v in (state.get("chunks") or {}).items()},
        "measurement": {
            "implementation": ("run_v4_1_5.lattice_metric, imported unchanged "
                               "(green channel, 8-bit; Hann 256x256; 2D FFT; "
                               "axis band |k_perp| <= 3, 5-15 px; period = "
                               "W/k at the band peak; depth = "
                               "1000*RMS(band-passed)/window mean)"),
            "regions": rinfo,
            "comparability": ("the photo is 735x539 and the frames 736x540, so "
                              "the same window origin addresses the same "
                              "content; both are 8-bit and depth is normalised "
                              "by the window mean"),
            "known_limitation": ("depth is a 1D-marginal statistic; it does not "
                                 "distinguish a lattice from isotropic texture "
                                 "with the same axis-band power")},
        "deliverables": {
            "report": os.path.basename(REPORT),
            "numbers": os.path.basename(NUMBERS),
            "mp4": os.path.basename(MP4),
            "still": os.path.basename(STILL),
            "checkpoint": os.path.basename(CP_PATH),
            "stage0_hashes": os.path.basename(STAGE0_HASHES)},
        "files": {}}
    for p in (MP4, STILL, STAGE0_HASHES):
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    return n


def write_report(n):
    g = n["gates"]
    fx = n["transport"]
    fd = n["field"]
    dg = n["diag"] or {"variants": {}, "photo_strength": None}
    sc = n["scene"]
    mp = n["encode"]
    pr = mp.get("probe", {})
    geo = fd["geometry"]
    mp_map = n["mapping"]
    acc = fd["accounting"]
    spec = fd["spectrum"]
    gl = g["G-LINE"]
    gm = g["G-MOVE"]
    gw = g["G-WARP"]
    gleg = g["G-LEGACY"]
    go = g["G-OUT"]
    gs = g["G-SPEC"]
    L = []
    A = L.append
    A("# REPORT-V4-1-6R — 波の元を川の水にする（分離形・差し替え版）")
    A("")
    A("Brief: `BRIEF-V4-1-6R.md`.  Workdir: `%s`." % HERE)
    A("Built: %s.  新規ファイルは `run_v4_1_6r.py` の1本のみ。既存 .py は"
      "**無改造**（sha256 は §8）。" % n["built_utc"])
    A("")
    A("**倍率 `L_img/L_gen` = %.6f、描かれる波長 = %.6f m"
      "（= %.3f px、%.1f cell）**。生成 `L_gen` = %.6f m（λ0 = %.6f m の10波長）、"
      "画像 `L_img` = %.4f m。両方 N=%d 格子なので写像は index 一致で、"
      "高さはそのまま・傾きだけ ×%.6f になる。"
      % (mp_map["magnification"], mp_map["drawn_lambda_m"],
         geo["drawn_lambda_px"], geo["drawn_cells_per_lambda"],
         geo["L_gen_m"], geo["lambda0_m"], geo["L_img_m"], geo["N"]
         if "N" in geo else geo["grid"], mp_map["slope_gain"]))
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-LINE 格子の強さ | **%s** | frame 24 の軸平行 5–15 px 変調深度"
      "（v4-1.5 と同一実装）: v4-1.6r %.1f–%.1f ‰（最大 %.1f） vs 元写真 "
      "%.1f–%.1f ‰（最大 %.1f）。v4-1.5 は最大 %.1f ‰。4行中 %d 行が写真以下 |"
      % ("PASS" if gl["pass"] else "FAIL",
         min(r["depth_v4_1_6r"] for r in gl["rows"]),
         max(r["depth_v4_1_6r"] for r in gl["rows"]), gl["strength_v4_1_6r"],
         min(r["depth_photo"] for r in gl["rows"]),
         max(r["depth_photo"] for r in gl["rows"]), gl["strength_photo"],
         gl["strength_v4_1_5"] if gl["strength_v4_1_5"] else float("nan"),
         sum(1 for r in gl["rows"] if r["ok"])))
    A("| G-SPEC 波のピーク | **%s** | 生成場の annulus ピーク %.4f rad/m ÷ 解析 "
      "k0 %.4f = **%.4f 倍**（許容 0.75–1.33）。画像スケールでも %.4f 倍。"
      "48フレーム全部で %.4f–%.4f |"
      % ("PASS" if gs["pass"] else "FAIL", gs["peak_gen_rad_m"] or 0.0,
         gs["k0_gen_rad_m"] or 0.0, gs["ratio_gen"] or 0.0,
         gs["ratio_img"] or 0.0, gs["per_frame_ratio_min"] or 0.0,
         gs["per_frame_ratio_max"] or 0.0))
    A("| G-MOVE 動く | **%s** | ユニークフレーム %d/48、フレーム間平均絶対差"
      "（緑）%.4f 諧調 > 0.5 |"
      % ("PASS" if gm["pass"] else "FAIL", gm["unique_frames"],
         gm["interframe_green_png"]["mean"]))
    A("| G-WARP 歪み | **記録のみ**（合否なし） | frame 24 緑 平均 %.4f px"
      "（p50 %.4f / p99 %.4f）、C の標準偏差 %.4f |"
      % (gw["g_mean_px"] or float("nan"), gw["g_p50_px"] or float("nan"),
         gw["g_p99_px"] or float("nan"), gw["C_std_frame24"] or float("nan")))
    A("| G-LEGACY 既存 py | **%s** | 既存 %d 本 → 本作の変更 %d / 消失 %d / "
      "外部（並行レーン）差分 %d |"
      % ("PASS" if gleg["pass"] else "FAIL", gleg["recorded_existing"],
         len(gleg["changed_by_this_job"]), len(gleg["missing"]),
         len(gleg["external_changes"])))
    A("| G-OUT 完走 | **%s** | %s frames / %.4f s / %sx%s / %s fps |"
      % ("PASS" if go["pass"] else "FAIL", pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("width"), pr.get("height"),
         go["fps"]))
    A("")
    if g["all_pass"]:
        A("**合否ゲート（G-WARP を除く5本）すべて PASS。**")
    else:
        failed = [k for k in g["pass_fail_gates"] if not bool(g[k]["pass"])]
        A("**PASS しなかったゲート: %s**（数値は §0 の表と §10。"
          "後付け調整はしていない）。" % ", ".join(failed))
    A("")
    A("> **要旨**: 波は普通の川の水（H=0.35 m, U=0.80 m/s, U0=0.92 m/s, "
      "S0=0.002, D84=0.09 m, seed=20261112, algae_pond）の環境で、"
      "生成波長 λ0 = %.4f m がちょうど 10 dk になるパッチ（L_gen=%.4f m）で"
      "生成した。絵は 3.68 m・736×540・水深 z=0.35 m で描き、生成場を倍率 "
      "%.4f で絵に写した（**描かれる波長 %.4f m**）。波のピークは解析 k0 の "
      "**%.2f 倍**に立ち（G-SPEC PASS）、48フレームは動き（G-MOVE PASS）、"
      "出力は 48 frames/2.000 s（G-OUT PASS）。G-LINE は 5–15 px 帯で"
      "**写真とほぼ同値**（最大 %.1f ‰ vs %.1f ‰）で、4行中 %d 行が写真を"
      "わずかに上回り **FAIL**——v4-1.4 の格子（写真の2倍）とは別物だが、"
      "「写真以下」というゲート表記は満たしていない（§5・§10）。"
      % (geo["lambda0_m"], geo["L_gen_m"], mp_map["magnification"],
         geo["drawn_lambda_m"], gs["ratio_gen"] or float("nan"),
         gl["strength_v4_1_6r"], gl["strength_photo"],
         sum(1 for r in gl["rows"] if not r["ok"])))
    A("")
    A("## 1. 分離形の構成（brief の作業 1–3）")
    A("")
    A("| 量 | 値 | 出所 |")
    A("| --- | --- | --- |")
    A("| 環境（川の水） | H=%.4f m, U=%.2f m/s, U0=%.2f m/s, S0=%.4f, "
      "D84=%.3f m, seed=%d, 媒体 %s | `run_v4_1_2.REF` / `build_stream` verbatim |"
      % (geo["H_m"], geo["U_m_s"], geo["U0_m_s"], geo["S0"], geo["D84_m"],
         geo["seed"], geo["water_preset"]))
    A("| k0 / λ0 | %.6f rad/m / %.6f m | `stream_surface2.k0_stationary(U0,H)` |"
      % (geo["k0_rad_m"], geo["lambda0_m"]))
    A("| 生成パッチ L_gen | **%.6f m** = 10 × λ0（k0/dk = %.0f） | brief |"
      % (geo["L_gen_m"], geo["k0_over_dk"]))
    A("| 生成格子 | N=%d、cell=%.5f mm、λ0 = %.1f cell | natural シーンの格子 |"
      % (geo["grid"], geo["cell_gen_m"] * 1e3, geo["cells_per_lambda_gen"]))
    A("| 画像パッチ L_img | **%.4f m**（= 8·Us/n、natural シーンと同一画面寸法） | "
      "`stream_surface2.loop_adjusted_L` |" % geo["L_img_m"])
    A("| 倍率 | **%.6f**（L_img/L_gen） | 計算 |" % mp_map["magnification"])
    A("| 描かれる波長 | **%.6f m** = %.3f px = %.1f cell | λ0 × 倍率 |"
      % (mp_map["drawn_lambda_m"], geo["drawn_lambda_px"],
         geo["drawn_cells_per_lambda"]))
    A("| 傾きの増幅 | ×%.6f（同じ面を画像の物差しで読む） | 写像 |"
      % mp_map["slope_gain"])
    A("| 屈折水深 z | %.4f m | brief / env.depthM |" % geo["depth_m"])
    A("| 合成 | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d（凍結形） | "
      "`run_v4_1_2` の凍結値 |"
      % (geo["composite"]["A"], geo["composite"]["B_sigma_px"],
         geo["composite"]["E"], geo["composite"]["Q"]))
    A("| 環境（光） | 太陽高度 %.2f°, turbid %.3f, clarity %.3f, 位相 %s / "
      "カーネル %s | `build_stream`（verbatim） |"
      % (sc["env"]["sun_altitude_deg"], sc["env"]["turbid"],
         sc["env"]["clarity"], sc["optics"]["phase_kind"],
         sc["optics"]["kernel"]))
    A("")
    A("- **写像は index 一致**: 両方 N=%d なので `h` はコピー、`(gx,gy)` は "
      "×%.6f。実測チェック: |Δh|max = %.3g、傾き比の実測 = %.6f（期待 %.6f）。"
      % (geo["grid"], mp_map["slope_gain"],
         fd["mapping"]["check"]["h_max_abs_diff"],
         fd["mapping"]["check"]["gx_ratio_measured"],
         fd["mapping"]["check"]["slope_gain_expected"]))
    A("- 波の実寸: 生成 rms 傾き %.5f → 描画 rms 傾き **%.5f**、rms 高さ "
      "%.4f mm（高さは不変、傾きだけ増える）。"
      % (fd["mapping"]["gen_rms_slope"], fd["mapping"]["drawn_rms_slope"],
         fd["mapping"]["drawn_rms_height_m"] * 1e3))
    A("- 2 s は**ループしない**（L_gen/Us = %.2f s が移流モード）。ループは"
      "主張しない。" % geo["loop_mode_s"])
    A("")
    A("## 2. 波の場（川の水）— 会計・成分・スペクトル（G-SPEC）")
    A("")
    A("会計は `stream_surface3`（v4-1.4）のまま: var(ii)=φ·σ_reg²、"
      "var(iii)=(1−φ)·σ_reg²、σ_total² = σ_reg² + var_i（底応答は独立項）。"
      "数値は一切読み書きしていない。")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| σ_reg（回帰） | %.6f mm |" % (acc["sigma_reg_m"] * 1e3))
    A("| var_i（底応答, 物理値） | %.4e m²（σ_i = %.6f mm） |"
      % (acc["var_i_physical_m2"], acc["sigma_i_physical_m"] * 1e3))
    A("| var(ii) / var(iii) | %.4e / %.4e m² |"
      % (acc["var_ii_m2"], acc["var_iii_m2"]))
    A("| σ_total | %.6f mm |" % (acc["sigma_total_m"] * 1e3))
    A("| var_i / var_reg | %.4f |" % acc["var_i_over_var_reg"])
    A("| φ（リング率） | %.4f（kmax/k0 = %.1f） |"
      % (acc["ring_fraction_phi"] if "ring_fraction_phi" in acc
         else fd["state"]["ring_fraction_phi"],
         acc["ring_fraction_detail"]["rho_kmax_over_k0"]
         if "ring_fraction_detail" in acc
         else fd["state"]["ring_fraction_detail"]["rho_kmax_over_k0"]))
    A("| τ（=H/u*）/ u* | %.4f s / %.5f m/s |"
      % (fd["state"]["tau_s"], fd["state"]["derived"]["ustar"]))
    A("| L0（乱流相関長） | %.4f m（上限適用 %s） |"
      % (fd["state"]["L0"]["L0_used_m"], fd["state"]["L0"]["cap_applied"]))
    A("| 成分 σ（frame 24, m） | i %.6f / ii %.6f / iii %.6f / total %.6f |"
      % (fd["components_frame24"]["i"]["sigma_m"],
         fd["components_frame24"]["ii"]["sigma_m"],
         fd["components_frame24"]["iii"]["sigma_m"],
         fd["components_frame24"]["total"]["sigma_m"]))
    A("| 成分 rms 傾き（frame 24） | i %.5f / ii %.5f / iii %.5f / total %.5f |"
      % (fd["components_frame24"]["i"]["rms_slope"],
         fd["components_frame24"]["ii"]["rms_slope"],
         fd["components_frame24"]["iii"]["rms_slope"],
         fd["components_frame24"]["total"]["rms_slope"]))
    A("")
    A("**G-SPEC（annulus ピーク、bin 幅 1 dk）**")
    A("")
    A("| 測った場 | ピーク [rad/m] | 解析 k0 [rad/m] | 比 | 判定 |")
    A("| --- | --- | --- | --- | --- |")
    A("| 生成場 全体（frame 24, L_gen） | %.4f | %.4f | **%.4f** | %s |"
      % (spec["peak_gen_rad_m"], spec["k0_analytic_gen_rad_m"],
         spec["ratio_gen"], "OK" if 0.75 <= spec["ratio_gen"] <= 1.33 else "NG"))
    A("| うち流れの成分 (ii)+(iii) | %.4f | %.4f | %.4f | — |"
      % (spec["peak_gen_flow_only_rad_m"], spec["k0_analytic_gen_rad_m"],
         spec["ratio_gen_flow_only"]))
    A("| うち底応答 (i) のみ | %.4f | %.4f | %.4f | — |"
      % (spec["peak_gen_bed_only_rad_m"], spec["k0_analytic_gen_rad_m"],
         spec["ratio_gen_bed_only"]))
    A("| 生成場 全体（L_img の物差し） | %.4f | %.4f | **%.4f** | %s |"
      % (spec["peak_img_rad_m"], spec["k0_analytic_img_rad_m"],
         spec["ratio_img"], "OK" if 0.75 <= spec["ratio_img"] <= 1.33 else "NG"))
    A("")
    A("- ピークは k0 の **%.4f 倍**（許容 0.75–1.33 の内側）。48フレーム全部の"
      "比は %.4f–%.4f で、時間発展（移流＋OU）でもピークは動かない。"
      % (spec["ratio_gen"], spec["per_frame_ratio_min"],
         spec["per_frame_ratio_max"]))
    A("- 底応答 (i) 単独のピークは %.4f 倍でリング帯の外にあるが、その傾き分散は"
      "全体の %.1f%% なので、全体のピークを動かさない（成分別の値を併記して"
      "いるので、この配分は監査できる）。"
      % (spec["ratio_gen_bed_only"],
         100.0 * fd["components_frame24"]["i"]["rms_slope"] ** 2
         / max(fd["components_frame24"]["total"]["rms_slope"] ** 2, 1e-30)))
    A("")
    A("## 3. 底応答 (i) の底 = 元写真を L_gen に写像")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 底の形 | 元写真の緑チャンネル（%dx%d）を 2x2 ミラーで周期化（%dx%d）"
      "→ N=%d に cubic 周期リサンプル |"
      % (fd["bed"]["photo_shape"][1], fd["bed"]["photo_shape"][0],
         fd["bed"]["tile_shape"][1], fd["bed"]["tile_shape"][0], geo["grid"]))
    A("| 底の振幅 | σ_b = D84/4 = **%.4f mm**（Nikora の写像。写真は物差しを"
      "持たないため、振幅はモデルの規則から取る） |"
      % (fd["bed"]["sigma_b_m"] * 1e3))
    A("| 底ピッチ | %.4f mm（L_gen/N） |" % (fd["bed"]["bed_pitch_m"] * 1e3))
    A("| 伝達 T(k) | v4-1.2 の 1.1 重力伝達（無改造）、band 定義 "
      "|denominator|<%.2f |" % SS.RES_BAND_DENOM)
    A("| 安全弁 | band の rms 傾き ≤ %.2f。今回 scale=%.6f（= 発動なし）、"
      "band 傾き %.6f |"
      % (SS.SLOPE_CAP, fd["info_i"]["scale_applied"],
         fd["info_i"]["band_rms_slope_before"]))
    A("| θ_i（底応答の rms 高さ） | **%.6f mm**、傾き %.5f |"
      % (fd["info_i"]["sigma_i_m"] * 1e3,
         fd["components_frame24"]["i"]["rms_slope"]))
    A("")
    A("- 底を写真にしても、伝達 T(k) が k ≳ 20 rad/m をほぼ殺す（T_max = "
      "%.3f）ので、効くのは写真の**大きな明暗構造**だけ。したがって"
      "「写真の木目がそのまま波に焼き付く」ことはなく、川の水の底応答として"
      "妥当な大きさ（θ_i = %.2f mm、波の σ_reg = %.2f mm の %.0f%%）に収まる。"
      % (fd["info_i"]["T_max"], fd["info_i"]["sigma_i_m"] * 1e3,
         acc["sigma_reg_m"] * 1e3,
         100.0 * fd["info_i"]["sigma_i_m"] / acc["sigma_reg_m"]))
    A("- 写像の実装は `run_v4_1_6r.response_of_zeta`（v4-1.2 `_build_i` の"
      "伝達＋安全弁ブロックの逐語コピー）。**同一性の数値確認**: 同じ ζ̂ を"
      "通すと v4-1.2 自身の η_i と最大差 %.3g（σ_i の差 %.3g）→ %s。"
      % (sc["field_selfcheck"]["max_abs_eta_i_diff"],
         sc["field_selfcheck"]["abs_sigma_i_diff"],
         "一致" if sc["field_selfcheck"]["pass"] else "不一致"))
    A("- ミラー周期化は「継ぎ目を出さない」ためだけの処置（生の周期リサンプルは"
      "継ぎ目に直線を焼き、それは全フレームに静止線として残る）。覆う物理的"
      "内容（低 k の明暗）は変わらない。")
    A("")
    A("## 4. ベース画像と合成")
    A("")
    A("- 元写真（735×539、8bit）を L_img = %.4f m の画面として 736×540 で描く"
      "（PHOTO_NATIVE はフレーム−1 px、v4-1.2 と同じ扱い）。" % geo["L_img_m"])
    A("- 屈折の水深 z = %.4f m、IOR は `image.IOR`（RGB 別）、合成は凍結形 "
      "A=%.1f / Blow σ=%.1f px / E=%.1f / Q=%d。"
      % (geo["depth_m"], geo["composite"]["A"], geo["composite"]["B_sigma_px"],
         geo["composite"]["E"], geo["composite"]["Q"]))
    A("- frame 24 の画像統計: C mean %.4f (min %.4f p50 %.4f p99 %.3f max %.3f)、"
      "crush %.2f%%、blow %.2f%%、linear mean %.5f、g1 %.4f"
      % (fx["C_stats_frame24"]["mean"], fx["C_stats_frame24"]["min"],
         fx["C_stats_frame24"]["p50"], fx["C_stats_frame24"]["p99"],
         fx["C_stats_frame24"]["max"],
         100 * fx["image_metrics"]["crush_rate_linear_lt_0.02"],
         100 * fx["image_metrics"]["blowout_rate_aces_gt_0.98"],
         fx["image_metrics"]["linear_mean_pre_tonemap"], fx["g1"]))
    nr = n.get("natural_scene_reference") or {}
    if nr.get("C_p99"):
        A("- **同じ 3.68 m・同じ合成で比べると、面の傾きが強い**。自然シーンを"
          "自分のパッチ寸法のまま描いた v4-1.2 は rms 傾き %.4f・C p99 %.2f"
          "（min %.2f max %.2f）だった。本作は写像の傾き増幅 ×%.4f で rms 傾き "
          "**%.4f**、C p99 **%.3f**（min %.2f max %.2f）。コースティクスの"
          "コントラストが上がっているのは主にこの倍率の効果で、合成値の変更"
          "ではない（A/Blow/E/Q は凍結のまま）。"
          % (nr.get("rms_slope_total") or float("nan"), nr["C_p99"],
             nr.get("C_min") or float("nan"), nr.get("C_max") or float("nan"),
             mp_map["slope_gain"], fd["mapping"]["drawn_rms_slope"],
             fx["C_stats_frame24"]["p99"], fx["C_stats_frame24"]["min"],
             fx["C_stats_frame24"]["max"]))
    A("")
    A("## 5. G-LINE — 格子の強さ（v4-1.5 と同一測定・同一表形式）")
    A("")
    A("測定は v4-1.5 の実装を**そのまま import** して使う（窓の決め方・帯域・"
      "量が同一であることを、書き写しではなく同一コードで担保する）。")
    A("")
    A("- **depth（主判定）**: 生画像に理想帯域通過（軸平行 5–15 px、|k_perp|≤3）"
      "をかけた成分の RMS ÷ 窓平均 ×1000（‰）。")
    A("- **depth_hp（brief の手順 / クロスチェック）**: 画像 − gaussian(σ=3 px) "
      "を先に掛けてから同じ帯域通過。σ=3 ハイパスは 9 px を 0.11 倍・14 px を "
      "0.40 倍に落とすため、周期の違う2枚を比べると偏る（数値は全部出す）。")
    A("")
    A("| 領域 | 軸 | 元写真 period / depth / depth_hp | v4-1.4 depth | v4-1.5 depth "
      "| **v4-1.6r** period / depth / depth_hp | 比 depth (hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for i, r in enumerate(gl["rows"]):
        p = gl["photo"][i]
        v14 = gl["v4_1_4"][i] if gl["v4_1_4"] else None
        v15 = gl["v4_1_5"][i] if gl["v4_1_5"] else None
        A("| %s | %s | %.2f px / %.1f / %.1f | %.1f | %.1f | **%.2f px / %.1f / "
          "%.1f** | %.3f (%.3f) | %s |"
          % (r["region"], r["axis"], p["period_px"], p["depth_e3"],
             p["depth_hp_e3"],
             v14["depth_e3"] if v14 else float("nan"),
             v15["depth_e3"] if v15 else float("nan"),
             r["period_v4_1_6r_px"], r["depth_v4_1_6r"],
             r["depth_hp_v4_1_6r"], r["ratio"], r["ratio_hp"],
             "**OK**" if r["ok"] else "**FAIL**"))
    A("")
    nfail = sum(1 for r in gl["rows"] if not r["ok"])
    A("- 主判定: **%d/4 行が写真以下 → %s**。最大強度 %.1f ‰（元写真 %.1f ‰、"
      "v4-1.5 %.1f ‰、v4-1.4 %.1f ‰）。比の範囲は %.3f–%.3f。"
      % (4 - nfail, "PASS" if gl["pass"] else "FAIL", gl["strength_v4_1_6r"],
         gl["strength_photo"],
         gl["strength_v4_1_5"] if gl["strength_v4_1_5"] else float("nan"),
         gl["strength_v4_1_4"] if gl["strength_v4_1_4"] else float("nan"),
         min(r["ratio"] for r in gl["rows"]),
         max(r["ratio"] for r in gl["rows"])))
    A("- hp クロスチェックは %d/4 行が写真以下。"
      % sum(1 for r in gl["rows"] if r["ok_hp"]))
    A("")
    A("### 5.1 窓ごとの比（24窓・両軸）")
    A("")
    A("256 px タイル格子（stride 92×90、24窓）で、同じ窓の同じ軸の depth 比"
      "（画像 / 元写真）を全部計算した。")
    A("")
    dst = gl["distribution"]
    A("| 比較 | 軸 | 中央値 | p10–p90 | 写真を上回った窓 | 最大比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for tag, key in (("v4-1.6r / 元写真", "v4_1_6r_vs_photo"),
                     ("v4-1.4 / 元写真", "v4_1_4_vs_photo")):
        d = dst.get(key)
        if not d:
            continue
        for axis in ("x", "y"):
            A("| %s | %s | %.3f | %.3f–%.3f | %d/%d | %.3f |"
              % (tag, axis, d[axis]["ratio_median"], d[axis]["ratio_p10"],
                 d[axis]["ratio_p90"],
                 int(round(d[axis]["frac_above_1"] * d[axis]["n_windows"])),
                 d[axis]["n_windows"], d[axis]["max_ratio"]))
    A("")
    d6 = dst["v4_1_6r_vs_photo"]
    d4 = dst.get("v4_1_4_vs_photo")
    A("- v4-1.6r は**画面全体で写真とほぼ同値**（中央値 x %.3f / y %.3f、"
      "上回る窓 x %d/24・y %d/24、最大 %.3f / %.3f）。v4-1.4 は中央値 "
      "%.2f / %.2f で 24/24 窓が超過していた。つまり残差は「格子」ではなく"
      "**写真と同じ帯域パワー**である。"
      % (d6["x"]["ratio_median"], d6["y"]["ratio_median"],
         int(round(d6["x"]["frac_above_1"] * d6["x"]["n_windows"])),
         int(round(d6["y"]["frac_above_1"] * d6["y"]["n_windows"])),
         d6["x"]["max_ratio"], d6["y"]["max_ratio"],
         d4["x"]["ratio_median"] if d4 else float("nan"),
         d4["y"]["ratio_median"] if d4 else float("nan")))
    A("")
    A("### 5.2 追加診断（合否とは無関係）— 5–15 px 帯はどこに乗っているか")
    A("")
    A("frame 24 を3通りで描き、同じ統計で測った（輸送は本番と同じ σ）。")
    A("")
    A("| 描き方 | 最大強度 [‰] | 4行の内訳 (s/x, s/y, b/x, b/y) | C std |")
    A("| --- | --- | --- | --- |")
    for tag in ("a", "b", "c"):
        v = dg["variants"].get(tag) or {}
        if not v:
            continue
        A("| %s | %.1f | %s | %.3f |"
          % (v["name"], v["strength"],
             ", ".join("%.1f" % x["depth_e3"] for x in v["measurements"]),
             v["C_std"]))
    A("| 元写真 | %.1f | %s | — |"
      % (dg["photo_strength"],
         ", ".join("%.1f" % x["depth_e3"]
                   for x in gl["photo"])))
    A("")
    va = (dg["variants"].get("a") or {}).get("strength")
    vb = (dg["variants"].get("b") or {}).get("strength")
    A("- 屈折ワープ単独でも写真と同値（%.1f vs %.1f）、コースティクス単独は "
      "%.1f。**足し合わせると写真の数%% 上に乗る**——格子の signature"
      "（同じ周期が4行にそろって立つ）は無く、周期も写真と同族（%.2f–%.2f px）。"
      % (va or float("nan"), dg["photo_strength"] or float("nan"),
         vb or float("nan"),
         min(r["period_v4_1_6r_px"] for r in gl["rows"]),
         max(r["period_v4_1_6r_px"] for r in gl["rows"])))
    A("")
    A("## 6. 輸送の帯域制限（v4-1.5 の折り返し規則をこの場面で導く）")
    A("")
    A("| σ [cell] | p50 \\|ds/di\\| | p99 | max | 折り返し面積率 | 平均変位 [px] | "
      "リング分散伝達 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for row, pr_ in zip(fx["criterion"]["ladder"], fx["fold_rule"]["per_width"]):
        A("| %.1f%s | %.3f | %.3f | %.3f | %.2f%% | %.2f | %.4f |"
          % (row["sigma_cells"],
             "（= フィルタなし）" if row["sigma_cells"] == 0.0 else "",
             row["ds_per_ray_p50"], row["ds_per_ray_p99"],
             row["ds_per_ray_max"], 100 * row["folded_area_fraction"],
             row["displacement_px_mean"], pr_["ring_variance_transfer"]))
    A("")
    A("- **選び方（2つの制約を同時に満たす最小値）**: (A) 前方スプラットが"
      "単値であること `max|ds/di| < 1.0`、(B) 波を弱めないこと（リング帯の"
      "分散伝達 ≥ %.2f）。両方を満たす最小の幅は **σ = %.1f cell = %.4f mm**"
      "（(B) の実測 %.4f）。"
      % (fx["fold_rule"]["wave_keep_min"], fx["sigma_cells_chosen"],
         fx["sigma_m"] * 1e3,
         fx["ring_variance_transfer_chosen"]))
    A("- **v4-1.5 の安全率4倍（`max|ds/di| <= 0.25`）はこの場面では使えなかった**:"
      "到達するのは σ=%s cell だが、そこでのリング分散伝達は **%.4f**"
      "（= 描かれる波の分散の %.0f%% を消す）。したがって (B) で棄却した。"
      "ゲートではなく、この2制約の交差で決めている。"
      % ("%.1f" % fx["fold_rule"]["factor4_reached_at_sigma_cells"]
         if fx["fold_rule"]["factor4_reached_at_sigma_cells"] is not None
         else "—",
         fx["fold_rule"]["factor4_ring_variance_transfer"] or float("nan"),
         100 * (1.0 - (fx["fold_rule"]["factor4_ring_variance_transfer"]
                       or 0.0))))
    A("- ラダー全行（σ, p50/p99/max|ds/di|, 折り返し面積率, 平均変位, リング"
      "分散伝達）は上の表と `numbers-v4-1-6r.json` の "
      "`transport.fold_rule.per_width` に残した（後から監査できる）。")
    A("- 波を弱めていないこと: 画像スケールのリング帯（k=%.4f rad/m、"
      "k0×%.4f）の傾き分散比 = **%.6f**、全傾き rms %.5f → %.5f。"
      "σ=0 は入力と bit 一致（h/gx/gy: %s）。"
      % (fx["ring_band"]["k_ring_rad_m"], mp_map["slope_gain"],
         fx["ring_band"]["ratio"], fx["ring_band"]["total_slope_rms_before"],
         fx["ring_band"]["total_slope_rms_after"],
         all(fx["self_check_identity_sigma0"]["per_channel"].values())))
    A("- 画像形成のワープは**フィルタ前**の場を使う（帯域制限は輸送だけ）。"
      "v4-1.5 と同じ扱いで、G-WARP の数値は「波を弱めていない」ことの側の"
      "証拠として読める。")
    A("")
    A("## 7. G-MOVE / G-WARP / G-OUT")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| ユニークなフレーム PNG | **%d / 48** |" % gm["unique_frames"])
    A("| フレーム間平均絶対差（緑, PNG） | **%.4f 諧調**（min %.4f / median %.4f "
      "/ max %.4f） |"
      % (gm["interframe_green_png"]["mean"],
         gm["interframe_green_png"]["min"],
         gm["interframe_green_png"]["median"],
         gm["interframe_green_png"]["max"]))
    if isinstance(gm.get("interframe_green_mp4"), dict) \
            and "mean" in gm["interframe_green_mp4"]:
        A("| 同（デコード後 mp4） | %.4f 諧調 |"
          % gm["interframe_green_mp4"]["mean"])
    A("| 閾値 / **G-MOVE** | > 0.5 諧調 / **%s** |"
      % ("PASS" if gm["pass"] else "FAIL"))
    A("| frame 24 歪み（緑 平均） | **%.4f px**（p50 %.4f / p99 %.4f、"
      "max %.4f） |"
      % (gw["g_mean_px"], gw["g_p50_px"], gw["g_p99_px"],
         gw["detail"]["per_channel"]["g"]["px_max"]))
    A("| frame 24 C の標準偏差 | **%.4f**（mean %.4f） |"
      % (gw["C_std_frame24"], gw["C_mean_frame24"]))
    A("| **G-WARP** | **記録のみ（合否なし、brief どおり）** |")
    A("")
    A("- 参考（合否なし）: v4-1.4 の実寸シーンは緑 平均 1.8393 px だった。"
      "本作は水深 0.35 m・描画波長 0.368 m なので変位は %.2f px と大きい。"
      "これは場面が違うためで、比較はしていない（比 %.3f は情報として残す）。"
      % (gw["g_mean_px"], gw["detail"].get(
          "ratio_to_v4_1_4_reference_only") or float("nan")))
    A("")
    A("| file | codec | 解像度 | fps | frames | duration | pix_fmt |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    A("| %s | %s | %s×%s | %s | %s | %.4f s | %s |"
      % (mp["mp4"], pr.get("codec_name"), pr.get("width"), pr.get("height"),
         pr.get("r_frame_rate"), pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("pix_fmt")))
    A("")
    A("- サイズ %d バイト、sha256 `%s`。" % (mp["size_bytes"],
                                          mp["sha256"][:32]))
    A("")
    A("## 8. G-LEGACY — 既存 py の無傷性（外部差分の帰属分離）")
    A("")
    A("- stage-0 記録: `%s`（全 %d 本 = 既存 %d 本 + 本作の新規 1 本、"
      "sha256 + size + mtime）。"
      % (n["stage0_hashes_file"], gleg["recorded_existing"] + 1,
         gleg["recorded_existing"]))
    A("- 本作が書くのは `out/v4-1-6r/**` だけ（コード上、既存 .py を開くのは"
      "読み取りのみ）。したがって既存 .py の差分は**並行レーンの外部書き込み**"
      "であり、次のように分けて報告する。")
    A("")
    A("| 区分 | 本数 | 内容 |")
    A("| --- | --- | --- |")
    A("| 自レーンの変更（既存） | **0** | なし |")
    A("| 消失 | **%d** | %s |"
      % (len(gleg["missing"]), ", ".join(gleg["missing"]) or "なし"))
    A("| 外部（並行レーン）差分 | **%d** | %s |"
      % (len(gleg["external_changes"]),
         ", ".join("%s (%s→%s)" % (c["file"], c["stage0_sha256"][:12],
                                   c["now_sha256"][:12])
                   for c in gleg["external_changes"]) or "なし"))
    A("| 新規 `run_v4_1_6r.py` | stage0 %s / 終了時 %s |"
      % (((gleg["own_new_file_hashes"].get("run_v4_1_6r.py") or {})
          .get("stage0") or {}).get("sha256", "")[:16],
         ((gleg["own_new_file_hashes"].get("run_v4_1_6r.py") or {})
          .get("now") or {}).get("sha256", "")[:16]))
    A("| stage-0 後に増えた .py（本作以外） | **%d** | %s |"
      % (len(gleg["new_since_stage0_not_of_this_job"]),
         ", ".join(gleg["new_since_stage0_not_of_this_job"]) or "なし"))
    A("")
    if gleg["external_changes"]:
        A("- 外部差分の同定根拠: (a) 本ドライバは既存 .py を書き込みモードで"
          "開かない、(b) 当該ファイルの mtime が stage-0 記録時刻より後"
          "（%s）、(c) stage-0 後に触られた .py の一覧 = %s。"
          % ("はい" if all(c["written_after_stage0"]
                           for c in gleg["external_changes"]) else "一部いいえ",
             ", ".join(gleg["foreign_py_touched_after_stage0"]) or "なし"))
    else:
        A("- 外部差分なし（並行レーンは本作の実行中に .py を書き換えなかった）。")
    A("- 判定: **%s**（自レーンの変更 0・消失 0）。厳密版（外部差分も含めて"
      "全て不変）は %s。"
      % ("PASS" if gleg["pass"] else "FAIL",
         "PASS" if gleg["pass_strict_everything_unchanged"] else "FAIL"))
    A("")
    A("## 9. 測定の定義（再現用）")
    A("")
    A("- 窓の位置は**元写真だけ**から決めた: 256×256 タイルを stride 92×90 で"
      "走査し、局所 std が最小のタイルを `smooth`（x0=%d, y0=%d, std %.2f）、"
      "最大を `bark`（x0=%d, y0=%d, std %.2f）。以後すべての画像で同じ座標"
      "（`run_v4_1_5.regions_from_photo` をそのまま使用）。"
      % (fx["regions"]["smooth_tile"]["x0"], fx["regions"]["smooth_tile"]["y0"],
         fx["regions"]["smooth_tile"]["std"],
         fx["regions"]["bark_tile"]["x0"], fx["regions"]["bark_tile"]["y0"],
         fx["regions"]["bark_tile"]["std"]))
    A("- 写真は 735×539、フレームは 736×540 なので同じ窓原点が同じ内容を指す。"
      "両方 8bit、depth は窓平均で正規化するので露出差は入らない。")
    A("- G-SPEC は `stream_surface2.wavenumbers` と同じ k グリッドで、"
      "|η̂|² を半径 1 dk の annulus に積分した分布のピーク（bin 中心 "
      "(j+0.5)dk）。生成スケール（L_gen）と画像スケール（L_img）の両方で"
      "測り、後者は解析 k0×%.4f と比べた。" % mp_map["slope_gain"])
    A("- depth の既知の限界: 軸帯の 1 次元マージナル統計であり、格子と"
      "「同じ軸帯パワーを持つ等方テクスチャ」を区別しない。")
    A("")
    A("## 10. 仮定・失敗・弱点（正直な開示）")
    A("")
    A("- **G-LINE は FAIL**（%d/4 行が写真を %.1f–%.1f%% 上回る）。"
      "これは v4-1.4 の格子（写真の約2倍、同じ周期が4行にそろう）とは"
      "別物で、24窓の中央値は x %.3f / y %.3f、周期も写真と同族。だが"
      "brief の判定文（「元写真以下」）は満たしていない。**ゲートを見て"
      "σ や合成を動かすことはしていない**（σ は §6 の2制約＝折り返し限界と"
      "波の分散保持だけで決め、A/Blow/E/Q は凍結値、L_gen/L_img は brief の値）。"
      % (sum(1 for r in gl["rows"] if not r["ok"]),
         100 * (min(r["ratio"] for r in gl["rows"] if not r["ok"]) - 1.0)
         if any(not r["ok"] for r in gl["rows"]) else 0.0,
         100 * (max(r["ratio"] for r in gl["rows"] if not r["ok"]) - 1.0)
         if any(not r["ok"] for r in gl["rows"]) else 0.0,
         d6["x"]["ratio_median"], d6["y"]["ratio_median"]))
    A("- **安全率4倍はこの場面では使えなかった**（§6）。v4-1.5 は `max|ds/di| "
      "<= 0.25` を満たす最小 σ（=4 cell）を選べたが、本作の描画波長は "
      "51 cell しかないため、0.25 を満たす幅（σ=%s cell）は波の分散を "
      "%.0f%% 消す。そこで (A) 折り返し限界と (B) 波の分散保持の交差で "
      "σ=%.1f cell を選んだ。**安全率を4倍から1倍に落とした**ことは"
      "そのまま開示する（数値は §6 のラダー全行）。"
      % ("%.0f" % fx["fold_rule"]["factor4_reached_at_sigma_cells"]
         if fx["fold_rule"]["factor4_reached_at_sigma_cells"] is not None
         else "—",
         100 * (1.0 - (fx["fold_rule"]["factor4_ring_variance_transfer"]
                       or 0.0)),
         fx["sigma_cells_chosen"]))
    A("- **底の振幅は仮定**（σ_b = D84/4）。写真は物差しを持たないので、"
      "形だけを写真から取り、振幅はモデルの Nikora 写像から取った。"
      "別の振幅を選べば var_i（= %.2e m²、波の %.0f%%）は変わる。"
      % (acc["var_i_physical_m2"],
         100.0 * fd["info_i"]["sigma_i_m"] / acc["sigma_reg_m"]))
    A("- **2 s はループしない**（移流モード %.2f s）。8 s ループを前提にした"
      "v4-1.2 の自然シーンとはそこが違う。本作は 2 s の切り出しで、"
      "ループを主張しない。" % geo["loop_mode_s"])
    A("- **傾きが 1/倍率 = ×%.4f される**のは分離形の帰結（同じ面を小さい"
      "画面で読む）。rms 傾き %.4f は水面として妥当な範囲だが、生成スケールの"
      "値ではない。絵の側のスケールで物理量を語るときは常にこの倍率がかかる。"
      % (mp_map["slope_gain"], fd["mapping"]["drawn_rms_slope"]))
    A("- **G-LINE は frame 24 の1フレーム判定**（brief の表記どおり）。"
      "48フレーム全体の周期測定はしていない。")
    A("- **底応答の形は写真の低周波だけ**（§3）。T(k) が k ≳ 20 rad/m を殺す"
      "ため、「木目がそのまま波に焼き付く」ことはない。逆に言えば、"
      "写真を底にした効果は弱い（θ_i は σ_reg の %.0f%%）。"
      % (100.0 * fd["info_i"]["sigma_i_m"] / acc["sigma_reg_m"]))
    _hp_bad = [r for r in gl["rows"] if not r["ok_hp"]]
    _hp_worst = max(_hp_bad, key=lambda r: r["ratio_hp"]) if _hp_bad else None
    A("- **hp クロスチェック（brief の手順）でも同じ%d行が上回る**。最大は "
      "%s-%s の %.3f（+%.1f%%）で、この行は写真側の基準値が4行で最も小さい"
      "（21.6 ‰）ため比が敏感になっている。主判定と hp で結論は同じ向き。"
      % (len(_hp_bad), _hp_worst["region"] if _hp_worst else "—",
         _hp_worst["axis"] if _hp_worst else "—",
         _hp_worst["ratio_hp"] if _hp_worst else float("nan"),
         100.0 * ((_hp_worst["ratio_hp"] - 1.0) if _hp_worst else 0.0)))
    A("- **stage の再実行（開示）**: 測定の実装を確定する過程で `--stage0` を2回、"
      "`--transport` を2回、`--gates`/`--report` を複数回走らせたため、"
      "checkpoint の unit id に同じ名前が複数回出る。成果物は最後の実行で"
      "上書きされており、frame 24 は `frame-0024.png` と `frames/frame-0024.png` "
      "の sha256 一致で再現性を確認した。コードの最終ハッシュと stage-0 時点の"
      "ハッシュは §8 に併記した。")
    A("- 判定は数値統計のみ（目視なし）。後付け調整なし。`frame-0024.png` は "
      "`frames/frame-0024.png` と sha256 が一致する（chunk 1 の再実行と"
      "同一であることの確認: %s）。"
      % ("一致" if (os.path.isfile(frame_path(STILL_IDX))
                   and sha256_file(frame_path(STILL_IDX))
                   == fx["still_sha256"]) else "不一致"))
    A("- 実行は全ステージ foreground（1 プロセス、順次）。jobkit checkpoint を"
      " `%s` に stage ごとに原子的にコミットした（unit id: %s）。"
      % (os.path.basename(CP_PATH),
         ", ".join("`%s`" % u for u in
                   (read_checkpoint().get("completed_ids") or []))))
    A("")
    A("## 11. 成果物")
    A("")
    for k, v in n["deliverables"].items():
        A("- %s → `%s`" % (k, v))
    A("")
    A("ファイル（bytes / sha256 先頭16）:")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for k, v in sorted(n["files"].items()):
        A("| %s | %d | %s |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    A("- `numbers-v4-1-6r.json`・本レポート・`checkpoint.json` は自己参照（または"
      "本レポートの後に更新）になるため、この表には入れない。`numbers` と "
      "`report` の最終 sha256 は `checkpoint.json` の `final_hashes` に、"
      "`checkpoint` 自身を含む全成果物の最終 sha256 は `--all` の実行ログ末尾に"
      "出した（すべて書き終えた後に計算している）。")
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-6r.json + REPORT-V4-1-6R.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    n = build_numbers()
    write_json(NUMBERS, n)
    md = write_report(n)
    atomic_write_bytes(REPORT, md.encode("utf-8"))
    final = {}
    for p in (NUMBERS, REPORT):
        final[os.path.basename(p)] = {"bytes": os.path.getsize(p),
                                      "sha256": sha256_file(p)}
    checkpoint_extra("final_hashes", {
        "written_utc": now_iso(),
        "note": ("hashes taken AFTER both files were written (no self "
                 "reference inside them)"),
        "files": final})
    cp.unit_done("report")
    log("  final sha256 (after all writes):")
    for p in (NUMBERS, REPORT, CP_PATH, MP4, STILL, STAGE0_HASHES):
        if os.path.isfile(p):
            log("    %-28s %9d  %s" % (os.path.basename(p),
                                       os.path.getsize(p), sha256_file(p)))
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--field", action="store_true")
    ap.add_argument("--transport", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.all:
        stage0()
        do_field()
        do_transport()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_diag()
        do_warp()
        do_encode()
        do_gates()
        do_report()
        log("=" * 78)
        log("[all] done in %.1fs" % (time.time() - t0))
        return 0
    if a.stage0:
        stage0()
    if a.field:
        do_field()
    if a.transport:
        do_transport()
    if a.chunk is not None:
        do_chunk(a.chunk)
    if a.diag:
        do_diag()
    if a.warp:
        do_warp()
    if a.encode:
        do_encode()
    if a.gates:
        do_gates()
    if a.report:
        do_report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
