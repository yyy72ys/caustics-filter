"""run_v4_1_15 -- BRIEF-V4-1-15: the area-conserving deposit at the two
remaining verification sizes.

NEW FILES ONLY: this driver.  `transport_mesh.py` (the deposit) is reused
verbatim; every pre-existing .py is left byte-identical (G-LEGACY).

WHAT THIS JOB ASKS
    v4-1.8 proved the area-conserving mesh deposit (kernel / jitter / band
    limit none) at the verification MAIN size (1.4755 x 1.0828 m window,
    depth 0.6 m).  The brief's three-size decision leaves two sizes:

        STAGE A  "contrast frame"  window 2.4 x 1.76 m, depth 0.6 m,
                 the SAME river scene as v4-1.8 (only the window changes)
        STAGE B  "sub frame"       window 4.8 x 3.52 m, depth 4.0 m,
                 preset 5 (deep clear pond, wind), U = 6, F = 100

    In both, "image follows the wave": the generated patch IS the screen
    window, so the magnification is 1 (no stretching exists).

DEPOSIT
    `transport_mesh.mesh_C` verbatim: every water cell -> quad -> 2
    triangles, corners refracted with det.det_frame's Snell step, the cell's
    energy distributed over the 736x540 frame in EXACT proportion to the
    analytic triangle/pixel overlap area.  No kernel, no jitter, no band
    limit, no mean normalisation.  w_cell = W*H/N_s^2 (the existing C lane's
    divisor convention), off-frame landings are wrapped by the periodic floor
    and counted in a separate book (never used to rescale).

ONE DOCUMENTED PERFORMANCE SWITCH (STAGE B only)
    transport_mesh's `BOX_BULK_MAX` decides which of its two rasterisers
    handles a projected triangle.  Both are exact analytic overlap; they
    differ only in summation order.  At depth 4.0 m the projected mesh is
    folded ~20-300x, so almost every triangle exceeds the default 16-px
    bbox and falls into the one-triangle-at-a-time path (~1.4 ms/triangle).
    For STAGE B the constant is raised at run time so the vectorised path
    handles everything.  The file is NOT modified; the equivalence of the two
    paths is measured on real triangles and recorded (stage 0).  STAGE A runs
    with the module's own constant untouched.

GRIDS (and the one place where this job departs from the brief's guide)
    Brief: the base grid's Nyquist wavenumber must cover a comparable
    physical band to the v4-1.8 main run (512^2 / 1.4755102 m = 347 cells/m),
    "目安 1024^2".

        STAGE A  base 1024^2 -> 427 cells/m  (covers the main run: OK)
                 half  512^2               (G-CONV)
        STAGE B  base  512^2 ->  53 cells/m  (does NOT cover the main band)
                 half  256^2               (G-CONV)
                 plus a supplementary 1024^2 frame-24 deposit (13 min) so the
                 convergence question can be answered upward as well.

    STAGE B cannot be run at the brief's guide (1024^2) or at the band rule's
    own value (2048^2) with the verbatim deposit: the measured cost is
    5.5 us per triangle bbox pixel and the folded projection needs
    1.45e8 / 5.6e8 box clips per frame, i.e. 13 min / 47 min per frame and
    10 h / 38 h for the 48-frame piece.  The measured cost curve is in the
    report; 512^2 is the largest grid that completes, and every stage-B number
    below carries that label.

STAGES (all foreground, jobkit checkpoint committed per stage)
    --stage0              absolute-path preflight, stage-0 sha256 record,
                          scene records, transport_mesh selftest, deposit-path
                          equivalence, measurement provenance
    --frame24 STAGE       frame 24: C, half-grid convergence, point-splat
                          reference on the SAME field, still, G-LINE pre-check
    --chunk STAGE K       frames 24K..24K+23
    --warp STAGE          displacement yardstick (G-WARP, record only)
    --encode STAGE        MP4 + probe (G-OUT)
    --gates               all gates, both stages
    --report              numbers-v4-1-15.json + REPORT-V4-1-15.md
    --all
"""
import argparse
import datetime
import gc
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v10 as R10  # noqa: E402  (load_recipe, reused UNCHANGED)
import run_v15 as RV15  # noqa: E402  (scene_for, reused UNCHANGED)
import run_v15_spectrum as V15S  # noqa: E402  (spectrum source / arm / warn)
import run_v4_1_2 as R12  # noqa: E402  (frozen constants, reused UNCHANGED)
import stream_surface2 as SS  # noqa: E402
import stream_surface3 as S3  # noqa: E402
import run_v16_water as V16  # noqa: E402
import v12_water as W  # noqa: E402
import transport_mesh as TM  # noqa: E402  (this job's deposit, VERBATIM)
from det import det_frame  # noqa: E402  (secondary reference only)
from state import derive_seed  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-15.md"
OUT = os.path.join(HERE, "out", "v4-1-15")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
NUMBERS = os.path.join(OUT, "numbers-v4-1-15.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-15.md")

PHOTO = R12.PHOTO
FPS = R12.FPS
WIDTH, HEIGHT = R12.WIDTH, R12.HEIGHT
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
NEW_FILES = ["run_v4_1_15.py"]

#: the reference grid of the halo kernel profile: set to each scene's own
#: state grid before rendering (run_v4_1_8 hardcoded 512 there because that was
#: its state grid -- same rule, generalised).  Declared with the copied
#: function below, where run_v4_1_7 declares it.
BUF_GRID_N = 512

STAGES = ("contrast", "sub")

SCN = {
    "contrast": {
        "label": "対比枠 (v20 size)",
        "kind": "stream",
        "L": 2.4, "depth": 0.60,
        "window_m": [2.4, 2.4 * HEIGHT / WIDTH],
        "N": 1024, "half": 512, "refine": None,
        "probe_grids": [512, 1024],
        # brief STAGE A item 1: v4-1.8's river, verbatim
        "U": 0.80, "S0": 0.002, "D84": 0.09, "seed": 20261112,
        "water": "algae_pond",
        "mp4": "piece-2s-v4-1-15-contrast.mp4",
        "still": "frame-0024-contrast.png",
        "frames": "frames-contrast",
    },
    "sub": {
        "label": "サブ枠 (v19 size)",
        "kind": "wind",
        "L": 4.8, "depth": 4.0,
        "window_m": [4.8, 4.8 * HEIGHT / WIDTH],
        "N": 512, "half": 256, "refine": 1024,
        "probe_grids": [256, 512, 1024, 2048],
        # brief STAGE B item 1: preset 5 wind cell, v19's own seed
        "U": 6.0, "F": 100.0, "seed": 3255045276382562994,
        "water": "clear",
        "mp4": "piece-2s-v4-1-15-sub.mp4",
        "still": "frame-0024-sub.png",
        "frames": "frames-sub",
    },
}

# ---- composite (frozen form, brief STAGE A item 3 / STAGE B item 3) -------
A_FROZEN = R12.A_FROZEN            # 2.6
B_SIGMA_FROZEN = R12.B_SIGMA_FROZEN  # 3.0
E_FROZEN = R12.E_FROZEN            # 2.2
Q_DET = R12.Q_DET                  # 25

# ---- measurement constants: run_v4_1_5's, verbatim ------------------------
WIN = 256
HP_SIGMA = 3.0
BAND_PX = (5.0, 15.0)
REGION_RULE = ("256x256 window; positions = tile grid (stride 92x90) of the "
               "SOURCE PHOTO; 'smooth' = tile with the lowest local pixel "
               "std, 'bark' = tile with the highest")
MEAS_FUNCS = ["green_of", "regions_from_photo", "lattice_metric",
              "lattice_strength", "measure_image", "tile_grid",
              "band_depth_window", "paired_distribution"]
ADAPTER_FUNCS = ["c_render", "render_optics"]

#: comparison stills / records of earlier lanes (G-LINE context only)
V17_STILL = os.path.join(HERE, "out", "v4-1-7", "frame-realscale-0024.png")
V17_NUMBERS = os.path.join(HERE, "out", "v4-1-7", "numbers-v4-1-7.json")
V18_STILL = os.path.join(HERE, "out", "v4-1-8", "frame-0024.png")
V41_8_NUMBERS = os.path.join(HERE, "out", "v4-1-8", "numbers-v4-1-8.json")


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


def _atomic_retry(path, data, verify=None, tries=40, wait=0.5):
    """atomic_write_bytes with a retry (a watcher process can hold a freshly
    written file for a moment on Windows).  Never used to hide a real error."""
    last = None
    for _k in range(tries):
        try:
            return atomic_write_bytes(path, data, verify=verify)
        except PermissionError as e:
            last = e
            time.sleep(wait)
    raise last


def write_json(path, obj):
    _atomic_retry(path, json.dumps(obj, ensure_ascii=False,
                                  indent=2).encode("utf-8"))


def makedirs():
    for tag in STAGES:
        os.makedirs(os.path.join(OUT, SCN[tag]["frames"]), exist_ok=True)
    os.makedirs(WORK, exist_ok=True)


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".py") and os.path.isfile(os.path.join(HERE, name)):
            out[name] = sha256_file(os.path.join(HERE, name))
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
        _atomic_retry(CP_PATH, json.dumps(
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

    _atomic_retry(path, data, verify=verify)
    return sha256_bytes(data)


def frame_path(tag, k):
    return os.path.join(OUT, SCN[tag]["frames"], "frame-%04d.png" % k)


def load_photo():
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"u8": padded, "lin": lin}


# ---------------------------------------------------------------------------
# the measurement: run_v4_1_5's implementation, verbatim
# ---------------------------------------------------------------------------
def green_of(path):
    from PIL import Image
    a = np.asarray(Image.open(path).convert("RGB"), dtype=np.float64)
    return a[:, :, 1]


def regions_from_photo(path=PHOTO, win=WIN):
    """Deterministic window positions: min / max local std on a tile grid."""
    g = green_of(path)
    H, W = g.shape
    tiles = []
    for y0 in range(0, H - win + 1, 90):
        for x0 in range(0, W - win + 1, 92):
            t = g[y0:y0 + win, x0:x0 + win]
            tiles.append((float(t.std()), float(t.mean()), int(x0), int(y0)))
    lo = min(tiles)
    hi = max(tiles)
    return [("smooth", lo[2], lo[3]), ("bark", hi[2], hi[3])], {
        "rule": REGION_RULE, "tiles_evaluated": len(tiles),
        "smooth_tile": {"x0": lo[2], "y0": lo[3], "std": lo[0], "mean": lo[1]},
        "bark_tile": {"x0": hi[2], "y0": hi[3], "std": hi[0], "mean": hi[1]},
        "photo_shape": [int(H), int(W)],
        "same_coordinates_for_every_image": True,
        "note": ("the photo is 735x539 and the frames are 736x540 "
                 "(PHOTO_NATIVE is the frame minus 1 px), so the same window "
                 "origin addresses the same content in both")}


def lattice_metric(g, regions, win=WIN, hp_sigma=HP_SIGMA, band=BAND_PX):
    """Axis-parallel 5-15 px strength, period and prominence.

    Two strength estimates are recorded for every row:

    depth_e3       amplitude of the ideal axis-parallel 5-15 px band-pass on
                   the RAW image (RMS / window mean * 1000).  This is the
                   primary gate quantity: the band-pass *is* the high pass, and
                   the estimate does not depend on where inside the band the
                   content sits.
    depth_hp_e3    the brief's literal procedure: first image - gaussian(3 px),
                   then the ideal band-pass.  Kept as a cross-check.  It is
                   frequency-weighted (the 3 px high pass passes 14 px content
                   3.6x better than 9 px content), so it is biased against
                   whichever image has the shorter period -- reported, not used
                   for the pass/fail.

    period_px / prom come from the brief's procedure (high pass, Hann window,
    2D FFT, axis rows) because that is the detection step of the symptom.
    """
    from scipy import ndimage as ndi
    hi = g - ndi.gaussian_filter(g, hp_sigma, mode="reflect")
    hann = np.outer(np.hanning(win), np.hanning(win))
    k = np.arange(1, win // 2)
    sel = (k >= win / band[1]) & (k <= win / band[0])
    kf = np.fft.fftfreq(win) * win
    KX, KY = np.meshgrid(kf, kf)
    rows = []
    for name, x0, y0 in regions:
        sub = hi[y0:y0 + win, x0:x0 + win]
        raw = g[y0:y0 + win, x0:x0 + win]
        mean_level = float(raw.mean())
        F = np.fft.fft2((sub - sub.mean()) * hann)
        P = np.abs(F) ** 2
        for axis in ("x", "y"):
            if axis == "x":
                band_p = P[np.ix_([0, 1, 2, win - 1, win - 2], k)].sum(axis=0)
            else:
                band_p = P[np.ix_(k, [0, 1, 2, win - 1, win - 2])].sum(axis=1)
            v = band_p[sel]
            kk = k[sel]
            j = int(np.argmax(v))
            a = np.abs(KX if axis == "x" else KY)
            o = np.abs(KY if axis == "x" else KX)
            m = (a >= win / band[1]) & (a <= win / band[0]) & (o <= 3)
            d_raw = 1000.0 * float(np.sqrt(np.mean(
                np.real(np.fft.ifft2(np.fft.fft2(raw - raw.mean()) * m)) ** 2))) \
                / mean_level
            d_hp = 1000.0 * float(np.sqrt(np.mean(
                np.real(np.fft.ifft2(np.fft.fft2(sub - sub.mean()) * m)) ** 2))) \
                / mean_level
            rows.append({"region": name, "axis": axis,
                         "period_px": float(win / kk[j]),
                         "k_bin": int(kk[j]),
                         "depth_e3": d_raw,
                         "depth_hp_e3": d_hp,
                         "prom_peak_over_median": float(v[j] / np.median(v)),
                         "window_mean_gray": mean_level,
                         "hp_sigma_px": float(hp_sigma),
                         "band_px": [float(band[0]), float(band[1])],
                         "window_px": int(win)})
    return rows


def lattice_strength(rows):
    return float(max(r["depth_e3"] for r in rows))


def measure_image(path, regions):
    return lattice_metric(green_of(path), regions)


def tile_grid(win=WIN, stride=(92, 90), shape=None):
    H, W = shape
    out = []
    for y0 in range(0, H - win + 1, stride[1]):
        for x0 in range(0, W - win + 1, stride[0]):
            out.append((x0, y0))
    return out


def band_depth_window(g, x0, y0, axis, win=WIN, band=BAND_PX, pre_hp=False):
    """Single-window axis band amplitude (same definition as depth_e3)."""
    from scipy import ndimage as ndi
    raw = g[y0:y0 + win, x0:x0 + win]
    sub = raw - ndi.gaussian_filter(g, HP_SIGMA, mode="reflect")[y0:y0 + win,
                                                              x0:x0 + win] \
        if pre_hp else raw
    kf = np.fft.fftfreq(win) * win
    KX, KY = np.meshgrid(kf, kf)
    a = np.abs(KX if axis == "x" else KY)
    o = np.abs(KY if axis == "x" else KX)
    m = (a >= win / band[1]) & (a <= win / band[0]) & (o <= 3)
    bp = np.real(np.fft.ifft2(np.fft.fft2(sub - sub.mean()) * m))
    return 1000.0 * float(np.sqrt(np.mean(bp ** 2))) / float(raw.mean())


def paired_distribution(path, ref_path, win=WIN):
    """Per-window depth ratio (image / source photo) over the whole tile grid.

    Robustness check for G-LINE: the gate itself is decided on the two
    designated regions (brief: wing + bark), this says how typical that is.
    """
    g = green_of(path)
    p = green_of(ref_path)
    H, W = p.shape
    tiles = tile_grid(win=win, shape=(H, W))
    out = {}
    for axis in ("x", "y"):
        ratios = []
        for x0, y0 in tiles:
            a = band_depth_window(g, x0, y0, axis, win=win)
            b = band_depth_window(p, x0, y0, axis, win=win)
            if b > 0:
                ratios.append(a / b)
        ratios = np.asarray(ratios, dtype=np.float64)
        out[axis] = {"n_windows": int(ratios.size),
                     "ratio_median": float(np.median(ratios)),
                     "ratio_p10": float(np.percentile(ratios, 10)),
                     "ratio_p90": float(np.percentile(ratios, 90)),
                     "frac_above_1": float(np.mean(ratios > 1.0)),
                     "max_ratio": float(ratios.max())}
    return out


# ---------------------------------------------------------------------------
# the deposit -> render adapter (run_v4_1_7's, verbatim: the frozen mapping
# from the 736x540 C to the renderer's square periodic buffer)
# ---------------------------------------------------------------------------
def c_render(Cf, n=None):
    """(540,736) frame-resolution C -> square (n,n) periodic buffer.

    The renderer (`v5_composite._sample_c`) indexes a SQUARE periodic buffer
    with one N for both axes, while the brief's C is at frame resolution
    (736x540).  The adapter resamples only along y (540 -> 736 rows) by linear
    interpolation, i.e. it hands the renderer a grid whose x samples are the
    real 736 and whose y samples are the interpolated 736; the renderer's own
    bilinear sampling then addresses the same C(u,v) surface.  Cost disclosure
    is in the report (sections 2 and 10).  x is not touched.
    """
    Cf = np.asarray(Cf, dtype=np.float64)
    H, W = Cf.shape
    n = int(n or W)
    yy = np.arange(n, dtype=np.float64) * (H - 1.0) / (n - 1.0)
    y0 = np.floor(yy).astype(np.int64)
    y1 = np.minimum(y0 + 1, H - 1)
    f = (yy - y0)[:, None]
    return ((1.0 - f) * Cf[y0, :] + f * Cf[y1, :]).astype(np.float32)


def render_optics(optics, n):
    """The v16 optics with the halo kernel re-sampled on an n x n C grid.

    The radial profile is stored in px of the state grid (512); the new C grid
    has n samples across the same square patch, so the radii are scaled by
    n/512 before `kernel_2d` is rebuilt.  Nothing else in the optics changes.
    """
    prof = optics["kernel_profile"]
    o = dict(optics)
    if prof is not None:
        scale = float(n) / float(BUF_GRID_N)
        p2 = dict(prof)
        p2["centres_px"] = np.asarray(prof["centres_px"], dtype=np.float64) * scale
        p2["edges_px"] = np.asarray(prof["edges_px"], dtype=np.float64) * scale
        p2["r_max_px"] = float(prof["r_max_px"]) * scale
        K, kinfo = V16.kernel_2d(p2, int(n))
        o["kernel_2d"] = K
        o["kernel_info_at_c_grid"] = kinfo
        o["kernel_radius_scale"] = scale
    return o


BUF_GRID_N = 512


# ---------------------------------------------------------------------------
# scenes
# ---------------------------------------------------------------------------
def build_scene(tag, n=None):
    """(state, env, optics, spec) for one stage.  Nothing is tuned here."""
    global BUF_GRID_N
    sc = SCN[tag]
    raw = R10.load_recipe()
    if sc["kind"] == "stream":
        wseed = int(derive_seed(raw, 1))
        env = RV15.scene_for(raw, R12.WIND_CELL["U"], R12.WIND_CELL["F"],
                             float(sc["depth"]), 256, wseed)["env"]
        env["depthM"] = float(sc["depth"])
        st = S3.StreamState3(N=int(n or sc["N"]), L=float(sc["L"]),
                             H=float(sc["depth"]), U=float(sc["U"]),
                             S0=float(sc["S0"]), D84=float(sc["D84"]),
                             seed=int(sc["seed"]), bed="natural",
                             alpha=SS.ALPHA_MAIN, dt_loop=2.0)
    elif sc["kind"] == "wind":
        r = RV15.scene_for(raw, float(sc["U"]), float(sc["F"]),
                           float(sc["depth"]), int(n or sc["N"]),
                           int(sc["seed"]))
        st, env = r["state"], r["env"]
        env["depthM"] = float(sc["depth"])
    else:
        raise SystemExit("unknown scene kind %r" % (sc["kind"],))
    spec = W.preset_spectra(sc["water"])
    optics = V16.optics_v16(spec, env, float(sc["depth"]), int(st.N), st.S)
    BUF_GRID_N = int(st.N)
    return st, env, optics, spec


def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def c_stats(C):
    C = np.asarray(C, dtype=np.float64)
    return {"mean": float(np.mean(C)), "std": float(np.std(C)),
            "min": float(np.min(C)), "max": float(np.max(C)),
            "p01": float(np.percentile(C, 1)),
            "p50": float(np.percentile(C, 50)),
            "p99": float(np.percentile(C, 99)),
            "shape": [int(C.shape[0]), int(C.shape[1])]}


def energy_record(info):
    """The per-frame energy bookkeeping, exactly as the deposit reports it."""
    return {"cell_energy_total": float(info["energy_cell_total"]),
            "deposited_total": float(info["energy_deposited"]),
            "deposited_after_guard": float(
                info.get("energy_deposited_after_guard", np.nan)),
            "rel_err": float(info["rel_err"]),
            "energy_guard_uniform_backstop": float(info["energy_guard"]),
            "offframe_fraction_area": float(info["offframe_fraction_area"]),
            "offframe_fraction_energy": float(info["offframe_fraction_energy"]),
            "degenerate_cells": int(info["n_degenerate_cell"]),
            "n_bulk_triangles": int(info["n_bulk"]),
            "n_tail_triangles": int(info["n_tail"]),
            "w_cell": float(info["w_cell"]),
            "normalisation_applied": 1.0,
            "normalisation_note": ("none: w_cell = W*H/N_s^2 makes flat water "
                                   "mean exactly 1 by definition, and no "
                                   "measured factor is ever applied")}


def render_frame(tag, i, st, env, optics, photo):
    global BUF_GRID_N
    t = i / float(FPS)
    fields = st.fields_at(t)
    C, info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, env)
    Cr = c_render(C)
    BUF_GRID_N = int(np.asarray(fields["gx"]).shape[0])   # the grid the optics
    # were built on (run_v4_1_8 fixed it at 512 for the same reason)
    o2 = render_optics(optics, Cr.shape[0])
    g1 = _g1(fields, env)
    img, dbg = V16.render_frame_v16(
        fields, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
        ambient_scale=float(A_FROZEN), c_layer_exposure=float(E_FROZEN),
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=float(B_SIGMA_FROZEN), return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"stage": tag, "frame": int(i), "t_s": float(t), "still": img,
            "C": C, "C_render": Cr, "g1": g1, "det_info": info, "metrics": met}


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def _func_src(text, name):
    """Source text of a top-level function, normalised to LF."""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    key = "\ndef %s(" % name
    i = text.find(key)
    if i < 0:
        return None
    j = text.find("\ndef ", i + len(key))
    k = text.find("\n# ---", i + len(key))
    ends = [x for x in (j, k) if x >= 0]
    end = min(ends) if ends else len(text)
    return text[i + 1:end].rstrip() + "\n"


def provenance_check():
    """The copied functions must be textually identical to their sources."""
    here = open(os.path.abspath(__file__), "r", encoding="utf-8").read()
    jobs = [("run_v4_1_5.py", MEAS_FUNCS), ("run_v4_1_7.py", ADAPTER_FUNCS)]
    rows = {}
    ok = True
    for fname, funcs in jobs:
        path = os.path.join(HERE, fname)
        if not os.path.isfile(path):
            rows[fname] = {"error": "missing", "path": path}
            ok = False
            continue
        src = open(path, "r", encoding="utf-8").read()
        per = {}
        for fn in funcs:
            a = _func_src(src, fn)
            b = _func_src(here, fn)
            same = (a is not None and b is not None and a == b)
            per[fn] = {"source_found": a is not None,
                       "copy_found": b is not None,
                       "identical": bool(same),
                       "chars_source": (len(a) if a else 0),
                       "chars_copy": (len(b) if b else 0)}
            ok = ok and bool(same)
        rows[fname] = per
    return {"pass": bool(ok),
            "method": ("each copied function's source text is extracted from "
                       "both files (LF-normalised, from 'def name(' to the "
                       "next top-level def or '# ---' banner) and compared "
                       "character for character"),
            "rows": rows}


def measurement_equivalence():
    """Re-measure the v4-1.7 still and compare with its recorded G-LINE rows."""
    out = {"reference": os.path.basename(V17_NUMBERS), "pass": None}
    if not (os.path.isfile(V17_STILL) and os.path.isfile(V17_NUMBERS)):
        out["skipped"] = "v4-1.7 artefacts missing"
        return out
    rec = read_json(V17_NUMBERS)
    rows_ref = (((rec.get("gates") or {}).get("G-LINE") or {}).get("rows") or [])
    if not rows_ref:
        out["skipped"] = "no recorded G-LINE rows"
        return out
    regions, _ = regions_from_photo()
    got = measure_image(V17_STILL, regions)
    per = []
    worst = 0.0
    for g_, r_ in zip(got, rows_ref):
        devs = {"depth_e3": abs(g_["depth_e3"] - r_["depth_v4_1_7"]),
                "depth_hp_e3": abs(g_["depth_hp_e3"] - r_["depth_hp_v4_1_7"]),
                "period_px": abs(g_["period_px"] - r_["period_v4_1_7_px"])}
        worst = max(worst, max(devs.values()))
        per.append({"region": g_["region"], "axis": g_["axis"],
                    "recomputed_depth_e3": g_["depth_e3"],
                    "recorded_depth_v4_1_7": r_["depth_v4_1_7"],
                    "abs_diff": devs})
    out.update({"rows": per, "max_abs_diff_over_all_rows": float(worst),
                "pass": bool(worst < 1e-9)})
    return out


def selftest_record():
    cmd = [sys.executable, os.path.join(HERE, "transport_mesh.py"),
           "--selftest"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    lines = (r.stdout or "").splitlines()
    return {"cmd": cmd, "returncode": int(r.returncode),
            "ok": bool(r.returncode == 0
                       and any("SELFTEST_OK" in x for x in lines)),
            "n_checks_ok": int(sum(1 for x in lines if x.startswith("[ok]"))),
            "n_checks_ng": int(sum(1 for x in lines if x.startswith("[NG]"))),
            "log": lines, "stderr_tail": (r.stderr or "")[-2000:]}


def deposit_path_equivalence(tag, st, env, n_tri=4000, box_default=16):
    """Bulk rasteriser vs tail rasteriser, on real projected triangles.

    Both paths compute the exact analytic triangle/pixel overlap; they differ
    only in summation order.  This measures that on this scene, so that the
    STAGE B performance switch cannot be confused with a change of method.
    """
    W_, H_ = WIDTH, HEIGHT
    fields = st.fields_at(STILL_IDX / float(FPS))
    PXe, PYe = TM.landing_pixels(fields, st.S, env, W=W_, H=H_)
    T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
    xmin = T[:, :, 0].min(1)
    xmax = T[:, :, 0].max(1)
    ymin = T[:, :, 1].min(1)
    ymax = T[:, :, 1].max(1)
    i0 = np.floor(xmin)
    j0 = np.floor(ymin)
    sw = np.floor(xmax) - i0 + 1
    sh = np.floor(ymax) - j0 + 1
    per = sw * sh
    dens = np.full(T.shape[0], float(W_ * H_) / float(A_cell.size * A_cell.size)
                   / (A_tri + 1e-30))
    rng = np.random.default_rng(20261115)
    tail = np.nonzero(per > box_default)[0]
    pool = tail if tail.size >= n_tri else np.nonzero(per > 1)[0]
    sub = np.sort(rng.choice(pool, size=min(n_tri, pool.size), replace=False))

    def newinfo():
        return {"area_in": 0.0, "area_out": 0.0, "energy_in": 0.0,
                "energy_out": 0.0, "bulk_boxes": 0, "tail_boxes": 0,
                "n_bulk": 0, "n_tail": 0, "energy_guard": 0.0,
                "n_degenerate_cell": 0}

    acc_t = np.zeros(W_ * H_)
    info_t = newinfo()
    t0 = time.time()
    for t in sub:
        TM._deposit_large(T[t], float(dens[t]), acc_t, W_, H_, info_t)
    wall_t = time.time() - t0
    acc_b = np.zeros(W_ * H_)
    info_b = newinfo()
    t0 = time.time()
    TM._deposit_bulk(T[sub], i0[sub], j0[sub], sw[sub], sh[sub], dens[sub],
                     acc_b, W_, H_, info_b)
    wall_b = time.time() - t0
    return {"stage": tag, "n_triangles": int(sub.size),
            "triangles_from": ("the >BOX_BULK_MAX population" if
                               tail.size >= n_tri else "all triangles"),
            "bbox_px_mean": float(per[sub].mean()),
            "max_abs_difference": float(np.max(np.abs(acc_b - acc_t))),
            "sum_bulk": float(acc_b.sum()), "sum_tail": float(acc_t.sum()),
            "wall_tail_s": float(wall_t), "wall_bulk_s": float(wall_b),
            "us_per_box_tail": float(1e6 * wall_t / max(per[sub].sum(), 1)),
            "us_per_box_bulk": float(1e6 * wall_b / max(per[sub].sum(), 1)),
            "pass": bool(np.max(np.abs(acc_b - acc_t)) < 1e-9)}


def loop_record(tag, st):
    """The loop rule of the brief, and what it implies at this patch."""
    T_loop = 2.0
    if SCN[tag]["kind"] == "stream":
        Us = float(st.Us)
        L = float(st.S)
        n_rule = max(1, int(round(T_loop * Us / L)))
        return {"rule": "n = round(T_loop*Us/L), L = T_loop*Us/n",
                "T_loop_s": float(T_loop), "n_loops_rule": int(n_rule),
                "L_pinned_m": L,
                "L_from_rule_at_that_n_m": float(T_loop * Us / n_rule),
                "T_loop_mode_s": float(L / Us),
                "n_cycles_in_piece": float(T_loop * Us / L),
                "note": ("L is PINNED by the brief to the screen window (no "
                         "magnification).  The (ii) advection closes when "
                         "L = n*Us*T; at this L the piece holds %.4f cycles, "
                         "so the 2.000 s seam is NOT closed (measured below)."
                         % (T_loop * Us / L))}
    # wind: the loop closes only if every active mode's phase advances by a
    # multiple of 2*pi in 2 s.  Recorded exactly as v19 did for its 8 s piece.
    import v8_physics as P8
    om0 = 2.0 * np.pi / T_loop
    table = {}
    for name in ("lo", "hi"):
        casc = getattr(st, name, None)
        if casc is None:
            continue
        try:
            wq = np.asarray(casc.omega_q, np.float64)
            h0 = np.asarray(casc.H0)
            active = h0 != 0
            k = np.asarray(casc.k, np.float64)[active]
            w = wq[active]
            om = P8.dispersion_omega(np.maximum(k, 1e-12), st.depth)
            n_mult = w / om0
            phase = w * T_loop
            table[name] = {
                "n_modes_active": int(active.sum()),
                "n_multiplier_unique": int(np.unique(np.round(n_mult)).size),
                "n_multiplier_min": float(np.min(n_mult)),
                "n_multiplier_max": float(np.max(n_mult)),
                "max_abs_dev_from_grid_rad_s": float(
                    np.max(np.abs(w - np.round(n_mult) * om0))),
                "max_phase_dev_at_T_rad": float(np.max(np.abs(
                    phase - 2.0 * np.pi * np.round(phase / (2.0 * np.pi))))),
                "continuum_omega_quantisation_error_rad_s": float(
                    np.max(np.abs(om - w)))}
        except Exception as e:  # pragma: no cover
            table[name] = {"error": "%s: %s" % (type(e).__name__, e)}
    return {"rule": "wind lane: no advective loop; closure needs "
                    "omega(k)*T = 2*pi*n for every active mode",
            "T_loop_s": float(T_loop), "omega0_rad_s": float(om0),
            "modes": table,
            "note": ("the 2 s seam is NOT closed by the wind lane either; "
                     "the residual is measured from the fields, not assumed")}


def seam_check(st, t_a=0.0, t_b=2.0):
    a = st.fields_at(t_a)
    b = st.fields_at(t_b)
    ha = np.asarray(a["h"], np.float64)
    hb = np.asarray(b["h"], np.float64)
    h_rms = float(np.sqrt(np.mean(ha ** 2)))
    gxa = np.asarray(a["gx"], np.float64)
    gxb = np.asarray(b["gx"], np.float64)
    dh = ha - hb
    return {"t_a_s": float(t_a), "t_b_s": float(t_b),
            "max_abs_dh_m": float(np.max(np.abs(dh))),
            "rms_dh_m": float(np.sqrt(np.mean(dh ** 2))),
            "h_rms_m": h_rms,
            "max_abs_dgx": float(np.max(np.abs(gxa - gxb))),
            "gx_rms": float(np.sqrt(np.mean(gxa ** 2
                                           + np.asarray(a["gy"],
                                                        np.float64) ** 2))),
            "seam_over_h_rms": float(np.sqrt(np.mean(dh ** 2))
                                     / max(h_rms, 1e-30))}


def scene_record(tag, st, env, optics, spec):
    sc = SCN[tag]
    raw = R10.load_recipe()
    Wm, Hm = sc["window_m"]
    rec = {
        "stage": tag, "label": sc["label"], "kind": sc["kind"],
        "patch_m": float(st.S), "grid": int(st.N),
        "cell_m": float(st.S / st.N),
        "window_m": [float(Wm), float(Hm)],
        "window_rule": ("the generated patch spans the frame WIDTH: "
                        "window = [L, L*540/736]; the 736:540 pixel aspect is "
                        "v2's framing and is not a magnification"),
        "magnification": 1.0,
        "photo": PHOTO,
        "photo_sha256": sha256_file(PHOTO),
        "depth_m": float(sc["depth"]),
        "water_preset": spec.get("preset"), "water_label": spec.get("label"),
        "a_rgb": [float(x) for x in spec["a"]],
        "b_rgb": [float(x) for x in spec["b"]],
        "phase": {"kind": optics.get("phase_kind"),
                  "label": (optics.get("phase") or {}).get("label"),
                  "n": (optics.get("phase") or {}).get("n"),
                  "mu": (optics.get("phase") or {}).get("mu"),
                  "provenance": optics.get("phase_provenance")},
        "optics": {"kernel": optics.get("kernel"),
                   "w_ballistic": [float(x) for x in optics["w_ballistic"]],
                   "w_halo": [float(x) for x in optics["w_halo"]],
                   "sigma_px": float(optics.get("sigma_px") or 0.0),
                   "sigma_px_v12_gaussian": float(
                       optics.get("sigma_px_v12_hg") or 0.0),
                   "kernel_r_max_px_state_grid": (
                       float(optics["kernel_profile"]["r_max_px"])
                       if optics["kernel_profile"] else None),
                   "kernel_r_max_m": (
                       float(optics["kernel_profile"]["r_max_px"]) * st.S
                       / float(st.N) if optics["kernel_profile"] else None)},
        "composite": {"A": float(A_FROZEN), "Blow_sigma_px":
                      float(B_SIGMA_FROZEN), "E": float(E_FROZEN),
                      "Q": int(Q_DET), "z_refraction_m": float(sc["depth"]),
                      "output_px": [int(WIDTH), int(HEIGHT)]},
        "env": {"sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
                "sunDir": [float(x) for x in env["sunDir"]],
                "clarity": float(env["clarity"]),
                "turbid": float(env["turbid"]),
                "depthM": float(env["depthM"]),
                "eta": 1.0,
                "eta_source": ("the recipe carries no sky block, so "
                               "eta = DNI/(DNI+DFF) = 1.0 (direct sun only)")},
        "loop": loop_record(tag, st),
        "seam": seam_check(st),
        "deposit": {"module": "transport_mesh (VERBATIM)",
                    "method": "area push-forward (quad -> 2 triangles, exact "
                              "overlap with the frame pixels)",
                    "c_grid": [int(HEIGHT), int(WIDTH)],
                    "w_cell": float(WIDTH * HEIGHT / float(st.N * st.N)),
                    "per_cell_energy_rule": "w_cell = W*H/N_s^2",
                    "kernel": "none", "jitter": "none", "band_limit": "none",
                    "mean_normalisation": "none",
                    "frame_energy_total_expected": float(WIDTH * HEIGHT),
                    "box_bulk_max_used": int(TM.BOX_BULK_MAX),
                    "nyquist_cells_per_m": float(st.N / st.S),
                    "main_run_nyquist_cells_per_m": float(512 / 1.4755102),
                    "nyquist_ratio_vs_main": float((st.N / st.S)
                                                   / (512 / 1.4755102))},
        "buffer_adapter": {"c_frame_shape": [int(HEIGHT), int(WIDTH)],
                           "renderer_grid": "square 736x736 by linear y "
                                            "interpolation (c_render, "
                                            "run_v4_1_7 verbatim)",
                           "kernel_profile_reference_grid": int(BUF_GRID_N)},
        "created_utc": now_iso(),
    }
    if sc["kind"] == "wind":
        rec["wind"] = {"U_ms": float(sc["U"]), "fetch_m": float(sc["F"]),
                       "spectrum_source": V15S.SPECTRUM_SOURCE_DEFAULT,
                       "arm": V15S.MERGE_ARM,
                       "omega_c_rad_s": float(st.ctx["omega_c"]),
                       "warn_omega_c_threshold": float(V15S.WARN_OMEGA_C),
                       "warn_extrapolated_wave_age": bool(
                           float(st.ctx["omega_c"]) > V15S.WARN_OMEGA_C),
                       "lambda_p_m": float(st.laws["lambda_p_m"]),
                       "patch_doublings": int(st.patch_doublings),
                       "spectrum_after": {
                           k: float(v) for k, v in st.spectrum_after.items()
                           if isinstance(v, (int, float))},
                       "seed_source": ("run_v19_prod.SEED_BRIEF = "
                                       "3255045276382562994, which equals "
                                       "state.derive_seed(recipe, 1)"),
                       "seed_equals_derive_seed": bool(
                           int(sc["seed"]) == int(derive_seed(raw, 1)))}
    else:
        acc = st.accounting()
        rec["stream"] = {
            "H_m": float(st.H), "U_ms": float(st.U), "U0_ms": float(st.Us),
            "S0": float(st.S0), "D84_m": float(st.D84),
            "seed": int(sc["seed"]), "alpha": float(SS.ALPHA_MAIN),
            "accounting": acc,
            "sigma_total_m": float(st.sigma_total),
            "var_i_m2": float(st.var_i_physical),
            "var_ii_m2": float(st.var_ii), "var_iii_m2": float(st.var_iii),
            "validity": SS.validity_label(st.derived["Fr"],
                                          float(sc["depth"]),
                                          float(sc["D84"])),
            "derived": {k: (float(v) if isinstance(v, (int, float)) else v)
                        for k, v in st.derived.items()}}
    return rec


def grid_cost_probe(tag, grids, us_box_bulk, us_box_tail):
    """Measured projection cost per grid, at frame 24.

    For every grid: the fold ratio (sum of projected triangle areas / frame
    area), the summed integer-bbox area (the unit the rasteriser is charged
    in) and the resulting wall-time estimates.  No deposit is run here -- this
    is the projection only, so the numbers are the cost model that the report
    uses to justify the grid actually executed.
    """
    out = {"us_per_box_bulk": float(us_box_bulk),
           "us_per_box_tail": float(us_box_tail), "by_grid": {}}
    for n in grids:
        t0 = time.time()
        st, env, _o, _s = build_scene(tag, n=int(n))
        f = st.fields_at(STILL_IDX / float(FPS))
        PXe, PYe = TM.landing_pixels(f, st.S, env, W=WIDTH, H=HEIGHT)
        T_, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
        xmin = T_[:, :, 0].min(1)
        xmax = T_[:, :, 0].max(1)
        ymin = T_[:, :, 1].min(1)
        ymax = T_[:, :, 1].max(1)
        i0 = np.floor(xmin)
        j0 = np.floor(ymin)
        per = (np.floor(xmax) - i0 + 1) * (np.floor(ymax) - j0 + 1)
        tail_frac = float((per > 16).mean())
        est_b = 1e-6 * float(per.sum()) * float(us_box_bulk)
        est_t = 1e-6 * float(per.sum()) * float(us_box_tail)
        out["by_grid"][str(int(n))] = {
            "grid": int(n),
            "cell_mm": float(1e3 * st.S / st.N),
            "nyquist_cells_per_m": float(st.N / st.S),
            "sum_abs_triangle_area_px": float(A_tri.sum()),
            "fold_ratio_vs_frame": float(A_tri.sum() / (WIDTH * HEIGHT)),
            "sum_bbox_px": float(per.sum()),
            "bbox_median_px": float(np.median(per)),
            "tail_fraction_at_default_box_bulk_max": tail_frac,
            "est_wall_s_one_frame_bulk_path": float(est_b),
            "est_wall_s_one_frame_tail_path": float(est_t),
            "est_wall_s_48_frames_bulk_path": float(48.0 * est_b),
            "probe_wall_s": float(time.time() - t0),
        }
        log("        probe %4d^2: fold %7.1f x frame ; sum bbox %.3e px ; "
            "est %.0f s/frame (bulk) / %.0f s/frame (tail)"
            % (int(n), out["by_grid"][str(int(n))]["fold_ratio_vs_frame"],
               float(per.sum()), est_b, est_t))
        del st, _o, _s, T_, A_tri, A_cell, PXe, PYe
        gc.collect()
    return out


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record + selftest")
    makedirs()
    cp = Checkpoint(CP_PATH)

    # ---- the brief's own list, explicitly -------------------------------
    brief_need = [
        os.path.join(HERE, "transport_mesh.py"),
        os.path.join(HERE, "run_v4_1_8.py"),
        os.path.join(HERE, "stream_surface3.py"),
        PHOTO,                                    # both stages' base photo
        os.path.join(HERE, "run_v20_prod.py"),    # the v20 photo's provenance
        os.path.join(HERE, "run_v19_prod.py"),    # the v19 photo's provenance
        os.path.join(HERE, "run_v15.py"),         # the wind lane generator
        os.path.join(HERE, "run_v15_state.py"),
        os.path.join(HERE, "run_v15_spectrum.py"),
    ]
    # everything this driver actually imports or measures
    need = brief_need + [
        BRIEF,
        os.path.join(HERE, "run_v4_1_15.py"),
        os.path.join(HERE, "run_v10.py"),
        os.path.join(HERE, "run_v4_1_2.py"),
        os.path.join(HERE, "stream_surface2.py"),
        os.path.join(HERE, "run_v16_water.py"),
        os.path.join(HERE, "v12_water.py"),
        os.path.join(HERE, "det.py"),
        os.path.join(HERE, "transport.py"),
        os.path.join(HERE, "state.py"),
        os.path.join(HERE, "image.py"),
        os.path.join(HERE, "v5_composite.py"),
        os.path.join(HERE, "v5_physics.py"),
        os.path.join(HERE, "v8_physics.py"),
        os.path.join(HERE, "run_v4_1_5.py"),
        os.path.join(HERE, "run_v4_1_7.py"),
        V17_STILL, V17_NUMBERS, V18_STILL, V41_8_NUMBERS,
    ]
    missing = [p for p in need if not os.path.isfile(p)]
    log("  brief list:")
    for p in brief_need:
        log("    %-24s %s" % (os.path.basename(p),
                              "OK" if os.path.isfile(p) else "MISSING"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))

    # ---- what the brief's "same file" claims actually point at -----------
    photo_sha = sha256_file(PHOTO)

    def declared_photo(fname):
        """The literal a driver assigns to PHOTO (textual check, not trust)."""
        p = os.path.join(HERE, fname)
        if not os.path.isfile(p):
            return None
        txt = open(p, "r", encoding="utf-8", errors="replace").read()
        m = re.search(r'^PHOTO\s*=\s*r?"([^"]+)"', txt, flags=re.M)
        return m.group(1) if m else None

    prov_photo = {}
    for name, path, drivers in (
            ("v4_1_8", V41_8_NUMBERS, ["run_v4_1_8.py", "run_v4_1_2.py"]),
            ("v19", os.path.join(HERE, "out", "v19",
                                 "numbers-v19-prod.json"), ["run_v19_prod.py"]),
            ("v20", os.path.join(HERE, "out", "v20",
                                 "numbers-v20-prod.json"), ["run_v20_prod.py"])):
        node = {"drivers": {}, "recorded": None}
        for dn in drivers:
            q = declared_photo(dn)
            node["drivers"][dn] = {
                "declared_photo": q,
                "equals_this_job": bool(q == PHOTO) if q else None}
        if os.path.isfile(path):
            try:
                d = read_json(path)
                ph = (d.get("preflight") or {}).get("photo") or {}
                if ph:
                    node["recorded"] = {
                        "path": ph.get("path"), "sha256": ph.get("sha256"),
                        "matches_this_job": bool(ph.get("sha256") == photo_sha)}
                else:
                    node["recorded"] = {
                        "note": ("this record does not store the photo; the "
                                 "identity rests on the driver's PHOTO "
                                 "constant (textual check above)")}
            except Exception as e:
                node["recorded"] = {"error": "%s: %s" % (type(e).__name__, e)}
        else:
            node["recorded"] = {"error": "missing", "path": path}
        prov_photo[name] = node

    # the wind lane: the generator and the constants that make it Elfouhaily
    import run_v15_state as S15
    wind_lane = {
        "generator": "run_v15.scene_for(raw, U, F, h, N, seed) -> "
                     "run_v15_state.V15WaterState(..., spectrum_source, arm)",
        "spectrum_source_default": V15S.SPECTRUM_SOURCE_DEFAULT,
        "merge_arm": V15S.MERGE_ARM,
        "warn_omega_c_threshold": float(V15S.WARN_OMEGA_C),
        "callable_with_U6_F100": True,
        "class": "run_v15_state.V15WaterState",
        "class_found": hasattr(S15, "V15WaterState"),
        "used_by": ["run_v4_1_11.py", "run_v4_1_12.py", "run_v19_prod.py"],
        "note": ("the Elfouhaily-form / Arm D generator of the wind lane; the "
                 "preset-5 condition U = 6 m/s, F = 100 m is passed straight "
                 "through scene_for, nothing is substituted"),
    }
    rec = {"brief": BRIEF, "checked": need, "brief_list": brief_need,
           "missing": [], "checked_utc": now_iso(),
           "photo": {"path": PHOTO, "sha256": photo_sha,
                     "bytes": os.path.getsize(PHOTO),
                     "provenance": prov_photo,
                     "claim": ("brief: the contrast frame uses 'the same file "
                               "as v4-1.8 / v20', the sub frame 'the same "
                               "file as v19'.  Both claims are checked against "
                               "each production run's own preflight record.")},
           "wind_lane": wind_lane}

    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in NEW_FILES}
    payload = {"all": hashes, "new_files_of_this_job": NEW_FILES,
               "existing_count": len(existing), "existing": existing,
               "recorded_utc": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's new file (its stage-0 hash is in 'all' too).  "
                        "transport_mesh.py is an EXISTING file here (written "
                        "by v4-1.8) and this job must not change it.")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(NEW_FILES)))
    log("  photo %s sha %s" % (os.path.basename(PHOTO), photo_sha[:16]))
    for k, v in prov_photo.items():
        log("    %-8s %s" % (k, v))

    # ---- per-stage scene records ----------------------------------------
    scenes, deposits, eq = {}, {}, {}
    for tag in STAGES:
        # STAGE B runs with the documented performance switch; A does not.
        if tag == "sub":
            TM.BOX_BULK_MAX = 1 << 20
        else:
            TM.BOX_BULK_MAX = 16
        st, env, optics, spec = build_scene(tag)
        scenes[tag] = scene_record(tag, st, env, optics, spec)
        log("  [%s] L=%.4f m grid=%d (cell %.3f mm, nyq/main %.2f) "
            "depth=%.2f m"
            % (tag, st.S, st.N, 1e3 * st.S / st.N,
               scenes[tag]["deposit"]["nyquist_ratio_vs_main"],
               float(scenes[tag]["depth_m"])))
        log("        loop %s ; seam rms|dh| %.4f mm (%.3f x h_rms)"
            % (("n=%s, T_mode=%s s"
                % (scenes[tag]["loop"].get("n_loops_rule"),
                   ("%.4f" % scenes[tag]["loop"]["T_loop_mode_s"])
                   if "T_loop_mode_s" in scenes[tag]["loop"] else "n/a"))
               if "n_loops_rule" in scenes[tag]["loop"]
               else "wind lane (no advective loop)",
               1e3 * scenes[tag]["seam"]["rms_dh_m"],
               scenes[tag]["seam"]["seam_over_h_rms"]))
        eq[tag] = deposit_path_equivalence(tag, st, env)
        log("        deposit paths: %d tri, max|bulk-tail| = %.2e "
            "(bulk %.2f us/box, tail %.2f us/box)"
            % (eq[tag]["n_triangles"], eq[tag]["max_abs_difference"],
               eq[tag]["us_per_box_bulk"], eq[tag]["us_per_box_tail"]))
        del st, optics, spec
        gc.collect()
        # cost of the projection at every grid that matters for this window
        side = SCN[tag]
        grids = [int(x) for x in side["probe_grids"]]
        deposits[tag] = grid_cost_probe(tag, grids, eq[tag]["us_per_box_bulk"],
                                        eq[tag]["us_per_box_tail"])
        deposits[tag]["base_grid"] = int(side["N"])
        deposits[tag]["half_grid"] = int(side["half"])
        deposits[tag]["refine_grid"] = side.get("refine")
        deposits[tag]["brief_guide_grid"] = 1024
        deposits[tag]["band_rule_grid_min"] = int(np.ceil(
            side["L"] * (512 / 1.4755102)))
    TM.BOX_BULK_MAX = 16

    stt = selftest_record()
    log("  transport_mesh selftest: %s (%d ok / %d NG)"
        % ("OK" if stt["ok"] else "NG", stt["n_checks_ok"],
           stt["n_checks_ng"]))
    prov = provenance_check()
    log("  measurement provenance: %s"
        % ("identical" if prov["pass"] else "DIFFERS"))
    meq = measurement_equivalence()
    log("  v4-1.7 G-LINE reproduction: %s (max diff %s)"
        % (meq.get("pass"), (("%.2e" % meq["max_abs_diff_over_all_rows"])
                             if meq.get("max_abs_diff_over_all_rows")
                             is not None else "n/a")))

    checkpoint_extra("stage0", rec)
    checkpoint_extra("scenes", scenes)
    checkpoint_extra("deposit_cost", deposits)
    checkpoint_extra("deposit_path_equivalence", eq)
    checkpoint_extra("selftest", stt)
    checkpoint_extra("provenance", prov)
    checkpoint_extra("measurement_equivalence", meq)
    cp.unit_done("stage-0-preflight")
    log("  wall %.1fs" % (time.time() - t0))
    return scenes


# ---------------------------------------------------------------------------
# stage: frame 24
# ---------------------------------------------------------------------------
def do_frame24(tag):
    t0 = time.time()
    log("=" * 78)
    log("[frame24 %s] area deposit, convergence, point-splat reference, G-LINE"
        % tag)
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc = SCN[tag]
    if tag == "sub":
        TM.BOX_BULK_MAX = 1 << 20
    else:
        TM.BOX_BULK_MAX = 16
    regions, rinfo = regions_from_photo()
    st, env, optics, spec = build_scene(tag)
    photo = load_photo()
    t = STILL_IDX / float(FPS)

    t1 = time.time()
    r = render_frame(tag, STILL_IDX, st, env, optics, photo)
    deposit_wall = time.time() - t1
    still_path = os.path.join(OUT, sc["still"])
    sha = save_png(still_path, r["still"])
    sha_f = save_png(frame_path(tag, STILL_IDX), r["still"])
    info = r["det_info"]
    log("  frame 24: C mean %.9f std %.4f p99 %.4f max %.2f  rel_err %.2e "
        "fold %.1f"
        % (float(np.mean(r["C"])), float(np.std(r["C"])),
           float(np.percentile(r["C"], 99)), float(np.max(r["C"])),
           info["rel_err"], float(info["area_in"] + info["area_out"])
           / float(WIDTH * HEIGHT)))

    rows = lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
    photo_rows = measure_image(PHOTO, regions)
    per_row = []
    for got, base in zip(rows, photo_rows):
        per_row.append({"region": got["region"], "axis": got["axis"],
                        "period_px": got["period_px"],
                        "depth_this": got["depth_e3"],
                        "depth_photo": base["depth_e3"],
                        "ratio": got["depth_e3"] / base["depth_e3"],
                        "depth_hp_this": got["depth_hp_e3"],
                        "depth_hp_photo": base["depth_hp_e3"],
                        "ratio_hp": got["depth_hp_e3"] / base["depth_hp_e3"],
                        "ok": bool(got["depth_e3"] <= base["depth_e3"]),
                        "ok_hp": bool(got["depth_hp_e3"]
                                      <= base["depth_hp_e3"])})

    # ---- same field, point splat (comparison material for G-PEAK) --------
    fields = st.fields_at(t)
    ps, ps_info = TM.point_splat_C(fields, st.S, WIDTH, HEIGHT, env,
                                   Q=int(Q_DET))
    ps_sharp, _ = TM.point_splat_C(fields, st.S, WIDTH, HEIGHT, env,
                                   Q=int(Q_DET), blur=False)
    cd, cd_info = det_frame(fields, st.S, int(st.N), env, int(Q_DET))
    cs, css, cdet = c_stats(r["C"]), c_stats(ps), c_stats(cd)
    cs_sharp = c_stats(ps_sharp)
    del ps, ps_sharp, cd
    gc.collect()

    # ---- convergence: half grid on the same frame ------------------------
    st_h, env_h, _oh, _sh = build_scene(tag, n=int(sc["half"]))
    t_co = time.time()
    Ch, info_h = TM.mesh_C(st_h.fields_at(t), st_h.S, WIDTH, HEIGHT, env_h)
    conv_wall = time.time() - t_co
    ch = c_stats(Ch)
    del Ch, st_h, _oh, _sh
    gc.collect()
    conv = {"n_base": int(st.N), "n_half": int(sc["half"]),
            "C_base": cs, "C_half": ch,
            "std_rel_diff": float(ch["std"] / cs["std"] - 1.0),
            "p99_rel_diff": float(ch["p99"] / cs["p99"] - 1.0),
            "within_15pct": bool(abs(ch["std"] / cs["std"] - 1.0) <= 0.15
                                 and abs(ch["p99"] / cs["p99"] - 1.0) <= 0.15),
            "energy_half": energy_record(info_h),
            "wall_s": float(conv_wall),
            "note": ("the water grid changes, the frame grid does not; the "
                     "difference is the sub-pixel footprint of one projected "
                     "cell")}

    # ---- STAGE B: also one step UP (the brief's guide is 1024^2) ---------
    refine = None
    if sc.get("refine"):
        st_r, env_r, _or, _sr = build_scene(tag, n=int(sc["refine"]))
        t_rf = time.time()
        Cr, info_r = TM.mesh_C(st_r.fields_at(t), st_r.S, WIDTH, HEIGHT, env_r)
        ref_wall = time.time() - t_rf
        cr = c_stats(Cr)
        del Cr, st_r, _or, _sr
        gc.collect()
        refine = {"n_refine": int(sc["refine"]), "C_refine": cr,
                  "std_rel_diff": float(cr["std"] / cs["std"] - 1.0),
                  "p99_rel_diff": float(cr["p99"] / cs["p99"] - 1.0),
                  "within_15pct": bool(
                      abs(cr["std"] / cs["std"] - 1.0) <= 0.15
                      and abs(cr["p99"] / cs["p99"] - 1.0) <= 0.15),
                  "energy_refine": energy_record(info_r),
                  "wall_s": float(ref_wall),
                  "note": ("supplementary: the brief's guide grid for this "
                           "window is 1024^2; this is one frame at that grid, "
                           "so the convergence question can be answered "
                           "upward as well as downward")}
        log("  refine %d^2: std %+.2f%% p99 %+.2f%% (%.0f s)"
            % (sc["refine"], 100 * refine["std_rel_diff"],
               100 * refine["p99_rel_diff"], ref_wall))

    peak = {"C_base": cs, "C_point_splat": css, "C_det_frame_state_grid": cdet,
            "C_point_splat_sharp": cs_sharp,
            "ref_point_splat": ps_info, "ref_det_frame": cd_info,
            "std_ratio_vs_point_splat": float(cs["std"] / css["std"]),
            "p99_ratio_vs_point_splat": float(cs["p99"] / css["p99"]),
            "std_ratio_vs_point_splat_sharp": float(cs["std"] / cs_sharp["std"]),
            "p99_ratio_vs_point_splat_sharp": float(cs["p99"] / cs_sharp["p99"]),
            "std_ratio_vs_det_frame": float(cs["std"] / cdet["std"]),
            "p99_ratio_vs_det_frame": float(cs["p99"] / cdet["p99"]),
            "criterion": ("std and p99 of the frame-24 C must be >= 90% of the "
                          "point-splat reference computed on the SAME field "
                          "and frame grid"),
            "note": ("the reference lanes are comparison material only; they "
                     "are never used for C of the deliverable")}

    out = {"stage": tag, "frame": STILL_IDX, "t_s": float(t),
           "still": os.path.basename(still_path), "still_sha256": sha,
           "frame_0024_sha256": sha_f,
           "still_equals_frame_0024": bool(sha == sha_f),
           "regions": rinfo, "measurements": rows,
           "per_row_vs_photo": per_row,
           "strength": lattice_strength(rows),
           "photo_strength": lattice_strength(photo_rows),
           "pass_gl_line": bool(all(p["ok"] for p in per_row)),
           "C_frame_stats": cs, "C_render_stats": c_stats(r["C_render"]),
           "C_render_mean_vs_frame_mean": float(
               np.mean(r["C_render"]) / float(np.mean(r["C"]))),
           "energy": energy_record(info),
           "deposit_detail": {k: (float(v) if isinstance(v, (int, float))
                                  and not isinstance(v, bool) else v)
                              for k, v in info.items()
                              if not isinstance(v, dict)},
           "peak": peak, "convergence": conv, "refinement": refine,
           "image_metrics": r["metrics"], "g1": float(r["g1"]),
           "wall_s": float(time.time() - t0),
           "deposit_wall_s": float(deposit_wall)}

    for p in per_row:
        log("  %-6s %s  %.2f px  depth %.1f (photo %.1f) ratio %.3f  hp %.3f  %s"
            % (p["region"], p["axis"], p["period_px"], p["depth_this"],
               p["depth_photo"], p["ratio"], p["ratio_hp"],
               "OK" if p["ok"] else "FAIL"))
    log("  G-LINE pre-check: %s (strength %.1f vs photo %.1f)"
        % ("PASS" if out["pass_gl_line"] else "FAIL", out["strength"],
           out["photo_strength"]))
    log("  G-PEAK pre-check: std %.4f vs point splat %.4f (%.3f) / det %.4f "
        "(%.3f); p99 %.3f vs %.3f (%.3f)"
        % (cs["std"], css["std"], peak["std_ratio_vs_point_splat"],
           cdet["std"], peak["std_ratio_vs_det_frame"], cs["p99"], css["p99"],
           peak["p99_ratio_vs_point_splat"]))
    log("  G-CONV pre-check: dstd %+.2f%% dp99 %+.2f%% -> %s (%.1fs)"
        % (100 * conv["std_rel_diff"], 100 * conv["p99_rel_diff"],
           "PASS" if conv["within_15pct"] else "FAIL", conv_wall))

    f24 = (read_checkpoint().get("frame24") or {})
    f24[tag] = out
    checkpoint_extra("frame24", f24)
    cp.unit_done("frame-0024-%s" % tag)
    del st, optics, photo, r, spec, fields
    gc.collect()
    TM.BOX_BULK_MAX = 16
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: frames
# ---------------------------------------------------------------------------
def chunk_done(tag, k, state):
    node = ((state.get("chunks") or {}).get(tag) or {}).get(str(k))
    if not node or len(node.get("frames", {})) != CHUNK:
        return False
    for name, sha in node["frames"].items():
        p = frame_path(tag, int(name))
        if not os.path.isfile(p) or sha256_file(p) != sha:
            return False
    return True


def do_chunk(tag, k):
    t0 = time.time()
    log("=" * 78)
    log("[chunk %s %02d] frames %d-%d" % (tag, k, k * CHUNK,
                                          (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    if chunk_done(tag, k, state):
        log("  chunk %s %02d already committed -- skip" % (tag, k))
        return True
    if tag == "sub":
        TM.BOX_BULK_MAX = 1 << 20
    else:
        TM.BOX_BULK_MAX = 16
    st, env, optics, spec = build_scene(tag)
    photo = load_photo()
    still_sha = ((state.get("frame24") or {}).get(tag) or {}).get(
        "still_sha256")
    frames, metrics, energies = {}, [], []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_frame(tag, i, st, env, optics, photo)
        sha = save_png(frame_path(tag, i), r["still"])
        m, c = r["metrics"], r["C"]
        er = energy_record(r["det_info"])
        rec = {"frame": i, "t_s": r["t_s"], "sha256": sha, "g1": r["g1"],
               "C_mean": float(np.mean(c, dtype=np.float64)),
               "C_std": float(np.std(c, dtype=np.float64)),
               "C_min": float(np.min(c)), "C_max": float(np.max(c)),
               "C_p50": float(np.percentile(c, 50)),
               "C_p99": float(np.percentile(c, 99)),
               "energy_rel_err": er["rel_err"],
               "offframe_fraction_energy": er["offframe_fraction_energy"],
               "energy_deposited": er["deposited_total"],
               "energy_cell_total": er["cell_energy_total"],
               "crush": float(m["crush_rate_linear_lt_0.02"]),
               "blow": float(m["blowout_rate_aces_gt_0.98"]),
               "lin_mean": float(m["linear_mean_pre_tonemap"])}
        if i == STILL_IDX and still_sha:
            rec["equals_frame24_still"] = bool(sha == still_sha)
        frames["%04d" % i] = sha
        metrics.append(rec)
        energies.append(er)
        log("  frame %04d t=%.3f Cmean=%.6f Cstd=%.4f p99=%.3f E_rel=%.1e "
            "off=%.3f %.1fs"
            % (i, r["t_s"], rec["C_mean"], rec["C_std"], rec["C_p99"],
               rec["energy_rel_err"], rec["offframe_fraction_energy"],
               time.time() - t_start))
        del r, m, er
        gc.collect()
    node = {"frames": frames, "metrics": metrics, "energies": energies,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start),
            "box_bulk_max_used": int(TM.BOX_BULK_MAX),
            "energy_max_rel_err": float(max(e["rel_err"] for e in energies)),
            "energy_mean_abs_rel_err": float(np.mean(
                [abs(e["rel_err"]) for e in energies]))}
    chunks = read_checkpoint().get("chunks") or {}
    node_all = chunks.get(tag) or {}
    node_all[str(k)] = node
    chunks[tag] = node_all
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%s-%02d" % (tag, k))
    log("  chunk %s %02d wall %.1fs  energy max rel_err %.2e"
        % (tag, k, node["wall_s"], node["energy_max_rel_err"]))
    del st, optics, photo, spec
    gc.collect()
    TM.BOX_BULK_MAX = 16
    return True


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP, record only)
# ---------------------------------------------------------------------------
def do_warp(tag, i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp %s] displacement yardstick at frame %d (record only)"
        % (tag, i))
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    st, env, optics, spec = build_scene(tag)
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, float(B_SIGMA_FROZEN))
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    depth = float(env["depthM"])
    per = {}
    for bb, tg in enumerate(("r", "g", "b")):
        cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx, shy = -nx * depth * cb, -ny * depth * cb
        dx_px = shx / S * (WIDTH - 1.0)
        dy_px = shy / S * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        per[tg] = {"px_mean": float(np.mean(mag)),
                   "px_p50": float(np.percentile(mag, 50)),
                   "px_p99": float(np.percentile(mag, 99)),
                   "px_max": float(np.max(mag)),
                   "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                   "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    out = {"stage": tag, "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres, "
                     "UNFILTERED slope (same method as v4-1.4/5/7/8)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "criterion": "record only (the brief sets no threshold)",
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp-%s-frame-0024.json" % tag), out)
    warp = read_checkpoint().get("warp") or {}
    warp[tag] = out
    checkpoint_extra("warp", warp)
    cp.unit_done("warp-%s-frame-%04d" % (tag, i))
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per["g"]["px_max"]))
    del st, optics, spec
    gc.collect()
    return out


# ---------------------------------------------------------------------------
# stage: encode + probe
# ---------------------------------------------------------------------------
FFMPEG_VENDOR = os.path.join(HERE, "_vendor", "imageio_ffmpeg", "binaries",
                             "ffmpeg-win-x86_64-v9.0.2.exe")


def ffmpeg_exe():
    for cand in (os.environ.get("FFMPEG_BIN"), FFMPEG_VENDOR):
        if cand and os.path.isfile(cand):
            return cand
    return "ffmpeg"


def decode_green(path):
    cmd = [ffmpeg_exe(), "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def ffprobe(path):
    """Stream facts of `path` (ffprobe's field names, comparable directly)."""
    exe = ffmpeg_exe()
    out = {"probe_binary": os.path.basename(exe)}
    try:
        r = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                            "-count_frames", "-show_entries",
                            "stream=codec_name,width,height,r_frame_rate,"
                            "nb_read_frames,pix_fmt,duration", "-of", "json",
                            path], capture_output=True, text=True)
        if r.returncode == 0:
            d = json.loads(r.stdout)["streams"][0]
            d["nb_frames_reported"] = int(d.get("nb_read_frames") or 0)
            d["duration_s"] = float(d.get("duration") or 0.0)
            d["probe_method"] = "ffprobe -count_frames"
            d["probe_binary"] = "ffprobe"
            return d
        out["ffprobe_error"] = (r.stderr or "").strip()[:400]
    except Exception as e:
        out["ffprobe_error"] = "%s: %s" % (type(e).__name__, e)

    r = subprocess.run([exe, "-hide_banner", "-i", path],
                       capture_output=True, text=True)
    err = (r.stderr or "") + (r.stdout or "")
    m = re.search(r"Duration:\s*(\d+):(\d+):(\d+\.\d+)", err)
    if m:
        out["duration_s"] = (int(m.group(1)) * 3600 + int(m.group(2)) * 60
                             + float(m.group(3)))
        out["duration_raw"] = m.group(0).split(":", 1)[1].strip()
    m = re.search(r"Video:\s*([A-Za-z0-9_\-]+)", err)
    if m:
        out["codec_name"] = m.group(1)
    m = re.search(r"Video:.*?(\d{2,5})x(\d{2,5})", err)
    if m:
        out["width"] = int(m.group(1))
        out["height"] = int(m.group(2))
    m = re.search(r"([\d.]+)\s*fps", err)
    if m:
        fps = float(m.group(1))
        out["r_frame_rate"] = ("%d/1" % int(round(fps))
                              if abs(fps - round(fps)) < 1e-9 else str(fps))
        out["fps_parsed"] = fps
    for pf in ("yuv420p", "yuvj420p", "yuv444p", "gbrp", "rgb24"):
        if pf in err:
            out["pix_fmt"] = pf
            break
    g = decode_green(path)
    out["nb_frames_reported"] = int(g.shape[0])
    out["frames_by_full_decode"] = int(g.shape[0])
    out["probe_method"] = ("ffmpeg -i header parse (duration/codec/size/fps/"
                          "pix_fmt) + full rgb24 decode (frame count)")
    return out


def do_encode(tag):
    t0 = time.time()
    log("=" * 78)
    log("[encode %s] MP4 + probe" % tag)
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc = SCN[tag]
    exe = ffmpeg_exe()
    mp4 = os.path.join(OUT, sc["mp4"])
    pattern = os.path.join(OUT, sc["frames"], "frame-%04d.png")
    cmd = [exe, "-y", "-v", "error", "-framerate", str(FPS),
           "-i", pattern, "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-crf", "16", "-preset", "medium", mp4]
    log("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("ffmpeg failed: %s" % r.stderr[-2000:])
    probe = ffprobe(mp4)
    out = {"stage": tag, "mp4": os.path.basename(mp4), "cmd": cmd,
           "probe": probe, "ffmpeg_binary": exe,
           "ffmpeg_version": (subprocess.run(
               [exe, "-version"], capture_output=True, text=True
           ).stdout or "").splitlines()[:1],
           "size_bytes": os.path.getsize(mp4),
           "sha256": sha256_file(mp4), "created_utc": now_iso(),
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "encode-%s.json" % tag), out)
    enc = read_checkpoint().get("encode") or {}
    enc[tag] = out
    checkpoint_extra("encode", enc)
    cp.unit_done("encode-%s" % tag)
    log("  %s  %sx%s %s fps frames=%s duration=%s pix_fmt=%s  %.2f MB  [%s]"
        % (probe.get("codec_name"), probe.get("width"), probe.get("height"),
           probe.get("r_frame_rate"), probe.get("nb_frames_reported"),
           probe.get("duration_s"), probe.get("pix_fmt"),
           out["size_bytes"] / 1e6, probe.get("probe_method")))
    return out


# ---------------------------------------------------------------------------
# stage: gates
# ---------------------------------------------------------------------------
def green_frames(tag):
    return [green_of(frame_path(tag, i)) for i in range(NFRAMES)]


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] both stages")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    photo_rows = measure_image(PHOTO, regions)
    gates = {}

    for tag in STAGES:
        sc = SCN[tag]
        f24 = (state.get("frame24") or {}).get(tag) or {}
        G = {}
        # ---- G-ENERGY ---------------------------------------------------
        chunks = (state.get("chunks") or {}).get(tag) or {}
        per_frame = []
        for k in sorted(chunks, key=int):
            for m in chunks[k].get("metrics", []):
                per_frame.append({"frame": m["frame"],
                                  "rel_err": m["energy_rel_err"],
                                  "offframe_fraction_energy":
                                      m["offframe_fraction_energy"],
                                  "deposited_total": m["energy_deposited"],
                                  "cell_energy_total": m["energy_cell_total"]})
        per_frame.sort(key=lambda x: x["frame"])
        worst = max((abs(x["rel_err"]) for x in per_frame), default=None)
        n_frames_e = len(per_frame)
        ofr = ([float(min(x["offframe_fraction_energy"] for x in per_frame)),
                float(max(x["offframe_fraction_energy"] for x in per_frame))]
               if per_frame else None)
        G["G-ENERGY"] = {
            "pass": bool(n_frames_e == NFRAMES and worst is not None
                         and worst < 1e-6),
            "frames_measured": n_frames_e, "threshold": 1e-6,
            "max_rel_err": worst,
            "max_rel_err_frame": (max(per_frame,
                                      key=lambda x: abs(x["rel_err"]))["frame"]
                                  if per_frame else None),
            "mean_rel_err": float(np.mean([abs(x["rel_err"])
                                           for x in per_frame]))
            if per_frame else None,
            "offframe_energy_fraction_range": ofr,
            "frame24_energy": f24.get("energy"),
            "frame24_deposit_detail": f24.get("deposit_detail"),
            "per_frame": per_frame,
            "normalisation_applied": 1.0,
            "criterion": ("for every frame the deposited total equals the "
                          "sum of the cell energies to a relative error "
                          "< 1e-6; the part of the projected footprint that "
                          "leaves the frame is counted separately (it is "
                          "wrapped back by the periodic rule) and never "
                          "hidden by a normalisation")}
        log("  [%s] G-ENERGY: %d frames, max rel_err %s -> %s ; off-frame "
            "energy share %s"
            % (tag, n_frames_e,
               ("%.2e" % worst) if worst is not None else "n/a",
               "PASS" if G["G-ENERGY"]["pass"] else "FAIL",
               ("%.3f-%.3f" % (ofr[0], ofr[1])) if ofr else "n/a"))

        # ---- G-LINE -----------------------------------------------------
        still = os.path.join(OUT, sc["still"])
        got = measure_image(still, regions)
        rows = []
        for g_, p_ in zip(got, photo_rows):
            rows.append({"region": g_["region"], "axis": g_["axis"],
                         "period_px": g_["period_px"],
                         "depth_this": g_["depth_e3"],
                         "depth_photo": p_["depth_e3"],
                         "ratio": g_["depth_e3"] / p_["depth_e3"],
                         "depth_hp_this": g_["depth_hp_e3"],
                         "depth_hp_photo": p_["depth_hp_e3"],
                         "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                         "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                         "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
        v18 = measure_image(V18_STILL, regions) if os.path.isfile(
            V18_STILL) else None
        G["G-LINE"] = {
            "pass": None,
            "pass_if_scored": bool(all(r["ok"] for r in rows)),
            "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
            "criterion": ("frame 24: the axis-parallel 5-15 px modulation "
                          "depth is <= the source photo's for both regions and "
                          "both axes; primary = unbiased band amplitude "
                          "(depth_e3).  DIAGNOSTIC ONLY / NOT A PASS-FAIL "
                          "GATE for this job: the brief says the single-"
                          "realization value is a reference for the "
                          "distribution diagnostic (G-LINE-stat) and must not "
                          "be used as a verdict."),
            "rows": rows, "measured": got, "photo": photo_rows,
            "v4_1_8": v18,
            "strength_this": lattice_strength(got),
            "strength_photo": lattice_strength(photo_rows),
            "strength_v4_1_8": lattice_strength(v18) if v18 else None,
            "distribution": {"this_vs_photo": paired_distribution(still, PHOTO),
                             "definition": ("per-window axis-parallel 5-15 px "
                                            "band amplitude ratio image/photo "
                                            "over the whole 256 px tile grid "
                                            "(stride 92x90, 24 windows)")},
            "regions": rinfo}
        for r_ in rows:
            log("  [%s] G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio "
                "%.3f  (hp %.3f)  %s"
                % (tag, r_["region"], r_["axis"], r_["period_px"],
                   r_["depth_this"], r_["depth_photo"], r_["ratio"],
                   r_["ratio_hp"], "OK" if r_["ok"] else "ABOVE"))

        # ---- G-PEAK -----------------------------------------------------
        pk = f24.get("peak") or {}
        sr = pk.get("std_ratio_vs_point_splat")
        pr = pk.get("p99_ratio_vs_point_splat")
        G["G-PEAK"] = {
            "pass": bool(sr is not None and pr is not None
                         and sr >= 0.9 and pr >= 0.9),
            "threshold": 0.9, "std_ratio": sr, "p99_ratio": pr,
            "std_ratio_sharp_reference": pk.get("std_ratio_vs_point_splat_sharp"),
            "p99_ratio_sharp_reference": pk.get("p99_ratio_vs_point_splat_sharp"),
            "std_ratio_det_frame": pk.get("std_ratio_vs_det_frame"),
            "p99_ratio_det_frame": pk.get("p99_ratio_vs_det_frame"),
            "C_base": pk.get("C_base"), "C_point_splat": pk.get("C_point_splat"),
            "C_det_frame_state_grid": pk.get("C_det_frame_state_grid"),
            "criterion": ("frame 24: std and p99 of C are >= 90% of the "
                          "point-splat reference on the same field and frame "
                          "grid"),
            "note": pk.get("note")}
        log("  [%s] G-PEAK: std ratio %s, p99 ratio %s (>= 0.90) -> %s"
            % (tag, ("%.3f" % sr) if sr else "n/a",
               ("%.3f" % pr) if pr else "n/a",
               "PASS" if G["G-PEAK"]["pass"] else "FAIL"))

        # ---- G-CONV -----------------------------------------------------
        cv = f24.get("convergence") or {}
        G["G-CONV"] = {
            "pass": bool(cv.get("within_15pct")),
            "threshold_rel": 0.15,
            "n_base": cv.get("n_base"), "n_half": cv.get("n_half"),
            "std_rel_diff": cv.get("std_rel_diff"),
            "p99_rel_diff": cv.get("p99_rel_diff"),
            "C_base": cv.get("C_base"), "C_half": cv.get("C_half"),
            "wall_s": cv.get("wall_s"),
            "criterion": ("frame 24: base grid vs half grid, |relative "
                          "change| of C std and p99 <= 15%"),
            "refinement": f24.get("refinement"),
            "note": cv.get("note")}
        log("  [%s] G-CONV: %s vs %s -> dstd %+.2f%%, dp99 %+.2f%% -> %s"
            % (tag, cv.get("n_base"), cv.get("n_half"),
               100 * (cv.get("std_rel_diff") or 0.0),
               100 * (cv.get("p99_rel_diff") or 0.0),
               "PASS" if G["G-CONV"]["pass"] else "FAIL"))
        rf = f24.get("refinement") or {}
        if rf:
            log("  [%s] G-CONV supplement: %s vs %s -> dstd %+.2f%%, "
                "dp99 %+.2f%%"
                % (tag, rf.get("n_refine"), cv.get("n_base"),
                   100 * (rf.get("std_rel_diff") or 0.0),
                   100 * (rf.get("p99_rel_diff") or 0.0)))

        # ---- G-MOVE -----------------------------------------------------
        frames = green_frames(tag)
        sha = [sha256_file(frame_path(tag, i)) for i in range(NFRAMES)]
        uniq = len(set(sha))
        st_ = interframe_stats(frames)
        mp4 = os.path.join(OUT, sc["mp4"])
        try:
            dg = decode_green(mp4)
            dec = interframe_stats(dg)
            dec["frames_decoded"] = int(dg.shape[0])
        except Exception as e:  # pragma: no cover
            dec = {"error": str(e)}
        G["G-MOVE"] = {
            "pass": bool(uniq == NFRAMES and st_["mean"] > 0.5),
            "unique_frames": int(uniq), "frames": NFRAMES,
            "interframe_green_png": st_, "interframe_green_mp4": dec,
            "threshold_gray": 0.5,
            "criterion": ("48/48 unique frame PNGs AND mean |d green| > 0.5 "
                          "gray")}
        log("  [%s] G-MOVE: unique %d/%d  mean|d green| %.4f -> %s"
            % (tag, uniq, NFRAMES, st_["mean"],
               "PASS" if G["G-MOVE"]["pass"] else "FAIL"))

        # ---- G-WARP -----------------------------------------------------
        warp = (state.get("warp") or {}).get(tag) or {}
        G["G-WARP"] = {"pass": None,
                       "criterion": "record only (no threshold)",
                       "g_mean_px": warp.get("g_mean_px"),
                       "g_p50_px": warp.get("g_p50_px"),
                       "g_p99_px": warp.get("g_p99_px"),
                       "per_channel": warp.get("per_channel"),
                       "depth_m": warp.get("depth_m")}
        log("  [%s] G-WARP: green mean %.3f px (record only)"
            % (tag, warp.get("g_mean_px") or float("nan")))

        # ---- G-OUT ------------------------------------------------------
        enc = (state.get("encode") or {}).get(tag) or {}
        pr_ = enc.get("probe") or {}
        fps_num = 0.0
        try:
            a_, b_ = str(pr_.get("r_frame_rate", "0/1")).split("/")
            fps_num = float(a_) / float(b_)
        except Exception:
            pass
        ok_out = (int(pr_.get("nb_frames_reported") or 0) == NFRAMES
                  and abs(float(pr_.get("duration_s") or 0) - 2.0) <= 0.05
                  and int(pr_.get("width") or 0) == WIDTH
                  and int(pr_.get("height") or 0) == HEIGHT
                  and abs(fps_num - FPS) <= 1e-6)
        G["G-OUT"] = {
            "pass": bool(ok_out), "probe": pr_, "fps": fps_num,
            "criterion": ("the MP4 contains 48 frames / 2.000 s / 736x540 / "
                          "24 fps (ffprobe's fields; obtained by header parse "
                          "plus a full decode when the ffprobe alias is "
                          "unreachable in this session)"),
            "probe_method": pr_.get("probe_method"),
            "probe_binary": pr_.get("probe_binary"),
            "size_bytes": enc.get("size_bytes"), "sha256": enc.get("sha256"),
            "mp4": enc.get("mp4")}
        log("  [%s] G-OUT: frames=%s duration=%s %sx%s fps=%s -> %s"
            % (tag, pr_.get("nb_frames_reported"), pr_.get("duration_s"),
               pr_.get("width"), pr_.get("height"), fps_num,
               "PASS" if ok_out else "FAIL"))
        G["order"] = ["G-ENERGY", "G-LINE", "G-PEAK", "G-CONV", "G-MOVE",
                      "G-WARP", "G-OUT"]
        G["scored"] = [g for g in G["order"] if G[g]["pass"] is not None]
        G["all_pass"] = all(bool(G[g]["pass"]) for g in G["scored"])
        gates[tag] = G

    # ---- G-LEGACY (global) ---------------------------------------------
    rec = read_json(STAGE0_HASHES)
    now = py_set()
    new_files = set(rec.get("new_files_of_this_job") or [])
    changed, missing, added = [], [], []
    for name, h in rec["existing"].items():
        if name not in now:
            missing.append(name)
        elif now[name] != h:
            changed.append(name)
    for name in now:
        if name not in rec["existing"] and name not in new_files:
            added.append(name)
    own = {n: {"stage0": rec["all"].get(n), "now": now.get(n),
               "changed": bool(rec["all"].get(n) != now.get(n))}
           for n in sorted(new_files)}
    imports = {"run_v10.py", "run_v15.py", "run_v15_state.py",
               "run_v15_spectrum.py", "run_v4_1_2.py", "stream_surface2.py",
               "stream_surface3.py", "run_v16_water.py", "v12_water.py",
               "det.py", "transport.py", "state.py", "image.py",
               "v5_composite.py", "v5_physics.py", "v8_physics.py",
               "transport_mesh.py"}
    ext = {}
    for name in changed:
        p = os.path.join(HERE, name)
        ext[name] = {
            "stage0_sha256": rec["all"].get(name), "now_sha256": now.get(name),
            "mtime": datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds"),
            "bytes_now": os.path.getsize(p),
            "imported_by_this_job": name in imports,
            "attribution": ("NOT this job: this driver never writes it; it "
                            "changed after the stage-0 record and is owned by "
                            "a parallel lane")}
    attributable = [n for n in changed if ext.get(n, {}).get(
        "imported_by_this_job")]
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "own_files_stable_during_run": bool(
            not any(v["changed"] for v in own.values())),
        "attributable_to_this_job": {"changed": attributable,
                                     "note": ("this job writes only its own "
                                              "new file")},
        "external_changes": ext,
        "transport_mesh_sha256_stage0": rec["all"].get("transport_mesh.py"),
        "transport_mesh_sha256_now": now.get("transport_mesh.py"),
        "criterion": ("every .py recorded at stage 0 (excluding this job's new "
                      "file) still has the same sha256; external "
                      "(parallel-lane) differences are reported with their "
                      "attribution split out")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d, added %d "
        "-> %s ; own file stable: %s"
        % (len(rec["existing"]), len(changed), len(missing), len(added),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL",
           gates["G-LEGACY"]["own_files_stable_during_run"]))

    gates["stages"] = list(STAGES)
    gates["all_pass_by_stage"] = {t: gates[t]["all_pass"] for t in STAGES}
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    for t in STAGES:
        log("  [%s] ALL: %s (%d scored gates)"
            % (t, "PASS" if gates[t]["all_pass"] else "FAIL",
               len(gates[t]["scored"])))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# stage: numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    n = {"brief": BRIEF, "workdir": HERE, "built_utc": now_iso(),
         "job": "BRIEF-V4-1-15",
         "stage0": state.get("stage0"),
         "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
         "scenes": state.get("scenes"),
         "deposit_cost": state.get("deposit_cost"),
         "deposit_path_equivalence": state.get("deposit_path_equivalence"),
         "selftest": state.get("selftest"),
         "provenance": state.get("provenance"),
         "measurement_equivalence": state.get("measurement_equivalence"),
         "frame24": state.get("frame24"),
         "warp": state.get("warp"),
         "encode": state.get("encode"),
         "gates": state.get("gates"),
         "chunks": {t: {k: {"wall_s": v.get("wall_s"),
                            "frames": sorted(v.get("frames", {})),
                            "energy_max_rel_err": v.get("energy_max_rel_err"),
                            "box_bulk_max_used": v.get("box_bulk_max_used"),
                            "metrics": v.get("metrics")}
                        for k, v in ((state.get("chunks") or {}).get(t)
                                     or {}).items()}
                   for t in STAGES},
         "sampling": {"frames_per_stage": NFRAMES, "fps": FPS,
                      "still_frame": STILL_IDX},
         "measurement": {
             "implementation": ("run_v4_1_5.lattice_metric, carried as a "
                                "verbatim text copy (textual identity + "
                                "numerical reproduction of the v4-1.7 record): "
                                "green channel, 8-bit; the brief's high pass = "
                                "image - gaussian(sigma=3 px); Hann 256x256; "
                                "2D FFT; axis band |k_perp| <= 3, 5-15 px; "
                                "period = W/k at the band peak; depth = "
                                "1000*RMS(band-passed)/mean"),
             "provenance": state.get("provenance"),
             "equivalence": state.get("measurement_equivalence"),
             "regions": rinfo,
             "comparability": ("the photo is 735x539 and the frames 736x540, "
                               "so the same window origin addresses the same "
                               "content; both are 8-bit and depth is "
                               "normalised by the window mean"),
             "known_limitation": ("depth is a 1D-marginal statistic; it does "
                                  "not distinguish a lattice from isotropic "
                                  "texture with the same axis-band power")},
         "deliverables": {"report": os.path.basename(REPORT),
                          "numbers": os.path.basename(NUMBERS),
                          "mp4": [SCN[t]["mp4"] for t in STAGES],
                          "still": [SCN[t]["still"] for t in STAGES],
                          "checkpoint": os.path.basename(CP_PATH),
                          "stage0_hashes": os.path.basename(STAGE0_HASHES)},
         "files": {}}
    paths = [CP_PATH, STAGE0_HASHES] + \
            [os.path.join(OUT, SCN[t]["mp4"]) for t in STAGES] + \
            [os.path.join(OUT, SCN[t]["still"]) for t in STAGES]
    for p in paths:
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p),
                "mtime": datetime.datetime.fromtimestamp(
                    os.path.getmtime(p)).astimezone().isoformat(
                        timespec="seconds")}
    n["still_equals_frame_0024"] = {
        t: bool(os.path.isfile(frame_path(t, STILL_IDX))
                and sha256_file(frame_path(t, STILL_IDX))
                == n["files"].get(SCN[t]["still"], {}).get("sha256"))
        for t in STAGES}
    n["transport_mesh_sha256"] = sha256_file(os.path.join(HERE,
                                                         "transport_mesh.py"))
    n["run_v4_1_15_sha256"] = sha256_file(os.path.abspath(__file__))
    return n


def _fmt(x, nd=1):
    if isinstance(x, list):
        return " / ".join(_fmt(v, nd) for v in x)
    try:
        return ("%%.%df" % nd) % float(x)
    except Exception:
        return str(x)


def _gateword(p):
    if p is None:
        return "**記録のみ**"
    return "**PASS**" if p else "**FAIL**"


def write_report(n):
    """REPORT-V4-1-15.md.  Every number printed here comes from
    numbers-v4-1-15.json (i.e. from the checkpoint), never from a constant."""
    g_all = n.get("gates") or {}
    scenes = n.get("scenes") or {}
    costs = n.get("deposit_cost") or {}
    eqv = n.get("deposit_path_equivalence") or {}
    f24 = n.get("frame24") or {}
    chunks = n.get("chunks") or {}
    warp = n.get("warp") or {}
    enc = n.get("encode") or {}
    st0 = n.get("stage0") or {}
    meas = n.get("measurement") or {}
    glue = g_all.get("G-LEGACY") or {}
    stt = n.get("selftest") or {}
    prov = n.get("provenance") or {}
    meq = n.get("measurement_equivalence") or {}
    names = {"contrast": "STAGE A 対比枠", "sub": "STAGE B サブ枠"}

    def gw(tag, key):
        return ((g_all.get(tag) or {}).get(key) or {}).get("pass")

    def fs(x, nd=2):
        try:
            return ("%%.%df" % nd) % float(x)
        except Exception:
            return "n/a"

    L = []
    A = L.append
    A("# REPORT-V4-1-15 — 面積保存の堆積を残り2つの検証サイズで確認")
    A("")
    A("Brief: `BRIEF-V4-1-15.md`（`%s`）  " % n.get("brief", ""))
    A("Workdir: `%s`  " % n.get("workdir", ""))
    A("生成: %s  " % n.get("built_utc", ""))
    A("新規ファイルは `run_v4_1_15.py` の 1 本のみ。`transport_mesh.py` は "
      "verbatim 再利用（既存 .py は一切変更していない。sha256 は "
      "`stage0-py-hashes.json` と第 9 節）。")
    A("")
    A("**このレポートは STAGE A（対比枠）と STAGE B（サブ枠）を並べて扱う。"
      "数値はすべて `numbers-v4-1-15.json`（＝checkpoint）から出している。"
      "目視判定はしていない。**")
    A("")

    # ---- 0. verdict ------------------------------------------------------
    A("## 0. 判定")
    A("")
    A("| ゲート | STAGE A（対比枠 2.4×1.76 m / 深さ 0.6 m） | "
      "STAGE B（サブ枠 4.8×3.52 m / 深さ 4.0 m） |")
    A("| --- | --- | --- |")
    for key, label in (("G-ENERGY", "G-ENERGY（堆積＝セル総和 < 1e-6）"),
                       ("G-PEAK", "G-PEAK（C の std・p99 ≥ 点スプラットの 90%）"),
                       ("G-CONV", "G-CONV（基準 vs 半分グリッド、差 ≤ 15%）"),
                       ("G-MOVE", "G-MOVE（48/48 ユニーク・平均差 > 0.5）"),
                       ("G-WARP", "G-WARP（歪み、記録のみ）"),
                       ("G-OUT", "G-OUT（48 frames / 2.000 s / 736×540 / 24 fps）")):
        A("| %s | %s | %s |"
          % (label, _gateword(gw("contrast", key)),
             _gateword(gw("sub", key))))
    A("| G-LINE（診断のみ・合否なし） | %s | %s |"
      % (_diag_line(g_all, "contrast", fs), _diag_line(g_all, "sub", fs)))
    A("| G-LEGACY（既存 .py の sha256） | %s | （両ステージ共通） |"
      % _gateword(glue.get("pass")))
    A("")
    A("合否対象（G-LINE と G-WARP を除く）の全通貨: "
      "STAGE A = **%s**、STAGE B = **%s**。"
      % ("PASS" if (g_all.get("all_pass_by_stage") or {}).get("contrast")
         else "FAIL",
         "PASS" if (g_all.get("all_pass_by_stage") or {}).get("sub")
         else "FAIL"))
    A("")

    # ---- 1. what was actually executed ----------------------------------
    A("## 1. 実行した形（ブリーフ固定値と、1 か所の逸脱）")
    A("")
    A("| 項目 | STAGE A 対比枠 | STAGE B サブ枠 |")
    A("| --- | --- | --- |")
    a, b = scenes.get("contrast") or {}, scenes.get("sub") or {}
    A("| 窓（生成パッチ＝画面） | %s × %s m | %s × %s m |"
      % (fs((a.get("window_m") or [0, 0])[0], 4),
         fs((a.get("window_m") or [0, 0])[1], 4),
         fs((b.get("window_m") or [0, 0])[0], 4),
         fs((b.get("window_m") or [0, 0])[1], 4)))
    A("| 場 | 川（stream_surface3, v4-1.8 と同一） | 風の池（preset 5, "
      "Elfouhaily-form / Arm D） |")
    A("| 深さ z（屈折） | %s m | %s m |"
      % (fs(a.get("depth_m"), 3), fs(b.get("depth_m"), 3)))
    A("| 基準グリッド | %s²（セル %s mm） | %s²（セル %s mm） |"
      % (a.get("grid"), fs(1000.0 * (a.get("cell_m") or 0), 3),
         b.get("grid"), fs(1000.0 * (b.get("cell_m") or 0), 3)))
    A("| G-CONV の半分グリッド | %s² | %s² |"
      % ((f24.get("contrast") or {}).get("convergence", {}).get("n_half"),
         (f24.get("sub") or {}).get("convergence", {}).get("n_half")))
    A("| Nyquist / 1.8 メイン（347 cells/m） | %s 倍 | %s 倍 |"
      % (fs((a.get("deposit") or {}).get("nyquist_ratio_vs_main"), 2),
         fs((b.get("deposit") or {}).get("nyquist_ratio_vs_main"), 2)))
    A("| 媒体 / 位相関数 | %s / %s (n=%s, μ=%s) | %s / %s (n=%s, μ=%s) |"
      % ((a.get("water_preset")), ((a.get("phase") or {}).get("label")),
         ((a.get("phase") or {}).get("n")), ((a.get("phase") or {}).get("mu")),
         (b.get("water_preset")), ((b.get("phase") or {}).get("label")),
         ((b.get("phase") or {}).get("n")), ((b.get("phase") or {}).get("mu"))))
    A("| 合成（凍結形） | A=%s / Blow σ=%s px / E=%s / Q=%s | 同左 |"
      % ((a.get("composite") or {}).get("A"),
         (a.get("composite") or {}).get("Blow_sigma_px"),
         (a.get("composite") or {}).get("E"),
         (a.get("composite") or {}).get("Q")))
    A("| 太陽 | alt %s° | alt %s° |"
      % (fs((a.get("env") or {}).get("sun_altitude_deg"), 2),
         fs((b.get("env") or {}).get("sun_altitude_deg"), 2)))
    A("| 出力 | 736×540 / 24 fps / 48 frames | 同左 |")
    A("")
    A("**倍率 1。** 生成パッチの一辺 L がそのまま画面の窓で、"
      "引き伸ばしは存在しない（L_img/L_gen = 1）。窓の高さは "
      "L·540/736 で、2.4 m 枠なら 1.7609 m（ブリーフの 1.76）、"
      "4.8 m 枠なら 3.5217 m（同 3.52）。736:540 の画素アスペクト非等方は "
      "v2 のフレーミングそのもの。")
    A("")
    A("### 1.1 STAGE B のグリッドだけがブリーフの目安から外れている（開示）")
    A("")
    bb = b.get("deposit") or {}
    A("ブリーフの規則は「基準グリッドの Nyquist 波数が 1.8 メイン実行"
      "（512² / 1.4755 m = 347 cells/m）と同等の物理帯域を覆うこと（目安 "
      "1024²）」。STAGE A は 1024² でこの規則を満たす（%s 倍 = %s cells/m）。"
      "STAGE B は **満たしていない**（%s 倍 = %s cells/m）。"
      % (fs((a.get("deposit") or {}).get("nyquist_ratio_vs_main"), 2),
         fs((a.get("deposit") or {}).get("nyquist_cells_per_m"), 0),
         fs(bb.get("nyquist_ratio_vs_main"), 2),
         fs(bb.get("nyquist_cells_per_m"), 0)))
    A("")
    A("規則どおりなら STAGE B は 2048²（1666² 以上が必要）、目安どおりなら "
      "1024² である。どちらも **verbatim の堆積では実行できない**ことが"
      "測って分かった（stage 0 で、各グリッドの frame 24 の投影だけを回して"
      "bbox 面積を数え、実測コストを掛けたもの）:")
    A("")
    A("| ステージ | グリッド | セル | Nyquist | 投影の重なり Σ\\|A\\|/画面 | "
      "Σ bbox px | 推定 1 フレーム | 推定 48 フレーム |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        cd = costs.get(t) or {}
        for key in sorted((cd.get("by_grid") or {}), key=int):
            d = cd["by_grid"][key]
            mark = ""
            if int(key) == int(cd.get("base_grid") or -1):
                mark = " **← 採用**"
            A("| %s | %s²%s | %.2f mm | %.0f /m | %.1f | %.2e | %.0f s | %.1f h |"
              % (names[t], key, mark, float(d.get("cell_mm") or 0),
                 float(d.get("nyquist_cells_per_m") or 0),
                 float(d.get("fold_ratio_vs_frame") or 0),
                 float(d.get("sum_bbox_px") or 0),
                 float(d.get("est_wall_s_one_frame_bulk_path") or 0),
                 float(d.get("est_wall_s_48_frames_bulk_path") or 0) / 3600.0))
    A("")
    A("推定は「Σ bbox × そのステージで実測した µs/box」（A = %.2f、B = %.2f）"
      "で、実測値そのものではない。"
      "512² を採用した理由は実行可能性だけで、物理的な良し悪しではない。"
      "**この 1 点が STAGE B の結論の適用範囲を狭める**（第 12 節）。"
      % (float((eqv.get("contrast") or {}).get("us_per_box_bulk") or 0.0),
         float((eqv.get("sub") or {}).get("us_per_box_bulk") or 0.0)))
    A("")
    A("STAGE B では 1024² の frame 24 を 1 枚だけ追加で回し、"
      "「基準 512² が 1024² に対してどれだけ動くか」を上向きにも測った"
      "（G-CONV の補足、第 5.1 節）。")
    A("")

    # ---- 1.2 wall time ---------------------------------------------------
    A("### 1.2 実行時間（すべて foreground、stage ごとに jobkit checkpoint）")
    A("")
    A("| ステージ | frame 24 | chunk 0 | chunk 1 | warp | encode |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        ch = chunks.get(t) or {}
        A("| %s | %s s | %s s | %s s | %s s | %s s |"
          % (names[t],
             fs((f24.get(t) or {}).get("wall_s"), 0),
             fs((ch.get("0") or {}).get("wall_s"), 0),
             fs((ch.get("1") or {}).get("wall_s"), 0),
             fs((warp.get(t) or {}).get("wall_s"), 0),
             fs((enc.get(t) or {}).get("wall_s"), 0)))
    A("")
    A("`checkpoint.json` のユニット履歴には、stage-0 が 2 回、gates と report が"
      "複数回入っている。理由は第 9 節の注記のとおり（stage-0 の後に"
      "レポート整形だけでドライバを編集したため、stage-0 をもう一度回して"
      "入力の同一性を再確認し、gates/report を再生成した）。"
      "**フレーム画像・MP4 は 1 回だけ生成しており、再生成していない**"
      "（chunk の sha256 が checkpoint に残っている）。")
    A("")

    # ---- 2. deposit ------------------------------------------------------
    A("## 2. 堆積（`transport_mesh.py`、verbatim）")
    A("")
    A("* **射影**: 水面セル（四辺形 → 2 三角形）の角を、既存輸送段と同じ Snell "
      "計算（単位太陽方向, η = 1/1.33434, 平面平行 z）で床へ落とす。セル角の "
      "h/gx/gy は 4 texel の双一次平均（周期）。")
    A("* **分配**: 射影三角形とフレーム画素の**重なり面積**に比例して "
      "セルのエネルギーを配る。面積は画素正方形への解析的クリップ"
      "（Sutherland–Hodgman）で**厳密**。カーネル・ジッタ・帯域制限は"
      "**一切なし**、平均 1 への強制正規化もなし。")
    A("* **セル重み**: `w_cell = W·H/N_s²`（既存 C レーンの divisor 規約の"
      "長方形フレーム版）。A: %s、B: %s。"
      % (fs((a.get("deposit") or {}).get("w_cell"), 6),
         fs((b.get("deposit") or {}).get("w_cell"), 6)))
    A("* **周期**: 合成段は C を mod 付き (u,v) で引くので、着地は "
      "`i mod W, j mod H`。unwrap 座標で枠外に出た分は別枠で数え、"
      "正規化には一切使わない。")
    A("* **セルフテスト**: 毎回 stage 0 で実行。**%s**（%s checks ok / %s ng）。"
      % ("PASS" if stt.get("ok") else "FAIL", stt.get("n_checks_ok"),
         stt.get("n_checks_ng")))
    A("* **測定の来歴**: G-LINE の統計は run_v4_1_5 の実装を逐語コピーし、"
      "ソーステキストを文字単位で比較（すべて一致 = %s）、さらに v4-1-7 の"
      "記録値を数値再現した（最大差 %s）。"
      % ("PASS" if prov.get("pass") else "FAIL",
         ("%.2e" % meq["max_abs_diff_over_all_rows"])
         if meq.get("max_abs_diff_over_all_rows") is not None else "n/a"))
    A("")
    A("### 2.1 1 つの性能スイッチ（STAGE B のみ、方法は変えない）")
    A("")
    A("`transport_mesh` の `BOX_BULK_MAX` は「整数 bbox が何画素までの三角形を"
      "ベクトル化経路で処理するか」を決めるだけの定数で、両経路とも面積は"
      "解析的（近似なし）、違いは和の順序だけである。深さ 4.0 m の STAGE B は"
      "投影が大きく折り返すため、既定値 16 では 9 割の三角形が 1 本ずつ処理される"
      "経路に落ちる。STAGE B ではこの定数を実行時に上げ、ファイルは触っていない。")
    A("")
    A("| ステージ | 経路 | 三角形数 | 平均 bbox | max\\|一括−逐次\\| | 一括 µs/box | 逐次 µs/box |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        e = eqv.get(t) or {}
        A("| %s | %s | %s | %s px | %.2e | %s | %s |"
          % (names[t], "一括のみ（BOX_BULK_MAX を上げた）" if t == "sub"
             else "既定（16、変更なし）", e.get("n_triangles"),
             fs(e.get("bbox_px_mean"), 1), e.get("max_abs_difference") or 0.0,
             fs(e.get("us_per_box_bulk"), 2), fs(e.get("us_per_box_tail"), 2)))
    A("")
    A("一括経路と逐次経路の差は **%.1e**（実三角形 %s 本で測定）。"
      "つまりこれは速度の都合であって、堆積の数学は同一である。"
      % (max([abs((eqv.get(t) or {}).get("max_abs_difference") or 0.0)
              for t in STAGES] or [0.0]),
         (eqv.get("sub") or {}).get("n_triangles")))
    A("")

    # ---- 3. G-ENERGY -----------------------------------------------------
    A("## 3. G-ENERGY — フレームごとの堆積会計")
    A("")
    A("`|堆積総量 − セルエネルギー総和| / セルエネルギー総和 < 1e-6`／ステージ。"
      "正規化で隠していない: 平均 1 は `w_cell = W·H/N_s²` の定義の帰結で、"
      "測定値による係数は掛けていない。")
    A("")
    A("| ステージ | 測ったフレーム | 最大相対誤差 | その frame | 平均\\|誤差\\| | 枠外（エネルギー）の範囲 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        ge = ((g_all.get(t) or {}).get("G-ENERGY") or {})
        ofr = ge.get("offframe_energy_fraction_range") or [None, None]
        A("| %s | %s | %.2e | %s | %.2e | %s 〜 %s |"
          % (names[t], ge.get("frames_measured"), float(ge.get("max_rel_err") or 0),
             ge.get("max_rel_err_frame"), float(ge.get("mean_rel_err") or 0),
             fs(ofr[0], 4) if ofr[0] is not None else "n/a",
             fs(ofr[1], 4) if ofr[1] is not None else "n/a"))
    A("")
    for t in STAGES:
        ge = ((g_all.get(t) or {}).get("G-ENERGY") or {})
        rows = ge.get("per_frame") or []
        A("### 3.%d %s — 48 フレームの表" % (1 if t == "contrast" else 2,
                                          names[t]))
        A("")
        A("| frame | 堆積総量 | セル総和 | 相対誤差 | 枠外(エネルギー) |")
        A("| --- | --- | --- | --- | --- |")
        for r in rows:
            A("| %d | %.3f | %.3f | %.2e | %.4f |"
              % (r.get("frame", -1), float(r.get("deposited_total") or 0),
                 float(r.get("cell_energy_total") or 0),
                 float(r.get("rel_err") or 0),
                 float(r.get("offframe_fraction_energy") or 0)))
        A("")

    # ---- 4. G-PEAK -------------------------------------------------------
    A("## 4. G-PEAK — ピークが削られていないか（同一フィールド・frame 24）")
    A("")
    A("| ステージ | C | mean | std | p99 | max |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        pk = ((f24.get(t) or {}).get("peak") or {})
        for key, lab in (("C_base", "面積堆積（本作, 736×540）"),
                         ("C_point_splat", "点スプラット（同一 736×540, Q=25）"),
                         ("C_det_frame_state_grid", "det_frame（状態グリッド）")):
            c = pk.get(key) or {}
            A("| %s | %s | %.6f | %.4f | %.3f | %.2f |"
              % (names[t], lab, float(c.get("mean") or 0),
                 float(c.get("std") or 0), float(c.get("p99") or 0),
                 float(c.get("max") or 0)))
    A("")
    for t in STAGES:
        pk = ((f24.get(t) or {}).get("peak") or {})
        A("* %s: std 比 **%s**（点スプラット）/ %s（det_frame）、"
          "p99 比 **%s** / %s。"
          % (names[t], fs(pk.get("std_ratio_vs_point_splat"), 3),
             fs(pk.get("std_ratio_vs_det_frame"), 3),
             fs(pk.get("p99_ratio_vs_point_splat"), 3),
             fs(pk.get("p99_ratio_vs_det_frame"), 3)))
    A("")
    A("参照は**本編には使わない**比較材料。面積堆積の射影は太陽円盤の中心方向 "
      "1 本で計算している（ブリーフの「滲み幅もジッタも要らない」に従う）。"
      "点スプラットは Q=25 本の画素フィルタなので、円盤平均を省いた分だけ"
      "ピーク比較が有利になりうる。**この効果は未測定**である（Q=1 と Q=25 の"
      "点スプラット比較は本作では取っていない）。深さ 4.0 m では円盤半径 "
      "0.00465 rad × 4.0 m ≈ 18.6 mm ≈ 2.9 px になるので、STAGE B の"
      "この非対称は STAGE A（0.6 m で ≈ 0.4 px）より大きい。")
    A("")

    # ---- 5. G-CONV -------------------------------------------------------
    A("## 5. G-CONV — 基準グリッド vs 半分グリッド（同じ 736×540 へ堆積）")
    A("")
    A("| ステージ | 水面 | C mean | C std | C p99 | C max |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        cv = ((f24.get(t) or {}).get("convergence") or {})
        for key, lab in (("C_base", "基準 %s²" % cv.get("n_base")),
                         ("C_half", "半分 %s²" % cv.get("n_half"))):
            c = cv.get(key) or {}
            A("| %s | %s | %.6f | %.4f | %.3f | %.2f |"
              % (names[t], lab, float(c.get("mean") or 0),
                 float(c.get("std") or 0), float(c.get("p99") or 0),
                 float(c.get("max") or 0)))
    A("")
    A("| ステージ | std の差 | p99 の差 | 判定（≤15%） |")
    A("| --- | --- | --- | --- |")
    for t in STAGES:
        cv = ((f24.get(t) or {}).get("convergence") or {})
        A("| %s | %+.2f%% | %+.2f%% | %s |"
          % (names[t], 100 * float(cv.get("std_rel_diff") or 0),
             100 * float(cv.get("p99_rel_diff") or 0),
             _gateword(cv.get("within_15pct"))))
    A("")
    A("### 5.1 STAGE B の上向きの確認（補足）")
    A("")
    rfb = ((f24.get("sub") or {}).get("refinement") or {})
    if rfb:
        A("基準 512² を 1024²（ブリーフの目安）と比べた 1 フレーム:")
        A("")
        A("| 水面 | C std | C p99 | 512² からの差 |")
        A("| --- | --- | --- | --- |")
        c0 = ((f24.get("sub") or {}).get("convergence") or {}).get("C_base") or {}
        cr = rfb.get("C_refine") or {}
        A("| 512²（基準） | %.4f | %.3f | — |"
          % (float(c0.get("std") or 0), float(c0.get("p99") or 0)))
        A("| 1024² | %.4f | %.3f | std %+.2f%%, p99 %+.2f%% |"
          % (float(cr.get("std") or 0), float(cr.get("p99") or 0),
             100 * float(rfb.get("std_rel_diff") or 0),
             100 * float(rfb.get("p99_rel_diff") or 0)))
        A("")
        rfb_std = 100.0 * float(rfb.get("std_rel_diff") or 0.0)
        if abs(rfb_std) > 15.0:
            A("512² と 1024² の std の差は **%+.1f%%**（15%% 超）。"
              "%s"
              "つまり「半分グリッドとの一致」は合格でも、"
              "**格子収束を意味しない**。この向きの差は第 12 節の限界に数える。"
              % (rfb_std,
                 "512² は 1024² より **分散が小さい**（細かくすると C が立つ）。"
                 if rfb_std < 0 else
                 "512² は 1024² より **分散が大きい**（粗くすると C が立つ）。"))
        else:
            A("512² と 1024² の std の差は **%+.1f%%** で 15%% 以内。"
              "下向き（256²）・上向き（1024²）の両方で統計が動かない。"
              % rfb_std)
        A("")

    # ---- 6. G-LINE -------------------------------------------------------
    A("## 6. G-LINE — 軸平行周期（5–15 px 帯）／診断のみ")
    A("")
    A("**この節は合否判定ではない。** ブリーフが明記するとおり、単一 "
      "realization の値は分布診断（G-LINE-stat）の参考値であり、"
      "合否には使わない。表の「元写真以下か」は参考列として出す。")
    A("")
    A("測定法は v4-1.5 と同一（緑・8 bit・写真で位置を固定した 256×256 窓・"
      "Hann・2D FFT・|x 軸は |ky| ≤ 3 で 5–15 px|）。主量は前段ハイパスなしの "
      "`depth_e3`、`depth_hp_e3` は 3 px ハイパス後の参考。")
    A("")
    for t in STAGES:
        gl = ((g_all.get(t) or {}).get("G-LINE") or {})
        A("### 6.%d %s（全体強度 %s、写真 %s、v4-1.8 は %s）"
          % (1 if t == "contrast" else 2, names[t],
             fs(gl.get("strength_this"), 1), fs(gl.get("strength_photo"), 1),
             fs(gl.get("strength_v4_1_8"), 1)))
        A("")
        A("| region | axis | period px | depth_e3 本作 | depth_e3 写真 | 比 | "
          "depth_hp 本作 / 写真 | 比(hp) | 参考 |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for r in (gl.get("rows") or []):
            A("| %s | %s | %.2f | %.2f | %.2f | %.3f | %.2f / %.2f | %.3f | %s |"
              % (r.get("region"), r.get("axis"), float(r.get("period_px") or 0),
                 float(r.get("depth_this") or 0),
                 float(r.get("depth_photo") or 0),
                 float(r.get("ratio") or 0),
                 float(r.get("depth_hp_this") or 0),
                 float(r.get("depth_hp_photo") or 0),
                 float(r.get("ratio_hp") or 0),
                 "写真以下" if r.get("ok") else "写真超"))
        A("")
        dist = ((gl.get("distribution") or {}).get("this_vs_photo") or {})
        if dist:
            A("全 24 窓の分布（本作/写真）: ")
            A("")
            A("| 軸 | 中央値 | p10 | p90 | 1 超の割合 | 最大 |")
            A("| --- | --- | --- | --- | --- | --- |")
            for axis in ("x", "y"):
                d = dist.get(axis) or {}
                if d:
                    A("| %s | %.3f | %.3f | %.3f | %.3f | %.3f |"
                      % (axis, float(d.get("ratio_median") or 0),
                         float(d.get("ratio_p10") or 0),
                         float(d.get("ratio_p90") or 0),
                         float(d.get("frac_above_1") or 0),
                         float(d.get("max_ratio") or 0)))
            A("")

    # ---- 7. G-MOVE + seam ------------------------------------------------
    A("## 7. G-MOVE とループの継ぎ目")
    A("")
    A("| ステージ | ユニーク | 平均\\|Δ緑\\| PNG | 最小 | 最大 | 平均\\|Δ緑\\| MP4 | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        gm = ((g_all.get(t) or {}).get("G-MOVE") or {})
        ip = gm.get("interframe_green_png") or {}
        im_ = gm.get("interframe_green_mp4") or {}
        A("| %s | %s/%s | %.2f | %.2f | %.2f | %.2f | %s |"
          % (names[t], gm.get("unique_frames"), gm.get("frames"),
             float(ip.get("mean") or 0), float(ip.get("min") or 0),
             float(ip.get("max") or 0), float(im_.get("mean") or 0),
             _gateword(gm.get("pass"))))
    A("")
    for t in STAGES:
        s = (scenes.get(t) or {})
        seam = s.get("seam") or {}
        loop = s.get("loop") or {}
        if not seam:
            continue
        A("* **%s の継ぎ目**（t=0 と t=2 s）: max|Δh| = %.4f m、RMS Δh = "
          "%.4f m（場の RMS %.4f m、比 %.2f）。"
          % (names[t], float(seam.get("max_abs_dh_m") or 0),
             float(seam.get("rms_dh_m") or 0),
             float(seam.get("h_rms_m") or 0),
             float(seam.get("seam_over_h_rms") or 0)))
        if t == "contrast":
            A("  ループ規則 n = %s、T_mode = L/Us = %.4f s（ピースはその "
              "%.3f 周期分）。**2 秒で厳密には閉じない**。"
              % (loop.get("n_loops_rule"), float(loop.get("T_loop_mode_s") or 0),
                 float(loop.get("n_cycles_in_piece") or 0)))
        else:
            modes = loop.get("modes") or {}
            lo = modes.get("lo") or {}
            hi = modes.get("hi") or {}
            A("  風レーンは移流ループを持たない。2 s で閉じるには全モードで "
              "ω(k)·T = 2πn が必要で、実測の最大位相ずれは lo %s rad / "
              "hi %s rad。**2 秒で厳密には閉じない**。"
              % (fs(lo.get("max_phase_dev_at_T_rad"), 3),
                 fs(hi.get("max_phase_dev_at_T_rad"), 3)))
    A("")

    # ---- 8. G-WARP -------------------------------------------------------
    A("## 8. G-WARP（合否なし・記録）")
    A("")
    A("| ステージ | 深さ | 緑 平均 px | p50 | p99 | max px |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in STAGES:
        wv = warp.get(t) or {}
        per = (wv.get("per_channel") or {}).get("g") or {}
        A("| %s | %s m | %.3f | %.3f | %.3f | %.3f |"
          % (names[t], fs(wv.get("depth_m"), 2), float(wv.get("g_mean_px") or 0),
             float(per.get("px_p50") or 0), float(per.get("px_p99") or 0),
             float(per.get("px_max") or 0)))
    A("")
    A("歪みレーンは既存のまま（屈折シフト `−n_xy·z·(IOR−1)/IOR` を画素へ換算）。"
      "同じ傾斜でも z が 6.7 倍違うので、STAGE B の歪みは STAGE A の約 6.7 倍になる。")
    A("")

    # ---- 9. G-LEGACY -----------------------------------------------------
    A("## 9. G-LEGACY — 既存 .py の同一性と並行レーンの差分")
    A("")
    ph = st0.get("photo") or {}
    A("STAGE 0 で確認した必須入力（ブリーフの 6 項目）:")
    A("")
    A("| 項目 | パス | 存在 |")
    A("| --- | --- | --- |")
    for p in (st0.get("brief_list") or []):
        A("| `%s` | `%s` | %s |"
          % (os.path.basename(p), p,
             "はい" if os.path.isfile(p) else "**いいえ**"))
    A("")
    A("**土台写真の同一性**（ブリーフの「同じファイル」要求の検証）: "
      "本ジョブが使うのは `%s`（sha256 `%s`）。"
      % (os.path.basename(ph.get("path") or ""), str(ph.get("sha256"))[:32]))
    A("")
    A("| 記録 | そのドライバが宣言する PHOTO | 本ジョブと同じ | "
      "その記録が保存している sha256 | 一致 |")
    A("| --- | --- | --- | --- | --- |")
    for k, v in (ph.get("provenance") or {}).items():
        drv = ", ".join("%s: %s" % (dn, ("一致" if dd.get("equals_this_job")
                                         else ("—（定数参照）"
                                               if dd.get("equals_this_job")
                                               is None else "不一致")))
                        for dn, dd in (v.get("drivers") or {}).items())
        rec = v.get("recorded") or {}
        A("| %s | %s | %s | `%s` | %s |"
          % (k, drv,
             "はい" if any(dd.get("equals_this_job")
                           for dd in (v.get("drivers") or {}).values())
             else "いいえ",
             str(rec.get("sha256"))[:24] if rec.get("sha256") else
             (rec.get("note") or rec.get("error") or "—"),
             ("—" if rec.get("matches_this_job") is None
              else rec.get("matches_this_job"))))
    A("")
    wl = st0.get("wind_lane") or {}
    A("**風レーンの生成器**: %s。spectrum_source = `%s`、arm = `%s`、"
      "Ω_c 警告閾値 = %s rad/s。プリセット5 の U = 6 / F = 100 は "
      "`scene_for` にそのまま渡している（代用なし）。"
      % (wl.get("generator"), wl.get("spectrum_source_default"),
         wl.get("merge_arm"), wl.get("warn_omega_c_threshold")))
    A("")
    A("stage-0 で記録した既存 .py %s 本: 変更 %s、欠落 %s。"
      "本ジョブが書くのは新規 1 本だけ。`transport_mesh.py` は stage-0 から "
      "%s（`%s`）。"
      % (glue.get("recorded_existing"), len(glue.get("changed") or []),
         len(glue.get("missing") or []),
         "バイト同一" if (glue.get("transport_mesh_sha256_stage0")
                          == glue.get("transport_mesh_sha256_now"))
         else "変化あり",
         str(glue.get("transport_mesh_sha256_now"))[:24]))
    A("")
    ext = glue.get("external_changes") or {}
    A("並行レーンの外部差分: %s。"
      % ("なし" if not ext else
         ", ".join("%s（mtime %s、本ジョブは import %s）"
                   % (k, v.get("mtime"),
                      "する" if v.get("imported_by_this_job") else "しない")
                   for k, v in ext.items())))
    A("")
    own = (glue.get("own_new_file_hashes") or {}).get("run_v4_1_15.py") or {}
    if own.get("changed"):
        A("**本ジョブ自身の新規ファイルの差分（正直な注記）**: `run_v4_1_15.py` は "
          "stage-0 の記録（`%s`）と現在（`%s`）でハッシュが違う。"
          "stage-0 の後に編集したのは **レポート生成部と、G-LINE を "
          "`pass` から `pass_if_scored` に変えた記録方法だけ**で、"
          "測定経路（`SCN` の定数、`build_scene`、`render_frame`、"
          "`TM.mesh_C` の呼び方、`BOX_BULK_MAX` の扱い）は stage-0 以降 "
          "1 文字も変えていない。48 フレーム・全ゲート・本レポートの数値は"
          "すべて最終ファイルで生成しており、stage-0 の記録だけが"
          "それより前の版を指している。"
          % (str(own.get("stage0"))[:16], str(own.get("now"))[:16]))
        A("")
    added = glue.get("new_since_stage0_not_of_this_job") or []
    if added:
        A("**stage-0 の後に作業ディレクトリへ増えた .py（本ジョブのものではない）**: "
          "%s。うち `run_v4_1_15_d15.py` は別ブリーフ"
          "（BRIEF-V4-1-15-D15、深さだけ 1.5 m に変えた 1 枚）のドライバで、"
          "本ドライバを import して読み、出力は `out/v4-1-15-extra/` に書く。"
          "`out/v4-1-15/` には書かない（本ジョブの成果物は混ざっていない）。"
          % ", ".join("`%s`" % a for a in added))
        A("")

    # ---- 10. G-OUT -------------------------------------------------------
    A("## 10. G-OUT")
    A("")
    A("| 項目 | STAGE A | STAGE B |")
    A("| --- | --- | --- |")
    pa = ((g_all.get("contrast") or {}).get("G-OUT") or {}).get("probe") or {}
    pb = ((g_all.get("sub") or {}).get("G-OUT") or {}).get("probe") or {}
    A("| codec / pix_fmt | %s / %s | %s / %s |"
      % (pa.get("codec_name"), pa.get("pix_fmt"),
         pb.get("codec_name"), pb.get("pix_fmt")))
    A("| 解像度 | %sx%s | %sx%s |"
      % (pa.get("width"), pa.get("height"), pb.get("width"), pb.get("height")))
    A("| frames / duration | %s / %s s | %s / %s s |"
      % (pa.get("nb_frames_reported"), pa.get("duration_s"),
         pb.get("nb_frames_reported"), pb.get("duration_s")))
    A("| fps | %s | %s |" % (pa.get("r_frame_rate"), pb.get("r_frame_rate")))
    A("| bytes | %s | %s |"
      % (((g_all.get("contrast") or {}).get("G-OUT") or {}).get("size_bytes"),
         ((g_all.get("sub") or {}).get("G-OUT") or {}).get("size_bytes")))
    A("| probe の方法 | %s | %s |"
      % (pa.get("probe_method"), pb.get("probe_method")))
    A("")
    A("**環境の注記**: このセッションでは PATH 上の `ffmpeg` / `ffprobe` は"
      "他ユーザープロファイル内の 0 バイト WinGet エイリアスで、CreateProcess が"
      "アクセス拒否（WinError 5）を返す。リポジトリ同梱の ffmpeg 9.0.2"
      "（`_vendor/imageio_ffmpeg`）でエンコードし、G-OUT の値は"
      "「ヘッダ解析 + 全フレームの rgb24 デコードによる枚数」で取得した。"
      "フィールド名は ffprobe のものをそのまま使っている。")
    A("")

    # ---- 11. deliverables ------------------------------------------------
    A("## 11. 成果物（すべて `out/v4-1-15/`）")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for name, rec in sorted((n.get("files") or {}).items()):
        A("| `%s` | %s | `%s` |" % (name, rec.get("bytes"),
                                   str(rec.get("sha256"))[:32]))
    A("| `numbers-v4-1-15.json` | （本文の数値の実体） | — |")
    A("")
    A("`frame-0024-*.png` と `frames-*/frame-0024.png` の同一性: "
      "A = %s, B = %s。`transport_mesh.py` sha256 `%s`、"
      "`run_v4_1_15.py` sha256 `%s`。"
      % ((n.get("still_equals_frame_0024") or {}).get("contrast"),
         (n.get("still_equals_frame_0024") or {}).get("sub"),
         str(n.get("transport_mesh_sha256"))[:32],
         str(n.get("run_v4_1_15_sha256"))[:32]))
    A("")
    A("注: `checkpoint.json` はこのレポートを書いた後に `report` ユニットを"
      "追記するので、上の sha256 はレポート生成時点のものである。"
      "`numbers-v4-1-15.json` は本文の数値の実体（機械可読）。")
    A("")

    # ---- 12. disclosures -------------------------------------------------
    A("## 12. 開示・限界")
    A("")
    A("1. **STAGE B のグリッドはブリーフの目安より粗い（512²）。** 目安 1024² は "
      "1 フレーム約 %.0f 分（48 フレームで約 %.0f 時間）、規則どおりの 2048² は "
      "1 フレーム約 %.0f 分（約 %.0f 時間）で、verbatim の堆積では実行できない。"
      "512² は「実行できる最大」であって「十分細かい」ではない。Nyquist は "
      "%.0f cells/m で 1.8 メイン実行の %.2f 倍しかない。したがって STAGE B の"
      "数値は **この窓サイズ・この深さ・この粗い格子での結果**であり、"
      "細かい格子での挙動を代表しない（第 5.1 節の上向き比較を併読）。"
      % (float(((costs.get("sub") or {}).get("by_grid") or {})
               .get("1024", {}).get("est_wall_s_one_frame_bulk_path") or 0)
         / 60.0,
         float(((costs.get("sub") or {}).get("by_grid") or {})
               .get("1024", {}).get("est_wall_s_48_frames_bulk_path") or 0)
         / 3600.0,
         float(((costs.get("sub") or {}).get("by_grid") or {})
               .get("2048", {}).get("est_wall_s_one_frame_bulk_path") or 0)
         / 60.0,
         float(((costs.get("sub") or {}).get("by_grid") or {})
               .get("2048", {}).get("est_wall_s_48_frames_bulk_path") or 0)
         / 3600.0,
         float(bb.get("nyquist_cells_per_m") or 0),
         float(bb.get("nyquist_ratio_vs_main") or 0)))
    A("2. **STAGE B の G-CONV は FAIL（下向き 512² vs 256² で +%.1f%%）。** "
      "上向き（512² vs 1024²）でも std が %+.1f%% 動く。つまり深さ 4.0 m の "
      "C 統計は格子依存で、**この構成の C の水準は格子に収束していない**。"
      "「半分グリッドと一致しない」ことは堆積の誤りではなく、"
      "折り返した写像を有限の格子で切っていることの帰結である（第 5.1 節）。"
      % (100.0 * float((((f24.get("sub") or {}).get("convergence") or {})
                        .get("std_rel_diff")) or 0),
         100.0 * float((((f24.get("sub") or {}).get("refinement") or {})
                        .get("std_rel_diff")) or 0)))
    A("3. **深さ 4.0 m では投影の折り返しが大きい。** frame 24 の "
      "Σ|A|/画面 = %s 倍（STAGE A は %s 倍）。C は平均 1 を保つが、"
      "局所の多重被覆は C > 1 として現れる（クリップも平滑化もしない）。"
      "STAGE B の C の p99 = %s、max = %s。"
      % (fs(((costs.get("sub") or {}).get("by_grid") or {})
            .get(str((costs.get("sub") or {}).get("base_grid")), {})
            .get("fold_ratio_vs_frame"), 1),
         fs(((costs.get("contrast") or {}).get("by_grid") or {})
            .get(str((costs.get("contrast") or {}).get("base_grid")), {})
            .get("fold_ratio_vs_frame"), 1),
         fs(((((f24.get("sub") or {}).get("peak") or {})
               .get("C_base") or {}).get("p99")), 2),
         fs(((((f24.get("sub") or {}).get("peak") or {})
               .get("C_base") or {}).get("max")), 2)))
    A("4. **円盤求積**: 射影の幾何は太陽円盤の中心方向 1 本。深さ 4.0 m では"
      "円盤半径が ≈ 2.9 px に相当するので、この省略は STAGE B で最も効く。"
      "G-PEAK の比較をどれだけ有利にしているかは未測定。")
    A("5. **BOX_BULK_MAX の実行時変更（STAGE B のみ）**: 方法は変えない"
      "（両経路とも解析的な重なり面積、差は和の順序）。実測差 %.1e。" %
      max([abs((eqv.get(t) or {}).get("max_abs_difference") or 0.0)
           for t in STAGES] or [0.0]))
    A("6. **画素アスペクト**: x と y で 1 m あたりの画素数が違う（736:540）。"
      "これは v2 のフレーミングで、倍率 1 は L_img/L_gen = 1 の意味。")
    A("7. **C の受け渡し**: 凍結合成器は (u,v) の正方 C を要求するため、"
      "736×540 の堆積を 736×736 へ線形 y 補間で再標本化している"
      "（純粋な再標本化、再正規化も追加のぼかしもしない）。")
    A("8. **ループ**: 両ステージとも 2 秒で厳密には閉じない（第 7 節の数値）。")
    A("9. **測定の限界**: depth は 1 次元周辺統計で、格子と等方的テクスチャを"
      "区別しない。C の統計は 1 フレーム・1 場のもので、シーンの一般性は"
      "主張しない。")
    A("")

    # ---- 13. what this shows --------------------------------------------
    A("## 13. この結果が言えること / 言えないこと")
    A("")
    A("* **言える**: 面積保存の堆積は、窓を 1.6 倍（2.4 m）に広げても、"
      "エネルギーを厳密に保存し（最大 %s）、ピークを点スプラット参照以上に"
      "残す（std 比 %s）。STAGE A の結果は v4-1.8 の結論を"
      "「窓だけを変えた枠」で再現している。"
      % (("%.1e" % float((((g_all.get("contrast") or {})
                           .get("G-ENERGY") or {}).get("max_rel_err")) or 0)),
         fs((((g_all.get("contrast") or {}).get("G-PEAK") or {})
             .get("std_ratio")), 2)))
    A("* **言える**: 深さ 4.0 m の風の池でも、堆積のエネルギー会計は同じ精度で"
      "閉じる（最大 %s）。"
      % ("%.1e" % float((((g_all.get("sub") or {}).get("G-ENERGY") or {})
                         .get("max_rel_err")) or 0)))
    A("* **言えない**: STAGE B が細かい格子でも同じ統計になること。"
      "基準 512² はブリーフの目安の半分で、上向き比較（1024²）では C の std が "
      "%+.1f%% 動く。G-CONV 自体も下向き（256²）で **%s** なので、"
      "深さ 4.0 m の C 統計は **格子に収束していない**。"
      "STAGE B で言えるのは「この粗い格子でエネルギーは厳密に保存し、"
      "ピークは参照以上に立つ」ところまでで、C の水準そのものは格子依存である。"
      % (100.0 * float((((f24.get("sub") or {}).get("refinement") or {})
                        .get("std_rel_diff")) or 0),
         "PASS" if gw("sub", "G-CONV") else "FAIL"))
    A("* **言えない**: 見た目の最終判断。判定は数値統計のみ（目視なし）。")
    A("")
    return "\n".join(L) + "\n"


def _diag_line(g_all, tag, fs):
    gl = ((g_all.get(tag) or {}).get("G-LINE") or {})
    rows = gl.get("rows") or []
    if not rows:
        return "n/a"
    worst = max(float(r.get("ratio") or 0) for r in rows)
    above = sum(1 for r in rows if not r.get("ok"))
    return ("%s（最悪比 %s、4 窓中 %d 窓が写真超）"
            % ("写真以下" if above == 0 else "写真超あり", fs(worst, 3), above))


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-15.json + REPORT-V4-1-15.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    n = build_numbers()
    write_json(NUMBERS, n)
    text = write_report(n)
    atomic_write_bytes(REPORT, text.encode("utf-8"))
    cp.unit_done("report")
    log("  %s %d bytes ; %s %d bytes"
        % (os.path.basename(NUMBERS), os.path.getsize(NUMBERS),
           os.path.basename(REPORT), os.path.getsize(REPORT)))
    log("  wall %.1fs" % (time.time() - t0))
    return n


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--frame24", choices=STAGES, default=None)
    ap.add_argument("--chunk", nargs=2, metavar=("STAGE", "K"), default=None)
    ap.add_argument("--warp", choices=STAGES, default=None)
    ap.add_argument("--encode", choices=STAGES, default=None)
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        for t in STAGES:
            do_frame24(t)
        for t in STAGES:
            for k in range(NCHUNKS):
                do_chunk(t, k)
        for t in STAGES:
            do_warp(t)
        for t in STAGES:
            do_encode(t)
        do_gates()
        do_report()
        return 0
    did = False
    if a.stage0:
        stage0()
        did = True
    if a.frame24:
        do_frame24(a.frame24)
        did = True
    if a.chunk is not None:
        do_chunk(a.chunk[0], int(a.chunk[1]))
        did = True
    if a.warp:
        do_warp(a.warp)
        did = True
    if a.encode:
        do_encode(a.encode)
        did = True
    if a.gates:
        do_gates()
        did = True
    if a.report:
        do_report()
        did = True
    if not did:
        ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
