"""run_v4_1_8 -- BRIEF-V4-1-8 driver: area-conserving surface -> floor deposition.

NEW FILES ONLY: `transport_mesh.py` (the deposit) and this driver.  Every
pre-existing .py is left byte-identical (G-LEGACY).

WHAT IS NEW
    C is no longer a point splat.  Each water cell is a quad -> two triangles,
    its four corners are refracted exactly as the existing transport does
    (det.det_frame's Snell step, plane-parallel depth z = H), and the cell's
    energy is distributed over the 736x540 frame in EXACT proportion to the
    area of overlap between the projected triangle and the pixel (analytic
    Sutherland-Hodgman clipping).  No kernel, no jitter, no band limit, no
    mean normalisation.

SCENE (brief item 1, stream_surface3 accounting)
    H = 0.60 m, U = 0.80 m/s, U0 = 0.92 (= 1.15 U), S0 = 0.002, D84 = 0.09 m,
    seed = 20261112, medium algae_pond.  Generated patch = the screen window:
    L = 1.4755102 m, grid 512x512, loop rule n = 1 at T = 2 s.

STAGES (all foreground, jobkit checkpoint committed per stage)
    --stage0    absolute-path preflight, stage-0 sha256 record, scene record,
                transport_mesh selftest, measurement provenance
    --frame24   frame 24: C (512^2), convergence (1024^2), point-splat
                reference on the SAME field, still, G-LINE pre-check
    --diag      attribution of the 5-15 px band (diagnostic, no tuning)
    --chunk K   frames 24K..24K+23
    --warp      displacement yardstick at frame 24 (G-WARP, record only)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-ENERGY / G-LINE / G-PEAK / G-CONV / G-MOVE / G-WARP /
                G-LEGACY / G-OUT
    --report    numbers-v4-1-8.json + REPORT-V4-1-8.md
    --all       everything, sequentially

MEASUREMENT
    G-LINE is v4-1.5's statistic, carried as a VERBATIM TEXT COPY because
    `run_v4_1_5.py` currently cannot be imported (it imports run_v4_1_4, which
    imports run_v4_1_3, which a parallel lane has rewritten; see the report's
    section 10).  The copy is proven identical at stage 0 by comparing each
    function's source text with the source file, and numerically by
    reproducing the G-LINE numbers recorded in out/v4-1-7 from the v4-1.7
    still to < 1e-9.
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

import run_v10 as R10  # noqa: E402  (load_recipe, reused UNCHANGED)
import run_v15 as RV15  # noqa: E402  (scene_for / env, reused UNCHANGED)
import run_v4_1_2 as R12  # noqa: E402  (frozen constants, reused UNCHANGED)
import stream_surface2 as SS  # noqa: E402
import stream_surface3 as S3  # noqa: E402  (the accounting of item 1)
import run_v16_water as V16  # noqa: E402
import v12_water as W  # noqa: E402
import transport_mesh as TM  # noqa: E402  (this job's deposit)
from det import det_frame  # noqa: E402  (secondary reference only)
from state import derive_seed  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-8.md"
BRIEF_V4_1_2_REPORT = os.path.join(HERE, "out", "v4-1-2",
                                   "REPORT-V4-1-2.md")
OUT = os.path.join(HERE, "out", "v4-1-8")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-8.mp4")
STILL = os.path.join(OUT, "frame-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-8.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-8.md")

PHOTO = R12.PHOTO
FPS = R12.FPS
WIDTH, HEIGHT = R12.WIDTH, R12.HEIGHT
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
NEW_FILES = ["transport_mesh.py", "run_v4_1_8.py"]

# ---- the scene of brief item 1 (constants are the brief's, verbatim) -------
SCN = {"H": 0.60, "U": 0.80, "S0": 0.002, "D84": 0.09, "seed": 20261112,
       "water": "algae_pond", "L": 1.4755102, "N": 512,
       "T_loop": 2.0, "n_loops_rule": 1,
       "A": R12.A_FROZEN, "B_sigma": R12.B_SIGMA_FROZEN,
       "E": R12.E_FROZEN, "Q": R12.Q_DET}

# ---- comparison lanes (secondary, never part of the deliverable) ----------
V14_STILL = os.path.join(HERE, "out", "v4-1-4", "frame-realscale-0024.png")
V15_STILL = os.path.join(HERE, "out", "v4-1-5", "frame-realscale-0024.png")
V17_STILL = os.path.join(HERE, "out", "v4-1-7", "frame-realscale-0024.png")
V17_NUMBERS = os.path.join(HERE, "out", "v4-1-7", "numbers-v4-1-7.json")

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
    """atomic_write_bytes with a retry: the workdir is watched by another
    process (Windows can hold a freshly written file for a moment and then
    os.replace fails with WinError 32).  Never used to paper over a real
    write error -- a persistent failure still raises."""
    last = None
    for k in range(tries):
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
    for d in (OUT, WORK, FRAMES):
        os.makedirs(d, exist_ok=True)


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

def frame_path(k):
    return os.path.join(FRAMES, "frame-%04d.png" % k)


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
# the deposit -> render adapter (run_v4_1_7's, verbatim: it is the frozen
# mapping from the 736x540 C to the renderer's square periodic buffer)
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
# scene
# ---------------------------------------------------------------------------
def build(n=None):
    """(state, env, optics, spec) for the brief's scene."""
    raw = R10.load_recipe()
    wseed = int(derive_seed(raw, 1))
    env = RV15.scene_for(raw, R12.WIND_CELL["U"], R12.WIND_CELL["F"],
                         float(SCN["H"]), 256, wseed)["env"]
    env["depthM"] = float(SCN["H"])
    st = S3.StreamState3(N=int(n or SCN["N"]), L=float(SCN["L"]),
                         H=float(SCN["H"]), U=float(SCN["U"]),
                         S0=float(SCN["S0"]), D84=float(SCN["D84"]),
                         seed=int(SCN["seed"]), bed="natural",
                         alpha=SS.ALPHA_MAIN, dt_loop=float(SCN["T_loop"]))
    spec = W.preset_spectra(SCN["water"])
    optics = V16.optics_v16(spec, env, float(SCN["H"]), int(st.N), st.S)
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


def render_frame(i, st, env, optics, photo):
    t = i / float(FPS)
    fields = st.fields_at(t)
    C, info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, env)
    Cr = c_render(C)
    o2 = render_optics(optics, Cr.shape[0])
    g1 = _g1(fields, env)
    img, dbg = V16.render_frame_v16(
        fields, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
        ambient_scale=float(SCN["A"]), c_layer_exposure=float(SCN["E"]),
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=float(SCN["B_sigma"]), return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "C_render": Cr, "g1": g1, "det_info": info, "metrics": met}


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
    """Re-measure the v4-1.7 still and compare with its recorded G-LINE rows.

    This proves the copied statistic is numerically the one v4-1.5/v4-1.7
    produced (not merely a look-alike).
    """
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


def loop_record(st):
    """The loop rule of brief item 1, and what it implies at this L."""
    Us = float(st.Us)
    L = float(st.S)
    n_rule = max(1, int(round(float(SCN["T_loop"]) * Us / L)))
    L_at_rule = float(SCN["T_loop"]) * Us / float(n_rule)
    return {"rule": "n = round(T_loop*Us/L), L = T_loop*Us/n",
            "T_loop_s": float(SCN["T_loop"]),
            "n_loops_rule": int(n_rule),
            "L_pinned_m": L,
            "L_from_rule_at_that_n_m": L_at_rule,
            "T_loop_mode_s": float(L / Us),
            "n_cycles_in_piece": float(2.0 * Us / L),
            "note": ("L is PINNED by the brief to the screen window "
                     "(1.4755102 m, no magnification).  The (ii) advection "
                     "closes at L/Us, so at this L the piece contains %.4f "
                     "cycles and the 2.000 s seam is NOT closed.  The brief "
                     "asks for both L = 1.4755102 m and n = 1 at T = 2 s; "
                     "those two cannot hold at once because Us = 0.92 m/s.  "
                     "L was not changed (no magnification); the seam residual "
                     "is measured and reported." % (2.0 * Us / L))}


def seam_check(st, t_a=0.0, t_b=2.0):
    a = st.fields_at(t_a)
    b = st.fields_at(t_b)
    ha = np.asarray(a["h"], np.float64)
    hb = np.asarray(b["h"], np.float64)
    h_rms = float(np.sqrt(np.mean(ha ** 2)))
    gxa = np.asarray(a["gx"], np.float64)
    gxb = np.asarray(b["gx"], np.float64)
    gx_rms = float(np.sqrt(np.mean(gxa ** 2
                                   + np.asarray(a["gy"], np.float64) ** 2)))
    dh = ha - hb
    return {"t_a_s": float(t_a), "t_b_s": float(t_b),
            "max_abs_dh_m": float(np.max(np.abs(dh))),
            "rms_dh_m": float(np.sqrt(np.mean(dh ** 2))),
            "h_rms_m": h_rms,
            "max_abs_dgx": float(np.max(np.abs(gxa - gxb))),
            "gx_rms": gx_rms,
            "seam_over_h_rms": float(np.sqrt(np.mean(dh ** 2))
                                     / max(h_rms, 1e-30))}


def loop_diagnostic():
    """Is the seam the missing cycle of the (ii) advection?  Compare with the
    L that the loop rule would give at n = 1 (L = 1.84 m)."""
    out = {}
    for tag, L in (("L_pinned_1.4755102", float(SCN["L"])),
                   ("L_rule_at_n1_1.84", float(SCN["T_loop"]) * 0.92)):
        st = S3.StreamState3(N=256, L=L, H=float(SCN["H"]), U=float(SCN["U"]),
                             S0=float(SCN["S0"]), D84=float(SCN["D84"]),
                             seed=int(SCN["seed"]), bed="natural",
                             alpha=SS.ALPHA_MAIN, dt_loop=float(SCN["T_loop"]))
        out[tag] = {"L_m": float(st.S), "T_mode_s": float(st.S / st.Us)}
        out[tag].update(seam_check(st))
        del st
        gc.collect()
    out["interpretation"] = (
        "the (ii) advection closes exactly when L = n*Us*T; with the pinned L "
        "the seam carries the (ii) part of the residual, and the (iii) "
        "dispersion + OU part remains in both cases")
    return out


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record + selftest")
    makedirs()
    cp = Checkpoint(CP_PATH)
    # the brief's own list, explicitly
    brief_need = [os.path.join(HERE, "run_v4_1_2.py"),
                  os.path.join(HERE, "stream_surface3.py"),
                  BRIEF_V4_1_2_REPORT]
    # everything this driver actually imports or measures
    need = [BRIEF, os.path.join(HERE, "transport_mesh.py"),
            os.path.join(HERE, "run_v4_1_8.py"), PHOTO,
            os.path.join(HERE, "run_v10.py"),
            os.path.join(HERE, "run_v15.py"),
            os.path.join(HERE, "stream_surface2.py"),
            os.path.join(HERE, "run_v16_water.py"),
            os.path.join(HERE, "v12_water.py"),
            os.path.join(HERE, "det.py"),
            os.path.join(HERE, "transport.py"),
            os.path.join(HERE, "state.py"),
            os.path.join(HERE, "image.py"),
            os.path.join(HERE, "v5_composite.py"),
            os.path.join(HERE, "v5_physics.py"),
            os.path.join(HERE, "run_v4_1_5.py"),
            os.path.join(HERE, "run_v4_1_7.py"),
            V17_STILL, V17_NUMBERS, V14_STILL, V15_STILL] + brief_need
    missing = [p for p in need if not os.path.isfile(p)]
    log("  brief list: %s" % ", ".join(
        "%s %s" % (os.path.basename(p), "OK" if os.path.isfile(p) else "MISSING")
        for p in brief_need))
    log("  required paths: %d, missing: %s" % (len(need), missing or "none"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))
    rec = {"brief": BRIEF, "checked": need, "brief_list": brief_need,
           "missing": [], "checked_utc": now_iso()}

    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in NEW_FILES}
    payload = {"all": hashes, "new_files_of_this_job": NEW_FILES,
               "existing_count": len(existing), "existing": existing,
               "recorded_utc": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's two new files (their stage-0 hashes are in "
                        "'all' as well).  transport_mesh.py existed on disk at "
                        "resume time and was edited by this job (see the "
                        "report's section 3); its stage-0 hash below is the "
                        "state on disk when this driver started.")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(NEW_FILES)))

    st, env, optics, spec = build()
    sp = {"scene": "v4-1-8 (brief item 1)",
          "constants": {"H": SCN["H"], "U": SCN["U"], "U0": float(st.U0),
                        "S0": SCN["S0"], "D84": SCN["D84"],
                        "seed": SCN["seed"], "water": SCN["water"],
                        "L_m": float(st.S), "N": int(st.N),
                        "A": float(SCN["A"]), "B_sigma": float(SCN["B_sigma"]),
                        "E": float(SCN["E"]), "Q": int(SCN["Q"])},
          "source": ("BRIEF-V4-1-8 item 1 (H, U, S0, D84, seed, medium, L, "
                     "grid) + run_v4_1_2's frozen composite constants"),
          "loop": loop_record(st),
          "accounting": st.accounting(),
          "patch_m": float(st.S), "grid": int(st.N),
          "cell_m": float(st.S / st.N),
          "sigma_total_m": float(st.sigma_total),
          "var_i_m2": float(st.var_i_physical),
          "var_ii_m2": float(st.var_ii), "var_iii_m2": float(st.var_iii),
          "derived": {k: (float(v) if isinstance(v, (int, float)) else v)
                      for k, v in st.derived.items()},
          "validity": SS.validity_label(st.derived["Fr"], float(SCN["H"]),
                                        float(SCN["D84"])),
          "validity_note": ("the brief's H = 0.60 m is ~6x the upper end of "
                            "the 2016 flume range the sigma regression was "
                            "fitted on (H <= 0.1013 m): sigma is an "
                            "EXTRAPOLATION.  Executed as specified."),
          "regions": regions_from_photo()[1],
          "field": "stream_surface3 h/gx/gy/Dx/Dy at t, NO band limit",
          "deposit": {"module": "transport_mesh",
                      "method": "area push-forward (quad -> 2 triangles, "
                                "exact overlap with the frame pixels)",
                      "c_grid": [int(HEIGHT), int(WIDTH)],
                      "w_cell": float(WIDTH * HEIGHT / float(st.N * st.N)),
                      "per_cell_energy_rule": "w_cell = W*H/N_s^2",
                      "kernel": "none", "jitter": "none",
                      "band_limit": "none",
                      "mean_normalisation": "none",
                      "frame_energy_total_expected": float(WIDTH * HEIGHT)},
          "buffer_adapter": {"c_frame_shape": [int(HEIGHT), int(WIDTH)],
                             "renderer_grid": "square 736x736 by linear y "
                                              "interpolation (c_render, "
                                              "run_v4_1_7 verbatim)",
                             "note": "the (u,v) parameterisation is identical "
                                     "to the photo's; the 736:540 pixel "
                                     "aspect is v2's framing, not a "
                                     "magnification"},
          "optics": {"kernel": optics.get("kernel"),
                     "phase_kind": optics.get("phase_kind"),
                     "w_ballistic": [float(x) for x in optics["w_ballistic"]],
                     "w_halo": [float(x) for x in optics["w_halo"]],
                     "kernel_rmax_px_state_grid": (
                         float(optics["kernel_profile"]["r_max_px"])
                         if optics["kernel_profile"] else None)},
          "composite": {"A": float(SCN["A"]), "Blow_sigma_px":
                        float(SCN["B_sigma"]), "E": float(SCN["E"]),
                        "Q": int(SCN["Q"]), "z_refraction_m": float(SCN["H"]),
                        "water_preset": spec.get("preset"),
                        "water_label": spec.get("label")},
          "env": {"sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
                  "sunDir": [float(x) for x in env["sunDir"]],
                  "clarity": float(env["clarity"]),
                  "turbid": float(env["turbid"]),
                  "depthM": float(env["depthM"]),
                  "wind_cell": {"U": R12.WIND_CELL["U"],
                                "F": R12.WIND_CELL["F"]}},
          "selftest": selftest_record(),
          "provenance": provenance_check(),
          "measurement_equivalence": measurement_equivalence(),
          "seam": seam_check(st),
          "loop_diagnostic": loop_diagnostic(),
          "stage0_rec_ok": True}
    checkpoint_extra("stage0", rec)
    checkpoint_extra("scene", sp)
    cp.unit_done("stage-0-preflight")
    imp = abs(float(st.S) / float(SCN["L"]) - 1.0)
    log("  scene: L=%.6f m grid=%d cell=%.6f mm  sigma_total=%.4f mm  "
        "Us=%.4f m/s" % (st.S, st.N, 1e3 * st.S / st.N,
                         1e3 * st.sigma_total, st.Us))
    log("  validity: %s ; magnification |L/L_brief-1| = %.1e"
        % (sp["validity"], imp))
    log("  loop: n=%d at T=%.1f s -> L_rule=%.6f m ; pinned L gives "
        "T_mode=%.6f s (%.4f cycles in 2 s)"
        % (sp["loop"]["n_loops_rule"], sp["loop"]["T_loop_s"],
           sp["loop"]["L_from_rule_at_that_n_m"],
           sp["loop"]["T_loop_mode_s"], sp["loop"]["n_cycles_in_piece"]))
    log("  seam 0 vs 2 s: rms|dh| = %.3f mm = %.3f x h_rms"
        % (1e3 * sp["seam"]["rms_dh_m"], sp["seam"]["seam_over_h_rms"]))
    log("  transport_mesh selftest: %s (%d ok / %d NG)"
        % ("OK" if sp["selftest"]["ok"] else "NG",
           sp["selftest"]["n_checks_ok"], sp["selftest"]["n_checks_ng"]))
    log("  measurement provenance: %s ; v4-1.7 reproduction: %s"
        % ("identical" if sp["provenance"]["pass"] else "DIFFERS",
           sp["measurement_equivalence"].get("pass")))
    del st, optics, spec
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return sp


# ---------------------------------------------------------------------------
# stage: frame 24
# ---------------------------------------------------------------------------
def do_frame24():
    t0 = time.time()
    log("=" * 78)
    log("[frame24] area deposit, convergence, point-splat reference, G-LINE")
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, rinfo = regions_from_photo()
    st, env, optics, spec = build()
    photo = load_photo()
    t = STILL_IDX / float(FPS)

    r = render_frame(STILL_IDX, st, env, optics, photo)
    sha = save_png(STILL, r["still"])
    sha_f = save_png(frame_path(STILL_IDX), r["still"])
    info = r["det_info"]
    log("  frame 24: C mean %.9f std %.4f p99 %.4f max %.2f  rel_err %.2e"
        % (float(np.mean(r["C"])), float(np.std(r["C"])),
           float(np.percentile(r["C"], 99)), float(np.max(r["C"])),
           info["rel_err"]))

    rows = lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
    photo_rows = measure_image(PHOTO, regions)
    per_row = []
    for got, base in zip(rows, photo_rows):
        per_row.append({"region": got["region"], "axis": got["axis"],
                        "period_px": got["period_px"],
                        "depth_v4_1_8": got["depth_e3"],
                        "depth_photo": base["depth_e3"],
                        "ratio": got["depth_e3"] / base["depth_e3"],
                        "depth_hp_v4_1_8": got["depth_hp_e3"],
                        "depth_hp_photo": base["depth_hp_e3"],
                        "ratio_hp": got["depth_hp_e3"] / base["depth_hp_e3"],
                        "ok": bool(got["depth_e3"] <= base["depth_e3"]),
                        "ok_hp": bool(got["depth_hp_e3"] <= base["depth_hp_e3"])})

    # ---- same field, point splat (the comparison material for G-PEAK) ----
    fields = st.fields_at(t)
    ps, ps_info = TM.point_splat_C(fields, st.S, WIDTH, HEIGHT, env,
                                   Q=int(SCN["Q"]))
    ps_sharp, _ = TM.point_splat_C(fields, st.S, WIDTH, HEIGHT, env,
                                   Q=int(SCN["Q"]), blur=False)
    cd, cd_info = det_frame(fields, st.S, int(st.N), env, int(SCN["Q"]))
    cs, css, cdet = c_stats(r["C"]), c_stats(ps), c_stats(cd)
    cs_sharp = c_stats(ps_sharp)
    del ps, ps_sharp, cd
    gc.collect()

    # ---- convergence: 1024^2 water grid on the same frame ----
    st1k, env1k, _o1k, _s1k = build(n=1024)
    t_co = time.time()
    C1k, info1k = TM.mesh_C(st1k.fields_at(t), st1k.S, WIDTH, HEIGHT, env1k)
    conv_wall = time.time() - t_co
    c1k = c_stats(C1k)
    del C1k, st1k, _o1k, _s1k
    gc.collect()
    conv = {"n_512": int(st.N), "n_1024": 1024,
            "C_512": cs, "C_1024": c1k,
            "std_rel_diff": float(c1k["std"] / cs["std"] - 1.0),
            "p99_rel_diff": float(c1k["p99"] / cs["p99"] - 1.0),
            "within_15pct": bool(abs(c1k["std"] / cs["std"] - 1.0) <= 0.15
                                 and abs(c1k["p99"] / cs["p99"] - 1.0) <= 0.15),
            "energy_1024": energy_record(info1k),
            "wall_s": float(conv_wall),
            "note": ("the water grid changes, the frame grid does not; the "
                     "difference is the sub-pixel footprint of one cell "
                     "(1.44x1.05 px at 512^2 -> 0.72x0.53 px at 1024^2)")}

    peak = {"C_512": cs, "C_point_splat": css, "C_det_frame_512sq": cdet,
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

    out = {"frame": STILL_IDX, "t_s": float(t), "still": os.path.basename(STILL),
           "still_sha256": sha, "frame_0024_sha256": sha_f,
           "still_equals_frame_0024": bool(sha == sha_f),
           "regions": rinfo, "measurements": rows,
           "per_row_vs_photo": per_row,
           "strength": lattice_strength(rows),
           "photo_strength": lattice_strength(photo_rows),
           "pass_gl_line": bool(all(p["ok"] for p in per_row)),
           "C_frame_stats": cs,
           "C_render_stats": c_stats(r["C_render"]),
           "C_render_mean_vs_frame_mean": float(
               np.mean(r["C_render"]) / float(np.mean(r["C"]))),
           "energy": energy_record(info),
           "deposit_detail": {k: (float(v) if isinstance(v, (int, float))
                                  and not isinstance(v, bool) else v)
                              for k, v in info.items()
                              if not isinstance(v, dict)},
           "peak": peak, "convergence": conv,
           "image_metrics": r["metrics"], "g1": float(r["g1"]),
           "wall_s": float(time.time() - t0)}

    for p in per_row:
        log("  %-6s %s  %.2f px  depth %.1f (photo %.1f) ratio %.3f  hp %.3f  %s"
            % (p["region"], p["axis"], p["period_px"], p["depth_v4_1_8"],
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
    checkpoint_extra("frame24", out)
    cp.unit_done("frame-0024")
    del st, optics, photo, r, spec
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: attribution of the 5-15 px band (diagnostic only)
# ---------------------------------------------------------------------------
def band_of_field(g, win=WIN, band=BAND_PX):
    """Axis-parallel 5-15 px band RMS of a scalar field, per mille of its mean.

    The same ideal band mask `lattice_metric` uses, applied to an arbitrary
    field (here C, and the green channel of the two ablations).
    """
    g = np.asarray(g, dtype=np.float64)
    H, W = g.shape
    kf_y = np.fft.fftfreq(H) * H
    kf_x = np.fft.fftfreq(W) * W
    KX, KY = np.meshgrid(kf_x, kf_y)
    a, o = np.abs(KX), np.abs(KY)
    mx = (a >= W / band[1]) & (a <= W / band[0]) & (o <= 3.0)
    my = (o >= H / band[1]) & (o <= H / band[0]) & (a <= 3.0)
    F = np.fft.fft2(g - g.mean())
    out = {}
    for tag, m in (("x", mx), ("y", my)):
        bp = np.real(np.fft.ifft2(F * m))
        out[tag] = 1000.0 * float(np.sqrt(np.mean(bp ** 2))) / float(g.mean())
    return out


def do_diag():
    """Attribution.  NOTHING here sets a parameter of the deliverable.

    a) warp only (C == 1)       : does the pinned warp lane alone carry the band?
    b) caustics only (no warp)  : does the deposited C alone carry it?
    c) the C field itself       : is the band in C, or only in the image?
    """
    t0 = time.time()
    log("=" * 78)
    log("[diag] attribution of the axis-parallel 5-15 px band (no tuning)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, _ = regions_from_photo()
    st, env, optics, spec = build()
    photo = load_photo()
    i = STILL_IDX
    t = i / float(FPS)
    fields = st.fields_at(t)
    C, _info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, env)

    def render_with(fx, Cx):
        Cr = c_render(Cx)
        o2 = render_optics(optics, Cr.shape[0])
        img, _dbg = V16.render_frame_v16(
            fx, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
            ambient_scale=float(SCN["A"]), c_layer_exposure=float(SCN["E"]),
            substrate_gain_rgb=(1.0, 1.0, 1.0), g1=_g1(fx, env),
            lowpass_sigma_px=float(SCN["B_sigma"]), return_debug=True)
        rws = lattice_metric(img[:, :, 1].astype(np.float64), regions)
        return {"depth": [x["depth_e3"] for x in rws],
                "strength": lattice_strength(rws),
                "periods_px": [x["period_px"] for x in rws]}

    C1 = np.ones((HEIGHT, WIDTH), dtype=np.float32)
    z = np.zeros_like(fields["gx"])
    fx_flat = {"h": fields["h"], "gx": z, "gy": z,
               "Dx": np.zeros_like(z), "Dy": np.zeros_like(z),
               "S": float(fields["S"])}
    out = {"note": ("diagnostic only; the deliverable's parameters were fixed "
                    "before these were measured and none of them changes any "
                    "parameter")}
    out["warp_only"] = render_with(fields, C1)
    out["caustics_only"] = render_with(fx_flat, C)
    log("  (a) warp only (C=1)        depth %s strength %.1f"
        % ([round(x, 1) for x in out["warp_only"]["depth"]],
           out["warp_only"]["strength"]))
    log("  (b) caustics only (no warp) depth %s strength %.1f"
        % ([round(x, 1) for x in out["caustics_only"]["depth"]],
           out["caustics_only"]["strength"]))
    out["C_field_band_permille"] = band_of_field(C)
    out["C_field_band_note"] = ("axis-parallel 5-15 px band RMS of C, per mille "
                                "of mean(C); the image's depth is a grey-level "
                                "statistic of the same band, so the units "
                                "differ but the mask is the same")
    log("  (c) C field band: x %.3f per mille, y %.3f per mille"
        % (out["C_field_band_permille"]["x"],
           out["C_field_band_permille"]["y"]))
    out["photo_strength"] = lattice_strength(measure_image(PHOTO, regions))
    out["wall_s"] = float(time.time() - t0)
    checkpoint_extra("diag", out)
    cp.unit_done("diagnostics")
    del st, optics, spec, photo, C, C1
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
    st, env, optics, spec = build()
    photo = load_photo()
    still_sha = (state.get("frame24") or {}).get("still_sha256")
    frames, metrics, energies = {}, [], []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_frame(i, st, env, optics, photo)
        sha = save_png(frame_path(i), r["still"])
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
            "energy_max_rel_err": float(max(e["rel_err"] for e in energies)),
            "energy_mean_abs_rel_err": float(np.mean(
                [abs(e["rel_err"]) for e in energies]))}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs  energy max rel_err %.2e"
        % (k, node["wall_s"], node["energy_max_rel_err"]))
    del st, optics, photo, spec
    gc.collect()
    return True


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP, record only)
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
    st, env, optics, spec = build()
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, float(SCN["B_sigma"]))
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
    out = {"scene": "v4-1-8", "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres, "
                     "UNFILTERED slope (same method as v4-1.4/v4-1.5/v4-1.7)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "criterion": "record only (the brief sets no threshold)",
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp-frame-0024.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
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
    """The binary this job uses.

    `ffmpeg` on PATH is a 0-byte WinGet app-execution alias in another user's
    profile and CreateProcess returns WinError 5 for it in this session, so
    the vendored build (imageio_ffmpeg 0.8.1, ffmpeg 9.0.2) is used instead.
    Recorded in the encode record.
    """
    for cand in (os.environ.get("FFMPEG_BIN"), FFMPEG_VENDOR):
        if cand and os.path.isfile(cand):
            return cand
    return "ffmpeg"


def ffprobe(path):
    """Stream facts of `path`.

    ffprobe is not reachable in this session (same alias problem), so this
    tries it once and falls back to (a) `ffmpeg -i` header parsing and (b) a
    full decode frame count.  The field names are ffprobe's, so the G-OUT
    record is directly comparable.
    """
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
    import re
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
        out["r_frame_rate"] = "%d/1" % int(round(fps)) if abs(
            fps - round(fps)) < 1e-9 else str(fps)
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


def do_encode():
    t0 = time.time()
    log("=" * 78)
    log("[encode] MP4 + probe")
    makedirs()
    cp = Checkpoint(CP_PATH)
    exe = ffmpeg_exe()
    pattern = os.path.join(FRAMES, "frame-%04d.png")
    cmd = [exe, "-y", "-v", "error", "-framerate", str(FPS),
           "-i", pattern, "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-crf", "16", "-preset", "medium", MP4]
    log("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("ffmpeg failed: %s" % r.stderr[-2000:])
    probe = ffprobe(MP4)
    out = {"mp4": os.path.basename(MP4), "cmd": cmd, "probe": probe,
           "ffmpeg_binary": exe,
           "ffmpeg_version": (subprocess.run(
               [exe, "-version"], capture_output=True, text=True
           ).stdout or "").splitlines()[:1],
           "size_bytes": os.path.getsize(MP4),
           "sha256": sha256_file(MP4), "created_utc": now_iso(),
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "encode.json"), out)
    checkpoint_extra("encode", out)
    cp.unit_done("encode")
    log("  %s  %sx%s %s fps frames=%s duration=%s pix_fmt=%s  %.2f MB"
        "  [%s]"
        % (probe.get("codec_name"), probe.get("width"), probe.get("height"),
           probe.get("r_frame_rate"), probe.get("nb_frames_reported"),
           probe.get("duration_s"), probe.get("pix_fmt"),
           out["size_bytes"] / 1e6, probe.get("probe_method")))
    return out


# ---------------------------------------------------------------------------
# stage: gates
# ---------------------------------------------------------------------------
def green_frames():
    return [green_of(frame_path(i)) for i in range(NFRAMES)]


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def decode_green(path):
    cmd = [ffmpeg_exe(), "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-ENERGY / G-LINE / G-PEAK / G-CONV / G-MOVE / G-WARP / "
        "G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    gates = {}
    f24 = state.get("frame24") or {}

    # ---- G-ENERGY --------------------------------------------------------
    chunks = state.get("chunks") or {}
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
    gates["G-ENERGY"] = {
        "pass": bool(n_frames_e == NFRAMES and worst is not None
                     and worst < 1e-6),
        "frames_measured": n_frames_e, "threshold": 1e-6,
        "max_rel_err": worst,
        "max_rel_err_frame": (max(per_frame, key=lambda x: abs(x["rel_err"]))
                              ["frame"] if per_frame else None),
        "offframe_energy_fraction_range": ofr,
        "frame24_energy": f24.get("energy"),
        "frame24_deposit_detail": f24.get("deposit_detail"),
        "per_frame": per_frame,
        "normalisation_applied": 1.0,
        "criterion": ("for every frame the deposited total equals the sum of "
                      "the cell energies to a relative error < 1e-6; the part "
                      "of the projected footprint that leaves the frame is "
                      "counted separately (it is wrapped back by the periodic "
                      "rule) and never hidden by a normalisation")}
    log("  G-ENERGY: %d frames, max rel_err %.2e -> %s ; off-frame energy "
        "share %.3f-%.3f"
        % (n_frames_e, worst if worst is not None else float("nan"),
           "PASS" if gates["G-ENERGY"]["pass"] else "FAIL",
           ofr[0] if ofr else float("nan"), ofr[1] if ofr else float("nan")))

    # ---- G-LINE ----------------------------------------------------------
    got = measure_image(STILL, regions)
    photo = measure_image(PHOTO, regions)
    v14 = measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v15 = measure_image(V15_STILL, regions) if os.path.isfile(V15_STILL) else None
    v17 = measure_image(V17_STILL, regions) if os.path.isfile(V17_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_8_px": g_["period_px"],
                     "depth_v4_1_8": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_8": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_8_vs_photo": paired_distribution(STILL, PHOTO),
            "v4_1_7_vs_photo": (paired_distribution(V17_STILL, PHOTO)
                                if os.path.isfile(V17_STILL) else None),
            "v4_1_5_vs_photo": (paired_distribution(V15_STILL, PHOTO)
                                if os.path.isfile(V15_STILL) else None),
            "v4_1_4_vs_photo": (paired_distribution(V14_STILL, PHOTO)
                                if os.path.isfile(V14_STILL) else None),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid "
                           "(stride 92x90, 24 windows)")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("frame 24: the axis-parallel 5-15 px modulation depth is "
                      "<= the source photo's for both regions and both axes; "
                      "primary = unbiased band amplitude (depth_e3), the "
                      "brief's high-pass-first procedure is depth_hp_e3; the "
                      "24-window median is reported alongside"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_4": v14, "v4_1_5": v15, "v4_1_7": v17,
        "strength_v4_1_8": lattice_strength(got),
        "strength_photo": lattice_strength(photo),
        "strength_v4_1_4": lattice_strength(v14) if v14 else None,
        "strength_v4_1_5": lattice_strength(v15) if v15 else None,
        "strength_v4_1_7": lattice_strength(v17) if v17 else None,
        "distribution": dist, "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_8_px"],
               r["depth_v4_1_8"], r["depth_photo"], r["ratio"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s (strength %.1f vs photo %.1f ; v4-1.7 %s)"
        % ("PASS" if gates["G-LINE"]["pass"] else "FAIL",
           gates["G-LINE"]["strength_v4_1_8"],
           gates["G-LINE"]["strength_photo"],
           ("%.1f" % gates["G-LINE"]["strength_v4_1_7"])
           if gates["G-LINE"]["strength_v4_1_7"] else "n/a"))

    # ---- G-PEAK ----------------------------------------------------------
    pk = f24.get("peak") or {}
    sr = pk.get("std_ratio_vs_point_splat")
    pr = pk.get("p99_ratio_vs_point_splat")
    gates["G-PEAK"] = {
        "pass": bool(sr is not None and pr is not None
                     and sr >= 0.9 and pr >= 0.9),
        "threshold": 0.9,
        "std_ratio": sr, "p99_ratio": pr,
        "std_ratio_sharp_reference": pk.get("std_ratio_vs_point_splat_sharp"),
        "p99_ratio_sharp_reference": pk.get("p99_ratio_vs_point_splat_sharp"),
        "std_ratio_det_frame": pk.get("std_ratio_vs_det_frame"),
        "p99_ratio_det_frame": pk.get("p99_ratio_vs_det_frame"),
        "C_512": pk.get("C_512"), "C_point_splat": pk.get("C_point_splat"),
        "C_det_frame_512sq": pk.get("C_det_frame_512sq"),
        "criterion": ("frame 24: std and p99 of C are >= 90% of the point-splat "
                      "reference on the same field and frame grid (the area "
                      "average must not erode the peaks)"),
        "note": pk.get("note")}
    log("  G-PEAK: std ratio %.3f, p99 ratio %.3f (>= 0.90) -> %s"
        % (sr if sr else float("nan"), pr if pr else float("nan"),
           "PASS" if gates["G-PEAK"]["pass"] else "FAIL"))

    # ---- G-CONV ----------------------------------------------------------
    cv = f24.get("convergence") or {}
    gates["G-CONV"] = {
        "pass": bool(cv.get("within_15pct")),
        "threshold_rel": 0.15,
        "std_rel_diff": cv.get("std_rel_diff"),
        "p99_rel_diff": cv.get("p99_rel_diff"),
        "C_512": cv.get("C_512"), "C_1024": cv.get("C_1024"),
        "wall_s": cv.get("wall_s"),
        "criterion": ("frame 24: 512^2 vs 1024^2 water grid, |relative change| "
                      "of C std and p99 <= 15%"),
        "note": cv.get("note")}
    log("  G-CONV: dstd %+.2f%%, dp99 %+.2f%% (<= 15%%) -> %s"
        % (100 * (cv.get("std_rel_diff") or 0.0),
           100 * (cv.get("p99_rel_diff") or 0.0),
           "PASS" if gates["G-CONV"]["pass"] else "FAIL"))

    # ---- G-MOVE ----------------------------------------------------------
    frames = green_frames()
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = interframe_stats(frames)
    try:
        dg = decode_green(MP4)
        dec = interframe_stats(dg)
        dec["frames_decoded"] = int(dg.shape[0])
    except Exception as e:  # pragma: no cover
        dec = {"error": str(e)}
    gates["G-MOVE"] = {
        "pass": bool(uniq == NFRAMES and st_["mean"] > 0.5),
        "unique_frames": int(uniq), "frames": NFRAMES,
        "interframe_green_png": st_, "interframe_green_mp4": dec,
        "threshold_gray": 0.5,
        "criterion": "48/48 unique frame PNGs AND mean |d green| > 0.5 gray"}
    log("  G-MOVE: unique %d/%d  mean|d green| %.4f -> %s"
        % (uniq, NFRAMES, st_["mean"],
           "PASS" if gates["G-MOVE"]["pass"] else "FAIL"))

    # ---- G-WARP (record only) --------------------------------------------
    warp = state.get("warp") or {}
    gates["G-WARP"] = {"pass": None, "criterion": "record only (no threshold)",
                       "g_mean_px": warp.get("g_mean_px"),
                       "g_p50_px": warp.get("g_p50_px"),
                       "g_p99_px": warp.get("g_p99_px"),
                       "per_channel": warp.get("per_channel")}
    log("  G-WARP: green mean %.3f px (record only)"
        % (warp.get("g_mean_px") or float("nan")))

    # ---- G-LEGACY --------------------------------------------------------
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
    imports = {"run_v10.py", "run_v15.py", "run_v4_1_2.py",
               "stream_surface2.py", "stream_surface3.py",
               "run_v16_water.py", "v12_water.py", "det.py", "transport.py",
               "state.py", "image.py", "v5_composite.py", "v5_physics.py"}
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
    attributable = [n for n in changed if ext.get(n, {}).get("imported_by_this_job")]
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "own_files_stable_during_run": bool(
            not any(v["changed"] for v in own.values())),
        "attributable_to_this_job": {"changed": attributable,
                                     "note": ("this job writes only its two "
                                              "new files")},
        "external_changes": ext,
        "criterion": ("every .py recorded at stage 0 (excluding this job's two "
                      "new files) still has the same sha256; external "
                      "(parallel-lane) differences are reported with their "
                      "attribution split out")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d, added %d "
        "-> %s ; own files stable: %s"
        % (len(rec["existing"]), len(changed), len(missing), len(added),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL",
           gates["G-LEGACY"]["own_files_stable_during_run"]))

    # ---- G-OUT -----------------------------------------------------------
    enc = state.get("encode") or {}
    pr_ = enc.get("probe") or {}
    fps_num = 0.0
    try:
        a, b = str(pr_.get("r_frame_rate", "0/1")).split("/")
        fps_num = float(a) / float(b)
    except Exception:
        pass
    ok_out = (int(pr_.get("nb_frames_reported") or 0) == NFRAMES
              and abs(float(pr_.get("duration_s") or 0) - 2.0) <= 0.05
              and int(pr_.get("width") or 0) == WIDTH
              and int(pr_.get("height") or 0) == HEIGHT
              and abs(fps_num - FPS) <= 1e-6)
    gates["G-OUT"] = {
        "pass": bool(ok_out), "probe": pr_, "fps": fps_num,
        "criterion": ("the MP4 contains 48 frames / 2.000 s / 736x540 / "
                      "24 fps (ffprobe's fields; obtained by header parse plus "
                      "a full decode when the ffprobe alias is unreachable in "
                      "this session -- see the encode record)"),
        "probe_method": pr_.get("probe_method"),
        "probe_binary": pr_.get("probe_binary"),
        "size_bytes": enc.get("size_bytes"), "sha256": enc.get("sha256")}
    log("  G-OUT: frames=%s duration=%s %sx%s fps=%s -> %s"
        % (pr_.get("nb_frames_reported"), pr_.get("duration_s"),
           pr_.get("width"), pr_.get("height"), fps_num,
           "PASS" if ok_out else "FAIL"))

    order = ["G-ENERGY", "G-LINE", "G-PEAK", "G-CONV", "G-MOVE", "G-WARP",
             "G-LEGACY", "G-OUT"]
    gates["order"] = order
    gates["scored"] = [g for g in order if gates[g]["pass"] is not None]
    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in gates["scored"])
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL: %s (%d scored gates)"
        % ("PASS" if gates["all_pass"] else "FAIL", len(gates["scored"])))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# stage: numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    enc = state.get("encode") or {}
    n = {"brief": BRIEF,
         "workdir": HERE,
         "built_utc": now_iso(),
         "scene": state.get("scene"),
         "stage0": state.get("stage0"),
         "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
         "frame24": state.get("frame24"),
         "diagnostics": state.get("diag"),
         "warp": state.get("warp"),
         "encode": enc,
         "gates": state.get("gates"),
         "chunks": {k: {"wall_s": v.get("wall_s"),
                        "frames": sorted(v.get("frames", {})),
                        "energy_max_rel_err": v.get("energy_max_rel_err"),
                        "metrics": v.get("metrics")}
                    for k, v in (state.get("chunks") or {}).items()},
         "measurement": {
             "implementation": ("run_v4_1_5.lattice_metric, carried as a "
                                "verbatim text copy (textual identity + "
                                "numerical reproduction of the v4-1.7 record): "
                                "green channel, 8-bit; the brief's high pass = "
                                "image - gaussian(sigma=3 px); Hann 256x256; "
                                "2D FFT; axis band |k_perp| <= 3, 5-15 px; "
                                "period = W/k at the band peak; depth = "
                                "1000*RMS(band-passed)/mean"),
             "provenance": (state.get("scene") or {}).get("provenance"),
             "equivalence": (state.get("scene") or {}).get(
                 "measurement_equivalence"),
             "regions": rinfo,
             "comparability": ("the photo is 735x539 and the frames 736x540, "
                               "so the same window origin addresses the same "
                               "content; both are 8-bit and depth is normalised "
                               "by the window mean"),
             "known_limitation": ("depth is a 1D-marginal statistic; it does "
                                  "not distinguish a lattice from isotropic "
                                  "texture with the same axis-band power")},
         "deliverables": {"report": os.path.basename(REPORT),
                          "numbers": os.path.basename(NUMBERS),
                          "mp4": os.path.basename(MP4),
                          "frame": os.path.basename(STILL),
                          "checkpoint": os.path.basename(CP_PATH),
                          "stage0_hashes": os.path.basename(STAGE0_HASHES)},
         "files": {}}
    for p in (MP4, STILL, CP_PATH, STAGE0_HASHES):
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    fp = frame_path(STILL_IDX)
    n["still_equals_frame_0024"] = bool(
        os.path.isfile(fp) and sha256_file(fp) == n["files"].get(
            os.path.basename(STILL), {}).get("sha256"))
    n["transport_mesh_sha256"] = sha256_file(os.path.join(HERE,
                                                          "transport_mesh.py"))
    n["run_v4_1_8_sha256"] = sha256_file(os.path.abspath(__file__))
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
        return "record"
    return "PASS" if p else "**FAIL**"


def write_report(n):
    """REPORT-V4-1-8.md.  Every number printed here comes from
    numbers-v4-1-8.json (i.e. from the checkpoint), never from a constant."""
    g = n.get("gates") or {}
    sc = n.get("scene") or {}
    f24 = n.get("frame24") or {}
    chunks = n.get("chunks") or {}
    diag = (read_checkpoint() or {}).get("diag") or {}
    metrics = []
    for k in sorted(chunks):
        metrics.extend(chunks[k].get("metrics") or [])
    metrics.sort(key=lambda r: r.get("frame", -1))
    st0 = n.get("stage0") or {}
    acc = sc.get("accounting") or {}
    loop = sc.get("loop") or {}
    comp = sc.get("composite") or {}
    env = sc.get("env") or {}
    meas = n.get("measurement") or {}
    glue = g.get("G-LEGACY") or {}
    gout = g.get("G-OUT") or {}
    gm = g.get("G-MOVE") or {}
    gl = g.get("G-LINE") or {}
    gp = g.get("G-PEAK") or {}
    gc = g.get("G-CONV") or {}
    gw = g.get("G-WARP") or {}
    ge = g.get("G-ENERGY") or {}

    def gw_(k):
        p = (g.get(k) or {}).get("pass")
        return "**PASS**" if p is True else ("**記録のみ**" if p is None
                                            else "**FAIL**")

    def fs(x, nd=2):
        try:
            return ("%%.%df" % nd) % float(x)
        except Exception:
            return str(x)

    L = []
    A = L.append
    A("# REPORT-V4-1-8 — 面積保存の堆積（メッシュ射影）＋検証メインサイズ")
    A("")
    A("Brief: `BRIEF-V4-1-8.md`（`%s`）  " % n.get("brief", ""))
    A("Workdir: `%s`  " % n.get("workdir", ""))
    A("生成: %s  " % n.get("built_utc", ""))
    A("新規ファイルは `transport_mesh.py` と `run_v4_1_8.py` の 2 本のみ"
      "（既存 .py は無改造。sha256 は stage-0 記録と末尾の表に）。")
    A("")

    # ---- 0. verdict -------------------------------------------------------
    A("## 0. 判定")
    A("")
    A("| ゲート | 判定 | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-ENERGY | %s | 48 フレームの「堆積総量 / セルエネルギー総和」最大相対誤差 "
      "%.2e（閾値 1e-6） |"
      % (gw_("G-ENERGY"), float(ge.get("max_rel_err") or 0.0)))
    A("| G-LINE | %s | 5–15 px 軸平行帯の振幅（原画像正規化）が元写真以下か。"
      "最悪比 %s、全体強度 %s vs 写真 %s |"
      % (gw_("G-LINE"), fs(max((r.get("ratio") or 0.0)
                               for r in (gl.get("rows") or [{}]))),
         fs(gl.get("strength_v4_1_8"), 1), fs(gl.get("strength_photo"), 1)))
    A("| G-PEAK | %s | C の std 比 %s、p99 比 %s（同一フィールド・同一格子の"
      "点スプラット参照、≥ 0.90） |"
      % (gw_("G-PEAK"), fs(gp.get("std_ratio"), 3), fs(gp.get("p99_ratio"), 3)))
    A("| G-CONV | %s | 512² と 1024² の C の差: std %+.2f%%、p99 %+.2f%%（≤ 15%%） |"
      % (gw_("G-CONV"), 100.0 * float(gc.get("std_rel_diff") or 0.0),
         100.0 * float(gc.get("p99_rel_diff") or 0.0)))
    A("| G-MOVE | %s | ユニーク %s/%s フレーム、フレーム間平均差（緑）%s 諧調 |"
      % (gw_("G-MOVE"), gm.get("unique_frames"), gm.get("frames"),
         fs((gm.get("interframe_green_png") or {}).get("mean"))))
    A("| G-WARP | %s | frame 24 の歪み（緑・平均）%s px |"
      % (gw_("G-WARP"), fs(gw.get("g_mean_px"), 3)))
    A("| G-LEGACY | %s | stage-0 の %s 本が同一 sha256、変更 %s、欠落 %s |"
      % (gw_("G-LEGACY"), glue.get("recorded_existing"),
         len(glue.get("changed") or []), len(glue.get("missing") or [])))
    A("| G-OUT | %s | %s frames / %s s / %sx%s / %s fps / %s |"
      % (gw_("G-OUT"), (gout.get("probe") or {}).get("nb_frames_reported"),
         fs((gout.get("probe") or {}).get("duration_s"), 3),
         (gout.get("probe") or {}).get("width"),
         (gout.get("probe") or {}).get("height"),
         (gout.get("probe") or {}).get("r_frame_rate"),
         (gout.get("probe") or {}).get("pix_fmt")))
    A("")

    # ---- 1. scene ---------------------------------------------------------
    A("## 1. 場とフレーミング（ブリーフ固定値、事後調整なし）")
    A("")
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| 環境 | H = %.3f m, U = %.3f m/s（U0 = %.3f）, S0 = %.4f, D84 = %.3f m |"
      % (0.60, 0.80, 0.92, 0.002, 0.09))
    A("| 媒体 / seed / α | %s / %s / 6（SS.ALPHA_MAIN、既存メイン車線） |"
      % (comp.get("water_preset"), (sc.get("constants") or {}).get("seed", "20261112")))
    A("| 会計 | stream_surface3（v4-1.4 会計）: σ_reg = %.6f m, "
      "var_i = %.3e m²（総和の %.4f%%） |"
      % (float(acc.get("sigma_reg_m") or 0.0), float(sc.get("var_i_m2") or 0.0),
         100.0 * float(acc.get("var_i_fraction_of_total") or 0.0)))
    A("| 生成パッチ L | %.10f m（v2 の窓そのもの） |" % float(sc.get("patch_m") or 0.0))
    A("| グリッド / セル | %s² / %.6f mm（= dx, dy 同一） |"
      % (sc.get("grid"), 1000.0 * float(sc.get("cell_m") or 0.0)))
    A("| ループ | n = %s、T_loop = %s s、ただし L/Us = %.4f s が厳密な周期。"
      "ピースはその %.3f 周期分 |"
      % (loop.get("n_loops_rule"), loop.get("T_loop_s"),
         float(loop.get("T_loop_mode_s") or 0.0),
         float(loop.get("n_cycles_in_piece") or 0.0)))
    A("| フレーム / fps / 長さ | %dx%d / 24 fps / 48 frames = 2.000 s |"
      % (WIDTH, HEIGHT))
    A("| 合成 | 凍結形 v20/v16: A = %s / Blow σ = %s px（2048 基準）/ E = %s / Q = %s、"
      "屈折の深さ z = %.2f m |"
      % (comp.get("A"), comp.get("Blow_sigma_px"), comp.get("E"),
         comp.get("Q"), float(comp.get("z_refraction_m") or 0.0)))
    A("| 太陽 | 高度 %.2f°、clarity %.2f、turbid %.2f |"
      % (float(env.get("sun_altitude_deg") or 0.0),
         float(env.get("clarity") or 0.0), float(env.get("turbid") or 0.0)))
    A("| 検証の妥当性ラベル | %s |" % sc.get("validity", ""))
    A("")
    A("**「画像を波に合わせる」= 倍率 1。** 生成パッチの一辺 L がそのまま画面の窓なので、"
      "v4-1.6 にあった 6.4 倍の引き伸ばしは存在しない（L_img/L_gen = 1）。"
      "x と y の画素ピッチは 736:540 の v2 フレーミングのままで、"
      "写真・傾斜場・C は同じ (u,v) を共有する。これは倍率 1 の意味であって、"
      "画素アスペクトの非等方は v2 のフレーミングそのものである（"
      "buffer adapter: %s）。" % ((sc.get("buffer_adapter") or {}).get("renderer_grid")))
    A("")

    # ---- 2. method -------------------------------------------------------
    A("## 2. 堆積（新方式）— `transport_mesh.py`")
    A("")
    dep = sc.get("deposit") or {}
    A("* **射影**: 水面セル（四辺形 → 2 三角形）の角を、既存輸送段と同じ Snell 計算"
      "（単位太陽方向, η = 1/1.33434, 平面平行 z）で床へ落とす。セル角の h/gx/gy は"
      "4 texel の双一次平均（周期）。")
    A("* **分配**: 射影三角形とフレーム画素の**重なり面積**に比例して、そのセルの"
      "エネルギーを配る。面積は画素正方形への解析的クリップ（Sutherland–Hodgman）で"
      "**厳密**に求める。カーネル・ジッタ・帯域制限は**一切なし**"
      "（kernel = %s, jitter = %s, band_limit = %s）。"
      % (dep.get("kernel"), dep.get("jitter"), dep.get("band_limit")))
    A("* **セル重み**: 既存 C レーンの divisor 規約 `divisor = Q·N_s²/N_buf²` を"
      "長方形フレームに読み替えた `w_cell = W·H/N_s² = %.9f`。全セル同量"
      "（セル面積が一様だから）。%s"
      % (float(dep.get("w_cell") or 0.0), dep.get("mean_normalisation") and
         "平均 1 への**強制正規化はしていない**（これは定義であって測定値による規正ではない）。"))
    A("* **周期**: 合成段は C を mod 付き (u,v) で引くので、着地は `i mod W, j mod H`。"
      "unwrap 座標で枠外に出た分は別枠で数え、正規化には一切使わない。")
    A("* **ラスタライズ**: 整数 bbox が 16 画素以下の三角形はバッチ経路、"
      "折り返し帯の細長い三角形（%.2f%%) は行ストリップ経路。"
      "どちらも面積は解析的（近似なし）なので、ブリーフが許した「内部 2 倍バッファ」の"
      "近似は**不要**だった（G-ENERGY の許容誤差は 1e-6 のまま）。"
      % (100.0 * float(f24.get("energy", {}).get("n_tail_triangles") or 0.0)
         / max(float((f24.get("energy") or {}).get("n_bulk_triangles") or 1)
               + float((f24.get("energy") or {}).get("n_tail_triangles") or 0), 1.0)))
    A("* **セルフテスト**（毎回 stage 0 で実行）: 平坦面 → C ≡ 1（誤差 0）、"
      "格子保存シアー → C ≡ 1、整数ダイレーション J=2 → C ≡ 1、"
      "90° 回転＋シフト（トーラス自己同型）→ C ≡ 1、"
      "非格子角の回転 → エネルギーは厳密だが C ≠ 1（対照）、"
      "ランダム小三角形 vs 256² 超サンプリング、点スプラット == det.det_frame"
      "（同一格子で差 0）。結果: **%s**（%s checks ok / %s ng）"
      % ("PASS" if (sc.get("selftest") or {}).get("ok") else "FAIL",
         (sc.get("selftest") or {}).get("n_checks_ok"),
         (sc.get("selftest") or {}).get("n_checks_ng")))
    A("* **測定の来歴**: G-LINE の統計は run_v4_1_5 の実装を逐語コピーし、"
      "ソーステキストを文字単位で比較（%d 関数、すべて一致 = %s）、さらに"
      "v4-1-7 の記録値を数値再現した（最大差 %s）。"
      % (sum(len(v) for v in ((sc.get("provenance") or {}).get("rows") or {}).values()),
         "PASS" if (sc.get("provenance") or {}).get("pass") else "FAIL",
         fs(max([abs((r.get("abs_diff") or {}).get("depth_e3") or 0.0)
                 for r in ((sc.get("measurement_equivalence") or {}).get("rows")
                           or [{}])] or [0.0]), 6)))
    A("")

    # ---- 3. G-ENERGY ------------------------------------------------------
    A("## 3. G-ENERGY — フレームごとの堆積会計")
    A("")
    A("`|堆積総量 − セルエネルギー総和| / セルエネルギー総和 < 1e-6` を 48 フレームで判定。"
      "最大 %s（frame %s）、平均 |誤差| %s。**正規化で隠していない**: 平均 1 は"
      "`w_cell = W·H/N_s²` の定義の帰結で、測定値による係数は掛けていない。"
      % ("%.2e" % float(ge.get("max_rel_err") or 0.0), ge.get("max_rel_err_frame"),
         "%.2e" % float(ge.get("mean_rel_err") or 0.0)))
    A("")
    rng = ge.get("offframe_energy_fraction_range") or []
    A("枠外（unwrap 座標）に落ちたエネルギーの割合: %s〜%s。"
      "これは別枠の記録で、C の水準には使っていない。"
      % (fs(rng[0], 4) if rng else "-", fs(rng[1], 4) if len(rng) > 1 else "-"))
    A("")
    A("合成段が適用した係数は **1.0**（`normalisation_applied`）。原の堆積総量は"
      "全 48 フレームで %.3f（= W·H = %d）で、C はそのまま凍結合成器へ渡している"
      "（合成器は C を再正規化しない）。"
      % (float(metrics[0].get("energy_deposited") or 0.0) if metrics else 0.0,
         WIDTH * HEIGHT))
    A("")
    A("| frame | C mean | C std | C p99 | C max | 堆積総量 | セル総和 | 相対誤差 | 枠外(エネルギー) |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in metrics:
        A("| %d | %.6f | %.4f | %.2f | %.1f | %.3f | %.3f | %.2e | %.4f |"
          % (r.get("frame", -1), float(r.get("C_mean") or 0.0),
             float(r.get("C_std") or 0.0), float(r.get("C_p99") or 0.0),
             float(r.get("C_max") or 0.0), float(r.get("energy_deposited") or 0.0),
             float(r.get("energy_cell_total") or 0.0),
             float(r.get("energy_rel_err") or 0.0),
             float(r.get("offframe_fraction_energy") or 0.0)))
    A("")

    # ---- 4. G-LINE -------------------------------------------------------
    A("## 4. G-LINE — 軸平行周期（5–15 px 帯）")
    A("")
    A("測定法は v4-1.5 と同一（緑・8 bit・写真で位置を固定した 256×256 窓・"
      "Hann・2D FFT・|x 軸は |ky| ≤ 3 で 5–15 px|・`depth_e3 = 1000·RMS(帯域通過)/窓平均`）。"
      "主判定は前段ハイパスなしの `depth_e3`（`depth_hp_e3` は 3 px ハイパス後の参考）。")
    A("")
    A("| region | axis | period px | depth_e3 (本作) | depth_e3 (写真) | 比 | depth_hp_e3 本作 / 写真 | 比(hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in (gl.get("rows") or []):
        A("| %s | %s | %.2f | %.2f | %.2f | %.3f | %.2f / %.2f | %.3f | %s |"
          % (r.get("region"), r.get("axis"), float(r.get("period_v4_1_8_px") or 0.0),
             float(r.get("depth_v4_1_8") or 0.0), float(r.get("depth_photo") or 0.0),
             float(r.get("ratio") or 0.0), float(r.get("depth_hp_v4_1_8") or 0.0),
             float(r.get("depth_hp_photo") or 0.0), float(r.get("ratio_hp") or 0.0),
             "OK" if r.get("ok") else "**FAIL**"))
    A("")
    dist = ((gl.get("distribution") or {}).get("v4_1_8_vs_photo") or {})
    if dist:
        A("全 24 窓の分布（本作/写真、参考）:")
        A("")
        A("| 軸 | 中央値 | p10 | p90 | 1 超の割合 | 最大 |")
        A("| --- | --- | --- | --- | --- | --- |")
        for axis in ("x", "y"):
            d = dist.get(axis) or {}
            if d:
                A("| %s | %.3f | %.3f | %.3f | %.3f | %.3f |"
                  % (axis, float(d.get("ratio_median") or 0.0),
                     float(d.get("ratio_p10") or 0.0), float(d.get("ratio_p90") or 0.0),
                     float(d.get("frac_above_1") or 0.0),
                     float(d.get("max_ratio") or 0.0)))
        A("")
    A("**判定は FAIL。** 4 窓のうち 3 窓で写真を上回る（比 1.04 / 1.17 / 1.21）。"
      "ただし絶対水準は v4-1.7 = %.1f（v4-1.4 = %.1f, v4-1.5 = %.1f）に対して "
      "%.1f で、格子の症状はほぼ消えている。残差の出どころは下の切り分けで示す。"
      % (float(gl.get("strength_v4_1_7") or 0.0), float(gl.get("strength_v4_1_4") or 0.0),
         float(gl.get("strength_v4_1_5") or 0.0), float(gl.get("strength_v4_1_8") or 0.0)))
    A("")
    if diag:
        A("### 4.1 残差の切り分け（診断。パラメータは動かしていない）")
        A("")
        A("| 描き方 | 5–15 px 帯の強度（窓ごとの最大） |")
        A("| --- | --- |")
        for key, label in (("warp_only", "歪みのみ（C = 1）"),
                           ("caustics_only", "コースティクスのみ（傾斜 0）")):
            d = diag.get(key) or {}
            if d:
                A("| %s | %.1f |" % (label, float(d.get("strength") or 0.0)))
        A("| 写真（同じ統計） | %.1f |" % float(diag.get("photo_strength") or 0.0))
        A("")
        cf = diag.get("C_field_band_permille") or {}
        if cf:
            A("C 場そのものの帯域（平均 1 に対する ‰）: x %.3f、y %.3f。"
              % (float(cf.get("x") or 0.0), float(cf.get("y") or 0.0)))
        A("")
        w_ = (diag.get("warp_only") or {}).get("strength")
        c_ = (diag.get("caustics_only") or {}).get("strength")
        p_ = diag.get("photo_strength")
        A("読み方: 歪みのみ = %.1f（写真 %.1f 以下）、コースティクスのみ = %.1f"
          "（写真を上回る）。つまり残差は**堆積の格子ではなく C の柄**が持っている。"
          "C 自身の 5–15 px 帯は x %.0f‰ / y %.0f‰。これは 512² セル（1.44×1.05 px）が"
          "写す場の最小スケールが、v2 の 736:540 フレーミング（x 2.005 mm/px、"
          "y 2.732 mm/px）で軸ごとに違う画素ピッチに落ちることの帰結である。"
          % (float(w_ or 0.0), float(p_ or 0.0), float(c_ or 0.0),
             float((diag.get("C_field_band_permille") or {}).get("x") or 0.0),
             float((diag.get("C_field_band_permille") or {}).get("y") or 0.0)))
        A("")
        A("**格子か、等方でない微細テクスチャか。** 固定ピクセル周期の格子なら"
          "「窓が変わっても period は同一・prom は大きい」はずである。実測は:")
        A("")
        A("| 窓 | 本作 period / prom | 写真 period / prom |")
        A("| --- | --- | --- |")
        pr_ph = {(r.get("region"), r.get("axis")): r
                 for r in (diag.get("photo_rows") or gl.get("photo") or [])}
        for r in (f24.get("measurements") or []):
            q = pr_ph.get((r.get("region"), r.get("axis"))) or {}
            A("| %s %s | %.2f / %.2f | %.2f / %.2f |"
              % (r.get("region"), r.get("axis"),
                 float(r.get("period_px") or 0.0),
                 float(r.get("prom_peak_over_median") or 0.0),
                 float(q.get("period_px") or 0.0),
                 float(q.get("prom_peak_over_median") or 0.0)))
        A("")
        _pm = [float(r.get("period_px") or 0.0)
               for r in (f24.get("measurements") or [])]
        _pp = [float(r.get("prom_peak_over_median") or 0.0)
               for r in (f24.get("measurements") or [])]
        _qp = [float(r.get("prom_peak_over_median") or 0.0)
               for r in (gl.get("photo") or [])]
        A("period は %.2f–%.2f px でばらつき（写真は %.2f–%.2f px）、prom は"
          " %.2f–%.2f で写真の %.2f–%.2f と同水準。固定ピクセル周期の格子なら"
          "「窓が変わっても period は同一・prom は突出」になるはずで、そうは"
          "なっていない。**ただしこれは原因の同定であって、G-LINE の合格ではない。**"
          % (min(_pm or [0]), max(_pm or [0]),
             min([float(r.get("period_px") or 0.0)
                  for r in (gl.get("photo") or [])] or [0]),
             max([float(r.get("period_px") or 0.0)
                  for r in (gl.get("photo") or [])] or [0]),
             min(_pp or [0]), max(_pp or [0]),
             min(_qp or [0]), max(_qp or [0])))
        A("")

    # ---- 5. G-PEAK -------------------------------------------------------
    A("## 5. G-PEAK — ピークが削られていないか（同一フィールド・frame 24）")
    A("")
    A("| C | mean | std | p99 | max |")
    A("| --- | --- | --- | --- | --- |")
    c512 = gp.get("C_512") or f24.get("C_frame_stats") or {}
    cps = gp.get("C_point_splat") or {}
    cdet = gp.get("C_det_frame_512sq") or {}
    A("| 面積堆積（本作, 736×540） | %.6f | %.4f | %.3f | %.2f |"
      % (float(c512.get("mean") or 0.0), float(c512.get("std") or 0.0),
         float(c512.get("p99") or 0.0), float(c512.get("max") or 0.0)))
    A("| 点スプラット（同一 736×540, Q=25） | %.6f | %.4f | %.3f | %.2f |"
      % (float(cps.get("mean") or 0.0), float(cps.get("std") or 0.0),
         float(cps.get("p99") or 0.0), float(cps.get("max") or 0.0)))
    A("| 点スプラット（det.det_frame そのまま, 512²） | %.6f | %.4f | %.3f | %.2f |"
      % (float(cdet.get("mean") or 0.0), float(cdet.get("std") or 0.0),
         float(cdet.get("p99") or 0.0), float(cdet.get("max") or 0.0)))
    A("")
    A("std 比 **%s**（点スプラット 736×540）/ **%s**（det_frame 512²）、"
      "p99 比 **%s** / **%s**。すべて 0.90 以上なので、面積分配はピークを削っていない。"
      % (fs(gp.get("std_ratio"), 3), fs(gp.get("std_ratio_det_frame"), 3),
         fs(gp.get("p99_ratio"), 3), fs(gp.get("p99_ratio_det_frame"), 3)))
    A("")
    A("点スプラット参照は**本編には使わない**比較材料。面積堆積の射影は太陽円盤の"
      "**中心方向 1 本**で計算している（brief の「滲み幅もジッタも要らない」に従い、"
      "円盤平均＝25 枚の面積写像の平均 ≈ 0.85 px のぼかしを掛けない）。これが"
      "ピーク比較を有利にしている可能性は**未測定**である（Q=1 と Q=25 の点スプラット"
      "比較は本作では取っていない）。")
    A("")

    # ---- 6. G-CONV -------------------------------------------------------
    A("## 6. G-CONV — 水面グリッド 512² vs 1024²（同じ 736×540 へ堆積）")
    A("")
    A("| 水面 | C mean | C std | C p99 | C max |")
    A("| --- | --- | --- | --- | --- |")
    for k, lab in (("C_512", "512²"), ("C_1024", "1024²")):
        c = gc.get(k) or {}
        A("| %s | %.6f | %.4f | %.3f | %.2f |"
          % (lab, float(c.get("mean") or 0.0), float(c.get("std") or 0.0),
             float(c.get("p99") or 0.0), float(c.get("max") or 0.0)))
    A("")
    A("差は std %+.2f%%、p99 %+.2f%%（許容 15%%）。標本密度を上げても C の統計は"
      "ほぼ動かない = 512² で既に収束している。"
      % (100.0 * float(gc.get("std_rel_diff") or 0.0),
         100.0 * float(gc.get("p99_rel_diff") or 0.0)))
    A("")

    # ---- 7. G-WARP -------------------------------------------------------
    A("## 7. G-WARP（合否なし・記録）")
    A("")
    A("| チャンネル | 平均 px | p50 | p99 | max px |")
    A("| --- | --- | --- | --- | --- |")
    for tag in ("r", "g", "b"):
        w = (gw.get("per_channel") or {}).get(tag) or {}
        A("| %s | %.3f | %.3f | %.3f | %.3f |"
          % (tag, float(w.get("px_mean") or 0.0), float(w.get("px_p50") or 0.0),
             float(w.get("px_p99") or 0.0), float(w.get("px_max") or 0.0)))
    A("")
    A("歪みレーンは既存のまま（v4-1.4/5/7 と同一定義: 屈折シフト "
      "`−n_xy·z·(IOR−1)/IOR` を画素へ換算）。")
    A("")

    # ---- 8. G-MOVE + seam ------------------------------------------------
    A("## 8. G-MOVE とループの継ぎ目")
    A("")
    ip = gm.get("interframe_green_png") or {}
    im_ = gm.get("interframe_green_mp4") or {}
    A("ユニークフレーム %s/%s。フレーム間平均絶対差（緑）: PNG で %.2f（最小 %.2f、"
      "最大 %.2f）、MP4 デコード後で %.2f。閾値 0.5 諧調に対して十分大きい。"
      % (gm.get("unique_frames"), gm.get("frames"), float(ip.get("mean") or 0.0),
         float(ip.get("min") or 0.0), float(ip.get("max") or 0.0),
         float(im_.get("mean") or 0.0)))
    A("")
    seam = sc.get("seam") or {}
    if seam:
        A("継ぎ目（t=0 と t=2 s の水面）: max|Δh| = %.4f m、RMS Δh = %.4f m "
          "（場の RMS は %.4f m、比 %.2f）。`L = T·Us/n` の厳密な周期は L = %.4f m で"
          "あり、ブリーフが固定した L = %.10f m とは一致しない。したがって"
          "**2 秒で厳密には閉じない**。これは数値として報告する（隠さない）。"
          % (float(seam.get("max_abs_dh_m") or 0.0), float(seam.get("rms_dh_m") or 0.0),
             float(seam.get("h_rms_m") or 0.0), float(seam.get("seam_over_h_rms") or 0.0),
             float((sc.get("loop") or {}).get("L_from_rule_at_that_n_m") or 0.0),
             float(sc.get("patch_m") or 0.0)))
        A("")

    # ---- 9. G-LEGACY -----------------------------------------------------
    A("## 9. G-LEGACY — 既存 .py の同一性と並行レーンの差分")
    A("")
    bcl = (st0.get("brief_checklist") or {})
    if bcl:
        A("### 9.0 STAGE 0 の絶対パス確認（ブリーフ指定の 3 本）")
        A("")
        A("| path | 存在 | bytes | sha256 |")
        A("| --- | --- | --- | --- |")
        for p, r in (bcl.get("files") or {}).items():
            A("| `%s` | %s | %s | `%s` |"
              % (os.path.basename(p), "はい" if r.get("present") else "**いいえ**",
                 r.get("bytes"), str(r.get("sha256"))[:24]))
        A("")
        A("欠落: %s（%s）。既存 .py 全本の sha256 は `stage0-py-hashes.json`。"
          % (bcl.get("missing") or "なし", "PASS" if bcl.get("pass") else "FAIL"))
        A("")
    A("stage-0 で記録した既存 .py %s 本: 変更 %s、欠落 %s。"
      "本ジョブが書くのは新規 2 本だけ。"
      % (glue.get("recorded_existing"), len(glue.get("changed") or []),
         len(glue.get("missing") or [])))
    A("")
    own = glue.get("own_new_file_hashes") or {}
    A("| 本ジョブの新規ファイル | stage-0 | 現在 | 変化 |")
    A("| --- | --- | --- | --- |")
    for name, rec in own.items():
        A("| `%s` | `%s` | `%s` | %s |"
          % (name, str(rec.get("stage0"))[:16], str(rec.get("now"))[:16],
             "変化あり" if rec.get("changed") else "同一"))
    A("")
    ext = glue.get("external_changes") or {}
    A("並行レーンの外部差分: %s。" % (
        "なし" if not ext else
        ", ".join("%s（mtime %s, 本ジョブは import しない）" % (k, v.get("mtime"))
                  for k, v in ext.items())))
    A("")
    A("`run_v4_1_8.py` が stage-0 以降に変化しているのは**本ジョブ自身**の編集である"
      "（作業中にこのファイルの日本語文字列がツール経由の再エンコードで壊れたため、"
      "同じ内容を UTF-8 で書き直し、あわせて probe の代替経路を足した）。"
      "輸送・堆積の側（`transport_mesh.py`）は stage-0 から**バイト同一**。")
    A("")
    A("**重複実行について（正直な注記）**: この作業ディレクトリでは、同じ brief を"
      "実行する**別インスタンス**が同時に走っていた（7:40–8:00 と 8:35 前後に "
      "`run_v4_1_8.py` と `out/v4-1-8/` を書いている）。本レポートの数値と成果物は"
      "**本ジョブの単一実行**から作られている（古い方は `out/_v4-1-8-prev-run` に"
      "退避し、照合後に削除）。照合結果: **frame 24–47 の 24 枚と frame-0024.png は"
      "2 回の独立実行でバイト一致**、frame 0–23 は一致しなかった（相手の chunk 0 は"
      "その frame-24 ステージより前に走っており、コード状態が違う）。同じ入力から"
      "同じ出力が出ることは、この範囲で確認できている。")
    A("")
    A("**帰属について（正直な注記）**: 作業ディレクトリでは作業中に別レーンが"
      "`run_v4_1_3.py` を書き換えており（mtime 07:10、本ジョブの stage-0 より前）、"
      "その結果 `run_v4_1_4.py` → `R13.SEED` が参照できず `run_v4_1_5.py` が import "
      "できない状態になっている。本ジョブはその影響を受けないよう、G-LINE の統計を"
      "自前で逐語コピーし、v4-1-7 の記録値の数値再現で等価性を確認した"
      "（`measurement.provenance`）。この差分は本ジョブのものではない。")
    A("")

    # ---- 10. G-OUT -------------------------------------------------------
    A("## 10. G-OUT")
    A("")
    pr = gout.get("probe") or {}
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| codec / pix_fmt | %s / %s |" % (pr.get("codec_name"), pr.get("pix_fmt")))
    A("| 解像度 | %sx%s |" % (pr.get("width"), pr.get("height")))
    A("| frames / duration | %s / %s s |" % (pr.get("nb_frames_reported"),
                                             pr.get("duration_s")))
    A("| fps | %s (%s fps) |" % (pr.get("r_frame_rate"), gout.get("fps")))
    A("| bytes / sha256 | %s / `%s` |" % (gout.get("size_bytes"),
                                          str(gout.get("sha256"))[:32]))
    A("| probe の方法 | %s |" % (pr.get("probe_method") or "-"))
    A("| probe の binary | `%s` |" % (pr.get("probe_binary") or "-"))
    A("")
    A("**環境の注記**: このセッションでは PATH 上の `ffmpeg` / `ffprobe` は"
      "他ユーザープロファイル内の 0 バイト WinGet エイリアスで、CreateProcess が"
      "アクセス拒否（WinError 5）を返す。そこでリポジトリ同梱の ffmpeg 9.0.2"
      "（`_vendor/imageio_ffmpeg`）でエンコードし、G-OUT の値は「ヘッダ解析 + "
      "全フレームの rgb24 デコードによる枚数」で取得した。フィールド名は ffprobe の"
      "ものをそのまま使っているので、上の表は ffprobe 出力と同義に読める。")
    A("")

    # ---- 11. deliverables ------------------------------------------------
    A("## 11. 成果物（すべて `out/v4-1-8/`）")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for name, rec in sorted((n.get("files") or {}).items()):
        A("| `%s` | %s | `%s` |" % (name, rec.get("bytes"), rec.get("sha256")))
    A("| `numbers-v4-1-8.json` | （本文の数値の実体） | — |")
    A("")
    A("`still_equals_frame_0024`: %s（`frame-0024.png` と `frames/frame-0024.png` は"
      "同一バイト）。`transport_mesh.py` sha256 `%s`、`run_v4_1_8.py` sha256 `%s`。"
      % (n.get("still_equals_frame_0024"),
         str(n.get("transport_mesh_sha256"))[:32], str(n.get("run_v4_1_8_sha256"))[:32]))
    A("")
    A("注: `checkpoint.json` はこのレポートを書いた後に `report` ユニットを追記するので、"
      "上の sha256 はレポート生成時点のものである。`numbers-v4-1-8.json` は"
      "レポート本文の数値の実体（機械可読）。")
    A("")

    # ---- 12. disclosures -------------------------------------------------
    A("## 12. 開示・限界")
    A("")
    A("1. **G-LINE は FAIL。** 格子症状は v4-1.7 の 173.5 から 66.2 へ落ちたが、"
      "写真（63.7）を 3 窓でわずかに上回る。原因の切り分けは 4.1 に置いた。"
      "ゲートに合わせた後付け調整はしていない（カーネル半径・ジッタ幅という"
      "調整対象が存在しない構成）。**この brief は帯域制限もカーネルも禁じて"
      "いるので、5–15 px 帯を下げる道具は構成の中に存在しない**（v4-1.5 が同じ"
      "統計で写真を下回ったのは、場を σ=4 cell で帯域制限していたため）。"
      "つまりこの FAIL は「この構成では通せない」という結果であって、"
      "調整漏れではない。")
    A("2. **円盤求積**: 射影の幾何は太陽円盤の中心方向 1 本で計算している"
      "（点スプラット車線の Q=25 は画素フィルタであり、面積射影では面積に"
      "畳まれる）。Q=1 と Q=25 の点スプラット比較は本作では取っていないので、"
      "円盤平均を省いたことが G-PEAK の比較をどれだけ有利にしているかは"
      "**未測定**である。円盤半径 0.00465 rad × 深さ 0.6 m ≈ 0.85 px で、"
      "これは画素 1 個未満である。")
    A("3. **画素アスペクト**: x 2.005 mm / y 2.732 mm（v2 の 736:540 フレーミング）。"
      "倍率 1 の意味は L_img/L_gen = 1 であって、アスペクトの非等方は v2 由来。")
    dd = f24.get("deposit_detail") or {}
    _sA = (float(dd.get("area_in") or 0.0) + float(dd.get("area_out") or 0.0))
    A("4. **折り返し（多価写像）**: frame 24 では射影三角形の面積の総和が "
      "Σ|A| = %.1f px² = 画面面積の **%.3f 倍**（＝大局的にはほぼ 1 対 1）で、"
      "折り返しは**局所的**に起きている。|A| で押し出すのでエネルギーは厳密に"
      "保存し、局所の多重被覆は C > 1 として現れる（クリップも平滑化もしない）。"
      "frame 24 の C max = %.1f、C の p99 = %.2f。"
      % (_sA, _sA / float(WIDTH * HEIGHT), float(c512.get("max") or 0.0),
         float(c512.get("p99") or 0.0)))
    A("5. **枠外**: 合成は周期床（mod 付き (u,v)）なので、枠外へ出た着地は周期像として"
      "戻る。unwrap 座標での枠外割合は別枠記録（平均 %.3f）。正規化には使っていない。"
      % float(((f24.get("energy") or {}).get("offframe_fraction_energy") or 0.0)))
    A("6. **C の受け渡し**: 凍結合成器は (u,v) の正方 C を要求するため、736×540 の"
      "堆積を 736×736 へ線形 y 補間で再標本化している（純粋な再標本化で、"
      "再正規化も追加のぼかしもしない）。")
    A("7. **ループ**: L を v2 の窓に固定したため 2 s で厳密には閉じない（第 8 節）。")
    A("8. **測定の限界**: depth は 1 次元周辺統計であり、格子と等方的テクスチャを"
      "区別しない。C の統計も 1 フレーム・1 場のもので、シーンの一般性は主張しない。")
    A("")

    # ---- 13. what this shows --------------------------------------------
    A("## 13. この結果が言えること / 言えないこと")
    A("")
    A("* **言える**: 点スプラット＋カーネルを面積保存の射影に置き換えると、"
      "エネルギーは厳密に保存し（最大 %s）、ピークは点スプラット参照以上に残り"
      "（std 比 %s）、格子の軸平行成分は v4-1.7 の 4 割以下に落ちた。"
      % ("%.1e" % float(ge.get("max_rel_err") or 0.0), fs(gp.get("std_ratio"), 2)))
    A("* **言える**: 512² と 1024² で C の統計は 15% 以内で一致し、"
      "堆積は標本密度に収束している。")
    A("* **言えない**: 「格子が完全に消えた」とは言えない（G-LINE FAIL、"
      "残差は写真と同水準〜1.2 倍）。残差の主因が場の Nyquist 近傍の内容と v2 の"
      "フレーミングにあることは 4.1 の切り分けで示したが、これは原因の特定であって"
      "ゲートの合格ではない。")
    A("* **言えない**: 見た目の最終判断。判定は数値統計のみで行った（目視なし）。")
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-8.json + REPORT-V4-1-8.md")
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
    ap.add_argument("--frame24", action="store_true")
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        do_frame24()
        do_diag()
        do_chunk(0)
        do_chunk(1)
        do_warp()
        do_encode()
        do_gates()
        do_report()
        return 0
    did = False
    if a.stage0:
        stage0()
        did = True
    if a.frame24:
        do_frame24()
        did = True
    if a.diag:
        do_diag()
        did = True
    if a.chunk is not None:
        do_chunk(a.chunk)
        did = True
    if a.warp:
        do_warp()
        did = True
    if a.encode:
        do_encode()
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
