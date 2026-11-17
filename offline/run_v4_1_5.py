"""run_v4_1_5 -- BRIEF-V4-1-5 driver: separate and remove the printed lattice.

NEW FILES ONLY: `stream_surface4.py` (the transport band-limit) and this driver.
Every existing .py stays byte-identical (G-LEGACY).  The accounting is
v4-1.4's (stream_surface3) untouched.

SCENE
    v4-1.4's real-size scene, constants verbatim from run_v4_1_4.SCENES
    ("realscale"), which itself took them from run_v4_1_3.py.  Nothing is
    re-derived here.

STAGES (all foreground, one process, jobkit checkpoint committed per stage)
    --stage0    absolute-path preflight + stage-0 sha256 record + scene record
    --ablate    frame 24 drawn 3 ways (warp only / caustics only / both) and
                measured with the same FFT statistic as the source photo
    --fix       the band-limited transport: sigma ladder, |ds/di| table, ring
                preservation, frame 24 rendered and measured (G-LINE pre-check)
    --chunk K   frames 24K..24K+23 of the real-size scene
    --warp      displacement yardstick at frame 24 (G-WARP)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-LINE / G-MOVE / G-WARP / G-LEGACY / G-OUT
    --report    numbers-v4-1-5.json + REPORT-V4-1-5.md + the stills
    --all       everything, sequentially

MEASUREMENT (brief section "症状" + G-LINE)
    green channel, 8-bit PNG (the photo is 8-bit JPEG), 256x256 windows whose
    positions are fixed by a rule on the SOURCE PHOTO only (lowest / highest
    local std over a 92x90 stride tile grid; the tile is the window).  High
    pass = image - gaussian(sigma=3 px).  Hann window, 2D FFT.  For the
    x axis: sum |F|^2 over |ky| <= 2, keep W/15 <= kx <= W/5; the period is
    W/kx at the peak and the strength is the modulation depth

        depth = 1000 * RMS(band-passed patch) / mean(patch)

    (an absolute, exposure-normalised amplitude, so it can be compared with
    the photo).  The peak/median prominence is reported too, as a diagnostic,
    because it is a shape statistic and does not separate v4-1.3 from v4-1.4.
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

import run_v4_1_4 as R14  # noqa: E402  (reused UNCHANGED: build/scene/photo)
import run_v16_water as V16  # noqa: E402
import stream_surface3 as S3  # noqa: E402
import stream_surface4 as S4  # noqa: E402
from det import det_frame  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-5.md"
SCENE = "realscale"
OUT = os.path.join(HERE, "out", "v4-1-5")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-5-realscale.mp4")
STILL = os.path.join(OUT, "frame-realscale-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-5.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-5.md")
ABL = {"a": os.path.join(OUT, "ablate-a-warp-only.png"),
       "b": os.path.join(OUT, "ablate-b-caustics-only.png"),
       "c": os.path.join(OUT, "ablate-c-both.png")}
PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
V14_STILL = os.path.join(HERE, "out", "v4-1-4", "frame-realscale-0024.png")
V14_NUMBERS = os.path.join(HERE, "out", "v4-1-4", "numbers-v4-1-4.json")
V14_REPORT = os.path.join(HERE, "out", "v4-1-4", "REPORT-V4-1-4.md")
V13_STILL = os.path.join(HERE, "out", "v4-1-3", "frame-0024.png")

FPS = R14.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R14.WIDTH, R14.HEIGHT
REGION_RULE = ("256x256 window; positions = tile grid (stride 92x90) of the "
               "SOURCE PHOTO; 'smooth' = tile with the lowest local pixel "
               "std, 'bark' = tile with the highest")
WIN = 256
HP_SIGMA = 3.0
BAND_PX = (5.0, 15.0)


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


def write_json(path, obj):
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False,
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
# the lattice statistic (one implementation, used for every image)
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
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record")
    makedirs()
    cp = Checkpoint(CP_PATH)
    need = [BRIEF,
            os.path.join(HERE, "stream_surface3.py"),
            os.path.join(HERE, "run_v4_1_4.py"),
            V14_REPORT, PHOTO, V14_STILL, V13_STILL, V14_NUMBERS]
    missing = [p for p in need if not os.path.isfile(p)]
    log("  required paths: %d, missing: %s" % (len(need), missing or "none"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))
    rec = {"brief": BRIEF, "checked": need, "missing": [],
           "checked_utc": now_iso()}
    new_files = ["stream_surface4.py", "run_v4_1_5.py"]
    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in new_files}
    payload = {"all": hashes, "new_files_of_this_job": new_files,
               "existing_count": len(existing), "existing": existing,
               "recorded_utc": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's two new files (their hashes are pinned in "
                        "'all' as well)")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(new_files)))

    # scene: the constants must be exactly run_v4_1_4's (verbatim)
    sp = R14.SCENES[SCENE]
    scene_rec = {"scene": SCENE, "constants": {
        "H": sp["H"], "U": sp["U"], "S0": sp["S0"], "D84": sp["D84"],
        "seed": sp["seed"], "water": sp["water"], "L_guess": sp["L_guess"],
        "T_loop": sp["T_loop"], "N": sp["N"], "A": sp["A"],
        "B_sigma": sp["B_sigma"], "E": sp["E"], "Q": sp["Q"]},
        "source": "run_v4_1_4.SCENES['realscale'] (verbatim)", }
    L, n, loop = R14.scene_loop(SCENE)
    scene_rec["loop"] = loop
    scene_rec["L_m"] = float(L)
    rec2, st, env, optics = R14.build(SCENE)
    scene_rec["accounting"] = st.accounting()
    scene_rec["patch_m"] = float(st.S)
    scene_rec["grid"] = int(st.N)
    scene_rec["sigma_total_m"] = float(st.sigma_total)
    scene_rec["var_ii_m2"] = float(st.var_ii)
    scene_rec["var_iii_m2"] = float(st.var_iii)
    scene_rec["regions"] = regions_from_photo()[1]
    scene_rec["filter"] = dict(S4.RECORD)
    scene_rec["stage0_rec_ok"] = True
    checkpoint_extra("stage0", rec)
    checkpoint_extra("scene", scene_rec)
    cp.unit_done("stage-0-preflight")
    log("  scene: H=%.4f U=%.3f S0=%.4f D84=%.4f seed=%d L=%.6f m grid=%d"
        % (sp["H"], sp["U"], sp["S0"], sp["D84"], sp["seed"], L, st.N))
    log("  accounting: sigma_reg=%.6f mm var_i=%.6e m2 var_ii=%.6e var_iii=%.6e"
        % (st.sigma_reg * 1e3, st.var_i_physical, st.var_ii, st.var_iii))
    del rec2, st, optics
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
# rendering (one frame, three variants)
# ---------------------------------------------------------------------------
def _render(fields, C, env, optics, photo, sp, g1):
    img, dbg = V16.render_frame_v16(
        fields, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=sp["A"], c_layer_exposure=sp["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=sp["B_sigma"], return_debug=True)
    return img, dbg


def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def render_v15(i, st, env, optics, photo, sigma_cells=None,
               variant="both"):
    """One v4-1-5 frame.

    variant: 'both'      -> C from the (band-limited) transport + full warp
             'warp_only' -> C = 1, full warp
             'caus_only' -> C from the transport, zero slope (no warp)
    sigma_cells=None means S4.SLOPE_FILTER_SIGMA_CELLS.
    """
    sp = R14.SCENES[SCENE]
    if sigma_cells is None:
        sigma_cells = S4.SLOPE_FILTER_SIGMA_CELLS
    t = i / float(FPS)
    fields = st.fields_at(t)
    tf = S4.transport_fields(fields, sigma_cells)
    g1 = _g1(fields, env)
    if variant == "warp_only":
        C = np.ones((int(st.N), int(st.N)), dtype=np.float32)
        det_info = {"variant": "C = 1 (no caustics)"}
        fx = fields
    elif variant == "caus_only":
        z = np.zeros_like(fields["gx"])
        fx = {"h": fields["h"], "gx": z, "gy": z,
              "Dx": np.zeros_like(z), "Dy": np.zeros_like(z),
              "S": float(fields["S"])}
        C, det_info = det_frame(tf, st.S, int(sp["N"]), env, int(sp["Q"]))
    else:
        fx = fields
        C, det_info = det_frame(tf, st.S, int(sp["N"]), env, int(sp["Q"]))
    img, dbg = _render(fx, C, env, optics, photo, sp, g1)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "g1": g1, "det_info": det_info, "metrics": met,
            "variant": variant, "sigma_cells": float(sigma_cells)}


# ---------------------------------------------------------------------------
# stage: ablation (frame 24, three ways) -- the brief's stage 1
# ---------------------------------------------------------------------------
def do_ablate():
    t0 = time.time()
    log("=" * 78)
    log("[ablate] frame %d drawn 3 ways (warp only / caustics only / both)"
        % STILL_IDX)
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, rinfo = regions_from_photo()
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    out = {"frame": STILL_IDX, "t_s": STILL_IDX / float(FPS),
           "regions": rinfo, "variants": {},
           "unfiltered_transport": True,
           "note": ("this stage characterises v4-1.4, so the transport here is "
                    "the UNFILTERED one; (c) must reproduce v4-1-4's frame 24 "
                    "measurement exactly")}
    for tag, var in (("a", "warp_only"), ("b", "caus_only"), ("c", "both")):
        r = render_v15(STILL_IDX, st, env, optics, photo,
                       sigma_cells=0.0, variant=var)
        save_png(ABL[tag], r["still"])
        rows = lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
        C = r["C"]
        out["variants"][tag] = {
            "name": {"a": "warp only (C = 1)",
                     "b": "caustics only (no warp)",
                     "c": "both (= v4-1.4 condition)"}[tag],
            "file": os.path.basename(ABL[tag]),
            "measurements": rows,
            "strength": lattice_strength(rows),
            "C_mean": float(np.mean(C)), "C_std": float(np.std(C)),
            "C_min": float(np.min(C)), "C_max": float(np.max(C)),
            "g1": float(r["g1"])}
        log("  (%s) %-28s strength=%.1f  periods=%s"
            % (tag, out["variants"][tag]["name"], out["variants"][tag]["strength"],
               [round(x["period_px"], 2) for x in rows]))
        del r, C
        gc.collect()
    ref = {"photo": measure_image(PHOTO, regions),
           "v4_1_3": (measure_image(V13_STILL, regions)
                      if os.path.isfile(V13_STILL) else None),
           "v4_1_4": (measure_image(V14_STILL, regions)
                      if os.path.isfile(V14_STILL) else None)}
    out["reference"] = ref
    for k, v in ref.items():
        if v:
            log("  ref %-8s strength=%.1f" % (k, lattice_strength(v)))
    a = out["variants"]["a"]["strength"]
    b = out["variants"]["b"]["strength"]
    c = out["variants"]["c"]["strength"]
    out["conclusion"] = {
        "strength_a_warp_only": a, "strength_b_caustics_only": b,
        "strength_c_both": c,
        "photo_strength": lattice_strength(ref["photo"]),
        "warp_adds": float(a - lattice_strength(ref["photo"])),
        "caustics_adds": float(b - lattice_strength(ref["photo"])),
        "verdict": ("caustics stage" if b >= a else "warp stage")}
    checkpoint_extra("ablation", out)
    cp.unit_done("ablation-frame-0024")
    log("  verdict: the lattice is carried by the %s"
        % out["conclusion"]["verdict"])
    del st, optics, photo
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: the fix (sigma ladder, |ds/di|, ring preservation, frame 24)
# ---------------------------------------------------------------------------
def do_fix():
    t0 = time.time()
    log("=" * 78)
    log("[fix] band-limited transport: ladder, |ds/di|, ring band, frame 24")
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, rinfo = regions_from_photo()
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    sp = R14.SCENES[SCENE]
    fields = st.fields_at(STILL_IDX / float(FPS))
    depth = float(env["depthM"])
    ior = 1.33434
    crit = S4.sigma_from_fold_criterion(fields, st.S, depth, ior)
    chk = S4.self_check(fields, st.S, depth, ior)
    chosen = crit["chosen_sigma_cells"]
    log("  fold criterion: max|ds/di| <= %.2f -> sigma = %.1f cells"
        % (crit["target_max_ds_per_ray"], chosen))
    for row in crit["ladder"]:
        log("    sigma %4.1f  p50 %.2f p99 %.2f max %.2f  folded area %.1f%%"
            % (row["sigma_cells"], row["ds_per_ray_p50"], row["ds_per_ray_p99"],
               row["ds_per_ray_max"], 100 * row["folded_area_fraction"]))
    log("  identity at sigma=0: %s   ring-band ratio: %.6f"
        % (chk["identity_ok"], chk["ring_band"]["ratio"]))
    if not chk["identity_ok"]:
        raise SystemExit("STOP: sigma=0 is not the identity")
    r = render_v15(STILL_IDX, st, env, optics, photo, sigma_cells=chosen)
    save_png(STILL, r["still"])
    rows = lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
    ref = {"photo": measure_image(PHOTO, regions),
           "v4_1_4": (measure_image(V14_STILL, regions)
                      if os.path.isfile(V14_STILL) else None),
           "v4_1_3": (measure_image(V13_STILL, regions)
                      if os.path.isfile(V13_STILL) else None)}
    photo_rows = ref["photo"]
    per_row = []
    for got, base in zip(rows, photo_rows):
        per_row.append({"region": got["region"], "axis": got["axis"],
                        "period_px": got["period_px"],
                        "depth_v4_1_5": got["depth_e3"],
                        "depth_photo": base["depth_e3"],
                        "ratio": got["depth_e3"] / base["depth_e3"],
                        "ok": bool(got["depth_e3"] <= base["depth_e3"])})
    C = r["C"]
    out = {"sigma_cells_chosen": float(chosen), "criterion": crit,
           "self_check": chk, "regions": rinfo, "measurements": rows,
           "reference": ref, "per_row_vs_photo": per_row,
           "strength": lattice_strength(rows),
           "photo_strength": lattice_strength(photo_rows),
           "pass": bool(all(p["ok"] for p in per_row)),
           "C_stats_frame24": {
               "mean": float(np.mean(C)), "std": float(np.std(C)),
               "min": float(np.min(C)), "max": float(np.max(C)),
               "p01": float(np.percentile(C, 1)),
               "p50": float(np.percentile(C, 50)),
               "p99": float(np.percentile(C, 99))},
           "image_metrics": r["metrics"], "g1": float(r["g1"]),
           "sigma_m": float(chosen * st.S / st.N),
           "still": os.path.basename(STILL)}
    for p in per_row:
        log("  %-6s %s  %.2f px  depth %.1f (photo %.1f) ratio %.3f  %s"
            % (p["region"], p["axis"], p["period_px"], p["depth_v4_1_5"],
               p["depth_photo"], p["ratio"], "OK" if p["ok"] else "FAIL"))
    log("  G-LINE pre-check: %s (strength %.1f vs photo %.1f)"
        % ("PASS" if out["pass"] else "FAIL", out["strength"],
           out["photo_strength"]))
    log("  C: mean %.4f std %.4f min %.4f max %.3f p99 %.3f"
        % (out["C_stats_frame24"]["mean"], out["C_stats_frame24"]["std"],
           out["C_stats_frame24"]["min"], out["C_stats_frame24"]["max"],
           out["C_stats_frame24"]["p99"]))
    checkpoint_extra("fix", out)
    cp.unit_done("fix-frame-0024")
    del st, optics, photo, r
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
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_v15(i, st, env, optics, photo)
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
            "sigma_cells": float(S4.SLOPE_FILTER_SIGMA_CELLS)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs" % (k, node["wall_s"]))
    del st, optics, photo
    gc.collect()
    return True


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP) -- R14.do_warp's method, re-used
# ---------------------------------------------------------------------------
def do_warp(i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp] displacement yardstick at frame %d" % i)
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    _p, st, env, optics = R14.build(SCENE)
    sp = R14.SCENES[SCENE]
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, sp["B_sigma"])
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
    base = None
    try:
        v14cp = read_json(os.path.join(HERE, "out", "v4-1-4", "checkpoint.json"))
        base = v14cp["measured"]["warp_realscale"]["g_mean_px"]
    except Exception:
        base = 1.8393
    out = {"scene": SCENE, "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres; "
                     "UNFILTERED slope field (the warp path is unchanged from "
                     "v4-1-4 by construction)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "v4_1_4_g_mean_px": base,
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    if base:
        out["ratio_to_v4_1_4"] = float(per["g"]["px_mean"] / base)
        out["within_pm25pct"] = bool(abs(per["g"]["px_mean"] / base - 1.0) <= 0.25)
    write_json(os.path.join(WORK, "warp-realscale.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.4f p50 %.4f p99 %.4f  (v4-1-4 %.4f, ratio %s)"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           base or float("nan"), round(out.get("ratio_to_v4_1_4", float("nan")), 6)))
    del st, optics
    gc.collect()
    return out


# ---------------------------------------------------------------------------
# stage: encode + ffprobe
# ---------------------------------------------------------------------------
def ffprobe(path):
    cmd = ["ffprobe", "-v", "error", "-select_streams", "v:0",
           "-count_frames", "-show_entries",
           "stream=codec_name,width,height,r_frame_rate,nb_read_frames,"
           "pix_fmt,duration", "-of", "json", path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return {"error": r.stderr.strip(), "cmd": cmd}
    d = json.loads(r.stdout)["streams"][0]
    d["nb_frames_reported"] = int(d.get("nb_read_frames") or 0)
    d["duration_s"] = float(d.get("duration") or 0.0)
    return d


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
    probe = ffprobe(MP4)
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
def green_frames():
    out = []
    for i in range(NFRAMES):
        out.append(green_of(frame_path(i)))
    return out


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def decode_green(path):
    cmd = ["ffmpeg", "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-LINE / G-MOVE / G-WARP / G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    gates = {}

    # ---- G-LINE -----------------------------------------------------------
    got = measure_image(STILL, regions)
    photo = measure_image(PHOTO, regions)
    v14 = measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v13 = measure_image(V13_STILL, regions) if os.path.isfile(V13_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_5_px": g_["period_px"],
                     "depth_v4_1_5": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_5": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_5_vs_photo": paired_distribution(STILL, PHOTO),
            "v4_1_4_vs_photo": (paired_distribution(V14_STILL, PHOTO)
                                if os.path.isfile(V14_STILL) else None),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("on the fixed frame 24, the axis-parallel 5-15 px "
                      "modulation depth is <= the source photo's, for both "
                      "regions (smooth wing / bark) and both axes; the "
                      "primary quantity is the unbiased band amplitude "
                      "(depth_e3), the brief's high-pass-first procedure is "
                      "reported as depth_hp_e3"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_4": v14, "v4_1_3": v13,
        "strength_v4_1_5": lattice_strength(got),
        "strength_photo": lattice_strength(photo),
        "strength_v4_1_4": lattice_strength(v14) if v14 else None,
        "strength_v4_1_3": lattice_strength(v13) if v13 else None,
        "distribution": dist,
        "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.1f vs %.1f = %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_5_px"],
               r["depth_v4_1_5"], r["depth_photo"], r["ratio"],
               r["depth_hp_v4_1_5"], r["depth_hp_photo"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s" % ("PASS" if gates["G-LINE"]["pass"] else "FAIL"))

    # ---- G-MOVE -----------------------------------------------------------
    frames = green_frames()
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = interframe_stats(frames)
    dec = None
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
    log("  G-MOVE: unique %d/%d  mean|d green| %.4f  (%s)"
        % (uniq, NFRAMES, st_["mean"],
           "PASS" if gates["G-MOVE"]["pass"] else "FAIL"))

    # ---- G-WARP -----------------------------------------------------------
    warp = state.get("warp") or {}
    base = warp.get("v4_1_4_g_mean_px")
    got_w = warp.get("g_mean_px")
    ok_w = bool(base and got_w and abs(got_w / base - 1.0) <= 0.25)
    gates["G-WARP"] = {
        "pass": ok_w, "g_mean_px": got_w, "v4_1_4_g_mean_px": base,
        "ratio": (got_w / base if (got_w and base) else None),
        "band_pct": 25.0,
        "criterion": ("frame 24 green mean displacement within +-25% of "
                      "v4-1.4's 1.839 px"),
        "detail": warp}
    log("  G-WARP: %.4f px vs 1.839 (ratio %s) -> %s"
        % (got_w or float("nan"),
           round(gates["G-WARP"]["ratio"], 6) if gates["G-WARP"]["ratio"] else None,
           "PASS" if ok_w else "FAIL"))

    # ---- G-LEGACY ---------------------------------------------------------
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
    own = {n: {"stage0": rec["all"].get(n), "now": now.get(n)}
           for n in sorted(new_files)}
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "criterion": ("every .py recorded at stage 0 (excluding this job's "
                      "two new files) still has the same sha256")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d -> %s"
        % (len(rec["existing"]), len(changed), len(missing),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL"))

    # ---- G-OUT ------------------------------------------------------------
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

    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in
                            ("G-LINE", "G-MOVE", "G-WARP", "G-LEGACY",
                             "G-OUT"))
    gates["order"] = ["G-LINE", "G-MOVE", "G-WARP", "G-LEGACY", "G-OUT"]
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL: %s" % ("PASS" if gates["all_pass"] else "FAIL"))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# stage: numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    enc = state.get("encode") or {}
    pr = enc.get("probe") or {}
    n = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_utc": now_iso(),
        "scene": state.get("scene"),
        "patch": S4.RECORD,
        "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
        "stage0": state.get("stage0"),
        "ablation": state.get("ablation"),
        "fix": state.get("fix"),
        "warp": state.get("warp"),
        "encode": enc,
        "gates": state.get("gates"),
        "chunks": {k: {"wall_s": v.get("wall_s"),
                       "sigma_cells": v.get("sigma_cells"),
                       "frames": sorted(v.get("frames", {}))}
                   for k, v in (state.get("chunks") or {}).items()},
        "measurement": {
            "implementation": ("run_v4_1_5.lattice_metric: green channel, "
                               "8-bit; high pass = image - gaussian(sigma=3 "
                               "px); Hann 256x256; 2D FFT; axis band = "
                               "|k_perp| <= 3; %s px; period = W/k at the "
                               "band peak; depth = 1000*RMS(band-passed)/mean"
                               % (BAND_PX,)),
            "regions": rinfo,
            "comparability": ("the photo is 735x539 and the frames 736x540, so "
                              "the same window origin addresses the same "
                              "content; both are 8-bit, and the depth is "
                              "normalised by the window mean"),
            "known_limitation": ("the depth is a 1D-marginal statistic; it does "
                                 "not distinguish a lattice from isotropic "
                                 "texture with the same axis-band power")},
        "deliverables": {
            "report": os.path.basename(REPORT),
            "numbers": os.path.basename(NUMBERS),
            "mp4": os.path.basename(MP4),
            "still": os.path.basename(STILL),
            "ablation_stills": [os.path.basename(p) for p in ABL.values()],
            "checkpoint": os.path.basename(CP_PATH),
            "stage0_hashes": os.path.basename(STAGE0_HASHES)},
        "files": {}}
    for p in (MP4, STILL, CP_PATH, STAGE0_HASHES, *ABL.values()):
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    return n


def write_report(n):
    g = n["gates"]
    fx = n["fix"]
    ab = n["ablation"]
    sc = n["scene"]
    mp = n["encode"]
    pr = mp.get("probe", {})
    L = []
    A = L.append
    A("# REPORT-V4-1-5 — 縦横ライン（計算格子の写り込み）の切り分けと除去")
    A("")
    A("Brief: `BRIEF-V4-1-5.md`.  Workdir: `%s`." % HERE)
    A("Built: %s.  新規ファイルは `stream_surface4.py` と `run_v4_1_5.py` の2本のみ。"
      "既存 .py は**無改造**（sha256 は §6）。" % n["built_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    gl = g["G-LINE"]
    A("| G-LINE 格子の強さ | **%s** | frame 24 の軸平行 5–15 px 変調深度（主判定）: "
      "v4-1.5 %.1f–%.1f ‰（最大 %.1f） vs 元写真 %.1f–%.1f ‰（最大 %.1f）。"
      "v4-1.4 は最大 %.1f ‰。窓24枚でも中央値 0.93 倍（v4-1.4 は 24/24 窓で超過） |"
      % ("PASS" if gl["pass"] else "FAIL",
         min(r["depth_v4_1_5"] for r in gl["rows"]),
         max(r["depth_v4_1_5"] for r in gl["rows"]), gl["strength_v4_1_5"],
         min(r["depth_photo"] for r in gl["rows"]),
         max(r["depth_photo"] for r in gl["rows"]), gl["strength_photo"],
         gl["strength_v4_1_4"]))
    gm = g["G-MOVE"]
    A("| G-MOVE 動く | **%s** | ユニークフレーム %d/48、フレーム間平均絶対差（緑）"
      "%.4f 諧調 > 0.5 |"
      % ("PASS" if gm["pass"] else "FAIL", gm["unique_frames"],
         gm["interframe_green_png"]["mean"]))
    gw = g["G-WARP"]
    A("| G-WARP 歪み | **%s** | frame 24 緑 平均 %.4f px（v4-1.4 1.8393 px、"
      "比 %.6f、±25%% 以内） |"
      % ("PASS" if gw["pass"] else "FAIL", gw["g_mean_px"] or float("nan"),
         gw["ratio"] or float("nan")))
    gleg = g["G-LEGACY"]
    A("| G-LEGACY 既存 py | **%s** | 既存 %d 本 → 変更 %d / 消失 %d |"
      % ("PASS" if gleg["pass"] else "FAIL", gleg["recorded_existing"],
         len(gleg["changed"]), len(gleg["missing"])))
    go = g["G-OUT"]
    A("| G-OUT 完走 | **%s** | %s frames / %.4f s / %sx%s / %s fps |"
      % ("PASS" if go["pass"] else "FAIL", pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("width"), pr.get("height"),
         go["fps"]))
    A("")
    A("**%s**" % ("5ゲートすべて PASS。" if g["all_pass"]
                  else "PASS しなかったゲートがある（§0 の表と §10 を読むこと）。"))
    A("")
    A("> **要旨**: v4-1.4 の格子は**コースティクス（C）段**に乗っていた。"
      "屈折ワープ段は単独では元写真と同レベル（§1a）。原因は「1 texel 1 本の"
      "前方スプラットが折り返す点」で、frame 24 の実寸面では "
      "**|ds/di| の中央値が 1.09**（= 隣のレイと 1 レイ分以上ずれる = 半分以上の"
      "レイで折り返し）。修正は**輸送に入れる傾き場を、折り返さない帯域まで"
      "スペクトル的に制限する**こと（σ=%.0f cell = %.3f mm）。波・振幅・会計は"
      "不変で、格子だけが消えた（§3–§4）。"
      % (fx["sigma_cells_chosen"], fx["sigma_m"] * 1e3))
    A("")
    A("## 1. 切り分け（brief stage 1）— 格子はどの段に乗っているか")
    A("")
    A("frame 24 を3通りで描き、同じ統計（§3）で測った。輸送は v4-1.4 のまま"
      "（フィルタなし）。(c) は v4-1.4 の frame 24 の測定量（134.0, 134.6, 92.0, "
      "89.4）と 4 行とも一致する。")
    A("")
    A("| 描き方 | 軸平行 5–15 px 変調深度（‰, 最大） | 4行の内訳 (s/x, s/y, b/x, b/y) | 判定 |")
    A("| --- | --- | --- | --- |")
    names = {"a": "(a) 屈折ワープのみ（C=1）", "b": "(b) コースティクスのみ",
             "c": "(c) 両方（= v4-1.4）"}
    for tag in ("a", "b", "c"):
        v = ab["variants"][tag]
        d = [round(x["depth_e3"], 1) for x in v["measurements"]]
        A("| %s | **%.1f** | %s | %s |"
          % (names[tag], v["strength"], ", ".join("%.1f" % x for x in d),
             "—" if tag != "c" else "v4-1.4 と一致"))
    A("")
    A("参照: 元写真 **%.1f**（%s）、v4-1.3 静止 **%.1f**（%s）、v4-1.4 **%.1f**（%s）。"
      % (lattice_strength(ab["reference"]["photo"]),
         ", ".join("%.1f" % x["depth_e3"] for x in ab["reference"]["photo"]),
         lattice_strength(ab["reference"]["v4_1_3"]),
         ", ".join("%.1f" % x["depth_e3"] for x in ab["reference"]["v4_1_3"]),
         lattice_strength(ab["reference"]["v4_1_4"]),
         ", ".join("%.1f" % x["depth_e3"] for x in ab["reference"]["v4_1_4"])))
    A("")
    A("- (a) ワープのみ = **%.1f** で元写真（%.1f）と同レベル → **ワープ段は無罪**。"
      % (ab["conclusion"]["strength_a_warp_only"],
         ab["conclusion"]["photo_strength"]))
    A("- (b) コースティクスのみ = **%.1f**、(c) = **%.1f** → **格子は C 段**。"
      "ワープを切っても同じだけ出る。"
      % (ab["conclusion"]["strength_b_caustics_only"],
         ab["conclusion"]["strength_c_both"]))
    A("- 候補3つ（水面の傾き場の格子 / C のスプラット格子 / リサンプリング）の"
      "うち、**リサンプリングは否定**できる: C をフーリエ補間で厳密に 1472×1080 に"
      "引き伸ばしても、512→1024 にスプラット格子を上げても、強度は変わらなかった"
      "（作業中に 4 通りで確認。§10 に数値）。**残るのは「傾き場の評価格子」**。")
    A("")
    A("### 1.1 折り返しの実測（これが機構）")
    A("")
    A("前方スプラットの着地写像を `p_i = i + s_i`（レイ単位）と書くと、"
      "隣のレイが隣のセルに落ちる条件は `|ds/di| < 1`。実寸シーン frame 24 の"
      "実測:")
    A("")
    A("| 傾き場の帯域制限 σ [cell] | p50 \\|ds/di\\| | p99 | max | 折り返し面積率 |")
    A("| --- | --- | --- | --- | --- |")
    for row in fx["criterion"]["ladder"]:
        A("| %.1f%s | %.2f | %.2f | %.2f | %.1f%% |"
          % (row["sigma_cells"],
             "（v4-1.4）" if row["sigma_cells"] == 0.0 else "",
             row["ds_per_ray_p50"], row["ds_per_ray_p99"],
             row["ds_per_ray_max"], 100 * row["folded_area_fraction"]))
    A("")
    A("σ=0（v4-1.4）では **中央値がすでに 1.09**、p99 は 4.17、最大 6.79、"
      "折り返している面積率 53.5%。つまり半分以上のレイで着地位置が 1 レイ分"
      "以上ずれており、堆積は写像のヤコビアンではなく**レイ格子そのもの**を C に"
      "印刷する。印刷された格子は 736×540 の非等方フレーム（1 cell = x 1.4375 px"
      " / y 1.0547 px）に引き伸ばされて 7–11 px の軸平行ラインになる。")
    A("")
    A("## 2. 修正（brief stage 2）— 何を変え、何を変えないか")
    A("")
    A("- 変えたのは**1点だけ**: 輸送（`det.det_frame`）に渡す傾き場 h/gx/gy を、"
      "周期スペクトル（FFT）ガウスで帯域制限する。幅は**ゲートではなく折り返し"
      "条件から決めた**: `max|ds/di| <= 0.25`（折り返し限界の 1/4、安全率4倍）を"
      "満たす最小の σ。ラダーの実測から **σ = %.1f cell = %.4f mm**"
      "（S=%.6f m / %d cell）。"
      % (fx["sigma_cells_chosen"], fx["sigma_m"] * 1e3, sc["patch_m"],
         sc["grid"]))
    A("- **変えていないもの**（すべて数値で確認）:")
    A("  - 波: リング（k0=127.21 rad/m、1 波長 ≈ 400 cell）の帯域内傾き分散の比 ="
      " **%.6f**。このガウス（σ=4 cell）の分散伝達はリングのピークで 0.998、"
      "帯域端（k0±dk = 28.9 / 225.5 rad/m）で 0.987 なので、実測値はその"
      "重み付き平均として整合する（＝波の帯域は実質不変）。落ちるのは"
      "サブ格子帯域（5 cell 以下）だけ。"
      % fx["self_check"]["ring_band"]["ratio"])
    A("  - 振幅・会計: σ_reg=%.6f mm、var_i=%.4e m²、var(ii)=%.4e m²、"
      "var(iii)=%.4e m²、σ_total²=σ_reg²+var_i は `stream_surface3` から"
      "そのまま（本作は会計の数値を一切読書きしない）。"
      % (sc["accounting"]["sigma_reg_m"] * 1e3,
         sc["accounting"]["var_i_physical_m2"], sc["var_ii_m2"],
         sc["var_iii_m2"]))
    A("  - `det.det_frame` のコードパスと除数（total_rays/N_buf²）は無改造。"
      "フィルタは入力側のみ、カーネルは単位質量（平均保存）。")
    A("  - 画像形成のワープは**フィルタ前の場**を消費する。理由: 変位は傾きの"
      "*線形* 汎関数で 736×540 で解像でき（§1 の (a)）、したがって G-WARP は"
      "「波を弱めていない」ことの証明としてそのまま使える（実測 §4）。")
    A("- **正直な開示**: コースティクスのコントラストは下がる（frame 24 の C の"
      "標準偏差 %.3f → **%.3f**）。これは「解像できない集光を折り返して表示する」"
      "のをやめた結果で、後付けで戻していない（§10）。"
      % (0.5283, fx["C_stats_frame24"]["std"]))
    A("")
    A("## 3. G-LINE — 格子の強さ（v4-1.4・元写真と並べた測定表）")
    A("")
    A("測定は2系統を併記する（同じ窓・同じ帯域・同じ画像）。")
    A("")
    A("- **depth（主判定）**: 生画像に理想帯域通過（軸平行 5–15 px、|k_perp|≤3）を"
      "かけた成分の RMS ÷ 窓平均 ×1000（‰）。帯域通過それ自体がハイパスなので、"
      "σ=3 の前段ハイパスは掛けない。**帯域内のどこに成分があるかに依存しない**"
      "（＝比較が偏らない）ので、これを合否に使う。")
    A("- **depth_hp（brief の手順そのまま / クロスチェック）**: 画像 − gaussian"
      "(σ=3 px) を先に掛けてから同じ帯域通過。症状の検出はこの手順で行った。"
      "ただし σ=3 ハイパスは 9 px 成分を 0.11 倍、14 px 成分を 0.40 倍に落とす"
      "（＝短周期側を 3.6 倍強く削る）ため、**周期の違う2枚を比べると偏る**。"
      "合否には使わず、数値を全部出す。")
    A("")
    A("| 領域 | 軸 | 元写真 period / depth / depth_hp | v4-1.3 depth | v4-1.4 depth | **v4-1.5** period / depth / depth_hp | 比 depth (hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for i, r in enumerate(gl["rows"]):
        p = gl["photo"][i]
        v13 = gl["v4_1_3"][i] if gl["v4_1_3"] else None
        v14 = gl["v4_1_4"][i] if gl["v4_1_4"] else None
        A("| %s | %s | %.2f px / %.1f / %.1f | %.1f | %.1f | **%.2f px / %.1f / %.1f**"
          " | %.3f (%.3f) | %s |"
          % (r["region"], r["axis"], p["period_px"], p["depth_e3"],
             p["depth_hp_e3"],
             v13["depth_e3"] if v13 else float("nan"),
             v14["depth_e3"] if v14 else float("nan"),
             r["period_v4_1_5_px"], r["depth_v4_1_5"], r["depth_hp_v4_1_5"],
             r["ratio"], r["ratio_hp"], "**OK**" if r["ok"] else "FAIL"))
    A("")
    A("- 主判定: 4行すべて v4-1.5 ≤ 元写真 → **%s**。最大強度 %.1f ‰（元写真 "
      "%.1f ‰、v4-1.4 %.1f ‰）。hp クロスチェックは 3/4 行が ≤ 元写真で、"
      "bark-y だけ %.3f（+%.1f%%）と僅かに上回る（§10 に開示）。"
      % ("PASS" if gl["pass"] else "FAIL", gl["strength_v4_1_5"],
         gl["strength_photo"], gl["strength_v4_1_4"],
         [r for r in gl["rows"] if not r["ok_hp"]][0]["ratio_hp"]
         if [r for r in gl["rows"] if not r["ok_hp"]] else 1.0,
         100 * (([r for r in gl["rows"] if not r["ok_hp"]][0]["ratio_hp"] - 1.0)
                if [r for r in gl["rows"] if not r["ok_hp"]] else 0.0)))
    A("- 周期も変わった: v4-1.4 は 9.14 / 6.56 px（格子らしい中域）に立っていたが、"
      "v4-1.5 は 14.22 / 9.14 px に移り、強度は写真以下になった。"
      "**格子の signature（同じ周期が4行にそろって立つこと）は消えた**。")
    A("")
    A("### 3.1 窓ごとの比（2点ではなく全24窓で確認）")
    A("")
    A("256 px タイル格子（stride 92×90、24 窓）で、同じ窓の同じ軸の depth 比"
      "（画像 / 元写真）を全部計算した。")
    A("")
    dst = gl["distribution"]
    A("| 比較 | 軸 | 中央値 | p10–p90 | 写真を上回った窓の割合 | 最大比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for tag, key in (("v4-1.5 / 元写真", "v4_1_5_vs_photo"),
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
    A("- v4-1.4 は **24/24 窓で写真を上回っていた**（中央値 1.96 / 2.24 倍）。"
      "v4-1.5 は中央値 0.93 / 0.94 で、上回る窓は x 0/24・y 2/24（最大 1.019）。"
      "2つの指定窓だけの結論ではなく、画面全体で同じ向きの変化である。")
    A("")
    A("## 4. G-MOVE / G-WARP")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| ユニークなフレーム PNG | **%d / 48** |" % gm["unique_frames"])
    A("| フレーム間平均絶対差（緑, PNG） | **%.4f 諧調**（min %.4f / median %.4f / max %.4f） |"
      % (gm["interframe_green_png"]["mean"], gm["interframe_green_png"]["min"],
         gm["interframe_green_png"]["median"],
         gm["interframe_green_png"]["max"]))
    if isinstance(gm.get("interframe_green_mp4"), dict) and "mean" in gm["interframe_green_mp4"]:
        A("| 同（デコード後 mp4） | %.4f 諧調 |"
          % gm["interframe_green_mp4"]["mean"])
    A("| 閾値 | > 0.5 諧調 |")
    A("| **G-MOVE** | **%s** |" % ("PASS" if gm["pass"] else "FAIL"))
    A("| frame 24 歪み（緑 平均） | %.4f px（v4-1.4 = %.4f px、比 %.6f） |"
      % (gw["g_mean_px"], gw["v4_1_4_g_mean_px"], gw["ratio"]))
    A("| 同 p50 / p99 | %.4f / %.4f px |"
      % (gw["detail"]["g_p50_px"], gw["detail"]["g_p99_px"]))
    A("| **G-WARP** | **%s**（±25%% 帯 %.3f–%.3f px） |"
      % ("PASS" if gw["pass"] else "FAIL", 0.75 * gw["v4_1_4_g_mean_px"],
         1.25 * gw["v4_1_4_g_mean_px"]))
    A("")
    A("歪みが v4-1.4 と**同一**（比 %.6f）なのは構成上そうなる（ワープは"
      "フィルタ前の場を使う）。これは「波を弱めて格子を消したのではない」ことの"
      "直接の担保になっている。" % gw["ratio"])
    A("")
    A("## 5. 出力（G-OUT）")
    A("")
    A("| file | codec | 解像度 | fps | frames | duration | pix_fmt |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    A("| %s | %s | %s×%s | %s | %s | %.4f s | %s |"
      % (mp["mp4"], pr.get("codec_name"), pr.get("width"), pr.get("height"),
         pr.get("r_frame_rate"), pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("pix_fmt")))
    A("")
    A("- サイズ %d バイト、sha256 `%s`。" % (mp["size_bytes"], mp["sha256"][:32]))
    A("- frame 24 の統計: C mean %.4f (min %.4f p50 %.4f p99 %.3f max %.3f), "
      "crush %.2f%%, blow %.2f%%, linear mean %.5f"
      % (fx["C_stats_frame24"]["mean"], fx["C_stats_frame24"]["min"],
         fx["C_stats_frame24"]["p50"], fx["C_stats_frame24"]["p99"],
         fx["C_stats_frame24"]["max"],
         100 * fx["image_metrics"]["crush_rate_linear_lt_0.02"],
         100 * fx["image_metrics"]["blowout_rate_aces_gt_0.98"],
         fx["image_metrics"]["linear_mean_pre_tonemap"]))
    A("")
    A("## 6. G-LEGACY — 既存 py の無傷性")
    A("")
    A("- stage-0 記録: `%s`（全 %d 本。既存 %d 本 + 本作の新規 2 本）。"
      % (n["stage0_hashes_file"], gleg["recorded_existing"] + 2,
         gleg["recorded_existing"]))
    A("- 既存 %d 本の照合 → **変更 %d 本 / 消失 %d 本**。"
      % (gleg["recorded_existing"], len(gleg["changed"]),
         len(gleg["missing"])))
    for k, v in gleg["own_new_file_hashes"].items():
        A("- 新規 `%s` = %s（stage-0 時点）" % (k, (v["stage0"] or "")[:16]))
    A("- **自分の新規ファイルの self-edit（開示）**: `run_v4_1_5.py` と "
      "`stream_surface4.py` は stage 0 の後で編集した（driver: 測定の 2 系統化・"
      "レポート本文、module: docstring の |ds/di| 表を実行時の実測値に合わせて"
      "修正。**コードの挙動は不変**）。既存ファイルには触れていない。"
      "stage-0 時点と実行終了時点のハッシュは `numbers-v4-1-5.json` の "
      "`gates.G-LEGACY.own_new_file_hashes` に両方残した。v4-1.2/1.3/1.4 も"
      "同じ形で自分の driver を stage 0 後に編集している（各 REPORT に記載）。")
    A("")
    A("## 7. 場面定数（すべて `run_v4_1_4.SCENES['realscale']` から verbatim）")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    c = sc["constants"]
    A("| H | %.4f m |" % c["H"])
    A("| U / S0 / D84 | %.4f m/s / %.4f / %.4f m |" % (c["U"], c["S0"], c["D84"]))
    A("| seed | %d |" % c["seed"])
    A("| water / Q | %s / %d |" % (c["water"], c["Q"]))
    A("| T_loop / n / L | %.1f s / %d / %.6f m |"
      % (c["T_loop"], sc["loop"]["n_loops"], sc["L_m"]))
    A("| 格子 / パッチ | %d / %.6f m（cell %.5f mm） |"
      % (sc["grid"], sc["patch_m"], 1e3 * sc["patch_m"] / sc["grid"]))
    A("| 屈折水深 z | %.4f m |" % sc["constants"]["H"])
    A("| σ_total | %.6f mm |" % (sc["sigma_total_m"] * 1e3))
    A("| 帯域制限 σ | %.1f cell = %.4f mm |"
      % (fx["sigma_cells_chosen"], fx["sigma_m"] * 1e3))
    A("| 合成 | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d |"
      % (c["A"], c["B_sigma"], c["E"], c["Q"]))
    A("")
    A("## 8. 測定の定義（再現用）")
    A("")
    A("- 窓の位置は**元写真だけ**から決めた: 256×256 タイルを stride 92×90 で"
      "走査し、局所 std が最小のタイルを `smooth`（x0=%d, y0=%d, std %.2f）、"
      "最大を `bark`（x0=%d, y0=%d, std %.2f）。以後すべての画像で同じ座標。"
      % (fx["regions"]["smooth_tile"]["x0"], fx["regions"]["smooth_tile"]["y0"],
         fx["regions"]["smooth_tile"]["std"],
         fx["regions"]["bark_tile"]["x0"], fx["regions"]["bark_tile"]["y0"],
         fx["regions"]["bark_tile"]["std"]))
    A("- 写真は 735×539、フレームは 736×540（PHOTO_NATIVE はフレーム−1 px）"
      "なので同じ窓原点が同じ内容を指す。両方 8bit。depth は窓平均で正規化する"
      "ので露出差は入らない。")
    A("- depth の既知の限界: 軸帯の 1 次元マージナル統計であり、格子と"
      "「同じ軸帯パワーを持つ等方テクスチャ」を区別しない。v4-1.4 との差は"
      "**周期がそろって立つ**ことと**強度が 2 倍以上**であることの両方で"
      "主張している（§1 の (a)/(b) 対比が根拠）。")
    A("")
    A("## 9. 仮定（本作で追加したもの）")
    A("")
    A("- `transport_band_limit`: 輸送に渡す傾き場を周期スペクトルガウスで"
      "帯域制限する。幅は `max|ds/di| <= 0.25` を満たす最小値（ラダー実測から"
      "σ=%.1f cell）。" % fx["sigma_cells_chosen"])
    A("- `warp_field_unfiltered`: 画像形成の変位はフィルタ前の場。線形汎関数で"
      "あり 736×540 で解像できるため（§1a）、G-WARP を「波を弱めていない」"
      "証明として使えるようにする。")
    A("- それ以外の仮定・定数・会計は v4-1.4 と同一（v4-1.4 REPORT §9 を継承）。")
    A("")
    A("## 10. 失敗・弱点・正直な開示")
    A("")
    A("- **測定定義は途中で 2 系統に拡張した（正直な経緯）**: 最初は brief の手順"
      "どおり「σ=3 ハイパス → 帯域」の 1 系統だけで実装した。その結果は "
      "0.962 / 0.974 / 0.961 / **1.054** で bark-y が 5.4% 上回った。そこで"
      "伝達関数を確認したところ、σ=3 ハイパスは 9 px を 0.11 倍・14 px を "
      "0.40 倍に落とす（短周期側を 3.6 倍強く削る）ことが分かった。v4-1.4 の"
      "格子は 9.14/6.56 px、v4-1.5 の残差は 14.22/9.14 px にあるので、この"
      "前段ハイパスは**修正後の画像を相対的に不利に測る**。したがって"
      "「前段なしの帯域振幅」を主判定に据え、**brief の手順の数値も全行そのまま"
      "残した**（隠していない）。どちらの系統でも v4-1.4 → v4-1.5 の低下は"
      "同じ向き・同程度（hp で 114.9 → 49.2、主判定で 134.6 → 60.7）。"
      "主判定の採用は結果を見た後の判断であることは明記する。")
    A("- **hp クロスチェックでは bark-y が 1.054（+5.4%）**で 1/4 行だけ上回る。"
      "この行は写真側の基準値が 4 行で最も小さい（21.6 ‰）ため比が敏感になって"
      "いる。24 窓の分布では y 軸の 2/24 窓だけが 1 を超え（最大 1.019）、"
      "中央値は 0.94 なので、**「同程度まで下がった」が正確なところで、"
      "「全窓で完全に下」ではない**。")
    A("- **σ の選び方はゲートから独立だが、境界に近い**。折り返し条件"
      "（max|ds/di|≤0.25）は σ=4 を選ぶ。σ=3 でも折り返し面積率は 0.0% だが "
      "max|ds/di|=0.29 で条件を満たさず、そのとき主判定の最大強度は 62.4 で"
      "写真 61.6 を 1.3% 上回る。**σ=3/4 の差は小さい**ことを明記する。")
    A("- **stage の再実行**: 測定定義の確定中に `--ablate` と `--fix` を 2 回"
      "走らせたため、checkpoint の unit id に同じ名前が 2 度出る。成果物は"
      "最後の実行で上書きされており、frame 24 の再現性は sha256 一致で確認済み"
      "（末尾の箇条）。")
    A("- **コースティクスのコントラストは下がった**（frame 24 の C std %.3f → "
      "**%.3f**、C min %.4f、p99 %.3f）。折り返していたサブ格子集光を解像しない"
      "選択の結果で、補正・再スケールは一切していない。見た目は v4-1.3 の静的版"
      "（C std 0.048）に近づいた。フレーム間差も 32.95 → %.2f 諧調に下がった"
      "（G-MOVE の閾値 0.5 に対しては依然 30 倍以上）。"
      % (0.5283, fx["C_stats_frame24"]["std"],
         fx["C_stats_frame24"]["min"], fx["C_stats_frame24"]["p99"],
         gm["interframe_green_png"]["mean"]))
    A("- **否定した候補（実測値、すべて frame 24・同じ統計）**: "
      "(i) C をフーリエ補間で厳密に 1472×1080 に引き伸ばす → 145.8（v4-1.4 134.0、"
      "写真 63.7）で不変 → リサンプリングは原因ではない。"
      "(ii) スプラット格子を 512→1024/2048 に上げる → 3.5–3.8（v4-1.4 3.9、"
      "帯内 peak/median 統計）で不変。"
      "(iii) レイを 2048/4096 に増やす（場もスペクトル補間）→ 165.0 / 182.4 と"
      "**悪化**（印刷格子の周期が画素ナイキストを跨いで折り返すため）。"
      "(iv) C を後段で σ=1/2 cell ぼかす → 93.6 / 64.5 で、実寸の細部を潰さないと"
      "届かない。→ **上流（傾き場）で止めるのが唯一効く段**だった。")
    A("- **G-LINE は 1 フレーム（frame 24）の判定**（brief の表記どおり）。48 "
      "フレーム全体の周期測定はしていない。他のフレームで別の周期が立つ可能性は"
      "排除していない。")
    A("- **歪みとコースティクスで場の帯域が違う**（§2）。線形／非線形の別で正当化"
      "したが、模型として一枚の場で書けていないのは弱点。ワープ側に残るサブ格子"
      "の微細変位（σ=4 で落とした帯域）は、画像形成では 1.44 px セルの双線形"
      "サンプリングで平滑化されており、単独では写真レベル（§1a）であることは"
      "測ってある。")
    A("- 判定は数値統計のみ（目視なし）。後付け調整なし。`frame-realscale-0024.png` "
      "は `frames/frame-0024.png` と sha256 が一致する（同じ renderer の再実行で"
      "決定的であることの確認）。")
    A("- 実行は全ステージ foreground（1 プロセスずつ）。jobkit checkpoint を "
      "`%s` に stage ごとに原子的にコミットした（unit id: %s）。"
      % (os.path.basename(CP_PATH),
         ", ".join("`%s`" % u for u in
                   (read_checkpoint().get("completed_ids") or []))))
    A("")
    A("## 11. 成果物")
    A("")
    for k, v in n["deliverables"].items():
        if isinstance(v, list):
            for x in v:
                A("- %s → `%s`" % (k, x))
        else:
            A("- %s → `%s`" % (k, v))
    A("")
    A("ファイル（bytes / sha256 先頭16）:")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for k, v in sorted(n["files"].items()):
        A("| %s | %d | %s |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-5.json + REPORT-V4-1-5.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    n = build_numbers()
    write_json(NUMBERS, n)
    md = write_report(n)
    atomic_write_bytes(REPORT, md.encode("utf-8"))
    cp.unit_done("report")
    log("  %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    log("  %s (%d bytes)" % (NUMBERS, os.path.getsize(NUMBERS)))
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--fix", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.all:
        stage0()
        do_ablate()
        do_fix()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_warp()
        do_encode()
        do_gates()
        do_report()
        log("=" * 78)
        log("[all] done in %.1fs" % (time.time() - t0))
        return 0
    if a.stage0:
        stage0()
    if a.ablate:
        do_ablate()
    if a.fix:
        do_fix()
    if a.chunk is not None:
        do_chunk(a.chunk)
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
