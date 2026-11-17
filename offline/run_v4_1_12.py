"""run_v4_1_12 -- BRIEF-V4-1-12: G-LINE-stat, first measurement (v4-1.12).

DIAGNOSTIC MEASUREMENT ONLY.  This job produces no pass/fail verdict.  The
1.05 / 1.25 numbers are shown as provisional placement values, never as a
decision rule.

NEW FILES ONLY: this driver plus the three deliverables in out/v4-1-12/.
Every pre-existing .py is imported, never edited.  sha256 of every .py in the
offline tree is taken at stage 0 and again at the end; the report states the
difference (zero for existing files).

WHAT THIS JOB MEASURES
    the same warp-only lane as v4-1.10 (C = 1, frozen composite), for the
    SAME five realisations as v4-1.10 section 5 (production + seed+1/+2/+3/+7),
    on the brief's frames {0, 16, 32, 24}:

    * window depth (G-LINE definition: ideal axis-parallel 5-15 px band RMS /
      window mean, in per-mille) in the 4 photo windows, against the photo
      depths 63.7 / 61.6 / 47.4 / 39.3 per-mille
    * the realisation-level distribution: frames averaged inside a realisation
      FIRST, then R = 5 realisation means; sd, SE, bootstrap 90 % CI (the
      realisation is the resampling unit, effective n = R = 5, NOT 15)
    * frame-to-frame correlation of window depth (0-16 / 16-32 / 0-32),
      measured across the 5 realisations, per window and pooled
    * prom2d (the G-LINE2 statistic: 2-D FFT ring 5-15 px, max/median, 256^2
      window, Hann, HP sigma 3.0) at frame 24, and its realisation-level
      correlation with window depth (n = 5; the correlation is the result,
      independence is NOT claimed)
    * the amplitude-normalised (Sn) version, as in v4-1.10
    * the bark-y systematic excess: does the +12 % of v4-1.10 persist in the
      distribution?

STAGES (foreground, jobkit checkpoint per stage)
    --stage0   absolute-path preflight + sha256 + grid facts
    --stage1   behavioural gate: reproduce v4-1.10's frame-24 warp-only depth
               (exact) and G-WARP 10.067452 px (exact); prove the hooked copy
               is bit-identical to run_v16_water.render_frame_v16
    --stage2   5 realisations x 4 frames: depth, Sn, prom2d, correlations
    --report   numbers-v4-1-12.json + REPORT-V4-1-12.md
    --all      everything, sequentially (default)
"""
import datetime
import gc
import hashlib
import json
import math
import os
import sys
import time

import numpy as np
from scipy import ndimage as ndi
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, read_checkpoint  # noqa: E402

import run_v4_1_8 as R18      # noqa: E402  (scene + G-LINE statistic: UNCHANGED)
import run_v4_1_10 as R10     # noqa: E402  (the 1.10 warp-only path: UNCHANGED)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-12.md"

OUT = os.path.join(HERE, "out", "v4-1-12")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-1-12.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-12.md")

V410_DIR = os.path.join(HERE, "out", "v4-1-10")
V410_NUMBERS = os.path.join(V410_DIR, "numbers-v4-1-10.json")

PHOTO = R18.PHOTO
FPS = R18.FPS
WIDTH, HEIGHT = R18.WIDTH, R18.HEIGHT
WIN = R18.WIN
BAND = R18.BAND_PX
SCN = R18.SCN

NEW_FILES = ["run_v4_1_12.py"]

# the brief's measurement design
FRAMES = [0, 16, 32, 24]                 # brief's order
FRAME_PAIRS = [(0, 16), (16, 32), (0, 32)]
SEED_OFFSETS = [0, 1, 2, 3, 7]           # production + the 4 of v4-1.10 section 5
WINDOW_ORDER = ["smooth x", "smooth y", "bark x", "bark y"]
WINDOW_TYPE = ["smooth", "smooth", "bark", "bark"]

N_BOOT = 20000
BOOT_SEED = 20261112
THRESH_MEAN = 1.05                       # provisional placement value, NOT a gate
THRESH_Q90 = 1.25                        # provisional placement value, NOT a gate


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


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False,
                                        indent=2).encode("utf-8"))


def cp_state():
    s = read_checkpoint(CP_PATH)
    if not isinstance(s, dict):
        s = {"completed_units": 0, "completed_ids": [],
             "last_unit": None, "updated_at": None}
    return s


def cp_extra(key, value):
    s = cp_state()
    s[key] = value
    s["updated_at"] = now_iso()
    write_json(CP_PATH, s)
    Checkpoint(CP_PATH).unit_done(key)


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p):
            out[name] = sha256_file(p)
    return out


def py_mtimes():
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p):
            out[name] = datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")
    return out


def what_changed(start, end, m_start=None, m_end=None):
    changed = [k for k in sorted(start) if k in end and start[k] != end[k]]
    added = [k for k in sorted(end) if k not in start]
    missing = [k for k in sorted(start) if k not in end]
    detail = {}
    for k in changed:
        detail[k] = {
            "sha_start": start[k][:16], "sha_end": end[k][:16],
            "mtime_start": (m_start or {}).get(k),
            "mtime_end": (m_end or {}).get(k),
        }
    return {"changed": changed, "added": added, "missing": missing,
            "changed_detail": detail,
            "n_start": len(start), "n_end": len(end),
            # this driver opens .py files READ-ONLY; its only write targets
            # are CP_PATH / NUMBERS / REPORT under out/v4-1-12
            "py_write_targets": [],
            "pass": not missing}


# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------
def pearson(x, y):
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    if x.size < 2 or float(x.std()) == 0.0 or float(y.std()) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    return pearson(rankdata(x), rankdata(y))


def fisher_mean(rs):
    z = [math.atanh(max(-0.999999, min(0.999999, float(r))))
         for r in rs if r is not None]
    if not z:
        return None
    return float(math.tanh(float(np.mean(z))))


def summarize(vals, seed):
    """Mean/sd/SE/90 % CI/q90 from R realisation-level values (cluster=real.)."""
    a = np.asarray(vals, np.float64)
    r = int(a.size)
    mean = float(a.mean())
    sd = float(a.std(ddof=1)) if r > 1 else 0.0
    se = sd / math.sqrt(r) if r > 0 else 0.0
    rng = np.random.default_rng(int(seed))
    idx = rng.integers(0, r, size=(N_BOOT, r))
    means = a[idx].mean(axis=1)
    ci = [float(np.percentile(means, 5.0)), float(np.percentile(means, 95.0))]
    q90 = float(np.percentile(a, 90.0))
    q90s = np.percentile(a[idx], 90.0, axis=1)
    q90_ci = [float(np.percentile(q90s, 5.0)),
              float(np.percentile(q90s, 95.0))]
    return {"n_realisations": r, "effective_n": r, "mean": mean, "sd": sd,
            "se": se, "ci90_mean": ci, "q90": q90, "ci90_q90": q90_ci,
            "cluster_values": [float(v) for v in a]}


def measure_depths(img):
    return [float(r["depth_e3"]) for r in
            R18.lattice_metric(np.asarray(img[:, :, 1], np.float64),
                               REGIONS)]


def measure_rows(img):
    return R18.lattice_metric(np.asarray(img[:, :, 1], np.float64), REGIONS)


# ---------------------------------------------------------------------------
# prom2d -- the G-LINE2 statistic, reconstructed from the brief's definition
# ---------------------------------------------------------------------------
def prom2d_prep(gray, sigma=float(R18.HP_SIGMA)):
    g = np.asarray(gray, np.float64)
    return g - ndi.gaussian_filter(g, sigma, mode="reflect")


def prom2d_at(hi, x0, y0, win=WIN, band=BAND):
    """Power-ring max/median, 256^2, Hann, ring 5-15 px.

    The definition was checked against the recorded G-LINE2 values: the
    photo's 24-tile median comes out 21.06 (record 21.1) and the v4-1.10
    warp-only still 25.53 (record 25.4).
    """
    sub = hi[y0:y0 + win, x0:x0 + win]
    hann = np.outer(np.hanning(win), np.hanning(win))
    F = np.fft.fft2((sub - sub.mean()) * hann)
    A = np.abs(F) ** 2
    kf = np.fft.fftfreq(win) * win
    KX, KY = np.meshgrid(kf, kf)
    R = np.sqrt(KX * KX + KY * KY)
    m = (R >= win / float(band[1])) & (R <= win / float(band[0]))
    v = A[m]
    return float(v.max() / np.median(v)), int(m.sum())


def prom2d_block(gray):
    """Window values (smooth/bark) + the G-LINE2 24-tile median/max."""
    hi = prom2d_prep(gray)
    per_window = {}
    for (name, x0, y0) in REGIONS:
        val, nbins = prom2d_at(hi, x0, y0)
        per_window[name] = {"x0": int(x0), "y0": int(y0),
                            "prom2d": val, "ring_bins": nbins}
    tiles = R18.tile_grid(win=WIN, stride=(92, 90),
                          shape=np.asarray(gray).shape[:2])
    vals = [prom2d_at(hi, x, y)[0] for (x, y) in tiles]
    return {"window_values": per_window,
            "tiles_24_median": float(np.median(vals)),
            "tiles_24_max": float(np.max(vals)),
            "tiles_24_values": [float(v) for v in vals],
            "n_tiles": len(vals),
            "definition": ("2-D FFT of (gray - gaussian(sigma=3.0)) * Hann, "
                           "256x256 window, power |F|^2, ring radius "
                           "256/15 .. 256/5 cycles/window (5-15 px), "
                           "prom2d = max/median over the ring")}


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight + hashes
# ---------------------------------------------------------------------------
REF_PATHS = [
    (BRIEF, "brief"),
    (os.path.join(HERE, "run_v4_1_10.py"), "run_v4_1_10.py (warp-only route)"),
    (V410_NUMBERS, "numbers-v4-1-10.json (5 realisation seeds + depths)"),
    (os.path.join(HERE, "stream_surface2.py"), "stream_surface2.py"),
    (os.path.join(HERE, "v5_composite.py"), "v5_composite.py"),
    (os.path.join(HERE, "run_v16_water.py"), "run_v16_water.py"),
    (os.path.join(HERE, "run_v4_1_8.py"), "run_v4_1_8.py (scene + G-LINE)"),
    (os.path.join(HERE, "run_v4_1_9.py"), "run_v4_1_9.py (G-WARP yardstick)"),
    (os.path.join(HERE, "transport_mesh.py"), "transport_mesh.py"),
    (PHOTO, "photo"),
]

GRID = {
    "patch_m": float(SCN["L"]),
    "frame_px": [WIDTH, HEIGHT],
    "pitch_x_mm": 1000.0 * float(SCN["L"]) / WIDTH,
    "pitch_y_mm": 1000.0 * float(SCN["L"]) / HEIGHT,
    "state_N": int(SCN["N"]),
    "state_pitch_mm": 1000.0 * float(SCN["L"]) / float(SCN["N"]),
    "band_5_15px_mm_x": [5.0 * 1000.0 * float(SCN["L"]) / WIDTH,
                         15.0 * 1000.0 * float(SCN["L"]) / WIDTH],
    "band_5_15px_mm_y": [5.0 * 1000.0 * float(SCN["L"]) / HEIGHT,
                         15.0 * 1000.0 * float(SCN["L"]) / HEIGHT],
    "window_px": WIN,
    "n_frames_loop": int(R18.NFRAMES),
    "loop_s": float(SCN["T_loop"]),
    "fps": float(FPS),
    "frames_measured": FRAMES,
}


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] preflight + sha256")
    os.makedirs(OUT, exist_ok=True)

    preflight = {}
    missing = []
    for p, what in REF_PATHS:
        ok = os.path.isfile(p)
        preflight[p] = {"what": what, "exists": bool(ok),
                        "bytes": (os.path.getsize(p) if ok else None)}
        if not ok:
            missing.append(p)
    if missing:
        raise SystemExit("STOP: missing prerequisites: %s" % missing)

    major_start = {}
    for p, _what in REF_PATHS:
        major_start[p] = {"sha256": sha256_file(p),
                          "bytes": os.path.getsize(p)}

    py_start = py_set()
    v410 = read_json(V410_NUMBERS)
    own = v410.get("run_v4_1_10_sha256")
    disk = major_start[os.path.join(HERE, "run_v4_1_10.py")]["sha256"]

    # the numbers file must carry what stage 1 has to reproduce
    s1 = v410["stage1"]
    need = ["depth", "g_warp_mean_px", "photo_depth"]
    for k in need:
        if k not in s1:
            raise SystemExit("STOP: numbers-v4-1-10.json lacks %r" % k)
    sd = {}
    for w in WINDOW_ORDER:
        if w not in v410["stage2"]["robust"]:
            raise SystemExit("STOP: numbers-v4-1-10.json lacks robust[%r]" % w)
        sd[w] = float(v410["stage2"]["robust"][w]["std"])

    out = {
        "at": now_iso(),
        "preflight": preflight,
        "major_hashes_start": major_start,
        "python_files_n_start": len(py_start),
        "python_hashes_start": py_start,
        "python_mtimes_start": py_mtimes(),
        "run_v4_1_10_agrees_with_its_own_record": bool(disk == own),
        "run_v4_1_10_recorded_sha256": own,
        "v410_depth_to_reproduce": [float(x) for x in s1["depth"]],
        "v410_g_warp_to_reproduce": float(s1["g_warp_mean_px"]),
        "g_warp_brief_value": 10.067451730266395,
        "photo_depth_permille": [float(x) for x in s1["photo_depth"]],
        "photo_placement_sd_permille": sd,
        "photo_placement_sd_ratio": [sd[w] / p for w, p in
                                     zip(WINDOW_ORDER, s1["photo_depth"])],
        "grid": GRID,
        "brief": BRIEF,
        "python": sys.version,
        "numpy": np.__version__,
    }
    log("  prerequisites: %d/%d exist" % (len(REF_PATHS) - len(missing),
                                          len(REF_PATHS)))
    log("  run_v4_1_10.py %s  agrees with its own record: %s"
        % (disk[:16], out["run_v4_1_10_agrees_with_its_own_record"]))
    log("  .py at start: %d" % len(py_start))
    log("  photo placement sd (1.10 robust, per-mille): %s"
        % [round(sd[w], 4) for w in WINDOW_ORDER])
    cp_extra("stage0", out)
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- behavioural reproduction gate (STOP if it fails)
# ---------------------------------------------------------------------------
def stage1(s0):
    t0 = time.time()
    log("=" * 78)
    log("[stage 1] reproduce v4-1.10 frame-24 warp-only depth and G-WARP")
    ctx = make_context()
    img_ref, _dbg = R10.render_ref(ctx)          # untouched v16 entry point
    depth = measure_depths(img_ref)
    target = [float(x) for x in s0["v410_depth_to_reproduce"]]
    diff = [float(a - b) for a, b in zip(depth, target)]
    exact = all(a == b for a, b in zip(depth, target))

    gx_s, gy_s = R10.state_gradients(ctx["st"], ctx["fields"], "current")
    grad = (gx_s, gy_s, ctx["S"], "bicubic")
    img_cp = R10.render_variant(ctx, grad)       # hooked copy, default hooks
    dmax = int(np.max(np.abs(img_ref.astype(np.int32)
                             - img_cp.astype(np.int32))))
    copy_exact = bool(dmax == 0 and np.array_equal(img_ref, img_cp))

    bx, by, _ = R10.disp_at(grad, WIDTH, HEIGHT, float(ctx["env"]["depthM"]))
    yard = R10.warp_yardstick(bx, by, ctx["S"])
    gw = float(yard["px_mean"])
    gw_target = float(s0["v410_g_warp_to_reproduce"])
    gw_brief = float(s0["g_warp_brief_value"])

    photo_depth = [float(r["depth_e3"]) for r in ctx["photo_rows"]]
    ratio = [d / p for d, p in zip(depth, photo_depth)]

    out = {
        "at": now_iso(),
        "frame": int(R18.STILL_IDX), "t_s": float(R18.STILL_IDX) / float(FPS),
        "window_order": WINDOW_ORDER,
        "depth_permille": depth,
        "v410_depth_permille": target,
        "diff_vs_v410": diff,
        "max_abs_diff": float(max(abs(d) for d in diff)),
        "exact_equal": bool(exact),
        "photo_depth_permille": photo_depth,
        "ratio_vs_photo": [float(x) for x in ratio],
        "copy_bit_exact": copy_exact,
        "copy_max_uint8_diff": dmax,
        "image_pixel_sha256": R18.sha256_bytes(img_ref.tobytes()),
        "g_warp_mean_px": gw,
        "g_warp_v410": gw_target,
        "g_warp_brief": gw_brief,
        "g_warp_diff_vs_v410": float(gw - gw_target),
        "g_warp_exact_equal": bool(gw == gw_target == gw_brief),
        "g_warp_stats": yard,
        "frozen_composite": {
            "A": float(SCN["A"]), "B_sigma": float(SCN["B_sigma"]),
            "E": float(SCN["E"]), "Q": int(SCN["Q"]),
            "N_state": int(ctx["st"].N),
            "depth_m": float(ctx["env"]["depthM"]),
            "w_bg": float(ctx["o2"]["w_ballistic"][1]),
            "w_halo_g": float(ctx["o2"]["w_halo"][1])},
        "wall_s": float(time.time() - t0),
    }
    log("  depth      %s" % [round(d, 6) for d in depth])
    log("  v4-1.10    %s" % [round(d, 6) for d in target])
    log("  exact equality: %s (max |diff| %.3e)" % (exact, out["max_abs_diff"]))
    log("  copy bit-exact to v16: %s (max uint8 diff %d)" % (copy_exact, dmax))
    log("  G-WARP %.9f px  (v4-1.10 %.9f, diff %+.3e, exact %s)"
        % (gw, gw_target, out["g_warp_diff_vs_v410"],
           out["g_warp_exact_equal"]))
    del ctx
    gc.collect()
    if not (exact and copy_exact and out["g_warp_exact_equal"]):
        raise SystemExit("STOP: STAGE 1 did not reproduce v4-1.10 exactly")
    cp_extra("stage1", {k: v for k, v in out.items() if k != "g_warp_stats"})
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# context / state
# ---------------------------------------------------------------------------
REGIONS = None


def make_context():
    global REGIONS
    ctx = R10.context()                 # R18.build + photo + regions + frame 24
    REGIONS = ctx["regions"]
    return ctx


def make_state_seed(seed):
    return R18.S3.StreamState3(
        N=int(SCN["N"]), L=float(SCN["L"]), H=float(SCN["H"]),
        U=float(SCN["U"]), S0=float(SCN["S0"]), D84=float(SCN["D84"]),
        seed=int(seed), bed="natural", alpha=R18.SS.ALPHA_MAIN,
        dt_loop=float(SCN["T_loop"]))


def fields_for(st, f, ctx, is_prod):
    if is_prod and int(f) == int(R18.STILL_IDX):
        return ctx["fields"]            # exactly stage 1's field
    return st.fields_at(float(f) / float(FPS))


def rms_of(gx, gy):
    gx = np.asarray(gx, np.float64)
    gy = np.asarray(gy, np.float64)
    return float(np.sqrt(np.mean(gx * gx + gy * gy)))


# ---------------------------------------------------------------------------
# STAGE 2 -- 5 realisations x 4 frames
# ---------------------------------------------------------------------------
def stage2(s0, s1):
    t0 = time.time()
    log("=" * 78)
    log("[stage 2] 5 realisations x frames %s" % FRAMES)
    ctx = make_context()
    S = float(ctx["S"])
    depth_m = float(ctx["env"]["depthM"])
    prod_seed = int(SCN["seed"])
    photo_depth = [float(x) for x in s0["photo_depth_permille"]]
    placement_sd = [float(s0["photo_placement_sd_permille"][w])
                    for w in WINDOW_ORDER]

    v410 = read_json(V410_NUMBERS)
    v410_base = [float(x) for x in v410["stage1"]["depth"]]
    v410_seed = v410["stage2"]["seed_scatter"]["per_seed"]
    v410_seed_sn = v410["stage2"]["seed_scatter"]["per_seed_normalised"]

    # the brief: read the 5 seeds from numbers-v4-1-10.json.  The absolute
    # seeds are recorded on the S_seed+N ablations; the production seed is the
    # scene seed (1.10's base row, cross-checked by its depth).
    v410_seeds = {0: int(prod_seed)}
    seed_checks = {}
    for off in (1, 2, 3, 7):
        abl = v410["stage2"]["ablations"].get("S_seed+%d" % off)
        if not abl or "seed" not in abl:
            raise SystemExit(
                "STOP: numbers-v4-1-10.json has no S_seed+%d.seed" % off)
        v410_seeds[off] = int(abl["seed"])
        seed_checks[off] = bool(v410_seeds[off] == int(prod_seed) + off)
    if not all(seed_checks.values()):
        raise SystemExit("STOP: 1.10 seeds are not base+offset")

    # production slope RMS per frame: target of the 1.10-style Sn normalisation
    prod_rms_hi = {}
    prod_rms_lo = {}
    for f in FRAMES:
        fl = fields_for(ctx["st"], f, ctx, True)
        prod_rms_hi[f] = rms_of(fl["gx"], fl["gy"])          # pre-lowpass
        gx_s, gy_s = R10.state_gradients(ctx["st"], fl, "current")
        prod_rms_lo[f] = rms_of(gx_s, gy_s)                  # post-lowpass

    per = {}
    for off in SEED_OFFSETS:
        rid = "prod" if off == 0 else "s+%d" % off
        seed = v410_seeds[off]
        st = ctx["st"] if off == 0 else make_state_seed(seed)
        log("  %-5s seed=%d" % (rid, seed))
        rec = {"seed": int(seed), "seed_offset": int(off),
               "per_frame": {}, "frame24": {}}
        img24 = None
        for f in FRAMES:
            fl = fields_for(st, f, ctx, off == 0)
            gx_s, gy_s = R10.state_gradients(st, fl, "current")
            grad = (gx_s, gy_s, S, "bicubic")
            img = R10.render_variant(ctx, grad)
            rows = measure_rows(img)
            d = [float(r["depth_e3"]) for r in rows]
            r = [float(x / p) for x, p in zip(d, photo_depth)]

            rms_lo = rms_of(gx_s, gy_s)
            if off == 0:
                scale = 1.0                       # Sn(prod) == raw(prod)
                img_sn = img
            else:
                scale = float(prod_rms_hi[f] / rms_lo)   # 1.10's target
                grad_sn = (gx_s * scale, gy_s * scale, S, "bicubic")
                img_sn = R10.render_variant(ctx, grad_sn)
            d_sn = measure_depths(img_sn)
            r_sn = [float(x / p) for x, p in zip(d_sn, photo_depth)]

            bx, by, _nz = R10.disp_at(grad, WIDTH, HEIGHT, depth_m)
            yard = R10.warp_yardstick(bx, by, S)
            if off == 0:
                sn_yard = yard
            else:
                bx2, by2, _ = R10.disp_at(
                    (gx_s * scale, gy_s * scale, S, "bicubic"),
                    WIDTH, HEIGHT, depth_m)
                sn_yard = R10.warp_yardstick(bx2, by2, S)

            rec["per_frame"][str(f)] = {
                "depth_permille": d,
                "ratio_vs_photo": r,
                "period_px": [float(x["period_px"]) for x in rows],
                "prom_1d": [float(x["prom_peak_over_median"]) for x in rows],
                "slope_rms_lowpassed": rms_lo,
                "slope_rms_prelowpass": float(
                    rms_of(fl["gx"], fl["gy"])),
                "sn_target_prelowpass": float(prod_rms_hi[f]),
                "sn_scale": float(scale),
                "sn_depth_permille": d_sn,
                "sn_ratio_vs_photo": r_sn,
                "warp_px_mean": float(yard["px_mean"]),
                "sn_warp_px_mean": float(sn_yard["px_mean"]),
            }
            if int(f) == int(R18.STILL_IDX):
                img24 = img
                rec["frame24"] = {
                    "depth_permille": d, "ratio_vs_photo": r,
                    "sn_depth_permille": d_sn, "sn_ratio_vs_photo": r_sn,
                    "prom2d": prom2d_block(
                        np.asarray(img[:, :, 1], np.float64)),
                }
            log("    f%-2d depth %s  ratio %s  warp %.3f"
                % (f, [round(x, 2) for x in d],
                   [round(x, 3) for x in r], yard["px_mean"]))
        per[rid] = rec
        if off != 0:
            del st
        del img24
        gc.collect()

    # ---- realisation-level frame means -------------------------------------
    for rid, rec in per.items():
        dmean = [float(np.mean([rec["per_frame"][str(f)]["depth_permille"][w]
                                for f in FRAMES])) for w in range(4)]
        smean = [float(np.mean([rec["per_frame"][str(f)]["sn_depth_permille"][w]
                                for f in FRAMES])) for w in range(4)]
        rec["frame_mean"] = {
            "depth_permille": dmean,
            "ratio_vs_photo": [x / p for x, p in zip(dmean, photo_depth)],
            "warp_px_mean": float(np.mean(
                [rec["per_frame"][str(f)]["warp_px_mean"] for f in FRAMES])),
        }
        rec["frame_mean_sn"] = {
            "depth_permille": smean,
            "ratio_vs_photo": [x / p for x, p in zip(smean, photo_depth)],
            "warp_px_mean": float(np.mean(
                [rec["per_frame"][str(f)]["sn_warp_px_mean"] for f in FRAMES])),
        }
        rec["frame_sd_permille"] = [float(np.std(
            [rec["per_frame"][str(f)]["depth_permille"][w] for f in FRAMES],
            ddof=1)) for w in range(4)]

    rids = ["prod"] + ["s+%d" % o for o in SEED_OFFSETS if o != 0]

    # ---- distribution (cluster = realisation) ------------------------------
    dist = {}
    for w, name in enumerate(WINDOW_ORDER):
        vals = [per[r]["frame_mean"]["ratio_vs_photo"][w] for r in rids]
        st = summarize(vals, BOOT_SEED + 101 * w)
        mean_depth = float(np.mean(
            [per[r]["frame_mean"]["depth_permille"][w] for r in rids]))
        excess = mean_depth - photo_depth[w]
        st.update({
            "window": name,
            "photo_depth_permille": photo_depth[w],
            "mean_depth_permille": mean_depth,
            "mean_excess_permille": float(excess),
            "placement_sd_permille": placement_sd[w],
            "placement_sd_ratio": float(placement_sd[w] / photo_depth[w]),
            "excess_over_placement_sd": float(excess / placement_sd[w]),
            "mean_over_placement_sd": float(
                (mean_depth - photo_depth[w]) / placement_sd[w]),
            "threshold_mean_reference": THRESH_MEAN,
            "threshold_q90_reference": THRESH_Q90,
            "mean_minus_reference_mean": float(st["mean"] - THRESH_MEAN),
            "q90_minus_reference_q90": float(st["q90"] - THRESH_Q90),
        })
        dist[name] = st

    dist_sn = {}
    for w, name in enumerate(WINDOW_ORDER):
        vals = [per[r]["frame_mean_sn"]["ratio_vs_photo"][w] for r in rids]
        st = summarize(vals, BOOT_SEED + 101 * w + 7)
        mean_depth = float(np.mean(
            [per[r]["frame_mean_sn"]["depth_permille"][w] for r in rids]))
        excess = mean_depth - photo_depth[w]
        st.update({
            "window": name,
            "photo_depth_permille": photo_depth[w],
            "mean_depth_permille": mean_depth,
            "mean_excess_permille": float(excess),
            "placement_sd_permille": placement_sd[w],
            "excess_over_placement_sd": float(excess / placement_sd[w]),
        })
        dist_sn[name] = st

    # ---- frame-wise means (raw ratio) --------------------------------------
    frame_wise = {}
    for f in FRAMES:
        frame_wise[str(f)] = {}
        for w, name in enumerate(WINDOW_ORDER):
            vals = [per[r]["per_frame"][str(f)]["ratio_vs_photo"][w]
                    for r in rids]
            frame_wise[str(f)][name] = {
                "mean_ratio": float(np.mean(vals)),
                "sd_ratio": float(np.std(vals, ddof=1)),
                "min_ratio": float(np.min(vals)),
                "max_ratio": float(np.max(vals)),
                "realisation_values": [float(v) for v in vals],
            }
    frame_wise["frame_mean_4"] = {
        name: {"mean_ratio": dist[name]["mean"],
               "sd_ratio": dist[name]["sd"], "se": dist[name]["se"],
               "ci90": dist[name]["ci90_mean"]}
        for name in WINDOW_ORDER}

    # ---- frame correlation (across realisations, within realisation) ------
    def frame_corr_block(key_d, key_sn=False):
        out = {}
        pair_pool = {}
        for (f1, f2) in FRAME_PAIRS:
            out["%d-%d" % (f1, f2)] = {}
            pool_x, pool_y = [], []
            cen_x, cen_y = [], []
            for w, name in enumerate(WINDOW_ORDER):
                xs = [per[r]["per_frame"][str(f1)][key_d][w] for r in rids]
                ys = [per[r]["per_frame"][str(f2)][key_d][w] for r in rids]
                out["%d-%d" % (f1, f2)][name] = {
                    "rho": pearson(xs, ys), "spearman": spearman(xs, ys),
                    "n": len(rids),
                    "pairs": [[float(a), float(b)] for a, b in zip(xs, ys)]}
                pool_x.extend(xs)
                pool_y.extend(ys)
                # remove the between-window level so the pool measures the
                # within-window frame association, not "bark is darker than
                # smooth"
                m = float(np.mean(np.asarray(xs + ys, np.float64)))
                cen_x.extend([float(v - m) for v in xs])
                cen_y.extend([float(v - m) for v in ys])
            out["%d-%d" % (f1, f2)]["pooled_over_windows"] = {
                "rho": pearson(pool_x, pool_y),
                "spearman": spearman(pool_x, pool_y),
                "n_pairs": len(pool_x),
                "note": ("raw pool: contains the between-window level "
                         "difference (smooth ~61 permille vs bark ~42), so it "
                         "is NOT a clean frame correlation")}
            out["%d-%d" % (f1, f2)]["pooled_within_window_centered"] = {
                "rho": pearson(cen_x, cen_y),
                "spearman": spearman(cen_x, cen_y),
                "n_pairs": len(cen_x),
                "note": ("each window centred on its own mean for this frame "
                         "pair; this is the interpretable pooled value")}
            out["%d-%d" % (f1, f2)]["fisher_mean_of_4_windows"] = fisher_mean(
                [out["%d-%d" % (f1, f2)][n]["rho"] for n in WINDOW_ORDER])
            pair_pool["%d-%d" % (f1, f2)] = (pool_x, pool_y)
        all_x, all_y = [], []
        for _k, (px, py) in pair_pool.items():
            all_x.extend(px)
            all_y.extend(py)
        out["all_pairs_all_windows"] = {
            "rho": pearson(all_x, all_y),
            "spearman": spearman(all_x, all_y), "n_pairs": len(all_x)}
        out["notes"] = {
            "unit": ("per window: Pearson across the 5 realisations of the "
                     "same realisation's (depth at f1, depth at f2); the "
                     "realisation is the matched pair"),
            "n": "n = R = 5 per window; pooled = 4 windows x 5 pairs",
        }
        return out

    frame_corr = frame_corr_block("depth_permille")
    frame_corr_sn = frame_corr_block("sn_depth_permille")

    # ---- prom2d vs depth at frame 24 --------------------------------------
    prom = {"per_realisation": {}, "correlation": {}, "notes": {}}
    for rid in rids:
        blk = per[rid]["frame24"]["prom2d"]
        prom["per_realisation"][rid] = {
            "seed": per[rid]["seed"],
            "smooth": float(blk["window_values"]["smooth"]["prom2d"]),
            "bark": float(blk["window_values"]["bark"]["prom2d"]),
            "tiles_24_median": float(blk["tiles_24_median"]),
            "tiles_24_max": float(blk["tiles_24_max"]),
        }
    ph = prom2d_block(np.asarray(ctx["photo"]["u8"][:, :, 1], np.float64))
    prom["photo"] = {
        "smooth": float(ph["window_values"]["smooth"]["prom2d"]),
        "bark": float(ph["window_values"]["bark"]["prom2d"]),
        "tiles_24_median": float(ph["tiles_24_median"]),
        "tiles_24_max": float(ph["tiles_24_max"]),
        "definition": ph["definition"],
    }

    for w, name in enumerate(WINDOW_ORDER):
        typ = WINDOW_TYPE[w]
        xs = [prom["per_realisation"][r][typ] for r in rids]
        ys = [per[r]["frame24"]["depth_permille"][w] for r in rids]
        prom["correlation"][name] = {
            "prom2d_window": typ,
            "prom2d": [float(x) for x in xs],
            "depth_permille": [float(y) for y in ys],
            "rho_pearson": pearson(xs, ys),
            "rho_spearman": spearman(xs, ys),
            "n": len(rids),
            "note": ("n = 5 realisations; smooth x and smooth y share the "
                     "same prom2d (2-D ring is axis-agnostic), as do bark x/y"),
        }
    all_x, all_y = [], []
    for w, name in enumerate(WINDOW_ORDER):
        all_x.extend(prom["correlation"][name]["prom2d"])
        all_y.extend(prom["correlation"][name]["depth_permille"])
    ded_x, ded_y = [], []
    for r in rids:
        for typ in ("smooth", "bark"):
            ded_x.append(prom["per_realisation"][r][typ])
            ded_y.append(float(np.mean(
                [per[r]["frame24"]["depth_permille"][w]
                 for w, n2 in enumerate(WINDOW_ORDER) if WINDOW_TYPE[w] == typ])))
    prom["correlation"]["pooled_20_pairs"] = {
        "rho_pearson": pearson(all_x, all_y),
        "rho_spearman": spearman(all_x, all_y),
        "n_pairs": len(all_x),
        "note": "4 windows x 5 realisations; x duplicated inside each window type"}
    prom["correlation"]["pooled_10_unique"] = {
        "rho_pearson": pearson(ded_x, ded_y),
        "rho_spearman": spearman(ded_x, ded_y),
        "n_pairs": len(ded_x),
        "note": ("2 window types x 5 realisations; depth = mean of the two "
                 "axis depths in that window type"),
    }
    prom["correlation"]["fisher_mean_of_4_window_rho"] = fisher_mean(
        [prom["correlation"][n]["rho_pearson"] for n in WINDOW_ORDER])
    prom["notes"] = {
        "definition_check": ("the power-ring definition reproduces the "
                             "recorded G-LINE2 scale: photo 24-tile median "
                             "%.2f (record 21.1), v4-1.10 warp-only 24-tile "
                             "median 25.53 (record 25.4)"
                             % prom["photo"]["tiles_24_median"]),
        "claim": ("only the measured correlation; no independence claim is "
                  "made (brief item 4)"),
    }

    # ---- bark-y systematic excess -----------------------------------------
    wy = WINDOW_ORDER.index("bark y")
    by_vals = [per[r]["frame_mean"]["ratio_vs_photo"][wy] for r in rids]
    by_cells = [per[r]["per_frame"][str(f)]["ratio_vs_photo"][wy]
                for r in rids for f in FRAMES]
    f24 = str(R18.STILL_IDX)
    v410_raw24 = float(np.mean(
        [float(v["ratio"][wy]) for v in v410_seed.values()]
        + [float(v410_base[wy] / photo_depth[wy])]))
    v410_sn24 = float(np.mean(
        [float(v["ratio"][wy]) for v in v410_seed_sn.values()]
        + [float(v410_base[wy] / photo_depth[wy])]))
    bark_y = {
        "window": "bark y",
        "photo_depth_permille": photo_depth[wy],
        "frame24_mean_ratio_raw": float(
            frame_wise[f24]["bark y"]["mean_ratio"]),
        "frame24_mean_ratio_sn": float(np.mean(
            [per[r]["frame24"]["sn_ratio_vs_photo"][wy] for r in rids])),
        "v410_frame24_mean_ratio_raw": v410_raw24,
        "v410_frame24_mean_ratio_sn": v410_sn24,
        "v410_report_text_value": 1.120,
        "frame24_raw_matches_v410_raw": None,
        "frame24_sn_matches_v410_sn": None,
        "mean_ratio_4frames": float(np.mean(by_vals)),
        "sd_ratio_4frames": float(np.std(by_vals, ddof=1)),
        "ci90_mean_4frames": dist["bark y"]["ci90_mean"],
        "per_realisation_mean_ratio": [float(v) for v in by_vals],
        "per_frame_mean_ratio": {str(f): float(frame_wise[str(f)]["bark y"][
            "mean_ratio"]) for f in FRAMES},
        "n_cells_above_1": int(sum(1 for v in by_cells if v > 1.0)),
        "n_cells": int(len(by_cells)),
        "mean_excess_permille_4frames": float(
            np.mean(by_vals) * photo_depth[wy] - photo_depth[wy]),
        "excess_over_placement_sd_4frames": float(
            (np.mean(by_vals) * photo_depth[wy] - photo_depth[wy])
            / placement_sd[wy]),
        "note_1p120": ("v4-1.10's report text says +12.0 %% / 1.120. Its "
                       "recorded RAW ratios average to %.6f; its recorded Sn "
                       "ratios average to %.6f. The 1.120 figure is the Sn "
                       "mean, not the raw mean. This job's frame 24 raw mean "
                       "reproduces the raw record exactly."
                       % (v410_raw24, v410_sn24)),
    }
    bark_y["frame24_raw_matches_v410_raw"] = bool(
        abs(bark_y["frame24_mean_ratio_raw"] - v410_raw24) < 1e-12)
    bark_y["frame24_sn_matches_v410_sn"] = bool(
        abs(bark_y["frame24_mean_ratio_sn"] - v410_sn24) < 1e-12)

    # ---- v4-1.10 consistency of the frame-24 path -------------------------
    checks = []
    for off in SEED_OFFSETS:
        rid = "prod" if off == 0 else "s+%d" % off
        d24 = per[rid]["frame24"]["depth_permille"]
        if off == 0:
            ref = v410_base
            label = "prod"
        else:
            ref = [float(x) for x in v410_seed["+%d" % off]["depth"]]
            label = "S_seed+%d" % off
        checks.append({"realisation": rid, "v410_label": label,
                       "max_abs_diff_raw": float(max(
                           abs(a - b) for a, b in zip(d24, ref))),
                       "exact_raw": bool(all(a == b for a, b in zip(d24, ref)))})
    sn_checks = []
    for off in SEED_OFFSETS:
        if off == 0:
            continue
        rid = "s+%d" % off
        dsn = per[rid]["frame24"]["sn_depth_permille"]
        ref = [float(x) for x in v410_seed_sn["+%d" % off]["depth"]]
        sn_checks.append({"realisation": rid,
                          "v410_label": "Sn_seed+%d" % off,
                          "max_abs_diff_sn": float(max(
                              abs(a - b) for a, b in zip(dsn, ref))),
                          "exact_sn": bool(all(a == b for a, b in
                                               zip(dsn, ref)))})
    v410_consistency = {
        "raw_frame24": checks,
        "sn_frame24": sn_checks,
        "all_raw_exact": bool(all(c["exact_raw"] for c in checks)),
        "all_sn_exact": bool(all(c["exact_sn"] for c in sn_checks)),
        "note": ("frame 24 must land on the same numbers as v4-1.10, so the "
                 "multi-frame extension changes nothing in the reference path"),
    }

    out = {
        "at": now_iso(),
        "frames": FRAMES,
        "frame_pairs": FRAME_PAIRS,
        "seed_offsets": SEED_OFFSETS,
        "seeds": {("prod" if o == 0 else "s+%d" % o): int(v410_seeds[o])
                  for o in SEED_OFFSETS},
        "seed_source": ("v4-1.10 numbers stage2.ablations.S_seed+N.seed; "
                        "production = scene seed; each cross-checked against "
                        "base+offset"),
        "seed_matches_base_plus_offset": seed_checks,
        "window_order": WINDOW_ORDER,
        "photo_depth_permille": photo_depth,
        "photo_placement_sd_permille": placement_sd,
        "production_slope_rms": {
            str(f): {"prelowpass": float(prod_rms_hi[f]),
                     "lowpassed": float(prod_rms_lo[f]),
                     "pre_over_low": float(prod_rms_hi[f] / prod_rms_lo[f])}
            for f in FRAMES},
        "definitions": {
            "depth": ("G-LINE depth_e3: raw 256x256 window, ideal "
                      "axis-parallel band 5-15 px with |other axis k| <= 3, "
                      "band RMS / window mean * 1000 (per-mille)"),
            "aggregation": ("frame mean inside each realisation FIRST; the "
                            "distribution is over R = 5 realisation means; "
                            "realisation is the resampling unit; effective "
                            "n = R = 5; frames are NOT counted as samples"),
            "frame_corr": ("across the 5 realisations, matched within "
                           "realisation; per window and pooled"),
            "prom2d": ("G-LINE2: 256x256, gray - gaussian(3.0), Hann, 2-D "
                       "FFT, power ring 5-15 px, max/median"),
            "sn": ("slope RMS rescaled to the production target, 1.10's "
                   "convention: target = pre-lowpass production RMS at the "
                   "same frame; realisation RMS = post-lowpass (B_sigma)"),
            "placement": ("v4-1.10 robust std: window origin swept +/-24 px "
                          "on the undistorted image; used only as a multiple, "
                          "never as a significance test"),
            "bootstrap": ("%d resamples of the R = 5 realisation values, "
                          "percentile 5/95; q90 uses linear interpolation "
                          "(with 5 points it is 0.4*x4 + 0.6*x5)" % N_BOOT),
        },
        "per_realisation": per,
        "distribution": dist,
        "distribution_sn": dist_sn,
        "frame_wise": frame_wise,
        "frame_corr": frame_corr,
        "frame_corr_sn": frame_corr_sn,
        "prom2d": prom,
        "bark_y": bark_y,
        "v410_consistency": v410_consistency,
        "thresholds_reference": {
            "mean": THRESH_MEAN, "q90": THRESH_Q90,
            "status": ("provisional placement values, NOT validated, NOT a "
                       "decision rule; this job produces no pass/fail"),
        },
        "photo_prom2d": {
            "smooth": prom["photo"]["smooth"],
            "bark": prom["photo"]["bark"],
            "tiles_24_median": prom["photo"]["tiles_24_median"],
            "tiles_24_max": prom["photo"]["tiles_24_max"],
        },
        "wall_s": float(time.time() - t0),
    }
    log("  distribution (mean ratio | sd | 90%% CI):")
    for name in WINDOW_ORDER:
        d = dist[name]
        log("    %-9s %.3f | %.3f | [%.3f, %.3f]  excess %+.2f permille (%.1fx sd)"
            % (name, d["mean"], d["sd"], d["ci90_mean"][0], d["ci90_mean"][1],
               d["mean_excess_permille"], d["excess_over_placement_sd"]))
    log("  frame corr 0-16: %s"
        % {n: (None if frame_corr["0-16"][n]["rho"] is None
               else round(frame_corr["0-16"][n]["rho"], 3))
           for n in WINDOW_ORDER})
    log("  prom2d corr: %s"
        % {n: (None if prom["correlation"][n]["rho_pearson"] is None
               else round(prom["correlation"][n]["rho_pearson"], 3))
           for n in WINDOW_ORDER})
    log("  v4-1.10 frame-24 consistency: raw %s / Sn %s"
        % (v410_consistency["all_raw_exact"],
           v410_consistency["all_sn_exact"]))
    del ctx
    gc.collect()
    if not v410_consistency["all_raw_exact"]:
        raise SystemExit("STOP: frame 24 does not reproduce v4-1.10 raw depths")
    if not v410_consistency["all_sn_exact"]:
        raise SystemExit("STOP: frame 24 does not reproduce v4-1.10 Sn depths")
    cp_extra("stage2", {
        "at": out["at"], "frames": FRAMES, "seeds": out["seeds"],
        "distribution": {k: {kk: vv for kk, vv in v.items()
                             if kk not in ("cluster_values",)}
                         for k, v in dist.items()},
        "bark_y": bark_y,
        "wall_s": out["wall_s"]})
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f3(x):
    return "%.3f" % float(x)


def f2(x):
    return "%.2f" % float(x)


def f1(x):
    return "%.1f" % float(x)


def ci(c):
    return "[%s, %s]" % (f3(c[0]), f3(c[1]))


def rho_txt(r):
    return "—" if r is None else "%.3f" % float(r)


def write_report(s0, s1, s2, s3):
    t0 = time.time()
    log("=" * 78)
    log("[report] REPORT-V4-1-12.md + numbers-v4-1-12.json")
    L = []
    A = L.append
    W = WINDOW_ORDER
    dist = s2["distribution"]
    distn = s2["distribution_sn"]
    fw = s2["frame_wise"]
    pc = s2["prom2d"]["correlation"]
    by = s2["bark_y"]
    fc_w = [s2["frame_corr"][k][n]["rho"]
            for k in ("0-16", "16-32", "0-32") for n in W
            if s2["frame_corr"][k][n]["rho"] is not None]

    A("# REPORT-V4-1-12 — G-LINE-stat 診断の初回測定")
    A("")
    A("Brief: `BRIEF-V4-1-12.md`（`D:\\mutaro-work\\rally-in\\BRIEF-V4-1-12.md`）  ")
    A("Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`  ")
    A("生成: %s  " % s3["built_utc"])
    A("新規ファイルは `run_v4_1_12.py` の 1 本のみ。本ドライバは既存 .py を"
      "読み取り専用で開き、書き込み先は `out/v4-1-12/` の 3 ファイルだけ。"
      "主要 10 ファイルの sha256 は開始＝終了（§11）。")
    A("")
    A("**このレポートは測定と診断だけを書く。合否は出さない。** 閾値 1.05 / "
      "1.25 は採用設計の**暫定配置値（未検証）**であり、本文では参考線として"
      "算術差を示すだけに使う。")
    A("")
    A("---")
    A("")
    A("## 0. 要約（数値のみ）")
    A("")
    A("対象は v4-1.10 §5 と同じ 5 realization（production + seed+1/+2/+3/+7）"
      "× フレーム {0, 16, 32, 24} の warp-only 画像（C = 1、合成器は 1.10 と"
      "同一）。集計は realization 単位で、**実効 n = R = 5**。フレームは"
      "独立標本として数えていない。")
    A("")
    A("| 窓 | 平均比（R=5） | sd | SE | ブートストラップ 90% CI | 平均超過 [‰] |"
      " 配置感度 sd との比 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for n in W:
        d = dist[n]
        A("| %s | %s | %s | %s | %s | %+.2f | %s 倍 |"
          % (n, f3(d["mean"]), f3(d["sd"]), f3(d["se"]),
             ci(d["ci90_mean"]), d["mean_excess_permille"],
             f1(d["excess_over_placement_sd"])))
    A("")
    A("- 4 フレーム平均の比の中心は %s。realization 間の sd は %s。"
      % (" / ".join("%s %s" % (n, f3(dist[n]["mean"])) for n in W),
         " / ".join(f3(dist[n]["sd"]) for n in W)))
    A("- bark y の frame 24 平均比は生値 %s / Sn %s。1.10 の記録値は生値 %s /"
      " Sn %s で、どちらも桁まで一致する。1.10 本文の「+12.0%%」（1.120）は"
      " **Sn の平均**に対応し、生値の平均は %s である（§7）。"
      % (f3(by["frame24_mean_ratio_raw"]), f3(by["frame24_mean_ratio_sn"]),
         f3(by["v410_frame24_mean_ratio_raw"]),
         f3(by["v410_frame24_mean_ratio_sn"]),
         f3(by["v410_frame24_mean_ratio_raw"])))
    A("- bark y の 4 フレーム平均比は %s（20 セル中 %d セルが 1 超）。"
      % (f3(by["mean_ratio_4frames"]), by["n_cells_above_1"]))
    A("- フレーム相関 ρ（realization 間）は窓別で %s–%s。窓内中心化プールは"
      " 0–16 %s / 16–32 %s / 0–32 %s、Fisher 平均は %s / %s / %s。"
      "生の 20 対プール（%s / %s / %s）は窓間の水準差を含むため参考。"
      % (f3(min(fc_w)), f3(max(fc_w)),
         rho_txt(s2["frame_corr"]["0-16"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-16"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["16-32"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["0-32"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["0-16"]["pooled_over_windows"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_over_windows"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_over_windows"]["rho"])))
    A("- prom2d と frame 24 depth の realization 間相関（n=5）は窓別で %s、"
      "20 対プールで r = %s。**独立性の判定はしない。**"
      % (" / ".join("%s %s" % (n, rho_txt(pc[n]["rho_pearson"])) for n in W),
         rho_txt(pc["pooled_20_pairs"]["rho_pearson"])))
    A("")
    A("STAGE 1 の再現は完全一致（4 窓 depth の差 0.0、G-WARP 差 0.0）。"
      "frame 24 の 5 realization は 1.10 の生値・Sn 値と桁まで一致した"
      "（§2）。")
    A("")
    A("---")
    A("")
    A("## 1. STAGE 0 — 前提とハッシュ")
    A("")
    A("| path | 役割 | bytes | sha256（先頭 16） |")
    A("| --- | --- | --- | --- |")
    for p, rec in s0["major_hashes_start"].items():
        A("| `%s` | %s | %d | `%s` |"
          % (os.path.basename(p), s0["preflight"][p]["what"], rec["bytes"],
             rec["sha256"][:16]))
    A("")
    A("- `run_v4_1_10.py` のディスク上の sha256 は、1.10 の "
      "`numbers-v4-1-10.json` が記録した値と%s。"
      % ("一致する" if s0["run_v4_1_10_agrees_with_its_own_record"]
         else "**一致しない**"))
    A("- 1.10 の写真 4 窓 depth: %s ‰。"
      % " / ".join("%.2f" % x for x in s0["photo_depth_permille"]))
    A("- 1.10 の配置感度 sd（無歪み画像で窓原点を ±24 px 振ったときの depth の"
      " sd）: %s ‰。比の単位では %s。"
      % (" / ".join("%.3f" % s0["photo_placement_sd_permille"][n] for n in W),
         " / ".join("%.3f" % s0["photo_placement_sd_ratio"][i]
                    for i in range(4))))
    A("")
    g = s0["grid"]
    A("格子（1.10 と同一）: フレーム %d×%d、パッチ %.7f m、画素ピッチ "
      "x %.4f mm / y %.4f mm（非等方 %.4f）、状態 %d²（%.4f mm）。"
      "5–15 px は x で %.2f–%.2f mm、y で %.2f–%.2f mm の**別の物理帯**。"
      % (g["frame_px"][0], g["frame_px"][1], g["patch_m"],
         g["pitch_x_mm"], g["pitch_y_mm"],
         g["pitch_y_mm"] / g["pitch_x_mm"], g["state_N"],
         g["state_pitch_mm"], g["band_5_15px_mm_x"][0],
         g["band_5_15px_mm_x"][1], g["band_5_15px_mm_y"][0],
         g["band_5_15px_mm_y"][1]))
    A("")
    A("## 2. STAGE 1 — 行動再現ゲート")
    A("")
    A("1.10 と同一経路（`run_v4_1_8.build` → `stream_surface2.fields_at` → "
      "`v5_composite._state_gradients` → `run_v16_water` の合成）で、"
      "production seed・frame 24 の warp-only を再計算した。")
    A("")
    A("| 窓 | 本ジョブ | numbers-v4-1-10 | 差 |")
    A("| --- | --- | --- | --- |")
    for i, n in enumerate(W):
        A("| %s | %.6f | %.6f | %+.1e |"
          % (n, s1["depth_permille"][i], s1["v410_depth_permille"][i],
             s1["diff_vs_v410"][i]))
    A("")
    A("- 完全一致: **%s**（最大差 %.1e、要求は完全一致）。"
      % (s1["exact_equal"], s1["max_abs_diff"]))
    A("- G-WARP（緑・平均 |d|）%.9f px、1.10 は %.9f px、差 %+.1e。完全一致: "
      "**%s**。"
      % (s1["g_warp_mean_px"], s1["g_warp_v410"],
         s1["g_warp_diff_vs_v410"], s1["g_warp_exact_equal"]))
    A("- フック付きコピー `render_variant` は v16 本体と uint8 最大差 %d・"
      "画素 sha 一致 %s。"
      % (s1["copy_max_uint8_diff"], s1["copy_bit_exact"]))
    A("")
    A("frame 24 の 5 realization も 1.10 の記録と突き合わせた（拡張が参照経路を"
      "変えていないことの確認）:")
    A("")
    A("| realization | 1.10 の行 | 生値の最大差 | 生値一致 | Sn の最大差 | Sn 一致 |")
    A("| --- | --- | --- | --- | --- | --- |")
    sn_map = {c["realisation"]: c for c in s2["v410_consistency"]["sn_frame24"]}
    for c in s2["v410_consistency"]["raw_frame24"]:
        sc = sn_map.get(c["realisation"])
        A("| %s | %s | %.1e | %s | %s | %s |"
          % (c["realisation"], c["v410_label"], c["max_abs_diff_raw"],
             c["exact_raw"],
             ("%.1e" % sc["max_abs_diff_sn"]) if sc else "—",
             sc["exact_sn"] if sc else "—"))
    A("")
    A("- 生値の全一致: **%s**、Sn の全一致: **%s**。"
      % (s2["v410_consistency"]["all_raw_exact"],
         s2["v410_consistency"]["all_sn_exact"]))
    A("")
    A("## 3. 方法（このジョブが測ったもの）")
    A("")
    A("### 3.1 対象")
    A("")
    A("| realization | seed | フレーム | レーン |")
    A("| --- | --- | --- | --- |")
    for rid, seed in s2["seeds"].items():
        A("| %s | %d | %s | warp-only（C = 1、合成器は 1.10 と同一） |"
          % (rid, seed, ", ".join(str(f) for f in s2["frames"])))
    A("")
    A("seed は 1.10 の numbers の `S_seed+N` から読み（base+offset と一致を"
      "確認）、+1/+2/+3/+7 の 5 本をそのまま使った。**無作為抽出ではない**"
      "（§9）。")
    A("")
    A("### 3.2 depth の定義")
    A("")
    A("G-LINE と同一。256×256 窓、生画像、理想的な軸平行 5–15 px 帯"
      "（直交軸は |k| ≤ 3 bin）、帯域 RMS / 窓平均 × 1000（‰）。窓位置は"
      "写真から決めた 4 窓（smooth x / smooth y / bark x / bark y）で、"
      "全画像で同じ座標。")
    A("")
    A("### 3.3 集計（realization 単位）")
    A("")
    A("各 realization で先にフレーム平均を取る（4 フレーム）。分布は R = 5 の"
      " realization 平均の列から作る。sd は realization 間の標本 sd、SE は "
      "sd/√R、90%% CI は realization を再抽出するブートストラップ"
      "（%d 回、5/95 パーセンタイル）。**実効 n = R = 5**。フレーム数を"
      " n に水増ししない。" % N_BOOT)
    A("")
    A("q90 は R=5 では補間値（5 点では 0.4×4 番目 + 0.6×5 番目）になる。"
      "参考線 1.25 との関係を見るために記録するが、推定としては薄い。")
    A("")
    A("### 3.4 フレーム相関 ρ")
    A("")
    A("realization 内で対応を取ったうえで、realization 間の相関を測る。"
      "窓 w・フレーム対 (f1, f2) について、5 realization の"
      " `(depth[r,f1,w], depth[r,f2,w])` の Pearson r。窓別の 4 本に加え、"
      "4 窓×5 対をプールした値と、4 窓の r の Fisher 平均を出す。"
      "相関は仮定せず実測値として報告する。")
    A("")
    A("### 3.5 prom2d（G-LINE2 の物差し）")
    A("")
    A("256×256 窓、`gray - gaussian(σ=3.0)`、Hann、2-D FFT、パワー |F|² の"
      "半径リング 256/15 … 256/5 cycles/window（5–15 px）の max/median。"
      "定義の再現確認: 本実装で写真の 24 窓中央値は **%.2f**（記録 21.1）、"
      "1.10 の warp-only は **25.53**（記録 25.4）。"
      % s2["photo_prom2d"]["tiles_24_median"])
    A("")
    A("窓別の prom2d は測定窓（smooth / bark）の 256² で測る。2-D リング統計"
      "なので軸を持たず、smooth x と smooth y は同じ prom2d を共有する"
      "（bark x / y も同様）。G-LINE2 のゲート集計（24 窓中央値）も併記するが、"
      "相関には窓別の値を使う。")
    A("")
    A("### 3.6 振幅正規化（Sn）")
    A("")
    A("1.10 と同方法: 同じ realization の傾斜場を、production の傾斜 RMS に"
      "スケールして再レンダする。目標は 1.10 と同じく**低域通過前**の"
      " production RMS（同時刻）。realization 側の RMS は実際に変位を作る"
      "**B_sigma 後**の傾斜 RMS。このため Sn の warp 平均は完全には揃わず、"
      "1.10 と同じく数 % 散らばる（§4.4）。production 自身は Sn = 生値。")
    A("")
    A("### 3.7 写真側の配置感度")
    A("")
    A("1.10 の配置感度 sd（上記）をそのまま使う。表現は"
      "「配置感度の約 ○ 倍」に限る。sd による有意差の言い方はしない。")
    A("")
    A("## 4. 結果 — 窓別 depth の分布")
    A("")
    A("### 4.1 主表（4 フレーム平均 → R=5 realization 平均の分布）")
    A("")
    A("| 窓 | 写真 [‰] | 平均 depth [‰] | 平均比 | sd（比） | SE（比） |"
      " 90% CI（平均比） | q90（比） | q90 の 90% CI |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for n in W:
        d = dist[n]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (n, f2(d["photo_depth_permille"]), f2(d["mean_depth_permille"]),
             f3(d["mean"]), f3(d["sd"]), f3(d["se"]), ci(d["ci90_mean"]),
             f3(d["q90"]), ci(d["ci90_q90"])))
    A("")
    A("実効 n = R = 5。CI は realization の再抽出による。")
    A("")
    A("### 4.2 realization 別の値（クラスタ＝realization）")
    A("")
    A("| 窓 | prod | s+1 | s+2 | s+3 | s+7 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for n in W:
        A("| %s | %s |" % (n, " | ".join(
            f3(v) for v in dist[n]["cluster_values"])))
    A("")
    A("各セルは 4 フレーム平均の比。分布の散らばりはこの 5 本の列から作った。")
    A("")
    A("### 4.3 フレーム別の平均比（realization 5 本の平均）")
    A("")
    A("| 窓 | frame 0 | frame 16 | frame 32 | frame 24 | 4 フレーム平均 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for n in W:
        A("| %s | %s | %s | %s | %s | %s |"
          % (n, f3(fw["0"][n]["mean_ratio"]), f3(fw["16"][n]["mean_ratio"]),
             f3(fw["32"][n]["mean_ratio"]), f3(fw["24"][n]["mean_ratio"]),
             f3(dist[n]["mean"])))
    A("")
    A("フレーム間の振れ（同じ窓の 4 列の範囲）は、realization 間の sd と"
      "同じ桁かそれ以上にある。フレームは独立標本ではないので、この振れを"
      " n の水増しには使わない。")
    A("")
    A("### 4.4 振幅正規化（Sn）版")
    A("")
    A("| 窓 | 平均比（Sn） | sd | SE | 90% CI | 平均超過 [‰] | 生値との差"
      "（平均比） |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for n in W:
        d = distn[n]
        A("| %s | %s | %s | %s | %s | %+.2f | %+.3f |"
          % (n, f3(d["mean"]), f3(d["sd"]), f3(d["se"]), ci(d["ci90_mean"]),
             d["mean_excess_permille"], d["mean"] - dist[n]["mean"]))
    A("")
    A("Sn でも 5 本の範囲は生値とほぼ同じ（下の warp 平均を参照）。"
      "したがって 4 窓の散らばりは、振幅の差ではなく"
      "**パターン（ realization ）の差**として残る。")
    A("")
    A("| realization | 生 warp 平均 [px] | Sn warp 平均 [px] | 比（Sn/生） |")
    A("| --- | --- | --- | --- |")
    for rid, rec in s2["per_realisation"].items():
        a = rec["frame_mean"]["warp_px_mean"]
        b = rec["frame_mean_sn"]["warp_px_mean"]
        A("| %s | %s | %s | %s |"
          % (rid, f3(a), f3(b), f3(b / a)))
    A("")
    A("production の傾斜 RMS（低域通過前 / 後、frame 24）: %s / %s"
      "（比 %.4f）。1.10 の Sn 目標は前の値なので、B_sigma 後の実振幅は"
      " realization ごとに数 %% ずれる。frame 24 だけを見ると生 warp は"
      " %s–%s px、Sn warp は %s–%s px（1.10 の Sn 記録 9.88–10.32 px と"
      "整合）。4 フレーム平均は上の表のとおりで、Sn は %s–%s px。"
      % ("%.4f" % s2["production_slope_rms"][str(R18.STILL_IDX)][
             "prelowpass"],
         "%.4f" % s2["production_slope_rms"][str(R18.STILL_IDX)][
             "lowpassed"],
         s2["production_slope_rms"][str(R18.STILL_IDX)]["pre_over_low"],
         f3(min(s2["per_realisation"][r]["per_frame"][str(R18.STILL_IDX)][
             "warp_px_mean"] for r in s2["per_realisation"])),
         f3(max(s2["per_realisation"][r]["per_frame"][str(R18.STILL_IDX)][
             "warp_px_mean"] for r in s2["per_realisation"])),
         f3(min(s2["per_realisation"][r]["per_frame"][str(R18.STILL_IDX)][
             "sn_warp_px_mean"] for r in s2["per_realisation"])),
         f3(max(s2["per_realisation"][r]["per_frame"][str(R18.STILL_IDX)][
             "sn_warp_px_mean"] for r in s2["per_realisation"])),
         f3(min(s2["per_realisation"][r]["frame_mean_sn"]["warp_px_mean"]
                for r in s2["per_realisation"])),
         f3(max(s2["per_realisation"][r]["frame_mean_sn"]["warp_px_mean"]
                for r in s2["per_realisation"]))))
    A("")
    A("## 5. フレーム相関 ρ（実測値）")
    A("")
    A("### 5.1 窓別（n = R = 5）")
    A("")
    A("| 窓 | 0–16 | 16–32 | 0–32 |")
    A("| --- | --- | --- | --- |")
    for n in W:
        A("| %s | %s | %s | %s |"
          % (n, rho_txt(s2["frame_corr"]["0-16"][n]["rho"]),
             rho_txt(s2["frame_corr"]["16-32"][n]["rho"]),
             rho_txt(s2["frame_corr"]["0-32"][n]["rho"])))
    A("")
    A("### 5.2 プール")
    A("")
    A("| 集計 | 0–16 | 16–32 | 0–32 |")
    A("| --- | --- | --- | --- |")
    A("| 窓内中心化プール（20 対） | %s | %s | %s |"
      % (rho_txt(s2["frame_corr"]["0-16"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_within_window_centered"]["rho"])))
    A("| 4 窓の Fisher 平均 | %s | %s | %s |"
      % (rho_txt(s2["frame_corr"]["0-16"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["16-32"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["0-32"]["fisher_mean_of_4_windows"])))
    A("| 生のプール（20 対、窓間水準差を含む・参考） | %s | %s | %s |"
      % (rho_txt(s2["frame_corr"]["0-16"]["pooled_over_windows"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_over_windows"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_over_windows"]["rho"])))
    A("")
    A("3 対×4 窓×5 realization の全 60 対をそのままプールすると r = %s だが、"
      "これは窓間の水準差をさらに強く含む（参考）。"
      % rho_txt(s2["frame_corr"]["all_pairs_all_windows"]["rho"]))
    A("")
    A("**プールの注意（重要）**: 生の 20 対プールが正に大きいのは、"
      "smooth 窓の depth（約 61‰）と bark 窓（約 42‰）という**窓間の水準差**が"
      "効いているためで、フレーム間で同じ realization を追う力ではない"
      "（Simpson 型）。窓内で中心化すると %s / %s / %s、4 窓の Fisher 平均は"
      " %s / %s / %s になる。フレーム相関として解釈できるのはこの 2 つ"
      "（窓別と窓内中心化）で、生のプールは brief の要求どおり記録するが"
      "解釈には使わない。"
      % (rho_txt(s2["frame_corr"]["0-16"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-16"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["16-32"]["fisher_mean_of_4_windows"]),
         rho_txt(s2["frame_corr"]["0-32"]["fisher_mean_of_4_windows"])))
    A("")
    A("読み方の注意: 各セルは 5 点の相関で、r の不確かさは大きい。窓別の値は"
      " %s–%s にわたり、同じ realization 内でも時刻によって帯の量の順位が"
      "変わることを示す。したがってフレームを独立な複製として n に数える"
      "ことはできない（実効 n = R = 5）。ここから合否は作らない。"
      % (f3(min(fc_w)), f3(max(fc_w))))
    A("")
    A("Sn 版（参考）: 窓内中心化プールは 0–16 %s / 16–32 %s / 0–32 %s。"
      % (rho_txt(s2["frame_corr_sn"]["0-16"][
          "pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr_sn"]["16-32"][
             "pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr_sn"]["0-32"][
             "pooled_within_window_centered"]["rho"])))
    A("")
    A("## 6. prom2d と depth の realization 間相関（G-LINE2 独立性チェック）")
    A("")
    A("### 6.1 prom2d の値（frame 24）")
    A("")
    A("| realization | smooth 窓 | bark 窓 | 24 窓中央値 | 24 窓最大 |")
    A("| --- | --- | --- | --- | --- |")
    for rid, v in s2["prom2d"]["per_realisation"].items():
        A("| %s | %s | %s | %s | %s |"
          % (rid, f2(v["smooth"]), f2(v["bark"]),
             f2(v["tiles_24_median"]), f2(v["tiles_24_max"])))
    A("| 写真 | %s | %s | %s | %s |"
      % (f2(s2["photo_prom2d"]["smooth"]), f2(s2["photo_prom2d"]["bark"]),
         f2(s2["photo_prom2d"]["tiles_24_median"]),
         f2(s2["photo_prom2d"]["tiles_24_max"])))
    A("")
    A("### 6.2 depth（frame 24）との相関")
    A("")
    A("| 窓 | prom2d の窓 | Pearson r | Spearman ρ | n |")
    A("| --- | --- | --- | --- | --- |")
    for n in W:
        c = pc[n]
        A("| %s | %s | %s | %s | %d |"
          % (n, c["prom2d_window"], rho_txt(c["rho_pearson"]),
             rho_txt(c["rho_spearman"]), c["n"]))
    A("| 全体（20 対プール） | — | %s | %s | %d |"
      % (rho_txt(pc["pooled_20_pairs"]["rho_pearson"]),
         rho_txt(pc["pooled_20_pairs"]["rho_spearman"]),
         pc["pooled_20_pairs"]["n_pairs"]))
    A("| 全体（窓内の 2 軸を平均した 10 対） | — | %s | %s | %d |"
      % (rho_txt(pc["pooled_10_unique"]["rho_pearson"]),
         rho_txt(pc["pooled_10_unique"]["rho_spearman"]),
         pc["pooled_10_unique"]["n_pairs"]))
    A("")
    A("窓別は **n = 5** の相関である。smooth x と smooth y は同じ prom2d を"
      "共有するので、その 2 行は独立ではない。20 対プールも x が窓内で"
      "重複する。**独立性の主張はしない。相関の実測値を出したところまでが"
      "本項の仕事**（brief item 4）。")
    A("")
    A("## 7. bark y の系統的超過は分布で残るか")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 写真 depth | %s ‰ |" % f2(by["photo_depth_permille"]))
    A("| 1.10 本文の表示 | 1.120（+12.0%） |")
    A("| 1.10 記録値の再計算（生 / Sn） | %s / %s |"
      % (f3(by["v410_frame24_mean_ratio_raw"]),
         f3(by["v410_frame24_mean_ratio_sn"])))
    A("| 本ジョブ frame 24（生 / Sn） | %s / %s（一致 %s / %s） |"
      % (f3(by["frame24_mean_ratio_raw"]),
         f3(by["frame24_mean_ratio_sn"]),
         by["frame24_raw_matches_v410_raw"],
         by["frame24_sn_matches_v410_sn"]))
    A("| 4 フレーム平均比（R=5、生） | %s |" % f3(by["mean_ratio_4frames"]))
    A("| その 90%% CI | %s |" % ci(by["ci90_mean_4frames"]))
    A("| realization 別（prod / +1 / +2 / +3 / +7） | %s |"
      % " / ".join(f3(v) for v in by["per_realisation_mean_ratio"]))
    A("| frame 別の 5 本平均（0 / 16 / 32 / 24） | %s |"
      % " / ".join(f3(by["per_frame_mean_ratio"][str(f)])
                   for f in s2["frames"]))
    A("| 1 超のセル（realization×frame） | %d / %d |"
      % (by["n_cells_above_1"], by["n_cells"]))
    A("| 平均超過（4 フレーム平均） | %+.2f ‰（配置感度の約 %s 倍） |"
      % (by["mean_excess_permille_4frames"],
         f1(by["excess_over_placement_sd_4frames"])))
    A("")
    A("**+12%% の内訳（1.10 との突き合わせ）**: 1.10 の本文は bark y を"
      "「+12.0%%（1.120）」と書くが、1.10 が残した realization 別の値を"
      "再計算すると、**生値の平均は %s（+11.0%%）、Sn 版の平均は %s"
      "（+12.0%%）**である。つまり 1.10 の 1.120 は Sn の平均で、生値は"
      " 1.110。本ジョブの frame 24 は生値 %s / Sn %s で両方を桁まで再現した。"
      "したがって「+12%%」を生値の水準として引用すると 1.0 ポイント"
      "過大になる。"
      % (f3(by["v410_frame24_mean_ratio_raw"]),
         f3(by["v410_frame24_mean_ratio_sn"]),
         f3(by["frame24_mean_ratio_raw"]),
         f3(by["frame24_mean_ratio_sn"])))
    A("")
    A("分布で見ると、bark y の上方向は残る。4 フレーム平均の中心は %s"
      "（写真の %+.1f%%）で、20 セル中 %d セルが 1 超、1 を下回る"
      " realization は %d 本。ただし大きさは基準で変わる: frame 24・生で"
      " +11.0%%、frame 24・Sn で +12.0%%、4 フレーム平均・生で %+.1f%%。"
      % (f3(by["mean_ratio_4frames"]),
         100.0 * (by["mean_ratio_4frames"] - 1.0),
         by["n_cells_above_1"],
         sum(1 for v in by["per_realisation_mean_ratio"] if v <= 1.0),
         100.0 * (by["mean_ratio_4frames"] - 1.0)))
    A("")
    A("他の 3 窓の中心は %s。smooth x / smooth y / bark x の中心は 1 の"
      "近傍ないし下にあり、bark y だけが分布の中心として上にある、という"
      " 1.10 §0.2 の形が 4 フレーム平均でも保たれている。"
      % " / ".join("%s %s" % (n, f3(dist[n]["mean"]))
                   for n in W if n != "bark y"))
    A("")
    A("## 8. 参考線（判定ではない）")
    A("")
    A("採用設計の 1.05 / 1.25 は**暫定配置値（未検証）**である。以下の表は"
      "算術差を示すだけで、合否・採否・「未解決」のいずれも本ジョブは"
      "出さない。")
    A("")
    A("| 窓 | 平均比 | 参考線 1.05 との差 | q90 | 参考線 1.25 との差 |")
    A("| --- | --- | --- | --- | --- |")
    for n in W:
        d = dist[n]
        A("| %s | %s | %+.3f | %s | %+.3f |"
          % (n, f3(d["mean"]), d["mean_minus_reference_mean"],
             f3(d["q90"]), d["q90_minus_reference_q90"]))
    A("")
    A("q90 は R=5 の補間値で、その 90%% CI は %s。"
      % " / ".join("%s %s" % (n, ci(dist[n]["ci90_q90"])) for n in W))
    A("")
    A("## 9. 限界・開示")
    A("")
    A("1. **R = 5 は薄い。** realization 間 sd の推定も、ブートストラップ CI も、"
      "5 クラスタから作っている（5⁵ = 3125 通りの再抽出しかない）。"
      "CI の幅は「データが許す幅」であって、モデルの realization 分布の"
      "幅ではない。")
    A("2. **写真は 1 枚。** 4 窓は 1 枚の写真から選んだ固定座標で、写真側の"
      "分布は無い。配置感度 sd は無歪み画像で窓原点を ±24 px 振ったときの"
      "統計の揺れであり、写真を複数枚撮ったときのばらつきではない。")
    A("3. **seed の選び方（+1/+2/+3/+7）は無作為抽出ではない。** production の"
      "近傍の 4 本で、モデルの realization 分布を代表する標本ではない。"
      "したがって「平均比」はこの 5 本の平均であって、モデルの期待値ではない。")
    A("4. **フレーム {0, 16, 32, 24} は brief の指定で、無作為でも等間隔でも"
      "ない。** ループは 48 フレーム 2 秒。4 フレーム平均はループ平均では"
      "ない。フレームは独立標本ではない（実効 n = 5）。")
    A("5. **相関は n = 5。** 窓別 r は 5 点で、1 点の移動で大きく変わる。"
      "プールしても窓内で x が重複する。独立性の判定はしない。")
    A("6. **Sn の振幅一致は厳密でない。** 1.10 と同じ目標（低域通過前の"
      " production RMS）を使うため、B_sigma 後の実振幅は数 % ずれる"
      "（§4.4 の warp 平均）。パターン差（±10–30%）に比べれば小さいが、"
      "「振幅を完全に揃えた」とは言えない。")
    A("7. **prom2d は再実装。** brief の定義（リング 5–15 px、max/median、"
      "256²、Hann、HP σ=3.0）を書き、記録値（写真 21.1、作品 25.4）との"
      "整合を確認した。ただし G-LINE2 の元実装そのものではない。窓別"
      " prom2d はゲートの 24 窓中央値とは別の量（テクスチャで選んだ窓）である。")
    A("8. **depth は 1 次元周辺統計。** 窓の位置・大きさ・帯の定義は固定で、"
      "その恣意性は残る。帯の方向性は軸平行の 2 方向しか見ていない。")
    A("9. **目視は使っていない。** 判断はすべて数値。")
    A("10. **閾値は参考線。** 1.05 / 1.25 は配置値で、統計的導出でも、"
      "新規独立 seed での検証済みでもない。本ジョブは合否を出さない。")
    A("11. **4 フレーム平均は frame 24 を含む。** frame 24 は 1.10 の参照 still"
      "でもあるが、重み付けはしていない。")
    A("12. **生のプール相関は窓間の水準差を含む。** smooth 窓（約 61‰）と"
      " bark 窓（約 42‰）をそのまま並べると、フレームと無関係な水準差で"
      "正の相関が出る（Simpson 型）。§5 では窓別・窓内中心化・Fisher 平均を"
      "主に読み、生のプールは記録として残す。")
    A("13. **1.10 の「+12%」は Sn 基準。** 1.10 の記録値から再計算すると、"
      "生値の bark y 平均は 1.110（+11.0%）、Sn 平均は 1.120（+12.0%）。"
      "本ジョブは両方を再現した。brief の +12% を生値の水準として使うと"
      " 1.0 ポイント過大になる（§7）。")
    A("")
    A("## 10. 言えること / 言えないこと（判定ではなく、数値の読み）")
    A("")
    A("- **言える**: 4 フレーム平均の比の中心は %s、realization 間 sd は %s。"
      % (" / ".join(f3(dist[n]["mean"]) for n in W),
         " / ".join(f3(dist[n]["sd"]) for n in W)))
    A("- **言える**: bark y の frame 24 平均比は 1.10 の生値 %s / Sn %s と"
      "桁まで一致し、4 フレーム平均の中心は %s。上方向は分布にしても残る。"
      "ただし 1.10 の「+12%%」は Sn 基準で、生値は +11.0%%。"
      % (f3(by["v410_frame24_mean_ratio_raw"]),
         f3(by["v410_frame24_mean_ratio_sn"]),
         f3(by["mean_ratio_4frames"])))
    A("- **言える**: フレーム相関は窓別で %s–%s、窓内中心化プールで %s / %s /"
      " %s。同じ realization 内でも時刻で帯の量の順位が変わるので、"
      "フレームを独立標本として n に数えない扱いと整合する。"
      % (f3(min(fc_w)), f3(max(fc_w)),
         rho_txt(s2["frame_corr"]["0-16"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["16-32"]["pooled_within_window_centered"]["rho"]),
         rho_txt(s2["frame_corr"]["0-32"]["pooled_within_window_centered"]["rho"])))
    A("- **言える**: prom2d と depth の realization 間相関は、窓別で %s、"
      "20 対プールで %s。窓別の値は §6.2 のとおり。"
      % (" / ".join(rho_txt(pc[n]["rho_pearson"]) for n in W),
         rho_txt(pc["pooled_20_pairs"]["rho_pearson"])))
    A("- **言えない**: この 5 本がモデルの realization 分布を代表するか。")
    A("- **言えない**: 写真 1 枚の 4 窓が「写真の分布」を代表するか。")
    A("- **言えない**: 1.05 / 1.25 が正しい配置か（未検証）。")
    A("- **言えない**: prom2d と depth が独立か（相関を測っただけ）。")
    A("")
    A("## 11. ゲートと無変更の記録")
    A("")
    A("本ジョブは診断なので、ゲートは「再現したこと」と「触っていないこと」"
      "の確認として記録する。")
    A("")
    A("| 確認 | 状態 | 数値 |")
    A("| --- | --- | --- |")
    A("| STAGE 1 再現（4 窓 depth） | 完全一致 | 最大差 %.1e |"
      % s1["max_abs_diff"])
    A("| STAGE 1 再現（G-WARP） | 完全一致 | %.9f px |" % s1["g_warp_mean_px"])
    A("| frame 24 の 1.10 整合（生値） | %s | 5 realization |"
      % s2["v410_consistency"]["all_raw_exact"])
    A("| frame 24 の 1.10 整合（Sn） | %s | 4 realization |"
      % s2["v410_consistency"]["all_sn_exact"])
    A("| 主要 10 ファイルの sha256 | %s | 開始＝終了 |"
      % ("変化なし" if s3["major_unchanged"] else "**変化あり**"))
    A("| 本ジョブが書き込みで開いた .py | 0 本 | 書き込み先は "
      "`out/v4-1-12/` の 3 ファイルのみ |")
    A("| 全 .py で開始/終了の sha が変わったもの | %d 本（本ジョブ外） | %s |"
      % (len(s3["py_legacy"]["changed"]),
         ", ".join("`%s`" % k for k in s3["py_legacy"]["changed"]) or "なし"))
    A("")
    A("### 11.1 sha256（開始 → 終了）")
    A("")
    A("| ファイル | 開始（先頭 16） | 終了（先頭 16） | 変化 |")
    A("| --- | --- | --- | --- |")
    for name in sorted(s3["major_hashes_end"]):
        a = s3["major_hashes_start"].get(name, {}).get("sha256", "—")[:16]
        b = s3["major_hashes_end"][name]["sha256"][:16]
        A("| `%s` | `%s` | `%s` | %s |"
          % (os.path.basename(name), a, b,
             "なし" if a == b else "**あり**"))
    A("")
    A("全 .py の変化（開始 %d 本 / 終了 %d 本）: %s。実行中に現れた .py: %s。"
      "欠落: %s。"
      % (s3["py_legacy"]["n_start"], s3["py_legacy"]["n_end"],
         ", ".join(s3["py_legacy"]["changed"]) or "なし",
         ", ".join(s3["py_legacy"]["added"]) or "なし",
         ", ".join(s3["py_legacy"]["missing"]) or "なし"))
    if s3["py_legacy"]["changed"] or s3["py_legacy"]["added"]:
        A("")
        if s3["py_legacy"]["changed"]:
            A("変化した .py の mtime（開始 → 終了）:")
            A("")
            A("| ファイル | sha 開始 → 終了 | mtime 開始 | mtime 終了 |")
            A("| --- | --- | --- | --- |")
            for k, d in s3["py_legacy"]["changed_detail"].items():
                A("| `%s` | `%s` → `%s` | %s | %s |"
                  % (k, d["sha_start"], d["sha_end"], d["mtime_start"],
                     d["mtime_end"]))
            A("")
        if s3["py_legacy"]["added"]:
            A("実行中に現れた .py の mtime: %s。"
              % ", ".join("`%s` %s" % (
                  k, s3["python_mtimes_end"].get(k))
                  for k in s3["py_legacy"]["added"]))
            A("")
        A("これらは本ジョブの書き込みではない。本ドライバは .py を読み取り専用で"
          "開き、書き込みモードで開くのは `out/v4-1-12/` の JSON/MD だけである。"
          "mtime が本ジョブの実行中なら、並行ジョブの成果物が更新・追加された"
          "ということ（本ジョブはそれらを import も write もしていない）。")
    A("")
    A("## 12. 成果物")
    A("")
    A("| ファイル | 役割 |")
    A("| --- | --- |")
    A("| `REPORT-V4-1-12.md` | 本レポート |")
    A("| `numbers-v4-1-12.json` | 数値の全量（分布・相関・realization 別・"
      "ハッシュ） |")
    A("| `checkpoint.json` | jobkit の段階記録（stage0 / stage1 / stage2 / "
      "report） |")
    A("")
    A("最終ファイルの bytes と sha256 は `checkpoint.json` の report ユニットに"
      "記録した。画像は保存していない（測定はメモリ上で完結）。")
    A("")
    text = "\n".join(L) + "\n"
    with open(REPORT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    log("  %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    log("  %s (%d bytes)" % (NUMBERS, os.path.getsize(NUMBERS)))
    log("  wall %.1fs" % (time.time() - t0))
    return {"report_bytes": os.path.getsize(REPORT),
            "numbers_bytes": os.path.getsize(NUMBERS)}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    log("BRIEF-V4-1-12 -- G-LINE-stat first measurement (diagnostic only)")
    s0 = stage0()
    s1 = stage1(s0)
    s2 = stage2(s0, s1)

    # end-of-job hashes and deliverable bookkeeping
    major_end = {p: {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
                 for p, _w in REF_PATHS}
    py_end = py_set()
    py_end_mtimes = py_mtimes()
    legacy = what_changed(s0["python_hashes_start"], py_end,
                          s0.get("python_mtimes_start"), py_end_mtimes)
    major_unchanged = all(
        major_end[p]["sha256"] == s0["major_hashes_start"][p]["sha256"]
        for p in major_end)
    built = now_iso()

    s3 = {
        "built_utc": built,
        "major_hashes_start": s0["major_hashes_start"],
        "major_hashes_end": major_end,
        "py_legacy": legacy,
        "major_unchanged": bool(major_unchanged),
        "python_files_n_end": len(py_end),
        "python_hashes_end": py_end,
        "python_mtimes_end": py_end_mtimes,
        "deliverables": {
            "numbers": "numbers-v4-1-12.json",
            "report": "REPORT-V4-1-12.md",
            "checkpoint": "checkpoint.json",
        },
    }

    numbers = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_utc": built,
        "new_files": NEW_FILES,
        "run_v4_1_12_sha256": sha256_file(os.path.abspath(__file__)),
        "scope": ("diagnostic measurement only; no pass/fail verdict; "
                  "1.05/1.25 are provisional placement values"),
        "stage0": s0,
        "stage1": s1,
        "stage2": s2,
        "stage3": s3,
        "files": {"numbers": "numbers-v4-1-12.json",
                  "report": "REPORT-V4-1-12.md",
                  "checkpoint": "checkpoint.json"},
    }
    write_json(NUMBERS, numbers)
    write_report(s0, s1, s2, s3)
    log("  %s (%d bytes)" % (NUMBERS, os.path.getsize(NUMBERS)))

    final = {}
    for p in (NUMBERS, REPORT):
        final[os.path.basename(p)] = {
            "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    cp_extra("report", {
        "at": now_iso(),
        "files": final,
        "checkpoint_bytes": (os.path.getsize(CP_PATH)
                             if os.path.isfile(CP_PATH) else None),
        "py_legacy": legacy,
        "note": ("diagnostic measurement only; no pass/fail; checkpoint's own "
                 "sha256 is not recorded (self-referential)"),
    })
    log("=" * 78)
    log("DONE  %.1fs" % (time.time() - t0))
    log("  %s" % REPORT)
    log("  %s" % NUMBERS)
    log("  %s" % CP_PATH)


if __name__ == "__main__":
    main()
