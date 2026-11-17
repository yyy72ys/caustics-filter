"""run_v4_1_13 -- BRIEF-V4-1-13: G-LINE-stat threshold validation (fresh seeds, R=10).

THIS IS A THRESHOLD VALIDATION, NOT A MODEL PASS/FAIL.  The job produces no
overall verdict; it reports per-window estimates and the three pre-registered
decision rules (i)(ii)(iii) as held / not held.

NEW FILES ONLY: this driver (plus out/v4-1-13/ deliverables).  run_v4_1_12.py
is imported unchanged and supplies the measurement pipeline (context, state,
warp-only render, G-LINE depth, prom2d, statistics helpers).

Pre-registration: D:\\mutaro-work\\rally-in\\PRE-REG-GLINE-STAT.md
  * seeds V1..V10 fixed and used in this order, no substitution
  * frames 0 / 16 / 32 are the validation frames; frame 24 is recorded for
    continuity only and is never used in a decision
  * per-window statistic = realisation representative value (mean of the three
    validation-frame ratios); mean / SE / 90 % CI computed over R realisations
    (effective n = R), CI by realisation-level cluster bootstrap, 5 000 draws
  * q90 reported as a reference statistic with its cluster-bootstrap CI width
  * thresholds 1.05 (mean) / 1.25 (q90) unchanged, never proposed anew
  * unresolved rule: a window whose mean CI straddles 1.05 gets +3
    realisations (cap R = 15); extra seeds are positions 11..15 of the same
    generator default_rng(20261114).integers(0, 2**63, ...)

STAGES (foreground, jobkit checkpoint per stage)
  --stage0   preflight + hashes + seed-sequence check + control reproduction
  --stage1   R = 10 fresh seeds x frames {0, 16, 32, 24}: depth, Sn, prom2d
  --stage2   realisation-level aggregation + decision rules + expansion
  --report   numbers-v4-1-13.json + REPORT-V4-1-13.md
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
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, read_checkpoint  # noqa: E402

import run_v4_1_12 as R12   # noqa: E402  (1.12 measurement implementation: UNCHANGED)
import run_v4_1_10 as R10   # noqa: E402  (warp-only route: UNCHANGED)
import run_v4_1_8 as R18    # noqa: E402  (scene + G-LINE statistic: UNCHANGED)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-13.md"
PREREG = r"D:\mutaro-work\rally-in\PRE-REG-GLINE-STAT.md"

OUT = os.path.join(HERE, "out", "v4-1-13")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-1-13.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-13.md")

V412_DIR = os.path.join(HERE, "out", "v4-1-12")
V412_NUMBERS = os.path.join(V412_DIR, "numbers-v4-1-12.json")
V412_REPORT = os.path.join(V412_DIR, "REPORT-V4-1-12.md")
V412_CP = os.path.join(V412_DIR, "checkpoint.json")

WINDOW_ORDER = R12.WINDOW_ORDER
WINDOW_TYPE = R12.WINDOW_TYPE

FRAMES_REF = [0, 16, 32, 24]      # brief: 24 measured for reference only
VALID_FRAMES = [0, 16, 32]        # the validation frames
FRAME_PAIRS = [(0, 16), (16, 32), (0, 32)]

N_BOOT = 5000                     # pre-reg: 5 000 cluster bootstrap draws
BOOT_SEED = 20261114              # same generator seed as the pre-registration
THRESH_MEAN = 1.05                # unchanged provisional placement value
THRESH_Q90 = 1.25                 # unchanged provisional placement value

SEED_GEN = 20261114
SEEDS_10 = [
    8760839803807278262,   # V1
    1533652889374406734,   # V2
    3629938633691225283,   # V3
    1808497988179510849,   # V4
    6792552148387992806,   # V5
    8651395248935177417,   # V6
    8482304938715652794,   # V7
    100472280562746898,    # V8
    8158281184839309514,   # V9
    2853851374193069971,   # V10
]

NEW_FILES = ["run_v4_1_13.py"]

REF_PATHS = [
    (BRIEF, "brief (1.13)"),
    (PREREG, "pre-registration (normative)"),
    (os.path.join(HERE, "run_v4_1_12.py"), "run_v4_1_12.py (1.12 measurement impl.)"),
    (V412_NUMBERS, "numbers-v4-1-12.json (control reference)"),
    (V412_REPORT, "REPORT-V4-1-12.md (1.12 report)"),
    (os.path.join(HERE, "run_v4_1_10.py"), "run_v4_1_10.py (warp-only route)"),
    (os.path.join(HERE, "run_v4_1_8.py"), "run_v4_1_8.py (scene + G-LINE)"),
    (os.path.join(HERE, "run_v4_1_9.py"), "run_v4_1_9.py (G-WARP yardstick)"),
    (os.path.join(HERE, "stream_surface2.py"), "stream_surface2.py"),
    (os.path.join(HERE, "v5_composite.py"), "v5_composite.py"),
    (os.path.join(HERE, "run_v16_water.py"), "run_v16_water.py"),
    (os.path.join(HERE, "transport_mesh.py"), "transport_mesh.py"),
    (R18.PHOTO, "photo"),
]

# the 4 paths STAGE 0 of the brief demands, on their own (STOP if any missing)
REQUIRED_4 = [
    os.path.join(HERE, "run_v4_1_12.py"),
    V412_NUMBERS,
    V412_REPORT,
    PREREG,
]


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
            "n_changed": len(changed), "n_added": len(added),
            "n_missing": len(missing),
            "n_start": len(start), "n_end": len(end),
            # this driver opens .py files READ-ONLY; its only write targets
            # are CP_PATH / NUMBERS / REPORT under out/v4-1-13
            "py_write_targets": []}


# ---------------------------------------------------------------------------
# statistics (same helpers as 1.12, N_BOOT = 5000 per the pre-reg)
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
    """Mean/sd/SE/90 % CI/q90 from R realisation-level values."""
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
            "q90_ci_width": float(q90_ci[1] - q90_ci[0]),
            "cluster_values": [float(v) for v in a]}


def rms_of(gx, gy):
    gx = np.asarray(gx, np.float64)
    gy = np.asarray(gy, np.float64)
    return float(np.sqrt(np.mean(gx * gx + gy * gy)))


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight, hashes, seed check, control reproduction
# ---------------------------------------------------------------------------
REQUIRED_MEAN_KEYS = ["depth_permille", "g_warp_mean_px", "photo_depth_permille"]


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] preflight + sha256 + seed check + control reproduction")
    os.makedirs(OUT, exist_ok=True)

    preflight = {}
    missing = []
    for p, what in REF_PATHS:
        ok = os.path.isfile(p)
        preflight[p] = {"what": what, "exists": bool(ok),
                        "bytes": (os.path.getsize(p) if ok else None)}
        if not ok:
            missing.append(p)
    required_missing = [p for p in REQUIRED_4 if not os.path.isfile(p)]
    if missing or required_missing:
        raise SystemExit("STOP: missing prerequisites: %s"
                         % (required_missing or missing))

    major_start = {p: {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
                   for p, _w in REF_PATHS}
    py_start = py_set()

    # the three 1.12 deliverables: explicit start fingerprints (integrity)
    v412_files_start = {}
    for p in (V412_NUMBERS, V412_REPORT, V412_CP):
        v412_files_start[os.path.basename(p)] = {
            "bytes": os.path.getsize(p), "sha256": sha256_file(p)}

    # --- seed sequence: must reproduce the pre-registered 10, in order -----
    rng = np.random.default_rng(SEED_GEN)
    gen15 = [int(x) for x in rng.integers(0, 2 ** 63, size=15)]
    seeds_match = [gen15[i] == SEEDS_10[i] for i in range(10)]
    if not all(seeds_match):
        raise SystemExit("STOP: seed sequence does not reproduce the pre-reg list")

    # --- 1.12 control reference -------------------------------------------
    v412 = read_json(V412_NUMBERS)
    s1 = v412.get("stage1", {})
    for k in REQUIRED_MEAN_KEYS:
        if k not in s1:
            raise SystemExit("STOP: numbers-v4-1-12.json stage1 lacks %r" % k)
    target_depth = [float(x) for x in s1["depth_permille"]]
    target_gw = float(s1["g_warp_mean_px"])
    target_sha = str(s1.get("image_pixel_sha256"))
    target_photo = [float(x) for x in s1["photo_depth_permille"]]

    v412_dist = v412["stage2"]["distribution"]
    v412_means = {w: float(v412_dist[w]["mean"]) for w in WINDOW_ORDER}
    v412_cis = {w: [float(x) for x in v412_dist[w]["ci90_mean"]]
                for w in WINDOW_ORDER}

    # minimal size/shape sanity of the 1.12 numbers before trusting them
    if len(target_depth) != 4 or len(target_photo) != 4:
        raise SystemExit("STOP: 1.12 control record has unexpected length")

    # --- control: 1.12 production seed, same pipeline, frame 24 -----------
    log("  control: reproducing 1.12 production seed frame 24 ...")
    ctx = R12.make_context()
    img_ref, _dbg = R10.render_ref(ctx)          # untouched v16 entry point
    depth = R12.measure_depths(img_ref)
    diff = [float(a - b) for a, b in zip(depth, target_depth)]
    exact = bool(all(a == b for a, b in zip(depth, target_depth)))

    gx_s, gy_s = R10.state_gradients(ctx["st"], ctx["fields"], "current")
    grad = (gx_s, gy_s, ctx["S"], "bicubic")
    img_cp = R10.render_variant(ctx, grad)
    dmax = int(np.max(np.abs(img_ref.astype(np.int32)
                             - img_cp.astype(np.int32))))
    copy_exact = bool(dmax == 0 and np.array_equal(img_ref, img_cp))
    img_sha = R18.sha256_bytes(img_ref.tobytes())

    bx, by, _ = R10.disp_at(grad, R12.WIDTH, R12.HEIGHT,
                            float(ctx["env"]["depthM"]))
    yard = R10.warp_yardstick(bx, by, ctx["S"])
    gw = float(yard["px_mean"])

    photo_depth = [float(r["depth_e3"]) for r in ctx["photo_rows"]]
    photo_match = bool(all(abs(a - b) < 1e-12
                           for a, b in zip(photo_depth, target_photo)))

    control = {
        "seed": int(R18.SCN["seed"]),
        "frame": int(R18.STILL_IDX),
        "window_order": WINDOW_ORDER,
        "depth_permille": depth,
        "v412_depth_permille": target_depth,
        "diff_vs_v412": diff,
        "max_abs_diff": float(max(abs(d) for d in diff)),
        "exact_equal": exact,
        "g_warp_mean_px": gw,
        "v412_g_warp_mean_px": target_gw,
        "g_warp_diff": float(gw - target_gw),
        "g_warp_exact_equal": bool(gw == target_gw),
        "image_pixel_sha256": img_sha,
        "v412_image_pixel_sha256": target_sha,
        "image_sha_equal": bool(img_sha == target_sha),
        "copy_bit_exact": copy_exact,
        "copy_max_uint8_diff": dmax,
        "photo_depth_permille": photo_depth,
        "v412_photo_depth_permille": target_photo,
        "photo_depth_equal": photo_match,
        "wall_s": float(time.time() - t0),
    }
    log("  control depth %s  (diff max %.3e, exact %s)"
        % ([round(x, 6) for x in depth], control["max_abs_diff"], exact))
    log("  control G-WARP %.9f px  (diff %+.3e, exact %s)"
        % (gw, control["g_warp_diff"], control["g_warp_exact_equal"]))
    log("  control pixel sha equal to 1.12 record: %s" % control["image_sha_equal"])
    del ctx, img_ref, img_cp
    gc.collect()

    if not (exact and control["g_warp_exact_equal"]):
        raise SystemExit("STOP: control did not reproduce 1.12 exactly")
    if not photo_match:
        raise SystemExit("STOP: photo depths differ from the 1.12 record")

    out = {
        "at": now_iso(),
        "preflight": preflight,
        "required_4": REQUIRED_4,
        "required_4_all_exist": True,
        "major_hashes_start": major_start,
        "v412_files_start": v412_files_start,
        "python_files_n_start": len(py_start),
        "python_hashes_start": py_start,
        "python_mtimes_start": py_mtimes(),
        "seed_generator": {
            "call": "numpy.random.default_rng(%d).integers(0, 2**63, size=15)"
                    % SEED_GEN,
            "numpy": np.__version__,
            "generated_15": gen15,
            "seeds_10_used": SEEDS_10,
            "matches_preregistered_10": bool(all(seeds_match)),
            "match_mask": seeds_match,
            "continuation_11_13": gen15[10:13],
            "continuation_14_15": gen15[13:15],
            "note": ("positions 11..15 are the continuation of the same "
                     "generator, exactly as the pre-reg / brief spell out"),
        },
        "v412_control_reference": {
            "source": V412_NUMBERS,
            "depth_permille": target_depth,
            "g_warp_mean_px": target_gw,
            "image_pixel_sha256": target_sha,
            "photo_depth_permille": target_photo,
            "v412_means_R5_4frames": v412_means,
            "v412_ci90_R5_4frames": v412_cis,
        },
        "control": control,
        "thresholds": {"mean": THRESH_MEAN, "q90": THRESH_Q90,
                       "status": "provisional placement values, unchanged"},
        "frames": {"measured": FRAMES_REF, "validation": VALID_FRAMES,
                   "reference_only": [24]},
        "brief": BRIEF,
        "prereg": PREREG,
        "python": sys.version,
    }
    log("  prerequisites: %d/%d exist (required 4: all)"
        % (len(REF_PATHS) - len(missing), len(REF_PATHS)))
    log("  seeds reproduce pre-reg list: %s  next3 %s"
        % (all(seeds_match), gen15[10:13]))
    log("  .py at start: %d" % len(py_start))

    cp_extra("stage0", {
        "at": out["at"],
        "preflight": preflight,
        "required_4_all_exist": True,
        "v412_files_start": v412_files_start,
        "seed_matches_preregistered_10": bool(all(seeds_match)),
        "control": {k: v for k, v in control.items() if k != "wall_s"},
        "n_python_start": len(py_start),
        "brief": BRIEF,
        "prereg": PREREG,
        "python": sys.version,
        "numpy": np.__version__,
    })
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- R = 10 fresh independent seeds
# ---------------------------------------------------------------------------
def prod_rms_targets(ctx):
    hi = {}
    lo = {}
    for f in FRAMES_REF:
        fl = R12.fields_for(ctx["st"], f, ctx, True)
        hi[f] = rms_of(fl["gx"], fl["gy"])                  # pre-lowpass
        gx_s, gy_s = R10.state_gradients(ctx["st"], fl, "current")
        lo[f] = rms_of(gx_s, gy_s)                          # post-lowpass
    return hi, lo


def measure_realisation(ctx, rid, seed, prod_hi, photo_depth, S, depth_m):
    st = R12.make_state_seed(seed)
    rec = {"realisation": rid, "seed": int(seed),
           "per_frame": {}, "frame24": None}
    img24 = None
    for f in FRAMES_REF:
        fl = st.fields_at(float(f) / float(R12.FPS))
        gx_s, gy_s = R10.state_gradients(st, fl, "current")
        grad = (gx_s, gy_s, S, "bicubic")
        img = R10.render_variant(ctx, grad)
        rows = R12.measure_rows(img)
        d = [float(r["depth_e3"]) for r in rows]
        r = [float(x / p) for x, p in zip(d, photo_depth)]

        rms_lo = rms_of(gx_s, gy_s)
        scale = float(prod_hi[f] / rms_lo)       # 1.10's Sn target
        grad_sn = (gx_s * scale, gy_s * scale, S, "bicubic")
        img_sn = R10.render_variant(ctx, grad_sn)
        d_sn = R12.measure_depths(img_sn)
        r_sn = [float(x / p) for x, p in zip(d_sn, photo_depth)]

        bx, by, _nz = R10.disp_at(grad, R12.WIDTH, R12.HEIGHT, depth_m)
        yard = R10.warp_yardstick(bx, by, S)

        rec["per_frame"][str(f)] = {
            "depth_permille": d,
            "ratio_vs_photo": r,
            "sn_scale": scale,
            "sn_depth_permille": d_sn,
            "sn_ratio_vs_photo": r_sn,
            "slope_rms_lowpassed": rms_lo,
            "sn_target_prelowpass": float(prod_hi[f]),
            "warp_px_mean": float(yard["px_mean"]),
            "period_px": [float(x["period_px"]) for x in rows],
            "prom_1d": [float(x["prom_peak_over_median"]) for x in rows],
        }
        if int(f) == 24:
            img24 = img
            rec["frame24"] = {
                "depth_permille": d,
                "ratio_vs_photo": r,
                "sn_depth_permille": d_sn,
                "sn_ratio_vs_photo": r_sn,
                "prom2d": R12.prom2d_block(
                    np.asarray(img[:, :, 1], np.float64)),
            }
        log("    %-4s f%-2d depth %s  ratio %s  warp %.3f"
            % (rid, f, [round(x, 2) for x in d],
               [round(x, 3) for x in r], yard["px_mean"]))
    del st, img24
    gc.collect()

    # representative value: mean of the three validation frames (raw ratio)
    rec["mean_validation"] = {
        "frames": VALID_FRAMES,
        "depth_permille": [float(np.mean(
            [rec["per_frame"][str(f)]["depth_permille"][w]
             for f in VALID_FRAMES])) for w in range(4)],
    }
    rec["mean_validation"]["ratio_vs_photo"] = [
        float(x / p) for x, p in zip(rec["mean_validation"]["depth_permille"],
                                     photo_depth)]
    rec["mean_validation_sn"] = {
        "frames": VALID_FRAMES,
        "depth_permille": [float(np.mean(
            [rec["per_frame"][str(f)]["sn_depth_permille"][w]
             for f in VALID_FRAMES])) for w in range(4)],
    }
    rec["mean_validation_sn"]["ratio_vs_photo"] = [
        float(x / p) for x, p in zip(
            rec["mean_validation_sn"]["depth_permille"], photo_depth)]
    # 4-frame mean incl. frame 24: continuity only, never a decision input
    rec["mean_4frame_reference"] = {
        "frames": FRAMES_REF,
        "depth_permille": [float(np.mean(
            [rec["per_frame"][str(f)]["depth_permille"][w]
             for f in FRAMES_REF])) for w in range(4)],
    }
    rec["mean_4frame_reference"]["ratio_vs_photo"] = [
        float(x / p) for x, p in zip(
            rec["mean_4frame_reference"]["depth_permille"], photo_depth)]
    return rec


def stage1(s0):
    t0 = time.time()
    log("=" * 78)
    log("[stage 1] R = 10 fresh seeds x frames %s" % FRAMES_REF)
    ctx = R12.make_context()
    S = float(ctx["S"])
    depth_m = float(ctx["env"]["depthM"])
    photo_depth = [float(x) for x in s0["v412_control_reference"]
                   ["photo_depth_permille"]]
    prod_hi, prod_lo = prod_rms_targets(ctx)

    per = {}
    for i, seed in enumerate(SEEDS_10, start=1):
        rid = "V%d" % i
        log("  %-4s seed=%d" % (rid, seed))
        per[rid] = measure_realisation(ctx, rid, seed, prod_hi,
                                       photo_depth, S, depth_m)

    photo_prom = None
    if "photo" in ctx:
        ph = R12.prom2d_block(np.asarray(ctx["photo"]["u8"][:, :, 1],
                                         np.float64))
        photo_prom = {
            "smooth": float(ph["window_values"]["smooth"]["prom2d"]),
            "bark": float(ph["window_values"]["bark"]["prom2d"]),
            "tiles_24_median": float(ph["tiles_24_median"]),
            "tiles_24_max": float(ph["tiles_24_max"]),
            "definition": ph["definition"],
        }

    out = {
        "at": now_iso(),
        "frames": FRAMES_REF,
        "validation_frames": VALID_FRAMES,
        "reference_only_frames": [24],
        "n_realisations": len(per),
        "window_order": WINDOW_ORDER,
        "photo_depth_permille": photo_depth,
        "production_slope_rms": {
            str(f): {"prelowpass": float(prod_hi[f]),
                     "lowpassed": float(prod_lo[f]),
                     "pre_over_low": float(prod_hi[f] / prod_lo[f])}
            for f in FRAMES_REF},
        "photo_prom2d": photo_prom,
        "per_realisation": per,
        "wall_s": float(time.time() - t0),
    }
    log("  photo prom2d: %s" % (None if photo_prom is None else
                                round(photo_prom["tiles_24_median"], 2)))
    del ctx
    gc.collect()
    cp_extra("stage1", {
        "at": out["at"], "frames": FRAMES_REF,
        "n_realisations": out["n_realisations"],
        "seeds": {rid: per[rid]["seed"] for rid in per},
        "effective_n": out["n_realisations"],
        "wall_s": out["wall_s"]})
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# STAGE 2 -- aggregation (realisation unit) + decision rules + expansion
# ---------------------------------------------------------------------------
def measure_extra(ctx, per, seed_positions, gen15, prod_hi, photo_depth,
                  S, depth_m, tag):
    """Append +N realisations from the pre-registered RNG continuation."""
    added = []
    for pos in seed_positions:
        seed = int(gen15[pos - 1])          # position is 1-based
        rid = "V%d" % pos
        log("  %-4s seed=%d (+%s)" % (rid, seed, tag))
        per[rid] = measure_realisation(ctx, rid, seed, prod_hi,
                                       photo_depth, S, depth_m)
        added.append({"realisation": rid, "seed": seed, "position": pos})
    return added


def aggregate(per, s0, rids, photo_prom):
    photo_depth = [float(x) for x in s0["v412_control_reference"]
                   ["photo_depth_permille"]]
    # 1.12 stored the placement sd on its own stage0; rebuild from record
    v412 = read_json(V412_NUMBERS)
    v412_sd = {w: float(v412["stage0"]["photo_placement_sd_permille"][w])
               for w in WINDOW_ORDER}
    placement_sd = [v412_sd[w] for w in WINDOW_ORDER]

    dist = {}
    dist_4 = {}
    dist_sn = {}
    for w, name in enumerate(WINDOW_ORDER):
        vals = [per[r]["mean_validation"]["ratio_vs_photo"][w] for r in rids]
        st = summarize(vals, BOOT_SEED + 101 * w)
        mean_depth = float(np.mean(
            [per[r]["mean_validation"]["depth_permille"][w] for r in rids]))
        excess = mean_depth - photo_depth[w]
        st.update({
            "window": name,
            "photo_depth_permille": photo_depth[w],
            "mean_depth_permille": mean_depth,
            "mean_excess_permille": float(excess),
            "placement_sd_permille": placement_sd[w],
            "placement_sd_ratio": float(placement_sd[w] / photo_depth[w]),
            "excess_over_placement_sd": float(excess / placement_sd[w]),
            "frames_used": VALID_FRAMES,
            "threshold_mean_reference": THRESH_MEAN,
            "threshold_q90_reference": THRESH_Q90,
            "mean_minus_reference_mean": float(st["mean"] - THRESH_MEAN),
            "q90_minus_reference_q90": float(st["q90"] - THRESH_Q90),
            "ci90_mean_straddles_1_05": bool(
                st["ci90_mean"][0] <= THRESH_MEAN <= st["ci90_mean"][1]),
        })
        dist[name] = st

        vals4 = [per[r]["mean_4frame_reference"]["ratio_vs_photo"][w]
                 for r in rids]
        st4 = summarize(vals4, BOOT_SEED + 101 * w + 3)
        st4["window"] = name
        st4["note"] = ("4-frame mean incl. frame 24: continuity only, "
                       "NOT a decision input")
        dist_4[name] = st4

        valsn = [per[r]["mean_validation_sn"]["ratio_vs_photo"][w]
                 for r in rids]
        stn = summarize(valsn, BOOT_SEED + 101 * w + 7)
        mean_depth_n = float(np.mean(
            [per[r]["mean_validation_sn"]["depth_permille"][w]
             for r in rids]))
        stn.update({
            "window": name, "photo_depth_permille": photo_depth[w],
            "mean_depth_permille": mean_depth_n,
            "mean_excess_permille": float(mean_depth_n - photo_depth[w]),
            "placement_sd_permille": placement_sd[w],
            "excess_over_placement_sd": float(
                (mean_depth_n - photo_depth[w]) / placement_sd[w]),
            "role": "secondary display only",
        })
        dist_sn[name] = stn

    # ---- frame correlation (within-window centred pool, as in 1.12) -----
    corr = {}
    for (f1, f2) in FRAME_PAIRS:
        key = "%d-%d" % (f1, f2)
        per_window = {}
        cen_x, cen_y = [], []
        raw_x, raw_y = [], []
        for w, name in enumerate(WINDOW_ORDER):
            xs = [per[r]["per_frame"][str(f1)]["depth_permille"][w]
                  for r in rids]
            ys = [per[r]["per_frame"][str(f2)]["depth_permille"][w]
                  for r in rids]
            per_window[name] = {
                "rho": pearson(xs, ys), "spearman": spearman(xs, ys),
                "n": len(rids),
                "pairs": [[float(a), float(b)] for a, b in zip(xs, ys)]}
            m = float(np.mean(np.asarray(xs + ys, np.float64)))
            cen_x.extend([float(v - m) for v in xs])
            cen_y.extend([float(v - m) for v in ys])
            raw_x.extend(xs)
            raw_y.extend(ys)
        corr[key] = {
            "per_window": per_window,
            "pooled_within_window_centered": {
                "rho": pearson(cen_x, cen_y),
                "spearman": spearman(cen_x, cen_y),
                "n_pairs": len(cen_x),
                "note": "each window centred on its own mean for this pair"},
            "fisher_mean_of_4_windows": fisher_mean(
                [per_window[n]["rho"] for n in WINDOW_ORDER]),
            "raw_pool_reference": {
                "rho": pearson(raw_x, raw_y), "n_pairs": len(raw_x),
                "note": ("contains the between-window level difference "
                         "(smooth ~61 vs bark ~42 permille); reference only")},
        }

    # ---- prom2d (G-LINE2 axis) vs depth ---------------------------------
    prom = {"per_realisation": {}, "correlation": {},
            "photo": photo_prom, "notes": {}}
    for rid in rids:
        blk = per[rid]["frame24"]["prom2d"]
        prom["per_realisation"][rid] = {
            "seed": per[rid]["seed"],
            "smooth": float(blk["window_values"]["smooth"]["prom2d"]),
            "bark": float(blk["window_values"]["bark"]["prom2d"]),
            "tiles_24_median": float(blk["tiles_24_median"]),
            "tiles_24_max": float(blk["tiles_24_max"]),
        }
    for w, name in enumerate(WINDOW_ORDER):
        typ = WINDOW_TYPE[w]
        xs = [prom["per_realisation"][r][typ] for r in rids]
        ys24 = [per[r]["frame24"]["depth_permille"][w] for r in rids]
        ys3 = [per[r]["mean_validation"]["depth_permille"][w] for r in rids]
        prom["correlation"][name] = {
            "prom2d_window": typ,
            "prom2d": [float(x) for x in xs],
            "depth_frame24_permille": [float(y) for y in ys24],
            "depth_validation_mean_permille": [float(y) for y in ys3],
            "rho_pearson_frame24": pearson(xs, ys24),
            "rho_spearman_frame24": spearman(xs, ys24),
            "rho_pearson_validation_mean": pearson(xs, ys3),
            "n": len(rids),
            "note": ("n = R realisations; smooth x/y share one prom2d and "
                     "bark x/y share one (2-D ring is axis-agnostic)"),
        }
    all_x, all_y = [], []
    ded_x, ded_y = [], []
    for w, name in enumerate(WINDOW_ORDER):
        all_x.extend(prom["correlation"][name]["prom2d"])
        all_y.extend(prom["correlation"][name]["depth_frame24_permille"])
    for r in rids:
        for typ in ("smooth", "bark"):
            ded_x.append(prom["per_realisation"][r][typ])
            ded_y.append(float(np.mean(
                [per[r]["frame24"]["depth_permille"][w]
                 for w, n2 in enumerate(WINDOW_ORDER)
                 if WINDOW_TYPE[w] == typ])))
    prom["correlation"]["pooled_all_pairs"] = {
        "rho_pearson": pearson(all_x, all_y),
        "rho_spearman": spearman(all_x, all_y),
        "n_pairs": len(all_x),
        "note": "4 windows x R; x duplicated inside each window type"}
    prom["correlation"]["pooled_unique_2xR"] = {
        "rho_pearson": pearson(ded_x, ded_y),
        "rho_spearman": spearman(ded_x, ded_y),
        "n_pairs": len(ded_x),
        "note": "2 window types x R; depth = mean of the two axis depths"}
    prom["correlation"]["fisher_mean_of_4_window_rho"] = fisher_mean(
        [prom["correlation"][n]["rho_pearson_frame24"]
         for n in WINDOW_ORDER])
    prom["notes"] = {
        "definition": ("2-D FFT of (gray - gaussian(sigma=3.0)) * Hann, "
                       "256x256, power |F|^2, ring radius 256/15..256/5, "
                       "prom2d = max/median"),
        "claim": "measured correlation only; no independence judgement",
    }

    # ---- decision rules --------------------------------------------------
    v412_means = {w: float(v412["stage2"]["distribution"][w]["mean"])
                  for w in WINDOW_ORDER}
    v412_ci = {w: [float(x) for x in
                   v412["stage2"]["distribution"][w]["ci90_mean"]]
               for w in WINDOW_ORDER}

    rule_i = {}
    for name in WINDOW_ORDER:
        new_mean = float(dist[name]["mean"])
        new_ci = dist[name]["ci90_mean"]
        ref_mean = float(v412_means[name])
        new_dir = "+" if new_mean > 1.0 else ("-" if new_mean < 1.0 else "0")
        ref_dir = "+" if ref_mean > 1.0 else ("-" if ref_mean < 1.0 else "0")
        contains = bool(new_ci[0] <= ref_mean <= new_ci[1])
        overlap = bool(not (new_ci[1] < v412_ci[name][0]
                            or new_ci[0] > v412_ci[name][1]))
        rule_i[name] = {
            "new_mean": new_mean, "new_ci90": new_ci,
            "new_direction": new_dir,
            "v412_mean_R5": ref_mean, "v412_ci90": v412_ci[name],
            "v412_direction": ref_dir,
            "direction_reproduced": bool(new_dir == ref_dir),
            "v412_mean_inside_new_ci": contains,
            "ci90_overlaps_v412_ci": overlap,
            "window_ok": bool(new_dir == ref_dir and (contains or overlap)),
            "v412_pre_reg_display": {"smooth x": 0.953, "smooth y": 0.999,
                                     "bark x": 0.932, "bark y": 1.056}[name],
        }
    rule_i_holds = bool(all(rule_i[n]["window_ok"] for n in WINDOW_ORDER))

    by = "bark y"
    rule_ii = {
        "window": by,
        "mean": float(dist[by]["mean"]),
        "ci90_mean": dist[by]["ci90_mean"],
        "lower_above_1": bool(dist[by]["ci90_mean"][0] > 1.0),
        "mean_above_1": bool(dist[by]["mean"] > 1.0),
        "excludes_1_upper_side": bool(
            dist[by]["ci90_mean"][0] > 1.0
            and dist[by]["ci90_mean"][1] > 1.0),
        "v412_mean_R5": float(v412_means[by]),
        "v412_ci90": v412_ci[by],
        "mean_excess_permille": float(dist[by]["mean_excess_permille"]),
        "excess_over_placement_sd": float(
            dist[by]["excess_over_placement_sd"]),
        "n_realisations_above_1": int(sum(
            1 for v in dist[by]["cluster_values"] if v > 1.0)),
        "n_realisations": int(len(dist[by]["cluster_values"])),
    }
    rule_ii_holds = bool(rule_ii["lower_above_1"])

    rule_iii = {}
    for name in WINDOW_ORDER:
        width = float(dist[name]["q90_ci_width"])
        rule_iii[name] = {"q90": float(dist[name]["q90"]),
                          "ci90_q90": dist[name]["ci90_q90"],
                          "width": width, "width_le_0_10": bool(width <= 0.10)}
    rule_iii_holds = bool(all(rule_iii[n]["width_le_0_10"]
                              for n in WINDOW_ORDER))

    unresolved = [n for n in WINDOW_ORDER
                  if dist[n]["ci90_mean_straddles_1_05"]]

    decision = {
        "rule_i": {"holds": rule_i_holds, "per_window": rule_i,
                   "test": ("same direction as 1.12 on all 4 windows AND "
                            "1.12 mean inside the new 90 % CI (or CIs "
                            "overlap)")},
        "rule_ii": {"holds": rule_ii_holds, "detail": rule_ii,
                    "test": "bark y 90 % CI lower bound > 1.0"},
        "rule_iii": {"holds": rule_iii_holds, "per_window": rule_iii,
                     "test": "q90 90 % CI width <= 0.10 on all 4 windows"},
        "both_i_and_ii": bool(rule_i_holds and rule_ii_holds),
        "threshold_status": ("provisional placement values; this job does "
                             "not upgrade or lower them"),
    }

    bark_y = {
        "window": by,
        "photo_depth_permille": photo_depth[WINDOW_ORDER.index(by)],
        "mean_ratio_validation_frames": float(dist[by]["mean"]),
        "sd_ratio": float(dist[by]["sd"]),
        "se": float(dist[by]["se"]),
        "ci90_mean": dist[by]["ci90_mean"],
        "q90": float(dist[by]["q90"]),
        "ci90_q90": dist[by]["ci90_q90"],
        "per_realisation_mean_ratio": {
            r: float(per[r]["mean_validation"]["ratio_vs_photo"]
                     [WINDOW_ORDER.index(by)]) for r in rids},
        "per_frame_mean_ratio": {
            str(f): float(np.mean(
                [per[r]["per_frame"][str(f)]["ratio_vs_photo"]
                 [WINDOW_ORDER.index(by)] for r in rids]))
            for f in FRAMES_REF},
        "n_cells_above_1": int(sum(
            1 for r in rids for f in VALID_FRAMES
            if per[r]["per_frame"][str(f)]["ratio_vs_photo"]
            [WINDOW_ORDER.index(by)] > 1.0)),
        "n_cells": int(len(rids) * len(VALID_FRAMES)),
        "mean_excess_permille": float(dist[by]["mean_excess_permille"]),
        "excess_over_placement_sd": float(dist[by]["excess_over_placement_sd"]),
        "frame24_mean_ratio_raw": float(np.mean(
            [per[r]["frame24"]["ratio_vs_photo"][WINDOW_ORDER.index(by)]
             for r in rids])),
        "frame24_mean_ratio_sn": float(np.mean(
            [per[r]["frame24"]["sn_ratio_vs_photo"]
             [WINDOW_ORDER.index(by)] for r in rids])),
        "v412_frame24_mean_ratio_raw": 1.110,
        "v412_frame24_mean_ratio_sn": 1.120,
    }

    frame_wise = {}
    for f in FRAMES_REF:
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

    return {
        "at": now_iso(),
        "n_realisations": len(rids),
        "realisation_ids": rids,
        "seeds": {r: int(per[r]["seed"]) for r in rids},
        "frames_validation": VALID_FRAMES,
        "frames_reference_only": [24],
        "photo_depth_permille": photo_depth,
        "photo_placement_sd_permille": placement_sd,
        "photo_placement_sd_ratio": [
            float(placement_sd[w] / photo_depth[w])
            for w in range(4)],
        "photo_placement_sd_ratio_brief": [0.035, 0.039, 0.099, 0.061],
        "distribution": dist,
        "distribution_4frame_reference": dist_4,
        "distribution_sn": dist_sn,
        "frame_corr": corr,
        "frame_wise": frame_wise,
        "prom2d": prom,
        "bark_y": bark_y,
        "decision": decision,
        "unresolved_windows": unresolved,
        "bootstrap": {"n_boot": N_BOOT, "seed": BOOT_SEED,
                      "unit": "realisation (cluster)",
                      "percentiles": [5.0, 95.0]},
        "definitions": {
            "representative_value": ("mean of the raw ratios at frames "
                                     "0/16/32 inside one realisation"),
            "distribution": ("mean / sd / SE / 90 % CI over the R "
                             "realisation values; effective n = R"),
            "q90": ("90th percentile of the R realisation values "
                    "(reference statistic); CI by the same cluster "
                    "bootstrap"),
            "rho": ("across realisations, matched within realisation; "
                    "pooled with each window centred on its own mean"),
            "sn": ("slope RMS rescaled to the production pre-lowpass "
                   "target, 1.10's convention; secondary display only"),
            "prom2d": ("G-LINE2 power-ring max/median (24 tile median) at "
                       "frame 24"),
        },
    }


def stage2(s0, s1):
    t0 = time.time()
    log("=" * 78)
    log("[stage 2] realisation-level aggregation + decision rules")
    per = s1["per_realisation"]
    rids = ["V%d" % i for i in range(1, len(per) + 1)]
    photo_prom = s1["photo_prom2d"]

    agg = aggregate(per, s0, rids, photo_prom)
    first_pass = {
        "n_realisations": agg["n_realisations"],
        "distribution": {
            n: {k: agg["distribution"][n][k] for k in (
                "mean", "sd", "se", "ci90_mean", "q90", "ci90_q90",
                "q90_ci_width", "ci90_mean_straddles_1_05")}
            for n in WINDOW_ORDER},
        "decision": {
            "rule_i_holds": agg["decision"]["rule_i"]["holds"],
            "rule_ii_holds": agg["decision"]["rule_ii"]["holds"],
            "rule_iii_holds": agg["decision"]["rule_iii"]["holds"],
            "unresolved_windows": list(agg["unresolved_windows"]),
        },
        "note": ("pre-registered primary n (R = 10) before the unresolved "
                 "rule was applied; the decisions in this report use the "
                 "final R after expansion"),
    }
    log("  distribution (mean | sd | 90% CI | q90 width):")
    for name in WINDOW_ORDER:
        d = agg["distribution"][name]
        log("    %-9s %.3f | %.3f | [%.3f, %.3f] | %.3f"
            % (name, d["mean"], d["sd"], d["ci90_mean"][0],
               d["ci90_mean"][1], d["q90_ci_width"]))
    log("  unresolved (CI straddles 1.05): %s" % (agg["unresolved_windows"] or "none"))

    expansion = []
    if agg["unresolved_windows"]:
        log("  pre-registered expansion triggered (+3, cap R = 15)")
        ctx = R12.make_context()
        S = float(ctx["S"])
        depth_m = float(ctx["env"]["depthM"])
        photo_depth = [float(x) for x in s0["v412_control_reference"]
                       ["photo_depth_permille"]]
        prod_hi, _prod_lo = prod_rms_targets(ctx)
        gen15 = s0["seed_generator"]["generated_15"]
        while agg["unresolved_windows"] and len(rids) < 15:
            n_add = min(3, 15 - len(rids))
            positions = list(range(len(rids) + 1, len(rids) + 1 + n_add))
            added = measure_extra(ctx, per, positions, gen15, prod_hi,
                                  photo_depth, S, depth_m,
                                  "expansion")
            rids = ["V%d" % i for i in range(1, len(per) + 1)]
            agg = aggregate(per, s0, rids, photo_prom)
            expansion.append({
                "added": added, "n_realisations": len(rids),
                "unresolved_after": agg["unresolved_windows"],
                "resolved": not agg["unresolved_windows"],
                "distribution_mean": {n: agg["distribution"][n]["mean"]
                                      for n in WINDOW_ORDER},
            })
            cp_extra("stage1_expansion_%d" % len(expansion), {
                "at": now_iso(), "n_realisations": len(rids),
                "added": added,
                "unresolved_after": agg["unresolved_windows"]})
            log("    -> R = %d, unresolved: %s"
                % (len(rids), agg["unresolved_windows"] or "none"))
        del ctx
        gc.collect()

    agg["first_pass_R10"] = first_pass
    agg["expansion_history"] = expansion
    agg["expansion_cap_reached"] = bool(
        agg["unresolved_windows"] and len(rids) >= 15)
    agg["wall_s"] = float(time.time() - t0)

    cp_extra("stage2", {
        "at": agg["at"],
        "n_realisations": agg["n_realisations"],
        "distribution": {k: {kk: vv for kk, vv in v.items()
                             if kk not in ("cluster_values",)}
                         for k, v in agg["distribution"].items()},
        "decision": {
            "rule_i_holds": agg["decision"]["rule_i"]["holds"],
            "rule_ii_holds": agg["decision"]["rule_ii"]["holds"],
            "rule_iii_holds": agg["decision"]["rule_iii"]["holds"],
        },
        "unresolved_windows": agg["unresolved_windows"],
        "expansion_history": expansion,
        "wall_s": agg["wall_s"],
    })
    log("  wall %.1fs" % (time.time() - t0))
    return agg


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def f3(x):
    return "%.3f" % float(x)


def f2(x):
    return "%.2f" % float(x)


def f4(x):
    return "%.4f" % float(x)


def ci(c):
    return "[%s, %s]" % (f3(c[0]), f3(c[1]))


def rho_txt(r):
    return "—" if r is None else "%.3f" % float(r)


def write_report(s0, s1, s2, s3):
    t0 = time.time()
    log("=" * 78)
    log("[report] REPORT-V4-1-13.md + numbers-v4-1-13.json")
    L = []
    A = L.append
    W = WINDOW_ORDER
    D = s2["distribution"]
    dec = s2["decision"]
    ctl = s0["control"]

    A("# REPORT-V4-1-13 — G-LINE-stat 閾値検証（新規独立 seed R=10 → 未解決により R=%d）"
      % s2["n_realisations"])
    A("")
    A("Brief: `BRIEF-V4-1-13.md`（`%s`）  " % s0["brief"])
    A("事前登録: `PRE-REG-GLINE-STAT.md`（`%s`、2026-11-14）  " % s0["prereg"])
    A("Workdir: `%s`  " % HERE)
    A("生成: %s" % s3["built_utc"])
    A("新規ファイルは `run_v4_1_13.py` の 1 本のみ。本ドライバは既存 .py を読み取り専用で開き、"
      "書き込み先は `out/v4-1-13/` の 3 ファイルだけ。")
    A("")
    A("**これは閾値の検証であって、モデルの合否判定ではない。合否は出さない。**  "
      "閾値 1.05（平均）/ 1.25（q90）の数値は動かさない。格上げの決定は本人が行う。")
    A("")
    A("---")
    A("")
    A("## 0. 要約（数値のみ）")
    A("")
    A("新規独立 seed 10 本（`default_rng(%d).integers(0, 2**63, size=10)` の V1–V10、"
      "事前登録の列をこの順で使用、代用なし）× 検証フレーム {0, 16, 32}。"
      "各 realization の代表値は 3 フレームの生比の平均。分布は realization 単位（実効 n = R）。"
      % SEED_GEN)
    A("")
    if s2["expansion_history"]:
        A("R=10 の時点で bark y の 90%% CI が参照線 1.05 をまたいだため、事前登録の"
          "未解決規則に従って V11–V15 を追加し、最終 R = %d で集計した（追加 seed は"
          "同じ生成器の 11–15 番目）。判定 (i)(ii)(iii) は最終 R の値で行った。"
          % s2["n_realisations"])
    else:
        A("R = %d で未解決の窓はなく、追加 realization は行っていない。"
          % s2["n_realisations"])
    A("")
    A("| 窓 | realization 平均（生比） | SE | 90% CI | q90 |")
    A("| --- | --- | --- | --- | --- |")
    for name in W:
        d = D[name]
        A("| %s | %s | %s | %s | %s |"
          % (name, f3(d["mean"]), f3(d["se"]), ci(d["ci90_mean"]),
             f3(d["q90"])))
    A("")
    A("- 判定規則 (i) 方向の再現: **%s**" % ("成立" if dec["rule_i"]["holds"]
                                            else "不成立"))
    A("- 判定規則 (ii) bark y の系統性: **%s**"
      % ("成立" if dec["rule_ii"]["holds"] else "不成立"))
    A("- 判定規則 (iii) q90 の復帰条件: **%s**"
      % ("成立" if dec["rule_iii"]["holds"] else "不成立"))
    if s2["unresolved_windows"]:
        A("- 未解決の窓（CI が 1.05 をまたぐ）: %s（§4）"
          % " / ".join(s2["unresolved_windows"]))
    else:
        A("- 未解決の窓（CI が 1.05 をまたぐ）: なし（§4）")
    A("")
    A("STAGE 0 の対照（1.12 本番 seed の frame 24 再現）は 4 窓 depth の最大差 "
      "%s、G-WARP 差 %+.3e px で完全一致した（§7.2）。" % ("0.0e+00"
                                                       if ctl["exact_equal"]
                                                       else "%.3e" % ctl["max_abs_diff"],
                                                       ctl["g_warp_diff"]))
    A("")
    A("---")
    A("")
    A("## 1. 数値表（窓ごと、realization 単位）")
    A("")
    A("### 1.1 主表 — 平均・SE・90% CI・q90・ρ")
    A("")
    A("| 窓 | 1.12 平均（R=5, 4f） | 本検証 平均（R=%d, 3f） | 差 | sd | SE | "
      "90%% CI | q90 | q90 の 90%% CI | q90 CI 幅 | ρ（窓内中心化プール 0–32） |"
      % s2["n_realisations"])
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for name in W:
        d = D[name]
        ref = s0["v412_control_reference"]["v412_means_R5_4frames"][name]
        A("| %s | %s | %s | %+.3f | %s | %s | %s | %s | %s | %s | %s |"
          % (name, f3(ref), f3(d["mean"]), d["mean"] - ref, f3(d["sd"]),
             f3(d["se"]), ci(d["ci90_mean"]), f3(d["q90"]),
             ci(d["ci90_q90"]), f3(d["q90_ci_width"]),
             rho_txt(s2["frame_corr"]["0-32"]["per_window"][name]["rho"])))
    A("")
    A("実効 n = R = %d。CI は realization を再抽出するクラスタ・ブートストラップ"
      "（%d 回、5/95 パーセンタイル、seed %d）。フレームは独立標本として数えていない。"
      % (s2["n_realisations"], N_BOOT, BOOT_SEED))
    A("")
    A("配置感度（1.10 §9.10 の窓原点シフト sd）との併記。配置感度は写真の分布ではない。")
    A("")
    A("| 窓 | 写真 depth [‰] | 平均 depth [‰] | 平均超過 [‰] | 配置感度 sd [‰] | 超過 / 配置感度 sd |")
    A("| --- | --- | --- | --- | --- | --- |")
    for w, name in enumerate(W):
        d = D[name]
        A("| %s | %s | %s | %+.2f | %s | %+.2f |"
          % (name, f2(d["photo_depth_permille"]),
             f2(d["mean_depth_permille"]), d["mean_excess_permille"],
             f2(d["placement_sd_permille"]), d["excess_over_placement_sd"]))
    A("")
    A("### 1.2 realization 別の代表値（生比、3 フレーム平均）")
    A("")
    A("V1–V10 は事前登録の検証標本。V11–V15 は未解決規則による追加分（生成器の 11–15 番目）。")
    A("")
    A("| realization | seed | " + " | ".join(W) + " |")
    A("| --- | --- | " + " | ".join(["---"] * 4) + " |")
    for rid in s2["realisation_ids"]:
        row = [f3(s2["distribution"][n]["cluster_values"][
            s2["realisation_ids"].index(rid)]) for n in W]
        A("| %s | %d | %s |" % (rid, s2["seeds"][rid], " | ".join(row)))
    A("")
    A("### 1.3 フレーム別の平均比（realization %d 本の平均）" % s2["n_realisations"])
    A("")
    A("| 窓 | frame 0 | frame 16 | frame 32 | frame 24（参考） | 3 フレーム平均 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for name in W:
        fw = s2["frame_wise"]
        A("| %s | %s | %s | %s | %s | %s |"
          % (name, f3(fw["0"][name]["mean_ratio"]),
             f3(fw["16"][name]["mean_ratio"]),
             f3(fw["32"][name]["mean_ratio"]),
             f3(fw["24"][name]["mean_ratio"]),
             f3(D[name]["mean"])))
    A("")
    A("frame 24 は参考記録のみで、判定には使っていない。4 フレーム平均（frame 24 を含む）は"
      "連続性の確認用に numbers へ別途記録した。")
    A("")
    A("### 1.4 フレーム間相関 ρ（実測、n = R = %d）" % s2["n_realisations"])
    A("")
    A("| 集計 | 0–16 | 16–32 | 0–32 |")
    A("| --- | --- | --- | --- |")
    A("| 窓内中心化プール（主） | %s | %s | %s |"
      % tuple(rho_txt(s2["frame_corr"][k]
                      ["pooled_within_window_centered"]["rho"])
              for k in ("0-16", "16-32", "0-32")))
    A("| 4 窓の Fisher 平均 | %s | %s | %s |"
      % tuple(rho_txt(s2["frame_corr"][k]["fisher_mean_of_4_windows"])
              for k in ("0-16", "16-32", "0-32")))
    A("| 生のプール（窓間水準差を含む・参考） | %s | %s | %s |"
      % tuple(rho_txt(s2["frame_corr"][k]["raw_pool_reference"]["rho"])
              for k in ("0-16", "16-32", "0-32")))
    A("")
    A("| 窓 | 0–16 | 16–32 | 0–32 |")
    A("| --- | --- | --- | --- |")
    for name in W:
        A("| %s | %s | %s | %s |"
          % (name,
             rho_txt(s2["frame_corr"]["0-16"]["per_window"][name]["rho"]),
             rho_txt(s2["frame_corr"]["16-32"]["per_window"][name]["rho"]),
             rho_txt(s2["frame_corr"]["0-32"]["per_window"][name]["rho"])))
    A("")
    A("### 1.5 振幅正規化（Sn）版 — 副表示のみ")
    A("")
    A("| 窓 | 平均比（Sn） | sd | SE | 90% CI | 生値との差 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for name in W:
        d = s2["distribution_sn"][name]
        A("| %s | %s | %s | %s | %s | %+.3f |"
          % (name, f3(d["mean"]), f3(d["sd"]), f3(d["se"]),
             ci(d["ci90_mean"]), d["mean"] - D[name]["mean"]))
    A("")
    A("### 1.6 G-LINE2 軸 — prom2d（24 窓中央値）と depth の realization 間相関")
    A("")
    A("| realization | seed | prom2d（24 窓中央値） | 24 窓最大 |")
    A("| --- | --- | --- | --- |")
    for rid in s2["realisation_ids"]:
        pr = s2["prom2d"]["per_realisation"][rid]
        A("| %s | %d | %s | %s |"
          % (rid, pr["seed"], f2(pr["tiles_24_median"]),
             f2(pr["tiles_24_max"])))
    if s2["prom2d"]["photo"]:
        ph = s2["prom2d"]["photo"]
        A("| 写真 | — | %s | %s |" % (f2(ph["tiles_24_median"]),
                                    f2(ph["tiles_24_max"])))
    A("")
    A("| 窓 | Pearson r（frame 24 depth） | Spearman ρ | Pearson r（3 フレーム平均 depth） | n |")
    A("| --- | --- | --- | --- | --- |")
    for name in W:
        c = s2["prom2d"]["correlation"][name]
        A("| %s | %s | %s | %s | %d |"
          % (name, rho_txt(c["rho_pearson_frame24"]),
             rho_txt(c["rho_spearman_frame24"]),
             rho_txt(c["rho_pearson_validation_mean"]), c["n"]))
    A("")
    A("プール（4 窓 × R）r = %s、窓型 2 × R の重複除去プール r = %s、"
      "4 窓の Fisher 平均 r = %s。**独立性の判定はしない**（DESIGN §3）。"
      % (rho_txt(s2["prom2d"]["correlation"]["pooled_all_pairs"]
                 ["rho_pearson"]),
         rho_txt(s2["prom2d"]["correlation"]["pooled_unique_2xR"]
                 ["rho_pearson"]),
         rho_txt(s2["prom2d"]["correlation"]
                 ["fisher_mean_of_4_window_rho"])))
    A("")
    A("---")
    A("")
    A("## 2. 事前登録の判定規則 (i)(ii)(iii)")
    A("")
    A("### (i) 方向の再現 — **%s**"
      % ("成立" if dec["rule_i"]["holds"] else "不成立"))
    A("")
    A("規則: 4 窓すべてで 1.12 と同じ方向（1 を上回る／下回る）を示し、かつ 1.12 の平均が"
      "新規 %d 本の 90%% CI の中に入る（または CI が重なる）。1 窓でも反転したら不成立。"
      % s2["n_realisations"])
    A("")
    A("| 窓 | 本検証 平均 | 方向 | 1.12 平均 | 方向 | 方向一致 | 1.12 が新 CI 内 | "
      "CI 重複 | 窓の判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for name in W:
        r = dec["rule_i"]["per_window"][name]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (name, f3(r["new_mean"]), r["new_direction"],
             f3(r["v412_mean_R5"]), r["v412_direction"],
             "一致" if r["direction_reproduced"] else "**反転**",
             "入る" if r["v412_mean_inside_new_ci"] else "外れる",
             "あり" if r["ci90_overlaps_v412_ci"] else "なし",
             "成立" if r["window_ok"] else "**不成立**"))
    A("")
    A("### (ii) bark y の系統性 — **%s**"
      % ("成立" if dec["rule_ii"]["holds"] else "不成立"))
    A("")
    r2 = dec["rule_ii"]["detail"]
    A("規則: bark y の realization 平均の 90% CI が 1.0 を上側に外すか。")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| bark y 平均比 | %s |" % f3(r2["mean"]))
    A("| 90%% CI | %s |" % ci(r2["ci90_mean"]))
    A("| CI 下限 > 1.0 | %s |" % ("はい" if r2["lower_above_1"] else "いいえ"))
    A("| 1 超の realization | %d / %d |"
      % (r2["n_realisations_above_1"], r2["n_realisations"]))
    A("| 平均超過 | %+.2f ‰（配置感度の %+.2f 倍） |"
      % (r2["mean_excess_permille"], r2["excess_over_placement_sd"]))
    A("| 1.12 平均（R=5, 4f） / 90%% CI | %s / %s |"
      % (f3(r2["v412_mean_R5"]), ci(r2["v412_ci90"])))
    A("")
    A("### (iii) q90 の復帰条件 — **%s**"
      % ("成立" if dec["rule_iii"]["holds"] else "不成立"))
    A("")
    A("規則: q90 のクラスタ・ブートストラップ 90% CI の幅が 4 窓すべてで ≤ 0.10 の"
      "場合に限り、q90 を判定量への復帰候補として記録する（復帰の決定は本人）。")
    A("")
    A("| 窓 | q90 | 90% CI | 幅 | ≤ 0.10 |")
    A("| --- | --- | --- | --- | --- |")
    for name in W:
        r3 = dec["rule_iii"]["per_window"][name]
        A("| %s | %s | %s | %s | %s |"
          % (name, f3(r3["q90"]), ci(r3["ci90_q90"]), f3(r3["width"]),
             "はい" if r3["width_le_0_10"] else "**いいえ**"))
    A("")
    A("---")
    A("")
    A("## 3. 閾値 1.05 / 1.25 への結論")
    A("")
    if dec["both_i_and_ii"]:
        A("(i) と (ii) が両方成立した。したがって**運用値の候補としての証拠が揃った**、"
          "とだけ書く。運用値への格上げは本人の決定であり、本ジョブは格上げを宣言しない。")
    else:
        if not dec["rule_i"]["holds"]:
            A("**(i) が不成立である。** 1.12 のパターンは seed 族の運だった可能性が高く、"
              "閾値は暫定配置値のままとする。設計の見直し（ラリーへ戻す）を提案する。")
            A("")
            A("設計見直しの提案（決定ではない）: どの窓で方向が反転したか・どの窓で 1.12 平均が"
              "新 CI の外に出たかを §2 の表で特定し、(a) 窓の再選択、(b) フレーム位相の追加、"
              "(c) realization 数の拡大、のどれを先に検討するかをラリーで決める。")
        if not dec["rule_ii"]["holds"]:
            A("**(ii) が不成立である。** bark y の系統的超過は確認できず、1.05 の参照線は "
              "bark y を分ける線として支持されない。閾値は暫定配置値のままとする。")
        if dec["rule_i"]["holds"] and not dec["rule_ii"]["holds"]:
            A("")
            A("(i) は成立したが (ii) が不成立のため、運用値の候補としての証拠は揃っていない。")
        if not dec["rule_i"]["holds"] and dec["rule_ii"]["holds"]:
            A("")
            A("(ii) は成立したが (i) が不成立のため、運用値の候補としての証拠は揃っていない。")
    A("")
    A("q90 については §2(iii) の判定どおり%s。閾値 1.25 の数値は動かしていない。"
      % ("、復帰候補として記録する" if dec["rule_iii"]["holds"]
         else "、参考統計のままとする"))
    A("")
    A("---")
    A("")
    A("## 4. 未解決の窓と追加 realization")
    A("")
    A("規則: realization 平均の CI が 1.05 の参照線をまたぐ窓があれば「未解決」とし、"
      "realization を +3 ずつ追加して再測定する（上限 R = 15）。")
    A("")
    if not s2["unresolved_windows"]:
        A("R = %d の時点で、4 窓すべての 90%% CI が 1.05 をまたがなかった。"
          "追加 realization は行っていない。" % s2["n_realisations"])
        A("")
        A("| 窓 | 90% CI | 1.05 をまたぐ |")
        A("| --- | --- | --- |")
        for name in W:
            A("| %s | %s | %s |" % (name, ci(D[name]["ci90_mean"]),
                                    "はい" if D[name]
                                    ["ci90_mean_straddles_1_05"] else "いいえ"))
    else:
        fp = s2["first_pass_R10"]
        A("R = %d の時点で未解決: %s。追加 realization を実行した。"
          % (fp["n_realisations"],
             " / ".join(fp["decision"]["unresolved_windows"])))
        A("")
        A("| 窓 | R=10 平均 | R=10 の 90%% CI | 1.05 をまたぐ | R=%d 平均 | R=%d の 90%% CI |"
          % (s2["n_realisations"], s2["n_realisations"]))
        A("| --- | --- | --- | --- | --- | --- |")
        for name in W:
            f10 = fp["distribution"][name]
            A("| %s | %s | %s | %s | %s | %s |"
              % (name, f3(f10["mean"]), ci(f10["ci90_mean"]),
                 "はい" if f10["ci90_mean_straddles_1_05"] else "いいえ",
                 f3(D[name]["mean"]), ci(D[name]["ci90_mean"])))
        A("")
        A("R=10 の判定: (i) %s / (ii) %s / (iii) %s。"
          % ("成立" if fp["decision"]["rule_i_holds"] else "不成立",
             "成立" if fp["decision"]["rule_ii_holds"] else "不成立",
             "成立" if fp["decision"]["rule_iii_holds"] else "不成立"))
        A("")
        for i, ex in enumerate(s2["expansion_history"], start=1):
            A("- 追加 %d: seeds %s → R = %d、未解決: %s"
              % (i, [a["seed"] for a in ex["added"]], ex["n_realisations"],
                 " / ".join(ex["unresolved_after"]) or "なし"))
        A("")
        if s2["expansion_cap_reached"]:
            A("**上限 R = 15 でも未解決が残った。** その旨をここに明記する。")
        A("")
        A("| 窓 | 90%% CI（最終 R=%d） | 1.05 をまたぐ |"
          % s2["n_realisations"])
        A("| --- | --- | --- |")
        for name in W:
            A("| %s | %s | %s |" % (name, ci(D[name]["ci90_mean"]),
                                    "はい" if D[name]
                                    ["ci90_mean_straddles_1_05"] else "いいえ"))
    A("")
    A("---")
    A("")
    A("## 5. 4 窓 AND の単一合否")
    A("")
    A("出さない（DESIGN §2 / 事前登録）。窓ごとの推定は §1 のとおりそのまま報告する。"
      "本レポートはモデルの合否も、4 窓を束ねた単一の判定も書かない。")
    A("")
    A("---")
    A("")
    A("## 6. 開示")
    A("")
    A("1. **prom2d は再実装。** 事前登録の定義（256²、`gray - gaussian(σ=3.0)`、Hann、"
      "2-D FFT、リング 5–15 px、max/median）を `run_v4_1_12.py` からそのまま使った。"
      "写真の 24 窓中央値は本実行で %s（記録 21.1）。G-LINE2 の元実装そのものではない。"
      % ("—" if not s2["prom2d"]["photo"]
         else f2(s2["prom2d"]["photo"]["tiles_24_median"])))
    A("2. **R は薄い。** 実効 n = R = %d。sd・CI・q90 は %d 個の realization 値から作った。"
      "CI の幅は「この標本が許す幅」であって、モデルの realization 分布の幅ではない。"
      % (s2["n_realisations"], s2["n_realisations"]))
    A("3. **seed の選び方。** 事前登録の `default_rng(%d).integers(0, 2**63, size=10)` の "
      "V1–V10 をこの順で使った（生成列が事前登録の 10 本と一致することを STAGE 0 で確認）。"
      "1.12 の 5 本（本番 +1/+2/+3/+7）は流用していない。無作為抽出ではなく、"
      "事前登録で固定された列である。" % SEED_GEN)
    A("4. **Sn は副表示のみ。** 主対象は raw 比。Sn の振幅一致は 1.10 の慣例"
      "（低域通過前の production RMS を目標）で、B_sigma 後の実振幅は数 % ずれる。"
      "判定はすべて raw 比で行った。")
    A("5. **frame 24 の扱い。** 参考記録として測ったが、判定・代表値・CI・q90・ρ の"
      "いずれにも使っていない。検証は 0 / 16 / 32 の 3 フレームのみ。")
    A("6. **窓・帯・depth の定義は 1.10 / 1.12 と同一で、測定中に変更していない。** "
      "窓の位置・サイズ・帯域（5–15 px）・depth の定義は凍結。")
    A("7. **配置感度は写真の分布ではない。** 1.10 §9.10 の窓原点シフト sd "
      "（比の単位で 3.5% / 3.9% / 9.9% / 6.1%）を超過の併記にだけ使った。"
      "有意差の言い方はしていない。")
    A("8. **目の権限は変わらない**（DESIGN §6）。本検証は数値の記録である。")
    A("")
    A("---")
    A("")
    A("## 7. 再現手順・無変更の記録")
    A("")
    A("### 7.1 再現手順")
    A("")
    A("```")
    A("cd %s" % HERE)
    A("python run_v4_1_13.py --all      # 全段: stage0 -> stage1 -> stage2 -> report")
    A("python run_v4_1_13.py --report   # 記録済み numbers からレポートだけ再生成")
    A("```")
    A("")
    A("`--stage0` / `--stage1` / `--stage2` / `--report` は内部の実行順で、通常の再現は "
      "`--all` の 1 回で完結する（stage0 は preflight + 対照、stage1 は V1–V10 の測定、"
      "stage2 は集計・未解決規則と追加分 V11–V15 の測定、report は数値と本レポート）。"
      "実行時間は約 5 分。")
    A("")
    A("### 7.2 STAGE 0 の対照（1.12 本番 seed, frame 24）")
    A("")
    A("| 窓 | 本ジョブ | numbers-v4-1-12 | 差 |")
    A("| --- | --- | --- | --- |")
    for w, name in enumerate(W):
        A("| %s | %.6f | %.6f | %+.1e |"
          % (name, ctl["depth_permille"][w],
             ctl["v412_depth_permille"][w], ctl["diff_vs_v412"][w]))
    A("")
    A("- 完全一致: **%s**（最大差 %.1e）。G-WARP %.9f px、1.12 記録 %.9f px、差 %+.1e。"
      % (ctl["exact_equal"], ctl["max_abs_diff"], ctl["g_warp_mean_px"],
         ctl["v412_g_warp_mean_px"], ctl["g_warp_diff"]))
    A("- 画素 sha256 の一致: **%s**。フック付きコピーの uint8 最大差 %d。"
      % (ctl["image_sha_equal"], ctl["copy_max_uint8_diff"]))
    A("- 対照は検証標本（R=%d）に含めていない。" % s2["n_realisations"])
    A("")
    A("### 7.3 G-LEGACY（既存ファイルの非改変）")
    A("")
    A("| 確認 | 数値 |")
    A("| --- | --- |")
    A("| .py の changed | %d |" % s3["py_legacy"]["n_changed"])
    A("| .py の added | %d |" % s3["py_legacy"]["n_added"])
    A("| .py の removed（missing） | %d |" % s3["py_legacy"]["n_missing"])
    A("| .py 本数（開始 / 終了） | %d / %d |"
      % (s3["py_legacy"]["n_start"], s3["py_legacy"]["n_end"]))
    A("| 主要参照ファイルの sha256 変化 | %s |"
      % ("なし" if s3["major_unchanged"] else "**あり**"))
    A("| 本ジョブが書き込みで開いた .py | 0 本 |")
    A("")
    if s3["py_legacy"]["changed"]:
        A("変化した .py（本ジョブ外の並行ジョブ）: %s"
          % ", ".join(s3["py_legacy"]["changed"]))
        A("")
    A("### 7.4 主要ファイルの sha256（開始 → 終了）")
    A("")
    A("| ファイル | 開始（先頭 16） | 終了（先頭 16） | 変化 |")
    A("| --- | --- | --- | --- |")
    for p, _w in REF_PATHS:
        a = s0["major_hashes_start"][p]["sha256"][:16]
        b = s3["major_hashes_end"][p]["sha256"][:16]
        A("| `%s` | `%s` | `%s` | %s |"
          % (os.path.basename(p), a, b, "なし" if a == b else "**あり**"))
    A("")
    A("1.12 の成果物 3 点（開始時の記録）: ")
    for k, v in s0["v412_files_start"].items():
        A("- `%s`: %d bytes, `%s`" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    A("**作業開始時の事故と復元（開示）**: 事前確認の途中で `run_v4_1_12.py` を誤って"
      "起動し、`out/v4-1-12/checkpoint.json` が書き換わった。`numbers-v4-1-12.json` が"
      "保持する stage0 / stage1 の元記録から復元し、復元後の bytes 32828 と sha256 "
      "`6b5a5a65065ff7db…` が 1.12 実行ログの記録（32828 / 6b5a5a65065ff7db）と一致する"
      "ことを確認した。1.12 の numbers とレポートは無改変（§7.4 の開始ハッシュ）。")
    A("")
    A("---")
    A("")
    A("## 8. 成果物")
    A("")
    A("| ファイル | 役割 |")
    A("| --- | --- |")
    A("| `REPORT-V4-1-13.md` | 本レポート |")
    A("| `numbers-v4-1-13.json` | 数値の全量（realization 別・分布・相関・prom2d・判定・ハッシュ） |")
    A("| `checkpoint.json` | jobkit の段階記録（stage0 / stage1 / stage2 / report） |")
    A("")
    A("画像は保存していない（測定はメモリ上で完結）。")
    A("")
    text = "\n".join(L) + "\n"
    os.makedirs(OUT, exist_ok=True)
    with open(REPORT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    log("  %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    log("  wall %.1fs" % (time.time() - t0))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run_report(s0, s1, s2):
    major_end = {p: {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
                 for p, _w in REF_PATHS}
    py_end = py_set()
    py_end_mtimes = py_mtimes()
    legacy = what_changed(s0["python_hashes_start"], py_end,
                          s0.get("python_mtimes_start"), py_end_mtimes)
    major_unchanged = all(
        major_end[p]["sha256"] == s0["major_hashes_start"][p]["sha256"]
        for p in major_end)

    v412_files_end = {}
    for p in (V412_NUMBERS, V412_REPORT, V412_CP):
        v412_files_end[os.path.basename(p)] = {
            "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    v412_unchanged = all(
        v412_files_end[k]["sha256"] == s0["v412_files_start"][k]["sha256"]
        for k in v412_files_end)

    s3 = {
        "built_utc": now_iso(),
        "major_hashes_start": s0["major_hashes_start"],
        "major_hashes_end": major_end,
        "major_unchanged": bool(major_unchanged),
        "v412_files_end": v412_files_end,
        "v412_unchanged": bool(v412_unchanged),
        "py_legacy": legacy,
        "python_files_n_end": len(py_end),
        "python_hashes_end": py_end,
        "python_mtimes_end": py_end_mtimes,
        "deliverables": {
            "numbers": "numbers-v4-1-13.json",
            "report": "REPORT-V4-1-13.md",
            "checkpoint": "checkpoint.json",
        },
    }

    numbers = {
        "brief": BRIEF,
        "prereg": PREREG,
        "workdir": HERE,
        "built_utc": s3["built_utc"],
        "new_files": NEW_FILES,
        "run_v4_1_13_sha256": sha256_file(os.path.abspath(__file__)),
        "scope": ("threshold validation only; no model pass/fail; 1.05/1.25 "
                  "unchanged provisional placement values"),
        "stage0": s0,
        "stage1": s1,
        "stage2": s2,
        "stage3": s3,
        "files": {"numbers": "numbers-v4-1-13.json",
                  "report": "REPORT-V4-1-13.md",
                  "checkpoint": "checkpoint.json"},
    }
    write_json(NUMBERS, numbers)
    write_report(s0, s1, s2, s3)

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
        "v412_unchanged": bool(v412_unchanged),
        "note": ("threshold validation only; no pass/fail; checkpoint's own "
                 "sha256 is not recorded (self-referential)"),
    })
    log("=" * 78)
    log("DONE")
    log("  %s" % REPORT)
    log("  %s" % NUMBERS)
    log("  %s" % CP_PATH)
    return {"report_bytes": os.path.getsize(REPORT),
            "numbers_bytes": os.path.getsize(NUMBERS)}


def main():
    stages = [a for a in sys.argv[1:] if a.startswith("--")]
    if "--all" in stages:
        stages = ["--stage0", "--stage1", "--stage2", "--report"]
    if not stages:
        stages = ["--stage0", "--stage1", "--stage2", "--report"]
    t0 = time.time()
    log("BRIEF-V4-1-13 -- G-LINE-stat threshold validation (R=10 fresh seeds)")

    s0 = s1 = s2 = None
    if "--stage0" in stages:
        s0 = stage0()
    if "--stage1" in stages:
        if s0 is None:
            s0 = read_json(NUMBERS)["stage0"]
        s1 = stage1(s0)
    if "--stage2" in stages:
        if s0 is None:
            s0 = read_json(NUMBERS)["stage0"]
        if s1 is None:
            s1 = read_json(NUMBERS)["stage1"]
        s2 = stage2(s0, s1)
    if "--report" in stages:
        if s0 is None:
            s0 = read_json(NUMBERS)["stage0"]
        if s1 is None:
            s1 = read_json(NUMBERS)["stage1"]
        if s2 is None:
            s2 = read_json(NUMBERS)["stage2"]
        run_report(s0, s1, s2)
    log("total wall %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
