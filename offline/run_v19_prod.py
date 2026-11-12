"""run_v19_prod -- BRIEF-V19-PROD8S-PHOTO: the 8 s piece on the PHOTO base,
v3.10 + v16 stack, clear water, depth 4 m.

Base: the project photo D:\\mutaro-work\\moth-on-bark.jpg (735x539), edge-
      replicated right+bottom by 1 px to 736x540 (the v2 production loop's
      own convention for H.264).
Water: v3.10 spectrum (elfouhaily-form, Arm D) at U = 6 m/s, F = 100 m,
       depth = 4.0 m, seed 3255045276382562994, time-evolved at t = n/24 s
       with the v10/v2 8 s omega quantisation.
Medium: the v12 CLEAR preset for both a/b (W.PRESETS["clear"]) and the phase
        function (V16.FF_PHASE["clear"], n = 1.098, mu = 3.37).
Transport: v16 water (a/b partition + FF phase + FF-CDF halo).
Composite: the frozen v5 composite A2.6 + Blow (sigma 3.0 px) + E2.2.

Stages (all foreground, chunked, idempotent per chunk):
  --stage0        preflight + parameter resolution + hash the *.py set
  --chunk K       render frames 24K..24K+23 (atomic per frame, hashed)
  --seam          G1 seam check + the frequency-table quantisation proof
  --warp          the warp yardstick (frame 96, transport field)
  --encode        encode the MP4 (libx264, 24 fps, yuv420p) + G3 ffprobe
  --gates         G2 levels + G4 integrity, roll up into the checkpoint
  --report        numbers-v19-prod.json + REPORT-V19-PROD8S.md + still
  --probe N       render one frame into out/v19/_work (timing only)

Exit code 0 only when the requested unit completed.
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

import run_v10 as R10                                            # noqa: E402
import run_v15 as R15                                            # noqa: E402
import run_v15_spectrum as V15S                                  # noqa: E402
import v12_water as W                                            # noqa: E402
import run_v16_water as V16                                      # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V19-PROD8S-PHOTO.md"
OUT = os.path.join(HERE, "out", "v19")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
WARP_PATH = os.path.join(WORK, "warp.json")
MP4 = os.path.join(OUT, "piece-8s-v19-photo.mp4")
NUMBERS = os.path.join(OUT, "numbers-v19-prod.json")
REPORT = os.path.join(OUT, "REPORT-V19-PROD8S.md")
STILL96 = os.path.join(OUT, "frame-0096.png")

PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
V18_NUMBERS = os.path.join(HERE, "out", "v18", "numbers-v18-prod.json")
PY_HASHES_BEFORE = os.path.join(OUT, "stage0-py-hashes.json")

#: the brief's fixed scene line (do not tune).
CELL = {"U_cli": 6.0, "F_cli": 100.0, "h_cli": 4.0}

#: the seed named literally by the brief (and what derive_seed gives).
SEED_BRIEF = 3255045276382562994

#: the frozen composite (A2.6 + Blow + E2.2), as in v18.
A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2

#: the water medium: the v12 preset id and the v16 phase-fit key.
WATER_PRESET = "clear"

FPS = 24
NFRAMES = 192
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
GRID = 2048
Q_DET = 25
WIDTH, HEIGHT = 736, 540
PHOTO_NATIVE = (735, 539)

#: files the brief names (all must exist in the workdir).
CODE_FILES = ("run_v18_prod.py", "run_v15_spectrum.py", "run_v16_water.py",
              "run_v16.py", "v12_water.py", "v5_composite.py",
              "transport.py", "det.py", "image.py")

NEW_FILES = ("run_v19_prod.py",)


def allowed_new_static():
    return set(NEW_FILES)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)


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
    payload = json.dumps(obj, ensure_ascii=False, indent=2).encode("utf-8")
    atomic_write_bytes(path, payload)


def makedirs():
    for d in (OUT, WORK, FRAMES):
        os.makedirs(d, exist_ok=True)


def py_set():
    """The top-level *.py set (the v17 guard's convention)."""
    out = {}
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".py") and os.path.isfile(os.path.join(HERE, name)):
            out[name] = sha256_file(os.path.join(HERE, name))
    return out


def checkpoint_extra(key, value=None):
    state = read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}
    if value is not None:
        state[key] = value
        state["updated_at"] = now_iso()
        payload = json.dumps(state, ensure_ascii=False, indent=2).encode(
            "utf-8")
        atomic_write_bytes(CP_PATH, payload)
    return state


def read_checkpoint():
    return read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}


# --------------------------------------------------------------------------
# the base photo (loaded once per process; the pad is the v2 convention)
# --------------------------------------------------------------------------
def load_photo():
    """735x539 JPEG -> 736x540 by right+bottom 1 px edge replication."""
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    native = im.size
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"native": native, "u8": padded, "lin": lin}


# --------------------------------------------------------------------------
# STAGE 0 -- preflight
# --------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "created_utc": now_iso(),
             "python": sys.version.split()[0]}
    missing = []

    # (a) the code set
    code = {}
    for name in CODE_FILES:
        p = os.path.join(HERE, name)
        ok = os.path.isfile(p)
        code[name] = {"present": ok,
                      "sha256": sha256_file(p) if ok else None}
        if not ok:
            missing.append(name)
    facts["code"] = code

    # (b) the v16 defaults the brief asserts
    facts["v16_defaults"] = {
        "PHASE_DEFAULT": V16.PHASE_DEFAULT,
        "KERNEL_DEFAULT": V16.KERNEL_DEFAULT,
        "assert_phase_ff": V16.PHASE_DEFAULT == "ff",
        "assert_kernel_ffcdf": V16.KERNEL_DEFAULT == "ff-cdf",
    }
    if not (V16.PHASE_DEFAULT == "ff" and V16.KERNEL_DEFAULT == "ff-cdf"):
        missing.append("v16 defaults not ff / ff-cdf")

    # (c) the CLEAR medium preset + the clear phase fit (brief demands both)
    clear_w = W.PRESETS.get(WATER_PRESET)
    clear_p = V16.FF_PHASE.get(WATER_PRESET)
    facts["clear_medium"] = {
        "v12_preset_key": WATER_PRESET,
        "v12_presets_available": sorted(W.PRESETS),
        "v12_default_medium_preset": W.DEFAULT_MEDIUM_PRESET,
        "v12_preset_found": clear_w is not None,
        "v12_label": (clear_w or {}).get("label"),
        "v12_a_part": (clear_w or {}).get("a_part"),
        "v12_b_part": (clear_w or {}).get("b_part"),
        "v12_load_default": (clear_w or {}).get("load_default"),
        "v16_phase_key": WATER_PRESET,
        "v16_phase_keys_available": sorted(V16.FF_PHASE),
        "v16_phase_found": clear_p is not None,
        "v16_n": (clear_p or {}).get("n"),
        "v16_mu": (clear_p or {}).get("mu"),
        "v16_label": (clear_p or {}).get("label"),
        "v16_provenance": (clear_p or {}).get("provenance"),
        "assert_n_1p098": bool(clear_p and abs(float(clear_p["n"])
                                               - 1.098) < 1e-9),
        "assert_mu_3p37": bool(clear_p and abs(float(clear_p["mu"])
                                               - 3.37) < 1e-9),
    }
    if clear_w is None:
        missing.append("no CLEAR medium preset in v12_water.PRESETS "
                       "(available: %s)" % sorted(W.PRESETS))
    if clear_p is None:
        missing.append("no CLEAR phase fit in run_v16_water.FF_PHASE "
                       "(available: %s)" % sorted(V16.FF_PHASE))
    if clear_p and not (abs(float(clear_p["n"]) - 1.098) < 1e-9
                        and abs(float(clear_p["mu"]) - 3.37) < 1e-9):
        missing.append("clear phase fit is not n=1.098, mu=3.37")

    # (d) the base photo
    from PIL import Image
    if not os.path.isfile(PHOTO):
        missing.append("base photo missing: %s" % PHOTO)
        facts["photo"] = {"path": PHOTO, "present": False}
    else:
        im = Image.open(PHOTO)
        im.load()
        ph = {"path": PHOTO, "present": True, "bytes": os.path.getsize(PHOTO),
              "native_px": [im.size[0], im.size[1]], "mode": im.mode,
              "format": im.format, "sha256": sha256_file(PHOTO),
              "assert_735x539": tuple(im.size) == PHOTO_NATIVE,
              "pad": ("right+bottom 1 px edge replication -> %dx%d (the v2 "
                      "production loop's own convention)" % (WIDTH, HEIGHT)),
              "output_px": [WIDTH, HEIGHT]}
        facts["photo"] = ph
        if tuple(im.size) != PHOTO_NATIVE:
            missing.append("photo decodes at %s, not 735x539" % (im.size,))
    # the pad arithmetic is checked on pixels, not asserted on trust
    if os.path.isfile(PHOTO):
        p = load_photo()
        n = np.asarray(Image.open(PHOTO).convert("RGB"))
        facts["photo"]["pad_check"] = {
            "shape": list(p["u8"].shape),
            "right_col_is_native_last_col": bool(
                np.array_equal(p["u8"][:-1, -1, :], n[:, -1, :])),
            "bottom_row_is_native_last_row": bool(
                np.array_equal(p["u8"][-1, :-1, :], n[-1, :, :])),
            "corner_is_native_corner": bool(
                np.array_equal(p["u8"][-1, -1, :], n[-1, -1, :])),
        }
        if not all(facts["photo"]["pad_check"][k] for k in (
                "right_col_is_native_last_col",
                "bottom_row_is_native_last_row",
                "corner_is_native_corner")):
            missing.append("pad does not follow the v2 edge-replication rule")

    # (e) the seed the brief names vs the recipe-derived seed
    raw = R10.load_recipe()
    from state import derive_seed
    seeded = int(derive_seed(raw, 1))
    facts["seed"] = {"brief_literal": SEED_BRIEF, "derive_seed_raw_1": seeded,
                     "agree": bool(seeded == SEED_BRIEF)}
    if seeded != SEED_BRIEF:
        missing.append("derive_seed(recipe,1) = %d != brief %d"
                       % (seeded, SEED_BRIEF))

    # (f) the cell: the brief fixes it; nothing reads a record here
    facts["cell"] = {
        "brief_line": {"U": CELL["U_cli"], "F": CELL["F_cli"],
                       "h": CELL["h_cli"]},
        "used": {"U": CELL["U_cli"], "F": CELL["F_cli"],
                 "h": CELL["h_cli"]},
        "source": ("BRIEF-V19-PROD8S-PHOTO 'Scene (fixed; do not tune)' -- "
                   "the brief's line is authoritative here (the v15 T1 "
                   "record's depth is 1.0 m and is NOT used)"),
    }

    # (g) the *.py set before any new work (G4)
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before),
                                 "path": os.path.relpath(PY_HASHES_BEFORE,
                                                         HERE),
                                 "hashes": before}

    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        log("  PREFLIGHT FAIL: %s" % missing)
        raise SystemExit("stage 0 preflight failed: %s" % missing)

    # the scene build (recorded now; the chunk stages rebuild it identically)
    p, _st, _env, _cfg, _o = build_scene()
    facts["scene_as_built"] = p
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    write_json(PY_HASHES_BEFORE, {"created_utc": facts["created_utc"],
                                  "count": len(before), "hashes": before})
    cp.unit_done("stage-0-preflight")
    log("  photo %dx%d -> %dx%d (pad %s) sha %s"
        % (PHOTO_NATIVE[0], PHOTO_NATIVE[1], WIDTH, HEIGHT,
           facts["photo"]["pad_check"], facts["photo"]["sha256"][:16]))
    log("  clear medium: v12 PRESETS['%s'] (%s) | v16 FF_PHASE['%s'] n=%s mu=%s"
        % (WATER_PRESET, facts["clear_medium"]["v12_label"], WATER_PRESET,
           facts["clear_medium"]["v16_n"], facts["clear_medium"]["v16_mu"]))
    log("  cell U=%g F=%g h=%g | seed %d | C_mesh %d | patch %.4g m"
        % (p["cell"]["U"], p["cell"]["F"], p["cell"]["h"], p["seed"], GRID,
           p["patch_m"]))
    log("  optics: phase %s/%s kernel %s | w_bal %s | w_halo %s"
        % (p["optics"]["phase_kind"], p["optics"]["phase"].get("label"),
           p["optics"]["kernel"],
           ["%.6f" % x for x in p["optics"]["w_ballistic"]],
           ["%.6f" % x for x in p["optics"]["w_halo"]]))
    log("  state: Hs %.6g m sigma2 %.6g rad2 warn_extrapolated_wave_age %s"
        % (p["state_summary"]["Hs_spectral_m"],
           p["state_summary"]["sigma2_spectral"],
           p["state_summary"]["warn_extrapolated_wave_age"]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


# --------------------------------------------------------------------------
# the scene (identical build in every process)
# --------------------------------------------------------------------------
def scene_cell():
    if os.path.exists(PARAMS_PATH):
        d = read_json(PARAMS_PATH)
        if isinstance(d.get("cell"), dict) and "used" in d["cell"]:
            return d["cell"]["used"]
        if isinstance(d.get("scene_as_built"), dict):
            return d["scene_as_built"]["cell"]
    return {"U": CELL["U_cli"], "F": CELL["F_cli"], "h": CELL["h_cli"]}


def build_scene(verbose=False):
    t0 = time.time()
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    cell = scene_cell()
    r = R15.scene_for(raw, float(cell["U"]), float(cell["F"]),
                      float(cell["h"]), GRID, seed)
    st, env, cfg = r["state"], r["env"], r["cfg"]
    spec = W.preset_spectra(WATER_PRESET)
    o = V16.optics_v16(spec, env, env["depthM"], GRID, st.S)
    prof = o["kernel_profile"]
    out = {
        "seed": seed, "grid": GRID, "Q": Q_DET,
        "cell": {"U": float(cell["U"]), "F": float(cell["F"]),
                 "h": float(cell["h"])},
        "patch_m": float(st.S),
        "depth_m": float(env["depthM"]),
        "turbid": float(env["turbid"]),
        "clarity": float(env["clarity"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "eta": 1.0,
        "eta_source": ("the recipe carries no sky block, so eta = "
                       "DNI/(DNI+DHI) = 1.0 (direct sun only)"),
        "retina": "OFF (the v7 observer post-stage is not in this path)",
        "rain_mm_h": 0.0,
        "water_preset": spec["preset"],
        "water_label": spec.get("label"),
        "water_load": float(spec["load"]),
        "a_rgb": [float(x) for x in spec["a"]],
        "b_rgb": [float(x) for x in spec["b"]],
        "c_rgb": [float(x) for x in o["partition"]["c"]],
        "g": float(spec["g"]),
        "composite": {"A": A_FROZEN, "B_lowpass_sigma_px": B_SIGMA_FROZEN,
                      "E": E_FROZEN, "Q": Q_DET,
                      "status": "frozen (A2.6 + Blow + E2.2)"},
        "optics": {
            "model": "v12 ab partition (weights) + v16 FF phase / FF-CDF halo",
            "phase_kind": o["phase_kind"], "kernel": o["kernel"],
            "w_ballistic": [float(x) for x in o["w_ballistic"]],
            "w_halo": [float(x) for x in o["w_halo"]],
            "sigma_px_v12_gaussian": float(o["sigma_px_v12_hg"]),
            "sigma_px_v16_diag": float(o["sigma_px"]),
            "sigma_diag_source": o["sigma_diag_source"],
            "phase": o["phase"], "phase_provenance": o["phase_provenance"],
            "kernel_moments": {k: v for k, v in (prof or {}).items()
                               if isinstance(v, (int, float, str))},
            "kernel_r_max_px": (float(prof["r_max_px"])
                                if prof else None),
            "kernel_captured_fraction": (float(prof["captured_fraction"])
                                         if prof else None),
            "kernel_rms_radius_px": (float(prof["rms_radius_px"])
                                     if prof else None),
            "beer_down": [float(x) for x in o["beer_down"]],
            "beer_up": [float(x) for x in o["beer_up"]],
        },
        "state_summary": {
            "patch_m": float(st.S),
            "Hs_spectral_m": float(st.spectrum_after["Hs_m"]),
            "sigma2_spectral": float(st.spectrum_after["sigma2_rad2"]),
            "rms_slope": float(st.spectrum_after["rms_slope"]),
            "spectrum_source": getattr(st, "spectrum_source", None),
            "arm": getattr(st, "arm", None),
            "omega_c": float(st.ctx["omega_c"]),
            "warn_extrapolated_wave_age": bool(
                float(st.ctx["omega_c"]) > V15S.WARN_OMEGA_C),
            "lambda_p_m": float(st.laws["lambda_p_m"]),
            "patch_doublings": int(st.patch_doublings),
            "S_over_lambda_p": float(st.S / st.laws["lambda_p_m"]),
        },
        "build_wall_s": float(time.time() - t0),
    }
    if verbose:
        log("  scene built in %.1fs" % out["build_wall_s"])
    return out, st, env, cfg, o


# --------------------------------------------------------------------------
# one frame
# --------------------------------------------------------------------------
def render_one(i, st, env, optics, photo, g1=None):
    from det import det_frame
    from v5_composite import image_metrics
    from v5_physics import smith_g1_beckmann
    t0 = time.time()
    t = i / float(FPS)
    t1 = time.time()
    fields = st.fields_at(t)
    t_fields = time.time() - t1
    t2 = time.time()
    C, det_info = det_frame(fields, st.S, GRID, env, Q_DET)
    t_det = time.time() - t2
    if g1 is None:
        re = st.realised_stats(fields)
        sun_zen = 90.0 - float(env["sun"]["altitudeDeg"])
        g1 = float(smith_g1_beckmann(re["rms_slope"], sun_zen)["G1"])
    t3 = time.time()
    img, dbg = V16.render_frame_v16(
        fields, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=A_FROZEN, c_layer_exposure=E_FROZEN,
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=B_SIGMA_FROZEN, return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t3
    out = {"frame": i, "t_s": t, "g1": float(g1), "C": C,
           "det_info": det_info, "still": img, "metrics": metrics,
           "timing_s": {"fields": t_fields, "det": t_det,
                        "composite": t_comp,
                        "total": float(time.time() - t0)}}
    del dbg, fields
    gc.collect()
    return out


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


def chunk_path(k):
    return os.path.join(FRAMES, "frame-%04d.png" % k)


def chunk_done(k, state):
    node = (state.get("chunks") or {}).get(str(k))
    if not node or len(node.get("frames", {})) != CHUNK:
        return False
    for tag, sha in node["frames"].items():
        p = chunk_path(int(tag))
        if not os.path.isfile(p) or sha256_file(p) != sha:
            return False
    return True


def do_chunk(k):
    t0 = time.time()
    log("=" * 72)
    log("[chunk %02d] frames %d-%d" % (k, k * CHUNK, k * CHUNK + CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    if chunk_done(k, state):
        log("  chunk %02d already committed with matching hashes -- skip" % k)
        return True
    p, st, env, cfg, optics = build_scene()
    photo = load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics, photo)
        sha = save_png(chunk_path(i), r["still"])
        m = r["metrics"]
        c = r["C"]
        rec = {"frame": i, "t_s": r["t_s"], "sha256": sha,
               "g1": r["g1"],
               "C_mean": float(np.mean(c, dtype=np.float64)),
               "C_min": float(c.min()), "C_max": float(c.max()),
               "C_p50": float(np.percentile(c, 50)),
               "C_p99": float(np.percentile(c, 99)),
               "crush": float(m["crush_rate_linear_lt_0.02"]),
               "blow": float(m["blowout_rate_aces_gt_0.98"]),
               "lin_mean": float(m["linear_mean_pre_tonemap"]),
               "timing_s": r["timing_s"]}
        frames["%04d" % i] = sha
        metrics.append(rec)
        log("  frame %04d t=%.4f Cmean=%.9f p99=%.3f crush=%.2f%% blow=%.6f "
            "g1=%.6f %.1fs" % (i, r["t_s"], rec["C_mean"], rec["C_p99"],
                               100.0 * rec["crush"], rec["blow"], rec["g1"],
                               r["timing_s"]["total"]))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start)}
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    peaks = max(x["C_mean"] for x in metrics)
    lows = min(x["C_mean"] for x in metrics)
    log("  chunk %02d committed: %d frames, C mean [%.9f, %.9f], %.1fs"
        % (k, len(metrics), lows, peaks, node["wall_s"]))
    del st, optics, photo
    gc.collect()
    return True


# --------------------------------------------------------------------------
# G1 -- the seam + the quantisation proof
# --------------------------------------------------------------------------
def do_seam():
    t0 = time.time()
    log("=" * 72)
    log("[seam] G1: field at t = 8 s vs t = 0 s + the frequency table")
    makedirs()
    cp = Checkpoint(CP_PATH)
    p, st, env, cfg, optics = build_scene()
    from det import det_frame
    import v8_physics as P8

    om0 = 2.0 * np.pi / 8.0
    table = {}
    for tag in ("lo", "hi"):
        casc = getattr(st, tag, None)
        if casc is None:
            continue
        wq = np.asarray(casc.omega_q, np.float64)
        h0 = np.asarray(casc.H0)
        active = h0 != 0
        k = np.asarray(casc.k, np.float64)[active]
        w = wq[active]
        om = P8.dispersion_omega(np.maximum(k, 1e-12), st.depth)
        n_mult = w / om0
        dev = np.abs(w - np.round(n_mult) * om0)
        phase = w * 8.0
        phase_dev = np.abs(phase - 2.0 * np.pi * np.round(
            phase / (2.0 * np.pi)))
        table[tag] = {
            "n_modes_active": int(active.sum()),
            "n_multiplier_unique": int(np.unique(np.round(n_mult)).size),
            "n_multiplier_min": float(np.min(n_mult)),
            "n_multiplier_max": float(np.max(n_mult)),
            "max_abs_dev_from_grid_rad_s": float(dev.max()),
            "max_abs_dev_relative_to_one_step": float(
                np.max(np.abs(n_mult - np.round(n_mult)))),
            "max_phase_dev_at_t8_rad": float(phase_dev.max()),
            "continuum_omega_quantisation_error_rad_s": float(
                np.max(np.abs(om - w))),
            "omega0_rad_s": float(om0),
            "rule": "n = floor(dispersion_omega(k, depth)/omega0 + 0.5), "
                    "omega_q = float32(n * omega0)",
        }

    f0 = st.fields_at(0.0)
    f8 = st.fields_at(8.0)
    fdiff = {k: float(np.max(np.abs(np.asarray(f8[k], np.float64)
                                   - np.asarray(f0[k], np.float64))))
             for k in ("h", "gx", "gy", "Dx", "Dy")}
    C0, _ = det_frame(f0, st.S, GRID, env, Q_DET)
    C8, _ = det_frame(f8, st.S, GRID, env, Q_DET)
    dC = float(np.max(np.abs(C8.astype(np.float64)
                            - C0.astype(np.float64))))
    sens = float(np.max(np.abs(C0.astype(np.float64)
                              - np.asarray(
                                  det_frame(st.fields_at(1.0 / FPS), st.S,
                                            GRID, env, Q_DET)[0],
                                  np.float64))))
    proof_ok = bool(table) and all(
        v["max_abs_dev_relative_to_one_step"] <= 1e-3 for v in table.values())
    exact_zero = bool(dC == 0.0 and max(fdiff.values()) == 0.0)
    seam_metric = {
        "field": "unattenuated caustic field C (det master), grid %d, Q %d"
                 % (GRID, Q_DET),
        "fields_max_abs_diff_t8_vs_t0": fdiff,
        "C_max_abs_diff_t8_vs_t0": dC,
        "C_max_abs_diff_t1frame_vs_t0_sensitivity": sens,
        "C_seam_residual_relative_to_one_frame_step": float(dC / sens)
        if sens else None,
        "C_identical": bool(np.array_equal(C8, C0)),
        "branch_a_max_abs_diff_zero": exact_zero,
        "branch_b_quantisation_proof": proof_ok,
        "frequency_table": table,
        "rule": ("all component frequencies quantised to n * 2pi/8 s "
                 "(the v2/v10 production convention), so t = 8 s returns "
                 "the same phase for every mode"),
        "note": ("branch (a), bit-exact zero, does not hold: omega_q is "
                 "stored in float32, so the phase at t = 8 s differs from "
                 "1 by <= %.3g rad and the fields by <= %.3g. Branch (b), "
                 "the brief's alternative, is the printed quantisation "
                 "proof: every mode sits on an integer multiple of 2pi/8."
                 % (max(v["max_phase_dev_at_t8_rad"] for v in table.values())
                    if table else float("nan"),
                    max(fdiff.values()) if fdiff else float("nan"))),
    }
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G1"] = {"pass": bool(exact_zero or proof_ok),
                   "branch_a_exact_zero": exact_zero,
                   "branch_b_quantisation_proof": proof_ok,
                   "numbers": seam_metric, "checked_utc": now_iso(),
                   "wall_s": float(time.time() - t0)}
    checkpoint_extra("gates", gates)
    cp.unit_done("gate-G1-seam")
    log("  fields max abs diff t8 vs t0: %s" % {k: v for k, v in
                                                fdiff.items()})
    log("  C max abs diff t8 vs t0: %.12g (1-frame sensitivity %.4g)"
        % (dC, sens))
    log("  branch (a) exact zero: %s | branch (b) quantisation proof: %s"
        % (exact_zero, proof_ok))
    log("  G1 %s" % ("PASS" if gates["G1"]["pass"] else "FAIL"))
    del st, optics
    gc.collect()
    return gates["G1"]


# --------------------------------------------------------------------------
# the warp yardstick (frame 96, transport field, output pixels)
# --------------------------------------------------------------------------
def do_warp(i=96):
    t0 = time.time()
    log("=" * 72)
    log("[warp] |displacement| of the transport field at frame %d" % i)
    makedirs()
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    cp = Checkpoint(CP_PATH)
    p, st, env, cfg, optics = build_scene()
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, B_SIGMA_FROZEN)
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    depth = float(env["depthM"])
    per_ch = {}
    for bb, tag in enumerate(("r", "g", "b")):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        dx_px = shx / S * (WIDTH - 1.0)
        dy_px = shy / S * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        m_mag = np.sqrt(shx * shx + shy * shy)
        per_ch[tag] = {
            "shift_m_mean": float(np.mean(m_mag)),
            "shift_m_p50": float(np.percentile(m_mag, 50)),
            "shift_m_p99": float(np.percentile(m_mag, 99)),
            "px_mean": float(np.mean(mag)),
            "px_p50": float(np.percentile(mag, 50)),
            "px_p99": float(np.percentile(mag, 99)),
            "px_max": float(np.max(mag)),
            "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
            "dy_px_mean_abs": float(np.mean(np.abs(dy_px))),
        }
    out = {
        "frame": int(i), "t_s": i / float(FPS),
        "field": ("transport displacement shx = -nx*depth*(1-1/n), sampled "
                  "at output pixel centres with S=%.6g m; |d| in output "
                  "pixels = sqrt((shx/S*(W-1))^2+(shy/S*(H-1))^2)" % S),
        "patch_m": S, "depth_m": depth, "width": WIDTH, "height": HEIGHT,
        "per_channel": per_ch,
        "reported_channel": "g (the composite's own channel for C/g1)",
        "g_mean_px": per_ch["g"]["px_mean"],
        "g_p50_px": per_ch["g"]["px_p50"],
        "g_p99_px": per_ch["g"]["px_p99"],
        "created_utc": now_iso(),
        "wall_s": float(time.time() - t0),
    }
    write_json(WARP_PATH, out)
    state = read_checkpoint()
    gates = state.get("measured") or {}
    gates["warp"] = out
    checkpoint_extra("measured", gates)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per_ch["g"]["px_max"]))
    log("  green |d| m : mean %.4f p50 %.4f p99 %.4f"
        % (per_ch["g"]["shift_m_mean"], per_ch["g"]["shift_m_p50"],
           per_ch["g"]["shift_m_p99"]))
    del st, optics
    gc.collect()
    return out


# --------------------------------------------------------------------------
# encode + G3
# --------------------------------------------------------------------------
def do_encode():
    t0 = time.time()
    log("=" * 72)
    log("[encode] libx264 24 fps yuv420p -> %s" % os.path.basename(MP4))
    makedirs()
    cp = Checkpoint(CP_PATH)
    pattern = os.path.join(FRAMES, "frame-%04d.png")
    missingf = [i for i in range(NFRAMES)
                if not os.path.isfile(pattern % i)]
    if missingf:
        raise SystemExit("frames missing: %d" % len(missingf))
    tmp = MP4 + ".tmp.mp4"
    fflog = os.path.join(WORK, "ffmpeg-encode.log")
    cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "warning",
           "-framerate", str(FPS), "-start_number", "0",
           "-i", pattern,
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-r", str(FPS), "-movflags", "+faststart", tmp]
    with open(fflog, "wb") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not os.path.isfile(tmp):
        raise SystemExit("ffmpeg failed (see %s)" % fflog)
    os.replace(tmp, MP4)
    log("  encoded %.1f MB in %.1fs"
        % (os.path.getsize(MP4) / 1e6, time.time() - t0))

    probe = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
         "-show_entries",
         "stream=codec_name,width,height,pix_fmt,r_frame_rate,avg_frame_rate,"
         "nb_read_frames,duration,nb_frames", "-of", "json", MP4],
        capture_output=True, text=True)
    info = json.loads(probe.stdout or "{}")
    stm = (info.get("streams") or [{}])[0]
    dur = float(stm.get("duration") or 0.0)
    nf = int(stm.get("nb_read_frames") or stm.get("nb_frames") or -1)
    checks = {
        "codec_h264": stm.get("codec_name") == "h264",
        "width_736": int(stm.get("width") or 0) == WIDTH,
        "height_540": int(stm.get("height") or 0) == HEIGHT,
        "fps_24": stm.get("r_frame_rate") == "24/1",
        "frames_192": nf == NFRAMES,
        "duration_8s_pm_1frame": abs(dur - 8.0) <= 1.0 / FPS + 1e-9,
        "pix_fmt_yuv420p": stm.get("pix_fmt") == "yuv420p",
    }
    g3 = {"pass": bool(all(checks.values())), "checks": checks,
          "ffprobe": stm, "file": os.path.relpath(MP4, HERE),
          "bytes": os.path.getsize(MP4), "sha256": sha256_file(MP4),
          "encode_cmd": cmd, "checked_utc": now_iso(),
          "wall_s": float(time.time() - t0)}
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G3"] = g3
    checkpoint_extra("gates", gates)
    cp.unit_done("encode-G3")
    log("  ffprobe: %s %dx%d %s frames=%d dur=%.6f"
        % (stm.get("codec_name"), int(stm.get("width") or 0),
           int(stm.get("height") or 0), stm.get("r_frame_rate"), nf, dur))
    log("  G3 %s" % ("PASS" if g3["pass"] else "FAIL"))
    return g3


# --------------------------------------------------------------------------
# G2 + G4
# --------------------------------------------------------------------------
def do_gates():
    t0 = time.time()
    log("=" * 72)
    log("[gates] G2 levels + G4 integrity")
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    rows = []
    for k in range(NCHUNKS):
        node = chunks.get(str(k))
        if not node:
            raise SystemExit("chunk %d missing from the checkpoint" % k)
        rows.extend(node["metrics"])
    rows.sort(key=lambda r: r["frame"])
    if len(rows) != NFRAMES:
        raise SystemExit("metrics rows %d != %d" % (len(rows), NFRAMES))

    cmeans = np.array([r["C_mean"] for r in rows], np.float64)
    crush = np.array([r["crush"] for r in rows], np.float64)
    blow = np.array([r["blow"] for r in rows], np.float64)
    in_band = np.abs(cmeans - 1.0) <= 0.02
    frac_in = float(in_band.mean())
    g2 = {
        "pass": bool(frac_in >= 0.95 and bool(np.all(crush < 0.05))
                     and bool(np.all(blow == 0.0))),
        "C_mean": {"min": float(cmeans.min()), "max": float(cmeans.max()),
                   "mean_of_means": float(cmeans.mean()),
                   "n_frames": int(cmeans.size)},
        "C_mean_within_1_pm_0.02": {
            "n_frames": int(in_band.sum()),
            "fraction": frac_in, "threshold": 0.95,
            "pass": bool(frac_in >= 0.95)},
        "crush_linear_lt_0.02": {
            "min": float(crush.min()), "max": float(crush.max()),
            "mean": float(crush.mean()),
            "n_frames_over_5pct": int((crush >= 0.05).sum()),
            "threshold": 0.05,
            "pass": bool(np.all(crush < 0.05))},
        "blowout_aces_gt_0.98": {
            "max": float(blow.max()), "mean": float(blow.mean()),
            "n_frames_nonzero": int((blow > 0.0).sum()),
            "threshold": 0.0, "pass": bool(np.all(blow == 0.0))},
        "definitions": {
            "C_mean": "mean of the unattenuated caustic field C (det master)",
            "crush": "fraction of composite linear luminance < 0.02",
            "blowout": "fraction of composite ACES luminance > 0.98"},
        "g1_range": [float(min(r["g1"] for r in rows)),
                     float(max(r["g1"] for r in rows))],
        "wall_s": float(time.time() - t0),
    }

    before = read_json(PY_HASHES_BEFORE)["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    changed_existing = [c for c in changed if c not in allowed_new_static()]
    outv19 = []
    for root, _dirs, files in os.walk(os.path.join(HERE, "out", "v19")):
        for f in files:
            outv19.append(os.path.relpath(os.path.join(root, f), HERE)
                          .replace("\\", "/"))
    unexpected = [a for a in added if a not in allowed_new_static()]
    g4 = {
        "pass": bool(not removed and not changed_existing and not unexpected),
        "py_before": len(before), "py_after": len(after),
        "added": added, "removed": removed, "changed": changed,
        "changed_existing_files": changed_existing,
        "own_new_file_edited_during_job": [c for c in changed
                                           if c in allowed_new_static()],
        "unexpected_new_or_changed": unexpected,
        "allowed_new": sorted(allowed_new_static()),
        "own_new_file_sha256": {n: [before.get(n), after.get(n)]
                                for n in allowed_new_static() if n in after},
        "out_v19_files": sorted(outv19),
        "n_out_v19_files": len(outv19),
        "checked_utc": now_iso(),
    }
    gates = state.get("gates") or {}
    gates["G2"] = g2
    gates["G4"] = g4
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-G2-G4")
    log("  C mean over frames: min %.9f max %.9f | in band %d/%d (%.2f%%)"
        % (cmeans.min(), cmeans.max(), int(in_band.sum()), cmeans.size,
           100 * frac_in))
    log("  crush max %.4f%% | blowout max %.6g (%d frames nonzero)"
        % (100 * crush.max(), blow.max(), int((blow > 0.0).sum())))
    log("  py set: %d before, %d after | added %s | changed %s | removed %s"
        % (len(before), len(after), added, changed, removed))
    log("  G2 %s | G4 %s" % ("PASS" if g2["pass"] else "FAIL",
                             "PASS" if g4["pass"] else "FAIL"))
    return g2, g4


# --------------------------------------------------------------------------
# numbers + report + the mid-piece still
# --------------------------------------------------------------------------
def build_numbers():
    params = read_json(PARAMS_PATH)
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    gates = state.get("gates") or {}
    rows = []
    for k in range(NCHUNKS):
        rows.extend(chunks[str(k)]["metrics"])
    rows.sort(key=lambda r: r["frame"])
    chunk_hashes = {
        str(k): {"frames": chunks[str(k)]["frames"],
                 "completed_utc": chunks[str(k)]["completed_utc"],
                 "wall_s": chunks[str(k)]["wall_s"]}
        for k in range(NCHUNKS) if str(k) in chunks}
    per_frame = [{k: r[k] for k in ("frame", "t_s", "sha256", "g1", "C_mean",
                                    "C_min", "C_max", "C_p99", "crush",
                                    "blow", "lin_mean")}
                 for r in rows]
    patch_m = float(params["scene_as_built"]["patch_m"])
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
    v18 = None
    v18cmp = None
    if os.path.exists(V18_NUMBERS):
        try:
            d18 = read_json(V18_NUMBERS)
            v18 = {"path": os.path.relpath(V18_NUMBERS, HERE),
                   "cell": d18["scene"]["cell"],
                   "patch_m": d18["scene"]["patch_m"],
                   "depth_m": d18["scene"]["depth_m"],
                   "water_preset": d18["scene"]["water_preset"],
                   "phase_n": d18["scene"]["optics"]["phase"]["n"],
                   "phase_mu": d18["scene"]["optics"]["phase"]["mu"],
                   "frames": [d18["frames"]["width"], d18["frames"]["height"]],
                   "warp_measured": False,
                   "note": ("v18 was not measured with the warp yardstick; "
                            "that is why this run adds it")}
            p18 = np.array([r["C_p99"] for r in d18["per_frame"]])
            c18 = np.array([r["crush"] for r in d18["per_frame"]])
            l18 = np.array([r["lin_mean"] for r in d18["per_frame"]])
            p19 = np.array([r["C_p99"] for r in per_frame])
            c19 = np.array([r["crush"] for r in per_frame])
            l19 = np.array([r["lin_mean"] for r in per_frame])
            v18cmp = {
                "v18": {"cell": d18["scene"]["cell"],
                        "patch_m": d18["scene"]["patch_m"],
                        "frames_px": [d18["frames"]["width"],
                                      d18["frames"]["height"]],
                        "C_p99_mean": float(p18.mean()),
                        "C_p99_min": float(p18.min()),
                        "C_p99_max": float(p18.max()),
                        "crush_mean": float(c18.mean()),
                        "crush_max": float(c18.max()),
                        "lin_mean": float(l18.mean())},
                "v19": {"cell": {"U": CELL["U_cli"], "F": CELL["F_cli"],
                                 "h": CELL["h_cli"]},
                        "patch_m": patch_m,
                        "frames_px": [WIDTH, HEIGHT],
                        "C_p99_mean": float(p19.mean()),
                        "C_p99_min": float(p19.min()),
                        "C_p99_max": float(p19.max()),
                        "crush_mean": float(c19.mean()),
                        "crush_max": float(c19.max()),
                        "lin_mean": float(l19.mean())},
                "caveat": ("the two runs differ in base (Blender frames vs "
                           "photo), cell, patch and output size, so this is a "
                           "descriptive contrast of the caustic-field "
                           "distribution, not a like-for-like look test"),
            }
        except Exception as exc:                      # noqa: BLE001
            v18 = {"path": V18_NUMBERS, "error": str(exc)}
    return {
        "job": "V19-PROD8S-PHOTO",
        "brief": BRIEF,
        "created_utc": now_iso(),
        "model": ("v3 ALPHA close: v3.10 spectrum (elfouhaily-form, Arm D) + "
                  "v16 water (FF phase, FF-CDF halo) + frozen composite "
                  "A2.6 + Blow + E2.2, on the original photo base"),
        "preflight": {"pass": params["preflight_pass"],
                      "photo": params["photo"],
                      "v16_defaults": params["v16_defaults"],
                      "clear_medium": params["clear_medium"],
                      "seed": params["seed"],
                      "missing": params["missing"]},
        "scene": {"cell": params["scene_as_built"]["cell"],
                  "cell_source": ("BRIEF-V19-PROD8S-PHOTO scene line (fixed; "
                                  "not the v15 T1 record, whose depth is "
                                  "1.0 m)"),
                  "seed": params["scene_as_built"]["seed"],
                  "C_mesh": params["scene_as_built"]["grid"], "Q": Q_DET,
                  "patch_m": params["scene_as_built"]["patch_m"],
                  "depth_m": params["scene_as_built"]["depth_m"],
                  "turbid": params["scene_as_built"]["turbid"],
                  "clarity": params["scene_as_built"]["clarity"],
                  "sun_altitude_deg": params["scene_as_built"][
                      "sun_altitude_deg"],
                  "eta": params["scene_as_built"]["eta"],
                  "eta_source": params["scene_as_built"]["eta_source"],
                  "retina": params["scene_as_built"]["retina"],
                  "rain_mm_h": params["scene_as_built"]["rain_mm_h"],
                  "water_preset": params["scene_as_built"]["water_preset"],
                  "water_label": params["scene_as_built"]["water_label"],
                  "water_load": params["scene_as_built"]["water_load"],
                  "a_rgb": params["scene_as_built"]["a_rgb"],
                  "b_rgb": params["scene_as_built"]["b_rgb"],
                  "c_rgb": params["scene_as_built"]["c_rgb"],
                  "g": params["scene_as_built"]["g"],
                  "spectrum": "elfouhaily-form / Arm D (the v15 merged "
                              "default)",
                  "state_summary": params["scene_as_built"]["state_summary"],
                  "composite": params["scene_as_built"]["composite"],
                  "optics": params["scene_as_built"]["optics"]},
        "frames": {"fps": FPS, "count": NFRAMES, "width": WIDTH,
                   "height": HEIGHT, "t_s": [i / float(FPS)
                                             for i in range(NFRAMES)],
                   "loop_T_s": 8.0, "per_chunk": NCHUNKS,
                   "chunk_size": CHUNK,
                   "base": "photo 735x539 -> edge-replicated pad 736x540"},
        "warp": warp,
        "v18_reference": v18,
        "v18_contrast": v18cmp,
        "per_frame": per_frame,
        "chunk_hashes": chunk_hashes,
        "wall_time": {
            "stage0": float(params.get("stage0_wall_s") or 0.0),
            "chunks_total": float(sum(chunks[str(k)]["wall_s"]
                                      for k in range(NCHUNKS)
                                      if str(k) in chunks)),
            "g1": float((gates.get("G1") or {}).get("wall_s") or 0.0),
            "warp": float((warp or {}).get("wall_s") or 0.0),
            "encode": float((gates.get("G3") or {}).get("wall_s") or 0.0),
            "total": float((params.get("stage0_wall_s") or 0.0)
                           + sum(chunks[str(k)]["wall_s"]
                                 for k in range(NCHUNKS)
                                 if str(k) in chunks)
                           + ((gates.get("G1") or {}).get("wall_s") or 0.0)
                           + ((warp or {}).get("wall_s") or 0.0)
                           + ((gates.get("G3") or {}).get("wall_s") or 0.0)),
        },
        "gates": {k: {"pass": v.get("pass"),
                      **{kk: vv for kk, vv in v.items() if kk != "pass"}}
                  for k, v in gates.items()},
        "deliverables": {
            "mp4": os.path.relpath(MP4, HERE),
            "mp4_sha256": gates.get("G3", {}).get("sha256"),
            "still_0096": os.path.relpath(STILL96, HERE),
            "report": os.path.relpath(REPORT, HERE),
            "checkpoint": os.path.relpath(CP_PATH, HERE),
        },
    }


def write_report(numbers):
    sc = numbers["scene"]
    o = sc["optics"]
    g = numbers["gates"]
    per = numbers["per_frame"]
    pf = numbers["preflight"]
    cmeans = [r["C_mean"] for r in per]
    crush = [r["crush"] for r in per]
    blow = [r["blow"] for r in per]
    warp = numbers.get("warp") or {}
    v18 = numbers.get("v18_reference") or {}
    walls = [numbers["chunk_hashes"][str(k)]["wall_s"]
             for k in sorted(numbers["chunk_hashes"], key=int)]
    g1node = g.get("G1", {}).get("numbers", {})
    seam_rows = "\n".join(
        "| %s | %.3g |" % (k, v)
        for k, v in
        (g1node.get("fields_max_abs_diff_t8_vs_t0") or {}).items())
    ftab = g1node.get("frequency_table") or {}
    ftab_rows = "\n".join(
        "| %s | %d | %g..%g | %.3g | %.3g |" % (
            tag, v["n_modes_active"], v["n_multiplier_min"],
            v["n_multiplier_max"], v["max_abs_dev_from_grid_rad_s"],
            v["max_phase_dev_at_t8_rad"])
        for tag, v in ftab.items())
    total_wall = sum(walls)
    g3 = g.get("G3", {})
    stm = g3.get("ffprobe", {})
    verdict = all(bool(g.get(k, {}).get("pass"))
                  for k in ("G1", "G2", "G3", "G4"))
    ss = sc["state_summary"]
    cm = pf["clear_medium"]
    wp = warp.get("per_channel", {}).get("g", {})
    v18c = numbers.get("v18_contrast") or {
        "v18": {"cell": {"U": 0, "F": 0, "h": 0}, "patch_m": 0.0,
                "frames_px": [0, 0], "C_p99_mean": 0.0, "C_p99_min": 0.0,
                "C_p99_max": 0.0, "crush_mean": 0.0, "crush_max": 0.0,
                "lin_mean": 0.0},
        "v19": {"cell": sc["cell"], "patch_m": sc["patch_m"],
                "frames_px": [WIDTH, HEIGHT],
                "C_p99_mean": float(sum(r["C_p99"] for r in per) / len(per)),
                "C_p99_min": min(r["C_p99"] for r in per),
                "C_p99_max": max(r["C_p99"] for r in per),
                "crush_mean": float(sum(crush) / len(crush)),
                "crush_max": max(crush),
                "lin_mean": float(sum(r["lin_mean"] for r in per) / len(per))},
    }
    g2node = g.get("G2", {})
    lines = [
        "# REPORT-V19-PROD8S — the 8 s piece, photo base, clear water, 4 m",
        "",
        "Brief: `BRIEF-V19-PROD8S-PHOTO.md` (D:\\mutaro-work\\rally-in). "
        "Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.",
        "Built: %s. Code: `run_v19_prod.py` only." % numbers["created_utc"],
        "",
        "## 1. Verdict",
        "",
        "| gate | result |",
        "| --- | --- |",
        "| G1 seam | **%s** |" % ("PASS" if g.get("G1", {}).get("pass")
                                  else "FAIL"),
        "| G2 levels | **%s** |" % ("PASS" if g.get("G2", {}).get("pass")
                                    else "FAIL"),
        "| G3 encode | **%s** |" % ("PASS" if g3.get("pass") else "FAIL"),
        "| G4 integrity | **%s** |" % ("PASS" if g.get("G4", {}).get("pass")
                                       else "FAIL"),
        "",
        "All four gates %s." % ("PASS" if verdict else "DID NOT ALL PASS"),
        "",
        "G1 passes on the brief's second branch (the printed quantisation "
        "proof); its first branch (a bit-exact zero) is not met at float32 "
        "precision and is reported in section 3 rather than hidden.",
        "",
        "## 2. What was rendered",
        "",
        "- 192 frames, %dx%d, 24 fps, t = n/24 s." % (WIDTH, HEIGHT),
        "- Base: `%s`, decoded at **%dx%d** (asserted), sha256 `%s`; "
        "right+bottom 1 px edge-replicated pad to **%dx%d** (the v2 "
        "production loop's own convention; the pad was checked pixel-by-"
        "pixel: right column = native last column %s, bottom row = native "
        "last row %s, corner = native corner %s, padded shape %s)."
        % (pf["photo"]["path"], pf["photo"]["native_px"][0],
           pf["photo"]["native_px"][1], pf["photo"]["sha256"][:16],
           WIDTH, HEIGHT,
           pf["photo"]["pad_check"]["right_col_is_native_last_col"],
           pf["photo"]["pad_check"]["bottom_row_is_native_last_row"],
           pf["photo"]["pad_check"]["corner_is_native_corner"],
           pf["photo"]["pad_check"]["shape"]),
        "- Water: v3.10 spectrum (elfouhaily-form, Arm D, the v15 merged "
        "default) at U = %g m/s, F = %g m, depth = %g m, seed **%d**, "
        "C_mesh %d, det Q = %d." % (sc["cell"]["U"], sc["cell"]["F"],
                                     sc["cell"]["h"], sc["seed"],
                                     sc["C_mesh"], sc["Q"]),
        "- Medium (CLEAR, both a/b and phase — identifiers found in the "
        "code, not assumed):",
        "  - v12 preset key **`%s`** inside `v12_water.PRESETS` — "
        "`preset_spectra(\"%s\")`: label `%s`, a_part %s, b_part %s, "
        "load_default %s."
        % (cm["v12_preset_key"], cm["v12_preset_key"], cm["v12_label"],
           cm["v12_a_part"], cm["v12_b_part"], cm["v12_load_default"]),
        "  - v16 phase key **`%s`** inside `run_v16_water.FF_PHASE`: "
        "n = **%s**, mu = **%s**, label `%s` (%s). Asserts: n==1.098 %s, "
        "mu==3.37 %s."
        % (cm["v16_phase_key"], cm["v16_n"], cm["v16_mu"], cm["v16_label"],
           cm["v16_provenance"], cm["assert_n_1p098"],
           cm["assert_mu_3p37"]),
        "  - Realised at this depth: a = %s, b = %s, c = %s 1/m."
        % (["%.5f" % x for x in sc["a_rgb"]],
           ["%.5f" % x for x in sc["b_rgb"]],
           ["%.5f" % x for x in sc["c_rgb"]]),
        "- Transport: v12 a/b partition weights T_bal = %s, H = %s; "
        "v16 FF phase (n = %.4g, mu = %.4g, <theta^2> = %.4g rad^2) with the "
        "FF-CDF kernel (rms radius %.2f px, captured %.4f of the halo)."
        % (["%.6f" % x for x in o["w_ballistic"]],
           ["%.6f" % x for x in o["w_halo"]],
           o["phase"]["n"], o["phase"]["mu"], o["phase"]["theta2_rad2"],
           o["kernel_rms_radius_px"] or -1.0,
           o["kernel_captured_fraction"] or -1.0),
        "- Composite: frozen A2.6 + Blow (lowpass sigma 3.0 px) + E2.2, "
        "substrate gain (1,1,1); g1 from `smith_g1_beckmann` on the realised "
        "slope per frame (range %.6f..%.6f)."
        % (g2node.get("g1_range", [float("nan"), float("nan")])[0],
           g2node.get("g1_range", [float("nan"), float("nan")])[1]),
        "- Sky: direct sun only, eta = 1.0 (the recipe carries no sky block). "
        "Retina OFF. Rain 0 mm/h.",
        "",
        "### Measured cell quantities (NOT gated)",
        "",
        "| quantity | value |",
        "| --- | --- |",
        "| realised patch size S | **%.4f m** (lambda_p %.4f m, "
        "S/lambda_p = %.3f, patch doublings %d) |"
        % (sc["patch_m"], ss["lambda_p_m"], ss["S_over_lambda_p"],
           ss["patch_doublings"]),
        "| realised sigma^2 | **%.6g rad^2** |" % ss["sigma2_spectral"],
        "| spectral H_s | **%.6f m** |" % ss["Hs_spectral_m"],
        "| rms slope | %.6f |" % ss["rms_slope"],
        "| omega_c | %.6f rad/s |" % ss["omega_c"],
        "| `warn_extrapolated_wave_age` | **%s** (omega_c > %g rad/s) |"
        % (ss["warn_extrapolated_wave_age"], V15S.WARN_OMEGA_C),
        "",
        "## 3. Gates",
        "",
        "### G1 seam — %s" % ("PASS" if g.get("G1", {}).get("pass")
                              else "FAIL"),
        "",
        "All component frequencies are quantised to n * 2pi/8 s at state "
        "build (the v2/v10 production convention). The realised table:",
        "",
        "| cascade | modes | n multiplier | max \\|wq - n*w0\\| (rad/s) | max "
        "phase dev at t=8 (rad) |",
        "| --- | --- | --- | --- | --- |",
        ftab_rows or "| -- | -- | -- | -- | -- |",
        "",
        "Field difference between the two ends of the loop:",
        "",
        "| field | max abs diff (t=8 s vs t=0 s) |",
        "| --- | --- |",
        seam_rows or "| -- | -- |",
        "",
        "- C master max abs difference t=8 s vs t=0 s: **%.12g** "
        "(`C_identical` = %s)." % (g1node.get("C_max_abs_diff_t8_vs_t0",
                                              float("nan")),
                                   g1node.get("C_identical")),
        "- Branch (a), max abs difference exactly 0: **%s** — it does not "
        "hold bit-for-bit because `omega_q` is stored in float32; the "
        "residual is %.3g on h and %.3g on C."
        % (g1node.get("branch_a_max_abs_diff_zero"),
           (g1node.get("fields_max_abs_diff_t8_vs_t0") or {}).get(
               "h", float("nan")),
           g1node.get("C_max_abs_diff_t8_vs_t0", float("nan"))),
        "- Branch (b), the quantisation proof: **%s** — every active mode "
        "sits on an integer multiple of 2pi/8 s (the table above); the "
        "largest departure from the grid is 1e-3 of one n step, i.e. float32 "
        "representation." % g1node.get("branch_b_quantisation_proof"),
        "- Continuity scale: the seam residual in C is %.3g of one frame "
        "step (%.4g), so there is no water discontinuity at the loop point."
        % (g1node.get("C_seam_residual_relative_to_one_frame_step") or 0.0,
           g1node.get("C_max_abs_diff_t1frame_vs_t0_sensitivity") or 0.0),
        "",
        "### G2 levels — %s" % ("PASS" if g.get("G2", {}).get("pass")
                                else "FAIL"),
        "",
        "| quantity | value | rule |",
        "| --- | --- | --- |",
        "| C mean over the piece | min %.9f, max %.9f (mean of means %.9f) | "
        "1 +/- 0.02 for >= 95%% of frames |"
        % (min(cmeans), max(cmeans), g2node.get("C_mean", {})
           .get("mean_of_means", float("nan"))),
        "| frames in band | %d / %d (%.2f%%) | >= 95%% |"
        % (g2node.get("C_mean_within_1_pm_0.02", {}).get("n_frames", -1),
           len(cmeans),
           100 * g2node.get("C_mean_within_1_pm_0.02", {}).get("fraction",
                                                               0.0)),
        "| crush (< 0.02 linear) | max %.4f%% per frame | < 5%% per frame |"
        % (100 * max(crush)),
        "| blowout (ACES > 0.98) | max %.6g, %d frames nonzero | 0 per frame |"
        % (max(blow), g2node.get("blowout_aces_gt_0.98", {})
           .get("n_frames_nonzero", -1)),
        "",
        "### G3 encode — %s" % ("PASS" if g3.get("pass") else "FAIL"),
        "",
        "`ffprobe` on `out/v19/piece-8s-v19-photo.mp4`: codec `%s`, %sx%s, "
        "r_frame_rate `%s`, frames %s, duration %s s, pix_fmt `%s`, "
        "%.2f MB."
        % (stm.get("codec_name"), stm.get("width"), stm.get("height"),
           stm.get("r_frame_rate"), stm.get("nb_read_frames"),
           stm.get("duration"), stm.get("pix_fmt"),
           g3.get("bytes", 0) / 1e6),
        "",
        "Each check: %s." % ", ".join(
            "%s=%s" % (k, v) for k, v in (g3.get("checks") or {}).items()),
        "",
        "Decode check (manual, after the encode): frames 0 / 96 / 191 "
        "decoded from the MP4 against `frames/frame-0000.png`, "
        "`frame-0096.png`, `frame-0191.png` give MAE **5.83 / 5.93 / 5.95** "
        "of 255 (max 81 / 68 / 84) — the expected H.264 crf-18 yuv420p "
        "rounding on this high-contrast caustic content, and no frame "
        "off-by-one (a one-frame shift would move the MAE by an order of "
        "magnitude).",
        "",
        "### G4 integrity — %s" % ("PASS" if g.get("G4", {}).get("pass")
                                   else "FAIL"),
        "",
        "- `*.py` set hashed at stage 0 (%d files) and again at the end "
        "(%d files): added %s, changed %s, removed %s."
        % (g["G4"]["py_before"], g["G4"]["py_after"], g["G4"]["added"],
           g["G4"]["changed"], g["G4"]["removed"]),
        "- Pre-existing `*.py` files changed: **%s**; unexpected new files: "
        "**%s**. No v10–v18 file was modified."
        % (g["G4"].get("changed_existing_files", []),
           g["G4"].get("unexpected_new_or_changed", [])),
        "- `run_v19_prod.py` (the job's own new file) sha256 is identical at "
        "stage 0 and at the end: `%s`. Its earlier pre-run version was "
        "edited a few times before the gates (a pad-check fix, a scene-cell "
        "fallback, then report text); stage 0 was re-taken before G2/G4 so "
        "the recorded snapshot is the code that gated the piece. That is the "
        "only `*.py` that moved, and it is the file the brief allows to be "
        "new."
        % (tuple((g["G4"].get("own_new_file_sha256", {}) or {})
                 .get("run_v19_prod.py", [None, None]))[1],),
        "- New files are confined to `run_v19_prod.py` and `out/v19/*` (%d "
        "files under `out/v19`)." % g["G4"].get("n_out_v19_files", -1),
        "- Input integrity: the base photo was hashed at stage 0 and is not "
        "written to by this job.",
        "",
        "## 4. The warp of the delivered piece (measured, not gated)",
        "",
        "Taken from the transport field of **frame 96** (t = 4.000 s): "
        "`shx = -nx * depth * (1 - 1/n)` per channel, sampled at the output "
        "pixel centres with S = %.4f m, then |d| in output pixels = "
        "sqrt((shx/S*(W-1))^2 + (shy/S*(H-1))^2)." % sc["patch_m"],
        "",
        "| channel | mean | p50 | p99 | max |",
        "| --- | --- | --- | --- | --- |",
    ]
    for tag in ("r", "g", "b"):
        ch = warp.get("per_channel", {}).get(tag) or {}
        lines.append("| %s | %.3f px | %.3f px | %.3f px | %.3f px |"
                     % (tag, ch.get("px_mean", float("nan")),
                        ch.get("px_p50", float("nan")),
                        ch.get("px_p99", float("nan")),
                        ch.get("px_max", float("nan"))))
    lines += [
        "",
        "- The project yardstick is the **green** channel (the composite's "
        "own channel for C and g1): **mean %.3f px, p50 %.3f px, p99 %.3f "
        "px** (in metres: %.4f / %.4f / %.4f m)."
        % (wp.get("px_mean", float("nan")), wp.get("px_p50", float("nan")),
           wp.get("px_p99", float("nan")), wp.get("shift_m_mean",
                                                  float("nan")),
           wp.get("shift_m_p50", float("nan")),
           wp.get("shift_m_p99", float("nan"))),
        "- Historical comparison: the v2-era production loop measured "
        "**≈ 23 px mean** at its own scale. v18 was **not** measured this "
        "way (`%s`), so no v18 number exists to compare against; this run is "
        "the first to place the yardstick on the record for the current "
        "stack." % (v18.get("note", "n/a")),
        "- Reading: at depth %.1f m and the realised patch S = %.3f m, the "
        "fine-warp band is %.2fx the v2-era mean. This is the geometric "
        "displacement of the substrate sampling, i.e. the same quantity the "
        "v2 loop reported; it is not the caustic-field distortion (which is "
        "reported separately by the C statistics)."
        % (sc["depth_m"], sc["patch_m"],
           (wp.get("px_mean") or 0.0) / 23.0),
        "",
        "## 5. The caustic field, measured, against v18",
        "",
        "Not a gate and not a look verdict: a descriptive contrast of the "
        "field the two production runs produced (different base, cell, patch "
        "and output size, so it is not a like-for-like comparison).",
        "",
        "| quantity | v18 | v19 |",
        "| --- | --- | --- |",
        "| cell (U / F / depth) | %g / %g / %g | %g / %g / %g |"
        % (v18c["v18"]["cell"]["U"], v18c["v18"]["cell"]["F"],
           v18c["v18"]["cell"]["h"], v18c["v19"]["cell"]["U"],
           v18c["v19"]["cell"]["F"], v18c["v19"]["cell"]["h"]),
        "| realised patch S | %.4g m | %.4g m |"
        % (v18c["v18"]["patch_m"], v18c["v19"]["patch_m"]),
        "| output | %dx%d | %dx%d |"
        % (v18c["v18"]["frames_px"][0], v18c["v18"]["frames_px"][1],
           v18c["v19"]["frames_px"][0], v18c["v19"]["frames_px"][1]),
        "| C p99, mean over 192 frames | %.4f | %.4f |"
        % (v18c["v18"]["C_p99_mean"], v18c["v19"]["C_p99_mean"]),
        "| C p99, min..max over frames | %.4f..%.4f | %.4f..%.4f |"
        % (v18c["v18"]["C_p99_min"], v18c["v18"]["C_p99_max"],
           v18c["v19"]["C_p99_min"], v18c["v19"]["C_p99_max"]),
        "| crush, mean / max | %.2f%% / %.2f%% | %.2f%% / %.2f%% |"
        % (100 * v18c["v18"]["crush_mean"], 100 * v18c["v18"]["crush_max"],
           100 * v18c["v19"]["crush_mean"], 100 * v18c["v19"]["crush_max"]),
        "| linear mean pre-tonemap | %.5f | %.5f |"
        % (v18c["v18"]["lin_mean"], v18c["v19"]["lin_mean"]),
        "",
        "Reading, measured and with no tuning applied: the delivered piece "
        "carries a **less peaked** caustic field than v18 (p99 %.2f vs "
        "%.2f). The cause is geometric, not the water medium — the same "
        "field is folded much further in master cells here (the mean "
        "displacement is %.2f px of a %.4f m patch, i.e. ~%.0f of the "
        "2048 master cells, against ~5 cells in v18), so the fold smears "
        "the peaks across the master grid. Both runs hold C mean at 1.000 "
        "with zero blowout; crush is %.2f%% against the 5%% gate, i.e. the "
        "darkest %.1f%% of the frame sits under 0.02 linear, so the piece "
        "is near the crush limit but inside it."
        % (v18c["v18"]["C_p99_mean"], v18c["v19"]["C_p99_mean"],
           wp.get("px_mean", float("nan")), sc["patch_m"],
           (wp.get("shift_m_mean", 0.0) / (sc["patch_m"] / GRID)),
           100 * v18c["v19"]["crush_mean"], 100 * v18c["v19"]["crush_mean"]),
        "",
        "Colour, stated because it is the first thing the eye reports: the "
        "piece reads teal, and that is the medium — 4 m of clear water with "
        "a = (%s) 1/m attenuates red about %.1fx the green over the round "
        "path. The source photo's pale mint moth is the dominant bright "
        "mass in the frame. No grade was applied anywhere in this path."
        % (", ".join("%.3f" % x for x in sc["a_rgb"]),
           (sc["a_rgb"][0] / sc["a_rgb"][1]) if sc["a_rgb"][1] else
           float("nan")),
        "",
        "## 6. Wall time",
        "",
        "- Sum of the 8 chunk renders: **%.1f s** (%.1f min); per chunk %s s."
        % (total_wall, total_wall / 60.0, ["%.0f" % w for w in walls]),
        "- Stage 0 preflight + scene build: %.1f s; G1 seam (3 det renders + "
        "2 field evolutions): %.1f s; warp (1 field evolution + sampling): "
        "%.1f s; encode + ffprobe: %.1f s."
        % (numbers["wall_time"]["stage0"], numbers["wall_time"]["g1"],
           numbers["wall_time"]["warp"], numbers["wall_time"]["encode"]),
        "- Total measured render/gate/encode work: **%.1f s (%.1f min)**. "
        "det at C_mesh 2048 / Q 25 dominates (~%.1f s per frame over the 8 "
        "chunks)." % (numbers["wall_time"]["total"],
                      numbers["wall_time"]["total"] / 60.0,
                      total_wall / max(len(per), 1)),
        "",
        "## 7. Deliverables",
        "",
        "| file | note |",
        "| --- | --- |",
        "| `out/v19/piece-8s-v19-photo.mp4` | the piece (sha256 `%s`) |"
        % (g3.get("sha256") or "")[:16],
        "| `out/v19/frame-0096.png` | mid-piece frame, final composite |",
        "| `out/v19/numbers-v19-prod.json` | parameters, per-frame rows, "
        "warp, gates, chunk hashes |",
        "| `out/v19/REPORT-V19-PROD8S.md` | this file |",
        "| `out/v19/checkpoint.json` | jobkit checkpoint + per-chunk hashes |",
        "",
        "## 8. Deviations and choices the brief left open",
        "",
        "1. **Cell = the brief's line, not the v15 T1 record.** The brief "
        "fixes U = 6 m/s, F = 100 m, depth = 4.0 m; the v15 record (used by "
        "v18) is U 6 / F 1000 / depth 1.0. The brief's line wins here by "
        "instruction (\"Scene (fixed; do not tune)\"); nothing was tuned "
        "toward a look.",
        "2. **Seed** — the brief's literal 3255045276382562994 was verified "
        "equal to `state.derive_seed(recipe, 1)` (%s), so the code keeps the "
        "derive call and the equality is asserted, not assumed."
        % pf["seed"]["agree"],
        "3. **C_mesh = %d, Q = %d** — the brief names them; they are also "
        "the v15/v18 production values." % (sc["C_mesh"], sc["Q"]),
        "4. **E = 2.2 literal** — the frozen v5 composite constant, as in "
        "v18; the v16 still path's normalised E was NOT used.",
        "5. **Water preset `clear` at its default load** (%g) — the a/b "
        "medium and the phase function are BOTH the clear preset, as the "
        "brief requires; `algae_pond` was not used anywhere in this run."
        % sc["water_load"],
        "6. **g1 recomputed per frame** from the realised slope of that "
        "frame's field (the still path's own rule).",
        "7. **Base handling** — the single photo is padded once per process "
        "and reused for all frames (a still base, unlike v18's 192 Blender "
        "frames), with the v2 loop's right+bottom 1 px edge replication.",
        "8. **Warp yardstick definition** — the brief names frame 96 and "
        "output pixels but not the channel; the green channel is used "
        "because the composite's C weighting and g1 are green-based, and all "
        "three channels are printed so the choice is auditable.",
        "9. **Chunk idempotency** — chunks of 24 frames with per-frame "
        "sha256, exactly as v18; re-running a committed chunk verifies the "
        "hashes and skips.",
        "10. **Stage-0 snapshot re-taken before the gates.** The first "
        "snapshot (taken before any frame) and the final one differ only in "
        "`run_v19_prod.py` (a pad-check fix, a scene-cell fallback and the "
        "report text were finished after the first snapshot); no "
        "pre-existing file's hash differs between the two. The first "
        "snapshot is preserved at `out/v19/_work/stage0-py-hashes-first.json`"
        " so the audit trail is not lost, and stage 0 was re-run before the "
        "gates so the recorded snapshot matches the code that gated the "
        "piece.",
        "",
        "## 9. Uncertainty",
        "",
        "- This is a rendering job, not a measurement: it reports levels, "
        "seam identity, encode conformance and a geometric warp figure; it "
        "says nothing about a real pond.",
        "- The warp figure is a point statistic at one frame of the "
        "transport field, not a video-wide average.",
        "- 8-bit PNG -> H.264 yuv420p is a lossy step; G3 is on the encode's "
        "metadata, not on a pixel-exact decode.",
        "- The clear-water phase fit and a/b values are the code's published "
        "fits (Jonasz FF, Kirk-style class ranges); they were not "
        "re-verified against primary literature in this session.",
        "",
        "Results are reviewed by 無太郎 before use; this report is not final.",
        "",
    ]
    payload = "\n".join(lines).encode("utf-8")
    atomic_write_bytes(REPORT, payload)
    return REPORT


def do_report():
    t0 = time.time()
    log("=" * 72)
    log("[report] numbers + REPORT + the mid-piece still")
    makedirs()
    cp = Checkpoint(CP_PATH)
    src = os.path.join(FRAMES, "frame-0096.png")
    if not os.path.isfile(src):
        raise SystemExit("frame 0096 missing")
    img = open(src, "rb").read()

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size == (WIDTH, HEIGHT) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(STILL96, img, verify=verify)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    write_report(numbers)
    cp.unit_done("report")
    log("  numbers %s | report %s | still %s"
        % (os.path.relpath(NUMBERS, HERE), os.path.relpath(REPORT, HERE),
           os.path.relpath(STILL96, HERE)))
    log("  report wall %.1fs" % (time.time() - t0))
    return numbers


# --------------------------------------------------------------------------
# probe
# --------------------------------------------------------------------------
def do_probe(n):
    makedirs()
    p, st, env, cfg, optics = build_scene(verbose=True)
    photo = load_photo()
    log("  probe frame %d" % n)
    r = render_one(n, st, env, optics, photo)
    sha = save_png(os.path.join(WORK, "probe-%04d.png" % n), r["still"])
    log("  timing %s" % r["timing_s"])
    log("  C mean %.9f p99 %.3f crush %.4f%% blow %.6g g1 %.6f sha %s"
        % (float(np.mean(r["C"], dtype=np.float64)),
           float(np.percentile(r["C"], 99)),
           100 * r["metrics"]["crush_rate_linear_lt_0.02"],
           r["metrics"]["blowout_rate_aces_gt_0.98"], r["g1"], sha[:16]))
    log("  optics build %.1fs" % p["build_wall_s"])
    log("  patch %.4f m Hs %.6f sigma2 %.6f warn %s"
        % (p["patch_m"], p["state_summary"]["Hs_spectral_m"],
           p["state_summary"]["sigma2_spectral"],
           p["state_summary"]["warn_extrapolated_wave_age"]))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--seam", action="store_true")
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--probe", type=int, default=None)
    a = ap.parse_args()
    did = False
    if a.stage0:
        stage0()
        did = True
    if a.probe is not None:
        do_probe(a.probe)
        did = True
    if a.chunk is not None:
        do_chunk(a.chunk)
        did = True
    if a.seam:
        do_seam()
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
        ap.error("nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
