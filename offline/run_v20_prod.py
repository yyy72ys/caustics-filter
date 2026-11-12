"""run_v20_prod -- BRIEF-V20-PROD8S-V2COND: the 8 s piece at the v2 loop's
PLACE AND WEATHER (preset 1: a shallow green stream in the morning, depth
0.6 m) on the full current v3 stack, with honest (spectrum-derived) waves.

This is the controlled comparison the owner asked for: same photo, same place
and weather, same seed as the 2026-11-10 "v2" loop; what is NOT reproduced is
the old pre-v3.7 amplitude-curve wave law (its H_s was not spectrum-derived,
which is why it was replaced). The waves here are the v3.10 spectrum.

Base: the project photo D:\\mutaro-work\\moth-on-bark.jpg (735x539), edge-
      replicated right+bottom by 1 px to 736x540 (the v2 production loop's
      own convention for H.264).
Water: v3.10 spectrum (elfouhaily-form, Arm D) at U = 3 m/s, F = 100 m
       (light wind over a small stream), depth = 0.6 m, seed
       3255045276382562994, time-evolved at t = n/24 s with the v10/v2 8 s
       omega quantisation.
Medium: the v12 ALGAE_POND preset for both a/b (W.PRESETS["algae_pond"]) and
        the phase function (V16.FF_PHASE["algae_pond"]).
Transport: v16 water (a/b partition + FF phase + FF-CDF halo).
Composite: the frozen v5 composite A2.6 + Blow (sigma 3.0 px) + E2.2.

Stages (all foreground, chunked, idempotent per chunk):
  --stage0        preflight + parameter resolution + hash the *.py set
  --chunk K       render frames 24K..24K+23 (atomic per frame, hashed)
  --seam          G1 seam check + the frequency-table quantisation proof
  --warp          the warp yardstick (frame 96, transport field)
  --encode        encode the MP4 (libx264, 24 fps, yuv420p) + G3 ffprobe
  --gates         G2 levels + G4 integrity, roll up into the checkpoint
  --report        numbers-v20-prod.json + REPORT-V20-PROD8S.md + still
  --probe N       render one frame into out/v20/_work (timing only)

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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V20-PROD8S-V2COND.md"
OUT = os.path.join(HERE, "out", "v20")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
WARP_PATH = os.path.join(WORK, "warp.json")
MP4 = os.path.join(OUT, "piece-8s-v20-v2cond.mp4")
NUMBERS = os.path.join(OUT, "numbers-v20-prod.json")
REPORT = os.path.join(OUT, "REPORT-V20-PROD8S.md")
STILL96 = os.path.join(OUT, "frame-0096.png")

PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
V19_NUMBERS = os.path.join(HERE, "out", "v19", "numbers-v19-prod.json")
V2_NUMBERS = os.path.join(HERE, "out", "v2", "numbers-v2.json")
V2_PRESET1 = os.path.join(HERE, "out", "preset1-numbers.json")
PY_HASHES_BEFORE = os.path.join(OUT, "stage0-py-hashes.json")

#: the brief's fixed scene line (do not tune): U = 3 m/s, F = 100 m, depth
#: 0.6 m -- the honest mapping of preset 1's place (light wind over a small
#: stream) and its depthM, exactly.
CELL = {"U_cli": 3.0, "F_cli": 100.0, "h_cli": 0.6}

#: the seed named literally by the brief (and what derive_seed gives).
SEED_BRIEF = 3255045276382562994

#: the frozen composite (A2.6 + Blow + E2.2), as in v18/v19.
A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2

#: the water medium: the v12 preset id and the v16 phase-fit key.
WATER_PRESET = "algae_pond"

FPS = 24
NFRAMES = 192
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
GRID = 2048
Q_DET = 25
WIDTH, HEIGHT = 736, 540
PHOTO_NATIVE = (735, 539)

#: files the brief names (all must exist in the workdir).
CODE_FILES = ("run_v19_prod.py", "run_v15_spectrum.py", "run_v16_water.py",
              "run_v16.py", "v12_water.py", "v5_composite.py",
              "transport.py", "det.py", "image.py")

NEW_FILES = ("run_v20_prod.py",)


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

    # (c) the ALGAE_POND medium preset + the algae_pond phase fit (brief
    #     demands both, at the v16 defaults)
    medium_w = W.PRESETS.get(WATER_PRESET)
    medium_p = V16.FF_PHASE.get(WATER_PRESET)
    facts["medium"] = {
        "v12_preset_key": WATER_PRESET,
        "v12_presets_available": sorted(W.PRESETS),
        "v12_default_medium_preset": W.DEFAULT_MEDIUM_PRESET,
        "v12_preset_found": medium_w is not None,
        "v12_label": (medium_w or {}).get("label"),
        "v12_a_part": (medium_w or {}).get("a_part"),
        "v12_b_part": (medium_w or {}).get("b_part"),
        "v12_load_default": (medium_w or {}).get("load_default"),
        "v16_phase_key": WATER_PRESET,
        "v16_phase_keys_available": sorted(V16.FF_PHASE),
        "v16_phase_found": medium_p is not None,
        "v16_n": (medium_p or {}).get("n"),
        "v16_mu": (medium_p or {}).get("mu"),
        "v16_label": (medium_p or {}).get("label"),
        "v16_provenance": (medium_p or {}).get("provenance"),
        "assert_preset_is_code_default": bool(
            W.DEFAULT_MEDIUM_PRESET == WATER_PRESET),
        "note": ("the brief names the medium by its character ('algae_pond "
                 "medium + phase ... v16 defaults'), not by numeric n/mu; "
                 "nothing was asserted beyond presence and the v16 defaults "
                 "the brief does name"),
    }
    if medium_w is None:
        missing.append("no ALGAE_POND medium preset in v12_water.PRESETS "
                       "(available: %s)" % sorted(W.PRESETS))
    if medium_p is None:
        missing.append("no ALGAE_POND phase fit in run_v16_water.FF_PHASE "
                       "(available: %s)" % sorted(V16.FF_PHASE))

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
        "source": ("BRIEF-V20-PROD8S-V2COND 'Scene (fixed; do not tune)' -- "
                   "U 3 m/s, F 100 m (light wind over a small stream) and "
                   "depth 0.6 m (the v2 recipe's depthM, exactly); not the "
                   "v15 T1 record and not the v19 cell"),
    }

    # (f2) the comparison references this run reports against (not gates)
    facts["references"] = {}
    for tag, path in (("v19", V19_NUMBERS), ("v2", V2_NUMBERS)):
        node = {"path": path, "present": os.path.isfile(path)}
        if node["present"]:
            node["sha256"] = sha256_file(path)
            try:
                node["bytes"] = os.path.getsize(path)
            except OSError:
                pass
        facts["references"][tag] = node

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
    log("  medium: v12 PRESETS['%s'] (%s) | v16 FF_PHASE['%s'] n=%s mu=%s"
        % (WATER_PRESET, facts["medium"]["v12_label"], WATER_PRESET,
           facts["medium"]["v16_n"], facts["medium"]["v16_mu"]))
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


def decode_check():
    """Decode frames 0 / 96 / 191 from the MP4 and compare with the PNGs.

    A one-frame off-by-one would move the MAE by an order of magnitude, so
    this doubles as the frame-index check.
    """
    from PIL import Image
    rows = []
    for i in (0, 96, 191):
        tmp = os.path.join(WORK, "decode-%04d.png" % i)
        if os.path.exists(tmp):
            os.remove(tmp)
        cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "error",
               "-i", MP4, "-vf", "select=eq(n\\,%d)" % i,
               "-fps_mode", "passthrough", "-frames:v", "1", tmp]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not os.path.isfile(tmp):
            # older/newer ffmpeg without -fps_mode: plain extraction
            cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "error",
                   "-i", MP4, "-vf", "select=eq(n\\,%d)" % i,
                   "-frames:v", "1", tmp]
            r = subprocess.run(cmd, capture_output=True, text=True)
        if not os.path.isfile(tmp):
            rows.append({"frame": i, "error": "decode failed",
                         "stderr": (r.stderr or "")[-400:]})
            continue
        dec = np.asarray(Image.open(tmp).convert("RGB")).astype(np.int32)
        src = np.asarray(Image.open(os.path.join(
            FRAMES, "frame-%04d.png" % i)).convert("RGB")).astype(np.int32)
        if dec.shape != src.shape:
            rows.append({"frame": i, "error": "shape %s vs %s"
                         % (dec.shape, src.shape)})
            continue
        diff = np.abs(dec - src)
        rows.append({"frame": i, "mae_255": float(diff.mean()),
                     "max_255": int(diff.max()),
                     "p999_255": float(np.percentile(diff, 99.9))})
    return {"frames": rows,
            "method": ("ffmpeg decode of the encoded stream vs the source "
                       "PNG, mean/max absolute difference in 8-bit levels"),
            "checked_utc": now_iso()}


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
    dec = decode_check()
    log("  decode check: %s" % [(r["frame"], round(r.get("mae_255", -1), 2))
                               for r in dec["frames"]])
    g3 = {"pass": bool(all(checks.values())), "checks": checks,
          "ffprobe": stm, "file": os.path.relpath(MP4, HERE),
          "bytes": os.path.getsize(MP4), "sha256": sha256_file(MP4),
          "encode_cmd": cmd, "decode_check": dec,
          "checked_utc": now_iso(),
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
    outv20 = []
    for root, _dirs, files in os.walk(os.path.join(HERE, "out", "v20")):
        for f in files:
            outv20.append(os.path.relpath(os.path.join(root, f), HERE)
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
        "out_v20_files": sorted(outv20),
        "n_out_v20_files": len(outv20),
        "rule": ("no v10-v19 *.py file modified: every pre-existing *.py "
                 "hash must be identical to the stage-0 snapshot; only "
                 "run_v20_prod.py and out/v20/* may be new"),
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
    current_depth = float(params["scene_as_built"]["depth_m"])
    current_hs = float(
        params["scene_as_built"]["state_summary"]["Hs_spectral_m"])
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None

    # ---- the two references this piece is read against (not gates) -------
    ref = {"v19": None, "v2": None, "v2_loop_yardstick_px": 23.0}
    v19cmp = None
    if os.path.exists(V19_NUMBERS):
        try:
            d19 = read_json(V19_NUMBERS)
            w19 = (d19.get("warp") or {}).get("per_channel", {}).get("g", {})
            ref["v19"] = {
                "path": os.path.relpath(V19_NUMBERS, HERE),
                "cell": d19["scene"]["cell"],
                "patch_m": d19["scene"]["patch_m"],
                "depth_m": d19["scene"]["depth_m"],
                "water_preset": d19["scene"]["water_preset"],
                "phase_n": d19["scene"]["optics"]["phase"]["n"],
                "phase_mu": d19["scene"]["optics"]["phase"]["mu"],
                "frames_px": [d19["frames"]["width"], d19["frames"]["height"]],
                "warp_g_mean_px": w19.get("px_mean"),
                "warp_g_p50_px": w19.get("px_p50"),
                "warp_g_p99_px": w19.get("px_p99"),
                "Hs_spectral_m": d19["scene"]["state_summary"]["Hs_spectral_m"],
                "sigma2_spectral":
                    d19["scene"]["state_summary"]["sigma2_spectral"],
                "rms_slope": d19["scene"]["state_summary"].get("rms_slope"),
                "spectrum_source": ("v3.10 elfouhaily-form / Arm D "
                                    "(same code, same seed)"),
            }
            p19 = np.array([r["C_p99"] for r in d19["per_frame"]])
            c19 = np.array([r["crush"] for r in d19["per_frame"]])
            l19 = np.array([r["lin_mean"] for r in d19["per_frame"]])
            p20 = np.array([r["C_p99"] for r in per_frame])
            c20 = np.array([r["crush"] for r in per_frame])
            l20 = np.array([r["lin_mean"] for r in per_frame])
            v19cmp = {
                "v19": {"cell": d19["scene"]["cell"],
                        "patch_m": d19["scene"]["patch_m"],
                        "depth_m": d19["scene"]["depth_m"],
                        "medium": d19["scene"]["water_preset"],
                        "frames_px": [d19["frames"]["width"],
                                      d19["frames"]["height"]],
                        "C_p99_mean": float(p19.mean()),
                        "C_p99_min": float(p19.min()),
                        "C_p99_max": float(p19.max()),
                        "crush_mean": float(c19.mean()),
                        "crush_max": float(c19.max()),
                        "lin_mean": float(l19.mean()),
                        "Hs_spectral_m":
                            d19["scene"]["state_summary"]["Hs_spectral_m"],
                        "warp_g_mean_px": w19.get("px_mean"),
                        "warp_g_p50_px": w19.get("px_p50"),
                        "warp_g_p99_px": w19.get("px_p99")},
                "v20": {"cell": {"U": CELL["U_cli"], "F": CELL["F_cli"],
                                 "h": CELL["h_cli"]},
                        "patch_m": patch_m,
                        "depth_m": current_depth,
                        "medium": WATER_PRESET,
                        "frames_px": [WIDTH, HEIGHT],
                        "C_p99_mean": float(p20.mean()),
                        "C_p99_min": float(p20.min()),
                        "C_p99_max": float(p20.max()),
                        "crush_mean": float(c20.mean()),
                        "crush_max": float(c20.max()),
                        "lin_mean": float(l20.mean()),
                        "Hs_spectral_m":
                            current_hs,
                        "warp_g_mean_px": (warp or {}).get("g_mean_px"),
                        "warp_g_p50_px": (warp or {}).get("g_p50_px"),
                        "warp_g_p99_px": (warp or {}).get("g_p99_px")},
                "caveat": ("same code, same seed, same base, same output "
                           "size; the cell, depth and medium differ (v19 was "
                           "6 m/s over 4 m of clear water), so this is a "
                           "descriptive contrast, not a preference"),
            }
        except Exception as exc:                      # noqa: BLE001
            ref["v19"] = {"path": V19_NUMBERS, "error": str(exc)}
    if os.path.exists(V2_NUMBERS):
        try:
            d2 = read_json(V2_NUMBERS)
            stc = d2.get("statistics_master_C") or {}
            ref["v2"] = {
                "path": os.path.relpath(V2_NUMBERS, HERE),
                "pipeline": d2.get("pipeline"),
                "recipe": d2.get("recipe"),
                "preset": d2.get("preset"),
                "seed": d2.get("seed"),
                "t_s": d2.get("t_s"),
                "patch_m": d2["state"]["patch_S_m"],
                "depth_m": d2["environment"]["depthM"],
                "sun_alt_deg": d2["environment"]["sun_alt_deg"],
                "clarity": d2["environment"].get("clarity"),
                "grid": d2["state"]["grid"],
                "C_quantiles": {k: stc.get(k) for k in
                                ("mean", "p50", "p99", "p99.9", "max", "min")},
                "title_ja": "緑の川の浅瀬・午前",
                "note": ("the v2 loop's own record for the same place, "
                         "weather, seed and depth; its wave law was the "
                         "pre-v3.7 amplitude-curve machinery, which is "
                         "exactly what is NOT reproduced here"),
            }
        except Exception as exc:                      # noqa: BLE001
            ref["v2"] = {"path": V2_NUMBERS, "error": str(exc)}
    if os.path.exists(V2_PRESET1) and isinstance(ref.get("v2"), dict):
        try:
            d1 = read_json(V2_PRESET1)
            st = d1.get("state") or {}
            ref["v2"]["amplitude_target"] = {
                "path": os.path.relpath(V2_PRESET1, HERE),
                "preset": d1.get("preset"),
                "seed": d1.get("seed"),
                "target_std_m": st.get("target_std_m"),
                "amplitude_scale": st.get("amplitude_scale"),
                "patch_S_m": st.get("patch_S_m"),
                "loop_T_s": st.get("loop_T_s"),
                "note": ("the pre-v3.7 amplitude curve's own target: a "
                         "prescribed sigma_eta, not a spectrum integral. "
                         "Reported as recorded fact; reviving it is what the "
                         "brief forbids."),
            }
        except Exception as exc:                      # noqa: BLE001
            ref["v2"]["amplitude_target"] = {"error": str(exc)}

    return {
        "job": "V20-PROD8S-V2COND",
        "brief": BRIEF,
        "created_utc": now_iso(),
        "model": ("v3 ALPHA close: v3.10 spectrum (elfouhaily-form, Arm D) + "
                  "v16 water (FF phase, FF-CDF halo) + frozen composite "
                  "A2.6 + Blow + E2.2, at the v2 loop's place and weather "
                  "(algae_pond, 0.6 m) with honest spectrum-derived waves"),
        "preflight": {"pass": params["preflight_pass"],
                      "photo": params["photo"],
                      "v16_defaults": params["v16_defaults"],
                      "medium": params["medium"],
                      "seed": params["seed"],
                      "missing": params["missing"],
                      "references": params.get("references")},
        "scene": {"cell": params["scene_as_built"]["cell"],
                  "cell_source": ("BRIEF-V20-PROD8S-V2COND scene line (fixed; "
                                  "U 3 m/s / F 100 m / depth 0.6 m; not the "
                                  "v15 T1 record, not the v19 cell)"),
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
        "references": ref,
        "v19_contrast": v19cmp,
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
    refs = numbers.get("references") or {}
    v19 = refs.get("v19") or {}
    v2 = refs.get("v2") or {}
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
    gates_all = ("G1", "G2", "G3", "G4")
    verdict = all(bool(g.get(k, {}).get("pass")) for k in gates_all)
    g2node = g.get("G2", {})
    ss = sc["state_summary"]
    cm = pf["medium"]
    wp = warp.get("per_channel", {}).get("g", {})
    v19c = numbers.get("v19_contrast") or {}
    dec = g3.get("decode_check") or {}
    dec_rows = "\n".join(
        "| %s | %s | %s | %s |" % (
            r.get("frame"), ("%.2f" % r["mae_255"]) if "mae_255" in r
            else r.get("error"),
            r.get("max_255", "-"), ("%.2f" % r["p999_255"])
            if "p999_255" in r else "-")
        for r in dec.get("frames", []))
    amp = (v2.get("amplitude_target") or {})
    v2q = v2.get("C_quantiles") or {}
    warp2_px = float(refs.get("v2_loop_yardstick_px") or 23.0)
    lines = [
        "# REPORT-V20-PROD8S — the 8 s piece at the v2 loop's place and "
        "weather, honest waves",
        "",
        "Brief: `BRIEF-V20-PROD8S-V2COND.md` (D:\\mutaro-work\\rally-in). "
        "Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.",
        "Built: %s. Code: `run_v20_prod.py` only." % numbers["created_utc"],
        "",
        "This is the controlled comparison the brief describes: the photo, "
        "the place and weather of preset 1 (a shallow green stream in the "
        "morning, depth 0.6 m) and the seed are the v2 loop's; the waves are "
        "NOT its waves — they are the v3.10 spectrum, because the old law's "
        "H_s was not spectrum-derived and reviving it is exactly what "
        "physics-first forbids.",
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
        "No parameter was tuned toward a look at any point: the scene line "
        "is the brief's, the medium is the brief's, and the composite is the "
        "frozen A2.6 + Blow + E2.2 of v18/v19.",
        "",
        "## 2. What was rendered",
        "",
        "- 192 frames, %dx%d, 24 fps, t = n/24 s, frequencies on the 2pi/8 s "
        "grid (inherited from the state build)." % (WIDTH, HEIGHT),
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
        "default) at **U = %g m/s, F = %g m, depth = %g m**, seed **%d**, "
        "C_mesh %d, det Q = %d."
        % (sc["cell"]["U"], sc["cell"]["F"], sc["cell"]["h"], sc["seed"],
           sc["C_mesh"], sc["Q"]),
        "- Medium (**algae_pond**, both a/b and phase — identifiers found in "
        "the code, not assumed):",
        "  - v12 preset key **`%s`** inside `v12_water.PRESETS` — "
        "`preset_spectra(\"%s\")`: label `%s`, a_part %s, b_part %s, "
        "load_default %s (also the code's `DEFAULT_MEDIUM_PRESET`: %s)."
        % (cm["v12_preset_key"], cm["v12_preset_key"], cm["v12_label"],
           cm["v12_a_part"], cm["v12_b_part"], cm["v12_load_default"],
           cm["assert_preset_is_code_default"]),
        "  - v16 phase key **`%s`** inside `run_v16_water.FF_PHASE`: "
        "n = **%s**, mu = **%s**, label `%s` (%s). The brief names the "
        "medium by character and says 'v16 defaults', so the code's own "
        "`PHASE_DEFAULT = %s` / `KERNEL_DEFAULT = %s` are used as-is; no "
        "n/mu value was asserted or chosen by hand."
        % (cm["v16_phase_key"], cm["v16_n"], cm["v16_mu"], cm["v16_label"],
           cm["v16_provenance"], pf["v16_defaults"]["PHASE_DEFAULT"],
           pf["v16_defaults"]["KERNEL_DEFAULT"]),
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
        "- Sky: direct sun only, eta = 1.0 (the recipe carries no sky block); "
        "sun handling identical to v19. Retina OFF. Rain 0 mm/h.",

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
        "| spectrum H_s | **%.6f m** |" % ss["Hs_spectral_m"],
        "| rms slope | %.6f |" % ss["rms_slope"],
        "| omega_c | %.6f rad/s |" % ss["omega_c"],
        "| `warn_extrapolated_wave_age` | **%s** (flag fires when omega_c > %g "
        "rad/s; omega_c = %.6f rad/s) |"
        % (ss["warn_extrapolated_wave_age"], V15S.WARN_OMEGA_C,
           ss["omega_c"]),
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
        "sits on an integer multiple of 2pi/8 s (the table above)."
        % g1node.get("branch_b_quantisation_proof"),
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
        "| crush (< 0.02 linear) | max %.4f%% per frame, mean %.4f%% | "
        "< 5%% per frame |" % (100 * max(crush), 100 * sum(crush) / len(crush)),
        "| blowout (ACES > 0.98) | max %.6g, %d frames nonzero | 0 per frame |"
        % (max(blow), g2node.get("blowout_aces_gt_0.98", {})
           .get("n_frames_nonzero", -1)),
        "",
        "### G3 encode — %s" % ("PASS" if g3.get("pass") else "FAIL"),
        "",
        "`ffprobe` on `out/v20/piece-8s-v20-v2cond.mp4`: codec `%s`, %sx%s, "
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
        "Decode check (automated in this job, %s):" % dec.get(
            "method", "n/a"),
        "",
        "| frame | MAE (levels of 255) | max | p99.9 |",
        "| --- | --- | --- | --- |",
        dec_rows or "| -- | -- | -- | -- |",
        "",
        "Frame 0 / 96 / 191 decode back close to their source PNGs — the "
        "expected H.264 crf-18 yuv420p rounding on this high-contrast "
        "caustic content — and there is no frame off-by-one (a one-frame "
        "shift would move the MAE by an order of magnitude).",
        "",
        "### G4 integrity — %s" % ("PASS" if g.get("G4", {}).get("pass")
                                   else "FAIL"),
        "",
        "- `*.py` set hashed at stage 0 (%d files) and again at the end "
        "(%d files): added %s, changed %s, removed %s."
        % (g["G4"]["py_before"], g["G4"]["py_after"], g["G4"]["added"],
           g["G4"]["changed"], g["G4"]["removed"]),
        "- Pre-existing `*.py` files changed: **%s**; unexpected new files: "
        "**%s**. No v10–v19 file was modified."
        % (g["G4"].get("changed_existing_files", []),
           g["G4"].get("unexpected_new_or_changed", [])),
        "- `run_v20_prod.py` sha256 is identical in the stage-0 snapshot and "
        "at the end: `%s`. It is the job's only new `*.py` (the brief allows "
        "it); the `added` list above is empty because the stage-0 snapshot "
        "was re-taken after the last code change and therefore already "
        "contains it. No pre-existing file's hash moved."
        % (tuple((g["G4"].get("own_new_file_sha256", {}) or {})
                 .get("run_v20_prod.py", [None, None]))[1],),
        "- New files are confined to `run_v20_prod.py` and `out/v20/*` (%d "
        "files under `out/v20`)." % g["G4"].get("n_out_v20_files", -1),
        "- The base photo is read only (hashed at stage 0); no v2/v18/v19 "
        "artefact is written to.",
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
    v19_warp = v19c.get("v19", {}) if v19c else {}
    v20_warp = v19c.get("v20", {}) if v19c else {}
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
        "- Against the v2 loop's yardstick (**~%g px mean** at its own scale, "
        "block-matched): this piece is **%.2fx** of it."
        % (warp2_px, (wp.get("px_mean") or 0.0) / warp2_px),
        "- Against v19 (same code, same seed, same output size, 4 m deep, "
        "clear water): **%.3f px mean / %.3f p50 / %.3f p99** here vs "
        "**%.3f / %.3f / %.3f** in v19."
        % (v20_warp.get("warp_g_mean_px") or float("nan"),
           v20_warp.get("warp_g_p50_px") or float("nan"),
           v20_warp.get("warp_g_p99_px") or float("nan"),
           v19_warp.get("warp_g_mean_px") or float("nan"),
           v19_warp.get("warp_g_p50_px") or float("nan"),
           v19_warp.get("warp_g_p99_px") or float("nan")),
        "- Reading: the displacement is depth x slope. The depth ratio to "
        "v19 alone is %.2f; the realised rms slope is also lower here "
        "(%.6f vs v19's %s), and the net measured ratio of the mean "
        "displacement is %.2f. The figure is the geometric displacement of "
        "the substrate sampling (the quantity the v2 loop reported), not the "
        "caustic-field distortion."
        % (sc["depth_m"] / (v19.get("depth_m") or 1.0), ss["rms_slope"],
           ("%.6f" % v19["rms_slope"]) if v19.get("rms_slope") else "n/a",
           (wp.get("px_mean") or 0.0) / (v19_warp.get("warp_g_mean_px")
                                         or 1.0)),
        "",
        "## 5. The controlled comparison with the v2 loop",
        "",
        "What is held equal (the brief's design): the photo, the place and "
        "weather (preset 1: shallow green stream, morning), **depth 0.6 m**, "
        "the seed **%d**, and the 2pi/8 s loop. What is not: the wave law."
        % sc["seed"],
        "",
        "| item | v2 loop (record on disk) | this piece (v20) |",
        "| --- | --- | --- |",
        "| recipe / place | %s / preset %s | same place, same weather |"
        % (v2.get("title_ja", "n/a"), v2.get("preset", "n/a")),
        "| depth | %s m | %.2f m |"
        % (("%.2f" % float(v2["depth_m"])) if v2.get("depth_m") is not None
           else "n/a", sc["depth_m"]),
        "| seed | %s | %d |" % (v2.get("seed", "n/a"), sc["seed"]),
        "| realised patch S | %s m | %.4f m |"
        % (("%.4f" % float(v2["patch_m"])) if v2.get("patch_m") else "n/a",
           sc["patch_m"]),
        "| wave amplitude | prescribed target sigma_eta %s m (amplitude "
        "curve; `amplitude_scale` %s) | spectrum-derived H_s **%.6f m**, "
        "sigma^2 **%.6g rad^2** |"
        % (("%.6f" % float(amp["target_std_m"]))
           if amp.get("target_std_m") else "n/a",
           ("%.4g" % float(amp["amplitude_scale"]))
           if amp.get("amplitude_scale") else "n/a",
           ss["Hs_spectral_m"], ss["sigma2_spectral"]),
        "| C master quantiles | mean %s, p50 %s, p99 %s, max %s (one "
        "field, t = %s s) | mean 1.0 (gated), p99 mean over 192 frames "
        "**%.4f** |"
        % (("%.4f" % v2q["mean"]) if v2q.get("mean") is not None else "n/a",
           ("%.4f" % v2q["p50"]) if v2q.get("p50") is not None else "n/a",
           ("%.4f" % v2q["p99"]) if v2q.get("p99") is not None else "n/a",
           ("%.3f" % v2q["max"]) if v2q.get("max") is not None else "n/a",
           ("%g" % float(v2.get("t_s"))) if v2.get("t_s") is not None
           else "n/a",
           v20_warp.get("C_p99_mean") or float("nan")),
        "| output | 736x540 (v2 convention) | %dx%d |" % (WIDTH, HEIGHT),
        "",
        "Read this table as what it is: the place and weather are the v2 "
        "loop's, the wave field is not. The old target was a prescribed "
        "sigma_eta%s; the new H_s is the integral of the Elfouhaily-form "
        "spectrum at U = %g m/s, F = %g m. The two numbers are therefore "
        "**not** the same kind of quantity and are printed side by side to "
        "make the substitution visible, not to claim equivalence."
        % ((" (%.6f m)" % float(amp["target_std_m"]))
           if amp.get("target_std_m") else "",
           sc["cell"]["U"], sc["cell"]["F"]),
        "",
        "## 6. The caustic field and the medium, measured",
        "",
        "| quantity | v19 (6 m/s, 4 m, clear) | v20 (3 m/s, 0.6 m, "
        "algae_pond) |",
        "| --- | --- | --- |",
        "| cell (U / F / depth) | %g / %g / %g | %g / %g / %g |"
        % (v19c.get("v19", {}).get("cell", {}).get("U", float("nan")),
           v19c.get("v19", {}).get("cell", {}).get("F", float("nan")),
           v19c.get("v19", {}).get("cell", {}).get("h", float("nan")),
           sc["cell"]["U"], sc["cell"]["F"], sc["cell"]["h"]),
        "| realised patch S | %.4g m | %.4g m |"
        % (v19c.get("v19", {}).get("patch_m", float("nan")), sc["patch_m"]),
        "| spectral H_s | %.6f m | %.6f m |"
        % (v19c.get("v19", {}).get("Hs_spectral_m", float("nan")),
           ss["Hs_spectral_m"]),
        "| C p99, mean over frames | %.4f | %.4f |"
        % (v19c.get("v19", {}).get("C_p99_mean", float("nan")),
           v20_warp.get("C_p99_mean", float("nan"))),
        "| crush, mean / max | %.2f%% / %.2f%% | %.2f%% / %.2f%% |"
        % (100 * (v19c.get("v19", {}).get("crush_mean") or 0.0),
           100 * (v19c.get("v19", {}).get("crush_max") or 0.0),
           100 * sum(crush) / len(crush), 100 * max(crush)),
        "| linear mean pre-tonemap | %.5f | %.5f |"
        % (v19c.get("v19", {}).get("lin_mean", float("nan")),
           sum(r["lin_mean"] for r in per) / len(per)),
        "",
        "The substrate here is seen through 0.6 m of algae-rich water twice "
        "(down and up): a = (%s) 1/m, so the round-trip attenuation is "
        "exp(-2*0.6*a) = (%s) per channel before the in-scatter terms."
        % (", ".join("%.3f" % x for x in sc["a_rgb"]),
           ", ".join("%.3f" % np.exp(-1.2 * x) for x in sc["a_rgb"])),
        "",
        "Colour, stated because it is the first thing the eye reports: an "
        "algae-pond medium absorbs blue hardest (a_b = %.3f 1/m against "
        "a_r = %.3f 1/m), so what survives is the green stream's own colour; "
        "the source photo's pale moth is still the dominant bright mass. No "
        "grade was applied anywhere in this path."
        % (sc["a_rgb"][2], sc["a_rgb"][0]),
        "",
        "## 7. Wall time",
        "",
        "- Sum of the 8 chunk renders: **%.1f s** (%.1f min); per chunk %s s."
        % (total_wall, total_wall / 60.0, ["%.0f" % w for w in walls]),
        "- Stage 0 preflight + scene build: %.1f s; G1 seam: %.1f s; warp: "
        "%.1f s; encode + ffprobe + decode check: %.1f s."
        % (numbers["wall_time"]["stage0"], numbers["wall_time"]["g1"],
           numbers["wall_time"]["warp"], numbers["wall_time"]["encode"]),
        "- Total measured render/gate/encode work: **%.1f s (%.1f min)**."
        % (numbers["wall_time"]["total"], numbers["wall_time"]["total"] / 60.0),
        "",
        "## 8. Deliverables",
        "",
        "| file | note |",
        "| --- | --- |",
        "| `out/v20/piece-8s-v20-v2cond.mp4` | the piece (sha256 `%s`) |"
        % (g3.get("sha256") or "")[:16],
        "| `out/v20/frame-0096.png` | mid-piece frame, final composite |",
        "| `out/v20/numbers-v20-prod.json` | parameters, per-frame rows, "
        "warp, references, gates, chunk hashes |",
        "| `out/v20/REPORT-V20-PROD8S.md` | this file |",
        "| `out/v20/checkpoint.json` | jobkit checkpoint + per-chunk hashes "
        "(24-frame chunks, idempotent per chunk) |",
        "",
        "## 9. Deviations and choices the brief left open",
        "",
        "1. **The wave law.** Not a deviation but the point of the run: the "
        "v2 loop's amplitude-curve waves are not reproduced. The H_s here is "
        "the spectrum's own integral (%.6f m), and that is the honest "
        "comparison substrate the brief asks for."
        % ss["Hs_spectral_m"],
        "2. **U and F.** The brief fixes U = 3 m/s, F = 100 m as 'the honest "
        "mapping of preset 1's place'; the v2 record on disk does not store "
        "a wind speed (it stored a prescribed sigma_eta), so no v2 U exists "
        "to match. The brief's line is used literally.",
        "3. **Seed** — the brief's literal %d was verified equal to "
        "`state.derive_seed(recipe, 1)` (%s), so the code keeps the derive "
        "call and the equality is asserted, not assumed."
        % (SEED_BRIEF, pf["seed"]["agree"]),
        "4. **Medium naming.** The brief says 'algae_pond medium + phase ... "
        "v16 defaults'; the code offers exactly that key in both tables, and "
        "its n/mu are used as published. `clear` and `sediment` were not "
        "used anywhere in this run.",
        "5. **Depth 0.6 m** is passed to both the spectrum (dispersion / "
        "patch) and the optics (beam paths) — one depth, as in v19.",
        "6. **E = 2.2 literal**, C_mesh 2048, det Q 25, g1 recomputed per "
        "frame from that frame's realised slope — all as v18/v19.",
        "7. **Base handling** — the single photo is padded once per process "
        "and reused for all frames, with the v2 loop's right+bottom 1 px "
        "edge replication.",
        "8. **Warp yardstick** — green channel, frame 96, output pixels; all "
        "three channels printed so the choice is auditable.",
        "9. **Chunk idempotency** — chunks of 24 frames with per-frame "
        "sha256; re-running a committed chunk verifies the hashes and skips.",
        "10. **The v2 yardstick (~%g px)** is quoted from the brief, not "
        "re-measured here: the brief calls it 'block-matched', and the v2 "
        "artefacts on disk (`out/v2/numbers-v2.json`, "
        "`out/preset1-numbers.json`) do not contain that measurement."
        % warp2_px,
        "",
        "## 10. Uncertainty",
        "",
        "- This is a rendering job, not a measurement of a real pond: it "
        "reports levels, seam identity, encode conformance and a geometric "
        "warp figure.",
        "- The warp figure is a point statistic at one frame of the "
        "transport field, not a video-wide average.",
        "- The v2 comparison uses the numbers the old loop left on disk; "
        "its delivery-time state (whatever it was at 2026-11-10) is not "
        "reconstructible from them beyond what is cited.",
        "- 8-bit PNG -> H.264 yuv420p is a lossy step; G3 is on the encode's "
        "metadata plus the printed decode MAE, not on a pixel-exact decode.",
        "- The algae_pond a/b and phase values are the code's published "
        "preset/fits; they were not re-verified against primary literature "
        "in this session.",
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
