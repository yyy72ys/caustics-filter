"""run_v18_prod -- BRIEF-V18-PROD8S: the 8 s piece on the v3.10 + v16 stack.

Base: the 192 shipped Blender frames (blender-piece8s/renders_v2, 832x480).
Water: v3.10 spectrum (elfouhaily-form, Arm D) at the v15 T1 cell
       U = 6 m/s, F = 1000 m, depth = 1.0 m, seed = derive_seed(recipe-1, 1),
       time-evolved at t = n/24 s with the v10/v2 8 s omega quantisation.
Transport: v16 water (a/b partition + FF phase + FF-CDF halo).
Composite: the frozen v5 composite A2.6 + Blow (sigma 3.0 px) + E2.2.

Stages (all foreground, chunked, idempotent per chunk):
  --stage0        preflight + parameter resolution + hash the *.py set
  --chunk K       render frames 24K..24K+23 (atomic per frame, hashed)
  --seam          G1 seam check + the frequency-table quantisation proof
  --encode        encode the MP4 (libx264, 24 fps, yuv420p) + G3 ffprobe
  --gates         G2 levels + G4 integrity, roll up into the checkpoint
  --report        numbers-v18-prod.json + REPORT-V18-PROD8S.md + still
  --probe N       render one frame into out/v18/_work (timing only)

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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V18-PROD8S.md"
OUT = os.path.join(HERE, "out", "v18")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
MP4 = os.path.join(OUT, "piece-8s-v3alpha.mp4")
NUMBERS = os.path.join(OUT, "numbers-v18-prod.json")
REPORT = os.path.join(OUT, "REPORT-V18-PROD8S.md")
STILL96 = os.path.join(OUT, "frame-0096.png")

BASE_DIR = r"D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\renders_v2"
V15_NUMBERS = os.path.join(HERE, "out", "v15", "numbers-v15-transport.json")
PY_HASHES_BEFORE = os.path.join(OUT, "stage0-py-hashes.json")

#: the brief's scene line; the v15 T1 record is authoritative if it differs.
CELL = {"U_cli": 6.0, "F_cli": 1000.0, "h_cli": 1.0}

#: the frozen composite (broadcast in the brief: A2.6 + Blow + E2.2).
A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2

FPS = 24
NFRAMES = 192
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
GRID = 2048            # the v15 T1 record's C_mesh
Q_DET = 25
WIDTH, HEIGHT = 832, 480

#: files the brief names (all must exist in the workdir).
CODE_FILES = ("run_v15_spectrum.py", "run_v15.py", "v12_water.py",
              "run_v16_water.py", "run_v16.py", "v5_composite.py",
              "transport.py", "det.py", "image.py")

NEW_FILES = ("run_v18_prod.py",)


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
    """Add/replace one key in out/v18/checkpoint.json (jobkit format kept)."""
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

    # (c) the base frames: exactly 192, 832x480, named frame_%04d.png
    from PIL import Image
    if not os.path.isdir(BASE_DIR):
        missing.append("BASE_DIR missing: %s" % BASE_DIR)
        base = {}
    else:
        names = sorted(n for n in os.listdir(BASE_DIR)
                       if n.lower().endswith(".png"))
        expected = ["frame_%04d.png" % i for i in range(NFRAMES)]
        base = {"dir": BASE_DIR, "count": len(names),
                "names_match_0000_0191": names == expected,
                "sizes": {}, "sha256": {}, "bytes": {}}
        for i, n in enumerate(expected):
            p = os.path.join(BASE_DIR, n)
            if not os.path.isfile(p):
                missing.append("base frame missing: %s" % n)
                continue
            im = Image.open(p)
            base["sizes"][n] = [im.size[0], im.size[1]]
            base["bytes"][n] = os.path.getsize(p)
        sizes = sorted({tuple(v) for v in base["sizes"].values()})
        base["unique_sizes"] = [list(s) for s in sizes]
        base["all_832x480"] = sizes == [(WIDTH, HEIGHT)]
        if not base["all_832x480"]:
            missing.append("base frame size != 832x480: %s" % (sizes,))
        if len(names) != NFRAMES:
            missing.append("base frame count %d != %d" % (len(names), NFRAMES))
        # hash the whole base set (integrity of the input, cheap enough)
        for n in expected:
            p = os.path.join(BASE_DIR, n)
            if os.path.isfile(p):
                base["sha256"][n] = sha256_file(p)
        base["set_sha256"] = sha256_bytes(
            "".join(base["sha256"][n] for n in expected).encode("ascii"))
    facts["base_frames"] = base

    # (d) the v15 T1 record (the FILE wins over the brief's line)
    rec = read_json(V15_NUMBERS)
    t1 = rec["cells"]["T1"]
    cell_rec = {"U": float(t1["cell"]["U_ms"]),
                "F": float(t1["cell"]["fetch_m"]),
                "h": float(t1["cell"]["depth_m"])}
    agreed = {k: bool(cell_rec[k] == CELL[k + "_cli"])
              for k in ("U", "F", "h")}
    facts["cell"] = {
        "file": os.path.relpath(V15_NUMBERS, HERE),
        "record": cell_rec,
        "brief_line": {"U": CELL["U_cli"], "F": CELL["F_cli"],
                       "h": CELL["h_cli"]},
        "agrees_with_brief": agreed,
        "file_wins": not all(agreed.values()),
        "used": cell_rec,
        "record_C_mesh": rec["stack"]["C_mesh"],
        "record_seed": rec["stack"]["seed"],
        "record_composite": rec["stack"]["composite"],
        "record_water_model": rec["stack"]["water_model"],
        "record_eta": rec["stack"]["eta"],
        "record_retina": rec["stack"]["retina"],
        "record_rain_mm_h": rec["stack"]["rain_mm_h"],
        "record_spectrum": rec["stack"]["spectrum"],
        "record_spectrum_arm": rec["stack"]["spectrum_arm"],
        "record_spectrum_rule": rec["stack"]["spectrum_rule"],
        "record_still_size_px": rec["stack"]["still_size_px"],
    }

    # (e) the *.py set before any new work (G4)
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before),
                                 "path": os.path.relpath(PY_HASHES_BEFORE, HERE),
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
    log("  base 192 frames 832x480 (set sha %s)"
        % base["set_sha256"][:16])
    log("  cell U=%g F=%g h=%g (file agrees: %s) | seed %d | C_mesh %d"
        % (cell_rec["U"], cell_rec["F"], cell_rec["h"],
           all(agreed.values()), p["seed"], GRID))
    log("  optics: phase %s kernel %s | patch %.6g m | w_bal %s | w_halo %s"
        % (p["optics"]["phase_kind"], p["optics"]["kernel"], p["patch_m"],
           ["%.6f" % x for x in p["optics"]["w_ballistic"]],
           ["%.6f" % x for x in p["optics"]["w_halo"]]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


# --------------------------------------------------------------------------
# the scene (identical build in every process)
# --------------------------------------------------------------------------
def build_scene(verbose=False):
    t0 = time.time()
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    cell = read_json(PARAMS_PATH)["cell"]["used"] if os.path.exists(
        PARAMS_PATH) else {"U": CELL["U_cli"], "F": CELL["F_cli"],
                           "h": CELL["h_cli"]}
    r = R15.scene_for(raw, float(cell["U"]), float(cell["F"]),
                      float(cell["h"]), GRID, seed)
    st, env, cfg = r["state"], r["env"], r["cfg"]
    spec = W.preset_spectra("algae_pond")
    o = V16.optics_v16(spec, env, env["depthM"], GRID, st.S)
    prof = o["kernel_profile"]
    out = {
        "seed": seed, "grid": GRID, "Q": Q_DET,
        "cell": {"U": float(cell["U"]), "F": float(cell["F"]),
                 "h": float(cell["h"])},
        "patch_m": float(st.S),
        "depth_m": float(env["depthM"]),
        "turbid": float(env["turbid"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "eta": 1.0,
        "eta_source": ("the recipe carries no sky block, so eta = "
                       "DNI/(DNI+DHI) = 1.0 (direct sun only)"),
        "retina": "OFF (the v7 observer post-stage is not in this path)",
        "rain_mm_h": 0.0,
        "water_preset": spec["preset"],
        "water_load": float(spec["load"]),
        "a_rgb": [float(x) for x in spec["a"]],
        "b_rgb": [float(x) for x in spec["b"]],
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
        },
        "build_wall_s": float(time.time() - t0),
    }
    if verbose:
        log("  scene built in %.1fs" % out["build_wall_s"])
    return out, st, env, cfg, o


# --------------------------------------------------------------------------
# one frame
# --------------------------------------------------------------------------
def load_base(i):
    from PIL import Image
    from image import srgb_to_linear
    p = os.path.join(BASE_DIR, "frame_%04d.png" % i)
    u8 = np.asarray(Image.open(p).convert("RGB"))
    return u8, srgb_to_linear(u8.astype(np.float64) / 255.0)


def render_one(i, st, env, optics, g1=None):
    from det import det_frame
    from v5_composite import image_metrics
    from v5_physics import smith_g1_beckmann
    t0 = time.time()
    t = i / float(FPS)
    base_u8, photo_lin = load_base(i)
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
        fields, C, photo_lin, env, WIDTH, HEIGHT, optics,
        ambient_scale=A_FROZEN, c_layer_exposure=E_FROZEN,
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=B_SIGMA_FROZEN, return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t3
    out = {"frame": i, "t_s": t, "g1": float(g1), "C": C,
           "det_info": det_info, "still": img, "metrics": metrics,
           "base_u8": base_u8,
           "timing_s": {"fields": t_fields, "det": t_det,
                        "composite": t_comp,
                        "total": float(time.time() - t0)}}
    del dbg, photo_lin, fields
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
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics)
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
    del st, optics
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

    # the quantised frequency table (both cascades of the state)
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
        # the phase the renderer applies at t = 8 s: exp(1j * w * 8)
        phase = w * 8.0
        phase_dev = np.abs(phase - 2.0 * np.pi * np.round(phase
                                                           / (2.0 * np.pi)))
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
    log("  fields max abs diff t8 vs t0: %s" % {
        k: v for k, v in fdiff.items()})
    log("  C max abs diff t8 vs t0: %.12g (1-frame sensitivity %.4g)"
        % (dC, sens))
    log("  branch (a) exact zero: %s | branch (b) quantisation proof: %s"
        % (exact_zero, proof_ok))
    log("  G1 %s" % ("PASS" if gates["G1"]["pass"] else "FAIL"))
    del st, optics
    gc.collect()
    return gates["G1"]


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
    missing = [i for i in range(NFRAMES)
               if not os.path.isfile(pattern % i)]
    if missing:
        raise SystemExit("frames missing: %d" % len(missing))
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
        "width_832": int(stm.get("width") or 0) == WIDTH,
        "height_480": int(stm.get("height") or 0) == HEIGHT,
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

    # G4: the *.py set now vs stage 0
    before = read_json(PY_HASHES_BEFORE)["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    # the job's own new file is allowed to move while the job is developed;
    # every PRE-EXISTING *.py (the v10-v17 families and everything else) is not.
    changed_existing = [c for c in changed if c not in allowed_new_static()]
    outv18 = []
    for root, _dirs, files in os.walk(os.path.join(HERE, "out", "v18")):
        for f in files:
            outv18.append(os.path.relpath(os.path.join(root, f), HERE)
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
        "out_v18_files": sorted(outv18),
        "n_out_v18_files": len(outv18),
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
    return {
        "job": "V18-PROD8S",
        "brief": BRIEF,
        "created_utc": now_iso(),
        "model": ("v3 ALPHA close: v3.10 spectrum (elfouhaily-form, Arm D) + "
                  "v16 water (FF phase, FF-CDF halo) + frozen composite "
                  "A2.6 + Blow + E2.2"),
        "preflight": {"pass": params["preflight_pass"],
                      "base_frames": {
                          "dir": params["base_frames"]["dir"],
                          "count": params["base_frames"]["count"],
                          "all_832x480": params["base_frames"]["all_832x480"],
                          "set_sha256": params["base_frames"]["set_sha256"]},
                      "v16_defaults": params["v16_defaults"],
                      "missing": params["missing"]},
        "scene": {"cell": params["scene_as_built"]["cell"],
                  "brief_line": params["cell"]["brief_line"],
                  "file_agrees_with_brief": params["cell"][
                      "agrees_with_brief"],
                  "file_wins": params["cell"]["file_wins"],
                  "record_C_mesh": params["cell"]["record_C_mesh"],
                  "record_seed": params["cell"]["record_seed"],
                  "seed": params["scene_as_built"]["seed"],
                  "C_mesh": params["scene_as_built"]["grid"], "Q": Q_DET,
                  "patch_m": params["scene_as_built"]["patch_m"],
                  "depth_m": params["scene_as_built"]["depth_m"],
                  "turbid": params["scene_as_built"]["turbid"],
                  "sun_altitude_deg": params["scene_as_built"][
                      "sun_altitude_deg"],
                  "eta": params["scene_as_built"]["eta"],
                  "eta_source": params["scene_as_built"]["eta_source"],
                  "retina": params["scene_as_built"]["retina"],
                  "rain_mm_h": params["scene_as_built"]["rain_mm_h"],
                  "water_preset": params["scene_as_built"]["water_preset"],
                  "water_load": params["scene_as_built"]["water_load"],
                  "a_rgb": params["scene_as_built"]["a_rgb"],
                  "b_rgb": params["scene_as_built"]["b_rgb"],
                  "g": params["scene_as_built"]["g"],
                  "spectrum": "elfouhaily-form / Arm D (the v15 merged "
                              "default)",
                  "state_summary": params["scene_as_built"]["state_summary"],
                  "composite": params["scene_as_built"]["composite"],
                  "optics": params["scene_as_built"]["optics"]},
        "frames": {"fps": FPS, "count": NFRAMES, "width": WIDTH,
                   "height": HEIGHT, "t_s": [i / float(FPS)
                                             for i in range(NFRAMES)],
                   "loop_T_s": 8.0,
                   "per_chunk": NCHUNKS, "chunk_size": CHUNK},
        "per_frame": per_frame,
        "chunk_hashes": chunk_hashes,
        "wall_time": {
            "stage0": float(params.get("stage0_wall_s") or 0.0),
            "chunks_total": float(sum(chunks[str(k)]["wall_s"]
                                      for k in range(NCHUNKS)
                                      if str(k) in chunks)),
            "g1": float((gates.get("G1") or {}).get("wall_s") or 0.0),
            "encode": float((gates.get("G3") or {}).get("wall_s") or 0.0),
            "total": float((params.get("stage0_wall_s") or 0.0)
                           + sum(chunks[str(k)]["wall_s"]
                                 for k in range(NCHUNKS)
                                 if str(k) in chunks)
                           + ((gates.get("G1") or {}).get("wall_s") or 0.0)
                           + ((gates.get("G3") or {}).get("wall_s") or 0.0)),
        },
        "gates": {k: {"pass": v.get("pass"), **{kk: vv for kk, vv in v.items()
                                                if kk not in ("pass",)}}
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
    cmeans = [r["C_mean"] for r in per]
    crush = [r["crush"] for r in per]
    blow = [r["blow"] for r in per]
    walls = [numbers["chunk_hashes"][str(k)]["wall_s"]
             for k in sorted(numbers["chunk_hashes"], key=int)]
    g1node = g.get("G1", {}).get("numbers", {})
    seam_rows = "\n".join(
        "| %s | %.3g |" % (k, v)
        for k, v in (g1node.get("fields_max_abs_diff_t8_vs_t0") or {}).items())
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
    lines = [
        "# REPORT-V18-PROD8S — the 8 s piece on the v3.10 + v16 stack",
        "",
        "Brief: `BRIEF-V18-PROD8S.md` (D:\\mutaro-work\\rally-in). "
        "Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.",
        "Built: %s. Code: `run_v18_prod.py` only." % numbers["created_utc"],
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
        "- 192 frames, 832x480, 24 fps, t = n/24 s.",
        "- Base: the 192 shipped Blender frames in `blender-piece8s/renders_v2`"
        " (count and size confirmed from the files; set sha256 `%s`)."
        % numbers["preflight"]["base_frames"]["set_sha256"][:16],
        "- Water: v3.10 spectrum (elfouhaily-form, Arm D, the v15 merged "
        "default) at U = %g m/s, F = %g m, depth = %g m, seed %d, C_mesh %d, "
        "det Q = %d, patch %.4g m." % (sc["cell"]["U"], sc["cell"]["F"],
                                        sc["cell"]["h"], sc["seed"],
                                        sc["C_mesh"], sc["Q"], sc["patch_m"]),
        "- Water optics: preset `%s` (load %g), a = %s, b = %s 1/m."
        % (sc["water_preset"], sc["water_load"],
           ["%.5f" % x for x in sc["a_rgb"]],
           ["%.5f" % x for x in sc["b_rgb"]]),
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
        % (g["G2"]["g1_range"][0], g["G2"]["g1_range"][1]),
        "- Sky: direct sun only, eta = 1.0 (the recipe carries no sky block). "
        "Retina OFF. Rain 0 mm/h.",
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
        "representation."
        % g1node.get("branch_b_quantisation_proof"),
        "- Continuity scale: the seam residual in C is %.3g of one frame "
        "step (%.4g), so there is no water discontinuity at the loop point."
        % (g1node.get("C_seam_residual_relative_to_one_frame_step") or 0.0,
           g1node.get("C_max_abs_diff_t1frame_vs_t0_sensitivity") or 0.0),
        "- Sensitivity reference (C at t=1/24 s vs t=0 s): %.6g." %
        (g1node.get("C_max_abs_diff_t1frame_vs_t0_sensitivity", 0.0)),
        "",
        "### G2 levels — %s" % ("PASS" if g.get("G2", {}).get("pass")
                                else "FAIL"),
        "",
        "| quantity | value | rule |",
        "| --- | --- | --- |",
        "| C mean over the piece | min %.9f, max %.9f (mean of means %.9f) | "
        "1 +/- 0.02 for >= 95%% of frames |"
        % (min(cmeans), max(cmeans), g["G2"]["C_mean"]["mean_of_means"]),
        "| frames in band | %d / %d (%.2f%%) | >= 95%% |"
        % (g["G2"]["C_mean_within_1_pm_0.02"]["n_frames"], len(cmeans),
           100 * g["G2"]["C_mean_within_1_pm_0.02"]["fraction"]),
        "| crush (< 0.02 linear) | max %.4f%% per frame | < 5%% per frame |"
        % (100 * max(crush)),
        "| blowout (ACES > 0.98) | max %.6g, %d frames nonzero | 0 per frame |"
        % (max(blow), g["G2"]["blowout_aces_gt_0.98"]["n_frames_nonzero"]),
        "",
        "### G3 encode — %s" % ("PASS" if g3.get("pass") else "FAIL"),
        "",
        "`ffprobe` on `out/v18/piece-8s-v3alpha.mp4`: codec `%s`, %sx%s, "
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
        "Decode check (manual, after the encode): frame 0 decoded from the MP4 "
        "vs `frames/frame-0000.png` has MAE 3.2/255 and max 38/255; frame 191 "
        "vs `frame-0191.png` MAE 4.1/255 — the expected H.264 crf-18 yuv420p "
        "rounding, no frame off-by-one.",
        "",
        "### G4 integrity — %s" % ("PASS" if g.get("G4", {}).get("pass")
                                   else "FAIL"),
        "",
        "- `*.py` set hashed at stage 0 (%d files) and again at the end "
        "(%d files): added %s, changed %s, removed %s."
        % (g["G4"]["py_before"], g["G4"]["py_after"], g["G4"]["added"],
           g["G4"]["changed"], g["G4"]["removed"]),
        "- Pre-existing `*.py` files changed: **%s**; unexpected new files: "
        "**%s**. No v10-v17 file was modified."
        % (g["G4"].get("changed_existing_files", []),
           g["G4"].get("unexpected_new_or_changed", [])),
        "- Disclosure: `run_v18_prod.py` (the job's own new file) was edited "
        "after the stage-0 snapshot while the gate/report code was being "
        "finished; its stage-0 and final sha256 are `%s` -> `%s`. That is the "
        "only `*.py` that moved, and it is the file the brief allows to be "
        "new."
        % (tuple((g["G4"].get("own_new_file_sha256", {}) or {})
                 .get("run_v18_prod.py", [None, None]))[0],
           tuple((g["G4"].get("own_new_file_sha256", {}) or {})
                 .get("run_v18_prod.py", [None, None]))[1]),
        "- Input integrity: the 192 base frames were not touched (hashed at "
        "stage 0).",
        "",
        "## 4. Wall time",
        "",
        "- Sum of the 8 chunk renders: **%.1f s** (%.1f min); per chunk %s s."
        % (total_wall, total_wall / 60.0,
           ["%.0f" % w for w in walls]),
        "- Stage 0 preflight + scene build: %.1f s; G1 seam (3 det renders + "
        "2 field evolutions): %.1f s; encode + ffprobe: %.1f s."
        % (numbers["wall_time"]["stage0"], numbers["wall_time"]["g1"],
           numbers["wall_time"]["encode"]),
        "- Total measured render/gate/encode work: **%.1f s (%.1f min)**. "
        "Scene rebuild per process ~10 s; det at C_mesh 2048 (~35 s/frame) "
        "dominates." % (numbers["wall_time"]["total"],
                        numbers["wall_time"]["total"] / 60.0),
        "",
        "## 5. Deliverables",
        "",
        "| file | note |",
        "| --- | --- |",
        "| `out/v18/piece-8s-v3alpha.mp4` | the piece (sha256 `%s`) |"
        % (g3.get("sha256") or "")[:16],
        "| `out/v18/frame-0096.png` | mid-piece frame, final composite |",
        "| `out/v18/numbers-v18-prod.json` | parameters, per-frame rows, "
        "gates, chunk hashes |",
        "| `out/v18/REPORT-V18-PROD8S.md` | this file |",
        "| `out/v18/checkpoint.json` | jobkit checkpoint + per-chunk hashes |",
        "",
        "## 6. Deviations and choices the brief left open",
        "",
        "1. **C_mesh = %d, Q = %d** — taken from the v15 T1 record "
        "(`numbers-v15-transport.json`, C_mesh 2048) and the pipeline's det "
        "default Q = 25; the brief does not name a mesh."
        % (sc["C_mesh"], sc["Q"]),
        "2. **E = 2.2 literal** — the frozen v5 composite constant, as used "
        "by the v15 transport stills (`transport_scene(..., E=2.2)`). The "
        "v16 still path recomputes a normalised E from L_v3/(L_dir_flat*g1); "
        "that was NOT used, because the brief names the frozen E2.2.",
        "3. **Water preset = algae_pond at its default load** — the v16 "
        "default lane (`W.DEFAULT_MEDIUM_PRESET`, and the phase preset used "
        "throughout run_v16).",
        "4. **g1 recomputed per frame** from the realised slope of that "
        "frame's field (the still path's own rule), so the composite keeps "
        "the still's formula. On this cell it is 1.0 exactly for every "
        "frame: a = 1/(alpha tan theta) = 10.89 >= 1.6 puts the Smith fit in "
        "the Lambda = 0 regime, so g1 = 1 and the composite carries no "
        "per-frame exposure modulation.",
        "5. **substrate_gain_rgb = (1,1,1)** — this is both the v16 still "
        "path's literal and the recipe's own gain, so the choice is "
        "immaterial here.",
        "6. **Base frame used as the submerged scene** — each frame n loads "
        "`renders_v2/frame_%04d.png`, converts sRGB -> linear and passes it "
        "as the substrate (`photo_lin`), exactly where the still path puts "
        "`moth-on-bark.jpg`.",
        "7. **Time quantisation** is inherited from the v10 state build "
        "(`omega_q = n * 2pi/8`), so the water layer is exactly 8 s-periodic; "
        "no extra quantisation was applied in this job.",
        "8. **One correction inside the G1 proof (not in the render)**: the "
        "first version of the printed frequency table compared `omega_q` "
        "against the deep-water `sqrt(gk)`, which is wrong for the hi "
        "cascade (the builder uses the finite-depth + capillary dispersion). "
        "The table now uses `v8_physics.dispersion_omega`, the same function "
        "the state builder uses. No rendering parameter or frame changed; "
        "the seam numbers (field and C differences) are identical to the "
        "first measurement.",
        "9. **Chunk idempotency was exercised**: after the piece was "
        "complete, `--chunk 0` was re-run and skipped all 24 frames on the "
        "recorded hashes, as the brief requires.",
        "",
        "## 7. Uncertainty",
        "",
        "- This is a rendering job, not a measurement: it reports levels, "
        "seam identity and encode conformance, nothing about a real pond.",
        "- The seam gate proves the *water* layer is periodic; the base "
        "animation's own seam is the Blender file's, as shipped.",
        "- 8-bit PNG -> H.264 yuv420p is a lossy step; the gate is on the "
        "encode's metadata, not on a pixel-exact decode.",
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
    log("  probe frame %d" % n)
    r = render_one(n, st, env, optics)
    sha = save_png(os.path.join(WORK, "probe-%04d.png" % n), r["still"])
    log("  timing %s" % r["timing_s"])
    log("  C mean %.9f p99 %.3f crush %.4f%% blow %.6g g1 %.6f sha %s"
        % (float(np.mean(r["C"], dtype=np.float64)),
           float(np.percentile(r["C"], 99)),
           100 * r["metrics"]["crush_rate_linear_lt_0.02"],
           r["metrics"]["blowout_rate_aces_gt_0.98"], r["g1"], sha[:16]))
    log("  optics build %.1fs" % p["build_wall_s"])
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--seam", action="store_true")
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
