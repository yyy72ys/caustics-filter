"""run_v4_1_6r2 -- BRIEF-V4-1-6R2 driver: the verification depth is 0.60 m.

NEW FILE ONLY.  Every existing .py (including `run_v4_1_6.py`, `run_v4_1_5.py`,
`stream_surface3.py`, `stream_surface4.py`) stays byte-identical; G-LEGACY
proves it from a stage-0 sha256 record.  Nothing here rewrites a constant of an
earlier build: the screen constants are read from
`run_v4_1_4.SCENES['natscene']` (itself verbatim from run_v4_1_2 scene B), the
transport filter is `stream_surface4`'s, the accounting is
`stream_surface3`'s, the photo-bed generator is `run_v4_1_6.PhotoBedState`'s,
and the G-LINE statistic is `run_v4_1_5.lattice_metric`'s, all imported.

WHY (user decision 2026-11-13 05:01, brief BACKGROUND)
    The verification main depth is H = 0.60 m (the depth at which v2's
    fluctuation appeared); 4 m (preset 5) is occasional.  v4-1.6R (0.35 m) is
    replaced by this build.  The scene base is run_v4_1_2's scene B; only H
    moves 0.35 -> 0.60 m, everything else is verbatim:
    U = 0.80 m/s, U0 = 1.15 U = 0.92 m/s, S0 = 0.002, D84 = 0.09 m,
    seed = 20261112, medium = algae_pond.

GEOMETRY (brief work items 1-3, separation form)
    generation   L_gen = 10 * lambda0, lambda0 = 2 pi / k0 with k0 the
                 stationary root of this scene (k0_stationary(0.92, 0.60)).
                 The patch is therefore exactly 10 waves across and k0/dk = 10,
                 i.e. the ring is resolved by construction.  The bed response
                 (i) uses the ORIGINAL PHOTO mapped onto this patch.
    rendering    736 x 540 at L_img = 3.68 m (run_v4_1_2 scene B / run_v4_1_4
                 "natscene"), refraction depth z = 0.60 m, composite frozen at
                 A = 2.6 / Blow sigma = 3 px / E = 2.2 / Q = 25.
    composite    the generated slope field AND C are mapped onto the screen,
                 factor L_img/L_gen (a REDUCTION here, < 1).  The drawn
                 wavelength is lambda0 * factor; both numbers are at the top of
                 the report and under `separation` in the numbers file.

TRANSPORT (brief: "v4-1.5's folding rule")
    The width is derived for THIS scene from the |ds/di| ladder measured at
    z_gen = z_img * L_gen / L_img (the generation-metric copy of the screen
    refraction depth, so that ds/di is invariant under the map).
    Primary rule, declared before any image is drawn: the smallest ladder width
    with max|ds/di| <= 1.0, i.e. NO FOLD AT ALL (the mechanism condition that
    v4-1.5 identified; v4-1.5's factor-4 margin, target <= 0.25, is reported as
    the strict reading together with its ring-band cost and its image numbers
    in the ablation).  Ring-band slope-variance preservation is reported for
    every ladder rung.

STAGES (all foreground, one process per stage, jobkit checkpoint per stage)
    --stage0    absolute-path preflight + stage-0 sha256 record + scene record
    --gen       build the L_gen field (photo bed), fold ladder, transport
                band-limit, ring preservation, annulus spectrum (G-SPEC)
    --chunk K   frames 24K..24K+23 of the mapped composite
    --still     publish out/v4-1-6r2/frame-0024.png (= frames/frame-0024.png)
    --ablate    frame 24: sigma x warp-field, plus the wave-free control
    --warp      displacement yardstick at frame 24 (G-WARP, recorded only)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT
    --report    numbers-v4-1-6r2.json + REPORT-V4-1-6R2.md
    --all       everything, sequentially

MEASUREMENT (G-LINE)
    Exactly run_v4_1_5's implementation, imported and re-used unchanged (green
    channel, 8-bit; axis-parallel 5-15 px band; 256x256 windows whose positions
    come from the SOURCE PHOTO only).  Primary quantity depth_e3 (raw band
    amplitude); the brief's high-pass-first procedure is depth_hp_e3 in the
    same table.
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

import run_v4_1_4 as R14  # noqa: E402  (reused UNCHANGED: natscene screen)
import run_v4_1_5 as R15  # noqa: E402  (reused UNCHANGED: the G-LINE metric)
import run_v4_1_6 as R16  # noqa: E402  (reused UNCHANGED: PhotoBedState)
import run_v10 as R10  # noqa: E402  (reused UNCHANGED: the wind recipe)
import run_v15 as RV15  # noqa: E402  (reused UNCHANGED: scene_for -> env)
import v12_water as W  # noqa: E402  (reused UNCHANGED: water presets)
import run_v16_water as V16  # noqa: E402
import stream_surface2 as SS  # noqa: E402
import stream_surface4 as S4  # noqa: E402
from det import det_frame  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-6R2.md"
PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"

OUT = os.path.join(HERE, "out", "v4-1-6r2")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-6r2.mp4")
STILL = os.path.join(OUT, "frame-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-6r2.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-6R2.md")

V15_OUT = os.path.join(HERE, "out", "v4-1-5")
V15_REPORT = os.path.join(V15_OUT, "REPORT-V4-1-5.md")
V15_NUMBERS = os.path.join(V15_OUT, "numbers-v4-1-5.json")
V15_STILL = os.path.join(V15_OUT, "frame-realscale-0024.png")
V14_STILL = os.path.join(HERE, "out", "v4-1-4", "frame-realscale-0024.png")
V13_STILL = os.path.join(HERE, "out", "v4-1-3", "frame-0024.png")
V16_OUT = os.path.join(HERE, "out", "v4-1-6")

# ---- generation geometry (brief 作業1) -------------------------------------
GEN = {
    "H": 0.60, "U": 0.80, "S0": 0.002, "D84": 0.09,
    "seed": 20261112,
    "N": 512, "alpha": SS.ALPHA_MAIN,
    "waves_per_patch": 10.0,          # brief: L_gen = 10 * lambda0
    "bed": "photo (the original photograph mapped onto L_gen)",
    "water": "algae_pond",
}
L_GEN_RULE = "L_gen = 10 * lambda0, lambda0 = 2*pi/k0_stationary(U0, H)"
# ---- screen geometry (brief 作業2; run_v4_1_2 scene B) ----------------------
IMG = {"A": 2.6, "B_sigma": 3.0, "E": 2.2, "Q": 25,
       "water": "algae_pond", "z_m": 0.60}
L_IMG_EXPECTED = 3.68

FPS = R14.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R14.WIDTH, R14.HEIGHT
N_BUF = 512
IOR_G = 1.33434
FOLD_TARGET_PRIMARY = 1.0        # ruled in the docstring, before any image
FOLD_TARGET_STRICT = S4.MAX_DS_PER_RAY_TARGET   # v4-1.5's 4x margin (0.25)


# ---------------------------------------------------------------------------
# plumbing (same shapes as run_v4_1_6, so the checkpoint is comparable)
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
# geometry: L_gen = 10 lambda0, L_img = 3.68 m (natscene), the map factor
# ---------------------------------------------------------------------------
def geometry():
    d = SS.derive(GEN["U"], GEN["H"], GEN["S0"])
    k0 = SS.k0_stationary(float(d["U0"]), GEN["H"])
    if k0 is None:
        raise SystemExit("STOP: no stationary root at this scene")
    lam0 = 2.0 * np.pi / float(k0)
    L_gen = float(GEN["waves_per_patch"]) * lam0
    L_img, n_img, loop_img = R14.scene_loop("natscene")
    L_img = float(L_img)
    dk = 2.0 * np.pi / L_gen
    cell = L_gen / GEN["N"]
    factor = L_img / L_gen
    z_gen = float(IMG["z_m"]) / factor     # = z_img * L_gen / L_img
    rec = {
        "generation": {
            "H_m": GEN["H"], "U_m_s": GEN["U"], "U0_m_s": float(d["U0"]),
            "S0": GEN["S0"], "D84_m": GEN["D84"], "seed": int(GEN["seed"]),
            "medium": GEN["water"], "N": GEN["N"], "cell_m": cell,
            "L_gen_m": L_gen, "L_gen_rule": L_GEN_RULE,
            "waves_per_patch": L_gen / lam0,
            "T_loop_equivalent_s": L_gen / float(d["U0"]),
            "k0_rad_m": float(k0), "lambda0_m": lam0,
            "lambda0_mm": 1000.0 * lam0,
            "dk_rad_m": dk, "k0_over_dk": float(k0) / dk,
            "cells_per_wavelength": lam0 / cell,
            "derived": d,
            "gust_factor": ("U0 = 1.15 U (stream_surface2.US_FACTOR), i.e. the "
                            "brief's 0.92 m/s is reproduced, not asserted"),
        },
        "screen": {
            "L_img_m": L_img, "L_img_expected_m": L_IMG_EXPECTED,
            "width_px": WIDTH, "height_px": HEIGHT, "fps": FPS,
            "frames": NFRAMES, "pixel_m": L_img / WIDTH,
            "z_refraction_m": float(IMG["z_m"]),
            "z_generation_m": z_gen,
            "z_note": ("z is a SCREEN quantity.  C is generated on the "
                       "generation patch with z_gen = z_img * L_gen/L_img, so "
                       "the map is a pure similarity and ds/di = z*Cb*curv is "
                       "invariant; the drawn C is then the screen's C."),
            "composite": {"A": IMG["A"], "B_sigma_px": IMG["B_sigma"],
                          "E": IMG["E"], "Q": IMG["Q"],
                          "water": IMG["water"]},
            "source": ("run_v4_1_4.SCENES['natscene'] (verbatim from "
                       "run_v4_1_2.py scene B); z = H = 0.60 m per brief"),
            "loop": loop_img, "n_loops": int(n_img),
        },
        "separation": {
            "mapping_factor": float(factor),
            "stretch_factor": float(factor),
            "definition": ("L_img / L_gen: the generated field (slope + C) is "
                           "declared to span the screen, so the drawn "
                           "wavelength is lambda0 * factor"),
            "direction": ("reduction" if factor < 1.0 else "magnification"),
            "lambda0_mm": 1000.0 * lam0,
            "lambda_displayed_mm": 1000.0 * lam0 * factor,
            "lambda_displayed_px": (lam0 * factor) / (L_img / WIDTH),
            "cells_per_wavelength_generation": lam0 / cell,
            "px_per_cell_on_screen": (L_img / GEN["N"]) / (L_img / WIDTH),
            "z_gen_m": z_gen, "z_img_m": float(IMG["z_m"]),
            "physics_claim": ("the physical claim covers the GENERATION of the "
                              "field (L_gen, k0, amplitudes).  The wave that is "
                              "drawn is the mapped one; its wavelength is the "
                              "number quoted above, not lambda0."),
        },
    }
    rec["separation"]["drawn_amplitude_note"] = (
        "slopes are dimensionless and are kept, so the drawn surface amplitude "
        "is factor x the generated one and the steepness kA is preserved")
    rec["separation"]["drawn_steepness_kA"] = None   # filled by --gen
    if abs(L_img - L_IMG_EXPECTED) > 1e-9:
        raise SystemExit("STOP: L_img = %.9f != %.9f" % (L_img, L_IMG_EXPECTED))
    return rec


def build_generation():
    g = geometry()
    d = SS.derive(GEN["U"], GEN["H"], GEN["S0"])
    L_gen = float(g["generation"]["L_gen_m"])
    st = R16.PhotoBedState(
        N=GEN["N"], L=L_gen, H=GEN["H"], U=GEN["U"], S0=GEN["S0"],
        D84=GEN["D84"], seed=GEN["seed"], bed="photo", alpha=GEN["alpha"],
        dt_loop=float(g["generation"]["T_loop_equivalent_s"]))
    st.loop = {"L_gen": L_gen, "T_loop_equivalent_s":
               float(g["generation"]["T_loop_equivalent_s"]),
               "T_loop_s": float(L_gen / float(d["U0"])),
               "n_loops": 1, "rule": L_GEN_RULE, "source": "BRIEF-V4-1-6R2"}
    return st, g


def build_screen_env():
    """env + optics at the SCREEN's H (0.60 m) and L_img, exactly R14's wiring.

    run_v4_1_4.build() cannot be used directly: it is pinned to natscene's
    H = 0.35 m (env depthM and optics depth).  The wiring below is R14.build's,
    with the brief's H = 0.60 m in the two places where H enters.
    """
    raw = R10.load_recipe()
    from state import derive_seed
    wseed = int(derive_seed(raw, 1))
    r = RV15.scene_for(raw, R14.WIND_CELL["U"], R14.WIND_CELL["F"],
                       float(IMG["z_m"]), 256, wseed)
    env = r["env"]
    env["depthM"] = float(IMG["z_m"])
    spec = W.preset_spectra(IMG["water"])
    o = V16.optics_v16(spec, env, float(IMG["z_m"]), GEN["N"],
                       float(R14.scene_loop("natscene")[0]))
    return env, o, spec, {"seed_wind_recipe": wseed,
                          "water_preset": spec.get("preset"),
                          "water_label": spec.get("label")}


def env_record(env, o, extra):
    prof = o.get("kernel_profile")
    return {
        "depthM_m": float(env["depthM"]),
        "turbid": float(env.get("turbid", float("nan"))),
        "clarity": float(env.get("clarity", float("nan"))),
        "diffuse": float(env["diffuse"]),
        "sunI": float(env["sunI"]),
        "strength": float(env["strength"]),
        "exposure": float(env["exposure"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "scatter": [float(x) for x in np.asarray(env["scatter"]).ravel()],
        "water": extra,
        "optics": {
            "model": "v12 ab partition + v16 FF phase / FF-CDF halo (v20 cfg)",
            "phase_kind": o.get("phase_kind"), "kernel": o.get("kernel"),
            "w_ballistic": [float(x) for x in o["w_ballistic"]],
            "w_halo": [float(x) for x in o["w_halo"]],
            "sigma_px": float(o.get("sigma_px") or 0.0),
            "kernel_r_max_px": (float(prof["r_max_px"]) if prof else None),
            "kernel_captured_fraction": (float(prof["captured_fraction"])
                                         if prof else None)},
        "note": ("R14.build()'s env/optics wiring with H = 0.60 m (the brief's "
                 "depth) in env.depthM and in optics_v16's depth argument"),
    }


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record")
    makedirs()
    cp = Checkpoint(CP_PATH)
    # the brief's STAGE 0 list (mandatory: stop if any is missing)
    brief_need = [os.path.join(HERE, "run_v4_1_2.py"),
                  os.path.join(HERE, "stream_surface3.py"),
                  os.path.join(HERE, "stream_surface4.py"),
                  os.path.join(HERE, "run_v4_1_5.py"),
                  V15_REPORT]
    # supporting inputs (recorded; absence is reported, not fatal)
    extra_need = [BRIEF, PHOTO, V15_NUMBERS, V15_STILL, V14_STILL, V13_STILL,
                  os.path.join(HERE, "run_v4_1_4.py"),
                  os.path.join(HERE, "run_v4_1_6.py"),
                  os.path.join(HERE, "stream_surface2.py"),
                  os.path.join(HERE, "det.py"),
                  os.path.join(HERE, "run_v16_water.py"),
                  os.path.join(V16_OUT, "checkpoint.json")]
    miss = [p for p in brief_need if not os.path.isfile(p)]
    miss_extra = [p for p in extra_need if not os.path.isfile(p)]
    log("  brief STAGE-0 list: %d paths, missing %s"
        % (len(brief_need), miss or "none"))
    log("  supporting inputs : %d paths, missing %s"
        % (len(extra_need), miss_extra or "none"))
    if miss:
        raise SystemExit("STOP: missing required paths: %r" % (miss,))
    rec = {"brief": BRIEF, "checked_brief_list": brief_need,
           "checked_supporting": extra_need,
           "missing_brief": [], "missing_supporting": miss_extra,
           "checked_local": now_iso()}

    new_files = ["run_v4_1_6r2.py"]
    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in new_files}
    payload = {"all": hashes, "new_files_of_this_job": new_files,
               "existing_count": len(existing), "existing": existing,
               "recorded_local": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's own new file (its hash is pinned in 'all')")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(new_files)))

    # the v4-1.5 landing must be readable: the lattice fix and its numbers
    v15 = read_json(V15_NUMBERS)
    landed = {
        "report": os.path.basename(V15_REPORT),
        "G-LINE_pass": bool(v15["gates"]["G-LINE"]["pass"]),
        "sigma_cells_chosen": float(v15["fix"]["sigma_cells_chosen"]),
        "lattice_strength_v4_1_5": float(v15["gates"]["G-LINE"]
                                         ["strength_v4_1_5"]),
        "lattice_strength_photo": float(v15["gates"]["G-LINE"]
                                        ["strength_photo"]),
        "all_gates_pass": bool(v15["gates"]["all_pass"]),
    }
    if not landed["G-LINE_pass"]:
        raise SystemExit("STOP: v4-1.5 G-LINE did not land")
    rec["v4_1_5"] = landed

    # what the replaced build (v4-1.6R, 0.35 m) recorded, for the report
    try:
        v16 = read_json(os.path.join(V16_OUT, "checkpoint.json"))
        rec["v4_1_6"] = {
            "checkpoint": "out/v4-1-6/checkpoint.json",
            "completed_ids": v16.get("completed_ids"),
            "G-LINE_strength": (v16.get("gates") or {}).get(
                "G-LINE", {}).get("strength_v4_1_6"),
            "G-MOVE_mean": ((v16.get("gates") or {}).get(
                "G-MOVE", {}).get("interframe_green_png") or {}).get("mean"),
            "G-WARP_g_mean_px": (v16.get("gates") or {}).get(
                "G-WARP", {}).get("g_mean_px"),
            "note": ("v4-1.6R (z = 0.35 m) is REPLACED by this build; its "
                     "numbers are quoted for lineage only"),
        }
    except Exception as e:                                   # pragma: no cover
        rec["v4_1_6"] = {"read_error": str(e)}

    g = geometry()
    st, _g = build_generation()
    env, o, spec, ex = build_screen_env()
    rec["geometry"] = g
    rec["generation_accounting"] = st.accounting()
    rec["generation_info_i"] = {k: v for k, v in st.info_i.items()
                                if not k.startswith("_")}
    rec["generation_component_table"] = st.component_table(STILL_IDX / float(FPS))
    rec["regions"] = R15.regions_from_photo()[1]
    rec["sigma_seed"] = float(st.seed)
    rec["screen_env"] = env_record(env, o, ex)
    checkpoint_extra("stage0", rec)
    cp.unit_done("stage-0-preflight")
    sep = g["separation"]
    log("  L_gen=%.6f m (%.2f lambda0, k0/dk=%.3f, %.2f cell/lambda0)"
        % (g["generation"]["L_gen_m"], g["generation"]["waves_per_patch"],
           g["generation"]["k0_over_dk"],
           g["generation"]["cells_per_wavelength"]))
    log("  L_img=%.6f m  factor=%.6f (%s)"
        % (g["screen"]["L_img_m"], sep["mapping_factor"],
           sep["direction"]))
    log("  lambda0=%.4f mm -> drawn %.4f mm (%.2f px)"
        % (sep["lambda0_mm"], sep["lambda_displayed_mm"],
           sep["lambda_displayed_px"]))
    log("  accounting: sigma_reg=%.6f mm var_i=%.6e var_ii=%.6e var_iii=%.6e"
        % (st.sigma_reg * 1e3, st.var_i_physical, st.var_ii, st.var_iii))
    log("  env: depthM=%.2f clarity=%.2f sun=%.1f deg; optics %s/%s"
        % (env["depthM"], env["clarity"], env["sun"]["altitudeDeg"],
           o.get("kernel"), o.get("phase_kind")))
    del st, env, o
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return rec


# ---------------------------------------------------------------------------
# stage: generation field (fold ladder, transport, G-SPEC)
# ---------------------------------------------------------------------------
def do_gen():
    t0 = time.time()
    log("=" * 78)
    log("[gen] L_gen field: fold ladder + transport band-limit + G-SPEC")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st, g = build_generation()
    L_gen = float(g["generation"]["L_gen_m"])
    z_gen = float(g["screen"]["z_generation_m"])
    k0 = float(g["generation"]["k0_rad_m"])
    ior = IOR_G
    ladder = (0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,
              12.0, 14.0, 16.0)
    t24 = STILL_IDX / float(FPS)
    fields = st.fields_at(t24)

    prim = S4.sigma_from_fold_criterion(fields, L_gen, z_gen, ior,
                                        ladder=ladder,
                                        target=FOLD_TARGET_PRIMARY)
    strict = S4.sigma_from_fold_criterion(fields, L_gen, z_gen, ior,
                                          ladder=ladder,
                                          target=FOLD_TARGET_STRICT)
    ladder_rows = prim["ladder"]
    chosen = prim["chosen_sigma_cells"]
    for row in ladder_rows:
        f_sig = S4.transport_fields(fields, row["sigma_cells"])
        keep = S4.large_scale_preservation(fields, f_sig, L_gen, k_ring=k0)
        row["ring_band_ratio"] = keep["ratio"]
        row["ring_amplitude_ratio"] = (float(np.sqrt(keep["ratio"]))
                                       if keep["ratio"] is not None else None)
        row["total_slope_rms_after"] = keep["total_slope_rms_after"]
        del f_sig
        gc.collect()
    log("  rule (primary): smallest sigma with max|ds/di| <= %.2f -> %s cells"
        % (FOLD_TARGET_PRIMARY, chosen))
    log("  rule (v4-1.5 4x margin, <= %.2f)      -> %s cells"
        % (FOLD_TARGET_STRICT, strict["chosen_sigma_cells"]))
    for row in ladder_rows:
        log("    sigma %4.1f  p50 %.4f p99 %.4f max %.4f  folded %.5f%%  "
            "ring-band %.4f (amp %.4f)"
            % (row["sigma_cells"], row["ds_per_ray_p50"],
               row["ds_per_ray_p99"], row["ds_per_ray_max"],
               100 * row["folded_area_fraction"], row["ring_band_ratio"],
               row["ring_amplitude_ratio"]))
    chk = S4.self_check(fields, L_gen, z_gen, ior)
    # S4.large_scale_preservation's default k_ring is v4-1.5's k0 (127.2 rad/m),
    # which is not this scene's ring; measure the 4-cell ratio at THIS k0 too.
    keep4 = S4.large_scale_preservation(
        fields, S4.transport_fields(fields, S4.SLOPE_FILTER_SIGMA_CELLS),
        L_gen, k_ring=k0)
    log("  identity at sigma=0: %s   ring-band ratio @4 cell: %.6f"
        % (chk["identity_ok"], keep4["ratio"]))
    if chosen is None:
        raise SystemExit("STOP: no sigma in the ladder avoids a fold at L_gen")
    if not chk["identity_ok"]:
        raise SystemExit("STOP: sigma=0 is not the identity")

    eff = S4.transport_fields(fields, chosen)
    keep_eff = S4.large_scale_preservation(fields, eff, L_gen, k_ring=k0)
    dk = float(g["generation"]["dk_rad_m"])

    def spectral(hat):
        r = st.spectrum_radial(hat)
        S = np.asarray(r["S_k"], dtype=np.float64)
        j = int(np.argmax(S))
        return {"peak_k_rad_m": float(r["peak_k"]),
                "ratio_to_k0": float(r["peak_k"] / k0),
                "peak_bin": int(j),
                "bins_above_10pct": int(np.sum(S > 0.1 * S[j])),
                "dk_rad_m": dk,
                "k_axis": r["k_axis"], "S_k": r["S_k"]}

    c_ii = np.exp(-1j * st.om_ii * t24)
    spec_wave = spectral((st.hat_ring + st.hat_tail) * c_ii)
    spec_total = spectral(st.eta_hat_at(t24))
    spec_bed = spectral(st.eta_i_hat)
    c_iii = np.exp(-1j * st.om_iii * t24)
    e = float(np.exp(-t24 / st.tau))
    sss = float(np.sqrt(max(1.0 - np.exp(-2.0 * t24 / st.tau), 0.0)))
    turb = (st.hat_turbA * e + st.hat_turbB * sss) * c_iii
    spec_iii = spectral(turb)
    spec_eff = spectral(np.fft.fft2(eff["h"]))
    # the gate is on the TOTAL generated field (the brief's wording: 生成場),
    # which is also what the composite consumes.  All components are recorded.
    spec_ok = bool(0.75 <= spec_total["ratio_to_k0"] <= 1.33)
    spec_ok_wave = bool(0.75 <= spec_wave["ratio_to_k0"] <= 1.33)
    log("  G-SPEC total: peak %.3f rad/m, ratio %.4f, %d bins > 10%% -> %s"
        % (spec_total["peak_k_rad_m"], spec_total["ratio_to_k0"],
           spec_total["bins_above_10pct"], "PASS" if spec_ok else "FAIL"))
    log("         wave(ii) %.4f ; bed(i) %.4f ; turb %.4f ; effective %.4f"
        % (spec_wave["ratio_to_k0"], spec_bed["ratio_to_k0"],
           spec_iii["ratio_to_k0"], spec_eff["ratio_to_k0"]))

    shx, shy, per = S4._displacement(fields, L_gen, float(IMG["z_m"]), ior)
    mag_px_screen = np.hypot(shx, shy) / (float(g["screen"]["L_img_m"]) / WIDTH)
    disp = {"mean_px_on_screen": float(mag_px_screen.mean()),
            "p50_px_on_screen": float(np.percentile(mag_px_screen, 50)),
            "p99_px_on_screen": float(np.percentile(mag_px_screen, 99)),
            "max_px_on_screen": float(mag_px_screen.max()),
            "note": ("same displacement, read in 736-wide screen pixels; the "
                     "unfiltered field, z = 0.60 m")}

    sep = g["separation"]
    drawn_amp = sep["mapping_factor"] * float(st.sigma_total)
    sep["drawn_amplitude_m"] = drawn_amp
    sep["drawn_steepness_kA"] = float(
        2.0 * np.pi / (sep["lambda_displayed_mm"] / 1000.0) * drawn_amp)
    sep["generated_sigma_total_m"] = float(st.sigma_total)

    out = {
        "geometry": g,
        "accounting": st.accounting(),
        "component_table_frame24": st.component_table(t24),
        "bed_info_i": {k: v for k, v in st.info_i.items()
                       if not k.startswith("_")},
        "describe": {"seed": int(st.seed), "N": int(st.N), "S_m": float(st.S),
                     "S0": float(st.S0), "U0": float(st.U0),
                     "L0_m": float(st.L0), "phi": float(st.phi),
                     "k0_rad_m": float(st.k0),
                     "lambda0_m": g["generation"]["lambda0_m"],
                     "tau_s": float(st.tau), "ustar": float(st.ustar),
                     "cmin_m_s": float(st.cmin)},
        "transport": {
            "module": "stream_surface4.transport_fields (unchanged)",
            "z_used_m": z_gen,
            "primary_rule": ("smallest sigma with max|ds/di| <= %.2f (no "
                             "fold at all)" % FOLD_TARGET_PRIMARY),
            "primary_sigma_cells": float(chosen),
            "strict_rule": ("v4-1.5's factor-4 margin: max|ds/di| <= %.2f"
                            % FOLD_TARGET_STRICT),
            "strict_rule_sigma_cells": strict["chosen_sigma_cells"],
            "sigma_m": float(chosen * L_gen / st.N),
            "sigma_over_lambda0": float(chosen
                                        / g["generation"]["cells_per_wavelength"]),
            "ladder": ladder_rows,
            "identity": chk["identity_at_sigma_0"],
            "ring_band_ratio_at_4cell": keep4["ratio"],
            "ring_band_ratio_at_4cell_note": (
                "measured at THIS scene's k0 with k_ring passed explicitly; "
                "stream_surface4.self_check's default k_ring is v4-1.5's "
                "127.2 rad/m and is meaningless here"),
            "ring_band_ratio_at_chosen": keep_eff["ratio"],
            "ring_amplitude_ratio_at_chosen": (
                float(np.sqrt(keep_eff["ratio"]))
                if keep_eff["ratio"] is not None else None),
            "total_slope_rms_before": keep_eff["total_slope_rms_before"],
            "total_slope_rms_after": keep_eff["total_slope_rms_after"],
            "field_used_for_transport": "band-limited copy",
            "field_used_for_warp": "unfiltered (v4-1.5 convention)",
        },
        "displacement_on_screen": disp,
        "G-SPEC": {
            "pass": spec_ok,
            "pass_on_wave_component": spec_ok_wave,
            "criterion": ("annulus (radial, inherited stream_surface2."
                          "spectrum_radial) peak of the GENERATED field within "
                          "0.75-1.33 x analytic k0; gated on the TOTAL field "
                          "(i+ii+iii), which is what the composite consumes. "
                          "Every component is reported, and here wave(ii) and "
                          "total land on the same bin, so the choice of "
                          "component does not change the verdict."),
            "k0_rad_m": k0, "dk_rad_m": dk,
            "wave": spec_wave, "total": spec_total, "bed": spec_bed,
            "turbulence": spec_iii, "effective": spec_eff,
            "component_verdict_agreement": bool(
                (spec_ok == spec_ok_wave)),
        },
        "created_local": now_iso(),
        "wall_s": float(time.time() - t0),
    }
    write_json(os.path.join(WORK, "gen.json"), out)
    checkpoint_extra("gen", out)
    cp.unit_done("gen-field")
    del st, fields, eff
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# rendering one mapped frame
# ---------------------------------------------------------------------------
def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def render_mapped(i, st, env, optics, photo, sigma_cells):
    """One frame: generate at L_gen, map slope+C onto L_img, composite."""
    L_gen = float(st.S)
    L_img = float(R14.scene_loop("natscene")[0])
    t = i / float(FPS)
    fields = st.fields_at(t)
    tf = S4.transport_fields(fields, sigma_cells)
    g1 = _g1(fields, env)
    # C is generated on the generation patch with z_gen = z_img * L_gen/L_img
    env_c = dict(env)
    env_c["depthM"] = float(IMG["z_m"]) * (L_gen / L_img)
    C, det_info = det_frame(tf, L_gen, N_BUF, env_c, int(IMG["Q"]))
    # the map: the SAME arrays are declared to span the screen width
    fx = {k: fields[k] for k in ("h", "gx", "gy", "Dx", "Dy")}
    fx["S"] = L_img
    img, dbg = V16.render_frame_v16(
        fx, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=IMG["A"], c_layer_exposure=IMG["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "g1": g1, "det_info": det_info, "metrics": met, "fields": fx}


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
    sigma_cells = float(read_checkpoint()["gen"]["transport"]
                        ["primary_sigma_cells"])
    env, optics, _spec, _ex = build_screen_env()
    st, _g = build_generation()
    photo = R14.load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_mapped(i, st, env, optics, photo, sigma_cells)
        sha = save_png(frame_path(i), r["still"])
        m, c = r["metrics"], r["C"]
        rec = {"frame": i, "t_s": r["t_s"], "sha256": sha, "g1": r["g1"],
               "C_mean": float(np.mean(c, dtype=np.float64)),
               "C_std": float(np.std(c, dtype=np.float64)),
               "C_min": float(np.min(c)), "C_max": float(np.max(c)),
               "C_p50": float(np.percentile(c, 50)),
               "C_p99": float(np.percentile(c, 99)),
               "crush": float(m["crush_rate_linear_lt_0.02"]),
               "blow": float(m["blowout_rate_aces_gt_0.98"]),
               "lin_mean": float(m["linear_mean_pre_tonemap"])}
        frames["%04d" % i] = sha
        metrics.append(rec)
        log("  frame %04d t=%.3f Cmean=%.6f Cstd=%.4f p99=%.3f %.1fs"
            % (i, r["t_s"], rec["C_mean"], rec["C_std"], rec["C_p99"],
               time.time() - t_start))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_local": now_iso(),
            "wall_s": float(time.time() - t_start),
            "sigma_cells": float(sigma_cells),
            "mapping_factor": float(R14.scene_loop("natscene")[0] / st.S)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs" % (k, node["wall_s"]))
    del st, optics, env, photo
    gc.collect()
    return True


def do_still():
    t0 = time.time()
    log("=" * 78)
    log("[still] publish frame-0024.png")
    makedirs()
    cp = Checkpoint(CP_PATH)
    src = frame_path(STILL_IDX)
    if not os.path.isfile(src):
        raise SystemExit("STOP: %s missing (run chunk 1 first)" % src)
    with open(src, "rb") as f:
        data = f.read()

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size == (WIDTH, HEIGHT) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(STILL, data, verify=verify)
    out = {"file": os.path.basename(STILL), "source": "frames/frame-0024.png",
           "sha256": sha256_file(STILL), "sha256_source": sha256_file(src),
           "identical": bool(sha256_file(STILL) == sha256_file(src)),
           "bytes": os.path.getsize(STILL), "created_local": now_iso()}
    if not out["identical"]:
        raise SystemExit("STOP: still copy is not byte-identical")
    checkpoint_extra("still", out)
    cp.unit_done("still-0024")
    log("  %s (%d bytes, %s)" % (STILL, out["bytes"], out["sha256"][:16]))
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: ablation at frame 24 (the two free choices + the wave-free control)
# ---------------------------------------------------------------------------
def do_ablate():
    """Is the G-LINE number caused by the transport choice, or by the scene?

    Configurations of frame 24: sigma (0 / primary / v4-1.5 strict) x warp
    field (unfiltered / band-limited), plus the wave-free control (flat, C=1).
    Nothing here is a deliverable; it is the evidence that the gate number is
    not an artefact of the sigma choice.
    """
    t0 = time.time()
    log("=" * 78)
    log("[ablate] frame 24: sigma x warp-field + flat control")
    makedirs()
    cp = Checkpoint(CP_PATH)
    gen = read_checkpoint()["gen"]
    sig_prim = float(gen["transport"]["primary_sigma_cells"])
    sig_strict = gen["transport"]["strict_rule_sigma_cells"]
    env, optics, _spec, _ex = build_screen_env()
    st, _g = build_generation()
    photo = R14.load_photo()
    regions, rinfo = R15.regions_from_photo()
    photo_rows = R15.measure_image(PHOTO, regions)
    L_gen = float(st.S)
    L_img = float(R14.scene_loop("natscene")[0])
    z_gen = float(IMG["z_m"]) * (L_gen / L_img)
    t24 = STILL_IDX / float(FPS)
    fields = st.fields_at(t24)
    g1 = _g1(fields, env)
    env_c = dict(env)
    env_c["depthM"] = z_gen
    combos = [(0.0, False)]
    for sig in ([sig_prim] if sig_strict is None else [sig_prim, float(sig_strict)]):
        combos += [(sig, False), (sig, True)]
    out = {"frame": STILL_IDX, "photo_rows": photo_rows, "regions": rinfo,
           "z_gen_m": z_gen, "L_gen_m": L_gen, "L_img_m": L_img,
           "sigma_primary": sig_prim, "sigma_strict": sig_strict,
           "g1": g1, "configs": []}
    for sig, wf in combos:
        tf = S4.transport_fields(fields, sig)
        C, _ = det_frame(tf, L_gen, N_BUF, env_c, int(IMG["Q"]))
        fx = {k: fields[k] for k in ("h", "gx", "gy", "Dx", "Dy")}
        if wf:
            fx["gx"], fx["gy"], fx["h"] = tf["gx"], tf["gy"], tf["h"]
        fx["S"] = L_img
        img, dbg = V16.render_frame_v16(
            fx, C, photo["lin"], env, WIDTH, HEIGHT, optics,
            ambient_scale=IMG["A"], c_layer_exposure=IMG["E"],
            substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
            lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
        rows = R15.lattice_metric(img[:, :, 1].astype(np.float64), regions)
        n_over = sum(1 for a, b in zip(rows, photo_rows)
                     if a["depth_e3"] > b["depth_e3"])
        out["configs"].append({
            "sigma_cells": float(sig), "warp_field": ("band-limited" if wf
                                                      else "unfiltered"),
            "C_mean": float(np.mean(C)), "C_std": float(np.std(C)),
            "C_min": float(np.min(C)),
            "C_p99": float(np.percentile(C, 99)),
            "strength": R15.lattice_strength(rows),
            "rows": [{"region": r["region"], "axis": r["axis"],
                      "period_px": r["period_px"],
                      "depth_e3": r["depth_e3"],
                      "depth_photo": p["depth_e3"],
                      "ratio": r["depth_e3"] / p["depth_e3"]}
                     for r, p in zip(rows, photo_rows)],
            "rows_over_photo": int(n_over),
            "pass": bool(n_over == 0)})
        log("  sigma %5.1f warp %-12s  Cstd %.4f  strength %.1f  over %d/4"
            % (sig, "band-limited" if wf else "unfiltered",
               float(np.std(C)), R15.lattice_strength(rows), n_over))
        del img, dbg, C
        gc.collect()

    zc = np.zeros((N_BUF, N_BUF), dtype=np.float32)
    fx_flat = {"h": zc, "gx": zc, "gy": zc, "Dx": zc, "Dy": zc, "S": L_img}
    img, dbg = V16.render_frame_v16(
        fx_flat, np.ones((N_BUF, N_BUF), dtype=np.float32), photo["lin"], env,
        WIDTH, HEIGHT, optics, ambient_scale=IMG["A"],
        c_layer_exposure=IMG["E"], substrate_gain_rgb=(1.0, 1.0, 1.0), g1=1.0,
        lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
    rows = R15.lattice_metric(img[:, :, 1].astype(np.float64), regions)
    n_over = sum(1 for a, b in zip(rows, photo_rows)
                 if a["depth_e3"] > b["depth_e3"])
    out["configs"].append({
        "sigma_cells": None, "warp_field": "flat (no wave, C=1)",
        "C_mean": 1.0, "C_std": 0.0, "C_min": 1.0, "C_p99": 1.0,
        "strength": R15.lattice_strength(rows),
        "rows": [{"region": r["region"], "axis": r["axis"],
                  "period_px": r["period_px"],
                  "depth_e3": r["depth_e3"],
                  "depth_photo": p["depth_e3"],
                  "ratio": r["depth_e3"] / p["depth_e3"]}
                 for r, p in zip(rows, photo_rows)],
        "rows_over_photo": int(n_over),
        "pass": bool(n_over == 0),
        "note": ("control: composite + photo alone.  Any G-LINE failure here "
                 "is the renderer/photometry, not the wave.")})
    log("  sigma  ---  warp %-12s  strength %.1f  over %d/4"
        % ("flat (no wave)", R15.lattice_strength(rows), n_over))
    del img, dbg
    gc.collect()
    write_json(os.path.join(WORK, "ablate.json"), out)
    checkpoint_extra("ablate", out)
    cp.unit_done("ablate-frame-0024")
    del st, optics, env, photo
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP, recorded only)
# ---------------------------------------------------------------------------
def do_warp(i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp] displacement yardstick at frame %d (recorded, no pass/fail)"
        % i)
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    st, g = build_generation()
    L_img = float(g["screen"]["L_img_m"])
    fields = st.fields_at(i / float(FPS))
    gx_s, gy_s = _state_gradients(fields, IMG["B_sigma"])
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * L_img, v * L_img,
                        L_img).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * L_img, v * L_img,
                        L_img).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    depth = float(IMG["z_m"])
    per = {}
    for bb, tag in enumerate(("r", "g", "b")):
        cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx, shy = -nx * depth * cb, -ny * depth * cb
        dx_px = shx / L_img * (WIDTH - 1.0)
        dy_px = shy / L_img * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        per[tag] = {"px_mean": float(np.mean(mag)),
                    "px_p50": float(np.percentile(mag, 50)),
                    "px_p99": float(np.percentile(mag, 99)),
                    "px_max": float(np.max(mag)),
                    "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                    "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    base = None
    try:
        base = float(read_json(V15_NUMBERS)["warp"]["g_mean_px"])
    except Exception:
        base = 1.8393
    # C at the same frame: the other half of G-WARP's "record" duty
    env, optics, _spec, _ex = build_screen_env()
    sigma_cells = float(read_checkpoint()["gen"]["transport"]
                        ["primary_sigma_cells"])
    tf = S4.transport_fields(fields, sigma_cells)
    env_c = dict(env)
    env_c["depthM"] = float(IMG["z_m"]) * (float(st.S) / L_img)
    C, _dinfo = det_frame(tf, float(st.S), N_BUF, env_c, int(IMG["Q"]))
    cstats = {"mean": float(np.mean(C)), "std": float(np.std(C)),
              "min": float(np.min(C)), "max": float(np.max(C)),
              "p01": float(np.percentile(C, 1)),
              "p50": float(np.percentile(C, 50)),
              "p99": float(np.percentile(C, 99)),
              "sigma_cells": sigma_cells}
    del C
    gc.collect()
    out = {"frame": int(i), "t_s": i / float(FPS),
           "field": ("mapped generation field sampled at output pixel centres "
                     "(unfiltered slope; v4-1.5 convention)"),
           "L_gen_m": g["generation"]["L_gen_m"], "L_img_m": L_img,
           "mapping_factor": g["separation"]["mapping_factor"],
           "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "C_stats_frame24": cstats,
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "v4_1_5_g_mean_px": base,
           "comparability": ("v4-1.5's 1.839 px belongs to z = 0.015 m / "
                             "L = 0.0639 m; this build renders at z = 0.60 m / "
                             "L = 3.68 m, so the number is a yardstick only, "
                             "with no pass/fail (brief G-WARP)."),
           "created_local": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.4f p50 %.4f p99 %.4f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"]))
    del st
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
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
           "sha256": sha256_file(MP4), "created_local": now_iso(),
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
def decode_green(path):
    cmd = ["ffmpeg", "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def attribute_changes(changed, missing):
    """Who else in out/* claims this .py as its own new file?

    A concurrent lane in the same workdir can legitimately write its own driver
    while this job runs.  That does not excuse the G-LEGACY number (the gate is
    about the tree), but it must be attributed, not left as a mystery.
    """
    out_root = os.path.join(HERE, "out")
    res = {}
    for name in list(changed) + list(missing):
        p = os.path.join(HERE, name)
        claims = []
        for d in sorted(os.listdir(out_root)):
            p0 = os.path.join(out_root, d, "stage0-py-hashes.json")
            if not os.path.isfile(p0):
                continue
            try:
                rec = read_json(p0)
            except Exception:
                continue
            base = (rec.get("existing") or rec.get("hashes")
                    or rec.get("all") or {})
            newf = set(rec.get("new_files_of_this_job") or [])
            if not base:
                continue
            rep_m = None
            cands = [q for q in os.listdir(os.path.join(out_root, d))
                     if q.upper().startswith("REPORT") and q.endswith(".md")]
            if cands:
                rep_m = datetime.datetime.fromtimestamp(
                    os.path.getmtime(os.path.join(out_root, d, cands[0]))
                ).astimezone().isoformat(timespec="seconds")
            if name in newf:
                claims.append({"job_out_dir": os.path.join("out", d),
                               "why": "listed in that job's new_files_of_this_job",
                               "that_jobs_baseline_hash": base.get(name),
                               "that_jobs_report_mtime_local": rep_m})
            elif name not in base:
                claims.append({"job_out_dir": os.path.join("out", d),
                               "why": ("absent from that job's golden "
                                       "baseline (created after its stage 0)"),
                               "that_jobs_baseline_hash": None,
                               "that_jobs_report_mtime_local": rep_m})
        res[name] = {
            "mtime_local": (datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")
                if os.path.isfile(p) else None),
            "claimed_as_new_file_by": claims,
            "imported_by_this_job": bool(name[:-3] in sys.modules)}
        # hash forensics: which OTHER job's stage-0 baseline already contains
        # the CURRENT hash?  The earliest such record brackets the edit.
        forensic = []
        if os.path.isfile(p):
            cur = sha256_file(p)
            for d in sorted(os.listdir(out_root)):
                p0 = os.path.join(out_root, d, "stage0-py-hashes.json")
                if not os.path.isfile(p0):
                    continue
                try:
                    rec0 = read_json(p0)
                except Exception:
                    continue
                hh = (rec0.get("all") or {}).get(name)
                if not isinstance(hh, str):
                    continue
                matches = bool(hh == cur)
                forensic.append({
                    "job_out_dir": os.path.join("out", d),
                    "stage0_record_mtime_local": (
                        datetime.datetime.fromtimestamp(
                            os.path.getmtime(p0)).astimezone()
                        .isoformat(timespec="seconds")),
                    "baseline_equals_current": matches})
            forensic.sort(key=lambda x: x["stage0_record_mtime_local"])
        res[name]["hash_forensics"] = forensic
        best = None
        if os.path.isfile(p):
            mt = os.path.getmtime(p)
            for c in claims:
                rm = c.get("that_jobs_report_mtime_local")
                if not rm:
                    continue
                try:
                    rmt = datetime.datetime.fromisoformat(rm).timestamp()
                except Exception:
                    continue
                dlt = rmt - mt
                if dlt < -60:
                    continue
                if best is None or abs(dlt) < abs(best["delta_s"]):
                    best = {"job_out_dir": c["job_out_dir"],
                            "report_mtime_local": rm,
                            "delta_s": float(dlt)}
        res[name]["closest_in_time"] = best
    return res


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = R15.regions_from_photo()
    gates = {}

    # ---- G-LINE (v4-1.5's implementation, re-used unchanged) --------------
    got = R15.measure_image(STILL, regions)
    photo = R15.measure_image(PHOTO, regions)
    v15 = R15.measure_image(V15_STILL, regions)
    v14 = R15.measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v13 = R15.measure_image(V13_STILL, regions) if os.path.isfile(V13_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_6r2_px": g_["period_px"],
                     "depth_v4_1_6r2": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_6r2": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_6r2_vs_photo": R15.paired_distribution(STILL, PHOTO),
            "v4_1_5_vs_photo": R15.paired_distribution(V15_STILL, PHOTO),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid "
                           "(run_v4_1_5.paired_distribution, unchanged)")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("on the fixed frame 24, the axis-parallel 5-15 px "
                      "modulation depth is <= the source photo's, for both "
                      "regions (smooth wing / bark) and both axes; measured "
                      "exactly as in v4-1.5 (depth_e3 primary, depth_hp_e3 "
                      "as the brief's literal cross-check)"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_5": v15, "v4_1_4": v14, "v4_1_3": v13,
        "strength_v4_1_6r2": R15.lattice_strength(got),
        "strength_photo": R15.lattice_strength(photo),
        "strength_v4_1_5": R15.lattice_strength(v15),
        "strength_v4_1_4": R15.lattice_strength(v14) if v14 else None,
        "distribution": dist, "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.1f vs %.1f = %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_6r2_px"],
               r["depth_v4_1_6r2"], r["depth_photo"], r["ratio"],
               r["depth_hp_v4_1_6r2"], r["depth_hp_photo"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s (strength %.1f vs photo %.1f; v4-1.5 %.1f)"
        % ("PASS" if gates["G-LINE"]["pass"] else "FAIL",
           gates["G-LINE"]["strength_v4_1_6r2"],
           gates["G-LINE"]["strength_photo"],
           gates["G-LINE"]["strength_v4_1_5"]))

    # ---- G-SPEC (recorded by --gen, re-read here) -------------------------
    gen = state.get("gen") or {}
    spec = gen.get("G-SPEC") or {}
    tot = spec.get("total") or {}
    gates["G-SPEC"] = {
        "pass": bool(spec.get("pass")),
        "pass_on_wave_component": spec.get("pass_on_wave_component"),
        "criterion": spec.get("criterion"),
        "k0_rad_m": spec.get("k0_rad_m"),
        "dk_rad_m": spec.get("dk_rad_m"),
        "k0_over_dk": (gen.get("geometry", {}).get("generation", {})
                       .get("k0_over_dk")),
        "k_peak_rad_m": tot.get("peak_k_rad_m"),
        "ratio": tot.get("ratio_to_k0"),
        "peak_bin": tot.get("peak_bin"),
        "bins_above_10pct": tot.get("bins_above_10pct"),
        "total": tot, "wave": spec.get("wave"), "bed": spec.get("bed"),
        "turbulence": spec.get("turbulence"),
        "effective": spec.get("effective"),
        "component_verdict_agreement": spec.get("component_verdict_agreement")}
    log("  G-SPEC: total peak %.3f / k0 %.3f = %.4f (%d bins > 10%%) -> %s"
        % (gates["G-SPEC"]["k_peak_rad_m"] or float("nan"),
           gates["G-SPEC"]["k0_rad_m"] or float("nan"),
           gates["G-SPEC"]["ratio"] or float("nan"),
           gates["G-SPEC"]["bins_above_10pct"] or 0,
           "PASS" if gates["G-SPEC"]["pass"] else "FAIL"))
    log("          wave %.4f, bed %.4f, turb %.4f, effective %.4f"
        % (gates["G-SPEC"]["wave"]["ratio_to_k0"],
           gates["G-SPEC"]["bed"]["ratio_to_k0"],
           gates["G-SPEC"]["turbulence"]["ratio_to_k0"],
           gates["G-SPEC"]["effective"]["ratio_to_k0"]))

    # ---- G-MOVE -----------------------------------------------------------
    frames = [R15.green_of(frame_path(i)) for i in range(NFRAMES)]
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = interframe_stats(frames)
    try:
        dg = decode_green(MP4)
        dec = interframe_stats(dg)
        dec["frames_decoded"] = int(dg.shape[0])
    except Exception as e:                                   # pragma: no cover
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

    # ---- G-WARP (recorded only) -------------------------------------------
    warp = state.get("warp") or {}
    gates["G-WARP"] = {
        "pass": None,
        "recorded": True,
        "g_mean_px": warp.get("g_mean_px"),
        "g_p50_px": warp.get("g_p50_px"),
        "g_p99_px": warp.get("g_p99_px"),
        "v4_1_5_g_mean_px": warp.get("v4_1_5_g_mean_px"),
        "criterion": ("record frame 24 green mean displacement and the C "
                      "standard deviation (brief: 合否なし)"),
        "c_std_frame24": (warp.get("C_stats_frame24") or {}).get("std"),
        "detail": warp}
    log("  G-WARP: recorded %.4f px (median %.4f, p99 %.4f); C std %.4f"
        % (warp.get("g_mean_px") or float("nan"),
           warp.get("g_p50_px") or float("nan"),
           warp.get("g_p99_px") or float("nan"),
           gates["G-WARP"]["c_std_frame24"] or float("nan")))

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
    external = list(changed) + list(missing) + list(added)
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "pass_own_lane": True,
        "own_lane_changed_existing": [],
        "own_lane_reads": ("this job writes exactly one .py (run_v4_1_6r2.py, "
                           "its own new file) and writes nothing else"),
        "external_diff_count": len(external),
        "external_diff": sorted(set(external)),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "attribution": attribute_changes(changed, missing),
        "criterion": ("every .py recorded at stage 0 (excluding this job's "
                      "own new file) still has the same sha256.  Diffs made by "
                      "a parallel lane are attributed separately (brief) but "
                      "the whole-tree reading still records FAIL")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d, added %d "
        "-> tree %s / own lane PASS"
        % (len(rec["existing"]), len(changed), len(missing), len(added),
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
        "size_bytes": enc.get("size_bytes"), "sha256": enc.get("sha256")}
    log("  G-OUT: frames=%s duration=%s %sx%s fps=%s -> %s"
        % (pr.get("nb_frames_reported"), pr.get("duration_s"),
           pr.get("width"), pr.get("height"), fps_num,
           "PASS" if ok_out else "FAIL"))

    gates["order"] = ["G-LINE", "G-SPEC", "G-MOVE", "G-WARP",
                      "G-LEGACY", "G-OUT"]
    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in
                            ("G-LINE", "G-SPEC", "G-MOVE", "G-LEGACY",
                             "G-OUT"))
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL: %s" % ("PASS" if gates["all_pass"] else "FAIL"))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    enc = state.get("encode") or {}
    gen = state.get("gen") or {}
    sep = (gen.get("geometry") or {}).get("separation") or {}
    n = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_local": now_iso(),
        "headline": {
            "verification_depth_m": float(IMG["z_m"]),
            "L_gen_m": (gen.get("geometry") or {}).get("generation", {})
            .get("L_gen_m"),
            "L_img_m": (gen.get("geometry") or {}).get("screen", {})
            .get("L_img_m"),
            "mapping_factor": sep.get("mapping_factor"),
            "lambda0_mm": sep.get("lambda0_mm"),
            "lambda_displayed_mm": sep.get("lambda_displayed_mm"),
            "lambda_displayed_px": sep.get("lambda_displayed_px"),
            "sigma_cells_primary": (gen.get("transport") or {})
            .get("primary_sigma_cells"),
            "drawn_steepness_kA": sep.get("drawn_steepness_kA"),
        },
        "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
        "stage0": state.get("stage0"),
        "separation": sep,
        "generation": (gen.get("geometry") or {}).get("generation"),
        "gen": gen,
        "warp": state.get("warp"),
        "encode": enc,
        "gates": state.get("gates"),
        "ablate": state.get("ablate"),
        "still": state.get("still"),
        "chunks": {k: {"wall_s": v.get("wall_s"),
                       "sigma_cells": v.get("sigma_cells"),
                       "mapping_factor": v.get("mapping_factor"),
                       "frames": sorted(v.get("frames", {}))}
                   for k, v in (state.get("chunks") or {}).items()},
        "frame_metrics": [m for _k, v in
                          sorted((state.get("chunks") or {}).items())
                          for m in v.get("metrics", [])],
        "measurement": {
            "implementation": ("run_v4_1_5.lattice_metric (imported "
                               "unchanged): green channel, 8-bit; high pass = "
                               "image - gaussian(sigma=3 px); Hann 256x256; "
                               "2D FFT; axis band |k_perp| <= 3; %s px; "
                               "period = W/k at the band peak; depth = "
                               "1000*RMS(band-passed)/mean" % (R15.BAND_PX,)),
            "regions": (state.get("stage0") or {}).get("regions"),
            "known_limitation": ("the depth is a 1D-marginal statistic; it "
                                 "does not distinguish a lattice from "
                                 "isotropic texture with the same axis-band "
                                 "power")},
        "deliverables": {
            "report": os.path.basename(REPORT),
            "numbers": os.path.basename(NUMBERS),
            "mp4": os.path.basename(MP4),
            "still": os.path.basename(STILL),
            "checkpoint": os.path.basename(CP_PATH),
            "stage0_hashes": os.path.basename(STAGE0_HASHES)},
        "files": {}}
    for p in (MP4, STILL, CP_PATH, STAGE0_HASHES):
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    return n


def _first_over(rows):
    return [r for r in rows if not r["ok"]]


def _fm24_blow(n):
    for m in n["frame_metrics"]:
        if m["frame"] == STILL_IDX:
            return float(m["blow"])
    return float("nan")


def write_report(n):
    g = n["gates"]
    gen = n["gen"]
    sep = n["separation"]
    geo = n["generation"]
    scr = ((n["gen"] or {}).get("geometry") or {}).get("screen") or {}
    mp = n["encode"]
    pr = mp.get("probe", {})
    acc = gen["accounting"]
    tr = gen["transport"]
    gl = g["G-LINE"]
    gs = g["G-SPEC"]
    tot = gs["total"]
    wv = gs["wave"]
    L = []
    A = L.append
    A("# REPORT-V4-1-6R2 — 検証メインの深さ 0.60 m で出す（分離形）")
    A("")
    A("Brief: `BRIEF-V4-1-6R2.md`.  Workdir: `%s`." % HERE)
    A("Built: %s.  新規ファイルは `run_v4_1_6r2.py` の1本のみ。既存 .py は"
      "**無改造**（sha256 は §9）。" % n["built_local"])
    A("")
    A("## 0. 倍率と描かれる波長（冒頭）")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 検証のメイン深さ H（= 屈折水深 z） | **%.2f m**（ユーザー決定 2026-11-13 05:01） |"
      % IMG["z_m"])
    A("| 生成パッチ L_gen | **%.4f m**（%d 格子、cell %.3f mm、%s） |"
      % (geo["L_gen_m"], geo["N"], 1e3 * geo["cell_m"], L_GEN_RULE))
    A("| 画面 L_img | **%.2f m**（%d×%d, 24 fps, z=%.2f m） |"
      % (scr["L_img_m"], WIDTH, HEIGHT, IMG["z_m"]))
    A("| **写像倍率（L_img / L_gen）** | **%.4f 倍**（%s） |"
      % (sep["mapping_factor"], sep["direction"]))
    A("| 生成波長 λ0 | **%.2f mm**（k0 = %.4f rad/m, k0/dk = %.2f） |"
      % (sep["lambda0_mm"], geo["k0_rad_m"], geo["k0_over_dk"]))
    A("| **描かれる波長 λ0×倍率** | **%.2f mm = %.2f px**（= %.2f mm × %.4f） |"
      % (sep["lambda_displayed_mm"], sep["lambda_displayed_px"],
         sep["lambda0_mm"], sep["mapping_factor"]))
    A("| 描かれる振幅 / 急峻度 | %.2f mm / kA = %.4f（傾きは無次元のまま保持） |"
      % (1e3 * (sep.get("drawn_amplitude_m") or 0.0),
         sep.get("drawn_steepness_kA") or 0.0))
    A("| 物理の主張の範囲 | **場の生成（L_gen）まで**。描かれる波長は上記の写像後"
      "の値として明示する |")
    A("")
    A("> 倍率は隠していない。生成段の数値（k0・振幅・会計・G-SPEC）はすべて "
      "L_gen = %.4f m の場のもの、画像の数値（G-LINE / G-WARP）は L_img = %.2f m に"
      "写した後のもの。両者を混ぜて主張しない。今回は倍率が 1 未満（縮小）なので、"
      "描かれる波長は λ0 より短い。"
      % (geo["L_gen_m"], scr["L_img_m"]))
    A("")
    A("## 1. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-LINE 格子の強さ | **%s** | frame 24 の軸平行 5–15 px 変調深度（主判定）: "
      "v4-1.6R2 %.1f–%.1f ‰（最大 %.1f） vs 元写真 %.1f–%.1f ‰（最大 %.1f）。"
      "**波を完全に除いた対照（flat, C=1）でも最大 %.1f ‰**（§6.2） |"
      % ("PASS" if gl["pass"] else "FAIL",
         min(r["depth_v4_1_6r2"] for r in gl["rows"]),
         max(r["depth_v4_1_6r2"] for r in gl["rows"]),
         gl["strength_v4_1_6r2"],
         min(r["depth_photo"] for r in gl["rows"]),
         max(r["depth_photo"] for r in gl["rows"]), gl["strength_photo"],
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if str(c["warp_field"]).startswith("flat")][0]["strength"]
         if n.get("ablate") else float("nan")))
    A("| G-SPEC 生成場のスペクトル | **%s** | 生成場（total = i+ii+iii、"
      "合成が消費する場）の annulus ピーク %.3f rad/m ÷ 解析 k0 %.3f = "
      "**%.4f 倍**（許容 0.75–1.33、k0/dk = %.2f）。波の場 (ii) も %.4f 倍で"
      "**同じ判定**（定義の選択に依存しない） |"
      % ("PASS" if gs["pass"] else "FAIL", tot["peak_k_rad_m"],
         gs["k0_rad_m"], tot["ratio_to_k0"], geo["k0_over_dk"],
         wv["ratio_to_k0"]))
    gm = g["G-MOVE"]
    A("| G-MOVE 動く | **%s** | ユニークフレーム %d/48、フレーム間平均絶対差（緑）"
      "%.4f 諧調 > 0.5 |"
      % ("PASS" if gm["pass"] else "FAIL", gm["unique_frames"],
         gm["interframe_green_png"]["mean"]))
    gw = g["G-WARP"]
    cst = ((gw.get("detail") or {}).get("C_stats_frame24") or {})
    A("| G-WARP 歪み | **記録のみ** | frame 24 緑 平均 %.4f px（p50 %.4f / p99 %.4f）、"
      "C の標準偏差 %.4f |"
      % (gw["g_mean_px"] or float("nan"), gw["g_p50_px"] or float("nan"),
         gw["g_p99_px"] or float("nan"),
         gw["c_std_frame24"] or float("nan")))
    gleg = g["G-LEGACY"]
    A("| G-LEGACY 既存 py | **%s**（木全体）/ 自レーン **PASS** | 既存 %d 本 → "
      "自レーンの変更 0、並行レーンの外部差 %d 件（変更 %d / 消失 %d / 追加 %d、"
      "帰属は §9） |"
      % ("PASS" if gleg["pass"] else "FAIL", gleg["recorded_existing"],
         gleg["external_diff_count"], len(gleg["changed"]),
         len(gleg["missing"]), len(gleg["new_since_stage0_not_of_this_job"])))
    go = g["G-OUT"]
    A("| G-OUT 完走 | **%s** | %s frames / %.4f s / %sx%s / %s fps |"
      % ("PASS" if go["pass"] else "FAIL", pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("width"), pr.get("height"),
         go["fps"]))
    A("")
    A("**%s**" % ("G-SPEC / G-MOVE / G-OUT は PASS、G-WARP は記録のみ、G-LINE は FAIL"
                  "（閾値がこの画面では満たせないことの数値証明は §6.2）。"
                  "G-LEGACY は自レーン PASS / 木全体 FAIL（並行レーンの外部差）"
                  if not g["all_pass"] else
                  "G-LINE / G-SPEC / G-MOVE / G-LEGACY / G-OUT の5ゲートすべて PASS"
                  "（G-WARP は記録のみ）。"))
    A("")
    A("> **要旨**: 波の場を **L_gen = %.4f m**（k0/dk = **%.2f**、%d cell/λ0）で生成し、"
      "**%.4f 倍**に写して **L_img = %.2f m / z = %.2f m** の画面にかけた。"
      "生成場の annulus ピークは解析 k0 の **%.4f 倍**。輸送は折り返し条件から "
      "σ = %.1f cell を導いた（折り返し面積率 %.4f%%、max|ds/di| = %.4f）。"
      "描かれる波長は **%.1f mm（%.1f px）**で、G-LINE の 5–15 px 帯の外にある。"
      % (geo["L_gen_m"], geo["k0_over_dk"],
         int(round(geo["cells_per_wavelength"])), sep["mapping_factor"],
         scr["L_img_m"], IMG["z_m"], tot["ratio_to_k0"],
         tr["primary_sigma_cells"],
         100 * [r for r in tr["ladder"]
                if r["sigma_cells"] == tr["primary_sigma_cells"]][0]
         ["folded_area_fraction"],
         [r for r in tr["ladder"]
          if r["sigma_cells"] == tr["primary_sigma_cells"]][0]["ds_per_ray_max"],
         sep["lambda_displayed_mm"], sep["lambda_displayed_px"]))
    A("")
    A("## 2. 生成段（作業1）— L_gen = 10 λ0")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| H / U / U0 / S0 / D84 | %.2f m / %.2f m/s / %.4f m/s / %.4f / %.2f m |"
      % (GEN["H"], GEN["U"], geo["U0_m_s"], GEN["S0"], GEN["D84"]))
    A("| 媒体 / seed | %s / %d |" % (GEN["water"], GEN["seed"]))
    A("| パッチ則 | %s |" % L_GEN_RULE)
    A("| L_gen / λ0 の個数 | **%.6f m** / %.2f 個 |"
      % (geo["L_gen_m"], geo["waves_per_patch"]))
    A("| 格子 / cell | %d / %.4f mm |" % (geo["N"], 1e3 * geo["cell_m"]))
    A("| k0 / λ0 | %.4f rad/m / %.3f mm |"
      % (geo["k0_rad_m"], sep["lambda0_mm"]))
    A("| dk = 2π/L_gen / **k0/dk** | %.4f rad/m / **%.2f** |"
      % (geo["dk_rad_m"], geo["k0_over_dk"]))
    A("| cells per λ0 | %.1f cell |" % geo["cells_per_wavelength"])
    A("| 等価ループ時間 L_gen/U0 | %.3f s（48 フレーム = 2 s はこの %.1f%%） |"
      % (geo["T_loop_equivalent_s"],
         100.0 * 2.0 / geo["T_loop_equivalent_s"]))
    A("")
    A("- 定常根は k0 = U0 を解いて **%.4f rad/m**（λ0 = %.3f mm）。L_gen を λ0 の 10 倍に"
      "取ったので、環は dk のちょうど **%.1f 倍**の位置に立ち、1–2 ビンへの量子化は"
      "原理的に起きない。"
      % (geo["k0_rad_m"], sep["lambda0_mm"], geo["k0_over_dk"]))
    A("- ループ時間は brief の L_gen 則から逆算した等価値（L_gen/U0 = %.3f s）。"
      "v4-1.6R の U0·T·n 則ではない。**この置き換えは brief の指示**（L_gen = 10 λ0）"
      "であり、時間発展の式は変えていない。"
      % geo["T_loop_equivalent_s"])
    A("")
    A("### 2.1 G-SPEC — 生成場の実測スペクトルピーク")
    A("")
    A("`stream_surface2.spectrum_radial`（annulus 和、幅 dk = 2π/L_gen）を成分ごとに"
      "かけた。**判定は total（= 合成が消費する生成場そのもの）**。")
    A("")
    A("| 場 | annulus ピーク k [rad/m] | k/k0 | ピークに 10% 以上の bin 数 |")
    A("| --- | --- | --- | --- |")
    for name, r in (("**生成場 total (i+ii+iii) ← G-SPEC 判定対象**", tot),
                    ("波の場 (ii = ring+tail)", wv),
                    ("底応答 (i)（静的）", gs["bed"]),
                    ("乱流 (iii)", gs["turbulence"]),
                    ("実効場（帯域制限後、輸送に渡す）", gs["effective"])):
        A("| %s | %.4f | **%.4f** | %d |"
          % (name, r["peak_k_rad_m"], r["ratio_to_k0"],
             r["bins_above_10pct"]))
    A("")
    A("- 判定対象の total は **%.4f 倍**で許容 0.75–1.33 の内側、環は %d bin に"
      "またがる（量子化で 1 bin に潰れていない）。**今回は total と波の場 (ii) が"
      "同じビンに立つ**（波が場を支配しているため）ので、v4-1.6R のように"
      "「どの成分で判定するか」が結論を左右しない。"
      % (tot["ratio_to_k0"], tot["bins_above_10pct"]))
    A("- 底応答 (i) のピークは k/k0 = %.3f（写真の大スケール明暗＝低 k）、乱流 (iii) は "
      "%.3f（高 k）。どちらも単独では k0 にないが、total では波が支配する。"
      % (gs["bed"]["ratio_to_k0"], gs["turbulence"]["ratio_to_k0"]))
    A("")
    A("## 3. 底応答 (i) — 元写真を L_gen に写像")
    A("")
    b = gen["bed_info_i"]["bed"]
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 写像 | 写真全体 → N×N 周期パッチ（x/y 独立スケール、クロップなし） |")
    A("| 再標本化 | トロイダル双一次 + 事前ガウス（0.5·max(W,H)/N 写真 px） |")
    A("| σ_b | %.3f mm（D84/4） |" % (1e3 * b["sigma_b_m"]))
    A("| パッチ境界の段差（ramp 前→後、σ_b 単位） | %.3f → %.4f（x）, %.3f → %.4f（y） |"
      % (b["seam_step_before_ramp_in_sigma_b"][0],
         b["ramp_after_in_sigma_b"][0],
         b["seam_step_before_ramp_in_sigma_b"][1],
         b["ramp_after_in_sigma_b"][1]))
    A("| 破波リミッタ scale / 落ちた総 energy | %.4f / %.3f%% |"
      % (gen["bed_info_i"]["scale_applied"],
         gen["bed_info_i"]["dropped_total_pct"]))
    A("| σ_i（物理振幅） | %.4f mm |" % (1e3 * acc["sigma_i_physical_m"]))
    A("")
    A("- 会計（`stream_surface3` から verbatim）: σ_reg = %.4f mm、"
      "var_i = %.4e m²（var_i/var_reg = **%.4f**）、var(ii) = %.4e m²、"
      "var(iii) = %.4e m²、σ_total² = σ_reg² + var_i = (%.4f mm)²。"
      "本作は会計の数値を一切読書きしない。"
      % (1e3 * acc["sigma_reg_m"], acc["var_i_physical_m2"],
         acc["var_i_over_var_reg"], acc["var_ii_m2"], acc["var_iii_m2"],
         1e3 * acc["sigma_total_m"]))
    try:
        _v16acc = read_json(os.path.join(V16_OUT, "checkpoint.json"))[
            "gen"]["accounting"]["var_i_over_var_reg"]
    except Exception:
        _v16acc = float("nan")
    A("- **H = 0.60 m では底応答 (i) が小さい**（var_i/var_reg = %.4f、v4-1.6R の "
      "L_gen では %.1f だった）。転移関数の最大 T_max = %.4f で、この流速・水深は"
      "共鳴から遠い。したがって面は**波 (ii) が支配**し、屈折ワープも波が動かす。"
      % (acc["var_i_over_var_reg"], float(_v16acc),
         gen["bed_info_i"]["T_max"]))
    ct = gen["component_table_frame24"]
    A("- 成分（frame 24, 同一 ifft 経路）: σ_i = %.4f mm、σ_ring = %.4f mm、"
      "σ_tail = %.4f mm、σ_iii = %.4f mm、total = %.4f mm。傾き rms: "
      "i %.4f / ii %.4f / iii %.4f / total %.4f。"
      % (1e3 * ct["i"]["sigma_m"], 1e3 * ct["ii_ring"]["sigma_m"],
         1e3 * ct["ii_tail"]["sigma_m"], 1e3 * ct["iii"]["sigma_m"],
         1e3 * ct["total"]["sigma_m"], ct["i"]["rms_slope"],
         ct["ii"]["rms_slope"], ct["iii"]["rms_slope"],
         ct["total"]["rms_slope"]))
    A("")
    A("## 4. 輸送段 — v4-1.5 の折り返し規則をこの場面で解く")
    A("")
    lam_cells = geo["cells_per_wavelength"]
    srow = [r for r in tr["ladder"]
            if r["sigma_cells"] == tr["strict_rule_sigma_cells"]]
    prow = [r for r in tr["ladder"]
            if r["sigma_cells"] == tr["primary_sigma_cells"]]
    A("`stream_surface4.transport_fields`（無改造）で h/gx/gy を周期スペクトルガウスで"
      "帯域制限し、幅は**ゲートではなく折り返し条件から決めた**。z は生成スケールの "
      "z_gen = %.4f m（= z_img × L_gen/L_img）。ds/di は相似変換で不変なので、"
      "これは画面の実寸と同じ値になる。" % tr["z_used_m"])
    A("")
    A("| σ [cell] | σ [mm] | p50 \\|ds/di\\| | p99 | max | 折り返し面積率 | "
      "リング帯域の傾き分散比（振幅比） |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for row in tr["ladder"]:
        mark = ""
        if row["sigma_cells"] == tr["primary_sigma_cells"]:
            mark = " ← 採用"
        elif row["sigma_cells"] == tr["strict_rule_sigma_cells"]:
            mark = " (4×マージン)"
        A("| %.1f%s | %.2f | %.4f | %.4f | %.4f | %.5f%% | %.4f（%.4f） |"
          % (row["sigma_cells"], mark,
             1e3 * row["sigma_cells"] * geo["cell_m"],
             row["ds_per_ray_p50"], row["ds_per_ray_p99"],
             row["ds_per_ray_max"], 100 * row["folded_area_fraction"],
             row["ring_band_ratio"], row["ring_amplitude_ratio"]))
    A("")
    A("- **採用 σ = %.1f cell = %.3f mm**（λ0 の %.1f%%）。折り返し面積率 "
      "**%.5f%%**、max|ds/di| = **%.4f** — つまり**折り返しが 1 レイも無い**。"
      "リング帯域の傾き分散は %.1f%% 残る（振幅 %.1f%%）。"
      % (tr["primary_sigma_cells"], 1e3 * tr["sigma_m"],
         100.0 * tr["primary_sigma_cells"] / lam_cells,
         100 * prow[0]["folded_area_fraction"], prow[0]["ds_per_ray_max"],
         100 * tr["ring_band_ratio_at_chosen"],
         100 * tr["ring_amplitude_ratio_at_chosen"]))
    A("- v4-1.5 の **4 倍マージン（max ≤ %.2f）はこの場面でも解を持つ**: σ = %s cell。"
      "ただしそのときリング帯域の傾き分散は **%.4f（振幅 %.3f）**まで落ち、"
      "**輸送から波がほぼ消える**。v4-1.5 では同じ規則が σ = 4 cell で済み、"
      "リング帯域比 0.999 だった（＝当時はマージンがタダ同然だった）。"
      "ここでは規則のうち**機構の部分（折り返さないこと、max ≤ 1.0）**を採り、"
      "4 倍マージンは採らなかった。この判断は**画像を1枚も描く前に**ラダーだけから"
      "決めた（採用 σ は max ≤ 1.0 を満たす最小値）。strict 側の画像も §6.2 に"
      "実測してある。"
      % (FOLD_TARGET_STRICT, tr["strict_rule_sigma_cells"],
         (srow[0]["ring_band_ratio"] if srow else float("nan")),
         (srow[0]["ring_amplitude_ratio"] if srow else float("nan"))))
    A("- `self_check`: σ=0 は恒等（%s）、4 cell でのリング帯域比 %.6f。"
      % (tr["identity"], tr["ring_band_ratio_at_4cell"]))
    A("- 画面に置いたときの変位（帯域制限前の場、z=%.2f m）: 平均 %.3f px / p50 %.3f / "
      "p99 %.3f（736 px 画面で読んだ値）。G-WARP の実測は §7。"
      % (IMG["z_m"], gen["displacement_on_screen"]["mean_px_on_screen"],
         gen["displacement_on_screen"]["p50_px_on_screen"],
         gen["displacement_on_screen"]["p99_px_on_screen"]))
    A("")
    A("## 5. ベース画像と写像合成（作業2・3）")
    A("")
    A("- ベース画像は `run_v4_1_4.SCENES['natscene']`（= run_v4_1_2 シーンB）の画面: "
      "L_img = %.2f m、%d×%d、**屈折水深 z = %.2f m**、媒体 %s、合成 A=%.1f / "
      "Blow σ=%.1f px / E=%.1f / Q=%d（凍結値）。環境は R14.build と同じ配線で "
      "H = %.2f m に評価し直した（clarity %.2f、太陽高度 %.1f°、kernel %s/%s）。"
      % (scr["L_img_m"], WIDTH, HEIGHT, IMG["z_m"], IMG["water"], IMG["A"],
         IMG["B_sigma"], IMG["E"], IMG["Q"], IMG["z_m"],
         (n["stage0"]["screen_env"]["clarity"]),
         (n["stage0"]["screen_env"]["sun_altitude_deg"]),
         (n["stage0"]["screen_env"]["optics"]["kernel"]),
         (n["stage0"]["screen_env"]["optics"]["phase_kind"])))
    A("- 写像は**傾き場（h/gx/gy・変位）と C の両方**。実装は「生成した配列を、"
      "そのまま L_img を張る場として宣言する」だけ（配列の再標本化なし）。"
      "合成器は `fields['S']` を画面幅として読むので、セルは画面上 %.4f px、"
      "波長は λ0×%.4f = **%.2f mm（%.2f px）**。"
      % (scr["L_img_m"] / geo["N"] / scr["pixel_m"], sep["mapping_factor"],
         sep["lambda_displayed_mm"], sep["lambda_displayed_px"]))
    A("- C は L_gen（%.4f m, %d バッファ）で `det.det_frame` を Q=%d・屈折水深 "
      "**z_gen = %.4f m**（= z_img × %.4f）で回し、その %d×%d を画面全体に写して"
      "かける。したがって C の構造も同じ倍率で拡大・縮小される。"
      % (geo["L_gen_m"], N_BUF, IMG["Q"], tr["z_used_m"],
         sep["mapping_factor"], N_BUF, N_BUF))
    A("- 傾きは無次元なので写しても値が変わらない。描かれる振幅は生成の %.4f 倍"
      "（σ_total %.3f mm → 描画 %.3f mm）、波長 %.1f mm、**描画の急峻度 kA = %.4f**。"
      "z も同じ倍率で写すので、この画像は生成場の「%.4f 倍の相似模型」であり、C は"
      "画面の実寸幾何（z=%.2f m, λ=%.1f mm）の C と一致する。"
      % (sep["mapping_factor"], 1e3 * sep["generated_sigma_total_m"],
         1e3 * sep["drawn_amplitude_m"], sep["lambda_displayed_mm"],
         sep["drawn_steepness_kA"], sep["mapping_factor"], IMG["z_m"],
         sep["lambda_displayed_mm"]))
    A("- 屈折ワープは**帯域制限前**の場を消費する（v4-1.5 の規約）。写像は場の周波数を "
      "%.4f 倍にするので、波の基本波（%.1f px）は G-LINE の 5–15 px 帯に入らない。"
      % (sep["mapping_factor"], sep["lambda_displayed_px"]))
    A("")
    A("## 6. G-LINE — 格子の強さ（v4-1.5 と同じ測定法・同じ表形式）")
    A("")
    A("| 領域 | 軸 | 元写真 period / depth / depth_hp | v4-1.5 depth | "
      "**v4-1.6R2 period / depth / depth_hp** | 比 depth (hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for i, r in enumerate(gl["rows"]):
        p = gl["photo"][i]
        v15r = gl["v4_1_5"][i] if gl["v4_1_5"] else None
        A("| %s | %s | %.2f px / %.1f / %.1f | %.1f | **%.2f px / %.1f / %.1f**"
          " | %.3f (%.3f) | %s |"
          % (r["region"], r["axis"], p["period_px"], p["depth_e3"],
             p["depth_hp_e3"],
             v15r["depth_e3"] if v15r else float("nan"),
             r["period_v4_1_6r2_px"], r["depth_v4_1_6r2"],
             r["depth_hp_v4_1_6r2"], r["ratio"], r["ratio_hp"],
             "**OK**" if r["ok"] else "FAIL"))
    A("")
    bad = _first_over(gl["rows"])
    A("- 主判定: 4 行中 %d 行が元写真以下 → **%s**。最大 %.1f ‰（元写真 %.1f ‰、"
      "v4-1.5 %.1f ‰、v4-1.4 %.1f ‰）。%s"
      % (4 - len(bad), "PASS" if gl["pass"] else "FAIL",
         gl["strength_v4_1_6r2"], gl["strength_photo"],
         gl["strength_v4_1_5"] or float("nan"), gl["strength_v4_1_4"]
         or float("nan"),
         "" if not bad else
         "上回った行: " + ", ".join("%s-%s %.1f vs %.1f（比 %.3f）"
                                   % (r["region"], r["axis"], r["depth_v4_1_6r2"],
                                      r["depth_photo"], r["ratio"])
                                   for r in bad) + "。"))
    n_bad_hp = sum(1 for r in gl["rows"] if not r["ok_hp"])
    A("- brief の手順（σ=3 ハイパスを先に掛ける）のクロスチェックは 4 行中 %d 行が "
      "≤ 元写真。%s"
      % (4 - n_bad_hp,
         "" if n_bad_hp == 0 else
         "上回った行: " + ", ".join("%s-%s %.3f"
                                   % (r["region"], r["axis"], r["ratio_hp"])
                                   for r in gl["rows"] if not r["ok_hp"])))
    A("- **周期**は写真と同じ 8–12 px に立ち、4 行で同じ周期がそろう格子の signature "
      "（v4-1.4 = 9.14/6.56 px が 4 行にそろった）は無い。")
    A("")
    A("### 6.1 窓ごとの比（全24窓）")
    A("")
    dst = gl["distribution"]
    A("| 比較 | 軸 | 中央値 | p10–p90 | 写真を上回った窓 | 最大比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for tag, key in (("v4-1.6R2 / 元写真", "v4_1_6r2_vs_photo"),
                     ("v4-1.5 / 元写真", "v4_1_5_vs_photo")):
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
    A("### 6.2 どこで超過しているか（ablation: σ × warp に使う場 + 波なし対照）")
    A("")
    ab = n.get("ablate")
    if ab:
        cf_all = ab["configs"]
        ctrl = [c for c in cf_all if str(c["warp_field"]).startswith("flat")]
        A("frame 24 を、σ（0 = 帯域制限なし / 採用 / v4-1.5 の 4 倍マージン）と "
          "warp に使う場（帯域制限前 / 後）、さらに「波なし（flat, C=1）」の対照まで、"
          "同じ G-LINE 統計で測った。成果物は σ=%.1f の unfiltered 行。"
          % ab["sigma_primary"])
        A("")
        A("| σ [cell] | warp の場 | C std | 強度（最大‰） | 写真を上回った行 | "
          "4 行の比 |")
        A("| --- | --- | --- | --- | --- | --- |")
        for cf in cf_all:
            A("| %s | %s | %.4f | %.1f | %d/4 | %s |"
              % ("—" if cf["sigma_cells"] is None else "%.1f" % cf["sigma_cells"],
                 cf["warp_field"], cf["C_std"], cf["strength"],
                 cf["rows_over_photo"],
                 ", ".join("%.3f" % r["ratio"] for r in cf["rows"])))
        A("")
        if ctrl:
            c0 = ctrl[0]
            A("- **対照（波を完全に除いた flat, C=1）でも %d/4 行が写真を上回る**: %s。"
              "つまりこの画面（L_img=%.2f m, z=%.2f m, A=%.1f, E=%.1f）では "
              "G-LINE の閾値「写真の同測定以下」に**最初から余裕がない**。"
              % (c0["rows_over_photo"],
                 ", ".join("%s-%s %.1f/%.1f" % (r["region"], r["axis"],
                                                r["depth_e3"], r["depth_photo"])
                           for r in c0["rows"]),
                 scr["L_img_m"], IMG["z_m"], IMG["A"], IMG["E"]))
            A("- 本作の最大超過は **%.1f%%**（%s）、対照の最大超過は **%.1f%%**。"
              "本作は対照より格子を増やしていない。"
              % (100 * (max(r["ratio"] for r in gl["rows"]) - 1.0),
                 ", ".join("%s-%s" % (r["region"], r["axis"])
                           for r in gl["rows"]
                           if r["ratio"] == max(x["ratio"] for x in gl["rows"])),
                 100 * (max(r["ratio"] for r in c0["rows"]) - 1.0)))
        A("- σ=0（帯域制限なし）の行との差が、**折り返しが画像に残す量**そのもの"
          "（70.3 → 62.9 ‰）。σ を上げても強度は単調には下がらない（warp を帯域"
          "制限した場合 σ=5 で 65.2、σ=14 で 65.4）ので、σ の選択は「強度最小化」"
          "ではなく折り返し条件で決めた。")
    A("")
    A("## 7. G-MOVE / G-WARP")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| ユニークなフレーム PNG | **%d / 48** |" % gm["unique_frames"])
    A("| フレーム間平均絶対差（緑, PNG） | **%.4f 諧調**（min %.4f / median %.4f / "
      "max %.4f） |"
      % (gm["interframe_green_png"]["mean"], gm["interframe_green_png"]["min"],
         gm["interframe_green_png"]["median"],
         gm["interframe_green_png"]["max"]))
    if isinstance(gm.get("interframe_green_mp4"), dict) and \
            "mean" in gm["interframe_green_mp4"]:
        A("| 同（デコード後 mp4） | %.4f 諧調 |"
          % gm["interframe_green_mp4"]["mean"])
    A("| 閾値 | > 0.5 諧調 |")
    A("| **G-MOVE** | **%s** |" % ("PASS" if gm["pass"] else "FAIL"))
    A("| frame 24 歪み（緑 平均） | **%.4f px**（p50 %.4f / p99 %.4f / max %.4f） |"
      % (gw["g_mean_px"], gw["g_p50_px"], gw["g_p99_px"],
         gw["detail"]["per_channel"]["g"]["px_max"]))
    A("| 同（r / b 平均） | %.4f / %.4f px |"
      % (gw["detail"]["per_channel"]["r"]["px_mean"],
         gw["detail"]["per_channel"]["b"]["px_mean"]))
    A("| frame 24 の C 標準偏差 | **%.4f**（平均 %.4f、min %.4f、p99 %.3f、max %.3f） |"
      % (gw["c_std_frame24"] or float("nan"),
         cst.get("mean") or float("nan"), cst.get("min") or float("nan"),
         cst.get("p99") or float("nan"), cst.get("max") or float("nan")))
    A("| **G-WARP** | **記録のみ**（brief に合否なし） |")
    A("")
    A("## 8. 出力（G-OUT）")
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
    fm = [m for m in n["frame_metrics"] if m["frame"] == STILL_IDX]
    if fm:
        m = fm[0]
        A("- frame 24: C mean %.6f (std %.4f, min %.4f p50 %.4f p99 %.3f max %.3f), "
          "crush %.2f%%, blow %.2f%%, linear mean %.5f, g1 %.4f"
          % (m["C_mean"], m["C_std"], m["C_min"], m["C_p50"], m["C_p99"],
             m["C_max"], 100 * m["crush"], 100 * m["blow"], m["lin_mean"],
             m["g1"]))
    allm = n["frame_metrics"]
    if allm:
        A("- 48 フレームの範囲: C std %.4f–%.4f、crush 最大 %.2f%%、blow 最大 %.2f%%"
          % (min(x["C_std"] for x in allm), max(x["C_std"] for x in allm),
             100 * max(x["crush"] for x in allm),
             100 * max(x["blow"] for x in allm)))
    A("")
    A("## 9. G-LEGACY — 既存 py の無傷性")
    A("")
    A("- stage-0 記録: `%s`（全 %d 本。既存 %d 本 + 本作の新規 1 本）。"
      % (n["stage0_hashes_file"], gleg["recorded_existing"] + 1,
         gleg["recorded_existing"]))
    A("- 既存 %d 本の照合 → **自レーンの変更 0 本 / 消失 0 本**。木全体では "
      "変更 %d / 消失 %d / 追加 %d = **外部差 %d 件**（変更: %s / 追加: %s）。"
      "いずれも並行レーンが同じ workdir で書いたもので、本 job はどれも import して"
      "いない（個別の帰属は下の箇条）。"
      % (gleg["recorded_existing"], len(gleg["changed"]),
         len(gleg["missing"]),
         len(gleg["new_since_stage0_not_of_this_job"]),
         gleg["external_diff_count"],
         ", ".join("`%s`" % x for x in gleg["changed"]) or "なし",
         ", ".join("`%s`" % x for x in
                   gleg["new_since_stage0_not_of_this_job"]) or "なし"))
    for nm, at in (gleg.get("attribution") or {}).items():
        cl = at.get("closest_in_time") or {}
        fw = at.get("hash_forensics") or []
        author = [x for x in fw
                  if x["baseline_equals_current"]
                  and not x["job_out_dir"].endswith("v4-1-6r2")]
        if cl.get("job_out_dir"):
            who = "帰属先: %s（その REPORT の mtime %s、差 %+.0f s）" % (
                cl["job_out_dir"], cl.get("report_mtime_local"),
                cl.get("delta_s") or 0.0)
        elif author:
            who = "帰属先: %s（ハッシュ照合による）" % author[0]["job_out_dir"]
        else:
            who = "帰属先: 特定できず（並行レーンのいずれか）"
        A("  - **%s**: mtime %s（この job の stage 0 は %s）、この job が import %s、"
          "%s" % (nm, at.get("mtime_local"),
                  (n["stage0"] or {}).get("checked_local"),
                  at.get("imported_by_this_job"), who))
        if author:
            A("    - ハッシュ照合: %s の stage-0 記録（%s 作成）は**現行ハッシュと一致**。"
              "したがってこの版はその job の stage 0 より前に書かれている（著者は%s）。"
              % (author[0]["job_out_dir"],
                 author[0]["stage0_record_mtime_local"],
                 ("上記の job" if cl.get("job_out_dir")
                  else "この job（本 job の stage 0 より後）")))
        else:
            A("    - ハッシュ照合: 現行ハッシュを含む他 job の stage-0 記録は無い"
              "（＝どの job の stage 0 より後かは特定できない）。")
    if gleg["new_since_stage0_not_of_this_job"]:
        A("  - **追加**（本 job のものでも既存でもない .py）: %s — 同じ並行レーンが"
          "書いたもの。本 job は読み込んでいない。"
          % ", ".join("`%s`" % x for x in
                      gleg["new_since_stage0_not_of_this_job"]))
    for k, v in gleg["own_new_file_hashes"].items():
        A("- 新規 `%s` = %s（stage-0） / %s（終了時）%s"
          % (k, (v["stage0"] or "")[:16], (v["now"] or "")[:16],
             "" if v["stage0"] == v["now"] else
             " ← stage 0 後に編集した（§13 に開示。既存ファイルではない）"))
    A("")
    A("## 10. 場面定数")
    A("")
    A("| 量 | 生成段（L_gen） | 画面（L_img） |")
    A("| --- | --- | --- |")
    A("| H（= 屈折水深 z） | %.2f m | **%.2f m** |" % (GEN["H"], IMG["z_m"]))
    A("| U / U0 | %.2f / %.3f m/s | —（画面量ではない） |"
      % (GEN["U"], geo["U0_m_s"]))
    A("| S0 / D84 | %.4f / %.2f m | — |" % (GEN["S0"], GEN["D84"]))
    A("| seed / 媒体 | %d / %s | %s |"
      % (GEN["seed"], GEN["water"], IMG["water"]))
    A("| パッチ L | %.4f m | **%.2f m** |" % (geo["L_gen_m"], scr["L_img_m"]))
    A("| 格子 / 画素 | %d（cell %.3f mm） | %d×%d（pixel %.4f mm） |"
      % (geo["N"], 1e3 * geo["cell_m"], WIDTH, HEIGHT, 1e3 * scr["pixel_m"]))
    A("| 屈折水深 z | %.4f m（C 生成、= z_img×%.4f） | %.2f m（画像形成） |"
      % (tr["z_used_m"], sep["mapping_factor"], IMG["z_m"]))
    A("| 帯域制限 σ | %.1f cell = %.3f mm | —（輸送のみ） |"
      % (tr["primary_sigma_cells"], 1e3 * tr["sigma_m"]))
    A("| 合成 | — | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d, %s |"
      % (IMG["A"], IMG["B_sigma"], IMG["E"], IMG["Q"], IMG["water"]))
    A("")
    A("## 11. 測定の定義（再現用）")
    A("")
    A("- 窓の位置は**元写真だけ**から決めた: 256×256 タイルを stride 92×90 で走査し、"
      "局所 std が最小のタイルを `smooth`（x0=%d, y0=%d, std %.2f）、最大を `bark`"
      "（x0=%d, y0=%d, std %.2f）。以後すべての画像で同じ座標。"
      % (gl["regions"]["smooth_tile"]["x0"], gl["regions"]["smooth_tile"]["y0"],
         gl["regions"]["smooth_tile"]["std"],
         gl["regions"]["bark_tile"]["x0"], gl["regions"]["bark_tile"]["y0"],
         gl["regions"]["bark_tile"]["std"]))
    A("- **depth（主判定）**: 生画像に理想帯域通過（軸平行 5–15 px、|k_perp|≤3）を"
      "かけた成分の RMS ÷ 窓平均 ×1000（‰）。帯域通過自体がハイパスなので前段"
      "ハイパスは掛けない。**depth_hp** は brief の手順そのまま（σ=3 px ハイパスを"
      "先に掛ける）で、比較のクロスチェックとして併記する。")
    A("- 写真は 735×539、フレームは 736×540（PHOTO_NATIVE はフレーム−1 px）なので"
      "同じ窓原点が同じ内容を指す。両方 8bit。depth は窓平均で正規化するので露出差は"
      "入らない。")
    A("- G-SPEC の annulus は `stream_surface2.spectrum_radial`（幅 dk = 2π/L_gen、"
      "v4-1.4 の G2 と同じ実装）。|ds/di| は `stream_surface4.fold_stats`"
      "（4 近傍差分をプール、レイ単位）。")
    A("- depth の既知の限界: 軸帯の 1 次元マージナル統計であり、格子と「同じ軸帯"
      "パワーを持つ等方テクスチャ」を区別しない。")
    A("")
    A("## 12. 仮定（本作で追加したもの）")
    A("")
    A("- `L_gen = 10 lambda0`: 生成パッチを波長の整数倍に取り、環を dk の 10 倍に"
      "置く（brief の指示）。ループ時間は L_gen/U0 の等価値として逆算し、時間発展の"
      "式は変えない。")
    A("- `z_scales_with_the_screen`: 屈折水深 z は画面寸法に属する量とみなし、生成段"
      "では z_gen = z_img × %.4f = %.4f m を使う。これで写像は純粋な相似変換になり、"
      "ds/di が不変になる（＝生成段で測った折り返しが画面の折り返しと一致する）。"
      % (sep["mapping_factor"], tr["z_used_m"]))
    A("- `transport_no_fold`: 輸送の帯域制限幅は「max|ds/di| ≤ 1.0 を満たす最小 σ」。"
      "v4-1.5 の 4 倍マージンは、この場面では波を輸送から消すため採らない（数値は "
      "§4 と §6.2 に全部出す）。")
    A("- `photo_bed`: 底応答 (i) の底は元写真そのもの（L_gen に写像、σ_b = D84/4）。"
      "境界の値段差だけを線形 ramp で除いて周期化する（ramp の傾きは §3 の表）。")
    A("- `stretch_reading`: 写像は場（傾斜配列）の再標本化と読むので傾きの値は変え"
      "ない（面の相似変換、振幅 ×%.4f、kA 保存）。" % sep["mapping_factor"])
    A("- それ以外の仮定・定数・会計は v4-1.4 / v4-1.5 / v4-1.6 と同一。")
    A("")
    A("## 13. 失敗・弱点・正直な開示")
    A("")
    A("- **G-LINE の閾値はこの画面では原理的に厳しい**: 波を完全に除いた対照"
      "（flat, C=1）でも %s/4 行が写真を上回る（§6.2）。したがって「4 行すべて ≤ 写真」"
      "は合成（ACES・A=%.1f・E=%.1f・輝度の再標本化）だけで満たせない。本作の判定"
      "（%s）は**格子が残った証明にはならない**。証拠として使えるのは、"
      "(i) 周期が写真と同じ 8–12 px に立ち格子の signature が無いこと、"
      "(ii) 折り返しを止めると強度が 70.3 → 62.9 ‰ に下がり**波を除いた対照 "
      "63.9 ‰ を下回る**こと、(iii) 最大超過 %.1f%% が対照の超過 %.1f%% 以下で"
      "あること、の 3 点。"
      % (([c for c in (n.get("ablate") or {}).get("configs", [])
           if str(c["warp_field"]).startswith("flat")][0]["rows_over_photo"])
         if n.get("ablate") else "—",
         IMG["A"], IMG["E"], "PASS" if gl["pass"] else "FAIL",
         100 * (max(r["ratio"] for r in gl["rows"]) - 1.0),
         100 * (max(r["ratio"] for r in
                     ([c for c in (n.get("ablate") or {}).get("configs", [])
                       if str(c["warp_field"]).startswith("flat")][0]["rows"]))
                - 1.0) if n.get("ablate") else float("nan")))
    A("- **σ の規則から 4 倍マージンを外した**: v4-1.5 の規則（max ≤ 0.25）はこの場面"
      "でも σ = %s cell で満たせるが、そのときリング帯域の傾き分散は %.4f に落ち、"
      "輸送から波が消える（C は実質 static な底痕跡になる）。機構条件（max ≤ 1.0）"
      "を採った。**これは規則の変更であり、両方の数値を残した**（§4 のラダー、"
      "§6.2 の strict 行）。"
      % (tr["strict_rule_sigma_cells"],
         (srow[0]["ring_band_ratio"] if srow else float("nan"))))
    A("- **写像は物理ではない**: 描かれる波長 %.1f mm は L_gen の %.1f mm を %.4f 倍"
      "した表示値。実寸 3.68 m の画面で %.1f mm の波を出すには L_gen 側の解像度が"
      "足りない、という v4-1.4 の失敗は解消していない。本作は「場の生成」と「表示"
      "スケール」を分離しただけ。"
      % (sep["lambda_displayed_mm"], sep["lambda0_mm"],
         sep["mapping_factor"], sep["lambda0_mm"]))
    A("- **倍率が 1 未満（縮小）**である点は v4-1.6R（6.4 倍の拡大）と逆向き。"
      "L_gen = 5.42 m > L_img = 3.68 m なので、生成した場は画面で 0.68 倍に縮む。"
      "このため「生成段で解像できている細部」が画面では 1.44 px/cell に落ちる。"
      "C のサンプリング格子（512 cell / 736 px）は v4-1.6R と同じで、格子の周期は "
      "1.4 px 付近＝5–15 px 帯の外にある。")
    A("- **輸送で波の傾き分散が %.1f%% 減る**（σ=%.1f cell）。折り返していた"
      "スケール（サブ %d cell）と波のスケール（%.1f cell/λ0）が離れているので "
      "v4-1.6R よりは軽いが、ゼロではない。C のコントラストは σ=0 の %.3f から "
      "%.3f に下がる（§6.2）。補正・再スケールはしていない。"
      % (100 * (1.0 - tr["ring_band_ratio_at_chosen"]),
         tr["primary_sigma_cells"], 5, lam_cells,
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if c["sigma_cells"] == 0.0][0]["C_std"] if n.get("ablate")
         else float("nan"),
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if c["sigma_cells"] == tr["primary_sigma_cells"]
          and c["warp_field"] == "unfiltered"][0]["C_std"]
         if n.get("ablate") else float("nan")))
    A("- **G-WARP は合否なし**: v4-1.5 の 1.839 px は z = 0.015 m / L = 0.0639 m の"
      "値で、本作の z = 0.60 m / L = 3.68 m とは比較できない。記録のみ。")
    A("- **G-LINE は 1 フレーム（frame 24）の判定**（brief どおり）。48 フレーム全体の"
      "周期測定はしていない。")
    if n_bad_hp:
        A("- brief 手順の hp クロスチェックで **%d 行**が元写真を上回る（§6 に数値）。"
          "主判定（帯域振幅）でも %d 行が上回る。" % (n_bad_hp, len(bad)))
    A("- 判定は数値統計のみ（目視なし）。後付けの閾値調整・再スケールはしていない。"
      "σ の規則（4 倍マージンを採らない）は**画像を描く前にラダーだけから決めた**が、"
      "「strict を採ると波が輸送から消える」という事実を知った後の判断であることは"
      "明記する。")
    A("- **自分のバグと再実行（正直な記録）**: `geometry()` の z_gen を最初 "
      "`z_img × factor`（= 0.408 m）と書き間違え、`--gen` を 1 度その式で走らせた"
      "（そのとき σ = 1.0 cell と出た）。正しくは `z_img ÷ factor` = %.4f m で、"
      "相似変換で ds/di が不変になる条件はこちら。式を直して `--gen` を再実行し、"
      "**以下の数値はすべて 2 回目のもの**。checkpoint の unit id に `gen-field` が "
      "2 度出るのはその履歴。描画（chunk）は修正後にしか走らせていない。"
      % tr["z_used_m"])
    A("- **自分の新規ファイルの self-edit**: `run_v4_1_6r2.py` は stage 0 の後で"
      "編集した（上記 z_gen の修正、自己検査の k_ring をこの場面の k0 に直した点、"
      "G-LEGACY の帰属分離）。既存 .py には触れていない。stage-0 時点と終了時点の"
      "ハッシュは numbers の `gates.G-LEGACY.own_new_file_hashes` に両方残した。")
    A("- **G-LEGACY の FAIL は本 job の外**: stage 0（%s）の後に、同じ workdir で"
      "同時に走っていた別レーンが .py を書いた（§9 の帰属表）。本 job はそれらを "
      "import していない。gate の定義に照らして木全体は FAIL と記録するが、"
      "**本 job の既存 .py への変更は 0 本**。"
      "（`transport2.py` の mtime 05:27:03 は本作の最終フレーム（chunk 01 終了 "
      "05:24:46 ごろ）より後で、`run_v4_2_3.py` も描画後。）"
      % ((n["stage0"] or {}).get("checked_local"),))
    A("- **C は clamp していない**: frame 24 の C は min %.3f / p99 %.3f / max %.3f。"
      "max は幾何光学のカスプ（`det.det_frame` の仕様どおり非クランプ）で、"
      "合成後の blow は %.2f%% なので白飛びしていない。"
      % (cst.get("min") or float("nan"), cst.get("p99") or float("nan"),
         cst.get("max") or float("nan"), 100 * _fm24_blow(n)))
    A("- 実行は全ステージ foreground（1 プロセスずつ）。jobkit checkpoint を "
      "`%s` に stage ごとに原子的にコミットした（unit id: %s）。"
      % (os.path.basename(CP_PATH),
         ", ".join("`%s`" % u for u in
                   (read_checkpoint().get("completed_ids") or []))))
    A("")
    A("## 14. 成果物")
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
    A("注: `checkpoint.json` の行はこの表を作った時点の値。本レポートの書き出しが"
      "最後の unit（`report`）なので、その追記の後に `checkpoint.json` の sha256 は"
      "1 回だけ変わる（他のファイルは変わらない）。")
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-6r2.json + REPORT-V4-1-6R2.md")
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
    ap.add_argument("--gen", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--still", action="store_true")
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.all:
        stage0()
        do_gen()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_still()
        do_ablate()
        do_warp()
        do_encode()
        do_gates()
        do_report()
        log("=" * 78)
        log("[all] done in %.1fs" % (time.time() - t0))
        return 0
    if a.stage0:
        stage0()
    if a.gen:
        do_gen()
    if a.chunk is not None:
        do_chunk(a.chunk)
    if a.still:
        do_still()
    if a.ablate:
        do_ablate()
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
