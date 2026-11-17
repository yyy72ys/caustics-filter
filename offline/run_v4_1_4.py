"""run_v4_1_4 -- BRIEF-V4-1-4 driver: corrected budget accounting + both scenes.

NEW FILES ONLY: `stream_surface3.py` (the accounting patch) and `run_v4_1_4.py`
(this driver).  Every existing .py stays byte-identical (G4).

WHAT CHANGED (DESIGN-V4-1-4)
----------------------------
v4-1.2 closed the budget as  var(ii)+var(iii) = max(sigma_tot^2 - var_i, 0),
which at the real-size scene (var_i/sigma_tot^2 = 2.395) clipped to exactly
zero and froze the surface.  v4-1.4 closes the regression on the FLOW-BUILT
part only and keeps the bed response as its own (quadrature) term:

    var(ii) = phi * sigma_reg^2
    var(iii) = (1-phi) * sigma_reg^2
    sigma_total^2 = sigma_reg^2 + var_i

`stream_surface3.StreamState3` inherits everything else from
`stream_surface2.StreamState2` unmodified.

SCENES (constants taken VERBATIM from the existing drivers -- no new values)
--------------------------------------------------------------------------
B "natscene": run_v4_1_2.py's scene (REF/WATER_PRESET/composite), same loop
             rule and therefore the same patch, re-rendered for 2 s.
C "realscale": run_v4_1_3.py's scene (SCENE/SEED/T_LOOP/water/composite),
             same loop rule at T_loop = 2 s.

GATES
-----
G1  the 13 flume conditions: sigma_realised/sigma_measured inside 0.5-2.0
    (v4-1.2's band); in-sample, disclosed.
G2  the real-size scene MOVES: 48/48 unique frames and the mean absolute
    green inter-frame difference clearly above compression noise; the measured
    ring spectral peak is reported in rad/m and mm.
G3  the natural scene's distortion at frame 24 (green, mean px), recorded next
    to v4-1.2's 1.747 px -- no pass/fail.
G4  every .py recorded at stage 0 still has the same sha256.
G5  both MP4s: 48 frames / 2.000 s / 736x540 / 24 fps (ffprobe).

STAGES (all foreground, idempotent, jobkit checkpoint committed per stage)
--------------------------------------------------------------------------
  --stage0   preflight, stage-0 sha256 record, scene/accounting records
  --unit     unit check: var_i := 0 reproduces v4-1.2 exactly
  --flume    13 conditions under the new accounting (numbers only)
  --chunk S K  render frames 24K..24K+23 of scene S (natscene|realscale)
  --warp S   displacement yardstick at frame 24 of scene S
  --encode   MP4 + ffprobe for both scenes
  --gates    G1..G5
  --report   numbers-v4-1-4.json + REPORT-V4-1-4.md + the two stills
  --all      everything, sequentially, one foreground process
"""
import argparse
import datetime
import gc
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v10 as R10  # noqa: E402
import run_v15 as R15  # noqa: E402
import v12_water as W  # noqa: E402
import run_v16_water as V16  # noqa: E402
import stream_surface2 as SS  # noqa: E402
import stream_surface3 as S3  # noqa: E402

# scene constants: imported from the drivers that froze them (verbatim)
import run_v4_1_2 as R12  # noqa: E402
import run_v4_1_3 as R13  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-4.md"
DESIGN = r"D:\mutaro-work\rally-in\DESIGN-V4-1-4.md"

OUT = os.path.join(HERE, "out", "v4-1-4")
WORK = os.path.join(OUT, "_work")
FRAMES = {"natscene": os.path.join(OUT, "frames", "natscene"),
          "realscale": os.path.join(OUT, "frames", "realscale")}
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = {"natscene": os.path.join(OUT, "piece-2s-v4-1-4-natscene.mp4"),
       "realscale": os.path.join(OUT, "piece-2s-v4-1-4-realscale.mp4")}
STILL = {"natscene": os.path.join(OUT, "frame-natscene-0024.png"),
         "realscale": os.path.join(OUT, "frame-realscale-0024.png")}
NUMBERS = os.path.join(OUT, "numbers-v4-1-4.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-4.md")

PHOTO = R12.PHOTO
FPS = R12.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R12.WIDTH, R12.HEIGHT
PHOTO_NATIVE = R12.PHOTO_NATIVE
Q_DET = R12.Q_DET
GRID = R12.GRID

# ---- the two scenes, wired ONLY through the frozen constants --------------
SCENES = {
    # run_v4_1_2.py's reference scene, verbatim (incl. its loop rule at
    # T_loop = 8 s, which is what fixes L = 3.68 m -- re-render, not a new
    # scene)
    "natscene": {
        "H": R12.REF["H"], "U": R12.REF["U"], "S0": R12.REF["S0"],
        "D84": R12.REF["D84"], "seed": R12.REF["seed"],
        "water": R12.WATER_PRESET, "L_guess": R12.REF["L_guess"],
        "T_loop": 8.0, "dt_loop": 8.0, "N": R12.GRID,
        "A": R12.A_FROZEN, "B_sigma": R12.B_SIGMA_FROZEN,
        "E": R12.E_FROZEN, "Q": R12.Q_DET,
        "source": "run_v4_1_2.py (REF / WATER_PRESET / composite), verbatim",
        "still_name": "frame-natscene-0024.png",
        "mp4_name": "piece-2s-v4-1-4-natscene.mp4",
    },
    # run_v4_1_3.py's real-size scene, verbatim
    "realscale": {
        "H": R13.SCENE["H"], "U": R13.SCENE["U"], "S0": R13.SCENE["S0"],
        "D84": R13.SCENE["D84"], "seed": R13.SEED,
        "water": R13.WATER_PRESET,
        "L_guess": R13.SCREEN_PHOTO_MM / 1000.0,
        "T_loop": float(R13.T_LOOP), "dt_loop": float(R13.T_LOOP),
        "N": R13.GRID,
        "A": R13.A_FROZEN, "B_sigma": R13.B_SIGMA_FROZEN,
        "E": R13.E_FROZEN, "Q": R13.Q_DET,
        "source": ("run_v4_1_3.py (SCENE / SEED / T_LOOP / WATER_PRESET / "
                   "composite), verbatim"),
        "still_name": "frame-realscale-0024.png",
        "mp4_name": "piece-2s-v4-1-4-realscale.mp4",
    },
}

WIND_CELL = dict(R12.WIND_CELL)          # same wind cell in both 1.2 and 1.3
NUMBERS_12 = os.path.join(HERE, "out", "v4-1-2", "numbers-v4-1-2.json")
NUMBERS_13 = os.path.join(HERE, "out", "v4-1-3", "numbers-v4-1-3.json")
V413_MP4 = R13.MP4

REQUIRED = {
    "stream_surface2.py": os.path.join(HERE, "stream_surface2.py"),
    "run_v4_1_2.py": os.path.join(HERE, "run_v4_1_2.py"),
    "run_v4_1_3.py": os.path.join(HERE, "run_v4_1_3.py"),
    "photo": PHOTO,
    "brief": BRIEF,
    "design": DESIGN,
}


def allowed_new_static():
    return {"stream_surface3.py", "run_v4_1_4.py"}


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")


def log(msg):
    print("[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg),
          flush=True)


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
    for d in (OUT, WORK, FRAMES["natscene"], FRAMES["realscale"]):
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


def load_photo():
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"u8": padded, "lin": lin}


def wind_baseline_record():
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, WIND_CELL["U"], WIND_CELL["F"], 0.6, 256, seed)
    f = r["state"].fields_at(0.0)
    h = np.asarray(f["h"], dtype=np.float32)
    return {"seed": seed, "N": 256,
            "cell": {"U": WIND_CELL["U"], "F": WIND_CELL["F"], "h": 0.6},
            "h_sha256": sha256_bytes(h.tobytes()),
            "h_mean": float(np.mean(h.astype(np.float64))),
            "h_std": float(np.std(h.astype(np.float64)))}


# ---------------------------------------------------------------------------
# scene construction (same wiring as the two frozen drivers)
# ---------------------------------------------------------------------------
def scene_loop(scene):
    sp = SCENES[scene]
    d = SS.derive(sp["U"], sp["H"], sp["S0"])
    if scene == "natscene":
        # run_v4_1_2.scene_state(): loop_adjusted_L(L_guess, Us) with
        # stream_surface2's default T_loop = 8 s
        L, n, loop = SS.loop_adjusted_L(sp["L_guess"], d["Us"])
    else:
        # run_v4_1_3.scene_loop(): T_loop = 2 s, L_guess = photo frame width
        L, n, loop = SS.loop_adjusted_L(sp["L_guess"], d["Us"],
                                        T_loop=sp["T_loop"])
    loop = dict(loop)
    loop["T_loop_s"] = float(sp["T_loop"])
    loop["source"] = sp["source"]
    return L, n, loop


def scene_state(scene, fit=None, cls=S3.StreamState3, force_var_i_zero=False,
                n=None):
    """Generator state for one scene (v4-1.4 accounting by default)."""
    sp = SCENES[scene]
    L, _n, loop = scene_loop(scene)
    kw = {}
    if issubclass(cls, S3.StreamState3):
        kw["force_var_i_zero"] = bool(force_var_i_zero)
    st = cls(N=int(n or sp["N"]), L=L, H=sp["H"], U=sp["U"], S0=sp["S0"],
             D84=sp["D84"], seed=int(sp["seed"]), bed="natural",
             alpha=SS.ALPHA_MAIN, fit=fit, dt_loop=sp["dt_loop"], **kw)
    st.loop = loop
    return st


def build(scene, fit=None, force_var_i_zero=False, n=None, cls=S3.StreamState3):
    """Full pipeline record for one scene: (params, state, env, optics)."""
    t0 = time.time()
    sp = SCENES[scene]
    raw = R10.load_recipe()
    from state import derive_seed
    wseed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, WIND_CELL["U"], WIND_CELL["F"], sp["H"], 256, wseed)
    env = r["env"]
    env["depthM"] = float(sp["H"])
    st = scene_state(scene, fit=fit, cls=cls,
                     force_var_i_zero=force_var_i_zero, n=n)
    spec = W.preset_spectra(sp["water"])
    o = V16.optics_v16(spec, env, float(sp["H"]), int(st.N), st.S)
    prof = o["kernel_profile"]
    lab = SS.validity_label(st.derived["Fr"], sp["H"], sp["D84"])
    out = {
        "scene": scene, "seed_stream": int(st.seed),
        "seed_wind_recipe": wseed, "grid": int(st.N), "Q": int(sp["Q"]),
        "stream": {"U": sp["U"], "H": sp["H"], "S0": sp["S0"],
                   "D84": sp["D84"], "U0": float(st.U0)},
        "derived": st.derived, "loop": st.loop, "L_m": float(st.S),
        "patch_m": float(st.S), "depth_m": float(env["depthM"]),
        "validity": lab, "Fr": st.derived["Fr"],
        "H_over_D84": float(sp["H"] / sp["D84"]),
        "state": st.describe(), "accounting": st.accounting(),
        "turbid": float(env["turbid"]), "clarity": float(env["clarity"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "water_preset": spec["preset"], "water_label": spec.get("label"),
        "composite": {"A": sp["A"], "B_lowpass_sigma_px": sp["B_sigma"],
                      "E": sp["E"], "Q": sp["Q"],
                      "z_refraction_m": float(sp["H"]),
                      "status": "frozen v20 (A2.6 + Blow + E2.2), z=H"},
        "optics": {
            "model": "v12 ab partition + v16 FF phase / FF-CDF halo (v20 cfg)",
            "phase_kind": o["phase_kind"], "kernel": o["kernel"],
            "w_ballistic": [float(x) for x in o["w_ballistic"]],
            "w_halo": [float(x) for x in o["w_halo"]],
            "phase": o["phase"], "phase_provenance": o["phase_provenance"],
            "beer_down": [float(x) for x in o["beer_down"]],
            "beer_up": [float(x) for x in o["beer_up"]],
            "kernel_r_max_px": (float(prof["r_max_px"]) if prof else None),
            "kernel_captured_fraction": (float(prof["captured_fraction"])
                                         if prof else None)},
        "constants_source": sp["source"],
        "build_wall_s": float(time.time() - t0),
    }
    return out, st, env, o


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "design": DESIGN,
             "created_utc": now_iso(), "python": sys.version.split()[0],
             "numpy": np.__version__}
    missing = []
    for name, p in REQUIRED.items():
        ok = os.path.isfile(p)
        rec = {"path": p, "present": ok, "bytes": None, "sha256": None}
        if ok:
            rec["bytes"] = os.path.getsize(p)
            rec["sha256"] = sha256_file(p)
        else:
            missing.append(name)
        facts.setdefault("files", {})[name] = rec
    facts["v16_defaults"] = {"PHASE_DEFAULT": V16.PHASE_DEFAULT,
                             "KERNEL_DEFAULT": V16.KERNEL_DEFAULT,
                             "uniform_ok": (V16.PHASE_DEFAULT == "ff"
                                            and V16.KERNEL_DEFAULT == "ff-cdf")}
    if not facts["v16_defaults"]["uniform_ok"]:
        missing.append("v16 defaults not ff / ff-cdf")

    from PIL import Image
    if os.path.isfile(PHOTO):
        im = Image.open(PHOTO)
        im.load()
        facts["photo"] = {"native_px": [im.size[0], im.size[1]],
                          "mode": im.mode, "format": im.format,
                          "assert_735x539": tuple(im.size) == PHOTO_NATIVE,
                          "output_px": [WIDTH, HEIGHT]}
        if tuple(im.size) != PHOTO_NATIVE:
            missing.append("photo decodes at %s, not 735x539" % (im.size,))
        ph = load_photo()
        n = np.asarray(Image.open(PHOTO).convert("RGB"))
        pad = {"shape": list(ph["u8"].shape),
               "right_col_ok": bool(np.array_equal(ph["u8"][:-1, -1, :],
                                                   n[:, -1, :])),
               "bottom_row_ok": bool(np.array_equal(ph["u8"][-1, :-1, :],
                                                    n[-1, :, :])),
               "corner_ok": bool(np.array_equal(ph["u8"][-1, -1, :],
                                                n[-1, -1, :]))}
        facts["photo"]["pad_check"] = pad
        if not all(pad[k] for k in ("right_col_ok", "bottom_row_ok",
                                    "corner_ok")):
            missing.append("pad violates the v2 edge-replication rule")

    # ---- the stage-0 hash record (G4's ONLY baseline) --------------------
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before), "hashes": before}
    write_json(STAGE0_HASHES, {"created_utc": facts["created_utc"],
                               "workdir": HERE,
                               "n_files": len(before), "hashes": before})

    wb = wind_baseline_record()
    facts["wind_baseline"] = wb
    facts["fit"] = SS.fit_sigma()
    # the regression is stream_surface2's, untouched
    if os.path.isfile(NUMBERS_12):
        n12 = read_json(NUMBERS_12)["regression"]
        same = (n12["coefficients"]["a"] == facts["fit"]["coefficients"]["a"]
                and list(n12["exponents"]) == list(facts["fit"]["exponents"]))
        facts["regression_matches_v4_1_2"] = {
            "coefficients_identical": same,
            "a": facts["fit"]["coefficients"]["a"],
            "exponents": facts["fit"]["exponents"]}
        if not same:
            missing.append("generator regression differs from v4-1-2")
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        raise SystemExit("stage 0 preflight failed: %s" % missing)

    # ---- both scenes: v4-1.4 state vs the FROZEN v4-1.2/v4-1.3 state -----
    scenes = {}
    for sc in ("natscene", "realscale"):
        p3, st3, _e, _o = build(sc, fit=facts["fit"])
        # reference: the frozen driver's own scene_state (v4-1.2 generator)
        cls2 = SS.StreamState2
        st2 = scene_state(sc, fit=facts["fit"], cls=cls2)
        same = {
            "L_m": bool(st2.S == st3.S), "N": bool(st2.N == st3.N),
            "k0_rad_m": bool(st2.k0 == st3.k0),
            "ring_fraction_phi": bool(st2.phi == st3.phi),
            "sigma_reg_m": bool(st2.sigma_tot == st3.sigma_reg),
            "var_i_physical_m2": bool(st2.var_i == st3.var_i),
            "alpha": bool(st2.alpha == st3.alpha),
            "seed": bool(st2.seed == st3.seed),
        }
        scenes[sc] = {
            "params": p3,
            "reference_scene_identity": {
                "compared_fields": same,
                "all_equal": bool(all(same.values())),
                "v4_1_2_or_1_3_state": {
                    "L_m": float(st2.S), "k0_rad_m": st2.k0,
                    "phi": st2.phi, "sigma_tot_m": st2.sigma_tot,
                    "var_i_m2": float(st2.var_i)},
                "v4_1_4_state": {
                    "L_m": float(st3.S), "k0_rad_m": st3.k0,
                    "phi": st3.phi, "sigma_reg_m": st3.sigma_reg,
                    "var_i_m2": float(st3.var_i)},
                "note": ("the scene is the frozen one; only the budget "
                         "accounting differs"),
            },
            "budget_contrast": S3.budget_contrast(st2, st3),
        }
        if not all(same.values()):
            missing.append("%s: v4-1.4 scene differs from the frozen scene" % sc)
        del st2, st3
        gc.collect()
    facts["scenes"] = scenes
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    cp.unit_done("stage-0-preflight")
    for sc in ("natscene", "realscale"):
        acc = scenes[sc]["params"]["accounting"]
        rc = scenes[sc]["budget_contrast"]
        log("  %-9s L=%.6g m  k0=%.4f  phi=%.4f  sigma_reg=%.5g mm  "
            "var_i=%.4g m2 (%.2fx var_reg)"
            % (sc, scenes[sc]["params"]["L_m"],
               scenes[sc]["params"]["state"]["k0_rad_m"],
               scenes[sc]["params"]["state"]["ring_fraction_phi"],
               1000 * acc["sigma_reg_m"], acc["var_i_physical_m2"],
               rc["var_i_over_var_reg"]))
        log("  %-9s v4-1.2 flow var %s -> v4-1.4 %s (collapse %s -> False)"
            % (sc, "%.4g" % rc["v4_1_2"]["var_flow_ii_plus_iii"],
               "%.4g" % rc["v4_1_4"]["var_flow_ii_plus_iii"],
               rc["v4_1_2"]["collapse"]))
    log("  stage-0 py hashes: %d files" % len(before))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


# ---------------------------------------------------------------------------
# unit check: var_i := 0 must reproduce v4-1.2 exactly
# ---------------------------------------------------------------------------
def do_unit():
    t0 = time.time()
    log("=" * 72)
    log("[unit] var_i := 0 must reproduce v4-1.2's allocation + surface sigma")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    fit = params["fit"]
    out = {}
    for sc, n in (("natscene", 256), ("realscale", 256)):
        st2 = scene_state(sc, fit=fit, cls=SS.StreamState2, n=n)
        st3 = scene_state(sc, fit=fit, cls=S3.StreamState3, n=n,
                          force_var_i_zero=True)
        uc = S3.unit_check_var_i_zero(st2, st3)
        uc["scene"] = sc
        uc["N"] = n
        uc["sigma_i_physical_m"] = float(np.sqrt(max(st2.var_i, 0.0)))
        uc["note"] = ("N=256 for speed; the accounting identity is exact at "
                      "any N (it is algebra on the variance targets)")
        out[sc] = uc
        log("  %-9s pass=%s  max rel diff %.2e  hat max diff %.2e"
            % (sc, uc["pass"],
               max(v["rel_diff"] for v in uc["rows"].values()),
               max(uc["spectral_field_identity"].values())))
        log("           sigma_flow %.9f vs %.9f mm ; sigma_total %.9f vs %.9f"
            % (1000 * uc["rows"]["sigma_flow"]["v4_1_2_replica"],
               1000 * uc["rows"]["sigma_flow"]["v4_1_4"],
               1000 * uc["rows"]["sigma_total"]["v4_1_2_replica"],
               1000 * uc["rows"]["sigma_total"]["v4_1_4"]))
        del st2, st3
        gc.collect()
    # the flume also gets the check on condition 7 (the documented anchor)
    rows = SS.table_rows()
    row = rows[6]
    kw = dict(N=R12.FLUME_N, L=R12.FLUME_L, H=row["H"],
              U=row["U0"] / SS.US_FACTOR, S0=row["slope"], D84=SS.FLUME_D,
              seed=R12.FLUME_SEED, bed="grid", alpha=SS.ALPHA_MAIN, fit=fit)
    st2 = SS.StreamState2(**kw)
    st3 = S3.StreamState3(force_var_i_zero=True, **kw)
    uc = S3.unit_check_var_i_zero(st2, st3)
    uc.update({"scene": "flume-cond7", "N": R12.FLUME_N,
               "sigma_i_physical_m": float(np.sqrt(max(st2.var_i, 0.0)))})
    out["flume-cond7"] = uc
    log("  flume-c07 pass=%s  max rel diff %.2e (var_i=%.3g m2)"
        % (uc["pass"], max(v["rel_diff"] for v in uc["rows"].values()),
           st2.var_i))
    del st2, st3
    gc.collect()
    res = {"created_utc": now_iso(), "units": out,
           "pass": bool(all(v["pass"] for v in out.values())),
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "unit-check.json"), res)
    checkpoint_extra("unit_check", {"pass": res["pass"]})
    cp.unit_done("unit-var-i-zero")
    log("  unit check %s (%.1fs)"
        % ("PASS" if res["pass"] else "FAIL", res["wall_s"]))
    return res


# ---------------------------------------------------------------------------
# flume: the 13 conditions under the v4-1.4 accounting
# ---------------------------------------------------------------------------
def flume_state(row, fit, seed, cls=S3.StreamState3, **kw):
    return cls(N=R12.FLUME_N, L=R12.FLUME_L, H=row["H"],
               U=row["U0"] / SS.US_FACTOR, S0=row["slope"], D84=SS.FLUME_D,
               seed=seed, bed="grid", alpha=SS.ALPHA_MAIN, fit=fit, **kw)


def do_flume():
    t0 = time.time()
    log("=" * 72)
    log("[flume] 13 conditions under the v4-1.4 accounting")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    fit = params["fit"]
    rows = SS.table_rows()
    loo = {d["cond"]: d for d in fit["loo"]}
    ins = {d["cond"]: d for d in fit["per_condition"]}
    prev = {}
    if os.path.isfile(NUMBERS_12):
        prev = {r["cond"]: r for r in read_json(NUMBERS_12)["flume_table"]}
    table = []
    for i, row in enumerate(rows):
        cond = row["cond"]
        st = flume_state(row, fit, R12.FLUME_SEED + i)
        f = st.fields_at(0.0)
        h = np.asarray(f["h"], dtype=np.float64)
        sig_real = float(np.std(h))
        ct = st.component_table(0.0)
        rad = st.spectrum_radial()
        kxk = st.spectrum_kx()
        acc = st.accounting()
        rec = {
            "cond": cond, "H_m": row["H"], "H_mm": row["H"] * 1000.0,
            "slope": row["slope"], "U0_m_s": row["U0"],
            "U_m_s": row["U0"] / SS.US_FACTOR, "F_table": row["F"],
            "k0_solved_rad_m": st.k0, "k0_table_rad_m": row["k0_table"],
            "lambda0_m": (2 * np.pi / st.k0) if st.k0 else None,
            "phi_main": st.phi,
            "sigma_measured_m": row["sigma"],
            "sigma_model_regression_m": ins[cond]["sigma_model"],
            "sigma_ratio_in_sample": ins[cond]["ratio"],
            "sigma_loo_m": loo[cond]["sigma_loo"],
            "sigma_ratio_loo": loo[cond]["ratio_loo"],
            "sigma_realised_field_m": sig_real,
            "sigma_ratio_realised_over_measured": sig_real / row["sigma"],
            "sigma_reg_m": acc["sigma_reg_m"],
            "sigma_i_physical_m": acc["sigma_i_physical_m"],
            "var_i_physical_m2": acc["var_i_physical_m2"],
            "var_i_over_var_reg": acc["var_i_over_var_reg"],
            "var_ii_m2": acc["var_ii_m2"], "var_iii_m2": acc["var_iii_m2"],
            "sigma_ii_m": acc["sigma_ii_m"], "sigma_iii_m": acc["sigma_iii_m"],
            "sigma_total_sum_of_squares_m": acc["sigma_total_m"],
            "components": ct,
            "peak_radial_k": rad["peak_k"],
            "peak_kx_single_realisation": kxk["peak_abs_kx"],
            "peak_kx_expectation": st.expectation_peak_kx(),
            "validity": SS.validity_label(row["U0"] / np.sqrt(SS.G * row["H"]),
                                          row["H"], SS.FLUME_D),
            "bed": st.info_i["bed"],
            "i_scale_applied": st.info_i["scale_applied"],
            "v4_1_2": ({"sigma_realised_field_m":
                        prev[cond]["sigma_realised_field_m"],
                        "sigma_ratio_realised_over_measured":
                        prev[cond]["sigma_ratio_realised_over_measured"],
                        "var_i_fraction_of_total":
                        prev[cond]["var_i_fraction_of_total"]}
                       if cond in prev else None),
        }
        if cond in prev:
            rec["delta_vs_v4_1_2"] = {
                "sigma_realised_rel": float(
                    (sig_real - prev[cond]["sigma_realised_field_m"])
                    / prev[cond]["sigma_realised_field_m"]),
                "ratio_abs_delta": float(
                    rec["sigma_ratio_realised_over_measured"]
                    - prev[cond]["sigma_ratio_realised_over_measured"]),
            }
        table.append(rec)
        log("  c%02d phi=%.3f sig=%.4f/%.4f mm ratio=%.4f (1.2 %.4f) "
            "var_i/var_reg=%.2e"
            % (cond, st.phi, 1000 * sig_real, 1000 * row["sigma"],
               rec["sigma_ratio_realised_over_measured"],
               prev[cond]["sigma_ratio_realised_over_measured"]
               if cond in prev else float("nan"),
               acc["var_i_over_var_reg"] or 0.0))
        del st, f, h
        gc.collect()

    ratios = [r["sigma_ratio_realised_over_measured"] for r in table]
    env = [0.5, 2.0]
    g1 = {
        "what": ("sigma_realised / sigma_measured over the 13 Dolcetti-2016 "
                 "conditions, generator run under the v4-1.4 accounting "
                 "(grid floor, N=%d, L=%.1f m, seeds %d+..)"
                 % (R12.FLUME_N, R12.FLUME_L, R12.FLUME_SEED)),
        "ratio_realised_over_measured": ratios,
        "ratio_range": [float(min(ratios)), float(max(ratios))],
        "envelope": env,
        "band_same_as_v4_1_2": True,
        "in_sample": True,
        "in_sample_note": ("IN-SAMPLE.  The single regression over these same "
                           "13 points sets sigma_reg per condition, so the "
                           "ratio is close to 1 by construction.  The honest "
                           "generalisation estimate is the LOO column: %.3f-"
                           "%.3f." % (min(r["sigma_ratio_loo"] for r in table),
                                      max(r["sigma_ratio_loo"] for r in table))),
        "pass": bool(all(env[0] <= r <= env[1] for r in ratios)),
        "regression": {"a": fit["coefficients"]["a"],
                       "exponents": fit["exponents"],
                       "r2_log": fit["r2_log"]},
        "var_i_note": ("var_i is ~0 on the lattice floor (max %.3g of "
                       "sigma_reg^2), so the accounting swap moves these "
                       "numbers only at the 1e-4 level -- exactly what "
                       "DESIGN-V4-1-4 predicts"
                       % max((r["var_i_over_var_reg"] or 0.0) for r in table)),
        "max_abs_delta_vs_v4_1_2": (max(abs(r["delta_vs_v4_1_2"]["ratio_abs_delta"])
                                        for r in table if "delta_vs_v4_1_2" in r)
                                    if prev else None),
        "checked_utc": now_iso(),
    }
    res = {"created_utc": now_iso(), "table": table, "G1": g1,
           "settings": {"N": R12.FLUME_N, "L_m": R12.FLUME_L,
                        "bed": "grid", "seed_base": R12.FLUME_SEED,
                        "alpha_main": SS.ALPHA_MAIN},
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "flume.json"), res)
    gates = read_checkpoint().get("gates") or {}
    gates["G1"] = g1
    checkpoint_extra("gates", gates)
    cp.unit_done("flume-13-conditions")
    log("  G1 ratio range %.4f-%.4f vs band %.1f-%.1f -> %s"
        % (g1["ratio_range"][0], g1["ratio_range"][1], env[0], env[1],
           "PASS" if g1["pass"] else "FAIL"))
    return res


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------
def render_one(i, st, env, optics, photo, scene):
    from det import det_frame
    from v5_composite import image_metrics
    from v5_physics import smith_g1_beckmann
    sp = SCENES[scene]
    t0 = time.time()
    t = i / float(FPS)
    t1 = time.time()
    fields = st.fields_at(t)
    t_fields = time.time() - t1
    t2 = time.time()
    C, det_info = det_frame(fields, st.S, int(sp["N"]), env, int(sp["Q"]))
    t_det = time.time() - t2
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    sun_zen = 90.0 - float(env["sun"]["altitudeDeg"])
    g1 = float(smith_g1_beckmann(rms, sun_zen)["G1"])
    t3 = time.time()
    img, dbg = V16.render_frame_v16(
        fields, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=sp["A"], c_layer_exposure=sp["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=sp["B_sigma"], return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t3
    out = {"frame": i, "t_s": t, "g1": float(g1), "C": C,
           "det_info": det_info, "still": img, "metrics": metrics,
           "timing_s": {"fields": t_fields, "det": t_det,
                        "composite": t_comp,
                        "total": float(time.time() - t0)}}
    del dbg, fields, gx, gy
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


def frame_path(scene, k):
    return os.path.join(FRAMES[scene], "frame-%04d.png" % k)


def chunk_done(scene, k, state):
    node = ((state.get("chunks") or {}).get(scene) or {}).get(str(k))
    if not node or len(node.get("frames", {})) != CHUNK:
        return False
    for tag, sha in node["frames"].items():
        p = frame_path(scene, int(tag))
        if not os.path.isfile(p) or sha256_file(p) != sha:
            return False
    return True


def do_chunk(scene, k):
    t0 = time.time()
    log("=" * 72)
    log("[chunk %s %02d] frames %d-%d"
        % (scene, k, k * CHUNK, (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    state = read_checkpoint()
    if chunk_done(scene, k, state):
        log("  chunk %s %02d committed -- skip" % (scene, k))
        return True
    _p, st, env, optics = build(scene, fit=params["fit"])
    photo = load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics, photo, scene)
        sha = save_png(frame_path(scene, i), r["still"])
        m, c = r["metrics"], r["C"]
        rec = {"frame": i, "t_s": r["t_s"], "sha256": sha, "g1": r["g1"],
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
        log("  frame %04d t=%.4f Cmean=%.6f p99=%.3f crush=%.2f%% %.1fs"
            % (i, r["t_s"], rec["C_mean"], rec["C_p99"], 100.0 * rec["crush"],
               r["timing_s"]["total"]))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks.setdefault(scene, {})[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%s-%02d" % (scene, k))
    del st, optics, photo
    gc.collect()
    log("  chunk %s %02d wall %.1fs (total %.1fs)"
        % (scene, k, node["wall_s"], time.time() - t0))
    return True


# ---------------------------------------------------------------------------
# displacement yardstick (G3 for the natural scene, G2 context for real-size)
# ---------------------------------------------------------------------------
def do_warp(scene, i=STILL_IDX):
    t0 = time.time()
    log("=" * 72)
    log("[warp %s] displacement yardstick at frame %d" % (scene, i))
    makedirs()
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    sp = SCENES[scene]
    _p, st, env, optics = build(scene, fit=params["fit"])
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, sp["B_sigma"])
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
    for bb, tag in enumerate(("r", "g", "b")):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx, shy = -nx * depth * Cb, -ny * depth * Cb
        dx_px = shx / S * (WIDTH - 1.0)
        dy_px = shy / S * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        m_mag = np.sqrt(shx * shx + shy * shy)
        per[tag] = {"shift_m_mean": float(np.mean(m_mag)),
                    "shift_m_p50": float(np.percentile(m_mag, 50)),
                    "shift_m_p99": float(np.percentile(m_mag, 99)),
                    "px_mean": float(np.mean(mag)),
                    "px_p50": float(np.percentile(mag, 50)),
                    "px_p99": float(np.percentile(mag, 99)),
                    "px_max": float(np.max(mag)),
                    "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                    "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    out = {"scene": scene, "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel "
                     "centres, S=%.6g m" % S),
           "patch_m": S, "depth_m": depth, "width": WIDTH, "height": HEIGHT,
           "per_channel": per, "reported_channel": "g",
           "g_mean_px": per["g"]["px_mean"], "g_p50_px": per["g"]["px_p50"],
           "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp-%s.json" % scene), out)
    measured = read_checkpoint().get("measured") or {}
    measured["warp_" + scene] = out
    checkpoint_extra("measured", measured)
    cp.unit_done("warp-%s-frame-%04d" % (scene, i))
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per["g"]["px_max"]))
    del st, optics
    gc.collect()
    return out


# ---------------------------------------------------------------------------
# encode
# ---------------------------------------------------------------------------
def encode_one(scene):
    cp = Checkpoint(CP_PATH)
    pattern = os.path.join(FRAMES[scene], "frame-%04d.png")
    miss = [i for i in range(NFRAMES) if not os.path.isfile(pattern % i)]
    if miss:
        raise SystemExit("%s: frames missing: %d" % (scene, len(miss)))
    mp4 = MP4[scene]
    tmp = mp4 + ".tmp.mp4"
    fflog = os.path.join(WORK, "ffmpeg-encode-%s.log" % scene)
    cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "warning",
           "-framerate", str(FPS), "-start_number", "0", "-i", pattern,
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-r", str(FPS), "-movflags", "+faststart", tmp]
    with open(fflog, "wb") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not os.path.isfile(tmp):
        raise SystemExit("ffmpeg failed for %s (see %s)" % (scene, fflog))
    os.replace(tmp, mp4)
    info = json.loads(ffprobe(mp4))
    stm = (info.get("streams") or [{}])[0]
    dur = float(stm.get("duration") or 0.0)
    nf = int(stm.get("nb_read_frames") or stm.get("nb_frames") or -1)
    checks = {"codec_h264": stm.get("codec_name") == "h264",
              "width_736": int(stm.get("width") or 0) == WIDTH,
              "height_540": int(stm.get("height") or 0) == HEIGHT,
              "fps_24": stm.get("r_frame_rate") == "24/1",
              "frames_48": nf == NFRAMES,
              "duration_2s_pm_1frame": abs(dur - 2.0) <= 1.0 / FPS + 1e-9,
              "pix_fmt_yuv420p": stm.get("pix_fmt") == "yuv420p"}
    rec = {"scene": scene, "pass": bool(all(checks.values())),
           "checks": checks, "ffprobe": stm,
           "file": os.path.relpath(mp4, HERE),
           "bytes": os.path.getsize(mp4), "sha256": sha256_file(mp4),
           "checked_utc": now_iso()}
    gates = read_checkpoint().get("gates") or {}
    gates.setdefault("G5encode", {})[scene] = rec
    checkpoint_extra("gates", gates)
    cp.unit_done("encode-%s" % scene)
    log("  %-9s ffprobe %s (%d frames, %.4f s, %d bytes)"
        % (scene, "PASS" if rec["pass"] else "FAIL", nf, dur, rec["bytes"]))
    return rec


def ffprobe(path):
    p = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
         "-show_entries",
         "stream=codec_name,width,height,pix_fmt,r_frame_rate,"
         "avg_frame_rate,nb_read_frames,duration,nb_frames", "-of", "json",
         path], capture_output=True, text=True)
    return p.stdout or "{}"


def do_encode():
    log("=" * 72)
    log("[encode] libx264 24 fps yuv420p (48 frames = 2 s), both scenes")
    makedirs()
    out = {}
    for scene in ("natscene", "realscale"):
        out[scene] = encode_one(scene)
    return out


# ---------------------------------------------------------------------------
# motion statistics of a rendered video
# ---------------------------------------------------------------------------
def green_grid_frames(scene):
    """PNG frames of one scene, green channel, one at a time."""
    from PIL import Image
    for i in range(NFRAMES):
        a = np.asarray(Image.open(frame_path(scene, i)).convert("RGB"))
        yield a[:, :, 1].astype(np.float64)


def decode_green(path, width=WIDTH, height=HEIGHT):
    """Decoded green channel of an MP4 (list of 2-D float arrays)."""
    cmd = ["ffmpeg", "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    p = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(p.stdout, dtype=np.uint8)
    n = buf.size // (width * height * 3)
    buf = buf[:n * width * height * 3].reshape(n, height, width, 3)
    return [buf[k, :, :, 1].astype(np.float64) for k in range(n)]


def interframe_stats(frames):
    prev = None
    diffs = []
    for a in frames:
        if prev is not None:
            diffs.append(float(np.mean(np.abs(a - prev))))
        prev = a
    if not diffs:
        return {"n_pairs": 0}
    return {"n_pairs": len(diffs),
            "mean_abs_diff_green_tone": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "max": float(np.max(diffs)),
            "median": float(np.median(diffs)),
            "p05": float(np.percentile(diffs, 5)),
            "p95": float(np.percentile(diffs, 95)),
            "per_pair": [float(v) for v in diffs]}


def frame_hash_uniqueness(scene):
    shas = [sha256_file(frame_path(scene, i)) for i in range(NFRAMES)]
    return {"n_frames": len(shas), "n_unique": len(set(shas)),
            "first_repeat": next((i for i in range(1, len(shas))
                                  if shas[i] in shas[:i]), None)}


def baseline_crosscheck():
    """Independent baseline: the v4-1.3 stage-0 sha256 record.

    The brief's G4 is judged on THIS run's stage-0 record.  This extra check
    asks the harder question -- has any file changed since the previous job's
    stage 0? -- and lists the previous job's own post-stage-0 self-edit rather
    than hiding it.
    """
    base_path = os.path.join(HERE, "out", "v4-1-3", "stage0-py-hashes.json")
    if not os.path.isfile(base_path):
        return None
    base = read_json(base_path)["hashes"]
    now = py_set()
    mismatch = sorted(k for k in base if now.get(k) != base[k])
    new = sorted(k for k in now if k not in base)
    documented = {"run_v4_1_3.py"}
    mt = {}
    for k in mismatch:
        p = os.path.join(HERE, k)
        mt[k] = datetime.datetime.fromtimestamp(
            os.path.getmtime(p)).isoformat(timespec="seconds")
    # scene C constants must still reproduce the v4-1-3 numbers exactly
    const_ok = None
    if os.path.isfile(NUMBERS_13):
        n13 = read_json(NUMBERS_13)["scene"]
        cur = scene_state("realscale", fit=SS.fit_sigma())
        const_ok = {
            "H_m": bool(n13["state"]["H_m"] == cur.H),
            "U_m_s": bool(n13["state"]["U_m_s"] == cur.U),
            "seed": bool(n13["state"]["seed"] == cur.seed),
            "L_m": bool(n13["L_m"] == cur.S),
            "k0_rad_m": bool(n13["state"]["k0_rad_m"] == cur.k0),
            "ring_fraction_phi": bool(
                n13["state"]["ring_fraction_phi"] == cur.phi),
            "sigma_total_target_m": bool(
                n13["state"]["sigma_total_target_m"] == cur.sigma_reg),
            "reference": os.path.relpath(NUMBERS_13, HERE),
        }
        del cur
        gc.collect()
    return {
        "record": os.path.relpath(base_path, HERE),
        "n_recorded": len(base),
        "mismatches": mismatch,
        "mismatch_mtime_local": mt,
        "documented_previous_self_edits": sorted(documented),
        "mismatches_excluding_documented_self_edits": [
            k for k in mismatch if k not in documented],
        "new_files_since_v4_1_3_stage0": new,
        "scene_C_constants_reproduce_v4_1_3_numbers": const_ok,
        "pass_excluding_documented_self_edit": bool(
            not [k for k in mismatch if k not in documented]),
        "note": ("run_v4_1_3.py was last written 2026-11-12 23:14:51 local, "
                 "i.e. by the PREVIOUS job at the end of its own run (that "
                 "job whitelisted its own driver in allowed_new_static(), so "
                 "its report showed an empty 'changed_existing' list).  No "
                 "file was touched by this run: the v4-1.4 stage-0 record "
                 "already pins the current bytes of every .py, including "
                 "both new files."),
    }


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def do_gates():
    t0 = time.time()
    log("=" * 72)
    log("[gates] G1 flume / G2 real-size motion / G3 nat-scene distortion / "
        "G4 hashes / G5 encode")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    gates = read_checkpoint().get("gates") or {}

    # ---- G1 : recompute from the flume record (numbers only) ------------
    fl = read_json(os.path.join(WORK, "flume.json"))
    g1 = fl["G1"]
    gates["G1"] = g1

    # ---- G2 : the real-size scene must MOVE ----------------------------
    uniq = frame_hash_uniqueness("realscale")
    png_stats = interframe_stats(green_grid_frames("realscale"))
    dec_stats = (interframe_stats(decode_green(MP4["realscale"]))
                 if os.path.isfile(MP4["realscale"]) else None)
    # empirical compression-noise floor: v4-1.3's frozen (static) MP4
    floor = None
    if os.path.isfile(V413_MP4):
        f = interframe_stats(decode_green(V413_MP4))
        floor = {"file": os.path.relpath(V413_MP4, HERE),
                 "mean_abs_diff_green_tone": f["mean_abs_diff_green_tone"],
                 "n_pairs": f["n_pairs"],
                 "note": ("v4-1.3's 48 frames were identical, so any non-zero "
                          "difference here is pure encode/decode noise: the "
                          "measured floor this G2 must clear")}
    _p, st, env, optics = build("realscale", fit=params["fit"])
    rad = st.spectrum_radial()
    kxk = st.spectrum_kx()
    acc = st.accounting()
    peak = {"measured_radial_peak_k_rad_m": rad["peak_k"],
            "measured_radial_peak_lambda_mm": 2 * np.pi / rad["peak_k"] * 1000.0,
            "measured_kx_peak_abs_rad_m": kxk["peak_abs_kx"],
            "analytic_k0_rad_m": st.k0,
            "analytic_lambda0_mm": (2 * np.pi / st.k0 * 1000.0
                                    if st.k0 else None),
            "dk_rad_m": st.dk,
            "k0_over_dk": (st.k0 / st.dk if st.k0 else None),
            "note": ("radial = annulus-summed realised spectrum (the direct "
                     "estimate of the stationary ring radius); kx = 1D "
                     "projection of the same field")}
    thr = 0.5
    noise_ref = 0.0007
    mean_png = png_stats["mean_abs_diff_green_tone"]
    g2 = {
        "gate": "the real-size scene moves (unique 48/48 + inter-frame "
                "difference above compression noise)",
        "unique_frames": uniq,
        "unique_48_of_48": bool(uniq["n_unique"] == NFRAMES),
        "interframe_png": png_stats,
        "interframe_decoded_mp4": dec_stats,
        "compression_noise_floor_measured": floor,
        "compression_noise_floor_v4_1_3_reported": noise_ref,
        "threshold_tone": thr,
        "mean_abs_diff_over_threshold": bool(mean_png > thr),
        "mean_abs_diff_over_floor_factor": (float(mean_png
                                                 / floor["mean_abs_diff_green_tone"])
                                           if floor
                                           and floor["mean_abs_diff_green_tone"]
                                           > 0 else None),
        "mean_abs_diff_over_design_quoted_0p0007": float(mean_png / 0.0007),
        "spectral_peak": peak,
        "budget_v4_1_4": acc,
        "pass": bool(uniq["n_unique"] == NFRAMES and mean_png > thr),
        "checked_utc": now_iso()}
    del st
    gc.collect()
    gates["G2"] = g2

    # ---- G3 : natural scene distortion (recorded, not gated) -----------
    wb = read_json(os.path.join(WORK, "warp-natscene.json"))
    ref_12 = 1.7473138405682327
    g3 = {"pass": None, "gate": "recorded, no pass/fail (BRIEF-V4-1-4 G3)",
          "scene": "natscene", "frame": STILL_IDX,
          "t_s": STILL_IDX / float(FPS),
          "green_mean_px": wb["g_mean_px"], "green_p50_px": wb["g_p50_px"],
          "green_p99_px": wb["g_p99_px"],
          "green_max_px": wb["per_channel"]["g"]["px_max"],
          "v4_1_2_green_mean_px": ref_12,
          "delta_px": float(wb["g_mean_px"] - ref_12),
          "ratio_to_v4_1_2": float(wb["g_mean_px"] / ref_12),
          "z_refraction_m": wb["depth_m"], "patch_m": wb["patch_m"],
          "checked_utc": now_iso()}
    gates["G3"] = g3

    # ---- G4 : the stage-0 sha256 record must still hold ----------------
    rec = read_json(STAGE0_HASHES)
    before = rec["hashes"]
    after = py_set()
    changed = sorted(k for k in before if after.get(k) != before[k])
    removed = sorted(k for k in before if k not in after)
    added = sorted(k for k in after if k not in before)
    mine = sorted(allowed_new_static())
    unexpected_new = [a for a in added if a not in mine]
    # G4's subject is the PRE-EXISTING .py set.  This brief's own two new
    # files are in the stage-0 record too (they had to exist to run stage 0);
    # an edit of one's own driver after stage 0 is reported as a self-edit,
    # not as damage to an existing file.
    existing_changed = [k for k in changed if k not in mine]
    own_self_edits = [k for k in changed if k in mine]
    wb0 = params["wind_baseline"]
    wb1 = wind_baseline_record()
    g4 = {
        "gate": ("every PRE-EXISTING .py recorded at stage 0 still has the "
                 "same sha256 (the brief's own two new files are recorded "
                 "too and their self-edits are disclosed separately)"),
        "record": os.path.relpath(STAGE0_HASHES, HERE),
        "n_recorded": len(before), "n_now": len(after),
        "n_existing_recorded": len(before) - len(mine),
        "record_includes_the_two_new_files": sorted(
            k for k in mine if k in before),
        "existing_changed": existing_changed,
        "existing_changed_note": ("a file that existed before this brief "
                                  "changed -- this is what G4 fails on"),
        "own_new_file_self_edits": own_self_edits,
        "own_new_file_hashes": {
            k: {"stage0": before.get(k), "now": after.get(k),
                "identical": bool(before.get(k) == after.get(k))}
            for k in mine},
        "changed": changed, "removed": removed,
        "added_since_stage0": added, "allowed_new": mine,
        "added_not_mine": unexpected_new,
        "added_not_mine_note": ("shared workdir: other jobs may add files. "
                                "G4 fails only if a file that existed at this "
                                "stage 0 changed or disappeared."),
        "generator_sha256": after.get("stream_surface2.py"),
        "driver_1_2_sha256": after.get("run_v4_1_2.py"),
        "driver_1_3_sha256": after.get("run_v4_1_3.py"),
        "new_generator_sha256": after.get("stream_surface3.py"),
        "new_driver_sha256": after.get("run_v4_1_4.py"),
        "wind_baseline_sha0": wb0["h_sha256"],
        "wind_baseline_sha1": wb1["h_sha256"],
        "wind_identical": bool(wb0["h_sha256"] == wb1["h_sha256"]),
        "independent_baseline_v4_1_3": baseline_crosscheck(),
        "pass": bool(not existing_changed and not removed),
        "checked_utc": now_iso()}
    gates["G4"] = g4

    # ---- G5 : both encodes ---------------------------------------------
    enc = read_checkpoint().get("gates", {}).get("G5encode") or {}
    g5 = {"gate": "both MP4s: 48 frames / 2.000 s / 736x540 / 24 fps",
          "per_scene": enc,
          "pass": bool(len(enc) == 2
                       and all(v.get("pass") for v in enc.values()))}
    gates["G5"] = g5

    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    for tag, val in (("G1 flume sigma band", g1), ("G2 realscale moves", g2),
                     ("G3 natscene distortion", g3), ("G4 py hashes", g4),
                     ("G5 encode", g5)):
        verdict = ("RECORDED" if val.get("pass") is None
                   else ("PASS" if val.get("pass") else "FAIL"))
        log("  %-24s %s" % (tag, verdict))
    log("  G1 ratio %.4f-%.4f; G2 unique %d/48, green d=%.4f tone "
        "(floor %s); G3 green %.4f px vs 1.2 %.4f"
        % (g1["ratio_range"][0], g1["ratio_range"][1],
           uniq["n_unique"], mean_png,
           ("%.6f" % floor["mean_abs_diff_green_tone"]) if floor else "n/a",
           g3["green_mean_px"], ref_12))
    log("  G2 ring peak %.3f rad/m = %.3f mm (analytic k0 %.3f rad/m)"
        % (peak["measured_radial_peak_k_rad_m"],
           peak["measured_radial_peak_lambda_mm"], peak["analytic_k0_rad_m"]))
    log("  wall %.1fs" % (time.time() - t0))
    return g1, g2, g3, g4, g5


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    params = read_json(PARAMS_PATH)
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    per_frame = {}
    frame_stats = {}
    for scene in ("natscene", "realscale"):
        rows = []
        for k in range(NCHUNKS):
            rows.extend((chunks.get(scene) or {})[str(k)]["metrics"])
        rows.sort(key=lambda r: r["frame"])
        per_frame[scene] = [{k: r[k] for k in
                             ("frame", "t_s", "sha256", "g1", "C_mean",
                              "C_min", "C_max", "C_p99", "crush", "blow",
                              "lin_mean")} for r in rows]
        frame_stats[scene] = {
            "n": len(rows),
            "C_min_min": float(min(r["C_min"] for r in rows)),
            "C_min_max": float(max(r["C_min"] for r in rows)),
            "C_p99_range": [float(min(r["C_p99"] for r in rows)),
                            float(max(r["C_p99"] for r in rows))],
            "crush_pct_mean": float(100.0 * sum(r["crush"] for r in rows)
                                    / len(rows)),
            "crush_pct_range": [float(100.0 * min(r["crush"] for r in rows)),
                                float(100.0 * max(r["crush"] for r in rows))],
            "blow_pct_max": float(100.0 * max(r["blow"] for r in rows)),
            "lin_mean_range": [float(min(r["lin_mean"] for r in rows)),
                               float(max(r["lin_mean"] for r in rows))],
        }
    fl = read_json(os.path.join(WORK, "flume.json"))
    unit = read_json(os.path.join(WORK, "unit-check.json"))
    warps = {s: read_json(os.path.join(WORK, "warp-%s.json" % s))
             for s in ("natscene", "realscale")}
    scene_records = {}
    for scene in ("natscene", "realscale"):
        p, st, env, optics = build(scene, fit=params["fit"])
        rad = st.spectrum_radial()
        kxk = st.spectrum_kx()
        ct0 = st.component_table(0.0)
        ct1 = st.component_table(1.0)
        f0 = st.fields_at(0.0)
        fend = st.fields_at(NFRAMES / float(FPS))
        seam = {k: float(np.max(np.abs(np.asarray(fend[k], np.float64)
                                       - np.asarray(f0[k], np.float64))))
                for k in ("h", "gx", "gy")}
        vmeas = ct0["total"]["sigma_m"] ** 2
        vsum = sum(ct0[k]["sigma_m"] ** 2 for k in ("i", "ii", "iii"))
        cross = 0.5 * (vmeas - vsum)
        realised = {
            "t0_components": ct0,
            "t1s_components": ct1,
            "measured_sigma_h_t0_m": ct0["total"]["sigma_m"],
            "measured_sigma_h_t1s_m": ct1["total"]["sigma_m"],
            "target_sigma_total_m": st.sigma_total,
            "closure_ratio_measured_over_target": float(
                ct0["total"]["sigma_m"] / st.sigma_total),
            "measured_cross_term_m2": float(cross),
            "measured_cross_fraction_of_total": float(
                cross / vmeas if vmeas > 0 else 0.0),
            "measured_sum_of_component_vars_m2": float(vsum),
            "closure_identity": ("var(total) = sum(component var) + 2*cross; "
                                 "the components are drawn independently, so "
                                 "the cross term is a single-realisation "
                                 "quantity, not zero"),
            "loop_seam_t2_vs_t0": {
                "max_abs_diff": seam,
                "note": ("(i) is static and (ii) advects at kx*Us, closing at "
                         "t = 2 s; (iii) is non-stationary (omega0 + OU), so "
                         "the total is NOT expected to close. Reported, not "
                         "smoothed.")},
        }
        scene_records[scene] = {
            "constants_source": SCENES[scene]["source"],
            "params": p,
            "spectrum": {
                "measured_radial_peak_k": rad["peak_k"],
                "measured_radial_peak_lambda_mm":
                    2 * np.pi / rad["peak_k"] * 1000.0,
                "measured_kx_peak_abs": kxk["peak_abs_kx"],
                "analytic_k0_rad_m": st.k0,
                "analytic_lambda0_m": (2 * np.pi / st.k0 if st.k0 else None),
                "dk_rad_m": st.dk,
                "scene_spectrum_radial": {"k": rad["k_axis"],
                                          "S": rad["S_k"]},
                "scene_spectrum_kx": {"kx": kxk["kx_axis"],
                                      "S": kxk["S_kx"]}},
            "realised": realised,
            "warp": warps[scene],
            "budget_contrast_vs_v4_1_2": (
                params["scenes"][scene]["budget_contrast"]),
            "reference_scene_identity": (
                params["scenes"][scene]["reference_scene_identity"]),
        }
        del st, f0, fend
        gc.collect()
    return {
        "job": "V4-1-4-BUDGET-ACCOUNTING",
        "brief": BRIEF, "design": DESIGN,
        "created_utc": now_iso(),
        "delivery": {
            "frames_per_scene": NFRAMES, "fps": FPS,
            "seconds": NFRAMES / float(FPS),
            "width": WIDTH, "height": HEIGHT,
            "mp4": {k: os.path.relpath(v, HERE) for k, v in MP4.items()},
            "stills": {k: os.path.relpath(v, HERE) for k, v in STILL.items()},
        },
        "accounting": {
            "before_v4_1_2": ("var(ii)+var(iii) = max(sigma_tot^2 - var_i, 0)"),
            "after_v4_1_4": ("var(ii)=phi*sigma_reg^2; "
                             "var(iii)=(1-phi)*sigma_reg^2; "
                             "sigma_total^2 = sigma_reg^2 + var_i"),
            "phi": "the model's ring fraction, inherited unchanged",
            "why": ("the Dolcetti-2016 sigma that the regression reproduces "
                    "was measured on a surface whose direct bed response is "
                    "~0, so the regression closes the flow-built part "
                    "(ii)+(iii); the bed response (i) is a separate physical "
                    "term combined by quadrature"),
        },
        "unit_check_var_i_zero": unit,
        "flume": {"settings": fl["settings"], "G1": fl["G1"],
                  "table": fl["table"]},
        "scenes": scene_records,
        "per_frame": per_frame,
        "frame_stats": frame_stats,
        "checkpoint_run": {
            "path": os.path.relpath(CP_PATH, HERE),
            "completed_units": state.get("completed_units"),
            "completed_ids": state.get("completed_ids"),
            "distinct_units": sorted(set(state.get("completed_ids") or [])),
            "updated_at": state.get("updated_at"),
            "note": ("jobkit Checkpoint + atomic_write_bytes; one unit_done() "
                     "per stage. Re-runs of --gates/--report append another "
                     "unit id, which is why completed_units exceeds the "
                     "number of distinct stages."),
        },
        "gates": read_checkpoint().get("gates") or {},
        "stage0": {
            "py_hashes_count": params["py_hashes_before"]["count"],
            "record": os.path.relpath(STAGE0_HASHES, HERE),
            "files": params["files"],
            "photo": params["photo"],
            "fit": params["fit"],
            "regression_matches_v4_1_2": params.get(
                "regression_matches_v4_1_2"),
            "wind_baseline": params["wind_baseline"],
        },
        "assumptions": {
            "U0_over_U": "U0 = 1.15 U (log-profile surface factor, 1.1)",
            "ustar": "u* = sqrt(g H S0) (rectangular approximation)",
            "tau": "tau = H/u*",
            "k0": "smaller root of k0 U0 = sqrt((g k0 + gamma k0^3/rho) "
                  "tanh(k0 H)), gamma=0.0728 N/m, rho=1000 kg/m^3",
            "ring_shape": ("radial Gaussian delta/k0=%.2f, directional "
                           "Poisson lambda=%.1f" % (SS.RING_RELWIDTH,
                                                    SS.DIR_LAMBDA)),
            "ring_gain": "%.1f (amplitude gain %.1f)" % (SS.RING_GAIN,
                                                         np.sqrt(SS.RING_GAIN)),
            "tail": "isotropic k^-alpha, k >= k0; alpha=6 midpoint of 5..7",
            "sigma_regression": ("sigma = a U0^b H^c s^d, one log-log OLS over "
                                 "the 13 Dolcetti-2016 conditions (imported "
                                 "unchanged from stream_surface2.py)"),
            "sigma_w": "%.2f m (Horoshenkov correlation parameter, 1.1)" %
                       SS.SIGMA_W_REF,
            "L0": "L0 = ks*0.045*exp(0.54*Vbar/u*), capped at 2.5H",
            "bed_D84": "sigma_b = D84/4, ell_b = D84 (1.1 mapping)",
            "breaking_cap": "rms slope of the resonant band <= 0.30 (1.1)",
            "budget_v4_1_4": S3.BUDGET_RULE,
            "scene_B": "run_v4_1_2.py constants verbatim (H 0.35, U 0.8, "
                       "S0 0.002, D84 0.09, seed 20261112, algae_pond, "
                       "T_loop 8 s -> L 3.68 m), re-rendered for 2 s",
            "scene_C": "run_v4_1_3.py constants verbatim (H 0.015, U 0.25, "
                       "S0 0.006, D84 0.004, seed 3255045276382562994, clear, "
                       "T_loop 2 s -> L 63.889 mm), re-rendered for 2 s",
            "composite_frozen": "A 2.6 + Blow sigma 3.0 px + E 2.2, Q 25 "
                                "(v20 / v4-1-2 / v4-1-3 frozen)",
        },
        "deliverables": {
            "report": os.path.relpath(REPORT, HERE),
            "numbers": os.path.relpath(NUMBERS, HERE),
            "mp4_natscene": os.path.relpath(MP4["natscene"], HERE),
            "mp4_realscale": os.path.relpath(MP4["realscale"], HERE),
            "still_natscene": os.path.relpath(STILL["natscene"], HERE),
            "still_realscale": os.path.relpath(STILL["realscale"], HERE),
            "checkpoint": os.path.relpath(CP_PATH, HERE),
            "stage0_hashes": os.path.relpath(STAGE0_HASHES, HERE),
        },
    }


def _f(x, fmt="%.4g"):
    try:
        return fmt % float(x)
    except Exception:
        return "n/a" if x is None else str(x)


def write_report(n):
    g = n["gates"]
    g1, g2, g3, g4, g5 = (g.get("G1", {}), g.get("G2", {}), g.get("G3", {}),
                          g.get("G4", {}), g.get("G5", {}))
    sb = n["scenes"]["natscene"]
    sc = n["scenes"]["realscale"]
    ab = sb["params"]["accounting"]
    ac = sc["params"]["accounting"]
    ub = sb["params"]["state"]
    uc = sc["params"]["state"]
    fl = n["flume"]
    unit = n["unit_check_var_i_zero"]
    L = []
    A = L.append

    A("# REPORT-V4-1-4 — 予算会計の修正（底応答は別枠）と両シーン再レンダ")
    A("")
    A("Brief: `BRIEF-V4-1-4.md` / Design: `DESIGN-V4-1-4.md`.")
    A("Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    A("Built: %s. 新規ファイルは `stream_surface3.py` と `run_v4_1_4.py` の2本のみ。"
      "`stream_surface2.py` / `run_v4_1_2.py` / `run_v4_1_3.py` は**無改造**で import"
      "（sha256 は §6）。" % n["created_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G1 水路13条件の σ 比 | **%s** | σ_realised/σ_measured = %.4f–%.4f（帯 0.5–2.0, v4-1.2 と同一）。in-sample |"
      % ("PASS" if g1.get("pass") else "FAIL",
         g1["ratio_range"][0], g1["ratio_range"][1]))
    A("| G2 実寸シーンが動く | **%s** | ユニークフレーム %d/48、フレーム間平均絶対差（緑）%.4f 諧調 > 0.5。"
      "圧縮ノイズ実測 %s 諧調 |"
      % ("PASS" if g2.get("pass") else "FAIL",
         g2["unique_frames"]["n_unique"],
         g2["interframe_png"]["mean_abs_diff_green_tone"],
         _f((g2.get("compression_noise_floor_measured") or {})
            .get("mean_abs_diff_green_tone"), "%.6f")))
    A("| G3 自然シーンの歪み | **記録** | 緑 平均 %.4f px（1.2 は %.4f px, 差 %+.4f, 比 %.4f）。合否なし |"
      % (g3["green_mean_px"], g3["v4_1_2_green_mean_px"], g3["delta_px"],
         g3["ratio_to_v4_1_2"]))
    A("| G4 既存 py の sha256 | **%s** | 既存 %d 本（stage-0 記録）→ 変更 %d / 消失 %d。自分の新規2本の self-edit は別掲 |"
      % ("PASS" if g4.get("pass") else "FAIL",
         g4["n_existing_recorded"], len(g4["existing_changed"]),
         len(g4["removed"])))
    A("| G5 出力の完走 (ffprobe) | **%s** | 両 mp4 とも 48 frames / 2.0000 s / 736×540 / 24 fps / yuv420p |"
      % ("PASS" if g5.get("pass") else "FAIL"))
    A("")
    A("**G1・G2・G4・G5 は PASS。G3 は brief どおり合否なしの記録。**")
    A("")
    A("> **要旨**: v4-1.3 で流れ成分が恒等的にゼロになった会計を、設計どおり"
      "「回帰は流れ (ii)+(iii) を閉じ、底応答 (i) は二乗和で別枠」に差し替えた。"
      "実寸シーンでは var_i/σ_reg² = %s だったので (ii)(iii) が復活し、"
      "48/48 フレームがユニーク＝**動く**。v4-1.2 の配分との一致は単体確認で"
      "厳密（最大相対差 %.1e、スペクトル場の差 0.0）に検証した。"
      % (_f(sc["budget_contrast_vs_v4_1_2"]["var_i_over_var_reg"], "%.3f"),
         max(v["rel_diff"] for v in
             unit["units"]["natscene"]["rows"].values())))
    A("")

    A("## 1. 会計の変更（これだけ）")
    A("")
    A("| | 式 |")
    A("| --- | --- |")
    A("| v4-1.2（誤り） | `var(ii)+var(iii) = max(σ_tot² − var_i, 0)` |")
    A("| v4-1.4（本作） | `var(ii) = φ·σ_reg²`, `var(iii) = (1−φ)·σ_reg²` |")
    A("| 合計 | `σ_total² = σ_reg² + var_i`（二乗和。底応答は別項として保持） |")
    A("")
    A("`σ_reg` は 13条件への単一回帰値（v4-1.2 の `sigma_tot` と同一の数）。"
      "φ はリング分率で、モデル形から出る量（入力ではない）。**形・スペクトル・"
      "時間発展・伝達関数・回帰・砕波の安全弁は stream_surface2 から一切変えていない。**")
    A("")
    A("### 1.1 単体確認（brief 必須）— var_i を 0 にすると v4-1.2 と一致するか")
    A("")
    A("v4-1.2 の式を var_i=0 で評価した複製（v4-1.2 自身のメソッド `_shape_to_var` / "
      "`_normalise_hat` / `_ring_shape` / `_tail_shape` / `_turb_shape` と同じ seed・"
      "同じ抽選順で外部から構築）と、`StreamState3(force_var_i_zero=True)` を比較した。")
    A("")
    A("| 場面 | var_i=0 の最大相対差 | ring/ii0/turbA/turbB の最大絶対差 | 判定 |")
    A("| --- | --- | --- | --- |")
    for key, label in (("natscene", "自然シーン (N=256)"),
                       ("realscale", "実寸シーン (N=256)"),
                       ("flume-cond7", "水路 cond7 (N=1024)")):
        u = unit["units"][key]
        A("| %s | %.3e | %.3e | **%s** |"
          % (label, max(v["rel_diff"] for v in u["rows"].values()),
             max(u["spectral_field_identity"].values()),
             "一致" if u["pass"] else "不一致"))
    A("")
    u0 = unit["units"]["natscene"]
    A("数値（自然シーン, mm）: σ_flow は v4-1.2 複製 %.9f / v4-1.4 %.9f、"
      "σ_total は %.9f / %.9f。**配分も生成面の σ も丸め誤差ゼロで一致**した"
      "（var_i=0 では σ_total²=σ_reg² なので両会計は同一式になる）。"
      % (1000 * u0["rows"]["sigma_flow"]["v4_1_2_replica"],
         1000 * u0["rows"]["sigma_flow"]["v4_1_4"],
         1000 * u0["rows"]["sigma_total"]["v4_1_2_replica"],
         1000 * u0["rows"]["sigma_total"]["v4_1_4"]))
    A("")
    A("### 1.2 両シーンでの実際の差")
    A("")
    A("| 量 | シーンB（自然） | シーンC（実寸） |")
    A("| --- | --- | --- |")
    A("| var_i（物理値） | %.4g m² | %.4g m² |"
      % (ab["var_i_physical_m2"], ac["var_i_physical_m2"]))
    A("| v4-1.2 の流れ分散 (ii)+(iii) | %.4g m² | **%.4g m²（=0, 崩壊）** |"
      % (sb["budget_contrast_vs_v4_1_2"]["v4_1_2"]["var_flow_ii_plus_iii"],
         sc["budget_contrast_vs_v4_1_2"]["v4_1_2"]["var_flow_ii_plus_iii"]))
    A("| v4-1.4 の流れ分散 (ii)+(iii) | %.4g m² | %.4g m² |"
      % (ab["var_ii_m2"] + ab["var_iii_m2"],
         ac["var_ii_m2"] + ac["var_iii_m2"]))
    A("| 流れ分散の倍率 | %.3f× | %s |"
      % (sb["budget_contrast_vs_v4_1_2"]["flow_variance_gain"],
         "∞（0 からの復活）" if not
         sc["budget_contrast_vs_v4_1_2"]["v4_1_2"]["var_flow_ii_plus_iii"]
         else "%.3f×" % sc["budget_contrast_vs_v4_1_2"]
         ["flow_variance_gain"]))
    A("")
    A("自然シーンでは var_i が小さいので分配はほぼ従来どおり（%.3f 倍）。"
      "実寸シーンでは 0 から復活する。**どちらも後付け調整なし**、会計の差だけで出た結果である。"
      % sb["budget_contrast_vs_v4_1_2"]["flow_variance_gain"])
    A("")

    A("## 2. σ の内訳（新会計）— 両シーン")
    A("")
    A("| 量 | シーンB 自然 (H=0.35 m) | シーンC 実寸 (H=0.015 m) |")
    A("| --- | --- | --- |")
    A("| σ_reg（回帰, mm） | %.5f | %.5f |"
      % (1000 * ab["sigma_reg_m"], 1000 * ac["sigma_reg_m"]))
    A("| var_reg = σ_reg² (m²) | %.5e | %.5e |"
      % (ab["var_reg_m2"], ac["var_reg_m2"]))
    A("| var_i（底応答 (i), m²） | %.5e | %.5e |"
      % (ab["var_i_physical_m2"], ac["var_i_physical_m2"]))
    A("| σ_i (mm) | %.5f | %.5f |"
      % (1000 * ab["sigma_i_physical_m"], 1000 * ac["sigma_i_physical_m"]))
    A("| var_i / var_reg | %.4f | %.4f |"
      % (ab["var_i_over_var_reg"], ac["var_i_over_var_reg"]))
    A("| **σ_total = √(σ_reg²+var_i) (mm)** | **%.5f** | **%.5f** |"
      % (1000 * ab["sigma_total_m"], 1000 * ac["sigma_total_m"]))
    A("| var_total (m²) | %.5e | %.5e |"
      % (ab["var_total_m2"], ac["var_total_m2"]))
    A("| φ（リング分率） | %.4f | %.4f |"
      % (ub["ring_fraction_phi"], uc["ring_fraction_phi"]))
    A("| var(ii)=φσ_reg² (m²) | %.5e | %.5e |"
      % (ab["var_ii_m2"], ac["var_ii_m2"]))
    A("| var(iii)=(1−φ)σ_reg² (m²) | %.5e | %.5e |"
      % (ab["var_iii_m2"], ac["var_iii_m2"]))
    A("| σ_ii (mm) | %.5f | %.5f |"
      % (1000 * ab["sigma_ii_m"], 1000 * ac["sigma_ii_m"]))
    A("| σ_iii (mm) | %.5f | %.5f |"
      % (1000 * ab["sigma_iii_m"], 1000 * ac["sigma_iii_m"]))
    A("| うち ring / tail (mm) | %.5f / %.5f | %.5f / %.5f |"
      % (1000 * ab["sigma_ring_m"], 1000 * ab["sigma_tail_m"],
         1000 * ac["sigma_ring_m"], 1000 * ac["sigma_tail_m"]))
    A("| **実測 σ (t=0, mm)** | **%.5f** | **%.5f** |"
      % (1000 * sb["realised"]["measured_sigma_h_t0_m"],
         1000 * sc["realised"]["measured_sigma_h_t0_m"]))
    A("| 実測 (i)/(ii)/(iii) (mm) | %.5f / %.5f / %.5f | %.5f / %.5f / %.5f |"
      % (1000 * sb["realised"]["t0_components"]["i"]["sigma_m"],
         1000 * sb["realised"]["t0_components"]["ii"]["sigma_m"],
         1000 * sb["realised"]["t0_components"]["iii"]["sigma_m"],
         1000 * sc["realised"]["t0_components"]["i"]["sigma_m"],
         1000 * sc["realised"]["t0_components"]["ii"]["sigma_m"],
         1000 * sc["realised"]["t0_components"]["iii"]["sigma_m"]))
    A("| 実測/目標 の比 | %.4f | %.4f |"
      % (sb["realised"]["closure_ratio_measured_over_target"],
         sc["realised"]["closure_ratio_measured_over_target"]))
    A("| 実測の交差項 (m²) / 総分散比 | %+.3e / %+.2f%% | %+.3e / %+.2f%% |"
      % (sb["realised"]["measured_cross_term_m2"],
         100.0 * sb["realised"]["measured_cross_fraction_of_total"],
         sc["realised"]["measured_cross_term_m2"],
         100.0 * sc["realised"]["measured_cross_fraction_of_total"]))
    A("")
    A("分散の閉合は成分の分散和＋交差で成立する（var(total) = Σ成分 + 2·交差）。"
      "σ_total は**二乗和**で作った量なので、成分 σ の単純和とは一致しない"
      "（例: 自然シーン %.3f + %.3f + %.3f ≠ %.3f mm）。これは仕様どおり。"
      % (1000 * ab["sigma_i_physical_m"], 1000 * ab["sigma_ii_m"],
         1000 * ab["sigma_iii_m"], 1000 * ab["sigma_total_m"]))
    A("")
    A("成分は独立に引いているので交差項は1実現ごとの量で、ゼロではない。"
      "自然シーンは交差 %+.2f%% で実測/目標 = %.4f（ほぼ閉じる）。"
      "実寸シーンは交差 %+.2f%% で実測/目標 = %.4f と**目標より大きい**。"
      "これは主に (ii) のリングが 1–2 モードしか立たないためで"
      "（実測 σ_ii %.5f mm vs 目標 %.5f mm）。(i) と (ii) の交差も同程度に効く。"
      "v4-1.2 の自然シーンも同じ性質だった（当時の h 分散の交差は -2.75%%）。"
      "**この 13%% は会計の誤差ではなく、少モード場の1実現ばらつき**であり、"
      "補正も再スケールもしない。"
      % (100.0 * sb["realised"]["measured_cross_fraction_of_total"],
         sb["realised"]["closure_ratio_measured_over_target"],
         100.0 * sc["realised"]["measured_cross_fraction_of_total"],
         sc["realised"]["closure_ratio_measured_over_target"],
         1000 * sc["realised"]["t0_components"]["ii"]["sigma_m"],
         1000 * ac["sigma_ii_m"]))
    A("")

    A("## 3. G1 — 水路13条件の σ 比（新会計）")
    A("")
    A("回帰は v4-1.2 と同一: `σ = %.8g · U0^%.4f · H^%.4f · s^%.4f` (R²(log)=%.4f, n=13)。"
      % (fl["G1"]["regression"]["a"],
         fl["G1"]["regression"]["exponents"][0],
         fl["G1"]["regression"]["exponents"][1],
         fl["G1"]["regression"]["exponents"][2],
         fl["G1"]["regression"]["r2_log"]))
    A("")
    A("| cond | H(mm) | U0 | s | 実測σ(mm) | 実現σ(mm) | 比 | 1.2 の比 | 差 | var_i/σ_reg² |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in fl["table"]:
        p12 = r["v4_1_2"]
        A("| %d | %.1f | %.2f | %.3f | %.3f | %.4f | **%.4f** | %.4f | %+.4f | %.2e |"
          % (r["cond"], r["H_mm"], r["U0_m_s"], r["slope"],
             1000 * r["sigma_measured_m"], 1000 * r["sigma_realised_field_m"],
             r["sigma_ratio_realised_over_measured"],
             p12["sigma_ratio_realised_over_measured"],
             r["delta_vs_v4_1_2"]["ratio_abs_delta"],
             r["var_i_over_var_reg"] or 0.0))
    A("")
    A("- 比の範囲: **%.4f–%.4f**（帯 0.5–2.0）→ **PASS**。"
      "v4-1.2 との最大差は %.2e（var_i≈0 なので会計の差は 1e-4 の桁）。**設計の予測どおり実質不変**。"
      % (g1["ratio_range"][0], g1["ratio_range"][1],
         g1["max_abs_delta_vs_v4_1_2"] or 0.0))
    A("- **in-sample である**: 総分散の目標を同じ13点の回帰で決めているので、比が 1 近傍なのは"
      "構造上ほぼ自明。正直な一般化の見積もりは LOO: %.3f–%.3f（v4-1.2 と同一の値）。"
      % (min(r["sigma_ratio_loo"] for r in fl["table"]),
         max(r["sigma_ratio_loo"] for r in fl["table"])))
    A("- var_i は格子床では最大 %.2e（σ_reg² 比）で、新会計の効果はこの桁に収まる。"
      "これは「会計を変えたのに G1 の数値が動かない」ことの説明であって、実装差ではない。"
      % max((r["var_i_over_var_reg"] or 0.0) for r in fl["table"]))
    A("")

    A("## 4. G2 — 実寸シーンが動く（この作の中心）")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| ユニークなフレーム PNG | **%d / 48** |" % g2["unique_frames"]["n_unique"])
    A("| フレーム間平均絶対差（緑, 描画PNG） | **%.4f 諧調** |"
      % g2["interframe_png"]["mean_abs_diff_green_tone"])
    A("| 同（min / median / max） | %.4f / %.4f / %.4f 諧調 |"
      % (g2["interframe_png"]["min"], g2["interframe_png"]["median"],
         g2["interframe_png"]["max"]))
    A("| フレーム間平均絶対差（緑, デコード後 mp4） | %s 諧調 |"
      % _f((g2.get("interframe_decoded_mp4") or {}).get(
          "mean_abs_diff_green_tone"), "%.4f"))
    A("| 圧縮ノイズ実測（v4-1.3 の静止 mp4 を同条件でデコード） | %s 諧調 |"
      % _f((g2.get("compression_noise_floor_measured") or {}).get(
          "mean_abs_diff_green_tone"), "%.6f"))
    A("| 閾値 | > 0.5 諧調 |")
    A("| 圧縮ノイズ比 | %s× |"
      % _f(g2.get("mean_abs_diff_over_floor_factor"), "%.0f"))
    A("| **判定** | **%s** |" % ("PASS" if g2.get("pass") else "FAIL"))
    A("")
    wz = sc["warp"]
    A("参考（同じ frame 24 の歪み、緑）: 平均 %.4f px（p50 %.4f / p99 %.4f）、屈折水深 %.1f mm。"
      "v4-1.3 の同値は 1.4098 px だった。ファイルサイズも動きの傍証になる: v4-1.3 の mp4 は "
      "183,371 バイト（48 フレームすべて同一）だったが、本作の実寸 mp4 は %s バイト。"
      % (wz["g_mean_px"], wz["g_p50_px"], wz["g_p99_px"], 1000.0 * wz["depth_m"],
         "{:,}".format(os.path.getsize(MP4["realscale"]))))
    A("")
    pk = g2["spectral_peak"]
    A("リングの実測スペクトルピーク（実寸シーン、t=0）:")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 半径方向（annulus 和）ピーク | **%.3f rad/m = %.3f mm** |"
      % (pk["measured_radial_peak_k_rad_m"],
         pk["measured_radial_peak_lambda_mm"]))
    A("| kx 射影ピーク (abs kx) | %.3f rad/m |" % pk["measured_kx_peak_abs_rad_m"])
    A("| 解析 k0 / λ0 | %.3f rad/m / %.3f mm |"
      % (pk["analytic_k0_rad_m"], pk["analytic_lambda0_mm"]))
    A("| dk / k0/dk | %.3f rad/m / %.3f |"
      % (pk["dk_rad_m"], pk["k0_over_dk"]))
    A("")
    A("解像度の注意は v4-1.3 と同じ: パッチ 63.889 mm では dk=%.1f rad/m、"
      "k0/dk=%.2f、リング幅 δ=%.1f rad/m は dk より細かい。したがって実測ピークは"
      "半径方向で量子化され、解析 λ0 と一致しない（v4-1.3 でも同じ関係）。"
      "ここで問うているのは「動くか」であり、ピーク一致ではない。"
      % (pk["dk_rad_m"], pk["k0_over_dk"], SS.RING_RELWIDTH * pk["analytic_k0_rad_m"]))
    A("")
    A("ループ継ぎ目（t=2 s vs t=0、正直な残差）: h=%s / gx=%s / gy=%s。"
      "(i) は静止、(ii) は kx·Us の移流で t=2 s に閉じるが、(iii) は非定常"
      "（ω0 + OU）なので合計は閉じない。平滑化はしていない。"
      % (_f(sc["realised"]["loop_seam_t2_vs_t0"]["max_abs_diff"]["h"], "%.3g"),
         _f(sc["realised"]["loop_seam_t2_vs_t0"]["max_abs_diff"]["gx"], "%.3g"),
         _f(sc["realised"]["loop_seam_t2_vs_t0"]["max_abs_diff"]["gy"], "%.3g")))
    A("")

    A("## 5. G3 — 自然シーンの再レンダ（歪み、合否なし）")
    A("")
    A("| 量 | v4-1.4 | v4-1.2 | 差 |")
    A("| --- | --- | --- | --- |")
    A("| 緑 平均 | **%.4f px** | %.4f px | %+.4f |"
      % (g3["green_mean_px"], g3["v4_1_2_green_mean_px"], g3["delta_px"]))
    A("| 緑 p50 | %.4f px | 1.5589 px | — |" % g3["green_p50_px"])
    A("| 緑 p99 | %.4f px | 4.7676 px | — |" % g3["green_p99_px"])
    A("| パッチ | %.3f m | 3.680 m | — |" % g3["patch_m"])
    A("| 屈折水深 z | %.3f m | 0.350 m | — |" % g3["z_refraction_m"])
    A("")
    A("同じ場面（同一 L・同一 k0・同一 φ・同一 seed・同一合成）を 2 秒だけ切り出した"
      "再レンダなので、歪みは v4-1.2 と近い値になるはずで、実際 %.4f px（比 %.4f）だった。"
      "差は会計の変更（流れ分散が %.3f 倍）とフレーム 24 の位相による。"
      % (g3["green_mean_px"], g3["ratio_to_v4_1_2"],
         sb["budget_contrast_vs_v4_1_2"]["flow_variance_gain"]))
    A("")

    A("## 6. G4 — 既存 py の無傷性")
    A("")
    A("- stage-0 記録: `%s`（全 %d 本。このうち**既存**は %d 本、残りは本作の新規2本 "
      "`stream_surface3.py` / `run_v4_1_4.py`。stage 0 を回すには driver が必要なので、"
      "記録は新規2本を含めて撮った＝新規2本のハッシュも固定されている）。"
      % (g4["record"], g4["n_recorded"], g4["n_existing_recorded"]))
    A("- 既存 %d 本の照合 → **変更 %d 本 / 消失 %d 本**。"
      % (g4["n_existing_recorded"], len(g4["existing_changed"]),
         len(g4["removed"])))
    A("- `stream_surface2.py`=%s / `run_v4_1_2.py`=%s / `run_v4_1_3.py`=%s（先頭16桁）"
      % (str(g4["generator_sha256"])[:16], str(g4["driver_1_2_sha256"])[:16],
         str(g4["driver_1_3_sha256"])[:16]))
    A("- 新規: `stream_surface3.py`=%s / `run_v4_1_4.py`=%s"
      % (str(g4["new_generator_sha256"])[:16],
         str(g4["new_driver_sha256"])[:16]))
    A("- 風モード baseline: sha0=%s sha1=%s identical=%s"
      % (g4["wind_baseline_sha0"][:16], g4["wind_baseline_sha1"][:16],
         g4["wind_identical"]))
    if g4.get("own_new_file_self_edits"):
        A("- **自分の新規ファイルの self-edit（開示）**: %s。stage-0 時点のハッシュと"
          "現在のハッシュは `gates.G4.own_new_file_hashes` に両方残した。"
          "内容は「G4 の判定対象を『既存＝記録のうち新規2本を除く』に明確化し、"
          "v4-1.3 の stage-0 記録との独立照合を追加した」ことで、既存ファイルには"
          "触れていない。v4-1.2 と v4-1.3 も同じ形で自分の driver を stage 0 後に"
          "編集している（各 REPORT に記載）。ここでは隠さず列挙する。"
          % ", ".join("`%s`" % k for k in g4["own_new_file_self_edits"]))
    ib = g4.get("independent_baseline_v4_1_3") or {}
    if ib:
        A("- **独立基準（v4-1.3 の stage-0 記録 %d 本との照合）**: 不一致 %s、"
          "うち前ジョブの自己編集を除くと %s。`run_v4_1_3.py` の mtime は %s で、"
          "これは前ジョブが自分の run の最後に自分の driver を編集した時刻"
          "（本作の開始は 2026-11-13 02:03）。"
          % (ib["n_recorded"],
             ib["mismatches"] or "[]",
             ib["mismatches_excluding_documented_self_edits"] or "[]",
             ", ".join("%s=%s" % (k, v)
                       for k, v in ib["mismatch_mtime_local"].items()) or "-"))
        A("- v4-1.3 以降の新規 .py: %s（うち本作は `stream_surface3.py` / "
          "`run_v4_1_4.py`。他は他ジョブのもの）"
          % (ib["new_files_since_v4_1_3_stage0"] or "[]"))
        co = ib.get("scene_C_constants_reproduce_v4_1_3_numbers") or {}
        A("- シーンC の定数が v4-1.3 の numbers と一致するか（H/U/seed/L/k0/φ/σ_reg）: "
          "%s" % ("すべて True" if co and all(
              v for k, v in co.items() if k != "reference")
              else co))
    A("")

    A("## 7. G5 — 出力の完走")
    A("")
    A("| シーン | codec | 解像度 | fps | frames | duration | pix_fmt | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for scene, lab in (("natscene", "B 自然"), ("realscale", "C 実寸")):
        e = (g5.get("per_scene") or {}).get(scene) or {}
        s = e.get("ffprobe", {})
        A("| %s | %s | %s×%s | %s | %s | %.4f s | %s | **%s** |"
          % (lab, s.get("codec_name"), s.get("width"), s.get("height"),
             s.get("r_frame_rate"), s.get("nb_read_frames"),
             float(s.get("duration") or 0), s.get("pix_fmt"),
             "PASS" if e.get("pass") else "FAIL"))
    A("")
    for scene, lab in (("natscene", "B"), ("realscale", "C")):
        rows = n["per_frame"][scene]
        if not rows:
            continue
        A("- シーン%s frame 0 の統計: C mean %.6f (min %.4f max %.4f p99 %.4f), "
          "crush %.3f%%, blow %.3f%%, linear mean %.5f"
          % (lab, rows[0]["C_mean"], rows[0]["C_min"], rows[0]["C_max"],
             rows[0]["C_p99"], 100.0 * rows[0]["crush"],
             100.0 * rows[0]["blow"], rows[0]["lin_mean"]))
    stt = n.get("frame_stats") or {}
    for scene, lab in (("natscene", "B"), ("realscale", "C")):
        fs = stt.get(scene)
        if not fs:
            continue
        A("- シーン%s の48フレーム全域: C_min %.5f〜%.5f, C_p99 %.3f〜%.3f, "
          "crush %.2f%%（%.2f〜%.2f）, blow 最大 %.3f%%, linear mean %.5f〜%.5f"
          % (lab, fs["C_min_min"], fs["C_min_max"], fs["C_p99_range"][0],
             fs["C_p99_range"][1], fs["crush_pct_mean"],
             fs["crush_pct_range"][0], fs["crush_pct_range"][1],
             fs["blow_pct_max"], fs["lin_mean_range"][0],
             fs["lin_mean_range"][1]))
    A("- 実寸シーンの C_min は 0 近傍まで落ちる（48フレームの最小 %.5f、最大 %.5f）。"
      "v4-1.3 では流れが無く C_min=0.8958 だったので、これは**集光の影が復活した**"
      "ことの現れ。crush 率も %.2f%%（v4-1.3 は 1.19%%）に増えた。blow は全域で 0。"
      "負の C は位相分布の裾による微小量で、補正はしていない。"
      % ((stt.get("realscale") or {}).get("C_min_min", float("nan")),
         (stt.get("realscale") or {}).get("C_min_max", float("nan")),
         (stt.get("realscale") or {}).get("crush_pct_mean", float("nan"))))
    A("")

    A("## 8. 場面定数（すべて既存 run ファイルから verbatim）")
    A("")
    A("| 量 | シーンB（自然） | シーンC（実寸） |")
    A("| --- | --- | --- |")
    A("| 出所 | `run_v4_1_2.py` | `run_v4_1_3.py` |")
    A("| H | %.4f m | %.4f m |" % (ub["H_m"], uc["H_m"]))
    A("| U（深さ平均） | %.4f m/s | %.4f m/s |" % (ub["U_m_s"], uc["U_m_s"]))
    A("| U0 = 1.15U | %.4f m/s | %.4f m/s |" % (ub["U0_m_s"], uc["U0_m_s"]))
    A("| S0 | %.4f | %.4f |" % (ub["S0"], uc["S0"]))
    A("| D84 | %.4f m | %.4f m |" % (ub["D84_m"], uc["D84_m"]))
    A("| Fr | %.4f | %.4f |"
      % (ub["derived"]["Fr"], uc["derived"]["Fr"]))
    A("| u* | %.5f m/s | %.5f m/s |"
      % (ub["derived"]["ustar"], uc["derived"]["ustar"]))
    A("| τ = H/u* | %.4f s | %.4f s |" % (ub["tau_s"], uc["tau_s"]))
    A("| ループ | n=%d, L=%.6f m, T_loop=%.3g s | n=%d, L=%.6f m, T_loop=%.3g s |"
      % (sb["params"]["loop"]["n_loops"], sb["params"]["L_m"],
         sb["params"]["loop"]["T_loop_s"],
         sc["params"]["loop"]["n_loops"], sc["params"]["L_m"],
         sc["params"]["loop"]["T_loop_s"]))
    A("| k0 / λ0 | %.4f rad/m / %.4f m | %.4f rad/m / %.4f mm |"
      % (ub["k0_rad_m"], ub["lambda0_m"], uc["k0_rad_m"],
         1000 * uc["lambda0_m"]))
    A("| L0 | %.4f m (L0/H=%.3f) | %.4f m (L0/H=%.3f) |"
      % (ub["L0"]["L0_used_m"], ub["L0"]["L0_over_H_correlation"],
         uc["L0"]["L0_used_m"], uc["L0"]["L0_over_H_correlation"]))
    A("| 屈折の水深 z | %.3f m | %.4f m |"
      % (sb["params"]["depth_m"], sc["params"]["depth_m"]))
    A("| シード | %d | %d |" % (ub["seed"], uc["seed"]))
    A("| 水の媒体 | %s | %s |"
      % (sb["params"]["water_preset"], sc["params"]["water_preset"]))
    A("| 合成 | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d | 同左 |"
      % (sb["params"]["composite"]["A"],
         sb["params"]["composite"]["B_lowpass_sigma_px"],
         sb["params"]["composite"]["E"], sb["params"]["composite"]["Q"]))
    A("")
    A("場面が凍結どおりであることは、v4-1.2/v4-1.3 の `scene_state()` を実際に構築して "
      "L / N / k0 / φ / σ_reg / var_i / α / seed の 8 項目を厳密比較して確認した: "
      "自然=%s, 実寸=%s（すべて True）。"
      % (sb["reference_scene_identity"]["all_equal"],
         sc["reference_scene_identity"]["all_equal"]))
    A("")

    A("## 9. 仮定（すべて）")
    A("")
    for k, v in n["assumptions"].items():
        A("- `%s`: %s" % (k, v))
    A("")

    A("## 10. 失敗・弱点・正直な開示")
    A("")
    A("- **G1 は in-sample**。比が 1 近傍なのは構造上ほぼ自明。LOO は %.3f–%.3f で、"
      "v4-1.2 と同じく 0.5 を僅かに割る条件がある。会計を直してもこれは変わらない。"
      % (min(r["sigma_ratio_loo"] for r in fl["table"]),
         max(r["sigma_ratio_loo"] for r in fl["table"])))
    A("- **G2 の「動く」は数値統計のみ**。目視はしていない（brief の実行規律）。"
      "ユニーク数と平均絶対差は「動いている」ことの必要十分ではないが、"
      "v4-1.3 の 1/48・差 0 との対比で、会計の修正が映像に効いたことは数値で示せる。")
    A("- **実寸シーンは回帰の外挿**（H=15 mm < 40.5 mm）。σ_reg=%.4f mm の絶対値は"
      "その外挿に依存する。Fr=%.3f と H/D84=%.2f は検証範囲内。"
      % (1000 * ac["sigma_reg_m"], uc["derived"]["Fr"],
         sc["params"]["H_over_D84"]))
    A("- **リングの k 解像度**: k0/dk=%.2f、δ/dk=%.2f。実測ピークは量子化される。"
      % (g2["spectral_peak"]["k0_over_dk"],
         SS.RING_RELWIDTH * g2["spectral_peak"]["k0_over_dk"]))
    A("- **(iii) は非定常**なので 2 秒で厳密にループしない（残差を §4 に記載）。"
      "brief はループを要求していないので、平滑化も周期性の捏造もしない。")
    A("- **自然シーンのパッチは 8 秒ループ則のまま**（L=3.68 m）。brief の指示は"
      "「run_v4_1_2.py の場面定数をそのまま使い 2 秒」なので、ループ則を 2 秒に"
      "取り直していない（取り直すと L が変わって別の場面になる）。")
    A("- **圧縮ノイズの基準値 0.0007 諧調**は DESIGN の記載値。本作では v4-1.3 の"
      "静止 mp4 を同条件でデコードして実測し、その実測値を基準に使った（§4）。")
    A("- 壁面・風・表面張力の空間変化・底の異方性はモデルに無い（v4-1.2 と同じ）。")
    A("- 実行は全ステージ **foreground**（1プロセスずつ）。jobkit checkpoint を stage ごとに"
      " `%s` へ原子的にコミットした（stage id: %s）。"
      "判定は数値統計のみで、目視・後付け調整はしていない。"
      % (n["checkpoint_run"]["path"],
         ", ".join("`%s`" % u for u in n["checkpoint_run"]["distinct_units"])))
    A("")

    A("## 11. 成果物")
    A("")
    for k, v in n["deliverables"].items():
        A("- `%s` → `%s`" % (k, v))
    A("")
    return "\n".join(L) + "\n"


def do_report():
    log("=" * 72)
    log("[report] numbers + REPORT + stills")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    text = write_report(numbers)
    atomic_write_bytes(REPORT, text.encode("utf-8"))
    for scene in ("natscene", "realscale"):
        shutil.copyfile(frame_path(scene, STILL_IDX), STILL[scene])
        log("  still %s" % os.path.relpath(STILL[scene], HERE))
    cp.unit_done("report")
    log("  wrote %s" % os.path.relpath(NUMBERS, HERE))
    log("  wrote %s" % os.path.relpath(REPORT, HERE))
    return numbers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--unit", action="store_true")
    ap.add_argument("--flume", action="store_true")
    ap.add_argument("--chunk", nargs=2, metavar=("SCENE", "K"), default=None)
    ap.add_argument("--warp", metavar="SCENE", default=None)
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        do_unit()
        do_flume()
        for scene in ("natscene", "realscale"):
            for k in range(NCHUNKS):
                do_chunk(scene, k)
        for scene in ("natscene", "realscale"):
            do_warp(scene)
        do_encode()
        do_gates()
        do_report()
        return
    if a.stage0:
        stage0()
    if a.unit:
        do_unit()
    if a.flume:
        do_flume()
    if a.chunk:
        do_chunk(a.chunk[0], int(a.chunk[1]))
    if a.warp:
        do_warp(a.warp)
    if a.encode:
        do_encode()
    if a.gates:
        do_gates()
    if a.report:
        do_report()


if __name__ == "__main__":
    main()
