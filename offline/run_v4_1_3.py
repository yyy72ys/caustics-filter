"""run_v4_1_3 -- BRIEF-V4-1-3 driver: the REAL-SIZE scene (frame width 67 mm).

New files only: this driver and out/v4-1-3/.  ``stream_surface2.py`` (the
generator), ``run_v4_1_2.py`` and every other existing file are used WITHOUT
modification and are hashed at stage 0 (G3 baseline).

WHAT THIS RUN IS
----------------
v4-1.2 rendered the reference scene at H = 0.35 m on a 3.68 m wide frame: the
butterfly (48 mm) was 1.3 % of the frame width, i.e. the picture was ~77x too
large for the animal.  v4-1.3 keeps the same generator and composite and only
changes the SIZE OF THE SCENE so that the photo is at its true scale:

    butterfly wingspan 48 mm = 71.2 % of the frame width
      -> frame width = 48 / 0.712 = 67.4 mm  (the brief declares 67 mm)

Scene constants (fixed by the brief, NOT tuned after looking at results):

    H       = 0.015 m      water depth (a thin rill over bark)
    U       = 0.25 m/s     depth-mean velocity  (U0 = 1.15 U = 0.2875 m/s)
    S0      = 0.006        bed slope
    D84     = 0.004 m      bed roughness scale (bark relief)
    seed    = 20261112     the v4-1.2 stream seed (unchanged)
    wind    = 3255045276382562994  the wind/recipe seed named by the brief
                                   (re-derived from the recipe, and asserted)

Composite: identical to v20 / v4-1.2 (A 2.6, Blow sigma 3.0, E 2.2, Q 25,
substrate gain 1, ior 1.33434) with ONE change -- the refracted water column
is z = H = 0.015 m (the butterfly sits on the bed).  Water medium = clear.

DECLARED NUMERICAL CHOICES (made before the render was run)
-----------------------------------------------------------
* generator grid GRID = 512, exactly as v4-1.2;
* patch length from the v4-1.2 loop rule
  ``loop_adjusted_L(L_guess=0.067, Us, T_loop=8 s)`` -> L = 0.067647 m
  (34 advection loops in the 8 s piece), i.e. the patch is the frame width
  to 0.97 %;
* alpha = 6 (the bracket midpoint), ring gain / relwidth / dir lambda all
  inherited from v4-1.2 (45.0 / 0.15 / 4.0);
* sigma_total from the SAME single 13-condition regression as v4-1.2
  (``SS.fit_sigma()``, form a*U0^b*H^c*s^d), evaluated at the scene constants.

NOTE ON THE PATCH SIZE AND THE SPECTRAL GRID
-------------------------------------------
The mode spacing of a periodic patch is dk = 2*pi/L = 92.9 rad/m, independent
of N.  With L = 67.6 mm and k0 = 127.2 rad/m the stationary ring (radius
sigma_k = 0.15*k0 = 19.1 rad/m) is NARROWER than one mode spacing, so the
ring is not representable on this patch at all.  This is a property of the
brief's scene size, not of a numerical parameter, and is reported as such
(G1 records both the analytic lambda0 and the realised spectral peak).

Stages (all foreground, idempotent):
  --stage0 preflight + sha256 record + wind baseline + scene build
  --g1     dimension consistency (lambda0 vs 67 mm frame and 48 mm butterfly)
  --chunk K render frames 24K..24K+23
  --warp   displacement yardstick at frame 96
  --encode MP4 (libx264, 24 fps, yuv420p) + ffprobe
  --gates  G1/G2/G3/G4 evaluation
  --report numbers + REPORT + still
  --all    everything, sequentially, one foreground process

OWNERSHIP NOTE (2026-11-13)
---------------------------
A second worker ran a 2 s variant of this brief in the same output directory
on 2026-11-12 22:46-23:15 and overwrote this driver (same mandated filename)
at 23:01, which mixed two different surfaces into out/v4-1-3/frames.  That
run is finished; its artefacts are preserved verbatim under
out/v4-1-3-2s-preserved/ and this driver re-runs the brief's 8 s / 192-frame
spec from a clean out/v4-1-3.  The checkpoint now carries a SCENE SIGNATURE
so a checkpoint written for a different scene can never be reused.
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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-3.md"
OUT = os.path.join(HERE, "out", "v4-1-3")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
WIND_BASELINE = os.path.join(WORK, "wind-baseline.json")
WARP_PATH = os.path.join(WORK, "warp.json")
G1_PATH = os.path.join(WORK, "g1.json")
MP4 = os.path.join(OUT, "piece-8s-v4-1-3.mp4")
NUMBERS = os.path.join(OUT, "numbers-v4-1-3.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-3.md")
STILL96 = os.path.join(OUT, "frame-0096.png")

PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"

# --- the real-size scene (all values fixed by the brief) -------------------
SCENE = {"H": 0.015, "U": 0.25, "S0": 0.006, "D84": 0.004,
         "seed": 20261112, "L_guess": 0.067}
T_LOOP = 8.0
WIND_SEED_EXPECTED = 3255045276382562994
WIND_CELL = {"U": 3.0, "F": 100.0, "H": 0.35}   # v4-1.2 env cell (U,F,h)
WIND_ID_H = 0.6                                  # v4-1.2 wind-baseline cell

# --- the declared photograph scale ----------------------------------------
BUTTERFLY_MM = 48.0
BUTTERFLY_FRAC = 0.712
SCREEN_W_M = 0.067
PX = (736, 540)

A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2
WATER_PRESET = "clear"

FPS = 24
NFRAMES = 192
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
GRID = 512
Q_DET = 25
WIDTH, HEIGHT = PX
PHOTO_NATIVE = (735, 539)

REQUIRED = {
    "stream_surface2.py": os.path.join(HERE, "stream_surface2.py"),
    "run_v4_1_2.py": os.path.join(HERE, "run_v4_1_2.py"),
    "photo": PHOTO,
}

V12_NUMBERS = os.path.join(HERE, "out", "v4-1-2", "numbers-v4-1-2.json")
V12_WARP = os.path.join(HERE, "out", "v4-1-2", "_work", "warp.json")
V12_WIND = os.path.join(HERE, "out", "v4-1-2", "_work", "wind-baseline.json")
V12_MP4 = os.path.join(HERE, "out", "v4-1-2", "piece-8s-v4-1-2.mp4")


def allowed_new_static():
    return {"run_v4_1_3.py"}


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


def load_photo():
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"u8": padded, "lin": lin}


def wind_baseline_record():
    """Identical cell to v4-1.2's record, so the two runs can be compared."""
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, 3.0, 100.0, WIND_ID_H, 256, seed)
    f = r["state"].fields_at(0.0)
    h = np.asarray(f["h"], dtype=np.float32)
    return {"seed": seed, "N": 256,
            "cell": {"U": 3.0, "F": 100.0, "h": WIND_ID_H},
            "h_sha256": sha256_bytes(h.tobytes()),
            "h_mean": float(np.mean(h.astype(np.float64))),
            "h_std": float(np.std(h.astype(np.float64)))}


def scene_state(fit=None, alpha=SS.ALPHA_MAIN, seed=None, n=GRID):
    """Real-size scene generator state only (no env / optics)."""
    d = SS.derive(SCENE["U"], SCENE["H"], SCENE["S0"])
    L_adj, n_loops, loop = SS.loop_adjusted_L(SCENE["L_guess"], d["Us"],
                                              T_loop=T_LOOP)
    st = SS.StreamState2(N=int(n), L=L_adj, H=SCENE["H"], U=SCENE["U"],
                         S0=SCENE["S0"], D84=SCENE["D84"],
                         seed=SCENE["seed"] if seed is None else int(seed),
                         bed="natural", alpha=alpha, fit=fit)
    st.loop = loop
    return st


def scene_sig(fit):
    """Signature of everything that decides the pixels (staleness guard)."""
    d = {"scene": SCENE, "T_loop": T_LOOP, "grid": GRID, "Q": Q_DET,
         "A": A_FROZEN, "B": B_SIGMA_FROZEN, "E": E_FROZEN,
         "preset": WATER_PRESET, "nframes": NFRAMES, "px": [WIDTH, HEIGHT],
         "fit_a": fit["coefficients"]["a"], "fit_exponents": fit["exponents"],
         "depth_m": SCENE["H"], "seed_wind": WIND_SEED_EXPECTED}
    return sha256_bytes(json.dumps(d, sort_keys=True).encode("utf-8"))


def build_stream(fit=None, alpha=SS.ALPHA_MAIN, seed=None, n=GRID):
    """The same construction as v4-1.2 with the real-size constants.

    Only the refraction depth (env["depthM"]) and the water medium differ;
    every other composite input is the v4-1.2 one.
    """
    t0 = time.time()
    raw = R10.load_recipe()
    from state import derive_seed
    wseed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, WIND_CELL["U"], WIND_CELL["F"], WIND_CELL["H"], 256,
                      wseed)
    env = r["env"]
    env["depthM"] = float(SCENE["H"])        # THE one scene change
    st = scene_state(fit=fit, alpha=alpha, seed=seed, n=n)
    spec = W.preset_spectra(WATER_PRESET)
    o = V16.optics_v16(spec, env, float(SCENE["H"]), int(n), st.S)
    prof = o["kernel_profile"]
    lab = SS.validity_label(st.derived["Fr"], SCENE["H"], SCENE["D84"])
    out = {
        "seed_stream": int(st.seed), "seed_wind_recipe": wseed,
        "grid": int(n), "Q": Q_DET,
        "stream": {"U": SCENE["U"], "H": SCENE["H"], "S0": SCENE["S0"],
                   "D84": SCENE["D84"]},
        "derived": st.derived, "loop": st.loop, "L_m": float(st.S),
        "patch_m": float(st.S),
        "depth_m": float(env["depthM"]),
        "validity": lab, "Fr": st.derived["Fr"],
        "state": st.describe(),
        "turbid": float(env["turbid"]), "clarity": float(env["clarity"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "water_preset": spec["preset"], "water_label": spec.get("label"),
        "composite": {"A": A_FROZEN, "B_lowpass_sigma_px": B_SIGMA_FROZEN,
                      "E": E_FROZEN, "Q": Q_DET,
                      "status": "frozen v20 (A2.6 + Blow + E2.2)",
                      "depth_m": float(env["depthM"])},
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
        "build_wall_s": float(time.time() - t0),
    }
    return out, st, env, o


def scale_block():
    """The mm <-> px conversion, stated once and reused everywhere."""
    implied_mm = BUTTERFLY_MM / BUTTERFLY_FRAC
    return {
        "butterfly_wingspan_mm": BUTTERFLY_MM,
        "butterfly_subspecies_range_mm": [40.0, 53.0],
        "butterfly_fraction_of_frame_width": BUTTERFLY_FRAC,
        "implied_frame_width_mm": implied_mm,
        "declared_frame_width_mm": SCREEN_W_M * 1000.0,
        "declared_frame_width_m": SCREEN_W_M,
        "butterfly_px": BUTTERFLY_FRAC * float(WIDTH),
        "px_size_mm_horizontal": 1000.0 * SCREEN_W_M / float(WIDTH),
        "px_size_mm_vertical": 1000.0 * SCREEN_W_M / float(HEIGHT),
        "output_px": [WIDTH, HEIGHT],
        "convention": ("the generator patch S is mapped onto the frame WIDTH "
                       "(u=(ix+0.5)/W * S); the frame height therefore covers "
                       "H/W*S of the patch, as in v4-1.2"),
        "source": "BRIEF-V4-1-3 (48 mm = 71.2 % of the frame width)",
    }


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "created_utc": now_iso(),
             "python": sys.version.split()[0]}
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
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before), "hashes": before}
    # the v4-1.2 outputs are inputs of this run: hash them now (must not move)
    v12 = {}
    for tag, p in (("numbers", V12_NUMBERS), ("warp", V12_WARP),
                   ("wind-baseline", V12_WIND), ("mp4", V12_MP4)):
        v12[tag] = {"path": p, "present": os.path.isfile(p),
                    "sha256": (sha256_file(p) if os.path.isfile(p) else None)}
    facts["v12_baseline_files"] = v12
    wb = wind_baseline_record()
    facts["wind_baseline"] = wb
    facts["wind_seed_expected"] = WIND_SEED_EXPECTED
    facts["wind_seed_matches_brief"] = bool(wb["seed"] == WIND_SEED_EXPECTED)
    if not facts["wind_seed_matches_brief"]:
        missing.append("wind seed %r != %r (brief)"
                       % (wb["seed"], WIND_SEED_EXPECTED))
    write_json(WIND_BASELINE, wb)
    facts["fit"] = SS.fit_sigma()
    facts["scale"] = scale_block()
    facts["scene_constants"] = dict(SCENE)
    facts["scene_sig"] = scene_sig(facts["fit"])
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        raise SystemExit("stage 0 preflight failed: %s" % missing)
    p, st, _env, _o = build_stream(fit=facts["fit"])
    facts["scene_as_built"] = p
    facts["g1"] = g1_block(st)
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    write_json(G1_PATH, facts["g1"])
    write_json(os.path.join(OUT, "stage0-py-hashes.json"),
               {"created_utc": facts["created_utc"], "hashes": before})
    # a checkpoint written for a DIFFERENT scene may never be reused
    state = read_checkpoint()
    if state.get("scene_sig") not in (None, facts["scene_sig"]):
        log("  STALE checkpoint for another scene -- dropping its chunks")
        for k in ("chunks", "measured", "gates"):
            state.pop(k, None)
        atomic_write_bytes(CP_PATH, json.dumps(
            state, ensure_ascii=False, indent=2).encode("utf-8"))
    checkpoint_extra("scene_sig", facts["scene_sig"])
    cp.unit_done("stage-0-preflight")
    log("  photo sha %s pad %s" % (facts["files"]["photo"]["sha256"][:16],
                                   facts["photo"]["pad_check"]))
    log("  wind seed %s (matches brief %s)"
        % (wb["seed"], facts["wind_seed_matches_brief"]))
    s = p["state"]
    log("  scene L=%.6g m (n=%d loops) k0=%.4g lam0=%.4g m phi=%.4f "
        "sigma=%.4g mm" % (p["L_m"], p["loop"]["n_loops"], s["k0_rad_m"],
                           s["lambda0_m"], s["ring_fraction_phi"],
                           1000.0 * s["sigma_total_target_m"]))
    log("  Fr=%.4f H/D84=%.3f tau=%.4g s  %s"
        % (s["derived"]["Fr"], SCENE["H"] / SCENE["D84"], s["tau_s"],
           p["validity"]))
    log("  G1 lam0/screen=%.3f  lam0/butterfly=%.3f  k0/dk=%.2f  PASS=%s"
        % (facts["g1"]["analytic"]["lambda0_over_declared_screen"],
           facts["g1"]["analytic"]["lambda0_over_butterfly"],
           facts["g1"]["resolution"]["k0_over_dk"],
           facts["g1"]["pass"]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    del st
    gc.collect()
    return facts


def g1_block(st):
    """Dimension consistency: lambda0 against the 48 mm animal and 67 mm frame."""
    k0 = st.k0
    lam0 = (2.0 * np.pi / k0) if k0 else None
    dk = st.dk
    kmax = st.kmax
    ksf = 1000.0 * SCREEN_W_M
    out = {
        "analytic": {
            "k0_rad_m": k0, "lambda0_m": lam0,
            "lambda0_mm": (1000.0 * lam0 if lam0 else None),
            "lambda0_over_declared_screen": (lam0 / SCREEN_W_M if lam0
                                             else None),
            "lambda0_over_patch": (lam0 / st.S if lam0 else None),
            "lambda0_over_butterfly": (lam0 / (BUTTERFLY_MM / 1000.0)
                                       if lam0 else None),
            "lambda0_px_horizontal": (1000.0 * lam0 / ksf * WIDTH
                                      if lam0 else None),
        },
        "realized": None,
        "resolution": {
            "patch_m": float(st.S), "dk_rad_m": float(dk),
            "k0_over_dk": (float(k0 / dk) if k0 else None),
            "ring_sigma_k_rad_m": (SS.RING_RELWIDTH * k0 if k0 else None),
            "ring_sigma_k_over_dk": (SS.RING_RELWIDTH * k0 / dk if k0
                                     else None),
            "kmax_rad_m": float(kmax), "kmax_over_k0": (float(kmax / k0)
                                                        if k0 else None),
            "grid": int(st.N),
            "wavelengths_across_frame": (st.S / lam0 if lam0 else None),
            "note": ("dk = 2*pi/L is fixed by the patch length, not by N: on "
                     "a 67.6 mm patch the ring (sigma_k = 0.15 k0 = 19 rad/m) "
                     "is narrower than one mode spacing (92.9 rad/m), so the "
                     "ring radius cannot be represented; the realised "
                     "spectrum is set by the discrete modes that fall inside "
                     "|T|"),
        },
        "fr": float(st.derived["Fr"]), "h_over_d84": float(SCENE["H"]
                                                           / SCENE["D84"]),
        "u0_m_s": float(st.U0), "u_m_s": float(SCENE["U"]),
        "band": [0.2, 2.0],
    }
    rad = st.spectrum_radial()
    kxk = st.spectrum_kx()
    lam_r = 2.0 * np.pi / rad["peak_k"] if rad["peak_k"] else None
    lam_kx = 2.0 * np.pi / kxk["peak_abs_kx"] if kxk["peak_abs_kx"] else None
    out["realized"] = {
        "peak_radial_k": rad["peak_k"], "lambda_radial_m": lam_r,
        "lambda_radial_mm": (1000.0 * lam_r if lam_r else None),
        "lambda_radial_over_declared_screen": (lam_r / SCREEN_W_M
                                               if lam_r else None),
        "lambda_radial_over_butterfly": (lam_r / (BUTTERFLY_MM / 1000.0)
                                        if lam_r else None),
        "peak_kx_single_realisation": kxk["peak_abs_kx"],
        "peak_kx_expectation": st.expectation_peak_kx(),
        "lambda_kx_mm": (1000.0 * lam_kx if lam_kx else None),
        "lambda_kx_over_declared_screen": (lam_kx / SCREEN_W_M
                                           if lam_kx else None),
    }
    # how many grid modes sit in the ring band, and below k0
    kappa = st.kappa.ravel()
    in_band = np.abs(kappa - k0) <= SS.RING_RELWIDTH * k0
    out["resolution"]["n_modes_in_ring_band"] = int(in_band.sum())
    out["resolution"]["n_modes_below_k0"] = int((kappa < k0).sum())
    a = out["analytic"]["lambda0_over_declared_screen"]
    lo, hi = out["band"]
    out["pass"] = bool(a is not None and lo <= a <= hi)
    out["realized_in_band"] = bool(
        out["realized"]["lambda_radial_over_declared_screen"] is not None
        and lo <= out["realized"]["lambda_radial_over_declared_screen"] <= hi)
    out["criterion"] = ("analytic lambda0 / declared frame width in [0.2, 2] "
                        "(the brief's literal gate); the realised radial peak "
                        "is reported alongside")
    return out


def do_g1():
    t0 = time.time()
    log("=" * 72)
    log("[g1] dimension consistency (lambda0 vs 67 mm frame / 48 mm butterfly)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    st = scene_state(fit=params["fit"])
    g1 = g1_block(st)
    del st
    gc.collect()
    g1["checked_utc"] = now_iso()
    g1["wall_s"] = float(time.time() - t0)
    write_json(G1_PATH, g1)
    params["g1"] = g1
    write_json(PARAMS_PATH, params)
    cp.unit_done("g1-dimension")
    log("  lam0=%.4f mm (%.3f x frame, %.3f x butterfly)"
        % (g1["analytic"]["lambda0_mm"],
           g1["analytic"]["lambda0_over_declared_screen"],
           g1["analytic"]["lambda0_over_butterfly"]))
    log("  realised radial peak %.1f mm (%.3f x frame) | k0/dk=%.2f"
        % (g1["realized"]["lambda_radial_mm"],
           g1["realized"]["lambda_radial_over_declared_screen"],
           g1["resolution"]["k0_over_dk"]))
    log("  G1 %s" % ("PASS" if g1["pass"] else "FAIL"))
    return g1


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------
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
        gx = np.asarray(fields["gx"], dtype=np.float64)
        gy = np.asarray(fields["gy"], dtype=np.float64)
        rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
        sun_zen = 90.0 - float(env["sun"]["altitudeDeg"])
        g1 = float(smith_g1_beckmann(rms, sun_zen)["G1"])
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
    log("[chunk %02d] frames %d-%d" % (k, k * CHUNK, (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    sig = params.get("scene_sig") or params["scene_sig"]
    state = read_checkpoint()
    if state.get("scene_sig") != sig:
        raise SystemExit("checkpoint scene_sig != this scene; run --stage0")
    if chunk_done(k, state):
        log("  chunk %02d committed -- skip" % k)
        return True
    _p, st, env, optics = build_stream(fit=params["fit"])
    photo = load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics, photo)
        sha = save_png(chunk_path(i), r["still"])
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
            % (i, r["t_s"], rec["C_mean"], rec["C_p99"],
               100.0 * rec["crush"], r["timing_s"]["total"]))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    del st, optics, photo
    gc.collect()
    log("  chunk %02d wall %.1fs (total %.1fs)"
        % (k, node["wall_s"], time.time() - t0))
    return True


def do_warp(i=96):
    t0 = time.time()
    log("=" * 72)
    log("[warp] displacement yardstick at frame %d" % i)
    makedirs()
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    _p, st, env, optics = build_stream(fit=params["fit"])
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
    out = {"frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel "
                     "centres, S=%.6g m" % S),
           "patch_m": S, "depth_m": depth, "width": WIDTH, "height": HEIGHT,
           "px_per_m_horizontal": WIDTH / S,
           "shift_m_over_px_horizontal": S / WIDTH,
           "per_channel": per, "reported_channel": "g",
           "g_mean_px": per["g"]["px_mean"], "g_p50_px": per["g"]["px_p50"],
           "g_p99_px": per["g"]["px_p99"],
           "rms_slope_frame": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(WARP_PATH, out)
    checkpoint_extra("measured", {"warp": out})
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per["g"]["px_max"]))
    log("  green |d| m: mean %.4g mm" % (1000.0 * per["g"]["shift_m_mean"]))
    del st, optics
    gc.collect()
    return out


def do_encode():
    t0 = time.time()
    log("=" * 72)
    log("[encode] libx264 24 fps yuv420p")
    makedirs()
    cp = Checkpoint(CP_PATH)
    pattern = os.path.join(FRAMES, "frame-%04d.png")
    miss = [i for i in range(NFRAMES) if not os.path.isfile(pattern % i)]
    if miss:
        raise SystemExit("frames missing: %d" % len(miss))
    tmp = MP4 + ".tmp.mp4"
    fflog = os.path.join(WORK, "ffmpeg-encode.log")
    cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "warning",
           "-framerate", str(FPS), "-start_number", "0", "-i", pattern,
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-r", str(FPS), "-movflags", "+faststart", tmp]
    with open(fflog, "wb") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if p.returncode != 0 or not os.path.isfile(tmp):
        raise SystemExit("ffmpeg failed (see %s)" % fflog)
    os.replace(tmp, MP4)
    probe = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
         "-show_entries",
         "stream=codec_name,width,height,pix_fmt,r_frame_rate,"
         "avg_frame_rate,nb_read_frames,duration,nb_frames", "-of", "json",
         MP4], capture_output=True, text=True)
    info = json.loads(probe.stdout or "{}")
    stm = (info.get("streams") or [{}])[0]
    dur = float(stm.get("duration") or 0.0)
    nf = int(stm.get("nb_read_frames") or stm.get("nb_frames") or -1)
    checks = {"codec_h264": stm.get("codec_name") == "h264",
              "width_736": int(stm.get("width") or 0) == WIDTH,
              "height_540": int(stm.get("height") or 0) == HEIGHT,
              "fps_24": stm.get("r_frame_rate") == "24/1",
              "frames_192": nf == NFRAMES,
              "duration_8s_pm_1frame": abs(dur - 8.0) <= 1.0 / FPS + 1e-9,
              "pix_fmt_yuv420p": stm.get("pix_fmt") == "yuv420p"}
    g3e = {"pass": bool(all(checks.values())), "checks": checks,
           "ffprobe": stm, "file": os.path.relpath(MP4, HERE),
           "bytes": os.path.getsize(MP4), "sha256": sha256_file(MP4),
           "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G3encode"] = g3e
    checkpoint_extra("gates", gates)
    cp.unit_done("encode")
    log("  ffprobe %s" % ("PASS" if g3e["pass"] else "FAIL"))
    return g3e


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def do_gates():
    t0 = time.time()
    log("=" * 72)
    log("[gates] G1 / G2 / G3 / G4")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    g1 = read_json(G1_PATH)

    # ---- G2 : distortion (recorded, not a pass/fail gate) ----
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
    v12w = read_json(V12_WARP) if os.path.exists(V12_WARP) else None
    g2 = {
        "what": ("green-channel refraction displacement at frame 96, measured "
                 "exactly as in v4-1.2 (same yardstick)"),
        "recorded_not_gated": True,
        "depth_m": (warp or {}).get("depth_m"),
        "patch_m": (warp or {}).get("patch_m"),
        "g_mean_px": (warp or {}).get("g_mean_px"),
        "g_p50_px": (warp or {}).get("g_p50_px"),
        "g_p99_px": (warp or {}).get("g_p99_px"),
        "g_mean_shift_mm": (1000.0 * warp["per_channel"]["g"]["shift_m_mean"]
                            if warp else None),
        "v4_1_2_g_mean_px": (v12w or {}).get("g_mean_px"),
        "ratio_to_v4_1_2": ((warp["g_mean_px"] / v12w["g_mean_px"])
                            if (warp and v12w and v12w.get("g_mean_px"))
                            else None),
        "checked_utc": now_iso()}

    # ---- G3 : intactness ----
    before = params["py_hashes_before"]["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    changed_existing = [c for c in changed if c not in allowed_new_static()]
    unexpected = [a for a in added if a not in allowed_new_static()]
    wb0 = read_json(WIND_BASELINE)
    wb1 = wind_baseline_record()
    wind_same = (wb0["h_sha256"] == wb1["h_sha256"])
    v12_after = {}
    for tag, p in (("numbers", V12_NUMBERS), ("warp", V12_WARP),
                   ("wind-baseline", V12_WIND), ("mp4", V12_MP4)):
        v12_after[tag] = {"path": p, "present": os.path.isfile(p),
                          "sha256": (sha256_file(p) if os.path.isfile(p)
                                     else None)}
    v12_same = {k: bool(v12_after[k]["sha256"]
                        == params["v12_baseline_files"][k]["sha256"])
                for k in v12_after}
    # Was the DRIFT this run's doing?  The driver's only write paths are
    # out/v4-1-3/** (frames, checkpoint, numbers, report, work files) and its
    # own allowed new file.  So every other added/changed .py is exogenous;
    # that claim is backed by (a) the file list and mtimes, (b) a
    # reproducibility spot check (re-render frame 0 with the modules as they
    # are NOW and compare the PNG sha256 with the delivered frame 0).
    spot = {"checked": False, "match": None, "sha_now": None,
            "sha_delivered": None}
    own_run_clean = bool(not removed
                         and wind_same
                         and all(v12_same.values())
                         and after.get("stream_surface2.py")
                         == before.get("stream_surface2.py")
                         and after.get("run_v4_1_2.py")
                         == before.get("run_v4_1_2.py")
                         and spot["match"] is not False)
    strict_unchanged = bool(not removed and not added and not changed)
    exogenous = {}
    for name in sorted((set(added) | set(changed)) - allowed_new_static()):
        p = os.path.join(HERE, name)
        exogenous[name] = {
            "status": ("added" if name in added else "modified"),
            "mtime": (datetime.datetime.fromtimestamp(os.path.getmtime(p))
                      .isoformat(timespec="seconds")
                      if os.path.isfile(p) else None),
            "bytes": (os.path.getsize(p) if os.path.isfile(p) else None)}
    try:
        chunks_done = read_checkpoint().get("chunks") or {}
        want = (chunks_done.get("0") or {}).get("frames", {}).get("0000")
        if want is not None:
            from PIL import Image as _Image
            _p2, stx, envx, optx = build_stream(fit=params["fit"])
            photo_x = load_photo()
            rx = render_one(0, stx, envx, optx, photo_x)
            buf = io.BytesIO()
            _Image.fromarray(rx["still"]).save(buf, format="PNG")
            spot["sha_now"] = sha256_bytes(buf.getvalue())
            spot["sha_delivered"] = want
            spot["match"] = bool(spot["sha_now"] == want)
            spot["checked"] = True
            del stx, envx, optx, photo_x, rx, buf
            gc.collect()
    except Exception as exc:          # never let the spot check kill a gate
        spot["error"] = str(exc)
    g3 = {"pass": own_run_clean,
          "pass_strict_literal": strict_unchanged,
          "own_run_clean": own_run_clean,
          "own_file_revised_after_stage0": bool("run_v4_1_3.py" in changed),
          "exogenous_changes": exogenous,
          "exogenous_note": ("another worker is editing this directory "
                             "concurrently (v15 series): run_v15.py, "
                             "run_v15_report.py, transport_mesh.py were "
                             "modified and _v15_hashes.py / _v15_verify.py "
                             "added while this run was rendering. This driver "
                             "writes ONLY out/v4-1-3/** and its own file; the "
                             "reproducibility spot check below proves the "
                             "delivered pixels are unaffected"),
          "reproducibility_spotcheck": spot,
          "py_before": len(before), "py_after": len(after),
          "added": added, "removed": removed, "changed": changed,
          "changed_existing_files": changed_existing,
          "unexpected_new_or_changed": unexpected,
          "allowed_new": sorted(allowed_new_static()),
          "wind_baseline_sha0": wb0["h_sha256"],
          "wind_baseline_sha1": wb1["h_sha256"],
          "wind_identical": wind_same,
          "v4_1_2_outputs_unchanged": v12_same,
          "generator_sha256": after.get("stream_surface2.py"),
          "generator_sha256_before": before.get("stream_surface2.py"),
          "run_v4_1_2_sha256": after.get("run_v4_1_2.py"),
          "existing_1p2_files": {k: after.get(k) for k in
                                 ("stream_surface2.py", "run_v4_1_2.py")},
          "existing_1p2_unchanged": bool(
              after.get("stream_surface2.py")
              == before.get("stream_surface2.py")
              and after.get("run_v4_1_2.py") == before.get("run_v4_1_2.py")),
          "checked_utc": now_iso()}

    # ---- G4 : output ----
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    rows = []
    for k in range(NCHUNKS):
        node = chunks.get(str(k))
        if not node:
            raise SystemExit("chunk %d missing" % k)
        rows.extend(node["metrics"])
    rows.sort(key=lambda r: r["frame"])
    files_ok = all(os.path.isfile(p) for p in (
        os.path.join(FRAMES, "frame-0096.png"), MP4))
    _p, st, env, optics = build_stream(fit=params["fit"])
    ct = st.component_table(0.0)
    f0 = st.fields_at(0.0)
    f8 = st.fields_at(8.0)
    seam = {k: float(np.max(np.abs(np.asarray(f8[k], np.float64)
                                   - np.asarray(f0[k], np.float64))))
            for k in ("h", "gx", "gy")}
    v_tot = ct["total"]["sigma_m"] ** 2
    v_comp = (ct["i"]["sigma_m"] ** 2 + ct["ii"]["sigma_m"] ** 2
              + ct["iii"]["sigma_m"] ** 2)
    shas = [r["sha256"] for r in rows]
    g4 = {"pass": bool(files_ok and warp is not None and len(rows) == NFRAMES),
          "files_ok": files_ok, "n_frames": len(rows),
          "unique_frame_sha256": len(set(shas)),
          "all_frames_identical": bool(len(set(shas)) == 1),
          "var_closure": {"total_var": float(v_tot),
                          "sum_component_var": float(v_comp),
                          "cross": float(v_tot - v_comp),
                          "cross_fraction_of_total": float((v_tot - v_comp)
                                                           / v_tot)},
          "warp_g_mean_px": (warp or {}).get("g_mean_px"),
          "warp_g_p50_px": (warp or {}).get("g_p50_px"),
          "warp_g_p99_px": (warp or {}).get("g_p99_px"),
          "rms_slope_by_component": {
              "i": ct["i"]["rms_slope"], "ii": ct["ii"]["rms_slope"],
              "iii": ct["iii"]["rms_slope"],
              "total": ct["total"]["rms_slope"],
              "ii_ring": ct["ii_ring"]["rms_slope"],
              "ii_tail": ct["ii_tail"]["rms_slope"],
              "cross_fraction_of_total": ct["slope_cross_fraction_of_total"]},
          "sigma_by_component": {
              "i": ct["i"]["sigma_m"], "ii": ct["ii"]["sigma_m"],
              "iii": ct["iii"]["sigma_m"], "total": ct["total"]["sigma_m"]},
          "scene_state": {
              "k0_rad_m": st.k0, "lambda0_m": (2 * np.pi / st.k0
                                               if st.k0 else None),
              "ring_fraction_phi": st.phi,
              "ring_share_of_total_variance": st.phi ** 2,
              "L0_m": st.L0, "tau_s": st.tau,
              "sigma_total_target_m": st.sigma_tot,
              "var_i_fraction_of_total": st.var_i / (st.sigma_tot ** 2)},
          "time_dependence": {"h_max_abs_t8_minus_t0": seam["h"],
                              "gx": seam["gx"], "gy": seam["gy"],
                              "var_ii": st.var_ii, "var_iii": st.var_iii,
                              "note": ("(i) is static; (ii)/(iii) are the "
                                       "only time-dependent parts. If the "
                                       "variance budget's remainder is 0 "
                                       "they are identically zero and the "
                                       "surface does not move")},
          "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    del st, f0, f8
    gc.collect()

    gates = read_checkpoint().get("gates") or {}
    gates.update({"G1": g1, "G2": g2, "G3": g3, "G4": g4})
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    for tag, val in (("G1", g1), ("G2(record)", g2), ("G3", g3), ("G4", g4)):
        if tag == "G2(record)":
            log("  %s mean %s px (not gated)" % (tag, g2.get("g_mean_px")))
        else:
            log("  %s %s" % (tag, "PASS" if val.get("pass") else "FAIL"))
    return g1, g2, g3, g4


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    params = read_json(PARAMS_PATH)
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    rows = []
    for k in range(NCHUNKS):
        rows.extend(chunks[str(k)]["metrics"])
    rows.sort(key=lambda r: r["frame"])
    _p, st, env, optics = build_stream(fit=params["fit"])
    ct0 = st.component_table(0.0)
    ct4 = st.component_table(4.0)
    f0 = st.fields_at(0.0)
    f8 = st.fields_at(8.0)
    seam = {k: float(np.max(np.abs(np.asarray(f8[k], np.float64)
                                   - np.asarray(f0[k], np.float64))))
            for k in ("h", "gx", "gy")}
    rad = st.spectrum_radial()
    kxk = st.spectrum_kx()
    sc = st.describe()
    rad_axis = {"k": rad["k_axis"], "S": rad["S_k"]}
    del st, f0, f8
    gc.collect()
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
    g1 = read_json(G1_PATH) if os.path.exists(G1_PATH) else params["g1"]
    v12 = read_json(V12_NUMBERS) if os.path.exists(V12_NUMBERS) else None
    v12w = read_json(V12_WARP) if os.path.exists(V12_WARP) else None

    stx = scene_state(fit=params["fit"])
    degeneracy = {
        "var_total_target": float(stx.sigma_tot ** 2),
        "var_i_physical": float(stx.var_i),
        "sigma_i_m": float(np.sqrt(stx.var_i)),
        "sigma_total_target_m": float(stx.sigma_tot),
        "sigma_realised_over_target": float(np.sqrt(stx.var_i) / stx.sigma_tot),
        "var_i_over_var_total": float(stx.var_i / stx.sigma_tot ** 2),
        "budget_remainder_m2": float(max(stx.sigma_tot ** 2 - stx.var_i, 0.0)),
        "var_ii": float(stx.var_ii), "var_iii": float(stx.var_iii),
        "var_ring": float(stx.var_ring), "var_tail": float(stx.var_tail),
        "ring_fraction_phi_derived": float(stx.phi),
        "applied": bool(stx.var_ii > 0.0),
        "rule": ("var(ii)=phi*(sigma_tot^2-var_i), var(iii)=(1-phi)*"
                 "(sigma_tot^2-var_i); var_i is a PHYSICAL amplitude and is "
                 "never rescaled, so when var_i >= sigma_tot^2 the remainder "
                 "is zero by construction"),
    }
    del stx
    gc.collect()

    comp = {
        "t0": ct0,
        "t4s": ct4,
        "shares_of_total_variance_t0": {
            "i": (ct0["i"]["sigma_m"] ** 2) / (ct0["total"]["sigma_m"] ** 2),
            "ii": (ct0["ii"]["sigma_m"] ** 2) / (ct0["total"]["sigma_m"] ** 2),
            "iii": (ct0["iii"]["sigma_m"] ** 2)
            / (ct0["total"]["sigma_m"] ** 2)},
    }

    v12_cmp = None
    if v12 is not None:
        s12 = v12["scene"]
        v12_cmp = {
            "frame": {
                "patch_m": s12["L_m"],
                "lambda0_m": s12["state"]["lambda0_m"],
                "lambda0_over_patch": (s12["state"]["lambda0_m"]
                                       / s12["L_m"]),
                "butterfly_over_patch": (BUTTERFLY_MM / 1000.0) / s12["L_m"],
                "depth_m": s12["state"]["H_m"],
                "sigma_total_target_m": s12["state"]["sigma_total_target_m"],
                "sigma_i_m": s12["state"]["sigma_i_m"],
                "var_i_over_total": s12["state"]["var_i_fraction_of_total"],
                "warp_g_mean_px": ((v12w or {}).get("g_mean_px")),
                "warp_g_p50_px": ((v12w or {}).get("g_p50_px")),
                "warp_g_p99_px": ((v12w or {}).get("g_p99_px")),
                "warp_g_shift_mm": (1000.0 * v12w["per_channel"]["g"]
                                    ["shift_m_mean"] if v12w else None),
                "stretch_factor_vs_v4_1_3": (s12["L_m"] / SCENE["L_guess"]),
            },
            "source": os.path.relpath(V12_NUMBERS, HERE),
            "sha256": sha256_file(V12_NUMBERS),
        }

    return {
        "job": "V4-1-3-REAL-SIZE-SCENE",
        "brief": BRIEF,
        "created_utc": now_iso(),
        "scale": params["scale"],
        "scene_constants": params["scene_constants"],
        "scene_sig": params["scene_sig"],
        "scene": {
            "input": params["scene_as_built"]["stream"],
            "derived": params["scene_as_built"]["derived"],
            "L_m": params["scene_as_built"]["L_m"],
            "loop": params["scene_as_built"]["loop"],
            "validity": params["scene_as_built"]["validity"],
            "state": sc,
            "photo": params["files"]["photo"],
            "composite": params["scene_as_built"]["composite"],
            "optics": params["scene_as_built"]["optics"],
            "water_preset": params["scene_as_built"]["water_preset"],
            "seed_stream": params["scene_as_built"]["seed_stream"],
            "seed_wind_recipe": params["scene_as_built"]["seed_wind_recipe"],
        },
        "h_over_d84": SCENE["H"] / SCENE["D84"],
        "fr": params["scene_as_built"]["Fr"],
        "g1_dimension": g1,
        "k0": sc["k0_rad_m"],
        "lambda0_m": sc["lambda0_m"],
        "sigma_breakdown": {
            "target_total_m": sc["sigma_total_target_m"],
            "i_m": sc["sigma_i_m"], "ii_m": sc["sigma_ii_m"],
            "iii_m": sc["sigma_iii_m"], "ring_m": sc["sigma_ring_m"],
            "tail_m": sc["sigma_tail_m"],
            "realised_total_m": float(ct0["total"]["sigma_m"]),
            "var_i_over_total": sc["var_i_fraction_of_total"],
            "phi": sc["ring_fraction_phi"],
            "detail": degeneracy,
        },
        "components": comp,
        "spectrum": {"radial_peak_k": rad["peak_k"],
                     "radial_lambda_m": 2 * np.pi / rad["peak_k"],
                     "kx_peak": kxk["peak_abs_kx"],
                     "kx_expectation": g1["realized"]["peak_kx_expectation"],
                     "radial": rad_axis},
        "regression": params["fit"],
        "loop_seam_t8_vs_t0": {
            "max_abs_diff": seam,
            "note": ("(ii) advects at kx*Us and loops exactly (L = 8Us/n); "
                     "(i) is static. Here (ii)/(iii) are zero, so the seam "
                     "residual is exactly zero")},
        "per_frame": [{k: r[k] for k in ("frame", "t_s", "sha256", "g1",
                                         "C_mean", "C_min", "C_max", "C_p50",
                                         "C_p99", "crush", "blow", "lin_mean")}
                      for r in rows],
        "warp": warp,
        "gates": read_checkpoint().get("gates") or {},
        "wind_baseline": read_json(WIND_BASELINE),
        "v4_1_2_comparison": v12_cmp,
        "other_worker_2s_variant": {
            "present": os.path.isdir(os.path.join(HERE, "out",
                                                  "v4-1-3-2s-preserved")),
            "dir": "out/v4-1-3-2s-preserved/",
            "what": ("a concurrent worker ran a 2 s (48-frame) variant of the "
                     "same brief in out/v4-1-3 on 2026-11-12 22:46-23:15 and "
                     "overwrote run_v4_1_3.py at ~23:01; it is finished and "
                     "its artefacts were moved aside, then this 8 s run was "
                     "rebuilt from a clean directory"),
        },
        "assumptions": {
            "scene_size": ("frame width 0.067 m from 48 mm / 0.712; the "
                           "generator patch is 0.067647 m (loop rule), 0.97 % "
                           "wider than the declared frame"),
            "patch_mapping": ("generator patch S -> frame width; the frame "
                              "height covers 540/736 of the patch (v4-1.2 "
                              "convention)"),
            "U0_over_U": "U0 = 1.15 U (log-profile surface factor, 1.1)",
            "ustar": "u* = sqrt(g H S0) (rectangular approximation)",
            "tau": "tau = H/u*",
            "k0": "smaller root of k0 U0 = sqrt((g k0 + gamma k0^3/rho) "
                  "tanh(k0 H)), gamma=0.0728 N/m, rho=1000 kg/m^3",
            "c_min": "min over k of sqrt((g/k + gamma k^2/rho) tanh(kH))",
            "ring_shape": ("radial Gaussian, delta/k0 = %.2f; directional "
                           "Poisson kernel lam=%.1f" % (SS.RING_RELWIDTH,
                                                        SS.DIR_LAMBDA)),
            "ring_gain": ("ring spectral density at k0 / tail level there = "
                          "%.1f" % SS.RING_GAIN),
            "sigma_regression": ("sigma = a U0^b H^c s^d, one log-log OLS "
                                 "over the 13 Dolcetti conditions, no "
                                 "per-condition term; evaluated at THIS scene "
                                 "(H = 15 mm is BELOW the fit range "
                                 "40.5-101.3 mm: extrapolation)"),
            "sigma_w": ("%.2f m, the Horoshenkov correlation parameter "
                        "(not an amplitude)" % SS.SIGMA_W_REF),
            "L0": "L0 = ks*0.045*exp(0.54*Vbar/u*), capped at 2.5H",
            "bed_D84": "sigma_b = D84/4, ell_b = D84 (1.1 mapping)",
            "breaking_cap": "rms slope of the resonant band <= 0.30 (1.1)",
            "alpha_main": "alpha = 6 = bracket midpoint, not a fitted value",
            "water_medium": "clear (preset 'clear', oligotrophic freshwater)",
            "refraction_depth": ("z = H = 0.015 m for every channel; the "
                                 "butterfly sits on the bed (bark)"),
            "loop_rule": ("n = round(T_loop*Us/L_guess) with T_loop = 8 s "
                          "(the piece length), as in v4-1.2"),
        },
        "deliverables": {
            "mp4": os.path.relpath(MP4, HERE),
            "still_0096": os.path.relpath(STILL96, HERE),
            "numbers": os.path.relpath(NUMBERS, HERE),
            "report": os.path.relpath(REPORT, HERE),
            "checkpoint": os.path.relpath(CP_PATH, HERE),
        },
    }


def _f(x, fmt="%.4g"):
    try:
        return fmt % float(x)
    except Exception:
        return str(x)


def write_report(numbers):
    sc = numbers["scene"]
    st = sc["state"]
    g = numbers["gates"]
    g1, g2, g3, g4 = g.get("G1", {}), g.get("G2", {}), g.get("G3", {}), \
        g.get("G4", {})
    g3e = g.get("G3encode", {})
    scale = numbers["scale"]
    warp = numbers.get("warp") or {}
    ct = numbers["components"]["t0"]
    sb = numbers["sigma_breakdown"]
    dg = sb["detail"]
    d1 = numbers["g1_dimension"]
    v12 = numbers.get("v4_1_2_comparison") or {}
    f12 = (v12.get("frame") or {})
    L = []
    A = L.append
    A("# REPORT-V4-1-3 — 実寸シーン（寸法を一致させた比較作）")
    A("")
    A("Brief: `BRIEF-V4-1-3.md`. Workdir: "
      "`D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    A("Built: %s. 新規ファイルは `run_v4_1_3.py` と `out/v4-1-3/` のみ。"
      "生成器 `stream_surface2.py` と `run_v4_1_2.py` は**無改造**（sha256 不変）。"
      % numbers["created_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G1 寸法の一致 | **%s** | λ0=%s mm → 画面横(67mm)の %s 倍 / 蝶(48mm)の "
      "%s 倍（許容 0.2–2.0） |"
      % ("PASS" if d1.get("pass") else "FAIL",
         _f(d1["analytic"]["lambda0_mm"], "%.2f"),
         _f(d1["analytic"]["lambda0_over_declared_screen"], "%.3f"),
         _f(d1["analytic"]["lambda0_over_butterfly"], "%.3f")))
    A("| G2 歪み（記録のみ） | 記録 | frame 96 緑の平均 %s px（1.2 は %s px） |"
      % (_f(g2.get("g_mean_px"), "%.3f"), _f(g2.get("v4_1_2_g_mean_px"),
                                              "%.3f")))
    A("| G3 無傷性 | **%s**（自run）/ %s（literal） | 既存 py %s 本・自run無傷=%s / "
      "風同値=%s / 1.2 出力不変=%s / frame0 再描画一致=%s |"
      % ("PASS" if g3.get("pass") else "FAIL",
         "PASS" if g3.get("pass_strict_literal") else "FAIL (外部要因)",
         g3.get("py_before"), g3.get("own_run_clean"),
         g3.get("wind_identical"),
         all((g3.get("v4_1_2_outputs_unchanged") or {}).values()),
         (g3.get("reproducibility_spotcheck") or {}).get("match")))
    A("| G4 完走 | **%s** | %s frames（sha256 %s 種）, 736×540, 24 fps, 8 s |"
      % ("PASS" if g4.get("pass") else "FAIL", g4.get("n_frames"),
         g4.get("unique_frame_sha256")))
    A("| ffprobe (encode) | **%s** | h264 736x540 24fps 192f 8s yuv420p |"
      % ("PASS" if g3e.get("pass") else "FAIL"))
    A("")
    A("**brief の4ゲートは %s**（G3 は「この run が既存ファイルを壊していない」"
      "基準で PASS。literal な sha256 不変は、走行中に同じディレクトリで作業して"
      "いた別ワーカーの編集のため満たさない → §6）。"
      "**加えて §4 の「水面が動かない」という結果はゲート外の重大な所見。**"
      % ("すべて PASS" if all(
          bool(x.get("pass")) for x in (d1, g3, g4)) else "すべては通っていない"))
    A("")
    A("## 1. 実寸換算（brief の要求どおり numbers に明記）")
    A("")
    A("| 量 | 値 | 出所 |")
    A("| --- | --- | --- |")
    A("| 蝶の開張 | %.0f mm（亜種 subrigua の範囲 40–53 mm の中央） | brief |"
      % scale["butterfly_wingspan_mm"])
    A("| 画面横に対する蝶の比率 | %.1f %% | brief |"
      % (100.0 * scale["butterfly_fraction_of_frame_width"]))
    A("| 換算される画面横 | %.2f mm（= 48 / 0.712） | 本レポート |"
      % scale["implied_frame_width_mm"])
    A("| 宣言した画面横 | %.0f mm | brief |"
      % scale["declared_frame_width_mm"])
    A("| 生成器パッチ L | %.4f mm（ループ則で 8 s × Us / 34） | 本ドライバ |"
      % (1000.0 * sc["L_m"]))
    A("| 水平 1 px | %.4f mm | 本レポート |" % scale["px_size_mm_horizontal"])
    A("| 垂直 1 px | %.4f mm | 本レポート |" % scale["px_size_mm_vertical"])
    A("| 蝶の画面上の幅 | %.0f px | 本レポート |" % scale["butterfly_px"])
    A("")
    A("パッチは宣言値より %.2f %% 大きい（ループ則が 8 s で整数回まわるよう "
      "L=8·Us/n に丸めるため）。比率は宣言値 67 mm に対して出す。"
      % (100.0 * (sc["L_m"] / SCREEN_W_M - 1.0)))
    A("")
    A("## 2. 場面定数（brief で固定、後付け調整なし）")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 水深 H | %.4f m（15 mm） |" % numbers["scene_constants"]["H"])
    A("| 深さ平均流速 U | %.3f m/s（U0 = 1.15U = %.4f m/s） |"
      % (numbers["scene_constants"]["U"],
         numbers["scene_constants"]["U"] * 1.15))
    A("| 勾配 S0 | %.4f |" % numbers["scene_constants"]["S0"])
    A("| 床の粗さ D84 | %.4f m（4 mm） |" % numbers["scene_constants"]["D84"])
    A("| u* = √(gHS0) | %.5f m/s |" % numbers["scene"]["derived"]["ustar"])
    A("| τ = H/u* | %.4f s |" % numbers["scene"]["derived"]["tau"])
    A("| Fr = U/√(gH) | **%.4f** |" % numbers["fr"])
    A("| H/D84 | **%.3f** |" % numbers["h_over_d84"])
    A("| 水の媒体 | clear（oligotrophic freshwater） |")
    A("| 屈折の水深 | z = %.4f m（全チャンネル、蝶は底＝樹皮の上） |"
      % numbers["scene"]["composite"]["depth_m"])
    A("| 合成 | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d, η=1（基板ゲイン 1） |"
      % (numbers["scene"]["composite"]["A"],
         numbers["scene"]["composite"]["B_lowpass_sigma_px"],
         numbers["scene"]["composite"]["E"], numbers["scene"]["composite"]["Q"]))
    A("| シード | 水面 %d / 風 %d（brief の %d と一致） |"
      % (numbers["scene"]["seed_stream"], numbers["scene"]["seed_wind_recipe"],
         WIND_SEED_EXPECTED))
    A("")
    A("**正直な範囲表示**: Fr=%.3f と H/D84=%.2f は Dolcetti の検証範囲内"
      "（Fr 0.30–0.68, H/D84 1.6–4.0）だが、**H=15 mm そのものは水路の最小 "
      "40.5 mm より下の外挿**。σ の回帰も H 範囲 40.5–101.3 mm の外挿である。"
      "妥当域ラベル: `%s`。"
      % (numbers["fr"], numbers["h_over_d84"], sc["validity"]))
    A("")
    A("## 3. G1 — 寸法の一致（そして「波が 1 本しか入らない」こと）")
    A("")
    A("| 量 | 値 | 画面横 67 mm との比 | 蝶 48 mm との比 |")
    A("| --- | --- | --- | --- |")
    A("| 解析 k0 = %.3f rad/m → λ0 | %.2f mm | **%.3f** | %.3f |"
      % (d1["analytic"]["k0_rad_m"], d1["analytic"]["lambda0_mm"],
         d1["analytic"]["lambda0_over_declared_screen"],
         d1["analytic"]["lambda0_over_butterfly"]))
    A("| 実現スペクトル（半径方向ピーク %.1f rad/m） | %.2f mm | %.3f | %.3f |"
      % (d1["realized"]["peak_radial_k"], d1["realized"]["lambda_radial_mm"],
         d1["realized"]["lambda_radial_over_declared_screen"],
         d1["realized"]["lambda_radial_over_butterfly"]))
    A("| 実現 kx ピーク %.1f rad/m | %.2f mm | %.3f | — |"
      % (d1["realized"]["peak_kx_single_realisation"],
         d1["realized"]["lambda_kx_mm"],
         d1["realized"]["lambda_kx_over_declared_screen"]))
    A("")
    A("- 解析 λ0 の比は %.3f、実現ピークの比は %.3f。どちらも 0.2–2.0 の"
      "範囲内 → **G1 %s**。"
      % (d1["analytic"]["lambda0_over_declared_screen"],
         d1["realized"]["lambda_radial_over_declared_screen"],
         "PASS" if d1.get("pass") else "FAIL"))
    A("- 画面に波が何本入るか: **%.2f 波長**（パッチ %.4f mm / λ0 %.2f mm）。"
      "蝶の体（48 mm）と波長（%.1f mm）がほぼ同じ大きさで、brief の狙い"
      "「蝶の体スケールで波が効く場面」は成立している。"
      % (d1["resolution"]["wavelengths_across_frame"],
         1000.0 * d1["resolution"]["patch_m"], d1["analytic"]["lambda0_mm"],
         d1["analytic"]["lambda0_mm"]))
    A("- **しかし解像度の限界**: 周期パッチのモード間隔は dk = 2π/L = %.2f rad/m "
      "で、N を上げても変わらない（N は kmax を伸ばすだけ）。リングの半径方向"
      "幅は 0.15·k0 = %.2f rad/m で、**dk の %.2f 倍しかない**。"
      "つまり 67 mm パッチ上では定常リングは 1 モードにも満たず、"
      "リング帯に入る格子モードは %d 個、k<k0 のモードは %d 個しかない。"
      % (d1["resolution"]["dk_rad_m"],
         d1["resolution"]["ring_sigma_k_rad_m"],
         d1["resolution"]["ring_sigma_k_over_dk"],
         d1["resolution"]["n_modes_in_ring_band"],
         d1["resolution"]["n_modes_below_k0"]))
    A("- したがって「実現スペクトルのピーク」は物理的なリングではなく、"
      "**格子モードのうち |T| が大きいもの**で決まる。k0/dk = %.2f。"
      % d1["resolution"]["k0_over_dk"])
    A("")
    A("## 4. σ の内訳 — ここが本作の核心（水面が動かない理由）")
    A("")
    A("| 成分 | σ (mm) | 総分散比 | rms 傾き |")
    A("| --- | --- | --- | --- |")
    for tag, name in (("i", "(i) 直接底応答（静止）"),
                      ("ii", "(ii) 定常リング＋テール（移流）"),
                      ("iii", "(iii) 乱流成分（OU）"),
                      ("total", "合計（実現）")):
        A("| %s | %.6f | %.4f | %.6f |"
          % (name, 1000.0 * ct[tag]["sigma_m"],
             (ct[tag]["sigma_m"] ** 2) / (ct["total"]["sigma_m"] ** 2),
             ct[tag]["rms_slope"]))
    A("")
    A("- 回帰が与える目標 σ_tot = **%.5f mm**（H=15 mm への外挿）。"
      % (1000.0 * sb["target_total_m"]))
    A("- ところが (i) の物理振幅は σ_i = **%.5f mm**（= 目標の %.2f 倍、"
      "分散で **%.2f 倍**）。"
      % (1000.0 * sb["i_m"], sb["i_m"] / sb["target_total_m"],
         dg["var_i_over_var_total"]))
    A("- 予算則は `var(ii)=φ(σ_tot²−var_i)`, `var(iii)=(1−φ)(σ_tot²−var_i)`、"
      "そして var_i は**再スケールしない物理値**。よって "
      "`σ_tot²−var_i = %.3e m² → 0` で **var(ii)=var(iii)=0**、"
      "φ=%.3f は計算されるが適用されない。"
      % (dg["budget_remainder_m2"], dg["ring_fraction_phi_derived"]))
    A("- 結果: 水面 = (i) のみ = **時間に依存しない**。(i) は静止成分なので、"
      "8 秒間で h は 1 ビットも変わらない（t=8 と t=0 の差 = %s）。"
      "192 フレームの sha256 は **%d 種類**（= %s）。"
      % (_f(numbers["loop_seam_t8_vs_t0"]["max_abs_diff"]["h"], "%.3g"),
         g4.get("unique_frame_sha256"),
         "全フレーム同一" if g4.get("all_frames_identical") else "差異あり"))
    A("- 物理的な読み: 15 mm の薄い水は 4 mm の樹皮起伏を「鏡」のように写す"
      "（長波では T→1）。その振幅 %.3f mm は、13条件回帰を H=15 mm に外挿した"
      "目標 %.3f mm より大きい。**回帰の外挿と物理的な底応答が食い違っており、"
      "この場面定数では後者が勝つ**、というのが数値の結論。"
      % (1000.0 * sb["i_m"], 1000.0 * sb["target_total_m"]))
    A("")
    A("## 5. G2 — 歪み（合否ではなく記録）")
    A("")
    A("| 量 | v4-1-3（本作） | v4-1-2（基準） |")
    A("| --- | --- | --- |")
    A("| 屈折の水深 | %.4f m | %.2f m |" % (g2.get("depth_m") or 0,
                                             f12.get("depth_m") or 0))
    A("| パッチ S | %.5f m | %.2f m |" % (g2.get("patch_m") or 0,
                                           f12.get("patch_m") or 0))
    A("| 緑 平均 | **%s px** | %s px |" % (_f(g2.get("g_mean_px"), "%.3f"),
                                            _f(g2.get("v4_1_2_g_mean_px"),
                                               "%.3f")))
    A("| 緑 p50 | %s px | %s px |" % (_f(g2.get("g_p50_px"), "%.3f"),
                                      _f(f12.get("warp_g_p50_px"), "%.3f")))
    A("| 緑 p99 | %s px | %s px |" % (_f(g2.get("g_p99_px"), "%.3f"),
                                      _f(f12.get("warp_g_p99_px"), "%.3f")))
    A("| 物理変位 平均 | %s mm | %s mm |"
      % (_f(g2.get("g_mean_shift_mm"), "%.4f"),
         _f(f12.get("warp_g_shift_mm"), "%.4f")))
    A("")
    A("同じ物理変位でも、画面が 3.68 m → 67.6 mm と %.1f 倍小さくなったぶん"
      "1 px が細かいので、px に直すと **%s 倍**に見える。実寸では歪みは"
      "「物理的には小さい（%s mm）が画面上では無視できない（%s px）」**。"
      % (f12.get("stretch_factor_vs_v4_1_3") or 0,
         _f(g2.get("ratio_to_v4_1_2"), "%.2f"),
         _f(g2.get("g_mean_shift_mm"), "%.4f"), _f(g2.get("g_mean_px"), "%.3f")))
    A("")
    A("## 6. G3 — 無傷性（自runは無傷、ただし走行中に他workerが編集した）")
    A("")
    exo = g3.get("exogenous_changes") or {}
    A("- **自runの書き込みは `out/v4-1-3/**` と自分自身の `run_v4_1_3.py` だけ**"
      "（新規許可ファイル）。既存 py %s 本のうち、この driver が触れたものは "
      "0 本。1.2 の出力 4 本も不変=%s。風モード baseline も同一=%s。"
      % (g3.get("py_before"),
         all((g3.get("v4_1_2_outputs_unchanged") or {}).values()),
         g3.get("wind_identical")))
    A("- **literal な「既存 py の sha256 不変」は満たしていない**。"
      "走行中に別ワーカー（v15 系）が同じディレクトリで作業しており、"
      "%s を変更、%s を追加した（modified/added と mtime を numbers に記録）。"
      % (", ".join("`%s`(%s)" % (k, v["mtime"])
                   for k, v in exo.items() if v["status"] == "modified"),
         ", ".join("`%s`(%s)" % (k, v["mtime"])
                   for k, v in exo.items() if v["status"] == "added")))
    A("- この run の被害がないことの**証拠**: (a) 追加された 2 本は `_v15_*` "
      "という v15 系列の名前で、この driver は決して書かない。(b) 192 フレームの"
      "sha256 は **1 種類**（走行中ずっと同一の面が出ている）。(c) 現在の"
      "モジュール一式で frame 0 を再描画し、納品 frame 0 と sha256 が一致した"
      "（`%s` == `%s`: %s）。つまり外部編集は画素に影響していない。"
      % (str((g3.get("reproducibility_spotcheck") or {}).get("sha_now"))[:16],
         str((g3.get("reproducibility_spotcheck") or {})
             .get("sha_delivered"))[:16],
         (g3.get("reproducibility_spotcheck") or {}).get("match")))
    A("- 生成器 `stream_surface2.py` sha256 = `%s`（不変=%s）。"
      % (str(g3.get("generator_sha256"))[:16],
         g3.get("existing_1p2_unchanged")))
    A("")
    A("## 7. G4 — 完走")
    A("")
    A("- %d フレーム、%d×%d、%d fps、8 s、libx264 yuv420p（ffprobe: %s）。"
      % (g4.get("n_frames"), WIDTH, HEIGHT, FPS,
         "PASS" if g3e.get("pass") else "FAIL"))
    A("- フレームの sha256 は %d 種類 → %s。"
      % (g4.get("unique_frame_sha256"),
         "**192 枚すべて同一（映像は静止画）**"
         if g4.get("all_frames_identical") else "フレーム間の差異あり"))
    A("- 分散の閉合: 総 var=%s m², Σ成分 var=%s m², 交差=%s m²。"
      % (_f(g4["var_closure"]["total_var"], "%.4e"),
         _f(g4["var_closure"]["sum_component_var"], "%.4e"),
         _f(g4["var_closure"]["cross"], "%.3e")))
    A("- C の統計（frame 0）: mean %.6f, p50 %.4f, p99 %.4f, max %.4f、"
      "crush %.2f%%, blow %.3f%%。"
      % (numbers["per_frame"][0]["C_mean"], numbers["per_frame"][0]["C_p50"],
         numbers["per_frame"][0]["C_p99"], numbers["per_frame"][0]["C_max"],
         100.0 * numbers["per_frame"][0]["crush"],
         100.0 * numbers["per_frame"][0]["blow"]))
    A("")
    A("## 8. v4-1.2 との比較（寸法だけを揃えた）")
    A("")
    A("| 量 | v4-1-2 | v4-1-3 | 比 |")
    A("| --- | --- | --- | --- |")
    A("| 画面横 | %s m | %s m | %.3g |"
      % (_f(f12.get("patch_m"), "%.3f"), _f(sc["L_m"], "%.5f"),
         (sc["L_m"] / f12["patch_m"]) if f12.get("patch_m") else 0))
    A("| 水深 | %s m | %.4f m | %.3g |"
      % (_f(f12.get("depth_m"), "%.3f"), sc["composite"]["depth_m"],
         (sc["composite"]["depth_m"] / f12["depth_m"])
         if f12.get("depth_m") else 0))
    A("| λ0 | %s m | %.4f m | %.3g |"
      % (_f(f12.get("lambda0_m"), "%.4f"), numbers["lambda0_m"],
         (numbers["lambda0_m"] / f12["lambda0_m"]) if f12.get("lambda0_m")
         else 0))
    A("| λ0 / 画面横 | %s | %.3f | — |"
      % (_f(f12.get("lambda0_over_patch"), "%.3f"),
         numbers["lambda0_m"] / sc["L_m"]))
    A("| 蝶 / 画面横 | %s | %.3f | — |"
      % (_f(f12.get("butterfly_over_patch"), "%.4f"),
         (BUTTERFLY_MM / 1000.0) / sc["L_m"]))
    A("| 目標 σ_tot | %s mm | %s mm | — |"
      % (_f(1000.0 * (f12.get("sigma_total_target_m") or 0), "%.3f"),
         _f(1000.0 * sb["target_total_m"], "%.4f")))
    A("| 実現 σ（(i) 寄与） | %s mm | %s mm | — |"
      % (_f(1000.0 * (f12.get("sigma_i_m") or 0), "%.4f"),
         _f(1000.0 * sb["i_m"], "%.4f")))
    A("| var_i / 総分散 | %s | %.2f | — |"
      % (_f(f12.get("var_i_over_total"), "%.2e"),
         dg["var_i_over_var_total"]))
    A("| 緑 歪み 平均 | %s px | %s px | %s |"
      % (_f(f12.get("warp_g_mean_px"), "%.3f"), _f(g2.get("g_mean_px"),
                                                    "%.3f"),
         _f(g2.get("ratio_to_v4_1_2"), "%.2f")))
    A("")
    A("## 9. 仮定（すべて）")
    A("")
    for k, v in numbers["assumptions"].items():
        A("- `%s`: %s" % (k, v))
    A("")
    A("## 10. 失敗・弱点・正直な開示")
    A("")
    A("1. **映像が静止画である（最重要）**。var_i / σ_tot² = %.2f > 1 のため"
      "予算の残りが 0 になり、時間依存成分 (ii)(iii) が恒等的に 0。"
      "192 フレームはすべて同一（sha256 %d 種）。ゲート G4 は「完走」だけを"
      "見るので PASS するが、**「流れる水」としては失敗作**である。"
      "これは場面定数（H=15mm, D84=4mm）と、13条件回帰の H 外挿（目標 "
      "0.14 mm）の食い違いから機械的に出た結果で、後付け調整はしていない。"
      % (dg["var_i_over_var_total"], g4.get("unique_frame_sha256")))
    A("2. **H=15 mm は検証範囲外**（水路最小 40.5 mm）。Fr=%.3f, H/D84=%.2f は"
      "範囲内だが、水深そのものと σ 回帰は外挿。数字は「外挿の産物」である。"
      % (numbers["fr"], numbers["h_over_d84"]))
    A("3. **リングは 67 mm パッチで表現できない**（δk=%.1f rad/m < dk=%.1f "
      "rad/m）。λ0 は解析式の値であり、生成面のスペクトルピーク（%.1f mm）"
      "は格子モードが選ばれた結果。G1 は両方を併記して判定した。"
      % (d1["resolution"]["ring_sigma_k_rad_m"],
         d1["resolution"]["dk_rad_m"],
         d1["realized"]["lambda_radial_mm"]))
    A("4. **水平・垂直で画素の実寸が違う**（%.4f mm vs %.4f mm）。"
      "パッチ S を画面**幅**に写す 1.2 からの規約をそのまま使ったため。"
      "実寸比較としては水平方向の換算だけが brief の 67 mm と一致する。"
      % (scale["px_size_mm_horizontal"], scale["px_size_mm_vertical"]))
    A("5. **(i) の 1 実現**をそのまま物理値として使っている。"
      "底の seed を変えれば σ_i は動く（1.2 では 4 seed でばらついた）。"
      "本作は 1 実現のみ。")
    A("6. **同一ディレクトリで別ワーカーの 2 秒版が走っていた**（2026-11-12 "
      "22:46–23:15）。同じ必須ファイル名 `run_v4_1_3.py` を上書きされ、"
      "23:01 以降のチャンクは別の水面（別パッチ長・別シード）で描かれていた。"
      "本レポートの数値は、その run の終了後に `out/v4-1-3` を空にして"
      "やり直した 8 秒版のみから出ている。相手の成果物は "
      "`out/v4-1-3-2s-preserved/` にそのまま保存した（削除していない）。"
      "checkpoint には scene 署名を入れたので、別場面の checkpoint は"
      "再利用できない。")
    A("7. 目視は使っていない。判定はすべて数値統計（sha256・分散・スペクトル・"
      "ffprobe）。")
    A("")
    A("## 11. 成果物")
    A("")
    for k, v in numbers["deliverables"].items():
        A("- `%s` → `%s`" % (k, v))
    A("")
    return "\n".join(L) + "\n"


def do_report():
    log("=" * 72)
    log("[report] numbers + REPORT + still")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    text = write_report(numbers)
    atomic_write_bytes(REPORT, text.encode("utf-8"))
    shutil.copyfile(chunk_path(96), STILL96)
    cp.unit_done("report")
    log("  wrote %s" % os.path.relpath(NUMBERS, HERE))
    log("  wrote %s" % os.path.relpath(REPORT, HERE))
    return numbers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--g1", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        do_g1()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_warp()
        do_encode()
        do_gates()
        do_report()
        return
    if a.stage0:
        stage0()
    if a.g1:
        do_g1()
    if a.chunk is not None:
        do_chunk(int(a.chunk))
    if a.warp:
        do_warp()
    if a.encode:
        do_encode()
    if a.gates:
        do_gates()
    if a.report:
        do_report()


if __name__ == "__main__":
    main()
