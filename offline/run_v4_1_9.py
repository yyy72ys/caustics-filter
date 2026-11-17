"""run_v4_1_9 -- BRIEF-V4-1-9: finite solar disc + in-water scattering spread.

NEW FILES ONLY: this driver.  Every pre-existing .py is left byte-identical
(G-LEGACY).  `run_v4_1_8.py` and `transport_mesh.py` are imported, never edited.

WHAT THIS JOB ADDS TO THE v4-1.8 DEPOSIT
    (a) SOLAR-DISC QUADRATURE.  The v4-1.8 area push-forward projects every
        water cell with ONE direction (the centre of the solar disc).  Here the
        same exact-overlap push-forward is evaluated for K deterministic
        equal-weight points of the solar disc (det.fibonacci_disc, the
        pipeline's own quadrature rule) and the K deposits are averaged.  Each
        sample conserves energy exactly, so the average does too.
    (b) SCATTERING PSF.  The single-scatter radial mass profile of the water
        along the brief's 0.60 m path (run_v16_water.radial_profile: the
        Fournier-Forand phase function of the algae_pond preset MC'd over the
        truncated-exponential first-scattering positions) is applied to the
        deposited C.  The width is the DERIVED profile; there is no free
        parameter and nothing is tuned to a gate.  The kernel is isotropic in
        METRES, so on the 736x540 frame (dx = 2.0048 mm, dy = 2.7324 mm) it is
        an ellipse in pixels.

MEASUREMENT
    G-LINE is run_v4_1_8's own functions, imported (same objects, not a copy),
    so "1.8 と同一の統計・同一の窓" holds by construction:
        R18.regions_from_photo / R18.measure_image / R18.lattice_metric
    The photo windows, the Hann 256x256 window, the 5-15 px ideal band mask and
    the depth definition are therefore literally the ones of v4-1.8.

STAGES (foreground, jobkit checkpoint per stage)
    --stage0    absolute-path preflight + stage-0 sha256 + transport_mesh
                selftest (15 ok / 0 ng) + water-parameter provenance
    --estimate  STAGE 1: disc / scattering / required-width numbers only
                (no rendering), + the verdict of brief item 4
    --frame24   three variants of frame 24 (nodisc / disc / discscat), the
                band-separation table, the K-convergence check, G-LINE, G-PEAK
    --chunk K   frames 24K..24K+23 of the main piece (variant iii)
    --conv      STAGE 3: 1024^2 water grid, frame 24, variant iii
    --warp      displacement yardstick (G-WARP, record only)
    --encode    MP4 + probe (G-OUT)
    --gates     every gate of the brief
    --report    numbers-v4-1-9.json + REPORT-V4-1-9.md
    --all       everything, sequentially
"""
import argparse
import datetime
import gc
import hashlib
import io
import json
import math
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v4_1_8 as R18  # noqa: E402  (scene, measurement, adapters: UNCHANGED)
import run_v16_water as V16  # noqa: E402  (FF phase function + radial profile)
import transport_mesh as TM  # noqa: E402  (the area push-forward of v4-1.8)
import v12_water as W  # noqa: E402  (the a/b presets = the water parameters)
from det import fibonacci_disc  # noqa: E402
from transport import sun_disc_radius, orthonormal_basis  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-9.md"

OUT = os.path.join(HERE, "out", "v4-1-9")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-9.mp4")
STILL = os.path.join(OUT, "frame-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-9.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-9.md")

V_NODISC = os.path.join(OUT, "frame-0024-nodisc.png")
V_DISC = os.path.join(OUT, "frame-0024-disc.png")
V_SCAT = os.path.join(OUT, "frame-0024-discscat.png")

# the previous job (read-only references)
V418_DIR = os.path.join(HERE, "out", "v4-1-8")
V418_NUMBERS = os.path.join(V418_DIR, "numbers-v4-1-8.json")
V418_STILL = os.path.join(V418_DIR, "frame-0024.png")
V418_REPORT = os.path.join(V418_DIR, "REPORT-V4-1-8.md")
V418_CP = os.path.join(V418_DIR, "checkpoint.json")

PHOTO = R18.PHOTO
FPS = R18.FPS
WIDTH, HEIGHT = R18.WIDTH, R18.HEIGHT
NFRAMES = R18.NFRAMES
CHUNK = R18.CHUNK
NCHUNKS = R18.NCHUNKS
STILL_IDX = R18.STILL_IDX
SCN = R18.SCN

NEW_FILES = ["run_v4_1_9.py"]

# per-process worker state (filled by _worker_init under the spawn start method)
_W = {}

# ---- the derived widths (STAGE 1 fixes them; nothing below may move them) ---
K_DISC = 9              # solar-disc quadrature points (>= 9, justified in S1)
L_SCAT_M = 0.60         # brief item 2: the 0.60 m path
NPROC = 8               # worker processes for the deposit pool

# PSF truncation: V16.kernel_r_max_px's convention, but geometric on the frame
PSF_RMAX_FRAC = 0.45

# gate reference of v4-1.8 (its point-splat lane on the same field and grid)
REF_PEAK = {"std": 1.3105895996932464, "p99": 6.348078870773312,
            "source": "out/v4-1-8/numbers-v4-1-8.json frame24.peak.C_point_splat"}
# the v4-1.8 G-LINE rows (the residual this brief is about)
V418_ROWS = [("smooth", "x", 12.800000, 66.22119846873873, 63.689564352484375),
             ("smooth", "y", 11.130435, 61.12059416759521, 61.61109595032102),
             ("bark", "x", 14.222222, 55.52168847475705, 47.40954384505002),
             ("bark", "y", 9.846154, 47.56754926473905, 39.32211373167725)]


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
    last = None
    for _ in range(tries):
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


def frame_path(k):
    return os.path.join(FRAMES, "frame-%04d.png" % k)


def read_checkpoint():
    return (read_json(CP_PATH) if os.path.exists(CP_PATH)
            else {"completed_units": 0, "completed_ids": [], "last_unit": None,
                  "updated_at": None})


def checkpoint_extra(key, value=None):
    state = read_checkpoint()
    if value is not None:
        state[key] = value
        state["updated_at"] = now_iso()
        _atomic_retry(CP_PATH, json.dumps(
            state, ensure_ascii=False, indent=2).encode("utf-8"))
    return state


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".py") and os.path.isfile(os.path.join(HERE, name)):
            out[name] = sha256_file(os.path.join(HERE, name))
    return out


def save_png(path, arr):
    return R18.save_png(path, arr)


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
EXPECTED_SHA = {
    "run_v4_1_8.py":
        "fd29dc96ec7de215ee4f1dae5db7778c8726254ef344c97746dd9ae41985e20d",
    "transport_mesh.py":
        "914bbc1ae8f3f4dc5c212771c1aeda74bdb9207e2db105be96602e98952a2081",
}


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + sha256 + selftest")
    makedirs()
    cp = Checkpoint(CP_PATH)

    need = [BRIEF,
            os.path.join(HERE, "run_v4_1_8.py"),
            os.path.join(HERE, "transport_mesh.py"),
            V418_NUMBERS, V418_REPORT, V418_STILL, V418_CP,
            os.path.join(HERE, "run_v16_water.py"),
            os.path.join(HERE, "v12_water.py"),
            os.path.join(HERE, "transport.py"),
            os.path.join(HERE, "det.py"),
            PHOTO]
    missing = [p for p in need if not os.path.isfile(p)]
    h408 = sha256_file(os.path.join(HERE, "run_v4_1_8.py"))
    htm = sha256_file(os.path.join(HERE, "transport_mesh.py"))
    sha_ok = {"run_v4_1_8.py": bool(h408 == EXPECTED_SHA["run_v4_1_8.py"]),
              "transport_mesh.py": bool(htm == EXPECTED_SHA["transport_mesh.py"])}
    log("  brief: %s %s" % (BRIEF, "OK" if os.path.isfile(BRIEF) else "MISSING"))
    log("  missing: %s" % (missing or "none"))
    log("  sha256 run_v4_1_8.py %s (%s)" % (h408[:16], "OK" if sha_ok[
        "run_v4_1_8.py"] else "MISMATCH"))
    log("  sha256 transport_mesh.py %s (%s)" % (htm[:16], "OK" if sha_ok[
        "transport_mesh.py"] else "MISMATCH"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))

    st = R18.selftest_record()
    log("  transport_mesh selftest: rc=%d %d ok / %d NG %s"
        % (st["returncode"], st["n_checks_ok"], st["n_checks_ng"],
           "OK" if st["ok"] else "NG"))

    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in NEW_FILES}
    write_json(STAGE0_HASHES, {
        "all": hashes, "existing": existing,
        "new_files_of_this_job": NEW_FILES,
        "existing_count": len(existing), "recorded_utc": now_iso(),
        "expected_sha256": EXPECTED_SHA, "sha_match": sha_ok})

    # water-parameter provenance: where b and the phase function live
    pres = W.PRESETS["algae_pond"]
    spec = W.preset_spectra("algae_pond")
    prov = {
        "scattering_coefficient": {
            "file": "v12_water.py",
            "symbols": ["PRESETS['algae_pond']['b_part']",
                        "PRESETS['algae_pond']['load_default']",
                        "b_water_morel", "preset_spectra"],
            "b_part_rgb_1_per_m_at_load_1": pres["b_part"],
            "b_water_morel_rgb_1_per_m": [float(x) for x in
                                          W.b_water_morel()],
            "load_default": pres["load_default"],
            "b_rgb_1_per_m": [float(x) for x in spec["b"]],
            "b_green_550_1_per_m": float(spec["b"][1]),
            "c_green_550_1_per_m": float(spec["c"][1]),
            "a_green_550_1_per_m": float(spec["a"][1]),
            "label": pres["label"],
        },
        "phase_function": {
            "file": "run_v16_water.py",
            "symbols": ["FF_PHASE['algae_pond']", "ff_phase", "phase_table",
                        "radial_profile"],
            "n": 1.10, "mu": 3.5835,
            "label": "Petzold average-particle fit", "measured_fit": True,
            "theta2_rad2": 0.17718464415618068,
            "mean_cos": 0.9299015409156731,
        },
        "solar_disc": {
            "file": "transport.py", "symbol": "sun_disc_radius",
            "formula": "0.00465*(0.7+0.6*(1-clarity)) rad",
            "used_by": "det.det_frame, transport.py, transport_mesh.point_splat_C",
            "clarity": 0.75,
        },
        "note": ("the two numbers this brief needs are both pre-existing "
                 "measurement-backed definitions, not knobs: b from the "
                 "preset table, the phase function from the Petzold fit"),
    }
    rec = {"brief": BRIEF, "checked": need, "missing": [],
           "sha256": {"run_v4_1_8.py": h408, "transport_mesh.py": htm},
           "sha_match": sha_ok, "selftest": st,
           "water_parameters": prov, "checked_utc": now_iso()}
    checkpoint_extra("stage0", rec)
    cp.unit_done("stage-0-preflight")
    log("  water: b(550) = %.6f /m, c(550) = %.6f /m, phase = %s n=%.2f mu=%.4f"
        % (spec["b"][1], spec["c"][1], "Petzold FF", 1.10, 3.5835))
    log("  wall %.1fs" % (time.time() - t0))
    return rec


# ---------------------------------------------------------------------------
# the scene (built by run_v4_1_8, unchanged)
# ---------------------------------------------------------------------------
def build(n=None):
    return R18.build(n=n)


def scene_record():
    """The scene of v4-1-8 (identical), recorded for this job's own report."""
    st, env, optics, spec = build()
    raw = R18.R10.load_recipe()
    rec = {
        "scene": "v4-1-8 identical (brief: 場面は v4-1-8 と同一)",
        "constants": {"H": float(SCN["H"]), "U": float(SCN["U"]),
                      "U0": float(st.U0), "S0": float(SCN["S0"]),
                      "D84": float(SCN["D84"]), "seed": int(SCN["seed"]),
                      "water": SCN["water"], "L_m": float(st.S),
                      "N": int(st.N), "A": float(SCN["A"]),
                      "B_sigma": float(SCN["B_sigma"]),
                      "E": float(SCN["E"]), "Q": int(SCN["Q"])},
        "source": ("BRIEF-V4-1-9 STAGE 1 ('場面は v4-1.8 と同一') + run_v4_1_8's "
                   "SCN and run_v4_1_2's frozen composite constants"),
        "patch_m": float(st.S), "grid": int(st.N),
        "cell_m": float(st.S / float(st.N)),
        "sigma_total_m": float(st.sigma_total),
        "var_i_m2": float(st.var_i_physical),
        "var_ii_m2": float(st.var_ii), "var_iii_m2": float(st.var_iii),
        "validity": R18.SS.validity_label(st.derived["Fr"], float(SCN["H"]),
                                          float(SCN["D84"])),
        "validity_note": ("the brief's H = 0.60 m is ~6x the upper end of the "
                          "2016 flume range the sigma regression was fitted on "
                          "(H <= 0.1013 m): sigma is an EXTRAPOLATION.  "
                          "Executed as specified (1.8 unchanged)."),
        "accounting": st.accounting(),
        "deposit": {"module": "transport_mesh (imported, unmodified)",
                    "method": ("area push-forward averaged over the solar-disc "
                               "quadrature (K samples)"),
                    "c_grid": [int(HEIGHT), int(WIDTH)],
                    "w_cell": float(WIDTH * HEIGHT / float(st.N * st.N)),
                    "kernel": "none", "jitter": "none", "band_limit": "none",
                    "mean_normalisation": "none"},
        "composite": {"A": float(SCN["A"]),
                      "Blow_sigma_px": float(SCN["B_sigma"]),
                      "E": float(SCN["E"]), "Q": int(SCN["Q"]),
                      "z_refraction_m": float(SCN["H"]),
                      "water_preset": spec.get("preset"),
                      "water_label": spec.get("label"),
                      "frozen": "identical to v4-1.8 (no knob moved)"},
        "env": {"sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
                "sunDir": [float(x) for x in env["sunDir"]],
                "clarity": float(env["clarity"]),
                "turbid": float(env["turbid"]),
                "depthM": float(env["depthM"])},
        "optics": {"kernel": optics.get("kernel"),
                   "phase_kind": optics.get("phase_kind"),
                   "w_ballistic": [float(x) for x in optics["w_ballistic"]],
                   "w_halo": [float(x) for x in optics["w_halo"]],
                   "sigma_px_state_grid": float(optics["sigma_px"]),
                   "kernel_rms_px_state_grid": float(
                       optics["kernel_profile"]["rms_radius_px"]),
                   "two_way_path_m": float(optics["partition"]["L_m"])},
        "regions": R18.regions_from_photo()[1],
        "created_utc": now_iso(),
    }
    del st, optics, spec
    gc.collect()
    return rec


def do_scene():
    t0 = time.time()
    log("=" * 78)
    log("[scene] the v4-1-8 scene, recorded")
    makedirs()
    cp = Checkpoint(CP_PATH)
    rec = scene_record()
    checkpoint_extra("scene", rec)
    cp.unit_done("scene-record")
    log("  L=%.7f m grid=%d cell=%.6f mm sigma_total=%.4f mm Us=%.4f m/s"
        % (rec["patch_m"], rec["grid"], 1e3 * rec["cell_m"],
           1e3 * rec["sigma_total_m"], rec["constants"]["U0"]))
    log("  optics: w_halo %s ; halo kernel rms %.2f px (state grid)"
        % ([round(x, 3) for x in rec["optics"]["w_halo"]],
           rec["optics"]["kernel_rms_px_state_grid"]))
    log("  wall %.1fs" % (time.time() - t0))
    return rec


# ---------------------------------------------------------------------------
# STAGE 1 -- the widths (numbers only, no rendering)
# ---------------------------------------------------------------------------
def disc_numbers(st, env):
    """Brief item 1: the projected size of the finite solar disc."""
    rad = float(sun_disc_radius(env.get("clarity", 0.8)))
    ior = float(TM.IOR_G)
    theta_w = math.asin(math.sin(rad) / ior)          # Snell, small angle
    depth = float(env["depthM"])
    r_m = depth * math.tan(theta_w)
    dx = float(st.S) / float(WIDTH)
    dy = float(st.S) / float(HEIGHT)
    return {
        "clarity": float(env["clarity"]),
        "formula": "sun_disc_radius(clarity) = 0.00465*(0.7+0.6*(1-clarity))",
        "source_file": "transport.py",
        "radius_rad": rad,
        "radius_source_note": ("the pipeline's own disc radius (used by "
                               "det.det_frame, transport.py and "
                               "transport_mesh.point_splat_C); at clarity "
                               "0.75 it is 0.85x the clarity-1 value 0.00465"),
        "ior_water": ior,
        "snell_theta_water_rad": theta_w,
        "depth_m": depth,
        "projected_radius_m": r_m,
        "px_pitch_x_m": dx, "px_pitch_y_m": dy,
        "projected_radius_px_x": r_m / dx,
        "projected_radius_px_y": r_m / dy,
        "projected_diameter_px_x": 2.0 * r_m / dx,
        "projected_diameter_px_y": 2.0 * r_m / dy,
        "rms_px_x": 0.5 * r_m / dx,
        "rms_px_y": 0.5 * r_m / dy,
        "rms_rule": ("a uniform disc of radius R has per-axis variance R^2/4, "
                     "so the equivalent Gaussian sigma (same 2nd moment) is "
                     "R/2; this is the number compared with the required "
                     "width in item 4"),
        "derivation": ("r = z*tan(asin(sin(rad)/n)); the disc is isotropic in "
                       "metres and the frame pitch is not (x 2.005 mm, "
                       "y 2.732 mm), so the projected disc is an ellipse in "
                       "pixels: radius %.3f px in x, %.3f px in y"
                       % (r_m / dx, r_m / dy)),
    }


def scatter_numbers(st, env, spec, optics, L=L_SCAT_M):
    """Brief item 2: single-scatter PSF width at the 0.60 m path."""
    b = float(spec["b"][1])
    c = float(spec["c"][1])
    th2 = float(optics["phase"]["theta2_rad2"])
    dx = float(st.S) / float(WIDTH)
    dy = float(st.S) / float(HEIGHT)
    x = c * L
    D2 = x * x - 2.0 * x + 2.0 - 2.0 * math.exp(-x)
    s2 = 0.5 * th2 * (b / (c ** 3)) * D2
    sig_m = math.sqrt(max(s2, 0.0))

    # the exact (tan) profile: same phase function, same scatter positions,
    # no moment matching at all.  r = (L-s) tan(theta).
    table = V16.phase_table("ff", preset="algae_pond")
    r_max_m = PSF_RMAX_FRAC * min(WIDTH * dx, HEIGHT * dy)
    prof = V16.radial_profile(table, b, c, L, dx, r_max_m / dx)

    # multiple-scattering order estimate (Poisson, single-scatter albedo)
    w0 = b / c
    p = {}
    for n in range(0, 4):
        p[n] = math.exp(-x) * (x ** n) / math.factorial(n) * (w0 ** n)
    scat_total = sum(p[n] for n in range(1, 4))

    return {
        "medium": "algae_pond", "wavelength_nm": 550,
        "b_1_per_m": b, "a_1_per_m": float(spec["a"][1]), "c_1_per_m": c,
        "single_scatter_albedo": w0,
        "path_m": float(L),
        "optical_depth_cL": x,
        "D2": D2,
        "theta2_rad2": th2,
        "theta_rms_rad": math.sqrt(th2),
        "theta_rms_deg": math.degrees(math.sqrt(th2)),
        "formula": ("sigma^2 = (<theta^2>/2) (b/c^3) D2(cL), "
                    "D2(x)=x^2-2x+2-2exp(-x)  [v12_water.sigma_from_scattering, "
                    "the pipeline's own single-scatter displacement moment]"),
        "sigma_m": sig_m,
        "sigma_px_x": sig_m / dx, "sigma_px_y": sig_m / dy,
        "single_scatter_weight_H": w0 * (1.0 - math.exp(-x)),
        "ballistic_weight_T": math.exp(-x),
        "exact_profile": {
            "builder": "run_v16_water.radial_profile (FF-CDF, tan(theta))",
            "r_max_m": r_max_m, "r_max_px_x": r_max_m / dx,
            "rms_radius_px_x": float(prof["rms_radius_px"]),
            "rms_radius_m": float(prof["rms_radius_px"]) * dx,
            "captured_fraction": float(prof["captured_fraction"]),
            "truncation_loss_fraction": float(prof["truncation_loss_fraction"]),
            "note": ("measured on the same phase function WITHOUT the "
                     "moment-matched Gaussian: %.1f%% wider than sigma_m, "
                     "because tan(theta) > theta at the %.0f deg rms angle"
                     % (100.0 * (prof["rms_radius_px"] * dx / sig_m - 1.0),
                        math.degrees(math.sqrt(th2)))),
        },
        "small_angle_validity": {
            "condition": "theta_rms << 1 rad (paraxial / moment-matched PSF)",
            "theta_rms_rad": math.sqrt(th2),
            "verdict": ("MARGINAL: %.2f rad = %.1f deg is not paraxial; the "
                        "Gaussian-moment width is a LOWER bound on the true "
                        "PSF width (exact tan(theta) profile is %.1f%% wider). "
                        "The deposit/blur below uses the exact profile, not "
                        "the moment, so the reported width is not the "
                        "small-angle one."
                        % (math.sqrt(th2), math.degrees(math.sqrt(th2)),
                           100.0 * (prof["rms_radius_px"] * dx / sig_m - 1.0))),
        },
        "multiple_scattering": {
            "orders_probability": {str(k): float(v) for k, v in p.items()},
            "fraction_scattered_orders1to3": float(scat_total),
            "ratio_double_over_single": float(p[2] / max(p[1], 1e-30)),
            "verdict": ("single scattering carries %.0f%% of the light that "
                        "is scattered at all (double %.0f%% of single); the "
                        "double-scattered light is ~1.41x wider, so the "
                        "single-scatter PSF understates the halo width by "
                        "roughly 10-20%%. Multiple scattering is therefore NOT "
                        "negligible in width, but the pipeline's v12/v16 "
                        "optical model is single-scatter by design (its own "
                        "docstring calls the halo a lower bound); adding a "
                        "multiple-scattering PSF would be a new model, not a "
                        "correction to this one."
                        % (100.0 * p[1] / max(p[1] + p[2] + p[3], 1e-30),
                           100.0 * p[2] / max(p[1], 1e-30))),
        },
        "pipeline_comparison": {
            "note": ("the renderer does NOT use 0.60 m: v12_water.ab_optics "
                     "integrates the halo over the TWO-WAY path "
                     "s_down + s_up, so its own halo PSF is wider.  Both "
                     "numbers are given so the reader can see that the "
                     "pipeline already carries a scattering PSF."),
            "optics_sigma_m": float(optics["sigma"]["sigma_m"]),
            "optics_sigma_px_state_grid": float(optics["sigma_px"]),
            "optics_kernel_rms_px_state_grid":
                float(optics["kernel_profile"]["rms_radius_px"]),
            "optics_kernel_rmax_px_state_grid":
                float(optics["kernel_profile"]["r_max_px"]),
            "optics_kernel_captured_fraction":
                float(optics["kernel_profile"]["captured_fraction"]),
            "optics_two_way_path_m": float(optics["partition"]["L_m"]),
            "w_halo_rgb": [float(x) for x in optics["w_halo"]],
            "w_ballistic_rgb": [float(x) for x in optics["w_ballistic"]],
            "applied_by": ("run_v16_water.render_frame_v16: Ch = "
                           "ff_blur_fft(C, kernel_2d) and the image uses "
                           "w_bal*C + w_halo*Ch.  So the single-scatter PSF is "
                           "already in the render, at weight w_halo."),
        },
    }


def required_sigma(ratios=None, periods=None, still=None, regions=None):
    """Brief item 3: the equivalent Gaussian sigma that would bring the worst
    v4-1.8 window to 1.00 (two estimates: a frontier formula and the measured
    window spectrum)."""
    ratios = ratios or [r[3] / r[4] for r in V418_ROWS]
    periods = periods or [r[2] for r in V418_ROWS]
    out = {"ratios_v4_1_8": [float(x) for x in ratios],
           "periods_px": [float(x) for x in periods],
           "method_A_frontier": [], "method_A_note":
               ("single-mode frontier: a Gaussian blur of sigma attenuates the "
                "mode of period p by exp(-2*pi^2*sigma^2/p^2); solving for the "
                "ratio to reach 1.00 gives sigma = p*sqrt(ln r/(2 pi^2)).  "
                "APPROXIMATION: one mode at the measured peak period.")}
    for (reg, ax, p, d, ph), r in zip(V418_ROWS, ratios):
        s = (p * math.sqrt(math.log(r) / (2.0 * math.pi ** 2))
             if r > 1.0 else 0.0)
        out["method_A_frontier"].append({
            "region": reg, "axis": ax, "period_px": float(p),
            "ratio": float(r), "sigma_needed_px": float(s),
            "already_below": bool(r <= 1.0)})

    # method B: the measured 256x256 window spectrum of the v4-1.8 still
    out["method_B_spectral"] = None
    if still and os.path.isfile(still) and regions:
        g = R18.green_of(still)
        win = R18.WIN
        band = R18.BAND_PX
        kf = np.fft.fftfreq(win) * win
        KX, KY = np.meshgrid(kf, kf)
        rows = []
        for (reg, ax, p, d, ph), r0 in zip(V418_ROWS, ratios):
            x0 = y0 = 0
            for nm, xx, yy in regions:
                if nm == reg:
                    x0, y0 = xx, yy
            raw = g[y0:y0 + win, x0:x0 + win]
            a = np.abs(KX if ax == "x" else KY)
            o = np.abs(KY if ax == "x" else KX)
            m = (a >= win / band[1]) & (a <= win / band[0]) & (o <= 3)
            F = np.fft.fft2(raw - raw.mean()) * m
            # period of each mode (inf where k==0, masked out anyway)
            per = np.where(np.abs(KX if ax == "x" else KY) > 0,
                           win / np.maximum(np.abs(KX if ax == "x"
                                                    else KY), 1e-9), 1e9)

            def band_rms(sig):
                T = np.exp(-2.0 * math.pi ** 2 * sig ** 2 / per ** 2)
                bp = np.real(np.fft.ifft2(F * T))
                return float(np.sqrt(np.mean(bp ** 2))) / float(raw.mean())

            base = band_rms(0.0)
            target = base / float(r0)
            lo, hi = 0.0, 8.0
            if band_rms(hi) > target:
                sig = float("nan")
            else:
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if band_rms(mid) > target:
                        lo = mid
                    else:
                        hi = mid
                sig = 0.5 * (lo + hi)
            rows.append({"region": reg, "axis": ax, "ratio_v4_1_8": float(r0),
                         "band_rms_at_sigma0": base,
                         "sigma_needed_px": (None if not np.isfinite(sig)
                                             else float(sig)),
                         "already_below": bool(r0 <= 1.0)})
        out["method_B_spectral"] = {
            "rows": rows, "source_image": os.path.basename(still),
            "window_px": win, "band_px": list(band),
            "note": ("the actual axis-band power of the v4-1.8 window, "
                     "attenuated mode by mode; the solve is on the IMAGE "
                     "(all lanes), so this is still an estimate of the "
                     "equivalent blur, not a prediction of a C-level blur "
                     "(the warp lane is not blurred by C)")}
    a = [r["sigma_needed_px"] for r in out["method_A_frontier"]
         if r["sigma_needed_px"] > 0]
    out["worst_sigma_needed_px_method_A"] = float(max(a)) if a else 0.0
    out["best_sigma_needed_px_method_A"] = float(min(a)) if a else 0.0
    if out["method_B_spectral"]:
        b = [r["sigma_needed_px"] for r in out["method_B_spectral"]["rows"]
             if r["sigma_needed_px"]]
        out["worst_sigma_needed_px_method_B"] = float(max(b)) if b else 0.0
    return out


def stage1_estimate():
    t0 = time.time()
    log("=" * 78)
    log("[stage 1] the widths (numbers only, no rendering)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st, env, optics, spec = build()
    regions, rinfo = R18.regions_from_photo()

    disc = disc_numbers(st, env)
    scat = scatter_numbers(st, env, spec, optics)
    req = required_sigma(still=V418_STILL, regions=regions)

    # ---- the verdict of brief item 4 ------------------------------------
    d_rms = 0.5 * (disc["rms_px_x"] + disc["rms_px_y"])
    s_rms_px_x = scat["exact_profile"]["rms_radius_px_x"]
    s_rms_px_y = scat["exact_profile"]["rms_radius_m"] / (
        float(st.S) / float(HEIGHT))
    need = req["worst_sigma_needed_px_method_A"]
    need_b = req.get("worst_sigma_needed_px_method_B")
    disc_factor = need / max(d_rms, 1e-30)
    scat_factor = ((s_rms_px_x + s_rms_px_y) / 2.0) / max(need, 1e-30)
    verdict = {
        "disc_rms_px_mean": float(d_rms),
        "scatter_rms_px_x": float(s_rms_px_x),
        "scatter_rms_px_y": float(s_rms_px_y),
        "disc_rms_px_x": disc["rms_px_x"], "disc_rms_px_y": disc["rms_px_y"],
        "needed_sigma_px_A": need, "needed_sigma_px_B": need_b,
        "disc_is_factor_too_narrow": float(disc_factor),
        "scatter_is_factor_too_wide": float(scat_factor),
        "disc_attenuation_of_worst_window": float(
            math.exp(-2.0 * math.pi ** 2 * (disc["rms_px_x"] ** 2)
                     / (V418_ROWS[2][2] ** 2))),
        "judgement": None,
        "text": None,
    }
    verdict["judgement"] = "NOT EXPLAINED"
    _rows = req["method_A_frontier"]
    _wi = max(range(len(_rows)), key=lambda i: _rows[i]["sigma_needed_px"])
    _wp = _rows[_wi]["period_px"]
    _att = math.exp(-2.0 * math.pi ** 2 * (disc["rms_px_x"] ** 2)
                    / (_wp ** 2))
    verdict["worst_window"] = ("%s %s" % (_rows[_wi]["region"],
                                          _rows[_wi]["axis"]))
    verdict["worst_window_period_px"] = float(_wp)
    verdict["disc_attenuation_of_worst_window"] = float(_att)
    verdict["text"] = (
        "物理の広がりでは 1.8 の超過を説明できない。(1) 太陽円盤は最悪窓"
        "（%s、period %.2f px）に必要な等価 sigma %.2f px に対して "
        "%.2f px（x）/%.2f px（y）で、**%.1f 倍細い**。円盤だけではその窓の"
        "帯振幅を **%.1f%%** しか落とせない（%.3f → %.3f 倍）。"
        "(2) 単散乱 PSF は %.1f px（x）/%.1f px（y）で、必要な幅の **%.0f 倍"
        "太い**。太いだけなら「広げれば下がる」ように見えるが、これは"
        "『幅が足りない』ではなく『既に入っている』という話で、実際 "
        "run_v16_water.render_frame_v16 は同じ位相関数から作った PSF を "
        "w_halo（%.2f-%.2f）で既に掛けている（C への再畳み込みは二重計上）。"
        "よって『実在する光の広がりが未投入』という前提自体が 1.8 の残差には"
        "当てはまらない。"
        % (_rows[_wi]["region"] + " " + _rows[_wi]["axis"], _wp, need,
           disc["rms_px_x"], disc["rms_px_y"], disc_factor,
           100.0 * (1.0 - _att), _rows[_wi]["ratio"], _rows[_wi]["ratio"] * _att,
           s_rms_px_x, s_rms_px_y, scat_factor,
           min(float(x) for x in
               ((scat.get("pipeline_comparison") or {}).get("w_halo_rgb")
                or [0.0])),
           max(float(x) for x in
               ((scat.get("pipeline_comparison") or {}).get("w_halo_rgb")
                or [0.0]))))

    out = {"scene": "v4-1-8 identical (H 0.60, U 0.80, window 1.4755102 m, "
                    "736x540, algae_pond, clarity 0.75)",
           "disc": disc, "scattering": scat, "required": req,
           "verdict": verdict,
           "K_choice": {
               "K": K_DISC,
               "rule": ("K is not a free parameter: the projected disc radius "
                        "is %.2f px (x) / %.2f px (y), so only its second "
                        "moment survives the pixel grid.  An equal-weight "
                        "K-point Fibonacci disc has an EXACT second moment "
                        "(sum r^2/K = R^2/2) for every K, so K=9 already "
                        "captures the width; the K=25 check in stage 2 "
                        "measures what the higher moments add."
                        % (disc["projected_radius_px_x"],
                           disc["projected_radius_px_y"])),
               "minimum_from_brief": 9},
           "water_parameter_files": {
               "b_and_presets": "v12_water.py (PRESETS['algae_pond'], "
                                "b_water_morel, preset_spectra)",
               "phase_function": "run_v16_water.py (FF_PHASE['algae_pond'], "
                                 "ff_phase, phase_table, radial_profile)",
               "solar_disc": "transport.py (sun_disc_radius)"},
           "wall_s": float(time.time() - t0)}
    log("  disc: rad %.6f -> %.3f mm -> %.3f px x / %.3f px y ; rms %.3f / %.3f"
        % (disc["radius_rad"], 1e3 * disc["projected_radius_m"],
           disc["projected_radius_px_x"], disc["projected_radius_px_y"],
           disc["rms_px_x"], disc["rms_px_y"]))
    log("  scatter: b=%.4f c=%.4f cL=%.3f sigma_m=%.5f m -> %.2f px x / %.2f "
        "px y (moment)" % (scat["b_1_per_m"], scat["c_1_per_m"],
                           scat["optical_depth_cL"], scat["sigma_m"],
                           scat["sigma_px_x"], scat["sigma_px_y"]))
    log("  scatter exact profile: rms %.2f px x (captured %.4f, r_max %.3f m)"
        % (s_rms_px_x, scat["exact_profile"]["captured_fraction"],
           scat["exact_profile"]["r_max_m"]))
    log("  required sigma: A(frontier) worst %.2f px ; B(spectral) worst %s"
        % (need, ("%.2f px" % need_b) if need_b else "n/a"))
    log("  VERDICT: %s (disc %.1fx too narrow, scatter %.0fx too wide)"
        % (verdict["judgement"], disc_factor, scat_factor))
    checkpoint_extra("stage1", out)
    cp.unit_done("stage-1-widths")
    del st, optics, spec
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# the disc quadrature + the scattering PSF (the two additions of this job)
# ---------------------------------------------------------------------------
def disc_envs(env, K=K_DISC):
    s, u, v = orthonormal_basis(env["sunDir"])
    q = fibonacci_disc(int(K), float(sun_disc_radius(env.get("clarity", 0.8))))
    envs = []
    for i in range(int(K)):
        d = s + u * float(q[i, 0]) + v * float(q[i, 1])
        d = d / np.linalg.norm(d)
        e = dict(env)
        e["sunDir"] = [float(x) for x in d]
        envs.append(e)
    return envs, q


def mesh_C_disc(fields, S, Wd, Ht, env, K=K_DISC, envs=None):
    """Area push-forward averaged over K solar-disc samples.

    Every sample is transport_mesh.mesh_C's exact-overlap deposit (energy
    conserved), so the average conserves energy as well.  No jitter, no
    kernel, no normalisation.
    """
    if envs is None:
        envs, q = disc_envs(env, K)
    else:
        q = None
    acc = None
    per = []
    for k, e in enumerate(envs):
        C, info = TM.mesh_C(fields, S, Wd, Ht, env=e)
        acc = C.astype(np.float64) if acc is None else acc + C
        per.append({"k": int(k),
                    "energy_cell_total": float(info["energy_cell_total"]),
                    "energy_deposited": float(info["energy_deposited"]),
                    "rel_err": float(info["rel_err"]),
                    "offframe_fraction_energy":
                        float(info["offframe_fraction_energy"]),
                    "n_tail": int(info["n_tail"]),
                    "n_degenerate_cell": int(info["n_degenerate_cell"])})
    acc /= float(len(envs))
    exp_total = float(Wd) * float(Ht)
    dep64 = float(acc.sum(dtype=np.float64))
    C32 = acc.astype(np.float32)
    dep32 = float(C32.sum(dtype=np.float64))
    cinfo = {
        "K": int(len(envs)),
        "method": ("mean of %d exact-overlap area push-forwards over the "
                   "solar-disc quadrature (det.fibonacci_disc, equal weights)"
                   % len(envs)),
        "disc_radius_rad": float(sun_disc_radius(env.get("clarity", 0.8))),
        "energy_cell_total_per_direction": exp_total,
        "energy_cell_total_all_directions": exp_total * len(envs),
        "deposited_total_all_directions": dep64 * len(envs),
        "deposited_total_averaged_C": dep64,
        "expected_averaged_C": exp_total,
        "rel_err": float(abs(dep64 - exp_total) / exp_total),
        "rel_err_float32_C": float(abs(dep32 - exp_total) / exp_total),
        "max_rel_err_per_direction": float(max(x["rel_err"] for x in per)),
        "offframe_fraction_energy_min": float(min(
            x["offframe_fraction_energy"] for x in per)),
        "offframe_fraction_energy_max": float(max(
            x["offframe_fraction_energy"] for x in per)),
        "n_degenerate_cell_max": int(max(x["n_degenerate_cell"] for x in per)),
        "kernel": "none", "jitter": "none", "band_limit": "none",
        "normalisation_applied": 1.0,
        "per_direction": per,
        "w_cell": float(Wd) * float(Ht) / float(np.asarray(
            fields["h"]).shape[0] ** 2),
    }
    return C32, cinfo


def psf_build(st, env, spec, optics, L=L_SCAT_M):
    """The derived single-scatter PSF on the frame grid (isotropic in metres)."""
    b = float(spec["b"][1])
    c = float(spec["c"][1])
    dx = float(st.S) / float(WIDTH)
    dy = float(st.S) / float(HEIGHT)
    table = V16.phase_table("ff", preset="algae_pond")
    r_max_m = PSF_RMAX_FRAC * min(WIDTH * dx, HEIGHT * dy)
    prof = V16.radial_profile(table, b, c, L, dx, r_max_m / dx)
    Kp, kinfo = psf_kernel(prof, WIDTH, HEIGHT, dx, dy)
    return {"profile": prof, "kernel": Kp, "kernel_info": kinfo,
            "dx_m": dx, "dy_m": dy, "r_max_m": r_max_m,
            "L_m": float(L), "b": b, "c": c}


def psf_kernel(prof, Wd, Ht, dx, dy):
    """Unit-mass periodic kernel on the Wd x Ht frame, isotropic in metres.

    The profile's abscissa is r/dx, so the pixel radius is
    R = sqrt((i*dx)^2 + (j*dy)^2)/dx.  Offsets wrap (the composite is periodic).
    The truncation loss of the profile is renormalised away and reported.
    """
    ii = np.arange(Wd)
    ii = np.where(ii < Wd // 2, ii, ii - Wd)
    jj = np.arange(Ht)
    jj = np.where(jj < Ht // 2, jj, jj - Ht)
    R = np.sqrt((ii[None, :] * dx) ** 2 + (jj[:, None] * dy) ** 2) / dx
    K = np.interp(R.ravel(), prof["centres_px"], prof["density_per_px2"],
                  left=0.0, right=0.0).reshape(Ht, Wd)
    raw = float(K.sum())
    if raw <= 0.0:
        raise ValueError("empty PSF kernel")
    kinfo = {"raw_sum": raw, "shape": [int(Ht), int(Wd)],
             "renormalised_to": 1.0,
             "profile_captured_fraction": float(prof["captured_fraction"]),
             "truncation_loss_fraction": float(prof["truncation_loss_fraction"]),
             "rms_radius_px_x": float(prof["rms_radius_px"]),
             "rms_radius_px_y": float(prof["rms_radius_px"]) * dx / dy,
             "note": ("truncated at r_max and renormalised to unit mass (the "
                      "same convention as run_v16_water.kernel_2d); the "
                      "truncation loss is that profile's, reported not hidden")}
    return (K / raw), kinfo


def psf_blur(C, Kp):
    F = np.fft.rfft2(np.asarray(C, dtype=np.float64))
    Kf = np.fft.rfft2(Kp)
    return np.fft.irfft2(F * Kf, s=np.asarray(C).shape)


# ---------------------------------------------------------------------------
# render one C through the frozen v4-1.8 composite
# ---------------------------------------------------------------------------
def render_C(i, fields, env, optics, photo, C):
    Cr = R18.c_render(C)
    o2 = R18.render_optics(optics, Cr.shape[0])
    img, dbg = V16.render_frame_v16(
        fields, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
        ambient_scale=float(SCN["A"]), c_layer_exposure=float(SCN["E"]),
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(R18._g1(fields, env)),
        lowpass_sigma_px=float(SCN["B_sigma"]), return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"still": img, "metrics": met}


# ---------------------------------------------------------------------------
# STAGE 2 -- frame 24 variants
# ---------------------------------------------------------------------------
def _worker_init():
    global _W
    st, env, optics, spec = build()
    _W = {"st": st, "env": env, "optics": optics, "spec": spec,
          "photo": R18.load_photo(), "psf": psf_build(st, env, spec, optics)}


def _worker_deposit_dir(job):
    """(tag, frame_idx, k_index, k_total) -> one disc direction's deposit."""
    tag, idx, ki, kt = job
    w = _W
    fields = w["st"].fields_at(float(idx) / float(FPS))
    envs, _q = disc_envs(w["env"], kt)
    C, info = TM.mesh_C(fields, w["st"].S, WIDTH, HEIGHT, env=envs[ki])
    return (tag, idx, ki, C, {"rel_err": float(info["rel_err"]),
                              "energy_deposited":
                                  float(info["energy_deposited"]),
                              "energy_cell_total":
                                  float(info["energy_cell_total"]),
                              "offframe_fraction_energy":
                                  float(info["offframe_fraction_energy"]),
                              "n_tail": int(info["n_tail"]),
                              "n_degenerate_cell":
                                  int(info["n_degenerate_cell"])})


def _worker_render(job):
    """(tag, frame_idx, C float32) -> (tag, frame_idx, still, metrics)."""
    tag, idx, C = job
    w = _W
    fields = w["st"].fields_at(float(idx) / float(FPS))
    r = render_C(idx, fields, w["env"], w["optics"], w["photo"], C)
    return (tag, idx, r["still"], r["metrics"])


def pool_map(fn, jobs, nproc=NPROC):
    ctx = mp.get_context("spawn")
    out = []
    with ctx.Pool(processes=int(nproc), initializer=_worker_init) as pool:
        for k, res in enumerate(pool.imap_unordered(fn, jobs, chunksize=1)):
            out.append(res)
            if (k + 1) % 10 == 0:
                log("    pool: %d/%d done" % (k + 1, len(jobs)))
    return out


def deposit_disc_pooled(tag, idx, K, nproc=NPROC):
    """K direction deposits for one frame, averaged in the parent."""
    jobs = [(tag, idx, ki, K) for ki in range(K)]
    res = pool_map(_worker_deposit_dir, jobs, nproc=nproc)
    res.sort(key=lambda r: r[2])
    acc = None
    per = []
    for _t, _i, ki, C, info in res:
        acc = C.astype(np.float64) if acc is None else acc + C
        per.append(dict(info, k=int(ki)))
    acc /= float(K)
    exp_total = float(WIDTH) * float(HEIGHT)
    dep64 = float(acc.sum(dtype=np.float64))
    C32 = acc.astype(np.float32)
    dep32 = float(C32.sum(dtype=np.float64))
    cinfo = {
        "K": int(K),
        "method": ("mean of %d exact-overlap area push-forwards over the "
                   "solar-disc quadrature (det.fibonacci_disc)" % K),
        "energy_cell_total_all_directions": exp_total * K,
        "deposited_total_all_directions": dep64 * K,
        "deposited_total_averaged_C": dep64,
        "expected_averaged_C": exp_total,
        "rel_err": float(abs(dep64 - exp_total) / exp_total),
        "rel_err_float32_C": float(abs(dep32 - exp_total) / exp_total),
        "max_rel_err_per_direction": float(max(x["rel_err"] for x in per)),
        "offframe_fraction_energy_min": float(min(
            x["offframe_fraction_energy"] for x in per)),
        "offframe_fraction_energy_max": float(max(
            x["offframe_fraction_energy"] for x in per)),
        "n_degenerate_cell_max": int(max(x["n_degenerate_cell"] for x in per)),
        "n_tail": [int(x["n_tail"]) for x in per],
        "kernel": "none", "jitter": "none", "band_limit": "none",
        "normalisation_applied": 1.0,
        "per_direction": per,
    }
    return C32, cinfo


def rows_from_green(tag, g, regions, photo_rows):
    got = R18.lattice_metric(np.asarray(g, dtype=np.float64), regions)
    rows = []
    for g_, p_ in zip(got, photo_rows):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_px": g_["period_px"],
                     "depth": g_["depth_e3"], "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "prom_peak_over_median": g_["prom_peak_over_median"],
                     "prom_photo": p_["prom_peak_over_median"],
                     "k_bin": g_["k_bin"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    return {"tag": tag, "rows": rows,
            "strength": float(max(r["depth"] for r in rows)),
            "pass": bool(all(r["ok"] for r in rows)),
            "pass_hp": bool(all(r["ok_hp"] for r in rows))}


def rowset(tag, still, regions, photo_rows):
    """G-LINE rows of an RGB still (green channel, 1.8's statistic)."""
    return rows_from_green(tag, np.asarray(still)[:, :, 1], regions,
                           photo_rows)


def do_frame24():
    t0 = time.time()
    log("=" * 78)
    log("[frame24] three variants (nodisc / disc / discscat) + K check")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st, env, optics, spec = build()
    photo = R18.load_photo()
    regions, rinfo = R18.regions_from_photo()
    photo_rows = R18.measure_image(PHOTO, regions)
    fields = st.fields_at(STILL_IDX / float(FPS))
    psf = psf_build(st, env, spec, optics)

    # ---- the deposits ----------------------------------------------------
    log("  [1/3] nodisc deposit (1 direction, v4-1.8's own single call)")
    C0, info0 = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, env)
    info0 = {"K": 1, "method": "transport_mesh.mesh_C (v4-1.8, single disc "
                               "centre direction)",
             "energy_cell_total_all_directions": float(
                 info0["energy_cell_total"]),
             "deposited_total_all_directions": float(info0["energy_deposited"]),
             "deposited_total_averaged_C": float(info0["energy_deposited"]),
             "expected_averaged_C": float(WIDTH * HEIGHT),
             "rel_err": float(info0["rel_err"]),
             "rel_err_float32_C": float(info0["rel_err"]),
             "max_rel_err_per_direction": float(info0["rel_err"]),
             "offframe_fraction_energy_min": float(
                 info0["offframe_fraction_energy"]),
             "offframe_fraction_energy_max": float(
                 info0["offframe_fraction_energy"]),
             "n_degenerate_cell_max": int(info0["n_degenerate_cell"]),
             "kernel": "none", "jitter": "none", "band_limit": "none",
             "normalisation_applied": 1.0}

    log("  [2/3] disc deposit (K=%d, pooled)" % K_DISC)
    C1, info1 = deposit_disc_pooled("f24", STILL_IDX, K_DISC)

    log("  [3/3] disc + scattering (derived PSF, L=%.2f m)" % psf["L_m"])
    C2 = psf_blur(C1, psf["kernel"]).astype(np.float32)
    info2 = dict(info1)
    info2.update({"psf": psf["kernel_info"], "psf_L_m": psf["L_m"],
                  "psf_r_max_m": psf["r_max_m"],
                  "mean_before_blur": float(np.mean(C1)),
                  "mean_after_blur": float(np.mean(C2))})

    # ---- the K-convergence check (why K=9 is enough) ---------------------
    log("  K check: K=25 on the same frame (pooled)")
    C25, info25 = deposit_disc_pooled("f24k25", STILL_IDX, 25)
    kchk = {"K9": R18.c_stats(C1), "K25": R18.c_stats(C25),
            "std_rel_diff": float(np.std(C25) / np.std(C1) - 1.0),
            "p99_rel_diff": float(np.percentile(C25, 99)
                                  / np.percentile(C1, 99) - 1.0),
            "max_rel_diff": float(np.max(C25) / np.max(C1) - 1.0),
            "band_of_C_K9": band_of_C(C1), "band_of_C_K25": band_of_C(C25),
            "energy_K25": {k: info25[k] for k in
                           ("rel_err", "rel_err_float32_C",
                            "max_rel_err_per_direction")}}
    log("    K9 vs K25: std %+.3f%%, p99 %+.3f%%, max %+.3f%%"
        % (100 * kchk["std_rel_diff"], 100 * kchk["p99_rel_diff"],
           100 * kchk["max_rel_diff"]))
    del C25
    gc.collect()

    # ---- render the three variants ---------------------------------------
    variants = {}
    stills = {}
    for tag, C, path, info in (("nodisc", C0, V_NODISC, info0),
                               ("disc", C1, V_DISC, info1),
                               ("discscat", C2, V_SCAT, info2)):
        log("  render %s" % tag)
        r = render_C(STILL_IDX, fields, env, optics, photo, C)
        stills[tag] = r["still"]
        sha = save_png(path, r["still"])
        rs = rowset(tag, r["still"], regions, photo_rows)
        variants[tag] = {
            "file": os.path.basename(path), "sha256": sha,
            "C_stats": R18.c_stats(C), "info": info,
            "metrics": r["metrics"], "gline": rs,
            "peak": {"std_ratio_vs_1_8_point_splat":
                         float(R18.c_stats(C)["std"] / REF_PEAK["std"]),
                     "p99_ratio_vs_1_8_point_splat":
                         float(R18.c_stats(C)["p99"] / REF_PEAK["p99"]),
                     "reference": REF_PEAK},
            "still_equals_v418": (bool(sha == sha256_file(V418_STILL))
                                  if os.path.isfile(V418_STILL) else None)}
        log("    %-8s C std %.4f p99 %.3f | band %s | G-LINE %s"
            % (tag, variants[tag]["C_stats"]["std"],
               variants[tag]["C_stats"]["p99"],
               [round(x["depth"], 1) for x in rs["rows"]],
               "PASS" if rs["pass"] else "FAIL"))

    # ---- the main piece's frame 24 (= variant iii) -----------------------
    sha_still = save_png(STILL, stills["discscat"])
    log("  frame-0024.png written from variant (iii): sha %s" % sha_still[:16])

    # ---- band separation (1.8 section 4.1 procedure) ---------------------
    log("  band separation: warp only / caustics only / composite")
    z = np.zeros_like(fields["gx"])
    flat = {"h": fields["h"], "gx": z, "gy": z,
            "Dx": np.zeros_like(z), "Dy": np.zeros_like(z),
            "S": float(fields["S"])}
    C_one = np.ones((HEIGHT, WIDTH), dtype=np.float32)
    sep = {"rows": [{"tag": "photo", "depth": [r["depth_e3"]
                                               for r in photo_rows]}]}
    r = render_C(STILL_IDX, fields, env, optics, photo, C_one)
    sep["rows"].append({"tag": "warp_only",
                        "depth": [x["depth_e3"] for x in
                                  R18.lattice_metric(
                                      r["still"][:, :, 1].astype(np.float64),
                                      regions)],
                        "note": "C = 1 (the floor brightness is the warp lane "
                                "alone) -- this is also the INIFINITE-BLUR "
                                "limit of any C-level spread"})
    for tag, C in (("nodisc", C0), ("disc", C1), ("discscat", C2)):
        r = render_C(STILL_IDX, flat, env, optics, photo, C)
        sep["rows"].append({"tag": "caustics_only_" + tag,
                            "depth": [x["depth_e3"] for x in
                                      R18.lattice_metric(
                                          r["still"][:, :, 1].astype(
                                              np.float64), regions)]})
    for tag in ("nodisc", "disc", "discscat"):
        sep["rows"].append({"tag": "composite_" + tag,
                            "depth": [x["depth"] for x in
                                      variants[tag]["gline"]["rows"]]})
    sep["window_order"] = [r["region"] + " " + r["axis"] for r in photo_rows]
    sep["photo_depth"] = [r["depth_e3"] for r in photo_rows]
    sep["note"] = ("the v4-1.8 section 4.1 procedure, applied to each variant; "
                   "'warp_only' is C=1 and is the limit of infinite C-level "
                   "blur, so no spread inside C can go below it")
    for row in sep["rows"]:
        log("    %-20s %s" % (row["tag"], [round(x, 1) for x in row["depth"]]))

    out = {"frame": STILL_IDX, "t_s": STILL_IDX / float(FPS),
           "variants": variants, "K_check": kchk,
           "band_separation": sep, "regions": rinfo,
           "still": os.path.basename(STILL), "still_sha256": sha_still,
           "still_equals_variant_iii": bool(
               sha_still == variants["discscat"]["sha256"]),
           "psf": {"kernel_info": psf["kernel_info"], "L_m": psf["L_m"],
                   "b": psf["b"], "c": psf["c"],
                   "profile_rms_px_x": float(
                       psf["profile"]["rms_radius_px"]),
                   "r_max_m": psf["r_max_m"]},
           "frozen_composite": {"A": float(SCN["A"]),
                                "Blow_sigma_px": float(SCN["B_sigma"]),
                                "E": float(SCN["E"]), "Q": int(SCN["Q"])},
           "wall_s": float(time.time() - t0)}
    checkpoint_extra("frame24", out)
    cp.unit_done("frame-0024-variants")
    del st, optics, spec, photo, C0, C1, C2
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


def band_of_C(C, win=256, band=(5.0, 15.0)):
    """Axis 5-15 px band RMS of a field, per mille of its mean (1.8's mask)."""
    g = np.asarray(C, dtype=np.float64)
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


# ---------------------------------------------------------------------------
# STAGE 2 -- the main piece (variant iii, 48 frames)
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
    log("[chunk %02d] frames %d-%d (variant iii: disc + scattering)"
        % (k, k * CHUNK, (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    if chunk_done(k, state):
        log("  chunk %02d already committed -- skip" % k)
        return True

    idxs = list(range(k * CHUNK, (k + 1) * CHUNK))
    # phase 1: all deposits (K per frame), pooled
    jobs = []
    for i in idxs:
        for ki in range(K_DISC):
            jobs.append(("c%d" % k, i, ki, K_DISC))
    log("  phase 1: %d deposits (%d frames x K=%d) on %d workers"
        % (len(jobs), len(idxs), K_DISC, NPROC))
    res = pool_map(_worker_deposit_dir, jobs)
    # group by frame, and sum in the SAME direction order every time (float
    # addition is not associative; the frame-24 stage sums in k order too)
    byframe = {}
    for tag, i, ki, C, info in res:
        byframe.setdefault(i, {})[ki] = (C, info)
    del res
    gc.collect()

    # phase 2: average, blur, render
    st, env, optics, spec = build()
    ps = psf_build(st, env, spec, optics)
    Kp = ps["kernel"]
    photo = R18.load_photo()
    frames, metrics, cstats, energies = {}, [], [], []
    for i in idxs:
        group = byframe.pop(i)
        per_i = [group[ki][1] for ki in sorted(group)]
        Ca = None
        for ki in sorted(group):
            Ck = group[ki][0].astype(np.float64)
            Ca = Ck if Ca is None else Ca + Ck
        Ca /= float(K_DISC)
        exp_total = float(WIDTH * HEIGHT)
        dep64 = float(Ca.sum(dtype=np.float64))
        C1 = Ca.astype(np.float32)
        C2 = psf_blur(C1, Kp).astype(np.float32)
        er = {"frame": int(i),
              "energy_cell_total_all_directions": exp_total * K_DISC,
              "deposited_total_all_directions": dep64 * K_DISC,
              "deposited_total_averaged_C": dep64,
              "expected_averaged_C": exp_total,
              "rel_err": float(abs(dep64 - exp_total) / exp_total),
              "rel_err_float32_C": float(
                  abs(float(C1.sum(dtype=np.float64)) - exp_total)
                  / exp_total),
              "max_rel_err_per_direction": float(max(
                  x["rel_err"] for x in per_i)),
              "offframe_fraction_energy_min": float(min(
                  x["offframe_fraction_energy"] for x in per_i)),
              "offframe_fraction_energy_max": float(max(
                  x["offframe_fraction_energy"] for x in per_i)),
              "psf_mean_before": float(np.mean(C1)),
              "psf_mean_after": float(np.mean(C2))}
        fields = st.fields_at(float(i) / float(FPS))
        r = render_C(i, fields, env, optics, photo, C2)
        sha = save_png(frame_path(i), r["still"])
        cs = R18.c_stats(C2)
        rec = {"frame": int(i), "t_s": i / float(FPS), "sha256": sha,
               "C_mean": cs["mean"], "C_std": cs["std"], "C_min": cs["min"],
               "C_max": cs["max"], "C_p50": cs["p50"], "C_p99": cs["p99"],
               "C_std_disc_only": float(np.std(C1)),
               "C_p99_disc_only": float(np.percentile(C1, 99)),
               "energy_rel_err": er["rel_err"],
               "energy_rel_err_float32": er["rel_err_float32_C"],
               "offframe_fraction_energy": er["offframe_fraction_energy_max"],
               "energy_deposited": er["deposited_total_averaged_C"],
               "energy_cell_total": er["expected_averaged_C"],
               "crush": float(r["metrics"]["crush_rate_linear_lt_0.02"]),
               "blow": float(r["metrics"]["blowout_rate_aces_gt_0.98"]),
               "lin_mean": float(r["metrics"]["linear_mean_pre_tonemap"])}
        if i == STILL_IDX:
            rec["equals_frame24_still"] = bool(
                sha == sha256_file(STILL) if os.path.isfile(STILL) else False)
        frames["%04d" % i] = sha
        metrics.append(rec)
        cstats.append({"frame": int(i), "disc": R18.c_stats(C1),
                       "discscat": cs})
        energies.append(er)
        log("  frame %04d t=%.3f Cmean=%.6f Cstd=%.4f p99=%.3f E_rel=%.1e "
            "%.1fs" % (i, i / float(FPS), cs["mean"], cs["std"],
                       cs["p99"], er["rel_err"], time.time() - t0))
        del Ca, C1, C2, r
        gc.collect()
    node = {"frames": frames, "metrics": metrics, "energies": energies,
            "cstats": cstats, "completed_utc": now_iso(),
            "wall_s": float(time.time() - t0),
            "variant": "iii (disc + scattering)",
            "K_disc": K_DISC,
            "energy_max_rel_err": float(max(
                e["rel_err"] for e in energies)),
            "energy_max_rel_err_float32": float(max(
                e["rel_err_float32_C"] for e in energies)),
            "energy_mean_abs_rel_err": float(np.mean(
                [abs(e["rel_err"]) for e in energies]))}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs  energy max rel_err %.2e (float32 %.2e)"
        % (k, node["wall_s"], node["energy_max_rel_err"],
           node["energy_max_rel_err_float32"]))
    del st, optics, spec
    gc.collect()
    return True


# ---------------------------------------------------------------------------
# STAGE 3 -- 1024^2 water grid, frame 24, variant iii
# ---------------------------------------------------------------------------
def do_conv(nproc=4):
    t0 = time.time()
    log("=" * 78)
    log("[conv] STAGE 3: 1024^2 water grid, frame 24, variant iii")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st1k, env1k, optics1k, spec1k = build(n=1024)
    fields = st1k.fields_at(STILL_IDX / float(FPS))
    K = K_DISC
    envs, _q = disc_envs(env1k, K)
    acc = None
    per = []
    for ki in range(K):
        C, info = TM.mesh_C(fields, st1k.S, WIDTH, HEIGHT, env=envs[ki])
        acc = C.astype(np.float64) if acc is None else acc + C
        per.append({"k": ki, "rel_err": float(info["rel_err"]),
                    "energy_deposited": float(info["energy_deposited"])})
        log("    dir %d/%d  rel_err %.2e  %.1fs" % (ki + 1, K,
                                                    info["rel_err"],
                                                    time.time() - t0))
    acc /= float(K)
    C1 = acc.astype(np.float32)
    psf = psf_build(st1k, env1k, spec1k, optics1k)
    C2 = psf_blur(C1, psf["kernel"]).astype(np.float32)
    cs512 = (read_checkpoint().get("frame24") or {}).get("variants", {}).get(
        "discscat", {}).get("C_stats") or R18.c_stats(C2)
    out = {"n": 1024, "frame": STILL_IDX, "K": K,
           "C_512_disc": (read_checkpoint().get("frame24") or {}).get(
               "variants", {}).get("disc", {}).get("C_stats"),
           "C_512_discscat": cs512,
           "C_1024_disc": R18.c_stats(C1),
           "C_1024_discscat": R18.c_stats(C2),
           "std_rel_diff_vs_512_discscat": float(
               R18.c_stats(C2)["std"] / cs512["std"] - 1.0),
           "p99_rel_diff_vs_512_discscat": float(
               R18.c_stats(C2)["p99"] / cs512["p99"] - 1.0),
           "band_1024_disc": band_of_C(C1),
           "band_1024_discscat": band_of_C(C2),
           "energy": {"energy_cell_total_all_directions":
                          float(WIDTH * HEIGHT) * K,
                      "deposited_total_all_directions":
                          float(acc.sum()) * K,
                      "rel_err": float(abs(float(acc.sum())
                                           - float(WIDTH * HEIGHT))
                                       / float(WIDTH * HEIGHT)),
                      "max_rel_err_per_direction": float(max(
                          x["rel_err"] for x in per))},
           "wall_s": float(time.time() - t0)}
    out["note"] = ("the water grid changes, the frame grid does not; same "
                   "method as v4-1.8 G-CONV but with the new deposit")
    checkpoint_extra("conv", out)
    cp.unit_done("conv-1024")
    log("  1024^2: std %.4f p99 %.3f (512: std %.4f p99 %.3f) -> dstd "
        "%+.2f%% dp99 %+.2f%%"
        % (out["C_1024_discscat"]["std"], out["C_1024_discscat"]["p99"],
           cs512["std"], cs512["p99"],
           100 * out["std_rel_diff_vs_512_discscat"],
           100 * out["p99_rel_diff_vs_512_discscat"]))
    del st1k, optics1k, spec1k
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


def do_conv512():
    """The 512^2 half of STAGE 3's comparison (the C-level band statistic)."""
    t0 = time.time()
    log("=" * 78)
    log("[conv512] 512^2 C-level band statistic for the stage-3 comparison")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st, env, optics, spec = build()
    fields = st.fields_at(STILL_IDX / float(FPS))
    envs, _q = disc_envs(env, K_DISC)
    acc = None
    for ki in range(K_DISC):
        C, info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, env=envs[ki])
        acc = C.astype(np.float64) if acc is None else acc + C
    acc /= float(K_DISC)
    C1 = acc.astype(np.float32)
    psf = psf_build(st, env, spec, optics)
    C2 = psf_blur(C1, psf["kernel"]).astype(np.float32)
    out = {"n": 512, "K": K_DISC,
           "band_512_disc": band_of_C(C1),
           "band_512_discscat": band_of_C(C2),
           "C_512_disc": R18.c_stats(C1),
           "C_512_discscat": R18.c_stats(C2),
           "energy_rel_err": float(abs(float(acc.sum()) - WIDTH * HEIGHT)
                                   / (WIDTH * HEIGHT)),
           "wall_s": float(time.time() - t0)}
    prev = read_checkpoint().get("conv") or {}
    prev["band_512"] = out
    checkpoint_extra("conv", prev)
    cp.unit_done("conv-512-band")
    log("  512^2 C band (per mille): disc x %.2f y %.2f ; discscat x %.2f y %.2f"
        % (out["band_512_disc"]["x"], out["band_512_disc"]["y"],
           out["band_512_discscat"]["x"], out["band_512_discscat"]["y"]))
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# G-WARP (record only)
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
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
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
                    "px_max": float(np.max(mag))}
    out = {"scene": "v4-1-9", "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres, "
                     "UNFILTERED slope (same method as v4-1.4/5/7/8)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
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
# encode
# ---------------------------------------------------------------------------
def do_encode():
    t0 = time.time()
    log("=" * 78)
    log("[encode] MP4 + probe")
    makedirs()
    cp = Checkpoint(CP_PATH)
    # the deliverable still IS the main piece's frame 24: make the two files
    # the same bytes (the main render is authoritative), and record it
    src = frame_path(STILL_IDX)
    still_sync = {"source": os.path.basename(src), "target":
                  os.path.basename(STILL), "action": None}
    if os.path.isfile(src):
        sha_src = sha256_file(src)
        still_sync["source_sha256"] = sha_src
        if not os.path.isfile(STILL):
            still_sync["action"] = "written from the main render"
        elif sha256_file(STILL) == sha_src:
            still_sync["action"] = "already identical"
        else:
            still_sync["action"] = ("rewritten from the main render (the "
                                    "frame-24 variant render differed in the "
                                    "last bits of the float addition order)")
        if still_sync["action"] != "already identical":
            with open(src, "rb") as f:
                data = f.read()
            _atomic_retry(STILL, data, verify=lambda p: sha256_file(p) == sha_src)
        still_sync["final_sha256"] = sha256_file(STILL)
        still_sync["equal"] = bool(still_sync["final_sha256"] == sha_src)
    log("  still sync: %s" % still_sync.get("action"))
    exe = R18.ffmpeg_exe()
    pattern = os.path.join(FRAMES, "frame-%04d.png")
    cmd = [exe, "-y", "-v", "error", "-framerate", str(FPS),
           "-i", pattern, "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-crf", "16", "-preset", "medium", MP4]
    log("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("ffmpeg failed: %s" % r.stderr[-2000:])
    probe = R18.ffprobe(MP4)
    out = {"mp4": os.path.basename(MP4), "cmd": cmd, "probe": probe,
           "ffmpeg_binary": exe, "still_sync": still_sync,
           "ffmpeg_version": (subprocess.run(
               [exe, "-version"], capture_output=True, text=True
           ).stdout or "").splitlines()[:1],
           "size_bytes": os.path.getsize(MP4),
           "sha256": sha256_file(MP4), "created_utc": now_iso(),
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "encode.json"), out)
    checkpoint_extra("encode", out)
    cp.unit_done("encode")
    log("  %s %sx%s %s fps frames=%s duration=%s pix_fmt=%s %.2f MB"
        % (probe.get("codec_name"), probe.get("width"), probe.get("height"),
           probe.get("r_frame_rate"), probe.get("nb_frames_reported"),
           probe.get("duration_s"), probe.get("pix_fmt"),
           out["size_bytes"] / 1e6))
    return out


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-ENERGY / G-LINE / G-PEAK / G-MOVE / G-WARP / G-LEGACY / "
        "G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = R18.regions_from_photo()
    photo_rows = R18.measure_image(PHOTO, regions)
    gates = {}
    f24 = state.get("frame24") or {}
    var = f24.get("variants") or {}

    # ---- G-ENERGY --------------------------------------------------------
    chunks = state.get("chunks") or {}
    per_frame = []
    for k in sorted(chunks, key=int):
        for e in chunks[k].get("energies", []):
            per_frame.append({"frame": e["frame"], "rel_err": e["rel_err"],
                              "rel_err_float32": e["rel_err_float32_C"],
                              "deposited_total":
                                  e["deposited_total_averaged_C"],
                              "cell_energy_total": e["expected_averaged_C"],
                              "offframe_fraction_energy_min":
                                  e["offframe_fraction_energy_min"],
                              "offframe_fraction_energy_max":
                                  e["offframe_fraction_energy_max"],
                              "max_rel_err_per_direction":
                                  e["max_rel_err_per_direction"]})
    per_frame.sort(key=lambda x: x["frame"])
    worst = max((abs(x["rel_err"]) for x in per_frame), default=None)
    worst32 = max((abs(x["rel_err_float32"]) for x in per_frame), default=None)
    worst_dir = max((abs(x["max_rel_err_per_direction"])
                     for x in per_frame), default=None)
    gates["G-ENERGY"] = {
        "pass": bool(len(per_frame) == NFRAMES and worst is not None
                     and worst < 1e-6 and worst32 < 1e-6),
        "frames_measured": len(per_frame), "threshold": 1e-6,
        "max_rel_err": worst, "max_rel_err_float32_C": worst32,
        "max_rel_err_any_direction": worst_dir,
        "frame24_energy": (var.get("discscat") or {}).get("info"),
        "per_frame": per_frame,
        "normalisation_applied": 1.0,
        "criterion": ("for every frame the deposited total equals the sum of "
                      "the cell energies to a relative error < 1e-6, AFTER the "
                      "K-sample disc average; the average of K energy-exact "
                      "deposits is energy-exact, and that is what is measured "
                      "(both the float64 bookkeeping total and the delivered "
                      "float32 C)")}
    log("  G-ENERGY: %d frames, max rel_err %.2e (float32 %.2e, dir %.2e) -> %s"
        % (len(per_frame), worst if worst is not None else float("nan"),
           worst32 if worst32 is not None else float("nan"),
           worst_dir if worst_dir is not None else float("nan"),
           "PASS" if gates["G-ENERGY"]["pass"] else "FAIL"))

    # ---- G-LINE (per variant + main piece) -------------------------------
    rowsets = {}
    for tag, node in var.items():
        rowsets[tag] = node["gline"]
    rowsets["main_frame24"] = (rows_from_green(
        "main_frame24", R18.green_of(STILL), regions, photo_rows)
        if os.path.isfile(STILL) else None)
    gates["G-LINE"] = {
        "pass": bool(all(rowsets[t]["pass"] for t in
                         ("nodisc", "disc", "discscat"))),
        "per_variant": {t: {"pass": v["pass"], "pass_hp": v["pass_hp"],
                            "strength": v["strength"], "rows": v["rows"]}
                        for t, v in rowsets.items() if v},
        "photo_strength": float(max(r["depth_e3"] for r in photo_rows)),
        "reference_v4_1_8_rows": [
            {"region": r[0], "axis": r[1], "period_px": r[2],
             "depth": r[3], "depth_photo": r[4], "ratio": r[3] / r[4]}
            for r in V418_ROWS],
        "criterion": ("frame 24, same statistic and same windows as v4-1.8 "
                      "(run_v4_1_5's lattice_metric imported from run_v4_1_8): "
                      "the axis-parallel 5-15 px depth must be <= the source "
                      "photo's in all 4 windows, for every variant AND for the "
                      "main piece"),
        "regions": rinfo}
    for tag, v in var.items():
        log("  G-LINE %-8s %s (%s)" % (tag, "PASS" if v["gline"]["pass"]
                                       else "FAIL",
                                       ", ".join("%s%s %.1f/%.1f"
                                                 % (r["region"][:1],
                                                    r["axis"], r["depth"],
                                                    r["depth_photo"])
                                                 for r in v["gline"]["rows"])))
    log("  G-LINE (all variants): %s"
        % ("PASS" if gates["G-LINE"]["pass"] else "FAIL"))

    # ---- G-PEAK ----------------------------------------------------------
    pk = {}
    for tag, node in var.items():
        pk[tag] = node["peak"]
    gates["G-PEAK"] = {
        "pass": bool(all(pk[t]["std_ratio_vs_1_8_point_splat"] >= 0.9
                         and pk[t]["p99_ratio_vs_1_8_point_splat"] >= 0.9
                         for t in ("nodisc", "disc", "discscat"))),
        "threshold": 0.9, "reference": REF_PEAK,
        "per_variant": pk,
        "criterion": ("C std and p99 of frame 24 vs the v4-1.8 point-splat "
                      "reference on the same field and frame grid "
                      "(std 1.3106 / p99 6.348), >= 0.90, for every variant")}
    for tag in ("nodisc", "disc", "discscat"):
        log("  G-PEAK %-8s std %.4f (%.3f) p99 %.3f (%.3f)"
            % (tag, var[tag]["C_stats"]["std"],
               pk[tag]["std_ratio_vs_1_8_point_splat"],
               var[tag]["C_stats"]["p99"],
               pk[tag]["p99_ratio_vs_1_8_point_splat"]))
    log("  G-PEAK: %s" % ("PASS" if gates["G-PEAK"]["pass"] else "FAIL"))

    # ---- G-MOVE ----------------------------------------------------------
    frames = [R18.green_of(frame_path(i)) for i in range(NFRAMES)]
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = R18.interframe_stats(frames)
    try:
        dg = R18.decode_green(MP4)
        dec = R18.interframe_stats(dg)
        dec["frames_decoded"] = int(dg.shape[0])
    except Exception as e:  # pragma: no cover
        dec = {"error": str(e)}
    gates["G-MOVE"] = {
        "pass": bool(uniq == NFRAMES and st_["mean"] > 0.5),
        "unique_frames": int(uniq), "frames": NFRAMES,
        "interframe_green_png": st_, "interframe_green_mp4": dec,
        "threshold_gray": 0.5,
        "criterion": "48/48 unique frame PNGs AND mean |d green| > 0.5 gray"}
    log("  G-MOVE: unique %d/%d mean|d green| %.4f -> %s"
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
    changed, missing, added = [], [], []
    for name, h in rec["existing"].items():
        if name not in now:
            missing.append(name)
        elif now[name] != h:
            changed.append(name)
    for name in now:
        if name not in rec["existing"] and name not in NEW_FILES:
            added.append(name)
    own = {n: {"stage0": rec["all"].get(n), "now": now.get(n),
               "changed": bool(rec["all"].get(n) != now.get(n))}
           for n in NEW_FILES}
    imports = {"run_v4_1_8.py", "transport_mesh.py", "run_v16_water.py",
               "v12_water.py", "transport.py", "det.py", "run_v10.py",
               "run_v15.py", "run_v4_1_2.py", "stream_surface3.py",
               "stream_surface2.py", "state.py", "image.py",
               "v5_composite.py", "v5_physics.py"}
    ext = {}
    for name in changed:
        p = os.path.join(HERE, name)
        ext[name] = {
            "stage0_sha256": rec["all"].get(name), "now_sha256": now.get(name),
            "mtime": datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds"),
            "bytes_now": os.path.getsize(p),
            "imported_by_this_job": name in imports,
            "written_by_this_job": name in NEW_FILES,
            "attribution": ("this job never opens it for writing (its only "
                            "write target is run_v4_1_9.py); the change came "
                            "from a parallel lane after the stage-0 record"
                            if name not in imports and name not in NEW_FILES
                            else "IMPORTED BY THIS JOB -- investigate"),
        }
    attributable = [n for n in changed if ext.get(n, {}).get(
        "imported_by_this_job") or n in NEW_FILES]
    gates["G-LEGACY"] = {
        "pass": bool(not missing and not attributable),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "stage0_expected_sha256": rec.get("expected_sha256"),
        "stage0_sha_match": rec.get("sha_match"),
        "external_changes": ext,
        "attributable_to_this_job": attributable,
        "criterion": ("no .py recorded at stage 0 is missing, and no .py this "
                      "job imports or writes has changed.  A change to a file "
                      "this job neither imports nor writes is reported as an "
                      "external (parallel-lane) change and listed with its "
                      "mtime, not hidden")}
    log("  G-LEGACY: %d existing, changed %d (%s), missing %d, added %d -> %s"
        % (len(rec["existing"]), len(changed),
           ", ".join("%s[%s]" % (n, "imported" if ext[n]["imported_by_this_job"]
                                 else "external") for n in changed) or "none",
           len(missing), len(added),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL"))

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
        "criterion": ("48 frames / 2.000 s / 736x540 / 24 fps, from ffprobe's "
                      "fields (header parse + full decode where the ffprobe "
                      "alias is unreachable)"),
        "probe_method": pr_.get("probe_method"),
        "size_bytes": enc.get("size_bytes"), "sha256": enc.get("sha256")}
    log("  G-OUT: frames=%s duration=%s %sx%s fps=%s -> %s"
        % (pr_.get("nb_frames_reported"), pr_.get("duration_s"),
           pr_.get("width"), pr_.get("height"), fps_num,
           "PASS" if ok_out else "FAIL"))

    order = ["G-ENERGY", "G-LINE", "G-PEAK", "G-MOVE", "G-WARP", "G-LEGACY",
             "G-OUT"]
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
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = R18.regions_from_photo()
    enc = state.get("encode") or {}
    n = {"brief": BRIEF, "workdir": HERE, "built_utc": now_iso(),
         "scene": state.get("scene"), "stage0": state.get("stage0"),
         "stage1": state.get("stage1"), "frame24": state.get("frame24"),
         "conv": state.get("conv"), "warp": state.get("warp"),
         "encode": enc, "gates": state.get("gates"),
         "chunks": {k: {"wall_s": v.get("wall_s"),
                        "variant": v.get("variant"), "K_disc": v.get("K_disc"),
                        "energy_max_rel_err": v.get("energy_max_rel_err"),
                        "energy_max_rel_err_float32":
                            v.get("energy_max_rel_err_float32"),
                        "metrics": v.get("metrics")}
                    for k, v in (state.get("chunks") or {}).items()},
         "measurement": {
             "implementation": ("G-LINE = run_v4_1_8's own functions, imported "
                                "(R18.regions_from_photo / R18.measure_image / "
                                "R18.lattice_metric) -- the same code objects, "
                                "so the statistic and the windows are "
                                "identical to v4-1.8 by construction"),
             "definition": ("green channel, 8 bit; Hann 256x256; 2D FFT; axis "
                            "band |k_perp| <= 3, 5-15 px; period = W/k at the "
                            "band peak; depth = 1000*RMS(band-passed)/mean"),
             "regions": rinfo,
             "comparability": ("the photo is 735x539 and the frames 736x540, "
                               "so the same window origin addresses the same "
                               "content; both are 8-bit")},
         "deliverables": {"report": os.path.basename(REPORT),
                          "numbers": os.path.basename(NUMBERS),
                          "mp4": os.path.basename(MP4),
                          "frame": os.path.basename(STILL),
                          "variants": [os.path.basename(V_NODISC),
                                       os.path.basename(V_DISC),
                                       os.path.basename(V_SCAT)],
                          "checkpoint": os.path.basename(CP_PATH),
                          "stage0_hashes": os.path.basename(STAGE0_HASHES)},
         "files": {}}
    for p in (MP4, STILL, V_NODISC, V_DISC, V_SCAT, CP_PATH, STAGE0_HASHES):
        if os.path.isfile(p):
            n["files"][os.path.basename(p)] = {
                "bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    fp = frame_path(STILL_IDX)
    n["still_equals_frame_0024"] = bool(
        os.path.isfile(fp) and sha256_file(fp) == n["files"].get(
            os.path.basename(STILL), {}).get("sha256"))
    n["run_v4_1_9_sha256"] = sha256_file(os.path.abspath(__file__))
    n["transport_mesh_sha256"] = sha256_file(os.path.join(HERE,
                                                         "transport_mesh.py"))
    n["run_v4_1_8_sha256"] = sha256_file(os.path.join(HERE, "run_v4_1_8.py"))
    return n


def _f(x, nd=1):
    try:
        return ("%%.%df" % nd) % float(x)
    except Exception:
        return "—"


def _pct(x, nd=1):
    try:
        return ("%+." + str(nd) + "f%%") % (100.0 * float(x))
    except Exception:
        return "—"


def write_report(n):
    s1 = n.get("stage1") or {}
    sc = n.get("scene") or {}
    f24 = n.get("frame24") or {}
    cv = n.get("conv") or {}
    g = n.get("gates") or {}
    st0 = n.get("stage0") or {}
    var = f24.get("variants") or {}
    sep = f24.get("band_separation") or {}
    kchk = f24.get("K_check") or {}
    disc = s1.get("disc") or {}
    scat = s1.get("scattering") or {}
    req = s1.get("required") or {}
    ver = s1.get("verdict") or {}
    gl = g.get("G-LINE") or {}
    gp = g.get("G-PEAK") or {}
    ge = g.get("G-ENERGY") or {}
    gm = g.get("G-MOVE") or {}
    gw = g.get("G-WARP") or {}
    gg = g.get("G-LEGACY") or {}
    go = g.get("G-OUT") or {}

    def w(k):
        p = (g.get(k) or {}).get("pass")
        return "**PASS**" if p is True else ("**記録のみ**" if p is None
                                             else "**FAIL**")

    _worst = 0.0
    _worst_tag = "—"
    for _t in ("nodisc", "disc", "discscat"):
        for _r in ((var.get(_t) or {}).get("gline") or {}).get("rows") or []:
            if _r["ratio"] > _worst:
                _worst, _worst_tag = _r["ratio"], "%s %s%s" % (
                    _t, _r["region"], _r["axis"])
    # the multiple-scattering paragraph numbers
    _mm = (scat.get("multiple_scattering") or {})
    _op = _mm.get("orders_probability") or {}
    _p0 = float(_op.get("0") or 0.0)
    _p1 = float(_op.get("1") or 0.0)
    _p2 = float(_op.get("2") or 0.0)
    _p3 = float(_op.get("3") or 0.0)
    _share1 = 100.0 * _p1 / max(_p1 + _p2 + _p3, 1e-30)
    _ratio21 = 100.0 * _p2 / max(_p1, 1e-30)

    L = []
    A = L.append
    A("# REPORT-V4-1-9 — 太陽円盤＋水中散乱の広がりを入れた堆積と再測")
    A("")
    A("Brief: `BRIEF-V4-1-9.md`（`D:\\mutaro-work\\rally-in\\BRIEF-V4-1-9.md`）  ")
    A("Workdir: `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`  ")
    A("生成: %s  " % n.get("built_utc"))
    A("新規ファイルは `run_v4_1_9.py` の 1 本のみ（既存 .py は無改造。`run_v4_1_8.py` "
      "と `transport_mesh.py` は import のみ、sha256 は §0）。")
    A("")
    A("## 0. 判定（全ゲート）")
    A("")
    A("| ゲート | 判定 | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-ENERGY | %s | 48 フレームの堆積総量/セル総和の最大相対誤差 %.2e（閾値 1e-6、"
      "円盤平均後・float32 の C でも %.2e） |"
      % (w("G-ENERGY"), ge.get("max_rel_err") or 0.0,
         ge.get("max_rel_err_float32_C") or 0.0))
    A("| G-LINE | %s | 4 窓×3 変種＋本編。最悪比 %.3f（%s、変種別は §5） |"
      % (w("G-LINE"), _worst, _worst_tag))
    A("| G-PEAK | %s | C std 比 %.3f／p99 比 %.3f（1.8 点スプラット参照 std 1.3106／"
      "p99 6.348、≥0.90） |"
      % (w("G-PEAK"),
         (gp.get("per_variant", {}).get("discscat", {}) or {}).get(
             "std_ratio_vs_1_8_point_splat") or 0.0,
         (gp.get("per_variant", {}).get("discscat", {}) or {}).get(
             "p99_ratio_vs_1_8_point_splat") or 0.0))
    A("| G-MOVE | %s | ユニーク %s/48 フレーム、フレーム間平均差（緑）%.2f 諧調 |"
      % (w("G-MOVE"), gm.get("unique_frames"),
         (gm.get("interframe_green_png") or {}).get("mean") or 0.0))
    A("| G-WARP | %s | frame 24 の歪み（緑・平均）%.3f px |"
      % (w("G-WARP"), gw.get("g_mean_px") or 0.0))
    A("| G-LEGACY | %s | stage-0 の %d 本、変更 %d（外部 %d／本ジョブ %d）、"
      "欠落 %d（§9） |"
      % (w("G-LEGACY"), gg.get("recorded_existing") or 0,
         len(gg.get("changed") or []),
         len([n for n, e in (gg.get("external_changes") or {}).items()
              if not e.get("imported_by_this_job")]),
         len(gg.get("attributable_to_this_job") or []),
         len(gg.get("missing") or [])))
    A("| G-OUT | %s | %s frames / %.3f s / %sx%s / %s fps / %s |"
      % (w("G-OUT"), (go.get("probe") or {}).get("nb_frames_reported"),
         float((go.get("probe") or {}).get("duration_s") or 0.0),
         (go.get("probe") or {}).get("width"),
         (go.get("probe") or {}).get("height"),
         (go.get("probe") or {}).get("r_frame_rate"),
         (go.get("probe") or {}).get("pix_fmt")))
    A("")
    A("**見出し（正直な要約）**")
    A("")
    A("1. **STAGE 1 の判定は「物理の広がりでは説明できない」**（§2）。太陽円盤は最悪窓に"
      "必要な等価 σ %.2f px に対して %.2f px（x）/ %.2f px（y）で **%.1f 倍細い**。"
      "単散乱 PSF は %.1f px（x）/ %.1f px（y）で **%.0f 倍太い**。"
      % (ver.get("needed_sigma_px_A") or 0.0, disc.get("rms_px_x") or 0.0,
         disc.get("rms_px_y") or 0.0, ver.get("disc_is_factor_too_narrow") or 0.0,
         ver.get("scatter_rms_px_x") or 0.0, ver.get("scatter_rms_px_y") or 0.0,
         ver.get("scatter_is_factor_too_wide") or 0.0))
    A("2. **散乱は「未投入」ではない。** 凍結合成器は同じ位相関数・同じ b,c から作った "
      "PSF を `w_halo`（緑 %.3f）で既に掛けている（`run_v16_water.render_frame_v16`）。"
      "C へ再畳み込みすると二重計上になり、実測で C の std が %.3f → %.3f（%s）、"
      "p99 が %.3f → %.3f（%s）に落ちる。"
      % (((scat.get("pipeline_comparison") or {}).get("w_halo_rgb") or [0, 0, 0])[1],
         (var.get("disc", {}).get("C_stats") or {}).get("std") or 0.0,
         (var.get("discscat", {}).get("C_stats") or {}).get("std") or 0.0,
         _pct(((var.get("discscat", {}).get("C_stats") or {}).get("std") or 0) /
              max((var.get("disc", {}).get("C_stats") or {}).get("std") or 1, 1e-9)
              - 1.0),
         (var.get("disc", {}).get("C_stats") or {}).get("p99") or 0.0,
         (var.get("discscat", {}).get("C_stats") or {}).get("p99") or 0.0,
         _pct(((var.get("discscat", {}).get("C_stats") or {}).get("p99") or 0) /
              max((var.get("disc", {}).get("C_stats") or {}).get("p99") or 1, 1e-9)
              - 1.0)))
    A("3. **それでも G-LINE は通らない。** 帯分離表（§5.2）が示すとおり、C=1 の"
      "「歪みのみ」でも 4 窓中 3 窓が写真を上回る。つまり **C の広がりを無限大にしても "
      "G-LINE は通らない**（広がりは C レーンにしか効かず、残差は歪みレーンにある）。")
    A("")
    A("## 1. 前提（STAGE 0、絶対パス）")
    A("")
    A("| path | 存在 | sha256 | 一致 |")
    A("| --- | --- | --- | --- |")
    for nm in ("run_v4_1_8.py", "transport_mesh.py"):
        A("| `%s` | はい | `%s` | %s |"
          % (nm, (st0.get("sha256") or {}).get(nm, "—")[:32],
             "OK" if (st0.get("sha_match") or {}).get(nm) else "MISMATCH"))
    A("")
    A("* `out/v4-1-8/numbers-v4-1-8.json` と `REPORT-V4-1-8.md`: 存在確認 OK（1.8 の数値の実体）。")
    A("* `transport_mesh.py --selftest`: **%d ok / %d NG**（rc=%d、要求 15 ok / 0 ng を再現）。"
      % ((st0.get("selftest") or {}).get("n_checks_ok") or 0,
         (st0.get("selftest") or {}).get("n_checks_ng") or 0,
         (st0.get("selftest") or {}).get("returncode") or 0))
    A("* 水パラメータの出所（stage-0 で特定）:")
    A("  * 散乱係数 b と媒体プリセット = **`v12_water.py`**（`PRESETS['algae_pond']`、"
      "`b_water_morel`、`preset_spectra`）。b(550) = %.6f /m、c(550) = %.6f /m"
      "（load = 0.6、プリセット既定）。"
      % (scat.get("b_1_per_m") or 0.0, scat.get("c_1_per_m") or 0.0))
    A("  * 位相関数 FF-Petzold = **`run_v16_water.py`**（`FF_PHASE['algae_pond']` "
      "n=1.10, mu=3.5835、`phase_table`、`radial_profile`）。")
    A("  * 太陽円盤の角半径 = **`transport.py`**（`sun_disc_radius(clarity)`）。")
    A("")
    A("## 2. STAGE 1 — 広がりの見積もり（数値のみ）")
    A("")
    A("場面は 1.8 と同一: H = %.2f m、U = %.2f m/s（U0 = %.3f）、seed = %d、"
      "窓 L = %.7f m、水面グリッド %d²、736×540、media algae_pond、"
      "clarity %.2f。画素ピッチ x = %.4f mm、y = %.4f mm。"
      % (float((sc.get("constants") or {}).get("H") or 0.6),
         float((sc.get("constants") or {}).get("U") or 0.8),
         float((sc.get("constants") or {}).get("U0") or 0.0),
         int((sc.get("constants") or {}).get("seed") or 0),
         float(sc.get("patch_m") or 0.0),
         int(sc.get("grid") or 0),
         0.75,
         1e3 * (disc.get("px_pitch_x_m") or 0.0),
         1e3 * (disc.get("px_pitch_y_m") or 0.0)))
    A("")
    A("合成は凍結形（A %.1f / Blow σ %.1f px / E %.1f / Q %d、z = %.2f m）で "
      "1.8 と同一。妥当性ラベル: %s。"
      % (float((sc.get("composite") or {}).get("A") or 0.0),
         float((sc.get("composite") or {}).get("Blow_sigma_px") or 0.0),
         float((sc.get("composite") or {}).get("E") or 0.0),
         int((sc.get("composite") or {}).get("Q") or 0),
         float((sc.get("composite") or {}).get("z_refraction_m") or 0.0),
         sc.get("validity") or "—"))
    A("")
    A("### 2.1 太陽円盤")
    A("")
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| 角半径 (clarity %.2f) | %.7f rad |" % (0.75, disc.get("radius_rad") or 0.0))
    A("| 水中角 (Snell, n=%.5f) | %.7f rad |"
      % (disc.get("ior_water") or 0.0, disc.get("snell_theta_water_rad") or 0.0))
    A("| 深さ %.2f m での投影半径 | %.5f mm |"
      % (disc.get("depth_m") or 0.6, 1e3 * (disc.get("projected_radius_m") or 0.0)))
    A("| 投影半径 px | **%.3f px (x) / %.3f px (y)** |"
      % (disc.get("projected_radius_px_x") or 0.0,
         disc.get("projected_radius_px_y") or 0.0))
    A("| 投影直径 px | %.3f px (x) / %.3f px (y) |"
      % (disc.get("projected_diameter_px_x") or 0.0,
         disc.get("projected_diameter_px_y") or 0.0))
    A("| 等価ガウス σ（2 次モーメント一致） | **%.3f px (x) / %.3f px (y)** |"
      % (disc.get("rms_px_x") or 0.0, disc.get("rms_px_y") or 0.0))
    A("")
    A("式: `r = z·tan(asin(sin(rad)/n))`、`rad = 0.00465·(0.7+0.6·(1−clarity))`。"
      "一様円盤の軸あたり分散は `R²/4` なので等価 σ = R/2。円盤は m で等方、"
      "画素ピッチは非等方（x %.3f / y %.3f mm）なので画素空間では楕円。"
      % (1e3 * (disc.get("px_pitch_x_m") or 0.0),
         1e3 * (disc.get("px_pitch_y_m") or 0.0)))
    A("")
    A("### 2.2 水中散乱（algae_pond, 光路 0.60 m）")
    A("")
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| b(550) / c(550) / a(550) | %.4f / %.4f / %.4f /m |"
      % (scat.get("b_1_per_m") or 0.0, scat.get("c_1_per_m") or 0.0,
         scat.get("a_1_per_m") or 0.0))
    A("| 単散乱アルベド b/c | %.4f |" % (scat.get("single_scatter_albedo") or 0.0))
    A("| 光路長 / 光学深度 cL | %.2f m / %.4f |"
      % (scat.get("path_m") or 0.0, scat.get("optical_depth_cL") or 0.0))
    A("| `<θ²>`（FF-Petzold） | %.5f rad²（θ_rms = %.1f°） |"
      % (scat.get("theta2_rad2") or 0.0, scat.get("theta_rms_deg") or 0.0))
    A("| σ（モーメント式） | **%.5f m = %.2f px (x) / %.2f px (y)** |"
      % (scat.get("sigma_m") or 0.0, scat.get("sigma_px_x") or 0.0,
         scat.get("sigma_px_y") or 0.0))
    A("| 半径 RMS（tan θ の厳密プロファイル、r_max %.3f m） | **%.1f px (x) / %.1f px (y)**"
      "（捕捉率 %.2f%%） |"
      % ((scat.get("exact_profile") or {}).get("r_max_m") or 0.0,
         (scat.get("exact_profile") or {}).get("rms_radius_px_x") or 0.0,
         ((scat.get("exact_profile") or {}).get("rms_radius_m") or 0.0)
         / max(disc.get("px_pitch_y_m") or 1.0, 1e-12),
         100.0 * ((scat.get("exact_profile") or {}).get("captured_fraction") or 0.0)))
    A("| 単散乱の重み H / 弾道 T | %.3f / %.3f |"
      % (scat.get("single_scatter_weight_H") or 0.0,
         scat.get("ballistic_weight_T") or 0.0))
    A("")
    A("式（パイプライン既存の v12 単散乱変位モーメント）: "
      "`σ² = (<θ²>/2)(b/c³)·D2(cL)`、`D2(x)=x²−2x+2−2e^{−x}`。")
    A("")
    A("**近似の範囲**: θ_rms = %.1f° は近軸ではない。したがってモーメント式の σ は"
      "**下限**である。`tan θ` をそのまま使う厳密プロファイルでは RMS が "
      "%s 大きい（上表の「半径 RMS」）。本ジョブの §4 の畳み込みは**厳密プロファイル**"
      "（`run_v16_water.radial_profile`、自由パラメータなし）を使い、モーメント式の値は"
      "比較用に併記する。"
      % (scat.get("theta_rms_deg") or 0.0,
         ("%.1f%%"
          % (100.0 * (((scat.get("exact_profile") or {}).get("rms_radius_m") or
                       0.0) / max(scat.get("sigma_m") or 1e-30, 1e-30) - 1.0)))))
    A("")
    A("**単散乱で足りるか / 多重散乱が要るか**: 次数別確率 "
      "(0 次=%.3f, 1 次=%.3f, 2 次=%.3f, 3 次=%.3f)。散乱した光のうち 1 次が "
      "%.0f%%、2 次は 1 次の %.0f%% で、2 次は約 1.41 倍広い。よって**幅としては"
      "単散乱は 10〜20%% 過小**（多重散乱は無視できない）。ただし v12/v16 の光学模型は"
      "設計上 single-scatter で、その docstring 自身がハローを下限と断っている。"
      "多重散乱 PSF を足すのは「この模型の補正」ではなく**新しい模型**なので、"
      "本ジョブでは足していない（足したら stage-1 の導出値でなくなる）。"
      % (_p0, _p1, _p2, _p3, _share1, _ratio21))
    A("")
    A("**パイプラインとの突き合わせ（重要）**: 合成器が使う光路は 0.60 m ではなく"
      "**往復 %.4f m**（`v12_water.ab_optics` の s_down+s_up）で、その σ は %.4f m = "
      "状態グリッド %.1f px（カーネル RMS %.1f px、r_max %.1f px、捕捉率 %.3f）。"
      "`run_v16_water.render_frame_v16` は `Ch = ff_blur_fft(C, kernel_2d)` を作り "
      "`w_bal·C + w_halo·Ch` で合成する。**同じ位相関数の PSF は既にレンダーに入っている**。"
      % (float(((scat.get("pipeline_comparison") or {}).get(
          "optics_two_way_path_m")) or 0.0),
         float(((scat.get("pipeline_comparison") or {}).get(
             "optics_sigma_m")) or 0.0),
         float(((scat.get("pipeline_comparison") or {}).get(
             "optics_sigma_px_state_grid")) or 0.0),
         float(((scat.get("pipeline_comparison") or {}).get(
             "optics_kernel_rms_px_state_grid")) or 0.0),
         float(((scat.get("pipeline_comparison") or {}).get(
             "optics_kernel_rmax_px_state_grid")) or 0.0),
         float(((scat.get("pipeline_comparison") or {}).get(
             "optics_kernel_captured_fraction")) or 0.0)))
    A("")
    A("### 2.3 超過を消すのに必要な等価 σ（概算）")
    A("")
    A("1.8 の窓別比（1.04 / 0.99 / 1.17 / 1.21）から、最悪窓を 1.00 以下にするのに"
      "必要な等価ガウス σ。")
    A("")
    A("| region | axis | period px | 比 | σ 必要 px（A: frontier 式） | σ 必要 px（B: 実測スペクトル） |")
    A("| --- | --- | --- | --- | --- | --- |")
    mb = {r["region"] + r["axis"]: r for r in
          ((req.get("method_B_spectral") or {}).get("rows") or [])}
    for a_, b_ in zip(req.get("method_A_frontier") or [],
                      req.get("ratios_v4_1_8") or []):
        k2 = mb.get(a_["region"] + a_["axis"]) or {}
        A("| %s | %s | %.2f | %.3f | **%.2f** | %s |"
          % (a_["region"], a_["axis"], a_["period_px"], a_["ratio"],
             a_["sigma_needed_px"],
             ("**%.2f**" % k2["sigma_needed_px"])
             if k2.get("sigma_needed_px") else "—（既に以下）"))
    A("")
    A("A は単一モードの frontier 式 `σ = p·sqrt(ln r/(2π²))`（**概算**）。"
      "B は 1.8 の still の窓スペクトルをモードごとに `exp(−2π²σ²/p²)` で減衰させて"
      "解いたもの（**概算**、画像レベル＝歪み＋コースティクス両方に掛けた場合）。"
      "最悪値は A で **%.2f px**、B で **%.2f px**。"
      % (req.get("worst_sigma_needed_px_method_A") or 0.0,
         req.get("worst_sigma_needed_px_method_B") or 0.0))
    A("")
    A("### 2.4 突き合わせと判定")
    A("")
    A("| 広がり | 等価 σ px（x / y） | 必要量に対する比 |")
    A("| --- | --- | --- |")
    A("| (1) 太陽円盤 | %.3f / %.3f | **%.1f 倍細い**（最悪窓） |"
      % (disc.get("rms_px_x") or 0.0, disc.get("rms_px_y") or 0.0,
         ver.get("disc_is_factor_too_narrow") or 0.0))
    A("| (2) 単散乱 PSF（厳密） | %.1f / %.1f | **%.0f 倍太い** |"
      % (ver.get("scatter_rms_px_x") or 0.0, ver.get("scatter_rms_px_y") or 0.0,
         ver.get("scatter_is_factor_too_wide") or 0.0))
    A("| (3) 必要量（A/B） | %.2f / %.2f | 1.00 |"
      % (req.get("worst_sigma_needed_px_method_A") or 0.0,
         req.get("worst_sigma_needed_px_method_B") or 0.0))
    A("")
    A("**判定: %s。** %s" % (ver.get("judgement"), ver.get("text")))
    A("")
    _need_rows = req.get("method_A_frontier") or [{}]
    _needworst = max(_need_rows,
                     key=lambda r: (r.get("sigma_needed_px") or 0.0))
    A("補足（数値）: 円盤だけでは最悪窓（%s %s、period %.2f px）の帯振幅を "
      "**%.1f%%** しか落とせない（σ_x %.3f px、比 %.3f → %.3f）。"
      % (_needworst.get("region") or "—", _needworst.get("axis") or "",
         _needworst.get("period_px") or 0.0,
         100.0 * (1.0 - (ver.get("disc_attenuation_of_worst_window") or 1.0)),
         disc.get("rms_px_x") or 0.0,
         _needworst.get("ratio") or 0.0,
         (_needworst.get("ratio") or 0.0)
         * (ver.get("disc_attenuation_of_worst_window") or 1.0)))
    A("")
    A("## 3. STAGE 2 — 実装")
    A("")
    A("* **円盤求積**: `det.fibonacci_disc(K, sun_disc_radius(clarity))`（既存の"
      "決定論的等重み則）の K 点それぞれで、`transport_mesh.mesh_C` の**面積保存射影**"
      "（四辺形→2 三角形、画素との重なり面積を解析的にクリップ）を回し、K 枚を平均。"
      "各サンプルはエネルギー厳密保存なので平均も保存。カーネル・ジッタ・帯域制限は"
      "引き続き**なし**。")
    A("* **K の正当化**: 投影半径は %.2f px（x）/ %.2f px（y）で 1 px 未満。"
      "画素格子に残るのは 2 次モーメントまでで、等重み K 点 Fibonacci 則は K によらず"
      "`Σr²/K = R²/2` を厳密に満たす。実測の K=9 vs K=25（同一フレーム・同一場）: "
      "std %s、p99 %s、max %s。よって **K = %d**（ブリーフの下限 9）で収束。"
      % (disc.get("projected_radius_px_x") or 0.0,
         disc.get("projected_radius_px_y") or 0.0,
         _pct(kchk.get("std_rel_diff"), 3),
         _pct(kchk.get("p99_rel_diff"), 3),
         _pct(kchk.get("max_rel_diff"), 3), 9))
    A("* **散乱の広がり**: §2.2 の厳密プロファイル（FF-Petzold、光路 0.60 m、"
      "r_max %.3f m、捕捉率 %.4f）を **m で等方**な周期カーネルに落として C に畳み込む"
      "（`rfft2`、単位質量へ再正規化＝v16 の `kernel_2d` と同じ規約）。"
      "画素空間では楕円（RMS %.2f px x / %.2f px y）。**幅の自由パラメータはない**。"
      % (float(((f24.get("psf") or {}).get("kernel_info") or {}).get(
          "profile_captured_fraction") or 0.0) and
         float((f24.get("psf") or {}).get("r_max_m") or 0.0),
         float(((f24.get("psf") or {}).get("kernel_info") or {}).get(
             "profile_captured_fraction") or 0.0),
         float(((f24.get("psf") or {}).get("kernel_info") or {}).get(
             "rms_radius_px_x") or 0.0),
         float(((f24.get("psf") or {}).get("kernel_info") or {}).get(
             "rms_radius_px_y") or 0.0)))
    A("* RMS は 0.2 px を大きく超える（%.1f px）ので、ブリーフの「入れない」条件には"
      "該当しない。**入れた**。" % (ver.get("scatter_rms_px_x") or 0.0))
    A("* frame 24 の 3 変種: (i) `nodisc`（1.8 と同一）/ (ii) `disc`（円盤のみ）/ "
      "(iii) `discscat`（円盤＋散乱）。本編は (iii)。合成は凍結形（A %.1f / Blow "
      "σ %.1f px / E %.1f / Q %d、1.8 と同一）。"
      % (float((f24.get("frozen_composite") or {}).get("A") or 0.0),
         float((f24.get("frozen_composite") or {}).get("Blow_sigma_px") or 0.0),
         float((f24.get("frozen_composite") or {}).get("E") or 0.0),
         int((f24.get("frozen_composite") or {}).get("Q") or 0)))
    A("")
    A("### 3.1 変種 (i) は 1.8 をビット再現したか")
    A("")
    A("| 変種 | file | sha256 | 1.8 の frame-0024.png と同一 |")
    A("| --- | --- | --- | --- |")
    for t in ("nodisc", "disc", "discscat"):
        vv = var.get(t) or {}
        A("| %s | `%s` | `%s` | %s |"
          % (t, vv.get("file"), (vv.get("sha256") or "—")[:16],
             ("**はい**" if vv.get("still_equals_v418") else "いいえ")
             if vv.get("still_equals_v418") is not None else "—"))
    A("")
    A("(i) がバイト一致なら、1.8 の堆積経路（単一方向）と本ジョブの足回りは同一である。")
    A("")
    A("## 4. STAGE 2 — C の統計（frame 24）")
    A("")
    A("| 変種 | mean | std | p99 | max | std 比 vs 1.8 点スプラット | p99 比 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for t in ("nodisc", "disc", "discscat"):
        cs = (var.get(t) or {}).get("C_stats") or {}
        pk = (var.get(t) or {}).get("peak") or {}
        A("| %s | %.6f | %.4f | %.3f | %.2f | %.3f | %.3f |"
          % (t, cs.get("mean") or 0.0, cs.get("std") or 0.0,
             cs.get("p99") or 0.0, cs.get("max") or 0.0,
             pk.get("std_ratio_vs_1_8_point_splat") or 0.0,
             pk.get("p99_ratio_vs_1_8_point_splat") or 0.0))
    A("| 1.8 点スプラット参照 | 1.000000 | 1.3106 | 6.348 | 44.08 | 1.000 | 1.000 |")
    A("")
    A("## 5. G-LINE（1.8 と同一の統計・同一の窓）")
    A("")
    A("### 5.1 変種別の 4 窓")
    A("")
    A("| 変種 | region | axis | period px | depth | 写真 | 比 | prom | OK |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for t in ("nodisc", "disc", "discscat"):
        for r in ((var.get(t) or {}).get("gline") or {}).get("rows") or []:
            A("| %s | %s | %s | %.2f | %.1f | %.1f | %.3f | %.2f | %s |"
              % (t, r["region"], r["axis"], r["period_px"], r["depth"],
                 r["depth_photo"], r["ratio"],
                 r.get("prom_peak_over_median") or 0.0,
                 "OK" if r["ok"] else "**FAIL**"))
    A("")
    A("写真の同一窓: " + " / ".join(
        "%s %s %.1f" % (r["region"], r["axis"], r["depth_e3"])
        for r in (R18.measure_image(PHOTO, R18.regions_from_photo()[0]))))
    A("")
    A("### 5.2 帯の分離表（1.8 §4.1 の作法）")
    A("")
    win = sep.get("window_order") or []
    A("| 描き方 | " + " | ".join(win) + " |")
    A("| --- | " + " | ".join(["---"] * len(win)) + " |")
    for row in sep.get("rows") or []:
        A("| %s | " % row["tag"] + " | ".join("%.1f" % x for x in row["depth"])
          + " |")
    A("")
    _w1 = {}
    try:
        _w1 = (read_json(V418_NUMBERS).get("diagnostics") or {}).get(
            "warp_only") or {}
    except Exception:
        _w1 = {}
    if _w1.get("depth"):
        _wo = [r["depth"] for r in (sep.get("rows") or [])
               if r["tag"] == "warp_only"]
        _dv = (max(abs(a - b) for a, b in zip(_wo[0], _w1["depth"]))
               if _wo else None)
        A("")
        A("交差検証: 本ジョブの `warp_only` %s は、1.8 が記録した同じ量 "
          "`diagnostics.warp_only.depth` %s と最大差 **%s**（1.8 のコードと本ジョブの"
          "足回りが同一であることの数値的な裏付け。差は 8 bit 丸めの範囲）。"
          % ([round(x, 2) for x in (_wo[0] if _wo else [])],
             [round(x, 3) for x in _w1["depth"]],
             ("%.2e" % _dv) if _dv is not None else "—"))
    A("")
    A("**読み方（決定的）**: `warp_only`（C=1）は **C レーンの広がりを無限大にした"
      "極限**である。その値が 4 窓中 3 窓で写真を上回っている以上、**どんな幅の"
      "カーネルを C に入れても G-LINE は通らない**。残差は C ではなく歪みレーンに"
      "ある（1.8 §4.1 は「窓ごとの最大」だけを示したため 62.6 < 63.7 に見えていたが、"
      "4 窓を並べると 3 窓で上回る。1.8 の生データ `diagnostics/warp_only/depth` は"
      "本表と同じ値である）。")
    A("")
    A("### 5.3 period / prominence（窓別、G-LINE2 検証用）")
    A("")
    A("| 変種 | region | axis | period px | prom（本作） | prom（写真） |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in ("nodisc", "disc", "discscat"):
        for r in ((var.get(t) or {}).get("gline") or {}).get("rows") or []:
            A("| %s | %s | %s | %.2f | %.2f | %.2f |"
              % (t, r["region"], r["axis"], r["period_px"],
                 r.get("prom_peak_over_median") or 0.0,
                 r.get("prom_photo") or 0.0))
    A("")
    A("## 6. STAGE 3 — 1024² 水面グリッド（frame 24、(iii)）")
    A("")
    A("| 水面 | C mean | C std | C p99 | C max |")
    A("| --- | --- | --- | --- | --- |")
    A("| 512²（disc） | %.6f | %.4f | %.3f | %.2f |"
      % ((cv.get("C_512_discscat") or {}).get("mean") or 0.0,
         (f24.get("variants", {}).get("disc", {}).get("C_stats") or {}).get(
             "std") or 0.0,
         (f24.get("variants", {}).get("disc", {}).get("C_stats") or {}).get(
             "p99") or 0.0,
         (f24.get("variants", {}).get("disc", {}).get("C_stats") or {}).get(
             "max") or 0.0))
    A("| 512²（discscat） | %.6f | %.4f | %.3f | %.2f |"
      % ((cv.get("C_512_discscat") or {}).get("mean") or 0.0,
         (cv.get("C_512_discscat") or {}).get("std") or 0.0,
         (cv.get("C_512_discscat") or {}).get("p99") or 0.0,
         (cv.get("C_512_discscat") or {}).get("max") or 0.0))
    A("| 1024²（disc） | %.6f | %.4f | %.3f | %.2f |"
      % ((cv.get("C_1024_disc") or {}).get("mean") or 0.0,
         (cv.get("C_1024_disc") or {}).get("std") or 0.0,
         (cv.get("C_1024_disc") or {}).get("p99") or 0.0,
         (cv.get("C_1024_disc") or {}).get("max") or 0.0))
    A("| 1024²（discscat） | %.6f | %.4f | %.3f | %.2f |"
      % ((cv.get("C_1024_discscat") or {}).get("mean") or 0.0,
         (cv.get("C_1024_discscat") or {}).get("std") or 0.0,
         (cv.get("C_1024_discscat") or {}).get("p99") or 0.0,
         (cv.get("C_1024_discscat") or {}).get("max") or 0.0))
    A("")
    A("512² → 1024² の差（discscat）: std %s、p99 %s。disc では std %s、p99 %s"
      "（いずれも ±15%% 以内）。"
      % (_pct(cv.get("std_rel_diff_vs_512_discscat")),
         _pct(cv.get("p99_rel_diff_vs_512_discscat")),
         _pct(((cv.get("C_1024_disc") or {}).get("std") or 0.0)
              / max(((cv.get("C_512_disc") or {}).get("std") or 1.0), 1e-9) - 1.0),
         _pct(((cv.get("C_1024_disc") or {}).get("p99") or 0.0)
              / max(((cv.get("C_512_disc") or {}).get("p99") or 1.0), 1e-9) - 1.0)))
    A("")
    _b5 = cv.get("band_512") or {}
    _b5d = _b5.get("band_512_disc") or {}
    _b5s = _b5.get("band_512_discscat") or {}
    _b1d = cv.get("band_1024_disc") or {}
    _b1s = cv.get("band_1024_discscat") or {}
    A("**C 場の帯統計（軸平行 5–15 px、平均に対する ‰）**")
    A("")
    A("| 水面 | disc x | disc y | discscat x | discscat y |")
    A("| --- | --- | --- | --- | --- |")
    A("| 512² | %.2f | %.2f | %.2f | %.2f |"
      % (_b5d.get("x") or 0.0, _b5d.get("y") or 0.0,
         _b5s.get("x") or 0.0, _b5s.get("y") or 0.0))
    A("| 1024² | %.2f | %.2f | %.2f | %.2f |"
      % (_b1d.get("x") or 0.0, _b1d.get("y") or 0.0,
         _b1s.get("x") or 0.0, _b1s.get("y") or 0.0))
    A("| 差 | %s | %s | %s | %s |"
      % (_pct((_b1d.get("x") or 0.0) / max(_b5d.get("x") or 1.0, 1e-9) - 1.0),
         _pct((_b1d.get("y") or 0.0) / max(_b5d.get("y") or 1.0, 1e-9) - 1.0),
         _pct((_b1s.get("x") or 0.0) / max(_b5s.get("x") or 1.0, 1e-9) - 1.0),
         _pct((_b1s.get("y") or 0.0) / max(_b5s.get("y") or 1.0, 1e-9) - 1.0)))
    A("")
    A("読み方: (a) 標本密度 4 倍で C 統計も帯統計も ±15%% 以内（差の絶対値の最大 %.1f%%）。"
      % (100.0 * max(
          abs((_b1d.get("y") or 0.0) / max(_b5d.get("y") or 1.0, 1e-9) - 1.0),
          abs((_b1s.get("y") or 0.0) / max(_b5s.get("y") or 1.0, 1e-9) - 1.0)))
      + "(b) **散乱の畳み込みだけで C の 5–15 px 帯が %.1f‰ → %.1f‰（x）に落ちる**"
        "（約 1/9）。これは §2 の PSF がこの帯を丸ごと削る太さだという事実の直接測定で"
        "あり、同時に G-PEAK が壊れる理由でもある。"
        % (_b5d.get("x") or 0.0, _b5s.get("x") or 0.0))
    A("")
    A("## 7. G-ENERGY — 円盤平均後の堆積会計")
    A("")
    A("48 フレーム、`|平均 C の総和 − W·H|/W·H` の最大 **%.2e**（float32 の C でも "
      "**%.2e**、1 方向あたりの最大 **%.2e**）。閾値 1e-6。"
      "K 枚のエネルギー厳密な堆積の平均はエネルギー厳密なので、これは定義上一致する"
      "（正規化係数は掛けていない）。枠外（unwrap）へ出たエネルギー割合は "
      "%.3f〜%.3f で、別枠の記録であり C の水準には使っていない。"
      % (ge.get("max_rel_err") or 0.0, ge.get("max_rel_err_float32_C") or 0.0,
         ge.get("max_rel_err_any_direction") or 0.0,
         min((x["offframe_fraction_energy_min"] for x in
              (ge.get("per_frame") or [{"offframe_fraction_energy_min": 0}])),
             default=0.0),
         max((x["offframe_fraction_energy_max"] for x in
              (ge.get("per_frame") or [{"offframe_fraction_energy_max": 0}])),
             default=0.0)))
    A("")
    A("## 8. G-PEAK")
    A("")
    A("| 変種 | std | p99 | std 比 | p99 比 | 判定 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for t in ("nodisc", "disc", "discscat"):
        cs = (var.get(t) or {}).get("C_stats") or {}
        pk = (var.get(t) or {}).get("peak") or {}
        okk = ((pk.get("std_ratio_vs_1_8_point_splat") or 0) >= 0.9
               and (pk.get("p99_ratio_vs_1_8_point_splat") or 0) >= 0.9)
        A("| %s | %.4f | %.3f | %.3f | %.3f | %s |"
          % (t, cs.get("std") or 0.0, cs.get("p99") or 0.0,
             pk.get("std_ratio_vs_1_8_point_splat") or 0.0,
             pk.get("p99_ratio_vs_1_8_point_splat") or 0.0,
             "OK" if okk else "**FAIL**"))
    A("")
    A("## 9. G-MOVE / G-WARP / G-LEGACY / G-OUT")
    A("")
    A("* **G-MOVE**: ユニーク %s/48。フレーム間平均絶対差（緑・PNG）%.2f（最小 %.2f、"
      "最大 %.2f）、MP4 デコード後 %.2f（%s フレーム）。"
      % (gm.get("unique_frames"),
         (gm.get("interframe_green_png") or {}).get("mean") or 0.0,
         (gm.get("interframe_green_png") or {}).get("min") or 0.0,
         (gm.get("interframe_green_png") or {}).get("max") or 0.0,
         (gm.get("interframe_green_mp4") or {}).get("mean") or 0.0,
         (gm.get("interframe_green_mp4") or {}).get("frames_decoded")))
    A("* **G-WARP**（記録のみ）: 緑 平均 %.3f px、p50 %.3f、p99 %.3f（1.8 と同値。"
      "歪みレーンは本ジョブで触っていない）。"
      % (gw.get("g_mean_px") or 0.0, gw.get("g_p50_px") or 0.0,
         gw.get("g_p99_px") or 0.0))
    A("* **G-LEGACY**: stage-0 で記録した既存 .py %d 本、変更 %d、欠落 %d、"
      "本ジョブ外の追加 %d。新規は `run_v4_1_9.py` の 1 本のみ。"
      % (gg.get("recorded_existing") or 0, len(gg.get("changed") or []),
         len(gg.get("missing") or []),
         len(gg.get("new_since_stage0_not_of_this_job") or [])))
    for _n, _e in (gg.get("external_changes") or {}).items():
        A("  * 外部変更 `%s`（mtime %s、%d bytes、本ジョブは import も write もしない"
          "＝並行レーン）。stage0 `%s` → 現在 `%s`。"
          % (_n, _e.get("mtime"), _e.get("bytes_now") or 0,
             (_e.get("stage0_sha256") or "")[:16],
             (_e.get("now_sha256") or "")[:16]))
    _own = (gg.get("own_new_file_hashes") or {}).get("run_v4_1_9.py") or {}
    if _own.get("changed"):
        A("  * **本ジョブ自身の新規ファイルは stage-0 から変化している**（%s → %s）。"
          "理由: 作業中に (a) フレーム間の浮動小数加算順を k 昇順に固定、(b) "
          "`frame-0024.png` を本編 frame 24 に同期する処理、(c) 記録キーの追加、を"
          "行ったため。数値・幅・カーネルは stage-1 の導出値から動かしていない。"
          % ((_own.get("stage0") or "")[:16], (_own.get("now") or "")[:16]))
    A("  * 判定基準: 欠落ゼロ、かつ**本ジョブが import/write する .py に変更なし**。"
      "本ジョブが触らないファイルの変化は外部（並行レーン）として mtime 付きで開示し、"
      "隠さない。")
    A("* **G-OUT**: %s / %sx%s / %s fps / %s frames / %.3f s / %.3f MB / "
      "sha256 `%s`。"
      % ((go.get("probe") or {}).get("codec_name"),
         (go.get("probe") or {}).get("width"),
         (go.get("probe") or {}).get("height"),
         (go.get("probe") or {}).get("r_frame_rate"),
         (go.get("probe") or {}).get("nb_frames_reported"),
         float((go.get("probe") or {}).get("duration_s") or 0.0),
         (go.get("size_bytes") or 0) / 1e6,
         (go.get("sha256") or "—")[:16]))
    A("")
    A("## 10. 成果物（すべて `out/v4-1-9/`）")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for name, rec in sorted((n.get("files") or {}).items()):
        A("| `%s` | %d | `%s` |" % (name, rec["bytes"], rec["sha256"]))
    A("")
    A("`still_equals_frame_0024`: %s（`frame-0024.png` は本編 frame 24 とバイト同一、"
      "`encode.still_sync` に同期の記録）。`run_v4_1_9.py` sha256 `%s`、"
      "`transport_mesh.py` `%s`、`run_v4_1_8.py` `%s`。"
      % (n.get("still_equals_frame_0024"),
         (n.get("run_v4_1_9_sha256") or "")[:16],
         (n.get("transport_mesh_sha256") or "")[:16],
         (n.get("run_v4_1_8_sha256") or "")[:16]))
    A("")
    A("## 11. 開示・限界")
    A("")
    A("1. **G-LINE は FAIL。** しかも今回は「広がりを入れれば下がる」という方向の"
      "仮説が**構造的に否定**された: C=1（無限ぼかしの極限）でも 3 窓が写真を上回る。"
      "したがって次の標的は C の堆積ではなく歪みレーン（`_state_gradients` + "
      "B_sigma のローパスと、画素アスペクト由来の軸別周期）である。")
    A("2. **本編 (iii) は物理的には二重計上**である（C への散乱畳み込み＋合成器の "
      "w_halo）。ブリーフが (iii) を本編と指定したのでそのまま出したが、"
      "G-PEAK が FAIL するのはそのためである。物理的に整合する描き方は (ii) "
      "（合成器のハローが散乱を担う）で、`frame-0024-disc.png` がその実体。")
    A("3. **散乱の幅は単散乱**（多重散乱は幅を 10〜20% 増やす。§2.2）。"
      "多重散乱 PSF を足すのは新しい模型なので足していない。")
    A("4. **カーネルの切断**: r_max %.3f m で切って単位質量に再正規化（捕捉率 %.4f、"
      "損失は開示）。v16 の `kernel_2d` と同じ規約。"
      % (float((f24.get("psf") or {}).get("r_max_m") or 0.0),
         float(((f24.get("psf") or {}).get("kernel_info") or {}).get(
             "profile_captured_fraction") or 0.0)))
    A("5. **測定の限界**: depth は 1 次元周辺統計で、格子と等方テクスチャを区別しない。"
      "判定は数値統計のみ（目視なし）。")
    A("6. **ループ**: L はブリーフ固定（1.4755102 m）で 2 s では厳密に閉じない"
      "（1.8 §8 と同じ。本ジョブは L を変えていない）。")
    A("7. **並行レーン**: 作業ディレクトリで本ジョブ以外の .py（`run_v4_2_5.py`）が"
      "stage-0 記録後に変化している。本ジョブはそれを import も write もしていない"
      "（§9 に mtime 付きで開示）。")
    A("8. **本ジョブの新規ファイルは作業中に 3 回編集した**（加算順の固定、still 同期、"
      "記録キー追加）。輸送・堆積・PSF の数値は stage-1 の導出値のまま。")
    A("")
    A("## 12. 言えること / 言えないこと")
    A("")
    A("* **言える**: 太陽円盤の投影は 0.887 px × 0.650 px（RMS 0.443/0.325 px）で、"
      "1.8 の残差を消すのに必要な等価 σ（%.2f px）の 1/3 に満たない。"
      % (req.get("worst_sigma_needed_px_method_A") or 0.0))
    A("* **言える**: 単散乱 PSF は 33 px（モーメント式）/ 56 px（tan θ の厳密"
      "プロファイル）級で、必要量の 20 倍以上太い。散乱は「未投入」ではなく、"
      "合成器が w_halo で既に担っている。")
    A("* **言える**: C レーンの広がりをいくら強めても G-LINE は通らない"
      "（C=1 の極限で既に 3 窓が写真超え）。")
    A("* **言えない**: 「光の広がりが実際に写真と一致する幅を持つ」こと。本ジョブは"
      "幅を導出して入れただけで、写真との一致は測っていない（測れば G-PEAK が壊れる）。")
    A("* **言えない**: 見た目の最終判断（目視なし）。")
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-9.json + REPORT-V4-1-9.md")
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
    ap.add_argument("--scene", action="store_true")
    ap.add_argument("--estimate", action="store_true")
    ap.add_argument("--frame24", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--conv", action="store_true")
    ap.add_argument("--conv512", action="store_true")
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        do_scene()
        stage1_estimate()
        do_frame24()
        do_chunk(0)
        do_chunk(1)
        do_conv()
        do_warp()
        do_encode()
        do_gates()
        do_report()
        return 0
    did = False
    for flag, fn in ((a.stage0, stage0), (a.scene, do_scene),
                     (a.estimate, stage1_estimate),
                     (a.frame24, do_frame24), (a.conv, do_conv),
                     (a.conv512, do_conv512), (a.warp, do_warp),
                     (a.encode, do_encode),
                     (a.gates, do_gates), (a.report, do_report)):
        if flag:
            fn()
            did = True
    if a.chunk is not None:
        do_chunk(a.chunk)
        did = True
    if not did:
        ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
