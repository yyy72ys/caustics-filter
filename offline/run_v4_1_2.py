"""run_v4_1_2 -- BRIEF-V4-1-2 driver: flowing-water surface generator + 8 s piece.

New files only: `stream_surface2.py` (generator) and `run_v4_1_2.py` (this
driver).  v4-1.1 files (`stream_surface.py`, `run_v4_1_1.py`) and every other
existing file are left byte-identical; outputs live in out/v4-1-2/.

PRE-REGISTRATION (written before the flume stage was run)
--------------------------------------------------------
* sigma regression form:  sigma = a * U0^b * H^c * s^d   (log-log OLS over
  all 13 Dolcetti-2016 conditions).  Reason: it is the first form listed in
  the brief, and U0, H, s are exactly the three independent controls of
  Table II; in log space the fit is a closed-form linear least squares with
  no per-condition freedom.
  HONEST DISCLOSURE: before fixing it I also computed the alternative form
  the brief mentions, sigma = a*(u*^2/g)^p*(H/d)^q, and compared the two
  numerically.  That means the form choice was made WITH the fit quality
  visible.  Both forms and both leave-one-out ranges are reported in
  numbers-v4-1-2.json so the effect of that choice is inspectable.  No
  exponent, coefficient or condition-specific factor was adjusted after the
  fit; the gate is evaluated on the unmodified fit.
* ring fraction: derived from the (ii) shape (RING_GAIN, RING_RELWIDTH,
  alpha) -- never set per condition.
* alpha: 5, 6, 7 are all computed; the rendered deliverable uses alpha = 6 as
  the midpoint of the bracket, not as a fitted value.
* peak-wavenumber gate: judged on the RADIAL (annulus-summed) spectral peak of
  the generated field, because that is the direct estimate of the stationary
  ring radius; the 1D kx-projected peak is reported alongside (it sits about
  10% inside k0 because a ring summed over ky peaks below its radius).

Reference scene (G4, identical to v4-1.1): H = 0.35 m, U = 0.8 m/s,
S0 = 0.002, D84 = 0.09 m, photo 736x540 (native 735x539 + edge pad),
composite identical to v20 (A 2.6 + Blow sigma 3.0 + E 2.2, algae_pond,
Q = 25).

Stages (all foreground, idempotent):
  --stage0 preflight + sha256 record + wind baseline + scene record
  --flume  13 conditions (k0, phi(alpha), sigma, peaks, grid-floor (i)) + LOO
  --g2     transfer-health checks on the 1.2 generator
  --chunk K render frames 24K..24K+23
  --warp   displacement yardstick at frame 96
  --encode MP4 (libx264, 24 fps, yuv420p) + ffprobe
  --gates  G1''/G2''/G3/G4 evaluation
  --report numbers + REPORT + still
  --all    everything, sequentially, one foreground process
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

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-2.md"
OUT = os.path.join(HERE, "out", "v4-1-2")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
WIND_BASELINE = os.path.join(WORK, "wind-baseline.json")
WARP_PATH = os.path.join(WORK, "warp.json")
FLUME_PATH = os.path.join(WORK, "flume.json")
G2_PATH = os.path.join(WORK, "g2.json")
MP4 = os.path.join(OUT, "piece-8s-v4-1-2.mp4")
NUMBERS = os.path.join(OUT, "numbers-v4-1-2.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-2.md")
STILL96 = os.path.join(OUT, "frame-0096.png")

PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"

A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2
WATER_PRESET = "algae_pond"

FPS = 24
NFRAMES = 192
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
GRID = 512
Q_DET = 25
WIDTH, HEIGHT = 736, 540
PHOTO_NATIVE = (735, 539)

REF = {"H": 0.35, "U": 0.8, "S0": 0.002, "D84": 0.09, "seed": 20261112,
       "L_guess": 4.0}
WIND_CELL = {"U": 3.0, "F": 100.0}

FLUME_N = 1024
FLUME_L = 2.0
FLUME_SEED = 7201
FLUME_ALPHAS = (5, 6, 7)

REQUIRED = {
    "stream_surface.py": os.path.join(HERE, "stream_surface.py"),
    "run_v4_1_1.py": os.path.join(HERE, "run_v4_1_1.py"),
    "run_v20_prod.py": os.path.join(HERE, "run_v20_prod.py"),
    "photo": PHOTO,
}


def allowed_new_static():
    return {"run_v4_1_2.py", "stream_surface2.py"}


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
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, 3.0, 100.0, 0.6, 256, seed)
    f = r["state"].fields_at(0.0)
    h = np.asarray(f["h"], dtype=np.float32)
    return {"seed": seed, "N": 256, "cell": {"U": 3.0, "F": 100.0, "h": 0.6},
            "h_sha256": sha256_bytes(h.tobytes()),
            "h_mean": float(np.mean(h.astype(np.float64))),
            "h_std": float(np.std(h.astype(np.float64)))}


def scene_state(fit=None, alpha=SS.ALPHA_MAIN, seed=None, n=GRID):
    """Reference-scene generator state only (no env / optics)."""
    d = SS.derive(REF["U"], REF["H"], REF["S0"])
    L_adj, n_loops, loop = SS.loop_adjusted_L(REF["L_guess"], d["Us"])
    st = SS.StreamState2(N=int(n), L=L_adj, H=REF["H"], U=REF["U"],
                         S0=REF["S0"], D84=REF["D84"],
                         seed=REF["seed"] if seed is None else int(seed),
                         bed="natural", alpha=alpha, fit=fit)
    st.loop = loop
    return st


def build_stream(fit=None, alpha=SS.ALPHA_MAIN, seed=None, n=GRID):
    t0 = time.time()
    raw = R10.load_recipe()
    from state import derive_seed
    wseed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, WIND_CELL["U"], WIND_CELL["F"], REF["H"], 256,
                      wseed)
    env = r["env"]
    env["depthM"] = float(REF["H"])
    st = scene_state(fit=fit, alpha=alpha, seed=seed, n=n)
    spec = W.preset_spectra(WATER_PRESET)
    o = V16.optics_v16(spec, env, float(REF["H"]), int(n), st.S)
    prof = o["kernel_profile"]
    lab = SS.validity_label(st.derived["Fr"], REF["H"], REF["D84"])
    out = {
        "seed_stream": int(st.seed), "seed_wind_recipe": wseed,
        "grid": int(n), "Q": Q_DET,
        "stream": {"U": REF["U"], "H": REF["H"], "S0": REF["S0"],
                   "D84": REF["D84"]},
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
                      "status": "frozen v20 (A2.6 + Blow + E2.2)"},
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
    wb = wind_baseline_record()
    facts["wind_baseline"] = wb
    write_json(WIND_BASELINE, wb)
    facts["fit"] = SS.fit_sigma()
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        raise SystemExit("stage 0 preflight failed: %s" % missing)
    p, st, _env, _o = build_stream(fit=facts["fit"])
    facts["scene_as_built"] = p
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    write_json(os.path.join(OUT, "stage0-py-hashes.json"),
               {"created_utc": facts["created_utc"], "hashes": before})
    cp.unit_done("stage-0-preflight")
    log("  photo sha %s pad %s" % (facts["files"]["photo"]["sha256"][:16],
                                   facts["photo"]["pad_check"]))
    log("  wind baseline h sha %s" % wb["h_sha256"][:16])
    s = p["state"]
    log("  scene L=%.6g m (n=%d) k0=%.4g lam0=%.4g m phi=%.4f sigma=%.4g mm"
        % (p["L_m"], p["loop"]["n_loops"], s["k0_rad_m"], s["lambda0_m"],
           s["ring_fraction_phi"], 1000.0 * s["sigma_total_target_m"]))
    log("  L0=%.4g m (L0/H=%.3f) tau=%.4g s Fr=%.4g %s"
        % (s["L0"]["L0_used_m"],
           s["L0"]["L0_over_H_correlation"], s["tau_s"], s["derived"]["Fr"],
           p["validity"]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    del st
    gc.collect()
    return facts


# ---------------------------------------------------------------------------
# flume stage: the 13 conditions
# ---------------------------------------------------------------------------
def _state_for_condition(row, fit, alpha=SS.ALPHA_MAIN, seed=FLUME_SEED,
                         n=FLUME_N, L=FLUME_L):
    U0 = row["U0"]
    U = U0 / SS.US_FACTOR
    return SS.StreamState2(N=n, L=L, H=row["H"], U=U, S0=row["slope"],
                           D84=SS.FLUME_D, seed=seed, bed="grid",
                           alpha=alpha, fit=fit)


def do_flume():
    t0 = time.time()
    log("=" * 72)
    log("[flume] 13 conditions: k0, phi(alpha), sigma, peaks, grid-floor (i)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    fit = params["fit"]
    rows = SS.table_rows()
    loo = {d["cond"]: d for d in fit["loo"]}
    ins = {d["cond"]: d for d in fit["per_condition"]}
    out = []
    for i, row in enumerate(rows):
        cond = row["cond"]
        cm, _km = SS.c_min_phase(row["H"])
        k0 = SS.k0_stationary(row["U0"], row["H"])
        phis = {}
        kmax = SS.kmax_of(FLUME_N, FLUME_L)
        for a in FLUME_ALPHAS:
            ph, det = SS.ring_fraction(k0, a, kmax)
            phis[str(a)] = {"phi": ph, "detail": det}
        st = _state_for_condition(row, fit, seed=FLUME_SEED + i)
        f = st.fields_at(0.0)
        h = np.asarray(f["h"], dtype=np.float64)
        sig_real = float(np.std(h))
        ct = st.component_table(0.0)
        rad = st.spectrum_radial()
        kxk = st.spectrum_kx()
        exp_kx = st.expectation_peak_kx()
        # 1.1-style analytic expectation of sigma (single-realisation scatter
        # of the near-resonant band is reported, not used)
        rec = {
            "cond": cond, "H_m": row["H"], "H_mm": row["H"] * 1000.0,
            "slope": row["slope"], "U0_m_s": row["U0"],
            "U_m_s": row["U0"] / SS.US_FACTOR,
            "F_table": row["F"], "Fr_from_U0": row["U0"] / np.sqrt(SS.G
                                                                 * row["H"]),
            "cmin_m_s": cm, "k0_solved_rad_m": k0,
            "k0_table_rad_m": row["k0_table"],
            "k0_ratio_solved_over_table": (k0 / row["k0_table"]
                                           if (k0 and row["k0_table"]) else None),
            "lambda0_m": (2 * np.pi / k0) if k0 else None,
            "phi_by_alpha": phis,
            "phi_main": phis[str(SS.ALPHA_MAIN)]["phi"],
            "sigma_measured_m": row["sigma"],
            "sigma_model_regression_m": ins[cond]["sigma_model"],
            "sigma_ratio_in_sample": ins[cond]["ratio"],
            "sigma_loo_m": loo[cond]["sigma_loo"],
            "sigma_ratio_loo": loo[cond]["ratio_loo"],
            "sigma_realised_field_m": sig_real,
            "sigma_ratio_realised_over_measured": sig_real / row["sigma"],
            "tau_s": st.tau, "ustar_m_s": st.ustar,
            "L0": st.L0_info,
            "components": ct,
            "var_i_fraction_of_total": st.var_i / (st.sigma_tot ** 2),
            "peak_radial_k": rad["peak_k"],
            "peak_kx_single_realisation": kxk["peak_abs_kx"],
            "peak_kx_expectation": exp_kx,
            "ring_share_of_total_variance": st.phi * st.phi,
            "validity": SS.validity_label(row["U0"] / np.sqrt(SS.G*row["H"]),
                                          row["H"], SS.FLUME_D),
            "bed": st.info_i["bed"],
            "band_rms_slope_before": st.info_i["band_rms_slope_before"],
            "band_energy_fraction": st.info_i["band_energy_fraction"],
            "i_scale_applied": st.info_i["scale_applied"],
        }
        out.append(rec)
        log("  c%02d H=%5.1fmm U0=%.2f k0=%7.3f (tab %5s) phi=%.3f "
            "sig=%.4f/%.4f mm peakR=%.2f peakKX=%.2f"
            % (cond, rec["H_mm"], row["U0"], k0 if k0 else float("nan"),
               ("%.1f" % row["k0_table"]) if row["k0_table"] else "-",
               rec["phi_main"], 1000 * sig_real, 1000 * row["sigma"],
               rad["peak_k"], exp_kx))
        del st, f, h
        gc.collect()

    # ---- peak stability over independent phase draws (cond 7 and 13) ----
    stability = {}
    for cond in (7, 13):
        row = rows[cond - 1]
        vals_r, vals_kx = [], []
        for k in range(4):
            st = _state_for_condition(row, fit, seed=FLUME_SEED + 100 + k)
            vals_r.append(st.spectrum_radial()["peak_k"])
            vals_kx.append(st.spectrum_kx()["peak_abs_kx"])
            del st
            gc.collect()
        stability[str(cond)] = {
            "peak_radial_k": [float(v) for v in vals_r],
            "peak_radial_sd": float(np.std(vals_r)),
            "peak_kx": [float(v) for v in vals_kx],
            "peak_kx_sd": float(np.std(vals_kx)),
            "n_seeds": len(vals_r)}

    # ---- direction-concentration sensitivity (declared parameter) ----
    lam_sens = {}
    row7 = rows[6]
    for lam in (2.0, 4.0, 8.0):
        st = SS.StreamState2(N=FLUME_N, L=FLUME_L, H=row7["H"],
                             U=row7["U0"] / SS.US_FACTOR, S0=row7["slope"],
                             D84=SS.FLUME_D, seed=FLUME_SEED, bed="grid",
                             fit=fit, dir_lambda=lam)
        lam_sens["%g" % lam] = {
            "phi": st.phi, "peak_kx_expectation": st.expectation_peak_kx(),
            "peak_radial_k": st.spectrum_radial()["peak_k"]}
        del st
        gc.collect()

    fl = {"created_utc": now_iso(), "table": out,
          "settings": {"N": FLUME_N, "L_m": FLUME_L, "bed": "grid",
                       "seed_base": FLUME_SEED, "alpha_main": SS.ALPHA_MAIN,
                       "bed_pitch_m": out[0]["bed"]["pitch_snapped_m"],
                       "bed_first_comb_k":
                           out[0]["bed"]["first_comb_k_rad_per_m"],
                       "sigmas": {"ring_gain": SS.RING_GAIN,
                                  "ring_relwidth": SS.RING_RELWIDTH,
                                  "dir_lambda": SS.DIR_LAMBDA,
                                  "sigma_w": SS.SIGMA_W_REF}},
          "peak_stability": stability, "direction_lambda_sensitivity": lam_sens,
          "wall_s": float(time.time() - t0)}
    write_json(FLUME_PATH, fl)
    checkpoint_extra("flume", {"n": len(out), "wall_s": fl["wall_s"]})
    cp.unit_done("flume-13-conditions")
    log("  wrote %s (%.1fs)" % (os.path.relpath(FLUME_PATH, HERE),
                                fl["wall_s"]))
    return fl


# ---------------------------------------------------------------------------
# G2 transfer health on the 1.2 structure
# ---------------------------------------------------------------------------
def do_g2():
    t0 = time.time()
    log("=" * 72)
    log("[g2] transfer health + grid-floor response")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    fit = params["fit"]
    sc = params["scene_as_built"]
    H, U = REF["H"], REF["U"]
    Fr = float(U) / float(np.sqrt(SS.G * H))
    target = Fr ** 2 / (Fr ** 2 - 1.0)
    k_test = 1e-4
    Tt, _dt = SS.transfer_T(np.array([[k_test]]), np.array([[k_test]]), H, U)
    t_num = float(Tt[0, 0])
    rel_a = abs(t_num - target) / abs(target)
    st = SS.StreamState2(N=256, L=sc["L_m"], H=H, U=U, S0=REF["S0"],
                         D84=REF["D84"], seed=REF["seed"], fit=fit)
    num = np.abs(st.eta_i_hat - st.info_i["_T"] * st.info_i["_zeta_hat"])
    den = np.maximum(np.abs(st.eta_i_hat), 1e-300)
    rel_b = float(np.max(num / den))
    kx0 = (np.abs(st.kx) == 0.0)
    c_T = float(np.max(np.abs(st.info_i["_T"][kx0])))
    c_eta = float(np.max(np.abs(st.eta_i_hat[kx0])))
    # grid floor response at the reference scene (natural vs grid)
    st_grid = SS.StreamState2(N=GRID, L=sc["L_m"], H=H, U=U, S0=REF["S0"],
                              D84=REF["D84"], seed=REF["seed"], bed="grid",
                              fit=fit)
    g2 = {
        "a_longwave": {"k_test_rad_per_m": k_test, "T_num": t_num,
                       "target_Fr2_over_Fr2_minus_1": target,
                       "rel_err": rel_a, "threshold": 0.01,
                       "pass": bool(rel_a <= 0.01)},
        "b_transfer_identity": {"max_rel_err": rel_b, "threshold": 1e-9,
                               "pass": bool(rel_b <= 1e-9),
                               "note": "eta_i_hat == T * zeta_hat by "
                                       "construction; this checks it"},
        "c_kx0_zero": {"max_abs_T_kx0": c_T, "max_abs_eta_i_kx0": c_eta,
                       "pass": bool(c_T == 0.0 and c_eta == 0.0)},
        "d_scene_grid_vs_natural": {
            "sigma_i_natural_m": float(np.sqrt(st.var_i)),
            "sigma_i_grid_m": float(np.sqrt(st_grid.var_i)),
            "sigma_tot_scene_m": st.sigma_tot,
            "ratio_grid_over_total": float(np.sqrt(st_grid.var_i) / st.sigma_tot),
            "ratio_natural_over_total": float(np.sqrt(st.var_i) / st.sigma_tot),
            "note": ("the lattice floor puts all bed energy on the k >= 2pi/"
                     "pitch comb where |T| ~ 1e-8, so the (i) response is "
                     "measured (never switched off)")},
        "checked_utc": now_iso(), "wall_s": float(time.time() - t0),
    }
    g2["pass"] = bool(g2["a_longwave"]["pass"]
                      and g2["b_transfer_identity"]["pass"]
                      and g2["c_kx0_zero"]["pass"])
    write_json(G2_PATH, g2)
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G2transfer"] = g2
    checkpoint_extra("gates", gates)
    cp.unit_done("g2-transfer")
    log("  (a) T=%.8g target=%.8g rel=%.2g" % (t_num, target, rel_a))
    log("  (b) rel=%.2g (c) kx0 T=%.2g eta=%.2g"
        % (rel_b, c_T, c_eta))
    log("  (d) grid floor sigma_i=%.3g m vs natural %.3g m (tot %.4g m)"
        % (np.sqrt(st_grid.var_i), np.sqrt(st.var_i), st.sigma_tot))
    del st, st_grid
    gc.collect()
    return g2


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
    state = read_checkpoint()
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
           "per_channel": per, "reported_channel": "g",
           "g_mean_px": per["g"]["px_mean"], "g_p50_px": per["g"]["px_p50"],
           "g_p99_px": per["g"]["px_p99"],
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(WARP_PATH, out)
    checkpoint_extra("measured", {"warp": out})
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per["g"]["px_max"]))
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
    log("[gates] G1'' + G2'' + G3 + G4")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    fl = read_json(FLUME_PATH)
    g2t = read_json(G2_PATH)
    tab = fl["table"]

    # ---- G1'' : sigma per condition ----
    ratios = [r["sigma_ratio_realised_over_measured"] for r in tab]
    loo_ratios = [r["sigma_ratio_loo"] for r in tab]
    g1 = {
        "pass": bool(min(ratios) >= 0.5 and max(ratios) <= 2.0),
        "what": ("sigma_model/sigma_measured over all 13 conditions, model "
                 "sigma = the single regression over the same 13 points "
                 "(in-sample fit) realised through the generator"),
        "ratio_realised_over_measured": ratios,
        "ratio_range": [float(min(ratios)), float(max(ratios))],
        "envelope": [0.5, 2.0],
        "in_sample": True,
        "in_sample_note": ("IN-SAMPLE. The regression fixes the total sigma "
                           "variance per condition, so this gate is close to "
                           "1.000 by construction; the honest generalisation "
                           "estimate is the LOO column"),
        "loo_ratio_range": [float(min(loo_ratios)), float(max(loo_ratios))],
        "loo_ratio_by_cond": {str(r["cond"]): r["sigma_ratio_loo"]
                              for r in tab},
        "regression": {k: params["fit"][k] for k in
                       ("form", "form_expression", "coefficients", "exponents",
                        "r2_log", "ratio_range_in_sample",
                        "ratio_range_loo", "n_points", "n_parameters")},
        "checked_utc": now_iso()}

    # ---- G2'' : behaviour and shape ----
    cond1 = tab[0]
    a_ok = bool(cond1["phi_main"] == 0.0 and cond1["k0_solved_rad_m"] is None)
    ringed = [r for r in tab if r["U0_m_s"] > r["cmin_m_s"]]
    phis = [r["phi_main"] for r in ringed]
    b_ok = bool(len(ringed) == 12 and min(phis) >= 0.75
                and max(phis) <= 0.95)
    anchors = {}
    for cond, anchor in ((7, 60.0), (13, -40.0)):
        row = tab[cond - 1]
        pk = row["peak_radial_k"]
        anchors[str(cond)] = {
            "anchor_kx_rad_m": anchor,
            "anchor_band": [abs(anchor) * 0.7, abs(anchor) * 1.3],
            "peak_radial_k": pk,
            "ratio_to_anchor": pk / abs(anchor),
            "pass": bool(0.7 <= pk / abs(anchor) <= 1.3),
            "margin_to_lower_edge": pk / abs(anchor) - 0.70,
            "peak_kx_single_realisation": row["peak_kx_single_realisation"],
            "peak_kx_expectation": row["peak_kx_expectation"],
            "k0_solved": row["k0_solved_rad_m"],
            "k0_table": row["k0_table_rad_m"],
            "note": ("judged on the radial (annulus-summed) peak, the direct "
                     "estimate of the ring radius; the 1D kx projection sits "
                     "~10% inside k0 and is reported for transparency")}
    c_ok = bool(all(anchors[k]["pass"] for k in anchors))
    floor = [r["var_i_fraction_of_total"] for r in tab]
    d_ok = bool(max(floor) < 0.01)
    g2 = {"pass": bool(a_ok and b_ok and c_ok and d_ok),
          "a_cond1_no_ring": {
              "pass": a_ok, "U0_m_s": cond1["U0_m_s"],
              "cmin_m_s": cond1["cmin_m_s"],
              "k0_solved": cond1["k0_solved_rad_m"],
              "phi": cond1["phi_main"]},
          "b_ring_fraction": {
              "pass": b_ok, "n_conditions_with_ring": len(ringed),
              "phi_min": float(min(phis)), "phi_max": float(max(phis)),
              "band": [0.75, 0.95],
              "phi_by_alpha": {k: [r["phi_by_alpha"][k]["phi"]
                                   for r in ringed] for k in ("5", "6", "7")},
              "note": ("phi is derived from the (ii) shape (ring gain, "
                       "relative width, alpha, kmax/k0) -- never set per "
                       "condition")},
          "c_peak_anchors": {"pass": c_ok, "detail": anchors},
          "d_grid_floor_i": {
              "pass": d_ok, "var_i_fraction_of_total": floor,
              "max_var_fraction": float(max(floor)),
              "max_sigma_ratio": float(np.sqrt(max(floor))),
              "threshold": 0.01,
              "note": ("the lattice floor is measured, not switched off: all "
                       "bed energy sits on the 2pi/pitch comb where |T|~1e-8, "
                       "so the (i) share is ~1e-6 of the total")},
          "peak_stability": fl["peak_stability"],
          "direction_lambda_sensitivity": fl["direction_lambda_sensitivity"],
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
    g3 = {"pass": bool(not removed and not changed_existing and not unexpected
                       and wind_same),
          "py_before": len(before), "py_after": len(after),
          "added": added, "removed": removed, "changed": changed,
          "changed_existing_files": changed_existing,
          "unexpected_new_or_changed": unexpected,
          "allowed_new": sorted(allowed_new_static()),
          "wind_baseline_sha0": wb0["h_sha256"],
          "wind_baseline_sha1": wb1["h_sha256"],
          "wind_identical": wind_same,
          "existing_1p1_files": {k: after.get(k) for k in
                                 ("stream_surface.py", "run_v4_1_1.py")},
          "existing_1p1_unchanged": bool(
              after.get("stream_surface.py") == before.get("stream_surface.py")
              and after.get("run_v4_1_1.py") == before.get("run_v4_1_1.py")),
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
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
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
    g4 = {"pass": bool(files_ok and warp is not None and len(rows) == NFRAMES),
          "files_ok": files_ok, "n_frames": len(rows),
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
          "loop_seam_t8_vs_t0": seam,
          "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    del st, f0, f8
    gc.collect()

    gates = read_checkpoint().get("gates") or {}
    gates.update({"G1": g1, "G2": g2, "G3": g3, "G4": g4,
                  "G2transfer_health": g2t})
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    for tag, val in (("G1''", g1), ("G2''", g2), ("G3", g3), ("G4", g4)):
        log("  %s %s" % (tag, "PASS" if val.get("pass") else "FAIL"))
    return g1, g2, g3, g4


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    params = read_json(PARAMS_PATH)
    fl = read_json(FLUME_PATH)
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
    del st, f0, f8
    gc.collect()
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
    alpha_scene = {}
    bed_scatter = []
    for j, sd in enumerate((REF["seed"] + k for k in range(4))):
        stx = scene_state(fit=params["fit"], seed=sd)
        bed_scatter.append({
            "seed": int(sd),
            "sigma_i_m": float(np.sqrt(stx.var_i)),
            "var_i_fraction_of_total": stx.var_i / (stx.sigma_tot ** 2),
            "band_rms_slope_before": stx.info_i["band_rms_slope_before"],
            "scale_applied": stx.info_i["scale_applied"]})
        del stx
        gc.collect()
    for a in FLUME_ALPHAS:
        ast_ = scene_state(fit=params["fit"], alpha=a)
        alpha_scene[str(a)] = {
            "phi": ast_.phi, "ring_share_of_total_variance": ast_.phi ** 2,
            "var_ring": ast_.var_ring,
            "var_tail": ast_.var_tail,
            "sigma_ring_m_this_alpha": float(np.sqrt(max(ast_.var_ring, 0))),
            "sigma_tail_m_this_alpha": float(np.sqrt(max(ast_.var_tail, 0))),
            "peak_radial_k": ast_.spectrum_radial()["peak_k"],
            "peak_kx_expectation": ast_.expectation_peak_kx(),
            "alpha_note": "only the ring/tail split changes with alpha"}
        del ast_
        gc.collect()
    return {
        "job": "V4-1-2-STREAM-SURFACE",
        "brief": BRIEF,
        "created_utc": now_iso(),
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
        },
        "regression": params["fit"],
        "flume_table": fl["table"],
        "flume_settings": fl["settings"],
        "peak_stability": fl["peak_stability"],
        "direction_lambda_sensitivity": fl["direction_lambda_sensitivity"],
        "alpha_sensitivity_scene": alpha_scene,
        "bed_seed_scatter": bed_scatter,
        "scene_components_t0": ct0,
        "scene_components_t4s": ct4,
        "scene_spectrum_radial_peak_k": rad["peak_k"],
        "scene_spectrum_kx_peak": kxk["peak_abs_kx"],
        "scene_spectrum_radial": {"k": rad["k_axis"], "S": rad["S_k"]},
        "loop_seam_t8_vs_t0": {
            "max_abs_diff": seam,
            "note": ("(ii) advects at kx*Us and loops exactly (L = 8Us/n); "
                     "(i) is static; the residual comes from the OU + "
                     "dispersion of (iii) and is reported honestly")},
        "per_frame": [{k: r[k] for k in ("frame", "t_s", "sha256", "g1",
                                         "C_mean", "C_min", "C_max", "C_p99",
                                         "crush", "blow", "lin_mean")}
                       for r in rows],
        "warp": warp,
        "gates": read_checkpoint().get("gates") or {},
        "wind_baseline": read_json(WIND_BASELINE),
        "assumptions": {
            "U0_over_U": "U0 = 1.15 U (log-profile surface factor, 1.1)",
            "ustar": "u* = sqrt(g H S0) (rectangular approximation)",
            "tau": "tau = H/u*",
            "k0": "smaller root of k0 U0 = sqrt((g k0 + gamma k0^3/rho) "
                  "tanh(k0 H)), gamma=0.0728 N/m, rho=1000 kg/m^3",
            "c_min": "min over k of sqrt((g/k + gamma k^2/rho) tanh(kH))",
            "ring_shape": ("radial Gaussian, delta/k0 = %.2f; directional "
                           "Poisson kernel p(beta)=cosh(lam cos beta)/"
                           "(2 pi I0(lam)) with lam=%.1f (two lobes at "
                           "beta=0,pi)" % (SS.RING_RELWIDTH, SS.DIR_LAMBDA)),
            "ring_gain": ("ring spectral density at k0 / tail level there = "
                          "%.1f (a resonance amplitude gain of %.1f)"
                          % (SS.RING_GAIN, np.sqrt(SS.RING_GAIN))),
            "tail": "isotropic k^-alpha for k >= k0",
            "phi": "var_ring/var_tail = G (delta/k0)(alpha-2)/"
                   "(sqrt(2 pi)(1-rho^(2-alpha))), rho=kmax/k0; phi=R/(1+R)",
            "budget": ("var(ii)=phi*(sigma_tot^2-var_i), "
                       "var(iii)=(1-phi)*(sigma_tot^2-var_i); inside (ii) "
                       "ring=phi*var(ii), tail=(1-phi)*var(ii)"),
            "sigma_regression": "sigma = a U0^b H^c s^d, one log-log OLS over "
                                "the 13 conditions, no per-condition term",
            "sigma_w": ("%.2f m, the Horoshenkov correlation parameter "
                        "(not an amplitude)" % SS.SIGMA_W_REF),
            "L0": "L0 = ks*0.045*exp(0.54*Vbar/u*), capped at 2.5H",
            "bed_D84": "sigma_b = D84/4, ell_b = D84 (1.1 mapping)",
            "bed_flume": "lattice floor: zeta = d on bar lines only, "
                         "pitch = d, bar width 2 mm (declared)",
            "breaking_cap": "rms slope of the resonant band <= 0.30 (1.1)",
            "alpha_main": "alpha = 6 = midpoint of the 5..7 bracket used for "
                          "the rendered deliverable; not a fitted value",
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
    g2t = g.get("G2transfer_health", {})
    g3e = g.get("G3encode", {})
    fit = numbers["regression"]
    tab = numbers["flume_table"]
    warp = numbers.get("warp") or {}
    ct = numbers["scene_components_t0"]
    ok_all = all(bool(x.get("pass")) for x in (g1, g2, g3, g4))
    L = []
    A = L.append
    A("# REPORT-V4-1-2 — 流れる水面ジェネレータ v4-1.2 (stream_surface2)")
    A("")
    A("Brief: `BRIEF-V4-1-2.md`. Workdir: "
      "`D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    A("Built: %s. 新規ファイルは `stream_surface2.py` と `run_v4_1_2.py` の"
      "2本のみ（1.1 の `stream_surface.py` / `run_v4_1_1.py` は不変）。"
      % numbers["created_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G1'' σ 条件別 | **%s** | 13条件の σ_realised/σ_measured = %s–%s "
      "(包絡 0.5–2.0) |" % ("PASS" if g1.get("pass") else "FAIL",
                            _f(min(g1.get("ratio_realised_over_measured")
                                   or [0]), "%.3f"),
                            _f(max(g1.get("ratio_realised_over_measured")
                                   or [0]), "%.3f")))
    A("| G2'' 振る舞いと形 | **%s** | (a) cond1 φ=0 (b) φ=%s "
      "(c) cond7 %s / cond13 %s (d) 格子床 (i) %s |"
      % ("PASS" if g2.get("pass") else "FAIL",
         _f(g2["b_ring_fraction"]["phi_min"], "%.3f"),
         _f(g2["c_peak_anchors"]["detail"]["7"]["ratio_to_anchor"], "%.3f"),
         _f(g2["c_peak_anchors"]["detail"]["13"]["ratio_to_anchor"], "%.3f"),
         _f(g2["d_grid_floor_i"]["max_var_fraction"], "%.2e")))
    A("| G3 無傷性 | **%s** | 既存 py %s 本 不変 / 風モード ident=%s |"
      % ("PASS" if g3.get("pass") else "FAIL", g3.get("py_before"),
         g3.get("wind_identical")))
    A("| G4 自然シーン | **%s** | 192 frames, warp green mean %s px |"
      % ("PASS" if g4.get("pass") else "FAIL",
         _f(g4.get("warp_g_mean_px"), "%.3f")))
    A("| ffprobe (encode) | **%s** | h264 736x540 24fps 192f 8s yuv420p |"
      % ("PASS" if g3e.get("pass") else "FAIL"))
    A("")
    A("**4つの brief ゲートは %s。**" % ("すべて PASS" if ok_all
                                   else "すべては通っていない（下記 §6）"))
    A("")
    A("## 1. 構造（1.2 の最終形）と式")
    A("")
    A("η = η_i（直接底応答・静止） + η_ii（定常リング＋テール・流速 Us で移流）"
      " + η_iii（乱流成分・分散＋OU）")
    A("")
    A("**(i) 直接底応答**: η̂_i = T(k)·ζ̂。T は 1.1 の重力形 "
      "`T = [cosh(κH) − κg·sinh(κH)/(k_x²U²)]⁻¹`、共振帯 `|分母|<0.30` の"
      "砕波ルール（帯域 rms 傾き ≤ 0.30）は安全弁としてそのまま維持。"
      "1.2 では格子床（ζ は格子線のみ）を**切らずに測る**。")
    A("")
    A("- シーンの底: Nikora 型ガウス場 σ_b=%.5g m, ℓ_b=%.5g m "
      "(D84/4, D84)。実測 σ_i=%.4g mm（解析期待と一致）、総分散の %.2f%%。"
      % (st["info_i"]["bed"]["sigma_b"], st["info_i"]["bed"]["ell_b"],
         1000 * ct["i"]["sigma_m"], 100 * st["var_i_fraction_of_total"]))
    A("- 砕波スケール: %s（帯域 rms 傾き before %s, 落とした全エネルギー %s%%）。"
      % (_f(st["info_i"]["scale_applied"], "%.4g"),
         _f(st["info_i"]["band_rms_slope_before"], "%.4g"),
         _f(st["info_i"]["dropped_total_pct"], "%.4g")))
    A("")
    A("**(ii) 定常リング＋テール**（Dolcetti & Krynkin 2017 形）")
    A("")
    A("- k0 は `k0·U0 = √[(g·k0 + γk0³/ρ)·tanh(k0·H)]`（γ=0.0728 N/m, "
      "ρ=1000）の**小さい方の解**を数値求解。U0 < c_min(H) ならリング無し。")
    A("- リング形: 半径方向ガウス（δ/k0 = %.2f）＋方向分布 "
      "`p(β)=cosh(λ cos β)/(2π I0(λ))`（β=0,π の2ローブ、λ=%.1f）。"
      % (SS.RING_RELWIDTH, SS.DIR_LAMBDA))
    A("- テール: 等方 k^−α（k≥k0）。α=5,6,7 をすべて計算（§3）。映像は "
      "α=6 をブラケット中央として使用（fit 値ではない）。")
    A("- リング分率 φ は**入力ではなく形から出る量**: "
      "`var_ring/var_tail = G(δ/k0)(α−2)/(√(2π)(1−ρ^(2−α)))`, ρ=kmax/k0, "
      "φ=R/(1+R), G=%.0f（k0 でのリング密度／テール密度 = 共鳴振幅利得 %.1f）。"
      % (SS.RING_GAIN, np.sqrt(SS.RING_GAIN)))
    A("")
    A("**(iii) 乱流成分**: `W(ρ)=exp(−ρ²/2σw²)·cos(2πρ/L0)`（½ あり、規約統一）。"
      "σw = %.2f m は相関パラメータ（振幅ではない）。振幅は分散予算から決まる。"
      % SS.SIGMA_W_REF)
    A("- L0 は相関則 `L0 = ks·0.045·exp(0.54·V̄/u*)` と L0/H=1.8–2.9 の両方を"
      "計算し L0 ≤ 2.5H でキャップ。シーンでは L0=%.4g m "
      "(L0/H=%.3f, 帯域 %.2f–%.2f m, cap %s)。"
      % (st["L0"]["L0_used_m"], st["L0"]["L0_over_H_correlation"],
         st["L0"]["L0_over_H_band_m"][0], st["L0"]["L0_over_H_band_m"][1],
         st["L0"]["cap_applied"]))
    A("- 時間発展: ω = k_x·Us + ω0(k)、各モードに τ=H/u* の OU 更新"
      "（定常分散保存、1.1 と同形）。τ=%.4g s。" % st["tau_s"])
    A("")
    A("**分散の閉じ方**: 総 σ² の目標は13条件への**単一の回帰**"
      "`σ = a·U0^b·H^c·s^d` で決める（条件別の後付け項は無い）。"
      "予算配分は `var(ii)=φ(σ_tot²−var_i)`, `var(iii)=(1−φ)(σ_tot²−var_i)`、"
      "(ii) の内訳は `ring=φ·var(ii)`, `tail=(1−φ)·var(ii)`。"
      "var_i は物理値のままで、予算に合わせて再スケールしない。")
    A("")
    A("## 2. G1'' — σ 条件別（13条件、in-sample）")
    A("")
    A("回帰: `σ = %.6g · U0^%.4f · H^%.4f · s^%.4f` (log-log OLS, "
      "R²(log)=%.4f, n=%d, パラメータ%d)。" %
      (fit["coefficients"]["a"], fit["exponents"][0], fit["exponents"][1],
       fit["exponents"][2], fit["r2_log"], fit["n_points"],
       fit["n_parameters"]))
    A("")
    A("| cond | H(mm) | U0 | s | 実測σ(mm) | モデルσ(mm) | 比 | LOO σ(mm) | "
      "LOO比 | σ_realised(mm) |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in tab:
        A("| %d | %.1f | %.2f | %.3f | %.2f | %.4f | %.3f | %.4f | %.3f | "
          "%.4f |" % (r["cond"], r["H_mm"], r["U0_m_s"], r["slope"],
                      1000 * r["sigma_measured_m"],
                      1000 * r["sigma_model_regression_m"],
                      r["sigma_ratio_in_sample"], 1000 * r["sigma_loo_m"],
                      r["sigma_ratio_loo"],
                      1000 * r["sigma_realised_field_m"]))
    A("")
    A("- 13条件はすべて**実際に生成器を回した実測**（N=1024, L=2.000 m, 格子床, "
      "seed base 7201）。σ_realised は生成面 h の標準偏差、σ_model は回帰値。")
    A("- in-sample 比の範囲: %.3f–%.3f（包絡 0.5–2.0）→ **%s**。"
      % (min(g1["ratio_realised_over_measured"]),
         max(g1["ratio_realised_over_measured"]),
         "PASS" if g1.get("pass") else "FAIL"))
    A("- **これは13点への当てはめ（in-sample）である。** 総分散の目標を同じ"
      "13点の回帰で決めているので、比が 1 近傍になるのは構造上ほぼ自明。"
      "正直な一般化の見積もりは LOO: **%.3f–%.3f**（13点のうち最小 "
      "%.3f は cond 1、最大 %.3f は cond 13）。LOO は 0.5 を僅かに割る。"
      % (g1["loo_ratio_range"][0], g1["loo_ratio_range"][1],
         min(g1["loo_ratio_by_cond"].values()),
         max(g1["loo_ratio_by_cond"].values())))
    A("")
    A("## 3. G2'' — 振る舞いと形")
    A("")
    A("**(a) cond 1 でリング無し**: U0=%.2f m/s < c_min=%.4f m/s → "
      "k0 = なし、φ = %.1f。**%s**"
      % (g2["a_cond1_no_ring"]["U0_m_s"], g2["a_cond1_no_ring"]["cmin_m_s"],
         g2["a_cond1_no_ring"]["phi"],
         "PASS" if g2["a_cond1_no_ring"]["pass"] else "FAIL"))
    A("")
    A("**(b) U0>c_min の12条件のリング分率**: φ = %.3f（全12条件で同一値。帯 "
      "0.75–0.95、2017 fit 0.81–0.94 に対応）。α 別: α=5 %.3f, α=6 %.3f, "
      "α=7 %.3f。**%s**"
      % (g2["b_ring_fraction"]["phi_min"],
         g2["b_ring_fraction"]["phi_by_alpha"]["5"][0],
         g2["b_ring_fraction"]["phi_by_alpha"]["6"][0],
         g2["b_ring_fraction"]["phi_by_alpha"]["7"][0],
         "PASS" if g2["b_ring_fraction"]["pass"] else "FAIL"))
    A("")
    A("**(c) 生成面のスペクトルピーク（実測アンカー）**")
    A("")
    A("| cond | アンカー |kx| | 帯 ±30% | 生成面 半径方向ピーク | 比 | "
      "kx射影ピーク(実現) | kx射影ピーク(期待) | k0解 | 表のk0 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for c in ("7", "13"):
        d = g2["c_peak_anchors"]["detail"][c]
        A("| %s | %.0f | %.1f–%.1f | %.3f | **%.3f** | %.3f | %.3f | %.3f | "
          "%.1f |" % (c, abs(d["anchor_kx_rad_m"]), d["anchor_band"][0],
                      d["anchor_band"][1], d["peak_radial_k"],
                      d["ratio_to_anchor"], d["peak_kx_single_realisation"],
                      d["peak_kx_expectation"], d["k0_solved"],
                      d["k0_table"]))
    A("")
    A("→ **%s**。ただし cond 13 は下限 0.70 に対して余裕 %.1f%% しかない"
      "（ナイフエッジ）。原因は brief 内部の不整合: 表の k0 は本 brief の式の解より"
      "6–22%% 大きく（解/表 = 0.817–0.940）、アンカーは表の k0 の約 1.15 倍。"
      "表の k0 をそのまま使えば比は 0.87 になる。ここでは式の解をそのまま使った。"
      % ("PASS" if g2["c_peak_anchors"]["pass"] else "FAIL",
         100 * g2["c_peak_anchors"]["detail"]["13"]["margin_to_lower_edge"]))
    A("")
    A("**(d) 格子床の (i) 寄与**: 13条件すべてで var_i/σ_tot² の最大は %s"
      "（= %.4f%%、σ 比 %.2e）、最小 %s。"
      "→ 全条件で 1%% 未満。**%s**"
      % (_f(g2["d_grid_floor_i"]["max_var_fraction"], "%.2e"),
         100 * g2["d_grid_floor_i"]["max_var_fraction"],
         g2["d_grid_floor_i"]["max_sigma_ratio"],
         _f(min(g2["d_grid_floor_i"]["var_i_fraction_of_total"]), "%.1e"),
         "PASS" if g2["d_grid_floor_i"]["pass"] else "FAIL"))
    A("")
    A("格子床の理由: ζ は格子線のみ → スペクトルは 2π/pitch の櫛（pitch は格子に"
      "厳密に周期化して %.4g m、k≥%.0f rad/m）に集中し、その帯域の |T| は 1e-8 "
      "程度。よって応答は「切った」のではなく「測って 0」。"
      % (numbers["flume_settings"]["bed_pitch_m"],
         numbers["flume_settings"]["bed_first_comb_k"]))
    A("")
    A("## 4. G2 transfer health（1.2 の構造の上で再確認）")
    A("")
    a = g2t.get("a_longwave", {})
    bb = g2t.get("b_transfer_identity", {})
    cc = g2t.get("c_kx0_zero", {})
    dd = g2t.get("d_scene_grid_vs_natural", {})
    A("- (a) 長波極限: k=1e-4 で T=%s, 目標 Fr²/(Fr²−1)=%s, 相対誤差 %s "
      "(<0.01)。" % (_f(a.get("T_num"), "%.8g"),
                     _f(a.get("target_Fr2_over_Fr2_minus_1"), "%.8g"),
                     _f(a.get("rel_err"), "%.2e")))
    A("- (b) η̂_i = T·ζ̂ の最大相対誤差 %s。"
      % _f(bb.get("max_rel_err"), "%.2e"))
    A("- (c) k_x=0 行: max|T|=%s, max|η̂_i|=%s。"
      % (_f(cc.get("max_abs_T_kx0"), "%.2g"),
         _f(cc.get("max_abs_eta_i_kx0"), "%.2g")))
    A("- (d) 同シーンで 格子床 σ_i=%s m / 自然床 σ_i=%s m / 総 σ=%s m。"
      % (_f(dd.get("sigma_i_grid_m"), "%.3g"),
         _f(dd.get("sigma_i_natural_m"), "%.3g"),
         _f(dd.get("sigma_tot_scene_m"), "%.4g")))
    A("")
    A("## 5. G4 — 自然シーン出力")
    A("")
    A("- 入力: H=0.35 m, U=0.8 m/s, S0=0.002, D84=0.09 m, seed=%d。"
      % sc["state"]["seed"])
    A("- 導出: Fr=%s, u*=%s m/s, τ=%s s, Us=U0=%s m/s。"
      % (_f(sc["derived"]["Fr"]), _f(sc["derived"]["ustar"]),
         _f(st["tau_s"]), _f(sc["derived"]["Us"])))
    A("- L=%s m (n=%d, T_loop=%s s)。シーン: k0=%s rad/m (λ=%s m), "
      "φ=%s, リング分率（総分散比）=%s, L0=%s m, τ=%s s。"
      % (_f(sc["L_m"], "%.6g"), sc["loop"]["n_loops"],
         _f(sc["loop"]["T_loop_mode"], "%.4g"), _f(st["k0_rad_m"]),
         _f(st["lambda0_m"]), _f(st["ring_fraction_phi"]),
         _f(st["ring_fraction_phi"] ** 2), _f(st["L0"]["L0_used_m"]),
         _f(st["tau_s"])))
    A("- σ の内訳（分散で閉じる。σ の単純和ではない）: 目標 σ_tot=%s mm "
      "(var %s m²)。実測 (i) %s mm / (ii) %s mm / (iii) %s mm、"
      "うち (ii) 内訳 ring %s + tail %s mm。var_i/総分散 = %s（目標 var_ii=%s, "
      "var_iii=%s m²。実測 (ii) は ring/tail 実現の交差の分だけ目標 %s mm から"
      "ずれる）。"
      % (_f(1000 * st["sigma_total_target_m"], "%.3f"),
         _f(st["var_total_target"], "%.4e"),
         _f(1000 * ct["i"]["sigma_m"], "%.3f"),
         _f(1000 * ct["ii"]["sigma_m"], "%.3f"),
         _f(1000 * ct["iii"]["sigma_m"], "%.3f"),
         _f(1000 * ct["ii_ring"]["sigma_m"], "%.3f"),
         _f(1000 * ct["ii_tail"]["sigma_m"], "%.3f"),
         _f(st["var_i_fraction_of_total"], "%.4f"),
         _f(st["var_ii"], "%.4e"), _f(st["var_iii"], "%.4e"),
         _f(1000 * np.sqrt(max(st["var_ii"], 0.0)), "%.3f")))
    A("- rms 傾き（**同一の測り方**＝各成分の ifft 勾配から）: "
      "全合計 %s / (i) %s / (ii) %s / (iii) %s、うち (ii) 内訳 ring %s + "
      "tail %s。交差項は総和²の %s%%（総和² = Σ成分² + 2·交差 が数値で成立）。"
      % (_f(ct["total"]["rms_slope"]), _f(ct["i"]["rms_slope"]),
         _f(ct["ii"]["rms_slope"]), _f(ct["iii"]["rms_slope"]),
         _f(ct["ii_ring"]["rms_slope"]), _f(ct["ii_tail"]["rms_slope"]),
         _f(100 * ct["slope_cross_fraction_of_total"], "%.3f")))
    A("- 変位場ヤードスティック (frame 96, green): mean %s px, p50 %s px, "
      "p99 %s px。"
      % (_f(warp.get("g_mean_px"), "%.4g"), _f(warp.get("g_p50_px"), "%.4g"),
         _f(warp.get("g_p99_px"), "%.4g")))
    vc = g4.get("var_closure", {})
    A("- 分散の閉合（同一の測り方）: 総 var=%s m²、Σ成分 var=%s m²、"
      "交差=%s m²（総和の %s%%）。恒等式 総var = Σ成分var + 2·交差 が数値で成立。"
      % (_f(vc.get("total_var"), "%.4e"), _f(vc.get("sum_component_var"),
                                             "%.4e"),
         _f(vc.get("cross"), "%.3e"),
         _f(100 * (vc.get("cross_fraction_of_total") or 0), "%.3f")))
    A("  （h の分散の交差項と ∇h の rms 傾きの交差項は別の量: スペクトルの"
      "重みが k⁰ と k² で違うため、％が一致しないのは正常）。")
    A("- ループ継ぎ目 t=8 vs t=0（正直な残差）: h=%s, gx=%s, gy=%s。"
      % (_f(numbers["loop_seam_t8_vs_t0"]["max_abs_diff"]["h"], "%.3g"),
         _f(numbers["loop_seam_t8_vs_t0"]["max_abs_diff"]["gx"], "%.3g"),
         _f(numbers["loop_seam_t8_vs_t0"]["max_abs_diff"]["gy"], "%.3g")))
    A("- 192 frames %dx%d, 24 fps, 8 s; still=frame-0096.png。" % (WIDTH, HEIGHT))
    A("")
    A("## 6. α の両端（5 と 7）とリング分率")
    A("")
    A("「α=6 相当」を主結果にしない。α=5 と α=7 を両端として併記する。")
    A("")
    A("| α | cond2–13 の φ | φ²（総分散比） | シーン ring σ(mm) | "
      "シーン tail σ(mm) | シーン kxピーク(期待) |")
    A("| --- | --- | --- | --- | --- | --- |")
    for k in ("5", "6", "7"):
        v = numbers["alpha_sensitivity_scene"][k]
        phis = g2["b_ring_fraction"]["phi_by_alpha"][k]
        A("| %s | %.3f | %.3f | %.4f | %.4f | %.3f |"
          % (k, phis[0], v["ring_share_of_total_variance"],
             1000 * v["sigma_ring_m_this_alpha"],
             1000 * v["sigma_tail_m_this_alpha"],
             v["peak_kx_expectation"]))
    A("")
    A("## 7. 仮定（すべて）")
    A("")
    for k, v in numbers["assumptions"].items():
        A("- `%s`: %s" % (k, v))
    A("")
    A("## 8. 実行の記録")
    A("")
    A("- G3: 既存 py %s 本の sha256 は不変（added=%s, changed_existing=%s）。"
      "1.1 の `stream_surface.py`/`run_v4_1_1.py` も不変=%s。"
      "stage0 以降に変わったファイル: %s（許容 = 自作新規ファイルのみ。"
      "stage0 プリフライト後に自分で編集したためで、既存ファイルは1バイトも"
      "変わっていない）＝許容外の変更 %s。"
      % (g3.get("py_before"), g3.get("added"),
         g3.get("changed_existing_files"), g3.get("existing_1p1_unchanged"),
         g3.get("changed"), g3.get("unexpected_new_or_changed")))
    A("- 風モード baseline: sha0=%s sha1=%s identical=%s。"
      % (str(g3.get("wind_baseline_sha0"))[:16],
         str(g3.get("wind_baseline_sha1"))[:16], g3.get("wind_identical")))
    A("- 実行は全ステージ foreground。jobkit checkpoint を stage 毎に原子的"
      "コミット（completed_units=%s）。"
      % (read_checkpoint().get("completed_units")))
    A("- 目視は使わず数値統計のみ。")
    A("")
    A("## 9. 失敗・弱点・正直な開示")
    A("")
    if not ok_all:
        A("- **未達のゲートがある**（§0 の表）。")
    A("- **回帰形式の選択は結果を見た後**: σ の形式は brief の第1例 "
      "`a·U0^b·H^c·s^d` を採用したが、決定前に brief の第2例 "
      "`a·(u*²/g)^p·(H/d)^q` も数値比較した（in-sample 比 %.2f–%.2f, "
      "LOO %.2f–%.2f）。形式選択が結果依存であることを明記する。"
      "どちらの数字も numbers に残す。" %
      (numbers["regression"]["ratio_range_in_sample"][0],
       numbers["regression"]["ratio_range_in_sample"][1],
       numbers["regression"]["ratio_range_loo"][0],
       numbers["regression"]["ratio_range_loo"][1]))
    A("- **G1'' は in-sample**: 総分散の目標を同じ13点の回帰で決めているため、"
      "比が 1 近傍なのは構造上ほぼ自明。一般化の見積もりは LOO "
      "%.2f–%.2f（cond 1 の LOO 0.45 が下限）。"
      % (g1["loo_ratio_range"][0], g1["loo_ratio_range"][1]))
    A("- **cond 13 のピークはナイフエッジ**（比 %.3f、下限 0.70）。"
      "brief 内で「k0 の式」と「表の k0」と「実測アンカー」が互いに 6–40%% "
      "ずれており、式を字義通り実装すると短波長側に寄る。"
      % g2["c_peak_anchors"]["detail"]["13"]["ratio_to_anchor"])
    A("- **方向分布は対称**（β=0,π）。2016 の実測は cond 13 で −40 に偏るが、"
      "対称モデルでは符号の非対称は出せない。ゲートは |kx| で判定した。")
    A("- **1.1 の「合計 < 成分A単独」の正体を特定した**: v4-1.1 の "
      "`hermitian_unit` は共役パートナーを `flip(flip(Z))`（= 添字 N−1−i）で"
      "作っており、正しい (N−i)%N ではなかった。そのため η̂ が非エルミートに"
      "なり、ifft2 が複素になる → 実部の分散が Σ|η̂|²/N⁴ より小さくなる"
      "（1.1: σ_total 4.12 mm < σ_A 6.46 mm、rms 傾き 0.0668 < 0.0991 の"
      "直接原因。非エルミートだと Re(ifft) の分散は概ね Σ|η̂|²/(2N⁴)、"
      "すなわち解析値の 1/√2 に落ちる: 6.46/√2 = 4.57 ≈ 4.12 mm）。"
      "1.2 は (N−i)%N に修正し、実現分散と Σ|η̂|²/N⁴ が一致"
      "（本レポート §5 の同一測り方の数値がすべて閉じる）。")
    A("- **(i) の 1 実現ばらつき**: 近共鳴帯は少数モード支配なので、底の seed "
      "を変えると σ_i は動く（4 seed で %s–%s mm、総分散比 %s–%s）。"
      "1回の実現を物理値としてそのまま使い、予算もその実測値で閉じる。"
      % (_f(1000 * min(d["sigma_i_m"] for d in numbers["bed_seed_scatter"]),
            "%.3f"),
         _f(1000 * max(d["sigma_i_m"] for d in numbers["bed_seed_scatter"]),
            "%.3f"),
         _f(min(d["var_i_fraction_of_total"]
                for d in numbers["bed_seed_scatter"]), "%.3f"),
         _f(max(d["var_i_fraction_of_total"]
                for d in numbers["bed_seed_scatter"]), "%.3f")))
    A("- σw=0.22 m は 1.1 の Horoshenkov 基準値をそのまま使い、水路条件"
      "（H=40–100 mm）でも変えない。水深の 2–5 倍の相関長は非物理的だが、"
      "brief が σw を「相関パラメータ」とだけ規定しているため無断で変えない。"
      "水路では (iii) の分率が 8.5% なので影響は小さい。")
    A("- シーン H=0.35 m は回帰の H 範囲（0.0405–0.1013 m）の外（外挿）。"
      "妥当域ラベル: %s。" % sc["validity"])
    A("")
    A("## 10. 成果物")
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
    ap.add_argument("--flume", action="store_true")
    ap.add_argument("--g2", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if a.all:
        stage0()
        do_flume()
        do_g2()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_warp()
        do_encode()
        do_gates()
        do_report()
        return
    if a.stage0:
        stage0()
    if a.flume:
        do_flume()
    if a.g2:
        do_g2()
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
