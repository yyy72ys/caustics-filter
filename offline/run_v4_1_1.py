"""run_v4_1_1 -- BRIEF-V4-1-1: flowing-water surface generator + 8 s piece.

New files only: stream_surface.py (generator), run_v4_1_1.py (this driver).
Existing files are never modified. Outputs live in out/v4-1-1/.

Reference scene (G4): H = 0.35 m, U = 0.8 m/s, S0 = 0.002 (stated value,
  gives u* = sqrt(gHS0) = 0.0829 m/s), D84 = 0.09 m, photo 736x540,
  composite identical to v20 (A2.6 + Blow sigma 3.0 + E2.2, v16 water
  a/b partition + FF phase + FF-CDF halo, algae_pond, Q = 25).

Transport/composite reuse: det.det_frame (Q = 25) + run_v16_water.optics_v16
  (algae_pond) + run_v16_water.render_frame_v16 with the frozen v20
  composite constants. The water SURFACE comes from stream_surface instead
  of the wind spectrum. Environment (sun/turbidity/weather) reuses the v20
  preset-1 weather (U_wind = 3 m/s, F = 100 m); stream hydraulics
  (U, H, S0) drive only the surface generator. depthM = H.

Stages (all foreground, idempotent):
  --stage0  preflight + sha256 record + wind baseline + scene record
  --g1      G1 flume reproduction (Dolcetti 2016 representative points)
  --g2      G2 transfer health
  --chunk K render frames 24K..24K+23
  --warp    displacement yardstick at frame 96
  --encode  MP4 (libx264, 24 fps, yuv420p) + ffprobe check
  --gates   G3 intactness + G4 output
  --report  numbers + REPORT-V4-1-1.md + still
  --all     everything above, sequentially, in one foreground process
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

import run_v10 as R10  # noqa: E402
import run_v15 as R15  # noqa: E402
import v12_water as W  # noqa: E402
import run_v16_water as V16  # noqa: E402
import stream_surface as SS  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-1.md"
OUT = os.path.join(HERE, "out", "v4-1-1")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
WIND_BASELINE = os.path.join(WORK, "wind-baseline.json")
WARP_PATH = os.path.join(WORK, "warp.json")
MP4 = os.path.join(OUT, "piece-8s-v4-1-1.mp4")
NUMBERS = os.path.join(OUT, "numbers-v4-1-1.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-1.md")
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

G1_POINTS = [{"H": 0.0405, "U": 0.19}, {"H": 0.07, "U": 0.40},
             {"H": 0.1013, "U": 0.60}]
G1_D = 0.0254
G1_S0 = 0.001
G1_N = 256
G1_L = 2.0
G1_SEED = 7001


def allowed_new_static():
    return {"run_v4_1_1.py", "stream_surface.py"}


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


def checkpoint_extra(key, value=None):
    state = read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}
    if value is not None:
        state[key] = value
        state["updated_at"] = now_iso()
        atomic_write_bytes(CP_PATH, json.dumps(
            state, ensure_ascii=False, indent=2).encode("utf-8"))
    return state


def read_checkpoint():
    return read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}


def load_photo():
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"u8": padded, "lin": lin}


def wind_baseline_record():
    """Wind-mode baseline for G3: v20-recipe env + wind state hash (N=256)."""
    raw = R10.load_recipe()
    from state import derive_seed
    seed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, 3.0, 100.0, 0.6, 256, seed)
    st = r["state"]
    f = st.fields_at(0.0)
    h = np.asarray(f["h"], dtype=np.float32)
    return {"seed": seed, "N": 256, "cell": {"U": 3.0, "F": 100.0, "h": 0.6},
            "h_sha256": sha256_bytes(h.tobytes()),
            "h_mean": float(np.mean(h.astype(np.float64))),
            "h_std": float(np.std(h.astype(np.float64)))}


def build_stream():
    """Reference stream scene: state + env + optics (identical every process)."""
    t0 = time.time()
    raw = R10.load_recipe()
    from state import derive_seed
    wseed = int(derive_seed(raw, 1))
    # env reuses v20 preset-1 weather; depth overridden to stream H below
    r = R15.scene_for(raw, WIND_CELL["U"], WIND_CELL["F"], REF["H"], 256,
                      wseed)
    env = r["env"]
    env["depthM"] = float(REF["H"])
    d = SS.derive(REF["U"], REF["H"], REF["S0"])
    L_adj, n_loops, loop = SS.loop_adjusted_L(REF["L_guess"], d["Us"])
    st = SS.StreamState(N=GRID, L=L_adj, H=REF["H"], U=REF["U"], S0=REF["S0"],
                        D84=REF["D84"], seed=REF["seed"])
    spec = W.preset_spectra(WATER_PRESET)
    o = V16.optics_v16(spec, env, float(REF["H"]), GRID, st.S)
    prof = o["kernel_profile"]
    out = {
        "seed_stream": REF["seed"], "seed_wind_recipe": wseed,
        "grid": GRID, "Q": Q_DET,
        "stream": {"U": REF["U"], "H": REF["H"], "S0": REF["S0"],
                   "D84": REF["D84"]},
        "derived": d,
        "loop": loop,
        "L_m": float(st.S),
        "patch_m": float(st.S),
        "depth_m": float(env["depthM"]),
        "bed_mapping": st.bed_mapping,
        "sigma_b": st.sigma_b, "ell_b": st.ell_b,
        "sigma_c_assumption": ("sigma_c = u*^2/g (velocity head); brief "
                               "gives no C amplitude, stated assumption"),
        "sigma_c": st.sigma_c,
        "sigma_w": st.sigma_w, "L0": st.L0,
        "tau_s": st.tau,
        "validity": SS.validity_label(st.Fr, REF["H"], REF["D84"]),
        "Fr": st.Fr,
        "breaking": {"band_def": st.Ainfo["res_band_def"],
                     "slope_cap": st.Ainfo["slope_cap"],
                     "scale_applied": st.Ainfo["scale_applied"],
                     "band_rms_before": st.Ainfo["band_rms_before"],
                     "band_energy_fraction":
                         st.Ainfo["band_energy_fraction"],
                     "dropped_band_pct": st.Ainfo["dropped_band_pct"],
                     "dropped_total_pct": st.Ainfo["dropped_total_pct"]},
        "turbid": float(env["turbid"]),
        "clarity": float(env["clarity"]),
        "sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
        "water_preset": spec["preset"],
        "water_label": spec.get("label"),
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
                                         if prof else None),
        },
        "build_wall_s": float(time.time() - t0),
    }
    return out, st, env, o


def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "created_utc": now_iso(),
             "python": sys.version.split()[0]}
    missing = []
    code = {}
    for name in ("run_v20_prod.py", "run_v16_water.py", "run_v19_prod.py"):
        p = os.path.join(HERE, name)
        ok = os.path.isfile(p)
        code[name] = {"present": ok, "sha256": sha256_file(p) if ok else None}
        if not ok:
            missing.append(name)
    facts["code"] = code
    facts["v16_defaults"] = {
        "PHASE_DEFAULT": V16.PHASE_DEFAULT,
        "KERNEL_DEFAULT": V16.KERNEL_DEFAULT,
        "assert_phase_ff": V16.PHASE_DEFAULT == "ff",
        "assert_kernel_ffcdf": V16.KERNEL_DEFAULT == "ff-cdf"}
    if not (V16.PHASE_DEFAULT == "ff" and V16.KERNEL_DEFAULT == "ff-cdf"):
        missing.append("v16 defaults not ff / ff-cdf")
    from PIL import Image
    if not os.path.isfile(PHOTO):
        missing.append("base photo missing: %s" % PHOTO)
        facts["photo"] = {"path": PHOTO, "present": False}
    else:
        im = Image.open(PHOTO)
        im.load()
        facts["photo"] = {"path": PHOTO, "present": True,
                          "bytes": os.path.getsize(PHOTO),
                          "native_px": [im.size[0], im.size[1]],
                          "mode": im.mode, "format": im.format,
                          "sha256": sha256_file(PHOTO),
                          "assert_735x539": tuple(im.size) == PHOTO_NATIVE,
                          "output_px": [WIDTH, HEIGHT]}
        if tuple(im.size) != PHOTO_NATIVE:
            missing.append("photo decodes at %s, not 735x539" % (im.size,))
    if os.path.isfile(PHOTO):
        p = load_photo()
        n = np.asarray(Image.open(PHOTO).convert("RGB"))
        facts["photo"]["pad_check"] = {
            "shape": list(p["u8"].shape),
            "right_col_ok": bool(np.array_equal(p["u8"][:-1, -1, :],
                                                n[:, -1, :])),
            "bottom_row_ok": bool(np.array_equal(p["u8"][-1, :-1, :],
                                                 n[-1, :, :])),
            "corner_ok": bool(np.array_equal(p["u8"][-1, -1, :],
                                             n[-1, -1, :]))}
        if not all(facts["photo"]["pad_check"][k]
                   for k in ("right_col_ok", "bottom_row_ok", "corner_ok")):
            missing.append("pad violates the v2 edge-replication rule")
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before), "hashes": before}
    wb = wind_baseline_record()
    facts["wind_baseline"] = wb
    write_json(WIND_BASELINE, wb)
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        raise SystemExit("stage 0 preflight failed: %s" % missing)
    p, _st, _env, _o = build_stream()
    facts["scene_as_built"] = p
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    write_json(os.path.join(OUT, "stage0-py-hashes.json"),
               {"created_utc": facts["created_utc"], "hashes": before})
    cp.unit_done("stage-0-preflight")
    log("  photo sha %s pad %s" % (facts["photo"]["sha256"][:16],
                                   facts["photo"]["pad_check"]))
    log("  wind baseline h sha %s" % wb["h_sha256"][:16])
    log("  stream L=%.5g m (n=%d) tau=%.4g s Fr=%.4g label=%s" % (
        p["L_m"], p["loop"]["n_loops"], p["tau_s"], p["Fr"], p["validity"]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


def _k0_of_state(st):
    eh = st.etaA_hat + st.A0
    N2 = st.N * st.N
    pw = np.abs(eh) ** 2
    pw[st.Cgrid["kappa"] == 0] = 0.0
    j = int(np.argmax(pw))
    return float(st.Cgrid["kappa"].ravel()[j])


def do_g1():
    t0 = time.time()
    log("=" * 72)
    log("[G1] Dolcetti-2016 flume reproduction (representative points)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    rows = []
    for i, pt in enumerate(G1_POINTS):
        m = SS.bed_mapping_sphere(G1_D)
        st = SS.StreamState(N=G1_N, L=G1_L, H=pt["H"], U=pt["U"], S0=G1_S0,
                            D84=0.0, sigma_b=m["sigma_b"], ell_b=m["ell_b"],
                            seed=G1_SEED + i)
        st.bed_mapping = dict(m, note="sphere bed d=25.4 mm mapping")
        f = st.fields_at(0.0)
        h = np.asarray(f["h"], dtype=np.float64)
        sigma_mm = float(np.std(h) * 1000.0)
        k0 = _k0_of_state(st)
        d = SS.derive(pt["U"], pt["H"], G1_S0)
        rows.append({
            "H_mm": pt["H"] * 1000.0, "U": pt["U"], "S0_assumed": G1_S0,
            "ustar": d["ustar"], "tau_s": d["tau"], "Fr": d["Fr"],
            "bed": {"d_mm": G1_D * 1000.0, **m},
            "sigma_eta_mm": sigma_mm,
            "sigma_envelope_mm": [0.05, 2.03],
            "sigma_pass": bool(0.05 <= sigma_mm <= 2.03),
            "k0_rad_per_m": k0,
            "k0_envelope": [31, 132],
            "k0_pass": bool(31 <= k0 <= 132),
            "validity": SS.validity_label(d["Fr"], pt["H"], G1_D),
        })
        log("  H=%.1fmm U=%.2f Fr=%.3f sigma=%.4gmm k0=%.4g %s" % (
            pt["H"] * 1000, pt["U"], d["Fr"], sigma_mm, k0,
            rows[-1]["validity"]))
        del st, f, h
        gc.collect()
    g1 = {"pass": bool(all(r["sigma_pass"] and r["k0_pass"] for r in rows)),
          "points": rows,
          "note": ("representative points only (brief allows); bed mapping "
                   "sigma_b=d/4, ell_b=d is a stated assumption; S0=0.001 "
                   "for tau/Us only (flume slope not in brief); no post-hoc "
                   "tuning -- FAIL is recorded honestly"),
          "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G1"] = g1
    checkpoint_extra("gates", gates)
    cp.unit_done("gate-G1-flume")
    log("  G1 %s" % ("PASS" if g1["pass"] else "FAIL"))
    return g1


def do_g2():
    t0 = time.time()
    log("=" * 72)
    log("[G2] transfer health")
    makedirs()
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
    sc = params["scene_as_built"]
    H, U = REF["H"], REF["U"]
    Fr = float(U) / float(np.sqrt(SS.G * H))
    target = Fr ** 2 / (Fr ** 2 - 1.0)
    # (a) long-wave numeric limit at k -> 0 along kx (ky = 0)
    k_test = 1e-4
    Tt, _dt = SS.transfer_T(np.array([[k_test]]), np.array([[k_test]]), H, U)
    t_num = float(Tt[0, 0])
    rel_a = abs(t_num - target) / abs(target)
    # grid smallest-mode value (informational: finite-kH correction)
    kx, ky, kap, dk = SS.wavenumbers(GRID, sc["L_m"])
    Tgrid, _dg = SS.transfer_T(kx, kap, H, U)
    m = (np.abs(ky) < 1e-12) & (np.abs(kx) > 0)
    kmin = float(np.min(np.abs(kx[m])))
    sel = (np.abs(np.abs(kx) - kmin) < 1e-12) & (np.abs(ky) < 1e-12)
    t_grid = float(np.mean(Tgrid[sel]))
    # (b) S_etaeta = |T|^2 S_bed identity: etaA_hat == T*zeta_hat exactly
    st = SS.StreamState(N=256, L=sc["L_m"], H=H, U=U, S0=REF["S0"],
                        D84=REF["D84"], seed=REF["seed"])
    prod = st.Ainfo["T"] * st.Ainfo["zeta_hat"] if "T" in st.Ainfo else None
    num = np.abs(st.etaA_hat - prod)
    den = np.maximum(np.abs(prod), 1e-300)
    rel_b = float(np.max(num / den))
    # (c) kx = 0 row is exactly 0
    kx0 = (np.abs(st.Cgrid["kx"]) == 0.0)
    c_T = float(np.max(np.abs(st.Ainfo["T"][kx0])))
    c_eta = float(np.max(np.abs(st.etaA_hat[kx0])))
    g2 = {
        "pass": bool(rel_a <= 0.01 and rel_b <= 1e-9 and c_T == 0.0
                      and c_eta == 0.0),
        "a_longwave": {"k_test_rad_per_m": k_test, "T_num": t_num,
                       "target_Fr2_over_Fr2_minus_1": target,
                       "rel_err": rel_a, "threshold": 0.01,
                       "grid_kmin_mode_T": t_grid,
                       "grid_kmin_note": ("finite-kH correction at the grid "
                                          "mode; the gate is on the k->0 "
                                          "numeric limit")},
        "b_transfer_identity": {
            "max_rel_err_etaA_vs_T_times_zeta": rel_b, "threshold": 1e-9,
            "meaning": ("S_etaeta = |T|^2 S_bed holds by construction; "
                        "per-mode power fluctuates as chi2 around it")},
        "c_kx0_zero": {"max_abs_T_kx0": c_T, "max_abs_etaA_kx0": c_eta},
        "checked_utc": now_iso(), "wall_s": float(time.time() - t0),
    }
    del st
    gc.collect()
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G2"] = g2
    checkpoint_extra("gates", gates)
    cp.unit_done("gate-G2-transfer")
    log("  (a) T=%.6g target=%.6g rel=%.2g" % (t_num, target, rel_a))
    log("  (b) identity rel=%.2g | (c) kx0 T=%.2g eta=%.2g" % (rel_b, c_T,
                                                              c_eta))
    log("  G2 %s" % ("PASS" if g2["pass"] else "FAIL"))
    return g2


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
    state = read_checkpoint()
    if chunk_done(k, state):
        log("  chunk %02d committed -- skip" % k)
        return True
    _p, st, env, optics = build_stream()
    photo = load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_one(i, st, env, optics, photo)
        sha = save_png(chunk_path(i), r["still"])
        m = r["metrics"]
        c = r["C"]
        gx = None
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
            "completed_utc": now_iso(), "wall_s": float(time.time() - t_start)}
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    del st, optics, photo
    gc.collect()
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
    _p, st, env, optics = build_stream()
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
    out = {"frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement, sampled at output pixel "
                     "centres, S=%.6g m" % S),
           "patch_m": S, "depth_m": depth, "width": WIDTH, "height": HEIGHT,
           "per_channel": per_ch, "reported_channel": "g",
           "g_mean_px": per_ch["g"]["px_mean"],
           "g_p50_px": per_ch["g"]["px_p50"],
           "g_p99_px": per_ch["g"]["px_p99"],
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(WARP_PATH, out)
    checkpoint_extra("measured", {"warp": out})
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           per_ch["g"]["px_max"]))
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
    missingf = [i for i in range(NFRAMES) if not os.path.isfile(pattern % i)]
    if missingf:
        raise SystemExit("frames missing: %d" % len(missingf))
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
          "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    state = read_checkpoint()
    gates = state.get("gates") or {}
    gates["G3encode"] = g3
    checkpoint_extra("gates", gates)
    cp.unit_done("encode")
    log("  G3encode(ffprobe) %s" % ("PASS" if g3["pass"] else "FAIL"))
    return g3


def do_gates():
    t0 = time.time()
    log("=" * 72)
    log("[gates] G3 intactness + G4 output")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS_PATH)
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
          "checked_utc": now_iso()}
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
        os.path.join(FRAMES, "frame-%04d.png" % 96), MP4))
    rms_finite = bool(np.all(np.isfinite([r["C_mean"] for r in rows])))
    _p, st, _env, _o = build_stream()
    rs = st.realised_stats()
    g4 = {"pass": bool(files_ok and rms_finite and warp is not None
                       and len(rows) == NFRAMES),
          "files_ok": files_ok, "n_frames": len(rows),
          "warp_g_mean_px": (warp or {}).get("g_mean_px"),
          "rms_slope_total_t0": rs["rms_slope"],
          "rms_slope_A": rs["rms_slope_A"],
          "rms_slope_C_t0": rs["rms_slope_C_t0"],
          "sigma_eta_A_m": rs["sigma_eta_A"],
          "sigma_eta_total_t0_m": rs["sigma_eta_total_t0"],
          "checked_utc": now_iso(), "wall_s": float(time.time() - t0)}
    del st
    gc.collect()
    gates = state.get("gates") or {}
    gates["G3"] = g3
    gates["G4"] = g4
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-G3-G4")
    log("  G3 %s (wind identical: %s)" % ("PASS" if g3["pass"] else "FAIL",
                                          wind_same))
    log("  G4 %s" % ("PASS" if g4["pass"] else "FAIL"))
    return g3, g4


def build_numbers():
    params = read_json(PARAMS_PATH)
    state = read_checkpoint()
    chunks = state.get("chunks") or {}
    gates = state.get("gates") or {}
    rows = []
    for k in range(NCHUNKS):
        rows.extend(chunks[str(k)]["metrics"])
    rows.sort(key=lambda r: r["frame"])
    warp = read_json(WARP_PATH) if os.path.exists(WARP_PATH) else None
    _p, st, _env, _o = build_stream()
    rs = st.realised_stats()
    # loop seam residual (honest: OU + dispersive part do not loop)
    f0 = st.fields_at(0.0)
    f8 = st.fields_at(8.0)
    seam = {k: float(np.max(np.abs(
        np.asarray(f8[k], np.float64) - np.asarray(f0[k], np.float64))))
        for k in ("h", "gx", "gy")}
    del st, f0, f8
    gc.collect()
    return {
        "job": "V4-1-1-STREAM-SURFACE",
        "brief": BRIEF,
        "created_utc": now_iso(),
        "scene": params["scene_as_built"],
        "realised": rs,
        "loop_seam_t8_vs_t0": {
            "max_abs_diff": seam,
            "note": ("advective part loops by construction (L = 8Us/n); "
                     "intrinsic-dispersive + OU parts do not; residual "
                     "reported honestly")},
        "per_frame": [{k: r[k] for k in ("frame", "t_s", "sha256", "g1",
                                         "C_mean", "C_min", "C_max", "C_p99",
                                         "crush", "blow", "lin_mean")}
                       for r in rows],
        "warp": warp,
        "gates": gates,
        "wind_baseline": read_json(WIND_BASELINE),
        "deliverables": {
            "mp4": os.path.relpath(MP4, HERE),
            "still_0096": os.path.relpath(STILL96, HERE),
            "report": os.path.relpath(REPORT, HERE),
            "checkpoint": os.path.relpath(CP_PATH, HERE),
        },
    }


def write_report(numbers):
    sc = numbers["scene"]
    g = numbers["gates"]
    per = numbers["per_frame"]
    warp = numbers.get("warp") or {}
    g1 = g.get("G1", {})
    g2 = g.get("G2", {})
    g3 = g.get("G3", {})
    g4 = g.get("G4", {})
    g3e = g.get("G3encode", {})
    verdict = all(bool(g.get(k, {}).get("pass"))
                  for k in ("G1", "G2", "G3", "G4"))
    L = []
    L.append("# REPORT-V4-1-1 — 流れる水面ジェネレータ (stream_surface)")
    L.append("")
    L.append("Brief: `BRIEF-V4-1-1.md`. Workdir: "
             "`D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    L.append("Built: %s. New files: `stream_surface.py`, `run_v4_1_1.py` "
             "only." % numbers["created_utc"])
    L.append("")
    L.append("## 0. 式と仮定（設計 v2 の再掲）")
    L.append("")
    L.append("- Fr = U/sqrt(gH), κ=|k|, u* = sqrt(g·H·S0)（矩形近似）。")
    L.append("- T(k) = [cosh(κH) − κg·sinh(κH)/(k_x²U²)]⁻¹。"
              "共振は (k_x·U)² = gκ·tanh(κH)。表面張力は T に入れない。")
    L.append("- k_x = 0 の行は T = 0。共振帯の定義: |分母| < %.3g。"
             "砕波ルール: 帯域の rms 傾き寄与が 0.30 を超えたら帯域を一様"
             "スケールして 0.30 に合わせる。" % SS.RES_BAND_DENOM)
    L.append("- D84 写像: sigma_b = D84/4, ell_b = D84。球床写像: "
             "sigma_b = d/4, ell_b = d。いずれも明記された仮定。")
    L.append("- 成分Bは独立に加えない。近共振帯を診断数値として切り分けて"
             "報告するのみ（二重計上なし）。")
    L.append("- 成分C: W(ρ)=exp(−ρ²/2σw²)·cos(2πρ/L0), σw=0.22 m, "
             "L0=0.17 m は物理 m としてそのまま使用（無断リスケールなし。"
             "L0/dx<8 または L<6σw なら異常終了）。振幅は brief に規定が"
             "ないため sigma_c = u*²/g（速度水頭）を仮定として置いた。")
    L.append("- 時間発展: ω = k_x·Us + ω0(k), Us = 1.15·U（対数則の表面係数"
             "の仮定）, ω0(k)=√((gk+σk³/ρ)·tanh(kH)), σ=0.0728, ρ=1000。"
             "各モード振幅に相関時間 τ=H/u* の OU 更新（定常分散保存）。"
             "τ = %.5g s。" % sc["tau_s"])
    L.append("- ループ整合: n = round(8·Us/L_guess), L = 8·Us/n。"
             "L_guess=%.3g → L=%.6g m (n=%d, T_loop=%.6g s)。"
             % (sc["loop"]["L_guess"], sc["L_m"], sc["loop"]["n_loops"],
                sc["loop"]["T_loop_mode"]))
    L.append("- 妥当域ラベル: %s (Fr=%.4g, H/D84=%.4g)。"
             % (sc["validity"], sc["Fr"], REF["H"] / REF["D84"]))
    L.append("")
    L.append("## 1. Verdict")
    L.append("")
    L.append("| gate | result |")
    L.append("| --- | --- |")
    L.append("| G1 flume | **%s** |" % ("PASS" if g1.get("pass") else "FAIL"))
    L.append("| G2 transfer | **%s** |"
             % ("PASS" if g2.get("pass") else "FAIL"))
    L.append("| G3 intactness | **%s** |"
             % ("PASS" if g3.get("pass") else "FAIL"))
    L.append("| G4 output | **%s** |" % ("PASS" if g4.get("pass") else "FAIL"))
    L.append("| ffprobe (encode) | **%s** |"
             % ("PASS" if g3e.get("pass") else "FAIL"))
    L.append("")
    L.append("All four brief gates %s."
             % ("PASS" if verdict else "DID NOT ALL PASS -- see below"))
    L.append("")
    L.append("## 2. G1 水路再現 (Dolcetti 2016, 球床 d=25.4 mm)")
    L.append("")
    L.append("代表点のみ（brief が代表点でよいと明記）。S0=0.001 は τ/Us "
             "用の仮定（brief に水路勾配の指定なし）。後付け調整なし。")
    L.append("")
    L.append("| H(mm) | U(m/s) | Fr | σ_η(mm) | 包絡 | k0(rad/m) | 包絡 | "
             "label |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in (g1.get("points") or []):
        L.append("| %.1f | %.2f | %.3f | %.4g | 0.05–2.03 | %.4g | 31–132 "
                 "| %s |" % (r["H_mm"], r["U"], r["Fr"], r["sigma_eta_mm"],
                              r["k0_rad_per_m"], r["validity"]))
    L.append("")
    L.append("## 3. G2 伝達の健全性")
    L.append("")
    a = g2.get("a_longwave", {})
    b = g2.get("b_transfer_identity", {})
    c = g2.get("c_kx0_zero", {})
    L.append("- (a) 長波数値極限: k=%.0e で T=%.8g, 目標 Fr²/(Fr²−1)=%.8g, "
             "相対誤差 %.3g (閾値 0.01)。"
             % (a.get("k_test_rad_per_m", 0), a.get("T_num", 0),
                a.get("target_Fr2_over_Fr2_minus_1", 0),
                a.get("rel_err", 9)))
    L.append("- (b) S_ηη=|T|²S_bed: η̂_A − T·ζ̂ の最大相対誤差 %.3g "
             "(閾値 1e-9)。" % b.get("max_rel_err_etaA_vs_T_times_zeta", 9))
    L.append("- (c) k_x=0 行: max|T|=%.3g, max|η̂|=%.3g。"
             % (c.get("max_abs_T_kx0", 9), c.get("max_abs_etaA_kx0", 9)))
    L.append("")
    L.append("## 4. G3 無傷性")
    L.append("")
    L.append("- 既存 py: before=%s after=%s, removed=%s, "
             "changed_existing=%s, unexpected=%s。"
             % (g3.get("py_before"), g3.get("py_after"), g3.get("removed"),
                g3.get("changed_existing_files"),
                g3.get("unexpected_new_or_changed")))
    L.append("- 風モード (run_v20_prod 系, N=256 baseline): sha0=%s sha1=%s "
             "identical=%s。"
             % (str(g3.get("wind_baseline_sha0"))[:16],
                str(g3.get("wind_baseline_sha1"))[:16],
                g3.get("wind_identical")))
    L.append("")
    L.append("## 5. G4 出力（基準シーン H=0.35, U=0.8, S0=0.002, D84=0.09）")
    L.append("")
    L.append("- u* = sqrt(9.81·0.35·0.002) = %.6g m/s, τ = H/u* = %.5g s, "
             "Us = 1.15U = %.4g m/s。"
             % (sc["derived"]["ustar"], sc["tau_s"], sc["derived"]["Us"]))
    L.append("- D84 写像: sigma_b=%.5g m, ell_b=%.5g m。"
             % (sc["sigma_b"], sc["ell_b"]))
    L.append("- 砕波スケール: %.4g (band rms before %.4g, band energy "
             "fraction %.4g, 落とした全エネルギー %.4g%%, 帯域内 %.4g%%)。"
             % (sc["breaking"]["scale_applied"],
                sc["breaking"]["band_rms_before"],
                sc["breaking"]["band_energy_fraction"],
                sc["breaking"]["dropped_total_pct"],
                sc["breaking"]["dropped_band_pct"]))
    L.append("- rms 傾き: 合計 %.5g (A %.5g + C(t0) %.5g)。"
             % (numbers["realised"]["rms_slope"],
                numbers["realised"]["rms_slope_A"],
                numbers["realised"]["rms_slope_C_t0"]))
    L.append("- 変位場ヤードスティック (frame 96, green): mean %.4g px, "
             "p50 %.4g px, p99 %.4g px。"
             % (warp.get("g_mean_px"), warp.get("g_p50_px"),
                warp.get("g_p99_px")))
    L.append("- ループ継ぎ目 t=8 vs t=0 (正直な残差): %s。"
             % numbers["loop_seam_t8_vs_t0"]["max_abs_diff"])
    L.append("- 192 frames %dx%d, 24 fps, 8 s; still=frame-0096.png。"
             % (WIDTH, HEIGHT))
    L.append("")
    L.append("## 6. 失敗の記録")
    L.append("")
    fails = []
    if not g1.get("pass"):
        fails.append("G1 FAIL: 少なくとも1代表点が包絡外（後付け調整なし）。")
    if not g2.get("pass"):
        fails.append("G2 FAIL: 伝達チェックのいずれかが閾値外。")
    if not g3.get("pass"):
        fails.append("G3 FAIL: 無傷性が崩れた。")
    if not g4.get("pass"):
        fails.append("G4 FAIL: 出力 deliverable が不完全。")
    if not fails:
        fails.append("なし（全ゲート PASS）。ただし G1 の S0=0.001 と C 振幅"
                     " sigma_c=u*²/g は仮定であり測定ではない。")
    for f in fails:
        L.append("- " + f)
    L.append("")
    return "\n".join(L) + "\n"


def do_report():
    import shutil
    log("=" * 72)
    log("[report] numbers + REPORT + still")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write(write_report(numbers))
    shutil.copyfile(chunk_path(96), STILL96)
    cp.unit_done("report")
    log("  wrote %s" % os.path.relpath(NUMBERS, HERE))
    log("  wrote %s" % os.path.relpath(REPORT, HERE))
    return numbers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--g1", action="store_true")
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
        do_g1()
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
    if a.g1:
        do_g1()
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
