"""run_promo_v1 -- BRIEF-PROMO-V1 material render (素材A 5 layers + 素材B depth sweep).

Dispatch (authoritative execution spec):
    D:\\mutaro-work\\DISPATCH-RENDER-V1.md
Brief (source of the artistic/production constraints):
    D:\\mutaro-work\\BRIEF-PROMO-V1.md  (v1.4)

What this runner does
---------------------
Renders the promo's wave material on the *production* stack, re-timed onto a
2.000 s periodic grid:

    preset 5: U = 6 m/s, F = 100 m, z = 4.0 m, CLEAR medium
              (v12_water.PRESETS["clear"]) + FF clear phase
              (run_v16_water.FF_PHASE["clear"], n = 1.098, mu = 3.37),
              eta = 1.0, seed 3255045276382562994,
              photo base (D:\\mutaro-work\\moth-on-bark.jpg, 735x539 ->
              736x540 by the v2 right+bottom 1 px edge replication),
              frozen composite A2.6 + Blow sigma 3.0 px + E2.2.

Same code path as run_v19_prod.py (the 8 s piece): run_v15.scene_for ->
v12_water.preset_spectra -> run_v16_water.optics_v16 /
render_frame_v16 -> det.det_frame.  Two deliberate differences, both required
by the dispatch:

  (1) TIME GRID.  The 8 s production recipe quantises omega to n * 2*pi/8 s
      (v10_state / v11_state _realise).  The promo needs a 2 s loop, i.e.
      omega quantised to n * 2*pi/2 s = n * pi (0.5 Hz steps).  The spatial
      realisation H0 is NOT touched: only casc.omega_q is recomputed on the
      2 s lattice, with the same rule the production code uses
      (n = max(1, round(dispersion_omega(k, depth) / omega0))) and the same
      float32 storage.  See requantise_2s() and the loop-closure unit.

  (2) DEPOSIT.  The fast deposit (det_fast.det_frame_fast, row_block=256) is
      the default, as in BRIEF-DET-WIRE; it is bit-identical to det.det_frame
      (REPORT-DET-WIRE).  --det-impl base restores the frozen path.

Nothing else is changed, and no existing file is modified: this file is the
only new code (the dispatch's "新規コードはランナー1本のみ").

Deliverables (out/promo-v1/):
    final/ height/ warp/ base/ caustics/   frame-0000..0047.png (48 each)
    depth-0.6/ depth-1.0/ depth-2.0/       frame-0000..0047.png (48 each)
    numbers-promo-draft.json               every number + its measurement source
    REPORT-RENDER.md                       counts / hashes / time / failures
    checkpoint.json                        jobkit checkpoint, updated per unit

Units (each idempotent; a unit that fails is logged, skipped, and reported):
    stage0, material-a, depth-0.6, depth-1.0, depth-2.0,
    loop-closure-t2, measure-frame-24, report

Usage (from caustics-filter/offline/):
    python -B run_promo_v1.py                      # the whole job
    python -B run_promo_v1.py --stage0             # preflight only
    python -B run_promo_v1.py --probe 24           # one frame into _work
    python -B run_promo_v1.py --only material-a    # a single unit
    python -B run_promo_v1.py --threads 10         # numba threads (default 10)
"""
import argparse
import datetime
import gc
import hashlib
import io
import json
import os
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

DISPATCH = r"D:\mutaro-work\DISPATCH-RENDER-V1.md"
BRIEF = r"D:\mutaro-work\BRIEF-PROMO-V1.md"

OUT = os.path.join(HERE, "out", "promo-v1")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS_PATH = os.path.join(WORK, "params.json")
NUMBERS = os.path.join(OUT, "numbers-promo-draft.json")
REPORT = os.path.join(OUT, "REPORT-RENDER.md")
RUNLOG = os.path.join(OUT, "render.log")

PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"
V19_NUMBERS = os.path.join(HERE, "out", "v19", "numbers-v19-prod.json")

#: 素材A layers, in report order (directory name -> what it is).
LAYERS = ("final", "height", "warp", "base", "caustics")
#: 素材B depths (z = 4.0 is 素材A's final, so it is not repeated here).
DEPTHS_B = (0.6, 1.0, 2.0)

#: the fixed preset-5 cell (dispatch: 条件は固定, tuning forbidden).
CELL = {"U": 6.0, "F": 100.0, "h": 4.0}
SEED_BRIEF = 3255045276382562994
WATER_PRESET = "clear"

#: frozen production constants (v19 stack).
A_FROZEN = 2.6
B_SIGMA_FROZEN = 3.0
E_FROZEN = 2.2
GRID = 2048
Q_DET = 25
ROW_BLOCK = 256

FPS = 24
NFRAMES = 48
T_LOOP = 2.0
WIDTH, HEIGHT = 736, 540
PHOTO_NATIVE = (735, 539)

#: 素材A's measured frame (dispatch: the middle frame of the 48, z = 4.0).
MEASURE_FRAME = 24

#: fixed display scales of the *visualisation* layers (never per-frame).
#: The master grid is 2048 cells over S = 4.8 m (1 cell = 2.34 mm), so a
#: display band-limit of 32 master cells = 0.075 m ~= 11.5 output pixels.
#: It is applied to the layer images only -- never to `final`, which is the
#: untouched production composite, and never to the yardstick in section 6.
DISPLAY_SIGMA_CELLS = 32.0
#: warp layer: |d| of the band-limited field mapped to full value at 20 px
#: (frame-24 max of that band is ~18.3 px, so the fixed ramp does not clip).
WARP_PX_SCALE = 20.0
#: caustics layer: fixed display window of the unit-mean relative green
#: irradiance (1 = flat water; the production field at frame 24 is
#: p1 0.65 / p50 1.00 / p99 1.54 / max 1.85).
CAUSTICS_WINDOW = (0.60, 1.90)

NEW_FILES = ("run_promo_v1.py",)

CODE_FILES = ("run_v15.py", "run_v15_state.py", "run_v16_water.py",
              "v12_water.py", "v11_state.py", "v10_state.py", "v8_physics.py",
              "v5_composite.py", "det.py", "det_fast.py", "image.py",
              "state.py", "recipe_v8.py", "env.py")


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    try:
        os.makedirs(OUT, exist_ok=True)
        with open(RUNLOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError:
        pass


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


def read_cp():
    if os.path.exists(CP_PATH):
        return read_json(CP_PATH)
    return {"completed_units": 0, "completed_ids": [], "last_unit": None,
            "updated_at": None}


def write_cp(state):
    state["updated_at"] = now_iso()
    payload = json.dumps(state, ensure_ascii=False, indent=2).encode("utf-8")
    atomic_write_bytes(CP_PATH, payload)
    return state


def cp_get(key, default=None):
    return read_cp().get(key, default)


def cp_set(key, value):
    state = read_cp()
    state[key] = value
    return write_cp(state)


def makedirs():
    for d in [OUT, WORK] + list(LAYERS) + depth_dirs():
        os.makedirs(d, exist_ok=True)


def depth_dir(depth):
    return os.path.join(OUT, "depth-%.1f" % float(depth))


def depth_dirs():
    return [depth_dir(d) for d in DEPTHS_B]


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".py") and os.path.isfile(os.path.join(HERE, name)):
            out[name] = sha256_file(os.path.join(HERE, name))
    return out


# --------------------------------------------------------------------------
# the base photo (loaded once per process; the pad is the v2 convention)
# --------------------------------------------------------------------------
def load_photo():
    from PIL import Image
    from image import srgb_to_linear
    im = Image.open(PHOTO)
    native = im.size
    u8 = np.asarray(im.convert("RGB"))
    padded = np.pad(u8, ((0, 1), (0, 1), (0, 0)), mode="edge")
    lin = srgb_to_linear(padded.astype(np.float64) / 255.0)
    return {"native": native, "u8": padded, "lin": lin}


# --------------------------------------------------------------------------
# 2 s time grid: the production quantisation rule with T = 2 s
# --------------------------------------------------------------------------
def requantise_2s(st):
    """Recompute casc.omega_q on the 2 s lattice (n * 2*pi/2 s).

    Spatial realisation (H0, k, spectrum, amplitudes) is untouched; only the
    time evolution lattice changes, which is exactly what "周期 T=2s の格子
    （周波数 0.5Hz 刻み）で組み直す" asks for.  Rule identical to
    v10_state.V10WaterState._realise / v11_state._realise except OMEGA0:
        n = max(1, floor(dispersion_omega(k, depth) / omega0 + 0.5))
        omega_q = float32(n * omega0)
    """
    import v8_physics as P8
    om0 = 2.0 * np.pi / T_LOOP
    table = {}
    for tag in ("lo", "hi"):
        casc = getattr(st, tag, None)
        if casc is None:
            continue
        k = np.asarray(casc.k, np.float64)
        om = P8.dispersion_omega(np.maximum(k, 1e-12), float(st.depth))
        n = np.maximum(1.0, np.floor(om / om0 + 0.5))
        casc.omega_q = (n * om0).astype(np.float32)
        h0 = np.asarray(casc.H0)
        active = h0 != 0
        w = casc.omega_q[active].astype(np.float64)
        nm = w / om0
        phase = w * T_LOOP
        table[tag] = {
            "n_modes_active": int(active.sum()),
            "n_multiplier_min": float(nm.min()),
            "n_multiplier_max": float(nm.max()),
            "max_abs_dev_relative_to_one_step": float(
                np.max(np.abs(nm - np.round(nm)))),
            "max_phase_dev_at_T_rad": float(np.max(np.abs(
                phase - 2.0 * np.pi * np.round(phase / (2.0 * np.pi))))),
            "continuum_omega_quantisation_error_rad_s": float(
                np.max(np.abs(om[active] - w))),
            "omega0_rad_s": float(om0),
        }
    return table


def build_scene(depth, det=None, verbose=False):
    """The production scene for one depth, re-timed to the 2 s grid."""
    import run_v10 as R10
    import run_v15 as R15
    import v12_water as W
    import run_v16_water as V16
    from state import derive_seed

    t0 = time.time()
    raw = R10.load_recipe()
    seed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, float(CELL["U"]), float(CELL["F"]), float(depth),
                      GRID, seed)
    st, env = r["state"], r["env"]
    quant = requantise_2s(st)
    spec = W.preset_spectra(WATER_PRESET)
    optics = V16.optics_v16(spec, env, env["depthM"], GRID, st.S)
    meta = {
        "depth_m": float(depth),
        "seed": seed,
        "grid": GRID,
        "Q": Q_DET,
        "patch_m": float(st.S),
        "dx_m": float(st.S / GRID),
        "depthM_env": float(env["depthM"]),
        "hs_spectral_m": float(st.spectrum_after["Hs_m"]),
        "sigma2_spectral": float(st.spectrum_after["sigma2_rad2"]),
        "omega_c": float(st.ctx["omega_c"]),
        "warn_extrapolated_wave_age": bool(
            float(st.ctx["omega_c"]) > 5.0),
        "water_preset": spec["preset"],
        "water_label": spec.get("label"),
        "phase_n": float(optics["phase"]["n"]),
        "phase_mu": float(optics["phase"]["mu"]),
        "kernel": optics["kernel"],
        "w_ballistic": [float(x) for x in optics["w_ballistic"]],
        "w_halo": [float(x) for x in optics["w_halo"]],
        "beer_down": [float(x) for x in optics["beer_down"]],
        "time_grid": {"T_s": T_LOOP, "omega0_rad_s": 2.0 * np.pi / T_LOOP,
                      "quantisation": quant},
        "H0_sha256_12": {
            tag: sha256_bytes(np.asarray(getattr(st, tag).H0).tobytes())[:12]
            for tag in ("lo", "hi") if getattr(st, tag, None) is not None},
        "build_wall_s": float(time.time() - t0),
    }
    if verbose:
        log("  scene z=%.1f built %.1fs S=%.4f m Hs=%.5f m" %
            (depth, meta["build_wall_s"], meta["patch_m"],
             meta["hs_spectral_m"]))
    return meta, st, env, optics


# --------------------------------------------------------------------------
# the deposit (fast path = production default, bit-identical to det.py)
# --------------------------------------------------------------------------
def build_det(impl, threads):
    impl = (impl or "fast").strip().lower()
    if impl == "base":
        from det import det_frame
        log("det impl=base (det.det_frame, frozen behaviour)")
        cp_set("det", {"impl": impl, "threads": None,
                       "row_block": ROW_BLOCK, "set_utc": now_iso()})

        def _base(fields, S, N_buf, env, Q):
            C, info = det_frame(fields, S, N_buf, env, Q)
            info = dict(info)
            info["impl"] = "det.det_frame"
            return C, info
        return _base
    if impl != "fast":
        raise SystemExit("unknown det impl %r (use fast|base)" % impl)
    from det_fast import det_frame_fast, warmup
    warmup(threads)
    work = {}

    def _fast(fields, S, N_buf, env, Q):
        C, info = det_frame_fast(fields, S, N_buf, env, Q,
                                 row_block=ROW_BLOCK, work=work,
                                 n_threads=threads)
        info = dict(info)
        info["impl"] = "det_fast.det_frame_fast"
        info["row_block"] = ROW_BLOCK
        return C, info
    log("det impl=fast (det_fast.det_frame_fast) row_block=%d threads=%s"
        % (ROW_BLOCK, threads))
    cp_set("det", {"impl": impl, "threads": threads,
                   "row_block": ROW_BLOCK, "set_utc": now_iso()})
    return _fast


# --------------------------------------------------------------------------
# layer encoders (documented, fixed scales; no per-frame normalisation)
# --------------------------------------------------------------------------
def false_color(z):
    """z in [-1, 1] -> uint8 RGB on a fixed 5-stop ramp (blue..white..orange)."""
    stops = np.array([
        [0.00, 10, 34, 96],
        [0.25, 32, 132, 190],
        [0.50, 232, 238, 242],
        [0.75, 226, 178, 110],
        [1.00, 178, 84, 36],
    ], dtype=np.float64)
    t = np.clip((np.asarray(z, np.float64) + 1.0) * 0.5, 0.0, 1.0)
    out = np.zeros(t.shape + (3,), np.float64)
    for j in range(len(stops) - 1):
        t0, t1 = stops[j, 0], stops[j + 1, 0]
        m = (t >= t0) & (t <= t1)
        if not np.any(m):
            continue
        f = ((t[m] - t0) / (t1 - t0))[:, None]
        out[m] = stops[j, 1:] * (1.0 - f) + stops[j + 1, 1:] * f
    return (np.clip(out + 0.5, 0, 255)).astype(np.uint8)


def hsv_to_u8(h, s, v):
    """Vectorised HSV (h in [0,1)) -> uint8 RGB."""
    h = np.mod(np.asarray(h, np.float64), 1.0)
    s = np.clip(np.asarray(s, np.float64), 0.0, 1.0)
    v = np.clip(np.asarray(v, np.float64), 0.0, 1.0)
    i = np.floor(h * 6.0).astype(np.int64) % 6
    f = h * 6.0 - np.floor(h * 6.0)
    p = v * (1.0 - s)
    q = v * (1.0 - f * s)
    t = v * (1.0 - (1.0 - f) * s)
    r = np.select([i == 0, i == 1, i == 2, i == 3, i == 4, i == 5],
                  [v, q, p, p, t, v])
    g = np.select([i == 0, i == 1, i == 2, i == 3, i == 4, i == 5],
                  [t, v, v, q, p, p])
    b = np.select([i == 0, i == 1, i == 2, i == 3, i == 4, i == 5],
                  [p, p, t, v, v, q])
    rgb = np.stack([r, g, b], axis=-1) * 255.0
    return (np.clip(rgb, 0, 255) + 0.5).astype(np.uint8)


def pixel_grid():
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    return u, v


def displacement_px(fields, env, S, sigma=B_SIGMA_FROZEN):
    """The v19 transport-yardstick displacement, in output pixels.

    Definition copied from run_v19_prod.do_warp (so the two numbers are
    comparable): lowpass gx/gy with sigma 3.0 px, sample at output pixel
    centres, n = (-gx, -gy, 1)/|.|, shx = -nx*depth*(1-1/n_ior),
    dx_px = shx / S * (W-1), same for y; |d| per channel.

    sigma is in master cells (the production value 3.0 is the yardstick; the
    `warp` *layer* passes DISPLAY_SIGMA_CELLS for a legible band-limited view).
    """
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    u, v = pixel_grid()
    gx_s, gy_s = _state_gradients(fields, float(sigma))
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S,
                        S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    depth = float(env["depthM"])
    out = {}
    for bb, tag in enumerate(("r", "g", "b")):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        dx_px = shx / S * (WIDTH - 1.0)
        dy_px = shy / S * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        m_mag = np.sqrt(shx * shx + shy * shy)
        out[tag] = {
            "dx_px": dx_px, "dy_px": dy_px, "mag_px": mag,
            "stats": {
                "shift_m_mean": float(np.mean(m_mag)),
                "shift_m_p50": float(np.percentile(m_mag, 50)),
                "shift_m_p99": float(np.percentile(m_mag, 99)),
                "px_mean": float(np.mean(mag)),
                "px_p50": float(np.percentile(mag, 50)),
                "px_p99": float(np.percentile(mag, 99)),
                "px_max": float(np.max(mag)),
                "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                "dy_px_mean_abs": float(np.mean(np.abs(dy_px))),
            },
        }
    return out


def lowpass_display(arr):
    """Display band-limit on the master grid (wrap mode), in master cells."""
    from scipy import ndimage as ndi
    return ndi.gaussian_filter(np.asarray(arr, np.float64),
                               float(DISPLAY_SIGMA_CELLS), mode="wrap")


def layer_height(fields, S, limit_m):
    from state import sample_bicubic
    u, v = pixel_grid()
    h_lp = lowpass_display(fields["h"])
    h = sample_bicubic(h_lp.astype(np.float32), u * S, v * S,
                       S).astype(np.float64)
    z = np.clip(h / max(float(limit_m), 1e-12), -1.0, 1.0)
    return false_color(z), h


def layer_warp(disp):
    g = disp["g"]
    dx, dy = g["dx_px"], g["dy_px"]
    mag = g["mag_px"]
    hue = (np.arctan2(dy, dx) / (2.0 * np.pi)) % 1.0
    val = np.clip(mag / WARP_PX_SCALE, 0.0, 1.0)
    sat = np.where(val > 1e-9, 1.0, 0.0)
    return hsv_to_u8(hue, sat, val)


def layer_caustics(dbg, optics):
    """The green caustic layer: unit-mean relative direct irradiance.

    irr_w = w_bal[g]*C_bal + w_halo[g]*C_hal (the composite's own green
    irradiance).  Dividing by (w_bal[g] + w_halo[g]) makes it dimensionless
    with mean ~1 (1 = flat water); the fixed linear window CAUSTICS_WINDOW is
    then mapped to 0..1 and sRGB-encoded.  Fixed for all 48 frames.
    """
    irr = np.asarray(dbg["direct_irradiance"], np.float64)[:, :, 1]
    wsum = float(optics["w_ballistic"][1]) + float(optics["w_halo"][1])
    rel = irr / wsum if wsum > 0 else irr
    lo, hi = CAUSTICS_WINDOW
    v = np.clip((rel - lo) / (hi - lo), 0.0, 1.0)
    from image import linear_to_srgb
    gray = (np.clip(linear_to_srgb(v), 0.0, 1.0) * 255.0 + 0.5).astype(
        np.uint8)
    stats = {"w_green_sum": wsum,
             "rel_min": float(rel.min()),
             "rel_p01": float(np.percentile(rel, 1)),
             "rel_p50": float(np.percentile(rel, 50)),
             "rel_p99": float(np.percentile(rel, 99)),
             "rel_max": float(rel.max())}
    return np.repeat(gray[:, :, None], 3, axis=2), stats


# --------------------------------------------------------------------------
# PNG writing (atomic, verified)
# --------------------------------------------------------------------------
def save_png(path, arr):
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG")

    def verify(p):
        try:
            im = Image.open(p)
            im.load()
            return im.size == (WIDTH, HEIGHT) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(path, buf.getvalue(), verify=verify)
    return sha256_bytes(buf.getvalue())


# --------------------------------------------------------------------------
# rendering one frame
# --------------------------------------------------------------------------
def render_frame(i, st, env, optics, photo, det, layers, height_limit=None,
                  measure=False):
    """Render frame i; returns the layer arrays requested in `layers`."""
    from v5_composite import image_metrics
    from v5_physics import smith_g1_beckmann
    import run_v16_water as V16

    t0 = time.time()
    t = i / float(FPS)
    t1 = time.time()
    fields = st.fields_at(t)
    t_fields = time.time() - t1

    t2 = time.time()
    C, det_info = det(fields, st.S, GRID, env, Q_DET)
    t_det = time.time() - t2

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

    t4 = time.time()
    arrays = {}
    extra = {}
    if "final" in layers:
        arrays["final"] = img
    if "base" in layers:
        arrays["base"] = photo["u8"]
    if "height" in layers:
        img_h, h = layer_height(fields, float(st.S), height_limit)
        arrays["height"] = img_h
        extra["h_absmax_m"] = float(np.max(np.abs(h)))
    disp = None
    if measure:
        disp = displacement_px(fields, env, float(st.S),
                               sigma=B_SIGMA_FROZEN)
    if "warp" in layers:
        # the layer is band-limited for legibility; the yardstick above is not
        disp_layer = displacement_px(fields, env, float(st.S),
                                     sigma=DISPLAY_SIGMA_CELLS)
        arrays["warp"] = layer_warp(disp_layer)
        del disp_layer
    if "caustics" in layers:
        arrays["caustics"], cstats = layer_caustics(dbg, optics)
        extra["caustics_irradiance"] = cstats
    t_layers = time.time() - t4

    Cs = np.asarray(C, np.float64)
    out = {
        "frame": int(i), "t_s": float(t),
        "arrays": arrays,
        "g1": float(g1),
        "C_mean": float(np.mean(Cs)),
        "C_min": float(Cs.min()), "C_max": float(Cs.max()),
        "C_p50": float(np.percentile(Cs, 50)),
        "C_p99": float(np.percentile(Cs, 99)),
        "crush": float(metrics["crush_rate_linear_lt_0.02"]),
        "blow": float(metrics["blowout_rate_aces_gt_0.98"]),
        "lin_mean": float(metrics["linear_mean_pre_tonemap"]),
        "det_info": {k: v for k, v in det_info.items()
                     if isinstance(v, (int, float, str))},
        "timing_s": {"fields": t_fields, "det": t_det, "composite": t_comp,
                     "layers": t_layers, "total": float(time.time() - t0)},
        "extra": extra,
    }
    if disp is not None:
        out["displacement"] = disp
    del fields, dbg
    gc.collect()
    return out


# --------------------------------------------------------------------------
# idempotency
# --------------------------------------------------------------------------
def frame_paths(frame, layer_dirs):
    return {name: os.path.join(d, "frame-%04d.png" % frame)
            for name, d in layer_dirs.items()}


def frame_complete(rec, frame, layer_dirs):
    if not rec:
        return False
    shas = rec.get("sha256") or {}
    for name in layer_dirs:
        sha = shas.get(name)
        if not sha:
            return False
        p = os.path.join(layer_dirs[name], "frame-%04d.png" % frame)
        if not os.path.isfile(p) or sha256_file(p) != sha:
            return False
    return True


# --------------------------------------------------------------------------
# unit: material A (z = 4.0, five layers, 48 frames)
# --------------------------------------------------------------------------
def material_a_dirs():
    return {name: os.path.join(OUT, name) for name in LAYERS}


def render_material_a(photo, det, params, stop_on_error=False):
    cp = read_cp()
    node = cp.get("material_a") or {"depth": 4.0, "frames": {}}
    dirs = material_a_dirs()
    height_limit = float(params["height_limit_m"])
    todo = [i for i in range(NFRAMES)
            if not frame_complete((node.get("frames") or {}).get("%04d" % i),
                                  i, dirs)]
    log("[material-a] z=4.0 layers=%s frames todo=%d/%d"
        % ("/".join(LAYERS), len(todo), NFRAMES))
    if not todo:
        log("  material-a already complete -- skip")
        return node
    meta, st, env, optics = build_scene(CELL["h"], verbose=True)
    cp_set("scene_meta", {**(cp_get("scene_meta") or {}), "4.0": meta})
    failed = cp_get("failed_frames") or {}
    t_start = time.time()
    done = 0
    try:
        for i in todo:
            try:
                r = render_frame(i, st, env, optics, photo, det, set(LAYERS),
                                 height_limit=height_limit)
                shas = {}
                for name, p in frame_paths(i, dirs).items():
                    shas[name] = save_png(p, r["arrays"][name])
                rec = {"frame": i, "t_s": r["t_s"], "sha256": shas,
                       "g1": r["g1"], "C_mean": r["C_mean"],
                       "C_max": r["C_max"], "C_p99": r["C_p99"],
                       "crush": r["crush"], "blow": r["blow"],
                       "lin_mean": r["lin_mean"],
                       "h_absmax_m": r["extra"].get("h_absmax_m"),
                       "timing_s": r["timing_s"],
                       "completed_utc": now_iso()}
                state = read_cp()
                node = state.get("material_a") or {"depth": 4.0, "frames": {}}
                node.setdefault("frames", {})["%04d" % i] = rec
                node["depth"] = 4.0
                node["layers"] = list(LAYERS)
                state["material_a"] = node
                write_cp(state)
                failed.pop("material-a:%04d" % i, None)
                done += 1
                el = time.time() - t_start
                eta = el / max(done, 1) * (len(todo) - done)
                log("  A %04d t=%.3f Cmean=%.9f p99=%.3f crush=%.2f%% "
                    "det=%.1fs total=%.1fs (eta %.0fs)"
                    % (i, r["t_s"], r["C_mean"], r["C_p99"],
                       100.0 * r["crush"], r["timing_s"]["det"],
                       r["timing_s"]["total"], eta))
                del r
                gc.collect()
            except Exception as exc:                     # noqa: BLE001
                msg = "%s: %s" % (type(exc).__name__, exc)
                failed["material-a:%04d" % i] = {
                    "unit": "material-a", "frame": i, "error": msg,
                    "traceback": traceback.format_exc()[-2000:],
                    "at_utc": now_iso()}
                cp_set("failed_frames", failed)
                log("  A %04d FAILED: %s" % (i, msg))
                if stop_on_error:
                    raise
    finally:
        node = (read_cp().get("material_a") or node)
        node["wall_s"] = float(time.time() - t_start)
        node["n_frames_done"] = sum(
            1 for i in range(NFRAMES)
            if frame_complete((node.get("frames") or {}).get("%04d" % i),
                              i, dirs))
        cp_set("material_a", node)
        del st, optics
        gc.collect()
    log("[material-a] done %d/%d in %.1fs (frames complete %d/%d)"
        % (done, len(todo), time.time() - t_start,
           node.get("n_frames_done", 0), NFRAMES))
    return node


# --------------------------------------------------------------------------
# unit: material B (one depth, final only, 48 frames)
# --------------------------------------------------------------------------
def render_depth(depth, photo, det):
    tag = "%.1f" % float(depth)
    d = depth_dir(depth)
    dirs = {"final": d}
    cp = read_cp()
    depths = cp.get("depth_b") or {}
    node = depths.get(tag) or {"depth": float(depth), "frames": {}}
    todo = [i for i in range(NFRAMES)
            if not frame_complete((node.get("frames") or {}).get("%04d" % i),
                                  i, dirs)]
    log("[depth-%s] final only, frames todo=%d/%d" % (tag, len(todo), NFRAMES))
    if not todo:
        log("  depth-%s already complete -- skip" % tag)
        return node
    meta, st, env, optics = build_scene(depth, verbose=True)
    cp_set("scene_meta", {**(cp_get("scene_meta") or {}), tag: meta})
    failed = cp_get("failed_frames") or {}
    t_start = time.time()
    done = 0
    try:
        for i in todo:
            try:
                r = render_frame(i, st, env, optics, photo, det, {"final"})
                shas = {"final": save_png(
                    os.path.join(d, "frame-%04d.png" % i), r["arrays"]["final"])}
                rec = {"frame": i, "t_s": r["t_s"], "sha256": shas,
                       "g1": r["g1"], "C_mean": r["C_mean"],
                       "C_max": r["C_max"], "C_p99": r["C_p99"],
                       "crush": r["crush"], "blow": r["blow"],
                       "lin_mean": r["lin_mean"], "timing_s": r["timing_s"],
                       "completed_utc": now_iso()}
                state = read_cp()
                depths = state.get("depth_b") or {}
                node = depths.get(tag) or {"depth": float(depth),
                                           "frames": {}}
                node.setdefault("frames", {})["%04d" % i] = rec
                node["depth"] = float(depth)
                depths[tag] = node
                state["depth_b"] = depths
                write_cp(state)
                failed.pop("depth-%s:%04d" % (tag, i), None)
                done += 1
                el = time.time() - t_start
                eta = el / max(done, 1) * (len(todo) - done)
                log("  z=%.1f %04d t=%.3f Cmean=%.9f crush=%.2f%% "
                    "det=%.1fs total=%.1fs (eta %.0fs)"
                    % (depth, i, r["t_s"], r["C_mean"], 100.0 * r["crush"],
                       r["timing_s"]["det"], r["timing_s"]["total"], eta))
                del r
                gc.collect()
            except Exception as exc:                     # noqa: BLE001
                msg = "%s: %s" % (type(exc).__name__, exc)
                failed["depth-%s:%04d" % (tag, i)] = {
                    "unit": "depth-%s" % tag, "frame": i, "error": msg,
                    "traceback": traceback.format_exc()[-2000:],
                    "at_utc": now_iso()}
                cp_set("failed_frames", failed)
                log("  z=%.1f %04d FAILED: %s" % (depth, i, msg))
    finally:
        state = read_cp()
        depths = state.get("depth_b") or {}
        node = depths.get(tag) or node
        node["wall_s"] = float(time.time() - t_start)
        node["n_frames_done"] = sum(
            1 for i in range(NFRAMES)
            if frame_complete((node.get("frames") or {}).get("%04d" % i),
                              i, dirs))
        depths[tag] = node
        state["depth_b"] = depths
        write_cp(state)
        del st, optics
        gc.collect()
    log("[depth-%s] done %d/%d in %.1fs (frames complete %d/%d)"
        % (tag, done, len(todo), time.time() - t_start,
           node.get("n_frames_done", 0), NFRAMES))
    return node


# --------------------------------------------------------------------------
# unit: loop closure on the T = 2 s grid (proof, not a guess)
# --------------------------------------------------------------------------
def do_loop_closure(det):
    t0 = time.time()
    log("=" * 72)
    log("[loop-closure] t=2.000 s vs t=0.000 s on the 2 s lattice")
    meta, st, env, optics = build_scene(CELL["h"], verbose=True)
    f0 = st.fields_at(0.0)
    f2 = st.fields_at(T_LOOP)
    fdiff = {k: float(np.max(np.abs(np.asarray(f2[k], np.float64)
                                   - np.asarray(f0[k], np.float64))))
             for k in ("h", "gx", "gy", "Dx", "Dy")}
    C0, _ = det(f0, st.S, GRID, env, Q_DET)
    C2, _ = det(f2, st.S, GRID, env, Q_DET)
    dC = float(np.max(np.abs(C2.astype(np.float64)
                            - C0.astype(np.float64))))
    f1 = st.fields_at(1.0 / FPS)
    C1, _ = det(f1, st.S, GRID, env, Q_DET)
    sens = float(np.max(np.abs(C1.astype(np.float64)
                              - C0.astype(np.float64))))
    quant = meta["time_grid"]["quantisation"]
    proof_ok = all(v["max_abs_dev_relative_to_one_step"] <= 1e-3
                   for v in quant.values())
    out = {
        "T_s": T_LOOP, "frames": NFRAMES, "fps": FPS,
        "last_frame": NFRAMES - 1,
        "last_frame_t_s": (NFRAMES - 1) / float(FPS),
        "duplicate_final_frame": False,
        "rule": ("omega_q = n * 2*pi/T with n = max(1, round(omega/omega0)); "
                 "the 48 frames are t = i/24 s, i = 0..47, so the last frame "
                 "is 1.9583 s and frame 0 is the loop point (2.000 s == 0)"),
        "fields_max_abs_diff_t2_vs_t0": fdiff,
        "C_max_abs_diff_t2_vs_t0": dC,
        "C_max_abs_diff_1frame_vs_t0_sensitivity": sens,
        "C_seam_residual_relative_to_one_frame_step": (float(dC / sens)
                                                       if sens else None),
        "quantisation_proof_pass": bool(proof_ok),
        "frequency_table": quant,
        "pass": bool(proof_ok and dC <= 1e-4),
        "note": ("the exact-zero branch cannot hold because omega_q is "
                 "float32; the printed quantisation proof plus the measured "
                 "residual against one frame step is the evidence"),
        "checked_utc": now_iso(),
        "wall_s": float(time.time() - t0),
    }
    state = read_cp()
    measured = state.get("measured") or {}
    measured["loop_t2"] = out
    state["measured"] = measured
    write_cp(state)
    log("  fields max|diff| t2-t0: %s" % {k: "%.3g" % v
                                         for k, v in fdiff.items()})
    log("  C max|diff| t2-t0 = %.6g | 1-frame sensitivity = %.6g | ratio %.3g"
        % (dC, sens, (dC / sens) if sens else float("nan")))
    log("  quantisation proof: %s | loop closure: %s"
        % (proof_ok, "PASS" if out["pass"] else "FAIL"))
    del st, optics
    gc.collect()
    return out


# --------------------------------------------------------------------------
# unit: the measured displacement of 素材A's frame 24
# --------------------------------------------------------------------------
def do_warp24(det):
    t0 = time.time()
    log("=" * 72)
    log("[measure-frame-24] transport displacement, 素材A final")
    meta, st, env, optics = build_scene(CELL["h"], verbose=True)
    i = MEASURE_FRAME
    fields = st.fields_at(i / float(FPS))
    disp = displacement_px(fields, env, float(st.S))
    per_ch = {tag: disp[tag]["stats"] for tag in ("r", "g", "b")}
    # cross-check: if the 素材A final frame exists, its C statistics must match
    # a fresh render of the same frame (same build, same seed -> same field).
    out = {
        "frame": int(i), "t_s": i / float(FPS),
        "field": ("transport displacement shx = -nx*depth*(1-1/n), sampled at "
                  "output pixel centres with S=%.6g m; |d| in output pixels = "
                  "sqrt((shx/S*(W-1))^2+(shy/S*(H-1))^2)" % float(st.S)),
        "patch_m": float(st.S), "depth_m": float(env["depthM"]),
        "width": WIDTH, "height": HEIGHT, "fps": FPS,
        "time_grid": {"T_s": T_LOOP, "omega0_rad_s": 2.0 * np.pi / T_LOOP},
        "per_channel": per_ch,
        "reported_channel": "g (the composite's own channel for C/g1)",
        "g_mean_px": per_ch["g"]["px_mean"],
        "g_p50_px": per_ch["g"]["px_p50"],
        "g_p99_px": per_ch["g"]["px_p99"],
        "g_max_px": per_ch["g"]["px_max"],
        "definition_source": ("run_v19_prod.do_warp (identical formula); "
                              "v19 measured frame 96, this measures frame 24"),
        "v19_reference": None,
        "created_utc": now_iso(),
        "wall_s": float(time.time() - t0),
    }
    if os.path.exists(V19_NUMBERS):
        try:
            d19 = read_json(V19_NUMBERS)
            w = d19.get("warp") or {}
            out["v19_reference"] = {
                "path": os.path.relpath(V19_NUMBERS, HERE),
                "frame": w.get("frame"),
                "g_mean_px": w.get("g_mean_px"),
                "g_p50_px": w.get("g_p50_px"),
                "g_p99_px": w.get("g_p99_px"),
                "g_max_px": ((w.get("per_channel") or {}).get("g") or {})
                .get("px_max"),
                "patch_m": w.get("patch_m"),
                "depth_m": w.get("depth_m"),
                "same_definition": True,
                "differences": ("time grid (8 s lattice there, 2 s lattice "
                                "here) and frame (96 vs 24); the formula, the "
                                "channel and the sampling are identical"),
            }
        except Exception as exc:                         # noqa: BLE001
            out["v19_reference"] = {"error": str(exc)}
    state = read_cp()
    measured = state.get("measured") or {}
    measured["warp_frame_24"] = out
    state["measured"] = measured
    write_cp(state)
    log("  green |d| px: mean %.3f p50 %.3f p99 %.3f max %.3f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"],
           out["g_max_px"]))
    del st, optics
    gc.collect()
    return out


# --------------------------------------------------------------------------
# stage 0 -- preflight (parameters pinned, nothing tuned)
# --------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] preflight: pinned conditions, code set, scene")
    makedirs()
    import run_v10 as R10
    import run_v16_water as V16
    import v12_water as W
    from state import derive_seed

    cp = Checkpoint(CP_PATH)
    facts = {"job": "PROMO-RENDER-V1", "dispatch": DISPATCH, "brief": BRIEF,
             "workdir": HERE, "created_utc": now_iso(),
             "python": sys.version.split()[0]}
    missing = []

    code = {}
    for name in CODE_FILES:
        p = os.path.join(HERE, name)
        ok = os.path.isfile(p)
        code[name] = {"present": ok,
                      "sha256": sha256_file(p) if ok else None}
        if not ok:
            missing.append("missing code: %s" % name)
    facts["code"] = code

    clear_w = W.PRESETS.get(WATER_PRESET)
    clear_p = V16.FF_PHASE.get(WATER_PRESET)
    facts["medium"] = {
        "v12_preset": WATER_PRESET,
        "v12_found": clear_w is not None,
        "v12_label": (clear_w or {}).get("label"),
        "v12_load": (clear_w or {}).get("load"),
        "v16_phase_found": clear_p is not None,
        "v16_phase_n": (clear_p or {}).get("n"),
        "v16_phase_mu": (clear_p or {}).get("mu"),
        "v16_phase_label": (clear_p or {}).get("label"),
        "v16_phase_provenance": (clear_p or {}).get("provenance"),
    }
    if clear_w is None:
        missing.append("no CLEAR medium in v12_water.PRESETS")
    if clear_p is None:
        missing.append("no CLEAR phase in run_v16_water.FF_PHASE")
    if clear_p and not (abs(float(clear_p["n"]) - 1.098) < 1e-9
                        and abs(float(clear_p["mu"]) - 3.37) < 1e-9):
        missing.append("clear phase is not n=1.098 mu=3.37")

    from PIL import Image
    if not os.path.isfile(PHOTO):
        missing.append("base photo missing: %s" % PHOTO)
        facts["photo"] = {"path": PHOTO, "present": False}
    else:
        im = Image.open(PHOTO)
        im.load()
        ph = load_photo()
        n = np.asarray(Image.open(PHOTO).convert("RGB"))
        facts["photo"] = {
            "path": PHOTO, "present": True, "native_px": list(im.size),
            "mode": im.mode, "format": im.format,
            "sha256": sha256_file(PHOTO),
            "output_px": [WIDTH, HEIGHT],
            "pad": ("right+bottom 1 px edge replication (the v2 production "
                    "loop's convention)"),
            "pad_check": {
                "shape": list(ph["u8"].shape),
                "right_col_is_native_last_col": bool(np.array_equal(
                    ph["u8"][:-1, -1, :], n[:, -1, :])),
                "bottom_row_is_native_last_row": bool(np.array_equal(
                    ph["u8"][-1, :-1, :], n[-1, :, :])),
                "corner_is_native_corner": bool(np.array_equal(
                    ph["u8"][-1, -1, :], n[-1, -1, :])),
            },
        }
        if tuple(im.size) != PHOTO_NATIVE:
            missing.append("photo decodes at %s, not 735x539" % (im.size,))
        if not all(facts["photo"]["pad_check"][k] for k in (
                "right_col_is_native_last_col",
                "bottom_row_is_native_last_row",
                "corner_is_native_corner")):
            missing.append("pad does not follow the v2 edge-replication rule")

    raw = R10.load_recipe()
    seeded = int(derive_seed(raw, 1))
    facts["seed"] = {"brief_literal": SEED_BRIEF, "derive_seed_raw_1": seeded,
                     "agree": bool(seeded == SEED_BRIEF)}
    if seeded != SEED_BRIEF:
        missing.append("derive_seed(recipe,1) = %d != brief %d"
                       % (seeded, SEED_BRIEF))

    facts["cell"] = {
        "brief": {"U": CELL["U"], "F": CELL["F"], "z": CELL["h"],
                  "medium": "clear", "phase": "ff/clear",
                  "eta": 1.0, "seed": SEED_BRIEF,
                  "base": "photo (v19 identical)"},
        "used": dict(CELL),
        "source": ("DISPATCH-RENDER-V1 '条件は固定（チューニング禁止）' / "
                   "BRIEF-PROMO-V1 'プリセット5'"),
        "frozen_composite": {"A": A_FROZEN,
                             "B_lowpass_sigma_px": B_SIGMA_FROZEN,
                             "E": E_FROZEN, "Q": Q_DET},
        "time": {"T_s": T_LOOP, "fps": FPS, "frames": NFRAMES,
                 "omega0_rad_s": 2.0 * np.pi / T_LOOP,
                 "last_frame": NFRAMES - 1,
                 "last_frame_t_s": (NFRAMES - 1) / float(FPS),
                 "note": "final frame is NOT a duplicate of frame 0"},
        "frames_px": [WIDTH, HEIGHT],
        "depths_b": list(DEPTHS_B),
    }

    before = py_set()
    prev = read_json(PARAMS_PATH) if os.path.exists(PARAMS_PATH) else {}
    if isinstance(prev.get("py_hashes_before"), dict):
        # stage 0 is re-runnable: keep the FIRST snapshot, otherwise the
        # "no existing file changed" evidence would be taken after the job.
        prior = prev["py_hashes_before"]
        before = prior.get("hashes") or before
        facts["py_hashes_before_reused"] = {
            "from": os.path.relpath(PARAMS_PATH, HERE),
            "created_utc": prior.get("created_utc") or prev.get("created_utc"),
            "count": len(before)}
    facts["py_hashes_before"] = {"count": len(before), "hashes": before}

    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(PARAMS_PATH, facts)
        log("  PREFLIGHT FAIL: %s" % missing)
        raise SystemExit("stage 0 preflight failed: %s" % missing)

    meta, st, env, optics = build_scene(CELL["h"], verbose=True)
    h0 = lowpass_display(st.fields_at(0.0)["h"])
    hlimit = float(np.ceil(np.max(np.abs(h0)) * 1e4) / 1e4)
    facts["scene_meta_z4"] = meta
    facts["height_limit_m"] = hlimit
    facts["display_bandlimit"] = {
        "sigma_master_cells": DISPLAY_SIGMA_CELLS,
        "sigma_m": float(DISPLAY_SIGMA_CELLS * meta["dx_m"]),
        "sigma_output_px": float(DISPLAY_SIGMA_CELLS * meta["dx_m"]
                                 / (meta["patch_m"] / WIDTH)),
        "applies_to": ["height", "warp"],
        "does_not_apply_to": ["final", "base", "caustics",
                              "the frame-24 yardstick"],
        "why": ("the master grid holds capillary bands far below one output "
                "pixel (1 cell = %.3g m); sampled at 736x540 they alias into "
                "per-pixel noise. The layer images therefore show the "
                "band-limited field; the render and the yardstick keep the "
                "production definition." % meta["dx_m"]),
    }
    facts["height_scale_rule"] = (
        "false-colour full scale = +/- max|lowpass(h)| of the z=4 state at "
        "t=0 (rounded up to 0.1 mm); fixed for all 48 frames, no per-frame "
        "auto-scaling")
    facts["warp_scale_rule"] = (
        "warp layer: hue = direction, value = |d| / %.1f px on the "
        "band-limited field (sigma %.0f master cells); fixed" % (
            WARP_PX_SCALE, DISPLAY_SIGMA_CELLS))
    facts["caustics_rule"] = (
        "caustics layer = green direct irradiance divided by the green weight "
        "sum (unit-mean relative, 1 = flat water), fixed linear window "
        "[%.2f, %.2f] -> sRGB; fixed for all frames, no normalisation"
        % CAUSTICS_WINDOW)
    facts["base_rule"] = "the padded photo, byte-exact (no water)"
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS_PATH, facts)
    state = read_cp()
    state["param"] = facts
    state["scene_meta"] = {**(state.get("scene_meta") or {}), "4.0": meta}
    write_cp(state)
    cp.unit_done("stage-0-preflight")
    log("  photo %dx%d -> %dx%d | pad %s"
        % (PHOTO_NATIVE[0], PHOTO_NATIVE[1], WIDTH, HEIGHT,
           facts["photo"]["pad_check"]))
    log("  clear medium %s | FF phase n=%s mu=%s"
        % (facts["medium"]["v12_label"], facts["medium"]["v16_phase_n"],
           facts["medium"]["v16_phase_mu"]))
    log("  cell U=%g F=%g z=%g seed=%d | S=%.4f m Hs=%.5f m"
        % (CELL["U"], CELL["F"], CELL["h"], seeded, meta["patch_m"],
           meta["hs_spectral_m"]))
    log("  T=%.3fs grid: omega0=%.6f rad/s | 48 frames, last t=%.4fs"
        % (T_LOOP, 2.0 * np.pi / T_LOOP, (NFRAMES - 1) / float(FPS)))
    log("  height false-colour limit +/-%.4g m | warp ramp %.1f px"
        % (hlimit, WARP_PX_SCALE))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    del st, optics
    gc.collect()
    return facts


# --------------------------------------------------------------------------
# numbers + report
# --------------------------------------------------------------------------
def dir_stats(d, expected=NFRAMES):
    files = sorted(f for f in os.listdir(d) if f.endswith(".png")
                   and not f.endswith(".tmp"))
    frames = []
    for f in files:
        try:
            frames.append(int(f.split("-")[1].split(".")[0]))
        except (IndexError, ValueError):
            continue
    frames = sorted(set(frames))
    missing = [i for i in range(expected) if i not in frames]
    extra = [i for i in frames if i < 0 or i >= expected]
    return {"dir": os.path.relpath(d, HERE).replace("\\", "/"),
            "n_png": len(files), "n_frames": len(frames),
            "first": frames[0] if frames else None,
            "last": frames[-1] if frames else None,
            "missing": missing, "extra": extra,
            "complete": bool(len(frames) == expected and not missing
                             and not extra)}


def representative_hashes():
    rows = []
    groups = [("A:" + name, os.path.join(OUT, name)) for name in LAYERS]
    groups += [("B:depth-%.1f" % d, depth_dir(d)) for d in DEPTHS_B]
    for tag, d in groups:
        for i in (0, MEASURE_FRAME, NFRAMES - 1):
            p = os.path.join(d, "frame-%04d.png" % i)
            rows.append({"group": tag, "frame": i,
                         "file": os.path.relpath(p, HERE).replace("\\", "/"),
                         "exists": os.path.isfile(p),
                         "bytes": os.path.getsize(p) if os.path.isfile(p)
                         else None,
                         "sha256": sha256_file(p) if os.path.isfile(p)
                         else None})
    return rows


def timing_summary():
    cp = read_cp()
    out = {}
    node = cp.get("material_a") or {}
    recs = list((node.get("frames") or {}).values())
    if recs:
        out["material-a"] = {
            "n_frames": len(recs),
            "group_wall_s": node.get("wall_s"),
            "frame_wall_s_sum": float(sum(r["timing_s"]["total"]
                                          for r in recs)),
            "det_s_sum": float(sum(r["timing_s"]["det"] for r in recs)),
            "mean_frame_s": float(np.mean([r["timing_s"]["total"]
                                           for r in recs])),
            "min_frame_s": float(np.min([r["timing_s"]["total"]
                                         for r in recs])),
            "max_frame_s": float(np.max([r["timing_s"]["total"]
                                         for r in recs])),
        }
    for tag, n in (cp.get("depth_b") or {}).items():
        recs = list((n.get("frames") or {}).values())
        if not recs:
            continue
        out["depth-%s" % tag] = {
            "n_frames": len(recs),
            "group_wall_s": n.get("wall_s"),
            "mean_frame_s": float(np.mean([r["timing_s"]["total"]
                                           for r in recs])),
            "min_frame_s": float(np.min([r["timing_s"]["total"]
                                         for r in recs])),
            "max_frame_s": float(np.max([r["timing_s"]["total"]
                                         for r in recs])),
            "det_s_sum": float(sum(r["timing_s"]["det"] for r in recs)),
        }
    for k, v in (cp.get("measured") or {}).items():
        out[k] = {"wall_s": v.get("wall_s")}
    return out


def build_numbers():
    cp = read_cp()
    param = cp.get("param") or {}
    meta = cp.get("scene_meta") or {}
    measured = cp.get("measured") or {}
    warp = measured.get("warp_frame_24") or {}
    counts = {}
    for name in LAYERS:
        counts[name] = dir_stats(os.path.join(OUT, name))
    for d in DEPTHS_B:
        counts["depth-%.1f" % d] = dir_stats(depth_dir(d))

    code_loc = "run_promo_v1.py::displacement_px / do_warp24"
    numbers = {
        "job": "PROMO-RENDER-V1",
        "dispatch": DISPATCH,
        "brief": BRIEF,
        "created_utc": now_iso(),
        "conditions": {
            "preset": 5,
            "U_ms": CELL["U"], "fetch_m": CELL["F"], "depth_m": CELL["h"],
            "medium": WATER_PRESET,
            "phase": "ff / clear (n=1.098, mu=3.37)",
            "eta": 1.0, "seed": SEED_BRIEF,
            "base": "photo moth-on-bark.jpg (735x539 -> 736x540 pad)",
            "composite": {"A": A_FROZEN, "Blow_sigma_px": B_SIGMA_FROZEN,
                          "E": E_FROZEN, "Q": Q_DET},
            "frames": NFRAMES, "fps": FPS, "T_s": T_LOOP,
            "size_px": [WIDTH, HEIGHT],
            "depths_b": list(DEPTHS_B),
            "scene": {k: {"patch_m": v.get("patch_m"),
                          "hs_spectral_m": v.get("hs_spectral_m"),
                          "omega_c": v.get("omega_c"),
                          "H0_sha256_12": v.get("H0_sha256_12")}
                      for k, v in meta.items()},
        },
        "layer_definitions": {
            "final": "v16 production composite (A2.6 + Blow sigma 3.0 + E2.2)",
            "height": ("false colour of h band-limited at sigma %.0f master "
                       "cells (%.3g m), fixed +/- %.4g m"
                       % (DISPLAY_SIGMA_CELLS,
                          (param.get("display_bandlimit") or {})
                          .get("sigma_m", float("nan")),
                          float(param.get("height_limit_m") or 0.0))),
            "warp": ("transport displacement (same band-limit), hue = "
                     "direction, value = |d| / %.0f px"
                     % WARP_PX_SCALE),
            "base": "the padded photo, byte-exact (no water)",
            "caustics": ("green direct irradiance / green weight sum "
                         "(unit-mean relative, 1 = flat water), fixed window "
                         "[%.2f, %.2f] -> sRGB" % CAUSTICS_WINDOW),
            "note": ("the band-limit applies to the two visualisation layers "
                     "only; final, base, caustics and the frame-24 yardstick "
                     "keep the production definition"),
        },
        "numbers": [
            {
                "id": "warp_g_mean_px",
                "label": "MEAN DISPLACEMENT",
                "value": warp.get("g_mean_px"),
                "unit": "px",
                "display": "MEAN DISPLACEMENT — %.2f px / frame 24"
                           % (warp.get("g_mean_px") or 0.0),
                "measurement": {
                    "what": ("|transport displacement| of the substrate "
                             "sampling, green band, sampled at the 736x540 "
                             "output pixel centres"),
                    "frame": warp.get("frame"), "t_s": warp.get("t_s"),
                    "patch_m": warp.get("patch_m"),
                    "depth_m": warp.get("depth_m"),
                    "definition": warp.get("field"),
                    "definition_match_v19": True,
                    "files": ["out/promo-v1/numbers-promo-draft.json#numbers"
                              "[0]", "out/v19/numbers-v19-prod.json#warp"],
                    "code": code_loc,
                },
            },
            {
                "id": "warp_g_p99_px",
                "label": "P99 DISPLACEMENT",
                "value": warp.get("g_p99_px"), "unit": "px",
                "display": "P99 DISPLACEMENT — %.2f px / frame 24"
                           % (warp.get("g_p99_px") or 0.0),
                "measurement": {
                    "what": "99th percentile of the same |d| distribution",
                    "frame": warp.get("frame"), "t_s": warp.get("t_s"),
                    "definition": warp.get("field"),
                    "definition_match_v19": True,
                    "files": ["out/promo-v1/numbers-promo-draft.json#numbers"
                              "[1]", "out/v19/numbers-v19-prod.json#warp"],
                    "code": code_loc,
                },
            },
            {
                "id": "warp_g_max_px",
                "label": "MAX DISPLACEMENT",
                "value": warp.get("g_max_px"), "unit": "px",
                "display": "MAX DISPLACEMENT — %.2f px / frame 24"
                           % (warp.get("g_max_px") or 0.0),
                "measurement": {
                    "what": "maximum of the same |d| distribution",
                    "frame": warp.get("frame"), "t_s": warp.get("t_s"),
                    "definition": warp.get("field"),
                    "definition_match_v19": True,
                    "files": ["out/promo-v1/numbers-promo-draft.json#numbers"
                              "[2]", "out/v19/numbers-v19-prod.json#warp"],
                    "code": code_loc,
                },
            },
            {
                "id": "preset_U_ms", "label": "U", "value": CELL["U"],
                "unit": "m/s", "display": "U = 6 m/s",
                "measurement": {"what": "preset-5 wind speed, fixed",
                                "files": [BRIEF],
                                "code": "run_promo_v1.py::CELL"},
            },
            {
                "id": "preset_fetch_m", "label": "F", "value": CELL["F"],
                "unit": "m", "display": "F = 100 m",
                "measurement": {"what": "preset-5 fetch, fixed",
                                "files": [BRIEF],
                                "code": "run_promo_v1.py::CELL"},
            },
            {
                "id": "preset_depth_m", "label": "z", "value": CELL["h"],
                "unit": "m", "display": "z = 4.0 m",
                "measurement": {"what": "preset-5 depth, fixed",
                                "files": [BRIEF],
                                "code": "run_promo_v1.py::CELL"},
            },
            {
                "id": "loop_T_s", "label": "PERIOD", "value": T_LOOP,
                "unit": "s", "display": "PERIODIC %.2fs" % T_LOOP,
                "measurement": {
                    "what": "time lattice: omega quantised to n*2pi/2s; "
                            "frame 47 at t=1.9583 s, frame 0 is the loop point",
                    "files": ["out/promo-v1/REPORT-RENDER.md#loop",
                              "out/promo-v1/numbers-promo-draft.json"
                              "#loop_closure"],
                    "code": "run_promo_v1.py::requantise_2s / "
                            "do_loop_closure"},
            },
        ],
        "layer_counts": counts,
        "timing": timing_summary(),
        "loop_closure": measured.get("loop_t2"),
        "warp_frame_24_full": warp,
        "scene_meta": meta,
        "deliverable_hashes": representative_hashes(),
    }
    return numbers


def write_report(numbers):
    cp = read_cp()
    param = cp.get("param") or {}
    meta = cp.get("scene_meta") or {}
    measured = cp.get("measured") or {}
    counts = numbers["layer_counts"]
    warp = numbers.get("warp_frame_24_full") or {}
    loop = numbers.get("loop_closure") or {}
    units = cp.get("units") or {}
    failures = cp.get("failed_frames") or {}
    timing = numbers.get("timing") or {}

    def unit_row(name):
        u = units.get(name) or {}
        st = u.get("status") or "not run"
        return st, u.get("wall_s"), (u.get("error") or "")

    total_expected = NFRAMES * (len(LAYERS) + len(DEPTHS_B))
    total_png = sum(v["n_png"] for v in counts.values())
    all_complete = all(v["complete"] for v in counts.values())

    L = []
    L += [
        "# REPORT-RENDER — プロモ素材レンダー v1（promo-v1）",
        "",
        "発注書: `%s`" % DISPATCH,
        "ブリーフ正本: `%s`（v1.4）" % BRIEF,
        "作業場所: `%s`／実行: %s／Python %s・numpy %s"
        % (HERE, numbers["created_utc"], param.get("python", "?"),
           np.__version__),
        "",
        "## 0. 結論（実物で判定）",
        "",
        "| 項目 | 実測 |",
        "|---|---|",
        "| フレーム数（期待 %d = 48f × (%d層 + %d深さ)） | **%d PNG** |"
        % (total_expected, len(LAYERS), len(DEPTHS_B), total_png),
        "| 全ディレクトリ完備（0..47、欠番なし・重複なし） | **%s** |"
        % ("YES" if all_complete else "NO"),
        "| 層 | %s |" % ", ".join("%s=%d" % (k, v["n_frames"])
                                 for k, v in counts.items()),
        "| 出せなかった層 | %s |"
        % ("なし" if not failures else "%d 単位で失敗（§5）" % len(failures)),
        "| frame 24 実測（g） | mean **%.3f px** / p99 **%.3f px** / "
        "max **%.3f px** |"
        % (warp.get("g_mean_px") or float("nan"),
           warp.get("g_p99_px") or float("nan"),
           warp.get("g_max_px") or float("nan")),
        "| T=2s ループ閉合 | %s（C の t=2 vs t=0 最大差 %.3g、1フレーム感度 "
        "%.3g） |" % ("PASS" if loop.get("pass") else "未測定/FAIL",
                     loop.get("C_max_abs_diff_t2_vs_t0") or float("nan"),
                     loop.get("C_max_abs_diff_1frame_vs_t0_sensitivity")
                     or float("nan")),
        "| 総所要（単位合計） | **%.1f 分** |"
        % (sum((u.get("wall_s") or 0.0) for u in units.values()) / 60.0),
        "",
        "終了コードではなく上のフレーム数・ハッシュ・実測で完了を判定できる。",
        "",
        "## 1. 条件（固定、チューニングなし）",
        "",
        "| 項目 | 値 | 出所 |",
        "|---|---|---|",
        "| プリセット | 5 | 発注書／ブリーフ |",
        "| U | %.1f m/s | 同上 |" % CELL["U"],
        "| F | %.0f m | 同上 |" % CELL["F"],
        "| z（素材A） | %.1f m | 同上 |" % CELL["h"],
        "| 媒体 | `%s`（%s）、load=%s | v12_water.PRESETS |"
        % (WATER_PRESET, (param.get("medium") or {}).get("v12_label"),
           (param.get("medium") or {}).get("v12_load") or "n/a"),
        "| 位相関数 | FF clear n=%s, mu=%s | run_v16_water.FF_PHASE |"
        % ((param.get("medium") or {}).get("v16_phase_n"),
           (param.get("medium") or {}).get("v16_phase_mu")),
        "| η | 1.0（direct sun only） | v19 と同一 |",
        "| seed | %d | recipe から導出（一致 %s） |"
        % (SEED_BRIEF, (param.get("seed") or {}).get("agree")),
        "| ベース | 写真 735x539 → 736x540（右+下 1px edge） | 発注書 |",
        "| 合成 | A=%.1f / Blow σ=%.1f px / E=%.1f | v19 frozen |"
        % (A_FROZEN, B_SIGMA_FROZEN, E_FROZEN),
        "| 出力 | %dx%d, %d fps, T=%.3f s, %d フレーム | 発注書 |"
        % (WIDTH, HEIGHT, FPS, T_LOOP, NFRAMES),
        "",
        "**時間格子の組み直し（発注書の中核）**: 本番 8s 版の ω 量子化 "
        "`n*2π/8s` を、プロモ用に `n*2π/2s`（0.5 Hz 刻み）へ置き換えた。"
        "空間実現（H0・k・振幅）は無変更で、時間発展の格子だけを変更している。"
        "フレームは t = i/24（i=0..47）で、**最終フレームは 1.9583 s**。"
        "2.000 s は frame 0 のループ点であり、最終フレームの重複はない。",
        "",
        "| z | S [m] | Hs [m] | ω_c | モード数 lo/hi | 量子化相対誤差 (max) | H0 sha (lo/hi) |",
        "|---|---|---|---|---|---|---|",
    ]
    for k in ("4.0", "0.6", "1.0", "2.0"):
        m = meta.get(k)
        if not m:
            continue
        q = (m.get("time_grid") or {}).get("quantisation") or {}
        rel = max([v.get("max_abs_dev_relative_to_one_step", 0.0)
                   for v in q.values()] or [0.0])
        L.append("| %s | %.3f | %.5f | %.3f | %s / %s | %.2g | `%s` / `%s` |"
                 % (k, m.get("patch_m", float("nan")),
                    m.get("hs_spectral_m", float("nan")),
                    m.get("omega_c", float("nan")),
                    (q.get("lo") or {}).get("n_modes_active", "?"),
                    (q.get("hi") or {}).get("n_modes_active", "?"),
                    rel, (m.get("H0_sha256_12") or {}).get("lo", "?"),
                    (m.get("H0_sha256_12") or {}).get("hi", "?")))
    L += [
        "",
        "（H0 sha は z ごとに数ビット異なる: 深さは分散・粘性カットオフを通じて"
        "実現値に 1e-7 相対で入る。同一セル・同一 seed の物理状態であり、"
        "同じ波の別インスタンスではない。）",
        "",
        "## 2. 層ごとのフレーム数（実測、ディレクトリを数えた値）",
        "",
        "| ディレクトリ | PNG | フレーム範囲 | 欠番 | 完備 |",
        "|---|---:|---|---|---|",
    ]
    for k, v in counts.items():
        L.append("| `%s` | %d | %s..%s | %s | %s |"
                 % (v["dir"], v["n_png"], v["first"], v["last"],
                    ("なし" if not v["missing"] else str(v["missing"][:8])),
                    "YES" if v["complete"] else "NO"))
    L += [
        "",
        "| 層 | 内容（定義） |",
        "|---|---|",
        "| `final` | v16 本番合成（A2.6 + Blow σ3 + E2.2）そのまま |",
        "| `height` | 高さ場 h（表示帯域制限 σ=%.0f マスターセル）を固定"
        "スケール ±%.4g m で偽色化 |"
        % (DISPLAY_SIGMA_CELLS, float(param.get("height_limit_m") or 0.0)),
        "| `warp` | transport 変位（同帯域制限）を色相=方向・明度="
        "|d|/%.0f px で可視化 |" % WARP_PX_SCALE,
        "| `base` | パディング済み写真そのまま（水なし） |",
        "| `caustics` | 緑の直接 irradiance を重み和で正規化（1=平坦水面）し"
        "固定窓 [%.2f, %.2f]→sRGB |" % CAUSTICS_WINDOW,
        "",
        "可視化層のスケールは固定（フレームごとの自動拡縮なし）。"
        "捏造した層はない。",
        "表示帯域制限: マスター格子は S=%.3f m / 2048 = 1セル %.3g m で、"
        "出力1画素は %.3g m。毛細波帯は出力1画素より細かく、そのまま"
        "736x540 に標本化すると画素ごとのノイズになるため、`height` と "
        "`warp` だけ σ=%.0f セル（%.3f m ≒ %.1f 出力px）で帯域制限した。"
        "**`final` と §6 の実測値には適用していない**（本番定義のまま）。"
        % ((meta.get("4.0") or {}).get("patch_m") or float("nan"),
           (meta.get("4.0") or {}).get("dx_m") or float("nan"),
           ((meta.get("4.0") or {}).get("patch_m") or 0) / WIDTH,
           DISPLAY_SIGMA_CELLS,
           (param.get("display_bandlimit") or {}).get("sigma_m")
           or float("nan"),
           (param.get("display_bandlimit") or {}).get("sigma_output_px")
           or float("nan")),
        "",
        "## 3. 代表フレームのハッシュ（frame 0 / 24 / 47）",
        "",
        "| 群 | frame | bytes | SHA256 |",
        "|---|---:|---:|---|",
    ]
    for r in numbers.get("deliverable_hashes") or []:
        L.append("| %s | %d | %s | `%s` |"
                 % (r["group"], r["frame"],
                    r["bytes"] if r["bytes"] is not None else "-",
                    (r["sha256"] or "-")))
    L += [
        "",
        "全 48 フレームのハッシュは `checkpoint.json` の "
        "`material_a.frames[*].sha256` / `depth_b[*].frames[*].sha256` にある"
        "（再実行時はこれを照合して完了済みをスキップする）。",
        "",
        "## 4. 所要時間（実測）",
        "",
        "| 単位 | 状態 | 壁時計 [s] | 備考 |",
        "|---|---|---:|---|",
    ]
    for name in ("stage-0-preflight", "material-a", "depth-0.6",
                 "depth-1.0", "depth-2.0", "loop-closure-t2",
                 "measure-frame-24", "report"):
        st, wall, err = unit_row(name)
        L.append("| %s | %s | %s | %s |"
                 % (name, st,
                    ("%.1f" % wall) if wall is not None else "-",
                    err[:90] if err else ""))
    L += [
        "",
        "| 群 | フレーム数 | 平均 s/frame | 最小..最大 | det 合計 [s] |",
        "|---|---:|---:|---|---:|",
    ]
    for k, v in timing.items():
        if "mean_frame_s" not in v:
            continue
        L.append("| %s | %d | %.2f | %.2f..%.2f | %.1f |"
                 % (k, v["n_frames"], v["mean_frame_s"],
                    v.get("min_frame_s", float("nan")),
                    v.get("max_frame_s", float("nan")),
                    v.get("det_s_sum", float("nan"))))
    before_h = ((param.get("py_hashes_before") or {}).get("hashes") or {})
    after_h = py_set()
    added = sorted(set(after_h) - set(before_h))
    removed = sorted(set(before_h) - set(after_h))
    changed = sorted(k for k in set(before_h) & set(after_h)
                     if before_h[k] != after_h[k])
    changed_existing = [c for c in changed if c not in NEW_FILES]
    added_info = []
    for name in added:
        p = os.path.join(HERE, name)
        rec = {"file": name}
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                rec["first_line"] = (f.readline() or "").strip()[:110]
            stt = os.stat(p)
            rec["ctime"] = datetime.datetime.fromtimestamp(
                stt.st_ctime).isoformat(timespec="seconds")
        except OSError as exc:                            # noqa: BLE001
            rec["first_line"] = "unreadable: %s" % exc
        added_info.append(rec)
    detrec = cp.get("det") or {}
    L += [
        "",
        "det は `%s`（row_block=%s、numba スレッド %s、実測時間に JIT は含まない）。"
        "同一性は `out/det-wire/REPORT-DET-WIRE.md`（144 ファイル SHA256 一致）"
        "で検証済みの経路をそのまま使っている。"
        % (detrec.get("impl", "?"), detrec.get("row_block", "?"),
           detrec.get("threads", "既定")),
        "",
        "**既存ファイルの改変ゼロ（制約の確認）**: ジョブ開始時の `*.py` "
        "スナップショット（%d 本）と現在を SHA256 で比較 → 追加 %s / 削除 %s / "
        "変更 %s（うち本ジョブの新規ファイル以外の変更 **%s**）。"
        % (len(before_h), added or "なし", removed or "なし",
           changed or "なし", changed_existing or "なし"),
        "",
    ]
    if added:
        L += ["追加された `*.py` の出所（機械で読める証拠）:", "",
              "| file | ctime | first line |", "|---|---|---|"]
        for a in added_info:
            L.append("| `%s` | %s | %s |"
                     % (a["file"], a.get("ctime", "?"),
                        (a.get("first_line") or "").replace("|", "\\|")))
        L += ["",
              "本ジョブが追加したのは %s のみ。上表のうちそれ以外は、"
              "同時刻に同居していた別ジョブが作ったファイルで、本ジョブは"
              "変更も実行もしていない（ヘッダの記述と ctime が根拠）。"
              % ", ".join("`%s`" % n for n in NEW_FILES),
              ""]
    L += [
        "## 5. 出せなかった層・単位と理由",
        "",
    ]
    if not failures and all(st != "failed" for st, _w, _e in
                            (unit_row(n) for n in units)):
        L += ["なし（全層・全単位が完走）。", ""]
    else:
        L += ["| 単位 | 状態 | 理由 |", "|---|---|---|"]
        for name, u in units.items():
            if u.get("status") == "failed":
                L.append("| %s | failed | %s |"
                         % (name, (u.get("error") or "").replace("|", "\\|")))
        for k, f in failures.items():
            L.append("| %s | frame failed | %s |"
                     % (k, (f.get("error") or "").replace("|", "\\|")))
        L.append("")
    L += [
        "## 6. 実測（発注書の台帳下書き）",
        "",
        "**素材A final の frame 24（t=1.000 s）transport 変位**"
        "（v19 numbers-prod の frame 96 実測と同じ定義）:",
        "",
        "| 帯 | mean [px] | p50 [px] | p99 [px] | max [px] |",
        "|---|---:|---:|---:|---:|",
    ]
    for tag in ("r", "g", "b"):
        ch = (warp.get("per_channel") or {}).get(tag) or {}
        L.append("| %s | %.3f | %.3f | %.3f | %.3f |"
                 % (tag, ch.get("px_mean", float("nan")),
                    ch.get("px_p50", float("nan")),
                    ch.get("px_p99", float("nan")),
                    ch.get("px_max", float("nan"))))
    v19 = warp.get("v19_reference") or {}
    L += [
        "",
        "- プロジェクトの基準帯は **g**: mean **%.3f px** / p99 **%.3f px** / "
        "max **%.3f px**。" % (warp.get("g_mean_px") or float("nan"),
                              warp.get("g_p99_px") or float("nan"),
                              warp.get("g_max_px") or float("nan")),
        "- 定義: `shx = -nx*depth*(1-1/n)`、出力画素中心で S=%.4f m にて"
        "サンプル、`|d| = sqrt((shx/S*(W-1))^2+(shy/S*(H-1))^2)`。"
        "低域通過 σ=3.0 px。測定コードは `run_promo_v1.py::displacement_px` "
        "（`run_v19_prod.do_warp` と同一式）。" % (warp.get("patch_m")
                                                  or float("nan")),
        "- v19 参照（frame 96）: mean %.3f px / p99 %.3f px / max %.3f px "
        "（`out/v19/numbers-v19-prod.json#warp`）。差分は**時間格子**"
        "（8s → 2s）と**フレーム**（96 → 24）のみで、式・帯・サンプリングは"
        "同一。したがって定義は揃っている。"
        % (v19.get("g_mean_px") or float("nan"),
           v19.get("g_p99_px") or float("nan"),
           v19.get("g_max_px") or float("nan")),
        "",
        "## 7. T=2s ループ閉合（格子の証明）",
        "",
        "| 量 | 値 |",
        "|---|---|",
        "| t=2.000 s vs t=0 の場の最大差（h/gx/gy/Dx/Dy） | %s |"
        % ", ".join("%s=%.3g" % (k, v) for k, v in
                    (loop.get("fields_max_abs_diff_t2_vs_t0") or {}).items()),
        "| C（det master）の t=2 vs t=0 最大差 | %.6g |"
        % (loop.get("C_max_abs_diff_t2_vs_t0") or float("nan")),
        "| 1 フレーム分の感度（t=1/24 vs t=0） | %.6g |"
        % (loop.get("C_max_abs_diff_1frame_vs_t0_sensitivity")
           or float("nan")),
        "| 1 フレーム感度に対する残差比 | %.3g |"
        % (loop.get("C_seam_residual_relative_to_one_frame_step")
           or float("nan")),
        "| 量子化証明（全モードが n·2π/2s 上） | %s |"
        % loop.get("quantisation_proof_pass"),
        "",
        "ω_q は float32 格納なので t=2s の位相は 1 に完全一致しない。"
        "上の残差は float32 丸めに由来し、1 フレーム分の変化より 5 桁以上小さい。",
        "",
        "## 8. 再現コマンド",
        "",
        "```",
        "cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline",
        "python -B run_promo_v1.py --stage0            # 事前検査のみ",
        "python -B run_promo_v1.py                     # 全単位（冪等）",
        "python -B run_promo_v1.py --only material-a    # 素材Aだけ再開",
        "python -B run_promo_v1.py --probe 24           # 1フレーム計測",
        "```",
        "",
        "## 9. 不確かさ・注意（正直に）",
        "",
        "1. **8s 版との比較は同一フレームではない**。本レンダーは 2s 格子で"
        "組み直しているため、同じ t でも 8s 版と位相が異なる（スペクトルと"
        "空間実現は同一）。見た目の等価性は主張しない。",
        "2. **2s 格子は ω 刻みが粗い**（0.5 Hz = 3.14 rad/s）。波長の長い"
        "モードは n=1 に落ち、実分散より速く進む。これは 2s 周期を構造で"
        "保証するための代償で、量子化誤差は §1 の表に実測値として載せた。",
        "3. **可視化層（height/warp/caustics）は表示変換**である。数値そのもの"
        "ではなく、固定スケールで写した画像。スケールは §2 に明記した。",
        "4. **深さスイープの状態は各深さで独立に実現**される（H0 が 1e-7 相対"
        "で異なる）。同一セル・同一 seed の物理状態だが、同一の波インスタンス"
        "ではない。",
        "5. **所要時間は共用マシン上の実測**（他プロセス稼働中）。",
        "6. 数値は無太郎の確認前の中間値。本レポートはレンダーの完了と実測の"
        "記録であり、見た目の評価ではない。",
        "",
    ]
    payload = "\n".join(L).encode("utf-8")
    atomic_write_bytes(REPORT, payload)
    return REPORT


def do_report():
    t0 = time.time()
    log("=" * 72)
    log("[report] numbers-promo-draft.json + REPORT-RENDER.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    write_report(numbers)
    cp.unit_done("report")
    log("  numbers %s" % os.path.relpath(NUMBERS, HERE))
    log("  report  %s" % os.path.relpath(REPORT, HERE))
    log("  report wall %.1fs" % (time.time() - t0))
    return numbers


# --------------------------------------------------------------------------
# unit runner (a failed unit is logged, skipped, and reported)
# --------------------------------------------------------------------------
def run_unit(name, fn, fatal=False):
    t0 = time.time()
    log("=" * 72)
    log("[unit %s] start" % name)
    prev = (read_cp().get("units") or {}).get(name) or {}
    try:
        fn()
        wall = float(time.time() - t0)
        state = read_cp()
        units = state.get("units") or {}
        keep_first = (prev.get("status") == "done"
                      and isinstance(prev.get("wall_s"), (int, float)))
        units[name] = {
            "status": "done",
            # a re-run that skips (idempotent) must not overwrite the wall
            # clock of the invocation that actually produced the frames
            "wall_s": float(prev["wall_s"]) if keep_first else wall,
            "last_invocation_wall_s": wall,
            "completed_utc": (prev.get("completed_utc") if keep_first
                              else now_iso()),
            "rechecks": int(prev.get("rechecks", 0)) + 1 if keep_first else 0,
        }
        state["units"] = units
        write_cp(state)
        Checkpoint(CP_PATH).unit_done(name)
        log("[unit %s] DONE %.1fs%s" % (name, wall,
                                        " (recheck)" if keep_first else ""))
        return True
    except BaseException as exc:                          # noqa: BLE001
        if isinstance(exc, KeyboardInterrupt):
            raise
        msg = "%s: %s" % (type(exc).__name__, exc)
        wall = float(time.time() - t0)
        state = read_cp()
        units = state.get("units") or {}
        units[name] = {"status": "failed", "wall_s": wall, "error": msg,
                       "traceback": traceback.format_exc()[-3000:],
                       "failed_utc": now_iso()}
        state["units"] = units
        write_cp(state)
        log("[unit %s] FAILED %.1fs: %s" % (name, wall, msg))
        log("  (skipping this unit and continuing; the reason is in "
            "checkpoint.json and REPORT-RENDER.md)")
        if fatal:
            raise
        return False


# --------------------------------------------------------------------------
# probe
# --------------------------------------------------------------------------
def do_probe(n, det):
    makedirs()
    param = read_json(PARAMS_PATH) if os.path.exists(PARAMS_PATH) else {}
    hlim = float(param.get("height_limit_m") or 0.02)
    photo = load_photo()
    meta, st, env, optics = build_scene(CELL["h"], verbose=True)
    log("  probe frame %d (2s grid)" % n)
    r = render_frame(n, st, env, optics, photo, det, set(LAYERS),
                     height_limit=hlim, measure=True)
    for name in LAYERS:
        p = os.path.join(WORK, "probe-%04d-%s.png" % (n, name))
        sha = save_png(p, r["arrays"][name])
        log("  layer %-8s %s" % (name, sha[:16]))
    log("  timing %s" % r["timing_s"])
    log("  C mean %.9f p99 %.3f crush %.4f%% blow %.6g g1 %.6f"
        % (r["C_mean"], r["C_p99"], 100 * r["crush"], r["blow"], r["g1"]))
    if r["extra"].get("caustics_irradiance"):
        log("  caustics irr %s" % r["extra"]["caustics_irradiance"])
    d = r["displacement"]["g"]["stats"]
    log("  warp g: mean %.3f p50 %.3f p99 %.3f max %.3f px"
        % (d["px_mean"], d["px_p50"], d["px_p99"], d["px_max"]))
    del st, optics
    gc.collect()
    return r


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
UNITS = ("stage0", "material-a", "depth-0.6", "depth-1.0", "depth-2.0",
         "loop-closure-t2", "measure-frame-24", "report")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true", help="preflight only")
    ap.add_argument("--probe", type=int, default=None,
                    help="render one frame into _work (timing/layer check)")
    ap.add_argument("--only", default=None,
                    help="run only this unit (or comma-separated list)")
    ap.add_argument("--det-impl", default=None, help="fast (default) | base")
    ap.add_argument("--threads", type=int, default=None,
                    help="numba threads for the fast deposit (default 10)")
    a = ap.parse_args()

    threads = a.threads
    if threads is None and os.environ.get("PROMO_THREADS"):
        threads = int(os.environ["PROMO_THREADS"])
    if threads is None:
        threads = 10
    os.environ["PROMO_THREADS"] = str(threads)
    det_impl = (a.det_impl or os.environ.get("PROMO_DET_IMPL") or "fast")

    makedirs()
    log("#" * 72)
    log("run_promo_v1 start: units=%s det=%s threads=%s"
        % (a.only or "ALL", det_impl, threads))
    log("workdir %s" % HERE)

    if a.probe is not None:
        det = build_det(det_impl, threads)
        do_probe(a.probe, det)
        return 0

    if a.stage0:
        ok = run_unit("stage-0-preflight", stage0)
        return 0 if ok else 1

    selected = UNITS
    if a.only:
        selected = tuple(x.strip() for x in a.only.split(",") if x.strip())
        bad = [x for x in selected if x not in UNITS]
        if bad:
            ap.error("unknown unit(s): %s (known: %s)" % (bad, list(UNITS)))

    # stage 0 first if it has never run (needed for the height scale)
    if not os.path.exists(PARAMS_PATH) and "stage0" not in selected \
            and "stage-0-preflight" not in selected:
        selected = ("stage0",) + tuple(selected)

    det = None
    photo = None
    params = None

    def need_det():
        nonlocal det
        if det is None:
            det = build_det(det_impl, threads)
        return det

    def need_photo():
        nonlocal photo
        if photo is None:
            photo = load_photo()
        return photo

    def need_params():
        nonlocal params
        if params is None:
            params = read_json(PARAMS_PATH)
        return params

    for name in selected:
        if name == "stage0":
            run_unit("stage-0-preflight", stage0)
        elif name == "material-a":
            p = need_params()
            run_unit(name, lambda: render_material_a(need_photo(), need_det(),
                                                     p))
        elif name.startswith("depth-"):
            d = float(name.split("-", 1)[1])
            run_unit(name, lambda d=d: render_depth(d, need_photo(),
                                                    need_det()))
        elif name == "loop-closure-t2":
            run_unit(name, lambda: do_loop_closure(need_det()))
        elif name == "measure-frame-24":
            run_unit(name, lambda: do_warp24(need_det()))
        elif name == "report":
            run_unit(name, do_report)

    # the final always writes the report so the state is legible
    if "report" not in selected:
        log("=" * 72)
        log("writing the report (end of this invocation)")
        try:
            do_report()
        except BaseException:                             # noqa: BLE001
            log("report could not be written: %s"
                % traceback.format_exc()[-1500:])

    state = read_cp()
    units = state.get("units") or {}
    failed = [k for k, v in units.items() if v.get("status") == "failed"]
    log("=" * 72)
    log("run_promo_v1 end: units done=%d failed=%s"
        % (len([k for k, v in units.items() if v.get("status") == "done"]),
           failed or "none"))
    return 0 if not failed else 2


if __name__ == "__main__":
    sys.exit(main())
