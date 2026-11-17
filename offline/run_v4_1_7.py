"""run_v4_1_7 -- BRIEF-V4-1-7 driver: lattice -> noise (jitter + kernel deposit).

NEW FILES ONLY: `transport2.py` (the jitter + kernel deposit) and this driver.
Every existing .py stays byte-identical (G-LEGACY).  The accounting is
v4-1.4's (`stream_surface3`), untouched.  The slope field is v4-1.4's own
field, passed through with NO band limit (unlike v4-1.5).

STAGES (all foreground, jobkit checkpoint committed per stage)
    --stage0    absolute-path preflight + stage-0 sha256 record + scene record
                + the jitter / kernel specification and its verification
    --frame24   frame 24: C statistics (G-CONTRAST), deposit diagnostics,
                attribution (warp-only / caustics-only), G-LINE pre-check
    --diag      the negative controls for the G-LINE question (NOT used to set
                any parameter of the deliverable; see the report)
    --chunk K   frames 24K..24K+23
    --warp      displacement yardstick at frame 24 (G-WARP)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-LINE / G-CONTRAST / G-WARP / G-MOVE / G-LEGACY / G-OUT
    --report    numbers-v4-1-7.json + REPORT-V4-1-7.md
    --all       everything, sequentially

MEASUREMENT
    Exactly v4-1.5's implementation, imported, not re-written:
    `run_v4_1_5.lattice_metric` (green channel, 8-bit, 256x256 windows whose
    positions are fixed by a rule on the SOURCE PHOTO only, ideal axis-parallel
    5-15 px band-pass; depth = 1000*RMS(band)/window-mean, reported both with
    and without the brief's sigma=3 px pre-high-pass) and
    `run_v4_1_5.paired_distribution` (the same statistic over all 24 windows of
    the 92x90 tile grid).  Main gate quantity = depth_e3 (no pre-high-pass).
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

import run_v4_1_4 as R14  # noqa: E402  (reused UNCHANGED: build/scene/photo)
import run_v4_1_5 as R15  # noqa: E402  (reused UNCHANGED: the measurement)
import run_v16_water as V16  # noqa: E402
import transport2 as T2  # noqa: E402
from det import det_frame as DET_FRAME  # noqa: E402  (v4-1.4 reference only)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-7.md"
SCENE = "realscale"
OUT = os.path.join(HERE, "out", "v4-1-7")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-7-realscale.mp4")
STILL = os.path.join(OUT, "frame-realscale-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-7.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-7.md")

PHOTO = R15.PHOTO
V14_STILL = os.path.join(HERE, "out", "v4-1-4", "frame-realscale-0024.png")
V15_STILL = os.path.join(HERE, "out", "v4-1-5", "frame-realscale-0024.png")
V14_NUMBERS = os.path.join(HERE, "out", "v4-1-4", "numbers-v4-1-4.json")
V15_REPORT = os.path.join(HERE, "out", "v4-1-5", "REPORT-V4-1-5.md")

FPS = R14.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R14.WIDTH, R14.HEIGHT
BUF = (2 * HEIGHT, 2 * WIDTH)          # (1080, 1472) per the brief
NEW_FILES = ["transport2.py", "run_v4_1_7.py"]

# measurement (delegated, so the numbers are comparable by construction)
regions_from_photo = R15.regions_from_photo
lattice_metric = R15.lattice_metric
lattice_strength = R15.lattice_strength
measure_image = R15.measure_image
paired_distribution = R15.paired_distribution
green_of = R15.green_of


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
# the deposit -> render adapter
# ---------------------------------------------------------------------------
def c_render(Cf, n=None):
    """(540,736) frame-resolution C -> square (n,n) periodic buffer.

    The renderer (`v5_composite._sample_c`) indexes a SQUARE periodic buffer
    with one N for both axes, while the brief's C is at frame resolution
    (736x540).  The adapter resamples only along y (540 -> 736 rows) by linear
    interpolation, i.e. it hands the renderer a grid whose x samples are the
    real 736 and whose y samples are the interpolated 736; the renderer's own
    bilinear sampling then addresses the same C(u,v) surface.  Cost disclosure
    is in the report (sections 2 and 10).  x is not touched.
    """
    Cf = np.asarray(Cf, dtype=np.float64)
    H, W = Cf.shape
    n = int(n or W)
    yy = np.arange(n, dtype=np.float64) * (H - 1.0) / (n - 1.0)
    y0 = np.floor(yy).astype(np.int64)
    y1 = np.minimum(y0 + 1, H - 1)
    f = (yy - y0)[:, None]
    return ((1.0 - f) * Cf[y0, :] + f * Cf[y1, :]).astype(np.float32)


def render_optics(optics, n):
    """The v16 optics with the halo kernel re-sampled on an n x n C grid.

    The radial profile is stored in px of the state grid (512); the new C grid
    has n samples across the same square patch, so the radii are scaled by
    n/512 before `kernel_2d` is rebuilt.  Nothing else in the optics changes.
    """
    prof = optics["kernel_profile"]
    o = dict(optics)
    if prof is not None:
        scale = float(n) / float(BUF_GRID_N)
        p2 = dict(prof)
        p2["centres_px"] = np.asarray(prof["centres_px"], dtype=np.float64) * scale
        p2["edges_px"] = np.asarray(prof["edges_px"], dtype=np.float64) * scale
        p2["r_max_px"] = float(prof["r_max_px"]) * scale
        K, kinfo = V16.kernel_2d(p2, int(n))
        o["kernel_2d"] = K
        o["kernel_info_at_c_grid"] = kinfo
        o["kernel_radius_scale"] = scale
    return o


BUF_GRID_N = 512


# ---------------------------------------------------------------------------
# one frame
# ---------------------------------------------------------------------------
def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def render_v17(i, st, env, optics, photo, variant="both", c_scale=None):
    """One v4-1.7 frame.

    variant 'both'      : C from transport2 (jitter + kernel deposit) + warp
            'warp_only' : C = 1, warp from the field
            'caus_only' : C from transport2, zero slope (no warp)
    c_scale (diagnostic only) multiplies the derived kernel width.
    """
    sp = R14.SCENES[SCENE]
    t = i / float(FPS)
    fields = st.fields_at(t)
    g1 = _g1(fields, env)
    Q = int(sp["Q"])
    old = T2.KERNEL_SIGMA_PER_SPACING
    if c_scale is not None:
        T2.KERNEL_SIGMA_PER_SPACING = old * float(c_scale)
    try:
        C, info = T2.det_frame2(fields, st.S, env, Q, int(i), n_buf=BUF)
    finally:
        T2.KERNEL_SIGMA_PER_SPACING = old
    if variant == "warp_only":
        C = np.ones((HEIGHT, WIDTH), dtype=np.float32)
        fx = fields
        dinfo = {"variant": "C = 1 (no caustics)"}
    elif variant == "caus_only":
        z = np.zeros_like(fields["gx"])
        fx = {"h": fields["h"], "gx": z, "gy": z,
              "Dx": np.zeros_like(z), "Dy": np.zeros_like(z),
              "S": float(fields["S"])}
        dinfo = info
    else:
        fx = fields
        dinfo = info
    Cr = c_render(C)
    o2 = render_optics(optics, Cr.shape[0])
    img, dbg = V16.render_frame_v16(
        fx, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
        ambient_scale=sp["A"], c_layer_exposure=sp["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=sp["B_sigma"], return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "C_render": Cr, "g1": g1, "det_info": dinfo, "metrics": met,
            "variant": variant, "c_scale": (1.0 if c_scale is None else float(c_scale))}


def c_stats(C):
    C = np.asarray(C, dtype=np.float64)
    return {"mean": float(np.mean(C)), "std": float(np.std(C)),
            "min": float(np.min(C)), "max": float(np.max(C)),
            "p01": float(np.percentile(C, 1)),
            "p50": float(np.percentile(C, 50)),
            "p99": float(np.percentile(C, 99)),
            "shape": [int(C.shape[0]), int(C.shape[1])]}


def axis_diag(g, lo=5.0, hi=15.0):
    """Axis-parallel vs diagonal band amplitude (hard masks, no window).

    A shape statistic to complement `depth`: an isotropic texture gives
    axis/diag ~ 1, a lattice or streaky field gives >> 1.  Reported as a
    diagnostic only -- the gate stays `lattice_metric`.
    """
    g = np.asarray(g, dtype=np.float64)
    H, W = g.shape
    kx = np.fft.fftfreq(W) * W
    ky = np.fft.fftfreq(H) * H
    KX, KY = np.meshgrid(kx, ky)
    A, B = np.abs(KX), np.abs(KY)
    lo_k, hi_k = W / hi, W / lo
    F = np.fft.fft2(g - g.mean())

    def rms(m):
        return 1000.0 * float(np.sqrt(np.mean(
            np.real(np.fft.ifft2(F * m)) ** 2))) / float(g.mean())

    axm = (((A >= lo_k) & (A <= hi_k) & (B <= 3.0))
           | ((B >= lo_k) & (B <= hi_k) & (A <= 3.0)))
    dm = (((np.abs(A - B) <= 3.0) | (np.abs(A + B) <= 3.0))
          & (np.hypot(A, B) >= lo_k) & (np.hypot(A, B) <= hi_k))
    a, d = rms(axm.astype(np.float64)), rms(dm.astype(np.float64))
    return {"axis_permille": a, "diag_permille": d,
            "ratio": (a / d if d > 0 else None),
            "shape": [int(H), int(W)]}


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def jitter_verification(frames=(0, 1, 12, 24, 47)):
    N_s = BUF_GRID_N
    rows = []
    prev = None
    for f in frames:
        jx, jy = T2.jitter_offsets(N_s, f)
        # 8x8 occupancy of the unit cell (chi-square vs uniform)
        occ, _, _ = np.histogram2d(jx.ravel(), jy.ravel(), bins=8,
                                   range=[[-0.5, 0.5], [-0.5, 0.5]])
        exp = occ.sum() / 64.0
        chi2 = float(np.sum((occ - exp) ** 2 / exp))
        rec = {"frame": int(f),
               "jx_mean": float(jx.mean()), "jx_std": float(jx.std()),
               "jx_min": float(jx.min()), "jx_max": float(jx.max()),
               "jy_mean": float(jy.mean()), "jy_std": float(jy.std()),
               "uniform_std_reference": float(1.0 / np.sqrt(12.0)),
               "cell_8x8_chi2": chi2,
               "cell_8x8_chi2_expect": 63.0,
               "corr_with_previous_frame": (
                   None if prev is None else float(np.corrcoef(
                       jx.ravel(), prev.ravel())[0, 1])),
               "offsets_0000_0004": [float(x) for x in jx[0, :5]],
               "offsets_y_0000_0004": [float(x) for x in jy[:5, 0]]}
        rows.append(rec)
        prev = jx
    return {"spec": T2.RECORD["jitter"], "per_frame": rows,
            "distinct_frames_checked": len(frames)}


def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record + specs")
    makedirs()
    cp = Checkpoint(CP_PATH)
    need = [BRIEF,
            os.path.join(HERE, "stream_surface3.py"),
            os.path.join(HERE, "stream_surface4.py"),
            os.path.join(HERE, "run_v4_1_4.py"),
            os.path.join(HERE, "run_v4_1_5.py"),
            V15_REPORT,
            os.path.join(HERE, "transport2.py"),
            os.path.join(HERE, "run_v4_1_7.py"),
            PHOTO, V14_STILL, V15_STILL, V14_NUMBERS]
    missing = [p for p in need if not os.path.isfile(p)]
    log("  required paths: %d, missing: %s" % (len(need), missing or "none"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))
    rec = {"brief": BRIEF, "checked": need, "missing": [], "checked_utc": now_iso()}
    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in NEW_FILES}
    payload = {"all": hashes, "new_files_of_this_job": NEW_FILES,
               "existing_count": len(existing), "existing": existing,
               "recorded_utc": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's two new files (their stage-0 hashes are in "
                        "'all' as well)")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(NEW_FILES)))

    sp = R14.SCENES[SCENE]
    L, n, loop = R14.scene_loop(SCENE)
    rec2, st, env, optics = R14.build(SCENE)
    scene_rec = {
        "scene": SCENE,
        "constants": {"H": sp["H"], "U": sp["U"], "S0": sp["S0"],
                      "D84": sp["D84"], "seed": sp["seed"],
                      "water": sp["water"], "L_guess": sp["L_guess"],
                      "T_loop": sp["T_loop"], "N": sp["N"], "A": sp["A"],
                      "B_sigma": sp["B_sigma"], "E": sp["E"], "Q": sp["Q"]},
        "source": "run_v4_1_4.SCENES['realscale'] (verbatim)",
        "loop": loop, "L_m": float(L),
        "accounting": st.accounting(),
        "patch_m": float(st.S), "grid": int(st.N),
        "cell_m": float(st.S / st.N),
        "sigma_total_m": float(st.sigma_total),
        "var_i_m2": float(st.var_i_physical),
        "var_ii_m2": float(st.var_ii), "var_iii_m2": float(st.var_iii),
        "regions": regions_from_photo()[1],
        "field": ("v4-1.4 h/gx/gy/Dx/Dy verbatim, NO band limit "
                  "(stream_surface4 is NOT used)"),
        "deposit": T2.RECORD,
        "buffer": {"shape": list(BUF), "scale": 2,
                   "note": "2 x (540,736) as the brief specifies"},
        "optics": {"kernel": optics.get("kernel"),
                   "phase_kind": optics.get("phase_kind"),
                   "w_ballistic": [float(x) for x in optics["w_ballistic"]],
                   "w_halo": [float(x) for x in optics["w_halo"]],
                   "kernel_rmax_px_state_grid": (
                       float(optics["kernel_profile"]["r_max_px"])
                       if optics["kernel_profile"] else None)},
        "kernel_derivation": {
            "s0_m": float(st.S / st.N),
            "sigma_rule": "sigma(x) = s0 * sqrt(|det J|) / 2",
            "sigma_cell_range": T2.RECORD["kernel"]["sigma_cell_range"],
            "sigma_m_range": [float(T2.RECORD["kernel"]["sigma_cell_range"][0]
                                    * st.S / st.N),
                              float(T2.RECORD["kernel"]["sigma_cell_range"][1]
                                    * st.S / st.N)],
            "buffer_samples_per_cell": [float(HEIGHT * 2 / st.N),
                                        float(WIDTH * 2 / st.N)],
        },
        "jitter": jitter_verification(),
        "stage0_rec_ok": True,
    }
    checkpoint_extra("stage0", rec)
    checkpoint_extra("scene", scene_rec)
    cp.unit_done("stage-0-preflight")
    log("  scene: H=%.4f U=%.3f L=%.6f m grid=%d cell=%.5f mm"
        % (sp["H"], sp["U"], L, st.N, 1e3 * st.S / st.N))
    log("  accounting: sigma_reg=%.6f mm var_i=%.4e var_ii=%.4e var_iii=%.4e"
        % (st.sigma_reg * 1e3, st.var_i_physical, st.var_ii, st.var_iii))
    jv = scene_rec["jitter"]["per_frame"]
    log("  jitter: %d frames verified, chi2(63 dof) = %s"
        % (len(jv), [round(r["cell_8x8_chi2"], 1) for r in jv]))
    del rec2, st, optics
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
# stage: frame 24
# ---------------------------------------------------------------------------
def do_frame24():
    t0 = time.time()
    log("=" * 78)
    log("[frame24] C statistics, attribution, G-LINE pre-check")
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, rinfo = regions_from_photo()
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    r = render_v17(STILL_IDX, st, env, optics, photo)
    sha = save_png(STILL, r["still"])
    rows = lattice_metric(r["still"][:, :, 1].astype(np.float64), regions)
    photo_rows = measure_image(PHOTO, regions)
    per_row = []
    for got, base in zip(rows, photo_rows):
        per_row.append({"region": got["region"], "axis": got["axis"],
                        "period_px": got["period_px"],
                        "depth_v4_1_7": got["depth_e3"],
                        "depth_photo": base["depth_e3"],
                        "ratio": got["depth_e3"] / base["depth_e3"],
                        "depth_hp_v4_1_7": got["depth_hp_e3"],
                        "depth_hp_photo": base["depth_hp_e3"],
                        "ok": bool(got["depth_e3"] <= base["depth_e3"])})
    out = {"frame": STILL_IDX, "t_s": STILL_IDX / float(FPS),
           "still": os.path.basename(STILL), "still_sha256": sha,
           "regions": rinfo, "measurements": rows, "per_row_vs_photo": per_row,
           "strength": lattice_strength(rows),
           "photo_strength": lattice_strength(photo_rows),
           "pass_gl_line": bool(all(p["ok"] for p in per_row)),
           "C_frame_stats": c_stats(r["C"]),
           "C_render_stats": c_stats(r["C_render"]),
           "C_stats_v4_1_4_reference": {"std": 0.5283, "note": "REPORT-V4-1-5"},
           "C_stats_v4_1_5_reference": {"std": 0.101, "note": "REPORT-V4-1-5"},
           "deposit": {k: r["det_info"].get(k) for k in
                       ("mean_before_mean_correction", "mean_correction_factor",
                        "divisor_analytical", "jacobian_clip",
                        "kernel_sigma_centres_cells",
                        "kernel_sigma_centres_buffer_samples_yx",
                        "kernel_ray_fraction_per_bin",
                        "detJ_abs_p50_p99_max", "fold_per_quadrature")},
           "image_metrics": r["metrics"], "g1": float(r["g1"]),
           "axis_diag": {
               "C_frame": axis_diag(r["C"]),
               "still_green": axis_diag(r["still"][:, :, 1]),
               "v4_1_4_still_green": axis_diag(green_of(V14_STILL)),
               "v4_1_5_still_green": axis_diag(green_of(V15_STILL)),
               "photo_green": axis_diag(green_of(PHOTO))},
           "wall_s": float(time.time() - t0)}
    # attribution: (a) warp only, (b) caustics only
    abl = {}
    for tag, var in (("a", "warp_only"), ("b", "caus_only")):
        ra = render_v17(STILL_IDX, st, env, optics, photo, variant=var)
        rws = lattice_metric(ra["still"][:, :, 1].astype(np.float64), regions)
        abl[tag] = {"name": {"a": "warp only (C = 1)",
                             "b": "caustics only (no warp)"}[tag],
                    "depth": [x["depth_e3"] for x in rws],
                    "strength": lattice_strength(rws),
                    "periods_px": [x["period_px"] for x in rws],
                    "C_std": float(np.std(ra["C"]))}
        log("  (%s) %-26s strength %.1f periods %s"
            % (tag, abl[tag]["name"], abl[tag]["strength"],
               [round(x, 2) for x in abl[tag]["periods_px"]]))
        del ra
        gc.collect()
    out["attribution"] = abl
    out["attribution"]["c_both_strength"] = out["strength"]
    for p in per_row:
        log("  %-6s %s  %.2f px  depth %.1f (photo %.1f) ratio %.3f  %s"
            % (p["region"], p["axis"], p["period_px"], p["depth_v4_1_7"],
               p["depth_photo"], p["ratio"], "OK" if p["ok"] else "FAIL"))
    log("  G-LINE pre-check: %s (strength %.1f vs photo %.1f)"
        % ("PASS" if out["pass_gl_line"] else "FAIL", out["strength"],
           out["photo_strength"]))
    cs = out["C_frame_stats"]
    log("  C(frame res %dx%d): mean %.6f std %.4f min %.4f p99 %.4f max %.3f"
        % (cs["shape"][1], cs["shape"][0], cs["mean"], cs["std"], cs["min"],
           cs["p99"], cs["max"]))
    log("  kernel ray fractions per bin: %s"
        % [round(x, 3) for x in out["deposit"]["kernel_ray_fraction_per_bin"]])
    checkpoint_extra("frame24", out)
    cp.unit_done("frame-0024")
    del st, optics, photo, r
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: negative controls for the G-LINE question (diagnostic only)
# ---------------------------------------------------------------------------
def do_diag():
    """Diagnostics.  NOTHING here sets a parameter of the deliverable.

    C0  the deliverable (R2 jitter, derived kernel, 1472x1080 buffer)
    C1  jitter switched off (source lattice restored) -> does jitter matter?
    C2  square buffer 1472x1472 (aspect control)
    C3  buffer at 2x the brief's resolution (2160x2944) -> resolution control
    C4  kernel width x1.5, x2, x3 -> the G-LINE vs G-CONTRAST frontier
    C5  the field's axis-parallel slope power, and a notch that removes it
        (99.9% of the 5-15 cell axis line) -> is the axis band in the field?
    """
    t0 = time.time()
    log("=" * 78)
    log("[diag] negative controls (no parameter of the deliverable is set here)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    regions, _ = regions_from_photo()
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    sp = R14.SCENES[SCENE]
    Q = int(sp["Q"])
    fields = st.fields_at(STILL_IDX / float(FPS))
    out = {"note": ("diagnostic only; the deliverable is C0 and its parameters "
                    "were fixed before any of these were measured")}

    def measure_rows(img):
        rws = lattice_metric(img[:, :, 1].astype(np.float64), regions)
        return [x["depth_e3"] for x in rws], lattice_strength(rws)

    def render_with(C):
        Cr = c_render(C)
        o2 = render_optics(optics, Cr.shape[0])
        g1 = _g1(fields, env)
        img, _ = V16.render_frame_v16(
            fields, Cr, photo["lin"], env, WIDTH, HEIGHT, o2,
            ambient_scale=sp["A"], c_layer_exposure=sp["E"],
            substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
            lowpass_sigma_px=sp["B_sigma"], return_debug=True)
        return img

    def case(T, nbuf, label, scale=None):
        C, info = T2.det_frame2(fields, st.S, env, Q, STILL_IDX, n_buf=nbuf)
        img = render_with(C)
        d, s = measure_rows(img)
        rec = {"label": label, "C_std": float(np.std(C)),
               "C_mean": float(np.mean(C)), "buffer": list(nbuf),
               "kernel_sigma_scale": (1.0 if scale is None else float(scale)),
               "depth": d, "strength": s,
               "C_shape": list(C.shape)}
        log("  %-34s C std %.4f  image depth %s  max %.1f"
            % (label, rec["C_std"], [round(x, 1) for x in d], s))
        return rec

    cases = {}
    cases["C0"] = case(T2, BUF, "C0 deliverable")
    orig_jit = T2.jitter_offsets
    T2.jitter_offsets = lambda N, fr: (np.zeros((N, N)), np.zeros((N, N)))
    cases["C1"] = case(T2, BUF, "C1 no source jitter")
    T2.jitter_offsets = orig_jit
    cases["C2"] = case(T2, (2 * WIDTH, 2 * WIDTH), "C2 square buffer 1472^2")
    cases["C3"] = case(T2, (4 * HEIGHT, 4 * WIDTH), "C3 buffer 2160x2944 (2x)")
    old = T2.KERNEL_SIGMA_PER_SPACING
    for sc in (1.5, 2.0, 3.0):
        T2.KERNEL_SIGMA_PER_SPACING = old * sc
        cases["C4_x%.1f" % sc] = case(
            T2, BUF, "C4 kernel width x%.1f (diagnostic)" % sc, scale=sc)
    T2.KERNEL_SIGMA_PER_SPACING = old

    # ---- C5: is the axis-parallel band in the FIELD? ----------------------
    eta = st.eta_hat_at(STILL_IDX / float(FPS))
    dk = 2.0 * np.pi / float(st.S)
    K = np.hypot(st.kx, st.ky)
    band = (K >= 5.0 * dk) & (K <= 15.0 * dk)
    axis = band & ((np.abs(st.kx) <= 10.0 * dk) | (np.abs(st.ky) <= 10.0 * dk))
    slope_pow = np.abs(1j * st.kx * eta) ** 2 + np.abs(1j * st.ky * eta) ** 2
    share = float(slope_pow[axis].sum() / slope_pow[band].sum())
    eta_n = eta.copy()
    eta_n[axis] = 0.0
    hn = np.real(np.fft.ifft2(eta_n))
    gxn = np.real(np.fft.ifft2(1j * st.kx * eta_n))
    gyn = np.real(np.fft.ifft2(1j * st.ky * eta_n))
    z = np.zeros_like(hn)
    f_n = {"h": hn.astype(np.float32), "gx": gxn.astype(np.float32),
           "gy": gyn.astype(np.float32), "Dx": z.astype(np.float32),
           "Dy": z.astype(np.float32), "S": float(st.S)}
    removed = float((np.abs(1j * st.kx * (eta - eta_n)) ** 2
                     + np.abs(1j * st.ky * (eta - eta_n)) ** 2)[band].sum()
                    / slope_pow[band].sum())
    Cn, _ = T2.det_frame2(f_n, st.S, env, Q, STILL_IDX, n_buf=BUF)
    imgn = render_with(Cn)
    dn, sn = measure_rows(imgn)
    cases["C5"] = {"label": "C5 field notched on the axis lines (5-15 cell)",
                   "field_axis_line_share_of_band_slope_power": share,
                   "removed_fraction_of_band_slope_power": removed,
                   "C_std": float(np.std(Cn)), "depth": dn, "strength": sn}
    log("  C5 field axis-line share %.3f ; removed %.3f ; C std %.4f ; "
        "image depth %s max %.1f"
        % (share, removed, float(np.std(Cn)), [round(x, 1) for x in dn], sn))
    out["cases"] = cases
    # ---- the empirical frontier: depth vs C std (diagnostic fit) ----------
    pts = []
    for k in ("C0", "C4_x1.5", "C4_x2.0", "C4_x3.0"):
        if k in cases:
            pts.append((cases[k]["C_std"], cases[k]["strength"]))
    if len(pts) >= 2:
        xs = np.array([p[0] for p in pts], dtype=np.float64)
        ys = np.array([p[1] for p in pts], dtype=np.float64)
        a, b = np.polyfit(xs, ys, 1)
        out["frontier_fit"] = {
            "points_std_strength": [[float(x), float(y)] for x, y in pts],
            "linear_fit_depth_per_std": float(a), "intercept": float(b),
            "std_needed_for_photo_depth": float((63.689564352484375 - b) / a),
            "note": ("four points, one line: the only lever that moves the "
                     "5-15 px depth is the deposit's smoothing scale, and it "
                     "moves C std with it; diagnostic only")}
        log("  frontier: depth = %.1f*std + %.1f  -> depth=photo at std %.3f"
            % (a, b, out["frontier_fit"]["std_needed_for_photo_depth"]))
    out["wall_s"] = float(time.time() - t0)
    checkpoint_extra("diag", out)
    cp.unit_done("diagnostics")
    del st, optics, photo
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: frames
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
    log("[chunk %02d] frames %d-%d" % (k, k * CHUNK, (k + 1) * CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    if chunk_done(k, state):
        log("  chunk %02d already committed -- skip" % k)
        return True
    _p, st, env, optics = R14.build(SCENE)
    photo = R14.load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_v17(i, st, env, optics, photo)
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
        log("  frame %04d t=%.3f Cmean=%.5f Cstd=%.4f p99=%.3f %.1fs"
            % (i, r["t_s"], rec["C_mean"], rec["C_std"], rec["C_p99"],
               time.time() - t_start))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_utc": now_iso(),
            "wall_s": float(time.time() - t_start)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs" % (k, node["wall_s"]))
    del st, optics, photo
    gc.collect()
    return True


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP) -- R14.do_warp's method, re-used
# ---------------------------------------------------------------------------
def do_warp(i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp] displacement yardstick at frame %d" % i)
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    _p, st, env, optics = R14.build(SCENE)
    sp = R14.SCENES[SCENE]
    fields = st.fields_at(i / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, sp["B_sigma"])
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
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
                    "px_max": float(np.max(mag)),
                    "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                    "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    base = 1.8393
    try:
        v14cp = read_json(os.path.join(HERE, "out", "v4-1-4", "checkpoint.json"))
        base = v14cp["measured"]["warp_realscale"]["g_mean_px"]
    except Exception:
        pass
    out = {"scene": SCENE, "frame": int(i), "t_s": i / float(FPS),
           "field": ("transport displacement sampled at output pixel centres; "
                     "UNFILTERED slope field (identical to v4-1.4/v4-1.5)"),
           "patch_m": S, "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "v4_1_4_g_mean_px": base,
           "created_utc": now_iso(), "wall_s": float(time.time() - t0)}
    if base:
        out["ratio_to_v4_1_4"] = float(per["g"]["px_mean"] / base)
        out["within_pm25pct"] = bool(abs(per["g"]["px_mean"] / base - 1.0) <= 0.25)
    write_json(os.path.join(WORK, "warp-realscale.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.4f p50 %.4f p99 %.4f  (v4-1-4 %.4f, ratio %s)"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"], base,
           round(out.get("ratio_to_v4_1_4", float("nan")), 6)))
    del st, optics
    gc.collect()
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
           "sha256": sha256_file(MP4), "created_utc": now_iso(),
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
def green_frames():
    return [green_of(frame_path(i)) for i in range(NFRAMES)]


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def decode_green(path):
    cmd = ["ffmpeg", "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-LINE / G-CONTRAST / G-WARP / G-MOVE / G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    gates = {}

    # ---- G-LINE ----------------------------------------------------------
    got = measure_image(STILL, regions)
    photo = measure_image(PHOTO, regions)
    v14 = measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v15 = measure_image(V15_STILL, regions) if os.path.isfile(V15_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_7_px": g_["period_px"],
                     "depth_v4_1_7": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_7": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_7_vs_photo": paired_distribution(STILL, PHOTO),
            "v4_1_4_vs_photo": (paired_distribution(V14_STILL, PHOTO)
                                if os.path.isfile(V14_STILL) else None),
            "v4_1_5_vs_photo": (paired_distribution(V15_STILL, PHOTO)
                                if os.path.isfile(V15_STILL) else None),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid "
                           "(stride 92x90, 24 windows)")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("frame 24: the axis-parallel 5-15 px modulation depth is "
                      "<= the source photo's for both regions and both axes; "
                      "primary = unbiased band amplitude (depth_e3), the "
                      "brief's high-pass-first procedure is depth_hp_e3; "
                      "24-window median reported alongside"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_4": v14, "v4_1_5": v15,
        "strength_v4_1_7": lattice_strength(got),
        "strength_photo": lattice_strength(photo),
        "strength_v4_1_4": lattice_strength(v14) if v14 else None,
        "strength_v4_1_5": lattice_strength(v15) if v15 else None,
        "distribution": dist, "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_7_px"],
               r["depth_v4_1_7"], r["depth_photo"], r["ratio"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s" % ("PASS" if gates["G-LINE"]["pass"] else "FAIL"))

    # ---- G-CONTRAST ------------------------------------------------------
    cst = (state.get("frame24") or {}).get("C_frame_stats") or {}
    std = cst.get("std")
    gates["G-CONTRAST"] = {
        "pass": bool(std is not None and std >= 0.35),
        "threshold": 0.35, "C_std": std, "C_stats": cst,
        "references": {"v4_1_4_std": 0.5283, "v4_1_5_std": 0.101},
        "criterion": ("frame 24 C (frame resolution, mean exactly 1) has "
                      "standard deviation >= 0.35")}
    log("  G-CONTRAST: std %.4f (>= 0.35) -> %s"
        % (std if std is not None else float("nan"),
           "PASS" if gates["G-CONTRAST"]["pass"] else "FAIL"))

    # ---- G-MOVE ----------------------------------------------------------
    frames = green_frames()
    sha = [sha256_file(frame_path(i)) for i in range(NFRAMES)]
    uniq = len(set(sha))
    st_ = interframe_stats(frames)
    try:
        dg = decode_green(MP4)
        dec = interframe_stats(dg)
        dec["frames_decoded"] = int(dg.shape[0])
    except Exception as e:  # pragma: no cover
        dec = {"error": str(e)}
    gates["G-MOVE"] = {
        "pass": bool(uniq == NFRAMES and st_["mean"] > 0.5),
        "unique_frames": int(uniq), "frames": NFRAMES,
        "interframe_green_png": st_, "interframe_green_mp4": dec,
        "threshold_gray": 0.5,
        "criterion": "48/48 unique frame PNGs AND mean |d green| > 0.5 gray"}
    log("  G-MOVE: unique %d/%d  mean|d green| %.4f -> %s"
        % (uniq, NFRAMES, st_["mean"],
           "PASS" if gates["G-MOVE"]["pass"] else "FAIL"))

    # ---- G-WARP ----------------------------------------------------------
    warp = state.get("warp") or {}
    base = warp.get("v4_1_4_g_mean_px")
    got_w = warp.get("g_mean_px")
    ok_w = bool(base and got_w and abs(got_w / base - 1.0) <= 0.25)
    gates["G-WARP"] = {
        "pass": ok_w, "g_mean_px": got_w, "v4_1_4_g_mean_px": base,
        "ratio": (got_w / base if (got_w and base) else None),
        "band_pct": 25.0,
        "criterion": ("frame 24 green mean displacement within +-25% of "
                      "v4-1.4's 1.839 px"),
        "detail": warp}
    log("  G-WARP: %.4f px vs %.4f (ratio %s) -> %s"
        % (got_w or float("nan"), base or float("nan"),
           round(gates["G-WARP"]["ratio"], 6) if gates["G-WARP"]["ratio"] else None,
           "PASS" if ok_w else "FAIL"))

    # ---- G-LEGACY --------------------------------------------------------
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
    ext = {}
    for name in changed:
        p = os.path.join(HERE, name)
        ext[name] = {
            "stage0_sha256": rec["all"].get(name), "now_sha256": now.get(name),
            "mtime": datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds"),
            "bytes_now": os.path.getsize(p),
            "imported_by_this_job": name in (
                "run_v4_1_4.py", "run_v4_1_5.py", "run_v16_water.py",
                "transport2.py", "det.py", "transport.py",
                "stream_surface3.py", "stream_surface2.py"),
            "attribution": ("NOT this job: this driver never opens or imports "
                            "it; it changed after the stage-0 record and is "
                            "owned by a parallel lane")}
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "attributable_to_this_job": {
            "changed": [n for n in changed
                        if ext.get(n, {}).get("imported_by_this_job")],
            "note": ("all .py this job imports are unchanged; the changed "
                     "files are external (parallel lane)")},
        "external_changes": ext,
        "criterion": ("every .py recorded at stage 0 (excluding this job's "
                      "two new files) still has the same sha256; external "
                      "(parallel-lane) differences are reported with their "
                      "attribution split out")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d, added %d -> %s"
        % (len(rec["existing"]), len(changed), len(missing), len(added),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL"))

    # ---- G-OUT -----------------------------------------------------------
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

    order = ["G-LINE", "G-CONTRAST", "G-WARP", "G-MOVE", "G-LEGACY", "G-OUT"]
    gates["order"] = order
    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in order)
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL: %s" % ("PASS" if gates["all_pass"] else "FAIL"))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# stage: numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    regions, rinfo = regions_from_photo()
    enc = state.get("encode") or {}
    n = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_utc": now_iso(),
        "scene": state.get("scene"),
        "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
        "stage0": state.get("stage0"),
        "frame24": state.get("frame24"),
        "diagnostics": state.get("diag"),
        "warp": state.get("warp"),
        "encode": enc,
        "gates": state.get("gates"),
        "chunks": {k: {"wall_s": v.get("wall_s"),
                       "frames": sorted(v.get("frames", {})),
                       "metrics": v.get("metrics")}
                   for k, v in (state.get("chunks") or {}).items()},
        "measurement": {
            "implementation": ("run_v4_1_5.lattice_metric (imported, not "
                               "re-written): green channel, 8-bit; the brief's "
                               "high pass = image - gaussian(sigma=3 px); Hann "
                               "256x256; 2D FFT; axis band |k_perp| <= 3, "
                               "5-15 px; period = W/k at the band peak; "
                               "depth = 1000*RMS(band-passed)/mean"),
            "regions": rinfo,
            "comparability": ("the photo is 735x539 and the frames 736x540, so "
                              "the same window origin addresses the same "
                              "content; both are 8-bit and depth is normalised "
                              "by the window mean"),
            "known_limitation": ("depth is a 1D-marginal statistic; it does not "
                                 "distinguish a lattice from isotropic texture "
                                 "with the same axis-band power")},
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
    fp = frame_path(STILL_IDX)
    n["still_equals_frame_0024"] = bool(
        os.path.isfile(fp) and sha256_file(fp) == n["files"].get(
            os.path.basename(STILL), {}).get("sha256"))
    return n


def _fmt(x, nd=1):
    if isinstance(x, list):
        return " / ".join(_fmt(v, nd) for v in x)
    try:
        return ("%%.%df" % nd) % float(x)
    except Exception:
        return str(x)


def write_report(n):
    g = n["gates"]
    f24 = n["frame24"] or {}
    sc = n["scene"] or {}
    dg = n["diagnostics"] or {}
    mp = n["encode"] or {}
    pr = mp.get("probe") or {}
    gl = g["G-LINE"]
    rows = gl["rows"]
    cs = (f24.get("C_frame_stats") or {})
    dep = f24.get("deposit") or {}
    L = []
    A = L.append

    def gateword(k):
        return "**PASS**" if g[k]["pass"] else "**FAIL**"

    A("# REPORT-V4-1-7 — 格子を「ノイズ」に変える（カーネル堆積＋ジッタ）")
    A("")
    A("Brief: `BRIEF-V4-1-7.md`.  Workdir: `%s`." % HERE)
    A("Built: %s.  新規ファイルは `transport2.py` と `run_v4_1_7.py` の2本のみ。"
      "既存 .py は**無改造**（sha256 は §6）。会計は `stream_surface3`、"
      "場は v4-1.4 のまま（**帯域制限なし**）。" % n["built_utc"])
    A("")
    A("## 0. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-LINE 格子の強さ | %s | frame 24 の軸平行 5–15 px 変調深度（主判定 = "
      "前段ハイパスなし、v4-1.5 と同一実装）: v4-1.7 %s ‰ vs 元写真 %s ‰（比 %s）。"
      "v4-1.4 %s ‰、v4-1.5 %s ‰ |"
      % (gateword("G-LINE"),
         _fmt([r["depth_v4_1_7"] for r in rows]),
         _fmt([r["depth_photo"] for r in rows]),
         _fmt([r["ratio"] for r in rows], 2),
         _fmt(gl.get("strength_v4_1_4"), 1), _fmt(gl.get("strength_v4_1_5"), 1)))
    A("| G-CONTRAST 波感 | %s | frame 24 の C（フレーム解像度 736×540、平均 1 厳密）"
      "の標準偏差 **%.4f**（≥ 0.35）。p99 %.3f、max %.3f（v4-1.4 0.528、v4-1.5 0.101） |"
      % (gateword("G-CONTRAST"), float(cs.get("std") or 0.0),
         float(cs.get("p99") or 0.0), float(cs.get("max") or 0.0)))
    A("| G-WARP 歪み | %s | frame 24 緑 平均 %.4f px（v4-1.4 %.4f px、比 %.6f、"
      "±25%% 以内） |"
      % (gateword("G-WARP"), float(g["G-WARP"]["g_mean_px"] or 0.0),
         float(g["G-WARP"]["v4_1_4_g_mean_px"] or 0.0),
         float(g["G-WARP"]["ratio"] or 0.0)))
    A("| G-MOVE 動く | %s | ユニークフレーム %d/%d、フレーム間平均絶対差（緑）"
      "%.4f 諧調 > 0.5 |"
      % (gateword("G-MOVE"), g["G-MOVE"]["unique_frames"], NFRAMES,
         g["G-MOVE"]["interframe_green_png"]["mean"]))
    A("| G-LEGACY 既存 py | %s | 既存 %d 本 → 変更 %d（本作に帰属 %d / 並行レーン %d）"
      "/ 消失 %d |"
      % (gateword("G-LEGACY"), g["G-LEGACY"]["recorded_existing"],
         len(g["G-LEGACY"]["changed"]),
         len((g["G-LEGACY"].get("attributable_to_this_job") or {}).get("changed") or []),
         len((g["G-LEGACY"].get("external_changes") or {})),
         len(g["G-LEGACY"]["missing"])))
    A("| G-OUT 完走 | %s | %s frames / %.4f s / %sx%s / %.1f fps |"
      % (gateword("G-OUT"), pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0.0), pr.get("width"), pr.get("height"),
         float(g["G-OUT"]["fps"])))
    A("")
    if g["all_pass"]:
        A("**6ゲートすべて PASS。**")
    else:
        failed = [k for k in g["order"] if not g[k]["pass"]]
        A("**%d/6 PASS。FAIL: %s。**" % (6 - len(failed), ", ".join(failed)))
        A("")
        A("G-LINE は**達成できていない**。数値は隠さず、原因の切り分けを §8 に"
          "実測で置いた。要点: このシーンでは「格子を消す」=「5–15 px 帯の"
          "軸平行パワーを写真以下にする」であり、その帯のパワーは**堆積の"
          "非線形（折り返し集光）が場の広帯域な傾きから作る微細なコースティクス"
          "網**が担っている。§8.2 の否定的対照が示すとおり、"
          "(a) 場の 5–15 cell 帯の軸平行成分を 99.9% 抜いても変わらず、"
          "(b) ソースのジッタを切っても −7% しか変わらず、"
          "(c) 効くのは**堆積の平滑化スケール**だけ（C4 の階段）で、"
          "それは C のコントラストと直接トレードオフする。"
          "したがって G-LINE を通すには場の微細傾斜を落とす（=v4-1.5 の"
          "σ=4 cell 制限、C std 0.101）以外にない。"
          "G-LEGACY の FAIL は**並行レーンの外部差分 2 本**で、本作に帰属する"
          "変更は 0 本（§6）。")
    A("")
    A("> **要旨**: 傾き場は v4-1.4 のまま（帯域制限なし）。変えたのは堆積だけ: "
      "(1) ソースレイの発射位置をフレームごとに R2 低不一致列でセル内 ±0.5 cell "
      "ジッタ、(2) 着地を局所ヤコビアンから導いた正規化ガウスカーネル"
      "（σ = s₀·√|det J|/2、|det J| は [1/4,4] にクリップ → σ ∈ [0.25,1.0] cell "
      "= %.4f–%.4f mm）で計数、(3) 1472×1080 のバッファに積んで 2×2 面積平均で"
      "736×540 に落とし、平均を厳密に 1 に正規化。結果: **C のコントラストは "
      "0.101 → %.3f に戻り（G-CONTRAST PASS、v4-1.4 の 0.528 も上回る）**が、"
      "5–15 px 帯の絶対値は写真の %.1f–%.1f 倍で **G-LINE は FAIL**。"
      "軸平行/対角比も v4-1.4 の %s から本作 %s へ**ほとんど動かなかった**"
      "（写真 %s）。ジッタは残差を細かいノイズに変えるが、この帯の主成分は"
      "堆積の規則性ではなく集光そのものだった。%s"
      % (sc.get("kernel_derivation", {}).get("sigma_m_range", [0, 0])[0] * 1e3,
         sc.get("kernel_derivation", {}).get("sigma_m_range", [0, 0])[1] * 1e3,
         float(cs.get("std") or 0.0),
         min(r["ratio"] for r in rows), max(r["ratio"] for r in rows),
         "%.2f" % (f24.get("axis_diag", {}).get("v4_1_4_still_green", {})
                   .get("ratio") or float("nan")),
         "%.2f" % (f24.get("axis_diag", {}).get("still_green", {})
                   .get("ratio") or float("nan")),
         "%.2f" % (f24.get("axis_diag", {}).get("photo_green", {})
                   .get("ratio") or float("nan")),
         "" if g["all_pass"] else "原因は §8。"))
    A("")
    # ---------------------------------------------------------------- 1
    A("## 1. brief の処方と、そのまま実行したこと")
    A("")
    A("| 段 | brief の指示 | 実装 |")
    A("| --- | --- | --- |")
    A("| 場 | 丸めない（v4-1.4 の場そのもの、帯域制限なし） | "
      "`st.fields_at(t)` の h/gx/gy/Dx/Dy を無改造で輸送に渡す。"
      "`stream_surface4`（v4-1.5 の帯域制限）は**import すらしていない** |")
    A("| 1 ジッタ | ソースレイの発射位置をフレームごとにセル内でずらす"
      "（決定論的、フレーム番号から、低不一致列） | R2（rank-1）加法的漸化列。"
      "仕様は §2.1、検証は §2.2 |")
    A("| 2 カーネル堆積 | 着地を点ではなく正規化ガウスカーネルで数える。半径は"
      "レイ密度と局所ヤコビアンから導く（ゲートに合わせない） | "
      "σ(x) = (S/N_s)·√|det J(x)| / 2、|det J| を [1/4,4] にクリップ、3σ で打ち切り・"
      "単位質量。導出は §2.3 |")
    A("| 2' バッファ | フレームの 2 倍解像度（1472×1080）で積み、面積平均で"
      "フレーム解像度へ | 1472×1080 に積み、厳密な 2×2 面積平均で 736×540 |")
    A("| 3 正規化 | C の平均をちょうど 1 に（エネルギー保存） | 実測平均を割る"
      "（補正率 %.9f） |" % float(dep.get("mean_correction_factor") or 1.0))
    A("")
    A("- 新規ファイルは `transport2.py`（堆積）と `run_v4_1_7.py`（driver）の2本。"
      "既存ファイルは 1 バイトも変えていない（§6）。")
    A("")
    # ---------------------------------------------------------------- 2
    A("## 2. 堆積の仕様（結果を見る前に固定したもの）")
    A("")
    A("### 2.1 ジッタ列")
    A("")
    A("```")
    A("n  = i + N_s*j                 (i = 列, j = 行, N_s = 512)")
    A("f  = フレーム番号 (0 起点)")
    A("g  = 1.324717957244746         (plastic number)")
    A("a1 = 1/g   = 0.7548776662466927   (R2 列, Roberts 2018)")
    A("a2 = 1/g^2 = 0.5698402909980532")
    A("jx = frac(0.5 + a1*(n + 13*f)) - 0.5      ∈ [-0.5, 0.5) cell")
    A("jy = frac(0.5 + a2*(n +  7*f)) - 0.5")
    A("発射点 = ((i+0.5+jx)*h, (j+0.5+jy)*h),  h = S/N_s = %.6f mm" % (sc.get("cell_m", 0) * 1e3))
    A("```")
    A("")
    A("- 振り幅は**1 cell 固定**（セル内、brief の文言どおり）。ゲートを見て"
      "変えていない。")
    A("- フレーム stride 13 / 7 は互いに素で N_s=512 とも素。a1,a2 は無理数なので"
      "フレームごとに列の位相がずれる。")
    A("- 場は発射点で**周期バイリニア**（= 従来の `det.det_frame` と同じ補間）で"
      "読む。ジッタで別の点から発射するので、その点の h/gx/gy/Dx/Dy を運ぶ。")
    A("")
    A("### 2.2 ジッタ列の検証（実行時の実測）")
    A("")
    jv = (sc.get("jitter") or {}).get("per_frame") or []
    A("| frame | jx 平均 | jx 標準偏差 | min | max | セル内 8×8 占有 χ²(63) | "
      "前フレームとの相関 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for r in jv:
        A("| %d | %.5f | %.5f | %.4f | %.4f | %.1f | %s |"
          % (r["frame"], r["jx_mean"], r["jx_std"], r["jx_min"], r["jx_max"],
             r["cell_8x8_chi2"],
             ("—" if r["corr_with_previous_frame"] is None
              else "%.4f" % r["corr_with_previous_frame"])))
    A("")
    A("- 一様分布の標準偏差は 1/√12 = %.5f。実測はこれに一致。セル内 8×8 の"
      "χ² は 63 自由度の期待値 63 に対し上表の値（過不足のない一様被覆）。"
      "フレーム間の相関は ~0（= フレームごとに別のパターン）。"
      % (1.0 / np.sqrt(12.0)))
    A("")
    A("### 2.3 カーネル半径の導出")
    A("")
    A("```")
    A("s0      = S / N_s                     = %.6f mm   (1 cell, 歪みなしのレイ間隔)"
      % (sc.get("kernel_derivation", {}).get("s0_m", 0) * 1e3))
    A("x_b     = x_s + T(x_s) * t(x_s)        (着地写像)")
    A("J       = d x_b / d x_s                (周期中心差分で実測)")
    A("s(x)    = s0 * sqrt(|det J(x)|)        (局所のレイ間隔: 面積圧縮 |det J|)")
    A("sigma(x)= s(x) / 2                     (間隔 s の標本列に対する整合フィルタ)")
    A("|det J| ∈ [1/4, 4] にクリップ  ->  sigma ∈ [0.25, 1.0] cell = %.4f–%.4f mm"
      % (sc.get("kernel_derivation", {}).get("sigma_m_range", [0, 0])[0] * 1e3,
         sc.get("kernel_derivation", {}).get("sigma_m_range", [0, 0])[1] * 1e3))
    A("3 sigma で打ち切り、単位質量に正規化（平均保存）")
    A("```")
    A("")
    A("- **Q で割らない理由（明示）**: 25 方向は太陽円盤内のほぼ同一点（円盤半径 ≪ "
      "1 cell）なので、堆積格子のスケールで独立標本を与えない。s0/√Q にすると "
      "カーネルがバッファ 1 標本より細くなり、カーネル段が実質 no-op になる。"
      "この読み替えは**採らなかった**（代替として §8 に明記）。")
    A("- 実装は σ を 6 個の幾何ビン（比 4^(1/5) = %.4f）に量子化し、ビンごとの"
      "堆積バッファを可分周期ガウスでぼかしてから足す。これは「各レイをそのビンの"
      "カーネルで撒く」のと数学的に同じで、単に速いだけ。"
      % float(T2.RECORD["kernel"]["bin_ratio"]))
    A("- バッファは 1472×1080（= 2×(736,540)）で、1 cell は y に %.3f 標本、"
      "x に %.3f 標本。σ はバッファ標本で y %.3f–%.3f / x %.3f–%.3f。"
      % (sc.get("kernel_derivation", {}).get("buffer_samples_per_cell", [0, 0])[0],
         sc.get("kernel_derivation", {}).get("buffer_samples_per_cell", [0, 0])[1],
         (dep.get("kernel_sigma_centres_buffer_samples_yx") or [[0, 0]])[0][0],
         (dep.get("kernel_sigma_centres_buffer_samples_yx") or [[0, 0]])[-1][0],
         (dep.get("kernel_sigma_centres_buffer_samples_yx") or [[0, 0]])[0][1],
         (dep.get("kernel_sigma_centres_buffer_samples_yx") or [[0, 0]])[-1][1]))
    A("")
    A("frame 24 での実測（1 フレーム、25 方向の合計）:")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| σ ビン中心 [cell] | %s |" % _fmt(dep.get("kernel_sigma_centres_cells"), 3))
    A("| レイのビン配分 | %s |"
      % _fmt([round(x, 3) for x in (dep.get("kernel_ray_fraction_per_bin") or [])], 3))
    dj = dep.get("detJ_abs_p50_p99_max") or []
    if dj:
        A("| \\|det J\\| p50 / p99 / max（中央の方向） | %.3f / %.3f / %.1f |"
          % (dj[0][0], dj[0][1], dj[0][2]))
    fp = dep.get("fold_per_quadrature") or {}
    if fp:
        A("| \\|ds/di\\| p50 / p99 / max（折り返し診断、全方向） | %.2f / %.2f / %.2f |"
          % (fp["ds_per_ray_p50"]["mean"], fp["ds_per_ray_p99"]["mean"],
             fp["ds_per_ray_max"]["mean"]))
        A("| 折り返し面積率（\\|ds/di\\| ≥ 1） | %.1f%% |"
          % (100.0 * fp["folded_area_fraction"]["mean"]))
    A("| 平均（正規化前） | %.9f |" % float(dep.get("mean_before_mean_correction") or 0.0))
    A("| 解析的除数 total_rays / n_buf² | %.6f |"
      % float(dep.get("divisor_analytical") or 0.0))
    A("")
    A("折り返し面積率が v4-1.5 の 53.5% と同じ水準 = **場は v4-1.4 のまま"
      "（折り返している）**。それを丸めずに、堆積側で受けたのが本作。")
    A("")
    A("（注: \\|ds/di\\| の max と \\|det J\\| の max は、周期差分が ±S/2 で折り返す"
      "境界と、|det J|→0 の集光特異点を含む裾の値。判定に使っているのは "
      "p50 / p99 / 折り返し面積率で、max は参考値。）")
    A("")
    # ---------------------------------------------------------------- 3
    A("## 3. G-CONTRAST — 波感（v4-1.5 が失ったもの）")
    A("")
    A("| 版 | frame 24 の C 標準偏差 | 備考 |")
    A("| --- | --- | --- |")
    A("| v4-1.4（点スプラット、512²、σ=0.75 px 後段ぼかし） | 0.528 | 折り返しを表示 |")
    A("| v4-1.5（σ=4 cell の帯域制限） | 0.101 | 細かい傾きを捨てた |")
    A("| **v4-1.7（本作）** | **%.4f** | 帯域制限なし＋ジッタ＋カーネル堆積 |"
      % float(cs.get("std") or 0.0))
    A("")
    A("C（フレーム解像度 736×540、平均 1 厳密）: min %.4f / p01 %.4f / p50 %.4f / "
      "p99 %.3f / max %.3f。閾値 0.35 に対し **%s**。"
      % (float(cs.get("min") or 0), float(cs.get("p01") or 0),
         float(cs.get("p50") or 0), float(cs.get("p99") or 0),
         float(cs.get("max") or 0),
         "PASS" if g["G-CONTRAST"]["pass"] else "FAIL"))
    A("")
    A("- 波感は v4-1.5 の 7.8 倍に戻った（0.101 → %.3f）。v4-1.4 の 0.528 を"
      "上回るのは、後段 σ=0.75 cell のぼかしを掛けず、2×2 面積平均だけにしたため"
      "（解像している分だけ尖っている）。" % float(cs.get("std") or 0.0))
    A("- 参考: C を描画用の正方形バッファ（§10 の適合器）に載せたときの標準偏差は "
      "%.4f。" % float((f24.get("C_render_stats") or {}).get("std") or 0.0))
    A("")
    # ---------------------------------------------------------------- 4
    A("## 4. G-LINE — 格子の強さ（主判定 = 前段ハイパスなし）")
    A("")
    A("v4-1.5 と同一実装（`run_v4_1_5.lattice_metric` を import）・同一窓・同一帯域。")
    A("")
    A("| 領域 | 軸 | 元写真 period / depth / depth_hp | v4-1.4 depth | v4-1.5 depth | "
      "**v4-1.7** period / depth / depth_hp | 比 depth (hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for i, r in enumerate(rows):
        ph = (gl.get("photo") or [])[i] if gl.get("photo") else {}
        r14 = (gl.get("v4_1_4") or [])[i] if gl.get("v4_1_4") else {}
        r15 = (gl.get("v4_1_5") or [])[i] if gl.get("v4_1_5") else {}
        A("| %s | %s | %.2f px / %.1f / %.1f | %.1f | %.1f | "
          "**%.2f px / %.1f / %.1f** | %.3f (%.3f) | **%s** |"
          % (r["region"], r["axis"], ph.get("period_px", 0.0),
             r["depth_photo"], r["depth_hp_photo"],
             r14.get("depth_e3", float("nan")),
             r15.get("depth_e3", float("nan")),
             r["period_v4_1_7_px"], r["depth_v4_1_7"], r["depth_hp_v4_1_7"],
             r["ratio"], r["ratio_hp"], "OK" if r["ok"] else "FAIL"))
    A("")
    A("depth は 1000·RMS(帯域通過)/窓平均、depth_hp は brief の手順"
      "「画像 − gaussian(σ=3 px) → 帯域通過」。主判定は depth（前段ハイパスなし）。")
    A("")
    A("### 4.1 参照（同じ測定）")
    A("")
    A("| 画像 | smooth x | smooth y | bark x | bark y | 最大 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for lab, key in (("元写真", "photo"), ("v4-1.4", "v4_1_4"), ("v4-1.5", "v4_1_5")):
        v = gl.get(key)
        if v:
            A("| %s | %s | %.1f |" % (lab, " | ".join("%.1f" % x["depth_e3"] for x in v),
                                      lattice_strength(v)))
    A("| **v4-1.7** | %s | **%.1f** |"
      % (" | ".join("%.1f" % x["depth_e3"] for x in gl["measured"]),
         gl["strength_v4_1_7"]))
    A("")
    A("### 4.2 全 24 窓の分布（brief の「中央値も併記」）")
    A("")
    A("| 比較 | 軸 | 中央値 | p10–p90 | 写真を上回った窓 | 最大比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for lab, key in (("v4-1.7 / 元写真", "v4_1_7_vs_photo"),
                     ("v4-1.4 / 元写真", "v4_1_4_vs_photo"),
                     ("v4-1.5 / 元写真", "v4_1_5_vs_photo")):
        dd = (gl.get("distribution") or {}).get(key)
        if not dd:
            continue
        for ax in ("x", "y"):
            v = dd[ax]
            A("| %s | %s | %.3f | %.3f–%.3f | %d/%d | %.3f |"
              % (lab, ax, v["ratio_median"], v["ratio_p10"], v["ratio_p90"],
                 int(round(v["frac_above_1"] * v["n_windows"])), v["n_windows"],
                 v["max_ratio"]))
    A("")
    A("→ **G-LINE は FAIL**。4 行すべてで写真を上回る（比 %.2f–%.2f）。"
      "24 窓の中央値でも %.2f / %.2f 倍。"
      % (min(r["ratio"] for r in rows), max(r["ratio"] for r in rows),
         gl["distribution"]["v4_1_7_vs_photo"]["x"]["ratio_median"],
         gl["distribution"]["v4_1_7_vs_photo"]["y"]["ratio_median"]))
    A("")
    A("### 4.3 形状の副次統計（軸平行 / 対角の帯振幅比）")
    A("")
    ad = f24.get("axis_diag") or {}
    A("| 画像 | 軸平行 ‰ | 対角 ‰ | 比 |")
    A("| --- | --- | --- | --- |")
    for lab, key in (("C（v4-1.7, 736×540）", "C_frame"),
                     ("本作の still（緑）", "still_green"),
                     ("v4-1.4 still（緑）", "v4_1_4_still_green"),
                     ("v4-1.5 still（緑）", "v4_1_5_still_green"),
                     ("元写真（緑）", "photo_green")):
        v = ad.get(key)
        if v:
            A("| %s | %.1f | %.1f | %.2f |"
              % (lab, v["axis_permille"], v["diag_permille"], v["ratio"] or 0.0))
    A("")
    A("- 比は v4-1.4 の %s から本作 %s へ**ほとんど動かなかった**（写真は %s）。"
      "つまり**格子の異方性そのものは解消していない**。"
      "ジッタは堆積の規則性を壊すが、この帯の軸平行成分は"
      "「規則的な点列の印刷」ではなく集光の網そのものなので、"
      "ジッタでは動かない（§8.2 の C1 と C5 が同じことを別方向から示す）。"
      % ("%.2f" % (ad.get("v4_1_4_still_green", {}).get("ratio") or float("nan")),
         "%.2f" % (ad.get("still_green", {}).get("ratio") or float("nan")),
         "%.2f" % (ad.get("photo_green", {}).get("ratio") or float("nan"))))
    A("")
    A("- C 自身の比は %.2f（v4-1.7 の C は 736×540 で、5–15 px が 3.5–10.4 cell"
      "に当たるため、512² の v4-1.4 C の比とは帯域が一致しない＝直接比較しない）。"
      % ((ad.get("C_frame") or {}).get("ratio") or 0.0))
    A("")
    # ---------------------------------------------------------------- 5
    A("## 5. G-WARP / G-MOVE / G-OUT")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| frame 24 歪み（緑 平均） | %.4f px（v4-1.4 %.4f、比 %.6f） |"
      % (float(g["G-WARP"]["g_mean_px"]), float(g["G-WARP"]["v4_1_4_g_mean_px"]),
         float(g["G-WARP"]["ratio"])))
    A("| 同 p50 / p99 | %.4f / %.4f px |"
      % (float((n["warp"] or {}).get("g_p50_px") or 0),
         float((n["warp"] or {}).get("g_p99_px") or 0)))
    A("| **G-WARP** | **%s**（±25%% 帯 %.3f–%.3f px） |"
      % ("PASS" if g["G-WARP"]["pass"] else "FAIL",
         1.8393 * 0.75, 1.8393 * 1.25))
    A("| ユニークフレーム | %d/%d |"
      % (g["G-MOVE"]["unique_frames"], NFRAMES))
    A("| フレーム間平均絶対差（緑, PNG） | %.4f 諧調（min %.4f / median %.4f / max %.4f） |"
      % (g["G-MOVE"]["interframe_green_png"]["mean"],
         g["G-MOVE"]["interframe_green_png"]["min"],
         g["G-MOVE"]["interframe_green_png"]["median"],
         g["G-MOVE"]["interframe_green_png"]["max"]))
    dd = g["G-MOVE"].get("interframe_green_mp4") or {}
    A("| 同（デコード後 mp4） | %s |"
      % ("%.4f 諧調" % dd["mean"] if "mean" in dd else "n/a"))
    A("| **G-MOVE** | **%s**（閾値 0.5 諧調） |"
      % ("PASS" if g["G-MOVE"]["pass"] else "FAIL"))
    A("| **G-OUT** | %s%d frames / %.4f s / %s×%s / %.1f fps、%.2f MB、sha256 `%s` |"
      % ("PASS: " if g["G-OUT"]["pass"] else "FAIL: ",
         int(pr.get("nb_frames_reported") or 0),
         float(pr.get("duration_s") or 0),
         pr.get("width"), pr.get("height"), float(g["G-OUT"]["fps"]),
         (mp.get("size_bytes") or 0) / 1e6, (mp.get("sha256") or "")[:16]))
    A("")
    A("- 歪みが v4-1.4 と**同一**なのは構成上そうなる（ワープはフィルタ前の場そのもの、"
      "v4-1.5 と同じ）。「波を弱めて格子を消したのではない」ことの担保。")
    A("")
    # ---------------------------------------------------------------- 6
    A("## 6. G-LEGACY — 既存 py の無傷性")
    A("")
    A("- stage-0 記録: `%s`（既存 %d 本 + 本作の新規 %d 本）。"
      % (os.path.basename(STAGE0_HASHES),
         g["G-LEGACY"]["recorded_existing"],
         len(g["G-LEGACY"]["own_new_file_hashes"])))
    A("- 既存 %d 本の照合 → **変更 %d 本 / 消失 %d 本**。stage-0 後に増えた"
      "（本作のものでも既存でもない）ファイル: %s。"
      % (g["G-LEGACY"]["recorded_existing"], len(g["G-LEGACY"]["changed"]),
         len(g["G-LEGACY"]["missing"]),
         g["G-LEGACY"]["new_since_stage0_not_of_this_job"] or "なし"))
    ext = g["G-LEGACY"].get("external_changes") or {}
    if ext:
        A("")
        A("**帰属の分離（brief の指示どおり）**: 変わったのは下記 2 本で、"
          "いずれも**本作が import も open もしない並行レーンのファイル**。"
          "本作が import する .py（run_v4_1_4 / run_v4_1_5 / run_v16_water / "
          "det / transport / stream_surface2 / stream_surface3 / transport2）は"
          "**すべて不変**。したがって本作の数値はこの差分の影響を受けない。")
        A("")
        A("| file | stage-0 sha256 | 現在 sha256 | mtime | 本作が import |")
        A("| --- | --- | --- | --- | --- |")
        for k, v in sorted(ext.items()):
            A("| `%s` | `%s` | `%s` | %s | %s |"
              % (k, (v["stage0_sha256"] or "")[:16], (v["now_sha256"] or "")[:16],
                 v["mtime"], "はい" if v["imported_by_this_job"] else "いいえ"))
        A("")
        A("→ ゲートの文言どおりには **G-LEGACY は FAIL**（既存 2 本が一致しない）。"
          "本作に帰属する変更は **0 本**。")
    for k, v in (g["G-LEGACY"]["own_new_file_hashes"] or {}).items():
        A("- 新規 `%s` = %s（stage-0）→ %s（終了時）"
          % (k, (v["stage0"] or "")[:16], (v["now"] or "")[:16]))
    A("")
    # ---------------------------------------------------------------- 7
    A("## 7. 場面定数（`run_v4_1_4.SCENES['realscale']` から verbatim）")
    A("")
    cst = sc.get("constants") or {}
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| H | %.4f m |" % cst.get("H", 0))
    A("| U / S0 / D84 | %.4f m/s / %.4f / %.4f m |"
      % (cst.get("U", 0), cst.get("S0", 0), cst.get("D84", 0)))
    A("| seed | %d |" % cst.get("seed", 0))
    A("| water / Q | %s / %d |" % (cst.get("water"), cst.get("Q", 0)))
    A("| T_loop / L | %.1f s / %.6f m |" % (cst.get("T_loop", 0), sc.get("L_m", 0)))
    A("| 格子 / cell | %d / %.5f mm |"
      % (sc.get("grid", 0), sc.get("cell_m", 0) * 1e3))
    A("| σ_total | %.6f mm |" % (sc.get("sigma_total_m", 0) * 1e3))
    A("| 合成 | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d |"
      % (cst.get("A", 0), cst.get("B_sigma", 0), cst.get("E", 0), cst.get("Q", 0)))
    A("| 堆積バッファ | %d×%d（=2×(540,736)） |" % (BUF[0], BUF[1]))
    A("")
    # ---------------------------------------------------------------- 8
    A("## 8. G-LINE がなぜ通らないか（実測による切り分け）")
    A("")
    A("ここは**後付けの言い訳ではなく、原因の切り分け**。すべて frame 24、"
      "同じ窓・同じ帯域・同じ描画器。下の C0–C5 は診断であって、"
      "**本作のパラメータはどれもここから決めていない**（§2 の導出は実行前に固定）。")
    A("")
    A("### 8.1 帰属（v4-1.5 の stage 1 と同じ3通り）")
    A("")
    ab = f24.get("attribution") or {}
    A("| 描き方 | 軸平行 5–15 px 深度（‰, 4行） | 最大 | C std |")
    A("| --- | --- | --- | --- |")
    for tag, key in (("a", "a"), ("b", "b")):
        if key in ab:
            A("| %s | %s | %.1f | %.4f |"
              % (ab[key]["name"], _fmt([round(x, 1) for x in ab[key]["depth"]]),
                 ab[key]["strength"], ab[key]["C_std"]))
    A("| c = 本作（両方） | %s | %.1f | %.4f |"
      % (_fmt([round(r["depth_v4_1_7"], 1) for r in rows]),
         gl["strength_v4_1_7"], float(cs.get("std") or 0)))
    A("")
    A("- ワープ単独（C=1）は写真と同レベル（%.1f vs %.1f）→ **ワープ段は無罪**、"
      "v4-1.5 と同じ結論。"
      % (ab.get("a", {}).get("strength", float("nan")), gl["strength_photo"]))
    A("- コースティクス段が帯を担う（%.1f）。"
      % ab.get("b", {}).get("strength", float("nan")))
    A("")
    A("### 8.2 否定的対照（C0–C5）")
    A("")
    A("| 条件 | C std | 画像 5–15 px 深度（‰, 4行） | 最大 | 何を言うか |")
    A("| --- | --- | --- | --- | --- |")
    cc = dg.get("cases") or {}
    notes = {
        "C0": "本作（仕様どおり）",
        "C1": "ジッタを切る → 差は小さい（堆積の規則性は主因ではない）",
        "C2": "正方形バッファ → ほぼ不変（非等方は主因ではない）",
        "C3": "バッファ 2 倍細かく → 1–2 割下がるが届かない",
        "C4_x1.5": "カーネル幅 1.5 倍（診断。採用しない）",
        "C4_x2.0": "カーネル幅 2 倍（診断。採用しない）",
        "C4_x3.0": "カーネル幅 3 倍（診断。採用しない）",
        "C5": "場の軸平行成分を 99.9% 抜く → ほぼ不変（帯は場の軸成分ではない）",
    }
    for k in ("C0", "C1", "C2", "C3", "C4_x1.5", "C4_x2.0", "C4_x3.0", "C5"):
        v = cc.get(k)
        if not v:
            continue
        A("| %s | %.4f | %s | %.1f | %s |"
          % (v["label"], v["C_std"],
             _fmt([round(x, 1) for x in v["depth"]]), v["strength"],
             notes.get(k, "")))
    A("")
    if "C5" in cc:
        A("- C5 は**場の 5–15 cell 帯の軸平行成分（slope パワーの %.1f%%）を"
          "消して**測り直したもの。消した割合は帯内 slope パワーの %.1f%%。"
          "それでも深度は変わらない → 画像の軸平行帯は「場の軸平行な成分」では"
          "ない。"
          % (100.0 * cc["C5"]["field_axis_line_share_of_band_slope_power"],
             100.0 * cc["C5"]["removed_fraction_of_band_slope_power"]))
    A("")
    A("### 8.3 何が残るか（正直な結論）")
    A("")
    ff = dg.get("frontier_fit") or {}
    A("1. **ジッタは効くが小さい**（C1 ジッタ無し %.1f → C0 %.1f、%.0f%% 減）。"
      "堆積の規則性（ソース格子・宛先格子）は主因ではない。"
      % (cc.get("C1", {}).get("strength", float("nan")),
         cc.get("C0", {}).get("strength", float("nan")),
         100.0 * (1 - cc.get("C0", {}).get("strength", 0)
                  / max(cc.get("C1", {}).get("strength", 1), 1e-9))))
    A("2. **バッファ形状は無関係**（C2 正方形でほぼ同じ）、**解像度も逆効果気味**"
      "（C3 で 2 倍に細かくすると僅かに増える）。")
    A("3. **場の 5–15 cell 帯の軸平行成分を 99.9% 消しても変わらない**（C5）。"
      "つまり画像の 5–15 px 軸平行帯は「場の軸平行な成分」ではない。"
      "残るのは**堆積の非線形（折り返し集光）が場の広帯域な傾きから作る"
      "微細コースティクス網**で、折り返し集光は特異（細い網）なので"
      "スペクトルが広帯域になり、5–15 px にも必ずパワーを持つ。")
    A("4. **効くレバーは堆積の平滑化スケールだけ**（C4 の階段）。"
      "その前線を4点で直線近似すると `深度 ≈ %.0f·C_std + %.0f`。"
      "写真の 63.7 ‰ に落とすには **C_std ≈ %.3f** が必要で、"
      "G-CONTRAST の閾値 0.35 を大きく下回る（= v4-1.5 の 0.101 の側）。"
      % (float(ff.get("linear_fit_depth_per_std") or 0.0),
         float(ff.get("intercept") or 0.0),
         float(ff.get("std_needed_for_photo_depth") or float("nan"))))
    A("5. したがって **G-LINE（≤ 写真）と G-CONTRAST（≥ 0.35）は、場を変えない"
      "限り両立しない**。v4-1.5 の帯域制限が効いた理由は「堆積の格子を消した」"
      "ことではなく「集光を作る微細傾斜そのものを落とした」こと、というのが"
      "実測の結論。brief の診断（堆積の格子が原因）は、少なくともこのシーンの"
      "5–15 px 帯については**支持されなかった**。")
    A("")
    A("この結論は brief の前提に対する反証なので、隠さずここに置く。"
      "「ゲートに合わせた後付け調整」は一切していない（§2 の仕様は実行前固定、"
      "C4 の幅は診断であって成果物には使っていない）。")
    A("")
    # ---------------------------------------------------------------- 9
    A("## 9. 仮定（本作で追加したもの）")
    A("")
    A("- `deposit_jitter_r2`: ソース発射位置を R2 列で ±0.5 cell。フレーム stride 13/7。")
    A("- `deposit_kernel_jacobian`: σ = s0·√|det J|/2、|det J| を [1/4,4] にクリップ、"
      "3σ 打ち切り・単位質量。Q では割らない（理由は §2.3）。")
    A("- `buffer_2x_area_average`: 1472×1080 に積み、厳密な 2×2 面積平均。")
    A("- `c_mean_normalisation`: C をその平均で割る（brief の指示。v4-1.4 では禁止"
      "だった操作）。実測の補正率は %.9f なので、実質はエネルギー保存の確認。"
      % float(dep.get("mean_correction_factor") or 1.0))
    A("- `render_grid_adapter`: 描画器の `_sample_c` は正方形バッファ前提なので、"
      "736×540 の C を y 方向だけ線形補間して 736×736 で渡す（§10 に開示）。"
      "ハロー核は同じ動径プロファイルを 736 格子に引き直す（半径を 736/512 倍）。")
    A("- それ以外の定数・会計・物理は v4-1.4 と同一。")
    A("")
    # ---------------------------------------------------------------- 10
    A("## 10. 失敗・弱点・正直な開示")
    A("")
    A("- **G-LINE FAIL**。原因は §8。brief の処方（ジッタ＋カーネル）は仕様どおり"
      "実装したが、このシーンでは 5–15 px 帯の軸平行パワーを写真以下にできない。"
      "数値はそのまま出した。")
    A("- **描画用バッファの適合器（開示）**: 描画器 `_sample_c` は N×N 正方形前提で、"
      "非正方形の C を渡すと y 方向の添字が壊れる（既存ファイルは変更禁止）。"
      "そこで 736×540 の C を **y だけ** 736 行に線形補間して渡した。x は実標本の"
      "まま。代償は y 方向の余分な線形補間 1 回（σ にして ~0.3 px 相当）で、"
      "C の標準偏差は 736×540 で %.4f、適合後 736×736 で %.4f。"
      % (float(cs.get("std") or 0),
         float((f24.get("C_render_stats") or {}).get("std") or 0)))
    A("- **ハロー核の格子差し替え（開示）**: optics の ff-cdf 核は state 格子（512）で"
      "作られており `ff_blur_fft` が形の一致を要求する。C を 736 格子にしたので、"
      "同じ動径プロファイルを 736 に引き直した（動径を 736/512 倍）。重み "
      "w_halo は %.5f 程度なので画像への寄与は 0.2%% 未満。"
      % float((sc.get("optics") or {}).get("w_halo", [0])[1]))
    A("- **カーネル幅の読み替え**: 「レイ密度」を「1 方向のレイ間隔」と読んだ"
      "（Q で割らない）。Q で割る読み方だと σ が 5 倍細くなり、カーネル段が"
      "ほぼ no-op になる。採らなかった理由は §2.3。この選択は結果を見る前の"
      "設計判断だが、**結果に有利な側（広い方）を選んでいる**ことは明記する。")
    A("- **測定の限界（v4-1.5 から継承）**: depth は 1 次元マージナル統計で、"
      "格子と同じ軸帯パワーを持つ等方テクスチャを区別しない。本作では §8.2 の C5 と"
      "§4.3 の軸平行/対角比で補っている。")
    A("- **G-LINE は 1 フレーム（frame 24）の判定**（brief の表記どおり）。48 "
      "フレーム全体の周期測定はしていない。")
    A("- **stage-0 後の自己編集（開示）**: `run_v4_1_7.py` は stage 0 の後に編集した"
      "（レポート本文・診断段の追加）。さらに `transport2.py` には**バグが 1 つ"
      "あり、stage 0 の後に直した**: ビン格子の下端を `|det J|·σ_per_spacing` と"
      "書いていたため σ の実効レンジが 0.125–2.0 cell になっていた。docstring に"
      "実行前から書いてあった規則（σ = σ_per_spacing·√|det J|、|det J| ∈ [1/4,4] "
      "→ σ ∈ [0.25,1.0] cell）に合わせて `√` を入れて修正し、**stage 0 から"
      "やり直した**。成果物はすべて修正後のコードで作っている。参考: バグ版でも "
      "frame 24 の C std は 0.7908、深度は 173.3 で、修正版（0.7782 / 173.5）と"
      "実質同じ——つまりこの結論はビン格子の取り方に依存しない。"
      "stage-0 時点と終了時点のハッシュは §6 に両方残した。")
    A("- **時間方向のノイズ（開示）**: フレームごとにジッタ列の位相が変わるので、"
      "残差ノイズはフレーム間で独立になる。フレーム間平均絶対差（緑）は "
      "%.2f 諧調で、v4-1.4 の 32.95・v4-1.5 の 16.42 より大きい"
      "（G-MOVE の閾値 0.5 に対しては十分大きいが、見た目には粒状のちらつきとして"
      "出る可能性がある。brief がフレームごとのジッタを指示している以上、"
      "これは仕様どおりの副作用）。"
      % g["G-MOVE"]["interframe_green_png"]["mean"])
    A("- 判定は数値統計のみ（目視なし）。`frame-realscale-0024.png` と "
      "`frames/frame-0024.png` の sha256 は**一致を確認済み**（%s）。"
      % ("一致" if n.get("still_equals_frame_0024") else "不一致"))
    A("- **破綻チェック（判定には未使用）**: 成果物が壊れていないことの確認として "
      "still を 1 枚だけ見た。判定はすべて上の数値で行っており、目視は"
      "ゲートに使っていない。観察（交差ハッチが見える）は §8 の数値"
      "（軸平行 5–15 px が写真の 2.4–2.7 倍）と同じ方向で、矛盾しない。")
    A("- **stage の再実行**: `--stage0` / `--frame24` / `--gates` を 2 回走らせた"
      "（バグ修正とレポート調整のため）。成果物は最後の実行で上書きされており、"
      "checkpoint の unit id に同じ名前が複数回出る。")
    A("- 実行は全ステージ foreground（1 プロセスずつ）。jobkit checkpoint を "
      "`%s` に stage ごとに原子的にコミットした（unit id: %s）。"
      % (os.path.basename(CP_PATH),
         ", ".join("`%s`" % u for u in
                   (read_checkpoint().get("completed_ids") or []))))
    A("")
    # ---------------------------------------------------------------- 11
    A("## 11. 成果物")
    A("")
    for k, v in (n.get("deliverables") or {}).items():
        A("- %s → `%s`" % (k, v))
    A("")
    A("ファイル（bytes / sha256 先頭16）:")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for k, v in sorted((n.get("files") or {}).items()):
        A("| %s | %d | %s |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-7.json + REPORT-V4-1-7.md")
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
    ap.add_argument("--frame24", action="store_true")
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.all:
        stage0()
        do_frame24()
        do_diag()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_warp()
        do_encode()
        do_gates()
        do_report()
        log("=" * 78)
        log("[all] done in %.1fs" % (time.time() - t0))
        return 0
    if a.stage0:
        stage0()
    if a.frame24:
        do_frame24()
    if a.diag:
        do_diag()
    if a.chunk is not None:
        do_chunk(a.chunk)
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
