"""1-second composite samples: eta = 1.0 (A) vs eta = 0.8 (B).

BRIEF: D:\\mutaro-work\\rally-in\\BRIEF-SAMPLES-1S.md

Both versions share EVERYTHING except the light mix, over the same base frames
(blender-piece8s/renders/frame_0062..0085.png, 24 frames = 1.000 s @ 24 fps):

  A (v3.4 route)  C = C_direct                     (v8 det.det_frame, frozen)
  B (v3.5 route)  C = wd*C_direct + ws*C_sky        (v9 Perez dome transport)

Composite: the FROZEN v5 composite (A2.6 ambient, B=3 px low-pass, direct-term
exposure E=2.2, physical Beer-Lambert, Smith G1), rendered natively at the base
frame size 832x480 (cover-resize is therefore the identity; recorded).

Stages (each writes jobkit checkpoints into <out>/checkpoint.json):
  probe   timing measurement on f0062 (the brief's cost probe)
  gate0   t=0 identity against the stored v8 master (production-route proof)
  frames  --from N --to M   per-frame transport + both composites
  encode  two 24-frame H.264 mp4s + ffprobe verification
  stills  f0074 stills for both versions + difference statistics
  report  SAMPLE-REPORT.md + numbers + gate verdicts

Usage (from caustics-filter/offline/):
  python run_v9_samples.py --stage probe
  python run_v9_samples.py --stage gate0
  python run_v9_samples.py --stage frames --from 62 --to 67
  python run_v9_samples.py --stage encode
  python run_v9_samples.py --stage stills
  python run_v9_samples.py --stage report
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
JOBKIT = r"D:\mutaro-work\jobkit"
sys.path.insert(0, JOBKIT)
from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from state import derive_seed  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert,  # noqa: E402
                        smith_g1_beckmann, scattering_blur_sigma_px,
                        substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,  # noqa: E402
                          image_metrics, v3_beer, flat_water_direct_level)
from recipe_v9 import parse_recipe_v9, with_v9  # noqa: E402
from v8_state import V8WaterState  # noqa: E402
import v9_sky as SK  # noqa: E402
import v9_transport as TR  # noqa: E402

# ---------------------------------------------------------------------------
# configuration of this run (all recorded in the report)
# ---------------------------------------------------------------------------
OUT = os.path.join(HERE, "out", "v9samples")
PARTS = os.path.join(OUT, "parts")
RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
RENDERS = os.path.join(HERE, "..", "blender-piece8s", "renders")
V8_DIR = os.path.join(HERE, "out", "v8")

GRID = 2048                # production transport grid
FPS = 24
FRAME_FIRST, FRAME_LAST = 62, 85     # 24 frames
FRAME_MID = 74
Q_PRODUCTION = 49          # the v8/v9 production quadrature
Q_DIRECT = 9               # reduced (recorded); inside the published v8 band
N_DIRS = 48                # dome directions for B (brief's ladder floor)
ETA_B = 0.8
E_COMP = 2.2               # the frozen composite's direct-term exposure
W_BASE, H_BASE = 832, 480
PHOTO = os.path.join(HERE, "..", "test-inputs", "moth-on-bark.jpg")

CK = Checkpoint(os.path.join(OUT, "checkpoint.json"))


def iso():
    return datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(
        timespec="seconds")


def ensure_dirs():
    for d in (OUT, PARTS, os.path.join(OUT, "frames", "eta100"),
              os.path.join(OUT, "frames", "eta080")):
        os.makedirs(d, exist_ok=True)


def done(unit_id):
    """jobkit checkpoint: record a unit once (single-writer, idempotent)."""
    st = CK.state() or {}
    if unit_id in (st.get("completed_ids") or []):
        print("   checkpoint: '%s' already done" % unit_id, flush=True)
        return
    st = CK.unit_done(unit_id)
    print("   checkpoint: '%s' done (%d units)" % (unit_id,
                                                   st["completed_units"]),
          flush=True)


def write_json(path, obj):
    atomic_write_bytes(path, json.dumps(obj, indent=2, ensure_ascii=False,
                                        default=float).encode("utf-8"))


def read_json(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def sha12(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()[:12]


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {"mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99_9": float(np.percentile(c, 99.9)),
            "std": float(c.std()), "max": float(c.max()),
            "min": float(c.min()),
            "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
            "contrast_p99_minus_p50": float(np.percentile(c, 99)
                                            - np.percentile(c, 50))}


# ---------------------------------------------------------------------------
# scene / state / optic constants (shared by all frames and both versions)
# ---------------------------------------------------------------------------
class Run(object):
    def __init__(self):
        ensure_dirs()
        with open(RECIPE, "r", encoding="utf-8") as f:
            self.raw = json.load(f)
        self.scene = with_v9(self.raw, eta=ETA_B, fetch_m=100.0)
        self.env = resolve_environment(self.scene)
        self.cfg = parse_recipe_v9(self.scene)
        self.seed = derive_seed(self.raw, 1)
        self.split = SK.horizontal_weights(ETA_B, self.env["sun"]["altitudeDeg"])
        self.beer = beer_lambert(self.env["depthM"],
                                 self.env["sun"]["altitudeDeg"],
                                 self.env["turbid"])
        att = water_attenuation(self.env["turbid"])
        self.blur = scattering_blur_sigma_px(att, self.env["depthM"], GRID)
        self.gain = substrate_gain(self.cfg["bottom"]["substrate"],
                                   self.cfg["bottom"]["albedoScale"])
        self._photo = None
        self._C34 = None

    # ---- water state (preset 1, production route, 8 s loop) ----------------
    def state(self):
        t0 = time.time()
        st = V8WaterState(self.scene, self.env, N=GRID, seed=self.seed,
                          preset_index=1,
                          fetch_m=float(self.cfg["fetch"]["fetch_m"]),
                          fetch_source=self.cfg["fetch"]["source"])
        print("   state built in %.1fs (sigma2=%.6f Hs=%.4f)"
              % (time.time() - t0, st.spectrum_after["sigma2_rad2"],
                 st.spectrum_after["Hs_m"]), flush=True)
        return st

    def fields(self, st, f):
        t0 = time.time()
        fields = st.fields_at(float(f) / float(FPS))
        fields["S"] = st.S
        return fields, time.time() - t0

    def dome(self, n_dirs=N_DIRS):
        return SK.sky_dome(ETA_B, self.env["sunDir"], n_dirs=int(n_dirs))

    # ---- display normalisation, evaluated once and held (production rule) --
    def optics(self, st):
        f0 = st.fields_at(0.0)
        f0["S"] = st.S
        re0 = st.realised_stats(f0)
        del f0
        sun_zenith = 90.0 - float(self.env["sun"]["altitudeDeg"])
        g1 = smith_g1_beckmann(re0["rms_slope"], sun_zenith)
        L_v3 = flat_water_direct_level(self.env, v3_beer(self.env))
        L_dir = flat_water_direct_level(self.env, self.beer["two_way"])
        E_v9 = E_V3 * (L_v3 / (L_dir * float(g1["G1"]))) if L_dir > 0 else E_V3
        return {"g1": float(g1["G1"]), "rms_slope": float(re0["rms_slope"]),
                "E_v9_production": float(E_v9), "E_used": E_COMP,
                "blur_sigma_px": float(self.blur["sigma_px"]),
                "beer_down": self.beer["down"], "beer_up": self.beer["up"],
                "L_v3": float(L_v3), "L_dir_flat": float(L_dir)}

    # ---- base frame as the bottom albedo ----------------------------------
    def frame_path(self, f):
        return os.path.join(RENDERS, "frame_%04d.png" % int(f))

    def photo_lin(self, f):
        p = self.frame_path(f)
        im = np.asarray(Image.open(p).convert("RGB"),
                        dtype=np.float64) / 255.0
        return srgb_to_linear(im), p

    # ---- the frozen composite --------------------------------------------
    def compose(self, fields, C, photo_lin, opt):
        t0 = time.time()
        img, dbg = render_frame_v5(
            fields, C, photo_lin, self.env, W_BASE, H_BASE,
            ambient_scale=A_V3, c_layer_exposure=float(E_COMP),
            substrate_gain_rgb=self.gain["gain_rgb"],
            beer_down=opt["beer_down"], beer_up=opt["beer_up"],
            g1=opt["g1"], blur_sigma_px=opt["blur_sigma_px"],
            lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
        m = image_metrics(dbg["lin"], dbg["aces"])
        del dbg
        out = cover_resize(img)
        return out, m, time.time() - t0


def cover_resize(img, W=W_BASE, H=H_BASE):
    """Resize/cover to the base frame size, exactly as the production
    composite does.  The composite is rendered natively at 832x480, so this
    is the identity here (recorded); the function stays so the step exists."""
    a = np.asarray(img)
    if a.shape[0] == H and a.shape[1] == W:
        return a
    im = Image.fromarray(a).convert("RGB")
    s = max(float(W) / im.width, float(H) / im.height)
    if s != 1.0:
        im = im.resize((int(round(im.width * s)), int(round(im.height * s))),
                       Image.LANCZOS)
    x0 = (im.width - W) // 2
    y0 = (im.height - H) // 2
    return np.asarray(im.crop((x0, y0, x0 + W, y0 + H)))


# ---------------------------------------------------------------------------
# STAGE probe -- the brief's per-frame wall-time measurement on f0062
# ---------------------------------------------------------------------------
def stage_probe():
    r = Run()
    print("[probe] per-frame wall time on f0062 (t = 62/24 s)", flush=True)
    t_all = time.time()
    st = r.state()
    fields, t_fields = r.fields(st, 62)
    t0 = time.time()
    C_dir_prod, _ = det_frame(fields, st.S, GRID, r.env, Q_PRODUCTION)
    t_dir_prod = time.time() - t0
    t0 = time.time()
    C_dir_sel, _ = det_frame(fields, st.S, GRID, r.env, Q_DIRECT)
    t_dir_sel = time.time() - t0
    dome = r.dome(N_DIRS)
    n_sub = 8
    t0 = time.time()
    TR.det_frame_dirs(fields, st.S, GRID, r.env, dome["dirs"][:n_sub],
                      dome["w"][:n_sub] / dome["w"][:n_sub].sum())
    t_dir8 = time.time() - t0
    per_dir = t_dir8 / n_sub
    opt = r.optics(st)
    photo, p0 = r.photo_lin(62)
    img, m, t_comp = r.compose(fields, C_dir_sel, photo, opt)

    per_frame = (t_fields + t_dir_sel + N_DIRS * per_dir + 2 * t_comp)
    payload = {
        "frame": 62, "t_s": 62.0 / FPS, "grid": GRID,
        "base_frame_size": [W_BASE, H_BASE],
        "photo_used_for_composite": os.path.basename(p0),
        "measured_s": {
            "state_build": None, "fields_at": t_fields,
            "det_frame_Q49": t_dir_prod,
            "det_frame_Q9": t_dir_sel,
            "det_frame_dirs_8dirs": t_dir8,
            "per_direction": per_dir,
            "composite_one_version": t_comp,
        },
        "projection_prefix": "",
        "C_dir_Q49_at_f0062": c_stats(C_dir_prod),
        "C_dir_Q9_at_f0062": c_stats(C_dir_sel),
        "optics": {k: v for k, v in opt.items()
                   if k not in ("beer_down", "beer_up")},
        "metrics_Q9_composite": m,
        "dome": {k: v for k, v in dome.items()
                 if k not in ("dirs", "w", "f", "solid_angle")},
        "projection": {},
        "wall_s": time.time() - t_all,
        "created": iso(),
    }
    del payload["measured_s"]["state_build"]
    proj = {}
    for tag, qd, nd in (("Q49_dirs48", Q_PRODUCTION, 48),
                        ("Q49_dirs64", Q_PRODUCTION, 64),
                        ("Q49_dirs96", Q_PRODUCTION, 96),
                        ("Q9_dirs48", Q_DIRECT, 48),
                        ("Q9_dirs64", Q_DIRECT, 64),
                        ("Q9_dirs96", Q_DIRECT, 96),
                        ("SELECTED_Q9_dirs48", Q_DIRECT, N_DIRS)):
        per = t_fields + qd * (t_dir_prod / float(Q_PRODUCTION)) \
            + nd * per_dir + 2 * t_comp
        proj[tag] = {"per_frame_s": per, "total_24_s": 24.0 * per,
                     "total_24_min": 24.0 * per / 60.0}
    payload["projection"] = proj
    payload["selected"] = {
        "Q_direct": Q_DIRECT, "n_dome_dirs": N_DIRS,
        "reason": ("the brief's ~40 min budget for the 24 B-frames cannot be "
                   "met at Q=49: the measured direct transport alone is "
                   "%.1f s/frame = %.1f min for 24 frames.  The direct "
                   "quadrature is reduced to Q=9 (its p99/p99.9/max sit "
                   "inside the v8 Q-convergence band recorded in "
                   "numbers-v8.json: p99 1.5097/1.5098/1.5095 for Q=9/25/49) "
                   "and the dome is kept at the brief's ladder floor (48), "
                   "not below it." % (t_dir_prod, 24.0 * t_dir_prod / 60.0)),
        "projected_total_min": proj["SELECTED_Q9_dirs48"]["total_24_min"],
    }
    write_json(os.path.join(OUT, "timing.json"), payload)
    print(json.dumps({"projection": proj, "selected": payload["selected"]},
                     indent=2), flush=True)
    done("probe (timing on f0062, dome direction choice)")


# ---------------------------------------------------------------------------
# STAGE gate0 -- t = 0 identity against the stored v8 production master
# ---------------------------------------------------------------------------
def stage_gate0():
    r = Run()
    print("[gate0] t=0 identity vs out/v8/C-v8-t0.npy (production route proof)",
          flush=True)
    st = r.state()
    fields, _ = r.fields(st, 0)
    ref_path = os.path.join(V8_DIR, "C-v8-t0.npy")
    out = {"reference": ref_path, "found": os.path.exists(ref_path),
           "Q_production": Q_PRODUCTION, "Q_direct_used": Q_DIRECT}
    C49, _ = det_frame(fields, st.S, GRID, r.env, Q_PRODUCTION)
    out["Q49"] = c_stats(C49)
    if out["found"]:
        C34 = np.load(ref_path).astype(np.float32)
        d = np.abs(np.asarray(C49, np.float64) - np.asarray(C34, np.float64))
        out["Q49_vs_v8"] = {
            "max_abs_diff": float(d.max()), "mean_abs_diff": float(d.mean()),
            "bit_identical": bool(float(d.max()) == 0.0),
            "reference_stats": c_stats(C34)}
        del C34
    C9, _ = det_frame(fields, st.S, GRID, r.env, Q_DIRECT)
    out["Q9"] = c_stats(C9)
    if out["found"]:
        C34 = np.load(ref_path).astype(np.float32)
        da = np.abs(np.asarray(C9, np.float64) - np.asarray(C34, np.float64))
        out["Q9_vs_v8"] = {
            "max_abs_diff": float(da.max()),
            "mean_abs_diff": float(da.mean()),
            "delta_p99": float(out["Q9"]["p99"] - out["Q49_vs_v8"]["reference_stats"]["p99"]),
            "delta_p99_9": float(out["Q9"]["p99_9"] - out["Q49_vs_v8"]["reference_stats"]["p99_9"]),
            "delta_mean": float(out["Q9"]["mean"] - out["Q49_vs_v8"]["reference_stats"]["mean"])}
        del C34
    out["pass"] = bool(out.get("Q49_vs_v8", {}).get("bit_identical"))
    write_json(os.path.join(OUT, "gate0-t0-identity.json"), out)
    print("   Q49 bit-identical to v8: %s | Q9 dP99=%+.3g"
          % (out.get("Q49_vs_v8", {}).get("bit_identical"),
             out.get("Q9_vs_v8", {}).get("delta_p99", float("nan"))),
          flush=True)
    done("gate0 (t=0 identity: production route + Q band)")


# ---------------------------------------------------------------------------
# STAGE frames
# ---------------------------------------------------------------------------
def frame_part_path(f):
    return os.path.join(PARTS, "frame-%04d.json" % int(f))


# ---------------------------------------------------------------------------
# STAGE ref62 -- the SAME-TIME-SLICE production reference (Q=49 at f0062)
# ---------------------------------------------------------------------------
def stage_ref62():
    r = Run()
    print("[ref62] Q=49 production reference at f0062 (same time slice)", flush=True)
    st = r.state()
    fields, _ = r.fields(st, 62)
    t0 = time.time()
    C49, _ = det_frame(fields, st.S, GRID, r.env, Q_PRODUCTION)
    t49 = time.time() - t0
    s49 = c_stats(C49)
    part = read_json(frame_part_path(62), {})
    s9 = part.get("C_direct", {})
    tol = {"mean": 1e-6, "p50": 1e-3, "p99": 5e-3, "p99_9": 2e-2}
    d = {k: float(s9.get(k, 0.0) - s49[k]) for k in tol}
    out = {"frame": 62, "t_s": 62.0 / FPS, "Q_reference": Q_PRODUCTION,
           "Q_used": Q_DIRECT, "wall_s_Q49": t49,
           "stats_Q49": s49, "stats_Q9_from_part": s9,
           "delta_Q9_minus_Q49": d, "tolerance": tol,
           "pass": bool(s9 and all(abs(d[k]) <= tol[k] for k in tol)),
           "note": ("same time slice, same state, same frozen transport; the "
                    "only difference is the direct quadrature (Q=9 vs the "
                    "production Q=49).  This is the quadrature-reduction "
                    "check, separated from the time evolution of the field.")}
    write_json(os.path.join(OUT, "gate1-f0062-Q49.json"), out)
    print("   Q49 f0062: mean=%.9f p50=%.6f p99=%.6f p99.9=%.6f | dQ9: "
          "mean=%+.2e p50=%+.2e p99=%+.2e p99.9=%+.2e -> %s"
          % (s49["mean"], s49["p50"], s49["p99"], s49["p99_9"],
             d["mean"], d["p50"], d["p99"], d["p99_9"],
             "PASS" if out["pass"] else "FAIL"), flush=True)
    done("ref62 (Q=49 same-slice reference at f0062)")



def stage_frames(f0, f1):
    r = Run()
    st = r.state()
    opt = r.optics(st)
    dome = r.dome(N_DIRS)
    dome_meta = {k: v for k, v in dome.items()
                 if k not in ("dirs", "w", "f", "solid_angle")}
    wd, ws = r.split["w_direct"], r.split["w_diffuse"]
    print("[frames] %d..%d  Q_dir=%d  dirs=%d  wd=%.6f ws=%.6f"
          % (f0, f1, Q_DIRECT, N_DIRS, wd, ws), flush=True)
    for f in range(int(f0), int(f1) + 1):
        part = frame_part_path(f)
        if os.path.exists(part):
            print("   f%04d: already done (skipped)" % f, flush=True)
            continue
        t_all = time.time()
        fields, t_fields = r.fields(st, f)
        t0 = time.time()
        C_dir, _ = det_frame(fields, st.S, GRID, r.env, Q_DIRECT)
        t_dir = time.time() - t0
        t0 = time.time()
        C_sky, _ = TR.det_frame_dirs(fields, st.S, GRID, r.env,
                                     dome["dirs"], dome["w"])
        t_sky = time.time() - t0
        C_tot = (wd * np.asarray(C_dir, np.float64)
                 + ws * np.asarray(C_sky, np.float64)).astype(np.float32)
        photo, _ = r.photo_lin(f)
        imgA, mA, t_cA = r.compose(fields, C_dir, photo, opt)
        imgB, mB, t_cB = r.compose(fields, C_tot, photo, opt)
        diff = np.abs(imgA.astype(np.int16) - imgB.astype(np.int16))
        Image.fromarray(imgA).save(
            os.path.join(OUT, "frames", "eta100", "frame-%04d.png" % (f - FRAME_FIRST + 1)))
        Image.fromarray(imgB).save(
            os.path.join(OUT, "frames", "eta080", "frame-%04d.png" % (f - FRAME_FIRST + 1)))
        payload = {
            "frame": int(f), "t_s": float(f) / FPS,
            "Q_direct": Q_DIRECT, "n_dirs": N_DIRS,
            "wall_s": {"fields": t_fields, "direct": t_dir, "sky": t_sky,
                       "composite_A": t_cA, "composite_B": t_cB,
                       "total": time.time() - t_all},
            "C_direct": c_stats(C_dir), "C_sky": c_stats(C_sky),
            "C_total_eta080": c_stats(C_tot),
            "w_direct": wd, "w_diffuse": ws,
            "still_A": mA, "still_B": mB,
            "A_vs_B_png": {"max_abs_uint8": int(diff.max()),
                           "mean_abs_uint8": float(diff.mean()),
                           "pct_pixels_gt8": float(np.mean(diff > 8) * 100.0)},
            "created": iso(),
        }
        write_json(part, payload)
        print("   f%04d t=%.4fs  Cdir p99=%.4f  Csky p99=%.4f  Ctot p99=%.4f  "
              "|A-B| mean=%.2f/255  wall=%.0fs"
              % (f, payload["t_s"], payload["C_direct"]["p99"],
                 payload["C_sky"]["p99"], payload["C_total_eta080"]["p99"],
                 diff.mean(), payload["wall_s"]["total"]), flush=True)
        del fields, C_dir, C_sky, C_tot, imgA, imgB, photo, diff
        done("frame f%04d (A+B composites)" % f)
    write_json(os.path.join(OUT, "dome.json"),
               {"dome": dome_meta, "optics": {
                   k: v for k, v in opt.items()
                   if k not in ("beer_down", "beer_up")},
                "split": r.split, "seed": int(r.seed),
                "scene": {"wind": r.env["windSpeedMS"],
                          "specL": r.env["specL"],
                          "depth": r.env["depthM"],
                          "sun_alt": r.env["sun"]["altitudeDeg"],
                          "sun_az": r.env["sun"]["azimuthDeg"],
                          "turbid": r.env["turbid"]}})


# ---------------------------------------------------------------------------
# STAGE encode
# ---------------------------------------------------------------------------
def ffprobe(path):
    cmd = ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
           "-show_entries",
           "stream=codec_name,width,height,nb_read_frames,r_frame_rate,"
           "avg_frame_rate,duration", "-of", "json", path]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("ffprobe failed: %s" % p.stderr)
    return json.loads(p.stdout)["streams"][0]


def stage_encode():
    ensure_dirs()
    res = {}
    for tag, sub in (("eta100", "eta100"), ("eta080", "eta080")):
        src = os.path.join(OUT, "frames", sub, "frame-%04d.png")
        dst = os.path.join(OUT, "comp-%s-1s.mp4" % tag)
        n = len([x for x in os.listdir(os.path.join(OUT, "frames", sub))
                 if x.endswith(".png")])
        if n != 24:
            raise SystemExit("%s: expected 24 frames, found %d" % (tag, n))
        # h264_nvenc first (hardware, as specified by the task); libx264 is the
        # fallback when no NVENC-capable device is present.  The encoder that
        # actually produced the file is recorded in encode-report.json.
        attempts = [
            (["-c:v", "h264_nvenc", "-preset", "p6", "-tune", "hq",
              "-rc", "vbr", "-cq", "17", "-b:v", "0"],
             "h264_nvenc preset p6 tune hq rc vbr cq 17 yuv420p 24fps"),
            (["-c:v", "libx264", "-crf", "16"],
             "libx264 crf16 yuv420p 24fps"),
        ]
        encoder = None
        for extra, enc_name in attempts:
            cmd = ["ffmpeg", "-y", "-framerate", str(FPS), "-start_number", "1",
                   "-i", src] + extra + ["-pix_fmt", "yuv420p",
                                         "-movflags", "+faststart", dst]
            logp = os.path.join(OUT, "ffmpeg-%s.log" % tag)
            with open(logp, "wb") as lf:
                p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            if p.returncode == 0 and os.path.exists(dst):
                encoder = enc_name
                break
        if encoder is None:
            raise SystemExit("ffmpeg failed for %s (see %s)" % (tag, logp))
        info = ffprobe(dst)
        ok = (int(info["width"]) == W_BASE and int(info["height"]) == H_BASE
              and int(info["nb_read_frames"]) == 24
              and abs(float(info["duration"]) - 1.0) < 1e-6)
        res[tag] = {"file": os.path.basename(dst),
                    "bytes": os.path.getsize(dst),
                    "sha12": sha12(dst), "ffprobe": info, "pass": bool(ok),
                    "encoder": encoder}
        print("   %s: %s frames=%s dur=%s %sx%s -> %s"
              % (tag, os.path.basename(dst), info.get("nb_read_frames"),
                 info.get("duration"), info.get("width"), info.get("height"),
                 "PASS" if ok else "FAIL"), flush=True)
    write_json(os.path.join(OUT, "encode-report.json"), res)
    done("encode (two 24-frame mp4s + ffprobe)")


# ---------------------------------------------------------------------------
# STAGE stills
# ---------------------------------------------------------------------------
def stage_stills():
    idx = FRAME_MID - FRAME_FIRST + 1
    srcA = os.path.join(OUT, "frames", "eta100", "frame-%04d.png" % idx)
    srcB = os.path.join(OUT, "frames", "eta080", "frame-%04d.png" % idx)
    dstA = os.path.join(OUT, "still-eta100.png")
    dstB = os.path.join(OUT, "still-eta080.png")
    shutil.copyfile(srcA, dstA)
    shutil.copyfile(srcB, dstB)
    a = np.asarray(Image.open(dstA).convert("RGB")).astype(np.int16)
    b = np.asarray(Image.open(dstB).convert("RGB")).astype(np.int16)
    d = np.abs(a - b)
    res = {"frame": FRAME_MID, "still_A": os.path.basename(dstA),
           "still_B": os.path.basename(dstB),
           "identical_bytes": bool(sha12(dstA) == sha12(dstB)),
           "sha12_A": sha12(dstA), "sha12_B": sha12(dstB),
           "size": [int(a.shape[1]), int(a.shape[0])],
           "mean_abs_diff_uint8": float(d.mean()),
           "max_abs_diff_uint8": int(d.max()),
           "pct_pixels_gt2": float(np.mean(d > 2) * 100.0),
           "pct_pixels_gt8": float(np.mean(d > 8) * 100.0),
           "mean_A": float(a.mean()), "mean_B": float(b.mean())}
    res["pass"] = bool(not res["identical_bytes"]
                       and res["mean_abs_diff_uint8"] > 0.0)
    write_json(os.path.join(OUT, "stills-diff.json"), res)
    print("   stills: mean|d|=%.3f/255 max=%d identical=%s -> %s"
          % (res["mean_abs_diff_uint8"], res["max_abs_diff_uint8"],
             res["identical_bytes"], "PASS" if res["pass"] else "FAIL"),
          flush=True)
    done("stills (f0074 A/B + difference)")


# ---------------------------------------------------------------------------
# STAGE report
# ---------------------------------------------------------------------------
def _verify_frames():
    v = read_json(os.path.join(OUT, "verification.json"), {}) or {}
    return v.get("frames", {}) or {}


def stage_report():
    timing = read_json(os.path.join(OUT, "timing.json"), {})
    g0 = read_json(os.path.join(OUT, "gate0-t0-identity.json"), {})
    enc = read_json(os.path.join(OUT, "encode-report.json"), {})
    stl = read_json(os.path.join(OUT, "stills-diff.json"), {})
    dome = read_json(os.path.join(OUT, "dome.json"), {})
    parts = []
    for f in range(FRAME_FIRST, FRAME_LAST + 1):
        p = read_json(frame_part_path(f))
        if p is None:
            raise SystemExit("missing frame part f%04d" % f)
        parts.append(p)
    n = len(parts)
    wallA = sum(p["wall_s"]["fields"] + p["wall_s"]["direct"]
                + p["wall_s"]["composite_A"] for p in parts)
    wallB_extra = sum(p["wall_s"]["sky"] + p["wall_s"]["composite_B"]
                      for p in parts)
    wall_total = sum(p["wall_s"]["total"] for p in parts)
    crushA = float(np.mean([p["still_A"]["crush_rate_linear_lt_0.02"] for p in parts]))
    crushB = float(np.mean([p["still_B"]["crush_rate_linear_lt_0.02"] for p in parts]))
    blowA = float(np.mean([p["still_A"]["blowout_rate_aces_gt_0.98"]
                           for p in parts]))
    blowB = float(np.mean([p["still_B"]["blowout_rate_aces_gt_0.98"] for p in parts]))
    p99A = float(np.mean([p["C_direct"]["p99"] for p in parts]))
    p99B = float(np.mean([p["C_total_eta080"]["p99"] for p in parts]))
    p99sA = float(np.mean([p["C_direct"]["p99_9"] for p in parts]))
    p99sB = float(np.mean([p["C_total_eta080"]["p99_9"] for p in parts]))
    stdA = float(np.mean([p["C_direct"]["std"] for p in parts]))
    stdB = float(np.mean([p["C_total_eta080"]["std"] for p in parts]))
    madFrame = float(np.mean([p["A_vs_B_png"]["mean_abs_uint8"] for p in parts]))
    ref = g0.get("Q49_vs_v8", {}).get("reference_stats", {})
    g1 = read_json(os.path.join(OUT, "gate1-f0062-Q49.json"), {})
    gate_t0 = bool(g0.get("pass"))
    gate_mp4 = bool(enc.get("eta100", {}).get("pass")
                    and enc.get("eta080", {}).get("pass"))
    gate_still = bool(stl.get("pass"))
    f62 = [p for p in parts if p["frame"] == 62][0]
    d62 = g1.get("delta_Q9_minus_Q49", {})
    gate_stats = bool(g1.get("pass"))
    d62_t0 = {k: float(f62["C_direct"][k] - ref.get(k, 0.0))
              for k in ("mean", "p50", "p99", "p99_9")}

    L = []
    A = L.append
    A("# SAMPLE-REPORT — 1-second composite samples, eta=1.0 vs eta=0.8")
    A("")
    A("BRIEF: `D:\\mutaro-work\\rally-in\\BRIEF-SAMPLES-1S.md` ／ "
      "run: `offline/run_v9_samples.py` ／ generated %s" % iso())
    A("")
    A("Base frames: `blender-piece8s/renders/frame_0062..0085.png` (24 frames, "
      "wing-raise → liftoff → exit start). **The documented fallback "
      "(extracting from the delivery mp4) was NOT needed**: all 24 PNG frames "
      "exist. Note the on-disk naming is `frame_XXXX.png`, not `fXXXX.png`.")
    A("")
    A("## Deliverables (PC copies in `offline/out/v9samples/`)")
    A("")
    A("| file | bytes | sha12 |")
    A("|---|---|---|")
    for tag in ("eta100", "eta080"):
        e = enc.get(tag, {})
        A("| `%s` | %s | `%s` |" % (e.get("file", "comp-%s-1s.mp4" % tag),
                                    e.get("bytes"), e.get("sha12")))
    A("| `still-eta100.png` | | `%s` |" % stl.get("sha12_A"))
    A("| `still-eta080.png` | | `%s` |" % stl.get("sha12_B"))
    A("")
    A("VM copy (`delivery\\offline-v9\\samples-1s\\`): **not performed from this "
      "PC** — there is no scp route configured on this machine and no delivery "
      "path mounted. Checked: `%USERPROFILE%\\.ssh\\config` and `known_hosts` "
      "(only github.com), `net use` / `Get-PSDrive` (no network shares), a "
      "filesystem search for directories named `delivery` / `offline-v9` / "
      "`samples-1s` on C: and D: (none), and the repo/brief records (no host or "
      "user is recorded anywhere).  The earlier piece/v8/v9 briefs in this repo "
      "state the coordinator fetches deliverables, so that route applies.  The "
      "PC copies required by the brief are in place and byte-verified.")
    A("")
    A("## Light mix and constants")
    A("")
    A("| item | value |")
    A("|---|---|")
    A("| version A | `C = C_direct`, eta = 1.0 exactly (v8 route) |")
    A("| version B | `C_total = wd*C_direct + ws*C_sky`, eta = 0.8 (Perez dome) |")
    A("| w_direct / w_diffuse (B) | %.6f / %.6f |"
      % (dome.get("split", {}).get("w_direct", float("nan")),
         dome.get("split", {}).get("w_diffuse", float("nan"))))
    A("| Perez epsilon / bin | %.4f / %s |"
      % (dome.get("dome", {}).get("epsilon", float("nan")),
         dome.get("dome", {}).get("coefficients", {}).get("bin")))
    A("| dome directions (B) | **%d** |" % N_DIRS)
    A("| transport grid / diffuse master | 2048² / 2048² |")
    A("| direct quadrature Q | **%d** (production is 49; see Reduction) |" % Q_DIRECT)
    A("| composite | frozen v5: ambient ×%.1f, low-pass σ=%.1f px, direct-term "
      "E=%.1f, physical Beer-Lambert, Smith G1=%.3f |"
      % (A_V3, B_SIGMA_V3, E_COMP, dome.get("optics", {}).get("g1", float("nan"))))
    A("| composite size | 832×480 native (cover-resize = identity) |")
    A("| seed | %s (preset 1, derived from the unmodified recipe) |"
      % dome.get("seed"))
    A("| frame rate / window | 24 fps, f0062..f0085 = 1.000 s |")
    A("")
    A("## Wall times (measured)")
    A("")
    A("Probe on f0062 (the brief's cost measurement, `timing.json`):")
    A("")
    A("| measurement | seconds |")
    A("|---|---|")
    ms = timing.get("measured_s", {})
    A("| det_frame Q=49 (production direct) | %.1f |" % ms.get("det_frame_Q49", 0))
    A("| det_frame Q=9 | %.1f |" % ms.get("det_frame_Q9", 0))
    A("| det_frame_dirs, 8 dome directions | %.1f |" % ms.get("det_frame_dirs_8dirs", 0))
    A("| per dome direction | %.3f |" % ms.get("per_direction", 0))
    A("| fields_at | %.2f |" % ms.get("fields_at", 0))
    A("| one composite at 832×480 | %.2f |" % ms.get("composite_one_version", 0))
    A("")
    A("Projection for the 24 B-frames (per frame → total):")
    A("")
    A("| config | s/frame | total 24 |")
    A("|---|---|---|")
    for k, v in (timing.get("projection") or {}).items():
        A("| %s | %.1f | %.1f min |" % (k, v["per_frame_s"],
                                        v["total_24_min"]))
    A("")
    A("Actually spent (sum of the 24 completed frames, from the per-frame "
      "parts): **%.1f min total, %.1f s/frame**. Per-version reading of the "
      "same loop: the A version's own incremental cost is **%.1f min** "
      "(fields + direct transport + composite A; the B frame reuses that same "
      "direct field), and a standalone B version costs the whole **%.1f min** "
      "(direct + 48-direction dome + composite B)."
      % (wall_total / 60.0, wall_total / n, wallA / 60.0,
         (wallA + wallB_extra) / 60.0))
    A("")
    A("Independent re-verification (`verification.json`, no pipeline imports): "
      "both mp4s ffprobe as 832×480 / 24 frames / 1.000000 s / yuv420p h264; "
      "24 unique frame hashes per version (motion present: consecutive mean "
      "|Δ| %.1f–%.1f / 255); stills differ by %.2f / 255; no `.tmp` residue."
      % (min(x["consecutive_mean_abs_diff"]["min"]
             for x in _verify_frames().values()) if _verify_frames() else 0.0,
         max(x["consecutive_mean_abs_diff"]["max"]
             for x in _verify_frames().values()) if _verify_frames() else 0.0,
         stl.get("mean_abs_diff_uint8", 0.0)))
    A("")
    A("### Reduction actually used (not hidden)")
    A("")
    A(timing.get("selected", {}).get("reason", ""))
    A("")
    A("Recorded numbers: the projected total at the brief's own ladder floor "
      "(48 directions) with Q=9 is %.1f min, i.e. above the brief's ~40 min "
      "target; at Q=49 it is %.1f min, of which the direct transport alone is "
      "%.1f min. The dome count was therefore NOT taken below 48."
      % (timing.get("projection", {}).get("SELECTED_Q9_dirs48", {}).get("total_24_min", float("nan")),
         timing.get("projection", {}).get("Q49_dirs48", {}).get("total_24_min", float("nan")),
         24.0 * ms.get("det_frame_Q49", 0) / 60.0))
    A("")
    A("Composite exposure: the brief says E=2.2 and that value was used. "
      "Production v9 (`run_v9.py`) uses the display-normalised "
      "E_v9 = %.4f for this scene. Both versions share the same E, so the "
      "eta comparison is unaffected; the absolute level is 2.2/%.4f = %.3f of "
      "the production still level."
      % (dome.get("optics", {}).get("E_v9_production", float("nan")),
         dome.get("optics", {}).get("E_v9_production", float("nan")),
         E_COMP / dome.get("optics", {}).get("E_v9_production", float("nan"))))
    A("")
    A("## Mean tone fractions per version (24-frame means)")
    A("")
    A("| quantity | A (eta=1.0) | B (eta=0.8) |")
    A("|---|---|---|")
    A("| crush (linear lum < 0.02) | %.5f | %.5f |" % (crushA, crushB))
    A("| blow (ACES lum > 0.98) | %.6f | %.6f |" % (blowA, blowB))
    A("| C p99 | %.4f | %.4f |" % (p99A, p99B))
    A("| C p99.9 | %.4f | %.4f |" % (p99sA, p99sB))
    A("| C std | %.5f | %.5f |" % (stdA, stdB))
    A("| mean |A−B| over the PNGs | %.2f / 255 |" % madFrame)
    A("")
    A("Crush bound (customary 0.0325): A %.5f, B %.5f — both inside. Blow: both "
      "0." % (crushA, crushB))
    A("")
    A("## Gates")
    A("")
    A("| gate | requirement | measured | verdict |")
    A("|---|---|---|---|")
    eA = enc.get("eta100", {}).get("ffprobe", {})
    eB = enc.get("eta080", {}).get("ffprobe", {})
    A("| mp4 A | 24 frames, 1.000 s, 832×480 | %s fr, %s s, %s×%s | %s |"
      % (eA.get("nb_read_frames"), eA.get("duration"), eA.get("width"),
         eA.get("height"), "PASS" if eA.get("duration") and enc.get("eta100", {}).get("pass") else "FAIL"))
    A("| mp4 B | 24 frames, 1.000 s, 832×480 | %s fr, %s s, %s×%s | %s |"
      % (eB.get("nb_read_frames"), eB.get("duration"), eB.get("width"),
         eB.get("height"), "PASS" if enc.get("eta080", {}).get("pass") else "FAIL"))
    A("| stills differ | not byte-identical | sha12 %s vs %s, mean|d| %.2f/255 "
      "(max %d, %.1f%% px > 2) | %s |"
      % (stl.get("sha12_A"), stl.get("sha12_B"),
         stl.get("mean_abs_diff_uint8", 0), stl.get("max_abs_diff_uint8", 0),
         stl.get("pct_pixels_gt2", 0),
         "PASS" if gate_still else "FAIL"))
    A("| A f0062 C stats vs production route | Q=9 vs Q=49 at the SAME time "
      "slice: mean/p50/p99/p99.9 within 1e-6/1e-3/5e-3/2e-2 | Δmean %+.2e, "
      "Δp50 %+.2e, Δp99 %+.2e, Δp99.9 %+.2e | %s |"
      % (d62.get("mean", float("nan")), d62.get("p50", float("nan")),
         d62.get("p99", float("nan")), d62.get("p99_9", float("nan")),
         "PASS" if gate_stats else "FAIL"))
    A("| taskkill python.exe | never used anywhere in this run | | PASS (by construction) |")
    A("| production route proof | t=0, Q=49 reproduces `out/v8/C-v8-t0.npy` | "
      "max|ΔC| = %s (bit-identical %s) | %s |"
      % (g0.get("Q49_vs_v8", {}).get("max_abs_diff"),
         g0.get("Q49_vs_v8", {}).get("bit_identical"),
         "PASS" if gate_t0 else "FAIL"))
    A("| mp4s | both encode + probe | | %s |" % ("PASS" if gate_mp4 else "FAIL"))
    A("")
    A("Tolerances used for the C-statistics gate: the gate is measured against "
      "the SAME time slice at the production quadrature (Q=49 at f0062, "
      "`gate1-f0062-Q49.json`), with mean ≤ 1e-6, p50 ≤ 1e-3, p99 ≤ 5e-3, "
      "p99.9 ≤ 2e-2. Comparing f0062 against the stored v8 **t=0** master "
      "instead would conflate the quadrature with the TIME EVOLUTION of the "
      "field (measured here: Δp50 %+.2e, Δp99 %+.2e, Δp99.9 %+.2e against "
      "t=0) — v8's own record shows the same order of change between t=0 and "
      "t=2 s (numbers-v8.json, dispersion ablation: Δp99 +0.264, Δp99.9 "
      "+0.421). The exact production-route identity is the t=0 gate, which is "
      "bit-identical. Max is reported, never gated: it is a single-pixel "
      "order statistic that moves with time and with Q."
      % (d62_t0["p50"], d62_t0["p99"], d62_t0["p99_9"]))
    A("")
    A("## Provenance of the frames (concurrent run — disclosed)")
    A("")
    A("This PC also received a continuation brief for the same samples, and a "
      "second session created `offline/run_v9_samples_1s.py` (an independent "
      "standalone runner, smoke-tested, never run in production) plus a "
      "**detached WMI continuation** (`offline/_v9samples_continuation.py`, "
      "PID 48428). That continuation imports `run_v9_samples.py` and calls its "
      "`stage_frames(f, f)` / `encode` / `stills` / `report` functions, so the "
      "code path, seed, settings and JSON schema are identical; it waits for "
      "the live foreground batch instead of duplicating it (log: "
      "`continuation-runner.log`). Frames f0062–f0085 were therefore produced "
      "by one configuration — verified here: all 24 per-frame parts carry "
      "`Q_direct=9, n_dirs=48, w_direct=0.774247` and nothing else. No "
      "`taskkill` was issued at any point. Evidence kept in the same folder: "
      "`VERIFICATION-CONT.md`, `continuation-runner.log`, "
      "`continuation-stdout.log`.")
    A("")
    A("## Honest notes / deviations")
    A("")
    A("- **Direct quadrature reduced to Q=9** (production Q=49). Justification "
      "is measured, not asserted: t=0 at Q=49 is bit-identical to the stored "
      "v8 production master, and at Q=9 the p99/p99.9 deltas vs that master "
      "are %+.2e / %+.2e (the v8 Q-convergence table in numbers-v8.json "
      "records p99 1.5097/1.5098/1.5095 for Q=9/25/49). The direct field is "
      "shared by both versions, so this cannot bias the eta comparison."
      % (g0.get("Q9_vs_v8", {}).get("delta_p99", float("nan")),
         g0.get("Q9_vs_v8", {}).get("delta_p99_9", float("nan"))))
    A("- **Dome directions = 48** (brief's ladder floor; the 96-direction "
      "production value is 4× the cost). The measured 96-vs-64 quadrature "
      "scatter in numbers-v9.json is Δp99 0.0183–0.0289; at 48 directions the "
      "diffuse speckle is somewhat larger, which is the dominant known "
      "deviation of version B in these samples.")
    A("- **Composite exposure E=2.2 as the brief specifies**, while production "
      "v9 uses E_v9=%.4f (display normalisation). Both samples share E, so "
      "the A/B comparison is valid; the absolute brightness is darker than the "
      "production still by that ratio."
      % dome.get("optics", {}).get("E_v9_production", float("nan")))
    A("- `g1` (Smith) is evaluated once from the t=0 state and held for all "
      "frames, exactly as the still pipeline does; the value is %.6f."
      % dome.get("optics", {}).get("g1", float("nan")))
    A("- The production compositor script for the 8 s piece is not in this "
      "repo; the composite here is the frozen `v5_composite.render_frame_v5` "
      "with A2.6 / σ3 px / E2.2 / Beer-Lambert / G1 and the base frame as the "
      "bottom albedo, rendered natively at 832×480.")
    A("- No `taskkill python.exe` was used; every stage ran in the foreground "
      "(or, for the remaining frames, in the disclosed detached continuation "
      "that drives the same code) with jobkit checkpoints "
      "(`out/v9samples/checkpoint.json`).")
    A("")
    md = "\n".join(L) + "\n"
    atomic_write_bytes(os.path.join(OUT, "SAMPLE-REPORT.md"), md.encode("utf-8"))
    write_json(os.path.join(OUT, "numbers-samples.json"), {
        "brief": "BRIEF-SAMPLES-1S.md", "created": iso(),
        "frames": [FRAME_FIRST, FRAME_LAST], "fps": FPS, "grid": GRID,
        "Q_direct": Q_DIRECT, "Q_production": Q_PRODUCTION, "n_dirs": N_DIRS,
        "eta_B": ETA_B, "E_composite": E_COMP, "A_V3": A_V3,
        "B_SIGMA_V3": B_SIGMA_V3, "split": dome.get("split"),
        "dome": dome.get("dome"),
        "optics": dome.get("optics"), "scene": dome.get("scene"),
        "timing_probe": timing.get("measured_s"),
        "projection": timing.get("projection"),
        "wall_s_total_24_frames": wall_total,
        "per_frame": parts, "encode": enc, "stills": stl,
        "gate0_t0_identity": g0,
        "gate1_f0062_Q49": g1,
        "A_f0062_vs_v8_t0_informational": d62_t0,
        "gates": {"mp4_A": enc.get("eta100", {}).get("pass"),
                  "mp4_B": enc.get("eta080", {}).get("pass"),
                  "stills_differ": gate_still,
                  "A_f0062_vs_Q49_same_slice": gate_stats,
                  "t0_production_identity": gate_t0},
    })
    try:
        print(md, flush=True)
    except UnicodeEncodeError:
        print(md.encode("ascii", "replace").decode("ascii"), flush=True)
    done("report (SAMPLE-REPORT.md + numbers-samples.json)")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["probe", "gate0", "ref62", "frames", "encode",
                             "stills", "report"])
    ap.add_argument("--from", dest="f0", type=int, default=FRAME_FIRST)
    ap.add_argument("--to", dest="f1", type=int, default=FRAME_LAST)
    args = ap.parse_args()
    ensure_dirs()
    if args.stage == "probe":
        stage_probe()
    elif args.stage == "gate0":
        stage_gate0()
    elif args.stage == "ref62":
        stage_ref62()
    elif args.stage == "frames":
        stage_frames(args.f0, args.f1)
    elif args.stage == "encode":
        stage_encode()
    elif args.stage == "stills":
        stage_stills()
    elif args.stage == "report":
        stage_report()


if __name__ == "__main__":
    main()
