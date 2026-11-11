"""Standalone runner for BRIEF-SAMPLES-1S (+ "-CONT" execution form).

1-second composite samples over the moth base frames f0062..f0085:

  version A  eta = 1.0  (direct only, the v8 route; C_total IS C_dir)
  version B  eta = 0.8  (v9 route: w_d C_dir + w_s C_sky, Perez dome with
                         48 directions, diffuse master on the 2048^2 grid,
                         which is the transport grid so no upsample is needed)

Both versions share everything except the light mix: the SAME frame's water
state, the SAME direct transport (`det.det_frame`, Q=49) and the SAME frozen
composite (`v5_composite.render_frame_v5`: A=2.6 ambient, lowpass sigma 3 px,
E on the direct term, Beer-Lambert on the direct path, turbidity blur,
substrate gain, G1).  The only difference between the two PNGs of a frame is
C_total.

Deterministic / idempotent:
  * a frame whose A and B PNGs both exist is skipped (also recorded in the
    jobkit checkpoint);
  * every PNG is written through jobkit's atomic_write_bytes (tmp + fsync +
    replace), so a killed process can never leave a half-written frame;
  * one line per frame is appended to progress.log and frames.jsonl, and the
    jobkit checkpoint.json is updated once per completed frame.

Products in <out>/ (default out/v9samples/):
  eta100/frame-XXXX.png   eta080/frame-XXXX.png
  still-eta100.png  still-eta080.png          (the mid-window frame f0074)
  checkpoint.json  progress.log  frames.jsonl  runner.log  runner.pid
  numbers-v9samples.json  SAMPLE-REPORT.md

Usage (from caustics-filter/offline/):
  python run_v9_samples_1s.py --frames 62-85
  python run_v9_samples_1s.py --smoke --frames 62 --out out/v9samples-smoke
Encoding to mp4 is NOT part of this script (a follow-up session does it).
"""
import argparse
import datetime
import gc
import io
import json
import os
import platform
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from env import resolve_environment                        # noqa: E402
from det import det_frame                                  # noqa: E402
from state import derive_seed                              # noqa: E402
from image import srgb_to_linear                           # noqa: E402
from v5_physics import (water_attenuation, beer_lambert,    # noqa: E402
                        smith_g1_beckmann,
                        scattering_blur_sigma_px,
                        substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3,           # noqa: E402
                          image_metrics)
from recipe_v9 import parse_recipe_v9, with_v9              # noqa: E402
from v8_state import V8WaterState                           # noqa: E402
import v9_sky as SK                                         # noqa: E402
import v9_transport as TR                                   # noqa: E402

# the frozen chunk-5 helpers are imported, not copied, so this runner cannot
# drift from the production route
from run_v9 import (scene_optics, combine, make_still,      # noqa: E402
                    c_stats, delta_stats, now_iso)

from checkpoint import Checkpoint, atomic_write_bytes       # noqa: E402

FPS = 24
BASE_W, BASE_H = 832, 480
REFERENCE_DIRECT_STATS = {          # recorded by the previous session (f0062)
    "mean": 1.0000002, "p99": 1.5482, "max": 2.1742,
    "tolerance": {"mean": 1e-5, "p99": 1e-2, "max": 1e-2},
    "source": ("BRIEF-SAMPLES-1S-CONT: the previous session's measured direct "
               "C statistics at the window; used only as a sanity check that "
               "the v8 direct route was actually used"),
}


def parse_frames(spec):
    out = []
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def png_bytes(img):
    buf = io.BytesIO()
    Image.fromarray(np.asarray(img, dtype=np.uint8)).save(buf, format="PNG")
    return buf.getvalue()


def flist(x):
    """numpy scalar/array -> plain list of floats (for json)."""
    return [float(v) for v in np.asarray(x, dtype=np.float64).ravel()]


def load_base_frame(photo_dir, f, W, H):
    """The moth base frame as linear albedo, resized/cropped to W x H."""
    path = os.path.join(photo_dir, "frame_%04d.png" % f)
    if not os.path.exists(path):
        raise SystemExit("base frame missing: %s" % path)
    im = Image.open(path).convert("RGB")
    if im.size != (W, H):
        # cover-fit (same aspect for the production 832x480 frames)
        sc = max(W / float(im.width), H / float(im.height))
        nw, nh = max(W, int(round(im.width * sc))), max(H, int(round(im.height * sc)))
        im = im.resize((nw, nh), Image.LANCZOS)
        x0 = (nw - W) // 2
        y0 = (nh - H) // 2
        im = im.crop((x0, y0, x0 + W, y0 + H))
    srgb = np.asarray(im, dtype=np.float64) / 255.0
    return srgb_to_linear(srgb), path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--frames", default="62-85")
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--dirs", type=int, default=48,
                    help="dome directions for version B (the CONT note fixes "
                         "48; the count is recorded in the report)")
    ap.add_argument("--diffuse-grid", type=int, default=2048)
    ap.add_argument("--eta-a", type=float, default=1.0)
    ap.add_argument("--eta-b", type=float, default=0.8)
    ap.add_argument("--fetch", type=float, default=100.0)
    ap.add_argument("--photo-dir", default="../blender-piece8s/renders")
    ap.add_argument("--out", default="out/v9samples")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--force", action="store_true",
                    help="re-render frames even if both PNGs already exist")
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--smoke", action="store_true",
                    help="256^2 / Q=9 / 8 dome directions (separate out dir)")
    ap.add_argument("--plan-only", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.grid = 256
        args.q = 9
        args.dirs = 8
        args.diffuse_grid = 256
        if args.out == "out/v9samples":
            args.out = "out/v9samples-smoke"

    out = args.out
    os.makedirs(out, exist_ok=True)
    dir_a = os.path.join(out, "eta100")
    dir_b = os.path.join(out, "eta080")
    os.makedirs(dir_a, exist_ok=True)
    os.makedirs(dir_b, exist_ok=True)
    runlog = os.path.join(out, "runner.log")
    prolog = os.path.join(out, "progress.log")
    jlpath = os.path.join(out, "frames.jsonl")
    ck = Checkpoint(os.path.join(out, "checkpoint.json"))

    # ---- no console under a WMI launch: everything goes to files ----
    logf = open(runlog, "a", encoding="utf-8")

    def log(msg):
        line = "[%s] %s" % (now_iso(), msg)
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()
        try:
            os.fsync(logf.fileno())
        except OSError:
            pass

    def prog(msg):
        with open(prolog, "a", encoding="utf-8") as f:
            f.write("%s %s\n" % (now_iso(), msg))
            f.flush()
            os.fsync(f.fileno())

    def jrec(rec):
        with open(jlpath, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False, default=float) + "\n")
            f.flush()
            os.fsync(f.fileno())

    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "runner.pid"), "w", encoding="utf-8") as f:
        f.write("%d\n" % os.getpid())

    frames = parse_frames(args.frames)
    if args.max_frames:
        frames = frames[:int(args.max_frames)]

    t_all = time.time()
    log("=" * 66)
    log("v9 1-second samples runner start  pid=%d  host=%s  python=%s"
        % (os.getpid(), platform.node(), platform.python_version()))
    log("argv: %s" % " ".join(sys.argv))
    log("frames=%s grid=%d q=%d dirs=%d diffuse_grid=%d etaA=%.3f etaB=%.3f"
        % (frames, args.grid, args.q, args.dirs, args.diffuse_grid,
           args.eta_a, args.eta_b))
    prog("runner_start pid=%d frames=%s" % (os.getpid(), args.frames))

    if args.plan_only:
        log("plan-only: exiting before any work")
        return

    # ==================================================================
    # scene / state / optics (identical to run_v9 STAGE 2)
    # ==================================================================
    with open(args.recipe, "r", encoding="utf-8") as f:
        raw = json.load(f)
    # the B recipe (eta = 0.8, manual) as written by run_v9; the seed is
    # derived from the UNMODIFIED recipe so the water realisation is the one
    # the production masters used
    scene = with_v9(raw, eta=float(args.eta_b), mode="manual",
                    fetch_m=float(args.fetch))
    env = resolve_environment(scene)
    cfg = parse_recipe_v9(scene)          # manual mode: no fetch, eta from file
    seed = (int(args.seed) if args.seed is not None else derive_seed(raw, 1))
    grid = int(args.grid)
    S_patch = 4.0 * float(env["specL"])
    log("recipe=%s schema=%s wind=%.2f depth=%.2f sun_alt=%.3f seed=%d"
        % (raw.get("name"), raw.get("schema"), env["windSpeedMS"],
           env["depthM"], env["sun"]["altitudeDeg"], seed))
    log("patch S=%.4f m  dx=%.4f mm  recipe eta=%.3f (%s)"
        % (S_patch, S_patch / grid * 1e3, cfg["sky"]["eta"],
           cfg["sky"]["source"]))

    t0 = time.time()
    st = V8WaterState(scene, env, N=grid, seed=seed, preset_index=args.preset,
                      fetch_m=float(cfg["fetch"]["fetch_m"]),
                      fetch_source=cfg["fetch"]["source"])
    t_state = time.time() - t0
    f0 = st.fields_at(0.0)
    f0["S"] = st.S
    re_state = st.realised_stats(f0)
    del f0
    gc.collect()

    att = water_attenuation(env["turbid"])
    opt = scene_optics(env, cfg, att, st, grid)
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1 = smith_g1_beckmann(re_state["rms_slope"], sun_zenith)
    G1 = float(g1["G1"])
    # the chunk-5 display normalisation (E_V3 = 2.2 renormalised so the
    # flat-water direct level equals the v3 definition); it is the same for
    # both versions, which is what makes the A/B comparison a light-mix one
    E_v9 = E_V3 * (opt["L_v3"] / (opt["L_dir_flat"] * G1))
    log("state %.1fs sigma2=%.6f Hs=%.4f rms_slope=%.4f G1=%.6f"
        % (t_state, re_state["sigma2_rad2"], re_state["Hs_m"],
           re_state["rms_slope"], G1))
    log("optics beer_down=%s beer_up=%s blur=%.4f px gain=%s"
        % ([round(x, 5) for x in opt["beer"]["down"]],
           [round(x, 5) for x in opt["beer"]["up"]],
           opt["blur"]["sigma_px"], [round(x, 4) for x in opt["gain"]["gain_rgb"]]))
    log("composite: A=%.2f B_sigma=%.2f E_v3=%.2f -> E used=%.6f"
        % (A_V3, B_SIGMA_V3, E_V3, E_v9))

    split_a = SK.horizontal_weights(args.eta_a, env["sun"]["altitudeDeg"])
    split_b = SK.horizontal_weights(args.eta_b, env["sun"]["altitudeDeg"])
    dome = SK.sky_dome(args.eta_b, env["sunDir"], n_dirs=int(args.dirs))
    log("split A eta=%.3f w_d=%.6f w_s=%.6f" % (args.eta_a,
                                                split_a["w_direct"],
                                                split_a["w_diffuse"]))
    log("split B eta=%.3f w_d=%.6f w_s=%.6f" % (args.eta_b,
                                                split_b["w_direct"],
                                                split_b["w_diffuse"]))
    log("dome B: dirs=%d eps=%.4f bin=%d w_eff=%.2f w_max=%.5f circ5=%.5f"
        % (dome["n_dirs"], dome["epsilon"], dome["coefficients"]["bin"],
           dome["w_effective_n"], dome["w_max"],
           dome["circumsolar_fraction_5deg"]))

    # ==================================================================
    # frames
    # ==================================================================
    done_ids = set()
    st_ck = ck.state()
    if st_ck and isinstance(st_ck.get("completed_ids"), list):
        done_ids = set(st_ck["completed_ids"])

    per_frame = []
    crush = {"A": [], "B": []}
    blow = {"A": [], "B": []}
    for f in frames:
        uid = "frame-%04d" % f
        pa = os.path.join(dir_a, "frame-%04d.png" % f)
        pb = os.path.join(dir_b, "frame-%04d.png" % f)
        have_a = os.path.exists(pa) and os.path.getsize(pa) > 0
        have_b = os.path.exists(pb) and os.path.getsize(pb) > 0
        if have_a and have_b and not args.force:
            log("frame %04d: skip (both PNGs exist)" % f)
            prog("frame=%04d status=skip (both PNGs exist)" % f)
            continue
        tv = time.time()
        t = f / float(FPS)
        try:
            t0 = time.time()
            fields = st.fields_at(t)
            fields["S"] = st.S
            t_fields = time.time() - t0

            photo_lin, base_path = load_base_frame(args.photo_dir, f,
                                                   BASE_W, BASE_H)
            PH, PW, _ = photo_lin.shape

            # ---- ONE direct transport shared by both versions ----
            t0 = time.time()
            C_dir, _ = det_frame(fields, st.S, grid, env, int(args.q))
            t_dir = time.time() - t0
            dir_stats = c_stats(C_dir)

            # version A: eta = 1 -> the sky weight is exactly 0, no array add
            C_A = combine(C_dir, None, split_a)

            # version B: the Perez dome on the same transport
            t0 = time.time()
            C_sky, _ = TR.det_frame_dirs(fields, st.S, int(args.diffuse_grid),
                                         env, dome["dirs"], dome["w"],
                                         tag="dome")
            t_sky = time.time() - t0
            if C_sky.shape[0] != grid:
                C_sky_up = TR.fft_upsample(C_sky, grid)
                upsampled = True
            else:
                C_sky_up = C_sky
                upsampled = False
            C_B = combine(C_dir, C_sky_up, split_b)
            tot_stats = c_stats(C_B)
            sky_stats = c_stats(C_sky)

            t0 = time.time()
            imgA, mA = make_still(fields, C_A, photo_lin, env, PW, PH,
                                  E_v9, G1, opt)
            t_compA = time.time() - t0
            t0 = time.time()
            imgB, mB = make_still(fields, C_B, photo_lin, env, PW, PH,
                                  E_v9, G1, opt)
            t_compB = time.time() - t0

            atomic_write_bytes(pa, png_bytes(imgA))
            atomic_write_bytes(pb, png_bytes(imgB))

            dA = imgA.astype(np.int16)
            dB = imgB.astype(np.int16)
            dd = np.abs(dA - dB)
            dt = time.time() - tv
            crush["A"].append(float(mA["crush_rate_linear_lt_0.02"]))
            crush["B"].append(float(mB["crush_rate_linear_lt_0.02"]))
            blow["A"].append(float(mA["blowout_rate_aces_gt_0.98"]))
            blow["B"].append(float(mB["blowout_rate_aces_gt_0.98"]))

            rec = {
                "frame": int(f), "t_s": float(t), "base_frame": base_path,
                "size": [int(PW), int(PH)], "grid": grid, "q": int(args.q),
                "dirs_b": int(args.dirs),
                "diffuse_grid": int(args.diffuse_grid),
                "upsampled": bool(upsampled),
                "timing_s": {"fields_at": t_fields, "det_frame": t_dir,
                             "dome_B": t_sky, "composite_A": t_compA,
                             "composite_B": t_compB, "frame_total": dt},
                "C_dir": dir_stats, "C_sky": sky_stats, "C_total_B": tot_stats,
                "metrics_A": {
                    "crush": float(mA["crush_rate_linear_lt_0.02"]),
                    "blow": float(mA["blowout_rate_aces_gt_0.98"]),
                    "linear_mean": float(mA["linear_mean_pre_tonemap"]),
                    "aces_max": float(mA["aces_luminance_max"])},
                "metrics_B": {
                    "crush": float(mB["crush_rate_linear_lt_0.02"]),
                    "blow": float(mB["blowout_rate_aces_gt_0.98"]),
                    "linear_mean": float(mB["linear_mean_pre_tonemap"]),
                    "aces_max": float(mB["aces_luminance_max"])},
                "image_AB": {"mean_abs_diff_uint8": float(dd.mean()),
                             "max_abs_diff_uint8": int(dd.max()),
                             "mean_A_uint8": float(imgA.mean()),
                             "mean_B_uint8": float(imgB.mean()),
                             "byte_identical": bool(dd.max() == 0)},
                "C_AB_delta": {k: float(tot_stats[k] - dir_stats[k])
                               for k in ("mean", "p99", "p99_9", "std",
                                         "max")},
            }
            per_frame.append(rec)
            jrec(rec)
            if uid not in done_ids:
                ck.unit_done(uid)
                done_ids.add(uid)
            prog("frame=%04d t=%.4fs rendered dt=%.1fs (dir %.0fs dome %.0fs) "
                 "Cmean=%.7f Cp99=%.4f Cmax=%.3f imgDiff=%.3f"
                 % (f, t, dt, t_dir, t_sky, dir_stats["mean"],
                    dir_stats["p99"], dir_stats["max"], dd.mean()))
            log("frame %04d done  %.1fs (fields %.1f / direct %.1f / dome %.1f "
                "/ comp %.1f+%.1f)  C_dir mean=%.7f p99=%.4f max=%.3f"
                % (f, dt, t_fields, t_dir, t_sky, t_compA, t_compB,
                   dir_stats["mean"], dir_stats["p99"], dir_stats["max"]))
            del fields, photo_lin, C_dir, C_A, C_sky, C_sky_up, C_B
            del imgA, imgB
            gc.collect()
        except Exception as exc:      # keep going: a 2 h detached run must not
            import traceback          # die on one bad frame
            tb = traceback.format_exc()
            log("frame %04d FAILED: %s" % (f, exc))
            log(tb)
            prog("frame=%04d status=error %s" % (f, exc))
            jrec({"frame": int(f), "status": "error", "error": str(exc),
                  "traceback": tb})
            gc.collect()
            continue

    wall = time.time() - t_all
    log("render phase finished: %d frames rendered/seen, wall %.1fs"
        % (len(per_frame), wall))

    # aggregate from the jsonl ledger so a resumed run reports the WHOLE
    # window, not only the frames this process happened to render
    all_recs = {}
    if os.path.exists(jlpath):
        with open(jlpath, "r", encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                if r.get("status") == "error" or "C_dir" not in r:
                    continue
                all_recs[int(r["frame"])] = r
    agg = [all_recs[k] for k in sorted(all_recs)]

    # ==================================================================
    # stills (mid-window frame) + report
    # ==================================================================
    mid = 74
    still_a = os.path.join(out, "still-eta100.png")
    still_b = os.path.join(out, "still-eta080.png")
    src_a = os.path.join(dir_a, "frame-%04d.png" % mid)
    src_b = os.path.join(dir_b, "frame-%04d.png" % mid)
    still_info = {}
    if os.path.exists(src_a) and os.path.exists(src_b):
        with open(src_a, "rb") as f:
            atomic_write_bytes(still_a, f.read())
        with open(src_b, "rb") as f:
            atomic_write_bytes(still_b, f.read())
        ia = np.asarray(Image.open(still_a).convert("RGB")).astype(np.int16)
        ib = np.asarray(Image.open(still_b).convert("RGB")).astype(np.int16)
        d = np.abs(ia - ib)
        still_info = {
            "frame": mid, "files": ["still-eta100.png", "still-eta080.png"],
            "mean_abs_diff_uint8": float(d.mean()),
            "max_abs_diff_uint8": int(d.max()),
            "byte_identical": bool(d.max() == 0 and
                                  open(src_a, "rb").read() ==
                                  open(src_b, "rb").read()),
        }
        log("stills: f%04d mean|dA-dB|=%.3f max=%d"
            % (mid, still_info["mean_abs_diff_uint8"],
               still_info["max_abs_diff_uint8"]))
    else:
        still_info = {"frame": mid, "missing": "frame %04d not rendered yet" % mid}
        log("stills: frame %04d not present (partial run)" % mid)

    # frame count actually on disk
    n_a = len([p for p in os.listdir(dir_a)
               if p.startswith("frame-") and p.endswith(".png")])
    n_b = len([p for p in os.listdir(dir_b)
               if p.startswith("frame-") and p.endswith(".png")])

    def mean_of(v):
        return float(np.mean(v)) if v else None

    # sanity: A/f0062 vs the previous session's recorded direct stats
    ref = dict(REFERENCE_DIRECT_STATS)
    tol = ref.pop("tolerance")
    ref.pop("source")
    sanity = {"reference": ref, "tolerance": tol, "measured": None,
              "pass": None}
    rec62 = all_recs.get(62)
    if rec62 is not None:
        m = {"mean": rec62["C_dir"]["mean"], "p99": rec62["C_dir"]["p99"],
             "max": rec62["C_dir"]["max"]}
        sanity["measured"] = m
        sanity["delta"] = {k: float(m[k] - ref[k]) for k in ref}
        sanity["pass"] = bool(all(abs(m[k] - ref[k]) <= tol[k] for k in ref))

    wall_per_frame = [r["timing_s"]["frame_total"] for r in agg]
    crush_agg = {"A": [r["metrics_A"]["crush"] for r in agg],
                 "B": [r["metrics_B"]["crush"] for r in agg]}
    blow_agg = {"A": [r["metrics_A"]["blow"] for r in agg],
                "B": [r["metrics_B"]["blow"] for r in agg]}
    numbers = {
        "job": "v9samples-1s",
        "brief": "BRIEF-SAMPLES-1S.md (execution form: BRIEF-SAMPLES-1S-CONT.md)",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "recipe": raw.get("name"), "preset": args.preset, "seed": int(seed),
        "frames": {"requested": frames, "fps": FPS,
                   "window_s": [frames[0] / float(FPS), frames[-1] / float(FPS)],
                   "rendered_now": [r["frame"] for r in per_frame],
                   "pngs_on_disk": {"eta100": n_a, "eta080": n_b}},
        "settings": {
            "grid": grid, "q": int(args.q),
            "dome_directions_B": int(args.dirs),
            "diffuse_grid": int(args.diffuse_grid),
            "diffuse_upsampled": bool(per_frame[0]["upsampled"])
            if per_frame else None,
            "eta_A": float(args.eta_a), "eta_B": float(args.eta_b),
            "split_A": split_a, "split_B": split_b,
            "dome_B": {k: dome[k] for k in
                       ("n_dirs", "epsilon", "circumsolar_fraction_5deg",
                        "w_effective_n", "w_max", "horizon_elev_deg")},
            "dome_B_coefficients": dome["coefficients"],
            "composite": {"A_ambient": A_V3, "B_sigma_px": B_SIGMA_V3,
                          "E_v3": E_V3, "E_used": float(E_v9),
                          "g1": G1, "blur_sigma_px": opt["blur"]["sigma_px"],
                          "beer_down": flist(opt["beer"]["down"]),
                          "beer_up": flist(opt["beer"]["up"]),
                          "substrate_gain_rgb": flist(opt["gain"]["gain_rgb"]),
                          "note": ("E is E_V3=2.2 renormalised by "
                                   "L_v3/(L_dir_flat*G1) exactly as run_v9 "
                                   "does (chunk-5 display level); G1 is "
                                   "1.0 at this scene, so the direct term "
                                   "keeps the frozen v3-v4 level")},
            "state": {"sigma2_rad2": re_state["sigma2_rad2"],
                      "Hs_m": re_state["Hs_m"],
                      "rms_slope": re_state["rms_slope"]},
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "depth_m": env["depthM"], "wind_ms": env["windSpeedMS"],
            "patch_S_m": S_patch,
        },
        "wall_s": {"this_process": wall,
                   "state_build": t_state,
                   "mean_per_frame": (float(np.mean(wall_per_frame))
                                      if agg else None),
                   "sum_per_frame_all_frames": (float(np.sum(wall_per_frame))
                                                if agg else None),
                   "projection_24_frames_s": (
                       float(np.mean(wall_per_frame) * 24)
                       if agg else None)},
        "tone": {"crush_mean_A": mean_of(crush_agg["A"]),
                 "crush_mean_B": mean_of(crush_agg["B"]),
                 "blow_mean_A": mean_of(blow_agg["A"]),
                 "blow_mean_B": mean_of(blow_agg["B"]),
                 "crush_bound": 0.0325, "blow_bound": 0.0},
        "stills": still_info,
        "gate_A_f0062_direct_matches_reference": sanity,
        "per_frame": agg,
        "fallbacks": {
            "base_frames": "blender-piece8s/renders/frame_%04d.png (all present)",
            "mp4_extraction": "not needed",
            "diffuse_master": ("computed per frame (no cached master for a "
                               "%d-direction dome at eta=%.2f)"
                               % (args.dirs, args.eta_b)),
        },
        "runner": {"pid": os.getpid(), "host": platform.node(),
                   "python": platform.python_version()},
    }
    num_path = os.path.join(out, "numbers-v9samples.json")
    atomic_write_bytes(num_path,
                       json.dumps(numbers, indent=2, ensure_ascii=False,
                                  default=float).encode("utf-8"))
    log("numbers written: %s" % num_path)

    # ---------------- SAMPLE-REPORT.md ----------------
    L = []
    a = L.append
    a("# 1-second composite samples — eta = 1.0 (A) vs eta = 0.8 (B)")
    a("")
    a("Brief: `BRIEF-SAMPLES-1S.md` / execution form "
      "`BRIEF-SAMPLES-1S-CONT.md`  ")
    a("Generated: %s by `offline/run_v9_samples_1s.py` (pid %d)"
      % (numbers["created_utc"], os.getpid()))
    a("")
    a("## 結果の要約")
    a("")
    a("- 対象フレーム: f%04d..f%04d（24 fps、%.3f s）"
      % (frames[0], frames[-1], frames[-1] / float(FPS)))
    a("- 今回は %d フレームをこの実行で処理（記録は窓の %d フレーム分）"
      % (len(per_frame), len(agg)))
    a("- 書き出し済み PNG: `eta100/` %d 枚 / `eta080/` %d 枚" % (n_a, n_b))
    a("- ドーム方向数（B）: **%d**（CONT 指定値。diffuse master は %d²）"
      % (args.dirs, args.diffuse_grid))
    a("- 合成は凍結のまま: A=%.1f, lowpass sigma %.1f px, E_v3=%.1f"
      % (A_V3, B_SIGMA_V3, E_V3))
    a("  （chunk-5 の表示レベルの正規化を適用した実効 E=%.6f、G1=%.6f）"
      % (E_v9, G1))
    a("")
    a("## 光の混合")
    a("")
    a("| 版 | eta | w_direct | w_diffuse |")
    a("|---|---|---|---|")
    a("| A | %.2f | %.6f | %.6f |" % (args.eta_a, split_a["w_direct"],
                                        split_a["w_diffuse"]))
    a("| B | %.2f | %.6f | %.6f |" % (args.eta_b, split_b["w_direct"],
                                        split_b["w_diffuse"]))
    a("")
    a("A は w_s が厳密に 0 なので配列加算そのものを行わない "
      "→ C_total は C_dir と同一（v8 経路）。")
    a("B は同じフレームの**同じ直達輸送**に、Perez ドーム %d 方向の拡散場を"
      "重み付きで加算する（v9 経路）。" % args.dirs)
    a("")
    a("## 所要時間（実測）")
    a("")
    if agg:
        tf = wall_per_frame
        a("| 指標 | 値 |")
        a("|---|---|")
        a("| フレームあたり平均 | %.1f s |" % float(np.mean(tf)))
        a("| うち直達輸送 (det_frame Q=%d) | %.1f s |"
          % (args.q, float(np.mean([r["timing_s"]["det_frame"]
                                    for r in agg]))))
        a("| うち拡散ドーム (%d 方向, %d²) | %.1f s |"
          % (args.dirs, args.diffuse_grid,
             float(np.mean([r["timing_s"]["dome_B"] for r in agg]))))
        a("| うち合成 2 枚 | %.1f s |"
          % float(np.mean([r["timing_s"]["composite_A"]
                           + r["timing_s"]["composite_B"]
                           for r in agg])))
        a("| 24 フレームへの外挿 | %.0f s (%.1f 分) |"
          % (float(np.mean(tf)) * 24, float(np.mean(tf)) * 24 / 60.0))
        a("")
        if int(args.dirs) in (48, 64, 96):
            a("前回セッションが測ったドームのみの実測（96: 306 s / 64: 238 s / "
              "48: 203 s 毎フレーム）と比較できる値である。")
    else:
        a("（新規フレームが無かったため、時間は未測定）")
    a("")
    a("## トーン（crush / blow の平均）")
    a("")
    a("| 版 | crush 平均 | blow 平均 | bound |")
    a("|---|---|---|---|")
    a("| A | %s | %s | crush<=0.0325 / blow=0 |"
      % (("%.6f" % mean_of(crush_agg["A"])) if crush_agg["A"] else "n/a",
         ("%.6f" % mean_of(blow_agg["A"])) if blow_agg["A"] else "n/a"))
    a("| B | %s | %s | 同上 |"
      % (("%.6f" % mean_of(crush_agg["B"])) if crush_agg["B"] else "n/a",
         ("%.6f" % mean_of(blow_agg["B"])) if blow_agg["B"] else "n/a"))
    a("")
    a("## ゲート")
    a("")
    if still_info.get("files"):
        a("- still-eta100.png / still-eta080.png（f%04d）は"
          " **%s**（mean|Δ|=%.3f, max=%d）"
          % (mid, "同一ではない" if not still_info["byte_identical"]
             else "同一（異常）", still_info["mean_abs_diff_uint8"],
             still_info["max_abs_diff_uint8"]))
    else:
        a("- stills: %s" % still_info.get("missing"))
    if sanity["pass"] is not None:
        a("- A 版 f0062 の C 統計 vs 前回記録 "
          "(mean %.7f / p99 %.4f / max %.4f): mean %.7f / p99 %.4f / max %.4f → **%s**"
          % (ref["mean"], ref["p99"], ref["max"],
             sanity["measured"]["mean"], sanity["measured"]["p99"],
             sanity["measured"]["max"],
             "PASS" if sanity["pass"] else "FAIL"))
    else:
        a("- A 版 f0062 の参照照合: 未測定（この実行では f0062 を描いていない）")
    a("- mp4 化は後続セッション（本スクリプトは PNG まで）")
    a("")
    a("## 成果物")
    a("")
    a("- `eta100/frame-XXXX.png`, `eta080/frame-XXXX.png`（832x480）")
    a("- `still-eta100.png`, `still-eta080.png`")
    a("- `checkpoint.json`（jobkit、フレーム単位）, `progress.log`, "
      "`frames.jsonl`, `runner.log`")
    a("- `numbers-v9samples.json`")
    a("")
    a("## 正直な記録")
    a("")
    a("- ドーム方向数は 96 ではなく **%d**（継続メモの指定）。"
      "eta 比較は測定型内で頑健だが、削減した事実はここに記録する。" % args.dirs)
    a("- 拡散マスターは毎フレーム計算（48 方向・2048² のキャッシュは無い）。")
    a("- 合成は凍結 A2.6+B3px に chunk-5 の表示正規化 E を適用"
      "（E_v3=2.2 を実効 %.6f に再正規化）。生の 2.2 を使うとフラット水面の"
      "直達レベルが v3/v4 定義より %.1f%%暗くなるため、production/chunk-5 と"
      "同じ式を採った。" % (E_v9, 100.0 * (1.0 - E_V3 / E_v9)))
    a("- A と B は同じフレームの同じ直達輸送・同じ合成・同じ E・同じ G1 を"
      "共有し、差は C_total の混合だけ。")
    a("")
    report_path = os.path.join(out, "SAMPLE-REPORT.md")
    atomic_write_bytes(report_path,
                       ("\n".join(L) + "\n").encode("utf-8"))
    log("report written: %s" % report_path)
    log("DONE wall=%.1fs frames=%d" % (time.time() - t_all, len(per_frame)))
    prog("runner_done wall=%.1fs" % (time.time() - t_all))
    logf.close()


if __name__ == "__main__":
    main()
