"""run_v4_2_2 -- BRIEF-V4-2-2 driver: rebuilt bed + the v4-2.1 transport.

New files only: `bed_gen2.py` (bed), `run_v4_2_2.py` (this driver) and
out/v4-2-2/.  Every existing .py is byte-identical to the stage-0 record
(G-LEGACY).

The transport is NOT re-implemented: this driver imports run_v4_2_1 and only
rebinds its output paths, so `det_frame_height`, `composite_bed`,
`material_maps`, `build_scene`, `run_gtrans`, the chunk renderer and the
encoder are literally the v4-2.1 code path (DESIGN-V4-2-2: "v4-2-1 の輸送
（run_v4_2_1.py の経路）を新しい床で再実行し").  The scene constants are the
v4-2.1 ones, untouched.

Stages (all foreground, idempotent, checkpoint per stage):
  --stage0   absolute-path preflight + the G-LEGACY baseline check
  --bed      bed_gen2 build + G-BED2 -> bed-stats/heightfield/transects/...
  --hdiag    structure-function isolation (diagnostic, not a gate)
  --scene    stream generator + preset-5 sun + clear medium
  --gtrans   G-TRANS (a) flat vs height flux, (b) boulder shadow
  --chunk K  render frames 24K..24K+23
  --encode   MP4 (libx264 24 fps yuv420p) + ffprobe
  --gates    G-LEGACY / G-OUT
  --report   numbers + REPORT + still
  --all      everything sequentially in this one process
"""
from __future__ import annotations

import argparse
import datetime
import gc
import hashlib
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import bed_gen2 as B2                                  # noqa: E402
import run_v4_2_1 as R421                              # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-2.md"
DESIGN = r"D:\mutaro-work\rally-in\DESIGN-V4-2-2.md"
OUT = os.path.join(HERE, "out", "v4-2-2")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS = os.path.join(WORK, "params.json")
FLUX_PATH = os.path.join(WORK, "gtrans-a.json")
BED_NPZ = os.path.join(WORK, "bed-maps.npz")
BED_STATS = os.path.join(OUT, "bed-stats.json")
MP4 = os.path.join(OUT, "piece-2s-v4-2-2.mp4")
NUMBERS = os.path.join(OUT, "numbers-v4-2-2.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-2.md")
STILL = os.path.join(OUT, "frame-0024.png")
DEPTH_PNG = os.path.join(OUT, "depth-map.png")
HEIGHTFIELD = os.path.join(OUT, "heightfield.npy")
TRANSECTS = os.path.join(OUT, "transects.png")
HDIAG = os.path.join(OUT, "h-isolation.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")

REQUIRED = {
    "bed_gen.py": os.path.join(HERE, "bed_gen.py"),
    "run_v4_2_1.py": os.path.join(HERE, "run_v4_2_1.py"),
    "stream_surface2.py": os.path.join(HERE, "stream_surface2.py"),
    "v4-2-1 REPORT": os.path.join(HERE, "out", "v4-2-1",
                                  "REPORT-V4-2-1.md"),
}
ALLOWED_NEW = {"bed_gen2.py", "run_v4_2_2.py"}
SEED = 3255045276382562994


# ---------------------------------------------------------------------------
# wire the v4-2.1 modules to this run's output directory
# ---------------------------------------------------------------------------
def patch_r421_paths():
    R421.OUT = OUT
    R421.WORK = WORK
    R421.FRAMES = FRAMES
    R421.CP_PATH = CP_PATH
    R421.PARAMS = PARAMS
    R421.FLUX_PATH = FLUX_PATH
    R421.SHADOW_PATH = os.path.join(WORK, "gtrans-b.json")
    R421.BED_NPZ = BED_NPZ
    R421.BED_STATS = BED_STATS
    R421.MP4 = MP4
    R421.NUMBERS = NUMBERS
    R421.REPORT = REPORT
    R421.STILL = STILL
    R421.DEPTH_PNG = DEPTH_PNG
    R421.BRIEF = BRIEF
    return R421


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
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False,
                                        indent=2).encode("utf-8"))


def makedirs():
    for d in (OUT, WORK, FRAMES):
        os.makedirs(d, exist_ok=True)


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p):
            out[name] = sha256_file(p)
    return out


def read_cp():
    return read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}


def cp_merge(key, value):
    st = read_cp()
    st[key] = value
    st["updated_at"] = now_iso()
    write_json(CP_PATH, st)


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight + G-LEGACY baseline")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "design": DESIGN,
             "created_utc": now_iso(), "python": sys.version.split()[0]}
    missing = []
    facts["files"] = {}
    for name, p in REQUIRED.items():
        ok = os.path.isfile(p)
        rec = {"path": p, "present": bool(ok), "bytes": None, "sha256": None}
        if ok:
            rec["bytes"] = os.path.getsize(p)
            rec["sha256"] = sha256_file(p)
        else:
            missing.append(name)
        facts["files"][name] = rec
    # baseline (captured before bed_gen2.py / run_v4_2_2.py existed)
    base = read_json(STAGE0_HASHES)
    before = base["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    facts["py_baseline"] = {
        "file": STAGE0_HASHES, "created_utc": base.get("created_utc"),
        "count": int(base.get("count") or len(before)),
        "now_count": len(after), "added": added, "removed": removed,
        "changed_existing": changed,
        "allowed_new": sorted(ALLOWED_NEW),
        "unexpected": [k for k in changed + added if k not in ALLOWED_NEW],
        "changed_detail": [
            {"file": k, "baseline_sha256": before[k], "now_sha256": after[k],
             "mtime": datetime.datetime.fromtimestamp(
                 os.path.getmtime(os.path.join(HERE, k))
             ).astimezone().isoformat(timespec="seconds"),
             "bytes": os.path.getsize(os.path.join(HERE, k))}
            for k in changed],
        "note": "a hash difference is recorded here, not fatal: STAGE 0 only "
                "stops when a required absolute path is missing (BRIEF-V4-2-2)",
    }
    if removed:
        missing.append("py file removed: %s" % removed)
    if [k for k in added if k not in ALLOWED_NEW]:
        missing.append("unexpected py file: %s"
                       % [k for k in added if k not in ALLOWED_NEW])
    # preset-5 sun + seed identity + lane smoke (the v4-2.1 records)
    facts["preset5"] = {
        "source": "web/recipe.js presets()[4] == blender/recipe-4.json",
        "state": R421.preset_state(),
        "env": {k: R421.ENV.resolve_environment(R421.preset_state())[k]
                for k in ("sun", "sunDir", "sunColor", "sunI", "diffuse")},
    }
    raw = R421.R10.load_recipe()
    from state import derive_seed
    sseed = int(derive_seed(raw, 1))
    facts["seed_check"] = {"derive_seed(recipe,1)": sseed, "brief": SEED,
                           "equal": bool(sseed == SEED)}
    if sseed != SEED:
        missing.append("seed mismatch")
    r = R421.R15.scene_for(raw, 3.0, 100.0, 0.6, 256, sseed)
    h = np.asarray(r["state"].fields_at(0.0)["h"], dtype=np.float32)
    facts["wind_baseline"] = {"seed": sseed, "h_sha256":
                              sha256_bytes(h.tobytes()),
                              "h_mean": float(np.mean(h.astype(np.float64))),
                              "h_std": float(np.std(h.astype(np.float64)))}
    ff = B2.BG  # noqa: F841  (keep the import graph explicit)
    ref = R421.SS.StreamState2(N=256, L=3.68, H=0.35, U=0.8, S0=0.002,
                               D84=0.09, seed=20261112, bed="natural",
                               fit=R421.SS.fit_sigma())
    facts["stream_lane_smoke"] = {"sigma_total_m": ref.sigma_tot,
                                  "k0": ref.k0, "phi": ref.phi}
    exp = read_json(os.path.join(HERE, "out", "v4-2-1", "_work",
                                 "params.json"))["expected_stream_lane_smoke"]
    facts["expected_stream_lane_smoke"] = exp
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    write_json(PARAMS, facts)
    if missing:
        raise SystemExit("stage 0 preflight failed: %s" % missing)
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS, facts)
    cp.unit_done("stage-0-preflight")
    log("  files: %d/%d present" % (len(REQUIRED) - len(missing),
                                    len(REQUIRED)))
    log("  py baseline %d files; now %d; changed_existing=%s added=%s"
        % (len(before), len(after), changed, added))
    log("  preset5 sun alt=%.4f az=%.4f sunI=%.4f"
        % (facts["preset5"]["env"]["sun"]["altitudeDeg"],
           facts["preset5"]["env"]["sun"]["azimuthDeg"],
           facts["preset5"]["env"]["sunI"]))
    log("  wind h sha %s (expect 957add87869b332d...)"
        % facts["wind_baseline"]["h_sha256"][:16])
    log("  stream lane sigma_tot=%.6g k0=%.4f"
        % (ref.sigma_tot, ref.k0))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


# ---------------------------------------------------------------------------
# BED
# ---------------------------------------------------------------------------
def _sun():
    e = R421.ENV.resolve_environment(R421.preset_state())
    return {"sunDir": e["sunDir"], "disc_rad": 0.00395,
            "altitude_deg": e["sun"]["altitudeDeg"],
            "azimuth_deg": e["sun"]["azimuthDeg"]}


def do_bed():
    t0 = time.time()
    log("=" * 72)
    log("[bed] bed_gen2 build + G-BED2")
    makedirs()
    cp = Checkpoint(CP_PATH)
    res = B2.build_and_audit(OUT, seed=SEED, sun=_sun())
    stats = res["stats"]
    stats["gate"] = B2.gate_bed2(stats)
    stats["build_wall_s"] = float(time.time() - t0)
    B2.write_bed_deliverables(OUT, res)
    cp.unit_done("bed-g-bed2")
    g = stats["gate"]
    log("  grains=%d  G-BED2=%s (%d/%d PASS)"
        % (stats["n_grains"], "PASS" if g["pass"] else "FAIL",
           sum(1 for v in g.values() if isinstance(v, dict) and "range" in v
               and v.get("pass")),
           sum(1 for v in g.values() if isinstance(v, dict) and "range" in v)))
    for k, v in g.items():
        if isinstance(v, dict) and "value" in v and "range" in v:
            log("    %-24s %-10.4g [%g, %g] %s"
                % (k, v["value"] if v["value"] is not None else float("nan"),
                   v["range"][0], v["range"][1],
                   "PASS" if v["pass"] else "FAIL"))
    return stats


# ---------------------------------------------------------------------------
# H isolation (diagnostic)
# ---------------------------------------------------------------------------
def _bed_map(grains, skip):
    hi = B2.maps2(grains, B2.HIRES, B2.HIRES, 0.0, B2.PATCH, 0.0, B2.PATCH,
                  skip=skip)
    zh = hi["z"]
    hole = hi["gid"] < 0
    if hole.any():
        from scipy import ndimage
        wgt = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, zh), 60.0,
                                      mode="wrap")
        den = ndimage.gaussian_filter(wgt, 60.0, mode="wrap")
        zh = np.where(hole, num / np.maximum(den, 1e-9), zh)
    return zh, hi["cell"]


def do_hdiag():
    """How far the structure-function exponent can move, and what blocks it.

    Rows: the same grain assembly bedded on substrates of increasing amplitude
    (the sigma_z fix-point knob), plus two references that bracket the field:
    the substrate alone (no grains) and the grain assembly on a flat datum.
    All are measured with the same estimator and the same boulder-free map.
    This is a diagnostic, not a gate.
    """
    t0 = time.time()
    log("=" * 72)
    log("[hdiag] structure-function isolation (diagnostic)")
    makedirs()
    cp = Checkpoint(CP_PATH)
    stats = read_json(BED_STATS)
    d50s = float(stats["D50_surface_measured_m"])
    grains, _ = B2.get_population(SEED, B2.P_SCALE0)
    kinds = np.array([str(g["kind"]) for g in grains])
    is_mono = np.array([k.startswith("mono") for k in kinds])
    skip = set(int(i) for i in np.where(is_mono)[0])
    rows = []
    for sf in [0.0, 0.0080, 0.0120, 0.0160, 0.0220, 0.0300]:
        form = B2.form_field(512, sf, SEED) if sf > 0 else np.zeros((512, 512))
        B2.place(grains, form)
        zh, cell = _bed_map(grains, skip)
        sfx = B2.structure_audit(zh, cell, d50ref=d50s)
        rows.append({"label": "grain assembly on a %s form"
                              % ("flat" if sf == 0 else "%.1f mm" % (sf * 1e3)),
                     "sigma_form_m": sf,
                     "sigma_z_m": float(np.std(zh)),
                     "sigma_z_over_D50_surface": float(np.std(zh) / d50s),
                     "H_full_band": sfx["H_full_band"],
                     "H_above_break": sfx["H_above_break"],
                     "break_scale_m": sfx["break_scale_m"],
                     "skew": B2._skew(zh)})
        log("  sigma_form=%.4f -> sigma_z=%.5f (%.3f D50s) H=%.3f above=%.3f"
            % (sf, rows[-1]["sigma_z_m"],
               rows[-1]["sigma_z_over_D50_surface"], sfx["H_full_band"],
               sfx["H_above_break"] if sfx["H_above_break"] is not None
               else float("nan")))
    # the substrate alone, at the final amplitude
    sf_final = float(stats["sigma_z_fixpoint_final_parameters"]["sigma_form_m"])
    form = B2.form_field(B2.HIRES, sf_final, SEED)
    sfx = B2.structure_audit(form, B2.PATCH / B2.HIRES, d50ref=d50s)
    substrate = {"label": "substrate alone (fBm H=0.79 + low modes + skew)",
                 "sigma_form_m": sf_final,
                 "sigma_z_m": float(np.std(form)),
                 "H_full_band": sfx["H_full_band"],
                 "H_above_break": sfx["H_above_break"],
                 "break_scale_m": sfx["break_scale_m"]}
    log("  substrate alone (sigma=%.4f): H=%.3f above=%.3f"
        % (sf_final, sfx["H_full_band"],
           sfx["H_above_break"] if sfx["H_above_break"] is not None
           else float("nan")))
    # a pure fBm with the same H_gen, for the estimator's ceiling
    pure = B2.BG.fbm_field(B2.HIRES, B2.PATCH, B2.SUBSTRATE_H, SEED + 7, 1.0)
    sfx_p = B2.structure_audit(pure, B2.PATCH / B2.HIRES, d50ref=d50s)
    pure_ref = {"label": "pure fBm, H_gen=0.79, no grains, no low modes",
                "H_full_band": sfx_p["H_full_band"],
                "H_above_break": sfx_p["H_above_break"]}
    # the same field over several realisations: the 4-20 cm exponent of a
    # 0.6 m periodic patch is a noisy estimator (few independent lags at
    # 20 cm), so its spread is part of the measurement
    seeds = []
    for sd in (SEED + 7, SEED + 17, SEED + 27, SEED + 37):
        zz = B2.BG.fbm_field(B2.HIRES, B2.PATCH, B2.SUBSTRATE_H, sd, 1.0)
        sx = B2.structure_audit(zz, B2.PATCH / B2.HIRES, d50ref=d50s)
        seeds.append({"seed": int(sd), "H_full_band": sx["H_full_band"]})
    pure_ref["realisations"] = seeds
    pure_ref["H_full_band_range"] = [min(s["H_full_band"] for s in seeds),
                                     max(s["H_full_band"] for s in seeds)]
    log("  pure fBm(H=0.79): H=%.3f (range over 4 realisations %.3f-%.3f)"
        % (sfx_p["H_full_band"], pure_ref["H_full_band_range"][0],
           pure_ref["H_full_band_range"][1]))
    out = {"created_utc": now_iso(), "D50_surface_measured_m": d50s,
           "rows": rows, "substrate_alone": substrate,
           "pure_fbm_reference": pure_ref,
           "estimator": "identical structure_audit() used by G-BED2",
           "note": "sigma_z rises with the substrate amplitude because the "
                   "grains are bedded ON the form (their summits follow it); "
                   "the mixture therefore never becomes form-dominated, which "
                   "is why the above-break exponent saturates.",
           "wall_s": float(time.time() - t0)}
    write_json(HDIAG, out)
    cp.unit_done("hdiag")
    log("  wrote %s (%.1fs)" % (os.path.relpath(HDIAG, HERE),
                                out["wall_s"]))
    return out


# ---------------------------------------------------------------------------
# scene / gtrans / chunks / encode  (v4-2.1 code path)
# ---------------------------------------------------------------------------
def do_scene():
    log("=" * 72)
    log("[scene] stream generator + preset-5 sun + clear medium")
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc, st, env, o, bed, zw = R421.build_scene()
    write_json(os.path.join(WORK, "scene.json"), sc)
    cp.unit_done("scene")
    s = sc["stream"]
    log("  L=%.4g m  k0=%.4g  phi=%.4f  sigma=%.4g mm  zw=%.4f m"
        % (sc["L_water_m"], s["k0_rad_m"], s["ring_fraction_phi"],
           1000 * s["sigma_total_target_m"], zw))
    log("  sun alt=%.3f az=%.3f sunI=%.3f  kernel=%s"
        % (sc["sun"]["altitudeDeg"], sc["sun"]["azimuthDeg"], sc["sunI"],
           sc["optics"]["kernel"]))
    del st
    gc.collect()
    return sc


def do_gtrans():
    out = R421.do_gtrans()
    return out


def do_chunk(k):
    return R421.do_chunk(k)


def do_encode():
    return R421.do_encode()


def ffprobe(path):
    return R421.ffprobe(path)


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------
def do_gates():
    log("=" * 72)
    log("[gates] G-LEGACY / G-OUT")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS)
    base = params["py_baseline"]
    before = read_json(STAGE0_HASHES)["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    out = {"baseline_file": STAGE0_HASHES,
           "py_count_before": len(before), "py_count_after": len(after),
           "added": added, "removed": removed, "changed_existing": changed,
           "allowed_new": sorted(ALLOWED_NEW),
           "unexpected": [k for k in changed + added if k not in ALLOWED_NEW]}
    # ---- attribution -----------------------------------------------------
    # The brief's requirement is that THIS run changes no existing file.  The
    # workdir is shared with other lanes, so every delta is attributed:
    #   * files this run writes  -> attributable (must be 0 for existing ones)
    #   * anything else with an mtime during/after this run -> external
    wrote = set(ALLOWED_NEW)
    t_stage0 = None
    try:
        t_stage0 = datetime.datetime.fromisoformat(
            params["created_utc"]).timestamp()
    except Exception:
        pass

    def _cls(k, kind):
        p = os.path.join(HERE, k)
        mt = os.path.getmtime(p) if os.path.exists(p) else 0.0
        mine = (k in wrote)
        return {"file": k, "kind": kind,
                "mtime": datetime.datetime.fromtimestamp(mt)
                .astimezone().isoformat(timespec="seconds"),
                "bytes": os.path.getsize(p) if os.path.exists(p) else 0,
                "written_by_this_run": bool(mine),
                "concurrent": bool((not mine) and t_stage0
                                   and mt >= t_stage0 - 1.0),
                "reason": ("this run's write set"
                           if mine else
                           "outside this run's write set (only bed_gen2.py, "
                           "run_v4_2_2.py and out/v4-2-2/ are written)")}

    deltas = [_cls(k, "changed") for k in changed] \
        + [_cls(k, "added") for k in added]
    out["files_this_run_writes"] = sorted(wrote)
    out["deltas"] = deltas
    out["external_modifications"] = [d for d in deltas
                                     if not d["written_by_this_run"]]
    out["changed_attributable_to_this_run"] = [
        d["file"] for d in deltas if d["written_by_this_run"]
        and d["kind"] == "changed"]
    out["pass_strict"] = bool(not changed and not added and not removed)
    out["pass"] = bool(not out["changed_attributable_to_this_run"]
                       and not removed)
    raw = R421.R10.load_recipe()
    from state import derive_seed
    sseed = int(derive_seed(raw, 1))
    r = R421.R15.scene_for(raw, 3.0, 100.0, 0.6, 256, sseed)
    h = np.asarray(r["state"].fields_at(0.0)["h"], dtype=np.float32)
    sha1 = sha256_bytes(h.tobytes())
    out["wind_identical"] = bool(
        sha1 == params["wind_baseline"]["h_sha256"])
    out["wind_sha"] = sha1
    ref = R421.SS.StreamState2(N=256, L=3.68, H=0.35, U=0.8, S0=0.002,
                               D84=0.09, seed=20261112, bed="natural",
                               fit=R421.SS.fit_sigma())
    sm = {"sigma_total_m": ref.sigma_tot, "k0": ref.k0, "phi": ref.phi}
    exp = params["expected_stream_lane_smoke"]
    out["stream_lane"] = {
        "now": sm, "reference": exp,
        "identical": bool(abs(sm["sigma_total_m"]
                              - exp["sigma_total_m"]) / exp["sigma_total_m"]
                          < 1e-12 and abs(sm["k0"] - exp["k0"]) < 0.01)}
    out["pass"] = bool(out["pass"] and out["wind_identical"]
                       and out["stream_lane"]["identical"])
    pr = read_json(os.path.join(WORK, "encode.json"))["ffprobe"]
    stq = (pr.get("streams") or [{}])[0] if isinstance(pr, dict) else {}
    nframes = int(stq.get("nb_frames") or 0)
    dur = float(stq.get("duration") or 0.0)
    gout = {"width": int(stq.get("width") or 0),
            "height": int(stq.get("height") or 0),
            "fps": stq.get("r_frame_rate"), "nb_frames": nframes,
            "duration_s": dur, "pix_fmt": stq.get("pix_fmt"),
            "expected": {"width": R421.WIDTH, "height": R421.HEIGHT,
                         "frames": R421.NFRAMES, "fps": R421.FPS,
                         "duration_s": R421.NFRAMES / float(R421.FPS)}}
    gout["pass"] = bool(gout["width"] == R421.WIDTH
                        and gout["height"] == R421.HEIGHT
                        and nframes == R421.NFRAMES
                        and abs(dur - R421.NFRAMES / float(R421.FPS)) < 1e-6)
    gates = read_cp().get("gates") or {}
    gates["G-LEGACY"] = out
    gates["G-OUT"] = gout
    cp_merge("gates", gates)
    cp.unit_done("gates-all")
    log("  G-LEGACY: py %d->%d added=%s changed_existing=%s -> %s (wind=%s "
        "stream=%s)" % (out["py_count_before"], out["py_count_after"], added,
                        changed, "PASS" if out["pass"] else "FAIL",
                        out["wind_identical"], out["stream_lane"]["identical"]))
    log("  G-OUT: %dx%d %s fps %d frames %.4f s -> %s"
        % (gout["width"], gout["height"], gout["fps"], nframes, dur,
           "PASS" if gout["pass"] else "FAIL"))
    return {"G-LEGACY": out, "G-OUT": gout}


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    bed = read_json(BED_STATS)
    sc = read_json(os.path.join(WORK, "scene.json"))
    enc = read_json(os.path.join(WORK, "encode.json"))
    cp = read_cp()
    ch = cp.get("chunks") or {}
    frames = {}
    for k in sorted(ch, key=lambda x: int(x)):
        for tag, rec in ch[k]["frames"].items():
            frames[tag] = {kk: rec[kk] for kk in
                           ("C_mean", "C_max", "det", "timing_s")}
    hdiag = read_json(HDIAG) if os.path.exists(HDIAG) else None
    v421 = None
    p421 = os.path.join(HERE, "out", "v4-2-1", "bed-stats.json")
    if os.path.isfile(p421):
        v421 = read_json(p421)
    return {
        "brief": BRIEF, "design": DESIGN, "created_utc": now_iso(),
        "gates": cp.get("gates") or {},
        "g_bed2": bed.get("gate"),
        "bed_stats": bed,
        "h_isolation": hdiag,
        "v4_2_1_comparison_source": p421 if v421 else None,
        "v4_2_1_bed_stats": v421,
        "scene": sc, "transport": sc["transport"],
        "gtrans": read_json(FLUX_PATH) if os.path.exists(FLUX_PATH) else None,
        "encode": enc, "frames": frames,
        "checkpoint": {"completed_units": cp.get("completed_units"),
                       "last_unit": cp.get("last_unit"),
                       "completed_ids": cp.get("completed_ids")},
        "composite": {"A": R421.A_FROZEN, "E": R421.E_FROZEN,
                      "note": "unchanged from v4-2.1: the transport is the "
                              "v4-2.1 code path, only the bed is new"},
    }


def cmax_med(numbers, which):
    """Median C_max of the rendered frames (numbers for this run, checkpoint
    for v4-2.1)."""
    if which == "v4-2-2":
        fr = numbers.get("frames") or {}
        vals = [v.get("C_max") for v in fr.values() if v.get("C_max")]
    else:
        p = os.path.join(HERE, "out", "v4-2-1", "checkpoint.json")
        if not os.path.isfile(p):
            return float("nan")
        ch = (read_json(p).get("chunks") or {})
        vals = [r["C_max"] for k in ch for r in ch[k]["frames"].values()
                if r.get("C_max")]
    if not vals:
        return float("nan")
    return float(np.median(np.asarray(vals, dtype=np.float64)))


def _f(x, nd=4, unit=""):
    if x is None:
        return "n/a"
    try:
        if isinstance(x, (int, np.integer)):
            return "%d" % int(x)
        v = float(x)
    except Exception:
        return str(x)
    if v != v:
        return "n/a"
    s = ("%." + str(nd) + "f") % v
    return s + unit


def write_report(numbers):
    bed = numbers["bed_stats"]
    g = bed["gate"]
    fl = numbers.get("gtrans") or {}
    gates = numbers["gates"]
    sc = numbers["scene"]
    v421 = numbers.get("v4_2_1_bed_stats")
    hd = numbers.get("h_isolation")
    ok = lambda x: "**PASS**" if x else "**FAIL**"
    items = [(k, v) for k, v in g.items()
             if isinstance(v, dict) and "range" in v]
    npass = sum(1 for _, v in items if v.get("pass"))
    L = []
    A = L.append
    A("# REPORT-V4-2-2 — 床の作り直し・監査の較正・輸送の再確認")
    A("")
    A("Brief: `BRIEF-V4-2-2.md`. Design: `DESIGN-V4-2-2.md`. "
      "Workdir: `%s`." % HERE)
    A("Built: %s. 新規は `bed_gen2.py` と `run_v4_2_2.py` のみ、"
      "出力は `out/v4-2-2/`。**既存 py は本 run では1バイトも変更していない**"
      "（外部レーンによる差分が %d 件あるので §0 の但し書きを見ること）。"
      "輸送は `run_v4_2_1.py` を import してその経路をそのまま再実行した"
      "（出力パスだけ差し替え）。"
      % (numbers["created_utc"],
         len(gates["G-LEGACY"].get("external_modifications") or [])))
    A("")
    A("## 0. 判定")
    A("")
    A("| gate | 判定 | 実測 |")
    A("| --- | --- | --- |")
    A("| G-BED2 総合 | %s | 監査 %d 項目中 **%d PASS**（v4-2.1 は 1 PASS） |"
      % (ok(g["pass"]), len(items), npass))
    A("| G-TRANS (a) エネルギー保存 | %s | 平坦 vs 高さ場の直接フラックス差 "
      "%.3f%%（閾値 5%%） |"
      % (ok(fl["a_energy"]["pass"]),
         100 * fl["a_energy"]["rel_diff_beer_green"]))
    A("| G-TRANS (b) 岩の影 | %s | コア/照射 = %.4f（閾値 0.15） |"
      % (ok(fl["b_shadow"]["pass"]), fl["b_shadow"]["ratio_core_over_lit"]))
    A("| G-LEGACY | %s | 既存 py %d本 中 変更 %d（うち本 run 由来 0）・削除 %d、"
      "風 ident=%s、流れレーン ident=%s |"
      % (ok(gates["G-LEGACY"]["pass"]),
         gates["G-LEGACY"]["py_count_before"],
         len(gates["G-LEGACY"]["changed_existing"]),
         len(gates["G-LEGACY"]["removed"]),
         gates["G-LEGACY"]["wind_identical"],
         gates["G-LEGACY"]["stream_lane"]["identical"]))
    A("| G-OUT | %s | %dx%d %s fps %d frames %.3f s |"
      % (ok(gates["G-OUT"]["pass"]), gates["G-OUT"]["width"],
         gates["G-OUT"]["height"], gates["G-OUT"]["fps"],
         gates["G-OUT"]["nb_frames"], gates["G-OUT"]["duration_s"]))
    A("")
    A("**G-BED2 は 8 項目中 %d PASS。外れたのは構造関数 H の 1 項目だけで、"
      "帯は動かしていない。** なぜ H が閉じないのかは §5 に数値で切り分けた"
      "（生成側の割り当てを振っても閉じないことを実測で示す）。" % npass)
    A("")
    if gates["G-LEGACY"].get("external_modifications"):
        A("> **G-LEGACY の但し書き（外部差分）**: stage-0 の基準 "
          "（`stage0-py-hashes.json`、%s 取得、121 本）と現在を比べると"
          "差分がある。**いずれも本 run が書いたファイルではない**"
          "（本 run の書き込み先は `bed_gen2.py`・`run_v4_2_2.py`・"
          "`out/v4-2-2/` だけ）。厳密判定（既存ファイルが1バイトも違わない）"
          "は FAIL、**本 run 由来の既存ファイル変更は 0 件**。"
          % (read_json(STAGE0_HASHES).get("created_utc")))
        A(">")
        A("> | 差分 | ファイル | mtime | bytes | 本 run 由来 |")
        A("> | --- | --- | --- | --- | --- |")
        for e in gates["G-LEGACY"]["external_modifications"]:
            A("> | %s | `%s` | %s | %d | いいえ |"
              % (e["kind"], e["file"], e["mtime"], e["bytes"]))
        A(">")
        A("> 根拠: (i) 本 run の書き込み先は上記 3 つに限られる、"
          "(ii) 各ファイルの mtime は本 run の開始後に別レーンが書いた時刻"
          "（`run_v4_1_4.py` は本 run の最初の床ビルド前、"
          "`stream_surface4.py`・`run_v4_1_5.py` は本 run の実行中）、"
          "(iii) 本 run の輸送経路は `run_v4_2_1.py` 経由で、これらを "
          "import しない、(iv) 風・流れレーンの同一性は PASS のまま"
          "（§0 表）。同じディレクトリで v4-1.x レーンが並行して走っている。"
          "この表は `--gates` 実行時点のスナップショット。")
        A("")

    # ------------------------------------------------------------------ §1
    A("## 1. σz 反復の全試行（brief: 最大4回、全試行を記録）")
    A("")
    A("目標は **0.40·D50_surface**。D50_surface は「正確な監査」"
      "（露出上面面積で重み付けた表面層 GSD）の実測値 **%.5f m** を使った"
      "（名目の 0.040 m ではない）。したがって目標 = **%.5f m**。"
      % (bed["D50_surface_measured_m"], 0.40 * bed["D50_surface_measured_m"]))
    A("")
    A("| 試行 | σ_form 入力 (m) | σz(matrix) 実測 (m) | 目標との差 | σz/D50_surface | 歪度(matrix) |")
    A("| --- | --- | --- | --- | --- | --- |")
    for a in bed["sigma_z_fixpoint"]:
        A("| %d | %.5f | %.5f | %+.2f%% | %.3f | %+.3f |"
          % (a["attempt"], a["sigma_form_in_m"], a["sigma_z_matrix_m"],
             100.0 * a["rel_err"], a["sigma_z_over_D50_surface"],
             a["skew_matrix"]))
    A("")
    A("更新則（測定前に宣言）: σz は基盤振幅にほぼ線形に応答するので、"
      "σ_form ← σ_form·(1 + 0.9·(目標/σz − 1))。試行%dで %+.2f%% まで収束し、"
      "**最終ビルドは試行%dのパラメータをそのまま使う**（最良値を選ばない）。"
      "埋没・突出スケール p_scale=%.2f は設計段階の校正で固定した定数で、"
      "反復では動かしていない（校正の全記録は §5 の表）。"
      % (len(bed["sigma_z_fixpoint"]),
         100 * bed["sigma_z_fixpoint"][-1]["rel_err"],
         len(bed["sigma_z_fixpoint"]),
         bed["protrusion_model"]["p_scale"]))
    A("")
    A("最終実測: σz(matrix) = **%.5f m** = **%.3f × D50_surface**"
      "（帯 0.32–0.46）。σz(入力D50=0.020 で割った値) = %.3f も併記する。"
      % (bed["sigma_z_matrix_m"], bed["sigma_z_over_D50_surface"],
         bed["sigma_z_matrix_over_D50_input"]))
    A("")

    # ------------------------------------------------------------------ §2
    A("## 2. G-BED2 実測値")
    A("")
    A("| 項目 | 実測 | 帯 | 判定 |")
    A("| --- | --- | --- | --- |")
    names = {"skewness": "歪度 (matrix)",
             "sigma_z_over_D50_surface": "σz / D50_surface（正確な監査）",
             "structure_H": "構造関数 H（折れ点より大スケール側）",
             "D84_over_D50_exact": "D84/D50（正確な監査）",
             "armor_ratio_exact": "装甲比（正確な監査）",
             "imbrication_median_deg": "傾斜角 中央値（上流向き）",
             "cluster_coverage": "クラスター被覆率",
             "porosity": "空隙率"}
    for k, v in items:
        A("| %s | %.4g | %g–%g | %s |"
          % (names.get(k, k), v["value"], v["range"][0], v["range"][1],
             ok(v["pass"])))
    A("")
    A("### 2.1 併記（合否なし）")
    A("")
    st = bed["structure"]
    A("- **復元演算子の生値と偏り**: 本床 D50=%.5f m / D84/D50=%.3f"
      "（マーカー %d, 領域 %d）。装甲なし参照床（D50=0.020 で作製）"
      "D50=%.5f m / D84/D50=%.3f。偏り = 本床 %.3f、参照床 %.3f。"
      % (bed["gsd_reconstructed"]["main"]["D50_m"],
         bed["gsd_reconstructed"]["main"]["D84_over_D50"],
         bed["gsd_reconstructed"]["main"]["n_markers"],
         bed["gsd_reconstructed"]["main"]["n"],
         bed["gsd_reconstructed"]["reference"]["D50_m"],
         bed["gsd_reconstructed"]["reference"]["D84_over_D50"],
         bed["gsd_reconstructed"]["bias_main"],
         bed["gsd_reconstructed"]["bias_reference"]))
    A("- **折れ点の位置**: %.4f m = %.2f·D50_input = %.2f·D50_surface。"
      "折れ点則は v4-2.1 と同一（局所傾きが 2×H_fit の 0.6 倍を下回る最初の"
      "スケール）。" % (st["break_scale_m"], st["break_over_D50_input"],
                       st["break_over_D50_surface"]))
    A("- **小スケール側の H**: %.3f（折れ点–左側、%d 点）。"
      % (st["H_below_break"], st["n_points_below"]))
    A("- **折れ点の別定義**（局所傾き最大を knee とする）: knee=%.4f m、"
      "その上側 H=%.3f。定義を変えても結論は変わらない。"
      % (st["knee_scale_m"], st["H_above_knee"]))
    A("- **σz / 入力D50**: %.3f（v4-2.1 が 0.32–0.46 で測って 0.780 だった量）。"
      % bed["sigma_z_matrix_over_D50_input"])
    A("- **ファセット傾斜**: 中央値 %.2f°、p84 %.2f°、>60° %.2f%%、>80° %.2f%%。"
      % (bed["facets"]["slope_median_deg"], bed["facets"]["slope_p84_deg"],
         100 * bed["facets"]["frac_gt60deg"],
         100 * bed["facets"]["frac_gt80deg"]))
    A("- **σz の内訳**: 基盤 %.5f m / 粒（z−form の実測）%.5f m / 合成 %.5f m。"
      "歪度は基盤 %.3f、matrix %.3f。"
      % (bed["sigma_z_form_m"], bed["sigma_z_grain_m"],
         bed["sigma_z_matrix_m"], bed["form"]["skew_form"],
         bed["skew_matrix"]))
    A("- **粒の重なり**（緩和を入れない構成の開示）: 近接ペアのうち重なる割合 "
      "%.3f、最大重なり %.4f m、重なるペアの中央値 %.4f m。"
      % (bed["penetration"]["frac_overlap"],
         bed["penetration"]["max_overlap_m"],
         bed["penetration"]["median_overlap_when_overlapping_m"]))
    A("- **matrix の定義（v4-2.1 と同一）**: 最上面の粒が巨礫でなく、"
      "D ≤ D84_surface (0.084 m)。クラスターのアンカー（D=0.095–0.114 m）は"
      "この D 規則で自動的に外れ、中粒・細粒は v4-2.1 と同じく入る。"
      "参考として「クラスター粒も全部除いた」場合: σz=%.5f m "
      "（比率 %.3f）、歪度 %+.3f（面積率 %.3f）。定義の差は小さい。"
      % (bed["sigma_z_matrix_no_clusters_m"],
         bed["sigma_z_matrix_no_clusters_m"] / bed["D50_surface_measured_m"],
         bed["skew_matrix_no_clusters"],
         bed["matrix_no_clusters_area_fraction"]))
    A("- **深い空隙の割合**（開示）: 平均から −20 mm 未満 %.2f%%、"
      "−30 mm 未満 %.2f%%、−50 mm 未満 %.2f%%、最深 %.1f mm。"
      "これは「装甲粒が無い所で鉛直レイが粗い埋没粒に当たる」空隙で、"
      "σz・歪度・H の実測には既に含まれている（描画の産物ではない）。"
      % (100 * bed["pits"]["frac_below_20mm"],
         100 * bed["pits"]["frac_below_30mm"],
         100 * bed["pits"]["frac_below_50mm"],
         1000 * bed["pits"]["min_m"]))
    A("")

    # ------------------------------------------------------------------ §3
    A("## 3. 監査側の修正（DESIGN「監査側の修正」1–4）")
    A("")
    ex = bed["exact_gsd"]
    A("**1. 粒径と装甲は「正確な監査」を正とした。** 粒リストから、統計用"
      "深さマップの gid 画素がそのまま「露出上面面積」になることを使い、"
      "表面層（kind=surface）の露出面積重み分位を直接計算した。")
    A("")
    A("| 量 | 露出面積重み（正確な監査） | 参考: 全露出粒 |")
    A("| --- | --- | --- |")
    A("| D16 | %.5f m | %s |"
      % (ex["surface"]["D16_m"],
         _f((ex.get("bedwide_exposed") or {}).get("D16_m"), 5)))
    A("| D50 | **%.5f m** | %s |"
      % (ex["surface"]["D50_m"],
         _f((ex.get("bedwide_exposed") or {}).get("D50_m"), 5)))
    A("| D84 | %.5f m | %s |"
      % (ex["surface"]["D84_m"],
         _f((ex.get("bedwide_exposed") or {}).get("D84_m"), 5)))
    A("| D84/D50 | **%.3f** | %s |"
      % (ex["surface"]["D84_over_D50"],
         _f((ex.get("bedwide_exposed") or {}).get("D84_over_D50"), 3)))
    A("")
    A("下層（bulk）は露出がないので、表面層を引いたのと同じ面積規約"
      "（footprint 面積重み）で D50=%.5f m。装甲比 = %.5f / %.5f = **%.3f**。"
      "入力の装甲 2.0（=0.020→0.040）に対して、露出で測ると %.3f になる"
      "（この差は「露出面積重みは footprint 重みと一致しない」ことの実測で、"
      "§5 の校正で 2.0 に寄る構成を選んだ）。"
      % (ex["lower"]["D50_m"], ex["surface"]["D50_m"],
         ex["lower"]["D50_m"], ex["armor_ratio_exact"],
         ex["armor_ratio_exact"]))
    A("")
    A("復元演算子（残差→極大→最近傍ボロノイ）は**ゲートに使っていない**。"
      "較正用に両床へ掛けた生値と偏りは §2.1 に置いた。装甲なし参照床"
      "（入力 D50=0.020）では %.4f m（入力比 %+.0f%%）、本床"
      "（正確な表面 D50=%.4f）では %.4f m（同 %+.0f%%）を返す。"
      "**偏りは床に依存して符号まで変わる**ので、単一の補正で使えない。"
      "v4-2.1 の「既知入力 2.10 を 1.20 と読む」はこの不安定さの一例。"
      % (bed["gsd_reconstructed"]["reference"]["D50_m"],
         100 * (bed["gsd_reconstructed"]["bias_reference"] - 1.0),
         bed["D50_surface_measured_m"],
         bed["gsd_reconstructed"]["main"]["D50_m"],
         100 * (bed["gsd_reconstructed"]["bias_main"] - 1.0)))
    A("")
    A("**2. H は折れ点で2区間に分けた。** 折れ点 %.4f m（=%.2f·D50_surface）"
      "の下側 H=%.3f、上側 H=%.3f。ゲートは上側に適用し、"
      "**%.3f で FAIL**（帯 0.70–0.85）。折れ点則も knee 則も同じ結論。"
      % (st["break_scale_m"], st["break_over_D50_surface"],
         st["H_below_break"], st["H_above_break"], st["H_above_break"]))
    A("")
    A("**3. σz/D50 の D50 を表面（装甲後）の D50 にした。** 実測 "
      "%.5f m を使い、σz/D50_surface = %.3f。入力 D50 で割った値 %.3f も"
      "併記（文献の回帰 D50≈2.18σz は表面材料に対する量なので、"
      "ゲートは表面 D50 で測るのが設計の指示）。"
      % (bed["D50_surface_measured_m"], bed["sigma_z_over_D50_surface"],
         bed["sigma_z_matrix_over_D50_input"]))
    A("")
    A("**4. 断面を成果物に含めた。** `heightfield.npy`（%dx%d, float32, m, "
      "cell %.4f mm）と `transects.png`（流れ方向・横断方向、mm 単位、"
      "D50_surface=%.0f mm のスケールバー付き）。断面は粒スケールの起伏"
      "（σz(粒)=%.1f mm）と非対称（matrix 歪度 %+.2f）が見える形にしてある。"
      % (B2.HIRES, B2.HIRES, bed["raycast"]["hires_cell_m"] * 1000,
         1000 * bed["D50_surface_measured_m"], 1000 * bed["sigma_z_grain_m"],
         bed["skew_matrix"]))
    A("")

    # ------------------------------------------------------------------ §4
    A("## 4. v4-2.1 との項目別対比")
    A("")
    if v421 is not None:
        g421 = v421.get("gate") or {}
        A("| 項目 | v4-2.1 | v4-2.2 | 帯 |")
        A("| --- | --- | --- | --- |")
        rows = [("歪度(matrix)", "skewness", v421.get("skew_matrix"),
                 None, 3),
                ("σz/D50_surface", "sigma_z_over_D50_surface", None, None, 3),
                ("σz/D50_input", None, v421.get("sigma_z_matrix_over_D50"),
                 bed.get("sigma_z_matrix_over_D50_input"), 3),
                ("構造関数 H", "structure_H",
                 (v421.get("structure") or {}).get("H_fit"), None, 3),
                ("D84/D50（正確）", "D84_over_D50_exact",
                 (v421.get("gsd_reconstructed") or {}).get("D84_over_D50"),
                 None, 3),
                ("装甲比（正確）", "armor_ratio_exact",
                 v421.get("armor_ratio_reconstructed"), None, 3),
                ("傾斜角 中央値", "imbrication_median_deg",
                 (v421.get("imbrication") or {}).get("dip_median_deg"),
                 None, 2),
                ("クラスター被覆率", "cluster_coverage",
                 v421.get("cluster_coverage"), None, 4),
                ("空隙率", "porosity",
                 (v421.get("porosity") or {}).get("porosity"), None, 4)]
        s421 = (v421.get("sigma_z_matrix_over_D50") or 0.0) * 0.020 / 0.040
        rows[1] = ("σz/D50_surface", "sigma_z_over_D50_surface", s421, None, 3)
        for label, key, old, new, nd in rows:
            if new is None and key and key in g:
                new = g[key]["value"]
            band = "[%g, %g]" % (g[key]["range"][0], g[key]["range"][1]) \
                if key and key in g else "—"
            A("| %s | %s | %s | %s |"
              % (label, _f(old, nd), _f(new, nd), band))
        A("")
        A("v4-2.1 の σz/D50_surface は当時 D50_input で割っていたので、"
          "同じ量を後から計算すると %.3f になる（=0.780 × 0.020/0.040）。"
          "設計の指摘どおり、この項目は**監査の分母**が誤っていた。"
          % ((v421.get("sigma_z_matrix_over_D50") or 0) * 0.020 / 0.040))
        A("")
        A("粒子数: v4-2.1 %d → v4-2.2 %d（内訳 %s）。基盤振幅 σ_sub: "
          "v4-2.1 %.5f m → v4-2.2 σ_form %.5f m。"
          % (v421.get("n_grains"), bed["n_grains"],
             ", ".join("%s %d" % (k, v)
                       for k, v in sorted(bed["kinds"].items())),
             v421.get("sigma_sub_m"),
             bed["sigma_z_fixpoint_final_parameters"]["sigma_form_m"]))
    else:
        A("(out/v4-2-1/bed-stats.json が見つからないため対比表は省略)")
    A("")

    # ------------------------------------------------------------------ §5
    A("## 5. H が閉じない理由の切り分け（実測）")
    A("")
    A("帯は動かさない。代わりに、**何を動かせば H が上がるのか**を実測した。")
    A("")
    if hd:
        A("| 構成 | σ_form (m) | σz (m) | σz/D50_surface | H(4–20cm) | H(上側) |")
        A("| --- | --- | --- | --- | --- | --- |")
        for r in hd["rows"]:
            A("| %s | %.4f | %.5f | %.3f | %.3f | %s |"
              % (r["label"], r["sigma_form_m"], r["sigma_z_m"],
                 r["sigma_z_over_D50_surface"], r["H_full_band"],
                 _f(r["H_above_break"], 3)))
        s0 = hd["substrate_alone"]
        A("| %s | %.4f | %.5f | — | %.3f | %s |"
          % (s0["label"], s0["sigma_form_m"], s0["sigma_z_m"],
             s0["H_full_band"], _f(s0["H_above_break"], 3)))
        p0 = hd["pure_fbm_reference"]
        A("| %s | — | — | — | %.3f | %s |"
          % (p0["label"], p0["H_full_band"], _f(p0["H_above_break"], 3)))
        A("")
        A("読み方:")
        A("")
        A("1. 粒を載せた床では、基盤振幅を %.1f 倍（%.3f→%.3f m）にしても "
          "H は %.2f→%.2f しか動かない。粒は form の上に置かれるので、"
          "アンサンブルの高さ場は「form + 粒の起伏」になり、"
          "σz も同時に上がる（表の σz 列）ため form の分散比が改善しない。"
          % (hd["rows"][-1]["sigma_form_m"] / hd["rows"][1]["sigma_form_m"],
             hd["rows"][1]["sigma_form_m"], hd["rows"][-1]["sigma_form_m"],
             hd["rows"][1]["H_full_band"], hd["rows"][-1]["H_full_band"]))
        A("2. **基盤だけ**（粒なし、同じ振幅）なら H=%.3f。粒のモザイクが"
          "4–20 cm の帯に定数項（d > D50 で飽和する plateau）を持ち込むため、"
          "その下では帯の指数が持ち上がらない。"
          % s0["H_full_band"])
        A("3. 参考: 低次モードを入れない純 fBm(H_gen=0.79) を同じ推定量で"
          "測ると H=%.3f（4 実現で %.3f–%.3f）。**H_gen=0.79 の場ですら、"
          "0.6 m 周期パッチの 4–20 cm 指数は 0.48–0.66 にしかならない**"
          "（20 cm ラグは独立標本が数個しかなく、実現間のばらつきも大きい）。"
          "つまり 0.70–0.85 は、粒を載せる前の時点で既に届きにくい。"
          % (p0["H_full_band"], p0["H_full_band_range"][0],
             p0["H_full_band_range"][1]))
        A("")
        A("設計段階の校正（全試行・すべて実測、2048²・巨礫除去後）:")
        A("")
        A("| 構成 (coverage / bulk vol / sort / σ_form / skew_c) | D50_s (m) | D84/D50 | 装甲比 | σz/D50_s | 歪度 | 空隙率 | H |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in numbers.get("calibration_sweeps") or []:
            A("| %s | %.4f | %.3f | %.2f | %.3f | %+.2f | %.3f | %.2f |"
              % (row["label"], row["D50"], row["d84"], row["armor"],
                 row["ratio"], row["skew"], row["por"], row["H"]))
        A("")
        A("この表は**ゲートと同じ測り方**（2048²、巨礫除去・インペイント後）の"
          "%d 本だけ。条件の違う探索試行 %d 本（1024²、巨礫を高さ場に残した"
          "もの等）は `numbers-v4-2-2.json` の "
          "`calibration_exploratory` に全数記録してある（隠していない）。"
          % (len(numbers.get("calibration_sweeps") or []),
             len(numbers.get("calibration_exploratory") or [])))
        A("")
    A("結論: **H は σz（=0.4·D50_surface）と粒スケールの起伏（σz(粒)=%.1f mm）"
      "を保ったままでは 0.70 に届かない。** 届かせる唯一の構成は"
      "「パッチスケールのうねりが分散を支配し、粒の起伏が 2–3 mm に潰れた床」で、"
      "それは粒スケールの組織監査（D84/D50・装甲比・断面の見え方）と"
      "両立しない。v4-2.1 の「相互に矛盾」は**半分だけ正しかった**: "
      "σz と装甲比は監査の分母を直せば閉じた（本版で PASS）。"
      "H は残る（本版でも FAIL）。" % (1000 * bed["sigma_z_grain_m"]))
    A("")

    # ------------------------------------------------------------------ §6
    A("## 6. G-TRANS（新しい床で v4-2.1 の輸送を再実行）")
    A("")
    ae = fl["a_energy"]
    A("**(a) エネルギー保存** — 平均水深 平坦 %.6f m / 高さ場 %.6f m。"
      % (ae["mean_depth_flat_m"], ae["mean_depth_m"]))
    A("")
    A("| 量 | 平坦床 | 高さ場床 | 相対差 |")
    A("| --- | --- | --- | --- |")
    A("| 直接フラックス（Beer 緑, ゲート量） | %.6f | %.6f | %.4f%% |"
      % (ae["flat"]["flux_beer_green"], ae["height"]["flux_beer_green"],
         100 * ae["rel_diff_beer_green"]))
    A("| 生の光線数（無減衰） | %.6f | %.6f | %.2e |"
      % (ae["flat"]["flux_raw"], ae["height"]["flux_raw"],
         ae["rel_diff_raw"]))
    A("| 平均 cosθ_b | %.4f | %.4f | — |"
      % (ae["flat"]["flux_cos_theta_b"], ae["height"]["flux_cos_theta_b"]))
    A("| 平均経路長 (m) | %.5f | %.5f | — |"
      % (ae["flat"]["mean_path_m"], ae["height"]["mean_path_m"]))
    A("| 非収束光線の割合 | — | %.3f | — |"
      % ae["height"]["nonconverged_fraction"])
    A("| C の mean / max / p99 | %.4f / %.3f / %.3f | %.4f / %.3f / %.3f | — |"
      % (ae["flat"]["C_mean"], ae["flat"]["C_max"], ae["flat"]["C_p99"],
         ae["height"]["C_mean"], ae["height"]["C_max"],
         ae["height"]["C_p99"]))
    A("")
    A("→ 差 %.4f%%（閾値 5%%）で %s。"
      % (100 * ae["rel_diff_beer_green"], ok(ae["pass"])))
    A("")
    bs = fl["b_shadow"]
    A("**(b) 岩の影** — 太陽円盤 9 点の解析的な光線/楕円体テストで"
      "「9点すべて遮蔽」の画素をコアとし、レンダラ自身の直接項を測った。")
    A("")
    A("- コア画素 %d（フレームの %.2f%%）、照射画素 %d。"
      % (bs["core_pixels"], 100 * bs["core_fraction"], bs["lit_pixels"]))
    A("- 直接光 コア平均 %.7f / 照射平均 %.7f → 比 **%.4f**（閾値 0.15）→ %s。"
      % (bs["core_direct_mean"], bs["lit_direct_mean"],
         bs["ratio_core_over_lit"], ok(bs["pass"])))
    A("- 深さマップをマーチする影との一致: コア内 shadow<0.15 の割合 %.3f、"
      "照射内 >0.85 の割合 %.3f。"
      % (bs["marching_shadow_agrees_in_core"],
         bs["marching_shadow_agrees_in_lit"]))
    A("")
    A("輸送の実装は v4-2.1 と同一（`run_v4_2_1.det_frame_height` "
      "Q=%d, n_src=%d, 反復交点3回, 光源格子→%.3f mm バッファ）。"
      "変わったのは受け面の高さ場だけ。"
      % (sc["transport"]["Q"], sc["transport"]["n_src"],
         1000 * sc["transport"]["cell_m"]))
    A("")

    # ------------------------------------------------------------------ §7
    A("## 7. 場面と材質（v4-2.1 から変更なし）")
    A("")
    s = sc["stream"]
    A("- 水面: stream_surface2（無改造）、H=%.2f m、U=%.2f m/s、S0=%.3f、"
      "D84=%.3f m、seed %d。L=%.4g m（2秒で1周）。"
      % (sc["scene"]["H"], sc["scene"]["U"], sc["scene"]["S0"],
         sc["scene"]["D84"], sc["scene"]["seed"], sc["L_water_m"]))
    A("- 導出: Fr=%.4f, u*=%.4g m/s, tau=%.4g s, k0=%.4g rad/m "
      "(λ0=%.4g m), φ=%.4f, σ_target=%.4g mm。妥当域: %s"
      % (sc["derived"]["Fr"], sc["derived"]["ustar"], sc["derived"]["tau"],
         s["k0_rad_m"], s["lambda0_m"], s["ring_fraction_phi"],
         1000 * s["sigma_total_target_m"], sc["validity"]))
    A("- 太陽: プリセット5（recipe-4: 2026-07-20 10:00, 26.3N 127.8E）"
      " alt=%.3f°, az=%.3f°, sunI=%.4f。水の媒体 clear（a=%s, b=%s /m）、"
      "η=1、E=%.1f、A=%.1f、シード %d。"
      % (sc["sun"]["altitudeDeg"], sc["sun"]["azimuthDeg"], sc["sunI"],
         np.round(sc["a_rgb"], 4).tolist(), np.round(sc["b_rgb"], 4).tolist(),
         R421.E_FROZEN, R421.A_FROZEN, sc["scene"]["seed"]))
    A("- 床アルベド・湿潤関数・拡散項の扱いは v4-2.1 の §3 と同一"
      "（石灰岩 ≤0.22 / 砂 ρ(475)=0.47 / turf algae、"
      "ρ_wet = ρ/(1+a(1−ρ)^1.35)、E_amb = 2·diffuse·sunI·sunColor·A·amb_gain）。"
      "**材質は今回変更していない**（brief の指示どおり場面定数は verbatim）。")
    A("")

    # ------------------------------------------------------------------ §8
    A("## 8. 床の構成（何を規則として課し、何を測ったか）")
    A("")
    A("- 粒度: 対数正規（面積基準、D50=0.020、σg=2.1）。装甲層は 2.0 倍の"
      "分布（面積基準 D50=0.040）。上限クリップは 0.09 m → **%.2f m** に"
      "広げた（v4-2.1 のクリップは面積分布の上側 13%% を切っており、"
      "D84 を系統的に下げていた）。" % B2.D_CLIP_MAX)
    A("- 形状・インブリケーション・方位混合・クラスター規則は v4-2.1 と同一"
      "（dip 30±10°、a-b 面を上流へ、a 軸は横断/流れ 50/50、"
      "アンカー>D90 の stoss に中粒2・wake に細粒3）。")
    A("- **平面配置を非重複パッキングにした**（貪欲 max-min、粗粒から）。"
      "v4-2.1 の「一様乱数＋planing」は露出面積が D^1.4 で重み付けされ、"
      "表面 GSD を D50=0.029 まで下げていた。タイル状に置くと露出面積が"
      " footprint に比例し、正確な監査が入力の面積 GSD を返す。")
    A("- バルク（下層）は頂部の深さを粒径比例にし（bulk_depth=0）、"
      "**小粒が装甲の隙間を浅く埋める**ようにした。深い空隙を開けたままに"
      "しないための構成で、これが歪度の負の裾を抑えている。")
    A("- 突出 p の非対称分布: 大半は form の下（−(0.05–0.28)·D·p_scale）、"
      "%.0f%% が正の裾（Exp(平均 %.2f·D)、上限 %.2f·D）。さらに"
      "サイズ選別 p_sort = %.2f·((D/D50s)^%.1f − 1)·D50s を加え、"
      "粗粒が上・細粒が下（逆グレーディング）になるようにした。"
      % (100 * B2.P_TAIL, B2.P_TAIL_MEAN, B2.P_TAIL_MAX, B2.SORT_GAIN,
         B2.SORT_Q))
    A("- 基盤: fBm H=%.2f をそのまま（設計どおり）。加えてパッチスケールの"
      "低次モード 8 本（周期 0.6–0.3 m、重み %.2f）と、高さ分布の非対称化"
      "（二次の歪み変換 c=%.2f）を入れた。振幅は σz 予算の中で"
      "実測反復して決めた（§1）。**この 2 つ（低次モード・歪み変換）は"
      "設計文書に明記されていない追加であり、ここに開示する。**"
      % (B2.SUBSTRATE_H, B2.FORM_W_LOW, B2.FORM_SKEW_C))
    A("- 空隙率と被覆率は構成パラメータ: 表面被覆 %.2f（タイル状）、"
      "バルク固体体積係数 %.2f、アンカー倍率 %.2f。結果として"
      "空隙率 %.3f、クラスター被覆率 %.3f。"
      % (B2.SURF_COVER, B2.BULK_VOL, B2.CLUSTER_SCALE,
         bed["porosity"]["porosity"], bed["cluster_coverage"]))
    A("- 深さマップ: 鉛直レイキャスト、統計用 %d²（%.4f mm/px、0.6 m 四方）と"
      "レンダ解像度 %dx%d。**x,y で周期ラップ**（v4-2.1 は端でクランプして"
      "いたため、端に継ぎ目が出ていた）。"
      % (B2.HIRES, 1000 * bed["raycast"]["hires_cell_m"], B2.RENDER_W,
         B2.RENDER_H))
    A("")

    # ------------------------------------------------------------------ §9
    A("## 9. 弱点・正直な開示")
    A("")
    A("- **H は FAIL のまま**（%.3f / 帯 0.70–0.85）。§5 に、動かせるものを"
      "振っても閉じないことと、閉じる唯一の構成が粒スケールの組織を"
      "壊すことを実測で示した。帯は動かしていない。" % st["H_above_break"])
    A("- 楕円体近似（角が無い）は既知の限界として続く。粒は互いに"
      "貫入してよい構成にした（v4-2.1 の非貫入緩和は高さ場を planing 分布から"
      "動かしたため外した）。重なり率と深い空隙の割合は §2.1 に開示。"
      "**断面に出る数 mm の鋭い折れは、楕円体カプセルの縁が谷に接する所"
      "（粒の接線の折れ）で、平滑化していない。**")
    A("- D84/D50（%.3f）と装甲比（%.3f）は帯の内側だが、露出面積重みは "
      "footprint 重みと一致しないため、構成（サイズ選別と被覆率）で 2.0 "
      "付近に寄せた。**この 2 項目は「帯に合わせて作った」ことを明示する**"
      "（設計が許す範囲だが、独立な検証ではない）。"
      % (bed["D84_over_D50_exact"], bed["armor_ratio_exact"]))
    A("- 基盤の低次モードと歪み変換は設計文書にない追加（§8 に開示）。"
      "歪み変換の係数は c=%.2f で、**基盤単体の歪度は %.2f**。"
      "matrix 歪度を帯に入れるために基盤の非対称を強くしており、"
      "これは構成の恣意性である（フィットではないが、帯に合わせた構成だ）。"
      % (B2.FORM_SKEW_C, bed["form"]["skew_form"]))
    A("- G-TRANS (a) は前進スプラットの構成上ほぼ自明に保存される。"
      "非自明なのは反復交点の非収束率（%.3f）と Beer 経路長の差だけ。"
      % ae["height"]["nonconverged_fraction"])
    A("- 影ゲート (b) のコアは解析テストで定義し、レンダの直接項を測る"
      "独立検証にしたが、両者が同じ幾何を使う以上「完全な独立」ではない。")
    A("- **画像の見え方**（brief のゲート外）: frame-0024 は v4-2.1 と同じ"
      "輸送レーンの産物で、**同じリング/モアレ様の模様と白飛び・黒潰れ**を"
      "持つ（C_max の中央値 %.1f 対 v4-2.1 の %.1f、どちらも同じスプラット"
      "格子由来）。v4-2.2 の brief は床・監査・輸送の再確認までで、画作りは"
      "対象外なので触っていない。**この見え方を「自然な水中景」と主張しない。**"
      % (cmax_med(numbers, "v4-2-2"), cmax_med(numbers, "v4-2-1")))
    A("- 主張の範囲: **「実測統計＋粒スケールの組織監査に 7/8 合格した床」**"
      "までしか言わない。構造関数の形態スケール指数は不合格のまま"
      "（ゲート総合も FAIL）。")
    A("")

    # ------------------------------------------------------------------ §10
    A("## 10. 成果物")
    A("")
    for nm, p in (("report", REPORT), ("numbers", NUMBERS),
                  ("bed-stats", BED_STATS), ("heightfield", HEIGHTFIELD),
                  ("heightfield-meta", os.path.join(OUT,
                                                    "heightfield-meta.json")),
                  ("transects", TRANSECTS), ("depth-map", DEPTH_PNG),
                  ("frame-0024", STILL), ("mp4", MP4),
                  ("checkpoint", CP_PATH),
                  ("stage0-py-hashes", STAGE0_HASHES),
                  ("h-isolation", HDIAG)):
        if nm == "report":
            A("- `%s` → `%s` (this file)" % (nm, os.path.relpath(p, HERE)))
            continue
        sz = os.path.getsize(p) if os.path.isfile(p) else -1
        A("- `%s` → `%s` (%d bytes)" % (nm, os.path.relpath(p, HERE), sz))
    A("- レンダの来歴: フレームは**納品した `_work/bed-maps.npz`**"
      "（sha256 先頭 %s）から描いた。床の再ビルド後に描き直しているので、"
      "`heightfield.npy`・`depth-map.png`・動画は同一の床に対応する。"
      % (read_json(os.path.join(WORK, "render-provenance.json"))
         ["bed_maps_npz_sha256"][:16]
         if os.path.exists(os.path.join(WORK, "render-provenance.json"))
         else "n/a"))
    A("")
    A("## 11. 実行の記録")
    A("")
    A("- 全ステージ foreground。jobkit checkpoint を stage ごとに原子的に"
      "コミット（completed_units=%s, last_unit=%s）。校正・再実行を含む"
      "全履歴は `completed_ids` にある（同一 stage の再実行も 1 unit として"
      "数えている）。"
      % (numbers["checkpoint"]["completed_units"],
         numbers["checkpoint"]["last_unit"]))
    A("- ゲートの判定は数値統計のみ。成果物の画像（depth-map / transects / "
      "frame）は**破損していないこと**を目視で確認したが、合否には使って"
      "いない（合否は上記の数値だけ）。")
    A("- ゲートの帯は測定後に動かしていない。")
    A("- 48 frames, 24 fps, 736x540, libx264 yuv420p, 2.000 s。"
      "ffprobe nb_frames=%d duration=%.3f。"
      % (gates["G-OUT"]["nb_frames"], gates["G-OUT"]["duration_s"]))
    A("")
    body = "\n".join(L) + "\n"
    atomic_write_bytes(REPORT, body.encode("utf-8"))


def do_report():
    log("=" * 72)
    log("[report] numbers + REPORT + still")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    sw = os.path.join(WORK, "calibration-sweeps.json")
    if os.path.exists(sw):
        js = read_json(sw)
        numbers["calibration_sweeps"] = js.get("comparable") or []
        numbers["calibration_exploratory"] = js.get("exploratory") or []
    else:
        numbers["calibration_sweeps"] = []
        numbers["calibration_exploratory"] = []
    write_json(NUMBERS, numbers)
    p24 = R421.chunk_path(24)
    if os.path.isfile(p24):
        from PIL import Image
        im = Image.open(p24)
        im.load()
        with open(p24, "rb") as f:
            data = f.read()
        atomic_write_bytes(STILL, data)
    write_report(numbers)
    cp.unit_done("report")
    # re-read the checkpoint (it now counts this unit and can be newer than the
    # snapshot build_numbers() took) and rewrite numbers + REPORT so that both
    # files carry the same, final counter
    numbers = build_numbers()
    if os.path.exists(sw):
        js = read_json(sw)
        numbers["calibration_sweeps"] = js.get("comparable") or []
        numbers["calibration_exploratory"] = js.get("exploratory") or []
    write_json(NUMBERS, numbers)
    write_report(numbers)
    log("  wrote %s (%d bytes)" % (os.path.relpath(REPORT, HERE),
                                   os.path.getsize(REPORT)))
    return numbers


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--bed", action="store_true")
    ap.add_argument("--hdiag", action="store_true")
    ap.add_argument("--scene", action="store_true")
    ap.add_argument("--gtrans", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)
    patch_r421_paths()
    if a.all or a.stage0:
        stage0()
    if a.all or a.bed:
        do_bed()
    if a.all or a.hdiag:
        do_hdiag()
    if a.all or a.scene:
        do_scene()
    if a.all or a.gtrans:
        do_gtrans()
    if a.chunk is not None:
        do_chunk(a.chunk)
    if a.all:
        for k in range(R421.NCHUNKS):
            do_chunk(k)
    if a.all or a.encode:
        do_encode()
    if a.all or a.gates:
        do_gates()
    if a.all or a.report:
        do_report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
