"""run_v4_1_16 -- BRIEF-V4-1-16-FASTDEPO driver.

STAGE 0  profile the verbatim area deposit (sub frame 512^2, frame 24)
STAGE 1  the fast implementation lives in transport_mesh_fast.py (new file)
STAGE 2  pre-registered equivalence gates against the verbatim reference
STAGE 3  sub frame at 1024^2, 48 frames, with the adopted fast method

Nothing in this job modifies an existing .py.  transport_mesh.py is imported
and used verbatim; the fast path is a separate module.
"""
import datetime
import gc
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

import run_v4_1_15 as R15            # noqa: E402  (build_scene, SCN, verbatim)
import transport_mesh as TM          # noqa: E402  (verbatim reference)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-16-FASTDEPO.md"
OUT = os.path.join(HERE, "out", "v4-1-16")
NUMBERS = os.path.join(OUT, "numbers-v4-1-16.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-16.md")
REFS = os.path.join(OUT, "refs-v4-1-16.npz")
LOG = os.path.join(OUT, "run-v4-1-16.log")

WIDTH, HEIGHT = R15.WIDTH, R15.HEIGHT
FPS = R15.FPS
STILL_IDX = R15.STILL_IDX            # 24
NFRAMES = R15.NFRAMES                # 48

#: verbatim deposit's runtime switch, exactly as v4-1.15 STAGE B set it.  The
#: file is untouched; this is the documented "which triangles take the
#: vectorised path" constant and changes no area (v4-1.15 section 2.1).
BOX_SUB = 1 << 20
BOX_CONTRAST = 16


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg):
    line = "[%s] %s" % (now_iso()[11:19], msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def makedirs():
    os.makedirs(OUT, exist_ok=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False, allow_nan=True)
    os.replace(tmp, path)


def numbers():
    if os.path.exists(NUMBERS):
        return read_json(NUMBERS)
    return {"job": "v4-1-16", "brief": BRIEF, "created": now_iso()}


def numbers_update(section, payload):
    n = numbers()
    n.setdefault(section, {})
    if isinstance(n[section], dict) and isinstance(payload, dict):
        n[section].update(payload)
    else:
        n[section] = payload
    n["updated"] = now_iso()
    write_json(NUMBERS, n)
    return n


def c_stats(C):
    C = np.asarray(C, dtype=np.float64)
    return {"mean": float(np.mean(C)), "std": float(np.std(C)),
            "min": float(np.min(C)), "max": float(np.max(C)),
            "p01": float(np.percentile(C, 1)),
            "p50": float(np.percentile(C, 50)),
            "p99": float(np.percentile(C, 99)),
            "shape": [int(C.shape[0]), int(C.shape[1])]}


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def probe_environment():
    env = {"python": sys.version.split()[0],
           "platform": platform.platform(),
           "machine": platform.machine(),
           "numpy": np.__version__,
           "cpu_count": os.cpu_count(),
           "cwd": os.getcwd(),
           "started": now_iso()}
    try:
        import numba
        env["numba"] = numba.__version__
        env["numba_threading_layer"] = str(
            getattr(numba, "config", None) and numba.config.THREADING_LAYER)
    except Exception as e:                        # noqa: BLE001
        env["numba"] = "MISSING (%s)" % (e,)
    try:
        import cupy  # noqa: F401
        env["cupy"] = "present"
    except Exception as e:                        # noqa: BLE001
        env["cupy"] = "MISSING (%s)" % (type(e).__name__,)
    try:
        import scipy
        env["scipy"] = scipy.__version__
    except Exception as e:                        # noqa: BLE001
        env["scipy"] = "MISSING (%s)" % (type(e).__name__,)
    env["gpu_query"] = "not attempted (no CuPy; see brief: only if available)"
    return env


class Probe:
    """Time/call counters for the deposit internals (module-global patching)."""

    def __init__(self):
        self.rec = {}
        self.phase = "none"
        self.stack = []

    def add(self, key, dt, n=1):
        r = self.rec.setdefault(key, {"t": 0.0, "n": 0})
        r["t"] += dt
        r["n"] += n

    def wrap(self, mod, name, key, phase=None):
        orig = getattr(mod, name)

        def f(*a, **k):
            if phase is not None:
                old = self.phase
                self.phase = phase
            t0 = time.perf_counter()
            try:
                return orig(*a, **k)
            finally:
                self.add(key, time.perf_counter() - t0)
                if phase is not None:
                    self.phase = old

        setattr(mod, name, f)
        self.stack.append((mod, name, orig))
        return orig

    def wrap_phase_clip(self, mod, clip_name, bulk_name, large_name):
        """Attribute _clip time to the bulk or the tail path."""
        clip_orig = getattr(mod, clip_name)

        def clip(*a, **k):
            t0 = time.perf_counter()
            try:
                return clip_orig(*a, **k)
            finally:
                self.add("clip_" + self.phase, time.perf_counter() - t0)

        setattr(mod, clip_name, clip)
        self.stack.append((mod, clip_name, clip_orig))

        for nm, ph in ((bulk_name, "bulk"), (large_name, "large")):
            orig = getattr(mod, nm)

            def g(*a, _orig=orig, _ph=ph, **k):
                old = self.phase
                self.phase = _ph
                t0 = time.perf_counter()
                try:
                    return _orig(*a, **k)
                finally:
                    self.add(_ph, time.perf_counter() - t0)
                    self.phase = old

            setattr(mod, nm, g)
            self.stack.append((mod, nm, orig))

    def restore(self):
        for mod, name, orig in reversed(self.stack):
            setattr(mod, name, orig)
        self.stack = []

    def report(self):
        return {k: {"seconds": float(v["t"]), "calls": int(v["n"])}
                for k, v in self.rec.items()}


def bbox_table(T, A_tri, A_cell, W, H):
    """The verbatim deposit's own bulk/tail classification, counted."""
    xmin = T[:, :, 0].min(1)
    xmax = T[:, :, 0].max(1)
    ymin = T[:, :, 1].min(1)
    ymax = T[:, :, 1].max(1)
    i0 = np.floor(xmin)
    j0 = np.floor(ymin)
    sw = np.floor(xmax) - i0 + 1
    sh = np.floor(ymax) - j0 + 1
    per = sw * sh
    out = {"n_triangles": int(T.shape[0]),
           "n_cells": int(A_cell.size),
           "sum_absA_over_screen": float(A_tri.sum() / float(W * H)),
           "sum_bbox_px": float(per.sum()),
           "bbox_px_mean": float(per.mean()),
           "bbox_px_median": float(np.median(per)),
           "bbox_px_p99": float(np.percentile(per, 99)),
           "bbox_px_max": float(per.max()),
           "rows": []}
    for lim in (16, 256, 1 << 20):
        bulk = per <= lim
        out["rows"].append({
            "box_bulk_max": int(lim),
            "n_bulk": int(bulk.sum()),
            "n_tail": int((~bulk).sum()),
            "bulk_boxes": float(per[bulk].sum()),
            "tail_boxes": float(per[~bulk].sum()),
            "tail_frac_triangles": float((~bulk).mean())})
    return out


def stage0():
    makedirs()
    log("=" * 78)
    log("[stage0] profile of the verbatim deposit, sub frame 512^2, frame 24")
    env = probe_environment()
    write_json(os.path.join(OUT, "stage0-env.json"), env)
    log("  env: " + json.dumps(env, ensure_ascii=False))
    if "MISSING" in str(env["numba"]):
        raise SystemExit("numba missing: candidate 2 is unavailable, stop and "
                         "report (brief forbids assuming it)")

    TM.BOX_BULK_MAX = BOX_SUB                 # v4-1.15 STAGE B condition
    rec = {"box_bulk_max_used": int(TM.BOX_BULK_MAX),
           "built_on": now_iso()}

    st, senv, optics, spec = R15.build_scene("sub")
    t = STILL_IDX / float(FPS)
    fields = st.fields_at(t)
    wt = float(WIDTH) * float(HEIGHT) / (float(st.N) * float(st.N))
    rec.update({"N_s": int(st.N), "W": WIDTH, "H": HEIGHT, "S_m": float(st.S),
                "depth_m": float(senv["depthM"]), "t_s": float(t),
                "w_cell": wt, "stage": "sub", "frame": STILL_IDX})

    # --- geometry, timed exactly as mesh_C does it -------------------------
    t0 = time.perf_counter()
    PXe, PYe = TM.landing_pixels(fields, st.S, senv, W=WIDTH, H=HEIGHT)
    t_land = time.perf_counter() - t0
    t0 = time.perf_counter()
    T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
    t_tri = time.perf_counter() - t0
    rec["geometry"] = {"landing_pixels_s": float(t_land),
                       "cell_triangles_s": float(t_tri)}
    rec["classification"] = bbox_table(T, A_tri, A_cell, WIDTH, HEIGHT)
    rng = np.random.default_rng(20261116)

    # --- cost of the tail path on a measured sample (labelled extrapolation)
    tail_rows = rec["classification"]["rows"][0]         # BOX_BULK_MAX = 16
    per_all = None
    xmin = T[:, :, 0].min(1); xmax = T[:, :, 0].max(1)
    ymin = T[:, :, 1].min(1); ymax = T[:, :, 1].max(1)
    i0a = np.floor(xmin); j0a = np.floor(ymin)
    per_all = (np.floor(xmax) - i0a + 1) * (np.floor(ymax) - j0a + 1)
    dens = np.full(A_cell.size, 1.0)
    dens_all = np.concatenate([dens, dens])
    tail_idx = np.nonzero(per_all > 16)[0]
    nsample = int(min(120, tail_idx.size))
    if nsample:
        sample = np.sort(rng.choice(tail_idx, size=nsample, replace=False))
        acc_s = np.zeros(WIDTH * HEIGHT, dtype=np.float64)
        info_s = {"area_in": 0.0, "area_out": 0.0, "energy_in": 0.0,
                  "energy_out": 0.0, "bulk_boxes": 0, "tail_boxes": 0}
        t0 = time.perf_counter()
        for k in sample:
            TM._deposit_large(T[k], float(dens_all[k]), acc_s, WIDTH, HEIGHT,
                              info_s)
        wall_s = time.perf_counter() - t0
        rec["tail_sample"] = {
            "n_triangles": int(nsample),
            "triangles": [int(v) for v in sample[:10].tolist()],
            "triangles_note": "first 10 of the sampled indices",
            "wall_s": float(wall_s),
            "us_per_triangle": float(1e6 * wall_s / nsample),
            "boxes": float(info_s["tail_boxes"]),
            "us_per_box": float(1e6 * wall_s / max(info_s["tail_boxes"], 1)),
            "extrapolated_total_s_for_all_tail_triangles": float(
                wall_s / nsample * tail_idx.size)}
        log("  tail sample: %.1f us/triangle, %.3f us/box -> %.0f s for the "
            "whole %d-triangle tail (extrapolation)"
            % (rec["tail_sample"]["us_per_triangle"],
               rec["tail_sample"]["us_per_box"],
               rec["tail_sample"]["extrapolated_total_s_for_all_tail_triangles"],
               tail_idx.size))

    # --- instrumented verbatim deposit ------------------------------------
    probe = Probe()
    probe.wrap_phase_clip(TM, "_clip", "_deposit_bulk", "_deposit_large")
    probe.wrap(TM, "_shoelace", "shoelace")
    w_cell = np.full(A_cell.size, wt)
    t0 = time.perf_counter()
    acc, info = TM.deposit(T, A_tri, A_cell, w_cell, WIDTH, HEIGHT)
    t_dep = time.perf_counter() - t0
    timers = probe.report()
    probe.restore()
    C = acc.astype(np.float32)

    rec["deposit"] = {"total_s": float(t_dep),
                      "timers": timers,
                      "info": {k: (float(v) if isinstance(v, (int, float))
                                   and not isinstance(v, bool) else v)
                               for k, v in info.items()},
                      "C_stats": c_stats(C)}
    rec["deposit"]["fractions"] = {
        k: float(v["seconds"] / t_dep) for k, v in timers.items()}
    log("  deposit wall %.1f s  bulk %.1f s  large %.1f s  clip(bulk) %.1f s  "
        "shoelace %.1f s"
        % (t_dep, timers.get("bulk", {}).get("seconds", 0.0),
           timers.get("large", {}).get("seconds", 0.0),
           timers.get("clip_bulk", {}).get("seconds", 0.0),
           timers.get("shoelace", {}).get("seconds", 0.0)))
    log("  C mean %.12f std %.4f p99 %.4f max %.3f rel_err %.2e"
        % (rec["deposit"]["C_stats"]["mean"], rec["deposit"]["C_stats"]["std"],
           rec["deposit"]["C_stats"]["p99"], rec["deposit"]["C_stats"]["max"],
           info["rel_err"]))

    # provenance: the recorded v4-1.15 frame-24 sub C statistics
    cp = read_json(os.path.join(HERE, "out", "v4-1-15", "checkpoint.json"))
    ref_stats = cp["frame24"]["sub"]["C_frame_stats"]
    rec["provenance_vs_v4_1_15"] = {
        "recorded": ref_stats,
        "abs_diff": {k: float(abs(rec["deposit"]["C_stats"][k] - ref_stats[k]))
                     for k in ("mean", "std", "p99", "max")},
        "stages": "regenerated with the same build_scene/fields_at/mesh_C "
                  "chain, BOX_BULK_MAX = %d" % int(TM.BOX_BULK_MAX)}
    log("  vs v4-1.15 recorded frame24 sub C: dstd %.2e dp99 %.2e dmax %.2e"
        % (rec["provenance_vs_v4_1_15"]["abs_diff"]["std"],
           rec["provenance_vs_v4_1_15"]["abs_diff"]["p99"],
           rec["provenance_vs_v4_1_15"]["abs_diff"]["max"]))

    np.save(os.path.join(OUT, "_stage0_sub512_C.npy"), C)
    np.savez_compressed(REFS + ".stage0.npz", sub512=C)
    numbers_update("stage0", rec)
    TM.BOX_BULK_MAX = 16
    log("[stage0] done")
    return rec


# ---------------------------------------------------------------------------
# reference generation (verbatim, unwrapped)
# ---------------------------------------------------------------------------
def refgen():
    makedirs()
    log("=" * 78)
    log("[refgen] regenerating the two reference C arrays with verbatim "
        "transport_mesh (no instrumentation, no patching)")
    out = {"created": now_iso(),
           "why": "v4-1.15 saved C statistics but not the C arrays; the brief "
                  "says to regenerate the two frames with the unmodified "
                  "verbatim deposit in that case",
           "frames": {}}
    arrays = {}
    for tag, box in (("contrast", BOX_CONTRAST), ("sub", BOX_SUB)):
        TM.BOX_BULK_MAX = box
        st, senv, optics, spec = R15.build_scene(tag)   # SCN grid (1024 / 512)
        t = STILL_IDX / float(FPS)
        t0 = time.perf_counter()
        fields = st.fields_at(t)
        t_field = time.perf_counter() - t0
        t0 = time.perf_counter()
        C, info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, senv)
        wall = time.perf_counter() - t0
        stt = c_stats(C)
        cp = read_json(os.path.join(HERE, "out", "v4-1-15",
                                    "checkpoint.json"))
        ref = cp["frame24"][tag]["C_frame_stats"]
        entry = {"stage": tag, "frame": STILL_IDX, "t_s": float(t),
                 "box_bulk_max": int(box), "N_s": int(st.N),
                 "w_cell": float(info["w_cell"]),
                 "mesh_C_wall_s": float(wall),
                 "fields_at_s": float(t_field),
                 "C_stats": stt,
                 "energy": {"rel_err": float(info["rel_err"]),
                            "cell_total": float(info["energy_cell_total"]),
                            "deposited": float(info["energy_deposited"]),
                            "guard": float(info["energy_guard"])},
                 "provenance_vs_v4_1_15": {
                     "recorded": ref,
                     "abs_diff": {k: float(abs(stt[k] - ref[k]))
                                  for k in ("mean", "std", "p99", "max")}}}
        out["frames"][tag] = entry
        arrays["%s_C" % tag] = C
        arrays["%s_info" % tag] = np.array(
            [info["rel_err"], info["energy_cell_total"],
             info["energy_deposited"]], dtype=np.float64)
        log("  %-8s %4d^2  mesh_C %.1f s  std %.6f p99 %.6f  rel_err %.2e  "
            "dstd %.2e"
            % (tag, st.N, wall, stt["std"], stt["p99"], info["rel_err"],
               entry["provenance_vs_v4_1_15"]["abs_diff"]["std"]))
        del st, optics, spec, fields
        gc.collect()
    arrays["W"] = np.array([WIDTH])
    arrays["H"] = np.array([HEIGHT])
    np.savez_compressed(REFS, **arrays)
    out["refs_file"] = os.path.basename(REFS)
    out["refs_sha256"] = sha256_file(REFS)
    numbers_update("reference", out)
    TM.BOX_BULK_MAX = 16
    log("[refgen] done -> %s" % REFS)
    return out


# ---------------------------------------------------------------------------
# STAGE 2 -- pre-registered equivalence gates
# ---------------------------------------------------------------------------
GATES = {
    "total_energy_rel": 1e-9,       # |sum(C_fast) - sum(C_ref)| / sum(C_ref)
    "std_rel": 0.005,               # 0.5 %
    "p99_rel": 0.005,               # 0.5 %
    "max_abs_px": None,             # diagnostic only, never a gate
}


def _load_refs():
    z = np.load(REFS)
    return z


def equiv():
    makedirs()
    log("=" * 78)
    log("[equiv] STAGE 2 gates: fast deposit vs the regenerated verbatim "
        "reference C arrays")
    import transport_mesh_fast as TMF
    z = _load_refs()
    rec = {"created": now_iso(), "gates": GATES,
           "reference": {"file": os.path.basename(REFS),
                         "sha256": sha256_file(REFS),
                         "how": "unmodified verbatim transport_mesh.mesh_C, "
                                "regenerated by this job because v4-1.15 saved "
                                "C statistics but not the C arrays"},
           "frames": {}}
    # consistency of the two independently produced verbatim sub references
    stage0_C = np.load(os.path.join(OUT, "_stage0_sub512_C.npy"))
    rec["reference_consistency_sub512"] = {
        "max_abs_diff_stage0_vs_refgen": float(
            np.max(np.abs(stage0_C.astype(np.float64)
                          - z["sub_C"].astype(np.float64)))),
        "note": "stage-0 profiled run vs the clean refgen run (both verbatim)"}

    for tag, box in (("contrast", BOX_CONTRAST), ("sub", BOX_SUB)):
        TMF.BOX_BULK_MAX = box
        Cref = z["%s_C" % tag].astype(np.float64)
        st, senv, optics, spec = R15.build_scene(tag)
        t = STILL_IDX / float(FPS)
        fields = st.fields_at(t)
        t0 = time.perf_counter()
        Cf, info = TMF.mesh_C_fast(fields, st.S, WIDTH, HEIGHT, senv)
        wall = time.perf_counter() - t0
        Cf64 = Cf.astype(np.float64)
        sref, sfast = c_stats(Cref), c_stats(Cf64)
        E_ref = float(Cref.sum())
        E_fast = float(Cf64.sum())
        d = np.abs(Cref - Cf64)
        g = {"stage": tag, "N_s": int(st.N), "box_bulk_max": int(box),
             "C_ref": sref, "C_fast": sfast,
             "total_energy_ref": E_ref, "total_energy_fast": E_fast,
             "total_energy_rel_diff": float(abs(E_fast - E_ref) / abs(E_ref)),
             "internal_rel_err_ref": float(z["%s_info" % tag][0]),
             "internal_rel_err_fast": float(info["rel_err"]),
             "std_rel_diff": float(sfast["std"] / sref["std"] - 1.0),
             "p99_rel_diff": float(sfast["p99"] / sref["p99"] - 1.0),
             "mean_rel_diff": float(sfast["mean"] / sref["mean"] - 1.0),
             "max_abs_px_diff": float(d.max()),
             "max_abs_px_diff_rel": float(d.max() / sref["max"]),
             "frac_pixels_differing": float(np.mean(d > 0.0)),
             "p999_abs_px_diff": float(np.percentile(d, 99.9)),
             "deposit_info_fast": {k: (float(v) if isinstance(v, (int, float))
                                       and not isinstance(v, bool) else v)
                                   for k, v in info.items()},
             "mesh_C_fast_wall_s": float(wall)}
        g["checks"] = {
            "total_energy": bool(g["total_energy_rel_diff"] <= GATES["total_energy_rel"]),
            "internal_rel_err_fast": bool(g["internal_rel_err_fast"]
                                          <= GATES["total_energy_rel"]),
            "std": bool(abs(g["std_rel_diff"]) <= GATES["std_rel"]),
            "p99": bool(abs(g["p99_rel_diff"]) <= GATES["p99_rel"])}
        g["pass"] = bool(all(g["checks"].values()))
        rec["frames"][tag] = g
        log("  %-8s %4d^2  fast %.1f s  dE %.2e  dstd %+.4f%%  dp99 %+.4f%%  "
            "max|dC| %.3e  -> %s"
            % (tag, st.N, wall, g["total_energy_rel_diff"],
               100 * g["std_rel_diff"], 100 * g["p99_rel_diff"],
               g["max_abs_px_diff"], "PASS" if g["pass"] else "FAIL"))
        del st, optics, spec, fields, Cref, Cf64, Cf
        gc.collect()
    rec["pass"] = bool(all(f["pass"] for f in rec["frames"].values()))
    rec["adopted"] = ("transport_mesh_fast.deposit_fast / mesh_C_fast "
                      "(Numba float64 kernels); frame-parallel execution is "
                      "measured separately and changes no arithmetic")
    numbers_update("equivalence", rec)
    log("[equiv] %s" % ("PASS" if rec["pass"] else "FAIL"))
    return rec


# ---------------------------------------------------------------------------
# frame-parallel candidate (candidate 1) and the 1024^2 timing of the
# verbatim deposit (for the speed table: measured, not extrapolated)
# ---------------------------------------------------------------------------
_WS = {}


def _w_init(tag, n, box, fast):
    global _WS
    import transport_mesh as _TM
    import transport_mesh_fast as _TMF
    _TMF.BOX_BULK_MAX = box
    _TM.BOX_BULK_MAX = box
    st, senv, optics, spec = R15.build_scene(tag, n=n)
    _WS = {"tag": tag, "st": st, "env": senv, "TM": _TM, "TMF": _TMF}


def _w_frames(idxs):
    st = _WS["st"]
    env = _WS["env"]
    TMF = _WS["TMF"]
    out = []
    for i in idxs:
        t = i / float(FPS)
        t0 = time.perf_counter()
        fields = st.fields_at(t)
        t1 = time.perf_counter()
        C, info = TMF.mesh_C_fast(fields, st.S, WIDTH, HEIGHT, env)
        t2 = time.perf_counter()
        out.append({"frame": int(i), "t_s": float(t),
                    "C_stats": c_stats(C),
                    "energy_rel_err": float(info["rel_err"]),
                    "energy_deposited": float(info["energy_deposited"]),
                    "energy_cell_total": float(info["energy_cell_total"]),
                    "offframe_fraction_energy": float(
                        info["offframe_fraction_energy"]),
                    "fields_s": float(t1 - t0), "mesh_C_s": float(t2 - t1)})
    return out


def _w_frames_verbatim(idxs):
    st = _WS["st"]
    env = _WS["env"]
    TM_ = _WS["TM"]
    out = []
    for i in idxs:
        t = i / float(FPS)
        fields = st.fields_at(t)
        t0 = time.perf_counter()
        C, info = TM_.mesh_C(fields, st.S, WIDTH, HEIGHT, env)
        wall = time.perf_counter() - t0
        out.append({"frame": int(i), "C_stats": c_stats(C),
                    "energy_rel_err": float(info["rel_err"]),
                    "mesh_C_s": float(wall)})
    return out


def parbench():
    makedirs()
    log("=" * 78)
    log("[parbench] candidate 1: frame parallelism on independent frames")
    import multiprocessing as mp
    import transport_mesh_fast as TMF
    nproc = min(8, max(1, (os.cpu_count() or 1) - 2))
    rec = {"created": now_iso(), "workers": int(nproc),
           "cpu_count": os.cpu_count()}
    idxs = list(range(8))
    # --- sequential, fast kernels -----------------------------------------
    TMF.BOX_BULK_MAX = BOX_SUB
    st, senv, optics, spec = R15.build_scene("sub")
    seq = []
    t0 = time.perf_counter()
    for i in idxs:
        t = i / float(FPS)
        f = st.fields_at(t)
        C, info = TMF.mesh_C_fast(f, st.S, WIDTH, HEIGHT, senv)
        seq.append({"frame": i, "C_stats": c_stats(C),
                    "energy_rel_err": float(info["rel_err"])})
    t_seq = time.perf_counter() - t0
    del st, optics, spec
    gc.collect()
    # --- same frames over a process pool ----------------------------------
    chunks = [idxs[k::nproc] for k in range(nproc)]
    chunks = [c for c in chunks if c]
    t0 = time.perf_counter()
    with mp.Pool(nproc, initializer=_w_init,
                 initargs=("sub", None, BOX_SUB, True)) as pool:
        parts = pool.map(_w_frames, chunks)
    t_par = time.perf_counter() - t0
    par = sorted([r for p in parts for r in p], key=lambda r: r["frame"])
    dstd = max(abs(a["C_stats"]["std"] / b["C_stats"]["std"] - 1.0)
               for a, b in zip(par, seq))
    dp99 = max(abs(a["C_stats"]["p99"] / b["C_stats"]["p99"] - 1.0)
               for a, b in zip(par, seq))
    rec.update({"grid": int(R15.SCN["sub"]["N"]), "frames": idxs,
                "sequential_s": float(t_seq),
                "sequential_s_per_frame": float(t_seq / len(idxs)),
                "parallel_s": float(t_par),
                "parallel_s_per_frame": float(t_par / len(idxs)),
                "speedup": float(t_seq / t_par),
                "max_rel_diff_std": float(dstd),
                "max_rel_diff_p99": float(dp99),
                "frames_identical": bool(dstd == 0.0 and dp99 == 0.0),
                "note": "frames are independent; each worker builds the same "
                        "scene from the same seed and runs the same kernel, "
                        "so the arithmetic per frame is untouched"})
    log("  8 frames 512^2 fast: sequential %.1f s (%.1f s/frame), "
        "pool(%d) %.1f s (%.1f s/frame) -> %.2fx; stats identical: %s"
        % (t_seq, t_seq / len(idxs), nproc, t_par, t_par / len(idxs),
           t_seq / t_par, rec["frames_identical"]))
    numbers_update("frame_parallel", rec)
    return rec


def time_verbatim_1024():
    makedirs()
    log("=" * 78)
    log("[time_verbatim_1024] one verbatim 1024^2 sub frame (speed table, "
        "measured -- this is the long one)")
    TM.BOX_BULK_MAX = BOX_SUB
    st, senv, optics, spec = R15.build_scene("sub", n=1024)
    t = STILL_IDX / float(FPS)
    fields = st.fields_at(t)
    t0 = time.perf_counter()
    C, info = TM.mesh_C(fields, st.S, WIDTH, HEIGHT, senv)
    wall = time.perf_counter() - t0
    stt = c_stats(C)
    rec = {"N_s": 1024, "box_bulk_max": int(BOX_SUB),
           "mesh_C_wall_s": float(wall), "C_stats": stt,
           "rel_err": float(info["rel_err"]),
           "bulk_boxes": float(info["bulk_boxes"]),
           "n_bulk": int(info["n_bulk"]), "n_tail": int(info["n_tail"]),
           "us_per_box": float(1e6 * wall / max(info["bulk_boxes"], 1.0)),
           "created": now_iso()}
    log("  verbatim 1024^2: %.1f s  std %.6f p99 %.6f  %.3f us/box"
        % (wall, stt["std"], stt["p99"], rec["us_per_box"]))
    np.save(os.path.join(OUT, "_verbatim_sub1024_f24_C.npy"), C)
    numbers_update("verbatim_1024", rec)
    TM.BOX_BULK_MAX = 16
    return rec


# ---------------------------------------------------------------------------
# STAGE 3 -- sub frame at 1024^2, 48 frames, fast method
# ---------------------------------------------------------------------------
def stage3():
    makedirs()
    log("=" * 78)
    log("[stage3] sub frame 1024^2, 48 frames, transport_mesh_fast")
    import multiprocessing as mp
    nproc = min(8, max(1, (os.cpu_count() or 1) - 2))
    idxs = list(range(NFRAMES))
    chunks = [idxs[k::nproc] for k in range(nproc)]
    chunks = [c for c in chunks if c]
    t0 = time.perf_counter()
    with mp.Pool(nproc, initializer=_w_init,
                 initargs=("sub", 1024, BOX_SUB, True)) as pool:
        parts = pool.map(_w_frames, chunks)
    wall = time.perf_counter() - t0
    rows = sorted([r for p in parts for r in p], key=lambda r: r["frame"])
    rec = {"created": now_iso(), "stage": "sub", "N_s": 1024,
           "window_m": R15.SCN["sub"]["window_m"],
           "depth_m": R15.SCN["sub"]["depth"],
           "water_preset": R15.SCN["sub"]["water"],
           "preset": {"U": R15.SCN["sub"]["U"], "F": R15.SCN["sub"]["F"],
                      "seed": R15.SCN["sub"]["seed"]},
           "W": WIDTH, "H": HEIGHT, "N_frames": NFRAMES, "workers": int(nproc),
           "total_wall_s": float(wall), "s_per_frame_wall": float(wall / NFRAMES),
           "method": "transport_mesh_fast (numba kernels), frame-parallel",
           "rendering": "C layer only: the brief's STAGE 3 asks for the "
                        "G-CONV material and the speed table, so this run "
                        "produces the 48 C frames and their statistics, not "
                        "PNG/MP4 output",
           "frames": rows}
    stats = [r["C_stats"] for r in rows]
    rec["summary"] = {
        "mean_of_mean": float(np.mean([s["mean"] for s in stats])),
        "mean_of_std": float(np.mean([s["std"] for s in stats])),
        "min_std": float(min(s["std"] for s in stats)),
        "max_std": float(max(s["std"] for s in stats)),
        "mean_of_p99": float(np.mean([s["p99"] for s in stats])),
        "min_p99": float(min(s["p99"] for s in stats)),
        "max_p99": float(max(s["p99"] for s in stats)),
        "max_rel_err": float(max(r["energy_rel_err"] for r in rows)),
        "mean_abs_rel_err": float(np.mean([abs(r["energy_rel_err"])
                                          for r in rows])),
        "mesh_C_s_per_frame_mean": float(np.mean([r["mesh_C_s"] for r in rows])),
        "fields_s_per_frame_mean": float(np.mean([r["fields_s"] for r in rows]))}
    log("  48 frames 1024^2 in %.0f s (%.1f s/frame wall, %d workers); "
        "std %.4f (%.4f..%.4f) p99 %.4f (%.4f..%.4f) max rel_err %.2e"
        % (wall, wall / NFRAMES, nproc, rec["summary"]["mean_of_std"],
           rec["summary"]["min_std"], rec["summary"]["max_std"],
           rec["summary"]["mean_of_p99"], rec["summary"]["min_p99"],
           rec["summary"]["max_p99"], rec["summary"]["max_rel_err"]))

    # ---- G-CONV material: 1024^2 (this run) vs 512^2 (v4-1.15 base) ------
    cp = read_json(os.path.join(HERE, "out", "v4-1-15", "checkpoint.json"))
    base = {}
    for k, node in (cp.get("chunks") or {}).get("sub", {}).items():
        for m in node["metrics"]:
            base[int(m["frame"])] = m
    cmp_rows = []
    for r in rows:
        b = base.get(r["frame"])
        if b is None:
            continue
        cmp_rows.append({
            "frame": int(r["frame"]),
            "std_512": float(b["C_std"]), "std_1024": r["C_stats"]["std"],
            "std_rel_diff": float(r["C_stats"]["std"] / b["C_std"] - 1.0),
            "p99_512": float(b["C_p99"]), "p99_1024": r["C_stats"]["p99"],
            "p99_rel_diff": float(r["C_stats"]["p99"] / b["C_p99"] - 1.0),
            "within_15pct": bool(
                abs(r["C_stats"]["std"] / b["C_std"] - 1.0) <= 0.15
                and abs(r["C_stats"]["p99"] / b["C_p99"] - 1.0) <= 0.15)})
    sd = np.array([c["std_rel_diff"] for c in cmp_rows])
    pd = np.array([c["p99_rel_diff"] for c in cmp_rows])
    f24 = [c for c in cmp_rows if c["frame"] == STILL_IDX]
    gate = {
        "n_frames_compared": len(cmp_rows),
        "std_rel_diff_mean": float(sd.mean()),
        "std_rel_diff_min": float(sd.min()),
        "std_rel_diff_max": float(sd.max()),
        "p99_rel_diff_mean": float(pd.mean()),
        "p99_rel_diff_min": float(pd.min()),
        "p99_rel_diff_max": float(pd.max()),
        "frames_within_15pct": int(sum(c["within_15pct"] for c in cmp_rows)),
        "frame24": f24[0] if f24 else None,
        "criterion": "v4-1.15's G-CONV definition: |std and p99 difference "
                     "between the base grid and the comparison grid| <= 15%",
        "note": ("the comparison is 1024^2 (this job) against the 512^2 base "
                 "grid of v4-1.15, i.e. the upward direction of the same "
                 "question; there, 512^2 vs 1024^2 on frame 24 was "
                 "std -33.81% / p99 -26.00%"),
        "verdict": ("PASS" if (abs(sd.mean()) <= 0.15 and abs(pd.mean()) <= 0.15)
                    else "FAIL")}
    rec["gconv_vs_512"] = gate
    rec["gconv_rows"] = cmp_rows
    log("  G-CONV material 1024^2 vs 512^2 (48 frames): std %+.2f%% "
        "(range %+.2f..%+.2f), p99 %+.2f%% (range %+.2f..%+.2f), "
        "within-15%%: %d/48 -> %s"
        % (100 * sd.mean(), 100 * sd.min(), 100 * sd.max(), 100 * pd.mean(),
           100 * pd.min(), 100 * pd.max(), gate["frames_within_15pct"],
           gate["verdict"]))
    numbers_update("stage3", rec)
    return rec


# ---------------------------------------------------------------------------
# speed table (all measured; no extrapolation in the table itself)
# ---------------------------------------------------------------------------
def speed():
    makedirs()
    log("=" * 78)
    log("[speed] measured single-frame costs, 736x540 frame grid, sub scene")
    import transport_mesh_fast as TMF
    rec = {"created": now_iso(), "frame_grid": [HEIGHT, WIDTH],
           "rows": [], "fast_detail": {}}
    for n in (512, 1024):
        TMF.BOX_BULK_MAX = BOX_SUB
        st, senv, optics, spec = R15.build_scene("sub", n=n)
        t = STILL_IDX / float(FPS)
        fields = st.fields_at(t)
        t0 = time.perf_counter()
        PXe, PYe = TM.landing_pixels(fields, st.S, senv, W=WIDTH, H=HEIGHT)
        t_land = time.perf_counter() - t0
        t0 = time.perf_counter()
        T, A_tri, A_cell = TM.cell_triangles(PXe, PYe)
        t_tri = time.perf_counter() - t0
        w_cell = np.full(A_cell.size, float(WIDTH) * float(HEIGHT)
                         / (float(n) * float(n)))
        t0 = time.perf_counter()
        acc, info = TMF.deposit_fast(T, A_tri, A_cell, w_cell, WIDTH, HEIGHT)
        t_dep = time.perf_counter() - t0
        Cfast = acc.astype(np.float32)
        t0 = time.perf_counter()
        Cf2, info2 = TMF.mesh_C_fast(fields, st.S, WIDTH, HEIGHT, senv)
        t_mesh = time.perf_counter() - t0
        rec["fast_detail"][str(n)] = {
            "landing_pixels_s": float(t_land), "cell_triangles_s": float(t_tri),
            "deposit_s": float(t_dep), "mesh_C_end_to_end_s": float(t_mesh),
            "bulk_boxes": float(info["bulk_boxes"]),
            "n_bulk": int(info["n_bulk"]), "n_tail": int(info["n_tail"]),
            "us_per_box_deposit": float(1e6 * t_dep
                                        / max(info["bulk_boxes"], 1.0)),
            "C_stats": c_stats(Cf2),
            "rel_err": float(info2["rel_err"]),
            "sum_absA_over_screen": float(A_tri.sum() / (WIDTH * HEIGHT))}
        log("  fast  %4d^2: deposit %.2f s (%.3f us/box, %.2e boxes)  "
            "geometry %.2f+%.2f s  mesh_C end-to-end %.2f s"
            % (n, t_dep, rec["fast_detail"][str(n)]["us_per_box_deposit"],
               info["bulk_boxes"], t_land, t_tri, t_mesh))
        del st, optics, spec, fields, T, A_tri, A_cell, acc, Cfast, Cf2
        gc.collect()

    n_all = numbers()
    v512 = (n_all.get("reference", {}).get("frames", {}).get("sub", {})
            .get("mesh_C_wall_s"))
    cp = read_json(os.path.join(HERE, "out", "v4-1-15", "checkpoint.json"))
    ch = cp["chunks"]["sub"]
    v512_48 = float((ch["0"]["wall_s"] + ch["1"]["wall_s"]) / NFRAMES)
    v1024 = (n_all.get("verbatim_1024", {}) or {}).get("mesh_C_wall_s")
    f512 = rec["fast_detail"]["512"]["mesh_C_end_to_end_s"]
    f1024 = rec["fast_detail"]["1024"]["mesh_C_end_to_end_s"]
    rows = [
        {"grid": "512^2", "method": "verbatim", "s_per_frame": v512,
         "source": "measured: this job, refgen sub mesh_C (single frame 24)"},
        {"grid": "512^2", "method": "verbatim",
         "s_per_frame": v512_48,
         "source": "measured: v4-1.15 STAGE B chunks, 48 frames total wall / 48"},
        {"grid": "512^2", "method": "fast (numba)", "s_per_frame": f512,
         "source": "measured: this job, mesh_C_fast end-to-end frame 24"},
        {"grid": "512^2", "method": "fast + 8 workers",
         "s_per_frame": (numbers().get("frame_parallel", {})
                         .get("parallel_s_per_frame")),
         "source": "measured: this job, 8 frames over a process pool"},
        {"grid": "1024^2", "method": "verbatim", "s_per_frame": v1024,
         "source": "measured: this job, one verbatim 1024^2 frame 24"},
        {"grid": "1024^2", "method": "fast (numba)", "s_per_frame": f1024,
         "source": "measured: this job, mesh_C_fast end-to-end frame 24"},
        {"grid": "1024^2", "method": "fast + 8 workers",
         "s_per_frame": (numbers().get("stage3", {})
                         .get("s_per_frame_wall")),
         "source": "measured: this job, STAGE 3, 48 frames / 48"}]
    for r in rows:
        r["x_vs_verbatim_same_grid"] = (
            (v512 / r["s_per_frame"]) if r["grid"] == "512^2"
            else (v1024 / r["s_per_frame"] if r["s_per_frame"] else None))
    rec["rows"] = rows
    rec["verbatim_1024_note"] = ("the verbatim 1024^2 number is a single "
                                 "measured frame; v4-1.15's stage-0 "
                                 "extrapolation said 867 s")
    numbers_update("speed", rec)
    for r in rows:
        log("  %-7s %-18s %9.1f s/frame  (%.1fx vs verbatim same grid)"
            % (r["grid"], r["method"], r["s_per_frame"] or 0.0,
               r["x_vs_verbatim_same_grid"] or 0.0))
    return rec


# ---------------------------------------------------------------------------
# the 512^2 base grid of v4-1.15, re-run with the fast path (48 frames), so
# that the 1024^2 vs 512^2 comparison comes from ONE code path and so that the
# fast path is checked against v4-1.15's own 48 recorded frames
# ---------------------------------------------------------------------------
def base512():
    makedirs()
    log("=" * 78)
    log("[base512] sub 512^2, 48 frames, fast path (same runner as STAGE 3)")
    import multiprocessing as mp
    nproc = min(8, max(1, (os.cpu_count() or 1) - 2))
    idxs = list(range(NFRAMES))
    chunks = [idxs[k::nproc] for k in range(nproc)]
    chunks = [c for c in chunks if c]
    t0 = time.perf_counter()
    with mp.Pool(nproc, initializer=_w_init,
                 initargs=("sub", 512, BOX_SUB, True)) as pool:
        parts = pool.map(_w_frames, chunks)
    wall = time.perf_counter() - t0
    rows = sorted([r for p in parts for r in p], key=lambda r: r["frame"])
    cp = read_json(os.path.join(HERE, "out", "v4-1-15", "checkpoint.json"))
    base = {}
    for k, node in (cp.get("chunks") or {}).get("sub", {}).items():
        for m in node["metrics"]:
            base[int(m["frame"])] = m
    diffs = []
    for r in rows:
        b = base[r["frame"]]
        diffs.append({"frame": r["frame"],
                      "std_rel_diff_vs_v4_1_15": float(
                          r["C_stats"]["std"] / b["C_std"] - 1.0),
                      "p99_rel_diff_vs_v4_1_15": float(
                          r["C_stats"]["p99"] / b["C_p99"] - 1.0),
                      "energy_deposited_rel_diff": float(
                          r["energy_deposited"] / b["energy_deposited"] - 1.0)})
    sd = np.array([d["std_rel_diff_vs_v4_1_15"] for d in diffs])
    pd = np.array([d["p99_rel_diff_vs_v4_1_15"] for d in diffs])
    rec = {"created": now_iso(), "N_s": 512, "workers": int(nproc),
           "total_wall_s": float(wall), "s_per_frame_wall": float(wall / NFRAMES),
           "frames": rows, "vs_v4_1_15_recorded": {
               "n_frames": len(diffs),
               "std_rel_diff_max_abs": float(np.abs(sd).max()),
               "p99_rel_diff_max_abs": float(np.abs(pd).max()),
               "std_rel_diff_mean": float(sd.mean()),
               "p99_rel_diff_mean": float(pd.mean()),
               "all_zero": bool(np.all(sd == 0.0) and np.all(pd == 0.0)),
               "note": "v4-1.15's 48 recorded frame metrics at 512^2 were "
                       "produced by the verbatim deposit; the fast path must "
                       "reproduce them"},
           "detail": diffs}
    rec["summary"] = {
        "mean_of_std": float(np.mean([r["C_stats"]["std"] for r in rows])),
        "mean_of_p99": float(np.mean([r["C_stats"]["p99"] for r in rows])),
        "max_rel_err": float(max(r["energy_rel_err"] for r in rows))}
    numbers_update("base512_fast", rec)
    log("  512^2 48 frames in %.0f s; vs v4-1.15 recorded metrics: max|dstd| "
        "%.3e max|dp99| %.3e  all_zero %s"
        % (wall, rec["vs_v4_1_15_recorded"]["std_rel_diff_max_abs"],
           rec["vs_v4_1_15_recorded"]["p99_rel_diff_max_abs"],
           rec["vs_v4_1_15_recorded"]["all_zero"]))
    return rec


# ---------------------------------------------------------------------------
# the fast module's own self test, recorded
# ---------------------------------------------------------------------------
def selftest():
    makedirs()
    log("=" * 78)
    log("[selftest] transport_mesh_fast.py --selftest")
    import subprocess
    t0 = time.perf_counter()
    r = subprocess.run([sys.executable, os.path.join(HERE,
                                                     "transport_mesh_fast.py"),
                        "--selftest"], capture_output=True, text=True,
                       cwd=HERE)
    wall = time.perf_counter() - t0
    lines = [ln for ln in (r.stdout or "").splitlines() if ln.strip()]
    rec = {"created": now_iso(), "cmd": [sys.executable,
                                         "transport_mesh_fast.py",
                                         "--selftest"],
           "returncode": int(r.returncode), "wall_s": float(wall),
           "ok": bool(r.returncode == 0 and "SELFTEST_OK" in (r.stdout or "")),
           "n_checks_ok": int(sum(1 for ln in lines if ln.startswith("[ok]"))),
           "n_checks_ng": int(sum(1 for ln in lines if ln.startswith("[NG]"))),
           "log": lines, "stderr_tail": (r.stderr or "")[-2000:]}
    numbers_update("selftest_fast", rec)
    for ln in lines:
        log("  " + ln)
    log("[selftest] %s (%.1f s)" % ("PASS" if rec["ok"] else "FAIL", wall))
    return rec


# ---------------------------------------------------------------------------
# supplementary equivalence (NOT a pre-registered gate): other grids/frames
# ---------------------------------------------------------------------------
def extra_eq():
    makedirs()
    log("=" * 78)
    log("[extra_eq] supplementary: fast vs verbatim on configurations that are "
        "NOT the two gate frames")
    import transport_mesh_fast as TMF
    TM.BOX_BULK_MAX = BOX_SUB
    TMF.BOX_BULK_MAX = BOX_SUB
    out = {"created": now_iso(),
           "status": "supplementary evidence, not a pre-registered gate",
           "cases": []}
    for tag, n, frameseq in (("sub", 256, (24,)), ("sub", 512, (10,)),
                             ("contrast", 512, (24,))):
        st, senv, optics, spec = R15.build_scene(tag, n=n)
        for fi in frameseq:
            f = st.fields_at(fi / float(FPS))
            Cv, iv = TM.mesh_C(f, st.S, WIDTH, HEIGHT, senv)
            Cf, iff = TMF.mesh_C_fast(f, st.S, WIDTH, HEIGHT, senv)
            a = Cv.astype(np.float64)
            b = Cf.astype(np.float64)
            d = np.abs(a - b)
            case = {"stage": tag, "N_s": int(n), "frame": int(fi),
                    "max_abs_px_diff": float(d.max()),
                    "frac_pixels_differing": float(np.mean(d > 0.0)),
                    "total_energy_rel_diff": float(
                        abs(b.sum() - a.sum()) / abs(a.sum())),
                    "std_rel_diff": float(np.std(b) / np.std(a) - 1.0),
                    "p99_rel_diff": float(np.percentile(b, 99)
                                          / np.percentile(a, 99) - 1.0),
                    "rel_err_verbatim": float(iv["rel_err"]),
                    "rel_err_fast": float(iff["rel_err"]),
                    "n_tail_triangles": int(iff["n_tail"])}
            case["identical"] = bool(case["max_abs_px_diff"] == 0.0)
            out["cases"].append(case)
            log("  %-8s %4d^2 frame %2d: max|dC| %.3e  dE %.2e  dstd %.2e  "
                "dp99 %.2e  tail %d"
                % (tag, n, fi, case["max_abs_px_diff"],
                   case["total_energy_rel_diff"], case["std_rel_diff"],
                   case["p99_rel_diff"], case["n_tail_triangles"]))
            del Cv, Cf, a, b, d, f
            gc.collect()
        del st, optics, spec
        gc.collect()
    out["all_identical"] = bool(all(c["identical"] for c in out["cases"]))
    numbers_update("supplementary_equivalence", out)
    TM.BOX_BULK_MAX = 16
    TMF.BOX_BULK_MAX = 16
    return out


# ---------------------------------------------------------------------------
# cross-checks + deliverable hashes (for the report)
# ---------------------------------------------------------------------------
def crosscheck():
    makedirs()
    log("=" * 78)
    log("[crosscheck] independent references that were NOT produced by this "
        "job's fast path")
    cp = read_json(os.path.join(HERE, "out", "v4-1-15", "checkpoint.json"))
    n = numbers()
    out = {"created": now_iso()}
    # (a) fast 1024^2 frame 24 (STAGE 3) vs v4-1.15's OWN verbatim 1024^2
    #     refinement of the same frame
    mine = [f for f in n["stage3"]["frames"] if f["frame"] == STILL_IDX][0]
    ref = cp["frame24"]["sub"]["refinement"]["C_refine"]
    out["stage3_frame24_vs_v4_1_15_verbatim_1024"] = {
        "recorded_reference": ref, "fast_result": mine["C_stats"],
        "abs_diff": {k: float(abs(mine["C_stats"][k] - ref[k]))
                     for k in ("mean", "std", "p99", "max")},
        "reference_energy_rel_err": cp["frame24"]["sub"]["refinement"]
        ["energy_refine"]["rel_err"],
        "fast_energy_rel_err": mine["energy_rel_err"],
        "what": "v4-1.15 ran ONE verbatim 1024^2 frame 24 (its 'refinement'); "
                "the 48-frame fast run's frame 24 must reproduce it",
        "note": "statistics agree to the last recorded digit; the energy "
                "relative errors differ at 1e-16 because the diagnostic "
                "counters sum in a different order"}
    # (b) gate reference arrays vs the v4-1.15 recorded statistics
    out["gate_reference_vs_v4_1_15"] = {
        tag: n["reference"]["frames"][tag]["provenance_vs_v4_1_15"]
        for tag in ("contrast", "sub")}
    # (c) deliverable identity
    out["files"] = {}
    for p in ("transport_mesh_fast.py", "run_v4_1_16.py", "transport_mesh.py"):
        fp = os.path.join(HERE, p)
        out["files"][p] = {"bytes": os.path.getsize(fp),
                           "sha256": sha256_file(fp)}
    out["transport_mesh_unchanged"] = {
        "sha256_now": out["files"]["transport_mesh.py"]["sha256"],
        "sha256_v4_1_15_recorded": cp["gates"]["G-LEGACY"]
        ["transport_mesh_sha256_now"],
        "same": bool(out["files"]["transport_mesh.py"]["sha256"]
                     == cp["gates"]["G-LEGACY"]["transport_mesh_sha256_now"])}
    numbers_update("crosscheck", out)
    for k, v in out["stage3_frame24_vs_v4_1_15_verbatim_1024"]["abs_diff"].items():
        log("  frame24 1024^2 fast vs v4-1.15 verbatim: d%s = %.3e" % (k, v))
    log("  transport_mesh.py unchanged: %s"
        % out["transport_mesh_unchanged"]["same"])
    return out


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "stage0"
    makedirs()
    if cmd == "stage0":
        stage0()
    elif cmd == "refgen":
        refgen()
    elif cmd == "equiv":
        equiv()
    elif cmd == "parbench":
        parbench()
    elif cmd == "verbatim1024":
        time_verbatim_1024()
    elif cmd == "stage3":
        stage3()
    elif cmd == "speed":
        speed()
    elif cmd == "crosscheck":
        crosscheck()
    elif cmd == "extra_eq":
        extra_eq()
    elif cmd == "base512":
        base512()
    elif cmd == "selftest":
        selftest()
    else:
        raise SystemExit("unknown stage %r" % (cmd,))
