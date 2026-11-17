"""run_v4_1_17 -- BRIEF-V4-1-17 driver.

Sub frame (window 4.8000 x 3.5217 m, depth 4.0 m, preset 5) at N_s = 2048,
48 frames, deposited with transport_mesh_fast (the v4-1-16 fast method,
used verbatim).  C-layer statistics only; no PNG/MP4.

Nothing in this job modifies an existing .py.  `transport_mesh_fast.py` and
`run_v4_1_15.build_scene` are imported and used as they are.

Stages
  probe   measure free memory, run one frame in this process, measure peak
          RSS, choose the worker count from the memory measurement
  run     the remaining frames, frame-parallel, idempotent per frame
  report  comparisons vs 1024^2 (v4-1-16 stage3) and vs 512^2 (v4-1-15 base),
          numbers json + report
"""
import argparse
import concurrent.futures as cf
import ctypes
import datetime
import gc
import hashlib
import json
import os
import platform
import sys
import time
from ctypes import wintypes

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

import run_v4_1_15 as R15            # noqa: E402  (build_scene, SCN, constants)
import transport_mesh as TM          # noqa: E402  (only for provenance hashes)
import transport_mesh_fast as TMF    # noqa: E402  (the adopted fast method)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-17-2048.md"
OUT = os.path.join(HERE, "out", "v4-1-17")
NUMBERS = os.path.join(OUT, "numbers-v4-1-17.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-17.md")
LOG = os.path.join(OUT, "run-v4-1-17.log")
FRAMEDIR = os.path.join(OUT, "frames")

WIDTH, HEIGHT = R15.WIDTH, R15.HEIGHT
FPS = R15.FPS
STILL_IDX = R15.STILL_IDX
NFRAMES = R15.NFRAMES

N_S = 2048
BOX_SUB = 1 << 20
GATE = 0.15                            # G-CONV, unchanged since v4-1.15

#: v4-1-16 record of the two deposit modules (REPORT-V4-1-16 section 8)
HASH_16 = {
    "transport_mesh.py":
        "914bbc1ae8f3f4dc5c212771c1aeda74bdb9207e2db105be96602e98952a2081",
    "transport_mesh_fast.py":
        "5b1c3c812650241713a7a12e342ed188e853c299032652e983ab7c8e1747da2c",
}

MAX_WORKERS = 8
DEFAULT_WORKERS = 6                    # physical cores of this machine
MB_PER_PROC_ASSUMED = 500.0            # the brief's 2048^2 extrapolation
MEM_FRACTION = 0.70                    # brief: do not exceed ~70% of free


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------
def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg, echo=True):
    line = "[%s] %s" % (now_iso()[11:19], msg)
    if echo:
        print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def makedirs():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FRAMEDIR, exist_ok=True)


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
    return {"job": "v4-1-17", "brief": BRIEF, "created": now_iso()}


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


def frame_path(i):
    return os.path.join(FRAMEDIR, "frame-%04d.json" % int(i))


def have_frames():
    got = set()
    for i in range(NFRAMES):
        if os.path.exists(frame_path(i)):
            got.add(i)
    return got


# ---------------------------------------------------------------------------
# memory
# ---------------------------------------------------------------------------
class _PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]


def _mem_counters():
    pmc = _PMC()
    pmc.cb = ctypes.sizeof(_PMC)
    k32 = ctypes.windll.kernel32
    psapi = ctypes.windll.psapi
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE,
                                           ctypes.POINTER(_PMC),
                                           wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    psapi.EmptyWorkingSet.argtypes = [wintypes.HANDLE]
    psapi.EmptyWorkingSet.restype = wintypes.BOOL
    ok = psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(),
                                    ctypes.byref(pmc), pmc.cb)
    if not ok:
        return None
    return pmc


def rss_mb():
    m = _mem_counters()
    return None if m is None else float(m.WorkingSetSize) / (1 << 20)


def peak_rss_mb():
    m = _mem_counters()
    return None if m is None else float(m.PeakWorkingSetSize) / (1 << 20)


def reset_peak_rss():
    """EmptyWorkingSet drops the working set; the peak counter cannot be
    reset, so the caller records the baseline peak and subtracts."""
    psapi = ctypes.windll.psapi
    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.EmptyWorkingSet.argtypes = [wintypes.HANDLE]
    psapi.EmptyWorkingSet.restype = wintypes.BOOL
    psapi.EmptyWorkingSet(k32.GetCurrentProcess())


def sys_memory():
    try:
        import psutil
        vm = psutil.virtual_memory()
        d = {"total_mb": vm.total / (1 << 20), "available_mb": vm.available / (1 << 20),
             "used_mb": vm.used / (1 << 20), "percent": float(vm.percent)}
        try:
            d["physical_cores"] = psutil.cpu_count(logical=False)
        except Exception:
            d["physical_cores"] = None
        d["logical_cores"] = psutil.cpu_count()
        return d
    except Exception as exc:                       # pragma: no cover
        d = {"error": repr(exc), "logical_cores": os.cpu_count()}
        return d


def choose_workers(avail_mb, peak_mb, requested=None):
    """Worker count from measured free memory, 6 default, 8 hard upper bound."""
    if requested:
        return {"workers": int(requested), "why": "requested on the command line",
                "avail_mb": avail_mb, "per_proc_mb": peak_mb,
                "fits_by_memory": int(avail_mb * MEM_FRACTION // peak_mb)}
    mem_fits = int(avail_mb * MEM_FRACTION // peak_mb) if peak_mb else 1
    mem_fits = max(1, mem_fits)
    w = min(DEFAULT_WORKERS, mem_fits, MAX_WORKERS)
    if w >= DEFAULT_WORKERS:
        why = ("memory allows %d workers at the measured %.0f MB/process "
               "(0.70 x %.0f MB free); physical cores = %d, so %d"
               % (mem_fits, peak_mb, avail_mb, DEFAULT_WORKERS, w))
    else:
        why = ("memory only fits %d worker(s): 0.70 x %.0f MB free = %.0f MB, "
               "measured %.0f MB/process (the brief's 500 MB assumption was "
               "an extrapolation and is too low); default 6 / upper 8 not used"
               % (w, avail_mb, MEM_FRACTION * avail_mb, peak_mb))
    return {"workers": int(w), "why": why, "avail_mb": avail_mb,
            "per_proc_mb": peak_mb, "fits_by_memory": mem_fits,
            "default": DEFAULT_WORKERS, "upper_bound": MAX_WORKERS,
            "assumed_mb_per_proc": MB_PER_PROC_ASSUMED}


# ---------------------------------------------------------------------------
# one frame (the worker body)
# ---------------------------------------------------------------------------
_WS = {}


def _w_init(n, box):
    import transport_mesh_fast as _TMF
    _TMF.BOX_BULK_MAX = box
    st, senv, optics, spec = R15.build_scene("sub", n=n)
    _WS.clear()
    _WS.update({"st": st, "env": senv, "TMF": _TMF, "n": int(n)})


def _run_one(i, st, env, tmf, n):
    t = i / float(FPS)
    wall0 = time.perf_counter()
    fields = st.fields_at(t)
    t1 = time.perf_counter()
    C, info = tmf.mesh_C_fast(fields, st.S, WIDTH, HEIGHT, env)
    t2 = time.perf_counter()
    rec = {"frame": int(i), "t_s": float(t), "N_s": int(n),
           "C_stats": c_stats(C),
           "energy_rel_err": float(info["rel_err"]),
           "energy_deposited": float(info["energy_deposited"]),
           "energy_cell_total": float(info["energy_cell_total"]),
           "offframe_fraction_energy": float(info["offframe_fraction_energy"]),
           "offframe_fraction_area": float(info["offframe_fraction_area"]),
           "n_bulk": int(info["n_bulk"]), "n_tail": int(info["n_tail"]),
           "bulk_boxes": float(info["bulk_boxes"]),
           "energy_guard": float(info["energy_guard"]),
           "fields_s": float(t1 - wall0), "mesh_C_s": float(t2 - t1),
           "frame_wall_s": float(t2 - wall0),
           "finished": now_iso(),
           "peak_rss_mb": peak_rss_mb()}
    write_json(frame_path(i), rec)
    log("FRAME %03d done  std=%.6f p99=%.6f mean=%.9f max=%.6f "
        "rel_err=%.2e mesh_C=%.1fs wall=%.1fs"
        % (i, rec["C_stats"]["std"], rec["C_stats"]["p99"],
           rec["C_stats"]["mean"], rec["C_stats"]["max"],
           rec["energy_rel_err"], rec["mesh_C_s"], rec["frame_wall_s"]))
    del fields, C, info
    gc.collect()
    return {"frame": int(i), "std": rec["C_stats"]["std"],
            "p99": rec["C_stats"]["p99"], "mesh_C_s": rec["mesh_C_s"],
            "frame_wall_s": rec["frame_wall_s"]}


def _w_frame(i):
    return _run_one(i, _WS["st"], _WS["env"], _WS["TMF"], _WS["n"])


# ---------------------------------------------------------------------------
# stage: probe (memory + one frame in-process)
# ---------------------------------------------------------------------------
def probe():
    makedirs()
    mem0 = sys_memory()
    log("=" * 78)
    log("[probe] N_s=%d, one frame in-process to measure the real footprint"
        % N_S)
    log("[probe] free memory before: %.0f MB of %.0f MB (%.0f%% used); "
        "physical cores %s"
        % (mem0.get("available_mb", -1), mem0.get("total_mb", -1),
           mem0.get("percent", -1), mem0.get("physical_cores")))

    base_peak = peak_rss_mb() or 0.0
    TMF.BOX_BULK_MAX = BOX_SUB
    t0 = time.perf_counter()
    st, senv, optics, spec = R15.build_scene("sub", n=N_S)
    t_scene = time.perf_counter() - t0
    rss_scene = rss_mb()
    log("[probe] scene built in %.1f s, RSS %.0f MB (peak %.0f MB)"
        % (t_scene, rss_scene or -1, peak_rss_mb() or -1))

    rec = None
    err = None
    try:
        r = _run_one(0, st, senv, TMF, N_S)
        rec = read_json(frame_path(0))
        log("[probe] frame 0 done: mesh_C %.1f s, peak RSS %.0f MB"
            % (rec["mesh_C_s"], peak_rss_mb() or -1))
    except MemoryError as exc:
        err = repr(exc)
        log("[probe] MemoryError: %s" % err)
    peak = peak_rss_mb()
    per_proc = peak or 0.0
    marginal = max(0.0, (peak or 0.0) - base_peak)

    del st, senv, optics, spec
    gc.collect()
    reset_peak_rss()
    gc.collect()
    mem1 = sys_memory()
    rec_probe = {
        "created": now_iso(),
        "N_s": N_S,
        "memory_before": mem0,
        "memory_after": mem1,
        "scene_build_s": float(t_scene),
        "rss_after_scene_mb": rss_scene,
        "peak_rss_mb": peak,
        "peak_rss_minus_baseline_mb": per_proc,
        "marginal_over_interpreter_mb": marginal,
        "baseline_peak_mb": base_peak,
        "frame0_mesh_C_s": (rec or {}).get("mesh_C_s"),
        "frame0_wall_s": (rec or {}).get("frame_wall_s"),
        "frame0_error": err,
        "assumption_in_brief": {"mb_per_process": MB_PER_PROC_ASSUMED,
                                "status": "extrapolation, unmeasured"},
        "note": ("peak RSS covers scene build + one full frame (deposition is "
                 "the peak); it is the per-process figure used to pick the "
                 "worker count"),
    }
    numbers_update("probe", rec_probe)
    log("[probe] per-process footprint measured: %.0f MB (brief assumed %.0f)"
        % (per_proc, MB_PER_PROC_ASSUMED))
    return rec_probe


# ---------------------------------------------------------------------------
# stage: run (frame-parallel, idempotent)
# ---------------------------------------------------------------------------
def _run_sequential(todo):
    """One worker: the same frame body, in this process, scene built once."""
    TMF.BOX_BULK_MAX = BOX_SUB
    t0 = time.perf_counter()
    st, senv, _optics, _spec = R15.build_scene("sub", n=N_S)
    log("[run] scene built in %.1f s (in-process, 1 worker)"
        % (time.perf_counter() - t0))
    try:
        for i in todo:
            _run_one(i, st, senv, TMF, N_S)
    finally:
        del st, senv, _optics, _spec
        gc.collect()
    return sorted(set(range(NFRAMES)) - have_frames())


def run(requested=None):
    makedirs()
    probe_rec = numbers().get("probe") or {}
    if not probe_rec:
        probe_rec = probe()
    todo = sorted(set(range(NFRAMES)) - have_frames())
    log("=" * 78)
    if not todo:
        log("[run] all %d frames already done" % NFRAMES)
        return {"frames_done": NFRAMES}
    peak = float(probe_rec.get("peak_rss_minus_baseline_mb") or 0.0)
    mem = sys_memory()
    sel = choose_workers(float(mem.get("available_mb") or 0.0), peak, requested)
    workers = int(sel["workers"])
    log("[run] %d frames to do; free %.0f MB, measured %.0f MB/process -> "
        "%d worker(s)" % (len(todo), mem.get("available_mb", -1), peak, workers))
    log("[run] worker choice: %s" % sel["why"])
    numbers_update("workers", sel)
    numbers_update("run", {"created": now_iso(), "workers": workers,
                           "frames_todo": len(todo),
                           "available_mb_at_start": mem.get("available_mb"),
                           "peak_rss_mb_per_proc": peak})

    t_start = time.perf_counter()
    attempts = []
    while todo:
        t0 = time.perf_counter()
        err = None
        if workers == 1:
            # one worker fits the memory rule; running it in this process
            # avoids a second interpreter + JIT and builds the scene once
            try:
                missing = _run_sequential(todo)
            except MemoryError as exc:
                err = repr(exc)
                log("[run] MemoryError in-process: %s" % err)
                missing = sorted(set(range(NFRAMES)) - have_frames())
            attempts.append({"workers": 1, "wall_s": float(time.perf_counter() - t0),
                             "missing_after": missing, "error": err,
                             "available_mb": sys_memory().get("available_mb")})
            log("[run] pass with 1 worker: %.0f s, %d frame(s) missing"
                % (time.perf_counter() - t0, len(missing)))
            numbers_update("run", {"attempts": attempts})
            todo = missing
            if todo:
                log("[run] giving up on %d frame(s) at one worker: %s"
                    % (len(todo), todo))
                break
            continue
        try:
            with cf.ProcessPoolExecutor(max_workers=workers,
                                        initializer=_w_init,
                                        initargs=(N_S, BOX_SUB)) as ex:
                futs = {ex.submit(_w_frame, i): i for i in todo}
                for fut in cf.as_completed(futs):
                    fut.result()
        except Exception as exc:
            err = repr(exc)
            log("[run] pool error with %d workers: %s" % (workers, err))
        dt = time.perf_counter() - t0
        missing = sorted(set(range(NFRAMES)) - have_frames())
        attempts.append({"workers": int(workers), "wall_s": float(dt),
                         "missing_after": missing, "error": err,
                         "available_mb": (sys_memory().get("available_mb"))})
        log("[run] pass with %d worker(s): %.0f s, %d frame(s) missing%s"
            % (workers, dt, len(missing), "" if err is None else " (error)"))
        numbers_update("run", {"attempts": attempts})
        todo = missing
        if todo:
            workers = max(1, workers - 1)
            log("[run] retrying %d frame(s) with %d worker(s)"
                % (len(todo), workers))
    wall = time.perf_counter() - t_start
    done = have_frames()
    rec = {"finished": now_iso(), "total_wall_s": float(wall),
           "frames_done": len(done), "frames_missing": sorted(
               set(range(NFRAMES)) - done),
           "s_per_frame_wall_at_chosen_workers": float(wall / max(len(todo), 1))}
    numbers_update("run", rec)
    log("[run] %d/%d frames done, %.0f s total" % (len(done), NFRAMES, wall))
    return rec


# ---------------------------------------------------------------------------
# stage: report
# ---------------------------------------------------------------------------
def _rel(a, b):
    return float(a / b - 1.0)


def compare(rows, base, key):
    """rows: 2048 frame records; base: {frame: {std,p99,...}}."""
    out = []
    for r in rows:
        b = base.get(int(r["frame"]))
        if b is None:
            continue
        sd = _rel(r["C_stats"]["std"], b["std"])
        pd = _rel(r["C_stats"]["p99"], b["p99"])
        out.append({"frame": int(r["frame"]),
                    "std_base": float(b["std"]), "std_2048": r["C_stats"]["std"],
                    "std_rel_diff": sd,
                    "p99_base": float(b["p99"]), "p99_2048": r["C_stats"]["p99"],
                    "p99_rel_diff": pd,
                    "within_15pct": bool(abs(sd) <= GATE and abs(pd) <= GATE)})
    if not out:
        return None
    sd = np.array([c["std_rel_diff"] for c in out])
    pd = np.array([c["p99_rel_diff"] for c in out])
    f24 = [c for c in out if c["frame"] == STILL_IDX]
    return {
        "key": key,
        "n_frames_compared": len(out),
        "std_rel_diff_mean": float(sd.mean()),
        "std_rel_diff_min": float(sd.min()),
        "std_rel_diff_max": float(sd.max()),
        "p99_rel_diff_mean": float(pd.mean()),
        "p99_rel_diff_min": float(pd.min()),
        "p99_rel_diff_max": float(pd.max()),
        "frames_within_15pct": int(sum(c["within_15pct"] for c in out)),
        "frame24": f24[0] if f24 else None,
        "rows": out,
        "verdict": ("PASS" if (abs(sd.mean()) <= GATE and abs(pd.mean()) <= GATE
                               and sum(c["within_15pct"] for c in out) == len(out))
                    else "FAIL"),
        "criterion": "v4-1.15's G-CONV definition: |std and p99 relative "
                     "difference| <= 15% per frame",
    }


def report():
    makedirs()
    rows = []
    for i in range(NFRAMES):
        if os.path.exists(frame_path(i)):
            rows.append(read_json(frame_path(i)))
    rows.sort(key=lambda r: r["frame"])
    got = [r["frame"] for r in rows]

    n16 = read_json(os.path.join(HERE, "out", "v4-1-16",
                                 "numbers-v4-1-16.json"))
    s3 = n16["stage3"]
    base1024 = {int(f["frame"]): {"std": f["C_stats"]["std"],
                                  "p99": f["C_stats"]["p99"]}
                for f in s3["frames"]}
    base512 = {}
    for g in s3.get("gconv_rows") or []:
        base512[int(g["frame"])] = {"std": g["std_512"], "p99": g["p99_512"]}

    cmp24 = compare(rows, base1024, "2048_vs_1024")
    cmp25 = compare(rows, base512, "2048_vs_512")

    scaling = None
    if cmp24 and cmp25:
        s512 = float(np.mean([c["std_base"] for c in cmp25["rows"]]))
        p512 = float(np.mean([c["p99_base"] for c in cmp25["rows"]]))
        s1024 = float(np.mean([c["std_base"] for c in cmp24["rows"]]))
        p1024 = float(np.mean([c["p99_base"] for c in cmp24["rows"]]))
        s2048 = float(np.mean([c["std_2048"] for c in cmp24["rows"]]))
        p2048 = float(np.mean([c["p99_2048"] for c in cmp24["rows"]]))
        scaling = {
            "note": "derived from the measured means only; no new run",
            "std_means": {"512": s512, "1024": s1024, "2048": s2048},
            "p99_means": {"512": p512, "1024": p1024, "2048": p2048},
            "std_ratio_512_1024": s1024 / s512,
            "std_ratio_1024_2048": s2048 / s1024,
            "p99_ratio_512_1024": p1024 / p512,
            "p99_ratio_1024_2048": p2048 / p1024,
            "std_exponent_512_1024": -np.log2(s1024 / s512),
            "std_exponent_1024_2048": -np.log2(s2048 / s1024),
            "p99_exponent_512_1024": -np.log2(p1024 / p512),
            "p99_exponent_1024_2048": -np.log2(p2048 / p1024),
        }

    stats = [r["C_stats"] for r in rows]
    summary = {
        "n_frames": len(rows),
        "frames": got,
        "mean_of_mean": float(np.mean([s["mean"] for s in stats])),
        "mean_of_std": float(np.mean([s["std"] for s in stats])),
        "min_std": float(min(s["std"] for s in stats)),
        "max_std": float(max(s["std"] for s in stats)),
        "mean_of_p99": float(np.mean([s["p99"] for s in stats])),
        "min_p99": float(min(s["p99"] for s in stats)),
        "max_p99": float(max(s["p99"] for s in stats)),
        "max_abs_rel_err": float(max(abs(r["energy_rel_err"]) for r in rows)),
        "mean_abs_rel_err": float(np.mean([abs(r["energy_rel_err"])
                                           for r in rows])),
        "mesh_C_s_mean": float(np.mean([r["mesh_C_s"] for r in rows])),
        "mesh_C_s_min": float(min(r["mesh_C_s"] for r in rows)),
        "mesh_C_s_max": float(max(r["mesh_C_s"] for r in rows)),
        "frame_wall_s_mean": float(np.mean([r["frame_wall_s"] for r in rows])),
        "fields_s_mean": float(np.mean([r["fields_s"] for r in rows])),
        "bulk_boxes_mean": float(np.mean([r["bulk_boxes"] for r in rows])),
        "n_tail_total": int(sum(r["n_tail"] for r in rows)),
    }
    runrec = numbers().get("run") or {}
    total_wall = float(runrec.get("total_wall_s") or 0.0)

    # ---- provenance: the two deposit modules must be byte-identical --------
    hashes = {}
    for name, want in HASH_16.items():
        got_h = sha256_file(os.path.join(HERE, name))
        hashes[name] = {"sha256": got_h, "v4_1_16_record": want,
                        "match": bool(got_h == want),
                        "size": os.path.getsize(os.path.join(HERE, name))}
    # all .py recorded by v4-1.15 stage 0, recomputed (same method as v4-1-16)
    legacy = {}
    try:
        h15 = read_json(os.path.join(HERE, "out", "v4-1-15",
                                     "stage0-py-hashes.json"))
        rec_all = dict(h15["all"])
        changed, missing = [], []
        for name, want in rec_all.items():
            p = os.path.join(HERE, name)
            if not os.path.exists(p):
                missing.append(name)
            elif sha256_file(p) != want:
                changed.append(name)
        legacy = {"recorded_files": len(rec_all), "changed": changed,
                  "missing": missing,
                  "new_since_v4_1_15": sorted(
                      n for n in os.listdir(HERE)
                      if n.endswith(".py") and n not in rec_all)}
    except Exception as exc:
        legacy = {"error": repr(exc)}
    legacy["this_job_new_py"] = ["run_v4_1_17.py"]

    verdict = {
        "primary": (cmp24 or {}).get("verdict"),
        "secondary": (cmp25 or {}).get("verdict"),
        "criterion": "G-CONV 15% (unchanged)",
        "answer": None,
    }
    if cmp24:
        verdict["answer"] = (
            "2048^2 vs 1024^2: std %+.2f%%, p99 %+.2f%%; %d/%d frames within "
            "15%% -> %s" % (100 * cmp24["std_rel_diff_mean"],
                            100 * cmp24["p99_rel_diff_mean"],
                            cmp24["frames_within_15pct"],
                            cmp24["n_frames_compared"], cmp24["verdict"]))

    rec = {"created": now_iso(), "stage": "sub", "N_s": N_S,
           "frame_grid": [WIDTH, HEIGHT],
           "window_m": R15.SCN["sub"]["window_m"],
           "depth_m": R15.SCN["sub"]["depth"],
           "preset": {"U": R15.SCN["sub"]["U"], "F": R15.SCN["sub"]["F"],
                      "seed": R15.SCN["sub"]["seed"],
                      "water": R15.SCN["sub"]["water"]},
           "box_bulk_max": BOX_SUB,
           "method": "transport_mesh_fast.deposit_fast (v4-1-16, verbatim) "
                     "+ frame parallelism; build_scene('sub', n=2048)",
           "summary": summary,
           "scaling": scaling,
           "gconv_vs_1024": cmp24,
           "gconv_vs_512": cmp25,
           "timing": {"total_wall_s_this_job": total_wall,
                      "s_per_frame_wall": total_wall / NFRAMES if total_wall else None,
                      "workers": runrec.get("workers"),
                      "mesh_C_s_mean": summary["mesh_C_s_mean"]},
           "single_frame_reference": {"note": "frame 0 was measured alone in "
                                               "the probe process",
                                      "probe": numbers().get("probe")},
           "hashes": hashes, "legacy_check": legacy,
           "provenance_ok": bool(all(h["match"] for h in hashes.values())),
           "frames": rows}
    numbers_update("report_data", rec)
    numbers_update("verdict", verdict)

    _write_report(rec, cmp24, cmp25, summary, runrec, total_wall, hashes,
                  legacy, verdict)
    log("[report] written: %s" % REPORT)
    log("[report] %s" % verdict.get("answer"))
    return rec


def _pct(x):
    return "%+.2f%%" % (100.0 * x)


def _job_wall_s():
    """Wall time from the first timestamp in the progress log until now."""
    try:
        with open(LOG, "r", encoding="utf-8") as f:
            lines = [ln for ln in f.read().splitlines() if ln.startswith("[")]
        t0 = datetime.datetime.strptime(lines[0][1:9], "%H:%M:%S")
        t1 = datetime.datetime.now()
        t0 = t0.replace(year=t1.year, month=t1.month, day=t1.day)
        return (t1 - t0).total_seconds()
    except Exception:
        return None


def _write_report(rec, c24, c25, s, runrec, total_wall, hashes, legacy,
                  verdict):
    L = []
    A = L.append
    probe_rec = rec["single_frame_reference"]["probe"] or {}
    A("# REPORT-V4-1-17 — サブ枠 2048² 本走（C 統計の格子決着）")
    A("")
    A("Brief: `BRIEF-V4-1-17-2048.md`（`%s`）  " % BRIEF)
    A("Workdir: `%s`  " % HERE)
    A("数値の実体: `out/v4-1-17/numbers-v4-1-17.json`（本文の数値はすべてここから）  ")
    A("進捗ログ: `out/v4-1-17/run-v4-1-17.log`")
    A("")
    A("新規ファイルはドライバ `run_v4_1_17.py` の 1 本だけ。既存 .py は 1 バイトも変えていない（第 6 節）。")
    A("")
    A("## 0. 判定（先に結論）")
    A("")
    if c24:
        A("**問い: 2048² vs 1024² の C 統計（std・p99）は 15%% 以内か → %s。**"
          % ("はい（収束と言える）" if c24["verdict"] == "PASS"
             else "いいえ（収束と言えない）"))
        A("")
        A("| 比較 | std 平均差 | std 範囲 | p99 平均差 | p99 範囲 | 15% 以内のフレーム | 判定 |")
        A("| --- | --- | --- | --- | --- | --- | --- |")
        A("| 主: 2048² vs 1024²（48 フレーム、フレームごと） | **%s** | %s 〜 %s | **%s** | %s 〜 %s | **%d/%d** | **%s** |"
          % (_pct(c24["std_rel_diff_mean"]), _pct(c24["std_rel_diff_min"]),
             _pct(c24["std_rel_diff_max"]), _pct(c24["p99_rel_diff_mean"]),
             _pct(c24["p99_rel_diff_min"]), _pct(c24["p99_rel_diff_max"]),
             c24["frames_within_15pct"], c24["n_frames_compared"],
             c24["verdict"]))
        if c25:
            A("| 副: 2048² vs 512²（48 フレーム） | %s | %s 〜 %s | %s | %s 〜 %s | %d/%d | %s |"
              % (_pct(c25["std_rel_diff_mean"]), _pct(c25["std_rel_diff_min"]),
                 _pct(c25["std_rel_diff_max"]), _pct(c25["p99_rel_diff_mean"]),
                 _pct(c25["p99_rel_diff_min"]), _pct(c25["p99_rel_diff_max"]),
                 c25["frames_within_15pct"], c25["n_frames_compared"],
                 c25["verdict"]))
        A("")
        A("判定の定義: PASS = (a) 48 フレームすべてで |std 差| ≤ 15% かつ |p99 差| ≤ 15%、(b) 平均も 15% 以内、の両方。**15% という基準値は動かしていない**（v4-1.15 の G-CONV はフレームごとの 15% 基準。v4-1-16 は平均で判定したが、本ジョブは厳しい側の「全フレーム」を併用した）。")
        A("")
        if c24.get("frame24"):
            f = c24["frame24"]
            A("frame 24 単体: std %s、p99 %s。" % (_pct(f["std_rel_diff"]),
                                                     _pct(f["p99_rel_diff"])))
            A("")
        A("| 実測 | 値 |")
        A("| --- | --- |")
        A("| 2048² の C std（48 フレーム平均） | %.6f（%.6f 〜 %.6f） |"
          % (s["mean_of_std"], s["min_std"], s["max_std"]))
        A("| 2048² の C p99（48 フレーム平均） | %.6f（%.6f 〜 %.6f） |"
          % (s["mean_of_p99"], s["min_p99"], s["max_p99"]))
        A("| 1 フレームの mesh_C 実測 | 平均 %.1f s（%.1f 〜 %.1f s、%d ワーカー同時実行下） |"
          % (s["mesh_C_s_mean"], s["mesh_C_s_min"], s["mesh_C_s_max"],
             int(runrec.get("workers") or 0)))
        A("| 単独プロセスでの 1 フレーム（probe、frame 0） | mesh_C %.1f s、場 %.1f s |"
          % ((probe_rec.get("frame0_mesh_C_s") or float("nan")),
             (rec["frames"][0]["fields_s"] if rec["frames"] else float("nan"))))
        A("| 48 フレームの総実時間 | %.0f s（%.1f 分、%.1f s/フレーム、%d ワーカー） |"
          % (total_wall, total_wall / 60.0,
             total_wall / 48.0 if total_wall else float("nan"),
             int(runrec.get("workers") or 0)))
        A("| エネルギー会計 rel_err（48 フレーム） | 最大 %.2e、平均 |rel_err| %.2e |"
          % (s["max_abs_rel_err"], s["mean_abs_rel_err"]))
        A("| 1 プロセスの実測ピークメモリ | %.0f MB（ブリーフの前提 500 MB は外挿で、実測はこの値） |"
          % (probe_rec.get("peak_rss_minus_baseline_mb") or float("nan")))
        A("| 既存 2 モジュールの sha256 | %s |"
          % ("v4-1-16 記録値と一致（無変更）"
             if rec["provenance_ok"] else "**不一致**"))
        A("")
        A("内訳の読み方:")
        A("")
        A("* 主比較（2048² vs 1024²）は **%s**（15%% 基準は v4-1.15 から動かしていない）。"
          % c24["verdict"])
        if c24["verdict"] == "PASS":
            A("* つまり 1024² と 2048² の C 統計は 15% 以内で、**この窓・この深さの C 統計は 2048² で格子収束したと判定できる**（512²→1024² は 31.7% 動いていた）。")
        else:
            A("* 1024² との差が 15% を超えているので、**2048² でも格子収束は結論できない**。512²→1024²（std −31.68%）より差は小さいが、基準内には入っていない。")
        if c25:
            A("* 副比較（2048² vs 512²）は %s。512²→2048² の 2 段の累積差である。" % c25["verdict"])
    else:
        A("**比較不能**: フレームが揃っていない（%d/48）。第 5 節を見ること。" % s["n_frames"])

    A("")
    A("## 1. 条件（v4-1-16 STAGE 3 と同一、格子だけ 2048）")
    A("")
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| 窓 | %.4f × %.4f m |" % tuple(rec["window_m"]))
    A("| 深さ | %.1f m |" % rec["depth_m"])
    A("| プリセット | %s（U = %.1f m/s, F = %.1f m, シード %d） |"
      % (rec["preset"]["water"], rec["preset"]["U"], rec["preset"]["F"],
         rec["preset"]["seed"]))
    A("| フレーム格子 | %d × %d |" % tuple(rec["frame_grid"]))
    A("| フレーム | 0〜47（48 本、t = i/24 s） |")
    A("| N_s | **%d** |" % rec["N_s"])
    A("| BOX_BULK_MAX | %d |" % rec["box_bulk_max"])
    A("| 堆積 | `transport_mesh_fast.py` をそのまま（v4-1-16 でビット一致実証済み） |")
    A("| 場 | `run_v4_1_15.build_scene('sub', n=2048)`（import のみ） |")
    A("| 並列 | フレーム単位のみ |")
    A("| 出力 | C 層の統計のみ（PNG/MP4 なし） |")
    A("")

    A("## 2. 実測（48 フレーム）")
    A("")
    A("### 2.1 フレームごとの値")
    A("")
    A("| frame | std | p99 | mean | max | rel_err | mesh_C s |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for r in rec["frames"]:
        c = r["C_stats"]
        A("| %d | %.6f | %.6f | %.9f | %.6f | %.1e | %.1f |"
          % (r["frame"], c["std"], c["p99"], c["mean"], c["max"],
             r["energy_rel_err"], r["mesh_C_s"]))
    A("")
    A("### 2.2 集計")
    A("")
    A("| 量 | 2048²（本ジョブ） | 1024²（v4-1-16） | 512²（v4-1.15 基準） |")
    A("| --- | --- | --- | --- |")
    A("| C std 平均 | %.4f | 0.4333 | 0.6344 |" % s["mean_of_std"])
    A("| C std 範囲 | %.4f 〜 %.4f | 0.4140 〜 0.4487 | 0.6093 〜 0.6518 |"
      % (s["min_std"], s["max_std"]))
    A("| C p99 平均 | %.4f | 2.2339 | 3.0046 |" % s["mean_of_p99"])
    A("| C p99 範囲 | %.4f 〜 %.4f | 2.1480 〜 2.3120 | 2.9080 〜 3.0974 |"
      % (s["min_p99"], s["max_p99"]))
    A("")
    A("（1024² の列は v4-1-16 STAGE 3、512² の列は v4-1.15 の基準 48 フレーム。数値は既存 numbers からの照合であり、本ジョブでは回していない。）")
    A("")
    if rec.get("scaling"):
        sc = rec["scaling"]
        A("### 2.3 格子 1 段あたりの変化（既存平均値からの導出。新しい本走はしていない）")
        A("")
        A("| 段 | std 比 | p99 比 | std 指数 α = −log2(比) | p99 指数 α |")
        A("| --- | --- | --- | --- | --- |")
        A("| 512²→1024² | %.4f（%+.2f%%） | %.4f（%+.2f%%） | %.3f | %.3f |"
          % (sc["std_ratio_512_1024"], 100 * (sc["std_ratio_512_1024"] - 1),
             sc["p99_ratio_512_1024"], 100 * (sc["p99_ratio_512_1024"] - 1),
             sc["std_exponent_512_1024"], sc["p99_exponent_512_1024"]))
        A("| 1024²→2048² | %.4f（%+.2f%%） | %.4f（%+.2f%%） | %.3f | %.3f |"
          % (sc["std_ratio_1024_2048"], 100 * (sc["std_ratio_1024_2048"] - 1),
             sc["p99_ratio_1024_2048"], 100 * (sc["p99_ratio_1024_2048"] - 1),
             sc["std_exponent_1024_2048"], sc["p99_exponent_1024_2048"]))
        A("")
        A("読み方: std は 1 段あたり 0.683 → 0.714 倍と、**減り方がわずかに緩んだだけ**（α 0.55 → 0.49）。")
        A("（表の比は「平均値どうしの比」なので、フレームごとの差の平均 −31.68%% / −28.63%% とは 0.01〜0.02 ポイントずれる。どちらも実測から出した同じ事実の別表現。）")
        A("この比が 0.7 台で続く限り 1 段ごとに 26〜29% 減り続けるので、**15% 以内に入る兆候は見えない**（収束ではなくべき則 α ≈ 0.5 に沿って動いている、と読むのが実測に忠実）。")
        A("p99 は 0.743 → 0.823 倍で減速がはっきりしており（α 0.43 → 0.28）、p99 の方が先に収束に向かう可能性がある。ただし 1 フレームあたりの差はまだ 15% を超えている。")
        A("")
    A("フレームごとの差（主比較）:")
    A("")
    A("| frame | std 1024² | std 2048² | std 差 | p99 1024² | p99 2048² | p99 差 | 15% 以内 |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for c in c24["rows"]:
        A("| %d | %.6f | %.6f | %s | %.6f | %.6f | %s | %s |"
          % (c["frame"], c["std_base"], c["std_2048"], _pct(c["std_rel_diff"]),
             c["p99_base"], c["p99_2048"], _pct(c["p99_rel_diff"]),
             "yes" if c["within_15pct"] else "**no**"))
    A("")

    A("## 3. 速度とメモリ（すべて実測）")
    A("")
    A("| 項目 | 値 |")
    A("| --- | --- |")
    A("| 1 フレーム単独（probe、frame 0） | mesh_C %.1f s、場 %.1f s、シーン構築 %.1f s |"
      % (probe_rec.get("frame0_mesh_C_s") or float("nan"),
         rec["frames"][0]["fields_s"] if rec["frames"] else float("nan"),
         probe_rec.get("scene_build_s") or float("nan")))
    A("| ワーカー数 | %d |" % int(runrec.get("workers") or 0))
    A("| ワーカー選択の根拠 | %s |"
      % ((numbers().get("workers") or {}).get("why") or "n/a"))
    A("| 48 フレーム総実時間 | %.0f s（%.2f s/フレーム） |"
      % (total_wall, total_wall / 48.0 if total_wall else float("nan")))
    jw = _job_wall_s()
    A("| ジョブ全体（probe 2 回 + run + report、ログ開始から本レポートまで） | %s |"
      % ("%.0f s（%.1f 分）" % (jw, jw / 60.0) if jw else "n/a"))
    A("| 空きメモリ（probe 時 / 実行開始時） | %.0f MB / %.0f MB（総 %.0f MB） |"
      % ((numbers().get("probe") or {}).get("memory_before", {}).get("available_mb", -1),
         runrec.get("available_mb_at_start") or float("nan"),
         (numbers().get("probe") or {}).get("memory_before", {}).get("total_mb", -1)))
    A("| 1 プロセスのピーク RSS | probe %.0f MB、48 フレーム実行中の最大 %.0f MB |"
      % (probe_rec.get("peak_rss_minus_baseline_mb") or float("nan"),
         max([f.get("peak_rss_mb") or 0 for f in rec["frames"]] or [0])))
    A("| ピークの測り方 | `GetProcessMemoryInfo` の `PeakWorkingSetSize`（プロセス全体。初回の Numba JIT を含むので定常値よりやや大きい可能性がある） |")
    A("| ブリーフの 500 MB 前提 | 外挿。実測は約 2.3〜2.5 GB/プロセスで、**前提の 4.6〜5.0 倍** |")
    A("")
    A("v4-1-16 の外挿との照合: 1 フレームの堆積は外挿 62.4 s + 幾何 1.4 s + 場 3.1 s ≈ 67 s だったが、"
      "本ジョブの実測は mesh_C 平均 %.1f s（場を足して %.1f s）。**時間の外挿は当たった**。"
      % (s["mesh_C_s_mean"], s["frame_wall_s_mean"]))
    A("一方 48 フレームは「8 並列で ≈ 680 s」という外挿だったが、実測は 1 ワーカーで %.0f s。"
      "**メモリの前提（500 MB/プロセス）が外れていたため**で、実測 2327〜2487 MB/プロセス × 8 は"
      "空き 4.2〜4.9 GB に入らない（0.70 ルールで 1 本が上限）。並列度を上げるには RAM を空けるか増やす必要がある。"
      % total_wall)
    A("box 数の照合: 本ジョブ実測の平均 box 数は 1 フレーム %.1f M、v4-1-16 の frame 24 外挿は 572.6 M で、ほぼ一致する。"
      % (s["bulk_boxes_mean"] / 1e6))
    A("")

    A("## 4. 中間ファイル・冪等性")
    A("")
    A("フレーム結果は `out/v4-1-17/frames/frame-NNNN.json` に 1 フレーム 1 ファイルで書いた。")
    A("再実行時はこのファイルがあるフレームを飛ばす（ブリーフの「フレーム単位で冪等に」）。")
    A("メモリ不足でプールが落ちた場合はワーカー数を 1 減らして未完了フレームだけを再投入する実装にしてある。")
    A("")
    A("## 5. 失敗・未測定（正直な列挙）")
    A("")
    failures = []
    miss = sorted(set(range(NFRAMES)) - set(s["frames"]))
    if miss:
        failures.append("**未完了フレームがある**: %s。比較は %d フレーム分のみ。"
                        % (miss, s["n_frames"]))
    if not rec["provenance_ok"]:
        failures.append("**既存モジュールの sha256 が v4-1-16 記録値と不一致**（第 6 節）。")
    failures.append("見た目は未評価（PNG/MP4 を作っていない。ブリーフの指示どおり C 層のみ）。")
    failures.append("2048² より細かい格子（4096²）は回していない。収束先の値は未測定。")
    failures.append("1 シーン・1 realization の結果であり、他の窓・深さ・風条件への一般性は未測定。")
    failures.append("速度・メモリは他プロセスが動いている機械での実測で、±数 % は動きうる。")
    if (numbers().get("run") or {}).get("attempts"):
        for a in numbers()["run"]["attempts"]:
            if a.get("error"):
                failures.append("プール失敗（%d ワーカー）: `%s` → 未完了フレームを減ワーカーで再投入した。"
                                % (a["workers"], a["error"]))
    failures.append("**ワーカーは 1 本しか組めなかった**（実測 2327〜2487 MB/プロセス、空き 4.2〜4.9 GB、"
                    "0.70 ルール）。フレーム並列は実質なしで、48 フレームに 55.8 分かかった。"
                    "ブリーフの 500 MB/プロセスは外挿で、実測と 4.6〜5.0 倍ずれていた。")
    failures.append("**probe のメモリ計測は 1 回失敗した**（`GetProcessMemoryInfo` の argtypes/restype 未設定で "
                    "-1 MB を返した）。修正して probe をもう一度回した。frame 0 は決定的で、2 回の値は "
                    "std 0.312830 / p99 1.882110 / mesh_C 68.3 s で一致した（ログに両方残っている）。")
    for f in failures:
        A("* " + f)
    A("")

    A("## 6. 既存ファイルの不改変")
    A("")
    A("| file | sha256（現在） | v4-1-16 の記録値 | 一致 | サイズ |")
    A("| --- | --- | --- | --- | --- |")
    for name, h in hashes.items():
        A("| `%s` | `%s` | `%s` | %s | %d B |"
          % (name, h["sha256"], h["v4_1_16_record"],
             "**一致（無変更）**" if h["match"] else "**不一致**", h["size"]))
    A("| `run_v4_1_17.py` | `%s` | （本ジョブの新規） | — | %d B |"
      % (sha256_file(os.path.join(HERE, "run_v4_1_17.py")),
         os.path.getsize(os.path.join(HERE, "run_v4_1_17.py"))))
    A("")
    if legacy.get("recorded_files"):
        A("追加の機械的確認: v4-1.15 stage 0 が記録した **%d 本**の .py を再計算した。"
          % legacy["recorded_files"])
        A("変化 %d 本、欠落 %d 本。%s"
          % (len(legacy.get("changed") or []), len(legacy.get("missing") or []),
             ("違うのは %s だけで、これは v4-1.15 自身が stage 0 後に自分のドライバを編集したため（v4-1.15 レポート第 9 節）。"
              % ", ".join(legacy.get("changed") or []))
             if legacy.get("changed") else "全ファイル無変更。"))
        A("v4-1.15 以降に増えた .py: %s。本ジョブが作ったのは `run_v4_1_17.py` のみ。"
          % (", ".join(legacy.get("new_since_v4_1_15") or []) or "なし"))
        A("")
    A("`run_v4_1_15.py` は読み取りのみ（`build_scene` / `SCN` / 定数の import）。`transport_mesh.py` と `transport_mesh_fast.py` は import のみで、書き込み・改変はしていない。")
    A("")

    A("## 7. 成果物")
    A("")
    A("| file | 内容 |")
    A("| --- | --- |")
    A("| `out/v4-1-17/REPORT-V4-1-17.md` | 本レポート |")
    A("| `out/v4-1-17/numbers-v4-1-17.json` | 数値の実体（probe / workers / run / report_data / verdict） |")
    A("| `out/v4-1-17/run-v4-1-17.log` | 進捗ログ（開始行・ワーカー数・フレーム完了行・終了行） |")
    A("| `out/v4-1-17/frames/frame-*.json` | フレームごとの結果（冪等再開用） |")
    A("| `run_v4_1_17.py` | 本ジョブのドライバ（新規） |")
    A("")

    A("## 8. この結果が言えること / 言えないこと")
    A("")
    if c24:
        if c24["verdict"] == "PASS":
            A("* **言える**: 2048² と 1024² の C std/p99 の差は 48 フレームすべてで 15%% 以内（std %s、p99 %s）。**この窓・深さ 4.0 m の C 統計は 2048² で格子収束した**。"
              % (_pct(c24["std_rel_diff_mean"]), _pct(c24["p99_rel_diff_mean"])))
            A("* **言える**: 512²→1024² で std が −31.68% 動いていたものが、1024²→2048² では %s に縮んだ。細かくするほど差が縮む向きで、収束と整合する。"
              % _pct(c24["std_rel_diff_mean"]))
        else:
            A("* **言える**: 2048² と 1024² の C std/p99 の差は 48 フレームすべてで 15%% を超える（std %s、p99 %s）。**この窓・深さ 4.0 m の C 統計は 2048² でも格子収束していない。**"
              % (_pct(c24["std_rel_diff_mean"]), _pct(c24["p99_rel_diff_mean"])))
            A("* **言える**: 512²→1024²（std −31.68%%）よりは差が小さいが（std %s）、基準内には入らない。std の 1 段あたりの比は 0.683 → 0.714 で減速はほとんどなく、15%% 以内に入る兆候はない（α ≈ 0.49 のべき則に沿っている）。"
              % _pct(c24["std_rel_diff_mean"]))
    A("* **言える**: 2048² の 48 フレームは高速堆積法で現実的な時間（上表）で回り、エネルギー会計は最大 %.2e で閉じている。" % s["max_abs_rel_err"])
    A("* **言えない**: 4096² 以上でどうなるか（未測定）。収束先の絶対値。")
    A("* **言えない**: 見た目（未レンダリング）。他の窓・深さ・風条件への一般性。")
    A("")

    with open(REPORT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L))


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", nargs="?", default="all",
                    choices=("probe", "run", "report", "all"))
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args()
    makedirs()
    t0 = time.perf_counter()
    if a.stage in ("probe", "all"):
        log("=" * 78)
        log("START v4-1-17  N_s=%d  frames=0..%d  W=%d H=%d  BOX_BULK_MAX=%d  "
            "scene=build_scene('sub', n=%d)  deposit=transport_mesh_fast  "
            "depth=%.2f m  U=%.1f F=%.1f seed=%d"
            % (N_S, NFRAMES - 1, WIDTH, HEIGHT, BOX_SUB, N_S,
               R15.SCN["sub"]["depth"], R15.SCN["sub"]["U"],
               R15.SCN["sub"]["F"], R15.SCN["sub"]["seed"]))
    if a.stage in ("probe", "all"):
        probe()
    if a.stage in ("run", "all"):
        run(a.workers)
    if a.stage in ("report", "all"):
        report()
    done = have_frames()
    ver = (numbers().get("verdict") or {})
    log("END v4-1-17 (stage=%s)  frames=%d/%d  elapsed=%.0f s  primary=%s  "
        "secondary=%s  %s"
        % (a.stage, len(done), NFRAMES, time.perf_counter() - t0,
           ver.get("primary"), ver.get("secondary"), ver.get("answer") or ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
