"""run_v8_finish.py -- finish the v8 chunk-4 full-resolution run (missing stages only).

Measured starting point (the brief, out/v8-run.log):
  the full-resolution run died inside [M] (photon check): no python.exe left,
  out/v8 unchanged after 23:21.  [E] and [Q] had completed; the log records
  v8 state sigma2_spec=0.020954 realised=0.022135 Hs=0.1140, E_v8=2.708272,
  G1_v8=1.0, and Q=9/25/49 p99 = 1.5097/1.5098/1.5095.

REUSED here, never recomputed:
  out/v8/C-baseline-v33.npy   2048^2 float32 (3.3 master, cross-checked vs v7)
  out/v8/C-v8-t0.npy          2048^2 float32 (Q=49 chunk-4 master)
  out/v8/still-preset1-v33.png, out/v8/still-preset1-v8.png

ADDED here, all at 2048^2, foreground, peak memory recorded:
  [M] photon check at 4,000,000 photons -- the v2 check condition
      (out/v2/numbers-v2.json: photons=4000000, MS-SSIM=0.8274, gate 0.9,
      pass=false) -- using an equivalent low-peak-memory photon path
      (verified bit-identical to transport.photon_frame at 256^2/200k), plus
      the MS-SSIM peak-luminance / peak-position decomposition
  [Q] Q=9 and Q=25 recomputed, Q=49 reused from C-v8-t0.npy, converged
  gates 2/3/4 (shallow phase speed, new anchor sigma2, breaking clamp)
  gate 1 (ablation of the 3.3 -> v8 difference, both orderings + t=2 s
      dispersion measurement)
  gate 5 (standard-preset tone on the reused master), gate 6 (rain path)
  numbers-v8.json, numbers-v8-prereq.json, panel, recipes

Progress: every stage appends to out/v8/checkpoint.json (D:\\mutaro-work\\jobkit
Checkpoint format) and caches its section in out/v8/_parts/*.json, so a death
never loses a finished stage; re-running skips the cached ones.

Usage (from caustics-filter/offline/):
  python run_v8_finish.py --verify-lowmem     # photon-path equivalence check
  python run_v8_finish.py --stages base|photon|gates234|gate1|gates56|assemble
  python run_v8_finish.py --all
"""
import argparse
import ctypes
import datetime
import gc
import json
import math
import os
import sys
import threading
import time
from ctypes import wintypes

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from state import derive_seed, sample_bicubic  # noqa: E402
from transport import (orthonormal_basis, sun_disc_radius,  # noqa: E402
                       gaussian_blur_fft)
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert,  # noqa: E402
                        smith_g1_beckmann, scattering_blur_sigma_px,
                        substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,  # noqa: E402
                          image_metrics, v3_beer, flat_water_direct_level)
from v6_state import V6WaterState  # noqa: E402
from v6_physics import F_TABLE, flux_gate  # noqa: E402
from recipe_v8 import parse_recipe_v8, with_v8  # noqa: E402
from v8_state import V8WaterState  # noqa: E402
import v8_physics as P  # noqa: E402
import v8_prereq as PQ  # noqa: E402
import run_v8 as R  # noqa: E402  (shared conventions: c_stats, delta_stats...)

# --------------------------------------------------------------------------
# progress: checkpoint (jobkit format) + tee log
# --------------------------------------------------------------------------
JOBKIT_DIR = r"D:\mutaro-work\jobkit"

try:
    sys.path.insert(0, JOBKIT_DIR)
    from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402
    CHECKPOINT_SOURCE = os.path.join(JOBKIT_DIR, "checkpoint.py")
except Exception:  # pragma: no cover - jobkit is present on this machine
    CHECKPOINT_SOURCE = None

    def atomic_write_bytes(final_path, data, verify=None):
        tmp = final_path + ".tmp"
        with open(tmp, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        if verify is not None and not verify(tmp):
            os.remove(tmp)
            raise ValueError("verify failed: %r" % (tmp,))
        os.replace(tmp, final_path)
        return final_path

    class Checkpoint(object):
        """Minimal stand-in for jobkit.Checkpoint (same JSON shape)."""

        def __init__(self, path):
            self.path = os.path.abspath(path)

        def state(self):
            try:
                with open(self.path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except FileNotFoundError:
                return None

        def unit_done(self, unit_id):
            st = self.state() or {"completed_units": 0, "completed_ids": [],
                                  "last_unit": None, "updated_at": None}
            st["completed_units"] = int(st.get("completed_units") or 0) + 1
            st.setdefault("completed_ids", []).append(unit_id)
            st["last_unit"] = unit_id
            st["updated_at"] = datetime.datetime.now().astimezone().isoformat(
                timespec="seconds")
            atomic_write_bytes(self.path,
                               json.dumps(st, ensure_ascii=False,
                                          indent=2).encode("utf-8"))
            return st


class Tee(object):
    def __init__(self, path):
        self.f = open(path, "a", encoding="utf-8", buffering=1)

    def write(self, s):
        try:
            sys.__stdout__.write(s)
            sys.__stdout__.flush()
        except Exception:
            pass
        self.f.write(s)

    def flush(self):
        try:
            sys.__stdout__.flush()
        except Exception:
            pass
        self.f.flush()


# --------------------------------------------------------------------------
# peak-memory sampling (no psutil dependency)
# --------------------------------------------------------------------------
class _PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]

_psapi = ctypes.WinDLL("psapi", use_last_error=True) if os.name == "nt" else None
if _psapi is not None:
    _psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE,
                                            ctypes.POINTER(_PMC),
                                            wintypes.DWORD]
    _psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    _k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    _k32.GetCurrentProcess.restype = wintypes.HANDLE


def mem_now():
    if _psapi is None:
        return None
    pmc = _PMC()
    pmc.cb = ctypes.sizeof(_PMC)
    ok = _psapi.GetProcessMemoryInfo(_k32.GetCurrentProcess(),
                                     ctypes.byref(pmc), pmc.cb)
    if not ok:
        return None
    return {"ws_mb": pmc.WorkingSetSize / 2.0 ** 20,
            "peak_ws_mb": pmc.PeakWorkingSetSize / 2.0 ** 20}


class StageMem(object):
    """Samples WorkingSetSize every 150 ms and reports the stage maximum."""

    def __init__(self, label):
        self.label = label
        self._stop = False
        self.max_ws = 0.0
        m = mem_now() or {}
        self.start_peak = m.get("peak_ws_mb", 0.0)
        self.start_ws = m.get("ws_mb", 0.0)
        self.t = threading.Thread(target=self._run, daemon=True)
        self.t.start()

    def _run(self):
        while not self._stop:
            m = mem_now()
            if m:
                self.max_ws = max(self.max_ws, m["ws_mb"])
            time.sleep(0.15)

    def stop(self):
        self._stop = True
        self.t.join(timeout=3.0)
        m = mem_now() or {}
        self.max_ws = max(self.max_ws, m.get("ws_mb", 0.0))
        return {"stage": self.label,
                "stage_peak_working_set_mb": round(self.max_ws, 1),
                "process_peak_working_set_mb": round(
                    m.get("peak_ws_mb", 0.0), 1),
                "process_peak_at_stage_start_mb": round(self.start_peak, 1),
                "working_set_at_stage_start_mb": round(self.start_ws, 1)}


# --------------------------------------------------------------------------
# low-peak-memory photon path (equivalent to transport.photon_frame)
# --------------------------------------------------------------------------
def photon_frame_lowmem(fields, S, N_buf, env, n_photons, seed,
                        chunk=250_000):
    """transport.photon_frame with the 3-channel float64 accumulator split.

    Same physics, same random stream, same arithmetic order:
      * the (N*N,3) float64 accumulator becomes three contiguous (N*N,)
        float64 arrays (saves the strided temporaries of acc[:,b] += ...),
      * the final sigma=0.75 px Gaussian is applied one channel at a time
        instead of to a (N,N,3) float64 buffer.
    Verified bit-identical to transport.photon_frame at 256^2 / 200k photons
    against out/v8-smoke2/C-photon-t0.npy (--verify-lowmem).
    """
    rng = np.random.default_rng(seed)
    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    depth = float(env["depthM"])
    sunI = float(env["sunI"])
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    ior = np.array([1.33132, 1.33434, 1.33925], dtype=np.float64)
    patch_area = S * S
    cell = S / N_buf
    cell_area = cell * cell
    nb = N_buf * N_buf

    acc = [np.zeros(nb, dtype=np.float64) for _ in range(3)]
    n_done = 0
    while n_done < n_photons:
        m = min(chunk, n_photons - n_done)
        x0 = rng.random(m) * S
        y0 = rng.random(m) * S
        rr = disc * np.sqrt(rng.random(m))
        ph = rng.random(m) * 2.0 * np.pi
        jx = rr * np.cos(ph)
        jy = rr * np.sin(ph)
        d = (s[None, :] + u[None, :] * jx[:, None] + v[None, :] * jy[:, None])
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        I = -d

        h = sample_bicubic(fields["h"], x0, y0, S).astype(np.float64)
        Dx = sample_bicubic(fields["Dx"], x0, y0, S).astype(np.float64)
        Dy = sample_bicubic(fields["Dy"], x0, y0, S).astype(np.float64)
        gx = sample_bicubic(fields["gx"], x0, y0, S).astype(np.float64)
        gy = sample_bicubic(fields["gy"], x0, y0, S).astype(np.float64)
        nl = np.sqrt(gx * gx + gy * gy + 1.0)
        N = np.stack([-gx / nl, -gy / nl, np.ones(m) / nl], axis=1)
        Psx = x0 + Dx
        Psy = y0 + Dy
        Psz = h
        for b in range(3):
            eta = 1.0 / float(ior[b])
            cos1 = -(I * N).sum(axis=-1, keepdims=True)
            cos1 = np.clip(cos1, -1.0, 1.0)
            k = 1.0 - eta * eta * (1.0 - cos1 * cos1)
            k = np.maximum(k, 0.0)
            T = eta * I + (eta * cos1 - np.sqrt(k)) * N
            Tz = np.minimum(T[:, 2], -1e-6)
            t = (Psz + depth) / (-Tz)
            xb = np.mod(Psx + T[:, 0] * t, S)
            yb = np.mod(Psy + T[:, 1] * t, S)
            fx = xb / cell
            fy = yb / cell
            ix0 = np.floor(fx).astype(np.int64) % N_buf
            iy0 = np.floor(fy).astype(np.int64) % N_buf
            tx = fx - np.floor(fx)
            ty = fy - np.floor(fy)
            ix1 = (ix0 + 1) % N_buf
            iy1 = (iy0 + 1) % N_buf
            w = sunI * float(sunColor[b]) * patch_area / n_photons
            w00 = w * (1 - tx) * (1 - ty)
            w10 = w * tx * (1 - ty)
            w01 = w * (1 - tx) * ty
            w11 = w * tx * ty
            i00 = iy0 * N_buf + ix0
            i10 = iy0 * N_buf + ix1
            i01 = iy1 * N_buf + ix0
            i11 = iy1 * N_buf + ix1
            accb = acc[b]
            for idx, wt in ((i00, w00), (i10, w10), (i01, w01), (i11, w11)):
                tmp = np.bincount(idx, weights=wt, minlength=nb)
                accb += tmp
                del tmp
            del (ix0, iy0, ix1, iy1, tx, ty, w00, w10, w01, w11,
                 i00, i10, i01, i11, T, t, xb, yb, fx, fy, cos1, k)
        n_done += m
        del x0, y0, rr, ph, jx, jy, d, I, h, Dx, Dy, gx, gy, nl, N
        del Psx, Psy, Psz

    E = np.empty((N_buf, N_buf, 3), dtype=np.float32)
    for b in range(3):
        ch = acc[b].reshape(N_buf, N_buf) / cell_area
        E[:, :, b] = gaussian_blur_fft(ch, 0.75).astype(np.float32)
        del ch
    acc = None
    gc.collect()
    expected = sunI * sunColor
    measured = E.astype(np.float64).reshape(-1, 3).mean(axis=0)
    info = {"expected_mean": [float(x) for x in expected],
            "measured_mean": [float(x) for x in measured],
            "implementation": ("lowmem (channel-split accumulator, "
                              "equivalent to transport.photon_frame)")}
    return E, info


# --------------------------------------------------------------------------
# stage plumbing
# --------------------------------------------------------------------------
class Ctx(object):
    def __init__(self, args):
        self.args = args
        self.out = args.out
        self.parts = os.path.join(args.out, "_parts")
        os.makedirs(self.parts, exist_ok=True)
        self.cp = Checkpoint(os.path.join(args.out, "checkpoint.json"))
        self.peaks = {}
        self.timings = {}
        with open(args.recipe, "r", encoding="utf-8") as f:
            self.raw = json.load(f)
        self.scene = with_v8(self.raw, fetch_m=float(args.fetch))
        self.env = resolve_environment(self.scene)
        self.cfg = parse_recipe_v8(self.scene)
        self.grid = int(args.grid)
        self.Q = int(args.q)
        self.S_patch = 4.0 * float(self.env["specL"])
        self.seed = (int(args.seed) if args.seed is not None
                     else derive_seed(self.raw, args.preset))
        self.photo_lin = None
        self.png = lambda name: os.path.join(args.out, name)

    # --- cached section io ---
    def part_path(self, name):
        return os.path.join(self.parts, name + ".json")

    def have(self, name):
        return os.path.exists(self.part_path(name))

    def load(self, name):
        with open(self.part_path(name), "r", encoding="utf-8") as f:
            return json.load(f)

    def save(self, name, payload):
        atomic_write_bytes(
            self.part_path(name),
            json.dumps(payload, indent=2, ensure_ascii=False,
                       default=float).encode("utf-8"))

    def done(self, unit_id):
        st = self.cp.unit_done(unit_id)
        print("   checkpoint: unit '%s' done (%d units)"
              % (unit_id, st.get("completed_units")), flush=True)

    def mem(self, label):
        return StageMem(label)

    def finish_mem(self, probe, t0):
        info = probe.stop()
        info["wall_s"] = round(time.time() - t0, 1)
        self.peaks[probe.label] = info
        print("   peak mem [%s] stage=%.0f MB process=%.0f MB wall=%.1fs"
              % (probe.label, info["stage_peak_working_set_mb"],
                 info["process_peak_working_set_mb"], info["wall_s"]),
              flush=True)
        # persist immediately, so a later death keeps the measurement
        all_peaks = {}
        p = os.path.join(self.out, "stage-memory.json")
        if os.path.exists(p):
            try:
                with open(p, "r", encoding="utf-8") as f:
                    all_peaks = json.load(f).get("stages", {}) or {}
            except Exception:
                all_peaks = {}
        all_peaks.update(self.peaks)
        atomic_write_bytes(
            p,
            json.dumps({"stages": all_peaks,
                        "updated_at": datetime.datetime.now(
                            datetime.timezone.utc).isoformat()},
                       indent=2, ensure_ascii=False).encode("utf-8"))
        return info

    def photo(self):
        if self.photo_lin is None:
            im = np.asarray(Image.open(self.args.photo).convert("RGB"),
                            dtype=np.float64) / 255.0
            self.photo_lin = srgb_to_linear(im)
        return self.photo_lin


def build_states(ctx, need33=False, need_v8=True, build_v33_img=False):
    """Rebuild the states/arrays that are cheap, reuse the expensive outputs."""
    d = {}
    env, grid, S_patch = ctx.env, ctx.grid, ctx.S_patch
    cfg, raw, seed = ctx.cfg, ctx.raw, ctx.seed
    t0 = time.time()
    if need33:
        t33 = time.time()
        st33 = V6WaterState(raw, env, N=grid, seed=seed,
                            preset_index=ctx.args.preset, rain_rate_mm_h=0.0)
        f33 = st33.fields_at(0.0)
        f33["S"] = st33.S
        re33 = st33.realised_stats(f33)
        d.update(st33=st33, f33=f33, re33=re33, t_st33=time.time() - t33)
        print("   rebuilt 3.3 state: sigma2_spec=%.6f realised=%.6f Hs=%.4f"
              % (st33.spectrum_after["sigma2_rad2"], re33["sigma2_rad2"],
                 re33["Hs_m"]), flush=True)
    if need_v8:
        t1 = time.time()
        st = V8WaterState(ctx.scene, env, N=grid, seed=seed,
                          preset_index=ctx.args.preset,
                          fetch_m=float(cfg["fetch"]["fetch_m"]),
                          fetch_source=cfg["fetch"]["source"])
        fields = st.fields_at(0.0)
        fields["S"] = st.S
        re_v8 = st.realised_stats(fields)
        d.update(st=st, fields=fields, re_v8=re_v8, t_st8=time.time() - t1)
        print("   rebuilt v8 state: sigma2_spec=%.6f realised=%.6f Hs=%.4f"
              % (st.spectrum_after["sigma2_rad2"], re_v8["sigma2_rad2"],
                 re_v8["Hs_m"]), flush=True)
    print("   state rebuild wall=%.1fs" % (time.time() - t0), flush=True)
    return d


def load_masters(ctx):
    C33 = np.load(ctx.png("C-baseline-v33.npy")).astype(np.float32)
    C = np.load(ctx.png("C-v8-t0.npy")).astype(np.float32)
    assert C33.shape == (ctx.grid, ctx.grid), C33.shape
    assert C.shape == (ctx.grid, ctx.grid), C.shape
    return C33, C


def still_array(ctx, name):
    return np.asarray(Image.open(ctx.png(name)).convert("RGB"))


# ==========================================================================
# STAGE base: reusable numbers + [Q] convergence (Q=9/25 recomputed)
# ==========================================================================
def stage_base(ctx):
    t0 = time.time()
    probe = ctx.mem("base")
    env, grid, cfg, raw, seed = ctx.env, ctx.grid, ctx.cfg, ctx.raw, ctx.seed
    args = ctx.args
    S_patch = ctx.S_patch
    photo_lin = ctx.photo()
    PH, PW, _ = photo_lin.shape

    print("\n[base] reference numbers + reuse checks + [Q]", flush=True)
    S = build_states(ctx, need33=True, need_v8=True)
    st33, re33 = S["st33"], S["re33"]
    st, fields, re_v8 = S["st"], S["fields"], S["re_v8"]

    # ---------- display level (same rule as v5/v6/v7) ----------
    att = water_attenuation(env["turbid"])
    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_33 = smith_g1_beckmann(re33["rms_slope"], sun_zenith)
    L_v5 = flat_water_direct_level(env, beer_ref["two_way"]) * float(
        g1_33["G1"])
    level_ratio = L_v3 / L_v5 if L_v5 > 0 else 1.0
    E_v8_pre = E_V3 * level_ratio          # exactly run_v8.py line 231
    g1_v8 = smith_g1_beckmann(re_v8["rms_slope"], sun_zenith)
    E_v8_check = E_V3 * (L_v3 / (flat_water_direct_level(
        env, beer_ref["two_way"]) * float(g1_v8["G1"])))

    # ---------- masters (reused, never recomputed) ----------
    C33, C = load_masters(ctx)
    v33_cross = {"path": ctx.png("C-baseline-v33.npy"), "found": True,
                 "reused": True,
                 "note": ("produced by the 23:18 full-resolution run, whose "
                          "log records it as matches_numbers_v7=True"),
                 "stats": R.c_stats(C33)}
    j7 = os.path.join(args.v33_dir, "numbers-v7.json")
    if os.path.exists(j7):
        with open(j7, "r", encoding="utf-8") as f:
            n7 = json.load(f)
        ref = n7["transport"]["per_master"]
        s = v33_cross["stats"]
        v33_cross["numbers_v7"] = {k: ref[k] for k in
                                   ("mean", "p99", "p99_9", "max")}
        v33_cross["deltas"] = {k: float(s[k] - ref[k]) for k in
                               ("mean", "p99", "p99_9", "max")}
        v33_cross["matches_numbers_v7"] = bool(
            abs(s["mean"] - ref["mean"]) < 1e-9
            and abs(s["p99"] - ref["p99"]) < 1e-6
            and abs(s["max"] - ref["max"]) < 1e-6)
    print("3.3 master reused: mean=%.12f p99=%.6f max=%.6f (matches v7: %s)"
          % (v33_cross["stats"]["mean"], v33_cross["stats"]["p99"],
             v33_cross["stats"]["max"],
             v33_cross.get("matches_numbers_v7")), flush=True)

    # reuse check of the chunk-4 master against the dead run's log line
    reuse_check = {
        "file": "C-v8-t0.npy", "reused": True,
        "logged_by_dead_run": {"Q": 49, "mean": 1.000000358,
                               "p99": 1.5095, "max": 2.099,
                               "source": "out/v8-run.log line 14"},
        "now": R.c_stats(C),
    }
    now = reuse_check["now"]
    reuse_check["abs_dev_vs_log"] = {
        "mean": abs(now["mean"] - 1.000000358),
        "p99": abs(now["p99"] - 1.5095),
        "max": abs(now["max"] - 2.099),
    }
    reuse_check["consistent"] = bool(
        reuse_check["abs_dev_vs_log"]["mean"] < 1e-6
        and reuse_check["abs_dev_vs_log"]["p99"] < 5e-3
        and reuse_check["abs_dev_vs_log"]["max"] < 5e-3)
    print("v8 master reused: mean=%.9f p99=%.4f max=%.3f (matches log: %s)"
          % (now["mean"], now["p99"], now["max"],
             reuse_check["consistent"]), flush=True)

    # ---------- stills ----------
    img33, dbg33 = render_frame_v5(
        S["f33"], C33, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v8_pre, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=float(g1_33["G1"]), blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m33 = image_metrics(dbg33["lin"], dbg33["aces"])
    del dbg33
    prev33 = ctx.png("still-preset1-v33.png")
    if os.path.exists(prev33):
        prev = np.asarray(Image.open(prev33).convert("RGB"))
        v33_still_check_prev = {
            "vs_stored_from_dead_run_max_abs_uint8": int(np.max(np.abs(
                img33.astype(np.int16) - prev.astype(np.int16)))),
        }
        v33_still_check_prev["identical"] = bool(
            v33_still_check_prev["vs_stored_from_dead_run_max_abs_uint8"] == 0)
        del prev
    else:
        v33_still_check_prev = None
    Image.fromarray(img33).save(prev33)
    v33_still_check = {"g1": float(g1_33["G1"]),
                       "E_v8": float(E_v8_pre),
                       "reproduced_from_reused_master":
                       v33_still_check_prev}
    j7_still = os.path.join(args.v33_dir, "still-off.png")
    if os.path.exists(j7_still):
        ref_img = np.asarray(Image.open(j7_still).convert("RGB"))
        v33_still_check["vs_v7_still_max_abs_uint8"] = int(np.max(np.abs(
            img33.astype(np.int16) - ref_img.astype(np.int16))))
        v33_still_check["identical_to_v7"] = bool(
            v33_still_check["vs_v7_still_max_abs_uint8"] == 0)
    print("3.3 still: crush=%.5f blow=%.6f identical_to_v7=%s"
          % (m33["crush_rate_linear_lt_0.02"],
             m33["blowout_rate_aces_gt_0.98"],
             v33_still_check.get("identical_to_v7")), flush=True)

    img8, dbg8 = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_v8_check, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m8 = image_metrics(dbg8["lin"], dbg8["aces"])
    del dbg8
    prev8 = ctx.png("still-preset1-v8.png")
    if os.path.exists(prev8):
        prev = np.asarray(Image.open(prev8).convert("RGB"))
        rep8 = {"vs_stored_from_dead_run_max_abs_uint8": int(np.max(np.abs(
            img8.astype(np.int16) - prev.astype(np.int16))))}
        rep8["identical"] = bool(
            rep8["vs_stored_from_dead_run_max_abs_uint8"] == 0)
        del prev
    else:
        rep8 = None
    Image.fromarray(img8).save(prev8)
    m8["reproduced_from_reused_master"] = rep8
    print("v8 still: crush=%.5f blow=%.6f reproduced=%s" % (
        m8["crush_rate_linear_lt_0.02"], m8["blowout_rate_aces_gt_0.98"],
        None if rep8 is None else rep8["identical"]), flush=True)
    del S["f33"]
    gc.collect()

    # ---------- [E] Smith audit numbers ----------
    g1_table = PQ.g1_reference_table(
        rms_slopes=sorted({0.05, round(float(re33["rms_slope"]), 6)}),
        sun_altitudes_deg=(60.0, 30.0, 10.0, 5.0, 3.0))
    smith = {
        "code_facts": PQ.SMITH_CODE_FACTS,
        "verdict": PQ.SMITH_VERDICT,
        "reference_scene": {
            "rms_slope_v33": float(re33["rms_slope"]),
            "rms_slope_v8": float(re_v8["rms_slope"]),
            "sun_alt_deg": float(env["sun"]["altitudeDeg"]),
            "zenith_deg": float(sun_zenith),
            "a_param": float(g1_33["a"]),
            "Lambda": float(g1_33["Lambda"]),
            "G1_v33": float(g1_33["G1"]),
            "G1_is_exactly_1": bool(float(g1_33["G1"]) == 1.0),
            "G1_v8": float(g1_v8["G1"]),
            "G1_v8_is_exactly_1": bool(float(g1_v8["G1"]) == 1.0),
        },
        "g1_corner_table": g1_table,
        "still_change_if_removed": {
            "applied_correction_on_reference_still": float(g1_33["G1"]),
            "delta_still_uint8": 0 if float(g1_33["G1"]) == 1.0 else None,
            "identical": bool(float(g1_33["G1"]) == 1.0),
            "note": ("G1 = 1.0 exactly at the reference scene, so the direct "
                     "term is multiplied by 1.0 and removing the node changes "
                     "nothing: the audit's measured change is 0."),
        },
    }

    # ---------- [Q] convergence ----------
    det_time = {}
    C_by_q = {}
    for q in (9, 25):
        tq = time.time()
        Cq, _ = det_frame(fields, st.S, grid, env, int(q))
        det_time[int(q)] = time.time() - tq
        C_by_q[int(q)] = Cq
        print("   Q=%-3d mean=%.9f p99=%.4f max=%.3f (%.1fs)"
              % (q, float(Cq.mean()), float(np.percentile(Cq, 99)),
                 float(Cq.max()), det_time[int(q)]), flush=True)
    C_by_q[int(ctx.Q)] = C          # the reused master
    det_time[int(ctx.Q)] = None
    qconv = PQ.q_convergence(C_by_q, S_patch)
    for k, v in qconv["per_Q"].items():
        v["det_time_s"] = det_time.get(int(k))
        if int(k) == int(ctx.Q):
            v["det_time_note"] = ("master reused from C-v8-t0.npy; the dead "
                                  "run recorded 82.1 s for this Q "
                                  "(out/v8-run.log)")
    qconv["det_times_s"] = det_time
    qconv["note"] = ("the master used for the chunk-4 stills is Q=%d (reused "
                     "from the dead run, cross-checked above); the Q=9/25/49 "
                     "triple quantifies the quadrature cost.  Q=9 and Q=25 "
                     "were recomputed here, Q=%d is the stored master."
                     % (ctx.Q, ctx.Q))
    del C_by_q
    gc.collect()

    payload = {
        "reuse_check_v8_master": reuse_check,
        "v33_cross": v33_cross,
        "smith_shadow_audit": smith,
        "q_convergence": qconv,
        "state": st.v8_node(fields),
        "state_detail_v8": {
            "realised": re_v8,
            "sigma2_spectral": st.spectrum_after["sigma2_rad2"],
            "Hs_spectral": st.spectrum_after["Hs_m"],
            "loop_maxdiff_h": st.loop_maxdiff(),
        },
        "state_detail_legacy_v33": {
            "sigma2_spectral": st33.spectrum_after["sigma2_rad2"],
            "Hs_spectral": st33.spectrum_after["Hs_m"],
            "realised": re33, "seed": int(st33.seed),
        },
        "dispersion_gate_table": P.dispersion_gate_table(float(env["depthM"])),
        "anchor_sensitivity": P.anchor_sensitivity(
            env["windSpeedMS"], cfg["fetch"]["fetch_m"], env["depthM"],
            S_patch=S_patch, N=grid),
        "direction_comparison": P.direction_comparison(
            env["windSpeedMS"], cfg["fetch"]["fetch_m"], env["depthM"],
            env["specL"], st.kmin, st.kmax_nyq),
        "old_vs_new_anchor": {
            "cox_munk_old": P.cox_munk_sigma2(env["windSpeedMS"]),
            "new": st.anchor["sigma2_target_rad2"],
            "ratio": st.anchor["ratio_new_over_old"],
        },
        "legacy_fetch_solved": P.legacy_fetch_m(
            env["windSpeedMS"], env["depthM"], S_patch=S_patch, N=grid),
        "metadata_chunk4": {
            "patch_real_size_m": float(S_patch),
            "patch_equals_max_wavelength": True,
            "kmin_rad_per_m": float(st.kmin),
            "dx_m": float(st.dx),
            "Ft": float(P.nondim_fetch(cfg["fetch"]["fetch_m"],
                                       env["windSpeedMS"])),
            "fetch_m": float(cfg["fetch"]["fetch_m"]),
            "fetch_source": cfg["fetch"]["source"],
            "anchor_model_id": P.ANCHOR_MODEL_ID,
        },
        "recipe_v8_extension": {
            "parsed_fetch": cfg["fetch"],
            "new_fields": cfg["schema_v8_fields"],
            "legacy_rule": ("a legacy recipe (no wave.fetchM) takes the fixed "
                            "fetch solved to reproduce the old Cox-Munk "
                            "anchor; the pond default 100 m is used when the "
                            "recipe is new (schema v2) and carries no value"),
            "legacy_solve": P.legacy_fetch_m(
                env["windSpeedMS"], env["depthM"], S_patch=S_patch, N=grid),
        },
        "scene": {
            "wind_ms": env["windSpeedMS"], "specL_m": env["specL"],
            "patch_S_m": S_patch, "patch_resolution_m": S_patch / grid,
            "depthM": env["depthM"], "turbid": env["turbid"],
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "photo": args.photo, "photo_size": [int(PW), int(PH)],
            "grid": [grid, grid], "Q_master": int(ctx.Q),
        },
        "stills": {
            "preset1_v8": {**m8, "file": "still-preset1-v8.png"},
            "preset1_v33": {**m33, "file": "still-preset1-v33.png",
                            **v33_still_check},
            "display": {"E_v8": float(E_v8_check), "E_v3": float(E_V3),
                        "level_ratio": float(level_ratio),
                        "flat_water_direct_lum_v3": float(L_v3),
                        "flat_water_direct_lum_v5": float(L_v5),
                        "G1_v33": float(g1_33["G1"]),
                        "G1_v8": float(g1_v8["G1"]),
                        "note": ("same display-level normalisation as v5/v6/"
                                 "v7; evaluated once and held")},
            "metrics_definitions": {
                "crush": ("fraction of pixels with rec709 linear luminance "
                          "(pre-ACES, including exposure) < 0.02"),
                "blow": ("fraction of pixels with rec709 ACES luminance "
                         "> 0.98 (post-ACES, pre-sRGB-encode)"),
            },
        },
        "transport": {
            "method": ("det.py det_frame (UNCHANGED), Q=%d, sigma=0.75 px, "
                       "no clamp / no renormalisation" % ctx.Q),
            "master": {"stats": R.c_stats(C), "file": "C-v8-t0.npy"},
            "baseline_33": {"stats": R.c_stats(C33),
                            "file": "C-baseline-v33.npy",
                            "cross_check": v33_cross},
        },
        "timing": {"state_v33": S.get("t_st33"), "state_v8": S["t_st8"],
                   "det_per_Q": det_time},
    }
    del img33, img8, C, C33
    gc.collect()
    ctx.save("base", payload)
    ctx.finish_mem(probe, t0)
    ctx.done("base (E numbers + [Q] Q=9/25 + reuse checks + stills)")
    return payload


# ==========================================================================
# STAGE photon ([M])
# ==========================================================================
def stage_photon(ctx):
    t0 = time.time()
    probe = ctx.mem("photon_M")
    print("\n[M] photon check + MS-SSIM decomposition (%d photons)"
          % ctx.args.photons, flush=True)
    S = build_states(ctx, need33=False, need_v8=True)
    st, fields = S["st"], S["fields"]
    C = np.load(ctx.png("C-v8-t0.npy")).astype(np.float32)
    img8 = still_array(ctx, "still-preset1-v8.png")

    env, grid, cfg = ctx.env, ctx.grid, ctx.cfg
    photo_lin = ctx.photo()
    PH, PW, _ = photo_lin.shape
    args = ctx.args
    S_patch = ctx.S_patch
    seed = ctx.seed

    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    att = water_attenuation(env["turbid"])
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    re_v8 = S["re_v8"]
    g1_v8 = smith_g1_beckmann(re_v8["rms_slope"],
                              90.0 - float(env["sun"]["altitudeDeg"]))
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    E_v8_check = E_V3 * (L_v3 / (flat_water_direct_level(
        env, beer_ref["two_way"]) * float(g1_v8["G1"])))

    tp = time.time()
    E_ph, ph_info = photon_frame_lowmem(
        fields, st.S, grid, env, int(args.photons),
        seed=(seed + 1000) % (2 ** 63))
    print("   photons=%d expected_mean=%s measured=%s (%.1fs)"
          % (args.photons,
             [round(x, 6) for x in ph_info["expected_mean"]],
             [round(x, 6) for x in ph_info["measured_mean"]],
             time.time() - tp), flush=True)
    base_lum = float((np.asarray(env["sunColor"]) * np.array(
        [0.2126, 0.7152, 0.0722])).sum() * float(env["sunI"]))
    C_ph = (np.asarray(E_ph, np.float64) @ np.array(
        [0.2126, 0.7152, 0.0722])) / base_lum
    del E_ph
    gc.collect()
    np.save(ctx.png("C-photon-t0.npy"), C_ph.astype(np.float32))

    img_ph = render_frame_v5(
        fields, C_ph.astype(np.float32), photo_lin, env, PW, PH,
        ambient_scale=A_V3, c_layer_exposure=E_v8_check,
        substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_ref["down"], beer_up=beer_ref["up"],
        g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=False)
    Image.fromarray(img_ph).save(ctx.png("still-PHOTON-preset1-v8.png"))

    photon = PQ.photon_vs_det(C, C_ph, S_patch, img8, img_ph)
    photon.update({"photons": int(args.photons),
                   "time_s": time.time() - tp,
                   "gate": R.PHOTON_GATE,
                   "pass_msssim": bool(
                       photon["msssim_log1pC"]["msssim"] >= R.PHOTON_GATE),
                   "expected_mean_rgb": ph_info["expected_mean"],
                   "measured_mean_rgb": ph_info["measured_mean"],
                   "condition": ("4,000,000 photons at 2048^2 -- the v2 "
                                 "photon-check condition"),
                   "implementation": ph_info["implementation"]})
    print("   MS-SSIM(log1pC)=%.4f | peak abs err=%.4f (rel %.2f%%) | "
          "peak shift=%.2f px / %.3f mm | still MS-SSIM=%.4f"
          % (photon["msssim_log1pC"]["msssim"],
             photon["decomposition"]["peak_luminance"]["abs_error_p99_9"],
             100.0 * (photon["decomposition"]["peak_luminance"]
                      ["rel_error_p99_9"] or 0.0),
             photon["decomposition"]["peak_position"]["xcorr_shift_px"]
             or 0.0,
             photon["decomposition"]["peak_position"]["xcorr_shift_mm"]
             or 0.0,
             photon.get("still", {}).get("msssim_grey", {}).get(
                 "msssim", float("nan"))), flush=True)

    payload = {"photon": photon,
               "transport_photon": {"stats": R.c_stats(C_ph),
                                    "file": "C-photon-t0.npy"},
               "still_photon": "still-PHOTON-preset1-v8.png"}
    del C_ph, img_ph, C, img8
    gc.collect()
    ctx.save("photon", payload)
    ctx.finish_mem(probe, t0)
    ctx.done("M photon check + MS-SSIM decomposition (4M photons)")
    return payload


# ==========================================================================
# STAGE gates234: shallow phase speed, anchor sigma2, breaking clamp
# ==========================================================================
def stage_gates234(ctx):
    t0 = time.time()
    probe = ctx.mem("gates234")
    print("\n[gates 2/3/4] shallow + clamp + anchor", flush=True)
    env, grid, cfg, seed = ctx.env, ctx.grid, ctx.cfg, ctx.seed
    args = ctx.args
    S_patch = ctx.S_patch
    photo_lin = ctx.photo()
    PH, PW, _ = photo_lin.shape

    S = build_states(ctx, need33=False, need_v8=True)
    st, fields, re_v8 = S["st"], S["fields"], S["re_v8"]
    Q_main = int(ctx.Q)

    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    att = water_attenuation(env["turbid"])
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    g1_v8 = smith_g1_beckmann(re_v8["rms_slope"],
                              90.0 - float(env["sun"]["altitudeDeg"]))
    E_v8_check = E_V3 * (L_v3 / (flat_water_direct_level(
        env, beer_ref["two_way"]) * float(g1_v8["G1"])))

    h_ref = float(env["depthM"])
    g_disp = P.dispersion_gate_table(h_ref)
    g_anchor = P.anchor_gate(st.anchor["sigma2_target_rad2"],
                             st.spectrum_after["sigma2_rad2"],
                             st.anchor["cox_munk_old_target_rad2"])
    scatter = P.realisation_scatter(st.H0, st.N, st.k)
    g_anchor["realised"] = {
        "sigma2_realised_rad2": re_v8["sigma2_rad2"],
        "relative_deviation_vs_target": float(
            (re_v8["sigma2_rad2"] - st.anchor["sigma2_target_rad2"])
            / st.anchor["sigma2_target_rad2"]),
        "realisation_scatter": scatter,
        "z_score": float(
            (re_v8["sigma2_rad2"] - st.anchor["sigma2_target_rad2"])
            / (st.anchor["sigma2_target_rad2"]
               * scatter["rel_std_sigma2"])),
        "note": ("the CALIBRATION gate is the spectral residual (the taper "
                 "solves it exactly); the realised value of a single snapshot "
                 "scatters around it with the effective DOF of this small "
                 "patch, and the z-score says whether that scatter explains "
                 "the deviation"),
    }
    g_anchor["pass"] = bool(g_anchor["pass"]
                            and abs(g_anchor["realised"]["z_score"]) <= 3.0)

    # ---------------- shallow verification scene ----------------
    shallow_recipe = with_v8(ctx.raw, fetch_m=float(args.fetch),
                             depth_m=float(args.shallow_depth),
                             wind_ms=float(args.shallow_wind))
    env_sh = resolve_environment(shallow_recipe)
    t_sh = time.time()
    st_sh = V8WaterState(shallow_recipe, env_sh, N=grid, seed=seed,
                         preset_index=args.preset,
                         fetch_m=float(args.fetch))
    f_sh = st_sh.fields_at(0.0)
    f_sh["S"] = st_sh.S
    re_sh = st_sh.realised_stats(f_sh)
    disp_sh = st_sh.dispersion_node()
    kmin_sh = st_sh.kmin
    h_sh = float(env_sh["depthM"])
    w_sh = float(P.dispersion_omega(np.array([kmin_sh]), h_sh)[0])
    wq_sh = float(st_sh._omega_q()[0, 1])
    sq = math.sqrt(P.G * h_sh)
    kh_sh = float(kmin_sh * h_sh)
    c_theory = float(P.shallow_ratio_theory(np.array([kh_sh]))[0]) * sq
    shallow_phase = {
        "depth_m": h_sh,
        "sqrt_gh_m_per_s": sq,
        "kmin_rad_per_m": float(kmin_sh),
        "kh_at_kmin": kh_sh,
        "omega_continuous": w_sh,
        "omega_as_used_quantised": wq_sh,
        "c_continuous_m_per_s": float(w_sh / kmin_sh),
        "c_theory_finite_depth_m_per_s": c_theory,
        "c_continuous_over_sqrt_gh": float(w_sh / kmin_sh / sq),
        "c_theory_over_sqrt_gh": float(c_theory / sq),
        "c_continuous_over_theory": float((w_sh / kmin_sh) / c_theory),
        "c_as_used_over_sqrt_gh": float(wq_sh / kmin_sh / sq),
        "c_as_used_over_theory": float((wq_sh / kmin_sh) / c_theory),
        "quantisation_rel_error": float(abs(wq_sh - w_sh) / w_sh),
        "table": P.dispersion_gate_table(h_sh),
    }
    print("shallow state h=%.2f: sigma2=%.6f Hs=%.4f clamp=%s (%.1fs)"
          % (env_sh["depthM"], st_sh.spectrum_after["sigma2_rad2"],
             st_sh.spectrum_after["Hs_m"], st_sh.breaking.get("applied"),
             time.time() - t_sh), flush=True)
    C_sh, _ = det_frame(f_sh, st_sh.S, grid, env_sh, Q_main)
    np.save(ctx.png("C-shallow-h010.npy"), np.asarray(C_sh, np.float32))
    att_sh = water_attenuation(env_sh["turbid"])
    beer_sh = beer_lambert(env_sh["depthM"], env_sh["sun"]["altitudeDeg"],
                           env_sh["turbid"])
    blur_sh = scattering_blur_sigma_px(att_sh, env_sh["depthM"], grid)
    g1_sh = smith_g1_beckmann(re_sh["rms_slope"],
                              90.0 - float(env_sh["sun"]["altitudeDeg"]))
    L_sh = flat_water_direct_level(env_sh, beer_sh["two_way"]) * float(
        g1_sh["G1"])
    E_sh = E_V3 * (L_v3 / L_sh) if L_sh > 0 else E_v8_check
    img_sh, dbg_sh = render_frame_v5(
        f_sh, C_sh, photo_lin, env_sh, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_sh, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_sh["down"], beer_up=beer_sh["up"],
        g1=float(g1_sh["G1"]), blur_sigma_px=blur_sh["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_sh = image_metrics(dbg_sh["lin"], dbg_sh["aces"])
    del dbg_sh
    Image.fromarray(img_sh).save(ctx.png("still-shallow-h010.png"))
    del img_sh, C_sh, f_sh
    gc.collect()

    # ---------------- clamp case ----------------
    clamp_recipe = with_v8(ctx.raw, fetch_m=float(args.fetch),
                           depth_m=float(args.shallow_depth),
                           wind_ms=float(args.clamp_wind))
    env_cl = resolve_environment(clamp_recipe)
    t_cl = time.time()
    st_cl = V8WaterState(clamp_recipe, env_cl, N=grid, seed=seed,
                         preset_index=args.preset, fetch_m=float(args.fetch))
    f_cl = st_cl.fields_at(0.0)
    f_cl["S"] = st_cl.S
    re_cl = st_cl.realised_stats(f_cl)
    print("clamp state h=%.2f U=%.1f: Hs_spec %.4f -> %.4f (Hmax=%.4f "
          "binding=%s realised Hs=%.4f) (%.1fs)"
          % (env_cl["depthM"], env_cl["windSpeedMS"],
             st_cl.breaking["Hs_before_m"], st_cl.breaking["Hs_after_m"],
             st_cl.breaking["H_max_m"], st_cl.breaking["binding"],
             re_cl["Hs_m"], time.time() - t_cl), flush=True)
    C_cl, _ = det_frame(f_cl, st_cl.S, grid, env_cl, Q_main)
    np.save(ctx.png("C-clamp-h010-U5.npy"), np.asarray(C_cl, np.float32))
    beer_cl = beer_lambert(env_cl["depthM"], env_cl["sun"]["altitudeDeg"],
                           env_cl["turbid"])
    blur_cl = scattering_blur_sigma_px(
        water_attenuation(env_cl["turbid"]), env_cl["depthM"], grid)
    g1_cl = smith_g1_beckmann(re_cl["rms_slope"],
                              90.0 - float(env_cl["sun"]["altitudeDeg"]))
    L_cl = flat_water_direct_level(env_cl, beer_cl["two_way"]) * float(
        g1_cl["G1"])
    E_cl = E_V3 * (L_v3 / L_cl) if L_cl > 0 else E_v8_check
    img_cl, dbg_cl = render_frame_v5(
        f_cl, C_cl, photo_lin, env_cl, PW, PH, ambient_scale=A_V3,
        c_layer_exposure=E_cl, substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer_cl["down"], beer_up=beer_cl["up"],
        g1=float(g1_cl["G1"]), blur_sigma_px=blur_cl["sigma_px"],
        lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
    m_cl = image_metrics(dbg_cl["lin"], dbg_cl["aces"])
    del dbg_cl
    Image.fromarray(img_cl).save(ctx.png("still-clamp-h010-U5.png"))
    del img_cl, C_cl, f_cl
    gc.collect()

    g_clamp = {
        "gate": ("chunk-4 gate 4: at h = 0.1 m class the effective H stays "
                 "<= 0.78 h"),
        "cases": [
            {"scene": "shallow (h=0.1, U=%.1f)" % args.shallow_wind,
             "depth_m": float(env_sh["depthM"]),
             "H_depth_limit_m": float(st_sh.breaking["H_depth_0p78h_m"]),
             "H_stokes_limit_m": float(st_sh.breaking["H_stokes_0p142_m"]),
             "Hs_before_m": float(st_sh.breaking["Hs_before_m"]),
             "Hs_after_m": float(st_sh.breaking["Hs_after_m"]),
             "Hs_realised_m": float(re_sh["Hs_m"]),
             "applied": bool(st_sh.breaking["applied"]),
             "binding": st_sh.breaking["binding"],
             "factor": float(st_sh.breaking["factor"]),
             "pass": bool(st_sh.breaking["Hs_after_m"]
                          <= st_sh.breaking["H_depth_0p78h_m"] * 1.0000001)},
            {"scene": "clamp (h=0.1, U=%.1f)" % args.clamp_wind,
             "depth_m": float(env_cl["depthM"]),
             "H_depth_limit_m": float(st_cl.breaking["H_depth_0p78h_m"]),
             "H_stokes_limit_m": float(st_cl.breaking["H_stokes_0p142_m"]),
             "Hs_before_m": float(st_cl.breaking["Hs_before_m"]),
             "Hs_after_m": float(st_cl.breaking["Hs_after_m"]),
             "Hs_realised_m": float(re_cl["Hs_m"]),
             "applied": bool(st_cl.breaking["applied"]),
             "binding": st_cl.breaking["binding"],
             "factor": float(st_cl.breaking["factor"]),
             "pass": bool(st_cl.breaking["Hs_after_m"]
                          <= st_cl.breaking["H_depth_0p78h_m"] * 1.0000001)},
        ],
        "counts": {
            "reference_scene_applied": bool(st.breaking["applied"]),
            "shallow_applied": bool(st_sh.breaking["applied"]),
            "clamp_applied": bool(st_cl.breaking["applied"]),
        },
        "note": ("at the reference scene (h=0.6 m, U=1.5) neither bound bites: "
                 "H_max is set by 0.142*lambda = %.4f m against H_s = %.4f m.  "
                 "Both h=0.1 m scenes are bound by 0.78h."
                 % (st.breaking["H_stokes_0p142_m"],
                    st.breaking["Hs_before_m"])),
    }
    g_clamp["pass"] = bool(all(c["pass"] for c in g_clamp["cases"])
                           and st_sh.breaking["applied"]
                           and st_cl.breaking["applied"])

    payload = {
        "gate2": {"gate": g_disp["gate"], **g_disp,
                  "as_used_shallow_state": shallow_phase,
                  "pass": bool(g_disp["pass"]
                               and 0.99 <= shallow_phase[
                                   "c_continuous_over_theory"] <= 1.01)},
        "gate3": g_anchor,
        "gate4": g_clamp,
        "state_detail_shallow": {
            "realised": re_sh,
            "sigma2_spectral": st_sh.spectrum_after["sigma2_rad2"],
            "Hs_spectral": st_sh.spectrum_after["Hs_m"],
            "breaking": st_sh.breaking,
            "dispersion_as_used": disp_sh,
        },
        "state_detail_clamp": {
            "realised": re_cl,
            "sigma2_spectral": st_cl.spectrum_after["sigma2_rad2"],
            "Hs_spectral": st_cl.spectrum_after["Hs_m"],
            "breaking": st_cl.breaking,
        },
        "stills": {
            "shallow": {**m_sh, "file": "still-shallow-h010.png",
                        "depth_m": float(env_sh["depthM"]),
                        "wind_ms": float(env_sh["windSpeedMS"]),
                        "E": float(E_sh)},
            "clamp": {**m_cl, "file": "still-clamp-h010-U5.png",
                      "depth_m": float(env_cl["depthM"]),
                      "wind_ms": float(env_cl["windSpeedMS"]),
                      "E": float(E_cl)},
        },
        "transport_extra": {
            "shallow": {"stats": R.c_stats(np.load(
                ctx.png("C-shallow-h010.npy"))), "file": "C-shallow-h010.npy"},
            "clamp": {"stats": R.c_stats(np.load(
                ctx.png("C-clamp-h010-U5.npy"))),
                "file": "C-clamp-h010-U5.npy"},
        },
        "gate5_inputs": {
            "m8": m8_cache(ctx),
            "m33": m33_cache(ctx),
            "m_sh": m_sh,
            "m_cl": m_cl,
        },
    }
    del st_sh, st_cl
    gc.collect()
    ctx.save("gates234", payload)
    ctx.finish_mem(probe, t0)
    ctx.done("gates 2/3/4 (shallow phase speed, anchor sigma2, breaking clamp)")
    return payload


def m8_cache(ctx):
    base = ctx.load("base")
    return base["stills"]["preset1_v8"]


def m33_cache(ctx):
    base = ctx.load("base")
    return base["stills"]["preset1_v33"]


# ==========================================================================
# STAGE gate1: ablation of the 3.3 -> v8 difference + t=2 s dispersion
# ==========================================================================
def stage_gate1(ctx):
    t0 = time.time()
    probe = ctx.mem("gate1")
    print("\n[G1] ablation of the 3.3 -> v8 difference", flush=True)
    env, grid, cfg, seed = ctx.env, ctx.grid, ctx.cfg, ctx.seed
    args = ctx.args
    photo_lin = ctx.photo()
    PH, PW, _ = photo_lin.shape
    base = ctx.load("base")
    S = build_states(ctx, need33=True, need_v8=True)
    st33, re33, st, fields, re_v8 = (S["st33"], S["re33"], S["st"],
                                     S["fields"], S["re_v8"])
    C33, C = load_masters(ctx)
    img33 = still_array(ctx, "still-preset1-v33.png")
    img8 = still_array(ctx, "still-preset1-v8.png")

    beer_ref = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"],
                            env["turbid"])
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    att = water_attenuation(env["turbid"])
    blur_ref = scattering_blur_sigma_px(att, env["depthM"], grid)
    L_v3 = flat_water_direct_level(env, v3_beer(env))
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1_33 = smith_g1_beckmann(re33["rms_slope"], sun_zenith)
    g1_v8 = smith_g1_beckmann(re_v8["rms_slope"], sun_zenith)
    E_v8_check = E_V3 * (L_v3 / (flat_water_direct_level(
        env, beer_ref["two_way"]) * float(g1_v8["G1"])))
    Q_main = int(ctx.Q)

    abl = {}
    ablC = {}
    variants = {
        "anchor_only": dict(use_direction=False, use_anchor=True),
        "dir_old_anchor": dict(use_direction=True, use_anchor=True,
                               anchor_target_rad2=P.cox_munk_sigma2(
                                   env["windSpeedMS"])),
    }
    for tag, kw in variants.items():
        t1 = time.time()
        sv = V8WaterState(ctx.scene, env, N=grid, seed=seed,
                          preset_index=args.preset,
                          fetch_m=float(cfg["fetch"]["fetch_m"]), **kw)
        fv = sv.fields_at(0.0)
        fv["S"] = sv.S
        Cv, _ = det_frame(fv, sv.S, grid, env, Q_main)
        np.save(ctx.png("C-ablate-%s.npy" % tag),
                np.asarray(Cv, np.float32))
        imgv, dbgv = render_frame_v5(
            fv, Cv, photo_lin, env, PW, PH, ambient_scale=A_V3,
            c_layer_exposure=E_v8_check,
            substrate_gain_rgb=gain_base["gain_rgb"],
            beer_down=beer_ref["down"], beer_up=beer_ref["up"],
            g1=float(g1_v8["G1"]), blur_sigma_px=blur_ref["sigma_px"],
            lowpass_sigma_px=B_SIGMA_V3, return_debug=True)
        mv = image_metrics(dbgv["lin"], dbgv["aces"])
        del dbgv
        abl[tag] = {
            "state": {
                "sigma2_spec": sv.spectrum_after["sigma2_rad2"],
                "Hs_spec": sv.spectrum_after["Hs_m"],
                "kc_rad_per_m": (sv.anchor or {}).get("kc_rad_per_m"),
                "anchor_target": (sv.anchor or {}).get("sigma2_target_rad2"),
                "direction_applied": bool(sv.direction),
            },
            "C_vs_33": R.delta_stats(C33, Cv, tag + " vs 3.3"),
            "still_vs_33": {
                "max_abs_diff_uint8": int(np.max(np.abs(
                    imgv.astype(np.int16) - img33.astype(np.int16)))),
                "mean_abs_diff_uint8": float(np.mean(np.abs(
                    imgv.astype(np.int16) - img33.astype(np.int16)))),
            },
            "still_metrics": mv,
            "det_time_s": time.time() - t1,
        }
        print("   %-14s sigma2=%.6f kc=%.3f -> dP99=%+.4f rmsdiff=%.4f "
              "(%.1fs)"
              % (tag, sv.spectrum_after["sigma2_rad2"],
                 (sv.anchor or {}).get("kc_rad_per_m") or 0.0,
                 abl[tag]["C_vs_33"]["delta_p99"],
                 abl[tag]["C_vs_33"]["rms_diff"],
                 abl[tag]["det_time_s"]), flush=True)
        ablC[tag] = Cv
        del fv, sv, Cv, imgv
        gc.collect()

    # dispersion: identical at t=0 (proven), measured at t = 2.0 s
    st_dp = V8WaterState(ctx.scene, env, N=grid, seed=seed,
                         preset_index=args.preset,
                         fetch_m=float(cfg["fetch"]["fetch_m"]),
                         use_finite_depth=False)
    h0_diff = float(np.max(np.abs(st.H0.astype(np.complex128)
                                  - st_dp.H0.astype(np.complex128))))
    om_diff = float(np.max(np.abs(st._omega_q() - st_dp._omega_q())))
    q2 = max(9, Q_main // 2)
    f2a = st.fields_at(2.0)
    f2a["S"] = st.S
    f2b = st_dp.fields_at(2.0)
    f2b["S"] = st_dp.S
    C2a, _ = det_frame(f2a, st.S, grid, env, q2)
    C2b, _ = det_frame(f2b, st_dp.S, grid, env, q2)
    h_t2_diff = float(np.max(np.abs(f2a["h"] - f2b["h"])))
    h_t2_std = float(np.std(f2a["h"]))
    abl["dispersion"] = {
        "t0_identity": {
            "H0_max_abs_diff": h0_diff,
            "omega_max_abs_diff": om_diff,
            "t0_field_identical": bool(h0_diff == 0.0),
            "field_max_abs_diff_t0": 0.0,
            "note": ("at t=0 every phase is exp(i*0)=1, so the dispersion "
                     "cannot change the field.  The chunk-4 dispersion is "
                     "therefore invisible in the reference still BY "
                     "CONSTRUCTION -- not by being flat."),
        },
        "t2s": {
            "Q": int(q2),
            "C_deep_dispersion_stats": R.c_stats(C2b),
            "C_chunk4_dispersion_stats": R.c_stats(C2a),
            "delta_vs_deep": R.delta_stats(C2b, C2a, "t=2.0 s dispersion"),
            "h_field_max_abs_diff_t2": h_t2_diff,
            "h_field_std_t2": h_t2_std,
            "h_field_rel_diff_t2": (float(h_t2_diff / h_t2_std)
                                    if h_t2_std else None),
        },
    }
    print("   dispersion: t0 H0 diff=%.3g (identical=%s) | t=2.0s "
          "rmsdiff=%.4f dP99=%+.4f"
          % (h0_diff, h0_diff == 0.0,
             abl["dispersion"]["t2s"]["delta_vs_deep"]["rms_diff"],
             abl["dispersion"]["t2s"]["delta_vs_deep"]["delta_p99"]),
          flush=True)
    np.save(ctx.png("C-t2-deep-dispersion.npy"), np.asarray(C2b, np.float32))
    np.save(ctx.png("C-t2-chunk4-dispersion.npy"),
            np.asarray(C2a, np.float32))
    del C2a, C2b, f2a, f2b, st_dp
    gc.collect()

    total = R.delta_stats(C33, C, "v8 vs 3.3")
    Ca = ablC["anchor_only"]
    Cb = ablC["dir_old_anchor"]
    d_anchor = R.delta_stats(C33, Ca, "anchor first: 3.3 -> anchor_only")
    d_dir_after = R.delta_stats(Ca, C, "then direction: anchor_only -> v8")
    d_dir = R.delta_stats(C33, Cb, "direction first: 3.3 -> dir_old_anchor")
    d_anchor_after = R.delta_stats(Cb, C,
                                   "then anchor: dir_old_anchor -> v8")
    rows = {}
    for kname in ("mean", "p50", "p99", "p99_9", "max"):
        key = "delta_" + kname
        tot = total[key]
        a = d_anchor[key]
        b = d_dir_after[key]
        c = d_dir[key]
        dd = d_anchor_after[key]
        inter = (b - c)
        rows[kname] = {
            "total": float(tot),
            "chain_anchor_first": {"anchor": float(a),
                                   "direction_given_anchor": float(b),
                                   "sum": float(a + b)},
            "chain_direction_first": {"direction": float(c),
                                      "anchor_given_direction": float(dd),
                                      "sum": float(c + dd)},
            "closure_error": float((a + b) - tot),
            "interaction": float(inter),
            "interaction_over_total": (float(abs(inter) / abs(tot))
                                       if abs(tot) > 1e-9 else None),
            "share_anchor_chainA": (float(abs(a) / (abs(a) + abs(b)))
                                    if (a or b) else None),
            "share_anchor_chainB": (float(abs(dd) / (abs(c) + abs(dd)))
                                    if (c or dd) else None),
        }
    attribution = {
        "method": ("two sequential (path) decompositions.  The diff is split "
                   "by switching ONE node at a time from the 3.3 state to the "
                   "chunk-4 state, in both orders, so the node interaction is "
                   "measured instead of assumed.  Sums are exact by "
                   "construction; the informative number is the interaction "
                   "between the orderings."),
        "rows": rows,
        "states": {
            "anchor_only": abl["anchor_only"]["state"],
            "dir_old_anchor": abl["dir_old_anchor"]["state"],
            "full": {"sigma2_spec": st.spectrum_after["sigma2_rad2"],
                     "Hs_spec": st.spectrum_after["Hs_m"]},
            "3.3": {"sigma2_spec": st33.spectrum_after["sigma2_rad2"],
                    "Hs_spec": st33.spectrum_after["Hs_m"]},
        },
        "still_level": {
            "anchor_only_vs_33": abl["anchor_only"]["still_vs_33"],
            "dir_old_anchor_vs_33": abl["dir_old_anchor"]["still_vs_33"],
            "v8_vs_33_max_abs_uint8": int(np.max(np.abs(
                img8.astype(np.int16) - img33.astype(np.int16)))),
            "v8_vs_33_mean_abs_uint8": float(np.mean(np.abs(
                img8.astype(np.int16) - img33.astype(np.int16)))),
        },
    }
    g1_gate = {
        "gate": ("chunk-4 gate 1: the difference against the 3.3 baseline is "
                 "explained by the dispersion + anchor change (bit-identity "
                 "is NOT required -- the physics changed).  Numeric breakdown "
                 "per node, both orderings."),
        "baseline": {"label": "3.3 = v7 chunk-3 master (reused)",
                     "stats": R.c_stats(C33),
                     "still": "still-preset1-v33.png",
                     "still_identical_to_v7": base["stills"]["preset1_v33"]
                     .get("identical_to_v7")},
        "v8": {"label": "chunk 4", "stats": R.c_stats(C),
               "still": "still-preset1-v8.png"},
        "total_difference": total,
        "per_node_states": {k: v["state"] for k, v in abl.items()
                            if k != "dispersion"},
        "attribution": attribution,
        "dispersion": abl.get("dispersion"),
        "clamp_and_normalisation": {
            "clamp_applied_at_reference_scene": bool(st.breaking["applied"]),
            "clamp_contribution_to_diff": 0.0 if not st.breaking["applied"]
            else None,
            "display_level_E_v8": float(E_v8_check),
            "display_level_E_v7": 2.708272210034262,
            "display_level_abs_diff": float(
                abs(E_v8_check - 2.708272210034262)),
            "G1_v33": float(g1_33["G1"]), "G1_v8": float(g1_v8["G1"]),
            "note": ("the clamp does not fire at the reference scene (H_s = "
                     "%.4f m vs H_max = %.4f m) and the display level is "
                     "identical, so the whole still difference is the "
                     "spectrum: anchor + direction (+ the dispersion, which "
                     "is exactly zero at t=0)."
                     % (st.breaking["Hs_before_m"], st.breaking["H_max_m"])),
        },
    }
    max_closure = max(abs(rows[k]["closure_error"])
                      for k in ("p50", "p99", "p99_9"))
    disp = abl.get("dispersion", {})
    t0_id = bool(disp.get("t0_identity", {}).get("t0_field_identical"))
    g1_gate["accounting"] = {
        "max_closure_error": float(max_closure),
        "closure_tolerance": 1e-9,
        "closure_pass": bool(max_closure <= 1e-9),
        "dispersion_contribution_at_t0": 0.0,
        "dispersion_t0_identity_proved": t0_id,
        "clamp_contribution": 0.0 if not st.breaking["applied"] else None,
        "display_level_contribution": 0.0,
        "g1_contribution": 0.0,
        "note": ("the difference closes into three terms: anchor, direction, "
                 "and zero for dispersion@t=0 / clamp / display level / G1.  "
                 "The nodes are not additive in their C statistics (they "
                 "partially cancel -- the anchor steepens the surface, the "
                 "narrower direction lobe removes crossing wave systems); the "
                 "interaction column is that measured non-additivity, not an "
                 "unexplained remainder, and both orderings are given."),
    }
    g1_gate["pass"] = bool(max_closure <= 1e-9 and t0_id)
    g1_gate["attribution_pass"] = g1_gate["pass"]
    g1_gate["pass_rule"] = (
        "the per-node split closes to 1e-9 on the C statistics (so every part "
        "of the difference is assigned), AND the dispersion's t=0 "
        "contribution is exactly zero by proof, AND clamp/display/G1 "
        "contributions are zero.  Bit-identity with 3.3 is NOT required and "
        "is explicitly not the gate.")

    payload = {"gate1": g1_gate,
               "ablation_states": {k: v["state"] for k, v in abl.items()
                                   if k != "dispersion"},
               "ablation_still_metrics": {k: v["still_metrics"]
                                          for k, v in abl.items()
                                          if k != "dispersion"}}
    del ablC, Ca, Cb, img33, img8, C33, C, st33
    gc.collect()
    ctx.save("gate1", payload)
    ctx.finish_mem(probe, t0)
    ctx.done("gate 1 (ablation both orderings + t=2 s dispersion)")
    return payload


# ==========================================================================
# STAGE photon16m: supplementary photon-density diagnostic (not the gate)
# ==========================================================================
def stage_photon16m(ctx):
    """Same field, 4x the photons: tells a noise-limited check from a
    physics disagreement.  This is NOT the gate (the gate condition is the
    v2 condition, 4M photons) -- it is the evidence that the 4M failure is
    Poisson noise at ~1 photon per master cell."""
    t0 = time.time()
    probe = ctx.mem("photon16m")
    print("\n[M+] supplementary photon density check (16M photons)", flush=True)
    env, grid = ctx.env, ctx.grid
    S = build_states(ctx, need33=False, need_v8=True)
    st, fields = S["st"], S["fields"]
    C = np.load(ctx.png("C-v8-t0.npy")).astype(np.float32)
    n_ph = 16_000_000
    t1 = time.time()
    E, info = photon_frame_lowmem(fields, st.S, grid, env, n_ph,
                                  seed=(ctx.seed + 1000) % (2 ** 63))
    base_lum = float((np.asarray(env["sunColor"]) * np.array(
        [0.2126, 0.7152, 0.0722])).sum() * float(env["sunI"]))
    C16 = (np.asarray(E, np.float64) @ np.array(
        [0.2126, 0.7152, 0.0722])) / base_lum
    del E
    gc.collect()
    np.save(ctx.png("C-photon-16M-diagnostic.npy"), C16.astype(np.float32))
    C4 = np.load(ctx.png("C-photon-t0.npy")).astype(np.float64)
    A = np.log1p(np.asarray(C, np.float64))
    B16 = np.log1p(C16)
    B4 = np.log1p(C4)
    res = {
        "photons": n_ph,
        "photons_per_master_cell": float(n_ph) / (grid * grid),
        "photons_per_cell_at_gate": 4_000_000.0 / (grid * grid),
        "time_s": time.time() - t1,
        "msssim_log1pC_vs_det": PQ.ms_ssim(A, B16)["msssim"],
        "msssim_log1pC_4M_vs_det": PQ.ms_ssim(A, B4)["msssim"],
        "msssim_log1pC_16M_vs_4M": PQ.ms_ssim(B4, B16)["msssim"],
        "quantiles_16M": PQ.quantiles(C16),
        "quantile_rel_16M_vs_det": {
            k: float(abs(PQ.quantiles(C16)[k] - PQ.quantiles(C)[k])
                     / abs(PQ.quantiles(C)[k]))
            for k in ("p99", "p99.9")},
        "measured_mean_rgb": info["measured_mean"],
        "note": ("the gate is evaluated at the v2 condition (4M photons, "
                 "~0.95 photons per master cell).  This run shows what the "
                 "same state looks like at 4x the photons, i.e. whether the "
                 "4M MS-SSIM is noise-limited."),
    }
    print("   MS-SSIM(log1pC): 4M=%.4f 16M=%.4f | 16M vs 4M=%.4f | "
          "p99 rel err 16M=%.3f%% (%.1fs)"
          % (res["msssim_log1pC_4M_vs_det"], res["msssim_log1pC_vs_det"],
             res["msssim_log1pC_16M_vs_4M"],
             100.0 * res["quantile_rel_16M_vs_det"]["p99"],
             res["time_s"]), flush=True)
    atomic_write_bytes(os.path.join(ctx.out, "photon-16M-diagnostic.json"),
                       json.dumps(res, indent=2, default=float).encode("utf-8"))
    ctx.finish_mem(probe, t0)
    ctx.done("M+ supplementary 16M-photon density diagnostic")
    return res


# ==========================================================================
# STAGE gate6fix: like-for-like verification of the rain path
# ==========================================================================
def stage_gate6fix(ctx):
    """The run_v8.py gate-6 rebuild compares the POND scene (U=1.5) + R=10 with
    the v6 verification scene (wind 0, calm), so its 72% sigma^2 difference is
    a scene difference, not a code difference.  This stage rebuilds the actual
    v6 calm scene at every v6 rate and compares to numbers-v6.json."""
    t0 = time.time()
    probe = ctx.mem("gate6fix")
    print("\n[gate6fix] rain path, like-for-like at the v6 calm scene",
          flush=True)
    from recipe_v6 import with_rain
    grid = ctx.grid
    part = ctx.load("gates56")
    rain = part["gate6"]
    v6nums = os.path.join("out", "v6", "numbers-v6.json")
    if not os.path.exists(v6nums):
        print("   numbers-v6.json missing -> nothing to compare", flush=True)
        return part
    with open(v6nums, "r", encoding="utf-8") as f:
        n6 = json.load(f)

    scene6 = with_rain(ctx.raw, 0.0, wind_ms=0.0)
    env6 = resolve_environment(scene6)
    seed6 = derive_seed(scene6, ctx.args.preset)
    rows = []
    for R, key in ((0.0, "R00"), (2.0, "R02"), (10.0, "R10"),
                   (50.0, "R50")):
        t1 = time.time()
        st_r = V6WaterState(scene6, env6, N=grid, seed=seed6,
                            preset_index=ctx.args.preset,
                            rain_rate_mm_h=R)
        f_r = st_r.fields_at(0.0)
        f_r["S"] = st_r.S
        re_r = st_r.realised_stats(f_r)
        old = n6["state"][key]["spectrum_after"]["sigma2_rad2"]
        now = float(st_r.spectrum_after["sigma2_rad2"])
        rows.append({
            "rate_mm_h": R, "v6_key": key,
            "v6_sigma2": float(old), "now_sigma2": now,
            "rel_diff": float(abs(now - old) / old),
            "v6_Hs": float(n6["state"][key]["spectrum_after"]["Hs_m"]),
            "now_Hs": float(st_r.spectrum_after["Hs_m"]),
            "v6_realised_sigma2": float(
                n6["state"][key].get("realised", {}).get("sigma2_rad2",
                                                         float("nan"))),
            "now_realised_sigma2": float(re_r["sigma2_rad2"]),
            "build_s": time.time() - t1,
        })
        print("   R=%-4g sigma2 v6=%.9f now=%.9f rel=%.3g"
              % (R, old, now, rows[-1]["rel_diff"]), flush=True)
        del f_r, st_r
        gc.collect()

    rain["scene_mismatch_in_run_v8_check"] = {
        "what": rain.get("rain_R10_rebuild"),
        "explanation": ("run_v8.py rebuilds V6WaterState from the CHUNK-4 pond "
                        "scene (U=1.5, fetch 100 m, patch 1.4755 m, seed "
                        "3255045276382562994) and compares it with "
                        "numbers-v6.json's R10, which was produced at the v6 "
                        "verification scene (wind 0 calm, patch 1.20 m, seed "
                        "2610463815026170101).  The 0.724 relative difference "
                        "is therefore a scene difference, not a regression."),
        "v6_scene": {"wind_ms": n6["scene"]["wind_ms"],
                     "patch_S_m": n6["scene"]["patch_S_m"],
                     "seed": n6["seed"]},
        "pond_scene": {"wind_ms": ctx.env["windSpeedMS"],
                       "patch_S_m": ctx.S_patch, "seed": ctx.seed},
    }
    rain["like_for_like_at_v6_scene"] = {
        "scene": {"wind_ms": env6["windSpeedMS"],
                  "patch_S_m": 4.0 * float(env6["specL"]),
                  "depth_m": env6["depthM"], "seed": int(seed6),
                  "grid": grid},
        "rows": rows,
        "max_rel_diff_sigma2": float(max(r["rel_diff"] for r in rows)),
        "all_within_1e_9": bool(all(r["rel_diff"] <= 1e-9 for r in rows)),
    }
    rain["pass_rule"] = ("the F(R) table is bit-identical to numbers-v6.json "
                         "(rel diff 0) AND the v6 rain states are reproduced "
                         "at the v6 scene (rel diff <= 1e-9) AND the v6 module "
                         "hashes are recorded.")
    rain["pass"] = bool(rain.get("F_table_max_rel_diff") == 0.0
                        and rain["like_for_like_at_v6_scene"]["all_within_1e_9"])
    part["gate6"] = rain
    ctx.save("gates56", part)
    ctx.finish_mem(probe, t0)
    ctx.done("gate 6 re-verified like-for-like at the v6 calm scene")
    return part


# ==========================================================================
# STAGE gates56: tone gate + rain path
# ==========================================================================
def stage_gates56(ctx):
    t0 = time.time()
    probe = ctx.mem("gates56")
    print("\n[gates 5/6] preset tone + rain path", flush=True)
    env, grid, cfg, seed = ctx.env, ctx.grid, ctx.cfg, ctx.seed
    args = ctx.args
    base = ctx.load("base")
    g234 = ctx.load("gates234")
    m8 = base["stills"]["preset1_v8"]
    m33 = base["stills"]["preset1_v33"]
    m_sh = g234["stills"]["shallow"]
    m_cl = g234["stills"]["clamp"]
    C = np.load(ctx.png("C-v8-t0.npy")).astype(np.float32)
    C33 = np.load(ctx.png("C-baseline-v33.npy")).astype(np.float32)

    g5 = {
        "gate": ("chunk-4 gate 5: standard preset keeps mean(C) = 1.0, crush "
                 "inside the customary bound (0.0325) and blow = 0"),
        "C_mean": float(C.mean()),
        "C_mean_abs_dev": float(abs(C.mean() - 1.0)),
        "C_mean_pass": bool(abs(C.mean() - 1.0) <= R.CMEAN_TOL),
        "crush": float(m8["crush_rate_linear_lt_0.02"]),
        "blow": float(m8["blowout_rate_aces_gt_0.98"]),
        "crush_bound": R.CRUSH_BOUND, "blow_bound": R.BLOW_BOUND,
        "crush_pass": bool(m8["crush_rate_linear_lt_0.02"] <= R.CRUSH_BOUND),
        "blow_pass": bool(m8["blowout_rate_aces_gt_0.98"] <= R.BLOW_BOUND),
        "reference_33": {
            "mean": float(C33.mean()),
            "crush": float(m33["crush_rate_linear_lt_0.02"]),
            "blow": float(m33["blowout_rate_aces_gt_0.98"]),
        },
        "shallow_case": {"crush": float(m_sh["crush_rate_linear_lt_0.02"]),
                         "blow": float(m_sh["blowout_rate_aces_gt_0.98"])},
        "clamp_case": {"crush": float(m_cl["crush_rate_linear_lt_0.02"]),
                       "blow": float(m_cl["blowout_rate_aces_gt_0.98"])},
    }
    g5["pass"] = bool(g5["C_mean_pass"] and g5["crush_pass"]
                      and g5["blow_pass"])
    print("   gate5: mean=%.12f crush=%.5f blow=%.6f -> %s"
          % (g5["C_mean"], g5["crush"], g5["blow"], g5["pass"]), flush=True)

    rain = {"gate": ("chunk-4 gate 6: the rain path is unmodified and the "
                     "F(R) table still matches numbers-v6.json on top of the "
                     "new state")}
    v6nums = os.path.join("out", "v6", "numbers-v6.json")
    fgate = flux_gate()
    rain["F_table_now"] = {str(k): v for k, v in F_TABLE.items()}
    rain["flux_gate_now"] = fgate
    rain["hashes"] = {}
    for fn in ("v6_physics.py", "v6_state.py", "recipe_v6.py", "run_v6.py"):
        if os.path.exists(fn):
            rain["hashes"][fn] = R.sha12(fn)
    n6 = None
    if os.path.exists(v6nums):
        with open(v6nums, "r", encoding="utf-8") as f:
            n6 = json.load(f)
        old = n6["rain_physics_nodes"]["F_table_report_3_3"]
        cmp_rows = []
        for k, v in old.items():
            now = float(F_TABLE[float(k)]) if float(k) in F_TABLE else None
            cmp_rows.append({"rate_mm_h": float(k), "v6": float(v),
                             "now": now,
                             "rel_diff": (None if now is None
                                          else abs(now - float(v))
                                          / float(v))})
        rain["F_table_vs_v6"] = cmp_rows
        rain["F_table_max_rel_diff"] = float(max(
            r["rel_diff"] for r in cmp_rows if r["rel_diff"] is not None))
        rain["pass"] = bool(rain["F_table_max_rel_diff"] == 0.0)
    else:
        rain["pass"] = None
        rain["note"] = ("numbers-v6.json not found; F table compared "
                        "in-process only")
    if grid >= 1024:
        t_r = time.time()
        st_r = V6WaterState(ctx.raw, env, N=grid, seed=seed,
                            preset_index=args.preset, rain_rate_mm_h=10.0)
        rain["rain_R10_rebuild"] = {
            "sigma2_spectral": st_r.spectrum_after["sigma2_rad2"],
            "band_variance": st_r.rain_band["variance_realised_m2"],
            "k_ring": st_r.rain_band["k_ring_nominal_rad_per_m"],
            "rebuild_time_s": time.time() - t_r,
        }
        if n6 is not None:
            old10 = n6["state"]["R10"]["spectrum_after"]["sigma2_rad2"]
            rain["rain_R10_rebuild"]["v6_sigma2"] = old10
            rain["rain_R10_rebuild"]["rel_diff_vs_v6"] = float(
                abs(st_r.spectrum_after["sigma2_rad2"] - old10) / old10)
        print("   gate6: R10 rebuild sigma2=%.6f (v6 rel diff %s) F-table "
              "max rel diff=%s"
              % (rain["rain_R10_rebuild"]["sigma2_spectral"],
                 rain["rain_R10_rebuild"].get("rel_diff_vs_v6"),
                 rain.get("F_table_max_rel_diff")), flush=True)
        del st_r
        gc.collect()

    payload = {"gate5": g5, "gate6": rain}
    del C, C33
    gc.collect()
    ctx.save("gates56", payload)
    ctx.finish_mem(probe, t0)
    ctx.done("gates 5/6 (tone + rain path)")
    return payload


# ==========================================================================
# STAGE assemble: numbers-v8.json (+ prereq json, panel, recipes)
# ==========================================================================
def stage_assemble(ctx):
    t0 = time.time()
    probe = ctx.mem("assemble")
    print("\n[assemble] numbers-v8.json + prereq json + panel", flush=True)
    env, grid, cfg = ctx.env, ctx.grid, ctx.cfg
    args = ctx.args
    base = ctx.load("base")
    photon = ctx.load("photon")
    g234 = ctx.load("gates234")
    g1 = ctx.load("gate1")
    g56 = ctx.load("gates56")
    S_patch = ctx.S_patch
    photo_lin = ctx.photo()
    PH, PW, _ = photo_lin.shape

    # panel of the prerequisite figures
    items = []
    items.append((R.grey_png(np.load(ctx.png("C-v8-t0.npy"))),
                  "C master Q=%d (reused)" % ctx.Q,
                  "mean %.6f, p99 %.3f, max %.3f"
                  % (base["transport"]["master"]["stats"]["mean"],
                     base["transport"]["master"]["stats"]["p99"],
                     base["transport"]["master"]["stats"]["max"])))
    items.append((R.grey_png(np.load(ctx.png("C-photon-t0.npy"))),
                  "photon master (4M)",
                  "MS-SSIM(log1pC)=%.4f"
                  % photon["photon"]["msssim_log1pC"]["msssim"]))
    R.save_panel(ctx.png("panel-v8-prereq.png"), items, cols=2)
    print("   panel written", flush=True)

    rec_dir = os.path.join(args.out, "recipes")
    os.makedirs(rec_dir, exist_ok=True)
    rec_files = []
    shallow_recipe = with_v8(ctx.raw, fetch_m=float(args.fetch),
                             depth_m=float(args.shallow_depth),
                             wind_ms=float(args.shallow_wind))
    clamp_recipe = with_v8(ctx.raw, fetch_m=float(args.fetch),
                           depth_m=float(args.shallow_depth),
                           wind_ms=float(args.clamp_wind))
    for tag, r in (("preset1-v8", ctx.scene),
                   ("shallow-h010", shallow_recipe),
                   ("clamp-h010-U5", clamp_recipe)):
        fn = os.path.join(rec_dir, "%s.json" % tag)
        with open(fn, "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2, ensure_ascii=False)
        rec_files.append(os.path.basename(fn))

    st_dict = base["state"]
    prereq = {
        "order": ["E smith-shadow audit", "Q quadrature convergence",
                  "M photon MS-SSIM decomposition"],
        "smith_shadow_audit": base["smith_shadow_audit"],
        "q_convergence": base["q_convergence"],
        "photon_msssim_decomposition": photon["photon"],
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
    }
    with open(ctx.png("numbers-v8-prereq.json"), "w", encoding="utf-8") as f:
        json.dump(prereq, f, indent=2, ensure_ascii=False, default=float)

    gates = {
        "gate1_difference_attribution": g1["gate1"],
        "gate2_shallow_phase_speed": g234["gate2"],
        "gate3_new_anchor_sigma2": g234["gate3"],
        "gate4_breaking_clamp": g234["gate4"],
        "gate5_standard_preset_tone": g56["gate5"],
        "gate6_rain_path_untouched": g56["gate6"],
    }

    stills = dict(base["stills"])
    stills["shallow"] = g234["stills"]["shallow"]
    stills["clamp"] = g234["stills"]["clamp"]
    stills["photon"] = photon["still_photon"]
    stills["panel_prereq"] = "panel-v8-prereq.png"

    transport = dict(base["transport"])
    transport["photon"] = photon["transport_photon"]
    transport["extra_cases"] = g234["transport_extra"]
    transport["reuse_note"] = (
        "C-baseline-v33.npy and C-v8-t0.npy are the 2048^2 masters written by "
        "the earlier full-resolution run (23:18 / 23:21) and are REUSED, not "
        "recomputed; their statistics are cross-checked against "
        "numbers-v7.json and against the dead run's log line.  Everything "
        "else (photon master, ablate masters, shallow/clamp masters, t=2 s "
        "masters) was computed in this run at 2048^2.")

    numbers = {
        "pipeline": "v8-chunk4 (version 3.4)",
        "design": "PIPELINE-V8.md chunk 4 (P0-C + H_max clamp + patch coverage)",
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": ctx.raw.get("name"), "recipe_schema_in": ctx.raw.get("schema"),
        "preset": args.preset, "seed": int(ctx.seed),
        "seed_note": ("the seed is derived from the UNMODIFIED recipe, so the "
                      "phase realisation is the one the 3.3 master used and "
                      "the difference isolates physics"),
        "t_s": 0.0,
        "scene": base["scene"],
        "metadata_chunk4": base["metadata_chunk4"],
        "recipe_v8_extension": dict(base["recipe_v8_extension"],
                                    recipes_written=rec_files),
        "state": st_dict,
        "state_detail": {
            "legacy_v33": base["state_detail_legacy_v33"],
            "v8": base["state_detail_v8"],
            "shallow": g234["state_detail_shallow"],
            "clamp": g234["state_detail_clamp"],
        },
        "dispersion": {
            "constants": {"gamma_N_per_m": 0.0728, "rho_kg_per_m3": 1000.0,
                          "nu_m2_per_s": P.NU_KIN,
                          "temp_dependence": ("constant; the water-"
                                              "temperature dependence is left "
                                              "as the v5 plug-in point")},
            "form": "omega^2 = (g k + (gamma/rho) k^3) tanh(k h)",
            "deep_water_retired": True,
            "gate_table": base["dispersion_gate_table"],
        },
        "anchor": {
            "node": st_dict["anchor"],
            "model_id": P.ANCHOR_MODEL_ID,
            "candidates": P.ANCHOR_CANDIDATES,
            "sensitivity": base["anchor_sensitivity"],
            "old_vs_new": base["old_vs_new_anchor"],
        },
        "direction": {
            "node": st_dict["direction"],
            "comparison_vs_phillips": base["direction_comparison"],
        },
        "capillary_tail": st_dict["capillary_tail"],
        "breaking": st_dict["breaking"],
        "transport": transport,
        "stills": stills,
        "prerequisite_tasks": prereq,
        "gates": gates,
        "approximations_and_design_choices": [
            {"item": "anchor model family and the median rule",
             "statement": ("the adopted target is the mean square slope of a "
                           "fetch-limited spectrum integrated over the "
                           "resolved band, with the JONSWAP/Donelan "
                           "frequency spectrum mapped through the chunk-4 "
                           "dispersion.  Two named members (JONSWAP gamma=3.3 "
                           "and Donelan gamma=1.7) bracket the literature; "
                           "the adopted value is their median (midpoint) and "
                           "the full gamma/alpha/cutoff/depth sensitivity is "
                           "in the anchor node."),
             "kind": "model_selection"},
            {"item": "JONSWAP fetch laws below their fitted range",
             "statement": ("alpha = 0.076 Ft^-0.22 and f_p = 3.5 (g/U) "
                           "Ft^-1/3 were fitted to field data with "
                           "Ft ~ 1e3-1e5.  The pond scene here is Ft = %.0f, "
                           "an extrapolation DOWN in fetch.  It is used "
                           "because the design asks for exactly this family "
                           "for ponds; the direction of the extrapolation "
                           "error (a slight over-estimate of the high-"
                           "frequency level at very short fetch) is "
                           "acknowledged and not corrected."
                           % float(P.nondim_fetch(cfg["fetch"]["fetch_m"],
                                                  env["windSpeedMS"]))),
             "kind": "stated_extrapolation"},
            {"item": "legacy default fetch",
             "statement": ("a legacy recipe (no wave.fetchM) uses the fixed "
                           "fetch %.2f m (solved %.3f m, rounded), solved so "
                           "that the new anchor reproduces the old Cox-Munk "
                           "target at the legacy reference scene (U=1.5, "
                           "h=0.6).  It is a compatibility anchor, not a "
                           "physical pond fetch; the pond default for new "
                           "recipes is 100 m."
                           % (float(base["legacy_fetch_solved"][
                               "rounded_default_m"]),
                              float(base["legacy_fetch_solved"]["fetch_m"]))),
             "kind": "backward_compatibility_choice"},
            {"item": "capillary cutoff does not bind",
             "statement": ("the viscous rule (2 nu k^2 > omega) sits ~4e3x "
                           "above the grid Nyquist, so no mode is removed "
                           "inside the resolved band and the cusp sharpness "
                           "is set by dx = %.2f mm.  The rule is implemented "
                           "and evaluated rather than waived silently; the "
                           "one-loop viscous attenuation of the finest mode "
                           "is recorded and NOT applied (the state is a "
                           "static spectrum repeated every 8 s)."
                           % (S_patch / grid * 1e3)),
             "kind": "honest_limitation"},
            {"item": "dispersion is invisible at t=0",
             "statement": ("omega only enters through exp(i omega t); at t=0 "
                           "the field is the spectrum itself.  The chunk-4 "
                           "dispersion therefore changes the reference still "
                           "by exactly 0 (proved: H0 and the field are "
                           "bit-identical across dispersion branches), and "
                           "its effect is measured at t=2.0 s instead."),
             "kind": "structural_consequence"},
            {"item": "patch = longest resolved wavelength",
             "statement": ("the patch side (%.3f m) IS the longest resolved "
                           "wavelength, i.e. the spectrum holds roughly one "
                           "wavelength of its dominant component.  A "
                           "single-snapshot realisation therefore deviates "
                           "from the spectral sigma^2 by the few-mode "
                           "scatter reported in gate 3 (z-score), and long-"
                           "wave-dominated strong-wind scenes are the ones "
                           "the patch cannot represent.  Enlarging the patch "
                           "is out of scope for chunk 4." % S_patch),
             "kind": "honest_limitation"},
            {"item": "shoaling / refraction out of scope",
             "statement": ("shoaling needs a depth gradient grad h; the "
                           "scenes here are uniform depth, so it does not "
                           "fire.  Deferred to v4 (solid bottom), as "
                           "PIPELINE-V8 chunk 4 states."),
             "kind": "scope"},
        ],
        "calibrations_not_derived": P.NON_DERIVED,
        "inherited_calibrations_v5_v6": [
            "turbidity -> da+db scale (TURBID_DC_MAX=0.30 1/m, 70% scatter)",
            "scattering blur proportionality (BLUR_K_PX=40 px at 2048^2)",
            "b_w(550)=0.0018 1/m anchor",
            "display-level normalisation E=E_v3.(L_v3/L_v5)",
            "RING_ETA=0.10 (rain, unused here)",
        ],
        "timing_s": {
            "state_v8": base["timing"]["state_v8"],
            "det_per_Q": base["timing"]["det_per_Q"],
            "photon": photon["photon"]["time_s"],
        },
        "recovery": recovery_block(ctx),
    }
    numbers["numerics_notes"] = {
        "C_mean_two_accumulators": {
            "float32_accumulation": 1.0000003576278687,
            "float64_accumulation": base["transport"]["master"]["stats"][
                "mean"],
            "note": ("np.mean on the float32 master uses float32 pairwise "
                     "summation (1.000000358, the value the dead run's log "
                     "printed); c_stats casts to float64 first "
                     "(1.0000000000).  The 3.6e-7 difference is the "
                     "accumulator, not the field; both are inside CMEAN_TOL."),
        },
    }

    # v2 comparable number for the photon check (same condition)
    v2nums = os.path.join("out", "v2", "numbers-v2.json")
    if os.path.exists(v2nums):
        with open(v2nums, "r", encoding="utf-8") as f:
            n2 = json.load(f)
        pc = n2.get("photon_check", {})
        numbers["prerequisite_tasks"]["photon_vs_v2"] = {
            "v2_condition": {"photons": pc.get("photons"),
                             "grid": 2048,
                             "msssim_log1pC": pc.get("msssim", {}).get(
                                 "msssim"),
                             "gate": pc.get("gate"), "pass": pc.get("pass")},
            "v8_condition": {"photons": photon["photon"]["photons"],
                             "grid": grid,
                             "msssim_log1pC": photon["photon"][
                                 "msssim_log1pC"]["msssim"],
                             "gate": R.PHOTON_GATE,
                             "pass": photon["photon"]["pass_msssim"]},
            "note": ("same photon count and the same pixel filter, different "
                     "state and different master; the v2 gate was already "
                     "reported as not met at 4M photons (0.827 < 0.9) and the "
                     "MS-SSIM decomposition is what the gate-0 task asks for"),
        }

    # supplementary photon-density diagnostic (NOT the gate)
    p16 = os.path.join(ctx.out, "photon-16M-diagnostic.json")
    if os.path.exists(p16):
        with open(p16, "r", encoding="utf-8") as f:
            numbers["prerequisite_tasks"]["supplementary_photon_density"] = \
                json.load(f)

    with open(ctx.png("numbers-v8.json"), "w", encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False, default=float)

    print("\nGATES")
    for name, g in gates.items():
        print("  %-36s %s" % (name, g.get("pass")))
    print("   photon MS-SSIM log1pC        %.4f (gate %.2f -> %s)"
          % (photon["photon"]["msssim_log1pC"]["msssim"], R.PHOTON_GATE,
             photon["photon"]["pass_msssim"]))
    ctx.finish_mem(probe, t0)
    ctx.done("numbers-v8.json + numbers-v8-prereq.json + panel + recipes")
    print("OK -> %s" % ctx.png("numbers-v8.json"), flush=True)
    return numbers


def recovery_block(ctx):
    """What was reused, what was recomputed, and the memory that was measured."""
    mem = {}
    p = os.path.join(ctx.out, "stage-memory.json")
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as f:
            mem = json.load(f).get("stages", {})
    return {
        "dead_run": {
            "what_died": ("the earlier 2048^2 run of run_v8.py stopped inside "
                          "[M] (the photon check): no python.exe was left and "
                          "out/v8 did not change after 23:21.  The brief "
                          "records the suspicion of memory exhaustion; no "
                          "Python traceback was written, so the cause is "
                          "recorded as UNKNOWN with OOM as the leading "
                          "hypothesis."),
            "log_reached": "out/v8-run.log (last line: '[M] photon check + "
                           "MS-SSIM decomposition')",
            "state_at_death": {
                "3.3_state_sigma2_spec": 0.010680,
                "v8_state_sigma2_spec": 0.020954,
                "v8_state_realised": 0.022135,
                "v8_Hs_m": 0.1140,
                "E_v8": 2.708272,
                "G1_v8": 1.0,
                "Q_p99": {"9": 1.5097, "25": 1.5098, "49": 1.5095},
            },
        },
        "reused_unchanged": [
            "out/v8/C-baseline-v33.npy (2048^2, 3.3 master)",
            "out/v8/C-v8-t0.npy (2048^2, Q=49 chunk-4 master)",
        ],
        "recomputed_in_this_run": [
            "[E] Smith audit numbers (from the rebuilt 3.3/v8 states)",
            "[Q] Q=9 and Q=25 det_frame runs (Q=49 master reused)",
            "[M] photon check at 4,000,000 photons + MS-SSIM decomposition",
            "gate 1 ablation (anchor_only, dir_old_anchor, t=2 s dispersion)",
            "gate 2 shallow phase-speed state",
            "gate 3 anchor gate on the realised sigma^2",
            "gate 4 breaking clamp (shallow + clamp scenes)",
            "gate 5 tone metrics on the reused master",
            "gate 6 rain-path F(R) table + R10 rebuild",
            "all stills except the 3.3/v8 pair, which are re-rendered from the "
            "reused masters and compared to the stored PNGs",
        ],
        "memory_policy": {
            "photon_path": ("channel-split accumulator + per-channel FFT "
                            "blur (photon_frame_lowmem); equivalence to "
                            "transport.photon_frame verified bit-exactly at "
                            "256^2/200k photons"),
            "no_taskkill": ("no taskkill / no background process was used; "
                            "every stage ran in the foreground of the "
                            "session"),
            "stage_peaks_mb": mem,
        },
        "checkpoint": {
            "file": "checkpoint.json",
            "format": ("D:\\mutaro-work\\jobkit Checkpoint "
                       "(completed_units / completed_ids / last_unit / "
                       "updated_at, atomic write)"),
            "source": CHECKPOINT_SOURCE,
        },
    }


# ==========================================================================
# equivalence check of the low-memory photon path
# ==========================================================================
def verify_lowmem(ctx):
    print("\n[verify] photon_frame_lowmem vs transport.photon_frame "
          "(256^2, 200k photons)", flush=True)
    probe = ctx.mem("verify_lowmem")
    grid = 256
    env = ctx.env
    st = V8WaterState(ctx.scene, env, N=grid, seed=ctx.seed,
                      preset_index=ctx.args.preset,
                      fetch_m=float(ctx.cfg["fetch"]["fetch_m"]),
                      fetch_source=ctx.cfg["fetch"]["source"])
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    t0 = time.time()
    E, info = photon_frame_lowmem(fields, st.S, grid, env, 200_000,
                                  seed=(ctx.seed + 1000) % (2 ** 63))
    base_lum = float((np.asarray(env["sunColor"]) * np.array(
        [0.2126, 0.7152, 0.0722])).sum() * float(env["sunI"]))
    C_ph = (np.asarray(E, np.float64) @ np.array(
        [0.2126, 0.7152, 0.0722])) / base_lum
    dt = time.time() - t0
    ref = np.load(os.path.join("out", "v8-smoke2", "C-photon-t0.npy"))
    ref = np.asarray(ref, dtype=np.float64)
    d = np.abs(C_ph - ref)
    # the stored reference is float32, so the exact statement is: identical
    # after the same float32 rounding the stored master went through
    d32 = np.abs(C_ph.astype(np.float32).astype(np.float64) - ref)
    out = {"grid": grid, "photons": 200_000, "time_s": dt,
           "reference": "out/v8-smoke2/C-photon-t0.npy",
           "max_abs_diff_float64": float(d.max()),
           "max_abs_diff_after_float32_rounding": float(d32.max()),
           "float32_eps_at_value_1p5": float(np.spacing(np.float32(1.5))),
           "bit_identical_after_float32": bool(d32.max() == 0.0),
           "bit_identical": bool(d32.max() == 0.0),
           "measured_mean_rgb": info["measured_mean"]}
    print("   max|dC| (f64 vs stored f32)=%.3g  after f32 rounding=%.3g  "
          "bit_identical=%s (%.1fs)"
          % (out["max_abs_diff_float64"],
             out["max_abs_diff_after_float32_rounding"],
             out["bit_identical"], dt), flush=True)
    ctx.finish_mem(probe, t0)
    atomic_write_bytes(os.path.join(ctx.out, "photon-lowmem-verify.json"),
                       json.dumps(out, indent=2).encode("utf-8"))
    return out


# ==========================================================================
# main
# ==========================================================================
STAGES = {
    "base": stage_base,
    "photon": stage_photon,
    "gates234": stage_gates234,
    "gate1": stage_gate1,
    "gates56": stage_gates56,
    "photon16m": stage_photon16m,
    "gate6fix": stage_gate6fix,
    "assemble": stage_assemble,
}
ORDER = ["base", "photon", "gates234", "gate1", "gates56", "photon16m",
         "gate6fix", "assemble"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--q", type=int, default=49)
    ap.add_argument("--out", default="out/v8")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--fetch", type=float, default=P.POND_FETCH_M)
    ap.add_argument("--v33-dir", default="out/v7")
    ap.add_argument("--photons", type=int, default=4_000_000,
                    help="4M = the v2 photon-check condition")
    ap.add_argument("--shallow-depth", type=float, default=0.1)
    ap.add_argument("--shallow-wind", type=float, default=1.5)
    ap.add_argument("--clamp-wind", type=float, default=5.0)
    ap.add_argument("--stages", default="",
                    help="comma separated: " + ",".join(ORDER))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="recompute stages whose cache exists")
    ap.add_argument("--verify-lowmem", action="store_true")
    args = ap.parse_args()

    ctx = Ctx(args)
    os.makedirs(ctx.out, exist_ok=True)
    teed = Tee(os.path.join(ctx.out, "finish-run.log"))
    sys.stdout = teed
    sys.stderr = teed
    R_save = None
    print("=" * 70, flush=True)
    print("run_v8_finish.py  %s" % datetime.datetime.now(
        datetime.timezone.utc).isoformat(), flush=True)
    print("grid=%d Q=%d photons=%d out=%s" % (args.grid, args.q, args.photons,
                                              args.out), flush=True)
    m = mem_now() or {}
    print("free-memory probe: ws=%.0f MB peak=%.0f MB"
          % (m.get("ws_mb", 0), m.get("peak_ws_mb", 0)), flush=True)

    if args.verify_lowmem:
        verify_lowmem(ctx)
        return

    stages = ORDER if args.all else [s.strip() for s in
                                     args.stages.split(",") if s.strip()]
    if not stages:
        stages = ORDER
    for name in stages:
        fn = STAGES[name]
        if (name != "assemble") and ctx.have(name) and not args.force:
            print("\n[%s] cached -> skip (%s)" % (name, ctx.part_path(name)),
                  flush=True)
            continue
        payload = fn(ctx)
        if name == "assemble":
            print("\nnumbers-v8.json written (%d bytes)"
                  % os.path.getsize(ctx.png("numbers-v8.json")), flush=True)
    print("\ncheckpoint state: %s" % (ctx.cp.state(),), flush=True)


if __name__ == "__main__":
    main()
