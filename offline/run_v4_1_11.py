"""run_v4_1_11 -- BRIEF-V4-1-11: the GF / nu groupiness diagnostic
(measurement only) for the four recorded scenes.

What this job is
----------------
It measures; it does not synthesise a new scene.  Every field it measures is
rebuilt deterministically from the recorded seed + the recorded generator
(run_v19_prod / run_v20_prod / run_prod_v4+state.py / run_v4_1_8+
stream_surface3), and STAGE 1 proves that the rebuilt field IS the recorded
one before a single groupiness number is quoted.

Measured quantities, per scene and per probe (16x16, on state grid nodes):
    height      eta(t)        -> Hilbert envelope -> SIWEH -> GF_eta, nu_eta
                                 (+ skewness / kurtosis)
    slope       |grad eta|(t) -> GF_slope  (and the gx / gy lanes separately)
    curvature   lap eta(t)    -> GF_curv, nu_curv

GF_S = std_W(A^2) / mean_W(A^2),   A = |Hilbert(s)|,
       W = circular boxcar of one peak period (sensitivity 0.5/1/2 T_p).
GF_A = std(A) / mean(A)  -- the envelope-ratio convention, reported alongside
       so that a reader can map to the literature's numbers.

nu = sqrt(m0*m2/m1^2 - 1)  (Longuet-Higgins width) from the discrete frequency
moments of the probe record; the theoretical (un-quantised omega) value is
reported next to it because this generator quantises omega to multiples of
2*pi/8 s.

Stages (all foreground, checkpoint written per stage):
    --stage0   preflight, generator hashes vs the recorded stage-0 tables
    --stage1   rebuild the four scenes + the reproduction gate (STOP gate)
    --stage2   the GF / nu sweep + sensitivities + controls
    --stage3   numbers-v4-1-11.json + REPORT-V4-1-11.md

Exit code 0 only if every gate passes.
"""
import argparse
import datetime
import gc
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-11.md"
OUT = os.path.join(HERE, "out", "v4-1-11")
WORK = os.path.join(OUT, "_work")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-1-11.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-11.md")
GEN = os.path.join(HERE, "run_v4_1_11.py")

FPS = 24
SEED_BRIEF = 3255045276382562994
OMEGA0 = 2.0 * math.pi / 8.0

NPROBE = 16                     # 16 x 16 probes, edges excluded
SMOOTH_TP = 1.0                 # primary SIWEH window, in peak periods
SMOOTH_LANES = (0.5, 1.0, 2.0)
N_ENSEMBLE = 400                # random-phase surrogate draws
SPACES = ("height", "slope", "curvature")

SCENES = [
    {"id": "WIND-A", "label": "WIND-A (v19)",
     "kind": "wind", "U": 6.0, "F": 100.0, "h": 4.0,
     "seed": SEED_BRIEF, "nframes": 192, "patch": 4.8,
     "gen": "run_v19_prod.py", "record": "out/v19/numbers-v19-prod.json",
     "recorded": {"Hs_spectral_m": 0.02340958139785273,
                  "sigma2_spectral": 0.020570696421876894,
                  "rms_slope": 0.14342488076298648}},
    {"id": "WIND-B", "label": "WIND-B (v20)",
     "kind": "wind", "U": 3.0, "F": 100.0, "h": 0.6,
     "seed": SEED_BRIEF, "nframes": 192, "patch": 2.4,
     "gen": "run_v20_prod.py", "record": "out/v20/numbers-v20-prod.json",
     "recorded": {"Hs_spectral_m": 0.0130268561297069,
                  "sigma2_spectral": 0.011486073753957202,
                  "rms_slope": 0.10717310182110622}},
    {"id": "WIND-C", "label": "WIND-C (v2 loop)",
     "kind": "v2", "U": None, "F": None, "h": 0.6,
     "seed": SEED_BRIEF, "nframes": 192, "patch": 1.4755102040816326,
     "gen": "run_prod_v4.py + state.py", "record": "out/prod/numbers-prod.json",
     "recorded": {"patch_S_m": 1.4755102040816326,
                  "dx_m": 0.0007204639668367347}},
    {"id": "STREAM", "label": "STREAM (v4-1.8)",
     "kind": "stream", "U": 0.8, "F": None, "h": 0.6,
     "seed": 20261112, "nframes": 48, "patch": 1.4755102,
     "gen": "run_v4_1_8.py + stream_surface3.py",
     "record": "out/v4-1-8/numbers-v4-1-8.json",
     "recorded": {"g_mean_px": 10.067451730266395,
                  "g_p50_px": 9.415957692929581,
                  "g_p99_px": 23.775464767253553}},
]

HASH_TABLES = ["out/v19/stage0-py-hashes.json",
               "out/v20/stage0-py-hashes.json",
               "out/v4-1-8/stage0-py-hashes.json"]

#: the generator files the brief names + everything the rebuilds import
CODE_FILES = ["run_v19_prod.py", "run_v20_prod.py", "run_prod_v4.py",
              "run_v4_1_8.py", "stream_surface3.py", "stream_surface2.py",
              "state.py", "env.py", "det.py", "image.py",
              "run_v10.py", "run_v15.py", "run_v15_state.py",
              "run_v15_spectrum.py", "v10_state.py", "v10_spectrum.py",
              "v11_state.py", "v11_spectrum.py", "v12_water.py",
              "run_v16_water.py", "v5_composite.py", "run_v14_state.py",
              "run_v14_spectrum.py", "v13_state.py", "transport.py",
              "recipe_v8.py", "v8_physics.py", "run_v4_1_2.py",
              "stream_surface4.py"]

BRIEF_FILE_ROWS = [("run_v19_prod.py", "WIND-A"), ("run_v20_prod.py", "WIND-B"),
                   ("run_prod_v4.py + state.py", "WIND-C"),
                   ("run_v4_1_8.py + stream_surface3.py", "STREAM")]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
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


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path, obj):
    payload = json.dumps(obj, ensure_ascii=False, indent=2,
                         allow_nan=True).encode("utf-8")
    atomic_write_bytes(path, payload)
    return path


def cp_extra(key, value=None):
    state = read_json(CP_PATH) if os.path.exists(CP_PATH) else {
        "completed_units": 0, "completed_ids": [], "last_unit": None,
        "updated_at": None}
    if value is not None:
        state[key] = value
        state["updated_at"] = now_iso()
        atomic_write_bytes(CP_PATH,
                           json.dumps(state, ensure_ascii=False,
                                      indent=2).encode("utf-8"))
    return state


def py_set():
    out = {}
    for name in sorted(os.listdir(HERE)):
        p = os.path.join(HERE, name)
        if name.endswith(".py") and os.path.isfile(p):
            out[name] = sha256_file(p)
    return out


def rel_err(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    d = max(abs(a), abs(b))
    return float(abs(a - b) / d) if d > 0 else float(abs(a - b))


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight + provenance of the generators
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] preflight: files, generator hashes vs the recorded tables")
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    cp = Checkpoint(CP_PATH)

    facts = {"workdir": HERE, "brief": BRIEF, "created_utc": now_iso(),
             "python": sys.version.split()[0], "numpy": np.__version__,
             "brief_exists": os.path.exists(BRIEF)}
    try:
        import scipy
        facts["scipy"] = scipy.__version__
    except Exception:                                    # pragma: no cover
        facts["scipy"] = None

    facts["code"] = {}
    missing = []
    for f in CODE_FILES:
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            facts["code"][f] = {"sha256": sha256_file(p), "exists": True}
        else:
            missing.append(f)
            facts["code"][f] = {"sha256": None, "exists": False}
    facts["code_missing"] = missing

    recorded = {}
    for tab in HASH_TABLES:
        p = os.path.join(HERE, tab)
        recorded[tab] = read_json(p) if os.path.exists(p) else None

    def table_hashes(d):
        if not d:
            return {}
        return d.get("hashes") or d.get("all") or {}

    prov_rows = []
    diff = []
    for f in CODE_FILES:
        row = {"file": f, "sha256_now": facts["code"][f]["sha256"]}
        verdicts = {}
        for tab, d in recorded.items():
            h = table_hashes(d).get(f)
            tag = ("not-listed" if h is None
                   else ("match" if h == row["sha256_now"] else "DIFF"))
            verdicts[tab] = {"sha256_recorded": h, "verdict": tag}
            if tag == "DIFF":
                diff.append({"file": f, "table": tab})
        row["recorded"] = verdicts
        prov_rows.append(row)
    # a job that edits itself records the FINAL hash in its own numbers file:
    # compare against that where the stage-0 table is only the opening snapshot.
    facts["self_hashes"] = {}
    for jid, nfile, key in (
            ("v4-1-8", "out/v4-1-8/numbers-v4-1-8.json", "run_v4_1_8_sha256"),
            ("v19", "out/v19/numbers-v19-prod.json", "run_v19_prod_sha256"),
            ("v20", "out/v20/numbers-v20-prod.json", "run_v20_prod_sha256")):
        p = os.path.join(HERE, nfile)
        if not os.path.exists(p):
            continue
        d = read_json(p)
        h = d.get(key)
        if h:
            facts["self_hashes"][jid] = {
                "file": key.replace("_sha256", "") + ".py",
                "recorded_final_sha256": h,
                "matches_now": bool(h == facts["code"].get(
                    key.replace("_sha256", "") + ".py", {}).get("sha256")),
                "source": nfile}

    facts["provenance_rows"] = prov_rows
    facts["hash_diff"] = diff

    need = ["out/v19/numbers-v19-prod.json", "out/v20/numbers-v20-prod.json",
            "out/prod/numbers-prod.json", "out/v4-1-8/numbers-v4-1-8.json",
            "out/v19/_work/params.json", "out/v20/_work/params.json",
            "out/prod/flow-0000.npy", "out/v4-1-8/_work/warp-frame-0024.json"]
    facts["records"] = {}
    for r in need:
        p = os.path.join(HERE, r)
        facts["records"][r] = ({"bytes": os.path.getsize(p)}
                               if os.path.exists(p) else False)
    facts["records_missing"] = [r for r in need
                                if not os.path.exists(os.path.join(HERE, r))]

    from state import derive_seed
    import run_v10 as R10
    raw = R10.load_recipe()
    derived = int(derive_seed(raw, 1))
    facts["seed"] = {"brief_literal": SEED_BRIEF, "derive_seed_raw_1": derived,
                     "agree": bool(derived == SEED_BRIEF)}

    from env import resolve_environment
    env2 = resolve_environment(raw)
    facts["v2_recipe_env"] = {k: env2.get(k) for k in
                             ("windSpeedMS", "depthM", "specL", "waveHeightM",
                              "clarity") if k in env2}
    facts["v2_recipe_note"] = (
        "BRIEF-V4-1-11's WIND-C row carries WIND-B's line (U=3, F=100, "
        "elfouhaily-form).  The RECORDED v2 loop (out/prod/numbers-prod.json, "
        "pipeline 'v4') was generated by run_prod_v4.py -> state.WaterState "
        "(Phillips spectrum) on recipe-1.json with windSpeedMS=%g, "
        "specL=%.10f and a sigma_eta target of %.10f m.  The recorded field is "
        "that one, so that is what STAGE 1 rebuilds; the report states the "
        "difference from the brief's table."
        % (env2.get("windSpeedMS", float("nan")),
           env2.get("specL", float("nan")),
           env2.get("waveHeightM", float("nan"))))

    facts["stage0_wall_s"] = float(time.time() - t0)
    cp_extra("stage0", facts)
    cp.unit_done("stage0-preflight")
    log("  code files: %d present, %d missing"
        % (len(CODE_FILES) - len(missing), len(missing)))
    log("  hash DIFF count %d: %s" % (len(diff),
                                      sorted({d["file"] for d in diff})))
    log("  seed %d == derive_seed(recipe,1)? %s"
        % (SEED_BRIEF, facts["seed"]["agree"]))
    log("  stage0 wall %.1fs" % facts["stage0_wall_s"])
    if missing or facts["records_missing"]:
        raise SystemExit("STOP: missing inputs %s / %s"
                         % (missing, facts["records_missing"]))
    return facts


# ---------------------------------------------------------------------------
# STAGE 1 -- rebuild + reproduction gate
# ---------------------------------------------------------------------------
def _dispersion(k, h):
    from v8_physics import dispersion_omega
    return dispersion_omega(np.maximum(np.asarray(k, dtype=np.float64), 1e-12),
                            float(h))


def _omega_to_k(w, h):
    w = float(w)
    lo, hi = 1e-6, 5.0e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if float(_dispersion(mid, h)) < w:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


class WindProvider:
    """h / gx / gy / lap on the scene grid, single precision, per frame.

    Replicates V10WaterState.fields_at for the three quantities the diagnostic
    needs (it skips Dx/Dy, which the diagnostic does not use) and adds the
    spectral Laplacian of the COMPOSITE scene field.  The hi cascade is folded
    exactly as _fold_hi does.
    """

    periodic = True

    def __init__(self, st):
        from scipy.fft import ifft2, fft2
        self._ifft2, self._fft2 = ifft2, fft2
        self.st = st
        self.S = float(st.S)
        self.N = int(st.lo.N)
        self.j = int(getattr(st, "cascade_j", 0))
        self.b = 2 ** self.j if self.j > 0 else 1
        self.casc = [c for c in (st.lo, st.hi) if c is not None]
        self.info = []
        for c in self.casc:
            n_idx = np.rint(c.omega_q.astype(np.float64) / OMEGA0
                            ).astype(np.int32)
            n_idx[~np.isfinite(c.omega_q)] = 0
            self.info.append({
                "n_idx": n_idx,
                "nmax": int(n_idx.max()) if n_idx.size else 1,
                "ikx": (1j * c.kx.astype(np.float64)).astype(np.complex64),
                "iky": (1j * c.ky.astype(np.float64)).astype(np.complex64),
                "H0": c.H0, "N": int(c.N),
                "fold": bool(c.tag == "hi" and self.b > 1),
                "active": int((c.H0 != 0).sum())})
        self.k2 = (st.lo.k.astype(np.float64) ** 2).astype(np.float32)

    def _fold(self, a):
        b = self.b
        if b == 1:
            return a
        n = a.shape[0] // b
        return np.tile(a.reshape(n, b, n, b).mean(axis=(1, 3)).astype(a.dtype),
                       (b, b))

    def fields(self, t):
        h = np.zeros((self.N, self.N), dtype=np.float32)
        gx = np.zeros_like(h)
        gy = np.zeros_like(h)
        for inf in self.info:
            lut = np.exp(1j * OMEGA0 * float(t)
                         * np.arange(inf["nmax"] + 1)).astype(np.complex64)
            H = inf["H0"] * np.take(lut, inf["n_idx"])
            del lut
            hh = self._ifft2(H).real.astype(np.float32)
            gxx = self._ifft2(H * inf["ikx"]).real.astype(np.float32)
            gyy = self._ifft2(H * inf["iky"]).real.astype(np.float32)
            del H
            if inf["fold"]:
                hh, gxx, gyy = self._fold(hh), self._fold(gxx), self._fold(gyy)
            h += hh
            gx += gxx
            gy += gyy
            del hh, gxx, gyy
        lap = self._ifft2(self._fft2(h) * self.k2).real.astype(np.float32)
        return {"h": h, "gx": gx, "gy": gy, "lap": lap, "S": self.S}


class V2Provider:
    """The v2 loop state (state.WaterState, Phillips).  Same four quantities."""

    periodic = True

    def __init__(self, st):
        from scipy.fft import ifft2, fft2
        self._ifft2, self._fft2 = ifft2, fft2
        self.S = float(st.S)
        self.N = int(st.N)
        self.n_idx = np.rint(st._omega_q() / OMEGA0).astype(np.int32)
        self.nmax = int(self.n_idx.max())
        self.H0 = st.H0
        self.ikx = (1j * st.kx.astype(np.float64)).astype(np.complex64)
        self.iky = (1j * st.ky.astype(np.float64)).astype(np.complex64)
        self.k2 = (st.k.astype(np.float64) ** 2).astype(np.float32)

    def fields(self, t):
        lut = np.exp(1j * OMEGA0 * float(t)
                     * np.arange(self.nmax + 1)).astype(np.complex64)
        H = self.H0 * np.take(lut, self.n_idx)
        del lut
        h = self._ifft2(H).real.astype(np.float32)
        gx = self._ifft2(H * self.ikx).real.astype(np.float32)
        gy = self._ifft2(H * self.iky).real.astype(np.float32)
        del H
        lap = self._ifft2(self._fft2(h) * self.k2).real.astype(np.float32)
        return {"h": h, "gx": gx, "gy": gy, "lap": lap, "S": self.S}


class StreamProvider:
    """v4-1.8 stream surface: the generator's own fields_at + spectral lap.

    The stream surface is NOT time periodic (its advection omega = kx*Us is not
    quantised), so the Hilbert transform is reflect padded (declared).
    """

    periodic = False

    def __init__(self, st):
        from scipy.fft import ifft2, fft2
        self._ifft2, self._fft2 = ifft2, fft2
        self.st = st
        self.S = float(st.S)
        self.N = int(st.N)
        self.k2 = (st.kappa.astype(np.float64) ** 2).astype(np.float32)

    def fields(self, t):
        f = self.st.fields_at(float(t))
        h = np.asarray(f["h"], dtype=np.float32)
        lap = self._ifft2(self._fft2(h) * self.k2).real.astype(np.float32)
        return {"h": h, "gx": np.asarray(f["gx"], dtype=np.float32),
                "gy": np.asarray(f["gy"], dtype=np.float32),
                "lap": lap, "S": self.S}


def probe_index(N, n=NPROBE):
    step = N // n
    return (np.arange(n) * step + step // 2).astype(np.int64)


def build_scene(sc):
    """(provider, state, env, facts)."""
    import run_v10 as R10
    t0 = time.time()
    raw = R10.load_recipe()
    facts = {"id": sc["id"], "kind": sc["kind"], "nframes": int(sc["nframes"])}
    if sc["kind"] == "wind":
        import run_v15 as R15
        r = R15.scene_for(raw, sc["U"], sc["F"], sc["h"], 2048, sc["seed"])
        st = r["state"]
        prov = WindProvider(st)
        env = r["env"]
        facts.update({
            "U": sc["U"], "F": sc["F"], "depth_m": sc["h"],
            "patch_m": float(st.S), "grid": int(st.lo.N),
            "cascade_j": int(st.cascade_j),
            "patch_doublings": int(st.patch_doublings),
            "k_x_eff": float(st.k_x_eff), "kmin": float(st.lo.kmin),
            "kmax_nyq_scene": float(st.lo.kmax_nyq), "k_cut": float(st.k_cut),
            "active_modes": int(sum(i["active"] for i in prov.info)),
            "state_summary": {k: v for k, v in st.spectrum_after.items()
                              if isinstance(v, (int, float))},
            "lambda_p_m": float(st.laws["lambda_p_m"]),
            "omega_c": float(st.ctx["omega_c"])})
    elif sc["kind"] == "v2":
        from env import resolve_environment
        from state import WaterState
        env = resolve_environment(raw)
        st = WaterState(raw, env, N=2048, seed=sc["seed"], preset_index=1)
        prov = V2Provider(st)
        facts.update({"U": float(env["windSpeedMS"]), "F": None,
                      "depth_m": float(env["depthM"]),
                      "patch_m": float(st.S), "grid": int(st.N),
                      "cascade_j": 0, "patch_doublings": None,
                      "kmin": float(st.kmin),
                      "kmax_nyq_scene": float(st.kmax_nyq),
                      "active_modes": int((st.H0 != 0).sum()),
                      "target_std_m": float(st.target_std),
                      "specL_m": float(st.specL),
                      "spectrum": "Phillips (state.WaterState)"})
    elif sc["kind"] == "stream":
        import run_v15 as RV15
        import run_v4_1_2 as R12
        import stream_surface2 as SS
        import stream_surface3 as S3
        from state import derive_seed
        wseed = int(derive_seed(raw, 1))
        env = RV15.scene_for(raw, R12.WIND_CELL["U"], R12.WIND_CELL["F"],
                             float(sc["h"]), 256, wseed)["env"]
        env["depthM"] = float(sc["h"])
        st = S3.StreamState3(N=512, L=float(sc["patch"]), H=float(sc["h"]),
                             U=float(sc["U"]), S0=0.002, D84=0.09,
                             seed=int(sc["seed"]), bed="natural",
                             alpha=SS.ALPHA_MAIN, dt_loop=2.0)
        prov = StreamProvider(st)
        facts.update({"U": float(sc["U"]), "U0": float(st.Us), "F": None,
                      "depth_m": float(sc["h"]), "patch_m": float(st.S),
                      "grid": int(st.N), "k0": float(st.k0 or 0.0),
                      "ring_radius_rad_per_m": float(st.k0 or 0.0),
                      "phi_ring": float(st.phi),
                      "sigma_total_m": float(st.sigma_total),
                      "sigma_reg_m": float(st.sigma_reg),
                      "T_seam_s": float(st.S / st.Us),
                      "accounting": st.accounting()})
    else:
        raise SystemExit("unknown scene kind %r" % (sc["kind"],))
    facts["build_wall_s"] = float(time.time() - t0)
    return prov, st, env, facts


def warp_green_px(st, env, frame=24):
    """The v4-1-8 G-WARP yardstick, verbatim (run_v4_1_8.do_warp)."""
    import run_v4_1_2 as R12
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    fields = st.fields_at(frame / float(FPS))
    S = float(st.S)
    gx_s, gy_s = _state_gradients(fields, float(R12.B_SIGMA_FROZEN))
    iy, ix = np.mgrid[0:R12.HEIGHT, 0:R12.WIDTH]
    u = (ix.astype(np.float64) + 0.5) / R12.WIDTH
    v = (iy.astype(np.float64) + 0.5) / R12.HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S, S
                        ).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S, S
                        ).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    depth = float(env["depthM"])
    cb = (float(IOR[1]) - 1.0) / float(IOR[1])
    shx, shy = -nx * depth * cb, -ny * depth * cb
    dx_px = shx / S * (R12.WIDTH - 1.0)
    dy_px = shy / S * (R12.HEIGHT - 1.0)
    mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
    return {"g_mean_px": float(np.mean(mag)),
            "g_p50_px": float(np.percentile(mag, 50)),
            "g_p99_px": float(np.percentile(mag, 99)),
            "g_max_px": float(np.max(mag))}


def stage1():
    t0 = time.time()
    log("=" * 78)
    log("[stage 1] rebuild from the recorded seed + the reproduction gate")
    cp = Checkpoint(CP_PATH)
    gates, inst, checks = {}, {}, {}
    for sc in SCENES:
        log("  building %s (%s)" % (sc["id"], sc["gen"]))
        prov, st, env, facts = build_scene(sc)
        gate = {"id": sc["id"], "rows": [], "pass": True,
                "generator": sc["gen"], "record": sc["record"]}
        if sc["kind"] == "wind":
            ss = st.spectrum_after
            keymap = {"Hs_spectral_m": "Hs_m",
                      "sigma2_spectral": "sigma2_rad2",
                      "rms_slope": "rms_slope"}
            for key, rec in sc["recorded"].items():
                got = float(ss[keymap.get(key, key)])
                e = rel_err(got, rec)
                gate["rows"].append({"quantity": key, "recorded": rec,
                                     "rebuilt": got, "rel_err": e,
                                     "source_field": keymap.get(key, key),
                                     "tol": 1e-3, "pass": bool(e <= 1e-3)})
            gate["source"] = ("st.spectrum_after (the spectral readout of the "
                              "realised mode set) vs " + sc["record"])
        elif sc["kind"] == "v2":
            for key, got in (("patch_S_m", float(st.S)),
                             ("dx_m", float(st.dx))):
                e = rel_err(got, sc["recorded"][key])
                gate["rows"].append({"quantity": key,
                                     "recorded": sc["recorded"][key],
                                     "rebuilt": got, "rel_err": e,
                                     "tol": 1e-3, "pass": bool(e <= 1e-3)})
            from det import flow_backward_master
            flow = flow_backward_master(st.fields_at(0.0), float(st.S), 2048,
                                        env)
            rec_flow = np.load(os.path.join(HERE, "out", "prod",
                                            "flow-0000.npy"))
            d = np.abs(flow.astype(np.float64) - rec_flow.astype(np.float64))
            mx = float(np.abs(rec_flow).max())
            e = float(d.max() / mx) if mx > 0 else float(d.max())
            gate["rows"].append({
                "quantity": "flow_backward_master(t=0) vs out/prod/flow-0000.npy",
                "recorded": "array 2048x2048x2 float32 (out/prod/flow-0000.npy)",
                "rebuilt": "array 2048x2048x2 float32",
                "max_abs_diff": float(d.max()),
                "max_abs_diff_over_max_abs_record": e,
                "rel_err": e, "tol": 1e-3, "pass": bool(e <= 1e-3)})
            del flow, rec_flow, d
            gate["source"] = ("recorded patch/dx (numbers-prod.json 'frozen') "
                              "+ the recorded transport flow array at t=0")
        elif sc["kind"] == "stream":
            w = warp_green_px(st, env, frame=24)
            for key, rec in sc["recorded"].items():
                got = float(w[key])
                e = rel_err(got, rec)
                gate["rows"].append({"quantity": key, "recorded": rec,
                                     "rebuilt": got, "rel_err": e,
                                     "tol": 1e-3, "pass": bool(e <= 1e-3)})
            gate["source"] = ("run_v4_1_8 do_warp at frame 24 (t=1.0 s) vs "
                              + sc["record"])
        gate["pass"] = all(r["pass"] for r in gate["rows"])
        gates[sc["id"]] = gate

        # cross-check of the lean field provider against the canonical one
        if sc["kind"] != "stream":
            t = 8.0 / 24.0
            f0 = st.fields_at(t)
            f1 = prov.fields(t)
            iy = probe_index(int(f0["h"].shape[0]))
            def at(a):
                return np.asarray(a, dtype=np.float64)[np.ix_(iy, iy)]
            d = {}
            for key in ("h", "gx", "gy"):
                a, b = at(f0[key]), at(f1[key])
                s = max(float(np.abs(a).std()), 1e-30)
                d[key] = {"max_abs_diff": float(np.abs(a - b).max()),
                          "rms_rel": float(np.sqrt(np.mean((a - b) ** 2)) / s)}
            h = np.asarray(f0["h"], dtype=np.float64)
            dx = float(st.S) / h.shape[0]
            fd = (np.roll(h, 1, 0) + np.roll(h, -1, 0) + np.roll(h, 1, 1)
                  + np.roll(h, -1, 1) - 4.0 * h) / (dx * dx)
            lap = at(f1["lap"])
            rec = fd[np.ix_(iy, iy)]
            d["lap_vs_fd"] = {
                "median_rel": float(np.median(np.abs(lap - rec)
                                              / max(np.abs(rec).std(), 1e-30))),
                "ratio_std": float(np.std(lap) / max(np.std(rec), 1e-30)),
                "note": ("spectral Laplacian of the sampled composite field vs "
                         "a 5-point finite difference of the generator's own "
                         "h; the gap is the FD discretisation, not a field "
                         "mismatch")}
            checks[sc["id"]] = d
            log("    provider check: h %.2e | gx %.2e | gy %.2e | lap std "
                "ratio %.4f" % (d["h"]["rms_rel"], d["gx"]["rms_rel"],
                                d["gy"]["rms_rel"], d["lap_vs_fd"]["ratio_std"]))
            del f0, f1
        inst[sc["id"]] = facts
        for r in gate["rows"]:
            log("    %-50s rel_err %.3e %s"
                % (str(r["quantity"])[:50], r["rel_err"],
                   "ok" if r["pass"] else "NG"))
        log("    gate %s: %s" % (sc["id"], "PASS" if gate["pass"] else "FAIL"))
        del prov, st, env
        gc.collect()

    ok = all(g["pass"] for g in gates.values()) and len(gates) == 4
    out = {"gates": gates, "provider_checks": checks, "pass": bool(ok),
           "instances": inst, "wall_s": float(time.time() - t0)}
    cp_extra("stage1", out)
    cp.unit_done("stage1-repro-gate")
    log("  STAGE 1 %s (wall %.1fs)" % ("PASS" if ok else "FAIL", out["wall_s"]))
    if not ok:
        raise SystemExit("STOP: the rebuild is not the recorded field, so no "
                         "measurement is meaningful")
    return out


# ---------------------------------------------------------------------------
# the measurement primitives
# ---------------------------------------------------------------------------
def hilbert_batch(x, pad="reflect"):
    """analytic signal along axis 1 (periodic or reflect-padded)."""
    x = np.asarray(x, dtype=np.float64)
    if pad == "reflect":
        T = x.shape[1]
        p = T // 2
        xp = np.concatenate([x[:, 1:p + 1][:, ::-1], x,
                             x[:, -p - 1:-1][:, ::-1]], axis=1)
    else:
        xp = x
    T = xp.shape[1]
    X = np.fft.fft(xp, axis=1)
    m = np.zeros(T, dtype=np.float64)
    m[0] = 1.0
    if T % 2 == 0:
        m[1:T // 2] = 2.0
        m[T // 2] = 1.0
    else:
        m[1:(T + 1) // 2] = 2.0
    z = np.fft.ifft(X * m[None, :], axis=1)
    if pad == "reflect":
        p = x.shape[1] // 2
        z = z[:, p:p + x.shape[1]]
    return z


def smooth_circ(x, win, kind="box"):
    """periodic smoothing along axis 1 (boxcar or Gaussian)."""
    T = x.shape[1]
    w = int(round(win))
    if w <= 1:
        return x.copy()
    if kind == "box":
        ker = np.zeros(T)
        ker[(np.arange(w) - w // 2) % T] = 1.0 / w
    else:
        sig = float(w) / 2.0
        i = np.arange(T)
        d = np.minimum(i, T - i)
        ker = np.exp(-0.5 * (d / sig) ** 2)
        ker /= ker.sum()
    X = np.fft.fft(x, axis=1)
    K = np.fft.fft(ker)
    return np.real(np.fft.ifft(X * K[None, :], axis=1))


def gf_batch(x, win, kind="box", pad="reflect"):
    """GF_S (SIWEH form) and GF_A (envelope-ratio form) per row."""
    A = np.abs(hilbert_batch(x, pad=pad))
    A2 = A * A
    S = smooth_circ(A2, win, kind)
    gf_s = np.std(S, axis=1, ddof=0) / np.mean(S, axis=1)
    gf_a = np.std(A, axis=1, ddof=0) / np.mean(A, axis=1)
    return gf_s, gf_a


def nu_batch(x, dt):
    """Longuet-Higgins width from the discrete one-sided moments."""
    x = np.asarray(x, dtype=np.float64)
    T = x.shape[1]
    X = np.fft.fft(x, axis=1)
    w = 2.0 * math.pi * np.fft.fftfreq(T, d=dt)
    P = (np.abs(X[:, 1:T // 2]) ** 2)
    wp = w[1:T // 2]
    m0 = P.sum(axis=1)
    m1 = (P * wp[None, :]).sum(axis=1)
    m2 = (P * wp[None, :] ** 2).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = m0 * m2 / (m1 * m1) - 1.0
    nu = np.sqrt(np.maximum(np.where(np.isfinite(v), v, 0.0), 0.0))
    return {"nu": nu, "omega_peak": wp[np.argmax(P, axis=1)]}


def skew_kurt(x):
    d = x - x.mean(axis=1)[:, None]
    s = d.std(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.mean(d ** 3, axis=1) / s ** 3,
                np.mean(d ** 4, axis=1) / s ** 4)


def q_stats(v):
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"median": None, "q25": None, "q75": None, "n": 0}
    return {"median": float(np.median(v)),
            "q25": float(np.percentile(v, 25)),
            "q75": float(np.percentile(v, 75)),
            "mean": float(np.mean(v)), "n": int(v.size)}


# ---------------------------------------------------------------------------
# STAGE 2 -- the GF / nu sweep
# ---------------------------------------------------------------------------
def spectrum_centroids(st, kind):
    """model-spectrum moments over the realised (pre-fold) mode set."""
    nyq = math.pi * FPS
    if kind == "stream":
        # forced pattern: (i) bed response is stationary, (ii) advects at kx*Us,
        # (iii) advects at kx*Us + omega0(kappa)
        import stream_surface2 as SS
        kk = st.kappa.astype(np.float64).ravel()
        kx = st.kx.astype(np.float64).ravel()
        om_i = np.zeros_like(kk)
        om_ii = kx * st.Us
        om_iii = kx * st.Us + SS.omega0_field(st.kappa, st.H).ravel()
        parts = [(np.abs(st.eta_i_hat) ** 2).ravel(), om_i,
                 (np.abs(st.hat_ii0) ** 2).ravel(), om_ii,
                 (np.abs(st.hat_turbA) ** 2).ravel(), om_iii]
        P = parts[0] + parts[2] + parts[4]
        out = {"k_centroid_height": math.sqrt(
                   float((P * kk ** 2).sum()) / max(float(P.sum()), 1e-30)),
               "k_centroid_slope": math.sqrt(
                   float((P * kk ** 4).sum())
                   / max(float((P * kk ** 2).sum()), 1e-30)),
               "k_centroid_curv": math.sqrt(
                   float((P * kk ** 6).sum())
                   / max(float((P * kk ** 4).sum()), 1e-30))}
        shares = {}
        for tag, jj in (("height", 0), ("slope", 2), ("curv", 4)):
            num = den = 0.0
            for a in range(0, len(parts), 2):
                wt = parts[a] * kk ** jj
                den += float(wt.sum())
                num += float(wt[np.abs(parts[a + 1]) > nyq].sum())
            shares[tag] = (num / den) if den > 0 else None
        out["aliased_share_above_temporal_nyquist"] = shares
        out["temporal_nyquist_rad_s"] = nyq
        out["note"] = ("stream: temporal frequency taken as kx*Us (ii) and "
                       "kx*Us+omega0(kappa) (iii); the bed response (i) is "
                       "stationary (omega = 0)")
        return out

    pairs = []
    if kind == "v2":
        k = st.k.astype(np.float64).ravel()
        P = (np.abs(st.H0.astype(np.complex128)) ** 2).ravel()
        m = P > 0
        pairs = [(k[m], P[m])]
    else:
        for c in (st.lo, st.hi):
            if c is None:
                continue
            k = c.k.astype(np.float64).ravel()
            P = (np.abs(c.H0.astype(np.complex128)) ** 2).ravel()
            m = P > 0
            pairs.append((k[m], P[m]))
    if not pairs:
        return {}
    k = np.concatenate([a for a, _ in pairs])
    P = np.concatenate([b for _, b in pairs])
    m = {j: float((P * k ** j).sum()) for j in range(7)}
    out = {"moments_k_m0_m6": {str(a): b for a, b in m.items()}}
    if m[0] > 0 and m[2] > 0:
        out["k_centroid_height"] = math.sqrt(m[2] / m[0])
    if m[2] > 0 and m[4] > 0:
        out["k_centroid_slope"] = math.sqrt(m[4] / m[2])
    if m[4] > 0 and m[6] > 0:
        out["k_centroid_curv"] = math.sqrt(m[6] / m[4])

    nyq = math.pi * FPS
    shares = {}
    for tag, jj in (("height", 0), ("slope", 2), ("curv", 4)):
        num = den = 0.0
        if kind == "v2":
            cs = [(st.k.astype(np.float64).ravel(),
                   st._omega_q().ravel(),
                   (np.abs(st.H0.astype(np.complex128)) ** 2).ravel())]
        else:
            cs = [(c.k.astype(np.float64).ravel(),
                   np.rint(c.omega_q.astype(np.float64).ravel() / OMEGA0)
                   * OMEGA0,
                   (np.abs(c.H0.astype(np.complex128)) ** 2).ravel())
                  for c in (st.lo, st.hi) if c is not None]
        for kk, w, P in cs:
            w = np.where(np.isfinite(w), w, 0.0)
            wt = P * kk ** jj
            den += float(wt.sum())
            num += float(wt[w > nyq].sum())
        shares[tag] = (num / den) if den > 0 else None
    out["aliased_share_above_temporal_nyquist"] = shares
    out["temporal_nyquist_rad_s"] = nyq
    return out


def measure_scene(sc, facts, prov, st):
    t0 = time.time()
    N = int(facts["grid"])
    nf = int(sc["nframes"])
    dt = 1.0 / FPS
    iy = probe_index(N)
    P = NPROBE * NPROBE
    series = {k: np.zeros((P, nf), dtype=np.float32)
              for k in ("h", "gx", "gy", "lap")}
    seam = None
    lattice = None
    for m in range(nf):
        f = prov.fields(m * dt)
        for k in ("h", "gx", "gy", "lap"):
            series[k][:, m] = f[k][np.ix_(iy, iy)].ravel()
        if m == 0:
            hg = np.asarray(f["h"], dtype=np.float64)
            on = hg[np.ix_(iy, iy)]
            offs = {}
            for off in (16, 32, 48, 80):
                j = (iy + off) % N
                offs[str(off)] = float(np.sqrt(np.var(hg[np.ix_(j, j)])))
            lattice = {
                "rms_full_grid": float(np.sqrt(np.var(hg))),
                "rms_probe_lattice": float(np.sqrt(np.var(on))),
                "rms_probe_over_full": float(np.sqrt(np.var(on))
                                             / max(np.sqrt(np.var(hg)), 1e-30)),
                "rms_shifted_lattices": offs,
                "note": ("the 16x16 probe lattice is a sparse sample of a field "
                         "that holds only S/lambda_p wavelengths across the "
                         "patch; the rms over the lattice differs from the full "
                         "grid by this much because of that spatial sampling, "
                         "the same finite-record effect that makes GF a random "
                         "variable.  GF is a normalised statistic.")}
        if m == 0 and facts.get("cascade_j", 0) > 0:
            h = f["h"]
            n2 = h.shape[0] // 2
            jj = np.abs(h[:, n2] - h[:, n2 - 1]).mean()
            ref = np.abs(h[:, n2 // 4 + 7] - h[:, n2 // 4 + 6]).mean()
            seam = {"mean_jump_at_tile_seam_m": float(jj),
                    "mean_neighbour_gradient_m": float(ref),
                    "ratio": float(jj / ref) if ref > 0 else None,
                    "seam_index": int(n2),
                    "note": ("the capillary tile is block-averaged then TILED, "
                             "so the composite h has a replication seam; the "
                             "probe stride is %d and the seam sits at %d"
                             % (N // NPROBE, n2))}
        del f
        if (m + 1) % 48 == 0:
            log("    %s frame %d/%d" % (sc["id"], m + 1, nf))
    gc.collect()

    eta = series["h"].astype(np.float64)
    gx = series["gx"].astype(np.float64)
    gy = series["gy"].astype(np.float64)
    lap = series["lap"].astype(np.float64)
    gmag = np.sqrt(gx * gx + gy * gy)

    pad = "periodic" if prov.periodic else "reflect"
    wp = float(np.median(nu_batch(eta, dt)["omega_peak"]))
    Tp = 2.0 * math.pi / wp if wp > 0 else float("nan")
    win = max(1, int(round(SMOOTH_TP * Tp / dt)))

    def block(x):
        gs, ga = gf_batch(x, win, "box", pad=pad)
        nu = nu_batch(x, dt)["nu"]
        return {"gf_S": q_stats(gs), "gf_A": q_stats(ga), "nu": q_stats(nu),
                "per_probe_gf_S": [float(v) for v in gs],
                "per_probe_nu": [float(v) for v in nu]}

    spaces = {"height": block(eta), "slope": block(gmag),
              "curvature": block(lap)}
    sk, ku = skew_kurt(eta)
    spaces["height"]["skew"] = q_stats(sk)
    spaces["height"]["kurt"] = q_stats(ku)
    spaces["slope_components"] = {"gx": block(gx), "gy": block(gy)}
    spaces["slope"]["nu_source"] = ("nu from the along-wind component gx; the "
                                    "gy lane is reported next to it")
    spaces["slope"]["definition"] = ("GF on |grad eta|(t); components gx, gy "
                                     "reported separately")

    amp = {"rms_eta_measured": float(np.sqrt(np.mean(eta ** 2))),
           "rms_slope_measured": float(np.sqrt(np.mean(gmag ** 2))),
           "rms_lap_measured": float(np.sqrt(np.mean(lap ** 2))),
           "rms_slope_x": float(np.sqrt(np.mean(gx ** 2))),
           "rms_slope_y": float(np.sqrt(np.mean(gy ** 2)))}
    if sc["kind"] == "wind":
        amp["rms_eta_spectral"] = float(math.sqrt(st.spectrum_after["m0_m2"]))
        amp["rms_slope_spectral"] = float(st.spectrum_after["rms_slope"])
    elif sc["kind"] == "v2":
        amp["rms_eta_spectral"] = float(st.target_std)
        amp["rms_slope_spectral"] = float(math.sqrt(
            (st.k ** 2 * (np.abs(st.H0) ** 2) / (2.0 * st.N ** 4)).sum()))

    cent = spectrum_centroids(st, sc["kind"])
    for tag, key in (("height", "k_centroid_height"),
                     ("slope", "k_centroid_slope"),
                     ("curvature", "k_centroid_curv")):
        if cent.get(key):
            cent["lambda_%s_m" % tag] = 2.0 * math.pi / cent[key]
            cent["cells_%s" % tag] = ((2.0 * math.pi / cent[key])
                                     / (float(facts["patch_m"]) / N))

    window = {"patch_m": float(facts["patch_m"]),
              "depth_m": float(facts["depth_m"]),
              "U_ms": facts.get("U"), "F_m": facts.get("F"),
              "grid": N, "dx_mm": float(facts["patch_m"]) / N * 1e3,
              "frames": nf, "fps": FPS, "record_s": nf * dt,
              "probes": P,
              "probe_layout": "%dx%d nodes, stride %d cells, edges excluded"
                              % (NPROBE, NPROBE, N // NPROBE),
              "T_peak_s": Tp, "omega_peak_rad_s": wp,
              "periods_in_record": (nf * dt) / Tp,
              "smooth_window_s": win * dt,
              "smooth_window_samples": int(win),
              "smooth_window_periods": win * dt / Tp,
              "boundary": pad,
              "lambda_p_m": facts.get("lambda_p_m")}
    return {"id": sc["id"], "label": sc["label"], "kind": sc["kind"],
            "window": window, "spaces": spaces, "amplitudes": amp,
            "spectral": cent, "seam": seam, "lattice": lattice, "pad": pad,
            "dt": dt, "win": win, "Tp": Tp, "series": series,
            "wall_s": float(time.time() - t0)}


def sensitivities(m):
    """window / probe-count / rate / boundary lanes (this scene, all spaces)."""
    dt, Tp, pad = m["dt"], m["Tp"], m["pad"]
    spaces = {"height": m["series"]["h"].astype(np.float64),
              "slope": np.sqrt(m["series"]["gx"].astype(np.float64) ** 2
                               + m["series"]["gy"].astype(np.float64) ** 2),
              "curvature": m["series"]["lap"].astype(np.float64)}
    out = {}
    lanes = {}
    for tag, x in spaces.items():
        row = {}
        for f in SMOOTH_LANES:
            w = max(1, int(round(f * Tp / dt)))
            for kind in ("box", "gauss"):
                gs, _ = gf_batch(x, w, kind, pad=pad)
                row["%s_%gp" % (kind, f)] = float(np.median(gs))
        lanes[tag] = row
    out["smoothing"] = {"lanes": lanes, "T_peak_s": Tp,
                        "note": ("median GF_S over the %d probes; the boxcar "
                                 "of %g T_p is the primary definition"
                                 % (x.shape[0], SMOOTH_TP))}

    prows = {}
    w = max(1, int(round(SMOOTH_TP * Tp / dt)))
    for g in (16, 8, 4):
        step = 16 // g
        sel = np.array([(a * step) * 16 + (b * step)
                        for a in range(g) for b in range(g)])
        prows[str(g)] = {}
        for tag, x in spaces.items():
            gs, _ = gf_batch(x[sel], w, "box", pad=pad)
            prows[str(g)][tag] = float(np.median(gs))
    out["probe_count"] = prows

    rates = {}
    for f in (1, 2, 3):
        d2 = dt * f
        w = max(1, int(round(SMOOTH_TP * Tp / d2)))
        rates[str(FPS // f)] = {}
        for tag, x in spaces.items():
            gs, _ = gf_batch(x[:, ::f], w, "box", pad=pad)
            rates[str(FPS // f)][tag] = float(np.median(gs))
    out["rate"] = {"by_fps": rates,
                   "note": ("the primary record is 24 fps (the scene's own "
                            "rate); 12 / 8 fps are decimations, shown only for "
                            "the rate sensitivity")}

    if pad == "periodic":
        wp = max(1, int(round(SMOOTH_TP * Tp / dt)))
        b = {}
        for tag, x in spaces.items():
            g1, _ = gf_batch(x, wp, "box", pad="periodic")
            g2, _ = gf_batch(x, wp, "box", pad="reflect")
            b[tag] = {"periodic": float(np.median(g1)),
                      "reflect": float(np.median(g2)),
                      "rel_diff": float(abs(np.median(g1) - np.median(g2))
                                        / max(abs(np.median(g1)), 1e-12))}
        out["boundary"] = b
        out["boundary_note"] = ("the wind record is exactly one 8 s period, so "
                                "periodic and reflect boundary treatments of "
                                "the Hilbert transform agree to this level")
    return out


def controls(m):
    """per-space random-phase ensemble band + the phase-aligned control."""
    dt, win, pad = m["dt"], m["win"], m["pad"]
    rng = np.random.default_rng(20261113)
    out = {}
    for tag in SPACES:
        x = {"height": m["series"]["h"].astype(np.float64),
             "slope": np.sqrt(m["series"]["gx"].astype(np.float64) ** 2
                              + m["series"]["gy"].astype(np.float64) ** 2),
             "curvature": m["series"]["lap"].astype(np.float64)}[tag]
        P, T = x.shape
        Pn = np.mean(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=0)
        pos = np.arange(1, T // 2)
        pw = np.maximum(Pn[pos], 0.0)

        def build(aligned, rr=None):
            S = np.zeros(T, dtype=np.complex128)
            if aligned:
                amp = np.sqrt(pw / 2.0)
                S[pos] = amp
                S[-pos] = amp
            else:
                a = (rr.standard_normal(pos.size)
                     + 1j * rr.standard_normal(pos.size))
                amp = np.sqrt(pw / 2.0) * a
                S[pos] = amp
                S[-pos] = np.conj(amp)
            return np.real(np.fft.ifft(S))

        sur = np.array([build(False, rng) for _ in range(N_ENSEMBLE)])
        gs_s, ga_s = gf_batch(sur, win, "box", pad=pad)
        gs_a, ga_a = gf_batch(build(True)[None, :], win, "box", pad=pad)
        band = np.percentile(gs_s, [5, 50, 95])
        pp = np.array(m["spaces"][tag]["per_probe_gf_S"])
        out[tag] = {
            "ensemble": {"n": N_ENSEMBLE,
                         "gf_S_p5": float(band[0]), "gf_S_p50": float(band[1]),
                         "gf_S_p95": float(band[2]),
                         "obs_median": float(np.median(pp)),
                         "probe_fraction_inside_p5_p95": float(
                             np.mean((pp >= band[0]) & (pp <= band[2]))),
                         "gf_A_p5": float(np.percentile(ga_s, 5)),
                         "gf_A_p50": float(np.percentile(ga_s, 50)),
                         "gf_A_p95": float(np.percentile(ga_s, 95))},
            "phase_aligned": {"gf_S": float(gs_a[0]),
                              "gf_A": float(ga_a[0])}}
    # instrument ordering test: use tones with an INTEGER number of cycles in
    # the record so that a single tone is exactly periodic (GF must be 0 then)
    T = m["series"]["h"].shape[1]
    Ts = T * dt
    tt = np.arange(T) * dt
    ncyc = max(2, int(round(Ts / m["Tp"])))
    f1, f2 = ncyc / Ts, (ncyc + 2) / Ts
    gs_i, _ = gf_batch(np.sin(2 * math.pi * f1 * tt)[None, :], win, "box",
                       pad=pad)
    gs_b, _ = gf_batch((np.sin(2 * math.pi * f1 * tt)
                        + np.sin(2 * math.pi * f2 * tt))[None, :],
                       win, "box", pad=pad)
    out["instrument_tests"] = {"pure_sine_gf_S": float(gs_i[0]),
                              "two_tone_1.12_gf_S": float(gs_b[0]),
                              "sine_cycles_in_record": int(ncyc),
                              "two_tone_hz": [float(f1), float(f2)],
                              "note": ("ordering test of the estimator: a "
                                       "single periodic sinusoid must give ~0, "
                                       "a beat pair must give a large GF")}
    # convention check: a narrowband Gaussian record (the Rayleigh-envelope
    # class a natural sea belongs to).  Its GF_A is the number that should be
    # compared with the literature's ~0.57, not GF_S.
    xh = m["series"]["h"].astype(np.float64)
    nc = int(np.argmax(np.mean(np.abs(np.fft.fft(xh, axis=1)) ** 2,
                               axis=0))) + 1
    nc = max(1, min(nc, T // 2 - 21))
    rng2 = np.random.default_rng(7)
    nb = np.zeros((24, T))
    for i in range(24):
        S = np.zeros(T, dtype=np.complex128)
        for kk in range(nc, nc + 20):
            a = rng2.standard_normal() + 1j * rng2.standard_normal()
            S[kk] = a
            S[-kk] = np.conj(a)
        nb[i] = np.real(np.fft.ifft(S))
    gs_nb, ga_nb = gf_batch(nb, win, "box", pad=pad)
    out["narrowband_gaussian"] = {
        "gf_S_median": float(np.median(gs_nb)),
        "gf_A_median": float(np.median(ga_nb)),
        "bins": [int(nc), int(nc + 19)], "draws": 24,
        "note": ("20 adjacent frequency bins, random phases, 24 draws: the "
                 "Rayleigh-envelope class a natural sea belongs to.  GF_A is "
                 "the convention to compare with the literature's ~0.57")}
    out["note"] = ("the ensemble is a random-phase surrogate built from the "
                   "MEASURED record-mean single-point power spectrum of each "
                   "space: it is the sampling distribution of GF for ONE record "
                   "of this length and spectral shape, not a new scene")
    return out


def stage2(s1):
    t0 = time.time()
    log("=" * 78)
    log("[stage 2] GF / nu sweep: 4 scenes x 3 spaces, 16x16 probes")
    cp = Checkpoint(CP_PATH)
    out = {"scenes": {}, "sensitivity": {}, "controls": {},
           "definitions": {
               "gf_S": ("std_W(A^2)/mean_W(A^2), A=|Hilbert(s)|, W = circular "
                        "boxcar of %g peak period (population std)" % SMOOTH_TP),
               "gf_A": "std(A)/mean(A) with the same A (envelope-ratio form)",
               "nu": ("sqrt(m0*m2/m1^2-1) from the discrete one-sided "
                      "frequency moments of the probe record"),
               "probe": ("16x16 state-grid nodes, stride N/16, edges excluded "
                         "(the samples the transport stage reads)"),
               "smoothing": ("SIWEH = smoothed instantaneous wave energy; the "
                             "window is one peak period of the record"),
               "boundary": ("periodic Hilbert for the wind lanes (the record is "
                            "exactly one 8 s period by construction); reflect "
                            "padding for the stream lane (its advection omega "
                            "is not quantised)")}}
    for sc in SCENES:
        log("  measuring %s (%d frames, %d probes)"
            % (sc["id"], sc["nframes"], NPROBE * NPROBE))
        prov, st, env, facts = build_scene(sc)
        m = measure_scene(sc, facts, prov, st)
        out["controls"][sc["id"]] = controls(m)
        out["sensitivity"][sc["id"]] = sensitivities(m)
        series = m.pop("series")
        m["facts"] = {k: v for k, v in facts.items() if k != "accounting"}
        out["scenes"][sc["id"]] = m
        iy = probe_index(int(facts["grid"]))
        dx = float(facts["patch_m"]) / int(facts["grid"])
        np.savez_compressed(
            os.path.join(OUT, "probes-%s.npz" % sc["id"].lower()),
            iy=iy, x_m=(iy + 0.5) * dx,
            t=np.arange(int(sc["nframes"])) / float(FPS),
            eta=series["h"], gx=series["gx"], gy=series["gy"],
            lap=series["lap"])
        log("    %s: T_peak %.4f s | periods in record %.2f | window %d "
            "samples (%.3f s)"
            % (sc["id"], m["Tp"], m["window"]["periods_in_record"],
               m["win"], m["window"]["smooth_window_s"]))
        for tag in SPACES:
            b = m["spaces"][tag]
            log("      %-9s GF_S %.4f [%.4f-%.4f]  nu %.4f  GF_A %.4f"
                % (tag, b["gf_S"]["median"], b["gf_S"]["q25"],
                   b["gf_S"]["q75"], b["nu"]["median"], b["gf_A"]["median"]))
        del m, series, prov, st, env
        gc.collect()
    out["wall_s"] = float(time.time() - t0)
    cp_extra("stage2", out)
    cp.unit_done("stage2-gf-nu-sweep")
    log("  stage 2 wall %.1fs" % out["wall_s"])
    return out


# ---------------------------------------------------------------------------
# STAGE 3 -- numbers JSON + report
# ---------------------------------------------------------------------------
def bfi_numbers(s1, s2):
    out = {}
    for sc in SCENES:
        if sc["id"] == "STREAM":
            out["STREAM"] = {"note": "forced pattern; BFI is not defined here"}
            continue
        inst = s1["instances"][sc["id"]]
        meas = s2["scenes"][sc["id"]]
        Hs = 4.0 * meas["amplitudes"]["rms_eta_measured"]
        wp = meas["window"]["omega_peak_rad_s"]
        kp = _omega_to_k(wp, float(inst["depth_m"]))
        lam = 2.0 * math.pi / kp
        ka = kp * Hs / 2.0
        nuf = max(float(meas["spaces"]["height"]["nu"]["median"]), 1e-9)
        out[sc["id"]] = {
            "Hs": float(Hs), "lambda_p": float(lam), "k_p": float(kp),
            "ka": float(ka), "domega_over_omega": float(nuf),
            "BFI": float(ka / nuf),
            "group_len_1_over_nu": float(1.0 / nuf),
            "group_len_m": float(lam / nuf),
            "definition": ("BFI = k_p*a/(domega/omega_p) with a = Hs/2 and "
                           "domega/omega = nu (the Longuet-Higgins width)"),
            "lambda_p_source": ("omega_peak of the probe record inverted "
                                "through the generator dispersion relation"),
            "lambda_p_recorded": inst.get("lambda_p_m")}
    return out


def build_numbers(s0, s1, s2):
    numbers = {
        "job": "v4-1-11 GF/nu groupiness diagnostic (measurement only)",
        "brief": BRIEF, "workdir": HERE, "created_utc": now_iso(),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "scipy": s0.get("scipy"),
        "generator": GEN, "generator_sha256": sha256_file(GEN),
        "constants": {"FPS": FPS, "NPROBE": NPROBE, "SMOOTH_TP": SMOOTH_TP,
                      "SMOOTH_LANES": list(SMOOTH_LANES),
                      "N_ENSEMBLE": N_ENSEMBLE, "OMEGA0": OMEGA0},
        "stage0": s0, "stage1": s1, "stage2": s2,
        "bfi": bfi_numbers(s1, s2), "gates": {}}
    numbers["gates"]["G-REPRO"] = {
        "pass": bool(s1["pass"]),
        "criterion": "relative error <= 1e-3 on every recorded scalar / array",
        "detail": {k: [{"quantity": r["quantity"], "rel_err": r["rel_err"]}
                       for r in g["rows"]] for k, g in s1["gates"].items()}}

    rows = []
    complete = True
    for sc in SCENES:
        scd = s2["scenes"][sc["id"]]
        for tag in SPACES:
            b = scd["spaces"][tag]
            w = scd["window"]
            ens = s2["controls"][sc["id"]][tag]["ensemble"]
            ok = (all(b["gf_S"][k] is not None
                      for k in ("median", "q25", "q75"))
                  and all(b["nu"][k] is not None
                          for k in ("median", "q25", "q75"))
                  and bool(w["probes"]) and w["smooth_window_s"] > 0)
            complete = complete and ok
            rows.append({"scene": sc["id"], "space": tag,
                         "gf_S": b["gf_S"], "gf_A": b["gf_A"], "nu": b["nu"],
                         "ensemble": ens, "window_summary": {
                             "patch_m": w["patch_m"], "depth_m": w["depth_m"],
                             "U_ms": w["U_ms"], "F_m": w["F_m"],
                             "grid": w["grid"], "dx_mm": w["dx_mm"],
                             "frames": w["frames"], "record_s": w["record_s"],
                             "probes": w["probes"],
                             "T_peak_s": w["T_peak_s"],
                             "periods_in_record": w["periods_in_record"],
                             "smooth_window_s": w["smooth_window_s"],
                             "smooth_window_samples": w["smooth_window_samples"],
                             "boundary": w["boundary"]},
                         "complete": bool(ok)})
    numbers["gates"]["G-DIAG"] = {
        "pass": bool(complete and len(rows) == 12),
        "criterion": ("GF+nu complete for 4 scenes x 3 spaces; window / probes "
                      "/ smoothing on every row; literature comparison present"),
        "rows": rows, "n_rows": len(rows)}
    return numbers


def write_report(numbers):
    S1, S2, B = numbers["stage1"], numbers["stage2"], numbers["bfi"]
    G = numbers["gates"]

    def f(x, n=4):
        return "n/a" if x is None else ("%.*f" % (n, x))

    def pct(x, n=1):
        return "n/a" if x is None else ("%.*f%%" % (n, 100.0 * x))

    def band(sid, tag):
        e = S2["controls"][sid][tag]["ensemble"]
        return "[%s–%s]" % (f(e["gf_S_p5"], 3), f(e["gf_S_p95"], 3)), e

    L = []
    A = L.append
    A("# REPORT-V4-1-11 — GF/ν 診断：波の群れの物差し（測定のみ）")
    A("")
    A("Brief: `BRIEF-V4-1-11.md`（`D:\\mutaro-work\\rally-in\\BRIEF-V4-1-11.md`）  ")
    A("Workdir: `%s`  " % numbers["workdir"])
    A("生成: %s  " % numbers["created_utc"])
    A("新規ファイルは `run_v4_1_11.py` の 1 本のみ。**本ジョブが参照・再生成に使う "
      "29 本の .py は無改造**（sha256 は stage0 と stage3 で一致。全 .py の sha256 は "
      "`numbers-v4-1-11.json` の `stage0.code`）。")
    A("")

    # ---- 0 ---------------------------------------------------------------
    A("## 0. 判定")
    A("")
    A("| ゲート | 判定 | 判定に使った量 |")
    A("| --- | --- | --- |")
    maxrel = max(r["rel_err"] for v in G["G-REPRO"]["detail"].values()
                 for r in v)
    A("| G-REPRO | **%s** | 4 場面の再生成値 vs 記録値の相対誤差 ≤ 1e-3"
      "（最大 %.2e） |" % ("PASS" if G["G-REPRO"]["pass"] else "FAIL", maxrel))
    A("| G-DIAG | **%s** | GF＋ν が 4 場面 × 3 空間で完全（%d 行）、"
      "各行に窓・プローブ数・平滑化、文献対比あり |"
      % ("PASS" if G["G-DIAG"]["pass"] else "FAIL", G["G-DIAG"]["n_rows"]))
    A("| G-OUT | **%s** | REPORT / numbers / checkpoint が `out/v4-1-11/`、"
      "既存 .py 無改造 = %s |"
      % ("PASS" if G["G-OUT"]["pass"] else "FAIL",
         G["G-OUT"]["existing_py_untouched"]))
    A("")

    A("**問い 1: うちの各場面の群性は自然なランダム位相の範囲にあるか。**")
    A("")
    A("| 場面 | GF_S（高さ、中央値） | 同じ場面のランダム位相帯 5–95% | 判定 |")
    A("| --- | --- | --- | --- |")
    for sid in ("WIND-A", "WIND-B", "WIND-C", "STREAM"):
        b, e = band(sid, "height")
        obs = S2["scenes"][sid]["spaces"]["height"]["gf_S"]["median"]
        A("| %s | %s | %s | %s |"
          % (sid, f(obs), b,
             "帯の内" if e["gf_S_p5"] <= obs <= e["gf_S_p95"] else "帯の外"))
    A("")
    A("3 つの風波場面はいずれも、**同じ分散スペクトル・同じ記録長・同じ推定量**で"
      "作ったランダム位相の代理アンサンブル（400 本）の 5–95%% 帯の内側にある"
      "（プローブ 256 点のうち帯に入る割合は WIND-A %s、WIND-B %s、WIND-C %s）。"
      "つまり **うちの場面の群性は「そのスペクトルがランダム位相で作られたときに"
      "期待される範囲」から外れていない**。カメラ側で群れを足す必要はない。"
      "（STREAM も表の上では帯の内側だが、その帯は強制パターンの代理なので"
      "解釈は §7 で分ける。）"
      % (pct(S2["controls"]["WIND-A"]["height"]["ensemble"]
             ["probe_fraction_inside_p5_p95"]),
         pct(S2["controls"]["WIND-B"]["height"]["ensemble"]
             ["probe_fraction_inside_p5_p95"]),
         pct(S2["controls"]["WIND-C"]["height"]["ensemble"]
             ["probe_fraction_inside_p5_p95"])))
    A("")
    A("**問い 2: 0.95 は目標か。** 目標ではない。同一の分散スペクトル・同一記録長で"
      "**位相だけを t=0 に揃える**と（NewWave 型の焦点。これは対照であって場面では"
      "ない）、GF_S は観測値の WIND-A %s 倍 / WIND-B %s 倍 / WIND-C %s 倍へ跳び、"
      "同じ場面のランダム位相帯の上限（WIND-A %s）を大きく超える。つまり"
      "**位相集中はこの推定量で明確に検出でき、うちの記録は検出されない側にいる**。"
      % (f(S2["controls"]["WIND-A"]["height"]["phase_aligned"]["gf_S"]
           / max(S2["scenes"]["WIND-A"]["spaces"]["height"]["gf_S"]["median"],
                 1e-12), 2),
         f(S2["controls"]["WIND-B"]["height"]["phase_aligned"]["gf_S"]
           / max(S2["scenes"]["WIND-B"]["spaces"]["height"]["gf_S"]["median"],
                 1e-12), 2),
         f(S2["controls"]["WIND-C"]["height"]["phase_aligned"]["gf_S"]
           / max(S2["scenes"]["WIND-C"]["spaces"]["height"]["gf_S"]["median"],
                 1e-12), 2),
         f(S2["controls"]["WIND-A"]["height"]["ensemble"]["gf_S_p95"], 3)))
    A("")
    A("注意: WIND-C の位相集中対照の絶対値は %s で、偶然 F&M の 0.95 とほぼ同じ"
      "数字になる。**これは無関係**である。WIND-C は窓が 1 波長しかない場なので、"
      "位相を揃えると SIWEH が「1 個の山」になり std/mean が 0.95 になるだけで"
      "あって、F&M の（多数の群れを持つ記録での）位相集中車線とは別物。"
      "倍数（%s 倍）で読むこと。"
      % (f(S2["controls"]["WIND-C"]["height"]["phase_aligned"]["gf_S"], 3),
         f(S2["controls"]["WIND-C"]["height"]["phase_aligned"]["gf_S"]
           / max(S2["scenes"]["WIND-C"]["spaces"]["height"]["gf_S"]["median"],
                 1e-12), 2)))
    A("")
    A("**Funke & Mansard の 0.95 という絶対値は 15 周期の記録では再現しない**。"
      "8 s に焦点が 1 回だけ入る記録では SIWEH が「ほぼ 0 ＋ 1 個の山」になり、"
      "std/mean は 1 を超える（WIND-A の対照で %s）。したがって 0.95 は目標としても"
      "比較尺度としても使えない。**使えるのは、同じ場・同じ記録長・同じ推定量で"
      "作ったランダム位相帯だけ**である（§6）。"
      % f(S2["controls"]["WIND-A"]["height"]["phase_aligned"]["gf_S"], 3))
    A("")
    spA = S2["scenes"]["WIND-A"]["spectral"]
    alA = spA["aliased_share_above_temporal_nyquist"]
    A("**問い 3: 勾配／曲率空間で何が加わったか。** 3 空間の GF・ν は同じ場面でも"
      "一致しない。微分階数を上げると統計を担う波数が上がる（WIND-A の重心 k: "
      "高さ %s / 勾配 %s / 曲率 %s rad/m、波長 %s / %s / %s mm）。曲率は"
      "格子分解能で %s セル/波長しかなく、しかも 24 fps の時間 Nyquist "
      "（%s rad/s）より上の帯に曲率分散の %s が落ちる。したがって"
      "**曲率の GF は「24 fps で記録された曲率時系列の包絡統計」であり、"
      "時間連続の曲率場の群性ではない**。高さ空間だけはエイリアス率が %s で、"
      "群性の比較にそのまま使える。"
      % (f(spA.get("k_centroid_height"), 0), f(spA.get("k_centroid_slope"), 0),
         f(spA.get("k_centroid_curv"), 0),
         f((spA.get("lambda_height_m") or 0) * 1e3, 1),
         f((spA.get("lambda_slope_m") or 0) * 1e3, 1),
         f((spA.get("lambda_curvature_m") or 0) * 1e3, 1),
         f(spA.get("cells_curvature"), 2),
         f(spA.get("temporal_nyquist_rad_s"), 1), pct(alA["curv"], 2),
         pct(alA["height"], 2)))
    A("")
    A("**問い 4: 測れたもの／測れなかったもの。**")
    A("")
    A("測れた: 4 場面 × 3 空間の GF_S・GF_A・ν（中央値＋四分位範囲、256 プローブ）、"
      "場面固有のランダム位相帯（＝小窓・短記録の標本ばらつきの実測）、位相集中"
      "との差、平滑化／プローブ数／時間率／境界処理への感度、各空間を担う波数と"
      "エイリアス帯の割合、実測振幅とスペクトル読みの差。")
    A("")
    A("測れなかった: (i) **8 s より長い記録での GF**。ω が 2π/8 s の整数倍に"
      "量子化されているので場は厳密に 8 s 周期で、記録を延ばしても同じ周期の"
      "繰り返しになる（新しい統計は増えない）。(ii) **時間連続の GF**。24 fps "
      "記録からは復元できず、4 倍のフレーム数が必要（本ジョブの範囲外、"
      "エイリアス率だけ定量化した）。(iii) **STREAM の群性**。2 s 記録は継ぎ目"
      "周期の %s 周期分しかなく、包絡の独立標本が 1〜2 個（§7）。(iv) **うねり"
      "場面**を測っていないので「風波が最大・古いうねりが最小」という順序は"
      "本ジョブでは検証できない（§6）。"
      % f(S2["scenes"]["STREAM"]["window"]["record_s"]
          / S1["instances"]["STREAM"]["T_seam_s"], 2))
    A("")

    # ---- 1 ---------------------------------------------------------------
    A("## 1. 再現ゲート（STAGE 1）— 測る前に「その場」であることを証明する")
    A("")
    A("| 場面 | 量 | 記録値 | 再生成値 | 相対誤差 | 判定 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for sid, gg in S1["gates"].items():
        for r in gg["rows"]:
            if isinstance(r.get("rebuilt"), (int, float)):
                A("| %s | `%s` | %s | %s | %.2e | %s |"
                  % (sid, r["quantity"], f(r["recorded"], 10),
                     f(r["rebuilt"], 10), r["rel_err"],
                     "ok" if r["pass"] else "NG"))
            else:
                A("| %s | `%s` | %s | %s | %.2e | %s |"
                  % (sid, r["quantity"], r["recorded"], r["rebuilt"],
                     r["rel_err"], "ok" if r["pass"] else "NG"))
    A("")
    A("WIND-C の `max_abs_diff` は **%s**（記録配列 `flow-0000.npy` との絶対差、"
      "float32 で完全一致）。STREAM の frame 24 緑歪み平均は **%.6f px** で記録"
      "（10.067451730266395）と一致。**再生成した場は記録された場そのもの**である。"
      % (f([r for r in S1["gates"]["WIND-C"]["rows"]
            if "max_abs_diff" in r][0]["max_abs_diff"], 1),
         S1["gates"]["STREAM"]["rows"][0]["rebuilt"]))
    A("")
    A("**生成スクリプトの来歴（実測 sha256）。**")
    A("")
    A("| ファイル | sha256（実測） | v19 記録 | v20 記録 | v4-1-8 記録 |")
    A("| --- | --- | --- | --- | --- |")
    rows = {r["file"]: r for r in numbers["stage0"]["provenance_rows"]}
    for fn in ("run_v19_prod.py", "run_v20_prod.py", "run_prod_v4.py",
               "state.py", "run_v4_1_8.py", "stream_surface3.py",
               "stream_surface2.py", "run_v15.py", "run_v15_spectrum.py"):
        r = rows[fn]
        def verdict(tab):
            v = r["recorded"].get(tab, {}).get("verdict", "—")
            if v == "DIFF" and fn == "run_v4_1_8.py":
                return "**自己編集**"
            return {"match": "一致", "DIFF": "**不一致**",
                    "not-listed": "未収載"}.get(v, v)
        A("| `%s` | `%s…` | %s | %s | %s |"
          % (fn, r["sha256_now"][:16],
             verdict("out/v19/stage0-py-hashes.json"),
             verdict("out/v20/stage0-py-hashes.json"),
             verdict("out/v4-1-8/stage0-py-hashes.json")))
    A("")
    n_diff = len({d["file"] for d in numbers["stage0"]["hash_diff"]})
    sh = numbers["stage0"].get("self_hashes", {})
    A("不一致は %d 本（%s）。" % (n_diff, ", ".join(
        "`%s`" % f for f in sorted({d["file"]
                                    for d in numbers["stage0"]["hash_diff"]}))))
    A("")
    A("* `run_v15.py` / `run_v15_spectrum.py` は v19/v20 の記録（2026-11-11/12）"
      "より後に更新されたため（v4-1-8 の記録 2026-11-13 とは一致する）。")
    A("* `run_v4_1_8.py` は**自分のジョブで自己編集された**ファイルで、stage-0 "
      "表の値は編集前のスナップショット。最終ハッシュは `numbers-v4-1-8.json` の "
      "`run_v4_1_8_sha256` に記録されており、現在のファイルと%s（`%s…`）。"
      % ("一致する" if sh.get("v4-1-8", {}).get("matches_now")
         else "一致しない",
         (sh.get("v4-1-8", {}).get("recorded_final_sha256") or "?")[:16]))
    A("")
    A("ブリーフの指示どおり **sha の字面ではなく行動ゲートで判定**し、"
      "WIND-A/B の 3 スカラーは相対誤差 0 で記録と完全一致、WIND-C は記録配列と"
      "完全一致、STREAM は歪み平均が一致した（＝いずれの更新も数値挙動に"
      "触れていない）。")
    A("")
    A("**測定量の生成器との一致（provider check, t = 8/24 s, プローブ点）。**")
    A("")
    A("| 場面 | h の rms 相対差 | gx | gy | 曲率: スペクトル/5点差分 の std 比 |")
    A("| --- | --- | --- | --- | --- |")
    for sid, pc in S1["provider_checks"].items():
        A("| %s | %.1e | %.1e | %.1e | %.4f |"
          % (sid, pc["h"]["rms_rel"], pc["gx"]["rms_rel"], pc["gy"]["rms_rel"],
             pc["lap_vs_fd"]["ratio_std"]))
    A("")
    A("h/gx/gy が 1e-7 台で一致するので、軽量生成は正規の `fields_at` と同じ場を"
      "見ている。曲率の比が 1 から離れるのは 5 点差分の離散化誤差（毛細管帯で "
      "k·dx が小でない）であって、場の不一致ではない。")
    A("")
    A("**ブリーフの表との食い違い（申告）。** ブリーフの WIND-C 行は WIND-B と"
      "同じ「U=3, F=100, elfouhaily-form」と書いているが、記録された v2 ループ"
      "（`out/prod/numbers-prod.json`, `pipeline: v4`）を生成したのは "
      "`run_prod_v4.py` → `state.WaterState`（**Phillips スペクトル**）で、"
      "recipe-1.json の実値は U = %s m/s、specL = %s m、σ_η 目標 = %s m。"
      "測定対象は「記録された場」なので、STAGE 1 は字面ではなくこの実際の生成器で"
      "再生成した。**WIND-C だけスペクトル型が違う**（elfouhaily-form/Arm D では"
      "ない）ことを明記する。"
      % (f(numbers["stage0"]["v2_recipe_env"].get("windSpeedMS"), 2),
         f(numbers["stage0"]["v2_recipe_env"].get("specL"), 8),
         f(numbers["stage0"]["v2_recipe_env"].get("waveHeightM"), 8)))
    A("")

    # ---- 2 ---------------------------------------------------------------
    A("## 2. 方法（数値で固定、事後調整なし）")
    A("")
    A("* **プローブ**: 状態格子のノード上に 16×16 = 256 点、stride N/16、縁を除く。"
      "輸送段が読むのと同じサンプルで、補間なし。")
    A("* **包絡**: 各プローブ時系列 s(t) の FFT ヒルベルト（解析信号）"
      "A(t) = |s + iH[s]|。境界は風波 3 場面が**周期**（ω が 2π/8 s の整数倍に"
      "量子化され、192 サンプルは厳密に 1 周期）、STREAM が**反射パディング**"
      "（移流 ω = kx·Us は量子化されていない）。")
    A("* **SIWEH**: A²(t) を **1 ピーク周期**の円形箱型窓で平滑化。"
      "GF_S = std(SIWEH)/mean(SIWEH)（母標準偏差）。")
    A("* **GF_A**: std(A)/mean(A)。同じ A の包絡比で、Rayleigh 包絡なら ≈0.42、"
      "SIWEH 型は ≈1.0。**流儀の違いが 0.42 と 1.0 の差**なので両方併記する"
      "（文献値と比べるときに必要な唯一の情報）。")
    A("* **ν**: プローブ記録の片側周波数モーメントから ν = √(m0·m2/m1² − 1)。"
      "勾配の ν は**風軸成分 gx**（gy は併記）。ω が量子化されているので"
      "理論値との差は §9 で開示。")
    A("* **時間**: 各場面の元ランと同じ 24 fps / 192 フレーム（風波）・48 フレーム"
      "（stream）。**間引きはしていない**。12/8 fps は感度のための参考車線。")
    A("* **集計**: 256 プローブの中央値と四分位範囲。単一値では語らない。")
    A("* **対照（場面ではない）**: 同一分散スペクトル・同一記録長のランダム位相"
      "代理 400 本（GF の標本分布）と、位相を t=0 に揃えた 1 本（位相集中の"
      "検出器）。どちらも 1 次元の記録で、新しい場面の生成・合成・動画化は"
      "行っていない。")
    A("* **推定量の検査**: 整数周期の純正弦 → GF_S = 0、2 音のビート → 0.69、"
      "狭帯域ガウス（Rayleigh 包絡クラス）→ GF_A 0.47 / GF_S 0.63。"
      "文献値と比べる前に、この 3 点で推定量の向きと流儀を固定した。")
    A("")

    # ---- 3 ---------------------------------------------------------------
    A("## 3. GF / ν 全表（場面 × 空間）")
    A("")
    A("GF_S・ν は 256 プローブの中央値 [q25–q75]。GF_A は同じ A の包絡比。"
      "「帯」は同じ場面・同じ空間のランダム位相代理（400 本）の 5–95%。"
      "文献の 0.57 と比べるのは GF_A の列（§6）。")
    A("")
    A("| 場面 | 空間 | GF_S | GF_A | ν | ランダム位相帯 5–95% | 内/外 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for row in G["G-DIAG"]["rows"]:
        sid = row["scene"]
        tag = {"height": "高さ η", "slope": "勾配 |∇η|",
               "curvature": "曲率 ∇²η"}[row["space"]]
        e = S2["controls"][sid][row["space"]]["ensemble"]
        A("| %s | %s | %s [%s–%s] | %s | %s [%s–%s] | [%s–%s] | %s |"
          % (sid, tag, f(row["gf_S"]["median"]), f(row["gf_S"]["q25"], 3),
             f(row["gf_S"]["q75"], 3), f(row["gf_A"]["median"]),
             f(row["nu"]["median"]), f(row["nu"]["q25"], 3),
             f(row["nu"]["q75"], 3), f(e["gf_S_p5"], 3), f(e["gf_S_p95"], 3),
             "内" if e["gf_S_p5"] <= row["gf_S"]["median"] <= e["gf_S_p95"]
             else "外"))
    A("")
    A("成分別の勾配（参考）:")
    A("")
    A("| 場面 | GF_S(gx) | GF_S(gy) | ν(gx) | ν(gy) |")
    A("| --- | --- | --- | --- | --- |")
    for sc in SCENES:
        sp = S2["scenes"][sc["id"]]["spaces"]["slope_components"]
        A("| %s | %s | %s | %s | %s |"
          % (sc["id"], f(sp["gx"]["gf_S"]["median"]),
             f(sp["gy"]["gf_S"]["median"]), f(sp["gx"]["nu"]["median"]),
             f(sp["gy"]["nu"]["median"])))
    A("")
    A("高さの歪度・尖度（束縛高調波の指紋）:")
    A("")
    A("| 場面 | skew 中央値 [q25–q75] | kurt 中央値 [q25–q75] | ガウス期待 |")
    A("| --- | --- | --- | --- |")
    for sc in SCENES:
        hb = S2["scenes"][sc["id"]]["spaces"]["height"]
        A("| %s | %s [%s–%s] | %s [%s–%s] | 0 / 3 |"
          % (sc["id"], f(hb["skew"]["median"], 3), f(hb["skew"]["q25"], 3),
             f(hb["skew"]["q75"], 3), f(hb["kurt"]["median"], 3),
             f(hb["kurt"]["q25"], 3), f(hb["kurt"]["q75"], 3)))
    A("")
    hA = S2["scenes"]["WIND-A"]["spaces"]["height"]
    hB = S2["scenes"]["WIND-B"]["spaces"]["height"]
    hC = S2["scenes"]["WIND-C"]["spaces"]["height"]
    A("**読み方**: 歪度はすべて 0 近傍（|中央値| < 0.02、IQR は ±0.1 以内）で、"
      "**束縛高調波は存在しない**（線形合成なので当然）。尖度は 3 より低い"
      "（WIND-A %s / WIND-B %s / WIND-C %s）。これは束縛波の欠如ではなく、"
      "8 s 記録に含まれる**有効モード数が少ない**ことによる sub-Gaussian 性で、"
      "M 個のランダム位相正弦の和なら kurt ≈ 3 − 3/M。逆算すると M_eff ≈ %s / "
      "%s / %s。**「ka が数 %% だから束縛波が出るはず」という議論はこの生成器には"
      "当てはまらない**（出ないのが正しい）。"
      % (f(hA["kurt"]["median"], 3), f(hB["kurt"]["median"], 3),
         f(hC["kurt"]["median"], 3),
         f(3.0 / max(3.0 - hA["kurt"]["median"], 1e-6), 1),
         f(3.0 / max(3.0 - hB["kurt"]["median"], 1e-6), 1),
         f(3.0 / max(3.0 - hC["kurt"]["median"], 1e-6), 1)))
    A("")

    # ---- 4 ---------------------------------------------------------------
    A("## 4. 窓とカバレッジ（各行が何を見ているか）")
    A("")
    A("| 場面 | パッチ | 深さ | U / F | 格子 / dx | フレーム / 秒 | T_p | 記録内の周期数 | 平滑化 | プローブ |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for sc in SCENES:
        w = S2["scenes"][sc["id"]]["window"]
        A("| %s | %s m | %s m | %s / %s | %d² / %s mm | %d / %s s | %s s | %.2f | %d サンプル (%s s) | %d (%d×%d) |"
          % (sc["id"], f(w["patch_m"], 4), f(w["depth_m"], 2),
             f(w["U_ms"], 2), f(w["F_m"], 0) if w["F_m"] else "—",
             w["grid"], f(w["dx_mm"], 3), w["frames"], f(w["record_s"], 3),
             f(w["T_peak_s"], 4), w["periods_in_record"],
             w["smooth_window_samples"], f(w["smooth_window_s"], 3),
             w["probes"], NPROBE, NPROBE))
    A("")
    A("パッチ / ピーク波長: " + "、".join(
        "%s %s" % (sc["id"],
                   ("S/λp = %s" % f(S2["scenes"][sc["id"]]["window"]
                                    ["patch_m"]
                                    / S2["scenes"][sc["id"]]["window"]
                                    ["lambda_p_m"], 2))
                   if S2["scenes"][sc["id"]]["window"]["lambda_p_m"]
                   else "λp 未記録")
        for sc in SCENES) + "。")
    A("")
    A("**小窓の注意。** WIND-C のパッチ %s m は、その場のエネルギーが最も強い"
      "波長（Phillips のピークは k_p = √2/specL、specL = %s m）に対して"
      "**1 波長程度**しかない。したがって WIND-C の GF は「有限窓に切り取られた"
      "低波数モード数本の干渉」の統計であって、開けた海の群性ではない。"
      "WIND-B も窓の中のピーク波は 9 個程度。**窓を明記しない GF の比較は"
      "意味を持たない**、というのがこの表の役割。"
      % (f(S2["scenes"]["WIND-C"]["window"]["patch_m"], 4),
         f(numbers["stage0"]["v2_recipe_env"].get("specL"), 4)))
    A("")

    # ---- 5 ---------------------------------------------------------------
    A("## 5. 高さ → 勾配 → 曲率：何が加わったか")
    A("")
    A("| 場面 | k 重心 高さ / 勾配 / 曲率 [rad/m] | λ 高さ / 勾配 / 曲率 [mm] | 曲率の格子分解能 [cells/λ] | 24fps Nyquist 超の分散率 高さ / 勾配 / 曲率 |")
    A("| --- | --- | --- | --- | --- |")
    for sc in SCENES:
        sp = S2["scenes"][sc["id"]]["spectral"]
        al = sp.get("aliased_share_above_temporal_nyquist") or {}
        A("| %s | %s / %s / %s | %s / %s / %s | %s | %s / %s / %s |"
          % (sc["id"], f(sp.get("k_centroid_height"), 0),
             f(sp.get("k_centroid_slope"), 0), f(sp.get("k_centroid_curv"), 0),
             f((sp.get("lambda_height_m") or 0) * 1e3, 2),
             f((sp.get("lambda_slope_m") or 0) * 1e3, 2),
             f((sp.get("lambda_curvature_m") or 0) * 1e3, 2),
             f(sp.get("cells_curvature"), 2),
             pct(al.get("height"), 2), pct(al.get("slope"), 2),
             pct(al.get("curv"), 2)))
    A("")
    A("実測振幅とスペクトル読みの差（測定量がどの場を見ているかの検算）:")
    A("")
    A("| 場面 | rms η 実測 / スペクトル | rms slope 実測 / スペクトル | rms ∇²η 実測 |")
    A("| --- | --- | --- | --- |")
    for sc in SCENES:
        am = S2["scenes"][sc["id"]]["amplitudes"]
        A("| %s | %s / %s | %s / %s | %s |"
          % (sc["id"], f(am["rms_eta_measured"], 6),
             f(am.get("rms_eta_spectral"), 6), f(am["rms_slope_measured"], 6),
             f(am.get("rms_slope_spectral"), 6),
             f(am["rms_lap_measured"], 5)))
    A("")
    lat = S2["scenes"]["WIND-A"].get("lattice") or {}
    if lat:
        am = S2["scenes"]["WIND-A"]["amplitudes"]
        b1 = S2["scenes"]["WIND-B"]["amplitudes"]
        rA = am["rms_eta_measured"] / max(am["rms_eta_spectral"], 1e-30)
        rB = b1["rms_eta_measured"] / max(b1["rms_eta_spectral"], 1e-30)
        rF = lat["rms_full_grid"] / max(am["rms_eta_spectral"], 1e-30)
        A("**実測 rms がスペクトル読みと数 %% ずれる理由（2 つ、どちらも開示）。** "
          "(a) スペクトル読み `m0` はモード集合の分散（アンサンブル期待値）で、"
          "1 回の実現場には `mean(z²)` の干渉項が残る。全面格子の rms ÷ 読みは "
          "WIND-A %s（WIND-B では −4.2%% の干渉項でこの差が説明できる）。"
          "(b) 16×16 プローブ格子は、パッチ内の波長数が 9〜18 個しかない場の疎な"
          "標本なので、プローブ rms ÷ 読み = WIND-A %s / WIND-B %s、プローブ rms ÷ "
          "全面 rms = %s（WIND-A）。格子を 16〜80 セルずらすと %s〜%s 倍に動く。"
          "**GF は規格化統計なのでこの振幅差は主項では消えるが、同じ有限標本効果が "
          "GF のばらつき（IQR と帯）の源**である。"
          % (f(rF, 3), f(rA, 3), f(rB, 3),
             f(lat["rms_probe_over_full"], 3),
             f(min(lat["rms_shifted_lattices"].values())
               / max(lat["rms_full_grid"], 1e-30), 3),
             f(max(lat["rms_shifted_lattices"].values())
               / max(lat["rms_full_grid"], 1e-30), 3)))
        A("")
    seam = S2["scenes"]["WIND-A"].get("seam")
    if seam:
        A("**毛細管タイルの折り返し（WIND-A）。** 高波数カスケードは 2×2 ブロック"
          "平均してタイル複製されるので、合成場に複製の継ぎ目が入る。継ぎ目での "
          "h の段差は平均 %s m、通常の隣接差は %s m（比 %s）。プローブは"
          "継ぎ目（index %d）を避けて配置した。曲率空間はこの折り返しの影響を"
          "最も強く受ける（ブロック平均は高波数ほど減衰させる）。"
          % (f(seam["mean_jump_at_tile_seam_m"], 7),
             f(seam["mean_neighbour_gradient_m"], 7), f(seam["ratio"], 2),
             seam["seam_index"]))
        A("")
    A("読み方: 高さ空間の GF は「うねりの群れ」を、勾配空間は「傾きの包絡」を、"
      "曲率空間は「曲がりの包絡」を見る。微分階数を上げると統計を担う波数が"
      "上がるため、同じ場面でも GF・ν の値は動く。**曲率の数値は 24 fps という"
      "記録の時間分解能と毛細管タイルの折り返しの影響を強く受ける**ので、"
      "コースティクスを決める曲率の群性の直接の物差しではなく、上限側の指標と"
      "して読む。")
    A("")

    # ---- 6 ---------------------------------------------------------------
    A("## 6. 文献との対比")
    A("")
    nbA = S2["controls"]["WIND-A"]["narrowband_gaussian"]
    A("* **Funke & Mansard**: 自然な合成海面 GF ≈ 0.57、位相集中 0.95。"
      "**絶対値の対応付けは本ジョブでは確定できない**（彼らの SIWEH 窓と記録長が"
      "不明で、こちらの窓は 1 T_p 固定）。確定できるのは 3 点: (i) 自然な海面が"
      "属する狭帯域ガウス記録（Rayleigh 包絡クラス）を同じ推定量・同じ記録長で"
      "測ると GF_A = %s / GF_S = %s、(ii) うちの高さの GF_A は WIND-A %s / "
      "WIND-B %s / WIND-C %s で、この Rayleigh クラスの値と数 %% 以内、"
      "(iii) 0.57 はうちの GF_A のランダム位相帯 [%s–%s] / [%s–%s] / [%s–%s] の"
      "内側。したがって**「自然な群性の範囲」という意味では、うちの 3 風波場面は"
      "文献の自然な海面と同じクラス**にいる。0.57 との厳密な一致は主張しない。"
      % (f(nbA["gf_A_median"], 3), f(nbA["gf_S_median"], 3),
         f(S2["scenes"]["WIND-A"]["spaces"]["height"]["gf_A"]["median"], 3),
         f(S2["scenes"]["WIND-B"]["spaces"]["height"]["gf_A"]["median"], 3),
         f(S2["scenes"]["WIND-C"]["spaces"]["height"]["gf_A"]["median"], 3),
         f(S2["controls"]["WIND-A"]["height"]["ensemble"]["gf_A_p5"], 3),
         f(S2["controls"]["WIND-A"]["height"]["ensemble"]["gf_A_p95"], 3),
         f(S2["controls"]["WIND-B"]["height"]["ensemble"]["gf_A_p5"], 3),
         f(S2["controls"]["WIND-B"]["height"]["ensemble"]["gf_A_p95"], 3),
         f(S2["controls"]["WIND-C"]["height"]["ensemble"]["gf_A_p5"], 3),
         f(S2["controls"]["WIND-C"]["height"]["ensemble"]["gf_A_p95"], 3)))
    A("* **位相集中 0.95 は目標でも尺度でもない。** 位相を揃えた対照は GF_S で "
      "%s / %s / %s と、ランダム位相帯の外へ大きく出る（帯の上限は WIND-A %s）。"
      "WIND-C の %s は F&M の 0.95 と数字が近いだけで無関係（§0 の注意）。"
      "この対照は「位相集中が起きていない」ことの検出器としてだけ使う。"
      "0.95 という数字は記録長と窓に依存するので、絶対値の比較には使わない。"
      % (f(S2["controls"]["WIND-A"]["height"]["phase_aligned"]["gf_S"], 3),
         f(S2["controls"]["WIND-B"]["height"]["phase_aligned"]["gf_S"], 3),
         f(S2["controls"]["WIND-C"]["height"]["phase_aligned"]["gf_S"], 3),
         f(S2["controls"]["WIND-A"]["height"]["ensemble"]["gf_S_p95"], 3),
         f(S2["controls"]["WIND-C"]["height"]["phase_aligned"]["gf_S"], 3)))
    A("* **SIWEH 型 GF_S の絶対値は 0.57 と比べない。** GF_S は同じ場で "
      "%s / %s / %s と GF_A より大きく、帯も [0.52–1.10] のように広い。"
      "窓を 2 T_p にすると 0.60、0.5 T_p にすると 0.80 に動く（§8）。"
      "**絶対値は推定量の流儀と窓の関数**であり、場の性質ではない。"
      % (f(S2["scenes"]["WIND-A"]["spaces"]["height"]["gf_S"]["median"], 3),
         f(S2["scenes"]["WIND-B"]["spaces"]["height"]["gf_S"]["median"], 3),
         f(S2["scenes"]["WIND-C"]["spaces"]["height"]["gf_S"]["median"], 3)))
    A("* **海洋解析（19 記録）**: 風波が GF の中央値で最大、古いうねりが最小、"
      "という順序が報告されている。うちの 3 風波場面はすべて風波（うねり成分"
      "なし）なので、この順序の中では上側の群に入る。**順序そのものは本ジョブ"
      "では検証できない**（うねり場面を測っていない）。")
    A("* **Chukchi 海**: 群性は急峻度・帯域幅・BFI と相関。うちの 4 場面の"
      "**高さ空間**では ν（帯域幅）と GF_S が単調な逆相関になる: ν = %s → GF %s、"
      "%s → %s、%s → %s、%s → %s。ν→0（狭帯域＝長い群れ）で GF が大きく、"
      "ν≳0.5（広帯域で溶ける）で小さい、という読みと整合する。"
      "ただし空間をまたぐと単調性は崩れる（同じ場面で勾配は ν が大きく GF が"
      "小さいが、WIND-C の曲率は ν 0.28 で GF 0.33 — 高さの 0.51 より小さい）。"
      "**ν は高さの帯域幅の物差しであって、勾配・曲率の群性の予測子ではない。**"
      % (f(S2["scenes"]["WIND-B"]["spaces"]["height"]["nu"]["median"], 3),
         f(S2["scenes"]["WIND-B"]["spaces"]["height"]["gf_S"]["median"], 3),
         f(S2["scenes"]["WIND-A"]["spaces"]["height"]["nu"]["median"], 3),
         f(S2["scenes"]["WIND-A"]["spaces"]["height"]["gf_S"]["median"], 3),
         f(S2["scenes"]["WIND-C"]["spaces"]["height"]["nu"]["median"], 3),
         f(S2["scenes"]["WIND-C"]["spaces"]["height"]["gf_S"]["median"], 3),
         f(S2["scenes"]["STREAM"]["spaces"]["height"]["nu"]["median"], 3),
         f(S2["scenes"]["STREAM"]["spaces"]["height"]["gf_S"]["median"], 3)))
    A("")
    A("**参考値（ゲートではない）: BFI と平均群長。**")
    A("")
    A("| 場面 | Hs [m] | λp（記録から逆算）[m] | λp（記録値）[m] | k_p·a | Δω/ω_p | BFI | 1/ν | λp/ν [m] |")
    A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for sc in SCENES:
        if sc["id"] == "STREAM":
            A("| STREAM | — | — | — | — | — | — | — | — （強制帯状パターン。BFI は定義しない） |")
            continue
        m = B[sc["id"]]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (sc["id"], f(m["Hs"], 5), f(m["lambda_p"], 4),
             f(m["lambda_p_recorded"], 4), f(m["ka"], 3),
             f(m["domega_over_omega"], 3), f(m["BFI"], 2),
             f(m["group_len_1_over_nu"], 2), f(m["group_len_m"], 3)))
    A("")
    A("BFI の定義は文献で O(1) 係数が動くので、`k_p·a`（a = Hs/2）と "
      "`Δω/ω_p`（= ν で近似）を両方出して、読者が任意の流儀で組み直せるように"
      "した。λp は記録の ω_peak を生成器の分散関係で逆算したもの（記録値と"
      "並記）。")
    A("")

    # ---- 7 ---------------------------------------------------------------
    A("## 7. STREAM は別枠")
    A("")
    A("stream 面はランダム位相の風波ではなく、**強制された帯状パターン**"
      "（定在リング k0 = %s rad/m ＋ 底面応答 ＋ 乱流成分）である。その「群れ」は"
      "包絡群ではなく定常パターンの構造なので、風波 3 場面と同一尺度で並べない。"
      % f(S1["instances"]["STREAM"].get("k0"), 3))
    A("")
    A("* 記録は %d フレーム / %s s。継ぎ目周期 T_seam = L/Us = %s s に対して "
      "**%s 周期**しかない（v4-1-8 の記録どおり）。"
      % (S2["scenes"]["STREAM"]["window"]["frames"],
         f(S2["scenes"]["STREAM"]["window"]["record_s"], 3),
         f(S1["instances"]["STREAM"]["T_seam_s"], 4),
         f(S2["scenes"]["STREAM"]["window"]["record_s"]
           / S1["instances"]["STREAM"]["T_seam_s"], 3)))
    A("* GF_S（高さ）= %s、ν = %s。**2 s 記録では包絡の独立標本が 1〜2 個しか"
      "取れない**ので、この数値は「stream 場の包絡統計の 1 標本」であって群性の"
      "推定値ではない。さらに T_p = %s s は 2 s 記録の周波数分解能（0.5 Hz 刻み）"
      "そのものなので、ピーク周期には ±0.25 s 程度の不確かさが乗る。"
      % (f(S2["scenes"]["STREAM"]["spaces"]["height"]["gf_S"]["median"], 4),
         f(S2["scenes"]["STREAM"]["spaces"]["height"]["nu"]["median"], 4),
         f(S2["scenes"]["STREAM"]["window"]["T_peak_s"], 3)))
    A("* §3 の表には同じ推定量を当てた値を載せてある（列は共通、解釈だけ分ける）。")
    A("* **曲率は特に読めない**: stream の曲率を担う波数は k ≈ %s rad/m"
      "（λ ≈ %s mm）で、格子分解能は %s セル/波長しかない。つまり曲率場は"
      "格子の Nyquist 近傍が支配的で、その 90.9%% は 24 fps の時間 Nyquist 超。"
      "**stream の曲率 GF は数値として report するが、物理量としては解釈しない**。"
      % (f((S2["scenes"]["STREAM"]["spectral"].get("k_centroid_curv") or 0), 0),
         f((S2["scenes"]["STREAM"]["spectral"].get("lambda_curvature_m") or 0)
           * 1e3, 2),
         f(S2["scenes"]["STREAM"]["spectral"].get("cells_curvature"), 2)))
    A("")

    # ---- 8 ---------------------------------------------------------------
    A("## 8. 感度（WIND-A で示す。全場面は numbers JSON）")
    A("")
    sen = S2["sensitivity"]["WIND-A"]
    A("平滑化窓（GF_S 中央値）:")
    A("")
    A("| 窓 | 高さ box | 高さ gauss | 勾配 box | 曲率 box |")
    A("| --- | --- | --- | --- | --- |")
    for kk in ("0.5p", "1p", "2p"):
        r = sen["smoothing"]["lanes"]
        A("| %s T_p | %s | %s | %s | %s |"
          % (kk, f(r["height"]["box_" + kk]), f(r["height"]["gauss_" + kk]),
             f(r["slope"]["box_" + kk]), f(r["curvature"]["box_" + kk])))
    A("")
    A("プローブ数 / 時間率 / 境界:")
    A("")
    A("| 変更 | 高さ | 勾配 | 曲率 |")
    A("| --- | --- | --- | --- |")
    v0 = {t: S2["scenes"]["WIND-A"]["spaces"][t]["gf_S"]["median"]
          for t in SPACES}
    A("| プローブ 16×16（本番） | %s | %s | %s |"
      % (f(v0["height"]), f(v0["slope"]), f(v0["curvature"])))
    for gname in ("8", "4"):
        v = sen["probe_count"][gname]
        A("| プローブ %s×%s | %s | %s | %s |"
          % (gname, gname, f(v["height"]), f(v["slope"]),
             f(v["curvature"])))
    for rate, v in sen["rate"]["by_fps"].items():
        A("| %s fps%s | %s | %s | %s |"
          % (rate, "（本番）" if rate == str(FPS) else "（間引き参考）",
             f(v["height"]), f(v["slope"]), f(v["curvature"])))
    if "boundary" in sen:
        bb = sen["boundary"]["height"]
        A("| 境界: 周期 vs 反射（高さ） | %s / %s（差 %s） | — | — |"
          % (f(bb["periodic"]), f(bb["reflect"]), pct(bb["rel_diff"])))
    A("")
    A("**読み方**: GF は平滑化窓と時間率に強く依存する（WIND-A の高さで "
      "2 T_p にすると %s、0.5 T_p では %s）。**絶対値ではなく同一条件内の比較**"
      "として読むべき、というのがこの感度表の結論。プローブ数を 8×8 に落としても"
      "高さの GF は %s（相対差 %s）で、配置の影響は小さい。"
      % (f(sen["smoothing"]["lanes"]["height"]["box_2p"]),
         f(sen["smoothing"]["lanes"]["height"]["box_0.5p"]),
         f(sen["probe_count"]["8"]["height"]),
         pct(abs(sen["probe_count"]["8"]["height"] - v0["height"])
             / max(v0["height"], 1e-12))))
    A("")
    A("**推定量そのものの検査**（各場面の T_p で、記録長に整数周期の音を使う）: "
      "純正弦の GF_S は WIND-A %s / WIND-B %s / WIND-C %s（0 なら平滑化と包絡の"
      "実装は正しい）。ビートを作る 2 音（WIND-A は %s / %s Hz）では "
      "%s / %s / %s。順序が正しいことを確認した。"
      % (f(S2["controls"]["WIND-A"]["instrument_tests"]["pure_sine_gf_S"], 5),
         f(S2["controls"]["WIND-B"]["instrument_tests"]["pure_sine_gf_S"], 5),
         f(S2["controls"]["WIND-C"]["instrument_tests"]["pure_sine_gf_S"], 5),
         f(S2["controls"]["WIND-A"]["instrument_tests"]["two_tone_hz"][0], 3),
         f(S2["controls"]["WIND-A"]["instrument_tests"]["two_tone_hz"][1], 3),
         f(S2["controls"]["WIND-A"]["instrument_tests"]["two_tone_1.12_gf_S"], 3),
         f(S2["controls"]["WIND-B"]["instrument_tests"]["two_tone_1.12_gf_S"], 3),
         f(S2["controls"]["WIND-C"]["instrument_tests"]["two_tone_1.12_gf_S"], 3)))
    A("")

    # ---- 9 ---------------------------------------------------------------
    A("## 9. 弱点・開示")
    A("")
    A("1. **8 s 記録は 1 周期**。ω が 2π/8 s の整数倍に量子化されているので場は"
      "厳密に 8 s 周期。記録を延ばしても同じ周期の繰り返しで新しい統計は増えない。"
      "GF の不確かさはプローブ間のばらつき（§3 の IQR）と代理アンサンブル"
      "（§6 の帯）でしか評価できない。")
    A("2. **小窓**。WIND-C はパッチ ≈ 1 波長、WIND-B は ≈ 9 波長。低波数モード数が"
      "少なく、GF の標本ばらつき自体が大きい（IQR と帯がその実測）。")
    A("3. **時間エイリアス**。24 fps の Nyquist は %s rad/s。毛細管帯の一部は"
      "これを超え、記録された時系列では折り返している（§5 の分散率）。高さでは"
      "小さいが曲率では大きい。**時間連続の GF は測っていない**（4 倍のフレーム"
      "数が必要）。"
      % f(S2["scenes"]["WIND-A"]["spectral"]["temporal_nyquist_rad_s"], 1))
    A("4. **毛細管タイルの折り返し**。WIND-A は高波数カスケードを 2×2 平均して"
      "タイル複製するため、合成場に継ぎ目と高波数減衰が入る（§5）。曲率空間は"
      "この影響を最も強く受ける。")
    A("5. **線形合成**。束縛高調波は入っていないので歪度は 0 近傍（§3）。尖度は"
      "有効モード数が少ないため 3 より低い（sub-Gaussian、§3）。")
    A("6. **GF の流儀依存**。GF_S と GF_A を併記した。同じ場で GF_A 0.47 / "
      "GF_S 0.72（WIND-A 高さ）と 1.5 倍違う。文献値と比べるときはどちらの流儀かを"
      "確認する必要がある。")
    A("7. **WIND-C の素性**。記録された v2 ループは Phillips スペクトル・"
      "U = %s m/s で、ブリーフの表の字面（U=3, elfouhaily）とは違う（§1）。"
      "測ったのは記録された場。"
      % f(numbers["stage0"]["v2_recipe_env"].get("windSpeedMS"), 2))
    A("8. **STREAM は 1 標本**（§7）。")
    A("9. **実現場の分散とスペクトル読みは一致しない**（§5）。`m0` はアンサンブル"
      "期待値で、1 回の実現場は mean(z²) 干渉項のぶん数 % ずれる。GF を"
      "「スペクトルの理論値」から計算してはいけない理由がこれ。")
    A("10. **ν の量子化**。ω の離散化幅は 0.785 rad/s で、低風速場面の帯域幅に"
      "対して粗い。ν は生成器のモード集合から直接計算した値でもあり、"
      "「波の群れの幅」としての解像度はこの刻みに制限される。")
    A("")

    # ---- 10 --------------------------------------------------------------
    A("## 10. 再現手順と成果物")
    A("")
    A("```")
    A("python run_v4_1_11.py")
    A("```")
    A("")
    A("| ファイル | 内容 |")
    A("| --- | --- |")
    A("| `REPORT-V4-1-11.md` | 本レポート |")
    A("| `numbers-v4-1-11.json` | 全数値（stage0/1/2、全行、プローブ別 GF・ν、"
      "感度、対照、BFI） |")
    A("| `checkpoint.json` | jobkit チェックポイント（段ごと、stage0/1/2/3） |")
    A("| `artifact-hashes.json` | 成果物の最終 sha256（自己参照を避けて別ファイル） |")
    A("| `probes-*.npz` | プローブ時系列 η / gx / gy / ∇²η（256 × フレーム、"
      "float32）と probe index・x[m] |")
    A("")
    A("既存 .py は無改造（stage0 の sha256 と stage3 再計算の sha256 が、参照する "
      "%d 本すべてで一致。`G-OUT.existing_py_untouched = %s`）。**フル η キューブは"
      "保存していない**（保存したのはプローブ時系列・スペクトルモーメント・集計表"
      "のみ）。"
      % (len(CODE_FILES), G["G-OUT"]["existing_py_untouched"]))
    A("")
    A("全数値の所在: `stage0`（ハッシュ・記録）、`stage1`（再現ゲート・"
      "provider check・場の諸元）、`stage2.scenes.*`（窓・空間別 GF/ν・"
      "プローブ別値・振幅・スペクトル重心）、`stage2.sensitivity.*`、"
      "`stage2.controls.*`（代理アンサンブル・位相集中・推定量検査）、`bfi`、"
      "`gates`。成果物の sha256 は `artifact-hashes.json`。")
    A("")
    body = "\n".join(L) + "\n"
    atomic_write_bytes(REPORT, body.encode("utf-8"))
    return REPORT


def stage3(s0, s1, s2):
    t0 = time.time()
    log("=" * 78)
    log("[stage 3] numbers + report")
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers(s0, s1, s2)

    # G-OUT is decided on the real state of the tree, BEFORE anything is
    # written: the .py files are not touched by writing the report, so the
    # verdict is known in advance and the report can be written once, final.
    after = py_set()
    untouched = all(after.get(f) == s0["code"][f]["sha256"]
                    for f in CODE_FILES if s0["code"][f]["exists"])
    files = {"probes": sorted([x for x in os.listdir(OUT)
                               if x.startswith("probes-")])}
    numbers["gates"]["G-OUT"] = {
        "pass": bool(untouched),
        "existing_py_untouched": bool(untouched),
        "files": files,
        "new_files_only": ["run_v4_1_11.py"],
        "note": ("self-referential hashes are avoided: the sha256 of the FINAL "
                 "REPORT-V4-1-11.md / numbers-v4-1-11.json / checkpoint.json is "
                 "recorded in checkpoint.json['artifact_hashes'], written after "
                 "the last write of every artefact")}
    write_json(NUMBERS, numbers)
    write_report(numbers)
    rep_sha = sha256_file(REPORT)
    numbers["gates"]["G-OUT"]["files"] = {
        "REPORT-V4-1-11.md": {"bytes": os.path.getsize(REPORT),
                              "sha256": rep_sha},
        "numbers-v4-1-11.json": {"bytes": None,
                                 "note": "own size/hash in artifact-hashes.json"},
        "checkpoint.json": {"bytes": None,
                            "note": "own size/hash in artifact-hashes.json"},
        "probes": files["probes"]}
    numbers["gates"]["G-OUT"]["pass"] = bool(
        untouched and os.path.getsize(REPORT) > 0
        and os.path.getsize(NUMBERS) > 0)
    write_json(NUMBERS, numbers)
    # final content is on disk now; the report does not depend on the byte
    # counts above, so it stays valid.  Record the final hashes in the
    # checkpoint (written last, so nothing goes stale).
    art = {"REPORT-V4-1-11.md": {"bytes": os.path.getsize(REPORT),
                                 "sha256": sha256_file(REPORT)},
           "numbers-v4-1-11.json": {"bytes": os.path.getsize(NUMBERS),
                                    "sha256": sha256_file(NUMBERS)},
           "probes": files["probes"]}
    cp_extra("artifact_hashes", art)
    cp_extra("gates", numbers["gates"])
    cp.unit_done("stage3-report")
    art["checkpoint.json"] = {"bytes": os.path.getsize(CP_PATH),
                              "sha256": sha256_file(CP_PATH),
                              "note": "recorded after this file's last write; "
                                      "the hash of checkpoint.json itself is "
                                      "not stored inside it"}
    # a small immutable sidecar makes every artefact verifiable without
    # re-running the job (and without a hash that contains its own file)
    write_json(os.path.join(OUT, "artifact-hashes.json"), art)
    log("  stage 3 wall %.1fs" % (time.time() - t0))
    return numbers["gates"]


def main():
    ap = argparse.ArgumentParser()
    for s in ("stage0", "stage1", "stage2", "stage3"):
        ap.add_argument("--" + s, action="store_true")
    args = ap.parse_args()
    only = any((args.stage0, args.stage1, args.stage2, args.stage3))

    t0 = time.time()
    s0 = s1 = s2 = None
    if not only or args.stage0:
        s0 = stage0()
    if not only or args.stage1:
        s1 = stage1()
    if not only or args.stage2:
        if s1 is None:
            s1 = read_json(CP_PATH)["stage1"]
        s2 = stage2(s1)
    if not only or args.stage3:
        if s0 is None:
            s0 = read_json(CP_PATH)["stage0"]
        if s1 is None:
            s1 = read_json(CP_PATH)["stage1"]
        if s2 is None:
            s2 = read_json(CP_PATH)["stage2"]
        gates = stage3(s0, s1, s2)
        log("  G-REPRO %s | G-DIAG %s | G-OUT %s"
            % tuple("PASS" if gates[k]["pass"] else "FAIL"
                    for k in ("G-REPRO", "G-DIAG", "G-OUT")))
        log("  report  %s" % REPORT)
        log("  numbers %s" % NUMBERS)
        log("total wall %.1fs" % (time.time() - t0))
        if not all(gates[k]["pass"] for k in ("G-REPRO", "G-DIAG", "G-OUT")):
            raise SystemExit(2)
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
