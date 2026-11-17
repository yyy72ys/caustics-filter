"""run_v4_1_15_d15 -- BRIEF-V4-1-15-D15: ONE still, STAGE B frame-0024 with
the water depth changed from 4.0 m to 1.5 m and nothing else changed.

NEW FILE ONLY: this driver.  `run_v4_1_15.py` (the v4-1.15 STAGE B driver) is
imported UNCHANGED and its `SCN["sub"]["depth"]` literal is the single number
touched -- and only in memory, for the duration of this process.  No existing
.py is edited, and nothing is written into out/v4-1-15/ (that directory belongs
to the still-running v4-1.15 job; it is READ only).

WHAT IS RENDERED
    stage      "sub"      (preset 5: Elfouhaily, U = 6 m/s, F = 100 m,
                           clear medium, FF clear phase n = 1.098 / mu = 3.37,
                           eta = 1)
    window     4.8 x 3.5217 m (frame WIDTH = patch, magnification 1)
    patch      4.8 m (patch_for rule -> 2 doublings of 1.2 m; same at both
               depths: lambda_p = 0.4317 m at U = 6, F = 100)
    grid       512^2 (STAGE B's own grid, cell 9.375 mm)
    seed       3255045276382562994 (= state.derive_seed(recipe, 1))
    field time 1.000 s = frame 24 at 24 fps
    composite  A = 2.6, Blow sigma = 3.0 px, E = 2.2 (frozen form)
    deposit    transport_mesh.mesh_C verbatim, BOX_BULK_MAX raised to 1<<20
               (the documented STAGE B performance switch, same value the
               STAGE B frame-24 still was produced with)
    delta      depth 4.0 m -> 1.5 m

WHY A CONTROL RUN COMES FIRST
    "Identical conditions" is not taken on trust.  Before the 1.5 m still is
    produced, the driver re-renders STAGE B's own frame 24 at depth 4.0 m
    through exactly the same call chain and compares the encoded PNG byte for
    byte with the still recorded by the v4-1.15 job
    (out/v4-1-15/frame-0024-sub.png).  If the hashes differ, the setup is NOT
    identical: the driver then stops, writes the reason into the record and
    produces no still (fail-closed, no fabrication).

ONE HONEST CAVEAT, MEASURED (not assumed)
    depth is a physical parameter of the wave field as well as of the optics.
    The seed and the complex draws are the same, but the mode-density
    amplitude (depth-limited Hs) and omega(k, h) are not; therefore the field
    at t = 1.0 s at 1.5 m is not bit-identical to the 4.0 m field.  That
    difference IS the requested change, and it is quantified in the record
    (phase-draw identity, amplitude ratio, active-mode mask, omega shift,
    h/slope differences at t = 0 and t = 1).  Every parameter that is NOT the
    depth is checked against the v4-1.15 STAGE B record and reported as
    identical.

OUTPUTS
    out/v4-1-15-extra/frame-0024-sub-d150.png   736 x 540 RGB
    out/v4-1-15-extra/record-d150.json

Usage (from caustics-filter/offline):
    python run_v4_1_15_d15.py
    python run_v4_1_15_d15.py --no-control     # skip the depth-4.0 control
"""
import argparse
import datetime
import hashlib
import io
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

import run_v4_1_15 as R15            # the STAGE B driver, imported UNCHANGED
import transport_mesh as TM          # the deposit (existing file, untouched)
import run_v10 as R10                # recipe loader (read only)
from state import derive_seed        # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-15-D15.md"
OUT = os.path.join(HERE, "out", "v4-1-15-extra")
PNG = os.path.join(OUT, "frame-0024-sub-d150.png")
REC = os.path.join(OUT, "record-d150.json")

REF_DIR = os.path.join(HERE, "out", "v4-1-15")      # READ ONLY (other job)
REF_CP = os.path.join(REF_DIR, "checkpoint.json")
REF_STILL = os.path.join(REF_DIR, "frame-0024-sub.png")

STAGE = "sub"
FRAME = 24
D_DELTA = 1.5
D_BASE = 4.0
BOX_BULK_STAGE_B = 1 << 20

#: what STAGE B was declared to be (text check on the imported module + on the
#: v4-1.15 record).  Anything else means "the setup cannot be reproduced
#: identically" and the driver stops.
EXPECT_SCN = {
    "kind": "wind", "L": 4.8, "depth": 4.0,
    "N": 512, "half": 256, "refine": 1024,
    "U": 6.0, "F": 100.0, "seed": 3255045276382562994,
    "water": "clear", "label": "サブ枠 (v19 size)",
}
EXPECT_COMPOSITE = {"A": 2.6, "Blow_sigma_px": 3.0, "E": 2.2, "Q": 25}
EXPECT_PHASE = {"kind": "ff", "n": 1.098, "mu": 3.37}

#: existing modules whose bytes must not change while this job runs
GUARD_FILES = ("run_v4_1_15.py", "transport_mesh.py", "run_v16_water.py",
               "v12_water.py", "v11_state.py", "v10_state.py",
               "run_v15_state.py", "run_v15_spectrum.py", "det.py",
               "v5_composite.py", "image.py", "state.py")

T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jnum(x):
    """numpy scalar -> plain python (JSON safe)."""
    if isinstance(x, (np.floating, float)):
        return float(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    return x


def jsonable(x):
    """Recursively numpy -> plain python (arrays become lists)."""
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return jnum(x)
    return x


def read_json_retry(path, tries=30, wait=0.5):
    last = None
    for _k in range(tries):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f), None
        except (PermissionError, OSError, ValueError) as e:
            last = e
            time.sleep(wait)
    return None, last


def write_json(path, obj):
    data = json.dumps(obj, ensure_ascii=False, indent=2,
                      allow_nan=True).encode("utf-8")
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    for _k in range(40):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.5)
    raise SystemExit("could not write %s" % path)


def png_bytes(arr):
    """Exactly R15.save_png's encoding, without touching the disk."""
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# reference (READ ONLY: out/v4-1-15 belongs to the running job)
# ---------------------------------------------------------------------------
def read_reference():
    ref = {"dir": REF_DIR, "writes_into_ref_dir": False}
    cp, err = read_json_retry(REF_CP)
    if cp is None:
        ref["checkpoint"] = {"error": "%s: %s" % (type(err).__name__, err)}
    else:
        f24 = ((cp.get("frame24") or {}).get(STAGE) or {})
        sc = ((cp.get("scenes") or {}).get(STAGE) or {})
        ref["checkpoint"] = {
            "updated_at": cp.get("updated_at"),
            "last_unit": cp.get("last_unit"),
            "completed_units": cp.get("completed_units"),
            "recorded_t_s": f24.get("t_s"),
            "recorded_still": f24.get("still"),
            "recorded_still_sha256": f24.get("still_sha256"),
            "recorded_C_frame_stats": f24.get("C_frame_stats"),
            "recorded_energy": f24.get("energy"),
            "recorded_deposit_wall_s": f24.get("deposit_wall_s"),
            "recorded_frame24_wall_s": f24.get("wall_s"),
            "scene": {k: sc.get(k) for k in
                      ("patch_m", "grid", "cell_m", "window_m", "depth_m",
                       "water_preset", "water_label", "a_rgb", "b_rgb",
                       "phase", "optics", "composite", "env", "loop",
                       "deposit", "buffer_adapter", "magnification",
                       "photo", "photo_sha256")},
            "wind": sc.get("wind"),
            "recorded_peak": f24.get("peak"),
        }
    if os.path.isfile(REF_STILL):
        try:
            ref["still_file"] = {"path": REF_STILL,
                                 "sha256": sha256_file(REF_STILL),
                                 "bytes": os.path.getsize(REF_STILL)}
        except OSError as e:
            ref["still_file"] = {"error": "%s: %s" % (type(e).__name__, e)}
    else:
        ref["still_file"] = {"error": "missing", "path": REF_STILL}
    # the reference the control is compared against: the delivered still file
    # if readable, else the checkpoint's recorded hash.
    cand = []
    if ref.get("still_file", {}).get("sha256"):
        cand.append(("out/v4-1-15/frame-0024-sub.png", ref["still_file"]["sha256"]))
    if ref.get("checkpoint", {}).get("recorded_still_sha256"):
        cand.append(("checkpoint frame24.sub.still_sha256",
                     ref["checkpoint"]["recorded_still_sha256"]))
    ref["control_reference"] = {"sources": [{"name": n, "sha256": s}
                                            for n, s in cand],
                                "sha256": cand[0][1] if cand else None,
                                "sources_agree": (len({s for _, s in cand}) <= 1
                                                  if cand else None)}
    return ref


def check_setup(ref):
    """Every non-depth parameter of STAGE B, checked.  Fail => stop."""
    sc = dict(R15.SCN[STAGE])
    rows = {}
    ok = True

    def eq(name, got, want):
        nonlocal ok
        same = bool(got == want)
        ok = ok and same
        rows[name] = {"stage_b_expected": want, "imported_driver": got,
                      "identical": same}

    for k, v in EXPECT_SCN.items():
        eq("SCN.%s" % k, sc.get(k), v)
    eq("A_FROZEN", float(R15.A_FROZEN), EXPECT_COMPOSITE["A"])
    eq("B_SIGMA_FROZEN", float(R15.B_SIGMA_FROZEN),
       EXPECT_COMPOSITE["Blow_sigma_px"])
    eq("E_FROZEN", float(R15.E_FROZEN), EXPECT_COMPOSITE["E"])
    eq("Q_DET", int(R15.Q_DET), EXPECT_COMPOSITE["Q"])
    eq("WIDTH", int(R15.WIDTH), 736)
    eq("HEIGHT", int(R15.HEIGHT), 540)
    eq("FPS", float(R15.FPS), 24.0)
    eq("STILL_IDX", int(R15.STILL_IDX), FRAME)

    rec = (ref.get("checkpoint") or {}).get("scene") or {}
    # the v4-1.15 record must describe the same STAGE B (if it is readable)
    if rec:
        rows["record.patch_m"] = {"stage_b_expected": float(sc["L"]),
                                  "record": rec.get("patch_m"),
                                  "identical": bool(rec.get("patch_m")
                                                    == float(sc["L"]))}
        rows["record.grid"] = {"stage_b_expected": int(sc["N"]),
                               "record": rec.get("grid"),
                               "identical": bool(rec.get("grid")
                                                 == int(sc["N"]))}
        rows["record.depth_m"] = {"stage_b_expected": float(sc["depth"]),
                                  "record": rec.get("depth_m"),
                                  "identical": bool(rec.get("depth_m")
                                                    == float(sc["depth"]))}
        rows["record.composite"] = {
            "stage_b_expected": {"A": float(R15.A_FROZEN),
                                 "Blow_sigma_px": float(R15.B_SIGMA_FROZEN),
                                 "E": float(R15.E_FROZEN),
                                 "Q": int(R15.Q_DET)},
            "record": {k: rec["composite"].get(k) for k in
                       ("A", "Blow_sigma_px", "E", "Q")},
            "identical": bool(
                rec["composite"].get("A") == float(R15.A_FROZEN)
                and rec["composite"].get("Blow_sigma_px")
                == float(R15.B_SIGMA_FROZEN)
                and rec["composite"].get("E") == float(R15.E_FROZEN)
                and rec["composite"].get("Q") == int(R15.Q_DET))}
        ph = rec.get("phase") or {}
        rows["record.phase"] = {
            "stage_b_expected": EXPECT_PHASE,
            "record": {k: ph.get(k) for k in ("kind", "n", "mu")},
            "identical": bool(ph.get("kind") == EXPECT_PHASE["kind"]
                              and ph.get("n") == EXPECT_PHASE["n"]
                              and ph.get("mu") == EXPECT_PHASE["mu"])}
        for name in ("record.patch_m", "record.grid", "record.depth_m",
                     "record.composite", "record.phase"):
            ok = ok and bool(rows[name]["identical"])
        seedrec = (ref["checkpoint"].get("wind") or {}).get("seed_source", "")
        rows["record.seed"] = {
            "stage_b_expected": int(sc["seed"]), "record": seedrec,
            "identical": bool(int(sc["seed"]) in (3255045276382562994,))}
        ok = ok and rows["record.seed"]["identical"]
    else:
        rows["record"] = {"note": "v4-1.15 checkpoint unreadable; "
                                  "module-level check only"}
    return {"pass": bool(ok), "rows": rows,
            "note": ("the depth literal is the ONLY value this job changes; it "
                     "is changed in memory after this check")}


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def guard_hashes():
    return {n: sha256_file(os.path.join(HERE, n)) for n in GUARD_FILES}


def state_node(st):
    cas = {}
    for name in ("lo", "hi"):
        c = getattr(st, name, None)
        if c is None:
            continue
        cas[name] = {"tag": c.tag, "N": int(c.N), "tile_m": float(c.tile_m),
                     "n_active_modes": int(np.count_nonzero(
                         np.asarray(c.H0) != 0)),
                     "k_build": jnum(c.stats.get("k_build")),
                     "viscous_cutoff": jnum(c.stats.get("viscous_cutoff")),
                     "m0_realised": jnum(c.stats.get("m0_realised")),
                     "sigma2_realised": jnum(c.stats.get("sigma2_realised"))}
    return {"class": type(st).__name__, "patch_m": float(st.S),
            "N": int(st.N), "depth_m": float(st.depth),
            "patch_doublings": int(st.patch_doublings),
            "lambda_p_m": float(st.laws["lambda_p_m"]),
            "Hs_depth_limited_m": float(st.laws.get("Hs_depth_limited_m",
                                                    float("nan"))),
            "spectrum_after": {k: float(v) for k, v in
                               st.spectrum_after.items()
                               if isinstance(v, (int, float))},
            "omega_c_rad_s": float(st.ctx["omega_c"]),
            "alpha_m_used": float(st.ctx["alpha_m_value"]),
            "cascades": cas}


def env_node(env):
    return {"sun_altitude_deg": float(env["sun"]["altitudeDeg"]),
            "sunDir": [float(x) for x in env["sunDir"]],
            "clarity": float(env["clarity"]),
            "turbid": float(env["turbid"]),
            "depthM": float(env["depthM"]),
            "windSpeedMS": float(env.get("windSpeedMS", float("nan")))}


def optics_node(optics):
    kp = optics.get("kernel_profile")
    return {"kernel": optics.get("kernel"),
            "phase_kind": optics.get("phase_kind"),
            "phase_label": (optics.get("phase") or {}).get("label"),
            "phase_n": (optics.get("phase") or {}).get("n"),
            "phase_mu": (optics.get("phase") or {}).get("mu"),
            "w_ballistic": [float(x) for x in optics["w_ballistic"]],
            "w_halo": [float(x) for x in optics["w_halo"]],
            "sigma_px_ff": float(optics.get("sigma_px") or 0.0),
            "sigma_px_v12_hg": float(optics.get("sigma_px_v12_hg") or 0.0),
            "kernel_r_max_px_state_grid": (float(kp["r_max_px"])
                                           if kp else None),
            "partition": jsonable(optics.get("partition") or {})}


def render(depth, label, keep_state=True):
    """One frame through STAGE B's own call chain, at `depth`."""
    log("[%s] depth %.3f m: build_scene (U=6, F=100, seed=%d, N=%d)"
        % (label, depth, int(R15.SCN[STAGE]["seed"]),
           int(R15.SCN[STAGE]["N"])))
    R15.SCN[STAGE]["depth"] = float(depth)
    TM.BOX_BULK_MAX = BOX_BULK_STAGE_B
    t0 = time.time()
    st, env, optics, spec = R15.build_scene(STAGE)
    t_build = time.time() - t0
    log("[%s] state built %.1fs: patch %.4f m, N %d, cell %.4f mm, "
        "lambda_p %.4f m, Hs %.4f m"
        % (label, t_build, st.S, st.N, 1e3 * st.S / st.N,
           float(st.laws["lambda_p_m"]),
           float(st.spectrum_after["Hs_m"])))
    photo = R15.load_photo()
    t1 = time.time()
    r = R15.render_frame(STAGE, FRAME, st, env, optics, photo)
    t_render = time.time() - t1
    info = r["det_info"]
    log("[%s] frame %d: t=%.3f s  C mean %.6f std %.6f p99 %.4f  "
        "E rel_err %.3e  wall %.1fs (deposit+composite)"
        % (label, FRAME, r["t_s"], float(np.mean(r["C"])),
           float(np.std(r["C"])), float(np.percentile(r["C"], 99)),
           float(info["rel_err"]), t_render))
    out = {"label": label, "depth_m": float(depth), "state": st, "env": env,
           "optics": optics, "spec": spec, "photo": photo, "frame": r,
           "t_build_s": float(t_build), "t_render_s": float(t_render),
           "box_bulk_max": int(TM.BOX_BULK_MAX),
           "png_bytes": png_bytes(r["still"])}
    out["png_sha256"] = sha256_bytes(out["png_bytes"])
    out["still_shape"] = [int(x) for x in np.asarray(r["still"]).shape]
    out["still_dtype"] = str(np.asarray(r["still"]).dtype)
    log("[%s] png sha256 %s" % (label, out["png_sha256"]))
    if not keep_state:
        out.pop("state")
    else:
        out["fields"] = {}
        for t in (0.0, float(r["t_s"])):
            f = st.fields_at(t)
            out["fields"][repr(t)] = {
                "h": np.asarray(f["h"], np.float64).copy(),
                "gx": np.asarray(f["gx"], np.float64).copy(),
                "gy": np.asarray(f["gy"], np.float64).copy()}
        out["H0"] = {}
        for name in ("lo", "hi"):
            c = getattr(st, name, None)
            if c is None:
                continue
            out["H0"][name] = {"H0": np.asarray(c.H0).copy(),
                               "omega_q": np.asarray(c.omega_q,
                                                     np.float64).copy()}
    return out


# ---------------------------------------------------------------------------
# comparisons
# ---------------------------------------------------------------------------
def _rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def field_cmp(a, b, t):
    fa = a["fields"][repr(t)]
    fb = b["fields"][repr(t)]
    ha, hb = fa["h"], fb["h"]
    dh = hb - ha
    ra, rb = _rms(ha), _rms(hb)
    corr = float(np.sum(ha * hb) / max(ra * rb * ha.size, 1e-300))
    sa = np.sqrt(fa["gx"] ** 2 + fa["gy"] ** 2)
    sb = np.sqrt(fb["gx"] ** 2 + fb["gy"] ** 2)
    return {"t_s": float(t),
            "h_rms_base_m": ra, "h_rms_d150_m": rb,
            "h_rms_ratio": float(rb / ra),
            "max_abs_dh_m": float(np.max(np.abs(dh))),
            "rms_dh_m": _rms(dh),
            "rms_dh_over_h_rms": float(_rms(dh) / max(ra, 1e-300)),
            "h_pearson_corr": corr,
            "slope_rms_base": _rms(sa), "slope_rms_d150": _rms(sb),
            "slope_rms_ratio": float(_rms(sb) / max(_rms(sa), 1e-300))}


def draw_cmp(a, b):
    """Are the complex draws the SAME draw (seed + phase)?"""
    out = {"note": ("H0 = scale*(a+ib)/sqrt2 with the same seeded (a+ib) for "
                    "both depths; the scale (mode density, i.e. the "
                    "depth-limited Hs) and the active mask depend on h")}
    for name in ("lo", "hi"):
        if name not in a["H0"] or name not in b["H0"]:
            continue
        Ha = a["H0"][name]["H0"]
        Hb = b["H0"][name]["H0"]
        oa = a["H0"][name]["omega_q"]
        ob = b["H0"][name]["omega_q"]
        ma, mb = Ha != 0, Hb != 0
        m = ma & mb
        node = {"n_active_base": int(ma.sum()),
                "n_active_d150": int(mb.sum()),
                "n_common": int(m.sum()),
                "n_only_base": int((ma & ~mb).sum()),
                "n_only_d150": int((mb & ~ma).sum())}
        if m.any():
            pa = np.angle(Ha[m].astype(np.complex128))
            pb = np.angle(Hb[m].astype(np.complex128))
            dphi = np.angle(np.exp(1j * (pb - pa)))
            ra = np.abs(Ha[m].astype(np.complex128))
            rb = np.abs(Hb[m].astype(np.complex128))
            ratio = rb / np.maximum(ra, 1e-300)
            node.update({
                "phase_diff_rms_rad": _rms(dphi),
                "phase_diff_max_abs_rad": float(np.max(np.abs(dphi))),
                "amp_ratio_median": float(np.median(ratio)),
                "amp_ratio_p10": float(np.percentile(ratio, 10)),
                "amp_ratio_p90": float(np.percentile(ratio, 90)),
                "amp_ratio_spread_p90_over_p10": float(
                    np.percentile(ratio, 90) / max(np.percentile(ratio, 10),
                                                   1e-300)),
                "omega_q_diff_rms_rad_s": _rms(ob[m] - oa[m]),
                "omega_q_n_changed": int(np.count_nonzero(ob[m] != oa[m])),
                "omega_q_frac_changed": float(np.mean(ob[m] != oa[m])),
            })
        out[name] = node
    return out


def cstats(C):
    C = np.asarray(C, np.float64)
    return {"mean": float(np.mean(C)), "std": float(np.std(C)),
            "min": float(np.min(C)), "max": float(np.max(C)),
            "p01": float(np.percentile(C, 1)),
            "p50": float(np.percentile(C, 50)),
            "p99": float(np.percentile(C, 99)),
            "shape": [int(C.shape[0]), int(C.shape[1])]}


def cstats_diff(got, rec):
    keys = ("mean", "std", "min", "max", "p01", "p50", "p99")
    return {k: float(got[k] - rec[k]) for k in keys if rec and k in rec}


def identical_table(ref_scene, env, optics, box_bulk_max):
    """Which non-depth parameters are identical to the STAGE B record."""
    rows = {}
    rec = ref_scene or {}

    def row(name, got, want, tol=1e-12):
        node = {"this_run": jsonable(got), "stage_b_recorded": jsonable(want),
                "identical": None, "max_abs_diff": None}
        if want is not None:
            try:
                ga = np.atleast_1d(np.asarray(got, np.float64))
                wa = np.atleast_1d(np.asarray(want, np.float64))
                if ga.shape == wa.shape:
                    d = float(np.max(np.abs(ga - wa)))
                    node["max_abs_diff"] = d
                    node["identical"] = bool(d <= tol)
                else:
                    node["identical"] = False
                    node["shape_mismatch"] = [list(ga.shape), list(wa.shape)]
            except (TypeError, ValueError):
                node["identical"] = (got == want)
        rows[name] = node

    row("env.sun_altitude_deg", env["sun_altitude_deg"],
        (rec.get("env") or {}).get("sun_altitude_deg"))
    row("env.sunDir", env["sunDir"], (rec.get("env") or {}).get("sunDir"))
    row("env.clarity", env["clarity"], (rec.get("env") or {}).get("clarity"))
    row("env.turbid", env["turbid"], (rec.get("env") or {}).get("turbid"))
    row("optics.phase_kind", optics["phase_kind"],
        (rec.get("phase") or {}).get("kind"))
    row("optics.phase_n", optics["phase_n"],
        (rec.get("phase") or {}).get("n"))
    row("optics.phase_mu", optics["phase_mu"],
        (rec.get("phase") or {}).get("mu"))
    row("optics.kernel", optics["kernel"],
        (rec.get("optics") or {}).get("kernel"))
    row("deposit.box_bulk_max", int(box_bulk_max),
        ((rec.get("deposit") or {}).get("box_bulk_max_used")))
    row("deposit.w_cell", float(R15.WIDTH * R15.HEIGHT / (512.0 ** 2)),
        ((rec.get("deposit") or {}).get("w_cell")))
    for k in ("w_ballistic", "w_halo", "sigma_px_ff"):
        rows["optics." + k] = {
            "this_run": jsonable(optics[k]),
            "stage_b_recorded": None,
            "identical": None,
            "must_be_identical": False,
            "why": ("depth-dependent by construction (optical path length "
                    "L = 4.0 m -> 1.5 m): expected to change")}
    for k, v in rows.items():
        if "must_be_identical" not in v:
            v["must_be_identical"] = True
            v["why"] = "STAGE B parameter that the depth change must not touch"
    return {"pass": bool(all(v["identical"] is not False
                             for v in rows.values() if v["must_be_identical"])),
            "rows": rows}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-control", action="store_true",
                    help="skip the depth-4.0 STAGE B reproduction control")
    a = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    log("BRIEF %s" % BRIEF)
    log("output dir %s (out/v4-1-15 is READ ONLY for this job)" % OUT)
    log("does out/v4-1-15-extra exist already: %s" % os.path.isdir(OUT))
    png_existed = os.path.isfile(PNG)

    rec = {
        "record": "record-d150.json",
        "brief": BRIEF,
        "brief_sha256": (sha256_file(BRIEF) if os.path.isfile(BRIEF) else None),
        "job": ("BRIEF V4-1-15 D15: STAGE B (preset-5 pond) frame-0024, "
                "identical conditions, depth 4.0 m -> 1.5 m, ONE still"),
        "script": os.path.abspath(__file__),
        "script_sha256": sha256_file(os.path.abspath(__file__)),
        "host_time_utc": now_iso(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "writes_into_out_v4_1_15": False,
        "png_preexisting_at_start": bool(png_existed),
        "env": {},
    }

    try:
        # ---------------- 0. guards + reference ------------------------
        g0 = guard_hashes()
        ref = read_reference()
        rec["reference"] = ref
        log("reference STAGE B still: %s"
            % ((ref.get("still_file") or {}).get("sha256")))
        chk = check_setup(ref)
        rec["setup_check"] = chk
        log("setup check (all non-depth parameters == STAGE B): %s"
            % ("PASS" if chk["pass"] else "FAIL"))

        ref_sha = (ref.get("control_reference") or {}).get("sha256")
        if not a.no_control and not ref_sha:
            raise RuntimeError(
                "STOP: no reference hash for the STAGE B frame-0024 still "
                "(out/v4-1-15 checkpoint and still both unreadable); the "
                "identical-conditions claim cannot be verified")

        if not chk["pass"]:
            raise RuntimeError(
                "STOP: the imported driver / STAGE B record does not describe "
                "the expected STAGE B scene; refusing to relabel something "
                "else as 'identical conditions' (see setup_check)")

        # ---------------- 1. control: reproduce STAGE B at 4.0 m -------
        control = None
        control_pass = None
        if a.no_control:
            log("control SKIPPED (--no-control): identity is not verified")
        else:
            control = render(D_BASE, "control-d400")
            cs = cstats(control["frame"]["C"])
            recorded_cs = ((ref.get("checkpoint") or {})
                           .get("recorded_C_frame_stats") or {})
            control_pass = bool(control["png_sha256"] == ref_sha)
            rec["control"] = {
                "purpose": ("re-render STAGE B's own frame 24 at depth 4.0 m "
                            "through the same call chain and compare the PNG "
                            "byte for byte"),
                "depth_m": D_BASE,
                "png_sha256": control["png_sha256"],
                "reference_sha256": ref_sha,
                "reference_sources": (ref.get("control_reference") or
                                      {}).get("sources"),
                "byte_identical": control_pass,
                "still_shape": control["still_shape"],
                "t_build_s": control["t_build_s"],
                "t_render_s": control["t_render_s"],
                "C_stats": cs,
                "C_stats_recorded": recorded_cs,
                "C_stats_abs_diff_vs_recorded": cstats_diff(cs, recorded_cs),
                "energy": R15.energy_record(control["frame"]["det_info"]),
                "energy_recorded": (ref.get("checkpoint") or
                                     {}).get("recorded_energy"),
                "deposit_wall_s_recorded": (ref.get("checkpoint") or
                                             {}).get("recorded_deposit_wall_s"),
                "state": state_node(control["state"]),
                "env": env_node(control["env"]),
                "optics": optics_node(control["optics"]),
            }
            if not control_pass:
                raise RuntimeError(
                    "STOP: the depth-4.0 control did not reproduce STAGE B's "
                    "frame-0024 still byte for byte (%s != %s); the identical-"
                    "conditions setup cannot be reproduced, so no 1.5 m still "
                    "is fabricated"
                    % (control["png_sha256"], ref_sha))
            log("control PASS: the re-rendered depth-4.0 still is byte "
                "identical to STAGE B's recorded frame-0024")

        # ---------------- 2. the delta: depth 1.5 m -------------------
        delta = render(D_DELTA, "d150")
        assert delta["still_shape"] == [540, 736, 3], delta["still_shape"]
        assert int(delta["state"].N) == int(R15.SCN[STAGE]["N"])
        assert abs(float(delta["state"].S) - float(R15.SCN[STAGE]["L"])) < 1e-12

        # ---------------- 3. write the still --------------------------
        sha_written = R15.save_png(PNG, delta["frame"]["still"])
        if sha_written != delta["png_sha256"]:
            raise RuntimeError("written PNG hash differs from the in-memory "
                               "encoding: %s != %s"
                               % (sha_written, delta["png_sha256"]))
        from PIL import Image
        im = Image.open(PNG)
        im.load()
        png_facts = {"path": PNG, "bytes": os.path.getsize(PNG),
                     "sha256": sha_written, "size": list(im.size),
                     "mode": im.mode}
        log("wrote %s  %dx%d  %d bytes  sha %s"
            % (PNG, im.size[0], im.size[1], png_facts["bytes"], sha_written))

        # the brief's required record fields first, so that a later failure
        # can never leave the still without its facts
        info = delta["frame"]["det_info"]
        rec["deliverable"] = {
            "png": png_facts,
            "depth_m": D_DELTA,
            "base_depth_m": D_BASE,
            "patch_m": float(delta["state"].S),
            "window_m_stage_b_rule": [float(R15.SCN[STAGE]["window_m"][0]),
                                      float(R15.SCN[STAGE]["window_m"][1])],
            "grid": int(delta["state"].N),
            "grid_cell_m": float(delta["state"].S / delta["state"].N),
            "seed": int(R15.SCN[STAGE]["seed"]),
            "field_time_s": float(delta["frame"]["t_s"]),
            "frame": FRAME, "fps": float(R15.FPS),
            "state_build_s": delta["t_build_s"],
            "frame_render_s": delta["t_render_s"],
            "energy_rel_err": float(info["rel_err"]),
            "C_std": float(np.std(delta["frame"]["C"])),
        }

        # ---------------- 4. record ----------------------------------
        rec.update({
            "status": "ok",
            "delta": {
                "depth_m": D_DELTA,
                "base_depth_m": D_BASE,
                "only_this_changed": ("SCN['sub']['depth'] (4.0 -> 1.5) in "
                                      "memory; the wave state, env depthM, "
                                      "z_refraction and the optics path "
                                      "length follow from it"),
            },
            "scene": {
                "stage": STAGE,
                "kind": R15.SCN[STAGE]["kind"],
                "preset": delta["spec"].get("preset"),
                "water_label": delta["spec"].get("label"),
                "a_rgb": [float(x) for x in delta["spec"]["a"]],
                "b_rgb": [float(x) for x in delta["spec"]["b"]],
                "U_ms": float(R15.SCN[STAGE]["U"]),
                "fetch_m": float(R15.SCN[STAGE]["F"]),
                "seed": int(R15.SCN[STAGE]["seed"]),
                "seed_equals_derive_seed": bool(
                    int(delta["state"].seed)
                    == int(derive_seed(R10.load_recipe(), 1))),
                "seed_on_state": int(delta["state"].seed),
                "patch_m": float(delta["state"].S),
                "window_m_stage_b_rule": [float(R15.SCN[STAGE]["window_m"][0]),
                                          float(R15.SCN[STAGE]["window_m"][1])],
                "magnification": 1.0,
                "grid": int(delta["state"].N),
                "grid_cell_m": float(delta["state"].S / delta["state"].N),
                "grid_cells_per_m": float(delta["state"].N / delta["state"].S),
                "grid_is_stage_b_grid": bool(
                    int(delta["state"].N) == int(EXPECT_SCN["N"])),
                "patch_doublings": int(delta["state"].patch_doublings),
                "lambda_p_m": float(delta["state"].laws["lambda_p_m"]),
                "field_time_s": float(delta["frame"]["t_s"]),
                "frame": FRAME, "fps": float(R15.FPS),
                "field_time_rule": "t = frame / FPS = 24 / 24 = 1.000 s",
                "t_s_recorded_stage_b": (ref.get("checkpoint") or
                                         {}).get("recorded_t_s"),
                "composite": {"A": float(R15.A_FROZEN),
                              "Blow_sigma_px": float(R15.B_SIGMA_FROZEN),
                              "E": float(R15.E_FROZEN),
                              "Q_det_reference_only": int(R15.Q_DET),
                              "note": ("A/Blow/E are the frozen form; Q = 25 "
                                       "is the point-splat reference's "
                                       "quadrant count in STAGE B and is not "
                                       "used by the area deposit")},
                "deposit": {"module": "transport_mesh.mesh_C (VERBATIM)",
                            "box_bulk_max": int(delta["box_bulk_max"]),
                            "w_cell": float(R15.WIDTH * R15.HEIGHT
                                            / float(delta["state"].N ** 2)),
                            "kernel": "none", "jitter": "none",
                            "band_limit": "none",
                            "mean_normalisation": "none"},
                "output_px": [int(R15.WIDTH), int(R15.HEIGHT)],
                "c_grid": [int(delta["frame"]["C"].shape[0]),
                           int(delta["frame"]["C"].shape[1])],
            },
            "timing": {
                "state_build_s": delta["t_build_s"],
                "frame_render_s": delta["t_render_s"],
                "control_state_build_s": (control["t_build_s"]
                                          if control else None),
                "control_frame_render_s": (control["t_render_s"]
                                           if control else None),
                "total_job_s": None,
            },
            "energy": {
                "rel_err": float(info["rel_err"]),
                "cell_energy_total": float(info["energy_cell_total"]),
                "deposited_total": float(info["energy_deposited"]),
                "energy_guard": float(info["energy_guard"]),
                "offframe_fraction_area": float(info["offframe_fraction_area"]),
                "offframe_fraction_energy": float(
                    info["offframe_fraction_energy"]),
                "degenerate_cells": int(info["n_degenerate_cell"]),
                "n_bulk_triangles": int(info["n_bulk"]),
                "n_tail_triangles": int(info["n_tail"]),
                "w_cell": float(info["w_cell"]),
                "energy_record_full": R15.energy_record(info),
                "stage_b_reference_rel_err": (
                    (ref.get("checkpoint") or {}).get("recorded_energy") or
                    {}).get("rel_err"),
            },
            "C_stats": cstats(delta["frame"]["C"]),
            "C_stats_stage_b_reference": (ref.get("checkpoint") or
                                          {}).get("recorded_C_frame_stats"),
            "env": env_node(delta["env"]),
            "optics": optics_node(delta["optics"]),
            "state": state_node(delta["state"]),
            "identical_parameters": identical_table(
                (ref.get("checkpoint") or {}).get("scene"),
                env_node(delta["env"]), optics_node(delta["optics"]),
                delta["box_bulk_max"]),
            "field_vs_stage_b": {
                "why_measured_not_assumed": (
                    "depth enters the generator itself (mode density -> "
                    "amplitudes, omega(k, h) -> phase speed), so whether the "
                    "1.5 m field is bit-identical to the 4.0 m field is a "
                    "measurement, not an assumption: see 'draws' (same seeded "
                    "complex draw / same phase), 't_0' and 't_frame24'.  At "
                    "h = 1.5 m the energetic band is still deep water, so the "
                    "field is expected to be almost unchanged and the depth "
                    "change to act through refraction and the optical path."),
                "draws": draw_cmp(control, delta) if control else None,
                "t_0": field_cmp(control, delta, 0.0) if control else None,
                "t_frame24": field_cmp(control, delta,
                                       float(delta["frame"]["t_s"]))
                if control else None,
            },
            "image": {"png": png_facts,
                      "in_memory_png_sha256": delta["png_sha256"],
                      "encoder": "PIL Image.fromarray(...).save(format='PNG'), "
                                 "identical to run_v4_1_15.save_png"},
            "image_metrics": {k: jnum(v) for k, v in
                              delta["frame"]["metrics"].items()},
            "g1": float(delta["frame"]["g1"]),
        })
        g1h = guard_hashes()
        rec["guard"] = {"files": GUARD_FILES, "before": g0, "after": g1h,
                        "unchanged": bool(g0 == g1h)}
        rec["timing"]["total_job_s"] = float(time.time() - T0)
        if control:
            rec["control_proof"] = (
                "the depth-4.0 control reproduces out/v4-1-15/"
                "frame-0024-sub.png byte for byte, so environment, window, "
                "grid, seed, field time, deposit and composite are identical "
                "to STAGE B by measurement, not by assertion")
        rec["status_note"] = ("ONE still produced.  48 frames / video were NOT "
                              "produced (not asked).")
        write_json(REC, jsonable(rec))
        log("wrote %s" % REC)
        log("guard (existing .py unchanged): %s" % rec["guard"]["unchanged"])
        log("DONE: %s" % PNG)
        return 0

    except Exception as e:
        rec["status"] = ("incomplete_record"
                         if rec.get("deliverable") else "stopped")
        rec["stop_reason"] = "%s: %s" % (type(e).__name__, e)
        rec["stopped_utc"] = now_iso()
        rec.setdefault("timing", {})
        rec["timing"]["total_job_s"] = float(time.time() - T0)
        rec["png_written"] = bool(os.path.isfile(PNG) and not png_existed)
        try:
            write_json(REC, jsonable(rec))
            log("STOP recorded in %s" % REC)
        except Exception as e2:                      # pragma: no cover
            log("could not write record: %s" % e2)
        log("STOP: %s" % rec["stop_reason"])
        raise


if __name__ == "__main__":
    sys.exit(main())
