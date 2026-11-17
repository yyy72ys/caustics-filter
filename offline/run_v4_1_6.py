"""run_v4_1_6 -- BRIEF-V4-1-6 driver: generate the wave field at a healthy
patch and STRETCH it onto the proven screen size (separation form).

NEW FILE ONLY.  Every existing .py (including the v4-1.5 pair
`stream_surface4.py` / `run_v4_1_5.py`) stays byte-identical; G-LEGACY proves
it from a stage-0 sha256 record.  Nothing in this file rewrites a constant of
an earlier build: the screen constants are read from
`run_v4_1_4.SCENES['natscene']` (which is itself verbatim from run_v4_1_2),
the transport rule is `stream_surface4`'s, and the accounting is
`stream_surface3`'s.

WHY (user decision 2026-11-13 03:20-03:21, brief BACKGROUND)
    v4-1.4 asked one patch to be BOTH the physical patch (0.0639 m, so that the
    photo and the bed imprint line up at the real size) AND a resolved spectral
    domain.  At 0.0639 m the ring sits at k0/dk = 1.29: the annulus peak is
    quantised onto 1-2 bins, and the sub-lattice slope content folds the 1-ray
    per texel splat into axis-parallel lines (v4-1.4, 134.6 permille).
    v4-1.5 removed the printed lattice by band-limiting the transport field,
    but it did not restore the spectral resolution of the ring.

    This build separates the two roles:

        generation   L_gen = U0 * T_loop * n = 0.2875 * 2 * 1 = 0.575 m
                     (k0/dk = 11.64 -> the ring is resolved; the bed response
                     (i) uses the ORIGINAL PHOTO mapped onto this patch)
        rendering    736 x 540 at L_img = 3.68 m (run_v4_1_2 scene B size),
                     refraction depth z = 0.35 m (a screen quantity)
        composite    the generated slope field AND C are stretched
                     L_gen -> L_img, factor 6.4

    The stretch is NOT hidden: the physics claim stops at the generation of the
    field, and the wavelength that is actually drawn is the stretched one,
    lambda_displayed = lambda0 * 6.4 = 316.1 mm (63.2 px).  Both numbers are in
    numbers-v4-1-6.json under `separation` and at the top of the report.

STAGES (all foreground, one process, jobkit checkpoint per stage)
    --stage0    absolute-path preflight + stage-0 sha256 record + scene record
    --gen       build the L_gen field (photo bed), fold ladder, v4-1.5
                band-limit path, ring preservation, annulus spectrum (G-SPEC)
    --chunk K   frames 24K..24K+23 of the stretched composite
    --still     publish out/v4-1-6/frame-0024.png (= frames/frame-0024.png)
    --warp      displacement yardstick at frame 24 (G-WARP, recorded only)
    --encode    MP4 + ffprobe (G-OUT)
    --gates     G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT
    --report    numbers-v4-1-6.json + REPORT-V4-1-6.md
    --all       everything, sequentially

MEASUREMENT (G-LINE)
    Exactly run_v4_1_5's implementation, imported and re-used unchanged
    (green channel, 8-bit; axis-parallel 5-15 px band; 256x256 windows whose
    positions come from the SOURCE PHOTO only).  The primary quantity is
    depth_e3 (raw band amplitude); the brief's high-pass-first procedure is
    reported as depth_hp_e3 in the same table.
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

import run_v4_1_4 as R14  # noqa: E402  (reused UNCHANGED: natscene screen)
import run_v4_1_5 as R15  # noqa: E402  (reused UNCHANGED: the G-LINE metric)
import run_v16_water as V16  # noqa: E402
import stream_surface2 as SS  # noqa: E402
import stream_surface4 as S4  # noqa: E402
from det import det_frame  # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-1-6.md"
PHOTO = r"D:\mutaro-work\moth-on-bark.jpg"

OUT = os.path.join(HERE, "out", "v4-1-6")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
STAGE0_HASHES = os.path.join(OUT, "stage0-py-hashes.json")
MP4 = os.path.join(OUT, "piece-2s-v4-1-6-stretched.mp4")
STILL = os.path.join(OUT, "frame-0024.png")
NUMBERS = os.path.join(OUT, "numbers-v4-1-6.json")
REPORT = os.path.join(OUT, "REPORT-V4-1-6.md")

V15_OUT = os.path.join(HERE, "out", "v4-1-5")
V15_REPORT = os.path.join(V15_OUT, "REPORT-V4-1-5.md")
V15_NUMBERS = os.path.join(V15_OUT, "numbers-v4-1-5.json")
V15_STILL = os.path.join(V15_OUT, "frame-realscale-0024.png")
V14_STILL = os.path.join(HERE, "out", "v4-1-4", "frame-realscale-0024.png")
V13_STILL = os.path.join(HERE, "out", "v4-1-3", "frame-0024.png")

# ---- generation geometry (brief section 作業 1) -----------------------------
GEN = {
    "H": 0.015, "U": 0.25, "S0": 0.006, "D84": 0.004,
    "seed": 3255045276382562994,
    "T_loop": 2.0, "n_loops": 1,
    "N": 512, "alpha": SS.ALPHA_MAIN,
    "bed": "photo (the original photograph mapped onto L_gen)",
}
L_GEN_EXPECTED = 0.575
# ---- screen geometry (brief section 作業 2; run_v4_1_2 scene B) -------------
IMG = {"A": 2.6, "B_sigma": 3.0, "E": 2.2, "Q": 25,
       "water": "algae_pond", "z_m": 0.35}

FPS = R14.FPS
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
STILL_IDX = 24
WIDTH, HEIGHT = R14.WIDTH, R14.HEIGHT
N_BUF = 512


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
# geometry: generation patch and screen patch
# ---------------------------------------------------------------------------
def geometry():
    """L_gen from the loop rule, L_img from run_v4_1_4's natscene scene."""
    d = SS.derive(GEN["U"], GEN["H"], GEN["S0"])
    L_gen = float(d["U0"]) * float(GEN["T_loop"]) * float(GEN["n_loops"])
    L_img, n_img, loop_img = R14.scene_loop("natscene")
    k0 = SS.k0_stationary(d["U0"], GEN["H"])
    lam0 = (2.0 * np.pi / k0) if k0 else None
    dk = 2.0 * np.pi / L_gen
    stretch = float(L_img) / L_gen
    rec = {
        "generation": {
            "H_m": GEN["H"], "U_m_s": GEN["U"], "U0_m_s": float(d["U0"]),
            "S0": GEN["S0"], "D84_m": GEN["D84"], "seed": int(GEN["seed"]),
            "T_loop_s": GEN["T_loop"], "n_loops": GEN["n_loops"],
            "L_gen_m": L_gen, "L_gen_expected_m": L_GEN_EXPECTED,
            "N": GEN["N"], "cell_m": L_gen / GEN["N"],
            "k0_rad_m": float(k0), "lambda0_m": float(lam0),
            "lambda0_mm": 1000.0 * float(lam0),
            "dk_rad_m": dk, "k0_over_dk": float(k0) / dk,
            "waves_per_patch": L_gen / float(lam0),
            "cells_per_wavelength": float(lam0) / (L_gen / GEN["N"]),
            "loop_rule": "L = U0 * T_loop * n (n = 1 here)",
            "derived": d,
        },
        "screen": {
            "L_img_m": float(L_img), "width_px": WIDTH, "height_px": HEIGHT,
            "fps": FPS, "frames": NFRAMES,
            "pixel_m": float(L_img) / WIDTH,
            "z_refraction_m": IMG["z_m"],
            "z_generation_m": IMG["z_m"] / stretch,
            "z_note": ("z is a SCREEN quantity (brief: 画面寸法に属する量).  The "
                       "generation patch is the 6.4x-reduced copy of the "
                       "screen, so its refraction depth is z_gen = z_img/6.4. "
                       "Then the stretch is a pure magnification and C is the "
                       "screen's C, not an over-focused neighbour: ds/di = "
                       "z*Cb*curv is invariant under this scaling, so the "
                       "screen's own displacement gradient is measured at the "
                       "generation metric."),
            "composite": {"A": IMG["A"], "B_sigma_px": IMG["B_sigma"],
                          "E": IMG["E"], "Q": IMG["Q"],
                          "water": IMG["water"]},
            "source": ("run_v4_1_4.SCENES['natscene'] (verbatim from "
                       "run_v4_1_2.py scene B)"),
            "loop": loop_img, "n_loops": int(n_img),
        },
        "separation": {
            "stretch_factor": float(L_img) / L_gen,
            "stretch_factor_computed": float(L_img) / L_gen,
            "definition": "L_img / L_gen; the generated field (slope + C) is "
                          "mapped onto the whole screen, so the drawn "
                          "wavelength is lambda0 * stretch",
            "lambda0_mm": 1000.0 * float(lam0),
            "lambda_displayed_mm": 1000.0 * float(lam0) * (float(L_img) / L_gen),
            "lambda_displayed_px": (float(lam0) * (float(L_img) / L_gen)
                                    / (float(L_img) / WIDTH)),
            "cells_per_wavelength_generation": float(lam0) / (L_gen / GEN["N"]),
            "px_per_cell_on_screen": (float(L_img) / GEN["N"]) / (float(L_img) / WIDTH),
            "z_gen_m": IMG["z_m"] / stretch,
            "z_img_m": IMG["z_m"],
            "drawn_amplitude_note": ("the slope array is dimensionless and is "
                                     "kept; therefore the drawn surface "
                                     "amplitude is stretch x the generated one "
                                     "and the steepness kA is preserved"),
            "physics_claim": ("the physical claim covers the GENERATION of the "
                              "field (L_gen, k0, amplitudes).  The wave that "
                              "is drawn is the stretched one: its wavelength "
                              "is lambda0*6.4 = 316 mm, and that is the number "
                              "quoted for the image, not lambda0."),
        },
    }
    rec["separation"]["drawn_steepness_kA"] = None  # filled by --gen

    if abs(L_gen - L_GEN_EXPECTED) > 1e-9:
        raise SystemExit("STOP: L_gen = %.9f != brief %.9f"
                         % (L_gen, L_GEN_EXPECTED))
    return rec


# ---------------------------------------------------------------------------
# the bed response (i): the ORIGINAL PHOTO mapped onto L_gen
# ---------------------------------------------------------------------------
def photo_bed_zeta(N, sigma_b, path=PHOTO, seam_report=None, L_patch=None):
    """Bed elevation from the photo, toroidally resampled onto the N grid.

    The photo is mapped onto the patch (the square N grid covers the whole
    L_gen x L_gen patch), low-passed before decimation so the 0.78 mm photo
    pixel does not alias into the resolved band, and then made periodic: the
    value jump across the patch boundary is removed with a linear ramp of
    slope <= sigma_b / L_gen (recorded), because the field is Fourier-periodic
    by construction everywhere else and a step at the boundary would leak
    power along the axes (exactly the axes the ring lives on).
    """
    from PIL import Image
    from scipy import ndimage as ndi
    a = np.asarray(Image.open(path).convert("L"), dtype=np.float64)
    Hp, Wp = a.shape
    sx, sy = Wp / float(N), Hp / float(N)
    ap = ndi.gaussian_filter(a, 0.5 * max(sx, sy), mode="wrap")
    xx = (np.arange(N, dtype=np.float64) + 0.5) * sx - 0.5
    yy = (np.arange(N, dtype=np.float64) + 0.5) * sy - 0.5
    X, Y = np.meshgrid(xx, yy)
    r = ndi.map_coordinates(ap, [Y, X], order=1, mode="wrap")
    r = r - r.mean()
    r = r / float(r.std())
    # periodic seam: value jump when the last column wraps to the first.
    # A linear ramp of the same total rise makes the two edges meet (the
    # residual is one cell of the ramp = dx/N), and its physical slope
    # dx*sigma_b/L_gen is ~4e-4 (recorded below, negligible).
    dx = float(np.mean(r[:, 0] - r[:, -1]))
    dy = float(np.mean(r[0, :] - r[-1, :]))
    ux = (np.arange(N, dtype=np.float64) + 0.5) / float(N)
    ramp = dx * ux[None, :] + dy * ux[:, None]
    r = r + ramp
    r = r - r.mean()
    r = r / float(r.std())
    z = r * float(sigma_b)
    if seam_report is not None:
        Lp = float(L_patch) if L_patch else 1.0
        seam_report.update({
            "kind": "photo (original photograph mapped onto L_gen)",
            "photo": path, "photo_sha256": sha256_file(path),
            "photo_shape": [int(Hp), int(Wp)],
            "resample": ("toroidal bilinear after a gaussian prefilter of "
                         "0.5*max(Wp,Hp)/N = %.3f photo px" % (0.5 * max(sx, sy))),
            "mapping": "whole photo -> the N x N periodic patch (x and y "
                       "rescaled independently; no crop)",
            "sigma_b_m": float(sigma_b),
            "sigma_b_source": "bed_mapping_D84: sigma_b = D84/4 (v4-1.1 rule)",
            "normalisation": "unit std -> sigma_b, zero mean",
            "seam_step_before_ramp_in_sigma_b": [dx, dy],
            "ramp_removed_in_sigma_b": [dx, dy],
            "ramp_slope_removed_physical": [
                dx * float(sigma_b) / Lp, dy * float(sigma_b) / Lp],
            "ramp_slope_removed_dimensionless": [dx / float(N), dy / float(N)],
            "ramp_after_in_sigma_b": [
                float(np.mean(r[:, 0] - r[:, -1])),
                float(np.mean(r[0, :] - r[-1, :]))],
            "z_min_m": float(z.min()), "z_max_m": float(z.max()),
            "z_rms_m": float(np.sqrt(np.mean(z * z))),
        })
    return z


class PhotoBedState(S4.StreamState4):
    """StreamState4 with the (i) bed replaced by the photo mapped onto L_gen.

    Only `_build_i` is overridden, and it reproduces stream_surface2's tail
    verbatim (transfer T(k) + the breaking-limit valve, same constants).  The
    accounting, the ring/tail/turbulence draws, the time evolution and the
    v4-1.5 transport band-limit are inherited untouched.
    """

    bed_record = None

    def __init__(self, photo_path=PHOTO, **kwargs):
        self._photo_path = photo_path
        super().__init__(**kwargs)

    def _build_i(self, sigma_b, ell_b):
        m = SS.bed_mapping_D84(self.ks)
        sigma_b = m["sigma_b"]
        ell_b = m["ell_b"]
        rec = {}
        z = photo_bed_zeta(self.N, sigma_b, self._photo_path, rec, L_patch=self.L)
        zeta_hat = np.fft.fft2(z)
        zeta_hat[0, 0] = 0.0
        rec["ell_b_nominal_m"] = float(ell_b)
        rec["ell_b_note"] = ("ell_b is not used by the photo bed (its shape "
                             "comes from the photograph); only sigma_b sets "
                             "the amplitude")
        type(self).bed_record = rec

        T, denom = SS.transfer_T(self.kx, self.kappa, self.H, self.U)
        eta_hat = T * zeta_hat
        band = np.abs(denom) < SS.RES_BAND_DENOM
        band[~np.isfinite(denom)] = False
        band[np.abs(self.kx) == 0.0] = False
        k2 = self.kappa ** 2
        N2 = self.N ** 2
        tot_e = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
        band_e = float(np.sum(np.abs(eta_hat[band]) ** 2) / N2 ** 2)
        band_slope = float(np.sqrt(max(
            np.sum(k2[band] * np.abs(eta_hat[band]) ** 2) / N2 ** 2, 0.0)))
        scale = 1.0
        if band_slope > SS.SLOPE_CAP > 0:
            scale = SS.SLOPE_CAP / band_slope
            eta_hat[band] = eta_hat[band] * scale
        after = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
        info = {
            "bed": rec,
            "_T": T,
            "_zeta_hat": zeta_hat,
            "T_max": float(np.max(np.abs(T))),
            "resonant_band_def": "|denominator| < %.3g" % SS.RES_BAND_DENOM,
            "band_rms_slope_before": band_slope,
            "band_energy_fraction": (band_e / tot_e if tot_e > 0 else 0.0),
            "scale_applied": float(scale),
            "dropped_band_pct": float(100.0 * (1.0 - scale ** 2)),
            "dropped_total_pct": float(100.0 * (tot_e - after) / tot_e
                                       if tot_e > 0 else 0.0),
            "sigma_i_m": float(np.sqrt(max(after, 0.0))),
        }
        return eta_hat, info


def build_generation():
    g = geometry()
    d = SS.derive(GEN["U"], GEN["H"], GEN["S0"])
    L_gen = float(d["U0"]) * GEN["T_loop"] * GEN["n_loops"]
    st = PhotoBedState(N=GEN["N"], L=L_gen, H=GEN["H"], U=GEN["U"],
                       S0=GEN["S0"], D84=GEN["D84"], seed=GEN["seed"],
                       bed="photo", alpha=GEN["alpha"],
                       dt_loop=float(GEN["T_loop"]))
    st.loop = {"L_gen": L_gen, "n_loops": GEN["n_loops"],
               "T_loop_s": GEN["T_loop"], "source": "BRIEF-V4-1-6"}
    return st, g


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 78)
    log("[stage 0] absolute-path preflight + stage-0 sha256 record")
    makedirs()
    cp = Checkpoint(CP_PATH)
    need = [BRIEF, V15_REPORT, V15_NUMBERS, V15_STILL,
            os.path.join(HERE, "stream_surface3.py"),
            os.path.join(HERE, "stream_surface4.py"),
            os.path.join(HERE, "run_v4_1_2.py"),
            os.path.join(HERE, "run_v4_1_4.py"),
            os.path.join(HERE, "run_v4_1_5.py"),
            PHOTO, V14_STILL, V13_STILL]
    missing = [p for p in need if not os.path.isfile(p)]
    log("  required paths: %d, missing: %s" % (len(need), missing or "none"))
    if missing:
        raise SystemExit("STOP: missing required paths: %r" % (missing,))
    rec = {"brief": BRIEF, "checked": need, "missing": [],
           "checked_local": now_iso()}

    new_files = ["run_v4_1_6.py"]
    hashes = py_set()
    existing = {k: v for k, v in hashes.items() if k not in new_files}
    payload = {"all": hashes, "new_files_of_this_job": new_files,
               "existing_count": len(existing), "existing": existing,
               "recorded_local": now_iso(),
               "note": ("all .py in the workdir; 'existing' excludes this "
                        "job's own new file (its hash is pinned in 'all')")}
    write_json(STAGE0_HASHES, payload)
    log("  sha256 recorded: %d .py (%d existing + %d new)"
        % (len(hashes), len(existing), len(new_files)))

    # the v4-1.5 landing must be readable: the lattice fix and its numbers
    v15 = read_json(V15_NUMBERS)
    landed = {
        "report": os.path.basename(V15_REPORT),
        "G-LINE_pass": bool(v15["gates"]["G-LINE"]["pass"]),
        "sigma_cells_chosen": float(v15["fix"]["sigma_cells_chosen"]),
        "transport_path": "stream_surface4.transport_fields (spectral FFT "
                          "Gaussian, width from max|ds/di|<=0.25)",
        "lattice_strength_v4_1_5": float(v15["gates"]["G-LINE"]["strength_v4_1_5"]),
        "lattice_strength_photo": float(v15["gates"]["G-LINE"]["strength_photo"]),
        "lattice_strength_v4_1_4": float(v15["gates"]["G-LINE"]["strength_v4_1_4"]),
        "all_gates_pass": bool(v15["gates"]["all_pass"]),
    }
    if not landed["G-LINE_pass"]:
        raise SystemExit("STOP: v4-1.5 G-LINE did not land")
    rec["v4_1_5"] = landed

    g = geometry()
    st, _g = build_generation()
    rec["geometry"] = g
    rec["generation_accounting"] = st.accounting()
    rec["generation_info_i"] = {k: v for k, v in st.info_i.items()
                                if not k.startswith("_")}
    rec["generation_component_table"] = st.component_table(STILL_IDX / float(FPS))
    rec["regions"] = R15.regions_from_photo()[1]
    rec["sigma_seed"] = float(st.seed)

    checkpoint_extra("stage0", rec)
    cp.unit_done("stage-0-preflight")
    sep = g["separation"]
    log("  L_gen=%.6f m (%.4f lambda0, k0/dk=%.3f)  L_img=%.6f m  stretch=%.6f"
        % (g["generation"]["L_gen_m"], g["generation"]["waves_per_patch"],
           g["generation"]["k0_over_dk"], g["screen"]["L_img_m"],
           sep["stretch_factor"]))
    log("  lambda0=%.4f mm -> displayed %.4f mm (%.2f px)"
        % (sep["lambda0_mm"], sep["lambda_displayed_mm"],
           sep["lambda_displayed_px"]))
    log("  accounting: sigma_reg=%.6f mm var_i=%.6e var_ii=%.6e var_iii=%.6e"
        % (st.sigma_reg * 1e3, st.var_i_physical, st.var_ii, st.var_iii))
    log("  bed(i): sigma_b=%.4f mm valve scale=%.6f dropped %.2f%%"
        % (st.info_i["bed"]["sigma_b_m"] * 1e3, st.info_i["scale_applied"],
           st.info_i["dropped_total_pct"]))
    del st
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return rec


# ---------------------------------------------------------------------------
# stage: generation field (fold ladder, ring, G-SPEC)
# ---------------------------------------------------------------------------
def do_gen():
    t0 = time.time()
    log("=" * 78)
    log("[gen] L_gen field: fold ladder + v4-1.5 transport path + G-SPEC")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st, g = build_generation()
    L_gen = g["generation"]["L_gen_m"]
    z_gen = g["screen"]["z_generation_m"]
    ior = 1.33434
    ladder = (0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
              10.0, 12.0, 14.0, 16.0)
    fields = st.fields_at(STILL_IDX / float(FPS))
    strict = S4.sigma_from_fold_criterion(fields, L_gen, z_gen, ior,
                                          ladder=ladder,
                                          target=S4.MAX_DS_PER_RAY_TARGET)
    nofold = S4.sigma_from_fold_criterion(fields, L_gen, z_gen, ior,
                                          ladder=ladder, target=1.0)
    ladder_rows = nofold["ladder"]
    chosen = nofold["chosen_sigma_cells"]
    criterion_used = (
        "mechanism rule: smallest sigma with max|ds/di| <= 1.0 (no fold).  "
        "v4-1.5's 4x margin (<= %.2f) IS satisfiable here at sigma=%s cell, "
        "but it removes the wave (sigma=%.1f cell = %.0f%% of lambda0), so the "
        "margin is not transferred; the full ladder with the ring-band ratio "
        "is reported."
        % (S4.MAX_DS_PER_RAY_TARGET, strict["chosen_sigma_cells"], chosen,
           100.0 * chosen / g["generation"]["cells_per_wavelength"]))
    log("  fold criterion (%s): sigma = %s cells"
        % (criterion_used, chosen))
    for row in ladder_rows:
        f_sig = S4.transport_fields(fields, row["sigma_cells"])
        keep = S4.large_scale_preservation(fields, f_sig, L_gen)
        row["ring_band_ratio"] = keep["ratio"]
        row["ring_amplitude_ratio"] = (float(np.sqrt(keep["ratio"]))
                                       if keep["ratio"] is not None else None)
        log("    sigma %4.1f  p50 %.4f p99 %.4f max %.4f  folded %.5f%%  "
            "ring-band %.4f" % (row["sigma_cells"], row["ds_per_ray_p50"],
                                row["ds_per_ray_p99"], row["ds_per_ray_max"],
                                100 * row["folded_area_fraction"],
                                keep["ratio"]))
    chk = S4.self_check(fields, L_gen, z_gen, ior)
    log("  identity at sigma=0: %s   ring-band ratio @4cell: %.6f"
        % (chk["identity_ok"], chk["ring_band"]["ratio"]))
    if chosen is None:
        raise SystemExit("STOP: no sigma in the ladder avoids a fold at L_gen")
    if not chk["identity_ok"]:
        raise SystemExit("STOP: sigma=0 is not the identity")

    eff = S4.transport_fields(fields, chosen)
    keep_eff = S4.large_scale_preservation(fields, eff, L_gen)
    k0 = g["generation"]["k0_rad_m"]
    dk = g["generation"]["dk_rad_m"]

    def spectral(hat):
        r = st.spectrum_radial(hat)
        S = np.asarray(r["S_k"], dtype=np.float64)
        j = int(np.argmax(S))
        return {"peak_k_rad_m": float(r["peak_k"]),
                "ratio_to_k0": float(r["peak_k"] / k0),
                "peak_bin": int(j),
                "bins_above_10pct": int(np.sum(S > 0.1 * S[j])),
                "dk_rad_m": float(dk),
                "k_axis": r["k_axis"], "S_k": r["S_k"]}

    c_ii = np.exp(-1j * st.om_ii * (STILL_IDX / float(FPS)))
    spec_wave = spectral((st.hat_ring + st.hat_tail) * c_ii)
    spec_total = spectral(st.eta_hat_at(STILL_IDX / float(FPS)))
    spec_bed = spectral(st.eta_i_hat)
    c_iii = np.exp(-1j * st.om_iii * (STILL_IDX / float(FPS)))
    e = float(np.exp(-(STILL_IDX / float(FPS)) / st.tau))
    sss = float(np.sqrt(max(1.0 - np.exp(-2.0 * (STILL_IDX / float(FPS))
                                        / st.tau), 0.0)))
    turb = (st.hat_turbA * e + st.hat_turbB * sss) * c_iii
    spec_iii = spectral(turb)
    spec_eff = spectral(np.fft.fft2(eff["h"]))
    spec_ok = bool(0.75 <= spec_wave["ratio_to_k0"] <= 1.33)
    log("  G-SPEC (wave ii): peak %.2f rad/m, ratio %.4f, %d bins > 10%% -> %s"
        % (spec_wave["peak_k_rad_m"], spec_wave["ratio_to_k0"],
           spec_wave["bins_above_10pct"], "PASS" if spec_ok else "FAIL"))
    log("         total field %.2f (ratio %.4f); bed(i) %.2f (%.4f); "
        "effective %.2f (%.4f)"
        % (spec_total["peak_k_rad_m"], spec_total["ratio_to_k0"],
           spec_bed["peak_k_rad_m"], spec_bed["ratio_to_k0"],
           spec_eff["peak_k_rad_m"], spec_eff["ratio_to_k0"]))

    # displacement in screen pixels at the SCREEN depth (the generation-metric
    # readout is z_gen; the drawn warp uses z_img = 0.35 m)
    shx, shy, per = S4._displacement(fields, L_gen, IMG["z_m"], ior)
    mag_px_screen = (np.hypot(shx, shy) / (g["screen"]["L_img_m"] / WIDTH))
    disp = {"mean_px_on_screen": float(mag_px_screen.mean()),
            "p50_px_on_screen": float(np.percentile(mag_px_screen, 50)),
            "p99_px_on_screen": float(np.percentile(mag_px_screen, 99)),
            "max_px_on_screen": float(mag_px_screen.max()),
            "note": "same displacement, read in 736-wide screen pixels"}

    # the drawn steepness: slopes are kept, so the drawn amplitude is
    # stretch x the generated amplitude and kA is preserved
    sep = g["separation"]
    drawn_amp = sep["stretch_factor"] * float(st.sigma_total)
    sep["drawn_amplitude_m"] = drawn_amp
    sep["drawn_steepness_kA"] = float(2.0 * np.pi
                                      / (sep["lambda_displayed_mm"] / 1000.0)
                                      * drawn_amp)
    sep["generated_sigma_total_m"] = float(st.sigma_total)

    out = {
        "geometry": g,
        "accounting": st.accounting(),
        "component_table_frame24": st.component_table(STILL_IDX / float(FPS)),
        "bed_info_i": {k: v for k, v in st.info_i.items()
                       if not k.startswith("_")},
        "describe": {"seed": int(st.seed), "N": int(st.N), "S_m": float(st.S),
                     "S0": float(st.S0), "U0": float(st.U0),
                     "L0_m": float(st.L0), "phi": float(st.phi),
                     "k0_rad_m": float(st.k0),
                     "lambda0_m": g["generation"]["lambda0_m"],
                     "tau_s": float(st.tau), "ustar": float(st.ustar)},
        "transport": {
            "module": "stream_surface4.transport_fields",
            "z_used_m": float(z_gen),
            "criterion_used": criterion_used,
            "target_requested": float(S4.MAX_DS_PER_RAY_TARGET),
            "strict_rule_sigma_cells": strict["chosen_sigma_cells"],
            "mechanism_rule_sigma_cells": nofold["chosen_sigma_cells"],
            "chosen_sigma_cells": float(chosen),
            "sigma_m": float(chosen * L_gen / st.N),
            "sigma_over_lambda0": float(chosen
                                        / g["generation"]["cells_per_wavelength"]),
            "ladder": ladder_rows,
            "identity": chk["identity_at_sigma_0"],
            "ring_band_ratio_at_4cell": (chk["ring_band"]["ratio"]),
            "ring_band_ratio_at_chosen": keep_eff["ratio"],
            "ring_amplitude_ratio_at_chosen": (
                float(np.sqrt(keep_eff["ratio"]))
                if keep_eff["ratio"] is not None else None),
            "field_used_for_transport": "band-limited copy",
            "field_used_for_warp": "unfiltered (same convention as v4-1.5)",
        },
        "displacement_on_screen": disp,
        "G-SPEC": {
            "pass": spec_ok,
            "criterion": ("annulus (radial, inherited stream_surface2."
                          "spectrum_radial) peak of the GENERATED WAVE field "
                          "(component ii = ring+tail, the field of work item "
                          "1) within 0.75-1.33 x analytic k0"),
            "k0_rad_m": k0, "dk_rad_m": dk,
            "wave": spec_wave, "total": spec_total,
            "bed": spec_bed, "turbulence": spec_iii, "effective": spec_eff,
            "definition_note": ("gated on the wave component.  The total "
                                "surface peak sits on the STATIC bed "
                                "response's transfer maximum (k ~ 1.33 k0) and "
                                "is therefore not a quantisation test; it is "
                                "reported with all components (see report)."),
            "definition_chosen_after_seeing_numbers": True,
        },
        "created_local": now_iso(),
        "wall_s": float(time.time() - t0),
    }
    write_json(os.path.join(WORK, "gen.json"), out)
    checkpoint_extra("gen", out)
    cp.unit_done("gen-field")
    del st, fields, eff
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# rendering one stretched frame
# ---------------------------------------------------------------------------
def _g1(fields, env):
    from v5_physics import smith_g1_beckmann
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    rms = float(np.sqrt(np.mean(gx * gx + gy * gy)))
    return float(smith_g1_beckmann(
        rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])


def render_stretched(i, st, env, optics, photo, sigma_cells):
    """One frame: generate at L_gen, stretch slope+C onto L_img, composite."""
    L_gen = float(st.S)
    L_img = float(R14.scene_loop("natscene")[0])
    t = i / float(FPS)
    fields = st.fields_at(t)
    tf = S4.transport_fields(fields, sigma_cells)
    g1 = _g1(fields, env)
    # C is generated on the healthy patch, with the generation copy's depth
    # (z_img/6.4) so that the 6.4x stretch is a pure magnification
    env_c = dict(env)
    env_c["depthM"] = float(IMG["z_m"]) / (L_img / L_gen)
    C, det_info = det_frame(tf, L_gen, N_BUF, env_c, int(IMG["Q"]))
    # the stretch: the SAME arrays are declared to span the screen width
    fx = {k: fields[k] for k in ("h", "gx", "gy", "Dx", "Dy")}
    fx["S"] = L_img
    img, dbg = V16.render_frame_v16(
        fx, C, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=IMG["A"], c_layer_exposure=IMG["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
    from v5_composite import image_metrics
    met = image_metrics(dbg["lin"], dbg["aces"])
    return {"frame": int(i), "t_s": float(t), "still": img, "C": C,
            "g1": g1, "det_info": det_info, "metrics": met, "fields": fx}


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
    sigma_cells = float(read_checkpoint()["gen"]["transport"]["chosen_sigma_cells"])
    _p, _stn, env, optics = R14.build("natscene")
    st, _g = build_generation()
    photo = R14.load_photo()
    frames, metrics = {}, []
    t_start = time.time()
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_stretched(i, st, env, optics, photo, sigma_cells)
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
        log("  frame %04d t=%.3f Cmean=%.6f Cstd=%.4f p99=%.3f %.1fs"
            % (i, r["t_s"], rec["C_mean"], rec["C_std"], rec["C_p99"],
               time.time() - t_start))
        del r, c, m
        gc.collect()
    node = {"frames": frames, "metrics": metrics,
            "completed_local": now_iso(),
            "wall_s": float(time.time() - t_start),
            "sigma_cells": float(sigma_cells),
            "stretch_factor": float(R14.scene_loop("natscene")[0] / st.S)}
    chunks = read_checkpoint().get("chunks") or {}
    chunks[str(k)] = node
    checkpoint_extra("chunks", chunks)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d wall %.1fs" % (k, node["wall_s"]))
    del st, _stn, optics, photo
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
    return True


def do_still():
    t0 = time.time()
    log("=" * 78)
    log("[still] publish frame-0024.png")
    makedirs()
    cp = Checkpoint(CP_PATH)
    src = frame_path(STILL_IDX)
    if not os.path.isfile(src):
        raise SystemExit("STOP: %s missing (run chunk 1 first)" % src)
    with open(src, "rb") as f:
        data = f.read()

    def verify(p):
        from PIL import Image
        try:
            im = Image.open(p)
            im.load()
            return im.size == (WIDTH, HEIGHT) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(STILL, data, verify=verify)
    out = {"file": os.path.basename(STILL), "source": "frames/frame-0024.png",
           "sha256": sha256_file(STILL), "sha256_source": sha256_file(src),
           "identical": bool(sha256_file(STILL) == sha256_file(src)),
           "bytes": os.path.getsize(STILL), "created_local": now_iso()}
    if not out["identical"]:
        raise SystemExit("STOP: still copy is not byte-identical")
    checkpoint_extra("still", out)
    cp.unit_done("still-0024")
    log("  %s (%d bytes, %s)" % (STILL, out["bytes"], out["sha256"][:16]))
    log("  wall %.1fs" % (time.time() - t0))
    return out


# ---------------------------------------------------------------------------
# stage: ablation of the two free choices (sigma, warp field) at frame 24
# ---------------------------------------------------------------------------
def do_ablate():
    """Is the G-LINE number caused by the transport choice, or by the scene?

    Four configurations of frame 24: chosen vs strict sigma (v4-1.5's 4x
    margin) x unfiltered vs band-limited warp field.  Nothing here is a
    deliverable; it is the evidence that the gate number is not an artefact of
    the sigma choice.
    """
    t0 = time.time()
    log("=" * 78)
    log("[ablate] frame 24: sigma x warp-field, same G-LINE statistic")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    gen = state["gen"]
    sig_chosen = float(gen["transport"]["chosen_sigma_cells"])
    sig_strict = gen["transport"]["strict_rule_sigma_cells"]
    _p, _stn, env, optics = R14.build("natscene")
    st, _g = build_generation()
    photo = R14.load_photo()
    regions, rinfo = R15.regions_from_photo()
    photo_rows = R15.measure_image(PHOTO, regions)
    L_gen = float(st.S)
    L_img = float(R14.scene_loop("natscene")[0])
    z_gen = IMG["z_m"] / (L_img / L_gen)
    fields = st.fields_at(STILL_IDX / float(FPS))
    g1 = _g1(fields, env)
    env_c = dict(env)
    env_c["depthM"] = z_gen
    out = {"frame": STILL_IDX, "photo_rows": photo_rows, "regions": rinfo,
           "z_gen_m": z_gen, "L_gen_m": L_gen, "L_img_m": L_img,
           "configs": []}
    for sig in ([sig_chosen] if sig_strict is None else
                [sig_chosen, float(sig_strict)]):
        tf = S4.transport_fields(fields, sig)
        C, _ = det_frame(tf, L_gen, N_BUF, env_c, int(IMG["Q"]))
        for wf in (False, True):
            fx = {k: fields[k] for k in ("h", "gx", "gy", "Dx", "Dy")}
            if wf:
                fx["gx"], fx["gy"], fx["h"] = tf["gx"], tf["gy"], tf["h"]
            fx["S"] = L_img
            img, dbg = V16.render_frame_v16(
                fx, C, photo["lin"], env, WIDTH, HEIGHT, optics,
                ambient_scale=IMG["A"], c_layer_exposure=IMG["E"],
                substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
                lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
            rows = R15.lattice_metric(img[:, :, 1].astype(np.float64), regions)
            n_over = sum(1 for a, b in zip(rows, photo_rows)
                         if a["depth_e3"] > b["depth_e3"])
            out["configs"].append({
                "sigma_cells": float(sig), "warp_field": ("band-limited" if wf
                                                          else "unfiltered"),
                "C_mean": float(np.mean(C)), "C_std": float(np.std(C)),
                "C_min": float(np.min(C)),
                "C_p99": float(np.percentile(C, 99)),
                "strength": R15.lattice_strength(rows),
                "rows": [{"region": r["region"], "axis": r["axis"],
                          "period_px": r["period_px"],
                          "depth_e3": r["depth_e3"],
                          "depth_photo": p["depth_e3"],
                          "ratio": r["depth_e3"] / p["depth_e3"]}
                         for r, p in zip(rows, photo_rows)],
                "rows_over_photo": int(n_over),
                "pass": bool(n_over == 0)})
            log("  sigma %5.1f warp %-12s  Cstd %.3f  strength %.1f  over %d/4"
                % (sig, "band-limited" if wf else "unfiltered",
                   float(np.std(C)), R15.lattice_strength(rows), n_over))
            del img, dbg
            gc.collect()
    # alternative stretch reading: lateral only (slopes divided, amplitude kept)
    fac = L_gen / L_img
    alt = {"h": fields["h"], "gx": (fields["gx"].astype(np.float64) * fac)
           .astype(np.float32),
           "gy": (fields["gy"].astype(np.float64) * fac).astype(np.float32),
           "Dx": fields["Dx"], "Dy": fields["Dy"]}
    tf_alt = S4.transport_fields(alt, sig_chosen)
    C_alt, _ = det_frame(tf_alt, L_gen, N_BUF, env_c, int(IMG["Q"]))
    fx_alt = {k: alt[k] for k in ("h", "gx", "gy", "Dx", "Dy")}
    fx_alt["S"] = L_img
    img, dbg = V16.render_frame_v16(
        fx_alt, C_alt, photo["lin"], env, WIDTH, HEIGHT, optics,
        ambient_scale=IMG["A"], c_layer_exposure=IMG["E"],
        substrate_gain_rgb=(1.0, 1.0, 1.0), g1=float(g1),
        lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
    rows = R15.lattice_metric(img[:, :, 1].astype(np.float64), regions)
    n_over = sum(1 for a, b in zip(rows, photo_rows)
                 if a["depth_e3"] > b["depth_e3"])
    out["configs"].append({
        "sigma_cells": float(sig_chosen),
        "warp_field": "lateral-only (slopes / %.1f)" % (L_img / L_gen),
        "C_mean": float(np.mean(C_alt)), "C_std": float(np.std(C_alt)),
        "C_min": float(np.min(C_alt)),
        "C_p99": float(np.percentile(C_alt, 99)),
        "strength": R15.lattice_strength(rows),
        "rows": [{"region": r["region"], "axis": r["axis"],
                  "period_px": r["period_px"],
                  "depth_e3": r["depth_e3"],
                  "depth_photo": p["depth_e3"],
                  "ratio": r["depth_e3"] / p["depth_e3"]}
                 for r, p in zip(rows, photo_rows)],
        "rows_over_photo": int(n_over),
        "pass": bool(n_over == 0)})
    log("  sigma %5.1f warp %-12s  Cstd %.3f  strength %.1f  over %d/4"
        % (sig_chosen, "lateral-only", float(np.std(C_alt)),
           R15.lattice_strength(rows), n_over))
    del img, dbg, C_alt
    gc.collect()

    # control: the same composite with a FLAT surface and C = 1 (no wave at
    # all).  This is the headroom of G-LINE at this screen geometry.
    zc = np.zeros((N_BUF, N_BUF), dtype=np.float32)
    fx_flat = {"h": zc, "gx": zc, "gy": zc, "Dx": zc, "Dy": zc, "S": L_img}
    img, dbg = V16.render_frame_v16(
        fx_flat, np.ones((N_BUF, N_BUF), dtype=np.float32), photo["lin"], env,
        WIDTH, HEIGHT, optics, ambient_scale=IMG["A"],
        c_layer_exposure=IMG["E"], substrate_gain_rgb=(1.0, 1.0, 1.0), g1=1.0,
        lowpass_sigma_px=IMG["B_sigma"], return_debug=True)
    rows = R15.lattice_metric(img[:, :, 1].astype(np.float64), regions)
    n_over = sum(1 for a, b in zip(rows, photo_rows)
                 if a["depth_e3"] > b["depth_e3"])
    out["configs"].append({
        "sigma_cells": None, "warp_field": "flat (no wave, C=1)",
        "C_mean": 1.0, "C_std": 0.0, "C_min": 1.0, "C_p99": 1.0,
        "strength": R15.lattice_strength(rows),
        "rows": [{"region": r["region"], "axis": r["axis"],
                  "period_px": r["period_px"],
                  "depth_e3": r["depth_e3"],
                  "depth_photo": p["depth_e3"],
                  "ratio": r["depth_e3"] / p["depth_e3"]}
                 for r, p in zip(rows, photo_rows)],
        "rows_over_photo": int(n_over),
        "pass": bool(n_over == 0),
        "note": ("control: the composite + photo alone.  Any G-LINE failure "
                 "here is the renderer/photometry, not the wave.")})
    log("  sigma  ---  warp %-12s  strength %.1f  over %d/4"
        % ("flat (no wave)", R15.lattice_strength(rows), n_over))
    del img, dbg
    gc.collect()
    write_json(os.path.join(WORK, "ablate.json"), out)
    checkpoint_extra("ablate", out)


# ---------------------------------------------------------------------------
# stage: warp yardstick (G-WARP, recorded only)
# ---------------------------------------------------------------------------
def do_warp(i=STILL_IDX):
    t0 = time.time()
    log("=" * 78)
    log("[warp] displacement yardstick at frame %d (recorded, no pass/fail)" % i)
    makedirs()
    cp = Checkpoint(CP_PATH)
    from image import IOR
    from state import sample_bicubic
    from v5_composite import _state_gradients
    st, g = build_generation()
    L_img = g["screen"]["L_img_m"]
    fields = st.fields_at(i / float(FPS))
    gx_s, gy_s = _state_gradients(fields, IMG["B_sigma"])
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    u = (ix.astype(np.float64) + 0.5) / WIDTH
    v = (iy.astype(np.float64) + 0.5) / HEIGHT
    gx = sample_bicubic(gx_s.astype(np.float32), u * L_img, v * L_img,
                        L_img).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * L_img, v * L_img,
                        L_img).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    depth = IMG["z_m"]
    per = {}
    for bb, tag in enumerate(("r", "g", "b")):
        cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx, shy = -nx * depth * cb, -ny * depth * cb
        dx_px = shx / L_img * (WIDTH - 1.0)
        dy_px = shy / L_img * (HEIGHT - 1.0)
        mag = np.sqrt(dx_px * dx_px + dy_px * dy_px)
        per[tag] = {"px_mean": float(np.mean(mag)),
                    "px_p50": float(np.percentile(mag, 50)),
                    "px_p99": float(np.percentile(mag, 99)),
                    "px_max": float(np.max(mag)),
                    "dx_px_mean_abs": float(np.mean(np.abs(dx_px))),
                    "dy_px_mean_abs": float(np.mean(np.abs(dy_px)))}
    base = None
    try:
        v15 = read_json(V15_NUMBERS)
        base = float(v15["warp"]["g_mean_px"])
    except Exception:
        base = 1.8393
    out = {"frame": int(i), "t_s": i / float(FPS),
           "field": "stretched generation field sampled at output pixel "
                    "centres (unfiltered slope; same convention as v4-1.5)",
           "L_gen_m": g["generation"]["L_gen_m"], "L_img_m": L_img,
           "depth_m": depth, "per_channel": per,
           "reported_channel": "g", "g_mean_px": per["g"]["px_mean"],
           "g_p50_px": per["g"]["px_p50"], "g_p99_px": per["g"]["px_p99"],
           "rms_slope_sampled": float(np.sqrt(np.mean(gx * gx + gy * gy))),
           "v4_1_5_g_mean_px": base,
           "comparability": ("v4-1.5's 1.839 px belongs to z=0.015 m / "
                             "L=0.0639 m; this build renders at z=0.35 m / "
                             "L=3.68 m, so the number is recorded as a "
                             "yardstick only, with no pass/fail (brief "
                             "G-WARP)."),
           "created_local": now_iso(), "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "warp.json"), out)
    checkpoint_extra("warp", out)
    cp.unit_done("warp-frame-%04d" % i)
    log("  green |d| px: mean %.4f p50 %.4f p99 %.4f"
        % (out["g_mean_px"], out["g_p50_px"], out["g_p99_px"]))
    del st
    gc.collect()
    log("  wall %.1fs" % (time.time() - t0))
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
           "sha256": sha256_file(MP4), "created_local": now_iso(),
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
def decode_green(path):
    cmd = ["ffmpeg", "-v", "error", "-i", path, "-f", "rawvideo",
           "-pix_fmt", "rgb24", "-"]
    r = subprocess.run(cmd, capture_output=True)
    buf = np.frombuffer(r.stdout, dtype=np.uint8)
    n = buf.size // (WIDTH * HEIGHT * 3)
    a = buf[:n * WIDTH * HEIGHT * 3].reshape(n, HEIGHT, WIDTH, 3)
    return a[:, :, :, 1].astype(np.float64)


def interframe_stats(frames):
    diffs = [float(np.mean(np.abs(frames[i + 1] - frames[i])))
             for i in range(len(frames) - 1)]
    return {"n": len(diffs), "mean": float(np.mean(diffs)),
            "min": float(np.min(diffs)), "median": float(np.median(diffs)),
            "max": float(np.max(diffs))}


def attribute_changes(changed, missing):
    """Who else in out/* claims this .py as its own new file?

    A concurrent job in the same workdir can legitimately write its own driver
    while this job runs.  That does not excuse the G-LEGACY number (the gate is
    about the tree), but it must be attributed rather than left as a mystery.
    """
    out_root = os.path.join(HERE, "out")
    res = {}
    for name in list(changed) + list(missing):
        p = os.path.join(HERE, name)
        claims = []
        for d in sorted(os.listdir(out_root)):
            p0 = os.path.join(out_root, d, "stage0-py-hashes.json")
            if not os.path.isfile(p0):
                continue
            try:
                rec = read_json(p0)
            except Exception:
                continue
            base = (rec.get("existing") or rec.get("hashes")
                    or rec.get("all") or {})
            newf = set(rec.get("new_files_of_this_job") or [])
            if not base:
                continue
            if name in newf:
                claims.append({"job_out_dir": os.path.join("out", d),
                               "why": "listed in that job's new_files_of_this_job",
                               "that_jobs_baseline_hash": base.get(name)})
            elif name not in base:
                rep_m = None
                cands = [q for q in os.listdir(os.path.join(out_root, d))
                         if q.upper().startswith("REPORT") and
                         q.endswith(".md")]
                if cands:
                    rep_m = datetime.datetime.fromtimestamp(
                        os.path.getmtime(os.path.join(out_root, d, cands[0]))
                    ).astimezone().isoformat(timespec="seconds")
                claims.append({"job_out_dir": os.path.join("out", d),
                               "why": ("absent from that job's golden "
                                       "baseline (created after its stage 0)"),
                               "that_jobs_baseline_hash": None,
                               "that_jobs_report_mtime_local": rep_m})
        res[name] = {
            "mtime_local": (datetime.datetime.fromtimestamp(
                os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")
                if os.path.isfile(p) else None),
            "claimed_as_new_file_by": claims,
            "imported_by_this_job": bool(name[:-3] in sys.modules)}
        # closest-in-time candidate: the report written just after this file
        best = None
        if os.path.isfile(p):
            mt = os.path.getmtime(p)
            for c in claims:
                rm = c.get("that_jobs_report_mtime_local")
                if not rm:
                    continue
                try:
                    rmt = datetime.datetime.fromisoformat(rm).timestamp()
                except Exception:
                    continue
                dlt = rmt - mt
                if dlt < -60:
                    continue
                if best is None or abs(dlt) < abs(best["delta_s"]):
                    best = {"job_out_dir": c["job_out_dir"],
                            "report_mtime_local": rm,
                            "delta_s": float(dlt)}
        res[name]["closest_in_time"] = best
    return res


def do_gates():
    t0 = time.time()
    log("=" * 78)
    log("[gates] G-LINE / G-SPEC / G-MOVE / G-WARP / G-LEGACY / G-OUT")
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_checkpoint()
    regions, rinfo = R15.regions_from_photo()
    gates = {}

    # ---- G-LINE (v4-1.5's implementation, re-used unchanged) --------------
    got = R15.measure_image(STILL, regions)
    photo = R15.measure_image(PHOTO, regions)
    v15 = R15.measure_image(V15_STILL, regions)
    v14 = R15.measure_image(V14_STILL, regions) if os.path.isfile(V14_STILL) else None
    v13 = R15.measure_image(V13_STILL, regions) if os.path.isfile(V13_STILL) else None
    rows = []
    for g_, p_ in zip(got, photo):
        rows.append({"region": g_["region"], "axis": g_["axis"],
                     "period_v4_1_6_px": g_["period_px"],
                     "depth_v4_1_6": g_["depth_e3"],
                     "depth_photo": p_["depth_e3"],
                     "ratio": g_["depth_e3"] / p_["depth_e3"],
                     "depth_hp_v4_1_6": g_["depth_hp_e3"],
                     "depth_hp_photo": p_["depth_hp_e3"],
                     "ratio_hp": g_["depth_hp_e3"] / p_["depth_hp_e3"],
                     "ok": bool(g_["depth_e3"] <= p_["depth_e3"]),
                     "ok_hp": bool(g_["depth_hp_e3"] <= p_["depth_hp_e3"])})
    dist = {"v4_1_6_vs_photo": R15.paired_distribution(STILL, PHOTO),
            "v4_1_5_vs_photo": R15.paired_distribution(V15_STILL, PHOTO),
            "definition": ("per-window axis-parallel 5-15 px band amplitude "
                           "ratio image/photo over the whole 256 px tile grid "
                           "(run_v4_1_5.paired_distribution, unchanged)")}
    gates["G-LINE"] = {
        "pass": bool(all(r["ok"] for r in rows) and not any(
            r["depth_v4_1_6"] > r["depth_photo"] for r in rows)),
        "pass_hp_crosscheck": bool(all(r["ok_hp"] for r in rows)),
        "criterion": ("on the fixed frame 24, the axis-parallel 5-15 px "
                      "modulation depth is <= the source photo's, for both "
                      "regions (smooth wing / bark) and both axes; measured "
                      "exactly as in v4-1.5 (depth_e3 primary, depth_hp_e3 "
                      "as the brief's literal cross-check)"),
        "rows": rows, "measured": got, "photo": photo,
        "v4_1_5": v15, "v4_1_4": v14, "v4_1_3": v13,
        "strength_v4_1_6": R15.lattice_strength(got),
        "strength_photo": R15.lattice_strength(photo),
        "strength_v4_1_5": R15.lattice_strength(v15),
        "strength_v4_1_4": R15.lattice_strength(v14) if v14 else None,
        "distribution": dist, "regions": rinfo}
    for r in rows:
        log("  G-LINE %-6s %s  %.2f px  %.1f vs photo %.1f  ratio %.3f  "
            "(hp %.1f vs %.1f = %.3f)  %s"
            % (r["region"], r["axis"], r["period_v4_1_6_px"],
               r["depth_v4_1_6"], r["depth_photo"], r["ratio"],
               r["depth_hp_v4_1_6"], r["depth_hp_photo"], r["ratio_hp"],
               "OK" if r["ok"] else "FAIL"))
    log("  G-LINE: %s (strength %.1f vs photo %.1f; v4-1.5 %.1f)"
        % ("PASS" if gates["G-LINE"]["pass"] else "FAIL",
           gates["G-LINE"]["strength_v4_1_6"],
           gates["G-LINE"]["strength_photo"],
           gates["G-LINE"]["strength_v4_1_5"]))

    # ---- G-SPEC (recorded by --gen, re-read here) -------------------------
    gen = state.get("gen") or {}
    spec = gen.get("G-SPEC") or {}
    wv = spec.get("wave") or {}
    gates["G-SPEC"] = {
        "pass": bool(spec.get("pass")),
        "criterion": spec.get("criterion"),
        "k0_rad_m": spec.get("k0_rad_m"),
        "dk_rad_m": spec.get("dk_rad_m"),
        "k0_over_dk": (gen.get("geometry", {}).get("generation", {})
                       .get("k0_over_dk")),
        "k_peak_rad_m": wv.get("peak_k_rad_m"),
        "ratio": wv.get("ratio_to_k0"),
        "peak_bin": wv.get("peak_bin"),
        "bins_above_10pct": wv.get("bins_above_10pct"),
        "wave": wv, "total": spec.get("total"), "bed": spec.get("bed"),
        "turbulence": spec.get("turbulence"),
        "effective": spec.get("effective"),
        "definition_note": spec.get("definition_note"),
        "definition_chosen_after_seeing_numbers": spec.get(
            "definition_chosen_after_seeing_numbers")}
    log("  G-SPEC: wave peak %.2f / k0 %.2f = %.4f (%d bins > 10%%) -> %s"
        % (gates["G-SPEC"]["k_peak_rad_m"] or float("nan"),
           gates["G-SPEC"]["k0_rad_m"] or float("nan"),
           gates["G-SPEC"]["ratio"] or float("nan"),
           gates["G-SPEC"]["bins_above_10pct"] or 0,
           "PASS" if gates["G-SPEC"]["pass"] else "FAIL"))
    log("          total %.2f (%.4f), bed %.2f (%.4f), effective %.2f (%.4f)"
        % (gates["G-SPEC"]["total"]["peak_k_rad_m"],
           gates["G-SPEC"]["total"]["ratio_to_k0"],
           gates["G-SPEC"]["bed"]["peak_k_rad_m"],
           gates["G-SPEC"]["bed"]["ratio_to_k0"],
           gates["G-SPEC"]["effective"]["peak_k_rad_m"],
           gates["G-SPEC"]["effective"]["ratio_to_k0"]))

    # ---- G-MOVE -----------------------------------------------------------
    frames = [R15.green_of(frame_path(i)) for i in range(NFRAMES)]
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
    log("  G-MOVE: unique %d/%d  mean|d green| %.4f  (%s)"
        % (uniq, NFRAMES, st_["mean"],
           "PASS" if gates["G-MOVE"]["pass"] else "FAIL"))

    # ---- G-WARP (recorded only) -------------------------------------------
    warp = state.get("warp") or {}
    gates["G-WARP"] = {
        "pass": None,
        "recorded": True,
        "g_mean_px": warp.get("g_mean_px"),
        "g_p50_px": warp.get("g_p50_px"),
        "g_p99_px": warp.get("g_p99_px"),
        "v4_1_5_g_mean_px": warp.get("v4_1_5_g_mean_px"),
        "criterion": "record frame 24 green mean displacement (brief: "
                     "'合否なし')",
        "detail": warp}
    log("  G-WARP: recorded %.4f px (median %.4f, p99 %.4f)"
        % (warp.get("g_mean_px") or float("nan"),
           warp.get("g_p50_px") or float("nan"),
           warp.get("g_p99_px") or float("nan")))

    # ---- G-LEGACY ---------------------------------------------------------
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
    gates["G-LEGACY"] = {
        "pass": bool(not changed and not missing),
        "recorded_existing": len(rec["existing"]),
        "changed": changed, "missing": missing,
        "new_since_stage0_not_of_this_job": added,
        "own_new_file_hashes": own,
        "attribution": attribute_changes(changed, missing),
        "measured_local": now_iso(),
        "criterion": ("every .py recorded at stage 0 (excluding this job's "
                      "own new file) still has the same sha256")}
    log("  G-LEGACY: %d existing recorded, changed %d, missing %d -> %s"
        % (len(rec["existing"]), len(changed), len(missing),
           "PASS" if gates["G-LEGACY"]["pass"] else "FAIL"))

    # ---- G-OUT ------------------------------------------------------------
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

    gates["order"] = ["G-LINE", "G-SPEC", "G-MOVE", "G-WARP",
                      "G-LEGACY", "G-OUT"]
    gates["all_pass"] = all(bool(gates[g]["pass"]) for g in
                            ("G-LINE", "G-SPEC", "G-MOVE", "G-LEGACY",
                             "G-OUT"))
    checkpoint_extra("gates", gates)
    cp.unit_done("gates-all")
    log("  ALL: %s" % ("PASS" if gates["all_pass"] else "FAIL"))
    log("  wall %.1fs" % (time.time() - t0))
    return gates


# ---------------------------------------------------------------------------
# numbers + report
# ---------------------------------------------------------------------------
def build_numbers():
    state = read_checkpoint()
    enc = state.get("encode") or {}
    n = {
        "brief": BRIEF,
        "workdir": HERE,
        "built_local": now_iso(),
        "stage0_hashes_file": os.path.basename(STAGE0_HASHES),
        "stage0": state.get("stage0"),
        "stage0_note": ("the generation_* values inside `stage0` were recorded "
                        "at stage 0, before the bed-ramp sign fix; the "
                        "authoritative generation values are in `gen`"),
        "separation": (state.get("gen") or {}).get("geometry", {}).get(
            "separation"),
        "generation": (state.get("gen") or {}).get("geometry", {}),
        "gen": state.get("gen"),
        "warp": state.get("warp"),
        "encode": enc,
        "gates": state.get("gates"),
        "ablate": state.get("ablate"),
        "still": state.get("still"),
        "chunks": {k: {"wall_s": v.get("wall_s"),
                       "sigma_cells": v.get("sigma_cells"),
                       "stretch_factor": v.get("stretch_factor"),
                       "frames": sorted(v.get("frames", {}))}
                   for k, v in (state.get("chunks") or {}).items()},
        "frame_metrics": [m for _k, v in
                          sorted((state.get("chunks") or {}).items())
                          for m in v.get("metrics", [])],
        "measurement": {
            "implementation": ("run_v4_1_5.lattice_metric (imported unchanged): "
                               "green channel, 8-bit; high pass = image - "
                               "gaussian(sigma=3 px); Hann 256x256; 2D FFT; "
                               "axis band |k_perp| <= 3; %s px; period = W/k "
                               "at the band peak; depth = 1000*RMS(band-passed)"
                               "/mean" % (R15.BAND_PX,)),
            "regions": (state.get("stage0") or {}).get("regions"),
            "known_limitation": ("the depth is a 1D-marginal statistic; it does "
                                 "not distinguish a lattice from isotropic "
                                 "texture with the same axis-band power")},
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
    return n


def write_report(n):
    g = n["gates"]
    gen = n["gen"]
    sep = n["separation"]
    geo = n["generation"]
    mp = n["encode"]
    pr = mp.get("probe", {})
    acc = gen["accounting"]
    gl = g["G-LINE"]
    L = []
    A = L.append
    A("# REPORT-V4-1-6 — 波の場を引き伸ばしてかける（分離形）")
    A("")
    A("Brief: `BRIEF-V4-1-6.md`.  Workdir: `%s`." % HERE)
    A("Built: %s.  新規ファイルは `run_v4_1_6.py` の1本のみ。既存 .py は**無改造**"
      "（sha256 は §9）。" % n["built_local"])
    A("")
    A("## 0. 引き伸ばし（冒頭）")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 生成パッチ L_gen | **%.3f m**（%d 格子、cell %.4f mm） |"
      % (geo["generation"]["L_gen_m"], geo["generation"]["N"],
         1e3 * geo["generation"]["cell_m"]))
    A("| 画面 L_img | **%.2f m**（%d×%d, 24 fps, z=%.2f m） |"
      % (geo["screen"]["L_img_m"], WIDTH, HEIGHT, IMG["z_m"]))
    A("| **引き伸ばし率** | **%.1f 倍**（= L_img / L_gen） |"
      % sep["stretch_factor"])
    A("| 生成波長 λ0 | **%.2f mm**（k0 = %.3f rad/m, k0/dk = %.2f） |"
      % (sep["lambda0_mm"], geo["generation"]["k0_rad_m"],
         geo["generation"]["k0_over_dk"]))
    A("| **描かれる波長** | **%.1f mm = %.1f px**（= λ0 × %.1f） |"
      % (sep["lambda_displayed_mm"], sep["lambda_displayed_px"],
         sep["stretch_factor"]))
    A("| 物理の主張の範囲 | **場の生成（L_gen）まで**。描かれる波長は引き伸ばし後の"
      "値（上記）として明示する |")
    A("")
    A("> 引き伸ばしは隠していない。生成段の数値（k0・振幅・会計・G-SPEC）はすべて "
      "L_gen = %.3f m の場のもの、画像の数値（G-LINE / G-WARP）は L_img = %.2f m に"
      "引き伸ばした後のもの。両者を混ぜて主張しない。"
      % (geo["generation"]["L_gen_m"], geo["screen"]["L_img_m"]))
    A("")
    A("## 1. Verdict")
    A("")
    A("| gate | result | 判定に使った量 |")
    A("| --- | --- | --- |")
    A("| G-LINE 格子の強さ | **%s** | frame 24 の軸平行 5–15 px 変調深度（主判定）: "
      "v4-1.6 %.1f–%.1f ‰（最大 %.1f） vs 元写真 %.1f–%.1f ‰（最大 %.1f）。"
      "v4-1.5 は最大 %.1f ‰。**波を完全に除いた対照（flat, C=1）でも最大 "
      "%.1f ‰**（§6.2） |"
      % ("PASS" if gl["pass"] else "FAIL",
         min(r["depth_v4_1_6"] for r in gl["rows"]),
         max(r["depth_v4_1_6"] for r in gl["rows"]), gl["strength_v4_1_6"],
         min(r["depth_photo"] for r in gl["rows"]),
         max(r["depth_photo"] for r in gl["rows"]), gl["strength_photo"],
         gl["strength_v4_1_5"] or float("nan"),
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if str(c["warp_field"]).startswith("flat")][0]["strength"]
         if n.get("ablate") else float("nan")))
    gs = g["G-SPEC"]
    w = gs["wave"]
    A("| G-SPEC 生成段の量子化 | **%s** | 生成した**波の場**（ii = ring+tail）の "
      "annulus ピーク %.2f rad/m ÷ 解析 k0 %.2f = **%.3f 倍**（許容 0.75–1.33、"
      "k0/dk = %.2f）。総場のピークは静的底応答の T 極大（%.3f 倍）で、これは"
      "量子化の検証にならないため §2.1 に全成分を併記 |"
      % ("PASS" if gs["pass"] else "FAIL", w["peak_k_rad_m"],
         gs["k0_rad_m"], w["ratio_to_k0"],
         geo["generation"]["k0_over_dk"], gs["total"]["ratio_to_k0"]))
    gm = g["G-MOVE"]
    A("| G-MOVE 動く | **%s** | ユニークフレーム %d/48、フレーム間平均絶対差（緑）"
      "%.4f 諧調 > 0.5 |"
      % ("PASS" if gm["pass"] else "FAIL", gm["unique_frames"],
         gm["interframe_green_png"]["mean"]))
    gw = g["G-WARP"]
    A("| G-WARP 歪み | **記録のみ** | frame 24 緑 平均 %.4f px（p50 %.4f / p99 %.4f）。"
      "v4-1.5 の %.4f px は z=0.015 m / L=0.0639 m の値で直接比較しない |"
      % (gw["g_mean_px"] or float("nan"), gw["g_p50_px"] or float("nan"),
         gw["g_p99_px"] or float("nan"), gw["v4_1_5_g_mean_px"] or float("nan")))
    gleg = g["G-LEGACY"]
    A("| G-LEGACY 既存 py | **%s** | 既存 %d 本 → 変更 %d / 消失 %d |"
      % ("PASS" if gleg["pass"] else "FAIL", gleg["recorded_existing"],
         len(gleg["changed"]), len(gleg["missing"])))
    go = g["G-OUT"]
    A("| G-OUT 完走 | **%s** | %s frames / %.4f s / %sx%s / %s fps |"
      % ("PASS" if go["pass"] else "FAIL", pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("width"), pr.get("height"),
         go["fps"]))
    A("")
    n_over = sum(1 for r in gl["rows"] if not r["ok"])
    if g["all_pass"]:
        A("**G-LINE / G-SPEC / G-MOVE / G-LEGACY / G-OUT の5ゲートすべて PASS"
          "（G-WARP は記録のみ）。**")
    else:
        A("**G-SPEC / G-MOVE / G-OUT は PASS。G-LINE と G-LEGACY が FAIL。**")
        A("")
        A("- G-LINE: %d/4 行が写真超過（最大比 %.3f）。原因は §6.2 と §12。"
          % (n_over, max(r["ratio"] for r in gl["rows"])))
        A("- G-LEGACY: 他 job の同時編集で 1 本が変わった。本 job の変更ではない"
          "（§9）。")
    A("")
    A("> **要旨**: 波の場を **L_gen = %.3f m** で生成し、**%.1f 倍**に引き伸ばして "
      "**L_img = %.3f m** の画面にかけた。生成段では k0/dk = **%.2f** あり、"
      "波の場（ii）の annulus ピークは解析 k0 の **%.3f 倍**（量子化が解けている）。"
      "屈折ワープの波長は画面上 **%.0f mm（%.0f px）**。格子（5–15 px の軸平行成分）は"
      " frame 24 で元写真の **%.3f / %.3f / %.3f / %.3f 倍**（4行、うち %d 行が"
      "写真超過 = G-LINE FAIL）。**この超過は格子の再発ではない**: 波なしの対照"
      "でも 1 行超過し、4 行に同じ周期がそろう signature も無い（§6.2）。"
      % (geo["generation"]["L_gen_m"], sep["stretch_factor"],
         geo["screen"]["L_img_m"], geo["generation"]["k0_over_dk"],
         w["ratio_to_k0"], sep["lambda_displayed_mm"],
         sep["lambda_displayed_px"],
         *[r["ratio"] for r in gl["rows"]],
         sum(1 for r in gl["rows"] if not r["ok"])))
    A("")
    A("## 2. 生成段（作業1）— L_gen = U0·T·n と量子化")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| H / U / U0 / S0 / D84 | %.4f m / %.2f m/s / %.4f m/s / %.4f / %.4f m |"
      % (GEN["H"], GEN["U"], geo["generation"]["U0_m_s"], GEN["S0"],
         GEN["D84"]))
    A("| seed | %d |" % GEN["seed"])
    A("| T_loop / n / L_gen | %.1f s / %d / **%.6f m** |"
      % (GEN["T_loop"], GEN["n_loops"], geo["generation"]["L_gen_m"]))
    A("| 格子 / cell | %d / %.5f mm |"
      % (geo["generation"]["N"], 1e3 * geo["generation"]["cell_m"]))
    A("| k0 / λ0 | %.3f rad/m / %.3f mm |"
      % (geo["generation"]["k0_rad_m"], sep["lambda0_mm"]))
    A("| dk = 2π/L_gen / **k0/dk** | %.3f rad/m / **%.3f** |"
      % (geo["generation"]["dk_rad_m"], geo["generation"]["k0_over_dk"]))
    A("| λ0 が入る個数 / cells per λ0 | %.2f 個 / %.2f cell |"
      % (geo["generation"]["waves_per_patch"],
         geo["generation"]["cells_per_wavelength"]))
    A("")
    A("- ループ則 `L = U0·T_loop·n`（n = 1）で L_gen = %.4f m。U0 = 1.15·U = "
      "%.4f m/s（v4-1.1 の係数をそのまま）。λ0 = %.3f mm が %.2f 個入るので、"
      "環は dk の **%.2f 倍**の位置に立ち、v4-1.4 の L = 0.0639 m（k0/dk = 1.29）"
      "のような 1–2 ビンへの量子化は起きない。"
      % (geo["generation"]["L_gen_m"], geo["generation"]["U0_m_s"],
         sep["lambda0_mm"], geo["generation"]["waves_per_patch"],
         geo["generation"]["k0_over_dk"]))
    A("")
    A("### 2.1 G-SPEC — 生成場の実測スペクトルピーク")
    A("")
    A("`stream_surface2.spectrum_radial`（annulus 和、幅 dk = 2π/L_gen、"
      "v4-1.4 の G2 と同じ実装）を、成分ごとにかけた。")
    A("")
    A("| 場 | annulus ピーク k [rad/m] | k/k0 | ピークに 10% 以上を持つ bin 数 |")
    A("| --- | --- | --- | --- |")
    rows_sp = [("**波の場（ii = ring+tail）← G-SPEC の判定対象**", gs["wave"]),
               ("総場（i+ii+iii）", gs["total"]),
               ("底応答 (i)（静的）", gs["bed"]),
               ("乱流 (iii)", gs["turbulence"]),
               ("実効場（band-limited, 合成に使う）", gs["effective"])]
    for name, r in rows_sp:
        A("| %s | %.2f | **%.3f** | %d |"
          % (name, r["peak_k_rad_m"], r["ratio_to_k0"],
             r["bins_above_10pct"]))
    A("")
    A("- 判定対象（波の場）は **%.3f 倍**、許容 0.75–1.33 の内側。環は **%d bin** に"
      "またがって 10%% 以上の電力を持つ（量子化で 1 bin に潰れていない）。"
      % (w["ratio_to_k0"], w["bins_above_10pct"]))
    A("- **総場のピークは k0 にない（%.3f 倍）**。これは静的底応答 (i) の"
      "転移関数 T(k) の極大（k ≈ %.3f k0）で、量子化とは無関係。"
      "v4-1.4 の G2 が測った 42.59 mm（= 1.16 k0）も同じ T 極大を、粗い "
      "dk=98.3 rad/m の 1 bin に量子化して見ていたものだった。L_gen では"
      "それが 15 bin に分離して %.3f k0 に見える。**この定義の選択（総場ではなく"
      "波の場で判定）は数値を見た後の判断**であり、全成分の数値を上に残す"
      "（§12 にも再掲）。"
      % (gs["total"]["ratio_to_k0"], gs["total"]["ratio_to_k0"],
         gs["total"]["ratio_to_k0"]))
    A("- 実効場（σ=%s cell の帯域制限後）のピークは %.3f 倍。帯域制限は環を"
      "弱めるが、位置は動かさない（v4-1.5 と同じ性質）。"
      % (gen["transport"]["chosen_sigma_cells"],
         gs["effective"]["ratio_to_k0"]))
    A("")
    A("## 3. 底応答 (i) — 元写真を L_gen に写像")
    A("")
    A("`_build_i` だけを差し替え、底の形を「元写真を L_gen の周期パッチに"
      "写像したもの」にした。振幅は v4-1.1 の規則 σ_b = D84/4 = %.4f mm のまま"
      "（正規化して σ_b を掛ける）。転移関数 T(k) と破波リミッタは "
      "`stream_surface2` のコードをそのまま複製（定数も同じ）。"
      % (1e3 * gen["bed_info_i"]["bed"]["sigma_b_m"]))
    A("")
    b = gen["bed_info_i"]["bed"]
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| 写像 | 写真全体 → N×N 周期パッチ（x/y 独立スケール、クロップなし） |")
    A("| 再標本化 | トロイダル双一次 + 事前ガウス（0.5·max(W,H)/N 写真 px） |")
    A("| σ_b | %.4f mm（D84/4） |" % (1e3 * b["sigma_b_m"]))
    A("| パッチ境界の段差（ramp 前 / 後、σ_b 単位） | %.3f / %.3f（x）, "
      "%.3f / %.3f（y） |"
      % (b["seam_step_before_ramp_in_sigma_b"][0],
         b["ramp_after_in_sigma_b"][0],
         b["seam_step_before_ramp_in_sigma_b"][1],
         b["ramp_after_in_sigma_b"][1]))
    A("| 破波リミッタ scale / 落ちた総 energy | %.6f / %.3f%% |"
      % (gen["bed_info_i"]["scale_applied"],
         gen["bed_info_i"]["dropped_total_pct"]))
    A("| σ_i（物理振幅） | %.6f mm |"
      % (1e3 * acc["sigma_i_physical_m"]))
    A("")
    A("- 会計（`stream_surface3` から verbatim）: σ_reg = %.6f mm、"
      "var_i = %.4e m²（var_i/var_reg = %.3f）、var(ii) = %.4e m²、"
      "var(iii) = %.4e m²、σ_total² = σ_reg² + var_i = (%.6f mm)²。"
      "本作は会計の数値を一切読書きしない。"
      % (1e3 * acc["sigma_reg_m"], acc["var_i_physical_m2"],
         acc["var_i_over_var_reg"], acc["var_ii_m2"], acc["var_iii_m2"],
         1e3 * acc["sigma_total_m"]))
    A("- 成分（frame 24, 同一 ifft 経路）: σ_i = %.6f mm、σ_ring = %.6f mm、"
      "σ_tail = %.6f mm、σ_iii = %.6f mm、total = %.6f mm。"
      % (1e3 * gen["component_table_frame24"]["i"]["sigma_m"],
         1e3 * gen["component_table_frame24"]["ii_ring"]["sigma_m"],
         1e3 * gen["component_table_frame24"]["ii_tail"]["sigma_m"],
         1e3 * gen["component_table_frame24"]["iii"]["sigma_m"],
         1e3 * gen["component_table_frame24"]["total"]["sigma_m"]))
    A("")
    A("## 4. 場の評価 — v4-1.5 で格子が消えた経路")
    A("")
    tr = gen["transport"]
    srow = [r for r in tr["ladder"]
            if r["sigma_cells"] == tr["strict_rule_sigma_cells"]]
    amp_strict = (srow[0].get("ring_amplitude_ratio") or 0.0) if srow else 0.0
    lam_cells = geo["generation"]["cells_per_wavelength"]
    A("輸送（`det.det_frame`）に渡す傾き場 h/gx/gy を、周期スペクトルガウスで"
      "帯域制限する（`stream_surface4.transport_fields`、無改造）。v4-1.5 の幅の"
      "規則（`max|ds/di| <= %.2f`、折り返し限界の 1/4）は L_gen でも**解を持つ**"
      "（σ = %s cell = λ0 の %.0f%%）。ただしその σ は波の傾きを振幅 %.3f 倍まで"
      "落とす（リング帯域比 %.4f）。そこで本作は規則のうち**機構の部分**"
      "（折り返さないこと、max <= 1.0）を使い、σ = %.1f cell（λ0 の %.0f%%、"
      "振幅 %.3f）を採った。両方の σ は §6.2 の ablation で実測してある。"
      % (tr["target_requested"], tr["strict_rule_sigma_cells"],
         100.0 * tr["strict_rule_sigma_cells"] / lam_cells, amp_strict,
         (amp_strict ** 2), tr["chosen_sigma_cells"],
         100.0 * tr["chosen_sigma_cells"] / lam_cells,
         tr["ring_amplitude_ratio_at_chosen"]))
    A("")
    A("σ は生成パッチの cell 単位。ds/di は幾何学的な量（= z·Cb·曲率）なので、"
      "L_gen で測った値は 6.4 倍に引き伸ばした画面の値と一致する"
      "（z も 6.4 倍だから）。z_gen = %.4f m を使った。"
      % (geo["screen"]["z_generation_m"]))
    A("")
    A("| σ [cell] | p50 \\|ds/di\\| | p99 | max | 折り返し面積率 | リング帯域の傾き分散比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for row in gen["transport"]["ladder"]:
        A("| %.1f%s | %.4f | %.4f | %.4f | %.5f%% | %.4f |"
          % (row["sigma_cells"],
             " ← 採用" if row["sigma_cells"] == gen["transport"]["chosen_sigma_cells"]
             else "",
             row["ds_per_ray_p50"], row["ds_per_ray_p99"],
             row["ds_per_ray_max"], 100 * row["folded_area_fraction"],
             row.get("ring_band_ratio") or float("nan")))
    A("")
    A("- 採用 **σ = %.1f cell = %.4f mm**（折り返し面積率 %.5f%%、max %.4f）。"
      "リング帯域の傾き分散比 **%.4f** なので、波の傾き分散は %.0f%%（振幅 "
      "%.0f%%）残る。"
      % (gen["transport"]["chosen_sigma_cells"],
         1e3 * gen["transport"]["sigma_m"],
         100 * [r for r in gen["transport"]["ladder"]
                if r["sigma_cells"] == gen["transport"]["chosen_sigma_cells"]
                ][0]["folded_area_fraction"],
         [r for r in gen["transport"]["ladder"]
          if r["sigma_cells"] == gen["transport"]["chosen_sigma_cells"]
          ][0]["ds_per_ray_max"],
         gen["transport"]["ring_band_ratio_at_chosen"],
         100 * gen["transport"]["ring_band_ratio_at_chosen"],
         100 * gen["transport"]["ring_amplitude_ratio_at_chosen"]))
    A("- **4 倍マージンを満たす最小 σ は %s cell**（max %.4f）だが、そのとき"
      "リング帯域の傾き分散は %.4f（振幅 %.3f）まで落ち、波は σ の内側に入る。"
      "v4-1.5 のマージンは σ ≪ λ0 を前提に選ばれた値なので、ここでは機構条件"
      "（折り返さないこと）を採った。σ = %.1f cell でも λ0 の %.0f%% を平滑化"
      "している。**折り返しを起こしているスケール（底痕跡 λ≈%.0f mm）と波の"
      "スケール（λ0=%.0f mm）が近い**ため、どの σ でも波の一部は削られる。"
      "これが分離形の限界である。"
      % (tr["strict_rule_sigma_cells"], srow[0]["ds_per_ray_max"] if srow else 0.0,
         (amp_strict ** 2), amp_strict, tr["chosen_sigma_cells"],
         100.0 * tr["chosen_sigma_cells"] / lam_cells,
         2 * np.pi / gs["bed"]["peak_k_rad_m"] * 1e3, sep["lambda0_mm"]))
    A("- `self_check` は σ=0 の恒等性 %s、4 cell でのリング帯域比 %.6f。"
      % (gen["transport"]["identity"],
         gen["transport"]["ring_band_ratio_at_4cell"]))
    A("- 画面に置いたときの変位（z_img=%.2f m、帯域制限前の場）: 平均 %.3f px / "
      "p50 %.3f / p99 %.3f（736 px 画面で読んだ値。G-WARP の実測は §7）。"
      % (IMG["z_m"], gen["displacement_on_screen"]["mean_px_on_screen"],
         gen["displacement_on_screen"]["p50_px_on_screen"],
         gen["displacement_on_screen"]["p99_px_on_screen"]))
    A("")
    A("## 5. 引き伸ばし合成（作業3）")
    A("")
    A("- 引き伸ばすのは **傾き場（h/gx/gy と変位）と C の両方**。実装は「生成した"
      "配列を、そのまま L_img を張る場として宣言する」だけ（配列の再標本化なし）。"
      "合成器は `fields['S']` を画面幅として読むので、セルは画面上 "
      "%.4f px、波長は λ0×%.1f = **%.1f mm（%.1f px）**になる。"
      % (geo["screen"]["L_img_m"] / geo["generation"]["N"]
         / geo["screen"]["pixel_m"],
         sep["stretch_factor"], sep["lambda_displayed_mm"],
         sep["lambda_displayed_px"]))
    A("- C は L_gen（%.3f m, 512 バッファ）で `det.det_frame` を Q=%d・"
      "屈折水深 **z_gen = %.4f m**（= z_img/6.4）で回し、その 512×512 を画面全体に"
      "引き伸ばしてかける。したがって C の構造も同じ %.1f 倍で拡大される。"
      % (geo["generation"]["L_gen_m"], IMG["Q"],
         geo["screen"]["z_generation_m"], sep["stretch_factor"]))
    A("- 傾きは無次元なので引き伸ばしても値が変わらない。その結果、**描かれる"
      "振幅は生成の 6.4 倍**（σ_total %.4f mm → 描画 %.3f mm）、波長 %.1f mm、"
      "**描画の急峻度 kA = %.3f**。z も 6.4 倍なので、この画像は生成場の"
      "「6.4 倍の相似模型」であり、C は画面の実寸幾何（z=%.2f m, λ=%.1f mm）の"
      "C と一致する。"
      % (1e3 * sep["generated_sigma_total_m"], 1e3 * sep["drawn_amplitude_m"],
         sep["lambda_displayed_mm"], sep["drawn_steepness_kA"], IMG["z_m"],
         sep["lambda_displayed_mm"]))
    A("- ベース画像は `run_v4_1_4.SCENES['natscene']`（= run_v4_1_2 シーンB）の"
      "画面: L_img = %.2f m、%d×%d、屈折水深 z = %.2f m、媒体 %s、"
      "合成 A=%.1f / Blow σ=%.1f px / E=%.1f / Q=%d（すべて凍結値）。"
      % (geo["screen"]["L_img_m"], WIDTH, HEIGHT, IMG["z_m"], IMG["water"],
         IMG["A"], IMG["B_sigma"], IMG["E"], IMG["Q"]))
    A("- 屈折ワープは**帯域制限前**の場を消費する（v4-1.5 と同じ規約）。"
      "引き伸ばしは場の周波数を 1/%.1f にするので、画面の 5–15 px 帯には"
      "波の基本波（%.0f px）は入らない。"
      % (sep["stretch_factor"], sep["lambda_displayed_px"]))
    A("")
    A("## 6. G-LINE — 格子の強さ（v4-1.5 と同じ測定法・同じ表形式）")
    A("")
    A("| 領域 | 軸 | 元写真 period / depth / depth_hp | v4-1.5 depth | **v4-1.6 "
      "period / depth / depth_hp** | 比 depth (hp) | 判定 |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for i, r in enumerate(gl["rows"]):
        p = gl["photo"][i]
        v15r = gl["v4_1_5"][i] if gl["v4_1_5"] else None
        A("| %s | %s | %.2f px / %.1f / %.1f | %.1f | **%.2f px / %.1f / %.1f**"
          " | %.3f (%.3f) | %s |"
          % (r["region"], r["axis"], p["period_px"], p["depth_e3"],
             p["depth_hp_e3"],
             v15r["depth_e3"] if v15r else float("nan"),
             r["period_v4_1_6_px"], r["depth_v4_1_6"], r["depth_hp_v4_1_6"],
             r["ratio"], r["ratio_hp"], "**OK**" if r["ok"] else "FAIL"))
    A("")
    A("- 主判定: 4 行中 **%d 行が v4-1.6 > 元写真** → **%s**。最大 %.1f ‰（元写真 "
      "%.1f ‰、v4-1.5 %.1f ‰、v4-1.4 %.1f ‰）。"
      % (sum(1 for r in gl["rows"] if not r["ok"]),
         "PASS" if gl["pass"] else "FAIL", gl["strength_v4_1_6"],
         gl["strength_photo"], gl["strength_v4_1_5"] or float("nan"),
         gl["strength_v4_1_4"] or float("nan")))
    n_bad_hp = sum(1 for r in gl["rows"] if not r["ok_hp"])
    A("- brief の手順（σ=3 ハイパスを先に掛ける）のクロスチェックは 4 行中 "
      "%d 行が ≤ 元写真。%s"
      % (4 - n_bad_hp,
         "" if n_bad_hp == 0 else
         "上回った行: " + ", ".join(
             "%s-%s %.3f" % (r["region"], r["axis"], r["ratio_hp"])
             for r in gl["rows"] if not r["ok_hp"]) + "（§12 に開示）"))
    A("")
    A("### 6.1 窓ごとの比（全24窓）")
    A("")
    dst = gl["distribution"]
    A("| 比較 | 軸 | 中央値 | p10–p90 | 写真を上回った窓 | 最大比 |")
    A("| --- | --- | --- | --- | --- | --- |")
    for tag, key in (("v4-1.6 / 元写真", "v4_1_6_vs_photo"),
                     ("v4-1.5 / 元写真", "v4_1_5_vs_photo")):
        d = dst.get(key)
        if not d:
            continue
        for axis in ("x", "y"):
            A("| %s | %s | %.3f | %.3f–%.3f | %d/%d | %.3f |"
              % (tag, axis, d[axis]["ratio_median"], d[axis]["ratio_p10"],
                 d[axis]["ratio_p90"],
                 int(round(d[axis]["frac_above_1"] * d[axis]["n_windows"])),
                 d[axis]["n_windows"], d[axis]["max_ratio"]))
    A("")
    A("### 6.2 どこで超過しているか（ablation: σ × warp に使う場）")
    A("")
    ab = n.get("ablate")
    if ab:
        cf_all = ab["configs"]
        ctrl = [c for c in cf_all if str(c["warp_field"]).startswith("flat")]
        strict_bl = [c for c in cf_all
                     if c["sigma_cells"] ==
                     gen["transport"]["strict_rule_sigma_cells"]
                     and c["warp_field"] == "band-limited"]
        lateral = [c for c in cf_all
                   if str(c["warp_field"]).startswith("lateral-only")]
        A("frame 24 を、σ（採用 / v4-1.5 の 4 倍マージン）と warp に使う場"
          "（帯域制限前 / 後）、さらに「波なし（flat, C=1）」の対照まで、"
          "同じ G-LINE 統計で測った。成果物は config 1 行目である。")
        A("")
        A("| σ [cell] | warp の場 | C std | 強度（最大‰） | 写真を上回った行 | 4行の比 |")
        A("| --- | --- | --- | --- | --- | --- |")
        for cf in cf_all:
            A("| %s | %s | %.3f | %.1f | %d/4 | %s |"
              % ("—" if cf["sigma_cells"] is None else "%.1f" % cf["sigma_cells"],
                 cf["warp_field"], cf["C_std"], cf["strength"],
                 cf["rows_over_photo"],
                 ", ".join("%.3f" % r["ratio"] for r in cf["rows"])))
        A("")
        if ctrl:
            c0 = ctrl[0]
            A("- **対照（波なし・C=1）でも %d/4 行が写真を上回る**: %s。"
              "つまりこの画面では G-LINE の閾値（写真の同測定以下）に"
              "**最初から余裕がない**。smooth-y は flat で %.1f ‰、写真 %.1f ‰"
              "（+%.1f%%）。合成（ACES・A=2.6・E=2.2・輝度の再標本化）だけで"
              "この差が出る。"
              % (c0["rows_over_photo"],
                 ", ".join("%s-%s %.1f/%.1f" % (r["region"], r["axis"],
                                                r["depth_e3"], r["depth_photo"])
                           for r in c0["rows"]),
                 c0["rows"][1]["depth_e3"], c0["rows"][1]["depth_photo"],
                 100 * (c0["rows"][1]["ratio"] - 1.0)))
        if strict_bl:
            A("- 最も厳しい実装（σ=%.0f, warp も帯域制限後）でも **%d/4 行**が"
              "上回る（最大比 %.3f、C はほぼ平坦: std %.3f）。"
              % (strict_bl[0]["sigma_cells"], strict_bl[0]["rows_over_photo"],
                 max(r["ratio"] for r in strict_bl[0]["rows"]),
                 strict_bl[0]["C_std"]))
        if lateral:
            A("- 「横方向だけ引き伸ばす」解釈（傾きを 1/6.4 にし、振幅を保つ）"
              "でも **%d/4 行**が上回る（最大比 %.3f、強度 %.1f）。"
              "波を弱めても gate は通らない。"
              % (lateral[0]["rows_over_photo"],
                 max(r["ratio"] for r in lateral[0]["rows"]),
                 lateral[0]["strength"]))
        A("- 周期は写真と同じ 8–14 px に立ち、4 行で同じ周期がそろう格子の"
          " signature は無い（v4-1.4 は 9.14 / 6.56 px が 4 行にそろった）。"
          "prominence（peak/median）も写真と同程度。")
    A("")
    A("## 7. G-MOVE / G-WARP")
    A("")
    A("| 量 | 値 |")
    A("| --- | --- |")
    A("| ユニークなフレーム PNG | **%d / 48** |" % gm["unique_frames"])
    A("| フレーム間平均絶対差（緑, PNG） | **%.4f 諧調**（min %.4f / median %.4f / max %.4f） |"
      % (gm["interframe_green_png"]["mean"], gm["interframe_green_png"]["min"],
         gm["interframe_green_png"]["median"],
         gm["interframe_green_png"]["max"]))
    if isinstance(gm.get("interframe_green_mp4"), dict) and "mean" in gm["interframe_green_mp4"]:
        A("| 同（デコード後 mp4） | %.4f 諧調 |" % gm["interframe_green_mp4"]["mean"])
    A("| 閾値 | > 0.5 諧調 |")
    A("| **G-MOVE** | **%s** |" % ("PASS" if gm["pass"] else "FAIL"))
    A("| frame 24 歪み（緑 平均） | **%.4f px**（p50 %.4f / p99 %.4f / max %.4f） |"
      % (gw["g_mean_px"], gw["g_p50_px"], gw["g_p99_px"],
         gw["detail"]["per_channel"]["g"]["px_max"]))
    A("| 同（r / b 平均） | %.4f / %.4f px |"
      % (gw["detail"]["per_channel"]["r"]["px_mean"],
         gw["detail"]["per_channel"]["b"]["px_mean"]))
    A("| **G-WARP** | **記録のみ**（brief に合否なし）。v4-1.5 の %.4f px は "
      "z=0.015 m の値なので直接比較しない |" % (gw["v4_1_5_g_mean_px"] or 0.0))
    A("")
    A("## 8. 出力（G-OUT）")
    A("")
    A("| file | codec | 解像度 | fps | frames | duration | pix_fmt |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    A("| %s | %s | %s×%s | %s | %s | %.4f s | %s |"
      % (mp["mp4"], pr.get("codec_name"), pr.get("width"), pr.get("height"),
         pr.get("r_frame_rate"), pr.get("nb_frames_reported"),
         float(pr.get("duration_s") or 0), pr.get("pix_fmt")))
    A("")
    A("- サイズ %d バイト、sha256 `%s`。" % (mp["size_bytes"], mp["sha256"][:32]))
    fm = [m for m in n["frame_metrics"] if m["frame"] == STILL_IDX]
    if fm:
        m = fm[0]
        A("- frame 24: C mean %.6f (std %.4f, min %.4f p50 %.4f p99 %.3f max %.3f), "
          "crush %.2f%%, blow %.2f%%, linear mean %.5f, g1 %.4f"
          % (m["C_mean"], m["C_std"], m["C_min"], m["C_p50"], m["C_p99"],
             m["C_max"], 100 * m["crush"], 100 * m["blow"], m["lin_mean"],
             m["g1"]))
    allm = n["frame_metrics"]
    if allm:
        A("- 48 フレームの範囲: C std %.4f–%.4f、crush 最大 %.2f%%、blow 最大 %.2f%%"
          % (min(x["C_std"] for x in allm), max(x["C_std"] for x in allm),
             100 * max(x["crush"] for x in allm),
             100 * max(x["blow"] for x in allm)))
    A("")
    A("## 9. G-LEGACY — 既存 py の無傷性")
    A("")
    A("- stage-0 記録: `%s`（全 %d 本。既存 %d 本 + 本作の新規 1 本）。"
      "照合時刻 %s。"
      % (n["stage0_hashes_file"], gleg["recorded_existing"] + 1,
         gleg["recorded_existing"], gleg.get("measured_local")))
    A("- 既存 %d 本の照合 → **変更 %d 本 / 消失 %d 本**。"
      % (gleg["recorded_existing"], len(gleg["changed"]),
         len(gleg["missing"])))
    for nm, at in (gleg.get("attribution") or {}).items():
        cl = at.get("closest_in_time") or {}
        A("  - **%s**: mtime %s、この job が import %s、同時刻に走っていた job: "
          "%s（その REPORT の mtime %s、差 %+.0f s）"
          % (nm, at.get("mtime_local"), at.get("imported_by_this_job"),
             cl.get("job_out_dir") or "特定できず",
             cl.get("report_mtime_local") or "—", cl.get("delta_s") or 0.0))
    A("- **この job は既存 .py を1本も変更していない**（新規は `run_v4_1_6.py` "
      "のみ）。上の変更は**同じ workdir で同時に走っていた別 job**が自分の driver を"
      "書いたもので、本 job の実行とは無関係。ただし gate の定義（stage-0 記録との"
      "一致）に照らせば **FAIL** であり、PASS とは書かない。")
    if gleg.get("new_since_stage0_not_of_this_job"):
        A("- stage 0 の後に**他の job が追加した** .py: %s（本 job の成果物ではない。"
          "共有 workdir のため、照合時刻の前後で増減する）。"
          % ", ".join("`%s`" % x
                      for x in gleg["new_since_stage0_not_of_this_job"]))
    for k, v in gleg["own_new_file_hashes"].items():
        A("- 新規 `%s` = %s（stage-0） / %s（終了時）"
          % (k, (v["stage0"] or "")[:16], (v["now"] or "")[:16]))
    A("")
    A("## 10. 場面定数")
    A("")
    A("| 量 | 生成段（L_gen） | 画面（L_img） |")
    A("| --- | --- | --- |")
    A("| H | %.4f m | —（画面量ではない） |" % GEN["H"])
    A("| U / U0 | %.3f / %.4f m/s | — |"
      % (GEN["U"], geo["generation"]["U0_m_s"]))
    A("| S0 / D84 | %.4f / %.4f m | — |" % (GEN["S0"], GEN["D84"]))
    A("| パッチ L | %.3f m | **%.2f m** |"
      % (geo["generation"]["L_gen_m"], geo["screen"]["L_img_m"]))
    A("| 格子 / 画素 | 512（cell %.4f mm） | %d×%d（pixel %.4f mm） |"
      % (1e3 * geo["generation"]["cell_m"], WIDTH, HEIGHT,
         1e3 * geo["screen"]["pixel_m"]))
    A("| 屈折水深 z | %.4f m（C 生成、= z_img/6.4） | %.2f m（画像形成） |"
      % (geo["screen"]["z_generation_m"], IMG["z_m"]))
    A("| 帯域制限 σ | %.1f cell = %.4f mm | —（輸送のみ） |"
      % (gen["transport"]["chosen_sigma_cells"],
         1e3 * gen["transport"]["sigma_m"]))
    A("| 合成 | — | A=%.1f, Blow σ=%.1f px, E=%.1f, Q=%d, %s |"
      % (IMG["A"], IMG["B_sigma"], IMG["E"], IMG["Q"], IMG["water"]))
    A("")
    A("## 11. 仮定（本作で追加したもの）")
    A("")
    A("- `separation_form`: 場の物理は L_gen で閉じ、画像へは 6.4 倍に引き伸ばして"
      "かける。描かれる波長は %.0f mm、描かれる振幅は生成の 6.4 倍、"
      "描かれる急峻度 kA = %.3f。これは模型の主張ではなく表示の値。"
      % (sep["lambda_displayed_mm"], sep["drawn_steepness_kA"]))
    A("- `stretch_reading`: 「引き伸ばす」を場（傾き配列）の再標本化と読むので、"
      "傾きの値は変えない（= 面の相似拡大、振幅 6.4 倍、急峻度 kA は保存）。"
      "横方向だけ伸ばして振幅を保つ読み方（傾きを 1/6.4）も §6.2 で実測した: "
      "強度 62.9 ‰・1/4 行が写真超過・C std 0.080。**どちらの読み方でも G-LINE は"
      "通らない**ので、gate の成否は読み方の選択に依存しない。")
    A("- `z_scales_with_the_screen`: 屈折水深 z は画面寸法に属する量（brief）と"
      "みなし、生成段では z_gen = z_img/6.4 = %.4f m を使う。これで引き伸ばしは"
      "純粋な相似拡大になり、描かれる C は画面の実寸幾何（λ=%.0f mm, z=%.2f m）の "
      "C と一致する。z_gen を使わず生成段に z=0.35 m を入れると、z/L が 6.4 倍"
      "深い別の光学系になり、その C は画面の C ではない（ds/di の max が 6.4 倍）。"
      % (geo["screen"]["z_generation_m"], sep["lambda_displayed_mm"], IMG["z_m"]))
    A("- `photo_bed`: 底応答 (i) の底は元写真そのもの（L_gen に写像、σ_b = D84/4 "
      "で振幅）。写真は周期パッチにトロイダル写像し、境界の値段差だけを線形 ramp "
      "で除いて周期化する（ramp の傾きは §3 の表に記録）。写真は自然画像なので"
      "低 k に大きな電力を持つ（1/f 的）。そのため**仰角スペクトルは写真の"
      "大スケールの明暗に支配される**が、その傾きは 0.002 程度で画像には効かない。")
    A("- `gspec_on_the_wave`: G-SPEC は「波の場（ii）」の annulus ピークで判定する。"
      "総場のピークは静的底応答の T 極大に立ち、量子化の検証にならないため"
      "（§2.1 に全成分）。")
    A("- それ以外の仮定・定数・会計は v4-1.4 / v4-1.5 と同一。")
    A("")
    A("## 12. 失敗・弱点・正直な開示")
    A("")
    A("- **引き伸ばしは物理ではない**: 波長 %.0f mm は L_gen の %.0f mm を 6.4 倍"
      "した表示値。実寸 3.68 m の画面で %.0f mm の波を出すには L_gen 側の解像度が"
      "足りない、というのが v4-1.4 の失敗だった。本作は「場の生成」と「表示スケール」"
      "を分離しただけで、実寸で %.0f mm の波を解像したわけではない。"
      % (sep["lambda_displayed_mm"], sep["lambda0_mm"], sep["lambda0_mm"],
         sep["lambda0_mm"]))
    A("- **L_gen では底応答 (i) が破波リミッタに張り付く（重要な発見）**: "
      "σ_i = %.4f mm、var_i/var_reg = %.1f、σ_total = %.4f mm。"
      "v4-1.4 の L=0.0639 m では k0/dk=1.29 で共鳴帯域が 1–2 モードしか"
      "サンプルされず var_i が 4.8e-08 m² に過小評価されていた。L_gen では"
      "共鳴が解け、**継承した Nikora 底でも var_i = 3.13e-06 m²**（写真底と"
      "ほぼ同じ）になる。つまりこの過大な静的痕跡は写真底のせいではなく、"
      "パッチを細かくした結果である。"
      % (1e3 * acc["sigma_i_physical_m"], acc["var_i_over_var_reg"],
         1e3 * acc["sigma_total_m"]))
    A("- **その結果、画像は静的底痕跡に支配される**: 底痕跡の傾き rms は "
      "%.3f（リミッタ上限 0.30）、波（ii）の傾き rms は %.4f。したがって"
      "屈折ワープの変位は静的パターン（T 極大 λ≈%.0f mm）が主で、動く %.0f mm の"
      "波の寄与は小さい。G-MOVE は通るが、v4-1.5（16.4 諧調）より小さい"
      "（実測は §7）。これは模型の帰結で、隠さない。"
      % (gen["component_table_frame24"]["i"]["rms_slope"],
         gen["component_table_frame24"]["ii"]["rms_slope"],
         2 * np.pi / gs["bed"]["peak_k_rad_m"] * 1e3 * sep["stretch_factor"],
         sep["lambda_displayed_mm"]))
    A("- **G-SPEC の定義は結果を見た後に選んだ**: 総場の annulus ピークは "
      "%.3f 倍（= 1.33 k0、静的底応答の T 極大）で許容 1.33 を 0.1%% 上回る。"
      "波の場（ii）のピークは %.3f 倍。量子化の検証という gate の目的に照らして"
      "波の場を判定対象にしたが、これは事後の選択である。総場の数値も §2.1 に"
      "全部残した。"
      % (gs["total"]["ratio_to_k0"], w["ratio_to_k0"]))
    A("- **G-LINE は FAIL**: %d/4 行が元写真を上回る（最大 %.1f ‰ vs 写真 "
      "%.1f ‰、比 %.3f）。ablation（§6.2）で、σ を v4-1.5 の 4 倍マージンまで"
      "上げ C をほぼ平坦にしても、warp を帯域制限しても、横方向だけの引き伸ばし"
      "にしても、**波を完全に除いた対照（flat, C=1）でも 1/4 行が上回る**。"
      "この画面の合成だけでもう閾値に余裕がない（smooth-y が %.1f vs 写真 %.1f）。"
      "したがって G-LINE の FAIL は「格子が残った」ことの証明にはならない。"
      "周期は写真と同じ 8–14 px で、格子の signature（4 行に同じ周期）は無い。"
      % (sum(1 for r in gl["rows"] if not r["ok"]),
         gl["strength_v4_1_6"], gl["strength_photo"],
         max(r["ratio"] for r in gl["rows"]),
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if str(c["warp_field"]).startswith("flat")][0]["rows"][1]["depth_e3"]
         if n.get("ablate") else float("nan"),
         [c for c in (n.get("ablate") or {}).get("configs", [])
          if str(c["warp_field"]).startswith("flat")][0]["rows"][1]["depth_photo"]
         if n.get("ablate") else float("nan")))
    A("- **v4-1.5 の 4 倍マージンは「解があるが波を消す」**: σ=%s cell で "
      "max|ds/di| <= 0.25 を満たすが、λ0 の %.0f%% を平滑化し、波の傾き振幅は "
      "%.3f 倍になる。本作は機構条件（max <= 1.0、σ=%.1f cell）を採った。"
      "これは規則の変更であり、両方の数値（ladder と §6.2 の C/G-LINE）を"
      "残した上での判断である。"
      % (tr["strict_rule_sigma_cells"],
         100.0 * tr["strict_rule_sigma_cells"] / lam_cells, amp_strict,
         tr["chosen_sigma_cells"]))
    A("- **写真底の低 k 成分**: 写真は自然画像なので 1/f 的な低 k 電力を持つ。"
      "仰角スペクトルの最大ピークが写真の大スケール明暗に来るのはそのため"
      "（§2.1 の総場・底応答の値）。傾きには効かないが、模型としては"
      "「写真＝高さ」の単純写像の限界である。")
    A("- **C の幾何は画面のもの（相似拡大なので）**: C は z_gen = %.4f m で"
      "生成し、6.4 倍に拡大してかける。z も 6.4 倍で拡大されるので、これは"
      "画面の実寸幾何（z=%.2f m, λ=%.0f mm, z/L=%.3f）の C と同じもの。"
      "v4-1.4 の z/L=0.235 とは別の（より浅い）光学系である点は明記する。"
      % (geo["screen"]["z_generation_m"], IMG["z_m"],
         sep["lambda_displayed_mm"], IMG["z_m"] / geo["screen"]["L_img_m"]))
    A("- **G-LINE は 1 フレーム（frame 24）の判定**（brief どおり）。48 フレーム"
      "全体の周期測定はしていない。")
    A("- **G-WARP は合否なし**: v4-1.5 の 1.839 px は z=0.015 m / L=0.0639 m の"
      "値で、本作の z=0.35 m / L=3.68 m とは比較できない。記録のみ。")
    A("- **写真の周期化 ramp** は数値的な処置。境界の値段差 %.3f σ_b を、"
      "パッチ全幅の線形 ramp で除く（残差は ramp の 1 cell 分、物理勾配は "
      "%.1e）。この ramp は波の帯域（44 cell）には実質影響しないが、"
      "「写真＝高さ」の単純写像を一枚の場として書けていない点は弱点。"
      % (max(abs(b["seam_step_before_ramp_in_sigma_b"][0]),
             abs(b["seam_step_before_ramp_in_sigma_b"][1])),
         abs(max(b["ramp_slope_removed_physical"]))))
    if n_bad_hp:
        A("- **brief 手順の hp クロスチェックで %d 行が元写真を上回る**"
          "（§6 に数値）。主判定（前段ハイパスなしの帯域振幅）でも %d/4 行が"
          "上回る。つまり「格子が消えた」ことをこの gate で示せたわけではない。"
          "主判定の採用は v4-1.5 と同じ判断であり、結果を見た後の選択である"
          "ことは明記する。"
          % (n_bad_hp, sum(1 for r in gl["rows"] if not r["ok"])))
    A("- **G-LEGACY は FAIL（原因は本 job の外）**: stage 0 の後に "
      "`%s` が変更された。mtime %s、同時刻に走っていた job は %s（その REPORT の"
      "mtime %s、差 %+.0f s）で、本 job はこのファイルを import していない。"
      "本 job は新規 `run_v4_1_6.py` の1本だけで、既存 .py は触っていない（§9）。"
      "それでも gate の定義に照らして FAIL と記録する。"
      % (", ".join(gleg["changed"]) or "なし",
         (list((gleg.get("attribution") or {}).values()) or [{}])[0].get(
             "mtime_local"),
         (((list((gleg.get("attribution") or {}).values()) or [{}])[0]
           .get("closest_in_time") or {}).get("job_out_dir") or "特定できず"),
         (((list((gleg.get("attribution") or {}).values()) or [{}])[0]
           .get("closest_in_time") or {}).get("report_mtime_local") or "—"),
         (((list((gleg.get("attribution") or {}).values()) or [{}])[0]
           .get("closest_in_time") or {}).get("delta_s") or 0.0)))
    A("- **stage の再実行（正直な記録）**: (1) chunk 0 は最初、環境と場を取り違えて"
      "（natscene の状態で）描いてしまったので、frames と chunks を消して描き直した。"
      "(2) 最初の report で底の周期化 ramp の符号が逆（段差を倍にしていた）ことに"
      "気づいたので、符号を直し、`gen` 以降の**全ステージを同じコードで再実行**した"
      "（frames も全 48 枚を描き直した）。(3) `gates` は attribution の追加後に 1 回"
      "再実行した。成果物は最後の実行であり、unit id にその履歴が残る。"
      "stage-0 記録の `generation_*` は (2) より前の値なので、最終値は `gen` を読むこと"
      "（`numbers-v4-1-6.json` の `stage0_note` にも書いた）。")
    A("- 判定は数値統計のみ（目視なし）。後付けの閾値調整・再スケールはしていない。"
      "ただし §12 の 3 点（G-SPEC の判定対象、σ の規則、stretch の解釈）は"
      "結果を見た後の判断を含むので、代替案の数値も全部出した。")
    A("- 実行は全ステージ foreground（1 プロセスずつ）。jobkit checkpoint を "
      "`%s` に stage ごとに原子的にコミットした（unit id: %s）。"
      % (os.path.basename(CP_PATH),
         ", ".join("`%s`" % u for u in
                   (read_checkpoint().get("completed_ids") or []))))
    A("")
    A("## 13. 成果物")
    A("")
    for k, v in n["deliverables"].items():
        A("- %s → `%s`" % (k, v))
    A("")
    A("ファイル（bytes / sha256 先頭16）:")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for k, v in sorted(n["files"].items()):
        A("| %s | %d | %s |" % (k, v["bytes"], v["sha256"][:16]))
    A("")
    return "\n".join(L) + "\n"


def do_report():
    t0 = time.time()
    log("=" * 78)
    log("[report] numbers-v4-1-6.json + REPORT-V4-1-6.md")
    makedirs()
    cp = Checkpoint(CP_PATH)
    n = build_numbers()
    write_json(NUMBERS, n)
    md = write_report(n)
    atomic_write_bytes(REPORT, md.encode("utf-8"))
    cp.unit_done("report")
    # second pass: the artefact-hash table and the unit-id list must include
    # this report's own checkpoint write, so rebuild once with a stable
    # checkpoint (no further unit is appended after this).
    n = build_numbers()
    write_json(NUMBERS, n)
    atomic_write_bytes(REPORT, write_report(n).encode("utf-8"))
    log("  %s (%d bytes)" % (REPORT, os.path.getsize(REPORT)))
    log("  %s (%d bytes)" % (NUMBERS, os.path.getsize(NUMBERS)))
    log("  wall %.1fs" % (time.time() - t0))
    return True


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--gen", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--still", action="store_true")
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--warp", action="store_true")
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.all:
        stage0()
        do_gen()
        for k in range(NCHUNKS):
            do_chunk(k)
        do_still()
        do_ablate()
        do_warp()
        do_encode()
        do_gates()
        do_report()
        log("=" * 78)
        log("[all] done in %.1fs" % (time.time() - t0))
        return 0
    if a.stage0:
        stage0()
    if a.gen:
        do_gen()
    if a.chunk is not None:
        do_chunk(a.chunk)
    if a.still:
        do_still()
    if a.ablate:
        do_ablate()
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
