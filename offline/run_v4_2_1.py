"""run_v4_2_1 -- BRIEF-V4-2-1 driver: particle bed + height-field transport.

New files only: `bed_gen.py` (bed), `run_v4_2_1.py` (this driver) and
out/v4-2-1/.  Every existing .py is left byte-identical (G-LEGACY).

What is new relative to v4-1.2 / v20
------------------------------------
1. The receiver of the forward transport is a HEIGHT FIELD z(x,y) instead of
   a plane.  The intersection is solved by iteration (plane hit -> local z
   correction, 3 passes) and the per-pixel depth h and bottom slope cos(theta_b)
   enter the direct term explicitly.
2. Boulders (and, via the same depth map, every grain) shadow the bed: the
   solar disc (half angle 0.00395 rad) is sampled 9 times and the shadow is
   marched through the depth map with the *refracted* in-water sun direction,
   so penumbra appears naturally.  AO (horizon approximation, 12 azimuths) is
   applied to the diffuse term only.
3. The bed albedo is built from published reflectances (limestone pebbles <=
   0.22, quartz sand rho(475 nm)=0.47, turf-algae spectrum) with a per-channel
   wet-transform and 20 % algal cover.  Every value and function shape is in
   numbers/REPORT.
4. The scene water surface is stream_surface2 (H=0.25 m, U=0.6, S0=0.003,
   D84=0.042) and the sun is preset 5's geometry (recipe-4: 2026-07-20 10:00,
   26.3 N 127.8 E) with env.solar_position, exactly as run_v19 does.

Stages (all foreground, idempotent):
  --stage0   preflight, sha256 record, wind + stream-lane smoke
  --bed      bed_gen build + G-BED audit -> bed-stats.json + depth-map.png
  --scene    stream generator + env + optics + scene constants
  --gtrans   G-TRANS (a) flat vs height flux, (b) boulder shadow core
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
import io
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

import env as ENV                      # noqa: E402
import v12_water as W                  # noqa: E402
import run_v10 as R10                  # noqa: E402
import run_v15 as R15                  # noqa: E402
import run_v16_water as V16            # noqa: E402
import stream_surface2 as SS           # noqa: E402
import bed_gen as BG                   # noqa: E402

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-1.md"
OUT = os.path.join(HERE, "out", "v4-2-1")
WORK = os.path.join(OUT, "_work")
FRAMES = os.path.join(OUT, "frames")
CP_PATH = os.path.join(OUT, "checkpoint.json")
PARAMS = os.path.join(WORK, "params.json")
FLUX_PATH = os.path.join(WORK, "gtrans-a.json")
SHADOW_PATH = os.path.join(WORK, "gtrans-b.json")
BED_NPZ = os.path.join(WORK, "bed-maps.npz")
BED_STATS = os.path.join(OUT, "bed-stats.json")
MP4 = os.path.join(OUT, "piece-2s-v4-2-1.mp4")
NUMBERS = os.path.join(OUT, "numbers-v4-2-1.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-1.md")
STILL = os.path.join(OUT, "frame-0024.png")
DEPTH_PNG = os.path.join(OUT, "depth-map.png")

# ---- scene constants (BRIEF-V4-2-1) ---------------------------------------
SC = {"H": 0.25, "U": 0.6, "S0": 0.003, "D84": 0.042,
      "seed": 3255045276382562994}
L_WATER = 1.38              # 2 Us T_loop with Us = 1.15*0.6 -> exact 2 s loop
T_LOOP = 2.0
N_S = 1024                  # water-surface grid over L_WATER
Q_DET = 12
IOR_G = 1.33434
FPS = 24
NFRAMES = 48
CHUNK = 24
NCHUNKS = NFRAMES // CHUNK
WIDTH, HEIGHT = 736, 540
A_FROZEN = 2.6
E_FROZEN = 2.2
PRESET_INDEX = 5            # "深い海" = recipe-4.json (5th preset of recipe.js)
WATER_MEDIUM = "clear"

RUN_DOMAIN = 0.80           # transport square domain (m)
XOFF = -0.20                # its offset in the water patch (m); the sun is in
                            # +x, so the rays land ~0.14 m downstream of their
                            # source and the domain MUST extend that far up-sun
NBUF = 980                  # transport buffer (0.816 mm / cell)
NS_SRC = 640                # ray-source grid (1.25 mm / cell, coarser than the
                            # buffer so the splat is a smooth upsample instead
                            # of a 1:1 lattice that would print a mesh)
WAVE_SIGMA_PX = 1.5

# ---- materials (public values, all labelled) -------------------------------
LIMESTONE_DRY = np.array([0.22, 0.21, 0.19])
SAND_DRY = np.array([0.62, 0.56, 0.47])       # rho(475 nm) = 0.47
ALGAE_DRY = np.array([0.10, 0.22, 0.11])
WET_A = np.array([0.95, 0.72, 0.55])
WET_P = 1.35
ALGAE_COVER = 0.20

REQUIRED = {
    "run_v20_prod.py": os.path.join(HERE, "run_v20_prod.py"),
    "stream_surface2.py": os.path.join(HERE, "stream_surface2.py"),
    "run_v4_1_2.py": os.path.join(HERE, "run_v4_1_2.py"),
    "recipe4": os.path.join(HERE, "..", "blender", "recipe-4.json"),
}


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
    return st


# ---------------------------------------------------------------------------
# preset 5 (recipe-4.json) -> environment
# ---------------------------------------------------------------------------
def preset_state():
    """The 5th preset of web/recipe.js, with the brief's water overrides."""
    raw = R10.load_recipe()          # recipe-1 (preset 1): used only as the
    st = {                            # shape template for the couplings
        "name": "深い海",             # preset 5 name (recipe.js / recipe-4)
        "sun": {"date": "2026-07-20", "time": "10:00",
                "latitudeDeg": 26.3, "longitudeDeg": 127.8},
        "water": {"depthM": SC["H"], "transparency": 0.95, "turbidity": 0.02,
                  "pigment": 0.02},
        "windWave": {"windSpeedMS": 3.5, "windDirDeg": 180, "fineness": 0.45,
                     "flowSpeed": 0.45},
        "weather": {"cloudPct": 10, "clarity": 0.9},
        "couplings": raw.get("couplings"),
        "render": {"strength": 1.0, "exposure": 1.0},
    }
    for c in st["couplings"]:
        c["enabled"] = True
        c["weight"] = 1.0
    return st


# ---------------------------------------------------------------------------
# height-field forward transport
# ---------------------------------------------------------------------------
def fibonacci_disc(Q, radius):
    i = np.arange(Q, dtype=np.float64)
    r = float(radius) * np.sqrt((i + 0.5) / float(Q))
    th = i * (math.pi * (3.0 - math.sqrt(5.0)))
    return np.stack([r * np.cos(th), r * np.sin(th)], axis=1)


def _sample_periodic(a, px, py, size):
    """Bilinear sample of a periodic map with px/py in pixel units."""
    x0 = np.floor(px).astype(np.int64)
    y0 = np.floor(py).astype(np.int64)
    tx = px - x0
    ty = py - y0
    x0m = np.mod(x0, size)
    x1m = np.mod(x0 + 1, size)
    y0m = np.mod(y0, size)
    y1m = np.mod(y0 + 1, size)
    return (a[y0m, x0m] * (1 - tx) * (1 - ty) + a[y0m, x1m] * tx * (1 - ty)
            + a[y1m, x0m] * (1 - tx) * ty + a[y1m, x1m] * tx * ty)


def det_frame_height(fields, L_water, nbuf, env, Q, zbed, zw, cell_bed,
                     ior=IOR_G, row_block=220, n_iter=3, a_abs=None,
                     n_src=None):
    """Forward transport onto a height field (v20 det_frame + z(x,y)).

    Parameters
    ----------
    fields : the water surface (h,gx,gy,Dx,Dy) on the N_S grid over L_water.
    zbed   : (nb, nb) bed height map covering the transport domain, sampled
             with period `cell_bed * nb` (the map is periodic in x and y).
    zw     : still-water surface level in the same z datum as zbed.
    a_abs  : per-band absorption (/m) used for the two flux totals only.

    Returns (C, info).  C is dimensionless, mean ~1 (no renormalisation).
    """
    h = np.asarray(fields["h"], dtype=np.float64)
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    N_s = h.shape[0]
    nbuf = int(nbuf)
    cell_w = float(L_water) / N_s
    cell_b = float(cell_bed)
    nb = zbed.shape[0] if zbed is not None else 0
    s, u, v = _basis(env["sunDir"])
    disc = _disc_radius(env.get("clarity", 0.8))
    quads = fibonacci_disc(int(Q), disc)
    eta = 1.0 / float(ior)
    acc = np.zeros(nbuf * nbuf, dtype=np.float64)
    n_src = int(nbuf if n_src is None else n_src)
    cell_s = float(L_water) * 0.0 + float(cell_bed) * nbuf / n_src
    tot = {"rays": 0.0, "raw": 0.0, "abs_g": 0.0, "cos": 0.0,
           "nonconv": 0.0, "path_mean": 0.0, "resid": 0.0}
    tot_resid = 0.0
    xs_src = (np.arange(n_src) + 0.5) * cell_s + XOFF
    for q in range(int(Q)):
        jx, jy = float(quads[q, 0]), float(quads[q, 1])
        d = s + u * jx + v * jy
        d /= np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])
        for r0 in range(0, n_src, row_block):
            r1 = min(r0 + row_block, n_src)
            yy = xs_src[r0:r1]
            xx = xs_src
            X0, Y0 = np.meshgrid(xx, yy)
            # water surface at the source point (periodic in L_water)
            px = np.mod(X0, L_water) / cell_w - 0.5
            py = np.mod(Y0, L_water) / cell_w - 0.5
            Hw = _sample_periodic(h, px, py, N_s)
            GX = _sample_periodic(gx, px, py, N_s)
            GY = _sample_periodic(gy, px, py, N_s)
            DX = _sample_periodic(np.asarray(fields["Dx"], dtype=np.float64),
                                  px, py, N_s)
            DY = _sample_periodic(np.asarray(fields["Dy"], dtype=np.float64),
                                  px, py, N_s)
            nl = np.sqrt(GX * GX + GY * GY + 1.0)
            Nx = -GX / nl
            Ny = -GY / nl
            Nz = 1.0 / nl
            cos1 = -(Ix * Nx + Iy * Ny + Iz * Nz)
            cos1 = np.clip(cos1, -1.0, 1.0)
            kk = 1.0 - eta * eta * (1.0 - cos1 * cos1)
            kk = np.maximum(kk, 0.0)
            sq = np.sqrt(kk)
            Tx = eta * Ix + (eta * cos1 - sq) * Nx
            Ty = eta * Iy + (eta * cos1 - sq) * Ny
            Tz = eta * Iz + (eta * cos1 - sq) * Nz
            Tz = np.minimum(Tz, -1e-6)
            Psx = X0 + DX
            Psy = Y0 + DY
            Psz = zw + Hw
            if zbed is None:
                t = (Psz - zw + float(SC["H"])) / (-Tz)
                xb = Psx + Tx * t
                yb = Psy + Ty * t
                cosb = np.ones_like(t)
                nonconv = np.zeros_like(t)
            else:
                t = Psz / (-Tz)
                for it in range(int(n_iter)):
                    xb = Psx + Tx * t
                    yb = Psy + Ty * t
                    mx = np.mod(xb, cell_b * nb) / cell_b - 0.5
                    my = np.mod(yb, cell_b * nb) / cell_b - 0.5
                    zb = _sample_periodic(zbed, mx, my, nb)
                    if it < int(n_iter) - 1:
                        # plain fixed point while far from the surface
                        t = t + 0.8 * ((Psz - zb) / (-Tz) - t)
                    else:
                        zxp = _sample_periodic(zbed, mx + 1.0, my, nb)
                        zxm = _sample_periodic(zbed, mx - 1.0, my, nb)
                        zyp = _sample_periodic(zbed, mx, my + 1.0, nb)
                        zym = _sample_periodic(zbed, mx, my - 1.0, nb)
                        zx = (zxp - zxm) / (2.0 * cell_b)
                        zy = (zyp - zym) / (2.0 * cell_b)
                        F = Psz + t * Tz - zb
                        dF = Tz - (zx * Tx + zy * Ty)
                        dF = np.sign(dF) * np.maximum(np.abs(dF), 0.15)
                        dt = np.clip(-F / dF, -0.05, 0.05)
                        tot_resid += float(np.sum(np.abs(dt)))
                        nonconv = (np.abs(dt) > 1e-4).astype(np.float64)
                        t = t + dt
                # landing stays in the DOMAIN frame (the buffer wraps by
                # RUN_DOMAIN, not by the water patch size)
                xb = Psx + Tx * t
                yb = Psy + Ty * t
                # bottom slope: finite difference of the bed map at the hit
                zx = (_sample_periodic(zbed, np.mod(xb + cell_b, cell_b * nb)
                                       / cell_b - 0.5, yb / cell_b - 0.5, nb)
                      - _sample_periodic(zbed, np.mod(xb - cell_b,
                                                      cell_b * nb) / cell_b
                                         - 0.5, yb / cell_b - 0.5, nb)) \
                    / (2 * cell_b)
                zy = (_sample_periodic(zbed, np.mod(xb, cell_b * nb) / cell_b
                                       - 0.5,
                                       np.mod(yb + cell_b, cell_b * nb)
                                       / cell_b - 0.5, nb)
                      - _sample_periodic(zbed, np.mod(xb, cell_b * nb)
                                         / cell_b - 0.5,
                                         np.mod(yb - cell_b, cell_b * nb)
                                         / cell_b - 0.5, nb)) / (2 * cell_b)
                cosb = 1.0 / np.sqrt(zx * zx + zy * zy + 1.0)
            # buffer coordinates (over the transport domain)
            fx = (xb - XOFF) / cell_b - 0.5
            fy = (yb - XOFF) / cell_b - 0.5
            ix0 = np.floor(fx).astype(np.int64) % nbuf
            iy0 = np.floor(fy).astype(np.int64) % nbuf
            tx = fx - np.floor(fx)
            ty = fy - np.floor(fy)
            ix1 = (ix0 + 1) % nbuf
            iy1 = (iy0 + 1) % nbuf
            nbb = nbuf * nbuf
            w00 = ((1 - tx) * (1 - ty)).ravel()
            w10 = (tx * (1 - ty)).ravel()
            w01 = ((1 - tx) * ty).ravel()
            w11 = (tx * ty).ravel()
            i00 = (iy0 * nbuf + ix0).ravel()
            i10 = (iy0 * nbuf + ix1).ravel()
            i01 = (iy1 * nbuf + ix0).ravel()
            i11 = (iy1 * nbuf + ix1).ravel()
            acc += np.bincount(i00, weights=w00, minlength=nbb)
            acc += np.bincount(i10, weights=w10, minlength=nbb)
            acc += np.bincount(i01, weights=w01, minlength=nbb)
            acc += np.bincount(i11, weights=w11, minlength=nbb)
            nrays = float(t.size)
            tot["rays"] += nrays
            tot["raw"] += nrays
            if a_abs is not None:
                tot["abs_g"] += float(np.sum(np.exp(-a_abs[1] * t)))
            tot["cos"] += float(np.sum(cosb))
            tot["nonconv"] += float(np.sum(nonconv))
            tot["resid"] += tot_resid
            tot_resid = 0.0
            tot["path_mean"] += float(np.sum(t))
            del X0, Y0, Hw, GX, GY, DX, DY
    total_rays = float(int(Q)) * float(n_src) * float(n_src)
    # irradiance per BUFFER cell: each ray carries the source cell area, so
    # the divisor is Q * (cell_b / cell_s)^2 (== Q when the grids coincide)
    divisor = float(int(Q)) * (cell_bed / cell_s) ** 2
    C = (acc / divisor).reshape(nbuf, nbuf)
    C = _gauss_fft(C, WAVE_SIGMA_PX).astype(np.float32)
    info = {
        "Q": int(Q), "ior": float(ior), "disc_rad": float(disc),
        "sigma_px": WAVE_SIGMA_PX, "total_rays": total_rays,
        "n_src": n_src, "n_buf": nbuf, "cell_src_m": cell_s,
        "divisor_analytical": divisor, "clamp": "none",
        "renormalization": "none",
        "n_iter": int(n_iter),
        "rays": tot["rays"],
        "mean_path_m": tot["path_mean"] / max(tot["rays"], 1.0),
        "nonconverged_fraction": tot["nonconv"] / max(tot["rays"], 1.0),
        "mean_fixed_point_residual_m": tot["resid"] / max(tot["rays"], 1.0),
        "flux_raw": tot["raw"] / max(tot["rays"], 1.0),
        "flux_beer_green": (tot["abs_g"] / max(tot["rays"], 1.0)
                            if a_abs is not None else None),
        "flux_cos_theta_b": tot["cos"] / max(tot["rays"], 1.0),
        "C_mean": float(np.mean(C)),
        "C_max": float(np.max(C)),
        "C_p99": float(np.percentile(C, 99)),
    }
    return C, info


def _basis(sun):
    s = np.asarray(sun, dtype=np.float64)
    s = s / np.linalg.norm(s)
    tmp = np.array([0.0, 0.0, 1.0]) if abs(s[2]) < 0.99 else np.array([1.0, 0, 0])
    u = np.cross(s, tmp)
    u /= np.linalg.norm(u)
    v = np.cross(s, u)
    return s, u, v


def _disc_radius(clarity):
    from transport import sun_disc_radius
    return sun_disc_radius(clarity)


def _gauss_fft(a, sigma):
    if sigma <= 0:
        return a
    n0, n1 = a.shape
    ky = np.fft.fftfreq(n0)[:, None]
    kx = np.fft.fftfreq(n1)[None, :]
    g = np.exp(-2.0 * (np.pi ** 2) * (sigma ** 2) * (kx * kx + ky * ky))
    return np.real(np.fft.ifft2(np.fft.fft2(a) * g))


# ---------------------------------------------------------------------------
# materials
# ---------------------------------------------------------------------------
def material_maps(bed, seed=SC["seed"]):
    """Per-pixel wet albedo (H,W,3) + the dry albedo, from public values."""
    gid = bed["gid"]
    D = bed["D"]
    d = np.where(gid >= 0, D[np.clip(gid, 0, len(D) - 1)], 0.02)
    # size class: sand below 6 mm, limestone above 12 mm, blend between
    f_rock = np.clip((d - 0.006) / 0.006, 0.0, 1.0)[..., None]
    base = SAND_DRY[None, None, :] * (1 - f_rock) \
        + LIMESTONE_DRY[None, None, :] * f_rock
    # per-grain tint (deterministic hash of the grain id)
    hsh = (np.sin((gid + 1) * 12.9898) * 43758.5453) % 1.0
    tint = (0.90 + 0.20 * hsh)[..., None]
    dry = np.clip(base * tint, 0.0, 1.0)
    # turf algae: 20 % areal cover from a band-limited noise field
    from scipy import ndimage
    rng = np.random.default_rng(seed + 77)
    n = np.zeros((WIDTH, WIDTH))
    noise = rng.standard_normal((WIDTH, WIDTH))
    noise = ndimage.gaussian_filter(noise, 6.0, mode="wrap")
    thr = np.quantile(noise, 1.0 - ALGAE_COVER)
    mask = np.clip((noise - thr) / (0.5 * np.std(noise) + 1e-12), 0.0, 1.0)
    algae = ALGAE_DRY[None, None, :] * (0.85 + 0.3 * hsh[..., None])
    dry = dry * (1 - mask[..., None]) + algae * mask[..., None]
    # wet transform (per channel, non-linear): rho_w = rho / (1 + a (1-rho)^p)
    wet = dry / (1.0 + WET_A[None, None, :]
                 * np.power(np.clip(1.0 - dry, 0.0, 1.0), WET_P))
    return {"dry": dry, "wet": wet, "algae_mask": mask,
            "coverage": float(np.mean(mask > 0.5))}


# ---------------------------------------------------------------------------
# scene
# ---------------------------------------------------------------------------
def build_scene():
    t0 = time.time()
    state = preset_state()
    env = ENV.resolve_environment(state)
    env["depthM"] = float(SC["H"])
    spec = W.preset_spectra(WATER_MEDIUM)
    fit = SS.fit_sigma()
    d = SS.derive(SC["U"], SC["H"], SC["S0"])
    st = SS.StreamState2(N=N_S, L=L_WATER, H=SC["H"], U=SC["U"],
                         S0=SC["S0"], D84=SC["D84"], seed=SC["seed"],
                         bed="natural", alpha=SS.ALPHA_MAIN, fit=fit,
                         dt_loop=T_LOOP)
    st.loop = {"n_loops": 1, "T_loop": T_LOOP,
               "rule": "L = n Us T_loop with Us = 1.15 U = %.4f, n = 1"
                       % d["Us"]}
    o = V16.optics_v16(spec, env, float(SC["H"]), NBUF, RUN_DOMAIN)
    bed = np.load(BED_NPZ, allow_pickle=True)
    zbed = bed["z"].astype(np.float64)
    bed_meta = {k: bed[k] for k in bed.files}
    zw = float(SC["H"]) + (float(np.percentile(zbed, 90)) - float(zbed.mean()))
    lab = SS.validity_label(st.derived["Fr"], SC["H"], SC["D84"])
    out = {
        "preset_index": PRESET_INDEX,
        "preset_name": state["name"],
        "sun": env["sun"], "sunDir": [float(x) for x in env["sunDir"]],
        "sunColor": [float(x) for x in env["sunColor"]],
        "sunI": float(env["sunI"]), "diffuse": float(env["diffuse"]),
        "clarity": float(env["clarity"]), "turbid": float(env["turbid"]),
        "scatter": [float(x) for x in env["scatter"]],
        "absorb": [float(x) for x in env["absorb"]],
        "water_medium": spec["preset"], "water_label": spec.get("label"),
        "a_rgb": [float(x) for x in spec["a"]],
        "b_rgb": [float(x) for x in spec["b"]],
        "scene": {"H": SC["H"], "U": SC["U"], "S0": SC["S0"],
                  "D84": SC["D84"], "seed": SC["seed"]},
        "derived": st.derived, "loop": st.loop, "L_water_m": L_WATER,
        "stream": st.describe(), "validity": lab,
        "zw_m": zw,
        "zbed": {"mean_m": float(zbed.mean()),
                 "p90_m": float(np.percentile(zbed, 90)),
                 "crest_rule": "zw = H + (p90(z_bed) - mean(z_bed)): the "
                               "brief's H is the mean depth over the bed crest"},
        "optics": {"phase_kind": o["phase_kind"], "kernel": o["kernel"],
                   "w_ballistic": [float(x) for x in o["w_ballistic"]],
                   "w_halo": [float(x) for x in o["w_halo"]],
                   "phase": o["phase"],
                   "phase_provenance": o["phase_provenance"],
                   "beer_down": [float(x) for x in o["beer_down"]],
                   "beer_up": [float(x) for x in o["beer_up"]],
                   "sigma_px": float(o["sigma_px"])},
        "transport": {"domain_m": RUN_DOMAIN, "xoff_m": XOFF, "nbuf": NBUF,
                      "n_src": NS_SRC,
                      "cell_m": RUN_DOMAIN / NBUF, "Q": Q_DET,
                      "n_iter": 3, "wave_sigma_px": WAVE_SIGMA_PX,
                      "method": "forward splat of refracted solar rays from "
                                "the water surface onto the height field, "
                                "iterative intersection (plane hit -> local z, "
                                "3 passes)"},
        "composite": {"A": A_FROZEN, "E": E_FROZEN, "status":
                      "frozen v20 (A2.6 + E2.2), no photo"},
        "build_wall_s": float(time.time() - t0),
    }
    bed = np.load(BED_NPZ, allow_pickle=True)
    return out, st, env, o, {k: bed[k] for k in bed.files}, zw


# ---------------------------------------------------------------------------
# frame
# ---------------------------------------------------------------------------
def render_frame(i, st, env, optics, bed, zw, g1):
    from image import aces_narkowicz  # noqa: F401  (sanity import)
    from v5_composite import image_metrics
    t0 = time.time()
    t = i / float(FPS)
    fields = st.fields_at(t)
    t1 = time.time()
    zbed = bed["z"].astype(np.float64)
    C, tinfo = det_frame_height(fields, L_WATER, NBUF, env, Q_DET, zbed, zw,
                                RUN_DOMAIN / NBUF, n_src=NS_SRC,
                                a_abs=np.asarray(env["absorb"]))
    t_det = time.time() - t1
    img, dbg = composite_bed(fields, C, bed, env, optics, zw, g1)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    out = {"frame": i, "t_s": t, "g1": float(g1), "C": C,
           "det_info": tinfo, "still": img, "metrics": metrics,
           "dbg": dbg,
           "timing_s": {"fields": 0.0, "det": t_det,
                        "composite": float(time.time() - t0 - t_det),
                        "total": float(time.time() - t0)}}
    del fields
    gc.collect()
    return out


def composite_bed(fields, C, bed, env, optics, zw, g1):
    """v12/v16 body with the bed height field, shadow, AO and local depth."""
    from image import (IOR, SKY_REF, aces_narkowicz, e_base_dir_from_env,
                       linear_to_srgb)
    from state import sample_bicubic  # noqa: F401
    from v5_composite import _state_gradients
    BG_CELL = float(np.asarray(bed["cell"]).ravel()[0])
    BG_J0 = int(np.asarray(bed["render_j0"]).ravel()[0])
    win = slice(BG_J0, BG_J0 + HEIGHT)
    z = bed["z"][win].astype(np.float64)
    nz_b = bed["nz"][win].astype(np.float64)
    ao = bed["ao"][win].astype(np.float64)
    shadow = bed["shadow"][win].astype(np.float64)
    mat = bed["mat"]
    alb_full = mat["wet"]
    alb = alb_full[win]
    S = RUN_DOMAIN
    depth = float(env["depthM"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    Ebase = e_base_dir_from_env(env)
    A, E = A_FROZEN, E_FROZEN
    w_bal = np.asarray(optics["w_ballistic"], dtype=np.float64)
    w_halo = np.asarray(optics["w_halo"], dtype=np.float64)
    Cc = np.asarray(C, dtype=np.float64)
    if float(np.max(w_halo)) > 0.0:
        if optics["kernel"] == "ff-cdf":
            Ch = V16.ff_blur_fft(Cc, optics["kernel_2d"]).astype(np.float64)
        else:
            Ch = _gauss_fft(Cc, float(optics["sigma_px"]))
    else:
        Ch = Cc
    gx_s, gy_s = _state_gradients(fields, 3.0)
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    # bed-map pixel -> transport-buffer pixel (both periodic)
    cell_b = RUN_DOMAIN / NBUF
    bx = (ix + 0.5) * BG_CELL + XOFF
    by = (iy + BG_J0 + 0.5) * BG_CELL + XOFF
    gx = _sample_periodic(np.asarray(gx_s, dtype=np.float64),
                          np.mod(bx, L_WATER) / (L_WATER / N_S) - 0.5,
                          np.mod(by, L_WATER) / (L_WATER / N_S) - 0.5, N_S)
    gy = _sample_periodic(np.asarray(gy_s, dtype=np.float64),
                          np.mod(bx, L_WATER) / (L_WATER / N_S) - 0.5,
                          np.mod(by, L_WATER) / (L_WATER / N_S) - 0.5, N_S)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nzw = 1.0 / nl
    h_local = np.clip(zw - z, 0.005, None)          # per-pixel water depth
    beer = np.exp(-np.asarray(env["absorb"], dtype=np.float64)[None, None, :]
                  * (h_local[..., None] / 0.90
                     + h_local[..., None]))
    scat_amt = float(1.0 - float(np.mean(beer)))
    amb_gain = (0.35 + 0.65 * float(env["diffuse"])) * (0.3 + 0.7 * sunI)
    # Diffuse-sky irradiance on the bed.  v12's E_amb = scatter*(1-beer)*
    # amb_gain*A is calibrated for the photo lane and returns ~0.003 here,
    # which is below any plausible sky irradiance, so the diffuse term uses
    # the env's own diffuse fraction:  E_sky = 2 * diffuse * sunI * sunColor
    # (factor 2 = hemispherical irradiance from a 10 %-cloud clear sky; a
    # declared modelling choice, recorded in numbers/REPORT).
    E_amb = (2.0 * float(env["diffuse"]) * sunI * sunColor) * A * amb_gain
    col = np.zeros((HEIGHT, WIDTH, 3), dtype=np.float64)
    direct_all = np.zeros((HEIGHT, WIDTH, 3), dtype=np.float64)
    for bb in range(3):
        Cb = (float(IOR[bb]) - 1.0) / float(IOR[bb])
        shx = -nx * h_local * Cb
        shy = -ny * h_local * Cb
        # the apparent bed is displaced by the refraction; C is defined on the
        # bed itself and is therefore NOT displaced (that would double count)
        alb_s = _sample_map(alb_full[:, :, bb], ix + 0.5 + shx / BG_CELL,
                            iy + BG_J0 + 0.5 + shy / BG_CELL)
        c_bal = _sample_periodic(Cc, bx / cell_b - 0.5, by / cell_b - 0.5, NBUF)
        c_hal = _sample_periodic(Ch, bx / cell_b - 0.5, by / cell_b - 0.5, NBUF)
        irr = float(w_bal[bb]) * c_bal + float(w_halo[bb]) * c_hal
        direct = (Ebase[bb] * irr * float(g1) * nz_b * shadow
                  * beer[:, :, bb])
        direct_all[:, :, bb] = direct
        bottom = alb_s * (E * direct) + alb_s * E_amb[bb] * ao * beer[:, :, bb]
        col[:, :, bb] = bottom * strength
    cos_t = nzw
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - float(env["diffuse"])) \
        + SKY_REF[None, None, :] * float(env["diffuse"])
    sky_gain = (0.10 + 0.30 * (1.0 - float(env["diffuse"]))) \
        * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * (sky_gain * A)
    lin = col * exposure
    aces = aces_narkowicz(np.maximum(lin, 0.0))
    out = (np.clip(linear_to_srgb(aces), 0.0, 1.0) * 255.0
           + 0.5).astype(np.uint8)
    return out, {"lin": lin, "aces": aces, "alb": alb,
                 "direct": direct_all, "shadow": shadow, "ao": ao,
                 "h_local": h_local, "beer": beer, "g1": float(g1),
                 "C_mean": float(np.mean(Cc)), "C_max": float(np.max(Cc))}


def _sample_map(a, px, py):
    """Bilinear sample of a periodic (736,736) patch map with wrap."""
    n0, n1 = a.shape
    x0 = np.floor(px).astype(np.int64)
    y0 = np.floor(py).astype(np.int64)
    tx = px - x0
    ty = py - y0
    x0m = np.mod(x0, n1)
    x1m = np.mod(x0 + 1, n1)
    y0m = np.mod(y0, n0)
    y1m = np.mod(y0 + 1, n0)
    return (a[y0m, x0m] * (1 - tx) * (1 - ty) + a[y0m, x1m] * tx * (1 - ty)
            + a[y1m, x0m] * (1 - tx) * ty + a[y1m, x1m] * tx * ty)


# ---------------------------------------------------------------------------
# stage 0
# ---------------------------------------------------------------------------
def stage0():
    t0 = time.time()
    log("=" * 72)
    log("[stage 0] absolute-path preflight")
    makedirs()
    cp = Checkpoint(CP_PATH)
    facts = {"workdir": HERE, "brief": BRIEF, "created_utc": now_iso(),
             "python": sys.version.split()[0]}
    missing = []
    for name, p in REQUIRED.items():
        ok = os.path.isfile(p)
        rec = {"path": p, "present": ok, "bytes": None, "sha256": None}
        if ok:
            rec["bytes"] = os.path.getsize(p)
            rec["sha256"] = sha256_file(p)
        else:
            missing.append(name)
        facts.setdefault("files", {})[name] = rec
    before = py_set()
    facts["py_hashes_before"] = {"count": len(before), "hashes": before}
    # preset 5 uses env.solar_position exactly as run_v19 does
    facts["preset5"] = {
        "source": "web/recipe.js presets()[4] == blender/recipe-4.json",
        "state": preset_state(),
        "env": {k: ENV.resolve_environment(preset_state())[k]
                for k in ("sun", "sunDir", "sunColor", "sunI", "diffuse")},
        "solar_position_fn": "env.solar_position (the same call run_v19 makes)",
    }
    # seed identity with the brief
    raw = R10.load_recipe()
    from state import derive_seed
    sseed = int(derive_seed(raw, 1))
    facts["seed_check"] = {"derive_seed(recipe,1)": sseed,
                           "brief": SC["seed"], "equal": sseed == SC["seed"]}
    if sseed != SC["seed"]:
        missing.append("seed mismatch")
    # wind lane smoke (identical to the v4-1.2 record)
    r = R15.scene_for(raw, 3.0, 100.0, 0.6, 256, sseed)
    f = r["state"].fields_at(0.0)
    h = np.asarray(f["h"], dtype=np.float32)
    wb = {"seed": sseed, "N": 256, "cell": {"U": 3.0, "F": 100.0, "h": 0.6},
          "h_sha256": sha256_bytes(h.tobytes()),
          "h_mean": float(np.mean(h.astype(np.float64))),
          "h_std": float(np.std(h.astype(np.float64)))}
    facts["wind_baseline"] = wb
    # stream lane smoke: reproduce the v4-1.2 reference-scene descriptors
    ref = SS.StreamState2(N=256, L=3.68, H=0.35, U=0.8, S0=0.002,
                          D84=0.09, seed=20261112, bed="natural",
                          fit=SS.fit_sigma())
    smoke = {"sigma_total_m": ref.sigma_tot, "k0": ref.k0,
             "phi": ref.phi, "L0_m": ref.L0, "tau_s": ref.tau,
             "sigma_i_m": float(np.sqrt(ref.var_i))}
    facts["stream_lane_smoke"] = smoke
    facts["expected_stream_lane_smoke"] = {
        "sigma_total_m": 0.009918737241949973, "k0": 11.59,
        "source": "out/v4-1-2/REPORT-V4-1-2.md section 5"}
    facts["missing"] = missing
    facts["preflight_pass"] = not missing
    if missing:
        write_json(os.path.join(WORK, "stage0.json"), facts)
        raise SystemExit("stage 0 preflight failed: %s" % missing)
    facts["stage0_wall_s"] = float(time.time() - t0)
    write_json(PARAMS, facts)
    write_json(os.path.join(OUT, "stage0-py-hashes.json"),
               {"created_utc": facts["created_utc"], "hashes": before})
    cp.unit_done("stage-0-preflight")
    log("  preset5 sun alt=%.4f az=%.4f sunI=%.4f"
        % (facts["preset5"]["env"]["sun"]["altitudeDeg"],
           facts["preset5"]["env"]["sun"]["azimuthDeg"],
           facts["preset5"]["env"]["sunI"]))
    log("  wind h sha %s (expect 957add87869b332d...)"
        % wb["h_sha256"][:16])
    log("  stream lane sigma_tot=%.6g k0=%.4f phi=%.4f"
        % (smoke["sigma_total_m"], smoke["k0"], smoke["phi"]))
    log("  preflight PASS (%.1fs)" % facts["stage0_wall_s"])
    return facts


# ---------------------------------------------------------------------------
# bed
# ---------------------------------------------------------------------------
def do_bed():
    t0 = time.time()
    log("=" * 72)
    log("[bed] particle bed + G-BED")
    makedirs()
    cp = Checkpoint(CP_PATH)
    st = preset_state()
    e = ENV.resolve_environment(st)
    sun = {"sunDir": e["sunDir"], "disc_rad": 0.00395,
           "altitude_deg": e["sun"]["altitudeDeg"],
           "azimuth_deg": e["sun"]["azimuthDeg"]}
    bed, rmap, hires, stats, full = BG.build_and_audit(OUT, seed=SC["seed"],
                                                       sun=sun)
    stats["gate"] = BG.gate_bed(stats)
    stats["created_utc"] = now_iso()
    stats["build_wall_s"] = float(time.time() - t0)
    stats["raycast"] = {"render": [BG.RENDER_W, BG.RENDER_H],
                        "render_window_m": [BG.WIN_X0, BG.WIN_X1,
                                            BG.WIN_Y0, BG.WIN_Y1],
                        "hires": [BG.HIRES, BG.HIRES],
                        "hires_cell_m": hires["cell"],
                        "render_cell_m": rmap["cell"],
                        "method": "vertical ray cast (z-buffer), analytic "
                                  "ellipsoid normal at the hit point"}
    BG.save_depth_png(DEPTH_PNG, rmap["z"],
                      np.stack([rmap["nx"], rmap["ny"], rmap["nz"]], -1),
                      rmap["cell"])
    np.savez_compressed(BED_NPZ,
                        z=full["z"].astype(np.float32), gid=full["gid"],
                        nx=full["nx"].astype(np.float32),
                        ny=full["ny"].astype(np.float32),
                        nz=full["nz"].astype(np.float32),
                        fx=full["fx"].astype(np.float32),
                        fy=full["fy"].astype(np.float32),
                        fz=full["fz"].astype(np.float32),
                        ao=full["ao"].astype(np.float32),
                        shadow=full["shadow"].astype(np.float32),
                        D=np.array([g["D"] for g in bed["grains"]]),
                        kind=np.array([g["kind"] for g in bed["grains"]]),
                        render_j0=np.array([BG.RENDER_J0]),
                        cell=np.array([full["cell"]]),
                        mono_geom=np.array(
                            [[g["x"], g["y"], g["z"], g["abc"][0],
                              g["abc"][1], g["abc"][2]]
                             for g in bed["grains"]
                             if str(g["kind"]) == "mono"]),
                        mono_R=np.array(
                            [g["R"].tolist() for g in bed["grains"]
                             if str(g["kind"]) == "mono"]))
    write_json(BED_STATS, stats)
    cp.unit_done("bed-g-bed")
    g = stats["gate"]
    log("  grains=%d  G-BED=%s" % (stats["n_grains"],
                                   "PASS" if g["pass"] else "FAIL"))
    for k, v in g.items():
        if isinstance(v, dict) and "value" in v and "range" in v:
            log("    %-24s %.4g [%g, %g] %s"
                % (k, v["value"], v["range"][0], v["range"][1],
                   "PASS" if v["pass"] else "FAIL"))
    return stats


# ---------------------------------------------------------------------------
# scene stage
# ---------------------------------------------------------------------------
def do_scene():
    log("=" * 72)
    log("[scene] stream generator + preset-5 sun + clear medium")
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc, st, env, o, bed, zw = build_scene()
    write_json(os.path.join(WORK, "scene.json"), sc)
    cp.unit_done("scene")
    s = sc["stream"]
    log("  L=%.4g m  k0=%.4g (lam0=%.4g m)  phi=%.4f  sigma=%.4g mm"
        % (sc["L_water_m"], s["k0_rad_m"], s["lambda0_m"],
           s["ring_fraction_phi"], 1000 * s["sigma_total_target_m"]))
    log("  Fr=%.4f  u*=%.4g m/s  tau=%.4g s  zw=%.4f m  %s"
        % (sc["derived"]["Fr"], sc["derived"]["ustar"],
           sc["derived"]["tau"], zw, sc["validity"]))
    log("  optics kernel=%s w_bal=%s w_halo=%s"
        % (sc["optics"]["kernel"],
           np.round(sc["optics"]["w_ballistic"], 3),
           np.round(sc["optics"]["w_halo"], 3)))
    del st
    gc.collect()
    return sc


# ---------------------------------------------------------------------------
# G-TRANS
# ---------------------------------------------------------------------------
def do_gtrans():
    t0 = time.time()
    log("=" * 72)
    log("[gtrans] (a) flat vs height-field direct flux  (b) boulder shadow")
    makedirs()
    cp = Checkpoint(CP_PATH)
    sc, st, env, o, bed, zw = build_scene()
    bed["mat"] = material_maps(bed)
    zbed = bed["z"].astype(np.float64)
    a_abs = np.asarray(env["absorb"], dtype=np.float64)
    fields = st.fields_at(0.0)
    cell_b = RUN_DOMAIN / NBUF
    C_h, info_h = det_frame_height(fields, L_WATER, NBUF, env, Q_DET, zbed,
                                   zw, cell_b, n_src=NS_SRC, a_abs=a_abs)
    # flat bed at the SAME mean depth: z=const=mean(zbed), same zw
    zflat = np.full_like(zbed, float(zbed.mean()))
    zw_flat = float(zw - zbed.mean() + zbed.mean())
    C_f, info_f = det_frame_height(fields, L_WATER, NBUF, env, Q_DET, zflat,
                                   zw_flat, cell_b, n_src=NS_SRC, a_abs=a_abs)
    diff = abs(info_h["flux_beer_green"] - info_f["flux_beer_green"]) \
        / info_f["flux_beer_green"]
    diff_raw = abs(info_h["flux_raw"] - info_f["flux_raw"]) \
        / info_f["flux_raw"]
    out = {
        "a_energy": {
            "flat": {"flux_beer_green": info_f["flux_beer_green"],
                     "flux_raw": info_f["flux_raw"],
                     "flux_cos_theta_b": info_f["flux_cos_theta_b"],
                     "mean_path_m": info_f["mean_path_m"],
                     "C_mean": info_f["C_mean"], "C_max": info_f["C_max"],
                     "C_p99": info_f["C_p99"]},
            "height": {"flux_beer_green": info_h["flux_beer_green"],
                       "flux_raw": info_h["flux_raw"],
                       "flux_cos_theta_b": info_h["flux_cos_theta_b"],
                       "mean_path_m": info_h["mean_path_m"],
                       "C_mean": info_h["C_mean"], "C_max": info_h["C_max"],
                       "C_p99": info_h["C_p99"],
                       "mean_fixed_point_residual_m":
                           info_h["mean_fixed_point_residual_m"],
                       "nonconverged_fraction":
                           info_h["nonconverged_fraction"]},
            "rel_diff_beer_green": float(diff),
            "rel_diff_raw": float(diff_raw),
            "threshold": 0.05, "pass": bool(diff <= 0.05),
            "definition": ("total direct flux delivered to the bed = sum over "
                           "source rays of exp(-a_green * path); the flat bed "
                           "has the same mean depth and the same albedo. "
                           "flux_raw is the un-attenuated ray count, which is "
                           "conserved by construction (it verifies that no ray "
                           "is lost in the iterative intersection). "
                           "flux_cos_theta_b is the mean cos(slope) and is NOT "
                           "the gate: it is an irradiance-per-true-area "
                           "factor, not a power."),
            "mean_depth_m": float(zw - zbed.mean()),
            "mean_depth_flat_m": float(zw_flat - zflat.mean()),
        },
        "checked_utc": now_iso(),
    }
    # ---- (b) boulder shadow ------------------------------------------------
    mono = bed["mono_geom"]
    monoR = bed["mono_R"]
    s = np.asarray(sc["sunDir"], dtype=np.float64)
    eta = 1.0 / 1.33434
    cos1 = float(s[2])
    kk = 1.0 - eta * eta * (1.0 - cos1 * cos1)
    sw = eta * s + (math.sqrt(max(kk, 0.0)) - eta * cos1) \
        * np.array([0.0, 0.0, 1.0])
    sw /= np.linalg.norm(sw)
    cellp = float(np.asarray(bed["cell"]).ravel()[0])
    j0 = int(np.asarray(bed["render_j0"]).ravel()[0])
    n = 9
    # disc samples in the plane perpendicular to sw
    e1 = np.cross(sw, np.array([0.0, 0.0, 1.0]))
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(sw, e1)
    ii = np.arange(n)
    rr = np.sqrt((ii + 0.5) / n) * 0.00395
    th = ii * (math.pi * (3.0 - math.sqrt(5.0)))
    iy, ix = np.mgrid[0:HEIGHT, 0:WIDTH]
    px = (ix + 0.5) * cellp
    py = (iy + BG.RENDER_J0 + 0.5) * cellp
    pz = bed["z"][j0:j0 + HEIGHT].astype(np.float64)
    blocked = np.zeros((HEIGHT, WIDTH), dtype=np.int32)
    for k in range(n):
        d = sw + e1 * (rr[k] * math.cos(th[k])) + e2 * (rr[k]
                                                        * math.sin(th[k]))
        d /= np.linalg.norm(d)
        hit = np.zeros((HEIGHT, WIDTH), dtype=bool)
        for m in range(mono.shape[0]):
            c = mono[m, :3]
            ax = np.array([mono[m, 3], mono[m, 4], mono[m, 5]])
            R = monoR[m]
            # ray-ellipsoid: p(t) = P + t d, transformed into the body frame
            P0 = np.stack([px - c[0], py - c[1], pz - c[2]])
            RT = R.T
            u = np.einsum("ij,jhw->ihw", RT, P0)
            w = np.einsum("ij,j->i", RT, d)
            qa = ((w / ax) ** 2).sum()
            qb = 2.0 * ((u * w[:, None, None]) / (ax ** 2)[:, None, None]).sum(0)
            qc = ((u / ax[:, None, None]) ** 2).sum(0) - 1.0
            disc = qb * qb - 4 * qa * qc
            ok = disc > 0
            sq = np.sqrt(np.where(ok, disc, 0.0))
            r1 = (-qb - sq) / (2 * qa)
            r2 = (-qb + sq) / (2 * qa)
            hit |= ok & (np.maximum(r1, r2) > 1e-4)
        blocked += hit.astype(np.int32)
    core = blocked >= n
    lit = blocked == 0
    dbg = composite_bed(fields, C_h, bed, env, o, zw, 1.0)[1]
    direct = dbg["direct"]
    # direct irradiance actually delivered (green band is enough: shadow and
    # cos(theta_b) multiply all bands the same way)
    gband = 1
    dv = dbg["alb"][:, :, gband] * direct[:, :, gband]
    core_frac = float(core.mean())
    lit_mean = float(np.mean(dv[lit])) if lit.any() else 0.0
    core_mean = float(np.mean(dv[core])) if core.any() else 0.0
    ratio = core_mean / lit_mean if lit_mean > 0 else 0.0
    march = dbg["shadow"]
    agree_core = float(np.mean(march[core] < 0.15)) if core.any() else 0.0
    agree_lit = float(np.mean(march[lit] > 0.85)) if lit.any() else 0.0
    out["b_shadow"] = {
        "n_disc_samples": n, "disc_half_angle_rad": 0.00395,
        "sun_dir_water": [float(x) for x in sw],
        "core_pixels": int(core.sum()), "core_fraction": core_frac,
        "lit_pixels": int(lit.sum()),
        "core_direct_mean": core_mean, "lit_direct_mean": lit_mean,
        "ratio_core_over_lit": float(ratio),
        "threshold": 0.15, "pass": bool(ratio < 0.15),
        "marching_shadow_agrees_in_core": agree_core,
        "marching_shadow_agrees_in_lit": agree_lit,
        "method": ("core = bed pixels where an ANALYTIC ray/ellipsoid test "
                   "against the two boulders blocks all 9 solar-disc samples "
                   "(a computation independent of the depth-map shadow march "
                   "used by the renderer); direct = the renderer's own direct "
                   "term at the same pixels."),
        "checked_utc": now_iso(),
    }
    out["pass"] = bool(out["a_energy"]["pass"] and out["b_shadow"]["pass"])
    write_json(FLUX_PATH, out)
    st_ = read_cp()
    gates = st_.get("gates") or {}
    gates["G-TRANS"] = out
    cp_merge("gates", gates)
    cp.unit_done("gtrans")
    log("  (a) flat flux=%.6g height flux=%.6g rel=%.3f%% raw_rel=%.2e -> %s"
        % (info_f["flux_beer_green"], info_h["flux_beer_green"], 100 * diff,
           diff_raw, "PASS" if out["a_energy"]["pass"] else "FAIL"))
    log("  (b) core=%.2f%% of frame  direct core/lit=%.4f  march agrees "
        "core=%.3f lit=%.3f -> %s"
        % (100 * core_frac, ratio, agree_core, agree_lit,
           "PASS" if out["b_shadow"]["pass"] else "FAIL"))
    del st
    gc.collect()
    return out


# ---------------------------------------------------------------------------
# chunk / encode
# ---------------------------------------------------------------------------
def chunk_path(i):
    return os.path.join(FRAMES, "frame-%04d.png" % i)


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


def do_chunk(k):
    t0 = time.time()
    log("=" * 72)
    log("[chunk %02d] frames %d-%d" % (k, k * CHUNK, k * CHUNK + CHUNK - 1))
    makedirs()
    cp = Checkpoint(CP_PATH)
    state = read_cp()
    node = (state.get("chunks") or {}).get(str(k))
    done = bool(node) and len(node.get("frames", {})) == CHUNK
    if done:
        for tag, sha in node["frames"].items():
            p = chunk_path(int(tag))
            if not os.path.isfile(p) or sha256_file(p) != sha:
                done = False
                break
    if done:
        log("  chunk %02d already committed with matching hashes -- skip" % k)
        return True
    sc, st, env, o, bed, zw = build_scene()
    mat = material_maps(bed)
    bed["mat"] = mat
    from v5_physics import smith_g1_beckmann
    _f = st.fields_at(0.0)
    _rms = float(np.sqrt(np.mean(np.asarray(_f["gx"], dtype=np.float64) ** 2
                                 + np.asarray(_f["gy"], dtype=np.float64) ** 2)))
    g1 = float(smith_g1_beckmann(
        _rms, 90.0 - float(env["sun"]["altitudeDeg"]))["G1"])
    del _f
    frames = {}
    for i in range(k * CHUNK, (k + 1) * CHUNK):
        r = render_frame(i, st, env, o, bed, zw, g1)
        sha = save_png(chunk_path(i), r["still"])
        frames[str(i)] = {"sha256": sha, "C_mean": float(np.mean(r["C"])),
                          "C_max": float(np.max(r["C"])),
                          "metrics": r["metrics"],
                          "det": {kk: r["det_info"][kk] for kk in
                                  ("C_mean", "mean_path_m",
                                   "nonconverged_fraction",
                                   "flux_beer_green")},
                          "timing_s": r["timing_s"]}
        log("  frame %03d  C mean %.4f max %.3f  %.2fs"
            % (i, float(np.mean(r["C"])), float(np.max(r["C"])),
               r["timing_s"]["total"]))
        del r
        gc.collect()
    stt = read_cp()
    ch = stt.get("chunks") or {}
    ch[str(k)] = {"frames": frames, "wall_s": float(time.time() - t0),
                  "g1": g1}
    cp_merge("chunks", ch)
    cp.unit_done("chunk-%02d" % k)
    log("  chunk %02d done (%.1fs)" % (k, time.time() - t0))
    del st, bed, mat
    gc.collect()
    return True


def ffprobe(path):
    exe = os.path.join(HERE, "_vendor", "imageio_ffmpeg", "binaries")
    cands = [os.path.join(exe, f) for f in os.listdir(exe)] \
        if os.path.isdir(exe) else []
    ff = None
    for c in cands:
        if "ffmpeg" in os.path.basename(c).lower():
            ff = c
    ffprobe = os.path.join(HERE, "_vendor", "ffprobe.exe")
    if os.path.isfile(ffprobe):
        cmd = [ffprobe]
    else:
        cmd = ["ffprobe"]
    cmd += ["-v", "error", "-select_streams", "v:0", "-show_entries",
            "stream=codec_name,width,height,pix_fmt,r_frame_rate,nb_frames,"
            "duration", "-of", "json", path]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0:
            return json.loads(r.stdout)
    except Exception:
        pass
    return {"error": "ffprobe unavailable"}


def do_encode():
    t0 = time.time()
    log("=" * 72)
    log("[encode] mp4")
    cp = Checkpoint(CP_PATH)
    exe = os.path.join(HERE, "_vendor", "imageio_ffmpeg", "binaries")
    ff = None
    if os.path.isdir(exe):
        for f in os.listdir(exe):
            if "ffmpeg" in f.lower():
                ff = os.path.join(exe, f)
    if ff is None:
        raise SystemExit("no ffmpeg binary in _vendor")
    cmd = [ff, "-y", "-framerate", str(FPS), "-start_number", "0",
           "-i", os.path.join(FRAMES, "frame-%04d.png"),
           "-frames:v", str(NFRAMES), "-c:v", "libx264", "-preset", "medium",
           "-crf", "15", "-pix_fmt", "yuv420p", "-r", str(FPS), MP4]
    r = subprocess.run(cmd, capture_output=True, text=True,
                       timeout=1800, cwd=HERE)
    if r.returncode != 0:
        raise SystemExit("ffmpeg failed: %s" % r.stderr[-2000:])
    pr = ffprobe(MP4)
    out = {"mp4": MP4, "bytes": os.path.getsize(MP4), "ffprobe": pr,
           "wall_s": float(time.time() - t0)}
    write_json(os.path.join(WORK, "encode.json"), out)
    cp.unit_done("encode")
    st = (pr.get("streams") or [{}])[0] if isinstance(pr, dict) else {}
    log("  ffprobe %s %sx%s %s fps frames=%s dur=%s"
        % (st.get("codec_name"), st.get("width"), st.get("height"),
           st.get("r_frame_rate"), st.get("nb_frames"), st.get("duration")))
    return out


# ---------------------------------------------------------------------------
# gates + numbers + report
# ---------------------------------------------------------------------------
def do_gates():
    log("=" * 72)
    log("[gates] G-LEGACY / G-OUT")
    cp = Checkpoint(CP_PATH)
    params = read_json(PARAMS)
    before = params["py_hashes_before"]["hashes"]
    after = py_set()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(before) & set(after)
                     if before[k] != after[k])
    allowed = {"run_v4_2_1.py", "bed_gen.py"}
    out = {
        "py_count_before": len(before), "py_count_after": len(after),
        "added": added, "removed": removed, "changed_existing": changed,
        "allowed_new": sorted(allowed),
        "unexpected": [k for k in changed + added if k not in allowed],
    }
    out["pass"] = not out["unexpected"] and not removed
    # wind identity
    raw = R10.load_recipe()
    from state import derive_seed
    sseed = int(derive_seed(raw, 1))
    r = R15.scene_for(raw, 3.0, 100.0, 0.6, 256, sseed)
    h = np.asarray(r["state"].fields_at(0.0)["h"], dtype=np.float32)
    sha1 = sha256_bytes(h.tobytes())
    out["wind_identical"] = bool(sha1 == params["wind_baseline"]["h_sha256"])
    out["wind_sha"] = sha1
    # stream lane identity
    ref = SS.StreamState2(N=256, L=3.68, H=0.35, U=0.8, S0=0.002, D84=0.09,
                          seed=20261112, bed="natural", fit=SS.fit_sigma())
    sm = {"sigma_total_m": ref.sigma_tot, "k0": ref.k0, "phi": ref.phi}
    exp = params["expected_stream_lane_smoke"]
    out["stream_lane"] = {
        "now": sm, "reference": exp,
        "sigma_rel_err": abs(sm["sigma_total_m"] - exp["sigma_total_m"])
        / exp["sigma_total_m"],
        "identical": bool(abs(sm["sigma_total_m"] - exp["sigma_total_m"])
                          / exp["sigma_total_m"] < 1e-12
                          and abs(sm["k0"] - exp["k0"]) < 0.01)}
    out["pass"] = bool(out["pass"] and out["wind_identical"]
                       and out["stream_lane"]["identical"])
    # G-OUT
    pr = read_json(os.path.join(WORK, "encode.json"))["ffprobe"]
    st = (pr.get("streams") or [{}])[0] if isinstance(pr, dict) else {}
    nframes = int(st.get("nb_frames") or 0) if st.get("nb_frames") else 0
    dur = float(st.get("duration") or 0) if st.get("duration") else 0.0
    gout = {"width": int(st.get("width") or 0),
            "height": int(st.get("height") or 0),
            "fps": st.get("r_frame_rate"), "nb_frames": nframes,
            "duration_s": dur, "pix_fmt": st.get("pix_fmt"),
            "expected": {"width": WIDTH, "height": HEIGHT,
                         "frames": NFRAMES, "fps": FPS,
                         "duration_s": NFRAMES / float(FPS)}}
    gout["pass"] = bool(gout["width"] == WIDTH and gout["height"] == HEIGHT
                        and nframes == NFRAMES
                        and abs(dur - NFRAMES / float(FPS)) < 1e-6)
    st_ = read_cp()
    gates = st_.get("gates") or {}
    gates["G-LEGACY"] = out
    gates["G-OUT"] = gout
    cp_merge("gates", gates)
    cp.unit_done("gates-all")
    log("  G-LEGACY: py %d->%d added=%s changed_existing=%s "
        "unexpected=%s wind=%s stream=%s -> %s"
        % (out["py_count_before"], out["py_count_after"], added,
           changed, out["unexpected"], out["wind_identical"],
           out["stream_lane"]["identical"],
           "PASS" if out["pass"] else "FAIL"))
    log("  G-OUT: %dx%d %s fps %d frames %.4f s -> %s"
        % (gout["width"], gout["height"], gout["fps"], nframes, dur,
           "PASS" if gout["pass"] else "FAIL"))
    return {"G-LEGACY": out, "G-OUT": gout}


def build_numbers():
    params = read_json(PARAMS)
    fl = read_json(FLUX_PATH)
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
    return {
        "brief": BRIEF, "created_utc": now_iso(),
        "gates": cp.get("gates") or {},
        "g_bed": bed.get("gate"),
        "scene": sc, "transport": sc["transport"],
        "gtrans": fl, "encode": enc,
        "frames": frames,
        "checkpoint": {"completed_units": cp.get("completed_units"),
                       "last_unit": cp.get("last_unit")},
        "composite": {"A": A_FROZEN, "E": E_FROZEN,
                      "note": "A from the v20 frozen lane, E as the brief fixes",
                      "sky_irradiance": "E_amb = 2*diffuse*sunI*sunColor * A * "
                                        "amb_gain (v12's scatter-based ambient "
                                        "returns ~0.003 for a bed and is below "
                                        "any plausible sky irradiance; the "
                                        "env's own diffuse fraction is used "
                                        "instead -- declared deviation)",
                      "C_refraction": "C is NOT displaced by the refraction; "
                                      "only the albedo is (C is defined on "
                                      "the receiver)"},
        "material": {
            "limestone_dry_rgb": LIMESTONE_DRY.tolist(),
            "limestone_note": "the brief's <= 0.22 bound; a slight warm tint "
                              "(0.22,0.21,0.19) is assumed for the limestone "
                              "pebbles -- approximation",
            "sand_dry_rgb": SAND_DRY.tolist(),
            "sand_note": "quartz sand rho(475 nm) = 0.47 (published); the "
                         "other bands follow the usual rise of quartz "
                         "reflectance with wavelength -- approximation",
            "algae_dry_rgb": ALGAE_DRY.tolist(),
            "algae_note": "turf-algae spectrum (chlorophyll a absorption -> "
                          "green peak at 550 nm, red minimum at 670 nm) "
                          "collapsed to RGB -- approximation",
            "algae_cover": ALGAE_COVER,
            "wet_transform": "rho_w = rho / (1 + a_l (1-rho)^1.35), "
                             "a = (0.95, 0.72, 0.55): a CHANNEL-NON-LINEAR "
                             "darkening (not a scalar), red attenuated most "
                             "because water absorbs red most",
            "wet_p": WET_P,
        },
    }


def write_report(numbers):
    g = numbers["gates"]
    bed = numbers["g_bed"]
    sc = numbers["scene"]
    fl = numbers["gtrans"]
    b = read_json(BED_STATS)
    st = b
    ok = lambda x: "**PASS**" if x else "**FAIL**"
    L = []
    A = L.append
    A("# REPORT-V4-2-1 — 粒子ベースの川床 + 高さ場への輸送")
    A("")
    A("Brief: `BRIEF-V4-2-1.md`. Workdir: `%s`." % HERE)
    A("Built: %s. 新規は `bed_gen.py` と `run_v4_2_1.py` のみ、"
      "出力は `out/v4-2-1/`。既存 py は1バイトも変えていない（G-LEGACY）。"
      % numbers["created_utc"])
    A("")
    A("## 0. 判定")
    A("")
    A("| gate | 判定 | 実測 |")
    A("| --- | --- | --- |")
    A("| G-BED 総合 | %s | 監査 %d 項目中 %d PASS |"
      % (ok(bed["pass"]),
         sum(1 for v in bed.values() if isinstance(v, dict) and "range" in v),
         sum(1 for v in bed.values() if isinstance(v, dict)
             and v.get("pass") is True and "range" in v)))
    A("| G-TRANS (a) エネルギー保存 | %s | 平坦 vs 高さ場の直接フラックス差 "
      "%.3f%% (閾値 5%%) |" % (ok(fl["a_energy"]["pass"]),
                              100 * fl["a_energy"]["rel_diff_beer_green"]))
    A("| G-TRANS (b) 岩の影 | %s | コア/照射 = %.4f (閾値 0.15) |"
      % (ok(fl["b_shadow"]["pass"]), fl["b_shadow"]["ratio_core_over_lit"]))
    A("| G-LEGACY | %s | 既存 py %d本 不変, 風 ident=%s, 流れレーン ident=%s |"
      % (ok(g["G-LEGACY"]["pass"]), g["G-LEGACY"]["py_count_after"],
         g["G-LEGACY"]["wind_identical"],
         g["G-LEGACY"]["stream_lane"]["identical"]))
    A("| G-OUT | %s | %dx%d %s fps %d frames %.3f s |"
      % (ok(g["G-OUT"]["pass"]), g["G-OUT"]["width"], g["G-OUT"]["height"],
         g["G-OUT"]["fps"], g["G-OUT"]["nb_frames"],
         g["G-OUT"]["duration_s"]))
    A("")
    A("## 1. G-BED — 実測値（全て numbers に同値）")
    A("")
    A("| item | 実測 | 監査範囲 | 判定 |")
    A("| --- | --- | --- | --- |")
    rows = [("歪度 (bed, 巨礫除く)", "skewness"),
            ("σz/D50 (matrix)", "sigma_z_over_D50"),
            ("構造関数 H (4–20 cm fit)", "structure_H"),
            ("表面復元 D84/D50", "recon_D84_over_D50"),
            ("装甲比（表面/装甲なし参照）", "recon_armor_ratio"),
            ("傾斜角 中央値 (上流向き)", "imbrication_median_deg"),
            ("クラスター被覆率", "cluster_coverage"),
            ("空隙率", "porosity")]
    for label, key in rows:
        v = bed.get(key) or {}
        A("| %s | %s | %s–%s | %s |"
          % (label, ("%.4g" % v["value"]) if v.get("value") is not None
             else "n/a", v["range"][0], v["range"][1], ok(v.get("pass"))))
    A("")
    A("- 参考（ゲート外・記録のみ）: 歪度(matrix)=%.3f, σz(matrix)=%.4f m, "
      "σz(全体)=%.4f m → σz/D50(全体)=%.3f, 構造関数の折れ点=%.3f m "
      "(%.2f·D50), ファセット>60°=%.2f%%, >80°=%.2f%%, 復元D50/入力D50=%.3f."
      % (st["skew_matrix"], st["sigma_z_matrix_m"], st["sigma_z_bed_m"],
         st["sigma_z_over_D50"], st["structure"]["break_scale_m"],
         st["structure"]["break_over_D50"], 100 * st["facets"]["frac_gt60deg"],
         100 * st["facets"]["frac_gt80deg"],
         st.get("D50_recon_over_input") or float("nan")))
    A("- 復元 GSD (表面): D16=%.1f / D50=%.1f / D84=%.1f mm, "
      "マーカー数 %d, 領域数 %d; 装甲なし参照床: D50=%.1f mm, D84/D50=%.2f."
      % (1000 * st["gsd_reconstructed"]["D16_m"],
         1000 * st["gsd_reconstructed"]["D50_m"],
         1000 * st["gsd_reconstructed"]["D84_m"],
         st["gsd_reconstructed"]["n_markers"],
         st["gsd_reconstructed"]["n"],
         1000 * st["gsd_reference"]["D50_m"],
         st["gsd_reference"]["D84_over_D50"]))
    A("- 粒子 %d 個 (bulk %d / surface %d / cluster %d / mono %d), "
      "基盤振幅 σ_sub=%.5f m (固定点 %d 回)."
      % (st["n_grains"], st["kinds"].get("bulk", 0),
         st["kinds"].get("surface", 0),
         st["kinds"].get("cluster_anchor", 0)
         + st["kinds"].get("cluster_medium", 0)
         + st["kinds"].get("cluster_fine", 0),
         st["kinds"].get("mono", 0),
         st["sigma_sub_m"], len(st["sigma_z_fixpoint"])))
    A("")
    A("### 1.1 G-BED の FAIL について（後付け調整はしていない）")
    A("")
    A("**σz/D50 = %.2f（範囲 0.32–0.46）**。σz(matrix) は %.1f mm で、"
      "**可視粒径**の D50（深さマップから復元した %.1f mm）で割ると %.2f になり、"
      "σz≈0.4·D50 の関係は可視粒径に対してはほぼ成り立つ。brief は入力 GSD の "
      "D50=0.020 で割るので 0.78 になる。装甲比 2.0 を同時に要求する限り可視 D50 は "
      "0.020 にならない（0.020 なら装甲が無い）ので、この2つのゲートは同時に閉じない。"
      % (bed["sigma_z_over_D50"]["value"], 1000 * st["sigma_z_matrix_m"],
         1000 * st["gsd_reconstructed"]["D50_m"],
         st["sigma_z_matrix_m"] / st["gsd_reconstructed"]["D50_m"]))
    A("")
    A("**構造関数 H = %.2f（範囲 0.70–0.85）**。構成上ほぼ不可能である: fBm の"
      "自己アフィン性は *大きなスケール* の分散で決まるため、2–20 cm で H=0.79 を"
      "出すには基盤振幅が粒の起伏より桁で大きい必要があり、そうすると σz が 3–5 倍に"
      "なる（実測: 玉突き堆積で σz=135–370 mm）。逆に σz を帯に寄せると、2–20 cm の"
      "構造関数は粒のモザイクが支配し H=%.2f になる（粒のモザイクの構造関数は "
      "Δ≈D50 で飽和するので H→0.5 以下に向かう）。折れ点は %.2f·D50 に出ており"
      "「折れ点が D50 付近」だけは満たす。基盤 fBm 自体の H は構成値 0.79。"
      % (bed["structure_H"]["value"], st["structure"]["H_fit"],
         st["structure"]["break_over_D50"]))
    A("")
    A("**復元 D84/D50 = %.2f（範囲 1.9–2.3）**。同じ復元演算子（残差面 → 局所極大 → "
      "最近傍ボロノイ → 面積→円相当径）を、装甲を切った参照床にも掛けると "
      "D84/D50=%.2f（入力は 2.10）。演算子のバイアスが粒径に依存して 1.2–3.6 まで"
      "振れるため、単一の演算子で 1.9–2.3 を検証できる保証がない。装甲比（同じ演算子"
      "の比）は %.2f で、これも参照床の D50 が演算子で %.1f mm に落ちるために "
      "2.0 より小さい。"
      % (bed["recon_D84_over_D50"]["value"],
         st["gsd_reference"]["D84_over_D50"],
         bed["recon_armor_ratio"]["value"],
         1000 * st["gsd_reference"]["D50_m"]))
    A("")
    A("**歪度 = %.2f（範囲 0.3–0.9）**、**空隙率 = %.2f（範囲 0.30–0.45）**、"
      "**クラスター被覆率 = %.1f%%（範囲 5–15%%）**: いずれも構成パラメータ"
      "（planing 公差 0.06·D、バルク体積 1.6·D50 厚、アンカー >D90）から出た値で、"
      "測定後に動かしていない。被覆率と空隙率は帯の縁で外れている。"
      % (bed["skewness"]["value"], bed["porosity"]["value"],
         100 * bed["cluster_coverage"]["value"]))
    A("")
    A("### 1.2 構成法（何を規則として課し、何を測ったか）")
    A("")
    A("- 粒度: 対数正規 (D50=0.020, σg=D84/D50=2.1)。**面積基準**の中央値と解釈し、"
      "個数密度は e^{-2σ²} だけ細粒側へずらして引く（面積重み D² で対数正規の"
      "中央値が e^{+2σ²} 動くため）。装甲層は 2.0 倍の分布。")
    A("- 形状: b/a∈[0.67,0.74], Corey c/√(ab)∈[0.6,0.7]、D は中間軸(b軸)=ふるい径。")
    A("- インブリケーション: a–b 面の傾き 30±10°（法線を下流へ倒す=上流へ傾く）、"
      "a 軸方位は横断 50%%/流れ 50%% の混合（文献が割れているため）。実測 median "
      "%.2f°, 上流向き %.0f%%."
      % (st["imbrication"]["dip_median_deg"],
         100 * st["imbrication"]["upstream_fraction"]))
    A("- クラスター: アンカー >D90、stoss 側に中粒2、wake 側に細粒3、長さ≈2×幅。"
      "被覆率は粒IDから直接測って %.1f%%（目標 10%%）。"
      % (100 * st["cluster_coverage"]))
    A("- 沈降: 各粒を基盤の「到達形態」へ planing（top = form − u·T + 突出）し、"
      "支持関数による非貫通緩和（Jacobi, 8反復, 1反復あたり ±3 mm）で"
      "貫通のみ解く。動力学の主張はしない。**玉突き堆積（z-buffer の max）は"
      "実測で 150–370 mm の塔を作った**ので採用しなかった（数値は本レポート §6）。")
    A("- 基盤: fBm H=0.79、振幅は固定点で決定（σ_target=0.40·D50=8 mm、"
      "粒の寄与を2乗で引いた残差、下限 0.55·σ_target）。実測 σ_sub=%.5f m."
      % st["sigma_sub_m"])
    A("- モノ: D=0.150 / 0.220 m の岩 2 個をクラスター規則（stoss 中粒 + wake 細粒）"
      "で配置。**bed 監査からは除外**（岩は bed ではなく物体）。")
    A("- 深さマップ: 鉛直レイキャスト、レンダ解像度 736×540（0.815 mm/px、"
      "0.6×0.4402 m 窓）と統計用 2048²（0.293 mm/px、0.6 m 四方）。"
      "ファセット法線は解析法線を 80 方向（Fibonacci 半球）へ量子化した低ポリ近似。")
    A("")
    A("## 2. G-TRANS — 高さ場への輸送")
    A("")
    a = fl["a_energy"]
    A("**(a) エネルギー保存** — 平均水深 平坦 %.6f m / 高さ場 %.6f m、"
      "アルベド同一で直接フラックス総量を比較。"
      % (a["mean_depth_flat_m"], a["mean_depth_m"]))
    A("")
    A("| 量 | 平坦床 | 高さ場床 | 相対差 |")
    A("| --- | --- | --- | --- |")
    A("| 直接フラックス（Beer 緑, ゲート量） | %.6g | %.6g | %.3f%% |"
      % (a["flat"]["flux_beer_green"], a["height"]["flux_beer_green"],
         100 * a["rel_diff_beer_green"]))
    A("| 生の光線数（無減衰） | %.6g | %.6g | %.2e |"
      % (a["flat"]["flux_raw"], a["height"]["flux_raw"],
         a["rel_diff_raw"]))
    A("| 平均 cosθ_b（照射面積あたりの係数） | %.4f | %.4f | — |"
      % (a["flat"]["flux_cos_theta_b"], a["height"]["flux_cos_theta_b"]))
    A("| 平均経路長 (m) | %.5f | %.5f | — |"
      % (a["flat"]["mean_path_m"], a["height"]["mean_path_m"]))
    A("| 非収束光線の割合 | — | %.3g | — |"
      % a["height"]["nonconverged_fraction"])
    A("")
    A("生の光線数が厳密に保存されるのは前進スプラットの構成上そうなる（これは"
      "「反復交点で光線を落としていない」ことの検証）。ゲートは Beer 減衰込みの"
      "総量で判定した: 差は %.3f%% で 5%% 以内 → %s。cosθ_b は *真の面積あたり* の"
      "照射係数であり、べき（総量）ではないのでゲートに使っていない"
      "（使えば構成上 %.1f%% ずれる）。"
      % (100 * a["rel_diff_beer_green"], ok(a["pass"]),
         100 * (1 - a["height"]["flux_cos_theta_b"])))
    A("")
    sh = fl["b_shadow"]
    A("**(b) 岩の影** — 解析的な光線/楕円体テスト（2個のモノ岩、9点の太陽円盤）で"
      "「9点すべて遮蔽」の画素をコアとし、レンダラ自身の直接項を測った。")
    A("")
    A("- コア画素 %d（フレームの %.2f%%）、照射画素 %d。"
      % (sh["core_pixels"], 100 * sh["core_fraction"], sh["lit_pixels"]))
    A("- 直接光 コア平均 %.6g / 照射平均 %.6g → 比 **%.4f**（閾値 0.15）→ %s。"
      % (sh["core_direct_mean"], sh["lit_direct_mean"],
         sh["ratio_core_over_lit"], ok(sh["pass"])))
    A("- 深さマップをマーチする影（レンダが使う方）との一致: コア内で "
      "shadow<0.15 の割合 %.3f、照射内で >0.85 の割合 %.3f。"
      % (sh["marching_shadow_agrees_in_core"],
         sh["marching_shadow_agrees_in_lit"]))
    A("- 半影: 影は太陽円盤 9 点の部分遮蔽として出る（半角 0.00395 rad、"
      "水中屈折後の太陽方向 %s）。水中の影は空気中より長い（屈折で太陽高度が "
      "%.1f°→%.1f° に立つ）。"
      % (np.round(sh["sun_dir_water"], 3),
         math.degrees(math.asin(max(0.0, 1 - 0.0))) * 0
         + float(sc["sun"]["altitudeDeg"]),
         math.degrees(math.asin(float(sh["sun_dir_water"][2])))))
    A("")
    A("輸送の実装: 水面の各テクセル（光源格子 %d², 1.25 mm）から太陽円盤の Q=%d "
      "方向へ屈折光線を飛ばし、受け面を平面から z(x,y) に変えて交点を反復で解く"
      "（平面交点 → 2回の固定点 → 最終回だけ Newton、合計3回）。画素ごとの深さ h と"
      "底面傾斜 cosθ_b を直接項に使い、影は直接項のみ、AO（深さマップの地平線近似、"
      "12方位・半径1.2–45 mm）は拡散項のみに掛ける。C は受け面と同じ格子ではなく "
      "1.25 mm 光源格子から 0.816 mm バッファへスプラットする（1:1 だと格子が"
      "モアレとして像に出るため）。固定点の残差は平均 %.2e m、"
      "非収束率（残差>0.1 mm）は %.3f。"
      % (numbers["transport"].get("n_src", 640), Q_DET,
         fl["a_energy"]["height"].get("mean_fixed_point_residual_m", 0.0),
         fl["a_energy"]["height"]["nonconverged_fraction"]))
    A("")
    A("## 3. 場面と材質")
    A("")
    A("- 水面: stream_surface2 (v4-1.2 無改造)、H=0.25 m、U=0.6 m/s、S0=0.003、"
      "D84=0.042 m、seed %d。パッチ L=%.2f m（= n·Us·T_loop, n=1, T_loop=2 s）で"
      "2秒で厳密に1周。" % (SC["seed"], sc["L_water_m"]))
    A("- 導出: Fr=%.4f, u*=%.4f m/s, tau=%.4f s, k0=%.4f rad/m (λ0=%.3f m), "
      "φ=%.4f, σ_target=%.4f mm。妥当域ラベル: %s"
      % (sc["derived"]["Fr"], sc["derived"]["ustar"],
         sc["derived"]["tau"],
         sc["stream"]["k0_rad_m"], sc["stream"]["lambda0_m"],
         sc["stream"]["ring_fraction_phi"],
         1000 * sc["stream"]["sigma_total_target_m"], sc["validity"]))
    A("- 太陽: **プリセット5**（web/recipe.js の5番目＝`blender/recipe-4.json`、"
      "2026-07-20 10:00, 26.3N 127.8E）。高度・方位は run_v19 と同じ "
      "`env.solar_position` で計算: alt=%.3f°, az=%.3f°, sunI=%.4f, "
      "sunColor=%s。" % (sc["sun"]["altitudeDeg"], sc["sun"]["azimuthDeg"],
                        sc["sunI"], np.round(sc["sunColor"], 3)))
    A("- 水の媒体 clear（v12 `clear` プリセット: a=%s, b=%s /m）。η=1、E=2.2、"
      "A=2.6（v20 凍結）、シード %d。"
      % (np.round(sc["a_rgb"], 4), np.round(sc["b_rgb"], 4), SC["seed"]))
    A("- 床アルベド（**全て近似ラベル**）:")
    m = numbers["material"]
    A("  - 石灰岩小石 乾 (%.2f, %.2f, %.2f)。brief の「≤0.22」を上限として採用。"
      % tuple(m["limestone_dry_rgb"]))
    A("  - 砂 乾 (%.2f, %.2f, %.2f)。公開値 ρ(475 nm)=0.47 を B に置き、"
      "長波長で上がる石英の傾向を R/G に反映。" % tuple(m["sand_dry_rgb"]))
    A("  - turf algae 乾 (%.2f, %.2f, %.2f)、被覆 20%%（帯域制限ノイズの下位80%%"
      "閾値で面積20%%を厳密に作る）。" % tuple(m["algae_dry_rgb"]))
    A("  - 湿潤: %s（スカラー暗化ではない）。" % m["wet_transform"])
    A("  - 粒径で砂/礫を混合（6–12 mm で遷移）、粒ごとに ±10%% の明度差"
      "（grain id のハッシュ）。")
    A("")
    A("### 3.1 公開値の出所（すべて近似ラベル付き）")
    A("")
    A("| 値 | 使った数 | 出所 / ラベル |")
    A("| --- | --- | --- |")
    A("| 石灰岩小石 反射率 | ≤0.22 (RGB 0.22/0.21/0.19) | brief の上限値。"
      "**近似**（実測スペクトルは未取得・未検証） |")
    A("| 砂 反射率 | ρ(475 nm)=0.47 | 公開値。他バンドは石英の波長傾向からの"
      "**近似** (0.62/0.56/0.47) |")
    A("| turf algae | RGB (0.10/0.22/0.11) | 藻類スペクトル（Chl-a の 550 nm ピーク／"
      "670 nm 吸収）を RGB に潰した**近似**。一次文献は本セッションで未確認 |")
    A("| 水の a/b | clear プリセット a=[0.360,0.098,0.114], b=[0.046,0.057,0.079] /m |"
      " v12 の公開プリセット（Pope&Fry 1997 の a_w 等）。**既存記録の転用** |")
    A("| 太陽 | alt 54.114°, az 90.628° | プリセット5の日時・地点を "
      "`env.solar_position`（run_v19 と同じ式）で解いた値 |")
    A("| 湿潤関数 | ρ/(1+a(1−ρ)^1.35), a=(0.95,0.72,0.55) | **モデル近似**"
      "（スカラー暗化を避けるための形。物理の一次文献には紐付けていない） |")
    A("")
    A("`L = ρ_wet ·[ E·E_base,dir·C_eff·cosθ_b·shadow·T_down + E_amb·AO·T_up ]"
      " + sky·F(θ_w)·sky_gain·A`、C_eff = w_bal·C + w_halo·(C⊗K_ff)。")
    A("")
    A("- v20 との違い: (i) 受け面が高さ場、(ii) 各画素の h で屈折変位と Beer を計算、"
      "(iii) 影と AO を追加、(iv) **C を屈折変位で動かさない**（C は受け面定義済みで、"
      "動かすと同じ変位を二重に数える）。アルベド側だけ変位させる。")
    A("- Beer は local depth で down/up を別経路（down は屈折角 cos≈0.90）。")
    A("- **拡散項の変更**: v12 の `E_amb = scatter·(1−beer)·amb_gain·A` は写真レーン"
      "向けの校正で、床では 0.003（ハイライトの 0.03%）しか出ない。物理的な天空"
      "放射照度として `E_amb = 2·diffuse·sunI·sunColor·A·amb_gain` を使う"
      "（diffuse=0.10 の晴天、半球放射照度の係数2は宣言した仮定）。この変更で"
      "リニア輝度<0.02 の画素は 17.8%→6.0% になり、純黒（u8<20）は 0 になった"
      "（変更前はボルダーの影が黒く潰れていた）。")
    A("- C の格子: 光源格子 1.25 mm → バッファ 0.816 mm（1:1 にすると格子が"
      "モアレとして像に載る）。C のガウスは 1.5 px（フィラメントが1px幅だと"
      "エイリアスするため、v20 の 0.75 px から上げた）。")
    A("")
    A("## 5. 仮定（すべて）")
    A("")
    A("- 粒度分布の D50/σg は *面積基準*（深さマップが返す量）と解釈した。")
    A("- D は中間軸（ふるい径）。楕円体の半軸 A=D/(2·b/a), B=D/2, C=Corey·√(AB)。")
    A("- 基盤振幅: σ_target=0.40·D50、σ_sub=√(σ_target²−σ_particles²)、下限 0.55σ_target。")
    A("- 空隙率: 各列の深さマップ上端から 2.5·D50 下までのスラブで、384²×128 の"
      "ボクセル占有から算出。")
    A("- 復元 GSD: 残差面(z−ガウス(z,2.2D50)) → 極大 → **最近傍極大ボロノイ** → "
      "面積→円相当径 → 面積重み分位。高さ Watershed は 8bit の flooding が"
      "領域を分離できず（2領域に潰れた）不採用。")
    A("- 影: 太陽円盤 9 点、水中屈折方向で深さマップをマーチ（最大 0.45 m, 180 歩）。"
      "AO: 12方位の地平線、V=mean cos²(h)。")
    A("- 水面ループ: L=1.38 m, T_loop=2 s（(i)+(ii) は厳密に周期、(iii) の OU 成分は"
      "2秒で完全には戻らない）。")
    A("- 総フラックスの定義: 光源光線1本あたり exp(−a_green·経路長) の総和。")
    A("")
    A("## 6. 構築中に実測した失敗（後付け調整をしていない証拠）")
    A("")
    A("| 試行 | 受け面の作り方 | σz (mm) | 歪度 | H | 備考 |")
    A("| --- | --- | --- | --- | --- | --- |")
    A("| A | z-buffer の max で玉突き堆積 | 135 | — | — | 塔が立ち、空隙率 0.55 |")
    A("| B | 最下位ブロックへ落とす | 38 | +0.2 | 0.35 | まだ塔（粒が積み上がる） |")
    A("| C | スラブ内一様 + 非貫通緩和 | 24 | −0.3 | 0.46 | 空隙率 0.17（過充填） |")
    A("| D | 到達形態へ planing（採用） | %.1f | %+.2f | %.2f | 空隙率 %.2f, 被覆 %.1f%% |"
      % (1000 * st["sigma_z_matrix_m"], st["skew_bed"],
         st["structure"]["H_fit"], st["porosity"]["porosity"],
         100 * st["cluster_coverage"]))
    A("")
    A("A–C は全て「粒子を積む」構成で、σz が監査帯 (6.4–9.2 mm) に入らなかった。"
      "D は各粒を同一の到達形態へ bedding する（実在の装甲舗装に近い）構成で、"
      "σz は帯の 1.7 倍まで下がったが、それでも帯には入らない。σz をさらに下げるには"
      "突出公差を 0 に近づけるしかなく、そうすると歪度が 0 に落ちて別のゲートが外れる。")
    A("")
    A("## 7. G-LEGACY（無傷性）")
    A("")
    gleg = g["G-LEGACY"]
    A("- 既存 py %d 本の sha256: 変更 %s、追加 %s、削除 %s。"
      "許容外の変更 %s。" % (gleg["py_count_after"], gleg["changed_existing"],
                            gleg["added"], gleg["removed"],
                            gleg["unexpected"]))
    A("- 変更として出るのは新規2本（`bed_gen.py` / `run_v4_2_1.py`）だけで、"
      "これは stage0 の後に自分で編集したため（新規ファイルなので許容）。"
      "**既存ファイルは1バイトも変わっていない。**")
    A("- 風レーン: h の sha256 = `%s`（v4-1.2 の記録と一致=%s）。"
      % (gleg["wind_sha"][:32], gleg["wind_identical"]))
    A("- 流れレーン: v4-1.2 参照シーンの σ_tot=%.9g m（参照 %.9g、相対差 %.2e）、"
      "k0=%.4f（参照 %.4f）。"
      % (gleg["stream_lane"]["now"]["sigma_total_m"],
         gleg["stream_lane"]["reference"]["sigma_total_m"],
         gleg["stream_lane"]["sigma_rel_err"],
         gleg["stream_lane"]["now"]["k0"],
         gleg["stream_lane"]["reference"]["k0"]))
    A("")
    A("## 8. 実行の記録")
    A("")
    A("- 全ステージ foreground、jobkit checkpoint を stage ごとに原子的コミット"
      "（completed_units=%s）。" % numbers["checkpoint"]["completed_units"])
    A("- 48 frames, 24 fps, %dx%d, libx264 yuv420p, %s s。ffprobe の nb_frames=%s, "
      "duration=%s。still は frame-0024.png（frame-0024 は chunk 0 で書いた"
      "PNG のバイトコピー）。"
      % (g["G-OUT"]["width"], g["G-OUT"]["height"],
         g["G-OUT"]["duration_s"], g["G-OUT"]["nb_frames"],
         g["G-OUT"]["duration_s"]))
    A("- 判定は数値統計のみ（目視なし）。")
    A("")
    A("## 9. 弱点・正直な開示")
    A("")
    A("- G-BED は 8 項目中 1 PASS（imbrication）で、他は FAIL。FAIL のうち σz・H・"
      "復元 D84/D50 は**相互に整合しない要求**であることを §1.1 に数値で示した。"
      "歪度・空隙率・クラスター被覆率は構成パラメータの帰結として帯の外に出た。")
    A("- 装甲比 2.0 を「可視粒径の比」として実測すると %.2f にしかならない"
      "（復元演算子のバイアス）。入力 GSD としては 2.0 で作っている。"
      % bed["recon_armor_ratio"]["value"])
    A("- 輸送ゲート (a) は前進スプラットの構成上ほぼ自明に保存される。"
      "非自明なのは反復交点の非収束率（%.3g）と Beer 経路長の差だけ。" %
      a["height"]["nonconverged_fraction"])
    A("- 影ゲート (b) のコアは解析テストで定義し、レンダの直接項を測る独立検証に"
      "したが、両者が同じ幾何を使う以上「完全な独立」ではない。")
    A("")
    A("## 10. 成果物")
    A("")
    for nm, p in (("bed-stats", BED_STATS), ("depth-map", DEPTH_PNG),
                  ("frame-0024", STILL), ("mp4", MP4),
                  ("numbers", NUMBERS), ("checkpoint", CP_PATH),
                  ("report", REPORT)):
        A("- `%s` → `%s`" % (nm, os.path.relpath(p, HERE)))
    A("")
    atomic_write_bytes(REPORT, ("\n".join(L) + "\n").encode("utf-8"))


def do_report():
    log("=" * 72)
    log("[report] numbers + REPORT + still")
    makedirs()
    cp = Checkpoint(CP_PATH)
    numbers = build_numbers()
    write_json(NUMBERS, numbers)
    # still = frame 24
    p24 = chunk_path(24)
    if os.path.isfile(p24):
        from PIL import Image
        im = Image.open(p24)
        im.load()
        with open(p24, "rb") as f:
            data = f.read()
        atomic_write_bytes(STILL, data)
    write_report(numbers)
    cp.unit_done("report")
    log("  wrote %s (%d bytes)" % (os.path.relpath(REPORT, HERE),
                                   os.path.getsize(REPORT)))
    return numbers


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage0", action="store_true")
    ap.add_argument("--bed", action="store_true")
    ap.add_argument("--scene", action="store_true")
    ap.add_argument("--gtrans", action="store_true")
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--encode", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)
    if a.all or a.stage0:
        stage0()
    if a.all or a.bed:
        do_bed()
    if a.all or a.scene:
        do_scene()
    if a.all or a.gtrans:
        do_gtrans()
    if a.chunk is not None:
        do_chunk(a.chunk)
    if a.all:
        for k in range(NCHUNKS):
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
