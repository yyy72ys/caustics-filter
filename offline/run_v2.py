"""Pipeline v2 runner (PIPELINE-V2, round 2).

Scope: preset 1, representative time t=0 (no full-frame production).
Products in <out>/ (default offline/out/v2/):
  C-t0.npy, flow-t0.npy, metadata.json, numbers-v2.json,
  preset1-still-v2.png

Gate tests recorded into numbers-v2.json:
  (a) angular-quadrature convergence Q=9 vs 25 (+49 if needed),
  (b) photon check (existing transport.py field -> same C definition,
      MS-SSIM of log(1+gamma*C), provisional gate >= 0.9),
  (c) resolution verdict (2048^2 area-mean downsample to 1024^2),
  (d) statistics of C (+ baseline-condition acceptance stats, OUTPUT-SPEC 7).
"""
import argparse
import datetime
import json
import os
import sys
import time

import numpy as np
from PIL import Image
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment, sun_color_from_altitude  # noqa: E402
from state import WaterState, derive_seed  # noqa: E402
from transport import photon_frame, sun_disc_radius  # noqa: E402
from det import det_frame, flow_backward_master, IOR_G  # noqa: E402
from image import (render_frame_v2, e_base_dir_from_env,  # noqa: E402
                   srgb_to_linear)

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)
GAMMA_DEFAULT = 1.0
MS_SSIM_WEIGHTS = (0.0448, 0.2856, 0.3001, 0.2363, 0.1333)


def to_lum(a):
    a = np.asarray(a, dtype=np.float64)
    return (a * REC709).sum(axis=-1)


def quantiles(a):
    a = np.asarray(a, dtype=np.float64).ravel()
    return {
        "mean": float(np.mean(a)),
        "p50": float(np.percentile(a, 50)),
        "p99": float(np.percentile(a, 99)),
        "p99.9": float(np.percentile(a, 99.9)),
        "max": float(np.max(a)),
        "min": float(np.min(a)),
    }


def log_rmse_rel(a, b):
    """RMSE(ln a - ln b) / RMS(ln b); a,b > 0 expected (guarded)."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    # NOTE: log only; the master C itself is never clamped. Tiny FFT-ring
    # negatives (<1e-9) are clipped for the logarithm domain.
    la = np.log(np.clip(a, 1e-12, None))
    lb = np.log(np.clip(b, 1e-12, None))
    rmse = float(np.sqrt(np.mean((la - lb) ** 2)))
    rms = float(np.sqrt(np.mean(lb ** 2)))
    med = float(np.median(np.abs(la - lb)))
    return {"rmse": rmse, "rms_ref": rms, "rel": rmse / rms if rms else 0.0,
            "median_abs": med}


def spectral_band_frac(f):
    """Fractional power of demeaned field f in 5 log k-bands (cycles/m)."""
    f = np.asarray(f, dtype=np.float64)
    N = f.shape[0]
    S = float(spectral_band_frac.S)
    F = np.fft.rfft2(f - f.mean())
    P = (F.real ** 2 + F.imag ** 2)
    fy = np.fft.fftfreq(N, d=S / N)[:, None]
    fx = np.fft.rfftfreq(N, d=S / N)[None, :]
    k = np.sqrt(fx ** 2 + fy ** 2)
    kmin, kmax = 1.0 / S, (N // 2) / S
    edges = np.logspace(np.log10(kmin), np.log10(kmax), 6)
    tot = float(P.sum())
    fracs = []
    for i in range(5):
        m = (k >= edges[i]) & (k < edges[i + 1] if i < 4 else k <= edges[i + 1])
        fracs.append(float(P[m].sum() / tot) if tot else 0.0)
    return {"edges_cycles_per_m": [float(e) for e in edges], "frac": fracs}


def spectral_reldiff(fa, fb):
    sa, sb = spectral_band_frac(fa), spectral_band_frac(fb)
    rel = []
    for x, y in zip(sa["frac"], sb["frac"]):
        rel.append(abs(x - y) / y if y else 0.0)
    return {"frac_a": sa["frac"], "frac_b": sb["frac"], "rel_per_band": rel,
            "edges": sa["edges_cycles_per_m"]}


def _ssim_cs(a, b, L):
    C1 = (0.01 * L) ** 2
    C2 = (0.03 * L) ** 2
    ma = ndimage.gaussian_filter(a, 1.5, mode="reflect")
    mb = ndimage.gaussian_filter(b, 1.5, mode="reflect")
    ma2 = ndimage.gaussian_filter(a * a, 1.5, mode="reflect")
    mb2 = ndimage.gaussian_filter(b * b, 1.5, mode="reflect")
    mab = ndimage.gaussian_filter(a * b, 1.5, mode="reflect")
    va = ma2 - ma * ma
    vb = mb2 - mb * mb
    vc = mab - ma * mb
    l = (2 * ma * mb + C1) / (ma * ma + mb * mb + C1)
    cs = (2 * vc + C2) / (va + vb + C2)
    return l, cs


def ms_ssim(a, b):
    """MS-SSIM (5 scales, Wang et al. weights) on float64 maps."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    L = float(max(a.max(), b.max()) - min(a.min(), b.min()))
    m = 1.0
    parts = []
    for i in range(5):
        l, cs = _ssim_cs(a, b, L)
        w = MS_SSIM_WEIGHTS[i]
        if i < 4:
            m *= float(np.mean(np.clip(cs, 0.0, None)) ** w)
            parts.append({"scale": i, "cs_mean": float(np.mean(cs))})
        else:
            m *= float(np.mean(np.clip(cs, 0.0, None)) ** w)
            m *= float(np.mean(np.clip(l, 0.0, None)) ** w)
            parts.append({"scale": i, "cs_mean": float(np.mean(cs)),
                          "l_mean": float(np.mean(l))})
        if i < 4:
            n0 = (a.shape[0] // 2) * 2
            n1 = (a.shape[1] // 2) * 2
            a = 0.25 * (a[:n0:2, :n1:2] + a[1:n0:2, :n1:2]
                        + a[:n0:2, 1:n1:2] + a[1:n0:2, 1:n1:2])
            b = 0.25 * (b[:n0:2, :n1:2] + b[1:n0:2, :n1:2]
                        + b[:n0:2, 1:n1:2] + b[1:n0:2, 1:n1:2])
    return {"msssim": float(m), "L": L, "weights": list(MS_SSIM_WEIGHTS),
            "scales": parts}


def ac_fwhm_px(f):
    """FWHM (px) of the normalized autocorrelation central peak."""
    f = np.asarray(f, dtype=np.float64)
    F = np.fft.fft2(f - f.mean())
    pw = (F.real ** 2 + F.imag ** 2)
    denom = float(pw.sum())
    if not np.isfinite(denom) or denom <= 0.0:
        return {"fwhm_px": float("nan"), "peak": float("nan"),
                "note": "degenerate (constant/non-finite) field"}
    ac = np.fft.fftshift(np.fft.ifft2(pw).real) / denom * pw.size
    N = f.shape[0]
    cy, cx = N // 2, N // 2
    prof = []
    rmax = N // 2
    yy, xx = np.mgrid[0:N, 0:N]
    rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    nb = min(256, rmax)
    edges = np.linspace(0.0, float(rmax), nb + 1)
    for i in range(nb):
        m = (rr >= edges[i]) & (rr < edges[i + 1])
        prof.append(float(ac[m].mean()) if m.any() else 0.0)
    prof = np.array(prof)
    half = 0.5 * prof[0]
    below = np.where(prof < half)[0]
    if len(below) == 0:
        return {"fwhm_px": float("nan"), "peak": float(prof[0])}
    j = int(below[0])
    if j == 0:
        r = 0.0
    else:
        t = (half - prof[j - 1]) / (prof[j] - prof[j - 1] + 1e-30)
        r = edges[j - 1] + t * (edges[j] - edges[j - 1])
    return {"fwhm_px": float(2 * r), "peak": float(prof[0])}


def count_components(mask):
    lab, n = ndimage.label(np.asarray(mask, dtype=bool),
                           structure=np.ones((3, 3), dtype=int))
    return int(n)


def area_mean_downsample(f, factor=2):
    f = np.asarray(f)
    H, W = f.shape[:2]
    H2, W2 = (H // factor) * factor, (W // factor) * factor
    f = f[:H2, :W2]
    if f.ndim == 2:
        return f.reshape(H2 // factor, factor, W2 // factor, factor).mean(axis=(1, 3))
    return f.reshape(H2 // factor, factor, W2 // factor, factor, -1).mean(axis=(1, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", default="../blender/recipe-1.json")
    ap.add_argument("--preset", type=int, default=1)
    ap.add_argument("--grid", type=int, default=2048)
    ap.add_argument("--photons", type=int, default=4000000)
    ap.add_argument("--qlevels", default="9,25")
    ap.add_argument("--gamma", type=float, default=GAMMA_DEFAULT)
    ap.add_argument("--out", default="out/v2")
    ap.add_argument("--photo", default="../test-inputs/moth-on-bark.jpg")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--no-baseline", action="store_true")
    args = ap.parse_args()

    t_all = time.time()
    os.makedirs(args.out, exist_ok=True)
    with open(args.recipe, "r", encoding="utf-8") as f:
        recipe = json.load(f)
    if recipe.get("schema") != "caustics-recipe/v1":
        raise SystemExit("recipe schema mismatch")
    env = resolve_environment(recipe)

    t0 = time.time()
    state = WaterState(recipe, env, N=args.grid, seed=args.seed,
                       preset_index=args.preset)
    t_state = time.time() - t0

    fields = state.fields_at(0.0)
    fields["S"] = state.S
    spectral_band_frac.S = state.S

    Ebase = e_base_dir_from_env(env)
    L_base = float((Ebase * REC709).sum())

    # ---- deterministic transport at requested quadrature levels ----
    qlevels = [int(x) for x in str(args.qlevels).split(",") if x.strip()]
    C_levels, det_infos, t_dets = {}, {}, {}
    for q in qlevels:
        tq = time.time()
        Cq, info = det_frame(fields, state.S, args.grid, env, q)
        t_dets[str(q)] = time.time() - tq
        C_levels[q] = Cq
        det_infos[str(q)] = info
        print("det Q=%d done mean=%.6f max=%.3f (%.1fs)"
              % (q, float(Cq.mean()), float(Cq.max()), t_dets[str(q)]),
              flush=True)

    def quad_compare(qa, qb):
        A, B = C_levels[qa], C_levels[qb]
        qa_s, qb_s = quantiles(A), quantiles(B)
        rel = {}
        for k in ("p50", "p99", "p99.9"):
            denom = abs(qb_s[k]) if qb_s[k] else 0.0
            rel[k] = abs(qa_s[k] - qb_s[k]) / denom if denom else 0.0
        lr = log_rmse_rel(A, B)
        sp = spectral_reldiff(A, B)
        return {"Qa": qa, "Qb": qb, "quantiles_a": qa_s, "quantiles_b": qb_s,
                "quantile_rel": rel, "log_rmse": lr, "spectral": sp}

    comparisons = {}
    for i in range(1, len(qlevels)):
        comparisons["Q%d_vs_Q%d" % (qlevels[i - 1], qlevels[i])] = \
            quad_compare(qlevels[i - 1], qlevels[i])

    TARGETS = {"log_rmse_rel": 0.01, "quantile_rel": 0.01,
               "spectral_band_rel": 0.02}

    def quad_pass(cmp):
        ok_log = cmp["log_rmse"]["rel"] <= TARGETS["log_rmse_rel"]
        ok_q = all(v <= TARGETS["quantile_rel"]
                   for v in cmp["quantile_rel"].values())
        ok_sp = all(v <= TARGETS["spectral_band_rel"]
                    for v in cmp["spectral"]["rel_per_band"])
        return {"log": bool(ok_log), "quantiles": bool(ok_q),
                "spectral": bool(ok_sp),
                "all": bool(ok_log and ok_q and ok_sp)}

    passes = {k: quad_pass(v) for k, v in comparisons.items()}

    # extend to Q=49 when the top comparison fails any target
    top_key = "Q%d_vs_Q%d" % (qlevels[-2], qlevels[-1]) if len(qlevels) > 1 else None
    if top_key is not None and not passes[top_key]["all"] and 49 not in qlevels:
        tq = time.time()
        C49, info49 = det_frame(fields, state.S, args.grid, env, 49)
        t_dets["49"] = time.time() - tq
        C_levels[49] = C49
        det_infos["49"] = info49
        print("det Q=49 done mean=%.6f max=%.3f (%.1fs)"
              % (float(C49.mean()), float(C49.max()), t_dets["49"]), flush=True)
        comparisons["Q%d_vs_Q49" % qlevels[-1]] = quad_compare(qlevels[-1], 49)
        passes["Q%d_vs_Q49" % qlevels[-1]] = quad_pass(
            comparisons["Q%d_vs_Q49" % qlevels[-1]])
        qlevels = qlevels + [49]

    # adopted = smallest Q whose comparison against the next level passes;
    # if nothing passes, best (max Q) with gate marked unmet.
    adopted = qlevels[-1]
    for i in range(len(qlevels) - 1):
        key = "Q%d_vs_Q%d" % (qlevels[i], qlevels[i + 1])
        if passes[key]["all"]:
            adopted = qlevels[i]
            break
    C = C_levels[adopted]

    # ---- photon check (existing transport.py, same condition) ----
    tp = time.time()
    E_ph, ph_info = photon_frame(fields, state.S, args.grid, env,
                                 args.photons,
                                 seed=(state.seed + 1000) % (2 ** 63))
    t_photon = time.time() - tp
    C_ph = (to_lum(E_ph) / L_base).astype(np.float64)
    C_det = C.astype(np.float64)
    gam = float(args.gamma)
    X_det, X_ph = np.log1p(gam * C_det), np.log1p(gam * C_ph)
    ms = ms_ssim(X_det, X_ph)
    ph_stats, det_stats = quantiles(C_ph), quantiles(C)
    ph_qrel = {}
    for k in ("p50", "p99", "p99.9"):
        denom = abs(ph_stats[k]) if ph_stats[k] else 0.0
        ph_qrel[k] = abs(det_stats[k] - ph_stats[k]) / denom if denom else 0.0
    photon_check = {
        "photons": args.photons, "gamma": gam,
        "msssim": ms, "gate": 0.9, "pass": bool(ms["msssim"] >= 0.9),
        "quantiles_det": det_stats, "quantiles_photon": ph_stats,
        "quantile_rel_det_vs_photon": ph_qrel,
        "spectral_det_vs_photon": spectral_reldiff(C_det, C_ph),
        "fwhm_det": ac_fwhm_px(C_det), "fwhm_photon": ac_fwhm_px(C_ph),
        "diff_X_mean_abs": float(np.mean(np.abs(X_det - X_ph))),
        "expected_mean_rgb": ph_info["expected_mean"],
        "measured_mean_rgb": ph_info["measured_mean"],
    }
    print("photon check MS-SSIM=%.4f gate=0.9 %s (%.1fs)"
          % (ms["msssim"], "PASS" if photon_check["pass"] else "FAIL", t_photon),
          flush=True)

    # ---- resolution verdict ----
    C1024 = area_mean_downsample(C_det, 2)
    s2048, s1024 = quantiles(C_det), quantiles(C1024)
    flux_rel = abs(s1024["mean"] - s2048["mean"]) / abs(s2048["mean"])
    p999_rel = (abs(s1024["p99.9"] - s2048["p99.9"]) / abs(s2048["p99.9"])
                if s2048["p99.9"] else 0.0)
    thr = float(np.percentile(C_det, 99))
    n2048_same = count_components(C_det > thr)
    n1024_same = count_components(C1024 > thr)
    n2048_own = count_components(C_det > float(np.percentile(C_det, 99)))
    n1024_own = count_components(C1024 > float(np.percentile(C1024, 99)))
    res = {"flux_rel_diff": float(flux_rel), "flux_gate": 0.001,
           "flux_pass": bool(flux_rel <= 0.001),
           "p999_rel_diff": float(p999_rel), "p999_gate": 0.05,
           "p999_pass": bool(p999_rel <= 0.05),
           "threshold_p99_det": thr,
           "components_2048_same_thr": n2048_same,
           "components_1024_same_thr": n1024_same,
           "components_2048_own_p99": n2048_own,
           "components_1024_own_p99": n1024_own,
           "topology_invariant": bool(n2048_same == n1024_same)}
    res["adopt_1024"] = bool(res["flux_pass"] and res["p999_pass"]
                             and res["topology_invariant"])
    res["master_resolution"] = 1024 if res["adopt_1024"] else args.grid
    print("resolution: flux=%.2e p999=%.3f topo %d->%d adopt=%s"
          % (flux_rel, p999_rel, n2048_same, n1024_same,
             res["master_resolution"]), flush=True)

    # ---- baseline-condition acceptance stats (OUTPUT-SPEC 7) ----
    baseline = None
    t_base = 0.0
    if not args.no_baseline and args.grid >= 1024:
        tb = time.time()
        b_recipe = json.loads(json.dumps(recipe))
        b_recipe["water"]["depthM"] = 1.0
        b_recipe["windWave"]["windSpeedMS"] = 3.0
        b_env = resolve_environment(b_recipe)
        b_env["specL"] = 2.0  # domain L = S = 8.0 m
        az = np.deg2rad(b_env["sun"]["azimuthDeg"])
        alt = np.deg2rad(45.0)
        sd = [np.sin(az) * np.cos(alt), np.cos(az) * np.cos(alt), np.sin(alt)]
        b_env["sunDir"] = sd
        b_env["sun"]["altitudeDeg"] = 45.0
        b_env["sunColor"] = list(sun_color_from_altitude(45.0))
        b_env["clarity"] = 0.55  # disc ~0.53 deg + weak diffuse cone
        b_env["depthM"] = 1.0
        b_state = WaterState(b_recipe, b_env, N=args.grid, seed=None,
                             preset_index=args.preset)
        b_fields = b_state.fields_at(0.0)
        b_fields["S"] = b_state.S
        Cb, b_info = det_frame(b_fields, b_state.S, args.grid, b_env, adopted)
        t_base = time.time() - tb
        bs = quantiles(Cb)
        baseline = {
            "condition": {"sun_alt_deg": 45.0, "depthM": 1.0,
                          "windSpeedMS": 3.0, "domain_L_m": 8.0,
                          "clarity": 0.55,
                          "disc_rad": float(sun_disc_radius(0.55))},
            "seed": int(b_state.seed), "specL_m": 2.0, "S_m": float(b_state.S),
            "Q": adopted, "det_info": b_info,
            "stats": bs, "p999_over_mean": (bs["p99.9"] / bs["mean"]
                                           if bs["mean"] else 0.0),
            "reference_p999_over_mean": [6.0, 12.0],
            "note": "6-12 is a condition-dependent reference, not a hard gate",
        }
        print("baseline Q=%d mean=%.4f p999/mean=%.2f (%.1fs)"
              % (adopted, bs["mean"], baseline["p999_over_mean"], t_base),
              flush=True)

    # ---- flow master (matches C master grid) ----
    flow_full = flow_backward_master(fields, state.S, args.grid, env)
    if res["master_resolution"] != args.grid:
        C_master = C1024.astype(np.float32)
        F_master = area_mean_downsample(flow_full, 2).astype(np.float32)
    else:
        C_master = C.astype(np.float32)
        F_master = flow_full

    # ---- v2 preview still ----
    photo = np.asarray(Image.open(args.photo).convert("RGB"),
                       dtype=np.float64) / 255.0
    photo_lin = srgb_to_linear(photo)
    PH, PW, _ = photo_lin.shape
    timg = time.time()
    still = render_frame_v2(fields, C_master, photo_lin, env, PW, PH)
    t_image = time.time() - timg
    Image.fromarray(still).save(os.path.join(args.out, "preset1-still-v2.png"))

    np.save(os.path.join(args.out, "C-t0.npy"), C_master)
    np.save(os.path.join(args.out, "flow-t0.npy"), F_master)

    s_master = quantiles(C_master)
    numbers = {
        "pipeline": "v2",
        "recipe": recipe.get("name"), "recipe_schema": recipe.get("schema"),
        "preset": args.preset, "seed": int(state.seed), "t_s": 0.0,
        "state": {
            "grid": [args.grid, args.grid], "patch_S_m": state.S,
            "dx_m": state.dx, "kmin_rad_per_m": float(state.kmin),
            "kmax_nyquist_rad_per_m": float(state.kmax_nyq),
            "band_decades": float(state.decades), "specL_m": float(state.specL),
            "choppy_lambda": 0.85, "loop_T_s": 8.0,
            "loop_maxdiff_h_m": float(state.loop_maxdiff()),
        },
        "transport_main": {
            "method": ("deterministic stratified forward mapping "
                       "(Jacobian-equivalent): 1 ray per surface texel per "
                       "solar-disc direction, Snell (G-band IOR %.5f), "
                       "bilinear splat + %.2fpx FFT Gaussian (same pixel "
                       "filter as photon path); all inverse images summed "
                       "structurally; no clamp / no renormalization / no "
                       "soft shoulder" % (IOR_G, 0.75)),
            "Q_levels": qlevels, "Q_adopted": adopted,
            "adopt_rule": ("smallest Q passing all targets vs next level; "
                           "max Q if none pass"),
            "targets": TARGETS,
            "comparisons": comparisons,
            "passes": passes,
            "per_Q": {str(q): {"stats": quantiles(C_levels[q]),
                               "det_info": det_infos[str(q)],
                               "time_s": t_dets[str(q)]} for q in qlevels},
            "E_base_dir_definition": ("flat-water bottom irradiance "
                                      "E_base,dir = sunI * sunColor; "
                                      "C = E_direct / E_base,dir"),
            "E_base_dir_rgb": [float(x) for x in Ebase],
            "E_base_dir_luminance": L_base,
        },
        "photon_check": photon_check,
        "resolution": res,
        "statistics_master_C": {
            **s_master,
            "p999_over_mean": (s_master["p99.9"] / s_master["mean"]
                               if s_master["mean"] else 0.0),
            "flux_ratio_vs_flat": float(s_master["mean"]),
            "spectral_frac": spectral_band_frac(
                np.asarray(C_master, dtype=np.float64)),
        },
        "baseline_condition": baseline,
        "environment": {
            "sun_alt_deg": env["sun"]["altitudeDeg"],
            "sun_az_deg": env["sun"]["azimuthDeg"],
            "sunI": env["sunI"], "sunColor": list(env["sunColor"]),
            "depthM": env["depthM"], "clarity": env.get("clarity"),
            "disc_rad": float(sun_disc_radius(env.get("clarity", 0.8))),
        },
        "timing_s": {"state_build": t_state,
                     **{"det_Q%s" % k: v for k, v in t_dets.items()},
                     "photon": t_photon, "baseline": t_base,
                     "image_v2": t_image, "wall_total": time.time() - t_all},
        "outputs": {
            "C": "C-t0.npy", "flow": "flow-t0.npy",
            "metadata": "metadata.json", "numbers": "numbers-v2.json",
            "still": "preset1-still-v2.png",
            "still_size": [int(PW), int(PH)],
            "C_shape": list(C_master.shape),
            "C_dtype": str(C_master.dtype),
            "flow_shape": list(F_master.shape),
            "flow_dtype": str(F_master.dtype),
        },
    }
    with open(os.path.join(args.out, "numbers-v2.json"), "w",
              encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)

    metadata = {
        "pipeline": "caustics-filter offline v2",
        "created_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "recipe": recipe.get("name"), "schema": recipe.get("schema"),
        "preset": args.preset, "seed": int(state.seed), "t_s": 0.0,
        "grids": {
            "surface_N": [args.grid, args.grid],
            "master_N": [int(res["master_resolution"]),
                         int(res["master_resolution"])],
            "patch_S_m": state.S, "dx_m": state.dx,
            "downsample": ("area mean 2x2" if res["adopt_1024"]
                           else "none (master = 2048^2)"),
        },
        "C": {
            "definition": ("dimensionless direct-sun illuminance ratio "
                           "C = E_direct / E_base,dir"),
            "E_base_dir_definition": ("mean bottom irradiance under "
                                      "identical conditions with a flat "
                                      "mean-plane water surface"),
            "E_base_dir_rgb": [float(x) for x in Ebase],
            "E_base_dir_luminance_rec709": L_base,
            "dtype": "float32", "channels": 1,
            "processing": "none: no clamp, no mean normalization, "
                          "no soft shoulder, no tone map",
            "quadrature": {"type": "deterministic Fibonacci disc",
                           "Q_adopted": adopted, "Q_tested": qlevels,
                           "disc_rad": float(sun_disc_radius(
                               env.get("clarity", 0.8)))},
            "pixel_filter": "bilinear splat + sigma=0.75px FFT Gaussian "
                            "(same as photon path)",
            "ior": {"note": "achromatic representative band", "G": IOR_G},
        },
        "flow_backward": {
            "definition": ("output-pixel -> input-photo reverse displacement "
                           "(G band), master-grid pixels"),
            "dtype": "float32", "channels": 2, "order": "[dx, dy]",
            "convention": ("origin = top-left corner of patch (u=0,v=0); "
                           "pixel centers at (i+0.5, j+0.5); +x = +u (right), "
                           "+y = +v (down, row index increasing); "
                           "sample = (i+0.5+dx, j+0.5+dy) in a master-sized "
                           "photo; forward detJ never mixed in"),
        },
        "software": {"det": "offline/det.py", "image": "offline/image.py "
                     "render_frame_v2 (OUTPUT-SPEC section 3 separation)"},
    }
    with open(os.path.join(args.out, "metadata.json"), "w",
              encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print("OK wall=%.1fs adoptedQ=%d msssim=%.4f master=%d" % (
        time.time() - t_all, adopted, ms["msssim"],
        res["master_resolution"]), flush=True)


if __name__ == "__main__":
    main()
