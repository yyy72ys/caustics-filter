"""v8 prerequisite-task analyses (PIPELINE-V8 "前提タスク" / chunk-4 gate 0).

Three items the rally left unverified, implemented as measurements rather
than opinions:

  E. SMITH SHADOW AUDIT -- read the direct-term path end to end and decide
     whether the deterministic transport (all preimage Jacobians) already
     resolves the same scale of occlusion that the statistical G1 removes.
     The answer is a *verdict plus code citations*, and a numeric test of
     what the term does on the reference still.

  Q.  SOLAR-DISC QUADRATURE CONVERGENCE -- C statistics (p50/p99/p99.9/max),
     mean, peak position and MS-SSIM for Q = 9/25/49.

  M.  MS-SSIM DECOMPOSITION vs the photon check -- split the difference into
     (a) peak-luminance absolute error and (b) peak-position shift, instead
     of reporting a single scalar that hides which one moved.

The MS-SSIM / autocorrelation helpers are copied verbatim from run_v2.py so
the numbers stay comparable with the frozen v2 photon-check convention.
"""
import math

import numpy as np
from scipy import ndimage

MS_SSIM_WEIGHTS = (0.0448, 0.2856, 0.3001, 0.2363, 0.1333)


# --------------------------------------------------------------------------
# shared image metrics (v2 convention)
# --------------------------------------------------------------------------
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
    """MS-SSIM (5 scales, Wang et al. weights) -- identical to run_v2.py."""
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
    """FWHM (px) of the normalised autocorrelation central peak (run_v2.py)."""
    f = np.asarray(f, dtype=np.float64)
    F = np.fft.fft2(f - f.mean())
    pw = (F.real ** 2 + F.imag ** 2)
    denom = float(pw.sum())
    if not np.isfinite(denom) or denom <= 0.0:
        return {"fwhm_px": float("nan"), "peak": float("nan")}
    ac = np.fft.fftshift(np.fft.ifft2(pw).real) / denom * pw.size
    N = f.shape[0]
    cy, cx = N // 2, N // 2
    rmax = N // 2
    yy, xx = np.mgrid[0:N, 0:N]
    rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    nb = min(256, rmax)
    edges = np.linspace(0.0, float(rmax), nb + 1)
    prof = np.zeros(nb)
    for i in range(nb):
        m = (rr >= edges[i]) & (rr < edges[i + 1])
        prof[i] = float(ac[m].mean()) if m.any() else 0.0
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


def quantiles(C):
    c = np.asarray(C, dtype=np.float64)
    return {"mean": float(c.mean()), "p50": float(np.percentile(c, 50)),
            "p99": float(np.percentile(c, 99)),
            "p99.9": float(np.percentile(c, 99.9)),
            "max": float(c.max()), "min": float(c.min())}


# --------------------------------------------------------------------------
# peak readouts / peak shift
# --------------------------------------------------------------------------
def peak_metrics(C, frac=1e-3):
    """Peak luminance and peak position of one C master."""
    c = np.asarray(C, dtype=np.float64)
    n = c.size
    k = max(1, int(round(frac * n)))
    flat = c.ravel()
    idx = np.argpartition(flat, -k)[-k:]
    thr = float(flat[idx].min())
    yy, xx = np.nonzero(c >= thr)
    w = c[yy, xx] - thr
    w = np.where(np.isfinite(w), np.maximum(w, 0.0), 0.0)
    cy = float((yy * w).sum() / w.sum()) if w.sum() > 0 else float(yy.mean())
    cx = float((xx * w).sum() / w.sum()) if w.sum() > 0 else float(xx.mean())
    ay, ax = np.unravel_index(int(np.argmax(c)), c.shape)
    return {
        "max": float(c.max()),
        "p99_9": float(np.percentile(c, 99.9)),
        "top_frac_mean": float(flat[idx].mean()),
        "top_frac": float(frac),
        "threshold": thr,
        "argmax_yx": [int(ay), int(ax)],
        "centroid_yx": [cy, cx],
        "n_pixels_above": int((c >= thr).sum()),
    }


def _xcorr_shift(a, b, half=24):
    """Sub-pixel shift (dy, dx) that best aligns b onto a (near their peaks)."""
    H, W = a.shape
    ay, ax = np.unravel_index(int(np.argmax(a)), a.shape)
    y0 = int(np.clip(ay - 96, 0, H - 1))
    y1 = int(np.clip(ay + 96, 1, H))
    x0 = int(np.clip(ax - 96, 0, W - 1))
    x1 = int(np.clip(ax + 96, 1, W))
    wa = a[y0:y1, x0:x1] - float(a[y0:y1, x0:x1].mean())
    best = None
    for dy in range(-half, half + 1):
        for dx in range(-half, half + 1):
            yb0, yb1 = y0 + dy, y1 + dy
            xb0, xb1 = x0 + dx, x1 + dx
            if yb0 < 0 or xb0 < 0 or yb1 > a.shape[0] or xb1 > a.shape[1]:
                continue
            wb = b[yb0:yb1, xb0:xb1] - float(b[yb0:yb1, xb0:xb1].mean())
            s = float((wa * wb).sum())
            if best is None or s > best[0]:
                best = (s, dy, dx)
    if best is None:
        return {"dy_px": None, "dx_px": None, "dist_px": None}
    return {"dy_px": float(best[1]), "dx_px": float(best[2]),
            "dist_px": float(math.hypot(best[1], best[2])),
            "corr": float(best[0])}


def _smooth(f, sigma=3.0):
    return ndimage.gaussian_filter(np.asarray(f, dtype=np.float64), sigma,
                                   mode="wrap")


def peak_shift(A, B, S_patch, label_a="det", label_b="photon", sigma=3.0):
    """Peak-position difference between two C masters (px and mm).

    Four readouts, because they answer different questions:
      argmax            -- raw single-pixel maximum (noise-sensitive)
      argmax_smoothed    -- maximum after a 3 px Gaussian (the caustic peak)
      centroid           -- mass centroid of the top-0.1% pixels (robust)
      xcorr              -- window cross-correlation shift (sub-pixel integer)
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    pa, pb = peak_metrics(A), peak_metrics(B)
    sa, sb = _smooth(A, sigma), _smooth(B, sigma)
    psa, psb = peak_metrics(sa), peak_metrics(sb)
    dy = pb["argmax_yx"][0] - pa["argmax_yx"][0]
    dx = pb["argmax_yx"][1] - pa["argmax_yx"][1]
    sdy = psb["argmax_yx"][0] - psa["argmax_yx"][0]
    sdx = psb["argmax_yx"][1] - psa["argmax_yx"][1]
    cdy = pb["centroid_yx"][0] - pa["centroid_yx"][0]
    cdx = pb["centroid_yx"][1] - pa["centroid_yx"][1]
    xc = _xcorr_shift(sa, sb)
    mm_per_px = float(S_patch) / A.shape[0]
    return {
        "a": label_a, "b": label_b,
        "argmax_a": pa["argmax_yx"], "argmax_b": pb["argmax_yx"],
        "argmax_shift_px": float(math.hypot(dy, dx)),
        "argmax_shift_mm": float(math.hypot(dy, dx) * mm_per_px),
        "smoothed_sigma_px": float(sigma),
        "argmax_smoothed_a": psa["argmax_yx"],
        "argmax_smoothed_b": psb["argmax_yx"],
        "argmax_smoothed_shift_px": float(math.hypot(sdy, sdx)),
        "argmax_smoothed_shift_mm": float(math.hypot(sdy, sdx) * mm_per_px),
        "centroid_shift_px": float(math.hypot(cdy, cdx)),
        "centroid_shift_mm": float(math.hypot(cdy, cdx) * mm_per_px),
        "xcorr_shift_px": xc.get("dist_px"),
        "xcorr_shift_mm": (None if xc.get("dist_px") is None
                           else float(xc["dist_px"] * mm_per_px)),
        "mm_per_master_px": mm_per_px,
        "note": ("the smoothed/centroid/xcorr readouts are the meaningful "
                 "ones for a noisy field; the raw argmax is reported for "
                 "completeness and is the one a single noisy pixel can own"),
    }


# --------------------------------------------------------------------------
# E. Smith shadow audit
# --------------------------------------------------------------------------
SMITH_CODE_FACTS = [
    {"file": "det.py", "function": "det_frame (lines 60-161)",
     "fact": ("each surface texel emits one ray per solar-disc direction and "
              "is forward-mapped onto the bottom buffer: "
              "xb = mod(Psx + Tx*t, S).  There is NO occlusion test anywhere "
              "-- no ray-ray blocking, no horizon test, no visibility mask."),
     "consequence": ("all preimages are summed structurally (that is the "
                     "Jacobian sum), but a surface point hidden behind a "
                     "crest still contributes its full weight.  Shadowing is "
                     "NOT resolved at any scale.")},
    {"file": "det.py", "function": "energy bookkeeping (lines 145-150)",
     "fact": ("divisor = total_rays / N_buf^2 is analytic; the buffer is "
              "never rescaled.  Summing all rays means the transport "
              "conserves energy exactly, so it cannot by itself produce the "
              "energy loss that shadowing causes."),
     "consequence": ("the resolved-scale energy deficit of a shadowed facet "
                     "is absent from C, i.e. the transport OVER-counts lit "
                     "area relative to a ray-exact shadow test.")},
    {"file": "transport.py", "function": "photon_frame (lines 72-152)",
     "fact": ("the photon path is the same single-bounce height-field map: "
              "one refraction, one bottom intersection, no re-entry test."),
     "consequence": "the two transports share the omission; it is not a det.py bug."},
    {"file": "v5_composite.py", "function": "render_frame_v5 (lines 196-202)",
     "fact": ("direct = E_base,dir * c * beer_down * beer_up * g1 with g1 a "
              "single scalar from smith_g1_beckmann(realised rms slope, sun "
              "zenith).  It multiplies the whole direct term exactly once."),
     "consequence": ("G1 is the only occlusion-like factor in the image, it "
                     "is statistical (Beckmann/Walter fit), and it is not "
                     "derived from the transported field.")},
    {"file": "v5_physics.py", "function": "smith_g1_beckmann (lines 218-242)",
     "fact": ("G1 = 1/(1+Lambda), a = 1/(alpha tan theta); for a >= 1.6 the "
              "fit returns Lambda = 0, G1 = 1 exactly."),
     "consequence": ("the term is a step-free switch to 1.0 in the "
                     "high-sun / low-slope corner, so it is inert there.")},
]

SMITH_VERDICT = {
    "verdict": "KEEP (maintain) with the scale separation stated explicitly",
    "double_counting": False,
    "reason": (
        "The deterministic transport resolves FOCUSING (the Jacobian sum over "
        "all preimages) but contains no occlusion test at all: no ray is ever "
        "blocked, and the energy divisor is analytic, so the transport cannot "
        "remove the light that a crest hides.  The statistical G1 removes "
        "exactly that unresolved masking-shadowing.  The two therefore act on "
        "disjoint physics (magnification vs visibility) and are not "
        "double-counted."),
    "why_not_remove": (
        "Removing G1 (Gemini's 'remove it all' option) would leave the "
        "pipeline with NO shadow term whatsoever and would make low-sun "
        "scenes brighter than physics allows; the claim that the forward "
        "Jacobian map already contains the shadow test is false for this "
        "code (see det.py citations above)."),
    "when_it_bites": (
        "G1 = 1 exactly whenever alpha*tan(theta) >= 1/1.6, i.e. high sun or "
        "gentle water.  It only acts in the low-sun / steep-wave corner, "
        "which is also the corner where the transport's missing occlusion is "
        "worst -- the two are complementary, not additive."),
    "minimal_fix_applied": "none (no double counting found -> no v8 change)",
    "measured_change_on_reference_still": (
        "G1 at the reference scene is exactly 1.0 (a >> 1.6) for both the "
        "v7/v6 state and the new v8 state, so the correction applied to the "
        "reference still is multiplication by 1.0: the still is bit-identical "
        "with and without the G1 node (verified, see the gate report)."),
}


def g1_reference_table(rms_slopes, sun_altitudes_deg):
    """G1 over the (slope, sun altitude) corner for the audit report."""
    from v5_physics import smith_g1_beckmann
    rows = []
    for a in rms_slopes:
        for alt in sun_altitudes_deg:
            z = 90.0 - float(alt)
            g = smith_g1_beckmann(a, z)
            rows.append({"rms_slope": float(a), "sun_alt_deg": float(alt),
                         "zenith_deg": float(z), "a_param": float(g["a"]),
                         "Lambda": float(g["Lambda"]), "G1": float(g["G1"])})
    return rows


# --------------------------------------------------------------------------
# Q. solar-disc convergence
# --------------------------------------------------------------------------
def q_convergence(C_by_q, S_patch):
    """Convergence of C statistics and peak position with the disc quadrature."""
    qs = sorted(C_by_q.keys())
    per = {}
    for q in qs:
        c = np.asarray(C_by_q[q], dtype=np.float64)
        per[str(q)] = {"stats": quantiles(c),
                       "peak": peak_metrics(c),
                       "det_time_s": None}
    pairs = []
    for i in range(len(qs) - 1):
        a, b = qs[i], qs[i + 1]
        ca = np.asarray(C_by_q[a], dtype=np.float64)
        cb = np.asarray(C_by_q[b], dtype=np.float64)
        qa, qb = per[str(a)]["stats"], per[str(b)]["stats"]
        rel = {}
        for k in ("p50", "p99", "p99.9", "max"):
            den = abs(qb[k]) if qb[k] else 0.0
            rel[k] = abs(qa[k] - qb[k]) / den if den else None
        xa, xb = np.log1p(ca), np.log1p(cb)
        pairs.append({
            "Qa": int(a), "Qb": int(b),
            "prev": str(a), "next": str(b),
            "quantile_rel": rel,
            "mean_rel": float(abs(qa["mean"] - qb["mean"])),
            "peak": peak_shift(ca, cb, S_patch, "Q%d" % a, "Q%d" % b),
            "peak_luminance_abs_err": float(abs(per[str(a)]["peak"]["max"]
                                                - per[str(b)]["peak"]["max"])),
            "msssim_log1pC": ms_ssim(xa, xb)["msssim"],
            "rms_diff": float(np.sqrt(np.mean((ca - cb) ** 2))),
        })
    return {"Q_levels": [int(q) for q in qs], "per_Q": per, "pairs": pairs}


# --------------------------------------------------------------------------
# M. photon check + decomposition
# --------------------------------------------------------------------------
def photon_vs_det(C_det, C_ph, S_patch, still_det=None, still_ph=None,
                  gamma=1.0):
    """MS-SSIM of det vs photon, decomposed into luminance and position."""
    A = np.asarray(C_det, dtype=np.float64)
    B = np.asarray(C_ph, dtype=np.float64)
    Xa, Xb = np.log1p(gamma * A), np.log1p(gamma * B)
    ms_C = ms_ssim(Xa, Xb)
    qa, qb = quantiles(A), quantiles(B)
    pa, pb = peak_metrics(A), peak_metrics(B)
    shift = peak_shift(A, B, S_patch, "det", "photon")
    peak_abs = float(abs(pa["max"] - pb["max"]))
    peak_abs_999 = float(abs(pa["p99_9"] - pb["p99_9"]))
    out = {
        "gamma": float(gamma),
        "msssim_log1pC": ms_C,
        "quantiles_det": qa, "quantiles_photon": qb,
        "quantile_rel": {k: (abs(qa[k] - qb[k]) / abs(qb[k]) if qb[k] else None)
                         for k in ("p50", "p99", "p99.9")},
        "fwhm_det": ac_fwhm_px(A), "fwhm_photon": ac_fwhm_px(B),
        "mean_abs_diff_log1pC": float(np.mean(np.abs(Xa - Xb))),
        "decomposition": {
            "peak_luminance": {
                "det_max": pa["max"], "photon_max": pb["max"],
                "abs_error": peak_abs,
                "rel_error": peak_abs / pb["max"] if pb["max"] else None,
                "det_p99_9": pa["p99_9"], "photon_p99_9": pb["p99_9"],
                "abs_error_p99_9": peak_abs_999,
                "rel_error_p99_9": (peak_abs_999 / pb["p99_9"]
                                    if pb["p99_9"] else None),
                "det_top0.1pct_mean": pa["top_frac_mean"],
                "photon_top0.1pct_mean": pb["top_frac_mean"],
                "note": ("single-pixel max is noise-dominated in the photon "
                         "field (Poisson); p99.9 and the top-0.1% mean are "
                         "the robust peak-luminance readouts"),
            },
            "peak_position": shift,
        },
    }
    if still_det is not None and still_ph is not None:
        sd = np.asarray(still_det, dtype=np.float64) / 255.0
        sp = np.asarray(still_ph, dtype=np.float64) / 255.0
        gd = sd.mean(axis=2)
        gp = sp.mean(axis=2)
        ms_s = ms_ssim(gd, gp)
        # the still's own brightest pixel is the PHOTO's brightest patch, not
        # a caustic, so the peak is read on a high-passed (caustic) version
        hd = gd - ndimage.gaussian_filter(gd, 12.0, mode="reflect")
        hp = gp - ndimage.gaussian_filter(gp, 12.0, mode="reflect")
        yd = np.unravel_index(int(np.argmax(hd)), hd.shape)
        yp = np.unravel_index(int(np.argmax(hp)), hp.shape)
        H, W = gd.shape
        mm_x = float(S_patch) / W
        mm_y = float(S_patch) / H
        out["still"] = {
            "msssim_grey": ms_s,
            "peak_definition": ("maximum of the grey high-pass (grey minus a "
                                "12 px Gaussian), i.e. the strongest caustic "
                                "filament, not the photo's brightest patch"),
            "peak_px_det": [int(yd[0]), int(yd[1])],
            "peak_px_photon": [int(yp[0]), int(yp[1])],
            "peak_shift_px": float(math.hypot(yd[0] - yp[0], yd[1] - yp[1])),
            "peak_shift_mm_xy": [float(abs(yd[1] - yp[1]) * mm_x),
                                 float(abs(yd[0] - yp[0]) * mm_y)],
            "peak_hp_value_det": float(hd.max()),
            "peak_hp_value_photon": float(hp.max()),
            "peak_hp_abs_err": float(abs(hd.max() - hp.max())),
            "peak_hp_rel_err": (float(abs(hd.max() - hp.max()) / hp.max())
                                if hp.max() else None),
            "mean_abs_diff_uint8_equiv": float(np.mean(np.abs(sd - sp)) * 255),
            "note": ("the still level is the deliverable level: the same "
                     "composite is rendered from the det master and from the "
                     "photon master, so the difference is transport only"),
        }
    return out
