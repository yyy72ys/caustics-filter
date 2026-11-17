"""transport2 -- BRIEF-V4-1-7: turn the printed lattice into fine noise.

NEW FILE.  Nothing else in the tree is touched: `det.py`, `transport.py`,
`stream_surface3.py`, `stream_surface4.py`, `run_v16_water.py` and every
driver stay byte-identical (G-LEGACY).  The accounting comes from
`stream_surface3`, the surface field is v4-1.4's *own* field (NO band limit,
unlike v4-1.5): the slope field is passed to the transport untouched.

WHY v4-1.5 WAS NOT ENOUGH (BRIEF-V4-1-7)
----------------------------------------
v4-1.5 killed the lattice by spectrally band-limiting the slope field to
sigma = 4 cells before the splat.  The lattice went away, but so did the fine
slopes -- i.e. the light that is actually focused (C std 0.528 -> 0.101).  The
diagnosis in the brief is that the *deposit* stage, not the field, is what
prints the lattice: a forward splat fires one ray per source texel from an
evenly spaced lattice and counts landings as points, so the residual deposit
error is a periodic function of the *destination* lattice, i.e. a lattice /
moire.  v4-1.4's |ds/di| table (REPORT-V4-1-5 section 1.1) shows |ds/di| has
median 1.09 at this scene, so the landing map folds almost everywhere and the
deposit is a print of the ray lattice.

THE FIX (deposit stage only -- the field is v4-1.4's, verbatim)
--------------------------------------------------------------
(1) RAY-POSITION JITTER.  The launch position of each source texel is shifted
    inside its own cell, per frame, by a deterministic 2-D low-discrepancy
    sequence.  The regularity of the printed lattice is destroyed; the residual
    deposit error becomes broadband noise instead of a periodic print.

    SPEC (fixed here, not tuned):
        n  = i + N_s * j                      (flat source-cell index)
        f  = frame index, 0-based integer
        g  = 1.324717957244746                (plastic number)
        a1 = 1/g   = 0.7548776662466927       (R2 sequence, Roberts 2018)
        a2 = 1/g^2 = 0.5698402909980532
        jx = frac(0.5 + a1 * (n + 13 * f)) - 0.5     in [-0.5, 0.5) cells
        jy = frac(0.5 + a2 * (n +  7 * f)) - 0.5
        launch = ((i + 0.5 + jx) * h, (j + 0.5 + jy) * h),  h = S / N_s
    The strides 13 and 7 are coprime with each other and with N_s = 512, so
    consecutive frames visit disjoint phases of the sequence; both are fixed
    a priori (they are the standard "rank-1 lattice" frame strides 1 and 2
    scaled by the co-prime offsets used in the R2 family).  The five field
    arrays (h, gx, gy, Dx, Dy) are sampled bilinearly (periodic) at the
    jittered point -- the ray carries the surface state of the point it is
    actually launched from, so the jitter does not bias the displacement.

(2) KERNEL DEPOSITION.  Landings are counted with a normalised Gaussian
    kernel, not as points.  The width is derived from the ray density and the
    LOCAL JACOBIAN of the landing map, never from a gate:

        nominal spacing in the bottom plane, one quadrature direction:
            s0 = S / N_s                                  (= 1 cell)
        landing map:  x_b = x_s + T(x_s) * t(x_s)
            J = d x_b / d x_s          (periodic central differences)
        local areal compression: |det J|
        local spacing:            s(x) = s0 * sqrt(|det J|)
        kernel width:             sigma(x) = s(x) / 2

    i.e. one Gaussian sigma per local inter-ray spacing (the matched filter
    for a sample sequence of spacing s has sigma = s/2: one sample per sigma).
    The single direction is used on purpose: the Q solar-disc directions are
    near-coincident copies (disc radius << one cell), so they do NOT supply
    independent samples at the scale of the deposit lattice -- they only reduce
    the variance, which the kernel already averages away.  Dividing s0 by
    sqrt(Q) instead would give a kernel 5x narrower than one buffer sample and
    would make the whole kernel stage a no-op; that alternative is recorded in
    the report as the rejected reading.

    Numerical guard, fixed a priori: |det J| is clipped to [1/4, 4], so
    sigma in [0.25, 1.0] cell.  Outside that band the deposit carries no
    information the frame can show (4x expansion = dark, detail-free; 4x
    compression = denser than the buffer's own sample grid, where a narrower
    kernel only adds shot noise).  The kernel is truncated at 3 sigma and
    renormalised to unit mass, so it is mean preserving.

    Implementation: sigma is quantised into K = 6 geometric bins over
    [0.25, 1.0] cell (ratio 4^(1/5) = 1.3195, each ray assigned to the nearest
    bin centre), the rays are deposited (bilinear, periodic) into one buffer
    per bin, and each buffer is convolved with the separable periodic Gaussian
    of its bin.  That is *exactly* splatting every ray with its own bin
    kernel; it is only cheaper.

(3) BUFFER = 2x FRAME RESOLUTION, THEN AREA AVERAGE.  The deposit buffer is
    1472 x 1080 = 2 x (736 x 540).  The frame-resolution C is the exact 2 x 2
    area average of that buffer, and is then corrected to mean exactly 1
    (energy conservation, as the brief requires; v4-1.4 forbade any rescaling,
    v4-1.7 requires it).

WHAT IS UNCHANGED
-----------------
* the field: v4-1.4's h/gx/gy/Dx/Dy untouched (no band limit, no filter).
* the accounting: sigma_reg / var_i / var(ii) / var(iii) from
  stream_surface3; this module never reads or writes a budget number.
* the physics: Snell refraction with the same representative IOR, the same
  Q-point Fibonacci solar disc, the same finite depth, the same analytic
  divisor (total_rays / n_buffer_pixels) before the mean correction.
* the image-formation warp: it consumes the field itself, as in v4-1.5, so
  G-WARP is still a "we did not weaken the wave" yardstick.
"""

import numpy as np

PATCH = "v4-1.7"

# ---- (1) jitter spec -------------------------------------------------------
PLASTIC = 1.324717957244746
JITTER_A1 = 1.0 / PLASTIC              # 0.7548776662466927
JITTER_A2 = 1.0 / (PLASTIC * PLASTIC)  # 0.5698402909980532
JITTER_STRIDE_X = 13
JITTER_STRIDE_Y = 7

# ---- (2) kernel spec -------------------------------------------------------
KERNEL_SIGMA_PER_SPACING = 0.5     # sigma = s/2
JAC_CLIP = (0.25, 4.0)             # |det J| clip -> sigma in [0.25, 1.0] cell
KERNEL_BINS = 6                    # geometric bins over [0.25, 1.0] cell
KERNEL_TRUNCATE_SIGMA = 3.0

# ---- (3) buffer spec -------------------------------------------------------
BUFFER_SCALE = 2                   # 2 x frame resolution

IOR_G = 1.33434
SIGMA_PX_LEGACY = 0.75             # det.SIGMA_PX (kept as a record only)

RECORD = {
    "patch": PATCH,
    "module": "transport2",
    "field": "v4-1.4 h/gx/gy/Dx/Dy verbatim (NO band limit)",
    "jitter": {
        "kind": "rank-1 (R2) additive recurrence, deterministic",
        "sequence": "R2 (Roberts 2018), g = plastic number",
        "a1": JITTER_A1, "a2": JITTER_A2,
        "stride_x_frames": JITTER_STRIDE_X, "stride_y_frames": JITTER_STRIDE_Y,
        "formula_x": "frac(0.5 + a1*(n + %d*f)) - 0.5" % JITTER_STRIDE_X,
        "formula_y": "frac(0.5 + a2*(n + %d*f)) - 0.5" % JITTER_STRIDE_Y,
        "n": "i + N_s*j (flat source-cell index)",
        "amplitude_cells": 1.0,
        "field_sampling": "periodic bilinear at the jittered launch point",
    },
    "kernel": {
        "shape": "normalised Gaussian, truncated at %.1f sigma" % KERNEL_TRUNCATE_SIGMA,
        "sigma_rule": "sigma = (S/N_s) * sqrt(|det J|) / 2, clipped",
        "sigma_per_spacing": KERNEL_SIGMA_PER_SPACING,
        "jacobian_clip": list(JAC_CLIP),
        "sigma_cell_range": [KERNEL_SIGMA_PER_SPACING * float(JAC_CLIP[0] ** 0.5),
                             KERNEL_SIGMA_PER_SPACING * float(JAC_CLIP[1] ** 0.5)],
        "bins": KERNEL_BINS,
        "bin_ratio": float((1.0 / 0.25) ** (1.0 / (KERNEL_BINS - 1))),
        "applied_as": ("separable periodic Gaussian on the binned deposit "
                       "(equivalent to per-ray splatting)"),
    },
    "buffer": {"scale": BUFFER_SCALE,
               "shape": [1080, 1472],
               "downsample": "2x2 area average -> frame resolution"},
    "normalisation": "C divided by its own mean (mean == 1 exactly)",
    "unchanged": ["v4-1.4 slope field", "stream_surface3 accounting",
                  "Snell physics / IOR / Q / Fibonacci disc",
                  "analytic divisor before the mean correction",
                  "image-formation warp field (unfiltered)"],
}


# ---------------------------------------------------------------------------
# (1) jitter
# ---------------------------------------------------------------------------
def jitter_offsets(N_s, frame):
    """Per-cell launch offsets (in cells) for one frame -- the spec above.

    Returns (jx, jy), each (N_s, N_s) float64 in [-0.5, 0.5).
    """
    i = np.arange(int(N_s), dtype=np.float64)
    n = i[None, :] + float(N_s) * i[:, None]
    f = float(frame)
    jx = np.mod(0.5 + JITTER_A1 * (n + JITTER_STRIDE_X * f), 1.0) - 0.5
    jy = np.mod(0.5 + JITTER_A2 * (n + JITTER_STRIDE_Y * f), 1.0) - 0.5
    return jx, jy


def _bilinear_periodic(F, ux, uy, N):
    """Sample F (N,N) at (ux,uy) given in cell units, periodic, bilinear."""
    x0 = np.floor(ux)
    y0 = np.floor(uy)
    fx = ux - x0
    fy = uy - y0
    x0 = x0.astype(np.int64) % N
    y0 = y0.astype(np.int64) % N
    x1 = (x0 + 1) % N
    y1 = (y0 + 1) % N
    return (F[y0, x0] * (1.0 - fx) * (1.0 - fy)
            + F[y0, x1] * fx * (1.0 - fy)
            + F[y1, x0] * (1.0 - fx) * fy
            + F[y1, x1] * fx * fy)


# ---------------------------------------------------------------------------
# optics helpers (re-used verbatim from transport.py / det.py)
# ---------------------------------------------------------------------------
def _basis_and_disc(env, Q, ior):
    from transport import orthonormal_basis, sun_disc_radius
    from det import fibonacci_disc
    s, u, v = orthonormal_basis(env["sunDir"])
    disc = sun_disc_radius(env.get("clarity", 0.8))
    quads = fibonacci_disc(int(Q), disc)
    return s, u, v, disc, quads, 1.0 / float(ior)


def _wrap_diff(a, S):
    """Periodic central difference along both axes (cell units)."""
    d = 0.5 * (np.roll(a, -1, axis=1) - np.roll(a, 1, axis=1))
    d = np.mod(d + 0.5 * S, S) - 0.5 * S
    return d


def _wrap_diff_y(a, S):
    d = 0.5 * (np.roll(a, -1, axis=0) - np.roll(a, 1, axis=0))
    d = np.mod(d + 0.5 * S, S) - 0.5 * S
    return d


def sigma_bins():
    """Bin centres (cells) and the log-spacing used to assign rays.

    sigma = KERNEL_SIGMA_PER_SPACING * sqrt(|det J|), so the clip
    |det J| in [1/4, 4] maps to sigma in [0.25, 1.0] cell.
    """
    lo = KERNEL_SIGMA_PER_SPACING * float(np.sqrt(JAC_CLIP[0]))   # 0.25 cell
    hi = KERNEL_SIGMA_PER_SPACING * float(np.sqrt(JAC_CLIP[1]))   # 1.00 cell
    K = int(KERNEL_BINS)
    ratio = (hi / lo) ** (1.0 / (K - 1))
    centres = lo * ratio ** np.arange(K, dtype=np.float64)
    return centres, ratio, lo, hi


# ---------------------------------------------------------------------------
# the transport
# ---------------------------------------------------------------------------
def det_frame2(fields, S, env, Q, frame_index, n_buf=(1080, 1472),
               ior=IOR_G, return_info=True):
    """Jittered, kernel-deposited forward transport for one time slice.

    fields : dict with h/gx/gy/Dx/Dy (N_s,N_s).  Used EXACTLY as given.
    S      : patch size (m).
    n_buf  : (rows, cols) of the deposit buffer (default 1080x1472).
    returns (C_frame, info) with C_frame (540,736) float32, mean exactly 1.
    """
    h = np.asarray(fields["h"], dtype=np.float64)
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    Dx = np.asarray(fields["Dx"], dtype=np.float64)
    Dy = np.asarray(fields["Dy"], dtype=np.float64)
    N_s = h.shape[0]
    S = float(S)
    Hb, Wb = int(n_buf[0]), int(n_buf[1])
    depth = float(env["depthM"])
    Q = int(Q)

    s, u, v, disc, quads, eta = _basis_and_disc(env, Q, ior)

    # ---- (1) jittered launch points + field sampled there ------------------
    jx, jy = jitter_offsets(N_s, int(frame_index))
    h_cell = S / float(N_s)
    ux = np.arange(N_s, dtype=np.float64)[None, :] + 0.5 + jx
    uy = np.arange(N_s, dtype=np.float64)[:, None] + 0.5 + jy
    hh = _bilinear_periodic(h, ux, uy, N_s)
    GGX = _bilinear_periodic(gx, ux, uy, N_s)
    GGY = _bilinear_periodic(gy, ux, uy, N_s)
    DDX = _bilinear_periodic(Dx, ux, uy, N_s)
    DDY = _bilinear_periodic(Dy, ux, uy, N_s)
    X0 = ux * h_cell
    Y0 = uy * h_cell

    nl = np.sqrt(GGX * GGX + GGY * GGY + 1.0)
    Nx = -GGX / nl
    Ny = -GGY / nl
    Nz = 1.0 / nl
    del GGX, GGY, nl

    # ---- (2) deposit geometry ---------------------------------------------
    centres, ratio, sig_lo, sig_hi = sigma_bins()
    K = centres.size
    nb = Hb * Wb
    acc = np.zeros(K * nb, dtype=np.float64)
    cell_x = S / float(Wb)
    cell_y = S / float(Hb)

    fold = {"ds_per_ray_p50": [], "ds_per_ray_p99": [], "ds_per_ray_max": [],
            "folded_area_fraction": []}
    detj_rows = []
    bin_counts = np.zeros(K, dtype=np.float64)

    for q in range(Q):
        jxq, jyq = float(quads[q, 0]), float(quads[q, 1])
        d = s + u * jxq + v * jyq
        d = d / np.linalg.norm(d)
        I = -d
        Ix, Iy, Iz = float(I[0]), float(I[1]), float(I[2])
        cos1 = -(Ix * Nx + Iy * Ny + Iz * Nz)
        cos1 = np.clip(cos1, -1.0, 1.0)
        kk = 1.0 - eta * eta * (1.0 - cos1 * cos1)
        kk = np.maximum(kk, 0.0)
        sq = np.sqrt(kk)
        Tx = eta * Ix + (eta * cos1 - sq) * Nx
        Ty = eta * Iy + (eta * cos1 - sq) * Ny
        Tz = eta * Iz + (eta * cos1 - sq) * Nz
        Tz = np.minimum(Tz, -1e-6)
        t = (hh + depth) / (-Tz)
        xb = np.mod(X0 + DDX + Tx * t, S)
        yb = np.mod(Y0 + DDY + Ty * t, S)

        # local Jacobian of the landing map (unit-mass panels, cell units)
        dxb_di = _wrap_diff(xb, S) / h_cell
        dyb_di = _wrap_diff(yb, S) / h_cell
        dxb_dj = _wrap_diff_y(xb, S) / h_cell
        dyb_dj = _wrap_diff_y(yb, S) / h_cell
        detJ = dxb_di * dyb_dj - dxb_dj * dyb_di
        aJ = np.abs(detJ)
        detj_rows.append((float(np.percentile(aJ, 50)),
                          float(np.percentile(aJ, 99)), float(aJ.max())))
        # fold diagnostic (source-space gradient, same definition as v4-1.5;
        # periodic differences so the wrap seam does not fake a 512-cell jump)
        dsp = np.concatenate([(np.abs(_wrap_diff(xb, S)) / h_cell).ravel(),
                              (np.abs(_wrap_diff_y(yb, S)) / h_cell).ravel()])
        fold["ds_per_ray_p50"].append(float(np.percentile(dsp, 50)))
        fold["ds_per_ray_p99"].append(float(np.percentile(dsp, 99)))
        fold["ds_per_ray_max"].append(float(dsp.max()))
        fold["folded_area_fraction"].append(float(np.mean(dsp >= 1.0)))
        del dsp

        # sigma -> nearest bin centre
        sig_cell = KERNEL_SIGMA_PER_SPACING * np.sqrt(
            np.clip(aJ, JAC_CLIP[0], JAC_CLIP[1]))
        kbin = np.rint(np.log(sig_cell / sig_lo) / np.log(ratio))
        kbin = np.clip(kbin, 0, K - 1).astype(np.int64)
        bin_counts += np.bincount(kbin.ravel(), minlength=K)
        del sig_cell, aJ, detJ, dxb_di, dyb_di, dxb_dj, dyb_dj

        # bilinear deposit at buffer resolution
        fx = xb / cell_x
        fy = yb / cell_y
        ix0 = np.floor(fx)
        iy0 = np.floor(fy)
        tx = fx - ix0
        ty = fy - iy0
        ix0 = ix0.astype(np.int64) % Wb
        iy0 = iy0.astype(np.int64) % Hb
        ix1 = (ix0 + 1) % Wb
        iy1 = (iy0 + 1) % Hb
        base = kbin * nb
        i00 = base + iy0 * Wb + ix0
        i10 = base + iy0 * Wb + ix1
        i01 = base + iy1 * Wb + ix0
        i11 = base + iy1 * Wb + ix1
        idx = np.concatenate([i00.ravel(), i10.ravel(), i01.ravel(), i11.ravel()])
        wts = np.concatenate([((1 - tx) * (1 - ty)).ravel(),
                              (tx * (1 - ty)).ravel(),
                              ((1 - tx) * ty).ravel(),
                              (tx * ty).ravel()])
        acc += np.bincount(idx, weights=wts, minlength=K * nb)
        del idx, wts, i00, i10, i01, i11, base, fx, fy, ix0, iy0, ix1, iy1
        del tx, ty, xb, yb, Tx, Ty, Tz, t, cos1, kk, sq

    # ---- (2b) the kernels -------------------------------------------------
    from scipy import ndimage as ndi
    buf = np.zeros((Hb, Wb), dtype=np.float64)
    sig_samples = []
    for k in range(K):
        blk = acc[k * nb:(k + 1) * nb].reshape(Hb, Wb)
        sg = centres[k] * np.array([Hb / float(N_s), Wb / float(N_s)])
        buf += ndi.gaussian_filter(blk, sigma=sg, mode="wrap",
                                   truncate=KERNEL_TRUNCATE_SIGMA)
        sig_samples.append([float(sg[0]), float(sg[1])])
        del blk
    del acc

    total_rays = float(Q) * float(N_s) * float(N_s)
    divisor = total_rays / float(nb)
    C_buf = buf / divisor
    mean_buf = float(np.mean(C_buf))

    # ---- (3) 2x2 area average -> frame resolution -------------------------
    C_frame = 0.25 * (C_buf[0::2, 0::2] + C_buf[0::2, 1::2]
                      + C_buf[1::2, 0::2] + C_buf[1::2, 1::2])
    mean_frame = float(np.mean(C_frame))
    C_frame = (C_frame / mean_frame) if mean_frame > 0 else C_frame
    C_frame = C_frame.astype(np.float32)

    info = {
        "Q": Q, "ior": float(ior), "disc_rad": float(disc),
        "N_s": int(N_s), "S_m": S, "h_cell_m": float(h_cell),
        "buffer_hw": [Hb, Wb], "buffer_scale": BUFFER_SCALE,
        "total_rays": total_rays, "divisor_analytical": divisor,
        "mean_before_mean_correction": mean_frame,
        "mean_buffer_before_correction": mean_buf,
        "mean_correction_factor": (1.0 / mean_frame if mean_frame else 1.0),
        "jitter_amplitude_cells": 1.0,
        "jitter_frame": int(frame_index),
        "kernel_bins": K,
        "kernel_sigma_centres_cells": [float(x) for x in centres],
        "kernel_sigma_centres_buffer_samples_yx": sig_samples,
        "kernel_ray_fraction_per_bin": [float(x / max(bin_counts.sum(), 1.0))
                                       for x in bin_counts],
        "kernel_truncate_sigma": KERNEL_TRUNCATE_SIGMA,
        "jacobian_clip": list(JAC_CLIP),
        "detJ_abs_p50_p99_max": detj_rows,
        "fold_per_quadrature": {
            k: {"min": float(np.min(v)), "max": float(np.max(v)),
                "mean": float(np.mean(v))} for k, v in fold.items()},
        "clamp": "kernel sigma clipped by the |det J| guard; C not clamped",
    }
    return C_frame, info
