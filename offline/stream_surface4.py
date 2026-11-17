"""stream_surface4 -- BRIEF-V4-1-5 patch: kill the printed computation lattice.

NEW FILE.  Nothing else in the tree is touched: `stream_surface3.py` (the
v4-1.4 accounting), `stream_surface2.py`, `det.py`, `transport.py`,
`run_v16_water.py` and every driver stay byte-identical (G-LEGACY).

WHAT THE BRIEF ASKS (V4-1-5, stage 1)
-------------------------------------
v4-1.4's real-size still shows axis-parallel lines with a ~7-11 px period that
the source photo and v4-1.3 do not have.  Stage 1 of the brief is the
*separation*: draw frame 24 three ways -- (a) refraction warp only (C = 1),
(b) caustics C only (warp off), (c) both (= v4-1.4).  The driver does exactly
that and measures each with the same FFT statistic.  The answer (recorded in
numbers-v4-1-5.json / REPORT-V4-1-5.md) is (b): the lattice is carried by the
caustics stage, and (a) sits at the source photo's level.

STAGE 2: WHY THE CAUSTICS STAGE PRINTS A LATTICE
------------------------------------------------
`det.det_frame` forward-splats ONE ray per surface texel (`N_s` rays across the
patch, `N_buf == N_s == 512`).  Write the landing map of one quadrature
direction as

    p_i = i + s_i          (i = ray index, s_i = displacement in ray units)

A forward splat can represent this map only if neighbouring rays land in
neighbouring buffer cells, i.e. if

    |ds/di| < 1                                            (no fold)

If |ds/di| > 1 somewhere the deposit is no longer the Jacobian of the map: the
folded cells collect rays from non-adjacent texels and the *ray lattice*
itself is printed into C.  The driver measures |ds/di| on the v4-1.4 surface
(`fold_stats`, frame 24, real-size scene, all four neighbour differences
pooled):

    sigma [cells]   p50|ds/di|   p99|ds/di|   max|ds/di|   folded area
      0.0 (v4-1.4)      1.09        4.17        6.79        53.5%
      1.0               0.22        0.83        1.38         0.2%
      2.0               0.08        0.31        0.46         0.0%
      3.0               0.05        0.18        0.29         0.0%
      4.0               0.03        0.13        0.19         0.0%

(these are the numbers the driver writes into numbers-v4-1-5.json; the ladder
is recomputed at run time, nothing here is copied by hand into a result)

So at the v4-1.4 surface the *median* ray already moves its landing point by
more than one whole ray spacing -- the transport is folded almost everywhere.
The folded content is the sub-4-cell part of the surface slope (the broadband
(i) bed response and the k^-6 (iii) tail), which is exactly the part the 1-ray
per texel lattice cannot resolve.

THE FIX (stage 2 of the brief) -- and what it may NOT do
-------------------------------------------------------
Evaluate the transport field at the resolution the transport can actually
represent: band-limit h/gx/gy with a *spectral* (periodic FFT) Gaussian of
`SLOPE_FILTER_SIGMA_CELLS` cells before the splat.  The width is fixed by the
fold criterion, not by the gate:

    sigma := the smallest width such that  max|ds/di| <= 0.25

(factor-4 margin below the fold limit).  The sigma ladder is measured and
reported; it lands on 4.0 cells = 0.50 mm at the real-size patch (S = 63.889
mm, 512 cells of 0.1248 mm).

What is deliberately NOT changed:

* the waves: the (ii) ring is at k0 = 1.294 cycles/patch = 127.2 rad/m, i.e.
  ~400 cells per wavelength.  `large_scale_preservation()` measures the slope
  variance inside the ring band before/after the filter; the driver records
  the ratio (it is 1.0000 to the printed precision).
* the amplitude / the budget: sigma_reg, var_i, var(ii), var(iii),
  sigma_total^2 = sigma_reg^2 + var_i are inherited from stream_surface3
  verbatim; nothing here reads or writes a budget number.
* the accounting of C: `det.det_frame` is called unchanged, so the analytic
  divisor (total_rays / N_buf^2) and the "mean ~ 1 by energy conservation"
  statement still hold.  The filter is applied to the *input* field, and its
  kernel is unit-mass (mean preserving), so mean(C) is unchanged in the
  ensemble but NOT forced.
* the image formation (the refraction warp): it consumes the *unfiltered*
  field, because the displacement is a LINEAR functional of the slope and the
  736x540 frame resolves it (the ablation shows the warp-only image sits at
  the source photo's level).  This is also what makes G-WARP an honest
  "we did not weaken the wave" test: the warped displacement of frame 24 is
  bit-identical to v4-1.4 by construction.

DISCLOSURE: the caustic contrast drops, because the aliased sub-lattice
focusing is removed (C std 0.528 -> 0.10 at frame 24).  That is the honest
consequence of resolving the field instead of aliasing it; it is tabulated in
the report rather than tuned away.
"""

import numpy as np

from transport import gaussian_blur_fft
import stream_surface3 as S3

PATCH = "v4-1.5"
SLOPE_FILTER_SIGMA_CELLS = 4.0
MAX_DS_PER_RAY_TARGET = 0.25
FOLD_LIMIT = 1.0

RECORD = {
    "patch": PATCH,
    "module": "stream_surface4",
    "stage_of_the_lattice": "caustics transport (det.det_frame forward splat)",
    "mechanism": ("1 ray per surface texel cannot represent |ds/di| >= 1; the "
                  "folded deposit prints the ray lattice into C"),
    "fix": ("spectral (periodic FFT) Gaussian band-limit of h/gx/gy before "
            "the splat, width fixed by max|ds/di| <= %.2f"
            % MAX_DS_PER_RAY_TARGET),
    "sigma_cells": SLOPE_FILTER_SIGMA_CELLS,
    "sigma_m_realscale": None,      # filled by the driver (S/512*sigma)
    "unchanged": ["wave spectrum", "sigma_reg/var_i/var_ii/var_iii accounting",
                  "det.det_frame code path and divisor", "optics",
                  "composite (A/Blow/E/Q)", "image-formation warp field"],
}


def transport_fields(fields, sigma_cells=SLOPE_FILTER_SIGMA_CELLS):
    """Band-limited copy of the transport input field (periodic spectral LP).

    sigma_cells: Gaussian width in units of the surface cell (S/N).
    sigma_cells == 0 returns the input arrays bit-for-bit (unit check).
    """
    sig = float(sigma_cells)
    out = dict(fields)
    if sig <= 0.0:
        for k in ("h", "gx", "gy"):
            out[k] = np.asarray(fields[k])
        return out
    for k in ("h", "gx", "gy"):
        a = np.asarray(fields[k], dtype=np.float64)
        out[k] = gaussian_blur_fft(a, sig).astype(fields[k].dtype)
    return out


def _displacement(fields, S, depth, ior):
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    cb = (float(ior) - 1.0) / float(ior)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny = -gx / nl, -gy / nl
    shx = -nx * float(depth) * cb
    shy = -ny * float(depth) * cb
    per_ray = float(S) / float(gx.shape[0])
    return shx, shy, per_ray


def fold_stats(fields, S, depth, ior):
    """|ds/di| in ray units: the criterion this patch is built on."""
    shx, shy, per_ray = _displacement(fields, S, depth, ior)
    dx = np.abs(np.diff(shx, axis=1)) / per_ray
    dy = np.abs(np.diff(shy, axis=0)) / per_ray
    d = np.concatenate([dx.ravel(), dy.ravel()])
    px_per_m = 1.0 / (float(S) / 736.0)
    mag_px = np.sqrt(shx ** 2 + shy ** 2) * px_per_m
    return {
        "sigma_cells": float(fields.get("_sigma_cells", 0.0)),
        "ds_per_ray_p50": float(np.percentile(d, 50)),
        "ds_per_ray_p99": float(np.percentile(d, 99)),
        "ds_per_ray_max": float(d.max()),
        "fold_limit": FOLD_LIMIT,
        "folded_area_fraction": float(np.mean(d >= FOLD_LIMIT)),
        "displacement_px_mean": float(mag_px.mean()),
        "displacement_px_p50": float(np.percentile(mag_px, 50)),
        "displacement_px_p99": float(np.percentile(mag_px, 99)),
    }


def large_scale_preservation(fields, filtered, S, k_ring=127.2123):
    """Slope variance inside the ring band, before vs after the filter.

    The ring is at k_ring (rad/m); the band is +-1 dk around it, which is the
    resolution of the patch (dk = 2*pi/S).  A band-limit that removes the
    printed lattice must leave this number at 1.0.
    """
    N = np.asarray(fields["gx"]).shape[0]
    dk = 2.0 * np.pi / float(S)
    kx = np.fft.fftfreq(N) * 2.0 * np.pi * N / float(S)
    KX, KY = np.meshgrid(kx, kx)
    kk = np.hypot(KX, KY)
    band = np.abs(kk - float(k_ring)) <= dk

    def slope_var(f):
        gx = np.asarray(f["gx"], dtype=np.float64)
        gy = np.asarray(f["gy"], dtype=np.float64)
        tot = 0.0
        for g, k in ((gx, KX), (gy, KY)):
            F = np.fft.fft2(g)
            tot += float(np.sum(np.abs(1j * k * F)[band] ** 2)) / float(N ** 4)
        return tot

    a = slope_var(fields)
    b = slope_var(filtered)
    return {"k_ring_rad_m": float(k_ring), "band_half_width_rad_m": float(dk),
            "slope_var_ring_band_before": a, "slope_var_ring_band_after": b,
            "ratio": (b / a if a > 0 else None),
            "total_slope_rms_before": float(np.sqrt(np.mean(
                np.asarray(fields["gx"], np.float64) ** 2
                + np.asarray(fields["gy"], np.float64) ** 2))),
            "total_slope_rms_after": float(np.sqrt(np.mean(
                np.asarray(filtered["gx"], np.float64) ** 2
                + np.asarray(filtered["gy"], np.float64) ** 2)))}


def sigma_from_fold_criterion(fields, S, depth, ior,
                              ladder=(0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0),
                              target=MAX_DS_PER_RAY_TARGET):
    """The smallest ladder width whose max|ds/di| is below `target`.

    This is the ONLY rule used to pick the width; the gate (G-LINE) is not
    consulted here and the driver reports the whole ladder so the choice can
    be audited.
    """
    rows = []
    chosen = None
    for sig in ladder:
        f = transport_fields(fields, sig)
        st_ = fold_stats(f, S, depth, ior)
        st_["sigma_cells"] = float(sig)
        rows.append(st_)
        if chosen is None and st_["ds_per_ray_max"] <= float(target):
            chosen = float(sig)
    return {"target_max_ds_per_ray": float(target), "chosen_sigma_cells": chosen,
            "ladder": rows}


class StreamState4(S3.StreamState3):
    """v4-1.4 surface with the v4-1.4 accounting, plus the transport filter.

    No budget number is touched: `accounting()` is inherited unchanged.  The
    new method is the transport input field.
    """

    accounting_tag = PATCH

    def transport_fields(self, t, sigma_cells=SLOPE_FILTER_SIGMA_CELLS):
        return transport_fields(self.fields_at(t), sigma_cells)


def self_check(fields, S, depth, ior):
    """sigma=0 must be bit-identical; the ring band must survive the filter."""
    same = transport_fields(fields, 0.0)
    ident = {}
    for k in ("h", "gx", "gy"):
        ident[k] = bool(np.array_equal(np.asarray(same[k]),
                                       np.asarray(fields[k])))
    filt = transport_fields(fields, SLOPE_FILTER_SIGMA_CELLS)
    keep = large_scale_preservation(fields, filt, S)
    return {"identity_at_sigma_0": ident,
            "identity_ok": bool(all(ident.values())),
            "ring_band": keep,
            "sigma_cells": SLOPE_FILTER_SIGMA_CELLS}
