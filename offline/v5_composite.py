"""v5 composite (PIPELINE-V5 chunk 1, report sections 1.3 / 2.1 / 2.2 / 2.3).

The OUTPUT-SPEC section 3 separation is kept unchanged:

    E_bottom(q,t) = E_amb(d) + E_base,dir(d) * C(q,t)

and the four new physics terms are inserted where the report says they
belong:

  1.3  Smith G1(sun) multiplies the DIRECT term (statistical approximation;
       ray-exact visibility is NOT implemented, and the two are never
       combined -- report section 1.3 / 7).
  2.1  Bottom albedo rho_b is a per-channel gain on the photo's linear
       values, applied to the direct-bottom term only (the ambient term
       arrived by a different path).
  2.2  Beer-Lambert exp(-c(lambda) s_down) on the way down and
       exp(-c(lambda) s_up) on the way back up replaces the old ad-hoc
       "absorb * path_len * 0.22" depth darkening.
  2.3  Turbidity scattering blur: C is convolved with a Gaussian whose
       sigma grows as sqrt(b*z).

The v3 composition constants (A = ambient scale, B = displacement lowpass,
C-layer exposure E) are kept as parameters so the baseline still is the v3
final composition under the same conditions.  The C-layer exposure is
renormalised once by the ratio of the flat-water two-way attenuation
(v3 ad-hoc / v5 physical) so the DISPLAY reference level is unchanged; that
is a display transform, not physics, and is recorded in numbers-v5.json.
"""
import numpy as np

from image import (IOR, SKY_REF, srgb_to_linear, linear_to_srgb,
                   aces_narkowicz, e_base_dir_from_env, render_frame_v3)
from state import sample_bicubic
from transport import gaussian_blur_fft

REC709 = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)

# Frozen v3-final composition constants (PIPELINE-V4 prod freeze).
A_V3 = 2.6
B_SIGMA_V3 = 3.0
E_V3 = 2.2


def to_lum(rgb):
    return (np.asarray(rgb, dtype=np.float64) * REC709).sum(axis=-1)


def image_metrics(lin, aces):
    lin = np.asarray(lin, dtype=np.float64)
    aces = np.asarray(aces, dtype=np.float64)
    lin_lum = to_lum(lin)
    aces_lum = to_lum(np.clip(aces, 0.0, 1.0))
    return {
        "linear_mean_pre_tonemap": float(np.mean(lin)),
        "linear_luminance_mean": float(np.mean(lin_lum)),
        "crush_rate_linear_lt_0.02": float(np.mean(lin_lum < 0.02)),
        "blowout_rate_aces_gt_0.98": float(np.mean(aces_lum > 0.98)),
        "linear_rgb_mean": [float(x) for x in lin.reshape(-1, 3).mean(0)],
        "aces_luminance_max": float(aces_lum.max()),
    }


def v3_beer(env):
    """The old ad-hoc depth darkening, verbatim from image.render_frame_v3."""
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    absorb = np.asarray(env["absorb"], dtype=np.float64)
    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    return np.exp(-absorb * path_len * 0.22)


def flat_water_direct_level(env, beer):
    """rec709 luminance of E_base,dir * beer (the flat-water direct level)."""
    return float((e_base_dir_from_env(env) * beer * REC709).sum())


def _smoothstep01(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def _wobbled_alb(photo_lin, u, v, shx, shy, S):
    PH, PW, _ = photo_lin.shape
    su = np.clip(u + shx / S, 0.0, 1.0) * (PW - 1.0)
    sv = np.clip(v + shy / S, 0.0, 1.0) * (PH - 1.0)
    x0 = np.floor(su).astype(np.int64)
    y0 = np.floor(sv).astype(np.int64)
    x1 = np.clip(x0 + 1, 0, PW - 1)
    y1 = np.clip(y0 + 1, 0, PH - 1)
    x0 = np.clip(x0, 0, PW - 1)
    y0 = np.clip(y0, 0, PH - 1)
    fx = su - np.floor(su)
    fy = sv - np.floor(sv)
    return (photo_lin[y0, x0] * (1 - fx)[:, :, None] * (1 - fy)[:, :, None]
            + photo_lin[y0, x1] * fx[:, :, None] * (1 - fy)[:, :, None]
            + photo_lin[y1, x0] * (1 - fx)[:, :, None] * fy[:, :, None]
            + photo_lin[y1, x1] * fx[:, :, None] * fy[:, :, None])


def _sample_c(Cc, u, v, shx, shy, S, N_buf):
    bu = np.mod(u + shx / S, 1.0) * N_buf - 0.5
    bv = np.mod(v + shy / S, 1.0) * N_buf - 0.5
    bx0 = np.floor(bu).astype(np.int64) % N_buf
    by0 = np.floor(bv).astype(np.int64) % N_buf
    bx1 = (bx0 + 1) % N_buf
    by1 = (by0 + 1) % N_buf
    fxb = bu - np.floor(bu)
    fyb = bv - np.floor(bv)
    return (Cc[by0, bx0] * (1 - fxb) * (1 - fyb)
            + Cc[by0, bx1] * fxb * (1 - fyb)
            + Cc[by1, bx0] * (1 - fxb) * fyb
            + Cc[by1, bx1] * fxb * fyb)


def _state_gradients(fields, lowpass_sigma_px):
    gx_in = np.asarray(fields["gx"], dtype=np.float64)
    gy_in = np.asarray(fields["gy"], dtype=np.float64)
    sig = float(lowpass_sigma_px)
    if sig > 0.0:
        from scipy import ndimage as _ndi
        n_state = gx_in.shape[0]
        sig_eff = sig * (float(n_state) / 2048.0)
        gx_s = _ndi.gaussian_filter(gx_in, sig_eff, mode="wrap")
        gy_s = _ndi.gaussian_filter(gy_in, sig_eff, mode="wrap")
    else:
        gx_s, gy_s = gx_in, gy_in
    return gx_s, gy_s


def render_frame_v5(fields, C, photo_lin, env, W, H,
                    ambient_scale=A_V3, c_layer_exposure=E_V3,
                    substrate_gain_rgb=(1.0, 1.0, 1.0),
                    beer_down=None, beer_up=None, g1=1.0,
                    blur_sigma_px=0.0, lowpass_sigma_px=B_SIGMA_V3,
                    glitter_gain=0.0, return_debug=False):
    """One v5 still.  Returns uint8 sRGB; debug dict with return_debug=True."""
    S = float(fields["S"])
    N_buf = int(np.asarray(C).shape[0])
    depth = float(env["depthM"])
    diffuse = float(env["diffuse"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    Ebase = e_base_dir_from_env(env)
    A = float(ambient_scale)
    E = float(c_layer_exposure)
    gain = np.asarray(substrate_gain_rgb, dtype=np.float64)

    # ---- 2.2 physical Beer-Lambert (replaces the ad-hoc absorb term) ----
    if beer_down is None or beer_up is None:
        from v5_physics import beer_lambert
        bl = beer_lambert(depth, env["sun"]["altitudeDeg"], env["turbid"])
        beer_down, beer_up = bl["down"], bl["up"]
    beer_down = np.asarray(beer_down, dtype=np.float64)
    beer_up = np.asarray(beer_up, dtype=np.float64)

    # ---- 2.3 turbidity scattering blur on C (mean preserving) ----
    Cc = np.asarray(C, dtype=np.float64)
    if blur_sigma_px and float(blur_sigma_px) > 0.0:
        Cc = gaussian_blur_fft(Cc, float(blur_sigma_px)).astype(np.float64)

    gx_s, gy_s = _state_gradients(fields, lowpass_sigma_px)
    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S
    gx = sample_bicubic(gx_s.astype(np.float32), px, py, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nz = 1.0 / nl

    # ---- ambient term (unchanged path; NOT modulated by C) ----
    # scat_amt is the water-scattering fraction.  It keeps the v3 definition
    # based on the v3 beer vector so the ambient level is bit-identical to the
    # v3 final composition; only the direct term's attenuation is physical.
    beer_v3 = v3_beer(env)
    scat_amt = float(1.0 - beer_v3.mean())
    amb_gain = (0.35 + 0.65 * diffuse) * (0.3 + 0.7 * sunI)
    E_amb = scatter * (scat_amt * amb_gain) * A

    PH, PW, _ = photo_lin.shape
    col = np.zeros((H, W, 3), dtype=np.float64)
    alb_all = np.zeros((H, W, 3), dtype=np.float64)
    direct_all = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        alb = _wobbled_alb(photo_lin[:, :, b:b + 1], u, v, shx, shy, S)[:, :, 0]
        alb_all[:, :, b] = alb
        c = _sample_c(Cc, u, v, shx, shy, S, N_buf)
        # 2.2 down + up, 1.3 Smith G1, 2.1 substrate gain on the direct term
        direct = (float(Ebase[b]) * c * float(beer_down[b])
                  * float(beer_up[b]) * float(g1))
        direct_all[:, :, b] = alb * float(gain[b]) * direct * strength * exposure
        bottom = alb * float(gain[b]) * (E * direct) + alb * float(E_amb[b])
        col[:, :, b] = bottom * strength

    cos_t = nz
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - diffuse) + SKY_REF[None, None, :] * diffuse
    sky_gain = (0.10 + 0.30 * (1.0 - diffuse)) * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * (sky_gain * A)

    if glitter_gain and float(glitter_gain) > 0.0:
        # 1.2 anisotropic glitter: sample the Cox-Munk (sigma_a^2, sigma_c^2)
        # Gaussian.  Off in the verification stills (the v3 offline composite
        # has no glitter term at all); implemented for the composite stage.
        from v5_physics import cox_munk_anisotropy
        an = cox_munk_anisotropy(float(env["windSpeedMS"]))
        s = np.asarray(env["sunDir"], dtype=np.float64)
        s = s / np.linalg.norm(s)
        Hv = s + np.array([0.0, 0.0, 1.0])
        Hv = Hv / np.linalg.norm(Hv)
        ndh = np.clip(-nx * Hv[0] - ny * Hv[1] + nz * Hv[2], 0.0, 1.0)
        # half-vector slope components (along/cross wind, wind along +x here)
        ma = np.sqrt(np.maximum(1.0 - ndh * ndh, 0.0))
        mc = np.zeros_like(ma)
        p = (np.exp(-0.5 * (ma * ma / max(an["sigma_a2"], 1e-9)
                            + mc * mc / max(an["sigma_c2"], 1e-9)))
             / (2.0 * np.pi * np.sqrt(an["sigma_a2"] * an["sigma_c2"])))
        spark = (np.power(ndh, 1400.0) * 2.0 + np.power(ndh, 180.0) * 0.35)
        col += sunColor[None, None, :] * (spark * p * float(glitter_gain))[:, :, None]

    lin = col * exposure
    aces = aces_narkowicz(np.maximum(lin, 0.0))
    out = (np.clip(linear_to_srgb(aces), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    if return_debug:
        return out, {"lin": lin, "aces": aces, "alb": alb_all,
                     "direct_lin": direct_all,
                     "beer_down": [float(x) for x in beer_down],
                     "beer_up": [float(x) for x in beer_up],
                     "g1": float(g1), "E": E, "A": A,
                     "blur_sigma_px": float(blur_sigma_px),
                     "substrate_gain_rgb": [float(x) for x in gain]}
    return out


def render_v3_final_repro(fields, C, photo_lin, env, W, H,
                          ambient_scale=A_V3, lowpass_sigma_px=B_SIGMA_V3,
                          c_layer_exposure=E_V3, return_debug=True):
    """Bit-level reproduction of the v3 final composition (v4 prod freeze).

    lin_E = lin_base + (E-1) * alb * direct_v3 * strength * exposure
    with direct_v3 = E_base,dir * C * beer_v3  (the old ad-hoc beer).
    """
    still, dbg = render_frame_v3(fields, C, photo_lin, env, W, H,
                                 ambient_scale=ambient_scale,
                                 lowpass_sigma_px=lowpass_sigma_px,
                                 return_debug=True)
    lin_base = dbg["lin"]
    S = float(fields["S"])
    N_buf = int(np.asarray(C).shape[0])
    depth = float(env["depthM"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    Ebase = e_base_dir_from_env(env)
    beer_v3 = v3_beer(env)
    Cc = np.asarray(C, dtype=np.float64)
    gx_s, gy_s = _state_gradients(fields, lowpass_sigma_px)
    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    gx = sample_bicubic(gx_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), u * S, v * S, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    direct_lin = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        alb = _wobbled_alb(photo_lin[:, :, b:b + 1], u, v, shx, shy, S)[:, :, 0]
        c = _sample_c(Cc, u, v, shx, shy, S, N_buf)
        direct_lin[:, :, b] = (alb * float(Ebase[b]) * c * float(beer_v3[b])
                               * strength * exposure)
    lin_E = lin_base + (float(c_layer_exposure) - 1.0) * direct_lin
    aces_E = aces_narkowicz(np.maximum(lin_E, 0.0))
    out = (np.clip(linear_to_srgb(aces_E), 0.0, 1.0) * 255.0
           + 0.5).astype(np.uint8)
    if return_debug:
        return out, {"lin": lin_E, "aces": aces_E, "lin_base": lin_base,
                     "direct_lin": direct_lin, "beer_v3": [float(x) for x in beer_v3]}
    return out
