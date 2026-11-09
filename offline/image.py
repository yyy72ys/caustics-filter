"""Image formation: same state -> photo wobble + bottom light + sky + ACES.

Per output pixel the shared water state (bicubic-sampled gradient) refracts
the view ray; the bottom sample point is shifted by
  shift_band = -N.xy * depth * (IOR_band-1)/IOR_band
(band-dependent IOR gives thin chromatic fringes).
Bottom light = albedo_linear * (E * beer + ambient), with Beer-Lambert
  beer_band = exp(-absorb_band * pathLen * 0.22),
  pathLen = (0.25+depth)*(1+turbid*1.5)   (same form as web/env.js),
and underwater-scattering ambient from the recipe scatter color.
Water surface reflection: Fresnel F = F0 + (1-F0)(1-cosT)^5, F0 = 0.02, with
sky gradient mix(sunColor, [0.62,0.70,0.76], diffuse).
ACES (Narkowicz approximation) + sRGB encode. The input photo is treated as
the bottom albedo (documented approximation).
"""
import numpy as np

IOR = np.array([1.33132, 1.33434, 1.33925], dtype=np.float64)
SKY_REF = np.array([0.62, 0.70, 0.76], dtype=np.float64)

# Dimensionless master ratio of the v2 main path (see det.py). Reference
# flat-water bottom irradiance per band, recorded in numbers/metadata.
E_BASE_DIR_DEFAULT = None  # resolved per-frame from env (sunI * sunColor)

from state import sample_bicubic, sample_bilinear  # noqa: E402


def srgb_to_linear(c):
    c = np.asarray(c, dtype=np.float64)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(c):
    c = np.clip(np.asarray(c, dtype=np.float64), 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c,
                    1.055 * (c ** (1.0 / 2.4)) - 0.055)


def aces_narkowicz(x):
    x = np.maximum(np.asarray(x, dtype=np.float64), 0.0)
    return (x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14)


def render_frame(fields, E, photo_lin, env, W, H):
    S = float(fields["S"])
    N_buf = E.shape[0]
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    diffuse = float(env["diffuse"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    absorb = np.asarray(env["absorb"], dtype=np.float64)
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)

    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S

    gx = sample_bicubic(fields["gx"], px, py, S).astype(np.float64)
    gy = sample_bicubic(fields["gy"], px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nz = 1.0 / nl

    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    beer = np.exp(-absorb * path_len * 0.22)  # (3,)
    scat_amt = float(1.0 - beer.mean())
    amb_gain = (0.35 + 0.65 * diffuse) * (0.3 + 0.7 * sunI)
    ambient = scatter * (scat_amt * amb_gain)  # (3,)

    PH, PW, _ = photo_lin.shape
    col = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
        # photo sample (clamped UV)
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
        alb = (photo_lin[y0, x0, b] * (1 - fx) * (1 - fy)
               + photo_lin[y0, x1, b] * fx * (1 - fy)
               + photo_lin[y1, x0, b] * (1 - fx) * fy
               + photo_lin[y1, x1, b] * fx * fy)
        # caustics buffer sample (periodic wrap, bilinear in texel space)
        bu = np.mod(u + shx / S, 1.0) * N_buf - 0.5
        bv = np.mod(v + shy / S, 1.0) * N_buf - 0.5
        bx0 = np.floor(bu).astype(np.int64) % N_buf
        by0 = np.floor(bv).astype(np.int64) % N_buf
        bx1 = (bx0 + 1) % N_buf
        by1 = (by0 + 1) % N_buf
        fxb = bu - np.floor(bu)
        fyb = bv - np.floor(bv)
        Eb = E[:, :, b].astype(np.float64)
        e = (Eb[by0, bx0] * (1 - fxb) * (1 - fyb)
             + Eb[by0, bx1] * fxb * (1 - fyb)
             + Eb[by1, bx0] * (1 - fxb) * fyb
             + Eb[by1, bx1] * fxb * fyb)
        col[:, :, b] = alb * (e * float(beer[b]) * strength + float(ambient[b]) * strength)

    cos_t = nz
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - diffuse) + SKY_REF[None, None, :] * diffuse
    sky_gain = (0.10 + 0.30 * (1.0 - diffuse)) * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * sky_gain

    col *= exposure
    col = aces_narkowicz(col)
    out = (np.clip(linear_to_srgb(col), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    return out


def e_base_dir_from_env(env):
    """Flat-water bottom reference irradiance per band (RGB).

    E_base,dir = sunI * sunColor. This is the denominator of the v2 master
    ratio C = E_direct / E_base,dir (its rec709 luminance is the scalar
    reference for the 1ch C). Entry surface reflectance is neglected, i.e.
    the reference is energy-conserving by construction, identical to the
    photon path's expected mean (transport.photon_frame info).
    """
    sunI = float(env["sunI"])
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    return sunI * sunColor


def render_frame_v2(fields, C, photo_lin, env, W, H):
    """Image formation with the OUTPUT-SPEC section 3 separation (v2).

      E_bottom(q,t) = E_amb(d) + E_base,dir(d) * C(q,t)
      L_camera = T_up * (rho/pi) * E_bottom + L_in-scatter + L_reflection

    C is the dimensionless 1ch master ratio (achromatic redistribution);
    the ambient term E_amb is NOT modulated by C (this was the v1 failure
    mode: the whole photo, ambient included, was scaled by a field of mean
    ~0.59). Downward transmission (Beer-Lambert), scattering ambient,
    Fresnel sky reflection, exposure and ACES are the conventional terms
    reused from render_frame (v1); only the composition order follows the
    spec. C is sampled with the same periodic bilinear wrap as v1.
    """
    S = float(fields["S"])
    N_buf = C.shape[0]
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    diffuse = float(env["diffuse"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    absorb = np.asarray(env["absorb"], dtype=np.float64)
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    Ebase = e_base_dir_from_env(env)  # (3,)

    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S

    gx = sample_bicubic(fields["gx"], px, py, S).astype(np.float64)
    gy = sample_bicubic(fields["gy"], px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nz = 1.0 / nl

    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    beer = np.exp(-absorb * path_len * 0.22)  # (3,)
    scat_amt = float(1.0 - beer.mean())
    amb_gain = (0.35 + 0.65 * diffuse) * (0.3 + 0.7 * sunI)
    E_amb = scatter * (scat_amt * amb_gain)  # (3,) ambient bottom irradiance

    Cc = np.asarray(C, dtype=np.float64)
    PH, PW, _ = photo_lin.shape
    col = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
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
        alb = (photo_lin[y0, x0, b] * (1 - fx) * (1 - fy)
               + photo_lin[y0, x1, b] * fx * (1 - fy)
               + photo_lin[y1, x0, b] * (1 - fx) * fy
               + photo_lin[y1, x1, b] * fx * fy)
        # master ratio sample (periodic wrap, bilinear in texel space)
        bu = np.mod(u + shx / S, 1.0) * N_buf - 0.5
        bv = np.mod(v + shy / S, 1.0) * N_buf - 0.5
        bx0 = np.floor(bu).astype(np.int64) % N_buf
        by0 = np.floor(bv).astype(np.int64) % N_buf
        bx1 = (bx0 + 1) % N_buf
        by1 = (by0 + 1) % N_buf
        fxb = bu - np.floor(bu)
        fyb = bv - np.floor(bv)
        c = (Cc[by0, bx0] * (1 - fxb) * (1 - fyb)
             + Cc[by0, bx1] * fxb * (1 - fyb)
             + Cc[by1, bx0] * (1 - fxb) * fyb
             + Cc[by1, bx1] * fxb * fyb)
        direct = float(Ebase[b]) * c * float(beer[b])
        bottom = direct + float(E_amb[b])  # ambient NOT modulated by C
        col[:, :, b] = alb * bottom * strength

    cos_t = nz
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - diffuse) + SKY_REF[None, None, :] * diffuse
    sky_gain = (0.10 + 0.30 * (1.0 - diffuse)) * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * sky_gain

    col *= exposure
    col = aces_narkowicz(col)
    out = (np.clip(linear_to_srgb(col), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    return out


def render_frame_v3(fields, C, photo_lin, env, W, H, ambient_scale=1.0,
                    lowpass_sigma_px=0.0, return_debug=False):
    """Image formation for the PIPELINE-V3 composition round (v2 state frozen).

    Same OUTPUT-SPEC section 3 separation as render_frame_v2, plus the two
    v3 variant axes (composition stage only; state generation and transport
    are untouched):

      A (ambient scale): E_amb and the camera-side additive term are
        multiplied by ``ambient_scale``. The direct term
        E_base,dir(d) * C(q,t) is bit-identical to v2. NOTE on the
        "in-scatter" wording: this offline path has no separate L_in-scatter
        buffer -- underwater scattering is folded into E_amb, and the only
        other C-independent additive is the Fresnel sky term (L_reflection).
        So A scales E_amb AND the sky-reflection additive; both are
        documented per variant. Direct light is never scaled.
      B (displacement band limit): when ``lowpass_sigma_px`` > 0, the state
        gradients (gx, gy) -- and identically h, though h is unused by image
        formation -- are smoothed with a periodic (wrap) Gaussian of that
        sigma (state-grid pixels; v3 spec value is 3.0 px at 2048^2) BEFORE
        the normal/shift computation. Albedo wobble and C sampling share the
        same smoothed shifts, so they stay coherent.

    Returns uint8 sRGB PNG pixels; with return_debug=True also returns a dict
    with lin (pre-ACES linear incl. exposure, float64 HxWx3), aces (post-ACES
    pre-sRGB-encode, float64 HxWx3) and alb (wobbled albedo linear, HxWx3).
    """
    S = float(fields["S"])
    N_buf = np.asarray(C).shape[0]
    depth = float(env["depthM"])
    turbid = float(env["turbid"])
    diffuse = float(env["diffuse"])
    sunI = float(env["sunI"])
    strength = float(env["strength"])
    exposure = float(env["exposure"])
    absorb = np.asarray(env["absorb"], dtype=np.float64)
    scatter = np.asarray(env["scatter"], dtype=np.float64)
    sunColor = np.asarray(env["sunColor"], dtype=np.float64)
    Ebase = e_base_dir_from_env(env)  # (3,)
    A = float(ambient_scale)

    gx_in = np.asarray(fields["gx"], dtype=np.float64)
    gy_in = np.asarray(fields["gy"], dtype=np.float64)
    sig = float(lowpass_sigma_px)
    if sig > 0.0:
        from scipy import ndimage as _ndi
        # Spec sigma is quoted at 2048^2; rescale if the state grid differs.
        n_state = gx_in.shape[0]
        sig_eff = sig * (float(n_state) / 2048.0)
        gx_s = _ndi.gaussian_filter(gx_in, sig_eff, mode="wrap")
        gy_s = _ndi.gaussian_filter(gy_in, sig_eff, mode="wrap")
    else:
        gx_s = gx_in
        gy_s = gy_in

    iy, ix = np.mgrid[0:H, 0:W]
    u = (ix.astype(np.float64) + 0.5) / W
    v = (iy.astype(np.float64) + 0.5) / H
    px = u * S
    py = v * S

    # NOTE: the lowpass acts in state space (smoothed state fields resampled
    # to the output lattice with the same bicubic path as v2), not in screen
    # space.
    gx = sample_bicubic(gx_s.astype(np.float32), px, py, S).astype(np.float64)
    gy = sample_bicubic(gy_s.astype(np.float32), px, py, S).astype(np.float64)
    nl = np.sqrt(gx * gx + gy * gy + 1.0)
    nx = -gx / nl
    ny = -gy / nl
    nz = 1.0 / nl

    path_len = (0.25 + depth) * (1.0 + turbid * 1.5)
    beer = np.exp(-absorb * path_len * 0.22)  # (3,)
    scat_amt = float(1.0 - beer.mean())
    amb_gain = (0.35 + 0.65 * diffuse) * (0.3 + 0.7 * sunI)
    E_amb = scatter * (scat_amt * amb_gain) * A  # (3,) scaled ambient

    Cc = np.asarray(C, dtype=np.float64)
    PH, PW, _ = photo_lin.shape
    col = np.zeros((H, W, 3), dtype=np.float64)
    alb_all = np.zeros((H, W, 3), dtype=np.float64)
    for b in range(3):
        Cb = (float(IOR[b]) - 1.0) / float(IOR[b])
        shx = -nx * depth * Cb
        shy = -ny * depth * Cb
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
        alb = (photo_lin[y0, x0, b] * (1 - fx) * (1 - fy)
               + photo_lin[y0, x1, b] * fx * (1 - fy)
               + photo_lin[y1, x0, b] * (1 - fx) * fy
               + photo_lin[y1, x1, b] * fx * fy)
        alb_all[:, :, b] = alb
        bu = np.mod(u + shx / S, 1.0) * N_buf - 0.5
        bv = np.mod(v + shy / S, 1.0) * N_buf - 0.5
        bx0 = np.floor(bu).astype(np.int64) % N_buf
        by0 = np.floor(bv).astype(np.int64) % N_buf
        bx1 = (bx0 + 1) % N_buf
        by1 = (by0 + 1) % N_buf
        fxb = bu - np.floor(bu)
        fyb = bv - np.floor(bv)
        c = (Cc[by0, bx0] * (1 - fxb) * (1 - fyb)
             + Cc[by0, bx1] * fxb * (1 - fyb)
             + Cc[by1, bx0] * (1 - fxb) * fyb
             + Cc[by1, bx1] * fxb * fyb)
        direct = float(Ebase[b]) * c * float(beer[b])  # NOT scaled by A
        bottom = direct + float(E_amb[b])  # ambient scaled by A
        col[:, :, b] = alb * bottom * strength

    cos_t = nz
    F = 0.02 + 0.98 * np.power(np.clip(1.0 - cos_t, 0.0, 1.0), 5.0)
    sky = sunColor[None, None, :] * (1.0 - diffuse) + SKY_REF[None, None, :] * diffuse
    sky_gain = (0.10 + 0.30 * (1.0 - diffuse)) * (0.25 + 0.75 * min(max(sunI + 0.25, 0.0), 1.0))
    col += sky * F[:, :, None] * (sky_gain * A)  # camera-side additive scaled by A

    lin = col * exposure
    aces = aces_narkowicz(np.maximum(lin, 0.0))
    out = (np.clip(linear_to_srgb(aces), 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    if return_debug:
        return out, {"lin": lin, "aces": aces, "alb": alb_all,
                     "ambient_scale": A, "lowpass_sigma_px": sig}
    return out
