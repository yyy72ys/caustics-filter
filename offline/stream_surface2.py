"""stream_surface2 -- BRIEF-V4-1-2 flowing-water surface generator (new file).

v4-1.1 (`stream_surface.py`) is left untouched; this is the v4-1.2 redesign.

STRUCTURE (the final form fixed by BRIEF-V4-1-2 + RALLY-V4-1-1 round-2)

(i) DIRECT BED RESPONSE (stationary)
    eta_i_hat = T(k) * zeta_hat, with the 1.1 gravity transfer
      T(k) = [cosh(kH) - k g sinh(kH)/(kx^2 U^2)]^-1,
    and the 1.1 breaking-limit safety valve kept (resonant band |denominator|
    < 0.30 uniformly scaled so its rms slope contribution <= 0.30; the
    dropped energy fraction is reported).
    The bed zeta is either
      - "grid": a lattice floor.  zeta is non-zero only on the lattice lines,
        so its spectrum is a comb at k >= 2*pi/pitch where |T| is
        astronomically small.  The response is therefore *measured* to be
        ~0 (never switched off by hand), or
      - "natural": a Nikora-type Gaussian random-phase field with
        S_bed(k) = sigma_b^2 ell_b^2/(2 pi) exp(-k^2 ell_b^2/2)
        (D84 mapping sigma_b = D84/4, ell_b = D84, as in 1.1).
    This component is STATIC (it is the steady surface imprint of the bed).

(ii) STATIONARY RING + TAIL (Dolcetti & Krynkin 2017 spectral model)
    k0 = smaller root of the still-water stationary condition
        k0 U0 = sqrt[(g k0 + gamma k0^3/rho) tanh(k0 H)]      (gamma=0.0728,
        rho=1000, U0 = surface velocity).
    If U0 < c_min(H) (the gravity-capillary minimum phase speed) there is no
    ring at all: k0 = None, ring fraction = 0.
    Ring shape: radial Gaussian of relative width RING_RELWIDTH around k0,
      directional Poisson-type factor  p(beta) = cosh(lam cos beta)/(2 pi I0(lam))
      (two lobes at beta = 0, pi; lam = DIR_LAMBDA), i.e. the Poisson kernel
      symmetrised to two lobes.
    Tail: isotropic k^-alpha for k >= k0 (alpha = 5, 6, 7 are all computed;
    the rendered deliverable uses alpha = 6, the midpoint of the bracket --
    NOT a fitted value).
    Ring fraction (derived, NOT an input):
      var_ring/var_tail = G * (delta/k0) * (alpha-2) / (sqrt(2 pi) (1-rho^(2-alpha)))
      with G = RING_GAIN (ring spectral density at k0 relative to the tail
      level there), delta = RING_RELWIDTH*k0, rho = kmax/k0, kmax = pi/dx.
      phi = R/(1+R).
    Time evolution: pure advection, omega = kx*Us, so the whole (ii) component
    loops exactly over 8 s when L = 8*Us/n.

(iii) TURBULENT COMPONENT (non-stationary)
    Horoshenkov-type autocovariance W(rho) = exp(-rho^2/2 sigma_w^2)
      cos(2 pi rho/L0)  (the 1/2 convention is kept), sigma_w = 0.22 m is a
      CORRELATION parameter, not an amplitude (1.1 convention).
    L0 from the correlation
      L0 = ks * 0.045 * exp(0.54 * Vbar/u*)       (Vbar = depth-mean velocity)
      and the alternative L0/H in 1.8..2.9 are both computed; L0 is capped at
      2.5 H.
    Spectrum = the periodic discrete correspondence of W on the grid
      (sample W, FFT, clip at 0), i.e. the 1.1 construction.
    Time evolution: omega = kx*Us + omega0(k) plus the OU update of every
      modal complex amplitude with tau = H/u* (stationary-variance preserving,
      the 1.1 form).

VARIANCE CLOSURE (no per-condition adjustment anywhere)
    sigma_tot = a * U0^b * H^c * s^d  (ONE log-log regression over the 13
    Dolcetti-2016 conditions; the form is fixed in advance -- see the driver).
    var(ii)  = phi * (sigma_tot^2 - var_i)
    var(iii) = (1-phi) * (sigma_tot^2 - var_i)
    var_i stays at its PHYSICAL value (it is never rescaled to fit a budget).
    Inside (ii): ring = phi*var(ii), tail = (1-phi)*var(ii).
    Every shape is generated from its model spectrum and then normalised by
    its realised variance, so "shape from the model, amplitude from the
    budget" holds literally.

U and U0
    Scene input is the depth-mean U; U0 = 1.15*U (the 1.1 log-profile
    assumption) is the surface velocity used by the stationary condition.

Conventions: periodic square patch Lx = Ly = L, Nx = Ny = N.  All spectral
quantities live on the numpy fftfreq grid; dk = 2 pi/L; the (0,0) mode is
zeroed so every field is zero-mean.
"""

import numpy as np

G = 9.81
SIGMA_K = 0.0728          # N/m  (surface tension, as fixed by the brief)
RHO = 1000.0              # kg/m^3
US_FACTOR = 1.15          # U0 = 1.15 U  (1.1 assumption, kept and declared)
SLOPE_CAP = 0.30          # breaking-limit safety valve for (i)
RES_BAND_DENOM = 0.30     # |denominator| < this  ->  resonant band

D84_SIGMA_DIVISOR = 4.0
D84_ELL_FACTOR = 1.0

RING_GAIN = 45.0          # ring spectral density at k0 / tail level there
RING_RELWIDTH = 0.15      # delta/k0
DIR_LAMBDA = 4.0          # Poisson-kernel concentration (two lobes at 0, pi)
ALPHA_BRACKET = (5, 6, 7)
ALPHA_MAIN = 6

SIGMA_W_REF = 0.22        # m, Horoshenkov correlation parameter (1.1 value)
L0_CORR_C = 0.045
L0_CORR_K = 0.54
L0_CAP_FACTOR = 2.5

# ---------------------------------------------------------------------------
# Dolcetti et al. (2016) Table II, verbatim from BRIEF-V4-1-2.
# columns: cond, H(mm), slope, U0(m/s), F, sigma(mm), k0_table(rad/m)
# ---------------------------------------------------------------------------
TABLE13 = [
    (1, 42.2, 0.001, 0.19, 0.30, 0.05, None),
    (2, 72.9, 0.001, 0.35, 0.41, 0.40, 89.7),
    (3, 101.0, 0.001, 0.41, 0.41, 0.50, 63.6),
    (4, 42.2, 0.002, 0.30, 0.47, 0.25, 131.5),
    (5, 101.3, 0.002, 0.49, 0.49, 1.79, 45.1),
    (6, 43.0, 0.003, 0.34, 0.52, 0.49, 99.7),
    (7, 73.1, 0.002, 0.46, 0.54, 1.21, 52.4),
    (8, 40.5, 0.004, 0.36, 0.57, 0.34, 89.5),
    (9, 43.4, 0.005, 0.40, 0.61, 0.46, 72.6),
    (10, 99.0, 0.003, 0.60, 0.61, 2.03, 31.4),
    (11, 72.4, 0.003, 0.54, 0.64, 1.17, 39.4),
    (12, 43.1, 0.006, 0.43, 0.66, 0.57, 63.8),
    (13, 73.2, 0.004, 0.58, 0.68, 1.10, 34.8),
]
FLUME_D = 0.0254          # mono-sphere bed diameter of the 2016 flume (1.1)
FLUME_MEASURED_PEAK = {7: 60.0, 13: -40.0}   # rad/m, on-axis kx anchors


def table_rows():
    """The 13 conditions as SI records."""
    out = []
    for (c, h_mm, s, u0, F, sg_mm, k0t) in TABLE13:
        out.append({"cond": c, "H": h_mm / 1000.0, "slope": s, "U0": u0,
                    "F": F, "sigma": sg_mm / 1000.0, "k0_table": k0t})
    return out


# ---------------------------------------------------------------------------
# basic hydraulics
# ---------------------------------------------------------------------------
def derive(U, H, S0):
    U = float(U)
    H = float(H)
    ustar = float(np.sqrt(G * H * float(S0))) if H > 0 and S0 > 0 else 0.0
    tau = float(H / ustar) if ustar > 0 else float("inf")
    return {"Fr": float(U / np.sqrt(G * H)), "ustar": ustar, "tau": tau,
            "Us": float(US_FACTOR * U), "U0": float(US_FACTOR * U),
            "formula_ustar": "u* = sqrt(g H S0) (rectangular approximation)",
            "formula_Us": "U0 = 1.15 U (log-profile surface factor, "
                          "inherited from v4-1.1)"}


def omega0(k, H):
    k = np.asarray(k, dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        w2 = (G * k + SIGMA_K * k ** 3 / RHO) * np.tanh(k * float(H))
    return np.sqrt(np.maximum(w2, 0.0))


def c_min_phase(H):
    """Minimum gravity-capillary phase speed omega/k (and the k where it sits)."""
    k = np.logspace(-2, 4, 40000)
    c = omega0(k, H) / k
    j = int(np.argmin(c))
    return float(c[j]), float(k[j])


def k0_stationary(U0, H, kmax=20000.0, n=400000):
    """Smaller (gravity-branch) root of k U0 = sqrt((gk+g k^3/rho) tanh(kH)).

    Returns None when U0 is at or below the minimum phase speed (no root).
    """
    U0 = float(U0)
    H = float(H)
    cm, _km = c_min_phase(H)
    if not np.isfinite(U0) or U0 <= cm:
        return None
    k = np.logspace(-3.0, np.log10(kmax), n)
    f = k * U0 - omega0(k, H)
    s = np.sign(f)
    idx = np.where(s[:-1] * s[1:] < 0.0)[0]
    if idx.size == 0:
        return None
    roots = []
    for j in idx[:4]:
        a, b = float(k[j]), float(k[j + 1])
        fa = a * U0 - float(omega0(a, H))
        for _ in range(200):
            m = 0.5 * (a + b)
            fm = m * U0 - float(omega0(m, H))
            if fa * fm <= 0.0:
                b = m
            else:
                a, fa = m, fm
        roots.append(0.5 * (a + b))
    roots.sort()
    return float(roots[0])


# ---------------------------------------------------------------------------
# sigma regression (single fit over the 13 conditions; form fixed in advance)
# ---------------------------------------------------------------------------
def _design(rows, form):
    X, y = [], []
    for r in rows:
        if form == "u0hs":
            X.append([1.0, np.log(r["U0"]), np.log(r["H"]),
                      np.log(r["slope"])])
        elif form == "ustar":
            u2g = G * r["H"] * r["slope"] / G          # = u*^2/g  (m)
            X.append([1.0, np.log(u2g), np.log(r["H"] / FLUME_D)])
        else:
            raise ValueError("unknown regression form %r" % (form,))
        y.append(np.log(r["sigma"]))
    return np.array(X), np.array(y)


def _predict(beta, form, r):
    if form == "u0hs":
        return float(np.exp(beta[0] + beta[1] * np.log(r["U0"])
                            + beta[2] * np.log(r["H"])
                            + beta[3] * np.log(r["slope"])))
    u2g = G * r["H"] * r["slope"] / G
    return float(np.exp(beta[0] + beta[1] * np.log(u2g)
                        + beta[2] * np.log(r["H"] / FLUME_D)))


def fit_sigma(rows=None, form="u0hs"):
    rows = table_rows() if rows is None else rows
    X, y = _design(rows, form)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    per = []
    for r in rows:
        p = _predict(beta, form, r)
        per.append({"cond": r["cond"], "sigma_measured": r["sigma"],
                    "sigma_model": p, "ratio": p / r["sigma"]})
    loo = []
    for i, r in enumerate(rows):
        sub = [q for j, q in enumerate(rows) if j != i]
        Xs, ys = _design(sub, form)
        bs, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
        p = _predict(bs, form, r)
        loo.append({"cond": r["cond"], "sigma_loo": p,
                    "ratio_loo": p / r["sigma"]})
    ratios = [d["ratio"] for d in per]
    lratios = [d["ratio_loo"] for d in loo]
    n = len(beta) - 1
    yhat = X @ beta
    resid = y - yhat
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "form": form,
        "form_expression": (
            "sigma = a * U0^b * H^c * s^d" if form == "u0hs"
            else "sigma = a * (u*^2/g)^p * (H/d)^q, d=0.0254 m"),
        "coefficients": {"a": float(np.exp(beta[0]))},
        "exponents": [float(x) for x in beta[1:]],
        "r2_log": float(1.0 - ss_res / ss_tot),
        "n_points": int(len(rows)), "n_parameters": int(n + 1),
        "per_condition": per, "loo": loo,
        "ratio_range_in_sample": [float(min(ratios)), float(max(ratios))],
        "ratio_range_loo": [float(min(lratios)), float(max(lratios))],
        "note": ("in-sample ratios are 1.000 by construction (the model sigma "
                 "IS the fitted sigma); the LOO column is the honest "
                 "generalisation estimate"),
    }


def sigma_from_fit(fit, U0, H, slope):
    b = [np.log(fit["coefficients"]["a"])] + list(fit["exponents"])
    return float(np.exp(b[0] + b[1] * np.log(U0) + b[2] * np.log(H)
                        + b[3] * np.log(slope)))


# ---------------------------------------------------------------------------
# grids / spectra
# ---------------------------------------------------------------------------
def wavenumbers(N, L):
    f = np.fft.fftfreq(int(N), d=float(L) / float(int(N)))
    k1 = f * 2.0 * np.pi
    kx, ky = np.meshgrid(k1, k1, indexing="xy")
    kappa = np.sqrt(kx * kx + ky * ky)
    return kx, ky, kappa, 2.0 * np.pi / float(L)


def kmax_of(N, L):
    """Largest |k| represented on the N x N grid of a patch of length L."""
    dk = 2.0 * np.pi / float(L)
    return float(np.sqrt(2.0) * (dk * (int(N) // 2)))


def hermitian_unit(N, rng):
    """Unit-variance Hermitian complex Gaussian field (CN(0,1), E|Z|^2=1).

    Note: the conjugation partner of index i is (-i) % N, NOT N-1-i.  Using
    np.flip here (as v4-1.1 did) silently breaks Hermiticity, which makes the
    inverse FFT complex and makes the realised variance disagree with
    sum|hat|^2/N^4.  Fixed in v4-1.2.
    """
    Z = (rng.standard_normal((N, N))
         + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    idx = (-np.arange(N)) % N
    Zr = Z[np.ix_(idx, idx)]
    return (Z + np.conj(Zr)) / np.sqrt(2.0)


def transfer_T(kx, kappa, H, U):
    """The 1.1 gravity transfer, unchanged (surface tension is NOT in T)."""
    kx = np.asarray(kx, dtype=np.float64)
    kappa = np.asarray(kappa, dtype=np.float64)
    kh = kappa * float(H)
    tanh = np.tanh(np.minimum(kh, 50.0))
    T = np.zeros_like(kh)
    nz = np.abs(kx) > 0.0
    safe = nz & (kh <= 50.0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cosh = np.cosh(np.minimum(kh, 50.0))
        denom = np.zeros_like(kh)
        f = kappa * G * tanh / np.maximum(kx * kx * float(U) ** 2, 1e-300)
        denom[safe] = cosh[safe] * (1.0 - f[safe])
        T[safe] = 1.0 / denom[safe]
    T[~nz] = 0.0
    with np.errstate(over="ignore"):
        denom_out = np.full_like(kh, np.inf)
        denom_out[safe] = denom[safe]
        denom_out[~nz] = np.inf
    return T, denom_out


def bed_spectrum(kappa, sigma_b, ell_b):
    kappa = np.asarray(kappa, dtype=np.float64)
    return (float(sigma_b) ** 2 * float(ell_b) ** 2 / (2.0 * np.pi)
            * np.exp(-0.5 * (kappa * float(ell_b)) ** 2))


def bed_mapping_D84(D84):
    return {"sigma_b": float(D84) / D84_SIGMA_DIVISOR,
            "ell_b": float(D84) * D84_ELL_FACTOR,
            "formula": "sigma_b = D84/4, ell_b = D84 (1.1 stated assumption)"}


def bed_mapping_sphere(d):
    return {"sigma_b": float(d) / 4.0, "ell_b": float(d),
            "formula": "sigma_b = d/4, ell_b = d (1.1 stated assumption)"}


def ring_fraction(k0, alpha, kmax, ring_gain=RING_GAIN,
                  relwidth=RING_RELWIDTH):
    """phi = var_ring/(var_ring+var_tail) derived from the model shape."""
    if k0 is None or k0 <= 0:
        return 0.0, {"R_over_T": 0.0, "reason": "no ring (U0 <= c_min)"}
    rho = float(kmax) / float(k0)
    delta = relwidth * float(k0)
    amp = 1.0 - rho ** (2.0 - alpha) if rho > 1.0 else 0.0
    if amp <= 0:
        amp = 1.0
    R = (ring_gain * (float(alpha) - 2.0) * (delta / float(k0))
         / (np.sqrt(2.0 * np.pi) * amp))
    phi = R / (1.0 + R)
    return float(phi), {"R_over_T": float(R), "rho_kmax_over_k0": rho,
                        "delta": delta, "alpha": float(alpha),
                        "ring_gain": ring_gain, "relwidth": relwidth}


def turb_L0(H, U, ustar, ks):
    """L0 from the correlation law + the L0/H=1.8..2.9 band, capped at 2.5H."""
    H = float(H)
    ks = float(ks)
    Vbar = float(U)
    if ustar > 0:
        arg = min(L0_CORR_K * Vbar / ustar, 700.0)
        L0c = ks * L0_CORR_C * float(np.exp(arg))
    else:
        L0c = float("nan")
    cap = L0_CAP_FACTOR * H
    L0 = float(min(L0c, cap))
    return {
        "L0_correlation_m": float(L0c),
        "L0_over_H_correlation": float(L0c / H) if np.isfinite(L0c) else None,
        "L0_over_H_expected_band": [1.8, 2.9],
        "L0_over_H_band_m": [1.8 * H, 2.9 * H],
        "cap_2p5H_m": float(cap),
        "L0_used_m": L0,
        "cap_applied": bool(L0c > cap),
        "formula": ("L0 = ks*0.045*exp(0.54*Vbar/u*), Vbar = depth-mean U; "
                    "capped at 2.5H"),
    }


def omega0_field(kappa, H):
    return omega0(kappa, H)


# ---------------------------------------------------------------------------
# component builders
# ---------------------------------------------------------------------------
def _grid_floor_zeta(N, L, pitch, half_width, amp):
    """Lattice floor: zeta = amp on the bar lines, 0 elsewhere (deterministic).

    The pitch is snapped to L/k with k a divisor of N, so the sampled pattern
    is EXACTLY periodic on the grid and its spectrum is a clean comb at
    2*pi*m/p (no window leakage -- leakage would fake low-k bed energy).
    """
    dx = L / float(N)
    divs = [d for d in range(4, int(N) + 1) if int(N) % d == 0]
    k = min(divs, key=lambda d: abs(L / float(d) - float(pitch)))
    p = L / float(k)
    m = max(1, int(round(float(half_width) / dx)))
    w = m * dx
    x = (np.arange(int(N)) + 0.5) * dx
    inbar = (x % p) < w
    z = np.where(inbar[:, None] | inbar[None, :], float(amp), 0.0)
    return z, {"pitch_snapped_m": p, "periods_per_patch": int(k),
               "cells_per_period": int(N) // int(k),
               "bar_width_m": float(w), "amplitude_m": float(amp),
               "first_comb_k_rad_per_m": 2.0 * np.pi / p,
               "pitch_target_m": float(pitch)}


class StreamState2:
    """v4-1.2 stream surface.  fields_at(t) returns (i)+(ii)+(iii)."""

    def __init__(self, N=512, L=4.0, H=0.35, U=0.8, S0=0.002, D84=0.09,
                 seed=1, sigma_b=None, ell_b=None, bed="natural",
                 alpha=ALPHA_MAIN, sigma_total=None, fit=None,
                 grid_pitch=None, grid_bar_halfwidth=0.002, grid_amp=None,
                 ring_gain=RING_GAIN, ring_relwidth=RING_RELWIDTH,
                 dir_lambda=DIR_LAMBDA, sigma_w=SIGMA_W_REF,
                 l0_override=None, dt_loop=8.0):
        self.N = int(N)
        self.L = float(L)
        self.S = float(L)
        self.H = float(H)
        self.U = float(U)
        self.S0 = float(S0)
        self.seed = int(seed)
        self.alpha = float(alpha)
        self.bed_kind = str(bed)
        d = derive(self.U, self.H, self.S0)
        self.derived = d
        self.Us = d["Us"]
        self.U0 = d["U0"]
        self.ustar = d["ustar"]
        self.tau = d["tau"]

        # --- stationary ring geometry -------------------------------------
        self.cmin, self.cmin_k = c_min_phase(self.H)
        self.k0 = k0_stationary(self.U0, self.H)
        self.kx, self.ky, self.kappa, self.dk = wavenumbers(self.N, self.L)
        self.kmax = kmax_of(self.N, self.L)
        self.phi, self.phi_info = ring_fraction(
            self.k0, self.alpha, self.kmax, ring_gain, ring_relwidth)
        self.ring_gain = float(ring_gain)
        self.ring_relwidth = float(ring_relwidth)
        self.dir_lambda = float(dir_lambda)

        # --- turbulent length scales --------------------------------------
        self.ks = float(D84)
        self.D84 = float(D84)
        l0 = turb_L0(self.H, self.U, self.ustar, self.ks)
        if l0_override is not None:
            l0["L0_used_m"] = float(l0_override)
            l0["override"] = True
        self.L0_info = l0
        self.L0 = float(l0["L0_used_m"])
        self.sigma_w = float(sigma_w)

        # --- variance budget ----------------------------------------------
        if sigma_total is None:
            fit = fit_sigma() if fit is None else fit
            self.sigma_tot = sigma_from_fit(fit, self.U0, self.H, self.S0)
        else:
            self.sigma_tot = float(sigma_total)
        self.fit_used = (fit is not None)

        # --- (i) direct bed response (PHYSICAL amplitude, never rescaled) --
        self.eta_i_hat, self.info_i = self._build_i(sigma_b, ell_b)
        self.var_i = self._var_of(self.eta_i_hat)

        # --- budget split --------------------------------------------------
        vtot = self.sigma_tot ** 2
        rem = max(vtot - self.var_i, 0.0)
        self.var_ii = self.phi * rem
        self.var_iii = (1.0 - self.phi) * rem
        self.var_ring = self.phi * self.var_ii
        self.var_tail = (1.0 - self.phi) * self.var_ii
        if self.k0 is None:
            # no ring: the whole (ii) share is the tail (alpha power law from
            # the lowest resolved wavenumber)
            self.var_ring = 0.0
            self.var_tail = self.var_ii

        rng = np.random.default_rng(self.seed)
        self.hat_ring = self._shape_to_var(
            self._ring_shape(), self.var_ring, rng, "ring")
        self.hat_tail = self._shape_to_var(
            self._tail_shape(), self.var_tail, rng, "tail")
        self.hat_ii0 = self.hat_ring + self.hat_tail
        self.var_ii_real = self._var_of(self.hat_ii0)

        # (iii) two independent draws with the SAME modal spectrum: the OU
        # update mixes them (A0 decaying in, B0 relaxing in) so the stationary
        # variance is preserved exactly.
        shape_iii = self._turb_shape()
        rngA = np.random.default_rng(self.seed + 0xC0A1)
        rngB = np.random.default_rng(self.seed + 0xB0B1)
        self.hat_turbA = self._normalise_hat(
            shape_iii * hermitian_unit(self.N, rngA), self.var_iii, "turbA")
        self.hat_turbB = self._normalise_hat(
            shape_iii * hermitian_unit(self.N, rngB), self.var_iii, "turbB")
        self.var_iii_real = self._var_of(self.hat_turbA)

        # --- time evolution ------------------------------------------------
        self.om_ii = self.kx * self.Us
        self.om_iii = self.kx * self.Us + omega0_field(self.kappa, self.H)
        self.dt_loop = float(dt_loop)

    # -- internals ---------------------------------------------------------
    def _var_of(self, hat):
        return float(np.mean(np.real(np.fft.ifft2(hat)) ** 2))

    def _normalise_hat(self, hat, var, name=""):
        v = self._var_of(hat)
        if var <= 0:
            return np.zeros_like(hat)
        if v <= 0:
            raise ValueError("degenerate %s shape (zero variance)" % name)
        return hat * np.sqrt(var / v)

    def _shape_to_var(self, shape_hat, var, rng, name):
        hat = shape_hat * hermitian_unit(self.N, rng)
        return self._normalise_hat(hat, var, name)

    def _build_i(self, sigma_b, ell_b):
        if self.bed_kind == "grid":
            pitch = (0.0254 if self.ks <= 0 else self.ks)
            amp = (self.ks if self.ks > 0 else pitch)
            z, gi = _grid_floor_zeta(self.N, self.L, pitch,
                                     0.002, amp)
            zeta_hat = np.fft.fft2(z)
            gi["kind"] = "grid (lattice floor: zeta non-zero on bar lines only)"
        else:
            if sigma_b is None or ell_b is None:
                m = bed_mapping_D84(self.ks)
                sigma_b, ell_b = m["sigma_b"], m["ell_b"]
            else:
                m = {"sigma_b": float(sigma_b), "ell_b": float(ell_b),
                     "formula": "explicit sigma_b/ell_b (caller-stated)"}
            rng = np.random.default_rng(self.seed + 0x5EED)
            S = bed_spectrum(self.kappa, sigma_b, ell_b)
            amp = float(self.N ** 2) * np.sqrt(np.maximum(S * self.dk ** 2, 0))
            zeta_hat = amp * hermitian_unit(self.N, rng)
            zeta_hat[0, 0] = 0.0
            gi = dict(m)
            gi["kind"] = "natural (Nikora-type Gaussian random-phase bed)"
        T, denom = transfer_T(self.kx, self.kappa, self.H, self.U)
        eta_hat = T * zeta_hat        # breaking-limit safety valve (1.1 rule, kept)
        band = np.abs(denom) < RES_BAND_DENOM
        band[~np.isfinite(denom)] = False
        band[np.abs(self.kx) == 0.0] = False
        k2 = self.kappa ** 2
        N2 = self.N ** 2
        tot_e = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
        band_e = float(np.sum(np.abs(eta_hat[band]) ** 2) / N2 ** 2)
        band_slope = float(np.sqrt(max(
            np.sum(k2[band] * np.abs(eta_hat[band]) ** 2) / N2 ** 2, 0.0)))
        scale = 1.0
        if band_slope > SLOPE_CAP > 0:
            scale = SLOPE_CAP / band_slope
            eta_hat[band] = eta_hat[band] * scale
        after = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
        info = {
            "bed": gi,
            "_T": T,
            "_zeta_hat": zeta_hat,
            "T_max": float(np.max(np.abs(T))),
            "resonant_band_def": "|denominator| < %.3g" % RES_BAND_DENOM,
            "band_rms_slope_before": band_slope,
            "band_energy_fraction": (band_e / tot_e if tot_e > 0 else 0.0),
            "scale_applied": float(scale),
            "dropped_band_pct": float(100.0 * (1.0 - scale ** 2)),
            "dropped_total_pct": float(100.0 * (tot_e - after) / tot_e
                                       if tot_e > 0 else 0.0),
            "sigma_i_m": float(np.sqrt(max(after, 0.0))),
        }
        return eta_hat, info

    def _ring_shape(self):
        if self.k0 is None:
            return np.zeros((self.N, self.N), dtype=np.complex128)
        k0 = self.k0
        delta = self.ring_relwidth * k0
        rad = np.exp(-0.5 * ((self.kappa - k0) / delta) ** 2)
        beta = np.arctan2(self.ky, self.kx)
        i0 = float(np.i0(self.dir_lambda))
        dirn = np.cosh(self.dir_lambda * np.cos(beta)) / (2.0 * np.pi * i0)
        S = rad * dirn
        N2 = float(self.N ** 2)
        return N2 * np.sqrt(np.maximum(S * self.dk ** 2, 0.0))

    def _tail_shape(self):
        klo = self.k0 if self.k0 is not None else 2.0 * self.dk
        S = np.zeros_like(self.kappa)
        m = self.kappa >= klo
        S[m] = self.kappa[m] ** (-self.alpha) / (2.0 * np.pi)
        N2 = float(self.N ** 2)
        return N2 * np.sqrt(np.maximum(S * self.dk ** 2, 0.0))

    def _turb_shape(self):
        """Horoshenkov W sampled on the grid, FFT'd, clipped at 0 (1.1 method)."""
        N, L = self.N, self.L
        dx = L / N
        ax = (np.arange(N) - N // 2) * dx
        X, Y = np.meshgrid(ax, ax, indexing="xy")
        r = np.sqrt(X * X + Y * Y)
        W = (np.exp(-r * r / (2.0 * self.sigma_w ** 2))
             * np.cos(2.0 * np.pi * r / self.L0))
        P = np.real(np.fft.fft2(np.fft.ifftshift(W))) * dx * dx
        S = np.clip(P, 0.0, None)
        S[0, 0] = 0.0
        N2 = float(N ** 2)
        return N2 * np.sqrt(np.maximum(S * self.dk ** 2, 0.0))

    # -- public ------------------------------------------------------------
    def eta_hat_at(self, t):
        t = float(t)
        tau = self.tau
        if np.isfinite(tau) and tau > 0:
            e = float(np.exp(-t / tau))
            s = float(np.sqrt(max(1.0 - np.exp(-2.0 * t / tau), 0.0)))
        else:
            e, s = 1.0, 0.0
        c_ii = np.exp(-1j * self.om_ii * t)
        c_iii = np.exp(-1j * self.om_iii * t)
        turb = (self.hat_turbA * e + self.hat_turbB * s) * c_iii
        return self.eta_i_hat + self.hat_ii0 * c_ii + turb

    def fields_at(self, t):
        eh = self.eta_hat_at(t)
        h = np.real(np.fft.ifft2(eh))
        gx = np.real(np.fft.ifft2(1j * self.kx * eh))
        gy = np.real(np.fft.ifft2(1j * self.ky * eh))
        z = np.zeros_like(h)
        return {"h": h.astype(np.float32), "gx": gx.astype(np.float32),
                "gy": gy.astype(np.float32), "Dx": z.astype(np.float32),
                "Dy": z.astype(np.float32), "S": float(self.S)}

    def measure(self, hat):
        """sigma and rms slope from ONE measurement path (ifft of the hat)."""
        h = np.real(np.fft.ifft2(hat))
        gx = np.real(np.fft.ifft2(1j * self.kx * hat))
        gy = np.real(np.fft.ifft2(1j * self.ky * hat))
        return {"sigma_m": float(np.sqrt(np.mean(h * h))),
                "rms_slope": float(np.sqrt(np.mean(gx * gx + gy * gy))),
                "rms_slope_x": float(np.sqrt(np.mean(gx * gx))),
                "rms_slope_y": float(np.sqrt(np.mean(gy * gy)))}

    def component_table(self, t=0.0):
        """Component-wise and total sigma/slope, all from the SAME path."""
        t = float(t)
        tau = self.tau
        if np.isfinite(tau) and tau > 0:
            e = float(np.exp(-t / tau))
            s = float(np.sqrt(max(1.0 - np.exp(-2.0 * t / tau), 0.0)))
        else:
            e, s = 1.0, 0.0
        c_ii = np.exp(-1j * self.om_ii * t)
        c_iii = np.exp(-1j * self.om_iii * t)
        h_i = self.eta_i_hat
        h_ii = self.hat_ii0 * c_ii
        h_iii = (self.hat_turbA * e + self.hat_turbB * s) * c_iii
        tot = h_i + h_ii + h_iii
        out = {"t_s": t,
               "i": self.measure(h_i), "ii": self.measure(h_ii),
               "iii": self.measure(h_iii), "total": self.measure(tot),
               "ii_ring": self.measure(self.hat_ring * c_ii),
               "ii_tail": self.measure(self.hat_tail * c_ii)}
        # consistency identity: total^2 = sum(comp^2) + 2*cross
        st = out["total"]["rms_slope"] ** 2
        sc = (out["i"]["rms_slope"] ** 2 + out["ii"]["rms_slope"] ** 2
              + out["iii"]["rms_slope"] ** 2)
        out["slope_cross_term"] = float(st - sc)
        out["slope_cross_fraction_of_total"] = float(
            (st - sc) / st if st > 0 else 0.0)
        out["note"] = ("every entry is measured from the same ifft-gradient "
                       "path, so total^2 = sum(component^2) + 2*cross holds "
                       "by construction")
        return out

    def spectrum_kx(self, hat=None):
        """1D ky-summed |eta_hat|^2 spectrum and its peak (sign-agnostic)."""
        if hat is None:
            hat = self.eta_hat_at(0.0)
        P = (np.abs(hat) ** 2) / (self.N ** 4)
        Sk = P.sum(axis=0)
        kxv = self.kx[0, :]
        m = np.abs(kxv) > 0
        j = int(np.argmax(Sk[m]))
        kpk = float(kxv[m][j])
        return {"kx_axis": [float(v) for v in kxv[m]],
                "S_kx": [float(v) for v in Sk[m]],
                "peak_kx": kpk, "peak_abs_kx": abs(kpk)}

    def spectrum_radial(self, hat=None):
        """Radial (annulus-summed) spectrum: the robust estimate of the ring
        radius.  Bins are one dk wide, so the peak lands on k0 within +-1 dk.
        """
        if hat is None:
            hat = self.eta_hat_at(0.0)
        P = (np.abs(hat) ** 2) / (self.N ** 4)
        kf = self.kappa.ravel()
        Pf = P.ravel()
        m = kf > 0
        kf, Pf = kf[m], Pf[m]
        nb = int(np.ceil(self.kmax / self.dk)) + 1
        idx = np.clip((kf / self.dk).astype(int), 0, nb - 1)
        Sk = np.bincount(idx, weights=Pf, minlength=nb)
        centres = (np.arange(nb) + 0.5) * self.dk
        j = int(np.argmax(Sk))
        return {"k_axis": [float(v) for v in centres],
                "S_k": [float(v) for v in Sk],
                "peak_k": float(centres[j])}

    def expectation_peak_kx(self, width_modes=0):
        """Peak of the 1D kx spectrum of the EXPECTED (phase-noise free)
        spectral amplitudes, i.e. the same quantity with |Z|^2 replaced by its
        mean.  Deterministic, so it is the noise-free counterpart of
        spectrum_kx().
        """
        hat = np.sqrt(np.abs(self.hat_ring) ** 2 + np.abs(self.hat_tail) ** 2
                      + np.abs(self.eta_i_hat) ** 2)
        return self.spectrum_kx(hat)["peak_abs_kx"]

    def describe(self):
        vtot = self.sigma_tot ** 2
        return {
            "N": self.N, "L_m": self.L, "H_m": self.H, "U_m_s": self.U,
            "U0_m_s": self.U0, "S0": self.S0, "D84_m": self.D84
            if hasattr(self, "D84") else self.ks,
            "seed": self.seed, "alpha": self.alpha,
            "derived": self.derived, "cmin_m_s": self.cmin,
            "cmin_k_rad_m": self.cmin_k,
            "k0_rad_m": self.k0,
            "lambda0_m": (2.0 * np.pi / self.k0) if self.k0 else None,
            "ring_fraction_phi": self.phi, "ring_fraction_detail": self.phi_info,
            "sigma_total_target_m": self.sigma_tot,
            "var_total_target": vtot,
            "var_i_physical": self.var_i,
            "var_i_fraction_of_total": (self.var_i / vtot if vtot > 0 else None),
            "var_ii": self.var_ii, "var_iii": self.var_iii,
            "var_ring": self.var_ring, "var_tail": self.var_tail,
            "sigma_i_m": float(np.sqrt(max(self.var_i, 0.0))),
            "sigma_ii_m": float(np.sqrt(max(self.var_ii, 0.0))),
            "sigma_iii_m": float(np.sqrt(max(self.var_iii, 0.0))),
            "sigma_ring_m": float(np.sqrt(max(self.var_ring, 0.0))),
            "sigma_tail_m": float(np.sqrt(max(self.var_tail, 0.0))),
            "L0": self.L0_info, "sigma_w_m": self.sigma_w,
            "tau_s": self.tau,
            "info_i": {k: v for k, v in self.info_i.items()
                       if not k.startswith("_")},
            "budget_rule": ("var(ii)=phi*(sigma_tot^2-var_i); "
                            "var(iii)=(1-phi)*(sigma_tot^2-var_i); "
                            "inside (ii): ring=phi*var(ii), tail=(1-phi)*var(ii)"),
        }


def validity_label(Fr, H, D84, H_fit_range=(0.0405, 0.1013)):
    outs = []
    if float(Fr) > 0.7:
        outs.append("Fr>0.7")
    if float(H) / float(D84) < 2.0:
        outs.append("H/D84<2")
    if not (H_fit_range[0] <= float(H) <= H_fit_range[1]):
        outs.append("H outside the 2016 fit range")
    return ("anchor out-of-domain (%s)" % ",".join(outs)) if outs else "in-domain"


def loop_adjusted_L(L_guess, Us, T_loop=8.0):
    n = max(1, int(round(T_loop * float(Us) / float(L_guess))))
    L = T_loop * float(Us) / float(n)
    return L, n, {"L_guess": float(L_guess), "L_adj": float(L),
                  "n_loops": int(n), "T_loop_mode": float(L / float(Us)),
                  "rule": "n = round(8 Us/L_guess), L = 8 Us/n"}
