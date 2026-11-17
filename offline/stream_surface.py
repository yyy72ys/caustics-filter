"""stream_surface -- BRIEF-V4-1-1 flowing-water surface generator (new file).

Design v2 (self-contained in BRIEF-V4-1-1.md):

Input: U (mean velocity m/s), H (mean depth m), S0 (slope -),
  bed spectrum (sigma_b rms m, ell_b correlation m; D84 mapping stated by
  caller), seed (int).

Derived: Fr = U/sqrt(gH), kappa = |k|, u* = sqrt(g H S0) (rectangular
  approximation; formula stated here and repeated in the report).

Component A (bed response, stationary):
  bed zeta = Gaussian random-phase field, spectrum S_bed defined by
  sigma_b and ell_b (isotropic Gaussian covariance):
    S_bed(k) = sigma_b^2 * ell_b^2 / (2*pi) * exp(-kappa^2 ell_b^2 / 2).
  eta_A_hat = T(k) * zeta_hat, with
    T(k) = [cosh(kH) - k*g*sinh(kH)/(kx^2 U^2)]^-1.
  Resonance is the gravity stationary condition (kx U)^2 = g k tanh(kH)
  (surface tension is NOT in T). kx = 0 rows: T = 0 (singularity avoidance).
  Resonant-band regularisation = breaking-limit rule: the resonant band
  (|denom| < RES_BAND_DENOM, defined below) contribution to surface slope,
  if its rms exceeds SLOPE_CAP = 0.30, is uniformly scaled to the cap.
  The scaled amount (dropped energy fraction %) is reported.

Component B: NOT added independently. The near-resonant band is diagnosed
  numerically only (band rms, band energy fraction) and reported separately.
  No double counting.

Component C (turbulence roughness, moving):
  Horoshenkov-type autocovariance W(rho) = exp(-rho^2/2σw^2) cos(2πρ/L0),
  reference σw = 0.22 m, L0 = 0.17 m. Mapping rule to scene dims (stated):
  σw, L0 are PHYSICAL metres used as-is (no rescaling) provided the patch
  resolves them (L/min(Nx,Ny)... at least 8 cells per L0 and patch >= 6σw);
  otherwise the run aborts with a clear error instead of silently rescaling.
  The spectrum is the discrete periodic correspondence of W: W sampled on
  the scene grid, FFT to power, clipped at 0, renormalised so that
  sum(S)*dkx*dky = sigma_c^2. sigma_c (rms of C) is NOT in the brief, so
  this module sets sigma_c = u*^2/g (velocity head) as a stated ASSUMPTION.
  Time evolution: each mode phase advances with ω = kx*Us + ω0(k),
  Us = surface velocity = US_FACTOR * U with US_FACTOR = 1.15 (log-profile
  assumption, stated), ω0(k) = sqrt((g k + σk^3/ρ) tanh(kH)) (gravity-
  capillary intrinsic frequency, σk = 0.0728 N/m, ρ = 1000 kg/m^3).
  Each mode complex amplitude additionally gets an OU update with
  correlation time τ = H/u* (stationary-variance preserving):
    A(t) = [A0 exp(-t/τ) + B sqrt(1-exp(-2t/τ))] exp(-iωt),
  A0, B independent draws with the same modal variance. τ is reported.

Loop consistency: advected component T_loop = Lx/Us must divide the 8 s
  loop an integer number of times. Lx (= Ly = L for the square det patch)
  is micro-adjusted: n = max(1, round(8 Us / L0_guess)), L = 8 Us / n.
  The adjustment is reported.

Validity label: Fr > 0.7 or H/D84 < 2 -> "anchor out-of-domain" label on
  numbers + report.

Conventions: periodic square patch Lx = Ly = L, Nx = Ny = N. Numpy FFT:
  eta(x) = Re(ifft2(eta_hat)), eta_hat = N^2 sqrt(S dk) Z (Z unit complex
  Gaussian, Hermitian-symmetrised), so that E|eta|^2 = sum S dk.
"""

import numpy as np

G = 9.81
SIGMA_K = 0.0728
RHO = 1000.0
US_FACTOR = 1.15
SLOPE_CAP = 0.30
RES_BAND_DENOM = 0.30

D84_SIGMA_DIVISOR = 4.0
D84_ELL_FACTOR = 1.0


def derive(U, H, S0):
    Fr = float(U) / float(np.sqrt(G * float(H)))
    ustar = float(np.sqrt(G * float(H) * float(S0)))
    tau = float(float(H) / ustar) if ustar > 0 else float("inf")
    Us = float(US_FACTOR * float(U))
    return {"Fr": Fr, "ustar": ustar, "tau": tau, "Us": Us,
            "formula_ustar": "u* = sqrt(g H S0) (rectangular approximation)",
            "formula_Us": "Us = 1.15 U (log-profile surface-factor assumption)"}


def bed_mapping_D84(D84):
    """D84 -> (sigma_b, ell_b) mapping (stated assumption)."""
    D84 = float(D84)
    return {"sigma_b": D84 / D84_SIGMA_DIVISOR,
            "ell_b": D84 * D84_ELL_FACTOR,
            "formula": "sigma_b = D84/4, ell_b = D84 (stated assumption)"}


def bed_mapping_sphere(d):
    """Mono-sphere bed d -> (sigma_b, ell_b) mapping (stated assumption)."""
    d = float(d)
    return {"sigma_b": d / 4.0, "ell_b": d,
            "formula": "sigma_b = d/4, ell_b = d (stated assumption)"}


def wavenumbers(N, L):
    f = np.fft.fftfreq(int(N), d=float(L) / float(int(N)))
    k1 = f * 2.0 * np.pi
    kx, ky = np.meshgrid(k1, k1, indexing="xy")
    kappa = np.sqrt(kx * kx + ky * ky)
    dkx = 2.0 * np.pi / float(L)
    return kx, ky, kappa, dkx


def hermitian_unit(N, rng):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
         / np.sqrt(2.0))
    # symmetrise: Zh(k) = (Z(k) + conj(Z(-k)))/sqrt2 -> Hermitian, CN(0,1)
    Zr = np.flip(np.flip(Z, axis=0), axis=1)
    Zh = (Z + np.conj(Zr)) / np.sqrt(2.0)
    return Zh


def transfer_T(kx, kappa, H, U):
    kx = np.asarray(kx, dtype=np.float64)
    kappa = np.asarray(kappa, dtype=np.float64)
    H = float(H)
    U = float(U)
    kh = kappa * H
    # stable: denom = cosh(kH) * (1 - k g tanh(kH)/(kx^2 U^2)); for kh > 50,
    # cosh overflows and T -> 0 regardless (denominator huge), so set T = 0.
    tanh = np.tanh(np.minimum(kh, 50.0))
    T = np.zeros_like(kh)
    nz = (np.abs(kx) > 0.0)
    safe = nz & (kh <= 50.0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cosh = np.cosh(np.minimum(kh, 50.0))
        denom = np.zeros_like(kh)
        f = (kappa * G * tanh / np.maximum(kx * kx * U * U, 1e-300))
        denom[safe] = cosh[safe] * (1.0 - f[safe])
        # kh > 50 (and kx != 0): denominator magnitude is astronomically
        # large -> T = 0 (already zero-initialised).
        T[safe] = 1.0 / denom[safe]
    T[~nz] = 0.0
    # reconstruct denom for diagnostics (inf where overflowed)
    with np.errstate(over="ignore"):
        denom_out = np.full_like(kh, np.inf)
        denom_out[safe] = denom[safe]
        denom_out[~nz] = np.inf
    return T, denom_out


def bed_spectrum(kappa, sigma_b, ell_b):
    kappa = np.asarray(kappa, dtype=np.float64)
    return (float(sigma_b) ** 2 * float(ell_b) ** 2 / (2.0 * np.pi)
            * np.exp(-0.5 * (kappa * float(ell_b)) ** 2))


def generate_zeta_hat(N, L, sigma_b, ell_b, rng):
    kx, ky, kappa, dk = wavenumbers(N, L)
    S = bed_spectrum(kappa, sigma_b, ell_b)
    N2 = int(N) * int(N)
    Zh = hermitian_unit(int(N), rng)
    amp = float(N2) * np.sqrt(np.maximum(S * dk * dk, 0.0))
    zeta_hat = amp * Zh
    # kx = 0 column of fftfreq grid: kx == 0 exactly for that column; the
    # zero response is enforced in T, zeta itself keeps kx=0 energy.
    return zeta_hat, {"S_bed": S, "kx": kx, "ky": ky, "kappa": kappa, "dk": dk}


def component_A(N, L, H, U, sigma_b, ell_b, rng):
    zeta_hat, grid = generate_zeta_hat(N, L, sigma_b, ell_b, rng)
    kx = grid["kx"]
    kappa = grid["kappa"]
    dk = grid["dk"]
    T, denom = transfer_T(kx, kappa, H, U)
    eta_hat = T * zeta_hat
    N2 = int(N) * int(N)
    # rms slope of resonant band before scaling
    band = np.abs(denom) < RES_BAND_DENOM
    band[~np.isfinite(denom)] = False
    # kx=0 row is T=0 by construction; exclude from band
    band[np.abs(kx) == 0.0] = False
    k2 = kappa * kappa
    tot_e = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
    band_e = float(np.sum(np.abs(eta_hat[band]) ** 2) / N2 ** 2)
    band_slope_var = float(np.sum(k2[band] * np.abs(eta_hat[band]) ** 2)
                           / N2 ** 2)
    band_rms = float(np.sqrt(max(band_slope_var, 0.0)))
    scale = 1.0
    if band_rms > SLOPE_CAP and band_rms > 0:
        scale = SLOPE_CAP / band_rms
        eta_hat[band] = eta_hat[band] * scale
    tot_e_after = float(np.sum(np.abs(eta_hat) ** 2) / N2 ** 2)
    dropped_total_pct = (100.0 * (tot_e - tot_e_after) / tot_e
                         if tot_e > 0 else 0.0)
    dropped_band_pct = (100.0 * (1.0 - scale ** 2) if band_e > 0 else 0.0)
    slope_var_after = float(np.sum(k2 * np.abs(eta_hat) ** 2) / N2 ** 2)
    info = {
        "T": T, "denom": denom, "zeta_hat": zeta_hat,
        "res_band_def": "|denom| < %.3g (denom = cosh(kH) - k g sinh(kH)/"
                        "(kx^2 U^2))" % RES_BAND_DENOM,
        "slope_cap": SLOPE_CAP,
        "band_rms_before": band_rms,
        "band_energy": band_e,
        "band_energy_fraction": (band_e / tot_e if tot_e > 0 else 0.0),
        "scale_applied": scale,
        "dropped_band_pct": float(dropped_band_pct),
        "dropped_total_pct": float(dropped_total_pct),
        "total_var_before": tot_e,
        "total_var_after": tot_e_after,
        "rms_slope_A": float(np.sqrt(max(slope_var_after, 0.0))),
        "sigma_eta_A": float(np.sqrt(max(tot_e_after, 0.0))),
    }
    return eta_hat, info


def horoshenkov_W_grid(N, L, sigma_w=0.22, L0=0.17):
    N = int(N)
    L = float(L)
    dx = L / N
    # resolution guard (mapping rule): abort rather than silently rescale
    if L0 / dx < 8.0:
        raise ValueError("L0 unresolved: L0/dx = %.2f < 8 (N=%d L=%.4g)" %
                         (L0 / dx, N, L))
    if L < 6.0 * float(sigma_w):
        raise ValueError("patch too small for W tail: L=%.4g < 6σw=%.4g" %
                         (L, 6.0 * float(sigma_w)))
    ax = (np.arange(N) - N // 2) * dx
    X, Y = np.meshgrid(ax, ax, indexing="xy")
    r = np.sqrt(X * X + Y * Y)
    W = np.exp(-r * r / (2.0 * float(sigma_w) ** 2)) * np.cos(
        2.0 * np.pi * r / float(L0))
    return W, dx


def component_C_spectrum(N, L, sigma_c, sigma_w=0.22, L0=0.17):
    Wg, dx = horoshenkov_W_grid(N, L, sigma_w, L0)
    # continuous-FT approx: P = FFT(W) dx dy, in fftfreq order (k=0 first).
    Praw = np.real(np.fft.fft2(np.fft.ifftshift(Wg))) * dx * dx
    # Praw is in fftfreq order already (k=0 first). Use as-is.
    Pclip = np.clip(Praw, 0.0, None)
    kx, ky, kappa, dk = wavenumbers(N, L)
    s = float(np.sum(Pclip) * dk * dk)
    if s <= 0:
        raise ValueError("Horoshenkov spectrum is empty after clipping")
    S = (float(sigma_c) ** 2) * Pclip / s
    return S, {"kx": kx, "ky": ky, "kappa": kappa, "dk": dk,
               "sigma_w": float(sigma_w), "L0": float(L0),
               "mapping": "σw,L0 are physical metres used as-is (no "
                          "rescaling); abort if L0/dx<8 or L<6σw"}


def omega0_field(kappa, H):
    k = np.asarray(kappa, dtype=np.float64)
    return np.sqrt(np.maximum((G * k + SIGMA_K * k ** 3 / RHO)
                              * np.tanh(k * float(H)), 0.0))


class StreamState:
    """Square periodic stream surface with fields_at(t) -> det-compatible."""

    def __init__(self, N=512, L=4.0, H=0.35, U=0.8, S0=0.002, D84=0.09,
                 sigma_b=None, ell_b=None, seed=1, sigma_w=0.22, L0=0.17):
        self.N = int(N)
        self.L = float(L)
        self.S = float(L)
        self.H = float(H)
        self.U = float(U)
        self.S0 = float(S0)
        self.seed = int(seed)
        d = derive(U, H, S0)
        self.Fr = d["Fr"]
        self.ustar = d["ustar"]
        self.tau = d["tau"]
        self.Us = d["Us"]
        self.D84 = float(D84)
        if sigma_b is None or ell_b is None:
            m = bed_mapping_D84(D84)
            sigma_b = m["sigma_b"]
            ell_b = m["ell_b"]
            self.bed_mapping = m
        else:
            self.bed_mapping = {
                "sigma_b": float(sigma_b), "ell_b": float(ell_b),
                "formula": "explicit sigma_b/ell_b (caller-stated)"}
        self.sigma_b = float(sigma_b)
        self.ell_b = float(ell_b)
        self.sigma_w = float(sigma_w)
        self.L0 = float(L0)
        # velocity-head assumption for C amplitude (stated)
        self.sigma_c = float(self.ustar ** 2 / G) if self.ustar > 0 else 0.0

        rng = np.random.default_rng(self.seed)
        self.etaA_hat, self.Ainfo = component_A(
            self.N, self.L, self.H, self.U, self.sigma_b, self.ell_b, rng)
        self.SC, self.Cgrid = component_C_spectrum(
            self.N, self.L, self.sigma_c, self.sigma_w, self.L0)
        kx = self.Cgrid["kx"]
        kappa = self.Cgrid["kappa"]
        dk = self.Cgrid["dk"]
        N2 = self.N * self.N
        rngA = np.random.default_rng(self.seed + 0xC0A1)
        rngB = np.random.default_rng(self.seed + 0xB0B1)
        self.A0 = (float(N2) * np.sqrt(np.maximum(self.SC * dk * dk, 0.0))
                   * hermitian_unit(self.N, rngA))
        self.B0 = (float(N2) * np.sqrt(np.maximum(self.SC * dk * dk, 0.0))
                   * hermitian_unit(self.N, rngB))
        self.om0 = omega0_field(kappa, self.H)
        self.om = kx * self.Us + self.om0

    def eta_hat_at(self, t):
        t = float(t)
        e = float(np.exp(-t / self.tau)) if np.isfinite(self.tau) else 1.0
        s = float(np.sqrt(max(1.0 - np.exp(-2.0 * t / self.tau), 0.0))) \
            if np.isfinite(self.tau) and self.tau > 0 else 0.0
        Ac = (self.A0 * e + self.B0 * s) * np.exp(-1j * self.om * t)
        return self.etaA_hat + Ac

    def fields_at(self, t):
        eh = self.eta_hat_at(t)
        h = np.real(np.fft.ifft2(eh)).astype(np.float64)
        kx = self.Cgrid["kx"]
        ky = self.Cgrid["ky"]
        gx = np.real(np.fft.ifft2(1j * kx * eh)).astype(np.float64)
        gy = np.real(np.fft.ifft2(1j * ky * eh)).astype(np.float64)
        z = np.zeros_like(h)
        return {"h": h.astype(np.float32), "gx": gx.astype(np.float32),
                "gy": gy.astype(np.float32), "Dx": z.astype(np.float32),
                "Dy": z.astype(np.float32), "S": float(self.S)}

    def realised_stats(self):
        f = self.fields_at(0.0)
        gx = np.asarray(f["gx"], dtype=np.float64)
        gy = np.asarray(f["gy"], dtype=np.float64)
        h = np.asarray(f["h"], dtype=np.float64)
        rms_slope = float(np.sqrt(np.mean(gx * gx + gy * gy)))
        N2 = self.N * self.N
        k2 = self.Cgrid["kappa"] ** 2
        rmsC0 = float(np.sqrt(np.sum(k2 * np.abs(self.A0) ** 2) / N2 ** 2))
        return {"rms_slope": rms_slope,
                "rms_slope_A": float(self.Ainfo["rms_slope_A"]),
                "rms_slope_C_t0": rmsC0,
                "sigma_eta_A": float(self.Ainfo["sigma_eta_A"]),
                "sigma_eta_total_t0": float(np.sqrt(np.mean(h * h)))}


def loop_adjusted_L(L_guess, Us, T_loop=8.0):
    n = max(1, int(round(T_loop * float(Us) / float(L_guess))))
    L = T_loop * float(Us) / float(n)
    return L, n, {"L_guess": float(L_guess), "L_adj": float(L),
                  "n_loops": int(n), "T_loop_mode": float(L / float(Us)),
                  "rule": "n = round(8 Us/L_guess), L = 8 Us/n"}


def validity_label(Fr, H, D84):
    outs = []
    if float(Fr) > 0.7:
        outs.append("Fr>0.7")
    if float(H) / float(D84) < 2.0:
        outs.append("H/D84<2")
    return ("anchor out-of-domain (%s)" % ",".join(outs)) if outs else "in-domain"
