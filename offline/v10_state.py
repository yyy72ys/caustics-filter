"""v10 state -- version 3.7 two-regime wave state (PIPELINE-V9 sections 4-6).

Wraps nothing: it is a new state in the v10 family.  It does NOT touch
v5..v9 model files; it imports the frozen dispersion / helpers from
v8_physics (the design says "existing") and the spectrum from v10_spectrum.

What it does
------------
1. AMPLITUDE (design section 4).  H_s comes from the fetch laws (JONSWAP
   family, depth-limited in the same Young-1988 family); ``waveHeightM`` is
   RETIRED as a source.  A recipe override is honoured and labelled
   ``hs_source = law|override``.  The breaking clamp min(0.142 lambda, 0.78 h)
   is a POST-CHECK: it scales the spectrum only if it binds, and records it.

2. PATCH (section 5).  The patch is the scene scale PATCH_BASE_M = 1.2 m
   (4 * 0.30 m, the v3.4 low-wind limit of env.py spec_L) and is NOT coupled
   to wind.  It is doubled in discrete steps only until it can hold the peak:
   S >= R_LAMBDA * lambda_p (R_LAMBDA = 8; the design's minimum is 4).  Every
   doubling is recorded, and the doubling is the convergence test of G6.

3. CASCADES (section 5).  The ONE spectrum S = S_g + S_h is realised on
   - a gravity grid (the scene grid: N points over S), and
   - a capillary grid (N points over S/2^j) when the scene grid cannot serve
     the capillary end (k_Nyq,scene < 1.1 k_cut_max).  The capillary tile is
     tiled 2^j times across the scene, so the composite field stays exactly
     scene-periodic.
   The split is a partition of unity in ln k: W_lo(k) + W_hi(k) = 1, centred
   on K_X = 100 rad/m (the centre of the design's 50-200 rad/m overlap band).
   No mode is carried twice, so the split cannot double count.
   Both grids quantise omega to the SAME 2 pi / 8 s grid with the same
   function.  Both draw phases from ONE RNG stream in a fixed order (the
   shared seed family), and the amplitudes of the surviving modes do not
   depend on k_cut, so the G6 cutoff sweep moves only the cutoff.

4. TIME (section 6).  The 8 s quantisation is untouched.  Its per-band error
   is measured (G7) on the realised mode set.

5. RAIN.  Out of scope for v3.7 (design section 8): R > 0 RAISES.  R = 0 is
   the only accepted value and the state is built with no rain branch at all,
   which is what makes gate G8(d) an equality rather than a comparison.

Units: k [rad/m], S [m^2 per rad/m], h/hs [m], omega [rad/s].
"""
import gc
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10

OMEGA0 = P8.OMEGA0
TAPER_NONE = None


def _logstep(k, k1, k2):
    """Raised-cosine partition in ln k: 0 below k1, 1 above k2."""
    k = np.maximum(np.asarray(k, dtype=np.float64), 1e-30)
    s = (np.log(k) - math.log(k1)) / (math.log(k2) - math.log(k1))
    s = np.clip(s, 0.0, 1.0)
    return s * s * (3.0 - 2.0 * s)


class _Cascade:
    """One spectral grid of the realisation."""

    def __init__(self, N, tile_m, tag):
        self.N = int(N)
        self.tile_m = float(tile_m)
        self.tag = tag
        self.dx = self.tile_m / self.N
        fx = np.fft.fftfreq(self.N, d=self.dx)
        kx1 = 2.0 * np.pi * fx
        kxg, kyg = np.meshgrid(kx1, kx1)
        self.kx = kxg.astype(np.float32)
        self.ky = kyg.astype(np.float32)
        k = np.sqrt(kxg.astype(np.float64) ** 2
                    + kyg.astype(np.float64) ** 2)
        del kxg, kyg
        self.k = k.astype(np.float32)
        self.kmin = 2.0 * math.pi / self.tile_m
        self.kmax_nyq = math.pi / self.dx
        self.dk = 2.0 * math.pi / self.tile_m
        self.H0 = None
        self.omega_q = None
        self.stats = None


class V10WaterState:
    """Version 3.7 environmental state: one spectrum, two regimes."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 fetch_m=None, height_m=None, k_cut=V10.KCUT_DEFAULT,
                 rain_rate_mm_h=0.0, patch_base_m=None, r_lambda=None,
                 gamma=V10.JONSWAP_GAMMA, theta_w=None):
        from state import derive_seed
        if float(rain_rate_mm_h) > 0.0:
            raise NotImplementedError(
                "v3.7 does not re-derive rain: rain must be 0 (PIPELINE-V9 "
                "section 8).  Refusing to build a state that would silently "
                "drop the rain band.")
        self.rain_rate = 0.0
        self.N = int(N)
        self.seed = int(seed) if seed is not None else int(
            derive_seed(recipe, preset_index))
        self.preset_index = int(preset_index)
        self.wind = float(recipe_env["windSpeedMS"])
        self.depth = float(recipe_env["depthM"])
        if fetch_m is None:
            raise RuntimeError("v10 state needs an explicit fetch (wave.fetchM)")
        self.fetch_m = float(fetch_m)
        wd = recipe_env.get("windDir", [1.0, 0.0])
        self.theta_w = (math.atan2(float(wd[1]), float(wd[0]))
                        if theta_w is None else float(theta_w))
        self.k_cut = float(k_cut)
        self.k_cut_max = max(V10.KCUT_SWEEP)
        base = float(patch_base_m if patch_base_m is not None
                     else V10.PATCH_BASE_M)
        rlam = float(r_lambda if r_lambda is not None else V10.R_LAMBDA)

        # ---- 1. spectrum assembly + amplitude ---------------------------
        self.ctx = V10.assemble(self.wind, self.fetch_m, self.depth,
                                hs_override=height_m, k_cut=self.k_cut,
                                gamma=gamma)
        self.laws = self.ctx["laws"]
        self.hs_source = self.ctx["hs_source"]

        # ---- 2. patch (scene scale, resolution-driven doublings) --------
        lam_p = float(self.laws["lambda_p_m"])
        pf = V10.patch_for(lam_p, base, rlam)
        self.S = float(pf["patch_m"])
        self.patch_doublings = int(pf["doublings"])
        self.patch_base_m = base
        self.r_lambda = rlam
        self.specL = self.S / 4.0
        self.dx = self.S / self.N
        self.kmin = 2.0 * math.pi / self.S
        self.kmax_nyq = math.pi / self.dx

        # ---- 3. cascade plan -------------------------------------------
        k_nyq_scene = self.kmax_nyq
        self.k_x_eff = float(V10.K_X)
        need_nyq = 1.1 * self.k_cut_max
        if k_nyq_scene >= need_nyq:
            j = 0
        else:
            j = int(math.ceil(math.log(need_nyq / k_nyq_scene, 2.0) - 1e-12))
            j = max(j, 1)
        # both grids must be able to REPRESENT their half of the partition:
        #   lo: k_Nyq,lo >= sqrt2 * K_X * R      hi: kmin,hi <= K_X/sqrt2 / R
        R_EDGE = 1.6
        kx_lo_max = k_nyq_scene / (math.sqrt(2.0) * R_EDGE)
        if j > 0:
            for _ in range(12):
                kmin_hi = 2.0 * math.pi / (self.S / float(2 ** j))
                kx_hi_min = math.sqrt(2.0) * R_EDGE * kmin_hi
                if kx_lo_max >= kx_hi_min:
                    break
                j += 1
            self.k_x_eff = float(min(V10.K_X, kx_lo_max))
            if self.k_x_eff < kx_hi_min:
                self.k_x_eff = float(0.5 * (kx_lo_max + kx_hi_min))
        self.k_x_eff = float(max(self.k_x_eff, 3.0 * self.kmin))
        self.cascade_j = j
        self.mode = "single" if j == 0 else "two_cascade"
        self.lo = _Cascade(self.N, self.S, "lo")
        self.hi = None
        self.W_hi_edges = [self.k_x_eff / math.sqrt(2.0),
                           self.k_x_eff * math.sqrt(2.0)]
        self._build()

    # ------------------------------------------------------------------
    # realisation
    # ------------------------------------------------------------------
    def _draw(self, rng, name):
        a = rng.standard_normal((self.N, self.N))
        b = rng.standard_normal((self.N, self.N))
        return a, b

    def _realise(self, casc, weight, rng, first):
        """Fill one cascade: H0, omega_q, stats.  weight is per-mode W."""
        kx = casc.kx.astype(np.float64)
        ky = casc.ky.astype(np.float64)
        k = casc.k
        E = V10.mode_density_2d(kx, ky, self.ctx, k_cut=None,
                                theta_w=self.theta_w)
        D = E * (casc.dk * casc.dk) * weight
        # --- termination of the MODE SET ---------------------------------
        # Modes are realised out to k_build so that the G6 cutoff sweep is a
        # pure mask test: the amplitude of every mode is independent of the
        # requested k_cut, which is applied afterwards.
        kvis = P8.viscous_cutoff(self.depth)
        k_build = min(self.k_cut_max, kvis, 0.95 * casc.kmax_nyq)
        D = np.where(k <= k_build, D, 0.0)
        D = np.where(np.isfinite(D), D, 0.0)
        del E
        a, b = self._draw(rng, casc.tag)
        scale = np.sqrt(2.0 * self.N ** 4 * np.maximum(D, 0.0))
        H0 = (scale * (a + 1j * b) / np.sqrt(2.0)).astype(np.complex64)
        del scale, a, b
        H0[0, 0] = 0.0
        # apply the REQUESTED cutoff (default 1000 rad/m)
        if self.k_cut < k_build:
            H0 = np.where(k <= self.k_cut, H0, np.complex64(0.0))
        casc.H0 = H0
        om = P8.dispersion_omega(np.maximum(k, 1e-12), self.depth)
        n = np.maximum(1.0, np.floor(om / OMEGA0 + 0.5))
        casc.omega_q = (n * OMEGA0).astype(np.float32)
        S_real = (np.abs(H0.astype(np.complex128)) ** 2) / (2.0 * self.N ** 4)
        sig = (k * k * D).sum()
        casc.stats = {
            "tag": casc.tag, "tile_m": casc.tile_m, "N": casc.N,
            "dx_mm": casc.dx * 1e3,
            "kmin": casc.kmin, "kmax_nyquist": casc.kmax_nyq,
            "k_build": float(k_build), "k_cut_applied": float(self.k_cut),
            "viscous_cutoff": float(kvis),
            "m0_target_precut": float(D.sum()),
            "sigma2_target_precut": float(sig),
            "m0_realised": float(S_real.sum()),
            "sigma2_realised": float((k * k * S_real).sum()),
            "n_modes_active": int((H0 != 0).sum()),
            "weight": ("W_hi" if casc.tag == "hi" else "W_lo"),
            "first_in_stream": bool(first),
        }
        # log-k binned cumulative m0 / sigma^2 versus k (for the G6 sweep)
        active = D > 0.0
        km = k[active]
        nb = 4000
        lo = max(float(km.min()) if km.size else 1e-6, 1e-9)
        hi = max(float(km.max()) if km.size else 1.0, lo * 1.0001)
        edges = np.logspace(math.log10(lo), math.log10(hi), nb + 1)
        idx = np.clip(np.digitize(km, edges) - 1, 0, nb - 1)
        b_m0 = np.bincount(idx, weights=D[active], minlength=nb)
        b_s2 = np.bincount(idx, weights=km * km * D[active], minlength=nb)
        b_m0r = np.bincount(idx, weights=S_real[active], minlength=nb)
        b_s2r = np.bincount(idx, weights=km * km * S_real[active],
                            minlength=nb)
        casc.stats["_cum"] = {
            "edges": edges,
            "bins_m0": b_m0, "bins_sigma2": b_s2,
            "m0": np.cumsum(b_m0), "sigma2": np.cumsum(b_s2),
            "m0_realised": np.cumsum(b_m0r),
            "sigma2_realised": np.cumsum(b_s2r),
        }
        del S_real, D, active, km
        return casc

    def _build(self):
        rng = np.random.default_rng(self.seed)
        if self.cascade_j == 0:
            weight_lo = np.ones((self.N, self.N), dtype=np.float64)
            self._realise(self.lo, weight_lo, rng, first=True)
            del weight_lo
        else:
            k2d = self.lo.k.astype(np.float64)
            w_hi = _logstep(k2d, self.W_hi_edges[0], self.W_hi_edges[1])
            w_lo = 1.0 - w_hi
            self._realise(self.lo, w_lo, rng, first=True)
            del w_lo, k2d
            self.hi = _Cascade(self.N, self.S / float(2 ** self.cascade_j),
                               "hi")
            k2d = self.hi.k.astype(np.float64)
            w_hi2 = _logstep(k2d, self.W_hi_edges[0], self.W_hi_edges[1])
            self._realise(self.hi, w_hi2, rng, first=False)
            del w_hi2, k2d
            del w_hi
        gc.collect()
        self._finish()

    def _finish(self):
        """Spectral readouts and the post-checks (no field arrays)."""
        lo = self.lo.stats
        hi = self.hi.stats if self.hi is not None else None
        m0 = lo["m0_realised"] + (hi["m0_realised"] if hi else 0.0)
        sig = lo["sigma2_realised"] + (hi["sigma2_realised"] if hi else 0.0)
        self.spectrum_after = {
            "m0_m2": float(m0), "Hs_m": float(4.0 * math.sqrt(max(m0, 0.0))),
            "sigma2_rad2": float(sig),
            "rms_slope": float(math.sqrt(max(sig, 0.0))),
            "note": ("spectral readout of the realised mode set (sum of the "
                     "drawn modal variances), before any spatial sampling"),
        }
        # breaking post-check
        lam_p = float(self.laws["lambda_p_m"])
        hmax = min(P8.BREAK_STEEPNESS * lam_p, P8.BREAK_DEPTH_FRAC * self.depth)
        hs = self.spectrum_after["Hs_m"]
        applied = bool(hs > hmax > 0.0)
        if applied:
            f = float(hmax / hs)
            self.lo.H0 = (self.lo.H0 * f).astype(np.complex64)
            if self.hi is not None:
                self.hi.H0 = (self.hi.H0 * f).astype(np.complex64)
            self.breaking_scale = f
            self.spectrum_after["m0_m2"] *= f * f
            self.spectrum_after["Hs_m"] *= f
            self.spectrum_after["sigma2_rad2"] *= f * f
            self.spectrum_after["rms_slope"] *= f
        else:
            self.breaking_scale = 1.0
        self.breaking = {
            "rule": "H_max = min(0.142 lambda_p, 0.78 h) (post-check)",
            "lambda_p_m": lam_p, "depth_m": self.depth,
            "H_stokes_0p142_m": float(P8.BREAK_STEEPNESS * lam_p),
            "H_stokes_1over7_m": float(P8.STOKES_STEEPNESS * lam_p),
            "H_depth_0p78h_m": float(P8.BREAK_DEPTH_FRAC * self.depth),
            "H_max_m": float(hmax),
            "Hs_before_m": float(hs),
            "applied": applied,
            "factor": float(self.breaking_scale),
            "Hs_after_m": float(hs * self.breaking_scale),
            "binding": ("stokes" if P8.BREAK_STEEPNESS * lam_p
                        <= P8.BREAK_DEPTH_FRAC * self.depth else "depth"),
        }
        # direction node + anisotropy
        self.direction = self._direction_node()
        # quantisation (G7)
        self.quantisation = self._quantisation_node()

    # ------------------------------------------------------------------
    # nodes
    # ------------------------------------------------------------------
    def _direction_node(self):
        kc = V10._k_grid(self.k_cut)
        an_c = V10.anisotropy_prediction(self.ctx, kc, self.kmin, self.k_cut)
        d = V10.eckv_delta(kc, self.depth, self.ctx["ustar"]["ustar_ms"],
                           self.ctx["k_p_rad_per_m"])
        return {
            "model": ("ECKV unified azimuthal function Phi = (1/2pi)(1 + "
                      "Delta cos 2(theta-theta_w)), Delta from Elfouhaily "
                      "et al. (1997) eq. 66-68"),
            "theta_w_deg": float(math.degrees(self.theta_w)),
            "a0": V10.ECKV_A0, "a_p": V10.ECKV_AP,
            "a_m": float(V10.ECKV_AM_K * self.ctx["ustar"]["ustar_ms"]
                         / V10.C_M),
            "delta_effective": float(an_c["delta_effective"]),
            "delta_min_over_band": float(np.min(d)),
            "delta_max": float(np.max(d)),
            "sigma2_along_predicted": float(an_c["sigma2_along"]),
            "sigma2_cross_predicted": float(an_c["sigma2_cross"]),
            "ratio_predicted": float(an_c["ratio_along_over_cross"]),
            "note": ("analytic on the assembled continuum (k ~ 1e-3 to "
                     "k_cut).  The realised field's own gx/gy ratio is "
                     "reported separately in state_node.realised and is the "
                     "grid-side check of the same number."),
        }

    def _quantisation_node(self):
        """Per-band max |domega|/omega of the 8 s loop (gate G7)."""
        bands = [("gt_1m", 0.0, 2 * math.pi),
                 ("0p1_1m", 2 * math.pi, 20 * math.pi),
                 ("1_10cm", 20 * math.pi, 200 * math.pi),
                 ("lt_1cm", 200 * math.pi, float("inf"))]
        acc = {b[0]: {"max_rel": 0.0, "n": 0} for b in bands}
        allmax = 0.0
        worst_k = None
        for casc in (self.lo, self.hi):
            if casc is None:
                continue
            k = casc.k.astype(np.float64)
            active = casc.H0 != 0
            om = P8.dispersion_omega(np.maximum(k, 1e-12), self.depth)
            rel = np.abs(casc.omega_q.astype(np.float64) - om) \
                / np.maximum(om, 1e-30)
            rel = np.where(active, rel, 0.0)
            for name, lo_k, hi_k in bands:
                m = active & (k >= lo_k) & (k < hi_k)
                if m.any():
                    r = float(rel[m].max())
                    acc[name]["max_rel"] = max(acc[name]["max_rel"], r)
                    acc[name]["n"] += int(m.sum())
            rmax = float(rel[active].max()) if active.any() else 0.0
            if rmax > allmax:
                i = int(np.argmax(np.where(active, rel, 0.0)))
                allmax = rmax
                worst_k = float(k.ravel()[i])
        return {
            "loop_quantum_rad_s": float(OMEGA0),
            "per_band_max_rel_error": acc,
            "max_rel_error": float(allmax),
            "k_at_max_rel_error_rad_per_m": worst_k,
            "note": ("omega is quantised to integer multiples of 2pi/8 s so "
                     "the 8 s loop closes exactly.  At t = 0 the field does "
                     "not depend on omega at all (exp(i*omega*0) = 1), so this "
                     "table only bounds the time-evolution error; the worst "
                     "mode is the longest one, as the design states"),
        }

    def sigma2_at_kcut(self, k_cut):
        """Spectral m0 / sigma^2 with the mask moved to k_cut (G6)."""
        m0 = 0.0
        sig = 0.0
        m0r = 0.0
        sigr = 0.0
        for casc in (self.lo, self.hi):
            if casc is None:
                continue
            c = casc.stats["_cum"]
            i = int(np.searchsorted(c["edges"], float(k_cut),
                                     side="right")) - 1
            i = min(i, len(c["m0"]) - 1)
            if i >= 0:
                m0 += float(c["m0"][i])
                sig += float(c["sigma2"][i])
                m0r += float(c["m0_realised"][i])
                sigr += float(c["sigma2_realised"][i])
        return {"k_cut": float(k_cut), "m0_target": m0, "sigma2_target": sig,
                "m0_realised": m0r, "sigma2_realised": sigr,
                "Hs_realised": float(4.0 * math.sqrt(max(m0r, 0.0)))}

    # ------------------------------------------------------------------
    # fields
    # ------------------------------------------------------------------
    def _cascade_fields(self, casc, t):
        wq = casc.omega_q
        if t == 0.0:
            H = casc.H0
        else:
            phase = np.exp(1j * wq.astype(np.float64) * float(t))
            H = (casc.H0 * phase.astype(np.complex64)).astype(np.complex64)
            del phase
        Hc = H.astype(np.complex128)
        h = np.real(np.fft.ifft2(Hc)).astype(np.float32)
        kxc = casc.kx.astype(np.complex128)
        kyc = casc.ky.astype(np.complex128)
        gx = np.real(np.fft.ifft2(1j * kxc * Hc)).astype(np.float32)
        out = {"h": h, "gx": gx}
        gy = np.real(np.fft.ifft2(1j * kyc * Hc)).astype(np.float32)
        out["gy"] = gy
        kmag = casc.k.astype(np.float64)
        kinv = np.zeros_like(kmag)
        nz = kmag > 0
        kinv[nz] = 1.0 / kmag[nz]
        lam = np.float32(0.85)
        Dx = np.real(np.fft.ifft2(-1j * kxc * kinv * lam * Hc)).astype(np.float32)
        Dy = np.real(np.fft.ifft2(-1j * kyc * kinv * lam * Hc)).astype(np.float32)
        out["Dx"] = Dx
        out["Dy"] = Dy
        del Hc
        return out

    def _fold_hi(self, f):
        """Block-average the capillary tile and tile it across the scene."""
        j = self.cascade_j
        if j == 0:
            return f
        b = 2 ** j
        n = self.N // b
        out = {}
        for key, arr in f.items():
            a = arr.astype(np.float64).reshape(n, b, n, b).mean(axis=(1, 3))
            out[key] = np.tile(a.astype(np.float32), (b, b))
            del a
        return out

    def fields_at(self, t):
        """Scene fields at time t: sum of the cascades (scene-periodic)."""
        lo = self._cascade_fields(self.lo, t)
        if self.hi is not None:
            hi = self._fold_hi(self._cascade_fields(self.hi, t))
            for key in lo:
                lo[key] = (lo[key] + hi[key]).astype(np.float32)
            del hi
            gc.collect()
        lo["S"] = float(self.S)
        return lo

    # ------------------------------------------------------------------
    # readouts
    # ------------------------------------------------------------------
    def frac_above(self, ks=(100.0, 370.0)):
        """Fraction of sigma^2 above k, from the mode sums (G2 validation)."""
        out = {}
        tot_r = self.spectrum_after["sigma2_rad2"]
        tot_t = self.sigma2_at_kcut(self.k_cut)["sigma2_target"]
        for kx in ks:
            hi_t = 0.0
            hi_r = 0.0
            for casc in (self.lo, self.hi):
                if casc is None:
                    continue
                c = casc.stats["_cum"]
                i = int(np.searchsorted(c["edges"], float(kx),
                                        side="right")) - 1
                i = min(i, len(c["sigma2"]) - 1)
                if i >= 0:
                    hi_t += float(c["sigma2"][-1]) - float(c["sigma2"][i])
                    hi_r += (float(c["sigma2_realised"][-1])
                             - float(c["sigma2_realised"][i]))
            out["above_%g" % kx] = {
                "target": (float(hi_t / tot_t) if tot_t > 0 else None),
                "realised": (float(hi_r / tot_r) if tot_r > 0 else None),
                "sigma2_above": float(hi_r),
            }
        return out

    def peak_readout(self):
        """Dominant wavelength: peak of the radially integrated m0 density."""
        best = None
        for casc in (self.lo, self.hi):
            if casc is None:
                continue
            c = casc.stats["_cum"]
            edges = c["edges"]
            dk = np.diff(edges)
            dens = c["bins_m0"] / np.maximum(dk, 1e-30)
            i = int(np.argmax(dens))
            if best is None or dens[i] > best["density"]:
                best = {"density": float(dens[i]),
                        "k_peak_rad_per_m": float(
                            math.sqrt(edges[i] * edges[i + 1])),
                        "cascade": casc.tag}
        kp = best["k_peak_rad_per_m"]
        return {"k_peak_rad_per_m": kp,
                "lambda_peak_m": float(2.0 * math.pi / kp),
                "variance_density_peak": best["density"],
                "peak_is_at_kmin": bool(kp <= self.kmin * 1.35),
                "carried_by_cascade": best["cascade"],
                "kmin_rad_per_m": float(self.kmin),
                "kmax_rad_per_m": float(self.kmax_nyq)}

    def realised_stats(self, fields):
        h = np.asarray(fields["h"], dtype=np.float64)
        gx = np.asarray(fields["gx"], dtype=np.float64)
        gy = np.asarray(fields["gy"], dtype=np.float64)
        m0 = float(np.mean(h * h))
        sig = float(np.mean(gx * gx + gy * gy))
        along = float(np.mean(gx * gx))
        cross = float(np.mean(gy * gy))
        return {"m0_m2": m0, "Hs_m": float(4.0 * math.sqrt(max(m0, 0.0))),
                "sigma2_rad2": sig, "rms_slope": float(math.sqrt(max(sig, 0.0))),
                "sigma2_along_rad2": along, "sigma2_cross_rad2": cross,
                "ratio_along_over_cross": (float(along / cross)
                                           if cross > 0 else None)}

    def band_shares(self):
        """Realised m0 / sigma^2 shares in the four bands (mode sums)."""
        bands = [("gt_1m", 0.0, 2 * math.pi),
                 ("0p1_1m", 2 * math.pi, 20 * math.pi),
                 ("1_10cm", 20 * math.pi, 200 * math.pi),
                 ("lt_1cm", 200 * math.pi, float("inf"))]
        out = {}
        tot_m0 = 0.0
        tot_s2 = 0.0
        per = {b[0]: {"m0": 0.0, "sigma2": 0.0} for b in bands}
        for casc in (self.lo, self.hi):
            if casc is None:
                continue
            k = casc.k.astype(np.float64)
            S = (np.abs(casc.H0.astype(np.complex128)) ** 2) / (2.0 * self.N ** 4)
            S = S * (self.breaking_scale ** 2)
            for name, lo_k, hi_k in bands:
                m = (k >= lo_k) & (k < hi_k)
                per[name]["m0"] += float(S[m].sum())
                per[name]["sigma2"] += float((k[m] ** 2 * S[m]).sum())
            tot_m0 += float(S.sum())
            tot_s2 += float((k * k * S).sum())
        for name, _, _ in bands:
            out[name] = {
                "m0": per[name]["m0"], "sigma2": per[name]["sigma2"],
                "m0_share": (per[name]["m0"] / tot_m0 if tot_m0 > 0 else None),
                "sigma2_share": (per[name]["sigma2"] / tot_s2
                                 if tot_s2 > 0 else None),
            }
        out["_total"] = {"m0": tot_m0, "sigma2": tot_s2}
        return out

    def state_node(self, fields=None):
        lo = self.lo.stats
        hi = self.hi.stats if self.hi is not None else None
        node = {
            "version": "3.7",
            "class": "V10WaterState",
            "mode": self.mode,
            "seed": self.seed,
            "patch": {
                "S_m": float(self.S), "base_m": float(self.patch_base_m),
                "doublings": int(self.patch_doublings),
                "r_lambda": float(self.r_lambda),
                "lambda_p_m": float(self.laws["lambda_p_m"]),
                "S_over_lambda_p": float(self.S / self.laws["lambda_p_m"]),
                "wind_independent_base": True,
                "note": ("the base is the scene scale 1.2 m and is not "
                         "coupled to wind; the doubling is a RESOLUTION "
                         "requirement (S >= r_lambda lambda_p) and is logged "
                         "per cell"),
            },
            "fetch": {"fetch_m": self.fetch_m,
                      "Ft": float(self.laws["Ft"]),
                      "f_p_Hz": float(self.laws["f_p_Hz"]),
                      "alpha": float(self.laws["alpha"]),
                      "t_min_s": float(self.laws["t_min_s"]),
                      "in_fitted_range_Ft_1e3_1e5":
                          bool(self.laws["in_fitted_range_Ft_1e3_1e5"])},
            "amplitude": {
                "Hs_law_m": float(self.laws["Hs_fetch_law_m"]),
                "Hs_depth_limited_m": float(self.laws["Hs_depth_limited_m"]),
                "depth_correction_rel":
                    float(self.laws["depth_correction_rel"]),
                "depth_correction_binds":
                    bool(self.laws["depth_correction_binds"]),
                "kp_times_h": float(self.laws["kp_times_h"]),
                "Hs_used_m": float(self.ctx["Hs_m"]),
                "hs_source": self.hs_source,
                "waveHeightM_used": False,
            },
            "spectrum": {
                "convention": "S(k), m0=int S dk, sigma2=int k^2 S dk",
                "S_g": "JONSWAP(omega; Hs, f_p) * domega/dk",
                "S_h": "Elfouhaily B_h(k; u*) / k^3",
                "additive_no_splice": True,
                "alpha_m": float(self.ctx["alpha_m_value"]),
                "alpha_m_raw": float(self.ctx["alpha_m"]["alpha_m_raw"]),
                "eckv_domain": self.ctx["eckv_domain"],
                "ustar_ms": float(self.ctx["ustar"]["ustar_ms"]),
                "u_over_cm": float(self.ctx["alpha_m"]["u_over_cm"]),
                "gamma": float(self.ctx["gamma"]),
                "k_cut_rad_per_m": float(self.k_cut),
                "k_cut_max_used": float(self.k_cut_max),
            },
            "cascades": {
                "n_cascades": 1 if self.hi is None else 2,
                "split_k_rad_per_m": float(self.k_x_eff),
                "split_source": ("K_X = 100 rad/m (centre of the 50-200 "
                                 "overlap band), moved only if a grid cannot "
                                 "represent its half"),
                "partition": ("W_lo + W_hi = 1, raised cosine in ln k between "
                              "%.3f and %.3f rad/m"
                              % (self.W_hi_edges[0], self.W_hi_edges[1])),
                "same_omega_grid": True,
                "shared_seed_family": ("one rng stream, gravity grid first, "
                                       "capillary grid second"),
                "lo": {k: v for k, v in lo.items() if k != "_cum"},
                "hi": (None if hi is None else
                       {k: v for k, v in hi.items() if k != "_cum"}),
                "tiled_capillary": bool(self.cascade_j > 0),
                "capillary_tiles_per_scene": int(2 ** self.cascade_j),
            },
            "breaking": self.breaking,
            "direction": self.direction,
            "quantisation": self.quantisation,
            "spectrum_after": self.spectrum_after,
            "bands": self.band_shares(),
        }
        if fields is not None:
            st = self.realised_stats(fields)
            node["realised"] = st
            node["peak"] = self.peak_readout()
            spec_s2 = self.spectrum_after["sigma2_rad2"]
            node["spectral_vs_spatial"] = {
                "Hs_spectral": self.spectrum_after["Hs_m"],
                "Hs_spatial": st["Hs_m"],
                "sigma2_spectral": spec_s2,
                "sigma2_spatial": st["sigma2_rad2"],
                "sigma2_resolved_fraction": (float(st["sigma2_rad2"] / spec_s2)
                                             if spec_s2 > 0 else None),
                "scene_nyquist_k": float(self.kmax_nyq),
                "note": ("the two readouts differ by (a) the few-mode "
                         "realisation scatter of this patch and (b) the scene "
                         "grid's own Nyquist, which at a large patch cannot "
                         "carry the capillary band into the FIELD.  The "
                         "spectral number is the physical surface statistic; "
                         "the spatial one is what the transport can see.  Both "
                         "are reported and the fraction is stated."),
            }
        return node
