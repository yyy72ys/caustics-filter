"""v6 water state: rain ring-wave band on top of the v5 calm state.

Wraps (does not modify) state.py / v5_state.py.  PIPELINE-V6 chunk 2 adds
exactly one new scene input, rain rate R [mm/h]; wind stays where the recipe
puts it (the verification scene is wind = 0, i.e. the Cox-Munk sigma^2 floor).

What is added to the spectral state (PIPELINE-V6 sections 1-4):

  1. An ISOTROPIC ring-wave band (report section 3.5 "spectral" method --
     no procedural ring texture, no per-drop particles), centred on
     k_ring(R) which is derived from the mass-median drop diameter D0(R)
     (bigger drops -> longer rings; report section 3.4).  The band's total
     variance is the physical energy budget

         sigma_ring^2 = eta * F(R) * Ebar(R) * tau / (rho g + sigma_s k^2)

     with F the impact flux (report 3.3), Ebar the flux-weighted mean impact
     kinetic energy (report 3.2), tau the ~150 ms impact-process lifetime
     (report 3.4) and eta the ONE calibrated constant of the chunk.

  2. Time evolution with the gravity-capillary dispersion
     omega^2 = g k + (sigma/rho) k^3 (report sections 3.4 / 6.1), quantised to
     multiples of 2*pi/8 s exactly as the existing state does, so the 8 s
     loop still closes bit-exactly.

  3. The second rain effect of report 3.4 / PIPELINE-V6 section 4: the
     sub-band (gravity) part of the spectrum is damped by a factor close to
     1 whose form is documented in v6_physics.gravity_band_damping and whose
     magnitude is a first-order rate estimate (no tuned constant).  In the
     windless scene it is almost a no-op, as designed.

R = 0 is a hard early return: the object is then bit-identical to a
v5_state.V5WaterState instance (same H0, same omega quantisation, no band,
no damping) -- that is the backward-compatibility gate.
"""
import numpy as np

from v5_state import V5WaterState
from v6_physics import (RING_ETA, RING_K_BAND, RING_SIGMA_K_FRAC,
                        gravity_capillary_omega, gravity_band_damping,
                        rain_state_summary, ring_band_variance,
                        ring_sigma_k, ring_wavenumber)

G = 9.81
OMEGA0 = 2.0 * np.pi / 8.0
RAIN_RNG_TAG = 0x6A17      # deterministic stream tag for the band phases


class V6WaterState(V5WaterState):
    """V5 calm state + a physically-budgeted isotropic ring-wave band."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 height_mode="wind", height_m=None, rain_rate_mm_h=0.0,
                 calm_calibrate=True, calm_target_sigma2=None,
                 eta=RING_ETA):
        self.rain_rate = max(0.0, float(rain_rate_mm_h))
        self.rain_eta = float(eta)
        self.rain = None
        self.rain_damping = None
        self.rain_band = None
        super().__init__(recipe, recipe_env, N=N, seed=seed,
                         preset_index=preset_index, height_mode=height_mode,
                         height_m=height_m, calm_calibrate=calm_calibrate,
                         calm_target_sigma2=calm_target_sigma2)
        self.calm_floor = self._spectrum_stats()
        if self.rain_rate > 0.0:
            self._apply_rain()

    # ------------------------------------------------------------------
    # dispersion
    # ------------------------------------------------------------------
    def _omega_q(self):
        """Quantised angular frequency.

        Without rain this is exactly the parent (deep-water sqrt(gk)) so the
        backward-compatibility gate is bit-exact.  With rain the full
        gravity-capillary dispersion is used, because the ring band sits on
        the crossover (report 3.4) where the capillary term carries as much
        as gravity; at the low-k wind band the capillary term is negligible.
        """
        if self.rain_rate <= 0.0:
            return super()._omega_q()
        k = np.maximum(self.k.astype(np.float64), 0.0)
        raw = gravity_capillary_omega(k)
        n = np.maximum(1.0, np.floor(raw / OMEGA0 + 0.5))
        return (n * OMEGA0).astype(np.float64)

    # ------------------------------------------------------------------
    # rain band
    # ------------------------------------------------------------------
    def _spectrum_density_of(self, H):
        return (np.abs(H.astype(np.complex128)) ** 2) / (2.0 * self.N ** 4)

    def _apply_rain(self):
        k = self.k.astype(np.float64)
        N = self.N

        # ---- (3) long-wave damping of the sub-band spectrum -------------
        S_now = self._spectrum_density_of(self.H0)
        low = k < RING_K_BAND[0]
        m0_low = float(S_now[low].sum())
        sig2_low = float((k[low] ** 2 * S_now[low]).sum())
        damp = gravity_band_damping(self.rain_rate, m0_low)
        w_low = 1.0 / (1.0 + (k / float(RING_K_BAND[0])) ** 4)
        gain_low = 1.0 / (1.0 + damp["gamma_per_s"] * damp["tau_s"] * w_low)
        self.H0 = (self.H0 * np.sqrt(gain_low)).astype(np.complex64)
        S_after = self._spectrum_density_of(self.H0)
        damp.update({
            "m0_lowk_before_m2": m0_low,
            "m0_lowk_after_m2": float(S_after[low].sum()),
            "sigma2_lowk_before_rad2": sig2_low,
            "sigma2_lowk_after_rad2": float((k[low] ** 2 * S_after[low]).sum()),
            "applied": True,
        })
        self.rain_damping = damp

        # ---- (1) isotropic ring band ------------------------------------
        rk = ring_wavenumber(self.rain_rate)
        k_ring = rk["k_ring_rad_per_m"]
        sk = ring_sigma_k(k_ring)
        var_target = ring_band_variance(self.rain_rate, eta=self.rain_eta)
        W = np.exp(-0.5 * ((k - k_ring) / sk) ** 2)
        tot_W = float(W.sum())
        if tot_W <= 0.0:
            raise RuntimeError("degenerate ring band weight")
        # E|H|^2 = A * W  ->  variance = sum(A W)/(2 N^4) = var_target
        A = float(var_target["ring_variance_m2"]) * (2.0 * N ** 4) / tot_W
        rng = np.random.default_rng(self.seed ^ RAIN_RNG_TAG)
        ga = rng.standard_normal((N, N))
        gb = rng.standard_normal((N, N))
        G = (ga + 1j * gb) / np.sqrt(2.0)
        H_ring = (G * np.sqrt(A * W)).astype(np.complex64)
        del ga, gb, G, W
        self.H0 = (self.H0 + H_ring).astype(np.complex64)

        S_ring = self._spectrum_density_of(H_ring)
        band_var = float(S_ring.sum())
        band_slope = float((k * k * S_ring).sum())
        in_band = (k >= RING_K_BAND[0]) & (k <= RING_K_BAND[1])
        frac_in = float(S_ring[in_band].sum() / band_var) if band_var > 0 else 0.0
        k_eff = float((k * S_ring).sum() / band_var) if band_var > 0 else 0.0
        k_eff_in = float((k[in_band] * S_ring[in_band]).sum()
                         / S_ring[in_band].sum()) if frac_in > 0 else 0.0
        # expected number of spectral modes carrying the band
        lam_band = float((S_ring > 1e-12 * max(S_ring.max(), 1e-30)).sum())
        self.rain_band = {
            "k_ring_nominal_rad_per_m": k_ring,
            "k_mean_effective_rad_per_m": k_eff,
            "k_mean_effective_in_band_rad_per_m": k_eff_in,
            "lambda_ring_m": rk["lambda_ring_m"],
            "sigma_k_rad_per_m": sk,
            "sigma_k_frac_of_band": RING_SIGMA_K_FRAC,
            "k_band_rad_per_m": list(RING_K_BAND),
            "variance_target_m2": var_target["ring_variance_m2"],
            "variance_realised_m2": band_var,
            "variance_realised_over_target": (band_var
                                              / var_target["ring_variance_m2"]),
            "slope_variance_realised_rad2": band_slope,
            "rms_slope_realised_rad": float(np.sqrt(band_slope)),
            "rms_height_realised_mm": float(np.sqrt(band_var) * 1e3),
            "variance_fraction_inside_band": frac_in,
            "modes_carrying_band": int(lam_band),
            "budget": var_target,
            "isotropic": True,
            "note": ("rings carry no direction: the bump is direction-uniform, "
                     "in contrast to the wind spectrum's cos^2 lobe "
                     "(PIPELINE-V6 section 2).  The Gaussian in k is not "
                     "hard-clipped: at the extreme rates the placement sits "
                     "on the band edge, so the tail crosses 310/630 "
                     "(containment fraction recorded)."),
        }
        self.rain = rain_state_summary(self.rain_rate, m0_low, eta=self.rain_eta)
        # the v5 attribute was captured before the band existed: refresh it so
        # every downstream readout sees the post-rain spectrum
        self.spectrum_after = self._spectrum_stats()
        return self.rain_band

    # ------------------------------------------------------------------
    # readouts
    # ------------------------------------------------------------------
    def rain_node(self, fields=None):
        """The rain node as reported in numbers-v6.json."""
        if self.rain_rate <= 0.0:
            return {"name": "rain.rate_mm_h", "value": 0.0, "active": False,
                    "note": "R=0 -> v5 behaviour exactly (no band, no damping)"}
        out = {"name": "rain.rate_mm_h", "value": float(self.rain_rate),
               "units": "mm/h", "active": True,
               "summary": self.rain, "band": self.rain_band,
               "damping": self.rain_damping}
        if fields is not None:
            st = self.realised_stats(fields)
            spec = self._spectrum_stats()
            out["realised"] = st
            out["sigma2_breakdown"] = {
                "calm_floor_spectral_rad2": self.calm_floor["sigma2_rad2"],
                "after_damping_spectral_rad2": None,
                "ring_band_spectral_rad2":
                    self.rain_band["slope_variance_realised_rad2"],
                "total_spectral_rad2": spec["sigma2_rad2"],
                "total_realised_rad2": st["sigma2_rad2"],
                "ring_share_of_total": (
                    self.rain_band["slope_variance_realised_rad2"]
                    / spec["sigma2_rad2"]) if spec["sigma2_rad2"] > 0 else None,
            }
            out["h_s_spectral_m"] = spec["Hs_m"]
            out["m0_spectral_m2"] = spec["m0_m2"]
        return out

    def band_variance_realised(self):
        """Ring-band height variance of the realised discrete spectrum."""
        if self.rain_band is None:
            return 0.0
        return float(self.rain_band["variance_realised_m2"])
