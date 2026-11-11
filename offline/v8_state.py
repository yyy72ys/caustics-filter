"""v8 water state -- chunk 4 (version 3.4): the wave spectrum, decided.

Wraps (does not modify) state.py / v5_state.py / v6_state.py.  The frozen v5,
v6 and v7 paths are untouched: this class is the only new state, and every
change it makes is switchable so the gate-1 breakdown can attribute the
difference to one physics node at a time.

What chunk 4 changes in the state (PIPELINE-V8.md chunk 4, items 1-6):

  1. DISPERSION  omega^2 = (g k + (gamma/rho) k^3) tanh(k h)
     -- the deep-water sqrt(gk) approximation is retired for the wind-driven
     background.  The 8 s loop quantisation is rebuilt on top of it
     (``_omega_q``).  The rain band path of v6 is left exactly as v6 had it
     (this chunk does not touch rain physics: gate 6).

  2. ANCHOR  the Cox-Munk (fully developed) sigma^2 target is replaced by
     the short-fetch slope model of ``v8_physics.anchor_sigma2``, driven by
     the scene fetch F (non-dimensional fetch Ft = gF/U^2).  The taper FORM
     is unchanged (g(k) = 1/(1+(k/kc)^n), n = 3, bisection) so the only
     difference against v5/v6/v7 is the target -- that is what makes gate 1
     an attribution instead of a guess.

  3. DIRECTION  D(theta) becomes an explicit state parameter
     (Mitsuyasu/Longuet-Higgins cos^{2s}), replacing the frozen
     cos^2 + 0.02 lobe.  The reweight is exact and phase-preserving: the
     random draw is untouched, only the per-mode amplitude changes, so the
     realisation keeps the seed of the reference state.

  4. CAPILLARY TAIL  the high-k cutoff rule (viscous damping 2 nu k^2 > omega)
     is applied as a mask and recorded (it does not bind inside the resolved
     band; the hierarchy is reported, not hidden).

  5. BREAKING CLAMP  H_max = min(0.142 lambda, 0.78 h) applied to the spectral
     H_s, with the Stokes 1/7 cap recorded alongside and the activation
     counted (``self.breaking``).

  6. Patch / fetch / model id metadata for numbers-v8.json.

Switches (all default True; used by run_v8.py for the ablation matrix):
  use_direction, use_anchor, use_finite_depth, use_breaking_clamp.
"""
import math

import numpy as np

from v6_state import V6WaterState
from v6_physics import RING_ETA
from state import OMEGA0
from v8_physics import (anchor_sigma2, apply_breaking_clamp, breaking_limits,
                        capillary_tail, capillary_crossover,
                        cos2s_norm, directional_factor, dispersion_domega_dk,
                        dispersion_omega, jonswap_fetch_relations,
                        patch_metadata, spreading_parameters, spreading_s,
                        viscous_cutoff)

TAPER_N = 3.0          # same taper form as v5 (only the target changes)
ANG_INT_OLD = None     # filled lazily: angular integral of cos^2 + 0.02


class V8WaterState(V6WaterState):
    """Chunk-4 spectrum: finite-depth + capillary dispersion, fetch anchor,
    explicit D(theta), capillary cutoff rule, breaking clamp."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 height_mode="wind", height_m=None, rain_rate_mm_h=0.0,
                 fetch_m=None, use_direction=True, use_anchor=True,
                 use_finite_depth=True, use_breaking_clamp=True,
                 fetch_source=None, anchor_target_rad2=None, eta=RING_ETA):
        self.v8_depth = float(recipe_env["depthM"])
        self.v8_wind = float(recipe_env["windSpeedMS"])
        self.v8_fetch = None if fetch_m is None else float(fetch_m)
        self.v8_fetch_source = fetch_source or "explicit"
        self.anchor_override = (None if anchor_target_rad2 is None
                                else float(anchor_target_rad2))
        self.use_direction = bool(use_direction)
        self.use_anchor = bool(use_anchor)
        self.use_finite_depth = bool(use_finite_depth)
        self.use_breaking_clamp = bool(use_breaking_clamp)
        self.direction = None
        self.anchor = None
        self.tail = None
        self.breaking = None

        # The parent builds the frozen Phillips + cos^2 spectrum with the
        # legacy amplitude normalisation.  No Cox-Munk taper is applied here
        # (calm_calibrate=False): the chunk-4 anchor replaces it below.
        super().__init__(recipe, recipe_env, N=N, seed=seed,
                         preset_index=preset_index, height_mode=height_mode,
                         height_m=height_m, rain_rate_mm_h=rain_rate_mm_h,
                         calm_calibrate=False, calm_target_sigma2=None,
                         eta=eta)
        self.v8_legacy_spectrum = self._spectrum_stats()

        if self.use_direction:
            self.direction = self._apply_direction()
        self.tail = self._apply_capillary_cutoff()
        if self.use_anchor and self.height_mode == "wind":
            self.anchor = self._apply_anchor()
        elif not self.use_anchor and self.height_mode == "wind":
            self.anchor = {
                "applied": False,
                "reason": "ablation: the legacy Cox-Munk target is not applied",
                "old_cox_munk_target_rad2": 0.003 + 5.12e-3 * self.v8_wind,
            }
        else:
            self.anchor = {
                "applied": False,
                "reason": ("height-driven mode: a global scalar would undo the "
                           "requested H_s (same rule as v5)"),
            }
        if self.use_breaking_clamp:
            self.breaking = self._apply_breaking()
        else:
            self.breaking = {"applied": False, "reason": "ablation: clamp off"}

        self.spectrum_after = self._spectrum_stats()

    # ------------------------------------------------------------------
    # 1. dispersion
    # ------------------------------------------------------------------
    def _omega_q(self):
        """Quantised omega under the new dispersion (8 s loop preserved).

        The rain branch is delegated to v6 unchanged (chunk 4 does not touch
        rain physics).  ``use_finite_depth=False`` restores the frozen
        deep-water sqrt(gk) for the gate-1 ablation.
        """
        if self.rain_rate > 0.0 or not self.use_finite_depth:
            return super()._omega_q()
        k = np.maximum(self.k.astype(np.float64), 0.0)
        raw = dispersion_omega(k, self.v8_depth)
        n = np.maximum(1.0, np.floor(raw / OMEGA0 + 0.5))
        return (n * OMEGA0).astype(np.float64)

    def dispersion_omega_raw(self):
        """Unquantised omega actually implied by the chosen branch."""
        k = np.maximum(self.k.astype(np.float64), 0.0)
        if self.rain_rate > 0.0:
            from v6_physics import gravity_capillary_omega
            return gravity_capillary_omega(k)
        if not self.use_finite_depth:
            return np.sqrt(9.81 * k)
        return dispersion_omega(k, self.v8_depth)

    def dispersion_node(self):
        raw = self.dispersion_omega_raw()
        q = self._omega_q()
        sl = self._spectrum_density() * (self.k.astype(np.float64) ** 2)
        rel = np.abs(q - raw) / np.maximum(raw, 1e-30)
        m = raw > 0.0
        return {
            "form": ("omega^2 = (g k + (gamma/rho) k^3) tanh(k h)"
                     if (self.use_finite_depth and self.rain_rate <= 0.0)
                     else "legacy deep-water sqrt(g k) (ablation)"),
            "gamma_N_per_m": 0.0728, "rho_kg_per_m3": 1000.0,
            "depth_m": float(self.v8_depth),
            "gamma_over_rho": 7.28e-5,
            "loop_quantum_rad_s": float(OMEGA0),
            "quantisation": {
                "max_rel_error": float(rel[m].max()) if m.any() else None,
                "slope_weighted_mean_rel_error": float((rel[m] * sl[m]).sum()
                                                       / sl[m].sum())
                if sl[m].sum() > 0 else None,
                "m0_weighted_mean_rel_error": float(
                    (rel[m] * self._spectrum_density()[m]).sum()
                    / self._spectrum_density()[m].sum()),
                "modes_below_OMEGA0": int((raw < OMEGA0).sum()),
                "n_modes": int(raw.size),
                "note": ("omega is quantised to integer multiples of "
                         "2*pi/8 s so the 8 s loop closes exactly.  The "
                         "t=0 field does not depend on omega at all (phases "
                         "are exp(i*omega*0)=1), so this table is the only "
                         "place the dispersion can differ at the reference "
                         "still; the time-evolution effect is measured "
                         "separately at t=2.0 s."),
            },
            "capillary_crossover": capillary_crossover(),
            "viscous_cutoff": capillary_tail(self.S, self.N, self.v8_depth),
            "deep_water_retired": bool(self.use_finite_depth
                                       and self.rain_rate <= 0.0),
        }

    # ------------------------------------------------------------------
    # 3. direction
    # ------------------------------------------------------------------
    def _apply_direction(self):
        kx = self.kx.astype(np.float64)
        ky = self.ky.astype(np.float64)
        k = self.k.astype(np.float64)
        wdx, wdy = 1.0, 0.0
        try:
            wdx, wdy = self._wind_dir
        except Exception:
            pass
        theta_w = math.atan2(float(wdy), float(wdx))

        if self.v8_fetch is None:
            raise RuntimeError("direction node needs a fetch")
        sp = spreading_parameters(self.v8_wind, self.v8_fetch)
        om = dispersion_omega(k, self.v8_depth)
        s = spreading_s(om, sp["omega_p_rad_s"], sp["s_max"])
        d_new = directional_factor(kx, ky, theta_w, s)

        theta = np.arctan2(ky, kx)
        c = np.cos(theta - theta_w)
        d_old_raw = np.where(c > 0.0, c * c, 0.0) + 0.02
        th_grid = np.linspace(-math.pi, math.pi, 4096, endpoint=False)
        int_old = float(np.trapz(np.cos(th_grid) ** 2 + 0.02, th_grid))
        ratio = d_new / (d_old_raw / int_old)
        ratio = np.where(np.isfinite(ratio), ratio, 0.0)

        before = self._spectrum_stats()
        self.H0 = (self.H0 * np.sqrt(ratio)).astype(np.complex64)
        mid = self._spectrum_stats()

        # keep the wave HEIGHT fixed: direction must only redistribute slope
        h = np.real(np.fft.ifft2(self.H0.astype(np.complex128)))
        std = float(np.std(h))
        behind = float(self.target_std) / std if std > 0 else 1.0
        self.H0 = (self.H0 * behind).astype(np.complex64)
        del h

        after = self._spectrum_stats()
        return {
            "applied": True,
            "model": ("cos^{2s}(dtheta/2), s = s_max (omega/omega_p)^m, "
                      "m=+5 below / -2.5 above the peak "
                      "(Mitsuyasu et al. 1975; Hasselmann et al. 1980)"),
            "s_max": float(sp["s_max"]), "U_over_cp": float(sp["U_over_cp"]),
            "omega_p_rad_s": float(sp["omega_p_rad_s"]),
            "s_min_over_band": float(s.min()), "s_max_over_band": float(s.max()),
            "old_factor": "cos^2(dtheta) + 0.02 (frozen state.py)",
            "old_angular_integral": int_old,
            "sigma2_before": before["sigma2_rad2"],
            "sigma2_after_reweight_unscaled": mid["sigma2_rad2"],
            "height_renormalisation": behind,
            "sigma2_after": after["sigma2_rad2"],
            "m0_before": before["m0_m2"], "m0_after": after["m0_m2"],
            "Hs_before": before["Hs_m"], "Hs_after": after["Hs_m"],
            "note": ("the random draw is untouched (same seed -> same phases); "
                     "only the per-mode amplitude is reweighted, and the "
                     "height variance is rescaled back to the frozen target "
                     "so the change is purely angular."),
        }

    # ------------------------------------------------------------------
    # 2. anchor
    # ------------------------------------------------------------------
    def _apply_anchor(self):
        if self.v8_fetch is None:
            raise RuntimeError("anchor node needs a fetch")
        node = anchor_sigma2(self.v8_wind, self.v8_fetch, self.v8_depth,
                             S_patch=self.S, N=self.N)
        target = float(node["sigma2_target_rad2"])
        if self.anchor_override is not None:
            target = float(self.anchor_override)
            node["target_override_rad2"] = target
            node["target_override_note"] = (
                "ablation: the v5/v6/v7 Cox-Munk target was used with the "
                "chunk-4 spectrum, so the direction node can be measured in "
                "the second ordering of the gate-1 attribution")
        S0 = self._spectrum_density()
        k = self.k.astype(np.float64)
        sig2_0 = float((k * k * S0).sum())
        if target >= sig2_0:
            node.update({"applied": False, "sigma2_before": sig2_0,
                         "reason": "target >= realised sigma^2; nothing to remove"})
            return node
        lo = self.kmin * 1e-3
        hi = self.kmax_nyq
        for _ in range(90):
            mid = math.sqrt(lo * hi)
            g = 1.0 / (1.0 + (k / mid) ** TAPER_N)
            s2 = float((k * k * S0 * g).sum())
            if s2 > target:
                hi = mid
            else:
                lo = mid
        kc = float(math.sqrt(lo * hi))
        g = 1.0 / (1.0 + (k / kc) ** TAPER_N)
        s2_cal = float((k * k * S0 * g).sum())
        m0_cal = float((S0 * g).sum())
        self.H0 = (self.H0 * np.sqrt(g)).astype(np.complex64)
        node.update({
            "applied": True,
            "sigma2_before": sig2_0,
            "sigma2_after_spectral": s2_cal,
            "sigma2_residual_rel": float(abs(s2_cal - target) / target),
            "m0_before": float(S0.sum()),
            "m0_after_spectral": m0_cal,
            "m0_removed_frac": float(1.0 - m0_cal / float(S0.sum())),
            "Hs_before": float(4.0 * math.sqrt(float(S0.sum()))),
            "Hs_after_spectral": float(4.0 * math.sqrt(m0_cal)),
            "taper": "g(k)=1/(1+(k/kc)^n)", "taper_n": TAPER_N,
            "kc_rad_per_m": kc, "lambda_c_m": float(2.0 * math.pi / kc),
            "kc_over_kmin": float(kc / self.kmin),
            "taper_form_unchanged_from_v5": True,
            "note": ("only the TARGET changed against v5/v6/v7; the taper "
                     "form and solver are identical, so the C difference is "
                     "attributable to the anchor"),
        })
        return node

    # ------------------------------------------------------------------
    # 4. capillary tail
    # ------------------------------------------------------------------
    def _apply_capillary_cutoff(self):
        k = self.k.astype(np.float64)
        kc = viscous_cutoff(self.v8_depth)
        mask = k < kc
        n_removed = int((~mask).sum())
        if n_removed:
            self.H0 = (self.H0 * mask).astype(np.complex64)
        node = capillary_tail(self.S, self.N, self.v8_depth)
        node["modes_removed_in_resolved_band"] = n_removed
        return node

    # ------------------------------------------------------------------
    # 5. breaking clamp
    # ------------------------------------------------------------------
    def _apply_breaking(self):
        self.H0, node = apply_breaking_clamp(
            self.H0, self.N, self.k, self.kmin, self.kmax_nyq,
            self.v8_depth)
        node["stokes_cap_1over7"] = {
            "applied": bool(node["Hs_before_m"] > node["H_stokes_1over7_m"])}
        return node

    # ------------------------------------------------------------------
    # readouts
    # ------------------------------------------------------------------
    def v8_node(self, fields=None):
        """Everything chunk 4 adds, in one dict for numbers-v8.json."""
        out = {
            "name": "wave-spectrum (chunk 4)",
            "fetch_m": self.v8_fetch,
            "fetch_source": self.v8_fetch_source,
            "nondim_fetch_Ft": (None if self.v8_fetch is None else
                                jonswap_fetch_relations(self.v8_wind,
                                                        self.v8_fetch)["Ft"]),
            "patch": patch_metadata(self.S, self.N, self.specL),
            "switches": {
                "use_direction": self.use_direction,
                "use_anchor": self.use_anchor,
                "use_finite_depth": self.use_finite_depth,
                "use_breaking_clamp": self.use_breaking_clamp,
            },
            "dispersion": self.dispersion_node(),
            "anchor": self.anchor,
            "direction": self.direction,
            "capillary_tail": self.tail,
            "breaking": self.breaking,
            "spectrum_legacy_before_chunk4": self.v8_legacy_spectrum,
            "spectrum_after": self.spectrum_after,
            "calm_floor_at_build": self.calm_floor if hasattr(self, "calm_floor")
            else None,
        }
        if fields is not None:
            st = self.realised_stats(fields)
            out["realised"] = st
            out["sigma2_realised_vs_target"] = (
                None if not self.anchor or not self.anchor.get("applied")
                else float((st["sigma2_rad2"] - self.anchor["sigma2_target_rad2"])
                           / self.anchor["sigma2_target_rad2"]))
        return out
