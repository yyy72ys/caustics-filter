"""v5 water state: H_s readout / height-driven mode + calmness (sigma^2).

Wraps the frozen Phillips-FFT state (state.py) without changing it.  Two
PIPELINE-V5 chunk-1 variables are added on top:

  1. Significant wave height H_s = 4*sqrt(m0).
     - "wind"  : the existing Phillips spectrum is kept; H_s is a READ-ONLY
                 derived node (numeric integral of the realised spectrum).
     - "fixed" : one scalar multiplier on the spectrum makes m0 = (H_s/4)^2
                 while the peak wavenumber stays where wind/fetch put it
                 (two-parameter Bretschneider-type form).  The Stokes limit
                 H/lambda <= 1/7 caps the realised height; the cap is
                 reported, never silently exceeded.

  2. Calmness = mean square slope sigma^2 = integral k^2 S(k) dk.
     - read out per frame from the FFT gradients (realised sigma^2),
       with the Cox-Munk reference 0.003 + 5.12e-3 U shown alongside.
     - in wind mode the high-wavenumber side of the spectrum is calibrated so
       the realised sigma^2 matches Cox-Munk.  The calibration is a smooth
       taper g(k) = 1/(1+(k/kc)^n) solved by bisection; the residual
       deviation inside the resolved band is reported in numbers-v5.json.
     - in height mode the calibration is NOT applied (a global scale would
       undo it); the realised sigma^2 deviation from Cox-Munk is reported.

Nothing in state.py is modified; the class subclasses WaterState and only
rescales H0 before the existing FFT path is used.
"""
import numpy as np

from state import WaterState
from v5_physics import (cox_munk_sigma2, cox_munk_anisotropy,
                        glitter_slope_sigmas, stokes_cap_hs)

TAPER_N = 3.0  # smoothness exponent of the calmness taper (documented)


class V5WaterState(WaterState):
    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 height_mode="wind", height_m=None, calm_calibrate=True,
                 calm_target_sigma2=None):
        super().__init__(recipe, recipe_env, N=N, seed=seed,
                         preset_index=preset_index)
        self._wind = float(recipe_env["windSpeedMS"])
        self._wind_dir = recipe_env.get("windDir", [1.0, 0.0])
        self.height_mode = str(height_mode)
        if self.height_mode not in ("wind", "fixed"):
            raise ValueError("wave.heightMode must be 'wind' or 'fixed'")
        self.height_m_requested = (None if height_m is None
                                   else float(height_m))
        self.calm_calibrate = bool(calm_calibrate)
        self.calm_target_sigma2 = calm_target_sigma2

        pre = self._spectrum_stats()
        self.spectrum_before = pre
        self.calm = None
        self.height = None

        if self.height_mode == "wind" and self.calm_calibrate:
            target = (float(calm_target_sigma2)
                      if calm_target_sigma2 is not None
                      else cox_munk_sigma2(recipe_env["windSpeedMS"]))
            self.calm = self._apply_calm_taper(target)
        elif self.height_mode == "wind":
            target = (float(calm_target_sigma2)
                      if calm_target_sigma2 is not None
                      else cox_munk_sigma2(recipe_env["windSpeedMS"]))
            self.calm = {"applied": False, "target_sigma2": target,
                         "reason": "calibration disabled"}

        if self.height_mode == "fixed":
            if self.height_m_requested is None or self.height_m_requested <= 0:
                raise ValueError("wave.height_m must be > 0 in fixed mode")
            self.height = self._apply_height(self.height_m_requested)

        self.spectrum_after = self._spectrum_stats()

    # -- helpers ---------------------------------------------------------
    def _spectrum_density(self):
        """Expected spectral density per mode: |H0|^2/(2 N^4) (m^2)."""
        N = self.N
        return (np.abs(self.H0.astype(np.complex128)) ** 2) / (2.0 * N ** 4)

    def _spectrum_stats(self):
        S = self._spectrum_density()
        k = self.k.astype(np.float64)
        m0 = float(S.sum())
        sig2 = float((k * k * S).sum())
        peak = self._spectral_peak(k, S)
        return {"m0_m2": m0, "Hs_m": float(4.0 * np.sqrt(m0)),
                "sigma2_rad2": sig2, "rms_slope": float(np.sqrt(sig2)),
                "peak": peak}

    def _spectral_peak(self, k, S):
        """Dominant wavelength: peak of the radially integrated spectrum.

        Bins are log-spaced over the resolved band; each bin sum is divided
        by its width in k, so the reported peak is the maximum of S(k) in
        k-space (not per log-k).
        """
        kf = k.ravel()
        Sf = S.ravel()
        nb = 400
        lo = max(self.kmin * 0.5, 1e-6)
        hi = self.kmax_nyq
        edges = np.logspace(np.log10(lo), np.log10(hi), nb + 1)
        idx = np.digitize(kf, edges) - 1
        ok = (idx >= 0) & (idx < nb)
        acc = np.bincount(idx[ok], weights=Sf[ok], minlength=nb)
        dk = np.diff(edges)
        dens = acc / np.maximum(dk, 1e-30)
        i = int(np.argmax(dens))
        kp = float(np.sqrt(edges[i] * edges[i + 1]))
        return {"k_peak_rad_per_m": kp,
                "lambda_peak_m": float(2.0 * np.pi / kp),
                "variance_density_peak": float(dens[i]),
                "kmin_rad_per_m": float(self.kmin),
                "kmax_rad_per_m": float(self.kmax_nyq),
                "peak_is_at_kmin": bool(kp <= self.kmin * 1.35)}

    def _apply_calm_taper(self, target):
        """Solve kc so that integral k^2 S/(1+(k/kc)^n) == target."""
        S0 = self._spectrum_density()
        k = self.k.astype(np.float64)
        sig2_0 = float((k * k * S0).sum())
        if not np.isfinite(target) or target <= 0.0:
            raise ValueError("calm target sigma^2 must be > 0")
        if target >= sig2_0:
            return {"applied": False, "target_sigma2": float(target),
                    "sigma2_before": sig2_0,
                    "reason": "target >= realised sigma^2; nothing to remove"}
        lo = self.kmin * 1e-3
        hi = self.kmax_nyq
        for _ in range(90):
            mid = np.sqrt(lo * hi)
            g = 1.0 / (1.0 + (k / mid) ** TAPER_N)
            s2 = float((k * k * S0 * g).sum())
            if s2 > target:
                hi = mid
            else:
                lo = mid
        kc = float(np.sqrt(lo * hi))
        g = 1.0 / (1.0 + (k / kc) ** TAPER_N)
        s2_cal = float((k * k * S0 * g).sum())
        m0_cal = float((S0 * g).sum())
        self.H0 = (self.H0 * np.sqrt(g)).astype(np.complex64)
        return {
            "applied": True,
            "target_sigma2": float(target),
            "sigma2_before": sig2_0,
            "sigma2_after_spectral": s2_cal,
            "sigma2_residual_rel": float(abs(s2_cal - target) / target),
            "m0_before": float(S0.sum()),
            "m0_after_spectral": m0_cal,
            "m0_removed_frac": float(1.0 - m0_cal / float(S0.sum())),
            "Hs_before": float(4.0 * np.sqrt(float(S0.sum()))),
            "Hs_after_spectral": float(4.0 * np.sqrt(m0_cal)),
            "taper": "g(k)=1/(1+(k/kc)^n)",
            "taper_n": TAPER_N,
            "kc_rad_per_m": kc,
            "lambda_c_m": float(2.0 * np.pi / kc),
            "kc_over_kmin": float(kc / self.kmin),
            "note": ("the Cox-Munk target removes most of the patch-scale "
                     "slope variance; m0 (and therefore H_s) necessarily "
                     "drops.  Reported, not hidden."),
        }

    def _apply_height(self, hs_target):
        """Global scalar multiplier so that m0 = (H_s/4)^2, then Stokes cap."""
        # By construction of WaterState, the realised std equals target_std.
        m0_cur = float(self.target_std) ** 2
        lam_p = self.spectrum_before["peak"]["lambda_peak_m"]
        hs_cap = stokes_cap_hs(lam_p)
        factor = (float(hs_target) / 4.0) / float(self.target_std)
        hs_raw = float(4.0 * np.sqrt(m0_cur * factor * factor))
        capped = bool(hs_raw > hs_cap)
        if capped:
            factor *= hs_cap / hs_raw
        self.H0 = (self.H0 * factor).astype(np.complex64)
        hs_real = float(4.0 * np.sqrt(m0_cur * factor * factor))
        return {
            "mode": "fixed", "height_m_requested": float(hs_target),
            "factor": float(factor),
            "Hs_before": float(4.0 * np.sqrt(m0_cur)),
            "Hs_uncapped_requested": hs_raw,
            "Hs_after": hs_real,
            "stokes_cap": {
                "applied": capped, "H_over_lambda_limit": 1.0 / 7.0,
                "lambda_peak_m": float(lam_p),
                "Hs_cap_m": float(hs_cap),
                "note": ("cap uses the dominant wavelength of the realised "
                         "spectrum; the cap is recorded, never silently "
                         "exceeded"),
            },
            "peak_kept": True,
        }

    # -- readouts --------------------------------------------------------
    def calmness_node(self, fields=None):
        """Node 'calmness' = sigma^2 with the Cox-Munk reference."""
        target = (self.calm_target_sigma2 if self.calm_target_sigma2
                  is not None else cox_munk_sigma2(self._wind))
        out = {
            "name": "calmness",
            "units": "rad^2",
            "definition": "sigma^2 = integral k^2 S(k) dk over the resolved "
                          "band (variance of the surface gradient)",
            "sigma2_spectral": self.spectrum_after["sigma2_rad2"],
            "rms_slope_rad": self.spectrum_after["rms_slope"],
            "rms_slope_deg": float(np.degrees(
                np.arctan(self.spectrum_after["rms_slope"]))),
            "cox_munk_reference": target,
            "cox_munk_anisotropy": cox_munk_anisotropy(self._wind),
            "glitter_slope_sigmas": glitter_slope_sigmas(
                self._wind, float(np.arctan2(self._wind_dir[1],
                                             self._wind_dir[0]))),
            "note": ("a 0-1 'calmness' display would be a labelled transform "
                     "of sigma^2, not physics (report section 7)"),
        }
        if fields is not None:
            st = self.realised_stats(fields)
            out["sigma2_realised"] = st["sigma2_rad2"]
            out["rms_slope_realised"] = st["rms_slope"]
            out["Hs_realised"] = st["Hs_m"]
            out["m0_realised"] = st["m0_m2"]
            out["deviation_vs_cox_munk_rel"] = float(
                (st["sigma2_rad2"] - target) / target) if target else None
        return out

    def realised_stats(self, fields):
        h = np.asarray(fields["h"], dtype=np.float64)
        gx = np.asarray(fields["gx"], dtype=np.float64)
        gy = np.asarray(fields["gy"], dtype=np.float64)
        m0 = float(np.mean(h * h))
        sig2 = float(np.mean(gx * gx + gy * gy))
        return {"m0_m2": m0, "Hs_m": float(4.0 * np.sqrt(max(m0, 0.0))),
                "sigma2_rad2": sig2,
                "rms_slope": float(np.sqrt(max(sig2, 0.0)))}

    def wave_node(self):
        U = self._wind
        return {
            "name": "H_s",
            "units": "m",
            "definition": "significant wave height H_s = 4*sqrt(m0), "
                          "m0 = integral S(k) dk of the realised spectrum",
            "mode": self.height_mode,
            "wind_speed_ms": float(U),
            "Hs_spectral": self.spectrum_after["Hs_m"],
            "m0_spectral": self.spectrum_after["m0_m2"],
            "spectrum_before": self.spectrum_before,
            "spectrum_after": self.spectrum_after,
            "height_driven": self.height,
            "calmness_calibration": self.calm,
            "stokes_limit": {"H_over_lambda": 1.0 / 7.0},
        }
