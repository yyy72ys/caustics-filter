"""v11 state -- version 3.8: a SOURCE SWITCH over the v10 state builder.

What this file is
-----------------
`V11WaterState` is `v10_state.V10WaterState` with exactly ONE thing switched:
the omnidirectional spectrum source.

    spectrum_source = "jonswap"          -> v3.7's two-regime spectrum, the
                                            legacy source, BYTE-IDENTICAL to
                                            v10 (ctx, H0, all readouts: the
                                            parent's code path is called)
    spectrum_source = "elfouhaily-form"  -> S(k) = (B_l + B_h)/k^3 with B_l
                                            from v11_spectrum (PIPELINE-V10
                                            section 3)

Everything downstream is inherited verbatim: the fetch-law amplitude route,
the scene patch and its resolution doublings, the partition-of-unity cascade
plan, the shared seed family, the 8 s omega quantisation, the breaking
post-check, the directional node and the band readouts.  v5-v10 model files
are imported and never edited; the three overridden methods are:

    _build          swap in the v11 ctx before the realisation
    _realise        the SAME realisation with the density call switched
    _direction_node the same estimator on the switched spectrum

`DEFAULT_SOURCE` is the merge gate: PIPELINE-V10 section 5 says the default
flips to "elfouhaily-form" only if every merge criterion holds in this
chunk's gate report.  The challenger is always available behind the switch.
"""
import gc
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10
import v10_state as S10
import v11_spectrum as V11

OMEGA0 = P8.OMEGA0

#: merge-gated default (PIPELINE-V10 section 5).  Stage 4 flips this line only
#: if the merge criteria all hold; the report records which way it went.
DEFAULT_SOURCE = "jonswap"


class V11WaterState(S10.V10WaterState):
    """Version 3.8 environmental state: one spectrum source, switched."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 fetch_m=None, height_m=None, k_cut=V10.KCUT_DEFAULT,
                 rain_rate_mm_h=0.0, patch_base_m=None, r_lambda=None,
                 gamma=V10.JONSWAP_GAMMA, theta_w=None,
                 spectrum_source=None):
        src = str(spectrum_source or DEFAULT_SOURCE)
        if src not in V11.SOURCES:
            raise ValueError("unknown spectrum_source %r" % (src,))
        self.spectrum_source = src
        self._v11_ctx = None
        if src != "jonswap":
            # built BEFORE the parent runs so that _build can swap it in
            # without a second, wasted realisation
            self._v11_ctx = V11.assemble(
                float(recipe_env["windSpeedMS"]), float(fetch_m),
                float(recipe_env["depthM"]), hs_override=height_m,
                k_cut=k_cut, gamma=gamma, source=src)
        super().__init__(recipe, recipe_env, N=N, seed=seed,
                         preset_index=preset_index, fetch_m=fetch_m,
                         height_m=height_m, k_cut=k_cut,
                         rain_rate_mm_h=rain_rate_mm_h,
                         patch_base_m=patch_base_m, r_lambda=r_lambda,
                         gamma=gamma, theta_w=theta_w)

    # ------------------------------------------------------------------
    def _build(self):
        if self._v11_ctx is not None:
            self.ctx = self._v11_ctx
            self.laws = self.ctx["laws"]
            self.hs_source = "predicted"
        super()._build()

    def _realise(self, casc, weight, rng, first):
        if self.spectrum_source == "jonswap":
            return super()._realise(casc, weight, rng, first)
        kx = casc.kx.astype(np.float64)
        ky = casc.ky.astype(np.float64)
        k = casc.k
        E = V11.mode_density_2d(kx, ky, self.ctx, k_cut=None,
                                theta_w=self.theta_w)
        D = E * (casc.dk * casc.dk) * weight
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
        if self.k_cut < k_build:
            H0 = np.where(k <= self.k_cut, H0, np.complex64(0.0))
        casc.H0 = H0
        om = P8.dispersion_omega(np.maximum(k, 1.0e-12), self.depth)
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

    def _direction_node(self):
        if self.spectrum_source == "jonswap":
            return super()._direction_node()
        kc = V11._k_grid(self.k_cut)
        an = V11.anisotropy_prediction(self.ctx, kc, self.kmin, self.k_cut)
        d = V10.eckv_delta(kc, self.depth, self.ctx["ustar"]["ustar_ms"],
                           self.ctx["k_p_rad_per_m"])
        return {
            "model": ("ECKV unified azimuthal function Phi = (1/2pi)(1 + "
                      "Delta cos 2(theta-theta_w)), Delta from Elfouhaily "
                      "et al. (1997) eq. 66-68, applied to the v3.8 total"),
            "theta_w_deg": float(math.degrees(self.theta_w)),
            "a0": V10.ECKV_A0, "a_p": V10.ECKV_AP,
            "a_m": float(V10.ECKV_AM_K * self.ctx["ustar"]["ustar_ms"]
                         / V10.C_M),
            "delta_effective": float(an["delta_effective"]),
            "delta_min_over_band": float(np.min(d)),
            "delta_max": float(np.max(d)),
            "sigma2_along_predicted": float(an["sigma2_along"]),
            "sigma2_cross_predicted": float(an["sigma2_cross"]),
            "ratio_predicted": float(an["ratio_along_over_cross"]),
            "note": ("analytic on the assembled continuum (k ~ 1e-3 to "
                     "k_cut); the realised field's own gx/gy ratio is "
                     "reported separately in state_node.realised"),
        }

    # ------------------------------------------------------------------
    def k1_node(self):
        """K1: predicted H_s against the fetch-law H_s, with the 35% flag."""
        hs_law = float(self.laws["Hs_depth_limited_m"])
        hs_pred = float(self.spectrum_after["Hs_m"])
        if self.spectrum_source == "jonswap":
            hs_pred_src = float(self.ctx["Hs_m"])
        else:
            hs_pred_src = float(self.ctx["Hs_m"])
        d = float(hs_pred / hs_law - 1.0)
        return {
            "gate": "K1 -- Delta H_s = H_s_pred - H_s_law",
            "Hs_law_m": hs_law,
            "Hs_law_form_0p0016_m": float(self.laws["Hs_fetch_law_m"]),
            "Hs_predicted_m": hs_pred,
            "Hs_predicted_assembled_m": hs_pred_src,
            "delta_rel": d,
            "delta_abs_m": float(hs_pred - hs_law),
            "flag_threshold_rel": 0.35,
            "flag": bool(abs(d) > 0.35),
            "hs_override_m": self.ctx.get("hs_override_m"),
            "hs_source": self.hs_source,
            "note": ("H_s under the Elfouhaily-form source is a PREDICTION "
                     "of the assembled spectrum (4 sqrt m0), never rescaled "
                     "toward the law; a recipe H_s override is carried as a "
                     "statement, not applied (PIPELINE-V10 section 4)"),
        }

    def v11_node(self):
        """The v3.8-specific readouts: B envelopes, domain labels, K1."""
        if self.spectrum_source == "jonswap":
            return {"spectrum_source": "jonswap",
                    "note": "legacy source: nothing v3.8-specific to report"}
        ctx = self.ctx
        khi = float(P8.viscous_cutoff(self.depth))
        k = V11._k_grid(khi, k_hi=khi)
        s = V11.s_parts(k, ctx)
        m = k >= self.kmin
        out = {
            "spectrum_source": self.spectrum_source,
            "label": ctx["label"],
            "sea_state": ctx["sea_state"],
            "gamma_bl": float(ctx["gamma_bl"]),
            "alpha_p": float(ctx["alpha_p"]),
            "sigma_bl": float(ctx["sigma_bl"]),
            "omega_c": float(ctx["omega_c"]),
            "eckv_domain": ctx["eckv_domain"],
            "alpha_m_domain": ctx.get("alpha_m_domain"),
            "B_l_max_full": float(np.max(s["B_l"])),
            "B_h_max_full": float(np.max(s["B_h"])),
            "B_total_max_full": float(np.max(s["B_l"] + s["B_h"])),
            "k_at_B_total_max": float(k[int(np.argmax(s["B_l"] + s["B_h"]))]),
            "sigma2_model_full": float(np.trapz(k[m] ** 2 * s["S"][m],
                                                k[m])),
            "m0_model_full": float(np.trapz(s["S"][m], k[m])),
            "Hs_model_full_m": float(4.0 * math.sqrt(
                max(float(np.trapz(s["S"][m], k[m])), 0.0))),
            "K1": self.k1_node(),
        }
        return out

    def state_node(self, fields=None):
        node = super().state_node(fields)
        node["spectrum_source"] = self.spectrum_source
        if self.spectrum_source != "jonswap":
            node["spectrum"] = {
                "convention": "S(k), m0=int S dk, sigma2=int k^2 S dk",
                "S_l": "Elfouhaily B_l(k; k_p, Omega_c from the fetch laws)/k^3",
                "S_h": "Elfouhaily B_h(k; u*)/k^3 (v3.7 part, unchanged)",
                "additive_no_splice": True,
                "alpha_m": node["spectrum"]["alpha_m"],
                "alpha_m_raw": node["spectrum"]["alpha_m_raw"],
                "alpha_p": float(self.ctx["alpha_p"]),
                "gamma_bl": float(self.ctx["gamma_bl"]),
                "sigma_bl": float(self.ctx["sigma_bl"]),
                "omega_c": float(self.ctx["omega_c"]),
                "eckv_domain": self.ctx["eckv_domain"],
                "alpha_m_domain": self.ctx.get("alpha_m_domain"),
                "ustar_ms": float(self.ctx["ustar"]["ustar_ms"]),
                "u_over_cm": float(self.ctx["alpha_m"]["u_over_cm"]),
                "k_cut_rad_per_m": float(self.k_cut),
                "k_cut_max_used": float(self.k_cut_max),
                "hs_source": "predicted",
                "hs_override_m": self.ctx.get("hs_override_m"),
            }
            node["v11"] = self.v11_node()
        return node


# --------------------------------------------------------------------------
# helpers used by the gates (K2 / K5)
# --------------------------------------------------------------------------
def height_quantiles(fields):
    """K2: quantiles of the built height field, plus its own m0 / H_s."""
    h = np.asarray(fields["h"], dtype=np.float64)
    m0 = float(np.mean(h * h))
    return {
        "p50_m": float(np.percentile(h, 50)),
        "p90_m": float(np.percentile(h, 90)),
        "p99_m": float(np.percentile(h, 99)),
        "p01_m": float(np.percentile(h, 1)),
        "min_m": float(h.min()), "max_m": float(h.max()),
        "mean_m": float(h.mean()), "std_m": float(h.std()),
        "m0_m2": m0, "Hs_m": float(4.0 * math.sqrt(max(m0, 0.0))),
        "skewness": float(np.mean(((h - h.mean()) / max(h.std(), 1e-30))
                                  ** 3)),
        "kurtosis": float(np.mean(((h - h.mean()) / max(h.std(), 1e-30))
                                  ** 4)),
    }


def slope_quantiles(fields):
    """The same for the two slope components (transport-visible)."""
    gx = np.asarray(fields["gx"], dtype=np.float64)
    gy = np.asarray(fields["gy"], dtype=np.float64)
    s = np.hypot(gx, gy)
    return {"p50": float(np.percentile(s, 50)),
            "p90": float(np.percentile(s, 90)),
            "p99": float(np.percentile(s, 99)),
            "mean": float(s.mean())}
