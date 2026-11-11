"""v13 state -- version 3.9: the v11 state builder with the labelled prior.

`V13WaterState` is `v11_state.V11WaterState` with exactly ONE thing changed:
`V11.assemble` is replaced by `V13.assemble`, i.e. the alpha_m resolution
(raw inside the ECKV fit domain, `max(raw, alpha_m_floor(cell))` below it with
the labelled prior).  Every other line of behaviour is inherited verbatim:

    _build / _realise / _direction_node   from V11WaterState
    everything else                        from V10WaterState

so the patch rule, the partition-of-unity cascade, the shared seed family, the
8 s omega quantisation, the breaking post-check and the transport-facing state
node are the untouched v11/v10 code.

`sigma2_h_target` is the prior's only setting.  It defaults to the design's
provisional 0.003 and is exposed so the mandatory sensitivity sweep can move it
over {0, 0.0015, 0.003, 0.006}.  A target of 0 reproduces the v3.8 hard zero
EXACTLY and is the sweep's control point, not a shipped configuration.

`DEFAULT_SOURCE` is the merge gate: PIPELINE-V11's rule flips it to
"elfouhaily-form" only if every criterion holds in this chunk's gate report.
Stage 4 writes that flip (or does not) into this file.
"""
import math

import numpy as np

import v8_physics as P8
import v10_spectrum as V10
import v10_state as S10
import v11_state as S11
import v13_spectrum as V13

OMEGA0 = P8.OMEGA0

#: merge-gated default (PIPELINE-V11's rule).  Stage 4 changes this line only
#: if the rule says so; the report records which way it went and why.
DEFAULT_SOURCE = "jonswap"

#: the design's provisional prior.
DEFAULT_SIGMA2_H_TARGET = V13.SIGMA2_H_TARGET_DEFAULT

V11_SOURCES = V13.SOURCES


class V13WaterState(S11.V11WaterState):
    """Version 3.9 environmental state: v3.8 source + labelled low-wind prior."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 fetch_m=None, height_m=None, k_cut=V10.KCUT_DEFAULT,
                 rain_rate_mm_h=0.0, patch_base_m=None, r_lambda=None,
                 gamma=V10.JONSWAP_GAMMA, theta_w=None,
                 spectrum_source=None,
                 sigma2_h_target=DEFAULT_SIGMA2_H_TARGET):
        src = str(spectrum_source or DEFAULT_SOURCE)
        if src not in V11_SOURCES:
            raise ValueError("unknown spectrum_source %r" % (src,))
        self.spectrum_source = src
        self.sigma2_h_target = float(sigma2_h_target)
        self._v11_ctx = None
        self._v13_ctx = None
        if src != "jonswap":
            # built BEFORE the parent runs so that _build can swap it in
            # without a second, wasted realisation (the v11 route)
            self._v11_ctx = V13.assemble(
                float(recipe_env["windSpeedMS"]), float(fetch_m),
                float(recipe_env["depthM"]), hs_override=height_m,
                k_cut=k_cut, gamma=gamma, source=src,
                sigma2_h_target=self.sigma2_h_target)
            self._v13_ctx = self._v11_ctx
        # the v10 constructor drives the build (v11 only added the switch)
        S10.V10WaterState.__init__(
            self, recipe, recipe_env, N=N, seed=seed,
            preset_index=preset_index, fetch_m=fetch_m, height_m=height_m,
            k_cut=k_cut, rain_rate_mm_h=rain_rate_mm_h,
            patch_base_m=patch_base_m, r_lambda=r_lambda,
            gamma=gamma, theta_w=theta_w)

    # ------------------------------------------------------------------
    def v13_node(self):
        """The v3.9-specific readouts: the alpha_m resolution and its work."""
        if self.spectrum_source == "jonswap":
            return {"spectrum_source": "jonswap",
                    "low_wind_prior": None,
                    "alpha_m_resolution": None,
                    "note": ("legacy source: the v3.9 prior lives in the "
                             "elfouhaily-form branch only, so nothing here "
                             "changes")}
        ctx = self.ctx
        res = ctx["alpha_m_resolution"]
        fp = V13.floor_proof(ctx)
        c1 = V13.c1_three_state(ctx)
        c1b = V13.c1b(ctx, fp)
        return {
            "spectrum_source": self.spectrum_source,
            "version": ctx["version"],
            "alpha_m_resolution": res,
            "low_wind_prior": ctx.get("low_wind_prior"),
            "retired_label_check": {
                "retired": V13.RETIRED_ALPHA_M_DOMAIN,
                "emitted": bool(ctx["alpha_m_domain"]
                                == V13.RETIRED_ALPHA_M_DOMAIN),
                "note": ("the v3.8 label eckv_domain=low_wind_floor is "
                         "retired; v3.9 uses alpha_m_domain"),
            },
            "floor_proof": fp,
            "C1": c1,
            "C1b": c1b,
        }

    def state_node(self, fields=None):
        node = S11.V11WaterState.state_node(self, fields)
        if self.spectrum_source != "jonswap":
            ctx = self.ctx
            node["version"] = "3.9"
            node["class"] = "V13WaterState"
            spec = node.get("spectrum", {})
            spec["alpha_m_used"] = float(ctx["alpha_m_value"])
            spec["alpha_m_floor"] = float(ctx["alpha_m_floor"])
            spec["alpha_m_domain"] = ctx["alpha_m_domain"]
            spec["u_over_fit_min"] = float(
                ctx["alpha_m_resolution"]["u_over_fit_min"])
            spec["in_eckv_fit_domain"] = bool(
                ctx["alpha_m_resolution"]["in_eckv_fit_domain"])
            spec["sigma2_h_target"] = float(ctx["sigma2_h_target"])
            spec["sigma2_h_total"] = float(ctx["sigma2_h_total"])
            spec["low_wind_prior"] = ctx.get("low_wind_prior")
            spec["retired_label_low_wind_floor"] = False
            node["spectrum"] = spec
            node["v13"] = self.v13_node()
        else:
            node["v13"] = self.v13_node()
        return node


# --------------------------------------------------------------------------
# helpers used by the gates (unchanged v11/v10 helpers, re-exported)
# --------------------------------------------------------------------------
def height_quantiles(fields):
    return S11.height_quantiles(fields)


def slope_quantiles(fields):
    return S11.slope_quantiles(fields)
