"""run_v14_state -- version 3.9 boundary redo: the three-arm state builder.

`V14WaterState` is `v13_state.V13WaterState` with ONE thing changed: the ctx is
built by `run_v14_spectrum.assemble(..., arm=...)` instead of `V13.assemble`.
Every other line of behaviour is inherited verbatim from
V11WaterState/V10WaterState (patch rule, cascade, shared seed family, 8 s
omega quantisation, breaking post-check, transport-facing state node).

`spectrum_source` defaults to v13's `DEFAULT_SOURCE` and is NOT touched by this
job (BRIEF-V14 line 21: "Set NO default change in this job").  `arm` is a
measurement switch, default "RAW" (v13 as shipped).
"""
import numpy as np

import v10_spectrum as V10
import v10_state as S10
import v11_spectrum as V11
import v11_state as S11
import v13_state as S13
import run_v14_spectrum as V14

#: read from v13, never re-declared: the merge-gated switch stays where v13
#: left it.  Recorded in every state node so a reader cannot mistake it.
DEFAULT_SOURCE = S13.DEFAULT_SOURCE
DEFAULT_SIGMA2_H_TARGET = V14.SIGMA2_H_TARGET_DEFAULT


class V14WaterState(S11.V11WaterState):
    """Version 3.9 boundary redo: v11 state + the arm's alpha_m resolution."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 fetch_m=None, height_m=None, k_cut=V10.KCUT_DEFAULT,
                 rain_rate_mm_h=0.0, patch_base_m=None, r_lambda=None,
                 gamma=V10.JONSWAP_GAMMA, theta_w=None,
                 spectrum_source=None, arm="RAW",
                 sigma2_h_target=DEFAULT_SIGMA2_H_TARGET):
        src = str(spectrum_source or DEFAULT_SOURCE)
        if src not in V11.SOURCES:
            raise ValueError("unknown spectrum_source %r" % (src,))
        if arm not in V14.ARMS and arm != V14.D_STRICT:
            raise ValueError("unknown arm %r" % (arm,))
        self.spectrum_source = src
        self.arm = arm
        self.sigma2_h_target = float(sigma2_h_target)
        self._v11_ctx = None
        if src != "jonswap":
            # built BEFORE the parent runs so `_build` can swap it in
            # without a second, wasted realisation (the v11/v13 route)
            self._v11_ctx = V14.assemble(
                float(recipe_env["windSpeedMS"]), float(fetch_m),
                float(recipe_env["depthM"]), arm=arm, hs_override=height_m,
                k_cut=k_cut, gamma=gamma, source=src,
                sigma2_h_target=self.sigma2_h_target)
            self._v14_ctx = self._v11_ctx
        # the v10 constructor drives the build
        S10.V10WaterState.__init__(
            self, recipe, recipe_env, N=N, seed=seed,
            preset_index=preset_index, fetch_m=fetch_m, height_m=height_m,
            k_cut=k_cut, rain_rate_mm_h=rain_rate_mm_h,
            patch_base_m=patch_base_m, r_lambda=r_lambda,
            gamma=gamma, theta_w=theta_w)

    # ------------------------------------------------------------------
    def v14_node(self):
        """The v14 readouts: the arm, its provenance and its own gates."""
        if self.spectrum_source == "jonswap":
            return {"spectrum_source": "jonswap", "arm": None,
                    "note": "legacy lane: no alpha_m resolution, no arm"}
        ctx = self.ctx
        res = ctx["alpha_m_resolution"]
        geo = V14.cell_geometry(ctx)
        mom = V14.moments(ctx)
        kcut = float(ctx["k_cut_rad_per_m"])
        k = V11._k_grid(kcut)
        sh = V14.shares(ctx, V11.integral_block(
            k, V11.s_parts(k, ctx), geo["kmin_rad_per_m"], kcut,
            float(ctx["k_p_rad_per_m"])))
        return {
            "spectrum_source": self.spectrum_source,
            "arm": self.arm,
            "arm_definition": V14.arm_definition(self.arm),
            "version": ctx["version"],
            "provenance": ctx["provenance"],
            "alpha_m_resolution": res,
            "low_wind_prior": ctx.get("low_wind_prior"),
            "C1": V14.gate_c1(sh),
            "C1b": V14.gate_c1b(sh, mom),
            "b_envelope": V14.b_envelope(ctx),
            "retired_label_check": {
                "retired": V14.V13.RETIRED_ALPHA_M_DOMAIN,
                "emitted": bool(ctx["alpha_m_domain"]
                                == V14.V13.RETIRED_ALPHA_M_DOMAIN)},
        }

    def state_node(self, fields=None):
        node = S11.V11WaterState.state_node(self, fields)
        ctx = getattr(self, "ctx", None)
        if self.spectrum_source != "jonswap" and ctx:
            node["version"] = "3.9+v14"
            node["class"] = "V14WaterState"
            spec = node.get("spectrum", {})
            spec["v14_arm"] = self.arm
            spec["provenance"] = ctx["provenance"]
            spec["alpha_m_used"] = float(ctx["alpha_m_value"])
            spec["alpha_m_raw"] = float(ctx["alpha_m_raw"])
            spec["alpha_m_floor"] = float(ctx["alpha_m_floor"])
            spec["alpha_m_domain"] = ctx["alpha_m_domain"]
            spec["top_up_acted"] = bool(
                ctx["alpha_m_resolution"]["top_up_acted"])
            spec["sigma2_h_target"] = float(ctx["sigma2_h_target"])
            spec["sigma2_h_total"] = float(ctx["sigma2_h_total"])
            spec["low_wind_prior"] = ctx.get("low_wind_prior")
            node["spectrum"] = spec
        node["v14"] = self.v14_node()
        return node


def height_quantiles(fields):
    return S13.height_quantiles(fields)


def slope_quantiles(fields):
    return S13.slope_quantiles(fields)


def _noop():                                    # pragma: no cover
    return np.zeros(1)
