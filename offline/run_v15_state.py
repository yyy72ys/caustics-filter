"""run_v15_state -- version 3.10 (chunk 11): the merged-default state builder.

`V15WaterState` is `run_v14_state.V14WaterState` with ONE thing changed: the
ctx is built by `run_v15_spectrum.assemble(..., arm=MERGE_ARM)` instead of
v14's RAW arm, and the spectrum source default is the v15 module's
`SPECTRUM_SOURCE_DEFAULT` ('elfouhaily-form').  Every other line of behaviour
is inherited verbatim from V14WaterState / V11WaterState / V10WaterState
(patch rule, cascade, shared seed family, 8 s omega quantisation, breaking
post-check, transport-facing state node).

`arm` is kept as a switch so the rollback/comparison arms can still be realised
(the merge default is run_v15_spectrum.MERGE_ARM = "D", the Arm-D rule set).
`spectrum_source="jonswap"` still takes the legacy lane untouched.
"""
import numpy as np

import v10_spectrum as V10
import v10_state as S10
import v11_spectrum as V11
import v11_state as S11
import v13_state as S13
import run_v14_state as S14
import run_v15_spectrum as V15

#: THE MERGE: this family's default source, read from the v15 spectrum module
#: and never re-declared here.
DEFAULT_SOURCE = V15.SPECTRUM_SOURCE_DEFAULT
DEFAULT_SIGMA2_H_TARGET = V15.SIGMA2_H_TARGET_DEFAULT
DEFAULT_ARM = V15.MERGE_ARM


class V15WaterState(S11.V11WaterState):
    """Version 3.10 merge: v11 state + the merged default spectrum."""

    def __init__(self, recipe, recipe_env, N=2048, seed=None, preset_index=1,
                 fetch_m=None, height_m=None, k_cut=V10.KCUT_DEFAULT,
                 rain_rate_mm_h=0.0, patch_base_m=None, r_lambda=None,
                 gamma=V10.JONSWAP_GAMMA, theta_w=None,
                 spectrum_source=None, arm=DEFAULT_ARM,
                 sigma2_h_target=DEFAULT_SIGMA2_H_TARGET):
        src = str(spectrum_source or DEFAULT_SOURCE)
        if src not in V15.SPECTRUM_SOURCES:
            raise ValueError("unknown spectrum_source %r" % (src,))
        if arm not in V15.ARMS and arm != V15.D_STRICT:
            raise ValueError("unknown arm %r" % (arm,))
        self.spectrum_source = src
        self.arm = arm
        self.sigma2_h_target = float(sigma2_h_target)
        self._v11_ctx = None
        if src != "jonswap":
            # built BEFORE the parent runs so `_build` can swap it in without a
            # second, wasted realisation (the v11/v13/v14 route)
            self._v11_ctx = V15.assemble(
                float(recipe_env["windSpeedMS"]), float(fetch_m),
                float(recipe_env["depthM"]), arm=arm, hs_override=height_m,
                k_cut=k_cut, gamma=gamma, source=src,
                sigma2_h_target=self.sigma2_h_target)
            self._v15_ctx = self._v11_ctx
        # the v10 constructor drives the build
        S10.V10WaterState.__init__(
            self, recipe, recipe_env, N=N, seed=seed,
            preset_index=preset_index, fetch_m=fetch_m, height_m=height_m,
            k_cut=k_cut, rain_rate_mm_h=rain_rate_mm_h,
            patch_base_m=patch_base_m, r_lambda=r_lambda,
            gamma=gamma, theta_w=theta_w)

    # ------------------------------------------------------------------
    def v15_node(self):
        """The v15 readouts: the merge provenance and its own gates."""
        if self.spectrum_source == "jonswap":
            return {"spectrum_source": "jonswap", "arm": None,
                    "note": ("rollback lane: no alpha_m resolution, no arm, "
                             "v11/v10 behaviour untouched")}
        ctx = self.ctx
        res = ctx["alpha_m_resolution"]
        geo = V15.cell_geometry(ctx)
        mom = V15.V14.moments(ctx)
        kcut = float(ctx["k_cut_rad_per_m"])
        k = V11._k_grid(kcut)
        sh = V15.V14.shares(ctx, V11.integral_block(
            k, V11.s_parts(k, ctx), geo["kmin_rad_per_m"], kcut,
            float(ctx["k_p_rad_per_m"])))
        prov = V15.provenance(ctx)
        return {
            "spectrum_source": self.spectrum_source,
            "spectrum_source_default": V15.SPECTRUM_SOURCE_DEFAULT,
            "merge_arm": V15.MERGE_ARM,
            "arm": self.arm,
            "arm_definition": V15.V14.arm_definition(self.arm),
            "version": ctx["version"],
            "provenance": ctx["provenance"],
            "alpha_m_resolution": res,
            "low_wind_prior": ctx.get("low_wind_prior"),
            "omega_c": prov["omega_c"],
            "warn_extrapolated_wave_age": prov["warn_extrapolated_wave_age"],
            "low_wind_cell": prov["low_wind_cell"],
            "low_wind_status": prov["low_wind_status"],
            "C1": V15.V14.gate_c1(sh),
            "C1b": V15.V14.gate_c1b(sh, mom),
            "conformance_a": V15.conformance_gravity(ctx),
            "conformance_c": V15.conformance_capillary(ctx),
            "observation_b": V15.observation_beta(ctx),
            "b_max_observation": V15.b_max_observation(ctx),
            "retired_clause": dict(V15.RETIRED_CLAUSE),
            "retired_label_check": {
                "retired": V15.V13.RETIRED_ALPHA_M_DOMAIN,
                "emitted": bool(ctx["alpha_m_domain"]
                                == V15.V13.RETIRED_ALPHA_M_DOMAIN)},
        }

    def state_node(self, fields=None):
        node = S11.V11WaterState.state_node(self, fields)
        ctx = getattr(self, "ctx", None)
        if self.spectrum_source != "jonswap" and ctx:
            node["version"] = "3.10+v15"
            node["class"] = "V15WaterState"
            spec = node.get("spectrum", {})
            spec["v15_arm"] = self.arm
            spec["spectrum_source_default"] = V15.SPECTRUM_SOURCE_DEFAULT
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
            spec["omega_c"] = float(ctx["omega_c"])
            spec["warn_extrapolated_wave_age"] = bool(
                float(ctx["omega_c"]) > V15.WARN_OMEGA_C)
            node["spectrum"] = spec
        node["v15"] = self.v15_node()
        return node


def height_quantiles(fields):
    return S13.height_quantiles(fields)


def slope_quantiles(fields):
    return S13.slope_quantiles(fields)


def _noop():                                    # pragma: no cover
    return np.zeros(1)
