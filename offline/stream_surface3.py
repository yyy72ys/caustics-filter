"""stream_surface3 -- BRIEF-V4-1-4 budget accounting patch on stream_surface2.

NEW FILE.  `stream_surface2.py` (v4-1.2) is imported and reused UNCHANGED; the
only thing this module does is replace the variance budget, exactly as
DESIGN-V4-1-4 states:

    var(ii)  = phi     * sigma_reg^2
    var(iii) = (1-phi) * sigma_reg^2
    (i)        is computed exactly as in v4-1.2 (physical amplitude, never
               rescaled), and its variance is added to the TOTAL as a sum of
               squares:
    sigma_total^2 = sigma_reg^2 + var_i

Why (v4-1.3 REPORT section 4): the sigma of Dolcetti 2016 that the regression
reproduces was measured on a water surface whose direct bed response (i) is
essentially zero (the 2016 flume sat on a mono-sphere lattice; v4-1.2 measured
var_i/sigma_tot^2 <= 3.14e-04 there).  So the regression must close the
FLOW-BUILT part (ii)+(iii), not the bed-inclusive total.  The bed response is a
separate physical quantity computed through the transfer function T(k), and
combining independent contributions by quadrature (sum of squares) is the
consistent accounting.  v4-1.2 instead did

    var(ii)+var(iii) = max(sigma_tot^2 - var_i, 0)

which at the real-size scene (var_i/sigma_tot^2 = 2.395) clipped to exactly
zero, freezing the surface.  That is the bug this module fixes.

Everything else -- shapes, ring fraction phi, tail, turbulence covariance, the
OU/advection time evolution, the transfer function, the 13-condition
regression, the breaking-limit safety valve, the measurement helpers -- is
inherited verbatim from stream_surface2.StreamState2.  Nothing is re-fitted and
no constant is invented here.

UNIT CHECK (required by the brief, run by the driver):
    setting var_i := 0 in the budget must reproduce v4-1.2's allocation
    var(ii)=phi*sigma_reg^2, var(iii)=(1-phi)*sigma_reg^2 AND the generated
    surface sigma bit-for-bit.  With var_i=0, sigma_total^2 = sigma_reg^2, so
    the two accountings coincide identically; the driver verifies it
    numerically against a replica built with v4-1.2's OWN methods/formula.
"""

import numpy as np

import stream_surface2 as SS

ACCOUNTING = "v4-1.4"
BUDGET_RULE = ("var(ii)=phi*sigma_reg^2; var(iii)=(1-phi)*sigma_reg^2; "
               "sigma_total^2 = sigma_reg^2 + var_i (sum of squares, the bed "
               "response is kept as its own term)")


class StreamState3(SS.StreamState2):
    """v4-1.2 surface generator with the v4-1.4 budget accounting.

    Constructor arguments are exactly StreamState2's, plus
    ``force_var_i_zero`` (unit-check only: treat var_i as 0 in the budget while
    still computing/keeping the physical (i) field).
    """

    def __init__(self, force_var_i_zero=False, **kwargs):
        # v4-1.2 build first: same shapes, same seeds, same transfer, same
        # regression.  Only the budget block is then redone below.
        super().__init__(**kwargs)
        self.force_var_i_zero = bool(force_var_i_zero)

        # `sigma_tot` (inherited name) IS the regression value sigma_reg.
        self.sigma_reg = float(self.sigma_tot)
        self.var_reg = float(self.sigma_reg ** 2)
        self.var_i_physical = float(self.var_i)          # from T(k)*zeta_hat
        self.var_i_budget = 0.0 if self.force_var_i_zero else self.var_i_physical
        self.var_total = float(self.var_reg + self.var_i_budget)
        self.sigma_total = float(np.sqrt(max(self.var_total, 0.0)))
        self._rebuild_budget()

    # -- the only changed physics ------------------------------------------
    def _rebuild_budget(self):
        v = self.var_reg
        self.var_ii = self.phi * v
        self.var_iii = (1.0 - self.phi) * v
        self.var_ring = self.phi * self.var_ii
        self.var_tail = (1.0 - self.phi) * self.var_ii
        if self.k0 is None:
            # no ring: the whole (ii) share is the tail, as in v4-1.2
            self.var_ring = 0.0
            self.var_tail = self.var_ii

        # Same construction order and the same RNG seeds as StreamState2's
        # __init__, so the shapes and phase draws are identical; only the
        # normalising variance differs.  Hence var_i=0 reproduces v4-1.2
        # bit-for-bit.
        rng = np.random.default_rng(self.seed)
        self.hat_ring = self._shape_to_var(
            self._ring_shape(), self.var_ring, rng, "ring")
        self.hat_tail = self._shape_to_var(
            self._tail_shape(), self.var_tail, rng, "tail")
        self.hat_ii0 = self.hat_ring + self.hat_tail
        self.var_ii_real = self._var_of(self.hat_ii0)

        shape_iii = self._turb_shape()
        rngA = np.random.default_rng(self.seed + 0xC0A1)
        rngB = np.random.default_rng(self.seed + 0xB0B1)
        self.hat_turbA = self._normalise_hat(
            shape_iii * SS.hermitian_unit(self.N, rngA), self.var_iii, "turbA")
        self.hat_turbB = self._normalise_hat(
            shape_iii * SS.hermitian_unit(self.N, rngB), self.var_iii, "turbB")
        self.var_iii_real = self._var_of(self.hat_turbA)

    # -- records ------------------------------------------------------------
    def accounting(self):
        d = {
            "accounting": ACCOUNTING,
            "budget_rule": BUDGET_RULE,
            "sigma_reg_m": self.sigma_reg,
            "var_reg_m2": self.var_reg,
            "sigma_total_m": self.sigma_total,
            "var_total_m2": self.var_total,
            "var_i_physical_m2": self.var_i_physical,
            "sigma_i_physical_m": float(np.sqrt(max(self.var_i_physical, 0.0))),
            "var_i_used_in_budget_m2": self.var_i_budget,
            "force_var_i_zero": self.force_var_i_zero,
            "var_i_over_var_reg": (self.var_i_budget / self.var_reg
                                   if self.var_reg > 0 else None),
            "var_i_fraction_of_total": (self.var_i_budget / self.var_total
                                        if self.var_total > 0 else None),
            "var_ii_m2": self.var_ii, "var_iii_m2": self.var_iii,
            "var_ring_m2": self.var_ring, "var_tail_m2": self.var_tail,
            "sigma_ii_m": float(np.sqrt(max(self.var_ii, 0.0))),
            "sigma_iii_m": float(np.sqrt(max(self.var_iii, 0.0))),
            "sigma_ring_m": float(np.sqrt(max(self.var_ring, 0.0))),
            "sigma_tail_m": float(np.sqrt(max(self.var_tail, 0.0))),
            "closure_check": {
                "var_total_minus_reg_minus_var_i": float(
                    self.var_total - self.var_reg - self.var_i_budget),
                "var_ii_plus_iii_over_var_reg": float(
                    (self.var_ii + self.var_iii) / self.var_reg
                    if self.var_reg > 0 else None),
            },
        }
        return d

    def describe(self):
        d = super().describe()
        acc = self.accounting()
        # the inherited describe() still speaks v4-1.2's vocabulary; make the
        # v4-1.4 quantities explicit and unambiguous.
        d.update({
            "accounting": ACCOUNTING,
            "budget_rule": BUDGET_RULE,
            "sigma_reg_m": self.sigma_reg,
            "var_reg_m2": self.var_reg,
            "var_total_target": self.var_total,
            "sigma_total_target_m": self.sigma_total,
            "var_i_physical": self.var_i_physical,
            "var_i_used_in_budget": self.var_i_budget,
            "force_var_i_zero": self.force_var_i_zero,
            "var_i_over_var_total": acc["var_i_fraction_of_total"],
        })
        return d


# ---------------------------------------------------------------------------
# unit check helper: a v4-1.2 replica whose budget is evaluated with var_i = 0
# (v4-1.2's OWN formula and OWN methods, called from outside the class so that
#  stream_surface2 itself stays untouched).
# ---------------------------------------------------------------------------
def v2_replica_var_i_zero(state2):
    """Rebuild a StreamState2's (ii)/(iii) fields with var_i forced to 0.

    Uses v4-1.2's formula  var(ii)=phi*(sigma_tot^2-var_i)  with var_i:=0 and
    v4-1.2's own _shape_to_var / _normalise_hat / _ring_shape / _tail_shape /
    _turb_shape / hermitian_unit, with the same RNG seeds and the same draw
    order.  Nothing here is specific to v4-1.4.
    """
    s2 = state2
    vtot = float(s2.sigma_tot ** 2)
    rem = vtot                      # <-- the "var_i forced to 0" premise
    var_ii = s2.phi * rem
    var_iii = (1.0 - s2.phi) * rem
    var_ring = s2.phi * var_ii
    var_tail = (1.0 - s2.phi) * var_ii
    if s2.k0 is None:
        var_ring = 0.0
        var_tail = var_ii
    rng = np.random.default_rng(s2.seed)
    hat_ring = s2._shape_to_var(s2._ring_shape(), var_ring, rng, "ring")
    hat_tail = s2._shape_to_var(s2._tail_shape(), var_tail, rng, "tail")
    hat_ii0 = hat_ring + hat_tail
    shape_iii = s2._turb_shape()
    rngA = np.random.default_rng(s2.seed + 0xC0A1)
    rngB = np.random.default_rng(s2.seed + 0xB0B1)
    hatA = s2._normalise_hat(
        shape_iii * SS.hermitian_unit(s2.N, rngA), var_iii, "turbA")
    hatB = s2._normalise_hat(
        shape_iii * SS.hermitian_unit(s2.N, rngB), var_iii, "turbB")
    flow = hat_ii0 + hatA            # t = 0: tau_iii = 1, so A + B with B at 0
    return {
        "var_i_used": 0.0,
        "var_ii": float(var_ii), "var_iii": float(var_iii),
        "var_ring": float(var_ring), "var_tail": float(var_tail),
        "hat_ring": hat_ring, "hat_tail": hat_tail, "hat_ii0": hat_ii0,
        "hat_turbA": hatA, "hat_turbB": hatB, "hat_flow": flow,
        "sigma_ii": float(np.sqrt(s2._var_of(hat_ii0))),
        "sigma_iii": float(np.sqrt(s2._var_of(hatA))),
        "sigma_flow": float(np.sqrt(s2._var_of(flow))),
        "sigma_total": float(np.sqrt(s2._var_of(flow + s2.eta_i_hat))),
    }


def unit_check_var_i_zero(state2, state3, rtol=1e-12):
    """Compare StreamState3(force_var_i_zero=True) with the v4-1.2 replica.

    Both must give the same (ii)/(iii) allocation and the same generated
    surface sigma.  Returns a JSON-able record.
    """
    ref = v2_replica_var_i_zero(state2)
    got = {
        "var_ii": float(state3.var_ii), "var_iii": float(state3.var_iii),
        "var_ring": float(state3.var_ring), "var_tail": float(state3.var_tail),
        "sigma_ii": float(np.sqrt(state3._var_of(state3.hat_ii0))),
        "sigma_iii": float(np.sqrt(state3._var_of(state3.hat_turbA))),
        "sigma_flow": float(np.sqrt(state3._var_of(
            state3.hat_ii0 + state3.hat_turbA))),
        "sigma_total": float(np.sqrt(state3._var_of(state3.eta_hat_at(0.0)))),
    }

    def rel(a, b):
        if a == b:
            return 0.0
        d = max(abs(a), abs(b))
        return float(abs(a - b) / d) if d > 0 else float(abs(a - b))

    rows = {}
    for k in ("var_ii", "var_iii", "var_ring", "var_tail",
              "sigma_ii", "sigma_iii", "sigma_flow", "sigma_total"):
        rows[k] = {"v4_1_2_replica": ref[k], "v4_1_4": got[k],
                   "rel_diff": rel(ref[k], got[k]),
                   "equal_within_rtol": bool(rel(ref[k], got[k]) <= rtol)}
    # spectral field identity (same shapes, same phase draws)
    hat_dev = {
        "hat_ring_max_abs_diff": float(np.max(np.abs(
            state3.hat_ring - ref["hat_ring"]))),
        "hat_ii0_max_abs_diff": float(np.max(np.abs(
            state3.hat_ii0 - ref["hat_ii0"]))),
        "hat_turbA_max_abs_diff": float(np.max(np.abs(
            state3.hat_turbA - ref["hat_turbA"]))),
        "hat_turbB_max_abs_diff": float(np.max(np.abs(
            state3.hat_turbB - ref["hat_turbB"]))),
    }
    ok = (all(v["equal_within_rtol"] for v in rows.values())
          and all(v == 0.0 for v in hat_dev.values()))
    return {"pass": bool(ok), "rtol": rtol, "rows": rows,
            "spectral_field_identity": hat_dev,
            "premise": ("var_i forced to 0: v4-1.2 uses "
                        "var(ii)=phi*(sigma_tot^2-var_i) with var_i=0, "
                        "v4-1.4 uses var(ii)=phi*sigma_reg^2 and "
                        "sigma_total^2=sigma_reg^2+var_i, so with var_i=0 the "
                        "two are the same formula")}


def budget_contrast(state2, state3):
    """What ACTUALLY differs between v4-1.2 and v4-1.4 at a real scene."""
    v2_var_ii = float(state2.var_ii)
    v2_var_iii = float(state2.var_iii)
    v2_flow = v2_var_ii + v2_var_iii
    return {
        "var_i_physical_m2": float(state2.var_i),
        "sigma_reg_m": float(state3.sigma_reg),
        "var_reg_m2": float(state3.var_reg),
        "var_i_over_var_reg": float(state2.var_i / state3.var_reg),
        "v4_1_2": {"var_ii": v2_var_ii, "var_iii": v2_var_iii,
                   "var_flow_ii_plus_iii": v2_flow,
                   "sigma_flow": float(np.sqrt(max(v2_flow, 0.0))),
                   "collapse": bool(v2_flow <= 0.0)},
        "v4_1_4": {"var_ii": float(state3.var_ii),
                   "var_iii": float(state3.var_iii),
                   "var_flow_ii_plus_iii": float(state3.var_ii
                                                 + state3.var_iii),
                   "sigma_flow": float(np.sqrt(max(state3.var_ii
                                                   + state3.var_iii, 0.0))),
                   "collapse": False},
        "flow_variance_gain": float((state3.var_ii + state3.var_iii)
                                    / v2_flow) if v2_flow > 0 else None,
    }
