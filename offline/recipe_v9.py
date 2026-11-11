"""caustics-recipe v9 extension -- sky block (PIPELINE-V8 chunk 5, item 4).

Adds ONE optional block on top of the v8 parser:

    sky: {
      "mode": "manual" | "real",
      "directFraction": 0.8,        # eta = DNI/(DNI+DHI); manual mode
      "turbidity": 3.0,             # Linke, kept as an input (default 3)
      "condition": "clear" | "overcast",   # real mode: which sample to use
      "real": { "eta": 0.8346, ... }        # already-resolved sample (offline)
    }

Rules (design addendum):
  * a recipe with NO sky block is LEGACY: eta = 1 exactly (pure direct), i.e.
    version 3.4 behaviour bit-for-bit.  This is the gate-1 requirement.
  * sky.mode = "manual" (or a bare directFraction): eta is the given number,
    validated to [0, 1].
  * sky.mode = "real": eta comes from the RESOLVED sample (v9_env).  The
    parser never touches the network; if the sample is missing it falls back
    to the design's representative value for the named condition and records
    that it did.
  * the chosen eta, its source and the turbidity are carried into the
    metadata so numbers-v9.json shows where each number came from.
"""
from recipe_v8 import parse_recipe_v8, with_v8
from v9_sky import FALLBACK_ETA, TURBIDITY_DEFAULT

SCHEMA_V9_FIELDS = ["sky.mode", "sky.directFraction", "sky.turbidity",
                    "sky.condition"]
SCHEMA_V9 = "caustics-recipe/v3"
SCHEMA_V9_NOTE = (
    "the frozen v5 parser accepts caustics-recipe/v1 and /v2 only, so a v3 "
    "recipe is normalised down to v2 for the frozen parsers before parsing.  "
    "The sky block is optional: a v2 (or v1) recipe without it is a legacy "
    "recipe with eta = 1.  This is recorded, not silent.")


def _normalise_schema(recipe):
    """v3 -> v2 for the frozen parsers; returns (recipe', record)."""
    sch = recipe.get("schema")
    if sch == SCHEMA_V9:
        r = dict(recipe)
        r["schema"] = "caustics-recipe/v2"
        return r, {"schema_in": sch,
                   "normalised_to": "caustics-recipe/v2",
                   "note": SCHEMA_V9_NOTE}
    return recipe, {"schema_in": sch, "normalised_to": sch}


def _sky_from_recipe(recipe):
    sky = recipe.get("sky")
    return sky if isinstance(sky, dict) else None


def parse_recipe_v9(recipe, resolved_sky=None):
    """v8 parse (unchanged) + the sky block."""
    src, schema_rec = _normalise_schema(recipe)
    cfg = dict(parse_recipe_v8(src))
    sky = _sky_from_recipe(recipe)

    turbidity = TURBIDITY_DEFAULT
    if sky and sky.get("turbidity") is not None:
        turbidity = float(sky["turbidity"])

    if sky is None:
        eta, mode, source = 1.0, "legacy", "legacy_no_sky_block"
        condition = None
        note = ("no sky block: pure direct, eta = 1 exactly.  This is the "
                "backward-compatibility rule and the gate-1 identity vs 3.4.")
    else:
        mode = str(sky.get("mode") or "manual").strip().lower()
        condition = sky.get("condition")
        if mode == "legacy":
            eta, source = 1.0, "explicit_legacy"
            note = "explicit legacy: eta = 1"
        elif mode == "manual":
            df = sky.get("directFraction")
            if df is None:
                raise SystemExit("sky.mode=manual needs sky.directFraction")
            eta = float(df)
            source = "manual"
            note = "manual mode: eta from sky.directFraction"
        elif mode == "real":
            real = sky.get("real")
            if isinstance(real, dict) and real.get("eta") is not None:
                eta = float(real["eta"])
                source = "real_resolved_inline"
                note = "real mode: eta from the inline resolved sample"
            else:
                cond = str(condition or "clear").strip().lower()
                val = None
                if isinstance(resolved_sky, dict):
                    ch = (resolved_sky.get("chosen") or {})
                    rec = ch.get(cond) or {}
                    if rec.get("eta") is not None:
                        val = float(rec["eta"])
                        source = ("real_resolved_cache:" +
                                  str(resolved_sky.get("source")))
                if val is None:
                    val = float(FALLBACK_ETA.get(cond, 0.8))
                    source = "fallback_representative"
                    note = ("real mode but no resolved sample: representative "
                            "value used and recorded (the design says do not "
                            "block on a failed fetch)")
                else:
                    note = "real mode: eta from the resolved sample"
                eta = val
        else:
            raise SystemExit("sky.mode must be manual|real|legacy (got %r)"
                             % (mode,))
        if not (0.0 <= eta <= 1.0):
            raise SystemExit("eta must be in [0,1] (got %r)" % (eta,))

    cfg["sky"] = {
        "mode": mode,
        "eta": float(eta),
        "directFraction": float(eta),
        "turbidity": float(turbidity),
        "condition": condition,
        "source": source,
        "note": note,
        "legacy_identity": bool(sky is None),
        "schema": schema_rec,
        "schema_fields": SCHEMA_V9_FIELDS,
        "resolution_split": {
            "dni_rel": float(eta), "dhi_rel": float(1.0 - eta),
            "definition": "eta = DNI/(DNI+DHI) with DNI+DHI = 1",
        },
    }
    return cfg


def with_v9(recipe, eta=None, mode=None, condition=None, turbidity=None,
            fetch_m=None, depth_m=None, wind_ms=None):
    """Deep copy of a recipe with an explicit sky block (and v8 options)."""
    import copy
    r = with_v8(recipe, fetch_m=fetch_m, depth_m=depth_m, wind_ms=wind_ms,
                schema="caustics-recipe/v3")
    if eta is not None or mode is not None or condition is not None:
        sky = r.setdefault("sky", {})
        if mode is not None:
            sky["mode"] = str(mode)
        if condition is not None:
            sky["condition"] = str(condition)
        if eta is not None:
            sky["directFraction"] = float(eta)
            if sky.get("mode") is None:
                sky["mode"] = "manual"
    if turbidity is not None:
        r.setdefault("sky", {})["turbidity"] = float(turbidity)
    return r
