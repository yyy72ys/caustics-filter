"""caustics-recipe v2 -> v6 rain field (PIPELINE-V6 chunk 2, section 5).

One optional block is added on top of the v5 parser:

    rain.rate_mm_h   number >= 0   (default 0)

R = 0 must read EXACTLY as v5 (PIPELINE-V6 section 5: "0 のとき v5 と完全同一
の挙動"), and a v1 recipe without any new field must read as before.  Existing
fields are not changed; parse_recipe_v5 does the v5 half of the job and is
called unchanged.
"""
from recipe_v5 import parse_recipe_v5

# Radar-based rain categories quoted in report section 3 (for the readout
# only, never for clamping the input).
RAIN_CATEGORIES_MM_H = {
    "light": (0.2, 2.0),
    "moderate": (2.0, 10.0),
    "heavy": (10.0, 30.0),
    "very_heavy": (30.0, 50.0),
    "intense": (50.0, 100.0),
}


def rain_category(rate_mm_h):
    R = float(rate_mm_h)
    if R <= 0.0:
        return "none"
    for name, (lo, hi) in RAIN_CATEGORIES_MM_H.items():
        if lo <= R <= hi:
            return name
    return "above_intense" if R > 100.0 else "below_light"


def parse_recipe_v6(recipe):
    """v5 parse (unchanged) + the rain block."""
    cfg = parse_recipe_v5(recipe)
    rain = recipe.get("rain") or {}
    if not isinstance(rain, dict):
        raise SystemExit("rain must be an object (got %r)" % (rain,))
    rate = rain.get("rate_mm_h", 0.0)
    try:
        rate = float(rate)
    except Exception:
        raise SystemExit("rain.rate_mm_h must be a number (got %r)"
                         % (rate,))
    if not (rate == rate) or rate in (float("inf"), float("-inf")):
        raise SystemExit("rain.rate_mm_h must be finite (got %r)" % (rate,))
    if rate < 0.0:
        raise SystemExit("rain.rate_mm_h must be >= 0 (got %r)" % (rate,))
    cfg = dict(cfg)
    cfg["rain"] = {
        "rate_mm_h": float(rate),
        "active": bool(rate > 0.0),
        "category": rain_category(rate),
        "legacy": (rate == 0.0 and not rain),
        "note": ("rate 0 == v5 behaviour exactly; the field is optional and "
                 "absent in v1 recipes"),
    }
    cfg["legacy_recipe"] = bool(cfg.get("legacy_recipe")
                                and not rain and rate == 0.0)
    cfg["schema_v6_fields"] = ["rain.rate_mm_h"]
    return cfg


def with_rain(recipe, rate_mm_h, wind_ms=0.0, schema="caustics-recipe/v2"):
    """Deep copy of a recipe with an explicit rain block (and optional wind)."""
    import copy
    r = copy.deepcopy(recipe)
    r["schema"] = schema
    r.setdefault("rain", {})["rate_mm_h"] = float(rate_mm_h)
    if wind_ms is not None:
        r.setdefault("windWave", {})["windSpeedMS"] = float(wind_ms)
    return r
