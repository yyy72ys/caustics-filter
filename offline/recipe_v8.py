"""caustics-recipe v1/v2 -> v8 fetch field (PIPELINE-V8 chunk 4, item 2).

One optional field is added on top of the v6 parser:

    wave.fetchM   number > 0   [m]   (the fetch the wind has blown over)

Rules (PIPELINE-V8 chunk 4):
  * new recipes: the pond default is 100 m (v8_physics.POND_FETCH_M);
  * LEGACY recipes (no fetch field anywhere): the default is the fixed value
    that is CLOSEST to the previous behaviour, i.e. the fetch whose
    short-fetch anchor reproduces the old Cox-Munk target at the legacy
    reference scene.  That value is solved once in
    ``v8_physics.legacy_fetch_m`` (4.60 m) and recorded in numbers-v8.json,
    so a legacy recipe does not silently jump to the new pond physics.
  * the chosen source is carried through to the metadata
    ("explicit" | "pond_default" | "legacy_compat").

A v1 recipe without the new block still reads (v5/v6 parse unchanged); the
only difference is which default fetch is used, and that choice is recorded.
"""
from recipe_v6 import parse_recipe_v6, with_rain
from v8_physics import LEGACY_FETCH_M, POND_FETCH_M

SCHEMA_V8_FIELDS = ["wave.fetchM"]


def _fetch_from_recipe(recipe):
    wave = recipe.get("wave") or {}
    wind = recipe.get("windWave") or {}
    for src, key in ((wave, "fetchM"), (wave, "fetch_m"),
                     (wind, "fetchM"), (wind, "fetch_m")):
        if isinstance(src, dict) and src.get(key) is not None:
            return float(src[key]), "explicit"
    return None, None


def parse_recipe_v8(recipe):
    """v6 parse (unchanged) + the fetch field and its default rule."""
    cfg = parse_recipe_v6(recipe)
    fetch, source = _fetch_from_recipe(recipe)
    if fetch is None:
        legacy = bool(cfg.get("legacy_recipe"))
        fetch = LEGACY_FETCH_M if legacy else POND_FETCH_M
        source = "legacy_compat" if legacy else "pond_default"
    if not (fetch > 0.0):
        raise SystemExit("wave.fetchM must be > 0 (got %r)" % (fetch,))
    cfg = dict(cfg)
    cfg["fetch"] = {
        "fetch_m": float(fetch),
        "source": source,
        "pond_default_m": float(POND_FETCH_M),
        "legacy_default_m": float(LEGACY_FETCH_M),
        "note": ("legacy recipes keep the fetch that reproduces the old "
                 "Cox-Munk anchor; explicit wave.fetchM (pond 100 m) selects "
                 "the chunk-4 short-fetch physics"),
    }
    cfg["schema_v8_fields"] = SCHEMA_V8_FIELDS
    return cfg


def with_v8(recipe, fetch_m=None, depth_m=None, wind_ms=None,
            schema="caustics-recipe/v2"):
    """Deep copy of a recipe with an explicit fetch (and optional depth/wind)."""
    import copy
    r = copy.deepcopy(recipe)
    r["schema"] = schema
    if fetch_m is not None:
        r.setdefault("wave", {})["fetchM"] = float(fetch_m)
    if depth_m is not None:
        r.setdefault("water", {})["depthM"] = float(depth_m)
    if wind_ms is not None:
        r.setdefault("windWave", {})["windSpeedMS"] = float(wind_ms)
    return r
