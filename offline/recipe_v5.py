"""caustics-recipe v1 -> v2 extension parsing (PIPELINE-V5 chunk 1).

The v2 recipe adds OPTIONAL fields; a v1 recipe without them must read exactly
as before ("旧レシピは新フィールド無し＝従来挙動で読めること").

    wave.heightMode   "wind" | "fixed"      (default "wind")
    wave.height_m     number                 (required when heightMode=fixed)
    bottom.substrate  "sand"|"rock"|"vegetation"|"mud"  (default absent=legacy)
    bottom.albedoScale number                (default 1.0)

For convenience the same wave.heightMode / wave.height_m keys are also
accepted inside the legacy "windWave" block; the top-level "wave" block wins.
"""
SCHEMA_V1 = "caustics-recipe/v1"
SCHEMA_V2 = "caustics-recipe/v2"
KNOWN_SUBSTRATES = ("sand", "rock", "vegetation", "mud")


def parse_recipe_v5(recipe):
    schema = recipe.get("schema")
    if schema not in (SCHEMA_V1, SCHEMA_V2):
        raise SystemExit("recipe schema mismatch: %r (expected %s or %s)"
                         % (schema, SCHEMA_V1, SCHEMA_V2))
    wave = recipe.get("wave") or {}
    wind = recipe.get("windWave") or {}
    mode = wave.get("heightMode", wind.get("heightMode", "wind"))
    mode = str(mode or "wind").strip().lower()
    if mode not in ("wind", "fixed"):
        raise SystemExit("wave.heightMode must be 'wind' or 'fixed' (got %r)"
                         % (mode,))
    height_m = wave.get("height_m", wind.get("height_m"))
    if mode == "fixed":
        if height_m is None:
            raise SystemExit("wave.heightMode='fixed' requires wave.height_m")
        height_m = float(height_m)
        if not (height_m > 0.0):
            raise SystemExit("wave.height_m must be > 0 (got %r)" % (height_m,))
    else:
        height_m = None

    bottom = recipe.get("bottom") or {}
    substrate = bottom.get("substrate")
    if substrate is not None:
        substrate = str(substrate).strip().lower()
        if substrate not in KNOWN_SUBSTRATES:
            raise SystemExit("bottom.substrate must be one of %s (got %r)"
                             % (list(KNOWN_SUBSTRATES), substrate))
    albedo_scale = float(bottom.get("albedoScale", 1.0))
    if not (albedo_scale >= 0.0):
        raise SystemExit("bottom.albedoScale must be >= 0")

    return {
        "schema": schema,
        "schema_is_v2": schema == SCHEMA_V2,
        "legacy_recipe": (schema == SCHEMA_V1 and not wave and not bottom),
        "wave": {"heightMode": mode, "height_m": height_m},
        "bottom": {"substrate": substrate, "albedoScale": albedo_scale},
    }
