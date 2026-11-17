"""presets_env -- the numbered ENVIRONMENT preset table (new, offline).

Brief: BRIEF-PRESET-TABLE.md (D:\\mutaro-work\\rally-in).

Why this module exists
----------------------
BRIEF-PRESET5 asked for preset 5 "深い澄んだ池・風" inside a numbered table
carrying `U / F / z / medium / phase n, mu / eta / spectrum`.  That table did
not exist in this tree: the only environment-preset list, `web/recipe.js`
`presets()`, is display-only and cannot carry `F`, a medium id, the FF
`n / mu`, `eta`, or the spectrum source, and substituting values into it was
forbidden.  This module is the missing table, in the brief's form.

The table is a TABLE, not a second physics engine.  Every physical number is
resolved by the existing modules, through the same code path the
BRIEF-PRESET5 reference resolution used:

    run_v15_spectrum.assemble(U, fetch, z, arm=..., source=...)
    v12_water.preset_spectra(medium)          -> medium preset / label / load
    run_v16_water.FF_PHASE[phase_key]         -> Fournier-Forand n, mu
    recipe_v9.parse_recipe_v9(recipe)         -> eta (sky block / legacy rule)

No value is rounded, substituted, or filled in by guesswork.  If a key is
absent, `get_preset_env` raises; it never invents an entry.

Schema of one entry (exactly these ten fields, no extras):

    {id: int, name: str, U: float, fetch: float, z: float,
     medium: str, phase_key: str, eta: float, spectrum: str, arm: str}

Only id = 5 is registered.  Ids 0-4 belong to the display-only web list and
are deliberately NOT copied here: this module must not invent the physical
values that list cannot carry.
"""

import json
import os

SCHEMA_FIELDS = ("id", "name", "U", "fetch", "z",
                 "medium", "phase_key", "eta", "spectrum", "arm")

#: The registered entries.  Keyed by the preset id, and the entry's own
#: "id" must agree with the key (checked at import time, below).
PRESETS_ENV = {
    5: {
        "id": 5,
        "name": "深い澄んだ池・風",
        "U": 6.0,
        "fetch": 100.0,
        "z": 4.0,
        "medium": "clear",
        "phase_key": "clear",
        "eta": 1.0,
        "spectrum": "elfouhaily-form",
        "arm": "D",
    },
}

#: The recipe used ONLY to resolve `eta` through the v9 sky-block parser.
#: It is a legacy recipe (no `sky` block), so `parse_recipe_v9` returns
#: eta = 1.0 with source "legacy_no_sky_block" -- the "sun only" condition of
#: this preset.  The path is the same one BRIEF-PRESET5's reference
#: resolution used.  The value reported by `resolve_preset_env` is the
#: PARSER's output; the table's `eta` is only the expectation.
ETA_RECIPE_REL = os.path.join("..", "blender", "recipe-1.json")


def _validate_entry(entry):
    """Fail closed on a malformed table row (no silent coercion)."""
    got = set(entry)
    want = set(SCHEMA_FIELDS)
    if got != want:
        raise ValueError("preset entry schema mismatch: missing %s, extra %s"
                         % (sorted(want - got), sorted(got - want)))
    if not isinstance(entry["id"], int):
        raise TypeError("preset id must be int (got %r)" % (entry["id"],))
    if not isinstance(entry["name"], str):
        raise TypeError("preset name must be str (got %r)" % (entry["name"],))
    for k in ("U", "fetch", "z", "eta"):
        if not isinstance(entry[k], float):
            raise TypeError("preset %r must be float (got %r)"
                            % (k, entry[k]))
    for k in ("medium", "phase_key", "spectrum", "arm"):
        if not isinstance(entry[k], str):
            raise TypeError("preset %r must be str (got %r)"
                            % (k, entry[k]))
    return entry


for _k, _v in PRESETS_ENV.items():
    _validate_entry(_v)
    if _v["id"] != _k:
        raise ValueError("preset key %r disagrees with entry id %r"
                         % (_k, _v["id"]))
del _k, _v


# --------------------------------------------------------------------------
# table access
# --------------------------------------------------------------------------
def preset_env_ids():
    """The registered ids, ascending.  Today: [5]."""
    return sorted(PRESETS_ENV)


def get_preset_env(preset_id):
    """One table entry, as a copy.  Raises KeyError for an unregistered id.

    There is deliberately no default and no nearest-match: an absent id is
    reported as absent, never filled in.
    """
    try:
        key = int(preset_id)
    except (TypeError, ValueError):
        raise KeyError("preset id must be an integer (got %r)" % (preset_id,))
    if key not in PRESETS_ENV:
        raise KeyError("no environment preset with id %r (registered: %s)"
                       % (preset_id, preset_env_ids()))
    return dict(PRESETS_ENV[key])


# --------------------------------------------------------------------------
# resolution (the existing stack, not new physics)
# --------------------------------------------------------------------------
def _eta_recipe_path():
    return os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     ETA_RECIPE_REL))


def resolve_preset_env(preset_id):
    """Resolve a table entry into the full parameter dict.

    Returns the table fields plus everything the production stack derives
    from them, at least::

        U, fetch, z, medium_preset, medium_label, medium_load,
        phase_n, phase_mu, eta, spectrum, arm, omega_c,
        warn_extrapolated_wave_age

    Every one of those comes from a call into the existing modules (see the
    module docstring for the four call sites).  Nothing here recomputes
    physics and nothing is rounded.
    """
    import run_v15_spectrum as V15S
    import v12_water as W
    import run_v16_water as V16
    import recipe_v9 as R9

    entry = get_preset_env(preset_id)

    # 1. sea state + spectrum: run_v15_spectrum.assemble (the v14/v13/v11
    #    assembly under the v3.10 merge: elfouhaily-form + Arm-D rules).
    ctx = V15S.assemble(entry["U"], entry["fetch"], entry["z"],
                        arm=entry["arm"], source=entry["spectrum"])

    # 2. water medium: v12's preset -> a/b spectra.
    medium = W.preset_spectra(entry["medium"])

    # 3. phase function: v16's published Fournier-Forand fits.
    phase = V16.FF_PHASE[entry["phase_key"]]

    # 4. eta: the v9 sky-block rule.  A legacy recipe (no sky block) gives
    #    eta = 1 exactly; the value below is the parser's, not the table's.
    recipe_path = _eta_recipe_path()
    if not os.path.isfile(recipe_path):
        raise FileNotFoundError(
            "eta recipe missing: %s (the eta path is recipe_v9."
            "parse_recipe_v9; it is not substituted)" % (recipe_path,))
    with open(recipe_path, "r", encoding="utf-8") as fh:
        recipe = json.load(fh)
    parsed = R9.parse_recipe_v9(recipe)
    sky = parsed["sky"]

    out = {
        # --- table fields -------------------------------------------------
        "id": entry["id"],
        "name": entry["name"],
        "U": float(entry["U"]),
        "fetch": float(entry["fetch"]),
        "z": float(entry["z"]),
        "medium": entry["medium"],
        "phase_key": entry["phase_key"],
        "spectrum": str(ctx["spectrum_source"]),
        "arm": str(ctx["merge_arm"]),
        "eta": float(sky["eta"]),
        # --- what the table asked for (for divergence reporting) ----------
        "eta_table": float(entry["eta"]),
        "spectrum_requested": entry["spectrum"],
        "arm_requested": entry["arm"],
        # --- sea state / spectrum (run_v15_spectrum.assemble) -------------
        "omega_c": float(ctx["omega_c"]),
        "warn_extrapolated_wave_age": bool(ctx["warn_extrapolated_wave_age"]),
        "warn_omega_c_threshold": float(V15S.WARN_OMEGA_C),
        "sigma2_h_total": float(ctx["sigma2_h_total"]),
        "spectrum_source_default": ctx["spectrum_source_default"],
        "merge_arm": ctx["merge_arm"],
        "v14_arm": ctx["v14_arm"],
        "cell_id": V15S.cell_id(entry["U"], entry["fetch"]),
        "depth_m": float(ctx["depth_m"]),
        "laws": {k: (float(v) if isinstance(v, (int, float)) else v)
                 for k, v in ctx["laws"].items()},
        "eckv_domain": ctx["eckv_domain"],
        "alpha_m_resolution": ctx["alpha_m_resolution"],
        "v15_provenance": ctx["v15_provenance"],
        # --- medium (v12_water.preset_spectra) ----------------------------
        "medium_preset": medium["preset"],
        "medium_label": medium["label"],
        "medium_load": float(medium["load"]),
        "medium_g": float(medium["g"]),
        "medium_cdom_a440": float(medium["cdom_a440"]),
        "a_rgb": [float(x) for x in medium["a"]],
        "b_rgb": [float(x) for x in medium["b"]],
        "c_rgb": [float(x) for x in medium["c"]],
        # --- phase (run_v16_water.FF_PHASE) -------------------------------
        "phase_n": float(phase["n"]),
        "phase_mu": float(phase["mu"]),
        "phase_label": phase["label"],
        "phase_provenance": phase["provenance"],
        "phase_assumption": bool(phase["phase_assumption"]),
        # --- eta provenance (recipe_v9.parse_recipe_v9) -------------------
        "eta_source": sky["source"],
        "eta_recipe": recipe_path,
        "eta_recipe_sha256": _sha256(recipe_path),
        "sky_mode": sky["mode"],
        "sun_only": bool(float(sky["eta"]) == 1.0),
    }
    return out


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    print(json.dumps(resolve_preset_env(5), ensure_ascii=False, indent=2))
