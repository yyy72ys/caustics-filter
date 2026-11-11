"""caustics-recipe v2 -> v7 observer block (report.md section 4).

One optional block is added on top of the v6 parser:

    observer.enable            bool    (default False = OFF, v6 behaviour)
    observer.peak_cd_m2        number > 0  (default 200: display peak white)
    observer.La_cd_m2          number > 0 or null (default null = image mean)
    observer.age_yr            number  (default 28.6; Watson-Yellott reference
                                       readout only, no effect on the MS node)
    observer.defocus_diopters  number >= 0 (default 0.0 = in focus)
    observer.pixels_per_degree number > 0  (default 60; viewing geometry,
                                       prerequisite of blur/CSF, report 4.3/4.4)
    observer.csf_strength      0..1    (default 0.0 = analysis only, no shaping)

Absent block (or enable=false) reads EXACTLY as v6 ("observer OFF" is a hard
bypass: the still is returned untouched, bit-identical -- gate 1).  A v1
recipe without any new field reads as before.  Existing fields are not
changed; parse_recipe_v6 does the v6 half of the job and is called unchanged.
"""

from recipe_v6 import parse_recipe_v6

OBSERVER_DEFAULTS = {
    "enable": False,
    "peak_cd_m2": 200.0,
    "La_cd_m2": None,
    "age_yr": 28.6,
    "defocus_diopters": 0.0,
    "pixels_per_degree": 60.0,
    "csf_strength": 0.0,
}


def _positive(name, value):
    v = float(value)
    if not (v == v) or v in (float("inf"), float("-inf")) or not (v > 0.0):
        raise SystemExit("%s must be a finite number > 0 (got %r)"
                         % (name, value))
    return v


def parse_recipe_v7(recipe):
    """v6 parse (unchanged) + the observer block."""
    cfg = parse_recipe_v6(recipe)
    obs = recipe.get("observer") or {}
    if not isinstance(obs, dict):
        raise SystemExit("observer must be an object (got %r)" % (obs,))

    enable = bool(obs.get("enable", False))
    peak = _positive("observer.peak_cd_m2",
                     obs.get("peak_cd_m2", OBSERVER_DEFAULTS["peak_cd_m2"]))
    La = obs.get("La_cd_m2", None)
    if La is not None:
        La = _positive("observer.La_cd_m2", La)
    try:
        age = float(obs.get("age_yr", OBSERVER_DEFAULTS["age_yr"]))
    except Exception:
        raise SystemExit("observer.age_yr must be a number (got %r)"
                         % (obs.get("age_yr"),))
    if not (age == age) or age in (float("inf"), float("-inf")):
        raise SystemExit("observer.age_yr must be finite (got %r)"
                         % (obs.get("age_yr"),))
    try:
        dD = float(obs.get("defocus_diopters",
                           OBSERVER_DEFAULTS["defocus_diopters"]))
    except Exception:
        raise SystemExit("observer.defocus_diopters must be a number (got %r)"
                         % (obs.get("defocus_diopters"),))
    if not (dD == dD) or not (dD >= 0.0) or dD in (float("inf"),):
        raise SystemExit("observer.defocus_diopters must be >= 0 (got %r)"
                         % (obs.get("defocus_diopters"),))
    ppd = _positive("observer.pixels_per_degree",
                    obs.get("pixels_per_degree",
                            OBSERVER_DEFAULTS["pixels_per_degree"]))
    try:
        csf = float(obs.get("csf_strength",
                            OBSERVER_DEFAULTS["csf_strength"]))
    except Exception:
        raise SystemExit("observer.csf_strength must be a number (got %r)"
                         % (obs.get("csf_strength"),))
    if not (0.0 <= csf <= 1.0) or not (csf == csf):
        raise SystemExit("observer.csf_strength must be in [0,1] (got %r)"
                         % (obs.get("csf_strength"),))

    cfg = dict(cfg)
    cfg["observer"] = {
        "enable": enable,
        "peak_cd_m2": float(peak),
        "La_cd_m2": (None if La is None else float(La)),
        "age_yr": float(age),
        "defocus_diopters": float(dD),
        "pixels_per_degree": float(ppd),
        "csf_strength": float(csf),
        "legacy": (not enable and not obs),
        "note": ("enable=false (default) is a hard bypass: the composite "
                 "still is returned untouched, bit-identical to the v6 "
                 "still (gate 1).  La=null means the image mean displayed "
                 "luminance (stated steady-state choice)."),
    }
    cfg["legacy_recipe"] = bool(cfg.get("legacy_recipe")
                                and not obs and not enable)
    cfg["schema_v7_fields"] = ["observer.enable", "observer.peak_cd_m2",
                               "observer.La_cd_m2", "observer.age_yr",
                               "observer.defocus_diopters",
                               "observer.pixels_per_degree",
                               "observer.csf_strength"]
    return cfg


def with_observer(recipe, enable=True, peak_cd_m2=200.0, La_cd_m2=None,
                  defocus_diopters=0.0, pixels_per_degree=60.0,
                  csf_strength=0.0, age_yr=28.6):
    """Deep copy of a recipe with an explicit observer block."""
    import copy
    r = copy.deepcopy(recipe)
    r["observer"] = {"enable": bool(enable),
                     "peak_cd_m2": float(peak_cd_m2),
                     "La_cd_m2": La_cd_m2,
                     "age_yr": float(age_yr),
                     "defocus_diopters": float(defocus_diopters),
                     "pixels_per_degree": float(pixels_per_degree),
                     "csf_strength": float(csf_strength)}
    return r
