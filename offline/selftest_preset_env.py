"""selftest_preset_env -- the 11 pre-registered checks for preset id 5.

Brief: BRIEF-PRESET-TABLE.md, section "自己テスト（検証値）".  The expected
values are the ones BRIEF-PRESET5 measured; they are NOT derived from this
module.  A mismatch is printed and reported as FAIL -- nothing is adjusted to
make a check pass.

Run:  python -B selftest_preset_env.py
Exit: 0 all PASS, 1 any FAIL.
"""

import sys

import presets_env

#: (label, key, expected, tolerance).  tolerance None = exact equality.
CHECKS = (
    ("1. U", "U", 6.0, None),
    ("2. fetch", "fetch", 100.0, None),
    ("3. z", "z", 4.0, None),
    ("4. medium_preset", "medium_preset", "clear", None),
    ("5. phase_n", "phase_n", 1.098, None),
    ("6. phase_mu", "phase_mu", 3.37, None),
    ("7. eta", "eta", 1.0, None),
    ("8. spectrum", "spectrum", "elfouhaily-form", None),
    ("9. arm", "arm", "D", None),
    ("10. omega_c", "omega_c", 7.296464998547599, 1e-6),
    ("11. warn_extrapolated_wave_age", "warn_extrapolated_wave_age", True, None),
)


def main():
    print("selftest_preset_env -- resolve_preset_env(5)")
    print("registered ids: %s" % (presets_env.preset_env_ids(),))
    entry = presets_env.get_preset_env(5)
    print("table entry   : %r" % (entry,))
    print("")

    res = presets_env.resolve_preset_env(5)

    print("%-34s %-10s %-26s %-26s %s"
          % ("check", "result", "measured", "expected", "delta"))
    print("-" * 108)

    fails = []
    for label, key, expected, tol in CHECKS:
        if key not in res:
            got = "<MISSING KEY>"
            ok = False
            delta = "n/a"
        else:
            got = res[key]
            if tol is None:
                ok = (got == expected) and (type(got) is type(expected))
                delta = "exact"
            else:
                ok = abs(float(got) - float(expected)) <= float(tol)
                delta = "%.3e (tol %g)" % (abs(float(got) - float(expected)),
                                           float(tol))
        if not ok:
            fails.append((label, got, expected))
        print("%-34s %-10s %-26r %-26r %s"
              % (label, "PASS" if ok else "FAIL", got, expected, delta))

    print("-" * 108)
    print("PASS %d / %d" % (len(CHECKS) - len(fails), len(CHECKS)))

    if fails:
        print("")
        print("MISMATCHES (reported, not adjusted):")
        for label, got, expected in fails:
            print("  %s: measured %r, expected %r" % (label, got, expected))
        return 1

    # context lines, printed only so the report can quote the provenance
    print("")
    print("context (not part of the 11):")
    for k in ("name", "medium", "medium_label", "medium_load", "phase_key",
              "phase_label", "eta_source", "eta_recipe_sha256", "cell_id",
              "sigma2_h_total", "eckv_domain", "warn_omega_c_threshold",
              "spectrum_source_default", "merge_arm", "v14_arm",
              "phase_assumption"):
        print("  %-24s %r" % (k, res[k]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
