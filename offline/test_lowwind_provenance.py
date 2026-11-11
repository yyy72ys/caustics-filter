"""test_lowwind_provenance -- v17 automated provenance test (NEW FILE).

Brief: BRIEF-V17-LOWWIND-STATUS.md, Stage 2.

What it checks
--------------
A. the settled wording in `lowwind_status` (status prefix, "not a gate pass",
   no "ADOPTED" in the heading or the text, the machine-readable dict, and the
   absence of the rejected phrase in the module source);
B. the four low-wind cells (U in {1.5, 3.0} x F in {100, 1000}) evaluated
   through the SAME entry point the v15 numbers file came from
   (`run_v15_spectrum.evaluate_cell`, depth 1.0 m as in
   `out/v15/numbers-v15-merge.json`):
     B1 the `low_wind_prior` label is present exactly when the cell is inside
        the prior zone (u* < c_m/e), and absent otherwise;
     B2 the provenance strings agree with the branch the cell is on
        (`low_wind_prior` marker inside the zone, `eckv_fit_domain_raw`
        outside it);
     B3 no low-wind cell is an unqualified PASS: every one is UNRESOLVED (or
        FAIL with a hard clause printed) and carries the unresolved label;
     B4 the recomputed record agrees with the stored v15 record (state,
        alpha_m_used, sigma2_h closure, provenance strings);
     B5 the brief's literal premise ("the prior was active in every low-wind
        cell") is reported as a DOCUMENTED DEVIATION where it does not hold
        (at U = 3.0, u* is above the ECKV fit crossing), never silently.

Field names as they actually exist in the v15 family (established by the
Stage 0 probe, not assumed here):
    record["provenance"]["provenance"]        -> arm label, e.g. "low_wind_absolute_gate"
    record["provenance"]["alpha_m_domain"]    -> "low_wind_prior" | "eckv_fit_domain"
    record["alpha"]["prior_active"]           -> bool
    record["alpha"]["low_wind_prior"]         -> label dict or None
    record["alpha"]["alpha_m_domain"]         -> same domain string
    record["provenance"]["low_wind_status"]   -> "UNRESOLVED_AT_BEST (U <= 3: ...)"

Exit code: 0 only if every non-deviation assertion passes.
"""
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))

try:                                   # keep the UTF-8 wording printable on cp932
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import lowwind_status as LS
import run_v15_spectrum as V15

STORED_PATH = os.path.join("out", "v15", "numbers-v15-merge.json")
OUT_PATH = os.path.join("out", "v17", "numbers-v17-lowwind.json")
DEPTH_M = 1.0                          # from the stored v15 records (cell.depth_m)
TOL = 1e-12                            # relative tolerance for the number cross-check

LOW_WIND_CELLS = ((1.5, 100.0), (1.5, 1000.0), (3.0, 100.0), (3.0, 1000.0))

RESULTS = []
DEVIATIONS = []


def check(aid, what, ok, expected=None, got=None, note=None):
    """Record one assertion and print it."""
    ok = bool(ok)
    RESULTS.append({"id": aid, "what": what, "pass": ok,
                    "expected": expected, "got": got, "note": note})
    line = "%-6s %-4s %s" % (aid, "PASS" if ok else "FAIL", what)
    print(line)
    if not ok:
        print("            expected: %r" % (expected,))
        print("            got     : %r" % (got,))
    if note:
        print("            note    : %s" % (note,))
    return ok


def rel_close(a, b, tol=TOL):
    if a is None or b is None:
        return a is b
    a, b = float(a), float(b)
    if a == b:
        return True
    return abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)


def label_paths(rec):
    """Every place the v15 family actually carries provenance / prior evidence."""
    prov = rec.get("provenance", {}) or {}
    alpha = rec.get("alpha", {}) or {}
    found = {
        "provenance.provenance": prov.get("provenance"),
        "provenance.alpha_m_domain": prov.get("alpha_m_domain"),
        "provenance.low_wind_cell": prov.get("low_wind_cell"),
        "provenance.low_wind_status": prov.get("low_wind_status"),
        "alpha.provenance": alpha.get("provenance"),
        "alpha.alpha_m_domain": alpha.get("alpha_m_domain"),
        "alpha.prior_active": alpha.get("prior_active"),
        "alpha.top_up_acted": alpha.get("top_up_acted"),
        "alpha.low_wind_prior": alpha.get("low_wind_prior"),
        "alpha.sigma2_h_target": alpha.get("sigma2_h_target"),
    }
    flat = json.dumps(found, ensure_ascii=False)
    found["_string_low_wind_prior_present_in"] = [
        k for k, v in found.items() if isinstance(v, str) and "low_wind_prior" in v]
    found["_string_low_wind_prior_anywhere"] = "low_wind_prior" in flat
    return found


def main():
    print("=" * 78)
    print("v17 low-wind provenance test")
    print("brief      : BRIEF-V17-LOWWIND-STATUS.md (Stage 2)")
    print("module     : lowwind_status (status_id=%s)" % LS.STATUS_ID)
    print("spectrum   : run_v15_spectrum.evaluate_cell (the v15 merge entry point)")
    print("stored ref : %s" % STORED_PATH)
    print("=" * 78)

    # ---------------- A. the settled wording ------------------------------
    print("\n-- A. settled wording -------------------------------------------")
    check("A1", 'STATUS_TEXT starts with "Status: UNRESOLVED"',
          LS.STATUS_TEXT.startswith("Status: UNRESOLVED"),
          "Status: UNRESOLVED...", LS.STATUS_TEXT[:40])
    check("A2", 'STATUS_TEXT contains "not a gate pass"',
          "not a gate pass" in LS.STATUS_TEXT,
          True, "not a gate pass" in LS.STATUS_TEXT)
    check("A3", '"ADOPTED" absent from STATUS_TEXT and STATUS_HEADING',
          ("ADOPTED" not in LS.STATUS_TEXT) and ("ADOPTED" not in LS.STATUS_HEADING),
          True,
          {"in_text": "ADOPTED" in LS.STATUS_TEXT,
           "in_heading": "ADOPTED" in LS.STATUS_HEADING})
    check("A4", 'heading leads with UNRESOLVED and says "not a gate pass"',
          LS.STATUS_HEADING.startswith("UNRESOLVED")
          and "not a gate pass" in LS.STATUS_HEADING,
          "UNRESOLVED ... (not validated; not a gate pass)", LS.STATUS_HEADING)

    d = LS.STATUS_DICT
    dict_ok = (d.get("status_id") == "UNRESOLVED_OPERATIONAL_PRIOR_RETAINED"
               and d.get("sigma2_h") == 0.003
               and list(d.get("sensitivity_interval", [])) == [0.0015, 0.006]
               and d.get("interval_kind") == "declared sensitivity interval (not an evidence bracket)"
               and d.get("provenance_label") == "low_wind_prior"
               and str(d.get("gate_C1_low_wind", "")).startswith("NOT PASSED")
               and "k > 200 rad/m" in str(d.get("supersession", "")))
    check("A5", "STATUS_DICT machine-readable fields", dict_ok,
          "status_id/sigma2_h/interval/interval_kind/label/gate_C1/supersession",
          d)

    blk = LS.status_block()
    block_ok = (all(k in blk for k in d)
                and blk.get("heading") == LS.STATUS_HEADING
                and blk.get("text") == LS.STATUS_TEXT)
    check("A6", "status_block() returns STATUS_DICT + heading + text", block_ok,
          "superset of STATUS_DICT with heading/text", sorted(blk.keys()))

    with open(LS.__file__, "r", encoding="utf-8") as fh:
        src = fh.read()
    check("A7", 'rejected phrase absent from the module source',
          "ADOPTED_WITH_BASIS" not in src,
          "0 occurrences", src.count("ADOPTED_WITH_BASIS"))

    # ---------------- load the stored v15 records -------------------------
    with open(STORED_PATH, "r", encoding="utf-8") as fh:
        stored_doc = json.load(fh)
    stored = {c["cell_id"]: c for c in stored_doc["cells"]}

    # ---------------- B. the four low-wind cells --------------------------
    cells_out = []
    for U, F in LOW_WIND_CELLS:
        cid = V15.cell_id(U, F)
        print("\n-- B. cell %s (U = %.1f m/s, F = %.0f m, h = %.1f m) --------"
              % (cid, U, F, DEPTH_M))
        rec = V15.evaluate_cell(float(U), float(F), DEPTH_M)
        st = stored[cid]
        found = label_paths(rec)

        ustar = float(rec["alpha"]["ustar_ms"])
        fit_min = float(rec["alpha"]["ustar_fit_min_ms"])
        in_domain = bool(rec["alpha"]["in_eckv_fit_domain"])
        expect_prior = not in_domain     # the prior zone is u* < c_m/e
        prior_active = bool(rec["alpha"]["prior_active"])
        label = rec["alpha"]["low_wind_prior"]

        # B1 label presence follows the branch
        b1 = (prior_active == expect_prior)
        if expect_prior:
            b1 = b1 and isinstance(label, dict) and label.get("sigma2_h_target") == 0.003
        else:
            b1 = b1 and label is None
        check("B1", "low_wind_prior label present exactly where the prior acts",
              b1,
              {"expect_prior": expect_prior, "label_sigma2_h": 0.003 if expect_prior else None},
              {"prior_active": prior_active, "low_wind_prior": label})

        # B2 provenance strings agree with the branch
        prov_arm = rec["provenance"]["provenance"]
        dom_prov = rec["provenance"]["alpha_m_domain"]
        dom_alpha = rec["alpha"]["alpha_m_domain"]
        if expect_prior:
            b2 = (dom_prov == "low_wind_prior" and dom_alpha == "low_wind_prior"
                  # the arm label the v15 family emits for the merge arm D
                  and prov_arm in ("low_wind_prior", "low_wind_absolute_gate"))
        else:
            b2 = (dom_prov == "eckv_fit_domain" and dom_alpha == "eckv_fit_domain"
                  and prov_arm == "eckv_fit_domain_raw")
        check("B2", "provenance strings match the cell's branch", b2,
              ("alpha_m_domain=low_wind_prior" if expect_prior
               else "alpha_m_domain=eckv_fit_domain, provenance=eckv_fit_domain_raw"),
              {"provenance.provenance": prov_arm,
               "provenance.alpha_m_domain": dom_prov,
               "alpha.alpha_m_domain": dom_alpha})

        # B3 never an unqualified PASS
        state = rec["state"]
        lw_status = str(rec["provenance"]["low_wind_status"])
        b3 = (state != "PASS"
              and bool(rec["low_wind_cell"]) is True
              and state in ("UNRESOLVED", "FAIL")
              and lw_status.startswith("UNRESOLVED_AT_BEST")
              and (state != "FAIL" or bool(rec["hard_clauses_fired"])))
        check("B3", "cell is not an unqualified PASS (UNRESOLVED/labelled)",
              b3,
              {"state": "UNRESOLVED|FAIL(hard clauses printed)",
               "low_wind_status": "UNRESOLVED_AT_BEST..."},
              {"state": state, "low_wind_cell": rec["low_wind_cell"],
               "low_wind_status": lw_status,
               "hard_clauses_fired": rec["hard_clauses_fired"]})

        # B4 agreement with the stored v15 record
        b4 = (state == st["state"]
              and rel_close(rec["alpha"]["alpha_m_used"], st["alpha"]["alpha_m_used"])
              and rel_close(rec["alpha"]["sigma2_h_total_closure"],
                            st["alpha"]["sigma2_h_total_closure"])
              and prov_arm == st["provenance"]["provenance"]
              and dom_prov == st["provenance"]["alpha_m_domain"])
        check("B4", "recomputed record agrees with out/v15/numbers-v15-merge.json",
              b4,
              {"state": st["state"],
               "alpha_m_used": st["alpha"]["alpha_m_used"],
               "sigma2_h_closure": st["alpha"]["sigma2_h_total_closure"]},
              {"state": state,
               "alpha_m_used": rec["alpha"]["alpha_m_used"],
               "sigma2_h_closure": rec["alpha"]["sigma2_h_total_closure"]})

        # B5 the brief's literal premise, as a documented deviation
        brief_premise = prior_active
        if not prior_active:
            dev = {
                "cell_id": cid,
                "brief_premise": "the low-wind prior was active in every low-wind cell",
                "measured": ("in_eckv_fit_domain = true: u* = %.6f m/s >= c_m/e = "
                             "%.6f m/s, so alpha_m is the raw ECKV value (%.6e) and no "
                             "floor is applied; the cell is still low-wind (U <= 3) and "
                             "still UNRESOLVED" % (ustar, fit_min,
                                                   rec["alpha"]["alpha_m_raw"])),
                "u_over_fit_min": rec["alpha"]["u_over_fit_min"],
                "documented": True,
            }
            DEVIATIONS.append(dev)
            print("DEVIATION  B5  %s: prior NOT active (u* = %.6f >= c_m/e = %.6f); "
                  "documented, not silently passed" % (cid, ustar, fit_min))
        else:
            print("INFO   B5  %s: prior active as the brief expected" % cid)

        cells_out.append({
            "cell_id": cid,
            "U_ms": float(U),
            "fetch_m": float(F),
            "depth_m": DEPTH_M,
            "state": state,
            "state_stored": st["state"],
            "low_wind_cell": bool(rec["low_wind_cell"]),
            "low_wind_status": lw_status,
            "ustar_ms": ustar,
            "ustar_fit_min_ms": fit_min,
            "u_over_fit_min": rec["alpha"]["u_over_fit_min"],
            "in_eckv_fit_domain": in_domain,
            "prior_expected_by_branch": expect_prior,
            "prior_active": prior_active,
            "prior_top_up_acted": bool(rec["alpha"]["top_up_acted"]),
            "alpha_m_raw": rec["alpha"]["alpha_m_raw"],
            "alpha_m_floor": rec["alpha"]["alpha_m_floor"],
            "alpha_m_used": rec["alpha"]["alpha_m_used"],
            "sigma2_h_target": rec["alpha"]["sigma2_h_target"],
            "sigma2_h_total_closure": rec["alpha"]["sigma2_h_total_closure"],
            "provenance_keys_found": found,
            "brief_literal_premise_holds": bool(brief_premise),
            "assertions": RESULTS[-4:],
        })

    # ---------------- verdict --------------------------------------------
    n_pass = sum(1 for r in RESULTS if r["pass"])
    n_fail = len(RESULTS) - n_pass
    verdict = "PASS" if n_fail == 0 else "FAIL"
    print("\n" + "=" * 78)
    print("assertions: %d PASS, %d FAIL, %d total" % (n_pass, n_fail, len(RESULTS)))
    print("documented deviations from the brief's literal premise: %d of 4 low-wind cells"
          % len(DEVIATIONS))
    for dv in DEVIATIONS:
        print("  - %s : %s" % (dv["cell_id"], dv["measured"]))
    print("VERDICT: %s" % verdict)
    print("=" * 78)

    out = {
        "stage": "stage2-lowwind-provenance-test",
        "brief": "BRIEF-V17-LOWWIND-STATUS.md",
        "module_under_test": "lowwind_status",
        "status_id": LS.STATUS_ID,
        "status_heading": LS.STATUS_HEADING,
        "status_text_sha256": hashlib.sha256(
            LS.STATUS_TEXT.encode("utf-8")).hexdigest(),
        "status_text_chars": len(LS.STATUS_TEXT),
        "spectrum_module": "run_v15_spectrum",
        "spectrum_entry_point": "evaluate_cell(U, F, h)",
        "depth_m": DEPTH_M,
        "stored_reference": STORED_PATH,
        "stored_reference_cell_states": {k: v["state"] for k, v in stored.items()},
        "provenance_field_names_actually_used": {
            "provenance.provenance": "arm label (low_wind_absolute_gate | eckv_fit_domain_raw)",
            "provenance.alpha_m_domain": "low_wind_prior | eckv_fit_domain",
            "alpha.prior_active": "bool",
            "alpha.low_wind_prior": "label dict or None",
            "alpha.alpha_m_domain": "low_wind_prior | eckv_fit_domain",
            "note": ("the brief expects a `low_wind_prior` PROVENANCE string; the v15 "
                     "family puts that string on alpha_m_domain and puts the arm label "
                     "on provenance.provenance"),
        },
        "cells": cells_out,
        "assertions": RESULTS,
        "assertion_summary": {"n_pass": n_pass, "n_fail": n_fail,
                              "n_total": len(RESULTS)},
        "documented_deviations": DEVIATIONS,
        "verdict": verdict,
        "exit_code": 0 if n_fail == 0 else 1,
    }
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("wrote %s" % OUT_PATH)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
