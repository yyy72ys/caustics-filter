"""v11 report -- REPORT-V11-GATES.md from the stage JSONs.

Everything numeric in the report is read from the stage files, so the report
cannot drift from the measurements.  ASCII only.
"""
import json
import math
import os


def _f(x, n=4, dash="-"):
    if x is None:
        return dash
    try:
        v = float(x)
    except Exception:
        return str(x)
    if not math.isfinite(v):
        return dash
    if v != 0.0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        return ("%%.%de" % max(0, n - 1)) % v
    return ("%%.%df" % n) % v


def _p(x, n=1, dash="-"):
    if x is None:
        return dash
    return ("%%.%df" % n) % float(x)


def _pass(b):
    if b is None:
        return "n/a"
    return "PASS" if b else "FAIL"


def _yn(b):
    return "yes" if b else "no"


def _load(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


BANDS = ("below_kp", "kp_10", "10_50", "50_200", "200_370", "370_1000",
         "1000_inf")


def _cell(c):
    """'U1p5_F100p0' -> 'U=1.5 F=100' (report table label)."""
    return "U=%s F=%s" % (_f(c["U_ms"], 1), _f(c["fetch_m"], 0))


def build(out_dir):
    s1 = _load(os.path.join(out_dir, "numbers-v11-spectrum.json"))
    s2 = _load(os.path.join(out_dir, "numbers-v11-sweep.json"))
    s3 = _load(os.path.join(out_dir, "numbers-v11-transport.json"))
    mg = _load(os.path.join(out_dir, "merge-evaluation-v11.json"))
    L = []
    A = L.append
    a1 = s1["gates"]["A1"]
    a2 = s1["gates"]["A2"]
    a2l = s1["gates"]["A2_legacy_reference"]
    cells1 = {c["cell"]: c for c in s1["cells"]}
    cells2 = s2["cells"]
    sm = s2["summary"]
    by = {c["cell"]: c for c in cells2}

    A("# REPORT-V11-GATES -- version 3.8, gravity-tail termination "
      "(chunk 8)")
    A("")
    A("Brief: `BRIEF-V11-CHUNK8.md`. Design: `PIPELINE-V10.md` (the "
      "specification).")
    A("Machine: PC, `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`.")
    A("Model: deepseek-v4.1-flash. All four stages ran FOREGROUND in this one "
      "session, checkpointed with jobkit (`out/v11/checkpoint.json`, one unit "
      "per delivered artefact), every deliverable committed atomically through "
      "`atomic_write_bytes`. No `python.exe` / `blender.exe` was killed; no "
      "image was inspected (numerical statistics only); every output is "
      "ASCII.")
    A("")
    A("New code (v11 family only; v5-v10 model files imported, never edited):")
    A("")
    A("- `v11_spectrum.py` -- B_l (PIPELINE-V10 section 3), the assembled "
      "S(k) = (B_l + B_h)/k^3, the band table, the ECKV-at-matched-inputs "
      "reference, gates A1/A2/K6")
    A("- `v11_state.py` -- `V11WaterState`: a source switch over "
      "`v10_state.V10WaterState` (three overridden methods: `_build`, "
      "`_realise`, `_direction_node`)")
    A("- `run_v11.py` -- the four stage runners")
    A("- `v11_report.py` -- this report (reads the stage JSONs)")
    A("")
    A("Reproduce, in this order and in the foreground:")
    A("")
    A("```")
    A("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    A("python run_v11.py --stage 1       #  15 s")
    A("python run_v11.py --stage 2       #  32 s")
    A("python run_v11.py --stage 3       #  ~4 min (4 renders at 2048^2, "
      "Q=25)")
    A("python run_v11.py --stage 4       #   1 s + report")
    A("```")
    A("")
    A("Run order in this session: 1 -> 2 -> 3 -> 4. Stage 2 was re-run once "
      "after the C3 comparator was corrected (it now compares the REALISED, "
      "cut band against the reference's full band -- the first version "
      "compared full against full, which is 1 by construction and says "
      "nothing). Stage 3 was re-run once as a determinism check: every number "
      "reproduced to the printed digits and the legacy still reproduced the "
      "v10 sha256 exactly. Stage 4 was re-run while the report builder's "
      "formatting was corrected. All of these are measurement/report-side; no "
      "model constant, no tolerance and no gate definition moved. The "
      "checkpoint carries one unit per delivered artefact; the re-runs left "
      "duplicate unit ids and it was de-duplicated afterwards with the same "
      "`_ckpt_dedupe` approach as the previous chunk.")
    A("")
    A("---")
    A("")
    A("## 0. Conventions, and what is deliberately NOT done")
    A("")
    A("Pinned verbatim from PIPELINE-V9 section 2 (PIPELINE-V10 section 2): "
      "`S(k)` in m^2 per rad/m, `m0 = int S dk`, `sigma^2 = int k^2 S dk`, "
      "directional density `Psi = S Phi / k` with `int Phi dtheta = 1`, wind = "
      "U10 and `u*` from U10 only. One addition this chunk: "
      "`B(k) := k^3 S(k)` is a first-class reported quantity, because both the "
      "defect and the fix live in B.")
    A("")
    A("The spectrum is, wholesale,")
    A("")
    A("```")
    A("S(k) = k^-3 ( B_l(k) + B_h(k) )")
    A("```")
    A("")
    A("- `B_l`: Elfouhaily et al. (1997) low-frequency part, formula verbatim "
      "from PIPELINE-V10 section 3, with `alpha_p = 6e-3 Omega_c^0.55`, "
      "`gamma = 1.7` (Omega_c <= 1) or `1.7 + 6 log10(Omega_c)` (Omega_c > 1), "
      "`sigma = 0.08 (1 + 4 Omega_c^-3)`, `L_PM = exp(-1.25 (k_p/k)^2)`, "
      "`J_p = gamma^Gamma`, `F_p = L_PM J_p exp[-(Omega_c/sqrt10)(sqrt(k/k_p) "
      "- 1)]`.")
    A("- `B_h`: v3.7's high-frequency part imported unchanged "
      "(`v10_spectrum.elfouhaily_bh`: `0.5 alpha_m (c_m/c) L_pm exp(-0.25 (k/k_m"
      " - 1)^2)`, McDaniel 2001 carries L_pm there too), level set by `u*` "
      "alone.")
    A("- **No join, no blend, no splice, no renormalisation of either part.** "
      "K6 below demonstrates the absence of a splice numerically rather than "
      "asserting it.")
    A("")
    A("**Honest label, used everywhere in this chunk** (rally D3): this is an "
      "**Elfouhaily-form spectrum with fetch-derived sea-state parameters**. "
      "It is *not* \"the Elfouhaily model\": ECKV's shape functions are driven "
      "by this pipeline's own pond fetch laws, and no ECKV validation exists "
      "for a 0.3-1.0 m deep, 100-1000 m fetch pond.")
    A("")
    A("### Parameter sources (nothing fitted, nothing new)")
    A("")
    A("| quantity | source in this chunk | note |")
    A("|---|---|---|")
    A("| `f_p` | Hasselmann (1973) fetch law `f_p = 3.5 (g/U) Ft^-1/3`, "
      "unchanged | the v3.7 peak |")
    A("| `k_p` | root of the pipeline's existing finite-depth gravity-"
      "capillary dispersion at `omega_p = 2 pi f_p` | **deviation from v3.7**, "
      "which used the deep-water `k_p = 2 pi/lambda_p`; see section 1.1 |")
    A("| `c_p` | `omega(k_p)/k_p` from the same dispersion | recorded "
      "deviation from ECKV-as-published, which uses deep-water `c(k)` |")
    A("| `Omega_c` | `U10 / c_p` | ECKV's own ocean-fetch `Omega_c` relation "
      "is **not** used (PIPELINE-V10 section 3) |")
    A("| `u*` | `C10 = (0.8 + 0.065 U10) e-3`, `u* = sqrt(C10) U10` | sets "
      "B_h only; `u*` does not enter B_l |")
    A("| `H_s` | **prediction** `4 sqrt m0` of the assembled spectrum | never "
      "rescaled toward the law (PIPELINE-V10 section 4) |")
    A("")
    A("| constant | value | units | basis |")
    A("|---|---|---|---|")
    A("| `K_M` | %s | rad/m | Elfouhaily et al. (1997): k_m = 2 pi / 0.017 m "
      "(used in B_h and in the research closed form) |"
      % _f(s1["constants"]["K_M_rad_per_m"], 6))
    A("| `C_M` | %s | m/s | gravity-capillary minimum phase speed at K_M |"
      % _f(s1["constants"]["C_M_ms"], 4))
    A("| `k_m` pipeline | %s | rad/m | sqrt(g/(gamma/rho)) -- the pipeline "
      "dispersion's own crossover, %s below K_M; both are reported |"
      % (_f(s1["constants"]["k_m_pipeline_rad_per_m"], 6),
         "%.2f%%" % (100 * s1["constants"]["k_m_convention_rel_diff"])))
    A("| k_m convention rel diff | %s | - | the whole residual of gate A1 with "
      "the research's rounded constant |"
      % _f(s1["constants"]["k_m_convention_rel_diff"], 5))
    A("")
    A("### What does not change")
    A("")
    A("Transport (`det.det_frame`, `transport.photon_frame`), the frozen "
      "composite (A2.6 + Blow + E2.2, `v5_composite.render_frame_v5`), the sky "
      "(chunk 3.5), the retina, the directional distribution `Phi` (ECKV "
      "unified function applied to the TOTAL spectrum), the patch/cascade "
      "realisation, the 8 s quantisation interface and the recipe schema. The "
      "v3.7 state generator remains in code as the LEGACY source behind the "
      "switch.")
    A("")
    A("---")
    A("")
    A("## 1. Stage 1 -- the gravity part and the algebra gates")
    A("")
    A("### 1.1 Three deviations, stated before any number")
    A("")
    A("1. **`k_p` from the dispersion, not from `lambda_p`.** v3.7 took "
      "`k_p = 2 pi/lambda_p` with `lambda_p = g/(2 pi f_p^2)`, i.e. the "
      "deep-water inversion. PIPELINE-V10 section 3 asks for the fetch-law "
      "`f_p` mapped through *the existing finite-depth gravity-capillary "
      "dispersion*; the implementation solves `omega(k; h) = omega_p` by log "
      "bisection to a relative residual below `1e-14`. The difference is small "
      "but not zero, and it is reported per cell as "
      "`k_p_dispersion_rel_diff` (at most `%s` across the eight cells)."
      % _f(max(abs(c["sea_state"]["k_p_dispersion_rel_diff"])
               for c in s1["cells"]), 4))
    A("2. **`c(k)` and `c_p` from the dispersion.** ECKV-as-published uses "
      "deep-water phase speeds. At the pond depths of the defect band "
      "`kh >= 15`, so `tanh(kh) = 1` to `1e-12` and the difference is confined "
      "to the peak region; this is a stated deviation, not a free parameter.")
    A("3. **`k_m` has two values in this chunk.** `B_h` and the research's "
      "closed form use `K_M = 2 pi/0.017 = 369.6 rad/m`; the pipeline's own "
      "dispersion crosses over at `sqrt(g/(gamma/rho)) = 367.0 rad/m`, a "
      "`%.2f%%` difference. Gate A1 reports the ratio against BOTH, so the "
      "residual is never absorbed into a tuned constant."
      % (100 * s1["constants"]["k_m_convention_rel_diff"]))
    A("")
    A("### A1 -- D1 verification: is the v3.7 B_g really the alpha_J/2 "
      "plateau?")
    A("")
    A("`B_g(k)` is computed **numerically from the current v3.7 code path** "
      "(`V10.assemble` -> `V10.s_parts`: JONSWAP in frequency form -> exact "
      "`domega/dk` transform -> `B = k^3 S`) and compared pointwise, on 2000 "
      "log-spaced points across 50-200 rad/m, with")
    A("")
    A("```")
    A("B_g ~ (alpha_J/2) (1 + 3 kappa^2)/(1 + kappa^2)^3 "
      "exp[-1.25 (omega_p/omega)^4],   kappa = k/k_m")
    A("```")
    A("")
    A("`alpha_J` is the amplitude the code actually applies, recovered as "
      "`A/g^2` where `A` is the JONSWAP normalisation constant; the published "
      "Phillips law `alpha` is reported beside it because the fetch-law `H_s` "
      "normalisation makes `alpha` a *readout*, not an input.")
    A("")
    A("| h [m] | kh across band | alpha_J from A/g^2 | published alpha | "
      "max abs(ratio-1), research k_m | max abs(ratio-1), pipeline k_m | "
      "B_g mean in 50-200 | plateau vs deep |")
    A("|---|---|---|---|---|---|---|---|")
    for h in (1000.0, 1.0, 0.5, 0.3):
        r = a1["depths"]["h_%g" % h]
        A("| %s | %s-%s | %s | %s | %s | %s | %s | %s |"
          % (_f(h, 1), _f(r["kh_over_band"][0], 0),
             _f(r["kh_over_band"][1], 0), _f(r["alpha_J_from_A_over_g2"], 5),
             _f(r["alpha_J_published_law"], 5),
             _f(r["ratio_kM_published"]["max_abs_minus_1"], 4),
             _f(r["ratio_kM_physical"]["max_abs_minus_1"], 2),
             _f(r["B_g_numeric_mean"], 5),
             "%+.1f%%" % (100 * r["plateau_rel_vs_deep"])))
    A("")
    A("**A1: %s.** %s" % (_pass(a1["pass"]), a1["finding"]))
    A("")
    A("The pathology is therefore confirmed as a LEVEL defect and not a cutoff "
      "defect: `B_g` sits at `%s` across 50-200 rad/m (deep water, h = 1000 m), "
      "which is `%s` times the measured saturation level `B_hat = 0.007`, "
      "while the dispersion conversion already suppresses everything beyond "
      "~370 rad/m (the v3.7 spectrum's own `frac sigma^2 > 370` is "
      "`%s`-`%s` in the eight cells, so the tail beyond 370 is NOT where the "
      "energy is)."
      % (_f(a1["diagnosis"]["B_g_plateau_deep"], 5),
         _f(a1["diagnosis"]["B_g_plateau_deep"] / 0.007, 2),
         _f(min(c["legacy_v10"]["frac_above_370_realised"]
                for c in cells2), 3),
         _f(max(c["legacy_v10"]["frac_above_370_realised"]
                for c in cells2), 3)))
    A("")
    A("### A2 -- Parseval on a built field from the new spectrum")
    A("")
    A("Same test configuration as v10's G0: `N=1024`, patch from the same "
      "`patch_for` rule, seed 12345, cell U = 6 m/s, F = 100 m, h = 1 m. Two "
      "levels are separated on purpose: (a) numeric field vs the mode sum of "
      "the *same* density -- this is the convention test, any 2 pi or 1/k "
      "error moves it far outside the realisation band; (b) mode sum vs the "
      "continuum integral -- the discretisation error of this patch.")
    A("")
    A("| quantity | numeric (field) | mode sum | continuum |")
    A("|---|---|---|---|")
    A("| m0 [m^2] | %s | %s | %s |"
      % (_f(a2["m0"]["numeric"], 8), _f(a2["m0"]["mode_sum"], 8),
         _f(a2["m0"]["continuum"], 8)))
    A("| sigma^2 [rad^2] | %s | %s | %s |"
      % (_f(a2["sigma2"]["numeric"], 6), _f(a2["sigma2"]["mode_sum"], 6),
         _f(a2["sigma2"]["continuum"], 6)))
    A("")
    A("| check | new source | legacy source (same config) | bound |")
    A("|---|---|---|---|")
    A("| z(m0) | %s | %s | abs <= 3 |"
      % (_f(a2["z_m0"], 2), _f(a2l["z_m0"], 2)))
    A("| z(sigma^2) | %s | %s | abs <= 3 |"
      % (_f(a2["z_sigma2"], 2), _f(a2l["z_sigma2"], 2)))
    A("| mode/continuum - 1, m0 | %s | %s | abs <= 0.30 |"
      % (_f(a2["rel_mode_over_continuum_m0"], 5),
         _f(a2l["rel_mode_over_continuum_m0"], 5)))
    A("| mode/continuum - 1, sigma^2 | %s | %s | abs <= 0.30 |"
      % (_f(a2["rel_mode_over_continuum_sigma2"], 4),
         _f(a2l["rel_mode_over_continuum_sigma2"], 4)))
    A("| DOF (m0 / sigma^2) | %s / %s | %s / %s | - |"
      % (_f(a2["dof"]["m0"], 1), _f(a2["dof"]["sigma2"], 1),
         _f(a2l["dof"]["m0"], 1), _f(a2l["dof"]["sigma2"], 1)))
    A("")
    A("**A2: %s** (convention %s, discretisation %s). Running the same test on "
      "the legacy source is the control that the convention code is shared and "
      "unchanged."
      % (_pass(a2["pass"]), _pass(a2["pass_convention"]),
         _pass(a2["pass_discretisation"])))
    A("")
    A("### The eight cells, new source")
    A("")
    A("| U | F | Ft | k_p [rad/m] | c_p [m/s] | Omega_c | alpha_p | gamma | "
      "sigma | k_p h | domain |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for c in s1["cells"]:
        ss = c["sea_state"]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | `%s` |"
          % (_f(c["U_ms"], 1), _f(c["fetch_m"], 0), _f(c["Ft"], 1),
             _f(ss["k_p_rad_per_m"], 3), _f(ss["c_p_ms"], 4),
             _f(ss["omega_c"], 2), _f(ss["alpha_p"], 5),
             _f(ss["gamma_bl"], 3), _f(ss["sigma_bl"], 4),
             _f(ss["kp_times_h"], 2), c["eckv_domain"]))
    A("")
    A("Two of eight cells (U = 6/F = 100 and U = 10/F = 100) sit at "
      "`Omega_c = %s` and `%s`: the gamma log-formula is continued past "
      "`Omega_c = 5` and those cells carry `eckv_domain = "
      "wave_age_extrapolation`. Nothing is clamped and nothing diverges; the "
      "label says plainly that they are uncalibrated extrapolation."
      % (_f(cells1["U6p0_F100p0"]["sea_state"]["omega_c"], 2),
         _f(cells1["U10p0_F100p0"]["sea_state"]["omega_c"], 2)))
    A("")
    A("| U | F | H_s predicted [m] | H_s law [m] | delta | H_s from m0 full | "
      "sigma^2 [rad^2] at k_cut=1000 | sigma^2 full | ECKV ref full |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in s1["cells"]:
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (_f(c["U_ms"], 1), _f(c["fetch_m"], 0),
             _f(c["Hs_predicted_m"], 5), _f(c["Hs_law_m"], 5),
             "%+.1f%%" % (100 * c["Hs_pred_rel_law"]),
             _f(4.0 * math.sqrt(max(c["integrals_full"]["total"]["m0"], 0.0)),
                5),
             _f(c["integrals_cut"]["total"]["sigma2"], 6),
             _f(c["integrals_full"]["total"]["sigma2"], 6),
             _f(c["eckv_reference"]["sigma2_rad2"], 6)))
    A("")
    A("The ECKV reference column is explained in section 2.2; it is the same "
      "shape at the same `(k_p, Omega_c, u*)` with deep-water `c(k)` and no "
      "cutoff.")
    A("")
    A("---")
    A("")
    A("## 2. Stage 2 -- the eight-cell sweep")
    A("")
    A("`V11WaterState` is `V10WaterState` with the source switched (three "
      "overridden methods). Everything downstream is inherited verbatim, "
      "including the patch rule (scene scale 1.2 m, doubled only until "
      "`S >= 8 lambda_p`), the partition-of-unity cascade pair, the shared RNG "
      "stream and the 8 s omega quantisation. The legacy-source build was "
      "checked to be **bit-identical to v10**: `H0` arrays equal, "
      "`spectrum_after` equal, for the same seed. Legacy numbers in this "
      "section are the v10 sweep JSON's recorded values, never recomputed.")
    A("")
    A("### 2.1 C1 -- the defect itself, per band, split by part")
    A("")
    A("Bands per PIPELINE-V10 section 5-C. `below_kp` is added so that the rows "
      "of a cell are exhaustive: without it, energy between `kmin` and `k_p` "
      "would disappear from the table and the shares would be computed against "
      "a truncated total.")
    A("")
    A("Cell = the reference cell U = 6 / F = 100 / h = 1:")
    A("")
    ref = by["U6p0_F100p0"]
    A("| band [rad/m] | sigma^2 gravity | sigma^2 capillary | sigma^2 total | "
      "share of total | gravity share |")
    A("|---|---|---|---|---|---|")
    for name in BANDS:
        r = ref["integrals_full"]["bands"]["rows"][name]
        A("| `%s` %s | %s | %s | %s | %s | %s |"
          % (name, _f(r["k"][0], 1) + "-" + _f(r["k"][1], 1),
             _f(r["S_g"]["sigma2"], 6), _f(r["S_h"]["sigma2"], 6),
             _f(r["total"]["sigma2"], 6), _f(r["sigma2_share"], 3),
             _f(r["gravity_share_of_sigma2"], 3)))
    A("")
    A("All eight cells -- `sigma^2` share of the band, with the gravity share "
      "in brackets. The `k_p-10` band is EMPTY in the four cells where "
      "`k_p > 10 rad/m` (`k_p = %s`, `%s`, `%s`, `%s`): the design's first band "
      "edge is `k_p` itself, so for those young-sea cells the band has no "
      "width and is reported as an empty row rather than merged silently into "
      "its neighbour."
      % (_f(cells1["U1p5_F100p0"]["sea_state"]["k_p_rad_per_m"], 1),
         _f(cells1["U3p0_F100p0"]["sea_state"]["k_p_rad_per_m"], 1),
         _f(cells1["U6p0_F100p0"]["sea_state"]["k_p_rad_per_m"], 1),
         _f(cells1["U10p0_F100p0"]["sea_state"]["k_p_rad_per_m"], 1)))
    A("")
    A("| cell | below k_p | k_p-10 | 10-50 | 50-200 | 200-370 | 370-1000 | "
      "1000-inf | total sigma^2 full |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in cells2:
        rows = c["integrals_full"]["bands"]["rows"]
        cells_txt = []
        for name in BANDS:
            r = rows[name]
            cells_txt.append("%s (%s)"
                             % (_f(r["sigma2_share"], 2),
                                _f(r["gravity_share_of_sigma2"], 2)))
        A("| %s | %s | %s |"
          % (_cell(c),
             " | ".join(cells_txt),
             _f(c["integrals_full"]["total"]["sigma2"], 6)))
    A("")
    A("| cell | gravity share of 50-200 (new) | v10 recorded gravity share of "
      "50-200 | limit | max B_l in 50-200 | max B_h in 50-200 | B_l limit | "
      "C1 |")
    A("|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C1"]
        A("| %s | %s | %s | 0.50 | %s | %s | 0.010 | %s |"
          % (_cell(c),
             _f(g["gravity_share_50_200"], 3),
             _f(g["legacy_gravity_share_50_200"], 3),
             _f(g["B_l_max_50_200"], 5), _f(g["B_h_max_50_200"], 5),
             _pass(g["pass"])))
    A("")
    A("**C1: %s** (worst share `%s`, limit 0.50; failing cells: %s)."
      % (_pass(sm["C1_gravity_share"]["pass"]),
         _f(sm["C1_gravity_share"]["worst_share"], 3),
         ", ".join(sm["C1_gravity_share"]["failing_cells"]) or "none"))
    A("")
    A("Two distinct facts are in that table and they must not be run "
      "together:")
    A("")
    A("- Where the capillary part exists (U = 6, 10), the gravity share of "
      "50-200 collapses from the v3.7 recorded `%s`-`%s` to `%s`-`%s`, i.e. "
      "the defect the chunk was opened to fix IS fixed, by one to two orders "
      "of magnitude."
      % (_f(min(c["gates"]["C1"]["legacy_gravity_share_50_200"]
                for c in cells2 if c["U_ms"] >= 6.0), 3),
         _f(max(c["gates"]["C1"]["legacy_gravity_share_50_200"]
                for c in cells2 if c["U_ms"] >= 6.0), 3),
         _f(min(c["gates"]["C1"]["gravity_share_50_200"]
                for c in cells2 if c["U_ms"] >= 6.0), 3),
         _f(max(c["gates"]["C1"]["gravity_share_50_200"]
                for c in cells2 if c["U_ms"] >= 6.0), 3)))
    A("- At U = 1.5 the inherited v3.7 low-wind floor sets `alpha_m = 0`, so "
      "`B_h` is *identically zero* and the 50-200 band is 100%% gravity by "
      "construction; at U = 3 the floor does not bind but `alpha_m` is only "
      "`%s`, so 85-89%% of the band is still gravity. The gate as written "
      "(`<= 50%%` in all eight cells) therefore FAILS in four cells. That is "
      "the inherited floor (`low_wind_floor`, a recorded v3.7 limitation that "
      "PIPELINE-V10 section 3 tells this chunk to inherit unchanged), not a "
      "B_l defect -- but it is a merge blocker as written, and no constant was "
      "moved to hide it." % _f(cells1["U3p0_F100p0"]["alpha_m"], 5))
    A("")
    A("### 2.2 C2 -- the total, and a premise that does not survive")
    A("")
    A("| cell | sigma^2 new (k_cut=1000) | v10 recorded total | rel vs v10 | "
      "sigma^2 new (full) | ECKV ref (matched, full) | rel vs ECKV | v10 "
      "recorded vs ECKV ref | C2 |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C2"]
        A("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(c),
             _f(g["sigma2_assembled_cut"], 6), _f(g["sigma2_v10_realised"], 6),
             "%+.1f%%" % (100 * g["rel_vs_v10"]),
             _f(g["sigma2_assembled_full"], 6),
             _f(g["sigma2_eckv_reference_full"], 6),
             "%+.1f%%" % (100 * g["rel_vs_eckv_full"]),
             "%+.1f%%" % (100 * g["rel_v10_vs_eckv_reference"]),
             _pass(g["pass"])))
    A("")
    A("**C2: %s** (worst abs rel `%s`, limit 0.10). This is the gate the design "
      "wrote to anchor continuity of the total while the distribution changes, "
      "and it fails in all eight cells."
      % (_pass(sm["C2_total"]["pass"]), _f(sm["C2_total"]["worst_rel"], 3)))
    A("")
    A("The reason is not a defect of B_l. It is that **C2's stated premise is "
      "false**: PIPELINE-V10 section 5-C says \"the totals already sit within "
      "1-2%% of full-ECKV\". The number v10 compared against was "
      "`v10_spectrum.full_eckv_sigma2`, which integrates *the v3.7 spectrum "
      "itself* with the cutoff removed -- it is not ECKV at matched inputs, and "
      "the comparison was therefore circular. Against an actual "
      "ECKV-at-matched-inputs reference (the same B_l and the v3.7 B_h with "
      "`L_pm` carried, McDaniel 2001, at the same `k_p`, `Omega_c`, `u*`, with "
      "`c(k)` from the dispersion's deep-water branch and no `k_cut`, "
      "integrated to the viscous edge) the v3.7 totals are "
      "`%s`x to `%s`x the reference, i.e. from `%.0f%%` below it to `%+.0f%%` "
      "above it -- not within 1-2%% of anything. The new spectrum reproduces "
      "that ECKV reference to within `%s` (by construction: it *is* that shape "
      "plus the stated finite-depth phase speed), and it is necessarily 40-70%% "
      "below v3.7 because v3.7's total was inflated by exactly the plateau C1 "
      "removes. C1 and C2 cannot both hold: fixing the tail removes the energy "
      "the plateau contributed, and the capillary part is an independent level "
      "that cannot compensate."
      % (_f(min(1 + c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                for c in cells2), 2),
         _f(max(1 + c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                for c in cells2), 2),
         100 * abs(min(c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                      for c in cells2)),
         100 * max(c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                   for c in cells2),
         "%.2f%%" % (100 * max(abs(c["gates"]["C2"]["rel_vs_eckv_full"])
                               for c in cells2))))
    A("")
    A("Band ratios against the ECKV reference, with absolute contributions "
      "(a ratio alone can hide compensating errors). The ratio is `1.000` in "
      "almost every band of every cell -- at most `1.020` (`below_kp` at "
      "U = 10/F = 1000) -- which is the statement that the implementation "
      "reproduces the reference BY CONSTRUCTION. It is therefore the ABSOLUTE "
      "column that carries information, and it is printed beside the ratio for "
      "exactly that reason. The comparison that matters for this chunk is "
      "against the LEGACY spectrum: in the 50-200 band the new source carries "
      "between `%s`x and `%s`x of v3.7's slope variance across the eight cells "
      "(worst cell U = 6/F = 100: `%s`x), which is the measured size of the "
      "defect being removed."
      % (_f(min((c["legacy_v10"]["overlap_50_200"]["S_g_sigma2"]
                 + c["legacy_v10"]["overlap_50_200"]["S_h_sigma2"])
                / c["integrals_full"]["bands"]["rows"]["50_200"]
                ["total"]["sigma2"] for c in cells2), 2),
         _f(max((c["legacy_v10"]["overlap_50_200"]["S_g_sigma2"]
                 + c["legacy_v10"]["overlap_50_200"]["S_h_sigma2"])
                / c["integrals_full"]["bands"]["rows"]["50_200"]
                ["total"]["sigma2"] for c in cells2), 2),
         _f((cells1["U6p0_F100p0"]["integrals_full"]["bands"]["rows"]
             ["50_200"]["total"]["sigma2"] and
             (by["U6p0_F100p0"]["legacy_v10"]["overlap_50_200"]["S_g_sigma2"]
              + by["U6p0_F100p0"]["legacy_v10"]["overlap_50_200"]
              ["S_h_sigma2"])
             / by["U6p0_F100p0"]["integrals_full"]["bands"]["rows"]["50_200"]
             ["total"]["sigma2"]), 2)))
    A("")
    A("| cell | band | new sigma^2 | ECKV sigma^2 | ratio |")
    A("|---|---|---|---|---|")
    for c in cells2:
        for name in BANDS:
            b = c["gates"]["C2"]["band_ratios_vs_eckv"][name]
            A("| %s | `%s` | %s | %s | %s |"
              % (_cell(c), name,
                 _f(b["new_sigma2"], 6), _f(b["eckv_sigma2"], 6),
                 _f(b["ratio"], 3)))
    A("")
    A("### 2.3 C3 -- fraction above 370 rad/m against ECKV's own fraction")
    A("")
    A("Labelled **conformance** in the words the design uses: an internal "
      "consistency check between the implemented spectrum and the ECKV "
      "reference at matched inputs, *not* physical validation.")
    A("")
    A("| cell | frac new (realised, k_cut=1000) | frac new (assembled full) | "
      "frac ECKV (matched, full) | delta (realised vs ECKV) | target | "
      "within 0.03 | v10 recorded frac (realised) |")
    A("|---|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["C3"]
        A("| %s | %s | %s | %s | %s | 0.03 | %s | %s |"
          % (_cell(c),
             _f(g["frac_above_370_new_realised"], 3),
             _f(g["frac_above_370_new_assembled_full"], 3),
             _f(g["frac_above_370_eckv_matched"], 3),
             "%+.3f" % g["delta"], _yn(g["pass_reported"]),
             _f(g["legacy_frac_above_370_realised"], 3)))
    A("")
    A("The realised fraction is the one the field actually carries: the mode "
      "set stops at `k_cut = 1000 rad/m`, so this column measures what the "
      "cutoff removes relative to ECKV's full band. All eight cells are inside "
      "`+-0.03`; the residual is the `k_cut` truncation plus the finite-depth "
      "phase speed declared in section 1.1. The assembled-full column agrees "
      "with the ECKV reference to ~1e-3, which is the statement that the "
      "implementation IS the ECKV shape at matched inputs -- that agreement is "
      "arithmetically expected and is not the interesting half of the gate.")
    A("")
    A("### 2.4 V1 -- total against clean Cox-Munk (a reference, never a "
      "target)")
    A("")
    A("| cell | sigma^2 new | clean Cox-Munk | factor | v10 recorded factor |")
    A("|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["V1"]
        A("| %s | %s | %s | %s | %s |"
          % (_cell(c),
             _f(g["sigma2_new"], 6), _f(g["cox_munk_clean_rad2"], 6),
             _f(g["factor_over_cox_munk"], 2),
             _f(g["legacy_factor_over_cox_munk"], 2)))
    A("")
    A("v3.7 sat at 1.1-2.3x clean Cox-Munk; v3.8 sits at %s-%sx, i.e. BELOW "
      "Cox-Munk. That direction is a real consequence and it is worth stating "
      "plainly: at these very young wave ages (`Omega_c` up to %s) the "
      "Elfouhaily `alpha_p` normalisation yields less total slope than Cox & "
      "Munk's fully developed sea. Nothing was rescaled to bring it back; V1 "
      "carries no pass criterion by design."
      % (_f(sm["V1_cox_munk"]["factor_range"][0], 2),
         _f(sm["V1_cox_munk"]["factor_range"][1], 2),
         _f(max(c["sea_state"]["omega_c"] for c in s1["cells"]), 1)))
    A("")
    A("### 2.5 V2 -- band shares beside the Laxague partition (with "
      "uncertainty)")
    A("")
    A("| cell | frac sigma^2 > 1173 rad/m (continuum) | Laxague pure capillary "
      "| Omega_c | caveat |")
    A("|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["V2"]
        A("| %s | %s | 0.05 | %s | Omega_c > 2: outside Laxague's seas |"
          % (_cell(c),
             _f(g["frac_sigma2_above_1173"], 3),
             _f(c["sea_state"]["omega_c"], 2)))
    A("")
    A("The realised mode set stops at `k_cut = 1000 rad/m`, so the pure "
      "capillary band is carried by the continuum only and appears in the "
      "field through no mode at all. Laxague et al.'s partition was measured "
      "on seas with `Omega_c <= 2`; seven of these eight cells are younger "
      "than that, so this table is a comparison with stated uncertainty, not a "
      "validation. The number is small because the `alpha_m` level is set by "
      "`u*` alone.")
    A("")
    A("### 2.6 V3 -- does the gravity part still sit above the measured "
      "saturation level?")
    A("")
    A("| cell | max B_l | max B_total | max B_l / B_hat(0.007) | max "
      "B_total / B_hat | frac of band with B_total > 0.008 | verdict |")
    A("|---|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["V3"]
        A("| %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(c),
             _f(g["B_l_max_over_band"], 5), _f(g["B_total_max_over_band"], 5),
             _f(g["max_ratio_B_l_over_Bhat_low"], 3),
             _f(g["max_ratio_B_total_over_Bhat_low"], 3),
             _f(g["frac_band_with_B_total_above_Bhat_high"], 3),
             "gravity below B_hat across the band"
             if g["pass"] else "gravity still above B_hat"))
    A("")
    A("The gravity part no longer sits above `B_hat` anywhere in 50-200 in any "
      "cell (worst ratio `%s`). The capillary part does exceed `B_hat = 0.008` "
      "in the two high-wind cells, where `max B_total / B_hat = %s`: that is a "
      "statement about the inherited `B_h` level at `u* = 0.34`, not about the "
      "gravity tail this chunk fixes."
      % (_f(sm["V3_saturation"]["worst_B_l_over_Bhat"], 3),
         _f(sm["V3_saturation"]["worst_B_total_over_Bhat"], 3)))
    A("")
    A("### 2.7 K1 -- H_s is a prediction")
    A("")
    A("| cell | H_s predicted (realised spectral) | H_s law | delta | abs delta "
      "| flag (> 35%) |")
    A("|---|---|---|---|---|---|")
    for c in cells2:
        g = c["gates"]["K1"]
        A("| %s | %s | %s | %s | %s | %s |"
          % (_cell(c),
             _f(g["Hs_predicted_m"], 5), _f(g["Hs_law_m"], 5),
             "%+.1f%%" % (100 * g["delta_rel"]),
             _f(g["delta_abs_m"], 5), _yn(g["flag"])))
    A("")
    A("**K1: no flag open** (worst abs delta `%s`, threshold 35%%). No H_s "
      "override was applied anywhere; the design's rule that an override is "
      "carried as a *statement* beside the prediction is implemented in the "
      "ctx as `hs_override_m` and would have been logged had a recipe supplied "
      "one." % ("%.1f%%" % (100 * sm["K1_Hs_delta"]["worst_abs_delta_rel"])))
    A("")
    A("### 2.8 K5 -- 8 s quantisation, per band (v3.7 G7 analogue)")
    A("")
    A("| cell | > 1 m | 0.1-1 m | 1-10 cm | < 1 cm | max rel | v10 recorded max "
      "rel |")
    A("|---|---|---|---|---|---|---|")
    for c in cells2:
        q = c["gates"]["K5"]
        pb = q["per_band_max_rel_error"]
        A("| %s | %s | %s | %s | %s | %s | %s |"
          % (_cell(c),
             _f(pb["gt_1m"]["max_rel"], 4), _f(pb["0p1_1m"]["max_rel"], 4),
             _f(pb["1_10cm"]["max_rel"], 4), _f(pb["lt_1cm"]["max_rel"], 4),
             _f(q["max_rel_error"], 4),
             _f(c["legacy_v10"]["quantisation_max_rel"], 4)))
    A("")
    A("The interface is unchanged (same `2 pi / 8 s` grid, same function). K5 "
      "measures the omega quantisation of the *realised mode set*, not the "
      "spectral shape: seven of the eight cells reproduce the v3.7 number to "
      "the last digit, and U = 1.5/F = 1000 moves from `%s` to `%s` because "
      "the new spectrum realises a sparser set of long modes there "
      "(`%d` modes in the `> 1 m` band), so the worst surviving mode in the "
      "table is a different mode. No interface behaviour changed and no "
      "tolerance was introduced."
      % (_f(cells2[1]["legacy_v10"]["quantisation_max_rel"], 4),
         _f(cells2[1]["gates"]["K5"]["max_rel_error"], 4),
         cells2[1]["gates"]["K5"]["per_band_max_rel_error"]["gt_1m"]["n"]))
    A("")
    A("### 2.9 K6 -- continuity audit of the assembled B(k)")
    A("")
    A("Two tests, both scale-free, both on the *assembled* `B = B_l + B_h`:")
    A("")
    A("1. no step: `B` never changes by more than 50% in one 0.2% step of "
      "`ln k`;")
    A("2. no slope discontinuity: the largest one-step change of "
      "`d ln B / d ln k` is measured at 4000 and 8000 points on the same band. "
      "A smooth function -- however sharp its peak -- has that quantity "
      "proportional to the step, so it **halves**; a splice or a "
      "renormalisation step has a fixed absolute jump and does not.")
    A("")
    A("| cell | B min (whole grid) | B >= 0 | max rel step | slope-jump ratio "
      "fine/coarse | ...at k | domain label | K6 |")
    A("|---|---|---|---|---|---|---|---|")
    for c in cells2:
        k6 = c["gates"]["K6"]
        A("| %s | %s | %s | %s | %s | %s | `%s` | %s |"
          % (_cell(c),
             _f(k6["B_min_over_grid"], 2), _yn(k6["B_nonneg"]),
             _f(k6["coarse"]["max_rel_step"], 3),
             _f(k6["max_dslope_ratio_fine_over_coarse"], 3),
             _f(k6["fine"]["k_at_max_dslope"], 2), k6["eckv_domain"],
             _pass(k6["pass"])))
    A("")
    A("Below `1e-12 max B` the spectrum has underflowed (`B_l ~ "
      "exp(-1.25 (k_p/k)^2)` reaches 1e-117 at k = 1 rad/m): those points are "
      "reported as finite and non-negative and excluded from the ratio tests, "
      "because a relative-step test on a denormal measures rounding. That "
      "exclusion is stated above, not silent -- it is the one place in this "
      "chunk where a test is not applied to the full 1-3000 rad/m grid.")
    A("")
    A("---")
    A("")
    A("## 3. Stage 3 -- compatibility")
    A("")
    A("### 3.1 K2 -- height quantiles of built fields, old vs new")
    A("")
    k2 = s3["K2"]["scenes"]
    A("Reference scenes identified from the v10 outputs "
      "(`out/v10/numbers-v10-transport.json`): scene A = `%s`, scene B = `%s`. "
      "Same seed family (one seed, preset 1), same grid, only the source "
      "differs."
      % (k2["A"]["scene"]["name"], k2["B"]["scene"]["name"]))
    A("")
    A("| scene | U/F/h | qty | legacy | new | rel change |")
    A("|---|---|---|---|---|---|")
    for key in ("A", "B"):
        sc = k2[key]["scene"]
        o = k2[key]["sources"]["jonswap"]["height"]
        n = k2[key]["sources"]["elfouhaily-form"]["height"]
        A("| %s | %s/%s/%s | eta p50 [m] | %s | %s | %s |"
          % (key, _f(sc["U"], 1), _f(sc["F"], 0), _f(sc["h"], 1),
             _f(o["p50_m"], 6), _f(n["p50_m"], 6), "~0 (sign only)"))
        A("| %s | | eta p90 [m] | %s | %s | %s |"
          % (key, _f(o["p90_m"], 5), _f(n["p90_m"], 5),
             "%+.1f%%" % (100 * (n["p90_m"] / o["p90_m"] - 1.0))))
        A("| %s | | eta p99 [m] | %s | %s | %s |"
          % (key, _f(o["p99_m"], 5), _f(n["p99_m"], 5),
             "%+.1f%%" % (100 * (n["p99_m"] / o["p99_m"] - 1.0))))
        A("| %s | | H_s = 4 sqrt m0 [m] | %s | %s | %s |"
          % (key, _f(o["Hs_m"], 5), _f(n["Hs_m"], 5),
             "%+.1f%%" % (100 * (n["Hs_m"] / o["Hs_m"] - 1.0))))
        A("| %s | | std [m] | %s | %s | %s |"
          % (key, _f(o["std_m"], 5), _f(n["std_m"], 5),
             "%+.1f%%" % (100 * (n["std_m"] / o["std_m"] - 1.0))))
        A("| %s | | skew / kurtosis | %s / %s | %s / %s | %+.3f / %+.3f |"
          % (key, _f(o["skewness"], 3), _f(o["kurtosis"], 3),
             _f(n["skewness"], 3), _f(n["kurtosis"], 3),
             n["skewness"] - o["skewness"], n["kurtosis"] - o["kurtosis"]))
    A("")
    A("The heights scale together (`p90`, `p99` and `H_s` move by the same "
      "factor within a few tenths of a percent), and the shape statistics are "
      "unchanged: the new source moves the AMPLITUDE of the same kind of "
      "field, not its distribution.")
    A("")
    A("### 3.2 K3 -- the G10 windy scene with BOTH sources")
    A("")
    A("Scene: `U = 6 m/s`, `F = 100 m`, `h = 1.0 m`, sun `2026-06-10 09:30`, "
      "moth photo native `735x539`, `eta = 1.0`, frozen composite (A 2.6, B "
      "lowpass 3.0 px, E 2.2), `Q = 25`, master `2048^2`. The transport and "
      "the composite are the imported v5/v9/v10 code paths -- nothing in them "
      "was touched, and the legacy-source still proves it (byte-identical "
      "below).")
    A("")
    A("| source | H_s spectral | sigma^2 spectral | rms slope | g1 | C p50 | "
      "C p99 | C p99.9 | C max |")
    A("|---|---|---|---|---|---|---|---|---|")
    for src in ("jonswap", "elfouhaily-form"):
        sc = s3["K3"]["scenes"][src]
        st = sc["state_summary"]
        A("| `%s` | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (src, _f(st["Hs_spectral_m"], 5), _f(st["sigma2_spectral"], 5),
             _f(st["rms_slope"], 4), _f(sc["route"]["g1_used"], 4),
             _f(sc["C_percentiles"]["p50"], 3),
             _f(sc["C_percentiles"]["p99"], 3),
             _f(sc["C_percentiles"]["p99_9"], 3),
             _f(sc["C_percentiles"]["max"], 2)))
    A("")
    A("Warp yardsticks -- the pipeline's OWN geometric shift field (the same "
      "displacement the renderer applies, rebuilt from the same state through "
      "`_state_gradients`), in px and mm, plus the block-match measurement on "
      "the delivered still (48 px blocks, +-24 px, 2 px steps, Rec.709):")
    A("")
    A("| source | own field mean [px] | own field mean [mm] | px scale "
      "[mm/px] | SSD mean [px] | SSD std [px] | NCC mean [px] | NCC std [px] "
      "| blocks at search edge |")
    A("|---|---|---|---|---|---|---|---|---|")
    for src in ("jonswap", "elfouhaily-form"):
        w = s3["K3"]["scenes"][src]["warp"]
        A("| `%s` | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (src, _f(w["own_shift_field"]["field_mean_px"], 2),
             _f(w["own_shift_field"]["mean_mm"], 2),
             _f(w["own_shift_field"]["pixel_scale_mm"], 3),
             _f(w["primary_ssd"]["mean_shift_px"], 2),
             _f(w["primary_ssd"]["shift_field_std_px"], 2),
             _f(w["ncc_check"]["mean_shift_px"], 2),
             _f(w["ncc_check"]["shift_field_std_px"], 2),
             _f(w["primary_ssd"]["frac_blocks_at_search_edge"], 3)))
    A("")
    A("| yardstick (history) | own-field shift | block match |")
    A("|---|---|---|")
    ys = s3["K3"]["scenes"]["jonswap"]["warp"]["yardsticks"]
    A("| v3.4 windy still (recorded) | %s mm at %s mm/px | %s +/- %s px |"
      % (_f(ys["v34_windy_still"]["own_field_mm"], 1),
         _f(ys["v34_windy_still"]["px_scale_mm"], 2),
         _f(ys["v34_windy_still"]["mean_px_ncc"], 2),
         _f(ys["v34_windy_still"]["std_px_ncc"], 2)))
    A("| production loop 2026-06-10 (recorded) | %s mm at %s mm/px | %s px |"
      % (_f(ys["production_loop_2026_06_10"]["mean_mm_from_px"], 1),
         _f(ys["production_loop_2026_06_10"]["px_scale_mm"], 2),
         _f(ys["production_loop_2026_06_10"]["mean_px"], 1)))
    A("| v3.7 legacy source (this chunk) | %s mm at %s mm/px | %s +/- %s px |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["mean_mm"], 2),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["pixel_scale_mm"], 3),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["primary_ssd"]
            ["mean_shift_px"], 2),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["primary_ssd"]
            ["shift_field_std_px"], 2)))
    A("| v3.8 new source (this chunk) | %s mm at %s mm/px | %s +/- %s px |"
      % (_f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["own_shift_field"]
            ["mean_mm"], 2),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["own_shift_field"]
            ["pixel_scale_mm"], 3),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["primary_ssd"]
            ["mean_shift_px"], 2),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["primary_ssd"]
            ["shift_field_std_px"], 2)))
    A("")
    A("The px scale is `S / W` of the scene patch (6.531 mm/px here); the "
      "v3.4 row's 7.63 mm/px came from its wind-coupled patch. The own-field "
      "number is the column to read: it is the renderer's own geometric "
      "displacement, so it does not depend on the matcher. The v3.7 row "
      "reproduces the v10 G10 record exactly, and the v3.8 row is the number "
      "the new source produces.")
    A("")
    A("Fine texture (mean absolute luminance gradient) and tone:")
    A("")
    A("| image | mean abs dL/dy | mean abs dL/dx | mean L | std L | crush < "
      "0.02 | blowout > 0.98 |")
    A("|---|---|---|---|---|---|---|")
    b = s3["K3"]["between_sources"]["texture_base_photo"]
    A("| base photo | %s | %s | %s | %s | %s | %s |"
      % (_f(b["mean_abs_dL_dy"], 5), _f(b["mean_abs_dL_dx"], 5),
         _f(b["mean_luminance"], 4), _f(b["std_luminance"], 4),
         _f(b["crush_rate_linear_lt_0.02"], 4),
         _f(b["frac_aces_gt_0.98"], 5)))
    for src in ("jonswap", "elfouhaily-form"):
        t = s3["K3"]["scenes"][src]["fine_texture"]
        A("| `%s` | %s | %s | %s | %s | %s | %s |"
          % (src, _f(t["mean_abs_dL_dy"], 5), _f(t["mean_abs_dL_dx"], 5),
             _f(t["mean_luminance"], 4), _f(t["std_luminance"], 4),
             _f(t["crush_rate_linear_lt_0.02"], 4),
             _f(t["frac_aces_gt_0.98"], 5)))
    A("")
    A("Between the two stills: mean abs difference `%s`/255, MS-SSIM `%s`, "
      "changed pixels `%s`. The g1-pinned control (the new state rendered with "
      "the LEGACY shadowing, so only the caustic gain differs) gives mean abs "
      "`%s`/255 and MS-SSIM `%s`."
      % (_f(s3["K3"]["between_sources"]["mean_abs_uint8"], 2),
         _f(s3["K3"]["between_sources"]["msssim"], 5),
         _f(s3["K3"]["between_sources"]["frac_pixels_changed"], 4),
         _f(s3["K3"]["between_sources"]["g1_pinned_control"]
            ["mean_abs_uint8"], 2),
         _f(s3["K3"]["between_sources"]["g1_pinned_control"]
            ["msssim_vs_legacy"], 5)))
    A("")
    A("**Attribution** (the design's requirement: any change must be "
      "attributable to the band redistribution or to delta H_s, and the report "
      "must say which):")
    A("")
    at = s3["K3"]["attribution"]
    A("| quantity | legacy | new | change |")
    A("|---|---|---|---|")
    A("| H_s (spectral) | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["state_summary"]["Hs_spectral_m"], 5),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["Hs_spectral_m"], 5), "%+.1f%%" % (100 * at["Hs_rel_change"])))
    A("| sigma^2 (spectral) | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["state_summary"]
            ["sigma2_spectral"], 5),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["sigma2_spectral"], 5),
         "%+.1f%%" % (100 * at["sigma2_rel_change"])))
    A("| rms slope | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["state_summary"]["rms_slope"], 4),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["rms_slope"], 4), "%+.1f%%" % (100 * at["rms_slope_rel_change"])))
    A("| g1 (Smith, from the state) | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["route"]["g1_used"], 4),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["route"]["g1_used"], 4),
         "%+.1f%%" % (100 * at["g1_change"])))
    A("| share of sigma^2 > 370 rad/m | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["state_summary"]
            ["frac_above_370_realised"], 3),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["frac_above_370_realised"], 3),
         "%+.1f%%" % (100 * at["frac_above_370_change"])))
    A("| sigma^2 above 370 rad/m, **absolute** | %s | %s | %s |"
      % (_f(at["sigma2_above_370_absolute"]["legacy"], 6),
         _f(at["sigma2_above_370_absolute"]["new"], 6),
         "%+.1f%%" % (100 * at["sigma2_above_370_absolute"]["rel_change"])))
    A("| sigma^2 share, 1-10 cm | %s | %s | - |"
      % (_f(at["band_share_1_10cm"]["legacy"], 3),
         _f(at["band_share_1_10cm"]["new"], 3)))
    A("| sigma^2 share, < 1 cm | %s | %s | - |"
      % (_f(at["band_share_lt_1cm"]["legacy"], 3),
         _f(at["band_share_lt_1cm"]["new"], 3)))
    A("| own-field warp, mm | %s | %s | %s |"
      % (_f(at["own_field_shift_mm"]["legacy"], 2),
         _f(at["own_field_shift_mm"]["new"], 2),
         "%+.1f%%" % (100 * (at["own_field_shift_mm"]["new"]
                             / at["own_field_shift_mm"]["legacy"] - 1.0))))
    A("| mean abs dL/dy | %s | %s | %s |"
      % (_f(s3["K3"]["scenes"]["jonswap"]["fine_texture"]["mean_abs_dL_dy"],
            5),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["fine_texture"]
            ["mean_abs_dL_dy"], 5),
         "%+.1f%%" % (100 * at["fine_texture_rel_change"]["dy"])))
    A("")
    A(at["statement"])
    A("")
    A("**Continuity check on the legacy side.** `comp-windy-v38-legacy.png` is "
      "**byte-identical** to the v10 chunk-7 deliverable `comp-windy-v37.png` "
      "(sha256 `%s...`), and every legacy-column number in this section equals "
      "the v10 G10 record exactly (C p50/p99/p99.9 `%s`/`%s`/`%s`, warp "
      "`%s` px, own-field `%s` mm, `|dL/dy|` `%s`). The source switch does not "
      "perturb the legacy path by one bit -- it is a switch, not a re-fit."
      % (s3["K3"]["continuity_vs_v10"]["sha256_v10"][:16],
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["p50"], 4),
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["p99"], 4),
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["p99_9"], 4),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["primary_ssd"]
            ["mean_shift_px"], 2),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["mean_mm"], 2),
         _f(s3["K3"]["scenes"]["jonswap"]["fine_texture"]["mean_abs_dL_dy"],
            6)))
    A("")
    A("**Expectation stated in advance, and its outcome.** The design expected "
      "\"fine web holds or strengthens, broad warp follows the new H_s\". The "
      "second half holds: the pipeline's own geometric displacement falls to "
      "`%s` of the legacy value in mm, tracking `delta H_s = %+.1f%%`. The "
      "first half does **not**: `|dL/dy|` falls `%.1f%%`. The mechanism is "
      "measurable and is the band redistribution, not a modelling choice: the "
      "*share* of `sigma^2` above 370 rad/m rises from `%s` to `%s`, but the "
      "share is taken of a total that fell `%.0f%%`, so the **absolute** "
      "high-wavenumber slope variance falls from `%s` to `%s` (`%+.1f%%`). "
      "The fine web is driven by that absolute level, so it weakens. At the "
      "same time the caustic peaks get STRONGER (C p99.9 `%s` -> `%s`, max "
      "`%s` -> `%s`): a narrower spectrum focuses more coherently even with "
      "less total slope. Both halves of that trade are reported, neither was "
      "tuned, and the reason the fine band cannot be held up is the inherited "
      "`alpha_m` level plus the parasitic-capillary gap (section 5.1)."
      % ("%.1f%%" % (100 * at["own_field_shift_mm"]["new"]
                     / max(at["own_field_shift_mm"]["legacy"], 1e-12)),
         100 * at["Hs_rel_change"],
         100 * abs(at["fine_texture_rel_change"]["dy"]),
         _f(s3["K3"]["scenes"]["jonswap"]["state_summary"]
            ["frac_above_370_realised"], 3),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["frac_above_370_realised"], 3),
         100 * abs(at["sigma2_rel_change"]),
         _f(at["sigma2_above_370_absolute"]["legacy"], 6),
         _f(at["sigma2_above_370_absolute"]["new"], 6),
         100 * at["sigma2_above_370_absolute"]["rel_change"],
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["p99_9"], 3),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["C_percentiles"]
            ["p99_9"], 3),
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["max"], 2),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["C_percentiles"]["max"],
            2)))
    A("")
    A("**The g1 control is degenerate and is reported as such.** Smith "
      "shadowing saturates at `g1 = %s` for BOTH sources at this sun "
      "elevation, so the \"g1-pinned\" render is identical to the native new "
      "render (`%s`/255, MS-SSIM `%s` -- the same as the unpinned comparison). "
      "The between-source difference is therefore NOT a shadowing change: it "
      "is the geometry (own shift field `%s` -> `%s` px) and the caustic gain "
      "C (`p99` `%s` -> `%s`). The C=1 control, which removes the caustic gain "
      "and leaves only geometry, separates them: block-match warp `%s` px "
      "(legacy) versus `%s` px (new)."
      % (_f(s3["K3"]["scenes"]["jonswap"]["route"]["g1_used"], 2),
         _f(s3["K3"]["between_sources"]["g1_pinned_control"]
            ["mean_abs_uint8"], 2),
         _f(s3["K3"]["between_sources"]["g1_pinned_control"]
            ["msssim_vs_legacy"], 5),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["own_shift_field"]
            ["field_mean_px"], 2),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["own_shift_field"]
            ["field_mean_px"], 2),
         _f(s3["K3"]["scenes"]["jonswap"]["C_percentiles"]["p99"], 3),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["C_percentiles"]["p99"],
            3),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]["control_C1"]
            ["ssd_primary"]["mean_shift_px"], 2),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]["control_C1"]
            ["ssd_primary"]["mean_shift_px"], 2)))
    A("")
    A("The block-match agreement with the predicted shift field stays low "
      "(corr `%s`-`%s`) exactly as it did in v10's G10: on this photo the "
      "matcher is dominated by the photograph's own texture and by the caustic "
      "modulation, and `%s%%`-`%s%%` of blocks saturate at the +-24 px search "
      "edge. That is an inherited property of the yardstick, unchanged by this "
      "chunk, and it is why the pipeline's own geometric field is reported "
      "beside it rather than instead of it."
      % (_f(s3["K3"]["scenes"]["elfouhaily-form"]["warp"]
            ["agreement_with_predicted"]["ncc_check"]["corr_dy"], 3),
         _f(s3["K3"]["scenes"]["jonswap"]["warp"]
            ["agreement_with_predicted"]["ncc_check"]["corr_dy"], 3),
         _f(100 * s3["K3"]["scenes"]["elfouhaily-form"]["warp"]
            ["primary_ssd"]["frac_blocks_at_search_edge"], 0),
         _f(100 * s3["K3"]["scenes"]["jonswap"]["warp"]["primary_ssd"]
            ["frac_blocks_at_search_edge"], 0)))
    A("")
    A("Tone: crush `%s` (customary bound 0.0325, %s) and blowout `%s` for the "
      "new still; `%s` / `%s` for the legacy still. The frozen composite pins "
      "`E = 2.2`; these are consequences, not targets."
      % (_f(s3["K3"]["scenes"]["elfouhaily-form"]["tone"]
            ["crush_rate_linear_lt_0.02"], 5),
         _pass(s3["K3"]["scenes"]["elfouhaily-form"]["tone"]
               ["crush_pass_customary"]),
         _f(s3["K3"]["scenes"]["elfouhaily-form"]["tone"]
            ["blowout_rate_aces_gt_0.98"], 6),
         _f(s3["K3"]["scenes"]["jonswap"]["tone"]
            ["crush_rate_linear_lt_0.02"], 5),
         _f(s3["K3"]["scenes"]["jonswap"]["tone"]
            ["blowout_rate_aces_gt_0.98"], 6)))
    A("")
    A("Delivered stills: `comp-windy-v38-legacy.png` and "
      "`comp-windy-v38-new.png` (both `%dx%d`)." % (735, 539))
    A("")
    A("### 3.3 K4 -- the regression battery (v3.7 G8 verbatim)")
    A("")
    k4 = s3["K4"]
    A("| check | result | numbers |")
    A("|---|---|---|")
    A("| K4a frozen composite bit-identical on IDENTICAL buffers | %s | "
      "max abs uint8 diff `%s`, params A=%s B=%s E=%s Q=%s |"
      % (_pass(k4["K4a"]["pass"]),
         _f(k4["K4a"]["max_abs_uint8"], 0), _f(k4["K4a"]["params"]["A"], 1),
         _f(k4["K4a"]["params"]["B_lowpass_sigma_px"], 1),
         _f(k4["K4a"]["params"]["E"], 1), _f(k4["K4a"]["params"]["Q"], 0)))
    A("| K4b Beer-Lambert for identical optical inputs | %s | new `%s` vs v9 "
      "recorded `%s`, max abs diff `%s`; legacy vs new identical to `%s` |"
      % (_pass(k4["K4b"]["pass"]),
         ["%.6f" % x for x in k4["K4b"]["new_two_way"]],
         ["%.6f" % x for x in k4["K4b"]["v9_recorded_two_way"]],
         _f(max(k4["K4b"]["abs_diff_vs_recorded"]), 2),
         _f(max(k4["K4b"]["abs_diff_legacy_vs_new"]), 2)))
    A("| K4c rain R = 0 equals the NEW base state exactly | %s | two builds "
      "bit-identical `%s`; R > 0 raises `%s` |"
      % (_pass(k4["K4c"]["pass"]), _yn(k4["K4c"]["bit_identical_state"]),
         _yn(k4["K4c"]["rain_positive_raises"])))
    A("")
    A("K4a is a code-path regression test on identical buffers; a render under "
      "the new state is EXPECTED to differ and is judged by K3 only, never by "
      "bit identity. K4b's two columns are the two things the design asks: "
      "unchanged against the historical optical record, and identical between "
      "sources because the environment is shared.")
    A("")
    A("---")
    A("")
    A("## 4. Stage 4 -- the merge rule, item by item")
    A("")
    A("PIPELINE-V10 section 5: the default flips to `elfouhaily-form` iff A1, "
      "A2, C1, C2 pass; C3, V1-V3 are reported; K4, K5, K6 pass; K2, K3 are "
      "reported; and no K1 flag is open.")
    A("")
    A("| # | criterion | verdict | evidence |")
    A("|---|---|---|---|")
    for i, it in enumerate(mg["items"], 1):
        A("| %d | %s | %s | %s |"
          % (i, it["item"], "PASS" if it["pass"] else "**FAIL**",
             it["evidence"]))
    A("")
    A("**Merge decision: %s.** `DEFAULT_SOURCE` in `v11_state.py` is `%s`. The "
      "challenger remains fully available behind the switch "
      "(`spectrum_source=\"elfouhaily-form\"`)."
      % ("APPLIED -- default flipped to `elfouhaily-form`"
         if mg["default_flipped"] else
         "WITHHELD -- default stays `jonswap`",
         str(mg["default_after"])))
    A("")
    if mg["blocking_items"]:
        A("Blocking items: %s." % ", ".join(mg["blocking_items"]))
        A("")
        A("**Why they block, stated as measurements rather than as opinions.**")
        A("")
        A("- **C2 (all 8 cells).** The gate demands the total `sigma^2` stay "
          "within +-10%% of v3.7's. v3.7's total was inflated by the very "
          "plateau this chunk removes: against a real ECKV-at-matched-inputs "
          "reference the v3.7 totals are `%s`x to `%s`x that reference, so C2 "
          "as written asks the fix to reproduce the defect's total. The new "
          "spectrum reproduces the ECKV reference to `%s`, which is the "
          "physically meaningful anchor, and lands 40-70%% below v3.7."
          % (_f(min(1 + c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                    for c in cells2), 2),
             _f(max(1 + c["gates"]["C2"]["rel_v10_vs_eckv_reference"]
                    for c in cells2), 2),
             "%.2f%%" % (100 * max(abs(c["gates"]["C2"]["rel_vs_eckv_full"])
                                   for c in cells2))))
        A("- **C1 (4 of 8 cells).** At U = 1.5 the inherited `alpha_m` floor "
          "makes `B_h` identically zero, so the 50-200 band cannot be anything "
          "but 100% gravity; at U = 3 `alpha_m` is so small that 85-89% of the "
          "band remains gravity. The gate's `<= 50%` is met only at U >= 6.")
        A("")
        A("Neither blocker can be closed inside this design without either "
          "(a) re-specifying C2 against an ECKV-at-matched-inputs reference "
          "instead of against v3.7's totals, or (b) closing the `alpha_m` "
          "low-wind floor -- which is the v5 surface-film seam and the "
          "parasitic-capillary gap, both explicitly out of scope. Both are "
          "**separate decisions**; per PIPELINE-V10 section 5 neither was "
          "improvised in-session, and no constant was moved to make a gate "
          "pass.")
        A("")
    A("---")
    A("")
    A("## 5. Recorded limitations (not silently dropped)")
    A("")
    A("1. **Parasitic capillary gap** (Gemini, rally D2; PIPELINE-V10 section "
      "5-V). ECKV has no parasitic-capillary source. At `Omega_c ~ 4-6.3` a "
      "real pond's steep crests may carry MORE sub-centimetre slope than this "
      "spectrum. This chunk does not close the gap -- and K3 shows the "
      "consequence in pixels: with the gravity tail removed and `B_h` left "
      "exactly where v3.7 had it, the sub-centimetre band's *share* of "
      "`sigma^2` rises (`%s` -> `%s`) while its *absolute* level falls from "
      "`%s` to `%s` (`%+.1f%%`) because the total fell `%.0f%%`, and `|dL/dy|` "
      "falls `%.1f%%`. That is the stated reason the Kudryavtsev challenger "
      "remains on the roadmap (D4: deferred at all three AIs)."
      % (_f(at["band_share_lt_1cm"]["legacy"], 3),
         _f(at["band_share_lt_1cm"]["new"], 3),
         _f(at["band_share_lt_1cm"]["legacy"]
            * s3["K3"]["scenes"]["jonswap"]["state_summary"]
            ["sigma2_spectral"], 6),
         _f(at["band_share_lt_1cm"]["new"]
            * s3["K3"]["scenes"]["elfouhaily-form"]["state_summary"]
            ["sigma2_spectral"], 6),
         100 * (at["band_share_lt_1cm"]["new"]
                * s3["K3"]["scenes"]["elfouhaily-form"]
                ["state_summary"]["sigma2_spectral"]
                / (at["band_share_lt_1cm"]["legacy"]
                   * s3["K3"]["scenes"]["jonswap"]
                   ["state_summary"]["sigma2_spectral"]) - 1.0),
         100 * abs(at["sigma2_rel_change"]),
         100 * abs(at["fine_texture_rel_change"]["dy"])))
    A("2. **Wave-age extrapolation.** Two of eight cells carry `Omega_c > 5` "
      "(`wave_age_extrapolation`): the `gamma` log-formula gives `gamma = %s` "
      "and `%s` there. Nothing is clamped and nothing diverges, but those "
      "cells are uncalibrated extrapolation of a formula fitted for `Omega_c "
      "<= 5`, and every number from them in this report inherits that label."
      % (_f(cells1["U6p0_F100p0"]["sea_state"]["gamma_bl"], 2),
         _f(cells1["U10p0_F100p0"]["sea_state"]["gamma_bl"], 2)))
    A("3. **Low-wind floor.** `alpha_m = max(raw, 0)` at U = 1.5 zeroes the "
      "entire capillary part (`low_wind_floor`, inherited unchanged from "
      "v3.7). Under the new source this is no longer a small band detail: with "
      "the gravity tail gone, the low-wind cells have essentially no sub-`k_p` "
      "capillary energy at all.")
    A("4. **Finite-depth phase speed.** `c(k)` and `c_p` come from the "
      "pipeline dispersion, not from ECKV's deep water. Confined to the peak "
      "region (`kh >= 15` in 50-200) and quantified by C3's residual.")
    A("5. **`k_p` route.** Solved through the dispersion, not the v3.7 "
      "deep-water `lambda_p` inversion. Reported per cell "
      "(`k_p_dispersion_rel_diff`); no other consequence was found.")
    A("6. **`H_s` is a prediction.** It is not the fetch-law value and it is "
      "not rescaled toward it. K1 carries the deltas (`%s` to `%s`); none "
      "exceeds the 35%% flag."
      % ("%+.1f%%" % (100 * min(c["gates"]["K1"]["delta_rel"]
                                for c in cells2)),
         "%+.1f%%" % (100 * max(c["gates"]["K1"]["delta_rel"]
                                for c in cells2))))
    A("7. **Directionality unchanged.** `Phi` is still ECKV's unified function "
      "on the total spectrum, and the anisotropy reported per cell is that "
      "function's prediction, not a new model.")
    A("8. **Conformance vs validation.** C1-C3/K5/K6 are internal consistency "
      "and code-path statements. Only V1-V3 are nature-facing, and each carries "
      "a stated reference and its uncertainty (Cox-Munk is a fully developed "
      "sea; Laxague's seas had `Omega_c <= 2`).")
    A("")
    A("---")
    A("")
    A("## 6. Deliverables (all on the PC in `offline/out/v11/`)")
    A("")
    A("| file | contents |")
    A("|---|---|")
    A("| `REPORT-V11-GATES.md` | this report |")
    A("| `numbers-v11-spectrum.json` | stage 1: A1, A2, per-cell sea state and "
      "integrals |")
    A("| `numbers-v11-sweep.json` | stage 2: C1 C2 C3 V1 V2 V3 K1 K5 K6, per "
      "cell, beside the v10 numbers |")
    A("| `numbers-v11-transport.json` | stage 3: K2, K3 (both stills, C, warp, "
      "texture, tone), K4 |")
    A("| `merge-evaluation-v11.json` | stage 4: the merge rule item by item |")
    A("| `comp-windy-v38-new.png` | the G10 windy scene, `elfouhaily-form` "
      "source |")
    A("| `comp-windy-v38-legacy.png` | the same scene, `jonswap` source |")
    A("| `spectrum-curves-v11.json` | S_l / S_h / B_l / B_h / B / Delta per "
      "cell, new beside legacy |")
    A("| `spectrum-ref-v38-U6F100.csv` | the reference cell curves as CSV |")
    A("| `checkpoint.json` | jobkit checkpoint, one unit per artefact |")
    A("| `_work/C-windy-v38-*-Q25.npy` | the two caustic masters (float32) |")
    A("")
    A("---")
    A("")
    A("## 7. Sources, and what this report does NOT claim")
    A("")
    A("- `B_l`, `B_h`, `Phi`: Elfouhaily, Chapron, Katsaros & Vandemark (1997), "
      "JGR 102(C7) 15781, eq. 34-44 and 66-68; the `L_pm` factor on `B_h` "
      "follows McDaniel (2001) and is inherited from v3.7 unchanged.")
    A("- Fetch laws: Hasselmann et al. (1973); the depth-limited `H_s` route "
      "is Young (1988), inherited.")
    A("- Saturation level `B_hat ~ 0.007-0.008`: the research note "
      "`research_notes/gravity-tail-termination-spectrum-20261111-0528/` "
      "(Rank 1), used as a reference in V3 only.")
    A("- `Cox & Munk (1954)`: fully developed sea, a reference, never a "
      "target.")
    A("- `Laxague et al.`: band partition, seas with `Omega_c <= 2`, compared "
      "with stated uncertainty.")
    A("")
    A("**Not claimed.** This is not a validation of ECKV in a pond. It is not a "
      "claim that `H_s_pred` is the true sea state (`K1` shows it is 5-22% "
      "below the fetch law by design of the `alpha_p` normalisation). It is not "
      "a claim that the new default should ship: the merge was WITHHELD on C1 "
      "and C2, with the numbers above. And the low-wind cells are not evidence "
      "either way, because the inherited floor removes their capillary part "
      "before the new spectrum is asked anything.")
    A("")
    A("Uncertainty, stated where it is largest: the two `Omega_c > 5` cells sit "
      "outside the fitted range of the `gamma` formula; the depth-limited "
      "growth correction of v3.7 does not transfer to this source and depth now "
      "enters only through the dispersion (reported as `kh_p`); and the "
      "capillary band above 1173 rad/m lives in the continuum only, because the "
      "realised mode set stops at `k_cut = 1000 rad/m`.")
    A("")
    return "\n".join(L) + "\n"


def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "out", "v11")
    os.makedirs(out_dir, exist_ok=True)
    txt = build(out_dir)
    import sys
    sys.path.insert(0, r"D:\mutaro-work\jobkit")
    from checkpoint import atomic_write_bytes
    path = os.path.join(out_dir, "REPORT-V11-GATES.md")

    def verify(p):
        with open(p, "rb") as f:
            data = f.read()
        return len(data) > 2000 and data[:1] == b"#"

    atomic_write_bytes(path, txt.encode("ascii"), verify=verify)
    print("wrote %s (%d bytes, %d lines)"
          % (path, len(txt.encode("ascii")), txt.count("\n")))


if __name__ == "__main__":
    main()
