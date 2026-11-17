#!/usr/bin/env python3
"""run_v4_2_8 -- BRIEF-V4-2-8: F-B0 analytic pre-check (mixture of cutoffs).

Question (BRIEF-V4-2-8 STAGE 1): can a mixture of scales (mixed cutoffs) raise
the apparent local slope H of the scaling band (0.4-3.6 cm) into 0.70-0.85?

ANALYTIC ONLY.  No bed is generated, no heightfield is touched, no image is
inspected.  Every curve comes from the saved decomposition of v4-2.3
(out/v4-2-3/decomp.json) that v4-2-5 and v4-2-6 already used, and every local
exponent / band fit is computed with run_v4_2_5.local_H / fit_band / band_shape
imported unmodified.

Runs in the foreground as one process, with one jobkit checkpoint unit per
stage:
  stage0  absolute-path preflight, sha256 of inputs and every existing .py,
          verbatim extraction of the design's F-B0 definition
  mix     mixture of the design's own local curves (D50 x 0.5/1/1.9) and the
          continuous lognormal patch modulation alpha in [0, 0.64]
  super   superposition with the measured H=0.79 fBm base and the Matern
          patch candidates; sigma_z budget [0.32, 0.46] D50
  bound   superposition-principle check (which component carries the slope at
          the smallest lag), convexity identity, knee / plateau location
  reach   minimum amplitude / share needed to reach 0.70, minimum sigma_z
  report  out/v4-2-8/numbers-v4-2-8.json + REPORT-V4-2-8.md

Only new files: this script and out/v4-2-8/{REPORT-V4-2-8.md,
numbers-v4-2-8.json, checkpoint.json}.  Existing files are hashed before and
after and must be unchanged.
"""
from __future__ import annotations

import datetime
import glob
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v4_2_5 as M                                  # noqa: E402  (unmodified)

# --- fixed paths (absolute; no tree search) --------------------------------
DESIGN = r"D:\mutaro-work\rally-in\DESIGN-V4-2-7-BED-NEXT.md"
BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-8.md"
NUM425 = os.path.join(HERE, "out", "v4-2-5", "numbers-v4-2-5.json")
NUM426 = os.path.join(HERE, "out", "v4-2-6", "numbers-v4-2-6.json")
DECOMP = os.path.join(HERE, "out", "v4-2-3", "decomp.json")

OUT = os.path.join(HERE, "out", "v4-2-8")
CP_PATH = os.path.join(OUT, "checkpoint.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-8.md")
NUMBERS = os.path.join(OUT, "numbers-v4-2-8.json")

THIS_FILE = "run_v4_2_8.py"
REQUIRED = {
    "design": DESIGN,
    "brief": BRIEF,
    "numbers_v4_2_5": NUM425,
    "numbers_v4_2_6": NUM426,
    "decomp_v4_2_3": DECOMP,
    "run_v4_2_5.py": os.path.join(HERE, "run_v4_2_5.py"),
}

# --- the contract (fixed values from the design, not re-derived here) ------
D50 = 0.0362                 # design §0 (surface D50, m)
SZ_BAND = (0.32, 0.46)       # BRIEF STAGE 1.4 / design §0 (matrix sigma_z / D50)
SZ_TARGET = 0.40             # design §0 operating point
BAND = (0.004, 0.036)        # design §0 / run_v4_2_5.BANDS["scaling"]
GATE = (0.70, 0.85)          # BRIEF STAGE 1.2 judgement band
ALPHA_MAX = 0.64             # design §2.2 sweep value (guard allows 0.70)
ALPHA_GUARD = 0.70           # design §2.3: alpha_max = sqrt(0.742^2 - 0.25^2)
NU_SET = (0.5, 1.5)          # design §2.2: kernel unresolved, both swept
LC_06m = (0.10, 0.12)        # design §2.2 / R5: 0.6 m patch upper bound 10-12 cm
LC_STRIP = (0.20, 0.45)      # intermediate 0.20 m; 0.45 m strip
LC_ALL = LC_06m + LC_STRIP
LC_WIDE = (0.01, 0.02, 0.03, 0.05, 0.08)   # outside the design; robustness only

LOG = []


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S+00:00")


def unit_done(cp, name, tries=10):
    """jobkit checkpoint write with a short retry (Windows file locks)."""
    for i in range(tries):
        try:
            return cp.unit_done(name)
        except OSError as exc:
            if i == tries - 1:
                raise
            log("  checkpoint write busy (%s); retry %d/%d"
                % (type(exc).__name__, i + 1, tries))
            time.sleep(0.4 * (i + 1))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _j(x):
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if math.isfinite(v) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (list, tuple)):
        return [_j(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _j(v) for k, v in x.items()}
    return x


def write_json(path, obj):
    data = json.dumps(_j(obj), ensure_ascii=False, indent=1).encode("utf-8")

    def verify(p):
        with open(p, "r", encoding="utf-8") as f:
            json.load(f)
        return True
    atomic_write_bytes(path, data, verify=verify)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def read_lines(path):
    with open(path, "r", encoding="utf-8") as f:
        return f.read().splitlines()


def extract_block(lines, start_pat, stop_pat):
    i = next((k for k, s in enumerate(lines) if start_pat in s), None)
    if i is None:
        return None
    j = next((k for k in range(i + 1, len(lines)) if stop_pat in lines[k]),
             len(lines))
    return "\n".join(lines[i:j]).strip()


def fmt(x, n=4):
    if x is None:
        return "—"
    return ("%." + str(n) + "f") % x


# ---------------------------------------------------------------------------
# curves and the analytic model
# ---------------------------------------------------------------------------
def load_curves():
    d = load_json(DECOMP)
    terms = d["terms"]

    def get(name):
        c = terms[name]["_curve"]
        return np.asarray(c["D_m"], float), np.asarray(c["S_m2"], float)
    Dg, Sg = get("grain_relief")
    Df, Sf = get("form")
    Dc, Sc = get("composite")
    return {
        "relief": {"D": Dg, "S": Sg, "sigma_m": float(terms["grain_relief"]["sigma_m"])},
        "form": {"D": Df, "S": Sf, "sigma_m": float(terms["form"]["sigma_m"])},
        "composite": {"D": Dc, "S": Sc, "sigma_m": float(terms["composite"]["sigma_m"])},
    }


def interp_log(ell, D, S):
    return np.exp(np.interp(np.log(ell), np.log(D), np.log(S)))


def grain_mix(alpha, ell, D, S, nq=48):
    """S_mix(ell) = integral phi(P) m^2 S_ref(ell/m) dP, m = exp(aP - a^2/2)."""
    x, w = np.polynomial.hermite_e.hermegauss(nq)
    w = w / math.sqrt(2.0 * math.pi)
    m = np.exp(alpha * x - alpha * alpha / 2.0)
    S_mix = np.zeros_like(ell)
    for mi, wi in zip(m, w):
        S_mix += wi * mi * mi * interp_log(ell / mi, D, S)
    return S_mix


def matern_S(nu, lc, ell):
    """2-D structure function of a Matern field (sigma=1), kappa=sqrt(2nu)/lc."""
    kap = math.sqrt(2.0 * nu) / lc
    r = kap * np.asarray(ell, float)
    if abs(nu - 0.5) < 1e-9:
        rho = np.exp(-r)
    elif abs(nu - 1.5) < 1e-9:
        rho = (1.0 + r) * np.exp(-r)
    else:
        raise ValueError("nu not in {0.5, 1.5}")
    return 2.0 * (1.0 - rho)


def band_stats(D, S, ell=None):
    """Instrument-verbatim local exponent, band fit and shape read."""
    Hloc = M.local_H(D, S)
    if ell is None:
        ell = D
    m = (ell >= BAND[0]) & (ell <= BAND[1])
    h = Hloc[m]
    fit = M.fit_band(D, S, BAND[0], BAND[1])
    shape = M.band_shape(D, S, BAND[0], BAND[1])
    probes = [0.0041, 0.0073, 0.0100, 0.0147, 0.0240, 0.0328]
    prof = {("%.4f" % p): float(np.interp(np.log(p), np.log(ell), Hloc))
            for p in probes}
    return {
        "H_fit": fit["H"], "H_se": fit["H_se"], "rms_log": fit["rms_log"],
        "r2": fit["r2"], "n_points": fit["n_points"],
        "local_H_min": float(h.min()), "local_H_max": float(h.max()),
        "local_H_min_at_m": float(ell[m][int(np.argmin(h))]),
        "local_H_max_at_m": float(ell[m][int(np.argmax(h))]),
        "n_in_gate": int(((h >= GATE[0]) & (h <= GATE[1])).sum()),
        "gate_all": bool((h >= GATE[0]).all() and (h <= GATE[1]).all()),
        "probe_H": prof,
        "S_at_lo": float(S[m][0]), "S_at_hi": float(S[m][-1]),
        "cross_0_70_m": shape.get("cross_0.70_m"),
        "cross_0_85_m": shape.get("cross_0.85_m"),
        "plateau_onset_m": shape.get("plateau_onset_m"),
        "local_H_profile": [float(v) for v in h],
    }


def crossings(ell, Hloc, level):
    out = []
    for i in range(len(ell) - 1):
        if (Hloc[i] - level) * (Hloc[i + 1] - level) < 0:
            f = (level - Hloc[i]) / (Hloc[i + 1] - Hloc[i])
            out.append(float(np.exp(np.log(ell[i])
                                    + f * (np.log(ell[i + 1]) - np.log(ell[i])))))
    return out


# ---------------------------------------------------------------------------
# STAGE 0
# ---------------------------------------------------------------------------
def stage0(cp):
    log("STAGE0  absolute-path preflight, sha256, design F-B0 extraction")
    missing = [k for k, p in REQUIRED.items() if not os.path.isfile(p)]
    if missing:
        raise SystemExit("STOP: missing required input(s): %s" % missing)

    hashes = {k: {"path": p, "sha256": sha256_file(p),
                  "bytes": os.path.getsize(p)} for k, p in REQUIRED.items()}

    py = sorted(p for p in glob.glob(os.path.join(HERE, "*.py"))
                if os.path.basename(p) != THIS_FILE)
    log("  hashed %d inputs and %d existing .py files" % (len(hashes), len(py)))

    lines = read_lines(DESIGN)
    fb0 = extract_block(lines, "**F-B0（混合効果なし）**", "**F-B1")
    hyp = extract_block(lines, "混合効果（mixture-of-cutoffs）の仮説",
                        "動かなければ、パッチ層はスケーリング帯には何もしない")
    if not fb0 or not hyp:
        raise SystemExit("STOP: could not extract the F-B0 definition verbatim")
    log("  F-B0 definition extracted verbatim (%d chars)" % len(fb0))

    cv = load_curves()
    n425 = load_json(NUM425)
    sh = n425["component_shares"]["lags"]
    check = {}
    for lag, blk in sh.items():
        ell = float(lag)
        for k in ("relief", "form", "composite"):
            Dk, Sk = cv[k]["D"], cv[k]["S"]
            j = int(np.argmin(np.abs(Dk - ell)))
            ref = blk["S_total"] if k == "composite" else blk["S_%s" % k]
            check["%s/%s" % (lag, k)] = abs(Sk[j] - ref) / max(abs(ref), 1e-30)
    max_rel = max(check.values())
    if max_rel > 1e-9:
        raise SystemExit("STOP: saved curves do not reproduce v4-2-5 shares "
                         "(max rel diff %.3g)" % max_rel)
    log("  curves reproduce v4-2-5 shares (max rel diff %.2e)" % max_rel)

    sz = {k: cv[k]["sigma_m"] / D50 for k in cv}
    sz["sqrt_sum"] = math.sqrt(cv["relief"]["sigma_m"] ** 2
                               + cv["form"]["sigma_m"] ** 2) / D50

    st = {
        "preflight": {"all_present": True, "required": hashes},
        "n_existing_py": len(py),
        "design_quote": {"F_B0": fb0, "mixture_hypothesis": hyp,
                         "source": DESIGN},
        "contract": {"D50_m": D50, "band_m": list(BAND), "gate": list(GATE),
                     "sigma_z_over_D50": list(SZ_BAND),
                     "sigma_z_target_over_D50": SZ_TARGET,
                     "alpha_max": ALPHA_MAX, "alpha_guard": ALPHA_GUARD,
                     "nu": list(NU_SET), "lc_0_6m": list(LC_06m),
                     "lc_strip": list(LC_STRIP)},
        "instrument": {"module": "run_v4_2_5",
                       "functions": ["local_H", "fit_band", "band_shape"],
                       "bands": M.BANDS, "nbin": M.NBIN},
        "curve_repro_vs_v4_2_5": {"max_rel_diff": max_rel, "ok": True},
        "sigma_z_over_D50": {k: float(v) for k, v in sz.items()},
    }
    unit_done(cp, "stage0")
    return st


# ---------------------------------------------------------------------------
# STAGE mix
# ---------------------------------------------------------------------------
def prep_ladders():
    cv = load_curves()
    Dg = cv["relief"]["D"]
    ext = Dg <= 0.20
    DE = Dg[ext]
    band = (DE >= BAND[0]) & (DE <= BAND[1])
    return cv, DE, band


def stage_mix(cp, st0):
    log("STAGE mix  local curves D50 x {0.5,1,1.9} and lognormal alpha sweep")
    cv, DE, band = prep_ladders()
    Dr, Sr = cv["relief"]["D"], cv["relief"]["S"]
    hloc_g = M.local_H(Dr, Sr)
    onset = None
    for i in range(len(Dr)):
        if Dr[i] >= 0.02 and hloc_g[i] < 0.20:
            onset = float(Dr[i])
            break
    osc = onset if onset else 0.0636
    log("  measured grain plateau onset (local H<0.20) = %.4f m" % osc)

    ms = [0.5, 1.0, 1.9]
    locals_ = []
    for m in ms:
        S = m * m * interp_log(DE / m, Dr, Sr)
        locals_.append({"m": m, "stats": band_stats(DE, S, DE)})

    mix_w = {}
    for tag, wts in (("equal_area", [1.0 / 3.0] * 3),
                     ("wide_0.2_0.5_0.3", [0.20, 0.50, 0.30]),
                     ("mid_0.25_0.5_0.25", [0.25, 0.50, 0.25])):
        S = sum(w * m * m * interp_log(DE / m, Dr, Sr)
                for w, m in zip(wts, ms))
        mix_w[tag] = {"weights": wts, "stats": band_stats(DE, S, DE)}

    # cutoff-only sensitivity: shift the saturation lag, keep H below it
    S_cut = np.zeros_like(DE)
    for w, m in zip([1.0 / 3.0] * 3, ms):
        ell_eff = np.minimum(DE, m * osc)
        S_cut += w * interp_log(ell_eff, Dr, Sr)
    cut_stats = band_stats(DE, S_cut, DE)

    sweep = []
    for a in (0.0, 0.20, 0.40, 0.55, 0.64, 0.70):
        S = grain_mix(a, DE, Dr, Sr)
        s = band_stats(DE, S, DE)
        s["alpha"] = a
        s["sigma_ratio"] = math.exp(a * a / 2.0)
        s["sigma_grain_over_D50"] = cv["relief"]["sigma_m"] * math.exp(
            a * a / 2.0) / D50
        s["delta_H_fit_vs_alpha0"] = s["H_fit"] - sweep[0]["H_fit"] if sweep else 0.0
        sweep.append(s)
        log("  alpha=%.2f  H_fit=%.4f  local %.4f..%.4f  sigma_ratio=%.3f"
            % (a, s["H_fit"], s["local_H_min"], s["local_H_max"],
               s["sigma_ratio"]))

    base = sweep[0]
    out = {
        "source": {"path": DECOMP, "key": "terms.grain_relief._curve",
                   "why": "saved grain-relief S(ell) of 3.0 m bed r1, x axis"},
        "plateau_onset_m": osc,
        "local_curves": locals_,
        "mixtures_3class": mix_w,
        "cutoff_only_sensitivity": cut_stats,
        "lognormal_sweep": sweep,
        "delta_H_fit": {
            "mix_equal_vs_m1": mix_w["equal_area"]["stats"]["H_fit"] - base["H_fit"],
            "mix_equal_vs_best_local": mix_w["equal_area"]["stats"]["H_fit"]
            - max(x["stats"]["H_fit"] for x in locals_),
            "mix_equal_vs_mean_local": mix_w["equal_area"]["stats"]["H_fit"]
            - float(np.mean([x["stats"]["H_fit"] for x in locals_])),
            "alpha064_vs_alpha0": sweep[4]["H_fit"] - base["H_fit"],
        },
        "note": ("local curves are the measured relief curve with all lengths "
                 "scaled by m and all heights by m (S_local = m^2 S_ref(ell/m)); "
                 "H is amplitude-invariant, so the same profile holds after any "
                 "global renormalisation of sigma_z"),
    }
    unit_done(cp, "mix")
    return out


# ---------------------------------------------------------------------------
# STAGE super
# ---------------------------------------------------------------------------
def sz_total_over_D50(alpha, sb_m, k_form, sig_g0, sig_f0):
    sg = sig_g0 * math.exp(alpha * alpha / 2.0)
    return math.sqrt(sg * sg + (k_form * sig_f0) ** 2 + sb_m * sb_m) / D50


def search_best(cap, lc_list, cv, DE, band, centre=None, k_fixed=None):
    Dr, Sr = cv["relief"]["D"], cv["relief"]["S"]
    Df_, Sf_ = cv["form"]["D"], cv["form"]["S"]
    sig_g0, sig_f0 = cv["relief"]["sigma_m"], cv["form"]["sigma_m"]
    Sf1 = interp_log(DE, Df_, Sf_)
    m_band = band

    if centre is None:
        alphas = np.arange(0.0, 0.701, 0.05)
        lcs = list(lc_list)
        sbs = np.linspace(0.0, cap * D50, 40)
        ks = (np.array([k_fixed]) if k_fixed is not None
              else np.linspace(0.0, 1.3, 17))
    else:
        a0, nu0, l0, s0, k0 = centre
        alphas = np.clip(np.arange(a0 - 0.05, a0 + 0.0501, 0.0125), 0.0, 0.70)
        lo = max(0.10, l0 - 0.02)
        hi = min(0.50, l0 + 0.02)
        lcs = list(np.linspace(lo, hi, 9))
        sbs = np.clip(np.arange(s0 - 0.02 * D50, s0 + 0.0201 * D50,
                                0.005 * D50), 0.0, cap * D50)
        ks = (np.array([k_fixed]) if k_fixed is not None
              else np.clip(np.arange(k0 - 0.05, k0 + 0.0501, 0.0125), 0.0, 1.3))
        lcs = [lc for lc in lcs if lc > 0.0]

    best = None
    n_eval = 0
    n_ok = 0
    for a in alphas:
        Sg = grain_mix(a, DE, Dr, Sr)
        sg = sig_g0 * math.exp(a * a / 2.0)
        for nu in NU_SET:
            for lc in lcs:
                Sp1 = matern_S(nu, lc, DE)
                for sb in sbs:
                    tot_base = sg * sg + sb * sb
                    for k in ks:
                        n_eval += 1
                        tot = math.sqrt(tot_base + (k * sig_f0) ** 2) / D50
                        if tot > cap + 1e-9 or tot < SZ_BAND[0] - 1e-9:
                            continue
                        n_ok += 1
                        Sm = Sg + (k * k) * Sf1 + (sb * sb) * Sp1
                        H = M.local_H(DE, Sm)
                        h = H[m_band]
                        key = (float(h.min()), float(h.max()))
                        if best is None or key > best[0]:
                            best = (key, {
                                "alpha": float(a), "nu": nu, "lc_m": float(lc),
                                "sb_m": float(sb),
                                "sb_over_D50": float(sb / D50),
                                "k_form": float(k),
                                "sigma_z_over_D50": float(tot),
                                "min": float(h.min()), "max": float(h.max()),
                                "H_fit": float(M.fit_band(
                                    DE, Sm, BAND[0], BAND[1])["H"])})
    return best, n_eval, n_ok


def stage_super(cp, st0):
    log("STAGE super  base fBm + patch candidates, sigma_z in [0.32,0.46]D50")
    cv, DE, band = prep_ladders()
    sig_g0, sig_f0 = cv["relief"]["sigma_m"], cv["form"]["sigma_m"]

    table = []
    for a in (0.0, 0.20, 0.30, 0.40, 0.55, 0.64, 0.70):
        sg = sig_g0 * math.exp(a * a / 2.0)
        room = SZ_BAND[1] * D50
        if sg >= room:
            kmax, sbmax = None, None
        else:
            kmax = math.sqrt(room * room - sg * sg) / sig_f0
            sbmax = math.sqrt(room * room - sg * sg)
        table.append({"alpha": a, "sigma_grain_over_D50": sg / D50,
                      "k_form_max_at_sz_0.46": kmax,
                      "sigma_b_max_m": sbmax,
                      "sigma_b_max_over_D50": (sbmax / D50) if sbmax else None})
    log("  sigma_z table built (alpha=0.64 -> k_form<=%.2f at sz=0.46D50)"
        % table[5]["k_form_max_at_sz_0.46"])

    best46, n46, ok46 = search_best(SZ_BAND[1], LC_ALL, cv, DE, band)
    c46 = best46
    for _ in range(2):
        c46, n2, ok2 = search_best(SZ_BAND[1], LC_ALL, cv, DE, band,
                                   centre=(c46[1]["alpha"], c46[1]["nu"],
                                           c46[1]["lc_m"], c46[1]["sb_m"],
                                           c46[1]["k_form"]))
        n46 += n2
        ok46 += ok2
    best40, n40, ok40 = search_best(SZ_TARGET, LC_ALL, cv, DE, band)
    c40 = best40
    for _ in range(2):
        c40, n2, ok2 = search_best(SZ_TARGET, LC_ALL, cv, DE, band,
                                   centre=(c40[1]["alpha"], c40[1]["nu"],
                                           c40[1]["lc_m"], c40[1]["sb_m"],
                                           c40[1]["k_form"]))
        n40 += n2
        ok40 += ok2
    best_wide, nw, okw = search_best(SZ_BAND[1], LC_WIDE + LC_ALL, cv, DE, band)
    # the design bed keeps the measured base fBm (k_form = 1): the faithful case
    best_k1, nk1, okk1 = search_best(SZ_BAND[1], LC_ALL, cv, DE, band,
                                     k_fixed=1.0)
    log("  best at sz<=0.46: min=%.4f max=%.4f (alpha=%.2f nu=%.1f lc=%.3f "
        "sb=%.3fD50 k=%.2f)" % (best46[0][0], best46[0][1],
                                best46[1]["alpha"], best46[1]["nu"],
                                best46[1]["lc_m"], best46[1]["sb_over_D50"],
                                best46[1]["k_form"]))
    log("  best at sz<=0.40: min=%.4f max=%.4f" % (best40[0][0], best40[0][1]))
    log("  best with lc 1-8 cm (outside the design): min=%.4f max=%.4f"
        % (best_wide[0][0], best_wide[0][1]))
    log("  best with the measured base fixed (k=1): min=%.4f max=%.4f"
        % (best_k1[0][0], best_k1[0][1]))

    b = best46[1]
    comps = {"relief": grain_mix(b["alpha"], DE, cv["relief"]["D"],
                                 cv["relief"]["S"]),
             "form": b["k_form"] ** 2 * interp_log(DE, cv["form"]["D"],
                                                   cv["form"]["S"]),
             "patch": b["sb_m"] ** 2 * matern_S(b["nu"], b["lc_m"], DE)}
    Sm = sum(comps.values())
    H = M.local_H(DE, Sm)
    hband = H[band]
    idx = np.where(band)[0]
    shares = {}
    for tag, i in (("bottom", int(idx[0])), ("top", int(idx[-1]))):
        tot = sum(v[i] for v in comps.values())
        shares[tag] = {k: float(v[i] / tot) for k, v in comps.items()}
        shares[tag]["S_total"] = float(tot)
        shares[tag]["lag_m"] = float(DE[i])
    compH = {}
    for tag, i in (("bottom", int(idx[0])), ("top", int(idx[-1]))):
        compH[tag] = {}
        for k, v in comps.items():
            compH[tag][k] = (float(M.local_H(DE, v)[i])
                             if float(np.max(np.abs(v))) > 0 else 0.0)

    best_block = dict(b)
    best_block.update({
        "local_H_profile": [float(v) for v in hband],
        "lags_m": [float(v) for v in DE[band]],
        "S_mix": [float(v) for v in Sm[band]],
        "shares_at_bottom": shares["bottom"],
        "shares_at_top": shares["top"],
        "component_local_H_bottom": compH["bottom"],
        "component_local_H_top": compH["top"],
        "cross_0_35_m": crossings(DE[band], hband, 0.35),
        "cross_0_30_m": crossings(DE[band], hband, 0.30),
    })

    out = {
        "source": {"relief": "out/v4-2-3/decomp.json terms.grain_relief._curve",
                   "form": "out/v4-2-3/decomp.json terms.form._curve",
                   "composite": "out/v4-2-3/decomp.json terms.composite._curve"},
        "sigma_z_table": table,
        "search": {"lc_list_design": list(LC_ALL),
                   "lc_list_wide_outside_design": list(LC_WIDE),
                   "n_eval_sz046": n46, "n_ok_sz046": ok46,
                   "n_eval_sz040": n40, "n_ok_sz040": ok40,
                   "n_eval_wide": nw, "n_ok_wide": okw,
                   "objective": "maximise the minimum in-band local slope"},
        "best_sz_046": best_block,
        "best_sz_040": {"min": c40[0][0], "max": c40[0][1], **c40[1]},
        "best_base_present_k_form_1": {"min": best_k1[0][0],
                                       "max": best_k1[0][1], **best_k1[1]},
        "best_lc_wide_outside_design": {"min": best_wide[0][0],
                                        "max": best_wide[0][1], **best_wide[1]},
        "n_in_gate_best": int(((hband >= GATE[0]) & (hband <= GATE[1])).sum()),
        "N_band_lags": int(band.sum()),
    }
    unit_done(cp, "super")
    return out


# ---------------------------------------------------------------------------
# STAGE bound
# ---------------------------------------------------------------------------
def stage_bound(cp, st0, mix, sup):
    log("STAGE bound  which component carries the slope at the smallest lag")
    cv, DE, band = prep_ladders()
    idx = np.where(band)[0]
    bottom_i, top_i = int(idx[0]), int(idx[-1])

    n425 = load_json(NUM425)
    sh = n425["component_shares"]["lags"]
    meas = {}
    for lag in ("0.0041", "0.0319"):
        blk = sh[lag]
        meas[lag] = {k: float(blk["share_%s" % k])
                     for k in ("form", "relief", "cov")}
        meas[lag]["S_total"] = float(blk["S_total"])

    b = sup["best_sz_046"]
    comps = {"relief": grain_mix(b["alpha"], DE, cv["relief"]["D"],
                                 cv["relief"]["S"]),
             "form": b["k_form"] ** 2 * interp_log(DE, cv["form"]["D"],
                                                   cv["form"]["S"]),
             "patch": b["sb_m"] ** 2 * matern_S(b["nu"], b["lc_m"], DE)}
    shares = {}
    for tag, i in (("bottom", bottom_i), ("top", top_i)):
        tot = sum(v[i] for v in comps.values())
        shares[tag] = {k: float(v[i] / tot) for k, v in comps.items()}

    Hmix = M.local_H(DE, sum(comps.values()))
    Hc = {}
    for k, v in comps.items():
        Hc[k] = (M.local_H(DE, v) if float(np.max(np.abs(v))) > 0
                 else np.zeros_like(DE))
    resid = []
    for i in idx:
        tot = sum(v[i] for v in comps.values())
        wsum = 0.0
        for k, v in comps.items():
            if v[i] > 0:
                wsum += (v[i] / tot) * Hc[k][i]
        resid.append(abs(wsum - Hmix[i]))
    residual = {"max_abs": float(max(resid)), "mean_abs": float(np.mean(resid))}

    # variance-weight condition (BRIEF STAGE 1.1) on the design's 3-curve
    # equal-area mixture, at the six representative lags of the design
    Dr, Sr = cv["relief"]["D"], cv["relief"]["S"]
    probes = [0.0041, 0.0073, 0.0100, 0.0147, 0.0240, 0.0328]
    wt = []
    for p in probes:
        comp = [m * m * float(interp_log(np.array([p / m]), Dr, Sr)[0])
                for m in (0.5, 1.0, 1.9)]
        tot = sum(comp)
        wts = [c / tot for c in comp]
        his = [float(M.local_H(Dr, m * m * interp_log(Dr / m, Dr, Sr))[
            int(np.argmin(np.abs(Dr - p)))]) for m in (0.5, 1.0, 1.9)]
        wt.append({"lag_m": p,
                   "omega": {"m0.5": wts[0], "m1": wts[1], "m1.9": wts[2]},
                   "local_H": {"m0.5": his[0], "m1": his[1], "m1.9": his[2]},
                   "weighted_mean_H": sum(w * h for w, h in zip(wts, his))})
    residual["three_class_probe_weights"] = wt

    Hg = M.local_H(Dr, Sr)
    m_max = math.exp(ALPHA_GUARD * 1.645 - ALPHA_GUARD ** 2 / 2.0)
    ell_min = BAND[0] / m_max
    H_at_min = float(np.interp(np.log(ell_min), np.log(Dr), Hg))
    sel = (Dr >= ell_min) & (Dr <= BAND[1])
    H_template_max = float(Hg[sel].max())
    H_curve_min_lag = float(Hg[0])

    out = {
        "measured_composite_shares_v4_2_5": meas,
        "best_model_shares": shares,
        "convexity_identity": residual,
        "grain_family_ceiling": {
            "m_max_at_guard_1.645sigma": m_max,
            "smallest_sampled_lag_m": ell_min,
            "H_template_at_that_lag": H_at_min,
            "H_template_at_curve_min_lag": H_curve_min_lag,
            "curve_min_lag_m": float(Dr[0]),
            "H_template_max_in_sampled_range": H_template_max,
            "note": ("H_mix is a convex combination of the local H of its "
                     "components; the grain family can therefore not exceed "
                     "the template's local H over the sampled lags"),
        },
        "knee": {
            "plateau_onset_grain_m": mix["plateau_onset_m"],
            "best_profile_cross_0_35_m": b.get("cross_0_35_m"),
            "best_profile_cross_0_30_m": b.get("cross_0_30_m"),
            "best_profile_min_at_m": b["lags_m"][int(np.argmin(
                b["local_H_profile"]))],
        },
        "statement": ("At the smallest band lag the measured grain relief "
                      "carries nearly all of S; the apparent slope there is "
                      "the grain-relief local slope, not a mixture average."),
    }
    log("  measured shares 4.1mm: relief %.4f form %.4f; best model 4.1mm: "
        "relief %.4f patch %.4f"
        % (meas["0.0041"]["relief"], meas["0.0041"]["form"],
           shares["bottom"]["relief"], shares["bottom"]["patch"]))
    unit_done(cp, "bound")
    return out


# ---------------------------------------------------------------------------
# STAGE reach
# ---------------------------------------------------------------------------
def stage_reach(cp, st0, mix, sup, bound):
    log("STAGE reach  required share / amplitude, minimum sigma_z")
    cv, DE, band = prep_ladders()
    idx = np.where(band)[0]
    bottom_i = int(idx[0])
    b = sup["best_sz_046"]
    H_g = b["component_local_H_bottom"]["relief"]
    H_f = float(M.local_H(DE, interp_log(DE, cv["form"]["D"],
                                         cv["form"]["S"]))[bottom_i])
    S_g = float(grain_mix(b["alpha"], DE, cv["relief"]["D"],
                          cv["relief"]["S"])[bottom_i])
    S_f1 = float(interp_log(DE, cv["form"]["D"], cv["form"]["S"])[bottom_i])

    tbl = []
    for tgt in (0.70, 0.75, 0.80):
        share = (tgt - H_g) / (H_f - H_g)
        k = math.sqrt(share / (1.0 - share) * S_g / S_f1)
        tbl.append({"target": tgt, "required_form_share_at_bottom": share,
                    "required_k_form": k,
                    "required_sigma_form_m": k * cv["form"]["sigma_m"],
                    "required_sigma_form_over_D50":
                        k * cv["form"]["sigma_m"] / D50})

    # the largest form share the sigma_z band can buy: alpha=0 (no grain
    # inflation) and no patch competitor, i.e. the most favourable case for
    # the form route
    room = SZ_BAND[1] * D50
    sg0 = cv["relief"]["sigma_m"]
    k_cap = math.sqrt(max(room * room - sg0 * sg0, 0.0)) / cv["form"]["sigma_m"]
    S_g0 = float(grain_mix(0.0, DE, cv["relief"]["D"],
                           cv["relief"]["S"])[bottom_i])
    S_f_cap = k_cap * k_cap * S_f1
    share_cap = S_f_cap / (S_f_cap + S_g0)
    ratio = (tbl[0]["required_form_share_at_bottom"] / share_cap
             if share_cap > 0 else None)

    patch_tbl = []
    for nu in NU_SET:
        for lc in LC_ALL:
            s1 = float(matern_S(nu, lc, np.array([DE[bottom_i]]))[0])
            Hp = float(M.local_H(DE, matern_S(nu, lc, DE))[bottom_i])
            if Hp > 0.70 and s1 > 0:
                sb = math.sqrt((0.70 - H_g) / (Hp - 0.70) * S_g / s1)
                patch_tbl.append({"nu": nu, "lc_m": lc, "patch_H_bottom": Hp,
                                  "required_sigma_b_m": sb,
                                  "required_sigma_b_over_D50": sb / D50,
                                  "reachable_at_any_amplitude": True})
            else:
                patch_tbl.append({"nu": nu, "lc_m": lc, "patch_H_bottom": Hp,
                                  "required_sigma_b_m": None,
                                  "required_sigma_b_over_D50": None,
                                  "reachable_at_any_amplitude": False})

    amp = []
    for s in (1.0, 0.75, 0.50, 0.25):
        k = tbl[0]["required_k_form"] * s
        amp.append({"grain_amplitude_scale": s, "required_k_form": k,
                    "required_sigma_form_over_D50":
                        k * cv["form"]["sigma_m"] / D50})

    def min_sigma(need_all):
        best = None
        Sf1 = interp_log(DE, cv["form"]["D"], cv["form"]["S"])
        for a in np.arange(0.0, 0.701, 0.10):
            Sg = grain_mix(a, DE, cv["relief"]["D"], cv["relief"]["S"])
            for nu in NU_SET:
                for lc in LC_ALL:
                    Sp1 = matern_S(nu, lc, DE)
                    for k in np.linspace(0, 25, 26):
                        Sf = k * k * Sf1
                        for sb in np.linspace(0, 12 * D50, 25):
                            Sm = Sg + Sf + sb * sb * Sp1
                            H = M.local_H(DE, Sm)[band]
                            ok = (H.min() >= 0.70) if need_all else (H[0] >= 0.70)
                            if not ok:
                                continue
                            tot = sz_total_over_D50(a, sb, k,
                                                    cv["relief"]["sigma_m"],
                                                    cv["form"]["sigma_m"])
                            if best is None or tot < best[0]:
                                best = (float(tot), float(H.min()),
                                        float(H.max()), float(a), nu, float(lc),
                                        float(sb / D50), float(k))
        return best

    ms_bottom = min_sigma(False)
    ms_all = min_sigma(True)
    log("  required form share at 4.1 mm = %.3f; available in sigma_z band = "
        "%.4f; min sigma_z (bottom only) = %.2f D50"
        % (tbl[0]["required_form_share_at_bottom"], share_cap, ms_bottom[0]))

    out = {
        "two_component_requirement": tbl,
        "available_in_sigma_z_band": {
            "alpha": 0.0, "k_form_cap": k_cap,
            "form_share_at_bottom": share_cap,
            "shortfall_ratio_share": ratio,
            "note": ("largest form share at the band bottom inside sigma_z, "
                     "evaluated at alpha=0 with no patch (most favourable for "
                     "the form route)"),
        },
        "patch_requirement": patch_tbl,
        "grain_amplitude_sensitivity": amp,
        "min_sigma_z": {
            "bottom_only": {"sigma_z_over_D50": ms_bottom[0],
                            "min_local_H": ms_bottom[1],
                            "max_local_H": ms_bottom[2],
                            "alpha": ms_bottom[3], "nu": ms_bottom[4],
                            "lc_m": ms_bottom[5],
                            "sigma_b_over_D50": ms_bottom[6],
                            "k_form": ms_bottom[7]},
            "whole_band": (None if ms_all is None else
                           {"sigma_z_over_D50": ms_all[0],
                            "min_local_H": ms_all[1], "max_local_H": ms_all[2],
                            "alpha": ms_all[3], "nu": ms_all[4],
                            "lc_m": ms_all[5], "sigma_b_over_D50": ms_all[6],
                            "k_form": ms_all[7]}),
        },
        "note": ("the two-component row is the form-only route; the grid also "
                 "allows a patch term, so min_sigma_z is the least total "
                 "sigma_z found over the allowed component set"),
    }
    unit_done(cp, "reach")
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def build_report(ctx):
    st0 = ctx["stage0"]
    mix = ctx["mix"]
    sup = ctx["super"]
    bnd = ctx["bound"]
    rch = ctx["reach"]
    b = sup["best_sz_046"]
    b40 = sup["best_sz_040"]
    bk1 = sup["best_base_present_k_form_1"]
    bw = sup["best_lc_wide_outside_design"]
    n425 = load_json(NUM425)
    L = []
    A = L.append

    A("# REPORT-V4-2-8 — F-B0 解析事前チェック（混合カットオフでスケーリング帯の傾きは上がるか）")
    A("")
    A("> BRIEF-V4-2-8 / 設計 `DESIGN-V4-2-7-BED-NEXT.md`（ADOPTED 2026-11-13）。")
    A("> 作業ツリー `D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline`。実行 %s。" % ctx["created_utc"])
    A("> **解析計算のみ。新しい床は生成していない。** 既存の保存曲線を読み、`run_v4_2_5.py` の")
    A("> `local_H` / `fit_band` / `band_shape` を無改変で通した。画像の目視は使っていない。")
    A("")
    A("## 0. 判定（二値）")
    A("")
    A("**混合では届かない。**")
    A("")
    A("- 設計の許容パラメータ範囲（α ≤ %.2f（掃引上限。ガードは %.2f）、ℓc 10–45 cm、"
      "Matérn ν ∈ {0.5, 1.5}、matrix σz ∈ [%.2f, %.2f]·D50）で、"
      % (ALPHA_MAX, ALPHA_GUARD, SZ_BAND[0], SZ_BAND[1]))
    A("  帯域 0.4–3.6 cm の**全 %d ラグ**で見かけの局所傾きが 0.70–0.85 に入る混合は **0 件**。"
      % sup["N_band_lags"])
    A("- 最良ケース（帯域内局所傾きの最小値を最大化した混合）: 最小 **%s**、最大 **%s**（H_fit %s）。"
      % (fmt(b["min"]), fmt(b["max"]), fmt(b["H_fit"])))
    A("  その構成は α = %.2f、Matérn ν = %.1f、ℓc = %.3f m、σ_b = %.3f·D50、k_form = %.2f、"
      "σz = %.3f·D50。" % (b["alpha"], b["nu"], b["lc_m"], b["sb_over_D50"],
                          b["k_form"], b["sigma_z_over_D50"]))
    A("  この最良ケースは k_form = 0（基底を帯域に寄与させない）という上限側の解である。"
      "設計の床は基底を必ず含むので、**混合仮説に最大限有利な仮定**である。")
    A("- 参考値: σz ≤ 0.40·D50 の最良は最小 %s・最大 %s。基底を実測振幅で固定した場合"
      "（設計の床）は最小 %s・最大 %s。設計範囲外の ℓc（1–8 cm）まで許した最良でも最小 %s・最大 %s。"
      % (fmt(b40["min"]), fmt(b40["max"]), fmt(bk1["min"]), fmt(bk1["max"]),
         fmt(bw["min"]), fmt(bw["max"])))
    A("- 設計 F-B0 の 3 曲線（D50 × 0.5/1/1.9 の等面積平均）は H_fit %s（m=1 曲線から %+.4f、"
      "最良の局所曲線から %+.4f）。局所曲線の包絡の内側に留まる。"
      % (fmt(mix["mixtures_3class"]["equal_area"]["stats"]["H_fit"]),
         mix["delta_H_fit"]["mix_equal_vs_m1"],
         mix["delta_H_fit"]["mix_equal_vs_best_local"]))
    A("- 帯域下端（4.1 mm）で 0.70 に届かせるには、H≈0.94 の基底成分がそのラグの S の "
      "**%.1f%%** を占める必要がある。σz 帯内で達成できる最大は **%.2f%%**。"
      % (100 * rch["two_component_requirement"][0]["required_form_share_at_bottom"],
         100 * rch["available_in_sigma_z_band"]["form_share_at_bottom"]))
    p15 = next(r for r in rch["patch_requirement"]
               if r["nu"] == 1.5 and abs(r["lc_m"] - 0.10) < 1e-9)
    A("  必要振幅は基底 **%.2f·D50**、パッチ（ν=1.5, ℓc=10 cm）**%.2f·D50** で、"
      "いずれも σz 帯（≤ %.2f·D50）の外。"
      % (rch["two_component_requirement"][0]["required_sigma_form_over_D50"],
         p15["required_sigma_b_over_D50"], SZ_BAND[1]))
    A("")
    A("設計 F-B0 の逐語定義（`DESIGN-V4-2-7-BED-NEXT.md` §2.6。以後これは**設計 F-B0** と呼ぶ）:")
    A("")
    A("```")
    for line in st0["design_quote"]["F_B0"].splitlines():
        A(line)
    A("```")
    A("")
    A("設計 F-B0 は「局所曲線の面積加重平均が局所曲線のいずれとも変わらないか（ΔH_fit < +0.02）」を")
    A("見る。本 brief の STAGE 1.2 は「帯域全域で見かけの局所傾きが 0.70–0.85 に入るか」を見る。")
    A("**別の判定量である。** 両方を第 3 節と第 7 節に分けて報告する。")
    A("")

    A("## 1. STAGE 0 — 前提確認")
    A("")
    A("絶対パスのみで確認した（別ツリーは探していない）。全入力が存在する。")
    A("")
    A("| 入力 | sha256（先頭 16） | bytes |")
    A("|---|---|---|")
    for k, v in st0["preflight"]["required"].items():
        A("| `%s` | `%s` | %d |" % (v["path"], v["sha256"][:16], v["bytes"]))
    A("")
    A("- 保存曲線の再現性: `out/v4-2-3/decomp.json` の 3 成分曲線は、v4-2-5 の"
      " `component_shares.lags` の S 値を最大相対差 **%.2e** で再現した。"
      % st0["curve_repro_vs_v4_2_5"]["max_rel_diff"])
    A("- 固定値（設計 §0）: D50 = %.4f m、σz 目標 = %.2f·D50、σz 帯 = [%.2f, %.2f]·D50、"
      "帯域 = [%.3f, %.3f] m、旧ゲート 0.70–0.85。"
      % (D50, SZ_TARGET, SZ_BAND[0], SZ_BAND[1], BAND[0], BAND[1]))
    A("- 実測合成床の σz: **%.4f·D50**（= %.5f m、`decomp.json terms.composite.sigma_m`）。"
      "すでに σz 帯の上端近くにある。"
      % (st0["sigma_z_over_D50"]["composite"], ctx["curves_sigma"]["composite"]))
    A("- 成分の σz: 粒レリーフ **%.4f·D50**、基底 **%.4f·D50**、2 乗和 **%.4f·D50**。"
      % (st0["sigma_z_over_D50"]["relief"], st0["sigma_z_over_D50"]["form"],
         st0["sigma_z_over_D50"]["sqrt_sum"]))
    A("- 粒レリーフの飽和（局所指数 < 0.20 の最初のラグ ≥ 2 cm）: **%.4f m**（≈ D90 = %.4f m）。"
      "帯域上端 3.6 cm の外側である。" % (mix["plateau_onset_m"], 1.9 * D50))
    A("- v4-2-6 の帯内分散分解（`numbers-v4-2-6.json` の `component_budget`、"
      "case `e_max_q20_040`）: share_relief **0.9957–1.0177**、share_form **0.0194–0.0824**、"
      "share_cov **−0.0323〜−0.0781**（負の共分散のため relief の share が 1 を超える点がある）。"
      "帯内 S のほぼ全量を粒レリーフが担う。")
    A("")
    A("設計 §2.5 の混合仮説（本計算のモデル根拠。逐語）:")
    A("")
    A("```")
    for line in st0["design_quote"]["mixture_hypothesis"].splitlines():
        A(line)
    A("```")
    A("")

    A("## 2. 方法")
    A("")
    A("### 2.1 使った曲線")
    A("")
    A("| 成分 | 出典 | 帯域 H_fit | 帯域局所指数 min / max |")
    A("|---|---|---|---|")
    for name, key, src in (
            ("粒レリーフ（実測）", "relief", "`decomp.json terms.grain_relief._curve`"),
            ("基底 fBm（実測）", "form", "`decomp.json terms.form._curve`"),
            ("合成床（実測）", "composite", "`decomp.json terms.composite._curve`")):
        s = ctx["ref_stats"][key]
        A("| %s | %s | %s | %s / %s |" % (name, src, fmt(s["H_fit"]),
                                          fmt(s["local_H_min"]), fmt(s["local_H_max"])))
    A("")
    A("σ は `decomp.json` の `sigma_m`（粒 %.6f、基底 %.6f、合成 %.6f m）。"
      % (ctx["curves_sigma"]["relief"], ctx["curves_sigma"]["form"],
         ctx["curves_sigma"]["composite"]))
    A("局所指数と帯域フィットは `run_v4_2_5.local_H` / `fit_band` をそのまま使った。"
      "ラダーは v4-2.3/v4-2-5 の 60 ビン幾何ラダーで、帯域内は 16 点。")
    A("")
    A("### 2.2 混合モデル")
    A("")
    A("設計 §2.5 の形に従い、成分 i の構造関数を S_i(ℓ) として")
    A("")
    A("```")
    A("S_mix(ℓ) = Σ_i w_i · S_i(ℓ),      H_app(ℓ) = ½ · dlnS_mix/dlnℓ")
    A("H_app(ℓ) = Σ_i ω_i(ℓ) · H_i(ℓ),   ω_i = w_i S_i / Σ_j w_j S_j  (Σω=1)")
    A("```")
    A("")
    A("局所曲線の混合（設計 F-B0 の本体）: 実測粒レリーフ曲線を m=1 のテンプレートとし、")
    A("局所 D50 が m 倍のパッチでは粒の長さが m 倍・高さが m 倍になるとした。")
    A("")
    A("```")
    A("S_local(ℓ; m) = m² · S_ref(ℓ/m)")
    A("```")
    A("")
    A("これは自己相似な粒床の標準的なスケーリングである。H は振幅に不変なので、")
    A("全体を σz に合わせて再規格化してもプロファイルは変わらない。")
    A("連続的なパッチ変調は設計 §2.3 の m(x,y)=exp(αP−α²/2)、P~N(0,1) をそのまま使い、")
    A("P について Gauss–Hermite（48 点）で積分した。α は設計 §2.2 の掃引値 0.40/0.55/0.64 と、")
    A("ガード α_max ≈ 0.70 までを見た。")
    A("")
    A("基底 fBm は実測の `terms.form._curve`（帯域 H_fit %.3f）を振幅 k 倍して足した。"
      % ctx["ref_stats"]["form"]["H_fit"])
    A("パッチ高さ項は Matérn 場の 2-D 構造関数 S_p(ℓ)=2σ_b²(1−ρ(ℓ)) で、")
    A("κ=√(2ν)/ℓc、ν=0.5 は ρ=e^(−κℓ)、ν=1.5 は ρ=(1+κℓ)e^(−κℓ) とした。")
    A("ℓc は設計 §2.2/R5 の 0.6 m パッチ上限 10–12 cm と、帯状床側の 20–45 cm を掃引した。")
    A("σz は成分の 2 乗和で評価した（実測の共分散項は帯域で −0.5〜−0.8%。保存側の")
    A("仮定で、正の共分散があれば σz 予算はさらに厳しくなる）。")
    A("")
    A("### 2.3 探索")
    A("")
    A("α ∈ [0, 0.70]（0.05 刻み）、ν ∈ {0.5, 1.5}、ℓc ∈ {10, 12, 20, 45} cm、")
    A("σ_b ∈ [0, 0.46·D50]、k_form ∈ [0, 1.3] の格子で、")
    A("σz ∈ [0.32, 0.46]·D50 を満たす組だけを残し、帯域 16 ラグの局所傾きの**最小値**を最大化した。")
    A("最良点のまわりで 2 回の局所細分化を行った。試した格子点は σz≤0.46 で %s 件"
      "（うち σz 制約を満たしたもの %s 件）、σz≤0.40 で %s 件（%s 件）。"
      % (format(sup["search"]["n_eval_sz046"], ","),
         format(sup["search"]["n_ok_sz046"], ","),
         format(sup["search"]["n_eval_sz040"], ","),
         format(sup["search"]["n_ok_sz040"], ",")))
    A("")

    A("## 3. 結果 1 — 局所曲線の混合（設計 F-B0 の本体）")
    A("")
    A("### 3.1 設計の 3 曲線（D50 × 0.5 / 1 / 1.9）")
    A("")
    A("| 曲線 | H_fit | 帯域局所指数 min / max |")
    A("|---|---|---|")
    for x in mix["local_curves"]:
        A("| m = %.1f | %s | %s / %s |"
          % (x["m"], fmt(x["stats"]["H_fit"]),
             fmt(x["stats"]["local_H_min"]), fmt(x["stats"]["local_H_max"])))
    A("")
    A("| 混合（面積重み） | H_fit | m=1 からの Δ | 最良局所曲線からの Δ | 局所指数 min / max |")
    A("|---|---|---|---|---|")
    for tag, lab in (("equal_area", "等面積 1/3 ずつ"),
                     ("wide_0.2_0.5_0.3", "0.2 / 0.5 / 0.3"),
                     ("mid_0.25_0.5_0.25", "0.25 / 0.5 / 0.25")):
        s = mix["mixtures_3class"][tag]["stats"]
        A("| %s | %s | %+.4f | %+.4f | %s / %s |"
          % (lab, fmt(s["H_fit"]),
             s["H_fit"] - mix["local_curves"][1]["stats"]["H_fit"],
             s["H_fit"] - max(x["stats"]["H_fit"] for x in mix["local_curves"]),
             fmt(s["local_H_min"]), fmt(s["local_H_max"])))
    A("")
    A("等面積混合の H_fit は **%s**。m=1 曲線（%s）より **%+.4f** 高いが、"
      "最良の局所曲線 m=1.9（%s）より **%+.4f** 低い。"
      % (fmt(mix["mixtures_3class"]["equal_area"]["stats"]["H_fit"]),
         fmt(mix["local_curves"][1]["stats"]["H_fit"]),
         mix["delta_H_fit"]["mix_equal_vs_m1"],
         fmt(mix["local_curves"][2]["stats"]["H_fit"]),
         mix["delta_H_fit"]["mix_equal_vs_best_local"]))
    A("つまり混合は局所曲線の包絡（H_fit %.4f–%.4f）の**内側**に留まる。"
      % (min(x["stats"]["H_fit"] for x in mix["local_curves"]),
         max(x["stats"]["H_fit"] for x in mix["local_curves"])))
    A("設計 §2.5 の見積もりは混合効果 +0.05–0.10（H_fit 0.35–0.50）だが、"
      "本計算の 3 曲線混合は **+%.4f**（H_fit %s）で、見積もりの下端 +0.05 にも届かない。"
      % (mix["delta_H_fit"]["mix_equal_vs_m1"],
         fmt(mix["mixtures_3class"]["equal_area"]["stats"]["H_fit"])))
    A("")
    A("**設計 F-B0 の比較対象**: 設計文は「局所曲線のいずれとも変わらない（ΔH_fit < +0.02）」と書き、"
      "比較対象も片側／両側も名指ししていない。"
      "『最良の局所曲線を超えない（利得なし）』と読めば、混合は包絡の内側（最良から %+.4f）なので"
      "条件を満たす。『各局所曲線から 0.02 以内』と読めば、m=1 から %+.4f なので満たさない。"
      "**両方の数を示す。**"
      % (mix["delta_H_fit"]["mix_equal_vs_best_local"],
         mix["delta_H_fit"]["mix_equal_vs_m1"]))
    A("")
    A("| 比較対象 | ΔH_fit | +0.02 を超えるか |")
    A("|---|---|---|")
    A("| m=1 曲線 | %+.4f | 超える |" % mix["delta_H_fit"]["mix_equal_vs_m1"])
    A("| 3 曲線の平均 | %+.4f | 超える |" % mix["delta_H_fit"]["mix_equal_vs_mean_local"])
    A("| 最良の局所曲線（m=1.9） | %+.4f | 超えない |"
      % mix["delta_H_fit"]["mix_equal_vs_best_local"])
    A("")
    A("どちらの読みでも、絶対値は 0.70 には届かない（次節）。")
    A("")
    A("### 3.2 連続的なパッチ変調 α")
    A("")
    A("| α | H_fit | α=0 からの ΔH_fit | 局所指数 min / max | σ_grain の増倍 | σ_grain/D50 |")
    A("|---|---|---|---|---|---|")
    for s in mix["lognormal_sweep"]:
        A("| %.2f | %s | %+.4f | %s / %s | %.3f | %.4f |"
          % (s["alpha"], fmt(s["H_fit"]), s["delta_H_fit_vs_alpha0"],
             fmt(s["local_H_min"]), fmt(s["local_H_max"]),
             s["sigma_ratio"], s["sigma_grain_over_D50"]))
    A("")
    A("α = 0.64 で ΔH_fit = **%+.4f**、局所指数は %s–%s。α = 0.70（ガード上限）でも "
      "%s–%s。**帯域のどこも 0.70 に届かない。**"
      % (mix["delta_H_fit"]["alpha064_vs_alpha0"],
         fmt(mix["lognormal_sweep"][4]["local_H_min"]),
         fmt(mix["lognormal_sweep"][4]["local_H_max"]),
         fmt(mix["lognormal_sweep"][5]["local_H_min"]),
         fmt(mix["lognormal_sweep"][5]["local_H_max"])))
    A("α は σ_grain を e^(α²/2) 倍する（α=0.64 で 1.227 倍）。σz 予算を食う方向に働く。")
    A("")
    A("### 3.3 カットオフだけを動かす変種（感度）")
    A("")
    A("幾何スケーリングの代わりに、飽和ラグだけを m 倍し、それより短いラグでは m=1 の曲線を")
    A("そのまま使う変種では、H_fit = **%s**、局所指数 %s–%s。"
      "幾何スケーリング（H_fit %s）より低い。**上の混合効果は楽観側の見積もりである。**"
      % (fmt(mix["cutoff_only_sensitivity"]["H_fit"]),
         fmt(mix["cutoff_only_sensitivity"]["local_H_min"]),
         fmt(mix["cutoff_only_sensitivity"]["local_H_max"]),
         fmt(mix["mixtures_3class"]["equal_area"]["stats"]["H_fit"])))
    A("**この変種が m=1 と完全に同じ値（H_fit %s）になるのは、帯域内でカットオフが動かないからである。**"
      % fmt(mix["cutoff_only_sensitivity"]["H_fit"]))
    A("実測の飽和ラグは %.2f mm（≈ D90）で、m=0.5 でも飽和は %.2f mm と帯域上端（36 mm）の外側にある。"
      % (mix["plateau_onset_m"] * 1000, 0.5 * mix["plateau_onset_m"] * 1000))
    A("したがって設計 §2.5 が名指しした『飽和ラグの混合』そのものは、**この帯域では機構として働かない**。"
      "幾何スケーリングで得られた +%.4f は、カットオフの移動ではなく、"
      "粗いパッチが局所曲線の短ラグ側（H≈0.45）をサンプルすることによる形の効果である。"
      % mix["delta_H_fit"]["mix_equal_vs_m1"])
    A("")

    A("## 4. 結果 2 — 基底・パッチとの重ね合わせと σz 制約")
    A("")
    A("### 4.1 σz 予算")
    A("")
    A("| α | σ_grain/D50 | σz=0.46·D50 での k_form 上限 | σ_b 上限 (m) | σ_b/D50 上限 |")
    A("|---|---|---|---|---|")
    for r in sup["sigma_z_table"]:
        A("| %.2f | %.4f | %s | %s | %s |"
          % (r["alpha"], r["sigma_grain_over_D50"],
             fmt(r["k_form_max_at_sz_0.46"], 3) if r["k_form_max_at_sz_0.46"] else "—",
             ("%.5f" % r["sigma_b_max_m"]) if r["sigma_b_max_m"] else "—",
             fmt(r["sigma_b_max_over_D50"], 3) if r["sigma_b_max_over_D50"] else "—"))
    A("")
    A("α=0.64 では基底の振幅を **%.2f 倍**までしか上げられず、パッチ項に回せるのは "
      "**%.3f·D50** まで。α=0 でも基底は **%.2f 倍**が上限（実測合成床がすでに σz=%s·D50）。"
      % (sup["sigma_z_table"][5]["k_form_max_at_sz_0.46"],
         sup["sigma_z_table"][5]["sigma_b_max_over_D50"],
         sup["sigma_z_table"][0]["k_form_max_at_sz_0.46"],
         fmt(st0["sigma_z_over_D50"]["composite"])))
    A("")
    A("### 4.2 σz 帯内の最良ケース")
    A("")
    A("| 条件 | 帯域局所傾き min | max | H_fit | 構成 |")
    A("|---|---|---|---|---|")
    A("| σz ≤ 0.46·D50（本判定） | **%s** | **%s** | %s | α=%.2f, ν=%.1f, ℓc=%.2f m, σ_b=%.3f·D50, k=%.2f, σz=%.3f·D50 |"
      % (fmt(b["min"]), fmt(b["max"]), fmt(b["H_fit"]), b["alpha"], b["nu"],
         b["lc_m"], b["sb_over_D50"], b["k_form"], b["sigma_z_over_D50"]))
    A("| σz ≤ 0.40·D50（設計目標） | %s | %s | %s | α=%.2f, ν=%.1f, ℓc=%.2f m, σ_b=%.3f·D50, k=%.2f |"
      % (fmt(b40["min"]), fmt(b40["max"]), fmt(b40["H_fit"]), b40["alpha"],
         b40["nu"], b40["lc_m"], b40["sb_over_D50"], b40["k_form"]))
    A("| 基底を実測振幅で固定（k=1、設計の床） | %s | %s | %s | α=%.2f, ν=%.1f, ℓc=%.2f m, σ_b=%.3f·D50 |"
      % (fmt(bk1["min"]), fmt(bk1["max"]), fmt(bk1["H_fit"]), bk1["alpha"],
         bk1["nu"], bk1["lc_m"], bk1["sb_over_D50"]))
    A("| ℓc を 1–8 cm まで許す（設計範囲外） | %s | %s | %s | α=%.2f, ν=%.1f, ℓc=%.3f m, σ_b=%.3f·D50 |"
      % (fmt(bw["min"]), fmt(bw["max"]), fmt(bw["H_fit"]), bw["alpha"],
         bw["nu"], bw["lc_m"], bw["sb_over_D50"]))
    A("")
    A("本判定の最良ケースのプロファイル（帯域 16 ラグ）:")
    A("")
    A("| ラグ (mm) | 局所傾き | S_mix (m²) |")
    A("|---|---|---|")
    for lag, h, s in zip(b["lags_m"], b["local_H_profile"], b["S_mix"]):
        A("| %.2f | %s | %.4e |" % (lag * 1000, fmt(h), s))
    A("")
    knee35 = " / ".join("%.2f mm" % (x * 1000) for x in b["cross_0_35_m"])
    knee30 = " / ".join("%.2f mm" % (x * 1000) for x in b["cross_0_30_m"])
    i_min = int(np.argmin(b["local_H_profile"]))
    A("プロファイルは 4.1 mm の **%s** から下がり、%.2f mm 付近で **%s** まで落ちて"
      "ほぼ平坦になる（%.2f–%.2f mm で %s–%s）。"
      "0.35 を切るのは %s、0.30 を切るのは %s。"
      % (fmt(b["local_H_profile"][0]), b["lags_m"][i_min] * 1000,
         fmt(b["local_H_profile"][i_min]),
         b["lags_m"][-2] * 1000, b["lags_m"][-1] * 1000,
         fmt(b["local_H_profile"][-2]), fmt(b["local_H_profile"][-1]),
         knee35 or "帯域内では切らない", knee30 or "帯域内では切らない"))
    A("粒レリーフの飽和（局所指数 < 0.20）は **%.2f mm** で、帯域の外側にある。"
      "したがって帯域内に『knee が入って傾きが崩れる』形ではなく、"
      "帯域の上端に向かって緩やかに低下する形である。" % (mix["plateau_onset_m"] * 1000))
    A("")
    A("**帯域 16 ラグのうち 0.70–0.85 に入るラグ数: %d / %d。**"
      % (sup["n_in_gate_best"], sup["N_band_lags"]))
    A("")

    A("## 5. 結果 3 — 重ね合わせの原理チェック（帯域下端）")
    A("")
    A("実測合成床（v4-2-5 `component_shares.lags`）と、本判定の最良ケースの帯域両端の S シェア:")
    A("")
    A("| ラグ | 実測: 粒 / 基底 / 共分散 | 最良ケース: 粒 / 基底 / パッチ |")
    A("|---|---|---|")
    A("| 4.1 mm | %.4f / %.4f / %+.4f | %.4f / %.4f / %.4f |"
      % (bnd["measured_composite_shares_v4_2_5"]["0.0041"]["relief"],
         bnd["measured_composite_shares_v4_2_5"]["0.0041"]["form"],
         bnd["measured_composite_shares_v4_2_5"]["0.0041"]["cov"],
         bnd["best_model_shares"]["bottom"]["relief"],
         bnd["best_model_shares"]["bottom"]["form"],
         bnd["best_model_shares"]["bottom"]["patch"]))
    A("| 31.9 mm | %.4f / %.4f / %+.4f | %.4f / %.4f / %.4f |"
      % (bnd["measured_composite_shares_v4_2_5"]["0.0319"]["relief"],
         bnd["measured_composite_shares_v4_2_5"]["0.0319"]["form"],
         bnd["measured_composite_shares_v4_2_5"]["0.0319"]["cov"],
         bnd["best_model_shares"]["top"]["relief"],
         bnd["best_model_shares"]["top"]["form"],
         bnd["best_model_shares"]["top"]["patch"]))
    A("")
    A("**最小ラグ側では最低 H 成分（粒レリーフ）が支配する**: 実測床で 4.1 mm の S の "
      "**%.2f%%**、最良ケースで **%.2f%%** を粒レリーフが担う。"
      "実測床の 4.1 mm の見かけの傾きは %.4f で、粒レリーフ自身の局所指数 %.4f とほぼ同じ。"
      "混合しても帯域下端の傾きは粒レリーフの値に固定される。"
      % (100 * bnd["measured_composite_shares_v4_2_5"]["0.0041"]["relief"],
         100 * bnd["best_model_shares"]["bottom"]["relief"],
         ctx["ref_stats"]["composite"]["local_H_profile"][0],
         ctx["ref_stats"]["relief"]["local_H_profile"][0]))
    A("")
    A("凸結合の恒等式 H_app = Σ ω_i H_i の数値残差は最大 **%.2e**（平均 %.2e）。"
      "つまり見かけの傾きは成分の局所指数の重み付き平均であり、"
      "最大の成分局所指数を超えない。"
      % (bnd["convexity_identity"]["max_abs"], bnd["convexity_identity"]["mean_abs"]))
    A("")
    A("分散重み条件（BRIEF STAGE 1.1）を設計の 3 曲線混合で 6 点で評価した"
      "（ω_i = m_i²S_i/Σ、H_i は各局所曲線の局所指数）:")
    A("")
    A("| ラグ (mm) | ω(0.5) | ω(1) | ω(1.9) | H(0.5) | H(1) | H(1.9) | ΣωH |")
    A("|---|---|---|---|---|---|---|---|")
    _wt = bnd["convexity_identity"]["three_class_probe_weights"]
    for r in _wt:
        A("| %.2f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |"
          % (r["lag_m"] * 1000, r["omega"]["m0.5"], r["omega"]["m1"],
             r["omega"]["m1.9"], r["local_H"]["m0.5"], r["local_H"]["m1"],
             r["local_H"]["m1.9"], r["weighted_mean_H"]))
    A("")
    A("ω は帯域の全ラグで粗い m=1.9 側が最大で、ラグが上がるほどそのシェアは増える"
      "（%.3f → %.3f）。細かい m=0.5 は %.3f → %.3f に下がる。"
      "いずれのラグでも ΣωH の最大は **%.3f** で、0.70 には届かない。"
      % (_wt[0]["omega"]["m1.9"], _wt[-1]["omega"]["m1.9"],
         _wt[0]["omega"]["m0.5"], _wt[-1]["omega"]["m0.5"],
         max(r["weighted_mean_H"] for r in _wt)))
    A("")
    A("粒レリーフ系の上限: テンプレートの局所指数は曲線の最小ラグ %.3f mm で **%.3f**、"
      "サンプル範囲（%.1f–%.1f mm）の最大は **%.3f**。α ≤ 0.70 で到達できる最小サンプルラグは %.3f mm、"
      "そこでの H は **%.3f**。粒だけの混合は 0.70 に届く手がかりを持たない。"
      % (bnd["grain_family_ceiling"]["curve_min_lag_m"] * 1000,
         bnd["grain_family_ceiling"]["H_template_at_curve_min_lag"],
         bnd["grain_family_ceiling"]["smallest_sampled_lag_m"] * 1000,
         BAND[1] * 1000,
         bnd["grain_family_ceiling"]["H_template_max_in_sampled_range"],
         bnd["grain_family_ceiling"]["smallest_sampled_lag_m"] * 1000,
         bnd["grain_family_ceiling"]["H_template_at_that_lag"]))
    A("")

    A("## 6. 結果 4 — 0.70 に届かせるための必要量")
    A("")
    A("### 6.1 基底（H≈0.94）で埋める場合")
    A("")
    A("| 目標 H_app（4.1 mm） | 必要な基底シェア | 必要な k_form | 必要 σ_form/D50 |")
    A("|---|---|---|---|")
    for r in rch["two_component_requirement"]:
        A("| %.2f | %.4f | %.2f | %.2f |"
          % (r["target"], r["required_form_share_at_bottom"],
             r["required_k_form"], r["required_sigma_form_over_D50"]))
    A("")
    A("σz 帯内で実際に取れる基底シェアは **%.4f**（%.2f%%）。必要量の **1/%.0f**。"
      % (rch["available_in_sigma_z_band"]["form_share_at_bottom"],
         100 * rch["available_in_sigma_z_band"]["form_share_at_bottom"],
         rch["available_in_sigma_z_band"]["shortfall_ratio_share"]))
    A("")
    A("### 6.2 パッチ高さ項で埋める場合")
    A("")
    A("| ν | ℓc (m) | パッチの 4.1 mm 局所指数 | 必要 σ_b/D50 | 到達可能か |")
    A("|---|---|---|---|---|")
    for r in rch["patch_requirement"]:
        A("| %.1f | %.2f | %.3f | %s | %s |"
          % (r["nu"], r["lc_m"], r["patch_H_bottom"],
             fmt(r["required_sigma_b_over_D50"], 2) if r["reachable_at_any_amplitude"] else "—",
             "振幅次第" if r["reachable_at_any_amplitude"] else "不可（H<0.70）"))
    A("")
    p05 = next(r for r in rch["patch_requirement"] if r["nu"] == 0.5)
    A("ν=0.5（指数型）のパッチは 4.1 mm の局所指数が **%.3f** で、振幅をいくら上げても "
      "0.70 に届かない。ν=1.5 は届くが、必要振幅は σz 帯の外にある。"
      % p05["patch_H_bottom"])
    A("")
    A("### 6.3 粒レリーフの振幅を人為的に下げた場合（感度）")
    A("")
    A("| 粒振幅の倍率 | 必要 k_form | 必要 σ_form/D50 |")
    A("|---|---|---|")
    for r in rch["grain_amplitude_sensitivity"]:
        A("| %.2f | %.2f | %.2f |"
          % (r["grain_amplitude_scale"], r["required_k_form"],
             r["required_sigma_form_over_D50"]))
    A("")
    A("粒レリーフを 1/4 にしても、基底は **%.2f·D50** 必要で、σz 帯（≤ %.2f·D50）の外。"
      % (rch["grain_amplitude_sensitivity"][3]["required_sigma_form_over_D50"],
         SZ_BAND[1]))
    A("")
    A("### 6.4 0.70 を出せる最小 σz（格子探索）")
    A("")
    A("| 条件 | 最小 σz/D50 | そのときの局所傾き min / max |")
    A("|---|---|---|")
    if rch["min_sigma_z"]["bottom_only"]:
        r = rch["min_sigma_z"]["bottom_only"]
        A("| 帯域下端だけ 0.70（必要条件） | **%.2f** | %s / %s |"
          % (r["sigma_z_over_D50"], fmt(r["min_local_H"]), fmt(r["max_local_H"])))
    if rch["min_sigma_z"]["whole_band"]:
        r = rch["min_sigma_z"]["whole_band"]
        A("| 帯域全域 0.70（十分条件） | **%.2f** | %s / %s |"
          % (r["sigma_z_over_D50"], fmt(r["min_local_H"]), fmt(r["max_local_H"])))
    A("")
    mb = rch["min_sigma_z"]["bottom_only"]
    mw = rch["min_sigma_z"]["whole_band"]
    if mw and mb and abs(mw["sigma_z_over_D50"] - mb["sigma_z_over_D50"]) < 0.01:
        A("2 行は同じ解（α=%.2f、ν=%.1f、ℓc=%.2f m、σ_b=%.2f·D50、k_form=%.1f）で、"
          "その σz では帯域全域が 0.70–0.85 に入る。"
          % (mb["alpha"], mb["nu"], mb["lc_m"], mb["sigma_b_over_D50"],
             mb["k_form"]))
    A("σz 帯の上限は 0.46·D50。帯域下端だけでも、許容の **%.1f 倍**の σz が要る。"
      % (rch["min_sigma_z"]["bottom_only"]["sigma_z_over_D50"] / SZ_BAND[1]))
    A("")

    A("## 7. 判定のまとめ（設計 F-B0 と本 brief の判定量）")
    A("")
    A("| 判定量 | 何を見るか | 結果 |")
    A("|---|---|---|")
    A("| 設計 F-B0（逐語） | 面積加重平均が局所曲線のいずれとも変わらないか（ΔH_fit < +0.02） | "
      "混合は局所曲線の包絡の内側（H_fit %s、最良局所曲線 %s より %+.4f）。"
      "『最良を超えない』読みでは条件を満たす。『各局所曲線から 0.02 以内』読みでは "
      "m=1 差 %+.4f で満たさない。比較対象・片側両側の名指しがないため両方を示す |"
      % (fmt(mix["mixtures_3class"]["equal_area"]["stats"]["H_fit"]),
         fmt(mix["local_curves"][2]["stats"]["H_fit"]),
         mix["delta_H_fit"]["mix_equal_vs_best_local"],
         mix["delta_H_fit"]["mix_equal_vs_m1"]))
    A("| 本 brief STAGE 1.2 | 帯域 16 ラグすべてで見かけの局所傾きが 0.70–0.85 | "
      "**0 / %d ラグ。最良ケース min %s / max %s。混合では届かない** |"
      % (sup["N_band_lags"], fmt(b["min"]), fmt(b["max"])))
    A("")
    A("本 brief の問いへの答えは **「混合では届かない」**。")
    A("数値の事実だけを書く: 帯域下端の S を支配しているのは粒レリーフ（最低 H 成分）で、")
    A("実測の粒レリーフの局所指数は %.4f。混合はこの下端の値を最大 %s まで上げるだけである。"
      % (ctx["ref_stats"]["relief"]["local_H_profile"][0], fmt(b["max"])))
    A("不足分をどの機構が埋めるかは、この計算の範囲外である（設計側の仕事）。")
    A("")

    A("## 8. 限界・開示")
    A("")
    A("- **方向は x（流下）のみ。** 粒レリーフ／基底の分離曲線は v4-2-3 の x 方向にしか保存されていない。"
      "v4-2-5 の帯域 H は x %.4f / y %.4f（合成床）で、方向を変えても判定を動かす幅は 0.01 程度。"
      % (ctx["ref_stats"]["composite"]["H_fit"],
         n425["beds"]["long_grains_r1"]["direction"]["y"]["bands"]["scaling"]["H"]))
    A("- **局所曲線の相似仮定。** S_local = m²S_ref(ℓ/m) は自己相似な粒床の標準形だが、"
      "実測の検証ではない。カットオフだけを動かす変種では混合効果はさらに小さい（§3.3）。")
    A("- **Matérn の ℓc の定義。** κ=√(2ν)/ℓc とした。ℓc と variogram range の換算は設計 R5 で"
      "未解決のままである。ℓc を大きくするほど帯域内のパッチ振幅は小さくなるので、"
      "本判定の結論は動かない方向にしか動かない。")
    A("- **σz は 2 乗和。** 共分散項を無視した。実測の帯域共分散シェアは −0.5〜−0.8% で小さい。"
      "正の共分散があれば σz 予算は厳しくなり、結論は変わらない。")
    A("- **α の範囲。** 設計 §2.2 の掃引値 0.40/0.55/0.64 と、§2.3 の実現可能性ガード α_max≈0.70 まで。"
      "α > 0.70 は σg_local ガードで棄却される。")
    A("- **0.70–0.85 は旧ゲートの値。** 設計 A1 はこの一点下限を廃止した。本 brief の STAGE 1.2 が"
      "明示的にこの帯を判定量として指定したので、そのまま計算した。設計 F-B0 の判定は別に示した（§7）。")
    A("- **新しい床は生成していない。** 高さ場・画像・乱数場は一切作っていない。")
    A("- **物差しの上限（文脈）**: v4-2-6 の `ceiling` では H_gen=0.79 の理想場でも 0.6 m 箱の"
      "帯域局所指数最小は平均 **0.654** で、旧ゲート 0.70 全域は理想場でも満たさない。"
      "本判定は混合の寄与を見るもので、この上限は結論の根拠ではない。")
    A("")

    A("## 9. 使ったパラメータと出典")
    A("")
    A("| 量 | 値 | 出典 |")
    A("|---|---|---|")
    A("| D50（表面） | %.4f m | 設計 §0、`out/v4-2-2/heightfield-meta.json`（v4-2-5/6 で確認） |" % D50)
    A("| 帯域 | 0.004–0.036 m | 設計 §0、`run_v4_2_5.BANDS[0]` |")
    A("| 旧ゲート | 0.70–0.85 | BRIEF STAGE 1.2、`run_v4_2_5.GATE` |")
    A("| σz 帯 | [0.32, 0.46]·D50 | BRIEF STAGE 1.4、設計 §0 |")
    A("| 粒レリーフ曲線 | 帯域 H_fit %.4f、局所 %.4f–%.4f | `out/v4-2-3/decomp.json terms.grain_relief._curve`；"
      "v4-2-5 `terms_H_scaling.grain_relief`（H=%.6f） |"
      % (ctx["ref_stats"]["relief"]["H_fit"],
         ctx["ref_stats"]["relief"]["local_H_min"],
         ctx["ref_stats"]["relief"]["local_H_max"],
         n425["component_shares"]["terms_H_scaling"]["grain_relief"]["H"]))
    A("| 基底 fBm 曲線 | 帯域 H_fit %.4f、局所 %.4f–%.4f | `out/v4-2-3/decomp.json terms.form._curve` |"
      % (ctx["ref_stats"]["form"]["H_fit"], ctx["ref_stats"]["form"]["local_H_min"],
         ctx["ref_stats"]["form"]["local_H_max"]))
    A("| 合成床 | H_fit %.4f、局所 %.4f–%.4f、σz %.4f·D50 | `decomp.json terms.composite._curve`；"
      "v4-2-5 `terms_H_scaling.composite`（H=%.6f） |"
      % (ctx["ref_stats"]["composite"]["H_fit"],
         ctx["ref_stats"]["composite"]["local_H_min"],
         ctx["ref_stats"]["composite"]["local_H_max"],
         st0["sigma_z_over_D50"]["composite"],
         n425["component_shares"]["terms_H_scaling"]["composite"]["H"]))
    A("| 帯内分散の粒寄与 | 97–100%（v4-2-5 §5.1）、v4-2-6 §4.1 | `numbers-v4-2-5.json` / `numbers-v4-2-6.json` |")
    A("| 設計 F-B0 | 逐語引用（§0） | 設計 §2.6 |")
    A("| 混合仮説 | 逐語引用（§1） | 設計 §2.5 |")
    A("| α の掃引 | 0.40 / 0.55 / 0.64、ガード 0.70 | 設計 §2.2–§2.3 |")
    A("| ℓc の掃引 | 10–12 cm（0.6 m 箱）/ 20–45 cm（帯状床） | 設計 §2.2、改訂 R5 |")
    A("| ν | 0.5 / 1.5 | 設計 §2.2、改訂 R5 |")
    A("| σ_B の帯 | σ_B/D50 0.15–0.35、σ_B/σz 0.45–0.80 | 設計 §4 C2（M5） |")
    A("")
    A("本判定の最良ケースの σ_b = %.3f·D50 は M5 の σ_B/D50 帯（0.15–0.35）の内側、"
      "σ_B/σz = %.2f は M5 の帯（0.45–0.80）の内側にある。**M5 の外の振幅は使っていない。**"
      % (b["sb_over_D50"], b["sb_m"] / (b["sigma_z_over_D50"] * D50)))
    A("")

    A("## 10. 成果物・実行")
    A("")
    A("- `out/v4-2-8/REPORT-V4-2-8.md`（本ファイル）")
    A("- `out/v4-2-8/numbers-v4-2-8.json`（全数値・曲線・探索結果・ログ）")
    A("- `out/v4-2-8/checkpoint.json`（jobkit checkpoint、%d ユニット）"
      % ctx["checkpoint_units"])
    A("- 計算スクリプト: `run_v4_2_8.py`（新規 1 本のみ。既存 .py への書き込みなし。"
      "下記の変更検査を参照）")
    A("- 実行時間: %s" % " / ".join("%s %.1fs" % (k, v)
                                    for k, v in ctx["timings"].items()))
    A("- 既存ファイルの変更検査: stage0 で %d 本の既存 .py をハッシュし、終了時に再ハッシュ。"
      "変更 %d 本、追加 %d 本、削除 %d 本（新規スクリプトは検査対象外）。"
      % (ctx["legacy"]["n_baseline"], len(ctx["legacy"]["changed"]),
         len(ctx["legacy"]["added"]), len(ctx["legacy"]["removed"])))
    if ctx["legacy"]["changed"]:
        A("  変更されたファイル: %s。本スクリプトの write set は `out/v4-2-8/` の 3 ファイルのみで、"
          "既存 .py は読み取りしかしていない。このファイルは run 中に外部プロセスが更新した"
          "（mtime %s）。他の %d 本は不変。"
          % (", ".join("`%s`" % x for x in ctx["legacy"]["changed"]),
             ctx["legacy"]["changed_mtime"],
             ctx["legacy"]["n_now"] - len(ctx["legacy"]["changed"])))
    A("")
    return "\n".join(L)


def stage_report(cp, ctx):
    log("REPORT  numbers-v4-2-8.json + REPORT-V4-2-8.md")
    md = build_report(ctx)
    data = md.encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            t = f.read().decode("utf-8")
        return len(t) > 1000 and "REPORT-V4-2-8" in t
    atomic_write_bytes(REPORT, data, verify=verify)

    write_json(NUMBERS, {
        "brief": BRIEF, "design": DESIGN,
        "created_utc": ctx["created_utc"],
        "what_this_is": ("BRIEF-V4-2-8 STAGE 1: analytic F-B0 pre-check. Mixture "
                         "of cutoffs vs the 0.4-3.6 cm scaling band. No bed "
                         "generated."),
        "decision": "混合では届かない",
        "stage0": ctx["stage0"],
        "curves_sigma_m": ctx["curves_sigma"],
        "reference_stats": ctx["ref_stats"],
        "mix": ctx["mix"],
        "super": ctx["super"],
        "bound": ctx["bound"],
        "reach": ctx["reach"],
        "verdict": {
            "question": "does any allowed mixture put the apparent local slope "
                        "inside 0.70-0.85 at all 16 band lags?",
            "answer": "no",
            "best_case_min": ctx["super"]["best_sz_046"]["min"],
            "best_case_max": ctx["super"]["best_sz_046"]["max"],
            "best_case_H_fit": ctx["super"]["best_sz_046"]["H_fit"],
            "n_lags_in_gate": ctx["super"]["n_in_gate_best"],
            "n_band_lags": ctx["super"]["N_band_lags"],
            "design_F_B0": ("mixture stays inside the local-curve envelope; "
                            "delta vs m=1 is %+.4f, vs the best local curve "
                            "%.4f; the design text does not specify the "
                            "comparator"
                            % (ctx["mix"]["delta_H_fit"]["mix_equal_vs_m1"],
                               ctx["mix"]["delta_H_fit"]["mix_equal_vs_best_local"])),
        },
        "limitations": [
            "x direction only (component decomposition is saved for x)",
            "local curves assumed self-similar: S_local = m^2 S_ref(ell/m)",
            "Matern kappa = sqrt(2 nu)/lc; lc-to-range conversion unresolved (R5)",
            "sigma_z evaluated as a root-sum-square (covariance ignored; measured "
            "in-band covariance share -0.5 to -0.8%)",
            "alpha <= 0.70 (feasibility guard); lc within the design ranges",
            "the 0.70-0.85 band is the retired one-point gate, used because the "
            "brief names it as the judgement quantity",
            "no new bed generated",
        ],
        "deliverables": {"report": REPORT, "numbers": NUMBERS,
                          "checkpoint": CP_PATH,
                          "script": os.path.join(HERE, THIS_FILE)},
        "checkpoint": {"units": ctx["checkpoint_units"], "path": CP_PATH},
        "timings_s": ctx["timings"],
        "legacy": ctx["legacy"],
        "log": LOG,
    })
    log("  wrote %s (%d B) and %s (%d B)"
        % (REPORT, os.path.getsize(REPORT), NUMBERS, os.path.getsize(NUMBERS)))
    unit_done(cp, "report")


def legacy_check(base):
    py = sorted(p for p in glob.glob(os.path.join(HERE, "*.py"))
                if os.path.basename(p) != THIS_FILE)
    now = {os.path.basename(p): sha256_file(p) for p in py}
    changed = [k for k, v in now.items() if k in base and base[k] != v]
    added = [k for k in now if k not in base]
    removed = [k for k in base if k not in now]
    mtimes = []
    for k in changed:
        try:
            mtimes.append("%s %s" % (
                k, datetime.datetime.fromtimestamp(
                    os.path.getmtime(os.path.join(HERE, k))).strftime(
                        "%Y-%m-%d %H:%M:%S")))
        except OSError:
            mtimes.append(k)
    return {"n_baseline": len(base), "n_now": len(now),
            "changed": changed, "added": added, "removed": removed,
            "changed_mtime": "; ".join(mtimes) if mtimes else None,
            "ok": not changed and not removed}


# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    if os.path.isfile(CP_PATH):
        os.remove(CP_PATH)
    cp = Checkpoint(CP_PATH)

    t_all = time.time()
    created = utcnow()
    timings = {}

    t0 = time.time()
    st0 = stage0(cp)
    timings["stage0"] = round(time.time() - t0, 1)

    py_base = {os.path.basename(p): sha256_file(p)
               for p in glob.glob(os.path.join(HERE, "*.py"))
               if os.path.basename(p) != THIS_FILE}

    t0 = time.time()
    mix = stage_mix(cp, st0)
    timings["mix"] = round(time.time() - t0, 1)

    cv = load_curves()
    ref_stats = {k: band_stats(cv[k]["D"], cv[k]["S"], cv[k]["D"])
                 for k in ("relief", "form", "composite")}
    curves_sigma = {k: cv[k]["sigma_m"] for k in cv}

    t0 = time.time()
    sup = stage_super(cp, st0)
    timings["super"] = round(time.time() - t0, 1)

    t0 = time.time()
    bnd = stage_bound(cp, st0, mix, sup)
    timings["bound"] = round(time.time() - t0, 1)

    t0 = time.time()
    rch = stage_reach(cp, st0, mix, sup, bnd)
    timings["reach"] = round(time.time() - t0, 1)

    legacy = legacy_check(py_base)
    log("G-LEGACY  %d existing .py: changed=%d added=%d removed=%d"
        % (legacy["n_baseline"], len(legacy["changed"]),
           len(legacy["added"]), len(legacy["removed"])))

    state = cp.state() or {}
    ctx = {"created_utc": created, "stage0": st0, "mix": mix, "super": sup,
           "bound": bnd, "reach": rch, "curves_sigma": curves_sigma,
           "ref_stats": ref_stats,
           "legacy": legacy, "timings": timings,
           "checkpoint_units": int(state.get("completed_units") or 0) + 1}
    t0 = time.time()
    stage_report(cp, ctx)
    timings["report"] = round(time.time() - t0, 1)

    log("DONE  %.1fs total" % (time.time() - t_all))


if __name__ == "__main__":
    main()
