"""Post-hoc diagnostics for v2 gate analysis (record-only, no gate logic).

Writes out/v2/numbers-v2-diag.json:
  D1: photon-vs-photon MS-SSIM (same condition, other seed) = noise ceiling.
  D2: G-band-only det-vs-photon MS-SSIM (removes chromatic-dispersion gap).
  D3: det-vs-photon MS-SSIM at 1024 (both area-meaned; noise suppression).
  D4: quadrature tail fractions for Q9v25 / Q25v49.
"""
import io
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env import resolve_environment  # noqa: E402
from state import WaterState  # noqa: E402
from transport import photon_frame  # noqa: E402
from image import e_base_dir_from_env  # noqa: E402
from run_v2 import (ms_ssim, area_mean_downsample, spectral_band_frac,  # noqa: E402
                    quantiles, to_lum, REC709)

OUT = "out/v2"


def main():
    t0 = time.time()
    with io.open("../blender/recipe-1.json", encoding="utf-8") as f:
        recipe = json.load(f)
    env = resolve_environment(recipe)
    state = WaterState(recipe, env, N=2048, seed=None, preset_index=1)
    fields = state.fields_at(0.0)
    fields["S"] = state.S
    spectral_band_frac.S = state.S
    Ebase = e_base_dir_from_env(env)
    L_base = float((Ebase * REC709).sum())

    C_det = np.load(os.path.join(OUT, "C-t0.npy")).astype(np.float64)
    gam = 1.0

    _, info1 = photon_frame(fields, state.S, 2048, env, 4000000,
                            seed=(state.seed + 1000) % (2 ** 63))
    E1, _ = photon_frame(fields, state.S, 2048, env, 4000000,
                         seed=(state.seed + 1000) % (2 ** 63))
    # NOTE: photon_frame has no seed-reset subtlety; rerun reproduces E1.
    E2, _ = photon_frame(fields, state.S, 2048, env, 4000000,
                         seed=(state.seed + 2000) % (2 ** 63))
    C1 = to_lum(E1) / L_base
    C2 = to_lum(E2) / L_base
    X1, X2 = np.log1p(gam * C1), np.log1p(gam * C2)
    d1 = {"ms_photon_vs_photon": ms_ssim(X1, X2),
          "quantiles_ph1": quantiles(C1), "quantiles_ph2": quantiles(C2)}

    Cg = (E1[:, :, 1] / float(Ebase[1])).astype(np.float64)
    d2 = {"ms_det_vs_photon_Gband": ms_ssim(np.log1p(gam * C_det),
                                            np.log1p(gam * Cg))}

    A1024 = area_mean_downsample(C_det, 2)
    B1024 = area_mean_downsample(C1, 2)
    d3 = {"ms_det_vs_photon_at1024": ms_ssim(np.log1p(gam * A1024),
                                             np.log1p(gam * B1024))}

    # quadrature tails (reload Q fields? recompute cheap proxy: use saved
    # C-t0 (Q49). Q9/Q25 not saved -> recompute not needed: tails were the
    # blocker; recompute Q9/Q25 quickly? ~56s. Do it for rigor.)
    from det import det_frame  # noqa: E402
    C9, _ = det_frame(fields, state.S, 2048, env, 9)
    C25, _ = det_frame(fields, state.S, 2048, env, 25)
    C49 = C_det
    tails = {}
    for name, A, B in (("Q9v25", C9.astype(np.float64), C25.astype(np.float64)),
                       ("Q25v49", C25.astype(np.float64), C49)):
        d = np.abs(np.log(np.clip(A, 1e-12, None))
                   - np.log(np.clip(B, 1e-12, None)))
        tails[name] = {"frac_gt_0.02": float((d > 0.02).mean()),
                       "frac_gt_0.10": float((d > 0.10).mean()),
                       "frac_gt_0.50": float((d > 0.50).mean()),
                       "p99_abs": float(np.percentile(d, 99)),
                       "p999_abs": float(np.percentile(d, 99.9))}
    diag = {"D1_noise_ceiling": d1, "D2_Gband": d2, "D3_at1024": d3,
            "D4_quadrature_tails": tails,
            "note": ("D1 reruns the check photon field (seed+1000) plus an "
                     "independent realization (seed+2000) at identical "
                     "condition (4M photons, same pixel filter)."),
            "wall_s": time.time() - t0}
    with io.open(os.path.join(OUT, "numbers-v2-diag.json"), "w",
                  encoding="utf-8") as f:
        json.dump(diag, f, indent=2, ensure_ascii=False)
    print("D1 ceiling msssim=%.4f" % d1["ms_photon_vs_photon"]["msssim"])
    print("D2 Gband  msssim=%.4f" % d2["ms_det_vs_photon_Gband"]["msssim"])
    print("D3 1024   msssim=%.4f" % d3["ms_det_vs_photon_at1024"]["msssim"])
    print("D4 tails:", json.dumps(tails))
    print("OK wall=%.1fs" % diag["wall_s"])


if __name__ == "__main__":
    main()
