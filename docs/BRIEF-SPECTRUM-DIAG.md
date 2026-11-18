# BRIEF — spectrum diagnosis: where do sigma^2 and H_s actually live?

Measurement only — do NOT change any model code. PC,
`D:\yohey2026\dsh\muse\caustics-filter\offline\` (v8 stack).
Model: deepseek-v4.1-flash. Foreground, jobkit checkpoints,
numerical output only.

## Background (measured facts, do not re-derive)
At wind 6.0 / fetch 100 m / z = 1.0 the v3.4 pipeline realised
sigma^2 = 0.0556 (anchor target met) but H_s = 0.332 m, while the
SMB fetch law for the same inputs gives H_s = 0.031 m — a 10x
disagreement. Hypothesis under test: the current spectral shape
places too much slope variance at the gravity peak and lacks the
capillary tail that, in real water, carries slope variance with
negligible height.

## Sweep
For each (wind, fetch) in {1.5, 3, 6, 10} x {100, 1000} (z = 1.0,
seed fixed to the production seed), build the v8 water state and
report:
- anchor sigma^2 target vs realised sigma^2
- integrated H_s from the spectrum vs the SMB fetch-law H_s
  (report the ratio)
- BAND DECOMPOSITION of both m0 (height variance) and sigma^2
  (slope variance) over wavelength bands: >1 m, 0.1-1 m,
  1-10 cm, <1 cm — each band's share (%) of m0 and of sigma^2
- lambda_p, patch size S, grid dx, Nyquist wavelength, and whether
  the capillary peak (k ~ 370 rad/m, lambda ~ 1.7 cm) is resolved
  by the grid at that patch size
- the Stokes/depth caps status

## Outputs (leave in `offline\out\spectrum-diag\`)
1. `spectrum-diag.json` — the full table
2. `REPORT-SPECTRUM-DIAG.md` — the table in readable form + answers
   to: (a) which band carries sigma^2 in the current shape, and
   which carries m0; (b) at which (wind, fetch) the H_s ratio vs
   SMB exceeds 2x; (c) what sigma^2 the SMB-consistent H_s would
   imply if the shape were kept (i.e. how much slope the current
   shape loses when height is forced right) — the size of the
   missing capillary contribution, per cell of the grid.
Coordinator fetches; do not push to the VM.
