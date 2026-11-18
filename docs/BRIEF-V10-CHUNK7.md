# BRIEF — chunk 7 = version 3.7 implementation (two-regime spectrum)

Design: `D:\mutaro-work\rally-in\PIPELINE-V9.md` (READ IT FIRST, in
full — it is the specification; this brief only adds mechanics).
Diagnosis evidence: PC `offline/out/spectrum-diag/` (report + JSON).
Working dir: `D:\yohey2026\dsh\muse\caustics-filter\offline`.
Model for this session: deepseek-v4.1-flash (project rule).

## Rules of the session
- Run FOREGROUND inside this session, stage by stage. Do NOT end the
  session until stage 4's deliverables are written. Use jobkit
  (D:\mutaro-work\jobkit) checkpoints per stage; atomic writes only.
- NEVER taskkill python.exe / blender.exe (parallel work + a LAN
  server share this PC). No see_image — numerical statistics only.
- ASCII in all output files. New code in the v10 family:
  v10_spectrum.py, v10_state.py, run_v10.py (do not edit v5–v9 model
  files; import from them where the design says "existing").
- If a gate fails for a PHYSICAL reason, record it openly in the
  report with the numbers (precedent: chunk 3.5 gate 3). If G9 fails
  by quadrature shredding, STOP at stage 3 and write the failure up —
  do not improvise a remedy.

## Stage 1 — spectrum module (design §2, §3)
Implement v10_spectrum.py: JONSWAP S_J(ω; H_s, f_p) from the fetch
laws (reuse the law family already in v8_physics where present);
conversion S_g(k) = S_J(ω(k))·dω/dk with the existing dispersion;
Elfouhaily B_h(k; u*) converted into the same S(k) convention and
added (S = S_g + S_h); ECKV Φ(k,θ); analytic integrals (m0, σ²,
band shares, per-part overlap 50–200 rad/m, anisotropy prediction).
Checks: G0 Parseval against a built test field; G4 positivity over
the 8-cell grid (U ∈ {1.5,3,6,10} × F ∈ {100,1000}, z=1.0) with the
low-wind α_m floor applied EXPLICITLY and logged (design §3).
Output: out/v10/stage1-spectrum.json + checkpoint.

## Stage 2 — state builder + sweep (design §4, §5, §6)
Implement v10_state.py: scene-fixed patch (recipe property);
two-cascade realisation of the single spectrum with a partition of
unity in ln k; both cascades quantise ω to the same 2π/8 s grid,
shared seed family; single-cascade mode when one grid provably
serves both ends (record the mode). waveHeightM is NOT used.
Sweep runner over the 8 cells producing the gate table:
G1 (H_s realised spectral+spatial vs law), G2 (σ² vs own integral;
band shares; fractions above k=100/370), G3 (σ² beside Cox–Munk
clean / retired 0.0548 / full-ECKV), G5 (anisotropy vs CM table),
G6 (k_cut sensitivity {800,1000,1200,1500}), G7 (per-band max
|Δω|/ω from loop quantisation), plus per-cell depth-correction
magnitude + t_min + hs_source + eckv_domain flags.
Output: out/v10/numbers-v10-sweep.json + checkpoint.

## Stage 3 — transport hookup + regressions (design §7 G8, G9)
Run the existing transport/composite on v10 states for the two
reference scenes (preset-1 family scene; the windy recipe below).
G8 battery: frozen composite on identical buffers bit-identical;
Beer–Lambert outputs for identical optical inputs unchanged vs the
v9 stack; η=1 golden still on the preset-1 scene compared against
the v9 η=1 route on ITS state (report the delta and its cause —
the state changed, the route must not have); rain R=0 equals the
v10 base state exactly.
G9: photon cross-check + C percentiles (mean, p50, p99, p99.9, max)
on both scenes; deliverable-still MS-SSIM vs photon render.
Output: out/v10/numbers-v10-transport.json + stills + checkpoint.

## Stage 4 — windy re-render + report (design §7 G10)
Re-render the windy-pond recipe under v3.7: moth photo
(D:\mutaro-work\moth-on-bark.jpg, native 735×539), wind 6.0 m/s,
fetch 100 m, depth 1.0 m, rain 0, η=1.0, sun 2026-06-10 09:30
(alt 59.03°, az 105.2°), frozen composite — same as
BRIEF-WINDY-STILL.md except the state generator is v10.
Measure the warp yardsticks on the result (block-matching
displacement vs the unwarped photo: 48 px blocks, ±24 px search,
2 px steps, luminance; AND the pipeline's own geometric shift
field, mean in px and in mm at the scene scale) and tabulate
against: v3.4 windy still (mean 10.50 px NCC / own-field 5.84 px
at 7.63 mm/px = 48.5 mm) and the production-loop yardstick
(23.2 px mean / 17.1 px std at ~1.77 mm/px).
Write out/v10/REPORT-V10-GATES.md: every gate G0–G10 with its
numbers and PASS/FAIL stated plainly, the preset-regression note
(design §4), and all recorded limitations.
Deliverables in out/v10/: REPORT-V10-GATES.md,
numbers-v10-sweep.json, numbers-v10-transport.json,
comp-windy-v37.png, spectrum plots as CSV or JSON (no images
needed beyond the still).
