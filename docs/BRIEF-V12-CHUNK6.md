# BRIEF-V12-CHUNK6: version 3.6 — in-water transport, absorption/scattering split

Design source: PIPELINE-V8.md, chunk 6 (rally-confirmed gap P0-B, user-approved
2026-11-10; implementation ordered by the user 2026-11-11). Read it first.

## What changes (and only this)
The water path currently attenuates with a per-channel Beer-Lambert factor driven
by a turbidity scalar, plus a depth-scaled blur. Replace that model with:
1. **a(lambda) and b(lambda) separated.** Retire the turbidity scalar. Presets pick
   the a/b ratio by particle type (e.g. clear / algae-rich pond / sediment). Every
   coefficient value must carry a literature source in the report; anything that is
   a modelling choice is labelled as such. No invented constants.
2. **Single scattering in the transport.** Phase function Henyey-Greenstein,
   g from the preset (freshwater forward scattering g ~ 0.85 is the reference
   candidate). The caustic softening currently done by the empirical blur
   sigma ~ sqrt(b*z) must be REPLACED by a form derived from the scattering
   integral (the blur may remain as the numerical representation, but sigma must
   come out of the physics, and the derivation must be written in the report).
3. **CDOM port only**: an additive slot in a(lambda) for CDOM absorption, default
   zero, contents are v5's domain — build the port, not the model.

## Base and boundaries
- New code family run_v12* (v12_water.py + run_v12.py), importing v5–v11 code
  unchanged; never edit earlier versions' files. Spectrum source default stays
  `jonswap` (the v3.8 merge was withheld); the spectrum is NOT this chunk's topic.
- Frozen composite (A2.6 + Blow + E2.2) unchanged. Retina OFF. Rain untouched.
- Recipe gains a water section (a/b source: particle preset or explicit spectra);
  legacy recipes without it must reproduce the current behaviour (b = 0 path).

## Gates (all must be measured, none asserted)
- G1: b = 0 reproduces the pre-3.6 result on the preset-1 scene (bit-identical
  target; if not bit-identical, state the measured difference and its cause).
- G2 (the rally's falsification condition): at the SAME extinction c = a + b,
  different a/b ratios must separate colour and contrast measurably (report the
  numbers; deliver the pair of stills).
- G3: increasing b at fixed a lowers caustic peaks and widens skirts monotonically
  (C percentiles + a skirt-width measure, at least 3 levels of b).
- G4: two-way attenuation (down to the bottom + reflected return) applies on both
  legs — verify numerically against the a/b model on the z = 0.6 vs z = 4.0 pair.

## Stages (jobkit checkpoints in offline/out/v12/, one unit per stage)
1. v12_water.py: a/b model, particle presets with sources, HG single-scatter
   softening derivation + formula-level checks.
2. run_v12.py wiring + G1 regression render (preset-1 scene).
3. G2/G3/G4 experiments + stills (same-c pair, b sweep, depth pair).
4. REPORT-V12-GATES.md + numbers-v12.json. Merge rule: the a/b water path becomes
   the default iff G1–G4 pass; otherwise it stays behind the recipe switch and
   the report says why.

## Rules
- Foreground execution; do not end the session until the stage-4 deliverables are
  written. NEVER taskkill python.exe/blender.exe. Write ONLY in offline/out/v12/
  and the run_v12*/v12_* files. Numerical image statistics only (no vision).
  ASCII output, atomic writes, failures and deviations recorded openly.
