# Freeze manifest — moth v1 / v3 series (re-assembled 2026-11-11)

This freeze was re-assembled on 2026-11-11 from the project records after the
original staging area was lost in a VM restart. The v3.1 commit (22593ba) is
the authentic original, recovered from the PC-side bundle. The two commits on
top (v3-series catch-up, moth v1) reproduce the originally planned freeze.

Excluded as regenerable: .npy caustic masters, .blend files (rebuildable from
the build scripts included here). Versions 3.2 and 3.3 were never tagged at
the time they were built; they are frozen by the catch-up commit, and the
v3.4 tag marks the top of the v3 series at freeze time. The v3.5 chunk (sky
split) and the 8 s piece base render landed after this freeze was specified
and will be frozen by a follow-up commit.

## v3 ALPHA follow-up freeze (assembled 2026-11-12; push approved by the owner 2026-11-12)

Per-version record of what was done: see `VERSION-HISTORY.md`.
Production follow-up (2026-11-12, after the ALPHA assembly): the three 8 s
pieces v18 (Blender moth base), v19 (photo base, preset 5 "deep clear pond")
and v20 (the v2 loop's place and weather on the current stack) — their
production code (offline/run_v18_prod.py, run_v19_prod.py, run_v20_prod.py)
and briefs (docs/BRIEF-V19-PROD8S-PHOTO.md, docs/BRIEF-V20-PROD8S-V2COND.md)
are included; the rendered videos themselves live in the project's delivery
records, not in this repo.

Contents: v3.5 (sky split, η=1 chosen for the piece), v3.6 water with the
2026-11-12 measurement correction (Fournier–Forand phase function, FF-CDF
kernel; HG g=0.85 refuted — see research/scattering-validation), v3.7 (kept as
rollback), v3.8 challenger, v3.9 low-wind floor as integrated by v3.10,
v3.10 (chunk 11 merge: elfouhaily-form spectrum is the default; region-wise
B-envelope gate replacing the retired universal B ≤ 0.012 clause — see
research/benvelope-cap-literature), the low-wind governance status
(docs/LOW-WIND-STATUS.md: UNRESOLVED — OPERATIONAL PRIOR RETAINED, decided
2026-11-12 after rally v42) with its code module (offline/lowwind_status.py)
and automated provenance test (offline/test_lowwind_provenance.py), and the
8 s piece base v2 (moth/piece-base-v2.mp4). Rain (3.2) remains a known-invalid
derivation, deferred to the v4 era by the owner's decision; v3 closes as
"v3 ALPHA".

Security note: this commit REMOVES offline/id_ed25519, a private SSH key that
was swept into the v3-series catch-up commit (1a09eed) when the offline tree
was copied wholesale. Removal from HEAD does not purge history. Rotation of
the key was considered and deferred by the owner's decision (2026-11-12);
the PC it authorizes is reachable only inside the owner's Tailscale network.

## v4 code freeze (assembled 2026-11-18; push approved by the owner 2026-11-18)

Contents: the v4-era production and measurement code of the offline tree,
synced from the production tree on the owner's PC (the same source as every
earlier freeze). New modules: the flowing-surface line (offline/stream_surface.py
through stream_surface4.py; drivers run_v4_1_1.py through run_v4_1_17.py), the
area-conserving mesh deposit (transport_mesh.py, transport_mesh_fast.py,
transport2.py), the bed line (bed_gen.py through bed_gen3.py; drivers
run_v4_2_1.py through run_v4_2_15.py, prep_v4_2_10.py), the deposit speedups
(det_fast.py, det_fast2.py, det_fast3.py), the environment preset table
(presets_env.py + selftest_preset_env.py, preset 5 registered), and the promo
material renderer (run_promo_v1.py). Updated: run_prod_v4.py (det_fast2 wired
in as the production deposit, bit-identical to the frozen det path, switchable
via DET_IMPL) and the post-freeze corrections to run_v15.py,
run_v15_report.py, run_v15_spectrum.py. docs/ROADMAP-V4.md and docs/PRESETS.md
are updated to their current versions; the per-version record is in
VERSION-HISTORY.md (v4 section).

Excluded: the v4 design drafts and rally records (the stage-1/stage-2 designs
are still in confirmation rallies and not yet adopted; they remain in the
project workspace), render outputs and .npy data (regenerable), and the PC
launch scripts (.cmd/.ps1/.bat) that earlier freezes also did not carry.

Security: before this push, the new tree and the full history were scanned for
secret material (key files, .env/credentials, private-key headers, token
patterns). No new secret material was found. The historical exposure of
offline/id_ed25519 in commit 1a09eed is unchanged and remains documented in
the v3 ALPHA section above.
