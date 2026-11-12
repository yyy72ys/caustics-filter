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
