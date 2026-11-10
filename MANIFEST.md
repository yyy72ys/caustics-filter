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
