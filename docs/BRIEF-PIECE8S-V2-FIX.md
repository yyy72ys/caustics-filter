# BRIEF — piece base v2: left-edge margin fix (f80–f96) + re-encode

Folder: `D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\`
Build: `build_piece_v2.py`, scene `piece-v2.blend`, frames in `renders-v2\`.
Read `PIECE-V2-LOG.md` §3 first — it contains the full audit and the
prescribed fix. This brief executes exactly that prescription.

## Defect (measured twice, independently)
At f86–f87 the flying (left) moth's wing is CUT by the left frame edge:
silhouette margin 3 px, moth-coloured pixels (119,158,153) on column 0.
Margin stays < 20 px through ~f94. The brief gate was margin >= 10 px
everywhere (design target ~20 px). Everything else passes: seam
pixel-identical, shadow offset/penumbra/opacity all monotone in height,
opacity never 0.

## Task
1. Shift the flight path in X by **+0.061 world (+17 px)** applied over
   f80–f96 with a smooth falloff to 0 outside that range (e.g. a raised
   cosine over f78–f98). Implement it the way the shipped build stores
   the path (PATH_X keys / fcurve offsets) — do NOT touch anything
   outside the falloff range, and do NOT touch the staying (right) moth,
   the shadow curves, or frames 0/191 (the pixel-exact seam must
   survive: verify f0/f191 keys unchanged after your edit).
2. Re-render ONLY frames 86–94 into `renders-v2\` (overwrite those
   PNGs; single-frame renders, NOT the full animation — the shipped
   render block is not per-frame idempotent).
3. Re-audit with the folder's own tools (`tools\v2_audit.py` or your
   equivalent): report the moth silhouette margin for every frame
   f80–f96; gate = min margin >= 10 px and NO moth-coloured pixels on
   column 0. Also re-verify the seam from the PNGs (f191 vs f0 must
   still be pixel-identical, 0 differing px).
4. Re-encode `piece-base-v2.mp4` from `renders-v2\` (h264, 832x480,
   24 fps, 192 frames). Verify with ffprobe: 192 frames, 8.000 s.
5. Append the fix record (what you changed, margin table f80–f96,
   seam re-check, ffprobe result) to `PIECE-V2-LOG.md` as a new §5.
6. Copy the re-encoded `piece-base-v2.mp4` to the VM:
   scp to hatch VM `~/workspace/caustics-filter/delivery/piece-8s-base-v2/`
   is NOT available from your side — instead just leave the mp4 in the
   folder and finish; the coordinator fetches it. (State clearly in
   your final message that the mp4 is ready in the folder.)

## Rules
- Foreground execution with jobkit checkpoints per step; do not end
  the session until steps 1–5 are written to disk.
- Numerical verification only (no vision tool available — do not try).
- Do NOT kill blender/python processes you did not start; never kill
  pythonw (the LAN server).
- If the path representation makes the falloff shift impossible without
  touching the seam, STOP and report instead of improvising.
