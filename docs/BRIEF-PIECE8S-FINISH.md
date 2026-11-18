# BRIEF — piece8s COMPLETION session (PC)

Working folder (PC): `D:\yohey2026\dsh\muse\caustics-filter\blender-piece8s\`.
The first worker session died on an SSH broken pipe at ~01:44. Everything it did is
preserved on disk. Do NOT redo finished stages. Read `PIECE-LOG.md` in full first —
it records the stage state, a path adjustment made during the build, and the
f168 Z = +0.015 note (the brief's own knot value; already adjudicated in the log).

## Verified state at handover (ground truth, coordinator-checked)
- checkpoint.json units done: copy, build_edit, f0_parity, render_start.
- tools/seam_verify.py already ran on piece.blend: **SEAM_ANIM_PASS** — 49 fcurves,
  max |value(f0) − value(f191)| = 0.0 EXACT; motion-log f0 == f191 (16 fields);
  right root Δ = 0.0 over all 192 rows; left root f175..f191 == REST exactly.
- The DETACHED render (WMI, pid via render.log) was at ~180/192 frames and blender.exe
  was still alive. assemble_piece.py, tools/seam_render.py, tools/pxscan.py,
  tools/poll_render.py already exist in the folder.

## Your work — only the remaining stages, in order
1. **Render completion**: confirm renders\ holds all 192 frames frame_0000..frame_0191
   (jobkit verify_units with the png verifier pattern from assemble.py: header
   832×480, size > 100000). If blender is still running, poll until it exits; if it
   died before 192, relaunch the same detached render command (Blender skips
   nothing — a full re-render of missing frames only if the frame range can be
   limited; otherwise rerun the whole render detached and say so in the log).
   Record render wall time (from render.log / file timestamps). Checkpoint `render_done`.
2. **Pixel gates** (tools/pxscan.py, already written — if its method differs from the
   proto rounds' "fully out = 0 px" method, state both): left-moth pixel counts at
   f80 (>0), f84 (==0), f120 (==0), the first frame with px>0 during the return
   (must be in [139,143] — scan f132..f150), and f168 sanity (root at rest from
   animation data; pixel count > 50% of f0's count). If a gate fails because of the
   path, do NOT silently re-tune and re-render the whole piece: report the failing
   numbers in the log and in RESULT-PIECE.md as measured, and flag them in your
   final message — the coordinator adjudicates.
3. **Assembly**: run assemble_piece.py → piece-base.mp4 (ffprobe-verified before
   atomic commit: 832×480, 24 fps, nb_read_frames == 192, duration 8.0 s, h264),
   encode-report.json, checkpoint `assemble`.
4. **Render-diff seam check**: run tools/seam_render.py (frame_0000 vs frame_0191:
   overall mean |Δ| on both 0–255 and 0–1 scales; moth mask from
   |frame − background-bark-832x480.png| > 12 dilated 3 px; mean |Δ| inside/outside
   mask; % pixels |Δ| > 8 in/out; GATE outside-mask mean ≤ 0.002 on the 0–1 scale).
   Append to seam-verify.json. Checkpoint `seam_verify`.
5. **Stills**: copy renders\frame_0000.png→f0.png, 0072→f72.png, 0084→f84.png,
   0120→f120.png, 0168→f168.png, 0191→f191.png (piece folder root).
6. **RESULT-PIECE.md** exactly as specified in BRIEF-PIECE8S-WORKER.md Step 6
   (timeline actuals, pixel-gate table, seam verification numbers, per-segment
   provenance labels, ffprobe block, render wall time, mp4 bytes, honest remains —
   including the hair/TAA shimmer note and the f168 Z note from PIECE-LOG).
   Checkpoint `report`.
7. Append a completion section to PIECE-LOG.md (what this session did, any incident).

Final message: the summary with key numbers (ffprobe line, seam diffs, pixel gates,
wall time). Do NOT copy deliverables to the VM — the coordinator fetches them.
Model: deepseek-v4.1-flash only (already fixed by the caller).
