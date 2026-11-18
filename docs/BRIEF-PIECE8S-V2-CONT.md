# CONTINUATION — 8 s piece base v2 (user review rework)

The previous session died on a transport error during the shadow-diagnosis
step and produced no v2 files. Read BRIEF-PIECE8S-V2.md in this folder for the
full spec (five user review points, new in-frame circuit choreography, gates,
deliverables). This note overrides two things:

1. The see_image / vision tool is UNAVAILABLE in your environment (it failed
   for every candidate). Do the shadow diagnosis NUMERICALLY instead: from
   v1's renders\frame_0072..frame_0084 and the moth's known animation state,
   measure the shadow blob (dark-region centroid offset from the moth centroid,
   core darkness ratio vs the perched moth's contact shadow, blob area vs
   height) with plain image statistics, and write the diagnosis into
   PIECE-V2-LOG.md from those numbers. Do not call see_image again.
2. Execution form: write build_piece_v2.py, launch the render DETACHED via
   WMI (Win32_Process.Create) with per-frame idempotent output to a NEW folder
   renders_v2\ (do NOT overwrite v1's renders\), verify from a NEW SSH session
   that the process is alive AND frames are advancing, report both checks,
   then END your session. Encoding/final verification happens in a follow-up
   session once frames complete. NEVER taskkill other jobs' processes.
