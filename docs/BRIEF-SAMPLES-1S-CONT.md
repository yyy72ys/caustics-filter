# CONTINUATION — 1 s composite samples (eta=1.0 vs eta=0.8)

The previous session died on a transport error AFTER measurement, before any
output files were written. Do not redo the measurement; the recorded numbers
from that session are: dome directions per-frame wall = 96: 306 s, 64: 238 s,
48: 203 s; direct C stats at the window (mean 1.0000002, p99 1.5482, max
2.1742 — matches the v8 route). Read BRIEF-SAMPLES-1S.md in this folder for
the full spec; this note overrides its execution form.

Execution form (mandatory):
1. Write ONE standalone runner script (Python) that computes, for base frames
   f0062..f0085: version A (eta=1, direct only, v8 route) and version B
   (eta=0.8, v9 route, **48 dome directions** — record this count in the
   report; diffuse master stays 2048^2), composites with the frozen composite,
   and writes per-frame PNGs idempotently (skip frames already done) plus a
   jobkit checkpoint.json update per frame, into offline\out\v9samples\.
   Reuse the same frame's direct transport for both versions where the code
   structure allows.
2. Launch the runner DETACHED via WMI (Win32_Process.Create), then verify from
   a NEW SSH session that the process is alive AND the output frame count is
   advancing. Report both checks in your final message, then END your session
   (the render continues without you).
3. Encoding to mp4 happens later (a follow-up session will encode + deliver
   once the frames are complete). Do NOT wait for the runner to finish.
NEVER taskkill python.exe. If WMI launch fails twice, fall back to foreground
execution and say so explicitly.
