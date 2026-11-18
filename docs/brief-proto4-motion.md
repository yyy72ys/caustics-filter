# BRIEF proto4-A: real moth motion research (YouTube frame-step analysis)

You are running headless on the user's PC. Execute this brief exactly, to completion, without stopping to ask questions. If a step fails, log the failure and continue with the next step/clip.

## Purpose
A Blender model of the emerald geometrid moth *Geometra papilionaria* needs MEASURED real-motion data to design its animation: perched micro-motion, the pre-takeoff sequence, and takeoff. Every number you report must trace to an actual measurement on downloaded footage (frame-difference series or a stated image measurement). NEVER invent numbers or fill gaps from general knowledge. If something cannot be measured from the footage, write exactly that.

## Working directory and checkpointing
- Work in `D:\yohey2026\dsh\muse\caustics-filter\moth-motion\` (create it). All outputs go there.
- Use the jobkit checkpoint module `D:\mutaro-work\jobkit\checkpoint.py` with `checkpoint.json` in the moth-motion dir. Units: setup, search, dl-<clipid>, analysis-<clipid>, report. Call unit_done after each unit so a dead session can resume without redoing finished units. On start, check the checkpoint and skip units already done.

## Steps
1. SETUP: Check yt-dlp (`yt-dlp --version`). If missing, create a venv at `moth-motion\venv` and `pip install yt-dlp` into it; use that yt-dlp afterwards. Check `ffmpeg -version` and `ffprobe -version` (they are installed). Check python has numpy and PIL (system python 3.10 does).
2. SEARCH: Use `yt-dlp --flat-playlist --print "%(id)s | %(title)s | %(duration)s | %(channel)s | %(webpage_url)s" "ytsearch8:<query>"` for each query:
   - `Geometra papilionaria`
   - `large emerald moth`
   - `emerald moth takeoff`
   - `moth takeoff 6000 fps`
   - `Ant Lab moths slow motion`
   - `geometrid moth slow motion`
   - `moth perched antennae movement`
   - `シャクガ 飛翔`
   Save ALL results (deduped by id) to `candidates.txt`.
3. SELECTION: Choose 2-4 clips to analyse. Priority order:
   (a) species is *Geometra papilionaria* or another emerald Geometrinae (e.g. a Nemoria emerald) — title/description evidence only;
   (b) a perched moth at rest with visible micro-movement;
   (c) a takeoff, slow-motion strongly preferred (the Ant Lab / Dr Adrian Smith high-fps moth takeoff videos are a known good lead — one includes a geometrid and an emerald);
   (d) moth large in frame, camera steady.
   Aim for at least one perched-behaviour clip and one takeoff clip. Write `selection.md`: chosen clips with URL, title, channel, duration, and the species-confidence + reason for each; also list near-miss candidates you rejected and why.
4. DOWNLOAD: For each chosen clip: `yt-dlp -f "bv*[height<=1080]+ba/b[height<=1080]" --merge-output-format mp4 -o "clips\%(id)s.%(ext)s" <url>`. If a download fails (403, bot check, sign-in required), append the clip id + exact error to `blockers.log` and try the next-best candidate from candidates.txt. Fallback if YouTube is entirely blocked: search Wikimedia Commons (commons.wikimedia.org API) and archive.org for *Geometra papilionaria* / emerald moth video files and download any real footage found; record what happened either way. Do NOT use browser cookies, account logins, or any Chrome profile. If no footage can be obtained at all, write that plainly in the report and finish the other steps that are still possible.
5. ANALYSIS per downloaded clip (ffmpeg + python numpy/PIL):
   - `ffprobe` the file: fps, duration, resolution, frame count -> save `meta-<id>.json`.
   - Extract frames as JPG into `frames-<id>\` (every frame for clips under ~3000 frames; every 2nd frame for longer ones — record the stride).
   - Frame-difference series: grayscale mean absolute difference between consecutive extracted frames. If the moth occupies a clear sub-region, ALSO compute the series restricted to a moth bounding box that you estimate and record. Save `diff-<id>.csv` (frame_index, time_s, diff_full, diff_mothbox if applicable).
   - Event detection: threshold = median + 5*MAD of the diff series (state the actual threshold value). Classify: short spikes (above threshold for 1-10 extracted frames, then back below) = candidate flicks/micro-movements; sustained elevation lasting until the moth leaves the frame or the clip ends = takeoff.
   - Perched behaviour: for each detected flick list start time (s), duration (frames and ms — convert using the clip's real fps and your extraction stride; for slow-motion clips ALSO give the real-time duration using the playback slowdown factor IF the title/description states the capture fps, otherwise state that real-time conversion is not possible), peak diff value (amplitude proxy), and the interval since the previous flick. List ALL intervals in order, then give mean/std/CV and state from those numbers whether the timing reads periodic or irregular.
   - Takeoff (from frame strips, step 6): describe the sequence in order with frame numbers: what moves first (antennae / wing raise / body shift / leg push), any hesitation (gap in frames between first detectable movement and liftoff), wing angle at the last perched frame and first airborne frame (measure from image geometry if the view allows; otherwise write "not measurable from this view"), left/right wing symmetry in the first wingbeat (qualitative, with frame refs), body pitch change (qualitative or measured, labelled which), first-wingbeat duration in frames.
   - Antennae: any visible antenna movement events with timestamps — qualitative from strips, labelled qualitative.
   - Two-moth reaction: if any analysed clip shows one moth flying off near another perched moth, record whether/when the perched one reacts (flick within 1 s after, etc.). If no such clip was analysed, state that.
6. STRIPS: For each takeoff clip, build montage strip images with ffmpeg tile or PIL: (a) the last ~1 s perched before takeoff (every 2nd or 3rd frame), (b) the takeoff moment frame-by-frame, (c) the first ~0.5 s of flight. Save as `strip-<id>-pre.png`, `strip-<id>-takeoff.png`, `strip-<id>-flight.png`. For perched clips, build one strip `strip-<id>-flicks.png` showing 2-3 frames around each of the 3 largest detected flicks. Keep each strip under ~4000 px wide.
7. REPORT: Write `analysis-report.md` in moth-motion with: chosen clips (URL/title/channel/duration/fps/species confidence), per-clip event tables, flick interval lists + CV + periodicity verdict, the pre-takeoff and takeoff narratives with frame refs, takeoff measurements, strip file names, all blockers/errors encountered, and a final section titled "NOT established by this footage" listing everything that could not be measured (be specific). Also write `events-summary.json` with the key numbers machine-readable (per clip: fps, slowdown factor if known, flicks [{t_start_s, dur_frames, peak_diff}], intervals_s, takeoff {first_move_frame, liftoff_frame, hesitation_frames, first_wingbeat_frames}).

## Rules
- ASCII file names only. No fabricated numbers, ever.
- Work only inside moth-motion (plus the rally-in brief you are reading). Do not touch other folders.
- Run to completion. Log blockers and continue; never wait for input.
