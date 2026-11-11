"""Pipeline v4 phase 2: edge-pad frames to 736x540 and encode H.264 MP4.

Reads caustics-filter/offline/out/prod/ frame-%04d.png (735x539 native),
edge-replicates right+bottom 1px -> temp pad-%04d.png (736x540), encodes
with vendored static ffmpeg (libx264, yuv420p, crf 18, 24fps), deletes temp,
updates numbers-prod.json (video section) and appends to progress.log.

Usage (from repo root):
  python caustics-filter/offline/make_video_v4.py
All logging to files; stdout minimal.
"""
import datetime
import glob
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "prod")
PAD = os.path.join(OUT, "padtmp")
FFMPEG = os.path.join(HERE, "_vendor", "imageio_ffmpeg", "binaries",
                      "ffmpeg-win-x86_64-v9.0.2.exe")
MP4_NAME = "loop8s-736x540.mp4"
FPS = 24


def main():
    t_all = time.time()
    runlog = os.path.join(OUT, "prod-run.log")
    prolog = os.path.join(OUT, "progress.log")

    def rlog(msg):
        with open(runlog, "a", encoding="utf-8") as f:
            f.write(msg + "\n")

    def prog(msg):
        with open(prolog, "a", encoding="utf-8") as f:
            f.write(msg + "\n")

    frames = sorted(glob.glob(os.path.join(OUT, "frame-*.png")))
    if len(frames) != 192:
        rlog("ERROR video phase needs 192 frames, found %d" % len(frames))
        raise SystemExit("need 192 frames, found %d" % len(frames))
    rlog("video phase start: 192 frames")

    os.makedirs(PAD, exist_ok=True)
    for src in frames:
        base = os.path.basename(src)  # frame-XXXX.png
        img = np.asarray(Image.open(src).convert("RGB"))
        if img.shape != (539, 735, 3):
            rlog("ERROR unexpected frame shape %s %s"
                 % (src, img.shape))
            raise SystemExit("unexpected frame shape")
        pad = np.pad(img, ((0, 1), (0, 1), (0, 0)), mode="edge")
        Image.fromarray(pad).save(os.path.join(PAD, "pad-" + base[6:]))
    rlog("padding done 735x539 -> 736x540 edge-replicate")
    # spot check: padded corner pixels equal source edge pixels
    chk = np.asarray(Image.open(
        os.path.join(PAD, "pad-0000.png")).convert("RGB"))
    src0 = np.asarray(Image.open(frames[0]).convert("RGB"))
    assert chk.shape == (540, 736, 3)
    assert (chk[539, 735] == src0[538, 734]).all()
    assert (chk[539, 0] == src0[538, 0]).all()
    assert (chk[0, 735] == src0[0, 734]).all()
    rlog("pad edge check OK")

    mp4_path = os.path.join(OUT, MP4_NAME)
    ff_log = os.path.join(OUT, "ffmpeg-encode.log")
    cmd = [FFMPEG, "-y", "-framerate", str(FPS), "-i",
           os.path.join(PAD, "pad-%04d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-movflags", "+faststart", mp4_path]
    with open(ff_log, "wb") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    rlog("ffmpeg returncode=%d" % p.returncode)
    if p.returncode != 0 or not os.path.exists(mp4_path):
        rlog("ERROR ffmpeg failed")
        raise SystemExit("ffmpeg failed")
    mp4_bytes = os.path.getsize(mp4_path)
    rlog("mp4 bytes=%d" % mp4_bytes)
    shutil.rmtree(PAD, ignore_errors=True)

    npath = os.path.join(OUT, "numbers-prod.json")
    with open(npath, "r", encoding="utf-8") as f:
        numbers = json.load(f)
    numbers["video"] = {
        "file": MP4_NAME,
        "codec": "libx264 (H.264/AVC)",
        "pix_fmt": "yuv420p",
        "size": [736, 540],
        "fps": FPS,
        "frames": 192,
        "duration_s": 8.0,
        "crf": 18,
        "padding": ("photo-native 735x539 (odd dims) right+bottom 1px "
                    "edge-replicated to 736x540 for H.264"),
        "file_bytes": mp4_bytes,
        "encoder": "imageio-ffmpeg static ffmpeg v9.0.2 (vendored)",
    }
    numbers["timing_s"]["video_encode"] = time.time() - t_all
    numbers["timing_s"]["wall_total"] = (numbers["timing_s"]
                                         .get("wall_render_total", 0.0)
                                         + time.time() - t_all)
    with open(npath, "w", encoding="utf-8") as f:
        json.dump(numbers, f, indent=2, ensure_ascii=False)
    prog("video %s done bytes=%d wall=%.0fs"
         % (MP4_NAME, mp4_bytes, time.time() - t_all))
    rlog("video phase done")


if __name__ == "__main__":
    main()
