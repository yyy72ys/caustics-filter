"""Check vendored ffmpeg for libx264 (file IO only, no stdout)."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "_vendor", "imageio_ffmpeg", "binaries",
                   "ffmpeg-win-x86_64-v9.0.2.exe")
LOG = os.path.join(HERE, "out", "ffmpeg-enc.log")
RES = os.path.join(HERE, "out", "ffmpeg-enc-res.txt")

with open(LOG, "wb") as lf:
    p = subprocess.run([BIN, "-hide_banner", "-encoders"], stdout=lf,
                       stderr=subprocess.STDOUT)
data = open(LOG, "rb").read().decode("utf-8", "replace")
hits = [ln for ln in data.splitlines() if "264" in ln or "mpeg4" in ln]
with open(RES, "w", encoding="utf-8") as f:
    f.write("returncode=%d\n" % p.returncode)
    f.write("\n".join(hits[:20]) + "\n")
