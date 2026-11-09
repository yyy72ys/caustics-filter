"""Fetch imageio-ffmpeg wheel from PyPI and unpack into _vendor (no stdout use)."""
import io
import json
import os
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
VENDOR = os.path.join(HERE, "_vendor")
os.makedirs(VENDOR, exist_ok=True)
LOG = os.path.join(HERE, "out", "vendor-install.log")


def log(msg):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def main():
    log("start vendor install")
    with urllib.request.urlopen("https://pypi.org/pypi/imageio-ffmpeg/json",
                                timeout=30) as r:
        meta = json.loads(r.read().decode("utf-8"))
    ver = meta["info"]["version"]
    log("latest version: " + ver)
    # pick win_amd64 wheel, prefer py3-none / cp310
    urls = [u for u in meta["urls"]
            if "win_amd64" in u["filename"] and u["filename"].endswith(".whl")]
    if not urls:
        log("ERROR: no win_amd64 wheel found")
        raise SystemExit(1)
    urls.sort(key=lambda u: (0 if "py3-none" in u["filename"] else 1,
                             u["filename"]))
    url = urls[0]["url"]
    log("wheel: " + urls[0]["filename"])
    with urllib.request.urlopen(url, timeout=120) as r:
        data = r.read()
    log("downloaded bytes: %d" % len(data))
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names = z.namelist()
        log("wheel members: %d" % len(names))
        for n in names:
            if n.endswith(".exe") or n.endswith("__init__.py") \
                    or n.endswith("version.py"):
                log("extract: " + n)
        z.extractall(VENDOR)
    log("OK vendor install")


if __name__ == "__main__":
    main()
