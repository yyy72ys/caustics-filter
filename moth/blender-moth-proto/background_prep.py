"""Produce background-bark-832x480.png (bark only, moth inpainted/quilting-filled).

Source photo: D:\\mutaro-work\\moth-on-bark.jpg (same file as butterfly-wan\\moth-on-bark.jpg).
Moth pixels are segmented by colour, then the hole is filled with min-error-cut
quilting from the photo's own bark/lichen patches (no external imagery, no
generative model).  Output is atomically committed (tmp -> verify -> os.replace)
and recorded in checkpoint.json.
"""
import os
import sys
import time

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from quilt import seg_moth, quilt, to_832x480  # noqa: E402
from checkpoint import Checkpoint  # noqa: E402

SRC = r"D:\mutaro-work\moth-on-bark.jpg"
OUT = os.path.join(HERE, "background-bark-832x480.png")
DBG = os.path.join(HERE, "_dbg")
TW, TH = 832, 480


def main():
    os.makedirs(DBG, exist_ok=True)
    t0 = time.time()
    img = cv2.imread(SRC, cv2.IMREAD_COLOR)
    if img is None:
        sys.exit("cannot read " + SRC)
    print("source", img.shape, "%.1f KB" % (os.path.getsize(SRC) / 1024.0))

    hole = seg_moth(img)
    print("moth hole %.1f%% of pixels" % (100.0 * hole.sum() / hole.size))

    bark = quilt(img, hole)
    crop = to_832x480(bark)

    cv2.imwrite(os.path.join(DBG, "bark_full_quilted.png"), bark)
    cv2.imwrite(os.path.join(DBG, "moth_mask.png"), (hole * 255).astype(np.uint8))

    # atomic commit: tmp -> verify -> os.replace
    tmp = OUT + ".tmp.png"
    cv2.imwrite(tmp, crop)
    chk = cv2.imread(tmp, cv2.IMREAD_UNCHANGED)
    assert chk is not None, "tmp unreadable"
    assert chk.shape[0] == TH and chk.shape[1] == TW, "bad size %r" % (chk.shape,)
    assert os.path.getsize(tmp) > 100000, "suspiciously small"
    os.replace(tmp, OUT)

    st = cv2.imread(OUT, cv2.IMREAD_COLOR)
    print("saved", OUT, st.shape, "%.1f KB" % (os.path.getsize(OUT) / 1024.0),
          "in %.1fs" % (time.time() - t0))

    cp = Checkpoint(os.path.join(HERE, "checkpoint.json"))
    cp.unit_done("background_ready")
    print("checkpoint:", cp.state())


if __name__ == "__main__":
    main()
