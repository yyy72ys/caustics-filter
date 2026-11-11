"""Detached CONTINUATION driver for the 1-second eta samples.

Why this file exists
--------------------
`BRIEF-SAMPLES-1S-CONT.md` requires the frame render to keep running after the
session that started it ends (WMI detached), and forbids degrading/duplicating
the deliverable.  When this session arrived, another session had already built
and started `offline/run_v9_samples.py` (probe -> gate0 -> frames 62-65) in the
SAME output root `offline/out/v9samples/`.  Two renderers with two different
recipes would have produced two contradictory frame sets (Q=49/E=2.708 vs
Q=9/E=2.2) and two conflicting checkpoint ledgers.

So this driver does NOT re-implement the pipeline and does NOT use a second
recipe: it imports the owner's own stage functions and finishes the owner's
own idempotency units (`out/v9samples/parts/frame-XXXX.json`,
`frames/eta100|eta080/frame-0001..0024.png`).  Whatever it writes is
byte-identical to what the owner would write, because it is the owner's code
with the owner's constants.

What it does
------------
1. For every frame f in 62..85 whose part file is missing:
     * if a live `run_v9_samples.py` batch covers f, wait for it (bounded) --
       the owner is already computing that frame, so taking it would be pure
       duplicate CPU;
     * otherwise render it with the owner's `stage_frames(f, f)`.
2. When all 24 parts exist, wait for the owner to finish encode/stills/report;
   if the owner is gone and those artifacts are missing, run the owner's
   `stage_encode` / `stage_stills` / `stage_report`.
3. Log everything to `out/v9samples/continuation-runner.log` (a WMI launch has
   no console) and exit.

It never kills anything, never invents numbers, and every frame it adds is the
production route of `run_v9_samples.py` (Q=9 direct, 48-direction Perez dome,
2048^2, frozen v5 composite with E=2.2) -- the constants recorded in the
owner's own timing.json / gate0-t0-identity.json.

Usage (from caustics-filter/offline/):
  python -u _v9samples_continuation.py            # finish the window
  python -u _v9samples_continuation.py --dry-run  # report the plan only
"""
import argparse
import datetime
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_v9_samples as R          # the owner's pipeline (same constants)

LOG_PATH = os.path.join(R.OUT, "continuation-runner.log")
OWNER_SCRIPT = "run_v9_samples.py"
WAIT_OWNER_S = 25 * 60              # max wait for a live owner batch frame
POLL_S = 30
DEADLINE_S = 8 * 3600               # hard stop (frames + tail)


def open_log():
    os.makedirs(R.OUT, exist_ok=True)
    f = open(LOG_PATH, "a", encoding="utf-8", buffering=1)
    return f


def iso():
    return datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(
        timespec="seconds")


def owner_batches():
    """Command lines of live run_v9_samples.py --stage frames batches."""
    out = []
    try:
        import psutil
    except Exception:
        return out
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            cl = p.info.get("cmdline") or []
            joined = " ".join(cl)
            if OWNER_SCRIPT not in joined or "frames" not in joined:
                continue
            f0 = f1 = None
            for i, tok in enumerate(cl):
                if tok == "--from" and i + 1 < len(cl):
                    f0 = int(cl[i + 1])
                if tok == "--to" and i + 1 < len(cl):
                    f1 = int(cl[i + 1])
            out.append({"pid": p.info["pid"], "from": f0, "to": f1,
                        "cmd": joined})
        except Exception:
            continue
    return out


def owner_alive():
    try:
        import psutil
    except Exception:
        return False
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            joined = " ".join(p.info.get("cmdline") or [])
            if OWNER_SCRIPT in joined and p.info["pid"] != os.getpid():
                return True
        except Exception:
            continue
    return False


def part_exists(f):
    return os.path.exists(R.frame_part_path(f))


def frame_pngs(sub):
    d = os.path.join(R.OUT, "frames", sub)
    if not os.path.isdir(d):
        return 0
    return len([x for x in os.listdir(d) if x.endswith(".png")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--frames", default="%d-%d" % (R.FRAME_FIRST, R.FRAME_LAST))
    args = ap.parse_args()

    lo, hi = [int(x) for x in args.frames.split("-")]
    frames = list(range(lo, hi + 1))
    logf = open_log()

    def log(msg):
        line = "[%s] %s" % (iso(), msg)
        print(line, flush=True)
        logf.write(line + "\n")

    R.ensure_dirs()
    log("=" * 68)
    log("continuation driver start pid=%d  frames=%s" % (os.getpid(), args.frames))
    log("owner pipeline constants: Q_DIRECT=%d N_DIRS=%d E_COMP=%.2f grid=%d"
        % (R.Q_DIRECT, R.N_DIRS, R.E_COMP, R.GRID))
    log("existing: parts=%d eta100=%d eta080=%d"
        % (len([1 for f in frames if part_exists(f)]),
           frame_pngs("eta100"), frame_pngs("eta080")))
    for b in owner_batches():
        log("live owner batch: pid=%d frames %s..%s" % (b["pid"], b["from"],
                                                        b["to"]))
    if args.dry_run:
        missing = [f for f in frames if not part_exists(f)]
        log("dry-run: missing frames = %s" % missing)
        return 0

    t_start = time.time()
    while True:
        missing = [f for f in frames if not part_exists(f)]
        if not missing:
            log("all %d frames present" % len(frames))
            break
        f = missing[0]
        if time.time() - t_start > DEADLINE_S:
            log("DEADLINE reached with frames still missing: %s" % missing)
            return 2
        batch = [b for b in owner_batches()
                 if b["from"] is not None and b["to"] is not None
                 and b["from"] <= f <= b["to"]]
        if batch:
            b = batch[0]
            log("f%04d is inside live owner batch pid=%d (%d..%d): waiting"
                % (f, b["pid"], b["from"], b["to"]))
            t_wait = time.time()
            while (part_exists(f) is False and time.time() - t_wait < WAIT_OWNER_S
                   and owner_alive()):
                time.sleep(POLL_S)
            if part_exists(f):
                log("f%04d produced by the owner; continuing" % f)
                continue
            log("f%04d not produced after %.0f s / owner gone: taking it"
                % (f, time.time() - t_wait))
        t0 = time.time()
        try:
            R.stage_frames(f, f)
            log("f%04d rendered by continuation driver in %.1f s (part=%s)"
                % (f, time.time() - t0,
                   os.path.exists(R.frame_part_path(f))))
        except Exception as exc:
            log("f%04d FAILED: %s" % (f, exc))
            log(traceback.format_exc())
            time.sleep(60)
            continue
        if not part_exists(f):
            log("f%04d: stage returned but no part file -- stopping to avoid a "
                "silent gap" % f)
            return 3

    # ---------------- tail: encode / stills / report ----------------
    log("tail: eta100=%d eta080=%d parts=%d"
        % (frame_pngs("eta100"), frame_pngs("eta080"),
           len([1 for f in frames if part_exists(f)])))
    if frame_pngs("eta100") < len(frames) or frame_pngs("eta080") < len(frames):
        log("tail skipped: PNG count short of %d" % len(frames))
        return 4

    # give the owner time to finish its own tail before stepping in
    t_wait = time.time()
    while owner_alive() and time.time() - t_wait < 900:
        if os.path.exists(os.path.join(R.OUT, "SAMPLE-REPORT.md")):
            log("owner finished the report; nothing to do")
            return 0
        time.sleep(POLL_S)

    enc = os.path.join(R.OUT, "encode-report.json")
    stl = os.path.join(R.OUT, "stills-diff.json")
    rep = os.path.join(R.OUT, "SAMPLE-REPORT.md")
    try:
        if not os.path.exists(enc):
            log("tail: encode")
            R.stage_encode()
        else:
            log("tail: encode-report.json present")
        if not os.path.exists(stl):
            log("tail: stills")
            R.stage_stills()
        else:
            log("tail: stills-diff.json present")
        if not os.path.exists(rep):
            log("tail: report")
            R.stage_report()
        else:
            log("tail: SAMPLE-REPORT.md present")
    except Exception as exc:
        log("tail FAILED: %s" % exc)
        log(traceback.format_exc())
        return 5
    log("DONE wall=%.1f min" % ((time.time() - t_start) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
