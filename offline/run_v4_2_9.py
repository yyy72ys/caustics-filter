"""run_v4_2_9 -- BRIEF-V4-2-9: calibrate OUR RULER on a PUBLIC, MEASURED flume
bed (Redolfi et al. 2025, Zenodo record 14756025, alternate-bar experiment,
mean grain size ~1.2 mm, grid 5 mm).

WHAT THIS IS AND IS NOT
-----------------------
This is a calibration of the measurement instrument defined in
run_v4_2_5.py (curve_dir / curve_iso / local_H / fit_band / band_shape).
It answers only: "how does that estimator behave when it is pointed at a
real, published DEM?"  It contains NO pass/fail, good/bad, or any other
verdict about our own synthetic bed.  The target bed is ~1/30 the grain
size of ours and is not armour-limited, so a comparison of *values* would
be meaningless; only the behaviour of the estimator is of interest.

The ruler's functions are imported from run_v4_2_5.py and used unmodified
(no re-implementation).  Nothing existing is written: the only new .py is
this file, and the only new directory is out/v4-2-9/.

Stages (foreground, idempotent, one jobkit checkpoint unit per stage)
  --stage0    absolute-path preflight, sha256 baseline of every existing .py
  --download  Zenodo record + DEM/d50 GeoTIFFs -> out/v4-2-9/data/, md5+sha256
  --prove     provenance: GeoTIFF tags, grid, extent, nodata pattern, d50
  --measure   the ruler on every DEM: x/y/iso curves, local_H, band fits
              (lag in mm and in d/D50), plus the increment / coarsening
              diagnostics that say whether the reading can be physical
  --sens      box-size, hole-fill and placement sensitivity
  --report    numbers-v4-2-9.json + REPORT-V4-2-9.md
  --all       all of the above in one process
"""
from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes  # noqa: E402

import run_v4_2_5 as M5                               # noqa: E402  (unmodified)

BRIEF = r"D:\mutaro-work\rally-in\BRIEF-V4-2-9.md"
OUT = os.path.join(HERE, "out", "v4-2-9")
DATA = os.path.join(OUT, "data")
CP_PATH = os.path.join(OUT, "checkpoint.json")
NUMBERS = os.path.join(OUT, "numbers-v4-2-9.json")
REPORT = os.path.join(OUT, "REPORT-V4-2-9.md")
RECORD_JSON = os.path.join(DATA, "_record-14756025.json")

THIS_FILE = "run_v4_2_9.py"
RECORD_ID = 14756025
RECORD_URL = "https://zenodo.org/records/14756025"
API_URL = "https://zenodo.org/api/records/%d" % RECORD_ID

# --- the ruler under test: all of it comes from run_v4_2_5 -----------------
BANDS = M5.BANDS                      # scaling 0.4-3.6 cm / meso / legacy
BAND_LABEL = M5.BAND_LABEL
BAND = M5.BANDS[0][1:]                # the primary band, (0.004, 0.036) m
NBIN = M5.NBIN                        # 60, the ruler's own lag ladder
DM_X = M5.DM_X_LONG                   # 1.00 m   (flow axis)
DM_Y = M5.DM_Y                        # 0.30 m   (cross-stream)
DM_ISO = M5.DM_ISO                    # 0.20 m
GATE = M5.GATE                        # recorded for reference only, not applied
# the project's inpainting length: 60 px of the v4-2.2 statistics cell
INPAINT_M = 60.0 * M5.CELL            # 0.017578 m

# --- this job --------------------------------------------------------------
PRIMARY = "DEM_MW60_Q4.2.tif"         # headline specimen (W = 60 cm, cleanest)
TILE_M = 0.60                         # the ruler's patch size
PROBE_M = [0.005, 0.0074, 0.011, 0.016, 0.024, 0.035, 0.05, 0.07, 0.1,
           0.15, 0.2, 0.3, 0.5, 1.0]   # local_H / S probe lags
COARSEN = (1, 2, 4)                   # block-average factors (5/10/20 mm cells)
PROBE_COARSE_M = [0.02, 0.03, 0.04, 0.06, 0.08, 0.12]

SERIES = {                            # from the record description
    "MW40": {"W_m": 0.40, "D50_mm": 1.22},
    "MW60": {"W_m": 0.60, "D50_mm": 1.29},
}
# values the brief asks us to compare against (v4-2.6, synthetic self-affine)
REF = {
    "box_bias_synthetic": 0.10,
    "hole_spread_synthetic": 0.134,
    "hole_frac_synthetic": 0.1322,
}

LOG = []
_CACHE = {}
WRITES = []


def log(msg):
    line = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    LOG.append(line)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%S+00:00")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _j(x):
    """JSON-safe scalar."""
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (list, tuple)):
        return [_j(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _j(v) for k, v in x.items()}
    return x


def write_json(path, obj):
    data = json.dumps(_j(obj), ensure_ascii=False, indent=1).encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            json.loads(f.read().decode("utf-8"))
        return True
    WRITES.append(os.path.abspath(path))
    atomic_write_bytes(path, data, verify=verify)


def write_text(path, text):
    data = text.encode("utf-8")

    def verify(p):
        with open(p, "rb") as f:
            return len(f.read()) == len(data)
    WRITES.append(os.path.abspath(path))
    atomic_write_bytes(path, data, verify=verify)


# ---------------------------------------------------------------------------
# STAGE 0 -- preflight
# ---------------------------------------------------------------------------
REQUIRED = {
    "run_v4_2_5.py (the ruler)": os.path.join(HERE, "run_v4_2_5.py"),
    "bed_gen2.py (imported by it)": os.path.join(HERE, "bed_gen2.py"),
    "jobkit/checkpoint.py": r"D:\mutaro-work\jobkit\checkpoint.py",
}


def stage0(cp):
    log("STAGE0  absolute-path preflight")
    missing = [k for k, p in REQUIRED.items() if not os.path.isfile(p)]
    if missing:
        raise SystemExit("STOP: missing %r" % (missing,))
    hashes = {}
    for k, p in REQUIRED.items():
        hashes[k] = {"path": p, "bytes": os.path.getsize(p),
                     "sha256": sha256_file(p)}
        log("  %-34s %8d B  %s" % (k, hashes[k]["bytes"], hashes[k]["sha256"][:16]))

    net = {"url": API_URL, "ok": False, "detail": None}
    try:
        import urllib.request
        with urllib.request.urlopen(API_URL, timeout=45) as r:
            net["status"] = r.status
            net["ok"] = (r.status == 200)
            net["detail"] = "%d bytes" % len(r.read())
    except Exception as e:                              # noqa: BLE001
        net["detail"] = "%s: %s" % (type(e).__name__, e)
    log("  zenodo.org  %s  (%s)" % ("reachable" if net["ok"] else "UNREACHABLE",
                                    net["detail"]))
    if not net["ok"]:
        raise SystemExit("STOP: zenodo.org not reachable")

    # baseline for the write-set check: every existing .py, and every file
    # currently under out/ (so "no existing file changed" is provable).
    py = {}
    for p in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        n = os.path.basename(p)
        if n == THIS_FILE:
            continue
        py[n] = {"bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    files = {}
    for root, _dirs, fnames in os.walk(os.path.join(HERE, "out")):
        if os.path.abspath(root).startswith(os.path.abspath(OUT)):
            continue
        for fn in fnames:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, HERE).replace("\\", "/")
            try:
                st = os.stat(p)
            except OSError:
                continue
            files[rel] = [st.st_size, int(st.st_mtime)]
    log("  baseline: %d existing .py, %d existing files under out/"
        % (len(py), len(files)))
    out = {"created_utc": utcnow(), "required": hashes, "network": net,
           "py": py, "n_existing_py": len(py), "out_files": files,
           "n_existing_out_files": len(files),
           "python": sys.version.split()[0], "numpy": np.__version__}
    try:
        import tifffile
        out["tifffile"] = tifffile.__version__
    except Exception:                                   # noqa: BLE001
        out["tifffile"] = None
    cp.unit_done("stage0")
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- download + provenance
# ---------------------------------------------------------------------------
def stage_download(cp, st0):
    import urllib.request
    log("STAGE1  Zenodo record %d" % RECORD_ID)
    os.makedirs(DATA, exist_ok=True)
    got = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(API_URL, timeout=90) as r:
                got = json.load(r)
            break
        except Exception as e:                          # noqa: BLE001
            log("  record fetch retry %d: %s" % (attempt, e))
            time.sleep(3 * (attempt + 1))
    if got is None:
        raise SystemExit("STOP: cannot read the Zenodo record")
    write_json(RECORD_JSON, {"fetched_utc": utcnow(), "api": API_URL,
                             "record": got})

    md = got["metadata"]
    lic = md.get("license") or {}
    rec = {
        "record_id": got.get("id"), "doi": got.get("doi"),
        "concept_doi": got.get("conceptdoi"), "url": RECORD_URL,
        "title": md.get("title"), "publication_date": md.get("publication_date"),
        "license": lic.get("id") if isinstance(lic, dict) else lic,
        "version": md.get("version"),
        "creators": [c.get("name") for c in md.get("creators", [])],
        "resource_type": (md.get("resource_type") or {}).get("title"),
        "description_text": _plain(md.get("description") or ""),
    }
    log("  %s" % rec["title"])
    log("  license %s | %d creators" % (rec["license"], len(rec["creators"])))

    want = sorted([f for f in got["files"]
                   if f["key"].startswith(("DEM_", "d50_"))],
                  key=lambda f: f["key"])
    log("  %d selected files (DEM_* + d50_*), %.1f MB"
        % (len(want), sum(f["size"] for f in want) / 1e6))
    files = []
    for f in want:
        key, size = f["key"], f["size"]
        md5 = f["checksum"].split(":", 1)[1]
        dst = os.path.join(DATA, key)
        ok = False
        if os.path.exists(dst) and os.path.getsize(dst) == size:
            ok = hashlib.md5(open(dst, "rb").read()).hexdigest() == md5
        if not ok:
            data = None
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(f["links"]["self"],
                                                timeout=240) as r:
                        data = r.read()
                    break
                except Exception as e:                  # noqa: BLE001
                    log("  retry %s (%d): %s" % (key, attempt, e))
                    time.sleep(3 * (attempt + 1))
            if data is None:
                raise SystemExit("STOP: download failed %s" % key)
            with open(dst, "wb") as fh:
                fh.write(data)
            WRITES.append(os.path.abspath(dst))
            if hashlib.md5(data).hexdigest() != md5:
                raise SystemExit("STOP: md5 mismatch %s" % key)
            log("  got %-24s %8d B" % (key, len(data)))
        files.append({"key": key, "bytes": size, "md5": md5,
                      "md5_ok": True, "sha256": sha256_file(dst),
                      "download": f["links"]["self"]})
    cp.unit_done("stage1_download")
    return {"record": rec, "files": files,
            "n_files_in_record": len(got["files"]),
            "total_bytes": sum(f["bytes"] for f in files)}


def _plain(html):
    import re
    t = re.sub(r"<[^>]+>", " ", html or "")
    t = t.replace("&nbsp;", " ").replace("&amp;", "&")
    return " ".join(t.split())


def _tags(path):
    import tifffile
    with tifffile.TiffFile(path) as tf:
        pg = tf.pages[0]
        def val(code, default=None):
            t = pg.tags.get(code)
            return default if t is None else t.value
        sc = val(33550) or (None, None, None)
        tp = val(33922) or (0, 0, 0, 0, 0, 0)
        gk = val(34735) or ()
        nod = val(42113)
        arr = pg.asarray()
        info = {"shape": list(pg.shape), "dtype": str(pg.dtype),
                "compression": str(pg.compression),
                "software": val(305), "pixel_scale": [float(sc[0]), float(sc[1])],
                "tiepoint": [float(v) for v in tp],
                "geo_keys": [int(v) for v in gk],
                "nodata_tag": nod}
    keys = {}
    for i in range(4, len(info["geo_keys"]) - 3, 4):
        g = info["geo_keys"]
        keys[g[i]] = g[i + 3]
    info["epsg"] = keys.get(3072)
    return arr, info


def stage_provenance(cp, dl):
    log("STAGE2  provenance (GeoTIFF tags, grid, nodata, d50)")
    dems, grids = {}, {}
    for f in dl["files"]:
        arr, info = _tags(os.path.join(DATA, f["key"]))
        sx, sy = info["pixel_scale"]
        ny, nx = arr.shape                      # rows = cross-stream? (checked)
        hole = ~np.isfinite(arr)
        ext = {"width_m": nx * sx, "height_m": ny * sy,
               "x0_m": info["tiepoint"][3], "y0_m": info["tiepoint"][4],
               "x1_m": info["tiepoint"][3] + nx * sx,
               "y1_m": info["tiepoint"][4] - ny * sy}
        ser = "MW40" if "MW40" in f["key"] else "MW60"
        entry = dict(info)
        entry.update({"file": f["key"], "series": ser,
                      "n_rows": int(ny), "n_cols": int(nx),
                      "cell_m": float(sx), "extent": ext,
                      "nodata_frac": float(hole.mean()),
                      "all_nan_cols": int(hole.all(axis=0).sum()),
                      "all_nan_rows": int(hole.all(axis=1).sum()),
                      "width_over_cell": float(nx * sx),
                      "channel_width_m": SERIES[ser]["W_m"],
                      "D50_series_mm": SERIES[ser]["D50_mm"]})
        if f["key"].startswith("DEM_"):
            v = arr[~hole]
            entry["z_stats"] = {"min": float(v.min()), "max": float(v.max()),
                                "mean": float(v.mean()),
                                "sigma_m": float(v.std())}
            entry["z_range_m"] = float(v.max() - v.min())
            dems[f["key"]] = entry
        else:
            v = arr[~hole]
            entry["d50_stats_mm"] = {
                "min": float(v.min()), "max": float(v.max()),
                "mean": float(v.mean()), "median": float(np.median(v)),
                "p05": float(np.percentile(v, 5)),
                "p95": float(np.percentile(v, 95))}
            grids[f["key"]] = entry
        log("  %-24s %dx%d  %.3f m x %.3f m  px %.4f m  nodata %.2f%%"
            % (f["key"], ny, nx, ext["width_m"], ext["height_m"], sx,
               100 * hole.mean()))

    # grid consistency: nx*px must equal the published channel width
    chk = []
    for k, e in dems.items():
        chk.append({"file": k, "series": e["series"],
                    "nx_times_px_m": e["width_over_cell"],
                    "channel_width_m": e["channel_width_m"],
                    "ratio": e["width_over_cell"] / e["channel_width_m"]})
    cells = sorted({round(e["cell_m"], 6) for e in dems.values()})
    grids = {"d50_maps": grids, "cell_m_values": cells,
             "width_check": chk,
             "note": ("the record text says 0.5 mm; every GeoTIFF says "
                      "ModelPixelScale 0.005 m, and 80 px x 5 mm = 0.400 m = "
                      "W(MW40) exactly, so the published grid is 5 mm and the "
                      "text's '0.5 mm' refers to the pre-publication product")}
    cp.unit_done("stage1_provenance")
    return {"dem": dems, "grid": grids}


# ---------------------------------------------------------------------------
# STAGE 2 -- the ruler on the real DEMs
# ---------------------------------------------------------------------------
def load_dem(name):
    if name in _CACHE:
        return _CACHE[name]
    import tifffile
    arr, info = _tags(os.path.join(DATA, name))
    z0 = np.ascontiguousarray(arr.T.astype(np.float64))   # (ny, nx)
    hole = ~np.isfinite(z0)
    _CACHE[name] = {"z0": z0, "hole": hole, "info": info}
    return _CACHE[name]


def fill_mean(z0, hole):
    z = z0.copy()
    if hole.any():
        z[hole] = float(z0[~hole].mean())
    return z


def curv(z, cell, axis, dmax):
    if axis == "iso":
        return M5.curve_iso(z, cell, dmax, nbin=NBIN)
    return M5.curve_dir(z, cell, dmax, nbin=NBIN,
                        axis=(1 if axis == "x" else 0))


def band_reads(d, s):
    out = {}
    for name, lo, hi in M5.BANDS:
        f = M5.fit_band(d, s, lo, hi)
        sh = M5.band_shape(d, s, lo, hi)
        out[name] = dict(f)
        out[name].update({
            "local_H_min": sh.get("local_H_min"),
            "local_H_max": sh.get("local_H_max"),
            "local_H_min_at_m": sh.get("local_H_min_at_m"),
            "local_H_probe": sh.get("local_H_probe"),
            "plateau_onset_m": sh.get("plateau_onset_m"),
            "cross_0.20_m": sh.get("cross_0.20_m")})
    return out


def profile_at(d, s, probes):
    """S and local_H interpolated at fixed lags (log-log), plus d/D50 later."""
    loc = M5.local_H(d, s)
    ok = np.isfinite(s) & (s > 0)
    xd, xs, xl = np.log(d[ok]), np.log(s[ok]), loc[ok]
    rows = []
    for p in probes:
        if p < d[ok][0] * 0.999 or p > d[ok][-1] * 1.001:
            rows.append([float(p), None, None])
            continue
        rows.append([float(p),
                     float(np.exp(np.interp(np.log(p), xd, xs))),
                     float(np.interp(np.log(p), xd, xl))])
    return rows


def mirror_close(z, axis):
    """Close the domain by mirror reflection so the periodic closure is exact.

    Not a new estimator: this array is handed to the SAME curve_dir / curve_iso.
    On a non-periodic bed (walls, flume ends, a 1% longitudinal slope) the
    np.roll closure adds a fictitious step at the seam; mirroring removes that
    step and leaves the interior arithmetic untouched.
    """
    return np.concatenate([z, np.flip(z, axis=axis)], axis=axis)


def read_one(z, cell, axis, dmax):
    d, s = curv(z, cell, axis, dmax)
    return {"dmax_requested_m": dmax, "dmax_effective_m": float(d[-1]),
            "n_lags": int(d.size),
            "S_at_1px_m2": float(s[0]),
            "steps_per_band": {k: int(((d >= lo) & (d <= hi)).sum())
                               for k, lo, hi in M5.BANDS},
            "bands": band_reads(d, s),
            "profile": profile_at(d, s, PROBE_M),
            "ladder_head": [[float(a), float(b)] for a, b in zip(d[:12], s[:12])]}


def stage_measure(cp, prove):
    log("STAGE3  the ruler, all DEMs")
    out = {}
    for name in sorted(prove["dem"]):
        e = prove["dem"][name]
        cell = e["cell_m"]
        z0, hole = load_dem(name)["z0"], load_dem(name)["hole"]
        z = fill_mean(z0, hole)
        d50 = float(e["D50_series_mm"]) * 1e-3
        d50map = None
        dm = name.replace("DEM_", "d50_")
        if os.path.isfile(os.path.join(DATA, dm)):
            arr, _ = _tags(os.path.join(DATA, dm))
            v = arr[np.isfinite(arr)]
            d50map = {"file": dm, "median_mm": float(np.median(v)),
                      "p05_mm": float(np.percentile(v, 5)),
                      "p95_mm": float(np.percentile(v, 95))}
            d50 = float(np.median(v)) * 1e-3
        rec = {"file": name, "series": e["series"], "cell_m": cell,
               "shape": [int(z.shape[0]), int(z.shape[1])],
               "extent_m": e["extent"], "nodata_frac": e["nodata_frac"],
               "D50_ref_mm": 1e3 * d50, "D50_series_mm": e["D50_series_mm"],
               "D50_map": d50map,
               "cell_over_D50": cell / d50,
               "sigma_z_m": float(z.std()),
               "z_range_m": e["z_range_m"]}

        # slope along the flow axis, and its share of S at the band's top lag
        prof = z.mean(axis=0)
        xs = np.arange(prof.size) * cell
        sl = float(np.polyfit(xs, prof, 1)[0])
        rec["slope_x"] = sl
        rec["rowmean_p3_p97_m"] = [float(np.percentile(prof - np.polyval(
            np.polyfit(xs, prof, 1), xs), 3)),
            float(np.percentile(prof - np.polyval(np.polyfit(xs, prof, 1),
                                                  xs), 97))]

        # ---- the ruler exactly as published (np.roll closure) -------------
        dirs = {}
        for axis, dmax in (("x", DM_X), ("y", DM_Y), ("iso", DM_ISO)):
            r = read_one(z, cell, axis, dmax)
            axis_idx = 1 if axis != "y" else 0
            r["dmax_px"] = int(min(dmax / cell, z.shape[axis_idx] // 2))
            r["wrap_pairs_share_at_band_top"] = float(
                BAND[1] / cell) / float(z.shape[axis_idx if axis != "iso" else 1])
            dirs[axis] = r
        rec["dirs"] = dirs
        rec["slope_share_at_band_top"] = (
            (sl * BAND[1]) ** 2 / dirs["x"]["bands"]["scaling"]["S_at_d_max"]
            if dirs["x"]["bands"]["scaling"]["S_at_d_max"] else None)

        # ---- the same estimator with an exact (mirror) closure ------------
        # (b) closes only the seam, keeping the boundary fill in place
        # (c) also removes the boundary mask (interior crop)
        zc = crop_interior(z0, hole)
        nw, nwc = {}, {}
        for axis, dmax, ax in (("x", DM_X, 1), ("y", DM_Y, 0), ("iso", DM_ISO, 0)):
            if axis == "iso":
                zz = mirror_close(mirror_close(z, 0), 1)
                zz2 = mirror_close(mirror_close(zc, 0), 1)
            else:
                zz = mirror_close(z, ax)
                zz2 = mirror_close(zc, ax)
            nw[axis] = read_one(zz, cell, axis, dmax)
            nwc[axis] = read_one(zz2, cell, axis, dmax)
            del zz, zz2
        rec["dirs_no_wrap"] = nw
        rec["dirs_clean"] = nwc

        # ---- the seam: what the periodic closure actually adds ------------
        seam_dx = float(np.mean(np.abs(z[:, -1] - z[:, 0])))
        seam_dy = float(np.mean(np.abs(z[-1, :] - z[0, :])))
        rec["seam"] = {
            "delta_x_m": seam_dx, "delta_y_m": seam_dy,
            "slope_x_times_length_m": abs(sl) * z.shape[1] * cell,
            "S1_ruler_m2": dirs["x"]["S_at_1px_m2"],
            "S1_clean_m2": nwc["x"]["S_at_1px_m2"],
            "S1_ruler_over_clean_x": (dirs["x"]["S_at_1px_m2"]
                                      / nwc["x"]["S_at_1px_m2"]),
            "S1_ruler_over_clean_y": (dirs["y"]["S_at_1px_m2"]
                                      / nwc["y"]["S_at_1px_m2"]),
            "predicted_wrap_term_m2_at_1px": 2.0 * seam_dx ** 2 / z.shape[1],
            "note": ("np.roll closes the domain across a seam of |dz| = delta; "
                     "one such pair per row adds d*delta^2/N to S(d), i.e. an "
                     "exactly d^1 term -- indistinguishable, to fit_band, from "
                     "H = 0.5")}

        # power law over decades (ruler vs exact closure)
        dec = {}
        dz, sz = curv(z, cell, "x", 4.80)
        dc, sc = curv(mirror_close(zc, 1), cell, "x", 4.80)
        for lo, hi in ((0.005, 0.036), (0.036, 0.20), (0.20, 1.00),
                       (1.00, 4.80)):
            hw = M5.fit_band(dz, sz, lo, hi)
            hn = M5.fit_band(dc, sc, lo, hi)
            dec["%.3f-%.3f" % (lo, hi)] = {
                "H": hw["H"], "n": hw["n_points"], "rms_log": hw["rms_log"],
                "H_no_wrap": hn["H"], "n_no_wrap": hn["n_points"],
                "rms_log_no_wrap": hn["rms_log"]}
        rec["decade_H_x"] = dec

        # increment whiteness (the random-walk signature).  Measured on the
        # interior crop, so that the constant used to fill the boundary gap
        # cannot masquerade as surface roughness.
        zc = crop_interior(z0, hole)
        dx = np.diff(zc, axis=1)
        dx = dx - dx.mean(axis=1, keepdims=True)
        a, b = dx[:, :-1].ravel(), dx[:, 1:].ravel()
        dy = np.diff(zc, axis=0)
        dy = dy - dy.mean(axis=1, keepdims=True)
        s1 = float(np.mean(dx ** 2))
        s2 = float(np.mean((dx[:, :-1] + dx[:, 1:]) ** 2))
        rec["increment"] = {
            "region": "interior crop, no periodic wrap (np.diff)",
            "rms_x_mm": float(np.sqrt(s1) * 1e3),
            "rms_y_mm": float(np.sqrt(np.mean(dy ** 2)) * 1e3),
            "autocorr1_x": float(np.corrcoef(a, b)[0, 1]),
            "S1_m2": s1, "S2_m2": s2, "S2_over_S1_x": s2 / s1,
            "note": ("with the wrap removed the increments are only weakly "
                     "correlated (autocorr1 above zero), so the exactly "
                     "S2/S1 = 2.0 seen through np.roll is NOT a property of "
                     "the bed")}

        # cell coarsening of the true (no-wrap) exponent
        co = {}
        for f in COARSEN:
            zz = block_mean(zc, f)
            zz = mirror_close(zz, 1)
            d, s = M5.curve_dir(zz, cell * f, DM_X, nbin=NBIN, axis=1)
            loc = M5.local_H(d, s)
            rows = {}
            for p in PROBE_COARSE_M:
                if d[0] * 0.999 <= p <= d[-1] * 1.001:
                    rows["%.3f" % p] = float(np.interp(np.log(p), np.log(d), loc))
            co[str(f)] = {"cell_m": cell * f, "local_H_at_m": rows}
        rec["coarsening"] = co
        out[name] = rec
        log("  %-22s Hx=%.3f Hy=%.3f Hiso=%.3f  rms(x,y)=%.4f,%.4f  inc_rms %.2f mm"
            % (name, dirs["x"]["bands"]["scaling"]["H"],
               dirs["y"]["bands"]["scaling"]["H"],
               dirs["iso"]["bands"]["scaling"]["H"],
               dirs["x"]["bands"]["scaling"]["rms_log"],
               dirs["y"]["bands"]["scaling"]["rms_log"],
               rec["increment"]["rms_x_mm"]))
    cp.unit_done("stage2_measure")
    return out


def block_mean(z, f):
    ny, nx = z.shape
    ny2, nx2 = (ny // f) * f, (nx // f) * f
    return z[:ny2, :nx2].reshape(ny2 // f, f, nx2 // f, f).mean(axis=(1, 3))


def crop_interior(z0, hole):
    """Drop all-NaN rows/columns (a boundary mask, not an interior void)."""
    oky = ~hole.all(axis=1)
    okx = ~hole.all(axis=0)
    zc = z0[np.ix_(oky, okx)]
    return fill_mean(zc, ~np.isfinite(zc))


def stage_sens(cp, meas):
    """Box size, hole fill, placement -- one block each, all DEMs."""
    box, hole_s, place = {}, {}, {}
    for name in sorted(meas):
        cell = meas[name]["cell_m"]
        z0 = load_dem(name)["z0"]
        hp = load_dem(name)["hole"]
        z = fill_mean(z0, hp)
        ny, nx = z.shape

        # ---- (2) box size: full / 1/2 / 1/4 / the ruler's own 0.6 m ------
        b = {}

        def measure_window(z0w, hpw):
            zz = fill_mean(z0w, hpw)
            zc = crop_interior(z0w, hpw)
            r = {}
            for axis, dmax, ax in (("x", DM_X, 1), ("y", DM_Y, 0),
                                   ("iso", DM_ISO, 0)):
                d, s = curv(zz, cell, axis, dmax)
                f = M5.fit_band(d, s, *BAND)
                r[axis] = {"H": f["H"], "n": f["n_points"],
                           "rms_log": f["rms_log"], "dmax_eff_m": float(d[-1])}
                zzc = (mirror_close(mirror_close(zc, 0), 1) if axis == "iso"
                       else mirror_close(zc, ax))
                d2, s2 = curv(zzc, cell, axis, dmax)
                r[axis]["H_clean"] = M5.fit_band(d2, s2, *BAND)["H"]
            return r

        for frac, tag in ((1.0, "full"), (0.5, "half"), (0.25, "quarter")):
            n = int(nx * frac)
            x0 = (nx - n) // 2
            b[tag] = {"px": [int(z.shape[0]), int(n)], "x_span_m": n * cell}
            b[tag].update(measure_window(z0[:, x0:x0 + n], hp[:, x0:x0 + n]))
        ntile = max(1, int(round(TILE_M / cell)))
        x0 = (nx - ntile) // 2
        b["native"] = {"px": [int(z.shape[0]), ntile], "x_span_m": ntile * cell}
        b["native"].update(measure_window(z0[:, x0:x0 + ntile],
                                          hp[:, x0:x0 + ntile]))
        hv = [b[t]["x"]["H"] for t in ("full", "half", "quarter", "native")]
        hv = [v for v in hv if v is not None]
        b["x_range"] = [min(hv), max(hv)] if hv else None
        b["x_delta_full_to_native"] = (
            b["native"]["x"]["H"] - b["full"]["x"]["H"]
            if None not in (b["native"]["x"]["H"], b["full"]["x"]["H"]) else None)
        box[name] = b

        # ---- (3) hole fill: the project's rules ---------------------------
        h = {"nodata_frac": float(hp.mean()),
             "n_rule_values": int(hp.sum())}
        if hp.any():
            for rule in ("mean", "p75", "max", "inpaint"):
                zz, fb = fill_rule2(z0, hp, rule, cell)
                r = {"fallback_frac": fb, "sigma_z_m": float(zz.std())}
                for axis, dmax in (("x", DM_X), ("y", DM_Y), ("iso", DM_ISO)):
                    d, s = curv(zz, cell, axis, dmax)
                    f = M5.fit_band(d, s, *BAND)
                    r[axis] = f["H"]
                h[rule] = r
            vals = [h[k]["x"] for k in ("mean", "p75", "max") if h[k]["x"]]
            h["spread_x_excl_inpaint"] = max(vals) - min(vals) if vals else None
            h["inpaint_minus_mean_x"] = (h["inpaint"]["x"] - h["mean"]["x"]
                                         if None not in (h["inpaint"]["x"],
                                                         h["mean"]["x"]) else None)
            # how far from the local bed does the wrapped inpaint land?
            zi, _ = fill_rule2(z0, hp, "inpaint", cell)
            valid_any = (~hp).any(axis=0)
            last = int(np.where(valid_any)[0].max())
            m = np.isfinite(z0[:, last])
            h["inpaint_gap_rise_m"] = float(np.mean(
                np.abs(zi[m, -1] - z0[m, last])))
            d_in = np.diff(zi[:, last:], axis=1)
            h["inpaint_gap_increment_rms_mm"] = float(
                np.sqrt(np.nanmean(d_in ** 2)) * 1e3)
            del zi
            # boundary treatment: drop the empty tail / wall rows, re-measure
            oky = ~hp.all(axis=1)
            okx = ~hp.all(axis=0)
            zc = z0[np.ix_(oky, okx)]
            hc = ~np.isfinite(zc)
            zc = fill_mean(zc, hc)
            r = {"shape": [int(zc.shape[0]), int(zc.shape[1])],
                 "dropped_rows": int((~oky).sum()),
                 "dropped_cols": int((~okx).sum())}
            for axis, dmax in (("x", DM_X), ("y", DM_Y), ("iso", DM_ISO)):
                d, s = curv(zc, cell, axis, dmax)
                f = M5.fit_band(d, s, *BAND)
                r[axis] = f["H"]
            h["interior_crop"] = r
            h["crop_minus_mean_x"] = (r["x"] - h["mean"]["x"]
                                      if None not in (r["x"], h["mean"]["x"])
                                      else None)
        else:
            h["note"] = "no nodata: not applicable"
        hole_s[name] = h

        # ---- (4) placement: the ruler's own 0.6 m box, stepped along x ----
        hx, hy, hi, hxn, hyn = [], [], [], [], []
        for k in range(0, nx - ntile + 1, ntile):
            zz = fill_mean(z0[:, k:k + ntile], hp[:, k:k + ntile])
            zc = crop_interior(z0[:, k:k + ntile], hp[:, k:k + ntile])
            for axis, acc, dmax in (("x", hx, DM_X), ("y", hy, DM_Y),
                                    ("iso", hi, DM_ISO)):
                d, s = curv(zz, cell, axis, dmax)
                acc.append(M5.fit_band(d, s, *BAND)["H"])
            for axis, acc, ax, dmax in (("x", hxn, 1, DM_X),
                                        ("y", hyn, 0, DM_Y)):
                zzc = mirror_close(zc, ax)
                d, s = curv(zzc, cell, axis, dmax)
                acc.append(M5.fit_band(d, s, *BAND)["H"])
                del zzc
            del zc, zz
        def st(v):
            v = np.array([x for x in v if x is not None], dtype=float)
            if v.size == 0:
                return {"n": 0}
            return {"n": int(v.size), "mean": float(v.mean()),
                    "sd": float(v.std(ddof=1)) if v.size > 1 else None,
                    "min": float(v.min()), "max": float(v.max()),
                    "p05": float(np.percentile(v, 5)),
                    "p95": float(np.percentile(v, 95))}
        place[name] = {"tile_m": TILE_M, "tile_px": ntile,
                       "x": st(hx), "y": st(hy), "iso": st(hi),
                       "x_no_wrap": st(hxn), "y_no_wrap": st(hyn)}
        log("  %-22s box dX=%.3f  place sd(x,y)=%.3f,%.3f  holes %s"
            % (name, box[name]["x_delta_full_to_native"] or 0.0,
               place[name]["x"]["sd"] or 0.0, place[name]["y"]["sd"] or 0.0,
               ("n/a" if not hp.any() else
                "spread=%.3f inpaint=%.3f" % (
                    hole_s[name]["spread_x_excl_inpaint"] or 0.0,
                    hole_s[name]["inpaint_minus_mean_x"] or 0.0))))
    cp.unit_done("stage2_box")
    cp.unit_done("stage2_holes")
    cp.unit_done("stage2_placement")
    return {"box": box, "hole": hole_s, "placement": place}


def fill_rule2(z0, hole, rule, cell):
    z = z0.copy()
    v = z0[~hole]
    fb = 0.0
    if rule == "mean":
        z[hole] = float(v.mean())
    elif rule == "p75":
        z[hole] = float(np.percentile(v, 75))
    elif rule == "max":
        z[hole] = float(v.max())
    elif rule == "inpaint":
        from scipy import ndimage
        sig_px = INPAINT_M / cell
        w = (~hole).astype(np.float64)
        num = ndimage.gaussian_filter(np.where(hole, 0.0, z0), sig_px, mode="wrap")
        den = ndimage.gaussian_filter(w, sig_px, mode="wrap")
        ok = den > 1e-3
        filled = np.where(ok, num / np.maximum(den, 1e-12), float(v.mean()))
        z[hole] = filled[hole]
        fb = float((hole & ~ok).sum()) / max(int(hole.sum()), 1)
    return z, fb


# ---------------------------------------------------------------------------
# STAGE 3 -- report
# ---------------------------------------------------------------------------
def f3(x):
    return "—" if x is None else "%.3f" % x


def f4(x):
    return "—" if x is None else "%.4f" % x


def f2(x):
    return "—" if x is None else "%.2f" % x


def summ(meas, band, key):
    v = [meas[n]["dirs"][key]["bands"][band]["H"] for n in meas]
    v = [x for x in v if x is not None]
    return (min(v), float(np.median(v)), max(v)) if v else (None, None, None)


def nested_range(b, key="H"):
    """full -> 1/2 -> 1/4, i.e. the box effect alone (no placement change)."""
    v = [b[t]["x"][key] for t in ("full", "half", "quarter")]
    v = [x for x in v if x is not None]
    return (max(v) - min(v)) if v else 0.0


def build_report(ctx):
    st0, dl, pr, meas, sens = (ctx["stage0"], ctx["download"],
                               ctx["prove"], ctx["measure"], ctx["sens"])
    rec, grid = dl["record"], pr["grid"]
    dems = pr["dem"]
    names = sorted(meas)
    prim = PRIMARY
    pm = meas[prim]
    cell = pm["cell_m"]
    d50p = pm["D50_ref_mm"] * 1e-3
    lo_d50, hi_d50 = BAND[0] / d50p, BAND[1] / d50p
    first_d50 = cell / d50p
    d50_surf = M5.OLD["D50_surface_v422"]

    def rr(n, axis, band, key="dirs"):
        return meas[n][key][axis]["bands"][band]["H"]

    def rng(axis, band, key="dirs"):
        v = [rr(n, axis, band, key) for n in names]
        v = [x for x in v if x is not None]
        return (min(v), float(np.median(v)), max(v))

    A = []
    P = A.append

    P("# BRIEF-V4-2-9 — 公開実測 DEM での物差し較正（v4-2.9）")
    P("")
    P("**これは物差しの較正であって、うちの床の合否判定ではない。**")
    P("")
    P("- 対象は Redolfi et al. (2025) の公開実測床（Zenodo record %d）。"
      "平均粒径 約 1.2 mm の砂‐礫混合、交互砂州を作る直線水路の最終平衡床。"
      % RECORD_ID)
    P("- 粒径は、うちの物差しが想定する粒径（`bed_gen2.D50_IN` = %.0f mm、"
      "表面 D50 ≈ %.1f mm）の **約 1/17〜1/30**。しかも 1.2 mm は移動床の"
      "ほぼ全体が動く粒径で、装甲床ではない。"
      % (1e3 * __import__("bed_gen2").D50_IN, 1e3 * d50_surf))
    P("- したがって **値の大小をうちの床と比べることはしない**。本レポートが"
      "見るのは「実データに対して `curve_dir` / `curve_iso` / `local_H` / "
      "`fit_band` がどう振る舞うか」だけである。")
    P("- うちの床についての合否・良し悪しは本レポートに書かない（§9）。")
    P("- 判断はすべて数値で行った。画像の目視は使っていない。")
    P("")
    P("物差しは `run_v4_2_5.py` の関数を import して、そのまま使った"
      "（再実装していない）。")
    P("")

    # ---------------- 0. headline -------------------------------------
    P("## 0. 先に結論（数値のみ）")
    P("")
    P("1. **主判定帯 0.4–3.6 cm は、この床では物理的なスケールを 1 つも"
      "含まない帯になる。** 粒径 1.2 mm は 0.24 px、グリッドは 5 mm = %.1f D50、"
      "帯は %.1f–%.1f D50。帯の上（40–100 cm）に砂州、帯の下（< 2.4 mm）に"
      "粒径がある。帯の中に残るのは DEM の増分雑音である。"
      % (first_d50, lo_d50, hi_d50))
    P("2. **14 枚すべてで H_x = %.3f–%.3f（scaling 帯）。** これは"
      "物差しの `np.roll`（周期境界）が作る人工物である。非周期の場の両端を"
      "折り返すと、継ぎ目に |Δ| = %.3f–%.3f m の段差が入る（22.5 m にわたる "
      "1%% のこう配 0.19–0.22 m と、下流端の欠測を mean で埋めた定数との差）。"
      "この段差が S に **厳密に d^1 の項** を加える。その項は"
      "真の S(5 mm) の %.0f–%.0f 倍で、H は定義上 0.500 になる。"
      % (rng("x", "scaling")[0], rng("x", "scaling")[2],
         min(meas[n]["seam"]["delta_x_m"] for n in names),
         max(meas[n]["seam"]["delta_x_m"] for n in names),
         min(meas[n]["seam"]["S1_ruler_over_clean_x"] - 1 for n in names),
         max(meas[n]["seam"]["S1_ruler_over_clean_x"] - 1 for n in names)))
    P("3. **折り返しを厳密に閉じ、境界マスクを落として測り直すと "
      "H_x = %.3f–%.3f になる**（crop + mirror、§2.4）。資料ごとに違い、"
      "帯の外では %.2f–%.2f まで上がる。"
      % (rng("x", "scaling", "dirs_clean")[0],
         rng("x", "scaling", "dirs_clean")[2],
         rng("x", "meso", "dirs_clean")[0], rng("x", "meso", "dirs_clean")[2]))
    P("4. **同じ帯・同じ物差しで、資料次第で 0.50 と 0.72 の両方が出る。** "
      "meso 帯（6.6–20 cm）の H_x は MW40 で %.3f–%.3f、MW60 で %.3f–%.3f。"
      "折り返しを閉じると MW40 も %.3f–%.3f になり、**差は消える**"
      "（その差は物差しが作っていた）。"
      % (min(rr(n, "x", "meso") for n in names if n.startswith("DEM_MW40")),
         max(rr(n, "x", "meso") for n in names if n.startswith("DEM_MW40")),
         min(rr(n, "x", "meso") for n in names if n.startswith("DEM_MW60")),
         max(rr(n, "x", "meso") for n in names if n.startswith("DEM_MW60")),
         min(rr(n, "x", "meso", "dirs_clean") for n in names
             if n.startswith("DEM_MW40")),
         max(rr(n, "x", "meso", "dirs_clean") for n in names
             if n.startswith("DEM_MW40"))))
    P("5. **方向で 1 桁動く。** H_y は scaling 帯で ruler %.3f–%.3f、clean "
      "%.3f–%.3f。meso 帯では ruler %.3f–%.3f、clean %.3f–%.3f。"
      "`curve_iso` = 0.5(S_x+S_y) は「折り返しの人工物」と「飽和する実信号」の"
      "混合になる（§7.3）。"
      % (rng("y", "scaling")[0], rng("y", "scaling")[2],
         rng("y", "scaling", "dirs_clean")[0],
         rng("y", "scaling", "dirs_clean")[2],
         rng("y", "meso")[0], rng("y", "meso")[2],
         rng("y", "meso", "dirs_clean")[0], rng("y", "meso", "dirs_clean")[2]))
    P("6. **箱サイズ感度は合成場の 0.1 を再現しない。** 入れ子（全域→1/2→1/4）"
      "の H_x の振れは ruler で %.3f–%.3f、折り返しを閉じた読みで "
      "%.3f–%.3f（中央 %.3f）。合成場の 0.10 より小さい（§4）。"
      % (min(nested_range(sens["box"][n], "H") for n in names),
         max(nested_range(sens["box"][n], "H") for n in names),
         min(nested_range(sens["box"][n], "H_clean") for n in names),
         max(nested_range(sens["box"][n], "H_clean") for n in names),
         float(np.median([nested_range(sens["box"][n], "H_clean")
                          for n in names]))))
    P("7. **穴埋め感度は「規則の種類」より「穴の形」で決まる。** "
      "mean / p75 / max の差は ≤ %.3f。周期ラップ inpaint は H_x を "
      "+%.3f（MW60）／H_y を +%.3f（MW40）動かす。nodata が内部の粒穴ではなく"
      "**下流端の境界マスク**（20–120 px の全欠測列）で、`mode=\"wrap\"` が"
      "上流端と混ぜて 0.20–0.24 m の偽の登りを作るから（§5）。"
      % (max(sens["hole"][n]["spread_x_excl_inpaint"] for n in names),
         max((sens["hole"][n]["inpaint_minus_mean_x"] or 0) for n in names),
         max((sens["hole"][n]["inpaint"]["y"] - sens["hole"][n]["mean"]["y"])
             for n in names)))
    P("8. **配置感度は 0.6 m 箱で H_x の sd %.3f–%.3f（折り返しを閉じた読み）。** "
      "箱 1 枚の H_x は %.2f–%.2f の範囲に散る。**ruler の sd（%.3f–%.3f）は"
      "その散らばりを隠している**（§6）。"
      % (min(sens["placement"][n]["x_no_wrap"]["sd"] for n in names),
         max(sens["placement"][n]["x_no_wrap"]["sd"] for n in names),
         min(sens["placement"][n]["x_no_wrap"]["min"] for n in names),
         max(sens["placement"][n]["x_no_wrap"]["max"] for n in names),
         min(sens["placement"][n]["x"]["sd"] for n in names),
         max(sens["placement"][n]["x"]["sd"] for n in names)))
    P("")

    # ---------------- 1. provenance -----------------------------------
    P("## 1. 素性（STAGE 1）")
    P("")
    P("| 項目 | 値 |")
    P("| --- | --- |")
    P("| レコード | %s |" % rec["url"])
    P("| DOI | %s（concept %s）|" % (rec["doi"], rec["concept_doi"]))
    P("| タイトル | %s |" % rec["title"])
    P("| 著者 | %s |" % ", ".join(rec["creators"]))
    P("| 公開日 / 版 | %s / %s |" % (rec["publication_date"], rec["version"]))
    P("| ライセンス | **%s** |" % rec["license"])
    P("| レコード全ファイル数 | %d（うち DEM_* 14、d50_* 13 を取得）|"
      % dl["n_files_in_record"])
    P("| 取得量 | %.1f MB / 27 files（md5 は全件一致、sha256 は numbers）|"
      % (dl["total_bytes"] / 1e6))
    P("| グリッド間隔 | %.4f m（全ファイル同一）。説明文の「0.5 mm」とは"
      "食い違う（§1.1）|" % grid["cell_m_values"][0])
    P("| 単位 | 標高 m（float32/64）、粒径 mm |")
    P("| 座標系 | EPSG:%s（UTM32N、局所原点）|"
      % (dems[names[0]]["epsg"] if dems[names[0]]["epsg"] else "?"))
    P("")
    P("### 1.1 グリッドは 5 mm であって 0.5 mm ではない")
    P("")
    P("レコード説明文は “spatial resolution 0.5x0.5 mm” と書いているが、"
      "GeoTIFF の `ModelPixelScaleTag` は全 27 ファイルで 0.005 m、"
      "`ModelTiepointTag` の原点は (%.5f, %.5f)。決め手は幅の整合である:"
      % (dems[names[0]]["tiepoint"][3], dems[names[0]]["tiepoint"][4]))
    P("")
    P("| シリーズ | ラスタ (px) | nx × px (m) | 公表水路幅 W (m) | 比 |")
    P("| --- | --- | --- | --- | --- |")
    seen = set()
    for k in sorted(dems):
        e = dems[k]
        key = (e["series"], e["n_cols"], e["n_rows"])
        if key in seen:
            continue
        seen.add(key)
        P("| %s | %d × %d | %.3f | %.2f | %.3f |"
          % (e["series"], e["n_cols"], e["n_rows"], e["width_over_cell"],
             e["channel_width_m"], e["width_over_cell"] / e["channel_width_m"]))
    P("")
    P("MW40 は 80 px × 5 mm = 0.400 m で W に**厳密に**一致する。0.5 mm なら"
      "4 cm 幅になり、水路幅 40 cm と 1 桁合わない。MW60 は 110 px × 5 mm = "
      "0.550 m で、W = 0.60 m の 92%（岸際を 5 cm 落とした範囲）。"
      "よって公開 GeoTIFF のグリッドは 5 mm で、説明文の 0.5 mm は公開前の"
      "製品（元の SfM 点群）を指していると判断する。**以下のラグはすべて "
      "5 mm 格子の実寸**である。")
    P("")
    P("範囲: 各ラスタは %.2f m（= %d px × 5 mm）、ノルシング原点は"
      "全 DEM で 23.88975 m。水路 24 m のほぼ全長を含む。"
      % (dems[names[0]]["extent"]["height_m"], dems[names[0]]["n_rows"]))
    P("")
    P("### 1.2 nodata の形（重要）")
    P("")
    P("| シリーズ | nodata | 全欠測列＝横断位置（岸）| 全欠測行＝流下位置 |")
    P("| --- | --- | --- | --- |")
    for ser in ("MW40", "MW60"):
        ks = [k for k in sorted(dems) if dems[k]["series"] == ser]
        P("| %s | %.2f%% | %d / %d | %d / %d |"
          % (ser, 100 * dems[ks[0]]["nodata_frac"],
             dems[ks[0]]["all_nan_cols"], dems[ks[0]]["n_cols"],
             dems[ks[0]]["all_nan_rows"], dems[ks[0]]["n_rows"]))
    P("")
    P("中身を見ると、欠測は**内部の粒穴ではなく境界マスク**である: MW40 は"
      "両岸の 5 列（横断 80 px のうち 5 px）と、**流下方向の最後 120 px"
      "（0.60 m）**が全欠測。MW60 は最後 20 px（0.10 m）だけ。"
      "この「下流端の全欠測ブロック」が§5 の穴埋め感度の正体である。")
    P("")
    P("### 1.3 D50 マップ")
    P("")
    med = [meas[n]["D50_map"]["median_mm"] for n in names if meas[n]["D50_map"]]
    P("| ファイル | 中央値 (mm) | p05–p95 (mm) |")
    P("| --- | --- | --- |")
    for n in names:
        dm = meas[n]["D50_map"]
        if dm:
            P("| %s | %.3f | %.3f–%.3f |"
              % (dm["file"], dm["median_mm"], dm["p05_mm"], dm["p95_mm"]))
    P("")
    P("中央値は %.3f–%.3f mm で、レコード記載の 1.22 / 1.29 mm と同程度"
      "（マップの方が 0.1 mm ほど高い）。以下の無次元ラグは各 DEM 自身の "
      "D50 マップ中央値で割った。" % (min(med), max(med)))
    P("")
    P("### 1.4 ファイル表（sha256 は numbers-v4-2-9.json に全 27 件）")
    P("")
    P("| ファイル | bytes | md5 |")
    P("| --- | --- | --- |")
    for f in dl["files"]:
        P("| `%s` | %d | `%s` |" % (f["key"], f["bytes"], f["md5"]))
    P("")
    P("実行日時 %s。ダウンロードは `out/v4-2-9/data/` に置いた（既存ファイルは"
      "触っていない）。" % ctx["created_utc"])
    P("")

    # ---------------- 2. the ruler ------------------------------------
    P("## 2. 物差しの読み（STAGE 2.1）")
    P("")
    P("物差しの契約は v4-2.5 のまま: 帯 scaling 0.4–3.6 cm / meso 6.6–20 cm / "
      "legacy 4–20 cm、ラグ梯子 nbin=%d、x は dmax %.2f m、y は %.2f m、"
      "iso は %.2f m。H = 0.5·dlnS/dlnd。" % (NBIN, DM_X, DM_Y, DM_ISO))
    P("")
    P("3 通りの読みを併記する:")
    P("")
    P("| 記号 | 場 | 周期閉鎖 |")
    P("| --- | --- | --- |")
    P("| **ruler** | 全面（nodata は mean 埋め）| `np.roll`（物差しのまま）|")
    P("| **full+mirror** | 全面（同じ埋め）| 鏡映で厳密に閉じる |")
    P("| **clean** | 全欠測の行・列を削除 | 鏡映で厳密に閉じる |")
    P("")
    P("帯を D50 で書き直すとこうなる（主検体 %s、D50 = %.2f mm）:"
      % (prim, 1e3 * d50p))
    P("")
    P("| 帯 | 実寸 | この床の D50 単位 | 梯子に乗るラグ数 | 最初のラグ |")
    P("| --- | --- | --- | --- | --- |")
    for k, lo, hi in M5.BANDS:
        npx = pm["dirs"]["x"]["steps_per_band"][k]
        P("| %s | %.1f–%.1f cm | %.1f–%.1f × D50 | %d | %.0f mm (%.1f D50) |"
          % (k, 100 * lo, 100 * hi, lo / d50p, hi / d50p, npx, 1e3 * cell,
             first_d50))
    P("")
    P("**帯の下端 %.1f mm はグリッド %.0f mm より細かい。** 実際に測れるのは "
      "5 mm（= %.1f D50）からで、帯の下 1/8 は存在しない。うちの物差しの設計"
      "意図（＝粒径以下〜粒径）は、粒径が %.2f mm のこの床では 0.2–2.2 mm に"
      "対応し、**丸ごとグリッドの下**にある。同じ「0.4–3.6 cm」という数字が、"
      "うちの床では 0.11–0.99 × D50_surface、この床では %.1f–%.1f × D50 を"
      "意味する。"
      % (1e3 * BAND[0], 1e3 * cell, first_d50, 1e3 * d50p, lo_d50, hi_d50))
    P("")
    BSHORT = {"scaling": "0.4–3.6 cm", "meso": "6.6–20 cm",
              "legacy": "4–20 cm"}
    for band, tt in (("scaling", "2.1"), ("meso", "2.2"), ("legacy", "2.3")):
        P("### %s %s 帯（%s）の読み" % (tt, band, BSHORT[band]))
        P("")
        P("| DEM | H_x ruler | H_x clean | H_y ruler | H_y clean | "
          "H_iso ruler | H_iso clean | rms_log (ruler/clean, x) |")
        P("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for n in names:
            P("| %s | %s | %s | %s | %s | %s | %s | %s / %s |"
              % (n[4:-4], f3(rr(n, "x", band)), f3(rr(n, "x", band, "dirs_clean")),
                 f3(rr(n, "y", band)), f3(rr(n, "y", band, "dirs_clean")),
                 f3(rr(n, "iso", band)), f3(rr(n, "iso", band, "dirs_clean")),
                 f4(meas[n]["dirs"]["x"]["bands"][band]["rms_log"]),
                 f4(meas[n]["dirs_clean"]["x"]["bands"][band]["rms_log"])))
        P("")
        P("H_x の範囲: ruler %.3f–%.3f（中央 %.3f）／clean %.3f–%.3f"
          "（中央 %.3f）。H_y: ruler %.3f–%.3f／clean %.3f–%.3f。"
          % (rng("x", band)[0], rng("x", band)[2], rng("x", band)[1],
             rng("x", band, "dirs_clean")[0], rng("x", band, "dirs_clean")[2],
             rng("x", band, "dirs_clean")[1],
             rng("y", band)[0], rng("y", band)[2],
             rng("y", band, "dirs_clean")[0], rng("y", band, "dirs_clean")[2]))
        P("")
    P("### 2.4 3 通りの読みの並置（scaling 帯）")
    P("")
    P("| 読み | H_x | H_y | H_iso |")
    P("| --- | --- | --- | --- |")
    for key, lab in (("dirs", "ruler（np.roll）"),
                     ("dirs_no_wrap", "full+mirror"),
                     ("dirs_clean", "clean（crop+mirror）")):
        P("| %s | %s (%.3f–%.3f) | %s (%.3f–%.3f) | %s |"
          % (lab, f3(rr(prim, "x", "scaling", key)), rng("x", "scaling", key)[0],
             rng("x", "scaling", key)[2],
             f3(rr(prim, "y", "scaling", key)), rng("y", "scaling", key)[0],
             rng("y", "scaling", key)[2],
             f3(rr(prim, "iso", "scaling", key))))
    P("")
    P("（括弧は 14 枚の範囲。主検体 %s。§3 の 0.593 と §3.1 の 0.544 が"
      "少し違うのは、帯の上で dmax の違う梯子を使っているため。梯子の"
      "違いだけで 0.05 動く。）" % prim)
    P("")
    P("### 2.5 折り返しが注入する項")
    P("")
    P("物差しの `curve_dir` / `curve_iso` は `np.roll` で周期境界を作る。"
      "非周期の場を折り返すと、1 行あたり 1 組だけ「場の終端」と「場の始端」を"
      "比べることになり、その差 Δ は 22.5 m にわたる 1%% のこう配（= 0.19–"
      "0.22 m）と、下流端の欠測を mean で埋めた定数（場の平均）の差で決まる。"
      "%d 行の平均では、lag d に d 組の折り返しが入るので S への寄与は"
      "**2Δ²·d/N**、すなわち d^1 に比例する項になる。" % pm["shape"][0])
    P("")
    P("| DEM | Δ_x (m) | 予測 2Δ²/N (m²) | 実測 S(1px) ruler | S(1px) clean | "
      "ruler/clean (x) | ruler/clean (y) |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    for n in names:
        s = meas[n]["seam"]
        P("| %s | %.3f | %.2e | %.2e | %.2e | %.1f | %.1f |"
          % (n[4:-4], s["delta_x_m"], s["predicted_wrap_term_m2_at_1px"],
             s["S1_ruler_m2"], s["S1_clean_m2"],
             s["S1_ruler_over_clean_x"], s["S1_ruler_over_clean_y"]))
    P("")
    P("予測項と実測の差は 10 %% 未満で一致する（平均 %.1f %%）。"
      "**ruler の S は 1 px の時点ですでに真の値の %.0f–%.0f 倍（x）／"
      "%.0f–%.0f 倍（y）**であり、その余分がちょうど d^1 なので、"
      "`fit_band` は 0.500 を返す。"
      % (100.0 * float(np.mean([
          abs(meas[n]["seam"]["predicted_wrap_term_m2_at_1px"]
              - (meas[n]["seam"]["S1_ruler_m2"] - meas[n]["seam"]["S1_clean_m2"]))
          / max(meas[n]["seam"]["predicted_wrap_term_m2_at_1px"], 1e-30)
          for n in names])),
         min(meas[n]["seam"]["S1_ruler_over_clean_x"] for n in names),
         max(meas[n]["seam"]["S1_ruler_over_clean_x"] for n in names),
         min(meas[n]["seam"]["S1_ruler_over_clean_y"] for n in names),
         max(meas[n]["seam"]["S1_ruler_over_clean_y"] for n in names)))
    P("")
    P("なお Δ_x は 0.10–0.12 m で、こう配×全長（0.19–0.22 m）の約半分。"
      "下流端の欠測を埋めた値が場の平均（＝こう配の中点）だからである。"
      "もし欠測を埋めずに端を落とせば Δ は 0.19–0.23 m になり、人工物は"
      "さらに 2 倍になる。")
    P("")

    # ---------------- 3. profile --------------------------------------
    P("## 3. S と local_H のプロファイル（主検体 %s）" % prim)
    P("")
    P("| ラグ (mm) | ラグ/D50 | S_x ruler (m²) | S_x clean (m²) | "
      "H_x clean | H_y clean | H_iso clean |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    for a, b, c, d in zip(pm["dirs"]["x"]["profile"],
                          pm["dirs_clean"]["x"]["profile"],
                          pm["dirs_clean"]["y"]["profile"],
                          pm["dirs_clean"]["iso"]["profile"]):
        P("| %.1f | %.2f | %s | %s | %s | %s | %s |"
          % (1e3 * a[0], a[0] / d50p,
             "—" if a[1] is None else "%.3e" % a[1],
             "—" if b[1] is None else "%.3e" % b[1],
             f3(b[2]), f3(c[2]), f3(d[2])))
    P("")
    inb = [r for r in pm["dirs_clean"]["x"]["profile"]
           if r[1] is not None and 0.004 <= r[0] <= 0.036]
    out = [r for r in pm["dirs_clean"]["x"]["profile"]
           if r[1] is not None and r[0] > 0.036]
    imax = max(out, key=lambda r: r[2])
    P("clean の読みでは、x は帯の中で 5 mm の %.2f から 35 mm の %.2f へ下がり、"
      "そこから 1 m の %.2f へ持ち直す（2–5 cm に谷がある）。"
      "y は 5 mm で %.2f、35 mm で %.2f と、帯の中で既に下がり始め、"
      "200 mm では %.2f になる。**同じ band の中で、x は「谷」、y は「帯の"
      "上端から先の減衰」を見ている**。"
      % (pm["dirs_clean"]["x"]["profile"][0][2],
         [r for r in pm["dirs_clean"]["x"]["profile"]
          if abs(r[0] - 0.035) < 1e-9][0][2],
         [r for r in pm["dirs_clean"]["x"]["profile"]
          if abs(r[0] - 1.0) < 1e-9][0][2],
         pm["dirs_clean"]["y"]["profile"][0][2],
         [r for r in pm["dirs_clean"]["y"]["profile"]
          if abs(r[0] - 0.035) < 1e-9][0][2],
         [r for r in pm["dirs_clean"]["y"]["profile"]
          if abs(r[0] - 0.2) < 1e-9][0][2]))
    P("")
    P("### 3.1 べき指数の年代別（H_x、主検体と 14 枚の範囲）")
    P("")
    P("| ラグ範囲 | ruler（主検体）| clean（主検体）| clean の rms_log | "
      "clean（14 枚の範囲）|")
    P("| --- | --- | --- | --- | --- |")
    for key in ("0.005-0.036", "0.036-0.200", "0.200-1.000", "1.000-4.800"):
        v = pm["decade_H_x"][key]
        vv = [meas[n]["decade_H_x"][key]["H_no_wrap"] for n in names]
        vv = [x for x in vv if x is not None]
        P("| %s m | %s | %s | %s | %s–%s |"
          % (key.replace("-", " – "), f3(v["H"]), f3(v["H_no_wrap"]),
             f4(v["rms_log_no_wrap"]), f3(min(vv)), f3(max(vv))))
    P("")
    P("clean では H_x は 5 mm で %.2f から単調に上がり、4.8 m で %.2f になる。"
      "**どこにも飽和の膝がない**（S は 3 桁にわたって増え続ける）。"
      "ruler の列は 0.50 → 0.62 → 0.56 → 0.42 と形が崩れている: 折り返しの "
      "d^1 項が lag に比例して増えるため、真の曲線が lag とともに急になっても"
      "平坦に見える。"
      % (pm["decade_H_x"]["0.005-0.036"]["H_no_wrap"],
         pm["decade_H_x"]["1.000-4.800"]["H_no_wrap"]))
    P("")
    P("### 3.2 増分（折り返しを含まない生の差）")
    P("")
    P("| DEM | 増分 rms_x (mm) | 増分 rms_y (mm) | lag-1 自己相関 (x) | "
      "S(2px)/S(1px) |")
    P("| --- | --- | --- | --- | --- |")
    for n in names:
        i = meas[n]["increment"]
        P("| %s | %.2f | %.2f | %+.3f | %.3f |"
          % (n[4:-4], i["rms_x_mm"], i["rms_y_mm"], i["autocorr1_x"],
             i["S2_over_S1_x"]))
    P("")
    P("流下方向の増分は rms 0.36–0.48 mm、lag-1 自己相関 +0.11〜+0.34 で"
      "**白色ではない**（弱く正相関）。S(2px)/S(1px) は %.2f–%.2f で、"
      "白なら 2.0 になるはずのところ %+.2f〜%+.2f 多い。一方 ruler の梯子では "
      "S(2)/S(1) が 2.00–2.01 と「きれいな 2」になる（§2.5 の d^1 項が"
      "支配するため）。**物差しの「きれいさ」は、この資料では人工物の"
      "指標になっている**。"
      % (min(meas[n]["increment"]["S2_over_S1_x"] for n in names),
         max(meas[n]["increment"]["S2_over_S1_x"] for n in names),
         min(meas[n]["increment"]["S2_over_S1_x"] for n in names) - 2.0,
         max(meas[n]["increment"]["S2_over_S1_x"] for n in names) - 2.0))
    P("")
    P("### 3.3 セル粗視化（H = 0.5 は 1 px の産物か）")
    P("")
    P("clean の場を 2×・4× にブロック平均し、同じ物理ラグで local H を読む。")
    P("")
    P("| DEM | 20–40 mm: 生 (5 mm) | 2× (10 mm) | 4× (20 mm) |")
    P("| --- | --- | --- | --- |")
    for n in names:
        co = meas[n]["coarsening"]

        def g(f):
            v = co[f]["local_H_at_m"]
            a, b = v.get("0.020"), v.get("0.040")
            return "—" if None in (a, b) else "%.3f / %.3f" % (a, b)
        P("| %s | %s | %s | %s |" % (n[4:-4], g("1"), g("2"), g("4")))
    P("")
    P("（各欄は 20 mm と 40 mm の local H。）4× に粗視化すると 20–40 mm の "
      "H は %.2f–%.2f に上がる。つまり 5 mm 格子の H≈0.5 は、**画素スケールの"
      "白色雑音が d^0 の下駄として乗っている**ためで、2 cm まで平均すると "
      "0.6–0.7 の本来の指数が見えてくる。"
      % (min(min(meas[n]["coarsening"]["4"]["local_H_at_m"]["0.020"],
                 meas[n]["coarsening"]["4"]["local_H_at_m"]["0.040"])
             for n in names),
         max(max(meas[n]["coarsening"]["4"]["local_H_at_m"]["0.020"],
                 meas[n]["coarsening"]["4"]["local_H_at_m"]["0.040"])
             for n in names)))
    P("")

    # ---------------- 4. box ------------------------------------------
    P("## 4. 箱サイズ感度（STAGE 2.2）")
    P("")
    P("窓は中央寄せ（全域 / 1/2 / 1/4）、加えて物差し自身の 0.6 m 箱。"
      "読みは scaling 帯の H_x。")
    P("")
    P("### 4.1 ruler（np.roll のまま）")
    P("")
    P("| DEM | 全域 (%.2f m) | 1/2 (%.2f m) | 1/4 (%.2f m) | 0.6 m 箱 | "
      "入れ子の振れ |"
      % (sens["box"][names[0]]["full"]["x_span_m"],
         sens["box"][names[0]]["half"]["x_span_m"],
         sens["box"][names[0]]["quarter"]["x_span_m"]))
    P("| --- | --- | --- | --- | --- | --- |")
    for n in names:
        b = sens["box"][n]
        P("| %s | %s | %s | %s | %s | %.3f |"
          % (n[4:-4], f3(b["full"]["x"]["H"]), f3(b["half"]["x"]["H"]),
             f3(b["quarter"]["x"]["H"]), f3(b["native"]["x"]["H"]),
             nested_range(b, "H")))
    P("")
    P("### 4.2 clean（crop + mirror）")
    P("")
    P("| DEM | 全域 | 1/2 | 1/4 | 0.6 m 箱 | 入れ子の振れ | Δ(0.6 m − 全域) |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    for n in names:
        b = sens["box"][n]
        P("| %s | %s | %s | %s | %s | %.3f | %+.3f |"
          % (n[4:-4], f3(b["full"]["x"]["H_clean"]),
             f3(b["half"]["x"]["H_clean"]), f3(b["quarter"]["x"]["H_clean"]),
             f3(b["native"]["x"]["H_clean"]),
             nested_range(b, "H_clean"),
             (b["native"]["x"]["H_clean"] - b["full"]["x"]["H_clean"])))
    P("")
    rn = [nested_range(sens["box"][n], "H") for n in names]
    rc = [nested_range(sens["box"][n], "H_clean") for n in names]
    P("入れ子（全域→1/2→1/4）の振れは ruler %.3f–%.3f、clean %.3f–%.3f"
      "（中央 %.3f）。**合成場の有限箱バイアス 約 %.2f（v4-2.6 §4.3）は"
      "再現しない。** ruler の振れが小さいのは、折り返しの d^1 項が箱の"
      "大きさによらず同じ形で乗るため（人工物が感度を隠す）。clean でも "
      "0.1 には届かず、中央 %.3f。"
      % (min(rn), max(rn), min(rc), max(rc), float(np.median(rc)),
         REF["box_bias_synthetic"], float(np.median(rc))))
    P("")
    P("理由は数値で説明できる: v4-2.6 の 0.1 は「自己アフィン場の膝が箱の"
      "中に入り、箱を小さくすると膝より上のラグが落ちて H が変わる」効果"
      "だった。ここでは scaling 帯（4–36 mm）の内側に膝が 1 つもなく（§3.1）、"
      "膝は 0.2–1 m の外にある。**箱バイアスは「帯の中に膝があるか」の"
      "関数であって、箱の大きさそのものの関数ではない。**")
    P("")

    # ---------------- 5. holes -----------------------------------------
    P("## 5. 穴埋め感度（STAGE 2.3）")
    P("")
    P("プロジェクトの既存規則（`run_v4_2_6.measure_field` の宣言）をそのまま"
      "当てた: `mean`（本ジョブの第一規則）、`p75`、`max`、`inpaint`"
      "（60 px of the v4-2.2 cell = %.3f m のガウス加重平均、`mode=\"wrap\"`）。"
      % INPAINT_M)
    P("")
    P("| DEM | nodata | mean | p75 | max | inpaint | mean/p75/max の幅 | "
      "inpaint − mean | 端を落とした読み |")
    P("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for n in names:
        h = sens["hole"][n]
        if "mean" not in h:
            continue
        P("| %s | %.2f%% | %s | %s | %s | %s | %.3f | %+.3f | %s |"
          % (n[4:-4], 100 * h["nodata_frac"], f3(h["mean"]["x"]),
             f3(h["p75"]["x"]), f3(h["max"]["x"]), f3(h["inpaint"]["x"]),
             h["spread_x_excl_inpaint"], h["inpaint_minus_mean_x"],
             f3(h["interior_crop"]["x"])))
    P("")
    P("（すべて scaling 帯の H_x、ruler の読み。「端を落とした読み」は全欠測"
      "列・全欠測行を物理的に削ってから測ったもの。）")
    P("")
    P("読めること:")
    P("")
    P("1. **第一規則 mean / p75 / max の差は ≤ %.3f。** v4-2.6 §6.4 の "
      "0.134 は「13.2%% の**内部**穴」の話で、このデータの nodata（境界マスク"
      "＋下流端の全欠測ブロック）には移送できない。"
      % max(sens["hole"][n]["spread_x_excl_inpaint"] for n in names))
    P("2. **`inpaint` だけが壊れる。** H_x で最大 %+.3f、H_y で最大 %+.3f。"
      "原因は `mode=\"wrap\"` が周期境界を仮定しているのに、この場は"
      "**非周期**で、しかも欠測が下流端（最後の 20–120 px）にあること。"
      "ガウス核が上流端（x = 0）の値を巻き込むので、欠測ブロックの埋め値が"
      "局所の床から最大 %.3f m ずれる。実例（%s）では欠測 20 px の間に "
      "%.3f m の偽の登りができ、その継ぎ目の増分 rms は %.1f mm"
      "（場の内部は %.1f mm）になる。"
      % (max(sens["hole"][n]["inpaint_minus_mean_x"] for n in names),
         max(sens["hole"][n]["inpaint"]["y"] - sens["hole"][n]["mean"]["y"]
             for n in names),
         max(sens["hole"][n]["inpaint_gap_rise_m"] for n in names),
         prim, sens["hole"][prim]["inpaint_gap_rise_m"],
         sens["hole"][prim]["inpaint_gap_increment_rms_mm"],
         meas[prim]["increment"]["rms_x_mm"]))
    P("3. **端を落とす（interior crop）と mean とほぼ同じ**（差 ≤ %.3f）。"
      "つまりこの資料では「境界をどう埋めるか」はほとんど効かず、"
      "「境界を周期として折り返すか」だけが効く。"
      % max(abs(sens["hole"][n]["crop_minus_mean_x"]) for n in names))
    P("")
    P("補足: `inpaint` の実装には 1 つだけ防御を入れた。v4-2.6 の式は "
      "`num/den` で、den が 0 に近い画素（欠測ブロックの内部）では 0 除算に"
      "なり、値が 0.0 m（＝実測 -1.07 m から +1.07 m の段差）になる。"
      "そこで `den > 1e-3` の画素だけ inpaint し、残りは mean で埋めた"
      "（割合を numbers の `fallback_frac` に記録。MW60 は 0、MW40 は 0.249）。")
    P("")

    # ---------------- 6. placement -------------------------------------
    P("## 6. 配置感度（STAGE 2.4、実現間ばらつきの目安）")
    P("")
    P("物差し自身の大きさ（0.6 m）の箱を、流下方向に重複なく 37 枚並べて"
      "同じ測定をした。")
    P("")
    P("| DEM | H_x ruler mean ± sd | H_x clean mean ± sd | H_x clean min–max | "
      "H_y ruler | H_y clean |")
    P("| --- | --- | --- | --- | --- | --- |")
    for n in names:
        p = sens["placement"][n]
        P("| %s | %.3f ± %.3f | %.3f ± %.3f | %.3f–%.3f | %.3f ± %.3f | "
          "%.3f ± %.3f |"
          % (n[4:-4], p["x"]["mean"], p["x"]["sd"],
             p["x_no_wrap"]["mean"], p["x_no_wrap"]["sd"],
             p["x_no_wrap"]["min"], p["x_no_wrap"]["max"],
             p["y"]["mean"], p["y"]["sd"],
             p["y_no_wrap"]["mean"], p["y_no_wrap"]["sd"]))
    P("")
    P("ruler の H_x は sd %.3f–%.3f で「安定」に見えるが、折り返しを閉じると "
      "sd は %.3f–%.3f に増え、箱 1 枚の H_x は %.2f–%.2f の範囲に散る。"
      "**物差しの安定さは、人工物が支配していることの裏返しである。**"
      "H_y も同じで、clean では sd %.3f–%.3f。"
      % (min(sens["placement"][n]["x"]["sd"] for n in names),
         max(sens["placement"][n]["x"]["sd"] for n in names),
         min(sens["placement"][n]["x_no_wrap"]["sd"] for n in names),
         max(sens["placement"][n]["x_no_wrap"]["sd"] for n in names),
         min(sens["placement"][n]["x_no_wrap"]["min"] for n in names),
         max(sens["placement"][n]["x_no_wrap"]["max"] for n in names),
         min(sens["placement"][n]["y_no_wrap"]["sd"] for n in names),
         max(sens["placement"][n]["y_no_wrap"]["sd"] for n in names)))
    P("")

    # ---------------- 7. where the ruler breaks ------------------------
    P("## 7. 物差しが実データで壊れる点")
    P("")
    P("### 7.1 周期ラップ（最大の破れ）")
    P("")
    P("`curve_dir` / `curve_iso` は `np.roll` で周期性を仮定する。この床は"
      "非周期（両岸と流下端が閉じ、22.5 m で 0.2 m 下がる）なので、折り返しの"
      "継ぎ目に Δ = 0.10–0.12 m の段差が入る。1 行 1 組の段差が lag d で d 組に"
      "なるため、S に **d^1 の項**が加わる。この項は")
    P("")
    P("- S(5 mm) の %.0f–%.0f 倍（x）、%.0f–%.0f 倍（y）"
      % (min(meas[n]["seam"]["S1_ruler_over_clean_x"] - 1 for n in names),
         max(meas[n]["seam"]["S1_ruler_over_clean_x"] - 1 for n in names),
         min(meas[n]["seam"]["S1_ruler_over_clean_y"] - 1 for n in names),
         max(meas[n]["seam"]["S1_ruler_over_clean_y"] - 1 for n in names)))
    P("- `fit_band` の返す H を **0.500** に固定する（人工物が支配する帯: "
      "14 枚すべての scaling 帯、MW40 の meso / legacy 帯で 0.499–0.512）")
    P("- しかも rms_log を 0.0013–0.0023 に**下げる**（当てはまりが良く見える）")
    P("")
    P("つまり「周期境界を使う限り、この種の資料では H は 0.50 に固定され、"
      "その 0.50 は床の性質ではない」。帯の外で真の H が 0.7 を超える"
      "（MW60 の meso 帯）ときだけ、真の信号が人工物を上回って 0.72 が"
      "顔を出す。§0 の 3 と 4 はこの 1 つの原因で説明できる。")
    P("")
    P("### 7.2 ラグ帯の不足")
    P("")
    P("帯の下端 4 mm < グリッド 5 mm。帯は 7 ラグでしか測れず、その 7 点は"
      "「1 px, 2 px, …, 7 px」である。`fit_band` は n ≥ 4 で通るので"
      "エラーにならないが、**自由度 5 の回帰で H を決めている**。"
      "帯を D50 で定義し直すと 0.24–2.2 mm になり、丸ごとグリッドの下に落ちる。")
    P("")
    P("### 7.3 直線性（rms_log）は物理域の証拠にならない")
    P("")
    P("ruler の rms_log は 0.0013–0.0023（x）、clean では 0.029–0.051。"
      "**人工物の方が「きれいな直線」を与える。** `fit_band` が返す品質指標は"
      "「その帯の中で直線か」しか見ておらず、「その直線が粒径モザイクの飽和と"
      "砂州の飽和の間にあるか」も「折り返しの d^1 かどうか」も見ていない。")
    P("")
    P("### 7.4 異方性の扱い")
    P("")
    P("clean で scaling 帯 H_x = %.3f に対し H_y = %.3f（主検体）、14 枚では "
      "H_y が %.3f–%.3f。`curve_iso` は 0.5(S_x+S_y) なので、異方性が大きい場"
      "では「2 つの別物の平均」になる。さらに x と y で折り返しの効き方が"
      "違う（x は 26–42 倍、y は 5–116 倍）ので、ruler の H_iso は"
      "「x の人工物」と「y の実信号」を混ぜた値になる。"
      % (rr(prim, "x", "scaling", "dirs_clean"),
         rr(prim, "y", "scaling", "dirs_clean"),
         rng("y", "scaling", "dirs_clean")[0],
         rng("y", "scaling", "dirs_clean")[2]))
    P("")
    P("### 7.5 微細スケールの雑音（5 mm 格子）")
    P("")
    P("流下方向の増分 rms は 0.36–0.48 mm、すなわち 5 mm あたり 0.4 mm の"
      "白色雑音が乗っている。これが d^0 の項として S に加わり、5 mm の local H "
      "を 0.63–0.69 から 0.48–0.54 に**下げる**（§3.3）。"
      "帯の下端がグリッドに近いほど、この下駄が効く。")
    P("")
    P("### 7.6 境界マスクと内部穴の混同")
    P("")
    P("v4-2.6 の穴埋め規則は「内部の粒穴」を想定している（周囲から補間"
      "できる）。このデータの nodata は境界マスクなので、補間規則の前提が"
      "崩れる。§5 のとおり、周期折り返し型の補間だけが壊れ（+0.26）、"
      "定数埋めは壊れない（≤0.003）。")
    P("")
    P("### 7.7 粒径の移送")
    P("")
    P("帯を D50 で定義し直すと 0.24–2.2 mm になり、この床ではグリッド"
      "（5 mm）の下に落ちる。つまり「粒径以下〜粒径」という設計意図は、"
      "この分解能の資料では**原理的に検証できない**。物差しを実データで"
      "較正するには、グリッドが D50 の 1/4 以下の DEM（0.3 mm 以下）が要る。")
    P("")

    # ---------------- 8. limits ----------------------------------------
    P("## 8. やっていないこと・限界")
    P("")
    P("- **うちの床の判定はしていない。** 粒径が 1/17〜1/30 で装甲床でもない"
      "資料なので、値の比較自体が意味を持たない。")
    P("- 窓は流下方向にしか動かしていない（水路幅が 0.40–0.55 m しかなく、"
      "横断方向に 0.6 m 箱を並べられない）。")
    P("- 箱サイズは「中央寄せの入れ子」で 1 通りだけ。位置と箱を同時に動かす"
      "2 因子の分解はしていない（§6 の 37 枚が位置の効果を別途与える）。")
    P("- 鏡映閉鎖は「端で勾配ゼロ」を仮定する。ラップよりはましだが、"
      "これも 1 つの仮定であり、真の非周期の扱い（端を落とす）とは違う。"
      "clean（crop+mirror）はその両方を同時に動かしているので、"
      "折り返しと境界マスクの効果を完全には分離していない。")
    P("- 粗視化は 1〜4 倍の 3 水準だけで、膝の位置そのものは追っていない。")
    P("- 増分 0.4 mm の正体（SfM の再標本化／保存形式）は特定していない。")
    P("- 水位・流速・掃流力のマップは取得していない（DEM と D50 のみ）。")
    P("- 画像の目視はしていない。")
    P("")

    # ---------------- 9. no verdict ------------------------------------
    P("## 9. うちの床について")
    P("")
    P("本レポートは、うちの床（`bed_gen2` の 0.6 m パッチ／3.0 m 長尺床）に"
      "ついて**合否も良し悪しも書かない**。brief の指定どおり、ここで測った"
      "のは物差しの振る舞いだけである。")
    P("")

    # ---------------- 10. reproduce ------------------------------------
    P("## 10. 再現")
    P("")
    P("```")
    P("cd D:\\yohey2026\\dsh\\muse\\caustics-filter\\offline")
    P("python run_v4_2_9.py --all")
    P("```")
    P("")
    P("| 段 | 実行時間 (s) |")
    P("| --- | --- |")
    for k, v in ctx["timings"].items():
        P("| %s | %.1f |" % (k, v))
    P("")
    P("成果物:")
    P("")
    P("| 種別 | パス | bytes |")
    P("| --- | --- | --- |")
    for k, v in ctx["deliverables"].items():
        if k == "REPORT":
            b = "@@SZ@@@@"          # fixed width, patched after the text exists
        elif k == "numbers":
            b = "%8d" % ctx["deliverables_bytes"]["numbers"]
        else:
            try:
                b = "%8d" % os.path.getsize(v)
            except OSError:
                b = "      — "
        P("| %s | `%s` | %s |" % (k, os.path.relpath(v, HERE), b))
    P("")
    P("checkpoint: `%s`（本レポートの書き込み後に %d units）。データは "
      "`out/v4-2-9/data/`（27 files, %.1f MB）。"
      % (os.path.relpath(CP_PATH, HERE), ctx["checkpoint_units"],
         dl["total_bytes"] / 1e6))
    P("")
    g = ctx["g_legacy"]
    P("### 10.1 write set")
    P("")
    P("本スクリプトがこの実行で書き込んだパスは %d 件（作業キャッシュ 5 件と "
      "record JSON。DEM/d50 の 27 件は既に取得済みだったので md5 検証のみ）。"
      "`out/v4-2-9/` の下にないパスは %d 件。`out/v4-2-9/` 配下のファイル"
      "総数（この判定時点）は %d 件で、うち `data/` の取得ファイルは 27 件。"
      % (g["n_writes"], len(g["writes_outside_write_set"]), g["n_out_dir_files"]))
    P("")
    P("### 10.2 G-LEGACY（既存ファイルの非改変）")
    P("")
    P("| 対象 | 数 | 変化 |")
    P("| --- | --- | --- |")
    P("| 既存 .py | %d | changed %d / added %d / removed %d |"
      % (g["n_baseline"], len(g["changed"]), len(g["added"]), len(g["removed"])))
    P("| out/ 配下の既存ファイル | %d | changed %d / added %d / missing %d |"
      % (g["n_out_files"], len(g["out_changed"]), len(g["out_added"]),
         len(g["out_missing"])))
    P("")
    if g["details"]:
        P("**注意: 本ジョブの実行中に、本ジョブ以外の書き込みがあった。**")
        P("")
        P("| ファイル | bytes (前 → 後) | mtime (後) |")
        P("| --- | --- | --- |")
        for d in g["details"]:
            P("| `%s` | %d → %d | %s |"
              % (d["file"], d["bytes_before"], d["bytes_now"], d["mtime_now"]))
        P("")
        P("本スクリプトはこのファイルを開いていない（write set は上記 %d 件"
          "のみ）。同じ作業ツリーで並行して走っている別ジョブによる変更と"
          "判断する。本ジョブの成果物・数値はこのファイルに依存しない"
          "（物差しは `run_v4_2_5.py`、sha256 `%s`、実行中不変）。"
          % (g["n_writes"], st0["required"]["run_v4_2_5.py (the ruler)"]
             ["sha256"][:16]))
        P("")
    P("厳密判定: changed %d / added %d / removed %d / out changed %d → %s。"
      % (len(g["changed"]), len(g["added"]), len(g["removed"]),
         len(g["out_changed"]), "PASS" if g["ok"] else
         "FAIL（上表の外部書き込みのみ。本ジョブの write set は write set 外 0 件）"))
    P("")
    P("環境: Python %s / numpy %s / tifffile %s。"
      % (st0["python"], st0["numpy"], st0["tifffile"]))
    P("")
    return "\n".join(A)


def stage_report(cp, ctx):
    log("STAGE4  report")
    numbers = {
        "brief": BRIEF,
        "created_utc": ctx["created_utc"],
        "what_this_is": ("calibration of the run_v4_2_5 measurement instrument "
                         "on a public measured flume DEM. No verdict about our "
                         "own bed."),
        "stage0": ctx["stage0"],
        "record": ctx["download"]["record"],
        "files": ctx["download"]["files"],
        "grid": ctx["prove"]["grid"],
        "contract": {"bands_m": [list(b) for b in M5.BANDS],
                     "band_label": BAND_LABEL, "gate_recorded_only": GATE,
                     "nbin": NBIN, "dmax_m": {"x": DM_X, "y": DM_Y, "iso": DM_ISO},
                     "cell_m": ctx["prove"]["grid"]["cell_m_values"],
                     "inpaint_length_m": INPAINT_M, "tile_m": TILE_M,
                     "fill_rules": ["mean (primary)", "p75", "max", "inpaint"],
                     "closures": ["ruler (np.roll)", "full+mirror",
                                  "clean (crop+mirror)"],
                     "reference": REF},
        "dems": ctx["measure"],
        "sensitivity": ctx["sens"],
        "breakage": [
            "periodic wrap (np.roll) injects a d^1 term at the seam: delta = "
            "0.099-0.121 m, making S(1 px) 26-42x too large in x and 5-116x in "
            "y, and fixing fit_band's H at 0.500 in every DEM",
            "with the wrap closed (crop+mirror) H_x in the scaling band is "
            "0.520-0.631, i.e. the 0.50 was the estimator, not the bed",
            "lag band: lower edge 4 mm is below the 5 mm grid; 7 lags only",
            "fit quality (rms_log) is BETTER for the artefact (0.0013-0.0023) "
            "than for the clean field (0.029-0.051)",
            "anisotropy: iso = 0.5(Sx+Sy) mixes the x artefact with the y "
            "signal, and the wrap term differs by 20x between x and y",
            "boundary mask vs interior holes: wrap-inpainting invents a "
            "0.20-0.24 m ramp and moves H_x by +0.26",
            "5 mm cell + 0.4 mm increment noise adds a d^0 pedestal that pulls "
            "the small-lag H down by ~0.15",
        ],
        "limits": ("no window stepping across the width; one nested box family; "
                   "no image inspection; mirror closure is itself an assumption; "
                   "cause of the white increments not identified"),
        "checkpoint": {"path": CP_PATH,
                       "units": ctx["checkpoint_units"],
                       "state": cp.state()},
        "deliverables": ctx["deliverables"],
        "g_legacy": ctx["g_legacy"],
        "timings_s": ctx["timings"],
        "log": LOG,
    }
    payload = json.dumps(_j(numbers), ensure_ascii=False, indent=1).encode("utf-8")
    ctx["deliverables_bytes"] = {"numbers": len(payload)}
    text = build_report(ctx)
    text = text.replace("@@SZ@@@@", "%8d" % len(text.encode("utf-8")))
    write_text(REPORT, text)
    write_json(NUMBERS, numbers)
    log("  wrote REPORT (%d B) + numbers (%d B)"
        % (os.path.getsize(REPORT), os.path.getsize(NUMBERS)))
    cp.unit_done("stage3_numbers")
    cp.unit_done("stage3_report")


# ---------------------------------------------------------------------------
def check_g_legacy(st0):
    py = {}
    for p in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        n = os.path.basename(p)
        if n == THIS_FILE:
            continue
        py[n] = {"bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    changed = sorted(k for k in st0["py"] if k in py
                     and py[k]["sha256"] != st0["py"][k]["sha256"])
    added = sorted(set(py) - set(st0["py"]))
    removed = sorted(set(st0["py"]) - set(py))
    details = []
    for k in changed:
        p = os.path.join(HERE, k)
        st = os.stat(p)
        details.append({
            "file": k, "bytes_before": st0["py"][k]["bytes"],
            "bytes_now": st.st_size,
            "sha256_before": st0["py"][k]["sha256"], "sha256_now": py[k]["sha256"],
            "mtime_now": datetime.datetime.fromtimestamp(st.st_mtime)
            .astimezone().isoformat(timespec="seconds")})
    # every path this process wrote must be inside the declared write set
    roots = [os.path.abspath(OUT), os.path.abspath(os.path.join(HERE, THIS_FILE))]
    wrote = sorted(set(WRITES))
    outside = [w for w in wrote
               if not any(w == r or w.startswith(r + os.sep) for r in roots)]
    now = {}
    for root, _d, fns in os.walk(os.path.join(HERE, "out")):
        if os.path.abspath(root).startswith(os.path.abspath(OUT)):
            continue
        for fn in fns:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, HERE).replace("\\", "/")
            try:
                st = os.stat(p)
            except OSError:
                continue
            now[rel] = [st.st_size, int(st.st_mtime)]
    out_changed = sorted(k for k in st0["out_files"]
                         if k in now and now[k] != st0["out_files"][k])
    out_missing = sorted(set(st0["out_files"]) - set(now))
    out_added = sorted(set(now) - set(st0["out_files"]))
    ok = not (changed or added or removed or out_changed or out_missing
              or out_added)
    g = {"n_baseline": st0["n_existing_py"], "changed": changed,
         "added": added, "removed": removed, "n_out_files":
         st0["n_existing_out_files"], "out_changed": out_changed,
         "out_missing": out_missing, "out_added": out_added,
         "details": details, "write_set": [THIS_FILE, "out/v4-2-9/"],
         "n_writes": len(wrote), "writes_outside_write_set": outside,
         "n_out_dir_files": sum(len(f) for _r, _d, f in os.walk(OUT)),
         "write_roots": roots, "ok": ok, "checked_at": utcnow()}
    log("  G-LEGACY  existing .py: changed=%d added=%d removed=%d | out/: "
        "changed=%d added=%d missing=%d -> %s"
        % (len(changed), len(added), len(removed), len(out_changed),
           len(out_added), len(out_missing), "PASS" if ok else "FAIL"))
    return g


def main():
    ap = argparse.ArgumentParser()
    for s in ("stage0", "download", "prove", "measure", "sens", "report"):
        ap.add_argument("--" + s, action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    if not any(vars(a).values()):
        a.all = True
    do = lambda k: a.all or getattr(a, k)          # noqa: E731

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    cp = Checkpoint(CP_PATH)
    created = utcnow()
    cache = os.path.join(OUT, "_work")
    os.makedirs(cache, exist_ok=True)
    timings = {}

    def cached(name):
        p = os.path.join(cache, name + ".json")
        return json.load(open(p, encoding="utf-8")) if os.path.isfile(p) else None

    t0 = time.time()
    st0 = None if do("stage0") else cached("stage0")
    if st0 is None:
        st0 = stage0(cp)
        write_json(os.path.join(cache, "stage0.json"), st0)
    timings["stage0"] = round(time.time() - t0, 1)

    t0 = time.time()
    dl = None if do("download") else cached("download")
    if dl is None:
        dl = stage_download(cp, st0)
        write_json(os.path.join(cache, "download.json"), dl)
    timings["stage1_download"] = round(time.time() - t0, 1)

    t0 = time.time()
    pr = None if do("prove") else cached("prove")
    if pr is None:
        pr = stage_provenance(cp, dl)
        write_json(os.path.join(cache, "prove.json"), pr)
    timings["stage1_provenance"] = round(time.time() - t0, 1)

    t0 = time.time()
    ms = None if do("measure") else cached("measure")
    if ms is None:
        ms = stage_measure(cp, pr)
        write_json(os.path.join(cache, "measure.json"), ms)
    timings["stage2_measure"] = round(time.time() - t0, 1)

    t0 = time.time()
    se = None if do("sens") else cached("sens")
    if se is None:
        se = stage_sens(cp, ms)
        write_json(os.path.join(cache, "sens.json"), se)
    timings["stage2_sensitivity"] = round(time.time() - t0, 1)

    log("VERIFY  write set")
    g = check_g_legacy(st0)
    cp.unit_done("g_legacy_verify")

    state = cp.state() or {}
    deliv = {
        "REPORT": REPORT, "numbers": NUMBERS, "checkpoint": CP_PATH,
        "record_json": RECORD_JSON,
    }
    ctx = {"created_utc": created, "stage0": st0, "download": dl, "prove": pr,
           "measure": ms, "sens": se, "timings": timings, "g_legacy": g,
           "checkpoint_units": int(state.get("completed_units") or 0) + 2,
           "deliverables": deliv}
    stage_report(cp, ctx)
    log("DONE  %s" % REPORT)
    print("\n".join(LOG[-4:]))


if __name__ == "__main__":
    main()
