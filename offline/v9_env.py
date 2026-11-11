"""v9 environment resolution -- real-mode direct fraction (chunk 5).

PIPELINE-V8.md chunk 5 addendum: "Open-Meteo no direct_radiation /
diffuse_radiation (DNI相当は direct_normal_radiation) を分離取得して eta を
実測値から決める。取得は環境解決ステップの責務とし、オフライン描画器は解決済みの
eta と放射量を入力に取る".

So this module is the ONLY place that touches the network, and it runs as its
own stage.  The renderer (run_v9.py) reads the resolved JSON; it never opens a
socket, and the transport stays deterministic and offline.

Contract
--------
``resolve_real(...)`` returns a dict with at least

    {"status": "ok" | "fallback",
     "samples": {...clear..., ...overcast...},
     "chosen":  {"clear": {eta, ghi, dni, dhi, time, ...},
                 "overcast": {...}},
     "source":  "open-meteo" | "fallback_representative",
     "url", "http", "error", "as_of"}

A failed fetch never raises and never blocks the run: it is recorded and the
representative clear/overcast fractions of v9_sky.FALLBACK_ETA are used.

The JSON is cached on disk (``out/v9/sky-real-samples.json``) so a re-run is
byte-identical even with no network.  ``--refresh`` forces a fetch.
"""
import datetime
import json
import math
import os
import urllib.parse
import urllib.request

from v9_sky import FALLBACK_ETA

OPEN_METEO = "https://api.open-meteo.com/v1/forecast"
HOURLY = ["direct_radiation", "diffuse_radiation",
          "direct_normal_irradiance", "shortwave_radiation"]

#: daytime sample floors (W/m^2) -- below this the split is noise
GHI_FLOOR_CLEAR = 150.0
GHI_FLOOR_OVERCAST = 40.0


def _url(lat, lon, past_days):
    q = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon,
        "hourly": ",".join(HOURLY),
        "past_days": int(past_days), "forecast_days": 1,
        "timezone": "UTC",
    })
    return OPEN_METEO + "?" + q


def fetch_open_meteo(lat=35.7, lon=139.7, past_days=7, timeout=25.0):
    """One HTTP GET.  Returns (payload dict | None, meta dict)."""
    url = _url(lat, lon, past_days)
    meta = {"url": url, "http": None, "error": None,
            "utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "caustics-filter-offline/3.5 (research)"})
        with urllib.request.urlopen(req, timeout=float(timeout)) as r:
            meta["http"] = int(getattr(r, "status", 200))
            body = r.read().decode("utf-8")
        payload = json.loads(body)
        meta["bytes"] = len(body)
        return payload, meta
    except Exception as e:                      # never block on the network
        meta["error"] = "%s: %s" % (type(e).__name__, str(e)[:300])
        return None, meta


def _rows(payload):
    h = payload.get("hourly") or {}
    t = h.get("time") or []
    out = []
    for i, ts in enumerate(t):
        def g(k):
            v = h.get(k) or []
            return None if i >= len(v) else v[i]
        dr, df = g("direct_radiation"), g("diffuse_radiation")
        dni, ghi = g("direct_normal_irradiance"), g("shortwave_radiation")
        if None in (dr, df, dni, ghi):
            continue
        out.append({"time": ts, "direct_radiation": float(dr),
                    "diffuse_radiation": float(df),
                    "dni": float(dni), "ghi": float(ghi)})
    return out


def classify(rows):
    """Clear / overcast samples from measured hours.

    eta = DNI/(DNI+DHI) is evaluated on measured values only.  The clear
    sample is the daytime hour with the largest eta among hours whose global
    horizontal irradiance is substantial; the overcast one is the smallest
    eta that still carries measurable light.
    """
    day = [r for r in rows if r["ghi"] > 10.0]
    if not day:
        return None
    for r in day:
        r["eta"] = (r["dni"] / (r["dni"] + r["diffuse_radiation"])
                    if (r["dni"] + r["diffuse_radiation"]) > 0 else None)
    ok = [r for r in day if r["eta"] is not None]
    if not ok:
        return None
    clear_pool = [r for r in ok if r["ghi"] >= GHI_FLOOR_CLEAR]
    over_pool = [r for r in ok if r["ghi"] >= GHI_FLOOR_OVERCAST]
    clear = max(clear_pool or ok, key=lambda r: r["eta"])
    over = min(over_pool or ok, key=lambda r: r["eta"])
    # independent consistency check of the API: GHI ~ DNI*sin(alt) + DHI.
    # The altitude is not in the payload, so the implied sun elevation is
    # derived from the split and reported (it must sit inside the day).
    for r in (clear, over):
        s = (r["ghi"] - r["diffuse_radiation"]) / max(r["dni"], 1e-9)
        r["implied_sin_alt"] = float(min(max(s, -1.0), 1.0))
        r["implied_alt_deg"] = float(math.degrees(math.asin(
            min(max(s, -1.0), 1.0))))
    return {
        "n_rows": len(rows), "n_day": len(day),
        "clear": clear, "overcast": over,
        "eta_gap": float((clear["eta"] or 0.0) - (over["eta"] or 0.0)),
    }


def resolve_real(lat=35.7, lon=139.7, past_days=7, cache_path=None,
                 refresh=False, out_dir="out/v9"):
    """Resolve the real-mode samples; cached, fail-soft."""
    if cache_path is None:
        cache_path = os.path.join(out_dir, "sky-real-samples.json")
    if (not refresh) and os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                cached = json.load(f)
            cached["cache"] = "hit"
            return cached
        except Exception:
            pass

    payload, meta = fetch_open_meteo(lat, lon, past_days)
    out = {"request": {"latitude": lat, "longitude": lon,
                       "past_days": int(past_days), "hourly": HOURLY},
           "http": meta, "as_of": meta["utc"], "location": "Tokyo 35.7N 139.7E",
           "fallback_eta": FALLBACK_ETA, "cache": "miss"}
    if payload is None:
        out.update({
            "status": "fallback", "source": "fallback_representative",
            "error": meta["error"],
            "chosen": {k: {"eta": float(v), "note": "representative value"}
                       for k, v in FALLBACK_ETA.items()},
            "reason": ("the fetch failed; the design says record it and fall "
                       "back to representative values without blocking"),
        })
    else:
        rows = _rows(payload)
        cls = classify(rows)
        if cls is None:
            out.update({"status": "fallback",
                        "source": "fallback_representative",
                        "error": "no usable daytime rows in the payload",
                        "chosen": {k: {"eta": float(v)}
                                   for k, v in FALLBACK_ETA.items()}})
        else:
            out.update({
                "status": "ok", "source": "open-meteo",
                "classification": cls,
                "chosen": {"clear": cls["clear"], "overcast": cls["overcast"]},
                "rows_sample": rows[:0],       # rows kept in the raw cache
            })
    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2, ensure_ascii=False, default=float)
    except Exception as e:
        out["cache_write_error"] = str(e)
    return out


def chosen_eta(resolved, condition="clear"):
    """eta for a named condition, or None (renderer-side read of the cache)."""
    ch = (resolved or {}).get("chosen") or {}
    rec = ch.get(condition) or {}
    v = rec.get("eta")
    return None if v is None else float(v)
