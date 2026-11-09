"""Environment resolution (port of caustics-filter/web/env.js resolveEnvironment).

Input: recipe dict (caustics-recipe/v1). Output: dict of derived physical
parameters used by the offline pipeline. Formulas mirror web/env.js so the
offline path stays consistent with the web front-end definition.
"""
import math


def _clamp(x, a, b):
    return max(a, min(b, x))


def _coupling_weight(couplings, cid):
    for c in couplings or []:
        if c.get("id") == cid:
            if not c.get("enabled", False):
                return 0.0
            return _clamp(float(c.get("weight", 1.0)), 0.0, 2.0)
    return 0.0


def solar_position(date_str, time_str, lat_deg, lon_deg):
    dp = str(date_str or "").split("-")
    tp = str(time_str or "").split(":")
    try:
        Y = int(dp[0])
    except Exception:
        Y = 2000
    try:
        M = int(dp[1])
    except Exception:
        M = 1
    try:
        D = int(dp[2])
    except Exception:
        D = 1
    try:
        hh = int(tp[0])
    except Exception:
        hh = 0
    try:
        mm = int(tp[1])
    except Exception:
        mm = 0
    zone_h = round((lon_deg or 0.0) / 15.0)
    utc_h = hh - zone_h + mm / 60.0
    a = math.floor((14 - M) / 12)
    y = Y + 4800 - a
    m = M + 12 * a - 3
    jdn = (D + math.floor((153 * m + 2) / 5) + 365 * y
           + math.floor(y / 4) - math.floor(y / 100)
           + math.floor(y / 400) - 32045)
    jd = jdn + (utc_h - 12.0) / 24.0
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360.0
    g = math.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = math.radians(L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g))
    eps = math.radians(23.439 - 0.0000004 * n)
    alpha = math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam))
    delta = math.asin(_clamp(math.sin(eps) * math.sin(lam), -1.0, 1.0))
    gmst = (18.697374558 + 24.06570982441908 * n) % 24.0
    lst = (gmst + (lon_deg or 0.0) / 15.0) % 24.0
    H = math.radians(lst * 15.0) - alpha
    lat = math.radians(lat_deg or 0.0)
    sin_alt = (math.sin(lat) * math.sin(delta)
               + math.cos(lat) * math.cos(delta) * math.cos(H))
    alt = math.asin(_clamp(sin_alt, -1.0, 1.0))
    az = math.atan2(math.sin(H),
                    math.cos(H) * math.sin(lat) - math.tan(delta) * math.cos(lat))
    az_deg = (math.degrees(az) + 180.0 + 360.0) % 360.0
    return {"altitudeDeg": math.degrees(alt), "azimuthDeg": az_deg}


def sun_color_from_altitude(alt_deg):
    day_amt = _clamp((alt_deg + 2.0) / 37.0, 0.0, 1.0)
    day_amt = day_amt * day_amt * (3.0 - 2.0 * day_amt)
    low = (1.0, 0.42, 0.13)
    high = (1.0, 0.97, 0.92)
    c = [low[i] + (high[i] - low[i]) * day_amt for i in range(3)]
    if alt_deg < -2.0:
        k = _clamp((alt_deg + 12.0) / 10.0, 0.0, 1.0)
        c = [0.10 * k + 0.02, 0.14 * k + 0.02, 0.22 * k + 0.03]
    return c


def resolve_environment(state):
    couplings = state.get("couplings") or []
    cloud = _clamp(float(state["weather"].get("cloudPct", 20)) / 100.0, 0.0, 1.0)
    clarity = state["weather"].get("clarity", 0.8)
    clarity = _clamp(float(clarity), 0.0, 1.0)
    wind = max(0.0, float(state["windWave"].get("windSpeedMS", 0.0)))
    wind_dir = math.radians(float(state["windWave"].get("windDirDeg", 0.0)))
    fineness = _clamp(float(state["windWave"].get("fineness", 0.5)), 0.0, 1.0)
    flow_spd = _clamp(float(state["windWave"].get("flowSpeed", 0.5)), 0.0, 1.0)
    depth = max(0.1, float(state["water"].get("depthM", 1.0)))
    transp = _clamp(float(state["water"].get("transparency", 0.7)), 0.0, 1.0)
    turbid = _clamp(float(state["water"].get("turbidity", 0.3)), 0.0, 1.0)
    pigment = _clamp(float(state["water"].get("pigment", 0.3)), 0.0, 1.0)

    sun = solar_position(state["sun"].get("date"), state["sun"].get("time"),
                         float(state["sun"].get("latitudeDeg", 35.7)),
                         float(state["sun"].get("longitudeDeg", 139.7)))
    sun_color = sun_color_from_altitude(sun["altitudeDeg"])
    up = _clamp(math.sin(math.radians(max(sun["altitudeDeg"], 0.0))), 0.0, 1.0)
    night_k = 1.0
    if sun["altitudeDeg"] < -2.0:
        night_k = _clamp((sun["altitudeDeg"] + 12.0) / 10.0, 0.05, 1.0)

    cloud_w = _coupling_weight(couplings, "cloud-sun")
    sun_i = (1.0 - cloud * 0.88 * min(cloud_w, 1.5)) * (0.35 + 0.65 * clarity)
    sun_i *= (0.15 + 0.85 * up)
    sun_i *= night_k
    sun_i = _clamp(sun_i, 0.0, 1.6)

    wind_w = _coupling_weight(couplings, "wind-wave")
    wave_amp = ((0.12 + pow(min(wind, 20.0) / 20.0, 0.8) * 1.5)
                * (0.25 + 0.75 * min(wind_w, 1.5)))
    if wind_w == 0.0:
        wave_amp = 0.35

    Lphys = (wind * wind) / 9.8
    spec_L = _clamp(0.30 + Lphys * 0.30, 0.30, 4.0)
    spec_amp = ((0.10 + pow(min(wind, 20.0) / 20.0, 0.8) * 1.30)
                * (0.25 + 0.75 * min(wind_w, 1.5)))
    if wind_w == 0.0:
        spec_amp = 0.30

    sc_w = _coupling_weight(couplings, "sun-caustics")
    caustics_i = sun_i * (0.35 + 0.65 * up) * 0.9 * min(max(sc_w, 0.0), 2.0)
    tb_w = _coupling_weight(couplings, "turbid-caustics")
    caustics_i *= (1.0 - turbid * 0.75 * min(tb_w, 1.5))
    caustics_i = _clamp(caustics_i, 0.0, 2.5)

    gl_w = _coupling_weight(couplings, "cloud-glitter")
    glitter_i = _clamp((0.4 + sun_i * 1.6) * night_k, 0.0, 2.5)
    diffuse = _clamp(cloud * min(max(gl_w, 0.0), 1.5), 0.0, 1.0)
    if gl_w == 0.0:
        diffuse = cloud * 0.5

    dt_w = _coupling_weight(couplings, "depth-tint")
    depth_shade = _clamp(math.log10(1.0 + depth) / math.log10(31.0), 0.0, 1.0)
    depth_shade *= (0.25 + 0.75 * (1.0 - transp * 0.6)) * min(max(dt_w, 0.0), 2.0)
    if dt_w == 0.0:
        depth_shade = 0.25
    depth_shade = _clamp(depth_shade, 0.0, 1.0)

    absorb_scale = (0.55 + 0.45 * (1.0 - transp)) * (1.0 + turbid * 1.2)
    absorb = [0.34 * absorb_scale + pigment * 0.10,
              0.057 * absorb_scale + pigment * 0.12,
              0.0092 * absorb_scale + pigment * 0.03 + (1.0 - transp) * 0.02]
    scatter = [0.03 + pigment * 0.10,
               0.06 + pigment * 0.12,
               0.05 + pigment * 0.03]

    is_preset0 = (state.get("name") == "真夏の沖縄・正午")
    caustics_kc = 1.0 if is_preset0 else 0.8
    dg_cap = 3.40 if is_preset0 else 4.20

    refract = 0.028 * (0.5 + wave_amp * 0.7)
    chroma = 0.10 + (1.0 - clarity) * 0.10

    flow_speed = 0.05 + flow_spd * 0.5 + wind * 0.008
    raw_fx = math.sin(wind_dir) * flow_speed * 0.35
    raw_fy = -math.cos(wind_dir) * flow_speed * 0.35
    flow_per_sec = [round(raw_fx * 8.0) / 8.0, round(raw_fy * 8.0) / 8.0]

    alt_r = math.radians(sun["altitudeDeg"])
    az_r = math.radians(sun["azimuthDeg"])
    sun_dir = [math.sin(az_r) * math.cos(alt_r),
               math.cos(az_r) * math.cos(alt_r),
               math.sin(alt_r)]
    if sun_dir[2] < 0.02:
        sun_dir[2] = 0.02
    sl = math.sqrt(sum(c * c for c in sun_dir))
    sun_dir = [c / sl for c in sun_dir]

    return {
        "sun": sun, "sunDir": sun_dir, "sunColor": sun_color,
        "sunI": sun_i, "diffuse": diffuse,
        "waveAmp": wave_amp,
        "specL": spec_L, "specAmp": spec_amp,
        "depthM": depth, "fine": fineness,
        "flowPerSec": flow_per_sec,
        "windDir": [math.sin(wind_dir), -math.cos(wind_dir)],
        "refract": refract, "chroma": chroma,
        "causticsI": caustics_i, "causticsKc": caustics_kc,
        "dgCap": dg_cap, "focusM": 0.8,
        "glitterI": glitter_i, "absorb": absorb, "scatter": scatter,
        "depthShade": depth_shade, "turbid": turbid,
        "clarity": clarity, "cloud": cloud, "pigment": pigment,
        "transparency": transp, "windSpeedMS": wind,
        "waveHeightM": wave_amp * 0.12,
        "strength": float(state.get("render", {}).get("strength", 1.0)),
        "exposure": float(state.get("render", {}).get("exposure", 1.0)),
    }
