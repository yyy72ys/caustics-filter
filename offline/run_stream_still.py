"""One composited still: v3.4 (v8) pipeline, eta = 1.0, gentle stream edge.

Task (D:\\mutaro-work\\rally-in\\BRIEF-STREAM-STILL.md):
  base frame  blender-piece8s/renders_v2/frame_0000.png   (832x480)
  caustics    v3.4 / v8 chunk-4 state (V8WaterState) + det_frame, master 2048^2
  composite   frozen production: A2.6-Blow + E2.2, lowpass 3.0 px, Q = 25
  light       eta = 1.0 -> direct sun only; the v9/Perez diffuse route is NOT
              used (a recipe with no sky block is eta = 1 exactly)
  environment gentle stream edge, morning (same sun/sky family as preset 1):
              depth 0.30 m, wind 2.5 m/s along the corridor, fetch 1000 m,
              rain 0, surface film none.
              H_s comes from the pipeline's own SMB/JONSWAP fetch computation.

Every step is checkpointed with D:\\mutaro-work\\jobkit (Checkpoint.unit_done)
and the three deliverables are committed atomically (atomic_write_bytes).

Deliverables -> out/stream-still/
  comp-stream-still.png   832x480 composited still
  numbers-stream.json     recipe + realised values + gates
  REPORT-STREAM.md        numbers table, pond-vs-stream (realised), gaps

Run (foreground, from caustics-filter/offline):
  python run_stream_still.py
"""
import gc
import hashlib
import io
import json
import math
import os
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\mutaro-work\jobkit")

from checkpoint import Checkpoint, atomic_write_bytes, verify_units  # noqa: E402

from env import resolve_environment  # noqa: E402
from det import det_frame  # noqa: E402
from state import derive_seed  # noqa: E402
from image import srgb_to_linear  # noqa: E402
from v5_physics import (water_attenuation, beer_lambert, smith_g1_beckmann,  # noqa: E402
                        scattering_blur_sigma_px, substrate_gain)
from v5_composite import (A_V3, B_SIGMA_V3, E_V3, render_frame_v5,  # noqa: E402
                          image_metrics)
from recipe_v8 import parse_recipe_v8, with_v8  # noqa: E402
from v8_state import V8WaterState  # noqa: E402
import v8_physics as P  # noqa: E402

RECIPE = os.path.join(HERE, "..", "blender", "recipe-1.json")
PHOTO = os.path.join(HERE, "..", "blender-piece8s", "renders_v2",
                     "frame_0000.png")
OUT = os.path.join(HERE, "out", "stream-still")
WORK = os.path.join(OUT, "_work")
CP = Checkpoint(os.path.join(OUT, "checkpoint.json"))

GRID = 2048
Q = 25
ETA = 1.0                 # direct only (the v8 route; no Perez diffuse)
A_COMP = A_V3             # 2.6  (frozen production ambient scale)
E_COMP = E_V3             # 2.2  (frozen production C-layer exposure)
B_LOWPASS = B_SIGMA_V3    # 3.0 px

# --- the brief's recipe numbers (stream edge) ---
STREAM = {"fetch_m": 1000.0, "wind_ms": 2.5, "depth_m": 0.30}
# --- preset 1 pond, for the side-by-side (same code, same seed) ---
POND = {"fetch_m": 100.0, "wind_ms": 1.5, "depth_m": 0.60}

CMEAN_TOL = 1e-3
CRUSH_BOUND = 0.02
BLOW_BOUND = 0.0


def log(msg):
    print(msg, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def coerce(obj):
    """numpy -> plain JSON types."""
    if isinstance(obj, dict):
        return {str(k): coerce(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [coerce(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return coerce(obj.tolist())
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


def c_stats(C):
    c = np.asarray(C, dtype=np.float64)
    return {
        "shape": list(c.shape),
        "mean": float(c.mean()),
        "p50": float(np.percentile(c, 50)),
        "p99": float(np.percentile(c, 99)),
        "p99_9": float(np.percentile(c, 99.9)),
        "max": float(c.max()),
        "min": float(c.min()),
        "mean_abs_dev_from_1": float(abs(c.mean() - 1.0)),
        "pass_mean_1": bool(abs(c.mean() - 1.0) <= CMEAN_TOL),
    }


def scene_for(raw, par, tag):
    return with_v8(raw, fetch_m=float(par["fetch_m"]),
                   depth_m=float(par["depth_m"]),
                   wind_ms=float(par["wind_ms"]))


def build(raw, scene, tag, grid, seed, photo_lin, PH, PW, C=None):
    """State -> (C master) -> still.  `C` lets the caller reuse a master."""
    t0 = time.time()
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    st = V8WaterState(scene, env, N=grid, seed=seed, preset_index=1,
                      fetch_m=float(cfg["fetch"]["fetch_m"]),
                      fetch_source=cfg["fetch"]["source"])
    fields = st.fields_at(0.0)
    fields["S"] = st.S
    re = st.realised_stats(fields)
    t_state = time.time() - t0
    log("  [%s] state %.1fs sigma2_spec=%.6f realised=%.6f Hs=%.4f"
        % (tag, t_state, st.spectrum_after["sigma2_rad2"],
           re["sigma2_rad2"], re["Hs_m"]))

    if C is None:
        t1 = time.time()
        C, info = det_frame(fields, st.S, grid, env, Q)
        t_det = time.time() - t1
        log("  [%s] det Q=%d %.1fs mean=%.9f p99=%.4f max=%.3f"
            % (tag, Q, t_det, float(C.mean()),
               float(np.percentile(C, 99)), float(C.max())))
    else:
        t_det = None
        info = {"reused": True, "Q": Q,
                "note": ("master reused from the earlier run of this same "
                         "script in this session (same seed, same Q=25); the "
                         "det stage was not re-timed")}

    att = water_attenuation(env["turbid"])
    beer = beer_lambert(env["depthM"], env["sun"]["altitudeDeg"], env["turbid"])
    blur = scattering_blur_sigma_px(att, env["depthM"], grid)
    gain_base = substrate_gain(cfg["bottom"]["substrate"],
                               cfg["bottom"]["albedoScale"])
    sun_zenith = 90.0 - float(env["sun"]["altitudeDeg"])
    g1 = smith_g1_beckmann(re["rms_slope"], sun_zenith)

    t2 = time.time()
    img, dbg = render_frame_v5(
        fields, C, photo_lin, env, PW, PH, ambient_scale=A_COMP,
        c_layer_exposure=E_COMP,
        substrate_gain_rgb=gain_base["gain_rgb"],
        beer_down=beer["down"], beer_up=beer["up"], g1=float(g1["G1"]),
        blur_sigma_px=blur["sigma_px"], lowpass_sigma_px=B_LOWPASS,
        return_debug=True)
    metrics = image_metrics(dbg["lin"], dbg["aces"])
    t_comp = time.time() - t2
    del dbg
    gc.collect()
    log("  [%s] composite %.1fs crush=%.5f blow=%.6f"
        % (tag, t_comp, metrics["crush_rate_linear_lt_0.02"],
           metrics["blowout_rate_aces_gt_0.98"]))

    return {
        "env": env, "cfg": cfg, "state": st, "fields": fields, "realised": re,
        "C": C, "C_stats": c_stats(C), "det_info": info,
        "still": img, "still_metrics": metrics,
        "beer": {"down": [float(x) for x in beer["down"]],
                 "up": [float(x) for x in beer["up"]],
                 "two_way": [float(x) for x in beer["two_way"]]},
        "blur_sigma_px": float(blur["sigma_px"]),
        "substrate_gain_rgb": list(gain_base["gain_rgb"]),
        "g1": float(g1["G1"]), "sun_zenith_deg": float(sun_zenith),
        "timing_s": {"state": t_state, "det": t_det, "composite": t_comp},
    }


def smb_node(env, fetch_m):
    rel = P.jonswap_fetch_relations(env["windSpeedMS"], fetch_m)
    sp = P.spreading_parameters(env["windSpeedMS"], fetch_m)
    return {
        "wind_ms": float(env["windSpeedMS"]), "fetch_m": float(fetch_m),
        "Ft": float(rel["Ft"]),
        "f_p_Hz": float(rel["f_p_Hz"]),
        "omega_p_rad_s": float(rel["omega_p_rad_s"]),
        "alpha": float(rel["alpha"]),
        "Hs_fetch_law_m": float(rel["Hs_fetch_law_m"]),
        "lambda_p_m": float(P.G / (2.0 * math.pi * rel["f_p_Hz"] ** 2)),
        "s_max_mitsuyasu": float(sp["s_max"]),
        "U_over_cp": float(sp["U_over_cp"]),
        "c_p_m_per_s": float(sp["c_p_m_per_s"]),
        "in_fitted_range_Ft_1e3_1e5": bool(1e3 <= rel["Ft"] <= 1e5),
        "relations": rel["relations"],
    }


def stream_gaps(env):
    wd = env.get("windDir", [1.0, 0.0])
    fl = env.get("flowPerSec", [0.0, 0.0])
    ang_w = math.degrees(math.atan2(wd[1], wd[0]))
    ang_f = (math.degrees(math.atan2(fl[1], fl[0]))
             if (fl[0] or fl[1]) else None)
    return {
        "wind_dir_vec": [float(wd[0]), float(wd[1])],
        "wind_dir_deg": float(ang_w),
        "flow_per_sec_vec": [float(fl[0]), float(fl[1])],
        "flow_axis_deg": (None if ang_f is None else float(ang_f)),
        "wind_minus_flow_deg": (None if ang_f is None
                                else float(ang_w - ang_f)),
        "flow_is_consumed_by_offline_composite": False,
        "flow_note": ("env.resolve_environment computes flowPerSec from "
                      "windWave.flowSpeed, but no offline module reads it "
                      "(grep: only env.py defines it), so the stream current "
                      "does not enter the caustics at all."),
    }


def fmt(x, n=4):
    if x is None:
        return "-"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return "%d" % x
    return ("%%.%df" % n) % float(x)


def main():
    t_all = time.time()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)

    with open(RECIPE, "r", encoding="utf-8") as f:
        raw = json.load(f)
    seed = derive_seed(raw, 1)
    log("base recipe=%s seed=%d" % (raw.get("name"), seed))

    base_np = np.asarray(Image.open(PHOTO).convert("RGB"))
    PH, PW, _ = base_np.shape
    photo_lin = srgb_to_linear(base_np.astype(np.float64) / 255.0)
    log("base frame %dx%d sha=%s" % (PW, PH, sha256_file(PHOTO)[:16]))

    # ==================================================================
    # 01 -- scene + SMB numbers (no transport)
    # ==================================================================
    scene = scene_for(raw, STREAM, "stream")
    env = resolve_environment(scene)
    cfg = parse_recipe_v8(scene)
    smb = smb_node(env, cfg["fetch"]["fetch_m"])
    log("[01] scene: wind=%.2f depth=%.2f fetch=%.1f Ft=%.1f f_p=%.4f Hz "
        "Hs_law=%.4f m s_max=%.2f"
        % (env["windSpeedMS"], env["depthM"], cfg["fetch"]["fetch_m"],
           smb["Ft"], smb["f_p_Hz"], smb["Hs_fetch_law_m"],
           smb["s_max_mitsuyasu"]))
    log("     sun alt=%.2f az=%.2f sunI=%.4f specL=%.4f S_patch=%.4f m"
        % (env["sun"]["altitudeDeg"], env["sun"]["azimuthDeg"], env["sunI"],
           env["specL"], 4.0 * env["specL"]))
    CP.unit_done("01-scene-smb")

    # ==================================================================
    # 02/03/04 -- stream state, master (Q=25), composite
    # ==================================================================
    cmaster_path = os.path.join(WORK, "C-stream-Q25.npy")
    C_reuse = None
    done = (CP.state() or {}).get("completed_ids") or []
    if "03-transport-Q25" in done and os.path.exists(cmaster_path):
        try:
            with open(cmaster_path, "rb") as f:
                buf = f.read()
            if len(buf) == GRID * GRID * 4:
                C_reuse = np.frombuffer(buf, dtype=np.float32).reshape(
                    GRID, GRID).copy()
                log("[02/03] reusing cached master %s" % cmaster_path)
        except Exception as exc:
            log("[02/03] cached master unusable (%s); recomputing" % exc)
    st = build(raw, scene, "stream", GRID, seed, photo_lin, PH, PW, C=C_reuse)
    CP.unit_done("02-state")
    if C_reuse is None:
        atomic_write_bytes(cmaster_path,
                           np.ascontiguousarray(st["C"]).tobytes())
    CP.unit_done("03-transport-Q25")

    png_bytes = io.BytesIO()
    Image.fromarray(st["still"]).save(png_bytes, format="PNG")
    png_bytes = png_bytes.getvalue()

    def verify_png(path):
        try:
            im = Image.open(path)
            im.load()
            return im.size == (PW, PH) and im.mode == "RGB"
        except Exception:
            return False

    atomic_write_bytes(os.path.join(OUT, "comp-stream-still.png"),
                       png_bytes, verify=verify_png)
    log("[04] wrote comp-stream-still.png (%dx%d, %d bytes)"
        % (PW, PH, len(png_bytes)))
    CP.unit_done("04-composite")

    # ==================================================================
    # 05 -- pond comparison with the SAME code and seed (realised values)
    # ==================================================================
    pond_scene = scene_for(raw, POND, "pond")
    pond = build(raw, pond_scene, "pond", GRID, seed, photo_lin, PH, PW)
    pond["smb"] = smb_node(pond["env"], pond["cfg"]["fetch"]["fetch_m"])
    log("[05] pond realised: Hs=%.4f sigma2=%.6f s_max=%.2f crush=%.5f"
        % (pond["realised"]["Hs_m"], pond["realised"]["sigma2_rad2"],
           pond["smb"]["s_max_mitsuyasu"],
           pond["still_metrics"]["crush_rate_linear_lt_0.02"]))
    del pond["C"], pond["fields"]
    gc.collect()
    CP.unit_done("05-pond-comparison")

    # ==================================================================
    # 06 -- gates + deliverables
    # ==================================================================
    still_png = np.asarray(Image.open(
        os.path.join(OUT, "comp-stream-still.png")).convert("RGB"))
    dif = np.abs(still_png.astype(np.int16) - base_np.astype(np.int16))
    gates = {
        "output_832x480": {"width": int(still_png.shape[1]),
                           "height": int(still_png.shape[0]),
                           "pass": bool(still_png.shape[:2] == (PH, PW))},
        "not_bit_identical_to_base": {
            "base_sha256": sha256_file(PHOTO),
            "still_sha256": sha256_file(
                os.path.join(OUT, "comp-stream-still.png")),
            "max_abs_diff_uint8": int(dif.max()),
            "mean_abs_diff_uint8": float(dif.mean()),
            "changed_pixel_fraction": float(np.mean(dif.max(axis=2) > 0)),
            "pass": bool(int(dif.max()) > 0),
        },
        "C_mean_1": {"mean": st["C_stats"]["mean"],
                     "abs_dev": st["C_stats"]["mean_abs_dev_from_1"],
                     "tolerance": CMEAN_TOL,
                     "pass": bool(st["C_stats"]["pass_mean_1"])},
        "crush_lt_2pct": {
            "crush": st["still_metrics"]["crush_rate_linear_lt_0.02"],
            "bound": CRUSH_BOUND,
            "pass": bool(st["still_metrics"]["crush_rate_linear_lt_0.02"]
                         < CRUSH_BOUND)},
        "blow_zero": {"blow": st["still_metrics"]
                      ["blowout_rate_aces_gt_0.98"],
                      "bound": BLOW_BOUND,
                      "pass": bool(st["still_metrics"]
                                   ["blowout_rate_aces_gt_0.98"]
                                   <= BLOW_BOUND)},
    }
    all_pass = all(g["pass"] for g in gates.values())
    log("[06] gates: %s (all=%s)"
        % ({k: v["pass"] for k, v in gates.items()}, all_pass))

    env_s, stt = st["env"], st["state"]
    gaps = stream_gaps(env_s)
    direction = stt.direction or {}

    numbers = {
        "task": "BRIEF-STREAM-STILL / one still, v3.4 (v8) pipeline, eta=1.0",
        "created_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat(),
        "pipeline": {
            "name": "v3.4 = v8 chunk 4",
            "runner": "run_stream_still.py (this run; imports run_v8's modules)",
            "state": "v8_state.V8WaterState",
            "transport": "det.det_frame (deterministic forward transport)",
            "composite": "v5_composite.render_frame_v5 (frozen production)",
            "eta": ETA,
            "eta_note": ("eta = 1.0 = direct sun only.  A recipe with no sky "
                         "block is eta = 1 exactly, so the v8 direct route is "
                         "used and the v9/Perez diffuse route is NOT invoked."),
            "master_grid": [GRID, GRID],
            "Q": Q,
            "A_ambient_scale": A_COMP,
            "B_lowpass_sigma_px": B_LOWPASS,
            "E_c_layer_exposure": E_COMP,
            "composite_note": ("E is the frozen production constant 2.2, not "
                               "the v8 display-level normalisation E_v8."),
        },
        "base_frame": {"path": PHOTO, "size": [int(PW), int(PH)],
                       "sha256": sha256_file(PHOTO)},
        "recipe_used": {
            "source": RECIPE,
            "name": raw.get("name"),
            "schema_out": scene.get("schema"),
            "fetch_m": float(cfg["fetch"]["fetch_m"]),
            "fetch_source": cfg["fetch"]["source"],
            "wind_ms": float(env_s["windSpeedMS"]),
            "wind_dir_deg": float(raw["windWave"]["windDirDeg"]),
            "depth_m": float(env_s["depthM"]),
            "rain_rate_mm_h": 0.0,
            "surface_film": "none",
            "sun": {"date": raw["sun"]["date"], "time": raw["sun"]["time"],
                    "lat_deg": raw["sun"]["latitudeDeg"],
                    "lon_deg": raw["sun"]["longitudeDeg"],
                    "altitude_deg": float(env_s["sun"]["altitudeDeg"]),
                    "azimuth_deg": float(env_s["sun"]["azimuthDeg"])},
            "weather": {"cloud_pct": raw["weather"]["cloudPct"],
                        "clarity": raw["weather"]["clarity"],
                        "diffuse_fraction": float(env_s["diffuse"])},
            "water": {"transparency": raw["water"]["transparency"],
                      "turbidity": raw["water"]["turbidity"],
                      "pigment": raw["water"]["pigment"]},
            "seed": int(stt.seed),
            "seed_note": ("derived from the UNMODIFIED recipe, so the phase "
                          "realisation is the one the piece's 3.x masters used"),
        },
        "smb_fetch_computation": smb,
        "realised": {
            "Hs_m": st["realised"]["Hs_m"],
            "Hs_target_from_anchor_m": stt.anchor.get(
                "Hs_after_spectral") if stt.anchor else None,
            "Hs_before_anchor_m": stt.anchor.get("Hs_before") if stt.anchor
            else None,
            "m0_m2": st["realised"]["m0_m2"],
            "sigma2_rad2": st["realised"]["sigma2_rad2"],
            "sigma2_spectral_rad2": stt.spectrum_after["sigma2_rad2"],
            "rms_slope_rad": st["realised"]["rms_slope"],
            "sigma2_target_rad2": (stt.anchor or {}).get(
                "sigma2_target_rad2"),
            "sigma2_relative_residual": (stt.anchor or {}).get(
                "sigma2_residual_rel"),
            "anchor_applied": (stt.anchor or {}).get("applied"),
            "taper_kc_rad_per_m": (stt.anchor or {}).get("kc_rad_per_m"),
            "Hs_consistency_check": {
                "integrated_Hs_m": stt.anchor["candidates"]["JONSWAP"]["Hs"]
                if stt.anchor and stt.anchor.get("candidates") else None,
                "jonswap_fetch_law_Hs_m": smb["Hs_fetch_law_m"],
            },
        },
        "direction": {
            "model": direction.get("model"),
            "s_max": direction.get("s_max"),
            "U_over_cp": direction.get("U_over_cp"),
            "s_min_over_band": direction.get("s_min_over_band"),
            "s_max_over_band": direction.get("s_max_over_band"),
            "half_width_deg_at_peak": (
                2.0 * math.degrees(math.acos(0.5 ** (1.0 /
                                                     (2.0 * direction["s_max"]))))
                if direction.get("s_max") else None),
            "sigma2_before": direction.get("sigma2_before"),
            "sigma2_after": direction.get("sigma2_after"),
            "height_renormalisation": direction.get("height_renormalisation"),
            "alignment": gaps,
            "narrowing_note": ("the code exposes NO free narrowing knob: "
                               "s_max = 11.5 (U10/c_p) is derived from the "
                               "fetch-limited peak.  Narrowing 'toward the "
                               "flow axis' was done by (a) putting the wind "
                               "along the flow axis and (b) taking the fetch "
                               "1000 m peak, which yields the narrowest lobe "
                               "this code honestly supports (s_max above)."),
        },
        "capillary_tail": {"k_viscous_rad_per_m":
                           stt.tail["k_viscous_rad_per_m"],
                           "k_nyquist_rad_per_m":
                           stt.tail["k_nyquist_rad_per_m"],
                           "kmin_rad_per_m": stt.tail["kmin_rad_per_m"],
                           "modes_removed_in_resolved_band":
                           stt.tail["modes_removed_in_resolved_band"]},
        "breaking": stt.breaking,
        "dispersion_quantisation":
            stt.dispersion_node()["quantisation"],
        "transport": {"Q": Q, "info": st["det_info"],
                      "C": st["C_stats"]},
        "still": {
            "file": "comp-stream-still.png",
            "size": [int(PW), int(PH)],
            "metrics": st["still_metrics"],
            "beer": st["beer"],
            "blur_sigma_px": st["blur_sigma_px"],
            "g1": st["g1"],
            "substrate_gain_rgb": st["substrate_gain_rgb"],
        },
        "gates": gates,
        "gates_all_pass": bool(all_pass),
        "pond_comparison_same_code_same_seed": {
            "pond_used": POND,
            "stream_used": STREAM,
            "labels_note": ("both columns are REALISED values produced by "
                            "this run; no preset labels are quoted"),
            "smb": pond["smb"],
            "realised": {
                "Hs_m": pond["realised"]["Hs_m"],
                "m0_m2": pond["realised"]["m0_m2"],
                "sigma2_rad2": pond["realised"]["sigma2_rad2"],
                "rms_slope_rad": pond["realised"]["rms_slope"],
                "sigma2_spectral_rad2":
                    pond["state"].spectrum_after["sigma2_rad2"],
                "sigma2_target_rad2": (pond["state"].anchor or {}).get(
                    "sigma2_target_rad2"),
                "anchor_applied": (pond["state"].anchor or {}).get("applied"),
                "kc_rad_per_m": (pond["state"].anchor or {}).get(
                    "kc_rad_per_m"),
                "s_max": pond["smb"]["s_max_mitsuyasu"],
                "half_width_deg_at_peak":
                    2.0 * math.degrees(math.acos(
                        0.5 ** (1.0 / (2.0 * pond["smb"]["s_max_mitsuyasu"])))),
                "Ft": pond["smb"]["Ft"],
                "lambda_p_m": pond["smb"]["lambda_p_m"],
                "C": pond["C_stats"],
                "crush": pond["still_metrics"]["crush_rate_linear_lt_0.02"],
                "blow": pond["still_metrics"]
                        ["blowout_rate_aces_gt_0.98"],
                "breaking_applied":
                    pond["state"].breaking.get("applied"),
            },
        },
        "honest_gaps": [
            {"kind": "still_cannot_show",
             "item": "current advection / Doppler shift of the spectrum",
             "statement": ("there is no advection or Doppler term in the "
                           "spectrum, and a still cannot show it anyway: "
                           "advection is motion (a phase drift at omega -> "
                           "omega + U.k), not a surface snapshot.")},
            {"kind": "genuinely_missing",
             "item": "current-induced steepening (U.k / omega of the waves)",
             "statement": ("genuinely absent from the code: nothing couples "
                           "the stream current to the wave spectrum, so the "
                           "opposing-flow steepening that makes a stream "
                           "edge crisp is NOT in this still.")},
            {"kind": "genuinely_missing",
             "item": "stream current in the composite",
             "statement": ("env.flowPerSec is computed but never read by any "
                           "offline module (env.py is its only occurrence), "
                           "so the recorded flow vector is display metadata "
                           "only.")},
            {"kind": "genuinely_missing",
             "item": "H_s is not produced by the SMB law in this code",
             "statement": (
                 "the brief asks for H_s out of the SMB/fetch computation, "
                 "but v8 only computes that number (%.4f m law, %.4f m "
                 "spectral integral) and does not use it as the height: the "
                 "field's amplitude is env.waveHeightM (recipe heuristic; "
                 "H_s %.4f m before the slope taper), the direction node "
                 "renormalises height, and the anchor taper then removes "
                 "slope.  Realised H_s = %.4f m is the field's own 4*sqrt(m0) "
                 "— %.2fx the SMB value.  Nothing was hand-set; the gap is "
                 "structural."
                 % (smb["Hs_fetch_law_m"],
                    stt.anchor["candidates"]["JONSWAP"]["Hs"]
                    if stt.anchor and stt.anchor.get("candidates") else
                    float("nan"),
                    float(stt.anchor.get("Hs_before") or float("nan")),
                    st["realised"]["Hs_m"],
                    st["realised"]["Hs_m"] / smb["Hs_fetch_law_m"]))},
            {"kind": "model_shape",
             "item": "realised spectral peak sits at the patch scale",
             "statement": ("the frozen Phillips shape peaks below the "
                           "resolved band (1/L = %.2f rad/m < kmin = %.2f), "
                           "so the realised spectrum peaks AT kmin "
                           "(breaking node: peak_is_at_kmin = %s) and the "
                           "dominant wavelength is the patch side %.2f m, "
                           "not the JONSWAP lambda_p %.2f m.  The fetch "
                           "anchor moves sigma^2, not the peak shape, so the "
                           "caustic scale is coarser than a real 1 km sea."
                           % (1.0 / float(env_s["specL"]),
                              float(stt.kmin),
                              str(bool(stt.breaking.get("peak_is_at_kmin"))),
                              4.0 * float(env_s["specL"]),
                              smb["lambda_p_m"]))},
            {"kind": "out_of_scope",
             "item": "shoaling / refraction over a depth gradient",
             "statement": ("uniform depth only; a stream edge has a cross-"
                           "channel gradient, which is out of scope for "
                           "chunk 4 (needs grad h).")},
            {"kind": "resolution_limit",
             "item": "spectral patch = 4*specL (longest resolved wavelength)",
             "statement": ("the patch side is %.3f m (set by the wind-driven "
                           "specL, not by the fetch); the fetch-1000 m peak "
                           "(%.3f m) is resolved, but the long-wave part of a "
                           "1 km sea is not representable."
                           % (4.0 * env_s["specL"], smb["lambda_p_m"]))},
            {"kind": "in_fitted_range",
             "item": "JONSWAP fetch laws",
             "statement": ("Ft = %.0f is INSIDE the fitted range 1e3-1e5, so "
                           "the stream numbers do not carry the downward "
                           "extrapolation the 100 m pond does (Ft = %.0f)."
                           % (smb["Ft"], pond["smb"]["Ft"]))},
        ],
        "timing_s": {"wall_total": time.time() - t_all},
        "checkpoint": CP.state(),
    }
    numbers = coerce(numbers)
    atomic_write_bytes(os.path.join(OUT, "numbers-stream.json"),
                       json.dumps(numbers, indent=2, ensure_ascii=False).encode(
                           "utf-8"))

    report = render_report(numbers)
    atomic_write_bytes(os.path.join(OUT, "REPORT-STREAM.md"),
                       report.encode("utf-8"))
    log("[06] wrote numbers-stream.json + REPORT-STREAM.md")
    CP.unit_done("06-gates-outputs")

    ok, bad = verify_units(OUT, "*.png", verify_png)
    ok_md, bad_md = verify_units(OUT, "*.md")
    ok_json, bad_json = verify_units(OUT, "*.json")
    print()
    print("checkpoint:", CP.state())
    print("verify: png ok=%d bad=%d | md ok=%d bad=%d | json ok=%d bad=%d"
          % (len(ok), len(bad), len(ok_md), len(bad_md),
             len(ok_json), len(bad_json)))
    print("gates_all_pass=%s wall=%.1fs" % (all_pass, time.time() - t_all))


def render_report(N):
    cmp = N["pond_comparison_same_code_same_seed"]
    P4 = cmp["smb"]
    pr = cmp["realised"]
    R = N["realised"]
    dr = N["direction"]
    st_c = N["transport"]["C"]
    sm = N["smb_fetch_computation"]
    g = N["gates"]

    L = []
    A = L.append
    A("# REPORT-STREAM — one still: v3.4 (v8) pipeline, gentle stream edge")
    A("")
    A("Base frame `blender-piece8s/renders_v2/frame_0000.png` (832x480), "
      "composited with the v3.4 (= v8 chunk-4) caustics pipeline at "
      "**eta = 1.0** (direct sun only; the v9/Perez diffuse route is not "
      "used). Composite frozen as production: **A2.6-Blow + E2.2**, master "
      "2048², Q = 25, σ_lowpass = 3.0 px, σ_blur = %s px, G1 = %s."
      % (fmt(N["still"]["blur_sigma_px"], 2), fmt(N["still"]["g1"], 6)))
    A("")
    A("Output `comp-stream-still.png` sha256 `%s`… (base `%s`…), "
      "max |Δ| = %d/255, %.3f%% of pixels changed."
      % (g["not_bit_identical_to_base"]["still_sha256"][:16],
         g["not_bit_identical_to_base"]["base_sha256"][:16],
         g["not_bit_identical_to_base"]["max_abs_diff_uint8"],
         100.0 * g["not_bit_identical_to_base"]["changed_pixel_fraction"]))
    A("")
    A("## 1. Numbers actually used")
    A("")
    A("| input | value | note |")
    A("|---|---|---|")
    A("| fetch | %s m | `wave.fetchM` explicit (pond default 100 m overridden) |"
      % fmt(N["recipe_used"]["fetch_m"], 1))
    A("| wind | %s m/s | along the corridor (`windDirDeg` %s → +x) |"
      % (fmt(N["recipe_used"]["wind_ms"], 2),
         fmt(N["recipe_used"]["wind_dir_deg"], 0)))
    A("| depth | %s m | shallow stream edge |"
      % fmt(N["recipe_used"]["depth_m"], 2))
    A("| rain | 0 mm/h | rain path untouched (F table not used) |")
    A("| surface film | none | no film slot in the code; nothing added |")
    A("| eta | 1.0 | legacy (no sky block) → direct only |")
    A("| sun | %s %s, %.1f°N %.1f°E → alt %s°, az %s° | preset-1 sun/sky family |"
      % (N["recipe_used"]["sun"]["date"], N["recipe_used"]["sun"]["time"],
         N["recipe_used"]["sun"]["lat_deg"],
         N["recipe_used"]["sun"]["lon_deg"],
         fmt(N["recipe_used"]["sun"]["altitude_deg"], 2),
         fmt(N["recipe_used"]["sun"]["azimuth_deg"], 2)))
    A("| water | turbidity %s, transparency %s, pigment %s | from recipe-1 |"
      % (N["recipe_used"]["water"]["turbidity"],
         N["recipe_used"]["water"]["transparency"],
         N["recipe_used"]["water"]["pigment"]))
    A("| grid / Q | 2048² / %d | det_frame, no clamp, no renormalisation |"
      % N["pipeline"]["Q"])
    A("")
    A("## 2. Realised values (this run)")
    A("")
    A("H_s is **computed, never set** — but the pipeline does not derive the "
      "height from the SMB law, and that gap is stated below rather than "
      "papered over: the SMB/JONSWAP value and the pipeline's own spectral "
      "integral are reported next to the realised field value.")
    A("")
    A("| quantity | value | where it comes from |")
    A("|---|---|---|")
    A("| Ft = gF/U² | %s | SMB/fetch input |" % fmt(sm["Ft"], 1))
    A("| f_p | %s Hz | JONSWAP `f_p = 3.5 (g/U) Ft^-1/3` |"
      % fmt(sm["f_p_Hz"], 4))
    A("| α (Phillips) | %s | `0.076 Ft^-0.22` |" % fmt(sm["alpha"], 5))
    A("| λ_p | %s m | from f_p |" % fmt(sm["lambda_p_m"], 3))
    A("| H_s (SMB fetch law) | %s m | `0.0016 Ft^0.5 U²/g` |"
      % fmt(sm["Hs_fetch_law_m"], 4))
    A("| H_s (pipeline spectral integral) | %s m | anchor candidates, "
      "JONSWAP/Donelan median |"
      % fmt((R["Hs_consistency_check"] or {}).get("integrated_Hs_m"), 4))
    A("| **H_s realised (4√m0)** | **%s m** | height field at t=0 |"
      % fmt(R["Hs_m"], 4))
    A("| H_s after the slope taper | %s m | before anchor: %s m |"
      % (fmt(R["Hs_target_from_anchor_m"], 4),
         fmt(R["Hs_before_anchor_m"], 4)))
    A("| σ² realised | %s rad² | mean(gx²+gy²) |" % fmt(R["sigma2_rad2"], 6))
    A("| σ² spectral (target) | %s rad² | anchor taper residual %.2e |"
      % (fmt(R["sigma2_target_rad2"], 6),
         float(R["sigma2_relative_residual"])))
    A("| rms slope | %s rad (%s°) | |"
      % (fmt(R["rms_slope_rad"], 4),
         fmt(math.degrees(math.atan(R["rms_slope_rad"])), 2)))
    A("| taper k_c | %s rad/m | λ_c = %s m |"
      % (fmt(R["taper_kc_rad_per_m"], 2),
         fmt(2.0 * math.pi / R["taper_kc_rad_per_m"], 3)
         if R["taper_kc_rad_per_m"] else "-"))
    A("| direction s_max | %s | Mitsuyasu `11.5 U/c_p` |"
      % fmt(dr["s_max"], 2))
    A("| angular half-width at peak | %s° | from s_max, no free knob |"
      % fmt(dr["half_width_deg_at_peak"], 1))
    A("| breaking clamp | %s | H_s %s m vs H_max %s m (%s binds) |"
      % ("applied" if N["breaking"].get("applied") else "not applied",
         fmt(N["breaking"].get("Hs_before_m"), 4),
         fmt(N["breaking"].get("H_max_m"), 4),
         N["breaking"].get("binding")))
    A("")
    A("C master (2048², Q=25): mean %s, p50 %s, p99 %s, p99.9 %s, max %s, "
      "min %s." % (fmt(st_c["mean"], 9), fmt(st_c["p50"], 4),
                   fmt(st_c["p99"], 4), fmt(st_c["p99_9"], 4),
                   fmt(st_c["max"], 4), fmt(st_c["min"], 4)))
    A("Still: crush %s %%, blow %s %%."
      % (fmt(100.0 * N["still"]["metrics"]["crush_rate_linear_lt_0.02"], 3),
         fmt(100.0 * N["still"]["metrics"]["blowout_rate_aces_gt_0.98"], 3)))
    A("")
    A("## 3. Pond (preset 1) vs stream — realised, same code, same seed")
    A("")
    A("Both columns were produced **in this run** by the same modules and "
      "the same recipe seed; only the three environment numbers differ. No "
      "preset label is quoted as a value.")
    A("")
    A("| realised quantity | pond (F 100 m, U 1.5, z 0.60) | stream (F 1000 m, U 2.5, z 0.30) | stream / pond |")
    A("|---|---|---|---|")
    rows = [
        ("Ft", P4["Ft"], sm["Ft"], 1),
        ("f_p [Hz]", P4["f_p_Hz"], sm["f_p_Hz"], 4),
        ("λ_p [m]", P4["lambda_p_m"], sm["lambda_p_m"], 3),
        ("s_max", P4["s_max_mitsuyasu"], sm["s_max_mitsuyasu"], 2),
        ("half-width at peak [°]",
         2.0 * math.degrees(math.acos(0.5 ** (1.0 / (2.0 * P4["s_max_mitsuyasu"])))),
         dr["half_width_deg_at_peak"], 1),
        ("H_s realised [m]", pr["Hs_m"], R["Hs_m"], 4),
        ("σ² spectral target [rad²]", pr["sigma2_target_rad2"],
         R["sigma2_target_rad2"], 6),
        ("σ² realised [rad²]", pr["sigma2_rad2"], R["sigma2_rad2"], 6),
        ("rms slope [rad]", pr["rms_slope_rad"], R["rms_slope_rad"], 4),
        ("k_c taper [rad/m]", pr["kc_rad_per_m"], R["taper_kc_rad_per_m"], 2),
        ("C mean", pr["C"]["mean"], st_c["mean"], 9),
        ("C p50", pr["C"]["p50"], st_c["p50"], 4),
        ("C p99", pr["C"]["p99"], st_c["p99"], 4),
        ("C p99.9", pr["C"]["p99_9"], st_c["p99_9"], 4),
        ("C max", pr["C"]["max"], st_c["max"], 4),
        ("crush [%]", 100.0 * pr["crush"],
         100.0 * N["still"]["metrics"]["crush_rate_linear_lt_0.02"], 3),
        ("blow [%]", 100.0 * pr["blow"],
         100.0 * N["still"]["metrics"]["blowout_rate_aces_gt_0.98"], 3),
    ]
    for name, a, b, n in rows:
        if a is None or b is None:
            ratio = "-"
        elif abs(float(a)) > 1e-12:
            ratio = "%.2f×" % (float(b) / float(a))
        else:
            ratio = "-"
        A("| %s | %s | %s | %s |" % (name, fmt(a, n), fmt(b, n), ratio))
    A("")
    A("Reading (measured, no mechanism claimed beyond these numbers): the "
      "stream is not simply 'more pond'. Its realised H_s is LARGER "
      "(%s m vs %s m) and its rms slope steeper (%s vs %s rad), yet its "
      "caustic tail is LOWER — C p99 %s vs %s (%.2fx), p99.9 %s vs %s, "
      "max %s vs %s — and it crushes less (%.3f%% vs %.3f%%). Under 0.30 m "
      "of water the longer, coarser field (λ_p %s m, peak at the patch "
      "scale, s_max %s → %s° lobe) focuses less sharply than the pond's "
      "short-wave field (λ_p %s m, %s° lobe)."
      % (fmt(R["Hs_m"], 4), fmt(pr["Hs_m"], 4),
         fmt(R["rms_slope_rad"], 4), fmt(pr["rms_slope_rad"], 4),
         fmt(st_c["p99"], 4), fmt(pr["C"]["p99"], 4),
         float(st_c["p99"]) / float(pr["C"]["p99"]),
         fmt(st_c["p99_9"], 4), fmt(pr["C"]["p99_9"], 4),
         fmt(st_c["max"], 4), fmt(pr["C"]["max"], 4),
         100.0 * N["still"]["metrics"]["crush_rate_linear_lt_0.02"],
         100.0 * pr["crush"],
         fmt(sm["lambda_p_m"], 3), fmt(sm["s_max_mitsuyasu"], 2),
         fmt(dr["half_width_deg_at_peak"], 1),
         fmt(P4["lambda_p_m"], 3),
         fmt(2.0 * math.degrees(math.acos(
             0.5 ** (1.0 / (2.0 * P4["s_max_mitsuyasu"])))), 1)))
    A("")
    A("## 4. Gates")
    A("")
    A("| gate | realised | pass |")
    A("|---|---|---|")
    A("| 832x480, opens | %dx%d RGB | %s |"
      % (g["output_832x480"]["width"], g["output_832x480"]["height"],
         fmt(g["output_832x480"]["pass"])))
    A("| not bit-identical to base | %d px differ (max Δ %d) | %s |"
      % (int(round(g["not_bit_identical_to_base"]["changed_pixel_fraction"]
                   * g["output_832x480"]["width"]
                   * g["output_832x480"]["height"])),
         g["not_bit_identical_to_base"]["max_abs_diff_uint8"],
         fmt(g["not_bit_identical_to_base"]["pass"])))
    A("| C mean = 1.0 | %s (|Δ| %.2e ≤ %s) | %s |"
      % (fmt(st_c["mean"], 9),
         float(g["C_mean_1"]["abs_dev"]),
         fmt(g["C_mean_1"]["tolerance"], 3),
         fmt(g["C_mean_1"]["pass"])))
    A("| crush < 2 %% | %s %% | %s |"
      % (fmt(100.0 * g["crush_lt_2pct"]["crush"], 3),
         fmt(g["crush_lt_2pct"]["pass"])))
    A("| blow = 0 | %s %% | %s |"
      % (fmt(100.0 * g["blow_zero"]["blow"], 3),
         fmt(g["blow_zero"]["pass"])))
    A("")
    A("All gates pass: **%s**." % fmt(N["gates_all_pass"]))
    A("")
    A("## 5. Honest gaps")
    A("")
    for item in N["honest_gaps"]:
        A("- **%s** (%s): %s" % (item["item"], item["kind"].replace("_", " "),
                                 item["statement"]))
    A("")
    A("## 6. Files")
    A("")
    A("- `comp-stream-still.png` — 832x480 composited still")
    A("- `numbers-stream.json` — recipe, realised values, gates, the pond "
      "side-by-side")
    A("- `REPORT-STREAM.md` — this file")
    A("- `checkpoint.json` — jobkit progress (%d unit_done events recorded at "
      "the time these files were written, retries included; the final "
      "`06-gates-outputs` unit lands immediately after them)"
      % (N.get("checkpoint") or {}).get("completed_units", 0))
    A("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
