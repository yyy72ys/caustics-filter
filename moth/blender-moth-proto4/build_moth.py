"""build_moth.py - PROTO2: procedural emerald-geometer moth pair + bark scene.

Changes vs proto1 (all aimed at the "vector illustration" look of v1):
  * wings: ONE fan-shaped mesh per side. The outline is traced from the
    reference photo (leading-edge curve + outer/inner margin curve), so the
    triangular forewing and the round hindwing read as a single silhouette,
    the way a resting geometer actually looks. Proto1 had 4 separate
    leaf-like parts.
  * posture: wings rest almost FLAT on the bark (3 deg, was 16 deg V-pose),
    camber reduced.
  * bands: the equal diagonal striping is gone. Two transverse bands plus a
    broken subterminal line, built as curves of near-constant span fraction
    (parallel to the outer margin), bowed, noise-wobbled and noise-broken.
  * colour: sampled from the photo (median wing green #6e9c90, dark mottle
    #5f826c, pale #83afa4) instead of the bright mint of proto1, with
    low-frequency mottling.
  * contact shadow: silhouette-shaped. The fan outline is rasterised into
    an alpha texture (scanline fill + box blur), so the shadow follows the
    moth's actual footprint instead of a generic radial ellipse.

Unchanged from proto1: bark plate (inpainted photo, unlit emission),
choreography (left moth departs by f24, right moth stays), Eevee 64
samples, motion blur, 832x480, 48 frames @ 24 fps, fixed camera,
body/legs/antennae construction.

Run examples
  python run_blender.py build_moth.py
  python run_blender.py build_moth.py --still 0 --out _dbg/t.png
  python run_blender.py build_moth.py --render
"""
import bpy
import bmesh
import json
import math
import os
import sys
from mathutils import Vector

HERE = os.path.dirname(os.path.abspath(__file__))
BG_IMG = os.path.join(HERE, "background-bark-832x480.png")
BLEND_OUT = os.path.join(HERE, "proto2.blend")

RES_X, RES_Y = 832, 480
FPS = 24
F_START, F_END = 0, 47
VIEW_W = 3.0
VIEW_H = VIEW_W * RES_Y / float(RES_X)


def fx(nx):
    return (nx - 0.5) * VIEW_W


def fy(ny):
    return (0.5 - ny) * VIEW_H


# ---------------------------------------------------------------- args
ARGV = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []


def arg(name, default=None, cast=str):
    if name in ARGV:
        return cast(ARGV[ARGV.index(name) + 1])
    return default


DO_RENDER = "--render" in ARGV
STILL = arg("--still", None, int)
STILL_OUT = arg("--out", os.path.join(HERE, "_dbg", "still.png"))
SUN_STRENGTH = arg("--sun", 0.66, float)
AMBIENT = arg("--ambient", 0.66, float)
TONE = arg("--tone", 1.0, float)
QUICK = "--quick" in ARGV
SILHOUETTE = "--silhouette" in ARGV
OBLIQUE_DEG = arg("--oblique", None, float)
NOSHADOW = "--noshadow" in ARGV
if STILL is None and "still" in ARGV:
    try:
        _si = ARGV.index("still")
        _nxt = ARGV[_si + 1] if _si + 1 < len(ARGV) else None
        if _nxt is not None and not _nxt.startswith("--"):
            if _nxt.endswith(".png"):
                STILL = 0
                STILL_OUT = _nxt if os.path.isabs(_nxt) else os.path.join(HERE, _nxt)
            else:
                try:
                    STILL = int(_nxt)
                except Exception:
                    STILL = 0
                if _si + 2 < len(ARGV) and ARGV[_si + 2].endswith(".png"):
                    _o = ARGV[_si + 2]
                    STILL_OUT = _o if os.path.isabs(_o) else os.path.join(HERE, _o)
    except Exception:
        pass


# ---------------------------------------------------------------- utils
def wipe():
    bpy.ops.wm.read_factory_settings(use_empty=True)


def srgb(hexstr):
    h = hexstr.lstrip("#")
    v = [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]
    return tuple((c / 12.92) if c <= 0.04045 else (((c + 0.055) / 1.055) ** 2.4) for c in v)


def new_mesh_obj(name, verts, faces, coll, uv=None, smooth=True):
    me = bpy.data.meshes.new(name)
    me.from_pydata([tuple(v) for v in verts], [], faces)
    me.update()
    me.validate(verbose=False)
    if uv is not None:
        uvl = me.uv_layers.new(name="UVMap")
        for li, loop in enumerate(me.loops):
            uvl.data[li].uv = uv[loop.vertex_index]
    if smooth:
        for p in me.polygons:
            p.use_smooth = True
    ob = bpy.data.objects.new(name, me)
    coll.objects.link(ob)
    return ob


def set_in(node, name, value):
    if name in node.inputs:
        node.inputs[name].default_value = value


def smooth_profile(points, n=160):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    out = []
    for i in range(n):
        t = i / float(n - 1)
        y = ys[-1] if t > xs[-1] else ys[0]
        for k in range(len(xs) - 1):
            if xs[k] <= t <= xs[k + 1]:
                f = (t - xs[k]) / max(1e-9, (xs[k + 1] - xs[k]))
                y = ys[k] * (1 - f) + ys[k + 1] * f
                break
        out.append(y)
    sm = out[:]
    for i in range(1, n - 1):
        sm[i] = 0.25 * out[i - 1] + 0.5 * out[i] + 0.25 * out[i + 1]
    return sm


def sample_curve(points, n):
    """points: [(t, x, y), ...] with t monotone -> (xs, ys) sampled at n rows."""
    xs = smooth_profile([(p[0], p[1]) for p in points], n)
    ys = smooth_profile([(p[0], p[2]) for p in points], n)
    return xs, ys


def fan_outline(n):
    """Sample LE/TE then displace TE along local normal for hindwing scallops."""
    lex, ley = sample_curve(FAN_LE, n)
    tex, tey = sample_curve(FAN_TE, n)
    tex = list(tex)
    tey = list(tey)
    for i in range(1, n - 1):
        t = i / float(n - 1)
        if 0.17 < t < 0.60:
            w = math.sin(math.pi * (t - 0.17) / 0.43) ** 0.7
            d = 0.0070 * math.sin(2 * math.pi * 7.0 * (t - 0.17) / 0.43 + 0.5) * w
            tx = tex[i + 1] - tex[i - 1]
            ty = tey[i + 1] - tey[i - 1]
            leng = math.hypot(tx, ty) or 1e-9
            nx, ny = -ty / leng, tx / leng
            dx, dy = lex[i] - tex[i], ley[i] - tey[i]
            if nx * dx + ny * dy < 0:
                nx, ny = -nx, -ny
            tex[i] += nx * d
            tey[i] += ny * d
    return lex, ley, tex, tey


def tube(name, path, radii, coll, segs=8, caps=True):
    pts = [Vector(p) for p in path]
    n = len(pts)
    verts, faces = [], []
    up = Vector((0, 0, 1))
    for i, p in enumerate(pts):
        if i == 0:
            tan = pts[1] - pts[0]
        elif i == n - 1:
            tan = pts[-1] - pts[-2]
        else:
            tan = pts[i + 1] - pts[i - 1]
        tan.normalize()
        ref = up if abs(tan.dot(up)) < 0.95 else Vector((1, 0, 0))
        a = tan.cross(ref).normalized()
        b = tan.cross(a).normalized()
        r = radii[i]
        for k in range(segs):
            ang = 2 * math.pi * k / segs
            verts.append(p + a * (math.cos(ang) * r) + b * (math.sin(ang) * r))
    for i in range(n - 1):
        for k in range(segs):
            k2 = (k + 1) % segs
            faces.append([i * segs + k, i * segs + k2, (i + 1) * segs + k2, (i + 1) * segs + k])
    if caps:
        faces.append(list(range(segs - 1, -1, -1)))
        base = (n - 1) * segs
        faces.append([base + k for k in range(segs)])
    return new_mesh_obj(name, verts, faces, coll)


# ---------------------------------------------------------------- materials
def mat_background():
    m = bpy.data.materials.new("BarkPlate")
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    emi = nt.nodes.new("ShaderNodeEmission")
    tex = nt.nodes.new("ShaderNodeTexImage")
    img = bpy.data.images.load(BG_IMG, check_existing=True)
    tex.image = img
    tex.extension = "EXTEND"
    emi.inputs["Strength"].default_value = 1.0
    nt.links.new(tex.outputs["Color"], emi.inputs["Color"])
    nt.links.new(emi.outputs["Emission"], out.inputs["Surface"])
    out.location = (400, 0)
    emi.location = (150, 0)
    tex.location = (-150, 0)
    return m


def mat_wing():
    """UV.x = span (0 base -> 1 apex), UV.y = chord (0 costa -> 1 dorsum).

    Palette sampled from the reference photo (median of green wing pixels):
    base #6e9c90, dark mottle #5f826c / #55755f, pale #83afa4.
    PROTO4 CHANGE 1: palette scaled by the measured photo/render gap
    (x1.24 R, x1.29 G, x1.37 B in sRGB) so the RENDERED wing median lands on
    the photo median (112,157,148); mottle amplitude 0.16 -> 0.25 and the
    grain colour multiply widened 0.94+0.11g -> 0.90+0.20g.
    """
    m = bpy.data.materials.new("WingProcedural2")
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    N, L = nt.nodes, nt.links.new

    out = N.new("ShaderNodeOutputMaterial"); out.location = (1700, 0)
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (1450, 0)

    tc = N.new("ShaderNodeTexCoord"); tc.location = (-1500, 0)
    sep = N.new("ShaderNodeSeparateXYZ"); sep.location = (-1300, 300)
    L(tc.outputs["UV"], sep.inputs["Vector"])
    t_ = sep.outputs["X"]   # span
    u_ = sep.outputs["Y"]   # chord

    def math_node(op, loc):
        n = N.new("ShaderNodeMath"); n.operation = op; n.location = loc
        return n

    def map_range(loc, value_out, from_min, from_max, to_min, to_max, smooth=True):
        n = N.new("ShaderNodeMapRange"); n.location = loc
        n.inputs["From Min"].default_value = from_min
        n.inputs["From Max"].default_value = from_max
        n.inputs["To Min"].default_value = to_min
        n.inputs["To Max"].default_value = to_max
        if smooth:
            try:
                n.interpolation_type = "SMOOTHSTEP"
            except Exception:
                pass
        L(value_out, n.inputs["Value"])
        return n.outputs["Result"]

    def mix_rgba(loc, factor_out, a_out, b_out, blend="MIX"):
        n = N.new("ShaderNodeMix"); n.data_type = "RGBA"; n.blend_type = blend
        n.location = loc
        if factor_out is not None:
            L(factor_out, n.inputs["Factor"])
        L(a_out, n.inputs[6])
        L(b_out, n.inputs[7])
        return n.outputs[2]

    def const_rgb(loc, hexstr):
        n = N.new("ShaderNodeRGB"); n.location = loc
        n.outputs[0].default_value = (*srgb(hexstr), 1.0)
        return n.outputs[0]

    # ---------------- mottled base colour (low-frequency patches)
    nzA = N.new("ShaderNodeTexNoise"); nzA.location = (-1300, -60)
    nzA.inputs["Scale"].default_value = 3.0
    nzA.inputs["Detail"].default_value = 4.0
    nzA.inputs["Roughness"].default_value = 0.6
    L(tc.outputs["UV"], nzA.inputs["Vector"])
    rampA = N.new("ShaderNodeValToRGB"); rampA.location = (-1080, -60)
    ra = rampA.color_ramp
    ra.elements[0].position = 0.18; ra.elements[0].color = (*srgb("#6e9b87"), 1)
    ra.elements[1].position = 0.52; ra.elements[1].color = (*srgb("#8ecccc"), 1)
    e = ra.elements.new(0.86); e.color = (*srgb("#aae8e9"), 1)
    L(nzA.outputs["Fac"], rampA.inputs["Fac"])
    col = rampA.outputs["Color"]

    # ---------------- darker outer margin
    margF = map_range((-1080, 420), t_, 0.87, 1.0, 0.0, 0.34)
    col = mix_rgba((-860, 300), margF, col, const_rgb((-860, 460), "#7cac99"))

    # ---------------- fore/hind crease + hindwing tint
    # boundary curve in UV: u_line(t) = 0.309 + 0.690 t
    line = math_node("MULTIPLY_ADD", (-1080, 640))
    L(t_, line.inputs[0]); line.inputs[1].default_value = 0.690; line.inputs[2].default_value = 0.309
    dsub = math_node("SUBTRACT", (-900, 640))
    L(u_, dsub.inputs[0]); L(line.outputs[0], dsub.inputs[1])
    hindF0 = map_range((-720, 700), dsub.outputs[0], 0.0, 0.07, 0.0, 1.0)
    hindF = math_node("MULTIPLY", (-560, 700))
    L(hindF0, hindF.inputs[0]); hindF.inputs[1].default_value = 0.38
    col = mix_rgba((-560, 420), hindF.outputs[0], col, const_rgb((-560, 560), "#a1e3e2"))
    dabs = math_node("ABSOLUTE", (-720, 560)); L(dsub.outputs[0], dabs.inputs[0])
    crease = map_range((-560, 560), dabs.outputs[0], 0.045, 0.0, 0.0, 0.16)
    oneMC = math_node("SUBTRACT", (-380, 560))
    oneMC.inputs[0].default_value = 1.0; L(crease, oneMC.inputs[1])
    colMul = N.new("ShaderNodeMix"); colMul.data_type = "RGBA"; colMul.blend_type = "MULTIPLY"
    colMul.location = (-380, 420); colMul.inputs["Factor"].default_value = 1.0
    L(col, colMul.inputs[6])
    greyRGB = N.new("ShaderNodeCombineColor"); greyRGB.location = (-260, 300)
    greyRGB.mode = "RGB"
    L(oneMC.outputs[0], greyRGB.inputs[0]); L(oneMC.outputs[0], greyRGB.inputs[1]); L(oneMC.outputs[0], greyRGB.inputs[2])
    L(greyRGB.outputs["Color"], colMul.inputs[7])
    col = colMul.outputs[2]

    # ---------------- veins: subtle lines of constant chord fraction
    tsin = math_node("SINE", (-900, -320)); 
    tmul = math_node("MULTIPLY", (-1080, -320)); L(t_, tmul.inputs[0]); tmul.inputs[1].default_value = 7.0
    L(tmul.outputs[0], tsin.inputs[0])
    vphase = math_node("MULTIPLY_ADD", (-720, -320))
    L(u_, vphase.inputs[0]); vphase.inputs[1].default_value = 2.0 * math.pi * 6.5
    vsin = math_node("SINE", (-540, -320))
    vadd = math_node("ADD", (-640, -380))
    L(vphase.outputs[0], vadd.inputs[0])
    tsc = math_node("MULTIPLY", (-640, -260)); L(tsin.outputs[0], tsc.inputs[0]); tsc.inputs[1].default_value = 0.5
    L(tsc.outputs[0], vadd.inputs[1])
    L(vadd.outputs[0], vsin.inputs[0])
    vabs = math_node("ABSOLUTE", (-380, -320)); L(vsin.outputs[0], vabs.inputs[0])
    vein = N.new("ShaderNodeValToRGB"); vein.location = (-200, -320)
    vr = vein.color_ramp
    vr.interpolation = "EASE"
    vr.elements[0].position = 0.0; vr.elements[0].color = (1, 1, 1, 1)
    vr.elements[1].position = 0.30; vr.elements[1].color = (0, 0, 0, 1)
    L(vabs.outputs[0], vein.inputs["Fac"])
    vfade = map_range((-200, -480), t_, 0.62, 1.0, 1.0, 0.3, smooth=False)
    veinAmt = math_node("MULTIPLY", (-40, -400))
    L(vein.outputs["Color"], veinAmt.inputs[0]); L(vfade, veinAmt.inputs[1])
    veinAmt2 = math_node("MULTIPLY", (60, -400))
    L(veinAmt.outputs[0], veinAmt2.inputs[0]); veinAmt2.inputs[1].default_value = 0.30
    col = mix_rgba((60, 300), veinAmt2.outputs[0], col, const_rgb((60, 460), "#739b87"))

    # ---------------- STAGE 2-2 discal (cell-end) marks (before bands)
    def _gauss_mark(tc_val, t0, wt, uc_val, u0, wu, amp, ybase):
        ts = math_node("SUBTRACT", (-1300, ybase)); L(tc_val, ts.inputs[0]); ts.inputs[1].default_value = t0
        td = math_node("DIVIDE", (-1120, ybase)); L(ts.outputs[0], td.inputs[0]); td.inputs[1].default_value = wt
        tp = math_node("POWER", (-960, ybase)); L(td.outputs[0], tp.inputs[0]); tp.inputs[1].default_value = 2.0
        us = math_node("SUBTRACT", (-1300, ybase - 90)); L(uc_val, us.inputs[0]); us.inputs[1].default_value = u0
        ud = math_node("DIVIDE", (-1120, ybase - 90)); L(us.outputs[0], ud.inputs[0]); ud.inputs[1].default_value = wu
        up = math_node("POWER", (-960, ybase - 90)); L(ud.outputs[0], up.inputs[0]); up.inputs[1].default_value = 2.0
        ad = math_node("ADD", (-800, ybase - 60)); L(tp.outputs[0], ad.inputs[0]); L(up.outputs[0], ad.inputs[1])
        ng = math_node("MULTIPLY", (-640, ybase - 60)); L(ad.outputs[0], ng.inputs[0]); ng.inputs[1].default_value = -1.0
        ex = math_node("EXPONENT", (-480, ybase - 60)); L(ng.outputs[0], ex.inputs[0])
        am = math_node("MULTIPLY", (-320, ybase - 60)); L(ex.outputs[0], am.inputs[0]); am.inputs[1].default_value = amp
        return am.outputs[0]
    _m1 = _gauss_mark(t_, 0.50, 0.05, u_, 0.38, 0.10, 0.16, 220)
    _m2 = _gauss_mark(t_, 0.34, 0.045, u_, 0.78, 0.09, 0.13, 20)
    col = mix_rgba((-160, 220), _m1, col, const_rgb((-160, 380), "#629274"))
    col = mix_rgba((0, 20), _m2, col, const_rgb((0, 180), "#629274"))

    # ---------------- transverse bands (bowed, wobbling, broken)
    # shared end fade along the chord
    ef1 = map_range((-1300, -700), u_, 0.0, 0.10, 0.0, 1.0)
    ef2 = map_range((-1300, -820), u_, 1.0, 0.90, 0.0, 1.0)
    efm = math_node("MULTIPLY", (-1120, -760)); L(ef1, efm.inputs[0]); L(ef2, efm.inputs[1])
    endfade = math_node("MULTIPLY_ADD", (-960, -760))
    L(efm.outputs[0], endfade.inputs[0]); endfade.inputs[1].default_value = 0.70; endfade.inputs[2].default_value = 0.30

    # ---------------- STAGE 2-1 lunulate shared factors L(u), W(u)
    _lu = math_node("MULTIPLY", (-1300, -560)); L(u_, _lu.inputs[0]); _lu.inputs[1].default_value = math.pi * 6.5
    _ls = math_node("SINE", (-1120, -560)); L(_lu.outputs[0], _ls.inputs[0])
    _la = math_node("ABSOLUTE", (-960, -560)); L(_ls.outputs[0], _la.inputs[0])
    _lp = math_node("POWER", (-800, -560)); L(_la.outputs[0], _lp.inputs[0]); _lp.inputs[1].default_value = 0.8
    _lun = math_node("MULTIPLY_ADD", (-640, -560)); L(_lp.outputs[0], _lun.inputs[0]); _lun.inputs[1].default_value = 0.75; _lun.inputs[2].default_value = 0.25
    lunL = _lun.outputs[0]
    _ws = math_node("SUBTRACT", (-1300, -660)); L(u_, _ws.inputs[0]); _ws.inputs[1].default_value = 0.48
    _wq = math_node("DIVIDE", (-1120, -660)); L(_ws.outputs[0], _wq.inputs[0]); _wq.inputs[1].default_value = 0.30
    _wp = math_node("POWER", (-960, -660)); L(_wq.outputs[0], _wp.inputs[0]); _wp.inputs[1].default_value = 2.0
    _wn = math_node("MULTIPLY", (-800, -660)); L(_wp.outputs[0], _wn.inputs[0]); _wn.inputs[1].default_value = -1.0
    _we = math_node("EXPONENT", (-640, -660)); L(_wn.outputs[0], _we.inputs[0])
    _ww = math_node("MULTIPLY_ADD", (-480, -660)); L(_we.outputs[0], _ww.inputs[0]); _ww.inputs[1].default_value = 0.47; _ww.inputs[2].default_value = 0.85
    lunW = _ww.outputs[0]

    def build_band(t0, bow, width, strength, seed, mask_lo, mask_hi, floor_, y0):
        # centre line tc(u) = t0 + bow*sin(pi u) + wobble + ripple
        bowN = math_node("SINE", (-960, y0))
        bm = math_node("MULTIPLY", (-1120, y0)); L(u_, bm.inputs[0]); bm.inputs[1].default_value = math.pi
        L(bm.outputs[0], bowN.inputs[0])
        bowS = math_node("MULTIPLY", (-800, y0)); L(bowN.outputs[0], bowS.inputs[0]); bowS.inputs[1].default_value = bow
        # wobble noise, seeded by UV offset
        off = N.new("ShaderNodeVectorMath"); off.operation = "ADD"; off.location = (-1300, y0 - 90)
        L(tc.outputs["UV"], off.inputs[0]); off.inputs[1].default_value = (seed * 0.37, seed * 0.91, 0.0)
        nzW = N.new("ShaderNodeTexNoise"); nzW.location = (-1120, y0 - 90)
        nzW.inputs["Scale"].default_value = 4.0; nzW.inputs["Detail"].default_value = 3.0
        L(off.outputs[0], nzW.inputs["Vector"])
        wob = math_node("MULTIPLY_ADD", (-960, y0 - 90))
        L(nzW.outputs["Fac"], wob.inputs[0]); wob.inputs[1].default_value = 0.055; wob.inputs[2].default_value = -0.0275
        # fine ripple
        rip = math_node("SINE", (-800, y0 - 90))
        rm = math_node("MULTIPLY_ADD", (-960, y0 - 180))
        L(u_, rm.inputs[0]); rm.inputs[1].default_value = 44.0; rm.inputs[2].default_value = seed * 13.0
        L(rm.outputs[0], rip.inputs[0])
        ripS = math_node("MULTIPLY", (-640, y0 - 90)); L(rip.outputs[0], ripS.inputs[0]); ripS.inputs[1].default_value = 0.0035
        cen = math_node("ADD", (-640, y0)); L(bowS.outputs[0], cen.inputs[0]); L(wob.outputs[0], cen.inputs[1])
        cen2 = math_node("ADD", (-480, y0)); L(cen.outputs[0], cen2.inputs[0]); L(ripS.outputs[0], cen2.inputs[1])
        cen3 = math_node("ADD", (-320, y0)); L(cen2.outputs[0], cen3.inputs[0]); cen3.inputs[1].default_value = t0
        dist = math_node("ABSOLUTE", (-160, y0))
        ds = math_node("SUBTRACT", (-240, y0)); L(t_, ds.inputs[0]); L(cen3.outputs[0], ds.inputs[1])
        L(ds.outputs[0], dist.inputs[0])
        distW = math_node("DIVIDE", (-80, y0)); L(dist.outputs[0], distW.inputs[0]); L(lunW, distW.inputs[1])
        prof = map_range((0, y0), distW.outputs[0], 0.0, width, 1.0, 0.0)
        # breakup mask
        sc = N.new("ShaderNodeVectorMath"); sc.operation = "MULTIPLY"; sc.location = (-1300, y0 - 300)
        L(tc.outputs["UV"], sc.inputs[0]); sc.inputs[1].default_value = (5.0, 16.0, 0.0)
        off2 = N.new("ShaderNodeVectorMath"); off2.operation = "ADD"; off2.location = (-1120, y0 - 300)
        L(sc.outputs[0], off2.inputs[0]); off2.inputs[1].default_value = (seed * 1.3, seed * 2.1, 0.0)
        nzB = N.new("ShaderNodeTexNoise"); nzB.location = (-960, y0 - 300)
        nzB.inputs["Scale"].default_value = 1.0; nzB.inputs["Detail"].default_value = 2.0
        L(off2.outputs[0], nzB.inputs["Vector"])
        mask = map_range((-800, y0 - 300), nzB.outputs["Fac"], mask_lo, mask_hi, 0.0, 1.0)
        mfac = math_node("MULTIPLY_ADD", (-640, y0 - 300))
        L(mask, mfac.inputs[0]); mfac.inputs[1].default_value = (1.0 - floor_); mfac.inputs[2].default_value = floor_
        amt = math_node("MULTIPLY", (160, y0)); L(prof, amt.inputs[0]); L(mfac.outputs[0], amt.inputs[1])
        amt2 = math_node("MULTIPLY", (300, y0)); L(amt.outputs[0], amt2.inputs[0]); L(endfade.outputs[0], amt2.inputs[1])
        amt3 = math_node("MULTIPLY", (420, y0)); L(amt2.outputs[0], amt3.inputs[0]); amt3.inputs[1].default_value = strength
        amt4 = math_node("MULTIPLY", (500, y0)); L(amt3.outputs[0], amt4.inputs[0]); L(lunL, amt4.inputs[1])
        return amt4.outputs[0]

    b1 = build_band(0.58, 0.060, 0.026, 0.38, 0.0, 0.38, 0.62, 0.15, 120)
    b2 = build_band(0.765, 0.075, 0.030, 0.62, 3.7, 0.38, 0.62, 0.15, -60)
    b3 = build_band(0.895, 0.060, 0.020, 0.34, 7.3, 0.45, 0.60, 0.05, -240)
    bsum = math_node("ADD", (560, 60)); L(b1, bsum.inputs[0]); L(b2, bsum.inputs[1])
    bsum2 = math_node("ADD", (680, 60)); L(bsum.outputs[0], bsum2.inputs[0]); L(b3, bsum2.inputs[1])
    bclamp = math_node("MINIMUM", (800, 60)); L(bsum2.outputs[0], bclamp.inputs[0]); bclamp.inputs[1].default_value = 1.0
    bamt = math_node("MULTIPLY", (920, 60)); L(bclamp.outputs[0], bamt.inputs[0]); bamt.inputs[1].default_value = 0.75
    col = mix_rgba((920, 300), bamt.outputs[0], col, const_rgb((920, 460), "#edf2e9"))

    # ---------------- fine grain + speckle (kept from proto1, subtler)
    # STAGE 2 note: n2 non-bump roles (colour ramp, roughness base) kept;
    # bump Height below is replaced by scale grain g (Change 2-3).
    n2 = N.new("ShaderNodeTexNoise"); n2.location = (-1300, 900)
    n2.inputs["Scale"].default_value = 55.0
    n2.inputs["Detail"].default_value = 6.0
    n2.inputs["Roughness"].default_value = 0.6
    L(tc.outputs["Object"], n2.inputs["Vector"])
    grain = N.new("ShaderNodeValToRGB"); grain.location = (-1080, 900)
    gr = grain.color_ramp
    gr.elements[0].position = 0.25; gr.elements[0].color = (0.90, 0.92, 0.90, 1)
    gr.elements[1].position = 0.75; gr.elements[1].color = (1.04, 1.03, 1.02, 1)
    L(n2.outputs["Fac"], grain.inputs["Fac"])
    gmix = N.new("ShaderNodeMix"); gmix.data_type = "RGBA"; gmix.blend_type = "MULTIPLY"
    gmix.location = (1120, 300); gmix.inputs["Factor"].default_value = 0.8
    L(col, gmix.inputs[6]); L(grain.outputs["Color"], gmix.inputs[7])
    col = gmix.outputs[2]

    vor = N.new("ShaderNodeTexVoronoi"); vor.location = (-1300, 1080)
    vor.inputs["Scale"].default_value = 120.0
    L(tc.outputs["Object"], vor.inputs["Vector"])
    speck = N.new("ShaderNodeValToRGB"); speck.location = (-1080, 1080)
    speck.color_ramp.elements[0].position = 0.02
    speck.color_ramp.elements[0].color = (0, 0, 0, 1)
    speck.color_ramp.elements[1].position = 0.22
    speck.color_ramp.elements[1].color = (0.5, 0.5, 0.5, 1)
    L(vor.outputs["Distance"], speck.inputs["Fac"])
    smix = N.new("ShaderNodeMix"); smix.data_type = "RGBA"; smix.blend_type = "SCREEN"
    smix.location = (1260, 300); smix.inputs["Factor"].default_value = 0.08
    L(col, smix.inputs[6]); L(speck.outputs["Color"], smix.inputs[7])
    col = smix.outputs[2]

    # ---------------- soft fresnel rim (STAGE 2-5: factor 0.08 -> 0.035)
    lw = N.new("ShaderNodeLayerWeight"); lw.location = (920, -560)
    lw.inputs["Blend"].default_value = 0.30
    rim = N.new("ShaderNodeValToRGB"); rim.location = (1080, -560)
    rr = rim.color_ramp
    rr.elements[0].position = 0.62; rr.elements[0].color = (0, 0, 0, 1)
    rr.elements[1].position = 1.0; rr.elements[1].color = (1, 1, 1, 1)
    L(lw.outputs["Facing"], rim.inputs["Fac"])
    rimw = math_node("MULTIPLY", (1240, -560)); L(rim.outputs["Color"], rimw.inputs[0]); rimw.inputs[1].default_value = 0.035
    col = mix_rgba((1400, 300), rimw.outputs[0], col, const_rgb((1240, -700), "#f2f8f0"))

    # ---------------- STAGE 2-3 scale grain g (anisotropic ~1/500 wingspan)
    _gmap = N.new("ShaderNodeMapping"); _gmap.location = (-1300, 700)
    _gmap.inputs["Scale"].default_value = (110.0, 200.0, 1.0)
    L(tc.outputs["UV"], _gmap.inputs["Vector"])
    _gn = N.new("ShaderNodeTexNoise"); _gn.location = (-1120, 700)
    _gn.inputs["Scale"].default_value = 1.0
    _gn.inputs["Detail"].default_value = 2.0
    _gn.inputs["Roughness"].default_value = 0.7
    L(_gmap.outputs["Vector"], _gn.inputs["Vector"])
    grainG = _gn.outputs["Fac"]

    # ---------------- STAGE 2-4 fade mottle + base-to-tip gradient (before grain multiply)
    _fmap = N.new("ShaderNodeMapping"); _fmap.location = (-1300, 500)
    _fmap.inputs["Scale"].default_value = (3.0, 2.0, 1.0)
    L(tc.outputs["UV"], _fmap.inputs["Vector"])
    _fn = N.new("ShaderNodeTexNoise"); _fn.location = (-1120, 500)
    L(_fmap.outputs["Vector"], _fn.inputs["Vector"])
    fadeS = map_range((-960, 500), _fn.outputs["Fac"], 0.55, 0.80, 0.0, 1.0)
    fade16 = math_node("MULTIPLY", (-800, 500)); L(fadeS, fade16.inputs[0]); fade16.inputs[1].default_value = 0.25
    hindW = math_node("MULTIPLY_ADD", (-800, 620)); L(hindF0, hindW.inputs[0]); hindW.inputs[1].default_value = 0.5; hindW.inputs[2].default_value = 1.0
    fadeF = math_node("MULTIPLY", (-640, 500)); L(fade16.outputs[0], fadeF.inputs[0]); L(hindW.outputs[0], fadeF.inputs[1])
    col = mix_rgba((-480, 500), fadeF.outputs[0], col, const_rgb((-480, 660), "#ccded2"))
    radS = map_range((-960, 340), t_, 0.16, 0.90, 0.0, 1.0)
    radF = math_node("MULTIPLY_ADD", (-800, 340)); L(radS, radF.inputs[0]); radF.inputs[1].default_value = 0.10; radF.inputs[2].default_value = 0.90
    radGrey = N.new("ShaderNodeCombineColor"); radGrey.location = (-640, 340)
    radGrey.mode = "RGB"
    L(radF.outputs[0], radGrey.inputs[0]); L(radF.outputs[0], radGrey.inputs[1]); L(radF.outputs[0], radGrey.inputs[2])
    radMul = N.new("ShaderNodeMix"); radMul.data_type = "RGBA"; radMul.blend_type = "MULTIPLY"
    radMul.location = (-480, 340); radMul.inputs["Factor"].default_value = 1.0
    L(col, radMul.inputs[6]); L(radGrey.outputs["Color"], radMul.inputs[7])
    col = radMul.outputs[2]

    # ---------------- STAGE 2-3 grain colour multiply (just before BSDF)
    grainF = math_node("MULTIPLY_ADD", (-640, 200)); L(grainG, grainF.inputs[0]); grainF.inputs[1].default_value = 0.30; grainF.inputs[2].default_value = 0.85
    grainGrey = N.new("ShaderNodeCombineColor"); grainGrey.location = (-480, 200)
    grainGrey.mode = "RGB"
    L(grainF.outputs[0], grainGrey.inputs[0]); L(grainF.outputs[0], grainGrey.inputs[1]); L(grainF.outputs[0], grainGrey.inputs[2])
    grainMul = N.new("ShaderNodeMix"); grainMul.data_type = "RGBA"; grainMul.blend_type = "MULTIPLY"
    grainMul.location = (-320, 300); grainMul.inputs["Factor"].default_value = 1.0
    L(col, grainMul.inputs[6]); L(grainGrey.outputs["Color"], grainMul.inputs[7])
    col = grainMul.outputs[2]

    # ---------- PROTO4 CHANGE 1b: fine grain octave (~1px at render scale)
    _gmap2 = N.new("ShaderNodeMapping"); _gmap2.location = (-1300, 1180)
    _gmap2.inputs["Scale"].default_value = (420.0, 760.0, 1.0)
    L(tc.outputs["UV"], _gmap2.inputs["Vector"])
    _gn2 = N.new("ShaderNodeTexNoise"); _gn2.location = (-1120, 1180)
    _gn2.inputs["Scale"].default_value = 1.0
    _gn2.inputs["Detail"].default_value = 2.0
    _gn2.inputs["Roughness"].default_value = 0.7
    L(_gmap2.outputs["Vector"], _gn2.inputs["Vector"])
    grainG2 = _gn2.outputs["Fac"]
    grainF2 = math_node("MULTIPLY_ADD", (-640, 60)); L(grainG2, grainF2.inputs[0]); grainF2.inputs[1].default_value = 0.13; grainF2.inputs[2].default_value = 0.935
    grainGrey2 = N.new("ShaderNodeCombineColor"); grainGrey2.location = (-480, 100)
    grainGrey2.mode = "RGB"
    L(grainF2.outputs[0], grainGrey2.inputs[0]); L(grainF2.outputs[0], grainGrey2.inputs[1]); L(grainF2.outputs[0], grainGrey2.inputs[2])
    grainMul2 = N.new("ShaderNodeMix"); grainMul2.data_type = "RGBA"; grainMul2.blend_type = "MULTIPLY"
    grainMul2.location = (-160, 300); grainMul2.inputs["Factor"].default_value = 1.0
    L(col, grainMul2.inputs[6]); L(grainGrey2.outputs["Color"], grainMul2.inputs[7])
    col = grainMul2.outputs[2]

    L(col, bsdf.inputs["Base Color"])

    # ---------------- STAGE 2-5 matte roughness 0.74-0.93 + STAGE 2-3 grain add
    rr2 = N.new("ShaderNodeMapRange"); rr2.location = (1080, -340)
    rr2.inputs["To Min"].default_value = 0.74
    rr2.inputs["To Max"].default_value = 0.93
    L(n2.outputs["Fac"], rr2.inputs["Value"])
    _gsub = math_node("SUBTRACT", (1240, -340)); L(grainG, _gsub.inputs[0]); _gsub.inputs[1].default_value = 0.5
    _gmul = math_node("MULTIPLY", (1360, -340)); L(_gsub.outputs[0], _gmul.inputs[0]); _gmul.inputs[1].default_value = 0.10
    _radd = math_node("ADD", (1480, -340)); L(rr2.outputs["Result"], _radd.inputs[0]); L(_gmul.outputs[0], _radd.inputs[1])
    L(_radd.outputs[0], bsdf.inputs["Roughness"])
    set_in(bsdf, "Specular IOR Level", 0.05)
    set_in(bsdf, "Sheen Weight", 0.10)
    set_in(bsdf, "Sheen Roughness", 0.55)
    if "Sheen Tint" in bsdf.inputs:
        bsdf.inputs["Sheen Tint"].default_value = (*srgb("#dff0da"), 1.0)

    # ---------------- STAGE 2-3 bump from g (replaces old n2 bump; no separate vein bump found)
    bump = N.new("ShaderNodeBump"); bump.location = (1180, -760)
    bump.inputs["Strength"].default_value = 0.26
    bump.inputs["Distance"].default_value = 0.0018
    _bmix = math_node("MULTIPLY_ADD", (1040, -760))
    L(grainG2, _bmix.inputs[0]); _bmix.inputs[1].default_value = 0.45
    _bmix2 = math_node("MULTIPLY", (920, -760)); L(grainG, _bmix2.inputs[0]); _bmix2.inputs[1].default_value = 0.55
    L(_bmix2.outputs[0], _bmix.inputs[2])
    L(_bmix.outputs[0], bump.inputs["Height"])
    L(bump.outputs["Normal"], bsdf.inputs["Normal"])

    # ---------------- STAGE 2-6 fringe alpha (grain-gated cutout past u=1)
    _au = map_range((920, -760), u_, 1.005, 1.035, 0.0, 1.0)
    _ah = math_node("MULTIPLY_ADD", (1080, -760)); L(grainG, _ah.inputs[0]); _ah.inputs[1].default_value = 0.5; _ah.inputs[2].default_value = 0.5
    _am = math_node("MULTIPLY", (1240, -760)); L(_au, _am.inputs[0]); L(_ah.outputs[0], _am.inputs[1])
    _a1 = math_node("SUBTRACT", (1360, -760)); _a1.inputs[0].default_value = 1.0; L(_am.outputs[0], _a1.inputs[1])
    L(_a1.outputs[0], bsdf.inputs["Alpha"])
    m.blend_method = "BLEND"
    if hasattr(m, "use_transparent_shadow"):
        m.use_transparent_shadow = True

    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


def mat_body(name, c_hex, rough=0.85, fur=0.5, scale=60.0, sheen=0.4, bump_dist=0.006):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    N, L = nt.nodes, nt.links.new
    out = N.new("ShaderNodeOutputMaterial"); out.location = (700, 0)
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (400, 0)
    tc = N.new("ShaderNodeTexCoord"); tc.location = (-800, 0)
    n = N.new("ShaderNodeTexNoise"); n.location = (-600, 0)
    n.inputs["Scale"].default_value = scale
    n.inputs["Detail"].default_value = 8.0
    n.inputs["Roughness"].default_value = 0.7
    L(tc.outputs["Object"], n.inputs["Vector"])
    n2 = N.new("ShaderNodeTexNoise"); n2.location = (-600, -260)
    n2.inputs["Scale"].default_value = scale * 5.0
    n2.inputs["Detail"].default_value = 4.0
    L(tc.outputs["Object"], n2.inputs["Vector"])
    ramp = N.new("ShaderNodeValToRGB"); ramp.location = (-380, 160)
    ramp.color_ramp.elements[0].position = 0.34
    ramp.color_ramp.elements[1].position = 0.72
    L(n.outputs["Fac"], ramp.inputs["Fac"])
    c = srgb(c_hex)
    dark = N.new("ShaderNodeRGB"); dark.location = (-380, -60)
    dark.outputs[0].default_value = (c[0] * 0.60, c[1] * 0.60, c[2] * 0.60, 1.0)
    light = N.new("ShaderNodeRGB"); light.location = (-380, -180)
    light.outputs[0].default_value = (min(1, c[0] * 1.18), min(1, c[1] * 1.18), min(1, c[2] * 1.18), 1.0)
    mix = N.new("ShaderNodeMix"); mix.data_type = "RGBA"; mix.location = (60, 120)
    L(ramp.outputs["Color"], mix.inputs["Factor"])
    L(dark.outputs[0], mix.inputs[6])
    L(light.outputs[0], mix.inputs[7])
    L(mix.outputs[2], bsdf.inputs["Base Color"])
    bsdf.inputs["Roughness"].default_value = rough
    set_in(bsdf, "Specular IOR Level", 0.15)
    set_in(bsdf, "Sheen Weight", sheen)
    bump = N.new("ShaderNodeBump"); bump.location = (180, -320)
    bump.inputs["Strength"].default_value = fur
    bump.inputs["Distance"].default_value = bump_dist
    L(n2.outputs["Fac"], bump.inputs["Height"])
    L(bump.outputs["Normal"], bsdf.inputs["Normal"])
    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


def mat_fur_shell(name, c_hex):
    # STAGE 2-7 fur shell: body colour lightened x1.1, alpha from Object noise > 0.45
    # PROTO4 CHANGE 2: finer tufts (noise 40 -> 70), lower cut (0.45 -> 0.40) and
    # a sheen step so the shell reads as fuzz at full view instead of a smooth skin.
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    m.blend_method = "BLEND"
    nt = m.node_tree
    nt.nodes.clear()
    N, L = nt.nodes, nt.links.new
    out = N.new("ShaderNodeOutputMaterial"); out.location = (700, 0)
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (400, 0)
    c = srgb(c_hex)
    bsdf.inputs["Base Color"].default_value = (min(1.0, c[0] * 1.1), min(1.0, c[1] * 1.1), min(1.0, c[2] * 1.1), 1.0)
    bsdf.inputs["Roughness"].default_value = 0.9
    set_in(bsdf, "Sheen Weight", 0.60)
    set_in(bsdf, "Sheen Roughness", 0.40)
    tc = N.new("ShaderNodeTexCoord"); tc.location = (-800, 0)
    nz = N.new("ShaderNodeTexNoise"); nz.location = (-600, 0)
    nz.inputs["Scale"].default_value = 90.0
    nz.inputs["Detail"].default_value = 2.0
    L(tc.outputs["Object"], nz.inputs["Vector"])
    gt = N.new("ShaderNodeMath"); gt.operation = "GREATER_THAN"; gt.location = (-300, -100)
    L(nz.outputs["Fac"], gt.inputs[0]); gt.inputs[1].default_value = 0.36
    L(gt.outputs[0], bsdf.inputs["Alpha"])
    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


# ---------------------------------------------------------------- contact shadow
# Footprint rect in moth-ROOT units (x span, y along body).
# STAGE 3-1: canvas covers root x in [-0.62, 0.62], y in [-0.42, 0.30].
SH_X0, SH_X1 = -0.62, 0.62
SH_Y0, SH_Y1 = -0.42, 0.30


def _scanline_fill(polys, W, H):
    """polys: list of [(x, y), ...] in local units -> HxW float mask (row 0 = SH_Y0)."""
    mask = [[0.0] * W for _ in range(H)]
    for poly in polys:
        for row in range(H):
            y = SH_Y0 + (row + 0.5) / H * (SH_Y1 - SH_Y0)
            xs = []
            n = len(poly)
            for i in range(n):
                x1, y1 = poly[i]
                x2, y2 = poly[(i + 1) % n]
                if (y1 <= y < y2) or (y2 <= y < y1):
                    xs.append(x1 + (y - y1) / (y2 - y1) * (x2 - x1))
            xs.sort()
            for k in range(0, len(xs) - 1, 2):
                c0 = int((xs[k] - SH_X0) / (SH_X1 - SH_X0) * W)
                c1 = int((xs[k + 1] - SH_X0) / (SH_X1 - SH_X0) * W)
                for c in range(max(0, c0), min(W, c1 + 1)):
                    mask[row][c] = 1.0
    return mask


def _box_blur(mask, W, H, radius):
    # horizontal pass (exact window sums; W,H are small)
    tmp = [[0.0] * W for _ in range(H)]
    for r in range(H):
        row = mask[r]
        for c in range(W):
            lo = max(0, c - radius)
            hi = min(W - 1, c + radius)
            tmp[r][c] = sum(row[lo:hi + 1]) / (hi - lo + 1)
    # vertical pass
    out = [[0.0] * W for _ in range(H)]
    for c in range(W):
        for r in range(H):
            lo = max(0, r - radius)
            hi = min(H - 1, r + radius)
            s = 0.0
            for rr in range(lo, hi + 1):
                s += tmp[rr][c]
            out[r][c] = s / (hi - lo + 1)
    return out


def _stroke_mask(polys, W, H, thickness=4.0):
    """Rasterise closed polylines (local units) as a thick line mask.

    Pixels whose centre is within thickness/2 (in px) of any segment = 1.0.
    Same local->pixel mapping as _scanline_fill (row 0 = SH_Y0, pixel (c, r)
    spans [c, c+1] x [r, r+1] in pixel space, centre at +0.5).
    """
    mask = [[0.0] * W for _ in range(H)]
    rad = thickness * 0.5
    sx = W / (SH_X1 - SH_X0)
    sy = H / (SH_Y1 - SH_Y0)
    for poly in polys:
        pts = [((x - SH_X0) * sx, (y - SH_Y0) * sy) for (x, y) in poly]
        n = len(pts)
        for i in range(n):
            x1, y1 = pts[i]
            x2, y2 = pts[(i + 1) % n]
            xlo = max(0, int(min(x1, x2) - rad - 1))
            xhi = min(W - 1, int(max(x1, x2) + rad + 1) + 1)
            ylo = max(0, int(min(y1, y2) - rad - 1))
            yhi = min(H - 1, int(max(y1, y2) + rad + 1) + 1)
            dx = x2 - x1
            dy = y2 - y1
            L2 = dx * dx + dy * dy
            for yy in range(ylo, yhi + 1):
                py = yy + 0.5
                for xx in range(xlo, xhi + 1):
                    px = xx + 0.5
                    if L2 < 1e-12:
                        d = math.hypot(px - x1, py - y1)
                    else:
                        t = ((px - x1) * dx + (py - y1) * dy) / L2
                        t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
                        d = math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))
                    if d <= rad:
                        mask[yy][xx] = 1.0
    return mask


def make_shadow_image(rot_deg=0.0, scale=1.0, tex_name="MothShadowTex", shift=None):
    """Rasterise both fan silhouettes + body ellipse into a soft alpha image.

    STAGE 3-1: rasterised in moth-ROOT coordinates (right fan += (+0.030,
    +0.012), left mirrored fan += (-0.030, +0.012)). Canvas covers root
    x in [-0.62, 0.62], y in [-0.42, 0.30], 512x296. The decal plane
    (build_shadow_decal, SH_* rect) covers exactly the same root rect so
    texture and plane align. Body ellipse is already in root coords.

    PROTO4 CHANGE 3 (contact core): the decal object itself sits at
    +SH_OFF (0.075,-0.075) in WORLD space, so a copy drawn inside the texture
    only lands back under the wing if it is shifted by -R(rot)^-1*SH_OFF/scale
    (the offset is applied after the decal's rotation/scale).  Each decal gets
    its own texture built with its own shift; the copy is drawn as a 3px
    stroke (blur r=1) plus a light fill, over the unchanged soft layers.
    Final alpha = min(0.96, 0.62*sharp + 0.45*wide + 0.95*stroke + 0.45*core
    + 0.60*contact_edge + 0.18*contact_fill).
    """
    W, H = 512, 296
    lex, ley, tex, tey = fan_outline(60)
    fan_polys = []
    for s in (1.0, -1.0):
        ox, oy = s * 0.030, 0.012
        poly = [(s * lex[i] + ox, ley[i] + oy) for i in range(60)]
        poly += [(s * tex[i] + ox, tey[i] + oy) for i in range(59, -1, -1)]
        fan_polys.append(poly)
    # body ellipse (thorax + abdomen + head)
    ep = []
    for i in range(28):
        a = 2 * math.pi * i / 28.0
        ep.append((0.070 * math.cos(a), -0.03 + 0.27 * math.sin(a)))
    polys = fan_polys + [ep]
    mask = _scanline_fill(polys, W, H)

    # contact copy shift (texture-local) that lands under the wing after the
    # decal's rotate(rot_deg) / scale / +SH_OFF transform
    if shift is None:
        th = math.radians(rot_deg)
        ox, oy = 0.075, -0.075
        vx = ox * math.cos(th) + oy * math.sin(th)
        vy = -ox * math.sin(th) + oy * math.cos(th)
        shift = (-vx / max(scale, 1e-6), -vy / max(scale, 1e-6))
    sx, sy = shift
    ct_polys = [[(x + sx, y + sy) for (x, y) in p] for p in polys]
    ct_mask = _scanline_fill(ct_polys, W, H)
    ct_edge = _stroke_mask(ct_polys, W, H, thickness=3.0)
    ct_b = _box_blur(ct_edge, W, H, 1)
    ct_fill = _box_blur(ct_mask, W, H, 3)

    sharp = _box_blur(mask, W, H, 2)
    wide = _box_blur(_box_blur(mask, W, H, 14), W, H, 14)
    stroke = _stroke_mask(fan_polys, W, H, thickness=5.0)
    stroke_b = _box_blur(stroke, W, H, 2)
    inner = _box_blur(mask, W, H, 2)
    px = []
    for r in range(H):
        for c in range(W):
            # tight dark core ~2 px inside the outline
            core = (inner[r][c] - 0.50) / 0.32
            core = 0.0 if core < 0.0 else (1.0 if core > 1.0 else core)
            a = (0.62 * sharp[r][c] + 0.45 * wide[r][c] + 0.85 * stroke_b[r][c]
                 + 0.22 * core + 0.22 * ct_b[r][c] + 0.08 * ct_fill[r][c])
            a = min(0.90, a)
            px.extend((0.0, 0.0, 0.0, a))
    img = bpy.data.images.new(tex_name, W, H, alpha=True)
    img.pixels.foreach_set(px)
    img.update()
    try:
        img.pack()
    except Exception:
        pass
    print("shadow tex %s contact shift=(%.4f, %.4f) rot=%.0f scale=%.3f" % (
        tex_name, sx, sy, rot_deg, scale))
    return img
    img = bpy.data.images.new("MothShadowTex", W, H, alpha=True)
    img.pixels.foreach_set(px)
    img.update()
    try:
        img.pack()
    except Exception:
        pass
    return img


def mat_shadow_silhouette(img):
    # Principled near-black with Alpha = texture alpha x opacity node
    # (same proven pattern as proto1's ShadowDecal material).
    m = bpy.data.materials.new("ShadowSilhouette")
    m.use_nodes = True
    m.blend_method = "BLEND"
    if hasattr(m, "surface_render_method"):
        m.surface_render_method = "BLENDED"
    if hasattr(m, "use_transparent_shadow"):
        m.use_transparent_shadow = True
    nt = m.node_tree
    nt.nodes.clear()
    N, L = nt.nodes, nt.links.new
    out = N.new("ShaderNodeOutputMaterial"); out.location = (700, 0)
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (420, 0)
    bsdf.inputs["Base Color"].default_value = (0.002, 0.003, 0.002, 1)
    bsdf.inputs["Roughness"].default_value = 1.0
    set_in(bsdf, "Specular IOR Level", 0.0)
    tex = N.new("ShaderNodeTexImage"); tex.location = (-200, 0)
    tex.image = img
    op = N.new("ShaderNodeValue"); op.location = (0, -260); op.label = "opacity"; op.name = "opacity"
    op.outputs[0].default_value = 1.0
    mul = N.new("ShaderNodeMath"); mul.operation = "MULTIPLY"; mul.location = (160, -120)
    L(tex.outputs["Alpha"], mul.inputs[0])
    L(op.outputs[0], mul.inputs[1])
    L(mul.outputs[0], bsdf.inputs["Alpha"])
    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


def build_shadow_decal(name, img, mat):
    verts = [(SH_X0, SH_Y0, 0.0012), (SH_X1, SH_Y0, 0.0012),
             (SH_X1, SH_Y1, 0.0012), (SH_X0, SH_Y1, 0.0012)]
    uv = [(0, 0), (1, 0), (1, 1), (0, 1)]
    ob = new_mesh_obj(name, verts, [(0, 1, 2, 3)], bpy.context.scene.collection, uv=uv, smooth=False)
    ob.data.materials.append(mat)
    ob.visible_shadow = False
    return ob


# ---------------------------------------------------------------- wing geometry
# Outline traced from the reference photo. t: 0 = wing base, 1 = apex.
# LE = costa (leading edge), TE = outer margin continuing into the inner
# margin of the hindwing. One continuous fan per side.
FAN_LE = [(0.00, 0.002, 0.052), (0.055, 0.003, 0.070), (0.11, 0.005, 0.086),
          (0.16, 0.006, 0.095), (0.286, 0.080, 0.130), (0.454, 0.180, 0.148),
          (0.622, 0.280, 0.134), (0.790, 0.375, 0.098), (0.916, 0.440, 0.048),
          (0.975, 0.492, 0.014), (1.00, 0.505, -0.002)]
FAN_TE = [(0.00, -0.030, -0.312), (0.055, -0.029, -0.305), (0.11, -0.024, -0.298),
          (0.16, 0.000, -0.295), (0.244, 0.030, -0.300), (0.345, 0.078, -0.318),
          (0.454, 0.132, -0.326), (0.580, 0.196, -0.352), (0.706, 0.235, -0.284),
          (0.815, 0.288, -0.215), (0.916, 0.380, -0.150), (0.975, 0.452, -0.076),
          (0.990, 0.497, -0.026), (1.00, 0.500, -0.012)]


def _sstep(a, b, x):
    f = min(1.0, max(0.0, (x - a) / (b - a)))
    return f * f * (3.0 - 2.0 * f)


def make_fan_wing(side, mir, name, coll, mat):
    sign = mir * side
    Nt, Nu = 64, 22
    # STAGE 2-6 fringe: 3 extra columns past the trailing edge
    FRINGE_U = (1.012, 1.024, 1.035)
    NCOL = (Nu + 1) + len(FRINGE_U)
    lex, ley, tex, tey = fan_outline(Nt + 1)
    verts, uvs, faces = [], [], []
    for i in range(Nt + 1):
        t = i / float(Nt)
        ty = tey[i]
        tx = tex[i]
        for j in range(NCOL):
            if j <= Nu:
                u = j / float(Nu)
                fringe_step = 0
            else:
                u = FRINGE_U[j - (Nu + 1)]
                fringe_step = j - Nu
            lx = lex[i] + (tx - lex[i]) * u
            ly = ley[i] + (ty - ley[i]) * u
            x = lx
            y = ly
            xr = side * (lx + 0.030)
            z = 0.008 * math.sin(math.pi * u) * (0.35 + 0.65 * t) \
                + 0.006 * math.exp(-(xr / 0.055) ** 2) * (1 - _sstep(0.02, 0.16, t)) \
                + 0.0035 * math.sin(2 * math.pi * (3.1 * lx + 2.3 * ly) + 1.7) * (0.3 + 0.7 * u) \
                - 0.010 * (u ** 6) * _sstep(0.25, 0.75, t) \
                * (0.55 + 0.45 * math.sin(2 * math.pi * 2.6 * t + 0.8)) \
                - 0.006 * t * t
            if fringe_step:
                z -= 0.0012 * fringe_step
            verts.append(Vector((sign * x, y, z)))
            uvs.append((t, u))
    for i in range(Nt):
        for j in range(NCOL - 1):
            a = i * NCOL + j
            b = a + 1
            c = a + NCOL
            d = c + 1
            if side > 0:
                faces.append([a, b, d, c])
            else:
                faces.append([a, c, d, b])
    # STAGE 3b-1 centreline half-patch: strip from the fan's t=0 chord to the
    # centreline, same mesh/material so it flaps with the wing.
    # Grid coords (pre-mirror): LE0/TE0 from FAN_LE/FAN_TE first points;
    # centreline grid x = -0.030 (mesh x_c = sign*-0.030; brief lists
    # -side*0.030 which equals this for mir=1 / the single-moth stills).
    _LE0x, _LE0y = 0.002, 0.052
    _TE0x, _TE0y = -0.030, -0.312
    _XC_G = -0.030
    _TP, _LIP = 0.02, 0.0015
    _patch_base = len(verts)
    for _k in range(11):
        _s = _k / 10.0
        _px = _LE0x + _s * (_TE0x - _LE0x)
        _py = _LE0y + _s * (_TE0y - _LE0y)
        _pts = [(_px, _py),
                ((_px + _XC_G) * 0.5, _py),
                (_XC_G, _py)]
        _uv = (0.03, 0.05 + 0.9 * _s)
        for (_lx, _ly) in _pts:
            _xr = side * (_lx + 0.030)
            _z = 0.008 * math.sin(math.pi * _s) * (0.35 + 0.65 * _TP) \
                + 0.006 * math.exp(-(_xr / 0.055) ** 2) * (1 - _sstep(0.02, 0.16, _TP)) \
                + 0.0035 * math.sin(2 * math.pi * (3.1 * _lx + 2.3 * _ly) + 1.7) * (0.3 + 0.7 * _s) \
                - 0.010 * (_s ** 6) * _sstep(0.25, 0.75, _TP) \
                * (0.55 + 0.45 * math.sin(2 * math.pi * 2.6 * _TP + 0.8)) \
                - 0.006 * _TP * _TP + _LIP
            verts.append(Vector((sign * _lx, _ly, _z)))
            uvs.append(_uv)
    for _k in range(10):
        _Ak = _patch_base + _k * 3
        _Bk = _Ak + 1
        _Ck = _Ak + 2
        _An = _Ak + 3
        _Bn = _Ak + 4
        _Cn = _Ak + 5
        if side > 0:
            faces.append([_Ak, _Bk, _Bn, _An])
            faces.append([_Bk, _Ck, _Cn, _Bn])
        else:
            faces.append([_Ak, _An, _Bn, _Bk])
            faces.append([_Bk, _Bn, _Cn, _Ck])
    ob = new_mesh_obj(name, verts, faces, coll, uv=uvs)
    ob.data.materials.append(mat)
    solid = ob.modifiers.new("Solidify", "SOLIDIFY")
    solid.thickness = 0.0016
    solid.offset = 0.0
    solid.use_even_offset = True
    return ob


# ---------------------------------------------------------------- moth
def build_moth(name, mirror, coll, mats):
    root = bpy.data.objects.new(name, None)
    root.empty_display_type = "PLAIN_AXES"
    root.empty_display_size = 0.08
    coll.objects.link(root)

    wing_objs = []
    for side in (1, -1):
        wname = "%s_wing_fan_%s" % (name, "R" if side > 0 else "L")
        w = make_fan_wing(side, mirror, wname, coll, mats["wing"])
        w.location = Vector((mirror * side * 0.030, 0.012, 0.015))
        w.parent = root
        wing_objs.append((w, "fan", side))

    body_parts = []
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=24, v_segments=16, radius=1.0)
    # PROTO4 CHANGE 2: lower + broaden the thorax and break the smooth oval
    # outline with a small seeded displacement (lumpy piliform mass, not an egg).
    # radii (0.062,0.080,0.050) z-centre 0.036 -> (0.072,0.086,0.042) z 0.030
    TH_SEED = 4.17
    TH_AMP = 0.13

    def _th_lump(x, y, z):
        s = TH_SEED
        return (0.60 * math.sin(17.3 * x + 1.1 + s) * math.sin(11.7 * y + 2.3 + s)
                * math.sin(19.1 * z + 0.4 + s)
                + 0.28 * math.sin(29.7 * x + 0.7 + 2 * s) * math.sin(19.3 * y + 1.9 + s)
                * math.sin(27.3 * z + 2.7 + 3 * s)
                + 0.12 * math.sin(41.1 * x + 2.1 + 3 * s) * math.sin(13.1 * y + 0.3 + 2 * s)
                * math.sin(37.3 * z + 1.7 + s))

    for v in bm.verts:
        ux, uy, uz = v.co.x, v.co.y, v.co.z
        d = 1.0 + TH_AMP * _th_lump(ux, uy, uz)
        v.co.x = ux * 0.072 * d
        v.co.y = uy * 0.086 * d + 0.028
        v.co.z = uz * 0.042 * d + 0.030
    me = bpy.data.meshes.new(name + "_thorax")
    bm.to_mesh(me)
    bm.free()
    for p in me.polygons:
        p.use_smooth = True
    th = bpy.data.objects.new(name + "_thorax", me)
    coll.objects.link(th)
    th.parent = root
    th.data.materials.append(mats["body"])
    body_parts.append(th)

    prof = []
    nseg = 10
    for i in range(nseg + 1):
        t = i / float(nseg)
        y = -0.020 - 0.215 * t
        r = (0.036 * (1.0 - t) ** 0.60 + 0.005) * 0.85
        r *= 1.0 + 0.05 * math.sin(2 * math.pi * 4.0 * t)
        prof.append((y, r))
    ring = 18
    verts, faces = [], []
    for i, (y, r) in enumerate(prof):
        t = i / float(nseg)
        z = 0.006 - 0.006 * t * t
        for k in range(ring):
            a = 2 * math.pi * k / ring
            verts.append((mirror * r * math.cos(a) * 0.80, y, z + r * math.sin(a) * 0.35))
    for i in range(len(prof) - 1):
        for k in range(ring):
            k2 = (k + 1) % ring
            a = i * ring + k
            b = i * ring + k2
            c = (i + 1) * ring + k2
            d = (i + 1) * ring + k
            faces.append([a, b, c, d])
    faces.append([k for k in range(ring - 1, -1, -1)])
    base = (len(prof) - 1) * ring
    faces.append([base + k for k in range(ring)])
    ab = new_mesh_obj(name + "_abdomen", verts, faces, coll)
    ab.data.materials.append(mats["abdomen"])
    ab.parent = root
    body_parts.append(ab)

    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=20, v_segments=12, radius=1.0)
    for v in bm.verts:
        v.co.x *= 0.036
        v.co.y *= 0.034
        v.co.z *= 0.032
        v.co.z += 0.042
        v.co.y += 0.100
    me = bpy.data.meshes.new(name + "_head")
    bm.to_mesh(me)
    bm.free()
    for p in me.polygons:
        p.use_smooth = True
    hd = bpy.data.objects.new(name + "_head", me)
    coll.objects.link(hd)
    hd.data.materials.append(mats["body"])
    hd.parent = root
    body_parts.append(hd)

    # STAGE 2-7 fur shells: duplicate ONLY thorax and head, scale 1.13 about own centres
    def _fur_shell(src_ob, shell_name):
        src_me = src_ob.data
        vs = [v.co.copy() for v in src_me.vertices]
        cx = sum(v.x for v in vs) / len(vs)
        cy = sum(v.y for v in vs) / len(vs)
        cz = sum(v.z for v in vs) / len(vs)
        sverts = [(cx + (v.x - cx) * 1.18, cy + (v.y - cy) * 1.18, cz + (v.z - cz) * 1.18) for v in vs]
        sfaces = [list(p.vertices) for p in src_me.polygons]
        sh = new_mesh_obj(shell_name, sverts, sfaces, coll)
        sh.data.materials.append(mats["fur"])
        sh.parent = root
        try:
            sh.visible_shadow = False
        except Exception:
            pass
        body_parts.append(sh)
        return sh
    _fur_shell(th, name + "_thorax_fur")
    _fur_shell(hd, name + "_head_fur")

    for side in (1, -1):
        bm = bmesh.new()
        bmesh.ops.create_uvsphere(bm, u_segments=12, v_segments=8, radius=1.0)
        for v in bm.verts:
            v.co *= 0.014
            v.co.x += mirror * side * 0.023
            v.co.y += 0.116
            v.co.z += 0.052
        me = bpy.data.meshes.new(name + "_eye")
        bm.to_mesh(me)
        bm.free()
        for p in me.polygons:
            p.use_smooth = True
        ey = bpy.data.objects.new(name + "_eye_%s" % ("R" if side > 0 else "L"), me)
        coll.objects.link(ey)
        ey.data.materials.append(mats["eye"])
        ey.parent = root
        body_parts.append(ey)

    leg_objs = []
    for side in (1, -1):
        for k, (y0, dx, dy) in enumerate([(0.076, 0.100, 0.140),
                                          (0.030, 0.115, 0.015),
                                          (-0.016, 0.108, -0.105)]):
            p0 = Vector((mirror * side * 0.026, y0, 0.020))
            p1 = Vector((mirror * side * (0.026 + dx * 0.55), y0 + dy * 0.45, 0.010))
            p2 = Vector((mirror * side * (0.026 + dx), y0 + dy, 0.001))
            leg = tube("%s_leg_%s%d" % (name, "R" if side > 0 else "L", k),
                       [p0, p1, p2], [0.0058, 0.0040, 0.0024], coll, segs=6)
            leg.data.materials.append(mats["leg"])
            leg.parent = root
            leg_objs.append(leg)

    for side in (1, -1):
        pts = []
        for i in range(9):
            t = i / 8.0
            pts.append(Vector((
                mirror * side * (0.014 + 0.055 * t),
                0.124 + 0.164 * t - 0.030 * t * t,
                0.058 + 0.055 * t - 0.010 * t * t)))
        radii = [0.0036 * 0.8 * (1.0 - 0.55 * (i / 8.0)) for i in range(9)]
        an = tube("%s_antenna_%s" % (name, "R" if side > 0 else "L"), pts, radii, coll, segs=6)
        an.data.materials.append(mats["leg"])
        an.parent = root
        leg_objs.append(an)

    return dict(root=root, wings=wing_objs, body=body_parts, legs=leg_objs)


# ---------------------------------------------------------------- build
def main():
    wipe()
    scene = bpy.context.scene
    coll = scene.collection

    scene.render.engine = "BLENDER_EEVEE_NEXT"
    scene.render.resolution_x = RES_X
    scene.render.resolution_y = RES_Y
    scene.render.resolution_percentage = 100
    scene.render.fps = FPS
    scene.frame_start = F_START
    scene.frame_end = F_END
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    scene.render.film_transparent = False
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGB"
    scene.render.image_settings.color_depth = "8"
    scene.render.image_settings.compression = 15
    ee = scene.eevee
    ee.taa_render_samples = 24 if QUICK else 64
    ee.use_gtao = True
    if hasattr(ee, "gtao_distance"):
        ee.gtao_distance = 0.10
    if hasattr(ee, "gtao_factor"):
        ee.gtao_factor = 0.9
    ee.use_shadows = True
    scene.render.use_motion_blur = True
    scene.render.motion_blur_shutter = 0.5

    mats = dict(
        wing=mat_wing(),
        body=mat_body("BodyMat", "#6f8f66", rough=0.88, fur=0.95, scale=70.0, sheen=0.85, bump_dist=0.009),
        abdomen=mat_body("AbdomenMat", "#6d8a64", rough=0.90, fur=0.55, scale=85.0),
        leg=mat_body("LegMat", "#4a4030", rough=0.66, fur=0.20, scale=120.0),
        eye=mat_body("EyeMat", "#191a12", rough=0.35, fur=0.0, scale=10.0),
        bark=mat_background(),
        fur=mat_fur_shell("FurShellMat", "#6f8f66"),
    )

    bw, bh = VIEW_W, VIEW_H
    verts = [(-bw / 2, -bh / 2, 0.0), (bw / 2, -bh / 2, 0.0),
             (bw / 2, bh / 2, 0.0), (-bw / 2, bh / 2, 0.0)]
    uv = [(0, 0), (1, 0), (1, 1), (0, 1)]
    me = bpy.data.meshes.new("BarkPlateMesh")
    me.from_pydata(verts, [], [(0, 1, 2, 3)])
    uvl = me.uv_layers.new(name="UVMap")
    for li, loop in enumerate(me.loops):
        uvl.data[li].uv = uv[loop.vertex_index]
    me.update()
    plate = bpy.data.objects.new("BarkPlate", me)
    coll.objects.link(plate)
    plate.data.materials.append(mats["bark"])
    plate.visible_shadow = False

    S = 1.63
    SCALE_L = 0.65 * S
    SCALE_R = 0.71 * S
    moth_l = build_moth("Moth_L", -1, coll, mats)
    moth_l["root"].location = Vector((fx(0.300), fy(0.520), 0.0))
    moth_l["root"].rotation_euler = (0.0, 0.0, math.radians(50.0))
    moth_l["root"].scale = (SCALE_L,) * 3
    REST_LOC_L = Vector(moth_l["root"].location)

    moth_r = build_moth("Moth_R", 1, coll, mats)
    moth_r["root"].location = Vector((fx(0.700), fy(0.490), 0.0))
    moth_r["root"].rotation_euler = (0.0, 0.0, math.radians(-60.0))
    moth_r["root"].scale = (SCALE_R,) * 3
    REST_LOC_R = Vector(moth_r["root"].location)

    # silhouette contact shadows, aligned with each moth's placement;
    # offset down-right in frame space (sun comes from the upper left).
    # PROTO4 CHANGE 3: one texture per decal so the contact-copy shift can
    # cancel that decal's world offset after its own rotation/scale.
    SH_OFF = Vector((0.075, -0.075, 0.0))
    sh_img_l = make_shadow_image(50.0, SCALE_L, "MothShadowTexL")
    sh_img_r = make_shadow_image(-60.0, SCALE_R, "MothShadowTexR")
    try:
        sh_img_r.filepath_raw = os.path.join(HERE, "_dbg", "shadow_tex.png")
        sh_img_r.file_format = "PNG"
        sh_img_r.save()
        print("shadow tex ->", sh_img_r.filepath_raw)
    except Exception as exc:
        print("shadow tex save failed:", exc)
    sh_mat_l = mat_shadow_silhouette(sh_img_l)
    sh_mat_r = mat_shadow_silhouette(sh_img_r)
    sh_l = build_shadow_decal("Shadow_L", sh_img_l, sh_mat_l)
    sh_l.location = REST_LOC_L + SH_OFF
    sh_l.rotation_euler = (0.0, 0.0, math.radians(50.0))
    sh_l.scale = (SCALE_L,) * 3
    sh_r = build_shadow_decal("Shadow_R", sh_img_r, sh_mat_r)
    sh_r.location = REST_LOC_R + SH_OFF
    sh_r.rotation_euler = (0.0, 0.0, math.radians(-60.0))
    sh_r.scale = (SCALE_R,) * 3

    cam_data = bpy.data.cameras.new("Cam")
    cam_data.lens = 85.0
    cam_data.sensor_width = 36.0
    cam_data.sensor_fit = "HORIZONTAL"
    cam = bpy.data.objects.new("Cam", cam_data)
    coll.objects.link(cam)
    dist = VIEW_W * cam_data.lens / cam_data.sensor_width
    cam.location = (0.0, 0.0, dist)
    cam.rotation_euler = (0.0, 0.0, 0.0)
    scene.camera = cam

    sun_data = bpy.data.lights.new("Sun", type="SUN")
    sun_data.energy = SUN_STRENGTH
    sun_data.angle = math.radians(12.0)
    sun_data.color = (1.0, 0.988, 0.96)
    sun = bpy.data.objects.new("Sun", sun_data)
    coll.objects.link(sun)
    sun.rotation_euler = (math.radians(-38.0), 0.0, math.radians(40.0))

    world = bpy.data.worlds.new("World")
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    bg.inputs["Color"].default_value = (0.66, 0.70, 0.64, 1.0)
    bg.inputs["Strength"].default_value = AMBIENT
    scene.world = world

    # ================================================ ANIMATION
    # PROTO4 STAGE C - choreography v2 (rules from MOTION-DOSSIER.md).
    #
    # Every animated value below is a PURE FUNCTION OF THE FRAME NUMBER: no
    # wall clock, no unseeded randomness, no frame-rate dependent sampling.
    # The irregular intervals were measured off real footage, but their exact
    # values are literals in MOTION_TABLE (SEED is documentation only).  All
    # keyframes are written once per frame and switched to LINEAR, so what is
    # keyed is exactly what is evaluated (motion blur only interpolates
    # between adjacent frames).
    #
    # The log written to _dbg/motion-log.json comes from these same functions
    # (see tools/motion_verify.py) - it is not a re-derivation.
    fps = FPS
    SEED = 20261110
    BASE_HZ = 9.5                      # proto3 flap_hz, kept as the base period
    BASE_T = fps / BASE_HZ             # 2.52632 frames per beat at 24 fps

    # ---------------------------------------------------- seeded event table
    # measured  = taken from the dossier's measurements
    # artistic  = explicitly invented (dossier "純粋な創作"), kept small
    MOTION_TABLE = {
        "seed": SEED,
        "meaning": "Every irregular value here is a literal; SEED is only the "
                   "label for the table (no runtime randomness anywhere).",
        "left": {
            "rest_deg": 3.0,                     # measured (proto2/3 rest pose)
            "settle": {"frame": 6, "dip_deg": 1.5, "kind": "artistic"},
            "antenna_twitch": {"peaks": [[3, 2.0]], "settle_frames": 2,
                               "kind": "artistic"},
            "raise": {"start": 8, "end": 15, "to_deg": 28.0,
                      "kind": "measured (raise phase before liftoff)"},
            "hold": {"from": 16, "to": 17,
                     "kind": "measured (short planted hold)"},
            "pitch_start": 17,                   # measured (small fwd pitch, planted)
            "liftoff": 18,                       # measured
            "beat_period_factors": [1.00, 0.90, 1.12, 0.94],
            "stroke_top_deg": 32.0,
            "twist_deg": 7.0,
            "first_downstroke_lead_frames": 1.0,  # artistic (dossier)
        },
        "right": {
            "breath_keys": [[0, 3.0], [6, 6.5], [14, 3.0], [23, 3.0],
                            [28, 6.5], [37, 3.0], [47, 3.0]],
            "antenna_both": {"peaks": [[5, 2.5]], "settle_frames": 3,
                             "kind": "artistic"},
            "wing_flick": {"peaks": [[15, 2.5]], "rise": 2, "fall": 2,
                           "kind": "artistic"},
            "antenna_right_only": {"peaks": [[22, 3.0]], "span": 4,
                                   "kind": "artistic"},
            "wing_raise_cycle": {"keys": [[30, 0.0], [34, 9.0], [36, 9.0],
                                          [39, 0.0]],
                                 "kind": "observed behaviour class (Hypagyrtis), "
                                         "staging artistic"},
        },
    }

    # ------------------------------------------------------------ helpers
    def _clamp01(x):
        return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)

    def _sstep(x):
        x = _clamp01(x)
        return x * x * (3.0 - 2.0 * x)

    def _piecewise(keys, f):
        """Smoothstep through explicit (frame, value) keys, clamped at the ends."""
        if f <= keys[0][0]:
            return float(keys[0][1])
        for i in range(len(keys) - 1):
            f0, v0 = keys[i]
            f1, v1 = keys[i + 1]
            if f0 <= f <= f1:
                return v0 + (v1 - v0) * _sstep((f - f0) / float(f1 - f0))
        return float(keys[-1][1])

    def _pchip(knots, flat_start=False, flat_end=False):
        """Monotone cubic (Fritsch-Carlson) through (frame, value) knots.

        Used for the exit PATH so positions travel smoothly through proto3's
        measured envelope without overshooting it.
        """
        xs = [float(k[0]) for k in knots]
        ys = [float(k[1]) for k in knots]
        n = len(xs)
        h = [xs[i + 1] - xs[i] for i in range(n - 1)]
        d = [(ys[i + 1] - ys[i]) / h[i] for i in range(n - 1)]
        m = [0.0] * n
        if n == 2:
            m = [d[0], d[0]]
        else:
            m[0] = d[0]
            m[-1] = d[-1]
            for i in range(1, n - 1):
                if d[i - 1] * d[i] <= 0.0:
                    m[i] = 0.0
                else:
                    w1 = 2.0 * h[i] + h[i - 1]
                    w2 = h[i] + 2.0 * h[i - 1]
                    m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
        if flat_start:
            m[0] = 0.0
        if flat_end:
            m[-1] = 0.0

        def ev(t):
            if t <= xs[0]:
                return ys[0]
            if t >= xs[-1]:
                return ys[-1]
            for i in range(n - 1):
                if xs[i] <= t <= xs[i + 1]:
                    s = (t - xs[i]) / h[i]
                    s2 = s * s
                    s3 = s2 * s
                    return ((2 * s3 - 3 * s2 + 1) * ys[i]
                            + (s3 - 2 * s2 + s) * h[i] * m[i]
                            + (-2 * s3 + 3 * s2) * ys[i + 1]
                            + (s3 - s2) * h[i] * m[i + 1])
            return ys[-1]
        return ev

    def flap_sign(mir, side):
        return -(mir * side)

    def set_wing(w, kind, side, mir, up_deg):
        s = flap_sign(mir, side)
        w.rotation_euler[1] = math.radians(s * up_deg)

    # ------------------------------------------------- irregular beat grid
    # "3-4 wingbeats with SEEDED irregular periods": the listed factors are
    # applied in order to the base period, starting at the liftoff frame.
    LIFT = MOTION_TABLE["left"]["liftoff"]
    BEAT_T = [float(LIFT)]
    for _k in MOTION_TABLE["left"]["beat_period_factors"]:
        BEAT_T.append(BEAT_T[-1] + BASE_T * _k)
    while BEAT_T[-1] < F_END + 4.0:
        BEAT_T.append(BEAT_T[-1] + BASE_T)      # keep flapping (off-screen)

    def stroke_phase(f):
        """Left wing stroke phase: integer = beat index, fraction = beat position."""
        if f <= BEAT_T[0]:
            return 0.0
        for i in range(len(BEAT_T) - 1):
            if BEAT_T[i] <= f <= BEAT_T[i + 1]:
                return i + (f - BEAT_T[i]) / (BEAT_T[i + 1] - BEAT_T[i])
        i = len(BEAT_T) - 2
        return i + (f - BEAT_T[i]) / (BEAT_T[i + 1] - BEAT_T[i])

    def lag_phase(f):
        """Phase of the lagging wing: the first downstroke starts 1 frame late
        and is compressed into the rest of the same beat, so the two wings are
        back in sync from the second beat onwards."""
        if f >= BEAT_T[1]:
            return stroke_phase(f)
        if f <= LIFT:
            return 0.0
        return (f - (LIFT + 1.0)) / (BEAT_T[1] - (LIFT + 1.0))

    # ---------------------------------------------------- LEFT moth posture
    L_UP = float(MOTION_TABLE["left"]["raise"]["to_deg"])
    STROKE_TOP = float(MOTION_TABLE["left"]["stroke_top_deg"])
    TWIST = float(MOTION_TABLE["left"]["twist_deg"])

    def left_up(f):
        """Left moth wing dihedral (deg) for the LEADING wing."""
        if f >= LIFT:
            return STROKE_TOP * math.cos(2.0 * math.pi * stroke_phase(f))
        if f <= 7:
            # rest 3 deg, plus the f6 settle dip (down 1 frame, back next)
            return 3.0 - (1.5 if f == 6 else 0.0)
        if f <= 15:
            return 3.0 + (L_UP - 3.0) * _sstep((f - 8) / 7.0)
        return L_UP                                   # f16-f17 planted hold

    def left_up_lag(f):
        """Same, for the wing that leads by 1 frame less (first stroke only)."""
        if f >= BEAT_T[1] or f < LIFT:
            return left_up(f)
        return L_UP * math.cos(2.0 * math.pi * lag_phase(f))

    def left_twist(f, phase):
        if f < LIFT:
            return 0.0
        return TWIST * math.cos(2.0 * math.pi * phase) * _clamp01((f - 17.0) / 3.0)

    # --------------------------------------------- LEFT moth exit envelope
    # Position keys 24 / 28 / 36 / 47 are proto3's own path values (captured
    # from the stage-B scene with tools/exitcap.py); f17-f18 hold at rest, so
    # the same exit envelope (and the same f24 framing) is preserved while the
    # wing/body behaviour is re-timed.  0..17 never translate.
    _rx, _ry = REST_LOC_L.x, REST_LOC_L.y
    PATH_X = _pchip([(17, _rx), (18, _rx), (24, -1.09158),
                     (28, -1.365), (36, -2.220), (47, -3.180)], flat_start=True)
    PATH_Y = _pchip([(17, _ry), (18, _ry), (24, 0.27882),
                     (28, 0.43269), (36, 0.84808), (47, 1.31538)], flat_start=True)
    PATH_Z = _pchip([(17, 0.0), (18, 0.015), (24, 0.60945),
                     (28, 0.760), (36, 1.020), (47, 1.220)], flat_start=True)
    PITCH_EV = _pchip([(16, 0.0), (17, 0.5), (18, 2.5), (24, 20.45),
                       (28, 22.0), (36, 20.0), (47, 18.0)])
    ROLL_EV = _pchip([(17, 0.0), (18, 0.0), (24, -26.22), (28, -30.0),
                      (36, -34.0), (47, -34.0)])
    YAW_EV = _pchip([(17, 0.0), (18, 0.0), (24, 13.255), (28, 16.0),
                     (36, 20.0), (47, 22.0)])

    # ---------------------------------------------- RIGHT moth (perched) --
    BREATH = MOTION_TABLE["right"]["breath_keys"]

    def right_wing_up(f):
        """Slow breathing (proto3's 3 -> 6.5 deg drift) + the seeded micro events.

        No periodic perched rhythm is added: the extra events are single,
        irregularly spaced and short (dossier: periodicity NOT measured).
        """
        v = _piecewise(BREATH, f)                       # breathing, continues
        v += _piecewise([[12, 0.0], [15, 2.5], [17, 0.0]], f)        # f13 flick
        v += _piecewise([[30, 0.0], [34, 9.0], [36, 9.0], [39, 0.0]], f)  # f31
        return v                                              # (0 deg elsewhere)

    def right_ant_both(f):
        return _piecewise([[4, 0.0], [5, 2.5], [8, 0.0]], f)

    def right_ant_right(f):
        return _piecewise([[20, 0.0], [22, 3.0], [24, 0.0]], f)

    def left_ant(f):
        return _piecewise([[2, 0.0], [3, 2.0], [5, 0.0]], f)

    # ------------------------------------------------------ antenna rig ----
    # The antennae are single tubes with their object origin at the moth
    # root origin, so a rotation about the object origin would swing the whole
    # tube.  Rotating about the tube's BASE POINT instead (pivot compensation
    # on the object location) is a pure animation transform - no geometry,
    # material or rig change.
    def _ant_pivot(ob):
        c = Vector((0.0, 0.0, 0.0))
        for i in range(6):                    # first tube ring == path[0]
            c += ob.data.vertices[i].co
        return c / 6.0

    def _set_antenna(ob, pivot, deg):
        a = math.radians(deg)
        ca, sa = math.cos(a), math.sin(a)
        ry = pivot.y * ca - pivot.z * sa
        rz = pivot.y * sa + pivot.z * ca
        ob.rotation_euler = (a, 0.0, 0.0)
        ob.location = Vector((0.0, pivot.y - ry, pivot.z - rz))

    ANT_L = [(o, _ant_pivot(o)) for o in moth_l["legs"] if "_antenna_" in o.name]
    ANT_R = [(o, _ant_pivot(o)) for o in moth_r["legs"] if "_antenna_" in o.name]
    L = moth_l["root"]
    R = moth_r["root"]

    # ------------------------------------------------- per-frame bake ------
    log_rows = []
    for f in range(F_START, F_END + 1):
        # ---- left moth
        up_lead = left_up(f)
        up_lag = left_up_lag(f)
        p_lead = stroke_phase(f) if f >= LIFT else 0.0
        p_lag = lag_phase(f) if f >= LIFT else 0.0
        tw_lead = left_twist(f, p_lead)
        tw_lag = left_twist(f, p_lag)
        ant_l = left_ant(f)
        loc = Vector((PATH_X(f), PATH_Y(f), PATH_Z(f)))
        pit, rol, yaw = PITCH_EV(f), ROLL_EV(f), YAW_EV(f)
        for w, kind, side in moth_l["wings"]:
            lead = w.name.endswith("_L")            # this wing starts the stroke
            set_wing(w, kind, side, -1, up_lead if lead else up_lag)
            w.rotation_euler[0] = math.radians(tw_lead if lead else tw_lag)
            w.keyframe_insert("rotation_euler", frame=f)
        L.location = loc
        L.keyframe_insert("location", frame=f)
        L.rotation_euler = (math.radians(pit), math.radians(rol),
                            math.radians(50.0 + yaw))
        L.keyframe_insert("rotation_euler", frame=f)
        for ob, piv in ANT_L:
            _set_antenna(ob, piv, ant_l)
            ob.keyframe_insert("rotation_euler", frame=f)
            ob.keyframe_insert("location", frame=f)

        # ---- right moth (never translates; feet planted)
        up_r = right_wing_up(f)
        ab = right_ant_both(f)
        ar = right_ant_right(f)
        for w, kind, side in moth_r["wings"]:
            set_wing(w, kind, side, 1, up_r)
            w.rotation_euler[0] = 0.0
            w.keyframe_insert("rotation_euler", frame=f)
        for ob, piv in ANT_R:
            _set_antenna(ob, piv, ar if ob.name.endswith("_R") else ab)
            ob.keyframe_insert("rotation_euler", frame=f)
            ob.keyframe_insert("location", frame=f)

        log_rows.append(dict(
            f=f,
            L_loc=[round(v, 6) for v in loc],
            L_pitch=round(pit, 4), L_roll=round(rol, 4), L_yaw=round(yaw, 4),
            L_wing_lead_deg=round(up_lead, 4),
            L_wing_lag_deg=round(up_lag, 4),
            L_twist_lead_deg=round(tw_lead, 4),
            L_twist_lag_deg=round(tw_lag, 4),
            L_ant_deg=round(ant_l, 4),
            R_wing_deg=round(up_r, 4),
            R_ant_both_deg=round(ab, 4),
            R_ant_right_deg=round(ar, 4),
            R_loc=[round(v, 6) for v in moth_r["root"].location],
        ))

    # dense keys -> LINEAR, so the evaluated curve is the maths above
    for _act in bpy.data.actions:
        for _fc in _act.fcurves:
            for _kp in _fc.keyframe_points:
                _kp.interpolation = "LINEAR"

    # left shadow: held while perched, then fades and drifts as the moth leaves
    op_l = sh_mat_l.node_tree.nodes["opacity"]
    for f, v in [(0, 0.95), (17, 0.95), (21, 0.62), (26, 0.30),
                 (32, 0.10), (40, 0.03), (47, 0.0)]:
        op_l.outputs[0].default_value = v
        op_l.outputs[0].keyframe_insert("default_value", frame=f)
    base_loc = REST_LOC_L + SH_OFF
    base_sc = SCALE_L
    for f, (dx, dy, sc) in {0: (0.0, 0.0, 1.0), 17: (0.0, 0.0, 1.0),
                            22: (-0.09, 0.05, 1.03), 28: (-0.20, 0.11, 1.08),
                            36: (-0.28, 0.15, 1.12), 47: (-0.34, 0.18, 1.15)}.items():
        sh_l.location = base_loc + Vector((dx, dy, 0.0))
        sh_l.keyframe_insert("location", frame=f)
        sh_l.scale = (base_sc * sc,) * 2 + (1.0,)
        sh_l.keyframe_insert("scale", frame=f)

    with open(os.path.join(HERE, "_dbg", "motion-log.json"), "w", encoding="utf-8") as _fh:
        json.dump(dict(seed=SEED, base_hz=BASE_HZ, base_period=BASE_T,
                       beat_starts=BEAT_T, table=MOTION_TABLE, frames=log_rows),
                  _fh, ensure_ascii=False, indent=1)
    print("motion log ->", os.path.join(HERE, "_dbg", "motion-log.json"))

    # --- PROTO3 STAGE1 debug views (Change 10) ---
    _single = SILHOUETTE or (OBLIQUE_DEG is not None)
    if _single:
        for _o in [moth_l["root"]] + [w for w, _, _ in moth_l["wings"]] + moth_l["body"] + moth_l["legs"]:
            try:
                _o.hide_render = True
                _o.hide_viewport = True
            except Exception:
                pass
        try:
            sh_l.hide_render = True
            sh_l.hide_viewport = True
        except Exception:
            pass
        try:
            R.rotation_euler = (0.0, 0.0, 0.0)
            R.location = REST_LOC_R
            R.scale = (SCALE_R,) * 3
        except Exception:
            pass
        for _w, _k, _s in moth_r["wings"]:
            try:
                _w.animation_data_clear()
            except Exception:
                pass
            _w.rotation_euler = (0.0, 0.0, 0.0)
        try:
            R.animation_data_clear()
        except Exception:
            pass
    if NOSHADOW:
        try:
            sh_l.hide_render = True
        except Exception:
            pass
        try:
            sh_r.hide_render = True
        except Exception:
            pass
    if SILHOUETTE:
        _wm = bpy.data.materials.new("WhiteEmit")
        _wm.use_nodes = True
        _nt = _wm.node_tree
        _nt.nodes.clear()
        _out = _nt.nodes.new("ShaderNodeOutputMaterial")
        _emi = _nt.nodes.new("ShaderNodeEmission")
        _emi.inputs["Color"].default_value = (1.0, 1.0, 1.0, 1.0)
        _emi.inputs["Strength"].default_value = 1.0
        _nt.links.new(_emi.outputs["Emission"], _out.inputs["Surface"])
        for _o in [w for w, _, _ in moth_r["wings"]] + moth_r["body"] + moth_r["legs"]:
            try:
                _o.data.materials.clear()
                _o.data.materials.append(_wm)
            except Exception:
                pass
        try:
            plate.hide_render = True
            plate.hide_viewport = True
        except Exception:
            pass
        try:
            sh_r.hide_render = True
            sh_r.hide_viewport = True
        except Exception:
            pass
        bg.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
        try:
            bg.inputs["Strength"].default_value = 1.0
        except Exception:
            pass
        _cd = bpy.data.cameras.new("TopOrtho")
        _cd.type = "ORTHO"
        _cd.ortho_scale = 3.4
        _co = bpy.data.objects.new("TopOrthoObj", _cd)
        coll.objects.link(_co)
        _co.location = (REST_LOC_R.x, REST_LOC_R.y, 5.0)
        _co.rotation_euler = (0.0, 0.0, 0.0)
        scene.camera = _co
    elif OBLIQUE_DEG is not None:
        try:
            plate.hide_render = True
            plate.hide_viewport = True
        except Exception:
            pass
        try:
            sh_r.hide_render = True
            sh_r.hide_viewport = True
        except Exception:
            pass
        bg.inputs["Color"].default_value = (0.35, 0.35, 0.35, 1.0)
        try:
            bg.inputs["Strength"].default_value = 1.0
        except Exception:
            pass
        _cx, _cy, _cz = REST_LOC_R.x, REST_LOC_R.y, 0.05
        _el = math.radians(float(OBLIQUE_DEG))
        _az = math.radians(-30.0)
        _d = 3.2
        _px = _cx + _d * math.cos(_el) * math.cos(_az)
        _py = _cy + _d * math.cos(_el) * math.sin(_az)
        _pz = _cz + _d * math.sin(_el)
        _cd2 = bpy.data.cameras.new("ObliqueCam")
        _cd2.type = "PERSP"
        _cd2.lens = 50.0
        _co2 = bpy.data.objects.new("ObliqueCamObj", _cd2)
        coll.objects.link(_co2)
        _co2.location = (_px, _py, _pz)
        _dir = Vector((_cx - _px, _cy - _py, _cz - _pz))
        _co2.rotation_euler = _dir.to_track_quat("-Z", "Y").to_euler()
        scene.camera = _co2

    scene.frame_set(0)
    img = bpy.data.images.get("background-bark-832x480.png")
    if img:
        try:
            img.pack()
        except Exception as exc:
            print("pack failed:", exc)
    bpy.ops.wm.save_as_mainfile(filepath=BLEND_OUT)
    print("saved", BLEND_OUT)

    params = dict(blender=bpy.app.version_string, sun=SUN_STRENGTH, ambient=AMBIENT,
                  view_w=VIEW_W, frames=[F_START, F_END], fps=FPS,
                  objects=len(bpy.data.objects), meshes=len(bpy.data.meshes),
                  wing_objects=len(moth_l["wings"]) + len(moth_r["wings"]),
                  flap_hz=BASE_HZ, stroke_top_deg=STROKE_TOP, proto=4,
                  choreography="proto4-C v2 (measured raise -> hold -> liftoff)",
                  seed=SEED,
                  beat_period_factors=MOTION_TABLE["left"]["beat_period_factors"],
                  beat_starts=[round(t, 4) for t in BEAT_T[:5]],
                  motion_log=os.path.join(HERE, "_dbg", "motion-log.json"))
    with open(os.path.join(HERE, "_dbg", "build_params.json"), "w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2)
    print("scene:", json.dumps(params))

    if STILL is not None:
        os.makedirs(os.path.dirname(STILL_OUT), exist_ok=True)
        scene.frame_set(STILL)
        scene.render.filepath = STILL_OUT
        bpy.ops.render.render(write_still=True)
        print("still ->", STILL_OUT)

    if DO_RENDER:
        out_dir = os.path.join(HERE, "renders")
        os.makedirs(out_dir, exist_ok=True)
        scene.render.filepath = os.path.join(out_dir, "frame_")
        scene.render.use_file_extension = True
        bpy.ops.render.render(animation=True)
        print("rendered into", out_dir)


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "_dbg"), exist_ok=True)
    main()
