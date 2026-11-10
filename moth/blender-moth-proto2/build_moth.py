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
    ra.elements[0].position = 0.18; ra.elements[0].color = (*srgb("#55755f"), 1)
    ra.elements[1].position = 0.52; ra.elements[1].color = (*srgb("#6e9c90"), 1)
    e = ra.elements.new(0.86); e.color = (*srgb("#83afa4"), 1)
    L(nzA.outputs["Fac"], rampA.inputs["Fac"])
    col = rampA.outputs["Color"]

    # ---------------- darker outer margin
    margF = map_range((-1080, 420), t_, 0.84, 1.0, 0.0, 0.34)
    col = mix_rgba((-860, 300), margF, col, const_rgb((-860, 460), "#5f826c"))

    # ---------------- fore/hind crease + hindwing tint
    # boundary curve in UV: u_line(t) = 0.42 + 0.58 t
    line = math_node("MULTIPLY_ADD", (-1080, 640))
    L(t_, line.inputs[0]); line.inputs[1].default_value = 0.58; line.inputs[2].default_value = 0.42
    dsub = math_node("SUBTRACT", (-900, 640))
    L(u_, dsub.inputs[0]); L(line.outputs[0], dsub.inputs[1])
    hindF0 = map_range((-720, 700), dsub.outputs[0], 0.0, 0.07, 0.0, 1.0)
    hindF = math_node("MULTIPLY", (-560, 700))
    L(hindF0, hindF.inputs[0]); hindF.inputs[1].default_value = 0.38
    col = mix_rgba((-560, 420), hindF.outputs[0], col, const_rgb((-560, 560), "#7cab9f"))
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
    vfade = map_range((-200, -480), t_, 0.55, 1.0, 1.0, 0.3, smooth=False)
    veinAmt = math_node("MULTIPLY", (-40, -400))
    L(vein.outputs["Color"], veinAmt.inputs[0]); L(vfade, veinAmt.inputs[1])
    veinAmt2 = math_node("MULTIPLY", (60, -400))
    L(veinAmt.outputs[0], veinAmt2.inputs[0]); veinAmt2.inputs[1].default_value = 0.30
    col = mix_rgba((60, 300), veinAmt2.outputs[0], col, const_rgb((60, 460), "#59755f"))

    # ---------------- transverse bands (bowed, wobbling, broken)
    # shared end fade along the chord
    ef1 = map_range((-1300, -700), u_, 0.0, 0.10, 0.0, 1.0)
    ef2 = map_range((-1300, -820), u_, 1.0, 0.90, 0.0, 1.0)
    efm = math_node("MULTIPLY", (-1120, -760)); L(ef1, efm.inputs[0]); L(ef2, efm.inputs[1])
    endfade = math_node("MULTIPLY_ADD", (-960, -760))
    L(efm.outputs[0], endfade.inputs[0]); endfade.inputs[1].default_value = 0.70; endfade.inputs[2].default_value = 0.30

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
        prof = map_range((0, y0), dist.outputs[0], 0.0, width, 1.0, 0.0)
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
        return amt3.outputs[0]

    b1 = build_band(0.50, 0.060, 0.026, 0.38, 0.0, 0.38, 0.62, 0.15, 120)
    b2 = build_band(0.72, 0.075, 0.030, 0.62, 3.7, 0.38, 0.62, 0.15, -60)
    b3 = build_band(0.875, 0.060, 0.020, 0.34, 7.3, 0.45, 0.60, 0.05, -240)
    bsum = math_node("ADD", (560, 60)); L(b1, bsum.inputs[0]); L(b2, bsum.inputs[1])
    bsum2 = math_node("ADD", (680, 60)); L(bsum.outputs[0], bsum2.inputs[0]); L(b3, bsum2.inputs[1])
    bclamp = math_node("MINIMUM", (800, 60)); L(bsum2.outputs[0], bclamp.inputs[0]); bclamp.inputs[1].default_value = 1.0
    bamt = math_node("MULTIPLY", (920, 60)); L(bclamp.outputs[0], bamt.inputs[0]); bamt.inputs[1].default_value = 0.75
    col = mix_rgba((920, 300), bamt.outputs[0], col, const_rgb((920, 460), "#edf2e9"))

    # ---------------- fine grain + speckle (kept from proto1, subtler)
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

    # ---------------- soft fresnel rim
    lw = N.new("ShaderNodeLayerWeight"); lw.location = (920, -560)
    lw.inputs["Blend"].default_value = 0.30
    rim = N.new("ShaderNodeValToRGB"); rim.location = (1080, -560)
    rr = rim.color_ramp
    rr.elements[0].position = 0.62; rr.elements[0].color = (0, 0, 0, 1)
    rr.elements[1].position = 1.0; rr.elements[1].color = (1, 1, 1, 1)
    L(lw.outputs["Facing"], rim.inputs["Fac"])
    rimw = math_node("MULTIPLY", (1240, -560)); L(rim.outputs["Color"], rimw.inputs[0]); rimw.inputs[1].default_value = 0.08
    col = mix_rgba((1400, 300), rimw.outputs[0], col, const_rgb((1240, -700), "#f2f8f0"))

    L(col, bsdf.inputs["Base Color"])

    rr2 = N.new("ShaderNodeMapRange"); rr2.location = (1080, -340)
    rr2.inputs["To Min"].default_value = 0.66
    rr2.inputs["To Max"].default_value = 0.85
    L(n2.outputs["Fac"], rr2.inputs["Value"])
    L(rr2.outputs["Result"], bsdf.inputs["Roughness"])
    set_in(bsdf, "Specular IOR Level", 0.12)
    set_in(bsdf, "Sheen Weight", 0.15)
    set_in(bsdf, "Sheen Roughness", 0.55)
    if "Sheen Tint" in bsdf.inputs:
        bsdf.inputs["Sheen Tint"].default_value = (*srgb("#dff0da"), 1.0)

    bump = N.new("ShaderNodeBump"); bump.location = (1180, -760)
    bump.inputs["Strength"].default_value = 0.15
    bump.inputs["Distance"].default_value = 0.002
    L(n2.outputs["Fac"], bump.inputs["Height"])
    L(bump.outputs["Normal"], bsdf.inputs["Normal"])

    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


def mat_body(name, c_hex, rough=0.85, fur=0.5, scale=60.0):
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
    set_in(bsdf, "Sheen Weight", 0.4)
    bump = N.new("ShaderNodeBump"); bump.location = (180, -320)
    bump.inputs["Strength"].default_value = fur
    bump.inputs["Distance"].default_value = 0.006
    L(n2.outputs["Fac"], bump.inputs["Height"])
    L(bump.outputs["Normal"], bsdf.inputs["Normal"])
    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


# ---------------------------------------------------------------- contact shadow
# Footprint rect in moth-local units (x span, y along body).
SH_X0, SH_X1 = -0.54, 0.54
SH_Y0, SH_Y1 = -0.36, 0.24


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


def make_shadow_image():
    """Rasterise both fan silhouettes + body ellipse into a soft alpha image."""
    W, H = 416, 336
    lex, ley = sample_curve(FAN_LE, 60)
    tex, tey = sample_curve(FAN_TE, 60)
    polys = []
    for s in (1.0, -1.0):
        poly = [(s * lex[i], ley[i]) for i in range(60)]
        poly += [(s * tex[i], tey[i]) for i in range(59, -1, -1)]
        polys.append(poly)
    # body ellipse (thorax + abdomen + head)
    ep = []
    for i in range(28):
        a = 2 * math.pi * i / 28.0
        ep.append((0.070 * math.cos(a), -0.03 + 0.27 * math.sin(a)))
    polys.append(ep)
    mask = _scanline_fill(polys, W, H)
    sharp = _box_blur(mask, W, H, 2)
    wide = _box_blur(_box_blur(mask, W, H, 14), W, H, 14)
    px = []
    for r in range(H):
        for c in range(W):
            # core is a strong darkening (~0.72), not a black cutout; wide soft halo
            a = 0.68 * sharp[r][c] + 0.50 * wide[r][c]
            a = min(0.72, a)
            px.extend((0.0, 0.0, 0.0, a))
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
FAN_LE = [(0.00, 0.006, 0.095), (0.15, 0.080, 0.130), (0.35, 0.180, 0.148),
          (0.55, 0.280, 0.140), (0.75, 0.375, 0.105), (0.90, 0.445, 0.055),
          (0.97, 0.487, 0.016), (1.00, 0.500, 0.000)]
FAN_TE = [(0.00, 0.000, -0.085), (0.10, 0.032, -0.150), (0.22, 0.080, -0.205),
          (0.35, 0.132, -0.245), (0.50, 0.196, -0.282), (0.65, 0.262, -0.262),
          (0.78, 0.322, -0.222), (0.90, 0.390, -0.158), (0.97, 0.450, -0.078),
          (0.99, 0.493, -0.028), (1.00, 0.496, -0.012)]


def _sstep(a, b, x):
    f = min(1.0, max(0.0, (x - a) / (b - a)))
    return f * f * (3.0 - 2.0 * f)


def make_fan_wing(side, mir, name, coll, mat):
    sign = mir * side
    Nt, Nu = 44, 18
    lex, ley = sample_curve(FAN_LE, Nt + 1)
    tex, tey = sample_curve(FAN_TE, Nt + 1)
    verts, uvs, faces = [], [], []
    for i in range(Nt + 1):
        t = i / float(Nt)
        # scalloped outer margin (only where TE is the outer margin)
        sc = _sstep(0.42, 0.80, t) * _sstep(1.0, 0.965, t)
        ty = tey[i] + 0.0075 * math.sin(2 * math.pi * 6.0 * t + 0.7) * sc
        tx = tex[i]
        for j in range(Nu + 1):
            u = j / float(Nu)
            x = lex[i] + (tx - lex[i]) * u
            y = ley[i] + (ty - ley[i]) * u
            z = 0.010 * math.sin(math.pi * u) * (0.35 + 0.65 * t) \
                + 0.0022 * math.sin(2 * math.pi * 2.5 * t + 3.0 * u) * t \
                - 0.012 * t * t
            verts.append(Vector((sign * x, y, z)))
            uvs.append((t, u))
    for i in range(Nt):
        for j in range(Nu):
            a = i * (Nu + 1) + j
            b = a + 1
            c = a + (Nu + 1)
            d = c + 1
            if side > 0:
                faces.append([a, b, d, c])
            else:
                faces.append([a, c, d, b])
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
        w.location = Vector((mirror * side * 0.030, 0.012, 0.026))
        w.parent = root
        wing_objs.append((w, "fan", side))

    body_parts = []
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=24, v_segments=16, radius=1.0)
    for v in bm.verts:
        v.co.x *= 0.058
        v.co.y *= 0.082
        v.co.z *= 0.055
        v.co.z += 0.040
        v.co.y += 0.028
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
        r = 0.036 * (1.0 - t) ** 0.60 + 0.005
        r *= 1.0 + 0.05 * math.sin(2 * math.pi * 4.0 * t)
        prof.append((y, r))
    ring = 18
    verts, faces = [], []
    for i, (y, r) in enumerate(prof):
        t = i / float(nseg)
        z = 0.042 - 0.024 * t * t
        for k in range(ring):
            a = 2 * math.pi * k / ring
            verts.append((mirror * r * math.cos(a) * 0.92, y, z + r * math.sin(a)))
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
        v.co.z += 0.046
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
            p0 = Vector((mirror * side * 0.026, y0, 0.038))
            p1 = Vector((mirror * side * (0.026 + dx * 0.55), y0 + dy * 0.45, 0.018))
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
                0.124 + 0.150 * t - 0.030 * t * t,
                0.058 + 0.055 * t - 0.010 * t * t)))
        radii = [0.0036 * (1.0 - 0.55 * (i / 8.0)) for i in range(9)]
        an = tube("%s_antenna_%s" % (name, "R" if side > 0 else "L"), pts, radii, coll, segs=6)
        an.data.materials.append(mats["leg"])
        an.parent = root
        leg_objs.append(an)
        for i in range(1, 7):
            t = i / 7.0
            p = Vector((
                mirror * side * (0.014 + 0.055 * t),
                0.124 + 0.150 * t - 0.030 * t * t,
                0.058 + 0.055 * t - 0.010 * t * t))
            q = p + Vector((mirror * side * 0.011, -0.019, 0.004))
            bar = tube("%s_pecten_%s%d" % (name, "R" if side > 0 else "L", i),
                       [p, q], [0.0020, 0.0011], coll, segs=5)
            bar.data.materials.append(mats["leg"])
            bar.parent = root
            leg_objs.append(bar)

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
        ee.gtao_factor = 0.6
    ee.use_shadows = True
    scene.render.use_motion_blur = True
    scene.render.motion_blur_shutter = 0.5

    mats = dict(
        wing=mat_wing(),
        body=mat_body("BodyMat", "#6f8f66", rough=0.88, fur=0.55, scale=70.0),
        abdomen=mat_body("AbdomenMat", "#6d8a64", rough=0.90, fur=0.55, scale=85.0),
        leg=mat_body("LegMat", "#4a4030", rough=0.66, fur=0.20, scale=120.0),
        eye=mat_body("EyeMat", "#191a12", rough=0.35, fur=0.0, scale=10.0),
        bark=mat_background(),
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
    # offset down-right in frame space (sun comes from the upper left)
    sh_img = make_shadow_image()
    try:
        sh_img.filepath_raw = os.path.join(HERE, "_dbg", "shadow_tex.png")
        sh_img.file_format = "PNG"
        sh_img.save()
        print("shadow tex ->", sh_img.filepath_raw)
    except Exception as exc:
        print("shadow tex save failed:", exc)
    sh_mat_l = mat_shadow_silhouette(sh_img)
    sh_mat_r = mat_shadow_silhouette(sh_img)
    SH_OFF = Vector((0.075, -0.075, 0.0))
    sh_l = build_shadow_decal("Shadow_L", sh_img, sh_mat_l)
    sh_l.location = REST_LOC_L + SH_OFF
    sh_l.rotation_euler = (0.0, 0.0, math.radians(50.0))
    sh_l.scale = (SCALE_L,) * 3
    sh_r = build_shadow_decal("Shadow_R", sh_img, sh_mat_r)
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
    fps = FPS

    def flap_sign(mir, side):
        return -(mir * side)

    def set_wing(w, kind, side, mir, up_deg):
        s = flap_sign(mir, side)
        w.rotation_euler[1] = math.radians(s * up_deg)

    # --- LEFT moth: rests FLAT (3 deg), then flaps away
    L = moth_l["root"]
    REST_L = 3.0
    for w, kind, side in moth_l["wings"]:
        set_wing(w, kind, side, -1, REST_L)
    for f in (0, 9):
        for w, kind, side in moth_l["wings"]:
            w.keyframe_insert("rotation_euler", index=1, frame=f)

    FLAP_HZ = 9.5
    AMP = 40.0
    for f in range(10, F_END + 1):
        ramp = min(1.0, max(0.0, (f - 10) / 3.0))
        amp = AMP * (0.25 + 0.75 * ramp)
        ph = 2 * math.pi * FLAP_HZ * (f - 10) / fps
        for w, kind, side in moth_l["wings"]:
            up = 6.0 * ramp + amp * math.sin(ph)
            set_wing(w, kind, side, -1, up)
            w.keyframe_insert("rotation_euler", index=1, frame=f)
            w.rotation_euler[0] = math.radians(7.0 * math.cos(ph) * ramp)
            w.keyframe_insert("rotation_euler", index=0, frame=f)

    LP = {
        0: (fx(0.300), fy(0.520), 0.000),
        9: (fx(0.300), fy(0.520), 0.010),
        14: (fx(0.285), fy(0.505), 0.130),
        20: (fx(0.205), fy(0.415), 0.430),
        28: (fx(0.045), fy(0.250), 0.760),
        36: (fx(-0.240), fy(0.010), 1.020),
        47: (fx(-0.560), fy(-0.260), 1.220),
    }
    for f, p in LP.items():
        L.location = Vector(p)
        L.keyframe_insert("location", frame=f)

    ROLL = {0: 0.0, 9: 0.0, 14: -6.0, 20: -20.0, 28: -30.0, 36: -34.0, 47: -34.0}
    PITCH = {0: 0.0, 9: 0.0, 14: 12.0, 20: 18.0, 28: 22.0, 36: 20.0, 47: 18.0}
    YAW = {0: 0.0, 9: 0.0, 14: 4.0, 20: 10.0, 28: 16.0, 36: 20.0, 47: 22.0}
    for f in LP:
        L.rotation_euler = (math.radians(PITCH[f]), math.radians(ROLL[f]), math.radians(50.0 + YAW[f]))
        L.keyframe_insert("rotation_euler", frame=f)

    # left shadow: fades and drifts as the moth leaves
    op_l = sh_mat_l.node_tree.nodes["opacity"]
    for f, v in [(0, 0.95), (9, 0.95), (14, 0.62), (22, 0.25), (30, 0.07), (47, 0.0)]:
        op_l.outputs[0].default_value = v
        op_l.outputs[0].keyframe_insert("default_value", frame=f)
    base_loc = REST_LOC_L + SH_OFF
    base_sc = SCALE_L
    for f, (dx, dy, sc) in {0: (0.0, 0.0, 1.0), 9: (0.0, 0.0, 1.0),
                            16: (-0.07, 0.04, 1.03), 24: (-0.18, 0.10, 1.08),
                            32: (-0.28, 0.15, 1.12), 47: (-0.34, 0.18, 1.15)}.items():
        sh_l.location = base_loc + Vector((dx, dy, 0.0))
        sh_l.keyframe_insert("location", frame=f)
        sh_l.scale = (base_sc * sc,) * 2 + (1.0,)
        sh_l.keyframe_insert("scale", frame=f)

    # --- RIGHT moth: stays put, slow shallow open/close around flat
    R = moth_r["root"]
    for w, kind, side in moth_r["wings"]:
        set_wing(w, kind, side, 1, 3.0)
    for f, up in [(0, 3.0), (6, 6.5), (14, 3.0), (23, 3.0), (28, 6.5), (37, 3.0), (47, 3.0)]:
        for w, kind, side in moth_r["wings"]:
            set_wing(w, kind, side, 1, up)
            w.keyframe_insert("rotation_euler", index=1, frame=f)

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
                  flap_hz=FLAP_HZ, flap_amp_deg=AMP, proto=2)
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
