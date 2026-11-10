"""build_moth.py - procedural emerald-geometer moth pair + bark scene (Blender 4.2 / Eevee).

Everything is generated here, no manual modelling:
  * wings: hand-authored outlines (span/chord profiles), cambered, hinged at the
    thorax as separate objects, UV'd so materials can run spanwise bands and
    chordwise veins. Forewing = triangular, hindwing = rounder (Geometridae).
  * body: thorax / abdomen (segmented lathe) / head / eyes / 6 legs / 2 antennae
    with a hint of pectination.
  * wing material: procedural only. Pale mint base (#a8d8b0), two wavy white
    bands, veins, a slightly darker margin, scale speckle, soft Fresnel rim.
  * shadows: black radial-gradient decals (procedural alpha) under each moth.
  * bark plate: the inpainted photo, unlit so it reproduces exactly.
  * fixed camera, one soft sun + weak ambient, 48 frames @ 24 fps, motion blur.

The photo is used ONLY as the flat bark plate behind the moths (and, before
this, as an outline reference). No photo pixel is mapped onto a moth.

Run examples
  python run_blender.py build_moth.py
  python run_blender.py build_moth.py --sun 1.3 --ambient 1.0 --still 0 --out _dbg/t.png
  python run_blender.py build_moth.py --render

Frame math: 24 fps, frames 0..47 (the brief's f0..f47).
World scale: frame width = 3.0 units -> 1 unit = 277.3 px at 832 px wide.
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
BLEND_OUT = os.path.join(HERE, "proto.blend")

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
    """UV.x = span (0 base -> 1 tip), UV.y = chord (0 costa -> 1 dorsum)."""
    m = bpy.data.materials.new("WingProcedural")
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    N, L = nt.nodes, nt.links.new

    out = N.new("ShaderNodeOutputMaterial"); out.location = (1500, 0)
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (1200, 0)

    tc = N.new("ShaderNodeTexCoord"); tc.location = (-1400, 0)
    sep = N.new("ShaderNodeSeparateXYZ"); sep.location = (-1200, 220)
    L(tc.outputs["UV"], sep.inputs["Vector"])
    t_ = sep.outputs["X"]   # span
    u_ = sep.outputs["Y"]   # chord

    def math_node(op, loc, a=None, b=None):
        n = N.new("ShaderNodeMath"); n.operation = op; n.location = loc
        return n

    # ---------------- two wavy white bands across the wing (vary along span)
    nz = N.new("ShaderNodeTexNoise"); nz.location = (-1200, -60)
    nz.inputs["Scale"].default_value = 5.0
    nz.inputs["Detail"].default_value = 3.0
    nz.inputs["Roughness"].default_value = 0.55
    L(tc.outputs["Object"], nz.inputs["Vector"])

    phase = math_node("MULTIPLY_ADD", (-1000, 60))
    L(nz.outputs["Fac"], phase.inputs[0])
    phase.inputs[1].default_value = 1.9
    phase.inputs[2].default_value = -0.95

    spanf = math_node("MULTIPLY_ADD", (-820, 220))
    L(t_, spanf.inputs[0])
    spanf.inputs[1].default_value = 2.0 * math.pi * 2.0
    spanf.inputs[2].default_value = 0.35

    ssum = math_node("ADD", (-640, 160))
    L(spanf.outputs[0], ssum.inputs[0])
    L(phase.outputs[0], ssum.inputs[1])
    ssin = math_node("SINE", (-460, 160)); L(ssum.outputs[0], ssin.inputs[0])
    sabs = math_node("ABSOLUTE", (-300, 160)); L(ssin.outputs[0], sabs.inputs[0])

    band = N.new("ShaderNodeValToRGB"); band.location = (-120, 160)
    cr = band.color_ramp
    cr.interpolation = "EASE"
    cr.elements[0].position = 0.0; cr.elements[0].color = (0, 0, 0, 1)
    cr.elements[1].position = 0.86; cr.elements[1].color = (0, 0, 0, 1)
    e = cr.elements.new(0.95); e.color = (0.75, 0.75, 0.75, 1)
    e = cr.elements.new(1.0); e.color = (1, 1, 1, 1)
    L(sabs.outputs[0], band.inputs["Fac"])

    # break the band up with fine noise so it is not a clean stripe
    bmod = math_node("MULTIPLY", (60, 160))
    L(band.outputs["Color"], bmod.inputs[0])
    bmod.inputs[1].default_value = 1.0
    bandmix = math_node("MULTIPLY_ADD", (240, 160))
    L(nz.outputs["Fac"], bandmix.inputs[0])
    bandmix.inputs[1].default_value = 0.5
    bandmix.inputs[2].default_value = 0.35
    bmul = math_node("MULTIPLY", (420, 160))
    L(bmod.outputs[0], bmul.inputs[0])
    L(bandmix.outputs[0], bmul.inputs[1])
    bcap = math_node("MULTIPLY", (500, 60))
    L(bmul.outputs[0], bcap.inputs[0])
    bcap.inputs[1].default_value = 0.62

    # ---------------- veins: fine darker lines running base -> margin
    vf = math_node("MULTIPLY", (-820, -240)); L(u_, vf.inputs[0])
    vf.inputs[1].default_value = 2.0 * math.pi * 7.0
    vs = math_node("SINE", (-640, -240)); L(vf.outputs[0], vs.inputs[0])
    va = math_node("ABSOLUTE", (-460, -240)); L(vs.outputs[0], va.inputs[0])
    vein = N.new("ShaderNodeValToRGB"); vein.location = (-280, -240)
    vr = vein.color_ramp
    vr.interpolation = "EASE"
    vr.elements[0].position = 0.0; vr.elements[0].color = (1, 1, 1, 1)
    vr.elements[1].position = 0.30; vr.elements[1].color = (0, 0, 0, 1)
    L(va.outputs[0], vein.inputs["Fac"])
    # veins fade out near the margin
    vfade = N.new("ShaderNodeMapRange"); vfade.location = (-280, -420)
    vfade.inputs["From Min"].default_value = 0.55
    vfade.inputs["From Max"].default_value = 1.0
    vfade.inputs["To Min"].default_value = 1.0
    vfade.inputs["To Max"].default_value = 0.25
    L(t_, vfade.inputs["Value"])

    # ---------------- darker outer margin
    marg = N.new("ShaderNodeMapRange"); marg.location = (-460, 400)
    marg.inputs["From Min"].default_value = 0.83
    marg.inputs["From Max"].default_value = 1.0
    marg.inputs["To Min"].default_value = 0.0
    marg.inputs["To Max"].default_value = 0.55
    L(t_, marg.inputs["Value"])

    # ---------------- scale speckle + grain
    vor = N.new("ShaderNodeTexVoronoi"); vor.location = (-1200, -420)
    vor.inputs["Scale"].default_value = 120.0
    L(tc.outputs["Object"], vor.inputs["Vector"])
    speck = N.new("ShaderNodeValToRGB"); speck.location = (-1000, -420)
    speck.color_ramp.elements[0].position = 0.02
    speck.color_ramp.elements[0].color = (0, 0, 0, 1)
    speck.color_ramp.elements[1].position = 0.22
    speck.color_ramp.elements[1].color = (0.5, 0.5, 0.5, 1)
    L(vor.outputs["Distance"], speck.inputs["Fac"])

    n2 = N.new("ShaderNodeTexNoise"); n2.location = (-1200, -640)
    n2.inputs["Scale"].default_value = 55.0
    n2.inputs["Detail"].default_value = 6.0
    n2.inputs["Roughness"].default_value = 0.6
    L(tc.outputs["Object"], n2.inputs["Vector"])
    grain = N.new("ShaderNodeValToRGB"); grain.location = (-1000, -640)
    gr = grain.color_ramp
    gr.elements[0].position = 0.25; gr.elements[0].color = (0.90, 0.92, 0.90, 1)
    gr.elements[1].position = 0.75; gr.elements[1].color = (1.04, 1.03, 1.02, 1)
    L(n2.outputs["Fac"], grain.inputs["Fac"])

    # ---------------- colour assembly
    base = N.new("ShaderNodeRGB"); base.location = (-120, 620)
    base.outputs[0].default_value = (*srgb("#8ccb9a"), 1.0)
    pale = N.new("ShaderNodeRGB"); pale.location = (-120, 500)
    pale.outputs[0].default_value = (*srgb("#dbe8d6"), 1.0)
    mix1 = N.new("ShaderNodeMix"); mix1.data_type = "RGBA"; mix1.location = (620, 480)
    L(bcap.outputs[0], mix1.inputs["Factor"])
    L(base.outputs[0], mix1.inputs[6])
    L(pale.outputs[0], mix1.inputs[7])

    # veins darken (green-brown, not black)
    veincol = N.new("ShaderNodeRGB"); veincol.location = (420, -180)
    veincol.outputs[0].default_value = (*srgb("#7b8f6d"), 1.0)
    vmix = N.new("ShaderNodeMix"); vmix.data_type = "RGBA"; vmix.location = (800, 380)
    L(vein.outputs["Color"], vmix.inputs["Factor"])
    L(mix1.outputs[2], vmix.inputs[6])
    L(veincol.outputs[0], vmix.inputs[7])
    vmul = math_node("MULTIPLY", (620, -300))
    L(vein.outputs["Color"], vmul.inputs[0])
    L(vfade.outputs["Result"], vmul.inputs[1])
    vsc = math_node("MULTIPLY", (760, -300))
    L(vmul.outputs[0], vsc.inputs[0])
    vsc.inputs[1].default_value = 0.55
    L(vsc.outputs[0], vmix.inputs["Factor"])

    # margin darkening
    mcol = N.new("ShaderNodeRGB"); mcol.location = (620, 260)
    mcol.outputs[0].default_value = (*srgb("#87a98c"), 1.0)
    mmix = N.new("ShaderNodeMix"); mmix.data_type = "RGBA"; mmix.location = (980, 320)
    L(marg.outputs["Result"], mmix.inputs["Factor"])
    L(vmix.outputs[2], mmix.inputs[6])
    L(mcol.outputs[0], mmix.inputs[7])

    # grain + speckle multiply
    gmix = N.new("ShaderNodeMix"); gmix.data_type = "RGBA"; gmix.blend_type = "MULTIPLY"
    gmix.location = (1160, 320)
    gmix.inputs["Factor"].default_value = 0.8
    L(mmix.outputs[2], gmix.inputs[6])
    L(grain.outputs["Color"], gmix.inputs[7])
    smix = N.new("ShaderNodeMix"); smix.data_type = "RGBA"; smix.blend_type = "SCREEN"
    smix.location = (1300, 320)
    smix.inputs["Factor"].default_value = 0.10
    L(gmix.outputs[2], smix.inputs[6])
    L(speck.outputs["Color"], smix.inputs[7])

    # fresnel rim (soft sheen suggesting a translucent edge)
    lw = N.new("ShaderNodeLayerWeight"); lw.location = (620, -560)
    lw.inputs["Blend"].default_value = 0.30
    rim = N.new("ShaderNodeValToRGB"); rim.location = (800, -560)
    rr = rim.color_ramp
    rr.elements[0].position = 0.62; rr.elements[0].color = (0, 0, 0, 1)
    rr.elements[1].position = 1.0; rr.elements[1].color = (1, 1, 1, 1)
    L(lw.outputs["Facing"], rim.inputs["Fac"])
    rimcol = N.new("ShaderNodeRGB"); rimcol.location = (800, -700)
    rimcol.outputs[0].default_value = (*srgb("#f4fbf2"), 1.0)
    rmix = N.new("ShaderNodeMix"); rmix.data_type = "RGBA"; rmix.location = (1400, 300)
    L(smix.outputs[2], rmix.inputs[6])
    L(rimcol.outputs[0], rmix.inputs[7])
    rimw = math_node("MULTIPLY", (1000, -560)); L(rim.outputs["Color"], rimw.inputs[0])
    rimw.inputs[1].default_value = 0.10
    L(rimw.outputs[0], rmix.inputs["Factor"])

    L(rmix.outputs[2], bsdf.inputs["Base Color"])

    rr2 = N.new("ShaderNodeMapRange"); rr2.location = (980, -340)
    rr2.inputs["To Min"].default_value = 0.68
    rr2.inputs["To Max"].default_value = 0.86
    L(n2.outputs["Fac"], rr2.inputs["Value"])
    L(rr2.outputs["Result"], bsdf.inputs["Roughness"])
    set_in(bsdf, "Specular IOR Level", 0.16)
    set_in(bsdf, "Sheen Weight", 0.22)
    set_in(bsdf, "Sheen Roughness", 0.55)
    if "Sheen Tint" in bsdf.inputs:
        bsdf.inputs["Sheen Tint"].default_value = (*srgb("#e6f4e2"), 1.0)

    bump = N.new("ShaderNodeBump"); bump.location = (1000, -760)
    bump.inputs["Strength"].default_value = 0.22
    bump.inputs["Distance"].default_value = 0.003
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


def mat_shadow():
    m = bpy.data.materials.new("ShadowDecal")
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
    bsdf = N.new("ShaderNodeBsdfPrincipled"); bsdf.location = (400, 0)
    bsdf.inputs["Base Color"].default_value = (0.003, 0.004, 0.003, 1)
    bsdf.inputs["Roughness"].default_value = 1.0
    set_in(bsdf, "Specular IOR Level", 0.0)
    tc = N.new("ShaderNodeTexCoord"); tc.location = (-760, 0)
    sub = N.new("ShaderNodeVectorMath"); sub.operation = "SUBTRACT"; sub.location = (-560, 0)
    L(tc.outputs["Generated"], sub.inputs[0])
    sub.inputs[1].default_value = (0.5, 0.5, 0.5)
    sc = N.new("ShaderNodeVectorMath"); sc.operation = "MULTIPLY"; sc.location = (-380, 0)
    L(sub.outputs[0], sc.inputs[0])
    sc.inputs[1].default_value = (2.0, 2.0, 0.0)
    ln = N.new("ShaderNodeVectorMath"); ln.operation = "LENGTH"; ln.location = (-200, 0)
    L(sc.outputs[0], ln.inputs[0])
    ramp = N.new("ShaderNodeValToRGB"); ramp.location = (0, 0)
    cr = ramp.color_ramp
    cr.interpolation = "EASE"
    cr.elements[0].position = 0.0; cr.elements[0].color = (1, 1, 1, 1)
    cr.elements[1].position = 0.38; cr.elements[1].color = (0.88, 0.88, 0.88, 1)
    e = cr.elements.new(0.64); e.color = (0.50, 0.50, 0.50, 1)
    e = cr.elements.new(0.86); e.color = (0.14, 0.14, 0.14, 1)
    e = cr.elements.new(1.0); e.color = (0.0, 0.0, 0.0, 1)
    L(ln.outputs["Value"], ramp.inputs["Fac"])
    op = N.new("ShaderNodeValue"); op.location = (0, -260); op.label = "opacity"; op.name = "opacity"
    op.outputs[0].default_value = 1.0
    mul = N.new("ShaderNodeMath"); mul.operation = "MULTIPLY"; mul.location = (320, -160)
    L(ramp.outputs["Color"], mul.inputs[0])
    L(op.outputs[0], mul.inputs[1])
    L(mul.outputs[0], bsdf.inputs["Alpha"])
    L(bsdf.outputs["BSDF"], out.inputs["Surface"])
    return m


# ---------------------------------------------------------------- wing geometry
WING_SPECS = {
    "fore": dict(L=0.440, root=0.012, camber=0.024,
                 ycen=[(0.0, 0.080), (0.25, 0.072), (0.55, 0.036), (0.80, -0.002), (1.0, -0.028)],
                 hw=[(0.0, 0.048), (0.15, 0.090), (0.35, 0.106), (0.55, 0.101),
                     (0.72, 0.077), (0.86, 0.045), (0.95, 0.018), (1.0, 0.004)]),
    "hind": dict(L=0.320, root=0.010, camber=0.020, sweep=-0.34,
                 ycen=[(0.0, 0.014), (0.3, -0.014), (0.6, -0.042), (1.0, -0.068)],
                 hw=[(0.0, 0.046), (0.25, 0.090), (0.45, 0.100), (0.65, 0.088),
                     (0.82, 0.058), (0.93, 0.026), (1.0, 0.004)]),
}


def make_wing(kind, side, mir, name, coll, mat):
    spec = WING_SPECS[kind]
    Nt, Nu = 24, 14
    yc = smooth_profile(spec["ycen"], Nt + 1)
    hw = smooth_profile(spec["hw"], Nt + 1)
    L_, camber = spec["L"], spec["camber"]
    sweep = spec.get("sweep", 0.0)
    verts, uvs, faces = [], [], []
    for i in range(Nt + 1):
        t = i / float(Nt)
        ang = sweep * t
        ca, sa = math.cos(ang), math.sin(ang)
        out = spec["root"] + L_ * t
        for j in range(Nu + 1):
            u = j / float(Nu)
            x = out * ca
            y = yc[i] + (2.0 * u - 1.0) * hw[i]
            if j == Nu:                       # scalloped trailing edge
                y += 0.0042 * math.sin(2 * math.pi * 4.5 * t) * t
            if j == 0:                        # very slight lead-edge wave
                y += 0.0016 * math.sin(2 * math.pi * 2.2 * t + 0.6) * t
            xr = x * ca - (y - yc[i]) * sa * 0.35
            yr = y + out * sa
            z = camber * math.sin(math.pi * u) * (0.22 + 0.78 * t) \
                + 0.010 * t * t \
                - 0.008 * (u - 0.5) * t
            verts.append(Vector((mir * side * xr, yr, z)))
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
    solid.thickness = 0.0018
    solid.offset = 0.0
    solid.use_even_offset = True
    return ob


def hinge_point(kind, mirror, side):
    if kind == "fore":
        return Vector((mirror * side * 0.042, 0.040, 0.028))
    return Vector((mirror * side * 0.038, -0.020, 0.020))


# ---------------------------------------------------------------- moth
def build_moth(name, mirror, coll, mats):
    root = bpy.data.objects.new(name, None)
    root.empty_display_type = "PLAIN_AXES"
    root.empty_display_size = 0.08
    coll.objects.link(root)

    wing_objs = []
    for side in (1, -1):
        for kind in ("fore", "hind"):
            wname = "%s_wing_%s_%s" % (name, kind, "R" if side > 0 else "L")
            w = make_wing(kind, side, mirror, wname, coll, mats["wing"])
            w.location = hinge_point(kind, mirror, side)
            w.parent = root
            wing_objs.append((w, kind, side))

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


def build_shadow(name, centre, rx, ry, coll, mat):
    nx, ny = 24, 24
    verts, faces = [], []
    for i in range(ny + 1):
        for j in range(nx + 1):
            u = j / float(nx)
            v = i / float(ny)
            verts.append((centre[0] + (u - 0.5) * rx, centre[1] + (v - 0.5) * ry, 0.0016))
    for i in range(ny):
        for j in range(nx):
            a = i * (nx + 1) + j
            faces.append([a, a + 1, a + nx + 2, a + nx + 1])
    ob = new_mesh_obj(name, verts, faces, coll, smooth=False)
    ob.data.materials.append(mat)
    return ob


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
        body=mat_body("BodyMat", "#93a878", rough=0.88, fur=0.55, scale=70.0),
        abdomen=mat_body("AbdomenMat", "#9aab7d", rough=0.90, fur=0.55, scale=85.0),
        leg=mat_body("LegMat", "#4a4030", rough=0.66, fur=0.20, scale=120.0),
        eye=mat_body("EyeMat", "#191a12", rough=0.35, fur=0.0, scale=10.0),
        bark=mat_background(),
        shadow=mat_shadow(),
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
    moth_l = build_moth("Moth_L", -1, coll, mats)
    moth_l["root"].location = Vector((fx(0.300), fy(0.520), 0.0))
    moth_l["root"].rotation_euler = (0.0, 0.0, math.radians(50.0))
    moth_l["root"].scale = (0.55 * S,) * 3

    moth_r = build_moth("Moth_R", 1, coll, mats)
    moth_r["root"].location = Vector((fx(0.700), fy(0.490), 0.0))
    moth_r["root"].rotation_euler = (0.0, 0.0, math.radians(-60.0))
    moth_r["root"].scale = (0.60 * S,) * 3

    sh_l = build_shadow("Shadow_L", (fx(0.300) + 0.095, fy(0.520) - 0.085), 1.30, 1.02, coll, mats["shadow"])
    sh_r = build_shadow("Shadow_R", (fx(0.700) + 0.095, fy(0.490) - 0.085), 1.34, 1.06, coll, mats["shadow"])

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
    # light travels toward +X/-Y/-Z => it comes from the upper left of frame
    sun.rotation_euler = (math.radians(-38.0), 0.0, math.radians(40.0))

    world = bpy.data.worlds.new("World")
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    bg.inputs["Color"].default_value = (0.66, 0.70, 0.64, 1.0)
    bg.inputs["Strength"].default_value = AMBIENT
    scene.world = world

    # ================================================ ANIMATION
    fps = FPS

    def key_all(ob, frames, interp=None):
        for f in frames:
            ob.keyframe_insert("rotation_euler", index=1, frame=f)

    def set_interp(ob, interp):
        ad = ob.animation_data
        if not ad or not ad.action:
            return
        for fc in ad.action.fcurves:
            for kp in fc.keyframe_points:
                kp.interpolation = interp

    def flap_sign(mir, side):
        return -(mir * side)

    def set_wing(w, kind, side, mir, up_deg):
        s = flap_sign(mir, side)
        a = up_deg if kind == "fore" else up_deg * 0.92
        w.rotation_euler[1] = math.radians(s * a)

    # --- LEFT moth
    L = moth_l["root"]
    REST_L = 16.0
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
            lag = 0.0 if kind == "fore" else -0.9
            up = 8.0 * ramp + amp * math.sin(ph + lag)
            set_wing(w, kind, side, -1, up)
            w.keyframe_insert("rotation_euler", index=1, frame=f)
            w.rotation_euler[0] = math.radians(7.0 * math.cos(ph + lag) * ramp)
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

    m_l = mats["shadow"].copy()
    m_l.name = "ShadowDecal_L"
    sh_l.data.materials[0] = m_l
    op = m_l.node_tree.nodes["opacity"]
    for f, v in [(0, 0.95), (9, 0.95), (14, 0.60), (22, 0.24), (30, 0.07), (47, 0.0)]:
        op.outputs[0].default_value = v
        op.outputs[0].keyframe_insert("default_value", frame=f)
    for f, (dx, dy, sc) in {0: (0.0, 0.0, 1.0), 9: (0.0, 0.0, 1.0),
                            16: (-0.18, -0.10, 1.10), 28: (-0.45, -0.26, 1.22),
                            47: (-0.70, -0.42, 1.26)}.items():
        sh_l.location = Vector((dx, dy, 0.0))
        sh_l.keyframe_insert("location", frame=f)
        sh_l.scale = (sc, sc, 1.0)
        sh_l.keyframe_insert("scale", frame=f)

    # --- RIGHT moth: stays put, two slow open/close cycles (+-8 deg)
    R = moth_r["root"]
    for w, kind, side in moth_r["wings"]:
        set_wing(w, kind, side, 1, 0.0)
    for f, up in [(0, 0.0), (5, 8.0), (13, 0.0), (22, 0.0), (27, 8.0), (36, 0.0), (47, 0.0)]:
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
                  flap_hz=FLAP_HZ, flap_amp_deg=AMP)
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
