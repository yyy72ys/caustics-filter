import os, subprocess, json, hashlib
OUT = r"D:\yohey2026\dsh\muse\caustics-filter\offline\out\v9samples"

def probe(p):
    cmd = ["ffprobe","-v","error","-count_frames","-select_streams","v:0",
           "-show_entries","stream=codec_name,width,height,r_frame_rate,avg_frame_rate,nb_read_frames,duration,pix_fmt",
           "-of","json",p]
    return json.loads(subprocess.run(cmd,capture_output=True,text=True).stdout)["streams"][0]

def sha12(p):
    return hashlib.sha256(open(p,"rb").read()).hexdigest()[:12]

rep = {}
for tag, name in (("eta100","comp-eta100-1s.mp4"),("eta080","comp-eta080-1s.mp4")):
    p = os.path.join(OUT,name)
    s = probe(p)
    nb = int(s.get("nb_read_frames") or 0)
    dur = float(s.get("duration") or 0)
    w,h = int(s["width"]), int(s["height"])
    ok = (s["codec_name"]=="h264" and w==832 and h==480 and s["pix_fmt"]=="yuv420p"
          and s["avg_frame_rate"]=="24/1" and nb==24 and abs(dur-1.0) < 1e-6)
    rep[tag] = {"file": name, "bytes": os.path.getsize(p), "sha12": sha12(p),
                "ffprobe": s, "pass": bool(ok),
                "encoder": "h264_nvenc preset p6 tune hq rc vbr cq 17 yuv420p 24fps"}
    print(tag, nb, dur, w, h, s["codec_name"], s["pix_fmt"], "PASS" if ok else "FAIL")

json.dump(rep, open(os.path.join(OUT,"encode-report.json"),"w"), indent=2)
print("written encode-report.json")
