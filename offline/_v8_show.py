"""Compact viewer for numbers-v8.json (chunk-4 reporting aid).

Usage (from offline/):
    python _v8_show.py out/v8/numbers-v8.json gates
    python _v8_show.py out/v8/numbers-v8.json prereq
    python _v8_show.py out/v8/numbers-v8.json keys
"""
import json
import sys


def walk(o, prefix="", depth=0, max_depth=3):
    if isinstance(o, dict):
        for k, v in o.items():
            p = prefix + "." + str(k) if prefix else str(k)
            if isinstance(v, (dict, list)) and depth < max_depth:
                if isinstance(v, list) and v and not isinstance(v[0], (dict, list)):
                    print("%-58s %s" % (p, str(v)[:110]))
                else:
                    print("%-58s <%s len=%d>" % (p, type(v).__name__, len(v)))
                    walk(v, p, depth + 1, max_depth)
            elif isinstance(v, dict):
                print("%-58s %s" % (p, json.dumps(v)[:110]))
            elif isinstance(v, list) and v and not isinstance(v[0], (dict, list)):
                print("%-58s %s" % (p, str(v)[:110]))
            else:
                print("%-58s %s" % (p, str(v)[:110]))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            walk(v, "%s[%d]" % (prefix, i), depth + 1, max_depth)


def main():
    path = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "keys"
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    if mode == "keys":
        walk(d, max_depth=2)
    elif mode == "gates":
        g = d.get("gates", d)
        for name, node in g.items():
            print("=" * 78)
            print(name, "PASS=", node.get("pass") if isinstance(node, dict) else None)
            walk(node, max_depth=2)
    elif mode == "prereq":
        p = d.get("prerequisite_tasks", d)
        walk(p, max_depth=2)
    elif mode == "all":
        walk(d, max_depth=3)
    else:
        node = d
        for part in mode.split("."):
            node = node[part]
        print(json.dumps(node, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
