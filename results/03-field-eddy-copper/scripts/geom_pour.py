import json, re

f = "/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
lines = open(f, encoding="utf-8", errors="replace").read().split("\n")
docs = []
cur = None
for l in lines:
    if not l.strip() or "||" not in l:
        continue
    h, p = l.split("||", 1)
    if p.endswith("|"):
        p = p[:-1]
    try:
        hj = json.loads(h)
    except Exception:
        continue
    t = hj.get("type")
    if t == "DOCHEAD":
        try:
            pj = json.loads(p)
        except Exception:
            pj = {}
        cur = {"type": pj.get("docType"), "uuid": pj.get("uuid"), "items": []}
        docs.append(cur)
    elif cur is not None:
        try:
            pj = json.loads(p)
        except Exception:
            pj = {"_raw": p[:300]}
        cur["items"].append((t, pj))

pcbs = [d for d in docs if d["type"] == "PCB"]

def parse_path(seq):
    """seq is the flat path list. Return list of (x,y) vertex points.
    Pattern: [x0,y0,'L',x1,y1,'ARC',angle,cx,cy,'L',...] or [x0,y0, x1,y1, ...]"""
    # flatten nested lists
    flat = []
    def fl(s):
        for e in s:
            if isinstance(e, list):
                fl(e)
            else:
                flat.append(e)
    fl(seq)
    seq = flat
    pts = []
    i = 0
    n = len(seq)
    # first two numbers
    while i < n:
        tok = seq[i]
        if isinstance(tok, str):
            if tok == "ARC":
                # next is angle, then a coordinate pair; skip angle
                i += 1
                if i < n and isinstance(seq[i], (int, float)):
                    i += 1  # skip angle
                continue
            else:
                i += 1
                continue
        else:
            # coordinate pair
            if i + 1 < n and isinstance(seq[i+1], (int, float)):
                pts.append((float(seq[i]), float(seq[i+1])))
                i += 2
            else:
                i += 1
    return pts

def all_pts(pj):
    out = []
    if "path" in pj:
        out += parse_path(pj["path"])
    if "pourFill" in pj:
        for sub in pj["pourFill"]:
            if "path" in sub:
                out += parse_path(sub["path"])
    return out

L2 = (245.0, -660.0)

for d in pcbs:
    print("##### PCB %s" % d["uuid"])
    for t, pj in d["items"]:
        if t == "POUR":
            pts = all_pts(pj)
            if not pts:
                continue
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            cover = (min(xs) <= L2[0] <= max(xs)) and (min(ys) <= L2[1] <= max(ys))
            print("  POUR %-10s layer=%s net=%-20s bbox x[%.1f,%.1f] y[%.1f,%.1f] coversL2=%s" % (
                pj.get("name"), pj.get("layerId"), pj.get("netName"),
                min(xs), max(xs), min(ys), max(ys), cover))
