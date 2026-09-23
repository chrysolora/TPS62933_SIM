import json

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

def flatten(seq):
    flat = []
    def fl(s):
        for e in s:
            if isinstance(e, list):
                fl(e)
            else:
                flat.append(e)
    fl(seq)
    return flat

def parse_path(seq):
    seq = flatten(seq)
    pts = []
    i = 0
    n = len(seq)
    while i < n:
        tok = seq[i]
        if isinstance(tok, str):
            if tok == "ARC":
                i += 1
                if i < n and isinstance(seq[i], (int, float)):
                    i += 1
                continue
            i += 1
            continue
        else:
            if i + 1 < n and isinstance(seq[i+1], (int, float)):
                pts.append((float(seq[i]), float(seq[i+1])))
                i += 2
            else:
                i += 1
    return pts

L2 = (245.0, -660.0)
# L2 6.5x6.5mm footprint box
H = 3.25
LB = (L2[0]-H, L2[1]-H, L2[0]+H, L2[1]+H)

for d in pcbs:
    print("##### PCB %s" % d["uuid"])
    hits = []
    for t, pj in d["items"]:
        if t == "POURED":
            pts = []
            for sub in pj.get("pourFill", []):
                pts += parse_path(sub.get("path", []))
            if not pts:
                continue
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            bx = (min(xs), min(ys), max(xs), max(ys))
            # overlap with L2 box?
            if not (bx[2] < LB[0] or bx[0] > LB[2] or bx[3] < LB[1] or bx[1] > LB[3]):
                hits.append((round(min(xs),1), round(min(ys),1), round(max(xs),1), round(max(ys),1), len(pts)))
    print("  POURED polygons overlapping L2 box(%s): %d" % (str(LB), len(hits)))
    for hh in hits[:20]:
        print("     bbox x[%.1f,%.1f] y[%.1f,%.1f] npts=%d" % hh)
