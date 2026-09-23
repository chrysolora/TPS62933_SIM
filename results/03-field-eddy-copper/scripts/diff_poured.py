import json, re, hashlib, sys

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
names = {p["uuid"]: p for p in pcbs}

def sigs(p):
    out = []
    for t, pj in p["items"]:
        if t in ("POURED", "POUR", "REGION"):
            s = json.dumps(pj, sort_keys=True)
            # canonicalize floats
            s = re.sub(r"-?\d+\.\d+", lambda m: ("%.4f" % float(m.group(0))), s)
            out.append((t, hashlib.md5(s.encode()).hexdigest()[:10], len(s), s[:180]))
    return out

for u, p in names.items():
    print("##### PCB %s  POURED/POUR/REGION count=%d" % (u, sum(1 for t,_ in p["items"] if t in ("POURED","POUR","REGION"))))

# compare region/poured sig sets
from collections import Counter
import itertools
uuids = list(names.keys())
for a, b in itertools.combinations(uuids, 2):
    ca = Counter(h for t,h,_,_ in sigs(names[a]))
    cb = Counter(h for t,h,_,_ in sigs(names[b]))
    only_a = ca - cb
    only_b = cb - ca
    print("### diff %s vs %s : only_in_A=%d only_in_B=%d" % (a[:8], b[:8], sum(only_a.values()), sum(only_b.values())))
    # print sample of differing entries with geometry size
    for src, cnt in (("A", only_a), ("B", only_b)):
        for h in list(cnt)[:6]:
            for t, hh, L, s in sigs(names[a if src=="A" else b]):
                if hh == h:
                    print("   %s %s len=%d %s" % (src, t, L, s))
                    break
