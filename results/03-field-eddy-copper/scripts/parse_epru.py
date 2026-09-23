import json, sys, collections

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
        cur = {"type": pj.get("docType"), "uuid": pj.get("uuid"),
               "title": pj.get("title"), "items": []}
        docs.append(cur)
    elif cur is not None:
        try:
            pj = json.loads(p)
        except Exception:
            pj = {"_raw": p[:200]}
        cur["items"].append((t, pj))

cnt = collections.Counter(d["type"] for d in docs)
print("DOC TYPES:", dict(cnt))
for d in docs:
    itc = collections.Counter(t for t, _ in d["items"])
    print("%-12s uuid=%s items=%d %s" % (d["type"], d["uuid"], len(d["items"]), dict(itc)))
