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
        cur = {"type": pj.get("docType"), "uuid": pj.get("uuid"), "items": []}
        docs.append(cur)
    elif cur is not None:
        try:
            pj = json.loads(p)
        except Exception:
            pj = {"_raw": p[:300]}
        cur["items"].append((t, pj))

pcbs = [d for d in docs if d["type"] == "PCB"]
what = sys.argv[1] if len(sys.argv) > 1 else "layer_phys"

for d in pcbs:
    print("##### PCB uuid=%s" % d["uuid"])
    for t, pj in d["items"]:
        if what == "layer_phys" and t == "LAYER_PHYS":
            print(json.dumps(pj, ensure_ascii=False))
        elif what == "pour":
            if t in ("POUR", "POURED"):
                s = json.dumps(pj, ensure_ascii=False)
                print(t, s[:600])
