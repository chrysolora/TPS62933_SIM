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

p = [d for d in docs if d["type"] == "PCB"][0]
# print first COMPONENT and first few ATTR
n = 0
for t, pj in p["items"]:
    if t == "COMPONENT" and n < 2:
        print("COMPONENT:", json.dumps(pj, ensure_ascii=False)[:700]); n += 1
n = 0
for t, pj in p["items"]:
    if t == "ATTR" and n < 6:
        print("ATTR:", json.dumps(pj, ensure_ascii=False)[:300]); n += 1
n = 0
for t, pj in p["items"]:
    if t == "LAYER" and n < 20:
        print("LAYER:", pj.get("layerType"), pj.get("layerName")); n += 1
n = 0
for t, pj in p["items"]:
    if t == "ELE_PLACEHOLDER" and n < 2:
        print("ELEPH:", json.dumps(pj, ensure_ascii=False)[:400]); n += 1
