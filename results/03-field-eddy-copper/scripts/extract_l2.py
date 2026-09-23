import json, re, sys

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

def num_extract(path):
    s = json.dumps(path)
    nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", s)]
    return nums

for d in pcbs:
    print("##### PCB uuid=%s" % d["uuid"])
    # designator -> parentId
    des = {}
    for t, pj in d["items"]:
        if t == "ATTR" and pj.get("key") == "Designator":
            des[pj.get("value")] = pj.get("parentId")
    # ELE_PLACEHOLDER positions keyed by id
    for t, pj in d["items"]:
        if t == "ELE_PLACEHOLDER":
            pid = pj.get("id")
            for dname, ppid in des.items():
                if ppid == pid:
                    print("  ELEPH for %s: %s" % (dname, json.dumps(pj, ensure_ascii=False)[:300]))
    # COMPONENT list with DeviceName
    for t, pj in d["items"]:
        if t == "COMPONENT":
            dn = pj["attrs"].get("DeviceName", "")
            if "ZEMS" in dn or "0650" in dn or "FXL" in dn:
                print("  COMP:", pj.get("x"), pj.get("y"), pj.get("angle"), pj.get("layerId"), dn[:120])
    # POUR bounding boxes on top layer
    for t, pj in d["items"]:
        if t == "POUR":
            path = pj.get("path")
            nums = num_extract(path)
            if not nums:
                continue
            xs = nums[0::2]; ys = nums[1::2]
            print("  POUR %-10s layer=%s net=%-22s bbox x[%.1f,%.1f] y[%.1f,%.1f]" % (
                pj.get("name"), pj.get("layerId"), pj.get("netName"),
                min(xs), max(xs), min(ys), max(ys)))
