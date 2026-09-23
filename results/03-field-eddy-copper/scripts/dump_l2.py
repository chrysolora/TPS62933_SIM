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

# L2 device uuid
L2DEV = "548c68cdcac279f4"

def bbox(path):
    xs, ys = [], []
    def walk(pp):
        if isinstance(pp, list):
            for e in pp:
                walk(e)
        elif isinstance(pp, (int, float)):
            pass
    # path is list of segments; points are numbers at even-ish positions
    it = str(path)
    # crude: extract all numbers, pair them
    import re
    nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", it)]
    pts = list(zip(nums[0::2], nums[1::2]))
    if not pts:
        return None
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))

for d in pcbs:
    print("##### PCB uuid=%s" % d["uuid"])
    # build ATTR designator map for components
    for t, pj in d["items"]:
        if t == "COMPONENT":
            s = json.dumps(pj, ensure_ascii=False)
            if L2DEV in s or "L2" in s:
                print("COMPONENT:", s[:400])
    # device -> designator via COMPONENT? print keys
    break
