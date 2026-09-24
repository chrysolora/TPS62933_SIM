import sys, math, os, json
sys.path.insert(0, "/mnt/raid10/sim-work/tps62933/fea3")
from parse_lib import load_pcbs, parse_path

IDX = int(sys.argv[1]); OUT = sys.argv[2]
F = "/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
MIL = 25.4e-6
L2 = (245.0, -660.0)
TB = 1.0e-3          # board thickness (ASSUMPTION, see report)
tcu = 0.035e-3

pcbs = load_pcbs(F); d = pcbs[IDX]

def bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), max(xs), min(ys), max(ys))

board = None
for t, pj in d["items"]:
    if t == "POLY" and pj.get("layerId") == 11:
        board = bbox(parse_path(pj["path"]))
bx0, bx1, by0, by1 = board
print("board mil x[%.2f,%.2f] y[%.2f,%.2f]" % board)

# --- top non-GND pours = holes in top GND plane (fea3 method) ---
tpours = []
for t, pj in d["items"]:
    if t == "POUR" and pj.get("layerId") == 1:
        pts = parse_path(pj["path"])
        if not pts: continue
        bb = bbox(pts)
        cb = (max(bb[0], bx0), min(bb[1], bx1), max(bb[2], by0), min(bb[3], by1))
        tpours.append((pj.get("name"), pj.get("netName"), cb))
print("top pours:", [(p[0], p[1]) for p in tpours])

# --- REGION(COPPER) keepouts by layer scope ---
def scope(lid):
    if lid == 1: return ("top",)
    if lid == 2: return ("bot",)
    if lid == 12: return ("top", "bot")   # MULTI
    return ()

rtop, rbot = [], []
for t, pj in d["items"]:
    if not t.upper().startswith("REGION"): continue
    if "COPPER" not in (pj.get("prohibitType") or []): continue
    pts = parse_path(pj.get("path"))
    if not pts: continue
    bb = bbox(pts); sc = scope(pj.get("layerId"))
    print("REGION layer=%s scope=%s bbox x[%.1f,%.1f] y[%.1f,%.1f]" % (pj.get("layerId"), sc, *bb))
    if "top" in sc: rtop.append((pj.get("name"), bb))
    if "bot" in sc: rbot.append((pj.get("name"), bb))

# convert to metres, origin at L2
def m(v): return v * MIL
LX, LY = m(L2[0]), m(L2[1])
def sx(v): return m(v) - LX
def sy(v): return m(v) - LY
bX0, bX1 = sx(bx0), sx(bx1)
bY0, bY1 = sy(by0), sy(by1)

MARGIN = 6e-3
ax0, ax1 = bX0 - MARGIN, bX1 + MARGIN
ay0, ay1 = bY0 - MARGIN, bY1 + MARGIN
az0, az1 = -8e-3, 12e-3
ins = 0.05e-3

def clip_box(cb):
    x0 = max(sx(cb[0]), bX0 + ins); x1 = min(sx(cb[1]), bX1 - ins)
    y0 = max(sy(cb[2]), bY0 + ins); y1 = min(sy(cb[3]), bY1 - ins)
    if x1 - x0 < 0.05e-3 or y1 - y0 < 0.05e-3: return None
    return (x0, y0, x1, y1)

def dedup(boxes):
    out = []; seen = set()
    for name, cb in boxes:
        b = clip_box(cb)
        if b is None: continue
        key = tuple(round(v, 9) for v in b)
        if key in seen: continue
        seen.add(key); out.append((name, *b))
    return out

top_voids = dedup([(nm, cb) for nm, net, cb in tpours if net != "GND"])
top_reg = dedup([("REG_" + str(nm), cb) for nm, cb in rtop])
bot_reg = dedup([("REG_" + str(nm), cb) for nm, cb in rbot])

boardA = (bX1 - bX0) * (bY1 - bY0)
def area(vs): return sum((x1 - x0) * (y1 - y0) for _, x0, y0, x1, y1 in vs)
print("top: pours_voids=%d regions=%d ; bot: regions=%d" % (len(top_voids), len(top_reg), len(bot_reg)))
print("board area %.4e m^2" % boardA)
print("top void area %.4e  (%.1f%%) ; top region %.4e ; bot region %.4e" % (
    area(top_voids), 100 * area(top_voids) / boardA, area(top_reg), area(bot_reg)))
top_cu = boardA - area(top_voids) - area(top_reg)
bot_cu = boardA - area(bot_reg)
print("top copper %.4e m^2 (%.1f%%) ; bot copper %.4e m^2 (%.1f%%)" % (
    top_cu, 100 * top_cu / boardA, bot_cu, 100 * bot_cu / boardA))

g = []
g.append('SetFactory("OpenCASCADE");')
g.append('// fea4 variant idx=%d ; units m, origin at L2 centre' % IDX)
g.append('lc_air=6.0e-3; lc_cu_far=0.9e-3; lc_cu_near=0.4e-3; lc_core=0.7e-3; lc_coil=0.5e-3;')
# air
g.append('Box(1) = {%r, %r, %r, %r, %r, %r};' % (ax0, ay0, az0, ax1 - ax0, ay1 - ay0, az1 - az0))
# copper top slab z=[-tcu,0]
g.append('Box(2) = {%r, %r, %r, %r, %r, %r};' % (bX0, bY0, -tcu, bX1 - bX0, bY1 - bY0, tcu))
# copper bottom slab z=[-TB-tcu,-TB]
g.append('Box(3) = {%r, %r, %r, %r, %r, %r};' % (bX0, bY0, -TB - tcu, bX1 - bX0, bY1 - bY0, tcu))

vid = 100
top_ids = []
for nm, x0, y0, x1, y1 in top_voids + top_reg:
    g.append('Box(%d) = {%r, %r, %r, %r, %r, %r};' % (vid, x0, y0, -0.5e-3, x1 - x0, y1 - y0, 1.0e-3))
    top_ids.append(vid); vid += 1
bot_ids = []
for nm, x0, y0, x1, y1 in bot_reg:
    g.append('Box(%d) = {%r, %r, %r, %r, %r, %r};' % (vid, x0, y0, -TB - 0.5e-3, x1 - x0, y1 - y0, 1.0e-3))
    bot_ids.append(vid); vid += 1

ct, cb = 2, 3
if top_ids:
    g.append('BooleanDifference(10) = { Volume{2}; Delete; }{ Volume{%s}; Delete; };' % ",".join(map(str, top_ids)))
    ct = 10
if bot_ids:
    g.append('BooleanDifference(11) = { Volume{3}; Delete; }{ Volume{%s}; Delete; };' % ",".join(map(str, bot_ids)))
    cb = 11

g.append('Physical Volume("CopperTop") = {%d};' % ct)
g.append('Physical Volume("CopperBot") = {%d};' % cb)
g.append('Box(200) = {%r, %r, 0.2e-3, 6.5e-3, 6.5e-3, 3.0e-3};' % (-3.25e-3, -3.25e-3))
g.append('Physical Volume("Core") = {200};')
g.append('Box(210) = {%r, %r, 0.9e-3, 5.2e-3, 5.2e-3, 0.6e-3};' % (-2.6e-3, -2.6e-3))
g.append('Box(211) = {%r, %r, 0.895e-3, 4.0e-3, 4.0e-3, 0.61e-3};' % (-2.0e-3, -2.0e-3))
g.append('BooleanDifference(212) = { Volume{210}; Delete; }{ Volume{211}; Delete; };')
g.append('Physical Volume("Coil") = {212};')
g.append('Physical Volume("Air") = {1};')
g.append('BooleanFragments{ Volume{1,%d,%d,200,212}; Delete; }{ }' % (ct, cb))
e = 0.2e-3
g.append('Physical Surface("Outer") = {')
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax0+e, ay1+e, az1+e))
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax1-e, ay0-e, az0-e, ax1+e, ay1+e, az1+e))
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax1+e, ay0+e, az1+e))
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay1-e, az0-e, ax1+e, ay1+e, az1+e))
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax1+e, ay1+e, az0+e))
g.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r}};' % (ax0-e, ay0-e, az1-e, ax1+e, ay1+e, az1+e))
g.append('Field[1]=Box; Field[1].VIn=lc_cu_near; Field[1].VOut=lc_air;')
g.append('Field[1].XMin=-10e-3; Field[1].XMax=10e-3; Field[1].YMin=-10e-3; Field[1].YMax=10e-3; Field[1].ZMin=-1.2e-3; Field[1].ZMax=1e-3;')
g.append('Field[2]=Box; Field[2].VIn=lc_cu_far; Field[2].VOut=lc_air;')
g.append('Field[2].XMin=%r; Field[2].XMax=%r; Field[2].YMin=%r; Field[2].YMax=%r; Field[2].ZMin=%r; Field[2].ZMax=0.5e-3;' % (bX0, bX1, bY0, bY1, -TB - tcu))
g.append('Field[3]=Box; Field[3].VIn=lc_core; Field[3].VOut=lc_air; Field[3].XMin=-4e-3; Field[3].XMax=4e-3; Field[3].YMin=-4e-3; Field[3].YMax=4e-3; Field[3].ZMin=0.0; Field[3].ZMax=3.5e-3;')
g.append('Field[4]=Min; Field[4].FieldsList={1,2,3}; Background Field=4;')
g.append('Mesh.MshFileVersion=2.2; Mesh.SaveAll=1;')
g.append('Mesh.CharacteristicLengthMax=lc_air; Mesh.CharacteristicLengthMin=lc_coil;')

os.makedirs(OUT, exist_ok=True)
open(os.path.join(OUT, "build.geo"), "w").write("\n".join(g) + "\n")
json.dump({"idx": IDX, "boardA": boardA, "top_cu": top_cu, "bot_cu": bot_cu,
           "top_voids": top_voids, "top_reg": top_reg, "bot_reg": bot_reg,
           "TB": TB, "tcu": tcu}, open(os.path.join(OUT, "geo_meta.json"), "w"), indent=1)
print("wrote", os.path.join(OUT, "build.geo"))
