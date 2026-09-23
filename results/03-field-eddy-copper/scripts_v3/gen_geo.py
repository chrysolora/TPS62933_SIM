import sys, math
sys.path.insert(0, ".")
from parse_lib import load_pcbs, parse_path

F = "/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
MIL = 25.4e-6          # mil -> m
L2 = (245.0, -660.0)   # L2 centre (mil)

pcbs = load_pcbs(F)
d = pcbs[0]

def bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), max(xs), min(ys), max(ys))

# board outline
board = None
for t, pj in d["items"]:
    if t == "POLY" and pj.get("layerId") == 11:
        board = bbox(parse_path(pj["path"]))
bx0, bx1, by0, by1 = board
print("board mil x[%.2f,%.2f] y[%.2f,%.2f]" % board)

# top-layer pours
top = []
for t, pj in d["items"]:
    if t == "POUR" and pj.get("layerId") == 1:
        pts = parse_path(pj["path"])
        if not pts: continue
        bb = bbox(pts)
        # clip to board
        cb = (max(bb[0], bx0), min(bb[1], bx1), max(bb[2], by0), min(bb[3], by1))
        top.append((pj.get("name"), pj.get("netName"), cb))
        print("  %-9s %-22s clipped x[%.1f,%.1f] y[%.1f,%.1f]" % (pj.get("name"), pj.get("netName"), *cb))

gnd = [p for p in top if p[1] == "GND"]
print("GND pours covering board:", [(p[0], p[2]) for p in gnd])
voids = [p for p in top if p[1] != "GND"]

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

# void rects (inset 0.05mm inside board to avoid coincident faces)
ins = 0.05e-3
vd = []
seen = set()
for name, net, cb in voids:
    x0 = max(sx(cb[0]), bX0 + ins); x1 = min(sx(cb[1]), bX1 - ins)
    y0 = max(sy(cb[2]), bY0 + ins); y1 = min(sy(cb[3]), bY1 - ins)
    if x1 - x0 < 0.1e-3 or y1 - y0 < 0.1e-3: continue
    key = (round(x0, 9), round(y0, 9), round(x1, 9), round(y1, 9))
    if key in seen: continue
    seen.add(key)
    vd.append((name, x0, y0, x1, y1))
print("voids:", len(vd))

tcu = 0.035e-3
geo = []
geo.append('SetFactory("OpenCASCADE");')
geo.append('// Auto-generated real-geometry model (board-fitted). Units m, origin at L2 centre.')
geo.append('lc_air=5.0e-3; lc_cu_far=0.6e-3; lc_cu_near=0.3e-3; lc_core=0.6e-3; lc_coil=0.4e-3;')
geo.append('// Air domain (board bbox + 6 mm margin)')
geo.append('Box(1) = {%r, %r, %r, %r, %r, %r};' % (ax0, ay0, az0, ax1-ax0, ay1-ay0, az1-az0))
geo.append('// Copper slab = board footprint, 35 um thick at z=[-35,0]um')
geo.append('Box(2) = {%r, %r, %r, %r, %r, %r};' % (bX0, bY0, -tcu, bX1-bX0, bY1-bY0, tcu))
# void boxes
vids = []
for i, (name, x0, y0, x1, y1) in enumerate(vd):
    vid = 100 + i
    vids.append(vid)
    geo.append('Box(%d) = {%r, %r, %r, %r, %r, %r}; // void %s' % (vid, x0, y0, -0.5e-3, x1-x0, y1-y0, 1.0e-3, name))
geo.append('BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{%s}; Delete; };' % ",".join(str(v) for v in vids))
geo.append('Physical Volume("Copper") = {3};')
# core + coil
geo.append('Box(200) = {%r, %r, %r, %r, %r, %r};' % (-3.25e-3, -3.25e-3, 0.2e-3, 6.5e-3, 6.5e-3, 3.0e-3))
geo.append('Physical Volume("Core") = {200};')
geo.append('Box(210) = {%r, %r, %r, %r, %r, %r};' % (-2.6e-3, -2.6e-3, 0.9e-3, 5.2e-3, 5.2e-3, 0.6e-3))
geo.append('Box(211) = {%r, %r, %r, %r, %r, %r};' % (-2.0e-3, -2.0e-3, 0.895e-3, 4.0e-3, 4.0e-3, 0.61e-3))
geo.append('BooleanDifference(212) = { Volume{210}; Delete; }{ Volume{211}; Delete; };')
geo.append('Physical Volume("Coil") = {212};')
geo.append('Physical Volume("Air") = {1};')
geo.append('BooleanFragments{ Volume{1,3,200,212}; Delete; }{ }')
# outer surfaces
e = 0.2e-3
geo.append('Physical Surface("Outer") = {')
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax0+e, ay1+e, az1+e))
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax1-e, ay0-e, az0-e, ax1+e, ay1+e, az1+e))
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax1+e, ay0+e, az1+e))
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay1-e, az0-e, ax1+e, ay1+e, az1+e))
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r},' % (ax0-e, ay0-e, az0-e, ax1+e, ay1+e, az0+e))
geo.append('  Surface In BoundingBox{%r,%r,%r,%r,%r,%r}};' % (ax0-e, ay0-e, az1-e, ax1+e, ay1+e, az1+e))
# refinement fields
geo.append('Field[1]=Box; Field[1].VIn=lc_cu_near; Field[1].VOut=lc_air;')
geo.append('Field[1].XMin=%r; Field[1].XMax=%r; Field[1].YMin=%r; Field[1].YMax=%r; Field[1].ZMin=-1e-3; Field[1].ZMax=1e-3;' % (-10e-3, 10e-3, -10e-3, 10e-3))
geo.append('Field[2]=Box; Field[2].VIn=lc_cu_far; Field[2].VOut=lc_air;')
geo.append('Field[2].XMin=%r; Field[2].XMax=%r; Field[2].YMin=%r; Field[2].YMax=%r; Field[2].ZMin=-0.5e-3; Field[2].ZMax=0.5e-3;' % (bX0, bX1, bY0, bY1))
geo.append('Field[3]=Box; Field[3].VIn=lc_core; Field[3].VOut=lc_air; Field[3].XMin=-4e-3; Field[3].XMax=4e-3; Field[3].YMin=-4e-3; Field[3].YMax=4e-3; Field[3].ZMin=0.0; Field[3].ZMax=3.5e-3;')
geo.append('Field[4]=Min; Field[4].FieldsList={1,2,3}; Background Field=4;')
geo.append('Mesh.MshFileVersion=2.2; Mesh.SaveAll=1;')
geo.append('Mesh.CharacteristicLengthMax=lc_air; Mesh.CharacteristicLengthMin=lc_coil;')

open("build.geo", "w").write("\n".join(geo) + "\n")
print("wrote build.geo with", len(vd), "voids")
print("board m x[%.5f,%.5f] y[%.5f,%.5f]" % (bX0, bX1, bY0, bY1))
print("air   m x[%.5f,%.5f] y[%.5f,%.5f] z[%.4f,%.4f]" % (ax0, ax1, ay0, ay1, az0, az1))
boardA = (bX1-bX0)*(bY1-bY0)
voidA = sum((x1-x0)*(y1-y0) for _, x0, y0, x1, y1 in vd)
print("copper area m^2: board %.4e - voids %.4e = %.4e (%.1f%% of board)" % (
    boardA, voidA, boardA-voidA, 100*(boardA-voidA)/boardA))
print("copper volume m^3 = %.4e" % ((boardA-voidA)*tcu))
