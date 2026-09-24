#!/usr/bin/env python3
"""TPS62933 board thermal 3x3 mesh (gmsh).
Geometry variants (copper keep-out regions, sourced from fea4 geo_meta):
  fullcu  : single-material board (no keep-out split)      -> 1 board body
  topcut  : top-layer copper under L2&L1 removed           -> board + 2 keep-out bodies
  dualcut : top+bottom copper under L2&L1 removed          -> board + 2 keep-out bodies
The keep-out bodies are full-thickness regions where the effective in-plane
conductivity is reduced (material assigned in make_case_3x3.py).
Board = anisotropic equivalent bulk (30um Cu + FR4). Components are 3D BODIES
with VOLUMETRIC heat source. Elmer needs 'Volumetric Heat Source'.
Usage: build_mesh_3x3.py <lcmax_mm> <lcmin_mm> <geo> <outname>
"""
import gmsh, sys, json
MM = 1e-3
lcmax = float(sys.argv[1]) * MM
lcmin = float(sys.argv[2]) * MM
geo   = sys.argv[3]
out   = sys.argv[4]

Lx, Ly, Lz = 26.50 * MM, 50.50 * MM, 1.49 * MM
z0 = -Lz / 2.0
# components: name, footprint_w(x), footprint_l(y), cx, cy, height(mm)
comps = [("U21", 2.10, 1.60,  8.25, 26.88, 0.60),
         ("D1",  4.30, 2.60, 14.10,  6.56, 2.20),
         ("L2",  6.50, 6.50,  6.22, 33.74, 3.00),
         ("L1",  4.20, 4.20, 16.38, 18.75, 2.00)]
terminals = [("TERM_IN",  3.81, 3.81,  4.95, 46.69),
             ("TERM_OUT", 3.81, 3.81,  5.00,  4.76)]
# copper keep-out rectangles (mm, board frame) derived from fea4 geo_meta.json
OX, OY = 6.223, 33.736   # fea4 origin (L2 centre) -> board frame offset (mm)
KEEPOUT = [("L2", OX-3.2999934, OY-3.499993,  OX+3.2999934, OY+3.499993),
           ("L1", OX+7.9838042, OY-17.272,    OX+12.3361958, OY-12.7098044)]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("board")
board = gmsh.model.occ.addBox(0, 0, z0, Lx, Ly, Lz)
extra = []
for (n, w, l, cx, cy, hh) in comps:
    extra.append(((cx-w/2)*MM, (cy-l/2)*MM, Lz/2.0, w*MM, l*MM, hh*MM))
ko_boxes = []
if geo != "fullcu":
    for (n, x0, y0, x1, y1) in KEEPOUT:
        extra.append((x0*MM, y0*MM, z0, (x1-x0)*MM, (y1-y0)*MM, Lz))
        ko_boxes.append((x0, y0, x1, y1))
comp_vols = [gmsh.model.occ.addBox(*b) for b in extra]
gmsh.model.occ.synchronize()
gmsh.model.occ.fragment([(3, board)], [(3, v) for v in comp_vols])
gmsh.model.occ.synchronize()

vols = [v for (d, v) in gmsh.model.getEntities(3)]
board_id = None
comp_ids = {}
ko_ids = []
comp_cent = {}
for v in vols:
    x, y, z = gmsh.model.occ.getCenterOfMass(3, v)
    X, Y = x*1e3, y*1e3
    if abs(z) < Lz/2.0 * 0.9:
        # board-type volume: board_main or keep-out
        isko = any(x0 <= X <= x1 and y0 <= Y <= y1 for (n, x0, y0, x1, y1) in KEEPOUT)
        if isko and geo != "fullcu":
            ko_ids.append(v)
        else:
            board_id = v
    else:
        comp_cent[v] = (x, y, z)
assert board_id is not None, "no board volume"
for (n, w, l, cx, cy, hh) in comps:
    best, bd = None, 1e18
    for v, (x, y, z) in comp_cent.items():
        d = (x - cx*MM)**2 + (y - cy*MM)**2
        if d < bd:
            bd, best = d, v
    comp_ids[n] = best
assert len(set(comp_ids.values())) == 4, comp_ids

# outer wetted surfaces = surfaces appearing once
from collections import Counter
allv = [board_id] + ko_ids + list(comp_ids.values())
cnt = Counter()
for v in allv:
    for (d, s) in gmsh.model.getBoundary([(3, v)], oriented=False):
        cnt[abs(s)] += 1
outer = sorted([s for s in cnt if cnt[s] == 1])

# terminal rects fragment on board top
tpad = [(2, gmsh.model.occ.addRectangle((cx-w/2)*MM, (cy-l/2)*MM, Lz/2.0, w*MM, l*MM))
        for (n, w, l, cx, cy) in terminals]
gmsh.model.occ.fragment([(3, board_id)], tpad)
gmsh.model.occ.synchronize()
allv = [board_id] + ko_ids + list(comp_ids.values())
cnt = Counter()
for v in allv:
    for (d, s) in gmsh.model.getBoundary([(3, v)], oriented=False):
        cnt[abs(s)] += 1
outer = sorted([s for s in cnt if cnt[s] == 1])
term_tags = {}
for (n, w, l, cx, cy) in terminals:
    best, bd = None, 1e18
    for s in outer:
        x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
        if abs(z - Lz/2.0) > 1e-9:
            continue
        d = (x - cx*MM)**2 + (y - cy*MM)**2
        if d < bd:
            bd, best = d, s
    term_tags[n] = best
assert len(set(term_tags.values())) == 2, term_tags
outer = [s for s in outer if s not in set(term_tags.values())]

gmsh.model.addPhysicalGroup(3, [board_id], 1)
for i, v in enumerate(ko_ids):
    gmsh.model.addPhysicalGroup(3, [v], 20 + i)
for i, (n, w, l, cx, cy, hh) in enumerate(comps):
    gmsh.model.addPhysicalGroup(3, [comp_ids[n]], 11 + i)
gmsh.model.addPhysicalGroup(2, outer, 1)
for i, (n, w, l, cx, cy) in enumerate(terminals):
    gmsh.model.addPhysicalGroup(2, [term_tags[n]], 200 + i)

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", outer)
gmsh.model.mesh.field.setNumber(1, "Sampling", 60)
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", lcmin)
gmsh.model.mesh.field.setNumber(2, "SizeMax", lcmax)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15e-3)
gmsh.model.mesh.field.setNumber(2, "DistMax", 3.0e-3)
gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)
gmsh.write(out + ".msh")
nn = len(gmsh.model.mesh.getNodes()[0])
nt = sum(len(t) for t in gmsh.model.mesh.getElements(3)[1])
meta = {"geo": geo, "lcmax_mm": lcmax*1e3, "lcmin_mm": lcmin*1e3, "nodes": nn, "tets": nt,
        "n_keepout": len(ko_ids), "keepout_rects_mm": KEEPOUT}
json.dump(meta, open(out + "_meta.json", "w"), indent=1)
print("MESH", out, meta)
gmsh.finalize()
