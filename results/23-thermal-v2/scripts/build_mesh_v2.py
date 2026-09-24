#!/usr/bin/env python3
"""TPS62933 board thermal mesh v2 (gmsh).
FIX vs v1: components are 3D BODIES with VOLUMETRIC heat generation, not
uniform surface-flux patches -> removes the log-singularity at flux-patch edges.
Board = anisotropic equivalent bulk. Component bodies carry package theta.
Physical groups (tags are read by make_case_v2.py from the meta json):
  volume 1 = board ; 11..14 = U21,D1,L2,L1 bodies
  surface 1 = all OUTER wetted surfaces (board outer + comp outer)
  surface 2.. = terminal pads (extra wire-fin path)
Usage: build_mesh_v2.py <lcmax_mm> <lcmin_mm> <outname>
"""
import gmsh, sys, json
MM = 1e-3
lcmax = (float(sys.argv[1]) if len(sys.argv) > 1 else 1.2) * MM
lcmin = (float(sys.argv[2]) if len(sys.argv) > 2 else 0.15) * MM
out   = sys.argv[3] if len(sys.argv) > 3 else "mesh_v2"
Lx, Ly, Lz = 26.50 * MM, 50.50 * MM, 1.49 * MM
z0 = -Lz / 2.0
# name, footprint_w(x), footprint_l(y), cx, cy, height
comps = [("U21", 2.10, 1.60,  8.25, 26.88, 0.60),
         ("D1",  4.30, 2.60, 14.10,  6.56, 2.20),
         ("L2",  6.50, 6.50,  6.22, 33.74, 3.00),
         ("L1",  4.20, 4.20, 16.38, 18.75, 2.00)]
# terminal pads (board-facts §3 J x2), small rect on top face for wire-fin path
terminals = [("TERM_IN",  3.81, 3.81,  4.95, 46.69),
             ("TERM_OUT", 3.81, 3.81,  5.00,  4.76)]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("board")
board = gmsh.model.occ.addBox(0, 0, z0, Lx, Ly, Lz)
comp_vols = []
for (n, w, l, cx, cy, hh) in comps:
    v = gmsh.model.occ.addBox((cx - w/2)*MM, (cy - l/2)*MM, Lz/2.0, w*MM, l*MM, hh*MM)
    comp_vols.append(v)
gmsh.model.occ.synchronize()
# fragment board with comp bodies so interfaces are conformal
gmsh.model.occ.fragment([(3, board)], [(3, v) for v in comp_vols])
gmsh.model.occ.synchronize()

vols = [v for (d, v) in gmsh.model.getEntities(3)]
assert len(vols) == 5, vols
# identify board vs comp volumes by centroid z (board centred at z=0)
board_id = None
comp_ids = {}
comp_cent = {}
for v in vols:
    x, y, z = gmsh.model.occ.getCenterOfMass(3, v)
    if abs(z) < Lz/2.0 * 0.9:
        board_id = v
    else:
        comp_cent[v] = (x, y, z)
assert board_id is not None
for (n, w, l, cx, cy, hh) in comps:
    best, bd = None, 1e18
    for v, (x, y, z) in comp_cent.items():
        d = (x - cx*MM)**2 + (y - cy*MM)**2
        if d < bd:
            bd, best = d, v
    comp_ids[n] = best
assert len(set(comp_ids.values())) == 4, comp_ids

# ---- outer wetted surfaces = surfaces of the 5 volumes that are NOT internal
all_surf = set()
for v in vols:
    for (d, s) in gmsh.model.getBoundary([(3, v)], oriented=False):
        all_surf.add(abs(s))
# internal = shared between two volumes -> appears twice
from collections import Counter
cnt = Counter()
for v in vols:
    for (d, s) in gmsh.model.getBoundary([(3, v)], oriented=False):
        cnt[abs(s)] += 1
outer = sorted([s for s in all_surf if cnt[s] == 1])

# top-face board surfaces for terminals (imprint terminal rects)
top_surfs = []
for s in outer:
    x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
    if abs(z - Lz/2.0) < 1e-9 and 0 < x < Lx and 0 < y < Ly:
        top_surfs.append((s, x, y))
# imprint terminals by fragmenting a thin disc? simpler: pick nearest existing
# top surface centre -> but board top is one big face; use occ fragment with rects
tpad = []
for (n, w, l, cx, cy) in terminals:
    tpad.append((2, gmsh.model.occ.addRectangle((cx - w/2)*MM, (cy - l/2)*MM, Lz/2.0, w*MM, l*MM)))
gmsh.model.occ.fragment([(3, board_id)], tpad)
gmsh.model.occ.synchronize()
vols2 = [v for (d, v) in gmsh.model.getEntities(3)]
# board_id may change; re-find
board_id = None
for v in vols2:
    x, y, z = gmsh.model.occ.getCenterOfMass(3, v)
    if abs(z) < Lz/2.0 * 0.9:
        board_id = v
assert board_id is not None
# re-collect ALL outer surfaces after fragment
cnt = Counter()
for v in vols2:
    for (d, s) in gmsh.model.getBoundary([(3, v)], oriented=False):
        cnt[abs(s)] += 1
outer = sorted([s for s in cnt if cnt[s] == 1])
# terminal surfaces = small top faces with centre near terminal centres
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

# terminal surfaces must NOT be in the generic outer-convection group
outer = [s for s in outer if s not in set(term_tags.values())]
assert len(outer) > 0

gmsh.model.addPhysicalGroup(3, [board_id], 1)
gmsh.model.addPhysicalGroup(2, outer, 1)
for i, (n, w, l, cx, cy, hh) in enumerate(comps):
    gmsh.model.addPhysicalGroup(3, [comp_ids[n]], 11 + i)
for i, (n, w, l, cx, cy) in enumerate(terminals):
    gmsh.model.addPhysicalGroup(2, [term_tags[n]], 200 + i)

# size field: refine near component bodies
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
tets = gmsh.model.mesh.getElements(3)[1]
nt = sum(len(t) for t in tets)
# volume + exterior-area of each device body
vol_body = {}
area_body = {}
for i, (n, w, l, cx, cy, hh) in enumerate(comps):
    vol_body[n] = w*MM*l*MM*hh*MM
    # exterior area of a box body: 2*(wl+lh+wh) minus bottom (bonded) = wl+2*(l*hh)+2*(w*hh)
    area_body[n] = (w*MM*l*MM) + 2*(l*MM*hh*MM) + 2*(w*MM*hh*MM)
meta = {"lcmax_mm": lcmax*1e3, "lcmin_mm": lcmin*1e3, "nodes": nn, "tets": nt,
        "board_id": 1, "comp_body_tags": {comps[i][0]: 11+i for i in range(len(comps))},
        "term_surf_tags": term_tags, "outer_surf_tag": 1,
        "comp_vol_m3": vol_body, "comp_extarea_m2": area_body,
        "footprint_area_m2": {comps[i][0]: comps[i][1]*MM*comps[i][2]*MM for i in range(len(comps))},
        "comp_heights_mm": {comps[i][0]: comps[i][5] for i in range(len(comps))}}
json.dump(meta, open(out + "_meta.json", "w"), indent=1)
print("MESH", meta)
gmsh.finalize()
