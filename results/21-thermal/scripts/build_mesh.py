#!/usr/bin/env python3
"""Board-level thermal mesh (gmsh). Board only; component FOOTPRINT patches are
imprinted on the top face and exported as separate physical surfaces so the
component power can be applied as a surface heat flux. SI units (meters).
Usage: build_mesh.py <lc_max_mm> <lc_min_mm> <outname>
"""
import gmsh, sys, json

lc_max = (float(sys.argv[1]) if len(sys.argv) > 1 else 1.5) * 1e-3
lc_min = (float(sys.argv[2]) if len(sys.argv) > 2 else 0.4) * 1e-3
out    = sys.argv[3] if len(sys.argv) > 3 else "board_msh"

MM = 1e-3
Lx, Ly, Lz = 26.50 * MM, 50.50 * MM, 1.49 * MM
z0 = -Lz / 2.0
comps = [
    ("D1",  4.3, 2.6, 14.10,  6.56),
    ("U21", 2.9, 1.6,  8.25, 26.88),
    ("L2",  7.0, 6.6,  6.22, 33.74),
    ("L1",  4.4, 4.2, 16.38, 18.75),
]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("board")
board = gmsh.model.occ.addBox(0, 0, z0, Lx, Ly, Lz)
pads = []
for (n, w, l, cx, cy) in comps:
    pads.append((2, gmsh.model.occ.addRectangle((cx - w/2)*MM, (cy - l/2)*MM, Lz/2.0, w*MM, l*MM)))
gmsh.model.occ.synchronize()
gmsh.model.occ.fragment([(3, board)], pads)
gmsh.model.occ.synchronize()

allv = [v for (d, v) in gmsh.model.getEntities(3)]
assert len(allv) == 1, allv
board_id = allv[0]

# top-face surfaces (z = Lz/2)
bnd = [abs(t) for (d, t) in gmsh.model.getBoundary([(3, board_id)], oriented=False)]
top = []
for s in bnd:
    x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
    if abs(z - Lz/2.0) < 1e-9:
        top.append((s, x, y))

fp = {}   # comp index -> surface tag
for i, (n, w, l, cx, cy) in enumerate(comps):
    ex, ey = cx*MM, cy*MM
    best, bd = None, 1e9
    for (s, x, y) in top:
        d = (x-ex)**2 + (y-ey)**2
        if d < bd:
            bd, best = d, s
    fp[i] = best
fpset = set(fp.values())
assert len(fpset) == len(comps), (fp, top)

ext = [s for s in bnd if s not in fpset]
gmsh.model.addPhysicalGroup(3, [board_id], 1)
gmsh.model.addPhysicalGroup(2, sorted(ext), 1)
for i in range(len(comps)):
    gmsh.model.addPhysicalGroup(2, [fp[i]], 100 + i)

# size field
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", sorted(fpset))
gmsh.model.mesh.field.setNumber(1, "Sampling", 50)
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_min)
gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_max)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.3e-3)
gmsh.model.mesh.field.setNumber(2, "DistMax", 4.0e-3)
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
meta = {"lc_max_mm": lc_max*1e3, "lc_min_mm": lc_min*1e3, "nodes": nn, "tets": nt,
        "board_id": board_id, "flux_boundary_tags": {comps[i][0]: 100+i for i in range(len(comps))},
        "areas_m2": {comps[i][0]: comps[i][1]*comps[i][2]*MM*MM for i in range(len(comps))}}
json.dump(meta, open(out + "_meta.json", "w"), indent=1)
print("MESH", meta)
gmsh.finalize()
