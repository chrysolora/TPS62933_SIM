import sys, numpy as np, collections, os, shutil

meshdir = sys.argv[1]; out = sys.argv[2]

# ---------- parse nodes ----------
xyz = {}
with open(meshdir + "/mesh.nodes") as f:
    for l in f:
        s = l.split()
        if len(s) == 5 and s[0].isdigit():
            xyz[int(s[0])] = (float(s[2]), float(s[3]), float(s[4]))
NT = {101: 2, 202: 3, 303: 3, 404: 4, 504: 4, 808: 8, 706: 6, 306: 3, 604: 6}

def read_blocks(path):
    """yield (first_line_tokens, [cont_lines]) with index into file."""
    lines = open(path).read().split("\n")
    out_lines = [l.split() for l in lines]
    i = 0
    blocks = []
    while i < len(out_lines):
        s = out_lines[i]
        if len(s) >= 3 and s[0].isdigit() and s[2].isdigit():
            et = int(s[2]); k = NT.get(et, 0); toks = s[3:]; j = i + 1
            while len(toks) < k and j < len(out_lines):
                toks += out_lines[j]; j += 1
            blocks.append((i, s, [out_lines[q] for q in range(i + 1, j)]))
            i = j
        else:
            i += 1
    return lines, out_lines, blocks

# ---------- bodies ----------
lines, ol, blocks = read_blocks(meshdir + "/mesh.elements")
bodynodes = collections.defaultdict(set)
for i, s, cont in blocks:
    b = int(s[1]); k = NT.get(int(s[2]), 0)
    toks = s[3:]
    for c in cont: toks += c
    for t in toks[:k]: bodynodes[b].add(int(t))

def classify(v):
    if v["z1"] < -0.5e-3: return "cobot"
    if v["z1"] <= 1e-6: return "cotop"
    if v["z0"] >= 0.8e-3 and v["z1"] <= 1.6e-3: return "coil"
    if v["z0"] >= 0.1e-3 and max(abs(v["x0"]), abs(v["x1"])) <= 3.3e-3: return "core"
    return "air"

role_old = {}
for b, ns in bodynodes.items():
    a = np.array([xyz[k] for k in ns])
    v = dict(n=len(ns), x0=a[:, 0].min(), x1=a[:, 0].max(), y0=a[:, 1].min(), y1=a[:, 1].max(),
             z0=a[:, 2].min(), z1=a[:, 2].max())
    r = classify(v); role_old[r] = b
    print("body %d -> %-6s n=%5d z[%.5f,%.5f]" % (b, r, v["n"], v["z0"], v["z1"]))

ORDER = ["cobot", "cotop", "coil", "air", "core"]
new_of_role = {r: i + 1 for i, r in enumerate(ORDER)}
bmap = {role_old[r]: new_of_role[r] for r in ORDER}

# rewrite elements body field
newlines = list(lines)
for i, s, cont in blocks:
    s2 = list(s); s2[1] = str(bmap[int(s[1])])
    newlines[i] = " ".join(s2)
open(meshdir + "/mesh.elements.tmp", "w").write("\n".join(newlines))
shutil.move(meshdir + "/mesh.elements.tmp", meshdir + "/mesh.elements")

# ---------- boundaries ----------
exists = os.path.exists(meshdir + "/mesh.boundary")
if exists:
    blines, bol, bblocks = read_blocks(meshdir + "/mesh.boundary")
    bnodes = collections.defaultdict(set)
    for i, s, cont in bblocks:
        bc = int(s[1]); et = int(s[4]); k = NT.get(et, 3)
        for t in s[-k:]: bnodes[bc].add(int(t))
    gn = np.array([xyz[k] for k in xyz])
    gx0, gx1 = gn[:, 0].min(), gn[:, 0].max()
    gy0, gy1 = gn[:, 1].min(), gn[:, 1].max()
    gz0, gz1 = gn[:, 2].min(), gn[:, 2].max()
    tol = 1e-9
    outer = []
    for bc, ns in bnodes.items():
        ok = True
        for k in ns:
            x, y, z = xyz[k]
            if not (abs(x - gx0) < tol or abs(x - gx1) < tol or abs(y - gy0) < tol or
                    abs(y - gy1) < tol or abs(z - gz0) < tol or abs(z - gz1) < tol):
                ok = False; break
        if ok: outer.append(bc)
    allbc = sorted(bnodes)
    brank = {b: i + 1 for i, b in enumerate(allbc)}
    newouter = sorted(brank[b] for b in outer)
    nl = list(blines)
    for i, s, cont in bblocks:
        s2 = list(s); s2[1] = str(brank[int(s[1])])
        nl[i] = " ".join(s2)
    shutil.move(meshdir + "/mesh.boundary", meshdir + "/mesh.boundary.orig")
    open(meshdir + "/mesh.boundary", "w").write("\n".join(nl))
    print("outer(orig)=%s -> new=%s of %d boundaries" % (sorted(outer), newouter, len(allbc)))
else:
    newouter = []

coil_id = new_of_role["coil"]

sif = """Header
  Mesh DB "." "mesh"
End
Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
End
Constants
  Permeability of Vacuum = 1.25664e-6
End
Body 1
  Target Bodies(1) = {cotop}
  Name = "CopperTop"
  Equation = 2
  Material = 1
End
Body 2
  Target Bodies(1) = {cobot}
  Name = "CopperBot"
  Equation = 2
  Material = 1
End
Body 3
  Target Bodies(1) = {coil}
  Name = "Coil"
  Equation = 1
  Material = 4
  Body Force = 1
End
Body 4
  Target Bodies(1) = {air}
  Name = "Air"
  Equation = 2
  Material = 2
End
Body 5
  Target Bodies(1) = {core}
  Name = "Core"
  Equation = 2
  Material = 3
End
Component 1
  Master Bodies(1) = {coil}
  Coil Closed = Logical True
  Coil Normal(3) = 0 0 1
  Desired Coil Current = 6.0
  Coil Type = String "Stranded"
End
Equation 1
  Name = "CoilEq"
  Active Solvers(3) = 1 2 3
End
Equation 2
  Name = "MagEq"
  Active Solvers(2) = 2 3
End
Solver 1
  Equation = "Coil"
  Procedure = "CoilSolver" "CoilSolver"
  Calculate Elemental Fields = Logical True
  Linear System Solver = Direct
  Linear System Direct Method = UMFPACK
End
Solver 2
  Equation = "MGDynamics"
  Procedure = "MagnetoDynamics" "WhitneyAVHarmonicSolver"
  Use Elemental CoilCurrent = Logical True
  Angular Frequency = 4925663.2
  Linear System Solver = Direct
  Linear System Direct Method = MUMPS
End
Solver 3
  Equation = "MGDynCalc"
  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
  Calculate Current Density = True
  Calculate Joule Heating = True
  Calculate Magnetic Field Strength = True
  Calculate Magnetic Flux Density = True
  Linear System Solver = Direct
  Linear System Direct Method = UMFPACK
End
Solver 4
  Exec Solver = After All
  Equation = "ResultOutput"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = case
  Output Format = vtu
  Ascii Output = Logical True
End
Material 1
  Name = "Copper"
  Relative Permeability = 1.0
  Electric Conductivity = 5.8e7
End
Material 2
  Name = "Air"
  Relative Permeability = 1.0
  Electric Conductivity = 0.0
End
Material 3
  Name = "Core"
  Relative Permeability = 40.0
  Electric Conductivity = 0.0
End
Material 4
  Name = "Coil"
  Relative Permeability = 1.0
  Electric Conductivity = 1.0
End
Body Force 1
  Current Density = 0.0
End
Boundary Condition 1
  Target Boundaries({nbc}) = {bc}
  AV = Real 0.0
End
""".format(cotop=new_of_role["cotop"], cobot=new_of_role["cobot"], coil=coil_id,
           air=new_of_role["air"], core=new_of_role["core"], nbc=len(newouter),
           bc=" ".join(map(str, newouter)))
open(out, "w").write(sif)
print("wrote", out, "body map:", bmap)
