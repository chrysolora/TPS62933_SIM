#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_sif4.py <meshdir> <epsr_slab> <tb>  -- Elmer sif; eps_r(z) via MATC.
usage: make_sif4.py mesh4_fullcu_air 4.3 0.00143
"""
import os,re,sys
EPS0=8.854187817e-12
SIF="""Header
  Mesh DB "." "%(mesh)s"
End
Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 3D
  Simulation Type = Steady state
  Steady State Max Iterations = 1
  Output File = "%(mesh)s.result"
  Post File = "%(mesh)s.vtu"
End
Constants
  Permittivity of Vacuum = %(eps0).9e
End
Body 1
  Equation = 1
  Material = 1
End
Material 1
  Relative Permittivity = Variable Coordinate
    Real MATC "1.0 + (%(epsr).6f-1.0)*(tx(2)<=0.0)*(tx(2)>=%(zlo).9g)"
End
Equation 1
  Active Solvers(1) = 1
End
Solver 2
  Exec Solver = After Simulation
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = "scalars.%(mesh)s.dat"
  Variable 1 = "electric energy density"
  Operator 1 = "volume integral"
  Variable 2 = "potential"
  Operator 2 = "volume"
End
Solver 1
  Equation = "StatElec"
  Procedure = "StatElecSolve" "StatElecSolver"
  Variable = Potential
  Calculate Electric Energy = Logical True
  Linear System Solver = "Direct"
  Linear System Direct Method = "MUMPS"
  Linear System Symmetric = True
  Linear System Row Equilibration = Logical True
End
Boundary Condition 1
  Target Boundaries(1) = %(sw)s
  Potential = 1.0
End
Boundary Condition 2
  Target Boundaries(1) = %(gnd)s
  Potential = 0.0
End
Boundary Condition 3
  Target Boundaries(1) = %(ref)s
  Potential = 0.0
End
%(topbnd)s"""
def names(md):
    d={}
    for m in re.finditer(r"\$\s*(\w+)\s*=\s*(\d+)", open(os.path.join(md,"mesh.names")).read()):
        d[m.group(1)]=int(m.group(2))
    return d
if __name__=="__main__":
    md=sys.argv[1].rstrip("/"); epsr=float(sys.argv[2]); tb=float(sys.argv[3])
    nm=names(md); print("names:",nm)
    top=""
    if "topgnd" in nm:
        top="\nBoundary Condition 4\n  Target Boundaries(1) = %d\n  Potential = 0.0\nEnd\n"%nm["topgnd"]
    s=SIF%dict(mesh=os.path.basename(md),eps0=EPS0,epsr=epsr,zlo=-tb,sw=nm["sw"],gnd=nm["gnd"],ref=nm["ref"],topbnd=top)
    out="case4_%s.sif"%os.path.basename(md); open(out,"w").write(s); print("wrote",out)
