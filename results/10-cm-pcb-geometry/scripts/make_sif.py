#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write Elmer case.sif for the electrostatic SW<->ref-plane capacitance."""
import os, sys

SIF = """Header
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
  Permittivity of Vacuum = 8.854187817e-12
End

Body 1
  Equation = 1
  Material = 1
End

Material 1
  Relative Permittivity = 1.0
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
"""

if __name__ == "__main__":
    mesh = sys.argv[1]
    sw, gnd, ref = sys.argv[2], sys.argv[3], sys.argv[4]
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "case_%s.sif" % mesh), "w").write(
        SIF % dict(mesh=mesh, sw=sw, gnd=gnd, ref=ref))
    print("case_%s.sif" % mesh)
