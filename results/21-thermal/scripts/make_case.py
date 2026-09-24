#!/usr/bin/env python3
"""Heat sources + case.sif for one (load, h, Tamb). Board-only model: component
power applied as surface heat flux [W/m^2] on imprinted top-face footprint patches.
Radiation: idealized grey-body to ambient (eps, E) -- on/off via arg.
Usage: make_case.py <Iout> <h> <Tamb_C> <casedir> [eps|nord]
"""
import sys, os, json

Vf, DCR_L2, DCR_L1 = 0.70, 0.085, 0.027
RdsHS, RdsLS = 0.076, 0.032
Vin, Vout = 24.0, 12.0
D = Vout / Vin
Iout_full, Iin_full, Ptot_full = 3.0, 1.70, 5.0
areas = {"D1": 4.3*2.6*1e-6, "U21": 2.9*1.6*1e-6, "L2": 7.0*6.6*1e-6, "L1": 4.4*4.2*1e-6}
bcTag = {"D1": 100, "U21": 101, "L2": 102, "L1": 103}
order = ["D1", "U21", "L2", "L1"]

def losses(Iout):
    Iin = Iin_full * (Iout / Iout_full)
    known_full = (Vf*Iin_full + Iout_full**2*DCR_L2 + Iin_full**2*DCR_L1
                  + Iout_full**2*(D*RdsHS + (1-D)*RdsLS))
    sw_full = Ptot_full - known_full
    p = {"D1": Vf*Iin, "L2": Iout**2*DCR_L2, "L1": Iin**2*DCR_L1,
         "U21_cond": Iout**2*(D*RdsHS + (1-D)*RdsLS)}
    p["U21_sw_E"] = sw_full * (Iout/Iout_full)
    return p, Iin

def write_case(Iout, h, Tamb_C, casedir, k_xy, k_z, eps=0.9):
    p, Iin = losses(Iout)
    P = {"D1": p["D1"], "U21": p["U21_cond"] + p["U21_sw_E"], "L2": p["L2"], "L1": p["L1"]}
    flux = {k: P[k]/areas[k] for k in order}
    Tamb = Tamb_C + 273.15
    s = []
    s.append("Header\n  CHECK KEYWORDS Warn\n  Mesh DB \".\" \"mesh\"\n  Include Path \"\"\n  Results Directory \"\"\nEnd\n")
    s.append("Simulation\n  Max Output Level = 3\n  Coordinate System = Cartesian\n"
             "  Simulation Type = Steady state\n  Steady State Max Iterations = 40\nEnd\n")
    s.append("Constants\n  Stefan Boltzmann = 5.67e-8\nEnd\n")
    s.append(f"Material 1\n  Name = \"BoardAniso\"\n  Heat Conductivity(3) = Real {k_xy:.5f} {k_xy:.5f} {k_z:.5f}\nEnd\n")
    s.append("Equation 1\n  Name = \"Heat\"\n  Active Solvers(1) = 1\nEnd\n")
    s.append("Solver 1\n  Equation = Heat Equation\n  Procedure = \"HeatSolve\" \"HeatSolver\"\n"
             "  Stabilize = True\n  Linear System Solver = Iterative\n"
             "  Linear System Iterative Method = BiCGStab\n  Linear System Max Iterations = 800\n"
             "  Linear System Convergence Tolerance = 1e-9\n"
             "  Steady State Convergence Tolerance = 1e-7\nEnd\n")
    s.append("Solver 2\n  Equation = Result Output\n  Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n"
             "  Output File Name = \"case\"\n  Output Format = Vtu\n  Exec Solver = Always\nEnd\n")
    s.append("Solver 3\n  Equation = SaveScalars\n  Procedure = \"SaveData\" \"SaveScalars\"\n"
             "  Filename = \"scalars.dat\"\n  Variable 1 = Temperature\n  Operator 1 = max\n"
             "  Variable 2 = Temperature\n  Operator 2 = min\n  Exec Solver = Always\nEnd\n")
    s.append("Body 1\n  Target Bodies(1) = 1\n  Name = \"Board\"\n  Equation = 1\n  Material = 1\nEnd\n")
    bc = (f"Boundary Condition 1\n  Target Boundaries(1) = 1\n  Name = \"Convection\"\n"
          f"  Heat Transfer Coefficient = Real {h}\n  External Temperature = Real {Tamb:.2f}\n")
    if eps is not None:
        bc += (f"  Radiation = Idealized\n  Emissivity = Real {eps}\n"
               f"  Radiation External Temperature = Real {Tamb:.2f}\n")
    bc += "End\n"
    s.append(bc)
    for i, n in enumerate(order):
        s.append(f"Boundary Condition {2+i}\n  Target Boundaries(1) = {bcTag[n]}\n"
                 f"  Name = \"flux_{n}\"\n  Heat Flux = Real {flux[n]:.4f}\nEnd\n")
    os.makedirs(casedir, exist_ok=True)
    open(os.path.join(casedir, "case.sif"), "w").write("\n".join(s))
    hs = {"load": {"Iout_A": Iout, "Iin_A": round(Iin, 4), "D": D},
          "powers_W": {k: round(P[k], 5) for k in order},
          "U21_split_W": {"cond_D": round(p["U21_cond"], 5), "sw_other_E": round(p["U21_sw_E"], 5)},
          "total_W": round(sum(P.values()), 5),
          "surface_flux_W_per_m2": {k: round(flux[k], 2) for k in order},
          "footprint_area_m2": areas,
          "sources": {
              "D1": {"model": "P=Vf*Iin", "Vf": "0.70V D (SS36 HJC ds)", "Iin": "S (board-facts §5)"},
              "U21_cond": {"model": "P=Iout^2*(D*RdsHS+(1-D)*RdsLS)", "Rds": "76/32mOhm D (TI SLUSEA4D)"},
              "U21_sw_other": {"model": "loss-budget remainder*(Iout/3A)", "level": "E",
                               "replaces": "DS switching loss (no Qg table in DS)", "impact": "see REPORT"},
              "L2": {"model": "P=Iout^2*DCR", "DCR": "85mOhm D", "note": "I_rms~Iout (ripple ignored)"},
              "L1": {"model": "P=Iin^2*DCR", "DCR": "27mOhm D"}},
          "k_xy_W_per_mK": k_xy, "k_z_W_per_mK": k_z, "h_W_per_m2K": h, "Tamb_C": Tamb_C,
          "emissivity": eps}
    json.dump(hs, open(os.path.join(casedir, "heat_sources.json"), "w"), indent=1)
    return hs

if __name__ == "__main__":
    k_xy = (2*400*0.030e-3*0.812 + 0.3*1.43e-3)/1.49e-3
    k_z = 1.49e-3/(2*0.030e-3/400.0 + 1.43e-3/0.3)
    rad = sys.argv[5] if len(sys.argv) > 5 else "0.9"
    eps = None if rad == "nord" else float(rad)
    info = write_case(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], k_xy, k_z, eps)
    print(json.dumps(info["powers_W"]), "total=", info["total_W"], "eps=", eps)
