#!/usr/bin/env python3
"""Build case.sif v2 for TPS62933 board thermal.
v2 model: board(aniso k) + 4 component BODIES with VOLUMETRIC heat source.
Each component body k chosen so its die->board thermal resistance = Rtheta.
Convection h(DeltaT) imposed by outer Python iteration (driver run_v2.py).
Wire/terminal fin path added as equivalent extra conductance at terminals.
Usage: make_case_v2.py <Iout> <h> <Tamb_C> <casedir> <outer_bc> <term_bc_in> <term_bc_out> <fin_on|fin_off> [eps|nord]
"""
import sys, os, json

Vf, DCR_L2, DCR_L1 = 0.70, 0.085, 0.027
RdsHS, RdsLS = 0.076, 0.032
Vin, Vout = 24.0, 12.0
D = Vout / Vin
Iout_full, Iin_full, Ptot_full = 3.0, 1.70, 5.0

# component geometry (must match build_mesh_v2.py)
GEOM = {"U21": (2.10, 1.60, 0.60), "D1": (4.30, 2.60, 2.20),
        "L2": (6.50, 6.50, 3.00), "L1": (4.20, 4.20, 2.00)}
# die->board thermal resistance (C/W) and source level
TH_SCALE = float(os.environ.get("THETA_SCALE", "1.0"))
_RTH0 = {"U21": (19.3, "D"), "D1": (30.0, "E"), "L2": (15.0, "E"), "L1": (20.0, "E")}
RTH = {k: (v[0]*TH_SCALE, v[1]) for k, v in _RTH0.items()}
BODY_TAG = {"U21": int(os.environ.get("BODY_U21","1")), "D1": int(os.environ.get("BODY_D1","1")),
            "L2": int(os.environ.get("BODY_L2","1")), "L1": int(os.environ.get("BODY_L1","1"))}
BOARD_TAG = int(os.environ.get("BODY_BOARD","1"))
order = ["U21", "D1", "L2", "L1"]

# terminal fin: 1 wire (0.75mm^2, 10cm Cu) per terminal. Fin conductance 2.328e-3 W/K each (see report).
G_FIN_PER_TERM = 2.328e-3
TERM_AREA = 3.81e-3 * 3.81e-3


def losses(Iout):
    Iin = Iin_full * (Iout / Iout_full)
    p = {"D1": Vf*Iin, "L2": Iout**2*DCR_L2, "L1": Iin**2*DCR_L1,
         "U21_cond": Iout**2*(D*RdsHS + (1-D)*RdsLS)}
    known_full = (Vf*Iin_full + Iout_full**2*DCR_L2 + Iin_full**2*DCR_L1
                  + Iout_full**2*(D*RdsHS + (1-D)*RdsLS))
    p["U21_sw_E"] = (Ptot_full - known_full) * (Iout/Iout_full)
    return p, Iin


def kcomp(name):
    w, l, h = GEOM[name]
    A = w*1e-3 * l*1e-3
    R = RTH[name][0]
    return h*1e-3 / (R * A)   # W/mK


def write_case(Iout, h, Tamb_C, casedir, outer_bc, term_in, term_out,
               fin_on=True, eps=0.9):
    p, Iin = losses(Iout)
    P = {"U21": p["U21_cond"] + p["U21_sw_E"], "D1": p["D1"], "L2": p["L2"], "L1": p["L1"]}
    Tamb = Tamb_C + 273.15
    k_xy = (2*400*0.030e-3*0.812 + 0.3*1.43e-3)/1.49e-3
    k_z = 1.49e-3/(2*0.030e-3/400.0 + 1.43e-3/0.3)
    s = []
    s.append("Header\n  CHECK KEYWORDS Warn\n  Mesh DB \".\" \"mesh\"\n  Include Path \"\"\n  Results Directory \"\"\nEnd\n")
    s.append("Simulation\n  Max Output Level = 3\n  Coordinate System = Cartesian\n"
             "  Simulation Type = Steady state\n  Steady State Max Iterations = 40\nEnd\n")
    s.append("Constants\n  Stefan Boltzmann = 5.67e-8\nEnd\n")
    s.append(f"Material 1\n  Name = \"BoardAniso\"\n  Heat Conductivity(3) = Real {k_xy:.5f} {k_xy:.5f} {k_z:.5f}\nEnd\n")
    for i, n in enumerate(order):
        s.append(f"Material {2+i}\n  Name = \"pkg_{n}\"\n  Heat Conductivity = Real {kcomp(n):.5f}\nEnd\n")
    s.append("Body Force 1\n  Name = \"qsrc_U21\"\n  Volumetric Heat Source = Real %.4e\nEnd\n" % (P["U21"] / (GEOM["U21"][0]*GEOM["U21"][1]*GEOM["U21"][2]*1e-9)))
    for i, n in enumerate(["D1", "L2", "L1"]):
        V = GEOM[n][0]*GEOM[n][1]*GEOM[n][2]*1e-9
        s.append("Body Force %d\n  Name = \"qsrc_%s\"\n  Volumetric Heat Source = Real %.4e\nEnd\n" % (2+i, n, P[n]/V))
    s.append("Equation 1\n  Name = \"Heat\"\n  Active Solvers(1) = 1\nEnd\n")
    s.append("Solver 1\n  Equation = Heat Equation\n  Procedure = \"HeatSolve\" \"HeatSolver\"\n"
             "  Stabilize = True\n  Linear System Solver = Iterative\n"
             "  Linear System Iterative Method = BiCGStab\n  Linear System Max Iterations = 1000\n"
             "  Linear System Convergence Tolerance = 1e-9\n"
             "  Nonlinear System Max Iterations = 1\n"
             "  Steady State Convergence Tolerance = 1e-8\nEnd\n")
    s.append("Solver 2\n  Equation = Result Output\n  Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n"
             "  Output File Name = \"case\"\n  Output Format = Vtu\n  Exec Solver = Always\nEnd\n")
    s.append("Solver 3\n  Equation = SaveScalars\n  Procedure = \"SaveData\" \"SaveScalars\"\n"
             "  Filename = \"scalars.dat\"\n  Variable 1 = Temperature\n  Operator 1 = max\n"
             "  Variable 2 = Temperature\n  Operator 2 = min\n  Exec Solver = Always\nEnd\n")
    s.append("Body 1\n  Target Bodies(1) = %d\n  Name = \"Board\"\n  Equation = 1\n  Material = 1\nEnd\n" % BOARD_TAG)
    for i, n in enumerate(order):
        s.append("Body %d\n  Target Bodies(1) = %d\n  Name = \"%s\"\n  Equation = 1\n  Material = %d\n  Body Force = %d\nEnd\n"
                 % (2+i, BODY_TAG[n], n, 2+i, 1+(["U21","D1","L2","L1"].index(n))))
    # outer convection + radiation
    bc = (f"Boundary Condition 1\n  Target Boundaries(1) = {outer_bc}\n  Name = \"Convection\"\n"
          f"  Heat Transfer Coefficient = Real {h}\n  External Temperature = Real {Tamb:.2f}\n")
    if eps is not None:
        bc += (f"  Radiation = Idealized\n  Emissivity = Real {eps}\n"
               f"  Radiation External Temperature = Real {Tamb:.2f}\n")
    bc += "End\n"
    s.append(bc)
    # terminals: convection + optional fin
    for j, (tb, nm) in enumerate([(term_in, "TERM_IN"), (term_out, "TERM_OUT")]):
        heff = h + (G_FIN_PER_TERM / TERM_AREA if fin_on else 0.0)
        tb_bc = (f"Boundary Condition {2+j}\n  Target Boundaries(1) = {tb}\n  Name = \"{nm}\"\n"
                 f"  Heat Transfer Coefficient = Real {heff:.4f}\n  External Temperature = Real {Tamb:.2f}\n")
        if eps is not None:
            tb_bc += (f"  Radiation = Idealized\n  Emissivity = Real {eps}\n"
                      f"  Radiation External Temperature = Real {Tamb:.2f}\n")
        tb_bc += "End\n"
        s.append(tb_bc)
    os.makedirs(casedir, exist_ok=True)
    open(os.path.join(casedir, "case.sif"), "w").write("\n".join(s))
    hs = {"load": {"Iout_A": Iout, "Iin_A": round(Iin, 4), "D": D},
          "powers_W": {k: round(P[k], 5) for k in order},
          "U21_split_W": {"cond_D": round(p["U21_cond"], 5), "sw_other_E": round(p["U21_sw_E"], 5)},
          "total_W": round(sum(P.values()), 5),
          "comp_effective_k_W_per_mK": {n: round(kcomp(n), 4) for n in order},
          "comp_Rtheta_C_per_W": {n: RTH[n][0] for n in order},
          "comp_Rtheta_level": {n: RTH[n][1] for n in order},
          "h_W_per_m2K": h, "Tamb_C": Tamb_C, "emissivity": eps,
          "fin_on": fin_on, "G_fin_per_terminal_W_per_K": G_FIN_PER_TERM,
          "sources": {
              "D1": {"model": "P=Vf*Iin", "Vf": "0.70V (SS36 HJC ds, D)", "Iin": "S"},
              "U21_cond": {"model": "P=Iout^2*(D*RdsHS+(1-D)*RdsLS)", "Rds": "76/32 mOhm (TI SLUSEA4D, D)"},
              "U21_sw_other": {"model": "loss-budget remainder", "level": "E",
                               "replaces": "DS switching loss", "impact": "see REPORT"},
              "L2": {"model": "P=Iout^2*DCR", "DCR": "85 mOhm (D)"},
              "L1": {"model": "P=Iin^2*DCR", "DCR": "27 mOhm (D)"}},
          "Rtheta_source": {"U21": "RthetaJB=19.3 C/W (TI SLUSEA4D Thermal Information, D)",
                            "D1": "~30 C/W (E: SMA diode, no DS theta; sensitivity +-50%)",
                            "L2": "~15 C/W (E: molded inductor, no DS theta; sensitivity +-50%)",
                            "L1": "~20 C/W (E: molded inductor, no DS theta; sensitivity +-50%)"}}
    json.dump(hs, open(os.path.join(casedir, "heat_sources.json"), "w"), indent=1)
    return hs


if __name__ == "__main__":
    a = sys.argv
    fin = (a[8] == "fin_on")
    eps = None if (len(a) > 9 and a[9] == "nord") else 0.9
    info = write_case(float(a[1]), float(a[2]), float(a[3]), a[4],
                      int(a[5]), int(a[6]), int(a[7]), fin, eps)
    print(json.dumps(info["powers_W"]), "total=", info["total_W"], "fin=", fin, "h=", a[2])
