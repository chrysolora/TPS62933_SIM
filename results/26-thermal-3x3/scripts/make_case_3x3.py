#!/usr/bin/env python3
"""Build case.sif for TPS62933 board thermal 3x3 study.
Board materials: BoardMain (background k_xy) + BoardCut (keep-out region, reduced k_xy).
4 component BODIES with VOLUMETRIC heat source (Elmer: 'Volumetric Heat Source').
Loss budget identical to thermal_final/make_case_final.py (no residual term).
Usage: make_case_3x3.py <Iout> <h> <Tamb_C> <casedir> <outer> <term_in> <term_out> <bodies.json>
Env: CUT_MODE (fullcu|topcut|dualcut), KZ (W/mK)
"""
import sys, os, json

# ---- level D parameters (sourced, same as thermal_final) ----
Vf      = 0.70
DCR_L2  = 0.085
DCR_L1  = 0.027
RdsHS   = 0.076
RdsLS   = 0.032
Vin, Vout = 24.0, 12.0
D       = Vout / Vin
Iout_full, Iin_full = 3.0, 1.70
PSW_FULL = 0.50   # E
GEOM = {"U21": (2.10, 1.60, 0.60), "D1": (4.30, 2.60, 2.20),
        "L2": (6.50, 6.50, 3.00), "L1": (4.20, 4.20, 2.00)}
_RTH0 = {"U21": (19.3, "D"), "D1": (30.0, "E"), "L2": (15.0, "E"), "L1": (20.0, "E")}
RTH = {k: v for k, v in _RTH0.items()}
order = ["U21", "D1", "L2", "L1"]
KZ = float(os.environ.get("KZ", "1.3612"))
VF = float(os.environ.get("VF", "0.70"))
PSW_ENV = float(os.environ.get("PSW", str(PSW_FULL)))
PSW_FULL = PSW_ENV
CUT_MODE = os.environ.get("CUT_MODE", "fullcu")

# ---- in-plane conductivity (30um Cu + FR4), coverage 0.812 (thermal_final convention) ----
K_CU, K_FR4 = 400.0, 0.30
TCU, T_FR4, T_TOT = 0.030e-3, 1.43e-3, 1.49e-3
COV = 0.812
KXY_BG = (2 * K_CU * TCU * COV + K_FR4 * T_FR4) / T_TOT
LAYER = K_CU * TCU * COV / T_TOT           # one copper layer contribution
if CUT_MODE == "fullcu":
    KXY_CUT = KXY_BG
elif CUT_MODE == "topcut":
    KXY_CUT = KXY_BG - LAYER
elif CUT_MODE == "dualcut":
    KXY_CUT = KXY_BG - 2 * LAYER
else:
    raise SystemExit("bad CUT_MODE " + CUT_MODE)

G_FIN_PER_TERM = 2.328e-3
TERM_AREA = 3.81e-3 * 3.81e-3


def losses(Iout):
    f = Iout / Iout_full
    Iin = Iin_full * f
    p = {"D1": VF*Iin, "L2": Iout**2*DCR_L2, "L1": Iin**2*DCR_L1,
         "U21_cond": Iout**2*(D*RdsHS + (1-D)*RdsLS), "U21_sw": PSW_FULL * f}
    return p, Iin


def kcomp(name):
    w, l, h = GEOM[name]
    A = w*1e-3 * l*1e-3
    return h*1e-3 / (RTH[name][0] * A)


def write_case(Iout, h, Tamb_C, casedir, outer_bc, term_in, term_out, bodies, fin_on=True, eps=0.9):
    p, Iin = losses(Iout)
    P = {"U21": p["U21_cond"] + p["U21_sw"], "D1": p["D1"], "L2": p["L2"], "L1": p["L1"]}
    Tamb = Tamb_C + 273.15
    s = []
    s.append("Header\n  CHECK KEYWORDS Warn\n  Mesh DB \".\" \"mesh\"\n  Include Path \"\"\n  Results Directory \"\"\nEnd\n")
    s.append("Simulation\n  Max Output Level = 3\n  Coordinate System = Cartesian\n"
             "  Simulation Type = Steady state\n  Steady State Max Iterations = 40\nEnd\n")
    s.append("Constants\n  Stefan Boltzmann = 5.67e-8\nEnd\n")
    s.append(f"Material 1\n  Name = \"BoardMain\"\n  Heat Conductivity(3) = Real {KXY_BG:.5f} {KXY_BG:.5f} {KZ:.5f}\nEnd\n")
    s.append(f"Material 2\n  Name = \"BoardCut\"\n  Heat Conductivity(3) = Real {KXY_CUT:.5f} {KXY_CUT:.5f} {KZ:.5f}\nEnd\n")
    for i, n in enumerate(order):
        s.append(f"Material {3+i}\n  Name = \"pkg_{n}\"\n  Heat Conductivity = Real {kcomp(n):.5f}\nEnd\n")
    s.append("Body Force 1\n  Name = \"qsrc_U21\"\n  Volumetric Heat Source = Real %.4e\nEnd\n"
             % (P["U21"] / (GEOM["U21"][0]*GEOM["U21"][1]*GEOM["U21"][2]*1e-9)))
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
    bidx = 1
    for bt in bodies["board_main"]:
        s.append("Body %d\n  Target Bodies(1) = %d\n  Name = \"BoardMain\"\n  Equation = 1\n  Material = 1\nEnd\n" % (bidx, bt))
        bidx += 1
    for bt in bodies["cutout"]:
        s.append("Body %d\n  Target Bodies(1) = %d\n  Name = \"BoardCut\"\n  Equation = 1\n  Material = 2\nEnd\n" % (bidx, bt))
        bidx += 1
    for i, n in enumerate(order):
        s.append("Body %d\n  Target Bodies(1) = %d\n  Name = \"%s\"\n  Equation = 1\n  Material = %d\n  Body Force = %d\nEnd\n"
                 % (bidx, bodies["comps"][n], n, 3+i, 1+i))
        bidx += 1
    bc = (f"Boundary Condition 1\n  Target Boundaries(1) = {outer_bc}\n  Name = \"Convection\"\n"
          f"  Heat Transfer Coefficient = Real {h}\n  External Temperature = Real {Tamb:.2f}\n")
    if eps is not None:
        bc += (f"  Radiation = Idealized\n  Emissivity = Real {eps}\n"
               f"  Radiation External Temperature = Real {Tamb:.2f}\n")
    bc += "End\n"
    s.append(bc)
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
    hs = {"cut_mode": CUT_MODE,
          "load": {"Iout_A": Iout, "Iin_A": round(Iin, 4), "D": D},
          "powers_W": {k: round(P[k], 5) for k in order},
          "U21_split_W": {"cond_D": round(p["U21_cond"], 5), "sw_E": round(p["U21_sw"], 5)},
          "total_W": round(sum(P.values()), 5),
          "efficiency_est": round(Vout*Iout/(Vout*Iout+sum(P.values())), 4),
          "comp_effective_k_W_per_mK": {n: round(kcomp(n), 4) for n in order},
          "comp_Rtheta_C_per_W": {n: RTH[n][0] for n in order},
          "comp_Rtheta_level": {n: RTH[n][1] for n in order},
          "board_k_xy_bg_W_mK": round(KXY_BG, 4),
          "board_k_xy_cut_W_mK": round(KXY_CUT, 4),
          "board_k_z_W_mK": KZ,
          "copper_coverage_layers": COV, "cu_layer_contribution_W_mK": round(LAYER, 4),
          "h_W_per_m2K": h, "Tamb_C": Tamb_C, "emissivity": eps,
          "fin_on": fin_on, "G_fin_per_terminal_W_per_K": G_FIN_PER_TERM,
          "sources": {
              "D1": {"model": "P=Vf*Iin", "Vf": "0.70V (SS36 ds, D)", "level": "D"},
              "U21_cond": {"model": "P=Iout^2*(D*RdsHS+(1-D)*RdsLS)", "Rds": "76/32 mOhm (TI SLUSEA4D, D)",
                           "level": "D"},
              "U21_sw": {"model": "P_sw=1/2*Vin*Iout*fsw*(tr+tf); fsw=805kHz; tr=tf=5ns -> 0.29W; "
                                  "nominal 0.50W incl Coss/gate", "level": "E", "range_W_fullload": [0.29, 0.70]},
              "L2": {"model": "P=Iout^2*DCR", "DCR": "85 mOhm (LCSC C41415624, D)", "level": "D"},
              "L1": {"model": "P=Iin^2*DCR", "DCR": "27 mOhm (LCSC C167203, D)", "level": "D"}},
          "Rtheta_source": {"U21": "RthetaJB=19.3 C/W (TI SLUSEA4D, D)",
                            "D1": "~30 C/W (E)", "L2": "~15 C/W (E)", "L1": "~20 C/W (E)"}}
    json.dump(hs, open(os.path.join(casedir, "heat_sources.json"), "w"), indent=1)
    return hs


if __name__ == "__main__":
    a = sys.argv
    bodies = json.load(open(a[8]))
    info = write_case(float(a[1]), float(a[2]), float(a[3]), a[4],
                      int(a[5]), int(a[6]), int(a[7]), bodies)
    print(json.dumps(info["powers_W"]), "total=", info["total_W"], "eta=", info["efficiency_est"],
          "cut_kxy=", info["board_k_xy_cut_W_mK"])
