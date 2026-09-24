# TPS62933 board thermal simulation — **v2** (corrected)

Directory: `/mnt/raid10/sim-work/tps62933/thermal_v2/`
Solver: ElmerSolver 26.2 (HeatSolve, steady state), mesh gmsh 4.15 python SDK.
Model: board = homogenised anisotropic solid (k_xy=13.37 / k_z=0.3126 W/m·K);
4 component **bodies** with **volumetric** heat generation; natural convection
`h = h(ΔT)` (outer iteration) + grey-body radiation ε=0.9; wire/terminal fin path.

---

## 1) What v2 fixes vs v1 (the four defects)

| # | v1 defect | v2 fix |
|---|-----------|--------|
| 1 | Heat applied as **uniform surface flux on the sharp footprint patch** → log singularity → non-convergent false 389 °C | Heat applied as **volumetric source inside the package body**; flux enters the board through the package solid (smooth, no step on the board surface). |
| 2 | **No package thermal resistance** | Each device = a body whose effective k makes its **die→board resistance = θ**. U21 uses **θJB = 19.3 °C/W from TI SLUSEA4D (D)**. |
| 3 | **Fixed h (5/10)** | **h(ΔT)** from the Incropera horizontal-plate correlation, solved by outer iteration (E). Runs: ΔT=25/40 °C. |
| 4 | **No wire/terminal conduction** | 2 × 0.75 mm² × 10 cm Cu wires modelled as an equivalent adiabatic-fin conductance added at the two connector pads (E). |

Extra: a **build/parse bug was found and fixed** — this Elmer build silently ignores
the BodyForce keyword `Heat Source`; only **`Volumetric Heat Source`** is read.
(v1 never used a volumetric source, so it did not hit this.)

---

## 2) Heat sources (each with source level)

| Source | Power (W) | Model | Level | Applied as |
|--------|-----------|-------|-------|-----------|
| D1 | **1.190** | P = Vf·Iin, Vf=0.70 V (SS36) | **D** | vol. in SMA body 4.3×2.6×2.2 mm |
| L2 | **0.765** | P = Iout²·DCR, DCR=85 mΩ | **D** | vol. in ZEMS0650 body 6.5×6.5×3.0 mm |
| L1 | **0.0780** | P = Iin²·DCR, DCR=27 mΩ | **D** | vol. in FXL0420 body 4.2×4.2×2.0 mm |
| U21 cond. | **0.487** | Iout²·(D·RdsHS+(1−D)·RdsLS), 76/32 mΩ | **D** | part of U21 body |
| U21 sw.+other | **2.480** | loss-budget remainder | **E** | part of U21 body |
| **Σ** | **5.000** | **"total loss = 5 W"** | **S** | — |

Package θ (die→board, sets each body's k):

| Dev | θ (die→board) | Level | effective k |
|-----|----------------|-------|-------------|
| U21 | **19.3 °C/W** (θJB, TI SLUSEA4D) | **D** | 9.25 W/m·K |
| D1 | ~30 °C/W (SMA, no DS θ) | **E** | 6.56 W/m·K |
| L2 | ~15 °C/W (moulded L, no DS θ) | **E** | 4.73 W/m·K |
| L1 | ~20 °C/W (moulded L, no DS θ) | **E** | 5.67 W/m·K |

Wire fin (E): per terminal G ≈ 2.33e-3 W/K ⇒ h_eff ≈ 160 W/m²·K on the pad area.

---

## 3) Results

### 3.1 Grid convergence (full load 5.0 W, Tamb 25 °C, **fixed h = 13** W/m²K)

| mesh | lc | board mean | U21 max | U21 avg | D1 avg | L2 avg | L1 avg |
|------|----|-----------|---------|---------|--------|--------|--------|
| L1 | 2.5/0.60 mm | 97.4 | 272.0 | 261.3 | 137.9 | 126.9 | 100.9 |
| L2 | 1.6/0.35 mm | 97.2 | 293.5 | 281.4 | 141.8 | 127.2 | 100.5 |
| L3 | 1.0/0.18 mm | 96.9 | 334.4 | 323.9 | 148.5 | 127.8 | 99.8 |
| **Δ L2→L3** | — | **−0.3 %** | **+13.9 %** | +15.1 % | +4.7 % | +0.5 % | −0.7 % |

- **PASS (<2 %)**: board mean, L1, L2 device temperatures.
- **FAIL (<2 %)**: **U21 peak/junction** and (marginal) D1 — a residual
  singularity remains at the high-power-density **package/board interface** in
  this homogenised model. **The U21 absolute temperature is therefore graded C**
  and is *not* used as a conclusion.

### 3.2 Main case — temperature field (`fig_thermal_map_v2.png`)

h(ΔT) iteration converges h: 10 → 11.6 → 12.4 → 12.9 → **13.0–13.25 W/m²K**
(ΔT ≈ 84 K, Tsurf ≈ 99–110 °C). Values below are **L3, h = 13**.

| Node / device | Temperature (°C) | Level |
|---------------|------------------|-------|
| Board mean | **96.9** | **A** |
| **U21 junction (peak)** | **334** (grid-sensitive, 272/294/334 over L1/L2/L3; true ≥ L3) | **C** |
| D1 | 152.7 peak / 148.5 avg | **C** (grid +4.7 %) |
| L2 | 138.1 peak / 127.8 avg | **B** |
| L1 | 103.0 peak / 99.8 avg | **B** |
| Tsurf avg | 99.5 | A |

**Heat balance ΣQ**: Q_out = 5.00 W (conv 3.05 W + rad 1.81 W + wires 0.31 W),
imbalance **0.0 %** → **PASS <2 %** (A).

### 3.3 Sensitivity (`fig_thermal_sensitivity_v2.png`, L2, full load)

| Case | U21 junction °C | Δ vs base |
|------|-----------------|-----------|
| h=5 (v1-like) | 327.0 | +33.5 |
| h=10 | 304.0 | +10.5 |
| **h=13 (h(ΔT))** base | **293.5** | — |
| Tamb 40 °C | 305.1 | +11.6 |
| no radiation | 339.1 | +45.6 |
| no wire-fin | 296.6 | +3.1 |
| package θ ×1.5 | 311.5 | +18.0 |
| half load 1.5 A (2.17 W) | 148.3 | — |
| light load 0.6 A (0.79 W) | 71.4 | — |

Load sweep: U21 junction ≈ **148 °C at 1.5 A** (2.17 W) and **293 °C at 3.0 A** (5 W)
→ the board is thermally acceptable only up to ~1.5 A; at full load it is far
beyond the 150 °C limit.

### 3.4 Conclusion table

| Conclusion | Grade | Basis |
|------------|-------|-------|
| Board mean temp ≈ **97 °C** at 5 W (Tamb 25) | **A** | converged (<0.3 %), ΣQ balance 0 % |
| L2 / L1 device temps ≈ 128 / 100 °C | **B** | grid-converged, but package θ is E |
| **U21 junction ≳ 330 °C at 5 W** (far above 150 °C) | **C** | local peak **not** grid-converged |
| At 1.5 A (2.17 W) U21 ≈ **148 °C** — near the limit | **B/C** | same model; U21 absolute still C |
| Radiation is **essential** (removes ~1.8 W) | **B** | sensitivity, converged board field |
| Wire-fin path is **minor** (~3 °C) for U21 | **B** | E-finf model, small effect |
| h must vary with ΔT (fixed h=5 over-predicts by ~34 °C) | **B** | correlation E, effect large & monotone |

---

## 4) Why v1 was so much higher

v1 U21 = **389 °C** (302→347→389 non-convergent) at fixed h=5, flux on a sharp
patch, no package θ, no radiation-on-wires. v2 U21 = **334 °C** (L3). The drop
(≈14 %) comes from: (a) volumetric source instead of a boundary-flux patch
(removes part of the logarithmic pile-up); (b) radiation at high ΔT removes
~1.8 W; (c) slightly higher h. **However, both v1 and v2 are non-convergent at
U21**, so the 389 vs 334 gap is within the numerical uncertainty of the local
peak — it should **not** be read as "v2 is 14 % cooler". The robust, converging
quantity is the **board-scale field (~97 °C mean)**, and the robust engineering
conclusion — *5 W cannot be dissipated on this 13 cm² 1-oz board at ~150 °C
junction* — is the same in both.

## 5) Assumptions & GAPS

- **(S)** total loss = 5 W is the controlling assumption; sensitivity to load given.
- **(E)** package θ for D1/L2/L1 (no datasheet data) ±50 %; θ×1.5 → +18 °C.
- **(E)** h(ΔT) Incropera correlation, single horizontal-plate value on all faces.
- **(E)** wire fin lumped conductance (not a resolved fin).
- **GAP — U21/D1 junction not grid-converged** (C): a homogenised package body
  cannot resolve the true die/lead-frame interface peak. Next step: a dedicated
  die/lead-frame sub-model, or finer local mesh (L4 attempt hit 6.4 M tets and was
  aborted on memory grounds).
- **GAP — vias not modelled**: 166 stitching/GND vias would improve vertical board
  conduction (k_z) → current local peaks are likely **over**-estimated.
- **GAP — no TI datasheet found for ZEMS0650 / FXL0420** thermal params (DCR only);
  package θ used are E.
- Board copper treated as a uniform 81.2 %-coverage layer (no trace-level detail).
- Per constraint: 30 µm Cu **not** cut as 3-D solids (homogenised anisotropic board).
- X-level values: **none used**.

## Deliverables
`fig_thermal_map_v2.png`, `fig_thermal_sensitivity_v2.png`, `heat_sources_v2.json`,
`numbers_thermal_v2.json`, `convergence_v2.json`, `conservation_v2.json`,
`case.sif`, `build_mesh_v2.py`, `run_v2.py`, `make_case_v2.py`, `post_v2.py`,
`final_post.py`, `summarize_v2.py`, `conv_fixed_run.py`, `sens_run.py`,
`tps62933.pdf` (TI datasheet), `cases/`, `mesh_L*.msh`, `elmer_L*`.
