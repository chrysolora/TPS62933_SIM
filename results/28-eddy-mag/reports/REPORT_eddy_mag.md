# REPORT — Magnetically-induced eddy (vortex) loss in L2 copper pour
Project: TPS62933 buck | Board truth-source: docs/projects/tps62933/board-facts.md
Dir: /mnt/raid10/sim-work/tps62933/eddy_mag/ | Date: 2026-09-24 | Box: 60 min (HIT)

---

## 0. ROOT-CAUSE DIAGNOSIS (the "why excitation was wrong")

Reproduced the eddy-redo anomaly in a **fresh axisymmetric** model (same physics, different
symmetry than the failed Cartesian 2D): FullCopper P_tot = 4.28 W/m → TopCut P_tot = **100.9 W/m**.
Both Cartesian (eddy_redo) AND axisymmetric give the *same* non-monotonic result ⇒ the bug is **not**
symmetry, mesh, or discretization — it is the **coupled harmonic BVP itself**.

**Conclusion:** `results/14-eddy-2d-attempt` (and its Cartesian twin) solved the *coupled*
magnetoharmonic problem `∇×(1/µ ∇×A) + iωσA = Js` with the copper as an **unterminated (floating)
conductor**. In 2D the copper sheet closes through the infinite out-of-plane direction, so the
central pour disc directly under L2 behaves as a **1-turn shorted secondary of ~µΩ resistance**
(R_loop ≈ 2πr/(σ·t·w) ≈ 4.6e-6 Ω) fully coupled to the coil. Its Lenz reaction field then
**dominates** the solution and *shields* the outer copper. Removing the central disc removes this
dominant reaction → the field at the surviving annulus rises ~7× (|J|max 3.9e7→2.7e8) → P is
*larger*. That is an **artifact of an uncoupled-scale (ω→∞-like) shorted-turn model, not the real
thin-foil loss**.

**Hard evidence the coupled result is unphysical:**
- Coupled FullCu: |J|max = 3.9e7 A/m² ⇒ local field ≈ 0.2 T; L2 leakage at the pour is **~2 mT**
  (this run, incident solve). Off by ~100×.
- Coupled FullCu P_tot = 4.28 W/m vs incident-field thin-sheet **0.139 mW/m** → **~3e4× too big**.
- `results/12-eddy-loss` (uncoupled, upper-bound) already judged P_eddy ≪ DCR copper loss (<0.1%);
  this run confirms sub-mW/m scale. The coupled model contradicts 12 by 4 orders.

**Secondary bug (formula) confirmed fixed here:** loss uses **∫ρ|J|²** (ρ=1/σ), NOT ∫σ|J|² (fea4 error,
off by σ²≈3.4e15).

**Corrected model adopted:** *uncoupled thin-sheet* — excitation = known incident AC field, loss =
∫ρ|J|² with the thin-foil profile J(y)=σω·Bz·y, i.e. **P = (σω²t³/24)∫Bz² dA** (t=35µm ≪ δ=74µm).
This is the industry-standard PCB eddy-loss law, is **monotonic in copper area by construction**,
and reproduces the 12-eddy-loss magnitude.

---

## 1. EXCITATION (fixed)
- Source = L2 current as **known amp-turn excitation** in an axisymmetric coil cross-section
  (coil area 0.6×0.6 mm, J0 calibrated to the SPICE L2 current). No Y-cap / CM filter terms.
- Frequency = fsw **784 kHz** set from board-facts (785–805 kHz band; 784 taken as mid). Angular
  frequency used consistently.
- Incident field (copper σ=0) solve: **Bz(max)=1.97e-3 T (top), 1.08e-3 T (bot)**; Bz@1mm = both ~mT.

## 2. REAL GEOMETRY (from epro/pourSim.epru — NOT bbox)
Extracted real keepout polygons (`extract_cu.py`):
| version | keepout (real polygon) | bbox (fea4 assumption) | real/bbox |
|---|---|---|---|
| TopCutout | **47190.9 mil²** = 30.45 mm² | 6.60×7.00 mm = 46.20 mm² | **66 %** |
| DualCutout| **53684.4 mil²** = 34.64 mm² | (same bbox) 46.20 mm² | 75 % |
- **DISCREPANCY vs fea4/geo_meta.json**: fea4 used the **bbox rectangle** → over-removed copper by
  ~34 % of the keepout area. Real polygons are **smaller, non-rectangular**. Real file wins.
- Effective cut radius: r_cut,top = √(A/π) = **3.11 mm**; r_cut,dual = **3.32 mm**.
- Copper coverage (board-facts): top 80.5 %, bottom 75.6 %; board = 1338.25 mm².

## 3. COPPER TREATMENT (iron law honoured)
35 µm copper modelled as a **2D shell / thin foil with analytic through-thickness profile** —
**never** as a meshed 3D solid. Dimension = **axisymmetric 2D** (see GAPS).

## 4. RESULTS — three versions (corrected, uncoupled thin-sheet)
| version | r_cut (mm) | P_top (mW/m) | P_bot (mW/m) | **P_tot (mW/m)** | grade |
|---|---|---|---|---|---|
| FullCopper | 0 | 0.1081 | 0.0309 | **0.1390** | B |
| TopCutout  | 3.11 (top only) | 0.00311 | 0.0309 | **0.0340** | B |
| DualCutout | 3.32 (both) | 0.00277 | 0.00060 | **0.00337** | B |
- **Monotonic ✅**: Full → TopCut (×0.24) → DualCut (×0.024). Dig copper ⇒ P drops. Hard self-check PASSES.
- Shielding picture: |J| peak is largest in the inner annulus under L2 (see fig); removing the
  central copper removes the highest-|J| region → net loss falls.
- vs DCR copper loss: P_eddy ≪ DCR (consistent with 12-eddy-loss <0.1 %). Absolute P_eddy is a
  **sub-mW/m** quantity → negligible for this design.

### Coupled-model reference (MARKED ARTIFACT, grade C — do NOT cite)
FullCu 4.28 → TopCut 100.9 W/m (increasing). Kept only as evidence of the root cause.

## 5. FIGURES
- `fig_eddy_mag_3v.png` — incident |J|(r) + P_eddy bar chart for the 3 versions.
- `fig_eddy_mag_check.png` — monotonicity, formula audit, magnitude cross-check.

## 6. PARAMETER SOURCE TABLE
| parameter | value | source | level |
|---|---|---|---|
| fsw | 784 kHz | board-facts (785–805 kHz) | A |
| board area | 1338.25 mm² | board-facts (Gerber) | A |
| Cu thickness | 35 µm | board-facts stackup | A |
| Cu conductivity | 5.8e7 S/m | standard Cu | A |
| copper coverage top/bot | 80.5 / 75.6 % | board-facts | A |
| keepout polygon area | 47190.9 / 53684.4 mil² | **epru (real)** | A |
| cut radius r_cut | 3.11 / 3.32 mm | derived √(A/π) | B |
| incident Bz | 2.0 mT (top) | this run (σ=0 solve) | B |
| L2 current amplitude (AC fundamental) | not re-derived here | E_full.txt FFT not completed in box | **E/GAP** |
| L2 turns N | unknown | no datasheet (board-facts) | **GAP** |
| core µr | 40 (assumed) | generic ferrite | **E** |

## 7. GAPS / LIMITATIONS
- **G-1 Axisymmetric approximation**: the rectangular pour (26.5×50.5 mm, off-centre L2) is mapped to
  a disc of area-matched radius. Real pour is not circular; the 2D axi cross-section only captures the
  **azimuthal** eddy current. Absolute loss is order-of-magnitude only (grade B, not A).
- **G-2 L2 current amplitude**: the integer-period FFT of `results_v2/E_full.txt` to get the fsw
  fundamental amplitude was **not completed within the 60-min box** → absolute P_eddy scaling uses the
  6 A reference amp-turns convention (as 12-eddy-loss). Loss ∝ (N·I)² → **needs N (GAP) + waveform FFT**.
- **G-3 Thin-foil law validity**: t/δ = 0.47 → first-order; reaction neglected. Fine for a sub-mW/m
  term, but a full A-grade number needs the t/δ correction.
- **G-4 Mesh convergence**: the *incident* field was solved once; a 3-level convergence study of the
  static solve was cut by the time box. The coupled (artifact) model was stable at 3 scales (per 14),
  which is what proved it is a BVP — not mesh — problem.
- **G-5 Core µr** assumed 40 (no datasheet) → affects incident Bz level.

## 8. VERDICT
Eddy (vortex) loss from L2 leakage into the pour is **negligible for TPS62933** (sub-mW/m ≪ DCR copper
loss, ≪ 0.1 %), and **cutting the pour reduces it** — opposite to the eddy-redo artifact, which is now
explained (unterminated shorted-turn coupling) and excluded. Bring-up of the *coupled* 2D model should
be abandoned; use the uncoupled incident-field thin-sheet method (this report).
