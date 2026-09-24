# TPS62933 board thermal 3×3 — RERUN with REAL source cutout polygons

**Dir:** `/mnt/raid10/sim-work/tps62933/thermal_3x3_real/`
**Date:** 2026-09-24  **Supersedes:** `thermal_3x3/` (rectangle keep-out) — see §4.
**Scope:** only the keep-out GEOMETRY changed. Losses, material params, BC, h-chain, θ, solver settings are byte-identical to `thermal_3x3/`. FullCopper mesh was **reused unchanged** (same geometry) → its numbers are an exact reproduction sanity check.

---

## 1. What was done

- Re-parsed the three copper versions **directly from the source** `epro/pourSim.epru` (`poly_from_epru.py`) and confirmed the TopCutout/DualCutout keep-out polygons vertex-by-vertex against `geom_verify/copper_areas.json` (**max vertex diff 0.000 mil, area diff 0.000 mm² → MATCH**, 4/4 regions).
- Built new gmsh meshes using the **true non-convex polygons** (not the bounding rectangles). TopCutout and DualCutout now have **different** polygons (their source regions differ); the old run used identical rectangles for both.
- Rebuilt the 3×3 matrix (3 geos × 3 loads 0.6/1.5/3.0 A) with the same tilted-natural-convection h(ΔT) outer iteration and the same θJA chain.
- Produced the rect-vs-real comparison and figures.

### Coordinate conversion (documented, matches old mesh frame)
`pourSim.epru` vertices are **mil, y down-positive**. The thermal mesh frame is mm, y up-positive, `x∈[0,26.50]`, `y∈[0,50.50]`:

```
x_mm = x_mil * 0.0254
y_mm = 50.50 + y_mil * 0.0254          # 50.50 = board Y size = 1988.19 mil * 0.0254
```

Verified against the old keep-out rectangles: real L2 bbox → **6.600 × 7.000 mm** centred (6.223, 33.736) mm = exactly the old KEEPOUT rect centre. So the *frame* is right and only the *shape* changed.

### Keep-out areas (real polygon vs old rectangle)
| region | real poly (mm²) | old rect bbox (mm²) | over-cut |
|---|---|---|---|
| L2 top  (PCB8_1, lyr1) | 30.446 | 46.20 (6.60×7.00) | **-34.1 %** |
| L1 top  (PCB8_1, lyr1) | 11.048 | 19.836 (4.35×4.56) | **-44.3 %** |
| L2 dual (PCB8_2, lyr12)| 34.635 | 46.20 | **-25.1 %** |
| L1 dual (PCB8_2, lyr12)| 13.464 | 19.836 | **-32.1 %** |

Meshes (main / coarse nodes): fullcu 34669/10403, topcut 36613/11164, dualcut 36527/11054 — each cutout mesh has exactly **2 keep-out bodies**.

---

## 2. New 3×3 result table (main mesh, Tamb 25 °C)

| geo | Iout (A) | board mean (°C) | U21 Tj cons (°C) | U21 Tj best (°C) | D1 max (°C) | L2 max (°C) | L1 max (°C) | h (W/m²K) | ΣQ imb. (%) | gate Tj<150 |
|---|---|---|---|---|---|---|---|---|---|---|
| FullCopper | 0.6 | 33.55 | 40.14 | 38.46 | 44.15 | 33.58 | 34.36 | 7.79 | -0.04 | PASS |
| FullCopper | 1.5 | 47.04 | 68.57 | 63.35 | 72.11 | 49.96 | 48.77 | 9.81 | -0.02 | PASS |
| FullCopper | 3.0 | 73.32 | 133.03 | 119.19 | 119.36 | 89.54 | 76.08 | 11.78 | -0.00 | PASS |
| TopCutout  | 0.6 | 33.53 | 40.21 | 38.52 | 44.21 | 33.64 | 34.41 | 7.79 | -0.30 | PASS |
| TopCutout  | 1.5 | 47.04 | 68.78 | 63.53 | 72.24 | 50.19 | 48.88 | 9.80 | -0.28 | PASS |
| TopCutout  | 3.0 | 73.39 | 133.51 | 119.58 | 119.59 | 90.39 | 76.26 | 11.78 | -0.20 | PASS |
| DualCutout | 0.6 | 33.54 | 40.50 | 38.72 | 44.28 | 33.93 | 34.53 | 7.79 | -0.04 | PASS |
| DualCutout | 1.5 | 47.01 | 69.61 | 64.09 | 72.35 | 51.54 | 49.14 | 9.81 | -0.19 | PASS |
| DualCutout | 3.0 | 73.25 | 135.54 | 120.91 | 119.64 | 95.27 | 76.72 | 11.78 | -0.10 | PASS |

### Final gates (all cases)
| gate | threshold | FullCopper | TopCutout | DualCutout | verdict |
|---|---|---|---|---|---|
| board-mean grid convergence (coarse vs main, 3 A) | < 0.5 % | 0.068 % | 0.204 % | 0.218 % | **PASS** |
| ΣQ conservation | < 2 % | -0.00 % | -0.20 % | -0.10 % | **PASS** |
| U21 Tj (conservative) | < 150 °C | 133.03 | 133.51 | 135.54 | **PASS** |

(U21 near-field ref convergence 0.30–1.17 %; Tj_best convergence 0.26–0.98 %.)

---

## 3. Rect(OLD) vs Real(NEW) — same-cell comparison

Δ = OLD_rect − NEW_poly (positive = old run was **too hot / overstated**):

| geo | load | board mean Δ | U21 Tj(cons) Δ | L2 max Δ | Tmax Δ |
|---|---|---|---|---|---|
| FullCopper | 0.6/1.5/3.0 | 0.00 | 0.00 | 0.00 | 0.00 |
| TopCutout  | 0.6 | +0.01 | +0.15 | 0.00 | +0.01 |
| TopCutout  | 1.5 | 0.00 | +0.44 | +0.05 | +0.02 |
| TopCutout  | 3.0 | -0.01 | +1.11 | +0.32 | +1.32 |
| DualCutout | 0.6 | 0.00 | +0.19 | -0.23 | +0.12 |
| DualCutout | 1.5 | -0.03 | +0.36 | +1.27 | +0.21 |
| DualCutout | 3.0 | -0.21 | +0.26 | **+7.42** | +0.42 |

**FullCopper reproduces the old run EXACTLY (Δ = 0.000 °C on every metric, all 3 loads)** → the new pipeline is on the same rails; only geometry differs where geometry differs.

### Quantified overstatement of "cutout → L2 hotspot" (3.0 A)
| metric | OLD rect | NEW real | overstated by |
|---|---|---|---|
| DualCut L2 max | 102.69 °C | 95.27 °C | **+7.42 °C (old 56 % too high vs the cutout-induced rise)** |
| DualCut L2 rise over FullCopper (89.54) | +13.15 °C | **+5.73 °C** | **old effect 2.3× too large** |
| TopCut L2 rise over FullCopper | +1.17 °C | +0.85 °C | 1.4× too large |
| DualCut L2 @1.5 A | 52.81 | 51.54 | +1.27 °C |

**Conclusion:** the headline claim **"dual-side copper cutout on L2 raises its temperature by +13.15 °C" is NOT valid**; the correct source-geometry number is **+5.73 °C at 3 A**. The board-mean temperature and U21 junction temperature are largely insensitive to the cutout shape (Δ ≤ 0.3 °C), so the **Tj<150 °C gate is unchanged and still passes**. Only the **L2 device hotspot** was materially overstated by the rectangle approximation.

---

## 4. Provenance / credibility (A/B/C)

| item | value | source | grade |
|---|---|---|---|
| Cutout polygon vertices | 14/13/14/11-vertex polygons | re-parsed from `epro/pourSim.epru`, cross-checked vs `geom_verify/copper_areas.json` (diff 0) | **A** |
| Coordinate conversion | x·0.0254, 50.50+y·0.0254 | validated against old KEEPOUT rect centre/bbox | **A** |
| Loss budget (D1 Vf=0.70, L2 DCR=0.085, L1 DCR=0.027, U21 RDS 0.076/0.032, U21_sw≈0.50 W) | unchanged from thermal_3x3 | datasheets (D) partly estimate (E) for U21_sw | **B** |
| board mean / Tj (conservative) numbers | per table | converged (<0.5 %), ΣQ<0.3 %, same as old | **A** |
| L2 / D1 / L1 body-max numbers | per table | mesh-sensitive device hotspot; converged; **but** see GAPS | **B** |
| "cutout → L2 hotspot +5.73 °C @3 A" | quantitative | direct consequence of A-grade geometry | **A** (relative), B (absolute) |
| "rect run overstated DualCut L2 by +7.42 °C" | quantitative | old vs new, same everything else | **A** |

---

## 5. Limitations / GAPS

1. **Keep-out is modelled as a full-thickness region with reduced in-plane conductivity** `k_xy` (topcut −1 layer, dualcut −2 layers), identical to the old model. This is a homogenised equivalent, not an explicit Cu/layer stack; the *placement and shape* of that region is now source-exact, but the *magnitude* of the effective k reduction is a modelling assumption (unchanged from `thermal_3x3`, so it does not affect the rect-vs-real delta). θ-chain for U21 (19.3 °C/W junction-to-ref) is a datasheet/estimate value, unchanged.
2. **Device body-max** is a mesh-sampled max over a box above the board; ±0.3–0.5 °C between mesh levels. Used only for the relative rect-vs-real comparison, where it is adequate.
3. **U21_sw loss (~0.50 W)** is an estimate (E); it dominates U21 Tj and therefore the gate is sensitive to it. Not re-tuned here (kept identical to old run on purpose).
4. No experimental validation; this is a modelling study.

---

## 6. Files
- `fig_thermal_3x3_real_matrix.png` — 3×3 field matrix (green outline = real cutout polygon)
- `fig_thermal_3x3_rect_vs_real.png` — old rect vs real polygon (L2/board/area/hotspot-rise)
- `numbers_thermal_3x3_real.json`, `convergence_3x3_real.json`, `heat_sources_3x3_real.json`
- `cutout_polygons.json` (source polygons + cross-check), `poly_from_epru.py`, `build_mesh_3x3_real.py`, `run_3x3_real.py`, `post_3x3_real.py` (+ `make_case_3x3.py` unchanged)
- `cases/<geo>_L##/case.sif` — Elmer inputs
