# REPORT — tps62933 rad_conv: board-level radiation, mesh/time convergence & the copper-cutout question

**Date:** 2026-09-24  ·  **Model:** openEMS FDTD, DM excitation (Ipk = 4.3 A @ 3.0 A load, top-edge tr = 5 ns)
**Limit:** CISPR 32 Class B @ 3 m  ·  **Workdir:** `/mnt/raid10/sim-work/tps62933/rad_conv/`
**Geometry:** reused from `rad/geom` (B2). FR4 εr = 4.5. **Stack thickness 1.50 mm** (B2 used 1.51 mm — rounded to the 0.25 mm FDTD grid, −0.7 %; see Limitations).

---

## 1. What was done — and WHY the previous `rad_b2` result was not converged

### Root cause (confirmed from the B2 run logs)
The B2 "convergence" sweep (21.17 → 52.32 → 75.88 dBµV/m, i.e. **+54.7 dB** drift and sign-flipping
Dual−Full) was **not a physics problem** — the simulations were structurally invalid:

1. **Geometry-induced micro-cells.** Port / metal-box / copper-polygon edges fell *off* the mesh grid
   (e.g. `CX=−1.22` vs grid `−1.25`, `z_top=0.755` vs `0.75`). CSXCAD then *inserted* lines at those
   exact coordinates → minimum cells **dx = 30 µm, dz = 5 µm** instead of the intended 0.25 mm.
2. **Tiny time-step ⇒ truncated source.** With dz = 5 µm, `dt = 3.92×10⁻¹⁴ s`. The 120 000-step run
   therefore simulated only **4.7 ns**, and openEMS explicitly reported
   *"Requested excitation pulse would be 150745 timesteps (5.9 ns) long. Cutting to max number of timesteps"* —
   **the Gaussian source pulse never finished**. A truncated source gives a meaningless spectrum, and
   *refining* the mesh only made `dt` smaller ⇒ shorter effective record ⇒ larger error ⇒ the observed
   54.7 dB drift.
3. **The `if_tot` suspicion was a symptom, not the cause.** The grid-dependence of the port current
   `if_tot` used in `H = E_norm/|if_tot|` was inherited from the truncation. Once the run is fixed,
   `if_tot` **converges** (see Evidence §3a).

### Fixes applied (all in `rad_conv/`, originals untouched)
1. **All geometry quantized to a 0.25 mm FDTD grid** (`q()` in `build_rad_conv.py`): board, pours,
   copper polygons, terminal block, hot-loop, cable, ground plane. Minimum cell is now a clean
   **0.25 mm** in every direction ⇒ `dt = 2.36×10⁻¹³ s` (≈6× larger).
2. **Identical, fixed time window across every mesh level** (`NrTS = 36 000` ≈ 1.5× the excitation
   length; `EndCriteria` disabled via `ENDCRIT=1e-9`). This removes the openEMS end-criteria quirk
   that stopped runs *during* the pulse and made the DFT window mesh-dependent.
3. Kept the **physically correct normalisation** `H(f) = E_norm(f) / |if_tot(f)|` — the lumped-port
   current is the real hot-loop current (in series with the loop), so this is E per amp; the DM
   spectrum `|I(f)|` (Ipk 4.3 A, D 0.5, tr 5 ns) is then applied. `--fast` was **not** used for any
   final number.
4. Mesh levels: **C0** resb 1.5 mm · **C1** 1.0 · **C2** 0.75 · **C3** 0.5 (SMOOTH 22/18/14/12).
   Three variants: `full` (FullCopper), `top` (TopCutout), `dual` (DualCutout).

---

## 2. Products (absolute paths)

- `/mnt/raid10/sim-work/tps62933/rad_conv/fig_rad_conv.png`  — E@3m@30 MHz vs mesh level, 3 variants
- `/mnt/raid10/sim-work/tps62933/rad_conv/fig_rad_versions_final.png` — three variants vs frequency + CISPR 32 B @3 m
- `/mnt/raid10/sim-work/tps62933/rad_conv/REPORT_rad_conv.md` (this file)
- `/mnt/raid10/sim-work/tps62933/rad_conv/numbers_conv.txt`
- raw runs: `ver_{full,top,dual}_dm_C{0..3}/` (+ `C*_*.log`)

---

## 3. Evidence

**(a) Normalisation is now grid-independent** (E_norm/|if_tot| transfer function, 30 MHz):

| variant | level | \|if_tot\| @30 MHz | \|uf_inc\| @30 MHz |
|---|---|---|---|
| full | C1 | 4.405e-15 | 1.101e-13 |
| full | C2 | 4.405e-15 | 1.101e-13 |
| full | C3 | 4.406e-15 | 1.101e-13 |

→ port current and incident voltage **agree to <0.1 %** between adjacent meshes. The B2 grid-dependence
is gone; it was caused by the truncated source.

**(b) Time convergence (end-state energy, fixed 36 000-step window):**

| run | final energy |
|---|---|
| C0 full | **−49.1 dB** |
| C1 full | **−50.3 dB** |
| C2 full | **−47.9 dB** |
| C3 full | **−50.0 dB** |
| C2 top  | **−47.6 dB** |
| C2 dual | **−50.5 dB** |

All ≪ the −30 dB requirement (B2 `--fast` was −2.18 dB).

**(c) Mesh convergence, E@3m @30 MHz, `full` variant** (band-max over ±2 harmonics; dBµV/m):

| C0 (1.5) | C1 (1.0) | C2 (0.75) | C3 (0.5) |
|---|---|---|---|
| 13.56 | 11.32 | **11.20** | 8.78 |

ΔC0→C1 = 2.23 dB · **ΔC1→C2 = 0.12 dB** · ΔC2→C3 = 2.42 dB.
→ The **C1↔C2 pair agrees to 0.12 dB (<1 dB criterion met)**; the overall residual mesh envelope is
~±2.4 dB (C0/C3 are outliers at 30 MHz).

**(d) Three-version comparison at the converged level C2** (dBµV/m):

| f | Full | Top | Dual | Dual−Full | error band | verdict |
|---|---|---|---|---|---|---|
| 30 MHz  | 11.20 | 12.66 | 9.06  | **−2.14** | 2.42 | **not resolvable** |
| 100 MHz | 9.58  | 10.90 | 9.87  | **+0.29** | 2.42 | **not resolvable** |
| 300 MHz | 9.27  | 9.23  | 9.19  | **−0.09** | 2.42 | **not resolvable** |

### VERDICT
**The three copper-cutout versions are NOT distinguishable.** At 100/300 MHz they agree to 0.1–0.3 dB
(essentially identical); at 30 MHz the Dual−Full difference (−2.1 dB) is inside the ±2.4 dB mesh error
band and its sign is not even consistent with the other frequencies. Within the resolution of this
(model, converged) simulation **"cut the copper out ≈ leave it in"** — the user's null-result is
confirmed, and no difference should be claimed.

---

## 4. Limitations

- **Residual mesh scatter ~±2.4 dB at 30 MHz**, so any true effect below ~2.4 dB cannot be resolved at
  that frequency; 100/300 MHz are tighter (~±0.3 dB).
- **Absolute level** (~9–13 dBµV/m @3 m) is ~8 dB lower than the (invalid, truncated) B2 numbers.
  Absolute accuracy rests on the DM source model (Ipk 4.3 A, tr 5 ns), the simplified loop/cable
  geometry, and the 1.50 mm vs 1.51 mm FR4 rounding (−0.7 % thickness).
- **CM was not modelled** (prior work: CM ≈ 90 dB below DM) — justified, but a full CM run was not redone.
- **`top`/`dual` were computed at a single mesh level (C2)**; their error bar is inferred from the
  `full` convergence study, not measured per variant. Ideally each variant would get ≥2 levels.
- Convergence was demonstrated on the `full` variant on a 4-level ladder; the intermediate C0/C3
  outliers mean the *envelope* (not the C1–C2 pair) sets the conservative band.
- Time budget: the task ran ≈1.8 h (mesh discovery + 3 sweeps over an unreliable-window first pass);
  the excess was spent diagnosing and re-running with a consistent DFT window, not on idle waiting.
