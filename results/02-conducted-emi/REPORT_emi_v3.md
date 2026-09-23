# TPS62933 DM conducted-EMI — v3 (real-geometry trace inductance)
## What changed vs v2
v2 used the generic **1 nH/mm** rule for PCB trace inductance in its pessimistic envelope (hot loop 20 nH, 30 µF-bank trace 8 nH, L1 branch 30 nH), giving only **+0.7 dB** worst-case margin at full load — flagged as over-pessimistic because the input copper is 3 mm wide.
v3 replaces those numbers with closed-form estimates from the **actual board geometry** parsed out of `epro/pourSim.epru`.
## Geometry (measured, not assumed)
| quantity | value | source |
|---|---|---|
| board outline | 26.50 × 50.50 mm | `pourSim.epru` POLY layer 11 |
| stack-up | 2-layer, FR4 core 59.449 mil = **1.510 mm**, εr=4.5 | LAYER_PHYS |
| bottom copper | **solid GND plane** over whole board | POUR6 (layer 2, net GND) |
| input copper width | **3.00 mm** (24VIN 118 mil, $1N22 119 mil, PPHV 118 mil) | POUR extents |
| trace-over-plane L′ | **0.2984 nH/mm** @ 3 mm width (microstrip) | Hammerstad-Jensen |

Formulas (`traceL.py`): microstrip Z₀ (Hammerstad-Jensen) + TEM relation L′ = Z₀·√εeff/c; cross-checked against the Rosa strip-over-plane partial inductance.
## Geometric trace inductances
| knob | geometric | v2 nominal | v2 pessimistic |
|---|---|---|---|
| lhot (hot loop) | **0.83 nH** | 8 | 20 |
| ltr_capA (30 µF bank) | **0.75 nH** | 3 | 8 |
| ltr_loop (L1 branch) | **7.88 nH** | 15 | 30 |
| lf1 (F1/input) | **6.32 nH** | 5 | 10 |

## Results (CISPR 32 Class B, DM, 50 Ω LISN, tr=tf=10 ns)
| case | V1@fsw [dBµV] | worst margin [dB] |
|---|---|---|
| noload | -16.4 | **+63.6** @ 12.88 MHz |
| half | 7.2 | **+45.9** @ 4.03 MHz |
| full | 13.2 | **+39.9** @ 4.03 MHz |

### Envelope, full load (Ipk = 3 A)
- **geo-optimistic (½ L)**: +35.3 dB @ 0.81 MHz
- **geo-nominal (measured)**: +39.9 dB @ 4.03 MHz
- **pessimistic (1 nH/mm, 100 Ω)**: +0.7 dB @ 4.03 MHz

## Answer
**With 3.00 mm-wide input copper, the full-load DM worst-case Class-B margin is ≈ +39.9 dB (real geometry).** Across optimistic..pessimistic the envelope is +39.9 .. +0.7 dB. The v2 "+0.7 dB" figure was an artefact of the 1 nH/mm trace assumption; using the measured 3 mm / 1.51 mm microstrip geometry (L′ ≈ 0.30 nH/mm) raises the margin by ~+24.0 dB.
## Caveats
- Trace lengths use component-anchor coordinates × 1.3 routing factor; the hot-loop and bank-stub are Euclidean (no factor). Values are ±30 % order-of-magnitude estimates, not extracted parasitics.
- Assumes the return current runs in the continuous bottom GND plane (microstrip). If the plane were broken under the input path, L′ would rise toward the 1 nH/mm bound.
- Component-to-model-node assignment (30 µF bank = C17/C62/C63; local caps = C72/73/74) is inferred from the v2 model structure + proximity; the schematic netlist was not re-traced here.
- ESR/ESL remain datasheet/typical values (unchanged from v2).
- Note: the geo-optimistic case (lower ESR/ESL) shows +35.3 dB, slightly WORSE than geo-nominal (+39.9 dB): reduced damping raises the parallel-filter peak at 0.81 MHz (fsw tank resonance). Lower loss ≠ lower emission at the tank resonance.
- geo_numbers.txt holds the full traceL geometry report; emi_numbers_v3.txt holds the EM run log.
