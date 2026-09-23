# REPORT — TPS62933 buck : frequency-domain DM conducted EMI (CISPR 32 approx.) + soft-start

**Design under study:** TPS62933F, 24 V → 12 V, Imax 3 A, fsw ≈ 805 kHz (RT = 26.7 kΩ).
**Simulated device model:** unencrypted *TPS62933P* transient model (different variant, no SS pin).
**Host:** node249 (`/mnt/raid10/sim-work/tps62933/emi/`), ngspice 42 + python3 (numpy 1.26 / scipy 1.11 / matplotlib 3.6, headless Agg).

---

## 1. Deliverables

| file | content |
|---|---|
| `emi_model.py` | complete, re-runnable frequency-domain evaluation script |
| `fig_emi_cispr.png` | DM conducted-emission spectrum (150 kHz–30 MHz) for 3 load cases + CISPR 32 Class B & A QP limits |
| `fig_startup.png` | analytic soft-start ramp (0→12 V, Tss = 14 ms) |
| `REPORT_emi.md` | this document |
| `emi_numbers.txt` | raw console log of all numbers quoted below |

Re-run command:
```bash
ssh node249
cd /mnt/raid10/sim-work/tps62933/emi
python3 emi_model.py          # ~5 s, no GUI, writes both PNGs + emi_numbers.txt
```

---

## 2. Why a frequency-domain method (not a re-run of the transient)

The real board input filter (rebuilt pin-by-pin from the epru schematic, confirmed) is

```
24VIN --F1--> PPHV --[C17||C62||C63 = 3x10uF]-- L1(1uH) --> PPHV_OUT_FILTER
                                   |                              |
                                 TVS D138                  C73||C74=2x10uF, C72=100nF,
                                                            R85(100m)+C16||C64 damping
                                                                   |
                                                              U21.VIN (TPS62933F)
U21.SW --L2(15uH)--> 12V  (C1 = output cap, value UNKNOWN)
```

Two hard constraints force the frequency-domain approach:

1. **The unencrypted model will not start** with the true >30 µF directly at the VIN pin
   (stage-1 had to shrink the input π to 20 µF–1 µH–10 µF; the real 3×10 µF bank never entered the transient).
2. **The `results_v2/E_*.txt` transient is not at a settled operating point.** Evidence:
   Vout is still drifting 10.6 → 11.5 V across the 80 µs run, and the measured switch duty
   (≈ 0.29 from both `v(SW)>Vin/2` and `d(i_L2)/dt` sign) is inconsistent with Vout ≈ 12 V in CCM
   (which would need D ≈ 0.5). The L2 current is still ~3.9 A instead of 3.0 A (output cap still charging).
   → the transient cannot be trusted as a *steady-state* EMI source.

What the transient **does** give reliably: the switching frequency.
`fsw` measured from the median of `i(L2)` slope-reversal periods = **805.0 kHz** (period 1.2423 µs),
consistent across the full and half-load cases and with the task's measured value.

So: **source = analytic trapezoid, filter = linear frequency-domain network.**
The data-derived spectrum is computed as an independent, documented cross-check (§6).

---

## 3. Method

### 3.1 DM noise source
The DM source is the pulsating current the converter draws at `U21.VIN` (the high-side switch current).
In CCM this is a **trapezoidal pulse train** at `fsw`:

* duty `D = Vout/Vin ≈ 0.50` (ideal buck; conduction losses push it slightly above 0.50 — ignored),
* peak `Ipk = I_L(DC) = I_load` (inductor average = load current in CCM),
* rise/fall `tr = tf = 15 ns` (**estimate** — not available from the model; see §5).

Cases: `full` Ipk = 3.0 A, `half` Ipk = 1.5 A, `noload` Ipk = 0.22 A, D = 0.15
(no-load runs PFM/burst — the trapezoid is only a coarse stand-in there, see §5).

Harmonics are obtained **numerically** (FFT of one synthesised pulse train, DC removed, Hann window →
amplitudes 2|X|/Σw). Using the numeric FFT instead of the closed-form sinc formula avoids any
normalisation/decade mistake. First-harmonic amplitudes: full 1.909 A, half 0.955 A, noload 0.064 A.

### 3.2 Input filter + LISN — transimpedance Z_t(f)
Linear 2-node nodal analysis at each harmonic frequency (nodes as in §2, fuse **F1 treated as an ideal short**):

* **node A** = LISN port = PPHV: `3×10µF` bank, LISN `50Ω` (via 1 µF coupling), LISN `50µH` to mains (AC gnd)
* **L1 = 1 µH** between A and B
* **node B** = PPHV_OUT_FILTER = U21.VIN (source node): `2×10µF`, `100nF`, damping `R85→100 mΩ + 2×10µF`

The noise current is injected at node B (unit current source), and we solve for the voltage at node A:

```
Z_t(f) = | V_LISN(f) / I_src(f) |
V_lisn(harmonic k) = I_src(k·fsw) · |Z_t(k·fsw)|        [V, rms]
EMI_dBuV = 20·log10( V_lisn / 1µV )
```

### 3.3 Limits (CISPR 32, QP, dBµV)
* Class B: 0.15–0.5 MHz 66→56 (linear in log f), 0.5–5 MHz 56, 5–30 MHz 60
* Class A: 0.15–0.5 MHz 79→73, 0.5–5 MHz 73, 5–30 MHz 73

### 3.4 Soft-start (analytic)
The datasheet-recommended relation is `Tss ≈ Css·Vref / Iss`. With the board value `Css = 100 nF`
and `Vref = 0.8 V`, the user-measured **Tss = 14 ms** is taken as authoritative; the implied
charge current is `Iss = Css·Vref/Tss = 0.571 µA… = 5.71 µA`. The curve drawn is a straight 0→12 V
ramp over 14 ms (soft-start ramps the reference, which maps linearly to Vout).

---

## 4. Results

### 4.1 DM spectrum vs CISPR 32 Class B  (`fig_emi_cispr.png`)

Transimpedance magnitude (filter + LISN): `|Z_t(0.805 MHz)| = 6.3 µΩ`, `|Z_t(5 MHz)| = 6.6 µΩ`,
`|Z_t(30 MHz)| = 19.7 µΩ`. The filter is a very effective DM shunt in the whole CISPR band.

| case | Iload | D | I_src(1st) | V_lisn(1st) | margin_min vs Class B | violations |
|---|---|---|---|---|---|---|
| noload | 0.02 A | 0.15 | 0.064 A | −7.9 dBµV | **+63.9 dB** @ 0.805 MHz | 0 |
| half   | 1.50 A | 0.50 | 0.955 A | 15.6 dBµV | **+40.4 dB** @ 0.805 MHz | 0 |
| full   | 3.00 A | 0.50 | 1.909 A | 21.6 dBµV | **+34.4 dB** @ 0.805 MHz | 0 |

* The spectrum is a **harmonic line spectrum at k·805 kHz** — the first in-band harmonic is at 0.805 MHz;
  150 k–805 k has no harmonic (no DM emission there).
* Peak of every curve is the **fundamental (0.805 MHz)** (largest harmonic of a D≈0.5 trapezoid),
  then the envelope rolls off — the expected shape.
* **All three DM cases pass CISPR 32 Class B with ≥ 34 dB margin** (Class A by ≥ 46 dB).

### 4.2 Resonances / filter features
* `L1 × C_series` resonance: `C_ser = 30µF ∥ 40.1µF = 17.2 µF`, `f = 1/(2π√(L1·C_ser)) = 38.4 kHz`
  — **below** the 150 kHz CISPR band start, so within the band the filter is capacitive/monotonic.
* In-band capacitor self-resonances (ESL-dominated): 30 µF bank ≈ 1.45 MHz, 40 µF bank ≈ 1.03 MHz,
  100 nF ≈ 15.9 MHz. Above these the caps turn inductive and `|Z_t|` rises — visible as the
  flattening/upturn of the curves toward 30 MHz.

### 4.3 Rough CM cross-check (order-of-magnitude only)
With an **assumed** stray capacitance `Cp = 5 pF` coupling the SW node to earth and the CM LISN
return impedance ≈ 25 Ω:

| harmonic | f | V_sw_h | I_cm | V_cm | Class B |
|---|---|---|---|---|---|
| 1 | 0.81 MHz | 15.3 V | 386 µA | 79.7 dBµV | 56 → **fail** |
| 3 | 2.42 MHz | 5.08 V | 386 µA | 79.7 dBµV | 56 → **fail** |
| 5 | 4.03 MHz | 3.04 V | 384 µA | 79.6 dBµV | 56 → **fail** |

CM scales **linearly with Cp**: pass needs Cp ≲ 0.5 pF, unrealistic on a real board.
→ **Engineering conclusion: DM is comfortably compliant; CM is expected to be the binding
constraint** and would set the need for CM mitigation (Y-caps / CM choke / shielded loop area).
*This is a crude estimate (see §5), not a validated CM model.*

---

## 5. Assumptions, estimates and limitations

**Assumptions / engineering estimates (explicitly flagged — values NOT given in the brief):**

| item | value used | status / impact |
|---|---|---|
| MLCC 10 µF/50 V ESR | 3 mΩ/cap | estimate; small effect at 805 kHz |
| MLCC 10 µF/50 V ESL | 1.2 nH/cap | estimate; sets the ~1–1.5 MHz bank SRF and the HF upturn |
| MLCC 100 nF ESR/ESL | 20 mΩ / 1.0 nH | estimate |
| L1 DCR / parasitic Cp | 30 mΩ / 1 pF | estimate; DCR negligible vs μΩ-level Z_t |
| switch rise/fall time | 15 ns | estimate; controls HF harmonic envelope (only matters ≫ 10 MHz here) |
| LISN representation | 50 Ω + 50 µH + 1 µF | **task-specified simplification** |
| fuse F1 | ideal short | ≈ 0 Ω in reality; adds a few nH → slightly *less* attenuation at HF |
| source Ipk | I_load (3.0 / 1.5 / 0.22 A) | nominal CCM value; the transient showed 3.9 A because it was *not settled* |

**Known limitations (do not treat as fact):**

1. **Model variant mismatch (P vs F).** The available model is TPS62933**P**; the board is TPS62933**F**.
   Dead-time, switch node voltage, and the control loop differ. The transient source is therefore
   *indicative only*; the EMI source used here is analytic, not device-extracted.
2. **DM only.** This is a differential-mode study. Common-mode is **not** modelled except for the
   crude order-of-magnitude note in §4.3. Real conducted compliance on a 24 V converter like this is
   normally CM-limited, so **absence of a DM violation does NOT imply EMC compliance.**
3. **`C1` (12 V output capacitor) value is UNKNOWN.** It does not enter the input-filter DM network
   (the input-side shunt path dominates at the LISN), so it does not change these numbers; it must
   **not** be assumed for output-ripple or loop-stability work.
4. **TVS D138 (SMF24A) off-state capacitance** is unknown and was **omitted** (not invented); it would
   add a small HF capacitance at PPHV, marginally changing the HF tail.
5. **Idealised, optimistic filter.** Ideal linear caps + estimated parasitics with no PCB trace
   inductance between stages. Real layout adds nH-level series inductance and higher effective ESR,
   which would **raise** DM emission above the values here. Treat the reported margins as an
   *upper bound*.
6. **Trapezoid source ignores switching ringing / gate-drive detail**; real waveforms have HF ringing
   that fills the trapezoid's spectral nulls (this is exactly why the data cross-check in §6 shows
   spikes at harmonics where the ideal trapezoid is ~0).
7. **The transient cross-check is not settled** (§2), so its harmonic amplitudes are only a sanity
   order-of-magnitude, not a validation.
8. **DM LISN impedance** is taken as a single 50 Ω per the brief. Strictly, for DM the two LISN 50 Ω
   appear in series (≈100 Ω), which would raise DM numbers by ~6 dB — still far from the limit.
9. **Soft-start is analytic**, not a simulation (the P model has no SS pin). The formula relation
   `Tss ≈ Css·Vref/Iss` is stated with the implied `Iss = 5.71 µA`; the **datasheet `Iss` was not
   available**, so this is an implication of the user's 14 ms, not a confirmed datasheet value.

---

## 6. Data cross-check (documented, indicative)

FFT of the reconstructed chopped input current `i(L2)·[v(SW)>Vin/2]` from `E_full.txt` / `E_half.txt`
(same uniform-grid + Hann normalisation) vs the analytic trapezoid:
median deviation **+1.9 dB (full)**, **+3.6 dB (half)** — i.e. the *levels* agree within a few dB,
but the data spread is huge (up to +150 dB at individual harmonics) because the ideal trapezoid has
spectral nulls that the non-ideal, unsettled waveform fills in. Conclusion: the analytic source is a
sound level estimate; the data cannot be used as a clean line spectrum for this board. The dotted
curves in `fig_emi_cispr.png` show this.

---

## 7. Bottom line

* fsw confirmed = **805 kHz**; DM spectrum peaks at fsw harmonics and rolls off as expected.
* The real input filter (30 µF | 1 µH | 20 µF+100nF+damping) gives **large DM margins
  (≥ 34 dB to CISPR 32 Class B) at all loads**; the main LC resonance (38 kHz) is below the band.
* **Caveat, emphasised:** this is an *idealised DM-only* result. The likely real-world problem is
  **common mode** (§4.3), and the DM numbers are an optimistic upper bound (§5.5–5.7).
* Soft-start: analytic 0→12 V ramp over 14 ms (`fig_startup.png`), consistent with `Css = 100 nF`.
