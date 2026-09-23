#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emi_model_v3.py -- TPS62933 DM conducted-EMI, VERSION 3
     = same nodal DM model as v2, but the PCB trace/loop inductances now come
       from the REAL GEOMETRY (emi_v3/traceL.py), not the generic 1 nH/mm rule.

Run: cd /mnt/raid10/sim-work/tps62933/emi_v3 && python3 emi_model_v3.py
Out: fig_emi_cispr_v3.png, emi_numbers_v3.txt, REPORT_emi_v3.md, geo_numbers.txt

The v2 model's "pessimistic" envelope used 1 nH/mm for the hot loop (20 nH) and
the 30uF bank trace (8 nH), giving only +0.7 dB at full load.  The board's input
copper is 3.00 mm wide (measured) over a solid GND plane 1.51 mm below -> the
geometric values are lhot=0.83 nH, ltr_capA=0.75 nH, ltr_loop=7.9 nH, lf1=6.3 nH.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceL                     # sibling module: geometric inductance

OUT = '/mnt/raid10/sim-work/tps62933/emi_v3/'
FSW = 805e3
G = traceL.main()                 # prints geo report, writes geo_numbers.txt, returns dict

# ---------------------------------------------------------------- model (from v2) ----
def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L
def capZ(f, C, esr, esl): return esr + Zl(f, esl) + Zc(f, C)

def trapezoid_harmonics(fsw, ipk, duty, tr, tf, kmax, nper=8, ppp=16384):
    T = 1.0 / fsw; N = nper * ppp; dt = (nper * T) / N
    t = np.arange(N) * dt; tt = t % T; on = duty * T
    x = np.zeros_like(t)
    if tr > 0:
        r = tt < min(tr, on); x[r] = ipk * tt[r] / tr
    fl = (tt >= tr) & (tt < on); x[fl] = ipk
    if tf > 0:
        fa = (tt >= on) & (tt < on + tf)
        x[fa] = ipk * (1 - (tt[fa] - on) / tf)
    x -= x.mean()
    w = np.hanning(N)
    X = np.fft.rfft(x * w); f = np.fft.rfftfreq(N, dt)
    amp = 2 * np.abs(X) / w.sum()
    ks = np.arange(1, kmax + 1); fh = ks * fsw
    idx = np.array([int(round(fr / f[1])) for fr in fh])
    return fh, amp[idx]

def solve(n, branches, inj):
    Y = np.zeros((n, n), complex); b = np.zeros(n, complex)
    for (i, j, Z) in branches:
        if abs(Z) < 1e-15: Z = 1e-12
        y = 1.0 / Z
        if i >= 0: Y[i, i] += y
        if j >= 0: Y[j, j] += y
        if i >= 0 and j >= 0: Y[i, j] -= y; Y[j, i] -= y
    b[inj] += 1.0
    return np.linalg.solve(Y, b)

def damping(P, fr):
    zc = capZ(fr, P['Ce'], P['esr_e'], P['esl_e'])
    return P['r85'] + 1.0 / (1.0 / zc + 1.0 / zc)

def Zt(fr, P):
    b = []
    b.append((0, -1, P['rs_dm'] + Zc(fr, P['c_lisn'])))
    b.append((0, -1, Zl(fr, P['l_mains'])))
    # nodes 0=A(LISN), 1=P(30uF bank), 2=V(VIN pin), 3=B(local caps)
    b.append((0, 1, P['rf1'] + Zl(fr, P['lf1'])))
    for _ in range(3):
        b.append((1, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'] + P['ltr_capA'])))
    b.append((1, 2, P['dcr_l1'] + Zl(fr, P['L1'] + P['ltr_loop'])))
    if P['cp_l1'] > 0: b.append((1, 2, Zc(fr, P['cp_l1'])))
    b.append((2, 3, Zl(fr, P['lhot'])))
    for _ in range(2):
        b.append((3, -1, capZ(fr, P['C10u'], P['esr_10u'], P['esl_10u'])))
    b.append((3, -1, capZ(fr, P['C100n'], P['esr_100n'], P['esl_100n'])))
    b.append((3, -1, damping(P, fr)))
    V = solve(4, b, 2)
    return abs(V[0])

def Zt_vec(f, P):
    return np.array([Zt(fr, P) for fr in np.atleast_1d(f)])

def cispr_B(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    66.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 56.0, 60.0))
def cispr_A(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    79.0 - 6.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 73.0, 73.0))

# ---------------------------------------------------------------- parameter sets ----
_BASE = dict(C10u=10e-6, C100n=100e-9, Ce=10e-6,
             L1=1e-6, dcr_l1=25e-3, cp_l1=4e-12,
             rs_dm=50.0, c_lisn=1e-6, l_mains=50e-6)

def P_v1repro():
    d = dict(_BASE)
    d.update(esr_10u=3e-3, esl_10u=1.2e-9, esr_100n=20e-3, esl_100n=1.0e-9,
             esr_e=0.0, esl_e=0.0, r85=0.1, dcr_l1=30e-3, cp_l1=1e-12,
             lf1=0.0, ltr_loop=0.0, lhot=0.0, ltr_capA=0.0, rf1=0.0)
    return d

def P_geo_nominal():
    """REAL GEOMETRY trace L (traceL.py) + nominal ESR/ESL, 50ohm LISN."""
    d = dict(_BASE)
    d.update(esr_10u=2e-3, esl_10u=1.0e-9, esr_100n=20e-3, esl_100n=0.8e-9,
             esr_e=4.0, esl_e=3.0e-9, r85=0.1,
             lf1=G['lf1'] * 1e-9, ltr_loop=G['ltr_loop'] * 1e-9,
             lhot=G['lhot'] * 1e-9, ltr_capA=G['ltr_capA'] * 1e-9, rf1=5e-3)
    return d

def P_geo_optimistic():
    """geometric L halved (shorter routing / ideal image return) + low parasitics."""
    d = dict(_BASE)
    d.update(esr_10u=1e-3, esl_10u=0.7e-9, esr_100n=15e-3, esl_100n=0.6e-9,
             esr_e=2.0, esl_e=2.0e-9, r85=0.1,
             lf1=0.5 * G['lf1'] * 1e-9, ltr_loop=0.5 * G['ltr_loop'] * 1e-9,
             lhot=0.5 * G['lhot'] * 1e-9, ltr_capA=0.5 * G['ltr_capA'] * 1e-9, rf1=3e-3)
    return d

def P_pessimistic():
    """v2 pessimistic: generic 1 nH/mm traces (hot loop 20nH, bank 8nH) + 100ohm LISN."""
    d = dict(_BASE)
    d.update(esr_10u=6e-3, esl_10u=1.5e-9, esr_100n=40e-3, esl_100n=1.2e-9,
             esr_e=8.0, esl_e=5.0e-9, r85=0.1, dcr_l1=30e-3, cp_l1=6e-12,
             lf1=10e-9, ltr_loop=30e-9, lhot=20e-9, ltr_capA=8e-9, rf1=10e-3,
             rs_dm=100.0)
    return d

CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]
KMAX = int(30e6 / FSW)
TR = TF = 10e-9

def spectrum(P, ipk, duty, tr=TR, tf=TF):
    f, a = trapezoid_harmonics(FSW, ipk, duty, tr, tf, KMAX)
    z = Zt_vec(f, P)
    return f, a, z, 20 * np.log10(np.maximum(a * z, 1e-30) / 1e-6)

def worst_margin(P, ipk, duty, tr=TR, tf=TF, lim=cispr_B):
    f, a, z, dbv = spectrum(P, ipk, duty, tr, tf)
    m = lim(f) - dbv
    j = int(np.argmin(m))
    return m[j], f[j], dbv[j], lim(f)[j]

lines = []
def log(s=''):
    print(s); lines.append(s)

# ================================================================ RUN ============
log('=' * 78)
log('TPS62933 DM conducted-EMI  model v3  (REAL-GEOMETRY trace L, fsw=%.1f kHz)' % (FSW / 1e3))
log('=' * 78)
log('')
log('[geometry] board 26.50 x 50.50 mm; 2-layer; FR4 core 1.510 mm, er=4.5;')
log('           bottom = solid GND plane; input copper width 3.00 mm (measured).')
log('           L\' = %.4f nH/mm for a 3 mm trace over a 1.51 mm plane.' % G['Lpm'])
log('           lhot=%.2f  ltr_capA=%.2f  ltr_loop=%.2f  lf1=%.2f nH  (v2 pess: 20/8/30/10)'
    % (G['lhot'], G['ltr_capA'], G['ltr_loop'], G['lf1']))

Pv1 = P_v1repro()
log('')
log('[validation] v1-repro |Zt|@805kHz = %.3g ohm  (v2 published 6.32e-06)'
    % Zt(805e3, Pv1))
log('             v1-repro |Zt|@30MHz  = %.3g ohm  (v2 published 1.97e-05)'
    % Zt(30e6, Pv1))

# --- three-case nominal run ---
Pgeo = P_geo_nominal()
log('')
log('--- NOMINAL (real geometry trace L), pi filter, 50ohm LISN, tr=tf=10ns ---')
log('%-7s %-6s %-5s %9s %11s %10s %9s' % ('case', 'Iload', 'D', 'I_src1[A]', '|Zt|1[ohm]', 'V1[dBuV]', 'marg1[dB]'))
nom = {}
for nm, il, ipk, duty, col in CASES:
    f, a, z, dbv = spectrum(Pgeo, ipk, duty)
    m = cispr_B(f) - dbv
    nom[nm] = dict(f=f, dbv=dbv, m=m, col=col, il=il)
    log('%-7s %-6.2f %-5.2f %9.3f %11.3g %10.1f %9.1f' % (nm, il, duty, a[0], z[0], dbv[0], m[0]))
log('')
log('Worst-case margin vs CISPR 32 Class B (min over 150k-30M), real geometry:')
for nm, il, ipk, duty, col in CASES:
    r = nom[nm]; j = int(np.argmin(r['m']))
    log('  %-7s margin_min = %6.1f dB @ %.2f MHz (V=%.1f dBuV, lim=%.1f)  viol=%d'
        % (nm, r['m'][j], r['f'][j] / 1e6, r['dbv'][j], cispr_B(r['f'])[j], int((r['m'] < 0).sum())))

# --- envelope full load: geometric / optimistic / pessimistic ---
log('')
log('--- ENVELOPE, full load (Ipk=3A, D=0.5) ---')
env = {}
for tag, P, tr in [('geo-optimistic', P_geo_optimistic(), 10e-9),
                   ('geo-nominal',    Pgeo, 10e-9),
                   ('v1-repro',       Pv1, 15e-9),
                   ('pessimistic',    P_pessimistic(), 20e-9)]:
    m, fj, vj, limj = worst_margin(P, 3.0, 0.50, tr, tr)
    env[tag] = (m, fj, vj)
    log('  %-15s worst margin = %6.1f dB @ %6.2f MHz (V=%6.1f dBuV, lim=%4.1f)  |Zt|@805k=%.3g ohm'
        % (tag, m, fj / 1e6, vj, limj, Zt(805e3, P)))
opt_m = max(v[0] for v in env.values()); pes_m = min(v[0] for v in env.values())
log('  => ENVELOPE: %+.1f dB .. %+.1f dB  (spread %.1f dB)' % (opt_m, pes_m, opt_m - pes_m))

# --- what changed vs v2 (attribution of the trace-L fix) ---
log('')
log('--- effect of the trace-L fix (full load, all else = v2 nominal) ---')
def clone(d): return dict(d)
Pv2nom = dict(Pgeo); Pv2nom.update(lf1=5e-9, ltr_loop=15e-9, lhot=8e-9, ltr_capA=3e-9)
m0, f0, *_ = worst_margin(Pv2nom, 3.0, 0.5)
m1, f1, *_ = worst_margin(Pgeo, 3.0, 0.5)
log('  v2-nominal trace L (1 nH/mm-ish): worst margin %+.1f dB @ %.2f MHz' % (m0, f0 / 1e6))
log('  v3 geometric trace L            : worst margin %+.1f dB @ %.2f MHz' % (m1, f1 / 1e6))
log('  => improvement from using real geometry: %+.1f dB' % (m1 - m0))
for lab, key in [('hot-loop', 'lhot'), ('30uF bank', 'ltr_capA'),
                 ('L1 branch', 'ltr_loop'), ('F1/input', 'lf1')]:
    Pt = clone(Pgeo); Pt[key] = 20e-9 if key != 'lf1' else 10e-9     # v2 pessimistic value
    mm_, ff, *_ = worst_margin(Pt, 3.0, 0.5)
    log('  set %-10s -> v2-pessimistic %4.0f nH : margin %+.1f dB (delta %+.1f)'
        % (lab, Pt[key] * 1e9, mm_, mm_ - m1))

# --- |Zt| checkpoints ---
fg = np.logspace(np.log10(150e3), np.log10(30e6), 3000)
ztg = Zt_vec(fg, Pgeo)
log('')
log('|Zt| geo-nominal: @150kHz=%.3g  @805kHz=%.3g  @5MHz=%.3g  @30MHz=%.3g ohm'
    % (np.interp(150e3, fg, ztg), np.interp(805e3, fg, ztg),
       np.interp(5e6, fg, ztg), np.interp(30e6, fg, ztg)))

FULL_GEO = nom['full']['m'].min()
log('')
log('*** ANSWER: with 3.00 mm input copper (real geometry), full-load DM worst')
log('    Class-B margin = %+.1f dB (geo-nominal); envelope %+.1f .. %+.1f dB. ***'
    % (FULL_GEO, opt_m, pes_m))

# ================================================================ FIGS ===========
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.6))
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
for ax in (ax1, ax2):
    ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B')
    ax.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.2, label='CISPR 32 Class A')
    ax.set_xlim(0.15, 30); ax.set_ylim(-30, 100)
    ax.grid(True, which='both', alpha=0.3); ax.set_xlabel('Frequency [MHz]')
for nm, il, ipk, duty, col in CASES:
    r = nom[nm]
    ax1.semilogx(r['f'] / 1e6, r['dbv'], '-o', ms=3.5, color=col,
                 label='%s Iload=%.2fA  min-marg=%+.1f dB' % (nm, il, r['m'].min()))
ax1.set_ylabel('DM noise on LISN [dB$\\mu$V]')
ax1.set_title('(a) three load cases, REAL-GEOMETRY trace L')
ax1.legend(loc='upper right', fontsize=8.5)
_, _, _, dbvi = spectrum(P_geo_optimistic(), 3.0, 0.5)
ax2.semilogx(nom['full']['f'] / 1e6, dbvi, '-o', ms=3, color='#2ca02c', label='geo-optimistic (L/2)')
ax2.semilogx(nom['full']['f'] / 1e6, nom['full']['dbv'], '-s', ms=3, color='#1f77b4', label='geo-nominal (measured geometry)')
_, _, _, dbvp = spectrum(P_pessimistic(), 3.0, 0.5, 20e-9, 20e-9)
ax2.semilogx(nom['full']['f'] / 1e6, dbvp, '-^', ms=3, color='#d62728', label='pessimistic (1 nH/mm, 100$\\Omega$)')
ax2.set_ylabel('DM noise on LISN [dB$\\mu$V]')
ax2.set_title('(b) full load: real-geometry vs pessimistic envelope')
ax2.legend(loc='upper right', fontsize=8.5)
fig.suptitle('TPS62933 DM conducted emission v3 (real-geometry trace L) vs CISPR 32', fontsize=12)
fig.tight_layout()
fig.savefig(OUT + 'fig_emi_cispr_v3.png', dpi=140)
log('')
log('Wrote fig_emi_cispr_v3.png')

with open(OUT + 'emi_numbers_v3.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
log('wrote emi_numbers_v3.txt')

# --- REPORT ---
rep = []
rep.append('# TPS62933 DM conducted-EMI — v3 (real-geometry trace inductance)\n')
rep.append('## What changed vs v2\n')
rep.append('v2 used the generic **1 nH/mm** rule for PCB trace inductance in its pessimistic '
           'envelope (hot loop 20 nH, 30 µF-bank trace 8 nH, L1 branch 30 nH), giving only '
           '**+0.7 dB** worst-case margin at full load — flagged as over-pessimistic because the '
           'input copper is 3 mm wide.\n')
rep.append('v3 replaces those numbers with closed-form estimates from the **actual board geometry** '
           'parsed out of `epro/pourSim.epru`.\n')
rep.append('## Geometry (measured, not assumed)\n')
rep.append('| quantity | value | source |\n|---|---|---|\n')
rep.append('| board outline | 26.50 × 50.50 mm | `pourSim.epru` POLY layer 11 |\n')
rep.append('| stack-up | 2-layer, FR4 core 59.449 mil = **1.510 mm**, εr=4.5 | LAYER_PHYS |\n')
rep.append('| bottom copper | **solid GND plane** over whole board | POUR6 (layer 2, net GND) |\n')
rep.append('| input copper width | **3.00 mm** (24VIN 118 mil, $1N22 119 mil, PPHV 118 mil) | POUR extents |\n')
rep.append('| trace-over-plane L′ | **%.4f nH/mm** @ 3 mm width (microstrip) | Hammerstad-Jensen |\n' % G['Lpm'])
rep.append('\nFormulas (`traceL.py`): microstrip Z₀ (Hammerstad-Jensen) + TEM relation '
           'L′ = Z₀·√εeff/c; cross-checked against the Rosa strip-over-plane partial inductance.\n')
rep.append('## Geometric trace inductances\n')
rep.append('| knob | geometric | v2 nominal | v2 pessimistic |\n|---|---|---|---|\n')
rep.append('| lhot (hot loop) | **%.2f nH** | 8 | 20 |\n' % G['lhot'])
rep.append('| ltr_capA (30 µF bank) | **%.2f nH** | 3 | 8 |\n' % G['ltr_capA'])
rep.append('| ltr_loop (L1 branch) | **%.2f nH** | 15 | 30 |\n' % G['ltr_loop'])
rep.append('| lf1 (F1/input) | **%.2f nH** | 5 | 10 |\n' % G['lf1'])
rep.append('\n## Results (CISPR 32 Class B, DM, 50 Ω LISN, tr=tf=10 ns)\n')
rep.append('| case | V1@fsw [dBµV] | worst margin [dB] |\n|---|---|---|\n')
for nm, il, ipk, duty, col in CASES:
    r = nom[nm]; j = int(np.argmin(r['m']))
    rep.append('| %s | %.1f | **%+.1f** @ %.2f MHz |\n' % (nm, r['dbv'][0], r['m'][j], r['f'][j] / 1e6))
rep.append('\n### Envelope, full load (Ipk = 3 A)\n')
for tag, P, tr in [('geo-optimistic (½ L)', P_geo_optimistic(), 10e-9),
                   ('geo-nominal (measured)', Pgeo, 10e-9),
                   ('pessimistic (1 nH/mm, 100 Ω)', P_pessimistic(), 20e-9)]:
    m, fj, vj, limj = worst_margin(P, 3.0, 0.50, tr, tr)
    rep.append('- **%s**: %+.1f dB @ %.2f MHz\n' % (tag, m, fj / 1e6))
rep.append('\n## Answer\n')
rep.append('**With 3.00 mm-wide input copper, the full-load DM worst-case Class-B margin is '
           '≈ %+.1f dB (real geometry).** Across optimistic..pessimistic the envelope is '
           '%+.1f .. %+.1f dB. The v2 "+0.7 dB" figure was an artefact of the 1 nH/mm trace '
           'assumption; using the measured 3 mm / 1.51 mm microstrip geometry '
           '(L′ ≈ %.2f nH/mm) raises the margin by ~%+.1f dB.\n'
           % (FULL_GEO, opt_m, pes_m, G['Lpm'], m1 - m0))
rep.append('## Caveats\n')
rep.append('- Trace lengths use component-anchor coordinates × 1.3 routing factor; the hot-loop '
           'and bank-stub are Euclidean (no factor). Values are ±30 % order-of-magnitude estimates, '
           'not extracted parasitics.\n')
rep.append('- Assumes the return current runs in the continuous bottom GND plane (microstrip). '
           'If the plane were broken under the input path, L′ would rise toward the 1 nH/mm bound.\n')
rep.append('- Component-to-model-node assignment (30 µF bank = C17/C62/C63; local caps = C72/73/74) '
           'is inferred from the v2 model structure + proximity; the schematic netlist was not re-traced here.\n')
rep.append('- ESR/ESL remain datasheet/typical values (unchanged from v2).\n')
with open(OUT + 'REPORT_emi_v3.md', 'w') as fh:
    fh.write(''.join(rep))
log('wrote REPORT_emi_v3.md')
log('done.')
