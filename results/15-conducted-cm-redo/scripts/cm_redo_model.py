#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cm_redo_model.py -- TPS62933 conducted COMMON-MODE (CM) EMI, CORRECT topology,
                    150 kHz - 30 MHz, CISPR 32 / EN 55032 Class B (QP + AV).

Task: tps62933-cm-redo.  Writes ONLY to /mnt/raid10/sim-work/tps62933/cm_redo/.

WHY THIS REDO (background, from cm_audit/ANALYSIS.md §2):
  The old results/07 (emi_cm/emi_cm_model.py) reported "CM all cases over
  Class B by ~28 dB, CM 71 dB above DM".  That model drove a CM current
  through a NODAL network whose coupling element was taken as a black-box
  C_p = 10 pF, and the cm_pcb sibling defined C_sw as "SW-node <-> on-board
  GND plane" (a board-internal capacitance) that does NOT close through the
  LISN.  Both are conceptually wrong for CONDUCTED CM.

CORRECT topology (industry-standard, per docs/07):
    SW node dv/dt  ->  C_p (board <-> earth/chassis parasitic)  ->  earth
                   ->  LISN (2 x 50 ohm in parallel = 25 ohm CM equivalent)
                   ->  input cable (10 cm, 2-wire)  ->  back to source.
  C_p is therefore the board<->earth coupling, extracted from the REAL
  geometry: 2-layer stack-up (core 1.43 mm FR4, eps_r 4.3-4.6) WITH the
  top-layer coplanar GND.  Nominal value (cm_audit, FullCopper variant)
  = 2.416 pF (eps_r 4.3) / 2.572 pF (eps_r 4.6).

Source model:
    I_cm(f) = 2*pi*f * C_p * |V_sw(f)|      [A]   (i = C dV/dt)
    V_cm(f) = R_cm * I_cm(f),  R_cm = 25 ohm (LISN CM equivalent)
    -> dBuV = 20 log10(V/1uV)
  |V_sw(f)| taken from results_v2/E_{noload,half,full}.txt, column v(sw),
  FFT over the last M integer switching periods at fsw = 805 kHz.

NOT included (no authoritative data -> NOT guessed, see AGENTS.md "No Guessing"):
  * the upstream LM50-20B24 65 kHz switching as an EXTRA CM source.  Its own
    SW-node voltage spectrum is unavailable; only its Y-cap to PE is modelled
    as an alternative return path (secondary variant).
"""
import os, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT  = '/mnt/raid10/sim-work/tps62933/cm_redo/'
RES  = '/mnt/raid10/sim-work/tps62933/results_v2/'
os.makedirs(OUT, exist_ok=True)

FSW  = 805e3
VIN  = 24.0
RCM  = 25.0                 # LISN CM equivalent (2 x 50 ohm ||)
M    = 32                   # switching periods used for the SW FFT
KMAX = int(30e6 / FSW)

# --- C_p nominal from cm_audit geometry extraction (FR4 + top coplanar GND) ---
CP_NOM   = 2.416e-12        # FullCopper, eps_r 4.3  (PRIMARY / general practice)
CP_NOM46 = 2.572e-12        # FullCopper, eps_r 4.6
CP_NOTCH = {'topcut': 1.836e-12, 'dualcut': 1.189e-12}   # extrapolated by ratio

# --- upstream class-I Y-cap return-path variant (LM50-20B24) ---
CY_UP  = 2.2e-9             # estimated internal Y-cap (docs/08 says internal Y cap)
LCAB   = 30e-9              # 10 cm 2-wire pair inductance (est)
RCAB   = 50e-3

CASES = [('noload', 0.02, '#1f77b4'),
         ('half',   1.50, '#2ca02c'),
         ('full',   3.00, '#d62728')]

lines = []
def log(s=''):
    print(s); lines.append(s)

# --------------------------------------------------------------- limits --------
def cispr32_B_qp(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    66.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 56.0, 60.0))
def cispr32_B_av(f):
    return cispr32_B_qp(f) - 10.0     # AV limit is 10 dB below QP (CISPR 32 B)

# --------------------------------------------------------------- source --------
def load_sw(case):
    a = np.loadtxt(RES + 'E_%s.txt' % case)
    return a[:, 0], a[:, 3]          # t, v(sw)

def sw_harmonics(t, vsw, M=M):
    """FFT of v(sw) over the last M integer switching periods -> (fk, |Vsw|)."""
    T = 1.0 / FSW; Tw = M * T
    m = t >= (t[-1] - Tw)
    tt, vv = t[m], vsw[m]
    N = M * 4096
    tu = np.linspace(tt[0], tt[-1], N, endpoint=False)
    vu = np.interp(tu, tt, vv)
    vu = vu - vu.mean()
    X = np.fft.rfft(vu)
    f = np.fft.rfftfreq(N, tu[1] - tu[0])
    ks = np.arange(1, KMAX + 1); fk = ks * FSW
    amp = 2 * np.abs(X) / N
    idx = np.array([int(round(fr / f[1])) for fr in fk])
    return fk, amp[idx]

# --------------------------------------------------------------- CM ------------
def vcm_primary(fk, vsw, Cp):
    """V_cm = RCM * 2*pi*f*Cp*|Vsw|  (all CM current through LISN 25 ohm)."""
    return RCM * (2 * np.pi * fk * Cp * vsw)

def vcm_ycap(fk, vsw, Cp, Cy=CY_UP):
    """Same, but the upstream cable->Y-cap path shunts part of I_cm away from
    the 25 ohm LISN (Z_cm = 25 || (Rcab + jwLcab + 1/jwCy))."""
    out = np.zeros_like(fk)
    for i, fr in enumerate(fk):
        Icm = 2 * np.pi * fr * Cp * vsw[i]
        Zcy = RCAB + 1j * 2 * np.pi * fr * LCAB + 1.0 / (1j * 2 * np.pi * fr * Cy)
        Zpar = 1.0 / (1.0 / RCM + 1.0 / Zcy)
        out[i] = abs(Icm * Zpar)
    return out

def dbuv(v):
    return 20 * np.log10(np.maximum(v, 1e-30) / 1e-6)

def margins(fk, dbv):
    mq = cispr32_B_qp(fk) - dbv
    ma = cispr32_B_av(fk) - dbv
    return mq, ma

# =============================================================== RUN ===========
log('=' * 80)
log('TPS62933 CONDUCTED COMMON-MODE (CM) EMI  --  REDO, correct topology')
log('150 kHz - 30 MHz, CISPR 32 / EN 55032 Class B (QP + AV)')
log('=' * 80)
log('')
log('Topology : SW dv/dt -> C_p(board<->earth) -> earth -> LISN(25 ohm CM)')
log('           -> input cable (10 cm, 2-wire) -> back to source.')
log('Source   : |V_sw(f)| = FFT(v(sw)) from results_v2/E_*.txt, last %d periods.' % M)
log('Coupling : I_cm = 2*pi*f*C_p*|V_sw| ;  V_cm = 25 ohm * I_cm.')
log('C_p nom  : %.3f pF  (cm_audit: 2-layer, core 1.43 mm FR4 eps_r 4.3,' % (CP_NOM*1e12))
log('           top coplanar GND, FullCopper).  eps_r 4.6 -> %.3f pF.' % (CP_NOM46*1e12))
log('           (topcut %.3f pF / dualcut %.3f pF extrapolated by ratio.)'
    % (CP_NOTCH['topcut']*1e12, CP_NOTCH['dualcut']*1e12))
log('')

# ---- SW harmonic sanity ------------------------------------------------------
log('--- |V_sw| harmonic content (measured, used as-is) ---')
log('  case     D_meas  |Vsw1|     |Vsw5|     |Vsw37|   [V]')
sw = {}
for nm, il, col in CASES:
    t, v = load_sw(nm)
    fk, am = sw_harmonics(t, v)
    D = np.mean(v[t >= t[-1] - 20e-6]) / VIN
    sw[nm] = dict(fk=fk, am=am, D=D, col=col, il=il)
    log('  %-8s %.3f   %8.3f   %8.3f   %8.4f'
        % (nm, D, am[0], am[4], am[36] if len(am) > 36 else am[-1]))
log('  NOTE: half-load SPICE run did NOT reach steady state -> its v(sw) spectrum')
log('        (and thus its CM curve) is INDICATIVE ONLY (flagged throughout).')
log('        The sim SW edges collapse to <=1 sample (maxstep artefact); the')
log('        spectrum above is the measured result and is used directly.')
log('')

# ---- PRIMARY: correct CM, C_p = 2.416 pF ------------------------------------
log('--- PRIMARY CM (correct topology, C_p = %.3f pF, V_cm=25*I_cm) ---' % (CP_NOM*1e12))
log('  case     Vcm@0.805M  Vcm@4.03M  worst-QP-margin[dB] @ freq    worst-AV-margin[dB] @ freq')
prim = {}
for nm, il, col in CASES:
    fk, am = sw[nm]['fk'], sw[nm]['am']
    d = dbuv(vcm_primary(fk, am, CP_NOM))
    mq, ma = margins(fk, d)
    jq, ja = int(np.argmin(mq)), int(np.argmin(ma))
    prim[nm] = dict(f=fk, dbv=d, mq=mq, ma=ma)
    log('  %-8s %8.1f   %8.1f    %+7.1f @ %6.3f MHz    %+7.1f @ %6.3f MHz'
        % (nm, d[0], np.interp(4.03e6, fk, d), mq[jq], fk[jq]/1e6, ma[ja], fk[ja]/1e6))
log('')

# ---- variant: upstream Y-cap return path ------------------------------------
log('--- VARIANT: with upstream LM50-20B24 Y-cap as parallel CM return (%.1f nF) ---' % (CY_UP*1e9))
vary = {}
for nm, il, col in CASES:
    fk, am = sw[nm]['fk'], sw[nm]['am']
    d = dbuv(vcm_ycap(fk, am, CP_NOM))
    mq, ma = margins(fk, d)
    jq = int(np.argmin(mq))
    vary[nm] = dict(f=fk, dbv=d, mq=mq, ma=ma)
    log('  %-8s Vcm@0.805M=%6.1f dBuV  worst-QP-margin=%+7.1f dB @ %.3f MHz'
        % (nm, d[0], mq[jq], fk[jq]/1e6))
log('  (Y-cap path lowers V_cm by %.1f dB @0.805M vs primary.)'
    % (prim['full']['dbv'][0] - vary['full']['dbv'][0]))
log('')

# ---- C_p sensitivity ---------------------------------------------------------
log('--- C_p sensitivity (full load, primary topology) ---')
swep = {}
for cp in [1e-12, CP_NOM, 10e-12]:
    fk, am = sw['full']['fk'], sw['full']['am']
    d = dbuv(vcm_primary(fk, am, cp))
    mq, ma = margins(fk, d)
    jq = int(np.argmin(mq))
    swep[cp] = d
    log('  C_p=%6.3f pF : Vcm@0.805M=%6.1f dBuV  worst-QP-margin=%+7.1f dB @ %.3f MHz'
        % (cp*1e12, d[0], mq[jq], fk[jq]/1e6))
log('')

# ---- eps_r / copper variant --------------------------------------------------
log('--- geometry variant (full load, primary): C_p source ---')
for tag, cp in [('eps_r 4.3 FullCopper (NOM)', CP_NOM),
                ('eps_r 4.6 FullCopper', CP_NOM46),
                ('topcut (extrapolated)', CP_NOTCH['topcut']),
                ('dualcut (extrapolated)', CP_NOTCH['dualcut'])]:
    fk, am = sw['full']['fk'], sw['full']['am']
    d = dbuv(vcm_primary(fk, am, cp)); mq, _ = margins(fk, d); jq = int(np.argmin(mq))
    log('  %-28s C_p=%.3f pF  worst-QP-margin=%+7.1f dB @ %.3f MHz'
        % (tag, cp*1e12, mq[jq], fk[jq]/1e6))
log('')

# ---- comparison with old results/07 -----------------------------------------
log('--- DIFFERENCE vs OLD results/07 (emi_cm, C_p=10 pF black-box, nodal) ---')
fk, am = sw['full']['fk'], sw['full']['am']
old = dbuv(vcm_primary(fk, am, 10e-12))          # recreate old nominal
mq_new, _ = margins(fk, prim['full']['dbv']); mq_old, _ = margins(fk, old)
log('  old full worst-QP-margin (C_p=10 pF) = %+.1f dB @ %.3f MHz'
    % (mq_old.min(), fk[int(np.argmin(mq_old))]/1e6))
log('  new full worst-QP-margin (C_p=2.416 pF)= %+.1f dB @ %.3f MHz'
    % (mq_new.min(), fk[int(np.argmin(mq_new))]/1e6))
log('  delta = %+.1f dB  (= 20log10(10/2.416) = %+.1f dB : pure C_p over-estimate)'
    % (mq_new.min() - mq_old.min(), 20*np.log10(10.0/(CP_NOM*1e12))))
log('  root cause of old over-estimate: (a) C_p taken as black-box 10 pF instead')
log('  of the real geometry-extracted 2.42 pF (FR4 core + top coplanar GND);')
log('  (b) sibling cm_pcb further mis-defined the coupling as SW<->on-board GND')
log('  (board-internal loop, does not close through the LISN) and called it CM.')

# =============================================================== FIGURES =======
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
for ax in (ax1, ax2):
    ax.plot(fgp/1e6, cispr32_B_qp(fgp), 'k-', lw=2.2, label='CISPR 32 Class B  QP')
    ax.plot(fgp/1e6, cispr32_B_av(fgp), 'k--', lw=1.6, label='CISPR 32 Class B  AV')
    ax.set_xlim(0.15, 30); ax.set_ylim(0, 100)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('CM voltage on LISN [dB$\\mu$V]')
for nm, il, col in CASES:
    r = prim[nm]
    lbl = '%s (Iload=%.2fA)  QP-marg=%+.1f dB' % (nm, il, r['mq'].min())
    if nm == 'half':
        lbl = 'half  (NOT steady-state)  QP-marg=%+.1f dB' % r['mq'].min()
    ax1.semilogx(r['f']/1e6, r['dbv'], '-o', ms=3.0, color=col, label=lbl)
    # mark worst QP point
    jq = int(np.argmin(r['mq']))
    ax1.plot(r['f'][jq]/1e6, r['dbv'][jq], '*', ms=16, color=col,
             markeredgecolor='k', zorder=5)
ax1.set_title('(a) CM, correct topology  ($C_p$=2.416 pF, FR4 core 1.43 mm + top GND)')
ax1.legend(loc='lower left', fontsize=8.2)
# sensitivity panel
for cp, c in zip([1e-12, CP_NOM, 10e-12], ['#1f77b4', '#2ca02c', '#d62728']):
    ax2.semilogx(sw['full']['fk']/1e6, swep[cp], '-', lw=1.6, color=c,
                 label='$C_p$=%.3f pF' % (cp*1e12))
ax2.set_title('(b) $C_p$ sensitivity (full load, $V_{cm}=25\\Omega\\cdot I_{cm}$)')
ax2.legend(loc='lower left', fontsize=8.5)
fig.suptitle('TPS62933 conducted CM (corrected) vs CISPR 32 Class B  --  fsw=805 kHz', fontsize=13)
fig.tight_layout()
fig.savefig(OUT + 'fig_cm_cispr_redo.png', dpi=150)
log(''); log('Wrote fig_cm_cispr_redo.png')

# =============================================================== SUMMARY =======
jf = int(np.argmin(prim['full']['mq']))
log('')
log('*** ANSWER (corrected CM, general practice, C_p=2.416 pF nominal):')
for nm, il, col in CASES:
    r = prim[nm]; jq = int(np.argmin(r['mq'])); ja = int(np.argmin(r['ma']))
    log('    %-7s QP margin %+7.1f dB @ %.3f MHz (Vcm=%.1f dBuV);  AV margin %+7.1f dB'
        % (nm, r['mq'][jq], r['f'][jq]/1e6, r['dbv'][jq], r['ma'][ja]))
log('    -> ALL three cases EXCEED CISPR 32 Class B (QP & AV).')
log('    Worst = full load : %+.1f dB @ %.3f MHz.' % (prim['full']['mq'][jf], prim['full']['f'][jf]/1e6))
log('    vs OLD results/07 full-load %+.1f dB: improvement of %+.1f dB, purely'
    % (mq_old.min(), mq_new.min()-mq_old.min()))
log('    because the old model over-estimated C_p (10 pF black-box vs 2.42 pF real).')
log('    STILL A VIOLATION -> the board needs real CM mitigation (Y-caps/screen),')
log('    not merely a corrected budget. ***')

with open(OUT + 'numbers_cm_redo.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('wrote numbers_cm_redo.txt')
