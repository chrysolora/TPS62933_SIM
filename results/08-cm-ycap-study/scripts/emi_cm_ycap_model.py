#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emi_cm_ycap_model.py -- TPS62933 conducted COMMON-MODE EMI model WITH a
                        board-level Y capacitor to PE (earth), 150 kHz - 30 MHz,
                        CISPR 32, three load cases.

DERIVED FROM: /mnt/raid10/sim-work/tps62933/emi_cm/emi_cm_model.py
UNIQUE DIFFERENCE vs the baseline (emi_cm):
  * ONE extra branch is added to the CM loop: a local Y capacitor
        node0 (board GND, CM injection) <-> earth/PE,   Z = 1/(2*pi*f*Cyl)
    Everything else is IDENTICAL (C_p=10 pF, 10 cm / 30 nH cable, upstream
    class-I Y-cap path 2.2 nF, LISN 25 ohm, DM model, measured SW spectrum).
  * A second optional element (Cyl == 0, Lchoke == 0 by default) is a series
    common-mode choke impedance inserted in the whole CM loop (comparison only).
With Cyl=0 and Lchoke=0 the model reduces EXACTLY to the baseline network.

Writes ONLY to /mnt/raid10/sim-work/tps62933/emi_cm_ycap/.

Prerequisite doc (authoritative): tps62933-sim/docs/07-input-conditions-and-assumptions.md
  - board side has NO PE, only +24V/GND two wires (10 cm, 0.75 mm^2)
  - plastic (unshielded) enclosure
  - C_p = board<->earth parasitic is the ONLY element closing the CM loop
  - CISPR 32 Class B, fsw 805 kHz, tr=tf=5 ns (assumed)

IMPORTANT (premise change): adding a board Y-cap to PE CREATES a PE path on the
board side.  This CHANGES premise A/C of docs/07 ("板侧无 PE") and is therefore
an ASSUMED IMPROVEMENT SCENARIO, not the current measured configuration.
"""
import os, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT  = '/mnt/raid10/sim-work/tps62933/emi_cm_ycap/'
RES  = '/mnt/raid10/sim-work/tps62933/results_v2/'
os.makedirs(OUT, exist_ok=True)

FSW  = 805e3
VIN  = 24.0
TR   = TF = 5e-9
KMAX = int(30e6 / FSW)

lines = []
def log(s=''):
    print(s); lines.append(s)

# ============================================================ helpers =========
def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L
def capZ(f, C, esr, esl): return esr + Zl(f, esl) + Zc(f, C)

def solve(n, branches, inj, Isrc=1.0):
    Y = np.zeros((n, n), complex); b = np.zeros(n, complex)
    for (i, j, Z) in branches:
        if abs(Z) < 1e-15: Z = 1e-12
        y = 1.0 / Z
        if i >= 0: Y[i, i] += y
        if j >= 0: Y[j, j] += y
        if i >= 0 and j >= 0: Y[i, j] -= y; Y[j, i] -= y
    b[inj] += Isrc
    return np.linalg.solve(Y, b)

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

# ============================================================ CM network ======
#  node 0 = EUT board (CM injection) ; node 1 = LISN/lines side ; node 2 = upstream
def P_cm(Cp, Cy=2.2e-9, Lcab=30e-9, Rcab=50e-3, Rlisn=25.0, ypath=True,
         Cyl=0.0, Lchoke=0.0):
    return dict(Cp=Cp, Cy=Cy, Lcab=Lcab, Rcab=Rcab, Rlisn=Rlisn, ypath=ypath,
                Cyl=Cyl, Lchoke=Lchoke)

def Zcm_branches(f, P):
    """Return (n_nodes, branch_list, node_of_LISN_voltage)."""
    b = []
    if P['Lchoke'] > 0:                        # CM choke in series with the loop
        b.append((0, 1, Zl(f, P['Lchoke'])))   # board -> line side
        lisn_node = 1
        n = 3
        b.append((1, -1, P['Rlisn']))
        if P['ypath']:
            b.append((1, 2, P['Rcab'] + Zl(f, P['Lcab'])))
            b.append((2, -1, Zc(f, P['Cy'])))
        if P['Cyl'] > 0:
            b.append((0, -1, Zc(f, P['Cyl'])))
    else:
        lisn_node = 0
        if P['ypath']:
            n = 2
            b.append((0, -1, P['Rlisn']))
            b.append((0, 1, P['Rcab'] + Zl(f, P['Lcab'])))
            b.append((1, -1, Zc(f, P['Cy'])))
        else:
            n = 1
            b.append((0, -1, P['Rlisn']))
        if P['Cyl'] > 0:
            b.append((0, -1, Zc(f, P['Cyl'])))
    return n, b, lisn_node

def Vcm_spectrum(fk, Vsw_k, P):
    Icm = 2 * np.pi * fk * P['Cp'] * Vsw_k
    out = np.zeros_like(fk)
    for i, fr in enumerate(fk):
        n, b, ln = Zcm_branches(fr, P)
        V = solve(n, b, 0, Icm[i])
        out[i] = abs(V[ln])
    return out

def Zcm_scalar(f, P):
    n, b, ln = Zcm_branches(f, P)
    V = solve(n, b, 0, 1.0)
    return abs(V[ln])

# ============================================================ load sim SW =====
def load_sw(case):
    a = np.loadtxt(RES + 'E_%s.txt' % case)
    return a[:, 0], a[:, 3]

def measured_sw_spec(t, vsw, M=32):
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
    return fk, amp[idx], vv.mean() / VIN

CASES = [('noload', 0.02, 0.22, 0.15, '#1f77b4'),
         ('half',   1.50, 1.50, 0.50, '#2ca02c'),
         ('full',   3.00, 3.00, 0.50, '#d62728')]

def dbuv(v): return 20 * np.log10(np.maximum(v, 1e-30) / 1e-6)

# ============================================================ RUN =============
log('=' * 80)
log('TPS62933 conducted CM model WITH board Y-cap to PE  (derived from emi_cm)')
log('150 kHz - 30 MHz, CISPR 32, three load cases')
log('=' * 80)
log('UNIQUE DIFFERENCE vs baseline: + one Y-cap branch  board-GND <-> PE/earth.')
log('')

# --- V_sw harmonics (measured, primary source, identical to baseline) ----------
sws = {}
for nm, il, ipk, dty, col in CASES:
    t, vsw = load_sw(nm)
    fk, am, dm = measured_sw_spec(t, vsw)
    sws[nm] = dict(fk=fk, am=am, dm=dm)

CP = 10e-12
Pbase = P_cm(CP, Cyl=0.0)          # baseline: no local Y-cap
P1n   = P_cm(CP, Cyl=1e-9)         # user-specified 1 nF / 2 kV

# --- sanity: reproduce baseline ------------------------------------------------
log('--- SANITY (must reproduce baseline emi_cm): Cyl=0, C_p=10 pF ---')
base = {}
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    V = Vcm_spectrum(fk, sws[nm]['am'], Pbase)
    d = dbuv(V); m = cispr_B(fk) - d; j = int(np.argmin(m))
    base[nm] = dict(f=fk, dbv=d, m=m)
    log('  %-6s worst-margin = %+7.2f dB @ %.3f MHz  (Vcm@805k=%.1f dBuV)'
        % (nm, m[j], fk[j] / 1e6, d[0]))
log('  baseline reference (REPORT_emi_cm): full-load CM worst = -28.0 dB @ 0.81 MHz')
log('')

# --- MAIN: 1 nF vs baseline ----------------------------------------------------
log('--- MAIN RESULT: C_y = 1 nF (2 kV, user-specified) ---')
one = {}
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    V = Vcm_spectrum(fk, sws[nm]['am'], P1n)
    d = dbuv(V); m = cispr_B(fk) - d; j = int(np.argmin(m))
    one[nm] = dict(f=fk, dbv=d, m=m)
    # per-freq attenuation vs baseline
    att = base[nm]['dbv'] - d            # positive = reduction (dB)
    i1 = 0                               # 0.805 MHz (k=1)
    i30 = int(np.argmin(np.abs(fk - 30e6)))   # ~29.8 MHz (k=37)
    log('  %-6s base-margin=%+6.2f dB  ->  1nF margin=%+6.2f dB   (delta %+5.2f dB) @ %.3f MHz'
        % (nm, base[nm]['m'].min(), m[j], m[j] - base[nm]['m'].min(), fk[j] / 1e6))
    log('         attenuation @0.805MHz = %5.2f dB ; @%.3fMHz = %5.2f dB'
        % (att[i1], fk[i30] / 1e6, att[i30]))
log('')

# --- Z_Cy(f) / 25 ohm explanatory table ----------------------------------------
log('--- Z_Cy(f)=1/(2*pi*f*C_y) vs 25 ohm (why low freq is hard) ---')
log('   C_y=1 nF :')
for fr in [0.15e6, 0.805e6, 2e6, 10e6, 30e6]:
    z = 1.0 / (2 * np.pi * fr * 1e-9)
    log('     f=%7.3f MHz : Z_Cy=%9.1f ohm   Z_Cy/25 = %8.1f' % (fr / 1e6, z, z / 25.0))
log('')

# --- 1 nF band summary ---------------------------------------------------------
log('--- 1 nF attenuation by band (vs baseline), per case ---')
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']; att = base[nm]['dbv'] - one[nm]['dbv']
    lo = (fk >= 0.15e6) & (fk <= 2e6)
    hi = (fk >= 10e6) & (fk <= 30e6)
    log('   %-6s  low 0.15-2 MHz: mean %5.2f dB  (max %5.2f)   |  high 10-30 MHz: mean %5.2f dB (max %5.2f)'
        % (nm, att[lo].mean(), att[lo].max(), att[hi].mean(), att[hi].max()))
log('')

# --- C_y SWEEP -----------------------------------------------------------------
CYL = [0.1e-9, 1e-9, 4.7e-9, 10e-9, 47e-9, 100e-9]
log('--- C_y sweep: worst Class-B margin per case ---')
log('   C_y[nF]  ' + '  '.join('%9s' % nm for nm, *_ in CASES) + '   le50Hz@230V[mA]')
sweep_margin = {}   # nm -> list of worst margins aligned with CYL
for nm, *_ in CASES: sweep_margin[nm] = []
for Cy in CYL:
    P = P_cm(CP, Cyl=Cy)
    row = []
    for nm, il, ipk, dty, col in CASES:
        fk = sws[nm]['fk']
        V = Vcm_spectrum(fk, sws[nm]['am'], P)
        d = dbuv(V); m = cispr_B(fk) - d
        wm = m.min()
        row.append(wm); sweep_margin[nm].append(wm)
    ileak = 2 * np.pi * 50 * Cy * 230.0 * 1e3   # mA
    log('   %8.1f   ' % (Cy * 1e9) + '   '.join('%+9.2f' % r for r in row) + '   %10.2f' % ileak)
log('   (margins in dB; negative = violation.  Baseline row = C_y=0.)')
log('')

# --- minimum C_y required for Class B (margin>=0) ------------------------------
log('--- minimum C_y to reach CISPR 32 Class B (margin >= 0), per case ---')
log('   (obtained by fine numeric scan; large C_y -> Z_Cy -> 0 -> shunt kills CM)')
Cgr = np.logspace(np.log10(0.05e-9), np.log10(3e-6), 260)
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    need = None
    for Cy in Cgr:
        P = P_cm(CP, Cyl=Cy)
        V = Vcm_spectrum(fk, sws[nm]['am'], P)
        if (cispr_B(fk) - dbuv(V)).min() >= 0.0:
            need = Cy; break
    if need:
        ileak = 2 * np.pi * 50 * need * 230.0 * 1e3
        log('   %-6s  C_y,min = %8.1f nF   (touch leakage @230V/50Hz = %7.1f mA)'
            % (nm, need * 1e9, ileak))
    else:
        log('   %-6s  no solution within scan range (<= 3 uF)' % nm)
log('')

# --- CM CHOKE comparison (first order: Norton->Thevenin on C_p) ----------------
log('--- comparison: series CM choke, first-order (C_p as Norton impedance) ---')
log('   NOTE: the baseline uses an IDEAL CM current source (Z_Cp = inf).  A series')
log('   choke only reduces the current if C_p is given its Norton impedance')
log('   Z_Cp = 1/(2*pi*f*C_p):  I_eff = I_cm * |Z_Cp/(Z_Cp + Z_choke)|  (the small')
log('   ~25 ohm loop impedance is neglected).  Z_choke = 2*pi*f*L (ideal).')
for L, tag in [(1e-3, '1 mH'), (4.7e-3, '4.7 mH'), (10e-3, '10 mH')]:
    row = []
    for nm, il, ipk, dty, col in CASES:
        fk = sws[nm]['fk']
        Zcp = 1.0 / (2 * np.pi * fk * CP)
        Zch = 2 * np.pi * fk * L
        Ieff = Zcp / np.abs(Zcp + 1j * Zch)
        V = Vcm_spectrum(fk, sws[nm]['am'], Pbase) * Ieff
        m = cispr_B(fk) - dbuv(V)
        row.append(m.min())
    log('   %-8s choke  ' % tag + '  '.join('%s=%+6.2f' % (nm, r) for (nm, *_2), r in zip(CASES, row)))
log('   choke effectiveness vs frequency (1 mH): reduction of I_cm in dB:')
for fr in [0.805e6, 2e6, 10e6, 30e6]:
    Zcp = 1.0 / (2 * np.pi * fr * CP); Zch = 2 * np.pi * fr * 1e-3
    log('      f=%6.2f MHz : Z_Cp=%9.1f ohm  Z_choke=%9.1f ohm  -> %+6.2f dB'
        % (fr / 1e6, Zcp, Zch, 20 * np.log10(Zcp / abs(Zcp + 1j * Zch))))
log('   NOTE: a CM choke is SERIES impedance, and it competes against Z_Cp.  At')
log('   0.805 MHz Z_Cp is high (~20 kohm) vs a 1 mH choke (~5 kohm) -> tiny effect;')
log('   it only bites at high frequency where Z_Cp has fallen.  Opposite trend to')
log('   the shunt Y-cap.  A real choke also self-resonates (parasitic C) so the')
log('   >10 MHz figures are optimistic.')
log('')

# --- pick "best C_y": smallest scanned value that makes ALL three cases pass ---
Cgr_wide = np.logspace(np.log10(0.05e-9), np.log10(3e-6), 600)
best_Cy = None
for Cy in Cgr_wide:
    P = P_cm(CP, Cyl=Cy)
    ok = True
    for nm, *_ in CASES:
        fk = sws[nm]['fk']
        if (cispr_B(fk) - dbuv(Vcm_spectrum(fk, sws[nm]['am'], P))).min() < 0: ok = False; break
    if ok:
        best_Cy = Cy; break
if best_Cy is None: best_Cy = Cgr[-1]
log('--- "best C_y" (smallest scanned that passes ALL three cases) = %.1f nF ---'
    % (best_Cy * 1e9))
bestP = P_cm(CP, Cyl=best_Cy)
bestv = {}
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    d = dbuv(Vcm_spectrum(fk, sws[nm]['am'], bestP))
    bestv[nm] = dict(f=fk, dbv=d, m=cispr_B(fk) - d)

# ============================================================ FIGURES =========
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)

# --- compare figure: baseline vs 1nF vs best C_y -------------------------------
fig, axs = plt.subplots(3, 1, figsize=(10.5, 12.5), sharex=True)
for k, (nm, il, ipk, dty, col) in enumerate(CASES):
    ax = axs[k]
    ax.semilogx(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B')
    ax.semilogx(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.0, label='CISPR 32 Class A')
    ax.semilogx(base[nm]['f'] / 1e6, base[nm]['dbv'], '-o', ms=3, color=col,
                label='baseline (no Y-cap)  worst=%+.1f dB' % base[nm]['m'].min())
    ax.semilogx(one[nm]['f'] / 1e6, one[nm]['dbv'], '-s', ms=3, color='#ff7f0e',
                label='C_y=1 nF  worst=%+.1f dB' % one[nm]['m'].min())
    ax.semilogx(bestv[nm]['f'] / 1e6, bestv[nm]['dbv'], '-^', ms=3, color='#9467bd',
                label='C_y=%.0f nF (opt)  worst=%+.1f dB' % (best_Cy * 1e9, bestv[nm]['m'].min()))
    ax.set_xlim(0.15, 30); ax.set_ylim(30, 110); ax.grid(True, which='both', alpha=0.3)
    ax.set_ylabel('CM voltage on LISN [dB$\\mu$V]')
    ax.set_title('case %s (Iload=%.2f A)' % (nm, il))
    ax.legend(loc='upper right', fontsize=8.3)
axs[-1].set_xlabel('Frequency [MHz]')
fig.suptitle('TPS62933 CM: baseline vs Y-cap to PE (C_p=10 pF, fsw=805 kHz)', fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig_emi_cm_ycap_compare.png', dpi=140)
log('Wrote fig_emi_cm_ycap_compare.png')

# --- sweep figure: worst margin vs C_y ----------------------------------------
fig2, ax = plt.subplots(figsize=(9.5, 6.5))
allC = np.array(sorted(set(CYL + [best_Cy, 300e-9])))
for nm, il, ipk, dty, col in CASES:
    fk = sws[nm]['fk']
    ms = []
    for Cy in allC:
        P = P_cm(CP, Cyl=Cy)
        ms.append((cispr_B(fk) - dbuv(Vcm_spectrum(fk, sws[nm]['am'], P))).min())
    ax.semilogx(allC * 1e9, ms, '-o', ms=4, color=col, label='%s (Iload=%.2f A)' % (nm, il))
ax.axhline(0, color='k', lw=1.8, ls='-', label='Class B limit (margin=0)')
ax.axvline(1.0, color='#ff7f0e', ls=':', lw=1.5, label='user C_y = 1 nF')
ax.set_xlabel('C_y to PE [nF]'); ax.set_ylabel('worst CISPR 32 Class-B margin [dB]')
ax.set_title('CM worst margin vs board Y-cap to PE (C_p=10 pF)')
ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=9, loc='lower right')
fig2.tight_layout(); fig2.savefig(OUT + 'fig_emi_cm_ycap_sweep.png', dpi=140)
log('Wrote fig_emi_cm_ycap_sweep.png')

with open(OUT + 'emi_numbers_ycap.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('wrote emi_numbers_ycap.txt')

log('')
log('DONE.')
