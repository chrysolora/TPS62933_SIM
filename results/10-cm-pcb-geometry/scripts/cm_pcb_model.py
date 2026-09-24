#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cm_pcb_model.py -- TPS62933 conducted CM, with the CM coupling capacitance
C_sw extracted from the REAL PCB copper geometry (Elmer StatElecSolver),
one value per pour variant (fullcu / topcut / dualcut).

Re-uses the CM circuit + CISPR-32 limits of emi_cm/emi_cm_model.py (read-only).
Writes ONLY to /mnt/raid10/sim-work/tps62933/cm_pcb/.

CM source: I_cm(f) = 2*pi*f*C_sw*|V_sw(f)|, V_sw from the MEASURED sim SW
waveform results_v2/E_{noload,half,full}.txt (col 4), integer-period FFT.
"""
import os, csv, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = '/mnt/raid10/sim-work/tps62933/cm_pcb/'
RES = '/mnt/raid10/sim-work/tps62933/results_v2/'
os.makedirs(OUT, exist_ok=True)

FSW = 805e3; VIN = 24.0; TR = TF = 5e-9
KMAX = int(30e6 / FSW)
VAR = ['fullcu', 'topcut', 'dualcut']
LBL = {'fullcu': 'FullCopper', 'topcut': 'TopCutout', 'dualcut': 'DualCutout'}
COL = {'fullcu': '#1f77b4', 'topcut': '#2ca02c', 'dualcut': '#d62728'}
CASES = [('noload', 0.02, '#1f77b4'), ('half', 1.50, '#2ca02c'), ('full', 3.00, '#d62728')]

lines = []
def log(s=''):
    print(s); lines.append(s)

def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L

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

def P_cm(Cp, Cy=2.2e-9, Lcab=30e-9, Rcab=50e-3, Rlisn=25.0, ypath=True):
    return dict(Cp=Cp, Cy=Cy, Lcab=Lcab, Rcab=Rcab, Rlisn=Rlisn, ypath=ypath)

def Zcm_branches(f, P):
    b = [(0, -1, P['Rlisn'])]
    if P['ypath']:
        b.append((0, 1, P['Rcab'] + Zl(f, P['Lcab'])))
        b.append((1, -1, Zc(f, P['Cy'])))
        return 2, b
    return 1, b

def Vcm_spectrum(fk, Vsw_k, P):
    Icm = 2 * np.pi * fk * P['Cp'] * Vsw_k
    out = np.zeros_like(fk)
    for i, fr in enumerate(fk):
        n, b = Zcm_branches(fr, P)
        V = solve(n, b, 0, Icm[i]); out[i] = abs(V[0])
    return out

def cispr_B(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    66.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 56.0, 60.0))

def load_sw(case):
    a = np.loadtxt(RES + 'E_%s.txt' % case)
    return a[:, 0], a[:, 3]

def measured_sw_spec(t, vsw, M=32):
    T = 1.0 / FSW; Tw = M * T
    m = t >= (t[-1] - Tw); tt, vv = t[m], vsw[m]
    N = M * 4096
    tu = np.linspace(tt[0], tt[-1], N, endpoint=False)
    vu = np.interp(tu, tt, vv); vu = vu - vu.mean()
    X = np.fft.rfft(vu); f = np.fft.rfftfreq(N, tu[1] - tu[0])
    ks = np.arange(1, KMAX + 1); fk = ks * FSW
    amp = 2 * np.abs(X) / N
    idx = np.array([int(round(fr / f[1])) for fr in fk])
    return fk, amp[idx], vv.mean() / VIN

def read_Csw():
    """C_sw [pF] at d=10cm nominal per variant from results.csv."""
    out = {}
    with open(OUT + 'results.csv') as fh:
        r = csv.DictReader(fh)
        for row in r:
            if row['C_pF'] in ('', None) or row['C_pF'] == 'FAIL':
                continue
            if abs(float(row['d_m']) - 0.1) < 1e-6 and abs(float(row['Lc']) - 0.0008) < 1e-12:
                out[row['variant']] = float(row['C_pF'])
    return out

if __name__ == '__main__':
    log('=' * 78)
    log('TPS62933 CM with PCB-derived C_sw (Elmer StatElec), CISPR-32 Class B')
    log('=' * 78)
    Csw = read_Csw()
    log('C_sw (d=10cm nominal) [pF]: ' + ', '.join('%s=%.4f' % (v, Csw.get(v, float('nan'))) for v in VAR))

    sws = {}
    for nm, il, col in CASES:
        t, vsw = load_sw(nm)
        fk, am, dm = measured_sw_spec(t, vsw)
        sws[nm] = dict(fk=fk, am=am, dm=dm)
        log('  SW case %-7s D_meas=%.3f  |Vsw1|=%.3f V  |Vsw5|=%.3f V'
            % (nm, dm, am[0], am[4]))

    # ---- per-variant CM ----
    res = {}
    for v in VAR:
        Cp = Csw[v] * 1e-12
        P = P_cm(Cp)
        d = {}
        for nm, il, col in CASES:
            fk = sws[nm]['fk']
            V = Vcm_spectrum(fk, sws[nm]['am'], P)
            dbv = 20 * np.log10(np.maximum(V, 1e-30) / 1e-6)
            m = cispr_B(fk) - dbv; j = int(np.argmin(m))
            d[nm] = dict(f=fk, dbv=dbv, m=m)
        res[v] = d
        log('%-8s C_sw=%.4f pF : ' % (v, Csw[v]) +
            '  '.join('%s worst-marg=%+.1f dB @%.2fMHz' % (nm, d[nm]['m'].min(),
                      d[nm]['f'][int(np.argmin(d[nm]['m']))] / 1e6) for nm, il, col in CASES))

    # ---- comparison to black-box 10 pF baseline ----
    log('')
    log('--- vs black-box baseline C_p=10 pF (emi_cm) ---')
    for nm, il, col in CASES:
        base = Vcm_spectrum(sws[nm]['fk'], sws[nm]['am'], P_cm(10e-12))
        dbb = 20 * np.log10(np.maximum(base, 1e-30) / 1e-6)
        log('  %-7s baseline worst-marg=%+.1f dB ; geometry gives (fullcu/topcut/dualcut) %s'
            % (nm, (cispr_B(sws[nm]['fk']) - dbb).min(),
               '/'.join('%+.1f' % res[v][nm]['m'].min() for v in VAR)))

    fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
    # ---- figure 1: three per-variant panels ----
    fig, axs = plt.subplots(1, 3, figsize=(19, 6.2), sharey=True)
    for ax, v in zip(axs, VAR):
        ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.0, label='CISPR 32 Class B (QP)')
        for nm, il, col in CASES:
            r = res[v][nm]
            ax.semilogx(r['f'] / 1e6, r['dbv'], '-o', ms=2.6, color=col,
                        label='%s  I=%.2fA  min-marg=%+.1f dB' % (nm, il, r['m'].min()))
        ax.set_xlim(0.15, 30); ax.set_ylim(0, 90); ax.grid(True, which='both', alpha=0.3)
        ax.set_xlabel('Frequency [MHz]')
        ax.set_title('%s   ($C_{sw}$=%.3f pF, d=10 cm)' % (LBL[v], Csw[v]))
        ax.legend(loc='upper right', fontsize=8)
    axs[0].set_ylabel('CM voltage on LISN [dB$\\mu$V]')
    fig.suptitle('TPS62933 CM vs CISPR 32 Class B -- C_sw from real PCB copper (Elmer)', fontsize=13)
    fig.tight_layout(); fig.savefig(OUT + 'fig_cm_per_variant.png', dpi=140)
    log('Wrote fig_cm_per_variant.png')

    # ---- figure 2: three-variant comparison, full load ----
    fig2, ax = plt.subplots(figsize=(10, 7))
    ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.2, label='CISPR 32 Class B (QP)')
    for v in VAR:
        r = res[v]['full']
        ax.semilogx(r['f'] / 1e6, r['dbv'], '-o', ms=3, color=COL[v],
                    label='%s (full, %.3f pF) min-marg=%+.1f dB' % (LBL[v], Csw[v], r['m'].min()))
    ax.set_xlim(0.15, 30); ax.set_ylim(0, 90); ax.grid(True, which='both', alpha=0.3)
    ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('CM voltage on LISN [dB$\\mu$V]')
    ax.set_title('Pour-variant comparison, full load (3.0A), d=10 cm')
    ax.legend(loc='upper right', fontsize=9)
    fig2.tight_layout(); fig2.savefig(OUT + 'fig_cm_variant_compare.png', dpi=140)
    log('Wrote fig_cm_variant_compare.png')

    with open(OUT + 'cm_pcb_numbers.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('wrote cm_pcb_numbers.txt')
