#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""post_conv.py -- radiation mesh-convergence check for the B2 'dump-copper' trend.

Compares E@3m(30 MHz) for FullCopper and DualCutout DM-only runs at three mesh levels:
  L0 = rad_b2 --fast          (resb 1.5, smooth 22, 318k cells)
  L1 = cm_audit resb 1.0/18   (425k cells)
  L2 = cm_audit resb 0.75/14  (588k cells)
Same source / same conditions (post_b2 normalization, Ipk=4.3A @3.0A load).
Writes ONLY to /mnt/raid10/sim-work/tps62933/cm_audit/.
"""
import os
import numpy as np

HERE = '/mnt/raid10/sim-work/tps62933/cm_audit'
RADB2 = '/mnt/raid10/sim-work/tps62933/rad_b2'

DM_FSW, DM_D, DM_TR = 805e3, 0.5, 5e-9

def harm(fmin, fmax, fsw, A, D, tr):
    N = 400001
    t = np.linspace(0, 1.0/fsw, N, endpoint=False)
    T = 1.0/fsw; top = max(D*T - 2*tr, 0.0)
    x = np.interp(t, [0, tr, tr+top, 2*tr+top, T], [0, A, A, 0, 0])
    X = np.abs(np.fft.rfft(x)/N)
    k = np.arange(1, int(fmax/fsw)+1); f = k*fsw
    m = (f >= fmin) & (f <= fmax)
    return f[m], X[k[m]]

def interp_cplx(H, fsrc, fdst):
    nf, a, b = H.shape
    Hf = H.reshape(nf, -1)
    out = np.empty((len(fdst), a*b), complex)
    for c in range(a*b):
        out[:, c] = np.interp(fdst, fsrc, Hf[:, c].real) + \
                    1j*np.interp(fdst, fsrc, Hf[:, c].imag)
    return out.reshape(len(fdst), a, b)

def dbmax(E):
    return 20*np.log10(np.maximum(np.abs(E).max(axis=(1, 2)), 1e-30)/1e-6)

def band(dB, FAP, t, half=2):
    i = int(np.argmin(np.abs(FAP - t)))
    return dB[max(0, i-half):i+half+1].max()

fh_dm, Ih_unit = harm(30e6, 1e9, DM_FSW, 1.0, DM_D, DM_TR)
FAP = fh_dm
IPK = 4.3                                   # 3.0 A load

LEVELS = [
    ('L0(fast 318k)', os.path.join(RADB2, 'ver_%s_dm')),
    ('L1(1.0/425k)',  os.path.join(HERE, 'ver_%s_dm_L1')),
    ('L2(0.75/588k)', os.path.join(HERE, 'ver_%s_dm_L2')),
]
T = [30e6, 100e6, 300e6]

lines = []
def log(s=''):
    print(s); lines.append(s)

log('=== B2 radiation mesh-convergence (DM only, 3.0A, Ipk=4.3A) ===')
log('E@3m [dBuV/m], band-max (+/-2 harmonics)')
E30 = {}
for v in ('full', 'dual'):
    log('')
    log('--- variant %s ---' % v)
    for lab, pat in LEVELS:
        p = pat % v
        fn = os.path.join(p, 'nf2ff.npz')
        if not os.path.exists(fn):
            log('  %-15s MISSING %s' % (lab, fn)); continue
        dd = np.load(fn); pd = np.load(os.path.join(p, 'port.npz'))
        Hdm = dd['E_norm'] / np.abs(pd['dm_if_tot'])[:, None, None]
        H = interp_cplx(Hdm, dd['freq'], FAP)
        Edm = H * (Ih_unit*IPK)[:, None, None] / 3.0
        dB = dbmax(Edm)
        row = [band(dB, FAP, t) for t in T]
        log('  %-15s  E@30M=%6.2f  E@100M=%6.2f  E@300M=%6.2f' % (lab, row[0], row[1], row[2]))
        E30.setdefault(v, []).append(row[0])

log('')
log('=== convergence vs trend ===')
for v in ('full', 'dual'):
    if len(E30.get(v, [])) >= 2:
        e = E30[v]
        log('  %-5s : levels E30 = %s   L0->L2 drift = %.2f dB' %
            (v, ' '.join('%.2f' % x for x in e), e[-1]-e[0]))
if all(len(E30.get(v, [])) >= 2 for v in ('full', 'dual')):
    for i in range(min(len(E30['full']), len(E30['dual']))):
        d = E30['dual'][i] - E30['full'][i]
        log('  level%d : Dual - Full = %+.2f dB' % (i, d))
    conv = max(abs(E30[v][-1]-E30[v][0]) for v in ('full', 'dual'))
    trend = abs(E30['dual'][-1] - E30['full'][-1])
    log('')
    log('  max self-convergence drift (L0->L2) = %.2f dB' % conv)
    log('  Dual-Full trend at finest level     = %.2f dB' % trend)
    log('  VERDICT: trend %s convergence drift -> %s' %
        ('>' if trend > conv else '<=',
         'resolvable' if trend > conv else 'NOT resolvable (within numerical noise)'))

open(os.path.join(HERE, 'conv_numbers.txt'), 'w').write('\n'.join(lines) + '\n')
print('wrote conv_numbers.txt')
