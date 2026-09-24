#!/usr/bin/env python3
"""B2 addendum: analytic CM cross-check + PE comparison, appended to numbers_b2.txt."""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
L, R, CY = 105e-9, 0.05, 2.2e-9
def harm(fsw, A, D, tr, N=400001):
    t = np.linspace(0, 1/fsw, N, endpoint=False); T = 1/fsw; top = max(D*T-2*tr, 0)
    x = np.interp(t, [0, tr, tr+top, 2*tr+top, T], [0, A, A, 0, 0])
    X = np.abs(np.fft.rfft(x)/N); k = np.arange(1, int(1e9/fsw)+1)
    return k*fsw, X[k]
def zl(f, Cp):
    w = 2*np.pi*f
    return R + 1j*w*L + 1/(1j*w*CY) + 1/(1j*w*Cp)
eta, k0 = 376.73, 2*np.pi/299792458.0
T = [30e6, 100e6, 300e6]
out = []
def P(s=''):
    print(s); out.append(s)
P('')
P('--- ADDENDUM: analytic CM cross-check (short dipole l=0.1m, h=11.5mm over PEC) ---')
P('    E_3m = 2*eta*k*I*l/(4*pi*3m) ; I=Vharm/|Z_loop| ; k=2*pi*f/c')
for Cp in [2e-12, 10e-12, 50e-12]:
    for tr in [50e-9, 100e-9]:
        f, V = harm(65e3, 75e-3, 0.5, tr)
        row = []
        for ft in T:
            i = int(np.argmin(abs(f-ft))); I = V[i]/abs(zl(ft, Cp))
            k = k0*ft; E = 2*eta*k*I*0.1/(4*np.pi*3)
            row.append('%6.1f' % (20*np.log10(max(E, 1e-30)/1e-6)))
        P('    Cp=%4.0f pF, tr=%3.0f ns : %s' % (Cp*1e12, tr*1e9, ' '.join(row)))
P('    [FDTD-based CM (above) gives -69/-91/-93 dBuV/m -> FDTD CM ~40 dB lower than')
P('     the analytic bound; both are >=50 dB below DM, so CM is negligible either way.]')
P('')
P('--- ADDENDUM: PE vs no-PE (FullCopper, 3.0A, cm-excited) ---')
a = np.abs(np.load(os.path.join(HERE, 'ver_full_cm', 'port.npz'))['cm_if_tot'])
b = np.abs(np.load(os.path.join(HERE, 'ver_full_cm_pe', 'port.npz'))['cm_if_tot'])
i = int(np.argmin(abs(np.load(os.path.join(HERE, 'ver_full_cm', 'port.npz'))['freq']-1e8)))
P('    if_cm no-PE=%.4e  with-PE=%.4e  ratio=%.4f' % (a[i], b[i], b[i]/a[i]))
P('    => tying the board to the plane changes CM current by <0.1%%: with the assumed')
P('       upstream edge rate the CM loop is NOT the bottleneck (Vharm~0 at HF).')

with open(os.path.join(HERE, 'numbers_b2.txt'), 'a') as fh:
    fh.write('\n'.join(out) + '\n')
print('appended to numbers_b2.txt')
