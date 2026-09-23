#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
traceL.py -- geometric PCB trace / loop inductance for the TPS62933 board.

Replaces the generic "1 nH/mm" rule of thumb with closed-form estimates that use
the REAL stack-up height and the REAL copper width parsed from the board file.

Stack-up (parsed from epro/pourSim.epru LAYER_PHYS):
    2-layer board; TOP copper / FR4 core / BOTTOM copper.
    core dielectric thickness  h = 59.449 mil = 1.510 mm ,  er = 4.5
    -> the input pours on TOP run over a CONTINUOUS BOTTOM GND plane (POUR6,
       net=GND, x[-60,1120] y[-2055,105] mil covers the whole board), i.e. a
       microstrip geometry with return on the plane 1.51 mm below.

Copper widths (parsed from POUR geometry, TOP layer):
    24VIN  pour x[210,328] -> 118 mil = 3.00 mm
    $1N22  pour x[402,521] -> 119 mil = 3.02 mm
    PPHV   pour x[589,707] -> 118 mil = 3.00 mm
    (user-stated "3 mm wide input copper" CONFIRMED by the CAD data)

FORMULAS
--------
Primary (microstrip, trace over a plane):
  Hammerstad-Jensen quasi-static microstrip impedance,
      u    = w/h
      eeff = (er+1)/2 + (er-1)/2 * (1 + 12/u)^-0.5                    [1]
      Z0   = (60/sqrt(eeff)) * ln( 6 + (2pi-6) exp[-(30.666/u)^0.7528] ) [1]
  A microstrip is a TEM line with phase velocity v = c/sqrt(eeff) and
      L' = Z0 / v = Z0 * sqrt(eeff) / c            [H per metre]        [2]
  and the *loop* inductance of a length l of that line is L = L' * l
  (the return current flows in the plane image ~ l below, so the loop length
   equals the trace length).                                          [2]

Cross-check (conductor over a ground plane, external partial inductance):
      L' = (mu0/2pi) * [ ln(2h/w) + 1/2 ]   (Rosa/Wheeler-type approx)  [3]

Reference for the generic rule of thumb it replaces: a trace whose return is
far away (>> w) sits near L' -> ~1 nH/mm.  For a 3 mm strip only 1.5 mm above a
plane, [2] gives ~0.36 nH/mm, i.e. the "1 nH/mm" envelope is ~3x pessimistic.
"""
import math

MU0 = 4e-7 * math.pi
C0 = 299792458.0
MIL = 25.4e-6


# --------------------------------------------------------------- formulas ----
def z0_microstrip(w, h, er):
    """Hammerstad-Jensen microstrip impedance [ohm]."""
    u = w / h
    eeff = (er + 1) / 2 + (er - 1) / 2 * (1 + 12 / u) ** -0.5
    a = (30.666 / u) ** 0.7528
    F = 6 + (2 * math.pi - 6) * math.exp(-a)
    z0 = 60 / math.sqrt(eeff) * math.log(F / u + math.sqrt(1 + (2 / u) ** 2))
    return z0, eeff


def L_per_m_microstrip(w, h, er):
    """Microstrip loop inductance per metre [H/m] (return = plane below)."""
    z0, eeff = z0_microstrip(w, h, er)
    return z0 * math.sqrt(eeff) / C0, z0, eeff


def L_per_m_over_plane(w, h):
    """External partial inductance of a strip over a plane [H/m] (cross-check)."""
    return MU0 / (2 * math.pi) * (math.log(2 * h / w) + 0.5)


def L_loop(l, Lp_m):
    """Loop inductance [H] of a trace of length l [m] with per-metre value."""
    return l * Lp_m


# --------------------------------------------------------------- geometry ----
# parsed from epro/pourSim.epru (see probe output).  mil units, screen y down.
BOARD = dict(x0=0.0, x1=1043.3, y0=-1958.2, y1=30.0)          # 26.50 x 50.50 mm
H_CORE = 59.449 * MIL            # 1.510 mm FR4 between TOP and BOTTOM copper
ER = 4.5
W_INPUT = 118.0 * MIL            # 3.00 mm  (measured pour width)

# component label anchor positions (mil), designator -> (x, y)
COMP = {
    'CN2': (196.85, -1800.71), 'F1': (210.0, -1465.0),
    'C17': (850.0, -1630.0), 'C62': (850.0, -1435.0), 'C63': (845.0, -1530.0),
    'C16': (865.0, -1120.0), 'C64': (865.0, -900.0),
    'C72': (425.0, -885.0), 'C73': (505.0, -885.0), 'C74': (600.0, -885.0),
    'L1': (645.0, -1250.0), 'U21': (325.0, -930.0), 'R85': (705.0, -1045.0),
}

# measured pour extents (mil) used for the input-path length
POUR_LEN = {'24VIN': (328.0 - 210.0, 1850.0 - 1485.0),
            '$1N22': (521.0 - 402.0, 1765.0 - 1485.0),
            'PPHV':  (707.0 - 589.0, 1766.0 - 1301.0)}

ROUTE_FACTOR = 1.3     # routed length / straight-line (90deg corners, pad entries)


def dist(a, b):
    """straight-line distance in mil."""
    return math.hypot(a[0] - b[0], a[1] - b[1])


def mm(v):  # mil -> mm of *length*
    return v * MIL * 1e3


def bank_centroid():
    xs = [COMP[k][0] for k in ('C17', 'C62', 'C63')]
    ys = [COMP[k][1] for k in ('C17', 'C62', 'C63')]
    return (sum(xs) / 3, sum(ys) / 3)


# --------------------------------------------------------------- segments ----
def measure():
    """Return a dict of named segment lengths in mm (routed estimate)."""
    Cc = bank_centroid()
    seg = {}
    # input path: connector -> fuse -> 30uF bank (node A->P)
    seg['conn_fuse_mil'] = dist(COMP['CN2'], COMP['F1'])
    seg['fuse_bank_mil'] = dist(COMP['F1'], Cc)
    # bank internal stub: worst cap to bank centroid  (-> ltr_capA)
    seg['bank_stub_mil'] = max(dist(COMP[k], Cc) for k in ('C17', 'C62', 'C63'))
    # L1 branch: bank -> L1 -> U21 VIN     (-> ltr_loop)
    seg['bank_L1_mil'] = dist(Cc, COMP['L1'])
    seg['L1_U21_mil'] = dist(COMP['L1'], COMP['U21'])
    # hot loop: U21 VIN pin -> nearest local decoupling cap, x2 for the loop
    near = min(('C72', 'C73', 'C74'), key=lambda k: dist(COMP['U21'], COMP[k]))
    seg['hot_oneway_mil'] = dist(COMP['U21'], COMP[near])
    seg['hot_loop_mil'] = 2 * dist(COMP['U21'], COMP[near])
    seg['hot_nearest_cap'] = near
    # input pour total path (24VIN + $1N22 + PPHV pour lengths, y-extents)
    seg['input_pour_mil'] = sum(POUR_LEN[k][1] for k in POUR_LEN)
    return seg


def main():
    out = []
    def log(s=''):
        print(s); out.append(s)

    log('=' * 74)
    log('traceL.py -- geometric trace / loop inductance (TPS62933 input stage)')
    log('=' * 74)
    log('stack-up: 2-layer, FR4 core h=%.3f mm, er=%.2f, bottom = solid GND plane' %
        (H_CORE * 1e3, ER))
    log('input copper width (measured from pours): %.2f mm' % (W_INPUT * 1e3))

    # --- per-mm inductance for candidate widths -----------------------------
    log('')
    log('inductance per mm of a TOP trace over the 1.51mm GND plane:')
    log('%8s %10s %10s %12s %12s' % ('w[mm]', 'Z0[ohm]', 'eeff', 'L\'[nH/mm]', 'overplane'))
    for wmm in (0.2, 0.5, 1.0, 2.0, 3.0):
        w = wmm * 1e-3
        Lp, z0, eeff = L_per_m_microstrip(w, H_CORE, ER)
        lo = L_per_m_over_plane(w, H_CORE)
        log('%8.2f %10.1f %10.3f %12.4f %12.4f' %
            (wmm, z0, eeff, Lp * 1e9 * 1e-3, lo * 1e9 * 1e-3))
    Lp3, z03, eeff3 = L_per_m_microstrip(W_INPUT, H_CORE, ER)
    Lpm = Lp3 * 1e9 * 1e-3          # nH/mm
    log('=> 3.00 mm copper: L\' = %.4f nH/mm  (generic 1 nH/mm is %.1fx higher)'
        % (Lpm, 1.0 / Lpm))

    # --- segment lengths + inductances -------------------------------------
    s = measure()
    log('')
    log('measured segments (mil -> mm; routed = straight x %.2f):' % ROUTE_FACTOR)

    def line(name, mil_len, fac=True):
        L = mm(mil_len) * (ROUTE_FACTOR if fac else 1.0)
        H = L * 1e-3 * Lp3 * 1e9       # nH
        return name, mm(mil_len), L, H

    segs = [
        line('lf1  input  F1->bank', s['fuse_bank_mil']),
        line('lf1  conn->F1', s['conn_fuse_mil']),
        line('ltr_capA bank cap stub', s['bank_stub_mil'], fac=False),
        line('ltr_loop bank->L1', s['bank_L1_mil']),
        line('ltr_loop L1->U21', s['L1_U21_mil']),
        line('lhot  U21->cap oneway', s['hot_oneway_mil'], fac=False),
        line('lhot  loop (x2)', s['hot_loop_mil'], fac=False),
    ]
    log('%-26s %9s %9s %10s' % ('segment', 'straight', 'routed', 'L_geo[nH]'))
    for nm, st, rt, H in segs:
        log('%-26s %8.2f  %8.2f  %9.2f' % (nm, st, rt, H))

    # --- assembled model parameters ----------------------------------------
    # map segments to the v2 model knobs
    def nH(mil_len, fac=True):
        return mm(mil_len) * (ROUTE_FACTOR if fac else 1.0) * 1e-3 * Lp3 * 1e9
    lf1 = nH(s['fuse_bank_mil'])
    ltr_loop = nH(s['bank_L1_mil']) + nH(s['L1_U21_mil'])
    ltr_capA = nH(s['bank_stub_mil'], fac=False)
    lhot = nH(s['hot_oneway_mil'], fac=False)
    log('')
    log('--- geometric values for the v2 model knobs (nominal/geo) ---')
    log('  lf1      (F1/input trace)      = %.2f nH' % lf1)
    log('  ltr_loop (L1-branch trace)     = %.2f nH' % ltr_loop)
    log('  ltr_capA (30uF bank trace)     = %.2f nH' % ltr_capA)
    log('  lhot     (hot-loop trace)      = %.2f nH' % lhot)
    log('  (v2 nominal had lf1=5, ltr_loop=15, ltr_capA=3, lhot=8 nH;')
    log('   v2 pessimistic had 10 / 30 / 8 / 20 nH -> clear over-estimation)')

    with open(OUT_DIR + 'geo_numbers.txt', 'w') as fh:
        fh.write('\n'.join(out) + '\n')

    # return dict for the EM model
    return dict(Lpm=Lpm, lf1=lf1, ltr_loop=ltr_loop, ltr_capA=ltr_capA, lhot=lhot,
                input_pour_mm=mm(s['input_pour_mil']), segs=segs)


import os
OUT_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'

if __name__ == '__main__':
    main()
