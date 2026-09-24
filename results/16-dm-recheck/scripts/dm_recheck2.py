#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dm_recheck2.py -- use the REAL SPICE IC input current as DM source; rerun worst-margin."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/raid10/sim-work/tps62933/dm_recheck/dm_recheck.py').read().split('# ============================================================================')[0])
OUT='/mnt/raid10/sim-work/tps62933/dm_recheck/'
FSW=805e3
lines=[]
def log(s=''): print(s); lines.append(s)

# real IC input current harmonics
d=np.loadtxt('/mnt/raid10/sim-work/tps62933/dm_recheck/src_check.txt')
t=d[:,0]; ic=d[:,3]
sel=(t>=15e-6)&(t<=45e-6); ts=t[sel]; ics=ic[sel]
tu=np.arange(ts[0],ts[-1],2e-9); icu=np.interp(tu,ts,ics)
x=icu-icu.mean(); n=len(x); w=np.hanning(n)
X=np.fft.rfft(x*w); fr=np.fft.rfftfreq(n,2e-9); a=2*np.abs(X)/w.sum()
KMAX=int(30e6/FSW); ks=np.arange(1,KMAX+1)
fh=ks*FSW
Ire=np.array([np.interp(k*FSW,fr,a) for k in ks])
f_an, a_an = trapezoid_harmonics(FSW,3.0,0.5,10e-9,10e-9,KMAX)

log('REAL-SOURCE DM CHECK (uses actual SPICE IC input current harmonics)')
log('  k   f[MHz]  I_analytic[A]  I_real[A]  diff[dB]')
for k in range(9):
    dif=20*np.log10(max(Ire[k],1e-12)/max(a_an[k],1e-12))
    log('  %2d  %6.2f   %10.4f  %9.4f  %+7.1f'%(k+1,fh[k]/1e6,a_an[k],Ire[k],dif))
log('  (even harmonics: analytic=0 by symmetry, real=nonzero -> source idealization)')

def dm_from_source(P, Iamp, tr=10e-9, lim=cispr_B, net='pi'):
    z=Zt_vec(net,fh,P)
    dbv=20*np.log10(np.maximum(Iamp*z,1e-30)/1e-6)
    m=lim(fh)-dbv; j=int(np.argmin(m))
    return m[j],fh[j],dbv[j],dbv,m

log('')
log('  WORST MARGIN (full load) using REAL source vs ANALYTIC source:')
for tag,P,tr,se in [('nominal',mk(P_nominal()),10e-9,'nom'),
                    ('pessimistic',mk(P_pessimistic()),20e-9,'pes')]:
    mr=dm_from_source(P,Ire,tr); ma=dm_from_source(P,a_an,tr)
    log('   %-12s analytic-src %+6.1f dB @%.2fMHz | real-src %+6.1f dB @%.2fMHz'
        %(tag,ma[0],ma[1]/1e6,mr[0],mr[1]/1e6))
# geom-worst combined
Pg=mk(P_nominal()); Pg.update(lhot=20e-9,ltr_capA=12e-9,esl_10u=2.0e-9)
mr=dm_from_source(Pg,Ire); ma=dm_from_source(Pg,a_an)
log('   geom-worst   analytic-src %+6.1f dB @%.2fMHz | real-src %+6.1f dB @%.2fMHz'
    %(ma[0],ma[1]/1e6,mr[0],mr[1]/1e6))

# ---- figure: 3 panels ----
fig,ax=plt.subplots(1,3,figsize=(16,5.4))
fq=np.logspace(np.log10(1.5e5),np.log10(3e7),400)
ax[0].semilogx(fh/1e6,20*np.log10(np.maximum(a_an,1e-12)/1e-6),'-o',ms=3,label='analytic trapezoid')
ax[0].semilogx(fh/1e6,20*np.log10(np.maximum(Ire,1e-12)/1e-6),'-s',ms=3,label='real SPICE i(VIN pin)')
ax[0].set_ylim(-20,140); ax[0].set_xlim(0.5,30)
ax[0].set_xlabel('MHz'); ax[0].set_ylabel('IC input current harmonic [dB$\\mu$A]')
ax[0].set_title('(a) DM source spectrum: analytic vs real'); ax[0].grid(alpha=.3,which='both'); ax[0].legend(fontsize=8)
_,_,_,dbv_r,_=dm_from_source(mk(P_nominal()),Ire)
_,_,_,dbv_a,_=dm_from_source(mk(P_nominal()),a_an)
_,_,_,dbv_rp,_=dm_from_source(mk(P_pessimistic()),Ire,20e-9)
ax[1].plot(fq/1e6,cispr_B(fq),'k-',lw=2,label='CISPR 32 B QP')
ax[1].plot(fq/1e6,cispr_B_avg(fq),'k--',lw=1.3,label='CISPR 32 B AV')
ax[1].semilogx(fh/1e6,dbv_a,'-o',ms=3,label='nominal (analytic src)')
ax[1].semilogx(fh/1e6,dbv_r,'-s',ms=3,label='nominal (real src)')
ax[1].semilogx(fh/1e6,dbv_rp,'-^',ms=3,color='r',label='pessimistic (real src)')
ax[1].set_xlim(0.5,30); ax[1].set_ylim(-10,90); ax[1].set_xlabel('MHz'); ax[1].set_ylabel('DM on LISN [dB$\\mu$V]')
ax[1].set_title('(b) DM spectrum, real vs analytic source'); ax[1].grid(alpha=.3,which='both'); ax[1].legend(fontsize=8)
# panel c: margin ranges
labels=['nominal\nanalytic','nominal\nreal','pess.\nanalytic','pess.\nreal','geom-worst\nreal','geom-worst\nanalytic']
vals=[dm_from_source(mk(P_nominal()),a_an)[0],dm_from_source(mk(P_nominal()),Ire)[0],
      dm_from_source(mk(P_pessimistic()),a_an,20e-9)[0],dm_from_source(mk(P_pessimistic()),Ire,20e-9)[0],
      dm_from_source(Pg,Ire)[0],dm_from_source(Pg,a_an)[0]]
cols=['#1f77b4']*2+['#ff7f0e']*2+['#d62728']*2
ax[2].bar(range(6),vals,color=cols); ax[2].axhline(0,color='k',lw=1.5)
ax[2].set_xticks(range(6)); ax[2].set_xticklabels(labels,fontsize=7.5,rotation=25)
ax[2].set_ylabel('worst Class-B QP margin [dB]'); ax[2].set_title('(c) margin summary')
ax[2].grid(alpha=.3,axis='y')
fig.tight_layout(); fig.savefig(OUT+'fig_dm_robustness_src.png',dpi=140)
log(''); log('Wrote fig_dm_robustness_src.png')
open(OUT+'numbers_dm_recheck_src.txt','w').write('\n'.join(lines)+'\n')
