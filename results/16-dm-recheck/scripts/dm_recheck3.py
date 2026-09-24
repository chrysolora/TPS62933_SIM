#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dm_recheck3.py -- final consolidated figure + numbers for DM robustness recheck."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/raid10/sim-work/tps62933/dm_recheck/dm_recheck.py').read().split('# ============================================================================')[0])
OUT='/mnt/raid10/sim-work/tps62933/dm_recheck/'
lines=[]
def log(s=''): print(s); lines.append(s)

Pn=mk(P_nominal()); Pp=mk(P_pessimistic())
def wm(P,tr=10e-9): return worst_margin('pi',P,3.0,0.5,tr,tr)[0]
def sweep(key,vals,tr=10e-9):
    return np.array([wm({**mk(Pn),key:v},tr) for v in vals])

# sweeps
v_lh=np.linspace(0,60e-9,25); y_lh=sweep('lhot',v_lh)
v_lc=np.linspace(0,40e-9,25); y_lc=sweep('ltr_capA',v_lc)
v_es=np.linspace(0.3e-9,3.0e-9,15); y_es=sweep('esl_10u',v_es)
v_tr=np.linspace(1e-9,30e-9,20); y_tr=np.array([wm(Pn,t) for t in v_tr])
v_li=np.array([25.,50.,75.,100.]); y_li=sweep('rs_dm',v_li)

# 2D map lhot x ltr_capA
lhg=np.linspace(0,45e-9,13); lcg=np.linspace(0,25e-9,13)
Z=np.zeros((len(lhg),len(lcg)))
for i,a in enumerate(lhg):
    for j,b in enumerate(lcg):
        Z[i,j]=wm({**mk(Pn),'lhot':a,'ltr_capA':b})

# real source
d=np.loadtxt(OUT+'src_check.txt'); t=d[:,0]; ic=d[:,3]
sel=(t>=15e-6)&(t<=45e-6); ts=t[sel]; ics=ic[sel]
tu=np.arange(ts[0],ts[-1],2e-9); icu=np.interp(tu,ts,ics)
x=icu-icu.mean(); n=len(x); w=np.hanning(n); X=np.fft.rfft(x*w)
fr=np.fft.rfftfreq(n,2e-9); aa=2*np.abs(X)/w.sum()
KMAX=int(30e6/FSW); ks=np.arange(1,KMAX+1); fh=ks*FSW
Ire=np.array([np.interp(k*FSW,fr,aa) for k in ks])
f_an,a_an=trapezoid_harmonics(FSW,3.0,0.5,10e-9,10e-9,KMAX)

# summary bars
def marg(P,I,tr=10e-9,lim=cispr_B,net='pi'):
    z=Zt_vec(net,fh,P); dbv=20*np.log10(np.maximum(I*z,1e-30)/1e-6)
    m=lim(fh)-dbv; j=int(np.argmin(m)); return m[j],fh[j]
Pg=mk(Pn); Pg.update(lhot=15e-9,ltr_capA=8e-9,esl_10u=2.0e-9)
bars_lab=['nominal QP','nominal AV','pessim QP','pessim AV','geomworst QP','geomworst AV']
bars_val=[marg(Pn,a_an)[0],marg(Pn,a_an,lim=cispr_B_avg)[0],
          marg(Pp,a_an,20e-9)[0],marg(Pp,a_an,20e-9,lim=cispr_B_avg)[0],
          marg(Pg,a_an)[0],marg(Pg,a_an,lim=cispr_B_avg)[0]]
bars_col=['#1f77b4','#9ecae1','#ff7f0e','#ffbb78','#d62728','#f4a3a3']

fig,axs=plt.subplots(2,4,figsize=(20,9)); axs=axs.ravel()
def pl(ax,x,y,xl,sc=1e9,note=True):
    ax.plot(x*sc,y,'-o',c='#1f77b4'); ax.axhline(0,color='r',ls='--',lw=1.2)
    ax.set_xlabel(xl); ax.set_ylabel('worst Class-B QP margin [dB]'); ax.grid(alpha=.3)
pl(axs[0],v_lh,y_lh,'hot-loop trace L [nH]')
for v,c,lb in [(8,'k','nom 8'),(20,'gray','pess 20'),(15,'g','geom-worst 15')]:
    axs[0].axvline(v,color=c,ls=':'); axs[0].text(v,30,lb,rotation=90,fontsize=7)
pl(axs[1],v_lc,y_lc,'30uF bank trace L [nH]')
for v,c in [(3,'k'),(8,'gray')]: axs[1].axvline(v,color=c,ls=':')
pl(axs[2],v_es,y_es,'10uF ESL [nH]')
pl(axs[3],v_tr,y_tr,'tr=tf [ns]')
pl(axs[4],v_li,y_li,'LISN DM Z [ohm]',sc=1.0)
axs[5].semilogx(fh/1e6,20*np.log10(np.maximum(a_an,1e-12)/1e-6),'-o',ms=3,label='analytic')
axs[5].semilogx(fh/1e6,20*np.log10(np.maximum(Ire,1e-12)/1e-6),'-s',ms=3,label='real SPICE')
axs[5].set_ylim(-20,140); axs[5].set_xlabel('MHz'); axs[5].set_ylabel('I$_{src}$ harmonic [dB$\\mu$A]')
axs[5].set_title('source: analytic vs real'); axs[5].grid(alpha=.3,which='both'); axs[5].legend(fontsize=8)
im=axs[6].contourf(lcg*1e9,lhg*1e9,Z,levels=np.linspace(-15,40,12),cmap='RdYlGn')
axs[6].contour(lcg*1e9,lhg*1e9,Z,levels=[0],colors='k',linewidths=2)
axs[6].plot(3,8,'k*',ms=15); axs[6].plot(8,20,'ks',ms=8)
axs[6].set_xlabel('30uF bank trace L [nH]'); axs[6].set_ylabel('hot-loop trace L [nH]')
axs[6].set_title('QP margin map (*nominal, sq=pess)'); plt.colorbar(im,ax=axs[6])
axs[7].bar(range(6),bars_val,color=bars_col); axs[7].axhline(0,color='k',lw=1.5)
axs[7].set_xticks(range(6)); axs[7].set_xticklabels(bars_lab,rotation=35,fontsize=7.5)
axs[7].set_ylabel('worst margin [dB]'); axs[7].set_title('summary (analytic src)'); axs[7].grid(alpha=.3,axis='y')
fig.suptitle('TPS62933 DM conducted-EMI robustness recheck (CISPR 32 Class B, 150k-30MHz, full load)',fontsize=13)
fig.tight_layout(); fig.savefig(OUT+'fig_dm_robustness.png',dpi=135)
log('Wrote fig_dm_robustness.png')
log('')
log('2D map: pass/fail boundary shown (0 dB). nominal(3,8)->%+.1f dB; pess(8,20)->%+.1f dB'%(Z[np.argmin(abs(lhg-8))][np.argmin(abs(lcg-3))],Z[np.argmin(abs(lhg-20))][np.argmin(abs(lcg-8))]))
log('frac of grid passing QP: %.0f%%'%(100*(Z>0).mean()))
log('geom-worst(15,8)=%+.1f dB'%Z[np.argmin(abs(lhg-15))][np.argmin(abs(lcg-8))])
open(OUT+'numbers_dm_recheck.txt','a').write('\n'.join(lines)+'\n')
